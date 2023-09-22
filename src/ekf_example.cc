// EKF related headers
#include <ekf_15dof_class.h>
// Autocoded FCS Libraries
#include <fcsModel.h>
#include <fcs_params.h>
// Matrix library
#include <Matrix/matrix_factorization_class.h>
// Navio 2 Utilities
#include "Navio/Common/MPU9250.h"
#include "Navio/Navio2/LSM9DS1.h"
#include "Navio/Common/Util.h"
#include <gps_utils.h>
#include <write_utils.h>
#include <imu_utils.h>
#include <rc_input_utils.h>
#include <pwm_output_utils.h>
// Standard C++ Libraries file, time and memory
#include <memory>

#include <chrono>
// Standard C++ Libraries for multi-threading
#include <thread>
#include <pthread.h>
#include <mutex>
#include <atomic>
#include <signal.h>
#include <sys/mman.h>
#include <iterator>

static fcsModel fcsModel_Obj;          // Instance of FCS model class

// To catch SIGINT
volatile sig_atomic_t sigint_flag = 0;

void sigint_handler(int sig){ // can be called asynchronously
  sigint_flag = 1; // set flag
}

int main(int argc, char *argv[]){
	system("sudo echo -1 >/proc/sys/kernel/sched_rt_runtime_us");

	GpsHelper gps_reader;

	WriteHelper data_writer("data_file.dat");

	ImuHelper imu_reader("mpu");
	imu_reader.InitializeImu();


	RcInputHelper rc_reader(8);
	rc_reader.InitializeRcInput();

	PwmOutputHelper pwm_writer(4);
	pwm_writer.InitializePwmOutput();

	vector<float> pwm_out_val(4, 0);

	float ned_pos_and_vel_meas[6];
	bool gps_meas_indices[6];
	float raw_lat_lon_alt[3];
	
	sched_param sch;
	int policy;	

	// Variables to set CPU affinity
	cpu_set_t cpuset;
  CPU_ZERO(&cpuset);
  CPU_SET(3, &cpuset);

	pthread_getschedparam(pthread_self(), &policy, &sch);
	sch.sched_priority = sched_get_priority_max(SCHED_FIFO);
	pthread_setschedparam(pthread_self(), SCHED_FIFO, &sch);
	int rc = pthread_setaffinity_np(pthread_self(),
                                    sizeof(cpu_set_t), &cpuset);    
	if (rc != 0) {
      std::cerr << "Error calling pthread_setaffinity_np on main(): " << rc << "\n";
    }

	// Register signals 
  	signal(SIGINT, sigint_handler); 

    // Variables to read data from the IMU
    // Accels
    float accel[3];
    //Gyros
    float gyro[3];
    //Mags
    float mag[3];

  // Initial State Variable
  MatrixInv<float> initial_state(15, 1);
    // Variable used in the measurement update of the EKF
	MatrixInv<float> sensor_meas(9, 1);
	// Sensor values used in the time propagation stage of the EKF
	MatrixInv<float> state_sensor_val(6, 1);
	// Q matrix of the EKF
	MatrixInv<float> process_noise_q(15, 15, "eye");
	// P matrix initial
	MatrixInv<float> initial_covariance_p(15, 15, "eye");
	// Angle process noise
	process_noise_q.Diag({0.000001, 0.000001, 0.000001, 0.000000001, 0.000000001, 0.000000001,
						  0.000000001, 0.000000001, 0.000000001, 0.000001, 0.000001, 0.000001,
						  0.000000001, 0.000000001, 0.000000001});
	// R matrix of the EKF
	MatrixInv<float> meas_noise_r(7, 7, "eye");
	meas_noise_r.Diag({0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01});
	// P init
	initial_covariance_p.Diag({30*DEG2RAD, 30*DEG2RAD, 30*DEG2RAD, 0.01, 0.01, 0.01, 100, 100, 100,
							   10, 10, 10, 0.1, 0.1, 0.1});

	// Variable to store current state at the end of each EKF run
	MatrixInv<float> current_state(15, 1);

	// NED to BODY DCM
	MatrixInv<float> c_ned2b;
	// NED to FEP DCM
	MatrixInv<float> c_ned2fep;

	// Array to indicate which measurement has been update
	bool meas_indices[] = {false,
						   false, false, false,
						   false, false, false};

	// Magnetic field declination in the Bay Area
	float magnetic_declination = 13.01*DEG2RAD;

	// Calibration variables for mag
	// Offset
	MatrixInv<float> mag_offset = {{17.8902136639429}, {35.5186740453011}, {-33.8067089624238}};
	// Rotation matrix for misalignment of the IMU axes
	MatrixInv<float> mag_a(3, 3, "eye");
	// Scale factor for each axis of the IMU
	MatrixInv<float> mag_scale = {{0.9652, 0, 0}, {0, 1.09, 0}, {0, 0, 0.9556}};	
	imu_reader.SetMagParams(magnetic_declination , mag_a, mag_offset, mag_scale);

	// Initial attitude and NED velocity
	float* init_att = imu_reader.ComputeInitialRollPitchAndYaw(200);

	gps_reader.InitializeGps(30);
	gps_reader.CreateGpsThread();
	data_writer.StartFileWriteThread();
	rc_reader.CreateRcInputReadThread();

	double* vned_init = gps_reader.GetInitNedVel();	

	for(size_t i_idx = 0; i_idx < 3; i_idx++){
		initial_state(i_idx) = init_att[i_idx];
		initial_state(i_idx + 9) = vned_init[i_idx];
	}

	// Create an 15 state EKF object
	Ekf15Dof<float> imu_gps_ekf(0.004, initial_state, process_noise_q, meas_noise_r, initial_covariance_p);

	// Make sure EKF always uses Magnetometer in the update step
	meas_indices[0] = true;

	// Variable to read imu data
	float* imu_data;

	//Variable to read rc input data
	int* rc_periods;

	//##############################################################
	// FCS Ctrl Params from fcs_params.h
	AssignFcsCtrlParams();
	// fcsModel input variable
	fcsModel::ExtU_fcsModel_T *ExtU_fcsModel_T_ =  new fcsModel::ExtU_fcsModel_T;
	
	// Initialize model
  fcsModel_Obj.initialize();
  // RC Cmd Input Variables
  busRcInCmds rcCmdsIn_;
  // Sensor inuput to the model
  busStateEstimate stateEstimate_;
  // From fcs_params.h
  ExtU_fcsModel_T_->ctrlParams = fcs_ctrl_params;

  // fcsModel output variable
  fcsModel::ExtY_fcsModel_T ExtY_fcsModel_T_;
  // Mutex to guard resource access to fcs outputs
	mutex fcs_out_mutex;
  {
  	unique_lock<mutex> fcs_out_lock(fcs_out_mutex);	
		ExtY_fcsModel_T_ = fcsModel_Obj.getExternalOutputs();
	}

	
  //##############################################################

  //write the initial state
	float zero_array[9] = {0};
	int rc_periods_ph[7];
	rc_periods_ph[0] = -1;
	rc_periods_ph[1] = -1;
	rc_periods_ph[2] = -1;
	rc_periods_ph[3] = -1;
	rc_periods_ph[4] = -1;
	rc_periods_ph[5] = -1;
	rc_periods_ph[6] = -1;
	data_writer.UpdateDataBuffer(0, 0, zero_array, sensor_meas, initial_state, gps_meas_indices, rc_periods_ph, ExtY_fcsModel_T_);

	// Loop counter
	size_t loop_count = 0;
	// initialize the duration
	chrono::microseconds delta (4000); 
	auto duration = chrono::duration_cast<chrono::microseconds> (delta);
	// start time
	auto time_point_start = chrono::high_resolution_clock::now();
	chrono::microseconds time_start_us = chrono::duration_cast<chrono::microseconds>(time_point_start.time_since_epoch());
	auto time_start_us_count = time_start_us.count();
	// loop
    while(1) {
    	auto loop_start = chrono::high_resolution_clock::now();
    	auto time_since_loop_start_us = chrono::duration_cast<chrono::microseconds>(loop_start.time_since_epoch()).count();
    	while (  time_since_loop_start_us < (time_start_us_count + 4000) ){
    		time_since_loop_start_us = chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now().time_since_epoch()).count();
    	}

    	// Make all the measurement flags corresponding to position and velocity false
			for( size_t idx_meas = 1; idx_meas < 7; idx_meas++ ){
    		meas_indices[idx_meas] = false;
    	}

    	// Read IMU data
	    imu_data = imu_reader.GetImuData();

	    // Assign required sensor value for time propagation of the ekf state and measurement update
	    for(size_t imu_idx = 0; imu_idx < 3; imu_idx++){
	    	state_sensor_val(imu_idx) = imu_data[imu_idx];
	    	state_sensor_val(imu_idx + 3) = imu_data[imu_idx + 3];
	    	sensor_meas(imu_idx) =  imu_data[imu_idx + 6];
	    }

	    gps_reader.GetRawLatLonAlt(raw_lat_lon_alt);
	    gps_reader.GetGpsNedPosAndVel(ned_pos_and_vel_meas, gps_meas_indices);
	    // Check if GPS position and velocity has been update and set appropriate flags
	    for(size_t idx_meas = 0; idx_meas < 6; idx_meas++){
	    	if(gps_meas_indices[idx_meas]){
	    		sensor_meas(idx_meas + 3) = ned_pos_and_vel_meas[idx_meas];
	    		meas_indices[idx_meas + 1] = true;
	    	}
	    }

	   	// Run one step of EKF
	    imu_gps_ekf.Run(state_sensor_val, sensor_meas, meas_indices);
	    // Get the state after EKF run
      current_state = imu_gps_ekf.GetCurrentState();

      //########################################
      // Assign the state values to the model input structure
      for(size_t idx = 0; idx < 3; idx++){
      	stateEstimate_.attitude_rad[idx] = current_state(idx);
      	stateEstimate_.bodyAngRates_radps[idx] = ( state_sensor_val(idx) - current_state(idx + 3) );
      	stateEstimate_.nedVel_mps[idx] = current_state(idx + 9);
      }

      // Get NED to Body DCM
      c_ned2b = GetDcm(current_state(0), current_state(1), current_state(2));

      // //GET NED to FEP DCM
      c_ned2fep = GetDcm(current_state(0), current_state(1), 0);

      for(size_t idx = 0; idx < 3; idx++){
      	for(size_t jdx = 0; jdx < 3; jdx++){
      		stateEstimate_.ned2BodyDcm_nd[idx*3 + jdx] = c_ned2b(jdx, idx);
      		stateEstimate_.ned2FepDcm_nd[idx*3 +  jdx] = c_ned2fep(jdx, idx);
      	}
      }


      stateEstimate_.geodeticPos.lat_rad = raw_lat_lon_alt[0];
      stateEstimate_.geodeticPos.lon_rad = raw_lat_lon_alt[1];
      stateEstimate_.geodeticPos.alt_m = current_state(8);

      ExtU_fcsModel_T_->stateEstimate = stateEstimate_;

      // Get RC Data
      if (remainder(loop_count, 5) == 0){
      	rc_periods = rc_reader.GetRcPeriods();
      	rcCmdsIn_.throttleCmd_nd = rc_periods[2];
      	rcCmdsIn_.joystickYCmd_nd = rc_periods[1];
      	rcCmdsIn_.joystickXCmd_nd = rc_periods[0];
      	rcCmdsIn_.joystickZCmd_nd = rc_periods[3];
      	rcCmdsIn_.rcSwitch1_nd = rc_periods[4];
      	rcCmdsIn_.rcSwitch2_nd = rc_periods[5];
      	rcCmdsIn_.rcSwitch3_nd = rc_periods[6];


      }

      ExtU_fcsModel_T_->rcCmdsIn = rcCmdsIn_;
      fcsModel_Obj.setExternalInputs(ExtU_fcsModel_T_);
      //########################################

      // Step the FCS model
      {
      	unique_lock<mutex> fcs_out_lock(fcs_out_mutex);
  			fcsModel_Obj.step();
  			ExtY_fcsModel_T_ = fcsModel_Obj.getExternalOutputs();
  		}


      if(remainder(loop_count, 5) == 0){
        	data_writer.UpdateDataBuffer(duration.count(), loop_count, imu_data, sensor_meas, current_state, gps_meas_indices, rc_periods, ExtY_fcsModel_T_);
	    }

	  //   if (remainder(loop_count, 50) == 0){
	  //   	printf("Roll [deg]: %+7.3f, Pitch[deg]: %+7.3f, Yaw[deg]: %+7.3f\n", current_state(0)*RAD2DEG, current_state(1)*RAD2DEG, current_state(2)*RAD2DEG);
	  //   	// printf("Pos N [m]: %+7.3f, Pos E [m]: %+7.3f, Pos D[m]: %+7.3f\n", current_state(6), current_state(7), current_state(8));
	  //   	// printf("Vel N [m]: %+7.3f, Vel E [m]: %+7.3f, Vel D[m]: %+7.3f\n", current_state(9), current_state(10), current_state(11));
	  //   	// printf("Throttle: %d, Roll: %d, Pitch: %d, Yaw: %d, Sw1: %d, Sw2: %d, Sw3: %d, State: %d\n", ExtU_fcsModel_T_->rcCmdsIn.throttleCmd_nd, ExtU_fcsModel_T_->rcCmdsIn.joystickXCmd_nd, 
	  //   	// 	ExtU_fcsModel_T_->rcCmdsIn.joystickYCmd_nd, ExtU_fcsModel_T_->rcCmdsIn.joystickZCmd_nd, ExtU_fcsModel_T_->rcCmdsIn.rcSwitch1_nd, ExtU_fcsModel_T_->rcCmdsIn.rcSwitch2_nd, 
	  //   	// 	ExtU_fcsModel_T_->rcCmdsIn.rcSwitch3_nd, ExtY_fcsModel_T_.fcsDebug.state);
	  //   	// // printf("%g, %g, %g, %g, %d\n",ExtY_fcsModel_T_.actuatorsCmds[0], ExtY_fcsModel_T_.actuatorsCmds[1], ExtY_fcsModel_T_.actuatorsCmds[2], ExtY_fcsModel_T_.actuatorsCmds[3], ExtY_fcsModel_T_.fcsDebug.state);
	  //   	printf("############################################\n");
		// }

		if(static_cast<uint8_t>(ExtY_fcsModel_T_.fcsDebug.state) != 0){
			if(rcCmdsIn_.throttleCmd_nd <= PWM_CHECK_MIN_THRESHOLD){
				pwm_out_val[0] = static_cast<float>(rcCmdsIn_.throttleCmd_nd*1.0);
				pwm_out_val[1] = static_cast<float>(rcCmdsIn_.throttleCmd_nd*1.0);
				pwm_out_val[2] = static_cast<float>(rcCmdsIn_.throttleCmd_nd*1.0);
				pwm_out_val[3] = static_cast<float>(rcCmdsIn_.throttleCmd_nd*1.0);
			}else{
		  	pwm_out_val[0] = max(PWM_CMD_MIN_THRESHOLD, min(PWM_CMD_MAX_THRESHOLD, ExtY_fcsModel_T_.actuatorsCmds[0]*RPM_TO_PWM_SCALE + PWM_MIN_THRESHOLD));
		  	pwm_out_val[1] = max(PWM_CMD_MIN_THRESHOLD, min(PWM_CMD_MAX_THRESHOLD, ExtY_fcsModel_T_.actuatorsCmds[3]*RPM_TO_PWM_SCALE + PWM_MIN_THRESHOLD));
		  	pwm_out_val[2] = max(PWM_CMD_MIN_THRESHOLD, min(PWM_CMD_MAX_THRESHOLD, ExtY_fcsModel_T_.actuatorsCmds[1]*RPM_TO_PWM_SCALE + PWM_MIN_THRESHOLD));
		  	pwm_out_val[3] = max(PWM_CMD_MIN_THRESHOLD, min(PWM_CMD_MAX_THRESHOLD, ExtY_fcsModel_T_.actuatorsCmds[2]*RPM_TO_PWM_SCALE + PWM_MIN_THRESHOLD));
			}
		}else{
				pwm_out_val[0] = ExtY_fcsModel_T_.actuatorsCmds[0]*1.0;
		  	pwm_out_val[1] = ExtY_fcsModel_T_.actuatorsCmds[3]*1.0;
		  	pwm_out_val[2] = ExtY_fcsModel_T_.actuatorsCmds[1]*1.0;
		  	pwm_out_val[3] = ExtY_fcsModel_T_.actuatorsCmds[2]*1.0;
		}
	  // if (remainder(loop_count, 50) == 0){
	  // 	printf("%g, %g, %g, %g\n", ExtY_fcsModel_T_.actuatorsCmds[0], ExtY_fcsModel_T_.actuatorsCmds[1], ExtY_fcsModel_T_.actuatorsCmds[2], ExtY_fcsModel_T_.actuatorsCmds[3]);
	  // 	// printf("%g, %g, %g, %g\n", pwm_out_val[0], pwm_out_val[1], pwm_out_val[2], pwm_out_val[3]);
	  // }

	  if( ExtY_fcsModel_T_.actuatorsCmds[0] > 1.0 && ExtY_fcsModel_T_.actuatorsCmds[1] > 1.0 && ExtY_fcsModel_T_.actuatorsCmds[2] > 1.0 &&
	  	ExtY_fcsModel_T_.actuatorsCmds[3] > 1.0 ){
	  	pwm_writer.SetPwmDutyCyle(pwm_out_val);
	}

		loop_count++;

	    // Get the stop time and compute the duration
	    auto loop_end = std::chrono::high_resolution_clock::now();

	    duration = chrono::duration_cast<chrono::microseconds>(loop_end - loop_start);
	    time_start_us_count = time_start_us_count + 4000;
	   	//cout<<loop_count<<"\n";
	   	if(sigint_flag == 1)
	   		break;
	}
	return 0;
}