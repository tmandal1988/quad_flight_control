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
#include <baro_utils.h>
#include <mahony_filter.h>
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

/**************************Flags to use EKF or/and Mahony*******************************/
bool use_mahony_filter = false;
bool use_ekf = true;
/**************************Flags to use EKF or/and Mahony*******************************/

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

	BaroHelper baro_reader;
	float baro_debug_data[9];
	array<float, 3> baro_accel{0};
	array<float, 3> baro_euler{0};
	baro_reader.SetKalmanFilterParams(array<float, 3> {100.0, 100.0, 100.0}, array<float, 3> {1e-6, 1e-6, 1e-6}, 
										  array<float, 2>{0.01, 0.09});
	baro_reader.StartBaroReader(1, 20);	
	 // Mutex to guard resource access to baro data
	mutex baro_out_mutex;

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

  /****************Variables to read data from the IMU*************************************/
  // Accels
  float accel[3];
  //Gyros
  float gyro[3];
  float gyro_offset[3];
  //Mags
  float mag[3];
  /****************Variables to read data from the IMU*************************************/

  /****************Variables to read data from the Baro*************************************/
  float baro_data[2];
  /****************Variables to read data from the Baro*************************************/

	//Quat
	float quat[4] = {0};
	quat[0] = 1;
	float mh_euler[3] = {0};

	/****************SAMPLE TIME VARIABLES*************************************/
	//Sample time
	float dt_s = 0.004;
	//Sample Step
	size_t dt_count = 4000;
	/****************SAMPLE TIME VARIABLES*************************************/


  // Initial State Variable for EKF
	MatrixInv<float> initial_state(15, 1);
	// Variable to store current state at the end of each EKF run
	MatrixInv<float> current_state(15, 1);
  
  // Variable used in the measurement update of the EKF
	MatrixInv<float> sensor_meas(20, 1);
	// Sensor values used in the time propagation stage of the EKF
	MatrixInv<float> state_sensor_val(6, 1);
	
	// P matrix initial
	MatrixInv<float> initial_covariance_p(15, 15, "eye");

	// Q matrix of the EKF
	// MatrixInv<float> process_noise_q(6, 6, "eye");
	// process_noise_q.Diag({1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5});
	MatrixInv<float> process_noise_q(15, 15, "eye");
	process_noise_q.Diag({0.000001, 0.000001, 0.000001, 1e-8, 1e-8, 1e-8, 1e-8, 1e-8,
    									  1e-8, 0.005, 0.005, 0.005, 1e-8, 1e-8, 1e-8});
	// R matrix of the EKF
	MatrixInv<float> meas_noise_r(7, 7, "eye");
	meas_noise_r.Diag({0.5, 10, 10, 10, 0.01, 0.01, 0.01});
	//meas_noise_r.Diag({0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001});
	// P init
	initial_covariance_p.Diag({30*DEG2RAD, 30*DEG2RAD, 30*DEG2RAD, 0.01, 0.01, 0.01, 100, 100, 100,
							   10, 10, 10, 0.1, 0.1, 0.1});	

	/****************SECONDARY FILTER DEBUG*************************************/
	MatrixInv<float> secondary_filter_debug(3, 1);
	/****************SECONDARY FILTER DEBUG*************************************/


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
	imu_reader.ComputeGyroOffset(100);
	imu_reader.GetGyroOffset(gyro_offset);
	float* init_att = imu_reader.ComputeInitialRollPitchAndYaw(200);
	float init_att_array[3];
	init_att_array[0] = init_att[0];
	init_att_array[1] = init_att[1];
	init_att_array[2] = init_att[2];

  GetQuatFromEuler(init_att_array, quat);

  MahonyFilter m_filt(0.0, 2.0, dt_s, quat);	


	bool gps_init_status = gps_reader.InitializeGps(30);

	if(!gps_init_status){
		use_mahony_filter = true;
		printf("EKF needs GPS, using Mahony filter for attitude computation. ONLY USE STABILIZE MODE\n");
	}
	gps_reader.CreateGpsThread();
	data_writer.StartFileWriteThread();
	rc_reader.CreateRcInputReadThread();

	double* vned_init = gps_reader.GetInitNedVel();	

	for(size_t i_idx = 0; i_idx < 3; i_idx++){
		initial_state(i_idx) = init_att[i_idx];
		initial_state(i_idx + 9) = vned_init[i_idx];
	}

	// Create an 15 state EKF object
	Ekf15Dof<float> imu_gps_ekf(dt_s, initial_state, process_noise_q, meas_noise_r, initial_covariance_p);

	// Make sure EKF always uses Magnetometer in the update step
	meas_indices[0] = true;

	// Variable to read imu data
	float* imu_data;

	//Variable to read rc input data
	int* rc_periods = new int[7];

	// Assign -1 to rc_periods to initialize
	for(size_t rc_idx = 0; rc_idx < 7; rc_idx++){
		rc_periods[rc_idx] = -1;
	}

	//##############################################################
	// fcsModel input variable
	fcsModel::ExtU_fcsModel_T *ExtU_fcsModel_T_ =  new fcsModel::ExtU_fcsModel_T;
	
	// Initialize model
  fcsModel_Obj.initialize();
  // RC Cmd Input Variables
  busRcInCmds rcCmdsIn_;
  // Sensor inuput to the model
  busStateEstimate stateEstimate_;
  // FCS Ctrl Params from fcs_params.h
  ExtU_fcsModel_T_->ctrlParams = AssignFcsCtrlParams();

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
	// int rc_periods_ph[7];
	// rc_periods_ph[0] = -1;
	// rc_periods_ph[1] = -1;
	// rc_periods_ph[2] = -1;
	// rc_periods_ph[3] = -1;
	// rc_periods_ph[4] = -1;
	// rc_periods_ph[5] = -1;
	// rc_periods_ph[6] = -1;
	data_writer.UpdateDataBuffer(0, 0, zero_array, sensor_meas, initial_state, secondary_filter_debug, gps_meas_indices, rc_periods, ExtY_fcsModel_T_);

	// Loop counter
	size_t loop_count = 0;
	size_t arm_loop_count = 0;
	//125 Hz flag
	bool one_twenty_five_hz_flag = false;
	// 50 Hz flag
	bool fifty_hz_flag = false;
	// 10 Hz flag
	bool ten_hz_flag = false;
	// initialize the duration
	chrono::microseconds delta (dt_count); 
	auto duration = chrono::duration_cast<chrono::microseconds> (delta);
	auto duration_count = duration.count();
	// start time
	auto time_point_start = chrono::high_resolution_clock::now();
	chrono::microseconds time_start_us = chrono::duration_cast<chrono::microseconds>(time_point_start.time_since_epoch());
	auto time_start_us_count = time_start_us.count();

	// Loop timers
	chrono::high_resolution_clock::time_point loop_start;
	chrono::high_resolution_clock::time_point loop_end;
	// long long time_since_loop_start_us;
	// long long loop_start_us;
	bool is_mtr_armed = false;

	// Intermediate trig variables for use in calculations
	float s_phi;
	float s_theta;
	float s_psi;

	float c_phi;
	float c_theta;
	float c_psi;

	// loop
    while(1) {
    	loop_count++;    	

    	/* Check if the current loop count is a multiple of 2 which will give a 125hz loop as main loop
    	runs at 250Hz
    	*/
    	if(loop_count % 2 == 0){
    		one_twenty_five_hz_flag = true;
    	}else{
    		one_twenty_five_hz_flag = false;
    	}

    	/* Check if the current loop count is a multiple of 5 which will give a 50hz loop as main loop
    	runs at 250Hz
    	*/
    	if(loop_count % 5 == 0){
    		fifty_hz_flag = true;
    	}else{
    		fifty_hz_flag = false;
    	}

    	/* Check if the current loop count is a multiple of 25 which will give a 10hz loop as main loop
    	runs at 250Hz
    	*/
    	if(loop_count % 25 == 0){
    		ten_hz_flag = true;
    	}else{
    		ten_hz_flag = false;
    	}

    	/* Get loop start time
    	*/
    	loop_start = chrono::high_resolution_clock::now();

    	// Make all the measurement flags corresponding to position and velocity false
			for( size_t idx_meas = 1; idx_meas < 7; idx_meas++ ){
    		meas_indices[idx_meas] = false;
    	}

    	// Read IMU data
	    imu_data = imu_reader.GetImuData();

	    // Read GPS data
	    gps_reader.GetRawLatLonAlt(raw_lat_lon_alt);
		  gps_reader.GetGpsNedPosAndVel(ned_pos_and_vel_meas, gps_meas_indices);

	    if (use_ekf){
		    // Assign required sensor value for time propagation of the ekf state and measurement update
		    for(size_t imu_idx = 0; imu_idx < 3; imu_idx++){
		    	state_sensor_val(imu_idx) = imu_data[imu_idx];
		    	state_sensor_val(imu_idx + 3) = imu_data[imu_idx + 3];
		    	sensor_meas(imu_idx) =  imu_data[imu_idx + 6];
		    }
		    
		    // Check if GPS position and velocity has been updated and set appropriate flags
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
    	}

    	if (use_mahony_filter){
    		m_filt.MahonyFilter9Dof(imu_data, gyro_offset);
    		m_filt.GetCurrentQuat(quat);
    		GetEulerFromQuat(quat, mh_euler);
    		secondary_filter_debug(0) = mh_euler[0];
    		secondary_filter_debug(1) = mh_euler[1];
    		secondary_filter_debug(2) = mh_euler[2];
    	}

    	// Get NED to Body DCM
      c_ned2b = GetDcm(current_state(0), current_state(1), current_state(2));

      // //GET NED to FEP DCM
      c_ned2fep = GetDcm(current_state(0), current_state(1), 0);

    	// Read Baro data at 50 Hz
    	if (fifty_hz_flag){
    			{
    				unique_lock<mutex> baro_out_lock(baro_out_mutex);
						baro_reader.GetBaroPressAndTemp(baro_data);
						baro_reader.GetBaroDebugData(baro_debug_data);
					}
					sensor_meas(9) = baro_data[0];
					sensor_meas(10) = baro_data[1];
					stateEstimate_.pressure_mbar = baro_data[0];
					stateEstimate_.temp_c = baro_data[1];					
					for(size_t idx_b = 0; idx_b < 9; idx_b++){
						sensor_meas(11 + idx_b) = baro_debug_data[idx_b];
					}	
					stateEstimate_.aglEst_m = baro_debug_data[6];
					stateEstimate_.climbRateEst_mps = baro_debug_data[7];		

					// Compute NED accels
					s_phi = sin(current_state(0));
					s_theta = sin(current_state(1));
					s_psi = sin(current_state(2));

					c_phi = cos(current_state(0));
					c_theta = cos(current_state(1));
					c_psi = cos(current_state(2));

					stateEstimate_.nedAccel_mps2[0] = imu_data[5]*(s_phi*s_psi + c_phi*c_psi*s_theta) - imu_data[4]*(c_phi*s_psi - 
																						c_psi*s_phi*s_theta) + imu_data[3]*c_psi*c_theta;
					stateEstimate_.nedAccel_mps2[1] = imu_data[4]*(c_phi*c_psi + s_phi*s_psi*s_theta) - imu_data[5]*(c_psi*s_phi - 
																						c_phi*s_psi*s_theta) + imu_data[3]*c_theta*s_psi;
					stateEstimate_.nedAccel_mps2[2] = -baro_debug_data[8];


					if(gps_init_status){
	    			baro_accel[0] = imu_data[3] - current_state(12);
	    			baro_accel[1] = imu_data[4] - current_state(13);
	    			baro_accel[2] = imu_data[5] - current_state(14);

	    			baro_euler[0] = current_state(0);
	    			baro_euler[1] = current_state(1);
	    			baro_euler[2] = current_state(2);
    			}else{      		

	    			baro_accel[0] = imu_data[3];
	    			baro_accel[1] = imu_data[4];
	    			baro_accel[2] = imu_data[5];

	    			baro_euler[0] = mh_euler[0];
	    			baro_euler[1] = mh_euler[1];
	    			baro_euler[2] = mh_euler[2];
    			}
    			{
    				unique_lock<mutex> baro_out_lock(baro_out_mutex);
    				baro_reader.SetBodyAccels(baro_accel);
    				baro_reader.SetEulerAngles(baro_euler);		
    			}				
    	}
    	
      //########################################
      // Assign the state values to the model input structure
      for(size_t idx = 0; idx < 3; idx++){
      	if(gps_init_status){
      		stateEstimate_.attitude_rad[idx] = current_state(idx);
      		stateEstimate_.bodyAngRates_radps[idx] = ( state_sensor_val(idx) - current_state(idx + 3) );
      	}else{      		
      		stateEstimate_.attitude_rad[idx] = mh_euler[idx];
      		stateEstimate_.bodyAngRates_radps[idx] = imu_data[idx] - gyro_offset[idx];
      	}
 
      	stateEstimate_.nedPos_m[idx] = current_state(idx + 6);
      	stateEstimate_.nedVel_mps[idx] = current_state(idx + 9);
      }

      for(size_t idx = 0; idx < 3; idx++){
      	for(size_t jdx = 0; jdx < 3; jdx++){
      		stateEstimate_.ned2BodyDcm_nd[idx*3 + jdx] = c_ned2b(jdx, idx);
      		stateEstimate_.ned2FepDcm_nd[idx*3 +  jdx] = c_ned2fep(jdx, idx);
      	}
      }


      stateEstimate_.geodeticPos.lat_rad = current_state(6);
      stateEstimate_.geodeticPos.lon_rad = current_state(7);
      stateEstimate_.geodeticPos.alt_m = current_state(8);

      ExtU_fcsModel_T_->stateEstimate = stateEstimate_;

      // Get RC Data
      if (fifty_hz_flag){
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

     if(ten_hz_flag){
        data_writer.UpdateDataBuffer(duration_count, loop_count, imu_data, sensor_meas, current_state, secondary_filter_debug, gps_meas_indices, rc_periods, ExtY_fcsModel_T_);
	   }

	    // if (remainder(loop_count, 50) == 0){
	    // 	printf("Roll [deg]: %+7.3f, Pitch[deg]: %+7.3f, Yaw[deg]: %+7.3f\n", current_state(0)*RAD2DEG, current_state(1)*RAD2DEG, current_state(2)*RAD2DEG);
	  // //   	// printf("Pos N [m]: %+7.3f, Pos E [m]: %+7.3f, Pos D[m]: %+7.3f\n", current_state(6), current_state(7), current_state(8));
	  // //   	// printf("Vel N [m]: %+7.3f, Vel E [m]: %+7.3f, Vel D[m]: %+7.3f\n", current_state(9), current_state(10), current_state(11));
	    	// printf("Throttle: %d, Roll: %d, Pitch: %d, Yaw: %d, Sw1: %d, Sw2: %d, Sw3: %d, State: %d, Flight Mode: %d\n", ExtU_fcsModel_T_->rcCmdsIn.throttleCmd_nd, ExtU_fcsModel_T_->rcCmdsIn.joystickXCmd_nd, 
	    	// 	ExtU_fcsModel_T_->rcCmdsIn.joystickYCmd_nd, ExtU_fcsModel_T_->rcCmdsIn.joystickZCmd_nd, ExtU_fcsModel_T_->rcCmdsIn.rcSwitch1_nd, ExtU_fcsModel_T_->rcCmdsIn.rcSwitch2_nd, 
	    	// 	ExtU_fcsModel_T_->rcCmdsIn.rcSwitch3_nd, ExtY_fcsModel_T_.fcsDebug.state, ExtY_fcsModel_T_.fcsDebug.flightMode);
	  // //   	// // printf("%g, %g, %g, %g, %d\n",ExtY_fcsModel_T_.actuatorsCmds[0], ExtY_fcsModel_T_.actuatorsCmds[1], ExtY_fcsModel_T_.actuatorsCmds[2], ExtY_fcsModel_T_.actuatorsCmds[3], ExtY_fcsModel_T_.fcsDebug.state);
	  //   	printf("Throttle: %d, Vz_Cmd: %+7.3f\n",ExtU_fcsModel_T_->rcCmdsIn.throttleCmd_nd, ExtY_fcsModel_T_.fcsDebug.outerLoopCtrlDebug.velCtrlDebug.cmd[2]);
	  //   	printf("############################################\n");
		// }

		if(static_cast<uint8_t>(ExtY_fcsModel_T_.fcsDebug.state) != 0){
				if(rcCmdsIn_.throttleCmd_nd <= PWM_CHECK_MIN_THRESHOLD){
					pwm_out_val[0] = static_cast<float>(rcCmdsIn_.throttleCmd_nd*1.0);
					pwm_out_val[1] = static_cast<float>(rcCmdsIn_.throttleCmd_nd*1.0);
					pwm_out_val[2] = static_cast<float>(rcCmdsIn_.throttleCmd_nd*1.0);
					pwm_out_val[3] = static_cast<float>(rcCmdsIn_.throttleCmd_nd*1.0);
				}else{
			  	// pwm_out_val[0] = max(PWM_CMD_MIN_THRESHOLD, min(PWM_CMD_MAX_THRESHOLD, ExtY_fcsModel_T_.actuatorsCmds[0]*RPM_TO_PWM_SCALE + PWM_MIN_THRESHOLD));
			  	// pwm_out_val[1] = max(PWM_CMD_MIN_THRESHOLD, min(PWM_CMD_MAX_THRESHOLD, ExtY_fcsModel_T_.actuatorsCmds[3]*RPM_TO_PWM_SCALE + PWM_MIN_THRESHOLD));
			  	// pwm_out_val[2] = max(PWM_CMD_MIN_THRESHOLD, min(PWM_CMD_MAX_THRESHOLD, ExtY_fcsModel_T_.actuatorsCmds[1]*RPM_TO_PWM_SCALE + PWM_MIN_THRESHOLD));
			  	// pwm_out_val[3] = max(PWM_CMD_MIN_THRESHOLD, min(PWM_CMD_MAX_THRESHOLD, ExtY_fcsModel_T_.actuatorsCmds[2]*RPM_TO_PWM_SCALE + PWM_MIN_THRESHOLD));
			  	pwm_out_val[0] = ExtY_fcsModel_T_.actuatorsPwmCmds[0];
			  	pwm_out_val[1] = ExtY_fcsModel_T_.actuatorsPwmCmds[3];
			  	pwm_out_val[2] = ExtY_fcsModel_T_.actuatorsPwmCmds[1];
			  	pwm_out_val[3] = ExtY_fcsModel_T_.actuatorsPwmCmds[2];
				}
		}else{
				pwm_out_val[0] = static_cast<float>(PWM_MIN_THRESHOLD*1.0);
		  	pwm_out_val[1] = static_cast<float>(PWM_MIN_THRESHOLD*1.0);
		  	pwm_out_val[2] = static_cast<float>(PWM_MIN_THRESHOLD*1.0);
		  	pwm_out_val[3] = static_cast<float>(PWM_MIN_THRESHOLD*1.0);
		}

	  pwm_writer.SetPwmDutyCyle(pwm_out_val);


    // Get the stop time and compute the duration
    loop_end = std::chrono::high_resolution_clock::now();

    duration_count = chrono::duration_cast<chrono::microseconds>(loop_end - loop_start).count();
    while(duration_count < dt_count){
    	duration_count = chrono::duration_cast<chrono::microseconds>(std::chrono::high_resolution_clock::now() - loop_start).count();
    }

   	if(sigint_flag == 1)
   		break;
	}
	return 0;
}
