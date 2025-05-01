// Autocoded FCS And State Estimator Libraries
#include <stateEstimatorEskfAutocode.h>
#include <fcsModel.h>
#include <fcs_params.h>
// Navio 2 Utilities
#include "Navio/Common/MPU9250.h"
#include "Navio/Navio2/LSM9DS1.h"
#include "Navio/Common/Util.h"
#include <gps_utils.h>
#include <write_utils.h>
#include <slow_loop_tasks.h>
#include <imu_utils.h>
#include <lidar_utils.h>
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
static stateEstimatorAutocode stateEstimator_Obj;          // Instance of State Esimator model class

// To catch SIGINT
volatile sig_atomic_t sigint_flag = 0;

void sigint_handler(int sig){ // can be called asynchronously
  sigint_flag = 1; // set flag
}

int main(int argc, char *argv[]){
	system("sudo echo -1 >/proc/sys/kernel/sched_rt_runtime_us");

	// GPS reader object
	GpsHelper gps_reader;

	// Data logger object
	WriteHelper data_writer("data_file.dat");

	// IMU reader object
	ImuHelper imu_reader("mpu");
	/*Initializes the hardware necessary to read data from MPU9250
	* This instance is run without any notch filter in the imu utility
	* as all filtering is handled in state estimator block now
	*/
	imu_reader.InitializeImu();

	// Baro reader object
	BaroHelper baro_reader;
	baro_reader.StartRawBaroReader(1, 20);

	// Lidar reader object
	LidarHelper lidar_reader;
	bool lidar_init_status = lidar_reader.InitializeLidar();
	bool is_lidar_valid = false;
	lidar_reader.CreateLidarThread();
	float lidar_range_m = 0;

	RcInputHelper rc_reader(8);
	rc_reader.InitializeRcInput();

	PwmOutputHelper pwm_writer(4);
	pwm_writer.InitializePwmOutput();

	vector<float> pwm_out_val(4, 0);

	// All the monitoring tasks such as monitoring battery voltage and current
	SlowLoopTasks slow_loop_tasks;
	slow_loop_tasks.Start();
	WriteHelper::SlowLoopTasksData slow_loop_tasks_data;

	float ned_pos_and_vel_meas[6];
	bool gps_meas_indices[6];
	float gps_raw_pos_and_vel[6];
	
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

  /****************Variables to read data from the Baro*************************************/
  float press_pa = 101325.0f, temp_c = 15.0f; // Start at mean sea level
  /****************Variables to read data from the Baro*************************************/

	/****************SAMPLE TIME VARIABLES*************************************/
	//Sample time
	float dt_s = 0.004;
	//Sample Step
	size_t dt_count = 4000;
	/****************SAMPLE TIME VARIABLES*************************************/

	// NED to BODY DCM
	MatrixInv<float> c_ned2b;
	// NED to FEP DCM
	MatrixInv<float> c_ned2fep;

	// Initialize
	bool gps_init_status = gps_reader.InitializeGps(60);
	bool is_gps_valid = false;

	gps_reader.CreateRawGpsThread();
	data_writer.StartFileWriteThread2();
	rc_reader.CreateRcInputReadThread();

	double* vned_init = gps_reader.GetInitNedVel();

	// Variable to read imu data
	float imu_raw_data[9] = {0};

	//Variable to read rc input data
	int* rc_periods = new int[7];

	// Assign -1 to rc_periods to initialize
	for(size_t rc_idx = 0; rc_idx < 7; rc_idx++){
		rc_periods[rc_idx] = -1;
	}

	//##############################################################
	// fcsModel input variable
	fcsModel::ExtU_fcsModel_T *ExtU_fcsModel_T_ =  new fcsModel::ExtU_fcsModel_T;
	// stateEstimator input variable
	stateEstimatorAutocode::ExtU_stateEstimatorEskfAutoco_T *ExtU_stateEstimatorAutocode_T_ =  
						new stateEstimatorAutocode::ExtU_stateEstimatorEskfAutoco_T;
	
	// Initialize model
  fcsModel_Obj.initialize();
  stateEstimator_Obj.initialize();

  // RC Cmd Input Variables
  busRcInCmds rcCmdsIn_;
  // Sensor inuput to the model
  busStateEstimate stateEstimate_;
  // FCS Ctrl Params from fcs_params.h
  ExtU_fcsModel_T_->ctrlParams = AssignFcsCtrlParams();

  // State estimator input variables
  busImuData imuData_;                
  busMagData magData_;                
  busGpsData gpsData_;                
  busBaroData baroData_;

  // fcsModel output variable
  fcsModel::ExtY_fcsModel_T ExtY_fcsModel_T_;

  // stateEstimator output variable
  stateEstimatorAutocode::ExtY_stateEstimatorEskfAutoco_T ExtY_stateEstimatorAutocode_T_;

  // Mutex to guard resource access to fcs outputs
	mutex fcs_out_mutex;
  // {
  // 	unique_lock<mutex> fcs_out_lock(fcs_out_mutex);	
	// 	ExtY_fcsModel_T_ = fcsModel_Obj.getExternalOutputs();
	// }

	// Loop counter
	size_t loop_count = 0;
	size_t arm_loop_count = 0;
	//125 Hz flag
	bool one_twenty_five_hz_flag = false;
	// 50 Hz flag
	bool fifty_hz_flag = false;
	// 10 Hz flag
	bool ten_hz_flag = false;
	// 25 Hz flag
	bool twenty_five_hz_flag = false;
	// initialize the duration
	chrono::microseconds delta (dt_count); 
	auto duration = chrono::duration_cast<chrono::microseconds> (delta);
	auto duration_count = duration.count();

	// Loop timers
	chrono::steady_clock::time_point loop_start;
	chrono::steady_clock::time_point loop_end;

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

    	/* Check if the current loop count is a multiple of 210 which will give a 25hz loop as main loop
    	runs at 250Hz
    	*/
    	if(loop_count % 10 == 0){
    		twenty_five_hz_flag = true;
    	}else{
    		twenty_five_hz_flag = false;
    	}

    	/* Get loop start time
    	*/
    	loop_start = chrono::steady_clock::now();

    	// Read IMU data
	    bool is_mag_valid = imu_reader.getRawImuData(imu_raw_data);
	    // Assign IMU data to state estimator input
	    for(size_t imuIdx = 0; imuIdx < 3; imuIdx++){
	    	ExtU_stateEstimatorAutocode_T_->imuData.bodyAccels_mps2[imuIdx] = imu_raw_data[imuIdx];
	    	ExtU_stateEstimatorAutocode_T_->imuData.bodyRates_radps[imuIdx] = imu_raw_data[imuIdx + 3];
	    	ExtU_stateEstimatorAutocode_T_->magData.bodyMagVector_uT[imuIdx] = imu_raw_data[imuIdx + 6];
	    }

	    ExtU_stateEstimatorAutocode_T_->magData.isMagDataValid = is_mag_valid;	  

	    // if(gps_init_status){
		    // Read GPS data
	    is_gps_valid = gps_reader.GetGpsRawPosAndVel(gps_raw_pos_and_vel);
	    for(size_t gpsIdx = 0; gpsIdx < 3; gpsIdx++){
	    	ExtU_stateEstimatorAutocode_T_->gpsData.latLonAlt[gpsIdx] = gps_raw_pos_and_vel[gpsIdx];
	    	ExtU_stateEstimatorAutocode_T_->gpsData.nedVel_mps[gpsIdx] = gps_raw_pos_and_vel[gpsIdx + 3];
	    }
	  	// }

	    ExtU_stateEstimatorAutocode_T_->gpsData.isGpsDataValid = is_gps_valid;
	    ExtU_stateEstimatorAutocode_T_->gpsData.isGpsInitialized = gps_init_status;

      bool is_baro_valid = baro_reader.GetRawPressAndTemp(press_pa, temp_c);
			ExtU_stateEstimatorAutocode_T_->baroData.pressure_pa = press_pa;
			ExtU_stateEstimatorAutocode_T_->baroData.temp_c = temp_c;
			ExtU_stateEstimatorAutocode_T_->baroData.isBaroDataValid = is_baro_valid;

			ExtU_stateEstimatorAutocode_T_->lidarData.isLidarInitialized = lidar_init_status;
			// is_lidar_valid = lidar_reader.GetLidarRange(lidar_range_m);
			if(lidar_init_status){
     			is_lidar_valid = lidar_reader.GetLidarRange(lidar_range_m);
     			ExtU_stateEstimatorAutocode_T_->lidarData.range_m = lidar_range_m;
     			ExtU_stateEstimatorAutocode_T_->lidarData.isLidarDataValid = is_lidar_valid;
     	}

			/****************** RUN STATE ESTIMATOR ****************************/
			stateEstimator_Obj.setExternalInputs(ExtU_stateEstimatorAutocode_T_);
			stateEstimator_Obj.step();
  		ExtY_stateEstimatorAutocode_T_ = stateEstimator_Obj.getExternalOutputs();
  		cout.flush();
  		if(ExtY_stateEstimatorAutocode_T_.stateEstimatorDebug.smMode == enumStateEstimateMode::INITIALIZE){
				cout<<"EKF INITIALIZATION PROGRESS: "<< ExtY_stateEstimatorAutocode_T_.stateEstimatorDebug.stateEstInitPct<<" %\r";
  		}

			// if(ten_hz_flag){
			// 	if(ExtY_stateEstimatorAutocode_T_.stateEstimatorDebug.smMode != enumStateEstimateMode::INITIALIZE){
			// 		printf("Roll : %.2f, Pitch: %.2f, Yaw: %.2f, SM Mode: %d\n", ExtY_stateEstimatorAutocode_T_.eulAng_rad[0]*RAD2DEG,
			// 			ExtY_stateEstimatorAutocode_T_.eulAng_rad[1]*RAD2DEG, ExtY_stateEstimatorAutocode_T_.eulAng_rad[2]*RAD2DEG, 
			// 			ExtY_stateEstimatorAutocode_T_.stateEstimatorDebug.smMode);
			// 	}
  		// }
  		/****************** RUN STATE ESTIMATOR ****************************/

  		//######################################## Set FCS Inputs ########################################
  		if(ExtY_stateEstimatorAutocode_T_.stateEstimatorDebug.smMode != enumStateEstimateMode::INITIALIZE){
  			/****************** Assign Values To StateEstimate field of the flight controller*/
  			for(size_t idx = 0; idx < 3; idx++){
      		stateEstimate_.attitude_rad[idx] = ExtY_stateEstimatorAutocode_T_.eulAng_rad[idx];
      		stateEstimate_.bodyAngRates_radps[idx] = ( imu_raw_data[idx + 3] - ExtY_stateEstimatorAutocode_T_.states[idx + 10] );
      		stateEstimate_.bodyAccels_mps2[idx] = ExtY_stateEstimatorAutocode_T_.bodyAccels_mps2[idx];
      		stateEstimate_.nedPos_m[idx] = ExtY_stateEstimatorAutocode_T_.states[idx + 4];
      		stateEstimate_.nedVel_mps[idx] = ExtY_stateEstimatorAutocode_T_.states[idx + 7];
      	}
  		
      	for(size_t dIdx = 0; dIdx < 9; dIdx++){
	  			stateEstimate_.ned2BodyDcm_nd[dIdx] = ExtY_stateEstimatorAutocode_T_.dcmNedToBody[dIdx];
	  			stateEstimate_.ned2FepDcm_nd[dIdx] = ExtY_stateEstimatorAutocode_T_.dcmNedToFep[dIdx];
	  		}
	  		stateEstimate_.aglEst_m = -ExtY_stateEstimatorAutocode_T_.states[6];
	  		stateEstimate_.climbRateEst_mps = -ExtY_stateEstimatorAutocode_T_.states[9];
	  		stateEstimate_.pressure_mbar = press_pa/100.0f;
	  		stateEstimate_.temp_c = temp_c;
	  		stateEstimate_.geodeticPos.lat_rad = gps_raw_pos_and_vel[0];
	  		stateEstimate_.geodeticPos.lon_rad = gps_raw_pos_and_vel[1];
	  		stateEstimate_.geodeticPos.alt_m = gps_raw_pos_and_vel[2];  

	      ExtU_fcsModel_T_->stateEstimate = stateEstimate_;
	      /****************** Assign Values To StateEstimate field of the flight controller*/

	      /****************** Assign Values To rcCmdsIn field of the flight controller*/
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
	      /****************** Assign Values To rcCmdsIn field of the flight controller*/
	      fcsModel_Obj.setExternalInputs(ExtU_fcsModel_T_);
	      //######################################## Set FCS Inputs ########################################

	      // Step the FCS model
	      {
	      	unique_lock<mutex> fcs_out_lock(fcs_out_mutex);
	  			fcsModel_Obj.step();
	  			ExtY_fcsModel_T_ = fcsModel_Obj.getExternalOutputs();
	  		}

	     if(fifty_hz_flag){
	     		slow_loop_tasks.GetSlowLoopTasksData(slow_loop_tasks_data);
	     	}
     	}

     	// if(ten_hz_flag){
			// 		// printf("SM Mode: %d, N Pos Cmd: %g, E Pos Cmd: %g, N Vel Cmd: %g, E Vel Cmd: %g\n", ExtY_stateEstimatorAutocode_T_.stateEstimatorDebug.smMode,
			// 		// 		ExtY_fcsModel_T_.fcsDebug.outerLoopCtrlDebug.posCtrlDebug.cmd[0], 
			// 		// 		ExtY_fcsModel_T_.fcsDebug.outerLoopCtrlDebug.posCtrlDebug.cmd[1],
			// 		// 		ExtY_fcsModel_T_.fcsDebug.outerLoopCtrlDebug.velCtrlDebug.cmd[0], 
			// 		// 		ExtY_fcsModel_T_.fcsDebug.outerLoopCtrlDebug.velCtrlDebug.cmd[1]);
     	// 		if(is_lidar_valid){
     	// 			printf("Lidar Valid and Lidar Range: %.2f\n", lidar_range_m);
     	// 		}else{
     	// 			printf("Lidar Invalid\n");
     	// 		}
			// 	}
     	std::array<double, 23> ekf_states_;
			for (size_t i = 0; i < 20; ++i) {
		    	ekf_states_[i] = static_cast<double>(ExtY_stateEstimatorAutocode_T_.states[i]);
			}
			ekf_states_[20] = 0.0f;
			ekf_states_[21] = 0.0f;
			ekf_states_[22] = 0.0f;
     data_writer.UpdateDataBuffer2(duration_count, loop_count, imu_raw_data, is_mag_valid, gps_raw_pos_and_vel,
     																is_gps_valid, press_pa, temp_c, is_baro_valid, lidar_range_m, is_lidar_valid,
     																ekf_states_, slow_loop_tasks_data, rc_periods, ExtY_fcsModel_T_);
	   
	    // if (fifty_hz_flag){
	    // 	printf("Roll [deg]: %+7.3f, Pitch[deg]: %+7.3f, Yaw[deg]: %+7.3f\n", current_state(0)*RAD2DEG, current_state(1)*RAD2DEG, current_state(2)*RAD2DEG);
	  // //   	// printf("Pos N [m]: %+7.3f, Pos E [m]: %+7.3f, Pos D[m]: %+7.3f\n", current_state(6), current_state(7), current_state(8));
	  // //   	// printf("Vel N [m]: %+7.3f, Vel E [m]: %+7.3f, Vel D[m]: %+7.3f\n", current_state(9), current_state(10), current_state(11));
	    	// printf("Throttle: %d, Roll: %d, Pitch: %d, Yaw: %d, Sw1: %d, Sw2: %d, Sw3: %d, State: %d, Flight Mode: %d\n", ExtU_fcsModel_T_->rcCmdsIn.throttleCmd_nd, ExtU_fcsModel_T_->rcCmdsIn.joystickXCmd_nd, 
	    		// ExtU_fcsModel_T_->rcCmdsIn.joystickYCmd_nd, ExtU_fcsModel_T_->rcCmdsIn.joystickZCmd_nd, ExtU_fcsModel_T_->rcCmdsIn.rcSwitch1_nd, ExtU_fcsModel_T_->rcCmdsIn.rcSwitch2_nd, 
	    		// ExtU_fcsModel_T_->rcCmdsIn.rcSwitch3_nd, ExtY_fcsModel_T_.fcsDebug.state, ExtY_fcsModel_T_.fcsDebug.flightMode);
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
    loop_end = std::chrono::steady_clock::now();

    duration_count = chrono::duration_cast<chrono::microseconds>(loop_end - loop_start).count();
    while(duration_count < dt_count){
    	duration_count = chrono::duration_cast<chrono::microseconds>(std::chrono::steady_clock::now() - loop_start).count();
    }

   	if(sigint_flag == 1)
   		break;
	}
	return 0;
}
