// EKF related headers
#include <ekf_15dof_class.h>
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

	float ned_pos_and_vel_meas[6];
	bool gps_meas_indices[6];
	
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

	gps_reader.InitializeGps();
	gps_reader.CreateGpsThread();
	data_writer.StartFileWriteThread();
	rc_reader.CreateRcInputReadThread();

	double* vned_init = gps_reader.GetInitNedVel();	

	for(size_t i_idx = 0; i_idx < 3; i_idx++){
		initial_state(i_idx) = init_att[i_idx];
		initial_state(i_idx + 9) = vned_init[i_idx];
	}

	//write the initial state
	float zero_array[9] = {0};
	data_writer.UpdateDataBuffer(0, 0, zero_array, sensor_meas, initial_state, gps_meas_indices);

	// Create an 15 state EKF object
	Ekf15Dof<float> imu_gps_ekf(0.004, initial_state, process_noise_q, meas_noise_r, initial_covariance_p);

	// Make sure EKF always uses Magnetometer in the update step
	meas_indices[0] = true;

	// Variable to read imu data
	float* imu_data;

	//Variable to read rc input data
	int* rc_periods;

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

        if(remainder(loop_count, 5) == 0){
        	data_writer.UpdateDataBuffer(duration.count(), loop_count, imu_data, sensor_meas, current_state, gps_meas_indices);
	    }

	    if (remainder(loop_count, 5) == 0){
	    	rc_periods = rc_reader.GetRcPeriods();
	    	// printf("Roll [deg]: %+7.3f, Pitch[deg]: %+7.3f, Yaw[deg]: %+7.3f\n", current_state(0)*RAD2DEG, current_state(1)*RAD2DEG, current_state(2)*RAD2DEG);
	    	// printf("Pos N [m]: %+7.3f, Pos E [m]: %+7.3f, Pos D[m]: %+7.3f\n", current_state(6), current_state(7), current_state(8));
	    	// printf("Vel N [m]: %+7.3f, Vel E [m]: %+7.3f, Vel D[m]: %+7.3f\n", current_state(9), current_state(10), current_state(11));
	    	printf("%d, %d, %d, %d, %d, %d, %d, %d\n", rc_periods[0], rc_periods[1], rc_periods[2], rc_periods[3], rc_periods[4], rc_periods[5], rc_periods[6], rc_periods[7]);
	    	printf("############################################\n");
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