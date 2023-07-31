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
// Standard C++ Libraries file, time and memory
#include <memory>
#include <pwm_output_utils.h>

#include<iostream>
#include <chrono>
// Standard C++ Libraries for multi-threading
#include <thread>
#include <pthread.h>
#include <mutex>
#include <atomic>
#include <signal.h>
#include <sys/mman.h>
#include <iterator>

int main(int argc, char *argv[]){
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

    RcInputHelper rc_reader(8);
	rc_reader.InitializeRcInput();
	rc_reader.CreateRcInputReadThread();

	PwmOutputHelper pwm_writer(4);
	pwm_writer.InitializePwmOutput();

	//Variable to read rc input data
	int* rc_periods;

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

		rc_periods = rc_reader.GetRcPeriods();
	  	int throttleCmd_nd = rc_periods[2];
	  	printf("Throttle: %d\n", throttleCmd_nd);

	  	vector<float> pwm_out_val(4, throttleCmd_nd*1.0);
	  	pwm_writer.SetPwmDutyCyle(pwm_out_val);
	  	// joystickYCmd_nd = rc_periods[1];
	  	// joystickXCmd_nd = rc_periods[0];
	  	// joystickZCmd_nd = rc_periods[3];
	  	// rcSwitch1_nd = rc_periods[4];
	  	// rcSwitch2_nd = rc_periods[5];
	  	// rcSwitch3_nd = rc_periods[6];
  	}

	return 0;
}