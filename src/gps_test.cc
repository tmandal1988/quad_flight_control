#include "Navio/Common/UbloxDriver.h"

#include <signal.h>
#include <chrono>

// To catch SIGINT
volatile sig_atomic_t sigint_flag = 0;

void sigint_handler(int sig){ // can be called asynstd::chronously
  sigint_flag = 1; // set flag
}

int main(int argc, char *argv[]){
	UbxDriver gps_m8n;

	int connection_status = gps_m8n.TestConnection();
	printf("Connection Status:%d\n", connection_status);

	UbxDriver::NavPvtData nav_pvt_data;

	gps_m8n.GetNavPvtData(nav_pvt_data);

	return 0;

	// initialize the duration
	size_t dt_count = 4000;
	std::chrono::microseconds delta (dt_count); 
	auto duration = std::chrono::duration_cast<std::chrono::microseconds> (delta);
	auto duration_count = duration.count();

	// Loop timers
	std::chrono::steady_clock::time_point loop_start;
	std::chrono::steady_clock::time_point loop_end;

	while(1){

	// 	/* Get loop start time
  //   */
  //   loop_start = std::chrono::steady_clock::now();

	// 	// gps_itow_ms = gps_reader.GetGpsNedPosVelAndTime(ned_pos_and_vel_meas, gps_meas_indices, pvt_times);
	// 	// if(gps_meas_indices[0] == 1){
	// 	// 	printf("GPS iTOW: %ld, Year: %d, Month: %d, Hour: %d, Min: %d, Sec: %d\n", gps_itow_ms, pvt_times[0], pvt_times[1],
	// 	// 				  pvt_times[2], pvt_times[3], pvt_times[4], pvt_times[5]);
	// 	// }

	// 	// Get the stop time and compute the duration
  //   // loop_end = std::std::chrono::steady_clock::now();

  //   // duration_count = std::chrono::duration_cast<std::chrono::microseconds>(loop_end - loop_start).count();
  //   // while(duration_count < dt_count){
  //   // 	duration_count = std::chrono::duration_cast<std::chrono::microseconds>(std::std::chrono::steady_clock::now() - loop_start).count();
  //   // }

    usleep(100);

		if(sigint_flag == 1)
   			break;
	}

	return 0;
}