#ifndef RCINPUTUTILS_H
#define RCINPUTUTILS_H

#include "Navio/Navio2/RCInput_Navio2.h"
#include<atomic>
#include<mutex>
#include<thread>
#include<pthread.h>
#include<unistd.h>
#include<iostream>

using namespace std;

class RcInputHelper{
	#define READ_FAILED -1

	public:
		// Constructors
		RcInputHelper(uint8_t num_channels);

		// Destructor
		~RcInputHelper();

		// To initialize the RC Input Channels
		void InitializeRcInput();

		// Start the GPS running thread
		void CreateRcInputReadThread();

		// Getter function
		int* GetRcPeriods();

	private:
		// Create a thread
		thread rc_input_thread_;
		// To set CPU affinity
		cpu_set_t cpuset_;
		sched_param sch_;
		int policy_;

		unique_ptr<RCInput> rc_input_ptr_;

		uint8_t num_channels_;
		// Array to store rc_periods
        std::unique_ptr<int[]> rc_periods_;
		
	    // Mutex to guard resource access between threads while running GpsReadLoop()
	    mutex rc_input_mutex_;

	    // Variable to indicate RcInputReadLoop() to stop
	    atomic<bool> stop_rc_input_read_loop_;

	    // Useful function that can be used to read RC Input data in a loop at a configured rate
		// This function can be passed to a thread to update RcInput Data
		void RcInputReadLoop();

		// Stops the above loop
		void StopRcInputReadLoop(){
			// Stop the RC Input reading loop
			stop_rc_input_read_loop_.store(true);
		}
};

#endif