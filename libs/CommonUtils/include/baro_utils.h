#ifndef BAROUTILS_H
#define BAROUTILS_H

#include<Navio/Common/MS5611.h>
#include<Navio/Common/Util.h>
#include<unistd.h>
#include<atomic>
#include<mutex>
#include<thread>
#include<pthread.h>
#include<iostream>

using namespace std;

class BaroHelper{
	public:
		// Constructors
		BaroHelper();

		// Destructor
		~BaroHelper();

		void StartBaroReader(int cpu_to_use, int32_t priority);
		void GetBaroPressAndTemp(float baro_data[]);

		bool is_baro_ready_;

	private:
		// Create a thread
		thread baro_reader_thread_;
		// To set CPU affinity
		cpu_set_t cpuset_;
		sched_param sch_;
		int policy_;

		// Mutex to guard resource access between threads while running GpsReadLoop()
	    mutex baro_data_mutex_;

		MS5611 barometer_;
		void BaroReadLoop();

		// Variable to indicate Baro Thread to stop
	    atomic<bool> stop_baro_read_thread_;
};


#endif