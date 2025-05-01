#ifndef SLOWLOOP_H
#define SLOWLOOP_H

#include<Navio/Navio2/ADC_Navio2.h>
#include<write_utils.h>

#include<iostream>


#include<atomic>
#include<mutex>
#include<thread>
#include<pthread.h>
#include<unistd.h>

using namespace std;

class SlowLoopTasks{
	public:
		SlowLoopTasks();
		~SlowLoopTasks();

		void Start();
		void RunSlowLoop();
		void GetSlowLoopTasksData(WriteHelper::SlowLoopTasksData &slow_loop_tasks_data);

	private:
		// Thread to write file
		thread slow_loop_tasks_thread_;

		// To set CPU affinity
		cpu_set_t cpuset_;
		sched_param sch_;
		int policy_;

		// For ADC
		std::unique_ptr <ADC> slow_loop_task_adc_;
		float volt_v_, current_amp_;

		// Variable to indicate slow tasks to stop
	    atomic<bool> stop_slow_loop_tasks_thread_{false};

	    // Mutex to guard resource access between threads while running Slow Tasks Loop
	    mutex slow_loop_tasks_mutex_;
};


#endif