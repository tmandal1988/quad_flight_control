#include "slow_loop_tasks.h"


SlowLoopTasks::SlowLoopTasks(){
	slow_loop_task_adc_ = std::unique_ptr <ADC>{ new ADC_Navio2() };
	slow_loop_task_adc_->initialize();
}


void SlowLoopTasks::Start(){
	slow_loop_tasks_thread_ =  thread(&SlowLoopTasks::RunSlowLoop, this);

	pthread_getschedparam(slow_loop_tasks_thread_.native_handle(), &policy_, &sch_);
	sch_.sched_priority = 8;
    pthread_setschedparam(slow_loop_tasks_thread_.native_handle(), SCHED_FIFO, &sch_);
    CPU_ZERO(&cpuset_);
    CPU_SET(1, &cpuset_);

    int rc = pthread_setaffinity_np(slow_loop_tasks_thread_.native_handle(),
                                    sizeof(cpuset_), &cpuset_);    
	if (rc != 0) {
      std::cerr << "Error calling pthread_setaffinity_np on slow loop tasks thread: " << rc << "\n";
    }
}

void SlowLoopTasks::GetSlowLoopTasksData(WriteHelper::SlowLoopTasksData &slow_loop_tasks_data){
	{
		unique_lock<mutex> slow_loop_tasks_data_lock(slow_loop_tasks_mutex_);
		slow_loop_tasks_data.volt_v_ = volt_v_;
		slow_loop_tasks_data.current_amp_ = current_amp_;
	}

}

void SlowLoopTasks::RunSlowLoop(){	
	while(1){
    	volt_v_ = slow_loop_task_adc_->read(2)/1000.0;
    	usleep(2000);
    	current_amp_ = slow_loop_task_adc_->read(3)/1000.0;
    	if(stop_slow_loop_tasks_thread_.load()){
				break;
		}

		usleep(48000);
	}
}

SlowLoopTasks::~SlowLoopTasks(){
	// Making sure that slow task loop is stopped
	stop_slow_loop_tasks_thread_.store(true);
	sleep(1);
	if(slow_loop_tasks_thread_.joinable())
		slow_loop_tasks_thread_.join();
}

