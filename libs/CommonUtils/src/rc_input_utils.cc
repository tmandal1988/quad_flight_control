#include "rc_input_utils.h"

//Constructor
RcInputHelper::RcInputHelper(uint8_t num_channels):
num_channels_(num_channels),
rc_periods_{new int[num_channels_]}{
	rc_input_ptr_ = unique_ptr<RCInput>{ new RCInput_Navio2() };
	for(size_t idx = 0; idx < num_channels_; idx++)
		rc_periods_[idx] = 0;
	// Variable to indicate RcInputReadLoop() to stop
	stop_rc_input_read_loop_.store(false);
}

void RcInputHelper::InitializeRcInput(){
	rc_input_ptr_->initialize();
}

void RcInputHelper::RcInputReadLoop(){
	// Keep reading the GPS data forever
	while(1){
		for(int idx = 0; idx < num_channels_; idx++){
            int temp_period = rc_input_ptr_->read(idx);
            {
            	unique_lock<mutex> rc_input_data_lock(rc_input_mutex_);
            	rc_periods_[idx] = temp_period;
            }                
            usleep(10000);
        }

		if(stop_rc_input_read_loop_.load()){
			break;
		}
		// 5 Hz read loop
		usleep(200000);
    }// While loop
}// RcInputReadLoop function

void RcInputHelper::CreateRcInputReadThread(){
	rc_input_thread_ = thread(&RcInputHelper::RcInputReadLoop, this);

    CPU_ZERO(&cpuset_);
    CPU_SET(2, &cpuset_);

	pthread_getschedparam(rc_input_thread_.native_handle(), &policy_, &sch_);
	sch_.sched_priority = 10;
	pthread_setschedparam(rc_input_thread_.native_handle(), SCHED_FIFO, &sch_);
	int rc = pthread_setaffinity_np(rc_input_thread_.native_handle(),
                                    sizeof(cpuset_), &cpuset_);    
	if (rc != 0) {
      std::cerr << "Error calling pthread_setaffinity_np on RC Input Read Thread: " << rc << "\n";
    }
}

int* RcInputHelper::GetRcPeriods(){
	return rc_periods_.get();
}

RcInputHelper::~RcInputHelper(){
	// Making sure that RC Input reading loop is stopped
	stop_rc_input_read_loop_.store(true);
	if(rc_input_thread_.joinable())
		rc_input_thread_.join();
}