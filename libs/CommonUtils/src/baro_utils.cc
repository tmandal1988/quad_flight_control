#include "baro_utils.h"

//Constructor
BaroHelper::BaroHelper():
is_baro_ready_(false){	
	// Variable to indicate Baro Thread to stop
	stop_baro_read_thread_.store(false);

	barometer_.initialize();
}

void BaroHelper::StartBaroReader(int cpu_to_use, int32_t priority){
	baro_reader_thread_ = thread(&BaroHelper::BaroReadLoop, this);

    CPU_ZERO(&cpuset_);
    CPU_SET(cpu_to_use, &cpuset_);

	pthread_getschedparam(baro_reader_thread_.native_handle(), &policy_, &sch_);
	sch_.sched_priority = priority;
	pthread_setschedparam(baro_reader_thread_.native_handle(), SCHED_FIFO, &sch_);
	int rc = pthread_setaffinity_np(baro_reader_thread_.native_handle(),
                                    sizeof(cpuset_), &cpuset_);    
	if (rc != 0) {
      std::cerr << "Error calling pthread_setaffinity_np on Baro Reader Thread: " << rc << "\n";
    }
}

void BaroHelper::BaroReadLoop(){
	{
        unique_lock<mutex> baro_data_lock(baro_data_mutex_);
		is_baro_ready_ = true;
	}
	// Keep reading the Baro data forever
	while(1){
		barometer_.refreshPressure();
		usleep(5000); // Waiting for pressure data ready

		barometer_.refreshTemperature();
        usleep(5000); // Waiting for temperature data ready
        barometer_.readTemperature();

        {
         	unique_lock<mutex> baro_data_lock(baro_data_mutex_);
        	barometer_.calculatePressureAndTemperature();
        }

        // 25 Hz read loop
		usleep(40000);

		if(stop_baro_read_thread_.load()){
				break;
		}
	}
}

void BaroHelper::GetBaroPressAndTemp(float baro_data[]){
	baro_data[0] = barometer_.getPressure();
	baro_data[1] = barometer_.getTemperature();

}

BaroHelper::~BaroHelper(){
	// Making sure that Baro read thread is stopped
	stop_baro_read_thread_.store(true);
	if(baro_reader_thread_.joinable())
		baro_reader_thread_.join();
}