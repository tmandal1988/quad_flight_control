#ifndef WRITEUTILS_H
#define WRITEUTILS_H

#include<Matrix/matrix_inv_class.h>
#include<fcs_out_data.h>

#include <fstream>
#include<atomic>
#include<mutex>
#include<thread>
#include<pthread.h>
#include<unistd.h> 

using namespace std;

class WriteHelper{
	#define MAX_BUFF_SIZE 1024
	public:
		WriteHelper(const string& file_name);
		~WriteHelper();

		// Thread to loop and data to file
		void StopWriteLoop();
		void UpdateDataBuffer(long long dt_ms, size_t count, float imu_data[9], MatrixInv<float> sensor_meas, MatrixInv<float> ekf_current_state, 
			bool gps_valid_flag[6], int* const& rc_periods, const FcsOutput &fcs_output);
		void StartFileWriteThread();
	private:
		// Thread to write file
		thread write_thread_;

		// To set CPU affinity
		cpu_set_t cpuset_;
		sched_param sch_;
		int policy_;

		// File object
		ofstream file_object_;

		// Mutex to guard resource access to write buffer
	    mutex data_mutex_;

		// Variable to store file name
		string file_name_;

		// Create two buffers to store the data to save before writing
		DataFields data_to_save1_[MAX_BUFF_SIZE];
		DataFields data_to_save2_[MAX_BUFF_SIZE];

		// Variable to indicate when the data buffers are full of data
	    atomic<bool> is_data_buff1_full_;
	    atomic<bool> is_data_buff2_full_;

	    // Variable to indicate wrte loop to stop
	    atomic<bool> stop_data_write_loop_;

	    size_t data_buff_idx1_;
	    size_t data_buff_idx2_;

	    // write to file
	    void WriteToFileLoop();

};


#endif