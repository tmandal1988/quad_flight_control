#ifndef WRITEUTILS_H
#define WRITEUTILS_H

#include<Matrix/matrix_inv_class.h>

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
			bool gps_valid_flag[6]);
		void StartFileWriteThread();
	private:
		// Struct defining what data to save
		struct data_fields{
			// Time it took for the EKF to run
			float dt_s;
			// GPS Time if available
			float time_of_week_ms;
			// Raw Sensor Data (Mag is calibrated and normalized)
			float imu_data[9];
			// GPS NED position and velocity if available;
			float ned_pos_m[3];
			float ned_vel_mps[3];
			// Current EKF State
			float ekf_current_state[16];
			// GPS Valid Flag
			float gps_valid_flag[6];
			// Debug Data
			// float ekf_state_jacobian[16*16];
			// // computed meas
			// float computed_meas[9];
			// // Debug Data
			// float ekf_meas_jacobian[9*16];
			// //Debug Data
			// float ekf_computed_meas[9];

		};

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
		data_fields data_to_save1_[MAX_BUFF_SIZE];
		data_fields data_to_save2_[MAX_BUFF_SIZE];

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