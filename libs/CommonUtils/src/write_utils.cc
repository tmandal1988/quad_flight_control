#include "write_utils.h"

WriteHelper::WriteHelper(const string& file_name):file_name_(file_name){
	is_data_buff1_full_.store(false);
	is_data_buff2_full_.store(false);

	data_buff_idx1_ = 0;
	data_buff_idx2_ = 0;

	file_object_.open(file_name_, ios::out | ios::binary);

	stop_data_write_loop_.store(false);
}

void WriteHelper::StopWriteLoop(){
	stop_data_write_loop_.store(true);
}

void WriteHelper::WriteToFileLoop(){
	const size_t bufsize = 1024 * 1024;
	unique_ptr<char[]> buf(new char[bufsize]);
	file_object_.rdbuf()->pubsetbuf(buf.get(), bufsize);
	while(1){
		if(is_data_buff1_full_.load()){
			file_object_.write((char*)data_to_save1_, 1024*sizeof(data_fields));
			is_data_buff1_full_.store(false);
		}else if(is_data_buff2_full_.load()){
			file_object_.write((char*)data_to_save2_, 1024*sizeof(data_fields));
			is_data_buff2_full_.store(false);
		}else{
			sleep(1);
		}
		if(stop_data_write_loop_.load())
			break;
	}// While loop

}

WriteHelper::~WriteHelper(){
	// Making sure that write loop is stopped
	stop_data_write_loop_.store(true);
	sleep(1);
	if(file_object_.is_open()){
		file_object_.close();
	}
	if(write_thread_.joinable())
		write_thread_.join();
}

void WriteHelper::UpdateDataBuffer(long long dt_ms, size_t count, float imu_data[9], 
	MatrixInv<float> sensor_meas, MatrixInv<float> ekf_current_state, bool gps_valid_flag[6]){
	if(data_buff_idx2_ == 0 && data_buff_idx1_ != MAX_BUFF_SIZE){
		{
			unique_lock<mutex> save_data_lock(data_mutex_);
			data_to_save1_[data_buff_idx1_].dt_s = static_cast< float >(dt_ms);
			data_to_save1_[data_buff_idx1_].time_of_week_ms = static_cast< float >(count);
			for(size_t imu_idx = 0; imu_idx < 3; imu_idx++){
				data_to_save1_[data_buff_idx1_].imu_data[imu_idx] = imu_data[imu_idx];
				data_to_save1_[data_buff_idx1_].imu_data[imu_idx + 3] = imu_data[imu_idx + 3];
				data_to_save1_[data_buff_idx1_].imu_data[imu_idx + 6] = imu_data[imu_idx + 6];
				data_to_save1_[data_buff_idx1_].ned_pos_m[imu_idx] = sensor_meas(imu_idx + 3);
				data_to_save1_[data_buff_idx1_].ned_vel_mps[imu_idx] = sensor_meas(imu_idx + 6);
			}

			for(size_t d_idx = 0; d_idx < ekf_current_state.get_nrows(); d_idx++){
				data_to_save1_[data_buff_idx1_].ekf_current_state[d_idx] = ekf_current_state(d_idx);
			}

			for(size_t g_idx = 0; g_idx < 6; g_idx++){
				data_to_save1_[data_buff_idx1_].gps_valid_flag[g_idx] = gps_valid_flag[g_idx];
			}

			/* To save ekf debugging messsages
			// size_t ekfdb_idx = 0;
			// for(size_t dbr_idx = 0; dbr_idx < ekf_state_jacobian.get_nrows(); dbr_idx++){
			// 	for(size_t dbc_idx = 0; dbc_idx < ekf_state_jacobian.get_ncols(); dbc_idx++){
			// 		data_to_save1_[data_buff_idx1_].ekf_state_jacobian[ekfdb_idx] = ekf_state_jacobian(dbr_idx, dbc_idx);
			// 		ekfdb_idx++;
			// 	}
			// }

			// for(size_t cm_idx = 0; cm_idx < 9; cm_idx++){
			// 	data_to_save1_[data_buff_idx1_].computed_meas[cm_idx] = computed_meas(cm_idx);
			// }

			// size_t ekfmj_idx = 0;
			// for(size_t dbr_idx = 0; dbr_idx < ekf_meas_jacobian.get_nrows(); dbr_idx++){
			// 	for(size_t dbc_idx = 0; dbc_idx < ekf_meas_jacobian.get_ncols(); dbc_idx++){
			// 		data_to_save1_[data_buff_idx1_].ekf_meas_jacobian[ekfmj_idx] = ekf_meas_jacobian(dbr_idx, dbc_idx);
			// 		ekfmj_idx++;
			// 	}
			// }

			// for(size_t tt_idx = 0; tt_idx < ekf_computed_meas.get_nrows(); tt_idx++){
			// 	data_to_save1_[data_buff_idx1_].ekf_computed_meas[tt_idx] = ekf_computed_meas(tt_idx);
			// }
			*/

			data_buff_idx1_++;
			if (data_buff_idx1_ == MAX_BUFF_SIZE){
				is_data_buff1_full_.store(true);
				data_buff_idx2_ = 0;
			}
		}
	}else if(data_buff_idx1_  == MAX_BUFF_SIZE){
		{
			unique_lock<mutex> save_data_lock(data_mutex_);
			data_to_save2_[data_buff_idx2_].dt_s = static_cast< float >(dt_ms);
			data_to_save2_[data_buff_idx2_].time_of_week_ms = static_cast< float >(count);
			for(size_t imu_idx = 0; imu_idx < 3; imu_idx++){
				data_to_save2_[data_buff_idx2_].imu_data[imu_idx] = imu_data[imu_idx];
				data_to_save2_[data_buff_idx2_].imu_data[imu_idx + 3] = imu_data[imu_idx + 3];
				data_to_save2_[data_buff_idx2_].imu_data[imu_idx + 6] = imu_data[imu_idx + 6];
				data_to_save2_[data_buff_idx2_].ned_pos_m[imu_idx] = sensor_meas(imu_idx + 3);
				data_to_save2_[data_buff_idx2_].ned_vel_mps[imu_idx] = sensor_meas(imu_idx + 6);
			}

			for(size_t d_idx = 0; d_idx < ekf_current_state.get_nrows(); d_idx++){
				data_to_save2_[data_buff_idx2_].ekf_current_state[d_idx] = ekf_current_state(d_idx);
			}

			for(size_t g_idx = 0; g_idx < 6; g_idx++){
				data_to_save2_[data_buff_idx2_].gps_valid_flag[g_idx] = gps_valid_flag[g_idx];
			}

			/* To save ekf debugging messsages
			// size_t ekfdb_idx = 0;
			// for(size_t dbr_idx = 0; dbr_idx < ekf_state_jacobian.get_nrows(); dbr_idx++){
			// 	for(size_t dbc_idx = 0; dbc_idx < ekf_state_jacobian.get_ncols(); dbc_idx++){
			// 		data_to_save2_[data_buff_idx2_].ekf_state_jacobian[ekfdb_idx] = ekf_state_jacobian(dbr_idx, dbc_idx);
			// 		ekfdb_idx++;
			// 	}
			// }

			// for(size_t cm_idx = 0; cm_idx < 9; cm_idx++){
			// 	data_to_save2_[data_buff_idx2_].computed_meas[cm_idx] = computed_meas(cm_idx);
			// }

			// size_t ekfmj_idx = 0;
			// for(size_t dbr_idx = 0; dbr_idx < ekf_meas_jacobian.get_nrows(); dbr_idx++){
			// 	for(size_t dbc_idx = 0; dbc_idx < ekf_meas_jacobian.get_ncols(); dbc_idx++){
			// 		data_to_save2_[data_buff_idx2_].ekf_meas_jacobian[ekfmj_idx] = ekf_meas_jacobian(dbr_idx, dbc_idx);
			// 		ekfmj_idx++;
			// 	}
			// }

			// for(size_t tt_idx = 0; tt_idx < ekf_computed_meas.get_nrows(); tt_idx++){
			// 	data_to_save2_[data_buff_idx2_].ekf_computed_meas[tt_idx] = ekf_computed_meas(tt_idx);
			// }
			*/

			data_buff_idx2_++;
			if (data_buff_idx2_ == MAX_BUFF_SIZE){
				is_data_buff2_full_.store(true);
				data_buff_idx1_ = 0;
				data_buff_idx2_ = 0;
			}
		}
	}
}

void WriteHelper::StartFileWriteThread(){
	write_thread_ =  thread(&WriteHelper::WriteToFileLoop, this);

	pthread_getschedparam(write_thread_.native_handle(), &policy_, &sch_);
	sch_.sched_priority = 5;
    pthread_setschedparam(write_thread_.native_handle(), SCHED_FIFO, &sch_);
    CPU_ZERO(&cpuset_);
    CPU_SET(1, &cpuset_);
    //CPU_SET(2, &cpuset_);

    int rc = pthread_setaffinity_np(write_thread_.native_handle(),
                                    sizeof(cpuset_), &cpuset_);    
	if (rc != 0) {
      std::cerr << "Error calling pthread_setaffinity_np on write thread: " << rc << "\n";
    }
}