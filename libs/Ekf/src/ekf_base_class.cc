#include "ekf_base_class.h"

template <typename T>
EkfBase<T>::EkfBase(size_t num_states, size_t num_meas, size_t num_states_sensor, T sample_time_s,
					MatrixInv<T> initial_state, MatrixInv<T> process_noise_q, MatrixInv<T> meas_noise_r, MatrixInv<T> initial_covariance_p,
					bool compute_q_each_iter, float process_noise_eps):
			 num_states_(num_states),
			 num_meas_(num_meas),
			 num_states_sensor_(num_states_sensor),
			 sample_time_s_(sample_time_s),
			 process_noise_q_(process_noise_q),
			 meas_noise_r_(meas_noise_r),
			 compute_q_each_iter_(compute_q_each_iter),
			 process_noise_eps_(process_noise_eps){	
		
		initial_state_ = initial_state;
		current_state_ = initial_state;
		computed_meas_ = MatrixInv<T> (num_meas_, 1);

		covariance_p_ = initial_covariance_p;

		state_jacobian_ = MatrixInv<T> (num_states_, num_states_);
		state_noise_jacobian_ = MatrixInv<T> (num_states_, num_states_, "eye");

		time_propagated_state_ = initial_state;
		meas_from_propogated_state_ = computed_meas_;
		meas_jacobian_ = MatrixInv<T> (num_meas_, num_states_);
		meas_noise_jacobian_ = MatrixInv<T> (num_meas_, num_meas_, "eye");

		kalman_eye_ = MatrixInv<T> (num_states_, num_states_, "eye");
		kalman_gain_seq_ = MatrixInv<T> (num_states_, 1);

		map_controls_to_state_ = MatrixInv<T> (num_states_, num_states_sensor_);
		process_noise_eps_matrix_ = MatrixInv<T> (num_states_, num_states_, "eye");
		process_noise_eps_matrix_ = process_noise_eps_matrix_*process_noise_eps_;
}

template <typename T>
EkfBase<T>::~EkfBase(){

}

template <typename T>
void EkfBase<T>::Run(const MatrixInv<T> &state_sensor_val, const MatrixInv<T> &meas_sensor_val, const bool meas_indices []){
	PropagateState(state_sensor_val);
	ComputeStateJacobian(state_sensor_val);

	GetMeas(meas_sensor_val, meas_indices);
	ComputeMeasJacobian(meas_sensor_val);

	//P = F*P*F' + L*Q*L';
	// if (compute_q_each_iter_){
	// 	ComputeControlToStateMap();
	// 	covariance_p_ = state_jacobian_*covariance_p_*state_jacobian_.Transpose() + map_controls_to_state_*process_noise_q_*map_controls_to_state_.Transpose() + process_noise_eps_matrix_;
	// }else{
	// 	covariance_p_ = state_jacobian_*covariance_p_*state_jacobian_.Transpose() + process_noise_q_;
	// }
	covariance_p_ = state_jacobian_*covariance_p_*state_jacobian_.Transpose() + process_noise_q_;

	current_state_ = time_propagated_state_;
	// sequentially update state with measurement
	for(size_t idx_r = 0; idx_r < num_meas_; idx_r++){
		// if (meas_indices[idx_r]){
			ComputeMeasFromState(idx_r);
			ComputeKalmanGainSequential(idx_r);
			MatrixInv<T> meas_jacobian_row = meas_jacobian_.GetRow(idx_r);
			current_state_ = current_state_ + kalman_gain_seq_*( computed_meas_(idx_r) - meas_from_propogated_state_(idx_r) );
			covariance_p_ = covariance_p_ - kalman_gain_seq_*meas_jacobian_row*covariance_p_;
		// }
	}
}

template <typename T>
void EkfBase<T>::Run(const MatrixInv<T> &state_sensor_val, MatrixInv<T> &meas_sensor_val, const bool meas_indices [], const long long int ins_dt_us, const long long int gps_time_us){
	ins_time_us_ += ins_dt_us;

	// printf("GPS TIME [us]: %lld, GPS Meas Indices: %d, %d\n", gps_time_us, meas_indices[0], meas_indices[1]);
	// Keep track of when sensor data are arriving, we are getting PVT from U-blox M8N so all the meas validity
	// flags are true
	if(meas_indices[1]){
		if(gps_prev_time_us_ == -1){
			// If this is the first time gps data was received set previous gps time to current time
			// set gps time between samples to zero. Also capture the ins time at the same time and set
			// ins time between gps update to 0
			gps_prev_time_us_ = gps_time_us;
			gps_dt_us_ = 0;
			ins_time_at_last_gps_update_us_ = ins_time_us_;
			ins_dt_bw_gps_update_us_ = 0;
		}else{
			// Find the time between last and current gps update and at the same time find out how far ins time has
			// progressed
			gps_dt_us_ = gps_time_us - gps_prev_time_us_;
			gps_prev_time_us_ = gps_time_us;
			ins_dt_bw_gps_update_us_  = ins_time_us_ - ins_time_at_last_gps_update_us_;
			ins_time_at_last_gps_update_us_ = ins_time_us_;
			printf("INS dt [us]: %lld, GPS dt [us]: %lld, 'GPS Meas Indices:%d, %d, %d, %d, %d, %d, %d\n", ins_dt_bw_gps_update_us_, gps_dt_us_, meas_indices[0], meas_indices[1], meas_indices[2],
					meas_indices[3], meas_indices[4], meas_indices[5], meas_indices[6]);
		}
	}else{
		ins_dt_bw_gps_update_us_  = 0;
		gps_dt_us_ = 1;
	}

	MatrixInv<T> meas_sensor_val_corrected = meas_sensor_val;

	// if ins time between gps updates is greater than time between gps updates then integrate gps position and velocity by the time
	// difference to make it current
	if(ins_dt_bw_gps_update_us_ > gps_dt_us_){
		PropagateState(state_sensor_val, meas_sensor_val_corrected, (ins_dt_bw_gps_update_us_ - gps_dt_us_));		
	}else{
		PropagateState(state_sensor_val);
	}

	ComputeStateJacobian(state_sensor_val);


	GetMeas(meas_sensor_val_corrected, meas_indices);
	ComputeMeasJacobian(meas_sensor_val_corrected);

	covariance_p_ = state_jacobian_*covariance_p_*state_jacobian_.Transpose() + process_noise_q_;

	current_state_ = time_propagated_state_;
	// sequentially update state with measurement
	for(size_t idx_r = 0; idx_r < num_meas_; idx_r++){
		if (meas_indices[idx_r]){
			ComputeMeasFromState(idx_r);
			ComputeKalmanGainSequential(idx_r);
			MatrixInv<T> meas_jacobian_row = meas_jacobian_.GetRow(idx_r);
			current_state_ = current_state_ + kalman_gain_seq_*( computed_meas_(idx_r) - meas_from_propogated_state_(idx_r) );
			covariance_p_ = covariance_p_ - kalman_gain_seq_*meas_jacobian_row*covariance_p_;
		}
	}

}

template <typename T>
inline void EkfBase<T>::ComputeKalmanGainSequential(size_t r_idx){
	kalman_gain_seq_ = covariance_p_*meas_jacobian_.Transpose().GetCol(r_idx);
	MatrixInv<T> k_denominator = meas_jacobian_.GetRow(r_idx)*kalman_gain_seq_ + meas_noise_r_(r_idx, r_idx);
	kalman_gain_seq_ = kalman_gain_seq_/k_denominator(0);
}

// Explicit template instantiation
template class EkfBase<float>;
template class EkfBase<double>;