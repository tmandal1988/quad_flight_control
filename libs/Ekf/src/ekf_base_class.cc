#include "ekf_base_class.h"

template <typename T>
EkfBase<T>::EkfBase(size_t num_states, size_t num_meas, size_t num_states_sensor, T sample_time_s,
					MatrixBase<T> initial_state, MatrixBase<T> process_noise_q, MatrixBase<T> meas_noise_r):
			 num_states_(num_states),
			 num_meas_(num_meas),
			 num_states_sensor_(num_states_sensor),
			 sample_time_s_(sample_time_s),
			 process_noise_q_(process_noise_q),
			 meas_noise_r_(meas_noise_r){	
		
		initial_state_ = initial_state;
		current_state_ = initial_state;
		current_state_sensor_val_ = MatrixBase<T> (num_states_sensor_, 1);
		current_meas_ = MatrixBase<T> (num_meas_, 1);

		covariance_p_ = MatrixBase<T> (num_states_, num_states_, "eye");

}

template <typename T>
EkfBase<T>::~EkfBase(){

}

template <typename T>
void EkfBase<T>::run(){
	MatrixBase<T> stateVal = GetStateFunction(current_state_, current_state_sensor_val_);
	current_state_ = stateVal;
}

// Explicit template instantiation
template class EkfBase<float>;
template class EkfBase<double>;