#include "ekf_base_class.h"

template <typename T>
EkfBase<T>::EkfBase(size_t num_states, size_t num_meas, size_t num_states_sensor, T sample_time_s,
					MatrixInv<T> initial_state, MatrixInv<T> process_noise_q, MatrixInv<T> meas_noise_r):
			 num_states_(num_states),
			 num_meas_(num_meas),
			 num_states_sensor_(num_states_sensor),
			 sample_time_s_(sample_time_s),
			 process_noise_q_(process_noise_q),
			 meas_noise_r_(meas_noise_r){	
		
		initial_state_ = initial_state;
		current_state_ = initial_state;
		computed_meas_ = MatrixInv<T> (num_meas_, 1);

		covariance_p_ = process_noise_q_;

		state_jacobian_ = MatrixInv<T> (num_states_, num_states_);
		state_noise_jacobian_ = MatrixInv<T> (num_states_, num_states_);

		time_propagated_state_ = initial_state;
		meas_from_propogated_state_ = computed_meas_;
		meas_jacobian_ = MatrixInv<T> (num_meas_, num_states_);
		meas_noise_jacobian_ = MatrixInv<T> (num_meas_, num_meas_, "eye");

		kalman_eye_ = MatrixInv<T> (num_states_, num_states_, "eye");
}

template <typename T>
EkfBase<T>::~EkfBase(){

}

template <typename T>
void EkfBase<T>::run(MatrixInv<T> state_sensor_val, MatrixInv<T> meas_sensor_val){
	PropagateState(current_state_, state_sensor_val);
	ComputeStateJacobian(current_state_, state_sensor_val);
	ComputeStateNoiseJacobian(current_state_);
	GetMeas(meas_sensor_val);
/*	ComputeMeasJacobian(MatrixInv<T> meas_sensor_val) = 0;
	ComputeMeasNoiseJacobian(MatrixInv<T> meas_sensor_val) = 0;*/
	ComputeMeasFromState(time_propagated_state_);

	//P = F*P*F' + L*Q*L';
	covariance_p_ = state_jacobian_*covariance_p_*state_jacobian_.Transpose() + state_noise_jacobian_*process_noise_q_*state_noise_jacobian_.Transpose();

	//Inverse part of Kalman gain H*P*H' + M*R*M'
	inv_part_of_kalman_gain_ = meas_jacobian_*covariance_p_*meas_jacobian_.Transpose() + meas_noise_jacobian_*meas_noise_r_*meas_noise_jacobian_.Transpose();
	// Calculate Kalman Gain P*H'*(H*P*H' + M*R*M')^-1
	kalman_gain_ = covariance_p_*meas_jacobian_.Transpose()*inv_part_of_kalman_gain_.Inverse();

	//update state x = x + K(y - h(x))
	current_state_ = time_propagated_state_ + kalman_gain_*(computed_meas_ - meas_from_propogated_state_);
	//update covariance P = (I - K*H)*P
	covariance_p_ = (kalman_eye_ - kalman_gain_*meas_jacobian_)*covariance_p_;
}

// Explicit template instantiation
template class EkfBase<float>;
template class EkfBase<double>;