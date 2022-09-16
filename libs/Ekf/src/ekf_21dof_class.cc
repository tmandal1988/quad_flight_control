#include "ekf_21dof_class.h"

template <typename T>
Ekf21Dof<T>::Ekf21Dof(T sample_time_s, MatrixBase<T> initial_state, MatrixBase<T> process_noise_q, MatrixBase<T> meas_noise_r):
			EkfBase<T>(NUM_STATES, NUM_MEAS, NUM_STATES_SENSOR, sample_time_s, initial_state, process_noise_q, meas_noise_r){	
}

template <typename T>
Ekf21Dof<T>::~Ekf21Dof(){

}

template <typename T>
MatrixBase<T> Ekf21Dof<T>::GetStateFunction(MatrixBase<T> previous_state, MatrixBase<T> current_states_meas){
	MatrixBase<T> A(this->num_states_, this->num_states_, "eye");
	return A;
}

template <typename T>
MatrixBase<T> Ekf21Dof<T>::GetStateJacobian(MatrixBase<T> previous_state, MatrixBase<T> current_states_meas){
	MatrixBase<T> A(this->num_states_, this->num_states_, "eye");
	return A;
}

template <typename T>
MatrixBase<T> Ekf21Dof<T>::GetStateNoiseJacobian(MatrixBase<T> previous_state, MatrixBase<T> current_states_meas){
	MatrixBase<T> A(this->num_states_, this->num_states_, "eye");
	return A;
}

template <typename T>
MatrixBase<T> Ekf21Dof<T>::GetMeasFunction(MatrixBase<T> current_meas){
	MatrixBase<T> A(this->num_states_, this->num_states_, "eye");
	return A;
}

template <typename T>
MatrixBase<T> Ekf21Dof<T>::GetMeasJacobian(MatrixBase<T> current_meas){
	MatrixBase<T> A(this->num_states_, this->num_states_, "eye");
	return A;
}

template <typename T>
MatrixBase<T> Ekf21Dof<T>::GetMeasNoiseJacobian(MatrixBase<T> current_meas){
	MatrixBase<T> A(this->num_states_, this->num_states_, "eye");
	return A;
}

// Explicit template instantiation
template class Ekf21Dof<float>;
template class Ekf21Dof<double>;