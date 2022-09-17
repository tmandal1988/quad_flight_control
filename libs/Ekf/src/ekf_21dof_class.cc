#include "ekf_21dof_class.h"

template <typename T>
Ekf21Dof<T>::Ekf21Dof(T sample_time_s, MatrixBase<T> initial_state, MatrixBase<T> process_noise_q, MatrixBase<T> meas_noise_r):
			EkfBase<T>(NUM_STATES, NUM_MEAS, NUM_STATES_SENSOR, sample_time_s, initial_state, process_noise_q, meas_noise_r){	


	T s_phi = 0;;
	T c_phi = 0;

	T s_theta = 0;
	T c_theta = 0;
	T t_theta = 0;
	T sc_theta = 0;

	T s_psi = 0;
	T c_psi = 0;
}

template <typename T>
Ekf21Dof<T>::~Ekf21Dof(){

}

template <typename T>
void Ekf21Dof<T>::PropagateState(MatrixBase<T> previous_state, MatrixBase<T> current_states_meas){
	// precalculate trignometric values
	s_phi = sin(previous_state(0));
	c_phi = cos(previous_state(0));

	s_theta = sin(previous_state(1));
	c_theta = cos(previous_state(1));
	t_theta = tan(previous_state(1));
	sc_theta = 1/c_theta;

	s_psi = sin(previous_state(2));
	c_psi = cos(previous_state(2));
	
	// propagate phi
	this->time_propagated_state_(0) = previous_state(0) + this->sample_time_s_*((current_states_meas(0) - previous_state(3)) + 
																		s_phi*t_theta*(current_states_meas(1) - previous_state(4)) +
																		c_phi*t_theta*(current_states_meas(2) - previous_state(5))) ;

	// propagate theta
	this->time_propagated_state_(1) = previous_state(1) + this->sample_time_s_*( c_phi*(current_states_meas(1) - previous_state(4)) - 
																		 s_phi*(current_states_meas(2) - previous_state(5)) );

	// propagate psi
	this->time_propagated_state_(2) = previous_state(2) + this->sample_time_s_*( s_phi*sc_theta*(current_states_meas(1) - previous_state(4)) - 
																		 c_phi*sc_theta*(current_states_meas(2) - previous_state(5)) );

	// propagate gyro bias
	this->time_propagated_state_(3) = previous_state(3);
	this->time_propagated_state_(4) = previous_state(4);
	this->time_propagated_state_(5) = previous_state(5);

	// propagate NED X position
	this->time_propagated_state_(6) = previous_state(6) + this->sample_time_s_*previous_state(9);
	// propagate NED Y position
	this->time_propagated_state_(7) = previous_state(7) + this->sample_time_s_*previous_state(10);
	// propagate NED Z position
	this->time_propagated_state_(8) = previous_state(8) + this->sample_time_s_*previous_state(11);

	// propagate NED X velocity
	this->time_propagated_state_(9) = previous_state(9) + this->sample_time_s_*( c_theta*c_psi*(current_states_meas(3) - previous_state(12)) +
																		 c_theta*s_psi*(current_states_meas(4) - previous_state(13)) -
																		 s_theta*(current_states_meas(5) - previous_state(14)) );
	// propagate NED Y velocity
	this->time_propagated_state_(10) = previous_state(10) + this->sample_time_s_*( (s_phi*s_theta*c_psi - c_phi*s_psi)*(current_states_meas(3) - previous_state(12)) +
																		   (s_phi*s_theta*s_psi + c_phi*c_psi)*(current_states_meas(4) - previous_state(13)) +
																		    s_phi*c_theta*(current_states_meas(5) - previous_state(14)) );
	// propagate NED Z velocity
	this->time_propagated_state_(11) = previous_state(11) + this->sample_time_s_*( (c_phi*s_theta*c_psi + s_phi*s_psi)*(current_states_meas(3) - previous_state(12)) +
																		   (c_phi*s_theta*s_psi - s_phi*c_psi)*(current_states_meas(4) - previous_state(13)) +
																		    c_phi*c_theta*(current_states_meas(5) - previous_state(14)) );

	// propagate accel bias
	this->time_propagated_state_(12) = previous_state(12);
	this->time_propagated_state_(13) = previous_state(13);
	this->time_propagated_state_(14) = previous_state(14);

	// filtered angular rates
	this->time_propagated_state_(15) = current_states_meas(0) - previous_state(3);
	this->time_propagated_state_(16) = current_states_meas(1) - previous_state(4);
	this->time_propagated_state_(17) = current_states_meas(2) - previous_state(5);

	//filtered accel
	this->time_propagated_state_(18) = current_states_meas(3) - previous_state(12);
	this->time_propagated_state_(19) = current_states_meas(4) - previous_state(13);
	this->time_propagated_state_(20) = current_states_meas(5) - previous_state(14);
}

template <typename T>
void Ekf21Dof<T>::ComputeStateJacobian(MatrixBase<T> previous_state, MatrixBase<T> current_states_meas){

	// 1st row
	this->state_jacobian_(0, 0) = 1 + this->sample_time_s_*( c_phi*t_theta*(current_states_meas(1) - previous_state(4)) -
															 s_phi*t_theta*(current_states_meas(2) - previous_state(5)) );
	this->state_jacobian_(0, 1) = this->sample_time_s_*(s_phi*sc_theta*sc_theta*(current_states_meas(1) - previous_state(4)) +
														c_phi*sc_theta*sc_theta*(current_states_meas(2) - previous_state(5)) );
	
	this->state_jacobian_(0, 3) = -this->sample_time_s_;
	this->state_jacobian_(0, 4) = -this->sample_time_s_*s_phi*t_theta;
	this->state_jacobian_(0, 5) = -this->sample_time_s_*c_phi*t_theta;

	this->state_jacobian_(0, 15) = this->sample_time_s_;
	this->state_jacobian_(0, 16) = this->sample_time_s_*s_phi*t_theta;
	this->state_jacobian_(0, 17) = this->sample_time_s_*c_phi*t_theta;

	//2nd row
	this->state_jacobian_(1, 0) = this->sample_time_s_*(-s_phi*(current_states_meas(1) - previous_state(4)) -
														 c_phi*(current_states_meas(2) - previous_state(5)) );
	this->state_jacobian_(1, 1) = 1;

	this->state_jacobian_(1, 4) = -this->sample_time_s_*c_phi;
	this->state_jacobian_(1, 5) = this->sample_time_s_*s_phi;

	this->state_jacobian_(1, 16) = this->sample_time_s_*c_phi;
	this->state_jacobian_(1, 17) = -this->sample_time_s_*s_phi;

	//3rd row
	this->state_jacobian_(2, 0) = this->sample_time_s_*( c_phi*sc_theta*(current_states_meas(1) - previous_state(4)) - 
														 s_phi*sc_theta*(current_states_meas(2) - previous_state(5)) );
	this->state_jacobian_(2, 1) = this->sample_time_s_*( s_phi*sc_theta*t_theta*(current_states_meas(1) - previous_state(4)) +
														 c_phi*sc_theta*t_theta*(current_states_meas(2) - previous_state(5)) );
	this->state_jacobian_(2, 2) = 1;

	this->state_jacobian_(2, 4) = -this->sample_time_s_*s_phi*sc_theta;
	this->state_jacobian_(2, 5) = -this->sample_time_s_*c_phi*sc_theta;

	this->state_jacobian_(2, 16) = this->sample_time_s_*s_phi*sc_theta;
	this->state_jacobian_(2, 17) = this->sample_time_s_*c_phi*sc_theta;

	//4th row
	this->state_jacobian_(3, 3) = 1;
	//5th row
	this->state_jacobian_(4, 4) = 1;
	//6th row
	this->state_jacobian_(5, 5) = 1;

	//7th row
	this->state_jacobian_(6, 6) = 1;
	this->state_jacobian_(6, 9) = this->sample_time_s_;

	//8th row
	this->state_jacobian_(7, 7) = 1;
	this->state_jacobian_(7, 10) = this->sample_time_s_;

	//9th row
	this->state_jacobian_(8, 8) = 1;
	this->state_jacobian_(8, 11) = this->sample_time_s_;

	//10th row
	this->state_jacobian_(9, 1) = this->sample_time_s_*( -s_theta*c_psi*(current_states_meas(3) - previous_state(12)) -
														 s_theta*s_psi*(current_states_meas(4) - previous_state(13)) -
														 c_theta*(current_states_meas(5) - previous_state(14)) );
	this->state_jacobian_(9, 2) = this->sample_time_s_*( -c_theta*s_psi*(current_states_meas(3) - previous_state(12)) + 
														 c_theta*c_psi*(current_states_meas(5) - previous_state(14)) );

	this->state_jacobian_(9, 9) = 1;

	this->state_jacobian_(9, 12) = -this->sample_time_s_*c_theta*c_psi;
	this->state_jacobian_(9, 13) = -this->sample_time_s_*c_theta*s_psi;
	this->state_jacobian_(9, 14) = this->sample_time_s_*s_theta;

	this->state_jacobian_(9, 18) = this->sample_time_s_*c_theta*c_psi;
	this->state_jacobian_(9, 19) = this->sample_time_s_*c_theta*s_psi;
	this->state_jacobian_(9, 20) = -this->sample_time_s_*s_theta;


	//11th row
	this->state_jacobian_(10, 0) = this->sample_time_s_*( (c_phi*s_theta*c_psi + s_phi*s_psi)*(current_states_meas(3) - previous_state(12)) + 
									 (c_phi*s_theta*s_psi - s_phi*c_psi)*(current_states_meas(4) - previous_state(13)) + 
									  c_phi*c_theta*(current_states_meas(5) - previous_state(14)) );
	this->state_jacobian_(10, 1) = this->sample_time_s_*( s_phi*c_theta*c_psi*(current_states_meas(3) - previous_state(12)) - 
														  s_phi*c_theta*s_psi*(current_states_meas(4) - previous_state(13)) -
														  s_phi*s_theta*(current_states_meas(5) - previous_state(14)) );
	this->state_jacobian_(10, 2) = this->sample_time_s_*( (-s_phi*s_theta*s_psi - c_phi*c_psi)*(current_states_meas(3) - previous_state(12)) +
														  (s_phi*s_theta*c_psi - c_phi*s_psi)*(current_states_meas(4) - previous_state(13)) );

	this->state_jacobian_(10, 10) = 1;

	this->state_jacobian_(10, 12) = this->sample_time_s_*( c_phi*s_psi - s_phi*s_theta*c_psi );
	this->state_jacobian_(10, 13) = this->sample_time_s_*( -s_phi*s_theta*s_psi - c_phi*c_psi );
	this->state_jacobian_(10, 14) = this->sample_time_s_*(-s_phi*c_theta);

	this->state_jacobian_(10, 18) = this->sample_time_s_*(s_phi*s_theta*c_psi - c_phi*s_psi);
	this->state_jacobian_(10, 19) = this->sample_time_s_*(s_phi*s_theta*s_psi + c_phi*c_psi);
	this->state_jacobian_(10, 20) = this->sample_time_s_*(s_phi*c_theta);

	//12th row
	this->state_jacobian_(11, 0) = this->sample_time_s_*( (-s_phi*s_theta*c_psi + c_phi*s_psi)*(current_states_meas(3) - previous_state(12)) +
														  (-s_phi*s_theta*s_psi - c_phi*c_psi)*(current_states_meas(4) - previous_state(13)) - 
														  s_phi*c_theta*(current_states_meas(5) - previous_state(14)) );
	this->state_jacobian_(11, 1) = this->sample_time_s_*( c_phi*c_theta*c_psi*(current_states_meas(3) - previous_state(12)) +
														  c_phi*c_theta*s_psi*(current_states_meas(4) - previous_state(13)) );
	this->state_jacobian_(11, 2) = this->sample_time_s_*( (-c_phi*s_theta*s_psi + s_phi*c_psi)*(current_states_meas(3) - previous_state(12)) +
														  (c_phi*s_theta*c_psi + s_phi*s_psi)*(current_states_meas(4) - previous_state(13)) );

	this->state_jacobian_(11, 11) = 1;
	this->state_jacobian_(11, 12) = this->sample_time_s_*(-c_phi*s_theta*c_psi - s_phi*s_psi);
	this->state_jacobian_(11, 13) = this->sample_time_s_*(s_phi*c_psi - c_phi*s_theta*s_psi);
	this->state_jacobian_(11, 14) = this->sample_time_s_*(-c_phi*c_theta);

	this->state_jacobian_(11, 18) = this->sample_time_s_*(c_phi*s_theta*c_psi + s_phi*s_psi);
	this->state_jacobian_(11, 19) = this->sample_time_s_*(-s_phi*c_psi + c_phi*s_theta*s_psi);
	this->state_jacobian_(11, 20) = this->sample_time_s_*(c_phi*c_theta);

	//13th row
	this->state_jacobian_(12, 12) = 1;
	//14th row
	this->state_jacobian_(13, 13) = 1;
	//15th row
	this->state_jacobian_(14, 14) = 1;

	//16th row
	this->state_jacobian_(15, 3) = -1;
	this->state_jacobian_(15, 15) = 1;
	//17th row
	this->state_jacobian_(16, 4) = -1;
	this->state_jacobian_(16, 16) = -1;
	//18th row
	this->state_jacobian_(17, 5) = -1;
	this->state_jacobian_(17, 17) = 1;
	//19th row
	this->state_jacobian_(18, 12) = -1;
	this->state_jacobian_(18, 18) = 1;
	//20th row
	this->state_jacobian_(19, 13) = -1;
	this->state_jacobian_(19, 19) = 1;
	//21th row
	this->state_jacobian_(20, 14) = -1;
	this->state_jacobian_(20, 20) = 1;

}

template <typename T>
MatrixBase<T> Ekf21Dof<T>::ComputeStateNoiseJacobian(MatrixBase<T> previous_state, MatrixBase<T> current_states_meas){
	MatrixBase<T> A(this->num_states_, this->num_states_, "eye");
	return A;
}

template <typename T>
MatrixBase<T> Ekf21Dof<T>::GetMeas(MatrixBase<T> current_meas){
	MatrixBase<T> A(this->num_states_, this->num_states_, "eye");
	return A;
}

template <typename T>
MatrixBase<T> Ekf21Dof<T>::ComputeMeasJacobian(MatrixBase<T> current_meas){
	MatrixBase<T> A(this->num_states_, this->num_states_, "eye");
	return A;
}

template <typename T>
MatrixBase<T> Ekf21Dof<T>::ComputeMeasNoiseJacobian(MatrixBase<T> current_meas){
	MatrixBase<T> A(this->num_states_, this->num_states_, "eye");
	return A;
}

// Explicit template instantiation
template class Ekf21Dof<float>;
template class Ekf21Dof<double>;