#include "ekf_15dof_class.h"

template <typename T>
Ekf15Dof<T>::Ekf15Dof(T sample_time_s, MatrixInv<T> initial_state, MatrixInv<T> process_noise_q, MatrixInv<T> meas_noise_r):
			EkfBase<T>(NUM_STATES, NUM_MEAS, NUM_STATES_SENSOR, sample_time_s, initial_state, process_noise_q, meas_noise_r){	

	this->meas_jacobian_(0, 2) = 1;
	this->meas_jacobian_(1, 6) = 1;
	this->meas_jacobian_(2, 7) = 1;
	this->meas_jacobian_(3, 8) = 1;
	this->meas_jacobian_(4, 9) = 1;
	this->meas_jacobian_(5, 10) = 1;
	this->meas_jacobian_(6, 11) = 1;


	T s_phi = 0;;
	T c_phi = 0;

	T s_theta = 0;
	T c_theta = 0;
	T t_theta = 0;
	T sc_theta = 0;

	T s_psi = 0;
	T c_psi = 0;

	g = MatrixInv<T>(3, 1);
	g(2) = 9.81;
}

template <typename T>
Ekf15Dof<T>::~Ekf15Dof(){

}

template <typename T>
void Ekf15Dof<T>::ComputeTrignometricValues(){
	// precalculate trignometric values
	s_phi = sin(this->current_state_(0));
	c_phi = cos(this->current_state_(0));

	s_theta = sin(this->current_state_(1));
	c_theta = cos(this->current_state_(1));
	t_theta = tan(this->current_state_(1));
	sc_theta = 1/c_theta;

	s_psi = sin(this->current_state_(2));
	c_psi = cos(this->current_state_(2));

}

template <typename T>
void Ekf15Dof<T>::PropagateState(MatrixInv<T> previous_state, MatrixInv<T> state_sensor_val){
	ComputeTrignometricValues();
	T dt = this->sample_time_s_;

	MatrixInv<T> propogated_euler_ang(3, 1);
	MatrixInv<T> propogated_velocity(3, 1);
	MatrixInv<T> body_rates_2_euler_rates(3, 3);
	MatrixInv<T> c_b2ned(3, 3);
	MatrixInv<T> angular_rate_vector(3, 1);
	MatrixInv<T> accel_vector(3, 1);

	body_rates_2_euler_rates(0, 0) = 1;
	body_rates_2_euler_rates(0, 1) = s_phi*t_theta;
	body_rates_2_euler_rates(0, 2) =  c_phi*t_theta;

	body_rates_2_euler_rates(1, 0) = 0;
	body_rates_2_euler_rates(1, 1) = c_phi;
	body_rates_2_euler_rates(1, 2) = -s_phi;

	body_rates_2_euler_rates(2, 0) = 0;
	body_rates_2_euler_rates(2, 1) = s_phi/c_theta;
	body_rates_2_euler_rates(2, 2) = c_phi/c_theta;

	angular_rate_vector(0) = state_sensor_val(0) - previous_state(3);
	angular_rate_vector(1) = state_sensor_val(1) - previous_state(4);
	angular_rate_vector(2) = state_sensor_val(2) - previous_state(5);

	propogated_euler_ang = body_rates_2_euler_rates*angular_rate_vector*dt;


	c_b2ned(0, 0) = c_psi*c_theta;
	c_b2ned(0, 1) = c_psi*s_phi*s_theta - c_phi*s_psi;
	c_b2ned(0, 2) = s_phi*s_psi + c_phi*c_psi*s_theta;

	c_b2ned(1, 0) = c_theta*s_psi;
	c_b2ned(1, 1) = c_phi*c_psi + s_phi*s_psi*s_theta;
	c_b2ned(1, 2) = c_phi*s_psi*s_theta - c_psi*s_phi;

	c_b2ned(2, 0) = -s_theta;
	c_b2ned(2, 1) = c_theta*s_phi;
	c_b2ned(2, 2) = c_phi*c_theta;

	accel_vector(0) = state_sensor_val(3) - previous_state(12);
	accel_vector(1) = state_sensor_val(4) - previous_state(13);
	accel_vector(2) = state_sensor_val(5) - previous_state(14);

	propogated_velocity = (c_b2ned*accel_vector + g)*dt;

	// propagate phi
	this->time_propagated_state_(0) = previous_state(0) + propogated_euler_ang(0);
	this->time_propagated_state_(1) = previous_state(1) + propogated_euler_ang(1);
	this->time_propagated_state_(2) = previous_state(2) + propogated_euler_ang(2);
	// this->time_propagated_state_(0) = previous_state(0) +  dt*( state_sensor_val(0) - previous_state(3) + c_phi*t_theta*( state_sensor_val(2) - previous_state(5) ) +
	// 															s_phi*t_theta*(state_sensor_val(1) - previous_state(4)) );

	// // propagate theta
	// this->time_propagated_state_(1) = previous_state(1)  + dt*( c_phi*(state_sensor_val(1) - previous_state(4)) - s_phi*(state_sensor_val(2) - previous_state(5)) );
	// // propagate psi
	// this->time_propagated_state_(2) = previous_state(2)  + dt*( (c_phi*(state_sensor_val(2) - previous_state(5)))/c_theta + (s_phi*(state_sensor_val(1) - previous_state(4)))/c_theta );

	// propagate gyro bias
	this->time_propagated_state_(3) = previous_state(3);
	this->time_propagated_state_(4) = previous_state(4);
	this->time_propagated_state_(5) = previous_state(5);

	// propagate NED X position
	this->time_propagated_state_(6) = previous_state(6) + dt*previous_state(9);
	// propagate NED Y position
	this->time_propagated_state_(7) = previous_state(7) + dt*previous_state(10);
	// propagate NED Z position
	this->time_propagated_state_(8) = previous_state(8) + dt*previous_state(11);

	// propagate NED X velocity
	this->time_propagated_state_(9) = previous_state(9) + propogated_velocity(0);
	this->time_propagated_state_(10) = previous_state(10) + propogated_velocity(1);
	this->time_propagated_state_(11) = previous_state(11) + propogated_velocity(2);
	// this->time_propagated_state_(9) = previous_state(9) + dt*((s_phi*s_psi + c_phi*c_psi*s_theta)*(state_sensor_val(5) - previous_state(14)) - 
	// 														  (c_phi*s_psi - c_psi*s_phi*s_theta)*(state_sensor_val(4) - previous_state(13)) + 
	// 														  c_psi*c_theta*(state_sensor_val(3) - previous_state(12)));
	// // propagate NED Y velocity
	// this->time_propagated_state_(10) = previous_state(10) + dt*((c_phi*c_psi + s_phi*s_psi*s_theta)*(state_sensor_val(4) - previous_state(13)) - 
	// 															(c_psi*s_phi - c_phi*s_psi*s_theta)*(state_sensor_val(5) - previous_state(14)) + 
	// 															c_theta*s_psi*(state_sensor_val(3) - previous_state(12)));
	// // propagate NED Z velocity
	// this->time_propagated_state_(11) = previous_state(11) + dt*(c_phi*c_theta*(state_sensor_val(5) - previous_state(14)) - 
	// 															s_theta*(state_sensor_val(3) - previous_state(12)) + 
	// 															c_theta*s_phi*(state_sensor_val(4) - previous_state(13)) + 9.81);

	// // propagate accel bias
	this->time_propagated_state_(12) = previous_state(12);
	this->time_propagated_state_(13) = previous_state(13);
	this->time_propagated_state_(14) = previous_state(14);
}

template <typename T>
void Ekf15Dof<T>::ComputeStateJacobian(MatrixInv<T> previous_state, MatrixInv<T> state_sensor_val){
	T dt = this->sample_time_s_;

	// 1st row
	this->state_jacobian_(0, 0) = 1 + dt*( c_phi*t_theta*(state_sensor_val(1) - previous_state(4)) - s_phi*t_theta*(state_sensor_val(2) - previous_state(5)) );
	this->state_jacobian_(0, 1) = dt*(s_phi*sc_theta*sc_theta*(state_sensor_val(1) - previous_state(4)) + c_phi*sc_theta*sc_theta*(state_sensor_val(2) - previous_state(5)) );
	
	this->state_jacobian_(0, 3) = -dt;
	this->state_jacobian_(0, 4) = -dt*s_phi*t_theta;
	this->state_jacobian_(0, 5) = -dt*c_phi*t_theta;

	//2nd row
	this->state_jacobian_(1, 0) = dt*(-s_phi*(state_sensor_val(1) - previous_state(4)) - c_phi*(state_sensor_val(2) - previous_state(5)) );
	this->state_jacobian_(1, 1) = 1;

	this->state_jacobian_(1, 4) = -dt*c_phi;
	this->state_jacobian_(1, 5) = dt*s_phi;

	//3rd row
	this->state_jacobian_(2, 0) = dt*( c_phi*sc_theta*(state_sensor_val(1) - previous_state(4)) - s_phi*sc_theta*(state_sensor_val(2) - previous_state(5)) );
	this->state_jacobian_(2, 1) = dt*( s_phi*sc_theta*t_theta*(state_sensor_val(1) - previous_state(4)) + c_phi*sc_theta*t_theta*(state_sensor_val(2) - previous_state(5)) );
	this->state_jacobian_(2, 2) = 1;

	this->state_jacobian_(2, 4) = -dt*s_phi*sc_theta;
	this->state_jacobian_(2, 5) = -dt*c_phi*sc_theta;

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
	this->state_jacobian_(9, 0) = dt*( (s_phi*s_psi + c_phi*s_theta*c_psi)*(state_sensor_val(4) - previous_state(13)) + 
                					   (c_phi*s_psi - s_phi*s_theta*c_psi)*(state_sensor_val(5) - previous_state(14)) );
	this->state_jacobian_(9, 1) = dt*( -s_theta*c_psi*(state_sensor_val(3) - previous_state(12)) + s_phi*c_theta*c_psi*(state_sensor_val(4) - previous_state(13)) + 
										c_phi*c_theta*c_psi*(state_sensor_val(5) - previous_state(14)) );
	this->state_jacobian_(9, 2) = dt*( -c_theta*s_psi*(state_sensor_val(3) - previous_state(12)) + (-c_phi*c_psi - s_phi*s_theta*s_psi)*(state_sensor_val(4) - previous_state(13)) +
										(s_phi*c_psi - c_phi*s_theta*s_psi)*(state_sensor_val(5) - previous_state(14)) );

	this->state_jacobian_(9, 9) = 1;

	this->state_jacobian_(9, 12) = -dt*c_theta*c_psi;
	this->state_jacobian_(9, 13) = -dt*(-c_phi*s_psi + s_phi*s_theta*c_psi);
	this->state_jacobian_(9, 14) = -dt*(s_phi*s_psi + c_phi*s_theta*c_psi);


	//11th row
	this->state_jacobian_(10, 0) = dt*( (-s_phi*c_psi + c_phi*s_theta*s_psi)*(state_sensor_val(4) - previous_state(13)) + 
    									 (-c_phi*c_psi - s_phi*s_theta*s_psi)*(state_sensor_val(5) - previous_state(14)) );
	this->state_jacobian_(10, 1) = dt*( -s_theta*s_psi*(state_sensor_val(3) - previous_state(12)) + s_phi*c_theta*s_psi*(state_sensor_val(4) - previous_state(13)) - 
   										 c_phi*c_theta*s_psi*(state_sensor_val(5) - previous_state(14)) );
	this->state_jacobian_(10, 2) = dt*( c_theta*c_psi*(state_sensor_val(3) - previous_state(12)) + 
    									(-c_phi*s_psi + s_phi*s_theta*c_psi)*(state_sensor_val(4) - previous_state(13)) + 
    									(s_phi*s_psi + c_phi*s_theta*c_psi)*(state_sensor_val(5) - previous_state(14)) );
	this->state_jacobian_(10, 10) = 1;

	this->state_jacobian_(10, 12) = -dt*c_theta*s_psi;
	this->state_jacobian_(10, 13) = -dt*(c_phi*c_psi + s_phi*s_theta*s_psi);
	this->state_jacobian_(10, 14) = -dt*(-s_phi*c_psi + c_phi*s_theta*s_psi);

	//12th row
	this->state_jacobian_(11, 0) = dt*( c_phi*c_theta*(state_sensor_val(4) - previous_state(13)) - s_phi*c_theta*(state_sensor_val(5) - previous_state(14)) );
	this->state_jacobian_(11, 1) = dt*( -c_theta*(state_sensor_val(3) - previous_state(12)) - s_phi*s_theta*(state_sensor_val(4) - previous_state(13)) - 
    									 c_phi*s_theta*(state_sensor_val(5) - previous_state(14)) );
	this->state_jacobian_(11, 11) = 1;
	this->state_jacobian_(11, 12) = dt*s_theta;
	this->state_jacobian_(11, 13) = -dt*s_phi*c_theta;
	this->state_jacobian_(11, 14) = -dt*(c_phi*c_theta);


	//13th row
	this->state_jacobian_(12, 12) = 1;
	//14th row
	this->state_jacobian_(13, 13) = 1;
	//15th row
	this->state_jacobian_(14, 14) = 1;
}

template <typename T>
void Ekf15Dof<T>::ComputeStateNoiseJacobian(MatrixInv<T> previous_state){
	//1st row
	this->state_noise_jacobian_(0, 15) = this->sample_time_s_;
	this->state_noise_jacobian_(0, 16) = this->sample_time_s_*s_phi*t_theta;
	this->state_noise_jacobian_(0, 17) = this->sample_time_s_*c_phi*t_theta;

	//2nd row
	this->state_noise_jacobian_(1, 16) = this->sample_time_s_*c_phi;
	this->state_noise_jacobian_(1, 17) = -this->sample_time_s_*s_phi;

	//3rd row
	this->state_noise_jacobian_(2, 16) = this->sample_time_s_*s_phi*sc_theta;
	this->state_noise_jacobian_(2, 17) = -this->sample_time_s_*c_phi*sc_theta;

	//4th row
	this->state_noise_jacobian_(3, 3) = 1;
	//5th row
	this->state_noise_jacobian_(4, 4) = 1;
	//6th row
	this->state_noise_jacobian_(5, 5) = 1;

	//7th row
	this->state_noise_jacobian_(6, 18) = pow(this->sample_time_s_, 2)*c_theta*c_psi;
	this->state_noise_jacobian_(6, 19) = pow(this->sample_time_s_, 2)*c_theta*s_psi;
	this->state_noise_jacobian_(6, 20) = -pow(this->sample_time_s_, 2)*s_theta;

	//8th row
	this->state_noise_jacobian_(7, 18) = pow(this->sample_time_s_, 2)*(s_phi*s_theta*c_psi - c_phi*s_psi);
	this->state_noise_jacobian_(7, 19) = pow(this->sample_time_s_, 2)*(s_phi*s_theta*s_psi + c_phi*c_psi);
	this->state_noise_jacobian_(7, 20) = pow(this->sample_time_s_, 2)*c_theta*s_phi;

	//9th row
	this->state_noise_jacobian_(8, 18) = pow(this->sample_time_s_, 2)*(c_phi*s_theta*c_psi + s_phi*s_psi);
	this->state_noise_jacobian_(8, 19) = pow(this->sample_time_s_, 2)*(c_phi*s_theta*s_psi - s_phi*c_psi);
	this->state_noise_jacobian_(8, 20) = pow(this->sample_time_s_, 2)*c_theta*c_phi;

	//10th row
	this->state_noise_jacobian_(9, 18) = this->sample_time_s_*c_theta*c_psi;
	this->state_noise_jacobian_(9, 19) = this->sample_time_s_*c_theta*s_psi;
	this->state_noise_jacobian_(9, 20) = -this->sample_time_s_*s_theta;

	//11th row
	this->state_noise_jacobian_(10, 18) = this->sample_time_s_*(s_phi*s_theta*c_psi - c_phi*s_psi);
	this->state_noise_jacobian_(10, 19) = this->sample_time_s_*(s_phi*s_theta*s_psi + c_phi*c_psi);
	this->state_noise_jacobian_(10, 20) = this->sample_time_s_*c_theta*s_phi;

	//12th row
	this->state_noise_jacobian_(11, 18) = this->sample_time_s_*(c_phi*s_theta*c_psi + s_phi*s_psi);
	this->state_noise_jacobian_(11, 19) = this->sample_time_s_*(c_phi*s_theta*s_psi - s_phi*c_psi);
	this->state_noise_jacobian_(11, 20) = this->sample_time_s_*c_theta*c_phi;

	//13th row
	this->state_noise_jacobian_(12, 12) = 1;
	//14th row
	this->state_noise_jacobian_(13, 13) = 1;
	//15th row
	this->state_noise_jacobian_(14, 14) = 1;
}

template <typename T>
void Ekf15Dof<T>::GetMeas(MatrixInv<T> meas_sensor_val){
	/* First 3 indices should be magnetometer x, y, z readings in body frame,
	next 3 indices should be GPS measured position translated into NED frame and last
	3 indices should be GPS measured NED velocity*/

	//Compute heading from magnetometer
	T m1 = meas_sensor_val(0);
	T m2 = meas_sensor_val(1);
	T m3 = meas_sensor_val(2);

	//Assume we are in bay area and use a declination 
	T magnetic_declination = 0;//13.01/57.29578;
	this->computed_meas_(0) = atan2( -m2*c_phi + m3*s_phi, m1*c_theta + (m2*s_phi + m3*c_phi)*s_theta ) + magnetic_declination;
	this->computed_meas_(1) = meas_sensor_val(3);//Px
	this->computed_meas_(2) = meas_sensor_val(4);//Py
	this->computed_meas_(3) = meas_sensor_val(5);//Pz
	this->computed_meas_(4) = meas_sensor_val(6);//Vx
	this->computed_meas_(5) = meas_sensor_val(7);//Vy
	this->computed_meas_(6) = meas_sensor_val(8);//Vz
}

template <typename T>
void Ekf15Dof<T>::ComputeMeasJacobian(MatrixInv<T> meas_sensor_val){
	// Do nothing
}

template <typename T>
void Ekf15Dof<T>::ComputeMeasNoiseJacobian(MatrixInv<T> meas_sensor_val){
	// Do Nothing
}

template <typename T>
void Ekf15Dof<T>::ComputeMeasFromState(MatrixInv<T> time_propagated_state){
	// In this case it is y = H*x
	this->meas_from_propogated_state_ = this->meas_jacobian_*this->time_propagated_state_;
}

// Explicit template instantiation
template class Ekf15Dof<float>;
template class Ekf15Dof<double>;