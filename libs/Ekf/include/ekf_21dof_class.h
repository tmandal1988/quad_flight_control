#ifndef EKF21DOF_H
#define EKF21DOF_H
#endif

#include<iostream>
#include<string>
#include <functional>

#include <ekf_base_class.h>

using namespace std;

#define NUM_STATES 			21
#define NUM_MEAS 			7

/* First three indices correspond to body angular rates and last 3 sensors corresponds to
accelerometers measurement without the gravity resolved in body frame*/
#define NUM_STATES_SENSOR 	6 

template <typename T>
class Ekf21Dof:public EkfBase<T>{
	public:
		// constructors
		Ekf21Dof(T sample_time_s, MatrixInv<T> initial_state, MatrixInv<T> process_noise_q, MatrixInv<T> meas_noise_r);

		// destructor
		~Ekf21Dof();	

	protected:
		// member functions
		void PropagateState(MatrixInv<T> previous_state, MatrixInv<T> state_sensor_val);
		void ComputeStateJacobian(MatrixInv<T> previous_state, MatrixInv<T> state_sensor_val);
		void ComputeStateNoiseJacobian(MatrixInv<T> previous_state);
		void GetMeas(MatrixInv<T> meas_sensor_val);
		void ComputeMeasJacobian(MatrixInv<T> meas_sensor_val);
		void ComputeMeasNoiseJacobian(MatrixInv<T> meas_sensor_val);
		void ComputeMeasFromState(MatrixInv<T> time_propagated_state);

		void ComputeTrignometricValues();

		// trignometric values
		T s_phi;
		T c_phi;

		T s_theta;
		T c_theta;
		T t_theta;
		T sc_theta;

		T s_psi;
		T c_psi;

			
};