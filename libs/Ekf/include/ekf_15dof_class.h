#ifndef EKF15DOF_H
#define EKF15DOF_H


#include<iostream>
#include<string>
#include <functional>

#include <ekf_base_class.h>

using namespace std;

#define NUM_STATES 			15
#define NUM_MEAS 			7

/* First three indices correspond to body angular rates and last 3 sensors corresponds to
accelerometers measurement without the gravity resolved in body frame*/
#define NUM_STATES_SENSOR 	6 

#define PI 					3.14159265358979
#define PIx2 				6.28318530717959
#define RAD2DEG			    57.2957795130823
#define DEG2RAD 			0.0174532925199433

template <typename T>
class Ekf15Dof:public EkfBase<T>{
	public:
		// constructors
		Ekf15Dof(T sample_time_s, MatrixInv<T> initial_state, MatrixInv<T> process_noise_q, MatrixInv<T> meas_noise_r);

		// destructor
		~Ekf15Dof();	

	protected:
		// member functions
		void PropagateState(MatrixInv<T> state_sensor_val);
		void ComputeStateJacobian(MatrixInv<T> state_sensor_val);
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

		MatrixInv<T> g;

			
};
#endif