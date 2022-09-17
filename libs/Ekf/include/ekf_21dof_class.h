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
		Ekf21Dof(T sample_time_s, MatrixBase<T> initial_state, MatrixBase<T> process_noise_q, MatrixBase<T> meas_noise_r);

		// destructor
		~Ekf21Dof();	

	protected:
		// member functions
		void PropagateState(MatrixBase<T> previous_state, MatrixBase<T> current_states_meas);
		void ComputeStateJacobian(MatrixBase<T> previous_state, MatrixBase<T> current_states_meas);
		MatrixBase<T> ComputeStateNoiseJacobian(MatrixBase<T> previous_state, MatrixBase<T> current_states_meas);
		MatrixBase<T> GetMeas(MatrixBase<T> current_meas);
		MatrixBase<T> ComputeMeasJacobian(MatrixBase<T> current_meas);
		MatrixBase<T> ComputeMeasNoiseJacobian(MatrixBase<T> current_meas);

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