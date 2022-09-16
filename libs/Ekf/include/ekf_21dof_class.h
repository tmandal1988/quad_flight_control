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
		MatrixBase<T> GetStateFunction(MatrixBase<T> previous_state, MatrixBase<T> current_states_meas);
		MatrixBase<T> GetStateJacobian(MatrixBase<T> previous_state, MatrixBase<T> current_states_meas);
		MatrixBase<T> GetStateNoiseJacobian(MatrixBase<T> previous_state, MatrixBase<T> current_states_meas);
		MatrixBase<T> GetMeasFunction(MatrixBase<T> current_meas);
		MatrixBase<T> GetMeasJacobian(MatrixBase<T> current_meas);
		MatrixBase<T> GetMeasNoiseJacobian(MatrixBase<T> current_meas);

			
};