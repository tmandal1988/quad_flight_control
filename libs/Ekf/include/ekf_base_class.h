#ifndef EKFBASE_H
#define EKFBASE_H
#endif

#include<iostream>
#include<string>
#include<functional>
#include<cmath>

#include <Matrix/matrix_inv_class.h>

using namespace std;

template <typename T>
class EkfBase{
	
	public:
		// constructors
		EkfBase(size_t num_states, size_t num_meas, size_t num_states_sensor, T sample_time_s, 
			MatrixInv<T> initial_state, MatrixInv<T> process_noise_q_, MatrixInv<T> meas_noise_r_);

		// destructor
		~EkfBase();

		void run(MatrixInv<T> state_sensor_val, MatrixInv<T> meas_sensor_val);
	protected:

		// EKF dimension and sample time
		const size_t num_states_;
		const size_t num_meas_;
		const size_t num_states_sensor_;
		const T sample_time_s_;


		// member functions
		virtual void PropagateState(MatrixInv<T> previous_state, MatrixInv<T> state_sensor_val) = 0;
		virtual void ComputeStateJacobian(MatrixInv<T> previous_state, MatrixInv<T> state_sensor_val) = 0;
		virtual void ComputeStateNoiseJacobian(MatrixInv<T> previous_state) = 0;
		virtual void GetMeas(MatrixInv<T> meas_sensor_val) = 0;
		virtual void ComputeMeasJacobian(MatrixInv<T> meas_sensor_val) = 0;
		virtual void ComputeMeasNoiseJacobian(MatrixInv<T> meas_sensor_val) = 0;
		virtual void ComputeMeasFromState(MatrixInv<T> time_propagated_state) = 0;

		

		// EKF variables
		MatrixInv<T> initial_state_;
		MatrixInv<T> current_state_;
		MatrixInv<T> computed_meas_;

		const MatrixInv<T> process_noise_q_;
		const MatrixInv<T> meas_noise_r_;
		MatrixInv<T> covariance_p_;

		// Variables to run EKF
		MatrixInv<T> time_propagated_state_;
		MatrixInv<T> meas_from_propogated_state_;
		MatrixInv<T> state_jacobian_;
		MatrixInv<T> state_noise_jacobian_;
		MatrixInv<T> meas_jacobian_;
		MatrixInv<T> meas_noise_jacobian_;
		MatrixInv<T> inv_part_of_kalman_gain_;
		MatrixInv<T> kalman_gain_;

		MatrixInv<T> kalman_eye_;
};