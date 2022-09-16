#ifndef EKFBASE_H
#define EKFBASE_H
#endif

#include<iostream>
#include<string>
#include <functional>

#include <Matrix/matrix_inv_class.h>

using namespace std;

template <typename T>
class EkfBase{
	
	public:
		// constructors
		EkfBase(size_t num_states, size_t num_meas, size_t num_states_sensor, T sample_time_s, 
			MatrixBase<T> initial_state, MatrixBase<T> process_noise_q_, MatrixBase<T> meas_noise_r_);

		// destructor
		~EkfBase();

		void run();

	protected:

		// EKF dimension and sample time
		const size_t num_states_;
		const size_t num_meas_;
		const size_t num_states_sensor_;
		const T sample_time_s_;


		// member functions
		virtual MatrixBase<T> GetStateFunction(MatrixBase<T> previous_state, MatrixBase<T> current_states_meas) = 0;
		virtual MatrixBase<T> GetStateJacobian(MatrixBase<T> previous_state, MatrixBase<T> current_states_meas) = 0;
		virtual MatrixBase<T> GetStateNoiseJacobian(MatrixBase<T> previous_state, MatrixBase<T> current_states_meas) = 0;
		virtual MatrixBase<T> GetMeasFunction(MatrixBase<T> current_meas) = 0;
		virtual MatrixBase<T> GetMeasJacobian(MatrixBase<T> current_meas) = 0;
		virtual MatrixBase<T> GetMeasNoiseJacobian(MatrixBase<T> current_meas) = 0;

		// EKF variables
		MatrixBase<T> initial_state_;
		MatrixBase<T> current_state_;
		MatrixBase<T> current_state_sensor_val_;
		MatrixBase<T> current_meas_;

		const MatrixBase<T> process_noise_q_;
		const MatrixBase<T> meas_noise_r_;
		MatrixBase<T> covariance_p_;


		
};