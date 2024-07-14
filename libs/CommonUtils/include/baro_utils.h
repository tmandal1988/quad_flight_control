#ifndef BAROUTILS_H
#define BAROUTILS_H

#include<Navio/Common/MS5611.h>
#include<Navio/Common/Util.h>
#include<constants.h>
#include<unistd.h>
#include<atomic>
#include<mutex>
#include<thread>
#include<pthread.h>
#include<iostream>
#include<vector>

using namespace std;

class BaroHelper{
	public:
		// Constructors
		BaroHelper(size_t init_num_samples = 50);

		// Destructor
		~BaroHelper();

		void StartBaroReader(int cpu_to_use, int32_t priority, float sample_time_s = 0.008);
		void GetBaroPressAndTemp(float baro_data[]);
		void GetAglAndClimbRateEst(float baro_data[]);
		void GetBaroDebugData(float baro_debug_data[]);
		void SetBodyAccels(const array<float, 3> &body_accel_mps2);
		void SetGpsVelAndAlt(const array<float, 2> &baro_gps_alt_and_climb_rate_mps, 
							 const array<bool, 2> & baro_gps_alt_and_climb_rate_flag);
		void SetGpsInitStatus(const bool gps_init_status);
		void SetEulerAngles(const array<float, 3> &euler_angles_rad);
		void SetComplimentaryFilterParams(const float accel_sigma, const float baro_sigma, 
										  const size_t zupt_length = 12, const float zupt_threshold_ = 0.5);
		void SetKalmanFilterParams(const array<float, 3> &p_cov_init, const array<float, 3> &proc_noise, 
								   const array<float, 4> &meas_noise, const size_t zupt_length = 12, 
								   const float zupt_threshold = 0.5);

		bool is_baro_ready_;

	private:
		// Number of samples for initialization of hinit
		size_t init_num_samples_;
		float h_init_m_;
		// Create a thread
		thread baro_reader_thread_;
		// To set CPU affinity
		cpu_set_t cpuset_;
		sched_param sch_;
		int policy_;

		// Mutex to guard resource access between threads while running GpsReadLoop()
	    mutex baro_data_mutex_;

		MS5611 barometer_;
		void BaroReadLoop();

		// Function to setup base station height
		void SetInitHeight();

		float GetHeightFromPressure();

		// Variables for Complimentary filter
		float sample_time_s_;
		array<float, 3> body_accel_mps2_{0};
		float gps_cr_mps{0};
		float gps_alt_m{0};
		bool gps_alt_new{false};
		bool gps_cr_new{false};
		bool gps_initialized{false};
		float baro_alt_gps_alt_offset{0};
		array<float, 3> euler_angles_rad_{0};
		array<float, 2> comp_kc_{NAN};

		float agl_est_m_;
		float climb_rate_est_mps_;
		size_t zupt_length_;
		float zupt_threshold_;
		bool use_comp_filter_;

		// Methods for altitude complimentary filter;
		void ComputeComplimentaryFilterGain(const float sigma_accel, const float sigma_baro);
		void BaroReadLoopCompFilter();

		// Methods to compute baro agl and vertical acceleration minus gravity
		float ComputeBaroAgl();
		float ComputeAzInNedWithoutGravity();

		float baro_debug_data_[9]{0};

		// Variables for Kalman Filter
		array< array<float, 3>, 3> proc_noise_{0};
		array< array<float, 4>, 4> meas_noise_{0};
		array< array<float, 3>, 3> p_cov_{0};
		array< array<float, 2>, 3> k_gain_{0};
		array< array<float, 2>, 3> k_gain_gps_{0};
		bool use_kalman_filter_;
		float veh_ned_az_est_mps2_;

		void BaroReadLoopKalmanFilter();

		// Variable to indicate Baro Thread to stop
	    atomic<bool> stop_baro_read_thread_;
};


#endif