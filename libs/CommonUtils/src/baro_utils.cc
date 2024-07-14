#include "baro_utils.h"

//Constructor
BaroHelper::BaroHelper(size_t init_num_samples):
init_num_samples_(init_num_samples),
is_baro_ready_(false),
use_comp_filter_(false),
use_kalman_filter_(false){	
	// Variable to indicate Baro Thread to stop
	stop_baro_read_thread_.store(false);

	barometer_.initialize();
}

void BaroHelper::StartBaroReader(int cpu_to_use, int32_t priority, float sample_time_s){
	sample_time_s_ = sample_time_s;

	SetInitHeight();
	agl_est_m_ = 0.0;
	climb_rate_est_mps_ = 0.0;
	veh_ned_az_est_mps2_ = 0.0;

	baro_reader_thread_ = thread(&BaroHelper::BaroReadLoop, this);

    CPU_ZERO(&cpuset_);
    CPU_SET(cpu_to_use, &cpuset_);

	pthread_getschedparam(baro_reader_thread_.native_handle(), &policy_, &sch_);
	sch_.sched_priority = priority;
	pthread_setschedparam(baro_reader_thread_.native_handle(), SCHED_FIFO, &sch_);
	int rc = pthread_setaffinity_np(baro_reader_thread_.native_handle(),
                                    sizeof(cpuset_), &cpuset_);    
	if (rc != 0) {
      std::cerr << "Error calling pthread_setaffinity_np on Baro Reader Thread: " << rc << "\n";
    }
}

void BaroHelper::BaroReadLoop(){
	if(use_kalman_filter_){
		BaroReadLoopKalmanFilter();
	}else{
		BaroReadLoopCompFilter();
	}
}

void BaroHelper::BaroReadLoopCompFilter(){
	{
        unique_lock<mutex> baro_data_lock(baro_data_mutex_);
		is_baro_ready_ = true;
	}
	float baro_agl = 0;
	float veh_ned_az_mps2 = 0;
	float eh_term = 0;

	float accel_eps = 1e-10;
	
	vector<float> zupt_array(zupt_length_, 0.0);
	size_t z_idx = 0;

	// Loop timers
	chrono::high_resolution_clock::time_point loop_start;
	chrono::high_resolution_clock::time_point loop_end;

	
	//Sample Step
	size_t dt_count = sample_time_s_*1e6;
	// initialize the duration
	chrono::microseconds delta (dt_count); 
	auto duration = chrono::duration_cast<chrono::microseconds> (delta);
	auto duration_count = duration.count();

	// Keep reading the Baro data forever
	while(1){
		// Get loop start time
		loop_start = chrono::high_resolution_clock::now();

		barometer_.refreshPressure();
		usleep(2500); // Waiting for pressure data ready
		barometer_.readPressure();

		barometer_.refreshTemperature();
        usleep(2500); // Waiting for temperature data ready
        barometer_.readTemperature();

        // {
        //  	unique_lock<mutex> baro_data_lock(baro_data_mutex_);
    	barometer_.calculatePressureAndTemperature();
        // }

        if (!( (abs(body_accel_mps2_[0]) < accel_eps) && (abs(body_accel_mps2_[1]) < accel_eps) &&
        	(abs(body_accel_mps2_[2]) < accel_eps) )){

	        // Run altitude complimentary filter
	        baro_agl = ComputeBaroAgl();
	        veh_ned_az_mps2 = ComputeAzInNedWithoutGravity();
	        eh_term = baro_agl - agl_est_m_;

	        float cest_tmp = sample_time_s_*eh_term;
	        float vz_tmp = sample_time_s_*veh_ned_az_mps2;
	        // {
	        //  	unique_lock<mutex> baro_data_lock(baro_data_mutex_);
        	agl_est_m_ += sample_time_s_*climb_rate_est_mps_ + comp_kc_[0]*cest_tmp + 
        			  0.5*sample_time_s_*comp_kc_[1]*cest_tmp + 0.5*sample_time_s_*vz_tmp;
        	climb_rate_est_mps_ += comp_kc_[1]*cest_tmp + vz_tmp;
	    	// }

	        // Do zero velocity update
	        zupt_array[z_idx] = abs(veh_ned_az_mps2);
	        size_t next_z_idx = (z_idx + 1) % zupt_length_;
	        z_idx = next_z_idx;        

	        size_t idx = 0;
	        for(idx = 0; idx < zupt_length_; idx++){
	        	if(zupt_array[idx] > zupt_threshold_)
	        		break;
	        }

	        if(idx == zupt_length_){
	        	// {
	         	// 	unique_lock<mutex> baro_data_lock(baro_data_mutex_);
	        	climb_rate_est_mps_ = 0.0;
	        	// }
	        }
    	}

		if(stop_baro_read_thread_.load()){
				break;
		}
		baro_debug_data_[0] = (float)duration_count;
		baro_debug_data_[1] = euler_angles_rad_[0];
		baro_debug_data_[2] = euler_angles_rad_[1];
		baro_debug_data_[3] = body_accel_mps2_[0];
		baro_debug_data_[4] = body_accel_mps2_[1];
		baro_debug_data_[5] = body_accel_mps2_[2];
		baro_debug_data_[6] = agl_est_m_;
		baro_debug_data_[7] = climb_rate_est_mps_;

		// Get the stop time and compute the duration
    	loop_end = std::chrono::high_resolution_clock::now();

    	duration_count = chrono::duration_cast<chrono::microseconds>(loop_end - loop_start).count();
    	while(duration_count < dt_count){
    		duration_count = chrono::duration_cast<chrono::microseconds>(std::chrono::high_resolution_clock::now() - loop_start).count();
    	}
	}
}

void BaroHelper::BaroReadLoopKalmanFilter(){
	{
        unique_lock<mutex> baro_data_lock(baro_data_mutex_);
		is_baro_ready_ = true;
	}
	float baro_agl = 0;
	float veh_ned_az_mps2 = 0;

	float accel_eps = 1e-10;
	
	vector<float> zupt_array(zupt_length_, 0.0);
	size_t z_idx = 0;

	// Loop timers
	chrono::high_resolution_clock::time_point loop_start;
	chrono::high_resolution_clock::time_point loop_end;

	
	//Sample Step
	size_t dt_count = sample_time_s_*1e6;
	// initialize the duration
	chrono::microseconds delta (dt_count); 
	auto duration = chrono::duration_cast<chrono::microseconds> (delta);
	auto duration_count = duration.count();

	array< array<float, 3>, 3> p_cov_tmp = {0};
	float agl_est_tmp_m = 0;
	float veh_ned_az_est_tmp_mps2 = 0;
	float climb_rate_est_tmp_mps = 0;

	// Keep reading the Baro data forever
	while(1){
		// Get loop start time
		loop_start = chrono::high_resolution_clock::now();

		barometer_.refreshPressure();
		usleep(2500); // Waiting for pressure data ready
		barometer_.readPressure();

		barometer_.refreshTemperature();
        usleep(2500); // Waiting for temperature data ready
        barometer_.readTemperature();

        // {
        //  	unique_lock<mutex> baro_data_lock(baro_data_mutex_);
        barometer_.calculatePressureAndTemperature();
        // }

        if (!( (abs(body_accel_mps2_[0]) < accel_eps) && (abs(body_accel_mps2_[1]) < accel_eps) &&
        	(abs(body_accel_mps2_[2]) < accel_eps) )){
        	float sample_time_s_sq = sample_time_s_*sample_time_s_;
	        // Compute baro alt and vertical acceleration minus gravity
	        baro_agl = ComputeBaroAgl();
	        veh_ned_az_mps2 = ComputeAzInNedWithoutGravity();

	        // propogate states, accel state uses constant dynamics
	        agl_est_m_ += (veh_ned_az_est_mps2_* sample_time_s_sq)*0.5 + climb_rate_est_mps_*sample_time_s_;
			climb_rate_est_mps_ += sample_time_s_*veh_ned_az_est_mps2_;

			// Update the covariance matrix
			// 1 row
			p_cov_tmp[0][0] = p_cov_[0][0] + proc_noise_[0][0] + p_cov_[2][0]*sample_time_s_ + sample_time_s_*(p_cov_[0][2] + 
							   p_cov_[2][2]*sample_time_s_ + (p_cov_[1][2]*sample_time_s_sq)*0.5) + (p_cov_[1][0]*sample_time_s_sq)*0.5 + 
							   (sample_time_s_sq*(p_cov_[0][1] + p_cov_[2][1]*sample_time_s_ + (p_cov_[1][1]*sample_time_s_sq)*0.5))*0.5;

			p_cov_tmp[0][1] = p_cov_[0][1] + p_cov_[2][1]*sample_time_s_ + (p_cov_[1][1]*sample_time_s_sq)*0.5;

			p_cov_tmp[0][2] = p_cov_[0][2] + p_cov_[2][2]*sample_time_s_ + sample_time_s_*(p_cov_[0][1] + p_cov_[2][1]*sample_time_s_ + 
							   (p_cov_[1][1]*sample_time_s_sq)*0.5) + (p_cov_[1][2]*sample_time_s_sq)*0.5;
			// 2 row
			p_cov_tmp[1][0] = p_cov_[1][0] + p_cov_[1][2]*sample_time_s_ + (p_cov_[1][1]*sample_time_s_sq)*0.5;
			p_cov_tmp[1][1] = p_cov_[1][1] + proc_noise_[1][1];
			p_cov_tmp[1][2] = p_cov_[1][2] + p_cov_[1][1]*sample_time_s_;
			// 3 row
			p_cov_tmp[2][0] = p_cov_[2][0] + (sample_time_s_sq*(p_cov_[2][1] + p_cov_[1][1]*sample_time_s_))*0.5 + p_cov_[1][0]*sample_time_s_ + 
							   sample_time_s_*(p_cov_[2][2] + p_cov_[1][2]*sample_time_s_);
			p_cov_tmp[2][1] = p_cov_[2][1] + p_cov_[1][1]*sample_time_s_;
			p_cov_tmp[2][2] = p_cov_[2][2] + proc_noise_[2][2] + p_cov_[1][2]*sample_time_s_ + sample_time_s_*(p_cov_[2][1] + p_cov_[1][1]*sample_time_s_);

			// Compute Kalman Gain
			// 1 row
			k_gain_[0][0] = (p_cov_tmp[0][0]*meas_noise_[1][1] + p_cov_tmp[0][0]*p_cov_tmp[1][1] - p_cov_tmp[0][1]*p_cov_tmp[1][0])/
							(p_cov_tmp[0][0]*meas_noise_[1][1] + p_cov_tmp[1][1]*meas_noise_[0][0] + meas_noise_[0][0]*meas_noise_[1][1] + 
							p_cov_tmp[0][0]*p_cov_tmp[1][1] - p_cov_tmp[0][1]*p_cov_tmp[1][0]);

			k_gain_[0][1] = (p_cov_tmp[0][1]*meas_noise_[0][0])/(p_cov_tmp[0][0]*meas_noise_[1][1] + p_cov_tmp[1][1]*meas_noise_[0][0] + 
							meas_noise_[0][0]*meas_noise_[1][1] + p_cov_tmp[0][0]*p_cov_tmp[1][1] - p_cov_tmp[0][1]*p_cov_tmp[1][0]);
			// 2 row
			k_gain_[1][0] = (p_cov_tmp[1][0]*meas_noise_[1][1])/(p_cov_tmp[0][0]*meas_noise_[1][1] + p_cov_tmp[1][1]*meas_noise_[0][0] + 
							meas_noise_[0][0]*meas_noise_[1][1] + p_cov_tmp[0][0]*p_cov_tmp[1][1] - p_cov_tmp[0][1]*p_cov_tmp[1][0]);

			k_gain_[1][1] = (p_cov_tmp[1][1]*meas_noise_[0][0] + p_cov_tmp[0][0]*p_cov_tmp[1][1] - p_cov_tmp[0][1]*p_cov_tmp[1][0])/
							(p_cov_tmp[0][0]*meas_noise_[1][1] + p_cov_tmp[1][1]*meas_noise_[0][0] + meas_noise_[0][0]*meas_noise_[1][1] + 
							p_cov_tmp[0][0]*p_cov_tmp[1][1] - p_cov_tmp[0][1]*p_cov_tmp[1][0]);
			// 3 row
			k_gain_[2][0] = (p_cov_tmp[2][0]*meas_noise_[1][1] - p_cov_tmp[1][0]*p_cov_tmp[2][1] + p_cov_tmp[1][1]*p_cov_tmp[2][0])/
							(p_cov_tmp[0][0]*meas_noise_[1][1] + p_cov_tmp[1][1]*meas_noise_[0][0] + meas_noise_[0][0]*meas_noise_[1][1] +
							p_cov_tmp[0][0]*p_cov_tmp[1][1] - p_cov_tmp[0][1]*p_cov_tmp[1][0]);

			k_gain_[2][1] = (p_cov_tmp[2][1]*meas_noise_[0][0] + p_cov_tmp[0][0]*p_cov_tmp[2][1] - p_cov_tmp[0][1]*p_cov_tmp[2][0])/
							(p_cov_tmp[0][0]*meas_noise_[1][1] + p_cov_tmp[1][1]*meas_noise_[0][0] + meas_noise_[0][0]*meas_noise_[1][1] + 
							p_cov_tmp[0][0]*p_cov_tmp[1][1] - p_cov_tmp[0][1]*p_cov_tmp[1][0]);

			// Update states
			// 1 row
			agl_est_tmp_m = agl_est_m_ + k_gain_[0][1]*(veh_ned_az_mps2 - veh_ned_az_est_mps2_) - k_gain_[0][0]*(agl_est_m_ - baro_agl);
			// 2 row
			veh_ned_az_est_tmp_mps2 = veh_ned_az_est_mps2_ + k_gain_[1][1]*(veh_ned_az_mps2 - veh_ned_az_est_mps2_) - k_gain_[1][0]*(agl_est_m_ - baro_agl);
			// 3 row
			climb_rate_est_tmp_mps = climb_rate_est_mps_ + k_gain_[2][1]*(veh_ned_az_mps2 - veh_ned_az_est_mps2_) - k_gain_[2][0]*(agl_est_m_ - baro_agl);

			agl_est_m_ = agl_est_tmp_m;
			veh_ned_az_est_mps2_ = veh_ned_az_est_tmp_mps2;
			climb_rate_est_mps_ = climb_rate_est_tmp_mps;

			// Update covariance Matrix
			// 1 row

			p_cov_[0][0] = - p_cov_tmp[0][0]*(k_gain_[0][0] - 1) - k_gain_[0][1]*p_cov_tmp[1][0];
			p_cov_[0][1] = - p_cov_tmp[0][1]*(k_gain_[0][0] - 1) - k_gain_[0][1]*p_cov_tmp[1][1];
			p_cov_[0][2] = - p_cov_tmp[0][2]*(k_gain_[0][0] - 1) - k_gain_[0][1]*p_cov_tmp[1][2];
			// 2 row
			p_cov_[1][0] = - p_cov_tmp[1][0]*(k_gain_[1][1] - 1) - k_gain_[1][0]*p_cov_tmp[0][0];
			p_cov_[1][1] = - p_cov_tmp[1][1]*(k_gain_[1][1] - 1) - k_gain_[1][0]*p_cov_tmp[0][1];
			p_cov_[1][2] = - p_cov_tmp[1][2]*(k_gain_[1][1] - 1) - k_gain_[1][0]*p_cov_tmp[0][2];
			// 3 row
			p_cov_[2][0] = p_cov_tmp[2][0] - k_gain_[2][0]*p_cov_tmp[0][0] - k_gain_[2][1]*p_cov_tmp[1][0];
			p_cov_[2][1] = p_cov_tmp[2][1] - k_gain_[2][0]*p_cov_tmp[0][1] - k_gain_[2][1]*p_cov_tmp[1][1];
			p_cov_[2][2] = p_cov_tmp[2][2] - k_gain_[2][0]*p_cov_tmp[0][2] - k_gain_[2][1]*p_cov_tmp[1][2];

			// If GPS Alt is available then do sequential update on it
			if(gps_alt_new){
				gps_alt_new = false;
				// Compute Kalman Gain
				// 1st row
				k_gain_gps_[0][0] = p_cov_[0][0]/(p_cov_[0][0] + meas_noise_[2][2]);
				// 2nd row
				k_gain_gps_[0][1] = p_cov_[1][0]/(p_cov_[0][0] + meas_noise_[2][2]);
				// 3rd row
				k_gain_gps_[0][2] = p_cov_[2][0]/(p_cov_[0][0] + meas_noise_[2][2]);

				// update states
				// agl
				agl_est_m_ = agl_est_m_ - k_gain_gps_[0][0]*(agl_est_m_ - gps_alt_m);
				// ned az
				veh_ned_az_est_mps2_ = veh_ned_az_est_mps2_ - k_gain_gps_[0][1]*(agl_est_m_ - gps_alt_m);
				// cr
				climb_rate_est_mps_ = climb_rate_est_mps_ - k_gain_gps_[0][2]*(agl_est_m_ - gps_alt_m);

				//Update covariance
				// 1 row
				p_cov_tmp[0][0] = p_cov_[0][0];
				p_cov_[0][0] = -p_cov_tmp[0][0]*(p_cov_tmp[0][0]/(p_cov_tmp[0][0] + meas_noise_[2][2]) - 1);

				p_cov_tmp[0][1] = p_cov_[0][1];
				p_cov_[0][1] = -p_cov_tmp[0][1]*(p_cov_tmp[0][0]/(p_cov_tmp[0][0] + meas_noise_[2][2]) - 1);

				p_cov_tmp[0][2] = p_cov_[0][2];
				p_cov_[0][2] = -p_cov_tmp[0][2]*(p_cov_tmp[0][0]/(p_cov_tmp[0][0] + meas_noise_[2][2]) - 1);

				// 2 row
				p_cov_tmp[1][0] = p_cov_[1][0];
				p_cov_[1][0] = p_cov_tmp[1][0] - (p_cov_tmp[0][0]*p_cov_tmp[1][0])/(p_cov_tmp[0][0] + meas_noise_[2][2]);

				p_cov_tmp[1][1] = p_cov_[1][1];
				p_cov_[1][1] = p_cov_tmp[1][1] - (p_cov_tmp[0][1]*p_cov_tmp[1][0])/(p_cov_tmp[0][0] + meas_noise_[2][2]);

				p_cov_tmp[1][2] = p_cov_[1][2];
				p_cov_[1][2] = p_cov_tmp[1][2] - (p_cov_tmp[0][2]*p_cov_tmp[1][0])/(p_cov_tmp[0][0] + meas_noise_[2][2]);

				// 3 row
				p_cov_tmp[2][0] = p_cov_[2][0];
				p_cov_[2][0] = p_cov_tmp[2][0] - (p_cov_tmp[0][0]*p_cov_tmp[2][0])/(p_cov_tmp[0][0] + meas_noise_[2][2]);

				p_cov_tmp[2][1] = p_cov_[2][1];
				p_cov_[2][1] = p_cov_tmp[2][1] - (p_cov_tmp[0][1]*p_cov_tmp[2][0])/(p_cov_tmp[0][0] + meas_noise_[2][2]);

				p_cov_tmp[2][2] = p_cov_[2][2];
				p_cov_[2][2] = p_cov_tmp[2][2] - (p_cov_tmp[0][2]*p_cov_tmp[2][0])/(p_cov_tmp[0][0] + meas_noise_[2][2]);

			}

			// // If GPS cr is available then do sequential update on it
			if(gps_cr_new){
				gps_cr_new = false;
				// 1 row
				k_gain_gps_[1][0] = p_cov_[0][2]/(p_cov_[2][2] + meas_noise_[3][3]);
				// 2 row
				k_gain_gps_[1][1] = p_cov_[1][2]/(p_cov_[2][2] + meas_noise_[3][3]);
				// 3 row
				k_gain_gps_[1][2] = p_cov_[2][2]/(p_cov_[2][2] + meas_noise_[3][3]);

				// update states
				// agl
				agl_est_m_ = agl_est_m_ - k_gain_gps_[1][0]*(climb_rate_est_mps_ - gps_cr_mps);
				// ned az
				veh_ned_az_est_tmp_mps2 = veh_ned_az_est_tmp_mps2 - k_gain_gps_[1][1]*(climb_rate_est_mps_ - gps_cr_mps);
				// cr
				climb_rate_est_mps_ = climb_rate_est_mps_ - k_gain_gps_[1][2]*(climb_rate_est_mps_ - gps_cr_mps);

				// 1 row
				p_cov_tmp[0][0] = p_cov_[0][0];
				p_cov_[0][0] = p_cov_tmp[0][0] - (p_cov_tmp[0][2]*p_cov_tmp[2][0])/(p_cov_tmp[2][2] + meas_noise_[3][3]);

				p_cov_tmp[0][1] = p_cov_[0][1];
				p_cov_[0][1] = p_cov_tmp[0][1] - (p_cov_tmp[0][2]*p_cov_tmp[2][1])/(p_cov_tmp[2][2] + meas_noise_[3][3]);

				p_cov_tmp[0][2] = p_cov_[0][2];
				p_cov_[0][2] = p_cov_tmp[0][2] - (p_cov_tmp[0][2]*p_cov_tmp[2][2])/(p_cov_tmp[2][2] + meas_noise_[3][3]);

				// 2 row
				p_cov_tmp[1][0] = p_cov_[1][0];
				p_cov_[1][0] = p_cov_tmp[1][0] - (p_cov_tmp[1][2]*p_cov_tmp[2][0])/(p_cov_tmp[2][2] + meas_noise_[3][3]);

				p_cov_tmp[1][1] = p_cov_[1][1];
				p_cov_[1][1] = p_cov_tmp[1][1] - (p_cov_tmp[1][2]*p_cov_tmp[2][1])/(p_cov_tmp[2][2] + meas_noise_[3][3]);

				p_cov_tmp[1][2] = p_cov_[1][2];
				p_cov_[1][2] = p_cov_tmp[1][2] - (p_cov_tmp[1][2]*p_cov_tmp[2][2])/(p_cov_tmp[2][2] + meas_noise_[3][3]);
				// 3 row
				p_cov_tmp[2][0] = p_cov_[2][0];
				p_cov_[2][0] = -p_cov_tmp[2][0]*(p_cov_tmp[2][2]/(p_cov_tmp[2][2] + meas_noise_[3][3]) - 1);

				p_cov_tmp[2][1] = p_cov_[2][1];
				p_cov_[2][1] = -p_cov_tmp[2][1]*(p_cov_tmp[2][2]/(p_cov_tmp[2][2] + meas_noise_[3][3]) - 1);

				p_cov_tmp[2][2] = p_cov_[2][2];
				p_cov_[2][2] = -p_cov_tmp[2][2]*(p_cov_tmp[2][2]/(p_cov_tmp[2][2] + meas_noise_[3][3]) - 1);
			}

	        // Do zero velocity update
	        zupt_array[z_idx] = abs(veh_ned_az_mps2);
	        size_t next_z_idx = (z_idx + 1) % zupt_length_;
	        z_idx = next_z_idx;        

	        size_t idx = 0;
	        for(idx = 0; idx < zupt_length_; idx++){
	        	if(zupt_array[idx] > zupt_threshold_)
	        		break;
	        }

	        if(idx == zupt_length_){
	        	// {
	         	// 	unique_lock<mutex> baro_data_lock(baro_data_mutex_);
        		climb_rate_est_mps_ = 0.0;
	        	// }
	        }
    	}

		if(stop_baro_read_thread_.load()){
				break;
		}
		baro_debug_data_[0] = (float)duration_count;
		baro_debug_data_[1] = euler_angles_rad_[0];
		baro_debug_data_[2] = euler_angles_rad_[1];
		baro_debug_data_[3] = body_accel_mps2_[0];
		baro_debug_data_[4] = body_accel_mps2_[1];
		baro_debug_data_[5] = body_accel_mps2_[2];
		baro_debug_data_[6] = agl_est_m_;
		baro_debug_data_[7] = climb_rate_est_mps_;
		baro_debug_data_[8] = veh_ned_az_est_mps2_;

		// Get the stop time and compute the duration
    	loop_end = std::chrono::high_resolution_clock::now();

    	duration_count = chrono::duration_cast<chrono::microseconds>(loop_end - loop_start).count();
    	while(duration_count < dt_count){
    		duration_count = chrono::duration_cast<chrono::microseconds>(std::chrono::high_resolution_clock::now() - loop_start).count();
    	}
	}
}

void BaroHelper::SetInitHeight(){
	float press_sum = 0;
	for(size_t idx = 0; idx < init_num_samples_; idx++){
		barometer_.refreshPressure();
		usleep(3500); // Waiting for pressure data ready
		barometer_.readPressure();

		barometer_.refreshTemperature();
        usleep(3500); // Waiting for temperature data ready
        barometer_.readTemperature();

        barometer_.calculatePressureAndTemperature();
        //100 Hz update
		usleep(3000);

		press_sum += barometer_.getPressure();
	}

	float mean_init_press_pa = press_sum*100.0/init_num_samples_;
	h_init_m_ = 44330.0*( 1 - pow( (mean_init_press_pa/SEA_LEVEL_PRESS_PA),0.19 ) );
}

void BaroHelper::SetGpsInitStatus(const bool gps_init_status){
	//If GPS is initialized then initial alt of NED GPS will be set to 0,
	//therefore intial height estimated by baro will be the offset between,
	//baro alt and gps alt
	if(gps_init_status){
		gps_initialized = true;
		baro_alt_gps_alt_offset = h_init_m_;
	}
}

float BaroHelper::ComputeBaroAgl(){
	return 44330.0*( 1 - pow( (barometer_.getPressure()*100.0/SEA_LEVEL_PRESS_PA),0.19 ) ) - h_init_m_;
}

float BaroHelper::ComputeAzInNedWithoutGravity(){
	float c_theta = cos(euler_angles_rad_[1]);
	// Convert body accel to NED accel
	float total_accel_ned_z = -sin(euler_angles_rad_[1])*body_accel_mps2_[0] + 
	sin(euler_angles_rad_[0])*c_theta*body_accel_mps2_[1] + 
	cos(euler_angles_rad_[0])*c_theta*body_accel_mps2_[2];

	// Remove gravity from total accel, -ve sign to make altitude increasing upwards
	return -(total_accel_ned_z + G_SI);
}

void BaroHelper::SetComplimentaryFilterParams(const float accel_sigma, const float baro_sigma, 
										  const size_t zupt_length, const float zupt_threshold){
	ComputeComplimentaryFilterGain(accel_sigma, baro_sigma);
	zupt_length_ = zupt_length;
	zupt_threshold_ = zupt_threshold;
	use_comp_filter_ = true;
}

void BaroHelper::SetKalmanFilterParams(const array<float, 3> &p_cov_init, const array<float, 3> &proc_noise, 
										  const array<float, 4> &meas_noise, const size_t zupt_length, const float zupt_threshold){
	p_cov_[0][0] = p_cov_init[0];
	p_cov_[1][1] = p_cov_init[1];
	p_cov_[2][2] = p_cov_init[2];

	proc_noise_[0][0] = proc_noise[0];
	proc_noise_[1][1] = proc_noise[1];
	proc_noise_[2][2] = proc_noise[2];

	meas_noise_[0][0] = meas_noise[0];
	meas_noise_[1][1] = meas_noise[1];
	meas_noise_[2][2] = meas_noise[2];
	meas_noise_[3][3] = meas_noise[3];


	zupt_length_ = zupt_length;
	zupt_threshold_ = zupt_threshold;
	use_kalman_filter_ = true;
}


void BaroHelper::ComputeComplimentaryFilterGain(const float sigma_accel, const float sigma_baro){
	comp_kc_[0] = sqrt(2 * sigma_accel/sigma_baro);
	comp_kc_[1] = sigma_accel/sigma_baro;
}

void BaroHelper::SetBodyAccels(const array<float, 3> &body_accel_mps2){
	 body_accel_mps2_ = body_accel_mps2;
}

void BaroHelper::SetGpsVelAndAlt(const array<float, 2> &baro_gps_alt_and_climb_rate_mps, 
	const array<bool, 2> & baro_gps_alt_and_climb_rate_flag){

	if(gps_initialized){
		gps_cr_mps = baro_gps_alt_and_climb_rate_mps[1];
		// Add offset between baro and gps
		gps_alt_m = baro_gps_alt_and_climb_rate_mps[0] + baro_alt_gps_alt_offset;
		gps_cr_new = baro_gps_alt_and_climb_rate_flag[1];
		gps_alt_new = baro_gps_alt_and_climb_rate_flag[0];
	}
}

void BaroHelper::SetEulerAngles(const array<float, 3> &euler_angles_rad){
	euler_angles_rad_ = euler_angles_rad;
}

void BaroHelper::GetBaroPressAndTemp(float baro_data[]){
	baro_data[0] = barometer_.getPressure();
	baro_data[1] = barometer_.getTemperature();

}

void BaroHelper::GetAglAndClimbRateEst(float baro_data[]){
	baro_data[0] = agl_est_m_;
	baro_data[1] = climb_rate_est_mps_;
}

void BaroHelper::GetBaroDebugData(float baro_debug_data[]){
	for(size_t idx = 0; idx < 9; idx++){
		baro_debug_data[idx] = baro_debug_data_[idx];
	}
}

BaroHelper::~BaroHelper(){
	// Making sure that Baro read thread is stopped
	stop_baro_read_thread_.store(true);
	if(baro_reader_thread_.joinable())
		baro_reader_thread_.join();
}
