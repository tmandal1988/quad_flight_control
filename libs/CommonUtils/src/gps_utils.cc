#include "gps_utils.h"

GpsHelper::GpsHelper(size_t n_valid_gps_count, size_t n_gps_meas_count, uint16_t time_in_ms_bw_samples):
n_valid_gps_count_(n_valid_gps_count),
n_gps_meas_count_(n_gps_meas_count),
time_in_ms_bw_samples_(time_in_ms_bw_samples),
ublox_check_flag(false){
    // Counter to keep track of how many valid llh position data was received from GPS during initialization
	gps_pos_count_ = 0;
    // Counter to keep track of how many valid ned velocity data was received from GPS during initialization
	gps_vel_count_ = 0;
    // Counter to keep track of how many valid 3D GPS fixes we got before we start capturing llh position and
    // ned velocity data to be used in EKF initialization
	gps_fix_count_ = 0;
    // Flag to indicate if GPS has a valid 3D fix or not
	gps_3d_fix_ = false;

	// Variables to store reference llh or the origin of NED frame
	lat_ref_ = 0;
	lon_ref_ = 0;
	height_ref_ = 0;

	// Variables to store initial NED velocity for EKF initialization
	vned_init_[0] = 0;
	vned_init_[1] = 0;
	vned_init_[2] = 0;

    // Variable to indicate GpsReadLoop() to stop
	stop_gps_read_loop_.store(false);

	sigint_flag = 0;

    /* Variable to indicate which GPS measurement has been updated
	first 3 indices are for position and last 3 indices are for velocities,
	this is shared between threads while running GpsReadLoop()
	*/
	// Set all the flags to false to show position and velocity data hasn't been received
	for( size_t idx_meas = 0; idx_meas < 6; idx_meas++ ){
		gps_meas_indices_[idx_meas] = false;
	}
}

bool GpsHelper::InitializeGps(float wait_duration_sec){
	bool gps_init_status = true;
	// Test GPS connection
	if(gps_.testConnection())
	{
		ublox_check_flag = true;
		printf("Ublox test OK\n");
		printf("\n");
		if ( !gps_.configureSolutionRate(time_in_ms_bw_samples_) )
		{
			printf("Setting new rate: FAILED\n");
		}
	}else{
		ublox_check_flag = false;
		printf("Ublox test failed\n");
		return false;
	}

	// initialize the duration
	// Recording the timestamp at the start of the code
    chrono::high_resolution_clock::time_point start_time;
    chrono::high_resolution_clock::time_point loop_end;
    chrono::microseconds delta (1000); 
	auto duration = chrono::duration_cast<chrono::microseconds> (delta);
    auto duration_count = duration.count();
    start_time = chrono::high_resolution_clock::now();

    // Run the loop till we have the desired number of llh position and NED velocity data for EKF initialization
	while( (gps_pos_count_ < n_gps_meas_count_) || (gps_vel_count_ < n_gps_meas_count_)){
		if  ( gps_3d_fix_ && ( gps_fix_count_ > n_valid_gps_count_ ) ){
			if(GetLlhPos() && gps_pos_count_ < n_gps_meas_count_){
				lat_ref_ += pos_data_[2];
				lon_ref_ += pos_data_[1];
				height_ref_ += pos_data_[3];
				gps_pos_count_++;
			}
		}

		if ( gps_3d_fix_ && ( gps_fix_count_ > n_valid_gps_count_ ) ){
			if(GetNedVel() && gps_vel_count_ < n_gps_meas_count_ ){
				vned_init_[0] += vel_data_[1];
				vned_init_[1] += vel_data_[2];
				vned_init_[2] += vel_data_[3];
				gps_vel_count_++;
			}
		}

		if( GetFixStatus() ){
			if ( (int)fix_data_[0] == GPS3DFIXID ){
				if (!gps_3d_fix_){
					cout<<"GPS 3D FIX OK\n"<<endl;
				}
				duration_count = 0;
				gps_3d_fix_ = true;
				// Increment the count that keeps track of how many valid 3D GPS fixes we have got so far
				if (gps_fix_count_ <= n_valid_gps_count_)
					gps_fix_count_++;
			}else{
				if (gps_3d_fix_){
					cout<<"GPS 3D FIX FAILED"<<endl;
				}
				gps_3d_fix_ = false;
				// Decrement the count that keeps track of how many valid 3D GPS fixes we have got so far if we loose the
				// 3D GPS fix during initialization
				if (gps_fix_count_ != 0)
					gps_fix_count_--;
			}
			
		}
		cout.flush();
		cout<<"GPS INITIALIZATION PROGRESS: "<< ( (float)(gps_pos_count_ + gps_vel_count_) * 100 )/(2*n_gps_meas_count_)<<" %\r";
		if(sigint_flag == 1 || (duration_count > (wait_duration_sec*1000000)) ){
			printf("GPS Initialization Failed\n");
			gps_init_status = false;
	   		break;
		}

	   	// Get the stop time and compute the duration
    	loop_end = std::chrono::high_resolution_clock::now();
    	duration_count = chrono::duration_cast<chrono::microseconds>(loop_end - start_time).count();
    }// While loop
    printf("\n");

    // Compute the initial llh position and initial NED velocity
    lat_ref_ = lat_ref_/n_gps_meas_count_;
    lon_ref_ = lon_ref_/n_gps_meas_count_;
    height_ref_ = height_ref_/n_gps_meas_count_;

    // Set the initial NED velocity
    vned_init_[0] = vned_init_[0]/n_gps_meas_count_;
    vned_init_[1] = vned_init_[1]/n_gps_meas_count_;
    vned_init_[2] = vned_init_[2]/n_gps_meas_count_;

    return gps_init_status;
}

void GpsHelper::GpsReadLoop(){
	if(ublox_check_flag){
		// Keep reading the GPS data forever
		while(1){
			// Set all the flags to false to show position and velocity data hasn't been received
			{
				unique_lock<mutex> gps_data_lock(gps_mutex_);
				for( size_t idx_meas = 0; idx_meas < 6; idx_meas++ ){
					gps_meas_indices_[idx_meas] = false;
				}
			}

	    	// // Check for llh position data if we have valid 3D GPS fix and number of 3D GPS fix is more than the required number of 3D GPS fixes
	    	if ( GetLlhPos() && gps_3d_fix_ && ( gps_fix_count_ > n_valid_gps_count_ ) )
			{
	        	// Convert the current llh position to NED position
				ned_pos_meas_ = Geodetic2Ned( pos_data_[2],
					pos_data_[1],
					pos_data_[3], 
					lat_ref_, 
					lon_ref_,
					height_ref_);
	        	// Set the flags to indicate llh position data is updated and assign data position data to the variables shared
	        	// between threads
				{
					unique_lock<mutex> gps_data_lock(gps_mutex_);
					for( size_t idx_meas = 0; idx_meas < 3; idx_meas++ ){
						gps_meas_indices_[idx_meas] = true;
						ned_pos_and_vel_meas_[idx_meas] = ned_pos_meas_(idx_meas);
					}
				}
			}

	    	// Check for NED velocity data if we have valid 3D GPS fix and number of 3D GPS fix is more than the required number of 3D GPS fixes
			if ( GetNedVel() && gps_3d_fix_ && ( gps_fix_count_ > n_valid_gps_count_ ) )
			{
	        	// Set the flags to indicate NED velocity data is updated and assign data position data to the variables shared
	        	// between threads
				unique_lock<mutex> gps_data_lock(gps_mutex_);
				for( size_t idx_meas = 3; idx_meas < 6; idx_meas++ ){
					gps_meas_indices_[idx_meas] = true;
					ned_pos_and_vel_meas_[idx_meas] = (float)(vel_data_[idx_meas - 2]);
				}
			}

	    	// Verify that we have valid 3D GPS fix
			if (GetFixStatus()){
				if ( (int)fix_data_[0] == GPS3DFIXID ){
					if (gps_fix_count_ <= n_valid_gps_count_)
						gps_fix_count_++;
					gps_3d_fix_ = true;
				}else{
					if (gps_fix_count_ != 0)
						gps_fix_count_--;
					gps_3d_fix_ = false;
				}
			}

			if(stop_gps_read_loop_.load()){
				break;
			}

			if(sigint_flag == 1)
		   		break;
	    }// While loop
	}
}// GpsReadLoop function

bool GpsHelper::GetFixStatus(){
	if ( gps_.decodeSingleMessage( Ublox::NAV_STATUS, fix_data_ ) == 1){
		return true;
	}else{
		return false;
	}
}

bool GpsHelper::GetNedVel(){
	if( gps_.decodeSingleMessage( Ublox::NAV_VELNED, vel_data_ ) == 1 ){
		vel_data_[1] = vel_data_[1]/100;
		vel_data_[2] = vel_data_[2]/100;
		vel_data_[3] = vel_data_[3]/100;
		return true;
	}else{
		return false;
	}
}

bool GpsHelper::GetLlhPos(){
	if( gps_.decodeSingleMessage( Ublox::NAV_POSLLH, pos_data_ ) == 1 ){
		pos_data_[1] = ( pos_data_[1] * DEG2RAD )/10000000;
		pos_data_[2] = ( pos_data_[2] * DEG2RAD )/10000000;
		pos_data_[3] = pos_data_[3]/1000;
		return true;
	}else{
		return false;
	}
}

void GpsHelper::CreateGpsThread(){
	gps_thread_ = thread(&GpsHelper::GpsReadLoop, this);

    CPU_ZERO(&cpuset_);
    CPU_SET(0, &cpuset_);

	pthread_getschedparam(gps_thread_.native_handle(), &policy_, &sch_);
	sch_.sched_priority = 5;
	pthread_setschedparam(gps_thread_.native_handle(), SCHED_FIFO, &sch_);
	int rc = pthread_setaffinity_np(gps_thread_.native_handle(),
                                    sizeof(cpuset_), &cpuset_);    
	if (rc != 0) {
      std::cerr << "Error calling pthread_setaffinity_np on GPS Thread: " << rc << "\n";
    }
}

GpsHelper::~GpsHelper(){
	// Making sure that GPS reading loop is stopped
	stop_gps_read_loop_.store(true);
	if(gps_thread_.joinable())
		gps_thread_.join();
}