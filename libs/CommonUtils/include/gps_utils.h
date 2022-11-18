#ifndef GPSUTILS_H
#define GPSUTILS_H

#include<Navio/Common/Ublox.h>
#include<Matrix/matrix_inv_class.h>
#include<constants.h>
#include<coordinate_transformation.h>

#include<vector>
#include<atomic>
#include<mutex>
#include<thread>
#include<pthread.h>

//https://content.u-blox.com/sites/default/files/products/documents/u-blox8-M8_ReceiverDescrProtSpec_UBX-13003221.pdf
using namespace std;

class GpsHelper{
	#define GPS3DFIXID 0x03 

	public:
		// Constructors
		GpsHelper(size_t n_valid_gps_count = 10, size_t n_gps_meas_count = 10, uint16_t time_in_ms_bw_samples = 100);

		// Destructor
		~GpsHelper();

		// To initialize the GPS
		void InitializeGps();

		// Useful function that can be used to read GPS data in a loop at the configured rate
		// This function can be passed to a thread to update position and velocity data in thread safe way
		void GpsReadLoop();

		// Start the GPS running thread
		void CreateGpsThread();

		// Get init NED Vel
		double* GetInitNedVel();

		// Get NED position and velocity
		void GetGpsNedPosAndVel(float (&ned_pos_and_vel_meas)[6], bool (&gps_meas_indices)[6]);

		// Stops the above loop
		void StopGpsReadLoop();

	private:
		// Create a thread
		thread gps_thread_;
		// To set CPU affinity
		cpu_set_t cpuset_;
		sched_param sch_;
		int policy_;
		// Create an Ublox object to read GPS data
		Ublox gps_;
		// Configuration parameter to set GPS sample rate
		uint16_t time_in_ms_bw_samples_;
		// Variable to store GPS fix quality data
	    vector<double> fix_data_;
	    // Counter to keep track of how many valid llh position data was received from GPS during initialization
	    size_t gps_pos_count_;
	    // Counter to keep track of how many valid ned velocity data was received from GPS during initialization
	    size_t gps_vel_count_;
	    // Counter to keep track of how many valid 3D GPS fixes we got before we start capturing llh position and
	    // ned velocity data to be used in EKF initialization
	    size_t gps_fix_count_;
	    // Number of valid llh position data to use for EKF initialization
	    size_t n_gps_meas_count_;
	    // Number of valid 3d GPS fixes required before capturing llh position and ned velocity data for EKF initialization
	    size_t n_valid_gps_count_;
	    // Flag to indicate if GPS has a valid 3D fix or not
	    bool gps_3d_fix_;

	    // Variables to store reference llh or the origin of NED frame
	    double lat_ref_;
	    double lon_ref_;
	    double height_ref_;

	    // Variables to store initial NED velocity for EKF initialization
	    double vned_init_[3];

	    // Measured NED position
	    MatrixInv<float> ned_pos_meas_{MatrixInv<float>(3, 1)};
	    // llh position data at each iteration received from the GPS
	    vector<double> pos_data_;
	    // NED velocity data at each iteration received from the GPS
	    vector<double> vel_data_;
	    // variable to store ned pos and ned velocity
	    float ned_pos_and_vel_meas_[6];

	    // Mutex to guard resource access between threads while running GpsReadLoop()
	    mutex gps_mutex_;

	    // Variable to indicate GpsReadLoop() to stop
	    atomic<bool> stop_gps_read_loop_;

	    /* Variable to indicate which GPS measurement has been updated
		first 3 indices are for position and last 3 indices are for velocities,
		this is shared between threads while running GpsReadLoop()
		*/
	    bool gps_meas_indices_[6];

	    // Function to get fix status
	    bool GetFixStatus();
	    // Function to get NED Velocity
	    bool GetNedVel();
	    // Function to get llh position
	    bool GetLlhPos();
};
#endif