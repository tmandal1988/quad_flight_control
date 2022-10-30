// EKF related headers
#include <ekf_15dof_class.h>
// Matrix library
#include <Matrix/matrix_factorization_class.h>
// Navio 2 Utilities
#include "Navio/Common/MPU9250.h"
#include "Navio/Navio2/LSM9DS1.h"
#include "Navio/Common/Util.h"
#include <Navio/Common/Ublox.h>
// Standard C++ Libraries file, time and memory
#include <memory>
#include <fstream>
#include <chrono>
// Standard C++ Libraries for multi-threading
#include <thread>
#include <mutex>

// Variable to store gps measured position and velocity, this is shared between threads
MatrixInv<float> ned_pos_and_vel_meas = {{0}, {0}, {0}, {0}, {0}, {0}};
// Variable to store initial velocity as measured by GPS during EKF initialization
float vned_init[] = {0, 0, 0};

// Variable to indicate if GPS has been initialized, this is shared between threads
bool is_gps_initialized = false;

/* Variable to indicate which GPS measurement has been updated
first 3 indices are for position and last 3 indices are for velocities,
this is shared between threads
*/
bool gps_meas_indices[] = {false, false, false, 
						   false, false, false};


// Mutex to guard resource access between threads
std::mutex gps_mutex;

// Function to construct an object of inertial sensor class for a given IMU type
std::unique_ptr <InertialSensor> get_inertial_sensor( std::string sensor_name)
{
    if (sensor_name == "mpu") {
        printf("Selected: MPU9250\n");
        auto ptr = std::unique_ptr <InertialSensor>{ new MPU9250() };
        return ptr;
    }
    else if (sensor_name == "lsm") {
        printf("Selected: LSM9DS1\n");
        auto ptr = std::unique_ptr <InertialSensor>{ new LSM9DS1() };
        return ptr;
    }
    else {
        return NULL;
    }
}

MatrixInv<float> Geodetic2Ecef(float lat, float lon, float height){
	/* This function converts Geodetic positions to ECEF positions
	Inputs:
		lat (float)		: Lattitude measured by GPS
		lon (float)		: Longitude measured by GPS
		height (float)	: Height above ellipsoid measured by GPS
	Outputs:
		MatrixInv<float>: 3x1 array of ECEF X, Y, Z coordinates
	*/

	// Square of Earth's eccentricity
	float ecc_sq = 0.0066943798522561;
	// Earth's semi measure axis
    float r_ea = 6378137.0;
    // Intermediate variable for coordinate conversion
    float n_lat = r_ea/sqrt( ( 1 - ecc_sq * pow(sin(lat), 2) ) );

    float c_lat = cos(lat);
    float c_lon = cos(lon);
    float s_lat = sin(lat);
    float s_lon = sin(lon);

    // Transform the coordinates and return
    MatrixInv<float> ecef_cord = { {(n_lat + height)*c_lat*c_lon}, {(n_lat + height)*c_lat*s_lon}, {(n_lat*(1 - ecc_sq) + height)*s_lat} };
    return ecef_cord;
}

MatrixInv<float> Geodetic2Ned(float lat, float lon, float height, float lat_ref, float lon_ref, float height_ref){
	/* This function converts Geodetic positions to NED positions
	Inputs:
		lat (float)			: Lattitude measured by GPS
		lon (float)			: Longitude measured by GPS
		height (float)		: Height above ellipsoid measured by GPS
		lat_ref (float) 	: Lattitude of the origin of the NED frame
		lon_ref (float) 	: Longitude of the origin of the NED frame
		height_ref (float)	: Height of the origin of the NED frame
	Outputs:
		MatrixInv<float> 	: 3x1 array of N, E, D coordinates
	*/

	// Convert the origin and current llh to ECEF coordinates
	MatrixInv<float> ecef_cord_ref = Geodetic2Ecef(lat_ref, lon_ref, height_ref);
	MatrixInv<float> ecef_cord = Geodetic2Ecef(lat, lon, height);

	float c_lat_ref = cos(lat_ref);
	float s_lat_ref = sin(lat_ref);

	float c_lon_ref = cos(lon_ref);
	float s_lon_ref = sin(lon_ref);

	// Create ECEF to NED transformation
    MatrixInv<float> ecef2ned = { {-s_lat_ref*c_lon_ref, -s_lat_ref*s_lon_ref, c_lat_ref}, 
    							  {-s_lon_ref, c_lon_ref, 0}, 
    							  {-c_lat_ref*c_lon_ref, -c_lat_ref*s_lon_ref, -s_lat_ref} };
   // Return NED coordinates
   	return ecef2ned*(ecef_cord - ecef_cord_ref);

}

void GetGpsData(){
	// Create an Ublox object to read GPS data
	Ublox gps;
	// Test GPS connection
    if(gps.testConnection())
    {
        printf("Ublox test OK\n");
        if (!gps.configureSolutionRate(100))
        {
            printf("Setting new rate: FAILED\n");
        }
    }

    // Variable to store GPS fix quality data
    std::vector<double> fix_data;
    // Counter to keep track of how many valid llh position data was received from GPS during initialization
    size_t gps_pos_count = 0;
    // Counter to keep track of how many valid ned velocity data was received from GPS during initialization
    size_t gps_vel_count = 0;
    // Counter to keep track of how many valid 3D GPS fixes we got before we start capturing llh position and
    // ned velocity data to be used in EKF initialization
    size_t gps_fix_count = 0;
    // Number of valid llh position data to use for EKF initialization
    size_t n_gps_meas_count = 10;
    // Number of valid 3d GPS fixes required before capturing llh position and ned velocity data for EKF initialization
    size_t n_valid_gps_count = 10;
    // Flag to indicate if GPS has a valid 3D fix or not
    bool gps_3d_fix = false;

    // Variables to store reference llh or the origin of NED frame
    double lat_ref = 0;
    double lon_ref = 0;
    double height_ref = 0;

    // Variables to store initial NED velocity for EKF initialization
    double vn_init = 0;
    double ve_init = 0;
    double vd_init = 0;

    // Measured NED position
    MatrixInv<float> ned_pos_meas;
    // llh position data at each iteration received from the GPS
    vector<double> pos_data;
    // NED velocity data at each iteration received from the GPS
    vector<double> vel_data;

    // Run the loop till we have the desired number of llh position and NED velocity data for EKF initialization
    while( (gps_pos_count < n_gps_meas_count) || (gps_vel_count < n_gps_meas_count) ){
    	
    	if  ( gps_3d_fix && ( gps_fix_count > n_valid_gps_count ) ) 
        {
        	// Check for llh position if we have valid GPS 3D fix required number of times and if the number of llh position
        	// received so far is less than the required number of llh position needed
        	if((gps.decodeSingleMessage(Ublox::NAV_POSLLH, pos_data) == 1) && (gps_pos_count < n_gps_meas_count) ){
        		// Get the lat, lon and height and sum it up
        		lat_ref += pos_data[2]/10000000;
        		lon_ref += pos_data[1]/10000000;
        		height_ref += pos_data[3]/1000;
        		gps_pos_count++;
        	}
        }  

        if ( gps_3d_fix && ( gps_fix_count > n_valid_gps_count) )
        {
        	// Check for NED velocity if we have valid GPS 3D fix required number of times and if the number of NED velocity
        	// received so far is less than the required number of NED velocity needed
        	if( (gps.decodeSingleMessage(Ublox::NAV_VELNED, vel_data) == 1) && (gps_vel_count < n_gps_meas_count) ){
        		// Get the NED velocity and sum it up
        		vn_init += vel_data[1]/100;
        		ve_init += vel_data[2]/100;
        		vd_init += vel_data[3]/100;
        		gps_vel_count++;
        	}
        }  
        if (gps.decodeSingleMessage(Ublox::NAV_STATUS, fix_data) == 1){
        	// Verify that we have 3D GPS fix
        	switch((int)fix_data[0]){
        		case 0x03:
        			if (!gps_3d_fix){
        				gps_mutex.lock();
        				cout<<"GPS 3D FIX OK"<<endl;
        				gps_mutex.unlock();
        			}
        			gps_3d_fix = true;
        			// Increment the count that keeps track of how many valid 3D GPS fixes we have got so far
        			if (gps_fix_count <= n_valid_gps_count)
        				gps_fix_count++;
        			break;
        		default:
        			if (gps_3d_fix){
        				gps_mutex.lock();
        				cout<<"GPS 3D FIX FAILED"<<endl;
        				gps_mutex.unlock();
        			}
        			// Decrement the count that keeps track of how many valid 3D GPS fixes we have got so far if we loose the
        			// 3D GPS fix during initialization
        			if (gps_fix_count != 0)
        				gps_fix_count--;
        			gps_3d_fix = false;
        	}
        }
        gps_mutex.lock();
        printf("GPS INITIALIZATION PROGRESS: %g%%\r", ( (float)(gps_pos_count + gps_vel_count) * 100 )/(2*n_gps_meas_count) );
        gps_mutex.unlock();
    }

  	// Compute the initial llh position and initial NED velocity
    lat_ref = (lat_ref * DEG2RAD)/n_gps_meas_count;
    lon_ref = (lon_ref * DEG2RAD)/n_gps_meas_count;
    height_ref = height_ref/n_gps_meas_count;

    gps_mutex.lock();
    vned_init[0] = vn_init/n_gps_meas_count;
    vned_init[1] = ve_init/n_gps_meas_count;
    vned_init[2] = vd_init/n_gps_meas_count;    
    // Set the flag that indicates that GPS has been initialized
    is_gps_initialized = true;
    gps_mutex.unlock();

    // Keep reading the GPS data forever
    while(1){
    	// Set all the flags to false to show position and velocity data hasn't been received
    	gps_mutex.lock();
    	for( size_t idx_meas = 0; idx_meas < 6; idx_meas++ ){
    		gps_meas_indices[idx_meas] = false;
    	}
    	gps_mutex.unlock();

    	// Check for llh position data if we have valid 3D GPS fix and number of 3D GPS fix is more than the required number of 3D GPS fixes
    	if ( gps.decodeSingleMessage(Ublox::NAV_POSLLH, pos_data) == 1 && gps_3d_fix && ( gps_fix_count > n_valid_gps_count ))
        	{
        		// Convert the current llh position to NED position
        		ned_pos_meas = Geodetic2Ned( (pos_data[2] * DEG2RAD)/10000000,
        											  (pos_data[1] * DEG2RAD)/10000000,
        											   pos_data[3]/1000, lat_ref, lon_ref,
        											   height_ref);

        		// Set the flags to indicate llh position data is updated and assign data position data to the variables shared
        		// between threads
        		gps_mutex.lock();
        		for( size_t idx_meas = 0; idx_meas < 3; idx_meas++ ){
    				gps_meas_indices[idx_meas] = true;
    				ned_pos_and_vel_meas(idx_meas) = ned_pos_meas(idx_meas);
    			}
    			gps_mutex.unlock();
        	}  

        	// Check for NED velocity data if we have valid 3D GPS fix and number of 3D GPS fix is more than the required number of 3D GPS fixes
        	if (gps.decodeSingleMessage(Ublox::NAV_VELNED, vel_data) == 1 && gps_3d_fix && ( gps_fix_count > n_valid_gps_count ))
        	{
        		// Set the flags to indicate NED velocity data is updated and assign data position data to the variables shared
        		// between threads
        		gps_mutex.lock();
        		for( size_t idx_meas = 3; idx_meas < 6; idx_meas++ ){
    				gps_meas_indices[idx_meas] = true;
    				ned_pos_and_vel_meas(idx_meas) = (float)(vel_data[idx_meas - 2]/100);
    			}
    			gps_mutex.unlock();
        	}  

        	// Verify that we have valid 3D GPS fix
        	if (gps.decodeSingleMessage(Ublox::NAV_STATUS, fix_data) == 1){
        		switch((int)fix_data[0]){
        			case 0x03:
        				if (gps_fix_count <= n_valid_gps_count)
        					gps_fix_count++;
        				gps_3d_fix = true;
        				break;
        			default:
        				if (gps_fix_count != 0)
        					gps_fix_count--;
        				gps_3d_fix = false;
        		}
        	}
    }

}

int main(int argc, char *argv[]){

	// Create IMU sensor object
	auto sensor = get_inertial_sensor("mpu");

	// Verify that correct IMU has been selected
	if (!sensor) {
        printf("Wrong sensor name. Select: mpu or lsm\n");
        return EXIT_FAILURE;
    }

    // Verify that the IMU is available
    if (!sensor->probe()) {
        printf("Sensor not enabled\n");
        return EXIT_FAILURE;
    }

    // Initialize the IMU
    sensor->initialize();

    // Create a seperate thread to received and update GPS data for the EKF
    thread gps_thread(GetGpsData);

    // Variables to read data from the IMU
    // Accels
    float ax, ay, az;
    //Gyros
    float gx, gy, gz;
    //Mags
    float mx, my, mz;

    // Initial State Variable
    MatrixInv<float> initial_state(15, 1);
    // Variable used in the measurement update of the EKF
	MatrixInv<float> sensor_meas(9, 1);
	// Sensor values used in the time propagation stage of the EKF
	MatrixInv<float> state_sensor_val(6, 1);
	// Q matrix of the EKF
	MatrixInv<float> process_noise_q(15, 15, "eye");
	// R matrix of the EKF
	MatrixInv<float> meas_noise_r(7, 7, "eye");

	// Variable to store calibrated mag data
	MatrixInv<float> temp_mag(3, 1);

	// Variable to store current state at the end of each EKF run
	MatrixInv<float> current_state(15, 1);

	// Array to indicate which measurement has been update
	bool meas_indices[] = {false,
						   false, false, false,
						   false, false, false};

	// Magnetic field declination in the Bay Area
	float magnetic_declination = 13.01*DEG2RAD;

	// Calibration variables for mag
	// Offset
	MatrixInv<float> mag_offset = {{17.8902136639429}, {35.5186740453011}, {-33.8067089624238}};
	// Rotation matrix for misalignment of the IMU axes
	MatrixInv<float> mag_a(3, 3, "eye");
	// Scale factor for each axis of the IMU
	MatrixInv<float> mag_scale = {{0.9652, 0, 0}, {0, 1.09, 0}, {0, 0, 0.9556}};

	// Compute initial roll, pitch, yaw but collecting and averaging data for 2 sec
	float sum_ax = 0;
	float sum_ay = 0;
	float sum_az = 0;
	float sum_mx = 0;
	float sum_my = 0;
	float sum_mz = 0;
	// Number of measurements to use for initialization
	size_t imu_count = 200;
	for (size_t idx_s = 0; idx_s < imu_count; idx_s++){
		// Read IMU data and sum them up
		sensor->update();
	    sensor->read_accelerometer(&ax, &ay, &az);
	    sensor->read_magnetometer(&mx, &my, &mz);

	    sum_ax += ay;
	    sum_ay += ax;
	    sum_az += -az;

	    sum_mx += mx;
	    sum_my += my;
	    sum_mz += mz;

		usleep(10);
	}

	// Compute average accel magnitude
	float avg_acc_mag = sqrt(pow((sum_ax/imu_count), 2) + pow((sum_ay/imu_count), 2) + pow((sum_az/imu_count), 2));
	// Set initial roll and pitch angle
	initial_state(0) = atan2(sum_ay/imu_count, avg_acc_mag);
	initial_state(1) = atan2(sum_ax/imu_count, avg_acc_mag);

	// Matrix to project 3D mag to 2D
	MatrixInv<float> mag2d_projection(2, 3);
	mag2d_projection(0, 0) = cos(initial_state(1));
	mag2d_projection(0, 1) = sin(initial_state(1))*sin(initial_state(0));
	mag2d_projection(0, 2) =  sin(initial_state(1))*cos(initial_state(0));

	mag2d_projection(1, 1)  = cos(initial_state(0));
	mag2d_projection(1, 2) = -sin(initial_state(0));

	// Compute average mag vector
	MatrixInv<float> mag_vector(3, 1);
	mag_vector(0) = sum_mx/imu_count;
	mag_vector(1) = sum_my/imu_count;
	mag_vector(2) = sum_mz/imu_count;
	// Correct for mag hard and soft iron errors
	mag_vector = mag_a*mag_scale*(mag_vector - mag_offset);
	// Normalize avg mag vector
	float avg_gaus_mag = sqrt(pow(mag_vector(0), 2) + pow(mag_vector(1), 2) + pow(mag_vector(2), 2));
	mag_vector = mag_vector/avg_gaus_mag;

	// Compute heading from mag readings
	MatrixInv<float> mag2d = mag2d_projection*mag_vector;
	initial_state(2) = atan2(-mag2d(1), mag2d(0)) + magnetic_declination;

	// Wait for gps to initialize
	cout<<"Waiting for 3D GPS Fix: "<<endl;
	while(!is_gps_initialized){usleep(1000000);}

	// Initial NED velocity
	initial_state(9) = vned_init[0];
	initial_state(10) = vned_init[1];
	initial_state(11) = vned_init[2];


	// Create an 15 state EKF object
	Ekf15Dof<float> imu_gps_ekf(0.01, initial_state, process_noise_q*0.00001, meas_noise_r*0.01);

	// file pointer to save data to if needed
    fstream fout;
  
 	// opens an existing csv file or creates a new file.
    fout.open("imu.csv", ios::out | ios::app);

    // Save the initial state in the first row
    fout << 0 << ", "
             << 0 << ", "
             << 0 << ", "
             << 0 << ", "
             << 0 << ", "
             << 0 << ","
             << 0 << ", "
             << 0 << ", "
             << 0 << ", "
             << initial_state(0) << ", "
             << initial_state(1) << ", "
             << initial_state(2) << ", "
             << 0 << ", "<< 0 << ", " << 0 << ", " << 0 << ", " << 0 << ", " << 0 << ", " << 0 << ", " << 0 << ", " << 0 << ", " << 0 << ", " << 0 << ", " << 0 << ", "
             << 0 << ", "<< 0 << ", " << 0
             << "\n";

	// Make sure EKF always uses Magnetometer in the update step
	meas_indices[0] = true;

	// Loop counter
	size_t loop_count = 0;
    while(1) {
    	// Get the current time
    	auto start = std::chrono::high_resolution_clock::now();
    	// Make all the measurement flags corresponding to position and velocity false
		for( size_t idx_meas = 1; idx_meas < 7; idx_meas++ ){
    		meas_indices[idx_meas] = false;
    	}

    	// Read IMU data
	    sensor->update();
	    sensor->read_accelerometer(&ax, &ay, &az);
	    sensor->read_gyroscope(&gx, &gy, &gz);
	    sensor->read_magnetometer(&mx, &my, &mz);

	    // Assign required sensor value for time propagation of the ekf state
	    state_sensor_val(0) = gy;
	    state_sensor_val(1) = gx;
	    state_sensor_val(2) = -gz;
	    state_sensor_val(3) = ay;
	    state_sensor_val(4) = ax;
	    state_sensor_val(5) = -az;

	    // Correct and assign mad data to the sensor measurement array
	    temp_mag(0) = mx;
	    temp_mag(1) = my;
	    temp_mag(2) = mz;
	    temp_mag = mag_scale*mag_a*(temp_mag - mag_offset);
	    sensor_meas(0) = temp_mag(0);
	    sensor_meas(1) = temp_mag(1);
	    sensor_meas(2) = temp_mag(2);

	    // Check if GPS position and velocity has been update and set appropriate flags
	    for(size_t idx_meas = 0; idx_meas < 6; idx_meas++){
	    	if(gps_meas_indices[idx_meas]){
	    		sensor_meas(idx_meas + 3) = ned_pos_and_vel_meas(idx_meas);
	    		meas_indices[idx_meas + 1] = true;
	    	}
	    }

	    // Run one step of EKF
	    imu_gps_ekf.Run(state_sensor_val, sensor_meas, meas_indices);
	    // Get the state after EKF run
        current_state = imu_gps_ekf.GetCurrentState();
	    if (remainder(loop_count, 50) == 0){
	    	printf("Roll [deg]: %+7.3f, Pitch[deg]: %+7.3f, Yaw[deg]: %+7.3f\n", current_state(0)*RAD2DEG, current_state(1)*RAD2DEG, current_state(2)*RAD2DEG);
	    	printf("Pos N [m]: %+7.3f, Pos E [m]: %+7.3f, Pos D[m]: %+7.3f\n", current_state(6), current_state(7), current_state(8));
	    	printf("Vel N [m]: %+7.3f, Vel E [m]: %+7.3f, Vel D[m]: %+7.3f\n", current_state(9), current_state(10), current_state(11));
	    	printf("############################################\n");
		}

		loop_count++;
		// This sleep number was tuned to get ~100Hz of loop execution time
	    usleep(7100);

	    // Get the stop time and compute the duration
	    auto stop = std::chrono::high_resolution_clock::now();

	    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
	}

	// Wait for GPS thread to finish
	gps_thread.join();


	return 0;
}