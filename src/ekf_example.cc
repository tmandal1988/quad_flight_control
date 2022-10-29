#include <ekf_15dof_class.h>
#include <Matrix/matrix_factorization_class.h>
#include "Navio/Common/MPU9250.h"
#include "Navio/Navio2/LSM9DS1.h"
#include "Navio/Common/Util.h"
#include <Navio/Common/Ublox.h>
#include <memory>
#include <fstream>
#include <chrono>
#include <thread>
#include <mutex>

MatrixInv<float> ned_pos_and_vel_meas = {{0}, {0}, {0}, {0}, {0}, {0}};
float vned_init[] = {0, 0, 0};

bool gps_meas_indices[] = {false, false, false, 
						   false, false, false};
bool is_gps_initialized = false;

std::mutex gps_mutex;

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
	float ecc_sq = 0.0066943798522561;
    float r_ea = 6378137.0;
    float n_lat = r_ea/sqrt( ( 1 - ecc_sq * pow(sin(lat), 2) ) );

    float c_lat = cos(lat);
    float c_lon = cos(lon);
    float s_lat = sin(lat);
    float s_lon = sin(lon);

    MatrixInv<float> ecef_cord = { {(n_lat + height)*c_lat*c_lon}, {(n_lat + height)*c_lat*s_lon}, {(n_lat*(1 - ecc_sq) + height)*s_lat} };
    return ecef_cord;
}

MatrixInv<float> Geodetic2Ned(float lat, float lon, float height, float lat_ref, float lon_ref, float height_ref){
	MatrixInv<float> ecef_cord_ref = Geodetic2Ecef(lat_ref, lon_ref, height_ref);
	MatrixInv<float> ecef_cord = Geodetic2Ecef(lat, lon, height);

	MatrixInv<float> llh = { {lat}, {lon}, {height}};
	MatrixInv<float> llh_ref = {{lat_ref}, {lon_ref}, {height_ref}};

	float c_lat_ref = cos(lat_ref);
	float s_lat_ref = sin(lat_ref);

	float c_lon_ref = cos(lon_ref);
	float s_lon_ref = sin(lon_ref);

    MatrixInv<float> ecef2ned = { {-s_lat_ref*c_lon_ref, -s_lat_ref*s_lon_ref, c_lat_ref}, 
    							  {-s_lon_ref, c_lon_ref, 0}, 
    							  {-c_lat_ref*c_lon_ref, -c_lat_ref*s_lon_ref, -s_lat_ref} };
   	return ecef2ned*(llh - llh_ref);

}

void GetGpsData(){
	Ublox gps;
    if(gps.testConnection())
    {
        printf("Ublox test OK\n");
        if (!gps.configureSolutionRate(100))
        {
            printf("Setting new rate: FAILED\n");
        }
    }

    std::vector<double> fix_data;
    size_t gps_pos_count = 0;
    size_t gps_vel_count = 0;
    size_t gps_fix_count = 0;
    size_t n_gps_meas_count = 10;
    size_t n_valid_gps_count = 10;
    bool gps_3d_fix = false;
    double lat_ref = 0;
    double lon_ref = 0;
    double height_ref = 0;

    double vn_init = 0;
    double ve_init = 0;
    double vd_init = 0;

    MatrixInv<float> ned_pos_meas;
    vector<double> pos_data;
    vector<double> vel_data;
    size_t count = 0;
    while( (gps_pos_count < n_gps_meas_count) || (gps_vel_count < n_gps_meas_count) ){
    	if  ( gps_3d_fix && ( gps_fix_count > n_valid_gps_count ) ) 
        {
        	if((gps.decodeSingleMessage(Ublox::NAV_POSLLH, pos_data) == 1) && (gps_pos_count < n_gps_meas_count) ){
        		lat_ref += pos_data[1]/10000000;
        		lon_ref += pos_data[2]/10000000;
        		height_ref += pos_data[3]/1000;
        		gps_pos_count++;
        	}
        }  

        if ( gps_3d_fix && ( gps_fix_count > n_valid_gps_count) )
        {
        	if( (gps.decodeSingleMessage(Ublox::NAV_VELNED, vel_data) == 1) && (gps_vel_count < n_gps_meas_count) ){
        		vn_init += vel_data[1]/100;
        		ve_init += vel_data[2]/100;
        		vd_init += vel_data[3]/100;
        		gps_vel_count++;
        	}
        }  
        if (gps.decodeSingleMessage(Ublox::NAV_STATUS, fix_data) == 1){
        	switch((int)fix_data[0]){
        		case 0x03:
        			if (!gps_3d_fix)
        				gps_mutex.lock();
        				cout<<"GPS 3D FIX OK"<<endl;
        				gps_mutex.unlock();
        			gps_3d_fix = true;
        			if (gps_fix_count <= n_valid_gps_count)
        				gps_fix_count++;
        			break;
        		default:
        			if (gps_3d_fix)
        				gps_mutex.lock();
        				cout<<"GPS 3D FIX FAILED"<<endl;
        				gps_mutex.unlock();
        			if (gps_fix_count != 0)
        				gps_fix_count--;
        			gps_3d_fix = false;
        	}
        }
        gps_mutex.lock();
        printf("GPS INITIALIZATION PROGRESS: %g%%\r", ( (float)(gps_pos_count + gps_vel_count) * 100 )/(2*n_gps_meas_count) );
        gps_mutex.unlock();
    }

    gps_mutex.lock();
    gps_mutex.unlock();
    lat_ref = (lat_ref * DEG2RAD)/n_gps_meas_count;
    lon_ref = (lon_ref * DEG2RAD)/n_gps_meas_count;
    height_ref = height_ref/n_gps_meas_count;

    gps_mutex.lock();
    vned_init[0] = vn_init/n_gps_meas_count;
    vned_init[1] = ve_init/n_gps_meas_count;
    vned_init[2] = vd_init/n_gps_meas_count;    
    is_gps_initialized = true;
    gps_mutex.unlock();

    while(1){
    	gps_mutex.lock();
    	for( size_t idx_meas = 0; idx_meas < 6; idx_meas++ ){
    		gps_meas_indices[idx_meas] = false;
    	}
    	gps_mutex.unlock();

    	if ( gps.decodeSingleMessage(Ublox::NAV_POSLLH, pos_data) == 1 && gps_3d_fix && ( gps_fix_count > n_valid_gps_count ))
        	{
        		ned_pos_meas = Geodetic2Ned( (pos_data[2] * DEG2RAD)/10000000,
        											  (pos_data[1] * DEG2RAD)/10000000,
        											   pos_data[3]/1000, lat_ref, lon_ref,
        											   height_ref);

        		gps_mutex.lock();
        		for( size_t idx_meas = 0; idx_meas < 3; idx_meas++ ){
    				gps_meas_indices[idx_meas] = true;
    				ned_pos_and_vel_meas(idx_meas) = ned_pos_meas(idx_meas);
    			}
    			gps_mutex.unlock();
        	}  

        	if (gps.decodeSingleMessage(Ublox::NAV_VELNED, vel_data) == 1 && gps_3d_fix && ( gps_fix_count > n_valid_gps_count ))
        	{
        		gps_mutex.lock();
        		for( size_t idx_meas = 3; idx_meas < 6; idx_meas++ ){
    				gps_meas_indices[idx_meas] = true;
    				ned_pos_and_vel_meas(idx_meas) = (float)(vel_data[idx_meas - 2]/100);
    			}
    			gps_mutex.unlock();
        	}  

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

	auto sensor = get_inertial_sensor("mpu");

	if (!sensor) {
        printf("Wrong sensor name. Select: mpu or lsm\n");
        return EXIT_FAILURE;
    }

    if (!sensor->probe()) {
        printf("Sensor not enabled\n");
        return EXIT_FAILURE;
    }
    sensor->initialize();

    thread gps_thread(GetGpsData);

    float ax, ay, az;
    float gx, gy, gz;
    float mx, my, mz;

    MatrixInv<float> initial_state(15, 1);
	MatrixInv<float> sensor_meas(9, 1);
	MatrixInv<float> state_sensor_val(6, 1);
	MatrixInv<float> process_noise_q(15, 15, "eye");
	MatrixInv<float> meas_noise_r(7, 7, "eye");

	
	MatrixInv<float> temp_mag(3, 1);

	MatrixInv<float> current_state(15, 1);
	MatrixInv<float> state_jacobian(1, 15);
	MatrixInv<float> computed_meas(7, 1);

	bool meas_indices[] = {false,
						   false, false, false,
						   false, false, false};

	//Compute initial heading
	float magnetic_declination = 13.01*DEG2RAD;

	MatrixInv<float> mag_offset = {{17.8902136639429}, {35.5186740453011}, {-33.8067089624238}};
	MatrixInv<float> mag_a(3, 3, "eye");
	MatrixInv<float> mag_scale = {{0.9652, 0, 0}, {0, 1.09, 0}, {0, 0, 0.9556}};


	size_t count = 0;
	int t_idx = atoi(argv[1]);
 	cout<<t_idx<<endl;

	// Compute initial roll, pitch, yaw but collecting and averaging data for 2 sec
	float sum_ax = 0;
	float sum_ay = 0;
	float sum_az = 0;
	float sum_mx = 0;
	float sum_my = 0;
	float sum_mz = 0;
	for (size_t idx_s = 0; idx_s < 200; idx_s++){
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

	float avg_acc_mag = sqrt(pow((sum_ax/200), 2) + pow((sum_ay/200), 2) + pow((sum_az/200), 2));
	initial_state(0) = atan2(sum_ay/200, avg_acc_mag);
	initial_state(1) = atan2(sum_ax/200, avg_acc_mag);

	MatrixInv<float> mag2d_projection(2, 3);
	mag2d_projection(0, 0) = cos(initial_state(1));
	mag2d_projection(0, 1) = sin(initial_state(1))*sin(initial_state(0));
	mag2d_projection(0, 2) =  sin(initial_state(1))*cos(initial_state(0));

	mag2d_projection(1, 1)  = cos(initial_state(0));
	mag2d_projection(1, 2) = -sin(initial_state(0));

	MatrixInv<float> mag_vector(3, 1);
	mag_vector(0) = sum_mx/200;
	mag_vector(1) = sum_my/200;
	mag_vector(2) = sum_mz/200;

	mag_vector = mag_a*mag_scale*(mag_vector - mag_offset);
	float avg_gaus_mag = sqrt(pow(mag_vector(0), 2) + pow(mag_vector(1), 2) + pow(mag_vector(2), 2));
	mag_vector = mag_vector/avg_gaus_mag;

	MatrixInv<float> mag2d = mag2d_projection*mag_vector;
	initial_state(2) = atan2(-mag2d(1), mag2d(0)) + magnetic_declination;

	// Wait for gps to initialize
	cout<<"Waiting for 3D GPS Fix: "<<endl;
	while(!is_gps_initialized){cout<<is_gps_initialized<<endl; usleep(1000000);}

	initial_state(9) = vned_init[0];
	initial_state(10) = vned_init[1];
	initial_state(11) = vned_init[2];

	printf("Initial Roll: %g [deg], Initial Pitch: %g [deg], Initial Yaw: %g [deg]\n", initial_state(0)*RAD2DEG, initial_state(1)*RAD2DEG, initial_state(2)*RAD2DEG);

	usleep(50);

	Ekf15Dof<float> imu_gps_ekf(0.01, initial_state, process_noise_q*0.00001, meas_noise_r*0.01);

	// // file pointer
    fstream fout;
  
 //    // opens an existing csv file or creates a new file.
    fout.open("imu.csv", ios::out | ios::app);

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

	double avg_duration = 0;	
	meas_indices[0] = true;
    while(1) {
    	auto start = std::chrono::high_resolution_clock::now();
		for( size_t idx_meas = 1; idx_meas < 7; idx_meas++ ){
    		meas_indices[idx_meas] = false;
    	}

	    sensor->update();
	    sensor->read_accelerometer(&ax, &ay, &az);
	    sensor->read_gyroscope(&gx, &gy, &gz);
	    sensor->read_magnetometer(&mx, &my, &mz);

	    state_sensor_val(0) = gy;
	    state_sensor_val(1) = gx;
	    state_sensor_val(2) = -gz;
	    state_sensor_val(3) = ay;
	    state_sensor_val(4) = ax;
	    state_sensor_val(5) = -az;

	    temp_mag(0) = mx;
	    temp_mag(1) = my;
	    temp_mag(2) = mz;
	    temp_mag = mag_scale*mag_a*(temp_mag - mag_offset);
	    sensor_meas(0) = temp_mag(0);
	    sensor_meas(1) = temp_mag(1);
	    sensor_meas(2) = temp_mag(2);

	    for(size_t idx_meas = 0; idx_meas < 6; idx_meas++){
	    	if(gps_meas_indices[idx_meas]){
	    		sensor_meas(idx_meas + 3) = ned_pos_and_vel_meas(idx_meas);
	    		meas_indices[idx_meas + 1] = true;
	    	}
	    }

	    imu_gps_ekf.Run(state_sensor_val, sensor_meas, meas_indices);
        current_state = imu_gps_ekf.GetCurrentState();
	    if (gps_meas_indices[1]){
	    	printf("Roll [deg]: %+7.3f, Pitch[deg]: %+7.3f, Yaw[deg]: %+7.3f\n", current_state(0)*RAD2DEG, current_state(1)*RAD2DEG, current_state(2)*RAD2DEG);
	    	printf("Pos N [m]: %+7.3f, Pos E [m]: %+7.3f, Pos D[m]: %+7.3f\n", current_state(6), current_state(7), current_state(8));
	    	printf("Vel N [m]: %+7.3f, Vel E [m]: %+7.3f, Vel D[m]: %+7.3f\n", current_state(9), current_state(10), current_state(11));
	    	printf("############################################\n");
		}

	    usleep(7100);

	    count++;
	    auto stop = std::chrono::high_resolution_clock::now();
	    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
	    avg_duration = avg_duration + duration.count();
	    //cout<<"duration.count(): "<<duration.count()<<endl;
	}
	cout<<"Avg Execution Time: "<<avg_duration/100<<endl;
	fout.flush();
	gps_thread.join();


	return 0;
}