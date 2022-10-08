#include <ekf_15dof_class.h>
#include <Matrix/matrix_factorization_class.h>
#include "Navio/Common/MPU9250.h"
#include "Navio/Navio2/LSM9DS1.h"
#include "Navio/Common/Util.h"
#include <memory>
#include <fstream>
/*#include <chrono>
*/
#define RAD2DEG 57.29578

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

int main(){

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

    float ax, ay, az;
    float gx, gy, gz;
    float mx, my, mz;

	MatrixInv<float> initial_state(15, 1);
	MatrixInv<float> sensor_meas(9, 1);
	MatrixInv<float> state_sensor_val(6, 1);
	MatrixInv<float> process_noise_q(15, 15, "eye");
	MatrixInv<float> meas_noise_r(7, 7, "eye");

	MatrixInv<float> current_state(15, 1);
	MatrixInv<float> state_jacobian(15, 15);

	//Compute initial heading
	float magnetic_declination = 13.01/RAD2DEG;
	/*sensor->read_magnetometer(&mx, &my, &mz);
	initial_state(2) = atan2( -my, mx) + magnetic_declination/RAD2DEG;
	printf("My = %+7.3f, Mx = %+7.3f\n", my, mx);
	printf("Initial Heading: %+7.3f\n", initial_state(2));*/
	// process_noise_q = process_noise_q*0.1;
	// process_noise_q(0, 0) = 0.01;
	// process_noise_q(1, 1) = 0.01;
	// process_noise_q(2, 2) = 0.01;

	// process_noise_q(6, 6) = 1;
	// process_noise_q(7, 7) = 1;
	// process_noise_q(8, 8) = 1;
	// process_noise_q(9, 9) = 1;
	// process_noise_q(10, 10) = 1;
	// process_noise_q(11, 11) = 1;

	// process_noise_q(12, 12) = 0.01;
	// process_noise_q(12, 12) = 0.01;
	// process_noise_q(12, 12) = 0.01;

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
	    sum_mz += -mz;

		usleep(10000);
	}

	float avg_acc_mag = sqrt(pow((sum_ax/200), 2) + pow((sum_ay/200), 2) + pow((sum_az/200), 2));
	float avg_gaus_mag = sqrt(pow((sum_mx/200), 2) + pow((sum_my/200), 2) + pow((sum_mz/200), 2));

	initial_state(0) = atan2(sum_ay/200, avg_acc_mag);
	initial_state(1) = atan2(sum_ax/200, avg_acc_mag);

	MatrixInv<float> mag2d_projection(2, 3);
	mag2d_projection(0, 0) = cos(initial_state(1));
	mag2d_projection(0, 1) = sin(initial_state(1))*sin(initial_state(0));
	mag2d_projection(0, 2) =  sin(initial_state(1))*cos(initial_state(0));

	mag2d_projection(1, 1)  = cos(initial_state(1));
	mag2d_projection(1, 2) = -sin(initial_state(0));

	MatrixInv<float> mag_vector(3, 1);
	mag_vector(0) = sum_mx/200/avg_gaus_mag;
	mag_vector(1) = sum_my/200/avg_gaus_mag;
	mag_vector(2) = sum_mz/200/avg_gaus_mag;

	MatrixInv<float> mag2d = mag2d_projection*mag_vector;
	initial_state(2) = -atan2(mag2d(1), mag2d(0)) + magnetic_declination;

	// initial_state(2) = atan2( -(sum_my/200/avg_gaus_mag)*cos(initial_state(0)) + (sum_mz/200/avg_gaus_mag)*sin(initial_state(0)), (sum_mx/200/avg_gaus_mag)*cos(initial_state(1)) + 
	// 	((sum_my/200/avg_gaus_mag)*sin(initial_state(0)) + (sum_mz/200/avg_gaus_mag)*cos(initial_state(0)))*sin(initial_state(1)) ) + magnetic_declination;

	printf("Initial Roll: %g [deg], Initial Pitch: %g [deg], Initial Yaw: %g [deg]\n", initial_state(0)*RAD2DEG, initial_state(1)*RAD2DEG, initial_state(2)*RAD2DEG);

	usleep(5000000);

	Ekf15Dof<float> imu_gps_ekf(0.01, initial_state, process_noise_q*0.00001, meas_noise_r*20);

	// // file pointer
    fstream fout;
  
 //    // opens an existing csv file or creates a new file.
    fout.open("imu.csv", ios::out | ios::app);

	size_t count = 0;
    while(count < 100000) {
/*    	std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
*/
	    sensor->update();
	    sensor->read_accelerometer(&ax, &ay, &az);
	    sensor->read_gyroscope(&gx, &gy, &gz);
	    sensor->read_magnetometer(&mx, &my, &mz);
	   //  //printf("Acc: %+7.3f %+7.3f %+7.3f  \n", ax, ay, az);
	   //  //printf("Gyr: %+8.3f %+8.3f %+8.3f  \n", gy, gx, -gz);
	   // //printf("Mag: %+7.3f %+7.3f %+7.3f\n", mx, my, mz);

	    // fout << gy << ", "
     //         << gx << ", "
     //         << -gz << ", "
     //         << ay << ", "
     //         << ax << ", "
     //         << az
     //         << "\n";
	    state_sensor_val(0) = gy;
	    state_sensor_val(1) = gx;
	    state_sensor_val(2) = -gz;
	    state_sensor_val(3) = ay;
	    state_sensor_val(4) = ax;
	    state_sensor_val(5) = -az;

	    sensor_meas(0) = mx;
	    sensor_meas(1) = my;
	    sensor_meas(2) = -mz;
	    imu_gps_ekf.Run(state_sensor_val, sensor_meas);
	    current_state = imu_gps_ekf.GetCurrentState();
	    state_jacobian = imu_gps_ekf.GetStateJacobian();
	    //current_state(0)*RAD2DEG

	    printf("Roll [deg]: %+7.3f, Pitch[deg]: %+7.3f, Yaw[deg]: %+7.3f\n", current_state(0)*RAD2DEG, current_state(1)*RAD2DEG, current_state(2)*RAD2DEG);
	    usleep(10000);
	      fout << gy << ", "
             << gx << ", "
             << -gz << ", "
             << ay << ", "
             << ax << ", "
             << az << ","
             << my << ", "
             << mx << ", "
             << -mz << ", "
             << current_state(0) << ", "
             << current_state(1) << ", "
             << current_state(2)
             << "\n";

 //         	this->state_jacobian_(1, 0) = dt*(-s_phi*(state_sensor_val(1) - previous_state(4)) - c_phi*(state_sensor_val(2) - previous_state(5)) );
	// this->state_jacobian_(1, 1) = 1;

	// this->state_jacobian_(1, 4) = -dt*c_phi;
	// this->state_jacobian_(1, 5) = dt*s_phi;
 //             this->state_jacobian_(2, 0) = dt*( c_phi*sc_theta*(state_sensor_val(1) - previous_state(4)) - s_phi*sc_theta*(state_sensor_val(2) - previous_state(5)) );
	// this->state_jacobian_(2, 1) = dt*( s_phi*sc_theta*t_theta*(state_sensor_val(1) - previous_state(4)) + c_phi*sc_theta*t_theta*(state_sensor_val(2) - previous_state(5)) );
	// this->state_jacobian_(2, 2) = 1;

	// this->state_jacobian_(2, 4) = -dt*s_phi*sc_theta;
	// this->state_jacobian_(2, 5) = -dt*c_phi*sc_theta;
	    /*std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
	    std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[µs]" << std::endl;*/
	    count++;
	}

	fout.flush();


	return 0;
}