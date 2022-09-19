#include <ekf_21dof_class.h>
#include <Matrix/matrix_factorization_class.h>
#include "Navio/Common/MPU9250.h"
#include "Navio/Navio2/LSM9DS1.h"
#include "Navio/Common/Util.h"
#include <memory>

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

	MatrixInv<float> initial_state(21, 1);
	MatrixInv<float> sensor_meas(9, 1);
	MatrixInv<float> state_sensor_val(6, 1);
	MatrixInv<float> process_noise_q(21, 21, "eye");
	MatrixInv<float> meas_noise_r(7, 7, "eye");

	MatrixInv<float> current_state(21, 1);

	//Compute initial heading
	float magnetic_declination = 13.01;
	sensor->read_magnetometer(&mx, &my, &mz);
	initial_state(2) = atan2( -my, mx) + magnetic_declination;

	Ekf21Dof<float> imu_gps_ekf(0.01, initial_state, process_noise_q*0.01, meas_noise_r*0.001);

    while(1) {
	    sensor->update();
	    sensor->read_accelerometer(&ax, &ay, &az);
	    sensor->read_gyroscope(&gx, &gy, &gz);
	    sensor->read_magnetometer(&mx, &my, &mz);
/*	    printf("Acc: %+7.3f %+7.3f %+7.3f  ", ax, ay, az);
	    printf("Gyr: %+8.3f %+8.3f %+8.3f  ", gx, gy, gz);
	    printf("Mag: %+7.3f %+7.3f %+7.3f\n", mx, my, mz);*/

	    state_sensor_val(0) = gx;
	    state_sensor_val(1) = gy;
	    state_sensor_val(2) = gz;
	    state_sensor_val(3) = ax;
	    state_sensor_val(4) = ay;
	    state_sensor_val(5) = az;

	    sensor_meas(0) = mx;
	    sensor_meas(1) = my;
	    sensor_meas(2) = mz;

	    imu_gps_ekf.Run(state_sensor_val, sensor_meas);
	    current_state = imu_gps_ekf.GetCurrentState();

	    printf("Roll [deg]: %+7.3f, Pitch[deg]: %+7.3f, Yaw[deg]: %+7.3f\n", current_state(0)*RAD2DEG, current_state(1)*RAD2DEG, current_state(2)*RAD2DEG);
	    usleep(10000);
	}


	return 0;
}