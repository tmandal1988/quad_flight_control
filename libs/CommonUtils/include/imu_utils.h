#ifndef IMUUTILS_H
#define IMUUTILS_H

// Navio 2 Utilities
#include "Navio/Common/MPU9250.h"
#include "Navio/Navio2/LSM9DS1.h"
#include<Matrix/matrix_inv_class.h>
#include<constants.h>

#include<cmath>
#include<memory>

using namespace std;

class ImuHelper{
	public:
		ImuHelper(const string& imu_name);
		~ImuHelper();
		void GetInertialSensor();
		void InitializeImu();
		void SetGyroOffset(size_t num_samples);
		float* ComputeInitialRollPitchAndYaw(size_t num_samples);
		MatrixInv<float> CorrectMagData(MatrixInv<float> mag_vector);
		MatrixInv<float> GetMag3DTo2DProj(float roll, float pitch);
		float* GetImuData();
		void SetMagParams(float mag_dec, MatrixInv<float> mag_a, MatrixInv<float> mag_offset, MatrixInv<float> mag_scale){
			MAG_DEC_ = mag_dec;
			MAG_A_ = mag_a;
			MAG_OFFSET_ = mag_offset;
			MAG_SCALE_ = mag_scale;
		}
		void MahonyFilterQuat(const float* imu_data, float twoKi, float twoKp, float sample_time_s, float (&quat)[4]);
		void MahonyFilter9Dof(const float* imu_data, float twoKi, float twoKp, float sample_time_s, float (&quat)[4]);
		void MahonyFilter6Dof(const float* imu_data, float twoKi, float twoKp, float sample_time_s, float (&quat)[4]);
		void GetEuler(float (&quat)[4], float (&mh_euler)[3]);
		float InvSqrt(float x);

	private:
		unique_ptr<InertialSensor> imu_sensor_;
		string imu_name_;

		// Variables to store data from the IMU
    	// Accels
		float accel_[3];
    	//Gyros
		float gyro_[3];
    	//Mags
		float mag_[3];

		// Gyro offsets
		float gyro_offset_[3];

		// Initial roll, pitch and yaw
		float init_att_[3];

		// IMU data array
		float imu_data_[9];

		// Mag parameters
		float MAG_DEC_;
		MatrixInv<float> MAG_A_;
		MatrixInv<float> MAG_SCALE_;
		MatrixInv<float> MAG_OFFSET_;

		void ReadRawImu();
		
};

#endif