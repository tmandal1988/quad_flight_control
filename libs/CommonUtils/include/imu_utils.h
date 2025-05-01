#ifndef IMUUTILS_H
#define IMUUTILS_H

// Navio 2 Utilities
#include "Navio/Common/MPU9250.h"
#include "Navio/Navio2/LSM9DS1.h"
#include<Matrix/matrix_inv_class.h>
#include<constants.h>
#include "notch_filter.h"

#include<cmath>
#include<memory>

using namespace std;

class ImuHelper{
	public:
		ImuHelper(const string& imu_name);
		~ImuHelper();
		void GetInertialSensor();
		void InitializeImu();
		void UpdateImuNotchFilterCoeffs(const array<float, 3> &notch_filter_num, const array<float, 3> &notch_filter_den);
		void EnableImuNotchFilters();
		void DisableImuNotchFilters();
		void SetAccelCalibParams(const array<array<float, 3>, 3> &accel_calib_m, const std::array<float, 3> &accel_calib_off);
		void ComputeGyroOffset(size_t num_samples);
		void GetGyroOffset(float (&gyro_offset)[3]);
		float* ComputeInitialRollPitchAndYaw(size_t num_samples);
		MatrixInv<float> CorrectMagData(MatrixInv<float> mag_vector);
		MatrixInv<float> GetMag3DTo2DProj(float roll, float pitch);
		float* GetImuData();
		bool getRawImuData(float (&imu_raw_data)[9]);
		void SetMagParams(float mag_dec, MatrixInv<float> mag_a, MatrixInv<float> mag_offset, MatrixInv<float> mag_scale){
			MAG_DEC_ = mag_dec;
			MAG_A_ = mag_a;
			MAG_OFFSET_ = mag_offset;
			MAG_SCALE_ = mag_scale;
		}

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
		float imu_raw_data_[9];
		float old_mag_[3] = {0.0f, 0.0f, 0.0f};
		bool is_mag_valid_ = false;

		// Mag parameters
		float MAG_DEC_;
		MatrixInv<float> MAG_A_;
		MatrixInv<float> MAG_SCALE_;
		MatrixInv<float> MAG_OFFSET_;

		array<float, 3> notch_filter_num_{1, 0, 0};
		array<float, 3> notch_filter_den_{0, 0, 0};

		array<array<float, 3>, 3 > accel_calib_m_{0};
		array<float, 3> accel_calib_off_{0};

		bool enable_notch_filter_{false};

		NotchFilter<float> gyro_accel_notch_filters_[6];

		void ReadRawImu();
		
};

#endif