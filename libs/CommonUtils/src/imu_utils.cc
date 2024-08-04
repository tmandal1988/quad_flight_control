#include "imu_utils.h"

ImuHelper::ImuHelper(const string& imu_name):imu_name_(imu_name){
}


ImuHelper::~ImuHelper(){

}

void ImuHelper::InitializeImu(){
	if (imu_name_ == "mpu") {
		printf("Selected: MPU9250\n");
		imu_sensor_ = unique_ptr <InertialSensor>{ new MPU9250() };
	}else if (imu_name_ == "lsm") {
		printf("Selected: LSM9DS1\n");
		imu_sensor_ = unique_ptr <InertialSensor>{ new LSM9DS1() };
	}else{
		std::cerr << "IMU type is not supported, valid types are lsm and mpu"<<endl;
	}

	imu_sensor_->initialize();

}

void ImuHelper::UpdateImuNotchFilterCoeffs(const array<float, 3> notch_filter_num, const array<float, 3> notch_filter_den){
	notch_filter_num_ = notch_filter_num;
	notch_filter_den_ = notch_filter_den;
	for(size_t idx = 0; idx < 6; idx++){
		gyro_accel_notch_filters_[idx].UpdateFilterCoeff(notch_filter_num_, notch_filter_den_);
	}
}

void ImuHelper::EnableImuNotchFilters(){
	enable_notch_filter_ = true;
}

void ImuHelper::DisableImuNotchFilters(){
	enable_notch_filter_ = false;
}

void ImuHelper::ComputeGyroOffset(size_t num_samples){
	//---------------------- Calculate the offset -----------------------------
    float offset[3] = {0.0, 0.0, 0.0};

    //-------------------------------------------------------------------------
    for(size_t idx = 0; idx < num_samples; idx++)
    {
        ReadRawImu();

        offset[0] += gyro_[0];
        offset[1] += gyro_[1];
        offset[2] += gyro_[2];

        usleep(5000);
    }
    offset[0]/=static_cast<float>(num_samples);
    offset[1]/=static_cast<float>(num_samples);
    offset[2]/=static_cast<float>(num_samples);

    gyro_offset_[0] = offset[0];
    gyro_offset_[1] = offset[1];
    gyro_offset_[2] = offset[2];
}

void ImuHelper::GetGyroOffset(float (&gyro_offset)[3]){
	gyro_offset[0] = gyro_offset_[0];
	gyro_offset[1] = gyro_offset_[1];
	gyro_offset[2] = gyro_offset_[2];
}

float* ImuHelper::ComputeInitialRollPitchAndYaw(size_t num_samples){
	// Variables to store running sum of accel and mag data
	float sum_accel[3] = {0, 0, 0};
	float sum_mag[3] = {0, 0, 0};

	for (size_t idx_s = 0; idx_s < num_samples; idx_s++){
		// Read IMU data and sum them up
		ReadRawImu();
		for(size_t imu_idx = 0; imu_idx < 3; imu_idx++){
			sum_accel[imu_idx] += accel_[imu_idx];
			sum_mag[imu_idx] += mag_[imu_idx];
		}
		usleep(5000);
	}

	// Compute average accel magnitude
	float avg_ax = sum_accel[0]/num_samples;
	float avg_ay = sum_accel[1]/num_samples;
	float avg_az = sum_accel[2]/num_samples;

	float avg_acc_mag = sqrt(pow(avg_ax, 2) + pow(avg_ay, 2) + pow(avg_az, 2));
	// Set initial roll and pitch angle
	init_att_[0] = atan2(avg_ay, avg_acc_mag);
	init_att_[1] = atan2(avg_ax, avg_acc_mag);

	// Compute average magn value and correct it
	float avg_mx = sum_mag[0]/num_samples;
	float avg_my = sum_mag[1]/num_samples;
	float avg_mz = sum_mag[2]/num_samples;

	MatrixInv<float> mag_vector =  { {avg_mx},
									 {avg_my},
									 {avg_mz} };

	// Correct magnetometer vector for errors
	mag_vector = CorrectMagData(mag_vector);

	// Compute heading from mag readings
	MatrixInv<float> mag2d = GetMag3DTo2DProj(init_att_[0], init_att_[1])*mag_vector;
	init_att_[2] = atan2(mag2d(1), mag2d(0)) + MAG_DEC_;

	return init_att_;
}

void ImuHelper::ReadRawImu(){
	imu_sensor_->update();
	imu_sensor_->read_accelerometer(accel_);
	imu_sensor_->read_gyroscope(gyro_);
	imu_sensor_->read_magnetometer(mag_);
}

MatrixInv<float> ImuHelper::CorrectMagData(MatrixInv<float> mag_vector){
	// Matrix to project 3D mag to 2D
	MatrixInv<float> corrected_mag;
	// Correct for mag hard and soft iron errors
	corrected_mag = MAG_A_*MAG_SCALE_*(mag_vector - MAG_OFFSET_);

	// Normalize the corrected mag vector
	float mag_norm = sqrt(pow(corrected_mag(0), 2) + pow(corrected_mag(1), 2) + pow(corrected_mag(2), 2));
	return corrected_mag/mag_norm;

}

MatrixInv<float> ImuHelper::GetMag3DTo2DProj(float roll, float pitch){
	MatrixInv<float> mag2d_projection = { {cos(pitch), sin(pitch)*sin(roll), sin(pitch)*cos(roll) },
										  {     0    , 			-cos(roll)  ,       sin(roll)     }   };
	return mag2d_projection;
}

float* ImuHelper::GetImuData(){
	ReadRawImu();

	MatrixInv<float> corrected_mag = CorrectMagData( MatrixInv<float>({ {mag_[0]}, 
																	  {mag_[1]},
																	  {mag_[2]} }) );
	for(size_t imu_idx = 0; imu_idx < 3; imu_idx++){
		if(enable_notch_filter_){
			imu_data_[imu_idx] = gyro_accel_notch_filters_[imu_idx].Filter(gyro_[imu_idx]);
			imu_data_[imu_idx + 3] = gyro_accel_notch_filters_[imu_idx + 3].Filter(accel_[imu_idx]);
		}else{
			imu_data_[imu_idx] = gyro_[imu_idx];
			imu_data_[imu_idx + 3] = accel_[imu_idx];
		}
		
		imu_data_[imu_idx + 6] = corrected_mag(imu_idx);
	}

	return imu_data_;
}