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

void ImuHelper::SetGyroOffset(size_t num_samples){
	//---------------------- Calculate the offset -----------------------------
    float offset[3] = {0.0, 0.0, 0.0};

    //-------------------------------------------------------------------------
    for(size_t idx = 0; idx < num_samples; idx++)
    {
        imu_sensor_->read_gyroscope(gyro_);

        offset[0] += gyro_[0];
        offset[1] += gyro_[1];
        offset[2] += gyro_[2];

        usleep(10000);
    }
    offset[0]/=static_cast<float>(num_samples);
    offset[1]/=static_cast<float>(num_samples);
    offset[2]/=static_cast<float>(num_samples);

    gyro_offset_[0] = offset[0];
    gyro_offset_[1] = offset[1];
    gyro_offset_[2] = offset[2];
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
		usleep(10);
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
	init_att_[2] = atan2(-mag2d(1), mag2d(0)) + MAG_DEC_;

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
										  {     0    , 			cos(roll)  ,       -sin(roll)     }   };
	return mag2d_projection;
}

float* ImuHelper::GetImuData(){
	ReadRawImu();

	MatrixInv<float> corrected_mag = CorrectMagData( MatrixInv<float>({ {mag_[0]}, 
																	  {mag_[1]},
																	  {mag_[2]} }) );
	for(size_t imu_idx = 0; imu_idx < 3; imu_idx++){
		imu_data_[imu_idx] = gyro_[imu_idx];
		imu_data_[imu_idx + 3] = accel_[imu_idx];
		imu_data_[imu_idx + 6] = corrected_mag(imu_idx);
	}

	return imu_data_;
}

void ImuHelper::MahonyFilterQuat(const float* imu_data, float twoKi, float twoKp, float sample_time_s, float (&quat)[4]){
	// Use Gyro and Accel only if magnetometer measurement invalid (avoids NaN in magnetometer normalisation)
    if((imu_data_[6] == 0.0f) && (imu_data_[7] == 0.0f) && (imu_data_[8] == 0.0f)) {
        MahonyFilter6Dof(imu_data, twoKi, twoKp, sample_time_s, quat);
        return;
    }else if(!((imu_data_[3] == 0.0f) && (imu_data_[4] == 0.0f) && (imu_data_[5] == 0.0f))){
    	MahonyFilter9Dof(imu_data_, twoKi, twoKp, sample_time_s, quat);
    	return;
    }
}

void ImuHelper::MahonyFilter9Dof(const float* imu_data, float twoKi, float twoKp, float sample_time_s, float (&quat)[4]){
	static float integralFBx = 0;
	static float integralFBy = 0;
	static float integralFBz = 0;
	/* Get the last time step quaternion
	*/
	float q0, q1, q2, q3;
	q0 = quat[0];
	q1 = quat[1];
	q2 = quat[2];
	q3 = quat[3];

	/* Pre process accel and gyro data
	*/
	float ax, ay, az;
	float gx, gy, gz;
	float mx, my, mz;

	ax = imu_data_[3]/G_SI;
	ay = imu_data_[4]/G_SI;
	az = imu_data_[5]/G_SI;

	gx = imu_data_[0] - gyro_offset_[0];
	gy = imu_data_[1] - gyro_offset_[1];
	gz = imu_data_[2] - gyro_offset_[2];

	mx = imu_data_[6];
	my = imu_data_[7];
	mz = imu_data_[8];

	// Normalise accelerometer measurement
    float recipNorm = InvSqrt(ax * ax + ay * ay + az * az);
    ax *= recipNorm;
    ay *= recipNorm;
    az *= recipNorm;

    // Normalise magnetometer measurement
    recipNorm = InvSqrt(mx * mx + my * my + mz * mz);
    mx *= recipNorm;
    my *= recipNorm;
    mz *= recipNorm;

    // Auxiliary variables to avoid repeated arithmetic
    float q0q0 = q0 * q0;
    float q0q1 = q0 * q1;
    float q0q2 = q0 * q2;
    float q0q3 = q0 * q3;
    float q1q1 = q1 * q1;
    float q1q2 = q1 * q2;
    float q1q3 = q1 * q3;
    float q2q2 = q2 * q2;
    float q2q3 = q2 * q3;
    float q3q3 = q3 * q3;

    // Reference direction of Earth's magnetic field
    float hx = 2.0f * (mx * (0.5f - q2q2 - q3q3) + my * (q1q2 - q0q3) + mz * (q1q3 + q0q2));
    float hy = 2.0f * (mx * (q1q2 + q0q3) + my * (0.5f - q1q1 - q3q3) + mz * (q2q3 - q0q1));
    float bx = sqrt(hx * hx + hy * hy);
    float bz = 2.0f * (mx * (q1q3 - q0q2) + my * (q2q3 + q0q1) + mz * (0.5f - q1q1 - q2q2));

    // Estimated direction of gravity and magnetic field
    float halfvx = q1q3 - q0q2;
    float halfvy = q0q1 + q2q3;
    float halfvz = q0q0 - 0.5f + q3q3;
    float halfwx = bx * (0.5f - q2q2 - q3q3) + bz * (q1q3 - q0q2);
    float halfwy = bx * (q1q2 - q0q3) + bz * (q0q1 + q2q3);
    float halfwz = bx * (q0q2 + q1q3) + bz * (0.5f - q1q1 - q2q2);

    // Error is sum of cross product between estimated direction and measured direction of field vectors
    float halfex = (ay * halfvz - az * halfvy) + (my * halfwz - mz * halfwy);
    float halfey = (az * halfvx - ax * halfvz) + (mz * halfwx - mx * halfwz);
    float halfez = (ax * halfvy - ay * halfvx) + (mx * halfwy - my * halfwx);

    // Compute and apply integral feedback if enabled
    if(twoKi > 0.0f) {
        integralFBx += twoKi * halfex * sample_time_s;	// integral error scaled by Ki
        integralFBy += twoKi * halfey * sample_time_s;
        integralFBz += twoKi * halfez * sample_time_s;
        gx += integralFBx;	// apply integral feedback
        gy += integralFBy;
        gz += integralFBz;
    }
    else {
        integralFBx = 0.0f;	// prevent integral windup
        integralFBy = 0.0f;
        integralFBz = 0.0f;
    }

    // Apply proportional feedback
    gx += twoKp * halfex;
    gy += twoKp * halfey;
    gz += twoKp * halfez;

    // Integrate rate of change of quaternion
    gx *= (0.5f * sample_time_s);		// pre-multiply common factors
    gy *= (0.5f * sample_time_s);
    gz *= (0.5f * sample_time_s);
    float qa = q0;
    float qb = q1;
    float qc = q2;
    q0 += (-qb * gx - qc * gy - q3 * gz);
    q1 += (qa * gx + qc * gz - q3 * gy);
    q2 += (qa * gy - qb * gz + q3 * gx);
    q3 += (qa * gz + qb * gy - qc * gx);

    // Normalise quaternion
    recipNorm = InvSqrt(q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3);
    q0 *= recipNorm;
    q1 *= recipNorm;
    q2 *= recipNorm;
    q3 *= recipNorm;

    quat[0] = q0;
    quat[1] = q1;
    quat[2] = q2;
    quat[3] = q3;
}

void ImuHelper::MahonyFilter6Dof(const float* imu_data, float twoKi, float twoKp, float sample_time_s, float (&quat)[4]){
	static float integralFBx = 0;
	static float integralFBy = 0;
	static float integralFBz = 0;

	/* Get the last time step quaternion
	*/
	float q0, q1, q2, q3;
	q0 = quat[0];
	q1 = quat[1];
	q2 = quat[2];
	q3 = quat[3];

	/* Pre process accel and gyro data
	*/
	float ax, ay, az;
	float gx, gy, gz;
	float mx, my, mz;

	ax = imu_data_[3];///G_SI;
	ay = imu_data_[4];///G_SI;
	az = imu_data_[5];///G_SI;

	gx = imu_data_[0] - gyro_offset_[0];
	gy = imu_data_[1] - gyro_offset_[1];
	gz = imu_data_[2] - gyro_offset_[2];

	// Normalise accelerometer measurement
    float recipNorm = InvSqrt(ax * ax + ay * ay + az * az);
    ax *= recipNorm;
    ay *= recipNorm;
    az *= recipNorm;

    // Estimated direction of gravity and vector perpendicular to magnetic flux
    float halfvx = q1 * q3 - q0 * q2;
    float halfvy = q0 * q1 + q2 * q3;
    float halfvz = q0 * q0 - 0.5f + q3 * q3;

    // Error is sum of cross product between estimated and measured direction of gravity
    float halfex = (ay * halfvz - az * halfvy);
    float halfey = (az * halfvx - ax * halfvz);
    float halfez = (ax * halfvy - ay * halfvx);

    // Compute and apply integral feedback if enabled
    if(twoKi > 0.0f) {
        integralFBx += twoKi * halfex * sample_time_s;	// integral error scaled by Ki
        integralFBy += twoKi * halfey * sample_time_s;
        integralFBz += twoKi * halfez * sample_time_s;
        gx += integralFBx;	// apply integral feedback
        gy += integralFBy;
        gz += integralFBz;
    }
    else {
        integralFBx = 0.0f;	// prevent integral windup
        integralFBy = 0.0f;
        integralFBz = 0.0f;
    }

    // Apply proportional feedback
    gx += twoKp * halfex;
    gy += twoKp * halfey;
    gz += twoKp * halfez;

    // Integrate rate of change of quaternion
    gx *= (0.5f * sample_time_s);		// pre-multiply common factors
    gy *= (0.5f * sample_time_s);
    gz *= (0.5f * sample_time_s);
    float qa = q0;
    float qb = q1;
    float qc = q2;
    q0 += (-qb * gx - qc * gy - q3 * gz);
    q1 += (qa * gx + qc * gz - q3 * gy);
    q2 += (qa * gy - qb * gz + q3 * gx);
    q3 += (qa * gz + qb * gy - qc * gx);

    // Normalise quaternion
    recipNorm = InvSqrt(q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3);
    q0 *= recipNorm;
    q1 *= recipNorm;
    q2 *= recipNorm;
    q3 *= recipNorm;
    quat[0] = q0;
    quat[1] = q1;
    quat[2] = q2;
    quat[3] = q3;
}

float ImuHelper::InvSqrt(float x)
{
    float halfx = 0.5f * x;
    float y = x;
    long i = *(long*)&y;
    i = 0x5f3759df - (i>>1);
    y = *(float*)&i;
    y = y * (1.5f - (halfx * y * y));
    return y;
}

void ImuHelper::GetEuler(float (&quat)[4], float (&mh_euler)[3]){
	float q0 = quat[0];
	float q1 = quat[1];
	float q2 = quat[2];
	float q3 = quat[3];
	mh_euler[0] = -atan2(2*(q0*q1+q2*q3), -1+2*(q1*q1+q2*q2)) * 180.0/PI;
   	mh_euler[1] = -asin(2*(q0*q2-q3*q1)) * 180.0/PI;
   	mh_euler[2] = atan2(2*(q0*q3+q1*q2), 1-2*(q2*q2+q3*q3)) * 180.0/PI;
}