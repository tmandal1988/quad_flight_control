#include "mahony_filter.h"
#include "coordinate_transformation.h"

MahonyFilter::MahonyFilter(float two_ki, float two_kp, float sample_time_s, const float (&quat)[4]):
							two_ki_(two_ki),
							two_kp_(two_kp),
							sample_time_s_(sample_time_s){

	q0_ = quat[0];
	q1_ = quat[1];
	q2_ = quat[2];
	q3_ = quat[3];	
	
	float recip_norm = InvSqrt(q0_ * q0_ + q1_ * q1_ + q2_ * q2_ + q3_ * q3_);
    q0_ *= recip_norm;
    q1_ *= recip_norm;
    q2_ *= recip_norm;
    q3_ *= recip_norm;
		
	integral_fbx_ = 0;
	integral_fby_ = 0;
	integral_fbz_ = 0;				
}

MahonyFilter::~MahonyFilter(){

}

void MahonyFilter::MahonyFilterQuat(const float* imu_data, const float (&gyro_offset)[3]){
	// Use Gyro and Accel only if magnetometer measurement invalid (avoids NaN in magnetometer normalisation)
    if((imu_data[6] == 0.0f) && (imu_data[7] == 0.0f) && (imu_data[8] == 0.0f)) {
        MahonyFilter6Dof(imu_data, gyro_offset);
        return;
    }else if(!((imu_data[3] == 0.0f) && (imu_data[4] == 0.0f) && (imu_data[5] == 0.0f))){
    	MahonyFilter9Dof(imu_data, gyro_offset);
    	return;
    }
}

void MahonyFilter::MahonyFilter6Dof(const float* imu_data, const float (&gyro_offset)[3]){
	/* Pre process accel and gyro data
	*/
	float ax, ay, az;
	float gx, gy, gz;

	ax = imu_data[3];///G_SI;
	ay = imu_data[4];///G_SI;
	az = imu_data[5];///G_SI;

	gx = imu_data[0] - gyro_offset[0];
	gy = imu_data[1] - gyro_offset[1];
	gz = imu_data[2] - gyro_offset[2];

	// Normalise accelerometer measurement
    float recip_norm = InvSqrt(ax * ax + ay * ay + az * az);
    ax *= recip_norm;
    ay *= recip_norm;
    az *= recip_norm;

    // Estimated direction of gravity and vector perpendicular to magnetic flux
    float halfvx = q1_ * q3_ - q0_ * q2_;
    float halfvy = q0_ * q1_ + q2_ * q3_;
    float halfvz = q0_ * q0_ - 0.5f + q3_ * q3_;

    // Error is sum of cross product between estimated and measured direction of gravity
    float halfex = (ay * halfvz - az * halfvy);
    float halfey = (az * halfvx - ax * halfvz);
    float halfez = (ax * halfvy - ay * halfvx);

    // Compute and apply integral feedback if enabled
    if(two_ki_ > 0.0f) {
        integral_fbx_ += two_ki_ * halfex * sample_time_s_;	// integral error scaled by Ki
        integral_fby_ += two_ki_ * halfey * sample_time_s_;
        integral_fbz_ += two_ki_ * halfez * sample_time_s_;
        gx += integral_fbx_;	// apply integral feedback
        gy += integral_fby_;
        gz += integral_fbz_;
    }
    else {
        integral_fbx_ = 0.0f;	// prevent integral windup
        integral_fby_ = 0.0f;
        integral_fbz_ = 0.0f;
    }

    // Apply proportional feedback
    gx += two_kp_ * halfex;
    gy += two_kp_ * halfey;
    gz += two_kp_ * halfez;

    // Integrate rate of change of quaternion
    gx *= (0.5f * sample_time_s_);		// pre-multiply common factors
    gy *= (0.5f * sample_time_s_);
    gz *= (0.5f * sample_time_s_);
    float qa = q0_;
    float qb = q1_;
    float qc = q2_;
    q0_ += (-qb * gx - qc * gy - q3_ * gz);
    q1_ += (qa * gx + qc * gz - q3_ * gy);
    q2_ += (qa * gy - qb * gz + q3_ * gx);
    q3_ += (qa * gz + qb * gy - qc * gx);

    // Normalise quaternion
    recip_norm = InvSqrt(q0_ * q0_ + q1_ * q1_ + q2_ * q2_ + q3_ * q3_);
    q0_ *= recip_norm;
    q1_ *= recip_norm;
    q2_ *= recip_norm;
    q3_ *= recip_norm;

    // float quat[4];
    // quat[0] = q0_;
    // quat[1] = q1_;
    // quat[2] = q2_;
    // quat[3] = q3_;

    // float euler[3];
    // GetEulerFromQuat(quat, euler);
    // printf("Roll [deg]: %+7.3f, Pitch[deg]: %+7.3f, Yaw[deg]: %+7.3f\n", euler[0] * 180/3.14, euler[1]*180/3.14, euler[2]*180/3.14);

}

void MahonyFilter::MahonyFilter9Dof(const float* imu_data, const float (&gyro_offset)[3]){
	/* Pre process accel and gyro data
	*/
	float ax, ay, az;
	float gx, gy, gz;
	float mx, my, mz;

	ax = imu_data[3];///G_SI;
	ay = imu_data[4];///G_SI;
	az = imu_data[5];///G_SI;

	gx = imu_data[0] - gyro_offset[0];
	gy = imu_data[1] - gyro_offset[1];
	gz = imu_data[2] - gyro_offset[2];

	mx = imu_data[6];
	my = imu_data[7];
	mz = imu_data[8];

	// Normalise accelerometer measurement
    float recip_norm = InvSqrt(ax * ax + ay * ay + az * az);
    ax *= recip_norm;
    ay *= recip_norm;
    az *= recip_norm;

    // Normalise magnetometer measurement
    recip_norm = InvSqrt(mx * mx + my * my + mz * mz);
    mx *= recip_norm;
    my *= recip_norm;
    mz *= recip_norm;

    // Auxiliary variables to avoid repeated arithmetic
    float q0q0 = q0_ * q0_;
    float q0q1 = q0_ * q1_;
    float q0q2 = q0_ * q2_;
    float q0q3 = q0_ * q3_;
    float q1q1 = q1_ * q1_;
    float q1q2 = q1_ * q2_;
    float q1q3 = q1_ * q3_;
    float q2q2 = q2_ * q2_;
    float q2q3 = q2_ * q3_;
    float q3q3 = q3_ * q3_;

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
    if(two_ki_ > 0.0f) {
        integral_fbx_ += two_ki_ * halfex * sample_time_s_;	// integral error scaled by Ki
        integral_fby_ += two_ki_ * halfey * sample_time_s_;
        integral_fbz_ += two_ki_ * halfez * sample_time_s_;
        gx += integral_fbx_;	// apply integral feedback
        gy += integral_fby_;
        gz += integral_fbz_;
    }
    else {
        integral_fbx_ = 0.0f;	// prevent integral windup
        integral_fby_ = 0.0f;
        integral_fbz_ = 0.0f;
    }

    // Apply proportional feedback
    gx += two_kp_ * halfex;
    gy += two_kp_ * halfey;
    gz += two_kp_ * halfez;

    // Integrate rate of change of quaternion
    gx *= (0.5f * sample_time_s_);		// pre-multiply common factors
    gy *= (0.5f * sample_time_s_);
    gz *= (0.5f * sample_time_s_);
    float qa = q0_;
    float qb = q1_;
    float qc = q2_;
    q0_ += (-qb * gx - qc * gy - q3_ * gz);
    q1_ += (qa * gx + qc * gz - q3_ * gy);
    q2_ += (qa * gy - qb * gz + q3_ * gx);
    q3_ += (qa * gz + qb * gy - qc * gx);

    // Normalise quaternion
    recip_norm = InvSqrt(q0_ * q0_ + q1_ * q1_ + q2_ * q2_ + q3_ * q3_);
    q0_ *= recip_norm;
    q1_ *= recip_norm;
    q2_ *= recip_norm;
    q3_ *= recip_norm;
}

float MahonyFilter::InvSqrt(float x)
{
    float halfx = 0.5f * x;
    float y = x;
    long i = *(long*)&y;
    i = 0x5f3759df - (i>>1);
    y = *(float*)&i;
    y = y * (1.5f - (halfx * y * y));
    return y;
}

void MahonyFilter::GetCurrentQuat(float (&quat)[4]){
	quat[0] = q0_;
	quat[1] = q1_;
	quat[2] = q2_;
	quat[3] = q3_;
}

