#ifndef MAHONYFILTER_H
#define MAHONYFILTER_H

// Navio 2 Utilities
#include<cmath>
#include "coordinate_transformation.h"

using namespace std;

class MahonyFilter{
	public:
		MahonyFilter(float two_ki, float two_kp, float sample_time_s, const float (&quat)[4]);
		~MahonyFilter();
		void MahonyFilterQuat(const float* imu_data, const float (&gyro_offset)[3]);
		void MahonyFilter9Dof(const float* imu_data, const float (&gyro_offset)[3]);
		void MahonyFilter6Dof(const float* imu_data, const float (&gyro_offset)[3]);
		void GetCurrentQuat(float (&quat)[4]);
		float InvSqrt(float x);

	private:
		float two_ki_;
		float two_kp_;
		float q0_; 
		float q1_; 
		float q2_; 
		float q3_;
		float sample_time_s_;
		float integral_fbx_;
		float integral_fby_;
		float integral_fbz_;
};

#endif