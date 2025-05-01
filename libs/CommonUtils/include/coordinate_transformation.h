#ifndef COORDTRANSFORM_H
#define COORDTRANSFORM_H

#include <Matrix/matrix_inv_class.h>
#include<cmath>
#include<array>

MatrixInv<float> Geodetic2Ecef(float lat, float lon, float height);
MatrixInv<float> Geodetic2Ned(float lat, float lon, float height, float lat_ref, float lon_ref, float height_ref);
MatrixInv<float> GetDcm(float roll, float pitch, float yaw);
void GetEulerFromQuat(const float (&quat)[4], float (&euler)[3]);
void GetQuatFromEuler(const float (&euler)[3], float (&quat)[4]);

template <typename T, size_t M, size_t N, size_t P>
std::array<std::array<T, P>, M> MatMult(const std::array<std::array<T, N>, M>& mat1, const std::array<std::array<T, P>, N>& mat2){
	std::array<std::array<T, P>, M> result;
	for(size_t idx = 0; idx < M; idx++){
		for(size_t jdx = 0; jdx < P; jdx++){
			result[idx][jdx] = 0;
			for(size_t kdx = 0; kdx < N; kdx++){
				result[idx][jdx] += mat1[idx][kdx]*mat2[kdx][jdx];
			}
		}
	}

	return result;
}

#endif