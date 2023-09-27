#ifndef COORDTRANSFORM_H
#define COORDTRANSFORM_H

#include <Matrix/matrix_inv_class.h>
#include<cmath>


MatrixInv<float> Geodetic2Ecef(float lat, float lon, float height);
MatrixInv<float> Geodetic2Ned(float lat, float lon, float height, float lat_ref, float lon_ref, float height_ref);
MatrixInv<float> GetDcm(float roll, float pitch, float yaw);
void GetEulerFromQuat(const float (&quat)[4], float (&euler)[3]);
void GetQuatFromEuler(const float (&euler)[3], float (&quat)[4]);
#endif