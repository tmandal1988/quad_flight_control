#include "coordinate_transformation.h"

MatrixInv<float> Geodetic2Ecef(float lat, float lon, float height){
	/* This function converts Geodetic positions to ECEF positions
	Inputs:
		lat (float)		: Lattitude measured by GPS
		lon (float)		: Longitude measured by GPS
		height (float)	: Height above ellipsoid measured by GPS
	Outputs:
		MatrixInv<float>: 3x1 array of ECEF X, Y, Z coordinates
	*/

	// Square of Earth's eccentricity
	float ecc_sq = 0.0066943798522561;
	// Earth's semi measure axis
    float r_ea = 6378137.0;
    // Intermediate variable for coordinate conversion
    float n_lat = r_ea/sqrt( ( 1 - ecc_sq * pow(sin(lat), 2) ) );

    float c_lat = cos(lat);
    float c_lon = cos(lon);
    float s_lat = sin(lat);
    float s_lon = sin(lon);

    // Transform the coordinates and return
    MatrixInv<float> ecef_cord = { {(n_lat + height)*c_lat*c_lon}, {(n_lat + height)*c_lat*s_lon}, {(n_lat*(1 - ecc_sq) + height)*s_lat} };
    return ecef_cord;
}

MatrixInv<float> Geodetic2Ned(float lat, float lon, float height, float lat_ref, float lon_ref, float height_ref){
	/* This function converts Geodetic positions to NED positions
	Inputs:
		lat (float)							: Lattitude measured by GPS
		lon (float)							: Longitude measured by GPS
		height (float)						: Height above ellipsoid measured by GPS
		lat_ref (float)						: Lattitude of the origin of the NED frame
		lon_ref (float)						: Longitude of the origin of the NED frame
		height_ref (float)					: Height above ellipsoid of the origin of the NED frame
	Outputs:
		MatrixInv<float> 	: 3x1 array of N, E, D coordinates
	*/

	// Convert current llh to ECEF coordinates
	MatrixInv<float> ecef_cord = Geodetic2Ecef(lat, lon, height);
	// Convert reference 
	MatrixInv<float> ecef_cord_ref = Geodetic2Ecef(lat_ref, lon_ref, height_ref);

	float c_lat_ref = cos(lat_ref);
	float s_lat_ref = sin(lat_ref);

	float c_lon_ref = cos(lon_ref);
	float s_lon_ref = sin(lon_ref);

	// Create ECEF to NED transformation
    MatrixInv<float> ecef2ned = { {-s_lat_ref*c_lon_ref, -s_lat_ref*s_lon_ref, c_lat_ref}, 
    							  {-s_lon_ref, c_lon_ref, 0}, 
    							  {-c_lat_ref*c_lon_ref, -c_lat_ref*s_lon_ref, -s_lat_ref} };
   // Return NED coordinates
   	return ecef2ned*(ecef_cord - ecef_cord_ref);
}

MatrixInv<float> GetDcm(float roll, float pitch, float yaw){
	/* This function computes the DCM for a given roll, pitch and yaw angle. The DCM takes from Inertial to Body coordinates.
	Inputs:
		roll (float) 						: Roll angle of the body w.r.t inertial frame
		pitch (float) 						: Pitch angle of the body w.r.t inertial frame
		yaw (float) 						: Yaw angle of the body w.r.t inertial frame
	Outputs:
		MatrixInv<float> 					: 3x3 DCM
	*/

	// precalculate trignometric values
	float s_phi = sin(roll);
	float c_phi = cos(roll);

	float s_theta = sin(pitch);
	float c_theta = cos(pitch);
	float t_theta = tan(pitch);
	float sc_theta = 1/c_theta;

	float s_psi = sin(yaw);
	float c_psi = cos(yaw);

	MatrixInv<float> c_ned2b = { {c_psi*c_theta, c_theta*s_psi, -s_theta},
							 {c_psi*s_phi*s_theta - c_phi*s_psi, c_phi*c_psi + s_phi*s_psi*s_theta, c_theta*s_phi},
							 {s_phi*s_psi + c_phi*c_psi*s_theta, c_phi*s_psi*s_theta - c_psi*s_phi, c_phi*c_theta} };

	return c_ned2b;
}

void GetEulerFromQuat(const float (&quat)[4], float (&euler)[3]){
	euler[0] = -atan2( 2*( quat[0]*quat[1] + quat[2]*quat[3] ), -1 + 2*( quat[1]*quat[1] + quat[2]*quat[2] ) );
   	euler[1] = -asin( 2*( quat[0]*quat[2] - quat[3]*quat[1] ) );
   	euler[2] = atan2( 2*( quat[0]*quat[3] + quat[1]*quat[2] ), 1-2*( quat[2]*quat[2] + quat[3]*quat[3] ) );
}

void GetQuatFromEuler(const float (&euler)[3], float (&quat)[4]){
	quat[0] = cos(euler[2]/2)*cos(euler[1]/2)*cos(euler[0]/2) + sin(euler[2]/2)*sin(euler[1]/2)*sin(euler[0]/2);
  	quat[1] = cos(euler[2]/2)*cos(euler[1]/2)*sin(euler[0]/2) - sin(euler[2]/2)*sin(euler[1]/2)*cos(euler[0]/2);
  	quat[2] = cos(euler[2]/2)*sin(euler[1]/2)*cos(euler[0]/2) + sin(euler[2]/2)*cos(euler[1]/2)*sin(euler[0]/2);
  	quat[3] = sin(euler[2]/2)*cos(euler[1]/2)*cos(euler[0]/2) - cos(euler[2]/2)*sin(euler[1]/2)*sin(euler[0]/2);
}