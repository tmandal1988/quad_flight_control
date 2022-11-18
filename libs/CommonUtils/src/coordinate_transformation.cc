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