#ifndef EXTMAG_H
#define EXTMAG_H

#include "Navio/Common/I2Cdev.h"
 #include <unistd.h>


#include<iostream>

#define DEV_ID 						   0x10

#define IST8310_RA_WAI		           0x00
#define IST8310_RA_ST1		           0x02
#define IST8310_RA_AVGCNTL		       0x41
#define IST8310_RA_PDCNTL		       0x42
#define IST8310_RA_CNTL1		       0x0A
#define IST8310_RA_STR  		       0x0C

#define IST8310_RA_DATAXL		       0x03

using namespace std;

class ExtMagDriver{
	public:
		ExtMagDriver();
		~ExtMagDriver();
		int TestConnection();
		int StartMeasurement();
		int GetMeasurement(int16_t mag_data []);

	private:
		uint8_t dev_addr_ = 0x0E;

};

#endif