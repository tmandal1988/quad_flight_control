#include "external_mag.h"

ExtMagDriver::ExtMagDriver(){

}

ExtMagDriver::~ExtMagDriver(){
	
}

int ExtMagDriver::TestConnection(){
	uint8_t data;
    int8_t status = I2Cdev::readByte(dev_addr_, IST8310_RA_WAI, &data);

    size_t mag_check_count = 0;
    while(mag_check_count < 5){
        if( (status > 0) && (data == DEV_ID) ){
            printf("External Mag ID check complete\n");
            break;
        }
        usleep(7000);
        mag_check_count++;
    }

    if(mag_check_count == 5){
        printf("External mag device ID check failed\n");
        return -1;
    }

    //Set Pulse Duration Control and Average Control Register,
    data = 0xC0;
    status = I2Cdev::writeByte(dev_addr_, IST8310_RA_PDCNTL, data);

    if( status < 0){
    	printf("Unable to set PDC register of the external_mag\n");
    	return -1;
    }

    data = 0x24;
    status = I2Cdev::writeByte(dev_addr_, IST8310_RA_AVGCNTL, data);

    if( status < 0){
    	printf("Unable to set AVG register of the external_mag\n");
    	return -1;
    }

    return 1;
}

int ExtMagDriver::StartMeasurement(){
	// Set up for single measurement mode
    uint8_t data = 0x01;
    int status = I2Cdev::writeByte(dev_addr_, IST8310_RA_CNTL1, data);

    return status;
}

int ExtMagDriver::GetMeasurement(int16_t mag_data []){
	uint8_t buffer[6];
	int status = I2Cdev::readBytes(dev_addr_, IST8310_RA_DATAXL, 6, buffer);
    mag_data[0] = (buffer[1] << 8) | buffer[0];
    mag_data[1] = (buffer[3] << 8) | buffer[2];
    mag_data[2] = (buffer[5] << 8) | buffer[4];

    return status;
}