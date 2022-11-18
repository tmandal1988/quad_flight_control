#ifndef _INERTIAL_SENSOR_H
#define _INERTIAL_SENSOR_H

class InertialSensor {
public:
    virtual bool initialize() = 0;
    virtual bool probe() = 0;
    virtual void update() = 0;

    float read_temperature() {return temperature;};
    void read_accelerometer(float (&accel)[3]) {accel[1] = _ax; accel[0] = _ay; accel[2] = -_az;};
    void read_gyroscope(float (&gyro)[3]) {gyro[1] = _gx; gyro[0] = _gy; gyro[2] = -_gz;};
    void read_magnetometer(float (&mag)[3]) {mag[0] = _mx; mag[1] = _my; mag[2] = _mz;};

protected:
    float temperature;
    float _ax, _ay, _az;
    float _gx, _gy, _gz;
    float _mx, _my, _mz;
};

#endif //_INERTIAL_SENSOR_H
