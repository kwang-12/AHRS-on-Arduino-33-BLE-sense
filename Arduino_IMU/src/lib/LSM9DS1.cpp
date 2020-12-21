#include "LSM9DS1.h"

#define LSM9DS1_ADDRESS 0x6b

#define LSM9DS1_WHO_AM_I 0x0f
#define LSM9DS1_CTRL_REG1_G 0x10
#define LSM9DS1_STATUS_REG 0x17
#define LSM9DS1_OUT_X_G 0x18
#define LSM9DS1_CTRL_REG6_XL 0x20
#define LSM9DS1_CTRL_REG8 0x22
#define LSM9DS1_OUT_X_XL 0x28

// magnetometer
#define LSM9DS1_ADDRESS_M 0x1e

#define LSM9DS1_CTRL_REG1_M 0x20
#define LSM9DS1_CTRL_REG2_M 0x21
#define LSM9DS1_CTRL_REG3_M 0x22
#define LSM9DS1_STATUS_REG_M 0x27
#define LSM9DS1_OUT_X_L_M 0x28

simpleLSM9DS1::simpleLSM9DS1(TwoWire &wire) { _wire = &wire; }

bool simpleLSM9DS1::ini()
{
    _wire->begin();
    _wire->setClock(400000);
    // reset IMU
    writeRegister(LSM9DS1_ADDRESS, LSM9DS1_CTRL_REG8, 0x05);
    writeRegister(LSM9DS1_ADDRESS_M, LSM9DS1_CTRL_REG2_M, 0x0c);
    delay(10);
    if (readRegister(LSM9DS1_ADDRESS, LSM9DS1_WHO_AM_I) != 0x68)
    {
        return false;
    }
    if (readRegister(LSM9DS1_ADDRESS_M, LSM9DS1_WHO_AM_I) != 0x3d)
    {
        return false;
    }
    // IMU settings
    writeRegister(LSM9DS1_ADDRESS, LSM9DS1_CTRL_REG1_G, 0xCB);  // 952 Hz, 500 dps, 100 Hz BW
    writeRegister(LSM9DS1_ADDRESS, LSM9DS1_CTRL_REG6_XL, 0xD0); // 952 Hz, 4G, 408 Hz BW
    writeRegister(LSM9DS1_ADDRESS_M, LSM9DS1_CTRL_REG1_M, 0xFE); // Temperature compensation enable, ultra-high performance, 80 Hz
    writeRegister(LSM9DS1_ADDRESS_M, LSM9DS1_CTRL_REG2_M, 0x00); // 4 Gauss
    writeRegister(LSM9DS1_ADDRESS_M, LSM9DS1_CTRL_REG3_M, 0x00); // Continuous conversion mode
    return true;
}

bool simpleLSM9DS1::readAccel(double &x, double &y, double &z)
{
    int16_t data[3];
    if (!readRegisters(LSM9DS1_ADDRESS, LSM9DS1_OUT_X_XL, (uint8_t *)data, sizeof(data)))
    {
        x = NAN;
        y = NAN;
        z = NAN;
        return false;
    }
    x = data[0] * 4.0 / 32768.0 * 9.8; // m/s2
    y = data[1] * 4.0 / 32768.0 * 9.8 * -1;
    z = data[2] * 4.0 / 32768.0 * 9.8;
    return true;
}

bool simpleLSM9DS1::readGyro(double &x, double &y, double &z)
{
    int16_t data[3];
    if (!readRegisters(LSM9DS1_ADDRESS, LSM9DS1_OUT_X_G, (uint8_t *)data, sizeof(data)))
    {
        x = NAN;
        y = NAN;
        z = NAN;
        return false;
    }
    x = data[0] * 0.0175 * 0.01745329; // rad/s
    y = data[1] * 0.0175 * 0.01745329 * -1;
    z = data[2] * 0.0175 * 0.01745329;
    return true;
}

bool simpleLSM9DS1::readMag(double &x, double &y, double &z)
{
    int16_t data[3];
    if (!readRegisters(LSM9DS1_ADDRESS_M, LSM9DS1_OUT_X_L_M, (uint8_t *)data, sizeof(data)))
    {
        x = NAN;
        y = NAN;
        z = NAN;
        return false;
    }
    x = data[0] * 0.014 * -1; // micro tesla
    y = data[1] * 0.014 * -1;
    z = data[2] * 0.014;
    return true;
}

int simpleLSM9DS1::readRegister(uint8_t slaveAddress, uint8_t address)
{
    _wire->beginTransmission(slaveAddress);
    _wire->write(address);
    if (_wire->endTransmission() != 0)
    {
        return -1;
    }

    if (_wire->requestFrom(slaveAddress, 1) != 1)
    {
        return -1;
    }

    return _wire->read();
}
int simpleLSM9DS1::readRegisters(uint8_t slaveAddress, uint8_t address, uint8_t *data, size_t length)
{
    _wire->beginTransmission(slaveAddress);
    _wire->write(0x80 | address);
    if (_wire->endTransmission(false) != 0)
    {
        return -1;
    }

    if (_wire->requestFrom(slaveAddress, length) != length)
    {
        return 0;
    }

    for (size_t i = 0; i < length; i++)
    {
        *data++ = _wire->read();
    }

    return 1;
}
int simpleLSM9DS1::writeRegister(uint8_t slaveAddress, uint8_t address, uint8_t value)
{
    _wire->beginTransmission(slaveAddress);
    _wire->write(address);
    _wire->write(value);
    if (_wire->endTransmission() != 0)
    {
        return 0;
    }

    return 1;
}

// simpleLSM9DS1 IMU_9DOF(Wire1);