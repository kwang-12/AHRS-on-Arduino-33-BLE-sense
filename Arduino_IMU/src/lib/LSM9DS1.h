#include <Arduino.h>
#include <Wire.h>
class simpleLSM9DS1
{
public:
    simpleLSM9DS1(TwoWire& wire);
    bool ini();
    
    bool readAccel(double& x, double& y, double& z);
    bool readGyro(double& x, double& y, double& z);
    bool readMag(double& x, double& y, double& z);

private:
    TwoWire* _wire;
    int readRegister(uint8_t slaveAddress, uint8_t address);
    int readRegisters(uint8_t slaveAddress, uint8_t address, uint8_t* data, size_t length);
    int writeRegister(uint8_t slaveAddress, uint8_t address, uint8_t value);
};