#include "src\lib\LSM9DS1.h"
#include "src\lib\matrix.h"
#include "src\lib\ekf.h"
#include <cmath>
#include <elapsedMillis.h>
#include <ArduinoBLE.h>

// timer
#define sample_mili 10
elapsedMillis timer = 0;

// BLE
BLEService imuService("5e85c96e-41d8-11eb-b378-0242ac130002");
BLEStringCharacteristic firstCharacteristic("5e85c96e-41d8-11eb-b378-0242ac130002", BLERead | BLENotify, 50);
char txBuffer[100];
bool is_firstPart = true;

// IMU
simpleLSM9DS1 IMU_9DOF(Wire1);
double xa, ya, za;
double xg, yg, zg;
double xm, ym, zm;
// EKF
Matrix P(4, 4);
Matrix Q(4, 4);
Matrix R(6, 6);
Matrix quaternion(4, 1);
Matrix Y(6, 1);
Matrix U(3, 1);
EKF ekf(quaternion, P, Q, R, double(sample_mili / 1000));

// hard iron bias - identified with motionCal.exe by PJRC
double mag_bias[3] = {8.648, 15.058, -89.069};

void setup()
{
    pinMode(LED_BUILTIN, OUTPUT);
    if (!BLE.begin())
    {
        Serial.println("Starting BLE failed");
        while (true)
        {
            delay(100);
        }
    }
    // Initialize BLE service
    BLE.setLocalName("IMU");
    BLE.setAdvertisedService(imuService);
    imuService.addCharacteristic(firstCharacteristic);
    BLE.addService(imuService);
    firstCharacteristic.writeValue(txBuffer);
    BLE.advertise();
    // Initialize IMU
    IMU_9DOF.ini();
    // read couple times to avoid nan
    IMU_9DOF.readAccel(xa, ya, za);
    IMU_9DOF.readGyro(xg, yg, zg);
    IMU_9DOF.readMag(xm, ym, zm);
    delay(100);
    IMU_9DOF.readAccel(xa, ya, za);
    IMU_9DOF.readGyro(xg, yg, zg);
    IMU_9DOF.readMag(xm, ym, zm);
    delay(100);

    // initialize EKF
    quaternion[0][0] = 1.0;
    P.setValueDiagonal(double(8));
    Q.setValueDiagonal(double(1e-6));
    R.setValueDiagonal(double(0.0001));
    ekf.reset(quaternion, P, Q, R);
}
void loop()
{
    // EKF
    if (timer >= sample_mili)
    {
        timer = 0;
        IMU_9DOF.readAccel(xa, ya, za);
        IMU_9DOF.readGyro(xg, yg, zg);
        IMU_9DOF.readMag(xm, ym, zm);
        if (isnan(xg))
        {
            xg = 0;
        }
        if (isnan(yg))
        {
            yg = 0;
        }
        if (isnan(zg))
        {
            zg = 0;
        }
        // measurement - accel and mag
        Y[0][0] = xa;
        Y[1][0] = ya;
        Y[2][0] = za;
        Y[3][0] = xm;
        Y[4][0] = ym;
        Y[5][0] = zm;
        // hard-iron bias
        Y[3][0] = Y[3][0] - mag_bias[0];
        Y[4][0] = Y[4][0] - mag_bias[1];
        Y[5][0] = Y[5][0] - mag_bias[2];
        // normalize
        double norm = sqrt(Y[0][0] * Y[0][0] + Y[1][0] * Y[1][0] + Y[2][0] * Y[2][0]);
        Y[0][0] = Y[0][0] / norm;
        Y[1][0] = Y[1][0] / norm;
        Y[2][0] = Y[2][0] / norm;
        norm = sqrt(Y[3][0] * Y[3][0] + Y[4][0] * Y[4][0] + Y[5][0] * Y[5][0]);
        Y[3][0] = Y[3][0] / norm;
        Y[4][0] = Y[4][0] / norm;
        Y[5][0] = Y[5][0] / norm;
        for (int tick = 0; tick < 6; tick++)
        {
            if (isnan(Y[tick][0]))
            {
                Serial.print(tick);
                Serial.print(" is nan ");
                Y[tick][0] = double(0.0);
            }
        }
        // if (!ekf.update(Y, U))
        if (!ekf.update(xg, yg, zg, Y))
        {
            quaternion.setValueAll(double(0.0));
            quaternion[0][0] = 1.0;
            ekf.reset(quaternion, P, Q, R);
        }
    }
    // BLE
    BLEDevice central = BLE.central();
    if (central)
    {
        if (central.connected())
        {
            // Light up built-in LED, send package
            digitalWrite(LED_BUILTIN, HIGH);
            quaternion = ekf.GetState();
            if (is_firstPart)
            {
                snprintf(txBuffer, sizeof(txBuffer) - 1, "x%.3fy%.3f\r\n", quaternion[1][0], quaternion[2][0]);
                firstCharacteristic.writeValue(txBuffer);
                is_firstPart = false;
            }
            else
            {
                snprintf(txBuffer, sizeof(txBuffer) - 1, "z%.3fw%.3f\r\n", quaternion[3][0], quaternion[4][0]);
                firstCharacteristic.writeValue(txBuffer);
                is_firstPart = true;
            }
        }
        digitalWrite(LED_BUILTIN, LOW);
    }
}
