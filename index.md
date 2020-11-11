## Welcome to GitHub Pages

This project aims to implement a kalman filter to estimate the orientation of the Arduino chip (Arduino 33 BLE sense) using the on-board 9DOF IMU (LSM9DS1).

Using euler angles, the orientation of the Arduino chip is defined by three independent angles: yaw, pitch, and roll relative to the world frame (inertial frame). 

### The reasons for using Arduino 33 BLE sense are: 

1. It includes a 9DOF (accelerometer+gyroscope+magnetometer) IMU which is essential to the filter estimation. 
The pitch and roll angle can be easily estimated with a complementary filter and a 6 DOF IMU (accelerometer+gyroscope). However, since a 6 DOF IMU doesn’t provide the necessary information for estimating yaw position. Even the slightest noise and biases in sensor reading can cause estimated yaw position to drift.
The additional magnetometer onboard the Arduino 33 BLE sense chip provides the critical information (earth magnetic field as a reference) for the filter to de-drift yaw position and output accurate orientation estimation.
2. The additional bluetooth functionality means that the Arduino chip can be deployed remotely and output orientation estimations via bluetooth protocol.
3. Minimal wiring.

