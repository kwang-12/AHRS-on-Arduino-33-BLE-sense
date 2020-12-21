# Wireless AHRS on Arduino 33 BLE sense

## Team Members
* Kezhuo Wang

## Links
### Project Repositories
[Github](https://github.com/kwang-12/AHRS-on-Arduino-33-BLE-sense)
### Final Presentation/DEMO Video
[Youtube](https://youtu.be/5BwvYRFV4Gc)

------------
## Project Introduction
#### Motivation
This project aims to implement an extended kalman filter to estimate the orientation of the Arduino chip (Arduino 33 BLE sense) using the on-board 9DOF IMU (LSM9DS1).

#### Goal
The goal is to enable wireless data transmission between the Arduino chip (read IMU data and estimate orientation) and the computer (data visualize and logging) via BLE.

#### Deliverable
1. Arduino code for estimating orientation and wireless data transmission via onboard BLE functionality.
2. App for data logging and visualizing on PC

#### Hardware used in this project
* Arduino 33 BLE sense
* 9v battery and connector
* PC

![device](/pics/device.png)

## Part 1: Orientation estimation with EKF
### Filter description
The complete EKF formulation is as below. 

![equation](/pics/1.png)

The filter has 4 states, 3 inputs and 6 measurments.
The 4 states are the 4 values from the quaternion that represented the system's orientation(q0, q1, q2, q3).
The 3 inputs are the angular rate measured by the gyroscope(p, q, r).
The 6 measurements are the 3 linear accelerations from the accelerometer(ax, ay, az) and 3 magnetic field strength from the magnetometer(Mx, My, Mz).
Magnetometer hard iron bias are: Mxb, Myb, Mzb. Sampling time is delta_t.

The prediction step in detail is:

![equation](/pics/2.png)

The measurement residuals are obtained:

![equation](/pics/3.png)

Jacobian H

![equation](/pics/4.png)

Jacobian F

![equation](/pics/5.png)

### Results
Logged orientation estimate vs estimation from smartphone
Estimation from smartphone is used as ground truth from evaluating implemented filter's performance. 
As it's shown in the figure below, the implemented filter is capable of producing reasonable results that roughly follows the ground truth.

![result](/pics/result.png)

## Part 2: Wireless data transmission via BLE
### Implementation
The 4 numbers from the quaternion are transmitted as 3 decimal place fractional values. Each number needs at least 5 or 6 byte depending on the sign. 
And since the maximum packet size for BLE is 20 bytes, the transmission package is divided into two parts: q0 and q1 are sent in one package, q2 and q3 are sent in another package. The receiving device will wait for all 4 numbers to be received before updating its stored quaternion and act accordingly.


## Part 3: Data logging and visualization on PC
### Description
A custom python app is written to handle logging and visualizing the data.

It features a wireframe render to visualize the current orientation of the device. The corresponding quaternion and euler angles are also displayed below.
The app allows logging quaternion data time history into csv file for post-processing.

![picture](/pics/python_app.png)


## Summary
Arduino 33 BLE sense comes with a pack full of sensors. The integration of bluetooth module and IMU made this project possible. Wireless data transmission minized the need for wiring and enabled the device to be deployed remotely. The device can be used for spatial localization, measure levelness and etc.
For future works, alternative communication methods should be explored to improve the wireless data transfer performance.

### Strength
* Can be deployed remotely.
* Accuracy: comparable performance against orientation filters on smartphones.
### Weakness
* Constant update requirement makes saving energy difficult.
* Requires a computer for data logging and visualization.
* Lack of interface makes it hard to interpret the state of device.
