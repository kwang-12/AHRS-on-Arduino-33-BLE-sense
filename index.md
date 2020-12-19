# Wireless AHRS on Arduino 33 BLE sense

## Team Members
* Kezhuo Wang

## Links
### Project Repositories
[Github](https://github.com/kwang-12/AHRS-on-Arduino-33-BLE-sense)
### Final Presentation/DEMO Video
Coming soon

------------
## Project Introduction
#### Motivation
This project aims to implement an extended kalman filter to estimate the orientation of the Arduino chip (Arduino 33 BLE sense) using the on-board 9DOF IMU (LSM9DS1).

#### Goal
The goal is to enable wireless data transmission between the Arduino chip (read IMU data and estimate orientation) and the computer (data visualize and logging) via BLE.

#### Deliverable
1. Arduino code for estimating orientation
2. Executable for data logging and visualizing on PC

#### Hardware used in this project
* Arduino 33 BLE sense
* 9v battery and connector
* PC

## Part 1: Orientation estimation with EKF
### Filter/code description
Coming soon
### Results
Coming soon

## Part 2: Wireless data transmission via BLE
### Communication/code description
Coming soon
### Results
Coming soon

## Part 3: Data logging and visualization on PC
### Python code description
Coming soon
### Results
Coming soon

## Summary
Arduino 33 BLE sense comes with a pack full of sensors. The integration of bluetooth module and IMU made this project possible. Wireless data transmission minized the need for wiring and enabled the device to be deployed remotely. The device can be used for spatial localization, measure levelness and etc.

### Strength
* Can be deployed remotely.
* Accuracy: comparable performance against orientation filters on smartphones.
* Flexibility: Arduino BLE 33 can do more since the IMU code excutes efficiencly at ~600microsecond per calculation cycle in a non-blocking manner.
### Weakness
* Constant update requirement makes saving energy difficult.
* Requires a computer for data inspection.
