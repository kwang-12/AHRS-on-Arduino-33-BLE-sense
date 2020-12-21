#pragma once
#include "matrix.h"

class EKF
{
public:
    EKF(Matrix &X, Matrix &P, Matrix &Q, Matrix &R, double sample_dt);
    void reset(Matrix &X, Matrix &P, Matrix &Q, Matrix &R);
    bool update(double p, double q, double r, Matrix &Y);
    Matrix GetState() { return Xm; }
private:
    double sample_dt = 0;
    Matrix Fm = Matrix(4, 4);
    Matrix Xm = Matrix(4, 1);
    Matrix P = Matrix(4, 4);
    Matrix H = Matrix(6, 4);
    Matrix Hx = Matrix(6, 1);
    Matrix Q = Matrix(4, 4);
    Matrix R = Matrix(6, 6);
    Matrix S = Matrix(6, 6);
    Matrix Gain = Matrix(4, 6);
};