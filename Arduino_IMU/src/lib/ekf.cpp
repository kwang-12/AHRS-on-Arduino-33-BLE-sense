#include "ekf.h"
EKF::EKF(Matrix &X, Matrix &P, Matrix &Xm, Matrix &R, double sample_dt)
{
    this->Xm = X.replicate();
    this->P = P.replicate();
    this->Xm = Xm.replicate();
    this->R = R.replicate();
    this->sample_dt = sample_dt;
}

void EKF::reset(Matrix &X, Matrix &P, Matrix &Q, Matrix &R)
{
    this->Xm = X.replicate();
    this->P = P.replicate();
    this->Q = Q.replicate();
    this->R = R.replicate();
}

bool EKF::update(double p, double q, double r, Matrix &Y)
{
    // Update F jacobian
    Fm[0][0] = 1.000;
    Fm[0][1] = -0.5 * p * sample_dt;
    Fm[0][2] = -0.5 * q * sample_dt;
    Fm[0][3] = -0.5 * r * sample_dt;

    Fm[1][0] = 0.5 * p * sample_dt;
    Fm[1][1] = 1.000;
    Fm[1][2] = 0.5 * r * sample_dt;
    Fm[1][3] = -0.5 * q * sample_dt;

    Fm[2][0] = 0.5 * q * sample_dt;
    Fm[2][1] = -0.5 * r * sample_dt;
    Fm[2][2] = 1.000;
    Fm[2][3] = 0.5 * p * sample_dt;

    Fm[3][0] = 0.5 * r * sample_dt;
    Fm[3][1] = 0.5 * q * sample_dt;
    Fm[3][2] = -0.5 * p * sample_dt;
    Fm[3][3] = 1.000;

    // Update priori
    Xm[0][0] += (0.5 * (-p * Xm[1][0] - q * Xm[2][0] - r * Xm[3][0])) * sample_dt;
    Xm[1][0] += (0.5 * (p * Xm[0][0] + r * Xm[2][0] - q * Xm[3][0])) * sample_dt;
    Xm[2][0] += (0.5 * (q * Xm[0][0] - r * Xm[1][0] + p * Xm[3][0])) * sample_dt;
    Xm[3][0] += (0.5 * (r * Xm[0][0] + q * Xm[1][0] - p * Xm[2][0])) * sample_dt;

    // Update priori covariance
    P = Fm * P * (Fm.transpose()) + Q;

    // Update H jacobian
    H[0][0] = -2 * Xm[2][0];
    H[0][1] = 2 * Xm[3][0];
    H[0][2] = -2 * Xm[0][0];
    H[0][3] = 2 * Xm[1][0];

    H[1][0] = 2 * Xm[1][0];
    H[1][1] = 2 * Xm[0][0];
    H[1][2] = 2 * Xm[3][0];
    H[1][3] = 2 * Xm[2][0];

    H[2][0] = 2 * Xm[0][0];
    H[2][1] = -2 * Xm[1][0];
    H[2][2] = -2 * Xm[2][0];
    H[2][3] = 2 * Xm[3][0];

    H[3][0] = 2 * Xm[0][0];
    H[3][1] = 2 * Xm[1][0];
    H[3][2] = -2 * Xm[2][0];
    H[3][3] = -2 * Xm[3][0];

    H[4][0] = -2 * Xm[3][0];
    H[4][1] = 2 * Xm[2][0];
    H[4][2] = 2 * Xm[1][0];
    H[4][3] = -2 * Xm[0][0];

    H[5][0] = 2 * Xm[2][0];
    H[5][1] = 2 * Xm[3][0];
    H[5][2] = 2 * Xm[0][0];
    H[5][3] = 2 * Xm[1][0];

    // Update gain
    S = (H * P * (H.transpose())) + R;
    Gain = P * (H.transpose()) * (S.inverse());

    // Map state to measurement space
    double xest00_2 = Xm[0][0] * Xm[0][0];
    double xest10_2 = Xm[1][0] * Xm[1][0];
    double xest20_2 = Xm[2][0] * Xm[2][0];
    double xest30_2 = Xm[3][0] * Xm[3][0];
    Hx[0][0] = 2 * Xm[1][0] * Xm[3][0] - 2 * Xm[0][0] * Xm[2][0];
    Hx[1][0] = 2 * Xm[2][0] * Xm[3][0] + 2 * Xm[0][0] * Xm[1][0];
    Hx[2][0] = xest00_2 - xest10_2 - xest20_2 + xest30_2;
    Hx[3][0] = xest00_2 + xest10_2 - xest20_2 - xest30_2;
    Hx[4][0] = 2 * (Xm[1][0] * Xm[2][0] - Xm[0][0] * Xm[3][0]);
    Hx[5][0] = 2 * (Xm[1][0] * Xm[3][0] + Xm[0][0] * Xm[2][0]);

    // Update posteriori
    Xm = Xm + (Gain * (Y - Hx));

    // Normalize
    if (!Xm.normalizeVector())
    {
        Xm[0][0] = double(1.0);
        Xm[1][0] = double(0.0);
        Xm[2][0] = double(0.0);
        Xm[3][0] = double(0.0);
    }
    Matrix I = Matrix(4, 4);
    I.setValueDiagonal(double(1.0));
    P = (I - (Gain * H)) * P;
    return true;
}
