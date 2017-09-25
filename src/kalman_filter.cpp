#include "kalman_filter.h"
#include "tools.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
  TODO:
    * predict the state
  */
    
    x_ = F_ * x_;
    P_ = (F_ * P_ * F_.transpose()) + Q_;
    
}


void KalmanFilter::Update(const VectorXd &z) {
   
    VectorXd y  = z - H_ * x_; // H_ is fixed in this case. No need to recompute it every timestep
    UpdateKalman(y);
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
    
    Tools  t = Tools();
    H_ = t.CalculateJacobian(x_);   // In Extended Kalman H is the Jacobian
    
    // Compute the converted state. Function is applied directly
    //
    
    VectorXd xc(3);
        xc << sqrt(x_[0]*x_[0]+x_[1]*x_[1]),
        atan2(x_[1],x_[0]),
        (x_[0]*x_[2]+x_[1]*x_[3])/sqrt(x_[0]*x_[0]+x_[1]*x_[1]);
    
    VectorXd y = z - xc;  // Compute Error
    
    // Normalize angle error. That is important if not error may be offseted by 2PI
    
    while (y[1] > M_PI){
        y[1] -= 2*M_PI;
    }
 
    while (y[1] < -M_PI){
        y[1] += 2*M_PI;
    }
    UpdateKalman(y);
 }

void KalmanFilter::UpdateKalman(VectorXd &error){
    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    
    MatrixXd Ht = H_.transpose();
    MatrixXd S = (H_ * P_ * Ht) + R_;
    MatrixXd K = P_ * Ht * S.inverse();
    
    x_ = x_ + (K * error);
    P_ = (I - K * H_) * P_;

}
