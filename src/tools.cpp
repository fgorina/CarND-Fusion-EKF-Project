#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
    
    VectorXd rmse(4);
    rmse << 0,0,0,0;
    
    if (estimations.size() == 0 || estimations.size() != ground_truth.size()){
        return rmse;
    }
    
    int n = estimations.size();
    
    //std::cout << "Size : " << estimations.size() << std::endl;
    
    //accumulate squared residuals
    for(int i=0; i < estimations.size(); ++i){
        VectorXd e = estimations[i];
        VectorXd g = ground_truth[i];
        
        VectorXd diff = e - g;
        VectorXd diffSqrd = diff.array() * diff.array();
        rmse += diffSqrd;
        
    }
    
    //calculate the mean
    // ... your code here
    rmse = rmse.array()/n;
    rmse = rmse.array().sqrt();
    //calculate the squared root
    // ... your code here
    
    //return the result
    return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
  TODO:
    * Calculate a Jacobian here.
  */
     MatrixXd Hj(3,4);
    //recover state parameters
    double px = x_state(0);
    double py = x_state(1);
    double vx = x_state(2);
    double vy = x_state(3);
    
    //pre-compute a set of terms to avoid repeated calculation
    double c1 = px*px+py*py;
    double c2 = sqrt(c1);
    double c3 = (c1*c2);
    
    //check division by zero
    if(fabs(c1) < 0.0001){
        cout << "CalculateJacobian () - Error - Division by Zero" << endl;
        return Hj;
    }
    
    //compute the Jacobian matrix
    Hj << (px/c2), (py/c2), 0, 0,
    -(py/c1), (px/c1), 0, 0,
    py*(vx*py - vy*px)/c3, px*(px*vy - py*vx)/c3, px/c2, py/c2;
    
    return Hj;
    
}
