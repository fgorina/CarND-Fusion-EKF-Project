#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;


/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
    is_initialized_ = false;
    
    previous_timestamp_ = 0;
    
    // initializing matrices
    R_laser_ = MatrixXd(2, 2);
    R_radar_ = MatrixXd(3, 3);
    H_laser_ = MatrixXd(2, 4);
    Hj_ = MatrixXd(3, 4);
    
    
    //measurement covariance matrix - laser
    R_laser_ << 0.0225, 0,
    0, 0.0225;
    
    //measurement covariance matrix - radar
    R_radar_ << 0.09, 0, 0,
    0, 0.0009, 0,
    0, 0, 0.09;
    
    /**
     TODO:
     * Finish initializing the FusionEKF.
     * Set the process and measurement noises
     */
    
    H_laser_ <<  1, 0, 0, 0,
        0, 1, 0, 0;
    
    Hj_ << 0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0;
    
    
    ekf_ = KalmanFilter();
    
    
}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
    
    
    /*****************************************************************************
     *  Initialization
     ****************************************************************************/
    if (!is_initialized_) {
        /**
         TODO:
         * Initialize the state ekf_.x_ with the first measurement.
         * Create the covariance matrix.
         * Remember: you'll need to convert radar from polar to cartesian coordinates.
         */
        // first measurement
        cout << "EKF: " << endl;
        ekf_.x_ = VectorXd(4);
        ekf_.x_ << 1, 1, 0, 0;
        
        ekf_.F_ = MatrixXd(4,4);
        ekf_.F_ << 1, 0, 1, 0,
            0, 1, 0, 1,
            0, 0, 1, 0,
            0, 0, 0, 1;

        
        previous_timestamp_ = measurement_pack.timestamp_;
        
        if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
            /**
             Convert radar from polar to cartesian coordinates and initialize state.
             */
            
            double r = measurement_pack.raw_measurements_[0];
            double theta = measurement_pack.raw_measurements_[1];
            
            ekf_.x_[0]  = r * cos(theta);
            ekf_.x_[1]  = r * sin(theta);
            ekf_.P_ = MatrixXd(4,4);
            
            ekf_.P_ << R_radar_(0,0)*cos(theta)+r*sin(theta)*R_radar_(1,1), 0, 0, 0,
            0, R_radar_(0,0)*sin(theta)+r*cos(theta)*R_radar_(1,1), 0, 0,
            0, 0, 10000, 0,
            0, 0, 0, 10000;

        }
        else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
            /**
             Initialize state.
             */
            
            ekf_.x_[0] = measurement_pack.raw_measurements_[0];
            ekf_.x_[1] = measurement_pack.raw_measurements_[1];
            ekf_.P_ = MatrixXd(4,4);
            
            ekf_.P_ << R_laser_(0,0), 0, 0, 0,
            0, R_laser_(1,1), 0, 0,
            0, 0, 10000, 0,
            0, 0, 0, 10000;

        }
        
        
        // done initializing, no need to predict or update
        is_initialized_ = true;
        return;
    }
    
    /*****************************************************************************
     *  Prediction
     ****************************************************************************/
    
    double delta = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
    
    double sigma2x = 9;     // Logically should came from measurement?
    double sigma2y = 9;
   
    ekf_.F_(0, 2) = delta;
    ekf_.F_(1, 3) = delta;
    
    ekf_.Q_ = MatrixXd(4, 4);
    
    ekf_.Q_ << pow(delta,4)*sigma2x/4, 0, pow(delta,3)*sigma2x/2, 0,
        0, pow(delta,4)*sigma2y/4, 0, pow(delta,3)*sigma2y/2,
        pow(delta,3)*sigma2x/2, 0, pow(delta,2)*sigma2x, 0,
        0, pow(delta,3)*sigma2y/2, 0, pow(delta,2)*sigma2y;
     
    ekf_.Predict();
    
    /*****************************************************************************
     *  Update
     ****************************************************************************/
    
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
        if (true){
            ekf_.R_ = R_radar_;
            ekf_.UpdateEKF(measurement_pack.raw_measurements_);
            Hj_ = ekf_.H_;
            previous_timestamp_ = measurement_pack.timestamp_;
        }
        
    } else {
        if (true){
            ekf_.H_ = H_laser_;
            ekf_.R_ = R_laser_;
            ekf_.Update(measurement_pack.raw_measurements_);
            H_laser_ = ekf_.H_;
            // Laser updates
            previous_timestamp_ = measurement_pack.timestamp_;
        }
    }
    
   
    
    // print the output
    cout << "x_ = " << ekf_.x_ << endl;
    //cout << "P_ = " << ekf_.P_ << endl;
}
