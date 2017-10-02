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
    
    previous_timestamp_ = 0.;
    
    /**
     TODO:
     * Finish initializing the FusionEKF.
     * Set the process and measurement noises
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
     */
    processNoise = VectorXd(2); // ax, ay
    processNoise << 9. , 9.;
    //processNoise << 5. , 5.;
    
}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
    
    /*****************************************************************************
     *  Initialization
     ****************************************************************************/
    double dt = 0.;
    if (!is_initialized_) {
        /**
         TODO:
         * Initialize the state ekf_.x_ with the first measurement.
         * Create the covariance matrix.
         * Remember: you'll need to convert radar from polar to cartesian coordinates.
         */
        previous_timestamp_= measurement_pack.timestamp_; // initial dt == 0.
        // first measurement
        
        ekf_.x() = VectorXd(4);
        ekf_.x() << 1, 1, 1, 1;
        //ekf_.x() << 0, 0, 0, 0;
        
        // KalmanFilter ekf_;
        // ekf_.x_ = measurement_pack.raw_measurements_;// this my not zero vx & vy
        
        //state covariance matrix P
        ekf_.P() = MatrixXd(4, 4);
        ekf_.P() <<     1, 0, 0,    0,
                        0, 1, 0,    0,
                        0, 0, 1000, 0,
                        0, 0, 0,    1000;
        // done initializing, no need to predict or update
        
        ekf_.initF(dt);
        ekf_.initQ(dt, processNoise);
        
        if (Tools::TESTING) {
            std::cout << "ProcessMeasurement-init-F&Q:" <<  ekf_.toString() << std::endl;
        }
        
        if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
            radarFilter.updateX(measurement_pack.raw_measurements_);
        } else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
            laserFilter.updateX(measurement_pack.raw_measurements_);
        }
        
        is_initialized_ = true;
        if (Tools::TESTING) {
            std::cout << "ProcessMeasurement-init-X:" <<  ekf_.toString() << std::endl;
        }
        return;
    }
    
    /*****************************************************************************
     *  Prediction
     ****************************************************************************/
    
    dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
    previous_timestamp_ = measurement_pack.timestamp_;
    
    /**
     TODO:
     * Update the state transition matrix F according to the new elapsed time.
     - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
     */
    
    /*****************************************************************************
     *  Update
     ****************************************************************************/
    
    /**
     TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
     */
    
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
        // Radar updates
        if (Tools::TESTING) {
            std::cout << "ProcessMeasurement-before radar predict:" <<  ekf_.toString() << std::endl;
        }
        radarFilter.Predict(dt, processNoise);
        if (Tools::TESTING) {
            std::cout << "ProcessMeasurement-after radar predict:" <<  ekf_.toString() << std::endl;
        }
        radarFilter.Update(measurement_pack.raw_measurements_, true);
        if (Tools::TESTING) {
            std::cout << "ProcessMeasurement-after radar update:" <<  ekf_.toString() << std::endl;
        }
    } else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
        // Laser updates
        if (Tools::TESTING) {
            std::cout << "ProcessMeasurement-before laser predict:" <<  ekf_.toString() << std::endl;
        }
        laserFilter.Predict(dt, processNoise);
        if (Tools::TESTING) {
            std::cout << "ProcessMeasurement-after laser predict:" <<  ekf_.toString() << std::endl;
        }
        laserFilter.Update(measurement_pack.raw_measurements_);
        if (Tools::TESTING) {
            std::cout << "ProcessMeasurement-after laser update:" <<  ekf_.toString() << std::endl;
        }
    }
    
    // print the output
    if (Tools::TESTING) {
        cout << "x_ = " << ekf_.x() << endl;
        cout << "P_ = " << ekf_.P() << endl;
    }
}

