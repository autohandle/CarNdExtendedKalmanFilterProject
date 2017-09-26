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

  // initializing matrices
  R_laser_ = MatrixXd(2, 2); // px, py
  R_radar_ = MatrixXd(3, 3); // rho, phi, rho dot
    
  H_laser_ = MatrixXd(2, 4);
  H_laser_ <<   1, 0, 0, 0,  // measure vector is: px py
                0, 1, 0, 0;
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ <<   0.0225, 0,
                0,      0.0225;

  //measurement covariance matrix - radar
  R_radar_ <<   0.09,   0,      0,
                0,      0.0009, 0,
                0,      0,      0.09;

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
    double dt = 0.;
    
    cout << "ProcessMeasurement-measurement_pack.sensor_type_:" << measurement_pack.sensor_type_ << ", dt:" << dt << "," << endl << "ekf_.x_:" << Tools::toString(ekf_.x_) << endl;
    cout << "ProcessMeasurement-is radar:" << (measurement_pack.sensor_type_==MeasurementPackage::RADAR) << endl;
    
    
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
        previous_timestamp_= measurement_pack.timestamp_; // initial dt == 0.
        // first measurement
        
        cout << "ProcessMeasurement-Initialization-EKF: " << endl;
        ekf_.x_ = VectorXd(4);
        ekf_.x_ << 1, 1, 1, 1;
        //ekf_.x_ << 0, 0, 0, 0;
        
        // KalmanFilter ekf_;
        // ekf_.x_ = measurement_pack.raw_measurements_;// this my not zero vx & vy
        
        //state covariance matrix P
        ekf_.P_ = MatrixXd(4, 4);
        ekf_.P_ <<  1, 0, 0,    0,
                    0, 1, 0,    0,
                    0, 0, 1000, 0,
                    0, 0, 0,    1000;
        // done initializing, no need to predict or update
        ekf_.F_ = Tools::makeF(dt, ekf_.x_);
        ekf_.Q_ = Tools::makeQ(dt, ekf_.x_, processNoise);
        if (Tools::TESTING) std::cout << "ProcessMeasurement-init-F:" <<  Tools::toString(ekf_.F_) << std::endl
            << "Q:" <<  Tools::toString(ekf_.Q_) << std::endl;
        
        //measurement matrix
        ekf_.H_ = MatrixXd(2, 4);
        ekf_.H_ <<  1, 0, 0, 0,
                    0, 1, 0, 0;
        
        //measurement covariance
        ekf_.R_ = MatrixXd(2, 2);
        ekf_.R_ <<  0.0225, 0,
                    0,      0.0225;
        
        ekf_.x_=Tools::updateX(ekf_.x_, measurement_pack.raw_measurements_);
        
        //if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
        //    FusionEKF::updateFromRadar(dt, measurement_pack.raw_measurements_);
        //} else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
        //    FusionEKF::updateFromLaser(dt, measurement_pack.raw_measurements_);
        //}
        
        is_initialized_ = true;
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
    
    ekf_.F_=Tools::updateF(dt, ekf_.F_); // linear F for both LASER & RADAR
    ekf_.Predict(); // same linear predict for both LASER & RADAR
    std::cout << "Tools::ProcessMeasurement-after predict-P:" << Tools::toString(ekf_.P_) << std::endl
                << "x:" <<Tools::toString(ekf_.x_) << std::endl;
    
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
        FusionEKF::updateFromRadar(dt, measurement_pack.raw_measurements_);
    } else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
        // Laser updates
        FusionEKF::updateFromLaser(dt, measurement_pack.raw_measurements_);
    }
    
    // print the output
    cout << "x_ = " << ekf_.x_ << endl;
    cout << "P_ = " << ekf_.P_ << endl;
}

void FusionEKF::updateFromLaser(double dt, const VectorXd &rawMeasurements) {
    assert(Tools::isNotZero(dt));
    ekf_.Q_=Tools::updateQ(dt, processNoise, ekf_.Q_);
    ekf_.H_=H_laser_;
    ekf_.R_=R_laser_;
    if (Tools::TESTING) std::cout << "updateFromLaser-F:" <<  Tools::toString(ekf_.F_) << std::endl
                                  << ", Q:" <<  Tools::toString(ekf_.Q_) << std::endl
                                  << ", H:" <<  Tools::toString(ekf_.H_) << std::endl
                                  << ", R:" <<  Tools::toString(ekf_.R_) << std::endl;
    
    ekf_.Update(rawMeasurements, false /* isRadar==false*/);// 2x for laser
    /**
     Initialize state.
     */
}

void FusionEKF::updateFromRadar(double dt, const VectorXd &rawMeasurements) {
    assert(Tools::isNotZero(dt));
    /**
     Convert radar from polar to cartesian coordinates and initialize state.
     */
    ekf_.Q_=Tools::updateQ(dt, processNoise, ekf_.Q_);
    ekf_.H_=Tools::CalculateHj(ekf_.x_);
    ekf_.R_=R_radar_;
    if (Tools::TESTING) std::cout << "updateFromRadar-F:" <<  Tools::toString(ekf_.F_) << std::endl
                                  << ", Q:" <<  Tools::toString(ekf_.Q_) << std::endl
                                  << ", H:" <<  Tools::toString(ekf_.H_) << std::endl
                                  << ", R:" <<  Tools::toString(ekf_.R_) << std::endl;    
    ekf_.Update(rawMeasurements, true /* isRadar==false*/);
    /**
     Initialize state.
     */
}
