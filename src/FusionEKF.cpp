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
  R_laser_ = MatrixXd(2, 2); // px, py
  R_radar_ = MatrixXd(3, 3); // rho, phi, rho dot
  H_laser_ = MatrixXd(2, 4);
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


}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {

  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
  previous_timestamp_ = measurement_pack.timestamp_;

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
      
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;
    
    // KalmanFilter ekf_;
    // ekf_.x_ = measurement_pack.raw_measurements_;// this my not zero vx & vy
      
    //state covariance matrix P
    ekf_.P_ = MatrixXd(4, 4);
    ekf_.P_ <<  1, 0, 0,    0,
                0, 1, 0,    0,
                0, 0, 1000, 0,
                0, 0, 0,    1000;
    // done initializing, no need to predict or update
    ekf_.F_ = Tools::makeF(0., ekf_.x_);
    ekf_.Q_ = Tools::makeQ(0., ekf_.x_, processNoise);
      
    //measurement matrix
    ekf_.H_ = MatrixXd(2, 4);
    ekf_.H_ <<  1, 0, 0, 0,
                0, 1, 0, 0;
      
    //measurement covariance
    ekf_.R_ = MatrixXd(2, 2);
    ekf_.R_ <<  0.0225, 0,
                0,      0.0225;
      
    ekf_.x_=Tools::updateX(ekf_.x_, measurement_pack.raw_measurements_);
    cout << "measurement_pack.sensor_type_:" << measurement_pack.sensor_type_ << ", dt:" << dt << "," << endl << "ekf_.x_:" << Tools::toString(ekf_.x_) << endl;
    cout << "is radar:" << (measurement_pack.sensor_type_==MeasurementPackage::RADAR) << endl;
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
        FusionEKF::updateFromRadar(dt, measurement_pack);
    } else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
        FusionEKF::updateFromLaser(dt, measurement_pack);
    }
      
      is_initialized_ = true;
      return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */

  ekf_.Predict();

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
    FusionEKF::updateFromRadar(dt, measurement_pack);
  } else {
    // Laser updates
    FusionEKF::updateFromLaser(dt, measurement_pack);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}

void FusionEKF::updateFromLaser(double dt, const MeasurementPackage &measurement_pack) {
    ekf_.F_=Tools::updateF(dt, &ekf_.F_);
    ekf_.Q_=Tools::updateQ(dt, measurement_pack.raw_measurements_, &ekf_.Q_);
    if (Tools::TESTING) std::cout << "updateFromLaser-F" <<  Tools::toString(ekf_.F_) << "Q" <<  Tools::toString(ekf_.Q_) << std::endl;
    
    ekf_.Predict();
    ekf_.Update(measurement_pack.raw_measurements_);
    //ekf_.Q_(100,100)=1;
    /**
     Initialize state.
     */
}

void FusionEKF::updateFromRadar(double dt, const MeasurementPackage &measurement_pack) {
    /**
     Convert radar from polar to cartesian coordinates and initialize state.
     */
    ekf_.F_=Tools::updateF(dt, &ekf_.F_);
    ekf_.Q_=Tools::updateQ(dt, measurement_pack.raw_measurements_, &ekf_.Q_);
    if (Tools::TESTING) std::cout << "updateFromLaser-F" <<  Tools::toString(ekf_.F_) << "Q" <<  Tools::toString(ekf_.Q_) << std::endl;
    
    ekf_.Predict();
    ekf_.Update(measurement_pack.raw_measurements_);
    //ekf_.Q_(100,100)=1;
    /**
     Initialize state.
     */
}
