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
  R_radar_ = MatrixXd(3, 3); // rho, phi, rhoDot
    
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
  //R_radar_ <<   0.09,   0,      0,
  //              0,      0.0006, 0,
  //              0,      0,      0.09;


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
    
    cout << "ProcessMeasurement-measurement_pack.sensor_type_:" << measurement_pack.sensor_type_ << "," << endl << "ekf_.x_:" << Tools::toString(ekf_.x()) << endl;
    cout << "ProcessMeasurement-is radar:" << (measurement_pack.sensor_type_==MeasurementPackage::RADAR) << endl;
    
    
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
        
        cout << "ProcessMeasurement-Initialization-EKF: " << endl;
        ekf_.x() = VectorXd(4);
        //ekf_.x_ << 1, 1, 1, 1;
        ekf_.x() << 0, 0, 0, 0;
        
        // KalmanFilter ekf_;
        // ekf_.x_ = measurement_pack.raw_measurements_;// this my not zero vx & vy
        
        //state covariance matrix P
        ekf_.P() = MatrixXd(4, 4);
        ekf_.P() <<  1, 0, 0,    0,
                    0, 1, 0,    0,
                    0, 0, 1000, 0,
                    0, 0, 0,    1000;
        // done initializing, no need to predict or update
        ekf_.F() = Tools::makeF(dt, ekf_.x());
        ekf_.Q() = Tools::makeQ(dt, ekf_.x(), processNoise);
        if (Tools::TESTING) std::cout << "ProcessMeasurement-init-F:" <<  Tools::toString(ekf_.F()) << std::endl
            << "Q:" <<  Tools::toString(ekf_.Q()) << std::endl;
        
        //measurement matrix
        ekf_.H() = MatrixXd(2, 4);
        ekf_.H() <<  1, 0, 0, 0,
                    0, 1, 0, 0;
        
        //measurement covariance
        ekf_.R() = MatrixXd(2, 2);
        ekf_.R() <<  0.0225, 0,
                    0,      0.0225;
        
        
        if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
            //FusionEKF::updateFromRadar(dt, measurement_pack.raw_measurements_);
            ekf_.x()=Tools::convertRadarToStateVector(measurement_pack.raw_measurements_);
        } else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
            ekf_.x()=Tools::updateX(ekf_.x(), measurement_pack.raw_measurements_);
            //FusionEKF::updateFromLaser(dt, measurement_pack.raw_measurements_);
        }
        
        is_initialized_ = true;
        if (true | Tools::TESTING) {
            std::cout << "ProcessMeasurement-init:" <<  ekf_.toString() << std::endl;
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
    
    //ekf_.Predict(dt, processNoise); // same linear predict for both LASER & RADAR
    
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
        FusionEKF::updateFromRadar(dt, measurement_pack.raw_measurements_);
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
        FusionEKF::updateFromLaser(dt, measurement_pack.raw_measurements_);
        if (Tools::TESTING) {
            std::cout << "ProcessMeasurement-after laser update:" <<  ekf_.toString() << std::endl;
        }
    }
    
    // print the output
    cout << "x_ = " << ekf_.x() << endl;
    cout << "P_ = " << ekf_.P() << endl;
}

void FusionEKF::updateFromLaser(double dt, const VectorXd &rawMeasurements) {
    assert(Tools::isNotZero(dt));
    //ekf_.Q_=Tools::updateQ(dt, processNoise, ekf_.Q_);
    ekf_.H()=H_laser_;
    ekf_.R()=R_laser_;
    if (Tools::TESTING) {
        if (Tools::TESTING) std::cout << "updateFromLaser-F:" <<  Tools::toString(ekf_.F()) << std::endl
                                  << ", Q:" <<  Tools::toString(ekf_.Q()) << std::endl
                                  << ", H:" <<  Tools::toString(ekf_.H()) << std::endl
                                  << ", R:" <<  Tools::toString(ekf_.R()) << std::endl;
    }
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
    //ekf_.Q_=Tools::updateQ(dt, processNoise, ekf_.Q_);
    ekf_.H()=Tools::CalculateHj(ekf_.x());
    ekf_.R()=R_radar_;
    if (Tools::TESTING) std::cout << "updateFromRadar-before-F:" <<  Tools::toString(ekf_.F()) << std::endl
                                  << ", P:" <<  Tools::toString(ekf_.P()) << std::endl
                                  << ", Q:" <<  Tools::toString(ekf_.Q()) << std::endl
                                  << ", H:" <<  Tools::toString(ekf_.H()) << std::endl
                                  << ", R:" <<  Tools::toString(ekf_.R()) << std::endl;
    ekf_.Update(rawMeasurements, true /* isRadar==false*/);
    if (Tools::TESTING) std::cout << "updateFromRadar-after-P:" <<  Tools::toString(ekf_.P()) << std::endl
        << ", x:" <<  Tools::toString(ekf_.x()) << std::endl;
    /**
     Initialize state.
     */
}
