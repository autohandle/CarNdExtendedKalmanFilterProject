#ifndef FusionEKF_H_
#define FusionEKF_H_

#include "measurement_package.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>
#include "kalman_filter.h"
#include "tools.h"

class FusionEKF {
public:
    /**
     * Constructor.
     */
    FusionEKF();
    
    /**
     * Destructor.
     */
    virtual ~FusionEKF();
    
    /**
     * Run the whole flow of the Kalman Filter from here.
     */
    void ProcessMeasurement(const MeasurementPackage &measurement_pack);
    void updateFromLaser(double dt, const VectorXd &x);
    void updateFromRadar(double dt, const VectorXd &x);
    
    /**
     * Kalman Filter update and prediction math lives in here.
     */
    KalmanFilterArrays anotherKalmanFilterArrays;
    KalmanFilterArrays kalmanFilterArrays;
    KalmanFilter ekf_ = KalmanFilter::KalmanFilter(kalmanFilterArrays);
    LaserFilter laserFilter = LaserFilter(kalmanFilterArrays);;
    RadarFilter radarFilter = RadarFilter(kalmanFilterArrays);;
    
    
private:
    // check whether the tracking toolbox was initialized or not (first measurement)
    bool is_initialized_;
    
    // previous timestamp
    long long previous_timestamp_;
    
    // tool object used to compute Jacobian and RMSE
    Tools tools;
    
    Eigen::MatrixXd R_laser_;
    Eigen::MatrixXd R_radar_;
    Eigen::MatrixXd H_laser_;
    Eigen::MatrixXd Hj_;
    
    Eigen::VectorXd processNoise;
    Eigen::VectorXd measurementNoiseRadar;
    Eigen::VectorXd measurementNoiseLidar;
    
};

#endif /* FusionEKF_H_ */
