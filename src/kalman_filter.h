#ifndef KALMAN_FILTER_H_
#define KALMAN_FILTER_H_
#include "Eigen/Dense"
#include <iostream>
#include "tools.h"

class KalmanFilterArrays {
public:
    
    // state vector
    Eigen::VectorXd x_;
    Eigen::VectorXd& x() {
        return x_;
    }
    
    // state covariance matrix
    Eigen::MatrixXd P_;
    Eigen::MatrixXd& P() {
        return P_;
    }
    
    // state transition matrix
    Eigen::MatrixXd F_;
    Eigen::MatrixXd& F() {
        return F_;
    }
    
    // process covariance matrix
    Eigen::MatrixXd Q_;
    Eigen::MatrixXd& Q() {
        return Q_;
    }
    
    // measurement matrix
    Eigen::MatrixXd H_;
    Eigen::MatrixXd& H() {
        return H_;
    }
    
    // measurement covariance matrix
    Eigen::MatrixXd R_;
    Eigen::MatrixXd& R() {
        return R_;
    }

    int ID() {
        return id;
    }
    /**
     * Constructor
     */
    KalmanFilterArrays() {
        id++;
    }
    
    std::string toString();

    void initQ(const float deltaT, const Eigen::VectorXd &theNoise);
    void initQ(const float deltaT, const VectorXd &theStateVector, const VectorXd &theNoise);
    static Eigen::MatrixXd makeQ(const float deltaT, const int theNumberOfMeasurements, const Eigen::VectorXd &theNoise);
    static Eigen::MatrixXd makeQ(const float deltaT, const Eigen::VectorXd &theMeasurements, const Eigen::VectorXd &theNoise);
    static Eigen::MatrixXd updateQ(const float deltaT, const Eigen::VectorXd &theNoise, Eigen::MatrixXd &Q);

    void initF(float deltaT);
    static Eigen::MatrixXd makeF(float deltaT, Eigen::VectorXd theMeasurements);
    static Eigen::MatrixXd makeF(float deltaT, int theNumberOfPositions);
    static Eigen::MatrixXd updateF(const float deltaT, Eigen::MatrixXd &F);

private:
    static int id;
};

class KalmanFilter {
public:
    // state vector
    Eigen::VectorXd& x() {
        return kalmanArrays().x();
    }
    
    // state covariance matrix
    Eigen::MatrixXd& P() {
        return kalmanArrays().P();
    }
    
    // state transition matrix
    Eigen::MatrixXd& F() {
        return kalmanArrays().F();
    }
    
    // process covariance matrix
    Eigen::MatrixXd& Q() {
        return kalmanArrays().Q();
    }
    
    // measurement matrix
    Eigen::MatrixXd& H() {
        return kalmanArrays().H();
    }
    
    // measurement covariance matrix
    Eigen::MatrixXd& R() {
        return kalmanArrays().R();
    }
    
    void initQ(const float deltaT, const VectorXd &theNoise);
    void initF(float deltaT);

    std::string toString() {
        return kalmanArrays().toString();
    }
    
private:
    
    KalmanFilterArrays& kalmanFilterArrays;
    KalmanFilterArrays& kalmanArrays() {
        return kalmanFilterArrays;
    }
    
public:
    /**
     * Constructor
     */
    KalmanFilter(KalmanFilterArrays &theKalmanFilterArrays);
    
    /**
     * Destructor
     */
    virtual ~KalmanFilter();
    
    /**
     * Init Initializes Kalman filter
     * @param x_in Initial state
     * @param P_in Initial state covariance
     * @param F_in Transition matrix
     * @param H_in Measurement matrix
     * @param R_in Measurement covariance matrix
     * @param Q_in Process covariance matrix
     */
    void Init(Eigen::VectorXd &x_in, Eigen::MatrixXd &P_in, Eigen::MatrixXd &F_in,
              Eigen::MatrixXd &H_in, Eigen::MatrixXd &R_in, Eigen::MatrixXd &Q_in);
    
    /**
     * Prediction Predicts the state and the state covariance
     * using the process model
     * @param delta_T Time between k and k+1 in s
     */
    virtual void Predict(const double deltaT, const Eigen::VectorXd &theProcessNoise);
    virtual void Predict();
    
    /**
     * Updates the state by using standard Kalman Filter equations
     * @param z The measurement at k+1
     */
    virtual void Update(const Eigen::VectorXd &z, const bool isRadar);
    
    virtual Eigen::VectorXd zPredicted(const Eigen::VectorXd &x);
    virtual Eigen::VectorXd normalize(Eigen::VectorXd &y);
    
    Eigen::VectorXd measurementUpdate(const Eigen::VectorXd &xPredicted, Eigen::MatrixXd &P, const Eigen::VectorXd &z, const Eigen::MatrixXd &H, const Eigen::MatrixXd &R, const Eigen::MatrixXd &I, const bool isRadar);
    
    /**
     * Updates the state by using Extended Kalman Filter equations
     * @param z The measurement at k+1
     */
    
    void UpdateEKF(const Eigen::VectorXd &z);
    
};

class ExtendedKalmanFilter: protected KalmanFilter {
public:
    
    ExtendedKalmanFilter(KalmanFilterArrays &theKalmanFilterArrays) : KalmanFilter::KalmanFilter(theKalmanFilterArrays) {
    }

    void Predict(const double deltaT, const Eigen::VectorXd &theProcessNoise) {
        KalmanFilter::Predict(deltaT, theProcessNoise);
    }
    
    virtual void Update(const Eigen::VectorXd &z, const bool isRadar);

    virtual Eigen::VectorXd zPredicted(const Eigen::VectorXd &x);
};

class LaserFilter: protected KalmanFilter {
private:
    MatrixXd R_;
    MatrixXd H_;
    
public:
    
    /**
     * Constructor
     */
    LaserFilter(KalmanFilterArrays &theKalmanFilterArrays);
    
    void Predict(const double deltaT, const Eigen::VectorXd &theProcessNoise) {
        KalmanFilter::Predict(deltaT, theProcessNoise);
    }

    void Update(const Eigen::VectorXd &z);

    void updateX(const VectorXd &theLaserMeasurement );

    Eigen::VectorXd zPredicted(const VectorXd &x);
};

class RadarFilter: protected ExtendedKalmanFilter {
private:
    MatrixXd R_;
    
    static MatrixXd CalculateHj(const VectorXd& x_state);

public:
    
    /**
     * Constructor
     */
    RadarFilter(KalmanFilterArrays &theKalmanFilterArrays);
    
    void Predict(const double deltaT, const Eigen::VectorXd &theProcessNoise) {
        ExtendedKalmanFilter::Predict(deltaT, theProcessNoise);
    }
    
    void Update(const Eigen::VectorXd &z, const bool isRadar);

    void updateX(const Eigen::VectorXd &theRadarMeasurement);
    Eigen::VectorXd normalize(Eigen::VectorXd &y);

};

#endif /* KALMAN_FILTER_H_ */
