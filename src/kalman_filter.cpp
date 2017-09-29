#include "kalman_filter.h"
#include "tools.h"
#include <iostream>
#include <sstream>


using Eigen::MatrixXd;
using Eigen::VectorXd;


std::string KalmanFilterArrays::toString() {
    std::ostringstream oss;
        oss << "KalmanFilterArrays::toString-P:" << std::endl
        << "x:" <<Tools::toString(x()) << std::endl
        << "P:" <<Tools::toString(P()) << std::endl
        << "F:" <<Tools::toString(F()) << std::endl
        << "Q:" <<Tools::toString(Q()) << std::endl
        << "H:" <<Tools::toString(F()) << std::endl
        << "R:" <<Tools::toString(Q()) << std::endl;
    return oss.str();
}

//KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x() = x_in;
  P() = P_in;
  F() = F_in;
  H() = H_in;
  R() = R_in;
  Q() = Q_in;
}

KalmanFilter::KalmanFilter(KalmanFilterArrays &theKalmanFilterArrays) : kalmanFilterArrays(theKalmanFilterArrays) {
    if (true || Tools::TESTING) {
        std::cout << "KalmanFilter::KalmanFilter-theKalmanFilterArrays:" << theKalmanFilterArrays.toString() << std::endl;
    }
}

void KalmanFilter::Predict(const double deltaT, const Eigen::VectorXd &theProcessNoise) {
    
    F()=Tools::updateF(deltaT, F());
    Q()=Tools::updateQ(deltaT, theProcessNoise, Q());
    
    Predict();
    if (Tools::TESTING) {
        std::cout << "KalmanFilter::Predict-after predict-P:" << std::endl
        << "x:" <<Tools::toString(x()) << std::endl
        << "P:" <<Tools::toString(P()) << std::endl
        << "F:" <<Tools::toString(F()) << std::endl
        << "Q:" <<Tools::toString(Q()) << std::endl;
    }
}

void KalmanFilter::Predict() {
    /**
     TODO:
     * predict the state
     */
    if (Tools::TESTING)
        std::cout << "Tools::predict-before-x=" << Tools::toString(x()) << ", F=" << Tools::toString(F()) << std::endl;
    x() = F() * x();
    MatrixXd Ft = F().transpose();
    P() = F() * P() * Ft + Q();
    if (Tools::TESTING) std::cout << "Tools::predict-after-F=" << Tools::toString(F()) << std::endl
        << "Q:" << Tools::toString(Q()) << std::endl
        << "P:" << Tools::toString(P()) << std::endl
        << "x:" << Tools::toString(x()) << std::endl;
}

void KalmanFilter::Update(const VectorXd &z, const bool isRadar) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
    int sizeOfStateSpace=x().size();
    MatrixXd I = MatrixXd::Identity(sizeOfStateSpace, sizeOfStateSpace);
    x()=Tools::measurementUpdate(x(), P(), z, H(), R(), I, isRadar);
}

void ExtendedKalmanFilter::Update(const VectorXd &z, const bool isRadar) {
    KalmanFilter::Update(z, isRadar);
}

//void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
//}
