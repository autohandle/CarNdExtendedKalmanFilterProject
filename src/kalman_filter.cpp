#include "kalman_filter.h"
#include "tools.h"
#include <iostream>
#include <sstream>


using Eigen::MatrixXd;
using Eigen::VectorXd;

int KalmanFilterArrays::id = 0;

std::string KalmanFilterArrays::toString() {
    std::ostringstream oss;
        oss << "KalmanFilterArrays::toString:" << std::endl
        << "id:" << id << std::endl
        << "x:" <<Tools::toString(x()) << std::endl
        << "P:" <<Tools::toString(P()) << std::endl
        << "F:" <<Tools::toString(F()) << std::endl
        << "Q:" <<Tools::toString(Q()) << std::endl
        << "H:" <<Tools::toString(F()) << std::endl
        << "R:" <<Tools::toString(Q()) << std::endl;
    return oss.str();
}

void KalmanFilterArrays::initQ(const float deltaT, const VectorXd &theNoise) {
    //assert(x().size()>0);
    Q()=KalmanFilterArrays::makeQ(deltaT, x(), theNoise);
    if (true || Tools::TESTING) std::cout << "KalmanFilterArrays::initQ:" <<  toString() << std::endl;
}

void KalmanFilterArrays::initQ(const float deltaT, const VectorXd &theStateVector, const VectorXd &theNoise) {
    assert(theStateVector.size()>0);
    Q()=KalmanFilterArrays::makeQ(deltaT, theStateVector, theNoise);
    if (true || Tools::TESTING) std::cout << "KalmanFilterArrays::initQ:" <<  toString() << std::endl;
}

MatrixXd KalmanFilterArrays::makeQ(const float deltaT, const VectorXd &theStateVector, const VectorXd &theNoise) {
    if (Tools::TESTING) std::cout << "Tools::makeQ-theStateVector" <<  Tools::toString(theStateVector) << std::endl;
    if (Tools::TESTING) std::cout << "Tools::makeQ-theNoise" <<  Tools::toString(theNoise) << std::endl;
    //assert(theMeasurements.size()==theNoise.size()); // assumes #(px,py,...) == #(vx,vy, ...)
    int stateVectorSize = theStateVector.size();
    if (Tools::TESTING) std::cout << "Tools::makeQ-stateVectorSize:" <<  stateVectorSize << std::endl;
    return makeQ(deltaT, stateVectorSize, theNoise);
}

MatrixXd KalmanFilterArrays::makeQ(const float deltaT, const int theStateVectorSize, const VectorXd &theNoise) {
    if (Tools::TESTING) std::cout << "makeQ-theStateVectorSize:" <<  theStateVectorSize << std::endl;
    // for every position there is a velocity
    MatrixXd Q = MatrixXd(theStateVectorSize, theStateVectorSize);
    return KalmanFilterArrays::updateQ(deltaT, theNoise, Q);
    MatrixXd x;
    return x;
}

MatrixXd KalmanFilterArrays::updateQ(const float deltaT, const VectorXd &theNoise, MatrixXd &Q) {
    // initialize Q
    Q.fill(0.);
    assert (deltaT<1) ;
    // initialialize top half of Q
    double quadConstant=(deltaT*deltaT*deltaT*deltaT)/4.;
    double cubeConstant=(deltaT*deltaT*deltaT)/2.;
    int numberOfStateVariables = Q.rows()/2;
    for (int row=0; row<numberOfStateVariables; row++) {
        Q(row, row)=quadConstant*theNoise(row);
        int column=numberOfStateVariables+row;
        Q(row, column)=cubeConstant*theNoise(row);
    }
    // initialialize bottom half of Q
    double squareConstant=(deltaT*deltaT);
    for (int row=0; row<numberOfStateVariables; row++) {
        int column=numberOfStateVariables+row;
        Q(column, row)=cubeConstant*theNoise(row);
        Q(column, column)=squareConstant*theNoise(row);
    }
    if (Tools::TESTING) std::cout << "KalmanFilterArrays::updateQ" <<  Tools::toString(Q) << std::endl;
    
    return Q;
}

void KalmanFilterArrays::initF(float deltaT) {
    F()=makeF(deltaT, x());
}

MatrixXd KalmanFilterArrays::makeF(float deltaT, VectorXd theMeasurements) {
    if (Tools::TESTING) std::cout << "KalmanFilterArrays::makeF-theMeasurements" <<  Tools::toString(theMeasurements) << std::endl;
    int numberOfMeasurements = theMeasurements.size();
    if (Tools::TESTING) std::cout << "KalmanFilterArrays::makeF-numberOfMeasurements:" <<  numberOfMeasurements << std::endl;
    return makeF(deltaT, numberOfMeasurements);
}

MatrixXd KalmanFilterArrays::makeF(float deltaT, int theNumberOfPositions) {
    if (Tools::TESTING) std::cout << "KalmanFilterArrays::makeF-theNumberOfPositions:" <<  theNumberOfPositions << std::endl;
    // for every position there is a velocity
    MatrixXd F = MatrixXd::Identity(theNumberOfPositions, theNumberOfPositions);
    
    F=KalmanFilterArrays::updateF(deltaT, F);
    if (Tools::TESTING) std::cout << "KalmanFilterArrays::makeF-f" <<  Tools::toString(F) << std::endl;
    return F;
}

MatrixXd KalmanFilterArrays::updateF(const float deltaT, MatrixXd &F) {
    if (Tools::TESTING) std::cout << "KalmanFilterArrays::updateF-deltaT:" <<  deltaT << std::endl;
    // every position measurement is adjust by its velocity
    int derivativeLocation=F.rows()/2;
    for (int row=0; row<derivativeLocation; row++) {
        int column=derivativeLocation+row;
        if (Tools::TESTING) std::cout << "KalmanFilterArrays::updateF-F[" <<  "row:" << row << "]=" << Tools::toString(F) << std::endl;
        F(row,column)=deltaT;
    }
    
    if (Tools::TESTING) std::cout << "KalmanFilterArrays::updateF-F" <<  Tools::toString(F) << std::endl;
    
    return F;
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
    
    F()=KalmanFilterArrays::updateF(deltaT, F());
    Q()=KalmanFilterArrays::updateQ(deltaT, theProcessNoise, Q());
    
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
    if (true || Tools::TESTING) std::cout << "KalmanFilter::Predict-before:" << toString() << std::endl;
    x() = F() * x();
    MatrixXd Ft = F().transpose();
    P() = F() * P() * Ft + Q();
    if (Tools::TESTING) std::cout << "KalmanFilter::Predict-after-F=" << Tools::toString(F()) << std::endl
        << "Q:" << Tools::toString(Q()) << std::endl
        << "P:" << Tools::toString(P()) << std::endl
        << "x:" << Tools::toString(x()) << std::endl;
}

void KalmanFilter::Update(const VectorXd &z, const bool isRadar) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
    if (Tools::TESTING) std::cout << "KalmanFilter::Update:" <<  toString() << std::endl;
    /**
     Initialize state.
     */
    int sizeOfStateSpace=x().size();
    MatrixXd I = MatrixXd::Identity(sizeOfStateSpace, sizeOfStateSpace);
    x()=measurementUpdate(x(), P(), z, H(), R(), I, isRadar);
}

Eigen::VectorXd KalmanFilter::zPredicted(const Eigen::VectorXd &x) {
    assert(1==0);
}

Eigen::VectorXd KalmanFilter::normalize(Eigen::VectorXd &y) {
    return y;
}

VectorXd KalmanFilter::measurementUpdate(const VectorXd &xPredicted, MatrixXd &P, const VectorXd &z, const MatrixXd &H, const MatrixXd &R, const MatrixXd &I, const bool isRadar) {
    /*
     * KF Measurement update step
     */
    if (Tools::TESTING) {
        std::cout << "KalmanFilter::measurementUpdate-xPredicted=" << Tools::toString(xPredicted) << std::endl
        << "H" << Tools::toString(H) << std::endl
        << "P" << Tools::toString(P) << std::endl
        << "R" << Tools::toString(R) << std::endl
        << "isRadar?" << isRadar << std::endl
        ;
    }
    VectorXd y = z-zPredicted(xPredicted/* 0:rho 1:phi 2:rho dot OR 0:px, 1:py, 2:vx, 3:vy*/);
    y = normalize(y);

    MatrixXd Ht = H.transpose();
    if (Tools::TESTING) {
        std::cout << "KalmanFilter::measurementUpdate:" << toString() << std::endl
        << "z:" << Tools::toString(z) << std::endl
        << "y:" << Tools::toString(y)
        ;
    }
    MatrixXd S = H * P * Ht + R;
    MatrixXd Si = S.inverse();
    MatrixXd K =  P * Ht * Si;
    
    //new state
    VectorXd xNewState = xPredicted + (K * y);
    P = (I - K * H) * P;
    if (Tools::TESTING) {
        std::cout << "Tools::measurementUpdate-new state-xNewState=" << Tools::toString(xNewState) << std::endl
        << "P=" << Tools::toString(P) << std::endl
        << "H=" << Tools::toString(H) << std::endl
        << "S=" << Tools::toString(S) << std::endl
        << "K=" << Tools::toString(K) << std::endl
        ;
    }
    return xNewState;
}

void KalmanFilter::initQ(const float deltaT, const VectorXd &theNoise) {
    kalmanArrays().initQ(deltaT, theNoise);
}

void KalmanFilter::initF(const float deltaT) {
    kalmanArrays().initF(deltaT);
}

void ExtendedKalmanFilter::Update(const VectorXd &z, const bool isRadar) {
    KalmanFilter::Update(z, isRadar);
}

VectorXd ExtendedKalmanFilter::zPredicted(const VectorXd &x) {
    double px = x[0];
    double py = x[1];
    double vx = x[2];
    double vy = x[3];
    
    double rho = sqrt(px*px+py*py);
    assert(Tools::isNotZero(rho));
    double phi = atan2(py,px);
    assert (abs(phi) <= M_PI);
    double rhoDot = (px*vx+py*vy)/rho;
    
    VectorXd zPredicted=VectorXd(3);
    zPredicted[0]=rho;
    zPredicted[1]=phi;
    zPredicted[2]=rhoDot;
    
    return zPredicted;
}

LaserFilter::LaserFilter(KalmanFilterArrays &theKalmanFilterArrays)  : KalmanFilter::KalmanFilter(theKalmanFilterArrays){
    R_ = MatrixXd(2, 2); // px, py
    R_ <<   0.0225, 0,
            0,      0.0225;
    
    H_ = MatrixXd(2, 4);
    H_ <<   1, 0, 0, 0,  // measure vector is: px py
            0, 1, 0, 0;
}

void LaserFilter::Update(const Eigen::VectorXd &z) {
    H()=H_;
    R()=R_;
    KalmanFilter::Update(z, false);
}

void LaserFilter::updateX(const VectorXd &theLaserMeasurement ) {// initialize x from laser measurement
    assert(x().size() >= theLaserMeasurement.size());
    for (int i=0; i<theLaserMeasurement.size(); i++) {
        x()(i)=theLaserMeasurement(i);
    }
}

VectorXd LaserFilter::zPredicted(const VectorXd &x) {
    return H()*x;
}

RadarFilter::RadarFilter(KalmanFilterArrays &theKalmanFilterArrays)  : ExtendedKalmanFilter::ExtendedKalmanFilter(theKalmanFilterArrays) {
    R_ = MatrixXd(3, 3); // rho, phi, rhoDot
    R_ <<   0.09,       0,      0,
            0,          0.0009, 0,
            0,          0,      0.09;
}

void RadarFilter::Update(const Eigen::VectorXd &z, const bool isRadar) {
    H()=CalculateHj(x());
    R()=R_;
    ExtendedKalmanFilter::Update(z, isRadar);
}


MatrixXd RadarFilter::CalculateHj(const VectorXd& x_state) {
    
    MatrixXd Hj(3,4);
    //recover state parameters
    float px = x_state(0);
    float py = x_state(1);
    float vx = x_state(2);
    float vy = x_state(3);
    
    //pre-compute a set of terms to avoid repeated calculation
    float c1 = px*px+py*py;
    float c2 = sqrt(c1);
    float c3 = (c1*c2);
    
    //check division by zero
    if(fabs(c1) < 0.0001){
        cout << "Tools::CalculateHj- () - Error - Division by Zero" << endl;
        return Hj;
    }
    
    //compute the Jacobian matrix
    Hj << (px/c2), (py/c2), 0, 0,
		  -(py/c1), (px/c1), 0, 0,
		  py*(vx*py - vy*px)/c3, px*(px*vy - py*vx)/c3, px/c2, py/c2;
    
    return Hj;
}

void RadarFilter::updateX(const VectorXd &theRadarMeasurement) {
    double rho = theRadarMeasurement[0];
    assert(Tools::isNotZero(rho));
    double phi = theRadarMeasurement[1];
    assert (abs(phi) <= M_PI);
    double rhoDot = theRadarMeasurement[2];
    double px = rho*sin(phi);
    double py = rho*cos(phi);
    double vx = rhoDot*sin(phi);
    double vy = rhoDot*sin(phi);
    
    VectorXd stateVector=VectorXd(4);
    x() << px, py, vx, vy;
}

Eigen::VectorXd RadarFilter::normalize(Eigen::VectorXd &y) {
    y(1)=fmod(y(1), M_PI);
    return y;
}

//void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
//}
