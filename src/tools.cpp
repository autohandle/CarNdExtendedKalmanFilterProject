#include <iostream>
#include "tools.h"
#include <math.h>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

MatrixXd Tools::CalculateHj(const VectorXd& x_state) {
    
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

/**
 TODO:
 * Calculate the RMSE here.
 */

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                       const vector<VectorXd> &ground_truth){
    
    VectorXd rmse(4);
    rmse << 0,0,0,0;
    
    // check the validity of the following inputs:
    //  * the estimation vector size should not be zero
    //  * the estimation vector size should equal ground truth vector size
    if(estimations.size() != ground_truth.size()
       || estimations.size() == 0){
        cout << "Tools::CalculateRMSE-Invalid estimation or ground_truth data" << endl;
        return rmse;
    }
    
    //accumulate squared residuals
    for(unsigned int i=0; i < estimations.size(); ++i){
        
        VectorXd residual = estimations[i] - ground_truth[i];
        
        //coefficient-wise multiplication
        residual = residual.array()*residual.array();
        rmse += residual;
    }
    
    //calculate the mean
    rmse = rmse/estimations.size();
    
    //calculate the squared root
    rmse = rmse.array().sqrt();
    
    //return the result
    return rmse;
}

void Tools::filter(VectorXd &x, MatrixXd &P, const MatrixXd &F, const MatrixXd &Q, const VectorXd &z, const MatrixXd &H, const MatrixXd &R, const MatrixXd &I, const bool isRadar) {
    //YOUR CODE HERE
    predict(x, P, F, Q);
    measurementUpdate(x, P, z, H, R, I, isRadar);
    std::cout << "Tools::filter-x=" << std::endl <<  x << std::endl;
    std::cout << "Tools::filter-P=" << std::endl <<  P << std::endl;
}

void Tools::predict(VectorXd &x, MatrixXd &P, const MatrixXd &F, const MatrixXd &Q ) {
    /*
     * KF Prediction step
     */
    std::cout << "Tools::predict-x=" << Tools::toString(x) << ", F=" << Tools::toString(F) << std::endl;
    x = F * x;
    MatrixXd Ft = F.transpose();
    std::cout << "Tools::predict-F=" << Tools::toString(F) << std::endl
                << "P:" << Tools::toString(P) << std::endl
                << "Q:" << Tools::toString(Q) << std::endl;
    P = F * P * Ft + Q;
}

void Tools::measurementUpdate(VectorXd &x, MatrixXd &P, const VectorXd &z, const MatrixXd &H, const MatrixXd &R, const MatrixXd &I, const bool isRadar) {
    /*
     * KF Measurement update step
     */
    if (TESTING) {
        std::cout << "Tools::measurementUpdate-x=" << Tools::toString(x) << std::endl
                << "H" << Tools::toString(H) << std::endl
                << "P" << Tools::toString(P) << std::endl
                << "R" << Tools::toString(R) << std::endl
                << "isRadar?" << isRadar << std::endl
        ;
    }
    VectorXd measurementPredicted = H * x;
    VectorXd y = z-measurementPredicted; // y[1] is phi for radar: 0:rho 1:phi 2:rho dot
    if (isRadar) {
        if (TESTING) {
            std::cout << "Tools::measurementUpdate-y=" << Tools::toString(y) << std::endl
            ;
        }
        y(1)=normalizeAngle(y(1));
        if (TESTING) {
            std::cout << "Tools::measurementUpdate-y=" << Tools::toString(y) << std::endl
            ;
        }
    }
    MatrixXd Ht = H.transpose();
    if (TESTING) {
        std::cout << "measurementPredicted=" << Tools::toString(measurementPredicted) << std::endl
                  << "z=" << Tools::toString(z) << std::endl
                  << "R=" << Tools::toString(R) << std::endl
        ;
    }
    MatrixXd S = H * P * Ht + R;
    MatrixXd Si = S.inverse();
    MatrixXd K =  P * Ht * Si;
    
    //new state
    x = x + (K * y);
    P = (I - K * H) * P;
    if (TESTING) {
        std::cout << "Tools::measurementUpdate-new state-x=" << Tools::toString(x) << std::endl
        << "P=" << Tools::toString(P) << std::endl
        ;
    }
    
}

VectorXd Tools::updateX(VectorXd &x, const VectorXd &theMeasurement ) {// initialize x from 1st measurement 
    assert(x.size() >= theMeasurement.size());
    for (int i=0; i<theMeasurement.size(); i++) {
        x(i)=theMeasurement(i);
    }
    return x;
}


MatrixXd Tools::makeF(float deltaT, VectorXd theMeasurements) {
    if (TESTING) std::cout << "Tools::makeF-theMeasurements" <<  Tools::toString(theMeasurements) << std::endl;
    int numberOfMeasurements = theMeasurements.size();
    if (TESTING) std::cout << "Tools::makeF-numberOfMeasurements:" <<  numberOfMeasurements << std::endl;
    return makeF(deltaT, numberOfMeasurements);
}

MatrixXd Tools::makeF(float deltaT, int theNumberOfPositions) {
    if (TESTING) std::cout << "Tools::makeF-theNumberOfPositions:" <<  theNumberOfPositions << std::endl;
    // for every position there is a velocity
    MatrixXd F = MatrixXd::Identity(theNumberOfPositions, theNumberOfPositions);
    
    F=updateF(deltaT, F);
    if (TESTING) std::cout << "Tools::makeF-f" <<  toString(F) << std::endl;
    return F;
}

MatrixXd Tools::updateF(const float deltaT, MatrixXd &F) {
    if (TESTING) std::cout << "Tools::updateF-deltaT:" <<  deltaT << std::endl;
    // every position measurement is adjust by its velocity
    int derivativeLocation=F.rows()/2;
    for (int row=0; row<derivativeLocation; row++) {
        int column=derivativeLocation+row;
        if (TESTING) std::cout << "Tools::updateF-F[" <<  "row:" << row << "]=" << toString(F) << std::endl;
        F(row,column)=deltaT;
    }
    
    if (TESTING) std::cout << "Tools::updateF-F" <<  toString(F) << std::endl;
    
    return F;
}

#include <sstream>

MatrixXd Tools::makeQ(const float deltaT, const VectorXd &theStateVector, const VectorXd &theNoise) {
    if (Tools::TESTING) std::cout << "Tools::makeQ-theStateVector" <<  Tools::toString(theStateVector) << std::endl;
    if (Tools::TESTING) std::cout << "Tools::makeQ-theNoise" <<  Tools::toString(theNoise) << std::endl;
    //assert(theMeasurements.size()==theNoise.size()); // assumes #(px,py,...) == #(vx,vy, ...)
    int stateVectorSize = theStateVector.size();
    if (TESTING) std::cout << "Tools::makeQ-stateVectorSize:" <<  stateVectorSize << std::endl;
    return makeQ(deltaT, stateVectorSize, theNoise);
}

MatrixXd Tools::makeQ(const float deltaT, const int theStateVectorSize, const VectorXd &theNoise) {
    if (Tools::TESTING) std::cout << "makeQ-theStateVectorSize:" <<  theStateVectorSize << std::endl;
    // for every position there is a velocity
    MatrixXd Q = MatrixXd(theStateVectorSize, theStateVectorSize);
    return Tools::updateQ(deltaT, theNoise, Q);
    MatrixXd x;
    return x;
}

MatrixXd Tools::updateQ(const float deltaT, const VectorXd &theNoise, MatrixXd &Q) {
    assert (deltaT<1) ;
    // initialize Q
    Q.fill(0.);
    if (Tools::TESTING) {
        std::cout << "updateQ-deltaT:" <<  deltaT << std::endl
        << "theNoise:" << Tools::toString(theNoise) << std::endl;
    }
    // initialialize top half of Q
    double quadConstant=(deltaT*deltaT*deltaT*deltaT)/4.;
    double cubeConstant=(deltaT*deltaT*deltaT)/2.;
    int numberOfStateVariables = Q.rows()/2;
    if (Tools::TESTING) {
        std::cout << "updateQ-cubeConstant:" <<  cubeConstant << std::endl
                    << "quadConstant:" <<  quadConstant << std::endl
                    << "numberOfStateVariables:" <<  numberOfStateVariables << std::endl
        ;

    }
    for (int row=0; row<numberOfStateVariables; row++) {
        Q(row, row)=quadConstant*theNoise(row);
        int column=numberOfStateVariables+row;
        Q(row, column)=cubeConstant*theNoise(row);
    }
    if (Tools::TESTING) std::cout << "updateQ-updateQ-top half-Q" <<  Tools::toString(Q) << std::endl;
    // initialialize bottom half of Q
    double squareConstant=(deltaT*deltaT);
    for (int row=0; row<numberOfStateVariables; row++) {
        int column=numberOfStateVariables+row;
        Q(column, row)=cubeConstant*theNoise(row);
        Q(column, column)=squareConstant*theNoise(row);
    }
    if (Tools::TESTING) std::cout << "updateQ-updateQ-bottom half-Q" <<  Tools::toString(Q) << std::endl;
    return Q;
}


bool Tools::areSame(VectorXd a, VectorXd b) {
    for (int r=0; r<a.rows(); r++) {
        if (!Tools::areSame(a[r],b[r])) {
            std::cout << std::endl
            << "a(" << r << "):" << a(r) << " != "
            << "b(" << r << "):" << b(r)
            << std::endl;
            return false;
        }
    }
    return true;
}

bool Tools::areSame(MatrixXd a, MatrixXd b) {
    for (int r=0; r<a.rows(); r++) {
        for (int c=0; c<a.cols(); c++) {
            if (!Tools::areSame(a(r,c),b(r,c))) {
                std::cout << std::endl
                << "a(" << r << "," << c << "):" << a(r,c) << " != "
                << "b(" << r << "," << c << "):" << b(r,c)
                << std::endl;
                return false;
            }
        }
    }
    return true;
}

static const double EPSILON=1.e-6;

bool Tools::areSame(double a, double b) {
    return fabs(a - b) < EPSILON;
}

bool Tools::isZero(double a) {
    return fabs(a) < EPSILON;
}

bool Tools::isNotZero(double a) {
    return !isZero(a);
}

std::string Tools::toString(VectorXd vector) {
    std::ostringstream ss;
    ss << "[" << vector.rows() << "x_]:";
    for (int r=0; r<vector.rows(); r++ ) {
        if ((r == 0) && (vector.rows() > 1)) {
            ss << std::endl;
        }
        if (r == 0) {
            ss << "[";
        } else {
            ss << "," << std::endl;
        }
        ss << vector(r);
    }
    ss << "]";
    //ss << "[" << toStringSize(vector.rows(), vector.cols()) << "]=" << vector;
    return ss.str();
}

std::string Tools::toString(MatrixXd matrix) {
    //std::cout << typeid(matrix(0)).name() << '\n';
    //std::cout << decltype(matrix) << '\n';
    std::ostringstream ss;
    ss << "[" << Tools::toStringSize(matrix.rows(), matrix.cols()) << "]=";
    for (int r=0; r<matrix.rows(); r++ ) {
        if ((r == 0) && (matrix.rows() > 1)) {
            ss << std::endl;
        }
        if (r != 0 ) ss << "," << std::endl;
        ss << "[" << matrix.row(r) << "]";
    }
    return ss.str();
}

std::string Tools::toStringSize(int rows, int columns) {
    std::ostringstream ss;
    ss << rows << "x" << columns;
    return ss.str();
}