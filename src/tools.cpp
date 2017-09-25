#include <iostream>
#include "tools.h"

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
        cout << "CalculateJacobian () - Error - Division by Zero" << endl;
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
        cout << "Invalid estimation or ground_truth data" << endl;
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

void Tools::filter(VectorXd &x, MatrixXd &P, MatrixXd F, MatrixXd Q, VectorXd z, MatrixXd H, MatrixXd R, MatrixXd I) {
    //YOUR CODE HERE
    predict(x, P, F, Q);
    measurementUpdate(x, P, z, H, R, I);
    std::cout << "x=" << std::endl <<  x << std::endl;
    std::cout << "P=" << std::endl <<  P << std::endl;
}

void Tools::predict(VectorXd &x, MatrixXd &P, MatrixXd F, MatrixXd Q ) {
    /*
     * KF Prediction step
     */
    std::cout << "x=" << Tools::toString(x) << ", F=" << Tools::toString(F) << std::endl;
    x = F * x;
    MatrixXd Ft = F.transpose();
    std::cout << "Ft=" << Tools::toString(Ft) << ", P=" << Tools::toString(P) << ", Q=" << Tools::toString(Q) << std::endl;
    P = F * P * Ft + Q;
}

void Tools::measurementUpdate(VectorXd &x, MatrixXd &P, VectorXd z, MatrixXd H, MatrixXd R, MatrixXd I) {
    /*
     * KF Measurement update step
     */
    std::cout << "measurementUpdate-x=" << Tools::toString(x) << "," << std::endl <<
        "z=" << Tools::toString(z) << std::endl << "," <<
        "H" << Tools::toString(H) << std::endl;
    VectorXd y = z - H * x;
    MatrixXd Ht = H.transpose();
    std::cout << "measurementUpdate-y=" << Tools::toString(y) << "," << std::endl <<
    "R=" << Tools::toString(R) << std::endl << "," <<
    "Ht" << Tools::toString(Ht) << std::endl;
    MatrixXd S = H * P * Ht + R;
    MatrixXd Si = S.inverse();
    MatrixXd K =  P * Ht * Si;
    
    //new state
    x = x + (K * y);
    std::cout << "measurementUpdate-I=" << Tools::toString(I) << "," << std::endl <<
    "K=" << Tools::toString(K) << std::endl;
    P = (I - K * H) * P;
    
}

VectorXd Tools::updateX(VectorXd &x, const VectorXd &theMeasurement ) {
    assert(x.size() >= theMeasurement.size());
    for (int i=0; i<theMeasurement.size(); i++) {
        x(i)=theMeasurement(i);
    }
    return x;
}


MatrixXd Tools::makeF(float deltaT, VectorXd theMeasurements) {
    if (TESTING) std::cout << "makeF-theMeasurements" <<  Tools::toString(theMeasurements) << std::endl;
    int numberOfMeasurements = theMeasurements.size();
    if (TESTING) std::cout << "makeF-numberOfMeasurements:" <<  numberOfMeasurements << std::endl;
    return makeF(deltaT, numberOfMeasurements);
}

MatrixXd Tools::makeF(float deltaT, int theNumberOfPositions) {
    if (TESTING) std::cout << "makeF-theNumberOfPositions:" <<  theNumberOfPositions << std::endl;
    // for every position there is a velocity
    MatrixXd f = MatrixXd::Identity(theNumberOfPositions, theNumberOfPositions);
    
    f=updateF(deltaT, &f);
    if (TESTING) std::cout << "makeF-f" <<  toString(f) << std::endl;
    return f;
}

MatrixXd Tools::updateF(float deltaT, MatrixXd *F) {
    if (TESTING) std::cout << "updateF-deltaT" <<  deltaT << std::endl;
    // every position measurement is adjust by its velocity
    int derivativeLocation=(*F).rows()/2;
    for (int row=0; row<derivativeLocation; row++) {
        int column=derivativeLocation+row;
        if (TESTING) std::cout << "updateF-F" <<  "row:" << row << "=" << toString(*F) << std::endl;
        (*F)(row,column)=deltaT;
    }
    
    if (TESTING) std::cout << "updateF-F" <<  toString(*F) << std::endl;
    
    return (*F);
}

#include <sstream>

MatrixXd Tools::makeQ(float deltaT, VectorXd theMeasurements, VectorXd theNoise) {
    if (Tools::TESTING) std::cout << "makeQ-theMeasurements" <<  Tools::toString(theMeasurements) << std::endl;
    if (Tools::TESTING) std::cout << "makeQ-theNoise" <<  Tools::toString(theNoise) << std::endl;
    //assert(theMeasurements.size()==theNoise.size()); // assumes #(px,py,...) == #(vx,vy, ...)
    int numberOfMeasurements = theMeasurements.size();
    if (TESTING) std::cout << "makeQ-numberOfMeasurements:" <<  numberOfMeasurements << std::endl;
    return makeQ(deltaT, numberOfMeasurements, theNoise);
}

MatrixXd Tools::makeQ(float deltaT, int theNumberOfMeasurements, VectorXd theNoise) {
    if (Tools::TESTING) std::cout << "makeQ-theNumberOfMeasurements:" <<  theNumberOfMeasurements << std::endl;
    // for every position there is a velocity
    MatrixXd Q = MatrixXd(theNumberOfMeasurements, theNumberOfMeasurements);
    
    return Tools::updateQ(deltaT, theNoise, &Q);
    MatrixXd x;
    return x;
}

MatrixXd Tools::updateQ(float deltaT, VectorXd theNoise, MatrixXd *Q) {
    // initialize Q
    (*Q).fill(0.);
    if (Tools::TESTING) std::cout << "updateQ-deltaT:" <<  deltaT << std::endl;
    // initialialize top half of Q
    double quadConstant=(deltaT*deltaT*deltaT*deltaT)/4.;
    double cubeConstant=(deltaT*deltaT*deltaT)/2.;
    int numberOfMeasurements = (*Q).rows()/2;
    for (int row=0; row<numberOfMeasurements; row++) {
        (*Q)(row, row)=quadConstant*theNoise(row);
        int column=numberOfMeasurements+row;
        (*Q)(row, column)=cubeConstant*theNoise(row);
    }
    if (Tools::TESTING) std::cout << "updateQ-updateQ-top half-Q" <<  Tools::toString(*Q) << std::endl;
    // initialialize bottom half of Q
    double squareConstant=(deltaT*deltaT);
    for (int row=0; row<numberOfMeasurements; row++) {
        int column=numberOfMeasurements+row;
        (*Q)(column, row)=cubeConstant*theNoise(row);
        (*Q)(column, column)=squareConstant*theNoise(row);
    }
    if (Tools::TESTING) std::cout << "updateQ-updateQ-bottom half-Q" <<  Tools::toString(*Q) << std::endl;
    return (*Q);
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

static const double EPSILON=0.0001;

bool Tools::areSame(double a, double b) {
    return fabs(a - b) < EPSILON;
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