#include <iostream>
#include "tools.h"
#include <math.h>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

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
    if (TESTING) std::cout << "Tools::predict-before-x=" << Tools::toString(x) << ", F=" << Tools::toString(F) << std::endl;
    x = F * x;
    MatrixXd Ft = F.transpose();
    P = F * P * Ft + Q;
    if (TESTING) std::cout << "Tools::predict-after-F=" << Tools::toString(F) << std::endl
        << "Q:" << Tools::toString(Q) << std::endl
        << "P:" << Tools::toString(P) << std::endl
        << "x:" << Tools::toString(x) << std::endl;
}

VectorXd Tools::zPredicted(const VectorXd &x) {
    double px = x[0];
    double py = x[1];
    double vx = x[2];
    double vy = x[3];
    
    double rho = sqrt(px*px+py*py);
    assert(isNotZero(rho));
    double phi = atan2(py,px);
    assert (abs(phi) <= M_PI);
    double rhoDot = (px*vx+py*vy)/rho;
    
    VectorXd zPredicted=VectorXd(3);
    zPredicted[0]=rho;
    zPredicted[1]=phi;
    zPredicted[2]=rhoDot;
    
    return zPredicted;
}

VectorXd Tools::measurementUpdate(const VectorXd &xPredicted, MatrixXd &P, const VectorXd &z, const MatrixXd &H, const MatrixXd &R, const MatrixXd &I, const bool isRadar) {
    /*
     * KF Measurement update step
     */
    if (TESTING) {
        std::cout << "Tools::measurementUpdate-xPredicted=" << Tools::toString(xPredicted) << std::endl
                << "H" << Tools::toString(H) << std::endl
                << "P" << Tools::toString(P) << std::endl
                << "R" << Tools::toString(R) << std::endl
                << "isRadar?" << isRadar << std::endl
        ;
    }
    VectorXd measurementPredicted;
    VectorXd y;

    if (isRadar) {
        measurementPredicted = zPredicted(xPredicted);// 0:rho 1:phi 2:rho dot
        y = z-measurementPredicted;
        y(1)=fmod(y(1), M_PI);
        if (TESTING) {
            std::cout << "Tools::measurementUpdate-radar-measurementPredicted=" << Tools::toString(measurementPredicted) << std::endl;
        }
        //if (TESTING) {
        //    std::cout << "Tools::measurementUpdate-y=" << Tools::toString(y) << std::endl
        //    ;
        //}
        //y(1)=normalizeAngle(y(1));
        //if (TESTING) {
        //    std::cout << "Tools::measurementUpdate-y=" << Tools::toString(y) << std::endl
        //    ;
        //}
    } else {
        measurementPredicted = H*xPredicted;
        y = z-measurementPredicted;
        if (TESTING) {
            std::cout << "Tools::measurementUpdate-lidar-measurementPredicted=" << Tools::toString(measurementPredicted) << std::endl;
        }
    }
    MatrixXd Ht = H.transpose();
    if (TESTING) {
        std::cout << "measurementPredicted=" << Tools::toString(measurementPredicted) << std::endl
                  << "z=" << Tools::toString(z) << std::endl
                  << "R=" << Tools::toString(R) << std::endl
                  << "y:" << Tools::toString(y) << std::endl;
        ;
    }
    MatrixXd S = H * P * Ht + R;
    MatrixXd Si = S.inverse();
    MatrixXd K =  P * Ht * Si;
    
    //new state
    VectorXd xNewState = xPredicted + (K * y);
    P = (I - K * H) * P;
    if (TESTING) {
        std::cout << "Tools::measurementUpdate-new state-xNewState=" << Tools::toString(xNewState) << std::endl
                  << "P=" << Tools::toString(P) << std::endl
                  << "H=" << Tools::toString(H) << std::endl
                  << "S=" << Tools::toString(S) << std::endl
                  << "K=" << Tools::toString(K) << std::endl
        ;
    }
    return xNewState;
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