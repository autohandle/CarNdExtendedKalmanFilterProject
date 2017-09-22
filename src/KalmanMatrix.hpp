//
//  MatrixPrinter.hpp
//  ExtendedKF
//
//  Created by David Scott on 9/22/17.
//
//

#ifndef KalmanMatrix_hpp
#define KalmanMatrix_hpp

#include <stdio.h>
#include <string>
#include <iostream>

#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;

class KalmanMatrix {
    
    
    public:
    
    static const bool TESTING=false;
    
    static MatrixXd makeQ(float deltaT, int theNumberOfMeasurements, VectorXd theNoise) {
        if (KalmanMatrix::TESTING) std::cout << "makeQ-theNumberOfMeasurements:" <<  theNumberOfMeasurements << std::endl;
        // for every position there is a velocity
        MatrixXd Q = MatrixXd(theNumberOfMeasurements*2, theNumberOfMeasurements*2);
        
        return KalmanMatrix::updateQ(deltaT, theNoise, &Q);
        MatrixXd x;
        return x;
    }
    
    static MatrixXd makeQ(float deltaT, VectorXd theMeasurements, VectorXd theNoise) {
        if (KalmanMatrix::TESTING) std::cout << "makeQ-theMeasurements" <<  toString(theMeasurements) << std::endl;
        if (KalmanMatrix::TESTING) std::cout << "makeQ-theNoise" <<  toString(theNoise) << std::endl;
        assert(theMeasurements.size()==theNoise.size());
        int numberOfMeasurements = theMeasurements.size();
        if (TESTING) std::cout << "makeQ-numberOfMeasurements:" <<  numberOfMeasurements << std::endl;
        return makeQ(deltaT, numberOfMeasurements, theNoise);
    }
    
    static MatrixXd updateQ(float deltaT, VectorXd theNoise, MatrixXd *Q) {
        // initialize Q
        (*Q).fill(0.);
        if (KalmanMatrix::TESTING) std::cout << "updateQ-deltaT:" <<  deltaT << std::endl;
        // initialialize top half of Q
        double quadConstant=(deltaT*deltaT*deltaT*deltaT)/4.;
        double cubeConstant=(deltaT*deltaT*deltaT)/2.;
        int numberOfMeasurements = (*Q).rows()/2;
        for (int row=0; row<numberOfMeasurements; row++) {
            (*Q)(row, row)=quadConstant*theNoise(row);
            int column=numberOfMeasurements+row;
            (*Q)(row, column)=cubeConstant*theNoise(row);
        }
        if (KalmanMatrix::TESTING) std::cout << "updateQ-updateQ-top half-Q" <<  KalmanMatrix::toString(*Q) << std::endl;
        // initialialize bottom half of Q
        double squareConstant=(deltaT*deltaT);
        for (int row=0; row<numberOfMeasurements; row++) {
            int column=numberOfMeasurements+row;
            (*Q)(column, row)=cubeConstant*theNoise(row);
            (*Q)(column, column)=squareConstant*theNoise(row);
        }
        if (KalmanMatrix::TESTING) std::cout << "updateQ-updateQ-bottom half-Q" <<  KalmanMatrix::toString(*Q) << std::endl;
        return (*Q);
    }
    
    #include <sstream>
    static std::string toString(VectorXd vector) {
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
    
    //#include <typeinfo>
    static std::string toString(MatrixXd matrix) {
        //std::cout << typeid(matrix(0)).name() << '\n';
        //std::cout << decltype(matrix) << '\n';
        std::ostringstream ss;
        ss << "[" << toStringSize(matrix.rows(), matrix.cols()) << "]=";
        for (int r=0; r<matrix.rows(); r++ ) {
            if ((r == 0) && (matrix.rows() > 1)) {
                ss << std::endl;
            }
            if (r != 0 ) ss << "," << std::endl;
            ss << "[" << matrix.row(r) << "]";
        }
        return ss.str();
    }

    static std::string toStringSize(int rows, int columns) {
        std::ostringstream ss;
        ss << rows << "x" << columns;
        return ss.str();
    }

    private:

};

#endif /* KalmanMatrix_hpp */
