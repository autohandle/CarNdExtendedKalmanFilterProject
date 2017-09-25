#ifndef TOOLS_H_
#define TOOLS_H_
#include <vector>
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

class Tools {
public:
  /**
  * Constructor.
  */
    
    Tools();

  /**
  * Destructor.
  */
  virtual ~Tools();

  /**
  * A helper method to calculate RMSE.
  */
  static VectorXd CalculateRMSE(const vector<VectorXd> &estimations, const vector<VectorXd> &ground_truth);

  /**
  * A helper method to calculate Jacobians.
  */
  static MatrixXd CalculateHj(const VectorXd& x_state);
    
  static void filter(VectorXd &x, MatrixXd &P, MatrixXd F, MatrixXd Q, VectorXd z, MatrixXd H, MatrixXd R, MatrixXd I);
  static void predict(VectorXd &x, MatrixXd &P, MatrixXd F, MatrixXd Q );
  static void measurementUpdate(VectorXd &x, MatrixXd &P, VectorXd z, MatrixXd H, MatrixXd R, MatrixXd I );

  static VectorXd updateX(VectorXd &x, const VectorXd &theMeasurement );

  static MatrixXd makeF(float deltaT, VectorXd theMeasurements);
  static MatrixXd makeF(float deltaT, int theNumberOfPositions);
  static MatrixXd updateF(float deltaT, MatrixXd *F);

  static MatrixXd makeQ(float deltaT, int theNumberOfMeasurements, VectorXd theNoise);
  static MatrixXd makeQ(float deltaT, VectorXd theMeasurements, VectorXd theNoise);
  static MatrixXd updateQ(float deltaT, VectorXd theNoise, MatrixXd *Q);

  static bool areSame(VectorXd a, VectorXd b);
  static bool areSame(MatrixXd a, MatrixXd b);
  static bool areSame(double a, double b);

  static const bool TESTING=true;

  static std::string toString(VectorXd vector);
  static std::string toString(MatrixXd matrix);
  static std::string toStringSize(int rows, int columns);
    
};

#endif /* TOOLS_H_ */
