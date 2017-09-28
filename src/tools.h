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
    
  static void filter(VectorXd &x, MatrixXd &P, const MatrixXd &F, const MatrixXd &Q, const VectorXd &z, const MatrixXd &H, const MatrixXd &R, const MatrixXd &I, const bool isRadar);
  static void predict(VectorXd &x, MatrixXd &P, const MatrixXd &F, const MatrixXd &Q );
  static VectorXd measurementUpdate(const VectorXd &x, MatrixXd &P, const VectorXd &z, const MatrixXd &H, const MatrixXd &R, const MatrixXd &I , const bool isRadar);

  static VectorXd updateX(VectorXd &x, const VectorXd &theMeasurement );
  static VectorXd zPredicted(const VectorXd &x);

  static MatrixXd makeF(float deltaT, VectorXd theMeasurements);
  static MatrixXd makeF(float deltaT, int theNumberOfPositions);
  static MatrixXd updateF(const float deltaT, MatrixXd &F);

  static MatrixXd makeQ(const float deltaT, const int theNumberOfMeasurements, const VectorXd &theNoise);
  static MatrixXd makeQ(const float deltaT, const VectorXd &theMeasurements, const VectorXd &theNoise);
  static MatrixXd updateQ(const float deltaT, const VectorXd &theNoise, MatrixXd &Q);

  static bool areSame(VectorXd a, VectorXd b);
  static bool areSame(MatrixXd a, MatrixXd b);
  static bool areSame(double a, double b);
  static bool isZero(double a);
  static bool isNotZero(double a);

  static double normalizeAngle(double angle) {
        double a = fmod(angle + M_PI, 2 * M_PI);
        return a >= 0 ? (a - M_PI) : (a + M_PI);
  }
    
  static const bool TESTING=false;

  static std::string toString(VectorXd vector);
  static std::string toString(MatrixXd matrix);
  static std::string toStringSize(int rows, int columns);
    
};

#endif /* TOOLS_H_ */
