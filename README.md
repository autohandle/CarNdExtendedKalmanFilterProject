# Extended Kalman Filters

#### Compiling
##### Code must compile without errors with cmake and make.

``` shell
    (carnd-term1) Software-Engineerings-MacBook-Pro:tmp david$ git clone https://github.com/autohandle/CarNdExtendedKalmanFilterProject.git
    Cloning into 'CarNdExtendedKalmanFilterProject'...
    remote: Counting objects: 813, done.
    remote: Compressing objects: 100% (537/537), done.
    remote: Total 813 (delta 278), reused 796 (delta 261), pack-reused 0
    Receiving objects: 100% (813/813), 2.56 MiB | 2.27 MiB/s, done.
    Resolving deltas: 100% (278/278), done.
    Checking connectivity... done.
    (carnd-term1) Software-Engineerings-MacBook-Pro:tmp david$ cd CarNdExtendedKalmanFilterProject/
    (carnd-term1) Software-Engineerings-MacBook-Pro:CarNdExtendedKalmanFilterProject david$ mkdir build
    (carnd-term1) Software-Engineerings-MacBook-Pro:CarNdExtendedKalmanFilterProject david$ cd build
    (carnd-term1) Software-Engineerings-MacBook-Pro:build david$ cmake ..
    -- The C compiler identification is AppleClang 7.3.0.7030031
    -- The CXX compiler identification is AppleClang 7.3.0.7030031
    -- Check for working C compiler: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc
    -- Check for working C compiler: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc -- works
    -- Detecting C compiler ABI info
    -- Detecting C compiler ABI info - done
    -- Detecting C compile features
    -- Detecting C compile features - done
    -- Check for working CXX compiler: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++
    -- Check for working CXX compiler: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ -- works
    -- Detecting CXX compiler ABI info
    -- Detecting CXX compiler ABI info - done
    -- Detecting CXX compile features
    -- Detecting CXX compile features - done
    -- Configuring done
    -- Generating done
    -- Build files have been written to: /tmp/CarNdExtendedKalmanFilterProject/build
    (carnd-term1) Software-Engineerings-MacBook-Pro:build david$ make
    Scanning dependencies of target ExtendedKF
    [ 20%] Building CXX object CMakeFiles/ExtendedKF.dir/src/main.cpp.o
    [ 40%] Building CXX object CMakeFiles/ExtendedKF.dir/src/tools.cpp.o
    [ 60%] Building CXX object CMakeFiles/ExtendedKF.dir/src/FusionEKF.cpp.o
    [ 80%] Building CXX object CMakeFiles/ExtendedKF.dir/src/kalman_filter.cpp.o
    [100%] Linking CXX executable ExtendedKF
    ld: warning: directory not found for option '-L/usr/local/Cellar/libuv/1.11.0/lib'
    [100%] Built target ExtendedKF
    (carnd-term1) Software-Engineerings-MacBook-Pro:build david$ ./ExtendedKF
    Listening to port 4567
```
#### Accuracy
##### RMSE less than [.11, .11, 0.52, 0.52]

Your algorithm will be run against Dataset 1 in the simulator which is the same as "data/obj_pose-laser-radar-synthetic-input.txt" in the repository. We'll collect the positions that your algorithm outputs and compare them to ground truth data. Your px, py, vx, and vy RMSE should be less than or equal to the values [.11, .11, 0.52, 0.52].

The RSME obtained for data set 1:
``` 
[0.0964425, 0.0852905, 0.415426, 0.431636]
```

The RSME obtained for data set 2:
```
[0.0722488, 0.0968095, 0.355655, 0.519296]
```

#### Follows the Correct Algorithm

##### follow the algorithm

The code is structured to execute the correct algrorithm. The `LaserFilter` is a `KalmanFilter`

```C++
class LaserFilter: protected KalmanFilter
```

While the `RadarFilter` is an `ExtendedKalmanFilter` and the `ExtendedKalmanFilter` is a `KalmanFilter`.

```C++
class ExtendedKalmanFilter: protected KalmanFilter
class RadarFilter: protected ExtendedKalmanFilter
```
The `KalmanFilter` implements the basic Predict
```C++
void KalmanFilter::Predict() {
    x() = F() * x();
    MatrixXd Ft = F().transpose();
    P() = F() * P() * Ft + Q();
}
```
used by both `RadarFilter` and `LaserFilter`.

`KalmanFilter` also implements the basic Update
``` C++
void KalmanFilter::Update(const VectorXd &z, const bool isRadar) {
    MatrixXd I = MatrixXd::Identity(sizeOfStateSpace, sizeOfStateSpace);
    x()=measurementUpdate(x(), P(), z, H(), R(), I, isRadar);
}

VectorXd KalmanFilter::measurementUpdate(const VectorXd &xPredicted, MatrixXd &P, const VectorXd &z, const MatrixXd &H, const MatrixXd &R, const MatrixXd &I, const bool isRadar) {
    VectorXd y = z-zPredicted(xPredicted/* 0:rho 1:phi 2:rho dot OR 0:px, 1:py, 2:vx, 3:vy*/);
    y = normalize(y);

    MatrixXd Ht = H.transpose();
    MatrixXd S = H * P * Ht + R;
    MatrixXd Si = S.inverse();
    MatrixXd K =  P * Ht * Si;
    
    VectorXd xNewState = xPredicted + (K * y);
    P = (I - K * H) * P;
    return xNewState;
}
```
However, `Update` virtualizes out `zPredicted` and `normalize`, so that the `RadarFilter` can override and convert the current state vector into polar coordinates:
```C++
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
```
and also `normalize` the resulting polar coordinates from the difference
``` C++
Eigen::VectorXd RadarFilter::normalize(Eigen::VectorXd &y) {
    y(1)=fmod(y(1), M_PI);
    return y;
}
```

##### first measurements to initialize

In `FusionEKF::ProcessMeasurement`, the program checks if has been initialized
``` C++
    if (!is_initialized_) {
```
then checks for the type of measurement in the `MeasurementPackage` before calling the correct initializer
``` C++
        if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
            radarFilter.updateX(measurement_pack.raw_measurements_);
        } else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
            laserFilter.updateX(measurement_pack.raw_measurements_);
        }
```
and then blocks further initialization calls:
``` C++
        is_initialized_ = true;
```

the `LaserFilter` simply copies the measurement vector into the state vector:
``` C++
void LaserFilter::updateX(const VectorXd &theLaserMeasurement ) {// initialize x from laser measurement
    assert(x().size() >= theLaserMeasurement.size());
    for (int i=0; i<theLaserMeasurement.size(); i++) {
        x()(i)=theLaserMeasurement(i);
    }
}
```
while the `RadarFilter` converts the measurement vector from polar coordinates before initializing the state vector:
``` C++
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
```

##### predict object position to the current timestep and then update

After the initiaizing measurement, all subsequent measurements call the `Predict` and then the `Update` step in that order

```C++
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
        radarFilter.Predict(dt, processNoise);
        radarFilter.Update(measurement_pack.raw_measurements_, true);
    } else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
        // Laser updates
        laserFilter.Predict(dt, processNoise);
        laserFilter.Update(measurement_pack.raw_measurements_);
    }
```

##### calls the correct measurement function

The `MeasurementPackage.sensor_type_` is used to call the corrent measurement function

```C++
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
        radarFilter.Predict(dt, processNoise);
        radarFilter.Update(measurement_pack.raw_measurements_, true);
    } else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
        // Laser updates
        laserFilter.Predict(dt, processNoise);
        laserFilter.Update(measurement_pack.raw_measurements_);
    }
```

#### Code Efficiency

##### Your algorithm should avoid unnecessary calculations.

The program was restructured to use object inheritance to improved comprehension and reduce the number of if-statement checks.

`FusionEKF` could have treated both `LaserFilter` and `RadarFilter` as a `KalmanFilter`, both pulled from a hash table based on `MeasurementPackage.sensor_type_`, but I didn't think that a hash table improved the program enough for the confusion it would have added.

#### Video Implementation

[Data Set 1 Video](https://s3.amazonaws.com/autohandle.com/video/VehicleDetectionDataSet1.mp4)

[Data Set 2 Video](https://s3.amazonaws.com/autohandle.com/video/VehicleDetectionDataSet1.mp4)

The video was created by using a [screen recording tool](http://www.idownloadblog.com/2016/02/26/how-to-record-part-of-mac-screen-quicktime/).
