#include <uWS/uWS.h>
#include <iostream>
#include "json.hpp"
#include <math.h>
#include "FusionEKF.h"
#include "tools.h"

using namespace std;

// for convenience
using json = nlohmann::json;

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
std::string hasData(std::string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("]");
  if (found_null != std::string::npos) {
    return "";
  }
  else if (b1 != std::string::npos && b2 != std::string::npos) {
    return s.substr(b1, b2 - b1 + 1);
  }
  return "";
}

VectorXd createEstimateVector(FusionEKF theFusionEKF) {
    double p_x = theFusionEKF.ekf_.x()(0);
    double p_y = theFusionEKF.ekf_.x()(1);
    double v1  = theFusionEKF.ekf_.x()(2);
    double v2 = theFusionEKF.ekf_.x()(3);

    VectorXd estimate(4);
    
    estimate(0) = p_x;
    estimate(1) = p_y;
    estimate(2) = v1;
    estimate(3) = v2;
    
    return estimate;
}

VectorXd createGroundTruthVector(std::string theSensorMeasurement) {
    
    MeasurementPackage measurementPackage;
    
    istringstream iss(theSensorMeasurement);
    
    // reads first element from the current line
    string sensorType;
    iss >> sensorType;
    
    if (sensorType.compare("L") == 0) {
        for (int toss=0; toss<3; toss++) {
            string measurement;
            iss >> measurement;
        }
        
    } else if (sensorType.compare("R") == 0) {
        for (int toss=0; toss<4; toss++) {
            string measurement;
            iss >> measurement;
        }
    }
    
    float x_gt;
    float y_gt;
    float vx_gt;
    float vy_gt;
    iss >> x_gt;
    iss >> y_gt;
    iss >> vx_gt;
    iss >> vy_gt;
    VectorXd gt_values(4);
    gt_values(0) = x_gt;
    gt_values(1) = y_gt;
    gt_values(2) = vx_gt;
    gt_values(3) = vy_gt;
    
    return gt_values;
}

MeasurementPackage createMeasurementPackage(std::string theSensorMeasurement) {
    
    MeasurementPackage measurementPackage;
    
    istringstream iss(theSensorMeasurement);
    long long timestamp;
    
    // reads first element from the current line
    string sensorType;
    iss >> sensorType;
    
    if (sensorType.compare("L") == 0) {
        
        measurementPackage.sensor_type_ = MeasurementPackage::LASER;
        measurementPackage.raw_measurements_ = VectorXd(2);
        float px;
        float py;
        iss >> px;
        iss >> py;
        measurementPackage.raw_measurements_ << px, py;
        iss >> timestamp;
        measurementPackage.timestamp_ = timestamp;
        
    } else if (sensorType.compare("R") == 0) {
        
        measurementPackage.sensor_type_ = MeasurementPackage::RADAR;
        measurementPackage.raw_measurements_ = VectorXd(3);
        float ro;
        float theta;
        float ro_dot;
        iss >> ro;
        iss >> theta;
        iss >> ro_dot;
        measurementPackage.raw_measurements_ << ro,theta, ro_dot;
        iss >> timestamp;
        measurementPackage.timestamp_ = timestamp;
    } else {
        measurementPackage.sensor_type_ = MeasurementPackage::INVALID;
    }
    
    return measurementPackage;
}

int runAsServer(FusionEKF fusionEKF) {
    
    // used to compute the RMSE later
    Tools tools;
    vector<VectorXd> estimations;
    vector<VectorXd> ground_truth;
    
    uWS::Hub h;
    
    h.onMessage([&fusionEKF,&tools,&estimations,&ground_truth](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length, uWS::OpCode opCode) {
        // "42" at the start of the message means there's a websocket message event.
        // The 4 signifies a websocket message
        // The 2 signifies a websocket event
        
        if (length && length > 2 && data[0] == '4' && data[1] == '2')
        {
            
            auto s = hasData(std::string(data));
            if (s != "") {
                
                auto j = json::parse(s);
                
                std::string event = j[0].get<std::string>();
                
                if (event == "telemetry") {
                    // j[1] is the data JSON object
                    
                    // sensor_measurment is just the straight data:
                    // "L 3.122427e-01 5.803398e-01 1477010443000000 6.000000e-01 6.000000e-01 5.199937e+00 0 0 6.911322e-03"
                    string sensor_measurment = j[1]["sensor_measurement"];
                    
                    MeasurementPackage meas_package;
                    istringstream iss(sensor_measurment);
                    long long timestamp;
                    
                    // reads first element from the current line
                    string sensor_type;
                    iss >> sensor_type;
                    
                    if (sensor_type.compare("L") == 0) {
                        meas_package.sensor_type_ = MeasurementPackage::LASER;
                        meas_package.raw_measurements_ = VectorXd(2);
                        float px;
                        float py;
                        iss >> px;
                        iss >> py;
                        meas_package.raw_measurements_ << px, py;
                        iss >> timestamp;
                        meas_package.timestamp_ = timestamp;
                    } else if (sensor_type.compare("R") == 0) {
                        
                        meas_package.sensor_type_ = MeasurementPackage::RADAR;
                        meas_package.raw_measurements_ = VectorXd(3);
                        float ro;
                        float theta;
                        float ro_dot;
                        iss >> ro;
                        iss >> theta;
                        iss >> ro_dot;
                        meas_package.raw_measurements_ << ro,theta, ro_dot;
                        iss >> timestamp;
                        meas_package.timestamp_ = timestamp;
                    }
                    float x_gt;
                    float y_gt;
                    float vx_gt;
                    float vy_gt;
                    iss >> x_gt;
                    iss >> y_gt;
                    iss >> vx_gt;
                    iss >> vy_gt;
                    VectorXd gt_values(4);
                    gt_values(0) = x_gt;
                    gt_values(1) = y_gt;
                    gt_values(2) = vx_gt;
                    gt_values(3) = vy_gt;
                    ground_truth.push_back(gt_values);
                    
                    //Call ProcessMeasurment(meas_package) for Kalman filter
                    fusionEKF.ProcessMeasurement(meas_package);
                    
                    //Push the current estimated x,y positon from the Kalman filter's state vector
                    
                    VectorXd estimate(4);
                    
                    double p_x = fusionEKF.ekf_.x()(0);
                    double p_y = fusionEKF.ekf_.x()(1);
                    double v1  = fusionEKF.ekf_.x()(2);
                    double v2 = fusionEKF.ekf_.x()(3);
                    
                    estimate(0) = p_x;
                    estimate(1) = p_y;
                    estimate(2) = v1;
                    estimate(3) = v2;
                    
                    estimations.push_back(estimate);
                    
                    VectorXd RMSE = tools.CalculateRMSE(estimations, ground_truth);
                    if (Tools::TESTING) cout<<"runAsServer-RMSE: <"<< Tools::toString(RMSE) << "\n";

                    json msgJson;
                    msgJson["estimate_x"] = p_x;
                    msgJson["estimate_y"] = p_y;
                    msgJson["rmse_x"] =  RMSE(0);
                    msgJson["rmse_y"] =  RMSE(1);
                    msgJson["rmse_vx"] = RMSE(2);
                    msgJson["rmse_vy"] = RMSE(3);
                    auto msg = "42[\"estimate_marker\"," + msgJson.dump() + "]";
                    // std::cout << msg << std::endl;
                    ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
                    
                }
            } else {
                
                std::string msg = "42[\"manual\",{}]";
                ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
            }
        }
        
    });
    
    // We don't need this since we're not using HTTP but if it's removed the program
    // doesn't compile :-(
    h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data, size_t, size_t) {
        const std::string s = "<h1>Hello world!</h1>";
        if (req.getUrl().valueLength == 1)
        {
            res->end(s.data(), s.length());
        }
        else
        {
            // i guess this should be done more gracefully?
            res->end(nullptr, 0);
        }
    });
    
    h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
        std::cout << "Connected Now2!!!" << std::endl;
    });
    
    h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code, char *message, size_t length) {
        ws.close();
        std::cout << "Disconnected" << std::endl;
    });
    
    int port = 4567;
    if (h.listen(port))
    {
        std::cout << "Listening to port " << port << std::endl;
    }
    else
    {
        std::cerr << "Failed to listen to port" << std::endl;
        return -1;
    }
    h.run();    return 0;
}

VectorXd compareRSME(VectorXd thePreviousRSME, VectorXd theCurrentRSME) {
    assert (thePreviousRSME.size()==theCurrentRSME.size());
    VectorXd comparison = VectorXd(4);
    for (int row=0; row<thePreviousRSME.size(); row++) {
        comparison(row)=theCurrentRSME(row)-thePreviousRSME(row);
    }
    return comparison;
}


int runAsFileProcessor(FusionEKF theFusionEKF, std::string theFileName) {
    
    ifstream measurementFile;
    measurementFile.open(theFileName, ios::in);
    if (Tools::TESTING) cout<<"runAsFileProcessor-theFileName: <"<< theFileName << ">, is_open? " << measurementFile.is_open() << "\n";
    
    int noise_ax = 5;
    int noise_ay = 5;
    int lineNumber = 0;
    VectorXd noise = VectorXd(2);
    VectorXd previousRSME = Eigen::VectorXd(4);
    previousRSME << 0,0,0,0;
    noise << noise_ax, noise_ay;
    if (Tools::TESTING) std::cout << "noise:" <<  noise << std::endl;
    
    if (measurementFile.is_open()) {
        vector<VectorXd> estimates;
        vector<VectorXd> groundTruthVector;
        if (Tools::TESTING) cout<<"runAsFileProcessor-theFileName: "<< theFileName << "\n";
        string measurementLine;
        while ( getline (measurementFile, measurementLine) ){
            lineNumber++;
            if (Tools::TESTING) {
                cout << "l-------------------------------------------" << "\n"
                << lineNumber << ": <" << measurementLine << ">\n"
                << "l-------------------------------------------" << "\n";
            }
            MeasurementPackage measurementPackage = createMeasurementPackage(measurementLine);
            if (measurementPackage.sensor_type_ == MeasurementPackage::RADAR || measurementPackage.sensor_type_ == MeasurementPackage::LASER) {
                theFusionEKF.ProcessMeasurement(measurementPackage);
                
                VectorXd groundTruthValues = createGroundTruthVector(measurementLine);
                if (Tools::TESTING) cout<<"runAsFileProcessor-groundTruthValues: <"<< Tools::toString(groundTruthValues) << "\n";
                groundTruthVector.push_back(groundTruthValues);
                VectorXd estimate = createEstimateVector(theFusionEKF);
                if (Tools::TESTING) cout<<"runAsFileProcessor-estimate: <"<< Tools::toString(estimate) << "\n";
                estimates.push_back(estimate);
                VectorXd rsme = Tools::CalculateRMSE(estimates, groundTruthVector);
                if (Tools::TESTING) cout<<"runAsFileProcessor-rsme: <"<< Tools::toString(rsme) << "\n";
                if (Tools::TESTING) {
                    cout << "r-------------------------------------------" << "\n"
                    << "<" << Tools::toString(compareRSME(previousRSME, rsme)) << ">\n"
                    << "-r------------------------------------------" << "\n";
                }
                previousRSME=rsme;
            } else {
                cout << "runAsFileProcessor-unknown-measurement_pack.sensor_type_:" << measurementPackage.sensor_type_ << endl;
            }
        }
    } else {
        cout<<"runAsFileProcessor-theFileName: <"<< theFileName << "> failed to open\n";
        return 1;
    }
    measurementFile.close();
    return 0;
}

int main(int argc, char *argv[] )
{
    // Create a Kalman Filter instance
    FusionEKF fusionEKF;

    int status=0;
    
    if (Tools::TESTING) cout<<"argc: "<< argc <<"\n";
    if (argc>1) {
        if (Tools::TESTING) cout<<"argv[0]: "<< argv[0] <<"\n";
        status=runAsFileProcessor(fusionEKF, argv[1]);
    } else {
        status=runAsServer(fusionEKF);
    }
  if (Tools::TESTING) cout<< "status: " << status <<"\n";
  return(status);
}
