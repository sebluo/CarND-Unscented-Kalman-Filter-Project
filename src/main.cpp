
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <stdlib.h>
#include "Eigen/Dense"
#include "ukf.h"
#include "ground_truth_package.h"
#include "measurement_package.h"
#include "tools.h"




using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;
#define DEBUG 0
#define TUNING_MODE 0

void check_arguments(int argc, char* argv[]) {
  string usage_instructions = "Usage instructions: ";
  usage_instructions += argv[0];
  usage_instructions += " path/to/input.txt output.txt";

  bool has_valid_args = false;

  // make sure the user has provided input and output files
  if (argc == 1) {
    cerr << usage_instructions << endl;
  } else if (argc == 2) {
    cerr << "Please include an output file.\n" << usage_instructions << endl;
  } else if (argc == 3) {
    has_valid_args = true;
  } else if (argc > 3) {
    cerr << "Too many arguments.\n" << usage_instructions << endl;
  }

  if (!has_valid_args) {
    exit(EXIT_FAILURE);
  }
}

void check_files(ifstream& in_file, string& in_name,
                 ofstream& out_file, string& out_name) {
  if (!in_file.is_open()) {
    cerr << "Cannot open input file: " << in_name << endl;
    exit(EXIT_FAILURE);
  }

  if (!out_file.is_open()) {
    cerr << "Cannot open output file: " << out_name << endl;
    exit(EXIT_FAILURE);
  }
}

int main(int argc, char* argv[]) {
    string in_file_name_;
    string out_file_name_;
	
	if(!DEBUG){


  check_arguments(argc, argv);

  in_file_name_ = argv[1];
        out_file_name_ = argv[2];
  }
  if(DEBUG){
  in_file_name_ = "../data/sample-laser-radar-measurement-data-1.txt";
   out_file_name_ = "../data/obj_pose-laser-radar-ekf-output-1.txt";
  }

  ifstream in_file_(in_file_name_.c_str(), ifstream::in);
 
  ofstream out_file_(out_file_name_.c_str(), ofstream::out);

  check_files(in_file_, in_file_name_, out_file_, out_file_name_);

  /**********************************************
   *  Set Measurements                          *
   **********************************************/

  vector<MeasurementPackage> measurement_pack_list;
  vector<GroundTruthPackage> gt_pack_list;

  string line;
  size_t number_of_radar_measurements=0;
  size_t number_of_lidar_measurements=0;

  // prep the measurement packages (each line represents a measurement at a
  // timestamp)
  while (getline(in_file_, line)) {
    string sensor_type;
    MeasurementPackage meas_package;
    GroundTruthPackage gt_package;
    istringstream iss(line);
    long long timestamp;

    // reads first element from the current line
    iss >> sensor_type;

    if (sensor_type.compare("L") == 0) {
      // laser measurement

      // read measurements at this timestamp
      meas_package.sensor_type_ = MeasurementPackage::LASER;
      meas_package.raw_measurements_ = VectorXd(2);
      float px;
      float py;
      iss >> px;
      iss >> py;
      meas_package.raw_measurements_ << px, py;
      iss >> timestamp;
      meas_package.timestamp_ = timestamp;
      measurement_pack_list.push_back(meas_package);
	  number_of_lidar_measurements++;
    } else if (sensor_type.compare("R") == 0) {
      // radar measurement

      // read measurements at this timestamp
      meas_package.sensor_type_ = MeasurementPackage::RADAR;
      meas_package.raw_measurements_ = VectorXd(3);
      float ro;
      float phi;
      float ro_dot;
      iss >> ro;
      iss >> phi;
      iss >> ro_dot;
      meas_package.raw_measurements_ << ro, phi, ro_dot;
      iss >> timestamp;
      meas_package.timestamp_ = timestamp;
      measurement_pack_list.push_back(meas_package);
	  number_of_radar_measurements++;
    }

      // read ground truth data to compare later
      float x_gt;
      float y_gt;
      float vx_gt;
      float vy_gt;
      iss >> x_gt;
      iss >> y_gt;
      iss >> vx_gt;
      iss >> vy_gt;
      gt_package.gt_values_ = VectorXd(4);
      gt_package.gt_values_ << x_gt, y_gt, vx_gt, vy_gt;
      gt_pack_list.push_back(gt_package);
  }

  // Create a UKF instance
  UKF ukf;

  // used to compute the RMSE later
  vector<VectorXd> estimations;
  vector<VectorXd> ground_truth;
    VectorXd ukf_x_cartesian_ = VectorXd(4);
    double x_estimate_;
    double y_estimate_;
    double vx_estimate_ ;
    double vy_estimate_ ;




    Tools tools;
  VectorXd rmse(4);

  double NIS_radar;
  double NIS_lidar;

  // start filtering from the second frame (the speed is unknown in the first
  // frame)

  size_t number_of_measurements = measurement_pack_list.size();

  size_t NIS_radar_outcount=0;
  size_t NIS_lidar_outcount=0;


        if(TUNING_MODE){
      // column names for output file
      double J_cost,J_max=0;
      int index=0;

      VectorXd rubric(6),select(2);
      if(DEBUG && in_file_name_ == "../data/sample-laser-radar-measurement-data-1.txt")
          rubric<<0.09,0.09,0.65,0.65,0.7,0.8;
      else if(DEBUG && in_file_name_ == "../data/sample-laser-radar-measurement-data-2.txt")
          rubric<<0.20,0.20,0.55,0.55,0.8,0.8;
      else
          rubric<<0.20,0.20,0.65,0.65,0.8,0.8;

      //out_file_ << "index" << "\t";
      out_file_ << "std_yawdd" << "\t";
      out_file_ << "std_a" << "\t";
      out_file_ << "rmse(0)" << "\t";
      out_file_ << "rmse(1)" << "\t";
      out_file_ << "rmse(2)" << "\t";
      out_file_ << "rmse(3)" << "\t";
      out_file_ << "NIS_radar" << "\t";
      out_file_ << "NIS_lidar" << "\n";
      //out_file_ << "J_cost" << "\n";


    select<<0.4,1.0;

    for(ukf.std_yawdd_=0.0;ukf.std_yawdd_<=1.2;ukf.std_yawdd_+=0.1){

      for(ukf.std_a_=0.0;ukf.std_a_<=6.0;ukf.std_a_+=0.1){
        NIS_radar_outcount=0;
        NIS_lidar_outcount=0;

        for (size_t k = 0; k < number_of_measurements; ++k) {
          ukf.ProcessMeasurement(measurement_pack_list[k]);

            x_estimate_ = ukf.x_(0);
            y_estimate_ = ukf.x_(1);
            vx_estimate_ = ukf.x_(2) * cos(ukf.x_(3));
            vy_estimate_ = ukf.x_(2) * sin(ukf.x_(3));

            ukf_x_cartesian_ << x_estimate_, y_estimate_, vx_estimate_, vy_estimate_;

            estimations.push_back(ukf_x_cartesian_);
            ground_truth.push_back(gt_pack_list[k].gt_values_);


            if (measurement_pack_list[k].sensor_type_ == MeasurementPackage::LASER && (ukf.NIS_laser_>6.0 || ukf.NIS_laser_<0.1))
                 NIS_lidar_outcount++;
            else if (measurement_pack_list[k].sensor_type_ == MeasurementPackage::RADAR && (ukf.NIS_radar_>7.81 || ukf.NIS_radar_<0.35))
                NIS_radar_outcount++;

        }
          NIS_radar=1.0-(NIS_radar_outcount+0.0)/number_of_radar_measurements;
          NIS_lidar=1.0-(NIS_lidar_outcount+0.0)/number_of_lidar_measurements;

          rmse=tools.CalculateRMSE(estimations, ground_truth);


          if(rmse(0)<rubric(0) && rmse(1)<rubric(1) && rmse(2)<rubric(2) && rmse(3)<rubric(3) && NIS_radar>rubric(4) && NIS_lidar>rubric(5)) {
          //if(rmse(0)<rubric(0) && rmse(1)<rubric(1) && rmse(2)<rubric(2) && rmse(3)<rubric(3) ) {
              //J_cost = 0.9 * (rubric(0)-rmse(0)+rubric(1)-rmse(1)+rubric(2)-rmse(2)+rubric(3)-rmse(3))+0.1*(-rubric(4)+NIS_radar -rubric(5)+NIS_lidar);
              J_cost = rubric(0)-rmse(0)+rubric(1)-rmse(1)+rubric(2)-rmse(2)+rubric(3)-rmse(3);

              out_file_ << ukf.std_yawdd_ << "\t";
              out_file_ << ukf.std_a_ << "\t";
              out_file_ << rmse(0) << "\t";
              out_file_ << rmse(1) << "\t";
              out_file_ << rmse(2) << "\t";
              out_file_ << rmse(3) << "\t";
              out_file_ << NIS_radar << "\t";
              out_file_ << NIS_lidar << "\n";
              //out_file_ << J_cost << "\n";

              if(J_cost>J_max){
                  J_max=J_cost;
                  select(0)= ukf.std_yawdd_;
                  select(1)= ukf.std_a_;

              }
          }

          estimations.clear();
          ground_truth.clear();
          //cout<<"the  estimations.size() is "<<estimations.size()<<endl;
          index++;
          cout<<"index "<<index<<"\t std_yawdd "<<ukf.std_yawdd_<<"\t std_a_ "<<ukf.std_a_ <<"\t NIS_radar "<<NIS_radar<<"\t NIS_lidar "<<NIS_lidar
              <<"\t rmse(3) "<<rmse(3)<<"\t rmse(2) "<<rmse(2)<<"\t rmse(1) "<<rmse(1)<<"\t rmse(0) "<<rmse(0)<<endl;

          ukf.is_initialized_=false;

      }

    }
      ukf.std_yawdd_=(double)select(0);
      ukf.std_a_=(double)select(1);
      cout<<"the selected std_yawdd_ is "<<(double)select(0)<<endl;
      cout<<"the selected std_a_ is "<<(double)select(1)<<endl;
  }
    else{
    // column names for output file
    out_file_ << "------------------------------------------------------------------------------------------------------------------" << "\n";
    out_file_ << "px" << "\t";
    out_file_ << "py" << "\t";
    out_file_ << "v" << "\t";
    out_file_ << "yaw_angle" << "\t";
    out_file_ << "yaw_rate" << "\t";
    out_file_ << "px_measured" << "\t";
    out_file_ << "py_measured" << "\t";
    out_file_ << "px_true" << "\t";
    out_file_ << "py_true" << "\t";
    out_file_ << "vx_true" << "\t";
    out_file_ << "vy_true" << "\t";
    out_file_ << "NIS" << "\n";

    NIS_radar_outcount=0;
    NIS_lidar_outcount=0;
    ukf.is_initialized_=false;

  for (size_t k = 0; k < number_of_measurements; ++k) {
    // Call the UKF-based fusion
      //long long measure_time=measurement_pack_list[k].timestamp_;
      //long long  previous_time=ukf.previous_timestamp_;
    ukf.ProcessMeasurement(measurement_pack_list[k]);

    // output the estimation
    out_file_ << ukf.x_(0) << "\t"; // pos1 - est
    out_file_ << ukf.x_(1) << "\t"; // pos2 - est
    out_file_ << ukf.x_(2) << "\t"; // vel_abs -est
    out_file_ << ukf.x_(3) << "\t"; // yaw_angle -est
    out_file_ << ukf.x_(4) << "\t"; // yaw_rate -est

    // output the measurements
    if (measurement_pack_list[k].sensor_type_ == MeasurementPackage::LASER) {
      // output the estimation

      // p1 - meas
      out_file_ << measurement_pack_list[k].raw_measurements_(0) << "\t";

      // p2 - meas
      out_file_ << measurement_pack_list[k].raw_measurements_(1) << "\t";
    } else if (measurement_pack_list[k].sensor_type_ == MeasurementPackage::RADAR) {
      // output the estimation in the cartesian coordinates
      float ro = measurement_pack_list[k].raw_measurements_(0);
      float phi = measurement_pack_list[k].raw_measurements_(1);
      out_file_ << ro * cos(phi) << "\t"; // p1_meas
      out_file_ << ro * sin(phi) << "\t"; // p2_meas
    }

    // output the ground truth packages
    out_file_ << gt_pack_list[k].gt_values_(0) << "\t";
    out_file_ << gt_pack_list[k].gt_values_(1) << "\t";
    out_file_ << gt_pack_list[k].gt_values_(2) << "\t";
    out_file_ << gt_pack_list[k].gt_values_(3) << "\t";

    // output the NIS values
    
    if (measurement_pack_list[k].sensor_type_ == MeasurementPackage::LASER) {
      out_file_ << ukf.NIS_laser_ << "\n";
	  
	  if(ukf.NIS_laser_>6.0 || ukf.NIS_laser_<0.1) NIS_lidar_outcount++;
	  
    } else if (measurement_pack_list[k].sensor_type_ == MeasurementPackage::RADAR) {
      out_file_ << ukf.NIS_radar_ << "\n";
	  
	  if(ukf.NIS_radar_>7.81 || ukf.NIS_radar_<0.35) NIS_radar_outcount++;
    }


    // convert ukf x vector to cartesian to compare to ground truth


     x_estimate_ = ukf.x_(0);
     y_estimate_ = ukf.x_(1);
     vx_estimate_ = ukf.x_(2) * cos(ukf.x_(3));
     vy_estimate_ = ukf.x_(2) * sin(ukf.x_(3));
    
    ukf_x_cartesian_ << x_estimate_, y_estimate_, vx_estimate_, vy_estimate_;
    
    estimations.push_back(ukf_x_cartesian_);
    ground_truth.push_back(gt_pack_list[k].gt_values_);

  }

   NIS_radar=1.0-(NIS_radar_outcount+0.0)/number_of_radar_measurements;
   NIS_lidar=1.0-(NIS_lidar_outcount+0.0)/number_of_lidar_measurements;

  // compute the accuracy (RMSE)

  rmse=tools.CalculateRMSE(estimations, ground_truth);

  cout << "Accuracy - RMSE:" << endl << rmse << endl;

if(DEBUG) {
  cout << "NIS_radar in range rate is " << NIS_radar << "\t with NIS_radar_outlier numbr is: " << NIS_radar_outcount
       << "\t over the total radar records " << number_of_radar_measurements << endl;

  cout << "NIS_lidar in range rate is " << NIS_lidar << "\t with NIS_lidar_outcount_ is: " << NIS_lidar_outcount
       << "\t over the total lidar records " << number_of_lidar_measurements << endl;

  //if(number_of_measurements !=number_of_radar_measurements+number_of_lidar_measurements)
  cout << "number_of_measurements==" << number_of_measurements << "\t number_of_radar_measurements=="
       << number_of_radar_measurements << "\t number_of_lidar_measurements==" << number_of_lidar_measurements << endl;
}
    // close files
  if (out_file_.is_open()) {
    out_file_.close();
  }

  if (in_file_.is_open()) {
    in_file_.close();
  }

  cout << "Done!" << endl;
  return 0;
}
}
