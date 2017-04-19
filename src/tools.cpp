#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;
using namespace std;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
  VectorXd rmse(4);
	rmse << 0,0,0,0;
	//input data validation check
	
  	if(estimations.size()==0){
	    cout<<"Error:the  estimations.size() is zero"<<endl;
	    return rmse;
	}
	if( estimations.size()!= ground_truth.size()) {
	   cout<<"Error:estimations.size()!= ground_truth.size()"<<endl;
	    return rmse;
	}


	// rmse calculate
	
	for(int i=0; i < estimations.size(); ++i){
		VectorXd e=estimations[i]-ground_truth[i];
		 e=e.array()*e.array();
		rmse+=e;
		
	}
	
	
    	//calculate the mean
	rmse = rmse/estimations.size();

	//calculate the squared root
	rmse = rmse.array().sqrt();
	//return the result
	return rmse;
}
