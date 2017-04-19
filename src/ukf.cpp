#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;
#define DEBUG 0

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
    // if this is false, laser measurements will be ignored (except during init)
    use_laser_ = true;

    // if this is false, radar measurements will be ignored (except during init)
    use_radar_ = true;

    n_x_ = 5;
    lambda_ = 3 - n_x_;
    n_aug_ = n_x_ + 2;

    // initial state vector
    x_ = VectorXd(n_x_);

    // initial covariance matrix
    P_ = MatrixXd(n_x_, n_x_);

    P_ <<   1, 0, 0, 0, 0,
            0, 1, 0, 0, 0,
            0, 0, 1, 0, 0,
            0, 0, 0, 1, 0,
            0, 0, 0, 0, 1;


    // Process noise standard deviation longitudinal acceleration in m/s^2
    std_a_ = 4.0;

    // Process noise standard deviation yaw acceleration in rad/s^2
    std_yawdd_ = 0.4;

    // Laser measurement noise standard deviation position1 in m
    std_laspx_ = 0.15;

    // Laser measurement noise standard deviation position2 in m
    std_laspy_ = 0.15;

    // Radar measurement noise standard deviation radius in m
    std_radr_ = 0.3;

    // Radar measurement noise standard deviation angle in rad
    std_radphi_ = 0.003;

    // Radar measurement noise standard deviation radius change in m/s
    std_radrd_ = 0.3;

    Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
    weights_ = VectorXd(2 * n_aug_ + 1);
    // set weights
    weights_(0) = lambda_ / (lambda_ + n_aug_);
    for (int i = 1; i < 2 * n_aug_ + 1; i++) {  //2n+1 weights
        weights_(i) = 0.5 / (n_aug_ + lambda_);;
    }

    R_Radar_ = MatrixXd(3, 3);
    R_Radar_ << std_radr_ * std_radr_, 0, 0,
            0, std_radphi_ * std_radphi_, 0,
            0, 0, std_radrd_ * std_radrd_;

    R_Laser_ = MatrixXd(2, 2);
    R_Laser_ << std_laspx_ * std_laspx_, 0,
            0, std_laspy_ * std_laspy_;

    is_initialized_ = false;
}

UKF::~UKF() {}

void UKF::init(MeasurementPackage meas_package) {
    //this function is to initialize x and P state

    //initialze the previous_timestamp
    previous_timestamp_ = meas_package.timestamp_;

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
        //Convert radar from polar to cartesian coordinates and initialize state.

        double rho = meas_package.raw_measurements_[0];
        double phi = meas_package.raw_measurements_[1];
        double rho_dot = meas_package.raw_measurements_[2];

        x_ << rho * cos(phi), rho * sin(phi), rho_dot, 0, 0;

    } else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
        // Laser coordinages

        double px = meas_package.raw_measurements_[0];
        double py = meas_package.raw_measurements_[1];

        if(px==0 & py ==0)
        {
            px = 0.001;
            py = 0.001;
        }
        x_ << px, py, 0, 0, 0;
    }

    is_initialized_ = true;
}

/**
 * The predict and update  session for each measurement from  either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
    if ((!use_laser_ )&& (meas_package.sensor_type_ == MeasurementPackage::LASER)) return;
    if ((!use_radar_ )&& (meas_package.sensor_type_ == MeasurementPackage::RADAR)) return;

    if (!is_initialized_) {
        init(meas_package);
        return;
    }

    //compute the time elapsed between the current and previous measurements
    double delta_t = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0;    //dt - expressed in seconds
    previous_timestamp_ = meas_package.timestamp_;
    /*
    while(delta_t>0.1) {
        const double dt=0.05;
        Prediction(dt);
        delta_t-=dt;
    }
     */
    Prediction(delta_t);

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
        UpdateRadar(meas_package);
    } else {
        UpdateLidar(meas_package);
    }

    // print the output
    if(DEBUG) {
        cout << "NIS_radar_ = " << NIS_radar_ << endl;
        cout << "NIS_laser_ = " << NIS_laser_ << endl;
        cout << endl;
    }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * no considering of either laser or radar measurement
 */
void UKF::Prediction(double delta_t) {
    MatrixXd Xsig_aug = GenerateSigmaPoints();
    PredictSigmaPoint(delta_t, Xsig_aug);

    //predicted state mean
    x_.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
        x_ = x_ + weights_(i) * Xsig_pred_.col(i);
    }

    //predicted state covariance matrix
    P_.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

        // state difference
        VectorXd x_diff = Xsig_pred_.col(i) - x_;
        //angle normalization
        while (x_diff(3) > M_PI) x_diff(3) -= 2. * M_PI;
        while (x_diff(3) < -M_PI) x_diff(3) += 2. * M_PI;

        P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
    }
}

MatrixXd UKF::GenerateSigmaPoints() {
    // generate the (n_aug_, 2 * n_aug_ + 1) dimension sigma points

    MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

    //create augmented mean vector
    VectorXd x_aug = VectorXd(n_aug_);

    //create augmented state covariance
    MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

    //create augmented mean state
    x_aug.head(n_x_) = x_;
    x_aug(n_x_) = 0;
    x_aug(n_x_ + 1) = 0;

    //create augmented covariance matrix
    P_aug.fill(0.0);
    P_aug.topLeftCorner(n_x_, n_x_) = P_;
    P_aug(n_x_, n_x_) = std_a_ * std_a_;
    P_aug(n_x_ + 1, n_x_ + 1) = std_yawdd_ * std_yawdd_;

    //create square root matrix
    MatrixXd L = P_aug.llt().matrixL();

    //create augmented sigma points
    Xsig_aug.col(0) = x_aug;
    for (int i = 0; i < n_aug_; i++) {
        Xsig_aug.col(i + 1) = x_aug + sqrt(lambda_ + n_aug_) * L.col(i);
        Xsig_aug.col(i + 1 + n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * L.col(i);
    }

    return Xsig_aug;
}

void UKF::PredictSigmaPoint(double delta_t, MatrixXd Xsig_aug) {
    //predict the x sigma points state according to the noise-augmented prvious x sigma points state and the eclapsed time delat_t
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {
        //avoid division by zero
        if (fabs(Xsig_aug(4, i)) > 0.001) {
            Xsig_pred_(0, i)  = Xsig_aug(0, i) + Xsig_aug(2, i) / Xsig_aug(4, i) * (sin(Xsig_aug(3, i) + Xsig_aug(4, i) * delta_t) - sin(Xsig_aug(3, i)));
            Xsig_pred_(1, i)  = Xsig_aug(1, i) + Xsig_aug(2, i) / Xsig_aug(4, i) * (cos(Xsig_aug(3, i)) - cos(Xsig_aug(3, i) + Xsig_aug(4, i) * delta_t));
        } else {
            Xsig_pred_(0, i)  = Xsig_aug(0, i) + Xsig_aug(2, i) * delta_t * cos(Xsig_aug(3, i));
            Xsig_pred_(1, i)  = Xsig_aug(1, i) + Xsig_aug(2, i) * delta_t * sin(Xsig_aug(3, i));
        }

        Xsig_pred_(2, i) = Xsig_aug(2, i);
        Xsig_pred_(3, i)= Xsig_aug(3, i)+ Xsig_aug(4, i) * delta_t;
        Xsig_pred_(4, i) = Xsig_aug(4, i);

        //add noise,calulate output sigma point prediction
        Xsig_pred_(0, i) +=  0.5 * Xsig_aug(5, i) * delta_t * delta_t * cos(Xsig_aug(3, i));
        Xsig_pred_(1, i) +=  0.5 * Xsig_aug(5, i) * delta_t * delta_t * sin(Xsig_aug(3, i));
        Xsig_pred_(2, i) +=  Xsig_aug(5, i) * delta_t;
        Xsig_pred_(3, i) +=  0.5 * Xsig_aug(6, i) * delta_t * delta_t;
        Xsig_pred_(4, i) +=  Xsig_aug(6, i) * delta_t;

    }
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
    //------------------------------------------------
    // MEASUREMENT SPACE PREDICTION

    //  laser measurement dimension
    int n_z = 2;

    //create matrix for sigma points in measurement space
    MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {
        Zsig(0, i) = Xsig_pred_(0, i);  //px
        Zsig(1, i) = Xsig_pred_(1, i);  //py
    }

    //mean predicted measurement
    VectorXd z_pred = VectorXd(n_z);
    z_pred.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {
        z_pred = z_pred + weights_(i) * Zsig.col(i);
    }

    //measurement covariance matrix S
    MatrixXd S = MatrixXd(n_z, n_z);
    S.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {
        //residual
        VectorXd z_diff = Zsig.col(i) - z_pred;

        //angle normalization
        while (z_diff(1) > M_PI) z_diff(1) -= 2. * M_PI;
        while (z_diff(1) < -M_PI) z_diff(1) += 2. * M_PI;

        S +=  weights_(i) * z_diff * z_diff.transpose();
    }

    //add measurement noise covariance matrix
    S +=  R_Laser_;

    //------------------------------------------------
    // UPDATE STATE

    //create matrix for cross correlation Tc
    MatrixXd Tc = MatrixXd(n_x_, n_z);

    //calculate cross correlation matrix
    Tc.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {

        //residual
        VectorXd z_diff = Zsig.col(i) - z_pred;
        //angle normalization
        while (z_diff(1) > M_PI) z_diff(1) -= 2. * M_PI;
        while (z_diff(1) < -M_PI) z_diff(1) += 2. * M_PI;

        // state difference
        VectorXd x_diff = Xsig_pred_.col(i) - x_;
        //angle normalization
        while (x_diff(3) > M_PI) x_diff(3) -= 2. * M_PI;
        while (x_diff(3) < -M_PI) x_diff(3) += 2. * M_PI;

        Tc += weights_(i) * x_diff * z_diff.transpose();
    }

    //Kalman gain K;
    MatrixXd K = Tc * S.inverse();

    VectorXd z = VectorXd(n_z);
    z << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1];

    //residual
    VectorXd z_diff = z - z_pred;

    //angle normalization
    while (z_diff(1) > M_PI) z_diff(1) -= 2. * M_PI;
    while (z_diff(1) < -M_PI) z_diff(1) += 2. * M_PI;

    //update state mean and covariance matrix
    x_ = x_ + K * z_diff;
    P_ = P_ - K * S * K.transpose();

    // update NIS: Norzmalized Innovation Squared
    NIS_laser_ = (meas_package.raw_measurements_ - z_pred).transpose() * S.inverse() *
                 (meas_package.raw_measurements_ - z_pred);
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
    //------------------------------------------------------------------------------------------------
    // MEASUREMENT SPACE PREDICTION

    int n_z = 3;

    //create matrix for sigma points in measurement space
    MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

    //transform sigma points into measurement space
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {

        // p_x   Xsig_pred_(0, i)
        // p_y  Xsig_pred_(1, i)
        // v    Xsig_pred_(2, i)
        // yaw  Xsig_pred_(03, i)

        // measurement model
        Zsig(0, i) = sqrt(Xsig_pred_(0, i) * Xsig_pred_(0, i) + Xsig_pred_(1, i) * Xsig_pred_(1, i));  //r
        Zsig(1, i) = atan2(Xsig_pred_(1, i), Xsig_pred_(0, i));                           //phi
        if ((fabs(Zsig(0, i)) > 0.001))
             Zsig(2, i) = (Xsig_pred_(0, i) * cos(Xsig_pred_(3, i)) * Xsig_pred_(2, i) + Xsig_pred_(1, i) * sin(Xsig_pred_(3, i)) * Xsig_pred_(2, i)) / Zsig(0, i);   //r_dot
        else Zsig(2, i) = 0;
    }

    //mean predicted measurement
    VectorXd z_pred = VectorXd(n_z);
    z_pred.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {
        z_pred = z_pred + weights_(i) * Zsig.col(i);
    }

    //measurement covariance matrix S
    MatrixXd S = MatrixXd(n_z, n_z);
    S.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {
        //residual
        VectorXd z_diff = Zsig.col(i) - z_pred;

        //angle normalization
        while (z_diff(1) > M_PI) z_diff(1) -= 2. * M_PI;
        while (z_diff(1) < -M_PI) z_diff(1) += 2. * M_PI;

        S +=  weights_(i) * z_diff * z_diff.transpose();
    }
    //add measurement noise covariance matrix

    S +=  R_Radar_;

    //------------------------------------------------
    // UPDATE STATE

    //create matrix for cross correlation Tc
    MatrixXd Tc = MatrixXd(n_x_, n_z);

    //calculate cross correlation matrix
    Tc.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

        //residual
        VectorXd z_diff = Zsig.col(i) - z_pred;
        //angle normalization
        while (z_diff(1) > M_PI) z_diff(1) -= 2. * M_PI;
        while (z_diff(1) < -M_PI) z_diff(1) += 2. * M_PI;

        // state difference
        VectorXd x_diff = Xsig_pred_.col(i) - x_;
        //angle normalization
        while (x_diff(3) > M_PI) x_diff(3) -= 2. * M_PI;
        while (x_diff(3) < -M_PI) x_diff(3) += 2. * M_PI;

        Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
    }

    //Kalman gain K;

    MatrixXd K = Tc * S.inverse();
    VectorXd z = VectorXd(n_z);

    z << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], meas_package.raw_measurements_[2];

    //residual
    VectorXd z_diff = z - z_pred;
    //angle normalization
    while (z_diff(1) > M_PI) z_diff(1) -= 2. * M_PI;
    while (z_diff(1) < -M_PI) z_diff(1) += 2. * M_PI;

    //update state mean and covariance matrix
    x_ = x_ + K * z_diff;
    P_ = P_ - K * S * K.transpose();

    // update NIS: Norzmalized Innovation Squared
    NIS_radar_ = (meas_package.raw_measurements_ - z_pred).transpose() * S.inverse() *
                 (meas_package.raw_measurements_ - z_pred);
}