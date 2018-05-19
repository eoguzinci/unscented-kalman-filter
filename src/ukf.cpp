#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
#include "measurement_package.h"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  /* PROCESS NOISE */
  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 2;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.3;
  
  /* MEASUREMENT NOISE */
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
  
  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  // at the beginning, not initialized
  is_initialized_ = false;

  // time is 0.0 at the beginning
  time_us_ = 0.0;

  // State dimension is 5 for CTRV model [px py v psi dpsi].trasnpose()
  n_x_ = 5;

  //Augmented state is concat(state,[n_ddx n_ddpsi])
  n_aug_ = 7;

  // width of sigma points
  lambda_ = 3- n_x_;

  // predicted sigma points
  Xsig_pred_ = MatrixXd(n_x_ , 2*n_aug_+1);

  // weights for sigma points
  weights_ = VectorXd(2*n_aug_+1);

}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  // check first if any sensor type is ignored, if yes do not count it
  if((meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) || (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_)){
    // INITIALIZE
    if(!is_initialized_)
    {
      x_.fill(0.0);
      P_.diagonal() << 1,1,1,1,1;

      time_us_ = meas_package.timestamp_;

      if(meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_) {
        x_(0) = meas_package.raw_measurements_(0);
        x_(1) = meas_package.raw_measurements_(1);
      }
      else if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_){
        float ro = meas_package.raw_measurements_(0);
        float phi = meas_package.raw_measurements_(1);
        x_(0) = ro * cos(phi);
        x_(1) = ro * sin(phi);
      }

      weights_(0) = lambda_ / (lambda_+n_aug_);
      for (int i = 1; i < 2*n_aug_+1; ++i) {
        weights_(i) = 0.5/(lambda_+n_aug_);
      }

      // So, it is initialized and won't go throught this part again.
      is_initialized_ = true;
      
      return;
    }


    //compute the time elapsed between the current and previous measurements
    float dt = (meas_package.timestamp_ - time_us_) / 1000000.0; //dt - expressed in seconds
    time_us_ = meas_package.timestamp_;

    Prediction(dt);

    if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      UpdateLidar(meas_package);
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      UpdateRadar(meas_package);
    }
  }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

  /* Generation of Sigma Points */

  //define spreading parameter
  double lambda = 3 - n_aug_;

  //create augmented mean vector
  VectorXd x_aug = VectorXd(7);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(7, 7);

  //create augmented sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  //create augmented mean state
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  cout << "PREDICTION STEP" << endl;
  cout<< "-------------" <<endl;

  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5,5) = std_a_*std_a_;
  P_aug(6,6) = std_yawdd_*std_yawdd_;

  //create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  //create augmented sigma points
  Xsig_aug.col(0)  = x_aug;
  for (int i = 0; i< n_aug_; i++){
    Xsig_aug.col(i+1)       = x_aug + sqrt(lambda+n_aug_) * L.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda+n_aug_) * L.col(i);
  }

  //create matrix with predicted sigma points as columns
  MatrixXd Xsig_pred = MatrixXd(n_x_, 2 * n_aug_ + 1);
  VectorXd Vsig(n_x_); VectorXd Asig(n_x_);

  // predict sigma points
  for (int i = 0; i < 2 * n_aug_ +1; ++i) {
    Xsig_pred.col(i) = Xsig_aug.col(i).head(n_x_);
    Vsig.fill(0.0);
    if (Xsig_aug(4)<0.0001){
      Vsig(0) = Xsig_aug(2,i)*cos(Xsig_aug(3,i))*delta_t;
      Vsig(1) = Xsig_aug(2,i)*sin(Xsig_aug(3,i))*delta_t;
    }
    else{
      Vsig(0) = Xsig_aug(2,i)/Xsig_aug(4,i)*(sin(Xsig_aug(3,i)+Xsig_aug(4,i)*delta_t)-sin(Xsig_aug(3,i)));
      Vsig(1) = Xsig_aug(2,i)/Xsig_aug(4,i)*(-cos(Xsig_aug(3,i)+Xsig_aug(4,i)*delta_t)+cos(Xsig_aug(3,i)));
      Vsig(3) = Xsig_aug(4,i)*delta_t;
    }
    Xsig_pred.col(i) += Vsig;

    Asig(0) = delta_t*delta_t*cos(Xsig_aug(3,i))*Xsig_aug(5,i)/2;
    Asig(1) = delta_t*delta_t*sin(Xsig_aug(3,i))*Xsig_aug(5,i)/2;
    Asig(2) = Xsig_aug(5,i)*delta_t;
    Asig(3) = delta_t*delta_t*Xsig_aug(6,i)/2;
    Asig(4) = delta_t*Xsig_aug(6,i);

    Xsig_pred.col(i) += Asig;
  }

  // predicted state mean
  x_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    x_ += weights_(i) * Xsig_pred.col(i);
  }

  //predicted state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
    // state difference
    VectorXd x_diff = Xsig_pred.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose() ;
  }

}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */

  // store the lidar measurements
  VectorXd z = meas_package.raw_measurements_;

  //set measurement dimension, radar can measure px and py
  int n_z_ = 2;

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z_, 2 * n_aug_ + 1);

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z_);

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z_,n_z_);

  // radar measurement transformation
  for (int j = 0; j < 2*n_aug_+1; ++j) {
    double px = Xsig_pred_(0,j);
    double py = Xsig_pred_(1,j);

    Zsig(0,j) = px;
    Zsig(1,j) = py;
  }

  //mean predicted measurement
  z_pred.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
      z_pred += weights_(i) * Zsig.col(i);
  }

  S.fill(0.0);
  MatrixXd Z_diff(n_z_,2*n_aug_+1);
  for (int i = 0; i < 2*n_aug_+1; ++i) {
    VectorXd z_diff = Zsig.col(i) - z_pred;
    Z_diff.col(i) = z_diff;
    S = S + weights_(i) * z_diff * z_diff.transpose() ;
  }

  // add measurement noise covariance matrix
  MatrixXd R(n_z_,n_z_);
  R.fill(0.0);
  R(0,0) = std_laspx_*std_laspx_;
  R(1,1) = std_laspy_*std_laspy_;

  S += R;

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z_);

  Tc.fill(0.0);
  for (int i = 0; i < 2*n_aug_+1; ++i) {
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    // covariance difference
    VectorXd z_diff = Z_diff.col(i);
    //angle normalization
    x_diff(3) = atan2(sin(x_diff(3)),cos(x_diff(3)));
    Tc = Tc + weights_(i) * x_diff * z_diff.transpose() ;
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = z - z_pred;

  //calculate NIS
  NIS_laser_ = z_diff.transpose() * S.inverse() * z_diff;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();

  // print the output
  cout << "x_ = " << x_ << endl;
  cout << "P_ = " << P_ << endl;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */

  // store the radar measurements
  VectorXd z = meas_package.raw_measurements_;

  //set measurement dimension, radar can measure r, phi, and r_dot
  int n_z_ = 3;

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z_, 2 * n_aug_ + 1);

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z_);

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z_,n_z_);

  cout << "RADAR UPDATE STEP" << endl;
  cout<< "-------------" <<endl;

  cout << "Xsig_pred_ =" << endl << Xsig_pred_ << endl; 

  // radar measurement transformation
  for (int j = 0; j < 2*n_aug_+1; ++j) {
    double px = Xsig_pred_(0,j);
    double py = Xsig_pred_(1,j);
    double v = Xsig_pred_(2,j);
    double psi = Xsig_pred_(3,j);

    Zsig(0,j) = sqrt(px*px+py*py);
    Zsig(1,j) = atan2(py,px);
    // prevent divide by 0
    if (fabs(Zsig(0,j)) > 1e-4) { 
      Zsig(2,j) = (px*cos(psi)*v + py*sin(psi)*v)/Zsig(0,j);
    }else{
      Zsig(2,j) = 0.0;
    }
    
  }

  cout << "ZSig =" << endl << Zsig << endl;

  //mean predicted measurement
  z_pred.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
      z_pred += weights_(i) * Zsig.col(i);
  }

  cout << "z_pred =" << endl << z_pred << endl;

  S.fill(0.0);
  MatrixXd Z_diff(n_z_,2*n_aug_+1);
  for (int i = 0; i < 2*n_aug_+1; ++i) {
    VectorXd z_diff = Zsig.col(i) - z_pred;
    // while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    // while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
    z_diff(1) = atan2(sin(z_diff(1)),cos(z_diff(1)));
    Z_diff.col(i) = z_diff;
    S = S + weights_(i) * z_diff * z_diff.transpose() ;
  }

  cout << "S = " << endl << S <<endl;

  // add measurement noise covariance matrix
  MatrixXd R(n_z_,n_z_);
  R.fill(0.0);
  R(0,0) = std_radr_*std_radr_;
  R(1,1) = std_radphi_*std_radphi_;
  R(2,2) = std_radrd_*std_radrd_;

  cout << "R =" <<endl << R << endl;

  S += R;

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z_);

  Tc.fill(0.0);
  for (int i = 0; i < 2*n_aug_+1; ++i) {
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    // covariance difference
    VectorXd z_diff = Z_diff.col(i);
    //angle normalization
    // while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    // while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
    // while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    // while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
    x_diff(3) = atan2(sin(x_diff(3)),cos(x_diff(3)));
    z_diff(1) = atan2(sin(z_diff(1)),cos(z_diff(1)));
    Tc = Tc + weights_(i) * x_diff * z_diff.transpose() ;
  }

  cout << "Tc" << endl << Tc << endl;

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  cout << "K" << endl << K << endl;

  //residual
  VectorXd z_diff = z - z_pred;

  //angle normalization
  // while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  // while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
  z_diff(1) = atan2(sin(z_diff(1)),cos(z_diff(1)));

  cout << "z_diff =" << endl << z_diff << endl;

  //calculate NIS
  NIS_radar_ = z_diff.transpose() * S.inverse() * z_diff;

  cout << "NIS_radar" <<endl<< NIS_radar_ << endl;

  // state update
  x_ = x_ + K*z_diff;

  // covariance update
  P_ = P_ - K*S*K.transpose();

  // print the output
  cout << "x_ = " << x_ << endl;
  cout << "P_ = " << P_ << endl;
}
