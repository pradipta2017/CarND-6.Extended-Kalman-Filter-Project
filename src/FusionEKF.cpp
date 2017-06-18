#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);
  Pfusion_ = MatrixXd(4, 4);
  Ffusion_ = MatrixXd(4, 4);
  Qfusion_ = MatrixXd(4, 4);
  xfusion_ = VectorXd(4);


  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */
  H_laser_ << 1, 0, 0, 0,
	  0, 1, 0, 0;

  xfusion_ << 1, 1, 1, 1;
  Pfusion_<< 1, 0, 0, 0,
	  0, 1, 0, 0,
	  0, 0, 100, 0,
	  0, 0, 0, 50;

  Ffusion_ << 1, 0, 0, 0,
	  0, 1, 0, 0,
	  0, 0, 1, 0,
	  0, 0, 0, 1;

  Qfusion_ << 0, 0, 0, 0,
	  0, 0, 0, 0,
	  0, 0, 0, 0,
	  0, 0, 0, 0;

  noise_ax = 9.0;
  noise_ay = 9.0;
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    cout << "EKF: " << endl;
	//kf_.x_ = VectorXd(4);
	//ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
		const float m1 = measurement_pack.raw_measurements_[0] * cos(measurement_pack.raw_measurements_[1]);
		const float m2 = measurement_pack.raw_measurements_[0] * sin(measurement_pack.raw_measurements_[1]);
		const double vx = measurement_pack.raw_measurements_[2] * cos(measurement_pack.raw_measurements_[1]);
		const double vy = measurement_pack.raw_measurements_[2] * sin(measurement_pack.raw_measurements_[1]);
		xfusion_ << m1, m2, vx, vy;
		ekf_.Init(xfusion_, Pfusion_, Ffusion_, Hj_, R_radar_, Qfusion_);
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
		xfusion_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
		ekf_.Init(xfusion_, Pfusion_, Ffusion_, H_laser_, R_laser_, Qfusion_);
    }

	
	previous_timestamp_ = measurement_pack.timestamp_;
    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */

  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0; //dt - expressed in seconds
  if (dt == 0) {
	  return;
  }
  previous_timestamp_ = measurement_pack.timestamp_;

  ekf_.F_ << 1, 0, dt, 0,
	  0, 1, 0, dt,
	  0, 0, 1, 0,
	  0, 0, 0, 1;

  //2. Set the process covariance matrix Q
  ekf_.Q_ << (pow(dt, 4) / 4)* noise_ax, 0, (pow(dt, 3) / 2)* noise_ax, 0,
	  0, (pow(dt, 4) / 4)* noise_ay, 0, (pow(dt, 3) / 2)* noise_ay,
	  (pow(dt, 3) / 2)* noise_ax, 0, pow(dt, 2)* noise_ax, 0,
	  0, (pow(dt, 3) / 2)* noise_ay, 0, pow(dt, 2)* noise_ay;

  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
	ekf_.R_ = R_radar_;
	ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
	ekf_.UpdateEKF(measurement_pack.raw_measurements_);

  } else {
    // Laser updates
    ekf_.H_ = H_laser_;
	ekf_.R_ = R_laser_;

	ekf_.Update(measurement_pack.raw_measurements_);

  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
