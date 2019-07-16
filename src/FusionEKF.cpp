#include "FusionEKF.h"
#include <iostream>
//#include "eigen3/Eigen/Dense"
#include "Eigen/Dense"
#include "tools.h"


using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

Tools tools;
MatrixXd CalculateJacobian(const VectorXd& x_state);


/**
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

	//measurement covariance matrix - laser
	R_laser_ << 0.0225, 0,
			  0, 0.0225;
	H_laser_ << 1,0,0,0,
				0,1,0,0;


	//measurement covariance matrix - radar
	R_radar_ << 0.09, 0, 0,
			  0, 0.0009, 0,
			  0, 0, 0.09;

	Hj_ << 1,0,0,0,
			0,1,0,0,
			0,0,1,0;

	ekf_.F_ = MatrixXd(4, 4);
	ekf_.F_ << 1, 0, 1, 0,
			  0, 1, 0, 1,
			  0, 0, 1, 0,
			  0, 0, 0, 1;

	ekf_.P_ = MatrixXd(4,4);
	ekf_.P_ << 1,0,0,0,
			  0,1,0,0,
			  0,0,1000,0,
			  0,0,0,1000;

	ekf_.Q_ = MatrixXd(4, 4);
	ekf_.P_ << 1,0,0,0,
			  0,1,0,0,
			  0,0,1,0,
			  0,0,0,1;

	noise_ax = 9.0;
	noise_ay = 9.0;
}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /**
   * Initialization
   */
	if (!is_initialized_) {
		/**
		 * Initialize the state ekf_.x_ with the first measurement.
		 * Create the covariance matrix.
		 */
		cout << "EKF: " << endl;
		ekf_.x_ = VectorXd(4);
		ekf_.x_ << 1, 1, 1, 1;

		if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
			//Convert radar from polar to cartesian coordinates 
			//         and initialize state.
			double rho = measurement_pack.raw_measurements_(0);
			double phi =  measurement_pack.raw_measurements_(1);
			double rho_dot =  measurement_pack.raw_measurements_(2);
			double M_x = rho * cos(phi) > 0.000001 ? 0.000001 : 0;
			double M_y = rho * sin(phi) > 0.000001 ? 0.000001 : 0;
			ekf_.x_ << M_x, M_y, 0, 0;  // x, y, vx, vy
		}
		else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
			// Initialize state.
			ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0; // x, y, vx, vy
		}
		// done initializing, no need to predict or update
		is_initialized_ = true;
		previous_timestamp_ = measurement_pack.timestamp_;
		return;
	}

  /**
   * Prediction
   */
	double dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
    previous_timestamp_ = measurement_pack.timestamp_;
    double dt2 = dt * dt;
    double dt3 = dt2 * dt;
    double dt4 = dt3 * dt;
    ekf_.F_(0, 2) = dt;
    ekf_.F_(1, 3) = dt;
    // 2. Set the process covariance matrix Q
    ekf_.Q_ = MatrixXd(4, 4);
    ekf_.Q_ << dt4/4*noise_ax, 0, dt3/2*noise_ax, 0,
            0, dt4/4*noise_ay, 0, dt3/2*noise_ay,
            dt3/2*noise_ax, 0, dt2*noise_ax, 0,
            0, dt3/2*noise_ay, 0, dt2*noise_ay;

	ekf_.Predict();
  /**
   * Update
   */

	if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
		Hj_ = tools.CalculateJacobian(ekf_.x_);
		ekf_.H_ = Hj_;
		ekf_.R_ = R_radar_;
		ekf_.UpdateEKF(measurement_pack.raw_measurements_);

  } else {
		ekf_.H_ = H_laser_;
		ekf_.R_ = R_laser_;
		ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "__________________" << endl;
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "-------------" << endl;
  cout << "P_ = " << ekf_.P_ << endl;
  cout << "__________________" << endl;

}
