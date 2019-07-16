#include "kalman_filter.h"
#include "tools.h"
using Eigen::MatrixXd;
using Eigen::VectorXd;
/* 
 * Please note that the Eigen library does not initialize 
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}




void KalmanFilter::Predict() {
	// predict the state
	x_ = F_ * x_;
	MatrixXd Ft = F_.transpose();
	P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
	// update the state by using Kalman Filter equations 
	VectorXd z_pred = H_ * x_;
	VectorXd y = z - z_pred;
	MatrixXd Ht = H_.transpose();
	MatrixXd S = H_ * P_ * Ht + R_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si;

	//new estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;

}



VectorXd CalculateCartesian2Polar(const VectorXd &x_state) {
	
	double px, py, vx, vy;
	px = x_state[0];
	py = x_state[1];
	vx = x_state[2];
	vy = x_state[3];

	double rho, phi, rho_dot;
	rho = sqrt(px*px + py*py);
	phi = atan2(py, px);

	if (fabs(rho) < 0.000001) {
		//cout << "CalculateCartasian2Polar - Error - Division by zero" << endl;
		rho_dot = 0.000001;
	}

	rho_dot = (px * vx + py * vy) / rho;
	VectorXd z_pred = VectorXd(3);
	z_pred << rho, phi, rho_dot;
	return z_pred;

}


void KalmanFilter::UpdateEKF(const VectorXd &z) {
	// update the state by using Extended Kalman Filter equations
	VectorXd z_pred = CalculateCartesian2Polar(x_);
	VectorXd y = z - z_pred;

	// normalize the angle between -pi to pi
	while(y(1) > M_PI){
		y(1) -= 2 * M_PI;
	}

	while(y(1) < -M_PI){
		y(1) += 2 * M_PI;
	}

	// following is exact the same as in the function of KalmanFilter::Update()
	MatrixXd Ht = H_.transpose();
	MatrixXd PHt = P_ * Ht;
	MatrixXd S = H_ * PHt + R_;
	MatrixXd Si = S.inverse();
	MatrixXd K = PHt * Si;

	//new estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;
}
