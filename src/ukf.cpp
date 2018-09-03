#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

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
    
    // set is_initialized to false;
    is_initialized_ = false;
    
    // set number of elements in the state vector.
    n_x_ = 5;
    
    // Set vector size for augmented state
    n_aug_ = 7;
    
    // Set vector size radar measurement
    n_radar = 3;
    
    // Set vector size lidar measurement
    n_lidar = 2;
    
    // calculate number of columns for sigma points matrix;
    n_sig_columns_ = n_aug_ * 2 + 1;
    
    // initial state vector
    x_ = VectorXd(n_x_);
    x_.fill(0.0);
    
    // initial covariance matrix
    P_ = MatrixXd(n_x_, n_x_);
    
    // Process noise standard deviation longitudinal acceleration in m/s^2
    std_a_ = 3;
    
    // Process noise standard deviation yaw acceleration in rad/s^2
    std_yawdd_ = 3;
    
    // set lamda value for spread of sigma points.
    lambda_ = 3 - n_aug_;
    
    p_init_ = 1;
    
    // init process co-varience.
    P_ = MatrixXd(n_x_, n_x_);
    P_ << p_init_, 0, 0, 0, 0,
        0, p_init_, 0, 0, 0,
        0, 0, p_init_, 0, 0,
        0, 0, 0, p_init_, 0,
        0, 0, 0, 0, p_init_;
    
    // set weights
    weights_ = VectorXd(n_sig_columns_);
    weights_(0) = lambda_ / (lambda_ + n_aug_);
    weights_.tail(n_sig_columns_ - 1).fill(0.5 / (n_aug_ + lambda_));

    // instantiate state sigma points matrix
    Xsig_pred_ = MatrixXd(n_x_, n_sig_columns_);
    
    // Linear measurement matrix for Lidar udpate
    H_ = MatrixXd(n_lidar, n_x_);
    H_ << 1, 0, 0, 0, 0,
    0, 1, 0, 0, 0;
    
    // Linear measurement covariance matrix for lidar
    R_laser_ = MatrixXd(2, 2);
    R_laser_ << std_laspx_*std_laspx_, 0,
        0, std_laspy*std_laspy;
    
    // Calculate radar noise matrix
    R_radar_ = MatrixXd(n_radar, n_radar);
    R_radar_ << std_radr_*std_radr_, 0, 0,
    0, std_radphi_*std_radphi_, 0,
    0, 0, std_radrd_*std_radrd_;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
    if (!IsMeasurementEnabled(meas_package.sensor_type_)) {
        cout << "skipping measurement with sensor type " << meas_package.sensor_type_ << endl;
        return;
    }
    
    if (!is_initialized_) {
        FirstMeasurement(meas_package);
        is_initialized_ = true;
        cout << "processed first message" << endl;
        return;
    }
    
    // Compute the time elapsed between the current and previous measurements in seconds
    double delta_t = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0;
    previous_timestamp_ = meas_package.timestamp_;
    
    Prediction(delta_t);
    
    if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
        UpdateLidar(meas_package);
    } else if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
        UpdateRadar(meas_package);
    }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
    // Create augmented mean vector to factor in noise
    VectorXd x_aug = VectorXd(n_aug_);
    x_aug.head(5) = x_;
    x_aug.tail(2).fill(0);
    
    // Create augmented covarience matrix
    MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
    P_aug.fill(0.0);
    P_aug.topLeftCorner(n_x_, n_x_) = P_;
    P_aug(5, 5) = std_a_*std_a_;
    P_aug(6,6) = std_yawdd_*std_yawdd_;
    
    // Create square root matrix
    MatrixXd L = P_aug.llt().matrixL();
    
    // Extract mean sigma points
    MatrixXd Xsig_aug = MatrixXd(n_aug_, n_sig_columns_);
    Xsig_aug.col(0) = x_aug;
    for (int i = 0; i< n_aug_; i++)
    {
        Xsig_aug.col(i+1)       = x_aug + sqrt(lambda_ + n_aug_) * L.col(i);
        Xsig_aug.col(i + 1 + n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * L.col(i);
    }
    
    PredictSigmaPoints(Xsig_aug, delta_t);
    
    // Predict state mean
    x_.fill(0.0);
    for (int i = 0; i < n_sig_columns_; i++) {
        x_ += weights_(i) * Xsig_pred_.col(i);
    }
    
    // Predict state covariance matrix
    P_.fill(0.0);
    for (int i = 0; i < n_sig_columns_; i++) {
        VectorXd x_diff = Xsig_pred_.col(i) - x_;
        x_diff(3) = NormalizeAngle(x_diff(3));
        
        P_ += weights_(i) * x_diff * x_diff.transpose();
    }
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
    VectorXd z_pred = H_ * x_;
    VectorXd y = meas_package.raw_measurements_ - z_pred;
    MatrixXd Ht = H_.transpose();
    MatrixXd S = H_ * P_ * Ht + R_laser_;
    MatrixXd Si = S.inverse();
    MatrixXd PHt = P_ * Ht;
    MatrixXd K = PHt * Si;
    
    // new estimate
    x_ = x_ + (K * y);
    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K * H_) * P_;
    
    // Calulate NIS
    VectorXd z_diff = z_pred - meas_package.raw_measurements_;
    float nis = z_diff.transpose() * S.inverse() * z_diff;
    
    cout << "NIS value for LIDAR: " << nis << endl;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
    //transform sigma points into measurement space
    //create matrix for sigma points in measurement space
    MatrixXd Zsig(n_radar, n_sig_columns_);
    for (int i = 0; i < n_sig_columns_; i++) {
        double px = Xsig_pred_(0, i);
        double py = Xsig_pred_(1, i);
        double v = Xsig_pred_(2, i);
        double yaw = Xsig_pred_(3, i);
        
        double rho = sqrt(px*px + py*py);
        
        Zsig(0, i) = rho;
        Zsig(1, i) = atan2(py, px);
        Zsig(2, i) = (px*cos(yaw)*v + py*sin(yaw)*v) / rho;
    }
    
    //mean predicted measurement
    VectorXd z_pred(n_radar);
    z_pred.fill(0.0);
    // Calculate mean predicted measurement
    for (int i = 0; i < n_sig_columns_; i++) {
        z_pred += weights_(i) * Zsig.col(i);
    }
    
    z_pred(1) = NormalizeAngle(z_pred(1));
    
    // Calculate measurement covariance matrix S
    MatrixXd S = MatrixXd(n_radar, n_radar);
    S.fill(0.0);
    for (int i = 0; i < n_sig_columns_; i++) {
        VectorXd z_diff = Zsig.col(i) - z_pred;
        
        z_diff(1) = NormalizeAngle(z_diff(1));
        
        S += weights_(i) * z_diff * z_diff.transpose();
    }
    
    S += R_radar_;
    
    // Create matrix for cross correlation Tc
    MatrixXd Tc(n_x_, n_radar);
    Tc.fill(0.0);
    
    // Calculate cross correlation matrix
    for (int i = 0; i < n_sig_columns_; i++) {
        VectorXd x_diff = Xsig_pred_.col(i) - x_;
        x_diff(3) = NormalizeAngle(x_diff(3));
        
        VectorXd z_diff = Zsig.col(i) - z_pred;
        z_diff(1) = NormalizeAngle(z_diff(1));
        
        Tc += weights_(i) * x_diff * z_diff.transpose();
    }
    
    // Calculate Kalman gain K;
    MatrixXd K = Tc * S.inverse();
    
    //update state mean and covariance matrix
    VectorXd z(3);
    z << meas_package.raw_measurements_(0), meas_package.raw_measurements_(1), meas_package.raw_measurements_(2);
    
    VectorXd z_diff = z - z_pred;
    z_diff(1) = NormalizeAngle(z_diff(1));
    
    x_ += K * z_diff;
    P_ -= K * S * K.transpose();
    
    // Calulate NIS
    float nis = z_diff.transpose() * S.inverse() * z_diff;
    
    cout << "NIS value for RADAR: " << nis << endl;
}

/**
 * Initializes the UKF on the first measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::FirstMeasurement(MeasurementPackage meas_package) {
    previous_timestamp_ = meas_package.timestamp_;
    if(meas_package.sensor_type_ == MeasurementPackage::LASER) {
        x_(0) = meas_package.raw_measurements_[0];
        x_(1) = meas_package.raw_measurements_[1];
    } else if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
        // Convert radar from polar to cartesian coordinates and initialize state.
        float rho = meas_package.raw_measurements_[0];
        float theta = meas_package.raw_measurements_[1];
        
        x_(0) = rho * cos(theta);
        x_(1) = rho * sin(theta);
    }
}

/**
 * Predicts sigma points from augmented sigma point matrix step k for step k+1
 * @param {MatrixXd} Xsig_aug
 * @param {MatrixXd} Xsig_aug
 */
void UKF::PredictSigmaPoints(MatrixXd& Xsig_aug, double delta_t) {
    // Predict sigma points k + 1
    for (int i = 0; i < n_sig_columns_; i++) {
        //extract values for better readability
        double p_x = Xsig_aug(0,i);
        double p_y = Xsig_aug(1,i);
        double v = Xsig_aug(2,i);
        double yaw = Xsig_aug(3,i);
        double yawd = Xsig_aug(4,i);
        double nu_a = Xsig_aug(5,i);
        double nu_yawdd = Xsig_aug(6,i);
        
        //predicted state values
        double px_p, py_p;
        
        //avoid division by zero
        if (fabs(yawd) > 0.001) {
            px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
            py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
        }
        else {
            px_p = p_x + v*delta_t*cos(yaw);
            py_p = p_y + v*delta_t*sin(yaw);
        }
        
        double v_p = v;
        double yaw_p = yaw + yawd*delta_t;
        double yawd_p = yawd;
        
        //add noise
        px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
        py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
        v_p = v_p + nu_a*delta_t;
        
        yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
        yawd_p = yawd_p + nu_yawdd*delta_t;
        
        //write predicted sigma point into right column
        Xsig_pred_(0,i) = px_p;
        Xsig_pred_(1,i) = py_p;
        Xsig_pred_(2,i) = v_p;
        Xsig_pred_(3,i) = yaw_p;
        Xsig_pred_(4,i) = yawd_p;
    }
}

/**
 * Predicts determines if measurement is supported baesd on implementation support and configuration
 * @returns {bool}
 */
bool UKF::IsMeasurementEnabled(MeasurementPackage::SensorType type) {
    switch (type) {
        case MeasurementPackage::SensorType::LASER:
            return use_laser_;
            break;
        case MeasurementPackage::SensorType::RADAR:
            return use_radar_;
            break;
        default:
            return false;
            break;
    }
}

float UKF::NormalizeAngle(double angle) {
    while (abs(angle) > M_PI) {
        if (angle > M_PI) {
            angle -= 2 * M_PI;
        } else {
            angle += 2 * M_PI;
        }
    }
    
    return angle;
}


