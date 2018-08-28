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
    
    // set number of elements in the augmented state vector to include noise.
    n_aug_ = 7;
    
    // calculate number of columns for sigma points matrix;
    n_sig_columns_ = n_aug_ * 2 + 1;
    
    // initial state vector
    x_ = VectorXd(n_x_);
    x_.fill(0.0);
    
    // initial covariance matrix
    P_ = MatrixXd(5, 5);
    
    // Process noise standard deviation longitudinal acceleration in m/s^2
    std_a_ = 30;
    
    // Process noise standard deviation yaw acceleration in rad/s^2
    std_yawdd_ = 30;
    
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
    weights_(0) = lambda_ / (lambda_ + n_aug_);
    weights_.tail(n_sig_columns_).fill(0.5 / (n_aug_ + lambda_));
    
    // instantiate state sigma points matrix
    Xsig_pred_ = MatrixXd(n_aug_, n_sig_columns_);
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
    if (!is_initialized_) {
        FirstMeasurement(meas_package);
        return;
    }
    
    Prediction(meas_package.timestamp_);
    
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
    for (int i = 0; i < n_sig_columns_; i++) {
        VectorXd x_diff = Xsig_pred_.col(i) - x_;
        P_ += weights_(i) * x_diff * x_diff.transpose();
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
}

/**
 * Initializes the UKF on the first measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::FirstMeasurement(MeasurementPackage meas_package) {
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
void UKF::PredictSigmaPoints(MatrixXd& Xsig_aug, long delta_t) {
    // Predict sigma points k + 1
    for (int i = 0; i < n_sig_columns_; i++) {
        //extract vars for readability.
        double px = Xsig_aug(0, i);
        double py = Xsig_aug(1, i);
        double vel = Xsig_aug(2, i);
        double yaw = Xsig_aug(3, i);
        double yawd = Xsig_aug(4, i);
        double nu_a = Xsig_aug(5, i);
        double nu_yawdd = Xsig_aug(6, i);
        
        //calc yaw rate with time delta
        double newYawRateDelta = yawd * delta_t;
        
        // calc vel with respect to yaw rate.
        double velYawd = vel/yawd;
        
        // square delta_t
        double deltaTSquared = delta_t * delta_t;
        
        //predicted state values
        double predPx;
        double predPy;
        
        if (abs(yawd) > .001) {
            predPx = px + velYawd * (sin(yaw + newYawRateDelta) - sin(yaw));
            predPy = py + velYawd * (cos(yaw) - cos(yaw + newYawRateDelta));
        } else {
            predPx = px + vel * cos(yaw) * delta_t;
            predPy = py + vel * sin(yaw) * delta_t;
        }
        
        double predVel = vel;
        double predYaw = yaw + yawd * delta_t;
        double predYawd = yawd;
        
        //add noise
        predPx = predPx + .5 * deltaTSquared * cos(yaw) * nu_a;
        predPy = predPy + .5 * deltaTSquared * sin(yaw) * nu_a;
        predVel = predVel + delta_t * nu_a;
        predYaw = predYaw + .5 * deltaTSquared * nu_yawdd;
        predYawd = predYawd + delta_t * nu_yawdd;
        
        //write column
        Xsig_pred_(0, i) = predPx;
        Xsig_pred_(1, i) = predPy;
        Xsig_pred_(2, i) = predVel;
        Xsig_pred_(3, i) = predYaw;
        Xsig_pred_(4, i) = predYawd;
    }
}
