#ifndef UKF_H
#define UKF_H

#include "measurement_package.h"
#include "nis_tracker.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF {
public:
    
    ///* initially set to false, set to true in first call of ProcessMeasurement
    bool is_initialized_;
    
    ///* if this is false, laser measurements will be ignored (except for init)
    bool use_laser_;
    
    ///* if this is false, radar measurements will be ignored (except for init)
    bool use_radar_;
    
    ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
    VectorXd x_;
    
    ///* state covariance matrix
    MatrixXd P_;
    
    ///* predicted sigma points matrix
    MatrixXd Xsig_pred_;
    
    ///* time when the state is true, in us
    long long time_us_;
    
    ///* Process noise standard deviation longitudinal acceleration in m/s^2
    double std_a_;
    
    ///* Process noise standard deviation yaw acceleration in rad/s^2
    double std_yawdd_;
    
    //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
    // Laser measurement noise standard deviation position1 in m
    static constexpr double std_laspx_ = 0.15;
    
    // Laser measurement noise standard deviation position2 in m
    static constexpr double std_laspy = 0.15;
    
    // Radar measurement noise standard deviation radius in m
    static constexpr double std_radr_ = 0.3;
    
    // Radar measurement noise standard deviation angle in rad
    static constexpr double std_radphi_ = 0.03;
    
    // Radar measurement noise standard deviation radius change in m/s
    static constexpr double std_radrd_ = 0.3;
    //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
    
    ///* Weights of sigma points
    VectorXd weights_;
    
    ///* State dimension
    int n_x_;
    
    ///* Augmented state dimension
    int n_aug_;
    
    ///* Number of sigma point columns
    int n_sig_columns_;
    
    ///* Vector size radar measurement
    int n_radar;
    
    ///* Vector size lidar measurement
    int n_lidar;
    
    ///* Sigma point spreading parameter
    double lambda_;
    
    ///* timestamp of previous time step.
    double previous_timestamp_;
    
    ///* Linear measurement matrix for linear measurements
    MatrixXd H_;
    
    ///* Laser noise matrix
    MatrixXd R_laser_;
    
    ///* Radar noise Matrix
    MatrixXd R_radar_;
    
    /**
     * Constructor
     */
    UKF();
    
    /**
     * Destructor
     */
    virtual ~UKF();
    
    /**
     * ProcessMeasurement
     * @param meas_package The latest measurement data of either radar or laser
     */
    void ProcessMeasurement(MeasurementPackage meas_package);
    
    /**
     * Prediction Predicts sigma points, the state, and the state covariance
     * matrix
     * @param delta_t Time between k and k+1 in s
     */
    void Prediction(double delta_t);
    
    /**
     * Updates the state and the state covariance matrix using a laser measurement
     * @param meas_package The measurement at k+1
     */
    void UpdateLidar(MeasurementPackage meas_package);
    
    /**
     * Updates the state and the state covariance matrix using a radar measurement
     * @param meas_package The measurement at k+1
     */
    void UpdateRadar(MeasurementPackage meas_package);
private:
    NisTracker nis_lidar;
    NisTracker nis_radar;
    /**
     * Sets the first measurement of the UKF.
     * @param meas_package The measurement
     */
    void FirstMeasurement(MeasurementPackage meas_package);
    /**
     * Predicts sigma points from augmented sigma point matrix step k for step k+1
     * @param {MatrixXd} Xsig_aug
     * @param {MatrixXd} Xsig_aug
     */
    void PredictSigmaPoints(MatrixXd& Xsig_aug, double delta_t);
    /**
     * Predicts determines if measurement is supported baesd on implementation support and configuration
     * @returns {bool}
     */
    bool IsMeasurementEnabled(MeasurementPackage::SensorType type);
    /**
     * Normalizes an angle between pi and -pi
     * @param {double} angle
     * @returns {float} new angle
    */
    double NormalizeAngle(double angle);
};

#endif /* UKF_H */
