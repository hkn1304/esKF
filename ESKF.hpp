#ifndef ESKF_HPP
#define ESKF_HPP

#include <iostream>
#include <fstream>
#include <sstream>

#include "math_utility.hpp"

using namespace Eigen;
using namespace std;


class ESKF 
{
  public:
    ESKF();

    /**
     * Set number of state and model order values of ESKF
     * 
     * @param number_of_state Number of state of ESKF
     * @param model_order Model order of ESKF
     * @param p_matrix_values p_att, p_vel, p_pos, p_wander, p_ba, p_bg
     * @param q_matrix_values q_att, q_vel, q_qos, q_wander, q_ba, q_bg
     * @param r_matrix_noise r_noise
     */
    void set_eskf(int number_of_state, int model_order, const VectorXd& p_matrix_values, const VectorXd& q_matrix_values, const double r_matrix_noise);

    /**
     * Time update section for ESKF
     * 
     * @param sampling_periode
     * @param latitude
     * @param height
    */
    void time_update(double sampling_periode, double latitude, double height);

    /**
     * Measurement update for ESKF
    */
    void measurement_update(void);

    /**
     * State space matrix
     * 
     * @param w_ie_w
     * @param l_b
     * @param nominal_state
     * @param f_ib_b
     */
    MatrixXd state_space_matrix(double latitude);

    /**
     * State transition
     * 
     * @param f_mat
     * @param Ts
     * @param model_order
     */
    MatrixXd state_transition(double Ts);

    /**
     * CL correction
     */
    void CL_correction(void);    

    // ---------------- PLOT VARIABLES ----------------
    int plot_iter=0;
    void plot_output(void);

    // ---------------------------------- Kalman Filter Initialization ----------------------------------
    // Number Of State "Attitude(3x1), Velocity(3x1), Pose(3x1), WSin(1x1), WCos(1x1), Ba(3x1), Bg(3x1)"
    MatrixXd h_matrix;    // State Selection Matrix
    MatrixXd f_matrix;    // State Space System Matrix
    MatrixXd phi_matrix;  // State Transition Matrix
    MatrixXd k_matrix;    // Kalman Gain Matrix
    VectorXd z_vector;    // Innovation matrix

    MatrixXd S_matrix;
    VectorXd S_mat_diag;

    // ---------------------------------------- P~Q~R Matrix Value --------------------------------------
    // High Kalman Gain:  if R >> P => P/(P+R) ~= K ~= 0 (Algorithm belief measurement)
    // Low  Kalman Gain:  if R << P => P/(P+R) ~= K ~= 1 (Algorithm belief kalman model)
    
    MatrixXd p_matrix;    // State Covariance Matrix
    double p_att, p_vel, p_pos, p_wander, p_ba, p_bg;

    MatrixXd q_matrix;    // Processes Noise Covariance Matrix
    double q_att, q_vel, q_pos, q_wander, q_ba, q_bg;

    MatrixXd r_matrix;    // Measurement Covariance Matrix
    double r_noise;

    // ------------------------------------- Algorithm Output Parameters --------------------------------
    double wander_ang;   //!< Wander angle [rad]

    VectorXd w_ie_w;   // Wander Frame Earth Angular Veocity Transform
    VectorXd b_a;      // Accel Dynamic Bias
    VectorXd b_g;      // Gyro  Dynamic Bias
    
    VectorXd nominal_state;   // Nominal State Vector (Measurement Value and Mechanization)
    VectorXd error_state;     // Error State Vector   (ESKF)
    VectorXd true_state;      // True State Vector = Nominal State Vector - Error State Vector
    VectorXd innovation;      // Kalman Filter Innovation Vector
    
    VectorXd f_ib_b;  // Acceleration Data
    VectorXd w_ib_b;  // Gyro Data

    //Plot variables
    std::vector<double> iterations;
    std::vector<double> values_attX;
    std::vector<double> values_attY;
    std::vector<double> values_attZ;
    std::vector<double> values_posX;
    std::vector<double> values_posY;
    std::vector<double> values_posZ;
    std::vector<double> values_velX;
    std::vector<double> values_velY;
    std::vector<double> values_velZ;
    std::vector<double> values_sinWander;
    std::vector<double> values_cosWander;
    std::vector<double> values_checkWander;
    std::vector<double> values_baX;
    std::vector<double> values_baY;
    std::vector<double> values_baZ;
    std::vector<double> values_bgX;
    std::vector<double> values_bgY;
    std::vector<double> values_bgZ;

  private:
    int model_order;      //!< State Transition Taylor Series Order  (Default Value of ESKF -> 3)
    int number_of_state;  //!< Algorithm Sate Vector Size            (Default Value of ESKF -> 17  "Attitude(3x1) Pose(3x1) Velocity(3x1) WanderSin(1x1) WanderCos(1x1) Ba(3x1) Bg(3x1)")
    int number_of_meas=3;
};

#endif // #ifndef ESKF_HPP