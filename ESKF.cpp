#include "ESKF.hpp"

#include "matplotlibcpp.h"
namespace plt = matplotlibcpp;


ESKF::ESKF() : wander_ang(0)
{
    cout << "ESKF has been created..." << endl;
}

void ESKF::set_eskf(int number_of_state, int model_order, const VectorXd& p_matrix_values, const VectorXd& q_matrix_values, const double r_noise)
{
    this->number_of_state = number_of_state;
    this->model_order = model_order;

    // ---------------------------- Kalman Filter Initialization -----------------------------
    h_matrix.setZero(number_of_meas, number_of_state);
    h_matrix.block<3,3>(0,6) = -Matrix3d::Identity();
    f_matrix.setZero(number_of_state, number_of_state);
    phi_matrix.setZero(number_of_state, number_of_state);
    k_matrix.setZero(number_of_state, number_of_meas);
    z_vector.setZero(number_of_meas);
    S_matrix.setZero(number_of_meas, number_of_meas);
    S_mat_diag.setZero(number_of_meas);

    // ---------------------------------- P~Q~R Matrix Value ----------------------------------
    p_matrix.setZero(number_of_state, number_of_state);
    q_matrix.setZero(number_of_state, number_of_state);
    r_matrix.setZero(number_of_meas, number_of_meas);

    p_matrix.block<3,3>(0,0)      = Matrix3d::Identity() * p_matrix_values(0);  // p_att;
    p_matrix.block<3,3>(3,3)      = Matrix3d::Identity() * p_matrix_values(1);  // p_vel;
    p_matrix.block<3,3>(6,6)      = Matrix3d::Identity() * p_matrix_values(2);  // p_pos;
    p_matrix.block<2,2>(9,9)      = Matrix2d::Identity() * p_matrix_values(3);  // p_wander;
    p_matrix.block<3,3>(11,11)    = Matrix3d::Identity() * p_matrix_values(4);  // p_ba;
    p_matrix.block<3,3>(14,14)    = Matrix3d::Identity() * p_matrix_values(5);  // p_bg;
    q_matrix.block<3,3>(0,0)      = Matrix3d::Identity() * q_matrix_values(0);  // q_att;
    q_matrix.block<3,3>(3,3)      = Matrix3d::Identity() * q_matrix_values(1);  // q_vel;
    q_matrix.block<3,3>(6,6)      = Matrix3d::Identity() * q_matrix_values(2);  // q_pos;
    q_matrix.block<2,2>(9,9)      = Matrix2d::Identity() * q_matrix_values(3);  // q_wander;
    q_matrix.block<3,3>(11,11)    = Matrix3d::Identity() * q_matrix_values(4);  // q_ba;
    q_matrix.block<3,3>(14,14)    = Matrix3d::Identity() * q_matrix_values(5);  // q_bg;
    r_matrix.block<3,3>(0,0)    = MatrixXd::Identity(3,3) * r_noise;

    w_ie_w.setZero(3);
    b_a.setZero(3);
    b_g.setZero(3);
    
    nominal_state.setZero(number_of_state);
    error_state.setZero(number_of_state);
    true_state.setZero(number_of_state);
    
    f_ib_b.setZero(3);
    w_ib_b.setZero(3);

    cout << "ESKF is set." << endl;
}

void ESKF::time_update(double sampling_periode, double latitude, double height)
{
    mechanization_wander(sampling_periode,
                         latitude,
                         height,
                         true_state,
                         nominal_state,
                         f_ib_b,
                         w_ib_b,
                         wander_ang,
                         w_ie_w);

    // Create state space matrix
    f_matrix = state_space_matrix(latitude);

        
    // P and Q updates
    phi_matrix = state_transition(sampling_periode);

    p_matrix = phi_matrix * p_matrix * phi_matrix.transpose() + q_matrix;
}

void ESKF::measurement_update(void)
{
    // Measurement Update (Correction) Stage
    S_matrix = (h_matrix * p_matrix * h_matrix.transpose() + r_matrix);  // Error Covariance Matrix

    S_mat_diag = S_matrix.diagonal();

    k_matrix = p_matrix * h_matrix.transpose() * S_matrix.inverse();
        
    // Calculate Innovation
    z_vector = -nominal_state.segment<3>(6); // write index 7 to 9
 
    error_state += k_matrix * z_vector;  

    p_matrix = (MatrixXd::Identity(17, 17) - k_matrix * h_matrix) * p_matrix * (MatrixXd::Identity(17, 17) - k_matrix * h_matrix).transpose() + k_matrix * r_matrix * k_matrix.transpose();

    CL_correction(/*eskf_.true_state, wander_ang, eskf_.nominal_state, eskf_.error_state*/);

    cout << "Wander Angle:\n" << wander_ang << endl;
    cout << "True State Vector: \n" << rad_to_deg_v(true_state.segment<3>(0)) << endl;

}

MatrixXd ESKF::state_space_matrix(double latitude_b)
{
    Eigen::Vector3d temp = nominal_state.segment<3>(0);
    Matrix3d c_b_w = Euler_to_CTM( temp);
    MatrixXd ss_mat(17, 17);
    ss_mat.block<17,17>(0,0)    = MatrixXd::Zero(17,17);
    ss_mat.block<3,3>(0,0)      = -skew_symmetric(w_ie_w);
    ss_mat.block<3,1>(0,9)      = Vector3d(0, w_ie * cos(latitude_b),0);
    ss_mat.block<3,1>(0,10)     = Vector3d(-w_ie * cos(latitude_b),0,0);
    ss_mat.block<3,3>(0,14)     = c_b_w;
    ss_mat.block<3,3>(3,0)      = -skew_symmetric(c_b_w * f_ib_b);
    ss_mat.block<3,3>(3,11)     = c_b_w;
    ss_mat.block<3,3>(6,3)      = Matrix3d::Identity();
 
    return ss_mat;
}

MatrixXd ESKF::state_transition(double Ts)
{
    int dim = f_matrix.rows();
    MatrixXd Phi = MatrixXd::Identity(dim, dim);

    switch (model_order) {
        case 1:
            Phi += (f_matrix * Ts);
            break;
        case 2:
            Phi += (f_matrix * Ts) + (0.5 * f_matrix * f_matrix * Ts * Ts);
            break;
        case 3:
            Phi += (f_matrix * Ts) + (0.5 * f_matrix * f_matrix * Ts * Ts) + ((1.0 / 6.0) * f_matrix * f_matrix * f_matrix * Ts * Ts * Ts);
            break;
        case 4:
            Phi += (f_matrix * Ts) + (0.5 * f_matrix * f_matrix * Ts * Ts) + ((1.0 / 6.0) * f_matrix * f_matrix * f_matrix * Ts * Ts * Ts) + ((1.0 / 24.0) * f_matrix * f_matrix * f_matrix * f_matrix * Ts * Ts * Ts * Ts);
            break;
        case 5:
            Phi += (f_matrix * Ts) + (0.5 * f_matrix * f_matrix * Ts * Ts) + ((1.0 / 6.0) * f_matrix * f_matrix * f_matrix * Ts * Ts * Ts) + ((1.0 / 24.0) * f_matrix * f_matrix * f_matrix * f_matrix * Ts * Ts * Ts * Ts) + ((1.0 / 120.0) * f_matrix * f_matrix * f_matrix * f_matrix * f_matrix * Ts * Ts * Ts * Ts * Ts);
            break;
        default:
            cerr << "Invalid model order." << endl;
            exit(1);
            break;
    }

    return Phi;
}

void ESKF::CL_correction(void)
{
    // Attitude update
    true_state.segment<3>(0) = CTM_to_Euler((Matrix3d::Identity() - skew_symmetric(error_state.segment<3>(0))) * Euler_to_CTM(nominal_state.segment<3>(0)));
 
    // Vel, pos, sinW, and cosW correction
    true_state.segment<3>(3) = nominal_state.segment<3>(3) - error_state.segment<3>(3); // Velocity correction
    true_state.segment<3>(6) = nominal_state.segment<3>(6) - error_state.segment<3>(6); // Position correction
 
    true_state.segment<2>(9) = nominal_state.segment<2>(9) - error_state.segment<2>(9); // Vel, pos, sinW and cosW correction
    // b_a correction
    true_state.segment<3>(11) = nominal_state.segment<3>(11) + error_state.segment<3>(11);
 
    // b_g correction
    true_state.segment<3>(14) = nominal_state.segment<3>(14) + error_state.segment<3>(14);

    // Wander Angle [Rad/Sn] Update
    wander_ang = atan2(true_state.coeff(9), true_state.coeff(10));
}

void ESKF::plot_output(void)
{
    std::vector<double> true_state_std(true_state.data(), true_state.data() + true_state.size());
    std::vector<double> Smat_diag_std(S_mat_diag.data(), S_mat_diag.data() + S_mat_diag.size());
    std::vector<double> z_vector_std(z_vector.data(), z_vector.data() + z_vector.size());
 
        // Algorithm State Vector Size   "Attitude(3x1) Velocity(3x1) Pose(3x1) WanderSin(1x1) WanderCos(1x1) Ba(3x1) Bg(3x1)"
        plot_iter++;
        iterations.push_back(plot_iter);
        values_attX.push_back(true_state_std[0]);
        values_attY.push_back(true_state_std[1]);
        values_attZ.push_back(true_state_std[2]);
        values_velX.push_back(true_state_std[3]);
        values_velY.push_back(true_state_std[4]);
        values_velZ.push_back(true_state_std[5]);
        values_posX.push_back(true_state_std[6]);
        values_posY.push_back(true_state_std[7]);
        values_posZ.push_back(true_state_std[8]);
        values_sinWander.push_back(true_state_std[9]);
        values_cosWander.push_back(true_state_std[10]);
        values_checkWander.push_back(sqrt(pow(true_state_std[9],2)+pow(true_state_std[10],2)));
        values_baX.push_back(true_state_std[11]);
        values_baY.push_back(true_state_std[12]);
        values_baZ.push_back(true_state_std[13]);
        values_bgX.push_back(true_state_std[14]);
        values_bgY.push_back(true_state_std[15]);
        values_bgZ.push_back(true_state_std[16]);


        //ONLINE PLOTTING
        if (plot_iter % 2000 == 0 || plot_iter == 0) { // Update the plot every xx iterations,
 
            plt::clf(); // Clear current figure
            plt::subplot(3, 3, 1);
            plt::named_plot("attX", iterations, values_attX, "r--");
            plt::named_plot("attY", iterations, values_attY, "g--");
            plt::named_plot("attZ", iterations, values_attZ, "b--");
            plt::legend();
 
            plt::subplot(3, 3, 2);
            plt::named_plot("velX", iterations, values_velX, "r--");
            plt::named_plot("velY", iterations, values_velY, "g--");
            plt::named_plot("velZ", iterations, values_velZ, "b--");
            plt::legend();
 
            plt::subplot(3, 3, 3);
            plt::named_plot("posX", iterations, values_posX, "r--");
            plt::named_plot("posY", iterations, values_posY, "g--");
            plt::named_plot("posZ", iterations, values_posZ, "b--");
            plt::legend();
 
            plt::subplot(3, 3, 4);
            plt::named_plot("sinWan", iterations, values_sinWander, "r--");
            plt::named_plot("cosWan", iterations, values_cosWander, "g--");
            plt::named_plot("Check", iterations, values_checkWander, "b--");
            plt::legend();
 
            plt::subplot(3, 3, 5);
            plt::named_plot("b_aX", iterations, values_baX, "r--");
            plt::named_plot("b_aY", iterations, values_baY, "g--");
            plt::named_plot("b_aZ", iterations, values_baZ, "b--");
            plt::legend();
 
            plt::subplot(3, 3, 6);
            plt::named_plot("b_gX", iterations, values_bgX, "r--");
            plt::named_plot("b_gY", iterations, values_bgY, "g--");
            plt::named_plot("b_gZ", iterations, values_bgZ, "b--");
            plt::legend();
 
            plt::pause(0.1); 
        } //end if plot
}