#pragma once

#include <Eigen/Core>
#include <Eigen/Dense>

namespace conmech
{
namespace stiffness_checker
{
/**
 * construct local element stiffness matrix in the element's local frame
 * length unit: meter, force unit: kN
 * @param L element's length
 * @param A cross section area
 * @param dim dimension, = 2 or 3
 * @param Jx torsion const around local x axis
 * @param Iy area moment of inertia around local y axis
 * @param Iz area moment of inertia around local y axis
 * @param E young's modulus
 * @param G shear modulus
 * @param mu poisson ratio
 * @param[out] K_loc returned (2*full_node_dof x 2*full_node_dof) local element stiffness matrix
 * in 2D case: 6*6 matrix (full_node_dof = 3)
 * in 3D case: 12*12 matrix (full_node_dof = 6)
 */
void createLocalStiffnessMatrix(const double& L, const double& A, const int& dim,
                                const double& Jx, const double& Iy, const double& Iz,
                                const double& E, const double& G, const double& mu,
                                Eigen::MatrixXd& K_loc);

/**
 * Compute global to local (3x3) rotation matrix
 * @param end_vert_u
 * @param end_vert_v
 * @param[out] rot_m
 * @param rot_y2x
 */
void getGlobal2LocalRotationMatrix(
    const Eigen::VectorXd & end_vert_u,
    const Eigen::VectorXd & end_vert_v,
    Eigen::Matrix3d& rot_m,
    const double& rot_y2x=0.0);

void getNodePoints(const Eigen::MatrixXd& Vertices, const int& end_u_id, const int& end_v_id, 
  Eigen::VectorXd& end_u, Eigen::VectorXd& end_v);

} // namespace stiffness_checker
} // namespace conmech
