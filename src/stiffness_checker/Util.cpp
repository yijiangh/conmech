#include <math.h>
#include <Eigen/Geometry>

#include "stiffness_checker/Util.h"

namespace conmech
{
namespace stiffness_checker
{

void createLocalStiffnessMatrix(const double &L, const double &A, const int &dim,
                                const double &Jx, const double &Iy, const double &Iz,
                                const double &E, const double &G, const double &mu,
                                Eigen::MatrixXd &K_eL)
{
  // TODO: add 2D case
  assert(3 == dim);

  switch(dim)
  {
    case 2:
      assert(false && "2D local stiffness not implemented.");
      return;
    case 3:
    {
      K_eL = Eigen::MatrixXd::Zero(12,12);

      // see: [Matrix Structural Analysis, McGuire et al., 2rd edition]
      // P73 - eq(4.34)
      Eigen::MatrixXd K_block(6,6);
      K_block.setZero();
      Eigen::VectorXd diag(6);

      // block_00 and block_11
      K_block(1,5) = 6*Iz / std::pow(L,2);
      K_block(2,4) = - 6*Iy / std::pow(L,2);
      K_block = K_block.eval() +  K_block.transpose().eval();
      K_eL.block<6,6>(0,0) = K_block;
      K_eL.block<6,6>(6,6) = -K_block;

      diag[0] = A/L;
      diag[1] = 12*Iz / std::pow(L,3);
      diag[2] = 12*Iy / std::pow(L,3);
      diag[3] = Jx / (2*(1+mu)*L);
      diag[4] = 4*Iy / L;
      diag[5] = 4*Iz / L;
      K_eL.block<6,6>(0,0) += Eigen::MatrixXd(diag.asDiagonal());
      K_eL.block<6,6>(6,6) += Eigen::MatrixXd(diag.asDiagonal());

      // block_01 and block_10
      K_block.setZero();

      K_block(1,5) = 6*Iz / std::pow(L,2);
      K_block(2,4) = - 6*Iy / std::pow(L,2);
      K_block = K_block.eval() - K_block.transpose().eval();
      K_eL.block<6,6>(0,6) = K_block;
      K_eL.block<6,6>(6,0) = -K_block;

      diag[0] = -A/L;
      diag[1] = -12*Iz / std::pow(L,3);
      diag[2] = -12*Iy / std::pow(L,3);
      diag[3] = -Jx / (2*(1+mu)*L);
      diag[4] = 2*Iy / L;
      diag[5] = 2*Iz / L;
      K_eL.block<6,6>(0,6) += Eigen::MatrixXd(diag.asDiagonal());
      K_eL.block<6,6>(6,0) += Eigen::MatrixXd(diag.asDiagonal());

      K_eL *= E;
    }
  }
}

/**
 * @brief Get the Global to Local Rotation Matrix object
 * Calculates a 3x3 matrix to tranform the global xyz axis to the element local axis.
 * TODO: add info on the axis convention, we have different conventions with compas_fea (abaqus)
 * 
 * The coordinate transformation matrix can be used to:
 *  - transform frame element end forces from the element (local) coordinate system
 *    to the structure (global) coordinate system
 *  - transfrom end displacements from the structural (global) coordinate system 
 *    to the element (local) coordinate system,
 *  - transform the frame element stiffness and mass matrices
 *    from element (local) coordinates to structral (global) coordinates.
 * Symbolically, the return matrix R = {local}_R_{global}
 * 
 * @param[in] end_vert_u 
 * @param[in] end_vert_v 
 * @param[out] rot_m 3x3 Eigen matrix, transforming global axis to local coordinate frame
 * @param[in] rot_y2x optional rotation of local y axis around the local x axis, defaults to zero
 */
void getGlobal2LocalRotationMatrix(
    const Eigen::VectorXd & end_vert_u,
    const Eigen::VectorXd & end_vert_v,
    Eigen::Matrix3d& rot_m,
    const double& rot_y2x)
{
  assert(end_vert_u.size() == end_vert_v.size() && "vert dimension not agree!");
  assert(end_vert_u.size() == 2 || end_vert_u.size() == 3);
  int dim = end_vert_u.size();

  // length of the element
  double L = (end_vert_v - end_vert_u).norm();
  // TODO: make tol as a common shared const
  assert(L < 1e6 && "vertices too close, might be duplicated pts.");

  // by convention, the new x axis is along the element's direction
  // directional cosine of the new x axis in the global world frame
  double c_x, c_y;
  c_x = (end_vert_v[0] - end_vert_u[0]) / L;
  c_y = (end_vert_v[1] - end_vert_u[1]) / L;

  Eigen::Matrix3d R = Eigen::Matrix3d::Zero();

  if (3 == dim)
  {
    double c_z = (end_vert_v[2] - end_vert_u[2]) / L;
    auto rot_axis = Eigen::AngleAxisd(rot_y2x, Eigen::Vector3d::UnitZ());

    if (abs(c_z) == 1.0)
    {
      // the element is parallel to global z axis
      // cross product is not defined, in this case
      // it's just a rotation about the global z axis
      // in x-y plane
      R(0, 2) = -c_z;
      R(1, 1) = 1;
      R(2, 0) = c_z;
    }
    else
    {
      // local x_axis = element's vector
      auto new_x = Eigen::Vector3d(c_x, c_y, c_z);

      // local y axis = cross product with global z axis
      Eigen::Vector3d new_y = -new_x.cross(Eigen::Vector3d::UnitZ());
      new_y.normalize();

      auto new_z = new_x.cross(new_y);

      R.block<3, 1>(0, 0) = new_x;
      R.block<3, 1>(0, 1) = new_y;
      R.block<3, 1>(0, 2) = new_z;
    }
    // This is essential!
    R = R * rot_axis;
    rot_m = R.transpose();
  }
  else
  {
    // 2D rotational matrix
    R(0,0) = c_x;
    R(0,1) = c_y;
    R(1,0) = -c_y;
    R(1,1) = c_x;
    R(2,2) = 1;

    auto rot_axis = Eigen::AngleAxisd(rot_y2x, Eigen::Vector3d::UnitZ());
    assert((R - rot_axis.toRotationMatrix()).norm() > 1e-3);

    rot_m = R;
  }
}

void getNodePoints(const Eigen::MatrixXd& Vertices, const int& end_u_id, const int& end_v_id, 
  Eigen::VectorXd& end_u, Eigen::VectorXd& end_v)
{
  end_u = Eigen::VectorXd(3);
  end_u << Vertices(end_u_id, 0), Vertices(end_u_id, 1), Vertices(end_u_id, 2);
  end_v = Eigen::VectorXd(3);
  end_v << Vertices(end_v_id, 0), Vertices(end_v_id, 1), Vertices(end_v_id, 2);
}

} // namespace stiffness_checker
} // namespace conmech
