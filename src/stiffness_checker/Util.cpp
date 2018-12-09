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

void getGlobal2LocalRotationMatrix(
    const Eigen::Vector3d& end_vert_u,
    const Eigen::Vector3d& end_vert_v,
    Eigen::Matrix3d& rot_m,
    const double& rot_y2x)
{
  // length of the element
  double L = (end_vert_v - end_vert_u).norm();

  // by convention, the new x axis is along the element's direction
  // directional cosine of the new x axis in the global world frame
  double c_x, c_y, c_z;
  c_x = (end_vert_v[0] - end_vert_u[0]) / L;
  c_y = (end_vert_v[1] - end_vert_u[1]) / L;
  c_z = (end_vert_v[2] - end_vert_u[2]) / L;

  rot_m = Eigen::Matrix3d::Zero();

  if (abs(c_z) == 1.0)
  {
    // the element is parallel to global z axis
    // cross product is not defined, in this case
    // it's just a rotation about the global z axis
    // in x-y plane
    rot_m(2,0) = 1.0;
    rot_m.block<2,2>(0,1) = Eigen::Rotation2Dd(rot_y2x).toRotationMatrix();
  }
  else
  {
    // local x_axis = element's vector
    auto new_x = Eigen::Vector3d(c_x,c_y,c_z);
    rot_m.block<1,3>(0,0) = new_x;

    // local y axis = cross product with global z axis
    auto new_y = Eigen::Vector3d::UnitZ().cross(new_x);
    Eigen::AngleAxisd aa(rot_y2x, new_x);
    new_y = aa * new_y;
    rot_m.block<1,3>(1,0) = new_y;

    rot_m.block<1,3>(2,0) = new_x.cross(new_y);
  }
}

} // namespace stiffness_checker
} // namespace conmech
