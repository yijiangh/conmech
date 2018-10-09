#include "stiffness_checker/Util.h"

namespace conmech
{
namespace stiffness_checker
{

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
    rot_m(0,2) = 1.0;
    rot_m.block<2,2>(1,0) = Eigen::Rotation2Dd(rot_y2x).toRotationMatrix();
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

//void CoordTrans::TransLocToGlob(
//    double t0, double t1, double t2, double t3, double t4,
//    double t5, double t6, double t7, double t8,
//    Eigen::MatrixXd	   &m, float r1, float r2)
//{
//  int     i, j, k;
//
//  Eigen::MatrixXd t(12, 12);
//  Eigen::MatrixXd mt(12, 12);
//
//  t.setZero();
//  mt.setZero();
//
//  for (i = 0; i < 4; i++)
//  {
//    t(3 * i, 3 * i) = t0;
//    t(3 * i, 3 * i + 1) = t1;
//    t(3 * i, 3 * i + 2) = t2;
//    t(3 * i + 1, 3 * i) = t3;
//    t(3 * i + 1, 3 * i + 1) = t4;
//    t(3 * i + 1, 3 * i + 2) = t5;
//    t(3 * i + 2, 3 * i) = t6;
//    t(3 * i + 2, 3 * i + 1) = t7;
//    t(3 * i + 2, 3 * i + 2) = t8;
//  }
//
//  /*  effect of finite node radius on coordinate transformation  is not supported now */
//
//  mt = t.transpose() * m * t;
//  m = mt;
//}

} // namespace stiffness_checker
} // namespace conmech
