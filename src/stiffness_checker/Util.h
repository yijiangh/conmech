#pragma once

#include <Eigen/Dense>

namespace conmech
{
namespace stiffness_checker
{

// Input: rot_y2x
void getGlobal2LocalRotationMatrix(
    const Eigen::Vector3d& end_vert_u,
    const Eigen::Vector3d& end_vert_v,
    Eigen::Matrix3d& rot_m,
    const double& rot_y2x=0.0);

} // namespace stiffness_checker
} // namespace conmech
