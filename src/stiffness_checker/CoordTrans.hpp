#pragma once

#include <vector>
#include <eigen3/Eigen/Dense>

#include "stiffness_checker/Vec.hpp"
#include "stiffness_checker/Frame.h"

/**
CoordTrans -  evaluate the 3D coordinate transformation coefficients
Default order of coordinate rotations...  typical for Y as the vertical axis
1. rotate about the global Z axis
2. rotate about the global Y axis
3. rotate about the local  x axis --- element 'roll'

If Zvert is defined as 1, then the order of coordinate rotations is typical
for Z as the vertical axis
1. rotate about the global Y axis
2. rotate about the global Z axis
3. rotate about the local  x axis --- element 'roll'

Q=TF;   U=TD;   T'T=I;   Q=kU;   TF=kTD;   T'TF=T'kTD;   T'kT = K;   F=KD
*/

namespace conmech
{
namespace stiffness_checker
{

class CoordTrans
{
 public:
  CoordTrans(){};
  ~CoordTrans(){};

  void CreateTransMatrix(
		  std::vector<Eigen::Vector3d>	xyz,
		  double L,			// length of the element(edge)
		  int n1, int n2,		// index fo endpoint of the element
		  double &t0, double &t1, double &t2, double &t3, double &t4,
		  double &t5, double &t6, double &t7, double &t8,
		  float p);

  void CreateTransMatrix(
		  trimesh::point u, trimesh::point v,
		  double &t0, double &t1, double &t2, double &t3, double &t4,
		  double &t5, double &t6, double &t7, double &t8,
		  float p);

  void TransLocToGlob(
		  double t0, double t1, double t2, double t3, double t4,
		  double t5, double t6, double t7, double t8,
		  Eigen::MatrixXd &m, float r1, float r2);
};

} // namespace stiffness_checker
} // namespace conmech
