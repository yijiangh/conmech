//
// Created by yijiangh on 9/24/18.
//

#ifndef CONMECH_EIGEN_DEMO_H
#define CONMECH_EIGEN_DEMO_H

#include <string>
#include <iostream>

#include <eigen3/Eigen/Dense>
#include "eigen_demo/eigen_demo.hpp"

void EigenSolveDemo::testEigen()
{
  using namespace std;
  using namespace Eigen;

  Matrix3f A;
  Vector3f b;
  A << 1,2,3,  4,5,6,  7,8,10;
  b << 3, 3, 4;
  cout << "Here is the matrix A:\n" << A << endl;
  cout << "Here is the vector b:\n" << b << endl;
  Vector3f x = A.colPivHouseholderQr().solve(b);
  cout << "The solution is:\n" << x << endl;
}

#endif //CONMECH_EIGEN_DEMO_H