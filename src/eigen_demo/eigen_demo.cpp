#include <string>
#include <iostream>

// #include <eigen3/Eigen/Dense>
#include <Eigen/Dense>
#include <eigen_demo/eigen_demo.h>
// #include "eigen_demo/eigen_demo_IO.hpp"

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

// void EigenSolveDemo::testJsonParse(const std::string& file_path)
// {
//   // parseJson(file_path);
// }
