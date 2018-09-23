//#include "eigen_demo/eigen_demo_IO.cpp"

#include <string>
#include <iostream>

#include <eigen3/Eigen/Dense>

class EigenSolveDemo
{
 public:
  EigenSolveDemo(){}

  void testEigen()
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

 private:
  std::string file_path;
};