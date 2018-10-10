//
// Created by yijiangh on 10/10/18.
//

#include <string>
#include <iostream>

#include <eigen3/Eigen/Dense>

void testMatrixPermutation()
{
  using namespace Eigen;

  Eigen::Matrix4d m;
  m << 1, 2, 3, 4,
      5, 6, 7, 8,
      9, 10, 11, 12,
      13, 14, 15, 16;
  std::cout << m << std::endl;

  Eigen::Vector4i p;
  p << 4,2,1,3;

  Eigen::PermutationMatrix<4, 4> P;
  P.indices() = p;

  auto P_mat = P.toDenseMatrix();
//  std::cout << P_mat << std::endl;
//  std::cout << P_mat.inverse() << std::endl;
  auto s = P_mat * m;
  std::cout << s << std::endl;

//  std::cout << "right perm:\n" << m * P_mat.inverse() << std::endl;
//
//  std::cout << "both-side perm:\n" << P_mat * m * P_mat.inverse() << std::endl;

//  Eigen::Transpositions<4> Tr(p);
//
//  std::cout << "Tr left perm:\n" << Tr * m << std::endl;
//  std::cout << "Tr right perm:\n" << m * Tr << std::endl;
//
//  std::cout << "Tr both-side perm:\n" << Tr * m * Tr << std::endl;
}

int main(int argc, char** argv)
{
  testMatrixPermutation();
  return 0;
}