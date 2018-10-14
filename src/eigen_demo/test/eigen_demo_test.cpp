//
// Created by yijiangh on 10/10/18.
//

#include <string>
#include <iostream>

#include <eigen3/Eigen/Core>
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
  p << 3,1,0,2;
  std::cout << p << std::endl;

//  Eigen::PermutationMatrix<Dynamic, Dynamic> P(4);
////  P.setIdentity();
//  P.indices() = p;
//  std::cout << P.toDenseMatrix() << std::endl;
//  for(int i=0; i<4; i++)
//  {
//    P.indices()[i] = p(i);
//  }

  Matrix4d P;
  P.setZero();
  for(int i=0; i<4; i++)
  {
    P(i, p(i)) = 1;
  }
  std::cout << P << std::endl;

  std::cout << "row perm: \n" << P*m << std::endl;
  std::cout << "perm right: \n" << m * P << std::endl;
  std::cout << "col perm: \n" << m * P.inverse() << std::endl;
  std::cout << "both side perm: \n" << (P * m * P.inverse()) << std::endl;

  Eigen::Transpositions<4> Tr(4);
  Tr.setIdentity();
  for(int i=0; i<4; i++)
  {
    Tr.indices()[i] = p(i);
  }

//  std::cout << "Tr left perm: \n" << Tr * m << std::endl;
//  std::cout << "Tr right perm: \n" << m * Tr.inverse() << std::endl;
//  std::cout << "Tr both-side perm:\n" << Tr * m * Tr.inverse() << std::endl;
}

int main(int argc, char** argv)
{
  testMatrixPermutation();
  return 0;
}