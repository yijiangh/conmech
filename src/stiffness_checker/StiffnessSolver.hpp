#pragma once

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>

#include "stiffness_checker/Timer.hpp"

using namespace std;

namespace conmech
{
namespace stiffness_checker
{

class StiffnessSolver
{
 private:
  typedef Eigen::SparseMatrix<double> SpMat;
  typedef Eigen::MatrixXd MX;
  typedef Eigen::Matrix3d M3;
  typedef Eigen::VectorXd VX;
  typedef Eigen::Vector3d V3;
  typedef Eigen::VectorXi VXi;
  typedef Eigen::MatrixXi MXi;

 public:
  StiffnessSolver();
  ~StiffnessSolver() {};

 public:
  /* solver I/O*/

  /*
  * SolveSystem  -  solve {F} =   [K]{D} via L D L' decomposition        2/Dec/2015
  *                 This override function is implemented for sovling Kqq * d_q = F_q,
  *                 where no partitioned LDLt is needed.
  * @param K   : stiffness matrix for the restrained frame
  * @param D   : displacement vector to be solved
  * @param F   : mechenical force(external load vector)
  * @param verbose	:  =1 for copious screenplay
  * @param info: <0 : not stiffness matrix positive definite
  * Note: This function use eigen library SimplicialLDLt module to solve the linear system.
  */
  bool SolveSystem(
      SpMat &K,
      VX &D,
      VX &F,
      int verbose,
      int &info
  );

  /*
  * SolveSystem  -  solve {F} =   [K]{D} via conjugate gradient with guess       30/Mar/2016
  *                 This override function is implemented for sovling Kqq * d_q = F_q,
  *                 where no partitioned LDLt is needed.
  * @param K   : stiffness matrix for the restrained frame
  * @param D   : displacement vector to be solved
  * @param F   : mechenical force(external load vector)
  * @param verbose	:  =1 for copious screenplay
  * @param info: <0 : not stiffness matrix positive definite
  * Note: This function use eigen library :ConjugateGradient module to solve the linear system.
  */
  bool SolveSystem(
      SpMat &K,
      VX &D,
      VX &F,
      VX &D0,
      int verbose,
      int &info
  );

  void Debug();

 public:
  /*
  * This LUDecomp module use Eigen library to solve the linear system
  */
  bool LUDecomp(
      MX &A,
      VX &x,
      VX &b
  );

 public:
  Timer compute_k_;
  Timer solve_d_;

  /* Timing Stat */
  bool detailed_timing_;
};

} // namespace stiffness_checker
} // namespace conmech