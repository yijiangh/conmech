#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "stiffness_checker/Timer.h"

using namespace std;
// TODO: clean up namespace
// TODO: no need for having a class

namespace conmech
{
namespace stiffness_checker
{

class StiffnessSolver
{
 public:
  StiffnessSolver() : timing_(false) {}
  ~StiffnessSolver() {}

 public:
  bool solveSystemLU(
    const Eigen::MatrixXd& A, const Eigen::VectorXd& b,
    Eigen::VectorXd& x);

  bool solveSparseSimplicialLDLT(
    const Eigen::SparseMatrix<double>& A, const Eigen::VectorXd& b,
    Eigen::VectorXd& x, const bool& verbose=false);

 public:
  Timer solve_timer_;

  /* Timing Stat */
  bool timing_;
};

} // namespace stiffness_checker
} // namespace conmech
