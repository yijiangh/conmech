#pragma once

#include <eigen3/Eigen/Dense>
#include "stiffness_checker/Timer.h"

using namespace std;

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

  bool solveSystemLU(const Eigen::MatrixXd& A, const Eigen::VectorXd& b, Eigen::VectorXd& x);

  //  bool SolveSystem(
//      SpMat &K,
//      VX &D,
//      VX &F,
//      int verbose,
//      int& info
//  );

 public:
  Timer solve_timer_;

  /* Timing Stat */
  bool timing_;
};

} // namespace stiffness_checker
} // namespace conmech