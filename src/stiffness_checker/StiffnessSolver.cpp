#include <Eigen/LU>

#include "stiffness_checker/StiffnessSolver.h"

namespace conmech
{
namespace stiffness_checker
{

//bool StiffnessSolver::SolveSystem(SpMat &K, VX &D, VX &F, int verbose, int &info)
//{
//  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
//
//  if (detailed_timing_)
//  {
//    compute_k_.Start();
//  }
//
//  solver.compute(K);
//
//  if (detailed_timing_)
//  {
//    compute_k_.Stop();
//  }
//
//  info = 0;
//
//  if (solver.info() != Eigen::Success)
//  {
//    fprintf(stderr, "SolverSystem(LDLT): Error in Decomposition!\n");
//    return false;
//  }
//
//  VX Diag = solver.vectorD();
//
//  for (int i = 0; i < Diag.size(); i++)
//  {
//    if (Diag[i] == 0.0)
//    {
//      fprintf(stderr, " SolveSystem(LDLT): zero found on diagonal ...\n");
//      fprintf(stderr, " d[%d] = %11.4e\n", i, Diag[i]);
//      return false;
//    }
//
//    if (Diag[i] < 0.0)
//    {
//      fprintf(stderr, " SolveSystem(LDLT): negative number found on diagonal ...\n");
//      fprintf(stderr, " d[%d] = %11.4e\n", i, Diag[i]);
//      info--;
//      return false;
//    }
//  }
//
//  if (info < 0)
//  {
//    fprintf(stderr, "Stiffness Matrix is not positive definite: %d negative elements\n", info);
//    fprintf(stderr, "found on decomp diagonal of K.\n");
//    fprintf(stderr, "The stucture may have mechanism and thus not stable in general\n");
//    fprintf(stderr, "Please Make sure that all six\n");
//    fprintf(stderr, "rigid body translations are restrained!\n");
//
//    return false;
//  }
//
//  if (detailed_timing_)
//  {
//    solve_d_.Start();
//  }
//
//  D = solver.solve(F);
//
//  if (detailed_timing_)
//  {
//    solve_d_.Stop();
//  }
//
//  if (solver.info() != Eigen::Success)
//  {
//    fprintf(stderr, "SolverSystem(LDLT): Error in Solving!\n");
//    return false;
//  }
//
//  return true;
//}

bool StiffnessSolver::solveSystemLU(
    const Eigen::MatrixXd& A, const Eigen::VectorXd& b, Eigen::VectorXd& x)
{
  if(timing_)
  {
    solve_timer_.Start();
  }

  x = A.fullPivLu().solve(b);

  if(timing_)
  {
    solve_timer_.Stop();
  }

  if ((A*x).isApprox(b))
  {
    return true;
  }
  else
  {
//    std::cout << "A is invertible? - " << A.fullPivLu().isInvertible() << std::endl;
    return false;
  }
}

} // namespace stiffness_checker
} // namespace conmech
