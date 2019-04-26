#include <Eigen/LU>
#include <Eigen/SparseCholesky>

#include "stiffness_checker/StiffnessSolver.h"

namespace conmech
{
namespace stiffness_checker
{

bool StiffnessSolver::solveSparseSimplicialLDLT(
    const Eigen::SparseMatrix<double>& A, const Eigen::VectorXd& b, Eigen::VectorXd& x)
{
  if(timing_)
  {
    solve_timer_.Start();
  }

  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
  solver.compute(A);

  if (solver.info() != Eigen::Success)
  {
    // fprintf(stderr, "SolverSystem(LDLT): Error in Decomposition!\n");
    return false;
  }

  int info = 0;
  auto Diag = solver.vectorD();
  for (int i = 0; i < Diag.size(); i++)
  {
    if (std::abs(Diag[i]) < 1e-300)
    {
      // fprintf(stderr, " SolveSystem(LDLT): zero found on diagonal ...\n");
      // fprintf(stderr, " d[%d] = %11.4e\n", i, Diag[i]);
      return false;
    }

    if (Diag[i] < -1e-300)
    {
      // fprintf(stderr, " SolveSystem(LDLT): negative number found on diagonal ...\n");
      // fprintf(stderr, " d[%d] = %11.4e\n", i, Diag[i]);
      info--;
      return false;
    }
  }

  if (info < 0)
  {
    // fprintf(stderr, "Stiffness Matrix is not positive definite: %d negative elements\n", info);
    // fprintf(stderr, "found on decomp diagonal of K.\n");
    // fprintf(stderr, "The stucture may have mechanism and thus not stable in general\n");
    // fprintf(stderr, "Please Make sure that all six\n");
    // fprintf(stderr, "rigid body translations are restrained!\n");
    return false;
  }

  x = solver.solve(b);
  if (solver.info() != Eigen::Success)
  {
    // fprintf(stderr, "SolverSystem(LDLT): Error in Solving!\n");
    return false;
  }

  if(timing_)
  {
    solve_timer_.Stop();
  }
  return true;
}

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
