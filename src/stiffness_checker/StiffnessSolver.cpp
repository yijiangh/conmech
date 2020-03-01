#include "stiffness_checker/SharedConst.h"
#include "stiffness_checker/StiffnessSolver.h"
#include <Eigen/LU>
#include <Eigen/SparseCholesky>
#include <iostream>

namespace conmech
{
namespace stiffness_checker
{

bool StiffnessSolver::solveSparseSimplicialLDLT(
    const Eigen::SparseMatrix<double>& A, const Eigen::VectorXd& b, Eigen::VectorXd& x, const bool& verbose)
{
  if(timing_)
  {
    solve_timer_.Start();
  }

  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
  solver.compute(A);

  if (solver.info() == Eigen::NumericalIssue)
  {
    if (verbose)
    {
      std::cerr << "SolverSystem(LDLT): Error in Decomposition!: " << solver.info() << std::endl;
    }
    int info = 0;
    auto Diag = solver.vectorD();
    for (int i = 0; i < Diag.size(); i++)
    {
      if (std::abs(Diag[i]) < 1e-300)
      {
        if (verbose)
        {
          std::cerr << " SolveSystem(LDLT): zero found on diagonal ..." << std::endl;
          std::cerr << " d[" << i << "] = " << Diag[i] << std::endl;
        }
      }
      if (Diag[i] < -1e-300)
      {
        if (verbose)
        {
          std::cerr << " SolveSystem(LDLT): negative number found on diagonal ..." << std::endl;
          std::cerr << " d[" << i << "] = " << Diag[i] << std::endl;
        }
        info--;
      }
    }
    if (info < 0)
    {
      if (verbose)
      {
        std::cerr << "Stiffness Matrix is not positive definite: " << info 
          << " negative elements found on decomp diagonal of K." << std::endl;
        std::cerr << "Matrix size:" << A.rows() << ", " << A.cols() << std::endl;
        std::cerr << "The stucture may have mechanism and thus not stable in general," << std::endl;
        std::cerr << "please Make sure that all six rigid body translations are restrained!" << std::endl;
      }
    }
    return false;
  }

  x = solver.solve(b);
  if (solver.info() != Eigen::Success)
  {
    if (verbose)
    {
      std::cerr << "SolverSystem(LDLT): Error in Solving!" << std::endl;
    }
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