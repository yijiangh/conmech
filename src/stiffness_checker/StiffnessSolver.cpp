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


////////////////////////////////////////////////////////////////////////////////
// File: crout_pivot.c                                                        //
// Routines:                                                                  //
//    Crout_LU_Decomposition_with_Pivoting                                    //
//    Crout_LU_with_Pivoting_Solve                                            //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  int Crout_LU_Decomposition_with_Pivoting(double *A, int pivot[], int n)   //
//                                                                            //
//  Description:                                                              //
//     This routine uses Crout's method to decompose a row interchanged       //
//     version of the n x n matrix A into a lower triangular matrix L and a   //
//     unit upper triangular matrix U such that A = LU.                       //
//     The matrices L and U replace the matrix A so that the original matrix  //
//     A is destroyed.                                                        //
//     Note!  In Crout's method the diagonal elements of U are 1 and are      //
//            not stored.                                                     //
//     Note!  The determinant of A is the product of the diagonal elements    //
//            of L.  (det A = det L * det U = det L).                         //
//     The LU decomposition is convenient when one needs to solve the linear  //
//     equation Ax = B for the vector x while the matrix A is fixed and the   //
//     vector B is varied.  The routine for solving the linear system Ax = B  //
//     after performing the LU decomposition for A is                         //
//                      Crout_LU_with_Pivoting_Solve.                         //
//     (see below).                                                           //
//                                                                            //
//     The Crout method with partial pivoting is: Determine the pivot row and //
//     interchange the current row with the pivot row, then assuming that     //
//     row k is the current row, k = 0, ..., n - 1 evaluate in order the      //
//     the following pair of expressions                                      //
//       L[i][k] = (A[i][k] - (L[i][0]*U[0][k] + . + L[i][k-1]*U[k-1][k]))    //
//                                 for i = k, ... , n-1,                      //
//       U[k][j] = A[k][j] - (L[k][0]*U[0][j] + ... + L[k][k-1]*U[k-1][j])    //
//                                                                  / L[k][k] //
//                                      for j = k+1, ... , n-1.               //
//       The matrix U forms the upper triangular matrix, and the matrix L     //
//       forms the lower triangular matrix.                                   //
//                                                                            //
//  Arguments:                                                                //
//     double *A       Pointer to the first element of the matrix A[n][n].    //
//     int    pivot[]  The i-th element is the pivot row interchanged with    //
//                     row i.                                                 //
//     int     n       The number of rows or columns of the matrix A.         //
//                                                                            //
//  Return Values:                                                            //
//     0  Success                                                             //
//    -1  Failure - The matrix A is singular.                                 //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     double A[N][N];                                                        //
//     int    pivot[N];                                                       //
//                                                                            //
//     (your code to intialize the matrix A)                                  //
//                                                                            //
//     err = Crout_LU_Decomposition_with_Pivoting(&A[0][0], pivot, N);        //
//     if (err < 0) printf(" Matrix A is singular\n");                        //
//     else { printf(" The LU decomposition of A is \n");                     //
//           ...                                                              //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //

#include <math.h>                                     // required for fabs()

int Crout_LU_Decomposition_with_Pivoting(double *A, int pivot[], int n)
{
   int row, i, j, k, p;
   double *p_k, *p_row, *p_col;
   double max;

//         For each row and column, k = 0, ..., n-1,

   for (k = 0, p_k = A; k < n; p_k += n, k++) {
 
//            find the pivot row

      pivot[k] = k;
      max = fabs( *(p_k + k) );
      for (j = k + 1, p_row = p_k + n; j < n; j++, p_row += n) {
         if ( max < fabs(*(p_row + k)) ) {
            max = fabs(*(p_row + k));
            pivot[k] = j;
            p_col = p_row;
         }
      }

//     and if the pivot row differs from the current row, then
//     interchange the two rows.
   
      if (pivot[k] != k)
         for (j = 0; j < n; j++) {
            max = *(p_k + j);
            *(p_k + j) = *(p_col + j);
            *(p_col + j) = max;
         }

//                and if the matrix is singular, return error

      if ( *(p_k + k) == 0.0 ) return -1;

//      otherwise find the upper triangular matrix elements for row k. 
 
      for (j = k+1; j < n; j++) {
         *(p_k + j) /= *(p_k + k);
      }

//            update remaining matrix

      for (i = k+1, p_row = p_k + n; i < n; p_row += n, i++)
         for (j = k+1; j < n; j++)
            *(p_row + j) -= *(p_row + k) * *(p_k + j);

   }
   return 0;
}


////////////////////////////////////////////////////////////////////////////////
//  int Crout_LU_with_Pivoting_Solve(double *LU, double B[], int pivot[],     //
//                                                        double x[], int n)  //
//                                                                            //
//  Description:                                                              //
//     This routine uses Crout's method to solve the linear equation Ax = B.  //
//     This routine is called after the matrix A has been decomposed into a   //
//     product of a lower triangular matrix L and a unit upper triangular     //
//     matrix U without pivoting.  The argument LU is a pointer to the matrix //
//     the superdiagonal part of which is U and the subdiagonal together with //
//     the diagonal part is L. (The diagonal part of U is 1 and is not        //
//     stored.)   The matrix A = LU.                                          //
//     The solution proceeds by solving the linear equation Ly = B for y and  //
//     subsequently solving the linear equation Ux = y for x.                 //
//                                                                            //
//  Arguments:                                                                //
//     double *LU      Pointer to the first element of the matrix whose       //
//                     elements form the lower and upper triangular matrix    //
//                     factors of A.                                          //
//     double *B       Pointer to the column vector, (n x 1) matrix, B.       //
//     int    pivot[]  The i-th element is the pivot row interchanged with    //
//                     row i.                                                 //
//     double *x       Solution to the equation Ax = B.                       //
//     int     n       The number of rows or columns of the matrix LU.        //
//                                                                            //
//  Return Values:                                                            //
//     0  Success                                                             //
//    -1  Failure - The matrix A is singular.                                 //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     double A[N][N], B[N], x[N];                                            //
//     int    pivot[N];                                                       //
//                                                                            //
//     (your code to create matrix A and column vector B)                     //
//     err = Crout_LU_Decomposition_with_Pivoting(&A[0][0], pivot, N);        //
//     if (err < 0) printf(" Matrix A is singular\n");                        //
//     else {                                                                 //
//        err = Crout_LU_with_Pivoting_Solve(&A[0][0], B, pivot, x, n);       //
//        if (err < 0) printf(" Matrix A is singular\n");                     //
//        else printf(" The solution is \n");                                 //
//           ...                                                              //
//     }                                                                      //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
int Crout_LU_with_Pivoting_Solve(double *LU, double B[], int pivot[], 
                                                            double x[], int n)
{
   int i, k;
   double *p_k;
   double dum;

//         Solve the linear equation Lx = B for x, where L is a lower
//         triangular matrix.                                      
   
   for (k = 0, p_k = LU; k < n; p_k += n, k++) {
      if (pivot[k] != k) {dum = B[k]; B[k] = B[pivot[k]]; B[pivot[k]] = dum; }
      x[k] = B[k];
      for (i = 0; i < k; i++) x[k] -= x[i] * *(p_k + i);
      x[k] /= *(p_k + k);
   }

//         Solve the linear equation Ux = y, where y is the solution
//         obtained above of Lx = B and U is an upper triangular matrix.
//         The diagonal part of the upper triangular part of the matrix is
//         assumed to be 1.0.

   for (k = n-1, p_k = LU + n*(n-1); k >= 0; k--, p_k -= n) {
      if (pivot[k] != k) {dum = B[k]; B[k] = B[pivot[k]]; B[pivot[k]] = dum; }
      for (i = k + 1; i < n; i++) x[k] -= x[i] * *(p_k + i);
      if (*(p_k + k) == 0.0) return -1;
   }
  
   return 0;
}