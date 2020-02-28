#include <stiffness_checker/SharedConst.h>
#include <Eigen/Dense>
#include <cstdlib>

namespace conmech_testing
{

template <
  typename DerivedA,
  typename DerivedB>
bool almost_equal(const DerivedA& a, const DerivedB& b, const double& zero=conmech::EPSILON)
{
  return std::abs(double(a) - double(b)) < zero;
}

template <
  typename DerivedA,
  typename DerivedB>
bool almost_equal_m(const Eigen::MatrixBase<DerivedA> & A,
                    const Eigen::MatrixBase<DerivedB> & B,
                    const double& zero=conmech::EPSILON)
{
  ASSERT (A.rows() == B.rows(), "#A:" + std::to_string(A.rows()) + ", #B:" + std::to_string(B.rows()));
  ASSERT (A.cols() == B.cols(), "#A:" + std::to_string(A.rows()) + ", #B:" + std::to_string(B.rows()));
  for(int i=0; i<A.rows(); i++)
  {
    for(int j=0; j<A.cols(); j++)
    {
      if (!almost_equal(A(i,j), B(i,j), zero))
      {
        return false;
      }
    }
  }
  return true;
}

} // ns conmech_testing