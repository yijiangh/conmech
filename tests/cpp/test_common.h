#pragma once

#include <stiffness_checker/SharedConst.h>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <catch2/catch.hpp>

#include <cctype>
#include <string>
#include <functional>
#include <algorithm>
#include <tuple>

namespace{
const std::string PathSeparator =
#ifdef _WIN32
"\\";
#else
"/";
#endif
}

namespace conmech_testing
{
  template<typename Param, typename Fun>
  void run_test_cases(const std::vector<Param> &params,  Fun test_case)
  {
    for(const auto &p : params)
      test_case(p);
  }

  // TODO: other cases
  inline std::vector<std::string> four_frame()
  {
    return
    {
      "four-frame.json",
    };
  };

  inline std::string data_path(std::string s)
  {
    // set in cmake
    return std::string(CONMECH_DATA_DIR) + "/" + s;
  };

  template <typename DerivedA, typename DerivedB>
  void assert_eq(
    const Eigen::MatrixBase<DerivedA> & A,
    const Eigen::MatrixBase<DerivedB> & B)
  {
    // Sizes should match
    REQUIRE(A.rows() == B.rows());
    REQUIRE(A.cols() == B.cols());
    for(int i = 0;i<A.rows();i++)
    {
      for(int j = 0;j<A.cols();j++)
      {
        // Create an ijv tuple to trick GoogleTest into printing (i,j) so we
        // know where the disagreement is.
        std::tuple<int,int,typename DerivedA::Scalar> Aijv {i,j,A(i,j)};
        std::tuple<int,int,typename DerivedB::Scalar> Bijv {i,j,B(i,j)};
        REQUIRE(Aijv == Bijv);
      }
    }
  }

  template <typename DerivedA, typename DerivedB>
  void assert_neq(
    const Eigen::MatrixBase<DerivedA> & A,
    const Eigen::MatrixBase<DerivedB> & B)
  {
    // Sizes should match
    REQUIRE(A.rows() == B.rows());
    REQUIRE(A.cols() == B.cols());
    bool all_equals = true;
    for(int i = 0;i<A.rows();i++)
    {
      for(int j = 0;j<A.cols();j++)
      {
        if (A(i,j) != B(i,j))
        {
          all_equals = false;
        }
      }
    }
    REQUIRE_FALSE(all_equals);
  }

  template <typename DerivedA, typename DerivedB>
  void assert_eq(
    const Eigen::SparseMatrix<DerivedA> & A,
    const Eigen::SparseMatrix<DerivedB> & B)
  {
    // Sizes should match
    REQUIRE(A.rows() == B.rows());
    REQUIRE(A.cols() == B.cols());
    Eigen::Matrix<long int,Eigen::Dynamic, 1> AI,AJ;
    Eigen::Matrix<long int,Eigen::Dynamic, 1> BI,BJ;
    Eigen::Matrix<DerivedA,Eigen::Dynamic, 1> AV;
    Eigen::Matrix<DerivedB,Eigen::Dynamic, 1> BV;
    // Assumes A and B are in same Major Ordering
    igl::find(A,AI,AJ,AV);
    igl::find(B,BI,BJ,BV);
    // This doesn't generalized to assert_near nicely, and it makes it hard to
    // tell which entries are different:
    assert_eq(AI,BI);
    assert_eq(AJ,BJ);
    assert_eq(AV,BV);
  }

  template <typename DerivedA, typename DerivedB, typename EpsType>
  void assert_near_m(
    const Eigen::MatrixBase<DerivedA> & A,
    const Eigen::MatrixBase<DerivedB> & B,
    const EpsType & eps)
  {
    // Sizes should match
    REQUIRE(A.rows() == B.rows());
    REQUIRE(A.cols() == B.cols());
    for(int i = 0;i<A.rows();i++)
    {
      for(int j = 0;j<A.cols();j++)
      {
        // Create an ijv tuple to trick GoogleTest into printing (i,j) so we
        // know where the disagreement is.
        //
        // Equivalent to ASSERT_NEAR(Aijv,Bijv)

        CAPTURE( i );
        CAPTURE( j );
        {
          // std::tuple<int,int,typename DerivedA::Scalar> Aijv {i,j,A(i,j)};
          // std::tuple<int,int,typename DerivedB::Scalar> Bijv {i,j,B(i,j)+eps};
          REQUIRE(A(i,j) < B(i,j)+eps);
        }
        {
          // std::tuple<int,int,typename DerivedA::Scalar> Aijv {i,j,A(i,j)+eps};
          // std::tuple<int,int,typename DerivedB::Scalar> Bijv {i,j,B(i,j)};
          REQUIRE(A(i,j)+eps > B(i,j));
        }
      }
    }
  }

  template <
    typename DerivedA,
    typename DerivedB,
    typename EpsType>
  void assert_near(const DerivedA& a, 
                   const DerivedB& b, 
                   const EpsType & eps)
  {
    CAPTURE( a );
    CAPTURE( b );
    REQUIRE(std::abs(double(a) - double(b)) < eps);
  }
} // ns conmech_testing