/*
* ==========================================================================
*		This file is part of the implementation of
*
*		<FrameFab: Robotic Fabrication of Frame Shapes>
*		Yijiang Huang, Juyong Zhang, Xin Hu, Guoxian Song, Zhongyuan Liu, Lei Yu, Ligang Liu
*		In ACM Transactions on Graphics (Proc. SIGGRAPH Asia 2016)
----------------------------------------------------------------------------
*		class:	Stiffness
*
*		Description:
*
*		Version:  2.0
*		Created:  Mar/23/2016
*		Updated:  Sep/2018
*
*		Author:  Xin Hu, Yijiang Huang, Guoxian Song
*		Company:  GCL@USTC
*		Citation:	Some part of this module is modified from frame3dd.c
*					stiffness matrix construction submodule.
*			Title:			Frame3dd source code
*							Static and dynamic structural analysis of 2D and 3D frames and trusses with
*							elastic and geometric stiffness.
*			Author:			Henri P. Gavin
*			Code Version:	20140514+
*			Availability:	http://frame3dd.sourceforge.net/
----------------------------------------------------------------------------
*		Copyright (C) 2016  Yijiang Huang, Xin Hu, Guoxian Song, Juyong Zhang
*		and Ligang Liu.
*
*		FrameFab is free software: you can redistribute it and/or modify
*		it under the terms of the GNU General Public License as published by
*		the Free Software Foundation, either version 3 of the License, or
*		(at your option) any later version.
*
*		FrameFab is distributed in the hope that it will be useful,
*		but WITHOUT ANY WARRANTY; without even the implied warranty of
*		MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*		GNU General Public License for more details.
*
*		You should have received a copy of the GNU General Public License
*		along with FrameFab.  If not, see <http://www.gnu.org/licenses/>.
* ==========================================================================
*/
#pragma once

#include "stiffness_checker/WireFrame.hpp"
#include "stiffness_checker/DualGraph.hpp"
#include "stiffness_checker/FiberPrintPARM.hpp"
#include "stiffness_checker/CoordTrans.hpp"
#include "stiffness_checker/GCommon.hpp"
#include "stiffness_checker/StiffnessIO.hpp"
#include "stiffness_checker/StiffnessSolver.hpp"
//#include "stiffness_checker/IllCondDetector.hpp"

namespace conmech
{
namespace stiffness_checker
{

class Stiffness
{
// private:
//  typedef Eigen::SparseMatrix<double> SpMat;
//  typedef Eigen::MatrixXd MX;
//  typedef Eigen::VectorXd VX;
//  typedef Eigen::VectorXi VXi;
//  typedef Eigen::MatrixXi MXi;
//  typedef trimesh::point point;

 public:
  Stiffness();
  Stiffness(DualGraph *ptr_dualgraph, FiberPrintPARM *ptr_parm,
			bool terminal_output = false);
  ~Stiffness();

 public:
  void Init();
  void CreateFe();
  void CreateF(VX *ptr_x = NULL);
  void CreateElasticK();
  void CreateGlobalK(VX *ptr_x = NULL);

  /* calculate D using LDLT */
  bool CalculateD(
		  VX &D,
		  VX *ptr_x = NULL,
		  bool cond_num = false,
		  int file_id = 0, string file_name = ""
  );

  /* calculate D using ConjugateGradient by Eigen */
  bool CalculateD(
		  VX &D,
		  VX &D0,                        // D0 is the last result
		  VX *ptr_x = NULL,
		  bool cond_num = false,
		  int file_id = 0, string file_name = ""
  );

  /* Check condition number */
  bool CheckIllCondition(IllCondDetector &stiff_inspector);
  bool CheckError(IllCondDetector &stiff_inspector, VX &D);

  /* Data I/O */
  SpMat *WeightedK()
  {
	  assert(&K_);
	  return &K_;
  }
  VX *WeightedF()
  {
	  assert(&F_);
	  return &F_;
  }

  MX eKe(int ei);            // ei: orig e id
  MX eKv(int ei);            // ei: orig e id
  VX Fe(int ei);                // ei: orig e id

  void PrintOutTimer();

 public:
  DualGraph *ptr_dualgraph_;
  FiberPrintPARM *ptr_parm_;

  StiffnessIO stiff_io_;
  StiffnessSolver stiff_solver_;

  CoordTrans transf_;

  SpMat K_;                        // x-Weighted global stiffness matrix, 6n*6n
  vector <MX> eK_;                    // elastic K, indexed by dual id
  VX F_;
  vector <VX> Fe_;

  int Ns_;

  // TODO: should be a list linked to elements
  double radius_; // radius of frame

  // TODO: do we need this?
  double nr_; // radius of node

  double density_;

  // gravitational acceleration
  double g_;

  // shear modulus
  double G_;

  // young's modulus
  double E_;

  // poisson ratio
  double v_;

  bool shear_;                    // 1 : shear deformation taken into consideration; 0 : not

  Timer create_fe_;
  Timer create_f_;
  Timer create_ek_;
  Timer create_k_;
  Timer check_ill_;
  Timer check_error_;

  bool terminal_output_;
};

} // namespace stiffness_checker
} // namespace conmech