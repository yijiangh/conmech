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

#include "stiffness_checker/Frame.h"
#include "stiffness_checker/DualGraph.hpp"
#include "stiffness_checker/StiffnessParm.hpp"
#include "stiffness_checker/CoordTrans.hpp"
#include "stiffness_checker/GCommon.hpp"
#include "stiffness_checker/StiffnessSolver.hpp"
#include "stiffness_checker/Vec.hpp"
//#include "stiffness_checker/IllCondDetector.hpp"

namespace conmech
{
namespace stiffness_checker
{

class Stiffness
{
 public:
  Stiffness();
  Stiffness(DualGraph *ptr_dualgraph, const StiffnessParm& parm,
            bool terminal_output = false);
  ~Stiffness();

 public:
  void Init();

  /* calculate D using LDLT */
  //
  bool CalculateD(
      Eigen::VectorXd &D,
      Eigen::VectorXd *ptr_x = NULL,
      bool cond_num = false);

  /* Check condition number */
//  bool CheckIllCondition(IllCondDetector &stiff_inspector);
//  bool CheckError(IllCondDetector &stiff_inspector, Eigen::VectorXd &D);

  // print out timing result on console
  void PrintOutTimer();

 private:
  void CreateFe();
  void CreateF(Eigen::VectorXd *ptr_x = NULL);
  void CreateElasticK();
  void CreateGlobalK(Eigen::VectorXd *ptr_x = NULL);

 protected:
  DualGraph* ptr_dualgraph_;
  StiffnessParm parm_;

  StiffnessSolver stiff_solver_;

  CoordTrans transf_;

  // x-Weighted global stiffness matrix, 6n*6n
  Eigen::SparseMatrix<double> K_;

  // elastic K, indexed by dual id
  std::vector<Eigen::MatrixXd> eK_;

  Eigen::VectorXd F_;
  std::vector<Eigen::VectorXd> Fe_;

  // number of wf nodes that is not fixed
  int num_free_wf_nodes;

  // TODO: do we need this?
  double nr_; // radius of node

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