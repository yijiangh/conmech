/*
* ==========================================================================
*		This file is part of the implementation of
*
*		<FrameFab: Robotic Fabrication of Frame Shapes>
*		Yijiang Huang, Juyong Zhang, Xin Hu, Guoxian Song, Zhongyuan Liu, Lei Yu, Ligang Liu
*		In ACM Transactions on Graphics (Proc. SIGGRAPH Asia 2016)
----------------------------------------------------------------------------
*		class:	FiberPrintPARM
*
*		Description:	FiberPrintPARM takes charge of computation related parameters configuration.
*
*		Version: 2.0
*		Created: Oct/10/2015
*		Updated: Aug/24/2016
*
*		Author:  Xin Hu, Yijiang Huang, Guoxian Song
*		Company:  GCL@USTC

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

#ifndef FRAMEFAB_TASK_SEQUENCE_PLANNER_FIBERPRIN_PARM_H
#define FRAMEFAB_TASK_SEQUENCE_PLANNER_FIBERPRIN_PARM_H

namespace conmech
{
namespace stiffness_checker
{

class StiffnessParm
{
 public:
  StiffnessParm(
		  double radius = 0.75,
		  double density = 1210 * 1e-12,
		  double g = -9806.3,
		  double youngs_modulus = 3457,
		  double shear_modulus = 1294,
		  double poisson_ratio = 0.335)
		  : radius_(radius), density_(density),
			g_(g),
			youngs_modulus_(youngs_modulus),
			shear_modulus_(shear_modulus),
			poisson_ratio_(poisson_ratio) {}

  ~StiffnessParm(){};

 public:
  // material & environment
  double radius_;
  double density_;
  double g_;
  double youngs_modulus_;
  double shear_modulus_;
  double poisson_ratio_;

};

} // ns stiffness_checker
} // ns conmech

#endif