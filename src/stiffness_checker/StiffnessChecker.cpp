//
// Created by yijiangh on 9/26/18.
//

#include <iostream>
#include "stiffness_checker/Stiffness.hpp"
#include "stiffness_checker/StiffnessIO.hpp"

#include "stiffness_checker/StiffnessChecker.hpp"

namespace conmech
{
namespace stiffness_checker
{

StiffnessChecker::StiffnessChecker(const std::string &json_file_path, bool verbose)
{
  terminal_output_ = verbose;
  stiff_solver_.detailed_timing_ = verbose;

  // construct dual graph
  this->ptr_dualgraph_ = new DualGraph();
  parseFrameJson(json_file_path, this->ptr_dualgraph_, this->parm_);

  this->nr_ = 0.0;
  this->shear_ = 0;

  Init();
}

StiffnessChecker::~StiffnessChecker()
{
  if(ptr_dualgraph_ != NULL)
  {
    if(ptr_dualgraph_->ptr_frame_ != NULL)
    {
      delete ptr_dualgraph_->ptr_frame_;
      ptr_dualgraph_->ptr_frame_ = NULL;
    }

    delete ptr_dualgraph_;
    ptr_dualgraph_ = NULL;
  }
}

bool StiffnessChecker::checkDeformation(const std::vector<int> &existing_element_ids)
{
  Eigen::VectorXd D;
  bool success = CalculateD(D, NULL, true, 0, "FiberTest");

  if (success)
  {
    std::cout << "[conmech::stiffness_checker] Maximal node deforamtion (mm): " << D.maxCoeff()
              << ", tolerance (mm): (N/A)" << std::endl;
  }
  else
  {
    std::cout << "[conmech::stiffness_checker] stiffness computation fails." << std::endl;
  }

  PrintOutTimer();
}

bool StiffnessChecker::checkDeformation(
    const std::vector<int> &existing_element_ids, std::vector<double> &nodes_deformation)
{

}

} // ns stiffness_checker
} // ns conmech
