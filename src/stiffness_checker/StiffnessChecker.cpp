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

  this->Init();
}

StiffnessChecker::~StiffnessChecker()
{
  // TODO: this garbage collection mechanism is very easy to crash the code...

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
  Eigen::VectorXd x;

  x.resize(ptr_dualgraph_->SizeOfVertList());
  x.setZero();

  for(int id : existing_element_ids)
  {
    x[id] = 1;
  }

  bool success = CalculateD(D, &x, terminal_output_);

  if(terminal_output_)
  {
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

  return success;
}

bool StiffnessChecker::checkDeformation(
    const std::vector<int> &existing_element_ids, std::vector<double> &nodes_deformation)
{
  std::cout << "not implemented at this moment." << std::endl;
}

} // ns stiffness_checker
} // ns conmech
