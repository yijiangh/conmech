//
// Created by yijiangh on 9/26/18.
//

#ifndef CONMECH_STIFFNESSCHECKER_H
#define CONMECH_STIFFNESSCHECKER_H

#include <stiffness_checker/Stiffness.hpp>

namespace conmech
{
namespace stiffness_checker
{

class StiffnessChecker : public Stiffness
{
 public:
  StiffnessChecker(const std::string &json_file_path, bool verbose = false);
  ~StiffnessChecker();

  bool checkDeformation(const std::vector<int> &existing_element_ids);
  bool checkDeformation(const std::vector<int> &existing_element_ids, std::vector<double> &nodes_deformation);
};

} // ns stiffness_checker
} // ns conmech

#endif //CONMECH_STIFFNESSCHECKER_H