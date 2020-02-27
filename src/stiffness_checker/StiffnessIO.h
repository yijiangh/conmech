#pragma once

#include <vector>
#include <Eigen/Dense>
#include "stiffness_checker/Material.h"

namespace conmech
{
namespace stiffness_checker
{
  bool parseFrameJson(Eigen::MatrixXd& V, Eigen::MatrixXi& E, 
                      Eigen::VectorXi& Fixities, std::vector<conmech::material::ConstantMaterial>& materials);

  /**
   * Parse material properties from a json file.
   * @param file_path target json file path
   * @param[out] a list of frame_parm
   * @return boolean flag for success
   */
  bool parseMaterialPropertiesJson(const std::string &file_path, std::vector<conmech::material::ConstantMaterial> &frame_parms);
  
  /**
   *
   * @param file_path target load case json file path
   * @param[out] Load a ((num_loaded_node) x (full_node_dof + 1)) Eigen double matrix
   * Here full_node_dof = 3 (dim 2) or 6 (dim 3)
   * @param[out] include_sw boolean flag for whether include self-weight or not
   * @return boolean flag of success
   */
  bool parseLoadCaseJson(const std::string &file_path, Eigen::MatrixXd& Load, bool& include_sw);

  /**
   * @brief write analysis output to a 
   * 
   * @param V : #V x dim matrix of vertex coordinates
   * @param E : #E x 2 matrix of indices of end points into V
   * @param node_displ 
   * @param fixities_reaction 
   * @param element_reaction 
   * @param file_path 
   * @return true 
   * @return false 
   */
  bool write_output_json(const Eigen::MatrixXd& V, const Eigen::MatrixXi& E,
                         const Eigen::MatrixXd &node_displ,
                         const Eigen::MatrixXd &fixities_reaction,
                         const Eigen::MatrixXd &element_reaction,
                         const std::string& file_path);

} // namespace stiffness_checker
} // namespace conmech
