#pragma once

#include "stiffness_checker/Material.h"
#include <nlohmann/json.hpp>
#include <vector>
#include <Eigen/Dense>

namespace conmech
{
namespace stiffness_checker
{
  /**
   * @brief parse frame data from a given json entry
   * 
   * @param[in] json_data 
   * @param[out] V : #V x 3 matrix of vertex coordinates
   * @param[out] E : #E x 2 matrix of indices of end points into V
   * @param[out] Fixities  : #FixedV x 7 (4) matrix of fixities spec, [node_id, bool fix or not]
   * @param[out] materials : #E std vector of materials
   * @return true 
   * @return false 
   */
  bool parseFrameJson(const nlohmann::json& json_data, 
                      Eigen::MatrixXd& V, Eigen::MatrixXi& E, Eigen::MatrixXi& Fixities, 
                      std::vector<conmech::material::Material>& materials);

  bool parseFrameJson(const std::string& file_path, 
                      Eigen::MatrixXd& V, Eigen::MatrixXi& E, Eigen::MatrixXi& Fixities, 
                      std::vector<conmech::material::Material>& materials);

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
