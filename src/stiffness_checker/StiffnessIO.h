#pragma once

#include <vector>
#include "stiffness_checker/StiffnessParm.h"
#include "stiffness_checker/Frame.h"

namespace conmech
{
namespace stiffness_checker
{
/**
 * Parse material properties from a json file.
 * @param file_path target json file path
 * @param[out] frame_parm
 * @return boolean flag for success
 */
bool parseMaterialPropertiesJson(const std::string &file_path, StiffnessParm &frame_parm);

/**
 *
 * @param file_path target load case json file path
 * @param[out] Load a ((num_loaded_node) x (full_node_dof + 1)) Eigen double matrix
 * Here full_node_dof = 3 (dim 2) or 6 (dim 3)
 * @param[out] include_sw boolean flag for whether include self-weight or not
 * @return boolean flag of success
 */
bool parseLoadCaseJson(const std::string &file_path, Eigen::MatrixXd& Load, bool& include_sw);

bool write_output_json(const Frame& frame,
    const Eigen::MatrixXd &node_displ,
    const Eigen::MatrixXd &fixities_reaction,
    const Eigen::MatrixXd &element_reaction,
    const std::string& file_path);

} // namespace stiffness_checker
} // namespace conmech
