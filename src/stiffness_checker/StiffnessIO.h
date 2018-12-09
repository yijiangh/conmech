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
 * @param file_path target json file path
 * @param[out] Load a ((#loaded_node) x (full_node_dof + 1)) Eigen double matrix
 * Here full_node_dof = 3 (dim 2) or 6 (dim 3)
 * @param[out] include_sw boolean flag for whether include self-weight or not
 * @return boolean flag of success
 */
bool parseLoadCaseJson(const std::string &file_path, Eigen::MatrixXd& Load, bool& include_sw);

/**
 * Create mesh data of deformed and undeformed mesh, use gnuplot
 * GNUPlot references:
 * http://svn.code.sourceforge.net/p/frame3dd/code/trunk/doc/Frame3DD-manual.html
 * https://people.duke.edu/~hpgavin/gnuplot.html
 * @param model_name model name (usually I simply use *.json), for naming the generated .plt and .dat files
 * @param file_path dir path for plt and dat files
 * @param frame
 * @param nodal_displ
 * @param exagg_static
 * @param draw_deform boolean flag for drawing deformed shape
 */
void createGnuPltStaticShape(
    const std::string& model_name,
    const std::string& file_path,
    const Frame& frame,
    const Eigen::MatrixXd& nodal_displ,
    double exagg_static, bool draw_deform);

/**
 * computes cubic deflection functions from end deflections
 * and end rotations.  Saves deflected shapes to a file.
 * These bent shapes are exact for mode-shapes, and for frames
 * loaded at their nodes.
*/
void createGnuPltCubicBentBeam(
    const Frame& frame,
    const Eigen::MatrixXd& nodal_displ,
    std::vector<Eigen::Vector3d>& beam,
    double exagg = 1.0);

} // namespace stiffness_checker
} // namespace conmech
