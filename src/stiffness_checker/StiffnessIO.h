#pragma once

#include <vector>
#include "stiffness_checker/StiffnessParm.h"
#include "stiffness_checker/Frame.h"

namespace conmech
{
namespace stiffness_checker
{
bool parseMaterialPropertiesJson(const std::string &file_path, StiffnessParm &frame_parm);

/**
 * Create mesh data of deformed and undeformed mesh, use gnuplot
 * GNUPlot references:
 * http://svn.code.sourceforge.net/p/frame3dd/code/trunk/doc/Frame3DD-manual.html
 * https://people.duke.edu/~hpgavin/gnuplot.html
 *
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
