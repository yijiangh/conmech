#pragma once

#include <vector>
#include "stiffness_checker/StiffnessParm.h"
#include "stiffness_checker/Frame.h"

namespace conmech
{
namespace stiffness_checker
{
bool parseMaterialPropertiesJson(const std::string &file_path, StiffnessParm &frame_parm);

void outputPath(const char *fname, char fullpath[], const int len, char *default_outdir, int verbose);

/**
 * Create mesh data of deformed and undeformed mesh, use gnuplot
*/
void createGnuPltStaticMesh(
    const std::string& file_path,
    const Frame& frame,
    const Eigen::MatrixXd& nodal_displ,
    double exagg_static, float scale);

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

/**
 * write input data to a .3dd file
*/
void writeFrame3ddData(
    const std::string &fpath,
    const Frame &frame,
    const StiffnessParm &parm,
    int verbose
);

} // namespace stiffness_checker
} // namespace conmech
