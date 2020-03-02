#pragma once
#include <nlohmann/json.hpp>

namespace conmech
{
namespace material
{
double computeShearModulus(const double& E, const double& poisson_ratio);

double computePoissonRatio(const double& E, const double& G);

struct Material
{
 public:
  // https://github.com/jpanetta/MeshFEM/blob/master/src/lib/MeshFEM/Materials.hh
  Material(){}
  Material(const std::string& file_path);
  Material(const nlohmann::json& config);

  void setFromFile(const std::string& materialFile);
  void setFromJson(const nlohmann::json& config);
  nlohmann::json getJson() const;

  /**
   * unit: kN/m^2
   */
  double youngs_modulus_;

  /**
   * unit: kN/m^2
   */
  double shear_modulus_;

  /**
   * unit: [] (unitless)
   * G = E / 2(1+mu),
   * where G is shear modulus, E is Young's modulus
   * mu is poission ratio
   */
  double poisson_ratio_;

  /**
   * unit: kN/m^2
   */
  double tensile_yeild_stress_;

  /**
   * material density, unit: kN/m^3
   */
  double density_;

  // TODO: this radius should be associated with geometry
  /**
   * raidus of the cross section
   * unit: meter
   */
  // double radius_;
  double cross_sec_area_;

  // torsion constant (around local x axis)
  // for solid circle: (1/2)pi*r^4, unit: m^4
  // see: https://en.wikipedia.org/wiki/Torsion_constant#Circle
  double Jx_;

  // area moment of inertia (bending about local y,z-axis)
  // assuming solid circular area of radius r, unit: m^4
  // see: https://en.wikipedia.org/wiki/List_of_area_moments_of_inertia
  // note this is slender rod of length L and Mass M, spinning around end
  double Iy_;
  double Iz_;
};


void parseMaterialPropertiesJson(const nlohmann::json& entry, Material& M);

} // ns stiffness_checker
} // ns conmech
