#include "Material.h"
#include <stdexcept>
#include <iostream>

namespace {
double convertModulusScale(const std::string &unit) {
  // TODO: this is super outdated! Change this!

  // 1 pa = 1 N/m^2 = 1e-4 N/cm^2 = 1e-6 N/mm^2
  // 1 Mpa = 1 N/mm^2

  if ("MPa" == unit)
  {
    return 1;
  }
  if ("GPa" == unit)
  {
    return 1e3;
  }
  if ("kN/cm2" == unit)
  {
    return 1e4;
  }
  if ("kN/m2" == unit)
  {
    // ! default unit inside stiffness_checker
    return 1;
  }
  // std::cout << "WARNING: unrecognized modulus unit in the input json file. "
  //              "Using MPa by default." << std::endl;
  return 1;
}

double convertDensityScale(const std::string &unit) {
  // TODO: this is super outdated! Change this!
  if ("kN/m3" == unit)
  {
    // ! default unit inside stiffness_checker
    return 1.0;
  }
  return 1.0;
}
} // anon util ns

namespace conmech{
namespace material{

// // https://github.com/jpanetta/MeshFEM/blob/master/src/lib/MeshFEM/Materials.cc
// void ConstantMaterial::setFromJson(const nlohmann::json &entry) {
//     std::string type = entry["type"];
//     if (type == "isotropic_material")        { parseIsotropic<_N>(entry, m_E);   }
//     else if (type == "isotropic")            { parseIsotropic<_N>(entry, m_E);   }
//     else if (type == "orthotropic_material") { parseOrthotropic<_N>(entry, m_E); }
//     else if (type == "orthotropic")          { parseOrthotropic<_N>(entry, m_E); }
//     else if (type == "symmetric_material")   { parseAnisotropic<_N>(entry, m_E); }
//     else if (type == "anisotropic")          { parseAnisotropic<_N>(entry, m_E); }
//     else { throw std::runtime_error("Invalid type."); }
// }

// void ConstantMaterial::setFromFile(const std::string &materialPath) {
//     std::ifstream is(materialPath);
//     if (!is.is_open()) {
//         throw std::runtime_error("Couldn't open material " + materialPath);
//     }
//     nlohmann::json entry;
//     is >> entry;
//     setFromJson(entry);
// }

// nlohmann::json ConstantMaterial::getJson() const {
//     nlohmann::json entry;

//     std::array<std::array<Real, flatLen(N)>, flatLen(N)> mat;
//     for (size_t i = 0; i < flatLen(N); ++i) {
//         for (size_t j = 0; j < flatLen(N); ++j) {
//             mat[i][j] = m_E.D(i, j);
//         }
//     }

//     entry["type"] = "anisotropic";
//     entry["material_matrix"] = mat;

//     return entry;
// }

double computeShearModulus(const double& E, const double& poisson_ratio) 
{
   return 0.5 * E / (1 + poisson_ratio); 
}

double computePoissonRatio(const double& E, const double& G)
{ 
  return (E / (2.0*G) - 1); 
}

void parseMaterialPropertiesJson(const nlohmann::json& entry, Material& mat)
{
  double unit_conversion;

  // TODO: check unit
  // kN/cm^2 -> kN/m^2
  unit_conversion = 1e4;
  mat.youngs_modulus_ = unit_conversion * double(entry["youngs_modulus"]);

  unit_conversion = 1e4;
  mat.shear_modulus_ = unit_conversion * double(entry["shear_modulus"]);

  // TODO: add tensile strength
  // kN/cm^2 -> kN/m^2
  // unit_conversion = 1e4;
  // frame_parm.tensile_yeild_stress_ = unit_conversion * entry["tensile_yeild_stress"].GetDouble();

  if (entry.contains("shear_modulus")) {
    unit_conversion = 1.0e4;
    mat.shear_modulus_ = unit_conversion * entry["shear_modulus"];
    mat.poisson_ratio_ = computePoissonRatio(mat.youngs_modulus_, mat.shear_modulus_);
    if (entry.contains("poisson_ratio")) {
      double tmp_poisson_ratio = entry["poisson_ratio"];
      // TODO: threhold made shared const
      if (std::abs(tmp_poisson_ratio - mat.poisson_ratio_) > 1e-5) {
          std::cout << "Warning: shear modulus and poission ratio not compatible! Diff: " << 
          std::abs(tmp_poisson_ratio - mat.poisson_ratio_) << ", given: " << tmp_poisson_ratio 
          << ", computed: " << mat.poisson_ratio_ << ", check formula between E, G, mu. Given E & G are used." << std::endl;
      }
    }
  } else {
    if (!entry.contains("poisson_ratio")){
      throw std::runtime_error("Both shear modulus and poisson_ratio not specified!");
    }
    mat.poisson_ratio_ = entry["poisson_ratio"];
    mat.shear_modulus_ = computeShearModulus(mat.youngs_modulus_, mat.poisson_ratio_);
  }

  // kN/m^3
  if (!entry.contains("density_unit") || !entry.contains("density")){
    throw std::runtime_error("Density not specified!");
  }
  unit_conversion = convertDensityScale(entry["density_unit"]);
  mat.density_ = unit_conversion * double(entry["density"]);

  // cm^2 -> m^2
  unit_conversion = 1.0e-4;
  if(!entry.contains("cross_sec_area")) {
    throw std::runtime_error("Cross section area property value not specified!");
  }
  mat.cross_sec_area_ = unit_conversion * double(entry["cross_sec_area"]);

  // cm^4 -> m^4
  if (!entry.contains("Jx") || 
      !entry.contains("Iy") ||
      !entry.contains("Iz")){
    throw std::runtime_error("Jx, Iy or Iz not specified!");
  }
  unit_conversion = 1.0e-8;
  mat.Jx_ = unit_conversion * double(entry["Jx"]);
  mat.Iy_ = unit_conversion * double(entry["Iy"]);
  mat.Iz_ = unit_conversion * double(entry["Iz"]);
}

} // end ns material
} // end ns conmech