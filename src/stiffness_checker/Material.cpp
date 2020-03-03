#include "stiffness_checker/Material.h"
#include "stiffness_checker/SharedConst.h"
#include <stdexcept>
#include <iostream>
#include <fstream>

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

double convertDensityScale(const std::string &unit) 
{
  // TODO: this is super outdated! Change this!
  if ("kN/m^3" == unit || "kN/m3" == unit)
  {
    // ! default unit inside stiffness_checker
    return 1.0;
  }
  return 1.0;
}
} // anon util ns

namespace conmech{
namespace material{

Material::Material(const std::string& file_path)
{
  setFromFile(file_path);
}

Material::Material(const nlohmann::json& entry)
{
  setFromJson(entry);
}

// https://github.com/jpanetta/MeshFEM/blob/master/src/lib/MeshFEM/Materials.cc
void Material::setFromJson(const nlohmann::json &entry) 
{
  double unit_conversion;

  // TODO: check unit
  // kN/cm^2 -> kN/m^2
  unit_conversion = 1e4;
  if (!entry.contains("youngs_modulus") || !entry.contains("youngs_modulus_unit"))
  {
    throw std::runtime_error("Young\'s modulus property value not specified!");
  }
  youngs_modulus_ = unit_conversion * double(entry["youngs_modulus"]);

  // TODO: add tensile strength
  // kN/cm^2 -> kN/m^2
  // unit_conversion = 1e4;
  // frame_parm.tensile_yeild_stress_ = unit_conversion * entry["tensile_yeild_stress"].GetDouble();

  if (entry.contains("shear_modulus")) 
  {
    unit_conversion = 1.0e4;
    shear_modulus_ = unit_conversion * double(entry["shear_modulus"]);
    poisson_ratio_ = computePoissonRatio(youngs_modulus_, shear_modulus_);
    if (entry.contains("poisson_ratio")) 
    {
      double tmp_poisson_ratio = double(entry["poisson_ratio"]);
      // TODO: threhold made shared const
      if (std::abs(tmp_poisson_ratio - poisson_ratio_) > 1e-4) {
          std::cerr << "Warning: shear modulus and poission ratio not compatible! Diff: " << 
          std::abs(tmp_poisson_ratio - poisson_ratio_) << ", given: " << tmp_poisson_ratio 
          << ", computed: " << poisson_ratio_ << ", check formula between E, G, mu. Given E & G are used." << std::endl;
      }
    }
  } else 
  {
    if (!entry.contains("poisson_ratio"))
    {
      throw std::runtime_error("Both shear modulus and poisson_ratio not specified!");
    }
    poisson_ratio_ = double(entry["poisson_ratio"]);
    shear_modulus_ = computeShearModulus(youngs_modulus_, poisson_ratio_);
  }

  // kN/m^3
  if (!entry.contains("density_unit") || !entry.contains("density"))
  {
    throw std::runtime_error("Density not specified!");
  }
  unit_conversion = convertDensityScale(entry["density_unit"]);
  density_ = unit_conversion * double(entry["density"]);

  // cm^2 -> m^2
  unit_conversion = 1.0e-4;
  if(!entry.contains("cross_sec_area")) 
  {
    throw std::runtime_error("Cross section area property value not specified!");
  }
  cross_sec_area_ = unit_conversion * double(entry["cross_sec_area"]);

  // cm^4 -> m^4
  if (!entry.contains("Jx") || 
      !entry.contains("Iy") ||
      !entry.contains("Iz"))
  {
    throw std::runtime_error("Jx, Iy or Iz not specified!");
  }
  unit_conversion = 1.0e-8;
  Jx_ = unit_conversion * double(entry["Jx"]);
  Iy_ = unit_conversion * double(entry["Iy"]);
  Iz_ = unit_conversion * double(entry["Iz"]);
}

void Material::setFromFile(const std::string &materialPath) {
    std::ifstream is(materialPath);
    if (!is.is_open()) {
        throw std::runtime_error("Couldn't open material: " + materialPath);
    }
    nlohmann::json entry;
    is >> entry;
    setFromJson(entry);
}

nlohmann::json Material::getJson() const {
  nlohmann::json entry;
  using namespace conmech;

  // TODO: assign material name
  entry["material_name"] = "exported material from conmech";
  entry["youngs_modulus"] = youngs_modulus_;
  entry["youngs_modulus_unit"] = PRESSURE_UNIT;
  entry["shear_modulus"] = shear_modulus_;
  entry["shear_modulus_unit"] = PRESSURE_UNIT;
  // entry["tensile_yeild_stress": 3.6,
  // entry["tensile_yeild_stress_unit": "kN/cm2",
  entry["density"] = density_;
  entry["density_unit"] = DENSITY_UNIT;
  entry["poisson_ratio"] = poisson_ratio_;

  entry["cross_sec_area"] = cross_sec_area_;
  entry["cross_sec_area_unit"] = AREA_UNIT;
  entry["Jx"] = Jx_;
  entry["Jx_unit"] = AREA_INERTIA_UNIT;
  entry["Iy"] = Iy_;
  entry["Iy_unit"] = AREA_INERTIA_UNIT;
  entry["Iz"] = Iz_;
  entry["Iz_unit"] = AREA_INERTIA_UNIT;

  return entry;
}

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
  mat.setFromJson(entry);
}

} // end ns material
} // end ns conmech