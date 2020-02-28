#include <stiffness_checker/Material.h>
#include <stiffness_checker/StiffnessIO.h>
#include <catch2/catch.hpp>
#include <fstream>
#include <nlohmann/json.hpp>
#include <Eigen/Dense>
#include "utils.h"

namespace
{
bool equal_PLA_material(const conmech::material::Material& m)
{
  using namespace conmech_testing;
  using namespace conmech::material;
  bool is_eq = true;

  // these readings are expected to be exact
  // kN/cm2 to kN/m2
  double pressure_unit = 1e4;
  is_eq &= almost_equal(m.youngs_modulus_, pressure_unit * 350);
  is_eq &= almost_equal(m.shear_modulus_, pressure_unit * 240);
  is_eq &= almost_equal(computePoissonRatio(m.youngs_modulus_, m.shear_modulus_), m.poisson_ratio_);

  // kN/m3
  double density_unit = 1;
  is_eq &= almost_equal(m.density_, density_unit * 12.2582);

  // cm^2 -> m^2
  double area_unit = 1e-4;
  is_eq &= almost_equal(m.cross_sec_area_, area_unit * 0.07068583470577035);

  // cm^4 -> m^4
  double inertia_unit = 1e-8;
  is_eq &= almost_equal(m.Jx_, inertia_unit * 0.0007952156404399163);
  is_eq &= almost_equal(m.Iy_, inertia_unit * 0.00039760782021995816);
  is_eq &= almost_equal(m.Iz_, inertia_unit * 0.00039760782021995816);
  return is_eq;
}
}

TEST_CASE("material deserialed from a json file", "conmech_io") 
{
  using namespace conmech_testing;
  using namespace conmech::material;
  const std::string test_frame_path = "C:\\Users\\yijiangh\\Documents\\pb_ws\\conmech\\tests\\assembly_instances\\extrusion\\four-frame.json";

  std::ifstream is(test_frame_path);
  if (!is.is_open()) {
      throw std::runtime_error("Couldn't open frame: " + test_frame_path);
  }
  nlohmann::json config;
  is >> config;

  Material m;
  parseMaterialPropertiesJson(config["material_properties"], m);
  REQUIRE(equal_PLA_material(m));
}

TEST_CASE("frame data deserialed from a json file", "conmech_io") 
{
  using namespace conmech_testing;
  using namespace conmech::material;
  using namespace conmech::stiffness_checker;
  const std::string test_frame_path = "C:\\Users\\yijiangh\\Documents\\pb_ws\\conmech\\tests\\assembly_instances\\extrusion\\four-frame.json";

  Eigen::MatrixXd V;
  Eigen::MatrixXi E;
  Eigen::MatrixXi Fixities;
  std::vector<Material> mats;

  parseFrameJson(test_frame_path, V, E, Fixities, mats);

  Eigen::MatrixXd realV(5,3);
  realV << 0, -20, -10,
           0, 20, -10,
           0, 20, 0,
           0, -20, 0,
           0, 0, 20;
  REQUIRE(almost_equal_m(V, realV));

  Eigen::MatrixXi realE(4,2);
  realE << 0, 3,
           1, 2,
           3, 4,
           2, 4;
  REQUIRE(almost_equal_m(E, realE));

  Eigen::MatrixXi realFixities(2, 7);
  realFixities << 0, 1, 1, 1, 1, 1, 1,
                  1, 1, 1, 1, 1, 1, 1;
  REQUIRE(almost_equal_m(Fixities, realFixities));

  for (const auto& m : mats)
  {
    REQUIRE(equal_PLA_material(m));
  }
}