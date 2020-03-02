#include <stiffness_checker/Material.h>
#include <stiffness_checker/StiffnessIO.h>
#include "test_common.h"
#include <catch2/catch.hpp>
#include <fstream>
#include <nlohmann/json.hpp>
#include <Eigen/Dense>

namespace
{
void assert_eq_PLA(const conmech::material::Material& m)
{
  using namespace conmech_testing;
  using namespace conmech::material;
  double eps = conmech::EPSILON;

  // these readings are expected to be exact
  // kN/cm2 to kN/m2
  double pressure_unit = 1e4;
  assert_near(m.youngs_modulus_, pressure_unit * 350, eps);
  assert_near(m.shear_modulus_, pressure_unit * 240, eps);
  assert_near(computePoissonRatio(m.youngs_modulus_, m.shear_modulus_), m.poisson_ratio_, eps);

  // kN/m3
  double density_unit = 1;
  assert_near(m.density_, density_unit * 12.2582, eps);

  // cm^2 -> m^2
  double area_unit = 1e-4;
  assert_near(m.cross_sec_area_, area_unit * 0.07068583470577035, eps);

  // cm^4 -> m^4
  double inertia_unit = 1e-8;
  assert_near(m.Jx_, inertia_unit * 0.0007952156404399163, eps);
  assert_near(m.Iy_, inertia_unit * 0.00039760782021995816, eps);
  assert_near(m.Iz_, inertia_unit * 0.00039760782021995816, eps);
}
}

TEST_CASE("material deserialed from a json file", "[io]") 
{
  using namespace conmech_testing;
  using namespace conmech::material;
  std::string test_frame_path = conmech_testing::data_path("four-frame.json");

  std::ifstream is(test_frame_path);
  if (!is.is_open()) {
      throw std::runtime_error("Couldn't open frame: " + test_frame_path);
  }
  nlohmann::json config;
  is >> config;

  Material m;
  parseMaterialPropertiesJson(config["material_properties"], m);
  assert_eq_PLA(m);

  m.setFromJson(config["material_properties"]);
  assert_eq_PLA(m);
}

TEST_CASE("frame data deserialed from a json file", "[io]") 
{
  using namespace conmech_testing;
  using namespace conmech::material;
  using namespace conmech::stiffness_checker;
  std::string test_frame_path = conmech_testing::data_path("four-frame.json");

  Eigen::MatrixXd V;
  Eigen::MatrixXi E;
  Eigen::MatrixXi Fixities;
  std::vector<Material> mats;
  double eps = conmech::EPSILON;

  parseFrameJson(test_frame_path, V, E, Fixities, mats);

  Eigen::MatrixXd realV(5,3);
  realV << 0, -20, -10,
           0, 20, -10,
           0, 20, 0,
           0, -20, 0,
           0, 0, 20;
  realV *= 0.001;
  assert_near_m(V, realV, eps);

  Eigen::MatrixXi realE(4,2);
  realE << 0, 3,
           1, 2,
           3, 4,
           2, 4;
  assert_near_m(E, realE, eps);

  Eigen::MatrixXi realFixities(2, 7);
  realFixities << 0, 1, 1, 1, 1, 1, 1,
                  1, 1, 1, 1, 1, 1, 1;
  assert_near_m(Fixities, realFixities, eps);

  for (const auto& m : mats)
  {
    assert_eq_PLA(m);
  }
}