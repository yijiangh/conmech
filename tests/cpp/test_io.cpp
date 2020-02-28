#include <catch2/catch.hpp>
#include <fstream>
#include <nlohmann/json.hpp>
#include <stiffness_checker/Material.h>

TEST_CASE("material deserialed from a json file", "conmech_material") {
  using namespace conmech::material;
  std::string test_frame_path = "C:\\Users\\yijiangh\\Documents\\pb_ws\\conmech\\tests\\assembly_instances\\extrusion\\four-frame.json";
  std::ifstream is(test_frame_path);
  if (!is.is_open()) {
      throw std::runtime_error("Couldn't open frame: " + test_frame_path);
  }
  nlohmann::json config;
  is >> config;

  Material m;
  parseMaterialPropertiesJson(config["material_properties"], m);

  // these readings are expected to be exact
  // kN/cm2 to kN/m2
  double pressure_unit = 1e4;
  REQUIRE(m.youngs_modulus_ == pressure_unit * 350);
  REQUIRE(m.shear_modulus_ == pressure_unit * 240);
  REQUIRE(std::abs(computePoissonRatio(m.youngs_modulus_, m.shear_modulus_) - m.poisson_ratio_) < 1e-6);

  // kN/m3
  double density_unit = 1;
  REQUIRE(m.density_ == density_unit * 12.2582);

  // cm^2 -> m^2
  double area_unit = 1e-4;
  REQUIRE(m.cross_sec_area_ == area_unit * 0.07068583470577035);

  // cm^4 -> m^4
  double inertia_unit = 1e-8;
  REQUIRE(m.Jx_ == inertia_unit * 0.0007952156404399163);
  REQUIRE(m.Iy_ == inertia_unit * 0.00039760782021995816);
  REQUIRE(m.Iz_ == inertia_unit * 0.00039760782021995816);
}