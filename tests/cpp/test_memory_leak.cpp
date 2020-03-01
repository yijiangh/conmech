#define _CRTDBG_MAP_ALLOC
// report filename / linenumber
// https://docs.microsoft.com/en-us/visualstudio/debugger/finding-memory-leaks-using-the-crt-library?view=vs-2019

#include <catch2/catch.hpp>
#include "test_common.h"
#include <stiffness_checker/Stiffness.h>
#include <stiffness_checker/StiffnessIO.h>
// #include <Eigen/Dense>
#include <string>
#include <iostream>
#include <fstream>
#include <ctime>
#include <nlohmann/json.hpp>

namespace
{
bool repTestStiffness(const std::string& test_frame_path, const bool& solve_exist_id, 
  const bool& verbose, const int& iter, const bool& re_init)
{
  assert(iter > 0);
  using namespace conmech::stiffness_checker;
  using namespace conmech_testing;
  using namespace conmech::material;

  std::cout << "Testing on file: " << test_frame_path << std::endl;
  Eigen::MatrixXd V;
  Eigen::MatrixXi E;
  Eigen::MatrixXi Fixities;
  std::vector<Material> mats;
  parseFrameJson(test_frame_path, V, E, Fixities, mats);

  Stiffness sf(V, E, Fixities, mats, verbose);

  // TODO: randomize element ids
  // partial structure ids
  // this one is for topopt-100
  std::vector<int> exist_e_ids{
    125, 126, 115, 122, 111, 108, 23, 22, 98, 75, 64, 34, 61, 65, 59, 60, 39, 36, 44, 67};

  std::clock_t c_start = std::clock();

  for (int i = 0; i < iter; i++) {
    std::cout << "==== iter #" << i << std::endl;

    if (i == 0 || re_init) {
      sf = Stiffness(V, E, Fixities, mats, verbose);
      bool include_sw = true;
      sf.setSelfWeightNodalLoad(include_sw);
    }

    bool success = false;
    if (!solve_exist_id) {
      std::cout << "full solve." << std::endl;
      success = sf.solve();
    }
    else {
      std::cout << "partial solve." << std::endl;
      success = sf.solve(exist_e_ids);
    }
    std::cout << "check result: " << success << std::endl;

    double trans_tol = sf.getTransTol();
    double rot_tol = sf.getRotTol();

    double max_trans, max_rot;
    int max_trans_vid, max_rot_vid;
    sf.getMaxNodalDeformation(max_trans, max_rot, max_trans_vid, max_rot_vid);
    std::cout << "max nodal deformation: Trans: " << max_trans << " / tol " << trans_tol << ", at node #" << max_trans_vid << std::endl;
    std::cout << "max nodal deformation: Rot: " << max_rot << " / tol " << rot_tol << ", at node #" << max_rot_vid << std::endl;

    double compliance;
    sf.getSolvedCompliance(compliance);
    REQUIRE(compliance > 0);
    std::cout << "compliance: " << compliance << std::endl;
  }

  std::cout << "=====================" << std::endl;
  std::cout << "Solve using existing ids: " << solve_exist_id << std::endl;
  std::cout << "Verbose: " << verbose << std::endl;
  std::cout << "Iter: " << iter << std::endl;
  std::cout << "Re-init: " << re_init << std::endl;

  std::clock_t c_end = std::clock();
  std::cout << "***** Test run " << iter << " times, avg run time: "
  << 1e3 * (c_end - c_start) / iter / CLOCKS_PER_SEC << " ms" << std::endl;

  return true;
}
} // util anon ns

TEST_CASE( "repetitive random solve memory leak check", "[memory_check]" ) 
{
  int run_iter = 10;
  bool verbose = false;

  SECTION("topopt-100 solve exist id w/o reinit")
  {
    bool solve_exist_id = true;
    bool reinit = false;
    std::string test_file_path = conmech_testing::data_path("topopt-100.json");
    REQUIRE(repTestStiffness(test_file_path, solve_exist_id, verbose, run_iter, reinit));
  }
  SECTION("topopt-100 solve exist id w reinit")
  {
    bool solve_exist_id = true;
    bool reinit = true;
    std::string test_file_path = conmech_testing::data_path("topopt-100.json");
    REQUIRE(repTestStiffness(test_file_path, solve_exist_id, verbose, run_iter, reinit));
  }
  SECTION("topopt-100 solve all w/o reinit")
  {
    bool solve_exist_id = false;
    bool reinit = false;
    std::string test_file_path = conmech_testing::data_path("topopt-100.json");
    REQUIRE(repTestStiffness(test_file_path, solve_exist_id, verbose, run_iter, reinit));
  }
  SECTION("topopt-100 solve all w reinit")
  {
    // void *testWhetherMemoryLeakDetectionWorks = malloc(1);
    bool solve_exist_id = false;
    bool reinit = true;
    std::string test_file_path = conmech_testing::data_path("topopt-100.json");
    REQUIRE(repTestStiffness(test_file_path, solve_exist_id, verbose, run_iter, reinit));
  }
}