#include <string>
#include <iostream>
#include "catch/catch.hpp"

#include <Eigen/Dense>
#include "stiffness_checker/Frame.h"
#include "stiffness_checker/Util.h"
#include "stiffness_checker/Stiffness.h"
#include "stiffness_checker/StiffnessIO.h"

namespace
{
const std::string PathSeparator =
#ifdef _WIN32
"\\";
#else
"/";
#endif

void repTestStiffness(const char* cfile_name,
  bool solve_exist_id, bool verbose, int iter, bool re_init)
{
  assert(iter > 0);
  using namespace conmech::stiffness_checker;

  std::string file_dir =
    "C:\\Users\\harry\\Documents\\pb-construction\\conmech\\src\\bindings\\pyconmech\\test\\assembly_instances\\extrusion\\";

  // set up shape file path
  std::string file_name = cfile_name;
  std::string file_path = file_dir + file_name;

  Stiffness sf(file_path, verbose);
  // partial structure ids
  // this one is for topopt-100
  std::vector<int> exist_e_ids{
    125, 126, 115, 122, 111, 108, 23, 22, 98, 75, 64, 34, 61, 65, 59, 60, 39, 36, 44, 67};

  assert(iter > 0);
  std::clock_t c_start = std::clock();

  for (int i = 0; i < iter; i++) {
    std::cout << "==== iter #" << i << std::endl;

    if (i == 0 || re_init) {
      sf = Stiffness(file_path, verbose);
      // bool include_sw = true;
      // sf.setSelfWeightNodalLoad(include_sw);
    }

    bool success = false;
    if (!solve_exist_id) {
      success = sf.solve();
    }
    else {
      success = sf.solve(exist_e_ids);
    }
    std::cout << "check result: " << success << std::endl;

    // Eigen::MatrixXd nD, fR, eR;
    // sf.getSolvedResults(nD, fR, eR, success);
    // std::cout << "get solved result: criteria pass : " << success << std::endl;
    // std::cout << "nodal disp: \n" << nD << std::endl;
    // std::cout << "fixities reaction: \n" << fR << std::endl;
    // std::cout << "element reaction: \n" << eR << std::endl;

    double trans_tol = sf.getTransTol();
    double rot_tol = sf.getRotTol();

    double max_trans, max_rot;
    int max_trans_vid, max_rot_vid;
    sf.getMaxNodalDeformation(max_trans, max_rot, max_trans_vid, max_rot_vid);
    std::cout << "max nodal deformation: Trans: " << max_trans << " / tol " << trans_tol << ", at node #" << max_trans_vid << std::endl;
    std::cout << "max nodal deformation: Rot: " << max_rot << " / tol " << rot_tol << ", at node #" << max_rot_vid << std::endl;

    double compliance;
    sf.getSolvedCompliance(compliance);
    std::cout << "compliance: " << compliance << std::endl;
  }

  std::cout << "=====================" << std::endl;
  std::cout << "Testing on file: " << cfile_name << std::endl;
  std::cout << "Solve using existing ids: " << solve_exist_id << std::endl;
  std::cout << "Verbose: " << verbose << std::endl;
  std::cout << "Iter: " << iter << std::endl;
  std::cout << "Re-init: " << re_init << std::endl;

  std::clock_t c_end = std::clock();
  std::cout << "***** Test run " << iter << " times, avg run time: "
  << 1e3 * (c_end - c_start) / iter / CLOCKS_PER_SEC << " ms" << std::endl;
}
} // util anon ns

// https://docs.microsoft.com/en-us/visualstudio/debugger/finding-memory-leaks-using-the-crt-library?view=vs-2019
// #if _MSC_VER
// 	void *testWhetherMemoryLeakDetectionWorks = malloc(1);
// #endif

TEST_CASE( "repetitive random solve memory leak check", "memory_check" ) {

    // REQUIRE( Factorial(0) == 1 );

}