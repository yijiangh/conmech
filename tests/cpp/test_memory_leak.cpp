#define _CRTDBG_MAP_ALLOC
// report filename / linenumber
// https://docs.microsoft.com/en-us/visualstudio/debugger/finding-memory-leaks-using-the-crt-library?view=vs-2019

#include <string>
#include <iostream>
#include "catch/catch.hpp"
#include <ctime>

// #include <Eigen/Dense>
#include <stiffness_checker/Frame.h>
#include <stiffness_checker/Stiffness.h>

#ifdef _MSC_VER
      #include <direct.h>
      #define GetCurrentDir _getcwd
      const std::string PathSep = "\\";
#else
     #include <unistd.h>
     #define GetCurrentDir getcwd
      const std::string PathSep = "/";
#endif

namespace
{

bool repTestStiffness(const std::string file_path, bool solve_exist_id, bool verbose, int iter, bool re_init)
{
  // https://github.com/ethz-asl/eigen_catkin/wiki/Eigen-Memory-Issues
  assert(iter > 0);
  using namespace conmech::stiffness_checker;

  std::cout << "Testing on file: " << file_path << std::endl;
  // Stiffness sf(file_path, verbose);

  // TODO: randomize element ids
  // partial structure ids
  // this one is for topopt-100
  std::vector<int> exist_e_ids{
    125, 126, 115, 122, 111, 108, 23, 22, 98, 75, 64, 34, 61, 65, 59, 60, 39, 36, 44, 67};

  // std::clock_t c_start = std::clock();

  // for (int i = 0; i < iter; i++) {
  //   std::cout << "==== iter #" << i << std::endl;

  //   if (i == 0 || re_init) {
  //     sf = Stiffness(file_path, verbose);
  //     bool include_sw = true;
  //     sf.setSelfWeightNodalLoad(include_sw);
  //   }

  //   bool success = false;
  //   if (!solve_exist_id) {
  //     std::cout << "full solve." << std::endl;
  //     success = sf.solve();
  //   }
  //   else {
  //     std::cout << "partial solve." << std::endl;
  //     success = sf.solve(exist_e_ids);
  //   }
  //   std::cout << "check result: " << success << std::endl;

  //   // TODO: compliance non-zero value check

  //   // Eigen::MatrixXd nD, fR, eR;
  //   // sf.getSolvedResults(nD, fR, eR, success);
  //   // std::cout << "get solved result: criteria pass : " << success << std::endl;
  //   // std::cout << "nodal disp: \n" << nD << std::endl;
  //   // std::cout << "fixities reaction: \n" << fR << std::endl;
  //   // std::cout << "element reaction: \n" << eR << std::endl;

  //   double trans_tol = sf.getTransTol();
  //   double rot_tol = sf.getRotTol();

  //   double max_trans, max_rot;
  //   int max_trans_vid, max_rot_vid;
  //   sf.getMaxNodalDeformation(max_trans, max_rot, max_trans_vid, max_rot_vid);
  //   std::cout << "max nodal deformation: Trans: " << max_trans << " / tol " << trans_tol << ", at node #" << max_trans_vid << std::endl;
  //   std::cout << "max nodal deformation: Rot: " << max_rot << " / tol " << rot_tol << ", at node #" << max_rot_vid << std::endl;

  //   double compliance;
  //   sf.getSolvedCompliance(compliance);
  //   std::cout << "compliance: " << compliance << std::endl;
  // }

  std::cout << "=====================" << std::endl;
  std::cout << "Solve using existing ids: " << solve_exist_id << std::endl;
  std::cout << "Verbose: " << verbose << std::endl;
  std::cout << "Iter: " << iter << std::endl;
  std::cout << "Re-init: " << re_init << std::endl;

  // std::clock_t c_end = std::clock();
  // std::cout << "***** Test run " << iter << " times, avg run time: "
  // << 1e3 * (c_end - c_start) / iter / CLOCKS_PER_SEC << " ms" << std::endl;

  return true;
}
} // util anon ns


std::string current_working_directory()
{
    char* dir_here = GetCurrentDir( 0, 0 ) ; // **** microsoft specific ****
    // char dir_here[129] = {0};
    // GetCurrentDir(dir_here, 128);
    std::string working_directory(dir_here) ;
    std::free(dir_here) ;
    return working_directory ;
}

TEST_CASE("frame init memory leak check", "memory_check") {
  std::string test_file_path = "C:\\Users\\yijiangh\\Documents\\pb_ws\\conmech\\tests\\assembly_instances\\extrusion\\four-frame.json";
  conmech::stiffness_checker::Frame frame;
  frame.loadFromJson(test_file_path);
}

// TEST_CASE( "repetitive random solve memory leak check", "memory_check" ) {
//     auto cwd = current_working_directory();
//     auto test_data_dir = cwd + PathSep + ".." + PathSep + "assembly_instances" + PathSep + "extrusion";

//     // pass in test_data file path from argument:
//     // https://stackoverflow.com/questions/25211110/find-external-test-file-for-unit-test-by-relative-path-c-cmake-guest
// ! this one even better: https://github.com/libigl/libigl/blob/master/tests/CMakeLists.txt#L33

//     std::string test_file_path = "C:\\Users\\yijiangh\\Documents\\pb_ws\\conmech\\tests\\assembly_instances\\extrusion\\four-frame.json";
//     // auto test_file_path = test_data_dir + PathSep + "fertility.json";
//     // auto test_file_path = test_data_dir + PathSep + "topopt-100_S1_03-14-2019_w_layer.json";

//     bool solve_exist_id = false;
//     bool verbose = false;
//     int run_iter = 1;
//     bool reinit = false;

//     REQUIRE( repTestStiffness(test_file_path, solve_exist_id, verbose, run_iter, reinit) == true );
// }