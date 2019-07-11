#include <string>
#include <iostream>
#include <sstream>
#include <chrono>
#include <ctime>

#include <Eigen/Dense>

#include "stiffness_checker/Frame.h"
#include "stiffness_checker/Util.h"
#include "stiffness_checker/Stiffness.h"
#include "stiffness_checker/StiffnessIO.h"

namespace
{
void testLocalGlobalTransf()
{
  auto end_u = Eigen::Vector3d(1,1,1);
  auto end_v = Eigen::Vector3d(1,1,1) + Eigen::Vector3d(3,5,8);
  Eigen::Matrix3d m;

  double pi = atan(1)*4;

  using namespace conmech::stiffness_checker;
  getGlobal2LocalRotationMatrix(end_u, end_v, m, 70 * pi/180);

  std::cout << m << std::endl;
}

void repTestStiffness(const char* cfile_name,
  bool solve_exist_id, bool verbose, int iter, bool re_init)
{
  assert(iter > 0);
  using namespace conmech::stiffness_checker;

  std::string file_dir =
    "/Users/yijiangh/Dropbox (MIT)/Projects/conmech/conmech/src/bindings/pyconmech/test/assembly_instances/extrusion/";
  std::string result_file_dir =
    "/Users/yijiangh/Dropbox (MIT)/Projects/conmech/conmech/src/stiffness_checker/test";

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

int main(int argc, char** argv)
{
  std::stringstream ss;
  bool use_exist_id = false;
  if (argc >= 3) {
    ss = std::stringstream(argv[2]);
    assert(ss >> std::boolalpha >> use_exist_id);
  }

  bool verbose = false;
  if (argc >= 4) {
    ss = std::stringstream(argv[3]);
    assert(ss >> std::boolalpha >> verbose);
  }

  int test_rep = 1;
  if (argc >= 5) {
    sscanf(argv[4], "%d", &test_rep);
  }

  bool reinit = false;
  if (argc >= 6) {
    ss = std::stringstream(argv[5]);
    assert(ss >> std::boolalpha >> reinit);
  }

  repTestStiffness(argv[1], use_exist_id, verbose, test_rep, reinit);

  return 0;
}
