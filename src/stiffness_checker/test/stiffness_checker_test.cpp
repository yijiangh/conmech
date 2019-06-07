#include <string>
#include <iostream>
#include <chrono>

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

void testStiffness(const char* cfile_name)
{
  using namespace conmech::stiffness_checker;

  std::string file_dir =
    "/Users/yijiangh/Dropbox (MIT)/Projects/conmech/conmech/src/bindings/pyconmech/test/assembly_instances/extrusion/";
  std::string result_file_dir =
    "/Users/yijiangh/Dropbox (MIT)/Projects/conmech/conmech/src/stiffness_checker/test";

  // set up shape file path
  std::string file_name = cfile_name;
  // "duck.json";
  // "tre_foil_knot.json";
  // "david.json";
  // "topopt-100_S1_03-14-2019_w_layer.json";
  // "sf-test_3-frame.json";
  //"simple_frame.json";
  // "tower_3D.json";
  std::string file_path = file_dir + file_name;
  Stiffness sf(file_path, false);

  // set up output json result
  std::string result_file_name = "conmech_cpp_test_result.json";
  // sf.setOutputJsonPath(result_file_dir, result_file_name);
  // sf.setOutputJson(true);

  // set up load case
  // std::string load_name = "tower_3D_load_case.json";
  // std::string load_file_path = file_dir + load_name;
  // Eigen::MatrixXd Load;
  // parseLoadCaseJson(load_file_path, Load, include_sw);
  // std::cout << "load parsed: \n" << Load << std::endl;

  // sf.setLoad(Load);

  bool include_sw = true;
  sf.setSelfWeightNodalLoad(include_sw);

  // partial structure ids
  // std::vector<int> exist_e_ids{0,1};
  // std::vector<int> exist_e_ids{
    // 125, 126, 115, 122, 111, 108, 23, 22, 98, 75, 64, 34, 61, 65, 59, 60, 39, 36, 44, 67};

  bool success = sf.solve();
  // bool success = sf.solve(exist_e_ids);
  std::cout << "stiffness check result: " << success << std::endl;

  Eigen::MatrixXd nD, fR, eR;
  sf.getSolvedResults(nD, fR, eR, success);
  std::cout << "get solved result: criteria pass : " << success << std::endl;
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

  double complaince;
  sf.getSolvedCompliance(complaince);
  std::cout << "complaince: " << complaince << std::endl;

  double exagg = 1.0;
  int disc = 2;

  bool draw_full_shape=false;
  Eigen::MatrixXd OB = sf.getOriginalShape(disc, draw_full_shape);
  // std::cout << "original beam polygonal points:\n" << OB << std::endl;

  Eigen::MatrixXd DB = sf.getDeformedShape(exagg, disc);
  // std::cout << "deformed beam polygonal points:\n" << DB << std::endl;
}

void repTestStiffness(const char* cfile_name, int iter, bool re_init)
{
  using namespace conmech::stiffness_checker;

  std::string file_dir =
    "/Users/yijiangh/Dropbox (MIT)/Projects/conmech/conmech/src/bindings/pyconmech/test/assembly_instances/extrusion/";
  std::string result_file_dir =
    "/Users/yijiangh/Dropbox (MIT)/Projects/conmech/conmech/src/stiffness_checker/test";

  // set up shape file path
  std::string file_name = cfile_name;
  // "topopt-100_S1_03-14-2019_w_layer.json";
  // "sf-test_3-frame.json";
  //"simple_frame.json";
  // "tower_3D.json";
  std::string file_path = file_dir + file_name;
  Stiffness sf(file_path, true);

  // set up output json result
  // std::string result_file_name = "conmech_cpp_test_result.json";
  // sf.setOutputJsonPath(result_file_dir, result_file_name);
  // sf.setOutputJson(true);

  bool include_sw = true;
  sf.setSelfWeightNodalLoad(include_sw);

  // partial structure ids
  // std::vector<int> exist_e_ids{0,1};
  std::vector<int> exist_e_ids{
    125, 126, 115, 122, 111, 108, 23, 22, 98, 75, 64, 34, 61, 65, 59, 60, 39, 36, 44, 67};

  // bool success = sf.solve();
  bool success = sf.solve(exist_e_ids);
  // std::cout << "stiffness check result: " << success << std::endl;

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
  //
  // double complaince;
  // sf.getSolvedCompliance(complaince);
  // std::cout << "complaince: " << complaince << std::endl;

}

} // util anon ns

int main(int argc, char** argv)
{
//  testLocalGlobalTransf();
  std::cout << "testing on file: " << argv[1] << std::endl;
  testStiffness(argv[1]);

  return 0;
}
