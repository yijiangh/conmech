#include <string>
#include <iostream>

#include <Eigen/Dense>

#include "stiffness_checker/Frame.h"
#include "stiffness_checker/Util.h"
#include "stiffness_checker/Stiffness.h"
#include "stiffness_checker/StiffnessIO.h"

namespace
{
void testGNUPlot()
{
  using namespace conmech::stiffness_checker;

  std::string file_dir =
    "/Users/yijiangh/Dropbox (MIT)/Projects/conmech/conmech/src/stiffness_checker/matlab/test/problem_instances/";

  std::string file_name = "sf-test_3-frame.json";
  std::string file_path = file_dir + file_name;

  std::string save_path = "/Users/yijiangh/Dropbox (MIT)/Projects/conmech/conmech/src/stiffness_checker/test";

  Frame f;

  f.loadFromJson(file_path);
  std::cout << "frame readed: vert size: " << f.sizeOfVertList() << ", element size: " << f.sizeOfElementList()
            << std::endl;
  std::cout << "fixed vert size: " << f.sizeOfFixedVert() << ", layer size: " << f.sizeOfLayer() << std::endl;

  Eigen::MatrixXd nodal_d;

  createGnuPltStaticShape(file_name, save_path, f, nodal_d, 1, 0);
}

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

void testStiffness()
{
  using namespace conmech::stiffness_checker;

  std::string file_dir =
    "/Users/yijiangh/Dropbox (MIT)/Projects/conmech/conmech/src/stiffness_checker/matlab/test/problem_instances/";

  std::string file_name = "sf-test_3-frame.json";
  std::string file_path = file_dir + file_name;

  Stiffness sf(file_path, true);

  std::string load_name = "sf-test_3-frame_load_case.json";
  std::string load_file_path = file_dir + load_name;
  bool include_sw = false;
  Eigen::MatrixXd Load;
  parseLoadCaseJson(load_file_path, Load, include_sw);

  std::vector<int> exist_e_ids;
  exist_e_ids.push_back(0);
//  exist_e_ids.push_back(1);

//  sf.setLoad(Load, true);
  sf.setSelfWeightNodalLoad();
  bool success = sf.solve();
// bool success = sf.solve(exist_e_ids);

  std::cout << "stiffness check result: " << success << std::endl;
}

} // util anon ns

int main(int argc, char** argv)
{
//  testLocalGlobalTransf();
  testStiffness();
//  testGNUPlot();

  return 0;
}