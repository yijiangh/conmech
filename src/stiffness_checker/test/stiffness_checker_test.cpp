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

void testLoadingFrame()
{
  std::string file_path_1 =
      "/home/yijiangh/Documents/assembly-instances/assembly_models/spatial_extrusion/sf-test_4-frame/sf-test_4-frameS1_10-13-2018.json";
//  std::string file_path_2 =
//      "/home/yijiangh/Documents/assembly-instances/assembly_models/spatial_extrusion/topopt-310/topopt-310_S1.0_09-17-2018.json";

  using namespace conmech::stiffness_checker;

  Frame f;

//  f.loadFromJson(file_path_2);
//  std::cout << "frame readed: vert size: " << f.sizeOfVertList() << ", element size: " << f.sizeOfElementList()
//            << std::endl;
//  std::cout << "fixed vert size: " << f.sizeOfFixedVert() << ", layer size: " << f.sizeOfLayer() << std::endl;

  f.loadFromJson(file_path_1);
  std::cout << "frame readed: vert size: " << f.sizeOfVertList() << ", element size: " << f.sizeOfElementList()
            << std::endl;
  std::cout << "fixed vert size: " << f.sizeOfFixedVert() << ", layer size: " << f.sizeOfLayer() << std::endl;


  for(int i=0; i<f.sizeOfVertList(); i++)
  {
    auto v = f.getVert(i);
    std::cout << "vert #" << i << " has degree" << v->degree() << ", ngbh elements:"
              << "grounded: " << v->isFixed() << std::endl;
    for(auto e : v->nghdElements())
    {
      std::cout << "E #" << e->id() << ", vert id: "
                << e->endVertU()->id() << ", " <<  e->endVertV()->id() << std::endl;
    }

    std::cout << "===========" << std::endl;
  }
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
  std::string file_path_1 =
      "/home/yijiangh/Documents/assembly-instances/assembly_models/spatial_extrusion/sf-test_4-frame/sf-test_3-frame.json";
  using namespace conmech::stiffness_checker;

  Stiffness sf(file_path_1, true);

  Eigen::MatrixXd ext_P(1,7);
  ext_P.setZero();
  ext_P(0,0) = 3; // node_id
  ext_P(0,3) = -500 * 1e3; // N

  std::vector<int> exist_e_ids;
  exist_e_ids.push_back(0);
//  exist_e_ids.push_back(1);

//  sf.setNodalLoad(ext_P, false);
  sf.setSelfWeightNodalLoad();
  bool success = sf.solve();
// bool success = sf.solve(exist_e_ids);

  std::cout << "stiffness check result: " << success << std::endl;
}

} // util anon ns

int main(int argc, char** argv)
{
//  testLoadingFrame();
//  testLocalGlobalTransf();
//  testStiffness();
  testGNUPlot();

  return 0;
}