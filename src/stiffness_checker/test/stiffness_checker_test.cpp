#include <string>
#include <iostream>

#include <eigen3/Eigen/Dense>
#include "stiffness_checker/Frame.h"
#include "stiffness_checker/Util.h"
#include "stiffness_checker/Stiffness.h"

namespace
{
void testLoadingFrame()
{
  std::string file_path_1 =
      "/home/yijiangh/Documents/assembly-instances/assembly_models/spatial_extrusion/voronoi/voronoi_S1_10-03-2018.json";
  std::string file_path_2 =
      "/home/yijiangh/Documents/assembly-instances/assembly_models/spatial_extrusion/topopt-310/topopt-310_S1.0_09-17-2018.json";

  using namespace conmech::stiffness_checker;

  Frame f;

  f.loadFromJson(file_path_2);
  std::cout << "frame readed: vert size: " << f.sizeOfVertList() << ", element size: " << f.sizeOfElementList()
            << std::endl;
  std::cout << "fixed vert size: " << f.sizeOfFixedVert() << ", layer size: " << f.sizeOfLayer() << std::endl;

  f.loadFromJson(file_path_1);
  std::cout << "frame readed: vert size: " << f.sizeOfVertList() << ", element size: " << f.sizeOfElementList()
            << std::endl;
  std::cout << "fixed vert size: " << f.sizeOfFixedVert() << ", layer size: " << f.sizeOfLayer() << std::endl;


  for(int i=0; i<f.sizeOfVertList(); i++)
  {
    auto v = f.getVert(i);
    std::cout << "vert #" << i << " has degree" << v->degree() << ", ngbh elements:" << std::endl;
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
      "/home/yijiangh/Documents/assembly-instances/assembly_models/spatial_extrusion/voronoi/voronoi_S1_10-03-2018.json";
  using namespace conmech::stiffness_checker;

  Stiffness sf(file_path_1, true);

  Eigen::MatrixXd empty_ext_P;
  sf.setNodalLoad(empty_ext_P, true);
  sf.solve();
}

} // util anon ns

int main(int argc, char** argv)
{
//  testLoadingFrame();
//  testLocalGlobalTransf();
  testStiffness();
  
  return 0;
}