#include <string>
#include <iostream>

//#include <stiffness_checker/StiffnessParm.hpp>
//#include <stiffness_checker/DualGraph.hpp>
//#include <stiffness_checker/StiffnessIO.hpp>
//#include <stiffness_checker/StiffnessChecker.hpp>

#include "stiffness_checker/Frame.h"

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
} // util anon ns

int main(int argc, char** argv)
{
  testLoadingFrame();
  return 0;
}