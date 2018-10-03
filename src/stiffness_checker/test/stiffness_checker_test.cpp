#include <string>

#include <stiffness_checker/StiffnessParm.hpp>
#include <stiffness_checker/DualGraph.hpp>
#include <stiffness_checker/StiffnessIO.hpp>
#include <stiffness_checker/StiffnessChecker.hpp>

int main(int argc, char** argv)
{
  std::string file_path =
      "/home/yijiangh/Documents/assembly-instances/assembly_models/spatial_extrusion/voronoi/voronoi_S1_10-03-2018.json";

  using namespace conmech::stiffness_checker;

  StiffnessChecker sc(file_path, true);

  std::vector<int> ids;
  ids.push_back(1);
  sc.checkDeformation(ids);

  return 0;
}