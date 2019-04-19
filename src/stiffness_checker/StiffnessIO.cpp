#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <stdio.h>

#include <rapidjson/document.h>
#include <rapidjson/filereadstream.h>
#include <rapidjson/filewritestream.h>
#include <rapidjson/prettywriter.h>

#include "stiffness_checker/Frame.h"
#include "stiffness_checker/StiffnessParm.h"
#include "stiffness_checker/StiffnessIO.h"

namespace
{

double convertLengthScale(const std::string &unit)
{
  if ("millimeter" == unit || "mm" == unit)
  {
    return 1;
  }
  if ("centimeter" == unit || "cm" == unit)
  {
    return 10;
  }
  if ("meter" == unit || "m" == unit)
  {
    return 1000;
  }
  if ("inch" == unit || "in" == unit)
  {
    return 25.4;
  }
  if ("foot" == unit || "ft" == unit)
  {
    return 304.8;
  }

  // default millimeter
  std::cout << "WARNING: unrecognized length unit in the input json file. Using millimeter by default." << std::endl;
  return 1;
}

double convertModulusScale(const std::string &unit)
{
  // 1 pa = 1 N/m^2 = 1e-4 N/cm^2 = 1e-6 N/mm^2
  // 1 Mpa = 1 N/mm^2

  if ("MPa" == unit)
  {
    return 1;
  }
  if ("GPa" == unit)
  {
    return 1e3;
  }

  // default MPa
  std::cout << "WARNING: unrecognized modulus unit in the input json file. "
               "Using MPa by default." << std::endl;
  return 1;
}

double convertDensityScale(const std::string &unit)
{
  if ("kg/m3" == unit)
  {
    return 1;
  }

  // default MPa
//  std::cout << "WARNING: unrecognized material density unit in the input json file. "
//      "Using kg/m3 by default." << std::endl;
  return 1;
}

} // anon util ns

namespace conmech
{
namespace stiffness_checker
{
bool parseMaterialPropertiesJson(const std::string &file_path, StiffnessParm &frame_parm)
{
  using namespace rapidjson;

  FILE *fp = fopen(file_path.c_str(), "r");

  assert(fp);

  char readBuffer[65536];
  FileReadStream is(fp, readBuffer, sizeof(readBuffer));

  Document document;

  if (document.ParseStream(is).HasParseError())
  {
    std::cout << "ERROR parsing the input json file!\n";
    return false;
  }

  fclose(fp);

  assert(document.HasMember("material_properties"));

  double unit_conversion;

  // TODO: check unit
  // kN/cm^2 -> kN/m^2
  assert(document["material_properties"].HasMember("youngs_modulus_unit"));
  assert(document["material_properties"].HasMember("youngs_modulus"));
//  unit_conversion = convertModulusScale(document["material_properties"]["youngs_modulus_unit"].GetString());
  unit_conversion = 1e4;
  frame_parm.youngs_modulus_ = unit_conversion * document["material_properties"]["youngs_modulus"].GetDouble();

  assert(document["material_properties"].HasMember("shear_modulus_unit"));
  assert(document["material_properties"].HasMember("shear_modulus"));
//  unit_conversion = convertModulusScale(document["material_properties"]["shear_modulus_unit"].GetString());
  unit_conversion = 1e4;
  frame_parm.shear_modulus_ = unit_conversion * document["material_properties"]["shear_modulus"].GetDouble();

  assert(document["material_properties"].HasMember("poisson_ratio"));
  frame_parm.poisson_ratio_ = document["material_properties"]["poisson_ratio"].GetDouble();

  // kN/m^3
  assert(document["material_properties"].HasMember("density_unit"));
  assert(document["material_properties"].HasMember("density"));
  unit_conversion = convertDensityScale(document["material_properties"]["density_unit"].GetString());
  frame_parm.density_ = unit_conversion * document["material_properties"]["density"].GetDouble();

  // cm -> m
  assert(document["material_properties"].HasMember("radius_unit"));
  assert(document["material_properties"].HasMember("radius"));
//  unit_conversion = convertDensityScale(document["material_properties"]["radius_unit"].GetString());
  unit_conversion = 1e-2;
  frame_parm.radius_ = unit_conversion * document["material_properties"]["radius"].GetDouble();

  return true;
}

bool parseLoadCaseJson(const std::string &file_path, Eigen::MatrixXd& Load, bool& include_sw)
{
  using namespace rapidjson;

  FILE *fp = fopen(file_path.c_str(), "r");

  assert(fp);

  char readbuffer[65536];
  rapidjson::FileReadStream is(fp, readbuffer, sizeof(readbuffer));

  rapidjson::Document document;

  if (document.ParseStream(is).HasParseError())
  {
    std::cout << "error parsing the input json file!\n";
    return false;
  }

  fclose(fp);

  assert(document.HasMember("dimension"));
  int dim = document["dimension"].GetInt();
  int node_full_dof = 0;
  if(3 == dim)
  {
    node_full_dof = 6;
  }
  else
  {
    node_full_dof = 3;
  }

  if(document.HasMember("include_self_weight"))
  {
    include_sw = document["include_self_weight"].GetBool();
  }
  else
  {
    include_sw = false;
  }

  // read point loads
  assert((include_sw || document.HasMember("point_load_list"))
  && "the load case must specify a load case, self-weight or point loads.");
  assert(document["point_load_list"].Size() > 0);
  int load_v_num  = document["point_load_list"].Size();
  Load = Eigen::MatrixXd::Zero(load_v_num, node_full_dof + 1);
  // std::cout << Load << std::endl;

  for(int i=0; i<load_v_num; i++)
  {
    const Value& p = document["point_load_list"][i];
    Load(i,0) = p["applied_node_id"].GetInt();

    if (3 == dim)
    {
      Load(i,1) = p["Fx"].GetDouble();
      Load(i,2) = p["Fy"].GetDouble();
      Load(i,3) = p["Fz"].GetDouble();
      Load(i,4) = p["Mx"].GetDouble();
      Load(i,5) = p["My"].GetDouble();
      Load(i,6) = p["Mz"].GetDouble();
    }
    else
    {
      Load(i,1) = p["Fx"].GetDouble();
      Load(i,2) = p["Fy"].GetDouble();
      Load(i,3) = p["Mz"].GetDouble();
    }
  }

  return true;
}

bool write_output_json(const Frame& frame,
    const Eigen::MatrixXd &node_displ,
    const Eigen::MatrixXd &fixities_reaction,
    const Eigen::MatrixXd &element_reaction,
    const std::string& file_path)
{
  using namespace rapidjson;
  // document is the root of a json message
  rapidjson::Document document;

  // define the document as an object rather than an array
  document.SetObject();

  // must pass an allocator when the object may need to allocate memory
  rapidjson::Document::AllocatorType& allocator = document.GetAllocator();

  int n_node = node_displ.rows();
  int n_element = element_reaction.rows();

  Value nodal_disp_container(rapidjson::kArrayType);
  for(int i=0; i<n_node; i++)
  {
    rapidjson::Value node_container(rapidjson::kObjectType);
    node_container.AddMember("node_id", int(node_displ(i,0)), allocator);

    auto pt = frame.getVertPosition(node_displ(i,0));
    rapidjson::Value node_pose(rapidjson::kArrayType);
    node_pose.PushBack(Value().SetDouble(pt[0]), allocator);
    node_pose.PushBack(Value().SetDouble(pt[1]), allocator);
    node_pose.PushBack(Value().SetDouble(pt[2]), allocator);

    rapidjson::Value node_dp(rapidjson::kArrayType);
    for(int j=1; j < node_displ.cols(); j++)
    {
      node_dp.PushBack(Value().SetDouble(node_displ(i,j)), allocator);
    }

    node_container.AddMember("node_pose", node_pose, allocator);
    node_container.AddMember("displacement", node_dp, allocator);
    nodal_disp_container.PushBack(node_container, allocator);
  }
  document.AddMember("node_displacement", nodal_disp_container, allocator);

  Value element_reaction_container(rapidjson::kArrayType);
  for(int i=0; i<n_element; i++)
  {
    rapidjson::Value element_container(rapidjson::kObjectType);
    auto e_id = element_reaction(i,0);
    element_container.AddMember("element_id", int(e_id), allocator);
    element_container.AddMember("node_u_id", frame.getElementEndVertU(e_id)->id(), allocator);
    element_container.AddMember("node_v_id", frame.getElementEndVertV(e_id)->id(), allocator);

    rapidjson::Value e_react(rapidjson::kArrayType);
    for(int j=1; j < element_reaction.cols(); j++)
    {
      e_react.PushBack(Value().SetDouble(element_reaction(i,j)), allocator);
    }
    element_container.AddMember("reaction", e_react, allocator);
    element_reaction_container.PushBack(element_container, allocator);
  }
  document.AddMember("element_reaction", element_reaction_container, allocator);

  int n_fix = fixities_reaction.rows();
  Value fixity_reaction_container(rapidjson::kArrayType);
  for(int i=0; i<n_fix; i++)
  {
    rapidjson::Value fix_container(rapidjson::kObjectType);
    fix_container.AddMember("node_id", int(fixities_reaction(i,0)), allocator);

    auto pt = frame.getVertPosition(fixities_reaction(i,0));
    rapidjson::Value node_pose(rapidjson::kArrayType);
    node_pose.PushBack(Value().SetDouble(pt[0]), allocator);
    node_pose.PushBack(Value().SetDouble(pt[1]), allocator);
    node_pose.PushBack(Value().SetDouble(pt[2]), allocator);

    rapidjson::Value node_reaction(rapidjson::kArrayType);
    for(int j=1; j < fixities_reaction.cols(); j++)
    {
      node_reaction.PushBack(Value().SetDouble(fixities_reaction(i,j)), allocator);
    }

    fix_container.AddMember("node_pose", node_pose, allocator);
    fix_container.AddMember("reaction", node_reaction, allocator);

    fixity_reaction_container.PushBack(fix_container, allocator);
  }
  document.AddMember("fixity_reaction", fixity_reaction_container, allocator);

  // output file to path
  std::string json_path = file_path;
  FILE *js_file = fopen(json_path.c_str(), "w+");
  if(NULL == js_file)
  {
    std::cout << "ERROR: invalid output file path!!!" << std::endl;
    return false;
  }

  char writeBuffer[65536];
  rapidjson::FileWriteStream os(js_file, writeBuffer, sizeof(writeBuffer));

  rapidjson::PrettyWriter<FileWriteStream> p_writer(os);
  document.Accept(p_writer);

  std::fclose(js_file);
  std::cout << "computed result file saved successfully!" << std::endl;
  return true;
}

} // ns stiffness_checker
} // ns conmech
