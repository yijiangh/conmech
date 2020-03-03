#include "stiffness_checker/SharedConst.h"
#include "stiffness_checker/Material.h"
#include "stiffness_checker/StiffnessIO.h"
#include <fstream>
#include <nlohmann/json.hpp>

namespace
{
double convertLengthScale(const nlohmann::json& unit)
{
  if ("millimeter" == unit || "mm" == unit)
  {
    return 0.001;
  }
  if ("centimeter" == unit || "cm" == unit)
  {
    return 0.01;
  }
  if ("meter" == unit || "m" == unit)
  {
    // ! default unit inside stiffness_checker
    return 1;
  }
  if ("inch" == unit || "in" == unit)
  {
    return 0.0254;
  }
  if ("foot" == unit || "ft" == unit)
  {
    return 0.3048;
  }

  // default millimeter
  std::cout << "WARNING: unrecognized length unit in the input json file. Using millimeter by default." << std::endl;
  return 1;
}
}

namespace conmech
{
namespace stiffness_checker
{
bool parseFrameJson(const std::string& file_path, 
                    Eigen::MatrixXd& V, Eigen::MatrixXi& E, 
                    Eigen::MatrixXi& Fixities, std::vector<conmech::material::Material>& materials)
{
  using json = nlohmann::json;

  std::ifstream is(file_path);
  if (!is.is_open()) {
    throw std::runtime_error("Couldn't open frame file: " + file_path);
  }
  nlohmann::json json_data;
  is >> json_data;
  
  return parseFrameJson(json_data, V, E, Fixities, materials);
}

bool parseFrameJson(const nlohmann::json& json_data, 
                    Eigen::MatrixXd& V, Eigen::MatrixXi& E, 
                    Eigen::MatrixXi& Fixities, std::vector<conmech::material::Material>& materials)
{
  using json = nlohmann::json;
  using namespace conmech::material;

	int dim = int(json_data["dimension"]);
	int node_full_dof = 0;
	if(3 == dim)
	{
	  node_full_dof = 6;
	}
	else
	{
	  node_full_dof = 3;
	}

  const int n_Node = json_data["node_list"].size();
  const int n_Element = json_data["element_list"].size();
	V = Eigen::MatrixXd::Zero(n_Node, 3);
	E = Eigen::MatrixXi::Zero(n_Element, 2);

  // parsing vertices position
  double length_scale = convertLengthScale(json_data["unit"]);
  std::vector<int> fixed_node_ids;
  json j_nodes = json_data["node_list"];
  for (json::iterator it = j_nodes.begin(); it != j_nodes.end(); ++it) 
  {
    json j_node = it.value();
    // we ignore the "node_id" entry here
    int v_id = int(it - j_nodes.begin());
    V(v_id, 0) = double(j_node["point"]["X"])*length_scale;
    V(v_id, 1) = double(j_node["point"]["Y"])*length_scale;
    V(v_id, 2) = double(j_node["point"]["Z"])*length_scale;
    if (int(j_node["is_grounded"]) == 1)
    {
       fixed_node_ids.push_back(v_id);
    }
  }

  // parse fixities
  const int n_VFix = fixed_node_ids.size();
	Fixities = Eigen::MatrixXi::Zero(n_VFix, 1+node_full_dof);
  for (int i=0; i<n_VFix; i++)
  {
    int fv_id = fixed_node_ids[i];
    Fixities(i,0) = fv_id;
    if (json_data["node_list"][fv_id].contains("fixities"))
    {
      json j_fix = json_data["node_list"][fv_id]["fixities"];
      if (int(j_fix.size()) != node_full_dof) 
      {
        throw std::runtime_error("Given fixities does not have enough dof specified: " + 
          std::to_string(j_fix.size()) + " | " + std::to_string(node_full_dof));
      }
      for(int di=0; di<node_full_dof; di++)
      {
        // shift to the right by one because fv_id
        Fixities(i, di+1) = int(j_fix[di]);
      }
    }
    else
    {
      // assume all fixed
      Fixities.block(i,1,1,node_full_dof) = Eigen::VectorXi::Constant(node_full_dof, 1).transpose();
    }
  }

  // parse elements
  json j_elements = json_data["element_list"];
  for (json::iterator it = j_elements.begin(); it != j_elements.end(); ++it) 
  {
    json j_element = it.value();
    // we ignore the "element_id" entry here
    int e_id = int(it - j_elements.begin());
    E(e_id, 0) = int(j_element["end_node_ids"][0]);
    E(e_id, 1) = int(j_element["end_node_ids"][1]);
  }

  // parse materials
  if (!json_data.contains("uniform_cross_section") || !json_data.contains("uniform_material_properties"))
  {
    throw std::runtime_error("Missing attributes in json: uniform_cross_section or uniform_material_properties!");
  }
  materials.clear();
  bool uniform = bool(json_data["uniform_cross_section"]) && bool(json_data["uniform_material_properties"]);
  for (int i=0; i<n_Element; i++)
  {
    Material m;
    if (uniform)
    {
      parseMaterialPropertiesJson(json_data["material_properties"], m);
    }
    else 
    {
      parseMaterialPropertiesJson(json_data["element_list"][i]["material_properties"], m);
    }
    materials.push_back(m);
  }
  return true;
}

bool parseLoadCaseJson(const std::string &file_path, Eigen::MatrixXd& Load, bool& include_sw)
{
  using json = nlohmann::json;
  std::ifstream is(file_path);
  if (!is.is_open()) {
    throw std::runtime_error("Couldn't open frame file: " + file_path);
  }
  nlohmann::json document;
  is >> document;

	int dim = int(document["dimension"]);
	int node_full_dof = 0;
	if(3 == dim)
	{
	  node_full_dof = 6;
	}
	else
	{
	  node_full_dof = 3;
	}

	if(document.contains("include_self_weight"))
	{
	  include_sw = bool(document["include_self_weight"]);
	}
	else
	{
	  include_sw = false;
	}

	int load_v_num  = int(document["point_load_list"].size());
	Load = Eigen::MatrixXd::Zero(load_v_num, node_full_dof + 1);

	for(int i=0; i<load_v_num; i++)
	{
	  json p = document["point_load_list"][i];
	  Load(i,0) = int(p["applied_node_id"]);

	  if (3 == dim)
	  {
	    Load(i,1) = double(p["Fx"]);
	    Load(i,2) = double(p["Fy"]);
	    Load(i,3) = double(p["Fz"]);
	    Load(i,4) = double(p["Mx"]);
	    Load(i,5) = double(p["My"]);
	    Load(i,6) = double(p["Mz"]);
	  }
	  else
	  {
	    Load(i,1) = double(p["Fx"]);
	    Load(i,2) = double(p["Fy"]);
	    Load(i,3) = double(p["Mz"]);
	  }
	}

  return true;
}

bool write_output_json(const Eigen::MatrixXd& V, const Eigen::MatrixXi& E,
                       const Eigen::MatrixXd &node_displ,
                       const Eigen::MatrixXd &fixities_reaction,
                       const Eigen::MatrixXd &element_reaction,
                       const std::string& file_path)
{
  using json = nlohmann::json;

  json json_data;
  json_data["length_unit"] = conmech::LENGTH_UNIT;
  json_data["rot_angle_unit"] = conmech::ROT_ANGLE_UNIT;
  json_data["force_unit"] = conmech::FORCE_UNIT;
  json_data["moment_unit"] = conmech::MOMENT_UNIT;

  int n_NodeDp = node_displ.rows();
  int n_ElementR = element_reaction.rows();

  json j_nodal_disp;
  for(int i=0; i<n_NodeDp; i++)
  {
    json j_node;
    j_node["node_id"] = int(node_displ(i,0));

    json j_node_pose;
    auto pt = V.block<1, 3>(int(node_displ(i,0)), 1);
    j_node_pose.push_back(pt[0]);
    j_node_pose.push_back(pt[1]);
    j_node_pose.push_back(pt[2]);

    json j_node_dp;
    for(int j=1; j<node_displ.cols(); j++)
    {
      j_node_dp.push_back(double(node_displ(i,j)));
    }
    j_node["node_pose"] = j_node_pose;
    j_node["displacement"] = j_node_dp;
    j_nodal_disp.push_back(j_node);
  }
  json_data["node_displacement"] = j_nodal_disp;

  json j_element_reaction;
  for(int i=0; i<n_ElementR; i++)
  {
    json j_element;
    int e_id = int(element_reaction(i,0));
    j_element["element_id"] = e_id;
    j_element["node_u_id"] = E(e_id, 0);
    j_element["node_v_id"] = E(e_id, 1);

    json e_react;
    for(int j=1; j < element_reaction.cols(); j++)
    {
      e_react.push_back(element_reaction(i,j));
    }
    j_element["reaction"] = e_react;
    j_element_reaction.push_back(j_element);
  }
  json_data["element_reaction"] = j_element_reaction;

  int n_fix = fixities_reaction.rows();
  json j_fixity_reaction;
  for(int i=0; i<n_fix; i++)
  {
    json j_fix;
    j_fix["node_id"] = int(fixities_reaction(i,0));

    auto pt = V.block<1, 3>(int(fixities_reaction(i,0)), 1);
    json j_node_pose;
    j_node_pose.push_back(pt[0]);
    j_node_pose.push_back(pt[1]);
    j_node_pose.push_back(pt[2]);

    json j_node_reaction;
    for(int j=1; j<fixities_reaction.cols(); j++)
    {
      j_node_reaction.push_back(fixities_reaction(i,j));
    }
    j_fix["node_pose"] = j_node_pose;
    j_fix["reaction"] = j_node_reaction;

    j_fixity_reaction.push_back(j_fix);
  }
  json_data["fixity_reaction"] = j_fixity_reaction;

  std::ofstream o("file_path");
  o << json_data << std::endl;
  // write prettified JSON to another file
  // std::setw(4) << 
  return true;
}

} // ns stiffness_checker
} // ns conmech
