#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <stdio.h>

// TODO: replace with: https://github.com/nlohmann/json
// https://nlohmann.github.io/json/
#include <rapidjson/document.h>
#include <rapidjson/filereadstream.h>
#include <rapidjson/filewritestream.h>
#include <rapidjson/prettywriter.h>
#include <rapidjson/cursorstreamwrapper.h>

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
  // TODO: this is super outdated! Change this!

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
  // TODO: this is super outdated! Change this!
  if ("kN/m3" == unit)
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
bool parseMaterialPropertiesJson(const std::string &file_path, std::vector<StiffnessParm> &frame_parms)
{
  using namespace rapidjson;

  FILE *fp = fopen(file_path.c_str(), "r");

  // assert(fp);
  // try {
  if(!fp) {
    throw std::runtime_error("Frame json file not found (material info should be included there)!\n");
  }
  // } catch (const std::runtime_error &e) {
  //   fprintf(stderr, "%s\n", e.what());
  //   fclose(fp);
  //   // return false;
  //   throw;
  // }

  // try {
	  char readBuffer[65536];
	  rapidjson::FileReadStream is(fp, readBuffer, sizeof(readBuffer));
	  rapidjson::MemoryPoolAllocator<> alloc;
		rapidjson::Document document(&alloc);

    CursorStreamWrapper<FileReadStream> csw(is);
		while(1) {
        document.ParseStream<rapidjson::kParseStopWhenDoneFlag,rapidjson::UTF8<>>(csw);
				if(document.HasParseError()) {
					break;
				}
    }
		fclose(fp);
    
    frame_parms.clear();

    if (!document.HasMember("uniform_cross_section") || !document.HasMember("uniform_material_properties")) { 
      throw std::runtime_error("Frame json file uniform cross sec / mat properties flag not specified!");
    }

    // TODO: make parsing uniform properties and element-wise properties the same function
    if (document["uniform_cross_section"].GetBool() && document["uniform_material_properties"].GetBool()) {
      StiffnessParm frame_parm;
      double unit_conversion;

      // TODO: check unit
      // kN/cm^2 -> kN/m^2
      if(!document["material_properties"].HasMember("youngs_modulus")) {
        throw std::runtime_error("Young's modulus property value not specified!");
      }
      unit_conversion = 1e4;
      frame_parm.youngs_modulus_ = unit_conversion * document["material_properties"]["youngs_modulus"].GetDouble();

      unit_conversion = 1e4;
      frame_parm.shear_modulus_ = unit_conversion * document["material_properties"]["shear_modulus"].GetDouble();

      // TODO: add tensile strength
      // kN/cm^2 -> kN/m^2
      // unit_conversion = 1e4;
      // frame_parm.tensile_yeild_stress_ = unit_conversion * document["material_properties"]["tensile_yeild_stress"].GetDouble();

      if (document["material_properties"].HasMember("shear_modulus")){
        unit_conversion = 1e4;
        frame_parm.shear_modulus_ = unit_conversion * document["material_properties"]["shear_modulus"].GetDouble();
        frame_parm.poisson_ratio_ = frame_parm.getPoissonRatio(frame_parm.youngs_modulus_, frame_parm.shear_modulus_);
        if (document["material_properties"].HasMember("poisson_ratio")) {
          double tmp_poisson_ratio = document["material_properties"]["poisson_ratio"].GetDouble();
          if (std::abs(tmp_poisson_ratio - frame_parm.poisson_ratio_) > 1e-5) {
              std::cout << "Warning: shear modulus and poission ratio not compatible! Diff: " <<  std::abs(tmp_poisson_ratio - frame_parm.poisson_ratio_) << ", check formula between E, G, mu. Given E & G are used." << std::endl;
            }
          }
        } else {
        if (!document["material_properties"].HasMember("poisson_ratio")){
          throw std::runtime_error("Both shear modulus and poisson_ratio not specified!");
        }
        frame_parm.poisson_ratio_ = document["material_properties"]["poisson_ratio"].GetDouble();
        frame_parm.shear_modulus_ = frame_parm.getShearModulus(frame_parm.youngs_modulus_, frame_parm.poisson_ratio_);
      }

      // kN/m^3
      if (!document["material_properties"].HasMember("density_unit") || !document["material_properties"].HasMember("density")){
        throw std::runtime_error("Density not specified!");
      }
      unit_conversion = convertDensityScale(document["material_properties"]["density_unit"].GetString());
      frame_parm.density_ = unit_conversion * document["material_properties"]["density"].GetDouble();

      // cm^2 -> m^2
      unit_conversion = 1e-4;
      if(!document["material_properties"].HasMember("cross_sec_area")) {
        throw std::runtime_error("Cross section area property value not specified!");
      }
      frame_parm.cross_sec_area_ = unit_conversion * document["material_properties"]["cross_sec_area"].GetDouble();

      // cm^4 -> m^4
      if (!document["material_properties"].HasMember("Jx") || 
          !document["material_properties"].HasMember("Iy") ||
          !document["material_properties"].HasMember("Iz")){
        throw std::runtime_error("Jx, Iy or Iz not specified!");
      }
      unit_conversion = 1e-8;
      frame_parm.Jx_ = unit_conversion * document["material_properties"]["Jx"].GetDouble();
      frame_parm.Iy_ = unit_conversion * document["material_properties"]["Iy"].GetDouble();
      frame_parm.Iz_ = unit_conversion * document["material_properties"]["Iz"].GetDouble();
      
      int element_num = document["element_list"].Size();
      for (int i=0; i<element_num; i++) {
        frame_parms.push_back(frame_parm);
      }
    }
    else {
      for(int i=0; i<document["element_list"].Size(); i++) {
        StiffnessParm frame_parm;
        double unit_conversion;

        // TODO: check unit
        // kN/cm^2 -> kN/m^2

        if(!document["element_list"][i]["material_properties"].HasMember("youngs_modulus")) {
          char msg[100];
          sprintf(msg, "Element %i: Young's modulus property value not specified!", i);
          throw std::runtime_error(msg);
        }
        unit_conversion = 1e4;
        frame_parm.youngs_modulus_ = unit_conversion * document["element_list"][i]["material_properties"]["youngs_modulus"].GetDouble();

        // kN/cm^2 -> kN/m^2
        // unit_conversion = 1e4;
        // frame_parm.tensile_yeild_stress_ = unit_conversion * document["material_properties"]["tensile_yeild_stress"].GetDouble();

        if (document["element_list"][i]["material_properties"].HasMember("shear_modulus")){
          unit_conversion = 1e4;
          frame_parm.shear_modulus_ = unit_conversion * document["element_list"][i]["material_properties"]["shear_modulus"].GetDouble();
          frame_parm.poisson_ratio_ = frame_parm.getPoissonRatio(frame_parm.youngs_modulus_, frame_parm.shear_modulus_);
          if (document["element_list"][i]["material_properties"].HasMember("poisson_ratio")) {
            double tmp_poisson_ratio = document["element_list"][i]["material_properties"]["poisson_ratio"].GetDouble();
            if (std::abs(tmp_poisson_ratio - frame_parm.poisson_ratio_) > 1e-5) {
              std::cout << "Warning: Element " << i << " shear modulus and poission ratio not compatible! Diff: " <<  std::abs(tmp_poisson_ratio - frame_parm.poisson_ratio_) << ", check formula between E, G, mu. Given E & G are used." << std::endl;
            }
          }
        } else {
          if (!document["element_list"][i]["material_properties"].HasMember("poisson_ratio")){
            char msg[100];
            sprintf(msg, "Element %i: Both shear modulus and poisson_ratio not specified!", i);
            throw std::runtime_error(msg);
          }
          frame_parm.poisson_ratio_ = document["element_list"][i]["material_properties"]["poisson_ratio"].GetDouble();
          frame_parm.shear_modulus_ = frame_parm.getShearModulus(frame_parm.youngs_modulus_, frame_parm.poisson_ratio_);
        }
        
        // kN/m^3
        if (!document["element_list"][i]["material_properties"].HasMember("density_unit") || !document["element_list"][i]["material_properties"].HasMember("density")){
          char msg[100];
          sprintf(msg, "Element %i: Density not specified!", i);
          throw std::runtime_error(msg);
        }
        unit_conversion = convertDensityScale(document["element_list"][i]["material_properties"]["density_unit"].GetString());
        frame_parm.density_ = unit_conversion * document["element_list"][i]["material_properties"]["density"].GetDouble();

        // cm^2 -> m^2
        unit_conversion = 1e-4;
        if(!document["element_list"][i]["material_properties"].HasMember("cross_sec_area")) {
          char msg[100];
          sprintf(msg, "Element %i: Cross section area property value not specified!", i);
          throw std::runtime_error(msg);
        }
        frame_parm.cross_sec_area_ = unit_conversion * document["element_list"][i]["material_properties"]["cross_sec_area"].GetDouble();
      
        // cm^4 -> m^4
        unit_conversion = 1e-8;
        if (!document["element_list"][i]["material_properties"].HasMember("Jx") || 
            !document["element_list"][i]["material_properties"].HasMember("Iy") ||
            !document["element_list"][i]["material_properties"].HasMember("Iz")){
          char msg[100];
          sprintf(msg, "Element %i: Jx, Iy or Iz not specified!", i);
          throw std::runtime_error(msg);
        }
        frame_parm.Jx_ = unit_conversion * document["element_list"][i]["material_properties"]["Jx"].GetDouble();
        frame_parm.Iy_ = unit_conversion * document["element_list"][i]["material_properties"]["Iy"].GetDouble();
        frame_parm.Iz_ = unit_conversion * document["element_list"][i]["material_properties"]["Iz"].GetDouble();

        frame_parms.push_back(frame_parm);
      }
    }
  // } catch (const std::runtime_error &e) {
	// 	fclose(fp);
  //   fprintf(stderr, "%s\n", e.what());
  //   throw;
  // }
  return true;
}

bool parseLoadCaseJson(const std::string &file_path, Eigen::MatrixXd& Load, bool& include_sw)
{
  using namespace rapidjson;

  FILE *fp = fopen(file_path.c_str(), "r");

  // assert(fp);
  try {
    if(!fp) {
      throw std::runtime_error("Load case json file not found!");
    }
  } catch (const std::runtime_error &e) {
    fprintf(stderr, "%s\n", e.what());
    fclose(fp);
    return false;
  }

  try {
	  char readBuffer[65536];
	  rapidjson::FileReadStream is(fp, readBuffer, sizeof(readBuffer));
	  rapidjson::MemoryPoolAllocator<> alloc;
		rapidjson::Document document(&alloc);

    CursorStreamWrapper<FileReadStream> csw(is);
		while(1) {
        document.ParseStream<rapidjson::kParseStopWhenDoneFlag,rapidjson::UTF8<>>(csw);
        // document.ParseStream(csw);

				if(document.HasParseError()) {
	    		// fprintf(stderr, "\nError(line %u, col %u, char %u): %s\n",
          //   	(unsigned)csw.GetLine(), // get line number
	        //     (unsigned)csw.GetColumn(), // get column number
	        //     (unsigned)document.GetErrorOffset(), // get number of chars (last parsing error)
	        //     GetParseError_En(document.GetParseError()));
					break;
				}
    }
		fclose(fp);
    // std::cout <<  "Size:" << alloc.Size() << "\n";
    // std::cout <<  "Capacity:" << alloc.Capacity() << "\n";

	  // assert(document.HasMember("dimension"));
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
	  // assert((include_sw || document.HasMember("point_load_list"))
	  // && "the load case must specify a load case, self-weight or point loads.");
	  // assert(document["point_load_list"].Size() > 0);
	  int load_v_num  = document["point_load_list"].Size();
	  Load = Eigen::MatrixXd::Zero(load_v_num, node_full_dof + 1);

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
  } catch (const std::runtime_error &e) {
		fclose(fp);
    fprintf(stderr, "%s\n", e.what());
    return false;
  }

  return true;
}

bool write_output_json(const Eigen::MatrixXd& V, const Eigen::MatrixXi& E,
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

  document.AddMember("length_unit", "meter", allocator);
  document.AddMember("rot_angle_unit", "rad", allocator);
  document.AddMember("force_unit", "kN", allocator);
  document.AddMember("moment_unit", "kN-m", allocator);

  Value nodal_disp_container(rapidjson::kArrayType);
  for(int i=0; i<n_node; i++)
  {
    rapidjson::Value node_container(rapidjson::kObjectType);
    node_container.AddMember("node_id", int(node_displ(i,0)), allocator);

    auto pt = V.block<1, 3>(int(node_displ(i,0)), 1);
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
    element_container.AddMember("node_u_id", E(i, 0), allocator);
    element_container.AddMember("node_v_id", E(i, 1), allocator);

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

    auto pt = V.block<1, 3>(int(fixities_reaction(i,0)), 1);
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
