#include <iostream>
#include <vector>
#include <string>
//#include <fstream>
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

  char readBuffer[65536];
  FileReadStream is(fp, readBuffer, sizeof(readBuffer));

  Document document;

  if (document.ParseStream(is).HasParseError())
  {
    std::cout << "ERROR parsing the input json file!\n";
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
  std::cout << Load << std::endl;

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
  std::cout << "path file saved successfully!" << std::endl;
  return true;
}

void createGnuPltStaticShape(
  const std::string &model_name,
  const std::string &plot_dir,
  const Frame &frame,
  const Eigen::MatrixXd &nodal_displ,
  double exagg_static, bool draw_deform)
{
  int nN = frame.sizeOfVertList();
  int nE = frame.sizeOfElementList();
  if (draw_deform)
  {
    assert(nodal_displ.rows() == nN && "nodal displacement size not covering all nodes!");
  }

  FILE *fpm = NULL;

  // flag for 3D shape, = ' ' or = '#'
  char D3 = ' ';

  // linewidth for original shape
  int lw = 1;
  double scale = 1;

  time_t now;

  std::string title = model_name;
  // write gnuplot plotting script commands

  //for (j = 0; j < nN; j++)
  //{
  //	// check for three-dimensional frame
  //	if (xyz[j][0] != 0.0) X = 1;
  //	if (xyz[j][1] != 0.0) Y = 1;
  //	if (xyz[j][2] != 0.0) Z = 1;
  //}
  //
  //if ((X && Y && Z) || D3_flag)
  //{
  //	D3 = ' '; D2 = '#';
  //}
  //else
  //{
  //	D3 = '#'; D2 = ' ';
  //}

  // file name for original shape and deformed shape data
  std::string plot_path = plot_dir + "/" + title + ".plt";
  std::string orig_shape_path = plot_dir + "/" + title + "_original" + ".dat";
  std::string deformed_shape_path = plot_dir + "/" + title + "deformed" + ".dat";

  // open plotting script file for writing
  if ((fpm = fopen(plot_path.c_str(), "w")) == NULL)
  {
    std::cout << "ERROR: cannot open gnuplot script file: " << plot_path << std::endl;
    exit(1);
  }

  // write header, plot-setup cmds, node label, and element label data
  // header & node number & element number labels
  fprintf(fpm, "# CONMECH FRAME STRUCTUAL ANALYSIS RESULTS");
  fprintf(fpm, "# %s\n", title.c_str());
  fprintf(fpm, "# %s", ctime(&now));
  fprintf(fpm, "# G N U P L O T   S C R I P T   F I L E \n");

  fprintf(fpm, "set autoscale\n");
  fprintf(fpm, "unset border\n");
  fprintf(fpm, "set pointsize 1.0\n");
  fprintf(fpm, "set xtics; set ytics; set ztics; \n");
  fprintf(fpm, "unset zeroaxis\n");
  fprintf(fpm, "unset key\n");
  fprintf(fpm, "unset label\n");
  fprintf(fpm, "set size ratio -1    # 1:1 2D axis scaling \n");
  fprintf(fpm, "# set view equal xyz # 1:1 3D axis scaling \n");

  fprintf(fpm, "# NODE NUMBER LABELS\n");
  for (int j = 0; j < nN; j++)
  {
    auto pos = frame.getVert(j)->position();
    fprintf(fpm, "set label ' %d' at %12.4e, %12.4e, %12.4e\n",
            j + 1, pos[0], pos[1], pos[2]);
  }

  fprintf(fpm, "# ELEMENT NUMBER LABELS\n");
  for (int m = 0; m < nE; m++)
  {
    auto pos_u = frame.getElementEndVertU(m)->position();
    auto pos_v = frame.getElementEndVertV(m)->position();

    double mx = 0.5 * (pos_u[0] + pos_v[0]);
    double my = 0.5 * (pos_u[1] + pos_v[1]);
    double mz = 0.5 * (pos_u[2] + pos_v[2]);
    fprintf(fpm, "set label ' %d' at %12.4e, %12.4e, %12.4e\n",
            m + 1, mx, my, mz);
  }

  // 3D plot setup commands
  fprintf(fpm, "%c set parametric\n", D3);
  fprintf(fpm, "%c set view 60, 70, %5.2f \n", D3, scale);
  fprintf(fpm, "%c set view equal xyz # 1:1 3D axis scaling \n", D3);
  fprintf(fpm, "%c unset key\n", D3);
  fprintf(fpm, "%c set xlabel 'x'\n", D3);
  fprintf(fpm, "%c set ylabel 'y'\n", D3);
  fprintf(fpm, "%c set zlabel 'z'\n", D3);
  //	 fprintf(fpm,"%c unset label\n", D3 );

  fprintf(fpm, "set title \"%s\"\n", title.c_str());

  if (draw_deform)
  {
    fprintf(fpm, "# deflection exaggeration: %.1f ", exagg_static);
  }

  fprintf(fpm, "unset clip; \nset clip one; set clip two\n");
//  fprintf(fpm, "set xyplane 0 \n"); // requires Gnuplot >= 4.6

  // 2D plot command
  //fprintf(fpm, "%c plot '%s' u 2:3 t 'undeformed mesh' w lp ",
  //	D2, meshpath);
  //if (!draw_deform) fprintf(fpm, "lw %d lt 1 pt 6 \n", lw);
  //else fprintf(fpm, "lw 1 lt 5 pt 6, '%s' u 1:2 t 'load case %d of %d' w l lw %d lt 3\n", meshfl, 1, 1, lw);

  // 3D plot command, draw original
  fprintf(fpm, "%c splot '%s' using 2:3:4 title 'original shape' with linespoints ",
          D3, orig_shape_path.c_str());

  if (!draw_deform)
  {
    // original shape
    fprintf(fpm, " linewidth %d linetype 1 pointtype 6 \n", lw);
  } else
  {
    // deformed shape
    fprintf(fpm,
            " linewidth 1 linetype 5 pointtype 6, '%s' using 1:2:3 title 'load case %d of %d' with lines linewidth %d linetype 3\n",
            deformed_shape_path.c_str(), 1, 1, lw);
  }

  if (draw_deform)
  {
    // wait until a carriage return is hit
    fprintf(fpm, "pause -1\n");
  }

  fclose(fpm);
  // end of writting gnuplot script file

  // write undeformed shape data
  // open the undeformed shape data file for writing
  if ((fpm = fopen(orig_shape_path.c_str(), "w")) == NULL)
  {
    std::cout << "ERROR: cannot open gnuplot undeformed shape data file: " << orig_shape_path << std::endl;
    exit(1);
  }

  fprintf(fpm, "# CONMECH FRAME STRUCTURAL ANALYSIS RESULTS\n");
  fprintf(fpm, "# %s\n", title.c_str());
  fprintf(fpm, "# %s", ctime(&now));
  fprintf(fpm, "# U N D E F O R M E D   S H A P E   D A T A   (global coordinates)\n");
  fprintf(fpm, "# Node        X            Y            Z \n");

  for (int m = 0; m < nE; m++)
  {
    const auto node_u = frame.getElementEndVertU(m);
    const auto node_v = frame.getElementEndVertV(m);
    auto pos_u = node_u->position();
    auto pos_v = node_v->position();

    double mx = 0.5 * (pos_u[0] + pos_v[0]);
    double my = 0.5 * (pos_u[1] + pos_v[1]);
    double mz = 0.5 * (pos_u[2] + pos_v[2]);

    fprintf(fpm, "# element %5d \n", m);
    fprintf(fpm, "%5d %12.4e %12.4e %12.4e \n",
            node_u->id(), pos_u[0], pos_u[1], pos_u[2]);
    fprintf(fpm, "%5d %12.4e %12.4e %12.4e",
            node_v->id(), pos_v[0], pos_v[1], pos_v[2]);
    fprintf(fpm, "\n\n\n");
  }
  fclose(fpm);

  if (!draw_deform)
  {
    return;  // no deformed shape
  }

  // write deformed shape data
  // open the deformed mesh data file for writing
  if ((fpm = fopen(deformed_shape_path.c_str(), "w")) == NULL)
  {
    std::cout << "ERROR: cannot open gnuplot deformed shape data file: " << deformed_shape_path << std::endl;
    exit(22);
  }

  fprintf(fpm, "# CONMECH FRAME STRUCTURAL ANALYSIS RESULTS\nmac");
  fprintf(fpm, "# %s\n", title.c_str());
  fprintf(fpm, "# %s", ctime(&now));
  fprintf(fpm, "# DEFORMED FRAME SHAPE DATA");
  fprintf(fpm, "#  deflection exaggeration: %.1f\n", exagg_static);
  fprintf(fpm, "#       X-dsp        Y-dsp        Z-dsp\n");

  for (int m = 0; m < nE; m++)
  {
    // write deformed shape data for each element
    fprintf(fpm, "# element %5d \n", m);

    std::vector<Eigen::Vector3d> deformed_e_pts;
//      GnuPltCubicBentBeam(deformed_e_pts, D, m, ptr_dualgraph, ptr_wf, exagg_static);

    for (int i = 0; i < deformed_e_pts.size(); i++)
    {
      fprintf(fpm, " %12.4e %12.4e %12.4e\n",
              deformed_e_pts[i][0],
              deformed_e_pts[i][1],
              deformed_e_pts[i][2]);
    }
    fprintf(fpm, "\n\n\n");
  }

  fclose(fpm);
}

/**
 * createGnuPltCubicBentBeam  -  computes cubic deflection functions from end deflections
 * and end rotations.  Saves deflected shapes to a file.  These bent shapes
 * are exact for mode-shapes, and for frames loaded at their nodes.
*/
void createGnuPltCubicBentBeam(
  const Frame &frame,
  const Eigen::MatrixXd &nodal_displ,
  std::vector<Eigen::Vector3d> &beam,
  double exagg)
{
//	double	t0, t1, t2, t3, t4, t5, t6, t7, t8, 	/* coord transf matrix entries	*/
//		u1, u2, u3, u4, u5, u6, u7, u8, u9, u10, u11, u12,
//		s, v, w, dX, dY, dZ;
//
//	WF_edge *ei = ptr_frame->GetEdge(ptr_dualgraph->e_orig_id(dual_i));
//	WF_edge *ej = ei->ppair_;
//	int dual_u = ptr_dualgraph->v_dual_id(ej->pvert_->ID());
//	int dual_v = ptr_dualgraph->v_dual_id(ei->pvert_->ID());
//
//	double L = ei->Length();
//
//	int	i1, i2;
//	int info;
//	char	errMsg[MAXL];
//
//	MX A(4,4);
//	VX a(4);
//	VX b(4);
//
//	trsf_.CreateTransMatrix(ej->pvert_->Position(), ei->pvert_->Position(),
//		t0, t1, t2, t3, t4, t5, t6, t7, t8, 0.0);
//
//	i1 = 6 * dual_u;	i2 = 6 * dual_v;
//
//	/* compute end deflections in local coordinates */
//
//	u1 = exagg*(t0*D[i1 + 0] + t1*D[i1 + 1] + t2*D[i1 + 2]);
//	u2 = exagg*(t3*D[i1 + 0] + t4*D[i1 + 1] + t5*D[i1 + 2]);
//	u3 = exagg*(t6*D[i1 + 0] + t7*D[i1 + 1] + t8*D[i1 + 2]);
//
//	u4 = exagg*(t0*D[i1 + 3] + t1*D[i1 + 4] + t2*D[i1 + 5]);
//	u5 = exagg*(t3*D[i1 + 3] + t4*D[i1 + 4] + t5*D[i1 + 5]);
//	u6 = exagg*(t6*D[i1 + 3] + t7*D[i1 + 4] + t8*D[i1 + 5]);
//
//	u7 = exagg*(t0*D[i2 + 0] + t1*D[i2 + 1] + t2*D[i2 + 2]);
//	u8 = exagg*(t3*D[i2 + 0] + t4*D[i2 + 1] + t5*D[i2 + 2]);
//	u9 = exagg*(t6*D[i2 + 0] + t7*D[i2 + 1] + t8*D[i2 + 2]);
//
//	u10 = exagg*(t0*D[i2 + 3] + t1*D[i2 + 4] + t2*D[i2 + 5]);
//	u11 = exagg*(t3*D[i2 + 3] + t4*D[i2 + 4] + t5*D[i2 + 5]);
//	u12 = exagg*(t6*D[i2 + 3] + t7*D[i2 + 4] + t8*D[i2 + 5]);
//
//	/* curve-fitting problem for a cubic polynomial */
//
//	a[0] = u2;		b[0] = u3;
//	a[1] = u8;   	b[1] = u9;
//	a[2] = u6;		b[2] = -u5;
//	a[3] = u12;		b[3] = -u11;
//
//	u7 += L;
//	A(0,0) = 1.0;   A(0,1) = u1;   A(0,2) = u1*u1;   A(0,3) = u1*u1*u1;
//	A(1,0) = 1.0;   A(1,1) = u7;   A(1,2) = u7*u7;   A(1,3) = u7*u7*u7;
//	A(2,0)= 0.0;    A(2,1) = 1.;   A(2,2) = 2.*u1;   A(2,3) = 3.*u1*u1;
//	A(3,0) = 0.0;   A(3,1) = 1.;   A(3,2) = 2.*u7;   A(3,3) = 3.*u7*u7;
//	u7 -= L;
//
//	VX xa(4);
//	info = solver_.LUDecomp(A, xa, a);		/* solve for cubic coef's */
//
//	if (!info)
//	{
//		sprintf(errMsg, " n1 = %d  n2 = %d  L = %e  u7 = %e \n", dual_u+1, dual_v+1, L, u7);
//		printf(errMsg);
//		exit(30);
//	}
//
//	VX xb(4);
//	info = solver_.LUDecomp(A, xb, b);		/* solve for cubic coef's */
//
//	if (!info)
//	{
//		sprintf(errMsg, " n1 = %d  n2 = %d  L = %e  u7 = %e \n", dual_u + 1, dual_v + 1, L, u7);
//		printf(errMsg);
//		exit(30);
//	}
//
//	// debug ... if deformed mesh exageration is too big, some elements
//	// may not be plotted.
//	//fprintf( fpm, "# u1=%e  L+u7=%e, dx = %e \n",
//	//				u1, fabs(L+u7), fabs(L+u7-u1)/10.0);
//	for (s = u1; fabs(s) <= 1.01*fabs(L + u7); s += fabs(L + u7 - u1) / 10.0)
//	{
//
//		/* deformed shape in local coordinates */
//		v = xa[0] + xa[1] * s + xa[2] * s*s + xa[3] * s*s*s;
//		w = xb[0] + xb[1] * s + xb[2] * s*s + xb[3] * s*s*s;
//
//		/* deformed shape in global coordinates */
//		dX = t0*s + t3*v + t6*w;
//		dY = t1*s + t4*v + t7*w;
//		dZ = t2*s + t5*v + t8*w;
//
//		beam.push_back(point(ej->pvert_->Position().x() + dX,
//			ej->pvert_->Position().y() + dY, ej->pvert_->Position().z() + dZ));
//		//fprintf(fpm, " %12.4e %12.4e %12.4e\n",);
//	}
//
//	return;
}

} // ns stiffness_checker
} // ns conmech
