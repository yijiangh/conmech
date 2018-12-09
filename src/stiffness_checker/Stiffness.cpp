#include <cmath>
#include <set>
#include <algorithm>

#include <Eigen/Dense>

#include "stiffness_checker/Util.h"
#include "stiffness_checker/StiffnessSolver.h"
#include "stiffness_checker/StiffnessIO.h"

#include "stiffness_checker/Stiffness.h"

namespace conmech
{
namespace stiffness_checker
{

Stiffness::Stiffness(Frame& frame, bool verbose, std::string model_type)
    : frame_(frame), verbose_(verbose), is_init_(false), transl_tol_(1.0), rot_tol_(5*(3.14/180))
{
  // frame sanity check
  int N_vert = frame_.sizeOfVertList();
  int N_element = frame_.sizeOfElementList();
  assert(N_vert>0 && N_element>0);

  stiff_solver_.timing_ = verbose;
  model_type_ = model_type;

  init();
}

Stiffness::Stiffness(const std::string &json_file_path, bool verbose, std::string model_type)
    : verbose_(verbose), is_init_(false), transl_tol_(1.0), rot_tol_(5*(3.14/180))
{
  verbose_ = verbose;
  stiff_solver_.timing_ = verbose;

  // parse frame
  frame_.loadFromJson(json_file_path);

  // parse material properties
  parseMaterialPropertiesJson(json_file_path, material_parm_);

  model_type_ = model_type;

  init();
}

bool Stiffness::init()
{
  // TODO: generalize to 2D
  dim_ = 3;

  // set up dimension, node_dof
  if (3 == dim_)
  {
    // 3D case, x, y, z, xx, yy, zz
    full_node_dof_ = 6;
    if (model_type_ == "frame")
    {
      node_dof_ = 6;
      xyz_dof_id_ = Eigen::VectorXi::LinSpaced(node_dof_*2,0,node_dof_*2-1);
      e_react_dof_id_ = Eigen::VectorXi::LinSpaced(node_dof_*2,0,node_dof_*2-1);
    }
    if (model_type_ == "truss")
    {
      node_dof_ = 3;
      xyz_dof_id_ = Eigen::VectorXi(6);
      xyz_dof_id_ << 0, 1, 2, 6, 7, 8;

      e_react_dof_id_ = Eigen::VectorXi(2);
      e_react_dof_id_ << 0, 6;
    }
  }
  else
  {
    // 2D case, x, y, theta
    full_node_dof_ = 3;
    if (model_type_ == "frame")
    {
      node_dof_ = 3;
      xyz_dof_id_ = Eigen::VectorXi::LinSpaced(node_dof_*2,0,node_dof_*2-1);
      e_react_dof_id_ = Eigen::VectorXi::LinSpaced(node_dof_*2,0,node_dof_*2-1);
    }
    if (model_type_ == "truss")
    {
      node_dof_ = 2;
      xyz_dof_id_ = Eigen::VectorXi(4);
      xyz_dof_id_ << 0, 1, 3, 4;

      e_react_dof_id_ = Eigen::VectorXi(2);
      e_react_dof_id_ << 0, 3;
    }
  }

  createElementStiffnessMatrixList();

  // create id_map_
  int dof = frame_.sizeOfVertList() * node_dof_;
  int N_element = frame_.sizeOfElementList();

  id_map_.resize(N_element, node_dof_*2);
  auto lin_sp_id = Eigen::VectorXi::LinSpaced(node_dof_,0,node_dof_-1); // 0,1,..,5

  for(int i=0; i<N_element; i++)
  {
    const auto e = frame_.getElement(i);
    const auto end_u = e->endVertU();
    const auto end_v = e->endVertV();

    id_map_.block<1,6>(i,0) = Eigen::VectorXi::Constant(node_dof_, node_dof_*end_u->id()) + lin_sp_id;
    id_map_.block<1,6>(i,node_dof_) = Eigen::VectorXi::Constant(node_dof_, node_dof_*end_v->id()) + lin_sp_id;
  }
//  std::cout << id_map_ << std::endl;

  if(verbose_)
  {
    std::cout << "Initialization done" << std::endl;
  }

  is_init_ = true;
}

void Stiffness::setNodalDisplacementTolerance(double transl_tol, double rot_tol)
{
  assert(transl_tol > 0 && rot_tol > 0 && "invalid tolerance: tolerance must be bigger than 0!");
  transl_tol_ = transl_tol;
  rot_tol_ = rot_tol;
}

void Stiffness::createSelfWeightNodalLoad(Eigen::VectorXd& self_weight_load_P)
{
  assert(is_init_);

  // refer [MSA McGuire et al.] P111
  // Loads between nodal points
  int N_vert = frame_.sizeOfVertList();
  int N_element = frame_.sizeOfElementList();
  self_weight_load_P.resize(N_vert*6);

  // assuming solid circular cross sec for now
  double Ax = M_PI * std::pow(material_parm_.radius_, 2);

  // TODO: change here
  // gravitional acc unit: m/s^2
  // in global frame
  double gz = 1;

  for (int i = 0; i < N_element; i++)
  {
    const auto e = frame_.getElement(i);
    const auto end_u = e->endVertU();
    const auto end_v = e->endVertV();

    double L = e->getLength();
    // uniform force density along the element
    // due to gravity
    double q = material_parm_.density_ * Ax * L * gz;
    // gravity in global frame
    Eigen::Vector3d P_G(0,0,q*L/2.0);

    Eigen::Matrix3d R_LG;
    getGlobal2LocalRotationMatrix(
        e->endVertU()->position(),
        e->endVertV()->position(), R_LG);
    // gravity in local frame
    Eigen::Vector3d P_L = R_LG * P_G;

    // force density in local frame
    double q_yL = P_L[1] / L;
    double q_zL = P_L[2] / L;

    // according to fixed-end force diagram
    // uniform load on a beam with both ends fixed ([MSA] p111.)
    // fixed-end fictitious reaction
    Eigen::VectorXd P_fL(6);
    P_fL.head(3) = P_L; // force
    P_fL[3] = 0.0; // no moment from axial component
    P_fL[4] = q_yL * std::pow(L,2) / 12;
    P_fL[5] = q_zL * std::pow(L,2) / 12;

    // transform back to global frame
    Eigen::VectorXd P_fG(6);
    auto R_GL = R_LG.transpose();
    P_fG.head(3) = R_GL * P_fL.head(3);
    P_fG.tail(3) = R_GL * P_fL.tail(3);

    self_weight_load_P.segment(end_u->id()*6, 6) += P_fG;

    // moment on the other end needs to be negated for sign consistency
    P_fG.tail(3) *= -1;

    self_weight_load_P.segment(end_v->id()*6, 6) += P_fG;
  }
}

void Stiffness::createExternalNodalLoad(
    const Eigen::MatrixXd& nodal_forces, Eigen::VectorXd& ext_load)
{
  assert(is_init_);
  assert(nodal_forces.cols() == node_dof_ + 1);

  int full_dof = node_dof_*frame_.sizeOfVertList();
  ext_load.resize(full_dof);

  for(int i=0; i < nodal_forces.rows(); i++)
  {
    int v_id = nodal_forces(i,0);
    assert(0<=v_id && v_id<frame_.sizeOfVertList());

    ext_load.segment(v_id*6,6) = nodal_forces.block<1,6>(i,1);
  }

  std::cout << "nodal_load" << ext_load << std::endl;
}

bool Stiffness::setLoad(const Eigen::MatrixXd& nodal_forces,
                             const bool& include_self_weight)
{
  assert(is_init_);

  int dof = node_dof_ * frame_.sizeOfVertList();
  bool is_empty_ext = nodal_forces.isZero(0) || nodal_forces.rows() == 0;

  assert((!include_self_weight && !is_empty_ext) && "No load is assigned.");

  Eigen::VectorXd sw_load(dof);
  sw_load.setZero();

  Eigen::VectorXd ext_load(dof);
  ext_load.setZero();

  nodal_load_P_.resize(dof);

  if(include_self_weight)
  {
    if(verbose_)
    {
      std::cout << "self-weight turned on." << std::endl;
    }

    createSelfWeightNodalLoad(sw_load);
  }

  if(!is_empty_ext)
  {
    createExternalNodalLoad(nodal_forces, ext_load);
  }
  else
  {
    if(verbose_)
    {
      std::cout << "no external load is assigned." << std::endl;
    }
  }

  nodal_load_P_ = sw_load + ext_load;
}

bool Stiffness::setSelfWeightNodalLoad()
{
  int dof = 6*frame_.sizeOfVertList();

  Eigen::VectorXd sw_load(dof);
  sw_load.setZero();

  nodal_load_P_.resize(dof);

  createSelfWeightNodalLoad(sw_load);
  nodal_load_P_ = sw_load;
}

void Stiffness::createElementStiffnessMatrixList()
{
  // check the complete element stiffness matrix
  // in local frame: [MSA McGuire et al.] P73

  // for all cross section's geometrical properties
  // cross section: m^2
  // Iz, Iy: m^4
  // J: m^4

  // element length L: m
  // E, G: kN/m^2
  // Force: kN
  // Moment: kN m

  // assuming all the elements have the same cross section
  // and all of them are solid circular shape
  double pi = atan(1)*4;

  // cross section area: assuming solid circular, unit: m^2
  double A = pi * std::pow(material_parm_.radius_,2);
//  double Asy = Ax * (6 + 12 * material_parm_.poisson_ratio_ + 6 * std::pow(material_parm_.poisson_ratio_,2))
//      / (7 + 12 * material_parm_.poisson_ratio_ + 4 * std::pow(material_parm_.poisson_ratio_,2));
//  double Asz = Asy;

  // torsion constant (around local x axis)
  // for solid circle: (1/2)pi*r^4, unit: mm^4
  // see: https://en.wikipedia.org/wiki/Torsion_constant#Circle
  double Jx = 0.5 * pi * std::pow(material_parm_.radius_, 4);
//  double Jx = (2.0/3.0) * pi * std::pow(material_parm_.radius_, 4);

  // area moment of inertia (bending about local y,z-axis)
  // assuming solid circular area of radius r, unit: m^4
  // see: https://en.wikipedia.org/wiki/List_of_area_moments_of_inertia
  // note this is slender rod of length L and Mass M, spinning around end
  double Iy = pi * std::pow(material_parm_.radius_,4) / 4;
  double Iz = Iy;

  // E,G: MPa; mu (poisson ratio): unitless
  double E = material_parm_.youngs_modulus_;
  double G = material_parm_.shear_modulus_;
  double mu = material_parm_.poisson_ratio_;

  assert(std::abs((G - E/(2*(1+mu)))/G) < 1e-3
      && "input poisson ratio not compatible with shear and Young's modulus!");

  int N_element = frame_.sizeOfElementList();
  assert(N_element > 0);

  using namespace std;

  element_K_list_.reserve(N_element);
  for (int i = 0; i < N_element; i++)
  {
    const auto e = frame_.getElement(i);
    const auto end_u = e->endVertU();
    const auto end_v = e->endVertV();

    // element length, unit: m
    double L = (end_v->position() - end_u->position()).norm();

    Eigen::Matrix3d R_LG;
    getGlobal2LocalRotationMatrix(end_u->position(), end_v->position(), R_LG);

    // element stiffness matrix in local frame
    Eigen::MatrixXd K_eL;
    createLocalStiffnessMatrix(L, A, dim_, Jx, Iy, Iz, E, G, mu, K_eL);

    // transform to global frame
    Eigen::MatrixXd R_LG_diag = Eigen::MatrixXd::Zero(full_node_dof_*2, full_node_dof_*2);

    // 4 block or 2 block
    for(int j=0;j<(full_node_dof_/3)*2;j++)
    {
      R_LG_diag.block<3, 3>(j*3, j*3) = R_LG;
    }

    // DO NOT combine transpose
    //https://libigl.github.io/matlab-to-eigen.html
    auto R_LG_diagT = R_LG_diag.transpose();
    auto K_eG = R_LG_diagT * K_eL * R_LG_diag;

    element_K_list_.push_back(K_eG);
  }
}

void Stiffness::createCompleteGlobalStiffnessMatrix(const std::vector<int>& exist_e_ids)
{
  assert(element_K_list_.size() > 0);
  assert(exist_e_ids.size()>0);
  assert(id_map_.rows() >= exist_e_ids.size());

  int total_dof = frame_.sizeOfVertList() * node_dof_;
  int N_element = frame_.sizeOfElementList();

  K_assembled_full_.resize(total_dof, total_dof);
  K_assembled_full_.setZero();
  for(const int e_id : exist_e_ids)
  {
    assert(e_id>=0 && e_id<frame_.sizeOfElementList());

    const auto K_e = element_K_list_[e_id];
    assert(K_e.rows() == 2*node_dof_ && K_e.cols() == 2*node_dof_);

    for(int i=0; i < 2*node_dof_; i++)
    {
      int row_id = id_map_(e_id,i);

      for(int j=0; j < 2*node_dof_; j++)
      {
        int col_id = id_map_(e_id,j);
        K_assembled_full_(row_id, col_id) += K_e(i,j);
      }
    }
  }
}

bool Stiffness::solve(
    const std::vector<int>& exist_element_ids,
    Eigen::MatrixXd& node_displ,
    Eigen::MatrixXd& fixities_reaction,
    Eigen::MatrixXd& element_reaction,
    const bool& cond_num)
{
  //todo
  using namespace std;

  assert(is_init_);

  int n_Element = exist_element_ids.size();
  int n_Exist_element = exist_element_ids.size();
  std::set<int> exist_node_ids;
  assert(n_Element>0 && n_Element<=frame_.sizeOfElementList());
  for(const int e_id : exist_element_ids)
  {
    assert(e_id>=0 && e_id<frame_.sizeOfElementList());
    int u_id = frame_.getElement(e_id)->endVertU()->id();
    int v_id = frame_.getElement(e_id)->endVertV()->id();
    exist_node_ids.insert(u_id);
    exist_node_ids.insert(v_id);
  }
  int n_Node = frame_.sizeOfVertList();

  if (verbose_)
  {
    create_k_.Start();
  }

  // turn exist_element_ids into a full dof free/fix/exist indicator
  // s[i] = 0 free, = 1 fixed, = -1 not exist
  std::vector<int> s(n_Node*6, -1);
  for(int e_id : exist_element_ids)
  {
    assert(e_id >= 0 && e_id < frame_.sizeOfElementList());
    auto e = frame_.getElement(e_id);

    std::vector<int> end_ids;
    end_ids.push_back(e->endVertU()->id());
    end_ids.push_back(e->endVertV()->id());

    for(auto id : end_ids)
    {
      if (frame_.getVert(id)->isFixed())
      {
        // TODO: assuming all dofs are fixed now
        for (int j = 0; j < 6; j++)
        {
          s[id * 6 + j] = 1;
        }
      }
      else
      {
        for (int j = 0; j < 6; j++)
        {
          s[id * 6 + j] = 0;
        }
      }
    }
  }

  // count supp_dof and res_dof to init K_slice
  int f_dof = std::count(s.begin(), s.end(), 0);
  int s_dof = std::count(s.begin(), s.end(), 1);
  int n_exist_dof = std::count(s.begin(), s.end(), -1);
  int dof = 6*n_Node;
  assert(s_dof + f_dof + n_exist_dof == dof);

  // mapping rearranged (R) dof_id into original (O) dof_id
  Eigen::VectorXi id_map_RO(dof);
  int f_tail = 0;
  int s_tail = f_dof;
  int n_exist_tail = f_dof + s_dof;
  for(int i=0; i<dof; i++)
  {
    if(-1 == s[i])
    {
      id_map_RO(n_exist_tail) = i;
      n_exist_tail++;
    }
    if(0 == s[i])
    {
      id_map_RO(f_tail) = i;
      f_tail++;
    }
    if(1 == s[i])
    {
      id_map_RO(s_tail) = i;
      s_tail++;
    }
  }

  Eigen::MatrixXd Perm(dof, dof);
  Perm.setZero();
  for(int i=0; i<dof; i++)
  {
    Perm(i, id_map_RO(i)) = 1;
  }
  auto Perm_inv = Perm.inverse();

  // permute the full stiffness matrix & carve the needed portion out
  createCompleteGlobalStiffnessMatrix(exist_element_ids);
  auto K_perm = Perm * K_assembled_full_ * Perm_inv;
  auto K_ff = K_perm.block(0,0,f_dof,f_dof);
  auto K_sf = K_perm.block(f_dof,0,s_dof,f_dof);

  Eigen::VectorXd P_perm = Perm * nodal_load_P_;

  auto P_ff = P_perm.segment(0, f_dof);

  if (verbose_)
  {
    create_k_.Stop();
  }

  Eigen::VectorXd U_ff(f_dof);
  if (!stiff_solver_.solveSystemLU(K_ff, P_ff, U_ff))
  {
    if(verbose_)
    {
      std::cout << "ERROR: Stiffness Solver fail!\n" << std::endl;
    }
    return false;
  }

  // reaction force
  Eigen::VectorXd R(dof);
  R.setZero();

  auto R_s = K_sf * U_ff;
  R.segment(f_dof, s_dof) = R_s;
  R = Perm_inv * R.eval();

  // assemble and permute the nodal displacement result back
  Eigen::VectorXd U(dof);
  U.setZero();

  U.segment(0, f_dof) = U_ff;
  U = Perm_inv * U.eval();

  // element internal reaction
  element_reaction.resize(n_Exist_element, 13);
  int fill_cnt = 0;
  for (const int e_id : exist_element_ids)
  {
    const auto e = frame_.getElement(e_id);
    const auto end_u = e->endVertU();
    const auto end_v = e->endVertV();
    assert(e_id == e->id());

    // transform to global frame
    Eigen::Matrix3d R_LG;
    getGlobal2LocalRotationMatrix(end_u->position(), end_v->position(), R_LG);

    Eigen::MatrixXd R_LG_diag(12, 12);
    R_LG_diag.setZero();
    for(int j=0;j<4;j++)
    {
      R_LG_diag.block<3, 3>(j*3, j*3) = R_LG;
    }

    Eigen::VectorXd Ue(12);
    Ue.segment(0,6) = U.segment(end_u->id()*6, 6);
    Ue.segment(6,6) = U.segment(end_v->id()*6, 6);

    element_reaction(fill_cnt,0) = e_id;
    element_reaction.block<1,12>(fill_cnt,1) = R_LG_diag * element_K_list_[e->id()] * Ue;

    fill_cnt++;
  }

  // assemble result into output formats
  // node displacement
  node_displ.resize(f_dof,3);
  fill_cnt=0;
  for(int i=0; i<n_Node; i++)
  {
    for(int j=0; j<6; j++)
    {
      if(0 == s[6*i+j])
      {
        assert(fill_cnt < f_dof);
        node_displ(fill_cnt,0) = i;
        node_displ(fill_cnt,1) = j;
        node_displ(fill_cnt,2) = U(6*i+j);

        fill_cnt++;
      }
    }
  }

  if(verbose_)
  {
    std::cout << "nodal disp (mm/rad):" << std::endl;
    std::cout << "node_id\tdof\tval" << std::endl;
    for(int i=0; i<node_displ.rows(); i++)
    {
      std::string dof_name;
      switch ((int)node_displ(i,1))
      {
        case 0: dof_name = "Tx";
          break;
        case 1: dof_name = "Ty";
          break;
        case 2: dof_name = "Tz";
          break;
        case 3: dof_name = "Rx";
          break;
        case 4: dof_name = "Ry";
          break;
        case 5: dof_name = "Rz";
          break;
      }
      std::cout << node_displ(i,0) << "\t"
                << node_displ(i,1) << "(" << dof_name << ")" << "\t"
                << node_displ(i,2) << std::endl;
    }
    std:cout << std::endl;
  }

  // fixities reaction
  fixities_reaction.resize(s_dof,3);
  fill_cnt=0;
  for(int i=0; i<n_Node; i++)
  {
    for(int j=0; j<6; j++)
    {
      if(1 == s[6*i+j])
      {
        assert(fill_cnt < s_dof);
        fixities_reaction(fill_cnt,0) = i;
        fixities_reaction(fill_cnt,1) = j;
        fixities_reaction(fill_cnt,2) = R(6*i+j);

        fill_cnt++;
      }
    }
  }

  if(verbose_)
  {
    std::cout << "fixities reaction (N, N-mm):" << std::endl;
    std::cout << "node_id\tdof\tval" << std::endl;
    for(int i=0; i<fixities_reaction.rows(); i++)
    {
      std::string dof_name;
      switch ((int)fixities_reaction(i,1))
      {
        case 0: dof_name = "Fx";
          break;
        case 1: dof_name = "Fy";
          break;
        case 2: dof_name = "Fz";
          break;
        case 3: dof_name = "Mx";
          break;
        case 4: dof_name = "My";
          break;
        case 5: dof_name = "Mz";
          break;
      }
      std::cout << fixities_reaction(i,0) << "\t"
                << fixities_reaction(i,1) << "(" << dof_name << ")" << "\t"
                << fixities_reaction(i,2) << std::endl;
    }
    std::cout << std::endl;
  }

  if(verbose_)
  {
    std::cout << "element reaction (N, N-mm):" << std::endl;
    std::cout << "e_id\tv_id\t" << std::endl;
    for(int i=0; i<element_reaction.rows(); i++)
    {
      std::cout << element_reaction(i,0) << "\t"
                << frame_.getElement(element_reaction(i,0))->endVertU()->id() << " ";

      std::cout << element_reaction.block<1,6>(i,1) << std::endl;

      std::cout << " " << "\t"
                << frame_.getElement(element_reaction(i,0))->endVertV()->id() << " ";

      std::cout << element_reaction.block<1,6>(i,7) << std::endl;
    }
    std::cout << std::endl;
  }

  if(verbose_)
  {
    printOutTimer();
  }

  // stiffness criteria check
  return checkStiffnessCriteria(node_displ, fixities_reaction, element_reaction);
}

bool Stiffness::solve(const std::vector<int>& exist_element_ids,
                      const bool& cond_num)
{
  Eigen::MatrixXd U,R,F;
  return solve(exist_element_ids, U, R, F, cond_num);
}

bool Stiffness::solve(Eigen::MatrixXd& node_displ,
                      Eigen::MatrixXd& fixities_reaction,
                      Eigen::MatrixXd& element_reaction,
                      const bool& cond_num)
{
  int nElement = frame_.sizeOfElementList();
  std::vector<int> all_e_ids;
  for(int i=0; i<nElement; i++)
  {
    all_e_ids.push_back(i);
  }
  return solve(
      all_e_ids, node_displ, fixities_reaction, element_reaction, cond_num);
}

bool Stiffness::solve(const bool& cond_num)
{
  Eigen::MatrixXd U,R,F;
  return solve(U, R, F, cond_num);
}

//bool Stiffness::CheckIllCondition(IllCondDetector &stiff_inspector)
//{
//	if (verbose_)
//	{
//		check_ill_.Start();
//	}
//
//	bool bSuccess = true;
//	double cond_num;
//	cond_num = stiff_inspector.ComputeCondNum();
//	printf("Condition Number = %9.3e\n", cond_num);
//	if (cond_num < MCOND_TOL)
//	{
//		if (verbose_)
//		{
//			printf(" < tol = %7.1e\n", MCOND_TOL);
//			printf(" * Acceptable Matrix! *\n");
//		}
//	}
//	else
//	{
//		printf(" > tol = %7.1e\n", MCOND_TOL);
//		printf(" * Ill Conditioned Stiffness Matrix! *\n");
//		printf("Press any key to exit...\n");
//		bSuccess = false;
//	}
//
//	if (verbose_)
//	{
//		check_ill_.Stop();
//	}
//
//	return bSuccess;
//}

bool Stiffness::checkStiffnessCriteria(const Eigen::MatrixXd& node_displ,
                                       const Eigen::MatrixXd& fixities_reaction,
                                       const Eigen::MatrixXd& element_reation)
{
  // stiffness check
  // nodal displacement check
  for(int i=0; i<node_displ.rows(); i++)
  {
    if(0<=node_displ(i,1) && node_displ(i,1)<=2)
    {
      if(std::abs(node_displ(i,2)) > transl_tol_)
      {
        if(verbose_)
        {
          std::cout << "node #" << node_displ(i,0)
                    << "translational dof #" << node_displ(i,1)
                    << " disp: " << node_displ(i,2)
                    << " > " << "tolerance " << transl_tol_ << std::endl;
        }
        return false;
      }
    }
    else
    {
      if(std::abs(node_displ(i,2)) > rot_tol_)
      {
        if(verbose_)
        {
          std::cout << "node #" << node_displ(i,0)
                    << "rotational dof #" << node_displ(i,1)
                    << " disp: " << node_displ(i,2)
                    << " > " << "tolerance " << rot_tol_ << std::endl;
        }
        return false;
      }
    }
  }

  // stability check
  // element reaction check
  // grounded element shouldn't be in tension
  for(int i=0; i<element_reation.rows(); i++)
  {
    int e_id = (int)element_reation(i,0);
    const auto e = frame_.getElement(e_id);
    if(e->endVertU()->isFixed() || e->endVertV()->isFixed())
    {
      // check tension
      if(element_reation(i,1) < 0)
      {
        assert(element_reation(i,7) > 0 && "element axial reaction sign not consistent!");

        if(verbose_)
        {
          std::cout << "grounded element #" << e_id
                    << " is in tension, the structure is not stable." << std::endl;
        }

        return false;
      }
    }
  }

  return true;
}


void Stiffness::printOutTimer()
{
  if (verbose_)
  {
    printf("***Stiffness timer result:\n");
    stiff_solver_.solve_timer_.Print("SolveK:");

    create_k_.Print("CreateGlobalK:");
    check_ill_.Print("CheckIllCond:");
  }
}

} // namespace stiffness_checker
} // namespace conmech
