#define _USE_MATH_DEFINES
#include <cstdlib>
#include <cmath>
#include <set>
#include <algorithm>
#include <numeric>
#include <chrono>
#include <ctime>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "stiffness_checker/Util.h"
#include "stiffness_checker/StiffnessSolver.h"
#include "stiffness_checker/StiffnessIO.h"

#include "stiffness_checker/Stiffness.h"

namespace{
const std::string PathSeparator =
#ifdef _WIN32
"\\";
#else
"/";
#endif
}

namespace conmech
{
namespace stiffness_checker
{

Stiffness::Stiffness(const Eigen::MatrixXd& V, const Eigen::MatrixXi& E, const Eigen::VectorXi& Fixities, 
                     const std::vector<std::string>& materials,
                     const bool& verbose, const std::string& model_type, const bool& output_json)
  : Vertices_(V), Elements_(E), Fixities_(Fixities),
    verbose_(verbose), is_init_(false), include_self_weight_load_(false),
    transl_tol_(1e-3), rot_tol_(3 * (3.14 / 180)), write_result_(output_json),
    has_stored_deformation_(false), stored_compliance_(-1.0), stored_alt_compliance_(-1.0)
{
  verbose_ = verbose;
  stiff_solver_.timing_ = verbose;

  // TODO: convert to parse json from strings
  // parseMaterialPropertiesJson(json_file_path, material_parms_);

  model_type_ = model_type;

  output_json_file_name_ = "";
  output_json_file_path_ = "";

  init();
}

Stiffness::~Stiffness()
{
}

bool Stiffness::init()
{
  dim_ = this->Vertices_.cols();

  // TODO: generalize to 2D
  assert(dim_ = 3 && "only support 3D structure now!");

  if (3 == dim_) {
    gravity_direction_ = Eigen::VectorXd(3);
    gravity_direction_ << 0, 0, -1;
  } else {
    gravity_direction_ = Eigen::VectorXd(2);
    gravity_direction_ << 0, -1;
  }

  // set up dimension, node_dof
  if (3 == dim_)
  {
    // 3D case, x, y, z, xx, yy, zz
    full_node_dof_ = 6;
    if (model_type_ == "frame")
    {
      node_dof_ = 6;
      xyz_dof_id_ = Eigen::VectorXi::LinSpaced(node_dof_ * 2, 0, node_dof_ * 2 - 1);
      e_react_dof_id_ = Eigen::VectorXi::LinSpaced(node_dof_ * 2, 0, node_dof_ * 2 - 1);
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
      xyz_dof_id_ = Eigen::VectorXi::LinSpaced(node_dof_ * 2, 0, node_dof_ * 2 - 1);
      e_react_dof_id_ = Eigen::VectorXi::LinSpaced(node_dof_ * 2, 0, node_dof_ * 2 - 1);
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
  // assert(xyz_dof_id_.size() == node_dof_ * 2);

  // init per-element material properties

  precomputeElementStiffnessMatrixList();
  precomputeElementSelfWeightLumpedLoad();

  // create id_map_
  int dof = nV() * node_dof_;
  int N_element = nE();
  int N_node = nV();

  id_map_.resize(N_element, node_dof_ * 2);
  v_id_map_.resize(N_node, node_dof_);
  auto lin_sp_id = Eigen::VectorXi::LinSpaced(node_dof_, 0, node_dof_ - 1); // 0,1,..,5

  for (int i = 0; i < N_node; i++)
  {
    v_id_map_.block(i, 0, 1, node_dof_) =
      (Eigen::VectorXi::Constant(node_dof_, node_dof_ * i) + lin_sp_id).transpose();
  }

  for (int i = 0; i < N_element; i++)
  {
    const auto end_u_id = Elements_(i, 0);
    const auto end_v_id = Elements_(i, 1);

    id_map_.block(i, 0, 1, node_dof_) =
      (Eigen::VectorXi::Constant(node_dof_, node_dof_ * end_u_id) + lin_sp_id).transpose();
    id_map_.block(i, node_dof_, 1, node_dof_) =
      (Eigen::VectorXi::Constant(node_dof_, node_dof_ * end_v_id) + lin_sp_id).transpose();
  }

  // create zero load
  nodal_load_P_ = Eigen::VectorXd::Zero(dof);

  if (verbose_)
  {
    std::cout << "Initialization done" << std::endl;
  }

  is_init_ = true;

  return true;
}

void Stiffness::setNodalDisplacementTolerance(double transl_tol, double rot_tol)
{
  if(transl_tol <= 0 || rot_tol <= 0) {
    throw std::runtime_error("invalid tolerance: tolerance must be bigger than 0!");
  }
  transl_tol_ = transl_tol;
  rot_tol_ = rot_tol;
}

void Stiffness::precomputeElementStiffnessMatrixList()
{
  using namespace std;
  // TODO: clean the units in the comments here
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

  int N_element = nE();

  element_K_list_.clear();
  rot_m_list_.clear();
  element_K_list_.reserve(N_element);
  rot_m_list_.reserve(N_element);

  for (int i = 0; i < N_element; i++)
  {
    int end_u_id = Elements_(i, 0);
    int end_v_id = Elements_(i, 1);
    auto end_u = Vertices_.block<1, 3>(end_u_id, 1);
    auto end_v = Vertices_.block<1, 3>(end_v_id, 1);

    // element length, unit: m
    double L = (end_u - end_v).norm();

    // material properties
    double A = material_parms_[i].cross_sec_area_;
    double Jx = material_parms_[i].Jx_;
    double Iy = material_parms_[i].Iy_;
    double Iz = material_parms_[i].Iz_;

    // E,G: MPa; mu (poisson ratio): unitless
    double E = material_parms_[i].youngs_modulus_;
    double mu = material_parms_[i].poisson_ratio_;
    double G = material_parms_[i].shear_modulus_;

    // assert(std::abs((material_parm_.shear_modulus_ - G) / G) < 1e-6
    //        && "input poisson ratio not compatible with shear and Young's modulus!");

    Eigen::Matrix3d R_LG;
    getGlobal2LocalRotationMatrix(end_u, end_v, R_LG);

    // element stiffness matrix in local frame
    Eigen::MatrixXd K_loc;
    createLocalStiffnessMatrix(L, A, dim_, Jx, Iy, Iz, E, G, mu, K_loc);

    // transform to global frame
    Eigen::MatrixXd R_LG_diag = Eigen::MatrixXd::Zero(full_node_dof_ * 2, full_node_dof_ * 2);

    // 4 block or 2 block
    for (int j = 0; j < (full_node_dof_ / 3) * 2; j++)
    {
      R_LG_diag.block<3, 3>(j * 3, j * 3) = R_LG;
    }

    Eigen::MatrixXd R_temp = Eigen::MatrixXd::Zero(e_react_dof_id_.size(), xyz_dof_id_.size());
    Eigen::MatrixXd K_loc_temp = Eigen::MatrixXd::Zero(e_react_dof_id_.size(), e_react_dof_id_.size());

    for (int k = 0; k < e_react_dof_id_.size(); k++)
    {
      for (int kk = 0; kk < xyz_dof_id_.size(); kk++)
      {
        // assert(e_react_dof_id_(k) < R_LG_diag.rows() && xyz_dof_id_(kk) < R_LG_diag.cols());
        R_temp(k, kk) = R_LG_diag(e_react_dof_id_(k), xyz_dof_id_(kk));
      }
    }
    R_LG_diag = R_temp;

    for (int k = 0; k < e_react_dof_id_.size(); k++)
    {
      for (int kk = 0; kk < e_react_dof_id_.size(); kk++)
      {
        // assert(e_react_dof_id_(k) < K_loc.rows() && e_react_dof_id_(kk) < K_loc.cols());
        K_loc_temp(k, kk) = K_loc(e_react_dof_id_(k), e_react_dof_id_(kk));
      }
    }
    K_loc = K_loc_temp;

    // DO NOT combine transpose
    //https://libigl.github.io/matlab-to-eigen.html
    auto R_LG_diagT = R_LG_diag.transpose();
    K_loc = R_LG_diagT * K_loc * R_LG_diag;

    element_K_list_.push_back(K_loc);
    rot_m_list_.push_back(R_LG_diag);
  }
}

void Stiffness::createCompleteGlobalStiffnessMatrix(const std::vector<int> &exist_e_ids)
{
  // std::clock_t c_start = std::clock();

  // assert(element_K_list_.size() > 0);
  // assert(exist_e_ids.size() > 0);
  // assert(id_map_.rows() >= exist_e_ids.size());

  int total_dof = nV() * node_dof_;
  int N_element = nE();

  K_assembled_full_.resize(total_dof, total_dof);
  std::vector<Eigen::Triplet<double>> K_triplets;

  for (const int e_id : exist_e_ids)
  {
    // assert(e_id >= 0 && e_id < nE());

    const auto K_e = element_K_list_[e_id];
    // assert(K_e.rows() == 2 * node_dof_ && K_e.cols() == 2 * node_dof_);

    for (int i = 0; i < 2 * node_dof_; i++)
    {
      int row_id = id_map_(e_id, i);

      for (int j = 0; j < 2 * node_dof_; j++)
      {
        int col_id = id_map_(e_id, j);
        K_triplets.push_back(Eigen::Triplet<double>(row_id, col_id, K_e(i, j)));
      }
    }
  } // end e_id
  K_assembled_full_.setFromTriplets(K_triplets.begin(), K_triplets.end());

  // std::clock_t c_end = std::clock();
  // std::cout << "create global stiffness matrix: " << 1e3 * (c_end - c_start) / CLOCKS_PER_SEC << " ms" << std::endl;
}

void Stiffness::setLoad(const Eigen::MatrixXd &nodal_forces)
{
  // assert(is_init_);

  int dof = node_dof_ * nV();
  bool is_empty_ext = nodal_forces.isZero(0) || nodal_forces.rows() == 0;

  if (!is_empty_ext)
  {
    createExternalNodalLoad(nodal_forces, nodal_load_P_);
  } else {
    throw std::runtime_error("No load is assigned.");
  }
}

void Stiffness::setUniformlyDistributedLoad(const Eigen::MatrixXd &element_load_density)
{
  // assert(is_init_);
  bool is_empty_ext = element_load_density.isZero(0) || element_load_density.rows() == 0;
  if (!is_empty_ext)
  {
    precomputeElementUniformlyDistributedLumpedLoad(element_load_density);
  } else {
    throw std::runtime_error("No load is assigned.");
  }
}

void Stiffness::setGravityDirection(const Eigen::VectorXd& gravity_direction) {
    if (gravity_direction.size() == dim_) {
      gravity_direction_ = gravity_direction;
      precomputeElementSelfWeightLumpedLoad();
    }
  }

void Stiffness::createExternalNodalLoad(
  const Eigen::MatrixXd &nodal_forces, Eigen::VectorXd &ext_load)
{
  // assert(is_init_);
  // assert(nodal_forces.cols() == node_dof_ + 1);

  int full_dof = node_dof_ * nV();
  ext_load.resize(full_dof);

  for (int i = 0; i < nodal_forces.rows(); i++)
  {
    int v_id = nodal_forces(i, 0);
    // assert(0 <= v_id && v_id < nV());

    ext_load.segment(v_id * 6, 6) = nodal_forces.block<1, 6>(i, 1);
  }
}

void Stiffness::computeLumpedUniformlyDistributedLoad(const Eigen::Vector3d &w_G, const Eigen::Matrix3d &R_LG, const double &Le, 
  Eigen::VectorXd &fixed_end_lumped_load)
{
  fixed_end_lumped_load = Eigen::VectorXd::Zero(2 * node_dof_);

  // node 0,1 force
  fixed_end_lumped_load.segment(0, 3) = - w_G * Le / 2.0;
  fixed_end_lumped_load.segment(6, 3) = - w_G * Le / 2.0;

  // transform global load density to local density
  Eigen::Vector3d w_l = - R_LG * w_G;

  // node 0, 1 local moment
  Eigen::Vector3d M_0_l(3);
  M_0_l << 0, - w_l(2)*std::pow(Le,2)/12.0, w_l(1)*std::pow(Le,2)/12.0;
  Eigen::Vector3d M_1_l = - M_0_l;

  // transform local moment back to global
  fixed_end_lumped_load.segment(3, 3) = R_LG.transpose() * M_0_l;
  fixed_end_lumped_load.segment(9, 3) = R_LG.transpose() * M_1_l;

  // equivalent nodal load P_e = - P_fixed_end
  fixed_end_lumped_load *= -1;
}

void Stiffness::precomputeElementUniformlyDistributedLumpedLoad(const Eigen::MatrixXd &element_load_density)
{
  // refer [MSA McGuire et al.] P111
  // Loads between nodal points
  int N_vert = nV();
  int N_element = nE();

  element_lumped_nload_list_.resize(N_element);
  for (int i = 0; i < N_element; i++)
  {
    element_lumped_nload_list_[i] = Eigen::VectorXd::Zero(6);
  }

  for (int i = 0; i < element_load_density.rows(); i++)
  {
    const int e_id = element_load_density(i, 0);
    int end_u_id = Elements_(e_id, 0);
    int end_v_id = Elements_(e_id, 1);
    auto end_u = Vertices_.block<1, 3>(end_u_id, 1);
    auto end_v = Vertices_.block<1, 3>(end_v_id, 1);
    double Le = (end_u - end_v).norm();

    Eigen::Vector3d w_g = element_load_density.block<1, 3>(i, 1);
    
    Eigen::Matrix3d R_LG;
    getGlobal2LocalRotationMatrix(end_u, end_v, R_LG);

    Eigen::VectorXd fixed_end_lumped_load;
    computeLumpedUniformlyDistributedLoad(w_g, R_LG, Le, fixed_end_lumped_load);

    element_lumped_nload_list_[e_id] = fixed_end_lumped_load;
  }
}

void Stiffness::createUniformlyDistributedLumpedLoad(const std::vector<int>& exist_e_ids, Eigen::VectorXd &ext_load)
{
  // solve-time calculation
  int N_vert = nV();
  ext_load = Eigen::VectorXd::Zero(N_vert * node_dof_);

  if (element_lumped_nload_list_.size() != nE()) {
    return;
  }

  // assert(is_init_);
  // assert(id_map_.cols() == 2 * node_dof_);

  for (const int e_id : exist_e_ids) {
    // assert(0 <= e_id && e_id < element_lumped_nload_list_.size());
    auto Qe = element_lumped_nload_list_[e_id];

    for (int j = 0; j < id_map_.cols(); j++) {
      ext_load[id_map_(e_id, j)] += Qe[j];
    }
  }
}

void Stiffness::precomputeElementSelfWeightLumpedLoad()
{
  // TODO: reuse element uniformly distributed load?
  // gravity - Precomputation
  // refer [MSA McGuire et al.] P111
  // Loads between nodal points
  int N_vert = nV();
  int N_element = nE();

  element_gravity_nload_list_.clear();
  element_gravity_nload_list_.resize(N_element);

  for (int i = 0; i < N_element; i++)
  {
    int end_u_id = Elements_(i, 0);
    int end_v_id = Elements_(i, 1);
    auto end_u = Vertices_.block<1, 3>(end_u_id, 1);
    auto end_v = Vertices_.block<1, 3>(end_v_id, 1);
    double Le = (end_u - end_v).norm();

    // uniform force density along the element
    // due to gravity
    // density kN / m^3 * m^2
    double q_sw = material_parms_[i].density_ * material_parms_[i].cross_sec_area_;

    if (model_type_ == "frame")
    {
      if (3 == dim_)
      {
        // Eigen::Vector3d w_g(3);
        // w_g << 0, 0, -q_sw;
        Eigen::Vector3d w_g = gravity_direction_ * q_sw;
        // std::cout << "gravity force: " << w_g.transpose() << std::endl;
    
        Eigen::Matrix3d R_LG;
        getGlobal2LocalRotationMatrix(end_u, end_v, R_LG);

        Eigen::VectorXd fixed_end_lumped_load;
        computeLumpedUniformlyDistributedLoad(w_g, R_LG, Le, fixed_end_lumped_load);

        element_gravity_nload_list_[i] = fixed_end_lumped_load;
      }
      else
      {
        // TODO
        assert(false && "2D frame gravity not implemented yet.");
      }
    }
    else
    {
      // TODO
      assert(false && "truss gravity not implemented yet.");
    }
  }
}

void Stiffness::createSelfWeightLumpedLoad(const std::vector<int> &exist_e_ids, Eigen::VectorXd &self_weight_load_P)
{
  // solve-time calculation
  // assert(is_init_);
  // assert(id_map_.cols() == 2 * node_dof_);

  int N_vert = nV();
  self_weight_load_P = Eigen::VectorXd::Zero(N_vert * node_dof_);

  for (const int e_id : exist_e_ids)
  {
    // assert(0 <= e_id && e_id < element_gravity_nload_list_.size());
    auto Qe = element_gravity_nload_list_[e_id];
    // assert(Qe.size() == id_map_.cols());

    for (int j = 0; j < id_map_.cols(); j++)
    {
      self_weight_load_P[id_map_(e_id, j)] += Qe[j];
    }
  }
}

bool Stiffness::solve(
  const std::vector<int> &exist_element_ids,
  Eigen::MatrixXd &node_displ,
  Eigen::MatrixXd &fixities_reaction,
  Eigen::MatrixXd &element_reaction,
  const bool &cond_num)
{
  int n_Element = nE();
  int n_Node = nV();

  if (verbose_)
  {
    create_k_.Start();
  }

  // start with assuming all dof does not exist (-1)
  int dof = node_dof_ * n_Node;
  Eigen::VectorXi full_f = Eigen::VectorXi::Constant(dof, -1);

  // try {
      if(!(exist_element_ids.size() > 0 && exist_element_ids.size() <= n_Element)) {
          throw std::invalid_argument("input existing ids not within range!");
      }
  // } catch (const std::invalid_argument& e) {
  //   fprintf(stderr, "%s", e.what());
  //   return false;
  // }

  std::set<int> sub_nodes_set;
  for (int e_id : exist_element_ids)
  {
    // assert(e_id >= 0 && e_id < nE());
    int end_u_id = Elements_(e_id, 0);
    int end_v_id = Elements_(e_id, 1);
    sub_nodes_set.insert(end_u_id);
    sub_nodes_set.insert(end_v_id);

    // turn the existing node's dofs to free(0)
    full_f.segment(node_dof_ * end_u_id, node_dof_) = Eigen::VectorXi::Constant(node_dof_, 0);
    full_f.segment(node_dof_ * end_v_id, node_dof_) = Eigen::VectorXi::Constant(node_dof_, 0);
  }

  // Assemble the list of fixed DOFs fixedList, a matrix of size
  // 1-by-(number of fixities)
  int n_SubFixedNode = 0;
  int n_Fixities = 0;
  for (int i = 0; i < Fixities_.rows(); i++)
  {
    int v_id = Fixities_(i, 0);

    if (sub_nodes_set.end() != sub_nodes_set.find(v_id))
    {
      full_f.segment(node_dof_ * v_id, node_dof_) = (Fixities_.block(i, 1, 1, node_dof_)).transpose();

      n_Fixities += full_f.segment(node_dof_ * v_id, node_dof_).sum();
      n_SubFixedNode++;
    }
  }

  if (0 == n_Fixities)
  {
    throw std::runtime_error("Not stable: At least one node needs to be fixed in the considered substructure!");
  }

  // count supp_dof and res_dof to init K_slice
  int n_Nexist = node_dof_ * (n_Node - sub_nodes_set.size());
  int n_Free = dof - n_Fixities - n_Nexist;

  // generate permute id map
  int free_tail = 0;
  int fix_tail = n_Free;
  int nexist_tail = n_Free + n_Fixities;
  Eigen::VectorXi id_map_RO = Eigen::VectorXi::LinSpaced(dof, 0, dof - 1);
  for (int i = 0; i < dof; i++)
  {
    if (0 == full_f[i])
    {
      id_map_RO[free_tail] = i;
      free_tail++;
    }
    if (1 == full_f[i])
    {
      id_map_RO[fix_tail] = i;
      fix_tail++;
    }
    if (-1 == full_f[i])
    {
      id_map_RO[nexist_tail] = i;
      nexist_tail++;
    }
  }

  // a row permuatation matrix (multiply left)
  Eigen::SparseMatrix<double> Perm(dof, dof);
  std::vector<Eigen::Triplet<double>> perm_triplets;
  for (int i = 0; i < dof; i++)
  {
    perm_triplets.push_back(Eigen::Triplet<double>(i, id_map_RO[i], 1));
  }
  Perm.setFromTriplets(perm_triplets.begin(), perm_triplets.end());
  // perm operator on the right (column) = inverse of perm operation on the left (row)
  auto Perm_T = Perm.transpose();

  // permute the full stiffness matrix & carve the needed portion out
  createCompleteGlobalStiffnessMatrix(exist_element_ids);

  auto K_perm = Perm * K_assembled_full_ * Perm_T;
  auto K_mm = K_perm.block(0, 0, n_Free, n_Free);
  auto K_fm = K_perm.block(n_Free, 0, n_Fixities, n_Free);

  Eigen::VectorXd nodal_load_P_tmp = nodal_load_P_;

  Eigen::VectorXd element_lumped_load;
  createUniformlyDistributedLumpedLoad(exist_element_ids, element_lumped_load);
  nodal_load_P_tmp += element_lumped_load;

  if (include_self_weight_load_)
  {
    Eigen::VectorXd load_sw;
    createSelfWeightLumpedLoad(exist_element_ids, load_sw);
    nodal_load_P_tmp += load_sw;
  }

  Eigen::VectorXd Q_perm = Perm * nodal_load_P_tmp;
  auto Q_m = Q_perm.segment(0, n_Free);
  auto Q_f = Q_perm.segment(n_Free, n_Fixities);

  if (verbose_)
  {
    create_k_.Stop();
  }

  Eigen::VectorXd U_m(n_Free);
  // TODO: tolerance here?
  if (nodal_load_P_tmp.isZero())
  {
    U_m.setZero();
  }
  else
  {
    // TODO: what's the best solve strategy?
    // TODO: Conjugate gradient, precondition? https://www.cs.cmu.edu/~baraff/papers/sig98.pdf
    if (!stiff_solver_.solveSparseSimplicialLDLT(K_mm, Q_m, U_m))
    {
      if (verbose_)
      {
        std::cout << "ERROR: Stiffness Solver fail! The sub-structure contains mechanism.\n" << std::endl;
      }
      return false;
    }
  }

  stored_compliance_ = 0.5 * U_m.dot(Q_m);

  // reaction force
  Eigen::VectorXd R(dof);
  R.setZero();

  auto R_f = K_fm * U_m - Q_f;
  R.segment(n_Free, n_Fixities) = R_f;
  R = Perm_T * R.eval();

//  std::cout << "R:\n" << R << std::endl;

  // displacement
  Eigen::VectorXd U(dof);
  U.setZero();
  U.segment(0, n_Free) = U_m;
  U = Perm_T * U.eval();

//  std::cout << "U:\n" << U << std::endl;

  // start raw computattion results conversion
  // element internal reaction
  element_reaction = Eigen::MatrixXd::Zero(exist_element_ids.size(), 1 + e_react_dof_id_.size());
  int cnt = 0;
  for (const int e_id : exist_element_ids)
  {
    Eigen::VectorXd Ue(2 * node_dof_);
    for (int k = 0; k < 2 * node_dof_; k++)
    {
      Ue[k] = U[id_map_(e_id, k)];
    }

    element_reaction(cnt, 0) = e_id;
    element_reaction.block(cnt, 1, 1, e_react_dof_id_.size()) =
      (rot_m_list_[e_id] * element_K_list_[e_id] * Ue).transpose();
    cnt++;
  }
  stored_element_reaction_ = element_reaction;

  if (verbose_)
  {
    std::cout << "element reaction (kN, kN-m):" << std::endl;
    std::cout << element_reaction << "\n" << std::endl;
  }

  // fixities reaction
  fixities_reaction = Eigen::MatrixXd::Zero(n_SubFixedNode, 1 + node_dof_);
  cnt = 0;
  for (int i = 0; i < Fixities_.rows(); i++)
  {
    int fix_node_id = Fixities_(i, 0);
    if (sub_nodes_set.end() != sub_nodes_set.find(fix_node_id))
    {
      fixities_reaction(cnt, 0) = fix_node_id;
      fixities_reaction.block(cnt, 1, 1, node_dof_) = (R.segment(fix_node_id * node_dof_, node_dof_)).transpose();
      cnt++;
    }
  }
  stored_fixities_reaction_ = fixities_reaction;

  if (verbose_)
  {
    std::cout << "fixities reaction (kN, kN-m):" << std::endl;
    std::cout << fixities_reaction << "\n" << std::endl;
  }

  // node displacement
  node_displ = Eigen::MatrixXd::Zero(sub_nodes_set.size(), 1 + node_dof_);
  cnt = 0;
  for (auto it = sub_nodes_set.begin(); it != sub_nodes_set.end(); ++it)
  {
    int v_id = *it;
    node_displ(cnt, 0) = v_id;
    node_displ.block(cnt, 1, 1, node_dof_) = (U.segment(v_id * node_dof_, node_dof_)).transpose();
    cnt++;
  }
  stored_nodal_deformation_ = node_displ;

  if (verbose_)
  {
    std::cout << "nodal disp (m/rad):" << std::endl;
    std::cout << node_displ << "\n" << std::endl;
  }

  if (verbose_)
  {
    printOutTimer();
  }

  // if (write_result_)
  // {
  //   // TODO: use os-dependent file separator here
  //   write_output_json(frame_, node_displ, fixities_reaction, element_reaction,
  //   output_json_file_path_ + PathSeparator + output_json_file_name_);
  // }

  has_stored_deformation_ = true;
  stored_existing_ids_ = exist_element_ids;

  // stiffness criteria check
  return checkStiffnessCriteria(node_displ, fixities_reaction, element_reaction);
} // end core solve function

// overloaded solve
bool Stiffness::solve(const std::vector<int> &exist_element_ids,
                      const bool &cond_num)
{
  Eigen::MatrixXd U, R, F;
  return solve(exist_element_ids, U, R, F, cond_num);
}

// overloaded solve
bool Stiffness::solve(Eigen::MatrixXd &node_displ,
                      Eigen::MatrixXd &fixities_reaction,
                      Eigen::MatrixXd &element_reaction,
                      const bool &cond_num)
{
  int nElement = nE();
  std::vector<int> all_e_ids;
  for (int i = 0; i < nElement; i++)
  {
    all_e_ids.push_back(i);
  }
  return solve(
    all_e_ids, node_displ, fixities_reaction, element_reaction, cond_num);
}

// overloaded solve
bool Stiffness::solve(const bool &cond_num)
{
  Eigen::MatrixXd U, R, F;
  return solve(U, R, F, cond_num);
}

bool Stiffness::getElementStiffnessMatrices(std::vector<Eigen::MatrixXd> &element_stiffness_mats)
{
  // try{
    if (!is_init_) {
      throw std::runtime_error("stiffness checker not inited yet.\n");
    }
    element_stiffness_mats = element_K_list_;
    return true;
  // }
  // catch (const std::runtime_error &e) {
  //   fprintf(stderr, "%s", e.what());
  //   return false;
  // }
}

bool Stiffness::getElementLocal2GlobalRotationMatrices(std::vector<Eigen::MatrixXd> &e_L2G_rot_mats)
{
  // try{
    if (!is_init_) {
      throw std::runtime_error("stiffness checker not inited yet.\n");
    }
    e_L2G_rot_mats = rot_m_list_;
    return true;
  // }
  // catch (const std::runtime_error &e) {
  //   fprintf(stderr, "%s", e.what());
  //   return false;
  // }
}

bool Stiffness::getSelfWeightNodalLoad(const std::vector<int>& exist_e_ids, Eigen::VectorXd& self_weight_load)
{
  // try {
      createSelfWeightLumpedLoad(exist_e_ids, self_weight_load);
      return true;
  // }
  // catch (const std::runtime_error &e) {
  //   fprintf(stderr, "%s", e.what());
  //   return false;
  // }
}

bool Stiffness::getUniformlyDistributedLumpedLoad(const std::vector<int>& exist_e_ids, Eigen::VectorXd& lumped_load)
{
  // try {
      createUniformlyDistributedLumpedLoad(exist_e_ids, lumped_load);
      return true;
  // }
  // catch (const std::runtime_error &e) {
  //   fprintf(stderr, "%s", e.what());
  //   return false;
  // }
}

int Stiffness::nE() const  
{
  // TODO: inline
  if (is_init_) {
    return this->Elements_.rows();
  }
  return 0;
}

int Stiffness::nV() const
{
  if (is_init_) {
    return this->Vertices_.rows();
  }
  return 0;
}

int Stiffness::nFixV() const
{
  if (is_init_) {
    return this->Fixities_.rows();
  }
  return 0;
}

int Stiffness::nFixDof() const
{
  if (is_init_) {
    return this->Vertices_.rows();
  }
  return 0;
}

int Stiffness::dim() const
{ 
  return this->dim_; 
}

bool Stiffness::getSolvedResults(Eigen::MatrixXd &node_displ,
                      Eigen::MatrixXd &fixities_reaction,
                      Eigen::MatrixXd &element_reaction,
                      bool &pass_criteria)
{
  if(!hasStoredResults()) {
    throw std::runtime_error("no stored result found.\n");
  }
  node_displ = stored_nodal_deformation_;
  fixities_reaction = stored_fixities_reaction_;
  element_reaction = stored_element_reaction_;
  pass_criteria = checkStiffnessCriteria(node_displ, fixities_reaction, element_reaction);
  return true;
}

bool Stiffness::getMaxNodalDeformation(double &max_trans, double &max_rot,
    int &max_trans_vid, int &max_rot_vid)
{
  if (!hasStoredResults()) {
    throw std::runtime_error("no stored result found.\n");
  }
  const int nNode = stored_nodal_deformation_.rows();
  int mt_i, mt_j, mr_i, mr_j;

  Eigen::VectorXd trans_norm(nNode);
  for(int i = 0; i < nNode; i++){
    for(int j = 1; j <= 3; j++){
      trans_norm[i] += std::pow(stored_nodal_deformation_(i, j), 2);
    }
    trans_norm[i] = sqrt(trans_norm[i]);
  }
  max_trans = trans_norm.maxCoeff(&mt_i);
  max_trans_vid = stored_nodal_deformation_(mt_i, 0);

  // max_trans = stored_nodal_deformation_.block(0, 1, nNode, 3).cwiseAbs().maxCoeff(&mt_i, &mt_j);
  // max_trans_vid = stored_nodal_deformation_(mt_i, 0);

  max_rot = stored_nodal_deformation_.block(0, 4, nNode, 3).cwiseAbs().maxCoeff(&mr_i, &mr_j);
  max_rot_vid = stored_nodal_deformation_(mr_j, 0);
  return true;
}

bool Stiffness::getSolvedCompliance(double &compliance)
{
  // try {
    if (!hasStoredResults()) {
      throw std::runtime_error("no stored result found.\n");
    }
    compliance = stored_compliance_;
    return true;
  // }
  // catch (const std::runtime_error &e) {
  //   compliance = 0;
  //   fprintf(stderr, "%s", e.what());
  //   return false;
  // }
}

bool Stiffness::getSolvedAltCompliance(double &compliance)
{
  // try {
    if (!hasStoredResults()) {
      throw std::runtime_error("no stored result found.\n");
    }
    compliance = stored_alt_compliance_;
    return true;
  // }
  // catch (const std::runtime_error &e) {
  //   compliance = 0;
  //   fprintf(stderr, "%s", e.what());
  //   return false;
  // }
}

Eigen::MatrixXd Stiffness::getOriginalShape(const int& disc, const bool& draw_full_shape)
{
  // assert(disc >= 1);
  Eigen::MatrixXd orig_beam;

  std::vector<int> draw_ids;
  if(!draw_full_shape && stored_existing_ids_.size() > 0)
  {
    draw_ids = stored_existing_ids_;
    orig_beam = Eigen::MatrixXd::Zero(stored_existing_ids_.size()*(disc+1), dim_);
  }
  else
  {
    int n_Element = nE();
    draw_ids = std::vector<int>(n_Element);
    std::iota(draw_ids.begin(), draw_ids.end(), 0);
    orig_beam = Eigen::MatrixXd::Zero(n_Element*(disc+1), dim_);
  }

  for(int i=0; i < draw_ids.size(); i++)
  {
    int end_u_id = Elements_(draw_ids[i], 0);
    int end_v_id = Elements_(draw_ids[i], 1);
    const auto end_u = Vertices_.block<1, 3>(end_u_id, 1);
    const auto end_v = Vertices_.block<1, 3>(end_v_id, 1);

    for(int k=0; k<disc+1; k++)
    {
      Eigen::VectorXd inter_pt = end_u + (double(k)/disc)*(end_v - end_u);
      orig_beam.block(i*(disc+1)+k,0,1,dim_) = inter_pt.transpose();
    }
  }

  return orig_beam;
}

Eigen::MatrixXd Stiffness::getDeformedShape(const double& exagg, const int& disc)
{
  if(!has_stored_deformation_ || stored_existing_ids_.empty())
  {
    throw std::runtime_error("Unable to compute deformed beam shape: no nodal deformation computed and stored!");
  }
  // assert(exagg >= 0 && disc >= 2);

  int nE = stored_existing_ids_.size();

  // TODO: incorporate 2D case
  Eigen::VectorXd d_end_u(node_dof_), d_end_v(node_dof_);
  Eigen::MatrixXd tmp_beam_disp;
  Eigen::MatrixXd beam_displ = Eigen::MatrixXd::Zero(nE*(disc+1), dim_);

  for(int i; i < stored_existing_ids_.size(); i++)
  {
    int end_u_id = Elements_(stored_existing_ids_[i], 0);
    int end_v_id = Elements_(stored_existing_ids_[i], 1);

    // TODO: check this node_dof and full_node_dof compatability
    int u_row_id, v_row_id = -1;

    for(int j=0; j<stored_nodal_deformation_.rows(); j++)
    {
      if(stored_nodal_deformation_(j,0) == end_u_id)
      {
        u_row_id = j;
      }
      if(stored_nodal_deformation_(j,0) == end_v_id)
      {
        v_row_id = j;
      }
    }
    if(u_row_id == -1 || v_row_id == -1) {
      throw std::runtime_error("existing node id not found in existing deformation matrix!");
    }

    d_end_u = stored_nodal_deformation_.row(u_row_id).segment(1, node_dof_);
    d_end_v = stored_nodal_deformation_.row(v_row_id).segment(1, node_dof_);

    auto end_u = Vertices_.block<1, 3>(end_u_id, 1);
    auto end_v = Vertices_.block<1, 3>(end_v_id, 1);
    computeCubicDeformedBeam(end_u, end_v,
                             d_end_u, d_end_v, exagg, disc,
                             tmp_beam_disp);

    beam_displ.block(i*(disc+1), 0, disc+1, dim_) = tmp_beam_disp;
  }

  return beam_displ;
}

bool Stiffness::computeCubicDeformedBeam(const Eigen::VectorXd& end_u, const Eigen::VectorXd& end_v,
                                         const Eigen::VectorXd& d_end_u, const Eigen::VectorXd& d_end_v,
                                         const double& exagg, const int& disc,
                                         Eigen::MatrixXd& BeamPolygon)
{
    // assert(exagg >= 0 && disc >= 2);
    double L = (end_v - end_u).norm();

    Eigen::Matrix3d R_LG;
    getGlobal2LocalRotationMatrix(end_u, end_v, R_LG);
    Eigen::MatrixXd R = Eigen::MatrixXd::Zero(full_node_dof_*2,full_node_dof_*2);
    for(int i=0; i<(full_node_dof_/3)*2; i++)
    {
      R.block(i*3, i*3, 3, 3) = R_LG;
    }

    // compute end deflection in local coordinates
    Eigen::VectorXd D_global(full_node_dof_*2), D_local(full_node_dof_*2);
    D_global << d_end_u, d_end_v;
    D_local = exagg * R * D_global;

    // curve-fitting problem for a cubic polynomial
    // coordinates in local coordinate frame (x-axis along linear element)
    double d_u_local_x = D_local[0];
    double d_v_local_x = D_local[full_node_dof_-1] + L;

    Eigen::MatrixXd A_m(4,4);
    A_m.block(0,0,1,4) << 1, d_u_local_x, std::pow(d_u_local_x,2), std::pow(d_u_local_x,3);
    A_m.block(1,0,1,4) << 1, d_v_local_x, std::pow(d_v_local_x,2), std::pow(d_v_local_x,3);
    A_m.block(2,0,1,4) << 0, 1,           2*d_u_local_x,   3*std::pow(d_u_local_x,2);
    A_m.block(3,0,1,4) << 0, 1,           2*d_v_local_x,   3*std::pow(d_v_local_x,2);
    auto A_M_decomp = A_m.colPivHouseholderQr();

    Eigen::VectorXd a_coeff, b_coeff;
    if(2 == dim_)
    {
      // TODO: 2d case
      assert(false && "not implemented!");
    }
    else
    {
      Eigen::VectorXd a_rhs(4);
      a_rhs << D_local[1], D_local[7], D_local[5], D_local[11];
      a_coeff = A_M_decomp.solve(a_rhs);

      Eigen::VectorXd b_rhs(4);
      b_rhs << D_local[2], D_local[8], -D_local[4], -D_local[10];
      b_coeff = A_M_decomp.solve(b_rhs);
    }

    // beam discretization, break the beam into disc polygons (thus, disc+1 vertices)
    // TODO: get this as input
    int n_poly_vert = disc+1;

    BeamPolygon = Eigen::MatrixXd::Zero(n_poly_vert, dim_);
    double d_u = D_local[0];
    double d_v = D_local[full_node_dof_];
    double i_x, i_y, i_z; // interpolated local x, y, z
    auto R_LG_decomp = R_LG.colPivHouseholderQr();
    Eigen::VectorXd idD_global(dim_);

    for(int i=0; i<n_poly_vert; i++)
    {
      i_x = d_u + i*(L + d_v - d_u)/disc;
      Eigen::VectorXd monomial_vec(4);
      monomial_vec << 1, i_x, std::pow(i_x,2), std::pow(i_x,3);

      i_y = a_coeff.dot(monomial_vec);
      if(3 == dim_)
      {
        i_z = b_coeff.dot(monomial_vec);

        Eigen::Vector3d idD_local;
        idD_local << i_x, i_y, i_z;
        idD_global = R_LG_decomp.solve(idD_local);
      }
      else
      {
        assert(false && "not implemented");
      }
      BeamPolygon.block(i,0,1,dim_) = (end_u + idD_global).transpose();
    }

    return true;
}

bool Stiffness::checkStiffnessCriteria(const Eigen::MatrixXd &node_displ,
                                       const Eigen::MatrixXd &fixities_reaction,
                                       const Eigen::MatrixXd &element_reation)
{
  // stiffness check
  // nodal displacement check
  const int nNode = node_displ.rows();

  Eigen::VectorXd trans_norm(nNode);
  for(int i = 0; i < nNode; i++){
    for(int j = 1; j <= 3; j++){
      trans_norm[i] += std::pow(node_displ(i, j), 2);
    }
    trans_norm[i] = sqrt(trans_norm[i]);
  }
  double max_trans = trans_norm.maxCoeff();
  if (max_trans > transl_tol_)
  {
    return false;
  }

  // const auto& max_trans = node_displ.block(0, 1, nNode, 3).cwiseAbs();
  // const auto& max_rot = node_displ.block(0, 4, nNode, 3).cwiseAbs();

  // if (verbose_)
  // {
  //   std::cout << "max translation deformation: " << max_trans.maxCoeff()
  //             << " / " << "tolerance " << transl_tol_ << std::endl;
  //   std::cout << "max rotational deformation: " << max_rot.maxCoeff()
  //             << " / " << "tolerance " << rot_tol_ << std::endl;
  // }
  // if (max_trans.maxCoeff() > transl_tol_)
  // {
  //   return false;
  // }
  // if (max_rot.maxCoeff() > rot_tol_)
  // {
  //   return false;
  // }

  // TODO: not well tested yet
  // stability check
  // element reaction check
  // grounded element shouldn't be in tension
  // double sec_mod = M_PI * std::pow(material_parm_.radius_, 3) / 4;
  // double area = M_PI * std::pow(material_parm_.radius_, 2);
  //
  // for (int i = 0; i < fixities_reaction.rows(); i++)
  // {
  //   int v_id = (int) fixities_reaction(i, 0);
  //   const auto vertF = frame_.getVert(v_id);
  //   const auto& nghdE = vertF->getNghdElement();
  //
  //   double eF_lx;
  //   for(const auto& e : nghdE) {
  //     if(e->isFixed()) {
  //       // find element reaction
  //       int e_result_id = -1;
  //
  //       for(int j = 0; j < element_reation.size(); j++) {
  //         if (e->id() == element_reation(j, 0)) {
  //           e_result_id = j;
  //           break;
  //         }
  //       }
  //
  //       if(e->endVertU()->id() == v_id) {
  //         eF_lx = element_reation(e_result_id, 1);
  //         // only want tensile
  //         if(eF_lx <= 0) {
  //           eF_lx = 0;
  //         }
  //       } else {
  //         eF_lx = element_reation(e_result_id, 7);
  //         // only want tensile
  //         if(eF_lx >= 0) {
  //           eF_lx = 0;
  //         }
  //       } // end pts in local axis
  //     } // e is fixed
  //   }
  //   std::cout << "result eM: " << eF_lx << std::endl;
  //
  //   // moment reaction at fixities
  //   Eigen::VectorXd fixM = fixities_reaction.block<1,2>(i, 4);
  //
  //   // tension/compression bending stress
  //   Eigen::VectorXd sigma_b = fixM / sec_mod;
  //   std::cout << "bending stess: " << sigma_b << std::endl;
  //
  //   // axial tensile stress
  //   Eigen::VectorXd sigma_a = (eF_lx / area) * Eigen::VectorXd::Ones(2);
  //   std::cout << "axial stess: " << sigma_a << std::endl;
  //
  //   // total tensile stress = bending tensile stress + axial tensile stress
  //   Eigen::VectorXd total_tensile = sigma_b + sigma_a;
  //   std::cout << "total stess: " << total_tensile << "\n/" << 0.1 * material_parm_.tensile_yeild_stress_ << std::endl;
  //
  //   // if  > 0.1 (tbd) * yield stress
  //   if(std::abs(total_tensile.maxCoeff()) > 0.1 * material_parm_.tensile_yeild_stress_) {
  //     return false;
  //   }
  // } // fixities_reaction

 // for (int i = 0; i < element_reation.rows(); i++)
 // {
 //   int e_id = (int) element_reation(i, 0);
 //   const auto e = frame_.getElement(e_id);
 //   if (e->endVertU()->isFixed() || e->endVertV()->isFixed())
 //   {
 //     // check tension
 //     if (element_reation(i, 1) < 0)
 //     {
 //       assert(element_reation(i, 7) > 0 && "element axial reaction sign not consistent!");
 //
 //       if (verbose_)
 //       {
 //         std::cout << "grounded element #" << e_id
 //                   << " is in tension, the structure is not stable." << std::endl;
 //       }
 //
 //       return false;
 //     }
 //   }
 // }

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
