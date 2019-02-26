  #include <cmath>
  #include <set>
  #include <algorithm>

  #include <Eigen/Dense>
  #include <Eigen/Sparse>

  #include "stiffness_checker/Util.h"
  #include "stiffness_checker/StiffnessSolver.h"
  #include "stiffness_checker/StiffnessIO.h"

  #include "stiffness_checker/Stiffness.h"

  namespace conmech
  {
  namespace stiffness_checker
  {

  // copy from an existing frame
  //Stiffness::Stiffness(Frame &frame, bool verbose, std::string model_type)
  //  : frame_(frame), verbose_(verbose), is_init_(false), transl_tol_(1.0), rot_tol_(5 * (3.14 / 180))
  //{
  //  // frame sanity check
  //  int N_vert = frame_.sizeOfVertList();
  //  int N_element = frame_.sizeOfElementList();
  //  assert(N_vert > 0 && N_element > 0);
  //
  //  stiff_solver_.timing_ = verbose;
  //  model_type_ = model_type;
  //
  //  init();
  //}

  Stiffness::Stiffness(const std::string& json_file_path, bool verbose, const std::string& model_type, bool output_json)
    : verbose_(verbose), is_init_(false), include_self_weight_load_(false),
      transl_tol_(1e-3), rot_tol_(1 * (3.14 / 180)), write_result_(output_json)
  {
    verbose_ = verbose;
    stiff_solver_.timing_ = verbose;

    // parse frame
    frame_.loadFromJson(json_file_path);

    // parse material properties
    parseMaterialPropertiesJson(json_file_path, material_parm_);

    model_type_ = model_type;

    output_json_file_name_ = "";
    output_json_file_path_ = "";

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
    assert(xyz_dof_id_.size() == node_dof_ * 2);

    createElementStiffnessMatrixList();
    createElementSelfWeightNodalLoad();

    // create id_map_
    int dof = frame_.sizeOfVertList() * node_dof_;
    int N_element = frame_.sizeOfElementList();

    id_map_.resize(N_element, node_dof_ * 2);
    auto lin_sp_id = Eigen::VectorXi::LinSpaced(node_dof_, 0, node_dof_ - 1); // 0,1,..,5

    for (int i = 0; i < N_element; i++)
    {
      const auto e = frame_.getElement(i);
      const auto end_u = e->endVertU();
      const auto end_v = e->endVertV();

      id_map_.block(i, 0, 1, node_dof_) =
        (Eigen::VectorXi::Constant(node_dof_, node_dof_ * end_u->id()) + lin_sp_id).transpose();
      id_map_.block(i, node_dof_, 1, node_dof_) =
        (Eigen::VectorXi::Constant(node_dof_, node_dof_ * end_v->id()) + lin_sp_id).transpose();
    }
   // std::cout << id_map_ << std::endl;

    // create fixities table
    int N_node = frame_.sizeOfVertList();
    assert(frame_.sizeOfFixedVert() && "at least one fixed vert in model");
    fixities_table_ = Eigen::MatrixXi::Zero(frame_.sizeOfFixedVert(), 1 + node_dof_);

    // TODO: this is ugly but..
    int cnt = 0;
    for (int i = 0; i < N_node; i++)
    {
      auto const v = frame_.getVert(i);
      // this should be true ...
      assert(v->id() == i);

      if (v->isFixed())
      {
        assert(v->fixities().size() == node_dof_);
        fixities_table_(cnt, 0) = v->id();
        const auto v_fs = v->fixities();
        for (int j = 1; j <= node_dof_; j++)
        {
          fixities_table_(cnt, j) = v_fs[j-1];
        }
        cnt++;
      }
    }
    // std::cout << fixities_table_ << std::endl;

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
    assert(transl_tol > 0 && rot_tol > 0 && "invalid tolerance: tolerance must be bigger than 0!");
    transl_tol_ = transl_tol;
    rot_tol_ = rot_tol;
  }

  // TODO: self-weight depends on the existing elements too!
  void Stiffness::createSelfWeightNodalLoad(const std::vector<int> &exist_e_ids, Eigen::VectorXd &self_weight_load_P)
  {
    assert(is_init_);
    assert(id_map_.cols() == 2 * node_dof_);

    int N_vert = frame_.sizeOfVertList();
    self_weight_load_P = Eigen::VectorXd::Zero(N_vert * node_dof_);

    // assuming solid circular cross sec for now
    double Ax = M_PI * std::pow(material_parm_.radius_, 2);

    for (const int e_id : exist_e_ids)
    {
      assert(0 <= e_id && e_id < element_gravity_nload_list_.size());
      auto Qe = element_gravity_nload_list_[e_id];

      assert(Qe.size() == id_map_.cols());

      for (int j = 0; j < id_map_.cols(); j++)
      {
        self_weight_load_P[id_map_(e_id, j)] += Qe[j];
      }
    }

  //  std::cout << self_weight_load_P << std::endl;
  }

  void Stiffness::createExternalNodalLoad(
    const Eigen::MatrixXd &nodal_forces, Eigen::VectorXd &ext_load)
  {
    assert(is_init_);
    assert(nodal_forces.cols() == node_dof_ + 1);

    int full_dof = node_dof_ * frame_.sizeOfVertList();
    ext_load.resize(full_dof);

    for (int i = 0; i < nodal_forces.rows(); i++)
    {
      int v_id = nodal_forces(i, 0);
      assert(0 <= v_id && v_id < frame_.sizeOfVertList());

      ext_load.segment(v_id * 6, 6) = nodal_forces.block<1, 6>(i, 1);
    }

  //  std::cout << "nodal_load" << ext_load << std::endl;
  }

  void Stiffness::setLoad(const Eigen::MatrixXd &nodal_forces)
  {
    assert(is_init_);

    int dof = node_dof_ * frame_.sizeOfVertList();
    bool is_empty_ext = nodal_forces.isZero(0) || nodal_forces.rows() == 0;

    assert(!is_empty_ext && "No load is assigned.");

    if (!is_empty_ext)
    {
      createExternalNodalLoad(nodal_forces, nodal_load_P_);
    }
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
    // cross section area: assuming solid circular, unit: m^2
    double A = M_PI * std::pow(material_parm_.radius_, 2);
  //  double Asy = Ax * (6 + 12 * material_parm_.poisson_ratio_ + 6 * std::pow(material_parm_.poisson_ratio_,2))
  //      / (7 + 12 * material_parm_.poisson_ratio_ + 4 * std::pow(material_parm_.poisson_ratio_,2));
  //  double Asz = Asy;

    // torsion constant (around local x axis)
    // for solid circle: (1/2)pi*r^4, unit: mm^4
    // see: https://en.wikipedia.org/wiki/Torsion_constant#Circle
    double Jx = 0.5 * M_PI * std::pow(material_parm_.radius_, 4);
  //  double Jx = (2.0/3.0) * pi * std::pow(material_parm_.radius_, 4);

    // area moment of inertia (bending about local y,z-axis)
    // assuming solid circular area of radius r, unit: m^4
    // see: https://en.wikipedia.org/wiki/List_of_area_moments_of_inertia
    // note this is slender rod of length L and Mass M, spinning around end
    double Iy = M_PI * std::pow(material_parm_.radius_, 4) / 4;
    double Iz = Iy;

    // E,G: MPa; mu (poisson ratio): unitless
    double E = material_parm_.youngs_modulus_;
    double G = material_parm_.shear_modulus_;
    double mu = material_parm_.poisson_ratio_;

    assert(std::abs((G - E / (2 * (1 + mu))) / G) < 1e-3
           && "input poisson ratio not compatible with shear and Young's modulus!");

    int N_element = frame_.sizeOfElementList();
    assert(N_element > 0);

    using namespace std;

    element_K_list_.reserve(N_element);
    rot_m_list_.reserve(N_element);

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
          assert(e_react_dof_id_(k) < R_LG_diag.rows() && xyz_dof_id_(kk) < R_LG_diag.cols());
          R_temp(k, kk) = R_LG_diag(e_react_dof_id_(k), xyz_dof_id_(kk));
        }
      }
      R_LG_diag = R_temp;

      for (int k = 0; k < e_react_dof_id_.size(); k++)
      {
        for (int kk = 0; kk < e_react_dof_id_.size(); kk++)
        {
          assert(e_react_dof_id_(k) < K_loc.rows() && e_react_dof_id_(kk) < K_loc.cols());
          K_loc_temp(k, kk) = K_loc(e_react_dof_id_(k), e_react_dof_id_(kk));
        }
      }
      K_loc = K_loc_temp;

      // DO NOT combine transpose
      //https://libigl.github.io/matlab-to-eigen.html
      auto R_LG_diagT = R_LG_diag.transpose();
      K_loc = R_LG_diagT * K_loc * R_LG_diag;

  //    std::cout <<  "K_loc:\n" << K_loc << std::endl;
  //    std::cout <<  "R_LG_diag:\n" << R_LG_diag << std::endl;

      element_K_list_.push_back(K_loc);
      rot_m_list_.push_back(R_LG_diag);
    }
  }

  void Stiffness::createElementSelfWeightNodalLoad()
  {
    // refer [MSA McGuire et al.] P111
    // Loads between nodal points
    int N_vert = frame_.sizeOfVertList();
    int N_element = frame_.sizeOfElementList();

    // assuming solid circular cross sec for now
    double Ax = M_PI * std::pow(material_parm_.radius_, 2);

    element_gravity_nload_list_.clear();
    element_gravity_nload_list_.reserve(N_element);

    for (int i = 0; i < N_element; i++)
    {
      const auto e = frame_.getElement(i);
      const auto end_u = e->endVertU();
      const auto end_v = e->endVertV();

      double Le = e->getLength();

      // uniform force density along the element
      // due to gravity
      double q_sw = material_parm_.density_ * Ax;

      if (model_type_ == "frame")
      {
        if (3 == dim_)
        {
          Eigen::Matrix3d R_eG;
          getGlobal2LocalRotationMatrix(
            e->endVertU()->position(),
            e->endVertV()->position(), R_eG);
          // gravity in local frame
          // https://github.mit.edu/yijiangh/frame3dd/blob/master/frame3dd_io.c#L905

          Eigen::VectorXd fixed_end_sw_load = Eigen::VectorXd::Zero(2 * node_dof_);

          // according to fixed-end force diagram
          // uniform load on a beam with both ends fixed ([MSA] p111.)
          // fixed-end fictitious reaction
          fixed_end_sw_load[2] = -q_sw * Le / 2.0;
          fixed_end_sw_load[3] =
            q_sw * std::pow(Le, 2) / 12.0 * ((-R_eG(1, 0) * R_eG(2, 2) + R_eG(1, 2) * R_eG(2, 0)) * (-1));
          fixed_end_sw_load[4] =
            q_sw * std::pow(Le, 2) / 12.0 * ((-R_eG(1, 1) * R_eG(2, 2) + R_eG(1, 2) * R_eG(2, 1)) * (-1));
          fixed_end_sw_load[5] = 0;

          fixed_end_sw_load[8] = -q_sw * Le / 2.0;
          fixed_end_sw_load[9] =
            -q_sw * std::pow(Le, 2) / 12.0 * ((-R_eG(1, 0) * R_eG(2, 2) + R_eG(1, 2) * R_eG(2, 0)) * (-1));
          fixed_end_sw_load[10] =
            -q_sw * std::pow(Le, 2) / 12.0 * ((-R_eG(1, 1) * R_eG(2, 2) + R_eG(1, 2) * R_eG(2, 1)) * (-1));
          fixed_end_sw_load[11] = 0;

          element_gravity_nload_list_.push_back(fixed_end_sw_load);
        }
        else
        {
          // TODO
          assert(false && "2D truss gravity not implemented yet.");
        }
      }
      else
      {
        // TODO
        assert(false && "truss gravity not implemented yet.");
      }
    }
  }

  void Stiffness::createCompleteGlobalStiffnessMatrix(const std::vector<int> &exist_e_ids)
  {
    assert(element_K_list_.size() > 0);
    assert(exist_e_ids.size() > 0);
    assert(id_map_.rows() >= exist_e_ids.size());

    int total_dof = frame_.sizeOfVertList() * node_dof_;
    int N_element = frame_.sizeOfElementList();

    K_assembled_full_.resize(total_dof, total_dof);
    std::vector<Eigen::Triplet<double>> K_triplets;

    for (const int e_id : exist_e_ids)
    {
      assert(e_id >= 0 && e_id < frame_.sizeOfElementList());

      const auto K_e = element_K_list_[e_id];
      assert(K_e.rows() == 2 * node_dof_ && K_e.cols() == 2 * node_dof_);

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
  }

  bool Stiffness::solve(
    const std::vector<int> &exist_element_ids,
    Eigen::MatrixXd &node_displ,
    Eigen::MatrixXd &fixities_reaction,
    Eigen::MatrixXd &element_reaction,
    const bool &cond_num)
  {
    assert(is_init_);

    int n_Element = frame_.sizeOfElementList();
    int n_Node = frame_.sizeOfVertList();

    if (verbose_)
    {
      create_k_.Start();
    }

    // start with assuming all dof does not exist (-1)
    int dof = node_dof_ * n_Node;
    Eigen::VectorXi full_f = Eigen::VectorXi::Constant(dof, -1);

    assert(exist_element_ids.size() > 0 && exist_element_ids.size() <= frame_.sizeOfElementList());
    std::set<int> sub_nodes_set;
    for (int e_id : exist_element_ids)
    {
      assert(e_id >= 0 && e_id < frame_.sizeOfElementList());
      auto e = frame_.getElement(e_id);
      sub_nodes_set.insert(e->endVertU()->id());
      sub_nodes_set.insert(e->endVertV()->id());

      // turn the existing node's dofs to free(0)
      full_f.segment(node_dof_ * e->endVertU()->id(), node_dof_) = Eigen::VectorXi::Constant(node_dof_, 0);
      full_f.segment(node_dof_ * e->endVertV()->id(), node_dof_) = Eigen::VectorXi::Constant(node_dof_, 0);
    }

    // Assemble the list of fixed DOFs fixedList, a matrix of size
    // 1-by-(number of fixities)
    int n_SubFixedNode = 0;

    int n_Fixities = 0;
    for (int i = 0; i < fixities_table_.rows(); i++)
    {
      int v_id = fixities_table_(i, 0);

      if (sub_nodes_set.end() != sub_nodes_set.find(v_id))
      {
        full_f.segment(node_dof_ * v_id, node_dof_) = (fixities_table_.block(i, 1, 1, node_dof_)).transpose();

        n_Fixities += full_f.segment(node_dof_ * v_id, node_dof_).sum();
        n_SubFixedNode++;
      }
    }

  //  std::cout << "full_f: " << full_f << std::endl;

    if (0 == n_Fixities)
    {
      if (verbose_)
      {
        std::cout << "Not stable: At least one node needs to be fixed in the considered substructure!" << std::endl;
      }
      return false;
    }

    // count supp_dof and res_dof to init K_slice
    int n_Nexist = node_dof_ * (frame_.sizeOfVertList() - sub_nodes_set.size());
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
  //  std::cout << "id_map_RO:\n" << id_map_RO << std::endl;

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

    if (include_self_weight_load_)
    {
      Eigen::VectorXd load_sw;
      createSelfWeightNodalLoad(exist_element_ids, load_sw);
      nodal_load_P_ += load_sw;
    }

    Eigen::VectorXd Q_perm = Perm * nodal_load_P_;

    auto Q_m = Q_perm.segment(0, n_Free);
    auto Q_f = Q_perm.segment(n_Free, n_Fixities);

  //  std::cout << "Q_m:\n" << Q_m << std::endl;

    if (verbose_)
    {
      create_k_.Stop();
    }

    Eigen::VectorXd U_m(n_Free);
    if (!stiff_solver_.solveSparseSimplicialLDLT(K_mm, Q_m, U_m))
    {
      if (verbose_)
      {
        std::cout << "ERROR: Stiffness Solver fail! The sub-structure contains mechanism.\n" << std::endl;
      }
      return false;
    }

  //  std::cout << "U_m:\n" << U_m << std::endl;

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
      element_reaction.block(cnt, 1, 1, e_react_dof_id_.size()) = (rot_m_list_[e_id] * element_K_list_[e_id] *
                                                                   Ue).transpose();
      cnt++;
    }

    if (verbose_)
    {
      std::cout << "element reaction (kN, kN-m):" << std::endl;
      std::cout << element_reaction << "\n" << std::endl;
    }

    // fixities reaction
    fixities_reaction = Eigen::MatrixXd::Zero(n_SubFixedNode, 1 + node_dof_);
    cnt = 0;
    for (int i = 0; i < fixities_table_.rows(); i++)
    {
      int fix_node_id = fixities_table_(i, 0);
      if (sub_nodes_set.end() != sub_nodes_set.find(fix_node_id))
      {
        fixities_reaction(cnt, 0) = fix_node_id;
        fixities_reaction.block(cnt, 1, 1, node_dof_) = (R.segment(fix_node_id * node_dof_, node_dof_)).transpose();
        cnt++;
      }
    }

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

    if (verbose_)
    {
      std::cout << "nodal disp (m/rad):" << std::endl;
      std::cout << node_displ << "\n" << std::endl;
    }

    if (verbose_)
    {
      printOutTimer();
    }

    if (write_result_)
    {
      write_output_json(frame_, node_displ, fixities_reaction, element_reaction,
      output_json_file_path_ + output_json_file_name_);
    }

    // stiffness criteria check
    return checkStiffnessCriteria(node_displ, fixities_reaction, element_reaction);
  }

  bool Stiffness::solve(const std::vector<int> &exist_element_ids,
                        const bool &cond_num)
  {
    Eigen::MatrixXd U, R, F;
    return solve(exist_element_ids, U, R, F, cond_num);
  }

  bool Stiffness::solve(Eigen::MatrixXd &node_displ,
                        Eigen::MatrixXd &fixities_reaction,
                        Eigen::MatrixXd &element_reaction,
                        const bool &cond_num)
  {
    int nElement = frame_.sizeOfElementList();
    std::vector<int> all_e_ids;
    for (int i = 0; i < nElement; i++)
    {
      all_e_ids.push_back(i);
    }
    return solve(
      all_e_ids, node_displ, fixities_reaction, element_reaction, cond_num);
  }

  bool Stiffness::solve(const bool &cond_num)
  {
    Eigen::MatrixXd U, R, F;
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

  bool Stiffness::checkStiffnessCriteria(const Eigen::MatrixXd &node_displ,
                                         const Eigen::MatrixXd &fixities_reaction,
                                         const Eigen::MatrixXd &element_reation)
  {
    // stiffness check
    // nodal displacement check
    for (int i = 0; i < node_displ.rows(); i++)
    {
      if (0 <= node_displ(i, 1) && node_displ(i, 1) <= 2)
      {
        if (std::abs(node_displ(i, 2)) > transl_tol_)
        {
          if (verbose_)
          {
            std::cout << "node #" << node_displ(i, 0)
                      << "translational dof #" << node_displ(i, 1)
                      << " disp: " << node_displ(i, 2)
                      << " > " << "tolerance " << transl_tol_ << std::endl;
          }
          return false;
        }
      }
      else
      {
        if (std::abs(node_displ(i, 2)) > rot_tol_)
        {
          if (verbose_)
          {
            std::cout << "node #" << node_displ(i, 0)
                      << "rotational dof #" << node_displ(i, 1)
                      << " disp: " << node_displ(i, 2)
                      << " > " << "tolerance " << rot_tol_ << std::endl;
          }
          return false;
        }
      }
    }

    // stability check
    // element reaction check
    // grounded element shouldn't be in tension
  //  for (int i = 0; i < element_reation.rows(); i++)
  //  {
  //    int e_id = (int) element_reation(i, 0);
  //    const auto e = frame_.getElement(e_id);
  //    if (e->endVertU()->isFixed() || e->endVertV()->isFixed())
  //    {
  //      // check tension
  //      if (element_reation(i, 1) < 0)
  //      {
  //        assert(element_reation(i, 7) > 0 && "element axial reaction sign not consistent!");
  //
  //        if (verbose_)
  //        {
  //          std::cout << "grounded element #" << e_id
  //                    << " is in tension, the structure is not stable." << std::endl;
  //        }
  //
  //        return false;
  //      }
  //    }
  //  }

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
