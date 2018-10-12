#include <cmath>
#include <algorithm>
#include "stiffness_checker/Util.h"
#include "stiffness_checker/StiffnessSolver.h"
#include "stiffness_checker/Stiffness.h"

namespace conmech
{
namespace stiffness_checker
{

Stiffness::Stiffness(Frame& frame, bool verbose)
    : frame_(frame), verbose_(verbose), is_init_(false)
{
  // frame sanity check
  int N_vert = frame_.sizeOfVertList();
  int N_element = frame_.sizeOfElementList();
  assert(N_vert>0 && N_element>0);

  init();
}

void Stiffness::createSelfWeightNodalLoad()
{
  // refer [MSA McGuire et al.] P111
  // Loads between nodal points
  int N_vert = frame_.sizeOfVertList();
  int N_element = frame_.sizeOfElementList();
  self_weight_load_P_ = Eigen::MatrixXd::Zero(N_vert, 6);

  // assuming solid circular cross sec for now
  double Ax = M_PI * std::pow(parm_.radius_, 2);

  // gravitional acc unit: m/s^2
  // in global frame
  double gz = parm_.g_;

  for (int i = 0; i < N_element; i++)
  {
    const auto e = frame_.getElement(i);
    const auto end_u = e->endVertU();
    const auto end_v = e->endVertV();

    double L = e->getLength();
    // uniform force density along the element
    // due to gravity
    double q = parm_.density_ * Ax * L * gz;
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

    self_weight_load_P_.block<1,6>(end_u->id(),0) += P_fG;

    // moment on the other end needs to be negated for sign consistency
    P_fG.tail(3) *= -1;
    self_weight_load_P_.block<1,6>(end_v->id(),0) += P_fG;
  }
}

void Stiffness::createExternalNodalLoad(const Eigen::MatrixXd& nodal_forces)
{
  // verify input
  for(int i=0; i < nodal_forces.rows(); i++)
  {
    int v_id = nodal_forces[i,0];
    assert(0<=v_id && v_id<N_vert);
  }

  ext_load_P_ = nodal_forces;
}

bool Stiffness::setNodalLoad(const Eigen::MatrixXd& nodal_forces,
                             const bool& include_self_weight)
{
  int dof = 6*frame_.sizeOfVertList();

  auto sw_load = Eigen::VectorXd::Zero(dof);
  auto ext_load = Eigen::VectorXd::Zero(dof);
  nodal_load_P_.resize(dof);

  if(include_self_weight)
  {
    createSelfWeightNodalLoad(sw_load);
  }
  createExternalNodalLoad(nodal_forces, ext_load);

  nodal_load_P_ = sw_load + ext_load;
}

void Stiffness::createElementStiffnessMatrixList()
{
  // check the complete element stiffness matrix
  // in local frame: [MSA McGuire et al.] P73
  // for unit used in this function, we conform with
  // the unit convention used in [MSA]'s example:

  // for all cross section's geometrical properties
  // cross section: mm^2
  // Iz, Iy: mm^4
  // J: mm^4

  // element length L: mm
  // E, G: GPa (1Pa = 1 N/m^2 = 1e-9 kN/mm^2)
  // Force: kN
  // Moment: kN m

  double pi = atan(1)*4;
  double A = pi * std::pow(parm_.radius_,2);
//  double Asy = Ax * (6 + 12 * parm_.poisson_ratio_ + 6 * std::pow(parm_.poisson_ratio_,2))
//      / (7 + 12 * parm_.poisson_ratio_ + 4 * std::pow(parm_.poisson_ratio_,2));
//  double Asz = Asy;

  // torsion constant (around local x axis)
  // see: https://en.wikipedia.org/wiki/Torsion_constant#Circle
  double Jx = 0.5 * pi * std::pow(parm_.radius_, 4);

  // area moment of inertia (bending about local y,z-axis)
  // assuming solid circular area of radius r
  // see: https://en.wikipedia.org/wiki/List_of_area_moments_of_inertia
  // note this is slender rod of length L and Mass M, spinning around end
  double Iy = F_PI * std::pow(parm_.radius_,4) / 4;
  double Iz = Iyy;

  double E = parm_.youngs_modulus_;
  double mu = parm_.poisson_ratio_;
  double G = E/(2*(1+mu));

  int N_element = frame_.sizeOfElementList();
  assert(N_element);

  element_K_list_.resize(N_element);
  for (int i = 0; i < N_element; i++)
  {
    const auto e = frame_.getElement(i);
    const auto end_u = e->endVertU();
    const auto end_v = e->endVertV();

    // element stiffness matrix in local frame
    Eigen::MatrixXd K_eL(12, 12);
    K_eL.setZero();

    Eigen::Matrix3d R_LG;
    getGlobal2LocalRotationMatrix(end_u->position(), end_v->position(), R_LG);

    // see: [Matrix Structural Analysis, McGuire et al., 2rd edition]
    // P73 - eq(4.34)
    Eigen::MatrixXd K_block(6,6);
    K_block.setZero();
    Eigen::VectorXd diag(6);

    // block_00 and block_11
    K_block(1,5) = 6*Iz / std::pow(L,2);
    K_block(2,4) = - 6*Izy / std::pow(L,2);
    K_block = K_block.eval() + K_block.transpose().eval();
    K_eL.block<6,6>(0,0) = K_block;
    K_eL.block<6,6>(6,6) = -K_block;

    diag[0] = A/L;
    diag[1] = 12*Iz / std::pow(L,3);
    diag[2] = 12*Iy / std::pow(L,3);
    diag[3] = Jx / (2*(1+mu)*L);
    diag[4] = 4*Iy / L;
    diag[5] = 4*Iz / L;
    K_eL.block<6,6>(0,0) = diag.asDiagonal();
    K_eL.block<6,6>(6,6) = diag.asDiagonal();

    // block_01 and block_10
    K_block.setZero();

    K_block(1,5) = 6*Iz / std::pow(L,2);
    K_block(2,4) = - 6*Izy / std::pow(L,2);
    K_block = K_block.eval() - K_block.transpose().eval();
    K_eL.block<6,6>(0,6) = K_block;
    K_eL.block<6,6>(6,0) = -K_block;

    diag[0] = -A/L;
    diag[1] = -12*Iz / std::pow(L,3);
    diag[2] = -12*Iy / std::pow(L,3);
    diag[3] = -Jx / (2*(1+mu)*L);
    diag[4] = 2*Iy / L;
    diag[5] = 2*Iz / L;
    K_eL.block<6,6>(0,6) = diag.asDiagonal();
    K_eL.block<6,6>(6,0) = diag.asDiagonal();

    K_eL *= E;

    assert(K_eL.isDiagonal());

    // transform to global frame
    Eigen::MatrixXd R_LG_diag(12, 12);
    R_LG_diag.setZero();
    for(i=0;i<4;i++)
    {
      R_LG_diag.block<3, 3>(i*3, i*3) = R_LG;
    }

    element_K_list_[i] = R_LG_diag.transpose() * K_eL * R_LG_diag;
  }
}

void Stiffness::createCompleteGlobalStiffnessMatrix()
{
  if (verbose_)
  {
    create_k_.Start();
  }

  assert(element_K_list_.size() > 0);

  int dof = frame_.sizeOfVertList() * 6;
  int N_element = frame_.sizeOfElementList();

  Eigen::MatrixXi id_map(N_element, 12);
  const Eigen::VectorXi lin_sp_id = Eigen::LinSpaced(6,0,5); // 0,1,..,5

  for(int i=0; i<N_element; i++)
  {
    const auto e = frame_.getElement(i);
    const auto end_u = e->endVertU();
    const auto end_v = e->endVertV();

    id_map.block<1,6>(i,0) = Eigen::VectorXi::Ones(6)*end_u->id() + lin_sp_id;
    id_map.block<1,6>(i,6) = Eigen::VectorXi::Ones(6)*end_v->id() + lin_sp_id;
  }

  K_assembled_full_.resize(dof, dof);
  K_assembled_full_.setZero();
  for(int i = 0; i < N_element; i++)
  {
    const auto K_e = element_K_list_[i];
    assert(K_e.rows() == 12 && K_e.cols() == 12);

    for(int j=0; j < 12; j++)
    {
      int row_id = id_map[i][j];

      for(int k=0; k < 12; k++)
      {
        int col_id = id_map[i][k];

        K_assembled_full_(row_id, col_id) += K_e(j,k);
      }
    }
  }

  if (verbose_)
  {
    create_k_.Stop();
  }
}

bool Stiffness::init()
{
  createElementStiffnessMatrixList();
  createCompleteGlobalStiffnessMatrix();

  is_init_ = true;
}

bool Stiffness::solve(
    const std::vector<int>& exist_element_ids,
    Eigen::MatrixXd& node_displ,
    Eigen::MatrixXd& fixities_reaction,
    Eigen::MatrixXd& element_reation,
    const bool& cond_num)
{
  assert(is_init_);

  int n_Node = frame_.sizeOfVertList();
  int n_Element = frame_.sizeOfElementList();

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
          s(id * 6 + j) = 1;
        }
      }
      else
      {
        for (int j = 0; j < 6; j++)
        {
          s(id * 6 + j) = 0;
        }
      }
    }
  }

  // count supp_dof and res_dof to init K_slice
  const int f_dof = std::count(s.begin(), s.end(), 0);
  const int s_dof = std::count(s.begin(), s.end(), 1);
  const int n_exist_dof = std::count(s.begin(), s.end(), -1);
  const int dof = 6*n_Node;

  assert(s_dof + f_dof + n_exist_dof == dof);
  assert(K_assembled_full_.rows() == dof && K_assembled_full_.cols() == dof);

  // mapping rearranged (R) dof_id into original (O) dof_id
  Eigen::VectorXi id_map_RO(dof);
  int f_tail = 0;
  int s_tail = f_dof;
  int n_exist_tail = f_dof + s_dof;
  for(int i=0; i<dof; i++)
  {
    if(-1 == s[i])
    {
      id_map_RO[n_exist_tail] = i;
      n_exist_tail++;
    }
    if(0 == s[i])
    {
      id_map_RO[f_tail] = i;
      f_tail++;
    }
    if(1 == s[i])
    {
      id_map_RO[s_tail] = i;
      s_tail++;
    }
  }

  Eigen::MatrixXd Perm(dof, dof);
  Perm.setZero();
  for(int i=0; i<dof; i++)
  {
    Perm(i, id_map_RO(i)) = 1;
  }

  // permute the full stiffness matrix & carve the needed portion out
  auto K_perm = Perm * K_assembled_full_ * Perm.inverse();
  auto K_ff = K_perm.block<f_dof, f_dof>(0,0);
  auto K_sf = K_perm.block<s_dof, f_dof>(f_dof,0);

  Eigen::VectorXd P_perm = Perm * nodal_load_P_;
  auto P_ff = P_perm.segment(0, f_dof);

  if (verbose_)
  {
    create_k_.Stop();
  }

  Eigen::VectorXd U_ff(f_dof);
  if (!stiff_solver_.solveSystemLU(K_ff, P_ff, U_ff))
  {
    std::cout << "ERROR: Stiffness Solver fail!\n" << std::endl;
    return false;
  }

  // reaction force
  auto R = Eigen::VectorXd::Zero(dof);
  auto R_s = K_sf * U_ff;
  R.segment(f_dof, f_dof + s_dof) = R_s;
  R = Perm.inverse() * R.eval();

  // assemble and permute the result back
  auto U = Eigen::VectorXd::Zero(dof);
  U.segment(0, f_dof) = U_ff;
  U = Perm.inverse() * U.eval();

  // TODO: get element internal reaction force

  // stiffness criteria check
//	if (!CheckError(stiff_inspector, U, R, int_F))
//	{
//		return false;
//	}

  // assemble result into formats
  // node displacement
  node_displ.resize(f_dof,3);
  int fill_cnt=0;
  for(int i=0; i<n_Node; i++)
  {
    for(int j=0; j<6; j++)
    {
      if(0 == s[6*i+j])
      {
        node_displ(fill_cnt,0) = i;
        node_displ(fill_cnt,1) = j;
        node_displ(fill_cnt,2) = U(6*i+j);

        fill_cnt++;
      }
    }
  }

  // fixities reaction
  fixities_reaction.resize(s_dof,3);
  int fill_cnt=0;
  for(int i=0; i<n_Node; i++)
  {
    for(int j=0; j<6; j++)
    {
      if(1 == s[6*i+j])
      {
        node_displ(fill_cnt,0) = i;
        node_displ(fill_cnt,1) = j;
        node_displ(fill_cnt,2) = R(6*i+j);

        fill_cnt++;
      }
    }
  }

  // element internal reaction


  return true;
}

bool Stiffness::solve(const std::vector<int>& exist_element_ids,
                      const bool& cond_num)
{
  Eigen::MatrixXd U,R,F;
  return solve(exist_element_ids, U, R, F, cond_num);
}

bool Stiffness::solve(Eigen::MatrixXd& node_displ,
                      Eigen::MatrixXd& fixities_reaction,
                      Eigen::MatrixXd& element_reation,
                      const bool& cond_num)
{
  int nElement = frame_.sizeOfElementList();
  std::vector<int> all_e_ids;
  for(int i=0; i<nElement; i++)
  {
    all_e_ids.push_back(i);
  }
  return solve(
      all_e_ids, node_displ, fixities_reaction, element_reation, cond_num);
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
//
//
//bool Stiffness::CheckError(IllCondDetector &stiff_inspector, Eigen::VectorXd &D)
//{
//	if (verbose_)
//	{
//		check_error_.Start();
//	}
//
//	bool bSuccess = true;
//	double error = stiff_inspector.EquilibriumError(K_, D, F_);
//	if (verbose_)
//	{
//		printf("Root Mean Square (RMS) equilibrium error = %9.3e\n", error);
//	}
//	if (error < STIFF_TOL)
//	{
//		if (verbose_)
//		{
//			printf(" < tol = %7.1e\n", STIFF_TOL);
//			printf(" * Converged *\n");
//		}
//	}
//	else
//	{
//		printf(" > tol = %7.1e\n", STIFF_TOL);
//		printf(" !! Not Converged !!\n");
//		printf("Press any key to exit...\n");
//		bSuccess = false;
//	}
//
//	if (verbose_)
//	{
//		check_error_.Stop();
//	}
//
//	return bSuccess;
//}


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
