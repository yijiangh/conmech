#include <cmath>
#include <algorithm>
#include "stiffness_checker/Util.h"
#include "stiffness_checker/Stiffness.h"

namespace conmech
{
namespace stiffness_checker
{

Stiffness::Stiffness(Frame& frame, bool verbose)
    : frame_(frame), verbose_(verbose)
{
  int N_vert = frame_.sizeOfVertList();
  int N_element = frame_.sizeOfElementList();
  assert(N_vert>0 && N_element>0);

  self_weight_load_P_ = Eigen::MatrixXd::Zero(N_vert, 6);
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

void Stiffness::createLocalStiffnessMatrixList()
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


void Stiffness::createGlobalStiffnessMatrix(const std::vector<int>& exist_element_ids,
                                            Eigen::Sparse<double>& K_assembled,
                                            Eigen::MatrixX2i& id_map)
{
  if (verbose_)
  {
    create_k_.Start();
  }

  assert(exist_element_ids.size()>0);

  // get existing elements' fixities
  std::vector<int> exist_fix_node_ids;
  assert(fixities_.rows() > 0 && fixities_.cols() == 7);

  for(int i=0; i<fixities_.rows(); i++)
  {
    if(std::find(exist_element_ids.begin(), exist_element_ids.end(), fixities_(i,0)))
    {
      exist_fix_node_ids.push_back(fixities_(i,0));
    }
  }

  std::vector<Eigen::Triplet<double>> K_list;

  K_.resize(6 * num_free_wf_nodes, 6 * num_free_wf_nodes);
  for (int i = 0; i < num_wf_edges; i++)
  {
    WF_edge *ei = ptr_frame->GetEdge(ptr_dualgraph_->e_orig_id(i));
    WF_edge *ej = ei->ppair_;
    int u = ej->pvert_->ID();
    int v = ei->pvert_->ID();
    int dual_u = ptr_dualgraph_->v_dual_id(u);
    int dual_v = ptr_dualgraph_->v_dual_id(v);

    if (dual_u < num_free_wf_nodes && dual_v < num_free_wf_nodes)
    {
      // dual_u and dual_v are both unrestrained node
      for (int k = 0; k < 6; k++)
      {
        for (int l = 0; l < 6; l++)
        {
          double tmp;

          tmp = (*ptr_x)[i] * eK_[i](k, l);
          if (fabs(tmp) > SPT_EPS)
          {
            K_list.push_back(Eigen::Triplet<double>(dual_u * 6 + k, dual_u * 6 + l, tmp));
          }

          tmp = (*ptr_x)[i] * eK_[i](k + 6, l + 6);
          if (fabs(tmp) > SPT_EPS)
          {
            K_list.push_back(Eigen::Triplet<double>(dual_v * 6 + k, dual_v * 6 + l, tmp));
          }

          tmp = (*ptr_x)[i] * eK_[i](k, l + 6);
          if (fabs(tmp) > SPT_EPS)
          {
            K_list.push_back(Eigen::Triplet<double>(dual_u * 6 + k, dual_v * 6 + l, tmp));
          }

          tmp = (*ptr_x)[i] * eK_[i](k + 6, l);
          if (fabs(tmp) > SPT_EPS)
          {
            K_list.push_back(Eigen::Triplet<double>(dual_v * 6 + k, dual_u * 6 + l, tmp));
          }
        }
      }
    }
    else
    if (dual_u < num_free_wf_nodes)
    {
      // dual_u is free while dual_v is restrained
      for (int k = 0; k < 6; k++)
      {
        for (int l = 0; l < 6; l++)
        {
          double tmp = (*ptr_x)[i] * eK_[i](k, l);
          if (fabs(tmp) > SPT_EPS)
          {
            K_list.push_back(Eigen::Triplet<double>(dual_u * 6 + k, dual_u * 6 + l, tmp));
          }
        }
      }
    }
    else
    if (dual_v < num_free_wf_nodes)
    {
      for (int k = 0; k < 6; k++)
      {
        for (int l = 0; l < 6; l++)
        {
          double tmp = (*ptr_x)[i] * eK_[i](k + 6, l + 6);
          if (fabs(tmp) > SPT_EPS)
          {
            K_list.push_back(Eigen::Triplet<double>(dual_v * 6 + k, dual_v * 6 + l, tmp));
          }
        }
      }
    }
  }
  K_.setFromTriplets(K_list.begin(), K_list.end());

  if (verbose_)
  {
    create_k_.Stop();
  }
}


bool Stiffness::CalculateD(
    Eigen::VectorXd &D, Eigen::VectorXd *ptr_x, bool cond_num)
{
  D.resize(6 * num_free_wf_nodes);
  D.setZero();

  CreateGlobalK(ptr_x);
  CreateF(ptr_x);

  /* Parameter for StiffnessSolver */
  int	info;
//	IllCondDetector	stiff_inspector(K_);

  /* --- Check Stiffness Matrix Condition Number --- */
//	if (cond_num)
//	{
//		if (!CheckIllCondition(stiff_inspector))
//		{
//			return false;
//		}
//	}

  /* --- Solving Process --- */
  if (verbose_)
  {
    fprintf(stdout, "Stiffness : Linear Elastic Analysis ... Element Gravity Loads\n");
  }

  if (!stiff_solver_.SolveSystem(K_, D, F_, verbose_, info))
  {
    cout << "Stiffness Solver fail!\n" << endl;
    return false;
  }

  /* --- Equilibrium Error Check --- */
//	if (!CheckError(stiff_inspector, D))
//	{
//		return false;
//	}

  return true;
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


void Stiffness::PrintOutTimer()
{
  if (verbose_)
  {
    printf("***Stiffness timer result:\n");
    stiff_solver_.compute_k_.Print("ComputeK:");
    stiff_solver_.solve_d_.Print("SolveD:");
    create_fe_.Print("CreateFe:");
    create_f_.Print("CreateF:");
    create_ek_.Print("CreateElasticK:");
    create_k_.Print("CreateGlobalK:");
    check_ill_.Print("CheckIllCond:");
    check_error_.Print("CheckError:");
  }
  else
  {
    printf("***Stiffness detailed timing turned off.\n");
  }
}

} // namespace stiffness_checker
} // namespace conmech
