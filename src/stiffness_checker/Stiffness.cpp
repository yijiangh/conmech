#include <cmath>
#include <stiffness_checker/Stiffness.h>

namespace conmech
{
namespace stiffness_checker
{

Stiffness::Stiffness()
{
  terminal_output_ = false;
  ptr_dualgraph_ = NULL;
}

Stiffness::Stiffness(
    DualGraph *ptr_dualgraph, const StiffnessParm& parm,
    bool terminal_output)
{
  /* in use */
  ptr_dualgraph_ = ptr_dualgraph;
  parm_ = parm;

  nr_ = 0.0;
  shear_ = 1;

  Init();

  terminal_output_ = terminal_output;
}


Stiffness::~Stiffness()
{
}


void Stiffness::Init()
{
  CreateFe();
  CreateElasticK();

  num_free_wf_nodes = ptr_dualgraph_->SizeOfFreeFace();
}


void Stiffness::CreateFe()
{
  if (terminal_output_)
  {
    create_fe_.Start();
  }

  WireFrame *ptr_frame = ptr_dualgraph_->ptr_frame_;
  int num_wf_edges = ptr_dualgraph_->SizeOfVertList();
//  int num_free_wf_nodes = ptr_dualgraph_->SizeOfFaceList();

  CoordTrans coord_trans;
  // assume solid cross section
  double Ax = M_PI * parm_.radius_ * parm_.radius_;

  // gravitional acc unit: m/s^2
  double gx = 0;
  double gy = 0;
  double gz = parm_.g_;

  Fe_.resize(num_wf_edges);
  for (int i = 0; i < num_wf_edges; i++)
  {
    Fe_[i].setZero();

    WF_edge *ei = ptr_frame->GetEdge(ptr_dualgraph_->e_orig_id(i));
    WF_edge *ej = ei->ppair_;
    int dual_u = ptr_dualgraph_->v_dual_id(ej->pvert_->ID());
    int dual_v = ptr_dualgraph_->v_dual_id(ei->pvert_->ID());

    double t0, t1, t2, t3, t4, t5, t6, t7, t8;
    coord_trans.CreateTransMatrix(ej->pvert_->Position(), ei->pvert_->Position(),
        t0, t1, t2, t3, t4, t5, t6, t7, t8, 0.0);

    Eigen::VectorXd Fei(12);
    Fei.setZero();

    // half of an element's weight, first three entrees are force
    double L = ei->Length(); // mm
    double self_mass = parm_.density_ * Ax * L * 1e-3; // kg

    Fei[0] = Fei[6] = self_mass * gx/ 2.0;
    Fei[1] = Fei[7] = self_mass * gy / 2.0;
    Fei[2] = Fei[8] = self_mass * gz / 2.0;

    // bending caused by self weight, second three entrees are moment
    Fei[3] = self_mass * L / 12.0 * ((-t3*t7 + t4*t6)*gy + (-t3*t8 + t5*t6)*gz);
    Fei[4] = self_mass * L / 12.0 * ((-t4*t6 + t3*t7)*gx + (-t4*t8 + t5*t7)*gz);
    Fei[5] = self_mass * L / 12.0 * ((-t5*t6 + t3*t8)*gx + (-t5*t7 + t4*t8)*gy);

    Fei[9]  = self_mass * L / 12.0 * ((t3*t7 - t4*t6)*gy + (t3*t8 - t5*t6)*gz);
    Fei[10] = self_mass  * L / 12.0 * ((t4*t6 - t3*t7)*gx + (t4*t8 - t5*t7)*gz);
    Fei[11] = self_mass  * L / 12.0 * ((t5*t6 - t3*t8)*gx + (t5*t7 - t4*t8)*gy);

    Fe_[i] = Fei;
  }

  if (terminal_output_)
  {
    create_fe_.Stop();
  }
}


void Stiffness::CreateF(Eigen::VectorXd *ptr_x)
{
  if (terminal_output_)
  {
    create_f_.Start();
  }

  /* Run only after CreadFe is done! */
  WireFrame *ptr_frame = ptr_dualgraph_->ptr_frame_;
  int num_wf_edge = ptr_dualgraph_->SizeOfVertList();

  Eigen::VectorXd x;
  if (ptr_x == NULL)
  {
    ptr_x = &x;
    ptr_x->resize(num_wf_edge);
    ptr_x->setOnes();
  }

  // only considering the unsupported ones
  F_.resize(6 * num_free_wf_nodes);
  F_.setZero();

  for (int i = 0; i < num_wf_edge; i++)
  {
    WF_edge *ei = ptr_frame->GetEdge(ptr_dualgraph_->e_orig_id(i));
    WF_edge *ej = ei->ppair_;

    int dual_u = ptr_dualgraph_->v_dual_id(ej->pvert_->ID());
    int dual_v = ptr_dualgraph_->v_dual_id(ei->pvert_->ID());

    for (int j = 0; j < 6; j++)
    {
      // only unrestrained node is added into stiffness equation
      if (dual_u < num_free_wf_nodes)
      {
        F_[dual_u * 6 + j] += (*ptr_x)[i] * Fe_[i][j];
      }
      if (dual_v < num_free_wf_nodes)
      {
        F_[dual_v * 6 + j] += (*ptr_x)[i] * Fe_[i][j + 6];
      }
    }
  }

  if (terminal_output_)
  {
    create_f_.Stop();
  }
}


void Stiffness::CreateElasticK()
{
  if (terminal_output_)
  {
    create_ek_.Start();
  }

  WireFrame *ptr_frame   = ptr_dualgraph_->ptr_frame_;
  std::vector<WF_edge*> wf_edge_list = *ptr_frame->GetEdgeList();
  std::vector<WF_vert*> wf_vert_list = *ptr_frame->GetVertList();

  int num_wf_edges = ptr_dualgraph_->SizeOfVertList();
//  int num_free_wf_nodes = ptr_dualgraph_->SizeOfFaceList();

  /* ref to Matrix Analysis of Strutures Aslam Kassimali, Section 8.2 Table 8.1*/
  double Ax = F_PI * std::pow(parm_.radius_,2);
  double Asy = Ax * (6 + 12 * parm_.poisson_ratio_ + 6 * std::pow(parm_.poisson_ratio_,2))
      / (7 + 12 * parm_.poisson_ratio_ + 4 * std::pow(parm_.poisson_ratio_,2));
  double Asz = Asy;

  /* torsion constant */
  double Jxx = 0.5 * F_PI * std::pow(parm_.radius_, 4);

  /* shear deformation constant */
  double Ksy = 0;
  double Ksz = 0;

  /* area moment of inertia (bending about local y,z-axis)
  * https://en.wikipedia.org/wiki/Bending (beam deflection equation)
  * https://en.wikipedia.org/wiki/List_of_area_moments_of_inertia (numerical value)
  * note this is slender rod of length L and Mass M, spinning around end
  */
  double Iyy = F_PI * std::pow(parm_.radius_,4) / 4;
  double Izz = Iyy;

  double t0, t1, t2, t3, t4, t5, t6, t7, t8;     /* coord transf matrix entries */

  double E = parm_.youngs_modulus_;
  double G = parm_.shear_modulus_;

  eK_.resize(num_wf_edges);
  for (int i = 0; i < num_wf_edges; i++)
  {
    eK_[i].setZero();

    //WF_edge *ei = ptr_frame->GetNeighborEdge(ptr_dualgraph_->v_orig_id(i));
    WF_edge *ei = wf_edge_list[ptr_dualgraph_->e_orig_id(i)];

    int u = ei->ppair_->pvert_->ID();
    int v = ei->pvert_->ID();
    double L = ei->Length();
    double Le = L - 2 * nr_;

    Eigen::MatrixXd eKuv(12, 12);
    eKuv.setZero();

    trimesh::point node_u = wf_vert_list[u]->Position();
    trimesh::point node_v = wf_vert_list[v]->Position();

    transf_.CreateTransMatrix(node_u, node_v, t0, t1, t2, t3, t4, t5, t6, t7, t8, 0);

    if (1 == shear_)
    {
      /* for circular cross-sections, the shape factor for shear(fs) = 1.2 (ref. Aslam's book) */
      double fs = 1.2;
      Ksy = 12. * E * Izz * fs / (G * Asy * Le * Le);
      Ksz = 12. * E * Iyy * fs / (G * Asz * Le * Le);
    }
    else
    {
      Ksy = 0;
      Ksz = 0;
    }

    int n1 = ptr_dualgraph_->v_dual_id(u);
    int n2 = ptr_dualgraph_->v_dual_id(v);

    eKuv(0,0) = eKuv(6,6) = E * Ax / Le;
    eKuv(1,1) = eKuv(7,7) = 12. * E * Izz / (Le * Le * Le * (1. + Ksy));
    eKuv(2,2) = eKuv(8,8) = 12. * E * Iyy / (Le * Le * Le * (1. + Ksz));
    eKuv(3,3) = eKuv(9,9) =  G * Jxx / Le;
    eKuv(4,4) = eKuv(10,10) = (4. + Ksz) * E * Iyy / (Le * (1. + Ksz));
    eKuv(5,5) = eKuv(11,11) = (4. + Ksy) * E * Izz / (Le * (1. + Ksy));

    eKuv(4,2) = eKuv(2,4) = -6. * E * Iyy / (Le * Le * (1. + Ksz));
    eKuv(5,1) = eKuv(1,5) = 6. * E * Izz / (Le * Le * (1. + Ksy));
    eKuv(6,0) = eKuv(0,6) = - eKuv(0,0);

    eKuv(11,7) = eKuv(7,11) = eKuv(7,5) = eKuv(5,7) = -eKuv(5,1);
    eKuv(10,8) = eKuv(8,10) = eKuv(8,4) = eKuv(4,8) = -eKuv(4,2);
    eKuv(9,3) = eKuv(3,9)   = - eKuv(3,3);
    eKuv(10, 2) = eKuv(2, 10) = eKuv(4, 2);
    eKuv(11, 1) = eKuv(1, 11) = eKuv(5, 1);

    eKuv(7,1) = eKuv(1,7) = -eKuv(1,1);
    eKuv(8,2) = eKuv(2,8) = -eKuv(2,2);
    eKuv(10,4) = eKuv(4,10) = (2. - Ksz) * E * Iyy / (Le * (1. + Ksz));
    eKuv(11,5) = eKuv(5,11) = (2. - Ksy) * E * Izz / (Le * (1. + Ksy));

    transf_.TransLocToGlob(t0, t1, t2, t3, t4, t5, t6, t7, t8, eKuv, 0, 0);

    // symmetry check
    for (int k = 0; k < 12; k++)
    {
      for (int l = k + 1; l < 12; l++)
      {
        if (eKuv(k, l) != eKuv(l, k))
        {
          if (fabs(eKuv(k,l) / eKuv(l,k) - 1.0) > SPT_EPS
              &&
                  (fabs(eKuv(k, l) / eKuv(k, k)) > SPT_EPS || fabs(eKuv(l, k) / eKuv(k, k)) > SPT_EPS)
              )
          {
            fprintf(stderr, "elastic_K: element stiffness matrix not symetric ...\n");
            fprintf(stderr, " ... k[%d][%d] = %15.6e \n", k, l, eKuv(k,l));
            fprintf(stderr, " ... k[%d][%d] = %15.6e   ", l, k, eKuv(l,k));
            fprintf(stderr, " ... relative error = %e \n", fabs(eKuv(k,l) / eKuv(l,k) - 1.0));
          }

          eKuv(k, l) = eKuv(l, k) = (eKuv(k, l) + eKuv(l, k)) / 2;
        }
      }
    }

    eK_[i] = eKuv;
  }

  if (terminal_output_)
  {
    create_ek_.Stop();
  }
}


void Stiffness::CreateGlobalK(Eigen::VectorXd *ptr_x)
{
  if (terminal_output_)
  {
    create_k_.Start();
  }

  WireFrame *ptr_frame = ptr_dualgraph_->ptr_frame_;
  int num_wf_edges = ptr_dualgraph_->SizeOfVertList();
  int num_free_wf_nodes = ptr_dualgraph_->SizeOfFaceList();

  Eigen::VectorXd x;
  if (ptr_x == NULL)
  {
    ptr_x = &x;
    ptr_x->resize(num_wf_edges);
    ptr_x->setOnes();
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

  if (terminal_output_)
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
  if (terminal_output_)
  {
    fprintf(stdout, "Stiffness : Linear Elastic Analysis ... Element Gravity Loads\n");
  }

  if (!stiff_solver_.SolveSystem(K_, D, F_, terminal_output_, info))
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
//	if (terminal_output_)
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
//		if (terminal_output_)
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
//	if (terminal_output_)
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
//	if (terminal_output_)
//	{
//		check_error_.Start();
//	}
//
//	bool bSuccess = true;
//	double error = stiff_inspector.EquilibriumError(K_, D, F_);
//	if (terminal_output_)
//	{
//		printf("Root Mean Square (RMS) equilibrium error = %9.3e\n", error);
//	}
//	if (error < STIFF_TOL)
//	{
//		if (terminal_output_)
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
//	if (terminal_output_)
//	{
//		check_error_.Stop();
//	}
//
//	return bSuccess;
//}


void Stiffness::PrintOutTimer()
{
  if (terminal_output_)
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
