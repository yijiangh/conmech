#include "stiffness_checker/Stiffness.h"
#include "stiffness_checker/Util.h"

namespace conmech
{
namespace stiffness_checker
{
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

  for(std::size_t i=0; i < draw_ids.size(); i++)
  {
    int end_u_id = Elements_(draw_ids[i], 0);
    int end_v_id = Elements_(draw_ids[i], 1);
    Eigen::VectorXd end_u, end_v;
    getNodePoints(Vertices_, end_u_id, end_v_id, end_u, end_v);

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

  for(std::size_t i=0; i<stored_existing_ids_.size(); i++)
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

    Eigen::VectorXd end_u, end_v;
    getNodePoints(Vertices_, end_u_id, end_v_id, end_u, end_v);

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
} // ns stiffness_checker
} // ns conmech