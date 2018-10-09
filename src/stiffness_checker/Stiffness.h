#pragma once

#include "stiffness_checker/Frame.h"
#include "stiffness_checker/StiffnessParm.h"
#include "stiffness_checker/Util.h"
#include "stiffness_checker/GCommon.h"
#include "stiffness_checker/StiffnessSolver.h"
//#include "stiffness_checker/IllCondDetector.hpp"

namespace conmech
{
namespace stiffness_checker
{

class Stiffness
{
 public:
  Stiffness();
  ~Stiffness();

 public:
  void init();

  /**
   * Construct external load (force, moment) on nodes.
   * @param nodal_forces
   * (n_desire x 7) matrix,
   * n_desire is the number of nodes that you wants to prescribe
   *    nodal_forces[i, :] = [node_id, Fx, Fy, Fz, Mx, My, Mz]
   *    described in global (world) frame. 0 <= node_id < N_vert.
   * @param include_self_weight
   *    boolean flag for include element's self weight.
   *    default to be false.
   * @return boolean success flag
   */
  bool setNodalLoad(const Eigen::MatrixXd& nodal_forces,
                    const bool& include_self_weight = false);

  /**
   * Construct infinitely rigid fixities on nodes.
   * @param fixities
   * (n_desire x 7) matrix,
   * n_desire is the number of nodes that you wants to prescribe
   *    fixities[i, :] = [node_id, is_Fx, is_Fy, is_Fz, is_Mx, is_My, is_Mz]
   *    where is_X is a boolean flag.
   * @return boolean success flag
   */
  bool setRigidFixities(const Eigen::Matrixi& fixities);

  /**
   * Compute nodal displacement given existing node's indices.
   * @param exist_node_ids
   * (n_desire x 1) int vector
   * n_desire is the number of nodes that you wants to prescribe
   *    exist_node_ids[i] = node_id,
   *    if the corresponding node exists in the (partial-)structure
   *    that you want to evaluate.
   * @param[out] node_displ
   * (n_exist x 7) double matrix
   *    node_displ[i] = [node_id, d_x, d_y, d_z, theta_x, theta_y, theta_z]
   *    where d_() is the translation displacement,
   *    theta_() is the rotational displ.
   *    See [MSA McGuire et al.] P75.
   * @param[out] fixities_reaction
   * (n_fix x 7) double matrix
   *    fixities_reaction[i] = [node_id, F_x, F_y, F_z, M_x, M_y, M_z]
   * @param[out] element_reation
   * (n_elements x 7) double matrix
   *    element_reaction[i] = [node_id, F_x, F_y, F_z, M_x, M_y, M_z]
   *    element's internal force described in global frame.
   * @param cond_num flag for including condition number checking, default true
   * @return boolean success flag
   */
  bool solve(
      const Eigen::VectorXi& exist_node_ids,
      Eigen::MatrixXd& node_displ,
      Eigen::MatrixXd& fixities_reaction,
      Eigen::MatrixXd& element_reation,
      const bool& cond_num = true);

  /**
   * Compute nodal displacement for the entire structure (assuming all elements exist).
   * @param[out] node_displ
   * (n_exist x 7) double matrix
   *    node_displ[i] = [node_id, d_x, d_y, d_z, theta_x, theta_y, theta_z]
   *    where d_() is the translation displacement,
   *    theta_() is the rotational displ.
   *    See [MSA McGuire et al.] P75.
   * @param[out] fixities_reaction
   * (n_fix x 7) double matrix
   *    fixities_reaction[i] = [node_id, F_x, F_y, F_z, M_x, M_y, M_z]
   * @param[out] element_reation
   * (n_elements x 7) double matrix
   *    element_reaction[i] = [node_id, F_x, F_y, F_z, M_x, M_y, M_z]
   *    element's internal force described in global frame.
   * @param cond_num flag for including condition number checking, default true
   * @return boolean success flag
   */
  bool solve(
      Eigen::MatrixXd& node_displ,
      Eigen::MatrixXd& fixities_reaction,
      Eigen::MatrixXd& element_reation,
      const bool& cond_num = true);

  /* Check condition number */
//  bool CheckIllCondition(IllCondDetector &stiff_inspector);
//  bool CheckError(IllCondDetector &stiff_inspector, Eigen::VectorXd &D);

  /**
   * print out timing result in console.
   */
  void PrintOutTimer();

 private:
  /**
   * create external nodal load to class var ext_load_P
   * @param nodal_forces
   * (n_desire x 7) matrix,
   * n_desire is the number of nodes that you wants to prescribe
   *    nodal_forces[i, :] = [node_id, Fx, Fy, Fz, Mx, My, Mz]
   *    described in global (world) frame. 0 <= node_id < N_vert.
   */
  void createExternalNodalLoad(const Eigen::MatrixXd& nodal_forces);

  /**
   * convert self-weight load (between nodal points) to nodal
   * loads and save it to the class var self_weight_load_P.
   */
  void createSelfWeightNodalLoad();

  void createLocalStiffnessMatrixList();

 protected:
  Frame frame_;
  StiffnessParm material_parm_;
  StiffnessSolver stiff_solver_;

  /**
   * a list of local stiffness matrix (TODO: in global or local frame?)
   * (N_all_element x (6x6)) list
   * Ne * 1 std::vector, Ne is number of elements
   * each entry is a 6*6 local stiffness matrix K
   */
  std::vector<Eigen::MatrixXd> local_K_list_;

  /**
   * external nodal load P
   * (N_all_node*6 x 1) double vector
   * @note the size of this vector is fixed to include
   * all the nodes in the structure. In solving process,
   * a sub-vector can be carved out for specific partial structure.
   */
  Eigen::VectorXd ext_load_P_;

  /**
   * self-weight nodal load P
   * (N_all_node*6 x 1) double vector
   * @note the size of this vector is fixed to include
   * all the nodes in the structure. In solving process,
   * a sub-vector can be carved out for specific partial structure.
   */
  Eigen::VectorXd self_weight_load_P_;

  Timer create_fe_;
  Timer create_f_;
  Timer create_ek_;
  Timer create_k_;
  Timer check_ill_;
  Timer check_error_;

  bool verbose_;
};

} // namespace stiffness_checker
} // namespace conmech