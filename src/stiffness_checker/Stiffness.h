#pragma once

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include "stiffness_checker/Frame.h"
#include "stiffness_checker/StiffnessParm.h"
#include "stiffness_checker/StiffnessSolver.h"
//#include "stiffness_checker/IllCondDetector.hpp"

namespace conmech
{
namespace stiffness_checker
{

class Stiffness
{
 public:
  Stiffness(Frame& frame, bool verbose=false);
  Stiffness(const std::string& json_file_path, bool verbose=false);

  ~Stiffness() {}

 public:

  bool init();

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
  // TODO: not implemented yet
  virtual bool setRigidFixities(const Eigen::MatrixXi& fixities) {}

  /**
   * Compute nodal displacement given existing node's indices.
   * @param exist_element_ids
   * (n_desire x 1) int std::vector
   * n_desire is the number of nodes that you wants to prescribe
   *    exist_element_ids[i] = node_id,
   *    if the corresponding node exists in the (partial-)structure
   *    that you want to evaluate.
   * @param[out] node_displ
   * (n_exist x 3) double matrix
   *    node_displ[i] = [node_id, dof_id, disp.]
   *    where dof_id in [0,5]
   * @param[out] fixities_reaction
   * (n_fix x 3) double matrix
   *    fixities_reaction[i] = [node_id, dof_id, R]
   * @param[out] element_reation
   * (n_elements x 7) double matrix
   *    element_reaction[i] = [element_id, F_x, F_y, F_z, M_x, M_y, M_z]
   *    element's internal force described in **local frame**.
   * @param cond_num flag for including condition number checking, default true
   * @return boolean success flag
   */
  bool solve(
      const std::vector<int>& exist_element_ids,
      Eigen::MatrixXd& node_displ,
      Eigen::MatrixXd& fixities_reaction,
      Eigen::MatrixXd& element_reation,
      const bool& cond_num = true);

  bool solve(
      const std::vector<int>& exist_element_ids,
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

  bool solve(
      const bool& cond_num = true);

  /* Check condition number */
//  bool CheckIllCondition(IllCondDetector &stiff_inspector);
//  bool CheckError(IllCondDetector &stiff_inspector, Eigen::VectorXd &D);

  /**
   * print out timing result in console.
   */
  void printOutTimer();

 private:
  /**
   * create external nodal load to class var ext_load_P
   * @param nodal_forces
   * (n_desire x 7) matrix,
   * n_desire is the number of nodes that you wants to prescribe
   *    nodal_forces[i, :] = [node_id, Fx, Fy, Fz, Mx, My, Mz]
   *    described in global (world) frame. 0 <= node_id < N_vert.
   */
  void createExternalNodalLoad(const Eigen::MatrixXd& nodal_forces, Eigen::VectorXd& ext_load);

  /**
   * convert self-weight load (between nodal points) to nodal
   * loads and save it to the class var self_weight_load_P.
   */
  void createSelfWeightNodalLoad(Eigen::VectorXd& self_weight_load);

  /**
   * set all grounded nodes in the input frame
   * to the same type of fixities
   * @param Tx flag for fixities restrain
   * @param Ty
   * @param Tz
   * @param Rx
   * @param Ry
   * @param Rz
   */
  void setGroundedNodesUniformFixities(bool Tx, bool Ty, bool Tz, bool Rx, bool Ry, bool Rz);

  /**
   * create a list of element stiffness matrix (in global frame),
   * and store it to the class var local_K_list_.
   */
  void createElementStiffnessMatrixList();

  /**
   * assemble global stiffness matrix from all elements,
   * i.e. (dof x dof) K_assembled, dof = n_Node*6
   */
  void createCompleteGlobalStiffnessMatrix();

 protected:
  Frame frame_;
  StiffnessParm material_parm_;
  StiffnessSolver stiff_solver_;

  Timer create_k_;
  Timer check_ill_;

  bool verbose_;

 private:
  /**
   * a list of element stiffness matrix in global frame
   * (N_all_element x (12x12)) list
   * Ne * 1 std::vector, Ne is number of elements
   * each entry is a 12*12 local stiffness matrix K
   */
  std::vector<Eigen::MatrixXd> element_K_list_;

  /**
   * Global assembled stiffness matrix for the entire structure
   * the matrix will be sliced to solve for partial-structure
   * by the getSlicesGlobalStiffnessMatrix method.
   */
  Eigen::MatrixXd K_assembled_full_;

  /**
   * external nodal load P
   * (dof) double vector
   */
  Eigen::VectorXd nodal_load_P_;

  /**
   * (N_fixities_node x 7) int matrix
   *    fixities[i] = [n_id, is_Tx, is_Ty, is_Tz, is_Rx, is_Ry, is_Rz]
   */
  Eigen::MatrixXi fixities_;

  bool is_init_;
};

} // namespace stiffness_checker
} // namespace conmech