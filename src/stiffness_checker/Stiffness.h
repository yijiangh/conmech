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
  ~Stiffness() {}

 public:
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
   * @param exist_element_ids
   * (n_desire x 1) int std::vector
   * n_desire is the number of nodes that you wants to prescribe
   *    exist_element_ids[i] = node_id,
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
      const std::vector<int>& exist_element_ids,
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
   * create a list of local element stiffness matrix,
   * and store it to the class var local_K_list_.
   */
  void createLocalStiffnessMatrixList();

  /**
   * assemble global stiffness matrix from the saved list of
   * stiffness matrix, given the existing elements' ids.
   * @param exist_element_ids
   * (n_desire x 1) int std::vector
   *    exist_element_ids[i] = element's id
   * @param[out] K_sp Sparse global stiffness matrix
   * (TODO: dimension?)
   *
   * @param[out] id_map
   */
  void createGlobalStiffnessMatrix(const std::vector<int>& exist_element_ids,
                                   Eigen::Sparse<double>& K_assembled,
                                   Eigen::VectorXi& id_map);

 protected:
  Frame frame_;
  StiffnessParm material_parm_;
  StiffnessSolver stiff_solver_;

  /**
   * a list of element stiffness matrix in global frame
   * (N_all_element x (12x12)) list
   * Ne * 1 std::vector, Ne is number of elements
   * each entry is a 12*12 local stiffness matrix K
   */
  std::vector<Eigen::MatrixXd> element_K_list_;

  /**
   * external nodal load P
   * (N_ext_loaded x 7) double matrix
   *    ext_load_P[i] = [node_id, (6x1)load vector]
   */
  Eigen::MatrixXd ext_load_P_;

  /**
   * self-weight nodal load P
   * (N_all_node x 6) double matrix
   * @note the size of this matrix is fixed to include
   * all the nodes in the structure. In solving process,
   * a sub-vector can be carved out for specific partial structure.
   */
  Eigen::MatrixXd self_weight_load_P_;

  /**
   * (N_fixities_node x 7) int matrix
   *    fixities[i] = [n_id, is_Tx, is_Ty, is_Tz, is_Rx, is_Ry, is_Rz]
   */
  Eigen::MatrixXi fixities_;

  Timer create_k_;
  Timer check_ill_;
  Timer check_error_;

  bool verbose_;
};

} // namespace stiffness_checker
} // namespace conmech