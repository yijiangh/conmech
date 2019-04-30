#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>
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
//  Stiffness(Frame &frame, bool verbose = false, std::string model_type = "frame");

  Stiffness(const std::string& json_file_path, bool verbose = false, const std::string& model_type = "frame", bool output_json = false);

  ~Stiffness() {}

public:
  /**
   * set maximal nodal translational and rotational tolerance
   * for stiffness checking criteria.
   */
  void setNodalDisplacementTolerance(double transl_tol, double rot_tol);

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
   */
  void setLoad(const Eigen::MatrixXd &nodal_forces);

  void setSelfWeightNodalLoad(bool include_sw) { include_self_weight_load_ = include_sw; }

  void setOutputJsonPath(const std::string& file_path, const std::string& file_name)
  {
    output_json_file_name_ = file_name;
    output_json_file_path_ = file_path;
  }

  void setOutputJson(const bool output_json) { write_result_ = output_json; }

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
    const std::vector<int> &exist_element_ids,
    Eigen::MatrixXd &node_displ,
    Eigen::MatrixXd &fixities_reaction,
    Eigen::MatrixXd &element_reation,
    const bool &cond_num = true);

  bool solve(
    const std::vector<int> &exist_element_ids,
    const bool &cond_num = true);

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
    Eigen::MatrixXd &node_displ,
    Eigen::MatrixXd &fixities_reaction,
    Eigen::MatrixXd &element_reation,
    const bool &cond_num = true);

  bool solve(const bool &cond_num = true);

  bool hasStoredResults() const { return has_stored_deformation_; }

  // TODO: pass by const reference or inmutatable
  bool getSolvedResults (Eigen::MatrixXd &node_displ,
                        Eigen::MatrixXd &fixities_reaction,
                        Eigen::MatrixXd &element_reaction,
                        bool &pass_criteria);

  bool getMaxNodalDeformation(double &max_trans, double &max_rot);

  double getTransTol() const { return transl_tol_; }
  double getRotTol() const { return rot_tol_; }

  Eigen::MatrixXd getOriginalShape(const int& disc=1, const bool& draw_full_shape=true);

  /** compute the cubic interpolated deformed shape
   * @param[in] deformation exaggeration ratio
   * @param[out] (nElement x 6) Eigen::MatrixXd DB
   *  DB[i,:] = (x_i, y_i, z_i, x_i, y_i, z_i)
   * @return success
   */
  Eigen::MatrixXd getDeformedShape(const double& exagg_ratio, const int& disc);

  bool computeCubicDeformedBeam(const Eigen::VectorXd& end_u, const Eigen::VectorXd& end_v,
                                const Eigen::VectorXd& d_end_u, const Eigen::VectorXd& d_end_v,
                                const double& exagg, const int& disc,
                                Eigen::MatrixXd& BeamPolygon);

  /* Check condition number */
//  bool CheckIllCondition(IllCondDetector &stiff_inspector);
//  bool CheckError(IllCondDetector &stiff_inspector, Eigen::VectorXd &D);

  /**
   * print out timing result in console.
   */
  void printOutTimer();

protected:

  bool init();

  virtual bool checkStiffnessCriteria(const Eigen::MatrixXd &node_displ,
                                      const Eigen::MatrixXd &fixities_reaction,
                                      const Eigen::MatrixXd &element_reation);

private:
  /**
   * create external nodal load to class var ext_load_P
   * @param nodal_forces
   * (n_loaded_nodes x (1+node_dof)) matrix,
   * n_loaded_nodes is the number of prescribed loaded nodes
   *    nodal_forces[i, :] = [node_id, Fx, Fy, Fz, Mx, My, Mz] (3D) or [node_id, Fx, Fy, Mz] (2D)
   *    described in global (world) frame. 0 <= node_id < N_vert.
   * @param[out] ext_load (dof x 1)
   */
  void createExternalNodalLoad(const Eigen::MatrixXd &nodal_forces, Eigen::VectorXd &ext_load);

  /**
   * create a list of element stiffness matrix (in global frame),
   * and store it to the class var local_K_list_.
   */
  void createElementStiffnessMatrixList();

  void createElementSelfWeightNodalLoad();

  /**
   * assemble global stiffness matrix from all elements,
   * i.e. (dof x dof) K_assembled, dof = n_Node*6
   */
  void createCompleteGlobalStiffnessMatrix(const std::vector<int> &exist_e_ids);

  /**
   * convert self-weight load (between nodal points) to nodal
   * loads and save it to the class var self_weight_load_P.
   */
  void createSelfWeightNodalLoad(const std::vector<int>& exist_e_ids, Eigen::VectorXd& self_weight_load);

protected:
  Frame frame_;
  StiffnessParm material_parm_;
  StiffnessSolver stiff_solver_;

  Timer create_k_;
  Timer check_ill_;

  double transl_tol_;
  double rot_tol_;

  bool verbose_;
  bool write_result_;

  std::string output_json_file_path_;
  std::string output_json_file_name_;

  /**
   * dimension of the model, 2 or 3
   */
  int dim_;

  /**
   * node degrees of freedom
   */
  int node_dof_;

  /**
   * node full degree of freedom (aka. frame)
   */
  int full_node_dof_;

  /**
   * elemental reaction dof indices among (0, ..., node_dof-1, node_dof, ..., 2*node_dof-1)
   * This is used for retrieving axial translational dof out the full dof that includes rotation.
   */
  Eigen::VectorXi e_react_dof_id_;

  /**
   * free dof indices among (0, ..., node_dof-1, node_dof, ..., 2*node_dof-1)
   * This is in particular used for retrieving translational dof out of the full dof that includes rotation.
   */
  Eigen::VectorXi xyz_dof_id_;

  /**
   * model type: truss or frame
   */
  std::string model_type_;

private:
  /**
   * a (N_element x (2*node_dof)) map
   * element id -> dof id map
   */
  Eigen::MatrixXi id_map_;

  /**
   * a ((num_fix_node) x node_dof_) eigen int matrix
   */
  Eigen::MatrixXi fixities_table_;

  /**
   * a list of element stiffness matrix in global frame
   * (N_all_element x (12x12)) list
   * Ne * 1 std::vector, Ne is number of elements
   * each entry is a 12*12 local stiffness matrix K
   * Note: length unit: mm, force unit: N
   */
  std::vector<Eigen::MatrixXd> element_K_list_;

  std::vector<Eigen::MatrixXd> rot_m_list_;

  std::vector<Eigen::VectorXd> element_gravity_nload_list_;

  /**
   * Global assembled stiffness matrix for the entire structure
   * the matrix will be sliced to solve for partial-structure
   * by the getSlicesGlobalStiffnessMatrix method.
   * Note: length unit: mm, force unit: N
   */
  Eigen::SparseMatrix<double> K_assembled_full_;

  /**
   * external nodal load P
   * (dof x 1) double vector
   * Note: force unit: kN
   */
  Eigen::VectorXd nodal_load_P_;

  /**
   * (N_fixities_node x 7) int matrix
   *    fixities[i] = [n_id, is_Tx, is_Ty, is_Tz, is_Rx, is_Ry, is_Rz]
   */
  Eigen::MatrixXi fixities_;

  /**
   * stored computation results
   */
  std::vector<int> stored_existing_ids_;

  Eigen::MatrixXd stored_nodal_deformation_;
  Eigen::MatrixXd stored_element_reaction_;
  Eigen::MatrixXd stored_fixities_reaction_;

  /**
   * boolean flag for if the model is inited (1) or not (0).
   */
  bool is_init_;

  bool include_self_weight_load_;

  bool has_stored_deformation_;
};

} // namespace stiffness_checker
} // namespace conmech
