#pragma once

#include "stiffness_checker/Material.h"
#include "stiffness_checker/StiffnessSolver.h"
#include <nlohmann/json.hpp>
#include <Eigen/Dense>
#include <Eigen/Sparse>
// TODO: check if needed: https://eigen.tuxfamily.org/dox/group__TopicStlContainers.html
// #include<Eigen/StdVector>

namespace conmech
{
namespace stiffness_checker
{
class Stiffness
{
public:
  /**
   * @brief Factory function - returned by value:
   * 
   * @param V 
   * @param E 
   * @param Fixities 
   * @param materials 
   * @param verbose 
   * @param model_type 
   * @param output_json 
   * @return Stiffness 
   */
  static Stiffness create(const Eigen::MatrixXd& V, const Eigen::MatrixXi& E, const Eigen::MatrixXi& Fixities,
                          const std::vector<conmech::material::Material>& materials,
                          const bool& verbose = false, const std::string& model_type = "frame", const bool& output_json = false);

  Stiffness(const std::string& file_path,
            const bool& verbose = false, const std::string& model_type = "frame", const bool& output_json = false);

  /**
   * @brief Construct a new Stiffness object
   * 
   * @param V : frame vertices coordinate list, n by 3 (2), where n is the number of vertices (#V)
   * @param E : element index list, m by 2, where m is the number of elements (#E)
   * @param BC : boundary condition fixities list, each row, vertex id, dof fixity indication, fn by 7 (4), TODO: dependent on frame or truss
   * @param materials : list of materials, parsed from json
   * @param verbose 
   * @param model_type : "frame" or "truss", "truss" unsupported now
   * @param output_json : if write analysis result to a json file
   */
  Stiffness(const Eigen::MatrixXd& V, const Eigen::MatrixXi& E, const Eigen::MatrixXi& Fixities,
            const std::vector<conmech::material::Material>& materials,
            const bool& verbose = false, const std::string& model_type = "frame", const bool& output_json = false);

  Stiffness(const Eigen::MatrixXd& V, const Eigen::MatrixXi& E, const Eigen::MatrixXi& Fixities,
            const std::vector<nlohmann::json>& material_jsons,
            const bool& verbose = false, const std::string& model_type = "frame", const bool& output_json = false);

  ~Stiffness();

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
   */
  void setLoad(const Eigen::MatrixXd &nodal_forces);

  /**
   * @brief Set the Uniformly Distributed Load
   * 
   * @param element_load_density 
   * (n_desired_elements x 4) matrix
   *  element_load_density[i, :] = [e_id, wx, wy, wz]
   *  described in global frame, wx, wy, wz are load densities in kN/m
   */
  void setUniformlyDistributedLoad(const Eigen::MatrixXd &element_load_density);

  void setSelfWeightNodalLoad(bool include_sw) { include_self_weight_load_ = include_sw;}
  void setGravityDirection(const Eigen::VectorXd& gravity_direction);

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

  // in global coordinate
  bool getElementStiffnessMatrices(std::vector<Eigen::MatrixXd> &element_stiffness_mats);
  bool getElementLocal2GlobalRotationMatrices(std::vector<Eigen::MatrixXd> &e_L2G_rot_mats);

  // for now, these load vectors contain all the nodes in the full structure
  // no index column
  bool getExternalNodalLoad(Eigen::VectorXd& ext_point_load) { ext_point_load = nodal_load_P_; return true; }
  bool getUniformlyDistributedLumpedLoad(const std::vector<int>& exist_e_ids, Eigen::VectorXd& lumped_load);
  bool getSelfWeightNodalLoad(const std::vector<int>& exist_e_ids, Eigen::VectorXd& self_weight_load);

  bool isIncludeSelfWeightLoad() const { return include_self_weight_load_; }

  bool hasStoredResults() const { return has_stored_deformation_; }

  // TODO: pass by const reference or inmutatable
  bool getSolvedResults (Eigen::MatrixXd &node_displ,
                        Eigen::MatrixXd &fixities_reaction,
                        Eigen::MatrixXd &element_reaction,
                        bool &pass_criteria);

  bool getMaxNodalDeformation(double &max_trans, double &max_rot,
                              int &max_trans_vid, int &max_rot_id);

  bool getSolvedCompliance(double &compliance);
  bool getSolvedAltCompliance(double &compliance);

  double getTransTol() const { return transl_tol_; }
  double getRotTol() const { return rot_tol_; }

  // TODO: avoid copy?
  Eigen::MatrixXi getElement2DofIdMap() const { return id_map_; }
  Eigen::MatrixXi getNode2DofIdMap() const { return v_id_map_; }

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

  /**
   * @brief total number of elements
   * 
   * @return int 
   */
  int nE() const;

  /**
   * @brief total number of vertices
   * 
   * @return int 
   */
  int nV() const;

  /**
   * @brief dimension of the structure
   * 
   * @return int const 
   */
  int dim() const;

  /**
   * @brief number of fixed vertices
   * 
   * @return int 
   */
  int nFixV() const;

  void computeLumpedUniformlyDistributedLoad(const Eigen::Vector3d &w_G, const Eigen::Matrix3d &R_LG, const double &Le, 
    Eigen::VectorXd &eq_nodal_load);

  /* Check condition number */
//  bool CheckIllCondition(IllCondDetector &stiff_inspector);
//  bool CheckError(IllCondDetector &stiff_inspector, Eigen::VectorXd &D);

  /**
   * print out timing result in console.
   */
  void printOutTimer();

protected:

  bool init(const Eigen::MatrixXd& V, const Eigen::MatrixXi& E, const Eigen::MatrixXi& Fixities, 
            const std::vector<conmech::material::Material>& materials,
            const bool& verbose, const std::string& model_type, const bool& output_json);

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
   * @brief Create a Uniformly Distributed Load
   * 
   * NOTE: this will erase previously assigned elemental loads!
   * 
   * @param element_load_density 
   * (n_loaded_elements x (1+3)) matrix,
   * n_loaded_elements is the number of prescribed loaded elements
   *    nodal_forces[i, :] = [element_id, wx, wy, wz] (3D)
   *    described in global (world) frame. 0 <= element_id < N_vert.
   * @param[out] ext_load lumped equivalent nodal load (full dofs)
   */
  void precomputeElementUniformlyDistributedLumpedLoad(const Eigen::MatrixXd &element_load_density);

  /**
   * @brief Create lumped nodal loads for external element-wise loads, this is called
   * when ``solve`` is called.
   * 
   * @param exist_e_ids : existing elements' indices
   * @param[out] ext_load lumped equivalent nodal load (full dofs)
   */
  void createUniformlyDistributedLumpedLoad(const std::vector<int>& exist_e_ids, Eigen::VectorXd &ext_load);


  /**
   * @brief Precompute lumped nodal loads for element self weight load.
   * 
   */
  void precomputeElementSelfWeightLumpedLoad();

  /**
   * convert self-weight load (between nodal points) to nodal
   * loads and save it to the class var self_weight_load_P.
   * 
   * This function is called whenever ``solve`` is called.
   */
  void createSelfWeightLumpedLoad(const std::vector<int>& exist_e_ids, Eigen::VectorXd& self_weight_load);

  /**
   * create a list of element stiffness matrix (in global frame),
   * and store it to the class var local_K_list_.
   */
  void precomputeElementStiffnessMatrixList();

  /**
   * assemble global stiffness matrix from all elements,
   * i.e. (dof x dof) K_assembled, dof = n_Node*6
   * This is called at solve-time.
   */
  void createCompleteGlobalStiffnessMatrix(const std::vector<int> &exist_e_ids);

protected:
  Eigen::MatrixXd Vertices_;
  Eigen::MatrixXi Elements_;
  /**
   * (N_fixities_node x 7) int matrix
   *    fixities[i] = [n_id, is_Tx, is_Ty, is_Tz, is_Rx, is_Ry, is_Rz]
   */
  Eigen::MatrixXi Fixities_;

  std::vector<conmech::material::Material> materials_;
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
   * a (N_node x (2*node_dof)) map
   * node id -> dof id map
   */
  Eigen::MatrixXi v_id_map_;

  // /**
  //  * a ((num_fix_node) x node_dof_) eigen int matrix
  //  */
  // Eigen::MatrixXi fixities_table_;

  /**
   * a list of element stiffness matrix in global frame
   * (N_all_element x (12x12)) list
   * Ne * 1 std::vector, Ne is number of elements
   * each entry is a 12*12 local stiffness matrix K
   * Note: length unit: mm, force unit: N
   */
  std::vector<Eigen::MatrixXd> element_K_list_;

  /**
   * @brief a list of 3x3 global-to-local coordinate transformation matrix
   * (N_all_element x (3x3)) list
   */
  std::vector<Eigen::MatrixXd> rot_m_list_;

  /**
   * @brief a list of (6 + 6) vectors in global frame
   * (N_all_element x (12)) list
   */
  std::vector<Eigen::VectorXd> element_lumped_nload_list_;

  /**
   * @brief a list of (6 + 6) vectors in global frame
   * 
   * Note: for now, we separate gravity and other element-wise loads
   * 
   * (N_all_element x (12)) list
   * 
   */
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
   * stored computation results
   */
  std::vector<int> stored_existing_ids_;

  Eigen::MatrixXd stored_nodal_deformation_;
  Eigen::MatrixXd stored_element_reaction_;
  Eigen::MatrixXd stored_fixities_reaction_;
  double stored_compliance_;

  /**
   * boolean flag for if the model is inited (1) or not (0).
   */
  bool is_init_;

  bool include_self_weight_load_;
  Eigen::VectorXd gravity_direction_;

  bool has_stored_deformation_;
};

} // namespace stiffness_checker
} // namespace conmech
