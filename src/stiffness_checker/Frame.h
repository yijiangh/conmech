#pragma once

#include <cassert>
#include <string>
#include <vector>
#include <memory>

#include <Eigen/Dense>

namespace conmech
{
namespace stiffness_checker
{

class FrameVert;
class FrameElement;

typedef std::shared_ptr<FrameVert> FrameVertPtr;
typedef std::shared_ptr<FrameElement> FrameElementPtr;

class FrameVert
{
 private:

 public:
  FrameVert()
      : id_(0), degree_(0) {}
  FrameVert(const Eigen::Vector3d& p)
      : position_(p),
        id_(0), degree_(0),
        b_fixed_(false) {}
  ~FrameVert() {}

 public:
  Eigen::Vector3d position() const { return position_; }
  int id() const { return id_; }
  int degree() const { return degree_; }
  const std::vector<FrameElementPtr> nghdElements() const { return p_nghd_element_list_; };

  // TODO: this "Fixed" is ambiguous, use "grounded" instead
  bool isFixed() const { return b_fixed_; }

  void setFixed(bool b_fixed) { b_fixed_ = b_fixed; }
  void setPosition(Eigen::Vector3d p) { position_ = p; }
  void setPosition(double x, double y, double z) { position_ = Eigen::Vector3d(x, y, z); }
  void setID(int id) { id_ = id; }
  void increaseDegree() { degree_++; }

  void addNghdElement(FrameElementPtr& e)
  {
    p_nghd_element_list_.push_back(e);
    degree_ = (int)p_nghd_element_list_.size();
  }
  void clearNghdElement()
  {
    p_nghd_element_list_.clear();
    degree_ = (int)p_nghd_element_list_.size();
  }

 private:
  std::vector<FrameElementPtr> p_nghd_element_list_;
  Eigen::Vector3d position_;

  int id_;
  int degree_;

  bool b_fixed_;
};

class FrameElement
{
 public:
  FrameElement()
      : id_(0), layer_(-1) {}
  ~FrameElement() {}

 public:
  int id() const { return id_; }
  int layer() const { return layer_; }

  void setID(int id) { id_ = id; }
  void setLayer(int layer) { layer_ = layer; }
  void setEndVertU(FrameVertPtr u) { pvert_u = u; }
  void setEndVertV(FrameVertPtr v) { pvert_v = v; }

  double getLength() const
  {
    if (pvert_u && pvert_v)
    {
      auto u = pvert_u->position();
      auto v = pvert_v->position();

      return (u-v).norm();
    }
    else
    {
      return 0;
    }
  }

  const FrameVertPtr endVertU() const
  {
    return pvert_u;
  }

  const FrameVertPtr endVertV() const
  {
    return pvert_v;
  }

 private:
  FrameVertPtr pvert_u;
  FrameVertPtr pvert_v;

  int id_;
  int layer_;
};

class Frame
{
 public:
  Frame(): layer_size_(0), fixed_vert_size_(0), unit_scale_(1.0) {}
  ~Frame() {}

 public:
  bool loadFromJson(const std::string& file_path);
  void clear();

  FrameVertPtr insertVertex(const Eigen::Vector3d& p);
  FrameElementPtr insertElement(FrameVertPtr u, FrameVertPtr v);

  void unify();

  inline void setUnitScale(double unit_scale) { unit_scale_ = unit_scale; }

  inline int sizeOfVertList() const { return (int)vert_list_.size(); }
  inline int sizeOfElementList() const { return (int)element_list_.size(); }
  inline int sizeOfFixedVert() const { return fixed_vert_size_; }
  inline int sizeOfLayer() const { return layer_size_; }

  const FrameVertPtr getVert(int vert_id) const;
  const std::vector<FrameElementPtr> getVertNeighborEdgeList(int v_id);
  Eigen::Vector3d getVertPosition(int v_id) const;
  int getVertDegree(int v_id) const;
  bool isVertFixed(int v_id) const;

  const FrameElementPtr getElement(int e_id) const;
  const FrameVertPtr getElementEndVertU(int e_id) const;
  const FrameVertPtr getElementEndVertV(int e_id) const;

  inline Eigen::Vector3d getCenterPos() const { return center_pos_; }
  inline Eigen::Vector3d getBaseCenterPos() const { return base_center_pos_; }

  inline double getUnitScale() const { return unit_scale_; }

 private:
  std::vector<FrameVertPtr> vert_list_;
  std::vector<FrameElementPtr> element_list_;

  int fixed_vert_size_;
  int layer_size_;

  double maxx_;
  double maxy_;
  double maxz_;
  double minx_;
  double miny_;
  double minz_;
  double basez_;

  Eigen::Vector3d center_pos_;
  Eigen::Vector3d base_center_pos_;

  /**
   * unit used inside Frame is millimeter.
   */
  double unit_scale_;
};

} // namespace stiffness_checker
} // namespace conmech