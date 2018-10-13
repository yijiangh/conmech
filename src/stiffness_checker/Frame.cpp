#include <iostream>
#include <rapidjson/document.h>
#include <rapidjson/filereadstream.h>

#include <eigen3/Eigen/Dense>

#include "stiffness_checker/Frame.h"

namespace
{
double convertUnitScale(const std::string& unit)
{
  if("millimeter" == unit || "mm" == unit)
  {
    return 1;
  }
  if("centimeter" == unit || "cm" == unit)
  {
    return 10;
  }
  if("meter" == unit || "m" == unit)
  {
    return 100;
  }
  if("inch" == unit || "in" == unit)
  {
    return 25.4;
  }
  if("foot" == unit || "ft" == unit)
  {
    return 304.8;
  }

  // default millimeter
  std::cout << "WARNING: unrecognized unit in the input json file. Using millimeter by default." << std::endl;
  return 1;
}
} // util anon ns

namespace conmech
{
namespace stiffness_checker
{

const FrameVertPtr Frame::getVert(int vert_id) const
{
  return (vert_id >= sizeOfVertList() || vert_id < 0) ? NULL : vert_list_[vert_id];
}

const FrameElementPtr Frame::getElement(int e_id) const
{
  return (e_id >= sizeOfElementList() || e_id < 0) ? NULL : element_list_[e_id];
}

const std::vector<FrameElementPtr> Frame::getVertNeighborEdgeList(int v_id)
{
  assert(v_id >= sizeOfVertList() || v_id < 0);
  return vert_list_[v_id]->nghdElements();
}

Eigen::Vector3d Frame::getVertPosition(int v_id) const
{
  assert(v_id < sizeOfVertList() && v_id >= 0);
  return vert_list_[v_id]->position();
}

int Frame::getVertDegree(int v_id) const
{
  assert(v_id < sizeOfVertList() && v_id >= 0);
  return vert_list_[v_id]->degree();
}

const FrameVertPtr Frame::getElementEndVertU(int e_id) const
{
  assert(e_id < sizeOfElementList() && e_id >= 0);
  return element_list_[e_id]->endVertU();
}

const FrameVertPtr Frame::getElementEndVertV(int e_id) const
{
  assert(e_id < sizeOfElementList() && e_id >= 0);
  return element_list_[e_id]->endVertV();
}

bool Frame::isVertFixed(int v_id) const
{
  assert(v_id < sizeOfVertList() && v_id >= 0);
  return vert_list_[v_id]->isFixed();
}

bool Frame::loadFromJson(const std::string &file_path)
{
  using namespace std;
  using namespace rapidjson;

  FILE *fp = fopen(file_path.c_str(), "r");
  assert(fp);

  char readBuffer[65536];
  FileReadStream is(fp, readBuffer, sizeof(readBuffer));

  Document document;
  if (document.ParseStream(is).HasParseError())
  {
    std::cout << "ERROR parsing the input json file!\n";
    fclose(fp);
    return false;
  }

  fclose(fp);

  // reset all existing vert and element list
  clear();

  if (!document.HasMember("unit"))
  {
    cout << "WARNING: no unit is specified - using millimeter by default" << endl;
  }
  else
  {
    // convert to millimeter
    unit_scale_ = convertUnitScale(document["unit"].GetString());
  }

  // read vertices
  assert(document.HasMember("node_list"));
  assert(document["node_list"].Size() > 0);

  for(int i=0; i<document["node_list"].Size(); i++)
  {
    const Value& p = document["node_list"][i]["point"];
    FrameVertPtr vert = insertVertex(
        Eigen::Vector3d(p["X"].GetDouble()*unit_scale_,
                        p["Y"].GetDouble()*unit_scale_,
                        p["Z"].GetDouble()*unit_scale_));

    if(document["node_list"][i]["is_grounded"].GetInt())
    {
      vert->setFixed(true);
      fixed_vert_size_++;
    }
  }

  // read edges (beams)
  assert(document.HasMember("element_list"));
  assert(document["element_list"].Size() > 0);

  for(int i=0; i<document["element_list"].Size(); i++)
  {
    int u = document["element_list"][i]["end_node_ids"][0].GetInt();
    int v = document["element_list"][i]["end_node_ids"][1].GetInt();
    int layer = document["element_list"][i]["layer_id"].GetInt();


    assert(0 <= u && u < vert_list_.size());
    assert(0 <= v && v < vert_list_.size());

    FrameElementPtr e = insertElement(vert_list_[u], vert_list_[v]);

    if (e != NULL)
    {
      assert(e->layer() == -1 && "Error: overwriting perviously assigned element's layer id");
      e->setLayer(layer);

      if(layer>layer_size_-1)
      {
        layer_size_++;
      }
    }
  }

  unify();
  return true;
}

void Frame::clear()
{
  for(auto e : element_list_)
  {
    e.reset();
  }
  for(auto v : vert_list_)
  {
    v.reset();
  }

  element_list_.clear();
  vert_list_.clear();

  layer_size_ = 0;
  fixed_vert_size_ = 0;
}

FrameVertPtr Frame::insertVertex(const Eigen::Vector3d& p)
{
  // detect duplication
  auto N = (int)sizeOfVertList();

  for (int i = 0; i < N; i++)
  {
    assert((p-vert_list_[i]->position()).norm() > 1e-5 && "duplicate point detected.");
  }

  FrameVertPtr pvert = FrameVertPtr(new FrameVert(p));
  vert_list_.push_back(pvert);
  pvert->setID(N);

  return pvert;
}

FrameElementPtr Frame::insertElement(FrameVertPtr u, FrameVertPtr v)
{
  // duplication detection
  auto m = (int)sizeOfElementList();

  for (int i = 0; i < m; i++)
  {
    auto u_prev = element_list_[i]->endVertU();
    auto v_prev = element_list_[i]->endVertV();

    bool match_1 = ((u_prev->position() - u->position()).norm() < 1e-5)
        && ((v_prev->position() - v->position()).norm() < 1e-5);

    bool match_2 = ((u_prev->position() - v->position()).norm() < 1e-5)
        && ((v_prev->position() - u->position()).norm() < 1e-5);

    assert((!match_1) && (!match_2) && "duplicate element detected.");
  }

  FrameElementPtr e = FrameElementPtr(new FrameElement());
  e->setEndVertU(u);
  e->setEndVertV(v);

  u->addNghdElement(e);
  v->addNghdElement(e);

  e->setID(m);

  element_list_.push_back(e);

  return e;
}

void Frame::unify()
{
  maxx_ = maxy_ = maxz_ = -1e20;
  minx_ = miny_ = minz_ = 1e20;
  basez_ = 1e20;

  layer_size_ = 0;

  auto N = sizeOfVertList();

  for (int i = 0; i < N; i++)
  {
    vert_list_[i]->setID(i);
    auto p = vert_list_[i]->position();

    // update min,max xyz
    if (!vert_list_[i]->isFixed())
    {
      if (p[0] > maxx_)
      {
        maxx_ = p[0];
      }
      if (p[1] > maxy_)
      {
        maxy_ = p[1];
      }
      if (p[2] > maxz_)
      {
        maxz_ = p[2];
      }
      if (p[0] < minx_)
      {
        minx_ = p[0];
      }
      if (p[1] < miny_)
      {
        miny_ = p[1];
      }
      if (p[2] < minz_)
      {
        minz_ = p[2];
      }
    }

    //set base_z for once
    if (p[2] < basez_)
    {
      basez_ = p[2];
    }
  }

  center_pos_ = Eigen::Vector3d((minx_ + maxx_) / 2, (miny_ + maxy_) / 2, (minz_ + maxz_) / 2);
  base_center_pos_ = Eigen::Vector3d((minx_ + maxx_) / 2, (miny_ + maxy_) / 2, basez_);
}

} // ns stiffness_checker
} // ns conmech