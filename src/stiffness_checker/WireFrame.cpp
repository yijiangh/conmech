#include <iostream>
#include "stiffness_checker/WireFrame.hpp"

#include <rapidjson/document.h>
#include <rapidjson/filereadstream.h>

namespace
{
double convertUnitScale(const std::string& unit)
{
  if("millimeter" == unit || "mm" == unit)
  {
    return 0.001;
  }
  if("centimeter" == unit || "cm" == unit)
  {
    return 0.01;
  }
  if("meter" == unit || "m" == unit)
  {
    return 1;
  }
  if("inch" == unit || "in" == unit)
  {
    return 0.0254;
  }
  if("foot" == unit || "ft" == unit)
  {
    return 0.3048;
  }

  // default millimeter
  std::cout << "WARNING: unrecognized unit in the input json file. Using millimeter by default." << std::endl;
  return 0.001;
}

} // util anon ns

namespace conmech
{
namespace stiffness_checker
{

WireFrame::WireFrame()
    : delta_tol_(1e-1), unify_size_(2.0), layer_size_(0), unit_scale_(1.0)
{
  pvert_list_ = new std::vector<WF_vert*>;
  pedge_list_ = new std::vector<WF_edge*>;
}

WireFrame::~WireFrame()
{
  int N = pvert_list_->size();
  for (int i = 0; i < N; i++)
  {
    delete (*pvert_list_)[i];
    (*pvert_list_)[i] = NULL;
  }
  delete pvert_list_;
  pvert_list_ = NULL;

  int M = pedge_list_->size();
  for (int i = 0; i < M; i++)
  {
    delete (*pedge_list_)[i];
    (*pedge_list_)[i] = NULL;
  }
  delete pedge_list_;
  pedge_list_ = NULL;
}

void WireFrame::LoadFromPWF(const char *path)
{
  /* assert that there is no replication in PWF */

  FILE *fp = fopen(path, "r");

  assert(fp);

  try
  {
    // read vertexes
    fseek(fp, 0, SEEK_SET);
    char pLine[512];
    char *tok;
    while (fgets(pLine, 512, fp))
    {
      if (pLine[0] == 'v' && pLine[1] == ' ')
      {
        Vec3f p;
        char tmp[128];
        tok = strtok(pLine, " ");
        for (int i = 0; i < 3; i++)
        {
          tok = strtok(NULL, " ");
          strcpy(tmp, tok);
          tmp[strcspn(tmp, " ")] = 0;
          p[i] = (float) atof(tmp);
        }

        p = point(p.x(), p.y(), p.z());
        InsertVertex(p);
      }
    }

    // read layer
    fseek(fp, 0, SEEK_SET);
    while (fgets(pLine, 512, fp))
    {
      if (pLine[0] == 'g' && pLine[1] == ' ')
      {
        tok = strtok(pLine, " ");

        char tmp[128];
        tok = strtok(NULL, " ");
        strcpy(tmp, tok);
        tmp[strcspn(tmp, " ")] = 0;
        int u = (int) atof(tmp) - 1;

        tok = strtok(NULL, " ");
        strcpy(tmp, tok);
        tmp[strcspn(tmp, " ")] = 0;
        int v = (int) atof(tmp) - 1;

        tok = strtok(NULL, " ");
        strcpy(tmp, tok);
        tmp[strcspn(tmp, " ")] = 0;
        int layer = (int) atof(tmp);

        if (!(0 <= u && u < pvert_list_->size()))
        {
          std::cout << "read layer: start node id overflow: " << u << "/" << pvert_list_->size() << std::endl;
          assert(0 <= u && u < pvert_list_->size());
        }
        if (!(0 <= v && v < pvert_list_->size()))
        {
          std::cout << "read layer: end node id overflow: " << v << "/" << pvert_list_->size() << std::endl;
          assert(0 <= v && v < pvert_list_->size());
        }

        WF_edge *e = InsertEdge((*pvert_list_)[u], (*pvert_list_)[v]);
        if (e != NULL)
        {
          if (e->Layer() != -1)
          {
            std::cout << "Overwrite previously set id! - prev layer id : " << e->Layer()
                      << ", current id: " << layer << std::endl;
            assert(e->Layer() == -1);
          }

          e->SetLayer(layer);
          e->ppair_->SetLayer(layer);
        }
      }
    }

    // read lines
    char c;
    int prev;
    int curv;
    fseek(fp, 0, SEEK_SET);
    while (c = fgetc(fp), c != EOF)
    {
      while (c != 'l' && c != EOF)
      {
        c = fgetc(fp);
      }

      if (c == '\n' || c == EOF || (c = fgetc(fp)) != ' ')
      {
        continue;
      }

      prev = -1;
      while (c != '\n' && c != EOF)
      {
        while (c = fgetc(fp), c != '\n' && c != EOF && !isdigit(c));

        if (c == '\n' || c == EOF)
        {
          break;
        }

        for (curv = 0; isdigit(c); c = fgetc(fp))
        {
          curv = curv * 10 + c - '0';
        }
        curv--;

        if (prev != -1)
        {
          if (!(0 <= prev && prev < pvert_list_->size()))
          {
            std::cout << "read lines: start node id overflow: " << prev << "/" << pvert_list_->size() << std::endl;
            assert(0 <= prev && prev < pvert_list_->size());
          }
          if (!(0 <= curv && curv < pvert_list_->size()))
          {
            std::cout << "read lines: end node id overflow: " << curv << "/" << pvert_list_->size() << std::endl;
            assert(0 <= curv && curv < pvert_list_->size());
          }

          InsertEdge((*pvert_list_)[prev], (*pvert_list_)[curv]);
        }

        prev = curv;
      }
    }

    // read pillar
    fseek(fp, 0, SEEK_SET);
    while (fgets(pLine, 512, fp))
    {
      if (pLine[0] == 'p' && pLine[1] == ' ')
      {
        tok = strtok(pLine, " ");

        char tmp[128];
        tok = strtok(NULL, " ");
        strcpy(tmp, tok);
        tmp[strcspn(tmp, " ")] = 0;
        int u = (int) atof(tmp) - 1;

        tok = strtok(NULL, " ");
        strcpy(tmp, tok);
        tmp[strcspn(tmp, " ")] = 0;
        int v = (int) atof(tmp) - 1;

        WF_vert *b = (*pvert_list_)[u];
        b->SetBase(true);

        WF_vert *f = (*pvert_list_)[v];
        f->SetFixed(true);

        WF_edge *e = InsertEdge(b, f);
        if (e != NULL)
        {
          e->SetPillar(true);
          e->ppair_->SetPillar(true);
        }
      }
    }

    // read ceiling
    fseek(fp, 0, SEEK_SET);
    while (fgets(pLine, 512, fp))
    {
      if (pLine[0] == 'c' && pLine[1] == ' ')
      {
        tok = strtok(pLine, " ");

        char tmp[128];
        tok = strtok(NULL, " ");
        strcpy(tmp, tok);
        tmp[strcspn(tmp, " ")] = 0;
        int u = (int) atof(tmp) - 1;

        tok = strtok(NULL, " ");
        strcpy(tmp, tok);
        tmp[strcspn(tmp, " ")] = 0;
        int v = (int) atof(tmp) - 1;

        WF_edge *e = InsertEdge((*pvert_list_)[u], (*pvert_list_)[v]);
        if (e != NULL)
        {
          e->SetCeiling(true);
          e->ppair_->SetCeiling(true);
        }
      }
    }

    Unify();
  }
  catch (...)
  {
    return;
  }

  fclose(fp);
}

bool WireFrame::LoadFromJson(const std::string &file_path)
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
    return false;
  }

  fclose(fp);

  if (!document.HasMember("unit"))
  {
    cout << "WARNING: no unit is specified - using millimeter by default" << endl;
  }
  else
  {
    // convert to meter (SI)
    unit_scale_ = convertUnitScale(document["unit"].GetString());
  }

  // read vertices
  assert(document.HasMember("node_list"));
  assert(document["node_list"].Size() > 0);

  for(int i=0; i<document["node_list"].Size(); i++)
  {
    const Value& p = document["node_list"][i]["point"];
    WF_vert* vert = InsertVertex(
        point(p["X"].GetDouble()*unit_scale_,
              p["Y"].GetDouble()*unit_scale_,
              p["Z"].GetDouble()*unit_scale_));

    if(document["node_list"][i]["is_grounded"].GetInt())
    {
      vert->SetFixed(true);
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

    if (!(0 <= u && u < pvert_list_->size()))
    {
      std::cout << "ERROR: read edge: start node id overflow: " << u << "/" << pvert_list_->size() << std::endl;
      assert(0 <= u && u < pvert_list_->size());
    }
    if (!(0 <= v && v < pvert_list_->size()))
    {
      std::cout << "ERROR read edge: end node id overflow: " << v << "/" << pvert_list_->size() << std::endl;
      assert(0 <= v && v < pvert_list_->size());
    }

    WF_edge *e = InsertEdge((*pvert_list_)[u], (*pvert_list_)[v]);
    if (e != NULL)
    {
      if (e->Layer() != -1)
      {
        std::cout << "ERROR: Overwrite previously set id! - prev layer id : " << e->Layer()
                  << ", current id: " << layer << std::endl;
        assert(e->Layer() == -1);
      }

      e->SetLayer(layer);
      e->ppair_->SetLayer(layer);

      if((*pvert_list_)[u]->isFixed() || (*pvert_list_)[v]->isFixed())
      {
        e->SetPillar(true);
        e->ppair_->SetPillar(true);
      }
    }
  }

  Unify();

  return true;
}

void WireFrame::WriteToPWF(
    bool bVert, bool bLine,
    bool bPillar, bool bCeiling,
    bool bCut, int min_layer, int max_layer,
    const char *path)
{
  using namespace std;

  if ((*pedge_list_)[0]->Layer() == -1 || !bCut)
  {
    min_layer = -(1 << 20);
    max_layer = (1 << 20);
  }

  FILE *fp = fopen(path, "wb+");
  int N = SizeOfVertList();
  int M = SizeOfEdgeList();

  vector<int> export_vert;
  vector<int> hash(N);
  fill(hash.begin(), hash.end(), -1);

  for (int i = 0; i < M; i++)
  {
    WF_edge *e1 = (*pedge_list_)[i];
    WF_edge *e2 = e1->ppair_;
    if (i < e2->ID() && e1->Layer() >= min_layer - 1 && e1->Layer() < max_layer)
    {
      if (e1->isPillar() && !bPillar)
      {
        continue;
      }

      int u = e2->pvert_->ID();
      int v = e1->pvert_->ID();
      if (hash[u] == -1)
      {
        export_vert.push_back(u);
        hash[u] = export_vert.size();
      }
      if (hash[v] == -1)
      {
        export_vert.push_back(v);
        hash[v] = export_vert.size();
      }
    }
  }

  int Ne = export_vert.size();
  if (bVert)
  {
    for (int i = 0; i < Ne; i++)
    {
      point p = (*pvert_list_)[export_vert[i]]->Position();
      fprintf(fp, "v %lf %lf %lf\n", p.x(), p.y(), p.z());
    }
  }

  if (bLine)
  {
    for (int i = 0; i < M; i++)
    {
      WF_edge *e1 = (*pedge_list_)[i];
      WF_edge *e2 = e1->ppair_;
      if (i < e2->ID() && !e1->isPillar() &&
          e1->Layer() >= min_layer - 1 && e1->Layer() < max_layer)
      {
        int u = e2->pvert_->ID();
        int v = e1->pvert_->ID();
        fprintf(fp, "l %d %d\n", hash[u], hash[v]);
      }
    }
  }

  if (bPillar)
  {
    for (int i = 0; i < M; i++)
    {
      WF_edge *e1 = (*pedge_list_)[i];
      WF_edge *e2 = e1->ppair_;
      if (i < e2->ID() && e1->isPillar()
          && e1->Layer() >= min_layer - 1 && e1->Layer() < max_layer)
      {
        int u = e2->pvert_->ID();
        int v = e1->pvert_->ID();
        assert(e2->pvert_->isBase());
        fprintf(fp, "p %d %d\n", hash[u], hash[v]);
      }
    }
  }

  if (bCeiling)
  {
    for (int i = 0; i < M; i++)
    {
      WF_edge *e1 = (*pedge_list_)[i];
      WF_edge *e2 = e1->ppair_;
      if (i < e2->ID() && e1->isCeiling()
          && e1->Layer() >= min_layer - 1 && e1->Layer() < max_layer)
      {
        int u = e2->pvert_->ID();
        int v = e1->pvert_->ID();
        fprintf(fp, "c %d %d\n", hash[u], hash[v]);
      }
    }
  }

  if (bCut)
  {
    for (int i = 0; i < M; i++)
    {
      WF_edge *e1 = (*pedge_list_)[i];
      WF_edge *e2 = e1->ppair_;
      if (i < e2->ID() && e1->Layer() != -1
          && e1->Layer() >= min_layer - 1 && e1->Layer() < max_layer)
      {
        int u = e2->pvert_->ID();
        int v = e1->pvert_->ID();
        fprintf(fp, "g %d %d %d\n", hash[u], hash[v], e1->Layer());
      }
    }
  }

  fclose(fp);
}

WF_vert *WireFrame::InsertVertex(Vec3f p)
{
  // detect duplication
  int N = SizeOfVertList();
  for (int i = 0; i < N; i++)
  {
    if (Dist(p, (*pvert_list_)[i]->Position()) < 1e-3)
    {
      std::cout << "duplicate point read! (" << p.x() << ", " << p.y() << ", " << p.z() << ")" << std::endl;
      return (*pvert_list_)[i];
    }
  }

  WF_vert *vert = new WF_vert(p);
  pvert_list_->push_back(vert);
  vert->SetID(N);
  return vert;
}

WF_edge *WireFrame::InsertEdge(WF_vert *u, WF_vert *v)
{
  // detect duplication
  WF_edge *e = u->pedge_;
  while (e != NULL)
  {
    if (e->pvert_ == v)
    {
      return e;
    }
    e = e->pnext_;
  }

  WF_edge *e1 = InsertOneWayEdge(u, v);
  WF_edge *e2 = InsertOneWayEdge(v, u);
  if (e1 != NULL)
  {
    e1->ppair_ = e2;
    e2->ppair_ = e1;
  }
  return e1;
}

WF_edge *WireFrame::InsertOneWayEdge(WF_vert *u, WF_vert *v)
{
  if (u == v)
  {
    return NULL;
  }

  WF_edge *edge = new WF_edge();
  edge->pvert_ = v;
  edge->pnext_ = u->pedge_;
  u->pedge_ = edge;
  u->IncreaseDegree();

  pedge_list_->push_back(edge);
  return edge;
}

void WireFrame::Unify()
{
  maxx_ = -1e20;
  maxy_ = -1e20;
  maxz_ = -1e20;
  minx_ = 1e20;
  miny_ = 1e20;
  minz_ = 1e20;
  basez_ = 1e20;

  fixed_vert_ = 0;
  base_vert_ = 0;
  pillar_size_ = 0;
  ceiling_size_ = 0;
  layer_size_ = 0;

  int N = SizeOfVertList();
  for (int i = 0; i < N; i++)
  {
    (*pvert_list_)[i]->SetID(i);
    point p = (*pvert_list_)[i]->Position();

    if (!(*pvert_list_)[i]->isFixed())
    {
      if (p.x() > maxx_)
      {
        maxx_ = p.x();
      }
      if (p.y() > maxy_)
      {
        maxy_ = p.y();
      }
      if (p.z() > maxz_)
      {
        maxz_ = p.z();
      }
      if (p.x() < minx_)
      {
        minx_ = p.x();
      }
      if (p.y() < miny_)
      {
        miny_ = p.y();
      }
      if (p.z() < minz_)
      {
        minz_ = p.z();
      }
    }

    //set base_z for once
    if (p.z() < basez_)
    {
      basez_ = p.z();
    }
  }

  int M = SizeOfEdgeList();
  for (int i = 0; i < M; i++)
  {
    (*pedge_list_)[i]->SetID(i);
    if ((*pedge_list_)[i]->isPillar())
    {
      pillar_size_++;
    }
    if ((*pedge_list_)[i]->isCeiling())
    {
      ceiling_size_++;
    }
    if ((*pedge_list_)[i]->Layer() > layer_size_)
    {
      layer_size_ = (*pedge_list_)[i]->Layer();
    }
  }

  if (pillar_size_ > 0)
  {
    // TODO: this assumption might be relaxed for future work
    // pillar is always in the first layer if exists
    if (0 == layer_size_)
    {
      layer_size_ = 1;
    }

    for (int i = 0; i < M; i++)
    {
      if (!(*pedge_list_)[i]->isPillar())
      {
        if (1 == layer_size_)
        {
          (*pedge_list_)[i]->SetLayer(1);
        }
        else
        {
          assert(0 != (*pedge_list_)[i]->Layer());
        }
      }
      else
      {
        (*pedge_list_)[i]->SetLayer(0);
      }
    }

    // add base
    layer_size_++;
  }

  float scaleX = maxx_ - minx_;
  float scaleY = maxy_ - miny_;
  float scaleZ = maxz_ - minz_;
  float scaleMax = scaleX;
  if (scaleMax < scaleY)
  {
    scaleMax = scaleY;
  }
  if (scaleMax < scaleZ)
  {
    scaleMax = scaleZ;
  }

  scaleV_ = unify_size_ / scaleMax;
  center_pos_ = point((minx_ + maxx_) / 2.f, (miny_ + maxy_) / 2.f, (minz_ + maxz_) / 2.f);
  base_center_pos_ = point((minx_ + maxx_) / 2.f, (miny_ + maxy_) / 2.f, basez_);

  for (size_t i = 0; i < N; i++)
  {
    (*pvert_list_)[i]->SetRenderPos(((*pvert_list_)[i]->Position() - center_pos_) * scaleV_);
    if ((*pvert_list_)[i]->isFixed())
    {
      fixed_vert_++;
    }
    if ((*pvert_list_)[i]->isBase())
    {
      base_vert_++;
    }
  }
}

} // ns stiffness_checker
} // ns conmech