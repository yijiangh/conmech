using System;
using System.Collections;
using System.Collections.Generic;

using Rhino;
using Rhino.Geometry;

using Grasshopper;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;

using System.IO;
using Karamba.Models;
using Karamba.Nodes;
using Karamba.CrossSections;
using Karamba.Elements;
using Karamba.Loads;
using Karamba.Materials;
using Karamba.Supports;
using Karamba.Geometry;
using Karamba.Joints;
using Newtonsoft.Json;


/// <summary>
/// This class will be instantiated on demand by the Script component.
/// </summary>
public class Script_Instance : GH_ScriptInstance
{
#region Utility functions
  /// <summary>Print a String to the [Out] Parameter of the Script component.</summary>
  /// <param name="text">String to print.</param>
  private void Print(string text) { /* Implementation hidden. */ }
  /// <summary>Print a formatted String to the [Out] Parameter of the Script component.</summary>
  /// <param name="format">String format.</param>
  /// <param name="args">Formatting parameters.</param>
  private void Print(string format, params object[] args) { /* Implementation hidden. */ }
  /// <summary>Print useful information about an object instance to the [Out] Parameter of the Script component. </summary>
  /// <param name="obj">Object instance to parse.</param>
  private void Reflect(object obj) { /* Implementation hidden. */ }
  /// <summary>Print the signatures of all the overloads of a specific method to the [Out] Parameter of the Script component. </summary>
  /// <param name="obj">Object instance to parse.</param>
  private void Reflect(object obj, string method_name) { /* Implementation hidden. */ }
#endregion

#region Members
  /// <summary>Gets the current Rhino document.</summary>
  private readonly RhinoDoc RhinoDocument;
  /// <summary>Gets the Grasshopper document that owns this script.</summary>
  private readonly GH_Document GrasshopperDocument;
  /// <summary>Gets the Grasshopper script component that owns this script.</summary>
  private readonly IGH_Component Component;
  /// <summary>
  /// Gets the current iteration count. The first call to RunScript() is associated with Iteration==0.
  /// Any subsequent call within the same solution will increment the Iteration count.
  /// </summary>
  private readonly int Iteration;
#endregion

  /// <summary>
  /// This procedure contains the user code. Input parameters are provided as regular arguments,
  /// Output parameters as ref arguments. You don't have to assign output parameters,
  /// they will have a default value.
  /// </summary>
  private void RunScript(System.Object _k_model, List<System.Object> _k_grounded_nodes, string model_name, string save_path, bool export, double limit_distance, bool format_json, bool reimport_from_json, ref object A)
  {
    var kmodel = _k_model as Model;
    if (kmodel == null) {
      throw new ArgumentException("The input is not of type karamba model!");
    }
    if (limit_distance == 0) {
      limit_distance = 0.001; // meter
    }
    // support fixities & gounded node tags
    var k_grounded = _k_grounded_nodes.ConvertAll(x => x as Support);

    // Karamba element tag rules:
    // [] (empty list): non is specified, use the default cross_sec and material
    // [""] : one material is specified, apply to ALL elements

    //    foreach(var bs in kmodel.ploads) {
    //      Print(bs.loadcase.ToString());
    //      //      foreach(var i in bs.elemIds){
    //      //        Print("{0}", i);
    //      //      }
    //      Print("---");
    //    }

    var c_model = new ConmechModel(ref kmodel, grounded_supports: k_grounded,
    model_name: model_name, limit_distance: limit_distance);

    // save to json
    string cmodel_file_path = String.Concat(save_path, Path.DirectorySeparatorChar,
      model_name, ".json");
    if(export) {
      JsonSerializer serializer = new JsonSerializer();
      serializer.NullValueHandling = NullValueHandling.Ignore;
      if(format_json) {
        serializer.Formatting = Newtonsoft.Json.Formatting.Indented;
      }
      // serialize shape to a json file
      using (StreamWriter sw = new StreamWriter(@cmodel_file_path))
      using (JsonWriter writer = new JsonTextWriter(sw))
      {
        serializer.Serialize(writer, c_model);
      }
      Print("Shape saved to {0}", cmodel_file_path);

      //      // serialize load case to a json file
      //      string load_file_path = String.Concat(save_path, Path.DirectorySeparatorChar,
      //        model_name, "_loads.json");
      //      using (StreamWriter sw = new StreamWriter(@load_file_path))
      //      using (JsonWriter writer = new JsonTextWriter(sw))
      //      {
      //        serializer.Serialize(writer, c_loadcase);
      //      }
      //      Print("Load saved to {0}", load_file_path);
    }

    // reconstruct karamba model
    if(reimport_from_json) {
      ConmechModel cmodel = JsonConvert.DeserializeObject<ConmechModel>(File.ReadAllText(@cmodel_file_path));

      Print(cmodel.model_name);
    }
  }

  // <Custom additional code> 
  // TODO compile these into a dll
  // TODO compile these into a dll
  public class ConmechFrameNode
  {
    public ConmechFrameNode() {}
    public ConmechFrameNode(Node k_node) {
      node_ind = -1;
      is_grounded = false;
      node_ind = k_node.ind;
      point = new List<double>(){k_node.pos.X, k_node.pos.Y, k_node.pos.Z};
    }

    public List<double> point { get; set; }
    public int node_ind { get; set; }

    // is_grounded marks if the node is a support node
    // used to differentiate the supports to make a 2D structure stable
    public bool is_grounded { get; set; }

    public Node ToKNode() {
      return new Node(this.node_ind, new Point3(this.point[0], this.point[1], this.point[2]));
    }
  }

  public class ConmechFrameElement
  {
    public ConmechFrameElement() {}
    public ConmechFrameElement(ModelTruss ke) {
      // BuilderBeam does not have thenode_inds attributes
      elem_tag = ke.id;
      elem_ind = ke.ind;
      end_node_inds = ke.node_inds as List<int>;
      bending_stiff = ke.bending_stiff;
      // default layer index
      layer_ind = -1;
    }
    public List<int> end_node_inds { get; set; }
    // element identifier
    public string elem_tag { get; set; }
    // zero-based index of element in model
    public int elem_ind { get; set; }
    // an extra tag identify for grouping elements, can be used in sequencing
    public int layer_ind { get; set; }
    public bool bending_stiff { get; set; }

    // TODO: export local coordinate frame and import into conmech

    public BuilderBeam ToKBuilderBeam() {
      return new BuilderBeam(this.end_node_inds[0], this.end_node_inds[1]);
    }
  }

  public class ConmechFemMaterial
  {
    public ConmechFemMaterial() {}
    public ConmechFemMaterial(FemMaterial km) {
      // TODO: read unit from Karamba
      E_unit = "kN/m2";
      G12_unit = "kN/m2";
      fy_unit = "kN/m2";
      density_unit = "kN/m3";

      type_name = km.typeName();
      if (type_name != "ISO") {
        throw new ArgumentException("orthotropic material not implemented!!");
      }
      family = km.family;
      name = km.name;
      //  material's identifiers of elements to which it shall be applied
      elem_tags = km.elemIds;
      E = km.E();
      G12 = km.G12();
      fy = km.fy();
      density = km.gamma();
    }

    public string family { get; set; }
    public string name { get; set; }
    public string type_name { get; set; }
    public string _source { get; set; }

    // list of element identifiers to which the cross section is to be attached.
    // applied to all elements if the list is empty
    public List<string> elem_tags { get; set; }

    public double E { get; set; }
    public string E_unit { get; set; }

    // in-plane shear modulus
    public double G12 { get; set; }
    public string G12_unit { get; set; }

    // transverse shear modulus
    //    public double G3 { get; set; }
    //    public string G3_unit { get; set; }

    // tensile yield stress
    public double fy { get; set; }
    public string fy_unit { get; set; }

    // gamma
    public double density { get; set; }
    public string density_unit { get; set; }
  }

  public class ConmechCrossSec {
    public ConmechCrossSec() {}
    public ConmechCrossSec(CroSec_Beam kcs) {
      // TODO: read unit from Karamba
      A_unit = "m2";
      Iy_unit = "m4";
      Iz_unit = "m4";
      Jx_unit = "m4";

      family = kcs.family;
      name = kcs.name;
      elem_tags = kcs.elemIds;
      A = kcs.A;
      Jx = kcs.Ipp;
      Iy = kcs.Iyy;
      Iz = kcs.Izz;
    }

    public string family { get; set; }
    public string name { get; set; }
    public string _source { get; set; }

    // list of element identifiers to which the cross section is to be attached.
    // applied to all elements if the list is empty
    public List<string> elem_tags { get; set; }

    public double A { get; set; }
    public string A_unit { get; set; }
    public string _A_note { get; set; }

    // around local axis
    public double Jx { get; set; }
    public string Jx_unit { get; set; }
    public string _Jx_note { get; set; }

    public double Iy { get; set; }
    public string Iy_unit { get; set; }
    public string _Iy_note { get; set; }

    public double Iz { get; set; }
    public string Iz_unit { get; set; }
    public string _Iz_note { get; set; }
  }

  public class ConmechJoint
  {
    public ConmechJoint() {}
    public ConmechJoint(Joint kj){
      c_conditions = new List<List<Nullable<double>>>();
      c_conditions.Add(kj.c_condition(0));
      c_conditions.Add(kj.c_condition(1));
      elem_tags = kj.elemIds;
    }
    public List<List<Nullable<double>>> c_conditions { get; set; }
    // element tags
    public List<string> elem_tags { get; set; }
  }

  public class ConmechSupport
  {
    public ConmechSupport() {}
    public ConmechSupport(Support ks){
      condition = ks.Condition as List<bool>;
      node_ind = ks.node_ind;
    }
    public List<bool> condition { get; set; }
    public int node_ind { get; set; }
  }

  public class ConmechLoad
  {
    public int loadcase { get; set; }
  }

  public class ConmechPointLoad : ConmechLoad
  {
    public ConmechPointLoad() {}
    public ConmechPointLoad(PointLoad kf) {
      force = new List<double>(){kf.force.X, kf.force.Y, kf.force.Z};
      force_unit = "kN";
      moment = new List<double>(){kf.moment.X, kf.moment.Y, kf.moment.Z};
      moment_unit = "kNm";
      node_ind = kf.node_ind;
      loadcase = kf.loadcase;
    }
    // force-vector (kN) in global coordinates
    public List<double> force { get; set; }
    public string force_unit { get; set; }
    // moment-vector (kN) in global coordinates
    public List<double> moment { get; set; }
    public string moment_unit { get; set; }
    public int node_ind { get; set; }
  }

  public class ConmechUnifDistLoad : ConmechLoad
  {
    public ConmechUnifDistLoad() {}
    public ConmechUnifDistLoad(UniformlyDistLoad  kf) {
      load_unit = "kN/m";
      load = new List<double>(){kf.Load.X, kf.Load.Y, kf.Load.Z};
      // element tags
      elem_tags = kf.beamIds;
      loadcase = kf.loadcase;

      // Q = new List<double>(){kf.Q.X, kf.Q.Y, kf.Q.Z};
    }
    // load-vector in the global coordinate
    public List<double> load { get; set; }
    public string load_unit { get; set; }
    public List<string> elem_tags { get; set; }

    //    // vector of distributed loads
    //    public List<double> Q { get; set; }
    //    public string Q_unit { get; set; }
    //    // orientation of distributed loads
    //    public int Q_orient { get; set; }
  }

  public class ConmechGravityLoad : ConmechLoad
  {
    public ConmechGravityLoad() {}
    public ConmechGravityLoad(GravityLoad kf) {
      force = new List<double>(){kf.force.X, kf.force.Y, kf.force.Z};
      force_unit = "kN";
      loadcase = kf.loadcase;
    }
    // vector of acceleration due to gravity.
    // A vector length of "1" corresponds to gravity as on earth.
    public List<double> force { get; set; }
    public string force_unit { get; set; }
  }

  public class ConmechLoadCase
  {
    public int lc_ind { get; set; }
    public List<ConmechPointLoad> ploads { get; set; }
    // TODO: only support uniformly dist element load for now
    public List<ConmechUnifDistLoad> eloads { get; set; }
    public ConmechGravityLoad gravity { get; set; }
  }

  public class ConmechModel
  {
    // -----------
    // constructor
    public ConmechModel() {}
    public ConmechModel(ref Model kmodel,
    string model_name = "", List<Support> grounded_supports = null, double limit_distance = 0.001) {
      this.model_name = model_name;
      this.generate_time = DateTime.Now.ToString();
      this.unit = "meter";
      this.element_num = kmodel.elems.Count;
      this.node_num = kmodel.nodes.Count;

      // initiate nodes with id and position
      this.nodes = kmodel.nodes.ConvertAll(kn => new ConmechFrameNode(kn));

      // element data
      this.elements = kmodel.elems.ConvertAll(ke => new ConmechFrameElement(ke as ModelTruss));

      // support data
      this.supports = kmodel.supports.ConvertAll(ks => new ConmechSupport(ks as Support));

      // material data
      this.materials = kmodel.materials.ConvertAll(km => new ConmechFemMaterial(km));

      // cross section data
      this.cross_secs = kmodel.crosecs.ConvertAll(kc => new ConmechCrossSec(kc as CroSec_Beam));

      // joint data
      this.joints = kmodel.joints.ConvertAll(kj => new ConmechJoint(kj));

      // tag grounded
      if (grounded_supports != null) {
        foreach(Support k_supp in kmodel.supports) {
          foreach(Support k_g in grounded_supports) {
            if(k_supp.position.DistanceTo(k_g.position) < limit_distance) {
              this.nodes[k_supp.node_ind].is_grounded = true;
            }
          }
        }
      }

      // save to attributes
      this.loadcases = new Dictionary<int, ConmechLoadCase>();
      for(int lc = 0; lc < kmodel.numLC; lc++){
        // load data
        var c_ploads = new List<ConmechPointLoad>();
        //          k_model.ploads.ConvertAll(kl => new ConmechPointLoad(kl as PointLoad));
        foreach(var pl in kmodel.ploads){
          if(pl.loadcase == lc) {
            c_ploads.Add(new ConmechPointLoad(pl as PointLoad));
          }
        }

        var c_eloads = new List<ConmechUnifDistLoad>();
        // k_model.eloads.ConvertAll(kl => new ConmechUnifDistLoad(kl as UniformlyDistLoad));
        foreach(var el in kmodel.eloads){
          if(el.loadcase == lc && el.GetType().Equals(typeof(UniformlyDistLoad))) {
            c_eloads.Add(new ConmechUnifDistLoad(el as UniformlyDistLoad));
          }
        }

        ConmechGravityLoad c_grav = null;
        Karamba.Loads.GravityLoad k_grav;
        kmodel.gravities.TryGetValue(lc, out k_grav);
        if(k_grav != null){
          c_grav = new ConmechGravityLoad(k_grav);
        }

        var c_loadcase = new ConmechLoadCase();
        c_loadcase.lc_ind = lc;
        c_loadcase.ploads = c_ploads;
        c_loadcase.eloads = c_eloads;
        c_loadcase.gravity = c_grav;

        this.loadcases[lc] = c_loadcase;
      }
    }

    // -----------
    // attributes
    public string model_name { get; set; }
    public string unit { get; set; }
    public string generate_time { get; set; }

    public int node_num { get; set; }
    public int element_num { get; set; }

    public List<ConmechFrameNode> nodes { get; set; }
    public List<ConmechFrameElement> elements { get; set; }

    public List<ConmechSupport> supports { get; set; }
    public List<ConmechJoint> joints { get; set; }

    public List<ConmechCrossSec> cross_secs { get; set; }
    public List<ConmechFemMaterial> materials { get; set; }

    public Dictionary<int, ConmechLoadCase> loadcases { get; set; }

    // -----------
    // model conversion
    public void ToKarambaModel(ref List<Karamba.Loads.Load> loads, out Model kmodel, double limit_dist=0.001) {
      var logger = new Karamba.Utilities.MessageLogger();
      var k3d = new KarambaCommon.Toolkit();

      double mass;
      Point3 cog; // center of gravity
      bool flag;
      string info;

      var builder_elems = this.elements.ConvertAll(kl => kl.ToKBuilderBeam());

      //      kmodel = k3d.Model.AssembleModel(, this.supports, ploads,
      //        out info, out mass, out cog, out info, out flag,
      //        joints: joints, points: points,
      //      crosecs: crosecs, materials = materials,
      //      limitDist = limitDist, );

      //      AssembleModel(
      //        IReadOnlyList<BuilderElement> elems,
      //        IReadOnlyList<Support> supports,
      //        List<Load > loads,
      //      out string info,
      //      out double mass,
      //      out Point3 cog,
      //      out string msg,
      //      out bool runtimeWarning,
      //
      //        IReadOnlyList<Joint > joints = null,
      //      IReadOnlyList<Point3> points = null,
      //        IReadOnlyList<ElemSet > beamsets = null,
      //      double limitDist = 0,005,
      //      IReadOnlyList<CroSec> crosecs = null,
      //        IReadOnlyList<FemMaterial > materials = null
      //      )

      // replace later
      kmodel = new Model();
    }
  }

  //  public void KarambaModelFromConmechModel(ref ConmechModel cmodel, out Model kmodel)
  //  {
  //
  //    //    var tmp_nodes = new List<Point3>();
  //    //    var beams = new List<BuilderBeam>();
  //    //    foreach(ConmechFrameElement celem in cmodel.elements) {
  //    //      var end_node_inds = celem.end_node_inds;
  //    //      var cnode0 = cmodel.nodes[end_node_inds[0]];
  //    //      var cnode1 = cmodel.nodes[end_node_inds[1]];
  //    //      var p0 = new Point3(cnode0.point[0], cnode0.point[1], cnode0.point[2]);
  //    //      var p1 = new Point3(cnode1.point[0], cnode1.point[1], cnode1.point[2]);
  //    //      var elems = k3d.Part.LineToBeam(new List<Line3>(){new Line3(p0, p1)}, new List<string>(){ celem.elem_tag },
  //    //        new List<CroSec>(), logger, out tmp_nodes, celem.bending_stiff);
  //    //      beams.AddRange(elems);
  //    //    }
  //
  //    //    var cond = new List<bool>(){ true, true, true, true, true, true};
  //    //    var support = k3d.Support.Support(0, cond);
  //    //    var supports = new List<Support>(){support};
  //    //
  //    //    var pload = k3d.Load.PointLoad(1, new Vector3(0, 0, -10), new Vector3());
  //    //    var ploads = new List<Load>(){pload};
  //
  //    kmodel = new Model();
  //    //    Model(List<Node>, List<FemMaterial>, List<CroSec>, List<ModelElement>,
  //    //    List < PointLoad>, List < PointMass>, List < MeshLoad>, List < ElementLoad>,
  //    //    List < Support>, List < ElemSet>, Dictionary < Int32, GravityLoad>,
  //    //    ItemSelector, IdManager, Int32)
  //  }
  // </Custom additional code> 
}