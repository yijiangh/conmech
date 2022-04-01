import os, sys

assert 'ironpython' in sys.version.lower() and os.name == 'nt', 'Karamba only available for IronPython on Windows now.'

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(THIS_DIR))
from io_base import Node, Element, Support, Joint, CrossSec, Material, Model, LoadCase, \
    PointLoad, UniformlyDistLoad, GravityLoad, get_defalut_element_id_format

###############################################

# https://ironpython.net/documentation/dotnet/dotnet.html
import clr
try:
    clr.AddReferenceToFileAndPath("C:\Program Files\Rhino 7\Plug-ins\Karamba.gha")
    clr.AddReferenceToFileAndPath("C:\Program Files\Rhino 7\Plug-ins\KarambaCommon.dll")
except:
    try:
        clr.AddReferenceToFileAndPath("C:\Program Files\Rhino 7\Plug-ins\Karamba\Karamba.gha")
        clr.AddReferenceToFileAndPath("C:\Program Files\Rhino 7\Plug-ins\Karamba\KarambaCommon.dll")
    except:
        raise RuntimeError("Karamba Plug-in path (where Karamba.gha and KarambaCommon.dll resides) not found!")

# https://scripting.karamba3d.com/
import feb # Karamba's C++ library (undocumented in the API)
import Karamba
import Karamba.Models.Model
from Karamba.Geometry import Point3, Line3, Vector3, Plane3
import Karamba.Nodes
import Karamba.Elements
import Karamba.CrossSections
import Karamba.Materials
import Karamba.Supports
import Karamba.Loads
import KarambaCommon
#
import System
from System import GC
from System.Collections.Generic import List
from System.Drawing import Color

####################################

class KarambaNode(Node):
    @classmethod
    def from_karamba(cls, knode, is_grounded=False):
        point = [knode.pos.X, knode.pos.Y, knode.pos.Z]
        return cls(point, knode.ind, is_grounded)

    def to_karamba(self, type='Point3'):
        if type == 'Point3':
            return Point3(*self.point)
        elif type == 'Node':
            return Karamba.Nodes.Node(self.node_ind, Point3(*self.point))
        else:
            raise ValueError

class KarambaElement(Element):
    @classmethod
    def from_karamba(cls, kelement):
        if kelement.id == '' and kelement.crosec:
            # * if element has attribute crosec, we create an element identifier for it 
            # https://www.karamba3d.com/help/2-2-0/html/652823a9-ab3d-f789-e9bc-d5e71a03d372.htm
            element_id = get_defalut_element_id_format(kelement.ind)
        else:
            element_id = kelement.id
        return cls(list(kelement.node_inds), kelement.ind, elem_tag=element_id,
            bending_stiff=kelement.bending_stiff)

    def to_karamba(self, type='builderbeam'):
        if type == 'builderbeam':
            bb = Karamba.Elements.BuilderBeam(*self.end_node_inds)
            bb.id = self.elem_tag
            bb.bending_stiff = self.bending_stiff
            return bb
        else:
            raise NotImplementedError

class KarambaSupport(Support):
    @classmethod
    def from_karamba(cls, ksupport):
        return cls(list(ksupport.Condition), ksupport.node_ind)

    def to_karamba(self):
        return Karamba.Supports.Support(self.node_ind, List[bool](list(map(bool, self.condition))), Plane3())

class KarambaJoint(Joint):
    @classmethod
    def from_karamba(cls, kjoint):
        assert kjoint.dofs == 6
        # condition is If null there is no joint for the corresponding DOF.
        return cls(list(kjoint.c_condition(0)) + list(kjoint.c_condition(1)), list(kjoint.elemIds))

    def to_karamba(self):
        jt = Karamba.Joints.Joint()
        jt.elemIds = List[str](self.elem_tags)
        jt.c = System.Array[System.Nullable[float]](self.c_conditions)
        return jt

class KarambaCrossSec(CrossSec):
    @classmethod
    def from_karamba(cls, kc):
        return cls(kc.A, kc.Ipp, kc.Iyy, kc.Izz, elem_tags=list(kc.elemIds), family=kc.family, name=kc.name)

    def to_karamba(self):
        beam_mod = Karamba.CrossSections.CroSec_BeamModifier()
        beam_mod.A = self.A
        beam_mod.Ipp = self.Jx
        beam_mod.Iyy = self.Iy
        beam_mod.Izz = self.Iz
        beam_mod.elemIds = List[str]([tag for tag in self.elem_tags if tag is not None])
        # beam_mod.clearElemIds()
        # for tag in self.elem_tags:
        #     if tag is not None:
        #         beam_mod.AddElemId(tag)
        beam_mod.family = self.family
        beam_mod.name = self.name
        return beam_mod

class KarambaMaterialIsotropic(Material):
    @classmethod
    def from_karamba(cls, km):
        assert km.typeName() == 'ISO', 'Input should be a Isotropic material!'
        return cls(km.E(), km.G12(), km.fc(), km.ft(), km.gamma(), elem_tags=list(km.elemIds),
            family=km.family, name=km.name, type_name=km.typeName(), G3=km.G3())

    def to_karamba(self):
        steel_alphaT = 1e-5 #1/deg
        # ! assuming equal compressive and transile strength now
        km = Karamba.Materials.FemMaterial_Isotrop(
	        self.family, self.name, self.E, self.G12, self.G3, self.density, self.ft, self.fc,
                Karamba.Materials.FemMaterial.FlowHypothesis.mises, steel_alphaT, None) #Color.FromName("SlateBlue")
        for tag in self.elem_tags:
            if tag is not None:
                km.AddBeamId(tag)
        return km

class KarambaPointLoad(PointLoad):
    @classmethod
    def from_karamba(cls, kf):
        return cls([kf.force.X, kf.force.Y, kf.force.Z],
                   [kf.moment.X, kf.moment.Y, kf.moment.Z], kf.node_ind, kf.LcName)

    def to_karamba(self, corotate=False):
        # corotate = true the pointload corotates with the node
        pl = Karamba.Loads.PointLoad(self.node_ind, Vector3(*self.force), Vector3(*self.moment), str(self.loadcase), corotate)
        return pl

class KarambaUniformlyDistLoad(UniformlyDistLoad):
    @classmethod
    def from_karamba(cls, kf):
        return cls([kf.Q.X, kf.Q.Y, kf.Q.Z],
                   [kf.Load.X, kf.Load.Y, kf.Load.Z], kf.beamIds, loadcase=kf.LcName)

    def to_karamba(self):
        # Karamba.Loads.LoadOrientation.global
        ul = Karamba.Loads.UniformlyDistLoad(List[str](self.elem_tags), Vector3(*self.q), 0, self.loadcase)
        return ul

class KarambaGravityLoad(GravityLoad):
    @classmethod
    def from_karamba(cls, kf):
        return cls([kf.force.X, kf.force.Y, kf.force.Z], kf.LcName)

    def to_karamba(self):
        return Karamba.Loads.GravityLoad(Vector3(*self.force), self.loadcase)

class KarambaLoadCase(LoadCase):
    @classmethod
    def from_loadcase(cls, loadcase):
        kpoint_loads = [KarambaPointLoad.from_data(pl.to_data()) for pl in loadcase.point_loads]
        kelem_loads = [KarambaUniformlyDistLoad.from_data(pl.to_data()) for pl in loadcase.uniform_element_loads]
        kgravity_load = None
        if loadcase.gravity_load is not None:
            kgravity_load = KarambaGravityLoad.from_data(loadcase.gravity_load.to_data())
        return cls(kpoint_loads, kelem_loads, kgravity_load)

    @classmethod
    def from_karamba(cls, km, lc=0):
        assert lc < km.numLC
        point_loads = []
        for pl in km.ploads:
            if pl.LcName == lc or pl.LcName == str(lc):
                point_loads.append(KarambaPointLoad.from_karamba(pl))

        uniform_element_loads = []
        for el in km.eloads:
            if (el.LcName == lc or el.LcName == str(lc)) and isinstance(el, Karamba.Loads.UniformlyDistLoad):
                uniform_element_loads.append(KarambaUniformlyDistLoad.from_karamba(el))

        gravity_load = None
        k_grav = None
        # k_grav = clr.Reference[Karamba.Loads.GravityLoad]()
        # km.gravities.TryGetValue(str(lc), k_grav)
        if km.gravities.ContainsKey(str(lc)):
            k_grav = km.gravities[str(lc)]
        if k_grav is not None:
            gravity_load = KarambaGravityLoad.from_karamba(k_grav)

        return cls(point_loads, uniform_element_loads, gravity_load)

    def to_karamba(self):
        in_loads = List[Karamba.Loads.Load]()
        in_loads.AddRange([pl.to_karamba() for pl in self.point_loads])
        in_loads.AddRange([pl.to_karamba() for pl in self.uniform_element_loads])
        if self.gravity_load is not None:
            in_loads.Add(self.gravity_load.to_karamba())
        return in_loads

####################################

class KarambaModel(Model):
    @classmethod
    def from_model(cls, model):
        knodes = [KarambaNode.from_data(n.to_data()) for n in model.nodes]
        kelements = [KarambaElement.from_data(e.to_data()) for e in model.elements]
        ksupports = [KarambaSupport.from_data(s.to_data()) for s in model.supports.values()]
        kjoints = [KarambaJoint.from_data(j.to_data()) for j in model.joints]
        kmaterials = [KarambaMaterialIsotropic.from_data(m.to_data()) for m in model.materials]
        kcrosssecs = [KarambaCrossSec.from_data(cs.to_data()) for cs in model.crosssecs]
        return cls(knodes, kelements, ksupports, kjoints, kmaterials, kcrosssecs,
            unit=model.unit, model_name=model.model_name)

    @classmethod
    def from_data(cls, data, verbose=False):
        model = super(KarambaModel, cls).from_data(data, verbose)
        return cls.from_model(model)

    @classmethod
    def from_karamba(cls, karamba_model, grounded_supports=None, model_name='', limit_distance=1e-6):
        # grounded_supports: list of Karamba.Supports.Support
        knodes = [KarambaNode.from_karamba(kn) for kn in karamba_model.nodes]
        kmaterials = [KarambaMaterialIsotropic.from_karamba(m) for m in karamba_model.materials]
        kcrosssecs = [KarambaCrossSec.from_karamba(cs) for cs in karamba_model.crosecs]

        kelements = [KarambaElement.from_karamba(e) for e in karamba_model.elems]
        for ke, e in zip(kelements, karamba_model.elems):
            if e.crosec and e.id == '':
                for kcs, cs in zip(kcrosssecs, karamba_model.crosecs):
                    if cs == e.crosec and ke.elem_tag not in kcs.elem_tags:
                        if '' in kcs.elem_tags:
                            kcs.elem_tags.pop(kcs.elem_tags.index(''))
                        kcs.elem_tags.append(ke.elem_tag)
                for kmat, mat in zip(kmaterials, karamba_model.materials):
                    if mat == e.crosec.material and ke.elem_tag not in kmat.elem_tags:
                        if '' in kmat.elem_tags:
                            kmat.elem_tags.pop(kmat.elem_tags.index(''))
                        kmat.elem_tags.append(ke.elem_tag)
 
        ksupports = [KarambaSupport.from_karamba(s) for s in karamba_model.supports]
        kjoints = [KarambaJoint.from_karamba(j) for j in karamba_model.joints]

        if grounded_supports is not None:
            for k_supp in karamba_model.supports:
                for k_g in grounded_supports:
                    if k_supp.position.DistanceTo(k_g.position) < limit_distance:
                      knodes[k_supp.node_ind].is_grounded = True

        return cls(knodes, kelements, ksupports, kjoints, kmaterials, kcrosssecs,
            unit="meter", model_name=model_name)

    def to_karamba(self, loadcase=None, limit_dist=1e-6):
        mass = clr.Reference[float]()
        cog = clr.Reference[Point3]()
        warning_flag = clr.Reference[bool]()
        info = clr.Reference[str]()
        msg = clr.Reference[str]()
        kmodel = clr.Reference[Karamba.Models.Model]()

        in_points = List[Point3]([n.to_karamba(type='Point3') for n in self.nodes])
        in_elems = List[Karamba.Elements.BuilderElement]([e.to_karamba(type='builderbeam') for e in self.elements])
        for e in in_elems:
            etag =  None if e.id not in self.crosssec_id_from_etag else e.id
            cs_id = self.crosssec_id_from_etag[etag]
            e.crosec = self.crosssecs[cs_id].to_karamba()
        in_supports = List[Karamba.Supports.Support]([s.to_karamba() for _, s in self.supports.items()])

        in_crosecs = List[Karamba.CrossSections.CroSec]([cs.to_karamba() for cs in self.crosssecs])
        in_materials = List[Karamba.Materials.FemMaterial]([m.to_karamba() for m in self.materials])
        in_joints = List[Karamba.Joints.Joint]([j.to_karamba() for j in self.joints])

        if loadcase:
            kloadcase = KarambaLoadCase.from_loadcase(loadcase)
            in_loads = kloadcase.to_karamba()
        else:
            in_loads = [Karamba.Loads.Load()]

        Karamba.Models.AssembleModel.solve(
        	in_points,
        	in_elems,
        	in_supports,
        	in_loads,
        	in_crosecs,
        	in_materials,
        	List[Karamba.Elements.ElemSet]([]),
        	in_joints,
        	limit_dist,
            # out
        	kmodel, info, mass, cog, msg, warning_flag)

        if warning_flag.Value:
            print('info {} | msg {}'.format(info.Value, msg.Value))

        # # calculate Th.I response
        # max_disp = clr.Reference[List[float]]()
        # out_g = clr.Reference[List[float]]()
        # out_comp = clr.Reference[List[float]]()
        # message = clr.Reference[str]()
        # model = k3d.Algorithms.AnalyzeThI(model, max_disp, out_g, out_comp, message);

        # ucf = UnitsConversionFactories.Conv();
        # cm = ucf.cm();
        # print("max disp: {} {}".format(cm.toUnit(max_disp.Value[0]), cm.unitB))

        return kmodel.Value
