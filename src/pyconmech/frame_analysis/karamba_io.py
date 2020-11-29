import os, sys

assert 'ironpython' in sys.version.lower() and os.name == 'nt', 'Karamba only available for IronPython on Windows now.'

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(THIS_DIR))
from io_base import Node, Element, Support, Joint, CrossSec, Material, Model, LoadCase, \
    PointLoad, UniformlyDistLoad, GravityLoad

###############################################

# https://ironpython.net/documentation/dotnet/dotnet.html
import clr
clr.AddReferenceToFileAndPath("C:\Program Files\Rhino 6\Plug-ins\Karamba.gha")
clr.AddReferenceToFileAndPath("C:\Program Files\Rhino 6\Plug-ins\KarambaCommon.dll")

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
from Karamba.Utilities import MessageLogger, UnitsConversionFactories
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
    def from_karamba(cls, kelement, is_grounded=False):
        return cls(kelement.node_inds, kelement.ind, elem_tag=kelement.id, 
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
        return Karamba.Supports.Support(self.node_ind, List[bool](self.condition), Plane3())

class KarambaJoint(Joint):
    @classmethod
    def from_karamba(cls, kjoint):
        return cls(list(kjoint.c_condition(0)) + list(kjoint.c_condition(1)), kjoint.elemIds)

    def to_karamba(self):
        return Karamba.Supports.Support(self.node_ind, List[bool](self.condition), Plane3())

class KarambaCrossSec(CrossSec):
    @classmethod
    def from_karamba(cls, kc):
        return cls(kc.A, kc.Ipp, kc.Iyy, kc.Izz, elem_tags=kc.elemIds, family=kc.family, name=kc.name)

    def to_karamba(self):
        beam_mod = Karamba.CrossSections.CroSec_BeamModifier()
        beam_mod.A = self.A
        beam_mod.Ipp = self.Jx
        beam_mod.Iyy = self.Iy
        beam_mod.Izz = self.Iz
        beam_mod.elemIds = List[str]([tag for tag in self.elem_tags if tag is not None])
        beam_mod.family = self.family
        beam_mod.name = self.name
        return beam_mod

class KarambaMaterialIsotropic(Material):
    @classmethod
    def from_karamba(cls, km):
        assert km.typeName() == 'ISO', 'Input should be a Isotropic material!'
        return cls(km.E(), km.G12(), km.fy(), km.gamma(), elem_tags=km.elemIds, 
            family=km.family, name=km.name, type_name=km.typeName(), G3=km.G3())

    def to_karamba(self):
        steel_alphaT = 1e-5 #1/deg
        km = Karamba.Materials.FemMaterial_Isotrop(
	        self.family, self.name, self.E, self.G12, self.G3, self.density, self.fy, steel_alphaT, None) #Color.FromName("SlateBlue")
        for tag in self.elem_tags:
            if tag is not None:
                km.AddBeamId(tag)
        return km

class KarambaPointLoad(PointLoad):
    @classmethod
    def from_karamba(cls, kpl):
        pass

    def to_karamba(self, corotate=False):
        # corotate = true the pointload corotates with the node
        pl = Karamba.Loads.PointLoad(self.node_ind, Vector3(*self.force), Vector3(*self.moment), corotate)
        pl.loadcase = self.loadcase
        return pl

class KarambaUniformlyDistLoad(UniformlyDistLoad):
    @classmethod
    def from_karamba(cls, kpl):
        pass

    def to_karamba(self):
        # Karamba.Loads.LoadOrientation.global
        ul = Karamba.Loads.UniformlyDistLoad(List[str](self.elem_tags), Vector3(*self.q), 0, self.loadcase)
        return ul

class KarambaGravityLoad(GravityLoad):
    @classmethod
    def from_karamba(cls, kpl):
        pass

    def to_karamba(self):
        return Karamba.Loads.GravityLoad(Vector3(*self.force), self.loadcase)

####################################

class KarambaModel(Model):
    @classmethod
    def from_model(cls, model):
        knodes = [KarambaNode.from_data(n.to_data()) for n in model.nodes]
        kelements = [KarambaElement.from_data(e.to_data()) for e in model.elements]
        ksupports = [KarambaSupport.from_data(s.to_data()) for s in model.supports.values()]
        kjoints = [KarambaJoint.from_data(j.to_data()) for j in model.joints.values()]
        kmaterials = [KarambaMaterialIsotropic.from_data(m.to_data()) for m in model.materials.values()]
        kcrosssecs = [KarambaCrossSec.from_data(cs.to_data()) for cs in model.crosssecs.values()]
        return cls(knodes, kelements, ksupports, kjoints, kmaterials, kcrosssecs, 
            unit=model.unit, model_name=model.model_name)

    @classmethod
    def from_data(cls, data, verbose=False):
        model = super(KarambaModel, cls).from_data(data, verbose)
        return cls.from_model(model)

    @classmethod
    def from_karamba(cls, km):
        pass

    def to_karamba(self, loadcase, limit_dist=1e-6):
        mass = clr.Reference[float]()
        cog = clr.Reference[Point3]()
        warning_flag = clr.Reference[bool]()
        info = clr.Reference[str]()
        msg = clr.Reference[str]()
        kmodel = clr.Reference[Karamba.Models.Model]()

        in_points = List[Point3]([n.to_karamba(type='Point3') for n in self.nodes])
        in_elems = List[Karamba.Elements.BuilderElement]([e.to_karamba(type='builderbeam') for e in self.elements])
        in_supports = List[Karamba.Supports.Support]([s.to_karamba() for _, s in self.supports.items()])

        # in_crosecs = List[Karamba.CrossSections.CroSec]([])
        # in_materials = List[Karamba.Materials.FemMaterial]([])
        # in_joints = List[Karamba.Joints.Joint]([])

        in_crosecs = List[Karamba.CrossSections.CroSec]([cs.to_karamba() for _, cs in self.crosssecs.items()])
        in_materials = List[Karamba.Materials.FemMaterial]([m.to_karamba() for _, m in self.materials.items()])
        in_joints = List[Karamba.Joints.Joint]([j.to_karamba() for _, j in self.joints.items()])

        in_loads = List[Karamba.Loads.Load]()
        kpoint_loads = [KarambaPointLoad.from_data(pl.to_data()) for pl in loadcase.point_loads]
        kelem_loads = [KarambaUniformlyDistLoad.from_data(pl.to_data()) for pl in loadcase.uniform_element_loads]
        in_loads.AddRange([pl.to_karamba() for pl in kpoint_loads])
        in_loads.AddRange([pl.to_karamba() for pl in kelem_loads])
        if loadcase.gravity_load is not None:
            kgravity_load = KarambaGravityLoad.from_data(loadcase.gravity_load.to_data())
            in_loads.Add(kgravity_load)

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