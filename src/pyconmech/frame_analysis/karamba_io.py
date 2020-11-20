import os, sys
THIS_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(THIS_DIR))
from io_base import Node, Element, Support, Joint

###############################################

# https://ironpython.net/documentation/dotnet/dotnet.html
import clr
clr.AddReferenceToFileAndPath("C:\Program Files\Rhino 6\Plug-ins\Karamba.gha")
clr.AddReferenceToFileAndPath("C:\Program Files\Rhino 6\Plug-ins\KarambaCommon.dll")

import feb # Karamba's C++ library (undocumented in the API)
import Karamba
import Karamba.Models.Model as KModel
import Karamba.Nodes
import Karamba.Elements
import Karamba.Elements.ModelTruss as KTruss
from Karamba.Geometry import Point3, Line3, Vector3, Plane3
import Karamba.CrossSections.CroSec as KCroSec
import Karamba.Supports.Support as KSupport
import Karamba.Loads.Load as KLoad
from Karamba.Utilities import MessageLogger, UnitsConversionFactories
import KarambaCommon
#
import System
from System import GC
from System.Collections.Generic import List

####################################

class KarambaNode(Node):
    @classmethod
    def from_karamba_node(cls, knode, is_grounded=False):
        point = [knode.pos.X, knode.pos.Y, knode.pos.Z]
        return cls(point, knode.ind, is_grounded)
    
    def to_karamba_node(self):
        return Karamba.Nodes.Node(self.node_ind, Point3(*self.point))

class KarambaElement(Element):
    @classmethod
    def from_karamba_element(cls, kelement, is_grounded=False):
        return cls(kelement.node_inds, kelement.ind, elem_tag=kelement.id, 
            bending_stiff=kelement.bending_stiff)

    def to_karamba_builder_beam(self):
        return Karamba.Elements.BuilderBeam(*self.end_node_inds)

class KarambaSupport(Support):
    @classmethod
    def from_karamba_support(cls, ksupport):
        return cls(list(ksupport.Condition), ksupport.node_ind)

    def to_karamba_support(self):
        return Karamba.Supports.Support(self.node_ind, List[bool](self.condition), Plane3())

class KarambaJoint(Joint):
    @classmethod
    def from_karamba_joint(cls, kjoint):
        return cls(list(kjoint.c_condition(0)) + list(kjoint.c_condition(1)), kjoint.elemIds)

    def to_karamba_support(self):
        return Karamba.Supports.Support(self.node_ind, List[bool](self.condition), Plane3())



####################################

# def model_from_karamba_model(karamba_model):
#     model = Model(nodes, elements, supports, joints, 
#         materials, crosssecs, unit='meter', model_name=None)