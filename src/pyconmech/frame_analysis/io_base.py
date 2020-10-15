import datetime

class Model(object):
    def __init__(self, nodes, elements, supports, joints, materials, crosssecs, model_name=None):
        self.model_name = model_name
        self.generate_time = str(datetime.datetime.now())
        self.nodes = nodes
        self.elements = elements
        self.supports = supports
        self.joints = joints
        self.crosssecs = crosssecs
        self.materials = materials

    @property
    def node_num(self):
        return len(self.nodes)

    @property
    def element_num(self):
        return len(self.elements)

    @classmethod
    def from_data(cls, data):
        from pyconmech.frame_analysis.frame_file_io import read_frame_data
        nodes, elements, supports, joints, materials, crosssecs, model_name, unit = read_frame_data(data)
        return cls(nodes, elements, supports, joints, materials, crosssecs, model_name=model_name)

    def to_data(self):
        from pyconmech.frame_analysis.frame_file_io import frame_to_data
        return frame_to_data(self.nodes, self.elements, self.supports, self.joints, self.materials, self.crosssecs, self.model_name)


class Node(object):
    def __init__(self, point, node_ind, is_grounded):
        self.point = point
        self.node_ind = node_ind
        self.is_grounded = is_grounded

    @classmethod
    def from_data(cls, data):
        return cls(data['point'], data['node_ind'], data['is_grounded'])

    def to_data(self):
        data = {'point' : self.point, 'node_ind' : self.node_ind, 'is_grounded' : self.is_grounded}
        return data

    def __repr__(self):
        return '{}(#{},{},Grd:{})'.format(self.__class__.__name__, self.node_ind, self.point, self.is_grounded)

class Element(object):
    def __init__(self, end_node_inds, elem_ind, elem_tag='', bending_stiff=True):
        assert end_node_inds[0] != end_node_inds[1], 'zero length element not allowed!'
        self.end_node_inds = end_node_inds
        self.elem_tag = elem_tag
        self.elem_ind = elem_ind
        self.bending_stiff = bending_stiff

    @classmethod
    def from_data(cls, data):
        return cls(data['end_node_inds'], data['elem_ind'], data['elem_tag'], data['bending_stiff'])

    def to_data(self):
        raise NotImplementedError()

    def __repr__(self):
        return '{}(#{}({}),{},Bend:{})'.format(self.__class__.__name__, self.elem_ind, self.elem_tag, self.end_node_inds, self.bending_stiff)

class Support(object):
    def __init__(self, condition, node_ind):
        self.condition = condition
        self.node_ind = node_ind

    @classmethod
    def from_data(cls, data):
        return cls(data['condition'], data['node_ind'])

    def to_data(self):
        raise NotImplementedError()

    def __repr__(self):
        return '{}(#{},{})'.format(self.__class__.__name__, self.node_ind, self.condition)

class Joint(object):
    def __init__(self, c_conditions, elem_tags):
        self.c_conditions = c_conditions
        self.elem_tags = elem_tags

    @classmethod
    def from_data(cls, data):
        return cls(data['c_conditions'], data['elem_tags'])

    def to_data(self):
        raise NotImplementedError()

# TODO: assumed to be converted to meter-based unit
class CrossSec(object):
    def __init__(self, A, Jx, Iy, Iz, elem_tags=None, family='unnamed', name='unnamed'):
        self.A = A
        self.Jx = Jx
        self.Iy = Iy
        self.Iz = Iz
        self.elem_tags = elem_tags if elem_tags else [None]
        self.family = family
        self.name = name

    @classmethod
    def from_data(cls, data):
        return cls(data['A'], data['Jx'], data['Iy'], data['Iz'], data['elem_tags'], data['family'], data['name'])

    def to_data(self):
        raise NotImplementedError()

    def __repr__(self):
        return 'family:{} name:{} area:{}[m2] Jx:{}[m4] Iy:{}[m4] Iz:{}[m4] applies to elements:{}'.format(
            self.family, self.name, self.A, self.Jx, self.Iy, self.Iz, self.elem_tags)

def mu2G(mu, E):
    """compute shear modulus from poisson ratio mu and Youngs modulus E
    """
    return E/(2*(1+mu))

def G2mu(G, E):
    """compute poisson ratio from shear modulus and Youngs modulus E
    """
    return E/(2*G)-1

# TODO: assumed to be converted to kN/m2, kN/m3
class Material(object):
    def __init__(self, E, G12, fy, density, elem_tags=None, family='unnamed', name='unnamed', type_name='ISO'):
        self.E = E
        self.G12 = G12
        self.mu = G2mu(G12, E)
        self.fy = fy
        self.density = density
        self.elem_tags = elem_tags if elem_tags else [None]
        self.family = family
        self.name = name
        self.type_name = type_name

    @classmethod
    def from_data(cls, data):
        return cls(data['E'], data['G12'], data['fy'], data['density'], data['elem_tags'], data['family'], data['name'], data['type_name'])

    def to_data(self):
        raise NotImplementedError()

    def __repr__(self):
        # G3:8076[kN/cm2]
        return 'Material: |{}| E:{}[kN/m2] mu:{} G12:{}[kN/m2] density:{}[kN/m3] fy:{}[kN/m2] applies to elements:{}'.format(
            self.family+'-'+self.name, self.E, self.mu, self.G12, self.density, self.fy, self.elem_tags)

class PointLoad(object):
    def __init__(self, force, moment, node_ind):
        self.force = force
        self.moment = moment
        self.node_ind = node_ind

    @classmethod
    def from_data(cls, data):
        return cls(data['force'], data['moment'], data['node_ind'])

    def to_data(self):
        raise NotImplementedError()

    def __repr__(self):
        return 'Pointload: node_ind {} | force {} | moment {}'.format(self.node_ind, self.force, self.moment)

class UniformlyDistLoad(object):
    def __init__(self, q, load, elem_tags):
        self.q = q
        # not sure if load is used in Karamba, we only use q here
        self.load = load
        self.elem_tags = elem_tags

    @classmethod
    def from_data(cls, data):
        return cls(data['q'], data['load'], data['elem_tags'])

    def to_data(self):
        raise NotImplementedError()

    def __repr__(self):
        return 'UniformlyDistLoad: element_tags {} | q {} | load {}'.format(self.elem_tags, self.q, self.load)

class GravityLoad(object):
    def __init__(self, force=[0,0,-1]):
        self.force = force

    @classmethod
    def from_data(cls, data):
        return cls(data['force'])

    def to_data(self):
        raise NotImplementedError()

    def __repr__(self):
        return 'GravityLoad: {}'.format(self.force)