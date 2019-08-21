"""
********************************************************************************
pyconmech.frame_analysis
********************************************************************************

.. currentmodule:: pyconmech.frame_analysis

Package with functionality to compute elastic deformation of frame structures.

Main Classes
--------------

.. autosummary::
    :toctree: generated/
    :nosignatures:

    StiffnessChecker


File IO
-------

Parsing / saving frame structures from json files 

.. autosummary::
    :toctree: generated/
    :nosignatures:

    read_frame_json

Result comparison
-----------------

"""


from __future__ import absolute_import

from .stiffness_checker import *
from .frame_file_io import *
# from .frame import *

__all__ = [name for name in dir() if not name.startswith('_')]