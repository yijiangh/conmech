"""
********************************************************************************
pyconmech.frame_analysis
********************************************************************************

.. currentmodule:: pyconmech.frame_analysis

Package with functionality to compute elastic deformation of frame structures.

StiffnessChecker
-----------------

.. autosummary::
    :toctree: generated/
    :nosignatures:

    StiffnessChecker


Data IO
-------

.. autosummary::
    :toctree: generated/
    :nosignatures:

    Model
    LoadCase

Result IO
----------

Parsing / saving analsis results from json files 

.. autosummary::
    :toctree: generated/
    :nosignatures:

    read_frame_analysis_result_json

"""

from __future__ import absolute_import

from .stiffness_checker import *
from .io_base import *
# from .result_compare_utils import *
# from .numpy_stiffness import *

__all__ = [name for name in dir() if not name.startswith('_')]