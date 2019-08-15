"""
********************************************************************************
pyconmech.frame_analysis
********************************************************************************

.. currentmodule:: pyconmech.frame_analysis

Package with functionality to compute elastic deformation of frame structures.

frame_analysis
--------------

.. autosummary::
    :toctree: generated/

    stiffness_checker

"""

from __future__ import absolute_import

from .stiffness_checker import *
from .frame_file_io import *
# from .frame import *

__all__ = [name for name in dir() if not name.startswith('_')]