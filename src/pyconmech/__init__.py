"""
********************************************************************************
pyconmech
********************************************************************************

.. currentmodule:: pyconmech

This library provides python wrappers for efficient evaluation of construction mechanics.

.. toctree::
    :maxdepth: 3

    pyconmech.frame_analysis
    pyconmech.database

"""


from __future__ import print_function

import os
import sys
import decimal

# from .cpp_stiffness_checker import *
from .frame_analysis import StiffnessChecker
from .database import *

from .__version__ import __author__, __author_email__, __copyright__, __description__, __license__, __title__, __url__, __version__


__all__ = [
    '__author__', '__author_email__', '__copyright__', '__description__', 
    '__license__', '__title__', '__url__', '__version__', 
    'raise_if_windows',
    'raise_if_not_windows',
    'raise_if_ironpython',
    'raise_if_not_ironpython',
]


def is_windows():
    """Check if the operating system is Windows.
    Returns
    -------
    bool
        True if the OS is Windows. False otherwise
    """
    return os.name == 'nt'
WINDOWS = is_windows()


def is_ironpython():
    """Check if the Python implementation is IronPython.
    Returns
    -------
    bool
        True if the implementation is IronPython. False otherwise
    """
    return 'ironpython' in sys.version.lower()
IPY = is_ironpython()


def raise_if_not_windows():
    if not WINDOWS:
        raise


def raise_if_windows():
    if WINDOWS:
        raise


def raise_if_not_ironpython():
    if not IPY:
        raise


def raise_if_ironpython():
    if IPY:
        raise