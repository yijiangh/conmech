=======
conmech
=======

.. start-badges

.. image:: https://travis-ci.com/yijiangh/conmech.svg?branch=master
    :target: https://travis-ci.com/yijiangh/conmech
    :alt: Travis CI

.. image:: https://ci.appveyor.com/api/projects/status/k0f10bas2fj4uqww?svg=true
    :target: https://ci.appveyor.com/project/yijiangh/conmech
    :alt: Appveyor CI

.. image:: https://readthedocs.org/projects/conmech/badge/?version=latest
    :target: https://conmech.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

.. image:: https://img.shields.io/github/license/yijiangh/conmech
    :target: ./LICENSE
    :alt: License MIT

.. image:: https://img.shields.io/badge/python-3.6|3.7-blue
    :target: https://pypi.org/project/pyconmech/
    :alt: PyPI - Python Version

.. .. image:: https://img.shields.io/badge/pypi-v0.3.1-orange
    :target: https://pypi.org/project/pyconmech/
    :alt: PyPI - Latest Release

.. end-badges

.. Write project description

**conmech** is a stiffness checker that performs elastic deformation analysis for 3D frame structures. 
It is designed for construction sequencing applications, which involves testing
the partially assembled structure (subset of element permutation) many times.

**conmech** is written in C++11 and wrapped friendly with Python via `pybind11 <https://github.com/pybind/pybind11>`_.

Installation
------------

::

  pip install pyconmech


Build from source
-----------------

Build python bindings
^^^^^^^^^^^^^^^^^^^^^

Prerequisites
"""""""""""""

*conmech* depends on `Eigen <http://eigen.tuxfamily.org/index.php?title=Main_Page>`_ for linear algebra and `rapidjson <https://github.com/Tencent/rapidjson>`_ for json (de-)serialization, both of which are shipped with conmech.

The following dependencies come from pybind11 for building the python wrappers.

**On Unix (Linux, OS X)**

* A compiler with C++11 support
* CMake >= 2.8.12

**On Windows**

* Visual Studio 2015 (required for all Python versions, see notes below)
* CMake >= 3.1

*conmech*'s python bindings are built with a CMake-based build system via pybind11.
**It is recommended (especially for Windows users) to test the environment with the** `cmake_example for pybind11 <https://github.com/pybind/cmake_example>`_ **before proceeding to build conmech.**

Then, clone this repository and pip install. Note the ``--recursive`` option which is needed for cloning the submodules:

::

  git clone --recursive https://github.com/yijiangh/conmech
  pip install ./conmech
  # try with '--user' if you encountered a sudo problem

Or for developers:

::

  git clone --recursive https://github.com/yijiangh/conmech
  cd conmech
  python setup.py sdist
  pip install --verbose dist/*.tar.gz

With the ``setup.py`` file included in the base folder, the pip install command will invoke CMake and build the pybind11 module as specified in CMakeLists.txt.

Build C++ code
^^^^^^^^^^^^^^

::

  mkdir build
  cd build
  cmake ..
  make -j4 # Unix

Or on Windows, replace the last line with

::

  cmake --build .
