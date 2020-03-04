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

.. image:: https://coveralls.io/repos/github/yijiangh/conmech/badge.svg?branch=master
    :target: https://coveralls.io/github/yijiangh/conmech?branch=master
    :alt: Coveralls

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

The following dependencies come from pybind11 for building the python wrappers.

**On Unix (Linux, OS X)**

* A compiler with C++11 support
* CMake >= 3.1

**On Windows**

* Visual Studio 2015 (required for all Python versions, see notes below)
* CMake >= 3.1

Then, clone this repository and pip install.

::

  cd conmech
  pip install .

With the ``setup.py`` file included in the base folder, the pip install command will invoke CMake and build the pybind11 module as specified in CMakeLists.txt.

**Note:**

*conmech*'s python bindings are built with a CMake-based build system via pybind11.
Take a look at `cmake_example for pybind11 <https://github.com/pybind/cmake_example>`_ 
if you want to learn more about this.

*conmech* depends on `Eigen <http://eigen.tuxfamily.org/index.php?title=Main_Page>`_ for linear algebra 
and `nlohmann::json <https://github.com/nlohmann/json>`_ 
for json (de-)serialization, both of which are handled automatically by cmake for downloading.

Build C++ code
^^^^^^^^^^^^^^

::

  mkdir build
  cd build
  cmake ..
  make -j2 # Unix

Or on Windows, replace the last line with

::

  cmake --build .
