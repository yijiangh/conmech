# conmech - an elastic analysis engine for 3D frame structures
[![Build Status](https://travis-ci.com/yijiangh/conmech.svg?branch=master)](https://travis-ci.com/yijiangh/conmech)
[![GitHub - License](https://img.shields.io/github/license/yijiangh/conmech)](./LICENSE)
[![PyPI - Python Version](https://img.shields.io/badge/python-2.5+|3.x-blue)](https://pypi.org/project/pyconmech/)
[![PyPI - Latest Release](https://img.shields.io/badge/pypi-v0.1.1-orange)](https://pypi.org/project/pyconmech/)
<!-- [![Build status](https://ci.appveyor.com/api/projects/status/k0f10bas2fj4uqww/branch/master?svg=true)](https://ci.appveyor.com/project/yijiangh/conmech/branch/master) -->

**conmech** is an open-source library to provide efficient stiffness checkers for architectural construction sequencing. It's written in C++11 and wrapped friendly with Python via [pybind11].

## Installation

```
pip install pyconmech
```

## Demo

(A cool gif should come here :satisfied:)

For examples of interactive usage in python (analysis for complete or partial structure in a construction sequence), see [stiffness_checker_test.ipynb](tests/notebook_demo/demo.ipynb).

## Build from source

### Build python bindings
#### Prerequisites
*conmech* depends on [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) for linear algebra and [rapidjson](https://github.com/Tencent/rapidjson) for json (de-)serialization, both of which are shipped with conmech.

The following dependencies come from [pybind11] for building the python wrappers.

**On Unix (Linux, OS X)**

* A compiler with C++11 support
* CMake >= 2.8.12

**On Windows**

* Visual Studio 2015 (required for all Python versions, see notes below)
* CMake >= 3.1

*conmech*'s python bindings are built with a CMake-based build system via [pybind11].
**It is recommended (especially for Windows users) to test the environment with the [cmake_example for pybind11](https://github.com/pybind/cmake_example) before proceeding to build conmech.**

Then, clone this repository and pip install. Note the `--recursive` option which is needed for cloning the submodules:

```bash
git clone --recursive https://github.com/yijiangh/conmech
pip install ./conmech
# try with '--user' if you encountered a sudo problem
```

Or for developers:

```bash
git clone --recursive https://github.com/yijiangh/conmech
cd conmech
python setup.py sdist
pip install --verbose dist/*.tar.gz
```

With the `setup.py` file included in the base folder, the pip install command will invoke CMake and build the pybind11 module as specified in CMakeLists.txt.

### Build C++ code

```bash
mkdir build
cd build
cmake ..
make -j4 # Unix
```
Or on Windows, replace the last line with
```
cmake --build .
```

[pybind11]: https://github.com/pybind/pybind11
[eigen]: http://eigen.tuxfamily.org/index.php?title=Main_Page
[BLAS]: https://www.netlib.org/blas/
[LAPACK]: http://www.netlib.org/lapack/
