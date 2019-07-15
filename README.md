# conmech - an elastic analysis engine for 3D frame structures

**conmech** is an open-source library to provide efficient stiffness checkers for architectural construction sequencing. It's written in C++11 and wrapped friendly with Python via [pybind11].

## Prerequisites
*conmech* depends on [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) for linear algebra and [rapidjson](https://github.com/Tencent/rapidjson) for json (de-)serialization, both of which are shipped with conmech.

The following dependencies come from [pybind11] for building the python wrappers.

**On Unix (Linux, OS X)**

* A compiler with C++11 support
* CMake >= 2.8.12

**On Windows**

* Visual Studio 2015 (required for all Python versions, see notes below)
* CMake >= 3.1

## Installation

*conmech*'s python bindings are built with a CMake-based build system via [pybind11].
**It is recommended (especially for Windows users) to test the environment with the [cmake_example for pybind11](https://github.com/pybind/cmake_example) before proceeding to build conmech.**


Just clone this repository and pip install. Note the `--recursive` option which is needed for cloning the submodules:

```bash
git clone --recursive https://github.com/yijiangh/conmech
pip install ./conmech --user
```

With the `setup.py` file included in the base folder, the pip install command will invoke CMake and build the pybind11 module as specified in CMakeLists.txt.

## Demo

(A cool gif should come here :satisfied:)

For examples of interactive usage in python (analysis for complete or partial structure in a construction sequence), see [stiffness_checker_test.ipynb](src/bindings/pyconmech/test/stiffness_checker_test.ipynb).

[pybind11]: https://github.com/pybind/pybind11
[eigen]: http://eigen.tuxfamily.org/index.php?title=Main_Page
[BLAS]: https://www.netlib.org/blas/
[LAPACK]: http://www.netlib.org/lapack/
