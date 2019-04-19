# conmech - a structural analysis computation engine for construction

**conmech** is an open-source library to provide efficient mechanics checkers for architectural construction sequencing. It's written in C++11 and wrapped friendly with Python via [pybind11].

## Prerequisites

**On Unix (Linux, OS X)**

* A compiler with C++11 support
* CMake >= 2.8.12
* Dependencies: [Eigen]

    ```bash
    sudo apt-get install libeigen3-dev
    ```

**On Windows**

Coming soon.

## Installation

Just clone this repository and pip install. Note the `--recursive` option which is needed for the pybind11 submodule:

```bash
git clone --recursive git@github.mit.edu:yijiangh/conmech.git
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
