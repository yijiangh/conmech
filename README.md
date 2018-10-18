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

## Test call

```python
import conmech_py as cm
json_path = <path_to_json/frame.json>
sc = cm.stiffness_checker(json_file_path = json_path, verbose = True)

# each row represents a nodal load
# entry 0 = node's id (specifed in the json fle),
# entry 1 - 7 = [Fx, Fy, Fz, Mx, My, Mz] (N, N-mm) in global coordinates.
ext_load = np.zeros([1,7])
ext_load[0,0] = 3
ext_load[0,3] = -500 * 1e3

sc.set_nodal_load(nodal_forces = ext_load, include_self_weight = False)
check_success = sc.solve()
print(check_success)

## partial structure check, the element ids here are compatible
## to the ones specified in the input json file
# existing_element_ids = [0,1]
# sc.solve(existing_element_ids)
```

[pybind11]: https://github.com/pybind/pybind11
[eigen]: http://eigen.tuxfamily.org/index.php?title=Main_Page
[BLAS]: https://www.netlib.org/blas/
[LAPACK]: http://www.netlib.org/lapack/
