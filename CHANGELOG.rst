
=========
Changelog
=========

All notable changes to this project will be documented in this file.

The format is based on `Keep a Changelog <https://keepachangelog.com/en/1.0.0/>`_
and this project adheres to `Semantic Versioning <https://semver.org/spec/v2.0.0.html>`_.

Unreleased
----------

**added**

* An example run script is added to the `examples`
* Added a direct call to `StiffnessChecker` in GHPython, via compas RPC proxy.
* Added `utils` module, providing function call access to stiffness checker computation.
* Added `StiffnessChecker.get_sigma_max_per_element`, only support rounded cross sec.

**changed**

* Changed using the default element tag `''` to assign a default, distinctive element tag for each element.
* Changed how we parse karamba model's crosssecs and materials if they are specified through element but directly as an input to the assemble component.
* Changed `stiffness_base`'s' `get_element_crosssec` and `get_element_material` function to raise error if e_tag is not assigned.
* The element reactions' sign is flipped, so equilibrium is expressed as `element reaction + external force + support reaction = 0`

**fixed**

- `karamba_io`'s import path problem fixed, try two default options, and raise error if neither works

0.6.0
----------

**added**

- Added `karamba_io` module to facilitate `Karamba3d`-`JSON` two way conversion.
- Added `Model` and `LoadCase` classes
- Added element information to `ValueError` when encountering zero-length line element
- Added print out attributes to `Node` and `Element`

**changed**

- Changed `StiffnessBase`'s `__init__` from inputting `nodes, elements, ...` to `Model`
- Moved `get_element_crossec` and `get_element_material` from `NumpyStiffness` to `StiffnessBase`
- Changed `StiffnessChecker`'s `set_loads` to take `LoadCase` object
- Changed `Joint`'s `c_conditions` attribute to a 12-entry list instead of a dictionary with two six-entry lists
- Karamba conversion is now done more formally with `karam_io` modules, without using an ad-hoc C# component.

**removed**

**fixed**


0.5.0
----------

**added**

- added the `numpy` engine
- added removing loads functionality by passing in `none`
- added base classes for input data: `io_base.node, element, support, joint, material, crosssec, pointload, etc.`

**changed**

- changed frame data format, use karamba exported format
- changed `stiffnesschecker`'s `__init__` function's arguments to take `io_base.*` data
- changed `stiffnesschecker`'s default behavior: not applying gravity load by default
- changed ``-dconmech_build_tests=off`` in ``setup_cmake_utils.py`` to disable cpp test building in ``python setup.py build``

**removed**

**fixed**

**changed**

0.4.0
-----------

**Added**

- Added ``StiffnessChecker`` class method, directly construct from frame data, without saving data to a temp json
- Added some initial cpp unit tests, test data fed in by CMake and tests organized by ``Catch2``

**Changed**

- Changed ``rapidjson`` to ``nlohmann::json``

**Removed**

- Removed the `Frame` data structure in Stiffness checker's cpp backend
- Removed all the git submodule and used CMake download external instead

**Fixed**

- Fixed the memory leak caused by the smart pointer cycle dependency in ``Frame``

0.3.1
----------

**Added**

- Added unit tests for `std::throw` in parsing material properties

0.3.0
----------

**Changed**

- Changed `try/catch` in the C++ file parsing to `std::throw` 

0.2.0
-----

**Changed**

- The original ``stiffness_checker`` extension module is wrapper as ``_stiffness_checker``.
  All the cpp modules are wrapper under a top-level python classes/functions, to give more
  flexibility.
- **API change**: ``stiffness_checker`` class is renamed to ``StiffnessChecker`` to conform
    to the class naming convention. All other APIs within this class are left unchanged.
- Delete ``radius`` entry from ``material_properties``.


**Added**

- documentation is hosted on readthedocs!
- add grasshopper examples - parse/save files, karamba comparsion, solve/get result in GH via ghpython-remote
- supports material / cross sectional properties for each element. 
- supports uniformly distributed load
- add gravity magnitude and direction

0.1.0
-----

Initial version
