
=========
Changelog
=========

All notable changes to this project will be documented in this file.

The format is based on `Keep a Changelog <https://keepachangelog.com/en/1.0.0/>`_
and this project adheres to `Semantic Versioning <https://semver.org/spec/v2.0.0.html>`_.

0.4.2
-----------

**Changed**

- Changed ``-DCONMECH_BUILD_TESTS=OFF`` in ``setup_cmake_utils.py`` to disable cpp test building in ``python setup.py build``

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
