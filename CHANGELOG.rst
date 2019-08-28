
=========
Changelog
=========

.. # with overline, for parts
.. * with overline, for chapters
.. =, for sections
.. -, for subsections
.. ^, for subsubsections
.. ", for paragraphs

0.1.2
-----

* Initial version

0.1.3
-----

Changed
^^^^^^^

- The original ``stiffness_checker`` extension module is wrapper as ``_stiffness_checker``.
  All the cpp modules are wrapper under a top-level python classes/functions, to give more
  flexibility.
- **API change**: ``stiffness_checker`` class is renamed to ``StiffnessChecker`` to conform
    to the class naming convention. All other APIs within this class are left unchanged.
- Delete ``radius`` entry from ``material_properties``.


Added
^^^^^

- documentation is hosted on readthedocs!
- add grasshopper examples - parse/save files, karamba comparsion, solve/get result in GH via ghpython-remote
- supports material / cross sectional properties for each element. 
- supports uniformly distributed load
