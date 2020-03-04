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

.. .. image:: https://img.shields.io/badge/pypi-v0.1.1-orange
    :target: https://pypi.org/project/pyconmech/
    :alt: PyPI - Latest Release
.. end-badges

.. Write project description

**conmech** is a stiffness checker that performs elastic deformation analysis for 3D frame structures. 
It is designed for construction sequencing applications, which involves testing
the partially assembled structure (subset of element permutation) many times.

**conmech** is written in C++11 and wrapped friendly with Python via `pybind11 <https://github.com/pybind/pybind11>`_.

========
Contents
========

.. toctree::
   :maxdepth: 1

   readme
   under_the_hood/index
   reference
   developer_notes/index
   contributing
   authors
   changelog

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
