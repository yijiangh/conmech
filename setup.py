#!/usr/bin/env python
# -*- encoding: utf-8 -*-
from __future__ import absolute_import, print_function

import os
import re
import io

from setuptools import setup, find_packages
from setuptools.command.build_ext import build_ext
from setup_cmake_utils import CMakeExtension, CMakeBuild


here = os.path.abspath(os.path.dirname(__file__))


def read(*names, **kwargs):
    return io.open(
        os.path.join(here, *names),
        encoding=kwargs.get('encoding', 'utf8')
    ).read()


about = {}
exec(read('src', 'pyconmech', '__version__.py'), about)

requirements = read('requirements.txt').split('\n')

ext_modules = [
    CMakeExtension('_pystiffness_checker'),
    ]


setup(
    name=about['__title__'],
    version=about['__version__'],
    license=about['__license__'],
    description=about['__description__'],
    author=about['__author__'],
    author_email=about['__author_email__'],
    url=about['__url__'],
    long_description='%s\n%s' % (
        re.compile('^.. start-badges.*^.. end-badges', re.M |
                   re.S).sub('', read('README.rst')),
        re.sub(':[a-z]+:`~?(.*?)`', r'``\1``', read('CHANGELOG.rst'))
    ),
    packages=find_packages('src'),
    package_dir={'': 'src'},
    package_data={'pyconmech': ['data/*.json']},
    ext_modules=ext_modules,
    cmdclass=dict(build_ext=CMakeBuild),
    zip_safe=False,
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering',
        "License :: OSI Approved :: MIT License",
        'Operating System :: Unix',
        'Operating System :: Microsoft :: Windows',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
    ],
    keywords=['3D frame analysis', 'Finite Element Analysis', 'Structural Analysis'],
    install_requires=requirements,
    # extras_require={},
    # entry_points={},
)
