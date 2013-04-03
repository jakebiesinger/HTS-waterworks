#!/usr/bin/env python

import sys, os
if not sys.version_info[0:2] >= (2,4):
    sys.stderr.write("Requires Python later than 2.4\n")
    sys.exit(1)

module_dependencies = open('requirements.txt').read().split()

from setuptools import setup
setup(
        name='HTS-waterworks',
        version='0.1',
        description='A python pipeline for analyzing high-throughput sequencing data',
        long_description=open('README.rst').read(),
        author='Jacob Biesinger',
        author_email='jake.biesinger@gmail.com',
        url='http://github.com/jakebiesinger/HTS-waterworks',
        install_requires = module_dependencies,
        setup_requires   = module_dependencies,
        classifiers=[
                    'Intended Audience :: End Users/Desktop',
                    'Intended Audience :: Developers',
                    'Intended Audience :: Science/Research',
                    'Intended Audience :: Information Technology',
                    'License :: OSI Approved :: MIT License',
                    'Programming Language :: Python',
                    'Topic :: Scientific/Engineering',
                    'Topic :: Scientific/Engineering :: Bio-Informatics',
                    'Topic :: System :: Distributed Computing',
                    'Topic :: Software Development :: Build Tools',
                    'Topic :: Software Development :: Build Tools',
                    'Topic :: Software Development :: Libraries',
                    'Environment :: Console',
                    ],
        license = "MIT",
        keywords = "pipeline parallel bioinformatics science",

        packages=['hts_waterworks'],
        package_dir={'hts_waterworks': 'hts_waterworks'},
        include_package_data = True,    # include everything in source control
     )

#
#  http://pypi.python.org/pypi
#  http://docs.python.org/distutils/packageindex.html
#   
# 
# 
# python setup.py register
# python setup.py sdist --format=gztar upload
