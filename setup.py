#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages
from pysashimi.version import __version__

setup(name='pysashimiplot',
      version=__version__,
      python_requires='>3.6',
      description='sashimiplot',
      author='Ran zhou',
      author_email='ranzhou1005@gmail.com',
      url='https://github.com/luguodexxx/pysashimi',
      # dependency_links=['https://pypi.douban.com/simple/'],
      license='MIT',
      classifiers=[
          'Development Status :: 5 - Production/Stable',
          'Intended Audience :: Science/Research',
          'Topic :: Scientific/Engineering :: Bio-Informatics',
          'License :: OSI Approved :: MIT License',
          'Programming Language :: Python :: 3',
      ],
      keywords='sashimiplot',
      packages=find_packages(),
      install_requires=[
          'xmltodict',
          'biothings_client',
          'pysam',
          'requests',
          'click',
          'numpy',
          'matplotlib',
          'loguru',
          'seaborn',
          'scipy'
      ],
      entry_points={
          'console_scripts': [
              'sashimiplot=pysashimi.runplot:cli'
          ],
      },
      )
