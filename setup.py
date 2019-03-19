from setuptools import setup
import sys
import os
import re

PACKAGENAME = 'sjoert'
packageDir = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          PACKAGENAME)

# Obtain the package version
# versionFile = os.path.join(packageDir, 'version.py')
# with open(versionFile, 'r') as f:
#           s = f.read()
# # Look up the string value assigned to __version__ in version.py using regexp
# versionRegExp = re.compile("__VERSION__ = \"(.*?)\"")
# # Assign to __version__
# __version__ =  versionRegExp.findall(s)[0]
# print(__version__)
__version__ = 1.0

# create requirements file
# setupDir = os.path.join(packageDir, '..', 'setup')
# genRequirements = os.path.join(setupDir, 'generate_requirements.py')
# print(genRequirements)


setup(# package information
      name=PACKAGENAME,
      version=__version__,
      description='collection of tools for astronomy',
      long_description=''' ''',
      # What code to include as packages
      packages=[PACKAGENAME],
      packagedir={PACKAGENAME: 'sjoert'},
      # What data to include as packages
      include_package_data=False
      # package_data={PACKAGENAME:['data/*.txt',
      #                            'filters/P48/*.dat',
      #                            'filters/SEDm/*.dat',
      #                            'filters/UVOT/*.txt',]}
      )