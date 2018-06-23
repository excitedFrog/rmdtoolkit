# Python 3.6.1

from setuptools import setup
from setuptools import find_packages
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

# Arguments marked as "Required" below must be included for upload to PyPI.
# Fields marked as "Optional" may be commented out.

setup(name='rmdtoolkit',
      version='0.9.0',
      description='A toolkit for LAMMPS and RAPTOR.',
      long_description=long_description,
      url='https://github.com/excitedFrog/rmdtoolkit',
      author='Jeff Li',
      classifiers=['Development Status :: 3 - Alpha',
                   'Intended Audience :: Science/Research',
                   'Topic :: Scientific/Engineering :: Chemistry',
                   'License :: OSI Approved :: MIT License',
                   'Programming Language :: Python :: 3'],
      keywords='chemistry, molecular_dynamics',
      packages=find_packages(),
      install_requires=['numpy', 'scipy', 'mendeleev', 'pandas', 'twilio', 'matplotlib'],  # TODO:Take care of this.
      extras_require={  # Optional
          'dev': ['check-manifest'],
          'test': ['coverage'],
      },
      package_data={  # Optional
          'sample': ['package_data.dat'],
      })
