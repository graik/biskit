## Setup version using setuptools to let autoinstallation with
## easy_install or pip install correctly treat the data files
## This setup.py has been adapted from the one shipped with django.
## It requires MANIFEST.in.

## building source distro : python setup.py sdist
## building windows distro: python setup.py bdist_wininst
## building rpm distro    : python setup.py bdist_rpm
## building egg distro    : python setup.py bdist_egg
## install source distro  : python setup.py install

## For custom installation folder use: --home=/where/i/want

from setuptools import setup
import os
import sys

# Compile the list of packages available, because distutils doesn't have
# an easy way to do this.
packages, data_files = [], []

root_dir = os.path.dirname(__file__)
# edit 26 May 2022: made a new variable rel_dir
rel_dir = os.path.relpath(__file__)
len_root_dir = len(root_dir)

biskit_dir = os.path.join(rel_dir, 'biskit')
doc_dir    = os.path.join(rel_dir, 'doc')
script_dir = os.path.join(rel_dir, 'scripts')

for dirpath, dirnames, filenames in os.walk(biskit_dir):
    # Ignore dirnames that start with '.'
    for i, dirname in enumerate(dirnames):
        if dirname.startswith('.'): del dirnames[i]

    if '__init__.py' in filenames:
        package = dirpath[len_root_dir:].lstrip('/').replace('/', '.')
        packages.append(package)
    else:
        data_files.append([dirpath, [os.path.join(dirpath, f) for f in filenames]])

## docs and scripts are moved from the root of the project into
## the package folder. That's why the separate treatment.
## First item in the data_files entry is the target folder, second item
## is the relative path in the svn project.
for dirpath, dirnames, filenames in os.walk( doc_dir ):
    # Ignore dirnames that start with '.'
    for i, dirname in enumerate(dirnames):
        if dirname.startswith('.'): del dirnames[i]
    data_files.append([os.path.join( 'biskit', dirpath ),
                       [os.path.join(dirpath, f) for f in filenames]])

for dirpath, dirnames, filenames in os.walk( script_dir ):
    # Ignore dirnames that start with '.'
    for i, dirname in enumerate(dirnames):
        if dirname.startswith('.'): del dirnames[i]
    data_files.append([os.path.join( 'biskit', dirpath ),
                       [os.path.join(dirpath, f) for f in filenames]])

# Small hack for working with bdist_wininst.
# See http://mail.python.org/pipermail/distutils-sig/2004-August/004134.html
if len(sys.argv) > 1 and sys.argv[1] == 'bdist_wininst':
    for file_info in data_files:
        file_info[0] = '/PURELIB/%s' % file_info[0]


long_description = \
 """Biskit is a modular, object-oriented Python library for structural
 bioinformatics research. It facilitates the manipulation and analysis
 of macromolecular structures, protein complexes, and molecular
 dynamics trajectories. For efficient number crunching, Biskit objects
 tightly integrate with numpy. Biskit also offers a platform for the
 integration of external programs and new algorithms into complex
 workflows. Calculations are thus often delegated to established
 programs like Xplor-NIH, Amber, Delphi, T-Coffee, TMAlign, Reduce and
 Modeller."""


setup(
    name = "biskit",
    version = "3.0.0",
    url = 'https://github.com/graik/biskit',
    download_url= 'https://github.com/graik/biskit/archive/refs/tags/v3.0.0.tar.gz',
    author = 'Raik Gruenberg, Johan Leckner and others',
    author_email = 'raik.gruenberg@crg.es',
    description = 'A Python platform for structural bioinformatics',
    long_description = long_description,
    provides=['biskit'],

    ## available on PyPi
    install_requires=['numpy', 'scipy', 'parmed'],
    packages = packages,
    include_package_data=True,
    data_files = data_files,
    scripts = ['scripts/bis.py'],

    classifiers= ['License :: OSI Approved :: GNU General Public License (GPL)',
                  'Topic :: Scientific/Engineering :: Bio-Informatics',
                  'Topic :: Scientific/Engineering :: Physics',
                  'Programming Language :: Python',
                  'Operating System :: OS Independent',
                  'Operating System :: POSIX',
                  'Operating System :: MacOS :: MacOS X',
                  'Intended Audience :: Science/Research',
                  'Development Status :: 5 - Production/Stable'
                  ]
)
