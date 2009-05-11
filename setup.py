## normal installation from source:     python setup.py install --exclude-scripts

## building source distro : python setup.py sdist

## building windows distro: python setup.py bdist_wininst

## building debian source package
## (requires stdeb: http://github.com/astraw/stdeb/tree/master)::
##
##   python -c "import stdeb; execfile('setup.py')" \
##          sdist_dsc --extra-cfg-file packaging/stdeb.cfg
##
##   cd deb_dist/
##
## building debian binary package
##   cd biskit-2.3/
##   cp ../../packaging/control debian/  ## fix control; edit if needed
##   cp ../../packaging/copyright debian/
##   dpkg-buildpackage -rfakeroot -uc -b
##   (this may require: sudo apt-get install python-all-dev)
##   (deb file located in deb_dist/)
##
## installing a local .deb file -- on Ubuntu:
##   gdebi python-biskit_2.3-1_all.deb
##
## installing a local .deb file -- other Debian:
## http://wiki.clug.org.za/wiki/How_do_I_install_a_.deb_file_I_downloaded_without_compromising_dependencies%3F

import sys

sys.path += ['Biskit']
from ez_setup import use_setuptools
use_setuptools()

from setuptools import setup, find_packages

from distutils.command.install import INSTALL_SCHEMES
import os

# Tell distutils to put the data_files in platform-specific installation
# locations. See here for an explanation:
# http://groups.google.com/group/comp.lang.python/browse_thread/thread/35ec7b2fed36eaec/2105ee4d9e8042cb
for scheme in INSTALL_SCHEMES.values():
    scheme['data'] = scheme['purelib']

try:
    root_dir = os.path.dirname(__file__)
except:
    root_dir = os.getcwd()

data_files = []
doc_dir    = os.path.join(root_dir, 'doc')
script_dir = os.path.join(root_dir, 'scripts')
## docs and scripts are moved from the root of the project into
## the package folder. That's why the separate treatment.
## First item in the data_files entry is the target folder, second item
## is the relative path in the svn project.
for dirpath, dirnames, filenames in os.walk( doc_dir ):
    # Ignore dirnames that start with '.'
    for i, dirname in enumerate(dirnames):
        if dirname.startswith('.'): del dirnames[i]
    data_files.append([os.path.join( 'Biskit', dirpath ),
                       [os.path.join(dirpath, f) for f in filenames]])

for dirpath, dirnames, filenames in os.walk( script_dir ):
    # Ignore dirnames that start with '.'
    for i, dirname in enumerate(dirnames):
        if dirname.startswith('.'): del dirnames[i]
    data_files.append([os.path.join( 'Biskit', dirpath ),
                       [os.path.join(dirpath, f) for f in filenames]])


setup(
    name = "biskit",
    version = "2.3",
    url = 'http://biskit.pasteur.fr',
    download_url = 'https://sourceforge.net/project/platformdownload.php?group_id=150896',
    author = 'Raik Gruenberg, Johan Leckner, and more',
    author_email = 'raik.gruenberg@crg.es',
    description = 'A Python platform for structural bioinformatics',

    ## available on PyPi: biopython
    install_requires=['biggles', 'scipy', 'biopython', 'ScientificPython',
                      'numpy' ],

    ## provide download links for packages not downloadable from PyPi
    dependency_links = \
    # biggles -- not available on PyPi
     ['http://biskit.pasteur.fr/install/mirror/biggles-1.6.5.tar.gz',

    # numpy, patched against issue #860 (overrides the PyPi download)
    # http://projects.scipy.org/numpy/ticket/860
      'http://biskit.pasteur.fr/install/mirror/numpy-1.3.0-patched.tar.gz',

    # ScientificPython -- registered on PyPi but download doesn't work
      'http://biskit.pasteur.fr/install/mirror/ScientificPython-2.8.tar.gz',

    # scipy -- registered on pypi but download doesn't work
      'http://biskit.pasteur.fr/install/mirror/scipy-0.7.0.tar.gz'],

    packages = find_packages(),
    include_package_data = True,
    data_files = data_files,


    scripts = ['scripts/Biskit/bispy'],

    classifiers= ['License :: OSI Approved :: GNU General Public License (GPL)',
                  'Topic :: Scientific/Engineering :: Bio-Informatics',
                  'Topic :: Scientific/Engineering :: Physics',
                  'Programming Language :: Python',
                  'Operating System :: POSIX',
                  'Operating System :: MacOS :: MacOS X',
                  'Intended Audience :: Science/Research',
                  'Development Status :: 5 - Production/Stable'
                  ]
)
