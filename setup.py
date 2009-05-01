## normal installation from source:     python setup.py install
##
## building source distro : python setup.py sdist
## building windows distro: python setup.py bdist_wininst
##
## building debian source package
## (requires stdeb: http://github.com/astraw/stdeb/tree/master)::
##
##   python -c "import stdeb; execfile('setup.py')" \
##          sdist_dsc --extra-cfg-file stdeb.cfg
##
##   cd deb_dist/
##
##   zip biskit-2.3-src.deb.zip *gz *dsc
##
##   You may need to edit deb_dist/biskit-2.3/debian/control
##

## building debian binary package
##   unzip biskit-2.3-src.deb.zip
##   dpkg-source -x biskit_2.3-1.dsc
##   cd biskit-2.3/
##   dpkg-buildpackage -rfakeroot -uc -b
##   (this may require: sudo apt-get install python-setuptools python-all-dev)
##
## installing a local .deb file -- on Ubuntu:
## gdebi python-biskit_2.3-1_all.deb
##
## installing a local .deb file -- other Debian:
## http://wiki.clug.org.za/wiki/How_do_I_install_a_.deb_file_I_downloaded_without_compromising_dependencies%3F

import sys

sys.path += ['Biskit']
from ez_setup import use_setuptools
use_setuptools()

from setuptools import setup, find_packages


setup(
    name = "biskit",
    version = "2.3",
    url = 'http://biskit.pasteur.fr',
    download_url = 'https://sourceforge.net/project/platformdownload.php?group_id=150896',
    author = 'Raik Gruenberg & Johan Leckner',
    author_email = 'raik.gruenberg@crg.es',
    description = 'A Python platform for structural bioinformatics',

    ## available on PyPi
    install_requires=['numpy', 'ScientificPython', 'biopython', 'scipy',
                      'biggles' ],

    ## not available on PyPi
    dependency_links = \
    # biggles
      ['http://downloads.sourceforge.net/biggles/python2-biggles-1.6.5.tar.gz',
    # scipy -- is registered on pypi but the download doesn't work
      'https://sourceforge.net/project/showfiles.php?group_id=27747'],

    packages = find_packages(),
    include_package_data = True,

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
