## Note, this setup.py has been adapted from the one shipped with django
from distutils.core import setup
from distutils.command.install import INSTALL_SCHEMES
import os
import sys

# Tell distutils to put the data_files in platform-specific installation
# locations. See here for an explanation:
# http://groups.google.com/group/comp.lang.python/browse_thread/thread/35ec7b2fed36eaec/2105ee4d9e8042cb
for scheme in INSTALL_SCHEMES.values():
    scheme['data'] = scheme['purelib']


# Compile the list of packages available, because distutils doesn't have
# an easy way to do this.
packages, data_files = [], []

root_dir = os.path.dirname(__file__)
len_root_dir = len(root_dir)

biskit_dir = os.path.join(root_dir, 'Biskit')
doc_dir    = os.path.join(root_dir, 'docs')
script_dir = os.path.join(root_dir, 'scripts')

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
    data_files.append([os.path.join( 'Biskit', dirpath ),
                       [os.path.join(dirpath, f) for f in filenames]])

for dirpath, dirnames, filenames in os.walk( script_dir ):
    # Ignore dirnames that start with '.'
    for i, dirname in enumerate(dirnames):
        if dirname.startswith('.'): del dirnames[i]
    data_files.append([os.path.join( 'Biskit', dirpath ),
                       [os.path.join(dirpath, f) for f in filenames]])

# Small hack for working with bdist_wininst.
# See http://mail.python.org/pipermail/distutils-sig/2004-August/004134.html
if len(sys.argv) > 1 and sys.argv[1] == 'bdist_wininst':
    for file_info in data_files:
        file_info[0] = '/PURELIB/%s' % file_info[0]



setup(
    name = "Biskit",
    version = "2.3",
    url = 'http://biskit.pasteur.fr',
    download_url = 'https://sourceforge.net/project/platformdownload.php?group_id=150896',
    author = 'Raik Gruenberg & Johan Leckner',
    author_email = 'raik.gruenberg@crg.es',
    description = 'A Python platform for structural bioinformatics',
    requires=['numpy', 'ScientificPython', 'scipy', 'biopython' ],

    packages = packages,
    data_files = data_files,
##     data_files = {'Biskit' : ['Biskit/data'] },
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
