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


def get_packages(target):
    packages = [
        dirpath.replace("/", ".").lstrip(".")
        for dirpath, _, filenames in os.walk(os.path.join(".", target))
        if "__init__.py" in filenames
    ]
    return [package for package in packages if len(package) > 0]


def get_files(target, prefix):
    data_files = [
        [
            os.path.join(prefix, dirpath),
            [os.path.join(dirpath, f) for f in filenames if "__init__.py" not in f],
        ]
        for dirpath, _, filenames in os.walk(target)
    ]
    return data_files


def get_data_files_and_packages():
    root_dir = os.path.dirname(__file__)

    biskit_dir = "biskit"
    doc_dir = "doc"
    script_dir = "scripts"

    cwd = os.getcwd()
    os.chdir(root_dir)
    packages = get_packages(biskit_dir)
    print(packages)
    data_files = get_files(biskit_dir, ".")
    data_files += get_files(doc_dir, biskit_dir)
    data_files += get_files(script_dir, biskit_dir)
    os.chdir(cwd)
    return packages, data_files


packages, data_files = get_data_files_and_packages()

# Small hack for working with bdist_wininst.
# See http://mail.python.org/pipermail/distutils-sig/2004-August/004134.html
if len(sys.argv) > 1 and sys.argv[1] == "bdist_wininst":
    for file_info in data_files:
        file_info[0] = "/PURELIB/%s" % file_info[0]


long_description = \
""" Biskit is a modular, object-oriented Python library for structural
 bioinformatics research. It facilitates the manipulation and analysis
 of macromolecular structures, protein complexes, and molecular
 dynamics trajectories. For efficient number crunching, Biskit objects
 tightly integrate with numpy. Biskit also offers a platform for the
 rapid integration of external programs and new algorithms into complex
 workflows. Calculations are thus often delegated to established
 programs like Xplor, Amber, Hex, Prosa, T-Coffee, TMAlign, Reduce and
 Modeller."""


setup(
    name="biskit",
    version="3.0.0a0",
    url="http://biskit.pasteur.fr",
    download_url="https://github.com/graik/biskit/archive/v3.0.0.tar.gz",
    author="Raik Gruenberg, Johan Leckner and others",
    author_email="raik.gruenberg@crg.es",
    description="A Python platform for structural bioinformatics",
    long_description=long_description,
    provides=["biskit"],
    ## available on PyPi
    install_requires=["numpy", "scipy", "biopython"],
    packages=packages,
    include_package_data=True,
    data_files=data_files,
    scripts=["scripts/bis.py"],
    classifiers=[
        "License :: OSI Approved :: GNU General Public License (GPL)",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Physics",
        "Programming Language :: Python",
        "Operating System :: OS Independent",
        "Operating System :: POSIX",
        "Operating System :: MacOS :: MacOS X",
        "Intended Audience :: Science/Research",
        "Development Status :: 5 - Production/Stable",
    ],
)
