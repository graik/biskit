biskit
==========
___a software platform for structural bioinformatics___

[![Build Status](https://travis-ci.org/graik/biskit.svg?branch=biskit3)](https://travis-ci.org/graik/biskit)
[![Coverage Status](https://coveralls.io/repos/github/graik/biskit/badge.svg?branch=biskit3&service=github)](https://coveralls.io/github/graik/biskit?branch=biskit3)

Please refer to 
            **http://biskit.pasteur.fr**
for installation and usage instructions, troubleshooting and
everything else!

Biskit is a modular, object-oriented python library that provides
intuitive classes for many typical tasks of structural bioinformatics
research. It facilitates the manipulation and analysis of
macromolecular structures, protein complexes, and molecular dynamics
trajectories. At the same time, Biskit offers a software platform for
the rapid integration of external programs and new algorithms into
complex structural bioinformatics workflows. Calculations are thus
often delegated to established programs like Xplor, Amber, Hex, Prosa,
and DelPhi; interfaces to further software can be added
easily. 

Python 3 Migration
-------------------

This is the Python 3 branch of Biskit. Migration is ongoing but the core functionalities and many of the wrappers have been migrated. The new version is found in (lower case) ```biskit```. All the modules that have **not** yet been migrated are in the old ```Biskit2``` folder. After installation, the new Python 3 Biskit is imported with ```import biskit```.

Short Installation Instructions
--------------------------------

___1. Install (plotting) dependencies___

The biskit library itself does not strictly need these and you can also safely install them later. However, biggles (https://biggles-plot.github.io/) is an excellent plotting library with a relatively intuitive syntax that is used throughout biskit and several unittests depend on it. gnuplot is wrapped by `biskit.gnuplot` and offers super-convenient quick and dirty line, scatter and histogram plots for rapid interactive data inspection. 

On Debian / Ubuntu:

  *  ```sh
     sudo apt-get install libplot-dev plotutils  ## needed for biggles compilation
     sudo apt-get install gnuplot ## program required by biskit.gnuplot
     ```

On Mac OS-X:

  * install Quartz (https://www.xquartz.org/)
  *  ```sh
     brew install plotutils --with-x11  # using homebrew
     ```

Then (all systems):

  *  ```
     pip install biggles
     ```

___2. Install biskit___

```sh
git clone https://github.com/graik/biskit.git biskit -b biskit3
pip install -r biskit/requirements.txt
pip install -e biskit
```
If not already available, this will also install numpy, scipy, and BioPython. Replace `git clone` by the appropriate `tar xvf *tgz` command to start from an official Biskit release bundle.

License
-------

Biskit is distributed under the GNU GPL version 3. See LICENSE.txt. Contact us, in case you prefer a different licensing model.
