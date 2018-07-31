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

```sh
git clone https://github.com/graik/biskit.git biskit -b biskit3
pip3 install -r biskit/requirements.txt
pip3 install -e biskit
```
Replace `git clone` by the appropriate `tar xvf *tgz` command to start from an official Biskit release bundle.

If not already available, this will install numpy, scipy, and BioPython. It will *not* however install the biggles plotting package, which is not strictly required but highly recommended.

___Install plotting libraries___

(1) **biggles** -- the main plotting library used by biskit, several unittests depend on it
 
  on Debian / Ubuntu:

   ```sh
   sudo apt-get install libplot-dev plotutils
   pip install biggles
   ```
 
  on Mac OSX:

   * install Quartz (https://www.xquartz.org/)
    
   ```sh
   brew install plotutils --with-x11 # homebrew
   pip install biggles
   ```

(2) **gnuplot** -- the `biskit.gnuplot` wrapper offers no-frills quick and dirty line and scatter plots (`plot()`, `scatter()`); 
  especially useful for rapid interactive inspection of data
  
  on Debian / Ubuntu: `sudo apt-get install gnuplot`
  
  on Mac OS-X: `brew install gnuplot --with-X11` (requires Quartz)

License
-------

Biskit is distributed under the GNU GPL version 3. See LICENSE.txt. Contact us, in case you prefer a different licensing model.
