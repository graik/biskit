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

The biskit library itself does not strictly need these and you can also safely install them later. However, **biggles** (https://biggles-plot.github.io/) is an excellent plotting library with a relatively intuitive syntax that is used throughout biskit and several unittests depend on it. **gnuplot** is wrapped by `biskit.gnuplot` and offers super-convenient quick and dirty line, scatter and histogram plots for rapid interactive data inspection. 

On Debian / Ubuntu:
  *  ```sh
     sudo apt-get install libplot-dev plotutils  ## needed for biggles compilation
     sudo apt-get install gnuplot ## program required by biskit.gnuplot
     pip install biggles
     ```

On Mac OS-X:
  * install Quartz (https://www.xquartz.org/)
  *  ```sh
     brew install plotutils --with-x11
     brew install gnuplot --with-x11
     pip install biggles
     ```


___2. Install biskit___

```sh
git clone https://github.com/graik/biskit.git biskit -b biskit3
pip install -r biskit/requirements.txt
pip install -e biskit
```
If not already available, this will also install numpy, scipy, and BioPython. Replace `git clone` by the appropriate `tar xvf *tgz` command to start from an official Biskit release bundle.

___3. Test your installation___

Biskit comes with a unittest suite that can be run using the `test.py` script that is part of the library. First you have to figure out where your biskit installation went. Typically it will be in something like `/usr/local/lib/python3.7/site-packages/biskit`. If you have no idea, open a python interpreter and ... :

   ```python
   >>> import biskit
   >>> biskit.__path__
   ['/usr/local/lib/python3.7/site-packages/biskit']
   ```
   
Now run the biskit test suite, *except* those tests that require external programs (`-e exe`) or are tagged as `old` or `fails`:
 
   ```sh
~> python3 /usr/local/lib/python3.7/site-packages/biskit/test.py -e exe old fails
collecting  'biskit'
collecting  'biskit.core'
collecting  'biskit.exe'
collecting  'biskit.md'
collecting  'biskit.dock'
collecting  'biskit.future'
amberResidueLibrary.Test.test_amberResidueLibrary  ............ ok
amberResidues.Test.test_amberPrepParser  ...................... ok
atomCharger.Test.test_atomcharger  ............................ ok  [ 0.56s]
biounit.Test.test_BioUnit  .................................... ok  [ 0.57s]
[... long list of further test results ...]
future.residue.Test.test_invalidIndex  ........................ ok

----------------------------------------------------------------------
Ran 91 tests in 30.048s

OK

The test log file has been saved to: '<stdout>'

Test Coverage:
=============

2 out of 70 modules had no test case:
	 biskit .   wlc
	 biskit .   core .   difflib_old

SUMMARY:
=======

A total of 91 tests from 68 modules were run.
   - 91 passed
   - 0 failed
DONE
   ```
Once you have installed third-party software such as Pymol, Delphi, Xplor-NIH, DSSP, surfaceRacer, etc, you can re-run the test without the -e exe option. If you want to test individual biskit wrappers for a given program, simply call the wrapping python module which will execute this particular test. For example, if you have just installed Pymol, you can now run the biskit.exe.pymoler test case to ensure biskit and Pymol are properly working together:

    ```sh
   ~> python3 /usr/local/lib/python3.7/site-packages/biskit/exe/pymoler.py
    ```
    
This should open a Pymoler window with a short MD movie.


License
-------

Biskit is distributed under the GNU GPL version 3. See LICENSE.txt. Contact us, in case you prefer a different licensing model.
