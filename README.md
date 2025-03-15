biskit
==========
___a software platform for structural bioinformatics___

[![Build Status](https://travis-ci.org/graik/biskit.svg?branch=master)](https://travis-ci.org/graik/biskit)
[![Coverage Status](https://coveralls.io/repos/github/graik/biskit/badge.svg?branch=master&service=github)](https://coveralls.io/github/graik/biskit?branch=master)

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

The core functionalities and many of the wrappers have been migrated from the old Biskit 2 branch. However, some smaller modules and most of the scripts didn't make it. All the modules that have **not** yet been migrated are in the old ```Biskit2``` folder (which cannot be imported any longer). After installation, the python 3 Biskit is imported with ```import biskit```.

Please refer to 
            **http://biskit.pasteur.fr** (note: http rather than https)
for much more detailed usage instructions based on the older Biskit 2 package. 

Short Installation Instructions
--------------------------------

___1. Install (plotting) dependencies___

The biskit library itself does not strictly need `biggles` and you can safely install it later. However, **biggles** (https://biggles-plot.github.io/) is still an excellent plotting library with a relatively intuitive syntax that is used throughout biskit and several unittests depend on it. If you are sitting on a linux / unix PC, you may as well get it set up right from the start. Unfortunately, biggles does not any longer seem to support installation in OSX. 

On Debian / Ubuntu:
  *  ```sh
     sudo apt-get install libplot-dev plotutils  ## needed for biggles compilation
     pip3 install biggles
     ```

This used to work on Mac OS-X but compilation fails:
  * install Quartz (https://www.xquartz.org/)
  *  ```sh
     brew install plotutils --with-x11
     pip3 install biggles
     ```

___2. Install biskit___

```sh
git clone https://github.com/graik/biskit.git biskit3
pip3 install -e biskit3
```
*Note:* The `-e` option will create an "editable" biskit installation where the git-controlled `biskit3` folder is not copied but sym-linked into your system's (or virtualenv) python `site-packages` folder [see explanation](http://codumentary.blogspot.com/2014/11/python-tip-of-year-pip-install-editable.html). *Alternatively*, run `python3 biskit3/setup.py install`
in order to create an actual copy of the biskit source code inside your Python 3 `site-packages` folder.

Depending on your environment, `pip install -e` will always install all the needed dependencies but it may or may not also create the link for putting biskit into the `$PYTHONPATH`. It seems to work within the virtualenv but not, e.g. in the OSX terminal. Test by changing to another folder, run `python` and run the following import command:
```py
>>> import biskit
```
Should this fail, you can instead add the biskit3 folder manually to your PYTHONPATH. Assuming you cloned the biskit project into `~/py/biskit3`, the following would work:
```sh
export PYTHONPATH=$PYTHONPATH:~/py/biskit3
```
Append this line to your `.bashrc` or `.bash_profile` file. zsh users should put it into `.zshenv`. 

___3. Test your installation___

Biskit comes with a unittest suite that can be run using the `test.py` script that is part of the library. First you have to figure out where your biskit installation went. If you used the `pip3 install -e` command above, biskit will still be in the same location where your `git clone` created it. Otherwise, it will be in something like `/usr/local/lib/python3.7/site-packages/biskit`. If you have no idea, open a python interpreter and ... :

  * ```python
    >>> import biskit
    >>> biskit.__path__
    ['/usr/local/lib/python3.13/site-packages/biskit']
    ```
   
Now run the biskit test suite, *except* those tests that require external programs (`-e exe`) or are tagged as `old` or `biggles`:
 
   ```sh
   ~> python3 biskit/test.py -e exe old biggles
   ```
Once you have installed third-party software such as Pymol, Delphi, Xplor-NIH, DSSP, surfaceRacer, etc, you can re-run the test without the -e exe option. If you want to test individual biskit wrappers for a given program, simply call the wrapping python module which will execute this particular test. For example, if you have just installed Pymol, you can now run the biskit.exe.pymoler test case to ensure biskit and Pymol are properly working together:

   ```sh
   ~> python3 biskit/exe/pymoler.py
   ```
    
This should open a Pymoler window with a short MD movie and will, generally, give you a more detailed test output.


License
-------

Biskit is distributed under the GNU GPL version 3. See LICENSE.txt. Contact us, in case you prefer a different licensing model.
