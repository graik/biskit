[![Build Status](https://travis-ci.org/graik/biskit.svg?branch=biskit3)](https://travis-ci.org/graik/biskit)
[![Coverage Status](https://coveralls.io/repos/github/graik/biskit/badge.svg?branch=biskit3)](https://coveralls.io/github/graik/biskit?branch=biskit3)
Biskit
==========
___a software platform for structural bioinformatics___

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
easily. Moreover, Biskit simplifies the parallelisation of
calculations via PVM (Parallel Virtual Machine).

Short Installation Instructions
--------------------------------

```sh
git clone https://github.com/graik/biskit.git biskit
pip install -r biskit/requirements.txt
pip install -e biskit
```
Replace `git clone` by the appropriate `tar xvf *tgz` command to start from an official Biskit release bundle.

*Note:* if not already available, this will automatically install numpy, scipy, and BioPython. It will *not* however install the biggles plotting package, which is not strictly required but highly recommended. Biggles is not available via `pip`. On Ubuntu, it can be installed with `sudo apt-get install python-pybiggles`. See http://biskit.pasteur.fr/install/short for step-by-step instructions on alternative installation methods.


License
-------

Biskit is distributed under the GNU GPL version 3. See LICENSE.txt. Contact us, in case you prefer a different licensing model.
