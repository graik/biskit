Biskit
==========
**a software platform for structural bioinformatics**

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

License
-------

Biskit is distributed under the GNU GPL version 3. See COPYING. Contact us, in case you prefer a different licensing model.
