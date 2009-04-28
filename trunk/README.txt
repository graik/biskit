Biskit -- a software platform for structural bioinformatics
===========================================================

Please refer to 
                   http://biskit.pasteur.fr
               (mirrored at http://biskit.sf.net)

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
and Hmmer; interfaces to further software can be added
easily. Moreover, Biskit simplifies the parallelisation of
calculations via PVM (Parallel Virtual Machine).

In this document:  * Release 2.3
                   * Release 2.2
		   * Release 2.1.0-beta
		   * Release 2.0.1
		   * Release 2.0.0
		   * Installation
		   * License
		   * Open issues

Release 2.3
-----------

Release 2.3 simplifies the installation of Biskit. The package can now
be installed system-wide with the standard 'python setup.py install'
command. This required a change to Biskit's directory layout: The
'test' and the 'external' folder were moved from the project root into
the Biskit python package. That means:

  old               ->    new
* biskit/external   ->    biskit/Biskit/data
* biskit/test       ->    biskit/Biskit/testdata

The idea here was to have all extra folders nicely bundled with the
python package rather than spreading them accross your system. This
also applies to the scripts and docs folder. These two remain at the
project root, but are copied into the package folder during system
installations. 

The advantage of this new setup is that we can now start creating
Debian, RPM and even Windows installation packages.

Other changes: 

* fixed Intervor wrapper

* some minor bug fixes

* some work-in-progress packages


Release 2.2
------------

Release 2.2 adds many bug fixes and a few new features to the
previous beta release:

* switch from CVS to SVN version control

* PDBModel: new report() and plot() function to get a fast overview
  over the content of a PDB

* PDBModel: PDBModel( '3TGI' ) will fetch the PDB entry 3TGI directly 
  from the ncbi database or a local mirror

* PDBModel: new maskDNA() function

* new CommandLine.py is a candidate to replace the old tools.cmdDict command
  line option parsing for scripts; with type-casting and test fixture
  built in

* Statistics/ROCalyzer module calculates ROC areas and their statistic 
  significance

* new Intervor module wraps the Voronoi interface calculator of Cazal et
  al (but Frederic hasn't got around to offer a intervor version for
  download yet... come on Frederic ;o)! )


Release 2.1.0-beta
------------------

Release 2.1.0 introduces a lot changes to the very core of Biskit:

* PDBModel overhaul -- the atom dictionaries are gone and all the
  information is now unified into `atoms` and `residues`
  ProfileCollections, take, compress etc. are **much** more efficient,
  profiles can be accessed in a more intuitive way, chain borders are
  tracked more consistently, ... The new PDBModel is, as usual,
  backward-compatible to pickles of previous versions.

* Migration from Numeric to Numpy -- plenty of bugs in the old Numeric
  have forced us to migrate earlier than anticipated. Unfortunately
  some third-party libraries like Scientific and Biopython still
  depend on Numeric -- so both libraries need to be installed. Also we
  needed to convert the biggles module to the new numpy library.

All tests I can currently perform (this excludes PVM, homology
modeling and some interfaces to external programs) are running through
without problems. But since the changes are very widespread and at the
center of the package, I expect that there are some bugs still
awaiting discovery.  Furthermore, I haven't yet converted the scripts
folder to the new numpy and, since we are right now lacking a good
test suite for this part of Biskit, many of the biskit/scripts will
probably be broken (albeit not difficult to fix). 

For this reason, we label this one as a "beta" release. Users who are
mostly interested in reproducing one of the existing workflows
(docking or homology modeling scripts) should for the moment stick to
release 2.0.1. However, for everyone else, the improvements should
outweight potential bugs and you spare yourself a later
migration. Please don't hesitate to report any problems!

Release 2.0.1
-------------

Release 2.0.1 introduces a new testing framework that is easier to
maintain and more convenient to use: The module Biskit.test
automatically collects and runs (unittest) test cases from across the
whole package and allows to skip or include tests depending on TAGS or
package. The Biskit.PVM package has been re-arranged: The compilation
of the built-in pypvm_core module used to be a notorious trouble
maker. We have now outsourced and re-packaged pypvm into a standard
distutils package that compiles and installs more smoothly. Besides,
many bugs were fixed (especially in the Biskit.Mod module) and both
code and documentation were polished in several places.

Release 2.0.0
-------------

This is the first 'official' release of Biskit. However, we are
working on (and using) this project since 2001 and what you are
looking at could most fairly be called generation 2 or 3 -- that's why
we are not starting at revision 1.0 here. Until now, we still have
been referring everybody to the CVS server for fetching the
package. This is still recommended if you want to do serious work with
it but the pre-packaged release may be more convenient for evaluation
purposes or simple tasks.


Installation
------------

Up-to-date installation instructions can be found online:

* http://biskit.pasteur.fr/install

(Some older instructions are still kept in docs/.)

The biskit python library (module) depends on other python libraries (numpy,
ScientificPython, biggles) which in turn depend on a few additional
libraries or programs (netcdf, gnuplot, plotutils). Some functions of 
Biskit moreover rely on scipy and on biopython, although we try to
keep this dependence optional.

The complete list is given here:
   http://biskit.pasteur.fr/install/short/

You can now install Biskit by running:
   python setup.py install
in the root of the project's directory. Alternatively,
you can simply copy or link the Biskit folder (the one containing
__init__.py) into your $PYTHONPATH.

More detailed instructions for installing the different third-party
packages can be found here:
   http://biskit.pasteur.fr/install/libraries

The installation of helper applications is described here:
   http://biskit.pasteur.fr/install/applist
   http://biskit.pasteur.fr/install/applications


License
-------

Biskit is distributed under the GNU GPL version 3. See license.txt!


Open issues
-----------

Please, check out the latest SVN snapshot and the list of open bugs at
http://sf.net/projects/biskit, before reporting any new problem!

The up-to-date list of open issues can be found at:

               http://biskit.pasteur.fr/doc/issues 

(or somewhere similar in case we still move it around. Use the search
field if the above link shouldn't work).


* most pickles in testdata/ are still using the older versions of
  PDBModel, this is partly intentional (ensuring backwards compability).

* not all scripts are tested

    The scripts in biskit/scripts are not included in any Test
    suites. Some of them haven't been used in quite a while. Please,
    report any problems to the bug tracker [Category `scripts`]!

* Amber* classes lack test cases

    We haven't yet written test cases for the AmberCrdEntropist,
    AmberEntropist, AmberParmBuilder. These classes are used for the
    entropy calculations in the Gruenberg et. al (2006) Structure
    paper. They should work with Amber 8. 

* the parallel versions of the Mod scripts may have bugs

    The non-parallel modeling (search_sequences,
    search_templates, clean_templates, align, model) should
    nevertheless work.

* API-Documentation issues

    The online API-documentation is usually lagging behind quite a
    bit. Also, epydoc gets confused by our practise to override the
    name of modules in the Biskit namespace with the classes they define.


For any other questions...
Did we mention the biskit home page yet ;-) ?
                   http://biskit.pasteur.fr

Good luck!
Johan & Raik
