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

In this document:  * Release 2.1.0-beta
		   * Release 2.0.1
		   * Release 2.0.0
		   * Installation
		   * License
		   * Open issues

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

The biskit python library (module) biskit/Biskit needs to be in your
$PYTHONPATH and depends on other python libraries (numpy,
Scientific, biggles) which in turn depend on a few additional
libraries or programs (netcdf, gnuplot, Numeric).

The complete list is given here:
   http://biskit.pasteur.fr/install/short/

More detailed instructions for installing the different third-party
packages can be found here:
   http://biskit.pasteur.fr/install/libraries

The installation of helper applications is described here:
   http://biskit.pasteur.fr/install/applist
   http://biskit.pasteur.fr/install/applications


License
-------

Biskit is distributed under the GNU GPL. See license.txt!


Open issues
-----------

Please, check out the latest CVS snapshot and the list of open bugs at
http://sf.net/projects/biskit, before reporting any new problem!

The up-to-date list of open issues can be found at:

               http://biskit.pasteur.fr/doc/issues 

(or somewhere similar in case we still move it around. Use the search
field if the above link shouldn't work).


* most pickles in test/ are still using the older versions of
  PDBModel, this is partly intentional (ensuring backwards compability).

* not all scripts are tested

    The scripts in biskit/scripts are not included in any Test
    suites. Some of them haven't been used in quite a while. Please,
    report any problems to the bug tracker [Category `scripts`]!

* Amber* classes lack test cases

    We haven't yet written test cases for the AmberCrdEntropist,
    AmberEntropist, AmberParmBuilder. These classes are used for the
    entropy calculations in the Gruenberg et. al (2006) Structure
    paper. They should work with Amber 8. I (Raik) plan to migrate
    them to Amber 9 and write a proper entropy howto sometime in 2007.

* the parallel versions of the Mod scripts may have bugs

    We have started working on some projects that will use the
    modeling pipeline.  Any problems will hopefully be sorted out
    along the way. The non-parallel modeling (search_sequences,
    search_templates, clean_templates, align, model) should
    nevertheless work.

* API-Documentation issues

    The online API-documentation is usually lacking behind quite a
    bit. Also, epydoc gets confused by our practise to override the
    name of modules in the Biskit namespace with the classes they define.


For any other questions...
Did we mention the biskit home page yet ;-) ?
                   http://biskit.pasteur.fr

Good luck!
Johan & Raik
