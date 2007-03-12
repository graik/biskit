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

In this document:  * Release 2.0.1
		   * Installation
		   * License
		   * Open issues


Release 2.0.1
-------------

Release 2.0.1 introduces a new testing framework that is easier to
maintain and more convenient to use: The module Biskit.test
automatically collects and runs (unittest) test cases from across the
whole package and allows to skip or include tests depending on TAGS or
by package. The Biskit.PVM package has been re-arranged: The compilation
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
$PYTHONPATH and depends on other python libraries (Numeric,
Scientific, biggles) which in turn depend on a few additional
libraries or programs (netcdf, gnuplot).

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


Open issues in release 2.0.1
----------------------------

Please, check out the latest CVS snapshot and the list of open bugs at
http://sf.net/projects/biskit, before reporting any new problem!

The up-to-date list of open issues can be found at:

               http://biskit.pasteur.fr/doc/issues 

(or somewhere similar in case we still move it around. Use the search
field if the above link shouldn't work).


* some test pickles still need to be re-created on 64 bit

    see bug 1624505

* not all scripts are tested

    The scripts in biskit/scripts are not included in any Test
    suites. Some of them haven't been used in quite a while. Please,
    report any problems to the bug tracker [Category `scripts`]!

* Some Amber* classes lack test cases

    We haven't yet written test cases for the AmberCrdEntropist,
    AmberEntropist, AmberParmBuilder. These classes are used for the
    entropy calculations in the Gruenberg et. al (2006) Structure
    paper. They should work with Amber 8. I (Raik) plan to migrate
    them to Amber 9 and write a proper entropy howto sometime in 2007.

* the parallel versions of the Mod scripts may have bugs

    We have started working on some projects that will use the
    modeling pipeline.  Any problems should be sorted out along the
    way. The non-parallel modeling (search_sequences,
    search_templates, clean_templates, align, model) should
    nevertheless work.

* host management is a bit awkward

    Biskit.PVM keeps track of available computer nodes in a settings
    file `~/.biskit/hosts.dat`.  Nodes are classified into groups
    (own, shared, others) to which one can assign different
    nice-levels. RAM and CPU-number can also be given and some scripts
    use this to e.g. exclude nodes that are not up to the task from
    memory hungry calculations. The resulting hosts.dat format is a
    bit awkward and not easy to maintain (or understand). There are
    three workarounds:
    
    - study it and use it
    - use the native pvm host management (see pvm documentation)
    - come up with your own solution and just pass your list to the master


For any other questions...
Did we mention the biskit home page yet ;-) ?
                   http://biskit.pasteur.fr

Good luck!
Johan & Raik
