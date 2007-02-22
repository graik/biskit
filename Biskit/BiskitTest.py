##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2006 Raik Gruenberg & Johan Leckner
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License as
## published by the Free Software Foundation; either version 2 of the
## License, or any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
## General Public License for more details.
##
## You find a copy of the GNU General Public License in the file
## license.txt along with this program; if not, write to the Free
## Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

## last $Author$
## last $Date$
## $Revision$

import unittest as U
import glob, types, re
import os.path

import Biskit
from Biskit.LogFile import StdLog, LogFile
import Biskit.tools as T

## System-wide test log 
TESTLOG = StdLog()

## categories
NORMAL = 0  ## standard test case
LONG   = 1  ## long running test case
PVM    = 2  ## depends on PVM
EXE    = 3  ## depends on external application


class BiskitTestError( Exception ):
    pass

class BiskitTest( U.TestCase):
    """
    Base class for Biskit test cases.
    BiskitTest adds some functionality over the standard L{unittest.TestCase}:

    * self.local reflects whether the Test is performed in the __main__ scope
      of the module it belongs to -- or as part of the whole test suite. 

    * Each test case can be classified by assigning flags to its
      static TAGS field -- this will be used by the test runner to
      filter out, e.g. very long tests or tests that require PVM.

    * self.log should capture all the output, by default it should go to the
      TESTLOG instance defined at the module level.

    Usage:
       Biskit test cases should be created by sub-classing BiskitTest and
       adding one or more test_* methods. setUp and tearDown can optionally
       be overriden but should be called in the overriding method.
    """
    ## categories for which this test case qualifies (class-wide)
    TAGS = [ NORMAL ]

    def setUp( self ):
        self.local =  self.__module__ == '__main__'
        self.log = TESTLOG

    def tearDown( self ):
        if self.local:
            globals().update( self.__dict__ )



class FilteredTestSuite( U.TestSuite ):
    """
    Collection of BiskitTests filtered by category tags.
    FilteredTestSuite silently ignores Test cases that are either

    * classified into any of the forbidden groups

    * not classified into any of the allowed groups

    The groups are defined in the L{Biskit.BiskitTest} module.

    By default (if initialized without parameters), all groups are allowed
    and no groups are forbidden.
    """

    def __init__( self, tests=(), allowed=[], forbidden=[] ):
        """
        @param tests: iterable of TestCases
        @type  tests: ( BiskitTest, )
        @param allowed: list of allowed tags
        @type  allowed: [ int ]
        @param forbidden : list of forbidden tags
        @type  forbidden : [ int ]
        """
        self._allowed   = allowed
        self._forbidden = forbidden
        U.TestSuite.__init__( self, tests=tests )


    def addTest( self, test ):
        """
        Add Biskit test case if it is matching the allowed and disallowed
        groups.
        @param test: test case
        @type  test: BiskitTest
        """
        assert isinstance( test, Biskit.BiskitTest.BiskitTest ), \
               'FilteredTestSuite only accepts BiskitTest instances not %r' \
               % test

        matches = [ g for g in test.TAGS if g in self._forbidden ]
        if len( matches ) > 0:
            return

        matches = [ g for g in test.TAGS if g in self._allowed ]

        if not self._allowed or (self._allowed and len(matches) > 0):
            U.TestSuite.addTest( self, test )
            

class BiskitTestLoader( object ):
    """
    A replacement for the unittest TestLoaders. It automatically
    collects all BiskitTests from a whole package (that means a
    folder with python files) rather than a single module.
    """

    def __init__( self, log=StdLog(),
                  allowed=[], forbidden=[], verbosity=2 ):

        self.allowed  = allowed
        self.forbidden= forbidden
        self.log = log
        self.verbosity = verbosity
        self.suite =  FilteredTestSuite( allowed=allowed, forbidden=forbidden )


    def modulesFromPath( self, path=T.projectRoot(), module='Biskit' ):
        """
        Import all python files of a package as modules. Sub-packages
        are ignored and have to be collected separately.
        @param path:  single search path for a package
        @type  path:  str
        @param module: name of the python package
        @type  module: str
        @return: list of imported python modules, see also L{__import__}
        @rtype : [ module ]
        @raise ImportError, if a python file cannot be imported
        """
        module_folder = module.replace('.', os.path.sep)
        
        files = glob.glob( os.path.join( path, module_folder,'*.py' ) )

        files = map( T.stripFilename, files )
        files = [ f for f in files if f[0] != '_' ]

        r = []

        for f in files:
            try:
                r += [ __import__( '.'.join([module, f]), globals(),
                                   None, [module]) ]
            except:
                pass  ## temporary // remove after testing

        return r


    def testsFromModules( self, modules ):
        """
        Extract all test cases from a list of python modules and add them to
        the internal test suite.
        @param modules: list of modules to be checked for BiskitTest classes
        @type  modules: [ module ]
        @param allowed: tags required for test to be considered, default: []
        @type  allowed: [ int ]
        @param forbidden: tags leading to the exclusion of test, default: []
        @type  forbidden: [ int ]
        @return: test cases (subject to the allowed and forbidden tags)
        @rtype: FilteredTestSuite
        """
        r = FilteredTestSuite( allowed=self.allowed, forbidden=self.forbidden )

        for m in modules:
            for i in m.__dict__.values():

                if type(i) is type \
                       and issubclass( i, Biskit.BiskitTest.BiskitTest )\
                       and i.__name__ != 'BiskitTest':
                    
                    suite = U.defaultTestLoader.loadTestsFromTestCase( i )
                    r.addTests( suite )

        return r


    def collectTests( self, path=T.projectRoot(), module='Biskit' ):
        """
        Add all BiskitTests found in a given module to the internal test suite.
        @param path:  single search path for a package
        @type  path:  str
        @param module: name of the python package
        @type  module: str
        """
        modules = self.modulesFromPath( path=path, module=module )
        self.suite.addTests( self.testsFromModules( modules ) )


    def report( self ):
        """
        Report how things went to stdout.
        """
        print '\nThe test log file has been saved to: %r'% self.log.fname
        total  = self.result.testsRun
        failed = len(self.result.failures) + len(self.result.errors)

        ## print a summary
        print '\nSUMMARY:\n=======\n'
        print 'A total of %i tests were run.' % total
        print '   - %i passed'% (total - failed)
        print '   - %i failed'% failed

        ## and a better message about which module failed 
        if failed:

            for test, ftrace in self.result.failures:
                print '      - failed: %s'% test.id()

            for test, ftrace in self.result.errors:
                print '      - error : %s'% test.id()



    def run( self ):

        runner = U.TextTestRunner( self.log.f(), verbosity=self.verbosity)
        self.result = runner.run( self.suite )

#########################
### Helper functions ####


def getOuterNamespace():
    """
    Fetch the namespace of the module/script running as __main__.
    @return: the namespace of the outermost calling stack frame
    @rtype: dict
    """
    import inspect

    try:
        frames = inspect.stack()
        f = frames[-1][0]   ## isolate outer-most calling frame
        r = f.f_globals
    finally:
        del frames, f

    return r


def extractTestCase( namespace ):
    """
    @return: first BisktTest found in given namespace
    @rtype: class
    """

    for i in namespace.values():

        if type(i) is type \
               and issubclass( i, Biskit.BiskitTest.BiskitTest )\
               and i.__name__ != 'BiskitTest':

            return i

    raise BiskitTestError, 'no BiskitTest class found in namespace'


def localTest( testclass=None, verbosity=2 ):
    """
    Perform the BiskitTest found in the scope of the calling module.
    After the test run, all fields of the BiskitTest instance are
    pushed into the global namespace so that they can be inspected in the
    interactive interpreter.
    
    @param testclass: BiskitTest-derived class [default: first one found]
    @type  testclass: class
    @param verbosity: verbosity level of TextTestRunner
    @type  verbosity: int

    @return: the test result object
    @rtype:  unittest.TestResult

    @raise BiskitTestError: if there is no BiskitTest-derived class defined
    """
    outer = getOuterNamespace()
    testclass = testclass or extractTestCase( outer )

    suite = U.TestLoader().loadTestsFromTestCase( testclass )
    runner= U.TextTestRunner(verbosity=verbosity)
    r = runner.run( suite )
    
    for t in suite._tests:
        outer.update( t.__dict__ )



if __name__ == '__main__':

    flog = LogFile( T.projectRoot()+'/test/test.log')

    l = BiskitTestLoader( allowed=[NORMAL],
                          forbidden=[PVM] )

    l.collectTests( module='Biskit')
    l.collectTests( module='Biskit.Mod' )
    l.run()
    l.report()


    print "DONE"
