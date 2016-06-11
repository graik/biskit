#!/usr/bin/env python
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2016 Raik Gruenberg & Johan Leckner
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License as
## published by the Free Software Foundation; either version 3 of the
## License, or any later version.
##
## This program is distributec!d in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
## General Public License for more details.
##
## You find a copy of the GNU General Public License in the file
## license.txt along with this program; if not, write to the Free
## Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

'''
Automatic Biskit module testing
===============================

Module tests are extremely important to keep the Biskit package
manageable and functional. Every Biskit module should contain a class
derrived from L{BiskitTest} (conventionally called C{Test}) that
comprises one or more test functions.  L{BiskitTestLoader} then
automatically extracts all BiskitTest child classes from the whole
package and bundles them into a L{FilteredTestSuite}.

Test cases that are L{LONG}, depend on L{PVM} or external programs
(L{EXE}) or are obsolete (L{OLD}) should be marked with the
appropriate tag(s) by overriding L{BiskitTest.TAGS} on the class
level. Test cases can then be filtered out or included into the
L{FilteredTestSuite} according to their TAGS.

Originally, testing code was simply in the __main__ section of each
module where we could execute it directly from emacs (or with python
-i) for interactive debugging. By comparison, unittest test fixtures
are less easy to execute stand-alone and intermediate variables remain
hidden within the test instance. The L{localTest}() removes this
hurdle and runs the test code of a single module as if it would be
executed directly in __main__. Simply putting the localTest() function
without parameters into the __main__ section of your module is
enough. Test.test_* methods should assign intermediate and final
results to self.|something| variables -- L{localTest} will then push
all self.* fields into the global namespace for interactive debugging.

Usage
=====

  By way of example, a Test case for MyModule would look like this::

    class MyClass:
        ...

    ### Module testing ###
    import Biskit.test as BT

    class Test(BT.BiskitTest):
        """Test MyModule"""

        TAGS = [ BT.LONG ]

        def test_veryLongComputation( self ):
            """MyModule.veryLongComputation test"""

            self.m = MyClass()
            self.result = self.m.veryLongComputation()

            if self.local:   ## only if the module is executed directly
                print self.result 

            self.assertEqual( self.result, 42, 'unexpected result' )


    if __name__ == '__main__':

        ## run Test and push self.* fields into global namespace
        BT.localTest( )

        print result  ## works thanks to some namespace magic in localTest


Note:
        - If TAG is not given, the test will have the default NORMAL tag.
        - Names of test functions **must** start with C{test_}.
        - The doc string of test_* will be reported as id of this test.
        - this module also acts as the script to collect and run the tests
'''

import unittest as U
import glob
import os.path

import Biskit
from LogFile import StdLog, LogFile
import tools as T

## categories
NORMAL = 0  ## standard test case
LONG   = 1  ## long running test case
PVM    = 2  ## depends on PVM
EXE    = 3  ## depends on external application
OLD    = 5  ## is obsolete
SCRIPT = 6  ## a script test case

class BiskitTestError( Exception ):
    pass

class BiskitTest( U.TestCase):
    """
    Base class for Biskit test cases.

    BiskitTest adds some functionality over the standard
    C{unittest.TestCase}:

      - self.local reflects whether the Test is performed in the
        __main__ scope of the module it belongs to -- or as part of
        the whole test suite.

      - Each test case can be classified by assigning flags to its
        static TAGS field -- this will be used by the test runner to
        filter out, e.g. very long tests or tests that require PVM.

      - self.log should capture all the output, by default it goes to
        the TESTLOG instance defined at the module level.

    Usage:
    ======

      Biskit test cases should be created by sub-classing BiskitTest
      and by adding one or more test_* methods (each method starting
      with 'test_' is treated as a separate test). The doc string of a
      test_* method becomes the id reported by the TextTestRunner.

      L{prepare} should be overriden for the definition of permanent
      input and temporary output files. L{cleanUp} should be overriden
      for clean up actions and to remove temporary files (use
      L{Biskit.tools.tryRemove}). L{cleanUp} is not called if
      BiskitTest is set into debugging mode (see L{BiskitTest.DEBUG}).

      The L{TAGS} field should be overriden to reflect any special categories
      of the test case (L{LONG}, L{PVM}, L{EXE} or L{OLD}).
    """
    #: categories for which this test case qualifies (class-wide)
    TAGS  = [ NORMAL ]

    #: debug mode -- don't delete temporary data
    DEBUG = False

    #: System-wide test log
    TESTLOG = StdLog()

    #: System-wide verbosity level
    VERBOSITY= 2

    def prepare(self):
        """Hook for test setup."""
        pass

    def cleanUp( self ):
        """Hook for post-test clean up; skipped if DEBUG==True."""
        pass

    def setUp( self ):
        self.local =  self.__module__ == '__main__'
        self.log =       getattr( self, 'log', BiskitTest.TESTLOG )
##      self.verbosity = getattr( self, 'verbosity', BiskitTest.VERBOSITY )
##      self.debugging = getattr( self, 'debugging', BiskitTest.DEBUG )
        self.prepare()

    def tearDown( self ):
        if not self.DEBUG:
            self.cleanUp()


class FilteredTestSuite( U.TestSuite ):
    """
    Collection of BiskitTests filtered by category tags.
    FilteredTestSuite silently ignores Test cases that are either

    * classified into any of the forbidden categories

    * not classified into any of the allowed categories

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
        assert isinstance( test, Biskit.test.BiskitTest ), \
               'FilteredTestSuite only accepts BiskitTest instances not %r' \
               % test

        matches = [ g for g in test.TAGS if g in self._forbidden ]
        if len( matches ) > 0:
            return

        matches = [ g for g in test.TAGS if g in self._allowed ]

        if not self._allowed or (self._allowed and len(matches) > 0):
            U.TestSuite.addTest( self, test )


import time

class PrettyTextTestResult( U.TextTestResult ):
    """
    Helper class for TextTestRunner.
    Only print either test description or doc-string.
    """

    
    def getDescription(self, test):
        s = test.id() + '  '
        ## remove leading package name
        if s.count('.') > 2:
            s = s[s.index('.')+1:]
        return s
    
    def startTest(self, test):
        super(U.TextTestResult, self).startTest(test)
        if self.showAll:
            desc = self.getDescription(test)
            self.stream.write(desc.ljust(60,'.'))
            self.stream.write("... ")
            self.stream.flush()
        self.startclock = time.clock()

    def addSuccess(self, test):
        ## super(U.TextTestResult, self).addSuccess(test)
        dt = time.clock() - self.startclock
        
        if self.showAll:
            if dt > 0.5:
                self.stream.writeln('ok  [%5.2fs]' % dt)
            else:
                self.stream.writeln('ok')

        elif self.dots:
            self.stream.write('.')
            self.stream.flush()


class SimpleTextTestRunner( U.TextTestRunner ):
    """
    Convince TextTestRunner to use the flushing text output rather
    than the default one.
    """

    def _makeResult(self):
        return PrettyTextTestResult(self.stream, self.descriptions,
                                    self.verbosity)


class BiskitTestLoader( object ):
    """
    A replacement for the unittest TestLoaders. It automatically
    collects all BiskitTests from a whole package (that means a
    folder with python files).
    """

    def __init__( self, log=StdLog(),
                  allowed=[], forbidden=[], verbosity=2, debug=False ):
        """
        @param log: log output target [default: L{Biskit.StdLog}]
        @type  log: Biskit.LogFile
        @param allowed: tags required for test to be considered, default: []
        @type  allowed: [ int ]
        @param forbidden: tags leading to the exclusion of test, default: []
        @type  forbidden: [ int ]
        @param verbosity: verbosity level for unittest.TextTestRunner
        @type  verbosity: int
        """

        self.allowed  = allowed
        self.forbidden= forbidden
        self.log = log
        self.verbosity = verbosity
        self.debugging = debug
        self.suite =  FilteredTestSuite( allowed=allowed, forbidden=forbidden )
        self.modules_untested = []  #: list of modules without test cases
        self.modules_tested = []    #: list of modules containing test cases
        self.result = U.TestResult() #: will hold test result after run()


    def modulesFromPath( self, path=T.projectRoot(), module='Biskit' ):
        """
        Import all python files of a package as modules. Sub-packages
        are ignored and have to be collected separately.

        @param path:  single search path for a package
        @type  path:  str
        @param module: name of the python package [default: Biskit]
        @type  module: str

        @return: list of imported python modules, see also __import__
        @rtype : [ module ]

        @raise ImportError: if a python file cannot be imported
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


    def addTestsFromModules( self, modules ):
        """
        Extract all test cases from a list of python modules and add them to
        the internal test suite.
        @param modules: list of modules to be checked for BiskitTest classes
        @type  modules: [ module ]
        """
        for m in modules:
            tested = 0

            for i in m.__dict__.values():

                if type(i) is type and \
                   issubclass( i, Biskit.test.BiskitTest ) and \
                   i.__name__ != 'BiskitTest':

                    suite = U.defaultTestLoader.loadTestsFromTestCase( i )
                    self.suite.addTests( suite )
                    tested = 1

            if tested:
                self.modules_tested += [m]
            else:
                self.modules_untested += [m]


    def collectTests( self, path=T.projectRoot(), module='Biskit' ):
        """
        Add all BiskitTests found in a given module to the internal test suite.
        @param path:  single search path for a package
        @type  path:  str
        @param module: name of the python package
        @type  module: str
        """
        modules = self.modulesFromPath( path=path, module=module )
        self.addTestsFromModules( modules )


    def __moduleNames(self, modules ):

        modules = [ m.__name__ for m in modules ]   ## extract name
        modules = [ m.replace('.',' .   ') for m in modules ]

        return modules


    def report( self ):
        """
        Report how things went to stdout.
        """
        print '\nThe test log file has been saved to: %r'% self.log.fname
        total  = self.result.testsRun
        failed = len(self.result.failures) + len(self.result.errors)

        m_tested  = len( self.modules_tested )
        m_untested= len( self.modules_untested)
        m_total  = m_tested + m_untested

        ## report coverage
        print '\nTest Coverage:\n=============\n'
        print '%i out of %i modules had no test case:' % (m_untested,m_total)
        for m in self.__moduleNames( self.modules_untested) :
            print '\t', m

        ## print a summary
        print '\nSUMMARY:\n=======\n'
        print 'A total of %i tests from %i modules were run.' %(total,m_tested)
        print '   - %i passed'% (total - failed)
        print '   - %i failed'% failed

        ## and a better message about which module failed 
        if failed:

            for test, ftrace in self.result.failures:
                print '      - failed: %s'% test.id()

            for test, ftrace in self.result.errors:
                print '      - error : %s'% test.id()


    def run( self, dry=False ):
        """
        @param dry: do not actually run the test but just set it up [False]
        @type  dry: bool
        """
        ## push global settings into test classes
        for testclass in self.suite:
            testclass.DEBUG = self.debugging
            testclass.VERBOSITY = self.verbosity
            testclass.TESTLOG = self.log

        runner = SimpleTextTestRunner(self.log.f(), verbosity=self.verbosity)
        if not dry:
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
        #print "DEBUG 1 ", inspect.getouterframes(f)
        #print
        #frames = inspect.getouterframes( inspect.currentframe() )
        #print "DEBUG 2 ", len( frames )
        #print 
        #print "DEBUG 2 "
        #for f in frames:
            #print f
            #print
    finally:
        del frames, f

    return r

def getCallingNamespace( callpattern='localTest' ):
    """
    Fetch the namespace of the module/script from which the localTest method
    is running. This is more general than the getOuterNamespace method and 
    also works with python -m <module> or programs that wrap the python call
    into their own methods.
    @param callpattern: str, name of calling method to look for ['localTest']
    @return: the namespace of the outermost calling stack frame
    @rtype: dict  
    """
    import inspect

    try:
        frames = inspect.getouterframes( inspect.currentframe() )

        i = 0
        while frames[i][3] != callpattern:
            i += 1
            if i >= len( frames ):
                raise BiskitTestError, \
                      'cannot find %s in current stack trace' % callpattern

        f = frames[i + 1][0]
        r = f.f_globals
    finally:
        del frames, f, i

    return r


def extractTestCases( namespace ):
    """
    @return: all BisktTest child classes found in given namespace
    @rtype: [ class ]
    @raise BiskitTestError: if there is no BiskitTest child
    """
    r =[]

    for i in namespace.values():

        if type(i) is type \
           and issubclass( i, Biskit.test.BiskitTest )\
           and i.__name__ != 'BiskitTest':

            r += [i]

    if not r:
        raise BiskitTestError, 'no BiskitTest class found in namespace'

    return r


def localTest( testclass=None, verbosity=BiskitTest.VERBOSITY,
               debug=BiskitTest.DEBUG, log=BiskitTest.TESTLOG ):
    """
    Perform the BiskitTest(s) found in the scope of the calling module.
    After the test run, all fields of the BiskitTest instance are
    pushed into the global namespace so that they can be inspected in the
    interactive interpreter. The BiskitTest instance itself is also
    put into the calling namespace as variable 'self', so that test code
    fragments referring to it can be executed interactively.

    @param testclass: BiskitTest-derived class [default: first one found]
    @type  testclass: class
    @param verbosity: verbosity level of TextTestRunner
    @type  verbosity: int
    @param debug: don't delete temporary files (skipp cleanUp) [0]
    @type  debug: int

    @return: the test result object
    @rtype:  unittest.TestResult

    @raise BiskitTestError: if there is no BiskitTest-derived class defined
    """
    ## get calling namespace
    outer = getCallingNamespace( callpattern='localTest' )
    if testclass:
        testclasses = [testclass]
    else:
        testclasses = extractTestCases( outer )

    suite = U.TestSuite()
    for test in testclasses:
        suite.addTests( U.TestLoader().loadTestsFromTestCase( test ) )

    for test in suite:
        test.DEBUG = debug
        test.VERBOSITY = verbosity
        test.TESTLOG = log

    runner= U.TextTestRunner(verbosity=verbosity)
    r = runner.run( suite )

    for t in suite._tests:
        outer.update( t.__dict__ )
        outer.update( {'self':t })

    return r

###############
## Mock test ##

class Test(BiskitTest):
    """Mock test, test doesn't test itself"""
    pass

########################
### Script functions ###

def _use( defaults ):
    print """
Run unittest tests for biskit.

    test.py [-i |include tag1 tag2..| -e |exclude tag1 tag2..|
             -p |package1 package2..|
             -v |verbosity| -log |log-file| -nox ]

    i    - include tags, only run tests with at least one of these tags   [All]
    e    - exclude tags, do not run tests labeled with one of these tags  [old]
         valid tags are:
             long   - long running test case
             pvm    - depends on PVM
             exe    - depends on external application
             old    - is obsolete
         (If no tags are given to -i this means all tests are included)

    p    - packages to test, e.g. Biskit Biskit.Dock Biskit.Mod           [All]
    v    - int, verbosity level, 3 switches on several graphical plots      [2]
    log  - path to logfile (overriden); empty -log means STDOUT        [STDOUT]
    nox  - suppress test plots                                          [False]
    dry  - do not actually run the test but just collect tests          [False]

Examples:

    * Run all but long or obsolete tests from Biskit and Biskit.Dock:
    test.py -e old long  -p Biskit Biskit.Dock

    * Run only PVM-dependent tests of the Biskit.Mod sub-package:
    test.py -i pvm  -p Biskit.Mod


Default options:
"""
    for key, value in defaults.items():
        print "\t-",key, "\t",value

    sys.exit(0)


def _str2tags( s ):
    """convert list of string options to list of valid TAGS"""
    try:
        r = [ x.upper() for x in s if x ] ## to list of uppercase str
        r = [ eval( x ) for x in r ]      ## to list of int
    except:
        EHandler.error('unrecognized tags: %r'%s)

    return r

def _convertOptions( o ):
    o['i'] = _str2tags( T.toList( o['i'] ) )
    o['e'] = _str2tags( T.toList( o['e'] ) )
    o['v'] = int( o['v'] )
    o['nox'] = ('nox' in o)
    o['dry'] = ('dry' in o)
    o['debug'] = ('debug' in o)
    if o['log']:
        o['log'] = LogFile( o['log'] )
    else:
        o['log'] = StdLog()
    o['p'] = T.toList( o['p'] )

if __name__ == '__main__':

    from Biskit import EHandler
    import sys


    defaults = {'i':'',
                'e':'old',
                'p':['Biskit', 'Biskit.Dock', 'Biskit.Mod', 'Biskit.PVM', 
                     'Biskit.Statistics'],
                'v':'2',
                'log': '', ##T.testRoot()+'/test.log',
                }

    o = T.cmdDict( defaults )

    if len( sys.argv ) == 1 and 'test.py' in sys.argv[0]:
        _use( defaults )

    _convertOptions( o )

    BiskitTest.VERBOSITY = o['v']
    BiskitTest.DEBUG = o['debug']

    l = BiskitTestLoader( allowed=o['i'], forbidden=o['e'],
                          verbosity=o['v'], log=o['log'], debug=o['debug'])


    for package in o['p']:
        print 'collecting ', repr( package )
        l.collectTests( module=package )

    l.run( dry=o['dry'] )
    l.report()

    print "DONE"
