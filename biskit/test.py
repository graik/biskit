#!/usr/bin/env python

##   Copyright 2014 - 2017 Raik Gruenberg, All Rights Reserved
'''
Automatic package testing
==========================

This is a fully self-contained extension to the standard PyUnit testing 
framework. The main features are:

* fully automatic collection of test cases from a package
* filtering and grouping of test cases by 'tags' and sub-packages
* test execution automatically kicks in if a module is run stand-alone
* test variables are pushed into global name space for interactive debugging
* the test module doubles as script for running the tests

Originally, our testing code was simply in the __main__ section of each
module where we could execute it directly from emacs (or with python
-i) for interactive debugging. By comparison, unittest test fixtures
are less easy to execute stand-alone and intermediate variables remain
hidden within the test instance. 

The :class:`localTest`() method removes this hurdle and runs the test code of a
single module as if it would be executed directly in __main__. Simply putting
the localTest() function without parameters into the __main__ section of your
module is enough. Your Test.test_* methods should assign intermediate and
final results to self.|something| variables -- :class:`localTest` will then push all
self.* fields into the global namespace for interactive debugging.

To get started, every module in your package should contain one or more classes
derrived from :class:`BiskitTest` (conventionally called C{Test}) that
each contains one or more test_* functions. :class:`BiskitTestLoader` then
automatically extracts all BiskitTest child classes from the whole
package and bundles them into a :class:`FilteredTestSuite`. Note, BiskitTest is 
derrived from the standard unittest.TestCase -- refer to the Python 
documentation for details on test writing.

Usage
=====

  By way of example, a Test case for MyModule would look like this::

    class MyClass:
        ...

    ### Module testing ###
    import testing

    class Test(testing.BiskitTest):
        """Test MyModule"""

        TAGS = [ testing.LONG ]

        def test_veryLongComputation( self ):
            """MyModule.veryLongComputation test"""

            self.m = MyClass()
            self.result = self.m.veryLongComputation()

            if self.local:   ## only if the module is executed directly
                print self.result 
                
            self.assertEqual( self.result, 42, 'unexpected result' )


    if __name__ == '__main__':

        ## run Test and push self.* fields into global namespace
        testing.localTest( )

        print result  ## works thanks to some namespace magic in localTest


Note:
        - If TAG is not given, the test will have the default NORMAL tag.
        - Names of test functions **must** start with C{test_}.
        - The doc string of test_* will be reported as id of this test.

This module also acts as the script to collect and run the tests. Run it without
arguments for help.
'''

import unittest as U
import glob
import os.path
import logging
import sys

## CONFIGURATION -- adapt the following values to your own Python package

#: list all packages from which test cases are collected by default (-p)
DEFAULT_PACKAGES = ['biskit', 'biskit.core', 'biskit.exe', 'biskit.md',
                    'biskit.dock', 'biskit.future']

#: tests with the following tags are excluded by default (override with -e)
DEFAULT_EXCLUDE  = ['old', 'extra', 'biggles']

## END OF CONFIGURATION

## categories
NORMAL = 0  ## standard test case
LONG   = 1  ## long running test case
PVM    = 2  ## depends on PVM
EXE    = 3  ## depends on external application
EXTRA  = 4  ## tests not routinely run
OLD    = 5  ## is obsolete
SCRIPT = 6  ## a script test case
FAILS  = 7  ## test known to currently fail (but procrastinated)
BIGGLES= 8  ## requires Biggles module

class BiskitTestError( Exception ):
    pass

########################################
### supporting file handling methods ###


def absfile( filename, resolveLinks=1 ):
    """
    Get absolute file path::
      - expand ~ to user home, change
      - expand ../../ to absolute path
      - resolve links
      - add working directory to unbound files ('ab.txt'->'/home/raik/ab.txt')

    :param filename: name of file
    :type  filename: str
    :param resolveLinks: eliminate any symbolic links (default: 1)
    :type  resolveLinks: 1|0
    
    :return: absolute path or filename
    :rtype: string

    :raise IOError: if a ~user part does not translate to an existing path
    """
    if not filename:
        return filename
    r = os.path.abspath( os.path.expanduser( filename ) )

    if '~' in r:
        raise IOError('Could not expand user home in %s' % filename)

    if resolveLinks:
        r = os.path.realpath( r )
    r = os.path.normpath(r)
    return r

def packageRoot():
    """
    The folder containing the parent python package.
    :return: str, absolute path of the root of current project
    """
    f = absfile(__file__)
    return absfile( os.path.join( os.path.split( f )[0], '..') )

def stripFilename( filename ):
    """
    Return filename without path and without ending.
    :param filename: str, name of file
    :return: str, base filename
    """
    name = os.path.basename( filename )      # remove path
    try:
        if name.find('.') != -1:
            name = name[: name.rfind('.') ]     # remove ending
    except:
        pass  ## just in case there is no ending to start with...
    return name


#########################
### Core test library ###

class LogHandle:
    """
    Thin wrapper around an open file handle that adds a `add()` and `add_nobreak`
    method to mirror biskit.LogFile
    """
    def __init__(self, handle=None):
        """
        :param handle: stream or file handle open for writing
        """
        self.handle = handle or sys.stdout
        self.name = handle.name
        
    def write(self, s):
        self.handle.write(s)
        
    def writeln(self, s):
        self.handle.write(s+'\n')
        
    def add(self, s):
        self.writeln(s)
        self.handle.flush()
        
    def add_nobreak(self,s):
        self.handle.write(s)
        self.handle.flush()
        
    def flush(self):
        return self.handle.flush()

class BiskitTest( U.TestCase):
    """
    New base class for test cases.

    BiskitTest adds some functionality over the standard
    C{unittest.TestCase}:

      - self.local reflects whether the Test is performed in the
        __main__ scope of the module it belongs to -- or as part of
        the whole test suite.

      - Each test case can be classified by assigning flags to its
        static TAGS field -- this will be used by the test runner to
        filter out, e.g. very long tests or tests that require PVM.

      - self.log (an open file handle) should capture all the output, 
        by default it reverts to BiskitTest.TESTLOG which opens STDOUT

    Usage:
    ======
    
      Test cases should be created by sub-classing BiskitTest
      and by adding one or more test_* methods (each method starting
      with 'test_' is treated as a separate test). The doc string of a
      test_* method becomes the id reported by the TextTestRunner.

      :class:`prepare` should be overriden for the definition of permanent
      input and temporary output files. :class:`cleanUp` should be overriden
      for clean up actions and to remove temporary files (use
      :class:`Biskit.tools.tryRemove`). :class:`cleanUp` is not called if
      BiskitTest is set into debugging mode (see :class:`BiskitTest.DEBUG`).

      The :class:`TAGS` field should be overriden to reflect any special categories
      of the test case (:class:`LONG`, :class:`PVM`, :class:`EXE` or :class:`OLD`).
    """
    #: categories for which this test case qualifies (class-wide)
    TAGS  = [ NORMAL ]

    #: debug mode -- don't delete temporary data
    DEBUG = False

    #: File handle for system-wide test log
    TESTLOG = sys.stdout
    if logging.getLogger().handlers:
        TESTLOG = logging.getLogger().handlers[0]

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
        self.log =  LogHandle(getattr( self, 'log', BiskitTest.TESTLOG ))
##      self.verbosity = getattr( self, 'verbosity', BiskitTest.VERBOSITY )
##      self.debugging = getattr( self, 'debugging', BiskitTest.DEBUG )
        self.prepare()

    def tearDown( self ):
        if not self.DEBUG:
            self.cleanUp()


def isTestClass( c ):
    """
    Checks whether the given class is derrived from the locally defined 
    TestCase extension. issubclass() would not work here because the class
    is defined in the same module.
    isTestClass( class ) <==> True if c is derrived from BiskitTest
    """
    if not issubclass( c, U.TestCase ):
        return False
    
    if c.__base__ and c.__base__.__name__ == 'BiskitTest':
        return True
    
    if c.__base__:
        return isTestClass( c.__base__ )
    
    return False

def isTestInstance( o ):
    return isTestClass( o.__class__ )
            
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
        :param tests: iterable of TestCases
        :type  tests: ( BiskitTest, )
        :param allowed: list of allowed tags
        :type  allowed: [ int ]
        :param forbidden : list of forbidden tags
        :type  forbidden : [ int ]
        """
        self._allowed   = allowed
        self._forbidden = forbidden
        U.TestSuite.__init__( self, tests=tests )


    def addTest( self, test ):
        """
        Add test case if it is matching the allowed and disallowed
        groups.
        :param test: test case
        :type  test: BiskitTest
        """
        assert isTestInstance( test ), \
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
        self.startclock = time.time()

    def addSuccess(self, test):
        ## super(U.TextTestResult, self).addSuccess(test)
        dt = time.time() - self.startclock
        
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

    def __init__( self, log=sys.stdout,
                  allowed=[], forbidden=[], verbosity=2, debug=False ):
        """
        :param log: log output target [default: STDOUT]
        :type  log: open file handle
        :param allowed: tags required for test to be considered, default: []
        :type  allowed: [ int ]
        :param forbidden: tags leading to the exclusion of test, default: []
        :type  forbidden: [ int ]
        :param verbosity: verbosity level for unittest.TextTestRunner
        :type  verbosity: int
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


    def modulesFromPath( self, path=packageRoot(), module='' ):
        """
        Import all python files of a package as modules. Sub-packages
        are ignored and have to be collected separately.

        :param path:  single search path for a package
        :type  path:  str
        :param module: name of the python package
        :type  module: str

        :return: list of imported python modules, see also __import__
        :rtype : [ module ]

        :raise ImportError: if a python file cannot be imported
        """
        module_folder = module.replace('.', os.path.sep)
        
        files = glob.glob( os.path.join( path, module_folder,'*.py' ) )

        files = map( stripFilename, files )
        files = [ f for f in files if f[0] != '_' ]
        
        r = []

        for f in files:
            try:
                r += [ __import__( '.'.join([module, f]), globals(),
                                   None, [module]) ]
            except Exception as why:
                if self.verbosity > 1:
                    logging.error( 'Import failure in %s: %r' % (f,why) )

        return r


    def addTestsFromModules( self, modules ):
        """
        Extract all test cases from a list of python modules and add them to
        the internal test suite.
        :param modules: list of modules to be checked for BiskitTest classes
        :type  modules: [ module ]
        """
        for m in modules:
            tested = 0

            for i in m.__dict__.values():
                
                if type(i) is type and isTestClass( i ):
                    
                    suite = U.defaultTestLoader.loadTestsFromTestCase( i )
                    self.suite.addTests( suite )
                    tested = 1

            if tested:
                self.modules_tested += [m]
            else:
                self.modules_untested += [m]


    def collectTests( self, path=packageRoot(), module='' ):
        """
        Add all BiskitTests found in a given module to the internal test suite.
        :param path:  single search path for a package
        :type  path:  str
        :param module: name of the python package
        :type  module: str
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
        print('\nThe test log file has been saved to: %r'% self.log.name)
        total  = self.result.testsRun
        failed = len(self.result.failures) + len(self.result.errors)

        m_tested  = len( self.modules_tested )
        m_untested= len( self.modules_untested)
        m_total  = m_tested + m_untested

        ## report coverage
        print('\nTest Coverage:\n=============\n')
        print('%i out of %i modules had no test case:' % (m_untested,m_total))
        for m in self.__moduleNames( self.modules_untested) :
            print('\t', m)

        ## print a summary
        print('\nSUMMARY:\n=======\n')
        print('A total of %i tests from %i modules were run.' %(total,m_tested))
        print('   - %i passed'% (total - failed))
        print('   - %i failed'% failed)

        ## and a better message about which module failed 
        if failed:

            for test, ftrace in self.result.failures:
                print('      - failed: %s'% test.id())

            for test, ftrace in self.result.errors:
                print('      - error : %s'% test.id())


    def run( self, dry=False ):
        """
        :param dry: do not actually run the test but just set it up [False]
        :type  dry: bool
        """
        ## push global settings into test classes
        for testclass in self.suite:
            testclass.DEBUG = self.debugging
            testclass.VERBOSITY = self.verbosity
            testclass.TESTLOG = self.log

        runner = SimpleTextTestRunner(self.log, verbosity=self.verbosity,
                                  descriptions=False)
        if not dry:
            self.result = runner.run( self.suite )

#########################
### Helper functions ####


def getOuterNamespace():
    """
    Fetch the namespace of the module/script running as __main__.
    :return: the namespace of the outermost calling stack frame
    :rtype: dict
    """
    import inspect

    try:
        frames = inspect.stack()
        f = frames[-1][0]   ## isolate outer-most calling frame
        r = f.f_globals
    finally:
        del frames, f

    return r


def extractTestCases( namespace ):
    """
    :return: all BisktTest child classes found in given namespace
    :rtype: [ class ]
    :raise BiskitTestError: if there is no BiskitTest child
    """
    r =[]

    for i in namespace.values():

        if type(i) is type and isTestClass( i ):

            r += [i]

    if not r:
        raise BiskitTestError('no BiskitTest class found in namespace')

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
    
    :param testclass: BiskitTest-derived class [default: first one found]
    :type  testclass: class
    :param verbosity: verbosity level of TextTestRunner
    :type  verbosity: int
    :param debug: don't delete temporary files (skipp cleanUp) [0]
    :type  debug: int

    :return: the test result object
    :rtype:  unittest.TestResult

    :raise BiskitTestError: if there is no BiskitTest-derived class defined
    """
    ## get calling namespace
    outer = getOuterNamespace()
    if testclass:
        testclasses = [testclass]
    else:
        testclasses = extractTestCases( outer )

    suite = U.TestSuite()
    for c in testclasses:
        suite.addTests( U.TestLoader().loadTestsFromTestCase(c) )

    all_tests = []  ## cache pointer to tests; suite will drop them after run
    for test in suite:
        test.DEBUG = debug
        test.VERBOSITY = verbosity
        test.TESTLOG = log
        all_tests += [ test ]

    runner= U.TextTestRunner(verbosity=verbosity)
    r = runner.run( suite )

    for t in all_tests:
        outer.update( t.__dict__ )
        outer.update( {'self':t })

    return r

###############
## Mock test ##

class Test(BiskitTest):
    """Mock test, test doesn't test itself"""
    pass

################################
### Script-related functions ###

def _use( defaults ):
    print("""
Run unittest tests for the whole package.

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

    p    - packages to test, e.g. mypackage mypackage.parser           [All]
    v    - int, verbosity level, 3 switches on several graphical plots      [2]
    log  - path to logfile (overriden); empty -log means STDOUT        [STDOUT]
    nox  - suppress test plots                                          [False]
    dry  - do not actually run the test but just collect tests          [False]
               
Examples:

    * Run all but long or obsolete tests from mypackage and mypackage.parser:
    test.py -e old long  -p mypackage mypackage.parser

    * Run only PVM-dependent tests of the mypackage.calc sub-package:
    test.py -i pvm -p mypackage.calc

        
Default options:
""")
    for key, value in defaults.items():
        print("\t-",key, "\t",value)
        
    sys.exit(0)
    
## quick and dirty command line argument parsing from Biskit.tools
def get_cmdDict(lst_cmd, dic_default):
    """
    Parse commandline options into dictionary of type C{ {<option> : <value>} }
    Options are recognised by a leading '-'.
    Error handling should be improved.
    
    Option C{ -x |file_name| } is interpreted as file with additional options.
    The key value pairs in lst_cmd replace key value pairs in the
    -x file and in dic_default.
    
    :param lst_cmd: list with the command line options::
                    e.g. ['-pdb', 'in1.pdb', 'in2.pdb', '-o', 'out.dat']
    :type  lst_cmd: [str]
    :param dic_default: dictionary with default options::
                        e.g. {'psf':'in.psf'}
    :type  dic_default: {str : str}

    :return: command dictionary::
             ala {'pdb':['in1.pdb', 'in2.pdb'], 'psf':'in.psf', 'o':'out.dat'}
    :rtype: {<option> : <value>}
    """
    dic_cmd = {}                     # return value
    try:

        for cmd in lst_cmd:
            if (cmd[0] == '-'):               # this entry is new option
                current_option = cmd[1:]      # take all but leading "-"
                dic_cmd[current_option] = ""  # make sure key exists even
                                              # w/o value
                counter = 0        # number of values for this option
            else:                  # this entry is value for latest option

                if counter < 1:
                    dic_cmd[current_option] = cmd

    # in case, several values follow after a "-xxx" option convert dictionary
    # entry into list and add all elements (until the next "-") to this list
                else:
                    if counter == 1:   # there is already a value assigned
    # convert to list
                        dic_cmd[current_option] = [dic_cmd[current_option]]
    # add value to list
                    dic_cmd[current_option] = dic_cmd[current_option] + [cmd]

                counter = counter + 1

    except (KeyError, UnboundLocalError) as why:
        raise UtilError("Can't resolve command line options.\n \tError:"+\
                  str(why))

    ## get extra options from external file
    try:
        if 'x' in dic_cmd:
            d = file2dic( dic_cmd['x'] )
            d.update( dic_cmd )
            dic_cmd = d
    except IOError:
        raise IOError("Error opening %s."% dic_cmd['x'])

    ## fill in missing default values
    dic_default.update( dic_cmd )
    dic_cmd = dic_default

    return dic_cmd

def cmdDict( defaultDic={} ):
    return get_cmdDict( sys.argv[1:], defaultDic )

  
def toList( o ):
    """toList(o) -> [o], or o,  if o is already a list"""
    if type( o ) != type( [] ):
        return [ o ]
    return o

def _str2tags( s ):
    """convert list of string options to list of valid TAGS"""
    try:
        r = [ x.upper() for x in s if x ] ## to list of uppercase str
        r = [ eval( x ) for x in r ]      ## to list of int
    except:
        logging.getLogger().error('unrecognized tags: %r'%s)

    return r

def _convertOptions( o ):
    o['i'] = _str2tags( toList( o['i'] ) )
    o['e'] = _str2tags( toList( o['e'] ) )
    o['v'] = int( o['v'] )
    o['nox'] = ('nox' in o)
    o['dry'] = ('dry' in o)
    o['debug'] = ('debug' in o)
    if o['log']:
        o['log'] = open( o['log'], 'a' )
    else:
        o['log'] = BiskitTest.TESTLOG
    o['p'] = toList( o['p'] )

############
### MAIN ###

if __name__ == '__main__':

    defaults = {'i':'',
                'e': DEFAULT_EXCLUDE,
                'p': DEFAULT_PACKAGES,
                'v':'2',
                'log': '', ##T.testRoot()+'/test.log',
                }

    o = cmdDict( defaults )

    if len( sys.argv ) == 1 and 'testing.py' in sys.argv[0]:
        _use( defaults )
        
    _convertOptions( o )

    BiskitTest.VERBOSITY = o['v']
    BiskitTest.DEBUG = o['debug']
    
    l = BiskitTestLoader( allowed=o['i'], forbidden=o['e'],
                          verbosity=o['v'], log=o['log'], debug=o['debug'])


    for package in o['p']:
        print('collecting ', repr( package ))
        l.collectTests( module=package )

    l.run( dry=o['dry'] )
    l.report()

    print("DONE")
    
    ## Travis-CI and others expect non-0 return on any error
    if l.result.wasSuccessful() and l.result.testsRun > 1:
        sys.exit(0)
    sys.exit(1)
