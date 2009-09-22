##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2009 Raik Gruenberg & Johan Leckner
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
##
##
## last $Author$
## last $Date$
## $Revision$

"""
Class for calls to external programs.
"""

import tempfile, os, time, subprocess

import Biskit.tools as t
import Biskit.settings as s
from Biskit.LogFile import StdLog
from Biskit.Errors import BiskitError
from Biskit.ExeConfigCache import ExeConfigCache
import Biskit as B

class RunError( BiskitError ):
    pass

class TemplateError( BiskitError ):
    pass

class Executor:
    """
    All calls of external programs should be done via this class or subclasses.

    Executor gets the necessary information about a program (binary,
    environment variables, etc) from ExeConfigCache, creates an input file
    or pipe from a template (if available) or an existing file, wrapps the
    program call into ssh and nice (if necessary), spawns an external process
    via subprocess.Popen, communicates the input file or string, waits for
    completion and collects the output file or string, and cleans up
    temporary files.

    There are two ways of using Executor
    ====================================
    
      1. (recommended) Create a subclass of Executor for a certain program
          call. Methods to override would be:

           - __init__  ... to set your own default values
                           (call parent __init__!)
           - prepare   ... called BEFORE program execution
           - finish    ... called AFTER successful program execution
           - isfailed  ... to detect the success status after program execution
           - fail      ... called if execution fails
           - cleanup   ... called AFTER everything else
                           (call parent cleanup!)
           
          Additionally, you should provide a simple program configuration file
          in biskit/Biskit/data/defaults/. See L{Biskit.ExeConfig} for
          details and examples!

      2.  Use Executor directly.
          An example is given in the __main__ section of this module.
          You first have to create an Executor instance with all the
          parameters, then call its run() method and collect the result.

          In the most simple cases this can be combined into one line:
         
          >>> out, error, returncode = Executor('ls', strict=0).run()

          strict=0 means, ExeConfig does not insist on an existing exe_ls.dat
          configuration file and instead looks for a program called 'ls' in the 
          search path.


    Templates
    =========
      Templates are files or strings that contain place holders like,
      for example::

          file_in=%(f_in)s
          file_out=%(f_out)s      

      At run time, Executor will create an input file or pipe from the
      template by replacing all place holders with values from its own
      fields. Let's assume, the above example is put into a file 'in.template'.

      >>> x = Executor( 'ls', template='in.template', f_in='in.dat')

      ... will then pass the following input to the ls program::
    
          file_in=in.dat
          file_out=/tmp/tmp1HYOvO

      However, the following input template will raise an error::
    
          file_in=%(f_in)s
          seed=%(seed)i
    
      ...because Executor doesn't have a 'seed' field. You could provide
      one by overwriting Executor.__init__. Alternatively, you can
      provide seed as a keyword to the original Executor.__init__:

      >>> x = Executor('ls', template='in.template',f_in='in.dat', seed=1.5)

      This works because Executor.__init__ puts all unknown key=value pairs
      into the object's name space and passes them on to the template.


    Communicating Input
    ===================

    Programs often expect scripts, commands or additional parameters
    from StdIn or from input files. Executor tries to support many
    scenarios -- which one is chosen mainly depends on the
    L{ExeConfig} `pipes` setting in exe_<program>.dat and on the
    `template` parameter given to Executor.__init__.  (Note: Executor
    loads the ExeConfig instance for the given program into its
    `self.exe` field.)

    Here is an overview over the different scenarios and how to
    activate them:

      1. B{ no input (default behaviour)}

        The program only needs command line parameters

        Condition:

          - template == None

      2. B{ input pipe from STDIN
        (== ``myprogram | 'some input string'``) }

        Condition:

          - exe.pipes == 1 / True
          - template != None ((or f_in points to existing file))

        Setup:

             1. `template` points to an existing file:

                 Executor reads the template file, completes it in
                 memory, and pushes it directly to the program.

             2. `template` points to string that doesn't look like a file name:

                 Executor completes the string in memory (using
                 `self.template % self.__dict__`) and pushes it
                 directly to the program. This is the fastest option
                 as it avoids file access alltogether.

             3. `template` == None but f_in points to an *existing* file:

                  Executor will read this file and push it unmodified to
                  the program via StdIn. (kind of an exception, if used at
                  all, f_in usual points to a *non-existing* file that
                  will receive the completed input file.)

      3. B{ input from file
        (== ``myprogram < input_file``) }

        Condition:

            - exe.pipes == 0 / False
            - template != None
            - push_inp == 1 / True (default)

        Setup:
          
           1. `template` points to an existing file:

                 Executor reads the template file, completes it in
                 memory, saves the completed file to disc (creating or
                 overriding self.f_in), opens the file and passes the
                 file handle to the program (instead of STDIN).

           2. `template` points to string that doesn't look like a file name:

                 Same as 3.1, except that the template is not read
                 from disc but directly taken from memory (see 2.2).

      4. B{ input file passed as argument to the program
        (== ``myprogram input_file``) }

        Condition:

          - exe.pipes == 0 / False

        Here it is up to you to provide the correct program
        argument.

        Setup:

           1. Use template completion:

                 The best option would be to set an explicit file name
                 for `f_in` and include this file name into  `args`, Example::

                   exe = ExeConfigCache.get('myprogram')
                   assert not exe.pipes 

                   x = Executor( 'myprogram', args='input.in', f_in='input.in',
                             template='/somewhere/input.template', cwd='/tmp' )

                 Executor creates your input file on the fly which is then
                 passed as first argument.

           2. Without template completion:

                 Similar, just that you don't give a template::

                   x = Executor( 'myprogram', args='input.in', f_in='input.in',
                             cwd='/tmp' )

                 It would then be up to you to provide the correct
                 input file in `/tmp/input.in`. You could override the
                 L{prepare()} hook method for creating it.

        There are other ways of doing the same thing.


    Look at L{generateInp()} to see what is actually going on. 


    References
    ==========

      - See also L{Biskit.IcmCad} and L{Biskit.Xplorer} for examples of
         how to overwrite and use Executor.

      - See also L{Biskit.ExeConfig} for a description of program
         configuration.
    """

    def __init__( self, name, args='', template=None, f_in=None, f_out=None,
                  f_err=None, strict=1, catch_out=1, push_inp=1, catch_err=0,
                  node=None, nice=0, cwd=None, log=None, debug=0,
                  verbose=None, validate=1, **kw ):

        """
        Create Executor. *name* must point to an existing program configuration
        unless *strict*=0. Executor will create a program input from
        the template and its own fields and put it into f_in. If f_in but
        no template is given, the unchanged f_in is used as input. If neither
        is given, the program is called without input. If a node is given,
        the process is wrapped in a ssh call. If *nice* != 0, the process
        is preceeded by nice. *cwd* specifies the working directory. By
        default, this setting is taken from the configuration file which
        defaults to the current working directory.

        @param name: program name (configured in .biskit/exe_name.dat)
        @type  name: str
        @param args: command line arguments
        @type  args: str
        @param template: template for input file -- this can be the template
                         itself or the path to a file containing it
                         (default: None)
        @type  template: str
        @param f_in: target for completed input file (default: None, discard)
        @type  f_in: str
        @param f_out: target file for program output (default: None, discard)
        @type  f_out: str
        @param f_err: target file for error messages (default: None, discard)
        @type  f_err: str
        @param strict: strict check of environment and configuration file
                       (default: 1)
        @type  strict: 1|0
        @param catch_out: catch output in file (f_out or temporary)
                          (default: 1)
        @type  catch_out: 1|0
        @param catch_err: catch errors in file (f_out or temporary)
                          (default: 1)
        @type  catch_err: 1|0
        @param push_inp: push input file to process via stdin ('< f_in') [1]
        @type  push_inp: 1|0
        @param node: host for calculation (None->no ssh) (default: None)
        @type  node: str
        @param nice: nice level (default: 0)
        @type  nice: int
        @param cwd: working directory, overwrites ExeConfig.cwd (default: None)
        @type  cwd: str
        @param log: execution log (None->STOUT) (default: None)
        @type  log: Biskit.LogFile
        @param debug: keep all temporary files (default: 0)
        @type  debug: 0|1
        @param verbose: print progress messages to log (default: log != STDOUT)
        @type  verbose: 0|1
        @param validate: validate binary and environment immedeatly (1=default)
                         or only before execution (0)
        @type validate: 1|0 
        @param kw: key=value pairs with values for template file or string
        @type  kw: key=value
        
        @raise ExeConfigError: if environment is not fit for running
                               the program
        """
        self.exe = ExeConfigCache.get( name, strict=strict )
        if validate:
            self.exe.validate()

        self.f_out = t.absfile( f_out )
        if not f_out and catch_out:
            self.f_out = tempfile.mktemp( '.out', name+'_' )

        self.f_err = t.absfile( f_err )
        if not f_err and catch_err:
            self.f_err = tempfile.mktemp( '.err', name+'_' )
                
        self.keep_out = f_out is not None
        self.keep_inp = f_in  is not None  #: do not clean up the input file
        self.catch_out = catch_out         #: capture STDOUT into file
        self.catch_err = catch_err         #: capture STDERR into file
        
        self.f_in = t.absfile( f_in )
        ## pre-define name of input file that will be generated later 
        if template and not f_in:
            self.f_in  = tempfile.mktemp('.inp', name+'_')

        self.push_inp = push_inp

        self.args = args
        self.template = template

        self.node  = node ## or os.uname()[1]
        self.nice  = nice
        self.debug = debug

        self.cwd = cwd or self.exe.cwd

        #: Log object for own messages
        self.log = log or StdLog()
        self.verbose = verbose
        if self.verbose is None:
            self.verbose = (log is not None)

        ## these are set by self.run():
        self.runTime = 0    #: time needed for last run
        self.output = None  #: STDOUT returned by process
        self.error = None   #: STDERR returned by process
        self.returncode = None #: int status returned by process
        self.pid = None     #: process ID

        self.result = None  #: set by self.finish()

        self.initVersion = self.version()

        self.__dict__.update( kw )


    def version( self ):
        """Version of class (at creation).
        @return: version
        @rtype: str
        """       
        return 'Executor $Revision$'


    def communicate( self, cmd, inp, bufsize=-1, executable=None,
                     stdin=None, stdout=None, stderr=None,
                     shell=0, env=None, cwd=None ):
        """
        Start and communicate with the new process. Called by execute().
        See subprocess.Popen() for a detailed description of the parameters!
        This method should work for pretty much any purpose but may fail for
        very long pipes (more than 100000 lines).
        
        @param inp: (for pipes) input sequence
        @type  inp: str
        @param cmd: command
        @type  cmd: str
        @param bufsize: see subprocess.Popen() (default: -1)
        @type  bufsize: int
        @param executable: see subprocess.Popen() (default: None)
        @type  executable: str
        @param stdin: subprocess.PIPE or file handle or None (default: None)
        @type  stdin: int|file|None
        @param stdout: subprocess.PIPE or file handle or None (default: None)
        @type  stdout: int|file|None
        @param stderr: subprocess.PIPE or file handle or None (default: None)
        @type  stderr: int|file|None
        @param shell: wrap process in shell; see subprocess.Popen()
                      (default: 0, use exe_*.dat configuration) 
        @type  shell: 1|0
        @param env: environment variables (default: None, use exe_*.dat config)
        @type  env: {str:str}
        @param cwd: working directory (default: None, means self.cwd)
        @type  cwd: str
        
        @return: output and error output
        @rtype: str, str
        
        @raise RunError: if OSError occurs during Popen or Popen.communicate
        """
        try:
            p = subprocess.Popen( cmd.split(),
                                  bufsize=bufsize, executable=executable,
                                  stdin=stdin, stdout=stdout, stderr=stderr,
                                  shell=shell or self.exe.shell,
                                  env=env or self.environment(), 
                                  cwd=cwd or self.cwd )

            self.pid = p.pid

            output, error = p.communicate( inp )
        
            self.returncode = p.returncode

        except OSError, e:
            raise RunError, \
                  "Couldn't run or communicate with external program: %r"\
                  % e.strerror

        return output, error


    def execute( self, inp=None ):
        """
        Run external command and block until it is finished.
        Called by L{ run() }.
        
        @param inp: input to be communicated via STDIN pipe (default: None)
        @type  inp: str

        @return: execution time in seconds
        @rtype: int
        
        @raise RunError: see communicate()
        """
        self.exe.validate() ##Check that binary and env variables are available
        
        start_time = time.time()

        cmd = self.command()
        
        shellexe = None
        if self.exe.shell and self.exe.shellexe:
            shellexe = self.exe.shellexe

        stdin = stdout = stderr = None

        if self.exe.pipes:
            stdin = subprocess.PIPE
            stdout= subprocess.PIPE
            stderr= subprocess.PIPE
        else:
            inp = None
            if self.f_in and self.push_inp:
                stdin = open( self.f_in )
            if self.f_out and self.catch_out:
                stdout= open( self.f_out, 'w' )
            if self.f_err and self.catch_err:    
                stderr= open( self.f_err, 'w' )
        
        if self.verbose:
            self.log.add('executing: %s' % cmd)
            self.log.add('in folder: %s' % self.cwd ) 
            self.log.add('input:  %r' % stdin )
            self.log.add('output: %r' % stdout )
            self.log.add('errors: %r' % stderr )
            self.log.add('wrapped: %r'% self.exe.shell )
            self.log.add('shell: %r'  % shellexe )
            self.log.add('environment: %r' % self.environment() )
            if self.exe.pipes and inp:
                self.log.add('%i byte of input pipe' % len(str(inp)))

        self.output, self.error = self.communicate( cmd, inp,
                            bufsize=-1, executable=shellexe, stdin=stdin,
                            stdout=stdout, stderr=stderr,
                            shell=self.exe.shell,
                            env=self.environment(), cwd=self.cwd )

        if self.exe.pipes and self.f_out:
            open( self.f_out, 'w').writelines( self.output )

        if self.verbose: self.log.add(".. finished.")

        return time.time() - start_time


    def run( self, inp_mirror=None ):
        """
        Run the callculation. This calls (in that order):
          - L{ prepare() },
          - L{ execute() },
          - L{ postProcess() },
          - L{ finish() } OR L{ fail() },
          - L{ cleanup() }
        
        @param inp_mirror: file name for formatted copy of inp file
                           (default: None) [not implemented]
        @type  inp_mirror: str

        @return: calculation result
        @rtype: any
        """
        try:
            self.prepare()

            self.inp  = self.generateInp()

            self.runTime = self.execute( inp=self.inp )

            self.postProcess()

        except MemoryError, why:
            try:
                self.fail()
            finally:
                self.cleanup()
            raise RunError, why

        try:
            if self.isFailed():
                self.fail()
            else:
                self.finish()
        finally:
            self.cleanup()

        return self.result


    def command( self ):
        """
        Compose command string from binary, arguments, nice, and node.
        Override (perhaps).
        
        @return: the command to execute
        @rtype: str
        """
        exe = t.absbinary( self.exe.bin )

        if self.args:
            exe = exe + ' ' + self.args

        str_nice = str_ssh = ''

        if self.nice != 0:
            str_nice = "%s -%i" % (s.nice_bin, self.nice)

        if self.node is not None:
            str_ssh  = "%s %s"  % (s.ssh_bin, self.node )

        cmd = "%s %s %s" % (str_ssh, str_nice, exe )
        cmd = cmd.strip()

        return cmd


    def environment( self ):
        """
        Setup the environment for the process. Override if needed.
        
        @return: environment dictionary
        @rtype: dict OR None
        """
        if not self.exe.replaceEnv:
            return None

        return self.exe.environment()


    def prepare( self ):
        """
        called before running external program, override!
        """
        pass


    def postProcess( self ):
        """
        called directly after running the external program, override!
        """
        pass


    def cleanup( self ):
        """
        Clean up after external program has finished (failed or not).
        Override, but call in child method!
        """
        if not self.keep_out and not self.debug and self.f_out:
            t.tryRemove( self.f_out )

        if not self.keep_inp and not self.debug:
            t.tryRemove( self.f_in )

        if self.f_err and not self.debug:
            t.tryRemove( self.f_err )


    def fail( self ):
        """
        Called if external program failed, override! 
        """
        B.EHandler.warning(\
            'Execution failed. (Override Executor.fail() to handle this!)')


    def finish( self ):
        """
        Called if external program finished successfully, override!
        """
        self.result = self.output, self.error, self.returncode


    def isFailed( self ):
        """
        Detect whether external program failed, override!
        """
        return 0


    def fillTemplate( self ):
        """
        Create complete input string from template with place holders.
        
        @return: input
        @rtype: str
        
        @raise TemplateError: if unknown option/place holder in template file
        """
        inp = self.template

        try:
            if os.path.isfile( inp ):
                inp = open( inp, 'r' ).read()
            return inp % self.__dict__

        except KeyError, why:
            s =  "Unknown option/place holder in template file."
            s += "\n  template file: " + str( self.template )
            s += "\n  Template asked for a option called " + str( why[0] )
            raise TemplateError, s


    def convertInput( self, inp):
        """
        Convert the input to a format used by the selected execution method.
        
        @param inp: path to existing input file or string with input
        @type  inp: str
        
        @return: input string if self.exe.pipes; file name otherwise
        @rtype: str
        """
        if self.exe.pipes:

            ## convert file to string
            if not inp and os.path.exists( self.f_in or '' ):
                return open( self.f_in, 'r' ).read()
            ## or return string unchanged
            return inp

        ## no pipes and no input string
        if inp is None:

            return inp

        ## else put input string into file
        assert self.f_in

        f = open( self.f_in, 'w')
        f.write(inp)
        f.close()
        return self.f_in


    def generateInp(self):
        """
        Prepare the program input (file or string) from a template (if
        any, file or string).

        @return: input file name OR (if pipes=1) content of input file OR None
        @rtype: str
        
        @raise TemplateError: if error while creating template file
        """
        try:
            inp = None

            if self.template:
                inp = self.fillTemplate()

            return self.convertInput( inp )

        except Exception, why:
            s =  "Error while creating input file from template."
            s += "\n  template file: " + str( self.template )
            s += "\n  why: " + str( why )
            s += "\n  Error:\n  " + t.lastError()
            raise TemplateError, s


#############
##  TESTING        
#############
import Biskit.test as BT

class Test(BT.BiskitTest):
    """Executor test"""

    TAGS = [ BT.EXE ]

    def prepare(self):
        import tempfile
        self.fout = tempfile.mktemp('_testexecutor.out')

    def cleanUp(self):
        t.tryRemove(self.fout)
    
    def test_Executor( self ):
        """Executor test (run emacs ~/.biskit/settings.cfg)"""
        ExeConfigCache.reset()

        self.x = ExeConfigCache.get( 'emacs', strict=0 )
        self.x.pipes = 1

        args = '.biskit/settings.cfg'
        if not self.local:
            args = '-kill ' + args

        self.e = Executor( 'emacs', args=args, strict=0,
                           f_in=None,
                           f_out=self.fout,
                           verbose=self.local, cwd=t.absfile('~') )
        
        self.r = self.e.run()

        if self.local:
            print 'Emacs was running for %.2f seconds'%self.e.runTime

        self.assert_( self.e.pid != None )

if __name__ == '__main__':

    BT.localTest()


