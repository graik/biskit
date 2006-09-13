##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2005 Raik Gruenberg & Johan Leckner
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

class RunError( BiskitError ):
    pass

class TemplateError( BiskitError ):
    pass

class Executor:
    """
    Executor
    ========
    
    All calls of external programs should be done via this class or subclasses.

    Executor gets the necessary information about a program (binary,
    environment variables, etc) from ExeConfigCache, creates an input file
    or pipe from a template (if available) or an existing file, wrapps the
    program call into ssh and nice (if necessary), spawns an external process
    via subprocess.Popen, communicates the input file or string, waits for
    completion and collects the output file or string, and cleans up
    temporary files.

    There are two ways of using Executor
    ------------------------------------
    
      1. (recommended) Create a subclass of Executor for a certain program
          call. Methods to override would be:

           - __init__  ... to set your own default values
                           (call parent __init__!)
           - prepare   ... called BEFORE program execution
           - cleanup   ... called AFTER program execution
                           (call parent cleanup!)
           - finish    ... called AFTER successful program execution
           - isfailed  ... to detect the success status after program execution
           - failed    ... called if execution fails
           
          Additionally, you should provide a simple program configuration file
          in biskit/external/defaults/. See L{Biskit.ExeConfig} for
          details and examples!

      2.  Use Executor directly.
          An example is given in the __main__ section of this module.
          You first have to create an Executor instance with all the
          parameters, then call its run() method and collect the result.

          In the most simple cases this can be combined into one line:
         
          >>> out, error, returncode = Executor('ls', strict=0).run()

          strict=0 means, ExeConfig does not insist on an existing exe_ls.dat
          file and instead looks for a program called 'ls' in the search path.


    Templates
    ---------
      Templates are files or strings that contain place holders like,
      for example:

      >>> file_in=%(f_in)s
      >>> file_out=%(f_out)s      

      At run time, Executor will create an input file or pipe from the
      template by replacing all place holders with values from its own
      fields. Let's assume, the above example is put into a file 'in.template'.

      >>> x = Executor( 'ls', template='in.template', f_in='in.dat')

      ... will then pass the following input to the ls program:
    
      >>> file_in=in.dat
      >>> file_out=/tmp/tmp1HYOvO

      However, the following input template will raise an error:
    
      >>> file_in=%(f_in)s
      >>> seed=%(seed)i
    
      ...because Executor doesn't have a 'seed' field. You could provide
      one by overwriting Executor.__init__. Alternatively, you can
      provide seed as a keyword to the original Executor.__init__:

      >>> x = Executor('ls', template='in.template',f_in='in.dat', seed=1.5)

      This works because Executor.__init__ puts all unknown key=value pairs
      into the object's name space and passes them on to the template.


    References
    ----------
      See also L{Biskit.IcmCad} for an Example of how to overwrite and
      use Executor.
      See also L{Biskit.ExeConfig} for a description of program configuration.
    """

    def __init__( self, name, args='', template=None, f_in=None, f_out=None,
                  f_err=None, strict=1, catch_out=1, push_inp=1, catch_err=0,
                  node=None, nice=0, cwd=None, log=None, debug=0,
                  verbose=0, **kw ):

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
        @param template: path to template for input file (default: None)
        @type  template: str
        @param f_in: path to complete input file (default: None, discard)
        @type  f_in: str
        @param f_out: target file for output (default: None, discard)
        @type  f_out: str
        @param strict: strict check of environment and configuration file
                       (default: 1)
        @type  strict: 1|0
        @param catch_out: catch output in file (f_out or temporary)
                          (default: 1)
        @type  catch_out: 1|0
        @param push_inp: push input file to process via stdin ('< f_in') [1]
        @type  push_inp: 1|0
        @param node: host for calculation (None->no ssh) (default: None)
        @type  node: str
        @param nice: nice level (default: 0)
        @type  nice: int
        @param cwd: working directory, overwrites ExeConfig.cwd (default: None)
        @type  cwd: str
        @param log: Biskit.LogFile, program log (None->STOUT) (default: None)
        @type  log: 
        @param debug: keep all temporary files (default: 0)
        @type  debug: 0|1
        @param verbose: print progress messages to log (log != STDOUT)
        @type  verbose: 0|1
        @param kw: key=value pairs with values for template file
        @type  kw: key=value
        
        @raise ExeConfigError: if environment is not fit for running
                               the program
        """
        self.exe = ExeConfigCache.get( name, strict=strict )
        self.exe.validate()

        self.f_out = t.absfile( f_out )
        if not f_out and catch_out:
            self.f_out = tempfile.mktemp( '.out' )

        self.f_err = t.absfile( f_err )
        if not f_err and catch_err:
            self.f_err = tempfile.mktemp( '.err' )
                
        self.keep_out  = f_out is not None
        self.catch_out = catch_out
        self.catch_err = catch_err
        
        self.f_in  = f_in  ## will be overridden by self.run()
        self.keep_inp = f_in is not None
        self.push_inp = push_inp

        self.args = args
        self.template = template

        self.node  = node ## or os.uname()[1]
        self.nice  = nice
        self.debug = debug

        self.cwd = cwd or self.exe.cwd

        ## Log object for own program messages
        self.log = log or StdLog()
        self.verbose = verbose
        if self.verbose is None:
            self.verbose = (log is not None)

        ## these are set by self.run():
        self.runTime = 0    ## time needed for last run
        self.output = None  ## STDOUT returned by process
        self.error = None   ## STDERR returned by process
        self.returncode = None ## int status returned by process

        self.result = None  ## set by self.finish()

        self.__dict__.update( kw )


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
        @param stdin: PIPE or file handle or None (default: None)
        @type  stdin: int|file|None
        @param stdout: PIPE or file handle or None (default: None)
        @type  stdout: int|file|None
        @param shell: wrap process in shell; see subprocess.Popen()
                      (default: 0) 
        @type  shell: 1|0
        @param env: environment variables (default: None)
        @type  env: {str:str}
        @param cwd: working directory (default: None)
        @type  cwd: str
        
        @return: output and error output
        @rtype: str, str
        
        @raise RunError: if OSError occurs during Popen or Popen.communicate
        """
        try:
            p = subprocess.Popen( cmd.split(),
                                  bufsize=bufsize, executable=executable,
                                  stdin=stdin, stdout=stdout, stderr=stderr,
                                  shell=self.exe.shell,
                                  env=self.environment(), cwd=self.cwd )

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
          - L{ finish() }/ L{ failed() },
          - L{ cleanup() }
        
        @param inp_mirror: file name for formatted copy of inp file
                           (default: None)
        @type  inp_mirror: str

        @return: calculation result
        @rtype: any
        """
        try:
            self.prepare()

            self.inp  = self.generateInp()

            self.runTime = self.execute( inp=self.inp )

        except RunError, why:
            try:
                self.failed()
            finally:
                self.cleanup()
            raise RunError, why

        try:
            if self.isFailed():
                self.failed()
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


    def failed( self ):
        """
        Called if external program failed, override! 
        """
        pass


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
            print '########', inp
            return inp

        ## no pipes and no input string
        if inp is None:

            return inp

        ## put input string into file
        self.f_in = self.f_in or tempfile.mktemp('_exec.inp')

        f = open( self.f_in, 'w')
        f.write(inp)
        f.close()
        return self.f_in


    def generateInp(self):
        """
        Replace formatstr place holders in inp by fields of this class.

        @return: input file name OR (if pipes=1) content of input file
        @rtype: str
        
        @raise TemplateError: if error while creating template file
        """
        try:
            inp = None

            if self.template:
                inp = self.fillTemplate()

            return self.convertInput( inp )

        except Exception, why:
            s =  "Error while creating template file."
            s += "\n  template file: " + str( self.template )
            s += "\n  why: " + str( why )
            s += "\n  Error:\n  " + t.lastError()
            raise TemplateError, s


#############
##  TESTING        
#############
        
class Test:
    """
    Test class
    """
    
    def run( self, local=0  ):
        """
        run function test

        @param local: transfer local variables to global and perform
                      other tasks only when run locally
        @type  local: 1|0
        
        @return: 1
        @rtype: int
        """
        ExeConfigCache.reset()

        x = ExeConfigCache.get( 'emacs', strict=0 )
        x.pipes = 1

        if local:
            e = Executor( 'emacs', args='.zshenv', strict=0,
                          f_in=None,
                          f_out=t.absfile('~/test.out'),
                          verbose=1, cwd=t.absfile('~') )
        else:
            e = Executor( 'emacs', args='-kill .zshenv', strict=0,
                          f_in=None,
                          f_out=t.absfile('~/test.out'),
                          verbose=0, cwd=t.absfile('~') )            
        
        r = e.run()

        if local:
            print 'Emacs were running for %.2f seconds'%e.runTime
            globals().update( locals() )
            
        return 1


    def expected_result( self ):
        """
        Precalculated result to check for consistent performance.

        @return: 1
        @rtype:  int
        """
        return 1
    

if __name__ == '__main__':

    test = Test()

    assert test.run( local=1 ) == test.expected_result()

    


