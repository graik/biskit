##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2005 Raik Gruenberg & Johan Leckner; All rights reserved
##

## last $Author$
## last $Date$
## $Revision$

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
    All calls of external programs should be done via this class or subclasses.

    Executor gets the necessary information about a program (binary,
    environment variables, etc) from ExeConfigCache, creates an input file
    or pipe from a template (if available) or an existing file, wrapps the
    program call into ssh and nice (if necessary), spawns an external process
    via subprocess.Popen, communicates the input file or string, waits for
    completion and collects the output file or string, and cleans up
    temporary files. There are two ways of using Executor:

    1) (recommended) Create a subclass of Executor for a certain program call.
    Methods to override would be:
    
    * __init__  ... to set your own default values (call parent __init__!)
    * prepare   ... called BEFORE program execution
    * cleanup   ... called AFTER program execution (call parent cleanup!)
    * finish    ... called AFTER successful program execution
    * isfailed  ... to detect the success status after program execution
    * failed    ... called if execution fails

    Additionally, you should provide a simple program configuration file in
    biskit/external/defaults/. See ExeConfig for details and examples!

    2) Use Executor directly.
    An example is given in the __main__ section of this module. You first have
    to create an Executor instance with all the parameters, then call its
    run() method and collect the result.

    In the most simple cases this can be combined into one line:
    out, error, returncode = Executor('ls', strict=0).run()

    strict=0 means, ExeConfig does not insist on an existing exe_ls.dat file
    and instead looks for a program called 'ls' in the search path.


    Templates:
    
    Templates are files or strings that contain place holders like,
    for example:

    '''
    file_in=%(f_in)s
    file_out=%(f_out)s
    '''

    At run time, Executor will create an input file or pipe from the
    template by replacing all place holders with values from its own
    fields. Let's assume, the above example is put into a file 'in.template'.
    
    x = Executor( 'ls', template='in.template', f_in='in.dat')

    ... will then pass the following input to the ls program:
    '''
    file_in=in.dat
    file_out=/tmp/tmp1HYOvO
    '''

    However, the following input template will raise an error:
    '''
    file_in=%(f_in)s
    seed=%(seed)i
    '''
    ...because Executor doesn't have a 'seed' field. You could provide
    one by overwriting Executor.__init__. Alternatively, you can
    provide seed as a keyword to the original Executor.__init__:

    x = Executor('ls', template='in.template',f_in='in.dat', seed=1.5)

    This works because Executor.__init__ puts all unknown key=value pairs
    into the object's name space and passes them on to the template.


    See also IcmCad for an Example of how to overwrite and use Executor.
    See also ExeConfig for a description of program configuration.
    """

    def __init__( self, name, args='', template=None, f_in=None, f_out=None,
                  strict=1, catch_out=1, push_inp=1, node=None, nice=0,
                  cwd=None,
                  log=None, debug=0, verbose=0, **kw ):
        
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
        
        name     - str, program name (configured in .biskit/exe_name.dat)
        args     - str, command line arguments                          []
        template - str, path to template for input file             [None]
        f_in     - str, path to complete input file      [None .. discard]
        f_out    - str, target file for output           [None .. discard]
        strict   - 1|0, strict check of environment and configuration file [1]
        catch_out- 1|0, catch output in file (f_out or temporary)      [1]
        push_inp - 1|0, push input file to process via stdin ('< f_in') [1] 
        node     - str, host for calculation (None->no ssh)         [None]
        nice     - int, nice level                                     [0]
        cwd      - str, working directory, overwrites ExeConfig.cwd [None]
        log      - Biskit.LogFile, program log (None->STOUT)        [None]
        debug    - 0|1, keep all temporary files                       [0]
        verbose  - 0|1, print progress messages to log     [log != STDOUT]
        **kw     - key=value pairs with values for template file
        !! ExeConfigError, if environment is not fit for running the program
        """
        self.exe = ExeConfigCache.get( name, strict=strict )
        self.exe.validate()
        
        self.f_out = t.absfile( f_out )
        if not f_out and catch_out:
            self.f_out = tempfile.mktemp( '.out' )
        self.keep_out = f_out is not None
        self.catch_out= catch_out

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
                     stdin=None, stdout=None, shell=0, env=None,
                     cwd=None ):
        """
        Start and communicate with the new process. Called by execute().
        See subprocess.Popen() for a detailed description of the parameters!
        This method should work for pretty much any purpose but may fail for
        very long pipes (more than 100000 lines).
        inp     - str, (for pipes) input sequence
        cmd     - str, command
        bufsize    - int, see subprocess.Popen()                           [-1]
        executable - str, see subprocess.Popen()                         [None]
        stdin      - int|file|None, PIPE or file handle or None          [None]
        stdout     - int|file|None, PIPE or file handle or None          [None]
        shell      - 1|0, wrap process in shell; see subprocess.Popen()     [0]
        env        - {str:str}, environment variables                    [None]
        cwd        - str, working directory                              [None]
        -> str, str - output and error output
        !! RunError, if OSError occurs during Popen or Popen.communicate
        """
        try:
            p = subprocess.Popen( cmd.split(),
                                  bufsize=bufsize, executable=executable,
                                  stdin=stdin, stdout=stdout,
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
        Run external command and block until it is finished. Called by run().
        inp  - str, input to be communicated via STDIN pipe    [None]
        !! RunError, see communicate()
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
            stderr= None

        if self.verbose:
            self.log.add('executing: %s' % cmd)
            self.log.add('in folder: %s' % self.cwd ) 
            self.log.add('input:  %r' % stdin )
            self.log.add('output: %r' % stdout )
            self.log.add('wrapped: %r'% self.exe.shell )
            self.log.add('shell: %r'  % shellexe )
            self.log.add('environment: %r' % self.environment() )
            if self.exe.pipes and inp:
                self.log.add('%i byte of input pipe' % len(str(inp)))

        self.output, self.error = self.communicate( cmd, inp,
                            bufsize=-1, executable=shellexe, stdin=stdin,
                            stdout=stdout, shell=self.exe.shell,
                            env=self.environment(), cwd=self.cwd )

        if self.exe.pipes and self.f_out:
            open( self.f_out, 'w').writelines( self.output )
        
        if self.verbose: self.log.add(".. finished.")

        return time.time() - start_time
       

    def run( self, inp_mirror=None ):
        """
        calls (in that order): prepare(), execute(),
        finish()/failed(), cleanup()
        inp_mirror - str, file name for formatted copy of inp file [None]
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
        -> str, the command to execute
        """
        exe  = t.absbinary( self.exe.bin )

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
        -> dict | None
        """
        if not self.exe.replaceEnv:
            return None
        
        return self.exe.environment()

    def prepare( self ):
        """called before running external program, override!"""
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


    def failed( self ):
        """Called if external program failed, override! """
        pass

    def finish( self ):
        """Called if external program finished successfully, override!"""
        self.result = self.output, self.error, self.returncode

    def isFailed( self ):
        """Detect whether external program failed, override!"""
        return 0

    def fillTemplate( self ):
        """
        Create complete input string from template with place holders.
        -> str, input
        !! TemplateError
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
        inp - str, path to existing input file or string with input
        -> str, input string if self.exe.pipes; file name otherwise
        """
        if self.exe.pipes:

            ## convert file to string
            if not inp and os.path.exists( self.f_in or '' ):    
                return open( self.f_in, 'r' ).read()

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

        -> input file name OR (if pipes=1) content of input file
        !! raise TemplateError
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


if __name__ == '__main__':

    ExeConfigCache.reset()

    x = ExeConfigCache.get( 'emacs', strict=0 )
    x.pipes = 1

##     e = Executor( 'emacs', args='.zshenv', strict=0, catch_out=0,
##                   verbose=1, cwd=t.absfile('~') )
    e = Executor( 'emacs', args='.zshenv', strict=0, node='magneto',
                  f_in=None,
                  f_out=t.absfile('~/test.out'),
                  verbose=1, cwd=t.absfile('~') )

    r = e.run()