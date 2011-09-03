#!/usr/bin/env python
##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2005 Raik Gruenberg & Johan Leckner
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License as
## published by the Free Software Foundation; either version 3 of the
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

import os

import tempfile

from Biskit.tools import *


o = {'bin': '/usr/local/bin/epydoc',
     't':projectRoot()+'/Biskit',
     'url':'http://biskit.pasteur.fr',
     'o':projectRoot()+'/docs/modules',
     'ver':3}


def _use():
    print """
    Create epyDoc HTML module documentation for Biskit

    bin  - str, path to epydoc bineries
    ver  - 2|3, epydoc version 2.x or 3.x
    t    - str, target folder/file
    url  - str, target url for top right link
    o    - str, output folder

Default options:
"""
    for key, value in o.items():
        print "\t-",key, "\t",value

    sys.exit(0)


GRAY = """
/* Body color */ 
body               { background: #ffffff; color: #000000; } 

/* Tables */ 
table.summary, table.details, table.index
                   { background: #dddddd; color: #000000; } 
tr.summary, tr.details, tr.index
                   { background: #999999; color: #000000;  
                     text-align: left; font-size: 120%; } 
tr.group           { background: #c0e0f8; color: #000000;
                     text-align: left; font-size: 120%;
                     font-style: italic; } 

/* Documentation page titles */
h2.module          { margin-top: 0.2em; }
h2.class           { margin-top: 0.2em; }

/* Headings */
h1.heading         { font-size: +140%; font-style: italic;
                     font-weight: bold; }
h2.heading         { font-size: +125%; font-style: italic;
                     font-weight: bold; }
h3.heading         { font-size: +110%; font-style: italic;
                     font-weight: normal; }

/* Base tree */
pre.base-tree      { font-size: 80%; margin: 0; }

/* Details Sections */
table.func-details { background: #cccccc; color: #000000;
                     border: 2px groove #c0d0d0;
                     padding: 0 1em 0 1em; margin: 0.4em 0 0 0; }
h3.func-detail     { background: transparent; color: #000000;
                     margin: 0 0 1em 0; }

table.var-details  { background: #cccccc; color: #000000;
                     border: 2px groove #c0d0d0;
                     padding: 0 1em 0 1em; margin: 0.4em 0 0 0; }
h3.var-details     { background: transparent; color: #000000;
                     margin: 0 0 1em 0; }

/* Function signatures */
.sig               { background: transparent; color: #000000;
                     font-weight: bold; }  
.sig-name          { background: transparent; color: #006080; }  
.sig-arg, .sig-kwarg, .sig-vararg
                   { background: transparent; color: #008060; }  
.sig-default       { background: transparent; color: #602000; }  
.summary-sig       { background: transparent; color: #000000; }  
.summary-sig-name  { background: transparent; color: #204080; }
.summary-sig-arg, .summary-sig-kwarg, .summary-sig-vararg
                   { background: transparent; color: #008060; }  

/* Doctest blocks */
.py-src            { background: transparent; color: #000000; }
.py-prompt         { background: transparent; color: #005050;
                     font-weight: bold;}
.py-string         { background: transparent; color: #006030; }
.py-comment        { background: transparent; color: #003060; }
.py-keyword        { background: transparent; color: #600000; }
.py-output         { background: transparent; color: #404040; }
pre.doctestblock   { background: #f4faff; color: #000000; 
                     padding: .5em; margin: 1em;
                     border: 1px solid #708890; }
table pre.doctestblock
                   { background: #dce4ec; color: #000000; 
                     padding: .5em; margin: 1em;
                     border: 1px solid #708890; }

/* Variable values */
pre.variable       { background: #dce4ec; color: #000000;
                     padding: .5em; margin: 0;
                     border: 1px solid #708890; }
.variable-linewrap { background: transparent; color: #604000; }
.variable-ellipsis { background: transparent; color: #604000; }
.variable-quote    { background: transparent; color: #604000; }
.re                { background: transparent; color: #000000; }
.re-char           { background: transparent; color: #006030; }
.re-op             { background: transparent; color: #600000; }
.re-group          { background: transparent; color: #003060; }
.re-ref            { background: transparent; color: #404040; }

/* Navigation bar */ 
table.navbar       { background: #bbbbbb; color: #0000ff;
                     border: 2px groove #c0d0d0; }
th.navbar          { background: #bbbbbb; color: #0000ff; } 
th.navselect       { background: #999999; color: #000000; } 
.nomargin          { margin: 0; }

/* Links */ 
a:link             { background: transparent; color: #0000ff; }  
a:visited          { background: transparent; color: #204080; }  
a.navbar:link      { background: transparent; color: #0000ff; 
                     text-decoration: none; }  
a.navbar:visited   { background: transparent; color: #204080; 
                     text-decoration: none; }  

/* Lists */
ul { margin-top: 0; }
"""

##################################################
# main function
##################################################
if __name__ == '__main__':

    if sys.argv[-1].upper() in ['H', '-H', 'HELP', '-HELP']:   
        _use()
        sys.exit(0)
    else:

        fname = tempfile.mktemp('.css', 'epydoc')
        f = open( fname, 'w' )
        f.writelines( GRAY )
        f.close()

        o = cmdDict( o )
        o['gray'] = fname
        o['dot'] = absbinary('dot')

        ## command for version 2.x of epydoc
        if int(o['ver'])==2:
            command = '%(bin)s --css white --private-css %(gray)s -u %(url)s -o %(o)s -n Biskit -t "doc_modules.html" -v %(t)s'%o

        ## command for version 3.x of epyd
        ##
        ## 3.aplha2 argument --show-sourcecode now need to be given explicitly
        if int(o['ver'])==3:
            command =  '%(bin)s -o %(o)s --css=white --url=%(url)s --name=Biskit --top=indices.html --show-sourcecode --graph=classtree --dotpath=%(dot)s --parse-only --verbose %(t)s'%o

    #  classtree, callgraph, umlclasstree

        os.system( command )
        os.unlink(fname)


#--dotpath=/usr/bin/dot  www.graphviz.org/ --graph all    
#--show-imports

## sudo rm -r /usr/local/lib/python2.4/site-packages/epydoc
## sudo rm /usr/local/bin/epydoc*
## cd src/
## sudo python setup.py install

#### COMMENTS FOR EPYDOC3
## in /usr/local/lib/python2.4/random.py i had to comment out three lines of the test

## epydoc -o modules --css=white --url=http://biskit.sf.net --name=Biskit --top=index.html --show-sourcecode --graph=classtree --graph=umlclasstree --dotpath=/usr/bin/dot --verbose --parse-only /shared_bin/biskit/Bisk


#####################
## EpyDoc 3 options

## usage: epydoc ACTION [options] NAMES...

## options:
##   --version             show program's version number and exit
##   -h, --help            show this help message and exit

##   Actions:
##     --html              Write HTML output.
##     --text              Write plaintext output. (not implemented yet)
##     --latex             Write LaTeX output.
##     --dvi               Write DVI output.
##     --ps                Write Postscript output.
##     --pdf               Write PDF output.
##     --check             Check completeness of docs.
##     --pickle            Write the documentation to a pickle file.

##   Options:
##     -o PATH, --output=PATH
##                         The output directory.  If PATH does not exist, then it
##                         will be created.
##     --inheritance=STYLE
##                         The format for showing inheritance objects.  STYLE
##                         should be one of: grouped, listed, included.
##     --docformat=NAME    The default markup language for docstrings.  Defaults
##                         to "epytext".
##     --css=STYLESHEET    The CSS stylesheet.  STYLESHEET can be either a
##                         builtin stylesheet or the name of a CSS file.
##     --name=NAME         The documented project's name (for the navigation
##                         bar).
##     --url=URL           The documented project's URL (for the navigation bar).
##     --navlink=HTML      HTML code for a navigation link to place in the
##                         navigation bar.
##     --top=PAGE          The "top" page for the HTML documentation.  PAGE can
##                         be a URL, the name of a module or class, or one of the
##                         special names "trees.html", "indices.html", or
##                         "help.html"
##     --help-file=FILE    An alternate help file.  FILE should contain the body
##                         of an HTML file -- navigation bars will be added to
##                         it.
##     --show-frames       Include frames in the HTML output. (default)
##     --no-frames         Do not include frames in the HTML output.
##     --show-private      Include private variables in the output. (default)
##     --no-private        Do not include private variables in the output.
##     --show-imports      List each module's imports.
##     --no-imports        Do not list each module's imports. (default)
##     -q, --quiet         Decrease the verbosity.
##     -v, --verbose       Increase the verbosity.
##     --debug             Show full tracebacks for internal errors.
##     --parse-only        Get all information from parsing (don't introspect)
##     --introspect-only   Get all information from introspecting (don't parse)
##     --profile-epydoc    Run the hotshot profiler on epydoc itself.  Output
##                         will be written to profile.out.
##     --dotpath=PATH      The path to the Graphviz 'dot' executable.
##     --config=FILE       A configuration file, specifying additional OPTIONS
##                         and/or NAMES.  This option may be repeated.
##     --graph=GRAPHTYPE   Include graphs of type GRAPHTYPE in the generated
##                         output.  Graphs are generated using the Graphviz dot
##                         executable.  If this executable is not on the path,
##                         then use --dotpath to specify its location.  This
##                         option may be repeated to include multiple graph types
##                         in the output.  GRAPHTYPE should be one of: all,
##                         classtree, callgraph, umlclasstree.
##     --graph-font=FONT   Specify the font used to generate Graphviz graphs.
##                         (e.g., helvetica or times).
##     --graph-font-size=SIZE
##                         Specify the font size used to generate Graphviz
##                         graphs, in points.
##     --separate-classes  When generating LaTeX or PDF output, list each class
##                         in its own section, instead of listing them under
##                         their containing module.
##     --show-sourcecode   Include source code with syntax highlighting in the
##                         HTML output.
##     --no-sourcecode     Do not include source code with syntax highlighting in
##                         the HTML output.
##     --pstat=FILE        A pstat output file, to be used in generating call
##                         graphs.

