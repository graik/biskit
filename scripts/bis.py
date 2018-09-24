#!/bin/sh
''''exec python3 -i -- "$0" ${1+"$@"} # '''

## above shebang hack calls python interpreter with interactive flag
## standard /usr/bin/env python -i doesn't get parsed in recent Linux 
## see: http://stackoverflow.com/questions/3306518/cannot-pass-an-argument-to-python-with-usr-bin-env-python

from biskit.tools import *
from biskit import *
import sys

## activate tab-completion
import rlcompleter, readline
readline.parse_and_bind('tab: complete')

def help():
    print("""
bispy [pickled_object or PDB_file]

Opens an interactive python prompt with all Biskit/* modules and Biskit.tools
functions imported into the default namespace. Tab-completion is enabled.

A single argument can point to: 
  - a pickled python object which is then unpickled and put into a variable x,
  - or a PDB file which is then loaded into a PDBModel object x
  - or a PDB ID to be fetched from the PDB database and loaded into PDBModel x

For example, all the following should work:
> bis.py mypickled_python.dic
> bis.py mystructure.pdb
> bis.py 3TGI
""")

args = sys.argv

help()

if len( args ) > 1:
    
    try:
        if (len(args[1]) > 3 and str.upper( args[1][-4:] ) == '.PDB') \
               or (len( args[1] ) ==4 ):
            x = PDBModel( args[1] )
        else:
            x = load( args[1] )
            
        print("Loaded %s and put it into variable x." % args[1])
        
    except:
        errWriteln( lastError() )
        errWriteln('%s does not seem to be a pickled python object.'\
                   % args[1] + '(and not a PDB file either.)')
