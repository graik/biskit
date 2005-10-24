## List of available cluster nodes
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
## Adapt this file and save it under biskit/Biskit/hosts.py !
##
## last $Date$
## last $Author$
"""
List of cluster computers.
Each host must be accessible via ssh w/o password.

nodes_* .. lists with one entry per computer
cpus_*  .. lists with one entry per CPU (usuallly that's the one used)

nodes/cpus_own   .. computers reserved for own use, highest priority
nodes/cpus_shared.. computers shared with others, medium priority
nodes/cpus_other .. computers mainly used by others, lowest priority

nodes/cpus_all   .. all computers in descending priority
nice_dic         .. nice value for each host
"""
import Biskit.mathUtils as MU
from tools import *
import os
import string

import ConfigParser

## Read configuration file
conf = ConfigParser.ConfigParser()
conf.read( os.path.expanduser('~/.biskit/hosts.dat' ) )

def getHosts( section, option ):
    """
    Strip comment from setting (separator #)
    Split list (separator ,)
    
    """
    if conf.has_option( section, option ):
    	setting = conf.get( section, option )
    	setting = string.split( setting, '#' )
	hosts = []
	for host in string.split( setting[0], ' ' ):
    		hosts += [ string.strip( host ) ] 
	return hosts
    return None


def getDict( section, option ):
    """
    Strip comment from setting (separator: #)
    Split entry (separator ,)
    Split dictionary (separator :)
    """
    dic = {}
    for i in getHosts( section, option ):
        i = string.split( i, ':' )
	dic[ string.strip(i[0]) ] =  float( i[1] )
    return dic
   
    
dual = []

##
## computers reserved for own use, highest priority
##
__dual   = getHosts( 'own_hosts', 'dual' ) # dual cpu computers reserved for own use, highest priority, separate with a blank space
__single = getHosts( 'own_hosts', 'single' ) # single cpu computers reserved for own use, highest priority, separate with space

nodes_own = __dual +__single
cpus_own = __dual * 2 + __single

dual += __dual


##
## computers shared with others, medium priority
##
__dual   = getHosts( 'shared_hosts', 'dual' ) # dual cpu computers shared with others, medium priority, separate with space
__single = getHosts( 'shared_hosts', 'single' ) # single cpu computers shared with others, medium priority, separate with space

dual += __dual

nodes_shared = __dual + __single
cpus_shared = __dual * 2 + __single


##
## computers mainly used by others, lowest priority
##
__dual   = getHosts( 'others_hosts', 'dual' ) # dual cpu computers mainly used by others, lowest priority, separate with space
__single = getHosts( 'others_hosts', 'single' ) # dual cpu computers mainly used by others, lowest priority, separate with space

dual += __dual

nodes_other =  __dual + __single
cpus_other = __dual * 2 + __single


##
## all computers in descending priority
##
nodes_all = nodes_own + nodes_shared + nodes_other
cpus_all = cpus_own + cpus_shared + cpus_other


##
## nice value for each computer
##
nice_dic = { 'default':19 }
for h in nodes_own:
    nice_dic[h] = 0
for h in cpus_shared:
    nice_dic[h] = 0
for h in cpus_other:
    nice_dic[h] = 17

##
## temporary fine tuning
##
own_nice = getDict( 'own_hosts', 'nice' ) # nice value for own computers, default 0, separate hosts with space and values with ":" 
shared_nice = getDict( 'shared_hosts', 'nice' ) # nice value shared computers, default 0, separate hosts with space and values with ":"
others_nice = getDict( 'others_hosts', 'nice' ) # nice value others computers, default 17, separate hosts with space and values with ":"

nice_dic.update(own_nice)
nice_dic.update(shared_nice)
nice_dic.update(others_nice)


##
## Memory available for each node, in GB
##
own_ram = getDict( 'own_hosts', 'ram' ) # installed RAM for own computers in GB, default values 0.5 for single and 1 for dual, separate hosts with space and values with ":" 
shared_ram = getDict( 'shared_hosts', 'ram' ) # installed RAM shared computers in GB, default values 0.5 for single and 1 for dual, separate hosts with space and values with ":" 
others_ram = getDict( 'others_hosts', 'ram' ) # installed RAM others computers in GB, default values 0.5 for single and 1 for dual, separate hosts with space and values with ":" 

ram_dic = { 'default':0.5 }
for h in dual:
    ram_dic[ h ] = 1.0
    
ram_dic.update(own_ram)
rma_dic.update(shared_ram)
ram_dic.update(others_ram)


##
## exclude list, temporarily remove some nodes
##
nodes_exclude = []
nodes_exclude += getHosts( 'own_hosts', 'exclude' ) # exclude list, temporarily remove some nodes, separate with space
nodes_exclude += getHosts( 'shared_hosts', 'exclude' ) # exclude list, temporarily remove some nodes, separate with space
nodes_exclude += getHosts( 'others_hosts', 'exclude' ) # exclude list, temporarily remove some nodes, separate with space
 

##
## Switch On / Off temporary removal of nodes 
## Set to 1 to take some nodes out of the list (i.e. don't use them)
## Set to 0 to let Biskit use all computers
##
REMOVE = 1

if REMOVE:
    for l in [ nodes_all, cpus_all, nodes_own, cpus_own,
               nodes_shared, cpus_shared, nodes_other, cpus_other ]:

        MU.removeFromList( l, nodes_exclude, all=1 )
