##
## Adapt this file and save it under ~/.biskit/hosts.py !
##
"""
The list of cluster computers available to PVM.
Each host must be accessible via ssh w/o password.

This is **NOT** the place to add your custom host list! Biskit.PVM.hosts
only serves as a mirror to read in the local host list configuration
from ~/.biskit/hosts.py or (if this doesn't exist) from
biskit/Biskit/data/defaults/hosts.py. It publishes the following fields
to the rest of Biskit:

-  nodes_* .. lists with one entry per computer
-  cpus_*  .. lists with one entry per CPU (usually that's the one used)

-  nodes/cpus_own   .. computers reserved for own use, highest priority
-  nodes/cpus_shared.. computers shared with others, medium priority
-  nodes/cpus_other .. computers mainly used by others, lowest priority
-  nodes/cpus_all   .. all the above computers in descending priority

-  nice_dic         .. nice value for each host
-  ram_dic          .. memory of each node (if set)

"""
import Biskit.tools as T
import Biskit as B
import os, user

__CFG_DEFAULT = T.dataRoot() + '/defaults/hosts.py'
__CFG_USER    = user.home + '/.biskit/hosts.py'

class HostsError( B.BiskitError ):
    """raised when there is a problem with the hosts.py configuration"""
    pass

##
#: computers reserved for own use, highest priority
nodes_own = []
cpus_own = []

##
#: computers shared with others, medium priority
nodes_shared = []
cpus_shared = []

##
#: computers mainly used by others, lowest priority
nodes_other =  []
cpus_other = []

##
#: all computers in descending priority
nodes_all = nodes_own + nodes_shared + nodes_other
cpus_all = cpus_own + cpus_shared + cpus_other

##
#: nice value for each computer
nice_dic = { 'default':19 }
for h in nodes_own:
    nice_dic[h] = 0
for h in nodes_shared:
    nice_dic[h] = 0
for h in nodes_other:
    nice_dic[h] = 17

##
#: Memory available for each node, in GB
ram_dic = { 'default':0.5 }

###########################
## Read local configuration

try:
    f = __CFG_DEFAULT

    if os.path.exists( __CFG_USER ):
        f = __CFG_USER

    execfile( f )

except Exception, why:
    raise HostsError( 'Cannot read host configuration from %s (%r)' % (f,why) )


################
## empty test ##

import Biskit.test as BT

class Test(BT.BiskitTest):
    """Mock test, hosts only contains data"""
    pass

## some cleanup
del __CFG_DEFAULT, __CFG_USER, T, B, os
