##
## Adapt this file and save it under ~/.biskit/hosts.py !
##
"""
The list of cluster computers available to PVM.
! Each host must be accessible via ssh w/o password !

Adapt this file (or rather a copy in ~/.biskit/hosts.py) to your needs.

The python code in this file will be executed in the scope of
Biskit.PVM.hosts.py.
"""

##
#: computers reserved for own use, highest priority
#
#  Example::
#     nodes_own = ['mynode1', 'mynode2' ]
#     cpus_own  = ['mynode1', 'mynode1', 'mynode2', 'mynode2' ]
#  That means: nodes_own should have one entry per machine, and cpus_own should
#              have one entry per available CPU.
#  You can shorten the above cpus_own line with normal python magic::
#     cpus_own  = nodes_own * 2
nodes_own = ['localhost']
cpus_own = ['localhost']

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
#
# Should be fairly obvious. You can assign nice values to individual nodes
# like this::
#     nice_dic['mynicenode'] = 17
# Or just use the below code to assign the same nice level to groups of nodes.
nice_dic = { 'default':19 }
for h in nodes_own:
    nice_dic[h] = 0
for h in nodes_shared:
    nice_dic[h] = 0
for h in nodes_other:
    nice_dic[h] = 17

##
#: Memory available for each node, in GB
#
#  Same principle as above. This dictionary is only very rarely used though.
ram_dic = { 'default':0.5 }
