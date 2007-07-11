## Template for modeller input file for multiple template modelling
## prepared to be filled out by Biskit/Mod/Modeller.py

from modeller import *
from modeller.automodel import *

env = environ()

a = automodel(env, alnfile='%(f_pir)s',
              knowns=( %(knowns_py)s ),
              sequence='%(target_id)s')

a.starting_model = %(starting_model)i
a.ending_model = %(ending_model)i

a.make()
