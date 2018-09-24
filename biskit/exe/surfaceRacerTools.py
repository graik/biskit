## numpy-oldnumeric calls replaced by custom script; 09/06/2016
## Automatically adapted for numpy-oldnumeric Mar 26, 2007 by alter_code1.py

##
## Biskit, a toolkit for the manipulation of macromolecular structures
## Copyright (C) 2004-2018 Raik Gruenberg & Johan Leckner
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
"""
Reference data for SurfaceRacer.
"""

import biskit.core.oldnumeric as N0

## dictionary with average accessabilities in 500 random
## GLY-XX-GLY peptides calculated with SurfaceRacer using
## a probe radius of 1.4 A and the Richards vdw radii set

## random coil molecular surface in a GLY-XXX-GLY peptide
ranMS = {'CYS': {'C'  :  1.969, 'CB' : 21.339, 'CA' :  8.077,
                 'O'  : 14.038, 'N'  :  5.713, 'SG' : 34.959},
         'GLN': {'C'  :  1.736, 'CB' : 18.742, 'CA' :  6.384,
                 'CG' : 20.672, 'O'  : 14.136, 'CD' :  2.825,
                 'N'  :  7.409, 'NE2': 23.326, 'OE1': 14.894},
         'ILE': {'C'  :  1.151, 'CB' :  7.933, 'CA' :  4.369,
                 'O'  : 14.185,  'N' :  5.428, 'CD1': 31.332,
                 'CG1': 17.690, 'CG2': 28.741},
         'SER': {'C'  :  2.336, 'OG' : 17.413, 'CB' : 26.667,
                 'CA' :  8.803, 'O'  : 13.894, 'N'  :  6.849},
         'VAL': {'C'  :  1.332, 'CB' :  7.877, 'CA' :  6.509,
                 'O'  : 13.741, 'N'  :  4.786, 'CG1': 31.049,
                 'CG2': 30.663},
         'LYS': {'C'  :  1.722, 'CB' : 18.077, 'CA' :  6.360,
                 'CG' : 16.832, 'CE' : 18.940, 'CD' : 16.329,
                 'NZ' : 32.635, 'O'  : 14.175, 'N'  :  7.007},
         'ASN': {'C'  :  2.050, 'CB' : 13.362, 'CA' : 15.158,
                 'CG' : 20.928, 'O'  :  6.352, 'N'  :  1.645,
                 'OD1': 14.146, 'ND2':  6.379 },
         'PRO': {'C'  : 15.067, 'CB' : 14.374, 'CA' :  5.090,
                 'CG' :  4.201, 'O'  : 15.123, 'CD' : 15.837,
                 'N'  :  1.779},
         'GLY': {'CA' : 20.259, 'C'  :  5.731, 'O'  : 21.192,
                 'N'  : 15.386},
         'THR': {'C'  : 19.678, 'CB' :  0.712, 'CA' :  1.439,
                 'OG1': 11.571, 'O'  :  8.170, 'N'  : 15.934,
                 'CG2': 14.310},
         'PHE': {'C'  :  6.163, 'CE1': 32.630, 'CB' :  1.976,
                 'CA' : 12.422, 'CG' : 21.216, 'O'  :  7.309,
                 'N'  :  1.825, 'CZ' : 14.120, 'CD1':  6.826,
                 'CD2': 15.219, 'CE2': 12.588},
         'ALA': {'CB' : 15.389, 'CA' : 15.277, 'C'  : 34.475,
                 'O'  :  9.672, 'N'  :  1.697},
         'HIS': {'C'  : 14.138, 'CE1':  7.174, 'CB' : 26.080,
                 'CA' :  3.136, 'CG' : 13.978, 'O'  :  8.065,
                 'N'  :  1.666, 'CD2': 14.372, 'ND1': 21.724,
                 'NE2':  8.415},
         'MET': {'C'  :  1.564, 'CB' : 14.093, 'CA' :  7.444,
                 'CG' : 15.702, 'CE' : 10.286, 'N'  : 15.480,
                 'O'  :  2.250, 'SD' : 19.502},
         'ASP': {'C'  :  7.110, 'CB' : 22.634, 'CA' : 14.257,
                 'CG' :  3.331, 'O'  : 14.934, 'N'  :  5.535,
                 'OD1': 14.057, 'OD2':  1.949},
         'LEU': {'C'  : 16.627, 'CB' :  5.280, 'CA' :  8.389,
                 'CG' : 14.299, 'O'  :  5.688, 'N'  : 30.423,
                 'CD1': 29.766, 'CD2':  1.777},
         'ARG': {'C'  : 18.112, 'CB' :  5.241, 'CA' : 16.730,
                 'CG' : 16.616, 'NE' : 13.828, 'O'  : 15.959,
                 'CD' :  1.912, 'CZ' : 21.727, 'NH1': 21.265,
                 'NH2':  7.178, 'N'  :  1.769},
         'TRP': {'C'  : 22.097, 'CZ2':  9.303, 'CB' :  4.131,
                 'CA' : 14.342, 'CG' :  5.499, 'O'  : 16.075,
                 'N'  : 15.575, 'CH2':  1.621, 'CE3': 22.125,
                 'CE2':  8.268, 'CD2':  2.468, 'CZ3': 13.935,
                 'NE1':  7.413, 'CD1': 15.173},
         'GLU': {'C'  : 22.925, 'CB' :  2.084, 'CA' : 15.389,
                 'CG' : 20.656, 'O'  : 21.058, 'CD' :  7.101,
                 'OE2':  2.090, 'N'  : 14.298, 'OE1':  6.603},
         'TYR': {'C'  :  5.766, 'CD2': 12.756, 'OH' : 13.124,
                 'CB' : 15.263, 'CA' :  1.525, 'CG' : 18.792,
                 'O'  :  6.204, 'N'  : 18.807, 'CZ' : 35.593,
                 'CD1':  7.103, 'CE1': 14.230, 'CE2': 14.968}}

## random coil solvent accessible surface in a GLY-XXX-GLY peptide
ranAS = {'CYS': {'C'  :  0.833, 'CB' : 34.358, 'CA' : 7.768,
                 'O'  : 29.027, 'N'  :  5.159, 'SG' : 76.311},
         'GLN': {'C'  :  0.694, 'CB' : 28.846, 'CA' :  3.896,
                 'CG' : 32.484, 'O'  : 28.716, 'CD' :  2.467,
                 'N'  :  7.736, 'NE2': 54.470, 'OE1': 32.637},
         'ILE': {'C'  :  0.336, 'CB' : 10.203, 'CA' :  2.117,
                 'O'  : 26.670, 'N'  :  4.723, 'CD1': 63.709,
                 'CG1': 32.543, 'CG2': 50.456},
         'SER': {'C'  :  1.003, 'OG' : 47.009, 'CB' : 44.666,
                 'CA' :  8.688, 'O'  : 29.511, 'N'  :  6.051},
         'VAL': {'C'  :  0.491, 'CB' :  9.140, 'CA' :  5.240,
                 'O'  : 23.612, 'N'  :  2.422, 'CG1': 65.475,
                 'CG2': 59.022},
         'LYS': {'C'  :  0.713, 'CB' : 27.111, 'CA' :  4.234,
                 'CG' : 24.319, 'CE' : 35.072, 'CD' : 23.182,
                 'NZ' : 76.704, 'O'  : 27.814, 'N'  :  7.147},
         'ASN': {'C'  :  0.671, 'CB' : 33.549, 'CA' :  6.662,
                 'CG' :  2.403, 'O'  : 27.881, 'N'  :  7.635,
                 'OD1': 36.200, 'ND2': 50.093},
         'PRO': {'C'  :  0.696, 'CB' : 34.966, 'CA' :  4.790,
                 'CG' : 46.846, 'O'  : 34.264, 'CD' : 32.028,
                 'N'  :  0.314},
         'GLY': {'CA' : 43.741, 'C'  :  2.298, 'O'  : 31.107,
                 'N'  : 12.702},
         'THR': {'C'  :  0.709, 'CB' : 15.180, 'CA' : 7.358,
                 'OG1': 38.230, 'O'  : 27.475, 'N'  : 4.937,
                 'CG2': 62.021},
         'PHE': {'C'  :  0.876, 'CE1': 33.974, 'CB' : 30.431,
                 'CA' :  4.790, 'CG' :  1.388, 'O'  : 28.205,
                 'N'  :  6.352, 'CZ' : 34.444, 'CD1': 19.752,
                 'CD2': 18.580, 'CE2': 32.812},
         'ALA': {'CB' : 70.070, 'CA' : 11.256, 'C'  :  0.717,
                 'O'  : 27.891, 'N'  :  7.782},
         'HIS': {'C'  :  0.923, 'CE1': 32.391, 'CB' : 31.658,
                 'CA' :  6.665, 'CG' :  1.280, 'O'  : 27.750,
                 'N'  :  7.365, 'CD2': 27.492, 'ND1': 15.030,
                 'NE2': 36.247},
         'MET': {'C'  :  0.590, 'CB' : 28.322, 'CA' :  3.961,
                 'CG' : 29.567, 'CE' : 74.891, 'N'  :  7.345,
                 'O'  : 27.658, 'SD' : 29.550},
         'ASP': {'C'  :  0.972, 'CB' : 37.245, 'CA' : 10.246,
                 'CG' :  3.665, 'O'  : 25.235, 'N'  :  3.251,
                 'OD1': 37.525, 'OD2': 33.852},
         'LEU': {'C'  :  0.793, 'CB' : 21.950, 'CA' :  2.852,
                 'CG' : 14.129, 'O'  : 27.978, 'N'  :  5.231,
                 'CD1': 61.814, 'CD2': 59.042},
         'ARG': {'C'  :  0.829, 'CB' : 28.313, 'CA' :  2.930,
                 'CG' : 27.099, 'NE' : 23.452, 'O'  : 28.419,
                 'CD' : 23.936, 'CZ' :  1.719, 'NH1': 50.063,
                 'NH2': 46.094, 'N'  :  7.109},
         'TRP': {'C'  :  0.856, 'CZ2': 31.924, 'CB' : 28.556,
                 'CA' :  3.339, 'CG' :  1.337, 'O'  : 28.269,
                 'N'  :  5.626, 'CH2': 34.844, 'CE3': 20.285,
                 'CE2':  4.521, 'CD2':  3.335, 'CZ3': 34.289,
                 'NE1': 32.834, 'CD1': 23.722},
         'GLU': {'C'  :  0.987, 'CB' : 27.019, 'CA' :  5.321,
                 'CG' : 41.443, 'O'  : 29.082, 'CD' :  3.502,
                 'OE2': 35.857, 'N'  :  4.529, 'OE1': 32.074},
         'TYR': {'C'  :  0.777, 'CD2': 20.732, 'OH' : 56.712,
                 'CB' : 29.172, 'CA' :  4.380, 'CG' :  1.517,
                 'O'  : 28.824, 'N'  :  5.930, 'CZ' :  6.352,
                 'CD1': 19.637, 'CE1': 30.502, 'CE2': 30.358}}

## C-terminal random coil molecular surface in a GLY-GLY-XXX peptide
ranMS_C = {'CYS': {'C'  :  4.901, 'OXT': 16.009, 'CB' : 21.505, 'CA' :  8.895,
                   'O'  : 15.937, 'N'  :  5.493, 'SG' : 32.718},
           'GLN': {'C'  :  4.268, 'OXT': 16.327, 'CB' : 18.696, 'CA' :  6.661,
                   'CG' : 21.508, 'O'  : 15.993, 'N'  :  6.289, 'CD' :  2.908,
                   'NE2': 23.350, 'OE1': 15.560},
           'ILE': {'C'  :  3.954, 'OXT': 15.895, 'CB' :  7.210, 'CA' :  6.836,
                   'O'  : 15.297, 'N'  :  5.218, 'CD1': 32.674, 'CG1': 17.185,
                   'CG2': 29.494  },
           'SER': {'C'  :  4.935, 'OG' : 17.498, 'CB' : 26.752, 'CA' :  9.621,
                   'O'  : 15.997, 'N'  :  6.296, 'OXT': 16.009  },
           'VAL': {'C'  :  3.886, 'OXT': 16.344, 'CB' :  8.237, 'CA' :  6.639,
                   'O'  : 16.111, 'N'  :  5.403, 'CG1': 31.003, 'CG2': 30.984},
           'LYS': {'C'  :  4.543, 'OXT': 15.863, 'CB' : 17.773, 'CA' :  9.434,
                   'CG' : 16.113, 'O'  : 15.696, 'N'  :  5.357, 'NZ' : 32.439,
                   'CE' : 18.742, 'CD' : 17.682},
           'TRP': {'C'  :  4.480, 'CD1': 13.511, 'CZ2': 15.023, 'CB' : 21.117,
                   'CA' :  7.578, 'CG' :  1.631, 'O'  : 16.137, 'N'  :  6.173,
                   'OXT': 16.517, 'CH2': 15.027, 'CE3': 13.489, 'CE2':  4.452,
                   'CD2':  3.832, 'CZ3': 14.883, 'NE1': 15.688},
           'PRO': {'C'  :  4.198, 'OXT': 16.105, 'CB' : 20.472, 'CA' : 10.419,
                   'CG' : 21.513, 'O'  : 16.161, 'CD' : 19.308, 'N'  :  0.553},
           'THR': {'C'  :  4.033, 'OXT': 16.336, 'CB' : 11.889, 'CA' :  7.479,
                   'OG1': 15.555, 'O'  : 16.001, 'N'  :  5.930, 'CG2': 33.396},
           'PHE': {'C'  :  4.941, 'CD2': 10.800, 'OXT': 16.072, 'CB' : 21.817,
                   'CA' :  7.746, 'CG' :  1.710, 'O'  : 15.997, 'N'  :  5.831,
                   'CZ' : 15.263, 'CD1': 12.185, 'CE1': 15.259, 'CE2': 14.726},
           'ALA': {'C'  :  4.785, 'OXT': 15.918, 'CB' : 34.751, 'CA' : 10.288,
                   'O'  : 15.844, 'N'  :  6.609},
           'GLY': {'OXT': 15.901, 'CA' : 26.945, 'C'  :  5.589, 'O'  : 15.783,
                   'N'  :  8.017},
           'HIS': {'C'  :  4.537, 'CD2': 13.867, 'OXT': 15.778, 'CB' : 22.103,
                   'CA' :  9.755, 'CG' :  1.686, 'O'  : 15.257, 'N'  :  5.629,
                   'CE1': 15.805, 'ND1': 12.331, 'NE2': 15.544  },
           'GLU': {'C'  :  4.787, 'OXT': 16.161, 'CB' : 20.204, 'CA' :  7.939,
                   'CG' : 20.793, 'O'  : 16.023, 'N'  :  5.478, 'OE2': 15.843,
                   'CD' :  4.832, 'OE1': 15.462},
           'LEU': {'C'  :  4.225, 'OXT': 16.433, 'CB' : 15.593, 'CA' :  6.733,
                   'CG' :  8.037, 'O'  : 16.122, 'N'  :  5.857, 'CD1': 31.131,
                   'CD2': 31.326},
           'ARG': {'C'  :  4.178, 'OXT': 16.028, 'CB' : 18.280, 'CA' :  9.543,
                   'CG' : 16.758, 'NE' : 20.127, 'O'  : 14.080, 'N'  :  4.504,
                   'CZ' :  2.192, 'NH1': 20.770, 'NH2': 22.158, 'CD' : 13.552},
           'ASP': {'C'  :  4.437, 'OXT': 16.092, 'CB' : 22.783, 'CA' :  9.494,
                   'CG' :  4.258, 'O'  : 16.423, 'N'  :  5.150, 'OD1': 16.141,
                   'OD2': 16.226},
           'ASN': {'C'  :  4.200, 'OXT': 16.330, 'CB' : 22.238, 'CA' :  8.738,
                   'CG' :  2.759, 'O'  : 16.268, 'N'  :  6.718, 'OD1': 15.585,
                   'ND2': 23.497},
           'TYR': {'C'  :  4.257, 'CE1': 15.574, 'OH' : 20.805, 'OXT': 16.007,
                   'CB' : 21.243, 'CA' :  9.453, 'CG' :  1.708, 'O'  : 16.290,
                   'N'  :  4.981, 'CZ' :  5.799, 'CD1': 10.416, 'CD2': 12.459,
                   'CE2': 15.491},
           'MET': {'C'  :  4.692, 'OXT': 16.239, 'CB' : 18.674, 'CA' :  7.639,
                   'CG' : 19.206, 'O'  : 15.957, 'N'  :  6.066, 'CE' : 36.053,
                   'SD' : 15.153  }}



## N-terminal random coil molecular surface in a XXX-GLY-GLY peptide
ranMS_N = {'CYS': {'C'  :  1.881, 'CB' : 20.689, 'CA' : 12.751, 'O'  : 12.456,
                   'N'  : 17.890, 'SG' : 33.579}, 
           'GLN': {'C'  :  1.838, 'CB' : 18.800, 'CA' : 12.632, 'CG' : 18.923,
                   'O'  : 12.560, 'CD' :  3.315, 'N'  : 17.352, 'NE2': 23.831,
                   'OE1': 13.464}, 
           'ILE': {'C'  :  1.681, 'CB' :  7.893, 'CA' :  8.675, 'O'  : 13.221,
                   'N'  : 18.819, 'CD1': 31.268, 'CG1': 17.636, 'CG2': 27.787},
           'SER': {'C'  :  2.303, 'OG' : 16.770, 'CB' : 26.889, 'CA' : 13.107,
                   'O'  : 13.152, 'N'  : 19.013}, 
           'VAL': {'C'  :  1.649, 'CB' :  8.327, 'CA' : 10.304, 'O'  : 13.305,
                   'N'  : 18.515, 'CG1': 29.069, 'CG2': 30.586}, 
           'LYS': {'C'  :  2.076, 'CB' : 16.542, 'CA' : 11.608, 'CG' : 17.298,
                   'CE' : 18.221, 'CD' : 17.673, 'NZ' : 32.602, 'O'  : 13.474,
                   'N'  : 19.247}, 
           'TRP': {'C'  :  1.788, 'CD1':  9.793, 'CZ2': 14.896, 'CB' : 19.663,
                   'CA' : 12.625, 'CG' :  1.542, 'O'  : 12.968, 'N'  : 18.599,
                   'CH2': 14.820, 'CE3': 13.676, 'CE2':  4.877, 'CD2':  4.279,
                   'CZ3': 14.890, 'NE1': 15.654},
           'PRO': {'C'  :  2.460, 'CB' : 21.196, 'CA' : 12.083, 'CG' : 21.410,
                   'O'  : 13.475, 'CD' : 24.855, 'N'  :  8.681 }, 
           'THR': {'C'  :  2.049, 'CB' : 11.750, 'CA' : 11.666, 'OG1': 15.173,
                   'O'  : 13.350, 'N'  : 18.585, 'CG2': 32.218},
           'PHE': {'C'  :  1.806, 'CD2': 12.542, 'CB' : 20.567, 'CA' : 12.435,
                   'CG' :  1.700, 'O'  : 12.991, 'N'  : 18.878, 'CZ' : 15.122,
                   'CD1':  9.868, 'CE1': 15.353, 'CE2': 15.245}, 
           'ALA': {'CB' : 34.147, 'CA' : 13.086, 'C'  :  2.120, 'O'  : 13.312,
                   'N'  : 19.125}, 
           'GLY': {'CA' : 29.487, 'C'  :  3.092, 'O'  : 13.362, 'N'  : 19.552},
           'HIS': {'C'  :  1.800, 'CD2': 13.812, 'CB' : 21.222, 'CA' : 12.204,
                   'CG' :  1.448, 'O'  : 12.229, 'N'  : 19.066, 'CE1': 15.538,
                   'ND1': 10.390, 'NE2': 15.168},
           'GLU': {'C'  :  2.148, 'CB' : 19.736, 'CA' : 11.607, 'CG' : 21.964,
                   'O'  : 13.334, 'CD' :  4.935, 'OE2': 15.359, 'N'  : 18.703,
                   'OE1': 15.546}, 
           'LEU': {'C'  :  1.846, 'CB' : 15.176, 'CA' : 12.338, 'CG' :  7.083,
                   'O'  : 12.648, 'N'  : 18.191, 'CD1': 29.800, 'CD2': 30.607},
           'ARG': {'C'  :  1.822, 'CB' : 16.419, 'CA' : 12.537, 'CG' : 16.382,
                   'NE' : 18.534, 'O'  : 12.524, 'CD' : 17.014, 'CZ' :  2.161,
                   'NH1': 22.420, 'NH2': 21.693, 'N'  : 18.034}, 
           'ASP': {'C'  :  1.943, 'CB' : 22.100, 'CA' : 13.127, 'CG' :  3.754,
                   'O'  : 13.528, 'N'  : 18.370, 'OD1': 14.891, 'OD2': 15.499},
           'ASN': {'C'  :  1.781, 'CB' : 21.494, 'CA' : 12.711, 'CG' :  2.215,
                   'O'  : 12.685, 'N'  : 18.523, 'OD1': 13.989, 'ND2': 21.674},
           'TYR': {'C'  :  1.891, 'CE1': 15.296, 'OH' : 20.825, 'CB' : 20.631,
                   'CA' : 12.481, 'CG' :  1.666, 'O'  : 12.892, 'N'  : 18.798,
                   'CZ' :  5.542, 'CD1': 12.374, 'CD2': 10.473, 'CE2': 14.935},
           'MET': {'C'  :  1.970, 'CB' : 17.536, 'CA' : 11.852, 'CG' : 18.881,
                   'CE' : 36.194, 'N'  : 19.261, 'O'  : 13.279, 'SD' : 15.618}}

## C-terminal random coil solvent accessible surface in a GLY-GLY-XXX peptide
ranAS_C = {'CYS': {'C'  :  3.576, 'OXT': 44.000, 'CB' : 35.823, 'CA' :  7.883,
                   'O'  : 40.390, 'N'  :  3.953, 'SG' : 63.386},
           'GLN': {'C'  :  2.560, 'OXT': 39.663, 'CB' : 23.902, 'CA' :  3.636,
                   'CG' : 35.848, 'O'  : 39.558, 'N'  :  5.031, 'CD' :  2.604,
                   'NE2': 52.310, 'OE1': 37.289},
           'ILE': {'C'  :  2.192, 'OXT': 38.990, 'CB' :  6.897, 'CA' :  4.449,
                   'O'  : 35.608, 'N'  :  3.047, 'CD1': 73.231, 'CG1': 27.316,
                   'CG2': 53.373},
           'SER': {'C'  :  3.562, 'OG' : 44.874, 'CB' : 45.772, 'CA' :  9.377,
                   'O'  : 40.102, 'N'  :  5.281, 'OXT': 44.024},
           'VAL': {'C'  :  1.947, 'OXT': 40.714, 'CB' : 10.749, 'CA' :  4.843,
                   'O'  : 38.524, 'N'  :  4.350, 'CG1': 64.304, 'CG2': 61.298},
           'LYS': {'C'  :  3.209, 'OXT': 43.059, 'CB' : 25.382, 'CA' :  9.294,
                   'CG' : 17.124, 'O'  : 37.874, 'N'  :  3.358, 'NZ' : 78.334,
                   'CE' : 34.282, 'CD' : 28.882},
           'TRP': {'C'  :  2.864, 'CD1': 23.471, 'CZ2': 30.132, 'CB' : 27.372,
                   'CA' :  4.926, 'CG' :  1.362, 'O'  : 39.445, 'N'  :  4.670,
                   'OXT': 42.192, 'CH2': 33.477, 'CE3': 21.297, 'CE2':  3.987,
                   'CD2':  3.165, 'CZ3': 32.837, 'NE1': 31.635},
           'PRO': {'C'  :  2.522, 'OXT': 42.757, 'CB' : 39.165, 'CA' : 12.754,
                   'CG' : 45.959, 'O'  : 36.644, 'CD' : 30.229, 'N'  :  0.162},
           'THR': {'C'  :  2.135, 'OXT': 41.079, 'CB' : 16.180, 'CA' : 5.510,
                   'OG1': 31.767, 'O'  : 39.334, 'N'  :  4.732, 'CG2': 69.879},
           'PHE': {'C'  :  3.491, 'CD2': 16.338, 'OXT': 43.822, 'CB' : 31.885,
                   'CA' :  5.300, 'CG' :  1.482, 'O'  : 40.849, 'N'  :  4.425,
                   'CZ' : 32.891, 'CD1': 20.269, 'CE1': 32.231, 'CE2': 27.941},
           'ALA': {'C'  :  3.470, 'OXT': 44.166, 'CB' : 69.869, 'CA' : 10.822,
                   'O'  : 40.456, 'N'  :  5.644 },
           'GLY': {'OXT': 44.099, 'CA' : 44.499, 'C'  :  5.688, 'O'  : 44.236,
                   'N'  : 10.787},
           'HIS': {'C'  :  3.235, 'CD2': 26.008, 'OXT': 42.132, 'CB' : 35.445,
                   'CA' :  9.836, 'CG' :  1.313, 'O'  : 34.817, 'N'  :  3.437,
                   'CE1': 36.648, 'ND1': 16.111, 'NE2': 36.711},
           'GLU': {'C'  :  3.233, 'OXT': 43.638, 'CB' : 25.879, 'CA' :  6.755,
                   'CG' : 28.174, 'O'  : 40.299, 'N'  :  4.009, 'OE2': 38.602,
                   'CD' :  4.560, 'OE1': 43.505},
           'LEU': {'C'  :  2.492, 'OXT': 41.188, 'CB' : 13.565, 'CA' :  3.909,
                   'CG' :  9.503, 'O'  : 38.884, 'N'  :  4.293,
                   'CD1': 62.247, 'CD2': 70.859},
           'ARG': {'C'  :  2.842, 'OXT': 43.324, 'CB' : 30.232, 'CA' :  9.378,
                   'CG' : 26.275, 'NE' : 31.446, 'O'  : 28.475, 'N'  :  2.424,
                   'CZ' :  2.091, 'NH1': 39.892, 'NH2': 58.200, 'CD' : 16.212},
           'ASP': {'C'  :  2.712, 'OXT': 42.871, 'CB' : 37.628, 'CA' :  9.681,
                   'CG' :  3.708, 'O'  : 35.649, 'N'  :  2.738, 'OD1': 35.467,
                   'OD2': 41.483},
           'ASN': {'C'  :  2.654, 'OXT': 40.748, 'CB' : 30.645, 'CA' :  6.482,
                   'CG' :  2.489, 'O'  : 38.499, 'N'  :  5.569, 'OD1': 38.650,
                   'ND2': 54.095},
           'TYR': {'C'  :  3.012, 'CE1': 25.434, 'OH' : 56.665, 'OXT': 43.093,
                   'CB' : 33.823, 'CA' :  9.198, 'CG' :  1.325, 'O'  : 35.242,
                   'N'  :  2.741, 'CZ' :  5.893, 'CD1': 11.890, 'CD2': 18.390,
                   'CE2': 30.470},
           'MET': {'C'  :  3.048, 'OXT': 43.014, 'CB' : 23.902, 'CA' :  5.848,
                   'CG' : 28.573, 'O'  : 39.775, 'N'  :  4.874, 'CE' : 76.639,
                   'SD' : 32.812}}



## N-terminal random coil solvent accessible surface in a XXX-GLY-GLY peptide
ranAS_N = {'CYS': {'C'  :  1.431, 'CB' : 38.159, 'CA' : 19.289, 'O'  : 23.438,
                   'N'  : 43.276, 'SG' : 67.559},
           'GLN': {'C'  :  1.649, 'CB' : 31.979, 'CA' : 19.263, 'CG' : 23.964,
                   'O'  : 21.682, 'CD' :  2.755, 'N'  : 38.957, 'NE2': 56.478,
                   'OE1': 31.267},
           'ILE': {'C'  :  1.024, 'CB' :  9.940, 'CA' :  8.458, 'O'  : 24.997,
                   'N'  : 42.475, 'CD1': 68.526, 'CG1': 33.941, 'CG2': 46.229},
           'SER': {'C'  :  1.629, 'OG' : 48.780, 'CB' : 42.463, 'CA' : 18.962,
                   'O'  : 26.451, 'N'  : 43.920},
           'VAL': {'C'  :  0.936, 'CB' : 10.647, 'CA' : 12.837, 'O'  : 25.260,
                   'N'  : 42.960, 'CG1': 52.342, 'CG2': 66.747},
           'LYS': {'C'  :  1.326, 'CB' : 17.475, 'CA' : 15.471, 'CG' : 28.459,
                   'CE' : 32.971, 'CD' : 32.126, 'NZ' : 75.829, 'O'  : 27.414,
                   'N'  : 45.390},
           'TRP': {'C'  :  1.715, 'CD1': 10.467, 'CZ2': 31.117, 'CB' : 32.211,
                   'CA' : 19.563, 'CG' :  1.157, 'O'  : 23.389, 'N'  : 38.572,
                   'CH2': 34.299, 'CE3': 21.185, 'CE2':  4.574, 'CD2':  3.138,
                   'CZ3': 34.416, 'NE1': 27.704},
           'PRO': {'C'  :  1.586, 'CB' : 36.124, 'CA' : 17.336, 'CG' : 47.576,
                   'O'  : 27.825, 'CD' : 51.994, 'N'  : 19.511},
           'THR': {'C'  :  1.375, 'CB' : 15.023, 'CA' : 15.341, 'OG1': 40.161,
                   'O'  : 26.875, 'N'  : 42.884, 'CG2': 61.218},
           'PHE': {'C'  :  1.579, 'CD2': 19.828, 'CB' : 34.960, 'CA' : 18.791,
                   'CG' :  1.281, 'O'  : 22.256, 'N'  : 42.396, 'CZ' : 34.159,
                   'CD1': 10.942, 'CE1': 28.856, 'CE2': 34.476},
           'ALA': {'CB' : 70.762, 'CA' : 19.700, 'C'  :  1.248, 'O'  : 27.698,
                   'N'  : 48.593},
           'GLY': {'CA' : 54.054, 'C'  :  2.369, 'O'  : 30.073, 'N'  : 53.879},
           'HIS': {'C'  :  1.443, 'CD2': 26.112, 'CB' : 35.612, 'CA' : 17.834,
                   'CG' :  1.116, 'O'  : 21.813, 'N'  : 45.116, 'CE1': 32.814,
                   'ND1': 12.827, 'NE2': 35.149},
           'GLU': {'C'  :  1.532, 'CB' : 21.289, 'CA' : 15.421, 'CG' : 34.752,
                   'O'  : 27.491, 'CD' :  4.877, 'OE2': 43.069, 'N'  : 43.177,
                   'OE1': 44.404},
           'LEU': {'C'  :  1.453, 'CB' : 22.380, 'CA' : 18.248, 'CG' :  7.839,
                   'O'  : 24.519, 'N'  : 42.514, 'CD1': 57.000, 'CD2': 66.283},
           'ARG': {'C'  :  1.558, 'CB' : 25.193, 'CA' : 18.837, 'CG' : 21.292,
                   'NE' : 23.924, 'O'  : 20.214, 'CD' : 24.711, 'CZ' :  2.057,
                   'NH1': 57.939, 'NH2': 50.567, 'N'  : 39.687},
           'ASP': {'C'  :  1.748, 'CB' : 38.818, 'CA' : 20.482, 'CG' :  3.560,
                   'O'  : 25.296, 'N'  : 42.730, 'OD1': 31.978, 'OD2': 36.750},
           'ASN': {'C'  :  1.576, 'CB' : 36.560, 'CA' : 19.196, 'CG' :  2.215,
                   'O'  : 22.760, 'N'  : 44.258, 'OD1': 29.512, 'ND2': 44.276},
           'TYR': {'C'  :  1.626, 'CE1': 30.412, 'OH' : 55.461, 'CB' : 34.783,
                   'CA' : 18.887, 'CG' :  1.309, 'O'  : 22.818, 'N'  : 42.020,
                   'CZ' :  5.631, 'CD1': 19.492, 'CD2': 12.154, 'CE2': 24.738},
           'MET': {'C'  :  1.458, 'CB' : 18.265, 'CA' : 16.140, 'CG' : 28.567,
                   'CE' : 79.176, 'N'  : 45.454, 'O'  : 25.118, 'SD' : 36.038}}



def relExposure( model, absSurf, key='AS', clip=1 ):
    """
    Calculate how exposed an atom is relative to the same
    atom in a GLY-XXX-GLY tripeptide, an approximation of
    the unfolded state.

    @param absSurf: Absolute MS OR AS values
    @type  absSurf: [float]
    @param key: MS or AS
    @type  key: MS|AS
    @param clip: clip values above 100% (default: 1)
    @type  clip: 1|0
    
    @return: rel - list of relative accessible surfaces
    @rtype: [float]
    """
    if not key=='MS' and not key=='AS':
        raise Exception('Incorrect key for relative exposiure: %s '%key)

    rel = []
    i=0

    ## loop over chains
    for j in range( model.lenChains()):
        c = model.takeChains([j])

        k=0
        cIdx = c.resIndex()
        ## and loop over atoms in chain
        for a in c.atoms.iterDicts():
            ## N-terminal residue
            if k < cIdx[1]:
                rel = __Nter( a, rel, absSurf, key, i )
            ## C-terminal residue
            if k >= cIdx[-1]:
                rel = __Cter( a, rel, absSurf, key, i )
            ## everything but N- and C termini
            if not k < cIdx[1] and not k >= cIdx[-1]:
                rel = __bulk( a, rel, absSurf, key, i )
            i+=1
            k+=1

    if clip:
        return  N0.clip( N0.array(rel), 0.0, 100.0 )
    else:
        return  N0.array(rel)


def __Nter( a, rel, absSurf, key, i ):
    """
    Get N-terminal relative exposures.

    @param rel: list in which relative exposures are collected
    @type  rel: [float]
    @param absSurf: Absolute MS OR AS values
    @type  absSurf: [float]
    @param key: 'MS' or 'AS'
    @type  key: str
    
    @return: rel - list of relative accessible surfaces with data from
                   N-terminal residues appended.
    @rtype: [float]
    """
    atom = a['name']
    resi = a['residue_name']
    ## don't want to divide with zero
    if absSurf[i] != 0:
        ## if solvent accessible surface
        if key == 'AS':
            rel += [ (absSurf[i] / ranAS_N[resi][atom]) *100 ]
        ## if molecular surface
        if key == 'MS':
            rel += [ (absSurf[i] / ranMS_N[resi][atom]) *100 ]
    else:
        rel += [0.0]
    return rel


def __Cter( a, rel, absSurf, key, i ):
    """
    Get C-terminal relative exposures

    @param rel: list in which relative exposures are collected
    @type  rel: [float]
    @param absSurf: Absolute MS OR AS values
    @type  absSurf: [float]
    @param key: MS or AS
    @type  key: MS|AS
    
    @return: rel - list of relative accessible surfaces with data from
                   C-terminal residues appended.
    @rtype: [float]    
    """
    atom = a['name']
    resi = a['residue_name']
    if absSurf[i] != 0:
        if key == 'AS':
            rel += [ (absSurf[i] / ranAS_C[resi][atom]) *100 ]
        if key == 'MS':
            rel += [ (absSurf[i] / ranMS_C[resi][atom]) *100 ]
    else:
        rel += [0.0]
    return rel


def __bulk( a, rel, absSurf, key, i ):
    """
    Get relative exposures for everything else.

    @param rel: list in which relative exposures are collected
    @type  rel: [float]
    @param absSurf: Absolute MS OR AS values
    @type  absSurf: [float]
    @param key: MS or AS
    @type  key: MS|AS
    
    @return: rel - list of relative accessible surfaces with data from
    '              none-terminal residues appended.
    @rtype: [float]   
    """
    atom = a['name']
    resi = a['residue_name']
    if absSurf[i] != 0:
        if key == 'AS':
            rel += [ (absSurf[i] / ranAS[resi][atom]) *100 ]
        if key == 'MS':
            rel += [ (absSurf[i] / ranMS[resi][atom]) *100 ]
    else:
        rel += [0.0]
    return rel



#############
##  TESTING        
#############
import biskit.test as BT
        
class Test(BT.BiskitTest):
    """Test case"""

    def test_surfaceRacerTools(self):
        """surfaceRacerTools test"""
        from biskit import PDBModel
        import biskit.tools as T
        
        ## load a structure
        self.m = PDBModel( T.testRoot()+'/lig/1A19.pdb' )
        self.m = self.m.compress( self.m.maskProtein() )
        self.m = self.m.compress( self.m.maskHeavy() )
        
        ## some fake surface data
        surf = N0.ones( self.m.lenAtoms()) * 10.0

        relExp = relExposure( self.m, surf )
        
##         if self.local:
##             globals().update( locals() )
            
        self.assertAlmostEqual( N0.sum(relExp), 44276.86085222386, 8 )
    
        
if __name__ == '__main__':

    BT.localTest()
