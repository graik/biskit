import biskit as B
import biskit.tools as T

rec = B.PDBModel('rec/1A2P_clean.pdb')
rec.saveAs('rec/1A2P.model')
recpsf = B.PCRModel(source=rec, fPsf='1A2P.psf', pdbCode='1A2P')

d = {1: recpsf}
T.dump(d, 'rec/1A2P_model.dic')


lig = B.PDBModel('lig/1A19_clean.pdb')
lig.saveAs('lig/1A19.model')
ligpsf = B.PCRModel(source=lig, fPsf='1A19.psf', pdbCode='1A19')

d = {1: ligpsf}
T.dump(d, 'lig/1A19_model.dic')


