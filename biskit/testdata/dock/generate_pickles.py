import biskit as B
import biskit.tools as T
import biskit.dock.hexparser as H

## generate PDBModel and model dictionary pickles for receptor and ligand
rec = B.PDBModel('rec/1A2P_clean.pdb')
rec.saveAs('rec/1A2P.model')
recpsf = B.PCRModel(source=rec, fPsf='1A2P.psf', pdbCode='1A2P')

rec_dic = {1: recpsf}
T.dump(rec_dic, 'rec/1A2P_model.dic')


lig = B.PDBModel('lig/1A19_clean.pdb')
lig.saveAs('lig/1A19.model')
ligpsf = B.PCRModel(source=lig, fPsf='1A19.psf', pdbCode='1A19')

lig_dic = {1: ligpsf}
T.dump(lig_dic, 'lig/1A19_model.dic')

## generate complex list from Hex docking result for this rec and lig
h = H.HexParser('hex/1A2P-1A19_hex.out', rec_dic, lig_dic)
c_lst = h.parseHex()

T.dump(c_lst, 'hex/complexes.cl')
