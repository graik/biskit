import biskit as B
import biskit.tools as T
import biskit.dock.hexparser as H

rec = B.XplorModel(source='rec/1A2P_clean.pdb', fPsf='1A2P.psf', pdbCode='1A2P')
rec.saveAs('rec/1A2P.model')

rec_dic = {1: rec}
T.dump(rec_dic, 'rec/1A2P_model.dic')


lig = B.XplorModel(source='lig/1A19_clean.pdb', fPsf='1A19.psf', pdbCode='1A19')
lig.saveAs('lig/1A19.model')

lig_dic = {1: lig}
T.dump(lig_dic, 'lig/1A19_model.dic')

## generate complex list from Hex docking result for this rec and lig
h = H.HexParser('hex/1A2P-1A19_hex5.out', rec_dic, lig_dic)
c_lst = h.parseHex()

T.dump(c_lst, 'hex/complexes.cl')
