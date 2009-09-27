import Biskit as B
import Biskit.tools as T
import numpy as N

def extract_unique_residues( m ):
    """
    @param m: PDBModel
    @return: dict, {'name':residue_PDBModel}
    """
    d = {}
    for res in m.resModels():
        name = res['residue_name'][0]
        if not name in d:
            d[name] = res
    return d

def normalize_residue( r ):
    """
    @param r: PDBModel of single residue
    @return: PDBModel
    """
    r['chain_id'] = ['A'] * len( r )
    r['segment_id'] = [''] * len( r )
    r['residue_number'] = [ 1 ] * len( r )
    r = r.centered()
    return r

def compile_standard_dict( fpdb, fout, fmask=None, skip_first=True,
                           take_chain=True):
    """
    @param fpdb: str, pdb file containing standard residues
    @param fout: str, file name for pickled dictionary
    @param fmask: func, method to mask standard residues
    """
    m = B.PDBModel( fpdb )
    if take_chain:
        m = m.takeChains( [0] )
    if skip_first:
        m.removeRes( 0 )
    if fmask is not None:
        m = m.compress( fmask(m) )

    d = extract_unique_residues( m )
    for name in d:
        d[name] = normalize_residue( d[name] )
    T.dump( d, fout )
    return d
    

if __name__ == '__main__':

    aa  = compile_standard_dict( '3TGI', 'amino_models.dict',
                                 fmask=B.PDBModel.maskProtein )
    dna = compile_standard_dict( '1d23', 'dna_models.dict',
                                 fmask=B.PDBModel.maskDNA )
    rna = compile_standard_dict( '2Q1R', 'rna_models.dict',
                                 fmask=B.PDBModel.maskRNA )
    xnp = compile_standard_dict( 'XNP.pdb', 'xnp_models.dict',
                                 fmask=None,
                                 skip_first=False, take_chain=False)


    ## write nucleotides to separate folder
##     all = xnp.values()
##     for m in all:
##         m.writePdb( 'mononuc/%s.pdb'%m['residue_name'][0] )
