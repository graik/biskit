import sys

import Biskit as B
import numpy as npy

def checkPresenceAll(test, reference):
    """
    This function checks that all items in test list are present 
    in reference list.
    Returns true if all elements of test are in reference.
    """
    for value in test:
        if value not in reference:
            return False
    return True

def takeProtein_noH(PDBModel, mask=False):
    m_H = npy.bitwise_not(PDBModel.maskH())
    m_prot = PDBModel.maskProtein() * m_H
    
    if mask:
        return m_prot
    else:
        return PDBModel.compress(m_prot)

def takeAmberNA(PDBModel, mask=False):
    """
    Take nucleic acid with Amber Residue naming
    If mask True, return the mask
    If false, return new PDBModel
    """
    m_H = npy.bitwise_not(PDBModel.maskH())
    maskNA = PDBModel.maskFrom('residue_name', ['RA','RU','RC','RG','RA3','RA5','RU3','RU5','RC3','RC5','RG3','RG5','DA','DA3','DA5','DG','DG3','DG5','DC','DC3','DC5','DT','DT3','DT5'])
    
    if mask:
        return m_H * maskNA
    else:
        return PDBModel.compress(m_H * maskNA)

def takeSoluteSolvent(PDBModel, noH=True, mask=True, names=None, soluteresnames=None):
    """ 
    This function will examinate a PDBModel instance
    and will try to separate the solute from the solvent.
    Then will return 2 masks corresponding to solute and solvent.
    
    It will try to extract Protein, RNA or DNA as solute
    then invert mask for the solvent.
    
    For purifying solvent, a list of residue names can be given
    under names argument.
    
    noH [0|1] Remove (True default) or not remove hydrogens
    mask[0|1] Return mask (True default) or return PDBModels (false)
    names: list or tuple with names of residues to study in the solvent
    soluteresnames: list or tuple with extra names to consider as solute. Useful when
        having some non standard residue names, glucids, fatty acids... etc.
    """
    # Mask Hydrogens and invert it
    if noH:
        m_H = npy.bitwise_not(PDBModel.maskH())
    else:
        m_H = True
        
    # Mask Protein + Nucleic Acid mask as solute
    m_prot = PDBModel.maskProtein()
    m_nucAcid = PDBModel.maskFrom('residue_name', ['RA','RU','RC','RG','RA3','RA5','RU3','RU5','RC3','RC5','RG3','RG5','DA','DA3','DA5','DG','DG3','DG5','DC','DC3','DC5','DT','DT3','DT5'])
    m_sys = m_prot + m_nucAcid
    
    if soluteresnames:
        m_extra = PDBModel.maskFrom('residue_name', soluteresnames)
        m_sys += m_extra
        
    m_solute = m_sys * m_H
    
    # Extract solvent depending on names list or 
    # make it inverse of the solute mask
    if names:
        if checkPresenceAll(names, npy.unique(PDBModel['residue_name'])):
            m_solvent = PDBModel.maskFrom('residue_name', names) * m_H
        else:
            sys.exit('ERROR: Bad residue name in arguments')
    else:
        m_solvent = npy.bitwise_not(m_sys) * m_H

    # Return mask or PDBModel instance
    if mask:
        return m_solute, m_solvent
    else:
        return PDBModel.compress(m_solute), PDBModel.compress(m_solvent)
    
def takeSolvent_noH(PDBModel, names=None, mask=False):
    """
    This function returns a compressed PDBModel or a Mask
    of the Solvent of a system (we consider the Solvent,
    all that is not the protein.
    Optionally one can give a list of names for the residues
    to extract from the solvent.
    """
    m_H = npy.bitwise_not(PDBModel.maskH())
    
    if names:
        # Check that all the names are correct
        if checkPresenceAll(names, npy.unique(PDBModel['residue_name'])) :
            m_solvent = PDBModel.maskFrom('residue_name', names) * m_H
        else:
            sys.exit('ERROR: Bad residue name in arguments')
    else:
        m_solvent = npy.bitwise_not(PDBModel.maskProtein()) * m_H
    
    if mask:
        return m_solvent
    
    else:
        return PDBModel.compress(m_solvent)