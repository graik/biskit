!------------------------------------------------------------------
! last $Date$
! last $Author$
!------------------------------------------------------------------

!------------------------------------------------------------------
! Xplor uses residue name CBX for N-methyl capping but normally
! this residue is called NME.
! Create residue called NME with identical topology
!------------------------------------------------------------------

RESIdue NME                 { can be linked   *     NME 
                                               \PEPT/     }
 GROUp
  ATOM N    TYPE=NH1   CHARge=-0.35   END
  ATOM H    TYPE=H     CHARge= 0.25   END
  ATOM CA   TYPE=CH3E  CHARge= 0.10   END

 BOND N    CA
 BOND N    H

 DONOr H    N

END {NME}


!------------------------------------------------------------------
! Acetyl and N-methyl capping residues should not receive terminal patch
!------------------------------------------------------------------

PRESidue NACE                { empty patch for N-teriminal ACE capping }
END {NACE}

PRESidue CCBX                { empty patch for C-teriminal CBX capping }
END {CCBX}

PRESidue CNME                { empty patch for C-teriminal NME capping }
END {CNME}
