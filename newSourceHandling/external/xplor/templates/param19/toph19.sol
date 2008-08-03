remarks  TOPH19.SOL
remarks  ==========
remarks  topology file for solvent molecules
remarks  water models available: TIP3P model

set echo=false end

!!the acceptor and donor terms are just for analysis
!!==================================================

{* default masses *}
MASS   HT     1.00800! TIPS3P water hydrogen
MASS   OT    15.99940 ! TIPS3P water oxygen
MASS   OM    15.99940 ! oxygen in O2 and CO
MASS   CM    12.01100 ! carbon in CO
MASS   CH2E  14.02700!   
MASS   CH3E  15.03500!   
MASS   OH1   15.99940! hydroxy oxygen
MASS   H      1.00800! hydrogen which can h-bond to neutral atom



AUTOGENERATE ANGLES=TRUE END   
 
!------------------------------------------------------------------

RESIdue TIP3       { TIPS3P WATER MODEL }
 GROUp
  ATOM OH2  TYPE=OT   CHARge= -0.834  END
  ATOM H1   TYPE=HT   CHARge=  0.417  END
  ATOM H2   TYPE=HT   CHARge=  0.417  END

 BOND OH2  H1
 BOND OH2  H2
! ACCEptor OH2  " "
! DONOr H1 OH2
! DONOr H2 OH2

END {* TIP3P *}

!------------------------------------------------------------------

RESIdue O2  {* oxygen *}
 GROUp
  ATOM O1  TYPE=OM   CHARGE= 0.021  END
  ATOM O2  TYPE=OM   CHARGE=-0.021  END
 BOND O1 O2
END {* O2 *}

!-------------------------------------------------------------------

RESIdue CO  {* carbon monoxide *}
 GROUp
  ATOM C  TYPE=CM  CHARGE=0.021   END
  ATOM O  TYPE=OM  CHARGE=-0.021  END
 BOND C  O
END {* CO *}

!-------------------------------------------------------------------

RESIdue ETH  {* ethylene *}
 GROUp
  ATOM C2   TYPE=CH3E    CHARge= 0.00  END
 GROUp
  ATOM C1   TYPE=CH2E    CHARge= 0.25  END  !#
  ATOM O    TYPE=OH1     CHARge=-0.65  END  !#
  ATOM H    TYPE=H       CHARge= 0.40  END  !#

 BOND C2   C1
 BOND C1   O
 BOND O    H

 DIHEdral C2   C1   O    H

 DONOr    H   O

 ACCEptor O  " "

END {* ETH *}

!----------------------------------------------------------------------

RESIdue COH
 GROUp
  ATOM C  TYPE=CH3E   CHARge= 0.25  END
  ATOM O  TYPE=OH1    CHARge=-0.65  END
  ATOM H  TYPE=H      CHARge=0.40   END
 BOND C  O      BOND O  H
 DONOr H  O
 ACCEptor O  " "
END  {* COH *}

!----------------------------------------------------------------------

set echo=true end
