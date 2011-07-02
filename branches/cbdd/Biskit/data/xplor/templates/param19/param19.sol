remarks   PARAM19.SOL (water parameters)
remarks   ===========
remarks   available: TIPS3P model

set echo=false end


{* TIP3P model *}
{* =========== *}
BOND HT   OT     450.0       0.9572    ! from TIPS3P geometry ( SHAKE w/PARAm)
BOND HT   HT       0.0       1.5139    ! from TIPS3P geometry ( SHAKE w/PARAm)
ANGLE HT   OT   HT      55.0     104.52    ! FROM TIPS3P geometry


!	We're gonna NBFIX the TIP3P water-water interactions
!	here to make them more like Jorgensen's.  The vdW parameters
!	specified above will be in effect, therefore, for ONLY
!	protein (read, protein OR nucleic acid)-water interactions.
!	OT-OT is exactly Jorgensen's; HT interactions are added
!	here.

{* for solute-water interactions *}
NONBONDED OT        0.1591   2.8509     0.1591   2.8509  !TIP3P water oxygen
NONBONDED HT        0.0498   1.4254     0.0498   1.4254  !TIP3P water hydrogen


{* for water-water interactions *}
!!-------------------EPS-------SIGMA---------EPS14---SIGMA14-----------------
!! NONBONDED OT      0.1521   3.1506       0.1521  3.1506 !TIP3P water oxygen
!! NONBONDED HT      0.04598  0.4000       0.04598 0.4000 !TIP3P water hydr.

!! THIS SHOULD BE THE LAST NBFIX !! ANY NEW ATOM TYPE WILL USE 
!! THESE PARAMETERS FOR MIXING RULES !!
!---------------A--------------B--------------A14-----------B14-----
nbfix ot  ot  581980.4948  595.0436396     581980.4948  595.0436396
nbfix ht  ht  3.085665E-06 7.533363E-04    3.085665E-06 7.533363E-04
nbfix ht  ot  327.8404792  10.47230620     327.8404792  10.47230620 

set echo=true end
