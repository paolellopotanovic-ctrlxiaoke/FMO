! MODULE MOD_SFORMAS
!>    @author  Joani Mato
!
!>    @brief   Contains old and new multiplicity in spin-flip ORMAS 
!
!>    @details  Stores a couple of multiplicities needed for an SF-ORMAS
!>              calculation. The initial high-spin and later low-spin CI. 
!    
!>    @param MULSCF The multiplicity of the underlying SCF calc.
!>    @param MULSFCI The multiplicity of the SF-CI calculation
!
MODULE mod_sformas

    IMPLICIT NONE

    INTEGER MULSCF
    INTEGER MULSFCI

!>    @author2 Katherine Ferreras
!
!>    @brief   Contains params needed to En calc: SF-ORMAS/PDFT with VVOs  
!
!>    @details params used to: screen for SF-ORMAS En calc in inputa.src.
!>             screen to proceed or block VVOs calc + substitution for 
!>             SF-ORMAS/PDFT in scflib. Will not allow En of SFORM + MP2.
!>             In vvos.src: NACVVO will stored sum of canonical occ. orbs
!>             and virtual val orbs. In ormas1.src: the correct state number
!>             is identified according to inputs in $CIDET group (kroot), and
!>             later used in mcpdft.src to print the correct state in output f.
!              
!>    @param SFOCI logical var: True for RUNTYP=ENER + CITYP=SFORM in inputa
!>    @param SFOVV logical var for VVOS in SF-ORMAS-PDFT calc stores input in
!>                 $SCF for SFOVVO in scflib
!>    @param NACVVO number of NACT when using VVOs in SF-ORMAS-PDFT calc
!>    @param KROOT state to calc DM1/2 from $CIDET. Use to label PDFT state

    LOGICAL :: SFOCI
    LOGICAL :: SFOVV
    INTEGER :: NACVVO
    INTEGER :: KROOT
   
END MODULE mod_sformas
