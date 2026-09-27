!  MODULE MODCCDDI
!>    @author  Taylor Harville
!
!>    @brief   Contains parameters for CCDDI COMMON Block in CC files
!
!>    @detail  Replaces the following common blocks:
!>    COMMON /CCDDI / D_VR,D_VR_IBT,D_VR_SYM,D_VR_BIT,D_VR_BI,
!>    *                D_VL,D_VL_IBT
!>    
!>    @param D_VR Address for (VO|VO) integral class for Distr. Mem.
!>    @param D_VR_IBT  Address for reordered (VO|VO) for Distr. Mem.
!>    @param D_VR_SYM  Address for symmtric (VO|VO) for Distr. Mem.
!>    @param D_VR_BIT Address for another reordered set of (VO|VO) for Distr. Mem.
!>    @param D_VR_VI Final address for last reordered set of (VO|VO) for Distr. Mem.
!>    @param D_VL Address for (VV|OO) for Distr. Mem.
!>    @param D_VL_IBT Address for reordered (VV|OO) for Distr. Mem.
 
MODULE modcc

    IMPLICIT NONE

    INTEGER :: D_VR 
    INTEGER :: D_VR_IBT
    INTEGER :: D_VR_SYM
    INTEGER :: D_VR_BIT
    INTEGER :: D_VR_BI
    INTEGER :: D_VL
    INTEGER :: D_VL_IBT

!  MODULE MODCCFILE
!>    @author  Taylor Harville
!
!>    @brief   Contains parameters for CCFILE COMMON Block in CC files
!
!>    @detail  Replaces the following common blocks:
!>      COMMON /CCFILE/ INTG,NT1,NT2,NT3,NVM,NVE,NFRLE,NRESF,NRESL
!>      INTEGER         INTG,NT1,NT2,NT3,NVM,NVE,NFRLE,NRESF,NRESL
!>
!>    @param INTG Associated with file open CCINTS
!>    @param NT1 " " CCT1AMP
!>    @param NT2 " " CCT2AMP
!>    @param NT3 " " CCT3AMP
!>    @param NVM " " CCVM (Integrals)
!>    @param NVE " " CCVE (Integrals)
!>    @param NFRLE " " CCDIIS
!>    @param NRESF " " CCREST
!>    @param NRESL " " ???

    INTEGER :: INTG,NT1,NT2,NT3,NVM,NVE,NFRLE,NRESF,NRESL

!  MODULE MODCCIKCT
!>    @author  Taylor Harville
!
!>    @brief   Contains parameters for CCPAR COMMON Block in CC files
!
!>    @detail  Replaces the following common blocks:
!>      COMMON /CCIKCT/ IKCUT
!>    @param IKCUT Value for IK cutoff

    INTEGER :: IKCUT

!>    @brief   Contains parameters for CCPAR COMMON Block in CC files
!
!>    @detail  Replaces the following common blocks:
!>    COMMON /CCPAR /  AMPTSH,METHCC,NCCTOT,NCCOCC,NCCFZC,NCCFZV,
!>    *                 MXCCIT,MXRLEIT,MWRDCC,ICCCNV,ICCRST,IDSKCC
!
!>    @param AMPTSH defines a threshold for eliminating small cluster with
!>           absolute values smaller than AMPTSH are set
!>    @param METHCC defines and holds type of CC method
!>    @param NCCTOT Total number of MOs
!>    @param NCCOCC Number of occupied MOs, including core
!>    @param NCCFZC Number core orbitals omitted from correlation
!>    @param NCCFZV Number of frozen external orbitals omitted
!>    @param MXCCIT CC max interations
!>    @param MXRLEIT CC max diis interations
!>    @param MWRDCC CC replicated memory in WORDS
!>    @param ICCCNV Convergence flag for CC
!>    @param ICCRST IREST flag for CC
!>    @param IDSKCC Disk flag for CC

    DOUBLE PRECISION :: AMPTSH
    INTEGER :: METHCC
    INTEGER :: NCCTOT
    INTEGER :: NCCOCC
    INTEGER :: NCCFZC
    INTEGER :: NCCFZV
    INTEGER :: MXCCIT
    INTEGER :: MXRLEIT
    INTEGER :: MWRDCC
    INTEGER :: ICCCNV
    INTEGER :: ICCRST
    INTEGER :: IDSKCC

!>    @brief   Contains parameters for CCRENO COMMON Block in CC files
!
!>    @detail  Replaces the following common blocks:
!>      COMMON /CCRENO/ OSS,ODS,ODD,OTS,OTD,OTT,ODS_S,ODS_D,ODS_T,
!>     *                OQS,OQDS,OQDD,OQTS,ESD,ETD,ETS,ETTM,ESD_TM
!>     Cannot figure out what the exact purpose for all these variables
!>     are, but they are some sort of variable used in the calculations
!>     of triples.
!>
!>    @param OSS
!>    @param ODS
!>    @param ODD
!>    @param OTS
!>    @param OTD
!>    @param OTT
!>    @param ODS_S
!>    @param ODS_D
!>    @param ODS_T
!>    @param OQS
!>    @param OQDS
!>    @param OQDD
!>    @param OQTS
!>    @param ESD
!>    @param ETD
!>    @param ETS
!>    @param ETTM
!>    @param ESD_TM

    DOUBLE PRECISION :: OSS,ODS,ODD,OTS,OTD,OTT,ODS_S,ODS_D,ODS_T
    DOUBLE PRECISION :: OQS,OQDS,OQDD,OQTS,ESD,ETD,ETS,ETTM,ESD_TM

!>    @brief   Contains parameters for CCPAR COMMON Block in CC files
!
!>    @detail  Replaces the following common blocks:
!>    COMMON /CCRLE / MXRLE,NRLE0,NRLE,IRLE,ITRLE
!>    @param MXRLE ???
!>    @param NRLE0 ???
!>    @param NRLE ???
!>    @param IRLE ???
!>    @param ITRLE ???

    INTEGER :: MXRLE,NRLE0,NRLE,IRLE,ITRLE

END MODULE modcc
