! MODULE MODMCPDFT
!>    @author  Laura Gagliardi Group
!
!>    @brief   Contains parameters for MCPDFT
!
!>    @detail  MCPDFT control variables
!
!     REVISION HISTORY:
!>    @date May, 2018 - Initial release
!
!>    @param   MCPDFTRUN - Logical that answers, Is this an MCPDFTRUN?
!>    @param   MCPDFTFUNCTIONAL - contains DFTYPE passed to DFTINP
MODULE modmcpdft

    IMPLICIT NONE

    LOGICAL :: MCPDFTRUN
    LOGICAL :: F_FLAG, T_FLAG
    DOUBLE PRECISION :: EPDFT,T_ALPHA,T_BETA
    DOUBLE PRECISION :: MCPDFTFUNCTIONAL,DFTEQ
    DOUBLE PRECISION :: PDFTYP
    DOUBLE PRECISION :: TPBE,FTPBE,TBLYP,FTBLYP,FTREVPBE,TREVPBE
    DOUBLE PRECISION :: eTPBE,eFTPBE,eTBLYP,eFTBLYP,eFTREVPBE,eTREVPBE
    DOUBLE PRECISION :: TOTELEA, TOTELEB, TOTELEAt, TOTELEBt
    DATA MCPDFTFUNCTIONAL /8HPBE     /
    integer :: DEbcou 
     
    DATA eTPBE /8HPBE     /
    DATA eFTPBE /8HPBE     /
    DATA eFTBLYP /8HBLYP    /
    DATA eTBLYP /8HBLYP    /
    DATA eTREVPBE /8HREVPBE  /
    DATA eFTREVPBE /8HREVPBE  /

    DATA TPBE /8HTPBE    /
    DATA FTPBE /8HFTPBE    /
    DATA FTBLYP /8HFTBLYP   /
    DATA TBLYP /8HTBLYP   /
    DATA TREVPBE /8HTREVPBE /
    DATA FTREVPBE /8HFTREVPBE/
END MODULE modmcpdft

! MODULE comm_MAXROOTS
!>    @author  Laura Gagliardi Group
!
!>    @brief   Contains parameters for MAXROOTS common block
!
!>    @detail  MCPDFT max roots
!
!     REVISION HISTORY:
!>    @date Sept 6, 2019 - Initial release
!
!>    @param   MXSTATE - maximum number of states to solve
MODULE comm_MAXROOTS
    IMPLICIT NONE
    INTEGER :: MXSTATE
END MODULE comm_MAXROOTS

! MODULE comm_RDM_MCPDFT
!>    @author  Laura Gagliardi Group
!
!>    @brief   Contains parameters for MCPDFT
!
!>    @detail  MCPDFT control variables
!
!     REVISION HISTORY:
!>    @date Sept 6, 2019 - Initial release
!
!>    @param   LD1A - Memory location for 1-RDM elements for active orbitals 
!>    @param   LXT4 - Memory location of 2-RDM elements for active orbitals  
!>    @param   NACTIVE - number of active orbitals 
!>    @param   NNZ - number of active+inactive orbitals
MODULE comm_RDM_MCPDFT
    IMPLICIT NONE
    INTEGER :: LD1A
    INTEGER :: LXT4
    INTEGER :: NACTIVE
    INTEGER :: NNZ
END MODULE comm_RDM_MCPDFT
