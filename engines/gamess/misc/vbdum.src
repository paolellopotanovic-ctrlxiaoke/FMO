! 24 MAR 07 - DUMMY FILE
!***********************************************************************
!     THIS IS A DUMMY MODULE TO SHORT CIRCUIT VB2000 OR XMVB CALLS,
!     WHEN VB2000 -OR- XMVB ARE NOT LINKED INTO GAMESS.  NOTE THAT
!     IT IS IMPOSSIBLE TO HOOK UP BOTH VALENCE BOND PROGRAMS AT THE
!     SAME TIME, SINCE THEY SHARE ONE (1) DUMMY FILE, AND ONE (1)
!     NAME FOR THE INTERFACE, NAMELY THIS ROUTINE'S NAME, VBGMS.
!
!     FOR THE REAL SOURCE CODE OF VB2000 (UP TO VERSION 3.0), SEE
!     PROF. RODRIGO BITZER, BITZER@IQ.UFRJ.BR     
!     FOR THE REAL SOURCE CODE OF XMVB, SEE PROFESSOR WEI WU,
!     XIAMEN UNIVERSITY, WEIWU@XMU.EDU.CN OR http://xmvb.org/
!
!     ORIGINAL FORM WRITTEN BY: BRIAN DUKE, NOVEMBER 2003
!     Updated by David Sousa, Sep 10, 2021
!***********************************************************************
!*MODULE VBDUM  *DECK VBGMS
!>
!> @brief   interface between VB2000 / XMVB and GAMESS 
!>
!> @author  Brian Duke
!>
!> @date September, 2021 - David Sousa
!> - adjusted to GAMESS new code standards
!> - included PRINT_NOS routine from parley.src
!>
      SUBROUTINE VBGMS(IOPT)
      IMPLICIT NONE
      INTEGER IOPT
      LOGICAL GOPARR,DSKWRK,MASWRK
      INTEGER IR,IW,IP,IS,IPK,IDAF,NAV,IODA
      INTEGER ME,MASTER,NPROC,IBTYP,IPTIM
      DOUBLE PRECISION SCFTYP,VBTYP,DFTYPE,TDDFTYP,CITYP,CCTYP
      INTEGER MPLEVL,MPCTYP
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /WFNOPT/ SCFTYP,VBTYP,DFTYPE,TDDFTYP,CITYP,CCTYP,
     *                MPLEVL,MPCTYP
!
!        USER SHOULD NOT SELECT VBTYP=VB2000 OR VBTYP=XMVB WITHOUT
!        HAVING LINKED THE REAL CODE FOR ONE -OR- THE OTHER!
!
!        IN THE REAL CODE, CALLING VBGMS(0) MUST READ THE INPUT DATA,
!        WHILE VBGMS(1) DOES THE ACTUAL VALENCE BOND COMPUTATION.
!
!        SINCE WE ALWAYS CALL VBGMS(0) FIRST, WE NEVER HIT THE -ELSE-.
!
      IF(IOPT.EQ.0) THEN
         IF(MASWRK) WRITE(IW,900) VBTYP,VBTYP
         CALL ABRT
      ELSE
         IF(MASWRK) WRITE(IW,*) 'VBGMS CALLED WITH IOPT=',IOPT
         CALL ABRT
      END IF
      RETURN
  900 FORMAT(1X,'--- ERROR: VBTYP=',A8,' REQUESTED IN $CONTRL ----'/
     *       1X,'BUT THE ',A8,' CODE WAS NOT LINKED INTO GAMESS')
      END

!*MODULE VBDUM  *DECK PRINT_NOS
!>
!> @brief   print natural orbitals of the VB function
!>           
!> @details Routine from vb2gms.src called in parley.src.
!>          It is included here just to not get
!>          'undefined reference' error when linking the GAMESS
!>          executable without VB2000     
!>
!> @author  David Sousa, September 2021
!>
      SUBROUTINE PRINT_NOS(VEC,OCCNO,M,L1)
      IMPLICIT NONE
      INTEGER M,L1
      DOUBLE PRECISION VEC(L1,L1),OCCNO(L1)
      RETURN
      END

!*MODULE VBDUM  *DECK TRUCON
!>
!> @brief   apply TRUDGE constraints to coordinates in VB2000
!>           
!> @details Routine from vb2gms.src called in trudge.src.
!>          It is included here just to not get
!>          'undefined reference' error when linking the GAMESS
!>          executable without VB2000     
!>
!> @author  Brian Duke
!>
      SUBROUTINE TRUCON
      RETURN
      END

!*MODULE VBDUM  *DECK mod_vb2000
!>
!> @brief   replace COMMON block /VBINTF/
!
!> @details Replaces the following common block:
!>          COMMON /VBINTF/ VBENGY,LPROP,NOSYMV,MAXOCC,LNOS,JSTEP,MGMS
!>          It is used in vb2000.src, vb2gms.src,
!>          gamess.src, inputa.src,parley.src, and statpt.src
!>           
!> @author  David Sousa, Sep 2021
!>
      module mod_vb2000
      implicit none
      double precision :: vbengy
      integer          :: lprop, nosymv, maxocc, lnos, jstep, mgms
      end module mod_vb2000
