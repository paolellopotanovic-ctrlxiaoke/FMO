
      MODULE EFP_LOGICAL

!> @brief logical variables that can turn on/off terms in EFP
!>
!> @author Peng XU
!> @details
!>   when MAKEFP:
!>    IDQDYN : if true, generate LMO dipole-quadrupole polarizability
!>    IQQDYN : if true, generate LMO quadrupole-quadrupole polarizability
!>    IDODYN : if true, generate LMO dipole-octopole polarizability
!>    IDD,IDQ,IQQ,IDO : if true, generate moleclar polarizabilities
!>    $MAKEFP DISP7=.true. ==> IDQDYN = .true.
!>    $MAKEFP aniso8=.true. ==> IDQDYN,IQQDYN,IDODYN = .true.
!>    $MAKEFP DISP8 = .true. ==> IQQDYN = .true.
!>    $MAKEFP MOLPOL=.true. ==> IDD, IDQ, IQQ, IDO = .true.
!>
!>   when doing EFP-EFP calculation:
!>    E7DISP : if true, calculate E7 dispersion (default false)
!>    E8DISP : if true, calculate isotropic E8 dispersion (default false)
!>    if do not want to calculate E7disp and/or E8disp, 
!>    use NODISP7 and/or NODISP8 in $EFRAG
!>    note that if NODISP8 is present, the old disp8=disp6/3
!>    approximation is used (needed if want to reproduce old
!>    calculations)
!>
!>    do_aniso_disp : if true, calculate anisotropic disperion 
!>
!>    for dispersion energy calculated from molecular polarizability
!>    DISP_GLOBAL ==> MOLPOL=.t. 
!>
!>
      LOGICAL :: MOLPOL= .false.
      LOGICAL :: IDD=.false.,IDQ=.false.,IQQ=.false.,IDO=.false.
      LOGICAL :: DISP7 = .true., DISP8 = .true., aniso8 = .false.
      LOGICAL :: IDDDYN = .true.
      LOGICAL :: IDQDYN = .false.
      LOGICAL :: IQQDYN = .false.
      LOGICAL :: IDODYN = .false.
      LOGICAL :: E7DISP = .false.
      LOGICAL :: E8DISP = .false.
!---PX: this logicals currently are not used, the codes are commented
!       out in locpol.src
      logical :: IQD =.false.
      logical :: IOD =.false.
! for EFP-EFP calculation, turn true if ANISO is present in input
      LOGICAL :: do_aniso_disp = .false. 
 
! for MAKEFP
      LOGICAL :: MAKEFP_DYN = .false.

      END MODULE EFP_LOGICAL

      MODULE DYNPOL_MOL
!> @brief molecular dynamic polarizability (standard and rotated)
!>
!> @author Peng Xu 
!> @date July 14, 2014
!>      - initial implementation, fixed dimension sizes
!> @date Sep, 2023  Peng Xu
!>      - change to allocatable arrays
!>
!> @details DYNDD : molecular dipole-dipole polarizability tensor
!>          DYNDQ : molecular dipole-quadrupole polarizability tensor
!>          DYNQQ : molecular quadrupole-quadrupole polarizability tensor
!>          DYNDO : molecular dipole-octopole polarizability tensor
!>          
!>          DYNDD_ROT : rotated dipole-dipole pol
!>          DYNDQ_ROT : rotated dipole-quadrupole pol
!>          DYNQQ_ROT : rotated quadrupole-quadrupole pol
!>          DYNDO_ROT : rotated dipole-octopole pol

      IMPLICIT NONE
! standard molecular polarizability are allocated and deallocated
! in subroutine EFINP
! the dimension of these polarizability, e.g. DYNQQ(3,3,3,3,MXIFRQ)

      double precision,allocatable :: DYNDD(:,:,:)
      double precision,allocatable :: DYNDQ(:,:,:,:)
      double precision,allocatable :: DYNDO(:,:,:,:,:)
      double precision,allocatable :: DYNQQ(:,:,:,:,:)

! rotated molecular polarizability are allocated in EFINP (efinp.src)
! and deallocated in EFPDYN (efdrvr.src) 
! the dimension of these polarizability, 
! e.g. DYNQQ_ROT(3,3,3,3,MXIFRQ,NFRG)

      double precision,allocatable,save :: DYNDD_ROT(:,:,:,:)
      double precision,allocatable,save :: DYNDQ_ROT(:,:,:,:,:)
      double precision,allocatable,save :: DYNDO_ROT(:,:,:,:,:,:)
      double precision,allocatable,save :: DYNQQ_ROT(:,:,:,:,:,:)


      END MODULE DYNPOL_MOL

      MODULE DYNPOL_DIST
!> @brief LMO dynamic polarizability tensors (standard and rotated)
!>
!> @author Peng Xu
!> @date July 14, 2014
!>      - initial implementation, fixed dimension sizes
!>
!> @date January 2016 -E.Guidez-
!>          DYNDQ_LMO_NEW: 2D array containing dipole-quadrupole
!>          polarizability needed for optimizations and MD simulations
!>
!> @date Sep, 2023  Peng Xu
!>      - change to allocatable arrays
!>
!> @details DYNDD_LMO : LMO dipole-dipole polarizability tensor
!>          DYNDQ_LMO : LMO dipole-quadrupole polarizability tensor
!>          DYNQQ_LMO : LMO quadrupole-quadrupole polarizability tensor
!>          DYNDO_LMO : LMO dipole-octopole polarizability tensor
!> 
!>          DYNDD_LMO_ROT : rotated LMO dipole-dipole pol
!>          DYNDQ_LMO_ROT : rotated LMO dipole-quadrupole pol
!>          DYNQQ_LMO_ROT : rotated LMO quadrupole-quadrupole pol
!>          DYNDO_LMO_ROT : rotated LMO dipole-octopole pol
!>          DER_DYNDD : Rotational derivatives of dipole-dipole polarizability
!>          DER_DYNDQ: Rotational derivative of dipole-quadrupole polarizability
!>         
!> @date January 2016 -E.Guidez-
!>          DYNDQ_LMO_NEW: 2D array containing dipole-quadrupole polarizability
!>                         Needed for optimizations and MD simulations

      use mx_limits, only: MXFRG,MXIFRQ
      IMPLICIT NONE
      integer, parameter :: MXLMO=50

      DOUBLE PRECISION :: DYNDQ_LMO_NEW(27,MXIFRQ*MXLMO*MXFRG),&
                          DER_DYNDD(3,9,MXIFRQ*MXLMO,MXFRG),&
                          DER_DYNDQ(3,27,MXIFRQ*MXLMO,MXFRG)

      double precision,allocatable :: DYNDD_LMO(:,:,:)
      double precision,allocatable :: DYNDQ_LMO(:,:,:,:)
      double precision,allocatable :: DYNDO_LMO(:,:,:,:,:)
      double precision,allocatable :: DYNQQ_LMO(:,:,:,:,:)

      double precision,allocatable,save :: DYNDD_LMO_ROT(:,:,:,:)
      double precision,allocatable,save :: DYNDQ_LMO_ROT(:,:,:,:,:)
      double precision,allocatable,save :: DYNDO_LMO_ROT(:,:,:,:,:,:)
      double precision,allocatable,save :: DYNQQ_LMO_ROT(:,:,:,:,:,:)

      END MODULE DYNPOL_DIST

      MODULE DYNPOL_Intermediate
!> @brief This module contains the intermediate array to store 
!>        various dynamic polarizability tensors
!>
!> @author Peng Xu
!>  - Sep, 2019
!> @details
!>    *_IN have compact index
!>    *_OUT have expanded index
!>    *_SFT refer to LMO polarizability with the expansion points
!>          being shifted from center of mass to LMO centroids
!>
!>    DD_IN  : dyn. mol. dip-dip polarizability (9,MXIFRQ) 
!>    DD_OUT : dyn. mol. dip-dip polarizability (3,3,MXIFRQ) 
!>    DQ_IN  : dyn. mol. dip-quad polarizability (18,MXIFRQ)
!>    DQ_OUT : dyn. mol. dip-quad polarizability (3,3,3,MXIFRQ)
!>    DO_IN  : dyn. mol. dip-octo polarizability (30,MXIFRQ)
!>    DO_OUT : dyn. mol. dip-octo polarizability (3,3,3,3,MXIFRQ)
!>    QQ_IN  : dyn. mol. quad-quad polarizability (36,MXIFRQ)
!>    QQ_OUT : dyn. mol. quad-quad polarizability (3,3,3,3,MXIFRQ)
!>    DDL    : dyn. LMO dip-dip polarizability (3,3,LNA,MXIFRQ)
!>    DQL_IN : dyn. LMO dip-quad polarizability (18,LNA,MXIFRQ)
!>    DQL_OUT: dyn. LMO dip-quad polarizability (3,3,3,LNA,MXIFRQ)
!>    DQL_SFT: dyn. LMO dip-quad polarizability (3,3,3,LNA,MXIFRQ)
!>    DOL_IN : dyn. LMO dip-octo polarizability (30,LNA,MXIFRQ)
!>    DOL_OUT: dyn. LMO dip-octo polarizability (3,3,3,3,LNA,MXIFRQ)
!>    DOL_SFT: dyn. LMO dip-octo polarizability (3,3,3,3,LNA,MXIFRQ)
!>    QQL_IN : dyn. LMO quad-quad polarizability (36,LNA,MXIFRQ)
!>    QQL_OUT: dyn. LMO quad-quad polarizability (3,3,3,3,LNA,MXIFRQ)
!>    QQL_SFT: dyn. LMO quad-quad polarizability (3,3,3,3,LNA,MXIFRQ)

      Implicit None
      double precision,allocatable :: & 
      DD_IN(:,:),DD_OUT(:,:,:),&
      DQ_IN(:,:),DQ_OUT(:,:,:,:),&
      DO_IN(:,:),DO_OUT(:,:,:,:,:),&
      QQ_IN(:,:),QQ_OUT(:,:,:,:,:)

      double precision,allocatable :: &
      DDL(:,:,:,:),&
      DQL_IN(:,:,:),&
      DQL_OUT(:,:,:,:,:),&
      DQL_SFT(:,:,:,:,:),&
      DOL_IN(:,:,:),&
      DOL_OUT(:,:,:,:,:,:),&
      DOL_SFT(:,:,:,:,:,:),&
      QQL_IN(:,:,:),&
      QQL_OUT(:,:,:,:,:,:),&
      QQL_SFT(:,:,:,:,:,:)
      
      END MODULE DYNPOL_Intermediate

      MODULE EFAIDIS
!> @brief This module contains global variables necessary for
!>         computing the EFP-ab initio dispersion energy
!>
!> @author Tony Smith
!> - Jan 10, 2013
!> @details
!>    SABEFPAI : THREE-DIMENSIONAL ARRAY OF DOUBLE PRECISION VALUES,
!>              allocated in SUBROUTINE EFPAIE in source file EFPAUL,
!>              deallocated in SUBROUTINE EFAIDISP in this module.
!>              Stores values of EFP-AI overlap integrals. Dimensions
!>              are (NFRG,NA-MMCORE,MMXMO), where NFRG = number of EFPs,
!>              NA-MMCORE = number of AI valence orbitals, MMXMO =
!>              number of dispersion points (centroids) per EFP
!>    IQMDISP : flag to turn on EFP-AI dispersion calculation
!>              POSSIBLE VALUES: 0 = off, non-zero = on.
!>              Set in SUBROUTINE EFINP in source file EFINP.
!>    IDISDMP : flag to specify type of damping function
!>              POSSIBLE VALUES: -1 = no damping, 1 = overlap-based,
!>                    all other values (usually 0) = Tang-Toennies
!>              Set in SUBROUTINE EFINP in source file EFINP.
!>    IDSGRD  : flag to skip the dispersion gradient calculation
!>              POSSIBLE VALUES: 0 = gradient off, non-zero = on.
!>              Set in SUBROUTINE EFINP in source file EFINP.
!>    MMCORE : number of core orbitals. Set in RHFCL in source
!>             file RHFUHF.
!>    DISPENRG : value of the EFP-AI dispersion energy
!>    MMXMO  :  number of dispersion points (centroids) per EFP
!>    MMXBF  :  number of basis functions
!
      DOUBLE PRECISION,ALLOCATABLE,SAVE :: SABEFPAI(:,:,:)
      INTEGER,SAVE :: IQMDISP,IDISDMP,IDSGRD,MMCORE
      DOUBLE PRECISION,SAVE :: DISPENRG,MMXMO,MMXBF
      END MODULE EFAIDIS

      MODULE EFAIEXREP
!> @brief This module contains global variables for
!>         the EFP-ab initio exchange repulsion energy
!>
!> @author Peng Xu 
!> - July 1, 2015
!> @details
!> EEXCH : EFP-ab initio exchange repulsion energy
!> @date October 2018, Peng Xu and Tosaporn Sattasathuchana
!> - add IOMODE: option for TEI of QM-EFP2 exchange repulsion

      DOUBLE PRECISION :: EFAITOL=1E-7
      DOUBLE PRECISION :: EEXCH 
      INTEGER :: IOMODE
      END MODULE EFAIEXREP

      MODULE EFAIELEC
!> @brief This module contains global variables for
!>         the EFP-ab initio electrostatic energy
!>
!> @author Tosaporn Sattasathuchana and Peng Xu
!> - October 2018
!> @details
!>   IEFAI_COUL : flag to specify type of damping function
!>            -1 = No damping (default for efp2)
!>             0 = GAUSSIAN damping (default for efp1)
!>             1 = EXPONENTIAL

      INTEGER :: IEFAI_COUL
!>   IEFAI_COULHI : damping options
!>            = 0 for QM-EFP damped for charge only
!>            = 1 for QM-EFP damped at higher order terms
      INTEGER :: IEFAI_COULHI
      END MODULE EFAIELEC

!------------ EFP COMMONBLOCK -----------------------------------
      MODULE comm_FRGINF
      use mx_limits, only: MXFRG,MXDPPT
      IMPLICIT NONE
      integer :: NMPTS(MXFRG),NMTTPT,IEFC,IEFD,IEFQ,IEFO
      integer :: NPPTS(MXFRG),NPTTPT,IEFP
      integer :: NRPTS(MXFRG),NRTTPT,IREP,ICHGP,NFRG
      integer :: NDPPTS(MXDPPT),NDPTTPT,IEFDP,LSTMPTS(MXFRG)
      integer :: NBSFN(MXFRG),NMXMO(MXFRG)

      END MODULE comm_FRGINF
!
      MODULE comm_FRGSTD
      use mx_limits, only: MXFRG,MXPT,MXIFRQ
      IMPLICIT NONE
      integer :: INLPR(4*MXPT),IKFR(MXPT),IKLR(MXPT),&
                 NFPNT(MXFRG)

      CHARACTER(len=8),dimension(MXPT) ::  MNAME,PNAME,RNAME,FNAME
      CHARACTER(len=8),dimension(MXIFRQ*MXPT) :: DPNAME

      double precision :: CORD(3,MXPT),PCORD(3,MXPT),&
                          POLT(9,MXPT),&
                          DPCORD(3,MXIFRQ*MXPT),&
                          DPOLT(9,MXIFRQ*MXPT),&
                          ACHG(2,MXPT),ADIP(3,MXPT),&
                          AQAD(6,MXPT),AOCT(10,MXPT),&
                          AMULT(MXPT),BMULT(MXPT),&
                          AMULT2(MXPT),BMULT2(MXPT),&
                          RCORD(3,MXPT),TCLPR(4*MXPT),&
                          TZLPR(4*MXPT),&
                          POLSCRN(MXPT),FCRD(3,MXPT)

      END MODULE comm_FRGSTD

      MODULE comm_DOTMP
      use mx_limits, only: MXPT
      IMPLICIT NONE 
      LOGICAL :: DOMTMP(MXPT),DODTMP(MXPT)
      LOGICAL :: DOQTMP(MXPT),DOOTMP(MXPT)
      END MODULE comm_DOTMP 
!
      MODULE comm_DOMULT
      use mx_limits, only: MXFGPT
      IMPLICIT NONE
      LOGICAL :: DOMONO(MXFGPT),DODIPO(MXFGPT)
      LOGICAL :: DOQUAD(MXFGPT),DOOCTU(MXFGPT)
      END MODULE comm_DOMULT
!
      MODULE comm_EFMULT
      use mx_limits, only: MXFGPT
      IMPLICIT NONE 
      CHARACTER(len=8),dimension(MXFGPT) :: FRGNME
      double precision :: EFC(3,MXFGPT),&
                          EFCHG(2,MXFGPT),&
                          EFATRM(MXFGPT),&
                          EFBTRM(MXFGPT),&
                          EFATRM2(MXFGPT),&
                          EFBTRM2(MXFGPT),&
                          EFDIP(3,MXFGPT),&
                          EFQAD(6,MXFGPT),&
                          EFOCT(10,MXFGPT)
      END MODULE comm_EFMULT
! 
      MODULE comm_EFNRGY
      IMPLICIT NONE
      LOGICAL :: EFSTATIC,MTD1,MTD2
      double precision :: ENUCRP,REPNUC0,ENEL,EFRMS,QEFEA0(5),&
                          QEFEAP(5)
      END MODULE comm_EFNRGY
!
      MODULE comm_EFPBAS
!>    @date Jan 2020  Peng Xu and Tosaporn Sattasathuchana
!>                  - add CMAXEF
!>    @param CMAXEF   maximum coefficient of s,p,d shells 
!>                    for EFP, used for screening in 
!>                    rotated axis integral scheme for QMEFP2

      use mx_limits, only: MXFGPT,MXDFG,&
         MXGEFP,MXSHEF
      IMPLICIT NONE
      integer :: KSTREF(MXSHEF,MXDFG),KATMEF(MXSHEF,MXDFG),&
                KTYPEF(MXSHEF,MXDFG),KNGEF(MXSHEF,MXDFG),&
                KLOCEF(MXSHEF,MXDFG),KMINEF(MXSHEF,MXDFG),&
                KMAXEF(MXSHEF,MXDFG),NSHLEF(MXDFG),&
                NGSSEF(MXDFG),&
                NATEF(MXDFG),NUMEF(MXDFG),NTPATM

      double precision :: EXEF(MXGEFP,MXDFG),&
                          CSEF(MXGEFP,MXDFG),&
                          CPEF(MXGEFP,MXDFG),&
                          CDEF(MXGEFP,MXDFG),&
                          CFEF(MXGEFP,MXDFG),&
                          CGEF(MXGEFP,MXDFG),&
                          PRNAME(MXFGPT),&
                          PRCORD(3,MXFGPT),&
                          EFZNUC(MXFGPT),&
                          CMAXEF(MXGEFP,MXDFG)

      END MODULE comm_EFPBAS
!
!
      MODULE comm_EFPFRC
      use mx_limits, only: MXFRG,MXFGPT
      IMPLICIT NONE
      LOGICAL :: IFRCPNT
      CHARACTER(len=8),dimension(MXFGPT) :: FRCNME
      integer :: NTFRCPN,NFRCPNT(MXFRG),LSTFPTS(MXFRG)
      double precision :: FRCCRD(3,MXFGPT)
      double precision :: FRCTRQ(6,MXFGPT)
      double precision :: FRCTRQX(6,MXFGPT)
      END MODULE comm_EFPFRC
!
      MODULE comm_EFPMP2
      IMPLICIT NONE
      integer :: MP2TYP
      END MODULE comm_EFPMP2
!
      MODULE comm_EFPOTD
      IMPLICIT NONE
      integer :: METHPOL,ISCCHG,MOVE,nfrag
      END MODULE comm_EFPOTD
!
      MODULE comm_FGRAD
      use mx_limits, only: MXFRG,MXFGPT
      IMPLICIT NONE
      double precision :: DEF(3,MXFGPT),DEFT(3,MXFRG)
      double precision :: TORQ(3,MXFRG)
      double precision :: EFCENT(3,MXFRG),FRGMAS(MXFRG)
      double precision :: FRGMI(6,MXFRG),ATORQ(3,MXFRG)
      END MODULE comm_FGRAD
!
!
      MODULE comm_FRGMSS
      use mx_limits, only: MXFGPT,MXPT
      IMPLICIT NONE
      double precision :: FPMASS(MXPT),FMASS(MXFGPT),&
                          FPNUC(MXPT),FGNUC(MXFGPT)
      END MODULE comm_FRGMSS
!
!
      MODULE comm_PAULMO
      use mx_limits, only: MXFRG,MXFGPT,MXPT,MXIFRQ,MXDFG,&
         MXDPPT,MXPAIRS,MXGEFP,MXSHEF,&
         MXCPUEFP
      IMPLICIT NONE
      integer :: NORBEF(MXFRG),NPBF(MXFRG),NTMO
      double precision :: CENTNM(MXFGPT),CENTCD(3,MXFGPT)
      END MODULE comm_PAULMO
!
!
      MODULE comm_DNTCUT
      use mx_limits, only: MXPAIRS
      IMPLICIT NONE
      logical :: LCTcut
      double precision :: CTcut2(MXPAIRS)
      END MODULE comm_DNTCUT 
!
!
      MODULE comm_EFPAR
      IMPLICIT NONE
      integer :: LOCNFG, NLOW,NUP
      END MODULE comm_EFPAR 
!
!
      MODULE comm_EFPAR2
      use mx_limits, only:MXCPUEFP
      IMPLICIT NONE
      integer :: NFRGAL(MXCPUEFP),NLOWAL(MXCPUEFP),NUPAL(MXCPUEFP),&
                 NPPTAL(MXCPUEFP)
      END MODULE comm_EFPAR2
!
!
      MODULE comm_EFPCT
      use mx_limits, only: MXDFG,MXAO
      IMPLICIT NONE
      INTEGER :: NTCTMO,NTAO
      INTEGER :: NOCC(MXDFG),NVIR(MXDFG),&
                 NAO(MXDFG),NCTMO(MXDFG)
      DOUBLE PRECISION :: CTFOK(MXAO,MXDFG)
      END MODULE comm_EFPCT
!
!
      MODULE comm_EFPGES
      use mx_limits, only: MXFRG
      IMPLICIT NONE
      double precision :: DESFRG(6,MXFRG)
      END MODULE comm_EFPGES
!
!
      MODULE comm_EFPIO
      use mx_limits, only: MXFRG
      IMPLICIT NONE
      integer :: MULMAT (MXFRG)
      END MODULE comm_EFPIO
!
!
      MODULE comm_EFPPAR
      use mx_limits, only: MXFGPT,MXPT,MXIFRQ,MXDFG,&
         MXDPPT,MXPAIRS,MXGEFP,MXSHEF,&
         MXCPUEFP
      IMPLICIT NONE
      CHARACTER(len=8),dimension(MXFGPT) ::  POLNAM
      CHARACTER(len=8),dimension(MXIFRQ*MXFGPT) :: DPOLNAM
      integer :: IPLSCR,INDSCR
      double precision :: ENO
      double precision :: EFP(3,MXFGPT),EFPOL(9,MXFGPT),&
                          EFDP(3,MXIFRQ*MXFGPT),&
                          EFDPOL(9,MXIFRQ*MXFGPT),&
                          DIND(3,MXFGPT),DINDD(3,MXFGPT),&
                          POLSCR(MXFGPT)
      END MODULE comm_EFPPAR
!
!
      MODULE comm_EFPPBC
      IMPLICIT NONE
      integer :: ISWF
      double precision :: XBOX,YBOX,ZBOX,XCUT,YCUT,ZCUT,RCUT,SWR1,SWR2
      END MODULE comm_EFPPBC
!
!
      MODULE comm_EFPSCR
      use mx_limits, only: MXFRG
      IMPLICIT NONE
      integer :: IELECS,IDISPS,NEW_POL
      double precision :: POLAB(MXFRG),B,SCROFF
      END MODULE comm_EFPSCR
!
!
      MODULE comm_FRGLJP
      use mx_limits, only: MXFGPT,MXPT,MXFRG
      IMPLICIT NONE
      character(len=8),dimension(MXPT) :: LJNAME
      character(len=8),dimension(MXFGPT) :: ELJNAM
      integer :: NLJPTS(MXFRG),NLJTTPT,ILJP,NMLJ(2)
      double precision :: CORDLJ(3,MXPT),ELJCRD(3,MXFGPT),&
                          POTLJ(2,MXPT),ELJPOT(2,MXFGPT)
      END MODULE comm_FRGLJP
!
!
      MODULE comm_FRGTYP
      use mx_limits, only: MXFRG,MXDFG
      IMPLICIT NONE
      integer :: NDFRG,ISET(MXFRG),NAMIDX(MXDFG)
      END MODULE comm_FRGTYP
!
!
      MODULE comm_IEACAL
      IMPLICIT NONE
      logical :: DOEFPIEA
      integer :: IEACAL,NFGIEA,NFRGIEA,NFRGI1,NFRGI2,IEADA,&
                 MXEFPT,MAXPOL,NPRIEA,NATFIEA,IEABDY,ITRLVL,&
                 LIEACNT,ieaold
      double precision :: REMPOL
      END MODULE comm_IEACAL
!
!
      MODULE comm_POLTMP
      use mx_limits, only: MXFRG,MXPT
      IMPLICIT NONE
      integer ::JNLPR(4*MXPT),JKFR(MXPT),JKLR(MXPT),&
                NDPTS(MXFRG),NDTTPT
      double precision :: TCLPD(4*MXPT),TZLPD(4*MXPT)
      END MODULE comm_POLTMP
!
!
      MODULE comm_PCMEFP
      IMPLICIT NONE
      double precision :: PB_MUL,PC_MUL,PB_POL,PC_POL,P_FF,P_NF,P_FN,&
                          P_I,P_IBIS,P_J,P_NUCC,P_NUCCBIS,ENPCM
      END MODULE comm_PCMEFP
!
!
      MODULE comm_FRGNRG
      IMPLICIT NONE
      double precision :: CCNRG,CDNRG,CQNRG,CONRG,DDNRG,DQNRG,QQNRG,&
                          ESNRG,REPNRG,POLNRG,DISNRG,DISNRG8,&
!                          EDISD,CHTNRG,EDISP7,EDISD6

                          EDISP,CHTNRG,EDISP7,EDISP6,&
                          EDISP8,EDISP8_aniso
      END MODULE comm_FRGNRG
!
!
      MODULE comm_EFPFMO
      IMPLICIT NONE
      integer :: IEFPFMO,ISWENUC
      double precision :: RUNEFP,REFFMO,REPNUCEFP,EFPEFP
      END MODULE comm_EFPFMO
!
!
      MODULE comm_EFPGCX
      use mx_limits, only: MXFRG
      IMPLICIT NONE
      double precision :: DCTFRG(6,MXFRG),DXRFRG(6,MXFRG)
      END MODULE comm_EFPGCX
!
!
      MODULE comm_EFPFM2
      IMPLICIT NONE
      integer :: LABFLDT,LDINDF2,LDINDDF2,LDELABFLD,&
                 LDELDIND,LDELDINDD,LALCNT2F
      double precision :: TOTEFPC,POLTOT
      END MODULE comm_EFPFM2
!
!
      MODULE comm_EFPTDG
      IMPLICIT NONE
      integer :: ITDEFG
      END MODULE comm_EFPTDG
!
!
      MODULE comm_EFPNFO
      IMPLICIT NONE
      INTEGER :: EFPMUL
      END MODULE comm_EFPNFO
      MODULE comm_EFPPRT
! @brief   This module contain keyword for printing control
!          for MAKEFP currently
! @author  Peng Xu
!         - Mar 2020
! @param   NPRTEF = 0 (default)
!          NPRTEF = -5 reduced printing
      IMPLICIT NONE
      INTEGER :: NPRTEF = 0
      END MODULE comm_EFPPRT    
!
      MODULE MAKEFP_CPHF
      IMPLICIT NONE
      logical :: cpmakefp=.false.
      logical :: makefp_print=.true.
      logical :: efmo_minout=.false.
      logical :: MAKEFP_DIRSCF=.true.
      logical :: makefp_print_default=.true.
      
      END MODULE MAKEFP_CPHF
!
!> @brief This routine is to fix the different declaration between
!> 1-D and 2-D matrix in FGRAD common block by converting 2-D AMAT to
!> 1-D.
!>
!> @author Tosaporn Sattasathuchana
!>         - March 2019
!> @details computing DDOT and finding maximum elements of the matrix
!>          this subroutine is used in GRADX routine.
!> @param   AMAT : 1-D matrix
!> @param   FMAX1 : maximum elements of the matrix
!> @param   DAMAT : dot product of AMAT*AMAT
!
      SUBROUTINE FMAXFGRAD(AMAT,NCF,FMAX1,DAMAT)
      use mx_limits, only: MXFRG
      IMPLICIT NONE
      double precision :: AMAT(3*MXFRG),FMAX1,DAMAT,DDOT
      integer :: NCF, IDAMAX
      FMAX1 = ABS(AMAT(IDAMAX(NCF,AMAT,1)))
      DAMAT = DDOT(NCF,AMAT,1,AMAT,1)
      RETURN
      END

      module makefp_timing
      implicit none
      private
      public makefp_timing_cpu
      contains
      subroutine makefp_timing_cpu(tstim)
      USE params, ONLY: intomp
      implicit none
      double precision, intent(inout) :: tstim
! MPI timing
      if(intomp.lt.1) then
          CALL TSECND(tstim)
      endif
! omp timing
!$      if(intomp.eq.1) then
!$          call makefp_timing_cpu_omp(tstim)
!$      endif
      return
      end subroutine makefp_timing_cpu
      
!$      subroutine makefp_timing_cpu_omp(tstim)
!$      USE omp_lib
!$      implicit none
!$      double precision, intent(inout) :: tstim
!$      tstim     = omp_get_wtime()
!$      return
!$      end subroutine makefp_timing_cpu_omp
      end module makefp_timing
