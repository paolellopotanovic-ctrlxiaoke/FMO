#ifdef __OFFLOAD

#define __cpu 0

#ifdef __INTEL_LLVM_COMPILER
#define __intel_offloading 1
#else
#define __intel_offloading 0
#endif

#ifdef USE_CUBLAS
#define __cublas 1
#else
#define __cublas 0
#endif

#ifdef USE_HIPBLAS
#define __hipblas 1
#else
#define __hipblas 0
#endif

#else

#define __cpu 1
#define __cublas 0
#define __hipblas 0
#define __intel_offloading 0

#endif

#define __PRINT_TIME_BUG 0
#define __PRINT_ESS 1
#define __PRINT_INT 0
#define __notes 0
#define __readRHF 0

#if __intel_offloading==1
   include "mkl_omp_offload.f90"
#endif


#ifndef __OFFLOAD
      subroutine GPUrimp2driver
         COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
         logical :: MASWRK
         IF(MASWRK) WRITE(*,*) "ri-mp2 energy need GMS_OPENMP_OFFLOAD=true in install.info"
         call ABRT
      end 
#else

      module blasHandles
#if __cublas==1
        use cublasf
        integer(c_int) :: cublas_return
        type(c_ptr)    :: cublas_handle
#endif
#if __hipblas==1
        use iso_c_binding
        use hipfort
        use hipfort_check
        use hipfort_hipblas
        type(c_ptr) :: hipblas_handle = c_null_ptr
        integer*4 ::  m ,n,k
        integer*4 :: lda, ldb, ldc
        double precision,parameter ::  alpha=1.0D00, beta=0.0D0
        type(c_ptr) :: da = c_null_ptr, db = c_null_ptr, dc = c_null_ptr
#endif
      end module blasHandles !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!>*module gpurimp2   *deck gpurimp2
!>
!>     @brief   driver for rimp2 energy on gpus
!>     @author  buu
!>     @date    june 7, 2020
!>
!>     @detail  driver for rimp2 energy on gpus.
!>
      SUBROUTINE GPURIMP2DRIVER

      use omp_lib

      use gpuGlobalRimp2Eng,only: MXAUXSH,MXAXGTOT
      use gpuGlobalRimp2Eng,only: MXSH,MXGTOT,MXATM,MXRT
      use gpuGlobalRimp2Eng,only: MaxAUXANG,MaxAUXBFX
      use gpuGlobalRimp2Eng,only: DAux,DAtm
      use gpuGlobalRimp2Eng,only: LOCSPH
      use gpuGlobalRimp2Eng,only: FlgSphAux,BUG
      use gpuGlobalRimp2Eng,only: MIXS,TLL
      use gpuGlobalRimp2Eng,only: LL1,LL2,LL3,LL4,D_B
      use gpuGlobalRimp2Eng,only: NGPU

      use shellProcessing,only: SORTSHELLS,contrCoeffs,deallocateShellProcessing
      use shellProcessing,only: n_S_shl,n_P_shl,n_D_shl,n_F_shl,n_G_shl
      use shellProcessing,only:n_aux_S_shl,n_aux_P_shl,n_aux_D_shl,n_aux_F_shl,n_aux_G_shl



      implicit none

      logical :: MASWRK, EXTCAB
      integer,parameter :: NORIMP=8,NOAUBF=2
      double precision,parameter :: P12=1.2D00
      integer :: DDI_NP,DDI_ME,DDI_NN,DDI_MY
      integer :: KQRIMP(NORIMP),KQAUXBF(NOAUBF)
      double precision :: QRIMP(NORIMP),QAUXBF(NOAUBF)

      double precision,allocatable :: DE_DDI(:,:)

      COMMON /AUXBAS/ EXAUX(MXAXGTOT),CAUXS(MXAXGTOT),CAUXP(MXAXGTOT),  &
                      CAUXD(MXAXGTOT),CAUXF(MXAXGTOT),CAUXG(MXAXGTOT),  &
                      CAUXH(MXAXGTOT),CAUXI(MXAXGTOT),                  &
                      KAUXST(MXAUXSH),KAUXAT(MXAUXSH),KAUXTY(MXAUXSH),  &
                      KAUXNG(MXAUXSH),KAUXLO(MXAUXSH),KAUXMI(MXAUXSH),  &
                      KAUXMX(MXAUXSH),NAUXSH
      double precision :: &
         EXAUX(MXAXGTOT),CAUXS(MXAXGTOT),CAUXP(MXAXGTOT), &
         CAUXD(MXAXGTOT),CAUXF(MXAXGTOT),CAUXG(MXAXGTOT),CAUXH(MXAXGTOT),CAUXI(MXAXGTOT)
      integer :: &
         KAUXST(MXAUXSH),KAUXAT(MXAUXSH),KAUXTY(MXAUXSH), &
         KAUXNG(MXAUXSH),KAUXLO(MXAUXSH),KAUXMI(MXAUXSH),KAUXMX(MXAUXSH),NAUXSH

      COMMON /NSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),      &
                      CF(MXGTOT),CG(MXGTOT),CH(MXGTOT),CI(MXGTOT),      &
                      KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),   &
                      KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
      double precision :: &
         EX(MXGTOT),CS(MXGTOT),CP(MXGTOT), &
         CD(MXGTOT),CF(MXGTOT),CG(MXGTOT),CH(MXGTOT),CI(MXGTOT)
      integer :: &
         KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH), &
         KNG(MXSH),KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL

      COMMON /BASSPH/ QMTTOL,ISPHER
      double precision :: QMTTOL
      integer :: ISPHER

      COMMON /ENRGMP/ EMP2,EMP3,EMP4,EMP2A
      double precision :: EMP2,EMP3,EMP4,EMP2A

      COMMON /ENRGYS/ ENUC,EELEC,ETOT,SZ,SZZ,ECORE,ESCF,EERD,E1,E2,     &
                      VEN,VEE,EPOT,EKIN,ESTATE(MXRT),STATN,EDFT(2),EDISP
      double precision :: ENUC,EELEC,ETOT,SZ,SZZ,ECORE,ESCF,EERD,E1,E2,     &
                      VEN,VEE,EPOT,EKIN,ESTATE(MXRT),STATN,EDFT(2),EDISP

      COMMON /IOFILE/ IR,IW,IP,IJK,IPK,IDAF,NAV,IODA(950)
      integer :: IR,IW,IP,IJK,IPK,IDAF,NAV,IODA(950)

      COMMON /MACHIN/ NWDVAR,MAXFM,MAXSM,LIMFM,LIMSM
      integer :: NWDVAR,MAXFM,MAXSM,LIMFM,LIMSM

      COMMON /MP2PAR/ OSPT,CODEMP,SCSPT,TOL,METHOD,NWDMP2,MEMPRI,MPPROP,&
                      NACORE,NBCORE,NOA,NOB,NO,NBF,NOMIT,MOCPHF,MAXITC
      double precision :: OSPT,CODEMP,SCSPT,TOL,METHOD
      integer :: NWDMP2,MEMPRI,MPPROP,NACORE,NBCORE,NOA,NOB, &
                  NO,NBF,NOMIT,MOCPHF,MAXITC
      COMMON /ONEELC/ EONE,E1A,E1B
      double precision :: EONE,E1A,E1B
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      integer :: ME,MASTER,NPROC,IBTYP,IPTIM
      logical :: DSKWRK, MASWRK, GOPARR
      COMMON /RESTAR/ TIMLIM,IREST,NREC,INTLOC,IST,JST,KST,LST      
      integer :: IREST,NREC,INTLOC,IST,JST,KST,LST
      double precision :: TIMLIM
      COMMON /RUNOPT/ RUNTYP,EXETYP,NEVALS,NGLEVL,NHLEVL
      double precision :: RUNTYP,EXETYP
      integer :: NEVALS,NGLEVL,NHLEVL
      COMMON /SCFWFN/ AROHF(3),BROHF(3),PACAVO(6),IACAVO,IUHFNO,ICUHF,MVOQ
      double precision :: AROHF(3),BROHF(3),PACAVO(6)
      integer :: IACAVO,IUHFNO,ICUHF,MVOQ

      COMMON /WFNOPT/ SCFTYP,VBTYP,DFTYPE,TDDFTYP,CITYP,CCTYP,          &
                      MPLEVL,MPCTYP
      double precision :: SCFTYP,VBTYP,DFTYPE,TDDFTYP,CITYP,CCTYP
      integer :: MPLEVL,MPCTYP
      COMMON /RIMPFI/ IFILV,IFILT2A,IFILT3A,IFILT2B,IFILT3B
      integer :: IFILV,IFILT2A,IFILT3A,IFILT2B,IFILT3B
      COMMON /OUTPUT/ NPRINT,ITOL,ICUT,NORMF,NORMP,NOPK
      integer :: NPRINT,ITOL,ICUT,NORMF,NORMP,NOPK
      COMMON /FUNCT / E,EG(3*MXATM)
      double precision :: E,EG(3*MXATM)
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,                    &
                      ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      integer :: NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,IAN
      double precision :: ZAN,C


double precision :: RHF,UHF,CHECK,RIMP2,AUXBAS,empty
      DATA RHF      /8HRHF     /
      DATA UHF      /8HUHF     /
      DATA CHECK    /8HCHECK   /

      DATA RIMP2    /8HRIMP2   /
      DATA AUXBAS   /8HAUXBAS  /
      DATA QRIMP    /8HIAUXBF  ,8HIVMTD   ,8HVTOL    ,8HSTOL    ,       &
                     8HOTHAUX  ,8HGOSMP   ,8HMEMSH   ,8HUSEDM   /
      DATA QAUXBF   /8HCABNAM  ,8HEXTCAB  /
      DATA KQRIMP   /1,1,3,3,0,0,1,0/
      DATA KQAUXBF  /5,0/
      data empty    /8H        /

      character(1) :: ADDR(1)
      integer :: SMP_NP,SMP_ME
#if __hipblas==1 
      ! hipfc requires "target" attribute.
      DOUBLE PRECISION,allocatable,target :: B32(:)
#else
      DOUBLE PRECISION,allocatable :: B32(:)
#endif

      LOGICAL :: NODEMEM
      logical :: isNodeMem
      double precision :: BDUM(1)

integer :: NBASIS,NOCCA,NVIRA,NORBA,NOCCB,NVIRB,NORBB,NCOR,NACT,NVIR,NORB,IAUXBF,IVMTD
double precision :: STOL,VTOL,CABNAM, CPUMEM,sizeFund,GOSMP,MEMSH,USEDM, w1
logical :: OTHAUX
integer :: NAUXCAT,NAUXBAS,NAUXSPH,MaxATMANG,MaxATMBFX,NANGM,MaxNGs,w0,NNODE,NAUXBASD,D_B_SMP,LB32,JRET,MYNODE

      ! LOGICAL :: LL1
      ! LOGICAL :: LL2
      ! LOGICAL :: LL3
      ! LOGICAL :: LL4

      ! double precision :: TLL(15,9,25)
      ! integer :: MIXS(2,25)

                          !!! ~ GO ~ !!!

! call RIMP2_reset_gpu()


      ! turn off debug flag
      BUG=.true.

      ! some mp2 vars
      NBASIS = NBF
      NOCCA  = NOA - NACORE
      NVIRA  = NO  - NOA
      NORBA  = NO  - NACORE
      NOCCB  = NOB - NBCORE
      NVIRB  = NO  - NOB
      NORBB  = NO  - NBCORE

      NBASIS = NBF

      NCOR = NACORE
      NACT = NOA - NACORE
      NVIR = NO - NOA
      NORB = NO


      NGPU = 2


      ! init and read $rimp2 group vars
      IAUXBF= 0
      IF(ISPHER.EQ.1.OR.ISPHER.EQ.0) IAUXBF=1
      ! OTHAUX= .FALSE.
      STOL  = 1.0D-6
      IVMTD = 0
      VTOL  = 1.0D-6
      MEMSH = 0

      FlgSphAux = .FALSE.
      IF(IAUXBF.EQ.1) FlgSphAux = .TRUE.

      CALL NAMEIO(IR,JRET,RIMP2,NORIMP,QRIMP,KQRIMP,                    &
          IAUXBF,IVMTD,VTOL,STOL,OTHAUX,GOSMP,MEMSH,USEDM,              &
          0,                                                            &
          0,0,0,0,0,  0,0,0,0,0,   0,0,0,0,0,                           &
          0,0,0,0,0,  0,0,0,0,0,   0,0,0,0,0,  0,0,0,0,0,               &
          0,0,0,0,0,  0,0,0,0,0,   0,0,0,0,0,  0,0,0,0,0)
      IF(JRET .EQ. 2) CALL ABRT


      ! ao and aux bases has to be both cartesian or spherical
      IF(IAUXBF.EQ.1 .AND. ISPHER.NE.1) CALL ABRT

      ! read $auxbas group
      CABNAM=EMPTY
      EXTCAB=.FALSE.
      CALL NAMEIO(IR,JRET,AUXBAS,NOAUBF,QAUXBF,KQAUXBF,                 &
          CABNAM,EXTCAB,                                                &
          0,0,0,0,0,  0,0,                                              &
          0,0,0,0,0,  0,0,0,0,0,   0,0,0,0,0,                           &
          0,0,0,0,0,  0,0,0,0,0,   0,0,0,0,0,  0,0,0,0,0,               &
          0,0,0,0,0,  0,0,0,0,0,   0,0,0,0,0,  0,0,0,0,0)
      IF(JRET .EQ. 2) CALL ABRT



      ! read aux basis
      CALL RICAUXBAS(CABNAM,EXTCAB,NAUXCAT)
      NAUXBAS = NAUXCAT

      ! read device group
      CALL DEVICEINP()
      
#if __PRINT_TIME_BUG==1
   WRITE(*,*) "debug_debug_debug entering SPHAUX"
#endif

      ! setup spherical aux basis
      IF(FlgSphAux) THEN
         ALLOCATE(LOCSPH(NAUXSH))
         CALL SPHAUX(LOCSPH,NAUXSPH)
         NAUXBAS = NAUXSPH
      ENDIF
#if __PRINT_ESS==1
   IF(MASWRK) write(*,'(A30,2I15)') "nums AO and AUX basis functions", &
                                     NBASIS,NAUXBAS
#endif  


#if __PRINT_TIME_BUG==1
   WRITE(*,*) "debug_debug_debug entering setting up some parameters"
#endif

      ! max ang momentum of ao and aux bases
      MaxAUXANG=MAXVAL(KAUXTY,1)-1
      MaxAUXBFX=(MaxAUXANG+1)*(MaxAUXANG+2)/2
      MaxATMANG=MAXVAL(KTYPE,1)-1
      MaxATMBFX=(MaxATMANG+1)*(MaxATMANG+2)/2
      IF(MaxAUXANG.GT.4.OR.MaxATMANG.GT.4) THEN
         WRITE(*,*) "AUX/AO BASIS WITH ANGULAR MOMENTUM", &
                    "GREATER THAN 4 (H AND ABOVE) HAS NOT", &
                    "BEEN PROGRAMED. THE CALCULATION WILL BE ABORTED"
         CALL ABRT
      ENDIF


      ! max basis functions in an ao and aux shell
      NANGM=MaxAUXBFX
      IF(MaxATMBFX.gt.MaxAUXBFX) NANGM=MaxATMBFX

      MaxNGs=MAXVAL(KAUXNG,1)
      ALLOCATE(DAux(NAUXSH,MaxNGs,28))
      MaxNGs=MAXVAL(KNG,1)
      ALLOCATE(DAtm(NSHELL,MaxNGs,28))

      
#if __PRINT_TIME_BUG==1      
   WRITE(*,*) "debug_debug_debug read density factor"
#endif

      ! contraction coefficients
      CALL DensFact_gpu()

      ! GET SPHERICAL MATRIX TRANSFORMATION
      CALL RIMP2CSTRM_gpu(MIXS,TLL)

      ! SORT AO AND AUX SHELLS
      CALL SORTSHELLS(NSHELL,KTYPE,KNG,.TRUE.)
#if __PRINT_ESS==1
   if(MASWRK) write(*,'(A30,5I10)') "num of S,P,D,F,G AO shells:", &
                                     n_S_shl,n_P_shl,n_D_shl,n_F_shl,n_G_shl
#endif

      CALL SORTSHELLS(NAUXSH,KAUXTY,KAUXNG,.FALSE.)
#if __PRINT_ESS==1
   if(MASWRK) write(*,'(A30,5I10)') "num of S,P,D,F,G AUX shells:", &
                                     n_aux_S_shl,n_aux_P_shl,n_aux_D_shl,n_aux_F_shl,n_aux_G_shl
#endif

      ! get contraction coefficients
      CALL contrCoeffs(EX,CS,CP,CD,CF,CG,CH,CI, &
                       KSTART,KATOM,KTYPE,KNG,  &
                       KLOC,KMIN,KMAX,NSHELL,   &
                       NORMF,NORMP,NOPK)

#if __PRINT_TIME_BUG==1
   IF(BUG) WRITE(*,*) "debug_debug_debug entering RIMP2_X_ENERGY_gpu"
#endif


#if __PRINT_ESS==1
   w0=omp_get_wtime()
#endif
      
! EACH MPI RANK IS ASSOCIATED WITH 1 GPU
! LL1
!    SERIAL MPI RUN.
!    B-MATRIX IS ALLOCATED ON NODE MEMORY.
! LL2
!    PARALLEL MPI RUN ON SINGLE NODE
!    B-MATRIX IS ALLOCATED ON SMP NODE MEMORY
! LL3
!    PARALLEL MPI RUN ACROSS MULTIPLE NODES
!    B-MATRIX SIZE IS SMALLER THAN NODE MEMORY
! LL4
!    PARALLEL MPI RUN ACROSS MULTIPLE NODES
!    B-MATRIX SIZE LARGER THAN SINGLE NODE MEMORY.

! number of nodes and mpi ranks
CALL DDI_NNODE(NNODE,MYNODE)

! CPU memory
CPUMEM = 500

! CPU storage for fundamental arrays
! B32(NAUXBAS,NVIR,NACT) + QVV(NVIR,NVIR) + BI(NAUXBAS,NVIR)
sizeFund = (NAUXBAS*NVIR*NACT + NVIR*NVIR + NAUXBAS*NVIR)*8.0/1073741824.0

#if __PRINT_ESS==1
IF(MASWRK) write(*,*) "total memory needed: ", sizeFund, "(GB)"
#endif

! check if a node memory can hold all fundamental arrays
isNodeMem = sizeFund .LT. CPUMEM

! separate some cases
LL1 = isNodeMem .and. (NPROC .eq. 1)
LL2 = isNodeMem .and. (NNODE .eq. 1)
LL3 = NNODE .gt. 1

! write some message about used memory models
#if __PRINT_ESS==1
IF(MASWRK) THEN
   IF(LL1) THEN
      write(IW,*) "LL1: 1 mpi rank, B32 is in allocatable node memory"
   ELSEIF(LL2) THEN
      write(IW,*) "LL2: 1 node, B32 is in SMP array"
   ELSEIF(LL3) THEN
      write(IW,*) "LL3: more than 1 nodes, B32 is in DDI array"
   ELSE
      write(IW,*) "wait, what's going on?"
   ENDIF
ENDIF

FLUSH(IW)
#endif

      ! SIZE=NAUXBAS*NVIR*NACT*8.0/(NPROC*1073741824.0)
      ! NODEMEM=500.0.GT.SIZE

      ! IF(MASWRK) THEN
      !    IF(NODEMEM) THEN
      !       WRITE(*,*) "NODE MEMORY CAN HANDLE MATRIX B"
      !    ELSE
      !       WRITE(*,*) "NODE MEMORY IS SMALLER THAN MATRIX B"
      !       WRITE(*,*) "THEREFORE NEED TO USE DDI ARRAY"
      !    ENDIF
      ! ENDIF

      ! LL1=NPROC.EQ.1
      ! LL2=NPROC.GT.1 .AND. NPROC.LE.6
      ! LL3=NPROC.GT.6 .AND. SIZE.LE.500.0
      ! LL4=NPROC.GT.6 .AND. SIZE.GT.500.0

      ! IF(MASWRK) THEN
      !    IF(LL1) write(*,*) "LL1: NPROC = 1"
      !    IF(LL2) write(*,*) "LL2: 1 < NPROC <= 6"
      !    IF(LL3) write(*,*) "LL3: NPROC > 6 && MEM < NODEMEM"
      !    IF(LL4) write(*,*) "LL4: NPROC > 6 && MEM > NODEMEM"
      ! ENDIF


      IF(LL1) THEN
         ALLOCATE(B32(NAUXBAS*NVIR*NACT))
         CALL RIMP2_X_ENERGY_gpu &
             (B32,NANGM,NAUXBASD,NAUXBAS,NAUXSH, &
              NCOR,NACT,NVIR,NORB,NBF, &
              STOL,VTOL,IVMTD,OTHAUX)
         DEALLOCATE(B32)
      ELSEIF(LL2) THEN
         CALL DDI_SMP_NPROC(SMP_NP,SMP_ME)
         CALL DDI_SMP_CREATE(NAUXBAS*NVIR*NACT,D_B_SMP)
         CALL DDI_SMP_OFFSET(D_B_SMP,ADDR,LB32)
         LB32=LB32+1
         CALL RIMP2_X_ENERGY_gpu &
             (ADDR(LB32),NANGM,NAUXBASD,NAUXBAS,NAUXSH,                    &
              NCOR,NACT,NVIR,NORB,NBF,                                     &
              STOL,VTOL,IVMTD,OTHAUX)
         CALL DDI_SMP_DESTROY(D_B_SMP)
      ELSEIF(LL3) THEN
         CALL DDI_CREATE(NAUXBAS*NVIR,NACT,D_B)
         CALL RIMP2_X_ENERGY_gpu &
             (BDUM,NANGM,NAUXBASD,NAUXBAS,NAUXSH, &
              NCOR,NACT,NVIR,NORB,NBF,            &
              STOL,VTOL,IVMTD,OTHAUX)
         CALL DDI_DESTROY(D_B)
      ENDIF

      ! IF(LL2.OR.LL3) THEN
      !    IF(LL3) CALL DDI_CREATE(NAUXBAS*NVIR,NACT,D_B)
      !    CALL DDI_SMP_NPROC(SMP_NP,SMP_ME)
      !    CALL DDI_SMP_CREATE(NAUXBAS*NVIR*NACT,D_B_SMP)
      !    CALL DDI_SMP_OFFSET(D_B_SMP,ADDR,LB32)
      !    LB32=LB32+1
      !    CALL RIMP2_X_ENERGY_gpu &
      !        (ADDR(LB32),NANGM,NAUXBASD,NAUXBAS,NAUXSH,                    &
      !         NCOR,NACT,NVIR,NORB,NBF,                                     &
      !         STOL,VTOL,IVMTD,OTHAUX)
      !    CALL DDI_SMP_DESTROY(D_B_SMP)
      !    IF(LL3) CALL DDI_DESTROY(D_B)
      ! ENDIF

call deallocateShellProcessing

#if __PRINT_ESS==1
w1=omp_get_wtime()
IF(BUG) WRITE(*,'(A40,I5,F10.1)') "TIME SUB TOTAL RIMP2_X_ENERGY_gpu",ME,w1-w0
#endif

      CALL DFINAL(1)
      IREST = 0
      
      IF(ISPHER.EQ.1 .OR. ISPHER.EQ.0) DEALLOCATE(LOCSPH)
      DEALLOCATE(DAux,DAtm)


      END
      !***********************************************************





!>*module gpurimp2   *deck rimp2_x_energy_gpu
!>
!>     @brief   rimp2 energy driver on gpus
!>     @author  buu
!>     @date    june 7, 2020
!>
!>     @detail  rimp2 energy driver on gpus
!>

      SUBROUTINE RIMP2_X_ENERGY_gpu &
                (B32,NANGM,NAUXBASD,NAUXBAS,NAUXSH,                     &
                 NCOR,NACT,NVIR,NORB,NBF,                               &
                 STOL,VTOL,IVMTD,OTHAUX)

      use omp_lib

      use gpuGlobalRimp2Eng,only: MXATM
      use gpuGlobalRimp2Eng,only: D_V
      use gpuGlobalRimp2Eng,only: BUG      
      use gpuGlobalRimp2Eng,only: MXAUXSH,MXAXGTOT, MXSH,MXGTOT, MXATM
      use gpuGlobalRimp2Eng,only: LL1,LL2,LL3,LL4,D_B

      implicit none

      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB, &
                      ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      integer :: NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,IAN
      double precision :: ZAN,C
      COMMON /IOFILE/ IR,IW,IP,IJK,IPK,IDAF,NAV,IODA(950)
      integer :: IR,IW,IP,IJK,IPK,IDAF,NAV,IODA(950)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      integer :: ME,MASTER,NPROC,IBTYP,IPTIM
      logical :: DSKWRK, MASWRK, GOPARR
      COMMON /AUXBAS/ EXAUX(MXAXGTOT),CAUXS(MXAXGTOT),CAUXP(MXAXGTOT),  &
                      CAUXD(MXAXGTOT),CAUXF(MXAXGTOT),CAUXG(MXAXGTOT),  &
                      CAUXH(MXAXGTOT),CAUXI(MXAXGTOT),                  &
                      KAUXST(MXAUXSH),KAUXAT(MXAUXSH),KAUXTY(MXAUXSH),  &
                      KAUXNG(MXAUXSH),KAUXLO(MXAUXSH),KAUXMI(MXAUXSH),  &
                      KAUXMX(MXAUXSH),NAUXSHX
      double precision :: &
         EXAUX(MXAXGTOT),CAUXS(MXAXGTOT),CAUXP(MXAXGTOT), &
         CAUXD(MXAXGTOT),CAUXF(MXAXGTOT),CAUXG(MXAXGTOT),CAUXH(MXAXGTOT),CAUXI(MXAXGTOT)
      integer :: &
         KAUXST(MXAUXSH),KAUXAT(MXAUXSH),KAUXTY(MXAUXSH), &
         KAUXNG(MXAUXSH),KAUXLO(MXAUXSH),KAUXMI(MXAUXSH),KAUXMX(MXAUXSH),NAUXSHX
      COMMON /NSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),      &
                      CF(MXGTOT),CG(MXGTOT),CH(MXGTOT),CI(MXGTOT),      &
                      KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),   &
                      KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
      double precision :: &
         EX(MXGTOT),CS(MXGTOT),CP(MXGTOT), &
         CD(MXGTOT),CF(MXGTOT),CG(MXGTOT),CH(MXGTOT),CI(MXGTOT)
      integer :: &
         KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH), &
         KNG(MXSH),KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
      COMMON /SPHERI/ PSHELL(3,3),DSHELL(6,6), &
                      FSHELL(10,10),GSHELL(15,15), &
                      PIHELL(3,3),DIHELL(6,6), &
                      FIHELL(10,10),GIHELL(15,15)
      double precision :: &
         PSHELL(3,3),DSHELL(6,6), &
         FSHELL(10,10),GSHELL(15,15), &
         PIHELL(3,3),DIHELL(6,6), &
         FIHELL(10,10),GIHELL(15,15)

integer :: NANGM,NCOR,NAUXSH,NOCC,NXO,NXV,NXN,NAUXBAS1,MXEX,ISH,MXEXAUX,NAUXBASD,NACT,NAUXSH,NVIR,NORB,NBF,NAUXBAS,IVMTD
double precision :: w0,w1,STOL,VTOL
logical :: OTHAUX



#if __hipblas==1
      double precision,allocatable,target,dimension(:,:) :: VDX, VEC
      double precision,allocatable,dimension(:,:) :: VXX
#else
      double precision,allocatable,dimension(:,:) :: VEC,VDX,VXX !, B32
#endif


      double precision,allocatable :: EIG(:)
      logical :: MASWRK

      ! double precision :: B32(NAUXBAS*NVIR,NACT)
      double precision :: B32(*)
      integer :: free,total
                        !///    INIT    //////



      ! sync ddi procs
      CALL DDI_SYNC(5123)

      ! some parameters
      NOCC=NCOR+NACT
      NORB=NCOR+NACT+NVIR

      NXO=NAUXBASD*NOCC
      NXV=NAUXBASD*NVIR
      NXN=NAUXBASD*NORB

! ! sort ao and aux shells
! ! sort AO shells
!       call SORTAOSHELLS(NSHELL,KTYPE,KNG)
!       write(*,*) "wwwwww",N_S_SHL,N_P_SHL,N_D_SHL



                        !///  V MATRIX  //////

#if __PRINT_TIME_BUG==1
   WRITE(*,*) "debug_debug_debug: RIMP2_VXX_gpu", ME
   CALL RIMP2_MEM_gpu(free,total)
#endif

#if __PRINT_ESS==1   
   w0=omp_get_wtime()
#endif

      ! FORM VDX
      ALLOCATE(VDX(NAUXBAS,NAUXBAS))
      CALL RIMP2_VXX_gpu(VDX, NAUXBAS)

#if __PRINT_ESS==1
   w1 = omp_get_wtime()
   WRITE(*,'(A40,I5,F10.1)') "TIME SUB TOTAL RIMP2_VXX_gpu ",ME, w1-w0
#endif

#if __PRINT_TIME_BUG==1   
   WRITE(*,*) "debug_debug_debug: RIMP2_CHOLESKY_DECOMPOSE_VXX_GPU", ME
   CALL RIMP2_MEM_gpu(free,total)
#endif



#if __PRINT_ESS==1   
   w0=omp_get_wtime()
#endif

      ! CHOLESKY DEC VXX
      NAUXBAS1=NAUXBAS
      CALL RIMP2_CHOLESKY_DECOMPOSE_VXX_GPU (VDX, NAUXBAS1)
      NAUXBASD = NAUXBAS1

#if __PRINT_ESS==1
   w1 = omp_get_wtime()
   WRITE(*,'(A40,I5,F10.1)') "TIME SUB RIMP2_CHOLESKY_DECOMPOSE_VXX ",ME, w1-w0
#endif


#if __PRINT_TIME_BUG==1
   CALL RIMP2_MEM_gpu(free,total)
#endif

 

               !///  SOME INPUT PREPARATION FOR I32  //////


      ! read MOs     
      ALLOCATE(VEC(NBF,NBF))
      CALL DAREAD(IDAF,IODA,VEC,NBF*NBF,15,0)
      CALL RIMP2_SYM_TRANS_MO_gpu(VEC,NBF,NCOR,NORB-NCOR)



      ! MAX NUMBER OF PRIMITIVE GAUSSIANS IN AO BASIS
      MXEX=0
      DO ISH=1,NSHELL
         MXEX=MXEX+KNG(ISH)
      ENDDO

      ! MAX NUMBER OF PRIMITIVE GAUSSIANS IN AUX BASIS
      MXEXAUX=0
      DO ISH=1,NAUXSH
         MXEXAUX=MXEXAUX+KAUXNG(ISH)
      ENDDO


      ! read MO VEC
      ALLOCATE(EIG(NBF))
      CALL DAREAD(IDAF,IODA,EIG,NBF,17,0)

! #if __readRHF==1

!       CALL writeRead &
!           (EIG, VEC,NCOR,NACT,NVIR,NBF, NAUXBAS,                         &
!            EXAUX,CAUXS,CAUXP,CAUXD,CAUXF,CAUXG,CAUXH,CAUXI,             &
!            KAUXST,KAUXAT,KAUXTY,KAUXNG,KAUXLO,KAUXMI,KAUXMX,NAUXSH,     &           
!            EX,CS,CP,CD,CF,CG,CH,CI,                                     &
!            KSTART,KATOM,KTYPE,KNG,KLOC,KMIN,KMAX,NSHELL,                &
!            MXEX,MXEXAUX)
! #endif

      ! IF(MASWRK) WRITE(*,'(A40,2I5)') "NUMBER OF PRIMITIVE AO AND AUX", MXEX,MXEXAUX

                 !///  FORM I32 in MO basis  //////

#if __PRINT_TIME_BUG==1
   WRITE(*,*) "debug_debug_debug: RIMP2_FORM_MO_I32_gpu", ME
   CALL RIMP2_MEM_gpu(free,total)
#endif

#if __PRINT_ESS==1   
   w0=omp_get_wtime()
#endif





      ! FORM I32 IN MO BASIS
      CALL RIMP2_FORM_MO_I32_gpu &
          (B32,VEC,NCOR,NACT,NVIR,NBF, NAUXBAS,                         &
           EXAUX,CAUXS,CAUXP,CAUXD,CAUXF,CAUXG,CAUXH,CAUXI,             &
           KAUXST,KAUXAT,KAUXTY,KAUXNG,KAUXLO,KAUXMI,KAUXMX,NAUXSH,     &           
           EX,CS,CP,CD,CF,CG,CH,CI,                                     &
           KSTART,KATOM,KTYPE,KNG,KLOC,KMIN,KMAX,NSHELL,                &
           PSHELL,DSHELL,FSHELL,GSHELL,PIHELL,DIHELL,FIHELL,GIHELL,     &
           MXEX,MXEXAUX)

      IF(LL2) THEN
         CALL DDI_SMP_SYNC()
      ENDIF



#if __PRINT_ESS==1
   w1=omp_get_wtime()
   WRITE(*,'(A40,I5,F10.1)') "TIME SUB TOTAL RIMP2_FORM_MO_I32_gpu",ME,w1-w0
#endif

#if __PRINT_TIME_BUG==1   
   CALL RIMP2_MEM_gpu(free,total)
#endif

                 !///  COMBINING VDX WITH I32  //////

#if __PRINT_ESS==1
   w0=omp_get_wtime()
#endif

#if __cpu != 1
      IF(LL1.OR.LL2) THEN
         CALL RIMP2_VDX_I32_cublas &
             (B32,EIG(NCOR+1),VDX,NAUXBAS,NAUXBASD, &
              NCOR,NACT,NVIR,NBF)
      ELSEIF(LL3) THEN
         CALL RIMP2_VDX_I32_DDI_cublas &
             (EIG(NCOR+1),VDX,NAUXBAS,NAUXBASD, &
              NCOR,NACT,NVIR,NBF)        
      ENDIF
#endif


#if __PRINT_ESS==1
   w1=omp_get_wtime()
   WRITE(*,'(A40,I5,F10.1)') "TIME SUB TOTAL RIMP2_VDX_I32_cublas",ME,w1-w0
#endif


#if __cpu==1      
      IF(LL1.OR.LL2.OR.LL3) THEN
         CALL RIMP2_VDX_I32_cpu &
             (B32,EIG(NCOR+1),VDX,NAUXBAS,NAUXBASD, &
              NCOR,NACT,NVIR,NBF)
      ENDIF
#endif


! #if __intel_offloading==1
!       IF(LL1.OR.LL2.OR.LL3) THEN
!          CALL RIMP2_VDX_I32_intel &
!              (B32,EIG(NCOR+1),VDX,NAUXBAS,NAUXBASD, &
!               NCOR,NACT,NVIR,NBF)
!       ENDIF
! #endif




      ! IF(LL3) THEN
      !    CALL DDI_SYNC(153)
      !    IF(MOD(ME+1,6).EQ.0) THEN
      !       CALL DDI_GET(D_B,1,NVIR*NAUXBAS,1,NACT,B32)
      !    ENDIF
      ! ENDIF

      ! CALL DDI_SMP_SYNC()

#if __PRINT_TIME_BUG==1
   w1=omp_get_wtime()
   WRITE(*,'(A40,I5,F10.1)') "TIME SUB TOTAL RIMP2_COMBINE_VDX_I32",ME,w1-w0
   ! CALL RIMP2_reset_gpu()
   CALL RIMP2_MEM_gpu(free,total)
#endif
deallocate(VDX)
deallocate(VEC)
deallocate(EIG)
      END
      ! ***************************************************************
 













!*MODULE RIMP2OMP   *DECK worksharedarray_gpu
!>
!>     @BRIEF   WORK SHARING AMONG NPROCS PROCESSES/THREADS
!>     @AUTHOR  BUU PHAM 
!>     @DATE    SEP 25, 2017
!>
!>     @DETAIL  A PORTION OF LOOP FROM Lstart TO Lend IS DIVIDED
!>              AMONG Nprocs (MPI PROCESSES OR THREADS). THE WORK
!>              SHRING IS STORED IN THE ISE ARRAY
!>
      SUBROUTINE WorkSharedArray_gpu(ISE,Nprocs,Lstart,Lend,LOMP)

      implicit none

      INTEGER,intent(in)  :: Nprocs
      INTEGER,intent(in)  :: Lstart, Lend

      LOGICAL,intent(out) :: LOMP
      INTEGER,intent(out) :: ISE(0:Nprocs-1,2)

      INTEGER :: NLOOPS,I,NCHUNK,NMOD

         
         ! ~~~ Start ~~~




      Nloops=Lend-Lstart+1
      IF(Nloops.lt.Nprocs) THEN
         ! number of ranks/threads > loops
         LOMP=.false.
         ISE = 0
         DO i=0, Nloops-1
            ISE(i,1)=Lstart+i
            ISE(i,2)=ISE(i,1)
         ENDDO
         RETURN
      ELSE
         ! number of ranks/threads < loops
         LOMP=.true.
         Nchunk=Nloops/Nprocs
         Nmod=mod(Nloops,Nprocs)

         ISE(0,1)=Lstart
         IF(Nmod.gt.0) then
            ISE(0,2)=Lstart+Nchunk
         ELSE
            ISE(0,2)=Lstart+Nchunk-1
         ENDIF
         Nmod=Nmod-1

         DO i=1, Nprocs-1
            ISE(i,1)=ISE(i-1,2)+1
            IF(Nmod.gt.0) then
               ISE(i,2)=ISE(i-1,2)+Nchunk+1
               Nmod=Nmod-1
            ELSE
               ISE(i,2)=ISE(i-1,2)+Nchunk
            ENDIF
         ENDDO
      ENDIF


      END
      !**********************************************************










      

!       SUBROUTINE RIMP2_AUXLOC_gpu(IauxLOC, KAUXSH,LSOE)
      
!       use gpuGlobalRimp2Eng,only: MXAXGTOT,MXAUXSH
!       use gpuGlobalRimp2Eng,only: MaxAUXBFX,MaxAUXANG
!       use gpuGlobalRimp2Eng,only: MaxATMANG,MaxATMBFX
!       use gpuGlobalRimp2Eng,only: FlgSphAux,LOCSPH



!       COMMON /AUXBAS/ EXAUX(MXAXGTOT),CAUXS(MXAXGTOT),CAUXP(MXAXGTOT),  &
!                       CAUXD(MXAXGTOT),CAUXF(MXAXGTOT),CAUXG(MXAXGTOT),  &
!                       CAUXH(MXAXGTOT),CAUXI(MXAXGTOT),                  &
!                       KAUXST(MXAUXSH),KAUXAT(MXAUXSH),KAUXTY(MXAUXSH),  &
!                       KAUXNG(MXAUXSH),KAUXLO(MXAUXSH),KAUXMI(MXAUXSH),  &
!                       KAUXMX(MXAUXSH),NAUXSH

!       double precision,allocatable :: TL(:)

!       logical :: LSOE


!       MXBSH=MaxAUXBFX
!       IF(MaxAUXANG < MaxATMANG) MXBSH=MaxATMBFX

!       LANGK = KAUXTY(KAUXSH)-1
!       MINK  = KAUXMI(KAUXSH)
!       MAXK  = KAUXMX(KAUXSH)
!       NK    = MAXK-MINK+1

!       ALLOCATE(TL(MXBSH*MXBSH))
!       IF(FlgSphAux) THEN
!             CALL RIMP2CSTRM(TL,LANGK,MINK,NK,MINKS,MAXKS)
!             LOCKS = LOCSPH(KAUXSH)-MINKS
!             NKS = MAXKS-MINKS+1
!             MINKT = MINKS
!             MAXKT = MAXKS
!             LOCKT = LOCKS
!       ELSE
!             MINKT = MINK
!             MAXKT = MAXK
!             LOCKT = LOCK
!       ENDIF

!       NKT=MAXKT-MINKT+1

!       IauxStart=LOCKT+MINKT
!       IauxEnd=IauxStart+NKT-1

!       IF(LSOE) THEN
!          IauxLOC=IauxStart
!       ELSE
!          IauxLOC=IauxEnd
!       ENDIF
! deallocate(TL)
!       END
!       ! *************************************************************









#if __cublas==1 || __hipblas==1
!>*module gpurimp2   *deck RIMP2_VDX_I32_gpu
!>
!>     @brief   combine VDX with I32
!>     @author  buu
!>     @date    june 7, 2020
!>
!>     @detail  combine VDX with I32
!>

      SUBROUTINE RIMP2_VDX_I32_cublas &
                (B32,EIG,VDX,NAUXBAS,NAUXBASD, &
                 NCOR,NACT,NVIR,NBF)

      use omp_lib
      use blasHandles
      use iso_c_binding
      use gpuGlobalRimp2Eng
      
      implicit none

      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      integer :: ME,MASTER,NPROC,IBTYP,IPTIM
      logical :: DSKWRK, MASWRK, GOPARR
      COMMON /IOFILE/ IR,IW,IP,IJK,IPK,IDAF,NAV,IODA(950)
      integer :: IR,IW,IP,IJK,IPK,IDAF,NAV,IODA(950)
      COMMON /ENRGYS/ ENUC,EELEC,ETOT,SZ,SZZ,ECORE,ESCF,EERD,E1,E2,     &
                      VEN,VEE,EPOT,EKIN,ESTATE(MXRT),STATN,EDFT(2),EDISP
      double precision :: ENUC,EELEC,ETOT,SZ,SZZ,ECORE,ESCF,EERD,E1,E2,     &
                      VEN,VEE,EPOT,EKIN,ESTATE(MXRT),STATN,EDFT(2),EDISP
      COMMON /MACHIN/ NWDVAR,MAXFM,MAXSM,LIMFM,LIMSM
      integer :: NWDVAR,MAXFM,MAXSM,LIMFM,LIMSM
      COMMON /ENRGMP/ EMP2,EMP3,EMP4,EMP2A
      double precision :: EMP2,EMP3,EMP4,EMP2A

      ! INPUT
      ! double precision :: B32(NAUXBAS*NVIR,NACT)
      double precision :: B32(NAUXBAS*NVIR,*)
      double precision :: VDX(NAUXBASD,NAUXBAS)
      double precision :: EIG(NACT+NVIR)

      ! LOCAL      
      double precision,allocatable :: QVV(:),T3(:),BI(:)
      integer,allocatable :: ISE(:,:)


      integer :: Nchunks
      logical :: GOPARR,DSKWRK,MASWRK

      integer :: SMP_NP,SMP_ME

      integer(c_int),dimension(1) :: cublasXt_deviceId
      integer(c_int) :: ndev,blockdim
      integer :: free,total
      
      
integer :: NAUXBAS,NAUXBASD,NCOR,NACT,NVIR,NBF,NOCC
integer :: LMPI_S,LMPI_E,MPI_S_E,MPI_S_E
integer :: NQVV,NSLICE,ISLICE, LGPU_S,LGPU_E
integer :: JACT,iQVV,IACT,II
double precision :: wddi0,wddi1
logical :: LOMP

      ! INIT CORRELATION ENERGY
      E2=0.0D00

      ! SMP
      CALL DDI_SMP_NPROC(SMP_NP,SMP_ME)

      ! MO ENERGY
      NOCC=NCOR+NACT

      ! GPU MEMORY
      CALL RIMP2_MEM_gpu(free,total)

      IF(LL3) ALLOCATE(BI(NVIR*NAUXBASD))


      ! MPI DISTRIBUTION OF OCC-OCC
      CALL RIMP2_TRAPE_DEC(LMPI_S,LMPI_E,ME,NPROC,NACT)
      MPI_S_E=LMPI_E-LMPI_S+1

      ! check if GPU memory if enough to hold an "easy" slice of B32
      IF(free .LT. NVIR*NAUXBAS*MPI_S_E*8) THEN
         WRITE(*,*) "***GPU MEMORY IS NOT ENOUGH TO HOLD A CPU SLICE OF B32***"
      ENDIF

         CALL DDI_SYNC(1562)

      ! set QVV size
      NQVV=MPI_S_E
      NSLICE=1
      DO II=1,MPI_S_E
         IF(NVIR*NVIR*NQVV*8 + NVIR*NAUXBAS*NQVV*8 + 2*NVIR*NAUXBAS*8 .LT. 0.8*free) EXIT
         NSLICE = NSLICE + 1
         NQVV = MPI_S_E/NSLICE + 1
      ENDDO
      IF(NQVV.LT.1) THEN
         write(*,*) "***not enough GPU memory for a single slice***"
         CALL ABRT
      ENDIF

      ! write(*,'(A30,5I10)') "NQVV,NSLICE,MPI_S_E,NACT",ME, NQVV, NSLICE,MPI_S_E,NACT


      ! allocate buffer arrays
      ALLOCATE(T3(NAUXBAS*NVIR))
      ALLOCATE(QVV(NVIR*NQVV*NVIR))


      ALLOCATE(ISE(0:NSLICE-1,2))
      CALL WorkSharedArray_gpu(ISE,NSLICE,LMPI_S,LMPI_E,LOMP)

#if __cublas==1
      ndev = 1
      blockdim = 2048
      cublasXt_deviceId(1) = 0
      cublas_return = cublasXtcreate(cublas_handle)
      cublas_return = cublasXtDeviceSelect(cublas_handle, ndev, cublasXt_deviceId)
      cublas_return = cublasXtSetBlockDim(cublas_handle, blockdim)
#endif
#if __hipblas==1
      call hipblasCheck(hipblasCreate(hipblas_handle))
#endif

#if __PRINT_TIME_BUG==1
   wI=0.0D00
   wII=0.0D00
   w0=omp_get_wtime()
#endif



!$omp target data map(to:EIG,VDX) map(alloc:QVV,T3)

      DO ISLICE=0,NSLICE-1     
         LGPU_S = ISE(ISLICE,1)
         LGPU_E = ISE(ISLICE,2)



wddi0 = omp_get_wtime()
      ! GET B32 if it is stored in distributed memory
      IF(LL3) THEN
         CALL DDI_GET(D_B,1,NVIR*NAUXBAS,LGPU_S,LGPU_E,B32(1:NVIR*NAUXBAS,LGPU_S:LGPU_E))
      ENDIF
wddi1 = omp_get_wtime()
#if __PRINT_ESS==1
write(*,'(A20,3I10,F15.3)') "DDI_GET time", ME, LGPU_S,LGPU_E, wddi1-wddi0
#endif

!$omp target data map(to:B32(1:NVIR*NAUXBAS,LGPU_S:LGPU_E))




#if __PRINT_TIME_BUG==1
   w1 = omp_get_wtime()
#endif

         DO JACT=LGPU_S,LGPU_E
            CALL RIMP2_B32_TRANSPOSE(B32(1,JACT),T3,NAUXBAS,NVIR)
#if __cublas==1 || __hipblas==1
            CALL RIMP2_FORM_B32_one(B32(1,JACT),T3,VDX,NAUXBASD,NVIR,NAUXBAS)
#elif __intel_offloading==1
            CALL RIMP2_FORM_B32_one_intel(B32(1,JACT),T3,VDX,NAUXBASD,NVIR,NAUXBAS)
#endif
            iQVV = JACT-LGPU_S+1
#if __cublas==1 || __hipblas==1
            CALL RIMP2_QVV_HEAD_cublas &
                (QVV,B32(1,LGPU_S),B32(1,JACT), &
                 NACT,NVIR,NAUXBASD,iQVV)
#elif __intel_offloading==1
            CALL RIMP2_QVV_HEAD_mkl &
                (QVV,B32(1,LGPU_S),B32(1,JACT), &
                 NACT,NVIR,NAUXBASD,iQVV)              
#endif
            CALL RIMP2_ENG_HEAD(E2,EIG,QVV,LGPU_S,JACT,NACT,NVIR,iQVV)

         ENDDO !JACT

#if __PRINT_TIME_BUG==1
   w2 = omp_get_wtime()
   wI = wI + (w2-w1)
#endif

         DO IACT=1,LGPU_S-1

            ! GET B32 if it is stored in distributed memory
            IF(LL3) THEN
               CALL DDI_GET(D_B,1,NVIR*NAUXBAS,IACT,IACT,BI)
               CALL RIMP2_TAIL &
                   (E2,BI,B32(1,LGPU_S),VDX,EIG,QVV,T3, &
                    NACT,NVIR,IACT,LGPU_S,LGPU_E,NQVV,NAUXBASD,NAUXBAS) 
            ELSE
               CALL RIMP2_TAIL &
                   (E2,B32(1,IACT),B32(1,LGPU_S),VDX,EIG,QVV,T3, &
                    NACT,NVIR,IACT,LGPU_S,LGPU_E,NQVV,NAUXBASD,NAUXBAS) 
            ENDIF


         ENDDO !IACT

#if __PRINT_TIME_BUG==1
   w3 = omp_get_wtime()
   wII = wII + (w3-w2)
#endif


!$omp end target data
      ENDDO
!$omp end target data

#if __cublas==1
cublas_return = cublasXtdestroy(cublas_handle)
#endif
#if __hipblas==1
call hipblasCheck(hipblasDestroy(hipblas_handle))
#endif


#if __PRINT_TIME_BUG==1
   w4=omp_get_wtime()
   write(*,'(A40,3F10.1)') "CHECK BALANCING:TOTAL:HEAD:TAIL", w4-w0, wI, wII
#endif


      ! correlation energy reduction
      CALL DDI_GSUMF(5999,E2,1)

      ! rhf + E(2)
      EMP2 = E2+ESCF


      ! print out result
      IF(MASWRK) WRITE(IW,100) ESCF,E2,EMP2


  100 FORMAT (/1X,'RHF RI MP2 ENERGY',/,                                &
             12X,'   ESCF=',1X,F20.10/                                  &
             12X,'   E(2)=',1X,F20.10/                                  &
             12X,' E(MP2)=',1X,F20.10 ,/)


! DEALLOCATE(ISE)
DEALLOCATE(ISE)
DEALLOCATE(T3)
if (allocated(BI)) then
   deallocate(BI)
endif
deallocate(QVV)
      END
#endif

#if __cublas==1 || __hipblas==1
!>*module gpurimp2   *deck RIMP2_VDX_I32_gpu
!>
!>     @brief   combine VDX with I32
!>     @author  buu
!>     @date    june 7, 2020
!>
!>     @detail  combine VDX with I32
!>

      SUBROUTINE RIMP2_VDX_I32_DDI_cublas &
                (EIG,VDX,NAUXBAS,NAUXBASD, &
                 NCOR,NACT,NVIR,NBF)

      
      use omp_lib
      use blasHandles
      use iso_c_binding
      use gpuGlobalRimp2Eng
      
      implicit none

      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      integer :: ME,MASTER,NPROC,IBTYP,IPTIM
      logical :: DSKWRK, MASWRK, GOPARR
      COMMON /IOFILE/ IR,IW,IP,IJK,IPK,IDAF,NAV,IODA(950)
      integer :: IR,IW,IP,IJK,IPK,IDAF,NAV,IODA(950)
      COMMON /ENRGYS/ ENUC,EELEC,ETOT,SZ,SZZ,ECORE,ESCF,EERD,E1,E2,     &
                      VEN,VEE,EPOT,EKIN,ESTATE(MXRT),STATN,EDFT(2),EDISP
      double precision :: ENUC,EELEC,ETOT,SZ,SZZ,ECORE,ESCF,EERD,E1,E2,     &
                      VEN,VEE,EPOT,EKIN,ESTATE(MXRT),STATN,EDFT(2),EDISP
      COMMON /ENRGMP/ EMP2,EMP3,EMP4,EMP2A
      double precision :: EMP2,EMP3,EMP4,EMP2A
      COMMON /MACHIN/ NWDVAR,MAXFM,MAXSM,LIMFM,LIMSM
      integer :: NWDVAR,MAXFM,MAXSM,LIMFM,LIMSM

integer :: NAUXBAS,NAUXBASD,NCOR,NACT,NVIR,NBF,NOCC
integer :: LMPI_S,LMPI_E,LGPU_S,LGPU_E,MPI_S_E,NQVV
integer :: NSLICE,ISLICE
integer :: JACT,IACT,JACT_new,iQVV
integer :: II,IMAX,NN
logical :: LOMP

      ! INPUT
      double precision :: VDX(NAUXBASD,NAUXBAS)
      double precision :: EIG(NACT+NVIR)

      ! LOCAL      
      double precision,allocatable :: QVV(:),T3(:),BI(:),B32(:,:)
      integer,allocatable :: ISE(:,:)


      integer :: Nchunks
      logical :: GOPARR,DSKWRK,MASWRK

      integer :: SMP_NP,SMP_ME

      integer(c_int),dimension(1) :: cublasXt_deviceId
      integer(c_int) :: ndev,blockdim
      integer :: free,total
      
      

      ! INIT CORRELATION ENERGY
      E2=0.0D00

      ! SMP
      CALL DDI_SMP_NPROC(SMP_NP,SMP_ME)

      ! MO ENERGY
      NOCC=NCOR+NACT

      ! GPU MEMORY
      CALL RIMP2_MEM_gpu(free,total)

      ALLOCATE(BI(NVIR*NAUXBASD))


      ! MPI DISTRIBUTION OF OCC-OCC
      CALL RIMP2_TRAPE_DEC(LMPI_S,LMPI_E,ME,NPROC,NACT)
      MPI_S_E=LMPI_E-LMPI_S+1

! write(*,'(A20,4I10)') "ME,LMPI_S,LMPI_E,NACT", ME,LMPI_S,LMPI_E,NACT

      ! check if GPU memory if enough to hold an "easy" slice of B32
      IF(free .LT. NVIR*NAUXBAS*MPI_S_E*8) THEN
         WRITE(*,*) "***GPU MEMORY IS NOT ENOUGH TO HOLD A CPU SLICE OF B32***"
      ENDIF

         CALL DDI_SYNC(1562)

IF(ME.GT.NACT) GOTO 911

      ! set QVV size
      NQVV=MPI_S_E
      NSLICE=1
      DO II=1,MPI_S_E
         IF(NVIR*NVIR*NQVV*8 + NVIR*NAUXBAS*NQVV*8 + 2*NVIR*NAUXBAS*8 .LT. 0.8*free) EXIT
         NSLICE = NSLICE + 1
         NQVV = MPI_S_E/NSLICE + 1
      ENDDO
      IF(NQVV.LT.1) THEN
         write(*,*) "not enough GPU memory"
      ENDIF

      ! write(*,'(A30,5I10)') "NQVV,NSLICE,MPI_S_E,NACT",ME, NQVV, NSLICE,MPI_S_E,NACT


      ! allocate buffer arrays
      ALLOCATE(T3(NAUXBAS*NVIR))
      ALLOCATE(QVV(NVIR*NQVV*NVIR))


      ALLOCATE(ISE(0:NSLICE-1,2))
      CALL WorkSharedArray_gpu(ISE,NSLICE,LMPI_S,LMPI_E,LOMP)

      IMAX = 1
      DO ISLICE=0,NSLICE-1     
         LGPU_S = ISE(ISLICE,1)
         LGPU_E = ISE(ISLICE,2)
         NN = LGPU_E - LGPU_S + 1
         IF(NN.GT.IMAX) IMAX = NN
      ENDDO

! write(*,*) "IMAX: ", IMAX

      ALLOCATE(B32(1:NVIR*NAUXBAS,IMAX))
! write(*,*) "afterwww"
#if __cublas==1
      ndev = 1
      blockdim = 2048
      cublasXt_deviceId(1) = 0
      cublas_return = cublasXtcreate(cublas_handle)
      cublas_return = cublasXtDeviceSelect(cublas_handle, ndev, cublasXt_deviceId)
      cublas_return = cublasXtSetBlockDim(cublas_handle, blockdim)
#endif
#if __hipblas==1
      call hipblasCheck(hipblasCreate(hipblas_handle))
#endif

#if __PRINT_TIME_BUG==1
   wI=0.0D00
   wII=0.0D00
   w0=omp_get_wtime()
#endif



!$omp target data map(to:EIG,VDX) map(alloc:QVV,T3)


      DO ISLICE=0,NSLICE-1     
         LGPU_S = ISE(ISLICE,1)
         LGPU_E = ISE(ISLICE,2)


! wddi0 = omp_get_wtime()
         CALL DDI_GET(D_B,1,NVIR*NAUXBAS,LGPU_S,LGPU_E,B32(1:NVIR*NAUXBAS,1:LGPU_E-LGPU_S+1))
! wddi1 = omp_get_wtime()
! write(*,'(A20,3I10,F15.3)') "DDI_GET time", ME, LGPU_S,LGPU_E, wddi1-wddi0


!$omp target data map(to:B32)




#if __PRINT_TIME_BUG==1
   w1 = omp_get_wtime()
#endif

         DO JACT=LGPU_S,LGPU_E
            JACT_new = JACT - LGPU_S + 1
            CALL RIMP2_B32_TRANSPOSE(B32(1,JACT_new),T3,NAUXBAS,NVIR)
            CALL RIMP2_FORM_B32_one(B32(1,JACT_new),T3,VDX,NAUXBASD,NVIR,NAUXBAS)
            iQVV = JACT-LGPU_S+1
            CALL RIMP2_QVV_HEAD_cublas &
                (QVV,B32,B32(1,JACT_new), &
                 NACT,NVIR,NAUXBASD,iQVV)
            CALL RIMP2_ENG_HEAD(E2,EIG,QVV,LGPU_S,JACT,NACT,NVIR,iQVV)

         ENDDO !JACT

#if __PRINT_TIME_BUG==1
   w2 = omp_get_wtime()
   wI = wI + (w2-w1)
#endif

         DO IACT=1,LGPU_S-1

            ! GET B32 if it is stored in distributed memory

            CALL DDI_GET(D_B,1,NVIR*NAUXBAS,IACT,IACT,BI)
            CALL RIMP2_TAIL &
                (E2,BI,B32,VDX,EIG,QVV,T3, &
                 NACT,NVIR,IACT,LGPU_S,LGPU_E,NQVV,NAUXBASD,NAUXBAS) 


         ENDDO !IACT

#if __PRINT_TIME_BUG==1
   w3 = omp_get_wtime()
   wII = wII + (w3-w2)
#endif


!$omp end target data
      ENDDO
!$omp end target data
#if __cublas==1
cublas_return = cublasXtdestroy(cublas_handle)
#endif
#if __hipblas==1
call hipblasCheck(hipblasDestroy(hipblas_handle))
#endif

#if __PRINT_TIME_BUG==1
   w4=omp_get_wtime()
   write(*,'(A40,3F10.1)') "CHECK BALANCING:TOTAL:HEAD:TAIL", w4-w0, wI, wII
#endif


      ! correlation energy reduction
      CALL DDI_GSUMF(5999,E2,1)

      ! rhf + E(2)
      EMP2 = E2+ESCF


      ! print out result
      IF(MASWRK) WRITE(IW,100) ESCF,E2,EMP2


  100 FORMAT (/1X,'RHF RI MP2 ENERGY',/,                                &
             12X,'   ESCF=',1X,F20.10/                                  &
             12X,'   E(2)=',1X,F20.10/                                  &
             12X,' E(MP2)=',1X,F20.10 ,/)


! DEALLOCATE(ISE)
DEALLOCATE(ISE)
DEALLOCATE(T3)
deallocate(BI)
deallocate(QVV)
deallocate(B32)
   911 CONTINUE
if (allocated(ISE)) then
   deallocate(ISE)
endif
if (allocated(T3)) then
   deallocate(T3)
endif
if (allocated(BI)) then
   deallocate(BI)
endif
if (allocated(QVV)) then
   deallocate(QVV)
endif
if (allocated(B32)) then
   deallocate(B32)
endif
      END


#endif




#if __cpu==1
!>*module gpurimp2   *deck RIMP2_VDX_I32_gpu
!>
!>     @brief   combine VDX with I32
!>     @author  buu
!>     @date    june 7, 2020
!>
!>     @detail  combine VDX with I32
!>

      SUBROUTINE RIMP2_VDX_I32_cpu &
                (B32,EIG,VDX,NAUXBAS,NAUXBASD, &
                 NCOR,NACT,NVIR,NBF)

      
      use omp_lib
      use cublasf
      use cublasHandles
      use iso_c_binding
      use gpuGlobalRimp2Eng
      
      implicit none

      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      integer :: ME,MASTER,NPROC,IBTYP,IPTIM
      logical :: DSKWRK, MASWRK, GOPARR
      COMMON /IOFILE/ IR,IW,IP,IJK,IPK,IDAF,NAV,IODA(950)
      integer :: IR,IW,IP,IJK,IPK,IDAF,NAV,IODA(950)
      COMMON /ENRGYS/ ENUC,EELEC,ETOT,SZ,SZZ,ECORE,ESCF,EERD,E1,E2,     &
                      VEN,VEE,EPOT,EKIN,ESTATE(MXRT),STATN,EDFT(2),EDISP
      double precision :: ENUC,EELEC,ETOT,SZ,SZZ,ECORE,ESCF,EERD,E1,E2,     &
                      VEN,VEE,EPOT,EKIN,ESTATE(MXRT),STATN,EDFT(2),EDISP
      COMMON /ENRGMP/ EMP2,EMP3,EMP4,EMP2A
      double precision :: EMP2,EMP3,EMP4,EMP2A
      COMMON /MACHIN/ NWDVAR,MAXFM,MAXSM,LIMFM,LIMSM
      integer :: NWDVAR,MAXFM,MAXSM,LIMFM,LIMSM

      ! INPUT
      double precision :: B32(NAUXBAS*NVIR,NACT)
      double precision :: VDX(NAUXBASD,NAUXBAS)
      double precision :: EIG(NACT+NVIR)

      ! LOCAL      
      double precision,allocatable :: QVV(:),T3(:),B32_GPU(:,:),BIJ(:)
      integer,allocatable :: ISE(:,:)


      integer :: Nchunks
      logical :: GOPARR,DSKWRK,MASWRK

      integer :: SMP_NP,SMP_ME

      double precision,allocatable,dimension(:,:) :: eij,eab
      

      ! INIT CORRELATION ENERGY
      E2=0.0D00

      ! SMP
      CALL DDI_SMP_NPROC(SMP_NP,SMP_ME)

      ! MO ENERGY
      NOCC=NCOR+NACT

      ! MPI DISTRIBUTION OF OCC-OCC
      CALL RIMP2_TRAPE_DEC(LMPI_S,LMPI_E,ME,NPROC,NACT)
      MPI_S_E=LMPI_E-LMPI_S+1

      ! GET B32
      IF(LL3) THEN
         CALL DDI_SYNC(1562)
         KK=CEILING(NVIR*NAUXBAS*NACT*8.0/2147483643.0)
         LL=CEILING(REAL(NACT)/REAL(KK))
         DO II=1,KK
            ISTART=(II-1)*LL+1
            IEND = II*LL
            IF(IEND.GT.NACT) IEND=NACT
            CALL DDI_GET(D_B,1,NVIR*NAUXBAS,ISTART,IEND,B32(1:NVIR*NAUXBAS,ISTART:IEND))
         ENDDO
      ENDIF

      ! BUFFERS
      ALLOCATE(T3(NAUXBAS*NVIR))
      NQVV=NACT/2
      ALLOCATE(QVV(NVIR*NQVV*NVIR))

      ! GPU BATCHING
free = 16914055168      
      NSLICE = CEILING(real(SIZEOF(B32(1:NVIR*NAUXBAS,LMPI_S:LMPI_E)))/(real(free)*0.8))
      ALLOCATE(ISE(0:NSLICE-1,2))
      CALL WorkSharedArray_gpu(ISE,NSLICE,LMPI_S,LMPI_E,LOMP)


MAXSLICE = 0
DO ISLICE=0,NSLICE-1
   LS=ISE(ISLICE,1)
   LE=ISE(ISLICE,2)
   IF(LE-LS+1.GT.MAXSLICE) MAXSLICE=LE-LS+1
ENDDO
ALLOCATE(B32_GPU(NVIR*NAUXBAS,MAXSLICE))
ALLOCATE(BIJ(NAUXBAS*NVIR))


      DO ISLICE=0,NSLICE-1     
         LGPU_S = ISE(ISLICE,1)
         LGPU_E = ISE(ISLICE,2)
         B32_GPU(1:NVIR*NAUXBAS,1:LGPU_E-LGPU_S+1) = B32(1:NVIR*NAUXBAS,LGPU_S:LGPU_E)
         DO JACT=LGPU_S,LGPU_E
            CALL RIMP2_B32_TRANSPOSE_cpu(B32(1,JACT),T3,NAUXBAS,NVIR)
            CALL RIMP2_FORM_B32_one_cpu &
                (B32_GPU(1,JACT-LGPU_S+1),T3, &
                 VDX,NAUXBASD,NVIR,NAUXBAS)
            DO IACTmod=1,(JACT-LGPU_S)/NQVV+1
               IACT = LGPU_S + (IACTmod-1)*NQVV
               IF(IACTmod*NQVV>JACT-LGPU_S+1) THEN
                  iQVV = JACT-LGPU_S+1-(IACTmod-1)*NQVV
               ELSE
                  iQVV = NQVV
               ENDIF
               CALL RIMP2_ENERGYIJ_cpu &
                   (E2, B32_GPU(1,IACT-LGPU_S+1),B32_GPU(1,JACT-LGPU_S+1), &
                    eij,eab,EIG, QVV,IACT,JACT,NACT, &
                    NVIR,NAUXBASD,iQVV)

            ENDDO
         ENDDO !JACT




         DO IACT=1,LGPU_S-1
            CALL RIMP2_B32_TRANSPOSE_cpu(B32(1,IACT),T3,NAUXBAS,NVIR)
            CALL RIMP2_FORM_B32_one_cpu(BIJ,T3,VDX,NAUXBASD,NVIR,NAUXBAS)
            DO JACTmod = 1,(LGPU_E-LGPU_S)/NQVV+1
               JACT = LGPU_S + (JACTmod-1)*NQVV
               IF(JACTmod*NQVV.GT.LGPU_E-LGPU_S+1) THEN
                  iQVV = LGPU_E-JACT+1-(JACTmod-1)*NQVV
               ELSE
                  iQVV = NQVV
               ENDIF
               CALL RIMP2_ENERGYIJ_II_cpu &
                   (E2, BIJ,B32_GPU(1,JACT-LGPU_S+1), &
                    eij,eab,EIG, QVV,IACT,JACT,NACT, &
                    NVIR,NAUXBASD,iQVV)
            ENDDO
         ENDDO
      ENDDO




      ! correlation energy reduction
      CALL DDI_GSUMF(5999,E2,1)

      ! rhf + E(2)
      EMP2 = E2+ESCF

      ! print out result
      IF(MASWRK) WRITE(IW,100) ESCF,E2,EMP2
  100 FORMAT (/1X,'RHF RI MP2 ENERGY',/,          &
             12X,'   ESCF=',1X,F20.10/            &
             12X,'   E(2)=',1X,F20.10/            &
             12X,' E(MP2)=',1X,F20.10 ,/)


! DEALLOCATE(ISE)
DEALLOCATE(ISE)
DEALLOCATE(T3)
deallocate(B32_GPU)
deallocate(BIJ)
      END
#endif









! !>*module gpurimp2   *deck RIMP2_VDX_I32_gpu
! !>
! !>     @brief   combine VDX with I32
! !>     @author  buu
! !>     @date    june 7, 2020
! !>
! !>     @detail  combine VDX with I32
! !>

!       SUBROUTINE RIMP2_VDX_I32_intel &
!                 (B32,EIG,VDX,NAUXBAS,NAUXBASD, &
!                  NCOR,NACT,NVIR,NBF)

      
!       use omp_lib
!       use cublasf
!       use cublasHandles
!       use iso_c_binding
!       use gpuGlobalRimp2Eng
      

!       COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
!       COMMON /IOFILE/ IR,IW,IP,IJK,IPK,IDAF,NAV,IODA(950)
!       COMMON /ENRGYS/ ENUCR,EELCT,ETOT,SZ,SZZ,ECORE,ESCF,EERD,E1,E2,    &
!                       VEN,VEE,EPOT,EKIN,ESTATE(MXRT),STATN,EDFT(2),EDISP
!       COMMON /ENRGMP/ EMP2,EMP3,EMP4,EMP2A
!       COMMON /MACHIN/ NWDVAR,MAXFM,MAXSM,LIMFM,LIMSM      

!       ! INPUT
!       double precision :: B32(NAUXBAS*NVIR,NACT)
!       double precision :: VDX(NAUXBASD,NAUXBAS)
!       double precision :: EIG(NACT+NVIR)

!       ! LOCAL      
!       double precision,allocatable :: QVV(:),T3(:)
!       integer,allocatable :: ISE(:,:)


!       integer :: Nchunks
!       logical :: GOPARR,DSKWRK,MASWRK

!       integer :: SMP_NP,SMP_ME

!       double precision,allocatable,dimension(:,:) :: eij,eab
      

!       ! INIT CORRELATION ENERGY
!       E2=0.0D00

!       ! SMP
!       CALL DDI_SMP_NPROC(SMP_NP,SMP_ME)

!       ! MO ENERGY
!       NOCC=NCOR+NACT

!       ! MPI DISTRIBUTION OF OCC-OCC
!       CALL RIMP2_TRAPE_DEC(LMPI_S,LMPI_E,ME,NPROC,NACT)
!       MPI_S_E=LMPI_E-LMPI_S+1

!       ! GET B32
!       IF(LL3) THEN
!          CALL DDI_SYNC(1562)
!          KK=CEILING(NVIR*NAUXBAS*NACT*8.0/2147483643.0)
!          LL=CEILING(REAL(NACT)/REAL(KK))
!          DO II=1,KK
!             ISTART=(II-1)*LL+1
!             IEND = II*LL
!             IF(IEND.GT.NACT) IEND=NACT
!             CALL DDI_GET(D_B,1,NVIR*NAUXBAS,ISTART,IEND,B32(1:NVIR*NAUXBAS,ISTART:IEND))
!          ENDDO
!       ENDIF

! #if __notes==1
!    write(*,*) "NOTES: make QVV chunk more flexible: NQVV=NACT/2"
! #endif
!       ! BUFFERS
!       ALLOCATE(T3(NAUXBAS*NVIR))
!       NQVV=NACT/2
!       ALLOCATE(QVV(NVIR*NQVV*NVIR))

!       ! GPU BATCHING
! #if __notes==1
!    write(*,*) "NOTES: make GPU more flexible: free = 16914055168"
! #endif
!       ! free = 16914055168
!       free = 6914055168
!       NSLICE = CEILING(real(SIZEOF(B32(1:NVIR*NAUXBAS,LMPI_S:LMPI_E)))/(real(free)*0.8))
!       ALLOCATE(ISE(0:NSLICE-1,2))
!       CALL WorkSharedArray_gpu(ISE,NSLICE,LMPI_S,LMPI_E,LOMP)



! !$omp target data map(to:EIG,VDX) map(alloc:QVV,T3)
!       DO ISLICE=0,NSLICE-1     
!          LGPU_S = ISE(ISLICE,1)
!          LGPU_E = ISE(ISLICE,2)
! !$omp target data map(to:B32(1:NVIR*NAUXBAS,LGPU_S:LGPU_E))
!          ! HEAD
!          DO JACT=LGPU_S,LGPU_E
!             CALL RIMP2_B32_TRANSPOSE(B32(1,JACT),T3,NAUXBAS,NVIR)
!             CALL RIMP2_FORM_B32_one_intel(B32(1,JACT),T3,VDX,NAUXBASD,NVIR,NAUXBAS)
!             DO IACTmod=1,(JACT-LGPU_S)/NQVV+1
!                IACT = LGPU_S + (IACTmod-1)*NQVV
!                IF(IACTmod*NQVV>JACT-LGPU_S+1) THEN
!                   iQVV = JACT-LGPU_S+1-(IACTmod-1)*NQVV
!                ELSE
!                   iQVV = NQVV
!                ENDIF
!                CALL RIMP2_ENERGYIJ_intel &
!                    (E2, B32(1,IACT),B32(1,JACT), &                   
!                     eij,eab,EIG, QVV,IACT,JACT,NACT, &
!                     NVIR,NAUXBASD,iQVV)
!             ENDDO
!          ENDDO !JACT
!          ! TAIL
!          DO IACT=1,LGPU_S-1
! !$omp target data map(to:B32(1:NVIR*NAUXBAS,IACT:IACT)) 
!             CALL RIMP2_B32_TRANSPOSE(B32(1,IACT),T3,NAUXBAS,NVIR)
!             CALL RIMP2_FORM_B32_one_intel(B32(1,IACT),T3,VDX,NAUXBASD,NVIR,NAUXBAS)
!             DO JACTmod = 1,(LGPU_E-LGPU_S)/NQVV+1
!                JACT = LGPU_S + (JACTmod-1)*NQVV
!                IF(JACTmod*NQVV.GT.LGPU_E-LGPU_S+1) THEN
!                   iQVV = LGPU_E-JACT+1-(JACTmod-1)*NQVV
!                ELSE
!                   iQVV = NQVV
!                ENDIF
!                CALL RIMP2_ENERGYIJ_II_intel &
!                    (E2, B32(1,IACT),B32(1,JACT), &
!                     eij,eab,EIG, QVV,IACT,JACT,NACT, &
!                     NVIR,NAUXBASD,iQVV)
!             ENDDO !JACTmod
! !$omp end target data 
!          ENDDO !IACT
! !$omp end target data         
!       ENDDO ! ISLICE
! !$omp end target data


!       ! correlation energy reduction
!       CALL DDI_GSUMF(5999,E2,1)

!       ! rhf + E(2)
!       EMP2 = E2+ESCF

!       ! print out result
!       IF(MASWRK) WRITE(IW,100) ESCF,E2,EMP2
!   100 FORMAT (/1X,'RHF RI MP2 ENERGY',/,                                &
!              12X,'   ESCF=',1X,F20.10/                                  &
!              12X,'   E(2)=',1X,F20.10/                                  &
!              12X,' E(MP2)=',1X,F20.10 ,/)


! ! DEALLOCATE(ISE)
! DEALLOCATE(ISE)
! DEALLOCATE(T3)

!       END
!       !***********************************************************







      SUBROUTINE RIMP2_FORM_B32 &
          (B32,T3,VDX, &
           NAUXBASD,NVIR,Nchunks,NAUXBAS, &
           LddiActStart,LddiActEnd)
      use omp_lib
      use blasHandles
      use iso_c_binding
      integer(c_int) :: NVIR,NAUXBAS,Nchunks,NAUXBASD
      double precision :: B32(NVIR*NAUXBAS,LddiActStart:LddiActEnd)
      double precision :: T3(NAUXBAS*NVIR*Nchunks)
      double precision :: VDX(NAUXBASD,NAUXBAS)

#if __cublas==1 || __hipblas==1

!$omp target data use_device_ptr(B32,VDX,T3)
#if __cublas==1
         cublas_return = cublasXtDgemm(cublas_handle, &
             CUBLAS_OP_N,CUBLAS_OP_N,                 &
             NAUXBASD,NVIR*Nchunks,NAUXBAS,       &
             1.0D00, VDX,NAUXBASD,                   &
                     T3,NAUXBASD,                    &
             0.0D00, B32(1,LddiActStart),NAUXBASD)
         cublas_return = cudaDeviceSynchronize()
#endif
#if __hipblas==1
        m=NAUXBASD
        n=NVIR*Nchunks
        k=NAUXBAS
        !alpha
        da=c_loc(VDX(1,1))
        lda=NAUXBASD
        db=c_loc(T3(1))
        ldb=NAUXBASD
        !beta
        dc=c_loc(B32(1,LddiActStart))
        ldc=NAUXBASD

        call hipblasCheck(hipblasDgemm(HIPBLAS_handle, &
          HIPBLAS_OP_N, & ! integer(kind(hipblas_op_n)), value
          HIPBLAS_OP_N, & ! integer(kind(hipblas_op_n)), value
          m, & ! integer(c_int), value
          n, & ! integer(c_int), value
          k, & ! integer(c_int), value
          alpha, & ! real(c_double)
          da, & ! type(c_ptr), value
          lda, & ! integer(c_int), value
          db, & ! type(c_ptr), value
          ldb, & ! integer(c_int), value
          beta, & ! real(c_double)
          dc, & ! type(c_ptr), value
          ldc & ! integer(c_int), value
        ))
        call hipCheck(hipDeviceSynchronize())
#endif

!$omp end target data

#endif

      END SUBROUTINE RIMP2_FORM_B32



#if __cublas==1 || __hipblas==1
      SUBROUTINE RIMP2_FORM_B32_one(B32,T3,VDX,NAUXBASD,NVIR,NAUXBAS)
      use omp_lib
      use blasHandles
      use iso_c_binding
      integer(c_int) :: NVIR,NAUXBAS,NAUXBASD
      double precision :: B32(NVIR*NAUXBAS)
      double precision :: T3(NAUXBAS*NVIR)
      double precision :: VDX(NAUXBASD,NAUXBAS)

            
            ! ~~~ Start ~~~

!$omp target data use_device_ptr(B32,VDX,T3)
#if __cublas==1
         cublas_return = cublasXtDgemm(cublas_handle, &
             CUBLAS_OP_N,CUBLAS_OP_N,                 &
             NAUXBASD,NVIR,NAUXBAS,                   &
             1.0D00, VDX,NAUXBASD,                    &
                     T3,NAUXBASD,                     &
             0.0D00, B32,NAUXBASD)
         cublas_return = cudaDeviceSynchronize()
#endif
#if __hipblas==1

        m=NAUXBASD
        n=NVIR
        k=NAUXBAS
        !alpha
        da=c_loc(VDX(1,1))
        lda=NAUXBASD
        db=c_loc(T3(1))
        ldb=NAUXBASD
        !beta
        dc=c_loc(B32(1))
        ldc=NAUXBASD

        call hipblasCheck(hipblasDgemm(HIPBLAS_handle, &
         HIPBLAS_OP_N, & ! integer(kind(hipblas_op_n)), value
         HIPBLAS_OP_N, & ! integer(kind(hipblas_op_n)), value
         m, & ! integer(c_int), value
         n, & ! integer(c_int), value
         k, & ! integer(c_int), value
         alpha, & ! real(c_double)
         da, & ! type(c_ptr), value
         lda, & ! integer(c_int), value
         db, & ! type(c_ptr), value
         ldb, & ! integer(c_int), value
         beta, & ! real(c_double)
         dc, & ! type(c_ptr), value
         ldc & ! integer(c_int), value
        ))
        call hipCheck(hipDeviceSynchronize())
#endif

!$omp end target data

      END
#endif      


#if __cpu==1
      SUBROUTINE RIMP2_FORM_B32_one_cpu &
          (B32,T3,VDX,NAUXBASD,NVIR,NAUXBAS)
      use omp_lib

      double precision :: B32(NVIR*NAUXBAS)
      double precision :: T3(NAUXBAS*NVIR)
      double precision :: VDX(NAUXBASD,NAUXBAS)
      integer :: NVIR,NAUXBAS,NAUXBASD

         
         ! ~~~ Start ~~~


      CALL DGEMM &
         ('N','N', &
           NAUXBASD,NVIR,NAUXBAS, &
           1.0D00,VDX,NAUXBASD, &
                  T3,NAUXBASD, &
           0.0D00,B32,NAUXBASD)



      END
#endif




#if __intel_offloading==1
      SUBROUTINE RIMP2_FORM_B32_one_intel &
          (B32,T3,VDX,NAUXBASD,NVIR,NAUXBAS)
      use omp_lib
      use onemkl_blas_omp_offload_ilp64
      ! use onemkl_blas_omp_offload_lp64

      double precision :: B32(NVIR*NAUXBAS)
      double precision :: T3(NAUXBAS*NVIR)
      double precision :: VDX(NAUXBASD,NAUXBAS)
      integer :: NVIR,NAUXBAS,NAUXBASD


!$omp target variant dispatch use_device_ptr(B32,VDX,T3)
      CALL DGEMM &
         ('N','N', &
           NAUXBASD,NVIR,NAUXBAS, &
           1.0D00,VDX,NAUXBASD, &
                  T3,NAUXBASD, &
           0.0D00,B32,NAUXBASD)
!$omp end target variant dispatch

      END SUBROUTINE RIMP2_FORM_B32_one_intel

#endif






      SUBROUTINE RIMP2_B32_TRANSPOSE &
                (B32,T3,NAUXBAS,NVIR)
      use omp_lib
      
      implicit none

      double precision :: B32(NVIR,NAUXBAS)
      double precision :: T3(NAUXBAS,NVIR)
      integer :: NAUXBAS,NVIR
      integer :: IACT,IAUX,IVIR,IDX

!$omp target teams distribute parallel do collapse(2) default(none) &
!$omp shared(B32,T3) &
!$omp shared(NAUXBAS,NVIR) &
!$omp private(IAUX,IVIR)
            DO IAUX=1,NAUXBAS
               DO IVIR=1,NVIR
                  T3(IAUX,IVIR) = B32(IVIR,IAUX)
               ENDDO
            ENDDO
!$omp end target teams distribute parallel do

      END SUBROUTINE RIMP2_B32_TRANSPOSE
      ! ----------------------------------------------------------------




      SUBROUTINE RIMP2_B32_TRANSPOSE_cpu &
                (B32,T3,NAUXBAS,NVIR)
      use omp_lib
      implicit none

      double precision :: B32(NVIR,NAUXBAS)
      double precision :: T3(NAUXBAS,NVIR)
      integer :: NAUXBAS,NVIR
      integer :: IACT,IAUX,IVIR,IDX

!!$omp target teams distribute parallel do collapse(2) default(none) &
!!$omp shared(B32,T3) &
!!$omp shared(NAUXBAS,NVIR) &
!!$omp private(IAUX,IVIR)
            DO IAUX=1,NAUXBAS
               DO IVIR=1,NVIR
                  T3(IAUX,IVIR) = B32(IVIR,IAUX)
               ENDDO
            ENDDO
!!$omp end target teams distribute parallel do

      END SUBROUTINE RIMP2_B32_TRANSPOSE_cpu
      ! ----------------------------------------------------------------






#if __intel_offloading==1
      SUBROUTINE RIMP2_I32_TRANSF_mkl &
                (I32,VEC,T11, &
                 NCOR,NACT,NVIR,NBASIS,NAUXBAS, &
                 NddiChunk,NddiChunk_max, &
                 LddiAuxStart,LddiAuxEnd)

      use omp_lib
      use onemkl_blas_omp_offload_ilp64
      ! use onemkl_blas_omp_offload_lp64

      implicit none

      double precision :: I32(NBASIS*NBASIS*NddiChunk)
      double precision :: T11(NBASIS*NACT,NddiChunk_max)
      double precision :: VEC(NBASIS,NBASIS)

      integer :: &
         NCOR,NACT,NVIR,NBASIS,NAUXBAS, &
         NddiChunk,NddiChunk_max, &
         LddiAuxStart,LddiAuxEnd, LL




!       ---------------
!!$omp target data map(tofrom:I32) map(to:VEC) map(alloc:T11)
!!$omp target variant dispatch use_device_ptr(I32,T11,VEC)
      CALL DGEMM &
         ('T','N', &
           NACT,NBASIS*NddiChunk, NBASIS, &
           1.0D00,VEC(1,1+NCOR),NBASIS, &
                  I32,NBASIS, &
           0.0D00,T11,NACT)

!!          -----------------
!!$omp end target variant dispatch
!!$omp end target data

      DO LL=LddiAuxStart,LddiAuxEnd
         ! -----------------
!!$omp target variant dispatch use_device_ptr(VEC,T11,I32)         
         CALL DGEMM &
            ('T','T', &
             NVIR,NACT,NBASIS, &
             1.0D00,VEC(1,1+NCOR+NACT),NBASIS, &
                    T11(1,LL-LddiAuxStart+1),NACT, &
             0.0D00,I32(NVIR*NACT*(LL-LddiAuxStart)+1),NVIR)
!!             --------------
!!$omp end target variant dispatch
      ENDDO
!!       ----------------
!!$omp end target data

      END
#endif





#if __cpu==1
      SUBROUTINE RIMP2_I32_TRANSF_cpu &
                (I32,VEC,T11, &
                 NCOR,NACT,NVIR,NBASIS,NAUXBAS, &
                 NddiChunk,NddiChunk_max, &
                 LddiAuxStart,LddiAuxEnd)

      use omp_lib

      implicit none

      double precision :: I32(NBASIS*NBASIS*NddiChunk)
      double precision :: T11(NBASIS*NACT,NddiChunk_max)
      double precision :: VEC(NBASIS,NBASIS)

      integer :: &
         NCOR,NACT,NVIR,NBASIS,NAUXBAS, &
         NddiChunk,NddiChunk_max, &
         LddiAuxStart,LddiAuxEnd, LL



      CALL DGEMM &
         ('T','N', &
           NACT,NBASIS*NddiChunk, NBASIS, &
           1.0D00,VEC(1,1+NCOR),NBASIS, &
                  I32,NBASIS, &
           0.0D00,T11,NACT)

      DO LL=LddiAuxStart,LddiAuxEnd
         CALL DGEMM &
            ('T','T', &
             NVIR,NACT,NBASIS, &
             1.0D00,VEC(1,1+NCOR+NACT),NBASIS, &
                    T11(1,LL-LddiAuxStart+1),NACT, &
             0.0D00,I32(NVIR*NACT*(LL-LddiAuxStart)+1),NVIR)
      ENDDO


      END
#endif


      SUBROUTINE RIMP2_B32_STORED_DDI &
                (I32,NACT,NVIR,LddiAuxStart,LddiAuxEnd,NddiChunk)

      use gpuGlobalRimp2Eng,only: D_B
      implicit none
 
      double precision :: I32(NVIR,NACT,NddiChunk)
      double precision,allocatable :: T3(:,:)

      integer :: NACT,NVIR,LddiAuxStart,LddiAuxEnd,NddiChunk
      integer :: LL,IBG,IED


      ALLOCATE(T3(NVIR,NACT))
      DO LL=LddiAuxStart,LddiAuxEnd
         IBG=NVIR*(LL-1)+1
         IED=NVIR*LL
         CALL DDI_PUT(D_B,IBG,IED,1,NACT,I32(1:NVIR,1:NACT,LL-LddiAuxStart+1))
      ENDDO

      deallocate(T3)
      END SUBROUTINE RIMP2_B32_STORED_DDI
      ! ----------------------------------------------------------------



      SUBROUTINE RIMP2_B32_STORED_NODE &
                (B32,I32,NACT,NVIR, &
                 NAUXBAS,LddiAuxStart,LddiAuxEnd,NddiChunk)

      implicit none
 
      double precision,intent(out) :: B32(NVIR,NAUXBAS,NACT)
      double precision,intent(in) :: I32(NVIR,NACT,NddiChunk)

      integer,intent(in) :: NACT,NVIR,NAUXBAS,LddiAuxStart,LddiAuxEnd,NddiChunk

      integer :: LL,IACT,IVIR
            

            ! ~~~ Start ~~~


!$omp parallel do collapse(3) default(none) &
!$omp shared(I32,B32) &
!$omp shared(LddiAuxStart,LddiAuxEnd,NVIR,NACT) &
!$omp private(LL,IACT,IVIR)
      DO LL=LddiAuxStart,LddiAuxEnd
         DO IACT=1,NACT
            DO IVIR=1,NVIR
               B32(IVIR,LL,IACT) = I32(IVIR,IACT,LL-LddiAuxStart+1)
            ENDDO
         ENDDO
      ENDDO
!$omp end parallel do

      END SUBROUTINE RIMP2_B32_STORED_NODE
      ! ----------------------------------------------------------------


      ! SUBROUTINE RIMP2_I32_SLICING_gpu &
      !           (NSLICES_I32,NddiChunk_max, &
      !            KAUXTY,KAUXMI,LOCSPH, &
      !            NACT,NBASIS,LddiAuxShEnd,LddiAuxShStart)

      ! use gpuGlobalRimp2Eng,only: I32WSA
      ! use gpuGlobalRimp2Eng,only: MIXS


      ! integer :: KAUXTY(*),KAUXMI(*),LOCSPH(*)
      ! logical :: LOOP,LOMP


      ! ! INIT NSLICES_I32
      ! NAuxShChunk=LddiAuxShEnd-LddiAuxShStart+1      
      ! ! SIZE=NBASIS*(NBASIS+NACT)*NAuxShChunk/536870912.0
      !       SIZE=NBASIS*(NBASIS+NACT)*NAuxShChunk/936870912.0
      ! NSLICES_I32=CEILING(SIZE)
      ! IF(NSLICES_I32.GT.NAuxShChunk) NSLICES_I32=NAuxShChunk
      
      ! LOOP=.TRUE.
      ! DO WHILE(LOOP)
      !    IF(ALLOCATED(I32WSA)) DEALLOCATE(I32WSA)
      !    ALLOCATE(I32WSA(NSLICES_I32,2))
      !    CALL WorkSharedArray_gpu &
      !        (I32WSA,NSLICES_I32,LddiAuxShStart,LddiAuxShEnd,LOMP)

      !    ! SLICE OF MAX SIZE: NddiChunk_max
      !    LShStart=I32WSA(1,1)      
      !    LShEnd=I32WSA(1,2)     
      !    LANGK = KAUXTY(LShEnd)-1
      !    MINK  = KAUXMI(LShEnd)
      !    INDEX = LANGK+MINK
      !    NKT = MIXS(2,INDEX)-MIXS(1,INDEX)+1
      !    LStart = LOCSPH(LShStart)
      !    LEnd = LOCSPH(LShEnd)+NKT-1
      !    NddiChunk_max = LEnd-LStart+1
      !    DO II=2,NSLICES_I32       
      !       LShStart=I32WSA(II,1)      
      !       LShEnd=I32WSA(II,2)     
      !       LANGK = KAUXTY(LShEnd)-1
      !       MINK  = KAUXMI(LShEnd)
      !       INDEX = LANGK+MINK
      !       NKT = MIXS(2,INDEX)-MIXS(1,INDEX)+1
      !       LStart = LOCSPH(LShStart)
      !       LEnd = LOCSPH(LShEnd)+NKT-1
      !       MaxTmp = LEnd-LStart+1
      !       IF(NddiChunk_max.LT.MaxTmp) NddiChunk_max=MaxTmp
      !    ENDDO

      !    ! TUNING NSLICES_I32
      !    SIZE=(NBASIS+NACT)*NBASIS*NddiChunk_max/134217728.0
      !    LOOP=.FALSE.
      !    IF(SIZE.GT.12.0) THEN
      !       LOOP=.TRUE.
      !       NSLICES_I32=NSLICES_I32+1
      !    ENDIF
      ! ENDDO

      ! END
      ! ! ----------------------------------------------------------------


      SUBROUTINE RIMP2_I32_SLICING_alter_gpu &
                (NSLICES_I32,NddiChunk_max, &
                 KAUXTY,KAUXMI,LOCSPH, &
                 NACT,NBASIS,LddiAuxShEnd,LddiAuxShStart)

      use gpuGlobalRimp2Eng,only: I32WSA
      use gpuGlobalRimp2Eng,only: MIXS

      implicit none

      integer :: NACT,NBASIS
      integer :: LddiAuxShEnd,LddiAuxShStart
      integer :: NddiChunk_max
      integer :: NSLICES_I32

      integer :: KAUXTY(*),KAUXMI(*),LOCSPH(*)

      integer :: NAuxShChunk

      double precision :: SIZE

      logical :: LOOP,LOMP
      integer :: total,free

      integer :: LShStart, LShEnd

integer :: LANGK,MINK,INDEX,NKT
integer :: LSTART,LEND
integer :: II,MAXTMP



      CALL RIMP2_MEM_gpu(free,total)

      ! INIT NSLICES_I32
      NAuxShChunk=LddiAuxShEnd-LddiAuxShStart+1

      ! I32(NBASIS*NBASIS*NAuxShChunk/Nchunk)
      ! T11(NBASIS*NACT*NAuxShChunk/Nchunk)
      ! 80% gpu mem
      ! 8*(NBASIS+NACT)*NBASIS*NAuxShChunk/Nchunk = free*0.8 
      ! Nchunk = (NBASIS+NACT)*NBASIS*NAuxShChunk/(free*10)

! (NBASIS*NBASIS*NAuxShChunk + NBASIS*NACT*NAuxShChunk)*8


      SIZE=10.0*real(NBASIS+NACT)*real(NBASIS)*real(NAuxShChunk)/real(free)
      NSLICES_I32=CEILING(SIZE)

! write(*,*) "NSLICES_I32: ",NSLICES_I32,SIZE,NBASIS,NACT,NAuxShChunk

      ! INIT NSLICES_I32
      IF(NSLICES_I32.GT.NAuxShChunk) NSLICES_I32=NAuxShChunk
      IF(ALLOCATED(I32WSA)) DEALLOCATE(I32WSA)
      ALLOCATE(I32WSA(NSLICES_I32,2))
      CALL WorkSharedArray_gpu(I32WSA,NSLICES_I32,LddiAuxShStart,LddiAuxShEnd,LOMP)

      
      LOOP=.TRUE.
      DO WHILE(LOOP)
         IF(ALLOCATED(I32WSA)) DEALLOCATE(I32WSA)
         ALLOCATE(I32WSA(NSLICES_I32,2))
         CALL WorkSharedArray_gpu &
             (I32WSA,NSLICES_I32,LddiAuxShStart,LddiAuxShEnd,LOMP)

         ! SLICE OF MAX SIZE: NddiChunk_max
         LShStart=I32WSA(1,1)      
         LShEnd=I32WSA(1,2)     
         LANGK = KAUXTY(LShEnd)-1
         MINK  = KAUXMI(LShEnd)
         INDEX = LANGK+MINK
         NKT = MIXS(2,INDEX)-MIXS(1,INDEX)+1
         LStart = LOCSPH(LShStart)
         LEnd = LOCSPH(LShEnd)+NKT-1
         NddiChunk_max = LEnd-LStart+1
         DO II=2,NSLICES_I32       
            LShStart=I32WSA(II,1)      
            LShEnd=I32WSA(II,2)     
            LANGK = KAUXTY(LShEnd)-1
            MINK  = KAUXMI(LShEnd)
            INDEX = LANGK+MINK
            NKT = MIXS(2,INDEX)-MIXS(1,INDEX)+1
            LStart = LOCSPH(LShStart)
            LEnd = LOCSPH(LShEnd)+NKT-1
            MaxTmp = LEnd-LStart+1
            IF(NddiChunk_max.LT.MaxTmp) NddiChunk_max=MaxTmp
         ENDDO

         ! TUNING NSLICES_I32
         SIZE=(NBASIS+NACT)*NBASIS*NddiChunk_max*8
         LOOP=.FALSE.
         IF(SIZE.GT.free*0.8) THEN
            LOOP=.TRUE.
            NSLICES_I32=NSLICES_I32+1
         ENDIF
      ENDDO


      END
      ! ----------------------------------------------------------------







#if __cublas==2
      subroutine Cholesky(V,lda,n)
      use omp_lib
      use iso_c_binding
      use cuda_cusolver
! parameters      
      integer*8,parameter :: cudaMemcpyDeviceToHost=0
      integer*8,parameter :: cudaMemcpyHostToDevice=1
      integer*4,parameter :: CUBLAS_OP_N=0
      integer*4,parameter :: CUBLAS_OP_T=1
      integer(c_int),parameter :: CUBLAS_FILL_MODE_LOWER=1
      integer(c_int),parameter :: CUBLAS_FILL_MODE_UPPER=0
      integer(c_int),parameter :: CUBLAS_DIAG_NON_UNIT = 0
      integer(c_int),parameter :: CUBLAS_DIAG_UNIT = 1
! V matrix
      real(c_double),target :: V(lda,n)
      real(c_double),allocatable :: WS(:),WS_i(:)
      integer(c_int) n,lda
! memory operation parameters
      integer(c_int) :: V_mem_stat
      type(c_ptr) :: d_V
      type(c_ptr) :: CPU_V_ptr
      integer(c_size_t) :: V_size
! buffer parameters
      integer(c_int),target :: Lwork, Lwork_i
      integer :: devInfo,devInfo_i
      integer*8 :: devInfo_size,Lwork_size,devInfo_size_i,Lwork_size_i
      integer*8 :: Workspace,Workspace_i
      type(c_ptr) :: d_Lwork,d_Lwork_i
      type(c_ptr) :: d_WS,d_WS_i
      type(c_ptr) :: d_devInfo,d_devInfo_i
      type(c_ptr)::CPU_Lwork_ptr, CPU_Lwork_ptr_i
! cusolver handle 
      type(c_ptr) :: cusolver_Hndl
! function status return  
      integer :: cusolver_stat
      integer :: WS_mem_stat, WS_mem_stat_i
      integer :: devInfo_mem_stat, devInfo_mem_stat_i
      integer :: Lwork_mem_stat, Lwork_mem_stat_i
! CPU information     
      type(c_ptr)::cpfre,cptot
      integer*8,target::free,total
      integer res

#if __cublas==1

      free = 0
      total = 0
      res = 1
      cpfre = c_loc(free)
      cptot = c_loc(total)
  
      res = cudaMemGetInfo(cpfre,cptot)
      ! write (*, '(A, I12)') "  free mem: ", free
      ! write (*, '(A, I12)') " total mem: ", total



! prepare V on host and device
      V_size = SIZEOF(V)
      CPU_V_ptr = C_LOC(V)
      ! V_mem_stat = cudaMalloc(d_V,V_size)
      ! V_mem_stat = cudaMemcpy(d_V,CPU_V_ptr,V_size,cudaMemcpyHostToDevice)

! prepare info arrays
      devInfo_size = SIZEOF(devInfo)
      Lwork_size = SIZEOF(Lwork)
      devInfo_size_i = SIZEOF(devInfo_i)
      Lwork_size_i = SIZEOF(Lwork_i)
      CPU_Lwork_ptr = C_Loc(Lwork)
      CPU_Lwork_ptr_i = C_Loc(Lwork_i)
      devInfo_mem_stat = cudaMalloc(d_devInfo,devInfo_size)
      Lwork_mem_stat   = cudaMalloc(d_Lwork,Lwork_size)
      devInfo_mem_stat_i = cudaMalloc(d_devInfo_i,devInfo_size_i)
      Lwork_mem_stat_i   = cudaMalloc(d_Lwork_i,Lwork_size_i)


      res = cudaMemGetInfo(cpfre,cptot)
      ! write (*, '(A, I12)') "  free mem: ", free
      ! write (*, '(A, I12)') " total mem: ", total

  

!$omp target data map(tofrom:V) 
! start cusolver handle
      cusolver_stat = cusolverDnCreate(cusolver_Hndl)  

!$omp target data use_device_ptr(V)
! Cholesky decomposition buffer
      cusolver_stat = cusolverDnDpotrf_bufferSize &
                     (cusolver_Hndl, &
                      CUBLAS_FILL_MODE_UPPER, &
                      n, &
                      V, &
                      lda, &
                      CPU_Lwork_ptr)
!$omp end target data
      allocate(WS(Lwork))

!$omp target data use_device_ptr(V)
! triangular matrix inversion buffer
      cusolver_stat = cusolverDnDtrtri_bufferSize &
                     (cusolver_Hndl, &
                      CUBLAS_FILL_MODE_UPPER, &
                      CUBLAS_DIAG_NON_UNIT, &
                      n, &
                      V, &
                      lda, &
                      CPU_Lwork_ptr_i)
!$omp end target data
      allocate(WS_i(Lwork_i))


!$omp target data map(alloc:WS,WS_i)

! Cholesky decomposition
!$omp target data use_device_ptr(V,WS)
      cusolver_stat = cusolverDnDpotrf &
                     (cusolver_Hndl, &
                      CUBLAS_FILL_MODE_UPPER, &
                      n, &
                      V, &
                      lda, &
                      WS, &
                      Lwork, &
                      d_devInfo)
      cusolver_stat = cudaDeviceSynchronize()      
!$omp end target data

! triangular matrix inversion
!$omp target data use_device_ptr(V,WS_i)
      cusolver_stat = cusolverDnDtrtri &
                     (cusolver_Hndl, &
                      CUBLAS_FILL_MODE_UPPER, &
                      CUBLAS_DIAG_NON_UNIT, &
                      n, &
                      V, &
                      lda, &
                      WS_i, &
                      Lwork_i, &
                      d_devInfo_i)
!$omp end target data

 
! close cusolver handle  
      cusolver_stat = cusolverDnDestroy(cusolver_Hndl)
!$omp end target data
!$omp end target data
 


#endif
deallocate(WS)
deallocate(WS_i)
      end subroutine Cholesky
#endif


      SUBROUTINE RIMP2_CHOLESKY_DECOMPOSE_VXX_GPU(V,NAUXBAS)
      use omp_lib
#if __intel_offloading==1
      use onemkl_blas_omp_offload_ilp64
      ! use onemkl_blas_omp_offload_lp64
#endif

      use iso_c_binding
      implicit none

      integer :: IAUX,JAUX, INFO

      integer(c_int) :: NAUXBAS
#if __hipblas==1
      real(c_double) :: V(NAUXBAS,NAUXBAS)
#else
      real(c_double),target :: V(NAUXBAS,NAUXBAS)
#endif

! #if __cublas==1
!       CALL Cholesky(V,NAUXBAS,NAUXBAS)
! #elif __hipblas==1
      ! CHOLESKY DECOMPOSION OF 2C-2E MATRIX V=L*LT
      CALL DPOTRF('L',NAUXBAS,V,NAUXBAS,INFO)
      ! DETERMINATION OF INVERSE OF CHOLESKY DECOMPOSED MATRIX L^(-1)
      CALL DTRTRI('L','N',NAUXBAS,V,NAUXBAS,INFO)
! #endif




! #if __cublas == 1 || __cpu == 1
#if __cpu == 1
      ! CHOLESKY DECOMPOSION OF 2C-2E MATRIX V=L*LT
      CALL DPOTRF('L',NAUXBAS,V,NAUXBAS,INFO)
      ! DETERMINATION OF INVERSE OF CHOLESKY DECOMPOSED MATRIX L^(-1)
      CALL DTRTRI('L','N',NAUXBAS,V,NAUXBAS,INFO)
#endif

! //////////////////////////
#if __intel_offloading==1
w0=omp_get_wtime()
! -------------
!!$omp target data map(tofrom:V) 
!!$omp target variant dispatch use_device_ptr(V)
      ! CHOLESKY DECOMPOSION OF 2C-2E MATRIX V=L*LT
      CALL DPOTRF('L',NAUXBAS,V,NAUXBAS,INFO)
!       ---------------
!!$omp end target variant dispatch
!!$omp end target data

!!$omp target data map(tofrom:V) 
!!$omp target variant dispatch use_device_ptr(V)
      ! DETERMINATION OF INVERSE OF CHOLESKY DECOMPOSED MATRIX L^(-1)     
      CALL DTRTRI('L','N',NAUXBAS,V,NAUXBAS,INFO)
!       ------------
!!$omp end target variant dispatch
!!$omp end target data
w1=omp_get_wtime()
write(*,*) "INTEL routine: DPOTRF & DTRTRI ",w1-w0
#endif

!!!!! ZERO OUT LOWER PART
      IF(NAUXBAS .GT. 1) THEN
         DO JAUX = 2, NAUXBAS
            DO IAUX = 1, JAUX-1
               V(IAUX,JAUX) = 0.0D00
            ENDDO
         ENDDO
      ENDIF


      END
      !***********************************************************

      SUBROUTINE RIMP2_MEM_gpu(free,total)
! looking up device mem info is slightly different for hipblas.
         use iso_c_binding
         use blasHandles
         use mod_device, only: freebytes
         USE COMM_PAR, ONLY: MASWRK
         USE COMM_IOFILE, ONLY: IW
         implicit none

         type(c_ptr)::cpfre,cptot
#if __hipblas==1
         integer(c_size_t) :: free,total,freeapi
         integer(kind(hipSuccess)) :: res
#else
  integer*8,target :: free,total,freeapi
        integer :: res
#endif
         free = 0
         total = 0
         res = 1
#if __cublas==1
         cpfre = c_loc(free)
         cptot = c_loc(total)
         res = cudaMemGetInfo(cpfre,cptot)
#endif

#if __hipblas==1
         res = hipMemGetInfo(free,total)
#endif

#if __cublas==1 || __hipblas==1
         freeapi = free
#if __PRINT_ESS==1
         if(maswrk) write(iw,'(55(1H*)/,A,":",1X,F8.3,1X,"GB")')  "DEVICE FREE  MEMORY (API)",freeapi/1073741824.d0
#endif
         if (freebytes > 0) then
            free = min(freebytes,freeapi)
#if __PRINT_ESS==1
            if(maswrk) write(iw,'(A,":",1X,F8.3,1X,"GB")')        "DEVICE FREE  MEMORY (INP)",free/1073741824.d0
#endif
         endif
#if __PRINT_ESS==1
         if(maswrk) write(iw,'(A,":",1X,F8.3,1X,"GB"/,55(1H*)/)') "DEVICE TOTAL MEMORY (API)",total/1073741824.d0
         flush(iw)
#endif
#endif

#if __cpu==1
         free = 16914055168
         total = 16914055168
#endif

#if __intel_offloading==1
         free = 16914055168
         total = 16914055168
#endif

#if __PRINT_TIME_BUG == 1
         write (*, '(A, I12)') "  free mem device: ", free
         write (*, '(A, I12)') " total mem: ", total
#endif
      END SUBROUTINE

      ! SUBROUTINE RIMP2_reset_gpu()
      !    use iso_c_binding
      !    use cublasf
      !    implicit none
      !    integer :: res
 
      !    res = cudaDeviceReset()
      ! END SUBROUTINE


#if __cublas==1 || __hipblas==1
   SUBROUTINE RIMP2_QVV_HEAD_cublas &
             (QVV,BI,BJ,NACT,NVIR,NAUXBASD,iQVV)


   use omp_lib
   use blasHandles
   use iso_c_binding

   implicit none

   ! input
   integer(c_int),intent(in) :: NACT,NVIR,NAUXBASD,iQVV
   double precision,intent(in) :: BI(NAUXBASD*NVIR*iQVV)
   double precision,intent(in) :: BJ(NAUXBASD*NVIR)

   ! output
   double precision,intent(out) :: QVV(NVIR,iQVV,NVIR)


       ! ~~~ start ~~~


! compute 4-2ERIs
!$omp target data use_device_ptr(BI,BJ,QVV)
#if __cublas==1
   cublas_return =  cublasXtDgemm  &
      (cublas_handle,CUBLAS_OP_T,CUBLAS_OP_N, &
       NVIR*iQVV,NVIR,NAUXBASD,               &
       1.0D00,BI,NAUXBASD,                    &
              BJ,NAUXBASD,                    &
       0.0D00,QVV,NVIR*iQVV)
   cublas_return = cudaDeviceSynchronize()
#endif
#if __hipblas==1

        m=NVIR*iQVV
        n=NVIR
        k=NAUXBASD
        !alpha
        da=c_loc(BI(1))
        lda=NAUXBASD
        db=c_loc(BJ(1))
        ldb=NAUXBASD
        !beta
        dc=c_loc(QVV(1,1,1))
        ldc=NVIR*iQVV

        call hipblasCheck(hipblasDgemm(HIPBLAS_handle, &
         HIPBLAS_OP_T, & ! integer(kind(hipblas_op_n)), value
         HIPBLAS_OP_N, & ! integer(kind(hipblas_op_n)), value
         m, & ! integer(c_int), value
         n, & ! integer(c_int), value
         k, & ! integer(c_int), value
         alpha, & ! real(c_double)
         da, & ! type(c_ptr), value
         lda, & ! integer(c_int), value
         db, & ! type(c_ptr), value
         ldb, & ! integer(c_int), value
         beta, & ! real(c_double)
         dc, & ! type(c_ptr), value
         ldc & ! integer(c_int), value
        ))
        call hipCheck(hipDeviceSynchronize())
#endif
!$omp end target data


   END SUBROUTINE RIMP2_QVV_HEAD_cublas
#endif




#if __intel_offloading==1
   SUBROUTINE RIMP2_QVV_HEAD_mkl &
             (QVV,BI,BJ,NACT,NVIR,NAUXBASD,iQVV)


   use omp_lib
   use onemkl_blas_omp_offload_ilp64

   implicit none

   ! output
   double precision,intent(out) :: QVV(NVIR,iQVV,NVIR)

   ! input
   double precision,intent(in) :: BI(NAUXBASD*NVIR*iQVV)
   double precision,intent(in) :: BJ(NAUXBASD*NVIR)
   integer,intent(in) :: NACT,NVIR,NAUXBASD,iQVV


       ! ~~~ start ~~~


! compute 4-2ERIs
!$omp target variant dispatch use_device_ptr(BI,BJ,QVV)
   CALL DGEMM &
      ('T','N', &
       NVIR*iQVV,NVIR,NAUXBASD,               &
       1.0D00,BI,NAUXBASD,                    &
              BJ,NAUXBASD,                    &
       0.0D00,QVV,NVIR*iQVV)
!$omp end target variant dispatch


   END SUBROUTINE RIMP2_QVV_HEAD_mkl
#endif



      SUBROUTINE RIMP2_ENERGYIJ_cpu &
                (E2,BI,BJ,eij,eab,EIG,                     &
                 QVV,IACT,JACT,NACT,NVIR,NAUXBASD,iQVV)

      use omp_lib
      
      implicit none

      ! output
      double precision :: E2

      ! input
      double precision :: BI(NAUXBASD*NVIR*iQVV)
      double precision :: BJ(NAUXBASD*NVIR)
      double precision :: eab(NVIR,NVIR)
      double precision :: eij(NACT,NACT)
      double precision :: EIG(NACT+NVIR)

      ! buffer
      double precision :: QVV(NVIR,iQVV,NVIR)

integer :: IACT,JACT,NACT,NVIR,NAUXBASD,iQVV, IC, IB, IA 
double precision ::  FAC,Tijab,Q_t,E2_t




      CALL DGEMM &
         ('T','N', &
           NVIR*iQVV,NVIR,NAUXBASD, &
           1.0D00,BI,NAUXBASD, &
                  BJ,NAUXBASD, &
           0.0D00,QVV,NVIR*iQVV)


!!$omp target teams distribute parallel do collapse(3) &
!!$omp map(to:QVV) reduction(+:E2) 
      DO IC=1,iQVV
         DO IB=1,NVIR
            DO IA=1,NVIR
               FAC=2.0D00
               IF(IACT+IC-1.EQ.JACT) FAC=1.0D00            
               Tijab=QVV(IA,IC,IB)/(EIG(IACT+IC-1)+EIG(JACT)-EIG(IA+NACT)-EIG(IB+NACT))
               Q_t=QVV(IA,IC,IB)+QVV(IA,IC,IB)
               E2_t=Tijab*(Q_t-QVV(IB,IC,IA))
               E2 = E2 + FAC*E2_t      
            ENDDO
         ENDDO
      ENDDO
!!$omp end target teams distribute parallel do 

      END





      SUBROUTINE RIMP2_ENERGYIJ_intel &
                (E2,BI,BJ,eij,eab,EIG, &
                 QVV,IACT,JACT,NACT,NVIR,NAUXBASD,iQVV)

      use omp_lib
#if __intel_offloading==1
      use onemkl_blas_omp_offload_ilp64
      ! use onemkl_blas_omp_offload_lp64
#endif
      implicit none

      ! output
      double precision :: E2

      ! input
      double precision :: BI(NAUXBASD*NVIR*iQVV)
      double precision :: BJ(NAUXBASD*NVIR)
      double precision :: eab(NVIR,NVIR)
      double precision :: eij(NACT,NACT)
      double precision :: EIG(NACT+NVIR)

      ! buffer
      double precision :: QVV(NVIR,iQVV,NVIR)

integer :: IACT,JACT,NACT,NVIR,NAUXBASD,iQVV

#if __intel_offloading==1


!$omp target data map(to:BI,BJ)
!$omp target variant dispatch use_device_ptr(BI,BJ,QVV)
      CALL DGEMM &
         ('T','N', &
           NVIR*iQVV,NVIR,NAUXBASD, &
           1.0D00,BI,NAUXBASD, &
                  BJ,NAUXBASD, &
           0.0D00,QVV,NVIR*iQVV)
!$omp end target variant dispatch


!$omp target teams distribute parallel do collapse(3) &
!$omp reduction(+:E2)
      DO IC=1,iQVV
         DO IB=1,NVIR
            DO IA=1,NVIR
               FAC=2.0D00
               IF(IACT+IC-1.EQ.JACT) FAC=1.0D00            
               Tijab=QVV(IA,IC,IB)/(EIG(IACT+IC-1)+EIG(JACT)-EIG(IA+NACT)-EIG(IB+NACT))
               Q_t=QVV(IA,IC,IB)+QVV(IA,IC,IB)
               E2_t=Tijab*(Q_t-QVV(IB,IC,IA))
               E2 = E2 + FAC*E2_t      
            ENDDO
         ENDDO
      ENDDO
!$omp end target teams distribute parallel do 
!$omp end target data



#endif
      END



#if __cublas==1 || __hipblas==1
      SUBROUTINE RIMP2_QVV_TAIL_cublas &
                (QVV,BI,BJ,NVIR,NAUXBASD,iQVV)

      use omp_lib
      use blasHandles
      use iso_c_binding

      implicit none

      ! input
      integer(c_int),intent(in) :: NVIR,NAUXBASD,iQVV
      double precision,intent(in) :: BI(NAUXBASD*NVIR)
      double precision,intent(in) :: BJ(NAUXBASD*NVIR*iQVV)

      ! output
      double precision,intent(out) :: QVV(NVIR,NVIR,iQVV)


         ! ~~~ start ~~~


! compute 4-2ERIs
!$omp target data use_device_ptr(BI,BJ,QVV)
#if __cublas==1
      cublas_return =  cublasXtDgemm  &
         (cublas_handle,CUBLAS_OP_T,CUBLAS_OP_N,   &
          NVIR,NVIR*iQVV,NAUXBASD,                 &
          1.0D00,BI,NAUXBASD,                      &
                 BJ,NAUXBASD,                      &
          0.0D00,QVV,NVIR)
      cublas_return = cudaDeviceSynchronize()
#endif
#if __hipblas==1
        m=NVIR
        n=NVIR*iQVV
        k=NAUXBASD
        !alpha
        da=c_loc(BI(1))
        lda=NAUXBASD
        db=c_loc(BJ(1))
        ldb=NAUXBASD
        !beta
        dc=c_loc(QVV(1,1,1))
        ldc=NVIR

        call hipblasCheck(hipblasDgemm(HIPBLAS_handle, &
         HIPBLAS_OP_T, & ! integer(kind(hipblas_op_n)), value
         HIPBLAS_OP_N, & ! integer(kind(hipblas_op_n)), value
         m, & ! integer(c_int), value
         n, & ! integer(c_int), value
         k, & ! integer(c_int), value
         alpha, & ! real(c_double)
         da, & ! type(c_ptr), value
         lda, & ! integer(c_int), value
         db, & ! type(c_ptr), value
         ldb, & ! integer(c_int), value
         beta, & ! real(c_double)
         dc, & ! type(c_ptr), value
         ldc & ! integer(c_int), value
        ))
        call hipCheck(hipDeviceSynchronize())
#endif
!$omp end target data


      END SUBROUTINE RIMP2_QVV_TAIL_cublas
#endif





#if __intel_offloading==1
      SUBROUTINE RIMP2_QVV_TAIL_mkl &
                (QVV,BI,BJ,NVIR,NAUXBASD,iQVV)

      use omp_lib
      use onemkl_blas_omp_offload_ilp64
      ! use onemkl_blas_omp_offload_lp64

      implicit none

      ! output
      double precision,intent(out) :: QVV(NVIR,NVIR,iQVV)

      ! input
      double precision,intent(in) :: BI(NAUXBASD*NVIR)
      double precision,intent(in) :: BJ(NAUXBASD*NVIR*iQVV)
      integer,intent(in) :: NVIR,NAUXBASD,iQVV


         ! ~~~ start ~~~


! compute 4-2ERIs

!$omp target variant dispatch use_device_ptr(BI,BJ,QVV)
      CALL DGEMM  &
         ('T','N',   &
          NVIR,NVIR*iQVV,NAUXBASD,                 &
          1.0D00,BI,NAUXBASD,                      &
                 BJ,NAUXBASD,                      &
          0.0D00,QVV,NVIR)
!$omp end target variant dispatch


      END SUBROUTINE RIMP2_QVV_TAIL_mkl
#endif





      SUBROUTINE RIMP2_ENERGYIJ_II_cpu &
                (E2,BI,BJ,eij,eab,EIG,                    &
                 QVV,IACT,JACT,NACT,NVIR,NAUXBASD,iQVV)

      use omp_lib


      implicit none

      ! output
      double precision :: E2

      ! input
      double precision :: BI(NAUXBASD*NVIR)
      double precision :: BJ(NAUXBASD*NVIR*iQVV)
      double precision :: eab(NVIR,NVIR)
      double precision :: eij(NACT,NACT)
      double precision :: EIG(NACT+NVIR)

      ! buffer
      double precision :: QVV(NVIR,NVIR,iQVV)

integer :: IACT,JACT,NACT,NVIR,NAUXBASD,iQVV
integer :: IA,IB,IC
double precision :: FAC, Tijab, Q_t, E2_t




      CALL DGEMM &
         ('T','N', & 
           NVIR,NVIR*iQVV,NAUXBASD, &
           1.0D00,BI,NAUXBASD, &
                  BJ,NAUXBASD, &
           0.0D00,QVV,NVIR)


!!$omp target teams distribute parallel do collapse(3) &
!!$omp map(to:QVV) reduction(+:E2) 
      DO IC=1,iQVV
         DO IB=1,NVIR
            DO IA=1,NVIR
               FAC=2.0D00
               IF(IACT.EQ.JACT+IC-1) FAC=1.0D00
               Tijab=QVV(IA,IB,IC)/(EIG(IACT)+EIG(JACT+IC-1)-EIG(IA+NACT)-EIG(IB+NACT))
               Q_t=QVV(IA,IB,IC)+QVV(IA,IB,IC)
               E2_t=Tijab*(Q_t-QVV(IB,IA,IC))
               E2=E2+FAC*E2_t
            ENDDO
         ENDDO
      ENDDO
!!$omp end target teams distribute parallel do



      END





      SUBROUTINE RIMP2_ENERGYIJ_II_intel &
                (E2,BI,BJ,eij,eab,EIG,                    &
                 QVV,IACT,JACT,NACT,NVIR,NAUXBASD,iQVV)

      use omp_lib
#if __intel_offloading==1
      use onemkl_blas_omp_offload_ilp64
      ! use onemkl_blas_omp_offload_lp64
#endif

      implicit none

      ! output
      double precision :: E2

      ! input
      double precision :: BI(NAUXBASD*NVIR)
      double precision :: BJ(NAUXBASD*NVIR*iQVV)
      double precision :: eab(NVIR,NVIR)
      double precision :: eij(NACT,NACT)
      double precision :: EIG(NACT+NVIR)

      ! buffer
      double precision :: QVV(NVIR,NVIR,iQVV)

integer :: IACT,JACT,NACT,NVIR,NAUXBASD,iQVV


#if __intel_offloading==1


!!$omp target data map(to:BI,BJ) 
!$omp target variant dispatch use_device_ptr(BI,BJ,QVV)
      CALL DGEMM &
         ('T','N', & 
           NVIR,NVIR*iQVV,NAUXBASD, &
           1.0D00,BI,NAUXBASD, &
                  BJ,NAUXBASD, &
           0.0D00,QVV,NVIR)
!$omp end target variant dispatch

!$omp target teams distribute parallel do collapse(3) &
!$omp reduction(+:E2) 
      DO IC=1,iQVV
         DO IB=1,NVIR
            DO IA=1,NVIR
               FAC=2.0D00
               IF(IACT.EQ.JACT+IC-1) FAC=1.0D00
               Tijab=QVV(IA,IB,IC)/(EIG(IACT)+EIG(JACT+IC-1)-EIG(IA+NACT)-EIG(IB+NACT))
               Q_t=QVV(IA,IB,IC)+QVV(IA,IB,IC)
               E2_t=Tijab*(Q_t-QVV(IB,IA,IC))
               E2=E2+FAC*E2_t
            ENDDO
         ENDDO
      ENDDO
!$omp end target teams distribute parallel do
!!$omp end target data

#endif

      END




      SUBROUTINE RIMP2_I32_TRANSF_cublas &
                (I32,VEC,T11, &
                 NCOR,NACT,NVIR,NBASIS,NAUXBAS, &
                 NddiChunk,NddiChunk_max, &
                 LddiAuxStart,LddiAuxEnd)

      use omp_lib
      use blasHandles
      use iso_c_binding

      use gpuGlobalRimp2Eng,only: zero,one

      implicit none

      integer :: NAUXBAS
      
      integer(c_int) :: &
         NCOR,NACT,NVIR,NBASIS, & 
         NddiChunk,NddiChunk_max, &
         LddiAuxStart,LddiAuxEnd

      double precision :: I32(NBASIS*NBASIS*NddiChunk)
      double precision :: T11(NBASIS*NACT,NddiChunk_max)
      double precision :: VEC(NBASIS,NBASIS)

      integer :: LL

        


        ! ~~~ start ~~~


#if __cublas==1
cublas_return = cublasCreate_v2(cublas_handle)
#endif
#if __hipblas==1
 call hipblasCheck(hipblasCreate(hipblas_handle))
#endif        

!$omp target data use_device_ptr(I32,VEC,T11)
#if __cublas==1
      cublas_return=cublasDgemm_v2 &
         (cublas_handle,CUBLAS_OP_T,CUBLAS_OP_N,   &
          NACT,NBASIS*NddiChunk, NBASIS,       &
          one, VEC(1,1+NCOR),NBASIS,           &
                  I32,NBASIS,                     &
          zero, T11,NACT)
      cublas_return=cudaDeviceSynchronize()
#endif
#if __hipblas==1
        m=NACT
        n=NBASIS*NddiChunk
        k=NBASIS
        !alpha
        da=c_loc(VEC(1,1+NCOR))
        lda=NBASIS
        db=c_loc(I32(1))
        ldb=NBASIS
        !beta
        dc=c_loc(T11(1,1))
        ldc=NACT

        call hipblasCheck(hipblasDgemm(HIPBLAS_handle, &
         HIPBLAS_OP_T, & ! integer(kind(hipblas_op_n)), value
         HIPBLAS_OP_N, & ! integer(kind(hipblas_op_n)), value
         m, & ! integer(c_int), value
         n, & ! integer(c_int), value
         k, & ! integer(c_int), value
         alpha, & ! real(c_double)
         da, & ! type(c_ptr), value
         lda, & ! integer(c_int), value
         db, & ! type(c_ptr), value
         ldb, & ! integer(c_int), value
         beta, & ! real(c_double)
         dc, & ! type(c_ptr), value
         ldc & ! integer(c_int), value
        ))
        call hipCheck(hipDeviceSynchronize())
#endif
      DO LL=LddiAuxStart,LddiAuxEnd
#if __cublas==1
         cublas_return=cublasDgemm_v2 &
            (cublas_handle,CUBLAS_OP_T,CUBLAS_OP_T,      &
             NVIR,NACT,NBASIS,                        &
             one,VEC(1,1+NCOR+NACT),NBASIS,         &
               T11(1,LL-LddiAuxStart+1),NACT,         &
             zero,I32(NVIR*NACT*(LL-LddiAuxStart)+1),NVIR)
         cublas_return=cudaDeviceSynchronize()
#endif
#if __hipblas==1
        m=NVIR
        n=NACT
        k=NBASIS
        !alpha
        da=c_loc(VEC(1,1+NCOR+NACT))
        lda=NBASIS
        db=c_loc(T11(1,LL-LddiAuxStart+1))
        ldb=NACT
        !beta
        dc=c_loc(I32(NVIR*NACT*(LL-LddiAuxStart)+1))
        ldc=NVIR

        call hipblasCheck(hipblasDgemm(HIPBLAS_handle, &
         HIPBLAS_OP_T, & ! integer(kind(hipblas_op_n)), value
         HIPBLAS_OP_T, & ! integer(kind(hipblas_op_n)), value
         m, & ! integer(c_int), value
         n, & ! integer(c_int), value
         k, & ! integer(c_int), value
         alpha, & ! real(c_double)
         da, & ! type(c_ptr), value
         lda, & ! integer(c_int), value
         db, & ! type(c_ptr), value
         ldb, & ! integer(c_int), value
         beta, & ! real(c_double)
         dc, & ! type(c_ptr), value
         ldc & ! integer(c_int), value
        ))
        call hipCheck(hipDeviceSynchronize())
#endif
      ENDDO
#if __cublas==1
cublas_return = cublasDestroy_v2(cublas_handle)
#endif
#if __hipblas==1
call hipblasCheck(hipblasDestroy(hipblas_handle))
#endif
!$omp end target data

#if __PRINT_TIME_BUG==1
   write(*,*) "NVIDIA routine: RIMP2_I32_TRANSF_gpu"
#endif


      END !****************************************************


SUBROUTINE RIMP2_TAIL &
          (E2,BI,BJ,VDX,EIG,QVV,T3, &
           NACT,NVIR,IACT,LGPU_S,LGPU_E, &
           NQVV,NAUXBASD,NAUXBAS)
   
   use omp_lib

   implicit none

   double precision,intent(inout) :: E2
   double precision,intent(in) :: BI(NAUXBASD*NVIR)
   double precision,intent(in) :: BJ(NAUXBASD*NVIR,NQVV)
   double precision,intent(in) :: VDX(NAUXBASD,NAUXBAS)
   double precision,intent(in) :: EIG(NACT+NVIR)

   integer,intent(in) :: NACT,NVIR,IACT,LGPU_S,LGPU_E,NQVV,NAUXBASD,NAUXBAS

   double precision :: QVV(NVIR*NQVV*NVIR)
   double precision :: T3(NVIR*NAUXBAS)

   integer :: JACTmod, iQVV, JACT

   double precision :: wtranspose,wb32one,wenergyb32ii,wengcontr
   double precision :: w0,w1,w2,w3,w4,w5

wtranspose = 0.0D00
wb32one = 0.0D00
wenergyb32ii = 0.0D00
wengcontr = 0.0D00

!$omp target data map(to:BI)
      CALL RIMP2_B32_TRANSPOSE(BI,T3,NAUXBAS,NVIR)
#if __cublas==1 || __hipblas==1
      CALL RIMP2_FORM_B32_one(BI,T3,VDX,NAUXBASD,NVIR,NAUXBAS)
#elif __intel_offloading==1
      CALL RIMP2_FORM_B32_one_intel(BI,T3,VDX,NAUXBASD,NVIR,NAUXBAS)      
#endif
   JACT = LGPU_S
   iQVV = LGPU_E - LGPU_S + 1
#if __cublas==1 || __hipblas==1
      CALL RIMP2_QVV_TAIL_cublas(QVV, BI,BJ,NVIR,NAUXBASD,iQVV)
#elif __intel_offloading==1
      CALL RIMP2_QVV_TAIL_mkl(QVV, BI,BJ,NVIR,NAUXBASD,iQVV)
#endif
!$omp end target data

   CALL RIMP2_ENG_TAIL(E2,EIG,QVV,IACT,JACT,NACT,NVIR,iQVV)





END SUBROUTINE RIMP2_TAIL



      SUBROUTINE RIMP2_ENG_HEAD &
                (E2,EIG,QVV,IACT,JACT,NACT,NVIR,iQVV)


      use omp_lib
      use blasHandles
      use iso_c_binding

      implicit none

      ! output
      double precision,intent(inout) :: E2

      ! input
      integer(c_int),intent(in) :: NACT,NVIR,iQVV
      double precision,intent(in) :: EIG(NACT+NVIR)

      ! buffer
      double precision,intent(in) :: QVV(NVIR,iQVV,NVIR)

      ! local
      double precision :: Tijab, Q_t,E2_t,FAC
      integer :: IACT,JACT,IC,IB,IA


          ! ~~~ start ~~~


! compute corr energy
!$omp target teams distribute reduction(+:E2) 
      DO IC=1,iQVV
         E2_t = 0.0D00
!$omp parallel do reduction(+:E2_t) collapse(2)
         DO IB=1,NVIR
            DO IA=1,NVIR            
               Tijab=QVV(IA,IC,IB)/(EIG(IACT+IC-1)+EIG(JACT)-EIG(IA+NACT)-EIG(IB+NACT))
               Q_t=QVV(IA,IC,IB)+QVV(IA,IC,IB)
               E2_t=E2_t + Tijab*(Q_t-QVV(IB,IC,IA))
            ENDDO
         ENDDO
!$omp end parallel do

         FAC=2.0D00
         IF(IACT+IC-1.EQ.JACT) FAC=1.0D00
         E2 = E2 + FAC*E2_t

      ENDDO
!$omp end target teams distribute

      END SUBROUTINE RIMP2_ENG_HEAD


      SUBROUTINE RIMP2_ENG_TAIL &
                (E2,EIG,QVV,IACT,JACT,NACT,NVIR,iQVV)

      use omp_lib
      use blasHandles
      use iso_c_binding

      implicit none

      ! output
      double precision :: E2

      ! input
      integer,intent(in) :: NACT,NVIR,iQVV
      double precision,intent(in) :: EIG(NACT+NVIR)
      
      ! buffer
      double precision,intent(in) :: QVV(NVIR,NVIR,iQVV)

      ! local
      integer :: IACT,JACT,IA,IB,IC
      double precision :: E2_T,Tijab,Q_t,FAC

         
         ! ~~~ start ~~~


! compute corr energy
!$omp target teams distribute reduction(+:E2) 
      DO IC=1,iQVV
         E2_t = 0.0D00
!$omp parallel do reduction(+:E2_t) collapse(2)
         DO IB=1,NVIR
            DO IA=1,NVIR
               Tijab=QVV(IA,IB,IC)/(EIG(IACT)+EIG(JACT+IC-1)-EIG(IA+NACT)-EIG(IB+NACT))
               Q_t=QVV(IA,IB,IC)+QVV(IA,IB,IC)
               E2_t=E2_t + Tijab*(Q_t-QVV(IB,IA,IC))
            ENDDO
         ENDDO
!$omp end parallel do

         FAC=2.0D00
         IF(IACT.EQ.JACT+IC-1) FAC=1.0D00
         E2 = E2 + FAC*E2_t

      ENDDO
!$omp end target teams distribute

      END




! !>*module rimp2grd   *deck RIMP2_X_FORM_MO_I32
! !>
! !>     @brief   form 3-2ERIs
! !>     @author  buu
! !>     @date    sept 8, 2018
! !>
! !>     @detail  form 3-2ERIs using Rys Quadrature method
! !>
!       SUBROUTINE writeRead &
!          (EIG, VEC,NCOR,NACT,NVIR,NBASIS,NAUXBAS,                            &
!           EXAUX,CAUXS,CAUXP,CAUXD,CAUXF,CAUXG,CAUXH,CAUXI,              &
!           KAUXST,KAUXAT,KAUXTY,KAUXNG,KAUXLO,KAUXMI,KAUXMX,NAUXSH,      &
!           EX,CS,CP,CD,CF,CG,CH,CI,                                      &
!           KSTART,KATOM,KTYPE,KNG,KLOC,KMIN,KMAX,NSHELL,                 &
!           PSHELL,DSHELL,FSHELL,GSHELL,PIHELL,DIHELL,FIHELL,GIHELL,      &
!           MXEX,MXEXAUX)
      
!       use gpuGlobalRimp2Eng,only: MXAUXSH,MXAXGTOT, MXSH,MXGTOT, MXATM
!       use gpuGlobalRimp2Eng,only: LX,LY,LZ
!       use gpuGlobalRimp2Eng,only: LOCSPH
!       use gpuGlobalRimp2Eng,only: MaxAUXBFX,MaxATMBFX,MaxAUXANG,MaxATMANG
!       use gpuGlobalRimp2Eng,only: MIXS,TLL





!       COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,                    &
!                       ZAN(MXATM),C(3,MXATM),IAN(MXATM)
!       COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
!       COMMON /OUTPUT/ NPRINT,ITOL,ICUT,NORMF,NORMP,NOPK

!       ! COMMON /AUXBAS/
!       double precision :: &
!          EXAUX(MXEXAUX),CAUXS(MXEXAUX),CAUXP(MXEXAUX), &
!          CAUXD(MXEXAUX),CAUXF(MXEXAUX),CAUXG(MXEXAUX),CAUXH(MXEXAUX),CAUXI(MXEXAUX)
!       integer :: &
!          KAUXST(NAUXSH),KAUXAT(NAUXSH),KAUXTY(NAUXSH), &
!          KAUXNG(NAUXSH),KAUXLO(NAUXSH),KAUXMI(NAUXSH),KAUXMX(NAUXSH),NAUXSH
!       ! COMMON /NSHEL /
!       double precision :: &
!          EX(MXEX),CS(MXEX),CP(MXEX), &
!          CD(MXEX),CF(MXEX),CG(MXEX),CH(MXEX),CI(MXEX)
!       integer :: &
!          KSTART(NSHELL),KATOM(NSHELL),KTYPE(NSHELL), &
!          KNG(NSHELL),KLOC(NSHELL),KMIN(NSHELL),KMAX(NSHELL),NSHELL

!       ! INPUT
!       double precision :: VEC(NBASIS,NBASIS)
!       double precision :: EIG(NBASIS)

!       logical :: MASWRK


!       character(80) :: filename


!       filename = "inp.inp"
!       irecl = 40


! IF(.TRUE.) THEN

! IF(MASWRK) THEN

!       open(unit=911, file=filename,status='new',form="unformatted", &
!            access='direct',recl=irecl,iostat=ierr,action='write')

!       write(911,iostat=ierr,rec=1) NCOR
!       write(911,iostat=ierr,rec=2) NACT
!       write(911,iostat=ierr,rec=3) NVIR
!       write(911,iostat=ierr,rec=4) NBF
!       write(911,iostat=ierr,rec=5) NAUXBAS

!       DO J=1,NBASIS
!          DO I=1,NBASIS
!             IREC = 5 + (J-1)*NBASIS + I
!             write(911,iostat=ierr,rec=IREC) VEC(I,J)
!          ENDDO
!       ENDDO

!       JREC=IREC
!       DO II=1,NBASIS
!          JREC = JREC + 1
!          write(911,iostat=ierr,rec=JREC) EIG(II)
!          write(*,*) EIG(II)
!       ENDDO

! ENDIF

! ELSE

!       open(unit=911, file=filename,status='old',form="unformatted", &
!            access='direct',recl=irecl,iostat=ierr,action='read')

!       read(911,iostat=ierr,rec=1) NCOR
!       read(911,iostat=ierr,rec=2) NACT
!       read(911,iostat=ierr,rec=3) NVIR
!       read(911,iostat=ierr,rec=4) NBF
!       read(911,iostat=ierr,rec=5) NAUXBAS

!       DO J=1,NBASIS
!          DO I=1,NBASIS
!             IREC = 5 + (J-1)*NBASIS + I
!             read(911,iostat=ierr,rec=IREC) VEC(I,J)
!          ENDDO
!       ENDDO

!       JREC=IREC
!       DO II=1,NBASIS
!          JREC = JREC + 1
!          read(911,iostat=ierr,rec=JREC) EIG(II)
!       ENDDO

! ENDIF


!       END SUBROUTINE writeRead


!*MODULE RIMP2GPU  *DECK DEVICEINP
!> @brief      Routine reads in device group
!>
!> @author     Sarom Leang
!>
!> @date    June 7, 2023
!>
!> @details    Routine reads in the $device group
!>
      SUBROUTINE DEVICEINP()

      USE CONSTANTS, ONLY: ZERO
      USE MOD_DEVICE, ONLY: FREEBYTES
      USE COMM_PAR, ONLY: MASWRK
      USE COMM_IOFILE, ONLY: IR,IW

      IMPLICIT NONE

      INTEGER, PARAMETER :: NNAM=1
      DOUBLE PRECISION   :: QNAM(NNAM)
      INTEGER            :: KQNAM(NNAM)
      INTEGER            :: JRET
      DOUBLE PRECISION   :: DEVICE, FREE

      DATA DEVICE/8HDEVICE  /
      DATA QNAM  /8HFREE    /
      DATA KQNAM /3/

!     INITIAL VALUES
      FREE=-1.0D0
      JRET=0

      CALL NAMEIO(IR,JRET,DEVICE,NNAM,QNAM,KQNAM,                       &
      FREE,0,0,0,                                                       &
       0,0,0,0,0,    0,0,0,0,0,   0,0,0,0,0,   0,0,0,0,0,               &
       0,0,0,0,0,    0,0,0,0,0,   0,0,0,0,0,   0,0,0,0,0,               &
       0,0,0,0,0,    0,0,0,0,0,   0,0,0,0,0,   0,0,0,0,0,               &
       0,0,0,0,0,    0)
      IF(JRET.EQ.2) THEN
         IF (MASWRK) WRITE (IW,*) 'ERROR READING $DEVICE GROUP'
         CALL ABRT
      END IF

      IF (FREE > ZERO) THEN
         FREEBYTES=INT(FREE*1073741824.D0)
         IF(MASWRK) WRITE(IW,9000) FREE
         IF(MASWRK) WRITE(IW,9001) FREEBYTES
         FLUSH(IW)
      ENDIF

9000  FORMAT(//5X,'$DEVICE OPTIONS'/5X,15(1H-)/,                        &
             1X,'FREE',4X,'=',F14.2,1X,"GB")
9001  FORMAT(1X,'FREE',4X,'=',I14,1X,"BYTES"/)

      END SUBROUTINE DEVICEINP

#endif





