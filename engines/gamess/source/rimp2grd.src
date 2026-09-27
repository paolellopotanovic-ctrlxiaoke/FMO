#define __bug 0

      MODULE Rimp2_Shared_Data

      use mx_limits, only: mxatm, mxrt, mxgtot, mxao, mxsh

      implicit double precision(a-h,o-z)

!!!!! PARAMETERS
      double precision,parameter ::         &
             ZERO=0.0D00                  , &
             HALF=0.5D00                  , &
             ONE=1.0D00                   , &
             TWO=2.0D00                   , &
             THREE=3.0D00                 , &
             FOUR=4.0D00                  , &
             SQRT3=1.73205080756888D00    , &
             SQRT5=2.23606797749979D00    , &
             SQRT7=2.64575131106459D00    , &
             SQRT9=3.0D00                 , &
             SQRT11=3.3166247903553998D00 , &
             PI252=34.986836655250D00

      integer,parameter :: MXAUXSH=10000,MXAXGTOT=20000

!!!!! ANG MOMENTUM COMPONENT
      integer :: LX(84),LY(84),LZ(84) ! angular momentum component

      DATA LX /   0,   1,   0,   0,   2,   0,   0,   1,   1,   0,       &
                  3,   0,   0,   2,   2,   1,   0,   1,   0,   1,       &
                  4,   0,   0,   3,   3,   1,   0,   1,   0,   2,       &
                  2,   0,   2,   1,   1,                                &
                  5,   0,   0,   4,   4,   1,   0,   1,   0,   3,       &
                  3,   2,   0,   2,   0,   3,   1,   1,   2,   2,       &
                  1,                                                    &
                  6,   0,   0,   5,   5,   1,   0,   1,   0,   4,       &
                  4,   2,   0,   2,   0,   4,   1,   1,   3,   3,       &
                  0,   3,   3,   2,   1,   2,   1,   2/

      DATA LY /   0,   0,   1,   0,   0,   2,   0,   1,   0,   1,       &
                  0,   3,   0,   1,   0,   2,   2,   0,   1,   1,       &
                  0,   4,   0,   1,   0,   3,   3,   0,   1,   2,       &
                  0,   2,   1,   2,   1,                                &
                  0,   5,   0,   1,   0,   4,   4,   0,   1,   2,       &
                  0,   3,   3,   0,   2,   1,   3,   1,   2,   1,       &
                  2,                                                    &
                  0,   6,   0,   1,   0,   5,   5,   0,   1,   2,       &
                  0,   4,   4,   0,   2,   1,   4,   1,   3,   0,       &
                  3,   2,   1,   3,   3,   1,   2,   2/

      DATA LZ /   0,   0,   0,   1,   0,   0,   2,   0,   1,   1,       &
                  0,   0,   3,   0,   1,   0,   1,   2,   2,   1,       &
                  0,   0,   4,   0,   1,   0,   1,   3,   3,   0,       &
                  2,   2,   1,   1,   2,                                &
                  0,   0,   5,   0,   1,   0,   1,   4,   4,   0,       &
                  2,   0,   2,   3,   3,   1,   1,   3,   1,   2,       &
                  2,                                                    &
                  0,   0,   6,   0,   1,   0,   1,   5,   5,   0,       &
                  2,   0,   2,   4,   4,   1,   1,   4,   0,   3,       &
                  3,   1,   2,   1,   2,   3,   3,   2/

!!!!! DEBUGING FLAG
      logical :: BUG ! rank load balancing

!!!!! DO FMO CALCULATIONS
      logical :: DOFMO ! DO FMO
      logical :: DOFMO1R ! DO FMO 1 RANK

!!!!! SPHERICAL BASIS FLAG
      logical :: FlgSPHAUX
      integer,allocatable :: LOCSPH(:)     ! Cartesian->Spherical Trans Matrix


!!!!! MAX ANGULAR MOMENTUM
      integer :: MaxAUXANG ! AUX
      integer :: MaxATMANG ! AO

!!!!! MAX BFS IN A SHELL
      integer :: MaxAUXBFX ! AUX
      integer :: MaxATMBFX ! AO

!!!!! CONTRACTION COEFFICIENT
      double precision,allocatable :: DAtm(:,:,:),DAux(:,:,:)

!!!!! DDI HANDDLES
      integer :: D_V    ! V-MATRIX (NAUXBAS,NAUXBAS)
      integer :: D_B    ! B-MATRIX (NAUXBAS,NVIR,NACT)
      integer :: D_GNAX ! partially transformed NON-SEPARABLE DENSITY (IN MO BASIS)

      integer :: D_VVAA   ! VXX x I32
      integer :: D_AIBJ_TMP


!!!!! SUBROUTINE RIMP2_XLOV
      double precision,allocatable :: XLOV_OMP(:,:),AOO(:,:),VV(:,:),BDV(:,:)
!$omp threadprivate(XLOV_OMP,AOO,VV,BDV)


!!!!! SUBROUTINE RIMP2_ACT_ACT
      double precision,save,allocatable ::                              &
                       PVV_OMP(:,:),POO_OMP(:,:),                       &
                       WOO_OMP(:,:),WVV_OMP(:,:),                       &
                       YDV_OMP(:,:),QVV(:,:),TVV(:,:)
!$omp threadprivate(PVV_OMP,YDV_OMP,POO_OMP,WOO_OMP,WVV_OMP,QVV,TVV)


!!!!! SUBROUTINE RIMP2_DLBMEM_AIBJ
      double precision,save,allocatable :: QOO(:,:),QOV(:,:)
!$omp threadprivate(QOO,QOV)


      double precision,allocatable,dimension(:) ::                      &
            X_WFN,X_XCH,X_GINT,X_FINT,X_SINT,X_IJKLG,X_DAB,             &
            X_CHRG, X_XY,X_XZ,X_YZ,X_X,X_Y,X_Z,X_SJ,X_SK,X_SL,X_GIJKL,  &
            X_GNKL,X_GNM,X_DIJ,X_DKL,X_B00,X_B01,X_B10,X_C00,X_D00,     &
            X_FXX,X_DIJSI,X_DIJSJ,X_DKLSK,X_DKLSL,X_ABV, X_CV,X_RW,     &
            X_AAI,X_AAJ,X_BBK,X_BBL,X_FI,X_FJ,X_FK,X_FL,X_SII,X_SJJ,    &
            X_SKK,X_SLL,X_SIJ,X_SIK,X_SIL,X_SJK,X_SJL,X_SKL
      integer,allocatable,dimension(:,:) :: X_NIJG


!$omp threadprivate(                                                    &
!$omp     X_GINT,X_FINT,X_SINT,X_IJKLG,X_DAB,                           &
!$omp     X_XY,X_XZ,X_YZ,X_X,X_Y,X_Z,X_SJ,X_SK,X_SL,X_GIJKL,            &
!$omp     X_GNKL,X_GNM,X_DIJ,X_DKL,X_B00,X_B01,X_B10,X_C00,X_D00,       &
!$omp     X_FXX,X_DIJSI,X_DIJSJ,X_DKLSK,X_DKLSL,X_ABV, X_CV,X_RW,       &
!$omp     X_AAI,X_AAJ,X_BBK,X_BBL,X_FI,X_FJ,X_FK,X_FL,X_SII,X_SJJ,      &
!$omp     X_SKK,X_SLL,X_SIJ,X_SIK,X_SIL,X_SJK,X_SJL,X_SKL)

      END MODULE



      MODULE RIMP2_TMP
!!!!! SUBROUTINE RIMP2_FORM_VI
      double precision,allocatable ::                                   &
            T1(:), GG(:),TL(:),XX1(:),XX2(:),VNX(:,:)
!$omp threadprivate(T1, GG,TL,XX1,XX2,VNX)


      integer,allocatable,dimension(:) :: NX_INIJG,NX_IIJKLG
!$omp threadprivate(NX_IIJKLG)
      double precision,allocatable,dimension(:) ::                      &
         X_IXCH  ,X_INIJG ,X_ICHRG ,X_IIJKLG,X_IDAB  ,X_IGINT ,X_IFINT ,&
         X_ISINT ,X_IGIJKL,X_IGNKL ,X_IGNM  ,X_IXY   ,X_IXZ   ,X_IYZ   ,&
         X_IX    ,X_IY    ,X_IZ    ,X_ISJ   ,X_ISK   ,X_ISL   ,X_IB00  ,&
         X_IB01  ,X_IB10  ,X_IC00  ,X_ID00  ,X_IF00  ,X_IDIJ  ,X_IDKL  ,&
         X_IDIJSI,X_IDIJSJ,X_IDKLSK,X_IDKLSL,X_IABV  ,X_ICV   ,X_IRW   ,&
         X_IAAI  ,X_IAAJ  ,X_IBBK  ,X_IBBL  ,X_IFI   ,X_IFJ   ,X_IFK   ,&
         X_IFL   ,X_ISII  ,X_ISJJ  ,X_ISKK  ,X_ISLL  ,X_ISIJ  ,X_ISIK  ,&
         X_ISIL  ,X_ISJK  ,X_ISJL  ,X_ISKL

!$omp threadprivate(                                                    &
!$omp                               X_IIJKLG,X_IDAB  ,X_IGINT ,X_IFINT ,&
!$omp    X_ISINT ,X_IGIJKL,X_IGNKL ,X_IGNM  ,X_IXY   ,X_IXZ   ,X_IYZ   ,&
!$omp    X_IX    ,X_IY    ,X_IZ    ,X_ISJ   ,X_ISK   ,X_ISL   ,X_IB00  ,&
!$omp    X_IB01  ,X_IB10  ,X_IC00  ,X_ID00  ,X_IF00  ,X_IDIJ  ,X_IDKL  ,&
!$omp    X_IDIJSI,X_IDIJSJ,X_IDKLSK,X_IDKLSL,X_IABV  ,X_ICV   ,X_IRW   ,&
!$omp    X_IAAI  ,X_IAAJ  ,X_IBBK  ,X_IBBL  ,X_IFI   ,X_IFJ   ,X_IFK   ,&
!$omp    X_IFL   ,X_ISII  ,X_ISJJ  ,X_ISKK  ,X_ISLL  ,X_ISIJ  ,X_ISIK  ,&
!$omp    X_ISIL  ,X_ISJK  ,X_ISJL  ,X_ISKL  )


      END MODULE




#ifndef _OPENMP
      subroutine rimp2grd
         implicit double precision(a-h,o-z)
         COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
         logical GOPARR,DSKWRK,MASWRK

         IF(MASWRK) WRITE(6,*) "ri-mp2 gradient need GMS_OPENMP=true in install.info"
         call ABRT
      end
#else

!>*module rimp2grd   *deck rimp2grd
!>
!>     @brief   read input for ri-mp2 grad
!>     @author  buu
!>     @date    sept 8, 2018
!>
!>     @detail  read input for ri-mp2 grad
!>
      SUBROUTINE rimp2grd
      use omp_lib
      use Rimp2_Shared_Data
      implicit double precision(a-h,o-z)

      logical :: GOPARR,DSKWRK,MASWRK
      logical :: USEDM, GOSMP, OTHAUX, EXTCAB
      integer,parameter :: NORIMP=8,NOAUBF=2
      double precision,parameter :: P12=1.2D00
      integer :: DDI_NP,DDI_ME,DDI_NN,DDI_MY
      integer :: KQRIMP(NORIMP),KQAUXBF(NOAUBF)
      double precision :: QRIMP(NORIMP),QAUXBF(NOAUBF),cabnam(2)

      double precision,allocatable :: DE_DDI(:,:)

      COMMON /AUXBAS/ EXAUX(MXAXGTOT),CAUXS(MXAXGTOT),CAUXP(MXAXGTOT),  &
                      CAUXD(MXAXGTOT),CAUXF(MXAXGTOT),CAUXG(MXAXGTOT),  &
                      CAUXH(MXAXGTOT),CAUXI(MXAXGTOT),                  &
                      KAUXST(MXAUXSH),KAUXAT(MXAUXSH),KAUXTY(MXAUXSH),  &
                      KAUXNG(MXAUXSH),KAUXLO(MXAUXSH),KAUXMI(MXAUXSH),  &
                      KAUXMX(MXAUXSH),NAUXSH
      COMMON /BASSPH/ QMTTOL,ISPHER
      COMMON /ENRGMP/ EMP2,EMP3,EMP4,EMP2A
      COMMON /ENRGYS/ ENUC,EELEC,ETOT,SZ,SZZ,ECORE,ESCF,EERD,E1,E2,     &
                      VEN,VEE,EPOT,EKIN,ESTATE(MXRT),STATN,EDFT(3),EDISP
      COMMON /IOFILE/ IR,IW,IP,IJK,IPK,IDAF,NAV,IODA(950)
      COMMON /MACHIN/ NWDVAR,MAXFM,MAXSM,LIMFM,LIMSM
      COMMON /MP2PAR/ OSPT,CODEMP,SCSPT,TOL,METHOD,NWDMP2,MEMPRI,MPPROP,&
                      NACORE,NBCORE,NOA,NOB,NO,NBF,NOMIT,MOCPHF,MAXITC
      COMMON /NSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),      &
                      CF(MXGTOT),CG(MXGTOT),CH(MXGTOT),CI(MXGTOT),      &
                      KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),   &
                      KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
      COMMON /ONEELC/ EONE,E1A,E1B
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /RESTAR/ TIMLIM,IREST,NREC,INTLOC,IST,JST,KST,LST
      COMMON /RUNOPT/ RUNTYP,EXETYP,NEVALS,NGLEVL,NHLEVL
      COMMON /SCFWFN/ AROHF(3),BROHF(3),PACAVO(6),IACAVO,IUHFNO,ICUHF,  &
                      MVOQ
      COMMON /WFNOPT/ SCFTYP,VBTYP,DFTYPE,TDDFTYP,CITYP,CCTYP,          &
                      MPLEVL,MPCTYP
      COMMON /RIMPFI/ IFILV,IFILT2A,IFILT3A,IFILT2B,IFILT3B
      COMMON /OUTPUT/ NPRINT,ITOL,ICUT,NORMF,NORMP,NOPK
      COMMON /GRAD  / DE(3*MXATM)
!$omp threadprivate(/GRAD  /)
      COMMON /FUNCT / E,EG(3*MXATM)
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,                    &
                      ZAN(MXATM),C(3,MXATM),IAN(MXATM)

      DATA RHF      /8HRHF     /
      DATA UHF      /8HUHF     /
      DATA CHECK    /8HCHECK   /

      DATA RIMP2    /8HRIMP2   /
      DATA AUXBAS   /8HAUXBAS  /
      DATA QRIMP    /8HIAUXBF  ,8HIVMTD   ,8HVTOL    ,8HSTOL    ,       &
                     8HOTHAUX  ,8HGOSMP   ,8HMEMSH   ,8HUSEDM   /
      DATA QAUXBF   /8HCABNAM  ,8HEXTCAB  /
      DATA KQRIMP   /1,1,3,3,0,0,1,0/
      DATA KQAUXBF  /25,0/
      data empty    /8H        /
      DATA G3MP2    /8HG3MP2   /


                          !!! ~ GO ~ !!!


!!!!! MP2 VARS
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


!!!!! INITIALIZE & READ $RIMP2 VARS
      IAUXBF= 0
      IF(ISPHER.EQ.1) IAUXBF=1
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

!!!!! AO & AUX ARE BOTH SPHERICAL/ OR NOT
      IF(IAUXBF.EQ.1 .AND. ISPHER.NE.1) CALL ABRT

!!!!! READ $AUXBAS VARS
      CABNAM(1)=EMPTY
      CABNAM(2)=EMPTY
      EXTCAB=.FALSE.
      CALL NAMEIO(IR,JRET,AUXBAS,NOAUBF,QAUXBF,KQAUXBF,                 &
          CABNAM,EXTCAB,                                                &
          0,0,0,0,0,  0,0,                                              &
          0,0,0,0,0,  0,0,0,0,0,   0,0,0,0,0,                           &
          0,0,0,0,0,  0,0,0,0,0,   0,0,0,0,0,  0,0,0,0,0,               &
          0,0,0,0,0,  0,0,0,0,0,   0,0,0,0,0,  0,0,0,0,0)
      IF(JRET .EQ. 2) CALL ABRT


!!!!! READ AUX BASIS
      call getauxbas(ibas)
! support g3mp2, get the first (ibas=1) aux bas in cabnam list
      if (runtyp.eq.g3mp2) then 
         ibas=1
         if(maswrk) write(*,911) cabnam(ibas),extcab
911                 format(1X,'rimp2 gradient auxiliary basis set',/,&
                           1X,'CABNAM=',A8,1X,'EXTCAB=',L2)
      endif

      CALL RICAUXBAS(CABNAM(ibas),EXTCAB,NAUXCAT)

      NAUXBAS = NAUXCAT


!!!!! SET UP SPHERICAL AUX BASIS SET
      IF(FlgSphAux) THEN
         ALLOCATE(LOCSPH(NAUXSH))
         CALL SPHAUX(LOCSPH,NAUXSPH)
         NAUXBAS = NAUXSPH
      END IF


!!!!! MAX ANG MOMENTUM OF AO AND AUX BASIS FUNCTIONS
      MaxAUXANG=MAXVAL(KAUXTY,1)-1
      MaxAUXBFX=(MaxAUXANG+1)*(MaxAUXANG+2)/2

      MaxATMANG=MAXVAL(KTYPE,1)-1
      MaxATMBFX=(MaxATMANG+1)*(MaxATMANG+2)/2


!!!!! MAX BASIS FUNCTIONS IN AN AO/AUX SHELL
      NANGM=MaxAUXBFX
      IF(MaxATMBFX.gt.MaxAUXBFX) NANGM=MaxATMBFX


!!!!! CONTRACTION COEFFICIENTS
      MaxNGs=MAXVAL(KAUXNG,1)
      IF (.NOT. ALLOCATED(DAux)) ALLOCATE(DAux(NAUXSH,MaxNGs,28))

      MaxNGs=MAXVAL(KNG,1)
      IF (.NOT. ALLOCATED(DAtm)) ALLOCATE(DAtm(NSHELL,MaxNGs,28))

      CALL DensFact()


!!!!! MAX ITERATION FOR Z-VECTOR SOLUTION
      ITERMX=MAXITC

#if __bug==1
      startw = omp_get_wtime()
#endif

!!!!! MAIN RI-MP2 GRADIENT DRIVER
      CALL RIMP2_X_GRADIENT                                             &
          (NANGM,NAUXBASD,NAUXBAS,NAUXSH,                               &
           NCOR,NACT,NVIR,NORB,NBF,                                     &
           STOL,VTOL,IVMTD,OTHAUX,ITERMX)

!!!!! COPY GRADIENT FROM GRAD TO FUNCT COMMON BLOCK
      FORALL(II=1:3*NAT)
         EG(II)=DE(II)
      ENDFORALL

      CALL DFINAL(1)
      IREST = 0

      IF(FlgSphAux)THEN
         DEALLOCATE(LOCSPH,DAux,DAtm)
      ENDIF

#if __bug==1
      endw = omp_get_wtime()
      WRITE(*,'(A40,I5,F10.1)') "SUB RIMP2_X_GRADIENT",ME,endw-startw
#endif

      RETURN

      END !***********************************************************




      SUBROUTINE rimp2_WSA(ISE,Nprocs,Lstart,Lend,LOMP)

      implicit double precision(a-h,o-z)

      integer,intent(in)  :: Nprocs
      integer,intent(in)  :: Lstart, Lend

      logical,intent(out) :: LOMP
      integer,intent(out) :: ISE(0:Nprocs-1,2)


      Nloops=Lend-Lstart+1
      IF(Nloops.lt.Nprocs) THEN
         LOMP=.false.
         RETURN
      ELSE
         LOMP=.true.
      ENDIF

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

      END !**********************************************************




!>*module rimp2grd   *deck RIMP2_X_GRADIENT
!>
!>     @brief   ri-mp2 grad driver
!>     @author  buu
!>     @date    sept 8, 2018
!>
!>     @detail  ri-mp2 grad driver
!>

      SUBROUTINE RIMP2_X_GRADIENT                                       &
                (NANGM,NAUXBASD,NAUXBAS,NAUXSH,                         &
                 NCOR,NACT,NVIR,NORB,NBF,                               &
                 STOL,VTOL,IVMTD,OTHAUX,ITERMX)

      use Rimp2_Shared_Data,only: D_V,D_B,D_GNAX,MXATM,BUG

      implicit double precision(a-h,o-z)

      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,                    &
                      ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /IOFILE/ IR,IW,IP,IJK,IPK,IDAF,NAV,IODA(950)
      COMMON /GRAD  / DE(3,MXATM)
!$omp threadprivate(/GRAD  /)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      logical :: GOPARR,DSKWRK,MASWRK
      logical :: ULRange
      logical :: OTHAUX


      double precision,allocatable ::                                   &
             PSCF(:,:),PAO(:),PMP2(:,:),POO(:,:),PVV(:,:),              &
             XLOV(:,:),ZOV(:,:),                                        &
             WAO(:),WMP2(:,:),WAA(:,:),WOO(:,:),WOV(:,:),WVV(:,:),      &
             EIG(:),VEC(:,:),VDX(:,:),VXX(:,:),                         &
             VI(:,:,:),GNNL(:,:,:),GNAL(:,:,:),                         &
             gXX(:,:),gXD(:,:),gDX(:,:), DE_DDI(:,:)






!///////////////////////////////////////////////////////////////////////
!/////                  INITIALIZATION                      ////////////
!///////////////////////////////////////////////////////////////////////


!!!!! DEBUGING FLAG
      BUG=.false.


!!!!! SYNCHRONIZE DDI PROCESSES
      CALL DDI_SYNC(5123)


!!!!! SOME PARAMETERS
      NOCC=NCOR+NACT
      NORB=NCOR+NACT+NVIR

      NXO=NAUXBASD*NOCC
      NXV=NAUXBASD*NVIR
      NXN=NAUXBASD*NORB


!///////////////////////////////////////////////////////////////////////
!/////                     FORM THE MATRIX B                ////////////
!///////////////////////////////////////////////////////////////////////


!!!!! FORM VDX
#if __bug==1
      startw = omp_get_wtime()
#endif  

      ALLOCATE(VDX(NAUXBAS,NAUXBAS))
      CALL RIMP2_FORM_VDX(VDX, STOL,VTOL,NAUXBAS,NAUXBASD,IVMTD,OTHAUX)
      CALL DDI_SYNC(5127)

#if __bug==1      
      endw = omp_get_wtime()
      WRITE(*,'(A40,I5,F10.1)') "TIME SUB RIMP2_FORM_VDX ",ME, endw-startw
#endif

!!!!! READ AND SYM TRANSFORM MO VECTORS
      ALLOCATE(VEC(NBF,NBF))
      CALL DAREAD(IDAF,IODA,VEC,NBF*NBF,15,0)
      CALL RIMP2_SYM_TRANS_MO(VEC,NBF,NCOR,NORB-NCOR)


!!!!! DDI_CREATE :: B32(NAUXBAS*NBF,NBF)
      CALL DDI_CREATE(NAUXBAS*NBF, NBF, D_B)


!!!!! FORM I32 AND DO MO-TRANSFORMATION
#if __bug==1
      startw = omp_get_wtime()
#endif

      CALL RIMP2_X_FORM_MO_I32(VEC,NBF,NACT,NVIR,NAUXBAS)
      CALL DDI_SYNC(5131)

#if __bug==1      
      endw = omp_get_wtime()
      WRITE(*,'(A40,I5,F10.1)') "TIME SUB RIMP2_X_FORM_MO_I32",ME, endw-startw
#endif

!!!!! COMBINE I32 WITH THE DECOMPOSED INVERSED V-MATRIX
!!!!! THE MATRIX B IS FORMED AND STORED IN DDI ARRAY (D_B)
#if __bug==1
      startw = omp_get_wtime()
#endif

      CALL RIMP2_COMBINE_VDX_I32(VDX, NAUXBAS,NAUXBASD,NBF)

#if __bug==1
      endw = omp_get_wtime()
      WRITE(*,'(A40,I5,F10.1)') "TIME SUB RIMP2_COMBINE_VDX_I32",ME, endw-startw
#endif



!///////////////////////////////////////////////////////////////////////
!/////                FORM (PART OF) DENSITY MATRICES        ///////////
!///////////////////////////////////////////////////////////////////////


!!!!! READ MO ENERGY
      ALLOCATE(EIG(NBF))
      CALL DAREAD(IDAF,IODA,EIG,NBF,17,0)


!!!!! ALLOCATE ARRAYS FOR DENSITY
!!!!! ENERGY-WEIGHTED DENSITY, LAGRANGIAN
      ALLOCATE(POO(NOCC,NOCC))
      ALLOCATE(PVV(NVIR,NVIR))
      ALLOCATE(WOO(NOCC,NOCC))
      ALLOCATE(WOV(NOCC,NVIR))
      ALLOCATE(WVV(NVIR,NVIR))
      ALLOCATE(XLOV(NOCC,NVIR))
      ALLOCATE(gDX(NAUXBASD,NAUXBAS))
      XLOV=0.0D00


!!!!! CREATE DDI ARRAY FOR 3-INDEX DENSITY
      CALL DDI_CREATE(NBF*NACT,NAUXBAS,D_GNAX)


!!!!! ACT-ACT OCCUPIED LOOP
#if __bug==1
      startw = omp_get_wtime()
#endif

      CALL RIMP2_ACT_ACT                                            &
          (gDX, POO,PVV, WOO,WOV,WVV,                               &
           XLOV,VDX,EIG,VEC, NANGM,NAUXBASD,NAUXBAS,                &
           NCOR,NACT,NOCC,NVIR,NORB,NBF)

#if __bug==1
      endw = omp_get_wtime()
      WRITE(*,'(A40,I5,F10.1)') "TIME SUB RIMP2_ACT_ACT",ME,endw-startw
#endif

!!!!! ON-NODE TEMPORARY GRADIENT ARRAY
      ALLOCATE(DE_DDI(3,NAT))
      DE_DDI=0.0D00




!///////////////////////////////////////////////////////////////////////
!/////           2-INDEX 2-PARTICLE GRADIENT CONTRIBUTION      /////////
!///////////////////////////////////////////////////////////////////////


!!!!! DDI WORK DISTRIBUTION :: NAUXSH LOOP
      CALL RIMP2_DDIAUXSH(LddiAuxShStart,LddiAuxShEnd, ME,NPROC,NAUXBAS)
      CALL RIMP2_AUXLOC(LddiAuxStart, LddiAuxShStart,.TRUE.)
      CALL RIMP2_AUXLOC(LddiAuxEnd, LddiAuxShEnd,.FALSE.)


!!!!! COMBINE 2-INDEX DENSITY AND THE V-MATRIX DERIVATIVE
#if __bug==1
      startw = omp_get_wtime()
#endif

      CALL RIMP2_D22GRD                                                 &
          (DE_DDI, gDX,VDX, NAUXBASD,NAUXBAS,                           &
           LddiAuxShStart,LddiAuxShEnd,                                 &
           LddiAuxStart,LddiAuxEnd)
      DEALLOCATE(gDX,VDX)

#if __bug==1
      endw = omp_get_wtime()
      WRITE(*,'(A40,I5,F10.1)') "TIME SUB RIMP2_D22GRD ",ME,endw-startw
#endif



!///////////////////////////////////////////////////////////////////////
!/////           3-INDEX 2-PARTICLE GRADIENT CONTRIBUTION      /////////
!///////////////////////////////////////////////////////////////////////


!!!!! ALLOCATE ON NODE ARRAYS
      startw = omp_get_wtime()
      ALLOCATE(GNNL(NBF,NBF,LddiAuxStart:LddiAuxEnd))
      ALLOCATE(GNAL(NBF,NACT,LddiAuxStart:LddiAuxEnd))


!!!!! CREATE DDI ARRAY
      CALL DDI_GET(D_GNAX, 1,NBF*NACT, LddiAuxStart,LddiAuxEnd, GNAL)


!!!!! FORM 3-INDEX DENSITY
      CALL RIMP2_32_DENS                                                &
          (GNNL, GNAL, VEC,                                             &
           LddiAuxShStart,LddiAuxShEnd,                                 &
           LddiAuxStart,LddiAuxEnd,                                     &
           NCOR, NACT,NBF,NAUXBAS,NAUXSH)
      DEALLOCATE(GNAL)

#if __bug==1
      endw = omp_get_wtime()
      WRITE(*,'(A40,I5,F10.1)') "TIME SUB RIMP2_32_DENS ",ME,endw-startw
#endif

!!!!! COMBINE 3-INDEX DENSITY AND 32ERI DERIVATIVE
      startw = omp_get_wtime()
      CALL RIMP2_D32GRD                                                 &
          (DE_DDI, GNNL,                                                &
           LddiAuxShStart,LddiAuxShEnd,                                 &
           LddiAuxStart,LddiAuxEnd,                                     &
           NBF, NCOR,NACT,NOCC,NVIR,NBF,NAUXBAS)
      DEALLOCATE(GNNL)

#if __bug==1      
      endw = omp_get_wtime()
      WRITE(*,'(A40,I5,F10.1)') "TIME SUB RIMP2_D32GRD ",ME,endw-startw
#endif

!!!!! CLEAN UP
      CALL DDI_DESTROY(D_GNAX)




!///////////////////////////////////////////////////////////////////////
!/////              COMPLETE DENSITY MATRICES                 //////////
!///////////////////////////////////////////////////////////////////////


!!!!! DENSITY REDUCTION FROM ALL DDI PROCESSES
      startw = omp_get_wtime()
      CALL DDI_GSUMF(6000,POO,NOCC*NOCC)
      CALL DDI_GSUMF(6001,PVV,NVIR*NVIR)
      CALL DDI_GSUMF(6003,WOO,NOCC*NOCC)
      CALL DDI_GSUMF(6005,WVV,NVIR*NVIR)

#if __bug==1
      endw = omp_get_wtime()
      WRITE(*,'(A40,I5,F10.1)') "SUB DDI_GSUMF(POO,PVV,WOO,WVV)",ME,endw-startw
#endif

!!!!! FORM WOO[II],WVV[II]
      startw = omp_get_wtime()
      CALL RIMP2_XLOV_VXX_WOO_WVV                                       &
          (XLOV,WOO,WVV, WOV,POO,PVV,EIG,                               &
           NAUXBASD,NAUXBAS, NCOR,NACT,NVIR,NORB)

#if __bug==1
      endw = omp_get_wtime()
      WRITE(*,'(A40,I5,F10.1)') "SUB RIMP2_XLOV_VXX_WOO_WVV ",ME,endw-startw
#endif

!!!!! REDUCING WOV FROM ALL DDI PROCESSES
      CALL DDI_GSUMF(6004,WOV,NOCC*NVIR)


!!!!! FORM XLOV
      startw = omp_get_wtime()
      CALL RIMP2_XLOV(XLOV, POO,PVV,VEC,EIG,NAUXBASD,NCOR,NACT,NVIR,NORB,NBF)

#if __bug==1
      endw = omp_get_wtime()
      WRITE(*,'(A40,I5,F10.1)') "SUB RIMP2_XLOV ",ME,endw-startw
#endif

!!!!! REDUCING XLOV FROM ALL DDI PROCESSES
      CALL DDI_GSUMF(6007,XLOV,NOCC*NVIR)


!!!!! SOLVE FOR Z-VECTOR (ZOV)
      startw = omp_get_wtime()
      ALLOCATE(ZOV(NOCC,NVIR))
      CALL RIMP2_ZVECTOR                                                &
          (ZOV, XLOV,EIG,ITERMX,NAUXBASD,NOCC,NVIR,NORB,NBF)

#if __bug==1
      endw = omp_get_wtime()
      WRITE(*,'(A40,I5,F10.1)') "SUB RIMP2_ZVECTOR ",ME,endw-startw
#endif

!!!!! FORM WOO AND WOV
      startw = omp_get_wtime()
      CALL RIMP2_WOO_WOV                                                &
          (WOO,WOV,POO,ZOV,PVV,VEC,EIG,NAUXBASD,NCOR,NACT,NVIR,NORB,NBF)

#if __bug==1
      endw = omp_get_wtime()
      WRITE(*,'(A40,I5,F10.1)') "SUB RIMP2_WOO_WOV ",ME,endw-startw
#endif

!!!!! CLEAN UP (DESTROY THE MATRIX B)
      CALL DDI_DESTROY(D_B)


!!!!! TRANSFORMING ENERGY-WEIGHTED DENSITY MATRIX INTO AO BASIS
!!!!! WMP2 <- WVV,WOV,WOO
      startw = omp_get_wtime()
      ALLOCATE(WMP2(NBF,NBF))
      CALL RIMP2_WMP2                                                   &
          (WMP2, WVV,WOV,WOO,VEC,NAUXBASD,NCOR,NACT,NVIR,NORB,NBF)
      DEALLOCATE(WVV,WOV,WOO)

#if __bug==1
      endw = omp_get_wtime()
      WRITE(*,'(A40,I5,F10.1)') "SUB RIMP2_WMP2 ",ME,endw-startw
#endif

!!!!! TRANSFORMING 1-PARTICLE DENSITY INTO THE AO BASIS
!!!!! PMP2 <- POO,ZOV,PVV
      startw = omp_get_wtime()
      ALLOCATE(PMP2(NBF,NBF))
      CALL RIMP2_PMP2(PMP2, POO,ZOV,PVV,VEC,NCOR,NACT,NVIR,NORB,NBF)
      DEALLOCATE(POO,ZOV,PVV)

#if __bug==1
      endw = omp_get_wtime()
      WRITE(*,'(A40,I5,F10.1)') "SUB RIMP2_PMP2 ",ME,endw-startw
#endif



!///////////////////////////////////////////////////////////////////////
!/////              1-ELECTRON GRADIENT CONTRIBUTION          //////////
!///////////////////////////////////////////////////////////////////////


!!!!! FORM PSCF, AND APPLY SYMMETRY FOR PMP2
      startw = omp_get_wtime()
      ALLOCATE(PSCF(NBF,NBF))
      CALL RIMP2_WAO_PAO(PMP2,PSCF,WMP2,EIG,VEC, NCOR,NACT,NVIR,NORB,NBF)
      DEALLOCATE(WMP2)

#if __bug==1
      endw = omp_get_wtime()
      WRITE(*,'(A40,I5,F10.1)') "SUB RIMP2_WAO_PAO ",ME,endw-startw
#endif

!!!!! 1-ELECTRON GRADIENT CONTRIBUTION
      startw = omp_get_wtime()
      CALL STVDER

#if __bug==1
      endw = omp_get_wtime()
      WRITE(*,'(A40,I5,F10.1)') "SUB STVDER ",ME,endw-startw
#endif

!!!!! GRADIENT REDUCTION
      CALL DDI_GSUMF(6013,DE_DDI,3*NAT)

      FORALL (II=1:3,JJ=1:NAT)
         DE(II,JJ)=DE(II,JJ)+DE_DDI(II,JJ)
      ENDFORALL




!///////////////////////////////////////////////////////////////////////
!/////           4-INDEX 2-ELECTRON GRADIENT CONTRIBUTION     //////////
!///////////////////////////////////////////////////////////////////////


      startw = omp_get_wtime()
      CALL RIMP2_JKDER

#if __bug==1
      endw = omp_get_wtime()
      WRITE(*,'(A40,I5,F10.1)') "SUB RIMP2_JKDER: ",ME, endw-startw
#endif

      END !!! done :P





      SUBROUTINE DensFact()

      use Rimp2_Shared_Data

      implicit double precision(a-h,o-z)

      COMMON /AUXBAS/ EXAUX(MXAXGTOT),CAUXS(MXAXGTOT),CAUXP(MXAXGTOT),  &
                      CAUXD(MXAXGTOT),CAUXF(MXAXGTOT),CAUXG(MXAXGTOT),  &
                      CAUXH(MXAXGTOT),CAUXI(MXAXGTOT),                  &
                      KAUXST(MXAUXSH),KAUXAT(MXAUXSH),KAUXTY(MXAUXSH),  &
                      KAUXNG(MXAUXSH),KAUXLO(MXAUXSH),KAUXMI(MXAUXSH),  &
                      KAUXMX(MXAUXSH),NAUXSH

      COMMON /NSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),      &
                      CF(MXGTOT),CG(MXGTOT),CH(MXGTOT),CI(MXGTOT),      &
                      KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),   &
                      KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL

      COMMON /OUTPUT/ NPRINT,ITOL,ICUT,NORMF,NORMP,NOPK

      logical NORM




      NORM=NORMF.NE.1 .OR. NORMP.NE.1

      DAux=0.0D00
      DO II=1,NAUXSH
         NGI   = KAUXNG(II)
         MINI  = KAUXMI(II)
         MAXI  = KAUXMX(II)
         I1    = KAUXST(II)
         CALL AuxDensFact(ii,EXAUX(I1),CAUXS(I1),CAUXP(I1),             &
              CAUXD(I1),CAUXF(I1),CAUXG(I1),CAUXH(I1),CAUXI(I1),        &
              NGI,MINI,MAXI,NORM)
      ENDDO

      DAtm=0.0D00
      DO ISH=1,NSHELL
         NGI   = KNG(ISH)
         MINI  = KMIN(ISH)
         MAXI  = KMAX(ISH)
         I1    = KSTART(ISH)

         CALL AtmDensFact(ish,EX(I1),CS(I1),CP(I1),CD(I1),CF(I1),       &
              CG(I1),CH(I1),CI(I1),NGI,MINI,MAXI,NORM)
      ENDDO

      END !***********************************************************



      SUBROUTINE AUXDENSFACT(II,EXI,CIS,CIP,CID,CIF,CIG,CIH,CII,        &
                 NGI,MINI,MAXI,NORM)

      use Rimp2_Shared_Data
      implicit double precision(a-h,o-z)

      logical :: norm
      double precision :: EXI(NGI),CIS(NGI),CIP(NGI),CID(NGI),CIF(NGI), &
                          CIG(NGI),CIH(NGI),CII(NGI)




      DO IG = 1, NGI
         AI  = EXI(IG)
         AIINV = ONE/AI

         ICNT = 0
         DO I = MINI, MAXI
            IF(I .EQ. 1) THEN
               DUM1 = CIS(IG)
            ELSE IF(I .LE.  4) THEN
               DUM1 = CIP(IG)
            ELSE IF(I .LE. 10) THEN
               DUM1 = CID(IG)
               IF(I .GE. 8 .AND. norm) DUM1 = DUM1*SQRT3
            ELSE IF(I .LE. 20) THEN
               DUM1 = CIF(IG)
               IF(I .GE. 14 .AND. norm) THEN
                  IF(I .LE. 19) THEN
                     DUM1 = DUM1*SQRT5
                  ELSE
                     DUM1 = DUM1*SQRT3*SQRT5
                  END IF
               END IF
            ELSE IF(I .LE. 35) THEN
               DUM1 = CIG(IG)
               IF(I .GE. 24 .AND. norm) THEN
                  IF(I .LE. 29) THEN
                     DUM1 = DUM1*SQRT7
                  ELSE IF(I .LE. 32) THEN
                     DUM1 = DUM1*SQRT7*SQRT5/SQRT3
                  ELSE
                     DUM1 = DUM1*SQRT7*SQRT5
                  END IF
               END IF
            ELSE IF(I .LE. 56) THEN
               DUM1 = CIH(IG)
               IF(I .GE. 39 .AND. norm) THEN
                  IF(I .LE. 44) THEN
                     DUM1 = DUM1*SQRT9
                  ELSE IF(I .LE. 50) THEN
                     DUM1 = DUM1*SQRT7*SQRT3
                  ELSE IF(I .LE. 53) THEN
                     DUM1 = DUM1*SQRT9*SQRT7
                  ELSE
                     DUM1 = DUM1*SQRT7*SQRT5*SQRT3
                  ENDIF
               ENDIF
            ELSE
               DUM1 = CII(IG)
               IF(I .GE. 60 .AND. norm) THEN
                  IF(I .LE. 65) THEN
                     DUM1 = DUM1*SQRT11
                  ELSE IF(I .LE. 71) THEN
                     DUM1 = DUM1*SQRT11*SQRT3
                  ELSE IF(I .LE. 74) THEN
                     DUM1 = DUM1*SQRT11*THREE
                  ELSE IF(I .LE. 77) THEN
                     DUM1 = DUM1*SQRT11*SQRT7*SQRT3/SQRT5
                  ELSE IF(I .LE. 83) THEN
                     DUM1 = DUM1*SQRT11*SQRT7*SQRT3
                  ELSE
                     DUM1 = DUM1*SQRT11*SQRT7*SQRT5
                  ENDIF
               END IF
            ENDIF
            ICNT = ICNT+1

            DAux(ii,IG,ICNT)=DUM1*AIINV

         END DO

      ENDDO

      END !***********************************************************




      SUBROUTINE ATMDENSFACT(ISH,EXI,CIS,CIP,CID,CIF,CIG,CIH,CII,       &
                 NGI,MINI,MAXI,NORMI)

      use Rimp2_Shared_Data

      implicit double precision(a-h,o-z)

      logical NORMI

      double precision EXI(NGI),CIS(NGI),CIP(NGI),CID(NGI),CIF(NGI),    &
                       CIG(NGI),CIH(NGI),CII(NGI)





      DO IG = 1, NGI

        icnt=0
         DO I = MINI, MAXI
            IF(I .EQ. 1) THEN
               DUM1 = CIS(IG)
            ELSE IF(I .LE.  4) THEN
               DUM1 = CIP(IG)
            ELSE IF(I .LE. 10) THEN
               DUM1 = CID(IG)
               IF(I .GE. 8 .AND. NORMI) DUM1 = DUM1*SQRT3
            ELSE IF(I .LE. 20) THEN
               DUM1 = CIF(IG)
               IF(I .GE. 14 .AND. NORMI) THEN
                  IF(I .LE. 19) THEN
                     DUM1 = DUM1*SQRT5
                  ELSE
                     DUM1 = DUM1*SQRT3*SQRT5
                  END IF
               END IF
            ELSE IF(I .LE. 35) THEN
               DUM1 = CIG(IG)
               IF(I .GE. 24 .AND. NORMI) THEN
                  IF(I .LE. 29) THEN
                     DUM1 = DUM1*SQRT7
                  ELSE IF(I .LE. 32) THEN
                     DUM1 = DUM1*SQRT7*SQRT5/SQRT3
                  ELSE
                     DUM1 = DUM1*SQRT7*SQRT5
                  END IF
               END IF
            ELSE IF(I .LE. 56) THEN
               DUM1 = CIH(IG)
               IF(I .GE. 39 .AND. NORMI) THEN
                  IF(I .LE. 44) THEN
                     DUM1 = DUM1*SQRT9
                  ELSE IF(I .LE. 50) THEN
                     DUM1 = DUM1*SQRT7*SQRT3
                  ELSE IF(I .LE. 53) THEN
                     DUM1 = DUM1*SQRT9*SQRT7
                  ELSE
                     DUM1 = DUM1*SQRT7*SQRT5*SQRT3
                  ENDIF
               ENDIF
            ELSE
               DUM1 = CII(IG)
               IF(I .GE. 60 .AND. NORMI) THEN
                  IF(I .LE. 65) THEN
                     DUM1 = DUM1*SQRT11
                  ELSE IF(I .LE. 71) THEN
                     DUM1 = DUM1*SQRT11*SQRT3
                  ELSE IF(I .LE. 74) THEN
                     DUM1 = DUM1*SQRT11*THREE
                  ELSE IF(I .LE. 77) THEN
                     DUM1 = DUM1*SQRT11*SQRT7*SQRT3/SQRT5
                  ELSE IF(I .LE. 83) THEN
                     DUM1 = DUM1*SQRT11*SQRT7*SQRT3
                  ELSE
                     DUM1 = DUM1*SQRT11*SQRT7*SQRT5
                  ENDIF
               END IF
            END IF

            icnt=icnt+1
            DAtm(ish,ig,icnt)=dum1

         END DO
      ENDDO

      END !***********************************************************



!>*module rimp2grd   *deck RIMP2_XLOV_VXX_WOO_WVV
!>
!>     @brief   form energy-weighted DM
!>     @author  buu
!>     @date    sept 8, 2018
!>
!>     @detail  form energy-weighted DM
!>

      SUBROUTINE RIMP2_XLOV_VXX_WOO_WVV                                 &
                (XLOV,WOO,WVV, WOV,POO,PVV,EIG,                         &
                 NAUXBASD,NAUXBAS, NCOR,NACT,NVIR,NORB)

      use omp_lib

      implicit double precision(a-h,o-z)

      double precision :: XLOV(NCOR+NACT,NVIR)
      double precision :: WOO(NCOR+NACT,NCOR+NACT)
      double precision :: WVV(NVIR,NVIR)

      double precision :: WOV(NCOR+NACT,NVIR)
      double precision :: POO(NCOR+NACT,NCOR+NACT)
      double precision :: PVV(NVIR,NVIR)
      double precision :: EIG(NORB)




!!!!! PARAMETERS
      NOCC=NCOR+NACT

!!!!! OMP THREADS
      Nthreads_env=omp_get_max_threads()
      Nthreads=Nthreads_env

#if __bug==1
      WRITE(*,'(A40,2I5)') "NUM_THREAD SUB RIMP2_XLOV_VXX_WOO_WVV", Nthreads, Nthreads_env
#endif

!$OMP PARALLEL NUM_THREADS(Nthreads)                                    &
!$omp default(none)                                                     &
!$omp shared(NOCC,NVIR,NAUXBAS,NAUXBASD)                                &
!$omp private(IX,II,JJ,IA,IB)                                           &
!$omp private(eab,eb,TMP)						&
!$omp shared(XLOV,WOV,WOO,POO,WVV,PVV,EIG)


!!!!! XLOV <- 0.5*WOV
!$omp workshare
      XLOV=XLOV-0.5D00*WOV
!$omp end workshare nowait


!!!!! WOO(i,j) <- POO(i,j) * (e(i)+e(j))
!$omp workshare
      FORALL(II=1:NOCC, JJ=1:NOCC)
         WOO(II,JJ)=WOO(II,JJ) - POO(II,JJ)*(EIG(II)+EIG(JJ))
      ENDFORALL
!$omp end workshare nowait


!!!!! WVV(a,b) <- PVV(a,b) * (e(a)+e(b))
!$omp do schedule(DYNAMIC)
      DO IB=1,NVIR
         eb=EIG(IB+NOCC)
         DO IA=1,NVIR
            eab=eb+EIG(IA+NOCC)
            WVV(IA,IB)=WVV(IA,IB) + PVV(IA,IB)*eab
         ENDDO
      ENDDO
!$omp end do

!$OMP END PARALLEL

      END !***********************************************************




!>*module rimp2grd   *deck RIMP2_XLOV
!>
!>     @brief   form Lagrangian
!>     @author  buu
!>     @date    sept 8, 2018
!>
!>     @detail  form Lagrangian
!>

      SUBROUTINE RIMP2_XLOV                                             &
                (XLOV, POO,PVV,VEC,EIG,                                 &
                 NAUXBASD,NCOR,NACT,NVIR,NORB,NBF)

      use omp_lib
      use Rimp2_Shared_Data,only : D_B
      use Rimp2_Shared_Data,only : XLOV_OMP,AOO,VV,BDV

      implicit double precision(a-h,o-z)

      COMMON /IOFILE/ IR,IW,IP,IJK,IPK,IDAF,NAV,IODA(950)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      logical :: GOPARR,DSKWRK,MASWRK


!!!!! OUTPUT
      double precision :: XLOV(NCOR+NACT,NVIR)


!!!!! INPUT
      double precision :: POO(NCOR+NACT,NCOR+NACT)
      double precision :: PVV(NVIR,NVIR)
      double precision :: VEC(NBF,NORB)
      double precision :: EIG(NBF)

!!!!! LOCAL DATA
      double precision,allocatable ::                                   &
         BDOO(:,:,:),BDOV(:,:,:),BDVV(:,:,:),VVV(:,:,:),BDVO(:,:,:),    &
         BDO(:,:),BI(:,:),BB(:,:)

      integer,allocatable :: ddiWSA(:,:),I_DLB_INDEX(:)

      logical :: lomp

!!!!! PARAMETERS
      NOCC=NCOR+NACT

      NDO=NAUXBASD*NOCC
      NDV=NAUXBASD*NVIR


!!!!! buffer zone for DBL
      N_DLB_BUFF=1 ! DLB buffer


!!!!! DDI WORK DISTRIBUTING FOR NOCC LOOP
      ALLOCATE(ddiWSA(0:NPROC-1,2))
      CALL rimp2_WSA(ddiWSA,NPROC,1,NOCC,LOMP)
      LddiOccStart=ddiWSA(ME,1)
      LddiOccEnd=ddiWSA(ME,2)


!!!!! DDI WORK DISTRIBUTING FOR NVIR LOOP
      CALL rimp2_WSA(ddiWSA,NPROC,1,NVIR,LOMP)
      LddiVirtStart=ddiWSA(ME,1)
      LddiVirtEnd=ddiWSA(ME,2)


!!!!! NUMBER OF THREADS FROM ENV VAR `OMP_NUM_THREADS`
      Nthreads_env=omp_get_max_threads()


!!!!! MAX NUMBER OF THREADS THAT CAN KEEP threadprivate DATA CONSISTENCY
      Nthreads=MIN(Nthreads_env,NOCC,NVIR,LddiVirtEnd-LddiVirtStart-N_DLB_BUFF)


!!!!! PRINT INFO ABOUT DATA CONSISTENCY
      IF(MASWRK) THEN
         IF(Nthreads_env.gt.Nthreads) THEN
            WRITE(IW,100) Nthreads_env, Nthreads
 100  FORMAT(/1X,'TO GUARANTEE DATA CONSISTENCY OF threadprivate ',/,   &
             1X,'IN SUBROUTINE RIMP2_XLOV THE NUM THREADS USED',/       &
             1X,'IN THIS REGION IS RESET FROM',1X,2I5,/)

             WRITE(IW,120)
 120  FORMAT(/1X,'THIS USUALLY HAPPENS WHEN A LARGE NUM THREADS/RANK',/,&
             1X,'IS USED FOR A SMALL MOLECULE. PLEASE ALSO CHECK IF',/  &
             1X,'YOU SET NUMBER OF THREADS/RANK = NUMBER OF',/   &
             1X,'THREADS/SOCKET OR NUMA NODE. IT IS USUALLY NOT',/ &
             1X,'EFFICIENT TO SET THREADS/RANK = THREADS/NODE'/)
         ENDIF
      ENDIF


!!!!! GET MATRIX B
      ALLOCATE(BDOO(NAUXBASD,NOCC,NOCC))
      ALLOCATE(BDOV(NAUXBASD,NOCC,LddiVirtStart:LddiVirtEnd))
      CALL DDI_GET(D_B, 1,NDO, 1,NOCC, BDOO)
      CALL DDI_GET(D_B, 1,NDO, LddiVirtStart+NOCC,LddiVirtEnd+NOCC, BDOV)


!!!!! TURN OFF DYNAMIC THREADS TO KEEP threadprivate DATA CONSISTENT
      CALL OMP_SET_DYNAMIC(.FALSE.)



!!!!! THREADPRIVATE ARRAYS
!$OMP PARALLEL NUM_THREADS(Nthreads)
      ALLOCATE(XLOV_OMP(NOCC,NVIR))
      ALLOCATE(AOO(NOCC,NOCC))
      XLOV_OMP=0.0D00
      AOO=0.0D00
!$OMP END PARALLEL



!///////////////////////////////////////////////////////////////////////
!////                                                            ///////
!////                  STEP 1: XLOV <- AOO*POO                   ///////
!////                                                            ///////
!//// - each rank is statically assigned a chunk of work         ///////
!////   >> (LddiVirtStart, LddiVirtEnd-N_DLB_BUFF)               ///////
!////                                                            ///////
!////                                                            ///////
!//// - the leaving work (>> N_DLB_BUFF on each rank)            ///////
!////   are dynamically distributed among ranks.                 ///////
!////                                                            ///////
!///////////////////////////////////////////////////////////////////////


!//// XLOV <- AOO*POO: static chunks

#if __bug==1
      startw = omp_get_wtime()
#endif

!$OMP PARALLEL NUM_THREADS(Nthreads)                                    &
!$omp default(none)                                                     &
!$omp shared(LddiVirtStart,LddiVirtEnd, N_DLB_BUFF)                     &
!$omp shared(POO,BDOO,BDOV, NVIR,NOCC,NAUXBASD)                         &
!$omp private(IA)

!$omp do schedule(DYNAMIC)
      ! DO IA=1,NVIR
      DO IA=LddiVirtStart,LddiVirtEnd-N_DLB_BUFF
         CALL RIMP2_XLOV__AOOPOO                                        &
             (XLOV_OMP, IA, POO,AOO, BDOO,BDOV(1,1,IA),                 &
              NAUXBASD,NOCC,NVIR)
      ENDDO
!$omp end do
!$OMP END PARALLEL

#if __bug==1
      endw_static = omp_get_wtime()
      WRITE(*,'(A40,I5,F10.1)') "XLOV_STATIC: AOO_POO",ME,endw_static-startw
#endif

!!!!! XLOV <- AOO*POO: DLB
      startw_dyn = omp_get_wtime()


!!!!! DBL buffer array
      ALLOCATE(I_DLB_INDEX(NVIR))
      CALL RIMP2_IDBLINDEX(I_DLB_INDEX,NVIR,N_DLB_BUFF,ddiWSA)

      ALLOCATE(BDO(NAUXBASD,NOCC))
      BDO = 0.0D00

      MINE = -1
      NEXT = -1
      DO IA=1,NVIR

         IF(I_DLB_INDEX(IA).EQ.0) CYCLE

!!!!! dynamic load balancing
         MINE = MINE + 1
         IF(MINE.GT.NEXT) CALL DDI_DLBNEXT(NEXT)
         IF(NEXT.NE.MINE) CYCLE

!!!!! rank with local data
         IRANK=NPROC-1
         DO IPROC=0,NPROC-2
            IF(IA.GE.ddiWSA(IPROC,2)-N_DLB_BUFF .AND. IA.LT.ddiWSA(IPROC+1,1)) THEN
               IRANK=IPROC
               EXIT
            ENDIF
         ENDDO

!!!!! XLOV_OMP (shared module)
         IF(ME.EQ.IRANK) THEN
            CALL RIMP2OMP_XLOV__AOOPOO                                  &
                (Nthreads, IA, POO, BDOO,BDOV(1,1,IA), NAUXBASD,NOCC,NVIR)
         ELSE
            CALL DDI_GET(D_B, 1,NDO, IA+NOCC,IA+NOCC, BDO)
            CALL RIMP2OMP_XLOV__AOOPOO                                  &
                (Nthreads,IA, POO, BDOO,BDO, NAUXBASD,NOCC,NVIR)
         ENDIF
      ENDDO
      CALL DDI_DLBRESET
      DEALLOCATE(BDOO,BDOV,BDO)

#if __bug==1
      endw = omp_get_wtime()
      WRITE(*,'(A40,I5,F10.1)') "XLOV_DYNAMIC: AOO_POO",ME,endw-startw_dyn
      WRITE(*,'(A40,I5,F10.1)') "XLOV_TOTAL: AOO_POO",ME,endw-startw
#endif


!///////////////////////////////////////////////////////////////////////
!////                                                            ///////
!////                  STEP 2: XLOV <- AVV*PVV                   ///////
!////                                                            ///////
!///////////////////////////////////////////////////////////////////////


#if __bug==1
      startw = omp_get_wtime()
#endif

!!!!! XLOV <- AVV*PVV : static chunk

      ALLOCATE(BDVV(NAUXBASD,NVIR,LddiVirtStart:LddiVirtEnd))
      ALLOCATE(BDVO(NAUXBASD,NVIR,LddiOccStart:LddiOccEnd))
      CALL DDI_GET(D_B, 1+NDO,NDV+NDO, LddiVirtStart+NOCC,LddiVirtEnd+NOCC, BDVV)
      CALL DDI_GET(D_B, 1+NDO,NDV+NDO, LddiOccStart,LddiOccEnd, BDVO)


!$OMP PARALLEL NUM_THREADS(Nthreads)                                    &
!$omp default(none)                                                     &
!$omp shared(N_DLB_BUFF)                                                &
!$omp shared(LddiVirtStart,LddiVirtEnd, LddiOccStart,LddiOccEnd)        &
!$omp shared(PVV,BDVV,BDVO, D_B,NDO,NDV, NOCC,NVIR,NAUXBASD)            &
!$omp private(II,IB)

      DEALLOCATE(AOO)
      ALLOCATE(VV(NVIR,NVIR))
      VV=0.0D00
      ALLOCATE(BDV(NAUXBASD,NVIR))

!$omp do schedule(DYNAMIC)
      DO II=1,NOCC
         IF(II .lt. LddiOccStart .or. II .gt. LddiOccEnd) THEN
!$omp critical
            CALL DDI_GET(D_B, 1+NDO,NDV+NDO, II,II, BDV)
!$omp end critical
            DO IB=LddiVirtStart,LddiVirtEnd - N_DLB_BUFF
               CALL RIMP2_XLOV__AVVPVV                                  &
                   (XLOV_OMP,IB,II,PVV,VV, BDVV(1,1,IB),BDV,            &
                    NAUXBASD,NOCC,NVIR)
            ENDDO
         ELSE
            DO IB=LddiVirtStart,LddiVirtEnd - N_DLB_BUFF
               CALL RIMP2_XLOV__AVVPVV                                  &
                   (XLOV_OMP,IB,II, PVV,VV, BDVV(1,1,IB),BDVO(1,1,II),  &
                    NAUXBASD,NOCC,NVIR)
            ENDDO
         ENDIF
      ENDDO
!$omp end do

!$OMP END PARALLEL

#if __bug==1
      endw_static = omp_get_wtime()
      WRITE(*,'(A40,I5,F10.1)') "XLOV_STATIC: AVV_PVV",ME,endw_static-startw
#endif

!!!!! XLOV <- AVV*PVV : DBL

      startw_dyn = omp_get_wtime()

      ALLOCATE(BI(NAUXBASD,NVIR))
      ALLOCATE(BB(NAUXBASD,NVIR))

      MINE=-1
      NEXT=-1
      DO IB=1,NVIR

         IF(I_DLB_INDEX(IB).EQ.0) CYCLE

         ! dynamic load balancing
         MINE = MINE + 1
         IF(MINE.GT.NEXT) CALL DDI_DLBNEXT(NEXT)
         IF(NEXT.NE.MINE) CYCLE

         ! rank with local data
         IRANK=NPROC-1
         DO IPROC=0,NPROC-2
            IF(IB.GE.ddiWSA(IPROC,2)-N_DLB_BUFF .AND. IB.LT.ddiWSA(IPROC+1,1)) THEN
               IRANK=IPROC
               EXIT
            ENDIF
         ENDDO

         IF(ME.NE.IRANK) CALL DDI_GET(D_B,1+NDO,NDV+NDO,IB+NOCC,IB+NOCC,BB)

         DO II=1,NOCC

            IF(II .GE. LddiOccStart .and. II .LE. LddiOccEnd) THEN
               CALL DCOPY(NAUXBASD*NVIR, BDVO(1,1,II),1, BI,1)
            ELSE
               CALL DDI_GET(D_B, 1+NDO,NDV+NDO, II,II, BI)
            ENDIF

            IF(ME.EQ.IRANK) THEN
               CALL RIMP2OMP_XLOV__AVVPVV                               &
                   (Nthreads, IB,II, PVV,VV, BDVV(1,1,IB),BI,NAUXBASD,NOCC,NVIR)
            ELSE
               CALL RIMP2OMP_XLOV__AVVPVV                               &
                   (Nthreads, IB,II, PVV,VV, BB,BI,NAUXBASD,NOCC,NVIR)
            ENDIF

         ENDDO !II

      ENDDO !IB
      CALL DDI_DLBRESET


!$OMP PARALLEL NUM_THREADS(Nthreads)
      DEALLOCATE(VV,BDV)
!$omp critical (XLOV_critical)
      XLOV=XLOV+XLOV_OMP
!$omp end critical (XLOV_critical)
      DEALLOCATE(XLOV_OMP)
!$OMP END PARALLEL

#if __bug==1
      endw = omp_get_wtime()
      WRITE(*,'(A40,I5,F10.1)') "XLOV_DYNAMIC: AVV_PVV",ME,endw-startw_dyn
      WRITE(*,'(A40,I5,F10.1)') "XLOV_TOTAL: AVV_PVV",ME,endw-startw
#endif

      END !***********************************************************





      SUBROUTINE RIMP2_XLOV__AOOPOO                                     &
                (XLOV, IA, POO,AOO, BDOO,BDO, NAUXBASD,NOCC,NVIR)

      implicit double precision(a-h,o-z)

      double precision :: XLOV(NOCC,NVIR)
      double precision :: POO(NOCC,NOCC)
      double precision :: AOO(NOCC,NOCC)

      double precision :: BDOO(NAUXBASD,NOCC,NOCC)
      double precision :: BDO(NAUXBASD,NOCC)


      DO II=1,NOCC

         CALL DGEMM('T','N', NOCC*NOCC,1,NAUXBASD,                      &
                     4.0D00, BDOO,NAUXBASD,                             &
                             BDO(1,II),NAUXBASD,                        &
                     0.0D00, AOO,NOCC*NOCC)

         CALL DGEMM('T','N', NOCC,NOCC,NAUXBASD,                        &
                    -2.0D00, BDO,NAUXBASD,                              &
                             BDOO(1,1,II),NAUXBASD,                     &
                     1.0D00, AOO, NOCC)

         TMP = DDOT(NOCC*NOCC, AOO,1, POO,1)
         XLOV(II,IA)=XLOV(II,IA) - TMP

      ENDDO

      END




      SUBROUTINE RIMP2OMP_XLOV__AOOPOO                                  &
                (Nthreads, IA, POO, BDOO,BDO, NAUXBASD,NOCC,NVIR)

      use omp_lib
      use Rimp2_Shared_Data,only : XLOV_OMP,AOO

      implicit double precision(a-h,o-z)


      ! double precision :: XLOV(NOCC,NVIR)
      double precision :: POO(NOCC,NOCC)

      double precision :: BDOO(NAUXBASD,NOCC,NOCC)
      double precision :: BDO(NAUXBASD,NOCC)

      integer,allocatable :: ompISE(:,:)

      logical :: lomp


      ! Nthreads = omp_get_max_threads()
      ! IF(Nthreads.gt.NOCC) Nthreads=NOCC

      ALLOCATE(ompISE(0:Nthreads-1,2))
      CALL WorkSharedArray(ompISE,Nthreads,1,NOCC,LOMP)


!$OMP PARALLEL NUM_THREADS(Nthreads)                                    &
!$omp default(none)                                                     &
!$omp shared(POO,BDOO,BDO, ompISE,IA,NVIR,NOCC,NAUXBASD)                &
!$omp private(II,TMP, IthreadID,Istart,Iend,Ichunk)

      IthreadID = omp_get_thread_num()
      Istart = ompISE(IthreadID,1)
      Iend = ompISE(IthreadID,2)

      DO II=Istart,Iend
         CALL DGEMM('T','N', NOCC*NOCC,1,NAUXBASD,                      &
                     4.0D00, BDOO,NAUXBASD,                             &
                             BDO(1,II),NAUXBASD,                        &
                     0.0D00, AOO,NOCC*NOCC)
         CALL DGEMM('T','N', NOCC,NOCC,NAUXBASD,                        &
                    -2.0D00, BDO,NAUXBASD,                              &
                             BDOO(1,1,II),NAUXBASD,                     &
                     1.0D00, AOO, NOCC)
         TMP = DDOT(NOCC*NOCC, AOO,1, POO,1)
         XLOV_OMP(II,IA)=XLOV_OMP(II,IA) - TMP
      ENDDO

!$OMP END PARALLEL

      END




      SUBROUTINE RIMP2_XLOV__AVVPVV                                     &
                (XLOV,IB,II, PVV,VV, BB,BI, NAUXBASD,NOCC,NVIR)

      ! use omp_lib
      ! use Rimp2_Shared_Data,only : XLOV,VV

      implicit double precision(a-h,o-z)

      double precision :: XLOV(NOCC,NVIR)
      double precision :: PVV(NVIR,NVIR)
      double precision :: VV(NVIR,NVIR)

      double precision :: BB(NAUXBASD,NVIR)
      double precision :: BI(NAUXBASD,NVIR)


      CALL DGEMM('T','N', NVIR,NVIR,NAUXBASD,                           &
                  1.0D00, BI,NAUXBASD,                                  &
                          BB,NAUXBASD,                                  &
                  0.0D00, VV,NVIR)

      TMP=DDOT(NVIR*NVIR, VV,1, PVV,1)
      XLOV(II,IB)=XLOV(II,IB)-TMP-TMP

      DO IC=1,NVIR
         pcb=PVV(IC,IB)
         DO IA=1,NVIR
            TMP=VV(IA,IC)*pcb !PVV(IC,IB)
            XLOV(II,IA)=XLOV(II,IA)+TMP+TMP+TMP+TMP
         ENDDO
      ENDDO


      END !*************************************************************





      SUBROUTINE RIMP2OMP_XLOV__AVVPVV                                  &
                (Nthreads, IB,II, PVV,VV, BB,BI, NAUXBASD,NOCC,NVIR)

      use omp_lib
      use Rimp2_Shared_Data,only : XLOV_OMP !,VV

      implicit double precision(a-h,o-z)

      ! double precision :: XLOV(NOCC,NVIR)
      double precision :: PVV(NVIR,NVIR)
      double precision :: VV(NVIR,NVIR)

      double precision :: BB(NAUXBASD,NVIR)
      double precision :: BI(NAUXBASD,NVIR)

      double precision,save :: TMP
!$omp threadprivate(TMP)

      ! Nthreads = omp_get_max_threads()
!$OMP PARALLEL NUM_THREADS(Nthreads)                                    &
!$omp shared(NVIR,NAUXBASD)                                             &
!$omp shared(BI,BB,VV,PVV)                                              &
!$omp private(pcb)

      TMP=0.0D00

!$omp do schedule(guided) ! firstprivate(TMP)
      DO IA=1,NVIR
         CALL DGEMM('T','N', NVIR,1,NAUXBASD,                           &
                     1.0D00, BI,NAUXBASD,                               &
                             BB(1,IA),NAUXBASD,                         &
                     0.0D00, VV(1,IA),NVIR)

         TMP=TMP+DDOT(NVIR, VV(1,IA),1, PVV(1,IA),1)
         ! XLOV_OMP(II,IB)=XLOV_OMP(II,IB)-TMP-TMP
      ENDDO
!$omp end do

      XLOV_OMP(II,IB)=XLOV_OMP(II,IB)-TMP-TMP

!$omp do schedule(guided)
      DO IC=1,NVIR
         pcb=PVV(IC,IB)
         DO IA=1,NVIR
            TMP=VV(IA,IC)*pcb !PVV(IC,IB)
            XLOV_OMP(II,IA)=XLOV_OMP(II,IA)+TMP+TMP+TMP+TMP
         ENDDO
      ENDDO
!$omp end do

!$OMP END PARALLEL

      END !*************************************************************


!>*module rimp2grd   *deck RIMP2_WMP2
!>
!>     @brief   transform energy-weighted DM to AO basis
!>     @author  buu
!>     @date    sept 8, 2018
!>
!>     @detail  transform energy-weighted DM to AO basis
!>

      SUBROUTINE RIMP2_WMP2                                             &
                (WMP2, WVV,WOV,WOO,VEC,NAUXBASD,NCOR,NACT,NVIR,NORB,NBF)
      use omp_lib

      implicit double precision(a-h,o-z)

      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      logical :: GOPARR,DSKWRK,MASWRK

!!!!! OUTPUT
      double precision :: WMP2(NBF,NBF)

!!!!! INPUT
      double precision :: WVV(NVIR,NVIR)
      double precision :: WOV(NCOR+NACT,NVIR)
      double precision :: WOO(NCOR+NACT,NCOR+NACT)
      double precision :: VEC(NBF,NORB)

      double precision,allocatable :: WB(:,:),WNO(:,:)

      integer,ALLOCATABLE :: ompISE(:,:)

      logical :: lomp

!!!!! PARAMETERS
      NOCC=NCOR+NACT


!!!!! WORKSHARE FOR NVIR LOOP
      ALLOCATE(ompISE(0:NPROC-1,2))
      CALL rimp2_WSA(ompISE,NPROC,1,NVIR,LOMP)
      LddiVirtStart=ompISE(ME,1)
      LddiVirtEnd=ompISE(ME,2)


!!!!! WORKSHARE FOR NBF LOOP
      CALL rimp2_WSA(ompISE,NPROC,1,NBF,LOMP)
      LddiNbfStart=ompISE(ME,1)
      LddiNbfEnd=ompISE(ME,2)


!!!!! BUFFER
      ALLOCATE(WNO(NBF,NOCC))
      WNO=0.0D00

      ALLOCATE(WB(NBF,NVIR))
      WB=0.0D00

      Nthreads=omp_get_max_threads()

!$OMP PARALLEL NUM_THREADS(Nthreads)                                    &
!$omp default(none)                                                     &
!$omp shared(NOCC,NBF,VEC,WOO,WB,WNO)                                   &
!$omp shared(LddiVirtStart,LddiVirtEnd)                                 &
!$omp private(II)                                                       &
!$omp shared(NVIR,WOV)                                                  &
!$omp private(IA)                                                       &
!$omp shared(WVV)


!!!!! WMP2 <- WOO
!$OMP DO SCHEDULE(guided)
      DO II=1,NOCC
         CALL DGEMM('N','N', NBF,1,NOCC,                                &
                     1.0D00, VEC,NBF,                                   &
                              WOO(1,II),NOCC,                           &
                     0.0D00, WNO(1,II),NBF)
      ENDDO
!$OMP END DO


!!!!! WMP2 <- WOV & WVV
!$OMP DO SCHEDULE(guided)
      ! DO IA=1,NVIR
      DO IA=LddiVirtStart,LddiVirtEnd
         CALL DGEMM('N','N', NBF,1,NOCC,                                &
                     1.0D00, VEC,NBF,                                   &
                              WOV(1,IA),NOCC,                           &
                     0.0D00, WB(1,IA),NBF)

         CALL DGEMM('N','N', NBF,1,NVIR,                                &
                     1.0D00, VEC(1,1+NOCC),NBF,                         &
                              WVV(1,IA),NVIR,                           &
                     1.0D00, WB(1,IA),NBF)
      ENDDO
!$OMP END DO
!$OMP END PARALLEL

!!!!! AllReduce:: WB(NBF,NVIR)
      CALL DDI_GSUMF(6011,WB,NBF*NVIR)


      WMP2=0.0D00

!$OMP PARALLEL DO SCHEDULE(guided) NUM_THREADS(Nthreads)                &
!$omp default(none)                                                     &
!$omp shared(WMP2,NOCC,NBF,VEC,WB,WNO)                                  &
!$omp shared(LddiNbfStart,LddiNbfEnd)                                   &
!$omp shared(NVIR)                                                      &
!$omp private(IP)


      ! DO IP=1,NBF
      DO IP=LddiNbfStart,LddiNbfEnd


!!!!! WMP2 <- WOO
         CALL DGEMM('N','T', NBF,1,NOCC,                                &
                     1.0D00, WNO,NBF,                                   &
                              VEC(IP,1),NBF,                            &
                     0.0D00, WMP2(1,IP),NBF)


!!!!! WMP2 <- WOV & WVV
         CALL DGEMM('N','T', NBF,1,NVIR,                                &
                     1.0D00, WB,NBF,                                    &
                              VEC(IP,1+NOCC),NBF,                       &
                     1.0D00, WMP2(1,IP),NBF)
      ENDDO

!$OMP END PARALLEL DO


!!!!! AllReduce:: WMP2(NBF,NBF)
      CALL DDI_GSUMF(6012,WMP2,NBF*NBF)

      END !***********************************************************



!>*module rimp2grd   *deck RIMP2_WAO_PAO
!>
!>     @brief   prepare DMs for 1-electron gradient
!>     @author  buu
!>     @date    sept 8, 2018
!>
!>     @detail  prepare DMs for 1-electron gradient
!>
      SUBROUTINE RIMP2_WAO_PAO                                          &
                (PMP2,PSCF, WMP2,EIG,VEC, NCOR,NACT,NVIR,NORB,NBF)

      implicit double precision(a-h,o-z)

      COMMON /IOFILE/ IR,IW,IP,IJK,IPK,IDAF,NAV,IODA(950)


      double precision PMP2(NBF,NBF)
      double precision PSCF(NBF,NBF)
      double precision WMP2(NBF,NBF)
      double precision EIG(NCOR+NACT)
      double precision VEC(NBF,NCOR+NACT)

      double precision,allocatable :: WAO(:),PAO(:)


!!!!! PARAMETERS
      NOCC=NCOR+NACT
      NBF3=(NBF*NBF+NBF)/2


!!!!! ALLOCATE ARRAYS
      ALLOCATE(WAO(NBF3))
      ALLOCATE(PAO(NBF3))

      IJ= 0
      DO I= 1,NBF
         DO J= 1,I
            IJ= IJ+1
            PAO(IJ)= PMP2(I,J)+PMP2(J,I)
            WAO(IJ)= WMP2(I,J)+WMP2(J,I)
         ENDDO
      ENDDO

      PAO = 0.5D00*PAO
      WAO = 0.5D00*WAO


!!!!! HF + MP2 ENERGY-WEIGHTED DENSITY
      IJ= 0
      DO I= 1,NBF
         DO J= 1,I
            IJ= IJ+1
            DUM= 0.0D00
            DO K= 1,NOCC
               DUM= DUM-EIG(K)*VEC(I,K)*VEC(J,K)
            ENDDO
            WAO(IJ)= WAO(IJ)+2.0D00*DUM
         ENDDO
      ENDDO

!!!!! SAVE WAO ON 309, PAO ON 307
      CALL DAWRIT(IDAF,IODA,WAO,NBF3,309,0)
      CALL DAWRIT(IDAF,IODA,PAO,NBF3,307,0)


!!!!! FORM PSCF
      CALL DGEMM('N','T',NBF,NBF,NOCC,                                  &
                  2.0D00, VEC,NBF,                                      &
                          VEC,NBF,                                      &
                  0.0D00, PSCF,NBF)


!!!!! TRANSFORM PSCF INTO TRIANGULAR FORM (WAO)
      IJ= 0
      DO I= 1,NBF
         DO J= 1,I
            IJ= IJ+1
            WAO(IJ)= PSCF(I,J)+PSCF(J,I)
         ENDDO
      ENDDO
      WAO = 0.5D00*WAO


!!!!! SAVE PSCF (IN TRIANGULAR FORM) INTO 308
      CALL DAWRIT(IDAF,IODA,WAO,NBF3,308,0)


!!!!! HF+MP2 DENSITY AND SAVE IT INTO 16
      PAO = PAO+WAO
      CALL DAWRIT(IDAF,IODA,PAO,NBF3,16,0)

      DEALLOCATE(WAO,PAO)

      END !***********************************************************





!>*module rimp2grd   *deck RIMP2_WOO_WOV
!>
!>     @brief   complete energy-weighted DM
!>     @author  buu
!>     @date    sept 8, 2018
!>
!>     @detail  complete energy-weighted DM
!>
      SUBROUTINE RIMP2_WOO_WOV                                          &
                (WOO,WOV, POO,ZOV,PVV,VEC,EIG,                          &
                 NAUXBASD,NCOR,NACT,NVIR,NORB,NBF)

      use omp_lib
      use Rimp2_Shared_Data,only : D_B

      implicit double precision(a-h,o-z)

      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      logical :: GOPARR,DSKWRK,MASWRK

      double precision :: WOO(NCOR+NACT,NCOR+NACT)
      double precision :: WOV(NCOR+NACT,NVIR)

      double precision :: POO(NCOR+NACT,NCOR+NACT)
      double precision :: ZOV(NCOR+NACT,NVIR)
      double precision :: PVV(NVIR,NVIR)

      double precision :: VEC(NBF,NORB)
      double precision :: EIG(NBF)

      double precision,allocatable :: BDVV(:,:,:),BDOO(:,:,:),BDOV(:,:,:),WOO_DDI(:,:)

      double precision,save,allocatable :: WOO_OMP(:,:),AOO(:,:),QOO(:,:),BDO(:,:),AVO(:,:)
!$omp threadprivate(WOO_OMP,AOO,QOO,BDO,AVO)

      integer,ALLOCATABLE :: ompISE(:,:)

      logical :: lomp

!!!!! PARAMETERS
      NOCC=NCOR+NACT

      NXC=NAUXBASD*NCOR
      NXA=NAUXBASD*NACT
      NXO=NAUXBASD*NOCC
      NXV=NAUXBASD*NVIR
      NXN=NAUXBASD*NORB


!!!!! WOV(i,a) <- ZOV(i,a) * e(i) * 2.0D00
      FORALL(II=1:NOCC,IA=1:NVIR)
         WOV(II,IA)=WOV(II,IA) - ZOV(II,IA)*(EIG(II)+EIG(II))
      ENDFORALL


!!!!! DDI WORK DISTRIBUTING FOR NOCC LOOP
      ALLOCATE(ompISE(0:NPROC-1,2))
      CALL rimp2_WSA(ompISE,NPROC,1,NOCC,LOMP)
      LddiOccStart=ompISE(ME,1)
      LddiOccEnd=ompISE(ME,2)


!!!!! DDI WORK DISTRIBUTING FOR NVIR LOOP
      CALL rimp2_WSA(ompISE,NPROC,1,NVIR,LOMP)
      LddiVirtStart=ompISE(ME,1)
      LddiVirtEnd=ompISE(ME,2)

      Nthreads=omp_get_max_threads()

      ALLOCATE(WOO_DDI(NOCC,NOCC))
      WOO_DDI=0.0D00

      ALLOCATE(BDOO(NAUXBASD,NOCC,NOCC))
      CALL DDI_GET(D_B, 1,NXO, 1,NOCC, BDOO)


!!!!! TURN OFF DYNAMIC THREADS
!!!!! TO KEEP threadprivate DATA CONSISTENT
      CALL OMP_SET_DYNAMIC(.FALSE.)

      startw = omp_get_wtime()
!$OMP PARALLEL NUM_THREADS(Nthreads)                                    &
!$omp default(none)                                                     &
!$omp shared(LddiOccStart,LddiOccEnd,NOCC,NAUXBASD)                     &
!$omp private(pkl,pklT,LL,KK)                                           &
!$omp shared(BDOO,POO,WOO_DDI)
      ALLOCATE(WOO_OMP(NOCC,NOCC))
      ALLOCATE(AOO(NOCC,NOCC))
      ALLOCATE(QOO(NOCC,NOCC))
      WOO_OMP=0.0D00
!$omp do schedule(guided)
      ! DO LL=1,NOCC
      DO LL=LddiOccStart,LddiOccEnd

!!!!! OFF DIAGONAL
         DO KK=1,LL-1
            CALL DGEMM('T','N', NOCC,NOCC,NAUXBASD,                     &
                        1.0D00, BDOO(1,1,LL),NAUXBASD,                  &
                                BDOO(1,1,KK),NAUXBASD,                  &
                        0.0D00, AOO,NOCC)
            QOO=transpose(AOO)
            CALL DGEMV('T', NOCC,NOCC,                                  &
                       -4.0D00, AOO,NOCC,                               &
                                POO(1,LL),1,                            &
                        1.0D00, WOO_OMP(1,KK),1)
            CALL DGEMV('T', NOCC,NOCC,                                  &
                       -4.0D00, QOO,NOCC,                               &
                                 POO(1,KK),1,                           &
                        1.0D00, WOO_OMP(1,LL),1)
            pkl=2.0D00*POO(KK,LL)
            pklT=2.0D00*POO(LL,KK)
            WOO_OMP=WOO_OMP+AOO*pkl+QOO*pklT
         ENDDO !KK

!!!!! DIAGONAL
         CALL DGEMM('T','N', NOCC,NOCC,NAUXBASD,                        &
                     1.0D00, BDOO(1,1,LL),NAUXBASD,                     &
                             BDOO(1,1,LL),NAUXBASD,                     &
                     0.0D00, AOO,NOCC)
         CALL DGEMV('T', NOCC,NOCC,                                     &
                    -4.0D00, AOO,NOCC,                                  &
                             POO(1,LL),1,                               &
                     1.0D00, WOO_OMP(1,KK),1)
         pkl=2.0D00*POO(KK,LL)
         WOO_OMP=WOO_OMP+AOO*pkl
      ENDDO !LL
!$omp end do
      DEALLOCATE(AOO)
!$OMP END PARALLEL

#if __bug==1
      endw = omp_get_wtime()
      WRITE(*,'(A40,I5,F10.1)') "RIMP2_WOO_WOV: AOO*POO ",ME,endw-startw
#endif

!!!!! PVV * AVV
      startw = omp_get_wtime()
      ALLOCATE(BDVV(NAUXBASD,NVIR,LddiVirtStart:LddiVirtEnd))
      CALL DDI_GET(D_B, 1+NXO,NXV+NXO, LddiVirtStart+NOCC,LddiVirtEnd+NOCC, BDVV)
!$OMP PARALLEL NUM_THREADS(Nthreads)                                    &
!$omp default(none)                                                     &
!$omp shared(LddiVirtStart,LddiVirtEnd,NXV,NXO,NOCC,NVIR,NAUXBASD)      &
!$omp shared(D_B,PVV,BDOO,BDVV)                                         &
!$omp private(IB,II)
      ALLOCATE(AVO(NVIR,NOCC))
!$omp do schedule(guided)
      ! DO IB=1,NVIR
      DO IB=LddiVirtStart,LddiVirtEnd
         DO II=1,NOCC
            CALL DGEMM('T','N', NVIR,NOCC,NAUXBASD,                     &
                        4.0D00, BDVV(1,1,IB),NAUXBASD,                  &
                                BDOO(1,1,II),NAUXBASD,                  &
                        0.0D00, AVO,NVIR)
            CALL DGEMV('T', NVIR,NOCC,                                  &
                        1.0D00, AVO,NVIR,                               &
                                PVV(1,IB),1,                            &
                        1.0D00, WOO_OMP(1,II),1)
         ENDDO !II
      ENDDO
!$omp end do
      DEALLOCATE(AVO)
!$OMP END PARALLEL
      DEALLOCATE(BDVV)

#if __bug==1
      endw = omp_get_wtime()
      WRITE(*,'(A40,I5,F10.1)') "RIMP2_WOO_WOV: AVO*PVV ",ME,endw-startw
#endif

!!!!! PVV*AVVOO
      startw = omp_get_wtime()

      IF(.FALSE.) THEN

      ALLOCATE(BDOV(NAUXBASD,NOCC,LddiVirtStart:LddiVirtEnd))
      CALL DDI_GET(D_B,1,NXO,LddiVirtStart+NOCC,LddiVirtEnd+NOCC,BDOV)
!$OMP PARALLEL NUM_THREADS(Nthreads)                                    &
!$omp default(none)                                                     &
!$omp shared(LddiVirtStart,LddiVirtEnd,NXO,NOCC,NVIR,NAUXBASD)          &
!$omp shared(D_B,PVV,BDOV)                                              &
!$omp private(pab,pbb)                                                  &
!$omp private(IA,IB,II)
      ALLOCATE(BDO(NAUXBASD,NOCC))
!$omp do schedule(guided)
      DO IB=LddiVirtStart,LddiVirtEnd
         DO IA=1,IB-1
            IF(IA.GE.LddiVirtStart .AND. IA.LE.LddiVirtEnd) THEN
               CALL DCOPY(NXO, BDOV(1,1,IA),1, BDO,1)
            ELSE
               CALL DDI_GET(D_B, 1,NXO, IA+NOCC,IA+NOCC, BDO)
            ENDIF
            pab=-2.0D00*PVV(IA,IB)
            CALL DGEMM('T','N', NOCC,NOCC,NAUXBASD,                     &
                        pab, BDO,NAUXBASD,                              &
                             BDOV(1,1,IB),NAUXBASD,                     &
                        0.0D00, QOO,NOCC)
            WOO_OMP=WOO_OMP+QOO+transpose(QOO)
         ENDDO !IA
         pbb=-PVV(IB,IB)-PVV(IB,IB)
         CALL DGEMM('T','N', NOCC,NOCC,NAUXBASD,                        &
                     pbb, BDOV(1,1,IB),NAUXBASD,                        &
                          BDOV(1,1,IB),NAUXBASD,                        &
                     0.0D00, QOO,NOCC)
         WOO_OMP=WOO_OMP+QOO
      ENDDO
!$omp end do
      DEALLOCATE(BDO)
!$OMP END PARALLEL
      DEALLOCATE(BDOV)


      ELSEIF(.TRUE.) THEN


      ALLOCATE(BDOV(NAUXBASD,NOCC,NVIR))
      CALL DDI_GET(D_B, 1,NXO, 1+NOCC,NVIR+NOCC, BDOV)
      NEXT = -1
      MINE = -1
      DO IB=1,NVIR
      ! DO IB=LddiVirtStart,LddiVirtEnd
         MINE = MINE + 1
         IF(MINE.GT.NEXT) CALL DDI_DLBNEXT(NEXT)
         IF(NEXT.NE.MINE) CYCLE
!$OMP PARALLEL NUM_THREADS(Nthreads)                                    &
!$omp default(none)                                                     &
!$omp shared(LddiVirtStart,LddiVirtEnd,NXO,NOCC,NVIR,NAUXBASD)          &
!$omp shared(D_B,PVV,BDOV)                                              &
!$omp private(pab,pbb)                                                  &
!$omp private(IA,II,JJ) &
!$omp shared(IB)
!$omp do schedule(guided)
         DO IA=1,IB-1
            pab=-PVV(IA,IB)-PVV(IA,IB)
            CALL DGEMM('T','N', NOCC,NOCC,NAUXBASD,                     &
                        pab, BDOV(1,1,IA),NAUXBASD,                     &
                             BDOV(1,1,IB),NAUXBASD,                     &
                        0.0D00, QOO,NOCC)
            WOO_OMP=WOO_OMP+QOO+transpose(QOO)
         ENDDO !IA
!$omp end do nowait
         pbb=-PVV(IB,IB)-PVV(IB,IB)
!$omp do schedule(guided)
         DO II=1,NOCC
            CALL DGEMM('T','N', NOCC,1,NAUXBASD,                        &
                        pbb, BDOV(1,1,IB),NAUXBASD,                     &
                             BDOV(1,II,IB),NAUXBASD,                    &
                        1.0D00, WOO_OMP(1,II),NOCC)
         ENDDO
!$omp end do
!$OMP END PARALLEL
      ENDDO
      CALL DDI_DLBRESET
      DEALLOCATE(BDOV)


      ELSEIF(.FALSE.) THEN


!    DDI WORK DISTRIBUTING FOR NVIR SYMMETRY-UTILIZED LOOP
      CALL RIMP2_SYMDIV (LddiVirtSymStart,LddiVirtSymEnd, ME,NPROC,NVIR)
      ALLOCATE(BDOV(NAUXBASD,NOCC,LddiVirtSymStart:LddiVirtSymEnd))
      CALL DDI_GET(D_B, 1,NXO, LddiVirtSymStart+NOCC,LddiVirtSymEnd+NOCC, BDOV)
!$OMP PARALLEL NUM_THREADS(Nthreads)                                    &
!$omp default(none)                                                     &
!$omp shared(LddiVirtSymStart,LddiVirtSymEnd)                           &
!$omp shared(NXO,NOCC,NVIR,NAUXBASD)                                    &
!$omp shared(D_B,PVV,BDOV)                                              &
!$omp private(pab,pbb)                                                  &
!$omp private(IA,IB,II)
      ALLOCATE(BDO(NAUXBASD,NOCC))
!$omp do schedule(guided)
      DO IB=LddiVirtSymStart,LddiVirtSymEnd
         DO IA=1,IB-1
            IF(IA.GE.LddiVirtSymStart .AND. IA.LE.LddiVirtSymEnd) THEN
               CALL DCOPY(NXO, BDOV(1,1,IA),1, BDO,1)
            ELSE
               CALL DDI_GET(D_B, 1,NXO, IA+NOCC,IA+NOCC, BDO)
            ENDIF
            pab=-2.0D00*PVV(IA,IB)
            CALL DGEMM('T','N', NOCC,NOCC,NAUXBASD,                     &
                        pab, BDO,NAUXBASD,                              &
                             BDOV(1,1,IB),NAUXBASD,                     &
                        0.0D00, QOO,NOCC)
            WOO_OMP=WOO_OMP+QOO+transpose(QOO)
         ENDDO !IA
         pbb=-PVV(IB,IB)-PVV(IB,IB)
         CALL DGEMM('T','N', NOCC,NOCC,NAUXBASD,                        &
                     pbb, BDOV(1,1,IB),NAUXBASD,                        &
                          BDOV(1,1,IB),NAUXBASD,                        &
                     0.0D00, QOO,NOCC)
         WOO_OMP=WOO_OMP+QOO
      ENDDO
!$omp end do
      DEALLOCATE(BDO)
!$OMP END PARALLEL
      DEALLOCATE(BDOV)

      ENDIF

#if __bug==1
      endw = omp_get_wtime()
      WRITE(*,'(A40,I5,F10.1)') "RIMP2_WOO_WOV: QOO*PVV ",ME, endw-startw
#endif



      startw = omp_get_wtime()
      ALLOCATE(BDOV(NAUXBASD,NOCC,LddiVirtStart:LddiVirtEnd))
      CALL DDI_GET(D_B, 1,NXO, LddiVirtStart+NOCC,LddiVirtEnd+NOCC, BDOV)
!    ZOV * AOV
!$OMP PARALLEL NUM_THREADS(Nthreads)                                    &
!$omp default(none)                                                     &
!$omp shared(NVIR,NOCC,NAUXBASD,ZOV)                                    &
!$omp shared(LddiVirtStart,LddiVirtEnd)                                 &
!$omp private(zka,I,J,IA,LL,KK)                                         &
!$omp shared(WOO_DDI,BDOV,BDOO)

!$OMP DO SCHEDULE(guided) !IA
      ! DO IA=1,NVIR
      DO IA=LddiVirtStart,LddiVirtEnd
         DO KK=1,NOCC
            CALL DGEMM('T','N', NOCC,NOCC,NAUXBASD,                     &
                        1.0D00, BDOV(1,1,IA),NAUXBASD,                  &
                                BDOO(1,1,KK),NAUXBASD,                  &
                        0.0D00, QOO,NOCC)
            CALL DGEMV('T', NOCC,NOCC,                                  &
                       -4.0D00, QOO,NOCC,                               &
                                ZOV(1,IA),1,                            &
                        1.0D00, WOO_OMP(1,KK),1)
            QOO=QOO+transpose(QOO)
            zka=ZOV(KK,IA)
            WOO_OMP=WOO_OMP+QOO*zka
         ENDDO !KK
      ENDDO !IA
!$omp end do !IA

!$omp critical
      WOO_DDI=WOO_DDI+WOO_OMP
!$omp end critical

      DEALLOCATE(WOO_OMP,QOO) !,AOO)

!$OMP END PARALLEL

#if __bug==1
      endw = omp_get_wtime()
      WRITE(*,'(A40,I5,F10.1)') "RIMP2_WOO_WOV: QOO*ZVO ",ME, endw-startw
#endif


      startw = omp_get_wtime()

!    ACCUMULATE DATA
      CALL DDI_GSUMF(6009,WOO_DDI,NOCC*NOCC)
      WOO=WOO+WOO_DDI

#if __bug==1
      endw = omp_get_wtime()
      WRITE(*,'(A40,I5,F10.1)') "RIMP2_WOO_WOV: REDUCTION ",ME, endw-startw
#endif

      END !***********************************************************



      logical FUNCTION ULRange(N1,N2,N)
            ULRange=.true.
            IF(N.LT.N1 .OR. N.GT.N2) ULRange=.FALSE.
      END !***********************************************************





!>*module rimp2grd   *deck RIMP2_PMP2
!>
!>     @brief   transform DM to AO basis
!>     @author  buu
!>     @date    sept 8, 2018
!>
!>     @detail  transform DM to AO basis
!>

      SUBROUTINE RIMP2_PMP2                                             &
                (PMP2, POO,ZOV,PVV, VEC,NCOR,NACT,NVIR,NORB,NBF)

      use omp_lib

      implicit double precision(a-h,o-z)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      logical :: GOPARR,DSKWRK,MASWRK

!!!!! OUTPUT
      double precision :: PMP2(NBF,NBF)


!!!!! INPUT
      double precision :: POO(NCOR+NACT,NCOR+NACT)
      double precision :: ZOV(NCOR+NACT,NVIR)
      double precision :: PVV(NVIR,NVIR)
      double precision :: VEC(NBF,NORB)


!!!!! LOCAL DATA
      double precision,allocatable :: PNO(:,:),PNV(:,:),PON(:,:)
      integer,allocatable :: ompISE(:,:),ompISEocc(:,:),ompISEvirt(:,:),&
                             ompISEnbf(:,:)

      logical :: lomp

!!!!! PARAMETERS
      NOCC=NCOR+NACT

!!!!! INITIALIZATION
      PMP2=0.0D00


!!!!! DDI WORKSHARE FOR NBF LOOP
      ALLOCATE(ompISE(0:NPROC-1,2))
      CALL rimp2_WSA(ompISE,NPROC,1,NBF,LOMP)
      LddiNbfStart=ompISE(ME,1)
      LddiNbfEnd=ompISE(ME,2)
      DEALLOCATE(ompISE)


!!!!! NUM THREADS FROM ENV SETTING
      Nthreads_env=omp_get_max_threads()


!!!!! MAX NUM THREADS USED TO GUARANTEE MEANINGFUL BOUNDS OF ARRAYS
      Nthreads=MIN(Nthreads_env,NOCC,NVIR,LddiNbfEnd-LddiNbfStart+1)



!!!!! OMP WORKSHARE FOR PAR REGIONS
      ALLOCATE(ompISEocc(0:Nthreads-1,2))
      ALLOCATE(ompISEvirt(0:Nthreads-1,2))
      ALLOCATE(ompISEnbf(0:Nthreads-1,2))

      CALL rimp2_WSA(ompISEocc,Nthreads,1,NOCC,LOMP)
      CALL rimp2_WSA(ompISEvirt,Nthreads,1,NVIR,LOMP)
      CALL rimp2_WSA(ompISEnbf,Nthreads,LddiNbfStart,LddiNbfEnd,LOMP)


!!!!! ALLOCATE AND INIT BUFF ARRAYS
      ALLOCATE(PNO(NBF,NOCC))
      ALLOCATE(PNV(NBF,NVIR))
      ALLOCATE(PON(NOCC,LddiNbfStart:LddiNbfEnd))

      PNO=0.0D00
      PNV=0.0D00
      PON=0.0D00


!$OMP PARALLEL NUM_THREADS(Nthreads)                                    &
!$omp default(none)                                                     &
!$omp shared(ompISEocc,ompISEvirt,ompISEnbf,Nthreads,LOMP)              &
!$omp private(IthreadID,LompOccStart,LompOccEnd)                        &
!$omp private(LompVirtStart,LompVirtEnd)                                &
!$omp private(LompNbfStart,LompNbfEnd)                                  &
!$omp shared(NOCC,NVIR,NBF)                                             &
!$omp shared(VEC,POO,PNO,PVV,PNV,ZOV)                                   &
!$omp private(II,IA,IP)                                                 &
!$omp shared(LddiNbfStart,LddiNbfEnd,PMP2,PON)

      IthreadID=omp_get_thread_num()

      LompOccStart=ompISEocc(IthreadID,1)
      LompOccEnd=ompISEocc(IthreadID,2)

      LompVirtStart=ompISEvirt(IthreadID,1)
      LompVirtEnd=ompISEvirt(IthreadID,2)

      LompNbfStart=ompISEnbf(IthreadID,1)
      LompNbfEnd=ompISEnbf(IthreadID,2)

!!!!! PMP2 <- POO
      ! DO II=1,NOCC
      DO II=LompOccStart,LompOccEnd
         CALL DGEMM('N','N', NBF,1,NOCC,                                &
                    -1.0D00, VEC,NBF,                                   &
                             POO(1,II),NOCC,                            &
                     0.0D00, PNO(1,II),NBF)

      ENDDO

!!!!! PMP2 <- PVV
      ! DO IA=1,NVIR
      DO IA=LompVirtStart,LompVirtEnd
         CALL DGEMM('N','N', NBF,1,NVIR,                                &
                     1.0D00, VEC(1,1+NOCC),NBF,                         &
                              PVV(1,IA),NVIR,                           &
                     0.0D00, PNV(1,IA),NBF)
      ENDDO

!$omp barrier

!!!!! PMP2 <- POO
      ! DO IP=1,NBF
      ! DO IP=LddiNbfStart,LddiNbfEnd
      DO IP=LompNbfStart,LompNbfEnd
         CALL DGEMM('N','T', NBF,1,NOCC,                                &
            1.0D00, PNO,NBF,                                            &
                     VEC(IP,1),NBF,                                     &
            0.0D00, PMP2(1,IP),NBF)
      ENDDO

!!!!! PMP2 <- PVV
      DO IP=LompNbfStart,LompNbfEnd
         CALL DGEMM('N','T', NBF,1,NVIR,                                &
                     1.0D00, PNV,NBF,                                   &
                              VEC(IP,1+NOCC),NBF,                       &
                     1.0D00, PMP2(1,IP),NBF)
      ENDDO

!!!!! PMP2 <- ZOV
      DO IP=LompNbfStart,LompNbfEnd
         CALL DGEMM('N','T', NOCC,1,NVIR,                               &
                     2.0D00, ZOV,NOCC,                                  &
                              VEC(IP,1+NOCC),NBF,                       &
                     0.0D00, PON(1,IP),NOCC)
      ENDDO

      DO IP=LompNbfStart,LompNbfEnd
         CALL DGEMM('N','N', NBF,1,NOCC,                                &
                     1.0D00, VEC,NBF,                                   &
                              PON(1,IP),NOCC,                           &
                    -2.0D00, PMP2(1,IP),NBF)

      ENDDO

!$OMP END PARALLEL



!!!!! REDUCTION :: PMP2(NBF,NBF)
      CALL DDI_GSUMF(6011,PMP2,NBF*NBF)


      END !***********************************************************




!>*module rimp2grd   *deck RIMP2_ZVECTOR
!>
!>     @brief   solve for Z-vector
!>     @author  buu
!>     @date    sept 8, 2018
!>
!>     @detail  solve for Z-vector
!>
      SUBROUTINE RIMP2_ZVECTOR                                          &
                (ZOV, XLOV,EIG, ITERMX,NAUXBASD,NOCC,NVIR,NORB,NBF)

      use omp_lib
      use Rimp2_Shared_Data,only: D_B,D_AIBJ_TMP

      implicit double PRECISION(A-H,O-Z)

      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /IOFILE/ IR,IW,IP,IJK,IPK,IDAF,NAV,IODA(950)
      COMMON /MACHSW/ KDIAG,ICORFL,IXDR,modio,mem10,lpnt10,mem10m
      COMMON /PCMPAR/ IPCM,NFT26,NFT27,IRPPCM,IEF,IP_F,NFMOPCM,IHET
      logical :: GOPARR,DSKWRK,MASWRK,masout
      integer :: KDIAG,ICORFL,IXDR,modio,mem10,lpnt10,mem10m

!!!!! OUTPUT
      double precision :: ZOV(NOCC,NVIR)


!!!!! INPUT
      double precision :: EIG(NBF)
      double precision :: XLOV(NOCC,NVIR)


!!!!! LOCAL
      double precision,allocatable ::                                   &
             GMRES1(:,:),GMRES2(:,:),GMVEC(:,:,:),GMTRN(:,:),           &
             GMWRK1(:,:),GMWRK2(:,:),GMFAC(:),AIBJ(:,:,:),              &
             AIBJMEM(:,:,:,:), AIBJ2(:,:,:),GMVEC_T(:,:,:),             &
             BDOV(:,:,:),BDOO(:,:,:),BDVV(:,:,:),                       &
             QOO(:,:),QOOT(:,:),BDV(:,:)

      integer,allocatable :: ddiWSA(:,:)
      integer :: D_AIBJ

      logical :: lomp




!!!!! PARAMETERS
      NDO=NAUXBASD*NOCC
      NDV=NAUXBASD*NVIR
      NDN=NAUXBASD*(NOCC+NVIR)
      masout=maswrk.and.iand(modio,16).eq.0

      THRESHOLD= 1.D-10


!!!!! DDI WORK DISTRIBUTING NVIR LOOP
      ALLOCATE(ddiWSA(0:NPROC-1,2))
      CALL rimp2_WSA(ddiWSA,NPROC,1,NVIR,LOMP)
      LddiVirtStart=ddiWSA(ME,1)
      LddiVirtEnd=ddiWSA(ME,2)


!!!!! INITIALIZE Z-VECTOR
      IF(.FALSE.) THEN
         ALLOCATE(GMVEC(NOCC,NVIR,0:ITERMX+1))
         ALLOCATE(AIBJ(NOCC,NOCC,NVIR))
      ELSEIF(.TRUE.) THEN
         ALLOCATE(AIBJMEM(NOCC,NOCC,NVIR,LddiVirtStart:LddiVirtEnd))
      ENDIF


      startw=omp_get_wtime()
      IF(.FALSE.) THEN
         IF(MASWRK) WRITE(IW,*) "(SYMMETRICALLY, EARN COMPUTATION; PAY IO, SYNC) FORM AIBJ"
         CALL DDI_CREATE(NOCC*NOCC*NVIR,NVIR,D_AIBJ)
         CALL RIMP2_SYM_AIBJ                                            &
             (D_AIBJ,GMVEC,XLOV,                                        &
              EIG, AIBJ,ITERMX,                                         &
              NAUXBASD,NOCC,NVIR,NORB,NBF)
      ELSEIF(.FALSE.) THEN
         IF(MASWRK) WRITE(IW,*) "(NON-SYMMETRICALLY, PAY COMPUTATION, EARN IO) FORM AIBJ"
         CALL DDI_CREATE(NOCC*NOCC*NVIR,NVIR,D_AIBJ)
         CALL RIMP2_NOSYM_AIBJ                                          &
             (D_AIBJ,GMVEC,XLOV,                                        &
              EIG, AIBJ,ITERMX,                                         &
              NAUXBASD,NOCC,NVIR,NORB,NBF)
      ELSEIF(.FALSE.) THEN
         IF(MASWRK) WRITE(IW,*) "(HIGH MEMORY REQUESTED) FORM AIBJ"
         CALL RIMP2_MEM_AIBJ                                            &
             (AIBJMEM, EIG,                                             &
              NAUXBASD,NOCC,NVIR,NORB,NBF,                              &
              LddiVirtStart,LddiVirtEnd)
      ELSEIF(.TRUE.) THEN
         IF(MASout) WRITE(IW,*) "(STATIC + DLB) FORM AIBJ"
         CALL RIMP2_DLBMEM_AIBJ                                         &
             (AIBJMEM, ddiWSA, EIG,                                     &
              NAUXBASD,NOCC,NVIR,NORB,NBF,                              &
              LddiVirtStart,LddiVirtEnd)
      ENDIF

#if __bug==1
      endw=omp_get_wtime()
      WRITE(*,'(A40,I5,F10.1)') "TIME AIBJ_TOTAL: ",ME, endw-startw
#endif


      ALLOCATE(GMVEC(NOCC,NVIR,0:ITERMX+1))

      CALL RIMP2_INIT_ZVEC                                              &
          (GMVEC, AIBJMEM,XLOV,EIG,                                     &
           NOCC,NVIR,NBF,ITERMX,                                        &
           LddiVirtStart,LddiVirtEnd)

!!!!! destroy ditributed array created inside RIMP2_DLBMEM_AIBJ
      CALL DDI_DESTROY(D_AIBJ_TMP)




!!!!! SECOND PARALLEL REGION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ALLOCATE(GMFAC(ITERMX))
      ALLOCATE(GMWRK2(ITERMX+1,ITERMX+1))
      ALLOCATE(GMRES1(ITERMX+1,ITERMX))
      ALLOCATE(GMRES2(ITERMX+1,ITERMX))
      ALLOCATE(GMTRN(ITERMX+1,ITERMX+1))
      ALLOCATE(GMWRK1(ITERMX+1,ITERMX+1))


      GMRES1 = 0.0D00
      GMRES2 = 0.0D00
      GMTRN  = 0.0D00
      GMWRK1 = 0.0D00


      FORALL(II=1:ITERMX+1)
         GMTRN(II,II)= 1.0D00
         GMWRK1(II,II)= 1.0D00
      ENDFORALL


      GAMMA= DDOT(NOCC*NVIR,GMVEC(1,1,1),1,GMVEC(1,1,1),1)
      GAMMA= SQRT(GAMMA)
      CALL DSCAL(NOCC*NVIR,1.0D00/GAMMA,GMVEC(1,1,1),1)


      Nthreads=omp_get_max_threads()


!!!!! START ITERATIONS
      DO ICL= 1,ITERMX

         IF(.FALSE.) THEN
            ! DO IB= 1,NVIR
            DO IB=LddiVirtStart,LddiVirtEnd
                  CALL DDI_GET(D_AIBJ, 1,NOCC*NOCC*NVIR, IB,IB, AIBJ)
                  CALL DGEMV('N',NOCC,NOCC*NVIR,                        &
                              1.0D00, AIBJ(1,1,1),NOCC,                 &
                                       GMVEC(1,1,ICL),1,                &
                              0.0D00, GMVEC(1,IB,ICL+1),1)
            ENDDO !IB
         ELSEIF(.TRUE.) THEN
!$OMP PARALLEL DO NUM_THREADS(Nthreads)                                 &
!$omp shared(LddiVirtStart,LddiVirtEnd)                                 &
!$omp shared(ICL,NOCC,NVIR)                                             &
!$omp shared(AIBJMEM,GMVEC)                                             &
!$omp private(IB)
            ! DO IB= 1,NVIR
            DO IB=LddiVirtStart,LddiVirtEnd
                  CALL DGEMV('N',NOCC,NOCC*NVIR,                        &
                              1.0D00, AIBJMEM(1,1,1,IB),NOCC,           &
                                      GMVEC(1,1,ICL),1,                 &
                              0.0D00, GMVEC(1,IB,ICL+1),1)
            ENDDO !IB
!$OMP END PARALLEL DO
         ENDIF

!       REDUCTION :: GMVEC
         CALL DDI_GSUMF(1002,GMVEC(1,1,ICL+1),NOCC*NVIR)


      ! CALL RIMP2_GMVEC(D_AIBJ,GMVEC,XLOV, EIG, AIBJ,ITERMX,       &
      !         NAUXBASD,NOCC,NVIR,NORB,NBF, &
      !         ICL,LddiVirtStart,LddiVirtEnd, &
      !         BDOV,BDOO,BDVV)



         DO II= 1,ICL
            GMRES1(II,ICL)=DDOT(NOCC*NVIR, GMVEC(1,1,II),1, GMVEC(1,1,ICL+1),1)
            CALL DAXPY(NOCC*NVIR,-GMRES1(II,ICL), GMVEC(1,1,II),1, GMVEC(1,1,ICL+1),1)
         ENDDO !II

         FAC1=DDOT(NOCC*NVIR,GMVEC(1,1,ICL+1),1, GMVEC(1,1,ICL+1),1)
         FAC1= SQRT(FAC1)
         GMRES1(ICL+1,ICL)= FAC1

         FAC1= 1.0D00/FAC1
         CALL DSCAL(NOCC*NVIR,FAC1,GMVEC(1,1,ICL+1),1)

         CALL DGEMM('N','N',ICL+1,ICL, ICL+1,                           &
                     1.0D00, GMTRN,ITERMX+1,                            &
                              GMRES1,ITERMX+1,                          &
                     0.0D00, GMRES2,ITERMX+1)

         FAC2 = GMRES2(ICL,ICL)*GMRES2(ICL,ICL) + GMRES2(ICL+1,ICL)*GMRES2(ICL+1,ICL)
         FAC2 = SQRT(FAC2)
         FAC2 = 1.0D00/FAC2


         WORK1 = GMRES2(ICL,ICL)*FAC2
         WORK2 = GMRES2(ICL+1,ICL)*FAC2

         GMWRK1(ICL,ICL) = WORK1
         GMWRK1(ICL+1,ICL) = -WORK2
         GMWRK1(ICL,ICL+1) = WORK2
         GMWRK1(ICL+1,ICL+1) = WORK1


         CALL DGEMM('N','N',ICL+1,ICL+1,ICL+1,                          &
                     1.0D00, GMWRK1,ITERMX+1,                           &
                             GMTRN,ITERMX+1,                            &
                     0.0D00, GMWRK2,ITERMX+1)

         FORALL(II=1:ICL+1, JJ=1:ICL+1)
            GMTRN(JJ,II) =GMWRK2(JJ,II)
         ENDFORALL


         GMWRK1(ICL,ICL)      = 1.0D00
         GMWRK1(ICL+1,ICL)    = 0.0D00
         GMWRK1(ICL,ICL+1)    = 0.0D00
         GMWRK1(ICL+1,ICL+1)  = 1.0D00

         NCYCLE= ICL

         IF(MASout) WRITE(IW,9100) ICL,GMTRN(ICL+1,1)*GAMMA
 9100 FORMAT(1X,'CYCLE',1X,I2,2X,'Z-VECTOR ERROR=',1X,E13.6)

!!!!! CHECK ERROR THRESHOLD
         IF(ABS(GMTRN(ICL+1,1)*GAMMA).LT.THRESHOLD) EXIT

      ENDDO !ICL




!!!!! THE FOURTH PARALLEL REGION

      CALL DGEMM('N','N',NCYCLE+1,NCYCLE,NCYCLE+1,                      &
                  1.0D00, GMTRN,ITERMX+1,                               &
                           GMRES1,ITERMX+1,                             &
                  0.0D00, GMRES2,ITERMX+1)


      DO ICYCLE= NCYCLE,1,-1
         QTRN= GMTRN(ICYCLE,1)*GAMMA
         DO II= ICYCLE+1,NCYCLE
            QTRN= QTRN-GMRES2(ICYCLE,II)*GMFAC(II)
         ENDDO
         GMFAC(ICYCLE)= QTRN/GMRES2(ICYCLE,ICYCLE)
      ENDDO


!!!!! Z-VECTOR
      DO IA=1,NVIR
         DO II=1,NOCC
            ZOV(II,IA)=GMVEC(II,IA,0)
         ENDDO
      ENDDO


      DO ICYCLE=1,NCYCLE
         DO IA=1,NVIR
            DO II=1,NOCC
               ZOV(II,IA)=ZOV(II,IA) + GMFAC(ICYCLE)*GMVEC(II,IA,ICYCLE)
            ENDDO
         ENDDO
      ENDDO

      IF(.FALSE.) CALL DDI_DESTROY(D_AIBJ)

      RETURN

      END !***********************************************************




      SUBROUTINE RIMP2_SYM_AIBJ (D_AIBJ,GMVEC,XLOV,EIG, AIBJ,ITERMX,    &
                 NAUXBASD,NOCC,NVIR,NORB,NBF)

      use omp_lib
      use Rimp2_Shared_Data,only: D_B

      implicit double precision(A-H,O-Z)

      logical :: GOPARR,DSKWRK,MASWRK

      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK

      double precision :: GMVEC(NOCC,NVIR,0:ITERMX+1)
      double precision :: XLOV(NOCC,NVIR)
      double precision :: EIG(NBF)
      double precision :: AIBJ(NOCC,NOCC,NVIR)

      double precision,allocatable ::                                   &
             AIBJ2(:,:,:),BDOV(:,:,:),BDOO(:,:,:),                      &
             BDVV(:,:,:),QOO(:,:),QOOT(:,:),BDV(:,:)

      integer,allocatable :: ompISE(:,:)

      integer :: D_AIBJ

      logical :: lomp



!!!!! PARAMETERS
      NDO=NAUXBASD*NOCC
      NDV=NAUXBASD*NVIR
      NDN=NAUXBASD*(NOCC+NVIR)


!!!!! DDI WORK DISTRIBUTING NVIR LOOP
      ALLOCATE(ompISE(0:NPROC-1,2))
      CALL rimp2_WSA(ompISE,NPROC,1,NVIR,LOMP)
      LddiVirtStart=ompISE(ME,1)
      LddiVirtEnd=ompISE(ME,2)
      DEALLOCATE(ompISE)


!!!!! DDI WORK DISTRIBUTING FOR NVIR SYMMETRY-UTILIZED LOOP
      CALL RIMP2_SYMDIV (LddiVirtSymStart,LddiVirtSymEnd, ME,NPROC,NVIR)


!!!!! COPY MO-ENERGY SCALED LAGRANGIAN TO Z-VECTOR
      ALLOCATE(BDOV(NAUXBASD,NOCC,NVIR))
      ALLOCATE(BDOO(NAUXBASD,NOCC,NOCC))
      ALLOCATE(AIBJ2(NOCC,NOCC,NVIR))

      CALL DDI_GET(D_B, 1,NDO, 1+NOCC,NVIR+NOCC, BDOV)
      CALL DDI_GET(D_B, 1,NDO, 1,NOCC, BDOO)


      Nthreads=omp_get_max_threads()


      Nchunk=LddiVirtEnd-LddiVirtStart+1
      IStartBDVV=LddiVirtSymStart
      IEN_DLB_BUFFDVV=LddiVirtSymStart+Nchunk
      ALLOCATE(BDVV(NAUXBASD,NVIR,IStartBDVV:IEN_DLB_BUFFDVV))
      CALL DDI_GET(D_B, 1+NDO,NDN, IStartBDVV+NOCC,IEN_DLB_BUFFDVV+NOCC, BDVV)

!!!!! AIJB = 4*(ai|bj) - (aj|bi)
      ! DO IB= 1,NVIR
      ! DO IB=LddiVirtStart,LddiVirtEnd
      DO IB=LddiVirtSymStart,LddiVirtSymEnd

!$OMP PARALLEL NUM_THREADS(Nthreads)                                    &
!$omp default(none)                                                     &
!$omp shared(LddiVirtSymStart,LddiVirtSymEnd,IStartBDVV,IEN_DLB_BUFFDVV)&
!$omp shared(IB,NDO,NDN,NOCC,NVIR,NBF,NAUXBASD)                         &
!$omp shared(AIBJ,AIBJ2,XLOV,GMVEC)                                     &
!$omp private(QOO,BDV,IA,I,J,IC,ebi)                                    &
!$omp shared(D_B,BDOV,BDOO,BDVV,EIG,eb)

         ALLOCATE(BDV(NAUXBASD,NVIR))
         ALLOCATE(QOO(NOCC,NOCC))

!$omp do schedule(guided)
!!!!! OFF-DIAGONAL ELEMENTS
         DO IA=1,IB !NVIR
!!!!! (ai|bj)
            CALL DGEMM('T','N', NOCC,NOCC,NAUXBASD,                     &
                        1.0D00, BDOV(1,1,IB),NAUXBASD,                  &
                                BDOV(1,1,IA),NAUXBASD,                  &
                        0.0D00, QOO,NOCC)

!!!!! UPPER PART :: OFF-DIAGONAL ELEMENTS
            IF(IB.LE.IEN_DLB_BUFFDVV) THEN
               CALL DCOPY(NAUXBASD*NVIR, BDVV(1,1,IB),1, BDV,1)
            ELSE
               CALL DDI_GET(D_B, 1+NDO,NDN, IB+NOCC,IB+NOCC, BDV)
            ENDIF

            CALL RIMP2_AIBJ(AIBJ(1,1,IA),QOO,                           &
                 BDOO,BDV(1,IA),EIG,IB, NOCC,NAUXBASD,NBF)

            IF(IA.LT.IB) THEN
               IF(IA.GE.IStartBDVV .AND. IA.LE.IEN_DLB_BUFFDVV) THEN
                  CALL DCOPY(NAUXBASD*NVIR, BDVV(1,1,IA),1, BDV,1)
               ELSE
                  CALL DDI_GET(D_B, 1+NDO,NDN, IA+NOCC,IA+NOCC, BDV)
               ENDIF

!!!!! LOWER PART :: OFF-DIAGONAL ELEMENTS
               CALL RIMP2_AIBJ(AIBJ2(1,1,IA),transpose(QOO),            &
                    BDOO,BDV(1,IB),EIG,IA, NOCC,NAUXBASD,NBF)
            ELSE
!!!!! ADD DELTA KRONECKER TO DIAGONAL ELEMENTS
               DO II=1,NOCC
                  AIBJ(II,II,IB)=AIBJ(II,II,IB)+1.0D00
               ENDDO
            ENDIF

         ENDDO
!$omp end do
      DEALLOCATE(BDV,QOO)
!$OMP END PARALLEL

!!!!! STORE LOWER PART
         IStart=NOCC*NOCC*(IB-1)+1
         IEnd=NOCC*NOCC*IB
         CALL DDI_PUT(D_AIBJ, IStart,IEnd,1,IB-1, AIBJ2)


!!!!! STORE UPPER PART
         IStartRow=1
         IEndRow=NOCC*NOCC*IB
         CALL DDI_PUT(D_AIBJ, IStartRow,IEndRow, IB,IB, AIBJ)

      ENDDO !IB

!!!!! DDI BARRIER
      CALL DDI_SYNC(9807)


      ! DO IB= 1,NVIR
      DO IB=LddiVirtStart,LddiVirtEnd
         CALL DCOPY(NOCC,XLOV(1,IB),1,GMVEC(1,IB,1),1)
         CALL DDI_GET(D_AIBJ, 1,NOCC*NOCC*NVIR, IB,IB, AIBJ)
         DO II=1,NOCC
            CALL DGEMV('N',1,NOCC*NVIR,                                 &
                       -1.0D00, AIBJ(II,1,1),NOCC,                      &
                                 GMVEC(1,1,0),1,                        &
                        1.0D00, GMVEC(II,IB,1),1)
         ENDDO

      ENDDO !IB

      END




      SUBROUTINE RIMP2_NOSYM_AIBJ                                       &
                (D_AIBJ,GMVEC,XLOV, EIG, AIBJ,ITERMX,                   &
                 NAUXBASD,NOCC,NVIR,NORB,NBF)

      use omp_lib
      use Rimp2_Shared_Data,only: D_B

      implicit double precision(A-H,O-Z)

      logical :: GOPARR,DSKWRK,MASWRK

      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK

      double precision :: GMVEC(NOCC,NVIR,0:ITERMX+1)
      double precision :: XLOV(NOCC,NVIR)
      double precision :: EIG(NBF)
      double precision :: AIBJ(NOCC,NOCC,NVIR)

      double precision,allocatable ::                                   &
            BDOV(:,:,:),BDOO(:,:,:),BDVV(:,:,:),QOO(:,:),BDV(:,:)

      integer,allocatable :: ompISE(:,:)

      integer :: D_AIBJ

      logical :: lomp

!!!!! PARAMETERS
      NDO=NAUXBASD*NOCC
      NDV=NAUXBASD*NVIR
      NDN=NAUXBASD*(NOCC+NVIR)


!!!!! DDI WORK DISTRIBUTING NVIR LOOP
      ALLOCATE(ompISE(0:NPROC-1,2))
      CALL rimp2_WSA(ompISE,NPROC,1,NVIR,LOMP)
      LddiVirtStart=ompISE(ME,1)
      LddiVirtEnd=ompISE(ME,2)
      DEALLOCATE(ompISE)


!!!!!  COPY MO-ENERGY SCALED LAGRANGIAN TO Z-VECTOR
      ALLOCATE(BDOV(NAUXBASD,NOCC,NVIR))
      ALLOCATE(BDOO(NAUXBASD,NOCC,NOCC))
      ALLOCATE(BDV(NAUXBASD,NVIR))
      CALL DDI_GET(D_B, 1,NDO, 1+NOCC,NVIR+NOCC, BDOV)
      CALL DDI_GET(D_B, 1,NDO, 1,NOCC, BDOO)


      Nthreads=omp_get_max_threads()

      NEXT = -1
      MINE = -1
      DO IB= 1,NVIR
      ! DO IB=LddiVirtStart,LddiVirtEnd

            MINE = MINE + 1
            IF(MINE.GT.NEXT) CALL DDI_DLBNEXT(NEXT)
            IF(NEXT.NE.MINE) CYCLE

         ! CALL DCOPY(NOCC,XLOV(1,IB),1,GMVEC(1,IB,1),1)
         CALL DDI_GET(D_B, 1+NDO,NDN, IB+NOCC,IB+NOCC, BDV)

         AIBJ=0.0D00

         eb=EIG(IB+NOCC)

!$OMP PARALLEL NUM_THREADS(Nthreads)                                    &
!$omp default(none)                                                     &
!$omp shared(IB,NOCC,NVIR,NAUXBASD)                                     &
!$omp shared(AIBJ,XLOV,GMVEC)                                           &
!$omp private(QOO,IA,I,J,IC,ebi)                                        &
!$omp shared(BDOV,BDOO,BDV,EIG,eb)


         ALLOCATE(QOO(NOCC,NOCC))
!$omp do schedule(guided)
         DO IA=1,NVIR

!!!!! AIJB = 4*(ai|bj) - (aj|bi)
            CALL DGEMM('T','N', NOCC,NOCC,NAUXBASD,                     &
                        1.0D00, BDOV(1,1,IB),NAUXBASD,                  &
                                BDOV(1,1,IA),NAUXBASD,                  &
                        0.0D00, QOO,NOCC)

            FORALL(J=1:NOCC,I=1:NOCC)
               AIBJ(I,J,IA)=4.0D00*QOO(I,J)-QOO(J,I)
            ENDFORALL


!!!!! AIJB <- -(ab|ij)
            CALL DGEMM('T','N', NOCC*NOCC,1, NAUXBASD,                  &
                       -1.0D00, BDOO,NAUXBASD,                          &
                                BDV(1,IA),NAUXBASD,                     &
                        1.0D00, AIBJ(1,1,IA),NOCC*NOCC)


!!!!! AIJB = AIBJ/[e(b)-e(i)]
            FORALL(J=1:NOCC,I=1:NOCC)
               AIBJ(I,J,IA)=AIBJ(I,J,IA)/(eb-EIG(I))
            ENDFORALL

         ENDDO !IA
!$omp end do

!!!!! AIJB(I,I) = AIJB(I,I) + 1.0D00
!$omp do schedule(guided)
         DO II=1,NOCC
            AIBJ(II,II,IB)=AIBJ(II,II,IB)+1.0D00
            CALL DGEMV('N',1,NOCC*NVIR,                                 &
                       -1.0D00, AIBJ(II,1,1),NOCC,                      &
                                GMVEC(1,1,0),1,                         &
                        1.0D00, GMVEC(II,IB,1),1)
         ENDDO
!$omp end do

         DEALLOCATE(QOO)

!$OMP END PARALLEL

!!!!! SAVE AIJB TO DDI ARRAY
         CALL DDI_PUT(D_AIBJ, 1, NOCC*NOCC*NVIR, IB,IB, AIBJ)

      ENDDO !IB


      CALL DDI_DLBRESET

      END





      SUBROUTINE RIMP2_MEM_AIBJ(AIBJMEM, EIG,                           &
                 NAUXBASD,NOCC,NVIR,NORB,NBF,                           &
                 LddiVirtStart,LddiVirtEnd)

      use omp_lib
      use Rimp2_Shared_Data,only: D_B

      implicit double precision(A-H,O-Z)

      logical :: GOPARR,DSKWRK,MASWRK

      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK


      double precision :: AIBJMEM(NOCC,NOCC,NVIR,LddiVirtStart:LddiVirtEnd)
      double precision :: EIG(NBF)

      double precision,allocatable ::BDOV(:,:,:),BDOO(:,:,:),           &
                                     BDVV(:,:,:)

      double precision,allocatable,save :: QOO(:,:),BDV(:,:)
!$omp threadprivate(QOO,BDV)

      integer,allocatable :: ompISE(:,:)

      integer :: D_AIBJMEM



!!!!! PARAMETERS
      NDO=NAUXBASD*NOCC
      NDV=NAUXBASD*NVIR
      NDN=NAUXBASD*(NOCC+NVIR)

!!!!! COPY MO-ENERGY SCALED LAGRANGIAN TO Z-VECTOR
      ALLOCATE(BDOV(NAUXBASD,NOCC,NVIR))
      ALLOCATE(BDOO(NAUXBASD,NOCC,NOCC))
      CALL DDI_GET(D_B, 1,NDO, 1+NOCC,NVIR+NOCC, BDOV)
      CALL DDI_GET(D_B, 1,NDO, 1,NOCC, BDOO)

      Nthreads=omp_get_max_threads()

      AIBJMEM=0.0D00

!$OMP PARALLEL NUM_THREADS(Nthreads)                                    &
!$omp default(none)                                                     &
!$omp shared(LddiVirtStart,LddiVirtEnd)                                 &
!$omp shared(D_B, NDO,NDN, NOCC,NVIR,NORB,NBF, NAUXBASD)                &
!$omp shared(AIBJMEM, BDOV,BDOO,EIG)                                    &
!$omp private(IB)

      ALLOCATE(QOO(NOCC,NOCC))
      ALLOCATE(BDV(NAUXBASD,NVIR))

!$omp do schedule(guided)
      DO IB=LddiVirtStart,LddiVirtEnd

         CALL DDI_GET(D_B, 1+NDO,NDN, IB+NOCC,IB+NOCC, BDV)

         CALL RIMP2_AIBJ_CORE                                           &
             (AIBJMEM(1,1,1,IB), IB,                                    &
              BDOV,BDOO,BDV,EIG,                                        &
              NAUXBASD,NOCC,NVIR,NORB,NBF)
      ENDDO !IB
!$omp end do

      DEALLOCATE(QOO,BDV)

!$OMP END PARALLEL

      END







      SUBROUTINE RIMP2_DLBMEM_AIBJ                                      &
                (AIBJMEM, ddiWSA, EIG, NAUXBASD,NOCC,NVIR,NORB,NBF,     &
                 LddiVirtStart,LddiVirtEnd)

      USE omp_lib
      USE Rimp2_Shared_Data,only: D_B,D_AIBJ_TMP,QOO,QOV

      implicit double precision(A-H,O-Z)

      logical :: GOPARR,DSKWRK,MASWRK,LIB

      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK


      double precision :: AIBJMEM(NOCC,NOCC,NVIR,LddiVirtStart:LddiVirtEnd)
      integer :: ddiWSA(0:NPROC-1,2)

      double precision :: EIG(NBF)

      double precision,allocatable ::BDOV(:,:,:),BDOO(:,:,:),           &
                                     BDVV(:,:,:),AIBJ(:,:,:),BDV(:,:)

      integer,allocatable :: I_DLB_INDEX(:)





!!!!! PARAMETERS
      NDO=NAUXBASD*NOCC
      NDV=NAUXBASD*NVIR
      NDN=NAUXBASD*(NOCC+NVIR)

      N_DLB_BUFF = 1  ! dynamic load balancing buffer
                      ! N_DLB_BUFF = 0 :: completely static work distribution
                      ! N_DLB_BUFF = LddiVirtEnd-LddiVirtStart+1 :: (almost completely) dynamic load balance



      Nthreads=omp_get_max_threads()

      startw=omp_get_wtime()

!!!!! GET B-MATRIX
      ALLOCATE(BDVV(NAUXBASD,NVIR,LddiVirtStart:LddiVirtEnd))
      ALLOCATE(BDOV(NAUXBASD,NOCC,NVIR))
      ALLOCATE(BDOO(NAUXBASD,NOCC,NOCC))

!$OMP PARALLEL NUM_THREADS(Nthreads) &
!$omp shared(BDVV,BDOV,BDOO) &
!$omp shared(LddiVirtStart,LddiVirtEnd) &
!$omp shared(NOCC,NVIR,NAUXBASD) &
!$omp private(IA,IB,IX,II,JJ)

      ! allocate threadprivate
      ALLOCATE(QOO(NOCC,NOCC))
      ALLOCATE(QOV(NOCC,NVIR))
      ! first touch
!$omp do schedule(DYNAMIC)
      DO IB=LddiVirtStart,LddiVirtEnd
         BDVV(1:NAUXBASD,1:NVIR,IB)=0.0D00
      ENDDO
!$omp end do

!$omp do schedule(DYNAMIC)
      DO IA=1,NVIR
         BDOV(1:NAUXBASD,1:NOCC,IA)=0.0D00
      ENDDO
!$omp end do

!$omp do schedule(DYNAMIC)
      DO JJ=1,NOCC
         BDOO(1:NAUXBASD,1:NOCC,JJ)=0.0D00
      ENDDO
!$omp end do
!$OMP END PARALLEL
      ! get B-MATRIX
      CALL DDI_GET(D_B, 1+NDO,NDN, LddiVirtStart+NOCC,LddiVirtEnd+NOCC, BDVV)
      CALL DDI_GET(D_B, 1,NDO, 1+NOCC,NVIR+NOCC, BDOV)
      CALL DDI_GET(D_B, 1,NDO, 1,NOCC, BDOO)

      AIBJMEM=0.0D00

! !/////////// on-node treatment for certain data chunk //////////////////
! !$OMP PARALLEL NUM_THREADS(Nthreads)                                    &
! !$omp default(none)                                                     &
! !$omp shared(LddiVirtStart,LddiVirtEnd, N_DLB_BUFF)                     &
! !$omp shared(D_B, NDO,NDN, NOCC,NVIR,NORB,NBF, NAUXBASD)                &
! !$omp shared(AIBJMEM, BDVV,BDOV,BDOO,EIG)                               &
! !$omp private(IB)
!       ALLOCATE(QOO(NOCC,NOCC))
! !$omp do schedule(guided)
!       ! DO IB=LddiVirtStart,LddiVirtEnd
!       DO IB=LddiVirtStart,LddiVirtEnd - N_DLB_BUFF

!          CALL RIMP2_AIBJ_CORE                                           &
!              (AIBJMEM(1,1,1,IB), IB,                                    &
!               BDOV,BDOO,BDVV(1,1,IB),EIG,                           &
!               NAUXBASD,NOCC,NVIR,NORB,NBF)

!       ENDDO !IB
! !$omp end do
!       ! DEALLOCATE(QOO)
! !$OMP END PARALLEL

      ! DO IB=LddiVirtStart,LddiVirtEnd
      DO IB=LddiVirtStart,LddiVirtEnd - N_DLB_BUFF

         CALL RIMP2OMP_AIBJ_CORE                                        &
             (AIBJMEM(1,1,1,IB), IB,                                    &
              BDOV,BDOO,BDVV(1,1,IB),EIG,                               &
              NAUXBASD,NOCC,NVIR,NORB,NBF)

      ENDDO !IB

#if __bug==1
      endw = omp_get_wtime()
      WRITE(*,'(A40,I5,F10.1)') "TIME AIBJ_STATIC: ",ME,endw-startw
#endif


      startw=omp_get_wtime()
!/////////// dynamic load balancing for the residual calculations //////
      ! DBL buffer array
      ALLOCATE(I_DLB_INDEX(NVIR))
      CALL RIMP2_IDBLINDEX(I_DLB_INDEX,NVIR,N_DLB_BUFF,ddiWSA)

      ! create ddi array for the dynamic buffer
      CALL DDI_CREATE(NOCC*NOCC*NVIR,N_DLB_BUFF*NPROC,D_AIBJ_TMP)

      ALLOCATE(AIBJ(NOCC,NOCC,NVIR))
      ALLOCATE(BDV(NAUXBASD,NVIR))

      ! dynamic load balancing for the buffer
      MINE = -1
      NEXT = -1
      DO IB=1,NVIR
         ! check if IB belong to the static chunk that already processed
         IF(I_DLB_INDEX(IB).EQ.0) CYCLE

         ! dynamic load balancing
         MINE = MINE + 1
         IF(MINE.GT.NEXT) CALL DDI_DLBNEXT(NEXT)
         IF(NEXT.NE.MINE) CYCLE

         ! rank where local data AIBJMEM(1,1,1,IB) resides
         IRANK=NPROC-1
         DO IPROC=0,NPROC-2
            IF(IB.GE.ddiWSA(IPROC,2)-N_DLB_BUFF .AND. IB.LT.ddiWSA(IPROC+1,1)) THEN
               IRANK=IPROC
               EXIT
            ENDIF
         ENDDO

         ! if IRANK .eq. ME => we can take advantage of local data
         IF(ME.EQ.IRANK) THEN
            CALL RIMP2OMP_AIBJ_CORE &
                (AIBJMEM(1,1,1,IB), IB, &
                 BDOV,BDOO,BDVV(1,1,IB),EIG, &
                 NAUXBASD,NOCC,NVIR,NORB,NBF)
            I_DLB_INDEX(IB) = 0
         ! or else boring get and put
         ELSE
            CALL DDI_GET(D_B, 1+NDO,NDN, IB+NOCC,IB+NOCC, BDV)
            CALL RIMP2OMP_AIBJ_CORE &
                (AIBJ, IB, &
                 BDOV,BDOO,BDV,EIG, &
                 NAUXBASD,NOCC,NVIR,NORB,NBF)
            CALL DDI_PUT(D_AIBJ_TMP, 1, NOCC*NOCC*NVIR, I_DLB_INDEX(IB),I_DLB_INDEX(IB), AIBJ)
         ENDIF

      ENDDO !IB
      CALL DDI_DLBRESET

      DEALLOCATE(BDV)

      ! fetch DDI data to AIBJMEM
      DO IB=LddiVirtEnd-N_DLB_BUFF,LddiVirtEnd
         IF(I_DLB_INDEX(IB).EQ.0) CYCLE
         CALL DDI_GET(D_AIBJ_TMP, 1,NOCC*NOCC*NVIR, I_DLB_INDEX(IB),I_DLB_INDEX(IB), AIBJMEM(1,1,1,IB))
      ENDDO

      ! to avoid synchronization,
      ! the D_AIBJ_TMP is destroy
      ! later when the opportunity comes

!$OMP PARALLEL NUM_THREADS(Nthreads)
      DEALLOCATE(QOO,QOV)
!$OMP END PARALLEL

#if __bug==1
      endw = omp_get_wtime()
      WRITE(*,'(A40,I5,F10.1)') "TIME AIBJ_DYNAMIC: ",ME, endw-startw
#endif

      END





      SUBROUTINE RIMP2_IDBLINDEX(I_DLB_INDEX,NLOOP,N_DLB_BUFF,ddiWSA)

      ! the dynamic load balancing buffer array

      implicit double precision(A-H,O-Z)

      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK

      integer :: I_DLB_INDEX(NLOOP)
      integer :: ddiWSA(0:NPROC-1,2)
      logical :: LIB
      logical :: GOPARR,DSKWRK,MASWRK

      I_DLB_INDEX=0
      ICOUNT=1
      DO IB=1,NLOOP
         ! check if IB belong to static part
         LIB=.FALSE.
         DO IPROC=0,NPROC-1
            IF(IB.GE.ddiWSA(IPROC,1) .AND. IB.LE.ddiWSA(IPROC,2)-N_DLB_BUFF) THEN
               LIB=.TRUE.
               EXIT
            ENDIF
         ENDDO
         IF(LIB) CYCLE
         I_DLB_INDEX(IB)=ICOUNT
         ICOUNT=ICOUNT+1
      ENDDO

      END




      SUBROUTINE RIMP2_AIBJ_CORE                                        &
                (AIBJ, IB, BDOV,BDOO,BDV,EIG,                           &
                 NAUXBASD,NOCC,NVIR,NORB,NBF)

      use omp_lib
      use Rimp2_Shared_Data,only: D_B
      USE Rimp2_Shared_Data,only: QOO

      implicit double precision(A-H,O-Z)

      logical :: GOPARR,DSKWRK,MASWRK

      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK


      double precision :: AIBJ(NOCC,NOCC,NVIR)

      double precision :: EIG(NBF)
      double precision :: BDOV(NAUXBASD,NOCC,NVIR)
      double precision :: BDOO(NAUXBASD,NOCC,NOCC)
      double precision :: BDV(NAUXBASD,NVIR)
      ! double precision :: QOO(NOCC,NOCC)


      eb=EIG(IB+NOCC)

      DO IA=1,NVIR
!!!!! AIJB = 4*(ai|bj) - (aj|bi)
         CALL DGEMM('T','N', NOCC,NOCC,NAUXBASD,                        &
                     1.0D00, BDOV(1,1,IB),NAUXBASD,                     &
                             BDOV(1,1,IA),NAUXBASD,                     &
                     0.0D00, QOO,NOCC)
         FORALL(J=1:NOCC,I=1:NOCC)
            AIBJ(I,J,IA)=4.0D00*QOO(I,J)-QOO(J,I)
         ENDFORALL
!!!!! AIJB <- -(ab|ij)
         CALL DGEMM('T','N', NOCC*NOCC,1, NAUXBASD,                     &
                    -1.0D00, BDOO,NAUXBASD,                             &
                             BDV(1,IA),NAUXBASD,                        &
                     1.0D00, AIBJ(1,1,IA),NOCC*NOCC)
!!!!! AIJB = AIBJ/[e(b)-e(i)]
         FORALL(J=1:NOCC,I=1:NOCC)
            AIBJ(I,J,IA)=AIBJ(I,J,IA)/(eb-EIG(I))
         ENDFORALL
      ENDDO !IA

!!!!! AIJB(I,I) = AIJB(I,I) + 1.0D00
      DO II=1,NOCC
         AIBJ(II,II,IB)=AIBJ(II,II,IB)+1.0D00
      ENDDO

      END








      SUBROUTINE RIMP2OMP_AIBJ_CORE                                     &
                (AIBJ, IB,BDOV,BDOO,BDV,EIG,NAUXBASD,NOCC,NVIR,NORB,NBF)

      USE omp_lib
      USE Rimp2_Shared_Data,only: D_B
      USE Rimp2_Shared_Data,only: QOO,QOV

      implicit double precision(A-H,O-Z)

      logical :: GOPARR,DSKWRK,MASWRK

      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK


      double precision :: AIBJ(NOCC,NOCC,NVIR)

      double precision :: EIG(NBF)
      double precision :: BDOV(NAUXBASD,NOCC,NVIR)
      double precision :: BDOO(NAUXBASD,NOCC,NOCC)
      double precision :: BDV(NAUXBASD,NVIR)

!       double precision,allocatable,save :: QOV(:,:)
! !$omp threadprivate(QOV)

      eb=EIG(IB+NOCC)
      Nthreads=omp_get_max_threads()
!$OMP PARALLEL NUM_THREADS(Nthreads)                                    &
!$omp default(none)                                                     &
!$omp shared(IB,NOCC,NVIR,NAUXBASD)                                     &
!$omp shared(AIBJ,BDOV,BDOO,BDV,EIG,eb)                                 &
!$omp private(IA,II,I,J)


      ! AIJB = 4*(ai|bj) - (aj|bi)
!$omp do schedule(dynamic)
      DO IA=1,NVIR
         CALL DGEMM('T','N', NOCC,NOCC,NAUXBASD,                        &
                     1.0D00, BDOV(1,1,IB),NAUXBASD,                     &
                             BDOV(1,1,IA),NAUXBASD,                     &
                     0.0D00, QOO,NOCC)

         FORALL(J=1:NOCC,I=1:NOCC)
            AIBJ(I,J,IA)=4.0D00*QOO(I,J)-QOO(J,I)
         ENDFORALL

      ENDDO !IA
!$omp end do


      ! AIJB <- -(ab|ij)
!$omp do schedule(dynamic)
      DO II=1,NOCC
         CALL DGEMM('T','N', NOCC,NVIR, NAUXBASD,                       &
                    -1.0D00, BDOO(1,1,II),NAUXBASD,                     &
                             BDV,NAUXBASD,                              &
                     0.0D00, QOV,NOCC)

         FORALL(J=1:NOCC,IA=1:NVIR)
            AIBJ(J,II,IA)=AIBJ(J,II,IA)+QOV(J,IA)
         ENDFORALL
      ENDDO !IA
!$omp end do ! <<< implicit barrier


      ! AIJB = AIBJ/[e(b)-e(i)]
!$omp do schedule(dynamic)
      DO IA=1,NVIR
         FORALL(J=1:NOCC,I=1:NOCC)
            AIBJ(I,J,IA)=AIBJ(I,J,IA)/(eb-EIG(I))
         ENDFORALL
      ENDDO
!$omp end do

      ! AIBJ + delta(ij,ab)
!$omp do schedule(dynamic)
      DO II=1,NOCC
         AIBJ(II,II,IB)=AIBJ(II,II,IB)+1.0D00
      ENDDO
!$omp end do

!$OMP END PARALLEL



      END








      SUBROUTINE RIMP2_INIT_ZVEC                                        &
                (GMVEC, AIBJMEM,XLOV,EIG,                               &
                 NOCC,NVIR,NBF,ITERMX,LddiVirtStart,LddiVirtEnd)

      use omp_lib

      implicit double precision(a-h,o-z)

      double precision :: GMVEC(NOCC,NVIR,0:ITERMX+1)
      double precision :: XLOV(NOCC,NVIR)
      double precision :: EIG(NBF)
      double precision :: AIBJMEM(NOCC,NOCC,NVIR,LddiVirtStart:LddiVirtEnd)


      GMVEC  = 0.0D00
      FORALL(IA=1:NVIR,II=1:NOCC)
         XLOV(II,IA)=XLOV(II,IA)/(EIG(IA+NOCC)-EIG(II))
      ENDFORALL
      CALL DCOPY(NOCC*NVIR,XLOV,1,GMVEC(1,1,0),1)

      Nthreads = omp_get_max_threads()

!$OMP PARALLEL DO NUM_THREADS(Nthreads)                                 &
!$omp shared(LddiVirtStart,LddiVirtEnd)                                 &
!$omp shared(NOCC,NVIR)                                                 &
!$omp shared(XLOV,GMVEC,AIBJMEM)                                        &
!$omp private(IB)

      DO IB=LddiVirtStart,LddiVirtEnd
         CALL DCOPY(NOCC,XLOV(1,IB),1,GMVEC(1,IB,1),1)
         CALL DGEMV('N',NOCC,NOCC*NVIR,                                 &
                    -1.0D00, AIBJMEM(1,1,1,IB),NOCC,                    &
                             GMVEC(1,1,0),1,                            &
                     1.0D00, GMVEC(1,IB,1),1)
      ENDDO
!$OMP END PARALLEL DO

      CALL DDI_GSUMF(1001,GMVEC(1,1,1),NOCC*NVIR)

      END








      SUBROUTINE RIMP2_SYMDIV (LddiSymStart,LddiSymEnd, IRANK,NPROC,NN)


      KK=NN*(NN+1)/(NPROC*2)

      LddiSymStart=0
      LddiSymEnd=0

      DO II=0,IRANK
         Nchunks=0
         IStart=LddiSymEnd+1
         LddiSymStart=IStart
         DO IPOINT=IStart,NN
               Nchunks=Nchunks+IPOINT
               LddiSymEnd=IPOINT
               IF(Nchunks .GE. KK) EXIT
         ENDDO !IPOINT
      ENDDO

      IF(IRANK .eq. NPROC-1) LddiSymEnd=NN

      END





      SUBROUTINE RIMP2_AIBJ(AOO,QOO,BDOO,BD,EIG,IB,NOCC,NAUXBASD,NBF)

      double precision :: AOO(NOCC,NOCC)
      double precision :: QOO(NOCC,NOCC)
      double precision :: BDOO(NAUXBASD,NOCC,NOCC)
      double precision :: BD(NAUXBASD)
      double precision :: EIG(NBF)

            AOO=4.0D00*QOO-transpose(QOO)

!          AIJB <- -(ab|ij)
            CALL DGEMM('T','N', NOCC*NOCC,1, NAUXBASD,                  &
                       -1.0D00, BDOO,NAUXBASD,                          &
                                 BD,NAUXBASD,                           &
                        1.0D00, AOO,NOCC*NOCC)
!          AIJB = AOO/[e()(i)]
            FORALL(I=1:NOCC,J=1:NOCC)
               AOO(I,J)=AOO(I,J)/(EIG(IB+NOCC)-EIG(I))
            ENDFORALL
      END



!>*module rimp2grd   *deck RIMP2_FORM_VDX
!>
!>     @brief   form matrix V
!>     @author  buu
!>     @date    sept 8, 2018
!>
!>     @detail  1) form matrix of 2-2ERIs
!>              2) decompose
!>              3) take inverse
!>
      SUBROUTINE RIMP2_FORM_VDX                                         &
                 (VXX, STOL,VTOL,NAUXBAS,NAUXBASD,IVMTD,OTHAUX)

      use Rimp2_Shared_Data

      implicit double precision(A-H,O-Z)

      COMMON /IOFILE/ IR,IW,IP,IJK,IPK,IDAF,NAV,IODA(950)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      DATA CHECK    /8HCHECK   /
      logical :: GOPARR,DSKWRK,MASWRK
      logical :: OTHAUX

!!!!! OUTPUT
      double precision VXX(NAUXBAS,NAUXBAS)


                          !!! ~ GO ~ !!!


!!!!! PARAMETERS
      NAUXBAS1 = NAUXBAS


! !!!!! DECOMPOSE VXX TO VDX
!       IF(IVMTD .EQ. 2) THEN
! !!!!! EIGENVALUE METHOD
!          CALL RIMP2_EIGEN_DECOMPOSE_VXX (VXX, VTOL,NAUXBAS1,NAUXBASD)
!       ELSE
! !!!!! CHOLESKY METHOD
!          CALL RIMP2_CHOLESKY_DECOMPOSE_VXX (VXX, NAUXBAS1)
!          NAUXBASD = NAUXBAS1
!       END IF


!!!!! FORM VXX
      CALL RIMP2_X_VXX(VXX, NAUXBAS)


!!!!! CHOLESKY DECOMPOSITION OF VXX
      CALL RIMP2_CHOLESKY_DECOMPOSE_VXX (VXX, NAUXBAS1)

      NAUXBASD = NAUXBAS1

      RETURN

      END !***********************************************************





      SUBROUTINE RIMP2_X_VXX(VXX, NAUXBAS)
      use Rimp2_Shared_Data
      use omp_lib
      implicit double PRECISION(A-H,O-Z)


      COMMON /AUXBAS/ EXAUX(MXAXGTOT),CAUXS(MXAXGTOT),CAUXP(MXAXGTOT),  &
                      CAUXD(MXAXGTOT),CAUXF(MXAXGTOT),CAUXG(MXAXGTOT),  &
                      CAUXH(MXAXGTOT),CAUXI(MXAXGTOT),                  &
                      KAUXST(MXAUXSH),KAUXAT(MXAUXSH),KAUXTY(MXAUXSH),  &
                      KAUXNG(MXAUXSH),KAUXLO(MXAUXSH),KAUXMI(MXAUXSH),  &
                      KAUXMX(MXAUXSH),NAUXSH
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,                    &
                      ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /OUTPUT/ NPRINT,ITOL,ICUT,NORMF,NORMP,NOPK

!!!!! OUTPUT
      double precision :: VXX(NAUXBAS,NAUXBAS)

!!!!! LOCAL
      integer,allocatable :: ompISE(:,:)
      double precision,allocatable :: TI(:),TJ(:),W(:),GHONDO(:)
      logical :: GOPARR,DSKWRK,MASWRK, IANDJ,LOMP
      double precision :: CC(3,2)


                          !!! ~ GO ~ !!!



!!!!! MAX BASIS FUNCTIONS
      MaxG=MaxAUXBFX*MaxAUXBFX

      ALLOCATE(TI(MaxG))
      ALLOCATE(TJ(MaxG))
      ALLOCATE(W(MaxG))
      ALLOCATE(GHONDO(MaxG))


!!!!! DDI WORK DISTRIBUTION NAUXSH LOOP
      ALLOCATE(ompISE(0:NPROC-1,2))
      call rimp2_WSA(ompISE,NPROC,1,NAUXSH,LOMP)
      LddiAuxShStart=ompISE(ME,1)
      LddiAuxShEnd=ompISE(ME,2)
      DEALLOCATE(ompISE)


!!!!! NUMBER OF THREADS
      Nthreads_env=omp_get_max_threads()
      Nthreads=MIN0(Nthreads_env,LddiAuxShEnd-LddiAuxShStart-1)

#if __bug==1
      WRITE(*,'(A40,2I5)') "NUM_THREAD SUB RIMP2_X_VXX", Nthreads, Nthreads_env
#endif

!!!!! OMP WORK DISTRIBUTION LddiAuxShStart:LddiAuxShEnd LOOP
      ALLOCATE(ompISE(0:Nthreads-1,2))
      call rimp2_WSA(ompISE,Nthreads,LddiAuxShStart,LddiAuxShEnd,LOMP)

      VXX=0.0D00

!!!!! FORM VXX
!$omp parallel NUM_THREADS(Nthreads)                                    &
!$omp shared(ompISE)                                                    &
!$omp shared(KAUXAT,C,KAUXTY,KAUXNG,KAUXMI,KAUXMX,KAUXST,KAUXLO)        &
!$omp shared(FlgSphAux,LOCSPH)                                          &
!$omp shared(EXAUX,CAUXI,CAUXH,CAUXG,CAUXF,CAUXD,CAUXP,CAUXS,NORM)      &
!$omp shared(NAUXBAS)                                                   &
!$omp shared(VXX)                                                       &
!$omp private(IthreadID,II,Istart,Iend)                                 &
!$omp private(IAT,CC,LANGI,NGI,MINI,MAXI,I1,LOCI,NI)                    &
!$omp private(MINIS,MAXIS,LOCIS,NIS)                                    &
!$omp private(JJ,JAT,LANGJ,NGJ,MINJ,MAXJ,J1,LOCJ,NJ,IANDJ,TI,TJ)        &
!$omp private(MINJS,MAXJS,LOCJS,NJS)                                    &
!$omp private(GHONDO,W)

      IthreadID=omp_get_thread_num()
      LompAuxShStart=ompISE(IthreadID,1)
      LompAuxShEnd=ompISE(IthreadID,2)

      ! DO II = 1, NAUXSH
      ! DO II=LddiAuxShStart,LddiAuxShEnd
      DO II=LompAuxShStart,LompAuxShEnd
         IAT = KAUXAT(II)
         CC(1,1) = C(1,IAT)
         CC(2,1) = C(2,IAT)
         CC(3,1) = C(3,IAT)
         LANGI = KAUXTY(II)-1
         NGI   = KAUXNG(II)
         MINI  = KAUXMI(II)
         MAXI  = KAUXMX(II)
         I1    = KAUXST(II)
         LOCI  = KAUXLO(II)-MINI
         NI    = MAXI-MINI+1

         IF(FlgSphAux) THEN
            CALL RIMP2CSTRM(TI,LANGI,MINI,NI,MINIS,MAXIS)
            LOCIS = LOCSPH(II)-MINIS
            NIS   = MAXIS-MINIS+1
         END IF

         DO JJ = 1, II

            JAT = KAUXAT(JJ)
            CC(1,2) = C(1,JAT)
            CC(2,2) = C(2,JAT)
            CC(3,2) = C(3,JAT)
            LANGJ = KAUXTY(JJ)-1
            NGJ   = KAUXNG(JJ)
            MINJ  = KAUXMI(JJ)
            MAXJ  = KAUXMX(JJ)
            J1    = KAUXST(JJ)
            LOCJ  = KAUXLO(JJ)-MINJ
            NJ = MAXJ-MINJ+1
            IANDJ = II .EQ. JJ

            IF(FlgSphAux) THEN
               CALL RIMP2CSTRM(TJ,LANGJ,MINJ,NJ,MINJS,MAXJS)
               LOCJS = LOCSPH(JJ)-MINJS
               NJS   = MAXJS-MINJS+1
            END IF


            CALL RIMP2_X_I22                                            &
                (GHONDO, II,JJ,CC,                                      &
                 EXAUX(I1),EXAUX(J1),                                   &
                 LANGI,LANGJ,NI,NJ,NGI,NGJ,MINI,MAXI,MINJ,MAXJ)

            CALL RIMP2_VXX                                              &
                (VXX, GHONDO, TI,TJ,W,NAUXBAS,                          &
                 LOCI,LOCJ,MINI,MAXI,MINJ,MAXJ,NI,NJ,                   &
                 LOCIS,LOCJS,MINIS,MAXIS,MINJS,MAXJS,NIS,NJS,           &
                 FlgSphAux,IANDJ)

         END DO
      END DO
!$OMP END PARALLEL

      CALL DDI_GSUMF(1504,VXX,NAUXBAS*NAUXBAS)

!!!!! CLEAN UP
      DEALLOCATE(TI,TJ,W,GHONDO,ompISE)

      RETURN

      END !***********************************************************














!
      SUBROUTINE RIMP2_X_I22                                            &
                (GG, II,JJ,CC,                                          &
                 EXI,EXJ,                                               &
                 LI,LJ,NI,NJ,NGI,NGJ,MINI,MAXI,MINJ,MAXJ)
!
      use Rimp2_Shared_Data
!
      implicit double precision(a-h,o-z)
!
      COMMON /ROOT  / XX,U(13),W(13),NROOTS
!$omp threadprivate(/ROOT  /)
!
      double precision GG(NI*NJ)
      double precision CC(3,2)
      double precision EXI(NGI),EXJ(NGJ)
      double precision DI(28),DJ(28),G(1911)
      double precision A1(3,13),B1(13),C1(13),A2(3,13),B2(13),F00(13)
!
      integer LLI(3,28),LLJ(3,28)


!    TOTAL ANG MOMENTUM, AND NUM OF ROOTS
      LTOT = LI+LJ
      NROOTS = LTOT/2+1
!
      RXIJ = CC(1,1)-CC(1,2)
      RYIJ = CC(2,1)-CC(2,2)
      RZIJ = CC(3,1)-CC(3,2)
      RIJ2 = RXIJ**2+RYIJ**2+RZIJ**2

!    ANG MOMENTUM COMPONENTS
      DO I = MINI, MAXI
         ICNT=(I-MINI)+1
         LLI(1,ICNT) = LX(I)
         LLI(2,ICNT) = LY(I)
         LLI(3,ICNT) = LZ(I)
      END DO
!
      DO J = MINJ, MAXJ
         JCNT=(J-MINJ)+1
         LLJ(1,JCNT) = LX(J)
         LLJ(2,JCNT) = LY(J)
         LLJ(3,JCNT) = LZ(J)
      END DO

!    ZERO OUT I22
      GG=0.0D00

!    I-SHELL
      DO IG = 1, NGI

         AI  = EXI(IG)
         AIINV = 1.0D00/AI

!       CONTRACTION COEFFICIENTS
         DO I=1,NI
            DI(I)=DAux(II,IG,I)
         ENDDO

!       J-SHELL
         DO JG = 1, NGJ
            AJ  = EXJ(JG)
            AJINV = 1.0D00/AJ

!          CONTRACTION COEFFICIENTS
            DO J=1,NJ
              DJ(J)=DAux(JJ,JG,J)
            ENDDO

            RHO = AI*AJ/(AI+AJ)
            XX = RHO*RIJ2

!          RYS ROOTS AND WEIGHTS
            SELECT CASE(NROOTS)
               CASE(:3)
                  CALL RT123
               CASE(4)
                  CALL ROOT4
               CASE(5)
                  CALL ROOT5
               CASE(6:)
                  CALL ROOT6
            END SELECT

!          IXYZ
            CALL I22XYZ(G, CC,U,W,A1,B1,C1,                             &
                 A2,B2,F00,AI,AJ,LI,LJ,NROOTS)

!          FORM I22
            CALL RIMP2_I22(GG, G, DI,DJ,LLI,LLJ,                        &
                 LI,LJ,NI,NJ,NROOTS)
         END DO
      END DO

      RETURN

      END !***********************************************************




!
      SUBROUTINE I22XYZ                                                 &
                (G, CC,U,W,                                             &
                 A1,B1,C1, A2,B2, F00,                                  &
                 AI,AJ, LI,LJ, NROOTS)

!
      implicit double precision(a-h,o-z)

!    OUTPUT
      double precision :: G(3,NROOTS,0:LI,0:LJ)

!    INPUT
      double precision :: CC(3,2),U(NROOTS),W(NROOTS)
      double precision :: A1(3,NROOTS),B1(NROOTS),C1(NROOTS)
      double precision :: A2(3,NROOTS),B2(NROOTS),F00(NROOTS)


!    PREFACTORS OF I22

      AIJINV = 1.0D00/(AI+AJ)
      RHO = AI*AJ*AIJINV

      AIINV  = 1.0D00/AI
      AJINV  = 1.0D00/AJ

      AIINV2 = 0.5D00*AIINV
      AJINV2 = 0.5D00*AJINV

      AIJINV2 = 0.5D00*AIJINV

      AIJINVSQ = SQRT(AIJINV)

      RXIJ = CC(1,2)-CC(1,1)
      RYIJ = CC(2,2)-CC(2,1)
      RZIJ = CC(3,2)-CC(3,1)

      DO IR = 1, NROOTS

         U2 = U(IR)*RHO
         T2 = U2/(RHO+U2)

         VAL1 = RHO*AIINV*T2
         VAL2 = RHO*AJINV*T2
         VAL3 = AIINV2-AIINV2*VAL1
         VAL4 = AJINV2-AJINV2*VAL2

         A1(1,IR) = VAL1*RXIJ
         A1(2,IR) = VAL1*RYIJ
         A1(3,IR) = VAL1*RZIJ

         A2(1,IR) = -VAL2*RXIJ
         A2(2,IR) = -VAL2*RYIJ
         A2(3,IR) = -VAL2*RZIJ

         B1(IR) = VAL3
         B2(IR) = VAL4

         C1(IR) = AIJINV2*T2

         F00(IR) = AIJINVSQ*W(IR)

      END DO

!    (s|s)
      DO IR = 1, NROOTS
         G(1,IR,0,0) = 1.0D00
         G(2,IR,0,0) = 1.0D00
         G(3,IR,0,0) = F00(IR)
      END DO

!    (p|s); (d|s), (f|s), (g|s)
      SELECT CASE(LI)
         CASE(1)
            CALL L10(G, LI,LJ, NROOTS,A1,F00)
         CASE(2:)
            CALL L20(G, LI,LJ, NROOTS,A1,B1,F00)
      END SELECT

      IF (LJ==0) RETURN

!    (s|p); (p|p); (d|p),(f|p),(g|p)
      SELECT CASE(LI)
         CASE(0)
            CALL L01(G, LI,LJ, NROOTS,A2,F00)
         CASE(1)
            CALL L11(G, LI,LJ, NROOTS,A1,A2,C1,F00)
         CASE(2:)
            CALL L21(G, LI,LJ, NROOTS,A1,A2,B1,C1,F00)
      END SELECT

      IF(LJ==1) RETURN

!    (s|d),(s|f),(s|g);
!    (p|d),(p|f),(p|g);
!    (d|d),(d|f),(d|g); (f|d),(f|f),(f|g); (g|d),(g|f),(g|g)

      SELECT CASE(LI)
         CASE(0)
            CALL L02(G, LI,LJ, NROOTS,A2,B2)
         CASE(1)
            CALL L12(G, LI,LJ, NROOTS,A1,A2,B2,C1)
         CASE(2:)
            CALL L22(G, LI,LJ, NROOTS,A1,A2,B1,B2,C1)
      END SELECT

      RETURN

      END !***********************************************************


!
      SUBROUTINE RIMP2_I22(GG, G, DI,DJ,LLI,LLJ,LI,LJ,NI,NJ,NROOTS)

      use Rimp2_Shared_Data

      implicit double PRECISION(A-H,O-Z)

      double precision GG(NI*NJ)

      double precision G(3,NROOTS,0:LI,0:LJ)
      double precision DI(NI),DJ(NJ)

      integer LLI(3,NI),LLJ(3,NJ)


      IJ = 0
      DO J = 1, NJ
         NJX = LLJ(1,J)
         NJY = LLJ(2,J)
         NJZ = LLJ(3,J)

         FACTOR=PI252*DJ(J)

         DO I = 1, NI
            NIX = LLI(1,I)
            NIY = LLI(2,I)
            NIZ = LLI(3,I)

            VAL = 0.0D00
            DO IR = 1, NROOTS
               VAL = VAL+G(1,IR,NIX,NJX)                                &
                        *G(2,IR,NIY,NJY)                                &
                        *G(3,IR,NIZ,NJZ)
            END DO

            IJ = IJ+1
            GG(IJ) = GG(IJ) + FACTOR*DI(I)*VAL

         END DO
      END DO

      RETURN

      END !***********************************************************


!
      SUBROUTINE RIMP2_VXX                                              &
          (VXX, GHONDO, TI,TJ,W,NAUXBAS,                                &
           LOCI,LOCJ,MINI,MAXI,MINJ,MAXJ,NI,NJ,                         &
           LOCIS,LOCJS,MINIS,MAXIS,MINJS,MAXJS,NIS,NJS,SPHAUX,IANDJ)

      implicit double PRECISION(A-H,O-Z)

      double precision VXX(NAUXBAS,NAUXBAS)

      double precision GHONDO(NI*NJ)
      double precision TI(NI*NI),TJ(NJ*NJ),W(NI*NJ)

      logical SPHAUX,IANDJ


!    CARTESIAN -> SPHERICAL
      IF(SPHAUX) THEN
         MINIT = MINIS
         MAXIT = MAXIS

         MINJT = MINJS
         MAXJT = MAXJS

         LOCIT = LOCIS
         LOCJT = LOCJS

         NIT = NIS

         CALL DGEMM('N','N', NI,NJS,NJ,      &
                     1.0D00, GHONDO,NI,      &
                              TJ,NJ,         &
                     0.0D00, W,NI)

         CALL DGEMM('T','N',  NIS,NJS,NI,    &
                     1.0D00, TI,NI,          &
                              W,NI,          &
                     0.0D00, GHONDO,NIS)

      ELSE
         MINIT = MINI
         MAXIT = MAXI

         MINJT = MINJ
         MAXJT = MAXJ

         LOCIT = LOCI
         LOCJT = LOCJ

         NIT = NI
      END IF

!    VXX <- I22
      JMAX = MAXJT

      DO I = MINIT, MAXIT
         ICNT=(I-MINIT)+1
         IAUX = LOCIT+I

         IF(IANDJ) JMAX = I
         DO J = MINJT, JMAX
            JCNT=(J-MINJT)+1
            JAUX = LOCJT+J
            IJ = ICNT+(JCNT-1)*NIT
            VXX(IAUX,JAUX) = GHONDO(IJ)
            VXX(JAUX,IAUX) = GHONDO(IJ)

         END DO
      END DO

      RETURN

      END !***********************************************************



!!!!! >> SAVED FOR LATER USE
      SUBROUTINE RIMP2_EIGEN_DECOMPOSE_VXX (V, VTOL,NAUXBAS,NAUXBASD)
      use Rimp2_Shared_Data, only : zero, one
      implicit double PRECISION(A-H,O-Z)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      double precision V(NAUXBAS,NAUXBAS)
      logical :: GOPARR,DSKWRK,MASWRK
      double precision,allocatable :: WRK1(:),TRI(:),WRK2(:),IWRK(:),SCR(:,:)
      ALLOCATE(WRK1(NAUXBAS*NAUXBAS))
      ALLOCATE(TRI((NAUXBAS*NAUXBAS+NAUXBAS)/2))
      ALLOCATE(WRK2(NAUXBAS))
      ALLOCATE(IWRK(NAUXBAS))
      ALLOCATE(SCR(NAUXBAS,8))
!!!!! DIAGONALIZE V MATRIX: VECTORS TO -V- AND EIGENVALUES TO -WRK2-
      CALL CPYSQT(V,TRI,NAUXBAS,1)
      CALL GLDIAG(NAUXBAS,NAUXBAS,NAUXBAS,TRI,SCR,WRK2,V,IERR,IWRK)
!!!!! MAKE V^(-1/2) MATRIX
      NAUXBASD = 0
      DO IAUX = NAUXBAS, 1, -1
         IF(WRK2(IAUX) .LT. VTOL) GO TO 100
         NAUXBASD = NAUXBASD+1
         WRK1(NAUXBASD) = WRK2(IAUX)
      END DO
 100  CONTINUE
      IVAL = NAUXBAS-NAUXBASD
      IF(MASWRK .AND. IVAL .GT. 0) THEN
         ICNT = 0
         DO I = 1, IVAL
            IF(WRK2(I) .LT. ZERO) ICNT = ICNT+1
         END DO
      END IF
      DO IAUX = 1, NAUXBASD
         WRK2(IAUX) = ONE/SQRT(WRK1(IAUX))
      END DO
      ICNT = 0
      DO JAUX = 1, NAUXBAS
         IAUX2 = NAUXBAS
         DO IAUX = 1, NAUXBASD
            ICNT = ICNT+1
            WRK1(ICNT) = V(JAUX,IAUX2)*WRK2(IAUX)
            IAUX2 = IAUX2-1
         END DO
      END DO
      CALL DCOPY(NAUXBAS*NAUXBASD,WRK1,1,V,1)
      DEALLOCATE(WRK1,TRI,WRK2,IWRK,SCR)
      RETURN
      END !***********************************************************



      SUBROUTINE RIMP2_CHOLESKY_DECOMPOSE_VXX(V,NAUXBAS)
      implicit double precision(A-H,O-Z)

      double precision :: V(NAUXBAS,NAUXBAS)


!!!!! CHOLESKY DECOMPOSION OF 2C-2E MATRIX V=L*LT
      CALL DPOTRF('L',NAUXBAS,V,NAUXBAS,INFO)


!!!!! DETERMINATION OF INVERSE OF CHOLESKY DECOMPOSED MATRIX L^(-1)
      CALL DTRTRI('L','N',NAUXBAS,V,NAUXBAS,INFO)


!!!!! ZERO OUT LOWER PART
      IF(NAUXBAS .GT. 1) THEN
         DO JAUX = 2, NAUXBAS
            DO IAUX = 1, JAUX-1
               V(IAUX,JAUX) = 0.0D00
            END DO
         END DO
      END IF

      RETURN

      END !***********************************************************


!>*module rimp2grd   *deck RIMP2_SYM_TRANS_MO
!>
!>     @brief   read MO coefficients
!>     @author  buu
!>     @date    sept 8, 2018
!>
!>     @detail  read MO coefficients
!>
      SUBROUTINE RIMP2_SYM_TRANS_MO(VEC,NBASIS,NACORE,NORBA)

      use Rimp2_Shared_Data, only : MXATM, MXAO
      use omp_lib

      implicit double precision(a-h,o-z)

      COMMON /IOFILE/ IR,IW,IP,IJK,IPK,IDAF,NAV,IODA(950)
      COMMON /IJPAIR/ IA(MXAO)
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,                    &
                      ZAN(MXATM),C(3,MXATM),IAN(MXATM)

      double precision VEC(NBASIS,NBASIS)

      double precision,allocatable :: MOLAB(:),MOIRP(:),                &
                              MODEG(:),Q(:),S(:),WRK(:)


!!!!! PARAMETERS
      NBF2 = NBASIS*(NBASIS+1)/2
      NBF3 = NBASIS*NBASIS
      NORB = NORBA+NACORE

!!!!! ALLOCATE ARRAYS
      ALLOCATE(S(NBF2))
      ALLOCATE(Q(NBF3))
      ALLOCATE(MOLAB(NBASIS))
      ALLOCATE(MOIRP(NBASIS))
      ALLOCATE(MODEG(NBASIS))
      ALLOCATE(WRK(NBASIS))


      CALL DAREAD(IDAF,IODA,S,NBF2,12,0)
      CALL DAREAD(IDAF,IODA,Q,NBF3,45,0)

!!!!! SYMMETRY TRANSFORM MO
      CALL TRFSYM(MOLAB,MOIRP,MODEG,Q,S,VEC,WRK,IA,NQMT,                &
                  NBASIS,NORB,NBASIS)

      END !***********************************************************



!>*module rimp2grd   *deck RIMP2_X_FORM_MO_I32
!>
!>     @brief   form 3-2ERIs
!>     @author  buu
!>     @date    sept 8, 2018
!>
!>     @detail  form 3-2ERIs using Rys Quadrature method
!>
      SUBROUTINE RIMP2_X_FORM_MO_I32 (VEC,NBASIS,NOCC,NVIR,NAUXBAS)

      use Rimp2_Shared_Data
      use omp_lib

      implicit double precision(a-h,o-z)

      COMMON /AUXBAS/ EXAUX(MXAXGTOT),CAUXS(MXAXGTOT),CAUXP(MXAXGTOT),  &
                      CAUXD(MXAXGTOT),CAUXF(MXAXGTOT),CAUXG(MXAXGTOT),  &
                      CAUXH(MXAXGTOT),CAUXI(MXAXGTOT),                  &
                      KAUXST(MXAUXSH),KAUXAT(MXAUXSH),KAUXTY(MXAUXSH),  &
                      KAUXNG(MXAUXSH),KAUXLO(MXAUXSH),KAUXMI(MXAUXSH),  &
                      KAUXMX(MXAUXSH),NAUXSH
      COMMON /NSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),      &
                      CF(MXGTOT),CG(MXGTOT),CH(MXGTOT),CI(MXGTOT),      &
                      KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),   &
                      KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,                    &
                      ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /OUTPUT/ NPRINT,ITOL,ICUT,NORMF,NORMP,NOPK


!!!!! INPUT
      double precision :: VEC(NBASIS,NBASIS)


!!!!! LOCAL
      double precision,allocatable :: T1(:),GG(:),TL(:),XX1(:),         &
                                      XX2(:), I32(:,:,:)
      integer,allocatable :: ddiWS(:,:)
      logical :: GOPARR,DSKWRK,MASWRK,KANDI, IANDJ,LOMP
      double precision :: CC(3,3)





!!!!! PARAMETERS
      MXBSH=MaxAUXBFX
      IF(MaxAUXANG < MaxATMANG) MXBSH=MaxATMBFX
      MaxT12=MAX(NBASIS*NBASIS, MXBSH*MXBSH*MXBSH)
      ISIZE = MAXVAL(KAUXNG)


!!!!! DDI WORK DISTRIBUTION :: NAUXSH LOOP
      ALLOCATE(ddiWS(0:NPROC-1,2))
      CALL rimp2_WSA(ddiWS,NPROC,1,NAUXSH,LOMP)
      LddiAuxShStart=ddiWS(ME,1)
      LddiAuxShEnd=ddiWS(ME,2)
      DEALLOCATE(ddiWS)

!!!!! ALTERNATIVELY :: NAUXBAS
      CALL RIMP2_AUXLOC(LddiAuxStart, LddiAuxShStart,.TRUE.)
      CALL RIMP2_AUXLOC(LddiAuxEnd, LddiAuxShEnd,.FALSE.)


!!!!! 3-INDEX 2-CENTER ERI
      ALLOCATE(I32(NBASIS,NBASIS,LddiAuxStart:LddiAuxEnd))


!!!!! OMP THREADS
      Nthreads_env=omp_get_max_threads()
      Nthreads=MIN0(Nthreads_env,LddiAuxShEnd-LddiAuxShStart-1)

#if __bug==1
      WRITE(*,'(A40,2I5)') "NUM_THREAD SUB RIMP2_X_FORM_MO_I32", Nthreads, Nthreads_env
#endif

!$OMP PARALLEL NUM_THREADS(NTHREADS)                                    &
!$omp default(none)                                                     &
!$omp shared(LddiAuxShStart,LddiAuxShEnd,Nthreads)                      &
!$omp shared(KAUXAT,C,EXAUX,KAUXTY,KAUXMI,KAUXMX,KAUXNG,KAUXST,KAUXLO)  &
!$omp shared(FlgSphAux,LOCSPH)                                          &
!$omp shared(NSHELL,NBASIS,EX,KATOM,KTYPE,KMIN,KMAX,KNG,KSTART,KLOC)    &
!$omp shared(MXBSH,ISIZE,MaxT12, NAUXSH,NAUXBAS)                        &
!$omp private(KAUXSH,KAT,CC,LANGK,MINK,MAXK,NK,NGK,L1,LOCK)             &
!$omp private(MINKS,MAXKS,LOCKS,NKS,MINKT,MAXKT,NKT,LOCKT)              &
!$omp private(IauxStart,IauxEnd)                                        &
!$omp private(ISH,IAT,LANGI,MINI,MAXI,NI,NGI,I1,LOCI,KANDI)             &
!$omp private(JSH,JAT,LANGJ,MINJ,MAXJ,NJ,NGJ,J1,LOCJ,IANDJ)             &
!$omp private(GG,TL,T1,XX1,XX2)                                         &
!$omp shared(VEC,I32)

      ALLOCATE(T1(MaxT12))
      ALLOCATE(GG(MaxT12))
      ALLOCATE(TL(MXBSH*MXBSH))
      ALLOCATE(XX1(ISIZE))
      ALLOCATE(XX2(ISIZE*28))

!!!!! K-SHELL (AUX)
!$omp do schedule(DYNAMIC)
      ! DO KAUXSH=1,NAUXSH
      DO KAUXSH=LddiAuxShStart,LddiAuxShEnd

            KAT = KAUXAT(KAUXSH)
            CC(1,3) = C(1,KAT)
            CC(2,3) = C(2,KAT)
            CC(3,3) = C(3,KAT)
            LANGK = KAUXTY(KAUXSH)-1
            MINK  = KAUXMI(KAUXSH)
            MAXK  = KAUXMX(KAUXSH)
            NK    = MAXK-MINK+1
            NGK   = KAUXNG(KAUXSH)
            L1    = KAUXST(KAUXSH)
            LOCK  = KAUXLO(KAUXSH)-MINK

            IF(FlgSphAux) THEN
                  CALL RIMP2CSTRM(TL,LANGK,MINK,NK,MINKS,MAXKS)
                  LOCKS = LOCSPH(KAUXSH)-MINKS
                  NKS = MAXKS-MINKS+1
                  MINKT = MINKS
                  MAXKT = MAXKS
                  LOCKT = LOCKS
            ELSE
                  MINKT = MINK
                  MAXKT = MAXK
                  LOCKT = LOCK
            ENDIF

            NKT=MAXKT-MINKT+1

            IauxStart=LOCKT+MINKT
            IauxEnd=IauxStart+NKT-1

!!!!! I-SHELL (AO)
            DO ISH = 1, NSHELL
                  IAT = KATOM(ISH)
                  CC(1,1) = C(1,IAT)
                  CC(2,1) = C(2,IAT)
                  CC(3,1) = C(3,IAT)
                  LANGI = KTYPE(ISH)-1
                  MINI  = KMIN(ISH)
                  MAXI  = KMAX(ISH)
                  NI    = KMAX(ISH)-KMIN(ISH)+1
                  NGI   = KNG(ISH)
                  I1    = KSTART(ISH)
                  LOCI  = KLOC(ISH)-MINI
                  KANDI = (KAUXSH==ISH)

!!!!! J-SHELL (AO)
                  DO JSH = 1, ISH
                        JAT = KATOM(JSH)
                        CC(1,2) = C(1,JAT)
                        CC(2,2) = C(2,JAT)
                        CC(3,2) = C(3,JAT)
                        LANGJ = KTYPE(JSH)-1
                        MINJ  = KMIN(JSH)
                        MAXJ  = KMAX(JSH)
                        NJ    = KMAX(JSH)-KMIN(JSH)+1
                        NGJ   = KNG(JSH)
                        J1    = KSTART(JSH)
                        LOCJ  = KLOC(JSH)-MINJ
                        IANDJ = (ISH==JSH)

                        CALL RIMP2_X_FORM_I32                           &
                            (GG, KAUXSH,ISH,JSH,                        &
                             CC, EX(I1),EX(J1),EXAUX(L1),               &
                             LANGI,LANGJ,LANGK,                         &
                             NI,NJ,NK,                                  &
                             NGI,NGJ,NGK,                               &
                             MINI,MAXI,MINJ,MAXJ,MINK,MAXK,             &
                             XX1,XX2)

                        CALL CP3CAUXINT(GG,I32(1,1,IauxStart),TL,T1,    &
                             NBASIS,LOCI,LOCJ,                          &
                             MINI,MAXI,MINJ,MAXJ,MINK,MAXK,             &
                             NI,NJ,NK,MINKS,MAXKS,                      &
                             NKS,FlgSphAux,IANDJ)

                  ENDDO
            ENDDO

            CALL RIMP2_MO_TRANSFORMED_I32                               &
                (I32(1,1,IauxStart),VEC,T1,                             &
                 LOCKT,MINKT,MAXKT,NBASIS,MXBSH)

      ENDDO
!$omp end do

      DEALLOCATE(T1,GG,TL,XX1,XX2)

!$OMP END PARALLEL


!!!!! DDI PUT :: I32 (NAUXBAS*NBF,NBF)
      DO LL=LddiAuxStart,LddiAuxEnd
         IBG=NBASIS*(LL-1)+1
         IED=NBASIS*LL
         CALL DDI_PUT(D_B, IBG,IED, 1,NBASIS, I32(1,1,LL))
      ENDDO


      RETURN

      END !***********************************************************





      SUBROUTINE RIMP2_AUXLOC(IauxLOC, KAUXSH,LSOE)

      use Rimp2_Shared_Data

      implicit double precision(a-h,o-z)

      COMMON /AUXBAS/ EXAUX(MXAXGTOT),CAUXS(MXAXGTOT),CAUXP(MXAXGTOT),  &
                      CAUXD(MXAXGTOT),CAUXF(MXAXGTOT),CAUXG(MXAXGTOT),  &
                      CAUXH(MXAXGTOT),CAUXI(MXAXGTOT),                  &
                      KAUXST(MXAUXSH),KAUXAT(MXAUXSH),KAUXTY(MXAUXSH),  &
                      KAUXNG(MXAUXSH),KAUXLO(MXAUXSH),KAUXMI(MXAUXSH),  &
                      KAUXMX(MXAUXSH),NAUXSH

      double precision,allocatable :: TL(:)

      logical :: LSOE


      MXBSH=MaxAUXBFX
      IF(MaxAUXANG < MaxATMANG) MXBSH=MaxATMBFX

      LANGK = KAUXTY(KAUXSH)-1
      MINK  = KAUXMI(KAUXSH)
      MAXK  = KAUXMX(KAUXSH)
      NK    = MAXK-MINK+1
      LOCK  = KAUXLO(KAUXSH)-MINK

      ALLOCATE(TL(MXBSH*MXBSH))
      IF(FlgSphAux) THEN
            CALL RIMP2CSTRM(TL,LANGK,MINK,NK,MINKS,MAXKS)
            LOCKS = LOCSPH(KAUXSH)-MINKS
            NKS = MAXKS-MINKS+1
            MINKT = MINKS
            MAXKT = MAXKS
            LOCKT = LOCKS
      ELSE
            MINKT = MINK
            MAXKT = MAXK
            LOCKT = LOCK
      ENDIF

      NKT=MAXKT-MINKT+1

      IauxStart=LOCKT+MINKT
      IauxEnd=IauxStart+NKT-1

      IF(LSOE) THEN
         IauxLOC=IauxStart
      ELSE
         IauxLOC=IauxEnd
      ENDIF

      END




      SUBROUTINE RIMP2_X_FORM_I32                                       &
                (GG, KAUXSH,ISH,JSH,                                    &
                 CC,EXI,EXJ,EXK,                                        &
                 LI,LJ,LK,                                              &
                 NI,NJ,NK,                                              &
                 NGI,NGJ,NGK,                                           &
                 MINI,MAXI,MINJ,MAXJ,MINK,MAXK,                         &
                 AKINVV,DK)

      use Rimp2_Shared_Data

      implicit double precision(a-h,o-z)

      integer,parameter :: MXLTOT=6,ITOLRI=10
      double precision,parameter :: RLN10=2.30258D00

      COMMON /ROOT  / XX,U(13),W(13),NROOTS
!$omp threadprivate(/ROOT  /)


!!!!! OUTPUT
      double precision :: GG(NI*NJ*NK)

!!!!! INPUT
      double precision :: CC(3,3), EXI(NGI),EXJ(NGJ),EXK(NGK),          &
                          AKINVV(NGK),DK(28,NGK),                       &
                          G2(2730),G3(10290),CAB(3,0:MXLTOT),           &
                          CNK(28),DB(784),P(3),PAV(3),CPV(3),           &
                          A1(3,13),B1(13),C1(13),                       &
                          A2(3,13),B2(13),F00(13)
      integer LLI(3,28),LLJ(3,28),LLK(3,28)



      TOL = ITOLRI*RLN10

      RAB2 = 0.0D00
      DO IA = 1, 3
         VAL = CC(IA,1)-CC(IA,2)
         RAB2 = RAB2+VAL*VAL
         DO J = 0, LJ
            CAB(IA,J) = VAL**J
         END DO
      END DO

      DO J = 0, LJ
         JJ = J*(J+1)/2+1
         CNK(JJ) = 1.0D00
         JJ = (J+1)*(J+2)/2
         CNK(JJ) = 1.0D00
         IF(J >= 2) THEN
            DO K = 1, J-1
               JJ = J*(J+1)/2+K+1
               JA = (J-1)*J/2+K
               JB = (J-1)*J/2+K+1
               CNK(JJ) = CNK(JA)+CNK(JB)
            END DO
         END IF
      END DO

      ICNT = 0
      DO I = MINI, MAXI
         ICNT = ICNT+1
         ! ICNT=(I-MINI)+1
         LLI(1,ICNT) = LX(I)
         LLI(2,ICNT) = LY(I)
         LLI(3,ICNT) = LZ(I)
      END DO

      JCNT = 0
      DO J = MINJ, MAXJ
         JCNT = JCNT+1
         ! JCNT=(J-MINJ)+1
         LLJ(1,JCNT) = LX(J)
         LLJ(2,JCNT) = LY(J)
         LLJ(3,JCNT) = LZ(J)
      END DO

      KCNT = 0
      DO K = MINK, MAXK
         KCNT = KCNT+1
         ! KCNT=(K-MINK)+1
         LLK(1,KCNT) = LX(K)
         LLK(2,KCNT) = LY(K)
         LLK(3,KCNT) = LZ(K)
      END DO

!!!!! DENSITY FACTOR
      DO KG = 1, NGK
         AKINVV(KG) = 1.0D00/EXK(KG)
         DO K=1,NK
            DK(K,KG)=DAux(KAUXSH,KG,K)
         ENDDO
      ENDDO

!!!!! TOTAL ANG MOM & ROOTS
      LB = LI+LJ
      LTOT = LB+LK
      NROOTS = LTOT/2+1

!!!!! ZERO OUT INTEGRAL ARRAY
      GG=0.0D00

!!!!! I-SHELL
      DO IG = 1, NGI
         AI = EXI(IG)

         TMP1 = AI*CC(1,1)
         TMP2 = AI*CC(2,1)
         TMP3 = AI*CC(3,1)

         AIRAB = AI*RAB2

!!!!! J-SHELL
         DO JG = 1, NGJ

            AJ = EXJ(JG)
            AB = AI+AJ
            ABINV = 1.0D00/AB

            DUM = AJ*ABINV*AIRAB

            IF (DUM > TOL) CYCLE

            EXPB = EXP(-DUM)

            P(1) = (TMP1+AJ*CC(1,2))*ABINV
            P(2) = (TMP2+AJ*CC(2,2))*ABINV
            P(3) = (TMP3+AJ*CC(3,2))*ABINV

            PAV(1) = P(1)-CC(1,1)
            PAV(2) = P(2)-CC(2,1)
            PAV(3) = P(3)-CC(3,1)

            CPV(1) = CC(1,3)-P(1)
            CPV(2) = CC(2,3)-P(2)
            CPV(3) = CC(3,3)-P(3)
            RR = CPV(1)*CPV(1)+CPV(2)*CPV(2)+CPV(3)*CPV(3)

            DUM0 = ABINV*EXPB

            IJ = 0
            DO I=1,NI
               DUM2 = DUM0*DAtm(ISH,IG,I)
               DO J=1,NJ
                  IJ = IJ+1
                  DB(IJ) = DUM2*DAtm(JSH,JG,J)
               ENDDO
            ENDDO !I

!!!!! K-SHELL
            DO KG = 1, NGK

               AK = EXK(KG)

               RHO = AB*AK/(AB+AK)
               XX = RHO*RR

!!!!! RYS ROOTS AND WEIGHT
               SELECT CASE(NROOTS)
                  CASE(:3)
                     CALL RT123
                  CASE(4)
                     CALL ROOT4
                  CASE(5)
                     CALL ROOT5
                  CASE(6:)
                     CALL ROOT6
               END SELECT

!!!!! (G2(XYZ,IR,LI+LJ,LK))
               CALL I32XYZ_1                                            &
                   (G2, U,W,PAV,CPV,                                    &
                    A1,B1,C1,                                           &
                    A2,B2,F00,                                          &
                    AB,AK,LB,LK,NROOTS)

!!!!! IXYZ G3(XYZ,IR,LI,LJ,LK)
               CALL I32XYZ_2(G3, G2, CAB,LI,LJ,LK,LB,NROOTS)

!!!!! FORM I32
               CALL FORM_I32                                            &
                   (GG, G3, DB,DK(1,KG),                                &
                    LLI,LLJ,LLK,                                        &
                    LI,LJ,LK,                                           &
                    NI,NJ,NK,NROOTS)

            END DO !KG
         END DO !JG
      END DO !IG

      RETURN

      END !***********************************************************





      SUBROUTINE I32XYZ_1(G, U,W,PA,CP,A1,B1,C1,                        &
                 A2,B2,F00,AI,AJ,LI,LJ,NROOTS)

      use Rimp2_Shared_Data, only : HALF, one

      implicit double precision(a-h,o-z)

      double precision G(3,NROOTS,0:LI,0:LJ)

      double precision U(NROOTS),W(NROOTS),PA(3),CP(3)
      double precision A1(3,NROOTS),B1(NROOTS),C1(NROOTS),              &
                       A2(3,NROOTS),B2(NROOTS),F00(NROOTS)



!!!!! HEAD UP: LI ~LB=LI+LJ IN THE SUB CALL I32XYZ_1

      AIJINV = 1.0D00/(AI+AJ)

      RHO = AI*AJ*AIJINV

      AIINV  = 1.0D00/AI
      AJINV  = 1.0D00/AJ

      AIINV2 = 0.5D00*AIINV
      AJINV2 = 0.5D00*AJINV
      AIJINV2 = 0.5D00*AIJINV
      AIJINVSQ = SQRT(AIJINV)

      DO IR = 1, NROOTS

         U2 = U(IR)*RHO
         T2 = U2/(RHO+U2)

         VAL1 = RHO*AIINV*T2
         VAL2 = RHO*AJINV*T2
         VAL3 = AIINV2-AIINV2*VAL1
         VAL4 = AJINV2-AJINV2*VAL2

         A1(1,IR) = PA(1)+VAL1*CP(1)
         A1(2,IR) = PA(2)+VAL1*CP(2)
         A1(3,IR) = PA(3)+VAL1*CP(3)

         A2(1,IR) = -VAL2*CP(1)
         A2(2,IR) = -VAL2*CP(2)
         A2(3,IR) = -VAL2*CP(3)

         B1(IR) = VAL3
         B2(IR) = VAL4

         C1(IR) = AIJINV2*T2

         F00(IR) = AIJINVSQ*W(IR)

      END DO

!!!!! --- (s|s) -------------------------------------------------------

      DO IR = 1, NROOTS
         G(1,IR,0,0) = 1.0D00
         G(2,IR,0,0) = 1.0D00
         G(3,IR,0,0) = F00(IR)
      END DO

!!!!! --- (p|s); (d|s), (f|s), (g|s) ----------------------------------

      SELECT CASE(LI)
         CASE(1)
            CALL L10(G, LI,LJ, NROOTS,A1,F00)
         CASE(2:)
            CALL L20(G, LI,LJ, NROOTS,A1,B1,F00)
      END SELECT

      IF (LJ==0) RETURN

!!!!! --- (s|p); (p|p); (d|p),(f|p),(g|p) -----------------------------

      SELECT CASE(LI)
         CASE(0)
            CALL L01(G, LI,LJ, NROOTS,A2,F00)
         CASE(1)
            CALL L11(G, LI,LJ, NROOTS,A1,A2,C1,F00)
         CASE(2:)
            CALL L21(G, LI,LJ, NROOTS,A1,A2,B1,C1,F00)
      END SELECT

      IF(LJ==1) RETURN


!!!!! --- (s|d),(s|f),(s|g); ---
!!!!! --- (p|d),(p|f),(p|g); ---
!!!!! --- (d|d),(d|f),(d|g); (f|d),(f|f),(f|g); (g|d),(g|f),(g|g) -----

      SELECT CASE(LI)
         CASE(0)
            CALL L02(G, LI,LJ, NROOTS,A2,B2)
         CASE(1)
            CALL L12(G, LI,LJ, NROOTS,A1,A2,B2,C1)
         CASE(2:)
            CALL L22(G, LI,LJ, NROOTS,A1,A2,B1,B2,C1)
      END SELECT

      RETURN

      END !***********************************************************



!
      SUBROUTINE L10(G, LI,LJ, NROOTS,A1,F00)

      implicit double precision(a-h,o-z)

      double precision G(3,NROOTS,0:LI,0:LJ)
      double precision A1(3,NROOTS)
      double precision F00(NROOTS)

!     --- IXYZ: (p|s) ---

      DO IR = 1, NROOTS
         G(1,IR,1,0) = A1(1,IR)
         G(2,IR,1,0) = A1(2,IR)
         G(3,IR,1,0) = A1(3,IR)*F00(IR)
      END DO

      END !***********************************************************


!
      SUBROUTINE L20(G, LI,LJ, NROOTS,A1,B1,F00)

      implicit double precision(a-h,o-z)

      double precision G(3,NROOTS,0:LI,0:LJ)
      double precision A1(3,NROOTS),B1(NROOTS)
      double precision F00(NROOTS)

!     --- IXYZ: (p|s) ---

      DO IR = 1, NROOTS
         G(1,IR,1,0) = A1(1,IR)
         G(2,IR,1,0) = A1(2,IR)
         G(3,IR,1,0) = A1(3,IR)*F00(IR)
      END DO

!     --- IXYZ: (d|s), (f|s), (g|s) ---

      DO I = 2, LI
         FORALL(IR=1:NROOTS, K=1:3)
            G(K,IR,I,0) = A1(K,IR)*G(K,IR,I-1,0)                        &
                     +(I-1)*B1(IR)*G(K,IR,I-2,0)
         ENDFORALL
      END DO

      END !***********************************************************


!
      SUBROUTINE L01(G, LI,LJ, NROOTS,A2,F00)

      implicit double precision(a-h,o-z)

      double precision G(3,NROOTS,0:LI,0:LJ)
      double precision A2(3,NROOTS)
      double precision F00(NROOTS)

!     --- IXYZ: (s|p) ---

      DO IR = 1, NROOTS
         G(1,IR,0,1) = A2(1,IR)
         G(2,IR,0,1) = A2(2,IR)
         G(3,IR,0,1) = A2(3,IR)*F00(IR)
      END DO

   END !***********************************************************


!
      SUBROUTINE L11(G, LI,LJ, NROOTS,A1,A2,C1,F00)

      implicit double precision(a-h,o-z)

      double precision G(3,NROOTS,0:LI,0:LJ)
      double precision A1(3,NROOTS),A2(3,NROOTS)
      double precision C1(NROOTS)
      double precision F00(NROOTS)

!     --- IXYZ: (s|p) ---

      DO IR = 1, NROOTS
         G(1,IR,0,1) = A2(1,IR)
         G(2,IR,0,1) = A2(2,IR)
         G(3,IR,0,1) = A2(3,IR)*F00(IR)
      END DO

!     --- IXYZ: (p|p) ---

      FORALL(K=1:3, IR=1:NROOTS)
         G(K,IR,1,1) = A1(K,IR)*G(K,IR,0,1)                             &
                        +C1(IR)*G(K,IR,0,0)
      ENDFORALL

      END !***********************************************************

!
      SUBROUTINE L21(G, LI,LJ, NROOTS,A1,A2,B1,C1,F00)

      implicit double precision(a-h,o-z)

      double precision G(3,NROOTS,0:LI,0:LJ)
      double precision A1(3,NROOTS),A2(3,NROOTS)
      double precision B1(NROOTS),C1(NROOTS)
      double precision F00(NROOTS)

!     --- IXYZ: (s|p) ---

      DO IR = 1, NROOTS
         G(1,IR,0,1) = A2(1,IR)
         G(2,IR,0,1) = A2(2,IR)
         G(3,IR,0,1) = A2(3,IR)*F00(IR)
      END DO

!     --- IXYZ: (p|p) ---

      FORALL(K=1:3, IR=1:NROOTS)
         G(K,IR,1,1) = A1(K,IR)*G(K,IR,0,1)                             &
                        +C1(IR)*G(K,IR,0,0)
      ENDFORALL

!     --- IXYZ: (d|p),(f|p),(g|p) ---

      DO I = 2, LI
         FORALL(K=1:3, IR=1:NROOTS)
            G(K,IR,I,1) = A1(K,IR)*G(K,IR,I-1,1)                        &
                     +(I-1)*B1(IR)*G(K,IR,I-2,1)                        &
                           +C1(IR)*G(K,IR,I-1,0)
         ENDFORALL
      END DO

      END !***********************************************************



      SUBROUTINE L02(G, LI,LJ, NROOTS,A2,B2)

      implicit double precision(a-h,o-z)

      double precision G(3,NROOTS,0:LI,0:LJ)
      double precision A2(3,NROOTS)
      double precision B2(NROOTS)

!     --- IXYZ: (s|d),(s|f),(s|g) ---

      DO J=2,LJ
         FORALL(K=1:3, IR=1:NROOTS)
            G(K,IR,0,J) = A2(K,IR)*G(K,IR,0,J-1)                        &
                     +(J-1)*B2(IR)*G(K,IR,0,J-2)
         ENDFORALL
      END DO

      END !***********************************************************


!
      SUBROUTINE L12(G, LI,LJ, NROOTS,A1,A2,B2,C1)

      implicit double precision(a-h,o-z)

      double precision G(3,NROOTS,0:LI,0:LJ)
      double precision A1(3,NROOTS),A2(3,NROOTS)
      double precision B2(NROOTS),C1(NROOTS)


!     --- IXYZ: (s|d),(s|f),(s|g) ---

      DO J=2,LJ
         FORALL(K=1:3, IR=1:NROOTS)
            G(K,IR,0,J) = A2(K,IR)*G(K,IR,0,J-1)                        &
                     +(J-1)*B2(IR)*G(K,IR,0,J-2)
         ENDFORALL
      END DO

!     --- IXYZ: (p|d),(p|f),(p|g) ---

      DO J = 2, LJ
         FORALL(K=1:3, IR=1:NROOTS)
            G(K,IR,1,J) = A1(K,IR)*G(K,IR,0,J)                          &
                         +J*C1(IR)*G(K,IR,0,J-1)
         ENDFORALL
      END DO

      END !***********************************************************


!
      SUBROUTINE L22(G, LI,LJ, NROOTS,A1,A2,B1,B2,C1)

      implicit double precision(a-h,o-z)

      double precision G(3,NROOTS,0:LI,0:LJ)
      double precision A1(3,NROOTS),A2(3,NROOTS)
      double precision B1(NROOTS),B2(NROOTS),C1(NROOTS)

!     --- IXYZ: (s|d),(s|f),(s|g) ---

      DO J=2,LJ
         FORALL(K=1:3, IR=1:NROOTS)
            G(K,IR,0,J) = A2(K,IR)*G(K,IR,0,J-1)                        &
                     +(J-1)*B2(IR)*G(K,IR,0,J-2)
         ENDFORALL
      END DO

!     --- IXYZ: (p|d),(p|f),(p|g) ---

      DO J = 2, LJ
         FORALL(K=1:3, IR=1:NROOTS)
            G(K,IR,1,J) = A1(K,IR)*G(K,IR,0,J)                          &
                         +J*C1(IR)*G(K,IR,0,J-1)
         ENDFORALL
      END DO

!     --- IXYZ: (d|d),(d|f),(d|g);
!               (f|d),(f|f),(f|g);
!               (g|d),(g|f),(g|g) ---

      DO J = 2, LJ
         DO I = 2, LI
            FORALL(K=1:3, IR=1:NROOTS)
               G(K,IR,I,J) = A1(K,IR)*G(K,IR,I-1,J)                     &
                       + (I-1)*B1(IR)*G(K,IR,I-2,J)                     &
                           + J*C1(IR)*G(K,IR,I-1,J-1)
            ENDFORALL
         END DO
      END DO

      END !***********************************************************





      SUBROUTINE I32XYZ_2(G3, G, CAB,LI,LJ,LK,LB,NROOTS)

      implicit double precision(a-h,o-z)

!!!!! OUTPUT
      double precision G3(3,NROOTS,0:LI,0:LJ,0:LK)

!!!!! INPUT
      double precision G(3,NROOTS,0:LB,0:LK)
      double precision CAB(3,0:LJ)



      DO KK = 0, LK
         DO JJ = 0, LJ-1

            FORALL(IA=1:3, IR=1:NROOTS, II=0:LI)
               G3(IA,IR,II,JJ,KK) = G(IA,IR,II,KK)
            ENDFORALL

            IF(LJ > 0) THEN
               DO II = 0, LI+LJ-JJ-1
                  FORALL(IA=1:3, IR=1:NROOTS)
                     G(IA,IR,II,KK) = CAB(IA,1)*G(IA,IR,II,KK)          &
                                               +G(IA,IR,II+1,KK)
                  ENDFORALL
               END DO
            END IF

         END DO !JJ
      END DO !KK


      FORALL(IA=1:3, IR=1:NROOTS, II=0:LI, KK=0:LK)
         G3(IA,IR,II,LJ,KK) = G(IA,IR,II,KK)
      ENDFORALL

      RETURN

      END !***********************************************************





      SUBROUTINE FORM_I32                                               &
                (GG, G3, DB,DK,                                         &
                 LLI,LLJ,LLK,                                           &
                 LI,LJ,LK,                                              &
                 NI,NJ,NK,NROOTS)

      implicit double precision(A-H,O-Z)

      double precision,parameter :: PI252=34.986836655250D00


!!!!! OUTPUT
      double precision :: GG(NI*NJ*NK)


!!!!! INPUT
      double precision :: G3(3,NROOTS,0:LI,0:LJ,0:LK)
      double precision :: DB(NI*NJ),DK(NK)

      integer LLI(3,NI),LLJ(3,NJ),LLK(3,NK)

!!!!! LOCAL DATA

      logical :: lomp

      IJ = 0
      IJK = 0
      DO I = 1, NI
         NIX = LLI(1,I)
         NIY = LLI(2,I)
         NIZ = LLI(3,I)
         DO J = 1, NJ
            NJX = LLJ(1,J)
            NJY = LLJ(2,J)
            NJZ = LLJ(3,J)
            IJ = IJ+1

            FACTOR=PI252*DB(IJ)

            DO K = 1, NK
               NKX = LLK(1,K)
               NKY = LLK(2,K)
               NKZ = LLK(3,K)
               VAL = 0.0D00
               DO IR = 1, NROOTS
                  VAL = VAL+G3(1,IR,NIX,NJX,NKX)                        &
                           *G3(2,IR,NIY,NJY,NKY)                        &
                           *G3(3,IR,NIZ,NJZ,NKZ)
               END DO
               IJK = IJK+1
               ! GG(IJK) = GG(IJK) + PI252*DB(IJ)*DK(K)*VAL
               GG(IJK) = GG(IJK) + FACTOR*DK(K)*VAL
            END DO
         END DO
      END DO

      RETURN

      END !***********************************************************



!>*module rimp2grd   *deck RIMP2_COMBINE_VDX_I32
!>
!>     @brief   combine 3-2ERI with VDX
!>     @author  buu
!>     @date    sept 8, 2018
!>
!>     @detail  combine 3-2ERI with VDX
!>

      SUBROUTINE RIMP2_COMBINE_VDX_I32(VDX, NAUXBAS,NAUXBASD,NBASIS)

      use Rimp2_Shared_Data
      use omp_lib

      implicit double PRECISION(A-H,O-Z)

      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK

!!!!! INPUT
      double precision :: VDX(NAUXBASD,NAUXBAS)

!!!!! LOCAL
      double precision,allocatable :: T2(:,:),T3(:,:,:)
      integer,allocatable :: ompISE(:,:)
      logical :: GOPARR,DSKWRK,MASWRK

      logical :: lomp


!!!!! DDI WORK DISTRIBUTING FOR NBASIS LOOP
      ALLOCATE(ompISE(0:NPROC-1,2))
      CALL rimp2_WSA(ompISE,NPROC,1,NBASIS,LOMP)
      LddiAoStart=ompISE(ME,1)
      LddiAoEnd=ompISE(ME,2)
      DEALLOCATE(ompISE)


!!!!! GET T3
      ALLOCATE(T3(NBASIS,NAUXBAS,LddiAoStart:LddiAoEnd))
      CALL DDI_GET(D_B, 1,NAUXBAS*NBASIS, LddiAoStart,LddiAoEnd, T3)


!!!!! OMP THREADS
      Nthreads_env=omp_get_max_threads()
      Nthreads=MIN0(Nthreads_env,LddiAoEnd-LddiAoStart-1)
      IF(BUG) WRITE(*,'(A40,2I5)') "NUM_THREAD SUB RIMP2_COMBINE_VDX_I32", Nthreads, Nthreads_env

!$OMP PARALLEL NUM_THREADS(NTHREADS)                                    &
!$omp default(none)                                                     &
!$omp shared(LddiAoStart,LddiAoEnd)                                     &
!$omp private(IthreadID,LompStart,LompEnd)                              &
!$omp shared(T3,NBASIS,NAUXBAS)                                         &
!$omp private(IBF,T2)                                                   &
!$omp shared(VDX,NAUXBASD)


      ALLOCATE(T2(NBASIS,NAUXBAS))
!$omp do schedule(guided)
      ! DO IBF=1,NBASIS
      DO IBF=LddiAoStart,LddiAoEnd

         CALL DCOPY(NBASIS*NAUXBAS, T3(1,1,IBF),1, T2,1)
         CALL DGEMM('N','T', NAUXBASD,NBASIS,NAUXBAS,                   &
              1.0D00, VDX,NAUXBASD,                                     &
                      T2,NBASIS,                                        &
              0.0D00, T3(1,1,IBF),NAUXBASD)


      ENDDO
!$omp end do
      DEALLOCATE(T2)

!$OMP END PARALLEL


!!!!! PUT V & MO-TRANSFORMED I32 BACK TO DDI ARRAYS
      CALL DDI_PUT(D_B, 1, NAUXBASD*NBASIS, LddiAoStart,LddiAoEnd, T3)


      RETURN

      END !***********************************************************






      SUBROUTINE RIMP2_MO_TRANSFORMED_I32                               &
                 (T0,VEC,T1,LOCLT,MINLT,MAXLT,NBASIS,MXBSH)

      use Rimp2_Shared_Data
      use omp_lib

      implicit double precision(a-h,o-z)

      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      logical :: GOPARR,DSKWRK,MASWRK

!!!!! INPUT
      double precision :: VEC(NBASIS,NBASIS)
      double precision :: T0(NBASIS,NBASIS,MXBSH)
      double precision :: T1(NBASIS*NBASIS)



      DO L=MINLT,MAXLT
            LL=L-MINLT+1

            CALL DGEMM('N','N', NBASIS,NBASIS, NBASIS,                  &
                 1.0D00, T0(1,1,LL),NBASIS,                             &
                          VEC,NBASIS,                                   &
                 0.0D00, T1,NBASIS)

            CALL DGEMM('T','N',NBASIS,NBASIS,NBASIS,                    &
                 1.0D00, VEC,NBASIS,                                    &
                          T1,NBASIS,                                    &
                 0.0D00, T0(1,1,LL),NBASIS)
      ENDDO


      END !***********************************************************



!>*module rimp2grd   *deck RIMP2_D22GRD
!>
!>     @brief   combine 2-2ERI gradient with 2-2PDM
!>     @author  buu
!>     @date    sept 8, 2018
!>
!>     @detail  combine 2-2ERI gradient with 2-2PDM
!>
      SUBROUTINE RIMP2_D22GRD                                           &
                (DE_DDI, gDX,VDX, NAUXBASD,NAUXBAS,                     &
                 LddiAuxShStart,LddiAuxShEnd,                           &
                 LddiAuxStart,LddiAuxEnd)

      use omp_lib
      use Rimp2_Shared_Data

      implicit double PRECISION(A-H,O-Z)

      COMMON /AUXBAS/ EXAUX(MXAXGTOT),CAUXS(MXAXGTOT),CAUXP(MXAXGTOT),  &
                      CAUXD(MXAXGTOT),CAUXF(MXAXGTOT),CAUXG(MXAXGTOT),  &
                      CAUXH(MXAXGTOT),CAUXI(MXAXGTOT),                  &
                      KAUXST(MXAUXSH),KAUXAT(MXAUXSH),KAUXTY(MXAUXSH),  &
                      KAUXNG(MXAUXSH),KAUXLO(MXAUXSH),KAUXMI(MXAUXSH),  &
                      KAUXMX(MXAUXSH),NAUXSH
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,                    &
                      ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK


!!!!! OUTPUT
      double precision :: DE_DDI(3,NAT)


!!!!! INPUT
      double precision,allocatable :: gXX(:,:)
      double precision :: gDX(NAUXBASD,NAUXBAS)
      double precision :: VDX(NAUXBASD,NAUXBAS)

!!!!! LOCAL DATA
      double precision,save,allocatable ::                              &
                D22(:,:),DE_OMP(:,:),TI(:),TJ(:),W(:)
!$omp threadprivate(D22,DE_OMP,TI,TJ,W)

      double precision :: CC(3,2)
      logical :: GOPARR,DSKWRK,MASWRK,LOMP,KANDL
      integer,allocatable :: ompISE(:,:)


!!!!! TURN OFF DYNAMIC THREADS
!!!!! TO KEEP threadprivate DATA CONSISTENT
      CALL OMP_SET_DYNAMIC(.FALSE.)


!!!!! OMP THREADS
      Nthreads_env=omp_get_max_threads()
      Nthreads=MIN0(Nthreads_env,NAUXBAS)
      IF(BUG) WRITE(*,'(A40,2I5)') "NUM_THREAD SUB RIMP2_D22GRD P1", Nthreads, Nthreads_env


!!!!! OMP WORK DISTRIBUTION
      ALLOCATE(ompISE(0:Nthreads-1,2))
      CALL WorkSharedArray(ompISE,Nthreads,1,NAUXBAS,LOMP)

      ALLOCATE(gXX(NAUXBAS,NAUXBAS))
      gXX=0.0D00

!$OMP PARALLEL NUM_THREADS(Nthreads)                                    &
!$omp private(IthreadID,Istart,Iend,Ichunk)                             &
!$omp shared(ompISE, NAUXBAS,NAUXBASD,VDX,gDX,gXX)

      IthreadID = omp_get_thread_num()
      Istart = ompISE(IthreadID,1)
      Iend = ompISE(IthreadID,2)
      Ichunk = Iend-Istart+1

! write(*,*) "Istart,Iend", Istart, Iend, Iend-Istart+1

      CALL DGEMM('T','N', NAUXBAS,Ichunk,NAUXBASD,                      &
                  4.0D00, VDX,NAUXBASD,                                 &
                          gDX(1,Istart),NAUXBASD,                       &
                  0.0D00, gXX(1,Istart),NAUXBAS)
!$OMP END PARALLEL

      call DDI_GSUMF(5463,gXX,NAUXBAS*NAUXBAS)


!!!!! OMP THREADS
      Nthreads_env=omp_get_max_threads()
      Nthreads=Nthreads_env
      IF(BUG) WRITE(*,'(A40,2I5)') "NUM_THREAD SUB RIMP2_D22GRD P2", Nthreads, Nthreads_env


!!!!! MAX BASIS FUNCTIONS IN A PAIR OF IJ AUX SHELL
      MaxG=MaxAuxBfx*MaxAuxBfx


!$OMP PARALLEL NUM_THREADS(Nthreads)                                    &
!$omp default(none)                                                     &
!$omp shared(DE_DDI,LddiAuxShStart,LddiAuxShEnd)                        &
!$omp shared(LddiAuxStart,LddiAuxEnd)                                   &
!$omp shared(EXAUX,KAUXAT,C,KAUXTY,KAUXNG,KAUXMI,KAUXMX,KAUXST,KAUXLO)  &
!$omp shared(NAUXSH,NAUXBAS,NAT)                                        &
!$omp shared(gXX,FlgSphAux,LOCSPH,MaxG)                                 &
!$omp private(LAUX,LAT,CC,LANGI,LANGIDER,NGI,MINI,MAXI,I1,LOCI,NI,KANDL)&
!$omp private(KAUX,KAT,LANGJ,NGJ,MINJ,MAXJ,J1,LOCJ,NJ)                  &
!$omp private(MINIS,MAXIS,LOCIS,NIS)                                    &
!$omp private(MINJS,MAXJS,LOCJS,NJS)


!!!!! ALLOCATE ARRAYS
      ALLOCATE(TI(MaxG))
      ALLOCATE(TJ(MaxG))
      ALLOCATE(W(MaxG))
      ALLOCATE(D22(MaxG,3))
      ALLOCATE(DE_OMP(3,NAT))

      DE_OMP=0.0D00

!$omp do schedule(DYNAMIC)

      ! DO LAUX = 1, NAUXSH
      DO LAUX = LddiAuxShStart,LddiAuxShEnd

         LAT = KAUXAT(LAUX)
         CC(1,1) = C(1,LAT)
         CC(2,1) = C(2,LAT)
         CC(3,1) = C(3,LAT)
         LANGI = KAUXTY(LAUX)-1
         NGI   = KAUXNG(LAUX)
         MINI  = KAUXMI(LAUX)
         MAXI  = KAUXMX(LAUX)
         I1    = KAUXST(LAUX)
         LOCI  = KAUXLO(LAUX)-MINI
         NI    = MAXI-MINI+1

         IF(FlgSphAux) THEN
            CALL RIMP2CSTRM(TI,LANGI,MINI,NI,MINIS,MAXIS)
            LOCIS = LOCSPH(LAUX)-MINIS
            NIS   = MAXIS-MINIS+1
         END IF

         ! DO KAUX = 1, NAUXSH
         DO KAUX = 1, LAUX

            KAT = KAUXAT(KAUX)

            IF(LAT.EQ.KAT) CYCLE

            CC(1,2) = C(1,KAT)
            CC(2,2) = C(2,KAT)
            CC(3,2) = C(3,KAT)
            LANGJ = KAUXTY(KAUX)-1
            NGJ   = KAUXNG(KAUX)
            MINJ  = KAUXMI(KAUX)
            MAXJ  = KAUXMX(KAUX)
            J1    = KAUXST(KAUX)
            LOCJ  = KAUXLO(KAUX)-MINJ
            NJ = MAXJ-MINJ+1
            KANDL=KAUX.EQ.LAUX

!!!!! CARTESIAN TO SPHERICAL AUX SHELL
            IF(FlgSphAux) THEN
               CALL RIMP2CSTRM(TJ,LANGJ,MINJ,NJ,MINJS,MAXJS)
               LOCJS = LOCSPH(KAUX)-MINJS
               NJS   = MAXJS-MINJS+1
            END IF

!!!!! INCREASE ANG MOMENTUM
            LANGIDER=LANGI+1

            CALL RIMP2_X_FORM_D22                                       &
                (D22, LAUX,KAUX,CC,                                     &
                 EXAUX(I1),EXAUX(J1),                                   &
                 LANGIDER,LANGJ,                                        &
                 NI,NJ,                                                 &
                 NGI,NGJ,                                               &
                 MINI,MAXI, MINJ,MAXJ)

            CALL RIMP2_D22_COMBINED_DENS                                &
                (DE_OMP,D22, gXX,                                       &
                 TI,TJ,W,NAUXBAS,                                       &
                 LAT,KAT,                                               &
                 LOCI,LOCJ, MINI,MAXI,MINJ,MAXJ, NI,NJ,                 &
                 LOCIS,LOCJS, MINIS,MAXIS,MINJS,MAXJS, NIS,NJS,         &
                 FlgSphAux,KANDL)

         END DO
      END DO

!$omp end do

!$omp critical(DE_DDI_critical)
      DE_DDI=DE_DDI+DE_OMP
!$omp end critical(DE_DDI_critical)

      DEALLOCATE(TI,TJ,W,D22,DE_OMP)

!$OMP END PARALLEL

      RETURN

      END !***********************************************************














!
      SUBROUTINE RIMP2_X_FORM_D22                                  &
                (D22, LAUX,KAUX,CC,                                     &
                 EXI,EXJ, LI,LJ, NI,NJ, NGI,NGJ,                        &
                 MINI,MAXI, MINJ,MAXJ)

      use omp_lib
      use Rimp2_Shared_Data

      implicit double precision(a-h,o-z)

      COMMON /ROOT  / XX,U(13),W(13),NROOTS
!$omp threadprivate(/ROOT  /)


!!!!! OUTPUT
      double precision :: D22(NI*NJ,3)


!!!!! INPUT
      double precision :: CC(3,2)
      double precision :: EXI(NGI),EXJ(NGJ)
      double precision :: DI(28),DJ(28),G(1911)
      double precision :: A1(3,13),B1(13),C1(13)
      double precision :: A2(3,13),B2(13),F00(13)


!!!!! LOCAL DATA
      integer :: LLI(3,28),LLJ(3,28)




!!!!! TOT ANG MOMENTUM & #ROOTS
      LTOT = LI+LJ
      NROOTS = LTOT/2+1


!!!!! KL SHELL CENTER DISTANCE
      RXIJ = CC(1,1)-CC(1,2)
      RYIJ = CC(2,1)-CC(2,2)
      RZIJ = CC(3,1)-CC(3,2)
      RIJ2 = RXIJ**2+RYIJ**2+RZIJ**2

!!!!! K&L SHELL ANG MOMENTUM
      DO I = MINI, MAXI
         ICNT=(I-MINI)+1
         LLI(1,ICNT) = LX(I)
         LLI(2,ICNT) = LY(I)
         LLI(3,ICNT) = LZ(I)
      END DO

      DO J = MINJ, MAXJ
         JCNT=(J-MINJ)+1
         LLJ(1,JCNT) = LX(J)
         LLJ(2,JCNT) = LY(J)
         LLJ(3,JCNT) = LZ(J)
      END DO

!!!!! ZERO OUT D22
      D22=0.0D00

!!!!! L-SHELL GAUSSIANS
      DO LG = 1, NGI
         AI  = EXI(LG)
         AIINV = 1.0D00/AI

!!!!! L-SHELL CONTRACTION COEFFICIENTS
         DO I=1,NI
            DI(I)=DAux(LAUX,LG,I)
         ENDDO

!!!!! K-SHELL GAUSSIANS
         DO KG = 1, NGJ
            AJ  = EXJ(KG)
            AJINV = 1.0D00/AJ

!!!!! J-SHELL CONTRACTION COEFFICIENTS
            DO J=1,NJ
              DJ(J)=DAux(KAUX,KG,J)
            ENDDO

            RHO = AI*AJ/(AI+AJ)
            XX = RHO*RIJ2

!!!!! RYS ROOTS AND WEIGHTS
            SELECT CASE(NROOTS)
               CASE(:3)
                  CALL RT123
               CASE(4)
                  CALL ROOT4
               CASE(5)
                  CALL ROOT5
               CASE(6:)
                  CALL ROOT6
            END SELECT

!!!!! IXYZ G(XYZ,IR,LI,LJ)
            CALL I22XYZ                                                 &
                (G, CC,U,W,                                             &
                 A1,B1,C1, A2,B2, F00,                                  &
                 AI,AJ, LI,LJ, NROOTS)

!!!!! FORM D22
            CALL RIMP2_FORM_D22                                         &
                (D22, G, AI,                                            &
                 DI,DJ,                                                 &
                 LLI,LLJ, LI,LJ, NI,NJ, NROOTS)
         END DO
      END DO

      RETURN

      END !***********************************************************





      SUBROUTINE RIMP2_FORM_D22                                         &
                (GG, G, AI,                                             &
                 DI,DJ,                                                 &
                 LLI,LLJ, LI,LJ, NI,NJ, NROOTS)
!
      use Rimp2_Shared_Data
!
      implicit double PRECISION(A-H,O-Z)

!    OUTPUT
      double precision :: GG(NI*NJ,3)

!    INPUT
      double precision :: G(3,NROOTS,0:LI,0:LJ)
      double precision :: DI(NI),DJ(NJ)
      integer :: LLI(3,NI),LLJ(3,NJ)

!    LOCAL DATA
      double precision :: VAL(3),VALXYZ(12,3),XYZ(12,3)


!
      AI2=2.0D00*AI

      IJ = 0
      DO J = 1, NJ
         NJX = LLJ(1,J)
         NJY = LLJ(2,J)
         NJZ = LLJ(3,J)
!
         FACXXX=PI252*DJ(J)
!
         DO I = 1, NI
            NIX = LLI(1,I)
            NIY = LLI(2,I)
            NIZ = LLI(3,I)
!
            FACTOR = FACXXX*DI(I)
!
            DO IR=1,NROOTS
               XYZ(IR,1)=G(2,IR,NIY,NJY)*G(3,IR,NIZ,NJZ)
               XYZ(IR,2)=G(1,IR,NIX,NJX)*G(3,IR,NIZ,NJZ)
               XYZ(IR,3)=G(1,IR,NIX,NJX)*G(2,IR,NIY,NJY)
            ENDDO

!
            VAL=0.0D00
            IJ = IJ+1

!          DX,DY,DZ
            DO II=1,3
               IF(LLI(II,I) == 0) THEN
                  DO IR=1,NROOTS
                     VALXYZ(IR,II) = AI2*G(II,IR,LLI(II,I)+1,LLJ(II,J))
                  ENDDO
               ELSE
                  DO IR=1,NROOTS
                     VALXYZ(IR,II) = AI2*G(II,IR,LLI(II,I)+1,LLJ(II,J)) &
                              -LLI(II,I)*G(II,IR,LLI(II,I)-1,LLJ(II,J))
                  ENDDO
               ENDIF
!
               DO IR=1,NROOTS
                  VAL(II) = VAL(II) + VALXYZ(IR,II)*XYZ(IR,II)
               ENDDO
               GG(IJ,II) = GG(IJ,II) + FACTOR*VAL(II)
            ENDDO !II

         END DO
      END DO

      RETURN

      END !***********************************************************


!
      SUBROUTINE RIMP2OMP_D22_COMBINED_DENS                             &
                (D22, gXX, TI,TJ,W,NAUXBAS,                             &
                 LAT,KAT,                                               &
                 LOCI,LOCJ, MINI,MAXI, MINJ,MAXJ, NI,NJ,                &
                 LOCIS,LOCJS, MINIS,MAXIS, MINJS,MAXJS, NIS,NJS,        &
                 SPHAUX)

!    OPENMP LIB
      use omp_lib

!    GLOBAL DATA
      use Rimp2_Shared_Data, only: MXATM
!
      implicit double precision(a-h,o-z)
!
      COMMON /GRAD  / DE(3,MXATM)
!$omp threadprivate(/GRAD  /)

!    INPUT
      double precision :: D22(NI*NJ,3)
      double precision :: gXX(NAUXBAS,NAUXBAS)
      double precision :: TI(NI*NI)
      double precision :: TJ(NJ*NJ)
      double precision :: W(NI*NJ)
      logical :: SPHAUX !,IANDJ

!    LOCAL DATA
      double precision :: VAL(3)
      double precision :: DEL(3),DEK(3)


                    !!!  \\\ (^^) /// !!!


!    CARTESIAN TO SPHERICAL AUX SHELL
      IF(SPHAUX) THEN
         MINIT = MINIS
         MAXIT = MAXIS
         MINJT = MINJS
         MAXJT = MAXJS
         LOCIT = LOCIS
         LOCJT = LOCJS
         NIT = NIS
         NJT = NJS
!
         DO II=1,3
            CALL DGEMM('N','N', NI,NJS, NJ,                             &
                 1.0D00, D22(1,II),NI,                                 &
                          TJ,NJ,                                        &
                 0.0D00, W,NI)
            CALL DGEMM('T','N', NIS,NJS, NI,                            &
                 1.0D00, TI,NI,                                        &
                          W,NI,                                         &
                 0.0D00, D22(1,II),NIS)
         ENDDO
!
      ELSE
         MINIT = MINI
         MAXIT = MAXI
         MINJT = MINJ
         MAXJT = MAXJ
         LOCIT = LOCI
         LOCJT = LOCJ
         NIT = NI
         NJT = NJ
      END IF


!    INITIALIZATION
      DEL=0.0D00
      DEK=0.0D00


!    D22 * DENSITY
      DO I = MINIT, MAXIT
         ICNT=(I-MINIT)+1
         IAUX = LOCIT+I

         DO J=MINJT,MAXJT
            JCNT=(J-MINJT)+1
            JAUX = LOCJT+J
            IJ = ICNT+(JCNT-1)*NIT

!          COPY DENSITY
            DEN = -gXX(JAUX,IAUX)

            FORALL(IXYZ=1:3)
               VAL(IXYZ) = DEN*D22(IJ,IXYZ)
               DEL(IXYZ) = DEL(IXYZ) + VAL(IXYZ)
               DEK(IXYZ) = DEK(IXYZ) - VAL(IXYZ)
            ENDFORALL

         END DO
      END DO

!$omp critical
      DO IXYZ=1,3
        DE(IXYZ,LAT) = DE(IXYZ,LAT) + DEL(IXYZ)
        DE(IXYZ,KAT) = DE(IXYZ,KAT) + DEK(IXYZ)
      ENDDO
!$omp end critical






      RETURN

      END !***********************************************************





      SUBROUTINE RIMP2_D22_COMBINED_DENS                                &
                (DE,D22, gXX,                                           &
                 TI,TJ,W,NAUXBAS,                                       &
                 LAT,KAT,                                               &
                 LOCI,LOCJ, MINI,MAXI, MINJ,MAXJ, NI,NJ,                &
                 LOCIS,LOCJS, MINIS,MAXIS, MINJS,MAXJS, NIS,NJS,        &
                 SPHAUX,KANDL)

      use omp_lib
      use Rimp2_Shared_Data, only: MXATM

      implicit double precision(a-h,o-z)

      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,                    &
                      ZAN(MXATM),C(3,MXATM),IAN(MXATM)


!!!!! OUTPUT
      double precision :: DE(3,NAT)


!!!!! INPUT
      double precision :: D22(NI*NJ,3)
      double precision :: gXX(NAUXBAS,NAUXBAS)
      double precision :: TI(NI*NI)
      double precision :: TJ(NJ*NJ)
      double precision :: W(NI*NJ)
      logical :: SPHAUX,KANDL


!!!!! LOCAL DATA
      double precision :: VAL(3)


!!!!! CARTESIAN TO SPHERICAL AUX SHELL
      IF(SPHAUX) THEN
         MINIT = MINIS
         MAXIT = MAXIS
         MINJT = MINJS
         MAXJT = MAXJS
         LOCIT = LOCIS
         LOCJT = LOCJS
         NIT = NIS
         NJT = NJS

         DO II=1,3
            CALL DGEMM('N','N', NI,NJS, NJ,                             &
                 1.0D00, D22(1,II),NI,                                  &
                         TJ,NJ,                                         &
                 0.0D00, W,NI)
            CALL DGEMM('T','N', NIS,NJS, NI,                            &
                 1.0D00, TI,NI,                                         &
                         W,NI,                                          &
                 0.0D00, D22(1,II),NIS)
         ENDDO

      ELSE
         MINIT = MINI
         MAXIT = MAXI
         MINJT = MINJ
         MAXJT = MAXJ
         LOCIT = LOCI
         LOCJT = LOCJ
         NIT = NI
         NJT = NJ
      END IF

!!!!! D22 * DENSITY
      JMAX=MAXJT
      DO I = MINIT, MAXIT
         ICNT=(I-MINIT)+1
         IAUX = LOCIT+I

         IF(KANDL) JMAX=I

         DO J=MINJT,MAXJT
            JCNT=(J-MINJT)+1
            JAUX = LOCJT+J
            IJ = ICNT+(JCNT-1)*NIT

!!!!! COPY DENSITY
            DEN = -gXX(JAUX,IAUX)
            ! IF(JAUX.EQ.IAUX) THEN
            !    DEN = -gXX(JAUX,IAUX)
            ! ELSE
            !    DEN = -gXX(JAUX,IAUX)-gXX(JAUX,IAUX)
            ! ENDIF

            FORALL(IXYZ=1:3)
               VAL(IXYZ) = DEN*D22(IJ,IXYZ)
               DE(IXYZ,LAT) = DE(IXYZ,LAT) + VAL(IXYZ)
               DE(IXYZ,KAT) = DE(IXYZ,KAT) - VAL(IXYZ)
            ENDFORALL

         END DO
      END DO

      RETURN

      END !***********************************************************




!>*module rimp2grd   *deck RIMP2_D32GRD
!>
!>     @brief   combine 3-2ERI gradient with 3-2PDM
!>     @author  buu
!>     @date    sept 8, 2018
!>
!>     @detail  combine 3-2ERI gradient with 3-2PDM
!>
      SUBROUTINE RIMP2_D32GRD                                           &
               (DE_DDI, GNNL,                                           &
                LddiAuxShStart,LddiAuxShEnd,                            &
                LddiAuxStart,LddiAuxEnd,                                &
                NBASIS, NCOR,NACT,NOCC,NVIR,NBF,NAUXBAS)

      use omp_lib
      use Rimp2_Shared_Data

      implicit double precision(a-h,o-z)

      COMMON /AUXBAS/ EXAUX(MXAXGTOT),CAUXS(MXAXGTOT),CAUXP(MXAXGTOT),  &
                      CAUXD(MXAXGTOT),CAUXF(MXAXGTOT),CAUXG(MXAXGTOT),  &
                      CAUXH(MXAXGTOT),CAUXI(MXAXGTOT),                  &
                      KAUXST(MXAUXSH),KAUXAT(MXAUXSH),KAUXTY(MXAUXSH),  &
                      KAUXNG(MXAUXSH),KAUXLO(MXAUXSH),KAUXMI(MXAUXSH),  &
                      KAUXMX(MXAUXSH),NAUXSH
      COMMON /NSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),      &
                      CF(MXGTOT),CG(MXGTOT),CH(MXGTOT),CI(MXGTOT),      &
                      KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),   &
                      KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,                    &
                      ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /OUTPUT/ NPRINT,ITOL,ICUT,NORMF,NORMP,NOPK

!!!!! OUTPUT
      double precision :: DE_DDI(3,NAT)

!!!!! INPUT
      double precision GNNL(NBF,NBF,LddiAuxStart:LddiAuxEnd)

!!!!! LOCAL DATA
      double precision,allocatable,save ::                              &
              T1(:),D32(:,:),TL(:),XX1(:),XX2(:),DE_OMP(:,:)
!$omp threadprivate(T1,D32,TL,XX1,XX2,DE_OMP)


      logical :: GOPARR,DSKWRK,MASWRK, IANDKAT, IANDJAT, IANDJ
      double precision CC(3,3)

      INTEGER,allocatable :: ompISE(:,:)

                     !!! //// ('') (oo) \\\\ !!!

      logical :: lomp

!!!!! MXBSH: MAX BASIS FUNCTIONS IN AN AO/AUX SHELL
      MXBSH=MaxAuxBfx
      IF(MaxAuxAng < MaxAtmAng) MXBSH=MaxATMBFX
      MaxG2=MXBSH*MXBSH
      MaxG3=MXBSH*MXBSH*MXBSH

      MaxT12=MAX(NBASIS*NBASIS, MXBSH*MXBSH*MXBSH)
      ISIZE = MAXVAL(KAUXNG)


!!!!! OMP THREADS
      Nthreads_env=omp_get_max_threads()
      Nthreads=MIN0(Nthreads_env,LddiAuxShEnd-LddiAuxShStart+1)

!!!!! OMP WORK DISTRIBUTING FOR LddiAuxShStart:LddiAuxShEnd LOOP
      ALLOCATE(ompISE(0:Nthreads-1,2))
      call rimp2_WSA(ompISE,Nthreads,LddiAuxShStart,LddiAuxShEnd,LOMP)
      IF(BUG) WRITE(*,'(A40,2I5)') "NUM_THREAD SUB RIMP2_D32GRD", Nthreads, Nthreads_env


!!!!! TURN OFF DYNAMIC THREADS TO KEEP threadprivate DATA CONSISTENT
      CALL OMP_SET_DYNAMIC(.FALSE.)


!!!!! START OMP REGION
!$OMP PARALLEL NUM_THREADS(Nthreads)                                    &
!$omp default(none)                                                     &
!$omp shared(DE_DDI)                                                    &
!$omp shared(ompISE,Nthreads)                                           &
!$omp shared(C,KAUXAT,KAUXTY,KAUXMI,KAUXMX,KAUXNG,KAUXST,KAUXLO)        &
!$omp shared(FlgSphAux,LOCSPH)                                          &
!$omp shared(NAT,NSHELL,KATOM,KTYPE,KMIN,KMAX,KNG,KSTART,KLOC)          &
!$omp shared(EX,EXAUX)                                                  &
!$omp shared(GNNL,NAUXBAS,NBF,MaxT12,MXBSH,ISIZE,MaxG2,MaxG3)           &
!$omp private(IthreadID,LompAuxShStart,LompAuxShEnd)                    &
!$omp private(KAUXSH,KAT,CC,LANGK,MINK,MAXK,NK,NGK,L1,LOCK)             &
!$omp private(NKS,MINKS,MAXKS,LOCKS,MINKT,MAXKT,LOCKT,NKT,ISC,IEC)      &
!$omp private(ISH,IAT,LANGI,MINI,MAXI,NI,NGI,I1,LOCI,IANDKAT)           &
!$omp private(JAT,LANGJ,MINJ,MAXJ,NJ,NGJ,J1,LOCJ)                       &
!$omp private(IANDJ,IANDJAT,LANGIDER,LANGJDER,LANGKDER)


!!!!! ALLOCATE ARRAYS
      ALLOCATE(D32(MaxG3,6))
      ALLOCATE(T1(MaxT12))
      ALLOCATE(TL(MaxG2))
      ALLOCATE(XX1(ISIZE))
      ALLOCATE(XX2(ISIZE*28))

      ALLOCATE(DE_OMP(3,NAT))
      DE_OMP=0.0D00

      IthreadID=omp_get_thread_num()
      LompAuxShStart=ompISE(IthreadID,1)
      LompAuxShEnd=ompISE(IthreadID,2)


!!!!! LOOP AUX &AO SHELLS
      ! DO KAUXSH=1,NAUXSH
      ! DO KAUXSH=LauxSHStart,LauxSHEnd
      ! DO KAUXSH=LddiAuxShStart,LddiAuxShEnd
      DO KAUXSH=LompAuxShStart,LompAuxShEnd

            KAT = KAUXAT(KAUXSH)
            CC(1,3) = C(1,KAT)
            CC(2,3) = C(2,KAT)
            CC(3,3) = C(3,KAT)
            LANGK = KAUXTY(KAUXSH)-1
            MINK  = KAUXMI(KAUXSH)
            MAXK  = KAUXMX(KAUXSH)
            NK    = MAXK-MINK+1
            NGK   = KAUXNG(KAUXSH)
            L1    = KAUXST(KAUXSH)
            LOCK  = KAUXLO(KAUXSH)-MINK

!!!!! CARTESIAN TO SPHERICAL (AUX)
            IF(FlgSphAux) THEN
                  CALL RIMP2_RIMP2CSTRM(TL,LANGK,MINK,NK,MINKS,MAXKS)
                  LOCKS = LOCSPH(KAUXSH)-MINKS
                  NKS = MAXKS-MINKS+1
                  MINKT = MINKS
                  MAXKT = MAXKS
                  LOCKT = LOCKS
            ELSE
                  MINKT = MINK
                  MAXKT = MAXK
                  LOCKT = LOCK
            ENDIF

            NKT=MAXKT-MINKT+1

            ISC = LOCKT+MINKT
            IEC = ISC+NKT-1

!!!!! I-SHELL
            DO ISH = 1, NSHELL
                  IAT = KATOM(ISH)
                  CC(1,1) = C(1,IAT)
                  CC(2,1) = C(2,IAT)
                  CC(3,1) = C(3,IAT)
                  LANGI = KTYPE(ISH)-1
                  MINI  = KMIN(ISH)
                  MAXI  = KMAX(ISH)
                  NI    = KMAX(ISH)-KMIN(ISH)+1
                  NGI   = KNG(ISH)
                  I1    = KSTART(ISH)
                  LOCI  = KLOC(ISH)-MINI
                  IANDKAT = (IAT.EQ.KAT)

!!!!! J-SHELL
                  ! DO JSH = 1, NSHELL
                  DO JSH = 1, ISH
                        JAT = KATOM(JSH)
                        CC(1,2) = C(1,JAT)
                        CC(2,2) = C(2,JAT)
                        CC(3,2) = C(3,JAT)
                        LANGJ = KTYPE(JSH)-1
                        MINJ  = KMIN(JSH)
                        MAXJ  = KMAX(JSH)
                        NJ    = KMAX(JSH)-KMIN(JSH)+1
                        NGJ   = KNG(JSH)
                        J1    = KSTART(JSH)
                        LOCJ  = KLOC(JSH)-MINJ
                        IANDJAT=IAT.EQ.JAT
                        IANDJ=ISH.EQ.JSH

                        IF(IANDKAT.AND.IANDJAT) CYCLE

!!!!! INCREASE ANGULAR MOMENTUM BY 1
                        LANGIDER=LANGI+1
                        IF(IANDJAT) THEN
                           LANGJDER=LANGJ+1
                           LANGKDER=LANGK
                        ELSE
                           LANGJDER=LANGJ
                           LANGKDER=LANGK+1
                        ENDIF

!!!!! FORM D32
                        CALL RIMP2_X_FORM_D32                           &
                            (D32, KAUXSH,ISH,JSH,                       &
                             IANDJAT,IANDKAT,                           &
                             CC, EX(I1),EX(J1),EXAUX(L1),               &
                             LANGI,LANGJ,LANGK,                         &
                             LANGIDER,LANGJDER,LANGKDER,                &
                             NI,NJ,NK,                                  &
                             NGI,NGJ,NGK,                               &
                             MINI,MAXI, MINJ,MAXJ, MINK,MAXK,           &
                             XX1,XX2)

!!!!! D32 * DENSITY
                        CALL RIMP2_D32_TIMES_DENS                       &
                            (DE_OMP,                                    &
                             D32,GNNL(1,1,ISC), IAT,JAT,KAT,            &
                             TL,T1,MaxG2,MaxG3,                         &
                             LOCI,LOCJ,LOCKT,                           &
                             MINI,MAXI, MINJ,MAXJ, MINK,MAXK,           &
                             NI,NJ,NK, NKT,                             &
                             MINKS,MAXKS,NKS,                           &
                             FlgSphAux,IANDJAT,IANDKAT,IANDJ,           &
                             NBF,NAUXBAS)
                  ENDDO
            ENDDO

      ENDDO

!$omp critical
      DE_DDI=DE_DDI+DE_OMP
!$omp end critical

      DEALLOCATE(D32,T1,TL,XX1,XX2,DE_OMP)

!$OMP END PARALLEL

      RETURN

      END !***********************************************************







      SUBROUTINE RIMP2_RIMP2CSTRM(TR,IT,MINI,NI,MINIS,MAXIS)

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      DOUBLE PRECISION TR(NI,NI)

      COMMON /SPHERI/ PSHELL(3,3),DSHELL(6,6), &
                      FSHELL(10,10),GSHELL(15,15), &
                      PIHELL(3,3),DIHELL(6,6), &
                      FIHELL(10,10),GIHELL(15,15)
      common /sphehi/ hshell(21,21),aishell(28,28), &
                      hihell(21,21),aiihell(28,28)

      IF(IT .EQ. 0) THEN
         MINIS = 1
         MAXIS = 1
         TR(1,1) = 1.0D00
      ELSE IF(IT .EQ. 1) THEN
         IF(MINI .EQ. 1) THEN
            MINIS = 1
            MAXIS = 4
            TR(1,1) = 1.0D00
            TR(2,1) = 0.0D00
            TR(3,1) = 0.0D00
            TR(4,1) = 0.0D00
            TR(1,2) = 0.0D00
            TR(1,3) = 0.0D00
            TR(1,4) = 0.0D00
            DO I = 2, 4
               DO J = 2, 4
                  TR(I,J) = PSHELL(I,J)
               END DO
            END DO
         ELSE
            MINIS = 2
            MAXIS = 4
            DO I = 1, 3
               DO J = 1, 3
                  TR(I,J) = PSHELL(I,J)
               END DO
            END DO
         END IF
      ELSE IF(IT .EQ. 2) THEN
         MINIS = 5
         MAXIS = 9
         DO I = 1, 6
            DO J = 1, 5
               TR(I,J) = DSHELL(I,J)
            END DO
         END DO
      ELSE IF(IT .EQ. 3) THEN
         MINIS = 10
         MAXIS = 16
         DO I = 1, 10
            DO J = 1, 7
               TR(I,J) = FSHELL(I,J)
            END DO
         END DO
      ELSE IF(IT .EQ. 4) THEN
         MINIS = 17
         MAXIS = 25
         DO I = 1, 15
            DO J = 1, 9
               TR(I,J) = GSHELL(I,J)
            END DO
         END DO
      ELSE IF(IT .EQ. 5) THEN
         MINIS = 26
         MAXIS = 36
         DO I = 1,21
            DO J = 1, 11
               TR(I,J) = HSHELL(I,J)
            END DO
         END DO
      ELSE IF(IT .EQ. 6) THEN
         MINIS = 37
         MAXIS = 49
         DO I = 1,28
            DO J = 1, 13
               TR(I,J) = AISHELL(I,J)
            END DO
         END DO
      ELSE
         WRITE(6,8000)
         CALL ABRT
      END IF
      RETURN
 8000 FORMAT(/1X,'*** ERROR ***'/ &
              1X,'RI BASIS WITH MOMENTUM LARGER THAN I', &
                 ' REQUESTED')

      END

!
      SUBROUTINE RIMP2_X_FORM_D32                                       &
                (D32, KAUXSH,ISH,JSH,                                   &
                 IANDJAT,IANDKAT,                                       &
                 CC, EXI,EXJ,EXK,                                       &
                 LI,LJ,LK,                                              &
                 LANGIDER,LANGJDER,LANGKDER,                            &
                 NI,NJ,NK,                                              &
                 NGI,NGJ,NGK,                                           &
                 MINI,MAXI, MINJ,MAXJ, MINK,MAXK,                       &
                 AKINVV,DK)
!    OPENMP LIB
      use omp_lib
!
      use Rimp2_Shared_Data
!
      implicit double precision(a-h,o-z)
!
      PARAMETER (MXLTOT=6)
      PARAMETER (ITOLRI=10)
      PARAMETER (RLN10=2.30258D00)
!
      COMMON /ROOT  / XX,U(13),W(13),NROOTS
!$omp threadprivate(/ROOT  /)

!
      double precision D32(NI*NJ*NK,6)
!
      double precision CC(3,3), EXI(NGI),EXJ(NGJ),EXK(NGK)
      double precision AKINVV(NGK),DK(28,NGK)
      double precision G2(2730),G3(10290),CAB(3,0:MXLTOT),              &
                       CNK(28),DB(784),P(3),PAV(3),CPV(3)
      double precision A1(3,13),B1(13),C1(13),                          &
                       A2(3,13),B2(13),F00(13)
      INTEGER LLI(3,28),LLJ(3,28),LLK(3,28)
      LOGICAL IANDJAT,IANDKAT

!
      TOL = ITOLRI*RLN10
!
      LB = LANGIDER + LANGJDER
      LTOT = LB + LANGKDER
      NROOTS = LTOT/2+1
!
      RAB2 = 0.0D00
      CAB = 0.0D00
      DO IA = 1, 3
         VAL = CC(IA,1)-CC(IA,2)
         RAB2 = RAB2+VAL*VAL
         DO J = 0, LJ
            CAB(IA,J) = VAL**J
         END DO
      END DO
!
      DO J = 0, LJ
         JJ = J*(J+1)/2+1
         CNK(JJ) = 1.0D00
         JJ = (J+1)*(J+2)/2
         CNK(JJ) = 1.0D00
         IF(J >= 2) THEN
            DO K = 1, J-1
               JJ = J*(J+1)/2+K+1
               JA = (J-1)*J/2+K
               JB = (J-1)*J/2+K+1
               CNK(JJ) = CNK(JA)+CNK(JB)
            END DO
         END IF
      END DO

!    ANGULAR MOMENTUM COMPONENT
      DO I = MINI, MAXI
         ICNT=(I-MINI)+1
         LLI(1,ICNT) = LX(I)
         LLI(2,ICNT) = LY(I)
         LLI(3,ICNT) = LZ(I)
      END DO
!
      DO J = MINJ, MAXJ
         JCNT=(J-MINJ)+1
         LLJ(1,JCNT) = LX(J)
         LLJ(2,JCNT) = LY(J)
         LLJ(3,JCNT) = LZ(J)
      END DO
!
      DO K = MINK, MAXK
         KCNT=(K-MINK)+1
         LLK(1,KCNT) = LX(K)
         LLK(2,KCNT) = LY(K)
         LLK(3,KCNT) = LZ(K)
      END DO

!    DENSITY FACTOR
      DO KG = 1, NGK
         AKINVV(KG) = 1.0D00/EXK(KG)
         DO K=1,NK
            DK(K,KG)=DAux(KAUXSH,KG,K)
         ENDDO
      ENDDO

!    ZERO OUT INTEGRAL ARRAY
      D32=0.0D00

!    LOOP OVER GAUSSIANS
      DO IG = 1, NGI ! I-SHELL (AO)

         AI = EXI(IG)

         AIRAB = AI*RAB2
!
         DO JG = 1, NGJ ! J-SHELL (AO)

            AJ = EXJ(JG)
            AB = AI+AJ
            ABINV = 1.0D00/AB

            DUM = AJ*ABINV*AIRAB

            IF (DUM > TOL) CYCLE

            EXPB = EXP(-DUM)
!
            FORALL(II=1:3)
               P(II) = (AI*CC(II,1)+AJ*CC(II,2))*ABINV
               PAV(II) = P(II)-CC(II,1)
               CPV(II) = CC(II,3)-P(II)
            ENDFORALL
!
            RR = CPV(1)*CPV(1)+CPV(2)*CPV(2)+CPV(3)*CPV(3)

!          COMBINE IJ-SHELL DENSITY
            DUM0 = ABINV*EXPB
            IJ = 0
            DO I=1,NI
               DUM2 = DUM0*DAtm(ISH,IG,I)
               DO J=1,NJ
                  IJ = IJ+1
                  DB(IJ) = DUM2*DAtm(JSH,JG,J)
               ENDDO
            ENDDO !I
!
            DO KG = 1, NGK ! K-SHELL (AUX)

               AK = EXK(KG)
               RHO = AB*AK/(AB+AK)
               XX = RHO*RR

!             RYS ROOTS AND WEIGHT
               SELECT CASE(NROOTS)
                  CASE(:3)
                     CALL RT123
                  CASE(4)
                     CALL ROOT4
                  CASE(5)
                     CALL ROOT5
                  CASE(6:)
                     CALL ROOT6
               END SELECT

!             IXYZ: G2(XYZ,IR,LI+LJ,LK) :: LB = LI+LJ
               CALL I32XYZ_1(G2, U,W,PAV,CPV,A1,B1,C1,                  &
                    A2,B2,F00,AB,AK,LB,LANGKDER,NROOTS)

!             IXYZ: G3(XYZ,IR,LI,LJ,LK) __
               CALL I32XYZ_2(G3, G2,                                    &
                    CAB,LANGIDER,LANGJDER,LANGKDER,LB,NROOTS)

!             IAT .EQ. KAT
               IF(IANDKAT) THEN
                  CALL FORM_IANDKAT_D32                                 &
                      (D32, G3, DB,DK(1,KG),                            &
                       AI,AJ,AK,                                        &
                       LLI,LLJ,LLK,                                     &
                       LANGIDER,LANGJDER,LANGKDER,                      &
                       NI,NJ,NK, NROOTS)

!             IAT .EQ. JAT
               ELSEIF(IANDJAT) THEN
                  CALL FORM_IANDJAT_D32                                 &
                      (D32, G3, DB,DK(1,KG),                            &
                       AI,AJ,AK,                                        &
                       LLI,LLJ,LLK,                                     &
                       LANGIDER,LANGJDER,LANGKDER,                      &
                       NI,NJ,NK, NROOTS)

!             IAT .NE. JAT .NE. KAT
               ELSE
                  CALL FORM_D32                                         &
                      (D32, G3, DB,DK(1,KG),                            &
                       AI,AJ,AK,                                        &
                       LLI,LLJ,LLK,                                     &
                       LANGIDER,LANGJDER,LANGKDER,                      &
                       NI,NJ,NK, NROOTS)
               ENDIF

            END DO !KG
         END DO !JG
      END DO !IG

      RETURN

      END !***********************************************************














!
      SUBROUTINE FORM_IANDJAT_D32(D32, G3,  DB,DK,                      &
                 AI,AJ,AK, LLI,LLJ,LLK,                                 &
                 LI,LJ,LK, NI,NJ,NK, NROOTS)
!
      implicit double PRECISION(A-H,O-Z)
!
      double precision D32(NI*NJ*NK,3)

      double precision G3(3,NROOTS,0:LI,0:LJ,0:LK)
      double precision DB(NI*NJ),DK(NK)

      integer LLI(3,NI),LLJ(3,NJ),LLK(3,NK)
!
      double precision VAL(3)
      double precision DIXYZ(18,3)
      double precision XYZ(18,3) ! HOFFULLY NROOTS <= 18

!
      PARAMETER (PI252=34.986836655250D00)

!
      AI2=2.0D00*AI
      AJ2=2.0D00*AJ
      AK2=2.0D00*AK

      IJ = 0
      IJK = 0
      DO I = 1, NI
         NIX = LLI(1,I)
         NIY = LLI(2,I)
         NIZ = LLI(3,I)
         DO J = 1, NJ
            NJX = LLJ(1,J)
            NJY = LLJ(2,J)
            NJZ = LLJ(3,J)
            IJ = IJ+1

            FACTOR=PI252*DB(IJ)

            DO K = 1, NK
               NKX = LLK(1,K)
               NKY = LLK(2,K)
               NKZ = LLK(3,K)

               DO IR=1,NROOTS
                  XYZ(IR,1)=G3(2,IR,NIY,NJY,NKY)*G3(3,IR,NIZ,NJZ,NKZ)
                  XYZ(IR,2)=G3(1,IR,NIX,NJX,NKX)*G3(3,IR,NIZ,NJZ,NKZ)
                  XYZ(IR,3)=G3(1,IR,NIX,NJX,NKX)*G3(2,IR,NIY,NJY,NKY)
               ENDDO

               IJK = IJK+1

               VAL=0.0D00
               DIXYZ=0.0D00

               DO II=1,3
!                IAT CENTER__
                  IF(LLI(II,I) == 0) THEN
                     DO IR=1,NROOTS
                        DIXYZ(IR,II)=                                   &
                           AI2*G3(II,IR,LLI(II,I)+1,LLJ(II,J),LLK(II,K))
                     ENDDO
                  ELSE
                     DO IR=1,NROOTS
                        DIXYZ(IR,II)=                                   &
                           AI2*G3(II,IR,LLI(II,I)+1,LLJ(II,J),LLK(II,K))&
                    -LLI(II,I)*G3(II,IR,LLI(II,I)-1,LLJ(II,J),LLK(II,K))
                     ENDDO
                  ENDIF

!                ADDING JAT CENTER DERIVATIVE __
                  IF(LLJ(II,J) == 0) THEN
                     DO IR=1,NROOTS
                        DIXYZ(IR,II) = DIXYZ(IR,II)                     &
                         + AJ2*G3(II,IR,LLI(II,I),LLJ(II,J)+1,LLK(II,K))
                     ENDDO
                  ELSE
                     DO IR=1,NROOTS
                        DIXYZ(IR,II) = DIXYZ(IR,II)                     &
                         + AJ2*G3(II,IR,LLI(II,I),LLJ(II,J)+1,LLK(II,K))&
                    -LLJ(II,J)*G3(II,IR,LLI(II,I),LLJ(II,J)-1,LLK(II,K))
                     ENDDO
                  ENDIF
!
                  DO IR=1,NROOTS
                     VAL(II) = VAL(II) + DIXYZ(IR,II)*XYZ(IR,II)
                  ENDDO
!
                  D32(IJK,II) = D32(IJK,II) + FACTOR*DK(K)*VAL(II)

               ENDDO !II

            END DO !NK
         END DO !NJ
      END DO !NI

      RETURN

      END !***********************************************************



!
      SUBROUTINE FORM_IANDKAT_D32(D32, G3,  DB,DK,                      &
                 AI,AJ,AK, LLI,LLJ,LLK,                                 &
                 LI,LJ,LK,NI,NJ,NK,NROOTS)
!
      implicit double PRECISION(A-H,O-Z)
!
      double precision D32(NI*NJ*NK,3)

      double precision G3(3,NROOTS,0:LI,0:LJ,0:LK)
      double precision DB(NI*NJ),DK(NK)

      integer LLI(3,NI),LLJ(3,NJ),LLK(3,NK)
!
      double precision VAL(3)
      double precision DIXYZ(18,3)
      double precision XYZ(18,3) ! HOFFULLY NROOTS <= 18

!
      PARAMETER (PI252=34.986836655250D00)

!
      AI2=2.0D00*AI
      AJ2=2.0D00*AJ
      AK2=2.0D00*AK

      IJ = 0
      IJK = 0
      DO I = 1, NI
         NIX = LLI(1,I)
         NIY = LLI(2,I)
         NIZ = LLI(3,I)
         DO J = 1, NJ
            NJX = LLJ(1,J)
            NJY = LLJ(2,J)
            NJZ = LLJ(3,J)
            IJ = IJ+1

            FACTOR=PI252*DB(IJ)

            DO K = 1, NK
               NKX = LLK(1,K)
               NKY = LLK(2,K)
               NKZ = LLK(3,K)

               DO IR=1,NROOTS
                  XYZ(IR,1)=G3(2,IR,NIY,NJY,NKY)*G3(3,IR,NIZ,NJZ,NKZ)
                  XYZ(IR,2)=G3(1,IR,NIX,NJX,NKX)*G3(3,IR,NIZ,NJZ,NKZ)
                  XYZ(IR,3)=G3(1,IR,NIX,NJX,NKX)*G3(2,IR,NIY,NJY,NKY)
               ENDDO

               IJK = IJK+1

               VAL=0.0D00
               DIXYZ=0.0D00

               DO II=1,3
!                IAT CENTER__
                  IF(LLI(II,I) == 0) THEN
                     DO IR=1,NROOTS
                        DIXYZ(IR,II)=                                   &
                           AI2*G3(II,IR,LLI(II,I)+1,LLJ(II,J),LLK(II,K))
                     ENDDO
                  ELSE
                     DO IR=1,NROOTS
                        DIXYZ(IR,II)=                                   &
                           AI2*G3(II,IR,LLI(II,I)+1,LLJ(II,J),LLK(II,K))&
                    -LLI(II,I)*G3(II,IR,LLI(II,I)-1,LLJ(II,J),LLK(II,K))
                     ENDDO
                  ENDIF

!                ADDING KAT CENTER DERIVATIVE __
                  IF(LLK(II,K) == 0) THEN
                     DO IR=1,NROOTS
                        DIXYZ(IR,II) = DIXYZ(IR,II)                     &
                         + AK2*G3(II,IR,LLI(II,I),LLJ(II,J),LLK(II,K)+1)
                     ENDDO
                  ELSE
                     DO IR=1,NROOTS
                        DIXYZ(IR,II) = DIXYZ(IR,II)                     &
                         + AK2*G3(II,IR,LLI(II,I),LLJ(II,J),LLK(II,K)+1)&
                    -LLK(II,K)*G3(II,IR,LLI(II,I),LLJ(II,J),LLK(II,K)-1)
                     ENDDO
                  ENDIF
!
                  DO IR=1,NROOTS
                     VAL(II)=VAL(II)+DIXYZ(IR,II)*XYZ(IR,II)
                  ENDDO
!
                  D32(IJK,II)=D32(IJK,II)+FACTOR*DK(K)*VAL(II)

               ENDDO !II

            END DO !NK
         END DO !NJ
      END DO !NI

      RETURN

      END !***********************************************************



!
      SUBROUTINE FORM_D32(D32, G3,  DB,DK,                              &
                 AI,AJ,AK, LLI,LLJ,LLK,                                 &
                 LI,LJ,LK,NI,NJ,NK,NROOTS)
!
      implicit double PRECISION(A-H,O-Z)
!
      double precision D32(NI*NJ*NK,3)

      double precision G3(3,NROOTS,0:LI,0:LJ,0:LK)
      double precision DB(NI*NJ),DK(NK)

      integer LLI(3,NI),LLJ(3,NJ),LLK(3,NK)
!
      double precision VAL(3)
      double precision DIXYZ(18,3)
      double precision XYZ(18,3) ! HOFFULLY NROOTS <= 18

!
      PARAMETER (PI252=34.986836655250D00)

!
      AI2=2.0D00*AI
      AJ2=2.0D00*AJ
      AK2=2.0D00*AK

      IJ = 0
      IJK = 0
      DO I = 1, NI
         NIX = LLI(1,I)
         NIY = LLI(2,I)
         NIZ = LLI(3,I)
         DO J = 1, NJ
            NJX = LLJ(1,J)
            NJY = LLJ(2,J)
            NJZ = LLJ(3,J)
            IJ = IJ+1

            FACTOR=PI252*DB(IJ)

            DO K = 1, NK
               NKX = LLK(1,K)
               NKY = LLK(2,K)
               NKZ = LLK(3,K)

               DO IR=1,NROOTS
                  XYZ(IR,1)=G3(2,IR,NIY,NJY,NKY)*G3(3,IR,NIZ,NJZ,NKZ)
                  XYZ(IR,2)=G3(1,IR,NIX,NJX,NKX)*G3(3,IR,NIZ,NJZ,NKZ)
                  XYZ(IR,3)=G3(1,IR,NIX,NJX,NKX)*G3(2,IR,NIY,NJY,NKY)
               ENDDO

               IJK = IJK+1

               VAL=0.0D00
               DIXYZ=0.0D00

!             IAT CENTER__
               DO II=1,3
                  IF(LLI(II,I) == 0) THEN
                     DO IR=1,NROOTS
                        DIXYZ(IR,II)=                                   &
                           AI2*G3(II,IR,LLI(II,I)+1,LLJ(II,J),LLK(II,K))
                     ENDDO
                  ELSE
                     DO IR=1,NROOTS
                        DIXYZ(IR,II)=                                   &
                           AI2*G3(II,IR,LLI(II,I)+1,LLJ(II,J),LLK(II,K))&
                    -LLI(II,I)*G3(II,IR,LLI(II,I)-1,LLJ(II,J),LLK(II,K))
                     ENDDO
                  ENDIF
!
                  DO IR=1,NROOTS
                     VAL(II)=VAL(II)+DIXYZ(IR,II)*XYZ(IR,II)
                  ENDDO
!
                  D32(IJK,II)=D32(IJK,II)+FACTOR*DK(K)*VAL(II)

               ENDDO

               VAL=0.0D00
               DIXYZ=0.0D00

!             KAT CENTER __
               DO II=1,3
                  IF(LLK(II,K) == 0) THEN
                     DO IR=1,NROOTS
                        DIXYZ(IR,II) = DIXYZ(IR,II)                     &
                         + AK2*G3(II,IR,LLI(II,I),LLJ(II,J),LLK(II,K)+1)
                     ENDDO
                  ELSE
                     DO IR=1,NROOTS
                        DIXYZ(IR,II) = DIXYZ(IR,II)                     &
                         + AK2*G3(II,IR,LLI(II,I),LLJ(II,J),LLK(II,K)+1)&
                    -LLK(II,K)*G3(II,IR,LLI(II,I),LLJ(II,J),LLK(II,K)-1)
                     ENDDO
                  ENDIF
!
                  DO IR=1,NROOTS
                     VAL(II)=VAL(II)+DIXYZ(IR,II)*XYZ(IR,II)
                  ENDDO
!
                  D32(IJK,II+3)=D32(IJK,II+3)+FACTOR*DK(K)*VAL(II)

               ENDDO !II

            END DO !NK
         END DO !NJ
      END DO !NI

      RETURN

      END !***********************************************************


!
      SUBROUTINE RIMP2OMP_D32_TIMES_DENS                                &
                (DE, D32,GNNL, IAT,JAT,KAT,                                 &
                 TL,W,MaxT12,MXBSH,                                     &
                 LOCI,LOCJ,LOCKT,                                       &
                 MINI,MAXI, MINJ,MAXJ, MINK,MAXK,                       &
                 NI,NJ,NK, NKT,                                         &
                 MINKS,MAXKS,NKS,                                       &
                 SPHAUX, IANDJAT,IANDKAT,                               &
                 NBF,NAUXBAS)

!    OPENMP LIB
      use omp_lib

!    GLOBAL DATA
      use Rimp2_Shared_Data,only: MXATM
!
      implicit double precision(A-H,O-Z)

!    GLOBAL DATA
      double precision :: DE(3,MXATM)
!    INPUT
      double precision :: D32(NI*NJ*NK,6)
      double precision :: GNNL(NBF,NBF,NKT)

      double precision :: TL(MXBSH*MXBSH)
      double precision :: W(MaxT12)
      logical :: SPHAUX, IANDJAT,IANDKAT

!    LOCAL DATA
      double precision :: DEIX,DEIY,DEIZ
      double precision :: DEJX,DEJY,DEJZ
      double precision :: DEKX,DEKY,DEKZ


                     !!! //// ('') (oo) (**) \\\\ !!!



!    NDIM=6 FOR THREE CENTERS; 3 FOR TWO CENTER DIFFERENTIATION
      NDIM=6
      IF(IANDJAT.OR.IANDKAT) NDIM=3

!    TRANSFORMING AUXILIARY SHELL TO SPHERICAL IF SPHAUX=.TRUE.
      IF(SPHAUX) THEN
         MINKT = MINKS
         MAXKT = MAXKS
         NKT = NKS
         NKJ = NKS*NJ
!
         DO II=1,NDIM
            CALL DGEMM('T','N', NKS,NI*NJ, NK,                          &
                 1.0D00, TL,NK,                                        &
                          D32(1,II),NK,                                 &
                 0.0D00, W,NKS)
            CALL DCOPY(NKJ*NI, W,1, D32(1,II),1)
         ENDDO
      ELSE
         MINKT = MINK
         MAXKT = MAXK
         NKT = NK
         NKJ = NK*NJ
      END IF


!    INITIALIZE TEMPORARY PRIVATE GRADIENT
      DEIX=0.0D00
      DEIY=0.0D00
      DEIZ=0.0D00

      DEJX=0.0D00
      DEJY=0.0D00
      DEJZ=0.0D00

      DEKX=0.0D00
      DEKY=0.0D00
      DEKZ=0.0D00


!    COMBINE THREE-INDEX DENSITY AND INTEGRAL GRADIENT
      DO LL = 1, NKT
         DO II = 1, NI
            IBF = (MINI+II-1) + LOCI
            DO JJ = 1, NJ
               JBF = (MINJ+JJ-1) + LOCJ
               LIJ = LL + (JJ-1)*NKT + (II-1)*NKJ

!             COPY DENSITY
               DEN = GNNL(JBF,IBF, LL)

!             IAT CENTER
               valIx = DEN*D32(LIJ,1)
               valIy = DEN*D32(LIJ,2)
               valIz = DEN*D32(LIJ,3)
!
               DEIX = DEIX + valIx
               DEIY = DEIY + valIy
               DEIZ = DEIZ + valIz

!             IAT==JAT :: use INVARIANCE FOR KAT CENTER
               IF(IANDJAT) THEN
                  DEKX = DEKX - valIx
                  DEKY = DEKY - valIy
                  DEKZ = DEKZ - valIz

!             IAT==KAT :: use INVARIANCE FOR JAT CENTER
               ELSEIF(IANDKAT) THEN
                  DEJX = DEJX - valIx
                  DEJY = DEJY - valIy
                  DEJZ = DEJZ - valIz

!             IAT != JAT != KAT :: ADD KAT DER CONTRIBUTION AND
!                                   use INVARIANCE FOR JAT CENTER
               ELSE
                  valKx = DEN*D32(LIJ,4)
                  valKy = DEN*D32(LIJ,5)
                  valKz = DEN*D32(LIJ,6)
!
                  DEKX = DEKX + valKx
                  DEKY = DEKY + valKy
                  DEKZ = DEKZ + valKz
!
                  DEJX = DEJX - valIx - valKx
                  DEJY = DEJY - valIy - valKy
                  DEJZ = DEJZ - valIz - valKz
               ENDIF
            END DO
         END DO
      END DO


!    COPY TEMPORARY PRIVATE GRAD VALUE TO THE SHARED
!     GRADIENT (DE) IN THE (/GRAD  /) COMMON BLOCK
!$omp critical
      DE(1,IAT) = DE(1,IAT) + DEIX
      DE(2,IAT) = DE(2,IAT) + DEIY
      DE(3,IAT) = DE(3,IAT) + DEIZ

      DE(1,JAT) = DE(1,JAT) + DEJX
      DE(2,JAT) = DE(2,JAT) + DEJY
      DE(3,JAT) = DE(3,JAT) + DEJZ

      DE(1,KAT) = DE(1,KAT) + DEKX
      DE(2,KAT) = DE(2,KAT) + DEKY
      DE(3,KAT) = DE(3,KAT) + DEKZ
!$omp end critical


      RETURN

      END !***********************************************************



!
      SUBROUTINE RIMP2_D32_TIMES_DENS                                   &
                (DE, &
                 D32,GNNL, IAT,JAT,KAT,                                 &
                 TL,W,MaxG2,MaxG3,                                      &
                 LOCI,LOCJ,LOCKT,                                       &
                 MINI,MAXI, MINJ,MAXJ, MINK,MAXK,                       &
                 NI,NJ,NK, NKT,                                         &
                 MINKS,MAXKS,NKS,                                       &
                 SPHAUX, IANDJAT,IANDKAT,IANDJ,                         &
                 NBF,NAUXBAS)

!
      use Rimp2_Shared_Data,only: MXATM
!
      implicit double precision(A-H,O-Z)
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,                    &
                      ZAN(MXATM),C(3,MXATM),IAN(MXATM)

!
      double precision :: DE(3,NAT)
!
      double precision :: D32(NI*NJ*NK,6)

      double precision :: TL(MaxG2),W(MaxG3)
      LOGICAL :: SPHAUX, IANDJAT,IANDKAT,IANDJ
      double precision :: GNNL(NBF,NBF,NKT)

      ! double precision,allocatable :: WTMP(:)


!
      NDIM=6
      IF(IANDJAT.OR.IANDKAT) NDIM=3

!    CARTESIAN TO SPHERICAL D32 (AUX SHELL)
      IF(SPHAUX) THEN
         MINKT = MINKS
         MAXKT = MAXKS
         NKT = NKS
         NKJ = NKS*NJ
!
         DO II=1,NDIM
            CALL DGEMM('T','N', NKS,NI*NJ, NK,                          &
                 1.0D00, TL,NK,                                        &
                          D32(1,II),NK,                                 &
                 0.0D00, W,NKS)
            CALL DCOPY(NKJ*NI, W,1, D32(1,II),1)
         ENDDO
      ELSE
         MINKT = MINK
         MAXKT = MAXK
         NKT = NK
         NKJ = NK*NJ
      END IF

!
      JJEND=NJ

      DO LL = 1, NKT
!
         DO II = 1, NI
            IBF = (MINI+II-1) + LOCI
!
            IF(IANDJ) JJEND=II
            DO JJ = 1, JJEND
               JBF = (MINJ+JJ-1) + LOCJ
!
               LIJ = LL + (JJ-1)*NKT + (II-1)*NKJ
!
               DEN = GNNL(JBF,IBF, LL)

!             IAT CENTER
               valIx = DEN*D32(LIJ,1)
               valIy = DEN*D32(LIJ,2)
               valIz = DEN*D32(LIJ,3)
!
               DE(1,IAT) = DE(1,IAT) + valIx
               DE(2,IAT) = DE(2,IAT) + valIy
               DE(3,IAT) = DE(3,IAT) + valIz

!             IAT==JAT :: use INVARIANCE FOR KAT CENTER
               IF(IANDJAT) THEN
                  DE(1,KAT) = DE(1,KAT) - valIx
                  DE(2,KAT) = DE(2,KAT) - valIy
                  DE(3,KAT) = DE(3,KAT) - valIz

!             IAT==KAT :: use INVARIANCE FOR JAT CENTER
               ELSEIF(IANDKAT) THEN
                  DE(1,JAT) = DE(1,JAT) - valIx
                  DE(2,JAT) = DE(2,JAT) - valIy
                  DE(3,JAT) = DE(3,JAT) - valIz

!             IAT != JAT != KAT :: ADD KAT DER CONTRIBUTION AND
!                                   use INVARIANCE FOR JAT CENTER
               ELSE
                  valKx = DEN*D32(LIJ,4)
                  valKy = DEN*D32(LIJ,5)
                  valKz = DEN*D32(LIJ,6)
!
                  DE(1,KAT) = DE(1,KAT) + valKx
                  DE(2,KAT) = DE(2,KAT) + valKy
                  DE(3,KAT) = DE(3,KAT) + valKz
!
                  DE(1,JAT) = DE(1,JAT) - valIx - valKx
                  DE(2,JAT) = DE(2,JAT) - valIy - valKy
                  DE(3,JAT) = DE(3,JAT) - valIz - valKz
               ENDIF

            END DO
         END DO
      END DO

      RETURN

      END !***********************************************************







      SUBROUTINE RIMP2_DDIAUXSH &
                 (LddiAuxShStart,LddiAuxShEnd, IRANK,NPROC,NAUXBAS)

      use Rimp2_Shared_Data

      implicit double precision(a-h,o-z)

      COMMON /AUXBAS/ EXAUX(MXAXGTOT),CAUXS(MXAXGTOT),CAUXP(MXAXGTOT),  &
                      CAUXD(MXAXGTOT),CAUXF(MXAXGTOT),CAUXG(MXAXGTOT),  &
                      CAUXH(MXAXGTOT),CAUXI(MXAXGTOT),                  &
                      KAUXST(MXAUXSH),KAUXAT(MXAUXSH),KAUXTY(MXAUXSH),  &
                      KAUXNG(MXAUXSH),KAUXLO(MXAUXSH),KAUXMI(MXAUXSH),  &
                      KAUXMX(MXAUXSH),NAUXSH

      double precision,allocatable :: TL(:)



      MXBSH=MaxAUXBFX
      IF(MaxAUXANG < MaxATMANG) MXBSH=MaxATMBFX
      ALLOCATE(TL(MXBSH*MXBSH))

      KK=NAUXBAS/NPROC

      LddiAuxShStart=0
      LddiAuxShEnd=0

      DO II=0,IRANK

         Nchunks=0

         IF (LddiAuxShEnd.EQ.NAUXSH) THEN
            LddiAuxShStart=LddiAuxShEnd
            EXIT
         ENDIF

         IStart=LddiAuxShEnd+1
         LddiAuxShStart=IStart
         DO KAUXSH=IStart,NAUXSH

               LANGK = KAUXTY(KAUXSH)-1
               MINK  = KAUXMI(KAUXSH)
               MAXK  = KAUXMX(KAUXSH)
               NK    = MAXK-MINK+1

               IF(FlgSphAux) THEN
                     CALL RIMP2CSTRM(TL,LANGK,MINK,NK,MINKS,MAXKS)
                     LOCKS = LOCSPH(KAUXSH)-MINKS
                     NKS = MAXKS-MINKS+1
                     MINKT = MINKS
                     MAXKT = MAXKS
                     LOCKT = LOCKS
               ELSE
                     MINKT = MINK
                     MAXKT = MAXK
                     LOCKT = LOCK
               ENDIF

               NKT=MAXKT-MINKT+1

               Nchunks=Nchunks+NKT

               LddiAuxShEnd=KAUXSH

               IF(Nchunks .GE. KK) EXIT

         ENDDO !KAUXSH


      ENDDO

      IF(IRANK .eq. NPROC-1) LddiAuxShEnd=NAUXSH


      END




!>*module rimp2grd   *deck RIMP2_32_DENS
!>
!>     @brief   form 3-2PDM
!>     @author  buu
!>     @date    sept 8, 2018
!>
!>     @detail  form 3-2PDM
!>

      SUBROUTINE RIMP2_32_DENS                                          &
                (GNNL, GNAL, VEC,                                       &
                 LddiAuxShStart,LddiAuxShEnd,                           &
                 LddiAuxStart,LddiAuxEnd,                               &
                 NCOR,NACT,NBF,NAUXBAS,NAUXSH)


      use omp_lib
      use Rimp2_Shared_Data,only: D_GNAX,BUG

      implicit double precision(a-h,o-z)


!!!!! OUTPUT
      double precision GNNL(NBF,NBF,LddiAuxStart:LddiAuxEnd)


!!!!! INPUT
      double precision :: GNAL(NBF,NACT,LddiAuxStart:LddiAuxEnd)
      double precision :: VEC(NBF,NBF)

!!!!! LOCAL DATA
      logical :: LOMP
      integer,allocatable :: ompISE(:,:)
      double precision,allocatable :: PMP2(:,:),GNN(:,:)

                     !!! //// ('') (oo) \\\\ !!!



!!!!! OMP THREADS
      Nthreads_env=omp_get_max_threads()
      Nthreads=Nthreads_env
      IF(BUG) WRITE(*,'(A40,2I5)') "NUM_THREAD SUB RIMP2_32_DENS", Nthreads, Nthreads_env


!!!!! START OMP PARALLEL REGION
!$OMP PARALLEL NUM_THREADS(Nthreads)                                    &
!$omp default(none)                                                     &
!$omp shared(LddiAuxStart,LddiAuxEnd)                                   &
!$omp shared(VEC,GNAL,NCOR,NACT)                                        &
!$omp shared(GNNL,NBF)                                                  &
!$omp private(KAUX,VI1,VI2,IM,IN)                                       &
!$omp private(GNN,IBF,JBF,TMP)

      ALLOCATE(GNN(NBF,NBF))

!$omp do schedule(DYNAMIC)
      DO KAUX =LddiAuxStart,LddiAuxEnd

!!!!! COMPLETE AO-TRANSFORMATION OF NON-SEPARABLE DENSITY
         CALL DGEMM('N','T', NBF,NBF, NACT,                             &
                     4.0D00, GNAL(1,1,KAUX),NBF,                        &
                             VEC(1,1+NCOR),NBF,                         &
                     0.0D00, GNNL(1,1,KAUX),NBF)

         DO IN=1,NBF
            DO IM=1,IN-1
               GNNL(IM,IN,KAUX)=GNNL(IM,IN,KAUX)+GNNL(IN,IM,KAUX)
            ENDDO
         ENDDO

      ENDDO
!$omp end do

      DEALLOCATE(GNN)

!$OMP END PARALLEL

      END !*************************************************************






!
      SUBROUTINE RIMP2_COPYI32                                          &
          (T0T, GHONDO,TL,W,NBASIS,                                     &
           LOCI,LOCJ,LOCKT,                                             &
           MINI,MAXI,MINJ,MAXJ,MINK,MAXK,                               &
           NI,NJ,NK, MINKS,MAXKS, NKS,                                  &
           SPHAUX,IANDJ,NAUXBAS)
!
      implicit double PRECISION(A-H,O-Z)
!
      double precision GHONDO(NI*NJ*NK)
      logical SPHAUX,IANDJ
      double PRECISION TL(*),W(*)
      double precision T0T(NAUXBAS,NBASIS,NI)
!
      PARAMETER (ZERO=0.0D00)
      PARAMETER (ONE=1.0D00)
!
      IF(SPHAUX) THEN
         MINKT = MINKS
         MAXKT = MAXKS
         NKT = NKS
         NKJ = NKS*NJ

         CALL DGEMM('T','N', NKS,NI*NJ, NK,                             &
              1.0D00, TL,NK, GHONDO,NK,                                &
              0.0D00, W,NKS)

         CALL DCOPY(NKJ*NI, W,1, GHONDO,1)
      ELSE
         MINKT = MINK
         MAXKT = MAXK
         NKT = NK
         NKJ = NK*NJ
      END IF

!



      ICNT = 0
      DO I = MINI, MAXI
         ICNT = ICNT+1
         IBF = LOCI+I

         JCNT = 0
         DO J = MINJ, MAXJ
            JCNT = JCNT+1
            JBF = LOCJ+J

            LL = 0
            DO L = MINKT, MAXKT
               LL = LL+1
               LAUX=LOCKT+L

               LIJ = LL+(JCNT-1)*NKT+(ICNT-1)*NKJ

               T0T(LAUX,JBF,ICNT) = GHONDO(LIJ)


            END DO
         END DO
      END DO



      ! LL = 0
      ! DO L = MINKT, MAXKT
      !    LL = LL+1
      !    LAUX=LOCKT+L

      !    ICNT = 0
      !    DO I = MINI, MAXI
      !       ICNT = ICNT+1
      !       IBF = LOCI+I

      !       JCNT = 0
      !       DO J = MINJ, MAXJ
      !          JCNT = JCNT+1
      !          JBF = LOCJ+J

      !          LIJ = LL+(JCNT-1)*NKT+(ICNT-1)*NKJ

      !          T0T(LAUX,JBF,ICNT) = GHONDO(LIJ)


      !       END DO
      !    END DO
      ! END DO

!
      RETURN
      END




      SUBROUTINE RIMP2_VVAA(BDVA, LddiActStart,LddiActEnd, &
                 NAUXBASD,NCOR,NACT,NVIR)

      use omp_lib
      use Rimp2_Shared_Data,only: D_B,D_VVAA

      implicit double precision(a-h,o-z)

      ! COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK

      double precision :: BDVA(NAUXBASD,NVIR,LddiActStart:LddiActEnd)

      double precision,allocatable :: VVA(:,:,:),BDV(:,:),BJ(:,:),VVA2(:,:,:)

      double precision,allocatable,save :: BI(:,:)
!$omp threadprivate(BI)

      integer,allocatable :: ompISE(:,:)


      NDO=NAUXBASD*(NCOR+NACT)
      NDV=NAUXBASD*NVIR
!    OMP THREAD
      Nthreads=omp_get_max_threads()

      ALLOCATE(BDV(NAUXBASD,NVIR))

      IStartRow=NVIR*NVIR*(LddiActStart-1)+1
      IEndRow=NVIR*NVIR*LddiActEnd

!       ALLOCATE(VVA(NVIR,NVIR,LddiActStart:LddiActEnd))

!       DO JACT=1,NACT
!          JJ=JACT+NCOR

!             CALL DDI_GET(D_B, 1+NDO,NDV+NDO, JJ,JJ, BDV)

! !$OMP PARALLEL DO NUM_THREADS(Nthreads)                                 &
! !$omp default(none)                                                     &
! !$omp shared(LddiActStart,LddiActEnd)                                   &
! !$omp private(II,IACT)                                                  &
! !$omp shared(NCOR,NVIR,NAUXBASD)                                        &
! !$omp shared(VVA,BDVA,BDV)

!             DO IACT=LddiActStart,LddiActEnd
!                II=IACT+NCOR
!                CALL DGEMM('T','N',NVIR,NVIR,NAUXBASD,                   &
!                            1.0D00, BDV,NAUXBASD,                        &
!                                    BDVA(1,1,IACT),NAUXBASD,             &
!                            0.0D00, VVA(1,1,IACT),NVIR)
!             ENDDO !JACT

! !$OMP END PARALLEL DO

!             CALL DDI_PUT(D_VVAA, IStartRow,IEndRow, JACT,JACT, VVA)

!       ENDDO !IACT
!       CALL DDI_SYNC(5678)




      NDO=NAUXBASD*(NCOR+NACT)
      NDV=NAUXBASD*NVIR
!    OMP THREAD
      Nthreads=omp_get_max_threads()
      ALLOCATE(BJ(NAUXBASD,NVIR))
      ALLOCATE(VVA(NVIR,NVIR,NACT))
      ALLOCATE(VVA2(NVIR,NVIR,NACT))

!$OMP PARALLEL
      ALLOCATE(BI(NAUXBASD,NVIR))
!$OMP END PARALLEL


      DO JACT=LddiActStart,LddiActEnd
         JJ=JACT+NCOR



!$OMP PARALLEL NUM_THREADS(Nthreads)                                    &
!$omp default(none)                                                     &
!$omp private(II,IACT,IA,IB)                                            &
!$omp shared(JACT,NCOR,NVIR,NAUXBASD,NDV,NDO)                           &
!$omp shared(BDVA,VVA,VVA2,D_B,BJ)

!$omp do schedule(guided)
         DO IACT=1,JACT ! NACT
            II=IACT+NCOR

            CALL DDI_GET(D_B, 1+NDO,NDV+NDO, II,II, BI)

            CALL DGEMM('T','N',NVIR,NVIR,NAUXBASD,       &
                        1.0D00, BDVA(1,1,JACT),NAUXBASD, &
                                BI,NAUXBASD,             &
                        0.0D00, VVA(1,1,IACT),NVIR)

         ENDDO !JACT
!$omp end do

!$omp do schedule(guided)
         DO IACT=1,JACT-1
            DO IB=1,NVIR
               DO IA=1,NVIR
                  VVA2(IA,IB,IACT)=VVA(IB,IA,IACT)
               ENDDO
            ENDDO
         ENDDO
!$omp end do

!$OMP END PARALLEL

         IStartRow=1
         IEndRow=NVIR*NVIR*JACT
         CALL DDI_PUT(D_VVAA, IStartRow,IEndRow, JACT,JACT, VVA)

         IStartRow=NVIR*NVIR*(JACT-1)+1
         IEndRow=NVIR*NVIR*JACT
         CALL DDI_PUT(D_VVAA, IStartRow,IEndRow, 1,JACT-1, VVA2)

      ENDDO !IACT

      CALL DDI_SYNC(5678)


!$OMP PARALLEL
      DEALLOCATE(BI)
!$OMP END PARALLEL

      DEALLOCATE(BJ,VVA,VVA2)



      END  !************************************************************



!
      SUBROUTINE RIMP2_PVV_YDV                                          &
                (PVV,YDV, BJ,QVV,TVV,EIG,                               &
                 IACT,JACT, NAUXBASD,NCOR,NACT,NVIR,NORB,NBF)

      implicit double precision(a-h,o-z)

      double precision :: PVV(NVIR,NVIR)
      double precision :: YDV(NAUXBASD,NVIR)

      double precision :: BJ(NAUXBASD,NVIR)
      double precision :: EIG(NACT+NVIR)

      double precision :: QVV(NVIR,NVIR)
      double precision :: TVV(NVIR,NVIR)




! !    TVV = 2*QVV - transpose(QVV)
!       TVV=QVV+QVV
!       TVV=TVV-transpose(QVV)

! !    TVV = TVV/Dijab
! !    QVV = QVV/Dijab
!       eij=EIG(IACT)+EIG(JACT)
!       DO IB=1,NVIR
!          eijB=eij-EIG(IB+NACT)
!          DO IA=1,NVIR
!             eijAB=eijB-EIG(IA+NACT)
!             TVV(IA,IB)=TVV(IA,IB)/eijAB
!             QVV(IA,IB)=QVV(IA,IB)/eijAB
!          ENDDO
!       ENDDO





!    TVV = TVV/Dijab
!    QVV = QVV/Dijab
      eij=EIG(IACT)+EIG(JACT)
      DO IB=1,NVIR
         eijB=eij-EIG(IB+NACT)
         DO IA=1,NVIR
            eijAB=eijB-EIG(IA+NACT)
            QVV(IA,IB)=QVV(IA,IB)/eijAB
         ENDDO
      ENDDO


!    TVV = 2*QVV - transpose(QVV)
      TVV=QVV+QVV
      TVV=TVV-transpose(QVV)



!    YDV <- BJ x transpose(TVV)
      CALL DGEMM('N','T', NAUXBASD,NVIR,NVIR,                           &
                  1.0D00, BJ,NAUXBASD,                                  &
                          TVV,NVIR,                                     &
                  1.0D00, YDV,NAUXBASD)


!    PVV <- QVV x transpose(TVV)
      CALL DGEMM('N','T', NVIR,NVIR,NVIR,                               &
                 -1.0D00, QVV,NVIR,                                     &
                          TVV,NVIR,                                     &
                  1.0D00, PVV,NVIR)

      END !***********************************************************





      SUBROUTINE RIMP2_YDV                                              &
                (YDV, BI,BJ,                                            &
                 QVV,TVV,EIG,IACT,JACT,                                 &
                 NAUXBASD,NCOR,NACT,NVIR,NORB,NBF)
!
      implicit double precision(a-h,o-z)
!
      double precision YDV(NAUXBASD,NVIR)
!
      double precision BI(NAUXBASD,NVIR)
      double precision BJ(NAUXBASD,NVIR)
      double precision QVV(NVIR,NVIR)
      double precision TVV(NVIR,NVIR)
      double precision EIG(NACT+NVIR)



!    ACTIVE OCCUPIED PAIRWISE MO ENERGY
      eij=EIG(IACT)+EIG(JACT)


!    QVV <- transpose(BI) x BJ
      CALL DGEMM('T','N',NVIR,NVIR,NAUXBASD,                            &
                  1.0D00, BI,NAUXBASD,                                  &
                           BJ,NAUXBASD,                                 &
                  0.0D00, QVV,NVIR)


!    TVV <- [2*QVV - transpose(QVV)] / Dijab
      TVV=QVV+QVV
      TVV=TVV-TRANSPOSE(QVV)
      FORALL(IA=1:NVIR,IB=1:NVIR)
         TVV(IA,IB)=TVV(IA,IB)/(eij-EIG(IA+NACT)-EIG(IB+NACT))
      ENDFORALL



!    YDV <- BJ x transpose(TVV)
      CALL DGEMM('N','T', NAUXBASD,NVIR,NVIR,                           &
                  1.0D00, BJ,NAUXBASD,                                  &
                           TVV,NVIR,                                    &
                  1.0D00, YDV,NAUXBASD)

      END !***********************************************************





!
      SUBROUTINE RIMP2_POO_WOO                                          &
                (POO,WOO, QVV,TVV,EIG,                                  &
                 IACT,JACT,NAUXBASD,NCOR,NACT,NVIR,NORB,NBF)
!
      USE Rimp2_Shared_Data,only: D_B,D_VVAA
!
      implicit double precision(a-h,o-z)

!    OUTPUT
      double precision POO(NCOR+NACT,NCOR+NACT)
      double precision WOO(NCOR+NACT,NCOR+NACT)

!    LOCAL DATA
      double precision QVV(NVIR,NVIR)
      double precision,allocatable :: VVA(:,:,:)
      double precision TVV(NVIR,NVIR)
      double precision EIG(NORB)


!    PARAMETERS
      NOCC=NCOR+NACT

      IOCC=IACT+NCOR



      ALLOCATE(VVA(NVIR,NVIR,NACT))
      CALL DDI_GET(D_VVAA, 1,NVIR*NVIR*NACT, JACT,JACT, VVA)


!    FORM WOO, POO
      DO IpACT=1,NACT
         IpOCC=IpACT+NCOR

!       POO <- (-1) * QVV/DVV * TVV
         eij=EIG(IpOCC)+EIG(JACT+NCOR)
         FORALL(IA=1:NVIR,IB=1:NVIR)
            QVV(IA,IB)=VVA(IB,IA,IpACT)/(eij-EIG(IA+NOCC)-EIG(IB+NOCC))
         ENDFORALL
         POO(IpOCC,IOCC)=POO(IpOCC,IOCC)-DDOT(NVIR*NVIR, QVV,1, TVV,1)

      ENDDO

      DEALLOCATE(VVA)




!       ALLOCATE(VVA(NVIR,NVIR,NACT))
!       CALL DDI_GET(D_VVAA, 1,NVIR*NVIR*NACT, JACT,JACT, VVA)

! !    FORM WOO, POO
!       DO IpACT=1,NACT
!          IpOCC=IpACT+NCOR

!          ! TMPWOO=0.0D00
!          TMPPOO=0.0D00
!          eij=EIG(IpOCC)+EIG(JACT+NCOR)
!          DO IB=1,NVIR
!             eijB=eij-EIG(IB+NOCC)
!             DO IA=1,NVIR
!                eijAB=eijB-EIG(IA+NOCC)
!                TMP=VVA(IB,IA,IpACT)
!                TMP=TMP*TVV(IA,IB)
!                ! TMPWOO=TMPWOO-TMP-TMP
!                TMPPOO=TMPPOO-TMP/eijAB
!             ENDDO
!          ENDDO

! ! !       WOO <- (-2) * QVV * TVV
! !          WOO(IpOCC,IOCC)=WOO(IpOCC,IOCC)+TMPWOO

! !       POO <- (-1) * QVV/DVV * TVV
!          POO(IpOCC,IOCC)=POO(IpOCC,IOCC)+TMPPOO

!       ENDDO

!       DEALLOCATE(VVA)


      END !***********************************************************





!
      SUBROUTINE RIMP2_WOO_POO_XLOV                                     &
                (WOO,POO,XLOV, BDVV,BDVC, YDV,EIG,                            &
                 IACT,NAUXBASD,NCOR,NACT,NVIR,NBF)
!
      implicit double precision(a-h,o-z)

!    OUTPUT
      double precision :: WOO(NCOR+NACT,NCOR+NACT)
      double precision :: POO(NCOR+NACT,NCOR+NACT)
      double precision :: XLOV(NCOR+NACT,NVIR)

!    INPUT
      double precision :: BDVV(NAUXBASD,NVIR,NVIR)
      double precision :: BDVC(NAUXBASD,NVIR,NCOR)
      double precision :: YDV(NAUXBASD,NVIR)
      double precision :: EIG(NCOR+NACT)



!    PARAMETERS
      NOCC=NCOR+NACT
      NXO=NAUXBASD*NOCC
      NXV=NAUXBASD*NVIR

!    XLOV(c,i) <- (2) * BC(R,a,c) * YDV(R,a,i)
      II=IACT+NCOR
      DO IA=1,NVIR
         TMP=DDOT(NAUXBASD*NVIR, BDVV(1,1,IA),1, YDV,1)
         XLOV(II,IA)=XLOV(II,IA) -TMP-TMP
      ENDDO




!    WOO(NACT,NCOR) = BI*YDV;
!    POO(NCOR,NACT) = BI*YDV/(ei - eI);
!    POO(NACT,NCOR) = POO(NCOR,NACT);
      DO ICOR=1,NCOR
         TMP=DDOT(NAUXBASD*NVIR, BDVC(1,1,ICOR),1, YDV,1)
         WOO(IACT+NCOR,ICOR)=-TMP-TMP
         POO(ICOR,II)=TMP/(EIG(II)-EIG(ICOR))
         POO(II,ICOR)=POO(ICOR,II)
      ENDDO


      END !***********************************************************









      SUBROUTINE RIMP2_XLOV_AA                                     &
                (WOO,POO, XLOV, BDVV,BDVC, YDV,EIG,                        &
                 IACT,NAUXBASD,NCOR,NACT,NVIR,NBF)
!
      implicit double precision(a-h,o-z)

!    OUTPUT
      double precision :: XLOV(NCOR+NACT,NVIR)

      double precision :: WOO(NCOR+NACT,NCOR+NACT)
      double precision :: POO(NCOR+NACT,NCOR+NACT)

!    INPUT
      double precision :: BDVV(NAUXBASD,NVIR,NVIR)
      double precision :: BDVC(NAUXBASD,NVIR,NCOR)
      double precision :: YDV(NAUXBASD,NVIR)
      double precision :: EIG(NCOR+NACT)



! !    PARAMETERS
!       NOCC=NCOR+NACT
!       NXO=NAUXBASD*NOCC
!       NXV=NAUXBASD*NVIR

      II=IACT+NCOR
      DO IA=1,NVIR
      ! DO IA=LddiVrtStart,LddiVrtEnd
         TMP=DDOT(NAUXBASD*NVIR, BDVV(1,1,IA),1, YDV,1)
         XLOV(II,IA)=XLOV(II,IA) -TMP-TMP
      ENDDO


! !    WOO(NACT,NCOR) = BI*YDV;
! !    POO(NCOR,NACT) = BI*YDV/(ei - eI);
! !    POO(NACT,NCOR) = POO(NCOR,NACT);
!       DO ICOR=1,NCOR
!          TMP=DDOT(NAUXBASD*NVIR, BDVC(1,1,ICOR),1, YDV,1)
!          WOO(IACT+NCOR,ICOR)=-TMP-TMP
!          POO(ICOR,II)=TMP/(EIG(II)-EIG(ICOR))
!          POO(II,ICOR)=POO(ICOR,II)
!       ENDDO


      END !***********************************************************





      SUBROUTINE RIMP2OMP_WOO_POO_XLOV(WOO_T,POO_T,XLOV, &
                 EIG,YDV_T, NAUXBASD,NCOR,NACT,NOCC,NVIR,IACT)

      use omp_lib

      USE Rimp2_Shared_Data,only: D_B
!
      implicit double precision(a-h,o-z)

      double precision :: WOO_T(NOCC,NOCC)
      double precision :: POO_T(NOCC,NOCC)
      double precision :: XLOV(NOCC,NVIR)

      double precision :: EIG(NOCC)
      double precision :: YDV_T(NAUXBASD,NVIR)

      double precision,allocatable :: BDV(:,:)

      II=IACT+NCOR
      NDV=NAUXBASD*NVIR
      NDO=NAUXBASD*NOCC

      Nthreads=omp_get_max_threads()


!    XLOV
!$OMP PARALLEL NUM_THREADS(Nthreads)                                    &
!$omp default(none)                                                     &
!$omp private(IV,KX,IA,IR,TMP)                                          &
!$omp shared(II,D_B,NDV,NDO,NAUXBASD,NCOR,NVIR,NOCC)                    &
!$omp private(BDV)                                                      &
!$omp shared(XLOV,WOO_T,POO_T,YDV_T,EIG)

      ALLOCATE(BDV(NAUXBASD,NVIR))

!$omp do schedule(guided)
!    XLOV
      DO IA=1,NVIR

         CALL DDI_GET(D_B, 1+NDO,NDV+NDO, IA+NOCC,IA+NOCC, BDV)

         ! TMP=0.0D00
         ! DO IV=1,NVIR
         !    DO KX=1,NAUXBASD
         !       TMP=TMP+BDV(KX,IV)*YDV_T(KX,IV)
         !    ENDDO
         ! ENDDO

         TMP = DDOT(NVIR*NAUXBASD, BDV,1, YDV_T,1)

         XLOV(II,IA)=XLOV(II,IA) -TMP -TMP
      ENDDO
!$omp end do nowait

!$omp do schedule(guided)
!    WOO,POO
      DO IR=1,NCOR

         CALL DDI_GET(D_B, 1+NDO,NDV+NDO, IR,IR, BDV)

         TMP=DDOT(NAUXBASD*NVIR,BDV,1,YDV_T,1)
         WOO_T(II,IR)=-TMP-TMP
         POO_T(IR,II)=POO_T(IR,II)+TMP/(EIG(II)-EIG(IR))
      ENDDO
!$omp end do



!$omp do schedule(guided)
!    WOO,POO
      ! DO IR=1,NCOR
      DO IR=NCOR+1,NOCC

         CALL DDI_GET(D_B, 1+NDO,NDV+NDO, IR,IR, BDV)

         TMP=DDOT(NAUXBASD*NVIR,BDV,1,YDV_T,1)
         WOO_T(II,IR)=-TMP-TMP
         ! POO_T(IR,II)=POO_T(IR,II)+TMP/(EIG(II)-EIG(IR))
      ENDDO
!$omp end do




      DEALLOCATE(BDV)

!$OMP END PARALLEL



!       ALLOCATE(BDV(NAUXBASD,NVIR))
! !    XLOV
!       DO IA=1,NVIR
!          CALL DDI_GET(D_B, 1+NDO,NDV+NDO, IA+NOCC,IA+NOCC, BDV)
!          TMP=DDOT(NAUXBASD*NVIR, BDV,1, YDV_T,1)
!          XLOV(II,IA)=XLOV(II,IA) -TMP-TMP
!       ENDDO
! !    WOO,POO
!       DO IR=1,NCOR
!          CALL DDI_GET(D_B, 1+NDO,NDV+NDO, IR,IR, BDV)
!          TMP=DDOT(NAUXBASD*NVIR,BDV,1,YDV_T,1)
!          WOO_T(II,IR)=-TMP-TMP
!          POO_T(IR,II)=POO_T(IR,II)+TMP/(EIG(II)-EIG(IR))
!       ENDDO
!       DEALLOCATE(BDV)

      END





      SUBROUTINE RIMP2_GNX_WOV_WVV                                      &
                (gDX, E2, WOV,WVV, GNX,                                 &
                 BDO,BDV,YDV,VDX, VEC, GXV,                             &
                 NAUXBASD,NAUXBAS,NCOR,NACT,NVIR,NBF,IACT)
      USE omp_lib
      USE Rimp2_Shared_Data,only: D_GNAX

      implicit double precision(a-h,o-z)


!!!!! OUTPUT
      double precision :: gDX(NAUXBASD,NAUXBAS)
      double precision :: E2
      double precision :: WOV(NCOR+NACT,NVIR)
      double precision :: WVV(NVIR,NVIR)


!!!!! DDI OUTPUT
      double precision :: GNX(NBF,NAUXBAS)


!!!!! INPUT
      double precision :: BDO(NAUXBASD,NCOR+NACT)
      double precision :: BDV(NAUXBASD,NVIR)
      double precision :: YDV(NAUXBASD,NVIR)
      double precision :: VDX(NAUXBASD,NAUXBAS)
      double precision :: VEC(NBF,NBF)


!!!!! LOCAL DATA
      double precision GXV(NAUXBAS,NVIR)


!!!!! PARAMETERS
      NOCC=NCOR+NACT


!**** CORRELATION ENERGY
      E2=E2+DDOT(NAUXBASD*NVIR, BDV,1, YDV,1)


      Nthreads=omp_get_max_threads()

!$OMP PARALLEL NUM_THREADS(Nthreads)                                    &
!$omp shared(NVIR,NOCC,NAUXBASD,BDO,YDV,WOV)                            &
!$omp private(IV)                                                       &
!$omp shared(VDX,GXV,BDV,WVV)

!$omp do schedule(guided)
      DO IV=1,NVIR
         CALL DGEMM('T','N', NOCC,1,NAUXBASD,                           &
                    -4.0D00, BDO,NAUXBASD,                              &
                             YDV(1,IV),NAUXBASD,                        &
                     1.0D00, WOV(1,IV),NOCC)
         CALL DGEMM('T','N', NVIR,1,NAUXBASD,                           &
                    -2.0D00, BDV,NAUXBASD,                              &
                             YDV(1,IV),NAUXBASD,                        &
                     1.0D00, WVV(1,IV),NVIR)
         CALL DGEMM('T','N',NAUXBAS,1,NAUXBASD,                         &
                     1.0D00, VDX,NAUXBASD,                              &
                             YDV(1,IV),NAUXBASD,                        &
                     0.0D00, GXV(1,IV),NAUXBAS)
      ENDDO
!$omp end do

!!!!! GNX = CNV x transpose(GXV)
!$omp do schedule(guided)
      DO IX=1,NAUXBAS
         CALL DGEMM('N','T', NBF,1, NVIR,                               &
                     1.0D00, VEC(1,1+NOCC),NBF,                         &
                             GXV(IX,1),NAUXBAS,                         &
                     0.0D00, GNX(1,IX),NBF)
      ENDDO
!$omp end do

!$OMP END PARALLEL

      CALL DGEMM('N','T', NAUXBASD,NAUXBAS,NVIR,                        &
                  1.0D00, BDV,NAUXBASD,                                 &
                          GXV,NAUXBAS,                                  &
                  1.0D00, gDX,NAUXBASD)

      END !***********************************************************




!>*module rimp2grd   *deck RIMP2_ACT_ACT
!>
!>     @brief   form density matrices
!>     @author  buu
!>     @date    sept 8, 2018
!>
!>     @detail  loop over active-active occ MOs
!>              to (partially) form density matrices
!>

      SUBROUTINE RIMP2_ACT_ACT                                          &
                (gDX, POO,PVV, WOO,WOV,WVV, XLOV,VDX,EIG,VEC,           &
                 NANGM,NAUXBASD,NAUXBAS, NCOR,NACT,NOCC,NVIR,NORB,NBF)

      USE omp_lib
      USE Rimp2_Shared_Data,only: D_B,D_GNAX, MXRT, PVV_OMP,POO_OMP,    &
                                  XLOV_OMP,WOO_OMP, WVV_OMP,YDV_OMP,TVV,&
                                  BUG

      implicit double precision(a-h,o-z)

      COMMON /IOFILE/ IR,IW,IP,IJK,IPK,IDAF,NAV,IODA(950)
      COMMON /ENRGYS/ ENUCR,EELCT,ETOT,SZ,SZZ,ECORE,ESCF,EERD,E1,E2,    &
                      VEN,VEE,EPOT,EKIN,ESTATE(MXRT),STATN,EDFT(3),EDISP
      COMMON /ENRGMP/ EMP2,EMP3,EMP4,EMP2A
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /MACHIN/ NWDVAR,MAXFM,MAXSM,LIMFM,LIMSM
      COMMON /FMOINF/ NFG,NLAYER,NATFMO,NBDFG,NAOTYP,NBODY,NSEGM
      COMMON /FMOOPT/ ESPSCA(9),RESPAP(2),RESPPC(2),RESDIM,RESTRI(4), &
                      RCORSD,RESPCT,CONVFG,CNVDMP,COROFF,RFLMO(4), &
                      ORSHFT,ORSHFT2,CNVAFO,ASCREEN(4),IXESP,MXITFG, &
                      NGUESS,NBSSE,MODORB,MODPAR,IRSTSTP,IRSTLAY,NPRFMO, &
                      NFMOPAL,MODPRP,MAXL1C,IPIEDA,MODGRD,MODESP,IVMUL, &
                      MODLMO,NOPDEN,MOFOCK,MODFD,MODFMM,NCENTM,NDUALB, &
                      ngab,modpan
      integer :: NFG,NLAYER,NATFMO,NBDFG,NAOTYP,NBODY,IXESP,MXITFG, &
                 NGUESS,NBSSE,MODORB,MODPAR,IRSTSTP,IRSTLAY,NPRFMO, &
                 NFMOPAL,MODPRP,MAXL1C,IPIEDA,MODGRD,MODESP,IVMUL,MODLMO &
             ,NOPDEN,MOFOCK,MODFD,MODFMM,NCENTM,NDUALB,ngab,modpan,NSEGM
      double precision :: ESPSCA,RESPAP,RESPPC,RESDIM,RESTRI,RCORSD, &
         RESPCT,CONVFG,CNVDMP,COROFF,RFLMO,ORSHFT,ORSHFT2,CNVAFO,ASCREEN



!!!!! OUTPUT
      double precision :: gDX(NAUXBASD,NAUXBAS)
      double precision :: POO(NOCC,NOCC)
      double precision :: PVV(NVIR,NVIR)
      double precision :: WOO(NOCC,NOCC)
      double precision :: WOV(NOCC,NVIR)
      double precision :: WVV(NVIR,NVIR)
      double precision :: XLOV(NOCC,NVIR)


!!!!! INPUT
      double precision :: VDX(NAUXBASD,NAUXBAS)
      double precision :: EIG(NBF)
      double precision :: VEC(NBF,NBF)


!!!!! LOCAL
      double precision,allocatable ::                                   &
                       YDV(:,:),GXV(:,:),GNX(:,:),                      &
                       BDO(:,:),BDOA(:,:,:),BDVA(:,:,:),BDVC(:,:,:),BDVV(:,:,:),    &
                       BDVA_T(:,:,:), eij(:,:), eab(:,:),QVVA(:,:,:)
      integer,allocatable :: ddiWSA(:,:),I_DLB_INDEX(:)
      logical :: GOPARR,DSKWRK,MASWRK,ULRange,LOMP,MASprt

      double precision, dimension(1) :: A_E2 ! for DDI_GSUMF with GCC 10


!!!!! TURN OFF DYNAMIC THREADS TO KEEP threadprivate DATA CONSISTENT
      CALL OMP_SET_DYNAMIC(.FALSE.)

      ! parameters
      NOCC=NCOR+NACT
      NORB=NCOR+NACT+NVIR

      NDO=NAUXBASD*NOCC
      NDV=NAUXBASD*NVIR
      NDN=NAUXBASD*NORB
      masprt=maswrk.and.(nfg.eq.0.or.iand(nprfmo,3).lt.2)


!!!!! SET BUFFER ZONE FOR 2 LAYER PARALLEL
      N_DLB_BUFF=0
      IF(NACT/NPROC .gt. 20) N_DLB_BUFF=1 ! DLB buffer

!!!!! OMP THREADS
      Nthreads_env=omp_get_max_threads()
      Nthreads=Nthreads_env
      IF(BUG) WRITE(*,'(A40,2I5)') "NUM_THREAD SUB RIMP2_ACT_ACT", Nthreads,Nthreads_env

      startw_other=omp_get_wtime()

!!!!! INIT
      gDX=0.0D00
      E2=0.0D00
      POO=0.0D00
      PVV=0.0D00
      WOO=0.0D00
      WOV=0.0D00
      WVV=0.0D00
      XLOV=0.0D00


!!!!! PAIRWISE MO ENERGY
      ALLOCATE(eab(NVIR,NVIR))
      FORALL(IA=1:NVIR,IB=1:NVIR)
         eab(IA,IB) = EIG(IA+NOCC) + EIG(IB+NOCC)
      ENDFORALL

      ALLOCATE(eij(NACT,NACT))
      FORALL(II=1:NACT,JJ=1:NACT)
         eij(II,JJ) = EIG(II+NCOR) + EIG(JJ+NCOR)
      ENDFORALL


!!!!! ddi local data chunk
      ALLOCATE(ddiWSA(0:NPROC-1,2))
      CALL rimp2_WSA(ddiWSA,NPROC,1,NACT,LOMP)
      LddiActStart=ddiWSA(ME,1)
      LddiActEnd=ddiWSA(ME,2)


!!!!! START threadprivate INITIALIZATION ////////////////////////////////

      ALLOCATE(BDVV(NAUXBASD,NVIR,NVIR))
      ALLOCATE(BDVA(NAUXBASD,NVIR,NACT))
      ALLOCATE(BDVC(NAUXBASD,NVIR,NCOR))
      ALLOCATE(BDOA(NAUXBASD,NOCC,LddiActStart:LddiActEnd))

!$OMP PARALLEL NUM_THREADS(Nthreads) &
!$omp shared(BDVV,BDVA,BDVC) &
!$omp shared(D_B,NDO,NDV, NCOR,NACT,NVIR,NAUXBASD) &
!$omp private(IR,II,IA,IB)

!!!!! first touch
      ! BDVV
!$omp do schedule(DYNAMIC)
      DO IB=1,NVIR
         BDVV(1:NAUXBASD,1:NVIR,IB) = 0.0D00
      ENDDO
!$omp end do nowait

      ! BDVA
!$omp do schedule(DYNAMIC)
      DO II=1,NACT
         BDVA(1:NAUXBASD,1:NVIR,II) = 0.0D00
      ENDDO
!$omp end do nowait

      ! BDVC
!$omp do schedule(DYNAMIC)
      DO IR=1,NCOR
         BDVC(1:NAUXBASD,1:NVIR,IR) = 0.0D00
      ENDDO
!$omp end do nowait

      ! BDOA
!$omp do schedule(DYNAMIC)
      DO II=LddiActStart,LddiActEnd
         BDOA(1:NAUXBASD,1:NOCC,II) = 0.0D00
      ENDDO
!$omp end do

!!!!! threadprivate ARRAYS
      ALLOCATE(PVV_OMP(NVIR,NVIR))
      ALLOCATE(POO_OMP(NOCC,NOCC))
      ALLOCATE(WOO_OMP(NOCC,NOCC))
      ALLOCATE(XLOV_OMP(NOCC,NVIR))
      ALLOCATE(YDV_OMP(NAUXBASD,NVIR))
      ALLOCATE(TVV(NVIR,NVIR))

      PVV_OMP=0.0D00
      POO_OMP=0.0D00
      WOO_OMP=0.0D00
      XLOV_OMP=0.0D00

!$OMP END PARALLEL


!!!!! GET MATRIX B
      CALL DDI_GET(D_B, 1+NDO,NDV+NDO, 1+NOCC,NVIR+NOCC, BDVV)
      CALL DDI_GET(D_B, 1+NDO,NDV+NDO, 1+NCOR,NACT+NCOR, BDVA)
      CALL DDI_GET(D_B, 1+NDO,NDV+NDO, 1,NCOR, BDVC)
      CALL DDI_GET(D_B, 1,NDO, LddiActStart+NCOR,LddiActEnd+NCOR, BDOA)

      endw_other=omp_get_wtime()
      IF(BUG) WRITE(*,'(A40,I5,F10.1)') "TIME ACT_ACT_OTHER: ",ME, endw_other-startw_other


!!!!! ON-NODE SHARED ARRAYS
      ALLOCATE(YDV(NAUXBASD,NVIR))
      ALLOCATE(GNX(NBF,NAUXBAS))
      ALLOCATE(GXV(NAUXBAS,NVIR))
      ALLOCATE(BDO(NAUXBASD,NOCC))
      ALLOCATE(QVVA(NVIR,NVIR,NACT))


      startw_static = omp_get_wtime()

      ! DO JACT=1,NACT
      DO JACT=LddiActStart,LddiActEnd !- N_DLB_BUFF
            JJ=JACT+NCOR

!!!!! YDV,PVV_OMP,POO_OMP
            YDV=0.0D00
            CALL RIMP2_IACT                                             &
                (YDV, BDVA,QVVA,eij,eab,                                &
                 JACT, NCOR,NACT,NVIR,NAUXBASD)

            CALL RIMP2_MEM_WOO_POO_XLOV                                 &
                (BDVV, BDVA,BDVC, EIG,YDV,                              &
                 NAUXBASD,NCOR,NACT,NOCC,NVIR,JACT)

!!!!! E2,WOV,WVV,GNX
            CALL RIMP2_GNX_WOV_WVV                                      &
                (gDX, E2, WOV,WVV, GNX,                                 &
                 BDOA(1,1,JACT),BDVA(1,1,JACT),YDV,VDX, VEC, GXV,       &
                 NAUXBASD,NAUXBAS,NCOR,NACT,NVIR,NBF,JACT)

!!!!! PUT GNX TO DDI MEM
            ISR = (JACT-1)*NBF + 1
            IER = ISR + NBF - 1
            CALL DDI_PUT(D_GNAX, ISR,IER, 1,NAUXBAS, GNX)

      ENDDO !JACT
      endw_static = omp_get_wtime()
      IF(BUG) WRITE(*,'(A40,I5,F10.1)') "TIME ACT_ACT_STATIC: ",ME, endw_static-startw_static




!       startw_dyn = omp_get_wtime()
!       ! DBL buffer array
!       ALLOCATE(I_DLB_INDEX(NACT))
!       CALL RIMP2_IDBLINDEX(I_DLB_INDEX,NACT,N_DLB_BUFF,ddiWSA)

!       MINE = -1
!       NEXT = -1

!       DO JACT=1,NACT
!             IF(I_DLB_INDEX(JACT).EQ.0) CYCLE

!             ! dynamic load balancing
!             MINE = MINE + 1
!             IF(MINE.GT.NEXT) CALL DDI_DLBNEXT(NEXT)
!             IF(NEXT.NE.MINE) CYCLE
!             ! rank with local data
!             IRANK=NPROC-1
!             DO IPROC=0,NPROC-2
!                IF(JACT.GE.ddiWSA(IPROC,2)-N_DLB_BUFF .AND. JACT.LT.ddiWSA(IPROC+1,1)) THEN
!                   IRANK=IPROC
!                   EXIT
!                ENDIF
!             ENDDO


!             JJ=JACT+NCOR
!             YDV=0.0D00


!             ! YDV,PVV_OMP,POO_OMP
!             CALL RIMP2_IACT                                             &
!                 (YDV, BDVA,QVVA,eij,eab,                                &
!                  JACT, NCOR,NACT,NVIR,NAUXBASD)
!             CALL RIMP2_MEM_WOO_POO_XLOV                                 &
!                 (BDVV, BDVA,BDVC, EIG,YDV,                              &
!                  NAUXBASD,NCOR,NACT,NOCC,NVIR,JACT)
!             ! CALL DDI_GET(D_B, 1,NDO, JJ,JJ, BDO)
! ! !          E2,WOV,WVV,GNX
! !             CALL RIMP2_GNX_WOV_WVV                                   &
! !                 (gDX, E2, WOV,WVV, GNX,                                 &
! !                  BDOA(1,1,JACT),BDVA(1,1,JACT),YDV,VDX, VEC, GXV,       &
! !                  NAUXBASD,NAUXBAS,NCOR,NACT,NVIR,NBF,JACT)
!             CALL DDI_GET(D_B, 1,NDO, JJ,JJ, BDO)
!             CALL RIMP2_GNX_WOV_WVV                                   &
!                 (gDX, E2, WOV,WVV, GNX,                                 &
!                  BDO,BDVA(1,1,JACT),YDV,VDX, VEC, GXV,       &
!                  NAUXBASD,NAUXBAS,NCOR,NACT,NVIR,NBF,JACT)


! !          PUT GNX TO DDI MEM
!             ISR = (JACT-1)*NBF + 1
!             IER = ISR + NBF - 1
!             CALL DDI_PUT(D_GNAX, ISR,IER, 1,NAUXBAS, GNX)
!       ENDDO !JACT
!       endw_dyn = omp_get_wtime()
!       WRITE(*,'(A40,I5,F10.1,I5)') "ACT_ACT/DYNAMIC: ", endw_dyn-startw_dyn, ME




      startw = omp_get_wtime()

!!!!! OMP REDUCTION :: PVV,POO,WOO
!$OMP PARALLEL NUM_THREADS(Nthreads) &
!$omp shared(PVV,POO,WOO)

!!!!! PVV
!$omp critical (PVV_critical)
      PVV=PVV+PVV_OMP
!$omp end critical (PVV_critical)

!!!!! POO
!$omp critical (POO_critical)
      POO=POO+POO_OMP
!$omp end critical (POO_critical)

!!!!! XLOV
!$omp critical (WOO_critical)
      XLOV=XLOV+XLOV_OMP
!$omp end critical (WOO_critical)

!!!!! WOO
!$omp critical (WOO_critical)
      WOO=WOO+WOO_OMP
!$omp end critical (WOO_critical)

      DEALLOCATE(WOO_OMP,XLOV_OMP,POO_OMP,PVV_OMP,TVV,YDV_OMP)

!$OMP END PARALLEL


      DO JJ=1+NCOR,NOCC
         DO II=1+NCOR,JJ-1
            POO(JJ,II)=POO(II,JJ)
         ENDDO
      ENDDO

!!!!! REDUCTION:: CORRELATION ENERGY
      A_E2(1) = E2
      CALL DDI_GSUMF(5999,A_E2,1)
      E2 = A_E2(1)

!!!!! RHF RI-MP2 ENERGY
      EMP2 = E2+ESCF

!!!!! PRINT OUT MP2 ENERGY
      IF(MASprt) WRITE(IW,100) ESCF,E2,EMP2
 100  FORMAT(/1X,'RHF RI MP2 ENERGY',/,                                 &
             12X,'   ESCF=',1X,F20.10/                                  &
             12X,'   E(2)=',1X,F20.10/                                  &
             12X,' E(MP2)=',1X,F20.10 ,/)

      endw = omp_get_wtime()
      IF(BUG) WRITE(*,'(A40,I5,F10.1)') "TIME ACT_ACT: OMP REDUCTION",ME, endw-startw



      END !***********************************************************






      SUBROUTINE RIMP2_IACT                                             &
                (YDV, BDVA,QVVA,eij,eab, JACT, NCOR,NACT,NVIR,NAUXBASD)

      USE omp_lib

      USE Rimp2_Shared_Data,only: PVV_OMP,YDV_OMP,POO_OMP,TVV

      implicit double precision(a-h,o-z)

      double precision :: YDV(NAUXBASD,NVIR)

      double precision :: BDVA(NAUXBASD,NVIR,NACT)
      double precision :: QVVA(NVIR,NVIR,NACT)

      double precision :: eab(NVIR,NVIR)
      double precision :: eij(NACT,NACT)




!!!!! TURN OFF DYNAMIC THREADS TO KEEP threadprivate DATA CONSISTENT
      CALL OMP_SET_DYNAMIC(.FALSE.)


!!!!! OMP THREADS
      Nthreads = omp_get_max_threads()


!$OMP PARALLEL NUM_THREADS(Nthreads)                                    &
!$omp default(none)                                                     &
!$omp shared(JACT, NCOR,NACT,NVIR,NAUXBASD)                             &
!$omp shared(YDV,BDVA,QVVA,eij,eab)                                     &
!$omp private(IACT,IA,IB)

!!!!! QVVA
!$omp do schedule(DYNAMIC)
      DO IACT=1, NACT
         CALL DGEMM('T','N', NVIR,NVIR,NAUXBASD,                        &
                     1.0D+00, BDVA(1,1,IACT),NAUXBASD,                  &
                              BDVA(1,1,JACT),NAUXBASD,                  &
                     0.0D+00, QVVA(1,1,IACT),NVIR)

         FORALL(IA=1:NVIR,IB=1:NVIR)
            QVVA(IA,IB,IACT)=QVVA(IA,IB,IACT)/(eij(IACT,JACT)-eab(IA,IB))
         ENDFORALL
      ENDDO
!$omp end do

      YDV_OMP=0.0D00

!$omp do schedule(DYNAMIC)
      DO IACT=1,NACT

!!!!! TVV, PVV_OMP, YDV_OMP
         CALL RIMP2_MEM_PVV_YDV                                         &
             (PVV_OMP,YDV_OMP,                                          &
              BDVA(1,1,IACT),QVVA(1,1,IACT),TVV,                        &
              NAUXBASD,NCOR,NACT,NVIR)

!!!!! POO_OMP
         CALL RIMP2_MEM_POO                                             &
             (POO_OMP, QVVA,TVV, JACT,IACT,                             &
              NAUXBASD,NCOR,NACT,NVIR)
      ENDDO !IACT
!$omp end do

!!!!! OMP REDUCTION :: YDV
!$omp critical (YDV_critical)
      YDV=YDV+YDV_OMP
!$omp end critical (YDV_critical)

!$OMP END PARALLEL

      END






      SUBROUTINE RIMP2_QVVA &
                (QVVA, JACT, BDVA,BJ,EIG,eab, NACT,NVIR,NAUXBASD)
      implicit none
      integer, intent(in) :: JACT, NACT, NVIR, NAUXBASD

      double precision :: QVVA(NVIR,NVIR,NACT)
      double precision :: BDVA(NAUXBASD,NVIR,NACT)
      double precision :: BJ(NAUXBASD,NVIR)
      double precision :: EIG(NACT+NVIR)
      double precision :: eab(NVIR,NVIR)
      double precision :: eij

      integer :: IACT, IA, IB

      DO IACT=1, NACT
         CALL DGEMM('T','N', NVIR,NVIR,NAUXBASD,        &
                     1.0D+00, BDVA(1,1,IACT),NAUXBASD,  &
                              BJ,NAUXBASD,              &
                     0.0D+00, QVVA(1,1,IACT),NVIR)

         eij=EIG(IACT)+EIG(JACT)
         FORALL(IA=1:NVIR,IB=1:NVIR)
            QVVA(IA,IB,IACT)=QVVA(IA,IB,IACT)/(eij-eab(IA,IB))
         ENDFORALL
      ENDDO

      END





!
      SUBROUTINE RIMP2_MEM_PVV_YDV                                      &
                (PVV,YDV, BI,QVV,TVV, NAUXBASD,NCOR,NACT,NVIR)

      implicit double precision(a-h,o-z)

      double precision :: PVV(NVIR,NVIR)
      double precision :: YDV(NAUXBASD,NVIR)

      double precision :: BI(NAUXBASD,NVIR)

      double precision :: QVV(NVIR,NVIR)
      double precision :: TVV(NVIR,NVIR)


!    TVV = 2*QVV - transpose(QVV)
      TVV = QVV + QVV
      TVV = TVV - transpose(QVV)

!    YDV <- BJ x transpose(TVV)
      CALL DGEMM('N','N', NAUXBASD,NVIR,NVIR,                           &
                  1.0D00, BI,NAUXBASD,                                  &
                          TVV,NVIR,                                     &
                  1.0D00, YDV,NAUXBASD)

!    PVV <- QVV x transpose(TVV)
      CALL DGEMM('N','T', NVIR,NVIR,NVIR,                               &
                 -1.0D00, QVV,NVIR,                                     &
                          TVV,NVIR,                                     &
                  1.0D00, PVV,NVIR)

      END !***********************************************************




!
      SUBROUTINE RIMP2_MEM_POO                                          &
                (POO, QVVA,TVV, JACT,IACT,                             &
                 NAUXBASD,NCOR,NACT,NVIR)

      implicit double precision(a-h,o-z)

!    OUTPUT
      double precision POO(NCOR+NACT,NCOR+NACT)

      double precision QVVA(NVIR,NVIR,NACT)
      double precision TVV(NVIR,NVIR)

!    PARAMETERS
      IOCC=IACT+NCOR

!    FORM WOO, POO
      DO IpACT=1, IACT ! NACT
         IpOCC=IpACT+NCOR
         POO(IpOCC,IOCC) = POO(IpOCC,IOCC) - DDOT(NVIR*NVIR, QVVA(1,1,IpACT),1, TVV,1)
      ENDDO


      END !***********************************************************





!
      SUBROUTINE RIMP2_MEM_POO_WOO                                          &
                (POO,WOO, BDVA,BJ, QVV,TVV,EIG,eab,                                  &
                 IACT,JACT,NAUXBASD,NCOR,NACT,NVIR,NORB,NBF)
!
      USE Rimp2_Shared_Data,only: D_B,D_VVAA
!
      implicit double precision(a-h,o-z)

!    OUTPUT
      double precision POO(NCOR+NACT,NCOR+NACT)
      double precision WOO(NCOR+NACT,NCOR+NACT)

      double precision BJ(NAUXBASD,NVIR)
      double precision BDVA(NAUXBASD,NVIR,NACT)

!    LOCAL DATA
      double precision QVV(NVIR,NVIR)
      double precision TVV(NVIR,NVIR)
      double precision EIG(NORB)
      double precision eab(NVIR,NVIR)

      ! double precision,allocatable :: QVVA(:,:,:)


!    PARAMETERS
      NOCC=NCOR+NACT
      IOCC=IACT+NCOR
      JOCC=JACT+NCOR

      POO(IOCC,IOCC) = POO(IOCC,IOCC) - DDOT(NVIR*NVIR, QVV,1, TVV,1)


!    FORM WOO, POO
      DO IpACT=1, IACT-1 ! NACT
         IpOCC=IpACT+NCOR

         CALL DGEMM('T','N', NVIR,NVIR,NAUXBASD,                        &
                     1.0D+00, BDVA(1,1,IpACT),NAUXBASD,                 &
                              BJ,NAUXBASD,                              &
                     0.0D+00, QVV,NVIR)


!       POO <- (-1) * QVV/DVV * TVV
         eij=EIG(IpOCC)+EIG(JOCC) ! EIG(JACT+NCOR)
         FORALL(IA=1:NVIR,IB=1:NVIR)
            QVV(IA,IB)=QVV(IA,IB)/(eij-eab(IA,IB))
         ENDFORALL
         POO(IpOCC,IOCC) = POO(IpOCC,IOCC) - DDOT(NVIR*NVIR, QVV,1, TVV,1)

      ENDDO


      END !***********************************************************





      SUBROUTINE RIMP2_MEM_WOO_POO_XLOV                                 &
                (BDVV,BDVA,BDVC, EIG,YDV,                               &
                 NAUXBASD,NCOR,NACT,NOCC,NVIR,IACT)

      USE omp_lib

      USE Rimp2_Shared_Data,only: D_B
      USE Rimp2_Shared_Data,only: XLOV_OMP,WOO_OMP,POO_OMP

      implicit double precision(a-h,o-z)


      double precision :: EIG(NOCC)
      double precision :: YDV(NAUXBASD,NVIR)

      double precision :: BDVV(NAUXBASD,NVIR,NVIR)
      double precision :: BDVA(NAUXBASD,NVIR,NACT)
      double precision :: BDVC(NAUXBASD,NVIR,NCOR)




      II=IACT+NCOR
      NDV=NAUXBASD*NVIR
      NDO=NAUXBASD*NOCC

      Nthreads=omp_get_max_threads()


!!!!! XLOV
!$OMP PARALLEL NUM_THREADS(Nthreads)                                    &
!$omp default(none)                                                     &
!$omp private(IA,IR,TMP)                                                &
!$omp shared(II,D_B,NAUXBASD,NCOR,NACT,NOCC,NVIR)                       &
!$omp shared(YDV,EIG) &
!$omp shared(BDVC,BDVA,BDVV)

!$omp do schedule(guided)
!!!!! XLOV_OMP
      DO IA=1,NVIR
         TMP = DDOT(NVIR*NAUXBASD, BDVV(1,1,IA),1, YDV,1)
         XLOV_OMP(II,IA)=XLOV_OMP(II,IA) -TMP -TMP
      ENDDO
!$omp end do nowait


!$omp do schedule(guided)
!!!!! WOO_OMP,POO_OMP
      DO IR=1,NCOR
         TMP=DDOT(NAUXBASD*NVIR,BDVC(1,1,IR),1,YDV,1)
         WOO_OMP(IR,II)=-TMP-TMP
         POO_OMP(IR,II)=POO_OMP(IR,II)+TMP/(EIG(II)-EIG(IR))
      ENDDO
!$omp end do nowait


!$omp do schedule(guided)
!!!!! WOO_OMP
      DO IIp=1+NCOR,NACT+NCOR
         TMP=DDOT(NAUXBASD*NVIR,BDVA(1,1,IIp-NCOR),1,YDV,1)
         WOO_OMP(II,IIp)=-TMP-TMP
      ENDDO
!$omp end do


!$OMP END PARALLEL


      END












      SUBROUTINE RIMP2_MEM_XLOV                                         &
                (XLOV, POO,PVV,VEC,EIG,                                 &
                 NAUXBASD,NCOR,NACT,NVIR,NORB,NBF)
!
      use omp_lib
      use Rimp2_Shared_Data,only : D_B,four,two,one,zero,BUG
!
      implicit double precision(a-h,o-z)
!
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      logical :: GOPARR,DSKWRK,MASWRK

!    OUTPUT
      double precision :: XLOV(NCOR+NACT,NVIR)

!    INPUT
      double precision :: POO(NCOR+NACT,NCOR+NACT)
      double precision :: PVV(NVIR,NVIR)
      double precision :: VEC(NBF,NORB)
      double precision :: EIG(NBF)

!    LOCAL DATA
      double precision,allocatable :: BDOO(:,:,:), BDOV(:,:,:),         &
                                      BDVO(:,:,:),BDVV(:,:,:),          &
                                      QOO(:,:),QVV(:,:),                &
                                      AOO(:,:),AVV(:,:),AVVT(:,:),      &
                                      AVVO(:,:,:),BDV(:,:),AVVV(:,:,:), &
                                      BA(:,:),BI(:,:),BDVV2(:,:,:)

      double precision,allocatable :: VVV(:,:,:),VVV_(:,:,:)

      double precision,save,allocatable :: XLOV_T(:,:)
!$omp threadprivate(XLOV_T)

      integer,allocatable :: ompISE(:,:),ompISEVirt(:,:)

      integer :: D_VVVO

      logical :: lomp

!    PARAMETERS
      NOCC=NCOR+NACT

      NDO=NAUXBASD*NOCC
      NDV=NAUXBASD*NVIR

!    DDI WORK DISTRIBUTING FOR NVIR LOOP
      ALLOCATE(ompISE(0:NPROC-1,2))
      CALL rimp2_WSA(ompISE,NPROC,1,NVIR,LOMP)
      LddiVirtStart=ompISE(ME,1)
      LddiVirtEnd=ompISE(ME,2)
      DEALLOCATE(ompISE)


      ALLOCATE(BDOO(NAUXBASD,NOCC,NOCC))
      ALLOCATE(BDOV(NAUXBASD,NOCC,LddiVirtStart:LddiVirtEnd))

      CALL DDI_GET(D_B, 1,NDO, 1,NOCC, BDOO)
      CALL DDI_GET(D_B, 1,NDO, LddiVirtStart+NOCC,LddiVirtEnd+NOCC, BDOV)


!    XLOV <- AOO*POO
      startw = omp_get_wtime()


!    OMP WORK DISTRIBUTING FOR LddiVirtStart:LddiVirtEnd LOOP
      Nthreads=omp_get_max_threads()
      Nchunk=LddiVirtEnd-LddiVirtStart+1
      IF(Nthreads.GT.Nchunk) Nthreads=Nchunk
      ALLOCATE(ompISE(0:Nthreads-1,2))
      CALL rimp2_WSA(ompISE,Nthreads,LddiVirtStart,LddiVirtEnd,LOMP)


!$OMP PARALLEL NUM_THREADS(Nthreads)                                    &
!$omp default(none)                                                     &
!$omp shared(ompISE)                                                    &
!$omp private(IthreadID,LompVirtStart,LompVirtEnd)                      &
!$omp private(IA,II)                                                    &
!$omp private(TMP,AOO)                                                  &
!$omp shared(NVIR,NOCC,NAUXBASD)                                        &
!$omp shared(XLOV,POO,BDOO,BDOV)

      ALLOCATE(AOO(NOCC,NOCC))

      IthreadID=omp_get_thread_num()
      LompVirtStart=ompISE(IthreadID,1)
      LompVirtEnd=ompISE(IthreadID,2)

      ! DO IA=1,NVIR
      ! DO IA=LddiVirtStart,LddiVirtEnd
      DO IA=LompVirtStart,LompVirtEnd
         DO II=1,NOCC
            CALL DGEMM('T','N', NOCC*NOCC,1,NAUXBASD,                   &
                        4.0D00, BDOO,NAUXBASD,                          &
                                 BDOV(1,II,IA),NAUXBASD,                &
                        0.0D00, AOO,NOCC*NOCC)
            CALL DGEMM('T','N', NOCC,NOCC,NAUXBASD,                     &
                       -2.0D00, BDOV(1,1,IA),NAUXBASD,                  &
                                 BDOO(1,1,II),NAUXBASD,                 &
                        1.0D00, AOO, NOCC)
            TMP = DDOT(NOCC*NOCC, AOO,1, POO,1)
            XLOV(II,IA)=XLOV(II,IA) - TMP
         ENDDO
      ENDDO

      DEALLOCATE(AOO)

!$OMP END PARALLEL

      DEALLOCATE(BDOO,BDOV)

      endw = omp_get_wtime()
      IF(BUG) WRITE(*,'(A40,I5,F10.1)') "XLOV: AOO_POO",ME,endw-startw





!    XLOV <- AVV*PVV
      startw = omp_get_wtime()

      ALLOCATE(BDVV(NAUXBASD,NVIR,LddiVirtStart:LddiVirtEnd))

      CALL DDI_GET(D_B, 1+NDO,NDV+NDO, LddiVirtStart+NOCC,LddiVirtEnd+NOCC, BDVV)

      Nthreads=omp_get_max_threads()


      ALLOCATE(BDV(NAUXBASD,NVIR))
      ALLOCATE(VVV(NVIR,NVIR,LddiVirtStart:LddiVirtEnd))

      DO II=1,NOCC

         CALL DDI_GET(D_B, 1+NDO,NDV+NDO, II,II, BDV)

!$OMP PARALLEL NUM_THREADS(Nthreads)                                    &
!$omp default(none)                                                     &
!$omp private(IA,IB,IC,TMP)                                             &
!$omp shared(II,LddiVirtStart,LddiVirtEnd)                              &
!$omp shared(NVIR,NAUXBASD)                                             &
!$omp shared(BDV,BDVV,VVV,PVV,XLOV)

!$omp do schedule(guided)
         DO IB=LddiVirtStart,LddiVirtEnd !DDI

            CALL DGEMM('T','N', NVIR,NVIR,NAUXBASD,                     &
                        1.0D00, BDV,NAUXBASD,                          &
                                 BDVV(1,1,IB),NAUXBASD,                 &
                        0.0D00, VVV(1,1,IB),NVIR)


            TMP=DDOT(NVIR*NVIR, VVV(1,1,IB),1, PVV,1)

            XLOV(II,IB)=XLOV(II,IB)-TMP-TMP

         ENDDO
!$omp end do


!$omp do schedule(guided)
         DO IA=1,NVIR
            TMP=0.0D00
            DO IB=LddiVirtStart,LddiVirtEnd !DDI
               DO IC=1,NVIR
                  TMP=TMP+VVV(IA,IC,IB)*PVV(IC,IB)
               ENDDO
            ENDDO
            XLOV(II,IA)=XLOV(II,IA)+TMP+TMP+TMP+TMP
         ENDDO
!$omp end do


!$OMP END PARALLEL

      ENDDO

      endw = omp_get_wtime()
      IF(BUG) WRITE(*,'(A40,I5,F10.1)') "XLOV: AVV_PVV",ME,endw-startw


      END !***********************************************************




      SUBROUTINE RIMP2SWAP(I,J)
         ITMP=I
         I=J
         J=ITMP
      END




      SUBROUTINE RIMP2_PRE_HONDO &
                (SFAC, ISH,JSH,KSH,LSH,NIJ,NIJ0,NKL,NKL0, &
                 INEW,JNEW,KNEW,LNEW, IIJJ,KKLL,X_INIJG,X_IIJKLG) !INIJG,IIJKLG)

      implicit double precision(a-h,o-z)

      integer :: X_IIJKLG(*) !,X_NIJG(*)
      integer :: X_INIJG(2,*)


      ! symmetry weight factor
      SFAC = 1.0D00
      IF(IIJJ.NE.KKLL) SFAC = 0.5D00
      IF(ISH.NE.JSH) SFAC = SFAC*0.5D00
      IF(KSH.NE.LSH) SFAC = SFAC*0.5D00
      IF(KSH.NE.ISH) SFAC = SFAC*2.0D00


      ! set pointers to the IJ& KL charge distribution
      CALL OEDRD(X_INIJG,NIJ,NIJ0,IIJJ)
      CALL OEDRD(X_INIJG,NKL,NKL0,KKLL)


      ! select differentiation center ??
      CALL JKDATM(INEW,JNEW,KNEW,LNEW)

      ! set indices for shell block ??
      CALL JKDSHL(INEW,JNEW,KNEW,LNEW)
      CALL JKDNDX(X_IIJKLG)

      END




      SUBROUTINE RIMP2_D42_INIT &
                (MP2,POPLE,HONDO, NANGM, &
                 CUTOFFDER,CUTOFF2,DABCUT, &
                 Q4,VTOL1,VTOL2,VTOLS,DTOL,RTOL)

      use mx_limits, only: mxgtot,mxsh
      implicit double precision(a-h,o-z)

      COMMON /NSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT), &
                      CF(MXGTOT),CG(MXGTOT),CH(MXGTOT),CI(MXGTOT), &
                      KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH), &
                      KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL !
      COMMON /INTOPT/ ISCHWZ,IECP,NECP,IEFLD
      COMMON /INTFIL/ NINTMX,NHEX,NTUPL,PACK2E,INTTYP,IGRDTYP
      COMMON /OUTPUT/ NPRINT,ITOL,ICUT,NORMF,NORMP,NOPK !
      COMMON /MP2PAR/ OSPT,CODEMP,SCSPT,TOL,METHOD,NWDMP2,MEMPRI,MPPROP,&
                      NACORE,NBCORE,NOA,NOB,NO,NBF,NOMIT,MOCPHF,MAXITC
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      logical :: GOPARR,DSKWRK,MASWRK
      ! Only IW is needed
      integer :: IR, IW, IP, IJK, IJKT, IDAF, NAV, IODA(950)
      COMMON /IOFILE/ IR, IW, IP, IJK, IJKT, IDAF, NAV, IODA

      logical :: MP2,POPLE,HONDO
      double precision,parameter :: RLN10=2.30258D+00

      CALL BASCHK(LMAX)
                    NANGM =  4
      IF(LMAX.EQ.2) NANGM =  6
      IF(LMAX.EQ.3) NANGM = 10
      IF(LMAX.EQ.4) NANGM = 15
      IF(LMAX.EQ.5) NANGM = 21
      IF(LMAX.EQ.6) NANGM = 28
      MAXG = NANGM**4


      ! find largeest shell type
      NBF3=(NBF*NBF+NBF)/2
      CALL BASCHK(MAXTYP)
      IF(MAXTYP.GT.5) THEN
         IF(MASWRK) WRITE(IW,*) 'RI-MP2 GRADIENT INTEGRALS PROGRAMMED FOR SPDFG ONLY'
         CALL ABRT
      END IF

      ! init
      MP2   = .TRUE.
      HONDO = .TRUE.
      POPLE = .TRUE.
      IF(IGRDTYP.EQ.2)          POPLE = .FALSE.
      IF(POPLE.AND.MAXTYP.LT.2) HONDO = .FALSE.

      ! Schwartz screening
      CUTOFFDER=1.0D-09
      IF(.NOT.POPLE) CUTOFFDER=CUTOFFDER/10.0D00
      CUTOFF2=CUTOFFDER*0.5D+00

      ! DABCUT
      ZBIG = 0.0D+00
      DO ISH=1,NSHELL
         I1=KSTART(ISH)
         I2=I1+KNG(ISH)-1
         DO IG=I1,I2
            IF(EX(IG).GT.ZBIG) ZBIG = EX(IG)
         ENDDO
      ENDDO
      DABCUT=1.0D-11
      IF(ZBIG.GT.1.0D+06) DABCUT = DABCUT/10.0D00
      IF(ZBIG.GT.1.0D+07) DABCUT = DABCUT/10.0D00

      ! Pople cut off
      Q4   = 1.0D00
      VTOL1= 1.0D-12
      VTOL2= 1.0D-12
      VTOLS= 1.0D-20
      DTOL = 10.0D00**(-ITOL)
      RTOL = RLN10*ITOL

      END





      SUBROUTINE RIMP2_D42DEN_INIT(XINT,PHF,PMP2,PTOT, NSH2,NBF3)

      COMMON /INTOPT/ ISCHWZ,IECP,NECP,IEFLD
      COMMON /IOFILE/ IR,IW,IP,IJK,IPK,IDAF,NAV,IODA(950)

      double precision :: XINT(NSH2)
      double precision :: PHF(NBF3)
      double precision :: PMP2(NBF3)
      double precision :: PTOT(NBF3)

      ! read Schwartz screening
      IF(ISCHWZ.EQ.1) THEN
         CALL DAREAD(IDAF,IODA,XINT,NSH2,54,0)
      ELSE
         CALL DCOPY(NSH2,1.0D01,0,XINT,1)
      END IF

      CALL DAREAD(IDAF,IODA,PHF,NBF3,308,0)
      CALL DAREAD(IDAF,IODA,PMP2,NBF3,307,0)
      PTOT = PHF + PMP2

      END





!>*module rimp2grd   *deck RIMP2_JKDER
!>
!>     @brief   4-2ERI gradient contribution
!>     @author  buu
!>     @date    sept 8, 2018
!>
!>     @detail  4-2ERI gradient contribution
!>

      SUBROUTINE RIMP2_JKDER

      USE omp_lib

      USE RIMP2_TMP,only:  NX_INIJG,NX_IIJKLG        ,&
         X_IXCH  ,X_INIJG ,X_ICHRG ,X_IIJKLG,X_IDAB  ,X_IGINT ,X_IFINT ,&
         X_ISINT ,X_IGIJKL,X_IGNKL ,X_IGNM  ,X_IXY   ,X_IXZ   ,X_IYZ   ,&
         X_IX    ,X_IY    ,X_IZ    ,X_ISJ   ,X_ISK   ,X_ISL   ,X_IB00  ,&
         X_IB01  ,X_IB10  ,X_IC00  ,X_ID00  ,X_IF00  ,X_IDIJ  ,X_IDKL  ,&
         X_IDIJSI,X_IDIJSJ,X_IDKLSK,X_IDKLSL,X_IABV  ,X_ICV   ,X_IRW   ,&
         X_IAAI  ,X_IAAJ  ,X_IBBK  ,X_IBBL  ,X_IFI   ,X_IFJ   ,X_IFK   ,&
         X_IFL   ,X_ISII  ,X_ISJJ  ,X_ISKK  ,X_ISLL  ,X_ISIJ  ,X_ISIK  ,&
         X_ISIL  ,X_ISJK  ,X_ISJL  ,X_ISKL

      use mx_limits,only:mxfrg,mxfgpt,mxdfg,mxdppt,mxgtot,mxsh,mxatm,&
      mxao,mxrt,mxnoro,mxgrid,mxfrz
      USE camdft, ONLY: CAMFLAG
      USE lrcdft, ONLY: LCFLAG, EMU, EMU2, LRFILE

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION METHMC,MD
!
      logical :: SKIPI,SKIPJ,SKIPK,SKIPL,NOCORE,NOOPEN
      logical :: HFSCF,UHFTYP,ROGVB,PACK2E,POPLE,HONDO
      logical :: GOPARR,DSKWRK,MASWRK,SLB,DLB
      logical :: SOME,OUT,DBG,MP2,MC,PER1,PER2,PER3,UMP2
      logical :: CANONC,FCORE,FORS,EKT,LINSER
      logical :: TDDFT,LRINT,SG1,PROJGRAD,LINEAR
      logical :: LCIS,LSFCIS,LSFDFT
!
!$    LOGICAL isdft,isgfdft,ischk,isgeo,ismd,isfmo
!
      double precision,parameter ::                                     &
         RLN10=2.30258D+00,TEN=10.0D+00,ONE=1.0D+00,TENM9=1.0D-09,      &
         TENM11=1.0D-11,TENM20=1.0D-20,PT5=0.5D+00,TENM12=1.0D-12



      COMMON /DETWFN/ WSTATE(MXRT),SPINS(MXRT),CRIT,PRTTOL,SDET,SZDET,  &
                      GRPDET,STSYM,GLIST,DWPARM,                        &
                      NFLGDM(MXRT),IWTS(MXRT),NCORSV,NCOR,NACT,NORBDT,  &
                      NADET,NBDET,KDET,KSTDET,IROOT,IPURES,MAXW1,NITER, &
                      MAXP,NCIDET,IGPDET,KSTSYM,NFTGCI,IDWEIGH,         &
                      fstate(mxrt),ifts(mxrt)
      COMMON /DFGRID/ DFTTHR,DFTGTHR,SWOFF,SW0,BSLRD(137),NDFTFG,       &
                      NRAD,NTHE,NPHI,NRAD0,NTHE0,NPHI0,                 &
                      NANGPT(MXGRID),NANGPT0(MXGRID),SG1,JANS
!     COMMON /DLT   / LAT,LBT,LCT,LDT
      COMMON /DSHLT / RTOL,DTOL,VTOL1,VTOL2,VTOLS,OUT,DBG
      COMMON /FFPARM/ NFFAT,NBOND,NANGL,NDIHR,NDIHB,NCMAP,NWAGG,        &
                      N1213J,N14J,NLKQMM,IDOCHG,IDOPOL,IDOLJ,IDOCMAP
      COMMON /FMCOM / X(1)
      COMMON /FMOINF/ NFG,NLAYER,NATFMO,NBDFG,NAOTYP,NBODY,NSEGM
      COMMON /FRGINF/ NMPTS(MXFRG),NMTTPT,IEFC,IEFD,IEFQ,IEFO,          &
                      NPPTS(MXFRG),NPTTPT,IEFP,                         &
                      NRPTS(MXFRG),NRTTPT,IREP,ICHGP,NFRG,              &
                      NDPPTS(MXDPPT),NDPTTPT,IEFDP,LSTMPTS(MXFRG),      &
                      NBSFN(MXFRG),NMXMO(MXFRG)

      COMMON /GUGWFN/ NFZC,NMCC,NDOC,NAOS,NBOS,NALP,NVAL,NNEXT,NFZV,    &
                      IFORS,IEXCIT,ICICI,NOIRR
      COMMON /IJPAIR/ IA(MXAO)
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,                    &
                      ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /INTFIL/ NINTMX,NHEX,NTUPL,PACK2E,INTTYP,IGRDTYP
      COMMON /INTOPT/ ISCHWZ,IECP,NECP,IEFLD
      COMMON /IOFILE/ IR,IW,IP,IJK,IPK,IDAF,NAV,IODA(950)
      COMMON /MASSES/ ZMASS(MXATM)
      COMMON /MCINP / METHMC,CISTEP,FINALCI,ACURCY,ENGTOL,DAMP,         &
                      MICIT,NWRDMC,NORBMC,NOROT(2,MXNORO),MOFRZ(MXFRZ), &
                      NPFLG(10),NOFO,MCFMO,IDIABAT,                     &
                      CANONC,FCORE,FORS,EKT,LINSER
      COMMON /MP2DM2/ NOC1,NOC2,ISTEP,NOC3,NOC4
      COMMON /MP2PAR/ OSPT,CODEMP,SCSPT,TOL,METHOD,NWDMP2,MEMPRI,MPPROP,&
                      NACORE,NBCORE,NOA,NOB,NORB,NBF,NOMIT,MOCPHF,MAXITC
      COMMON /NSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),      &
                      CF(MXGTOT),CG(MXGTOT),CH(MXGTOT),CI(MXGTOT),      &
                      KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),   &
                      KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
      COMMON /OUTPUT/ NPRINT,ITOL,ICUT,NORMF,NORMP,NOPK
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /PAULMO/ CENTNM(MXFGPT),CENTCD(3,MXFGPT),NORBEF(MXFRG),    &
                      NPBF(MXFRG),NTMO
      COMMON /PCMPAR/ IPCM,NFT26,NFT27,IRPPCM,IEF,IP_F,NFMOPCM,IHET
      COMMON /RESTAR/ TIMLIM,IREST,NREC,INTLOC,IST,JST,KST,LST
      COMMON /RUNOPT/ RUNTYP,EXETYP,NEVALS,NGLEVL,NHLEVL
      COMMON /SHLBAS/ MAXTYP,MAXNUM
      COMMON /SIMDAT/ NACC,NREJ,IGOMIN,NRPA,IBWM,NACCT,NREJT,NRPAT,     &
                      NPRTGO,IDPUNC,IGOFLG
      COMMON /SYMTRY/ MAPSHL(MXSH,48),MAPCTR(MXATM,48),                 &
                      T(432),INVT(48),NT

      COMMON /WFNOPT/ SCFTYP,VBTYP,DFTYPE,TDDFTYP,CITYP,CCTYP,          &
                      MPLEVL,MPCTYP
      COMMON /ZMAT  / NZMAT,NZVAR,NVAR,NSYMC,LINEAR

      COMMON/DERMEM/IWFN,IXCH,INIJG,IGINT,IFINT,ISINT,IIJKLG,           &
       IDAB,ICHRG,IXY,IXZ,IYZ,IX,IY,IZ,ISJ,ISK,ISL,IGIJKL,IGNKL,IGNM,   &
       IDIJ,IDKL,IB00,IB01,IB10,IC00,ID00,IF00,                         &
       IAAI,IAAJ,IBBK,IBBL,IFI,IFJ,IFK,IFL,                             &
       ISII,ISJJ,ISKK,ISLL,ISIJ,ISIK,ISIL,ISJK,ISJL,ISKL,               &
       IDIJSI,IDIJSJ,IDKLSK,IDKLSL,IABV,ICV,IRW

      COMMON /DERMEMALLOCATE/ &
         N_IWFN  ,N_IXCH  ,N_IDAB  ,N_INIJG ,N_IGINT ,N_IFINT ,N_ISINT ,&
         N_IIJKLG,N_ICHRG ,N_IXY   ,N_IXZ   ,N_IYZ   ,N_IX    ,N_IY    ,&
         N_IZ    ,N_ISJ   ,N_ISK   ,N_ISL   ,N_IGIJKL,N_IGNKL ,N_IGNM  ,&
         N_IDIJ  ,N_IDKL  ,N_IB00  ,N_IB01  ,N_IB10  ,N_IC00  ,N_ID00  ,&
         N_IF00  ,N_IDIJSI,N_IDIJSJ,N_IDKLSK,N_IDKLSL,N_IABV  ,N_ICV   ,&
         N_IRW   ,N_IAAI  ,N_IAAJ  ,N_IBBK  ,N_IBBL  ,N_IFI   ,N_IFJ   ,&
         N_IFK   ,N_IFL   ,N_ISII  ,N_ISJJ  ,N_ISKK  ,N_ISLL  ,N_ISIJ  ,&
         N_ISIK  ,N_ISIL  ,N_ISJK  ,N_ISJL  ,N_ISKL

      COMMON /DERSKP/ IIAT,JJAT,KKAT,LLAT,SKIPI,SKIPJ,SKIPK,SKIPL
!$omp threadprivate(/DERSKP/)
      COMMON /DFTPAR/ DFTTYP(20),EXENA,EXENB,EXENC,IDFT34,NAUXFUN,NAUXSHL
!$omp threadprivate(/DFTPAR/)
      COMMON /DSHLNO/ LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,              &
                      MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL,          &
                      NIJ,IJD,KLD,IJ,KL
!$omp threadprivate(/DSHLNO/)
      COMMON /GRAD  / DE(3,MXATM)
!$omp threadprivate(/GRAD/)
      COMMON /NLRCF / LRINT
!$omp threadprivate(/NLRCF /)
      COMMON /TMVALS/ TI,TX,TIM
!$omp threadprivate(/TMVALS/)
!$omp threadprivate(/DERMEM/)
      COMMON/DERPAR/IDER,JDER,KDER,LDER,NDER
!$omp threadprivate(/DERPAR/)

      COMMON/INDD80/IMAX,JMAX,KKKMAX,LMAX
!$omp threadprivate(/INDD80/)

      integer,dimension(48),save :: M0,M1,M2,M3
!$omp threadprivate(M0,M1,M2,M3)




      integer :: LENSHL(5)
      DATA LENSHL/1,4,10,20,35/
      DATA CHECK,GRD2,DEBUG/8HCHECK   ,8HGRD2    ,8HDEBUG   /
      DATA RHF,UHF,ROHF,GVB                                             &
          /8HRHF     ,8HUHF     ,8HROHF    ,8HGVB     /
      DATA RMC/8HMCSCF    /
      DATA GUGA,ALDET,GENCI,ORMAS,GMCCI                                 &
             /8HGUGA    ,8HALDET   ,8HGENCI   ,8HORMAS   ,8HGMCCI   /
      DATA CIS/8HCIS     /,RNONE/8HNONE    /
      DATA OPT,SADPT,MD/8HOPTIMIZE,8HSADPOINT,8HMD      /
      DATA SFCIS/8HSFCIS   /
      DATA EXCITE,SPNFLP/8HEXCITE  ,8HSPNFLP  /
      DOUBLE PRECISION MRSF
      DATA MRSF/8HMRSF    /

      logical :: I_AND_J,I_AND_K,I_AND_JK,I_AND_L,I_AND_JKL


      double precision,allocatable :: &
         PHF(:),PMP2(:),PTOT(:),DE_DDI(:,:),DE_OMP(:,:)



      ! some parameters
      MP2=.TRUE.
      L1 = NUM !number of basis functions
      L2 = (NUM*NUM+NUM)/2
      L3 = NUM*NUM
      NSH2=(NSHELL*NSHELL+NSHELL)/2

      ! number of threads
      Nthreads = omp_get_max_threads()

      ! turn off dynamic threads
      CALL OMP_SET_DYNAMIC(.FALSE.)


      ! initilize Rys/Pople var
      HONDO = .TRUE.
      POPLE = IGRDTYP.NE.2

      ! integral cutoffs
      CALL RIMP2_INTCUTOFF &
          (CUTOFF,CUTOFF2,DABCUT,                                       &
           VTOL1,VTOL2,VTOLS,DTOL,RTOL, ITOL,POPLE,KSTART,EX,KNG,NSHELL)


      ! setup COMMON /IJPAIR/ IA(MXAO)
      FORALL(I=1:NUM)
         IA(I) = (I*I-I)/2
      END FORALL

      ! find largest shell
      CALL BASCHK(MAXTYP)
      IF(MAXTYP.GT.5)  THEN
         IF(MASWRK) WRITE(IW,*) 'GRADIENT INTEGRALS FOR SPDFG ONLY'
         CALL ABRT
      END IF

      ! set Rys to false if only SP integrals
      IF (POPLE.AND.MAXTYP.LT.2) HONDO = .FALSE.

      ! set memory pointers LENGTH=1; LOADFM=1
      CALL RIMP2_JKDMEM(1,1,IADDR,1,MINXYZ,MAXXYZ,MINVEC,POPLE,MP2)

      ! densities
      ALLOCATE(PHF(L2))    !scf density
      ALLOCATE(PMP2(L2))   !mp2 density
      ALLOCATE(PTOT(L2))   !scf+mp2 density
      !exchange screening integral
      ALLOCATE(X_IXCH(N_IXCH))
      CALL RIMP2_READSTUFF(PHF,PMP2,PTOT,X_IXCH, N_IXCH,L2)

      ! 1-electron charge distribution
      ALLOCATE(NX_INIJG(N_INIJG))
      ALLOCATE(X_ICHRG(N_ICHRG))
      IF(HONDO) CALL OEDHND(NX_INIJG,X_ICHRG)



      ! allocate (threadprivate) arrays for Rys
      CALL RIMP2_D42ARRAYS


      ! copy 1-electron gradient to DE_DDI
      ALLOCATE(DE_DDI(3,NAT))
      DE_DDI(1:3,1:NAT) = (1.0D00/NPROC) * DE(1:3,1:NAT)

      ! zero out gradient vector
!$OMP PARALLEL NUM_THREADS(Nthreads) shared(NAT)
      DE(1:3,1:NAT) = 0.0D00
!$OMP END PARALLEL

      ! initilize the number of
      ! zero integral derivative blocks
      IISKIP1= 0
      IISKIP2= 0
      IDID = 0


      ! square DTOL for Rys package
      DTOL = DTOL*DTOL
      NC=1
      LDF=1

      NEXT = -1
      MINE = -1

      ! I-SHELL
      DO II = 1,NSHELL
         IATM= KATOM(II)

         ! CALL RIMP2_M0(M0,II)

         ! go parallel
         MINE = MINE + 1
         IF (MINE.GT.NEXT) CALL DDI_DLBNEXT(NEXT)
         IF (NEXT.NE.MINE) CYCLE

      ! ! turn off dynamic threads
      ! CALL OMP_SET_DYNAMIC(.FALSE.)


!$OMP PARALLEL NUM_THREADS(Nthreads)                                    &
!$omp default(none)                                                     &
!$omp shared(IA,II,IATM,HONDO)                                          &
!$omp shared(KATOM,PHF,PMP2,PTOT, NX_INIJG,X_IXCH,X_ICHRG)              &
!$omp shared(CUTOFF,CUTOFF2, DABCUT)                                    &
!$omp shared(L1,L2,MAXXYZ,MINVEC)                                       &
!$omp copyin(M0,NDER)                                                   &
!$omp private(JJ,KK,JATM,KATM, I_AND_J,I_AND_K,I_AND_JK)                &
!$omp private(IIJJ,NIJ0)                                                &
!!!!
!$omp shared(IGRDTYP,KTYPE,NBF)                                         &
!$omp private(MAXLL,LATM,I_AND_L,I_AND_JKL,Q4)                          &
!$omp private(KLKL,INEW,JNEW,KNEW,LNEW)                                 &
!$omp private(IAT,JAT,KAT,LAT)                                          &
!$omp private(POPLE,IJIJ,GMAX,DABMAX,JTYPE)                             &
!$omp private(ISH,JSH,KSH,LSH)                                          &
!$omp private(NKL,NKL0)                                                 &
!$omp private(FC,NC,DF,LDF,DDA)


!!!!! J-SHELL
!$omp do schedule(DYNAMIC)
         DO JJ=1,II
            JATM = KATOM(JJ)
            I_AND_J = IATM.EQ.JATM


!!!!! SET POINTERS TO CHARGE DISTRIBUTION
            IF(HONDO) THEN
               IIJJ=IA(MAX(II,JJ))+MIN(II,JJ)
               CALL OEDRD(NX_INIJG,NIJ,NIJ0,IIJJ)
               IF(NIJ.EQ.0) CYCLE
            END IF

!!!!! K-SHELL
            DO KK = 1,II
               KATM = KATOM(KK)
               I_AND_K = IATM.EQ.KATM
               I_AND_JK = I_AND_J.AND.I_AND_K

               CALL RIMP2_M3(M3,KK,II)

               MAXLL = KK
               IF (KK.EQ.II) MAXLL = JJ

               DO LL = 1,MAXLL

                  LATM = KATOM(LL)
                  I_AND_L = IATM.EQ.LATM
                  I_AND_JKL = I_AND_JK.AND.I_AND_L

!!!!! 1-CENTER INTEGRAL DERIVATIVE = 0.0D00
                  IF(I_AND_JKL) CYCLE

                  Q4=1.0D00

!!!!! CHECK IF POPLE/RYS PACKAGE IS USED
                  POPLE = .TRUE.
                  IF(IGRDTYP.EQ.2)   POPLE=.FALSE.
                  IF(KTYPE(II).GT.2) POPLE=.FALSE.
                  IF(KTYPE(JJ).GT.2) POPLE=.FALSE.
                  IF(KTYPE(KK).GT.2) POPLE=.FALSE.
                  IF(KTYPE(LL).GT.2) POPLE=.FALSE.

                  IJIJ=IA(MAX(II,JJ)) + MIN(II,JJ)
                  KLKL=IA(MAX(KK,LL)) + MIN(KK,LL)
                  GMAX = X_IXCH(IJIJ)*X_IXCH(KLKL)

!!!!! EXCHANGE INTEGRAL SCREENING
!!!!! JCC 12 (1991) 1058-1064

                  IF (GMAX.LT.CUTOFF) CYCLE

                  ISH=II
                  JSH=JJ
                  KSH=KK
                  LSH=LL

!!!!! POPLE PACKAGE
                  IF (POPLE) THEN
                     CALL RIMP2_PRE_POPLE                               &
                         (INEW,JNEW,KNEW,LNEW,                          &
                          IAT,JAT,KAT,LAT, JTYPE,                       &
                          ISH,JSH,KSH,LSH,                              &
                          Q4,KTYPE,KATOM)

!!!!! POPLE DENSITY STORED
!!!!! IN COMMON /GSPG80/ E(4,4,4,4)
                     CALL RIMP2_POPLE_DEN                               &
                         (INEW,JNEW,KNEW,LNEW, PHF,PMP2,PTOT,           &
                          DABMAX,L1,L2,Q4,POPLE)

!!!!! FINE SCREENING: INTEGRAL * DENSITY
                     IF(DABMAX*GMAX.LT.CUTOFF2) CYCLE

                     CALL JKDG80                                        &
                         (ONE,DABMAX,INEW,JNEW,KNEW,LNEW,               &
                          JTYPE,IAT,JAT,KAT,LAT)
!!!!! RYS PACKAGE
                  ELSE

                     CALL RIMP2_PRE_RYS                                 &
                         (INEW,JNEW,KNEW,LNEW, NKL,NKL0,                &
                          KK,LL, ISH,JSH,KSH,LSH,                       &
                          NX_INIJG,NX_IIJKLG)

                     IF(NKL.EQ.0) CYCLE
                     IF(SKIPI.AND.SKIPJ.AND.SKIPK.AND.SKIPL) CYCLE

!!!!! RYS DENSITY STORED IN X_IDAB
                     CALL RIMP2_RYS_DEN                                 &
                         (INEW,JNEW,KNEW,LNEW, PHF,PMP2,PTOT,           &
                          X_IDAB,DABMAX,L1,L2,Q4,POPLE)

!!!!! FINE SCREENING: INTEGRAL * DENSITY
                     IF(DABMAX*GMAX.LT.CUTOFF2) CYCLE

                     CALL JKDSPD                                        &
                         (NIJ0,NKL,NKL0,                                &
                          X_ICHRG ,X_IGINT ,X_IFINT ,X_ISINT ,NX_IIJKLG,&
                          X_IGIJKL,X_IGNKL ,X_IGNM  ,X_IXY   ,X_IXZ    ,&
                          X_IYZ   ,X_IX    ,X_IY    ,X_IZ    ,X_ISJ    ,&
                          X_ISK   ,X_ISL   ,X_IB00  ,X_IB01  ,X_IB10   ,&
                          X_IC00  ,X_ID00  ,X_IF00  ,X_IDIJ  ,X_IDKL   ,&
                          X_IDIJSI,X_IDIJSJ,X_IDKLSK,X_IDKLSL,X_IABV   ,&
                          X_ICV   ,X_IRW   ,X_IAAI  ,X_IAAJ  ,X_IBBK   ,&
                          X_IBBL  ,X_IFI   ,X_IFJ   ,X_IFK   ,X_IFL    ,&
                          X_ISII  ,X_ISJJ  ,X_ISKK  ,X_ISLL  ,X_ISIJ   ,&
                          X_ISIK  ,X_ISIL  ,X_ISJK  ,X_ISJL  ,X_ISKL   ,&
                          X_IDAB  , &
                          MAXXYZ,FC,NC,DF,LDF,NBF,DDA,Q4,MINVEC,DABCUT, &
                          DABMAX,ONE)
                  END IF
               ENDDO

            ENDDO
         ENDDO
!$omp end do
!$OMP END PARALLEL
      ENDDO

!!!!! reset DLB
      CALL DDI_DLBRESET


      ALLOCATE(DE_OMP(3,NAT))
      DE_OMP = 0.0D00

!$OMP PARALLEL NUM_THREADS(Nthreads)

!$omp critical(DE_OMP_critical)
      DE_OMP(1:3,1:NAT) = DE_OMP(1:3,1:NAT) + DE(1:3,1:NAT)
!$omp end critical(DE_OMP_critical)

      DEALLOCATE( X_IGIJKL,X_IYZ   ,X_ISK   ,X_IC00  ,X_IDIJSI,X_ICV   ,&
         X_IBBL  ,X_ISII  ,X_ISIK  ,X_IDAB  ,X_IGINT ,X_IGNKL ,X_IX    ,&
         X_ISL   ,X_ID00  ,X_IDIJSJ,X_IRW   ,X_IFI   ,X_ISJJ  ,X_ISIL  ,&
         X_IFINT ,X_IGNM  ,X_IY    ,X_IB00  ,X_IF00  ,X_IDKLSK,X_IAAI  ,&
         X_IFJ   ,X_ISKK  ,X_ISJK  ,X_ISINT ,X_IXY   ,X_IZ    ,X_IB01  ,&
         X_IDIJ  ,X_IDKLSL,X_IAAJ  ,X_IFK   ,X_ISLL  ,X_ISJL  ,X_IXZ   ,&
         X_ISJ   ,X_IB10  ,X_IDKL  ,X_IABV  ,X_IBBK  ,X_IFL   ,X_ISIJ  ,&
         X_ISKL  ,NX_IIJKLG)

!$OMP END PARALLEL

      ! DEALLOCATE(PHF,PMP2,PTOT,X_IXCH,NX_INIJG,X_ICHRG)
      DEALLOCATE(X_IXCH,NX_INIJG,X_ICHRG)

      DE(1:3,1:NAT) = DE_DDI(1:3,1:NAT) + DE_OMP(1:3,1:NAT)
      CALL DDI_GSUMF(1600,DE,3*NAT)

      END





      SUBROUTINE RIMP2_READSTUFF(PHF,PMP2,PTOT,X_IXCH, N_IXCH,NBF2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      double precision :: PHF(NBF2),PMP2(NBF2),PTOT(NBF2)
      double precision :: X_IXCH(N_IXCH)

      COMMON /IOFILE/ IR,IW,IP,IJK,IPK,IDAF,NAV,IODA(950)
      COMMON /MP2DM2/ NOC1,NOC2,ISTEP,NOC3,NOC4
      COMMON /WFNOPT/ SCFTYP,VBTYP,DFTYPE,TDDFTYP,CITYP,CCTYP,          &
                      MPLEVL,MPCTYP
      COMMON /INTOPT/ ISCHWZ,IECP,NECP,IEFLD

      DATA NFT16/16/
      DATA RHF,ROHF/8HRHF     ,8HROHF    /


      ! read scf/mp2 density
      CALL DAREAD(IDAF,IODA,PHF ,NBF2,308,0)
      CALL DAREAD(IDAF,IODA,PMP2,NBF2,307,0)
      ! scf+mp2 density
      PTOT = PHF + PMP2


      ! read exchange screening integral
      IF(ISCHWZ.EQ.1) THEN
         CALL DAREAD(IDAF,IODA,X_IXCH,N_IXCH,54,0)
      ELSE
         X_IXCH = 1.0D00
      END IF

      RETURN

      END !*************************************************************




!
!     -----  AT THIS POINT IT IS GOOD TO REMEMBER THAT    -----
!            -MAXTYP- = HIGHEST SHELL ANGULAR MOMENTUM
!            -MAXFUN- = NUMBER OF FUNCTIONS WITH ANGULAR
!                       MOMENTUM LESS OR EQUAL TO -MAXTYP-
!            -MAXNUM- = NUMBER OF FUNCTIONS WITH ANGULAR
!                       MOMENTUM         EQUAL TO -MAXTYP-
!            -MAXXYZ- = MAXIMUM NUMBER OF PRIMITIVE INTEGRALS
!                       THAT CAN BE HANDLED IN ONE -VECTOR-
!            -NUMXYZ- = ACTUAL MAXIMUM LENGTH OF ONE -VECTOR-
!            -MAXXYZ- = IT IS NUMXYZ/3 . SINCE THE X, Y, AND Z
!                       COMPONENTS ARE TREATED AS A SINGLE VECTOR,
!                       -MAXXYZ- IS THE NUMBER OF (PRIMITIVE-ROOTS)
!                       COMBINATIONS WHICH CAN BE TREATED IN ONE
!                       VECTOR. FOR -SSSS- INTEGRALS WHICH REQUIRE
!                       ONE RYS ROOT, MAXXYZ HAPPENS TO COINCIDE WITH
!                       THE NUMBER OF PRIMITIVE INTEGRALS TREATED IN
!                       ONE VECTOR. FOR -DDDD- INTEGRALS WHICH
!                       REQUIRE FIVE RYS ROOTS, THE NUMBER OF PRIMITIVE
!                       INTEGRALS TREATED IN ONE VECTOR IS -MAXXYZ-/5 .
!

!
!     X(IWFN  ) = WAVEFUNCTION DATA
!     X(IXCH  ) = EXCHANGE INTEGRAL THRESHOLD
!     X(INIJG ) = CHARGE DISTRIBUTION POINTERS
!     X(IGINT ) = ELECTRON REPULSION INTEGRALS
!     X(IFINT ) = FIRST DERIVATIVE INTEGRALS
!     X(ISINT ) = SECOND DERIVATIVE INTEGRALS
!     X(IIJKLG) = INDICES
!     X(IDAB  ) = DENSITY ARRAY -DAB-
!     X(ICHRG ) = CHARGE DISTRIBUTION PARAMETERS
!     X(IXY   ) = TEMPORARY ARRAY WHEN FORMING DERIVATIVE INTEGRALS
!     X(IXZ   ) = TEMPORARY ARRAY WHEN FORMING DERIVATIVE INTEGRALS
!     X(IYZ   ) = TEMPORARY ARRAY WHEN FORMING DERIVATIVE INTEGRALS
!     X(IX    ) = TEMPORARY ARRAY WHEN FORMING DERIVATIVE INTEGRALS
!     X(IY    ) = TEMPORARY ARRAY WHEN FORMING DERIVATIVE INTEGRALS
!     X(IZ    ) = TEMPORARY ARRAY WHEN FORMING DERIVATIVE INTEGRALS
!     X(ISJ   ) = TEMPORARY ARRAY WHEN -SP- SHELLS
!     X(ISK   ) = TEMPORARY ARRAY WHEN -SP- SHELLS
!     X(ISL   ) = TEMPORARY ARRAY WHEN -SP- SHELLS
!     X(IGIJKL) = ( 2-D , 4 CENTERS ) INTEGRALS
!     X(IGNKL ) = ( 2-D , 3 CENTERS ) INTEGRALS
!     X(IGNM  ) = ( 2-D , 2 CENTERS ) INTEGRALS
!     X(IDIJ  ) = CONTRACTION DENSITY FOR -IJ- CHARGE DISTRIBUTION
!     X(IDKL  ) = CONTRACTION DENSITY FOR -KL- CHARGE DISTRIBUTION
!     X(IB00  ) = -B00-
!     X(IB01  ) = -B01-
!     X(IB10  ) = -B10-
!     X(IC00  ) = -C00-
!     X(ID00  ) = -D00-
!     X(IFXX  ) = -F00-
!     X(IDIJSI) = SCALING FACTOR FOR -S- FUNCTION OF AN -SP- II SHELL
!     X(IDIJSJ) = SCALING FACTOR FOR -S- FUNCTION OF AN -SP- JJ SHELL
!     X(IDKLSK) = SCALING FACTOR FOR -S- FUNCTION OF AN -SP- KK SHELL
!     X(IDKLSL) = SCALING FACTOR FOR -S- FUNCTION OF AN -SP- LL SHELL
!     X(IABV  ) = -AB- VECTOR FOR PRIMITIVE INTEGRALS
!     X(ICV   ) = -CV- VECTOR FOR PRIMITIVE INTEGRALS
!     X(IRW   ) = -RW- VECTOR FOR RYS ROOTS AND WEIGHTS
!     X(IAAI  ) = EXPONENT FOR DERIVATIVE OF II SHELL
!     X(IAAJ  ) = EXPONENT FOR DERIVATIVE OF JJ SHELL
!     X(IBBK  ) = EXPONENT FOR DERIVATIVE OF KK SHELL
!     X(IBBL  ) = EXPONENT FOR DERIVATIVE OF LL SHELL
!     X(IFI   ) = FIRST DERIVATIVE WRT. II OF ( 2-D , 4 CENTERS ) INT.
!     X(IFJ   ) = FIRST DERIVATIVE WRT. JJ OF ( 2-D , 4 CENTERS ) INT.
!     X(IFK   ) = FIRST DERIVATIVE WRT. KK OF ( 2-D , 4 CENTERS ) INT.
!     X(IFL   ) = FIRST DERIVATIVE WRT. LL OF ( 2-D , 4 CENTERS ) INT.
!     X(ISII  ) = SECOND DER. WRT II AND II OF ( 2-D , 4 CENTERS ) INT.
!     X(ISJJ  ) = SECOND DER. WRT JJ AND JJ OF ( 2-D , 4 CENTERS ) INT.
!     X(ISKK  ) = SECOND DER. WRT KK AND KK OF ( 2-D , 4 CENTERS ) INT.
!     X(ISLL  ) = SECOND DER. WRT LL AND LL OF ( 2-D , 4 CENTERS ) INT.
!     X(ISIJ  ) = SECOND DER. WRT II AND JJ OF ( 2-D , 4 CENTERS ) INT.
!     X(ISIK  ) = SECOND DER. WRT II AND KK OF ( 2-D , 4 CENTERS ) INT.
!     X(ISIL  ) = SECOND DER. WRT II AND LL OF ( 2-D , 4 CENTERS ) INT.
!     X(ISJK  ) = SECOND DER. WRT JJ AND KK OF ( 2-D , 4 CENTERS ) INT.
!     X(ISJL  ) = SECOND DER. WRT JJ AND LL OF ( 2-D , 4 CENTERS ) INT.
!     X(ISKL  ) = SECOND DER. WRT KK AND LL OF ( 2-D , 4 CENTERS ) INT.
!


      SUBROUTINE RIMP2_JKDMEM                                           &
                (MDER,LOADFM,IADDR,LENGTH,                              &
                 MINXYZ,MAXXYZ,MINVEC,POPLE,MP2)

      use mx_limits, only: mxgtot,mxsh

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      integer,parameter :: LENVEC=255
      double precision,parameter ::                                     &
         ONE=1.0D+00,SQRT3=1.73205080756888D+00,                        &
         SQRT5=2.23606797749979D+00,SQRT7=2.64575131106459D+00



      COMMON /DSHLT / RTOL,DTOL,VTOL1,VTOL2,VTOLS,OUT,DBG
      COMMON /IOFILE/ IR,IW,IP,IJK,IPK,IDAF,NAV,IODA(950)
      COMMON /MACHIN/ NWDVAR,MAXFM,MAXSM,LIMFM,LIMSM
      COMMON /NSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),      &
                      CF(MXGTOT),CG(MXGTOT),CH(MXGTOT),CI(MXGTOT),      &
                      KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),   &
                      KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
      COMMON /OUTPUT/ NPRINT,ITOL,ICUT,NORMF,NORMP,NOPK
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON/SHLBAS/MAXTYP,MAXNUM
      COMMON /SHLNRM/ PNRM(84)
      COMMON/DERPAR/IDER,JDER,KDER,LDER,NDER
!$omp threadprivate(/DERPAR/)
      COMMON /DERMEMALLOCATE/ &
         N_IWFN  ,N_IXCH  ,N_IDAB  ,N_INIJG ,N_IGINT ,N_IFINT ,N_ISINT ,&
         N_IIJKLG,N_ICHRG ,N_IXY   ,N_IXZ   ,N_IYZ   ,N_IX    ,N_IY    ,&
         N_IZ    ,N_ISJ   ,N_ISK   ,N_ISL   ,N_IGIJKL,N_IGNKL ,N_IGNM  ,&
         N_IDIJ  ,N_IDKL  ,N_IB00  ,N_IB01  ,N_IB10  ,N_IC00  ,N_ID00  ,&
         N_IF00  ,N_IDIJSI,N_IDIJSJ,N_IDKLSK,N_IDKLSL,N_IABV  ,N_ICV   ,&
         N_IRW   ,N_IAAI  ,N_IAAJ  ,N_IBBK  ,N_IBBL  ,N_IFI   ,N_IFJ   ,&
         N_IFK   ,N_IFL   ,N_ISII  ,N_ISJJ  ,N_ISKK  ,N_ISLL  ,N_ISIJ  ,&
         N_ISIK  ,N_ISIL  ,N_ISJK  ,N_ISJL  ,N_ISKL

      logical :: SOME,OUT,DBG,NORM,POPLE,MP2,GOPARR,DSKWRK,MASWRK
      integer :: LENSHL(5)

      DATA LENSHL /1,4,10,20,35/






      SOME = MASWRK  .AND.  NPRINT.NE.-5

      ! available memory (X)
      CALL GOTFM(NGOTMX)

      ! ! This is a dirty parameter.
      ! ! Ideally, there should be a
      ! ! function to get available memory
      ! NGOTMX = 100000000


      MAXVEC = LENVEC/3
      MINVEC = IGETGRDVECLEN(MAXVEC)
!
!    ----- THIS PROGRAM FOR DERIVATIVES -----
!
      NDER=MDER
      IF(NDER.NE.1.AND.NDER.NE.2) NDER=0
      NDER0=0
      NDER1=0
      NDER2=0
      IF(OUT              ) NDER0=1
      IF(OUT.AND.NDER.GE.1) NDER1=12
      IF(OUT.AND.NDER.EQ.2) NDER2=78
      IF(        NDER.EQ.2) NDER1=12


      ! max angular momentum
      CALL BASCHK(LMAX)
      MAXTYP=LMAX+1
      MINXYZ=(4*MAXTYP -2 +NDER)/2

      ! number of primitive charge distribution
      NIJG=0
      DO II=1,NSHELL
        DO JJ=1,II
          NIJG=NIJG+KNG(II)*KNG(JJ)
         ENDDO
      ENDDO

      ! normalization constants
      MAXFUN=LENSHL(MAXTYP)
      DO I=1,MAXFUN
         PNRM(I)=ONE
      ENDDO

      NORM=NORMF.NE.1.OR.NORMP.NE.1
      IF(NORM) THEN
         SQRT53=SQRT5/SQRT3
         DO I=1,MAXFUN
           IF (I.EQ.1.OR.I.EQ.2.OR.I.EQ.5.OR.I.EQ.11.OR.I.EQ.21) THEN
              FI = ONE
           ELSE IF (I.EQ.8.OR.I.EQ.20.OR.I.EQ.33) THEN
              FI=FI*SQRT3
           ELSE IF (I.EQ.14) THEN
              FI=FI*SQRT5
           ELSE IF (I.EQ.24) THEN
              FI=FI*SQRT7
           ELSE IF (I.EQ.30) THEN
              FI=FI*SQRT53
           END IF
           PNRM(I)=FI
         ENDDO
      ENDIF
!
!        ALWAYS ALLOCATE -L- SHELLS IF -P- IS HIGHEST ANGULAR MOMENTUM
!
      IDUM=MAXTYP-1
      MAXNUM=((IDUM+1)*(IDUM+2))/2
      IF(MAXNUM.LT.4) MAXNUM=4
      NGIJKL=(MAXNUM**4)

      ! max angular momentum after differentiated
      MODTYP=MAXTYP+NDER

      ! vector length and core pointers
      LVAR=0
      LFIX=LENGTH
      LFIX=LFIX  +(NSHELL*(NSHELL+1))/2

      ! memory for first derivate -SP- functions
      IWFN  = LOADFM + 1
      IXCH  = IWFN + LENGTH
      IDAB  = IXCH + (NSHELL*(NSHELL+1))/2
      ILAST = IDAB + NGIJKL
      INEED = ILAST- LOADFM
      IADDR = ILAST

      ! for allocatable arrays
      N_IWFN = LENGTH
      N_IXCH = (NSHELL*(NSHELL+1))/2
      N_IDAB = NGIJKL


      IF (POPLE.AND.MAXTYP.LT.3.AND.NDER.EQ.1) RETURN
 !        THEN
 !         WRITE(IW,9995) IWFN,IXCH,INEED,NGOTMX
 ! 9995 FORMAT(' SPECIAL -SP- ROUTINES USED.',                            &
 !             ' IWFN,IXCH,INEED,MAXFM = ',4I8)
 !         RETURN
 !      ENDIF

      LFIX=LFIX+( (NSHELL*(NSHELL+1))/2 )*2
      LFIX=LFIX+NGIJKL*(NDER0+NDER1+NDER2)
      LFIX=LFIX+NGIJKL* 4
      LFIX=LFIX+NGIJKL
      LFIX=LFIX+NIJG*15
      LVAR=     ( MODTYP**2       * MODTYP**2       )*3
      LVAR=LVAR+( MODTYP**2       *(MODTYP+MODTYP-1))*3
      LVAR=LVAR+((MODTYP+MODTYP-1)*(MODTYP+MODTYP-1))*3
      LVAR=LVAR+( MODTYP**2                         )*3
      LVAR=LVAR+((MODTYP+MODTYP-1)                  )*3
      LVAR=LVAR+((MODTYP+MODTYP-1)* 3               )*3
      LVAR=LVAR+(  3                                )*3
      LVAR=LVAR+(  9                                )
      LVAR=LVAR+(  4                                )
      LVAR=LVAR+(  5                                )
      LVAR=LVAR+( 18                                )
      LVAR=LVAR+(  2                                )
      LVAR=LVAR+(  4                                )*3
      LVAR=LVAR+( MODTYP**2       * MODTYP**2       )*3*14

 !      MAXXYZ=(NGOTMX-LFIX-1)/LVAR
 !      IF(MP2) MAXXYZ = MIN(MAXXYZ,2*MINXYZ)
 !      MINMEM = (MINXYZ*LVAR)+1+LFIX
 !      IF(MAXXYZ.LT.MINXYZ) THEN
 !         IF(MASWRK) WRITE(IW,9998) MINMEM
 ! 9998 FORMAT(/1X,'NOT ENOUGH MEMORY FOR THE TWO-ELECTRON GRADIENT'/     &
 !              1X,'YOU WILL NEED AT LEAST ',I10,' WORDS.')
 !         CALL ABRT
 !      END IF

      IF(MP2) MAXXYZ = 2*MINXYZ



      IF(.NOT.MP2  .AND.  SOME) WRITE(IW,9994) MINMEM
 9994 FORMAT(/,' THE MINIMUM MEMORY REQUIRED FOR THIS RUN IS ',I10,     &
               ' WORDS.')
      IF(MAXXYZ.GT.MAXVEC) MAXXYZ=MAXVEC


      N_IXCH   =  (NSHELL*(NSHELL+1))/2
      N_INIJG  = (((NSHELL*(NSHELL+1))/2)*2)/NWDVAR
      N_IGINT  =   NGIJKL*NDER0
      N_IFINT  =   NGIJKL*NDER1
      N_ISINT  =   NGIJKL*NDER2
      N_IIJKLG =   NGIJKL*4/NWDVAR
      N_IDAB   =   NGIJKL
      N_ICHRG  =   NIJG  *15
      N_IXY    = (  1                                )*MAXXYZ
      N_IXZ    = (  1                                )*MAXXYZ
      N_IYZ    = (  1                                )*MAXXYZ
      N_IX     = (  1                                )*MAXXYZ
      N_IY     = (  1                                )*MAXXYZ
      N_IZ     = (  1                                )*MAXXYZ
      N_ISJ    = (  1                                )*MAXXYZ
      N_ISK    = (  1                                )*MAXXYZ
      N_ISL    = (  1                                )*MAXXYZ
      N_IGIJKL = ( MODTYP**2       * MODTYP**2       )*MAXXYZ*3
      N_IGNKL  = ( MODTYP**2       *(MODTYP+MODTYP-1))*MAXXYZ*3
      N_IGNM   = ((MODTYP+MODTYP-1)*(MODTYP+MODTYP-1))*MAXXYZ*3
      N_IDIJ   = ( MODTYP**2                         )*MAXXYZ*3
      N_IDKL   = ((MODTYP+MODTYP-1)                  )*MAXXYZ*3
      N_IB00   = ((MODTYP+MODTYP-1)                  )*MAXXYZ*3
      N_IB01   = ((MODTYP+MODTYP-1)                  )*MAXXYZ*3
      N_IB10   = ((MODTYP+MODTYP-1)                  )*MAXXYZ*3
      N_IC00   = (  1                                )*MAXXYZ*3
      N_ID00   = (  1                                )*MAXXYZ*3
      N_IF00   = (  1                                )*MAXXYZ*3
      N_IDIJSI = (  1                                )*MAXXYZ
      N_IDIJSJ = (  1                                )*MAXXYZ
      N_IDKLSK = (  1                                )*MAXXYZ
      N_IDKLSL = (  1                                )*MAXXYZ
      N_IABV   = (  5                                )*MAXXYZ
      N_ICV    = ( 18                                )*MAXXYZ
      N_IRW    = (  2                                )*MAXXYZ
      N_IAAI   = (  1                                )*MAXXYZ*3
      N_IAAJ   = (  1                                )*MAXXYZ*3
      N_IBBK   = (  1                                )*MAXXYZ*3
      N_IBBL   = (  1                                )*MAXXYZ*3
      N_IFI    = ( MODTYP**2       * MODTYP**2       )*MAXXYZ*3
      N_IFJ    = ( MODTYP**2       * MODTYP**2       )*MAXXYZ*3
      N_IFK    = ( MODTYP**2       * MODTYP**2       )*MAXXYZ*3
      N_IFL    = ( MODTYP**2       * MODTYP**2       )*MAXXYZ*3
      N_ISII   = ( MODTYP**2       * MODTYP**2       )*MAXXYZ*3
      N_ISJJ   = ( MODTYP**2       * MODTYP**2       )*MAXXYZ*3
      N_ISKK   = ( MODTYP**2       * MODTYP**2       )*MAXXYZ*3
      N_ISLL   = ( MODTYP**2       * MODTYP**2       )*MAXXYZ*3
      N_ISIJ   = ( MODTYP**2       * MODTYP**2       )*MAXXYZ*3
      N_ISIK   = ( MODTYP**2       * MODTYP**2       )*MAXXYZ*3
      N_ISIL   = ( MODTYP**2       * MODTYP**2       )*MAXXYZ*3
      N_ISJK   = ( MODTYP**2       * MODTYP**2       )*MAXXYZ*3
      N_ISJL   = ( MODTYP**2       * MODTYP**2       )*MAXXYZ*3
      N_ISKL   = ( MODTYP**2       * MODTYP**2       )*MAXXYZ*3


      END




      SUBROUTINE RIMP2_M0(M0,II)
      USE Rimp2_Shared_Data,only: MXSH,MXATM
      implicit double precision(a-h,o-z)
      integer :: M0(48)
      COMMON /SYMTRY/ MAPSHL(MXSH,48),MAPCTR(MXATM,48),                 &
                      T(432),INVT(48),NT

      DO IT = 1,NT
         ID = MAPSHL(II,IT)
         IF (ID.GT.II) CYCLE
         M0(IT) = ID
      ENDDO

      END


      SUBROUTINE RIMP2_M1M2(M1,M2, JJ, M0,II)
      USE Rimp2_Shared_Data,only: MXSH,MXATM
      implicit double precision(a-h,o-z)
      integer :: M0(48),M1(48),M2(48)
      COMMON /SYMTRY/ MAPSHL(MXSH,48),MAPCTR(MXATM,48),                 &
                      T(432),INVT(48),NT

      DO IT = 1,NT
         JD = MAPSHL(JJ,IT)
         IF (JD.GT.II) CYCLE
         ID = M0(IT)
         IF (ID.LT.JD) THEN
            ND = ID
            ID = JD
            JD = ND
         ENDIF
         IF (ID.EQ.II .AND. JD.GT.JJ) CYCLE
         M1(IT) = ID
         M2(IT) = JD
      ENDDO
      END


      SUBROUTINE RIMP2_M3(M3,KK,II)
      USE Rimp2_Shared_Data,only: MXSH,MXATM
      implicit double precision(a-h,o-z)
      integer :: M3(48)
      COMMON /SYMTRY/ MAPSHL(MXSH,48),MAPCTR(MXATM,48),                 &
                      T(432),INVT(48),NT

      DO IT = 1,NT
         KD = MAPSHL(KK,IT)
         IF (KD.GT.II) CYCLE
         M3(IT) = KD
      ENDDO
      END


      SUBROUTINE RIMP2_D42ARRAYS
      USE omp_lib
      USE RIMP2_TMP,only:  NX_INIJG,NX_IIJKLG        ,&
         X_IXCH  ,X_INIJG ,X_ICHRG ,X_IIJKLG,X_IDAB  ,X_IGINT ,X_IFINT ,&
         X_ISINT ,X_IGIJKL,X_IGNKL ,X_IGNM  ,X_IXY   ,X_IXZ   ,X_IYZ   ,&
         X_IX    ,X_IY    ,X_IZ    ,X_ISJ   ,X_ISK   ,X_ISL   ,X_IB00  ,&
         X_IB01  ,X_IB10  ,X_IC00  ,X_ID00  ,X_IF00  ,X_IDIJ  ,X_IDKL  ,&
         X_IDIJSI,X_IDIJSJ,X_IDKLSK,X_IDKLSL,X_IABV  ,X_ICV   ,X_IRW   ,&
         X_IAAI  ,X_IAAJ  ,X_IBBK  ,X_IBBL  ,X_IFI   ,X_IFJ   ,X_IFK   ,&
         X_IFL   ,X_ISII  ,X_ISJJ  ,X_ISKK  ,X_ISLL  ,X_ISIJ  ,X_ISIK  ,&
         X_ISIL  ,X_ISJK  ,X_ISJL  ,X_ISKL

      implicit double precision(a-h,o-z)

      COMMON /DERMEMALLOCATE/ &
         N_IWFN  ,N_IXCH  ,N_IDAB  ,N_INIJG ,N_IGINT ,N_IFINT ,N_ISINT ,&
         N_IIJKLG,N_ICHRG ,N_IXY   ,N_IXZ   ,N_IYZ   ,N_IX    ,N_IY    ,&
         N_IZ    ,N_ISJ   ,N_ISK   ,N_ISL   ,N_IGIJKL,N_IGNKL ,N_IGNM  ,&
         N_IDIJ  ,N_IDKL  ,N_IB00  ,N_IB01  ,N_IB10  ,N_IC00  ,N_ID00  ,&
         N_IF00  ,N_IDIJSI,N_IDIJSJ,N_IDKLSK,N_IDKLSL,N_IABV  ,N_ICV   ,&
         N_IRW   ,N_IAAI  ,N_IAAJ  ,N_IBBK  ,N_IBBL  ,N_IFI   ,N_IFJ   ,&
         N_IFK   ,N_IFL   ,N_ISII  ,N_ISJJ  ,N_ISKK  ,N_ISLL  ,N_ISIJ  ,&
         N_ISIK  ,N_ISIL  ,N_ISJK  ,N_ISJL  ,N_ISKL


      Nthreads = omp_get_max_threads()

!$OMP PARALLEL NUM_THREADS(Nthreads) &
!$omp shared( &
!$omp    N_IWFN  ,N_IXCH  ,N_IDAB  ,N_INIJG ,N_IGINT ,N_IFINT ,N_ISINT ,&
!$omp    N_IIJKLG,N_ICHRG ,N_IXY   ,N_IXZ   ,N_IYZ   ,N_IX    ,N_IY    ,&
!$omp    N_IZ    ,N_ISJ   ,N_ISK   ,N_ISL   ,N_IGIJKL,N_IGNKL ,N_IGNM  ,&
!$omp    N_IDIJ  ,N_IDKL  ,N_IB00  ,N_IB01  ,N_IB10  ,N_IC00  ,N_ID00  ,&
!$omp    N_IF00  ,N_IDIJSI,N_IDIJSJ,N_IDKLSK,N_IDKLSL,N_IABV  ,N_ICV   ,&
!$omp    N_IRW   ,N_IAAI  ,N_IAAJ  ,N_IBBK  ,N_IBBL  ,N_IFI   ,N_IFJ   ,&
!$omp    N_IFK   ,N_IFL   ,N_ISII  ,N_ISJJ  ,N_ISKK  ,N_ISLL  ,N_ISIJ  ,&
!$omp    N_ISIK  ,N_ISIL  ,N_ISJK  ,N_ISJL  ,N_ISKL )

      ALLOCATE(X_IGIJKL(N_IGIJKL))
      ALLOCATE(X_IYZ   (N_IYZ   ))
      ALLOCATE(X_ISK   (N_ISK   ))
      ALLOCATE(X_IC00  (N_IC00  ))
      ALLOCATE(X_IDIJSI(N_IDIJSI))
      ALLOCATE(X_ICV   (N_ICV   ))
      ALLOCATE(X_IBBL  (N_IBBL  ))
      ALLOCATE(X_ISII  (N_ISII  ))
      ALLOCATE(X_ISIK  (N_ISIK  ))
      ALLOCATE(X_IDAB  (N_IDAB  ))

      ALLOCATE(X_IGINT (N_IGINT ))
      ALLOCATE(X_IGNKL (N_IGNKL ))
      ALLOCATE(X_IX    (N_IX    ))
      ALLOCATE(X_ISL   (N_ISL   ))
      ALLOCATE(X_ID00  (N_ID00  ))
      ALLOCATE(X_IDIJSJ(N_IDIJSJ))
      ALLOCATE(X_IRW   (N_IRW   ))
      ALLOCATE(X_IFI   (N_IFI   ))
      ALLOCATE(X_ISJJ  (N_ISJJ  ))
      ALLOCATE(X_ISIL  (N_ISIL  ))

      ALLOCATE(X_IFINT (N_IFINT ))
      ALLOCATE(X_IGNM  (N_IGNM  ))
      ALLOCATE(X_IY    (N_IY    ))
      ALLOCATE(X_IB00  (N_IB00  ))
      ALLOCATE(X_IF00  (N_IF00  ))
      ALLOCATE(X_IDKLSK(N_IDKLSK))
      ALLOCATE(X_IAAI  (N_IAAI  ))
      ALLOCATE(X_IFJ   (N_IFJ   ))
      ALLOCATE(X_ISKK  (N_ISKK  ))
      ALLOCATE(X_ISJK  (N_ISJK  ))

      ALLOCATE(X_ISINT (N_ISINT ))
      ALLOCATE(X_IXY   (N_IXY   ))
      ALLOCATE(X_IZ    (N_IZ    ))
      ALLOCATE(X_IB01  (N_IB01  ))
      ALLOCATE(X_IDIJ  (N_IDIJ  ))
      ALLOCATE(X_IDKLSL(N_IDKLSL))
      ALLOCATE(X_IAAJ  (N_IAAJ  ))
      ALLOCATE(X_IFK   (N_IFK   ))
      ALLOCATE(X_ISLL  (N_ISLL  ))
      ALLOCATE(X_ISJL  (N_ISJL  ))

      ALLOCATE(NX_IIJKLG(N_IIJKLG))
      ALLOCATE(X_IXZ   (N_IXZ   ))
      ALLOCATE(X_ISJ   (N_ISJ   ))
      ALLOCATE(X_IB10  (N_IB10  ))
      ALLOCATE(X_IDKL  (N_IDKL  ))
      ALLOCATE(X_IABV  (N_IABV  ))
      ALLOCATE(X_IBBK  (N_IBBK  ))
      ALLOCATE(X_IFL   (N_IFL   ))
      ALLOCATE(X_ISIJ  (N_ISIJ  ))
      ALLOCATE(X_ISKL  (N_ISKL  ))

!$OMP END PARALLEL

      END



      SUBROUTINE RIMP2_Q4(Q4, M1,M2,M3,II,JJ,KK,LL)

      USE Rimp2_Shared_Data,only: MXSH,MXATM

      implicit double precision(a-h,o-z)

      COMMON /SYMTRY/ MAPSHL(MXSH,48),MAPCTR(MXATM,48),                 &
                      T(432),INVT(48),NT

      integer :: M1(48),M2(48),M3(48)

      N4 = 0

      DO IT = 1,NT
         LD = MAPSHL(LL,IT)
         IF (LD.GT.II) CYCLE
         KD = M3(IT)

         IF (KD.LT.LD) CALL RIMP2SWAP(KD,LD)

         ID = M1(IT)
         JD = M2(IT)
         IF (ID.NE.II .AND. KD.NE.II) CYCLE
         IF (KD.GE.ID) THEN
            IF (KD.NE.ID .OR. LD.GT.JD) THEN
               CALL RIMP2SWAP(ID,KD)
               CALL RIMP2SWAP(JD,LD)
            ENDIF
         ENDIF

         IF(JD.NE.JJ) CYCLE
         IF(KD.NE.KK) CYCLE
         IF(LD.NE.LL) CYCLE

         N4 = N4+1
      ENDDO

      Q4 = NT
      Q4 = Q4 / N4

      END



      SUBROUTINE RIMP2_PRE_POPLE &
                (INEW,JNEW,KNEW,LNEW, &
                 IAT,JAT,KAT,LAT, JTYPE, &
                 ISH,JSH,KSH,LSH, &
                 Q4,KTYPE,KATOM)
      USE omp_lib
      implicit double precision(a-h,o-z)

      double precision,parameter :: &
         ZER=0.0D+00,PT5=0.5D+00,P25=.25D+00,F04=4.0D+00

      COMMON/INDD80/IMAX,JMAX,KKKMAX,LMAX
!$omp threadprivate(/INDD80/)

      integer :: KTYPE(*),KATOM(*)



      AX1=PT5
      IF(ISH.NE.JSH) AX1=AX1+AX1
      IF(KSH.NE.LSH) AX1=AX1+AX1
      IF(ISH.NE.KSH.OR.JSH.NE.LSH) AX1=AX1+AX1
      Q4 = Q4*AX1
      INEW=ISH
      JNEW=JSH
      KNEW=KSH
      LNEW=LSH
      IMAX=KTYPE(INEW)-1
      JMAX=KTYPE(JNEW)-1
      KKKMAX=KTYPE(KNEW)-1
      LMAX=KTYPE(LNEW)-1

      IF (IMAX.LT.JMAX) CALL RIMP2SWAP(INEW,JNEW)
      IF (KKKMAX.LT.LMAX) CALL RIMP2SWAP(KNEW,LNEW)
      IF ((IMAX+JMAX).LT.(KKKMAX+LMAX)) THEN
        CALL RIMP2SWAP(INEW,KNEW)
        CALL RIMP2SWAP(JNEW,LNEW)
      END IF

      IMAX=3*(KTYPE(INEW)-1)+1
      JMAX=3*(KTYPE(JNEW)-1)+1
      KKKMAX=3*(KTYPE(KNEW)-1)+1
      LMAX=3*(KTYPE(LNEW)-1)+1
      JTYPE=(IMAX+JMAX+KKKMAX+KKKMAX+LMAX-2)/3

      IAT = KATOM(INEW)
      JAT = KATOM(JNEW)
      KAT = KATOM(KNEW)
      LAT = KATOM(LNEW)

      END






      SUBROUTINE RIMP2_POPLE_DEN &
                (II,JJ,KK,LL, &
                 PHF,PMP2,PTOT, &
                 DABMAX,L1,L2,Q4,POPLE)

      USE omp_lib
      use mx_limits, only: mxsh,mxgtot,mxao

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      double precision :: PHF(L2),PMP2(L2),PTOT(L2)

      logical :: pople

      COMMON /IJPAIR/ IA(MXAO)
      COMMON /IOFILE/ IR,IW,IP,IJK,IPK,IDAF,NAV,IODA(950)
      COMMON /INDD80/ LA,LB,LC,LD
!$omp threadprivate(/INDD80/)
      COMMON /MP2DM2/ NOC1,NOC2,ISTEP,NOC3,NOC4
      COMMON /MP2PAR/ OSPT,CODEMP,SCSPT,TOL,METHOD,NWORD,MEMPRI,MPPROP, &
                      NACORE,NBCORE,NOA,NOB,NORB,NBF,NOMIT,MOCPHF,MAXITC
      LOGICAL SOME,OUT,DBUG
      COMMON /MP2PRT/ SOME,OUT,DBUG
      COMMON /NSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),      &
                      CF(MXGTOT),CCG(MXGTOT),CH(MXGTOT),CI(MXGTOT),     &
                      KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),   &
                      KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
      LOGICAL IIEQJJ,KKEQLL,IJEQKL,IJGTKL,IJLTKL
      COMMON /SHLEQU/ IIEQJJ,KKEQLL,IJEQKL,IJGTKL,IJLTKL
!$omp threadprivate(/SHLEQU/)
      COMMON /SHLLMN/ IGXYZ(4,35),JGXYZ(4,35),KGXYZ(4,35),LGXYZ(4,35)
!$omp threadprivate(/SHLLMN/)
      COMMON /SHLNRM/ PNRM(84)
      COMMON /GSPG80/ E(4,4,4,4)
!$omp threadprivate(/GSPG80/)
!
      double precision,parameter :: &
         ZER=0.0D+00,PT5=0.5D+00,P25=.25D+00,F04=4.0D+00
!
!     ----- FORM TWO-PARTICLE DENSITY MATRIX FOR MP2 GRADIENT -----
!
      DABMAX= ZER
      MINI= KMIN(II)
      MINJ= KMIN(JJ)
      MINK= KMIN(KK)
      MINL= KMIN(LL)
      LOCI= KLOC(II)-MINI
      LOCJ= KLOC(JJ)-MINJ
      LOCK= KLOC(KK)-MINK
      LOCL= KLOC(LL)-MINL
      NOC = NOC1-1

      DO L=1,LD
         NNU= LOCL+L
         DO K=1,LC
            NMU= LOCK+K
            MUNU=IA(MAX(NMU,NNU))+MIN(NMU,NNU)
            DO J=1,LB
               NSI= LOCJ+J
               DO I=1,LA
                  NLA= LOCI+I
                  LASI=IA(MAX(NLA,NSI))+MIN(NLA,NSI)

                  DF1= ZER

                  MUSI= IA(MAX(NMU,NSI)) + MIN(NMU,NSI)
                  LANU= IA(MAX(NLA,NNU)) + MIN(NLA,NNU)
                  MULA= IA(MAX(NMU,NLA)) + MIN(NMU,NLA)
                  NUSI= IA(MAX(NNU,NSI)) + MIN(NNU,NSI)

                  DF1= DF1+PTOT(MUNU)*PHF(LASI) + PHF(MUNU)*PMP2(LASI)
                  DQ1=     PTOT(MUSI)*PHF(LANU) + PHF(MUSI)*PMP2(LANU) &
                          +PTOT(MULA)*PHF(NUSI) + PHF(MULA)*PMP2(NUSI)
                  DF1= DF1-P25*DQ1

                  DF1= DF1*Q4
                  IF(DABMAX.LT. ABS(DF1)) DABMAX = ABS(DF1)
                  E(I,J,K,L)= DF1
               ENDDO
            ENDDO
         ENDDO
      ENDDO


      END






      SUBROUTINE RIMP2_RYS_DEN &
                (II,JJ,KK,LL, &
                 PHF,PMP2,PTOT, &
                 DAB,DABMAX,L1,L2,Q4,POPLE)

      use mx_limits, only: mxsh,mxgtot,mxao

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      double precision :: PHF(L2),PMP2(L2),PTOT(L2), DAB(*)
      LOGICAL POPLE


      COMMON /IJPAIR/ IA(MXAO)
      COMMON /IOFILE/ IR,IW,IP,IJK,IPK,IDAF,NAV,IODA(950)
      COMMON /INDD80/ LA,LB,LC,LD
!$omp threadprivate(/INDD80/)
      COMMON /MP2DM2/ NOC1,NOC2,ISTEP,NOC3,NOC4
      COMMON /MP2PAR/ OSPT,CODEMP,SCSPT,TOL,METHOD,NWORD,MEMPRI,MPPROP, &
                      NACORE,NBCORE,NOA,NOB,NORB,NBF,NOMIT,MOCPHF,MAXITC
      LOGICAL SOME,OUT,DBUG
      COMMON /MP2PRT/ SOME,OUT,DBUG
      COMMON /NSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),      &
                      CF(MXGTOT),CCG(MXGTOT),CH(MXGTOT),CI(MXGTOT),     &
                      KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),   &
                      KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
      LOGICAL IIEQJJ,KKEQLL,IJEQKL,IJGTKL,IJLTKL
      COMMON /SHLEQU/ IIEQJJ,KKEQLL,IJEQKL,IJGTKL,IJLTKL
!$omp threadprivate(/SHLEQU/)
      COMMON /SHLLMN/ IGXYZ(4,35),JGXYZ(4,35),KGXYZ(4,35),LGXYZ(4,35)
!$omp threadprivate(/SHLLMN/)
      COMMON /SHLNRM/ PNRM(84)


      double precision,parameter :: &
         ZER=0.0D+00,PT5=0.5D+00,P25=.25D+00,F04=4.0D+00


      DABMAX= ZER
      MINI= KMIN(II)
      MINJ= KMIN(JJ)
      MINK= KMIN(KK)
      MINL= KMIN(LL)
      LOCI= KLOC(II)-MINI
      LOCJ= KLOC(JJ)-MINJ
      LOCK= KLOC(KK)-MINK
      LOCL= KLOC(LL)-MINL
      NOC = NOC1-1


      MAXI= KMAX(II)
      MAXJ= KMAX(JJ)
      MAXK= KMAX(KK)
      MAXL= KMAX(LL)
      DO I=MINI,MAXI
         P1I = PNRM(I)
         NLA = LOCI+I
         JMAX= MAXJ
         IF(IIEQJJ) JMAX= I
         DO J=MINJ,JMAX
            P2J = P1I*PNRM(J)
            NSI = LOCJ+J
            LASI= IA(MAX0(NLA,NSI))+MIN0(NLA,NSI)
            KKMAX=MAXK
            IF(IJEQKL) KKMAX= I
            DO K=MINK,KKMAX
               P3K = P2J*PNRM(K)
               NMU = LOCK+K
               LMAX= MAXL
               IF(KKEQLL) LMAX= K
               IF(IJEQKL .AND. K.EQ.I) LMAX= J
               DO L=MINL,LMAX
                  P4L = P3K*PNRM(L)
                  NNU = LOCL+L
                  MUNU= IA(MAX0(NMU,NNU))+MIN0(NMU,NNU)

                  DF1= ZER

                  MUSI= IA(MAX0(NMU,NSI))+MIN0(NMU,NSI)
                  LANU= IA(MAX0(NLA,NNU))+MIN0(NLA,NNU)
                  MULA= IA(MAX0(NMU,NLA))+MIN0(NMU,NLA)
                  NUSI= IA(MAX0(NNU,NSI))+MIN0(NNU,NSI)

                  DF1= DF1+PTOT(MUNU)*PHF(LASI) + PHF(MUNU)*PMP2(LASI)
                  DQ1=     PTOT(MUSI)*PHF(LANU) + PHF(MUSI)*PMP2(LANU) &
                          +PTOT(MULA)*PHF(NUSI) + PHF(MULA)*PMP2(NUSI)
                  DF1= DF1-P25*DQ1


                  DF1= DF1*F04
                  IF(NMU .EQ.NNU ) DF1= DF1*PT5
                  IF(NLA .EQ.NSI ) DF1= DF1*PT5
                  IF(MUNU.EQ.LASI) DF1= DF1*PT5

                  IF(DABMAX.LT. ABS(DF1)) DABMAX= ABS(DF1)
                  IJKL=IGXYZ(1,I)+JGXYZ(1,J)+KGXYZ(1,K)+LGXYZ(1,L)
                  DAB(IJKL)= DF1*P4L
               ENDDO
            ENDDO
         ENDDO
      ENDDO


      RETURN

      END



      SUBROUTINE RIMP2_PRE_RYS &
                (INEW,JNEW,KNEW,LNEW, NKL,NKL0, &
                 KK,LL, ISH,JSH,KSH,LSH, &
                 NX_INIJG,NX_IIJKLG)

      use mx_limits, only: mxao

      implicit double precision(A-H,O-Z)

      COMMON /IJPAIR/ IA(MXAO)
      COMMON /DERSKP/ IIAT,JJAT,KKAT,LLAT,SKIPI,SKIPJ,SKIPK,SKIPL
!$omp threadprivate(/DERSKP/)

      integer :: NX_INIJG(*),NX_IIJKLG(*)

      ! set pointers for -KL- charge distribution
      KKLL=IA(MAX0(KK,LL))+MIN0(KK,LL)
      CALL OEDRD(NX_INIJG,NKL,NKL0,KKLL)

      ! select center for derivative
      CALL JKDATM(ISH,JSH,KSH,LSH)

      ! set indices for shell block
      CALL JKDSHL(ISH,JSH,KSH,LSH)
      CALL JKDNDX(NX_IIJKLG)

      INEW = ISH
      JNEW = JSH
      KNEW = KSH
      LNEW = LSH

      END




      SUBROUTINE RIMP2_INTCUTOFF                                        &
                (CUTOFF,CUTOFF2,DABCUT,                                 &
                 VTOL1,VTOL2,VTOLS,DTOL,RTOL,                           &
                 ITOL,POPLE,KSTART,EX,KNG,NSHELL)

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      integer :: KSTART(*)
      integer :: KNG(*)
      double precision :: EX(*)
      logical :: POPLE

      ! Schwarz screening cutoff
      CUTOFF=1.0D-09
      IF(.NOT.POPLE) CUTOFF=CUTOFF/10.0D00
      CUTOFF2=CUTOFF/2.0D00

      ! two particle density cutoff
      ZBIG = 0.0D00
      DO ISH=1,NSHELL
         I1=KSTART(ISH)
         I2=I1+KNG(ISH)-1
         DO IG=I1,I2
            IF(EX(IG).GT.ZBIG) ZBIG = EX(IG)
         ENDDO
      ENDDO
      DABCUT=1.0D-11
      IF(ZBIG.GT.1.0D06) DABCUT = DABCUT/10.0D00
      IF(ZBIG.GT.1.0D07) DABCUT = DABCUT/10.0D00

      ! Pople cutoff
      VTOL1 = 1.0D-12
      VTOL2 = 1.0D-12
      VTOLS = 1.0D-20
      DTOL = 10.0D00**(-ITOL)
      RTOL = 2.30258D00*ITOL

      END


#endif
