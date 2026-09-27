#define __bug 0

#ifndef _OPENMP
      subroutine OMPRIMP2DRIVER
         COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
         logical :: MASWRK
         IF(MASWRK) WRITE(*,*) "ri-mp2 energy need GMS_OPENMP=true in install.info"
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
      SUBROUTINE OMPRIMP2DRIVER
      use omp_lib
      use Rimp2_Shared_Data
      implicit double precision(a-h,o-z)

      integer,parameter :: NORIMP=8,NOAUBF=2
      double precision,parameter :: P12=1.2D00

      logical :: GOPARR,DSKWRK,MASWRK
      logical :: USEDM, GOSMP, OTHAUX, EXTCAB
      integer :: DDI_NP,DDI_ME,DDI_NN,DDI_MY
      integer :: KQRIMP(NORIMP),KQAUXBF(NOAUBF)
      double precision :: QRIMP(NORIMP),QAUXBF(NOAUBF),CABNAM(2)

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

!!!!! TURN ON/OFF DEBUG FLAG
      BUG=.false.

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
      IF(ISPHER.EQ.1 .OR. ISPHER.EQ.0) IAUXBF=1
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
! support g3mp2, get the second (ibas=2) aux bas in cabnam list
      if (runtyp.eq.g3mp2) then
         ibas=2
         if(maswrk) write(*,911) cabnam(ibas),extcab
911                 format(1X,'rimp2 energy auxiliary basis set',/,&
                           1X,'CABNAM=',A8,1X,'EXTCAB=',L2)
      endif

      CALL RICAUXBAS(CABNAM(ibas),EXTCAB,NAUXCAT)

      NAUXBAS = NAUXCAT


!!!!! SET UP SPHERICAL AUX BASIS SET
#if __bug==1
   WRITE(*,*) "wwww: SPHAUX"
#endif
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
      ALLOCATE(DAux(NAUXSH,MaxNGs,28))
      MaxNGs=MAXVAL(KNG,1)
      ALLOCATE(DAtm(NSHELL,MaxNGs,28))
#if __bug==1
   WRITE(*,*) "wwww: DensFact", ME
   startw = omp_get_wtime()
#endif
      CALL DensFact()


!!!!! MAIN RI-MP2 GRADIENT DRIVER
#if __bug==1
   WRITE(*,*) "wwww: RIMP2_X_ENERGY", ME
#endif
      CALL RIMP2_X_ENERGY                                               &
          (NANGM,NAUXBASD,NAUXBAS,NAUXSH,                               &
           NCOR,NACT,NVIR,NORB,NBF,                                     &
           STOL,VTOL,IVMTD,OTHAUX)


      CALL DFINAL(1)
      IREST = 0

      IF(ISPHER.EQ.1 .OR. ISPHER.EQ.0) DEALLOCATE(LOCSPH)
      DEALLOCATE(DAux,DAtm)

#if __bug==1
      endw = omp_get_wtime()
      WRITE(*,'(A40,I5,F10.1)') "TIME SUB RIMP2_X_ENERGY",ME,endw-startw
#endif

      RETURN

      END !***********************************************************





!>*module rimp2grd   *deck RIMP2_X_ENERGY
!>
!>     @brief   ri-mp2 grad driver
!>     @author  buu
!>     @date    sept 8, 2018
!>
!>     @detail  ri-mp2 grad driver
!>

      SUBROUTINE RIMP2_X_ENERGY                                         &
                (NANGM,NAUXBASD,NAUXBAS,NAUXSH,                         &
                 NCOR,NACT,NVIR,NORB,NBF,                               &
                 STOL,VTOL,IVMTD,OTHAUX)

      use Rimp2_Shared_Data,only: D_V,D_B,D_GNAX,MXATM,BUG,DOFMO,DOFMO1R

      implicit double precision(a-h,o-z)

      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,                    &
                      ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /IOFILE/ IR,IW,IP,IJK,IPK,IDAF,NAV,IODA(950)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK

      logical :: ULRange
      logical :: GOPARR,DSKWRK,MASWRK
      logical :: OTHAUX

      double precision,allocatable,dimension(:) :: EIG
      double precision,allocatable,dimension(:,:) :: VEC,VDX,VXX, B32



                        !/////////////////////
                        !///    INIT    //////
                        !/////////////////////


!!!!! DEBUGING FLAG
      DOFMO=.TRUE.
      DOFMO1R=DOFMO.AND.(NPROC.EQ.1)


!!!!! SYNCHRONIZE DDI PROCESSES
      CALL DDI_SYNC(5123)


!!!!! SOME PARAMETERS
      NOCC=NCOR+NACT
      NORB=NCOR+NACT+NVIR

      NXO=NAUXBASD*NOCC
      NXV=NAUXBASD*NVIR
      NXN=NAUXBASD*NORB


                        !/////////////////////
                        !///  B MATRIX  //////
                        !/////////////////////


!!!!! FORM VDX
#if __bug==1
   WRITE(*,*) "wwww: RIMP2_1R_X_VXX"
   startw = omp_get_wtime()
#endif

      ALLOCATE(VDX(NAUXBAS,NAUXBAS))
      ! CALL RIMP2_FORM_VDX(VDX, STOL,VTOL,NAUXBAS,NAUXBASD,IVMTD,OTHAUX)
      ! FORM VXX
      CALL RIMP2_1R_X_VXX(VDX, NAUXBAS)

      ! CHOLESKY DEC VXX
#if __bug==1
   WRITE(*,*) "wwww: RIMP2_CHOLESKY_DECOMPOSE_VXX"
   startw = omp_get_wtime()
#endif

      NAUXBAS1=NAUXBAS
      CALL RIMP2_CHOLESKY_DECOMPOSE_VXX (VDX, NAUXBAS1)
      NAUXBASD = NAUXBAS1

      ! CALL DDI_SYNC(5127)

#if __bug
   endw = omp_get_wtime()
   WRITE(*,'(A40,I5,F10.1)') "TIME SUB RIMP2_FORM_VDX ",ME, endw-startw
#endif


!!!!! READ AND SYM TRANSFORM MO VECTORS
      ALLOCATE(VEC(NBF,NBF))
      CALL DAREAD(IDAF,IODA,VEC,NBF*NBF,15,0)
      CALL RIMP2_SYM_TRANS_MO(VEC,NBF,NCOR,NORB-NCOR)


!!!!! FORM I32 AND DO MO-TRANSFORMATION
#if __bug==1
      startw = omp_get_wtime()
#endif
      IF(NPROC.EQ.1) THEN
#if __bug==1
   WRITE(*,*) "wwww: RIMP2_1R_FAT_X_FORM_MO_I32"
#endif
         ! store B32 on node
         ALLOCATE(B32(NAUXBAS*NVIR,NACT))
         CALL RIMP2_1R_FAT_X_FORM_MO_I32(B32,VEC,NCOR,NACT,NVIR,NBF, NAUXBAS)
      ELSE
      ! the second option can be storing B32 on SMP memory, but I am lazy now
#if __bug==1
   WRITE(*,*) "wwww: RIMP2_FAT_X_FORM_MO_I32"
#endif
         ! store B32 in ddi array
         CALL DDI_CREATE(NAUXBAS*NVIR,NACT, D_B)
         CALL RIMP2_FAT_X_FORM_MO_I32(VEC,NCOR,NACT,NVIR,NBF, NAUXBAS)
      ! ELSE
      !    IF(BUG) WRITE(*,*) "wwww entering RIMP2_X_FORM_MO_I32"
      !    ! CALL RIMP2_X_FORM_MO_I32(VEC,NBF,NACT,NVIR,NAUXBAS)
      ENDIF
#if __bug==1
   endw = omp_get_wtime()
   WRITE(*,'(A40,I5,F10.1)') "TIME SUB RIMP2_X_FORM_MO_I32",ME, endw-startw
#endif
      CALL DDI_SYNC(5131)




!!!!! COMBINE I32 WITH THE DECOMPOSED INVERSED V-MATRIX
!!!!! THE MATRIX B IS FORMED AND STORED IN DDI ARRAY (D_B)
      startw = omp_get_wtime()

      IF(NPROC.EQ.1) THEN
#if __bug==1
   WRITE(*,*) "wwww: RIMP2_1R_FAT_COMBINE_VDX_I32"
#endif
         CALL RIMP2_1R_FAT_COMBINE_VDX_I32(B32, VDX,NAUXBAS,NAUXBASD,NACT,NVIR,NBF)
      ELSE
#if __bug==1
   WRITE(*,*) "wwww: RIMP2_FAT_COMBINE_VDX_I32"
#endif
         CALL RIMP2_FAT_COMBINE_VDX_I32(VDX,NAUXBAS,NAUXBASD,NACT,NVIR,NBF)
      ! ELSE
      !    IF(BUG) WRITE(*,*) "wwww entering RIMP2_COMBINE_VDX_I32"
      !    CALL RIMP2_COMBINE_VDX_I32(VDX, NAUXBAS,NAUXBASD,NBF)
      ENDIF
#if __bug==1
   endw = omp_get_wtime()
   WRITE(*,'(A40,I5,F10.1)') "TIME SUB RIMP2_COMBINE_VDX_I32",ME, endw-startw
#endif

!!!!! READ MO ENERGY
#if __bug==1
   WRITE(*,*) "wwww: DAREAD"
#endif
      ALLOCATE(EIG(NBF))
      CALL DAREAD(IDAF,IODA,EIG,NBF,17,0)


                        !/////////////////////
                        !///   ENERGY    /////
                        !/////////////////////


!!!!! ACT-ACT OCCUPIED LOOP
#if __bug==1
   WRITE(*,*) "computing energy"
   startw = omp_get_wtime()
#endif

      IF(NPROC.EQ.1) THEN
#if __bug==1
   WRITE(*,*) "wwww: RIMP2_1R_FAT_ENERGY"
#endif
         CALL RIMP2_1R_FAT_ENERGY                                       &
             (B32, EIG,VEC, NANGM,NAUXBASD,NAUXBAS,                     &
              NCOR,NACT,NOCC,NVIR,NORB,NBF)
      ELSEIF(DOFMO) THEN
#if __bug==1
   WRITE(*,*) "wwww: RIMP2_FAT_ENERGY"
#endif
         CALL RIMP2_FAT_ENERGY                                          &
             (EIG,VEC, NANGM,NAUXBASD,NAUXBAS,                          &
              NCOR,NACT,NOCC,NVIR,NORB,NBF)
         CALL DDI_DESTROY(D_B)
      ! ELSE
      !    CALL RIMP2_ENERGY                                              &
      !        (EIG,VEC, NANGM,NAUXBASD,NAUXBAS,                          &
      !         NCOR,NACT,NOCC,NVIR,NORB,NBF)
      ENDIF
#if __bug==1
   endw = omp_get_wtime()
   WRITE(*,'(A40,I5,F10.1)') "TIME SUB RIMP2_FAT_ENERGY",ME,endw-startw
#endif


      END !!! done :P







!>*module rimp2grd   *deck RIMP2_ENERGY
!>
!>     @brief   form density matrices
!>     @author  buu
!>     @date    sept 8, 2018
!>
!>     @detail  loop over active-active occ MOs
!>              to (partially) form density matrices
!>

      SUBROUTINE RIMP2_ENERGY                                          &
                (EIG,VEC, NANGM,NAUXBASD,NAUXBAS,           &
                 NCOR,NACT,NOCC,NVIR,NORB,NBF)

      USE omp_lib
      USE Rimp2_Shared_Data,only: D_B,BUG,MXRT

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      COMMON /IOFILE/ IR,IW,IP,IJK,IPK,IDAF,NAV,IODA(950)
      COMMON /ENRGYS/ ENUCR,EELCT,ETOT,SZ,SZZ,ECORE,ESCF,EERD,E1,E2,    &
                      VEN,VEE,EPOT,EKIN,ESTATE(MXRT),STATN,EDFT(3),EDISP
      COMMON /ENRGMP/ EMP2,EMP3,EMP4,EMP2A
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /MACHIN/ NWDVAR,MAXFM,MAXSM,LIMFM,LIMSM


!!!!! INPUT
      double precision :: VDX(NAUXBASD,NAUXBAS)
      double precision :: EIG(NBF)
      double precision :: VEC(NBF,NBF)

!!!!! LOCAL
      double precision,allocatable,dimension(:,:) :: YDV,eij,eab,QVV,TVV
      double precision,allocatable,dimension(:,:,:) :: BDVA

      integer,allocatable :: ddiWSA(:,:),I_DLB_INDEX(:)
      logical :: ULRange,LOMP
      logical :: GOPARR,DSKWRK,MASWRK


                     !!!!!!!!!!!!!!!!!!!!!!!!!!!
                     !!!!!      INIT       !!!!!
                     !!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!! TURN OFF DYNAMIC THREADS TO KEEP threadprivate DATA CONSISTENT
      CALL OMP_SET_DYNAMIC(.FALSE.)

      ! parameters
      NOCC=NCOR+NACT
      NORB=NCOR+NACT+NVIR

      NDO=NAUXBASD*NOCC
      NDV=NAUXBASD*NVIR
      NDN=NAUXBASD*NORB


!!!!! SET BUFFER ZONE FOR 2 LAYER PARALLEL
      N_DLB_BUFF=0
      IF(NACT/NPROC .gt. 20) N_DLB_BUFF=1 ! DLB buffer

!!!!! OMP THREADS
      Nthreads_env=omp_get_max_threads()
      Nthreads=Nthreads_env
      IF(BUG) WRITE(*,'(A40,2I5)') "NUM_THREAD SUB RIMP2_ENERGY", Nthreads,Nthreads_env

!!!!! INIT
      E2=0.0D00

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

!!!!! GET MATRIX B
      ALLOCATE(BDVA(NAUXBASD,NVIR,NACT))
      CALL DDI_GET(D_B, 1+NDO,NDV+NDO, 1+NCOR,NACT+NCOR, BDVA)



                     !!!!!!!!!!!!!!!!!!!!!!!!!!!
                     !!!!!      MAIN       !!!!!
                     !!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!! BUFFERS
      ALLOCATE(YDV(NAUXBASD,NVIR))
      ALLOCATE(QVV(NVIR,NVIR))

      startw_static = omp_get_wtime()

      ! DO JACT=1,NACT
      DO JACT=LddiActStart,LddiActEnd !- N_DLB_BUFF
            JJ=JACT+NCOR

            YDV=0.0D00

            DO IACT=1, NACT
               CALL DGEMM('T','N', NVIR,NVIR,NAUXBASD,                        &
                           1.0D+00, BDVA(1,1,IACT),NAUXBASD,                  &
                                    BDVA(1,1,JACT),NAUXBASD,                  &
                           0.0D+00, QVV,NVIR)

               FORALL(IA=1:NVIR,IB=1:NVIR)
                  QVV(IA,IB)=QVV(IA,IB)/(eij(IACT,JACT)-eab(IA,IB))
               ENDFORALL

               TVV = QVV + QVV
               TVV = TVV - transpose(QVV)

!    YDV <- BJ x transpose(TVV)
               CALL DGEMM('N','N', NAUXBASD,NVIR,NVIR,                           &
                           1.0D00, BDVA(1,1,IACT),NAUXBASD,                                  &
                                   TVV,NVIR,                                     &
                           1.0D00, YDV,NAUXBASD)
            ENDDO
!**** CORRELATION ENERGY
               E2=E2+DDOT(NAUXBASD*NVIR, BDVA(1,1,JACT),1, YDV,1)

      ENDDO !JACT
      endw_static = omp_get_wtime()
      IF(BUG) WRITE(*,'(A40,I5,F10.1)') "TIME ACT_ACT_STATIC: ",ME, endw_static-startw_static


                     !!!!!!!!!!!!!!!!!!!!!!!!!!!
                     !!!!!      CLOSE      !!!!!
                     !!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!! REDUCTION:: CORRELATION ENERGY
      CALL DDI_GSUMF(5999,E2,1)

!!!!! RHF RI-MP2 ENERGY
      EMP2 = E2+ESCF

!!!!! PRINT OUT MP2 ENERGY
      IF(MASWRK) WRITE(IW,100) ESCF,E2,EMP2
 100  FORMAT(/1X,'RHF RI MP2 ENERGY',/,                                 &
             12X,'   ESCF=',1X,F20.10/                                  &
             12X,'   E(2)=',1X,F20.10/                                  &
             12X,' E(MP2)=',1X,F20.10 ,/)

      endw = omp_get_wtime()
      IF(BUG) WRITE(*,'(A40,I5,F10.1)') "TIME ACT_ACT: OMP REDUCTION",ME, endw-startw



      END !***********************************************************






!>*module rimp2grd   *deck RIMP2_X_FORM_MO_I32
!>
!>     @brief   form 3-2ERIs
!>     @author  buu
!>     @date    sept 8, 2018
!>
!>     @detail  form 3-2ERIs using Rys Quadrature method
!>
      SUBROUTINE RIMP2_FAT_X_FORM_MO_I32(VEC,NCOR,NACT,NVIR,NBASIS,NAUXBAS)

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
      double precision,allocatable,dimension(:) :: T1,GG,TL,XX1,XX2
      double precision,allocatable,dimension(:,:) :: B32
      double precision,allocatable,dimension(:,:,:) :: I32

      integer,allocatable,dimension(:,:) :: ddiWS
      logical :: KANDI, IANDJ,LOMP
      logical :: GOPARR,DSKWRK,MASWRK
      double precision :: CC(3,3)

      double precision,save :: w21,w32,w43
!$omp threadprivate(w21,w32,w43)


                     !!!!!!!!!!!!!!!!!!!!!!!!!!!
                     !!!!!      INIT       !!!!!
                     !!!!!!!!!!!!!!!!!!!!!!!!!!!





      IF(ME.GT.NAUXSH-1) GOTO 120

!!!!! NUMBER OF BFXs IN A SHELL
      MXBSH=MaxAUXBFX
      IF(MaxAUXANG < MaxATMANG) MXBSH=MaxATMBFX
      MaxT12=MAX(NBASIS*NBASIS, MXBSH*MXBSH*MXBSH)
      ISIZE = MAXVAL(KAUXNG)


!!!!! DDI WORK DISTRIBUTION :: NAUXSH LOOP
      IF(.FALSE.) THEN
         ALLOCATE(ddiWS(0:NPROC-1,2))
         CALL WorkSharedArray(ddiWS,NPROC,1,NAUXSH,LOMP)
         LddiAuxShStart=ddiWS(ME,1)
         LddiAuxShEnd=ddiWS(ME,2)
         DEALLOCATE(ddiWS)
      ELSE
         ! this work distribution seems to be better
         CALL RIMP2_AUXSHAUX &
         (LddiAuxShStart,LddiAuxShEnd,LddiAuxStart,LddiAuxEnd, NAUXBAS)
      ENDIF

      ! WRITE(*,'(A7,6I5)') "SPLITTTTT",LddiAuxShStart,LddiAuxShEnd,LddiAuxShEnd-LddiAuxShStart+1, &
      !                                 LddiAuxStart,LddiAuxEnd, LddiAuxEnd-LddiAuxStart+1

!!!!! B32 BUFFER
      IBG=NVIR*(LddiAuxStart-1)+1
      IED=NVIR*LddiAuxEnd
      ALLOCATE(B32(IBG:IED,1:NACT))

!!!!! OMP THREADS
      Nthreads_env=omp_get_max_threads()
      Nthreads=MIN0(Nthreads_env,LddiAuxShEnd-LddiAuxShStart+1)

#if __bug==1
      WRITE(*,'(A40,2I5)') "NUM_THREAD SUB RIMP2_FAT_X_FORM_MO_I32", Nthreads, Nthreads_env
#endif


                     !!!!!!!!!!!!!!!!!!!!!!!!!!!
                     !!!!!      MAIN       !!!!!
                     !!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! w21=0.0
      ! w1=omp_get_wtime()

!$OMP PARALLEL NUM_THREADS(NTHREADS)                                    &
!$omp default(none)                                                     &
!$omp shared(LddiAuxShStart,LddiAuxShEnd,Nthreads,MXBSH,ISIZE,MaxT12,   &
!$omp        LddiAuxStart,LddiAuxEnd,ME, &
!$omp        KAUXAT,C,EXAUX,KAUXTY,KAUXMI,KAUXMX,KAUXNG,KAUXST,KAUXLO,  &
!$omp        NAUXSH,NAUXBAS, FlgSphAux,LOCSPH,MaxAUXBFX,                &
!$omp        NSHELL,EX,KATOM,KTYPE,KMIN,KMAX,KNG,KSTART,KLOC,           &
!$omp        NCOR,NACT,NVIR,NBASIS,VEC,B32, D_B)                        &
!$omp private(KAUXSH,KAT,CC,LANGK,MINK,MAXK,NK,NGK,L1,LOCK,             &
!$omp         MINKS,MAXKS,LOCKS,NKS,MINKT,MAXKT,NKT,LOCKT,              &
!$omp         ISH,IAT,LANGI,MINI,MAXI,NI,NGI,I1,LOCI,KANDI,             &
!$omp         JSH,JAT,LANGJ,MINJ,MAXJ,NJ,NGJ,J1,LOCJ,IANDJ,             &
!$omp         IauxStart,IauxEnd,IBG,IED,                                &
!$omp         I32,GG,TL,T1,XX1,XX2)


!!!!! 3-2ERI ARRAY
      ALLOCATE(I32(NBASIS,NBASIS,MaxAUXBFX))

!!!!! BUFFERS
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

                        CALL fmoRIMP2_X_FORM_I32                        &
                            (GG, KAUXSH,ISH,JSH,                        &
                             CC, EX(I1),EX(J1),EXAUX(L1),               &
                             LANGI,LANGJ,LANGK,                         &
                             NI,NJ,NK,                                  &
                             NGI,NGJ,NGK,                               &
                             MINI,MAXI,MINJ,MAXJ,MINK,MAXK,             &
                             XX1,XX2)


                        CALL CP3CAUXINT(GG,I32,TL,T1,                   &
                             NBASIS,LOCI,LOCJ,                          &
                             MINI,MAXI,MINJ,MAXJ,MINK,MAXK,             &
                             NI,NJ,NK,MINKS,MAXKS,                      &
                             NKS,FlgSphAux,IANDJ)

                  ENDDO
            ENDDO

!!!!! MO TRANSFORMED I32

            CALL RIMP2_FAT_MO_TRANSFORMED_I32                           &
                (I32,VEC,T1,                                            &
                 LOCKT,MINKT,MAXKT,NCOR,NACT,NVIR,NBASIS,NKT)

!!!!! DDI PUT I32

            DO LL=IauxStart,IauxEnd
               IBG=NVIR*(LL-1)+1
               IED=NVIR*LL
               CALL DCOPY(NVIR*NACT,I32(1,1,LL-IauxStart+1),1,B32(IBG:IED,1:NACT),1)
            ENDDO

      ENDDO

!$omp end do

      DEALLOCATE(I32, T1,GG,TL,XX1,XX2)

!$OMP END PARALLEL

      ! w2=omp_get_wtime()
      ! w21=w21 + w2-w1
      ! WRITE(*,'(A10,F5.1,I5)') "kkkkkkkk", w21,ME


      ! w65=0.0
      ! w5=omp_get_wtime()

      IBG=NVIR*(LddiAuxStart-1)+1
      IED=NVIR*LddiAuxEnd
      CALL DDI_PUT(D_B, IBG,IED, 1,NACT, B32(IBG:IED,1:NACT))

      ! w6=omp_get_wtime()
      ! w65=w65+ w6-w5
      ! WRITE(*,'(A10,F5.1,I5)') "jjjjjjjjj", w65,ME


  120 CONTINUE


      END !***********************************************************




      SUBROUTINE fmoRIMP2_X_FORM_I32                                    &
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

      double precision :: TMP(3),VAL(3)
      ! double precision,allocatable :: DUM2(:)
      ! ! integer,allocatable :: JLJ(:)

      ! allocate(DUM2(1:NI))

      ! ! ALLOCATE(JLJ(0:LJ))
      ! ! FORALL(J=0:LJ)
      ! !    JLJ(J)=J
      ! ! ENDFORALL

      TOL = ITOLRI*RLN10

      VAL(1:3) = CC(1:3,1)-CC(1:3,2)

      RAB2 = VAL(1)*VAL(1)
      RAB2 = RAB2 + VAL(2)*VAL(2)
      RAB2 = RAB2 + VAL(3)*VAL(3)

      DO J = 0, LJ
         CAB(1:3,J) = VAL(1:3)**J
      ENDDO

      LLI(1,1:NI) = LX(MINI:MAXI)
      LLI(2,1:NI) = LY(MINI:MAXI)
      LLI(3,1:NI) = LZ(MINI:MAXI)

      LLJ(1,1:NJ) = LX(MINJ:MAXJ)
      LLJ(2,1:NJ) = LY(MINJ:MAXJ)
      LLJ(3,1:NJ) = LZ(MINJ:MAXJ)

      LLK(1,1:NK) = LX(MINK:MAXK)
      LLK(2,1:NK) = LY(MINK:MAXK)
      LLK(3,1:NK) = LZ(MINK:MAXK)


!!!!! DENSITY FACTOR
      AKINVV(1:NGK) = 1.0D00/EXK(1:NGK)
      DO K=1,NK
         DK(K,1:NGK)=DAux(KAUXSH,1:NGK,K)
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

         TMP(1:3) = AI*CC(1:3,1)

         AIRAB = AI*RAB2

!!!!! J-SHELL
         DO JG = 1, NGJ

            AJ = EXJ(JG)
            AB = AI+AJ
            ABINV = 1.0D00/AB

            DUM = AJ*ABINV*AIRAB

            IF (DUM > TOL) CYCLE

            EXPB = EXP(-DUM)

            P(1:3) = (TMP(1:3)+AJ*CC(1:3,2))*ABINV

            PAV(1:3) = P(1:3)-CC(1:3,1)

            CPV(1:3) = CC(1:3,3)-P(1:3)

            RR = CPV(1)*CPV(1)
            RR = RR + CPV(2)*CPV(2)
            RR = RR + CPV(3)*CPV(3)

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






!>*module rimp2grd   *deck RIMP2_X_FORM_MO_I32
!>
!>     @brief   form 3-2ERIs
!>     @author  buu
!>     @date    sept 8, 2018
!>
!>     @detail  form 3-2ERIs using Rys Quadrature method
!>
      SUBROUTINE RIMP2_1R_FAT_X_FORM_MO_I32(B32,VEC,NCOR,NACT,NVIR,NBASIS,NAUXBAS)

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


!!!!! OUTPUT
      double precision :: B32(NAUXBAS*NVIR,NACT)


!!!!! INPUT
      double precision :: VEC(NBASIS,NBASIS)


!!!!! LOCAL
      double precision,allocatable,dimension(:) :: T1,GG,TL,XX1,XX2
      double precision,allocatable,dimension(:,:,:) :: I32

      integer,allocatable,dimension(:,:) :: ddiWS
      logical :: KANDI, IANDJ,LOMP
      logical :: GOPARR,DSKWRK,MASWRK
      double precision :: CC(3,3)


                     !!!!!!!!!!!!!!!!!!!!!!!!!!!
                     !!!!!      INIT       !!!!!
                     !!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!! NUMBER OF BFXs IN A SHELL
      MXBSH=MaxAUXBFX
      IF(MaxAUXANG < MaxATMANG) MXBSH=MaxATMBFX
      MaxT12=MAX(NBASIS*NBASIS, MXBSH*MXBSH*MXBSH)
      ISIZE = MAXVAL(KAUXNG)

!!!!! OMP THREADS
      Nthreads=omp_get_max_threads()
      IF(BUG) WRITE(*,'(A40,2I5)') "NUM_THREAD SUB RIMP2_X_FORM_MO_I32", Nthreads


                     !!!!!!!!!!!!!!!!!!!!!!!!!!!
                     !!!!!      MAIN       !!!!!
                     !!!!!!!!!!!!!!!!!!!!!!!!!!!



!$OMP PARALLEL NUM_THREADS(NTHREADS)                                    &
!$omp default(none)                                                     &
!$omp shared(LddiAuxShStart,LddiAuxShEnd,Nthreads,MXBSH,ISIZE,MaxT12,   &
!$omp        KAUXAT,C,EXAUX,KAUXTY,KAUXMI,KAUXMX,KAUXNG,KAUXST,KAUXLO,  &
!$omp        NAUXSH,NAUXBAS, FlgSphAux,LOCSPH,MaxAUXBFX,                &
!$omp        NSHELL,EX,KATOM,KTYPE,KMIN,KMAX,KNG,KSTART,KLOC,           &
!$omp        NCOR,NACT,NVIR,NBASIS,VEC,B32)                             &
!$omp private(KAUXSH,KAT,CC,LANGK,MINK,MAXK,NK,NGK,L1,LOCK,             &
!$omp         MINKS,MAXKS,LOCKS,NKS,MINKT,MAXKT,NKT,LOCKT,              &
!$omp         ISH,IAT,LANGI,MINI,MAXI,NI,NGI,I1,LOCI,KANDI,             &
!$omp         JSH,JAT,LANGJ,MINJ,MAXJ,NJ,NGJ,J1,LOCJ,IANDJ,             &
!$omp         IauxStart,IauxEnd,IBG,IED,                                &
!$omp         I32,GG,TL,T1,XX1,XX2)

!!!!! 3-2ERI ARRAY
      ALLOCATE(I32(NBASIS,NBASIS,MaxAUXBFX))

!!!!! BUFFERS
      ALLOCATE(T1(MaxT12))
      ALLOCATE(GG(MaxT12))
      ALLOCATE(TL(MXBSH*MXBSH))
      ALLOCATE(XX1(ISIZE))
      ALLOCATE(XX2(ISIZE*28))

!!!!! K-SHELL (AUX)
!$omp do schedule(DYNAMIC)

      DO KAUXSH=1,NAUXSH

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

                        CALL CP3CAUXINT(GG,I32,TL,T1,                   &
                             NBASIS,LOCI,LOCJ,                          &
                             MINI,MAXI,MINJ,MAXJ,MINK,MAXK,             &
                             NI,NJ,NK,MINKS,MAXKS,                      &
                             NKS,FlgSphAux,IANDJ)

                  ENDDO
            ENDDO

!!!!! MO TRANSFORMED I32

            CALL RIMP2_FAT_MO_TRANSFORMED_I32                           &
                (I32,VEC,T1,                                            &
                 LOCKT,MINKT,MAXKT,NCOR,NACT,NVIR,NBASIS,NKT)

!!!!! DDI PUT I32

            DO LL=IauxStart,IauxEnd
               IBG=NVIR*(LL-1)+1
               IED=NVIR*LL
               CALL DCOPY(NVIR*NACT, I32(1,1,LL-IauxStart+1),1,         &
                                     B32(IBG:IED,1:NACT),1)
            ENDDO

      ENDDO

!$omp end do

      DEALLOCATE(I32, T1,GG,TL,XX1,XX2)

!$OMP END PARALLEL


      END !***********************************************************






!>*module rimp2grd   *deck RIMP2_X_FORM_MO_I32
!>
!>     @brief   form 3-2ERIs
!>     @author  buu
!>     @date    sept 8, 2018
!>
!>     @detail  form 3-2ERIs using Rys Quadrature method
!>
      SUBROUTINE RIMP2_X_FORM_MO_I32_test(B32,VEC,NCOR,NACT,NVIR,NBASIS,NAUXBAS)

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


!!!!! OUTPUT
      double precision :: B32(NAUXBAS*NVIR,NACT)


!!!!! INPUT
      double precision :: VEC(NBASIS,NBASIS)


!!!!! LOCAL
      double precision,allocatable,dimension(:) :: T1,GG,TL,XX1,XX2
      double precision,allocatable,dimension(:,:,:) :: I32

      integer,allocatable,dimension(:,:) :: ddiWS
      logical :: KANDI, IANDJ,LOMP
      logical :: GOPARR,DSKWRK,MASWRK
      double precision :: CC(3,3)


                     !!!!!!!!!!!!!!!!!!!!!!!!!!!
                     !!!!!      INIT       !!!!!
                     !!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!! NUMBER OF BFXs IN A SHELL
      MXBSH=MaxAUXBFX
      IF(MaxAUXANG < MaxATMANG) MXBSH=MaxATMBFX
      MaxT12=MAX(NBASIS*NBASIS, MXBSH*MXBSH*MXBSH)
      ISIZE = MAXVAL(KAUXNG)


!!!!! DDI WORK DISTRIBUTION :: NAUXSH LOOP
      ! ALLOCATE(ddiWS(0:NPROC-1,2))
      ! CALL WorkSharedArray(ddiWS,NPROC,1,NAUXSH,LOMP)
      ! LddiAuxShStart=ddiWS(ME,1)
      ! LddiAuxShEnd=ddiWS(ME,2)
      ! DEALLOCATE(ddiWS)

      CALL RIMP2_DDIAUXSH(LddiAuxShStart,LddiAuxShEnd, ME,NPROC,NAUXBAS)



!!!!! ALTERNATIVELY :: NAUXBAS
      CALL RIMP2_AUXLOC(LddiAuxStart, LddiAuxShStart,.TRUE.)
      CALL RIMP2_AUXLOC(LddiAuxEnd, LddiAuxShEnd,.FALSE.)


!!!!! OMP THREADS
      Nthreads_env=omp_get_max_threads()
      Nthreads=MIN0(Nthreads_env,LddiAuxShEnd-LddiAuxShStart+1)
      IF(BUG) WRITE(*,'(A40,2I5)') "NUM_THREAD SUB RIMP2_X_FORM_MO_I32", Nthreads, Nthreads_env



                     !!!!!!!!!!!!!!!!!!!!!!!!!!!
                     !!!!!      MAIN       !!!!!
                     !!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!! 3-2ERI ARRAY
      ALLOCATE(I32(NBASIS,NBASIS,MaxAUXBFX))

!!!!! BUFFERS
      ALLOCATE(T1(MaxT12))
      ALLOCATE(GG(MaxT12))
      ALLOCATE(TL(MXBSH*MXBSH))
      ALLOCATE(XX1(ISIZE))
      ALLOCATE(XX2(ISIZE*28))

!!!!! K-SHELL (AUX)
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

                        CALL CP3CAUXINT(GG,I32,TL,T1,                   &
                             NBASIS,LOCI,LOCJ,                          &
                             MINI,MAXI,MINJ,MAXJ,MINK,MAXK,             &
                             NI,NJ,NK,MINKS,MAXKS,                      &
                             NKS,FlgSphAux,IANDJ)

                  ENDDO
            ENDDO

!!!!! MO TRANSFORMED I32

            CALL RIMP2_FAT_MO_TRANSFORMED_I32                           &
                (I32,VEC,T1,                                            &
                 LOCKT,MINKT,MAXKT,NCOR,NACT,NVIR,NBASIS,NKT)

!!!!! DDI PUT I32

            DO LL=IauxStart,IauxEnd
               IBG=NVIR*(LL-1)+1
               IED=NVIR*LL
               CALL DCOPY(NVIR*NACT, I32(1,1,LL-IauxStart+1),1,B32(IBG:IED,1:NACT),1)
            ENDDO

      ENDDO

      DEALLOCATE(I32, T1,GG,TL,XX1,XX2)


      IF(NPROC.GT.1) THEN
         IBG=NVIR*(LddiAuxStart-1)+1
         IED=NVIR*LddiAuxEnd
         CALL DDI_PUT(D_B, IBG,IED, 1,NACT,B32(IBG:IED,1:NACT))
      ENDIF



      END !***********************************************************




      SUBROUTINE RIMP2_FAT_MO_TRANSFORMED_I32                           &
                 (I32,VEC,T1,LOCLT,MINLT,MAXLT,                         &
                  NCOR,NACT,NVIR,NBASIS,NKT)

      use Rimp2_Shared_Data
      use omp_lib

      implicit double precision(a-h,o-z)

      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      logical :: GOPARR,DSKWRK,MASWRK

!!!!! OUTPUT
      double precision :: I32(NBASIS,NBASIS,NKT)

!!!!! INPUT
      double precision :: VEC(NBASIS,NBASIS)
      double precision :: T1(NBASIS*NBASIS)



      DO L=MINLT,MAXLT
            LL=L-MINLT+1

            CALL DGEMM('N','N', NBASIS,NACT, NBASIS,                    &
                 1.0D00, I32(1,1,LL),NBASIS,                            &
                         VEC(1,1+NCOR),NBASIS,                          &
                 0.0D00, T1,NBASIS)

            CALL DGEMM('T','N',NVIR,NACT,NBASIS,                        &
                 1.0D00, VEC(1,1+NCOR+NACT),NBASIS,                     &
                         T1,NBASIS,                                     &
                 0.0D00, I32(1,1,LL),NVIR)
      ENDDO


      END !***********************************************************








!>*module rimp2grd   *deck RIMP2_COMBINE_VDX_I32
!>
!>     @brief   combine 3-2ERI with VDX
!>     @author  buu
!>     @date    sept 8, 2018
!>
!>     @detail  combine 3-2ERI with VDX
!>

      SUBROUTINE RIMP2_FAT_COMBINE_VDX_I32(VDX, NAUXBAS,NAUXBASD,NACT,NVIR,NBASIS)

      use Rimp2_Shared_Data
      use omp_lib

      implicit double PRECISION(A-H,O-Z)

      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK

! !!!!! OUTPUT
!       double precision :: B32_1(NAUXBAS*NVIR,NACT)

!!!!! INPUT
      double precision :: VDX(NAUXBASD,NAUXBAS)

!!!!! LOCAL
      double precision,allocatable,save,dimension(:,:) :: T2
!$omp threadprivate(T2)
      double precision,allocatable,dimension(:,:) :: T3,B32

      integer,allocatable :: ompISE(:,:)
      logical :: GOPARR,DSKWRK,MASWRK,LOMP




! !!!!! DDI WORK DISTRIBUTING FOR NACT LOOP
!       ALLOCATE(ompISE(0:NPROC-1,2))
!       CALL WorkSharedArray(ompISE,NPROC,1,NACT,LOMP)
!       LddiActStart=ompISE(ME,1)
!       LddiActEnd=ompISE(ME,2)
!       DEALLOCATE(ompISE)


!!!!! ddi work distribution based on
!!!!! local portion of DDI arrays
      CALL DDI_DISTRIB(D_B,ME,ILO,IHI,LddiActStart,LddiActEnd)

      ! CALL DDI_NDISTRIB(D_B,ME,ILO,IHI,JLO,JHI)
      ! WRITE(*,*) "xxxxx", ME,LddiActStart,LddiActEnd,LddiActEnd-LddiActStart+1,NACT



!!!!! omp threads
      Nthreads_env=omp_get_max_threads()
      Nthreads=MIN0(Nthreads_env,LddiActEnd-LddiActStart+1)
      IF(BUG) WRITE(*,'(A40,2I5)') "NUM_THREAD SUB RIMP2_FAT_COMBINE_VDX_I32", Nthreads, Nthreads_env


!!!!! buffers
      ALLOCATE(T3(NAUXBAS*NVIR,LddiActStart:LddiActEnd))


      ! CALL DDI_DISTRIB(D_B,ME,ILO,IHI,LddiActStart,LddiActEnd)

      ALLOCATE(B32(NAUXBAS*NVIR,LddiActStart:LddiActEnd)) !JLO:JHI))

      ! ! CALL DDI_NDISTRIB(D_B,ME,ILO,IHI,JLO,JHI)
      ! WRITE(*,*) "xxxxx", ME,JLO,JHI,JHI-JLO+1,NACT



! !!!!! RIP EXTRA RANKS
!       CALL DDI_SYNC(1456)
!       IF(ME.GT.NACT-1) GOTO 120


!!!!! GET B32
      ! IF(NPROC.GT.1) CALL DDI_GET(D_B, 1,NAUXBAS*NVIR, 1,NACT, B32)

      CALL DDI_GET(D_B, 1,NAUXBAS*NVIR, LddiActStart,LddiActEnd, B32)


!!!!! RIP EXTRA RANKS
      CALL DDI_SYNC(1456)
      IF(ME.GT.NACT-1) GOTO 120






!$OMP PARALLEL NUM_THREADS(NTHREADS)                                    &
!$omp default(none)                                                     &
!$omp shared(LddiActStart,LddiActEnd,                                   &
!$omp        NACT, NVIR,NAUXBAS,NAUXBASD,                               &
!$omp        B32,VDX)                                                   &
!$omp private(IACT)

      ALLOCATE(T2(NVIR,NAUXBAS))

!$omp do schedule(DYNAMIC)
      ! DO IACT=1,NACT
      DO IACT=LddiActStart,LddiActEnd

         CALL DCOPY(NAUXBAS*NVIR, B32(1,IACT),1, T2,1)

         CALL DGEMM('N','T', NAUXBASD,NVIR,NAUXBAS,                     &
              1.0D00, VDX,NAUXBASD,                                     &
                      T2,NVIR,                                          &
              0.0D00, B32(1:NAUXBASD*NVIR,IACT:IACT),NAUXBASD)
      ENDDO
!$omp end do

      DEALLOCATE(T2)

!$OMP END PARALLEL

  120 CONTINUE
      IF(NPROC.GT.1) THEN
         CALL DDI_PUT(D_B, 1,NAUXBASD*NVIR, LddiActStart,LddiActEnd,B32(1,LddiActStart))
      ENDIF

      DEALLOCATE(B32)


      END !***********************************************************









!>*module rimp2grd   *deck RIMP2_COMBINE_VDX_I32
!>
!>     @brief   combine 3-2ERI with VDX
!>     @author  buu
!>     @date    sept 8, 2018
!>
!>     @detail  combine 3-2ERI with VDX
!>

      SUBROUTINE RIMP2_1R_FAT_COMBINE_VDX_I32(B32, VDX, NAUXBAS,NAUXBASD,NACT,NVIR,NBASIS)

      use Rimp2_Shared_Data
      use omp_lib

      implicit double PRECISION(A-H,O-Z)

      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK

!!!!! OUTPUT
      double precision :: B32(NAUXBAS*NVIR,NACT)

!!!!! INPUT
      double precision :: VDX(NAUXBASD,NAUXBAS)

!!!!! LOCAL
      double precision,allocatable,save,dimension(:,:) :: T2
!$omp threadprivate(T2)

      integer,allocatable :: ompISE(:,:)
      logical :: GOPARR,DSKWRK,MASWRK,LOMP



!!!!! OMP THREADS
      Nthreads_env=omp_get_max_threads()
      Nthreads=MIN0(Nthreads_env,NACT)
      IF(BUG) WRITE(*,'(A40,2I5)') "NUM_THREAD SUB RIMP2_FAT_COMBINE_VDX_I32", Nthreads, Nthreads_env


!$OMP PARALLEL NUM_THREADS(NTHREADS)                                    &
!$omp default(none)                                                     &
!$omp shared(LddiActStart,LddiActEnd,                                   &
!$omp        NACT, NVIR,NAUXBAS,NAUXBASD,                               &
!$omp        B32,VDX)                                                   &
!$omp private(IACT)

      ALLOCATE(T2(NVIR,NAUXBAS))

!$omp do schedule(DYNAMIC)
      DO IACT=1,NACT
         CALL DCOPY(NAUXBAS*NVIR, B32(1,IACT),1, T2,1)
         CALL DGEMM('N','T', NAUXBASD,NVIR,NAUXBAS,                     &
              1.0D00, VDX,NAUXBASD,                                     &
                      T2,NVIR,                                          &
              0.0D00, B32(1,IACT),NAUXBASD)
      ENDDO
!$omp end do

      DEALLOCATE(T2)

!$OMP END PARALLEL




      END !***********************************************************









!>*module rimp2grd   *deck RIMP2_ENERGY
!>
!>     @brief   form density matrices
!>     @author  buu
!>     @date    sept 8, 2018
!>
!>     @detail  loop over active-active occ MOs
!>              to (partially) form density matrices
!>

      SUBROUTINE RIMP2_FAT_ENERGY                                       &
                (EIG,VEC, NANGM,NAUXBASD,NAUXBAS,                   &
                 NCOR,NACT,NOCC,NVIR,NORB,NBF)

      USE omp_lib
      USE Rimp2_Shared_Data,only: D_B,BUG,MXRT

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      COMMON /IOFILE/ IR,IW,IP,IJK,IPK,IDAF,NAV,IODA(950)
      COMMON /ENRGYS/ ENUCR,EELCT,ETOT,SZ,SZZ,ECORE,ESCF,EERD,E1,E2,    &
                      VEN,VEE,EPOT,EKIN,ESTATE(MXRT),STATN,EDFT(3),EDISP
      COMMON /ENRGMP/ EMP2,EMP3,EMP4,EMP2A
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /MACHIN/ NWDVAR,MAXFM,MAXSM,LIMFM,LIMSM


!!!!! INPUT
      double precision :: EIG(NBF)
      double precision :: VEC(NBF,NBF)

!!!!! LOCAL
      double precision,allocatable,dimension(:,:) :: eij,eab,B32_J,B32
      double precision,allocatable,dimension(:,:),save :: QVV,BI
!$omp threadprivate(QVV,BI)
      double precision,save :: E2_omp
!$omp threadprivate(E2_omp)

      integer,allocatable :: TRI(:,:)
      logical :: ULRange,LOMP
      logical :: GOPARR,DSKWRK,MASWRK


                     !!!!!!!!!!!!!!!!!!!!!!!!!!!
                     !!!!!      INIT       !!!!!
                     !!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!! parameters
      NOCC=NCOR+NACT
      NORB=NCOR+NACT+NVIR

      NDO=NAUXBASD*NOCC
      NDV=NAUXBASD*NVIR
      NDN=NAUXBASD*NORB

      E2=0.0D00

!!!!! make sure B32 has completely
!!!!! been put into the ddi array
      CALL DDI_SYNC(156)

!!!!! RIP EXTRA RANKS
      IF(ME.GT.NACT-1) THEN
         WRITE(*,*) "RANK SKIPPED", ME
         GOTO 120
      ENDIF


!!!!! occupied & virtual MO energy pairs
      ALLOCATE(eab(NVIR,NVIR))
      DO IB=1,NVIR
         DO IA=1,IB
            eab(IA,IB) = EIG(IA+NOCC) + EIG(IB+NOCC)
            eab(IB,IA) = eab(IA,IB)
         ENDDO
      ENDDO

      ALLOCATE(eij(NACT,NACT))
      DO JJ=1,NACT
         DO II=1,JJ
            eij(II,JJ)=EIG(II+NCOR) + EIG(JJ+NCOR)
         ENDDO
      ENDDO


!!!!! TURN OFF DYNAMIC THREADS TO
!!!!! KEEP threadprivate DATA CONSISTENT
      CALL OMP_SET_DYNAMIC(.FALSE.)

!!!!! omp threads
      Nthreads=omp_get_max_threads()
      IF(BUG) WRITE(*,'(A40,2I5)') "NUM_THREAD SUB RIMP2_FAT_ENERGY", Nthreads

!!!!! trapezoidal decomposition act-act occ pairs
      IF(NACT.LT.NPROC) THEN
         NSLICES = NACT
      ELSE
         NSLICES = NPROC
      ENDIF
      CALL RIMP2_TRAPE_DEC(LddiActStart,LddiActEnd,ME,NSLICES,NACT)
      ALLOCATE(B32_J(NAUXBASD*NVIR,LddiActStart:LddiActEnd))


!!!!! tmp fusing array TRI
      ALLOCATE(TRI(NACT*(NACT+1)/2,2))
      LddiActStart_t=LddiActStart
      IF(MASWRK) LddiActStart_t=LddiActStart+1
      ISIZE=0
      DO JACT=LddiActStart_t,LddiActEnd
         DO IACT=1,JACT-1
            ISIZE=ISIZE+1
            TRI(ISIZE,1)=JACT
            TRI(ISIZE,2)=IACT
         ENDDO
      ENDDO



!!!!! check available node memory for B32
      IERR = 1
      NSPLIT=0
      DO II=1,NACT
         ALLOCATE(B32(NAUXBASD*NVIR,NACT-NSPLIT),STAT=IERR)
         IF(IERR.GT.0) THEN
            NSPLIT=NSPLIT+10
            DEALLOCATE(B32)
         ELSE
            ! if the whole B32 fits process mem goto
            ! label 300 and call RIMP2_ENERGY_FUSED
            IF(NSPLIT.EQ.0) THEN

      ! w21=0.0D00
      ! w1=omp_get_wtime()

               ! WRITE(*,'(A34,I4)') "THE WHOLE B32 IS COPIED TO PROC", ME
               DEALLOCATE(B32_J)
               CALL DDI_GET(D_B, 1,NAUXBASD*NVIR, 1,NACT-NSPLIT, B32)

      ! w2=omp_get_wtime()

      ! w21=w21 + w2-w1
      ! WRITE(*,'(A10,F5.1,5I5)') "hhhhhhh", w21,ME,LddiActStart,LddiActEnd,LddiActEnd-LddiActStart+1,ISIZE

               GOTO 300
            ELSE
               WRITE(*,'(F5.2,A38,I4)') real(NACT-NSPLIT)/real(NACT), "OF B32 WOULD BE COPIED INTO RANK", ME
               CALL DDI_GET(D_B, 1,NAUXBAS*NVIR, LddiActStart,LddiActEnd, B32_J)
               EXIT
            ENDIF
         ENDIF
      ENDDO



        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!! ONLY PART of B32 loaded on node mem !!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        WRITE(*,*) "ppppppppp", ME

!$OMP PARALLEL NUM_THREADS(Nthreads)                                    &
!$omp default(none)                                                     &
!$omp shared(LddiActStart,LddiActEnd, Istart,Iend,                      &
!$omp        NACT,NVIR,NAUXBASD,                                        &
!$omp        B32,B32_J,eij,eab,E2,ME,MASWRK,D_B,NSPLIT,ISIZE,TRI)       &
!$omp private(KK,IACT,JACT) &
!$omp private(w43,w44,w45,w46,IthreadID,LOMP)

      ALLOCATE(BI(NAUXBASD,NVIR))
      ALLOCATE(QVV(NVIR,NVIR))

      E2_omp = 0.0D00

!!!!! IACT .ne. JACT
!$omp do schedule(DYNAMIC)
      DO KK=1,ISIZE

         JACT=TRI(KK,1)
         IACT=TRI(KK,2)

         IF(JACT-1.LE.NACT-NSPLIT) THEN
            CALL DCOPY(NVIR*NAUXBASD,B32(1,IACT),1,BI,1)
         ELSEIF(IACT.GE.LddiActStart .AND. IACT.LE.LddiActEnd) THEN
            CALL DCOPY(NVIR*NAUXBASD,B32_J(1,IACT),1,BI,1)
         ELSE
            CALL DDI_GET(D_B, 1,NAUXBASD*NVIR, IACT,IACT, BI)
         ENDIF

         CALL RIMP2_ENERGYIJ                                            &
            (E2_omp, BI,B32_J(1,JACT),                                  &
             eij(IACT,JACT),eab, QVV,                                   &
             IACT,JACT,2.0D00,NVIR,NAUXBASD)

      ENDDO !KK
!$omp end do nowait


!!!!! IACT .eq. JACT
!$omp do schedule(DYNAMIC)
      DO JACT=LddiActStart,LddiActEnd
         CALL RIMP2_ENERGYIJ                                            &
             (E2_omp, B32_J(1,JACT),B32_J(1,JACT),                      &
              eij(JACT,JACT),eab, QVV,                                  &
              JACT,JACT,1.0D00,NVIR,NAUXBASD)
      ENDDO !JACT
!$omp end do nowait

!!!!! reduce E2
!$omp atomic
      E2 = E2 + E2_omp

      DEALLOCATE(BI)
      DEALLOCATE(QVV)

!$OMP END PARALLEL

      GOTO 450



                  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  !!!!! whole B32 loaded on node mem !!!!!
                  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  300 CONTINUE
      ! w43=0.0D00
      ! w3=omp_get_wtime()

      CALL RIMP2_ENERGY_FUSED &
          (E2,TRI,ISIZE,B32, eij,eab,LddiActStart,LddiActEnd,NAUXBASD,NACT,NVIR)

      w4=omp_get_wtime()
      w43=w43 + w4-w3
      ! WRITE(*,'(A10,F5.1,I5)') "oooooooooo", w43,ME

                       !!!!!!!!!!!!!!!!!!!!!!!!!!!
                       !!!!!      CLOSE      !!!!!
                       !!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!! skip RIMP2_ENERGY_FUSED
  450 CONTINUE


!!!!! end point
!!!!! extra ranks jump down here
  120 CONTINUE


!!!!! REDUCTION:: CORRELATION ENERGY
      CALL DDI_GSUMF(5999,E2,1)


!!!!! RHF RI-MP2 ENERGY
      EMP2 = E2+ESCF


!!!!! PRINT OUT MP2 ENERGY
      IF(MASWRK) WRITE(IW,100) ESCF,E2,EMP2
  100 FORMAT (/1X,'RHF RI MP2 ENERGY',/,                                &
             12X,'   ESCF=',1X,F20.10/                                  &
             12X,'   E(2)=',1X,F20.10/                                  &
             12X,' E(MP2)=',1X,F20.10 ,/)



      END !***********************************************************





      SUBROUTINE RIMP2_ENERGY_FUSED(E2,TRI,ISIZE, B32, eij,eab,LddiActStart,LddiActEnd,NAUXBASD,NACT,NVIR)

      use omp_lib
      use Rimp2_Shared_Data,only : BUG,D_B

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK

!!!!! OUTPUT
      double precision :: E2

!!!!! INPUT
      integer :: TRI(NACT*(NACT+1)/2,2)
      double precision :: eij(NACT,NACT)
      double precision :: eab(NVIR,NVIR)
      double precision :: B32(NAUXBASD*NVIR,NACT)

!!!!! LOCAL DATA
      logical :: GOPARR,DSKWRK,MASWRK
      double precision,save :: E2_omp
!$omp threadprivate(E2_omp)
      double precision,allocatable,save :: QVV(:,:)
!$omp threadprivate(QVV)





!!!!! TURN OFF DYNAMIC THREADS TO
!!!!! KEEP threadprivate DATA CONSISTENT
      CALL OMP_SET_DYNAMIC(.FALSE.)


!!!!! omp threads
      Nthreads=omp_get_max_threads()
      IF(BUG) WRITE(*,'(A40,2I5)') "NUM_THREAD SUB RIMP2_ENERGY_FUSED", Nthreads


!$OMP PARALLEL NUM_THREADS(Nthreads)                                    &
!$omp default(none)                                                     &
!$omp shared(LddiActStart,LddiActEnd, Istart,Iend,                      &
!$omp        NACT,NVIR,NAUXBASD,                                        &
!$omp        B32,eij,eab,E2,TRI,ISIZE,MASWRK)                           &
!$omp private(IACT,JACT,w43,w44)

      E2_omp = 0.0D00

      ALLOCATE(QVV(NVIR,NVIR))

!!!!! IACT .ne. JACT
!$omp do schedule(DYNAMIC)
      DO KK=1,ISIZE

         JACT=TRI(KK,1)
         IACT=TRI(KK,2)

         CALL RIMP2_ENERGYIJ                                            &
             (E2_omp, B32(1,IACT),B32(1,JACT),                          &
              eij(IACT,JACT),eab, QVV,                                  &
              IACT,JACT,2.0D00,NVIR,NAUXBASD)

      ENDDO !KK
!$omp end do nowait


!!!!! IACT .eq. JACT
!$omp do schedule(DYNAMIC)
      DO JACT=LddiActStart,LddiActEnd
         CALL RIMP2_ENERGYIJ                                            &
             (E2_omp, B32(1,JACT),B32(1,JACT),                          &
              eij(JACT,JACT),eab, QVV,                                  &
              JACT,JACT,1.0D00,NVIR,NAUXBASD)
      ENDDO !JACT
!$omp end do nowait

!$omp atomic
      E2 = E2 + E2_omp

      DEALLOCATE(QVV)

!$OMP END PARALLEL

      END !*************************************************************




      ! SUBROUTINE RIMP2_TRAPE_DEC(Istart,Iend,N)


      ! COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK

      ! logical :: GOPARR,DSKWRK,MASWRK

      ! IF(N .LE. NPROC) THEN
      !    Istart = ME + 1
      !    Iend = ME + 1
      !    RETURN
      ! ENDIF

      ! NTOT = (N*(N+1)) / NPROC
      ! Istart = 1
      ! DO II = 0,ME
      !    TMP = NTOT + Istart*Istart - Istart
      !    End = (SQRT(4*TMP + 1.0) - 1.0)/2.0
      !    Iend = nint(End)
      !    IF(II.LT.ME) Istart = Iend + 1
      ! ENDDO

      ! END !*************************************************************

      SUBROUTINE RIMP2_TRAPE_DEC(Istart,Iend,ME,NSlices,N)

      implicit none

      INTEGER,INTENT(IN) :: N ! outer loop bound of the triangular loop
      INTEGER,INTENT(IN) :: NSlices ! slicing triangular loop into Nslices
      INTEGER,INTENT(IN) :: ME ! the ME-th slice
      INTEGER,INTENT(OUT) :: Istart,Iend ! start and end of ME-th slice

      ! some local variables
      INTEGER,allocatable :: WSA(:,:)
      INTEGER :: KK,I,J,NITER,NRES
      INTEGER :: NChunksize
      INTEGER :: NChunks
      INTEGER :: II

      ! not satisfying work distribution

      ALLOCATE(WSA(0:NSlices-1, 2))

      ! when the outer loop bound is smaller than nproc
      IF (N .LE. NSlices) THEN
         DO II = 0, N-1
            WSA(II,1) = II + 1
            WSA(II,2) = II + 1
         ENDDO
         Istart = WSA(ME,1)
         Iend = WSA(ME,2)
         NChunks = N
         RETURN
      ENDIF

      ! when the outer loop bound is greater than nproc
      NChunksize = (N*(N+1)) / (2*NSlices)

      KK = 0
      WSA(0,1) = 1
      DO J = 1, N
         NITER = J*(J-1)/2
         DO I = 1, J
            NITER = NITER + 1

            IF (MOD(NITER, NChunksize) .EQ. 0) THEN
               NRES =  (NSlices-KK) - (N-J+1)
               IF(NRES .GE. 0) THEN
                  WSA(KK,2) = J - NRES
                  GOTO 123
               ELSE
                  WSA(KK,2) = J
               ENDIF

               KK = KK + 1

               WSA(KK,1) = J + 1

            ENDIF
         ENDDO
      ENDDO
      WSA(NSlices-1,2) = N
      GOTO 456


123 CONTINUE
      DO J = WSA(KK,2) + 1, N
         KK = KK + 1
         WSA(KK,1) = J
         WSA(KK,2) = J
      ENDDO

456 CONTINUE
      NChunks = NSlices
      Istart = WSA(ME,1)
      Iend = WSA(ME,2)

      END ! ***************************************************************



      SUBROUTINE RIMP2_ENERGYIJ                                         &
                (E2, BI,BJ,eij,eab,                                     &
                 QVV, IACT,JACT,FAC,NVIR,NAUXBASD)

      implicit double precision(a-h,o-z)

!!!!! OUTPUT
      double precision :: E2

!!!!! INPUT
      double precision :: BI(NAUXBASD*NVIR)
      double precision :: BJ(NAUXBASD*NVIR)
      double precision :: eab(NVIR,NVIR)
      double precision :: eij

!!!!! BUFFER
      double precision :: QVV(NVIR,NVIR)

      double precision,allocatable,dimension(:,:) :: QVV_t


      ALLOCATE(QVV_t(NVIR,NVIR))


      ! CALL DGEMM('T','N', NVIR,NVIR,NAUXBASD,                           &
      !             1.0D+00, BI,NAUXBASD,                                 &
      !                      BJ,NAUXBASD,                                 &
      !             0.0D+00, QVV_t,NVIR)
      ! E2_t = 0.0D00
      ! DO IB=1,NVIR
      !    DO IA=1,NVIR
      !       T2 = QVV_t(IA,IB)/(eij-eab(IA,IB))
      !       TMP = QVV_t(IA,IB)+QVV_t(IA,IB)
      !       TMP = TMP-QVV_t(IB,IA)
      !       E2_t = E2_t + T2*TMP
      !    ENDDO
      ! ENDDO
      ! E2 = E2 + FAC*E2_t



      CALL DGEMM('T','N', NVIR,NVIR,NAUXBASD,                           &
                  1.0D00, BI,NAUXBASD,                                  &
                          BJ,NAUXBASD,                                  &
                  0.0D00, QVV_t,NVIR)

      QVV = transpose(QVV_t)
      QVV = QVV_t + QVV_t - QVV
      QVV_t = QVV_t/(eij - eab)
      E2_t = DDOT(NVIR*NVIR,QVV,1,QVV_t,1)
      E2 = E2 + FAC*E2_t


      END








!>*module rimp2grd   *deck RIMP2_ENERGY
!>
!>     @brief   form density matrices
!>     @author  buu
!>     @date    sept 8, 2018
!>
!>     @detail  loop over active-active occ MOs
!>              to (partially) form density matrices
!>

      SUBROUTINE RIMP2_1R_FAT_ENERGY                                    &
                (B32,EIG,VEC, NANGM,NAUXBASD,NAUXBAS,                   &
                 NCOR,NACT,NOCC,NVIR,NORB,NBF)

      USE omp_lib
      USE Rimp2_Shared_Data,only: D_B,BUG,MXRT

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      COMMON /IOFILE/ IR,IW,IP,IJK,IPK,IDAF,NAV,IODA(950)
      COMMON /ENRGYS/ ENUCR,EELCT,ETOT,SZ,SZZ,ECORE,ESCF,EERD,E1,E2,    &
                      VEN,VEE,EPOT,EKIN,ESTATE(MXRT),STATN,EDFT(3),EDISP
      COMMON /ENRGMP/ EMP2,EMP3,EMP4,EMP2A
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /MACHIN/ NWDVAR,MAXFM,MAXSM,LIMFM,LIMSM


!!!!! INPUT
      double precision :: B32(NAUXBASD*NVIR,NACT)
      double precision :: VDX(NAUXBASD,NAUXBAS)
      double precision :: EIG(NBF)
      double precision :: VEC(NBF,NBF)

!!!!! LOCAL
      double precision,allocatable,dimension(:,:) :: eij,eab

      integer,allocatable :: ddiWSA(:,:),I_DLB_INDEX(:),TRI(:,:)
      logical :: ULRange,LOMP
      logical :: GOPARR,DSKWRK,MASWRK




                     !!!!!!!!!!!!!!!!!!!!!!!!!!!
                     !!!!!      INIT       !!!!!
                     !!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!! TURN OFF DYNAMIC THREADS TO KEEP threadprivate DATA CONSISTENT
      CALL OMP_SET_DYNAMIC(.FALSE.)

      ! parameters
      NOCC=NCOR+NACT
      NORB=NCOR+NACT+NVIR

      NDO=NAUXBASD*NOCC
      NDV=NAUXBASD*NVIR
      NDN=NAUXBASD*NORB


!!!!! OMP THREADS
      Nthreads_env=omp_get_max_threads()
      Nthreads=Nthreads_env
      IF(BUG) WRITE(*,'(A40,2I5)') "NUM_THREAD SUB RIMP2_1R_FAT_ENERGY", Nthreads,Nthreads_env

!!!!! INIT
      E2=0.0D00

!!!!! occupied & virtual MO energy pairs
      ALLOCATE(eab(NVIR,NVIR))
      DO IB=1,NVIR
         DO IA=1,IB
            eab(IA,IB) = EIG(IA+NOCC) + EIG(IB+NOCC)
            eab(IB,IA) = eab(IA,IB)
         ENDDO
      ENDDO

      ALLOCATE(eij(NACT,NACT))
      DO JJ=1,NACT
         DO II=1,JJ
            eij(II,JJ)=EIG(II+NCOR) + EIG(JJ+NCOR)
         ENDDO
      ENDDO


      ALLOCATE(TRI(NACT*(NACT+1)/2,2))
      ISIZE=0
      DO JACT=2,NACT
         DO IACT=1,JACT-1
            ISIZE=ISIZE+1
            TRI(ISIZE,1)=JACT
            TRI(ISIZE,2)=IACT
         ENDDO
      ENDDO


                     !!!!!!!!!!!!!!!!!!!!!!!!!!!
                     !!!!!      MAIN       !!!!!
                     !!!!!!!!!!!!!!!!!!!!!!!!!!!


      CALL RIMP2_ENERGY_FUSED(E2,TRI,ISIZE, B32, eij,eab, &
        1,NACT,NAUXBASD,NACT,NVIR)


                     !!!!!!!!!!!!!!!!!!!!!!!!!!!
                     !!!!!      CLOSE      !!!!!
                     !!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! DEALLOCATE(TRI)

!!!!! RHF RI-MP2 ENERGY
      EMP2 = E2+ESCF



!!!!! PRINT OUT MP2 ENERGY
      IF(MASWRK) WRITE(IW,100) ESCF,E2,EMP2
  100 FORMAT (/1X,'RHF RI MP2 ENERGY',/,                                &
             12X,'   ESCF=',1X,F20.10/                                  &
             12X,'   E(2)=',1X,F20.10/                                  &
             12X,' E(MP2)=',1X,F20.10 ,/)



      END !***********************************************************






! !>*module rimp2grd   *deck RIMP2_ENERGY
! !>
! !>     @brief   form density matrices
! !>     @author  buu
! !>     @date    sept 8, 2018
! !>
! !>     @detail  loop over active-active occ MOs
! !>              to (partially) form density matrices
! !>

!       SUBROUTINE RIMP2_FAT_ENERGY___reserved                            &
!                 (B32,EIG,VEC, NANGM,NAUXBASD,NAUXBAS,                   &
!                  NCOR,NACT,NOCC,NVIR,NORB,NBF)

!       USE omp_lib
!       USE Rimp2_Shared_Data,only: D_B,BUG,MXRT

!       IMPLICIT DOUBLE PRECISION(A-H,O-Z)

!       COMMON /IOFILE/ IR,IW,IP,IJK,IPK,IDAF,NAV,IODA(950)
!       COMMON /ENRGYS/ ENUCR,EELCT,ETOT,SZ,SZZ,ECORE,ESCF,EERD,E1,E2,    &
!                       VEN,VEE,EPOT,EKIN,ESTATE(MXRT),STATN,EDFT(3),EDISP
!       COMMON /ENRGMP/ EMP2,EMP3,EMP4,EMP2A
!       COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
!       COMMON /MACHIN/ NWDVAR,MAXFM,MAXSM,LIMFM,LIMSM


! !!!!! INPUT
!       double precision :: B32(NAUXBASD*NVIR,NACT)
!       double precision :: VDX(NAUXBASD,NAUXBAS)
!       double precision :: EIG(NBF)
!       double precision :: VEC(NBF,NBF)

! !!!!! LOCAL
!       double precision,allocatable,dimension(:,:) :: YDV,eij,eab,QVV,TVV

!       integer,allocatable :: ddiWSA(:,:),I_DLB_INDEX(:)
!       logical :: MASWRK,ULRange,LOMP


!                      !!!!!!!!!!!!!!!!!!!!!!!!!!!
!                      !!!!!      INIT       !!!!!
!                      !!!!!!!!!!!!!!!!!!!!!!!!!!!


! !!!!! TURN OFF DYNAMIC THREADS TO KEEP threadprivate DATA CONSISTENT
!       CALL OMP_SET_DYNAMIC(.FALSE.)

!       ! parameters
!       NOCC=NCOR+NACT
!       NORB=NCOR+NACT+NVIR

!       NDO=NAUXBASD*NOCC
!       NDV=NAUXBASD*NVIR
!       NDN=NAUXBASD*NORB


! !!!!! OMP THREADS
!       Nthreads_env=omp_get_max_threads()
!       Nthreads=Nthreads_env
!       IF(BUG) WRITE(*,'(A40,2I5)') "NUM_THREAD SUB RIMP2_ENERGY", Nthreads,Nthreads_env

! !!!!! INIT
!       E2=0.0D00

! !!!!! occupied MO energy pair
!       ALLOCATE(eab(NVIR,NVIR))
!       FORALL(IA=1:NVIR,IB=1:NVIR)
!          eab(IA,IB) = EIG(IA+NOCC) + EIG(IB+NOCC)
!       ENDFORALL

! !!!!! virtual MO energy pair
!       ALLOCATE(eij(NACT,NACT))
!       FORALL(II=1:NACT,JJ=1:NACT)
!          eij(II,JJ) = EIG(II+NCOR) + EIG(JJ+NCOR)
!       ENDFORALL

! !!!!! ddi local data chunk
!       ALLOCATE(ddiWSA(0:NPROC-1,2))
!       CALL rimp2_WSA(ddiWSA,NPROC,1,NACT,LOMP)
!       LddiActStart=ddiWSA(ME,1)
!       LddiActEnd=ddiWSA(ME,2)


! !!!!! BUFFERS
!       ALLOCATE(YDV(NAUXBASD,NVIR))
!       ALLOCATE(QVV(NVIR,NVIR))


!                      !!!!!!!!!!!!!!!!!!!!!!!!!!!
!                      !!!!!      MAIN       !!!!!
!                      !!!!!!!!!!!!!!!!!!!!!!!!!!!




!       startw_static = omp_get_wtime()

!       ! DO JACT=1,NACT
!       DO JACT=LddiActStart,LddiActEnd !- N_DLB_BUFF
!             JJ=JACT+NCOR

!             YDV=0.0D00

!             DO IACT=1, NACT
!                CALL DGEMM('T','N', NVIR,NVIR,NAUXBASD,                  &
!                            1.0D+00, B32(1,IACT),NAUXBASD,               &
!                                     B32(1,JACT),NAUXBASD,               &
!                            0.0D+00, QVV,NVIR)

!                FORALL(IA=1:NVIR,IB=1:NVIR)
!                   QVV(IA,IB)=QVV(IA,IB)/(eij(IACT,JACT)-eab(IA,IB))
!                ENDFORALL

!                TVV = QVV + QVV
!                TVV = TVV - transpose(QVV)

! !!!!! YDV <- BJ x transpose(TVV)
!                CALL DGEMM('N','N', NAUXBASD,NVIR,NVIR,                  &
!                            1.0D00, B32(1,IACT),NAUXBASD,                &
!                                    TVV,NVIR,                            &
!                            1.0D00, YDV,NAUXBASD)
!             ENDDO

! !!!!! CORRELATION ENERGY
!             E2=E2+DDOT(NAUXBASD*NVIR, B32(1,JACT),1, YDV,1)

!       ENDDO !JACT
!       endw_static = omp_get_wtime()
!       IF(BUG) WRITE(*,'(A40,I5,F10.1)') "TIME ACT_ACT_STATIC: ",ME, endw_static-startw_static


!                      !!!!!!!!!!!!!!!!!!!!!!!!!!!
!                      !!!!!      CLOSE      !!!!!
!                      !!!!!!!!!!!!!!!!!!!!!!!!!!!


! !!!!! REDUCTION:: CORRELATION ENERGY
!       CALL DDI_GSUMF(5999,E2,1)

! !!!!! RHF RI-MP2 ENERGY
!       EMP2 = E2+ESCF

! !!!!! PRINT OUT MP2 ENERGY
!       IF(MASWRK) WRITE(IW,100) ESCF,E2,EMP2
!  100  FORMAT(/1X,'RHF RI MP2 ENERGY',/,                                 &
!              12X,'   ESCF=',1X,F20.10/                                  &
!              12X,'   E(2)=',1X,F20.10/                                  &
!              12X,' E(MP2)=',1X,F20.10 ,/)

!       endw = omp_get_wtime()
!       IF(BUG) WRITE(*,'(A40,I5,F10.1)') "TIME ACT_ACT: OMP REDUCTION",ME, endw-startw



!       END !***********************************************************



      SUBROUTINE RIMP2_1R_X_VXX(VXX, NAUXBAS) !1R = 1 RANK
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
      double precision,allocatable,dimension(:) :: TI,TJ,W,GHONDO

      integer,allocatable,dimension(:,:) :: ompISE
      logical :: GOPARR,DSKWRK,MASWRK, IANDJ,LOMP
      double precision :: CC(3,2)


                          !!! ~ GO ~ !!!

      VXX=0.0D00

!!!!! MAX BASIS FUNCTIONS
      MaxG=MaxAUXBFX*MaxAUXBFX

      ALLOCATE(TI(MaxG))
      ALLOCATE(TJ(MaxG))
      ALLOCATE(W(MaxG))
      ALLOCATE(GHONDO(MaxG))


      ! evaluate VXX by 1 process

      Nthreads=omp_get_max_threads()

!!!!! FORM VXX
!$omp parallel NUM_THREADS(Nthreads)                                    &
!$omp shared(VXX, ompISE,                                               &
!$omp        KAUXAT,C,KAUXTY,KAUXNG,KAUXMI,KAUXMX,KAUXST,KAUXLO,        &
!$omp        EXAUX,CAUXI,CAUXH,CAUXG,CAUXF,CAUXD,CAUXP,CAUXS,NORM,      &
!$omp        FlgSphAux,LOCSPH,NAUXBAS)                                  &
!$omp private(GHONDO,W, IthreadID,II,Istart,Iend,                       &
!$omp         IAT,CC,LANGI,NGI,MINI,MAXI,I1,LOCI,NI,                    &
!$omp         MINIS,MAXIS,LOCIS,NIS,                                    &
!$omp         JJ,JAT,LANGJ,NGJ,MINJ,MAXJ,J1,LOCJ,NJ,IANDJ,TI,TJ,        &
!$omp         MINJS,MAXJS,LOCJS,NJS)

!$omp do schedule(DYNAMIC)
      DO II = 1, NAUXSH
      ! DO II=LddiAuxShStart,LddiAuxShEnd
      ! DO II=LompAuxShStart,LompAuxShEnd
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

!$omp end do

!$OMP END PARALLEL

  ! 120 CONTINUE
  !     CALL DDI_GSUMF(1504,VXX,NAUXBAS*NAUXBAS)

! !!!!! CLEAN UP
!       ! DEALLOCATE(TI,TJ,W,GHONDO,ompISE)
!       DEALLOCATE(TI)
!       DEALLOCATE(TJ)
!       DEALLOCATE(W)
!       DEALLOCATE(GHONDO)



      IF(BUG) WRITE(*,*) "wwww exiting RIMP2_X_VXX", ME


      RETURN

      END !***********************************************************





      SUBROUTINE RIMP2_AUXSHAUX &
                 (LddiAuxShStart,LddiAuxShEnd, LddiAuxStart,LddiAuxEnd, NAUXBAS)

      use Rimp2_Shared_Data

      implicit double precision(a-h,o-z)

      COMMON /AUXBAS/ EXAUX(MXAXGTOT),CAUXS(MXAXGTOT),CAUXP(MXAXGTOT),  &
                      CAUXD(MXAXGTOT),CAUXF(MXAXGTOT),CAUXG(MXAXGTOT),  &
                      CAUXH(MXAXGTOT),CAUXI(MXAXGTOT),                  &
                      KAUXST(MXAUXSH),KAUXAT(MXAUXSH),KAUXTY(MXAUXSH),  &
                      KAUXNG(MXAUXSH),KAUXLO(MXAUXSH),KAUXMI(MXAUXSH),  &
                      KAUXMX(MXAUXSH),NAUXSH
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK

      double precision,allocatable,dimension(:) :: TL
      integer,allocatable,dimension(:) :: NSIZE
      integer,allocatable,dimension(:,:) :: ISHELL

      logical :: GOPARR,DSKWRK,MASWRK


      ! divide AUX shells based on AUX

                     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                     !!!!!  ARRAY OF AUX SHELLS !!!!!!
                     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! num aux basis fnx per shell
      ! stored in NSIZE(NAUXSH)
      MXBSH=MaxAUXBFX
      IF(MaxAUXANG < MaxATMANG) MXBSH=MaxATMBFX
      ALLOCATE(TL(MXBSH*MXBSH))
      ALLOCATE(NSIZE(NAUXSH))
      DO KAUXSH=1,NAUXSH
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
         NSIZE(KAUXSH)=MAXKT-MINKT+1
      ENDDO !KAUXSH


                     !!!!!!!!!!!!!!!!!!!
                     !!!!!  INIT  !!!!!!
                     !!!!!!!!!!!!!!!!!!!

      ! shells
      LddiAuxShStart=0
      LddiAuxShEnd=0

      ! basis functions
      LddiAuxStart=0
      LddiAuxEnd=0


                     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                     !!!!!  NPROCS .GE. NAUXSH  !!!!!!
                     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! if num of ddi procs larger or equal to aux shells
      IF(NPROC.GE.NAUXSH) THEN
         IF(ME.GE.NAUXSH) RETURN
         ! shell
         LddiAuxShStart=ME+1
         LddiAuxShEnd=ME+1
         ! basis functions
         LddiAuxStart = SUM(NSIZE(1:LddiAuxShStart)) - NSIZE(LddiAuxShStart) + 1
         LddiAuxEnd = SUM(NSIZE(1:LddiAuxShEnd))
         RETURN
      ENDIF


                     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                     !!!!!!  NPROC .LT. NAUXSH  !!!!!!
                     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      ! estimiate num aux basis func per process
      chunks = real(NAUXBAS)/real(NPROC)


      ALLOCATE(ISHELL(0:NPROC-1,2))
      ! ALLOCATE(IBAFNX(0:NPROC-1,2))

      ISHELL=99999
      ISHELL(0,1) = 1


      Nchunks = 0

      DO II=0,NPROC-1
         Nchunks = Nchunks + chunks
         NSUM=0
         DO KK=1,NAUXSH
            NSUM=NSUM+NSIZE(KK)
            IF(NSUM.GE.Nchunks) THEN
               IF(NSUM-Nchunks .LE. Nchunks+NSIZE(KK)-NSUM) THEN
                  ISHELL(II,2) = KK
               ELSE
                  ISHELL(II,2) = KK-1
               ENDIF
               EXIT
            ENDIF
         ENDDO

         IF(ISHELL(II,1).GE.NAUXSH) THEN
            ISHELL(II,2)=ISHELL(II,1)
            EXIT
         ENDIF

         IF(ISHELL(II,1).GE.ISHELL(II,2)) ISHELL(II,2)=ISHELL(II,1)

         IF(II.LT.NPROC-1) ISHELL(II+1,1)=ISHELL(II,2)+1

      ENDDO

      IF(II.LT.NPROC-1) THEN
         LL=0
         DO KK=II,0,-1
            LL=LL+ISHELL(KK,2)-ISHELL(KK,1)
            IF(LL.GE.NPROC-1-II) EXIT
         ENDDO

         NEXTRA=LL-(NPROC-1-II)
         ISHELL(KK,2)=ISHELL(KK,1)+NEXTRA ! (LL-NPROC+II)

         DO MM=KK+1,NPROC-1
            ISHELL(MM,1) = ISHELL(MM-1,2)+1
            ISHELL(MM,2) = ISHELL(MM,1)
         ENDDO

         GOTO 120
      ENDIF

      NEXTRA=NAUXSH-ISHELL(NPROC-1,2) + 1

      IF(NEXTRA.GT.0) THEN
         KK=0
         DO II=NPROC-NEXTRA+1,NPROC-1
            KK=KK+1
            ISHELL(II,2)=ISHELL(II,2)+KK
            IF(II.LT.NPROC-1) ISHELL(II+1,1)=ISHELL(II,2)+1
         ENDDO
      ENDIF


  120 CONTINUE
      LddiAuxShStart=ISHELL(ME,1)
      LddiAuxShEnd=ISHELL(ME,2)

      LddiAuxStart=SUM(NSIZE(1:LddiAuxShStart))-NSIZE(LddiAuxShStart)+1
      LddiAuxEnd=SUM(NSIZE(1:LddiAuxShEnd))


      ! WRITE(*,'(A9,4I5,A5,3I5,F7.2,4I5)') "vvvvvv333", &
      !      LddiAuxShStart,LddiAuxShEnd,LddiAuxShEnd-LddiAuxShStart+1,NAUXSH, &
      !      "##33", &
      !      LddiAuxStart,LddiAuxEnd,LddiAuxEnd-LddiAuxStart+1,chunks, &
      !      NAUXBAS,SUM(NSIZE), NPROC,ME


      ! DEALLOCATE(TL,NSIZE,ISHELL)


      END



!*MODULE RIMP2OMP   *DECK WorkSharedArray
!>
!>     @BRIEF   WORK SHARING AMONG NPROCS PROCESSES/THREADS
!>     @AUTHOR  BUU PHAM
!>     @DATE    SEP 25, 2017
!>
!>     @DETAIL  A PORTION OF LOOP FROM Lstart TO Lend IS DIVIDED
!>              AMONG Nprocs (MPI PROCESSES OR THREADS). THE WORK
!>              SHRING IS STORED IN THE ISE ARRAY
!>
      SUBROUTINE WorkSharedArray(ISE,Nprocs,Lstart,Lend,LOMP)

      IMPLICIT DOUBLE precision(a-h,o-z)

      INTEGER,intent(in)  :: Nprocs
      INTEGER,intent(in)  :: Lstart, Lend

      LOGICAL,intent(out) :: LOMP
      INTEGER,intent(out) :: ISE(0:Nprocs-1,2)


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


      END !**********************************************************


#endif
