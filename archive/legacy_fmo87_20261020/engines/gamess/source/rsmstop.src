C*MODULE RSM_STOP  *DECK RSMCC_MOD
C>
C>    @brief   module to store the hole and particle orbital energy array
C>             detemined in gas phase
C>
C>    @author  Daisuke Yokogawa
C>
      MODULE RSMCC_MOD
      IMPLICIT NONE

      DOUBLE PRECISION, ALLOCATABLE :: OEH_ORG(:), OEP_ORG(:)

      END MODULE

C*MODULE RSM_STOP
C>
C>    @brief   module for RISM I/O
C>
C>    @author  Daisuke Yokogawa
C>
      MODULE RISM_IO
      IMPLICIT NONE

      INTEGER, PARAMETER :: IRSDAF=26, IRSSQF=27
      INTEGER, PARAMETER :: MAXIO = 1000, IRECLN = 4090
C
C     RISM INPUT, OUTPUT, and PUNCH FILES (GLOBAL)
C
      INTEGER            :: IR, IW, IP
C
C     RISM DICTIONARY FILE (GLOBAL)
C
      INTEGER            :: IRECST,IFILEN(MAXIO),IODARS(MAXIO)
C
      CONTAINS

C*MODULE RSM_STOP  *DECK RSM_SETIO
C>
C>    @brief   read i/o data from dictionary file
C>
C>    @author  Daisuke Yokogawa
C>
        SUBROUTINE RSM_SETIO

          INTEGER :: I0
C
C         IR, IW, IP
C
          CALL RSM_DARE(I0,1,501,0)

          IP = MOD(I0,100)
          I0 = I0/100
          IW = MOD(I0,100)
          IR = I0/100

        END SUBROUTINE
C*MODULE RSM_STOP  *DECK RSM_PUTIO
C>
C>    @brief   save i/o data from dictionary file
C>
C>    @author  Daisuke Yokogawa
C>
        SUBROUTINE RSM_PUTIO

          INTEGER :: IRIWIP
C
C         IR, IW, IP
C
          IRIWIP = 10000*IR + 100*IW  + IP
          CALL RSM_DAWR(IRIWIP,1,501,0)
        END SUBROUTINE

      END MODULE

C*MODULE RISM_MOD
C>
C>    @brief   module for RISM calculation
C>
C>    @author  Daisuke Yokogawa
C>
      MODULE RISM_MOD
      IMPLICIT NONE
C
C     INFORMATION ABOUT SOLUTE and SOLVENT
C
C     NTAB = NGRID (0 to NGRID-1)
C          = NRAD
C
      INTEGER, SAVE        :: NU, NV, NTAB, NSHEL0A       ! G1
C
C     PARAMETER FOR FFT
C
      INTEGER, SAVE        :: IGRD_TYP                    ! G1
C
C     PARAMETERS OF SOLVENT
C
      DOUBLE PRECISION     :: BETA                        ! G2
C
C     PARAMETERS FOR RISM
C
      INTEGER              :: NRISM
      DOUBLE PRECISION     :: ALPA, CLTYP                 ! G1
      INTEGER              :: IRSM_DAMP
C
C     PARAMETER FOR RISM-SCF
C
      INTEGER              :: IREST                       ! G1
C
C     FREE ENERGYIES
C
      DOUBLE PRECISION     :: FHNC, FGF, FKH, FPSE2, FPSE3, EUU, FRISM  !
      DOUBLE PRECISION     :: EUV_E2, EUV_E8, E0, ESTATE(10)
      DOUBLE PRECISION     :: EBOX(22)
C
C     PARAMETER FOR PARALLEL CALCULATION (GLOBAL)
C
      LOGICAL, SAVE        :: GOPARR, MASWRK, DSKWRK
      INTEGER              :: NPROC, ME, MASTER
C
C     RISM PART
C
      INTEGER              :: JDIIS, NEXP, MXDIIS, LJRULE ! G1
      INTEGER, ALLOCATABLE :: INDU(:), INDV(:), IJUV(:,:)
      DOUBLE PRECISION     :: DELEX, XMN, THR, RMAX       ! G1

      DOUBLE PRECISION, ALLOCATABLE :: QU(:),QV(:)
      DOUBLE PRECISION, ALLOCATABLE :: ET0(:,:),FQS(:,:),FQLK(:,:),
     *                                 CSB(:,:),FNEP(:,:),WU(:,:,:),
     *                                 HWV(:,:,:),DENV(:),DHK(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: RTAB(:), RKTAB(:)
C
C     LOOKUP TABLE
C
      INTEGER              :: IDISP,  IREP
      INTEGER              :: NP_DIS, NP_REP
      DOUBLE PRECISION              :: EU_DIS, EU_REP
      DOUBLE PRECISION, ALLOCATABLE :: TAB_DIS(:,:,:), TAB_REP(:,:,:),
     *                                  DELSI_DIS(:,:), DELSI_REP(:,:)
C
C     BRIDGE
C
      DOUBLE PRECISION, SAVE        :: FREF, FPT, FCH0             ! G2
      DOUBLE PRECISION, ALLOCATABLE :: HREF(:,:),BRG(:,:),FRLJ(:),
     *                                  HKREF(:,:)
C
C     INTEGRATION
C
      DOUBLE PRECISION              :: RMX_SYMP, RKMX_SYMP         ! G2
C
C     MCQDPT
C
      INTEGER              :: IROOT
C
C     TDDFT
C
      DOUBLE PRECISION              :: EDIFMX
C
      CONTAINS
C*MODULE RISM_MOD  *DECK RSM_ISTIEND
C>
C>    @brief   give starting and ending points for parallel calculation
C>
C>    @author  Daisuke Yokogawa
C>
C>    @param   NUM: the number of loop to be parallelized
C>             ISTRT: starting point
C>             IEND: ending point
C>             NDIM: the size of loop
C>
        SUBROUTINE RSM_ISTIEND(NUM,ISTRT,IEND,NDIM)

        INTEGER, INTENT(IN   ) :: NUM
        INTEGER, INTENT(  OUT) :: ISTRT, IEND, NDIM

        INTEGER NMOD, NDIV

        NMOD=MOD(NUM,NPROC)
        NDIV=FLOOR(REAL(NUM/NPROC,8))

        IF(ME.LT.NMOD) THEN
           ISTRT=ME*(NDIV+1)+1
           IEND =ISTRT+NDIV
           NDIM =NDIV+1
        ELSEIF(ME.EQ.NMOD) THEN
           ISTRT=ME*(NDIV+1)+1
           IEND =ISTRT+NDIV-1
           NDIM =NDIV
        ELSE
           ISTRT=ME*NDIV+NMOD+1
           IEND =ISTRT+NDIV-1
           NDIM =NDIV
        ENDIF

        RETURN
        END SUBROUTINE
C*MODULE RISM_MOD  *DECK RSM_PUTPAR
C>
C>    @brief   save parameters for RISM calculation
C>
C>    @author  Daisuke Yokogawa
C>
        SUBROUTINE RSM_PUTPAR

        INTEGER  :: I0, I1

        I0 = 100*NU    + NV
        I1 = 10000*LJRULE + 1000*IRSM_DAMP + 100*IREST + JDIIS

        CALL RSM_DAWR(I0       ,    1,551,0)
        CALL RSM_DAWR(NTAB     ,    1,552,0)
        CALL RSM_DAWR(NSHEL0A  ,    1,553,0)
        CALL RSM_DAWR(IGRD_TYP ,    1,554,0)
        CALL RSM_DAWR(ALPA     ,    1,555,0)
        CALL RSM_DAWR(CLTYP    ,    1,556,0)
        CALL RSM_DAWR(I1       ,    1,557,0)
        CALL RSM_DAWR(NEXP     ,    1,558,0)
        CALL RSM_DAWR(MXDIIS   ,    1,559,0)
        CALL RSM_DAWR(DELEX    ,    1,560,0)
        CALL RSM_DAWR(XMN      ,    1,561,0)
        CALL RSM_DAWR(THR      ,    1,562,0)
        CALL RSM_DAWR(BETA     ,    1,563,0)
        CALL RSM_DAWR(FREF     ,    1,564,0)
        CALL RSM_DAWR(FPT      ,    1,565,0)
        CALL RSM_DAWR(FCH0     ,    1,566,0)
        CALL RSM_DAWR(RMX_SYMP ,    1,567,0)
        CALL RSM_DAWR(RKMX_SYMP,    1,568,0)
        CALL RSM_DAWR(IROOT    ,    1,570,0)
        CALL RSM_DAWR(EDIFMX   ,    1,571,0)
        CALL RSM_DAWR(RMAX     ,    1,573,0)

        END SUBROUTINE
C*MODULE RISM_MOD  *DECK RSM_PUTPAR
C>
C>    @brief   read parameters for RISM calculation
C>
C>    @author  Daisuke Yokogawa
C>
        SUBROUTINE RSM_SETPAR

        INTEGER  :: I0, I1

        CALL RSM_DARE(I0       ,    1,551,0)
        CALL RSM_DARE(NTAB     ,    1,552,0)
        CALL RSM_DARE(NSHEL0A  ,    1,553,0)
        CALL RSM_DARE(IGRD_TYP ,    1,554,0)
        CALL RSM_DARE(ALPA     ,    1,555,0)
        CALL RSM_DARE(CLTYP    ,    1,556,0)
        CALL RSM_DARE(I1       ,    1,557,0)
        CALL RSM_DARE(NEXP     ,    1,558,0)
        CALL RSM_DARE(MXDIIS   ,    1,559,0)
        CALL RSM_DARE(DELEX    ,    1,560,0)
        CALL RSM_DARE(XMN      ,    1,561,0)
        CALL RSM_DARE(THR      ,    1,562,0)
        CALL RSM_DARE(BETA     ,    1,563,0)
        CALL RSM_DARE(FREF     ,    1,564,0)
        CALL RSM_DARE(FPT      ,    1,565,0)
        CALL RSM_DARE(FCH0     ,    1,566,0)
        CALL RSM_DARE(RMX_SYMP ,    1,567,0)
        CALL RSM_DARE(RKMX_SYMP,    1,568,0)
        CALL RSM_DARE(IROOT    ,    1,570,0)
        CALL RSM_DARE(EDIFMX   ,    1,571,0)
        CALL RSM_DARE(RMAX     ,    1,573,0)

        NV = MOD(I0,100)
        NU = I0/100

        LJRULE    = I1/10000
        I1        = I1 - 10000*LJRULE
        IRSM_DAMP = I1/1000
        I1        = I1 - 1000*IRSM_DAMP
        IREST     = I1/100
        JDIIS     = I1 - 100 *IREST

        END SUBROUTINE
C*MODULE RISM_MOD  *DECK RSM_PUTPAR
C>
C>    @brief   read matrix data for RISM calculation
C>
C>    @author  Daisuke Yokogawa
C>
        SUBROUTINE RSM_SETVAR

        INTEGER :: L01, L02, L03, L05, L11, L12, L13, L14, L21, L22, L23

        L01 = NU
        L02 = NV
        L03 = NTAB
        L05 = NTAB
        L11 = NU*NU
        L12 = NV*NV
        L13 = NU*NV
        L14 = (NSHEL0A+1)*(NSHEL0A+1)
        L21 = NTAB*NV*NV
        L22 = NTAB*NU*NV
        L23 = NTAB*NU*NU

        ALLOCATE(ET0(NTAB,L13),
     *           FQS(NTAB,L13),FQLK(NTAB,L13),CSB(NTAB,L13),
     *           FNEP(NTAB,L13),WU(NTAB,NU,NU),QU(L01),
     *           QV(L02),HWV(NTAB,NV,NV),DENV(L02),DHK(NTAB,L13),
     *           RTAB(L03),RKTAB(L03),
     *           HREF(NTAB,L13),BRG(NTAB,L13),FRLJ(NU),HKREF(NTAB,L13))
        ALLOCATE(INDV(L13),IJUV(NU,NV),INDU(L13))
C
C       PEAD STORED DATA
C
        CALL RSM_DARE(DENV,L02, 4,0)
        CALL RSM_DARE(HWV ,L21,11,0)
        CALL RSM_DARE(ET0 ,L22,16,0)
        CALL RSM_DARE(CSB ,L22,17,0)
        CALL RSM_DARE(FQLK,L22,18,0)
        CALL RSM_DARE(INDU,L13,26,0)
        CALL RSM_DARE(INDV,L13,27,0)
        CALL RSM_DARE(IJUV,L13,28,0)

        CALL RSM_DARE(QV   ,L02,  3,0)
        CALL RSM_DARE(RTAB ,L05,  8,0)
        CALL RSM_DARE(RKTAB,L05,  9,0)
        CALL RSM_DARE(FNEP ,L22, 12,0)
        CALL RSM_DARE(BRG  ,L22, 20,0)
        CALL RSM_DARE(HREF ,L22, 21,0)
        CALL RSM_DARE(WU   ,L23, 47,0)
!       CALL RSM_DARE(FQS  ,L22, 97,0)

        END SUBROUTINE
C*MODULE RISM_MOD  *DECK RSM_PUTPAR
C>
C>    @brief   save matrix data for RISM calculation
C>
C>    @author  Daisuke Yokogawa
C>
        SUBROUTINE RSM_PUTVAR

        INTEGER :: L22

        L22 = NTAB*NU*NV
C
C       STORE DATA
C
        CALL RSM_DAWR(ET0 ,L22,16,0)
        CALL RSM_DAWR(CSB ,L22,17,0)
        CALL RSM_DAWR(FQLK,L22,18,0)
C       CALL RSM_DAWR(FQS ,L22,19,0)
C
        DEALLOCATE(ET0,FQS,FQLK,CSB,FNEP,WU,QU,QV,HWV,DENV,DHK,RTAB,
     *             RKTAB,HREF,BRG,FRLJ,HKREF,INDV,IJUV,INDU)

        END SUBROUTINE
C*MODULE RISM_MOD  *DECK RSM_PUTENE
C>
C>    @brief   save solvation free energies
C>
C>    @author  Daisuke Yokogawa
C>
        SUBROUTINE RSM_PUTENE

        INTEGER :: I

        EBOX( 1) = FHNC
        EBOX( 2) = FGF
        EBOX( 3) = FKH
        EBOX( 4) = EUU
        EBOX( 5) = FRISM
        EBOX( 6) = EUV_E2
        EBOX( 7) = EUV_E8
        EBOX( 8) = E0
        EBOX( 9) = FPSE2
        EBOX(10) = FPSE3

        DO I=1, 10
          EBOX(I+10) = ESTATE(I)
        ENDDO

        CALL RSM_DAWR(EBOX ,20,569,0)

        END SUBROUTINE
C*MODULE RISM_MOD  *DECK RSM_SETENE
C>
C>    @brief   read solvation free energies
C>
C>    @author  Daisuke Yokogawa
C>
        SUBROUTINE RSM_SETENE

        INTEGER :: I

        CALL RSM_DARE(EBOX ,20,569,0)

        FHNC   = EBOX( 1)
        FGF    = EBOX( 2)
        FKH    = EBOX( 3)
        EUU    = EBOX( 4)
        FRISM  = EBOX( 5)
        EUV_E2 = EBOX( 6)
        EUV_E8 = EBOX( 7)
        E0     = EBOX( 8)
        FPSE2  = EBOX( 9)
        FPSE3  = EBOX(10)

        DO I=1, 10
          ESTATE(I) = EBOX(I+10)
        ENDDO

        END SUBROUTINE
C
      END MODULE
C*MODULE RISMUN_MOD
C>
C>    @brief   module for grid information of uniform-grid RISM
C>
C>    @author  Daisuke Yokogawa
C>
      MODULE RISMUN_MOD
      IMPLICIT NONE
C
C     GRID DATA
C
      DOUBLE PRECISION, SAVE        :: DELTR, DELTK

      END MODULE
C*MODULE RISMLN_MOD
C>
C>    @brief   module for grid information of log-grid RISM
C>
C>    @author  Daisuke Yokogawa
C>
      MODULE RISMLN_MOD
      IMPLICIT NONE
C
C
C
      INTEGER             :: NRAD,MRAD,MXPNT
      DOUBLE PRECISION    :: DRHO,RHOM
C
      DOUBLE PRECISION, ALLOCATABLE :: R32(:)

      END MODULE
C*MODULE ABSMOD
C>
C>    @brief   module for auxiliary basis set
C>
C>    @author  Daisuke Yokogawa
C>
      MODULE ABSMOD
      IMPLICIT NONE
C
      INTEGER              :: NRAD, NMUESP, NMUREP, IFITU              ! G1
      DOUBLE PRECISION     :: XMNFIT, XMUMX                            ! G1
C
      INTEGER, SAVE        :: MXG, MXS, NSHEL0A, NQMT0A, JABS,         ! G1
     *                        MXANG_ABS

      INTEGER, ALLOCATABLE :: KSTRT0A(:),KATM0A(:),KTYP0A(:),KNG0A(:), ! G1
     *                        KLOC0A(:),KMIN0A(:),KMAX0A(:)
      DOUBLE PRECISION, ALLOCATABLE :: EX0A(:),COF0A(:)                         ! G1

      INTEGER, ALLOCATABLE :: KATM_QMT0A(:)                            ! G1
      DOUBLE PRECISION, SAVE        :: QMTTOL                                   ! G1
      DOUBLE PRECISION, ALLOCATABLE :: TQ(:,:),FMAT(:,:)                        ! G1
C
      DOUBLE PRECISION              :: TRN1(9,10),TRN2(10,9),TRN3(10,10)        ! G1
C
      CONTAINS
C*MODULE ABSMOD  *DECK RSM_SETABS
C>
C>    @brief   read information about ABS
C>
C>    @author  Daisuke Yokogawa
C>
        SUBROUTINE RSM_SETABS

        INTEGER :: I0

        CALL RSM_DARE(NRAD   ,     1         ,596,0)
        CALL RSM_DARE(I0     ,     1         ,597,0)
        CALL RSM_DARE(IFITU  ,     1         ,598,0)
        CALL RSM_DARE(XMNFIT ,     1         ,599,0)
        CALL RSM_DARE(XMUMX  ,     1         ,600,0)
        CALL RSM_DARE(MXG    ,     1         ,576,0)
        CALL RSM_DARE(MXS    ,     1         ,577,0)
        CALL RSM_DARE(NSHEL0A,     1         ,578,0)
        CALL RSM_DARE(NQMT0A ,     1         ,579,0)
        CALL RSM_DARE(JABS   ,     1         ,580,0)
        CALL RSM_DARE(QMTTOL ,     1         ,581,0)

        NMUREP = MOD(I0,100)
        NMUESP = I0/100

        ALLOCATE(KSTRT0A(MXS),KATM0A(MXS),KTYP0A(MXS),KNG0A(MXS),
     *           KLOC0A(MXS),KMIN0A(MXS),KMAX0A(MXS))
        ALLOCATE(EX0A(MXG),COF0A(MXG))
        ALLOCATE(KATM_QMT0A(MXS))
        ALLOCATE(TQ(NSHEL0A,NSHEL0A),FMAT(NQMT0A,NRAD))

        CALL RSM_DARE(KSTRT0A   ,MXS            ,582,0)
        CALL RSM_DARE(KATM0A    ,MXS            ,583,0)
        CALL RSM_DARE(KTYP0A    ,MXS            ,584,0)
        CALL RSM_DARE(KNG0A     ,MXS            ,585,0)
        CALL RSM_DARE(KLOC0A    ,MXS            ,586,0)
        CALL RSM_DARE(KMIN0A    ,MXS            ,587,0)
        CALL RSM_DARE(KMAX0A    ,MXS            ,588,0)
        CALL RSM_DARE(EX0A      ,MXG            ,589,0)
        CALL RSM_DARE(COF0A     ,MXG            ,590,0)
        CALL RSM_DARE(KATM_QMT0A,MXS            ,591,0)
        CALL RSM_DARE(TQ        ,NSHEL0A*NSHEL0A,592,0)
        CALL RSM_DARE(TRN1      ,90             ,593,0)
        CALL RSM_DARE(TRN2      ,90             ,594,0)
        CALL RSM_DARE(TRN3      ,100            ,595,0)
        CALL RSM_DARE(FMAT      ,NQMT0A*NRAD    ,601,0)
        CALL RSM_DARE(MXANG_ABS ,1              ,602,0)

        END SUBROUTINE
C*MODULE ABSMOD  *DECK RSM_PUTABS
C>
C>    @brief   save information about ABS
C>
C>    @author  Daisuke Yokogawa
C>
        SUBROUTINE RSM_PUTABS(I1,I2,I3,I4,I5,I6,I7,I8,V1,V2,V3,V4,V5,V6,
     *                        V7)

        INTEGER,          INTENT(IN   ) :: I1(*),I2(*),I3(*),I4(*),
     *                                     I5(*),I6(*),I7(*),I8(*)
        DOUBLE PRECISION, INTENT(IN   ) :: V1(*),V2(*),V3(*),V4(*),
     *                                     V5(*),V6(*),V7(*)

        INTEGER :: I0

        I0 = 100*NMUESP + NMUREP

        CALL RSM_DAWR(NRAD   ,     1         ,596,0)
        CALL RSM_DAWR(I0     ,     1         ,597,0)
        CALL RSM_DAWR(IFITU  ,     1         ,598,0)
        CALL RSM_DAWR(XMNFIT ,     1         ,599,0)
        CALL RSM_DAWR(XMUMX  ,     1         ,600,0)
        CALL RSM_DAWR(MXG    ,     1         ,576,0)
        CALL RSM_DAWR(MXS    ,     1         ,577,0)
        CALL RSM_DAWR(NSHEL0A,     1         ,578,0)
        CALL RSM_DAWR(NQMT0A ,     1         ,579,0)
        CALL RSM_DAWR(JABS   ,     1         ,580,0)
        CALL RSM_DAWR(QMTTOL ,     1         ,581,0)
        CALL RSM_DAWR(I1     ,MXS            ,582,0)
        CALL RSM_DAWR(I2     ,MXS            ,583,0)
        CALL RSM_DAWR(I3     ,MXS            ,584,0)
        CALL RSM_DAWR(I4     ,MXS            ,585,0)
        CALL RSM_DAWR(I5     ,MXS            ,586,0)
        CALL RSM_DAWR(I6     ,MXS            ,587,0)
        CALL RSM_DAWR(I7     ,MXS            ,588,0)
        CALL RSM_DAWR(V1     ,MXG            ,589,0)
        CALL RSM_DAWR(V2     ,MXG            ,590,0)
        CALL RSM_DAWR(I8     ,MXS            ,591,0)
        CALL RSM_DAWR(V3     ,NSHEL0A*NSHEL0A,592,0)
        CALL RSM_DAWR(V4     ,90             ,593,0)
        CALL RSM_DAWR(V5     ,90             ,594,0)
        CALL RSM_DAWR(V6     ,100            ,595,0)
        CALL RSM_DAWR(V7     ,NQMT0A*NRAD    ,601,0)
        CALL RSM_DAWR(MXANG_ABS,    1        ,602,0)

        END SUBROUTINE
C*MODULE ABSMOD  *DECK RSM_PUTABS
C>
C>    @brief   deallocate memories of ABS
C>
C>    @author  Daisuke Yokogawa
C>
        SUBROUTINE RSM_CLOSEABS

        DEALLOCATE(KSTRT0A,KATM0A,KTYP0A,KNG0A,KLOC0A,KMIN0A,KMAX0A,
     *             EX0A,COF0A,KATM_QMT0A,TQ,FMAT)
        END SUBROUTINE
C
      END MODULE
C*MODULE RSMSED
C>
C>    @brief   module for RISM-SCF calculation
C>
C>    @author  Daisuke Yokogawa
C>
      MODULE RSMSED
      IMPLICIT NONE
C
C     RISM-SCF-CALCULATION (GLOBAL)
C
      INTEGER, SAVE          :: ICNTRL(18)
      INTEGER, SAVE          :: IRISM        ! G2
C
C     PRINT OPTION
C
      INTEGER, SAVE          :: NPRINT       ! G2
C
C     QM CALCULATION TYPE  (GLOBAL)
C
C     ----- CC calculation
C
      INTEGER, SAVE          :: ICC          ! G2
C
C     ----- MCSCF calculation
C
      INTEGER, SAVE          :: IMCSCF              ! G2
C
C     ----- DFT calculation
C
      INTEGER, SAVE          :: IDFT, ITDDFT, IZINC ! G2
      INTEGER, SAVE          :: ITER_RSMTD          ! G2
      DOUBLE PRECISION, SAVE :: TDTYP               ! G2
C
C     ----- CONICAL INTERSECTIOn SERCH
C
      INTEGER, SAVE          :: ICONCL              ! G2
C
C     RISM-SCF CYCLE
C
      DOUBLE PRECISION, SAVE :: RCNTRL(8)
      INTEGER, SAVE          :: NSERCH            ! G2
      INTEGER, SAVE          :: ITRRSM, IEXP      ! G2
      DOUBLE PRECISION, SAVE :: FEXP, ETOT, EINT1, EINT2, EINT3, G1, G3 ! G2
      INTEGER, SAVE          :: ITDCONV
      INTEGER                :: ISKIP_RSM
      INTEGER                :: IRSM_COUNT
C
C     SPECTROSCOPIC CALCULATION
C
      INTEGER, SAVE        :: ISPEC        ! G2
C
C     ----- EQUILIBRIUM/NONEQUILIBRIUM CALCULATION
C
      INTEGER, SAVE        :: IEQ          ! G2
C
C     ----- INITIAL AND FINAL STATES
C
      INTEGER, SAVE        :: NROOT_I, MUL_I
C
C     D correction
C
      INTEGER, SAVE        :: IDCORR       ! G2
C
C     /RSINFA/
C
      INTEGER              :: INFOMAT(9)
      INTEGER              :: NAT, ICH, MUL, NUM, NQMT, NE, NA, NB,
     *                        NROOT  ! G1
      INTEGER, ALLOCATABLE :: IAN(:),IZCORE(:)                      ! G1
      DOUBLE PRECISION, ALLOCATABLE :: ZAN(:),C(:,:),ANAM(:)        ! G1
C
C     RISM-SCF PART
C
      DOUBLE PRECISION, ALLOCATABLE :: QMT(:),VPR(:),DMT(:),VCR(:)

      CONTAINS
C*MODULE RSMSED  *DECK RSM_SCF_SETVAR
C>
C>    @brief   read matrix data for RISM-SCF calculation
C>
C>    @author  Daisuke Yokogawa
C>
        SUBROUTINE RSM_SCF_SETVAR

        USE ABSMOD, ONLY : NQMT0A

        INTEGER :: L03, L14

        L03 = NQMT0A+1
        L14 = (NQMT0A+1)*(NQMT0A+1)

        ALLOCATE(QMT(L14),VPR(L03),DMT(L03),VCR(L03))

        CALL RSM_DARE(VPR ,L03, 6,0)
        CALL RSM_DARE(DMT ,L03,24,0)
        CALL RSM_DARE(QMT ,L14,41,0)

        END SUBROUTINE
C*MODULE RSMSED  *DECK RSM_SCF_PUTVAR
C>
C>    @brief   save matrix data for RISM-SCF calculation
C>
C>    @author  Daisuke Yokogawa
C>
        SUBROUTINE RSM_SCF_PUTVAR

        USE ABSMOD, ONLY : NQMT0A

        INTEGER :: L03, L14

        L03 = NQMT0A+1
        L14 = (NQMT0A+1)*(NQMT0A+1)
        CALL RSM_DAWR(QMT   ,  L14, 41,0)
        CALL RSM_DAWR(VPR   ,  L03,  6,0)
        CALL RSM_DAWR(DMT   ,  L03, 24,0)

        DEALLOCATE(QMT,VPR,DMT,VCR)

        END SUBROUTINE
C*MODULE RSMSED  *DECK RSM_SCF_SETPAR
C>
C>    @brief   read parameters for RISM-SCF calculation
C>
C>    @author  Daisuke Yokogawa
C>
        SUBROUTINE RSM_SCF_SETPAR

        CALL RSM_DARE(ICNTRL,   18,621,0)
        CALL RSM_DARE(RCNTRL,    8,624,0)

        IRISM     = ICNTRL(1)
        ICC       = ICNTRL(2)
        IDFT      = ICNTRL(3)
        ITDDFT    = ICNTRL(4)
        IDCORR    = ICNTRL(5)
        IEQ       = ICNTRL(6)
        ITRRSM    = ICNTRL(7)
        IEXP      = ICNTRL(8)
        NSERCH    = ICNTRL(9)
        ITER_RSMTD= ICNTRL(10)
        NPRINT    = ICNTRL(11)
        ISPEC     = ICNTRL(12)
        IZINC     = ICNTRL(13)
        ISKIP_RSM = ICNTRL(14)
        IRSM_COUNT= ICNTRL(15)
        ITDCONV   = ICNTRL(16)
        ICONCL    = ICNTRL(17)
        IMCSCF    = ICNTRL(18)

        FEXP      = RCNTRL(1)
        ETOT      = RCNTRL(2)
        EINT1     = RCNTRL(3)
        EINT2     = RCNTRL(4)
        EINT3     = RCNTRL(5)
        TDTYP     = RCNTRL(6)
        G1        = RCNTRL(7)
        G3        = RCNTRL(8)

        END SUBROUTINE
C*MODULE RSMSED  *DECK RSM_SCF_PUTPAR
C>
C>    @brief   save parameters for RISM-SCF calculation
C>
C>    @author  Daisuke Yokogawa
C>
        SUBROUTINE RSM_SCF_PUTPAR

        ICNTRL( 1) = IRISM
        ICNTRL( 2) = ICC
        ICNTRL( 3) = IDFT
        ICNTRL( 4) = ITDDFT
        ICNTRL( 5) = IDCORR
        ICNTRL( 6) = IEQ
        ICNTRL( 7) = ITRRSM
        ICNTRL( 8) = IEXP
        ICNTRL( 9) = NSERCH
        ICNTRL(10) = ITER_RSMTD
        ICNTRL(11) = NPRINT
        ICNTRL(12) = ISPEC
        ICNTRL(13) = IZINC
        ICNTRL(14) = ISKIP_RSM
        ICNTRL(15) = IRSM_COUNT
        ICNTRL(16) = ITDCONV
        ICNTRL(17) = ICONCL
        ICNTRL(18) = IMCSCF

        RCNTRL( 1) = FEXP
        RCNTRL( 2) = ETOT
        RCNTRL( 3) = EINT1
        RCNTRL( 4) = EINT2
        RCNTRL( 5) = EINT3
        RCNTRL( 6) = TDTYP
        RCNTRL( 7) = G1
        RCNTRL( 8) = G3

        CALL RSM_DAWR(ICNTRL,   18,621,0)
        CALL RSM_DAWR(RCNTRL,    8,624,0)

        END SUBROUTINE
C*MODULE RSMSED  *DECK RSM_SCF_SETINF1
C>
C>    @brief   read information of molecule for RISM-SCF calculation
C>
C>    @author  Daisuke Yokogawa
C>
        SUBROUTINE RSM_SCF_SETINF1
C
        CALL RSM_DARE(INFOMAT,9,511,0)
        NAT  = INFOMAT(1)
        ICH  = INFOMAT(2)
        MUL  = INFOMAT(3)
        NUM  = INFOMAT(4)
        NQMT = INFOMAT(5)
        NE   = INFOMAT(6)
        NA   = INFOMAT(7)
        NB   = INFOMAT(8)
        NROOT= INFOMAT(9)

        ALLOCATE(IAN(NAT),IZCORE(NAT),ZAN(NAT),ANAM(NAT))

        CALL RSM_DARE(IAN   ,  NAT,518,0)  ! IAN
        CALL RSM_DARE(IZCORE,  NAT,519,0)  ! IZCORE
        CALL RSM_DARE(ZAN   ,  NAT,520,0)  ! ZAN
        CALL RSM_DARE(ANAM  ,  NAT,522,0)
C
        END SUBROUTINE
C*MODULE RSMSED  *DECK RSM_SCF_PUTINF1
C>
C>    @brief   save information of molecule for RISM-SCF calculation
C>
C>    @author  Daisuke Yokogawa
C>
C>    @param   I1: IAN
C>             I2: IZCORE
C>             R1: ZAN
C>             R2: ANAM
C>
        SUBROUTINE RSM_SCF_PUTINF1(I1,I2,R1,R2)
C
        INTEGER,          INTENT(IN   ) :: I1(*), I2(*)
        DOUBLE PRECISION, INTENT(IN   ) :: R1(*), R2(*)
C
        INTEGER :: I0

        INFOMAT(1) = NAT
        INFOMAT(2) = ICH
        INFOMAT(3) = MUL
        INFOMAT(4) = NUM
        INFOMAT(5) = NQMT
        INFOMAT(6) = NE
        INFOMAT(7) = NA
        INFOMAT(8) = NB
        INFOMAT(9) = NROOT

        CALL RSM_DAWR(INFOMAT, 9,511,0)
        CALL RSM_DAWR(I1  ,  NAT,518,0)  ! IAN
        CALL RSM_DAWR(I2  ,  NAT,519,0)  ! IZCORE
        CALL RSM_DAWR(R1  ,  NAT,520,0)  ! ZAN
        CALL RSM_DAWR(R2  ,  NAT,522,0)  ! ANAM

        END SUBROUTINE
C*MODULE RSMSED  *DECK RSM_SCF_CLOSEINF1
C>
C>    @brief   deallocate memolies for INF1
C>
C>    @author  Daisuke Yokogawa
C>
        SUBROUTINE RSM_SCF_CLOSEINF1
C
        DEALLOCATE(IAN,IZCORE,ZAN,ANAM)

        END SUBROUTINE
C*MODULE RSMSED  *DECK RSM_SCF_SETINF2
C>
C>    @brief   read information of molecule for RISM-SCF calculation
C>
C>    @author  Daisuke Yokogawa
C>
        SUBROUTINE RSM_SCF_SETINF2
C
        ALLOCATE(C(3,NAT))
C
        CALL RSM_DARE(C,3*NAT,521,0)  ! C

        END SUBROUTINE
C*MODULE RSMSED  *DECK RSM_SCF_PUTINF2
C>
C>    @brief   save information of molecule for RISM-SCF calculation
C>
C>    @author  Daisuke Yokogawa
C>
        SUBROUTINE RSM_SCF_PUTINF2(R2)
C
        DOUBLE PRECISION, INTENT(IN   ) :: R2(*)
C
        CALL RSM_DAWR(R2  ,3*NAT,521,0)  ! C

        END SUBROUTINE
C*MODULE RSMSED  *DECK RSM_SCF_CLOSEINF2
C>
C>    @brief   deallocate memories for INF2
C>
C>    @author  Daisuke Yokogawa
C>
        SUBROUTINE RSM_SCF_CLOSEINF2
C
        DEALLOCATE(C)

        END SUBROUTINE
C
      END MODULE
C*MODULE RSM_BSMOD
C>
C>    @brief   module for basis set of QM calculations
C>
C>    @author  Daisuke Yokogawa
C>
      MODULE RSM_BSMOD
      IMPLICIT NONE
C
      INTEGER              :: MXANG, IPURED, IPUREF, NSHELL, NPRIM
      LOGICAL              :: NORM

      INTEGER, ALLOCATABLE :: KSTART(:),KATOM(:),KTYPE(:),KNG(:),
     *                        KLOC(:),KMIN(:),KMAX(:)

      DOUBLE PRECISION, ALLOCATABLE :: EX(:),CS(:),CP(:),CD(:),CF(:),
     *                                 CG(:),CH(:),CI(:)

      CONTAINS
C*MODULE RSM_BSMOD  *DECK RSM_SCF_SETBS
C>
C>    @brief   read information about basis set
C>
C>    @author  Daisuke Yokogawa
C>
        SUBROUTINE RSM_SCF_SETBS

        INTEGER :: I0

        CALL RSM_DARE(I0    ,     1,531,0)
        CALL RSM_DARE(NSHELL,     1,532,0)
        CALL RSM_DARE(NPRIM ,     1,533,0)
        IPUREF = MOD(I0,100)
        I0     = I0/100
        IPURED = MOD(I0,100)
        I0     = I0/100
        MXANG  = MOD(I0,100)
        I0     = I0/100
        IF(I0 == 1) THEN
          NORM = .TRUE.
        ELSE
          NORM = .FALSE.
        ENDIF

        ALLOCATE(KSTART(NSHELL),KATOM(NSHELL),KTYPE(NSHELL),KNG(NSHELL),
     *           KLOC(NSHELL),KMIN(NSHELL),KMAX(NSHELL))
        ALLOCATE(EX(NPRIM),CS(NPRIM),CP(NPRIM),CD(NPRIM),CF(NPRIM),
     *           CG(NPRIM),CH(NPRIM),CI(NPRIM))

        CALL RSM_DARE(KSTART,NSHELL,534,0)
        CALL RSM_DARE(KATOM ,NSHELL,535,0)
        CALL RSM_DARE(KTYPE ,NSHELL,536,0)
        CALL RSM_DARE(KNG   ,NSHELL,537,0)
        CALL RSM_DARE(KLOC  ,NSHELL,538,0)
        CALL RSM_DARE(KMIN  ,NSHELL,539,0)
        CALL RSM_DARE(KMAX  ,NSHELL,540,0)
        CALL RSM_DARE(EX    ,NPRIM ,541,0)
        CALL RSM_DARE(CS    ,NPRIM ,542,0)
        CALL RSM_DARE(CP    ,NPRIM ,543,0)
        CALL RSM_DARE(CD    ,NPRIM ,544,0)
        CALL RSM_DARE(CF    ,NPRIM ,545,0)
        CALL RSM_DARE(CG    ,NPRIM ,546,0)
        CALL RSM_DARE(CH    ,NPRIM ,547,0)
        CALL RSM_DARE(CI    ,NPRIM ,548,0)

        END SUBROUTINE
C*MODULE RSM_BSMOD  *DECK RSM_SCF_PUTBS
C>
C>    @brief   save information about basis set
C>
C>    @author  Daisuke Yokogawa
C>
C>    @param   I1: KSTART
C>             I2: KATOM
C>             I3: KTYPE
C>             I4: KNG
C>             I5: KLOC
C>             I6: KMIN
C>             I7: KMAX
C>             V1: EX
C>             V2: CS
C>             V3: CP
C>             V4: CD
C>             V5: CF
C>             V6: CG
C>             V7: CH
C>             V8: CI
C>
        SUBROUTINE RSM_SCF_PUTBS(I1,I2,I3,I4,I5,I6,I7,V1,V2,V3,V4,V5,V6,
     *                           V7,V8)

        INTEGER,          INTENT(IN   ) :: I1(*),I2(*),I3(*),I4(*),
     *                                     I5(*),I6(*),I7(*)
        DOUBLE PRECISION, INTENT(IN   ) :: V1(*),V2(*),V3(*),V4(*),
     *                                     V5(*),V6(*),V7(*),V8(*)

        INTEGER :: I0, INORM

        IF(NORM) THEN
          INORM = 1
        ELSE
          INORM = 0
        ENDIF

        I0 = 1000000*INORM + 10000*MXANG + 100*IPURED + IPUREF

        CALL RSM_DAWR(I0    ,     1,531,0)
        CALL RSM_DAWR(NSHELL,     1,532,0)
        CALL RSM_DAWR(NPRIM ,     1,533,0)
        CALL RSM_DAWR(I1    ,NSHELL,534,0)
        CALL RSM_DAWR(I2    ,NSHELL,535,0)
        CALL RSM_DAWR(I3    ,NSHELL,536,0)
        CALL RSM_DAWR(I4    ,NSHELL,537,0)
        CALL RSM_DAWR(I5    ,NSHELL,538,0)
        CALL RSM_DAWR(I6    ,NSHELL,539,0)
        CALL RSM_DAWR(I7    ,NSHELL,540,0)
        CALL RSM_DAWR(V1    ,NPRIM ,541,0)
        CALL RSM_DAWR(V2    ,NPRIM ,542,0)
        CALL RSM_DAWR(V3    ,NPRIM ,543,0)
        CALL RSM_DAWR(V4    ,NPRIM ,544,0)
        CALL RSM_DAWR(V5    ,NPRIM ,545,0)
        CALL RSM_DAWR(V6    ,NPRIM ,546,0)
        CALL RSM_DAWR(V7    ,NPRIM ,547,0)
        CALL RSM_DAWR(V8    ,NPRIM ,548,0)

        END SUBROUTINE
C*MODULE RSM_BSMOD  *DECK RSM_SCF_CLOSEBS
C>
C>    @brief   deallocate memories for basis set
C>
C>    @author  Daisuke Yokogawa
C>
        SUBROUTINE RSM_SCF_CLOSEBS

        DEALLOCATE(KSTART,KATOM,KTYPE,KNG,KLOC,KMIN,KMAX)
        DEALLOCATE(EX,CS,CP,CD,CF,CG,CH,CI)

        END SUBROUTINE

      END MODULE

C*MODULE RSM_STOP  *DECK RSMSTP
      SUBROUTINE RSMSTP
C
      IMPLICIT NONE
C
      WRITE(6,*) 'RISM-SCF-cSED RUNS ARE NOT ENABLED, PLEASE RE-COMPILE'
      CALL ABRT
C
      RETURN
      END
C
C*MODULE RSM_GMS_GMSPAR  *DECK RSM_PARIN
C>
C>    @brief   set parameters for parallel calculation
C>
C>    @author  Daisuke Yokogawa
C>
      SUBROUTINE RSM_PARIN

      USE RISM_MOD, ONLY : GOPARR, MASWRK, NPROC, ME, MASTER, DSKWRK

      IMPLICIT NONE

      LOGICAL GOPARR_GMS,DSKWRK_GMS,MASWRK_GMS
      INTEGER ME_GMS,MASTER_GMS,NPROC_GMS,IBTYP,IPTIM
      COMMON /PAR   / ME_GMS,MASTER_GMS,NPROC_GMS,IBTYP,IPTIM,
     *                GOPARR_GMS,DSKWRK_GMS,MASWRK_GMS

      ME     = ME_GMS
      MASTER = MASTER_GMS
      NPROC  = NPROC_GMS

      GOPARR = GOPARR_GMS
      DSKWRK = DSKWRK_GMS
      MASWRK = MASWRK_GMS

      END SUBROUTINE
C*MODULE RSM_GMS_GMSPAR  *DECK RSM_BCAST
C>
C>    @brief   broadcast the data
C>
C>    @author  Daisuke Yokogawa
C>
      SUBROUTINE RSM_BCAST(MSGTAG,TYPE,BUFF,LEN,FROM)
      IMPLICIT NONE

      INTEGER,     INTENT(IN   ) :: MSGTAG, BUFF(*), LEN, FROM
      CHARACTER*1, INTENT(IN   ) :: TYPE

      CALL DDI_BCAST( MSGTAG, TYPE, BUFF, LEN, FROM )

      END SUBROUTINE

C*MODULE RSM_GMS_GETDAT  *DECK RSM_GETDAT0
C>
C>    @brief   get the information about IO
C>
C>    @author  Daisuke Yokogawa
C>
      SUBROUTINE RSM_GETDAT0
C
      USE RISM_IO, ONLY : IR, IW, IP
      USE RSMSED,  ONLY : IRISM

      IMPLICIT NONE

      INTEGER IR_GMS,IW_GMS,IP_GMS,IJK,IJKT,IDAF,NAV,IODA
      COMMON /IOFILE/ IR_GMS,IW_GMS,IP_GMS,IJK,IJKT,IDAF,NAV,IODA(950)

      INTEGER NEVALS,NGLEVL,NHLEVL
      DOUBLE PRECISION RUNTYP,EXETYP,CHECK
      COMMON /RUNOPT/ RUNTYP,EXETYP,NEVALS,NGLEVL,NHLEVL
      DATA CHECK/8HCHECK   /
C
      IR = IR_GMS
      IW = IW_GMS
      IP = IP_GMS
C
      IRISM = 100
      CALL SEQREW(IR)
      CALL FNDGRP(IR,' $RISMUN',IRISM)
      IF(IRISM /= 0) THEN
        CALL SEQREW(IR)
        CALL FNDGRP(IR,' $RISMLN',IRISM)
      ENDIF
C
C     RISM DOES NOT SUPPORT "EXETYP=CHECK"
C
      IF(EXETYP .EQ. CHECK) IRISM = 100

      END SUBROUTINE
C*MODULE RSM_IOS  *DECK RSM_OPDA
C>
C>    @brief   open dictionary file for RISM-SCF calculation
C>
C>    @author  Daisuke Yokogawa
C>
C>    @param   IFIRST: is it required for initialization? (0: no, 1: yes)
C>

      SUBROUTINE RSM_OPDA(IFIRST)
C
      USE RISM_MOD, ONLY : MASWRK
      USE RISM_IO,  ONLY : IRSDAF, IODARS, MAXIO, IRECLN, IRECST, IFILEN
      USE RSMSED,    ONLY : IRISM
C
C     ---- OPEN RANDOM ACCESS FILE ----
C
      IMPLICIT NONE
C
      INTEGER, INTENT(IN   ) :: IFIRST
C
      INTEGER :: I
C
      CHARACTER*256 FILENM
C
      IF(IRISM.EQ.1) RETURN
      IF (MASWRK) THEN
C
        CALL GMS_GETENV('RISM01',FILENM)
        OPEN (UNIT=IRSDAF, FILE=FILENM, STATUS='UNKNOWN',
     *        ACCESS='DIRECT', FORM='UNFORMATTED',
     *        RECL=8*IRECLN)
C
      END IF
C
      IF(IFIRST == 1) THEN
        IRECST = 1
        DO I = 1,MAXIO
          IODARS(I) = -1
        ENDDO
        IRECST = IRECST + 1
        IF(MASWRK) WRITE(UNIT=IRSDAF, REC=1) IRECST,IODARS,IFILEN
      ELSE
        IF(MASWRK) READ(UNIT=IRSDAF, REC=1) IRECST,IODARS,IFILEN
      ENDIF

      RETURN
      END
C*MODULE RSM_IOS  *DECK RSM_CLDA
C>
C>    @brief   close dictionary file
C>
C>    @author  Daisuke Yokogawa
C>
C>
      SUBROUTINE RSM_CLDA

      USE RISM_IO,  ONLY : IRSDAF
      USE RSMSED,    ONLY : IRISM

      IMPLICIT NONE

      IF(IRISM.EQ.1) RETURN
      CLOSE(IRSDAF)

      RETURN
      END
C*MODULE RSM_IOS  *DECK RSM_DAWR
C>
C>    @brief   save data in dictionary file
C>
C>    @author  Daisuke Yokogawa
C>
C>    @param   V: data array
C>             LEN: the length of the data
C>             NREC: the recording position
C>             IDTYP: not used
C>
      SUBROUTINE RSM_DAWR(V,LEN,NREC,IDTYP)
C
      USE RISM_MOD, ONLY : MASWRK
      USE RISM_IO,  ONLY : IRSDAF, IODARS, IFILEN, IRECST, IRECLN
      USE RSMSED,    ONLY : IRISM
C
      IMPLICIT NONE
C
      INTEGER,          INTENT(IN   ) :: LEN, NREC, IDTYP
      DOUBLE PRECISION, INTENT(IN   ) :: V(LEN)

      LOGICAL                :: NEWREC
      INTEGER                :: N, IST, NS, NSP, LENT, LENW, IF
C
C         WRITE A LOGICAL RECORD ON THE DAF DICTIONARY FILE
C         A LOGICAL RECORD MAY SPAN SEVERAL PHYSICAL RECORDS
C
      IF(IRISM.EQ.1) RETURN
      N = IODARS(NREC)
      IF (N .GT. 0 .AND. LEN .NE. IFILEN(NREC)) GO TO 800
      NEWREC = .FALSE.
C
C     NEW DATA
C
      IF (N <= 0) THEN
        IODARS(NREC) = IRECST
        IFILEN(NREC) = LEN
        NEWREC       = .TRUE.
        IRECST       = IRECST + (LEN-1)/IRECLN + 1
        N            = IODARS(NREC)
      ENDIF
      IST = -IRECLN + 1
      NS = N
      LENT = LEN
  120 CONTINUE
         IST = IST + IRECLN
         IF = IST + LENT - 1
         IF ((IF-IST+1) .GT. IRECLN) IF = IST+IRECLN-1
         NSP = NS
         LENW = IF - IST + 1
         CALL RSM_DAWRT(V(IST),LENW,IRSDAF,NSP)
         LENT = LENT - IRECLN
         NS = NS + 1
         N = NS
      IF (LENT .GE. 1) GO TO 120

      IF (NEWREC .AND. MASWRK)
     *      WRITE(UNIT=IRSDAF,REC=1) IRECST,IODARS,IFILEN
      RETURN
C
  800 CONTINUE
      IF (MASWRK) WRITE (6,9008) NREC,LEN,IFILEN(NREC)
      CALL ABRT
      RETURN
C
 9008 FORMAT(1X,'RSM_DAWR HAS REQUESTED A RECORD WITH LENGTH',
     *       1X,'DIFFERENT THAN BEFORE - ABORT FORCED.'/
     *       1X,'RISM01 RECORD ',I5,' NEW LENGTH =',I5,
     *          ' OLD LENGTH =',I5)
      END
C*MODULE RSM_IOS  *DECK RSM_DARE
C>
C>    @brief   read data in dictionary file
C>
C>    @author  Daisuke Yokogawa
C>
C>    @param   V: data array
C>             LEN: the length of the data
C>             NREC: the recording position
C>             IDTYP: not used
C>
      SUBROUTINE RSM_DARE(V,LEN,NREC,IDTYP)
C
      USE RISM_MOD, ONLY : MASWRK
      USE RISM_IO,  ONLY : IRSDAF, IODARS, IFILEN, IRECST, IRECLN
      USE RSMSED,    ONLY : IRISM
C
      IMPLICIT NONE

      INTEGER,          INTENT(IN   ) :: LEN, NREC, IDTYP
      DOUBLE PRECISION, INTENT(  OUT) :: V(LEN)
C
      INTEGER                :: N, IS, NS, NSP, LENT, LENW, IF
C
      IF(IRISM.EQ.1) RETURN
      N = IODARS(NREC)
      IF(N.EQ.-1) GO TO 800
      IS = -IRECLN + 1
      NS = N
      LENT = LEN

  100 CONTINUE
         IS = IS + IRECLN
         IF = IS + LENT - 1
         IF ((IF-IS+1) .GT. IRECLN) IF = IS + IRECLN - 1
         NSP = NS
         LENW = IF - IS + 1
         CALL RSM_DARD(V(IS),LENW,IRSDAF,NSP)
         LENT = LENT - IRECLN
         NS = NS + 1
         N = NS
      IF (LENT .GE. 1) GO TO 100
      RETURN
C
  800 CONTINUE
      IF (MASWRK) WRITE(6,9000) NREC,LEN
      CALL ABRT
      RETURN
C
 9000 FORMAT(1X,'*** ERROR ***, ATTEMPT TO READ A -RISM01- RECORD',
     *         ' THAT WAS NEVER WRITTEN.'/1X,'NREC,LEN=',I5,I10)
      END

C*MODULE RSM_IOS  *DECK RSM_DAWRT
C>
C>    @brief   save data in dictionary file (called by RSM_DAWR)
C>
C>    @author  Daisuke Yokogawa
C>
C>    @param   V: data array
C>             LEN: the length of the data
C>             IDAF: the file number of the dictionary file
C>             NS: record number
C>
      SUBROUTINE RSM_DAWRT(V,LEN,IDAF,NS)

      USE RISM_MOD, ONLY : MASWRK
      USE RSMSED,    ONLY : IRISM

      IMPLICIT NONE

      INTEGER,          INTENT(IN   ) :: LEN, IDAF, NS
      DOUBLE PRECISION, INTENT(IN   ) :: V(LEN)
C
C     ----- WRITE A PHYSICAL RECORD ON THE DAF -----
C
      IF(IRISM.EQ.1) RETURN
      IF (MASWRK) WRITE (UNIT=IDAF, REC=NS) V
      RETURN
      END

C*MODULE RSM_IOS  *DECK RSM_DARD
C>
C>    @brief   read data in dictionary file (called by RSM_DARE)
C>
C>    @author  Daisuke Yokogawa
C>
C>    @param   V: data array
C>             LEN: the length of the data
C>             IDAF: the file number of the dictionary file
C>             NS: record number
C>
      SUBROUTINE RSM_DARD(V,LEN,IDAF,NS)

      USE RISM_MOD, ONLY : MASWRK, GOPARR, MASTER
      USE RSMSED,    ONLY : IRISM

      IMPLICIT NONE

      INTEGER,          INTENT(IN   ) :: LEN, IDAF, NS
      DOUBLE PRECISION, INTENT(  OUT) :: V(LEN)

      INTEGER                :: ITYP
C
C       READ A PHYSICAL RECORD FROM THE DAF
C
      IF(IRISM.EQ.1) RETURN
      IF (MASWRK) READ (UNIT=IDAF, REC=NS) V
      ITYP = 100 + IDAF
      IF (GOPARR) CALL RSM_BCAST(ITYP,'F',V,LEN,MASTER)
      RETURN
      END
C*MODULE RSM_IOS  *DECK RSM_SQOPN
C>
C>    @brief   open a sequential file for RISM-SCF calculation
C>
C>    @author  Daisuke Yokogawa
C>
C>    @param   IUNIT: file number
C>             FNAME: file name
C>             FSTAT: file status
C>             FMT: file format
C>
      SUBROUTINE RSM_SQOPN(IUNIT,FNAME,FSTAT,FMT)
C
      USE RISM_MOD, ONLY : MASWRK
      USE RSMSED,    ONLY : IRISM
C
      IMPLICIT NONE
C
      INTEGER,       INTENT(IN   ) :: IUNIT
      CHARACTER*(*), INTENT(IN   ) :: FNAME, FSTAT, FMT
C     PARAMETER (MXUNIT=299)
C
      CHARACTER*256 :: FILENM
      CHARACTER*1   :: NULL
C
      INTEGER       :: KOL
C
c     DO NOT OPEN/CLOSE RISM FILES IF NOT DOING RISM
      IF(IRISM.EQ.1) RETURN
      FILENM=' '
      IF (MASWRK) CALL GMS_GETENV(FNAME,FILENM)
      IF (.NOT.MASWRK) RETURN

      NULL = CHAR(0)
      DO KOL=1,256
        IF(FILENM(KOL:KOL).EQ.' '  .OR.
     *     FILENM(KOL:KOL).EQ.NULL) GO TO 2
      ENDDO
      KOL=257
    2 CONTINUE
      IF(KOL.EQ.1) THEN
        WRITE(6,3) FNAME
        CALL ABRT
      END IF
      KOL=KOL-1
      OPEN(UNIT=IUNIT, FILE=FILENM(1:KOL), STATUS=FSTAT,
     *        ACCESS='SEQUENTIAL', FORM=FMT, ERR=4)
      RETURN
C
C         ERROR HANDLING (E.G. FOR ERICFMT, ...)
C
    3 FORMAT(1X,'YOU MUST ASSIGN GENERIC NAME ',A,' WITH A SETENV.')
C
    4 CONTINUE
      IF(FSTAT.EQ.'OLD') THEN
         IF(MASWRK) WRITE(6,5) 'PRE-EXISTING',FNAME,FILENM(1:KOL)
      ELSE
         IF(MASWRK) WRITE(6,5) 'NEW',FNAME,FILENM(1:KOL)
      ENDIF
      CALL ABRT
    5 FORMAT(//1X,'ERROR OPENING ',A,' FILE ',A,','/
     *         1X,'ASSIGNED TO EXPLICIT FILE NAME ',A,','/
     *         1X,'PLEASE CHECK THE -SETENV- FILE ASSIGNMENTS',
     *            ' IN YOUR -RUNGMS- SCRIPT.')
C
      RETURN
      END
C*MODULE RSM_IOS  *DECK RSM_SQOPN
C>
C>    @brief   close the sequential file for RISM-SCF calculation
C>
C>    @author  Daisuke Yokogawa
C>
C>    @param   IUNIT: file number
C>
      SUBROUTINE RSM_SQCLS(IUNIT)
      USE RSMSED,    ONLY : IRISM
      IMPLICIT NONE

      INTEGER, INTENT(IN   ) :: IUNIT

c     DO NOT OPEN/CLOSE RISM FILES IF NOT DOING RISM
      IF(IRISM.EQ.1) RETURN
      CLOSE(IUNIT)

      END

C*MODULE RSM_STOP  *DECK RSM_SCF_STRT1
C>
C>    @brief   store the parameters that are independent of geometries
C>             for RISM-SCF calculation
C>
C>    @author  Daisuke Yokogawa
C>
      SUBROUTINE RSM_SCF_STRT1

      USE RSMSED,    ONLY : IRISM, RSM_SCF_PUTPAR
      USE RISM_IO,   ONLY : IRSSQF, RSM_PUTIO

      IMPLICIT NONE
C
C     CHECK RISM CALCULATION
C
C
      CALL RSM_PARIN
      CALL RSM_OPDA(1)
      CALL RSM_SQOPN(IRSSQF,'RISM02','UNKNOWN','UNFORMATTED')
      CALL RSM_GETDAT0

      CALL RSM_PUTIO
      CALL RSM_SCF_PUTPAR
      CALL RSM_CLDA
      CALL RSM_SQCLS(IRSSQF)

      IF(IRISM /= 0) THEN
        RETURN
      ENDIF

      CALL RSMSTP

      RETURN
      END
C*MODULE RSM_STOP  *DECK RSM_SCF_STRT2
C>
C>    @brief   store the parameters (ex. coordinates) for RISM-SCF calculation
C>
C>    @author  Daisuke Yokogawa
C>
      SUBROUTINE RSM_SCF_STRT2

      USE RISM_IO,  ONLY : RSM_SETIO, IRSSQF
      USE RSMSED,   ONLY : IRISM, RSM_SCF_SETPAR

      IMPLICIT NONE
C
      CALL RSM_PARIN
      CALL RSM_OPDA(0)
      CALL RSM_SQOPN(IRSSQF,'RISM02','UNKNOWN','UNFORMATTED')
      CALL RSM_SETIO
      CALL RSM_SCF_SETPAR

      CALL RSM_CLDA
      CALL RSM_SQCLS(IRSSQF)

      IF(IRISM /= 0) THEN
        RETURN
      ENDIF

      CALL RSMSTP
C
      RETURN
      END
C*MODULE RSM_STOP  *DECK RSM_SCF_RHF
C>
C>    @brief   main subroutine for RISM-SCF calculation with RHF
C>
C>    @author  Daisuke Yokogawa
C>
      SUBROUTINE RSM_SCF_RHF(ITER,DIFF,FA,DA)

      USE RISM_IO,  ONLY : RSM_SETIO, IRSSQF
      USE RSMSED,   ONLY : IRISM, RSM_SCF_SETPAR

      IMPLICIT NONE

      INTEGER,          INTENT(IN   ) :: ITER
      DOUBLE PRECISION, INTENT(IN   ) :: DIFF, DA(*)
      DOUBLE PRECISION, INTENT(  OUT) :: FA(*)
C
C     SET MODULE DATA
C
C
      CALL RSM_PARIN
      CALL RSM_OPDA(0)
      CALL RSM_SQOPN(IRSSQF,'RISM02','UNKNOWN','UNFORMATTED')
      CALL RSM_SETIO
      CALL RSM_SCF_SETPAR

      CALL RSM_CLDA
      CALL RSM_SQCLS(IRSSQF)

      IF(IRISM /= 0) THEN
        RETURN
      ENDIF

      CALL RSMSTP

      RETURN
      END
C*MODULE RSM_STOP  *DECK RSM_SCF_UHF
C>
C>    @brief   main subroutine for RISM-SCF calculation with UHF
C>
C>    @author  Daisuke Yokogawa
C>
      SUBROUTINE RSM_SCF_UHF(ITER,DIFF,FA,FB,DA,DB,WRK)

      USE RISM_IO,  ONLY : RSM_SETIO, IRSSQF
      USE RSMSED,   ONLY : IRISM, RSM_SCF_SETPAR

      IMPLICIT NONE
C
      INTEGER,          INTENT(IN   ) :: ITER
      DOUBLE PRECISION, INTENT(IN   ) :: DIFF, DA(*), DB(*)
      DOUBLE PRECISION, INTENT(  OUT) :: FA(*), FB(*), WRK(*)
C
C     SET MODULE DATA
C
      CALL RSM_PARIN
      CALL RSM_OPDA(0)
      CALL RSM_SQOPN(IRSSQF,'RISM02','UNKNOWN','UNFORMATTED')
      CALL RSM_SETIO
      CALL RSM_SCF_SETPAR

      CALL RSM_CLDA
      CALL RSM_SQCLS(IRSSQF)

      IF(IRISM /= 0) THEN
        RETURN
      ENDIF

      CALL RSMSTP

      RETURN
      END
C*MODULE RSM_STOP  *DECK RSM_SCF_CC
C>
C>    @brief   main subroutine for RISM-SCF calculation with CC
C>
C>    @author  Daisuke Yokogawa
C>
      SUBROUTINE RSM_SCF_CC(ITER,NO,NU,IFC,FH,FP,FPH,OEH,OEP,T1)

      USE RISM_IO,  ONLY : RSM_SETIO, IRSSQF
      USE RSMSED,   ONLY : IRISM, RSM_SCF_SETPAR

      IMPLICIT NONE

      INTEGER,          INTENT(IN   ) :: ITER, NO, NU, IFC
      DOUBLE PRECISION, INTENT(INOUT) ::FH(*),FP(*),FPH(*),OEH(*),OEP(*)
      DOUBLE PRECISION, INTENT(  OUT) :: T1(*)
C
      CALL RSM_PARIN
      CALL RSM_OPDA(0)
      CALL RSM_SQOPN(IRSSQF,'RISM02','UNKNOWN','UNFORMATTED')
      CALL RSM_SETIO
      CALL RSM_SCF_SETPAR

      CALL RSM_CLDA
      CALL RSM_SQCLS(IRSSQF)

      IF(IRISM /= 0) THEN
        RETURN
      ENDIF

      CALL RSMSTP

      RETURN
      END
C*MODULE RSM_GMS_MAIN  *DECK RSM_SCF_CI
C>
C>    @brief   main subroutine for RISM-SCF calculation with CI
C>
C>    @author  Daisuke Yokogawa
C>
      SUBROUTINE RSM_SCF_CI(NRNFG,NPFLG)

      USE RISM_IO,  ONLY : RSM_SETIO, IRSSQF
      USE RSMSED,   ONLY : IRISM, RSM_SCF_SETPAR

      IMPLICIT NONE
C
      INTEGER, INTENT(INOUT) :: NRNFG(10), NPFLG(10)
C
C     RISM PART
C
      CALL RSM_PARIN
      CALL RSM_OPDA(0)
      CALL RSM_SQOPN(IRSSQF,'RISM02','UNKNOWN','UNFORMATTED')
      CALL RSM_SETIO
      CALL RSM_SCF_SETPAR

      CALL RSM_CLDA
      CALL RSM_SQCLS(IRSSQF)

      IF(IRISM /= 0) THEN
        RETURN
      ENDIF

      CALL RSMSTP

      END SUBROUTINE
C*MODULE RSM_STOP  *DECK RSM_SCF_MCQDPT
C>
C>    @brief   main subroutine for RISM-SCF calculation with MCQDPT
C>
C>    @author  Daisuke Yokogawa
C>
      SUBROUTINE RSM_SCF_MCQDPT(ITYP,LUNFT0,NDOUB,NMOACT,NSTATE,NCSF,
     *                             NMO,MAXCSF,EVAL,AVECOE)

      USE RISM_IO,  ONLY : RSM_SETIO, IRSSQF
      USE RSMSED,   ONLY : IRISM, RSM_SCF_SETPAR

      IMPLICIT NONE

      INTEGER,          INTENT(IN   ) :: ITYP,LUNFT0,NDOUB,NMOACT,
     *                                   NSTATE,NCSF,NMO,MAXCSF
      DOUBLE PRECISION, INTENT(IN   ) :: AVECOE(*)
      DOUBLE PRECISION, INTENT(  OUT) :: EVAL(*)
C
C     RISM PART
C
      CALL RSM_PARIN
      CALL RSM_OPDA(0)
      CALL RSM_SQOPN(IRSSQF,'RISM02','UNKNOWN','UNFORMATTED')
      CALL RSM_SETIO
      CALL RSM_SCF_SETPAR

      CALL RSM_CLDA
      CALL RSM_SQCLS(IRSSQF)

      IF(IRISM /= 0) THEN
        RETURN
      ENDIF

      CALL RSMSTP

      END SUBROUTINE
C*MODULE RSM_STOP  *DECK RSM_SCF_MCSCF1
C>
C>    @brief   main subroutine for RISM-SCF calculation with MCSCF
C>             (determine the memory size)
C>
C>    @author  Daisuke Yokogawa
C>
      SUBROUTINE RSM_SCF_MCSCF1(LOADFM,LAST,NEEDRSM,L2,L3,LD,LH1,LC,
     *                            MORBS,MCORBS)

      USE RISM_IO,  ONLY : RSM_SETIO, IRSSQF
      USE RSMSED,   ONLY : IRISM, RSM_SCF_SETPAR

      IMPLICIT NONE
C
      INTEGER, INTENT(IN   ) :: LOADFM, L2, L3, MORBS, MCORBS
      INTEGER, INTENT(  OUT) :: LD, LH1, LC, LAST, NEEDRSM
C
C     SET MODULE DATA
C
      CALL RSM_PARIN
      CALL RSM_OPDA(0)
      CALL RSM_SQOPN(IRSSQF,'RISM02','UNKNOWN','UNFORMATTED')
      CALL RSM_SETIO
      CALL RSM_SCF_SETPAR

      CALL RSM_CLDA
      CALL RSM_SQCLS(IRSSQF)

      RETURN
      END
C*MODULE RSM_STOP  *DECK RSM_SCF_MCSCF2
C>
C>    @brief   main subroutine for RISM-SCF calculation with MCSCF
C>
C>    @author  Daisuke Yokogawa
C>
      SUBROUTINE RSM_SCF_MCSCF2(MCDENAO,ITER,DIFF,FA,DA)

      USE RISM_IO,  ONLY : RSM_SETIO, IRSSQF
      USE RSMSED,   ONLY : IRISM, ITRRSM, RSM_SCF_SETPAR

      IMPLICIT NONE

      LOGICAL,          INTENT(  OUT) :: MCDENAO
      INTEGER,          INTENT(IN   ) :: ITER
      DOUBLE PRECISION, INTENT(IN   ) :: DA(*)
      DOUBLE PRECISION, INTENT(  OUT) :: DIFF, FA(*)
C
C     SET MODULE DATA
C
      CALL RSM_PARIN
      CALL RSM_OPDA(0)
      CALL RSM_SQOPN(IRSSQF,'RISM02','UNKNOWN','UNFORMATTED')
      CALL RSM_SETIO
      CALL RSM_SCF_SETPAR

      CALL RSM_CLDA
      CALL RSM_SQCLS(IRSSQF)

      IF(IRISM /= 0) THEN
        RETURN
      ENDIF

      CALL RSMSTP

      RETURN
      END
C*MODULE RSM_STOP  *DECK RSM_SCF_TD
C>
C>    @brief   main subroutine for RISM-SCF calculation with TD
C>
C>    @author  Daisuke Yokogawa
C>
C>
      SUBROUTINE RSM_SCF_TD

      USE RISM_IO,  ONLY : RSM_SETIO, IRSSQF
      USE RSMSED,   ONLY : IRISM, ITDDFT, RSM_SCF_SETPAR

      IMPLICIT NONE
C
C     RISM PART
C
      CALL RSM_PARIN
      CALL RSM_OPDA(0)
      CALL RSM_SQOPN(IRSSQF,'RISM02','UNKNOWN','UNFORMATTED')
      CALL RSM_SETIO
      CALL RSM_SCF_SETPAR

      CALL RSM_CLDA
      CALL RSM_SQCLS(IRSSQF)

      IF(IRISM /= 0 .or. ITDDFT == 0) RETURN

      CALL RSMSTP

      RETURN
      END
C*MODULE RSM_STOP  *DECK RSM_SCF_TRNSTN
C>
C>    @brief   main subroutine for RISM-SCF calculation with
C>             transition moments calculation
C>
C>    @author  Daisuke Yokogawa
C>
C>
      SUBROUTINE RSM_SCF_TRNSTN

      USE RISM_IO,  ONLY : RSM_SETIO
      USE RSMSED,   ONLY : IRISM, RSM_SCF_SETPAR

      IMPLICIT NONE
C
C
C     SET MODULE DATA
C
      CALL RSM_PARIN
      CALL RSM_OPDA(0)
      CALL RSM_SETIO
      CALL RSM_SCF_SETPAR

      CALL RSM_CLDA

      IF(IRISM /= 0) THEN
        RETURN
      ENDIF

      CALL RSMSTP

      RETURN
      END
C*MODULE RSM_STOP  *DECK RSM_SCFEND_HF
C>
C>    @brief   finalize RISM-SCF calculation with HF
C>
C>    @author  Daisuke Yokogawa
C>
      SUBROUTINE RSM_SCFEND_HF(ETOT,JOBTYP)

      USE RISM_IO,  ONLY : RSM_SETIO, IRSSQF
      USE RSMSED,   ONLY : IRISM, RSM_SCF_SETPAR

      IMPLICIT NONE

      INTEGER,          INTENT(IN   ) :: JOBTYP
      DOUBLE PRECISION, INTENT(INOUT) :: ETOT

      CALL RSM_PARIN
      CALL RSM_OPDA(0)
      CALL RSM_SQOPN(IRSSQF,'RISM02','UNKNOWN','UNFORMATTED')
      CALL RSM_SETIO
      CALL RSM_SCF_SETPAR

      CALL RSM_CLDA
      CALL RSM_SQCLS(IRSSQF)

      IF(IRISM /= 0) THEN
        RETURN
      ENDIF

      CALL RSMSTP

      RETURN
      END
C*MODULE RSM_STOP  *DECK RSM_SCFEND_CC
C>
C>    @brief   finalize RISM-SCF calculation with CC
C>
C>    @author  Daisuke Yokogawa
C>
      SUBROUTINE RSM_SCFEND_CC(ENRG)

      USE RISM_IO,  ONLY : RSM_SETIO, IRSSQF
      USE RSMSED,   ONLY : IRISM, RSM_SCF_SETPAR

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(IN   ) :: ENRG

      CALL RSM_PARIN
      CALL RSM_OPDA(0)
      CALL RSM_SQOPN(IRSSQF,'RISM02','UNKNOWN','UNFORMATTED')
      CALL RSM_SETIO
      CALL RSM_SCF_SETPAR

      CALL RSM_CLDA
      CALL RSM_SQCLS(IRSSQF)

      IF(IRISM /= 0) THEN
        RETURN
      ENDIF

      CALL RSMSTP

      RETURN
      END
C*MODULE RSM_STOP  *DECK RSM_SCFEND_MP2
C>
C>    @brief   finalize RISM-SCF calculation with MP2
C>
C>    @author  Daisuke Yokogawa
C>
      SUBROUTINE RSM_SCFEND_MP2(ETOT)

      USE RISM_IO,  ONLY : RSM_SETIO, IRSSQF
      USE RSMSED,   ONLY : IRISM, ICC, RSM_SCF_SETPAR
      USE RISM_MOD, ONLY : EUU, EUV_E8, RSM_SETENE, RSM_PUTENE, FRISM

      IMPLICIT NONE
C
      DOUBLE PRECISION, INTENT(IN   ) :: ETOT
C
      CALL RSM_PARIN
      CALL RSM_OPDA(0)
      CALL RSM_SQOPN(IRSSQF,'RISM02','UNKNOWN','UNFORMATTED')
      CALL RSM_SETIO
      CALL RSM_SCF_SETPAR

      CALL RSM_CLDA
      CALL RSM_SQCLS(IRSSQF)

      IF(IRISM /= 0) THEN
        RETURN
      ENDIF

      CALL RSMSTP

      RETURN
      END
C*MODULE RSM_STOP  *DECK RSM_SCFEND_TD
C>
C>    @brief   finalize RISM-SCF calculation with TD
C>
C>    @author  Daisuke Yokogawa
C>
      SUBROUTINE RSM_SCFEND_TD

      USE RISM_IO,  ONLY : RSM_SETIO, IRSSQF
      USE RSMSED,   ONLY : IRISM, RSM_SCF_SETPAR

      IMPLICIT NONE

      CALL RSM_PARIN
      CALL RSM_OPDA(0)
      CALL RSM_SQOPN(IRSSQF,'RISM02','UNKNOWN','UNFORMATTED')
      CALL RSM_SETIO
      CALL RSM_SCF_SETPAR

      CALL RSM_CLDA
      CALL RSM_SQCLS(IRSSQF)

      IF(IRISM /= 0) THEN
        RETURN
      ENDIF

      CALL RSMSTP

      RETURN
      END
C*MODULE RSM_STOP  *DECK RSM_SCFEND_MCQDPT2
C>
C>    @brief   finalize RISM-SCF calculation with MCQDPT2
C>
C>    @author  Daisuke Yokogawa
C>
      SUBROUTINE RSM_SCFEND_MCQDPT2(NSTATE,EMCQDPT2)

      USE RISM_IO,  ONLY : RSM_SETIO, IRSSQF
      USE RSMSED,   ONLY : IRISM, RSM_SCF_SETPAR

      IMPLICIT NONE

      INTEGER, INTENT(IN   ) :: NSTATE
      REAL*8,  INTENT(IN   ) :: EMCQDPT2(*)

      CALL RSM_PARIN
      CALL RSM_OPDA(0)
      CALL RSM_SQOPN(IRSSQF,'RISM02','UNKNOWN','UNFORMATTED')
      CALL RSM_SETIO
      CALL RSM_SCF_SETPAR

      CALL RSM_CLDA
      CALL RSM_SQCLS(IRSSQF)

      IF(IRISM /= 0) THEN
        RETURN
      ENDIF

      CALL RSMSTP

      RETURN
      END
C*MODULE RSM_STOP  *DECK RSM_SCF_MCSCF3
C>
C>    @brief   get energy of the (averaged) state in MCSCF iterations
C>
C>    @author  Daisuke Yokogawa
C>
      SUBROUTINE RSM_SCF_MCSCF3(ITER,MXRT,E,ESTATE)

      USE RISM_IO,  ONLY : RSM_SETIO, IRSSQF
      USE RSMSED,   ONLY : IRISM, RSM_SCF_SETPAR

      IMPLICIT NONE

      INTEGER,          INTENT(IN   ) :: ITER, MXRT
      DOUBLE PRECISION, INTENT(IN   ) :: ESTATE(*)
      DOUBLE PRECISION, INTENT(INOUT) :: E
C
      CALL RSM_PARIN
      CALL RSM_OPDA(0)
      CALL RSM_SQOPN(IRSSQF,'RISM02','UNKNOWN','UNFORMATTED')
      CALL RSM_SETIO
      CALL RSM_SCF_SETPAR

      CALL RSM_CLDA
      CALL RSM_SQCLS(IRSSQF)

      IF(IRISM /= 0) THEN
        RETURN
      ENDIF

      CALL RSMSTP

      RETURN
      END
C*MODULE RSM_STOP  *DECK RSM_SCFEND_MCSCF
C>
C>    @brief   finalize RISM-SCF calculation with MCSCF
C>
C>    @author  Daisuke Yokogawa
C>
      SUBROUTINE RSM_SCFEND_MCSCF(MASWRK,L1,L2,HCORE,ESTATE,NEEDRSM)

      USE RISM_IO,  ONLY : RSM_SETIO, IRSSQF
      USE RSMSED,   ONLY : IRISM, RSM_SCF_SETPAR

      IMPLICIT NONE

      LOGICAL,          INTENT(IN   ) :: MASWRK
      INTEGER,          INTENT(IN   ) :: L1, L2, NEEDRSM
      DOUBLE PRECISION, INTENT(INOUT) :: ESTATE(*), HCORE(*)
C
      CALL RSM_PARIN
      CALL RSM_OPDA(0)
      CALL RSM_SQOPN(IRSSQF,'RISM02','UNKNOWN','UNFORMATTED')
      CALL RSM_SETIO
      CALL RSM_SCF_SETPAR
C
      CALL RSM_CLDA
      CALL RSM_SQCLS(IRSSQF)
C
      RETURN
      END
C*MODULE RSM_STOP  *DECK RSM_SCFEND_TRNSTN
C>
C>    @brief   finalize RISM-SCF calculation with transition moments calculation
C>
C>    @author  Daisuke Yokogawa
C>
      SUBROUTINE RSM_SCFEND_TRNSTN

      USE RISM_IO,  ONLY : RSM_SETIO, IRSSQF
      USE RSMSED,   ONLY : IRISM, RSM_SCF_SETPAR

      IMPLICIT NONE

      CALL RSM_PARIN
      CALL RSM_OPDA(0)
      CALL RSM_SETIO
      CALL RSM_SCF_SETPAR

      CALL RSM_CLDA

      IF(IRISM /= 0) THEN
        RETURN
      ENDIF

      CALL RSMSTP

      RETURN
      END
C*MODULE RSM_STOP  *DECK RSM_SCFEND
C>
C>    @brief   finalize RISM-SCF calculation
C>
C>    @author  Daisuke Yokogawa
C>
      SUBROUTINE RSM_SCFEND

      USE RISM_IO,  ONLY : RSM_SETIO, IRSSQF
      USE RSMSED,   ONLY : IRISM, RSM_SCF_SETPAR

      IMPLICIT NONE
C
      CALL RSM_PARIN
      CALL RSM_OPDA(0)
      CALL RSM_SQOPN(IRSSQF,'RISM02','UNKNOWN','UNFORMATTED')
      CALL RSM_SETIO
      CALL RSM_SCF_SETPAR

      CALL RSM_CLDA
      CALL RSM_SQCLS(IRSSQF)

      IF(IRISM /= 0) THEN
        RETURN
      ENDIF

      CALL RSMSTP

      RETURN
      END
C*MODULE RSM_STOP  *DECK RSM_SCFGRD_1
C>
C>    @brief   turn off PROJGRAD
C>
C>    @author  Daisuke Yokogawa
C>
      SUBROUTINE RSM_SCFGRD_1(PROJGRAD)

      USE RISM_IO,  ONLY : RSM_SETIO, IRSSQF
      USE RSMSED,   ONLY : IRISM, RSM_SCF_SETPAR

      IMPLICIT NONE

      LOGICAL, INTENT(  OUT) :: PROJGRAD

      CALL RSM_PARIN
      CALL RSM_OPDA(0)
      CALL RSM_SETIO
      CALL RSM_SCF_SETPAR

      CALL RSM_CLDA

      IF(IRISM /= 0) THEN
        RETURN
      ENDIF

      CALL RSMSTP

      RETURN
      END
C*MODULE RSM_STOP  *DECK RSM_SCFGRD_2
C>
C>    @brief   compute gradients in solution
C>
C>    @author  Daisuke Yokogawa
C>
      SUBROUTINE RSM_SCFGRD_2

      USE RISM_IO,   ONLY : RSM_SETIO, IRSSQF
      USE RSMSED,    ONLY : IRISM, RSM_SCF_SETPAR

      IMPLICIT NONE
C
      INTEGER :: NEVALS, NGLEVL, NHLEVL
      DOUBLE PRECISION :: RUNTYP, EXETYP
      COMMON /RUNOPT/ RUNTYP,EXETYP,NEVALS,NGLEVL,NHLEVL
C
      DOUBLE PRECISION :: IMS, CHECK
      DATA IMS/8HIMS     /
      DATA CHECK  /8HCHECK   /
C
      IF(EXETYP == CHECK) RETURN
C
C     SET MODULE DATA
C
      CALL RSM_PARIN
      CALL RSM_OPDA(0)
      CALL RSM_SQOPN(IRSSQF,'RISM02','UNKNOWN','UNFORMATTED')
      CALL RSM_SETIO
      CALL RSM_SCF_SETPAR

      CALL RSM_CLDA
      CALL RSM_SQCLS(IRSSQF)

      IF(IRISM /= 0) THEN
        RETURN
      ENDIF

      CALL RSMSTP

      RETURN
      END
C*MODULE RSM_STOP  *DECK RSM_SETIRISM
C>
C>    @brief   set IRISM
C>
C>    @author  Daisuke Yokogawa
C>
C>    @param   IRISM_GMS: 0: RISM-SCF calculation, 1: gas phase
C>
      SUBROUTINE RSM_SETIRISM(IRISM_GMS)

      USE RSMSED,   ONLY : IRISM, RSM_SCF_SETPAR

      IMPLICIT NONE
C
      INTEGER, INTENT(  OUT) :: IRISM_GMS

      CALL RSM_PARIN
      CALL RSM_OPDA(0)
      CALL RSM_SCF_SETPAR

      IRISM_GMS = IRISM

      CALL RSM_CLDA

      RETURN
      END
C*MODULE RSM_STOP  *DECK RSM_CCINI
C>
C>    @brief   store the hole and particle orbital energies determied in gas
C>             phase in OEH_ORG and OEP_ORG, respectively
C>
C>    @author  Daisuke Yokogawa
C>
      SUBROUTINE RSM_CCINI(NO,NU,OEH,OEP)

      USE RSMSED,    ONLY : IRISM
      USE RSMCC_MOD, ONLY : OEH_ORG, OEP_ORG

      IMPLICIT NONE
C
      INTEGER,          INTENT(IN   ) :: NO, NU
      DOUBLE PRECISION, INTENT(IN   ) :: OEH(*), OEP(*)
C
      RETURN
      END SUBROUTINE
C*MODULE RSM_STOP  *DECK RSM_CCFIN
C>
C>    @brief   return the hole and particle orbital energies stored in
C>             OEH_ORG and OEP_ORG
C>
C>    @author  Daisuke Yokogawa
C>
      SUBROUTINE RSM_CCFIN(NO,NU,OEH,OEP)

      USE RSMSED,    ONLY : IRISM
      USE RSMCC_MOD, ONLY : OEH_ORG, OEP_ORG

      IMPLICIT NONE
C
      INTEGER,          INTENT(IN   ) :: NO, NU
      DOUBLE PRECISION, INTENT(  OUT) :: OEH(*), OEP(*)
C
      RETURN
      END SUBROUTINE
C*MODULE RSM_STOP  *DECK RSM_CPHF3
C>
C>    @brief   modify WAX and AA in CPHF calculations and complete F^{SOLV}
C>
C>    @author  Daisuke Yokogawa
C>
      SUBROUTINE RSM_CPHF3(NOCC,NVIR,IA,WAX,AA)

      USE RISM_IO,  ONLY : RSM_SETIO, IRSSQF
      USE RSMSED,   ONLY : IRISM, RSM_SCF_SETPAR

      IMPLICIT NONE

      INTEGER,          INTENT(IN   ) :: NOCC, NVIR, IA(*)
      DOUBLE PRECISION, INTENT(INOUT) :: WAX(NOCC,NVIR,*),
     *                                   AA((NOCC*NOCC+NOCC)/2,*)
C
C     SET MODULE DATA
C
      CALL RSM_PARIN
      CALL RSM_OPDA(0)
      CALL RSM_SQOPN(IRSSQF,'RISM02','UNKNOWN','UNFORMATTED')
      CALL RSM_SETIO
      CALL RSM_SCF_SETPAR

      CALL RSM_CLDA
      CALL RSM_SQCLS(IRSSQF)

      IF(IRISM /= 0) THEN
        RETURN
      ENDIF

      CALL RSMSTP

      RETURN
      END
C*MODULE RSM_STOP  *DECK RSM_CPHF4
C>
C>    @brief   compute tQ * v * Q * (trial vector)
C>             (see: the last term in eq 20 in JCP, 155, 204102 (2021))
C>
C>    @author  Daisuke Yokogawa
C>
      SUBROUTINE RSM_CPHF4(ITYP,PMN,FMN,NXYZF,ERR,ICONV)

      USE RISM_IO,  ONLY : RSM_SETIO, IRSSQF
      USE RSMSED,   ONLY : IRISM, RSM_SCF_SETPAR

      IMPLICIT NONE
C
      INTEGER,          INTENT(IN   ) :: ITYP, NXYZF
      INTEGER,          INTENT(  OUT) :: ICONV
      DOUBLE PRECISION, INTENT(IN   ) :: PMN(*), ERR
      DOUBLE PRECISION, INTENT(  OUT) :: FMN(*)
C
C     SET MODULE DATA
C
      CALL RSM_PARIN
      CALL RSM_OPDA(0)
      CALL RSM_SQOPN(IRSSQF,'RISM02','UNKNOWN','UNFORMATTED')
      CALL RSM_SETIO
      CALL RSM_SCF_SETPAR

      CALL RSM_CLDA
      CALL RSM_SQCLS(IRSSQF)

      IF(IRISM /= 0) THEN
        RETURN
      ENDIF

      CALL RSMSTP

      RETURN
      END
C*MODULE RSM_STOP  *DECK RSM_CPHF5
C>
C>    @brief   add the solvent contribution to EG and FCM
C>
C>    @author  Daisuke Yokogawa
C>
      SUBROUTINE RSM_CPHF5(NAT,FCM_GMS)

      USE RISM_IO,  ONLY : RSM_SETIO, IRSSQF
      USE RSMSED,   ONLY : IRISM, RSM_SCF_SETPAR

      IMPLICIT NONE
C
      INTEGER,          INTENT(IN   ) :: NAT
      DOUBLE PRECISION, INTENT(INOUT) :: FCM_GMS(*)
C
C     SET MODULE DATA
C
      CALL RSM_PARIN
      CALL RSM_OPDA(0)
      CALL RSM_SQOPN(IRSSQF,'RISM02','UNKNOWN','UNFORMATTED')
      CALL RSM_SETIO
      CALL RSM_SCF_SETPAR

      CALL RSM_CLDA
      CALL RSM_SQCLS(IRSSQF)

      IF(IRISM /= 0) THEN
        RETURN
      ENDIF

      CALL RSMSTP

      RETURN
      END
C*MODULE RSM_STOP  *DECK RSM_CPHF7
C>
C>    @brief   read and write the response vector
C>
C>    @author  Daisuke Yokogawa
C>
      SUBROUTINE RSM_CPHF7(MODE,NNXYZ,NOCC,NVIR,NROT,YA,VEC1)

      USE RISM_IO,  ONLY : RSM_SETIO, IRSSQF
      USE RSMSED,   ONLY : IRISM, RSM_SCF_SETPAR

      IMPLICIT NONE

      INTEGER,          INTENT(IN   ) :: MODE, NNXYZ, NOCC, NVIR, NROT
      DOUBLE PRECISION, INTENT(INOUT) :: YA(*), VEC1(*)
C
C     SET MODULE DATA
C
      CALL RSM_PARIN
      CALL RSM_OPDA(0)
      CALL RSM_SQOPN(IRSSQF,'RISM02','UNKNOWN','UNFORMATTED')
      CALL RSM_SETIO
      CALL RSM_SCF_SETPAR

      CALL RSM_CLDA
      CALL RSM_SQCLS(IRSSQF)

      IF(IRISM /= 0) THEN
        RETURN
      ENDIF

      CALL RSMSTP

      RETURN
      END
C*MODULE RSM_STOP  *DECK RSM_CPHF8
C>
C>    @brief   modify the L.H.S. of Z-vector equation
C>
C>    @author  Daisuke Yokogawa
C>
      SUBROUTINE RSM_CPHF8(PMN,XLAI,CMO,NOCC,NVIR,NBF)

      USE RISM_IO,  ONLY : RSM_SETIO, IRSSQF
      USE RSMSED,   ONLY : IRISM, RSM_SCF_SETPAR

      IMPLICIT NONE

      INTEGER,          INTENT(IN   ) :: NOCC, NVIR, NBF
      DOUBLE PRECISION, INTENT(IN   ) :: PMN(*), CMO(NBF,*)
      DOUBLE PRECISION, INTENT(INOUT) :: XLAI(NOCC,NVIR)
C
C     SET MODULE DATA
C
      CALL RSM_PARIN
      CALL RSM_OPDA(0)
      CALL RSM_SQOPN(IRSSQF,'RISM02','UNKNOWN','UNFORMATTED')
      CALL RSM_SETIO
      CALL RSM_SCF_SETPAR

      CALL RSM_CLDA
      CALL RSM_SQCLS(IRSSQF)

      IF(IRISM /= 0) THEN
        RETURN
      ENDIF

      CALL RSMSTP

      RETURN
      END
C*MODULE RSM_STOP  *DECK RSM_CPHF9
C>
C>    @brief   modify the A matrix in the R.H.S. of Z-vector equation
C>
C>    @author  Daisuke Yokogawa
C>
      SUBROUTINE RSM_CPHF9(GMVEC0,GMVEC1,NOCC,NVIR,NBF,CMO,ENG,FACTOR)

      USE RISM_IO,  ONLY : RSM_SETIO, IRSSQF
      USE RSMSED,   ONLY : IRISM, RSM_SCF_SETPAR

      IMPLICIT NONE
C
      INTEGER,          INTENT(IN   ) :: NOCC, NVIR, NBF
      DOUBLE PRECISION, INTENT(IN   ) :: GMVEC0(NOCC,NVIR),CMO(NBF,NBF),
     *                                   ENG(*),FACTOR
      DOUBLE PRECISION, INTENT(INOUT) :: GMVEC1(NOCC,NVIR)
C
C     SET MODULE DATA
C
      CALL RSM_PARIN
      CALL RSM_OPDA(0)
      CALL RSM_SQOPN(IRSSQF,'RISM02','UNKNOWN','UNFORMATTED')
      CALL RSM_SETIO
      CALL RSM_SCF_SETPAR

      CALL RSM_CLDA
      CALL RSM_SQCLS(IRSSQF)

      IF(IRISM /= 0) THEN
        RETURN
      ENDIF

      CALL RSMSTP

      RETURN
      END
C*MODULE RSM_STOP  *DECK RSM_CPHF10
C>
C>    @brief   modify W_{ij}[2] for gradient calculation
C>
C>    @author  Daisuke Yokogawa
C>
      SUBROUTINE RSM_CPHF10(PMN,WIJ,NOCC,NBF,CMO)

      USE RISM_IO,  ONLY : RSM_SETIO, IRSSQF
      USE RSMSED,   ONLY : IRISM, RSM_SCF_SETPAR

      IMPLICIT NONE
C
      INTEGER,          INTENT(IN   ) :: NOCC, NBF
      DOUBLE PRECISION, INTENT(IN   ) :: PMN(NBF,NBF), CMO(NBF*NBF)
      DOUBLE PRECISION, INTENT(INOUT) :: WIJ(NOCC,*)
C
C     SET MODULE DATA
C
      CALL RSM_PARIN
      CALL RSM_OPDA(0)
      CALL RSM_SQOPN(IRSSQF,'RISM02','UNKNOWN','UNFORMATTED')
      CALL RSM_SETIO
      CALL RSM_SCF_SETPAR

      CALL RSM_CLDA
      CALL RSM_SQCLS(IRSSQF)

      IF(IRISM /= 0) THEN
        RETURN
      ENDIF

      CALL RSMSTP

      RETURN
      END
C*MODULE RSM_STOP  *DECK RSM_NMR_DRV
C>
C>    @brief   prepare dY/dB
C>             (see: eq 34 in JCP, 152, 194102 (2020))
C>
C>    @author  Kosuke Imamura
C>
      SUBROUTINE RSM_NMR_DRV

      USE RISM_IO,  ONLY : RSM_SETIO, IRSSQF
      USE RSMSED,   ONLY : IRISM, RSM_SCF_SETPAR

      IMPLICIT NONE
C
C     SET MODULE DATA
C
      CALL RSM_PARIN
      CALL RSM_OPDA(0)
      CALL RSM_SQOPN(IRSSQF,'RISM02','UNKNOWN','UNFORMATTED')
      CALL RSM_SETIO
      CALL RSM_SCF_SETPAR

      CALL RSM_CLDA
      CALL RSM_SQCLS(IRSSQF)

      IF(IRISM /= 0) THEN
        RETURN
      ENDIF

      CALL RSMSTP

      RETURN
      END
C*MODULE RSM_STOP  *DECK RSM_NMR_1
C>
C>    @brief   fetch dS/dB
C>
C>    @author  Kosuke Imamura
C>
      SUBROUTINE RSM_NMR_1(NCART,S10,L2S)

      USE RISM_IO,  ONLY : RSM_SETIO, IRSSQF
      USE RSMSED,   ONLY : IRISM, RSM_SCF_SETPAR

      IMPLICIT NONE
C
      INTEGER,          INTENT(IN   ) :: NCART, L2S
      DOUBLE PRECISION, INTENT(IN   ) :: S10(*)
C
C     SET MODULE DATA
C
      CALL RSM_PARIN
      CALL RSM_OPDA(0)
      CALL RSM_SQOPN(IRSSQF,'RISM02','UNKNOWN','UNFORMATTED')
      CALL RSM_SETIO
      CALL RSM_SCF_SETPAR

      CALL RSM_CLDA
      CALL RSM_SQCLS(IRSSQF)

      IF(IRISM /= 0) THEN
        RETURN
      ENDIF

      CALL RSMSTP

      RETURN
      END
C*MODULE RMS_STOP  *DECK RSM_NMR_2
C>
C>    @brief   fetch dF/dB in the gas-phase and add solvation term
C>
C>    @author  Kosuke Imamura
C>
      SUBROUTINE RSM_NMR_2(F10,L2S)

      USE RISM_IO,  ONLY : RSM_SETIO, IRSSQF
      USE RSMSED,   ONLY : IRISM, RSM_SCF_SETPAR

      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(INOUT) :: F10(*)
      INTEGER,          INTENT(IN   ) :: L2S
C
C     SET MODULE DATA
C
      CALL RSM_PARIN
      CALL RSM_OPDA(0)
      CALL RSM_SQOPN(IRSSQF,'RISM02','UNKNOWN','UNFORMATTED')
      CALL RSM_SETIO
      CALL RSM_SCF_SETPAR

      CALL RSM_CLDA
      CALL RSM_SQCLS(IRSSQF)

      IF(IRISM /= 0) THEN
        RETURN
      ENDIF

      CALL RSMSTP

      RETURN
      END
C*MODULE RSM_STOP  *DECK RSM_GETDAT4
C>
C>    @brief   switch the TD-DFT conditions
C>
C>    @author  Daisuke Yokogawa
C>
      SUBROUTINE RSM_GETDAT4(NRUN)
C
      USE RISM_IO,   ONLY : RSM_SETIO, IRSSQF
      USE RSMSED,    ONLY : IRISM, RSM_SCF_SETPAR
C
      IMPLICIT NONE
C
      INTEGER, INTENT(IN   ) :: NRUN
C
C     CHECK RISM CALCULATION
C
      CALL RSM_PARIN
      CALL RSM_OPDA(0)
      CALL RSM_SQOPN(IRSSQF,'RISM02','UNKNOWN','UNFORMATTED')
      CALL RSM_SETIO
      CALL RSM_SCF_SETPAR

      CALL RSM_CLDA
      CALL RSM_SQCLS(IRSSQF)

      IF(IRISM /= 0) THEN
        RETURN
      ENDIF

      CALL RSMSTP

      RETURN
      END
C*MODULE RSM_STOP  *DECK RSM_MODENG
C>
C>    @brief   modify the nergies for the conical intersection search
C>
C>    @author  Daisuke Yokogawa
C>
      SUBROUTINE RSM_MODENG(IXSTAT,NAT,TDDFTYP,ENG1,ENG2,EG1,EG2)
C
      USE RISM_IO,   ONLY : RSM_SETIO, IRSSQF
      USE RSMSED,    ONLY : IRISM, RSM_SCF_SETPAR
C
      IMPLICIT NONE
C
      INTEGER,          INTENT(IN   ) :: IXSTAT(*), NAT
      DOUBLE PRECISION, INTENT(IN   ) :: TDDFTYP
      DOUBLE PRECISION, INTENT(  OUT) :: ENG1, ENG2, EG1(*), EG2(*)
C
C     CHECK RISM CALCULATION
C
      CALL RSM_PARIN
      CALL RSM_OPDA(0)
      CALL RSM_SQOPN(IRSSQF,'RISM02','UNKNOWN','UNFORMATTED')
      CALL RSM_SETIO
      CALL RSM_SCF_SETPAR

      CALL RSM_CLDA
      CALL RSM_SQCLS(IRSSQF)

      IF(IRISM /= 0) THEN
        RETURN
      ENDIF

      CALL RSMSTP

      RETURN
      END
C*MODULE RSM_STOP  *DECK RSM_CONIADDLOOP
C>
C>    @brief   modify the the conical intersection search loop
C>
C>    @author  Daisuke Yokogawa
C>
      SUBROUTINE RSM_CONIADDLOOP(EXETYP,METHOD,NAT,NPRT,NPRTHS,NEG2CL,
     *                        IXROOT)

      USE RISM_IO,   ONLY : RSM_SETIO, IRSSQF
      USE RSMSED,    ONLY : IRISM, RSM_SCF_SETPAR

      IMPLICIT NONE
C
      DOUBLE PRECISION, INTENT(IN   ) :: EXETYP, METHOD
      INTEGER,          INTENT(IN   ) :: NAT, NPRT, NPRTHS, IXROOT(2)
      INTEGER,          INTENT(INOUT) :: NEG2CL
C
C
C     CHECK RISM CALCULATION
C
      CALL RSM_PARIN
      CALL RSM_OPDA(0)
      CALL RSM_SQOPN(IRSSQF,'RISM02','UNKNOWN','UNFORMATTED')
      CALL RSM_SETIO
      CALL RSM_SCF_SETPAR

      CALL RSM_CLDA
      CALL RSM_SQCLS(IRSSQF)
      IF(IRISM /= 0) RETURN

      CALL RSMSTP

      RETURN
      END
