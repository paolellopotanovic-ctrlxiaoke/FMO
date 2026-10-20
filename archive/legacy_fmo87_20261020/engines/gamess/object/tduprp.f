C*MODULE TDUPRP *DECK TDUPRP
C>
C>    @brief   Calculate properties for unrelaxed R-TDDFT densities
C>
C>    @details Inspired by property calculations in TDGRAD and also
C>             MULKEN. Here, the unrelaxed densities are used.
C>
C>    @author Christian Friedl
C>
C>    @date   Apr, 2022 - Christian Friedl
C>    - Added subroutine
C>    @date   Nov, 2022 - Christian Friedl
C>    - Memory is always initialized before usage
C>    - Zero FED coupling is handled properly
C>    - Logicals are used instead of bitfields
C>    - NPRINT=-5 is handled properly
C>    @date   Nov, 2022 - Christian Friedl
C>    - Enabled FCD printing
C>    - Added Multi-FED-FCD
C>
      SUBROUTINE TDUPRP
C
      USE MX_LIMITS, ONLY: MXATM,MXRT,MXAO
      USE CONSTANTS, ONLY: ZERO,ONE,TWO,THREE
      IMPLICIT NONE
C     Declarations for common blocks
      DOUBLE PRECISION QMTTOL,X,ZAN,C,
     *                 CNVTOL,PFREQ,SPCP,
     *                 RUNTYP,EXETYP,espscf,e0scf,emp2s
      INTEGER ISPHER,IA,NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,IAN,
     *        MODTD,JANST,NRADT,NTHET,NPHIT,NLEBT,NSTAT,NTRIAL,MAXVEC,
     *        NTHST,IRECTD,ITDFG,ITDPRP,NONEQR,MULTD,MTHST,IFEDAT,
     *        IR,IW,IP,IS,IPK,IDAF,NAV,IODA,
     *        NPRINT,ITOL,ICUT,NORMF,NORMP,NOPK,
     *        ME,MASTER,NPROC,IBTYP,IPTIM,
     *        NEVALS,NGLEVL,NHLEVL,nfg,nlayer,natfmo,nbdfg,naotyp,nbody,
     *        IDAFMO,icurfg,jcurfg,kcurfg,icurlay,icurunt,nat1e,ncursh,
     *        ngau,icurpop,ifmostp,moncor,needr,modrst,norbproj,nunesp,
     *        iskipesp,IESDPPC,idoprop,mp2run,icurit,idmfmo,iddfmo,
     *        iddcur,nddleft,ivmfmo,nzmtfmo,ifmobas,itmfmo,nsegm
      LOGICAL DFTBFL,SCC,SRSCC,DFTB3,DAMPXH,LCDFTB
      LOGICAL TRIPLET,SG1T,TAMMD,TPA,ALPHKWD,BETAKWD
      LOGICAL MREKT,MRDEA
      LOGICAL GOPARR,DSKWRK,MASWRK
C     Other declarations
      DOUBLE PRECISION, PARAMETER :: HALF = 0.5D+00
      DOUBLE PRECISION, PARAMETER :: SQRT_TWO_HALF = SQRT(TWO)*HALF
      DOUBLE PRECISION, PARAMETER :: TOEV = 27.21138386D+00
      DOUBLE PRECISION CHECK,X12D,X12A,DX12,DX1,DX2
      DOUBLE PRECISION E1,E2,COUPL,EL1,EL2,TANTH,TAN2TH
      DOUBLE PRECISION SINSQTH,COSSQTH
C     DOUBLE PRECISION SINCOSTH
      DOUBLE PRECISION Q12D,Q12A,DQ1,DQ2,DQ12
      DOUBLE PRECISION Q0D,Q0A,DQ0,Q02D,Q02A,DQ02
      DOUBLE PRECISION DX,DY,DZ,OS
      INTEGER NDSR,LX,NOCC,NVIR,L1,L2,L3,L7,L7MAX
      INTEGER LAST,LOADFM,NEED,NEED1,NEED2
      INTEGER LINXOV,LVROALL,LVLOALL,LTDEN,LDEN,LZAN,LV,LSCR
      INTEGER LOVLP,LMULPOP,LMULCH,LWRK1,LWRK2,LWRK3,LSHMAT
      INTEGER LX12D,LX12A,LQ12D,LQ12A,LEE,LQ02D,LQ02A
      INTEGER LBLE1,LBLE2,LBCT1,LBCT2,LHLE1,LHLE2,LHCT1,LHCT2
      INTEGER LBLE1D,LBLE2D,LBCT1D,LBCT2D,LHFIN
      INTEGER LVRO,LVLO,LXMAT,LQMAT,LDMAT,LW
      INTEGER LOS,LSCR2,LAX,LAY,LAZ,LTXYZ,LTXYZ2
      INTEGER LTXYZMF,LTXYZ2MF,LX12DL,LX12AL,LQ12DL,LQ12AL
      INTEGER LDONI,LACCI
      INTEGER I,J,ITMP,IST,JST,ITF,NTHSTSV
      INTEGER IFINV,IFTSP,INFO,ILECT,ILE2LE1,ICT2CT1
      INTEGER ICURFGST,JCURFGST,NDON,NACC
      LOGICAL TDPRP,TRNSD,UNRLX,ALTRN,FED,DOMULTIST
      INTEGER, EXTERNAL :: IXFTCH
C
      COMMON /BASSPH/ QMTTOL,ISPHER
      COMMON /DFTB  / DFTBFL,SCC,SRSCC,DFTB3,DAMPXH,LCDFTB
      COMMON /FMCOM / X(1)
      common /fmoinf/ nfg,nlayer,natfmo,nbdfg,naotyp,nbody,nsegm
      common /fmorun/ espscf,e0scf(2),emp2s,IDAFMO,icurfg,jcurfg,kcurfg,
     *                icurlay,icurunt,nat1e,ncursh,ngau,icurpop,ifmostp,
     *                moncor,needr,modrst,norbproj,nunesp,iskipesp,
     *                IESDPPC,idoprop,mp2run,icurit,idmfmo,iddfmo,
     *                iddcur,nddleft,ivmfmo,nzmtfmo,ifmobas,itmfmo(2)
      COMMON /IJPAIR/ IA(MXAO)
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /INFOTD/ CNVTOL,PFREQ(2),MODTD,
     *                JANST,NRADT,NTHET,NPHIT,NLEBT,
     *                NSTAT,NTRIAL,MAXVEC,NTHST,IRECTD,ITDFG,ITDPRP,
     *                TRIPLET,SG1T,NONEQR,TAMMD,TPA,ALPHKWD,BETAKWD,
     *                SPCP(3),MULTD,MREKT,MRDEA,MTHST,IFEDAT(4)
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON /OUTPUT/ NPRINT,ITOL,ICUT,NORMF,NORMP,NOPK
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /RUNOPT/ RUNTYP,EXETYP,NEVALS,NGLEVL,NHLEVL
C
      DATA CHECK/8HCHECK   /
C
      IF(TRIPLET) THEN
        IF(MASWRK) THEN
          WRITE(IW,*) 'TRIPLET IS NOT SUPPORTED FOR TDDFT PROPERTIES'
        END IF
        CALL ABRT
      END IF
C
      IF (NTHST.GT.MXRT) THEN
        IF(MASWRK) THEN
          WRITE(IW,*) 'TDDFT PROPERTY CALCULATION IS ALLOWED ONLY',
     *                'FOR IROOT.LE.',MXRT
        END IF
        CALL ABRT
      END IF
      IF(NPRINT.EQ.-5) CALL AOLIM
C
C     --- DEFINE LOGICALS
C
      TDPRP = IAND(ITDPRP,1).NE.0
      TRNSD = IAND(ITDPRP,2).NE.0
      UNRLX = IAND(ITDPRP,4).NE.0
      ALTRN = IAND(ITDPRP,8).NE.0
C     In FMO, do FED only for dimers.
      FED = (IAND(ITDPRP,16).NE.0.AND.NFG.EQ.0.).OR.
     *    (IAND(ITDPRP,16).NE.0.AND.NFG.NE.0.AND.IFMOSTP.EQ.4)
C     This is set by subroutine FEDFMOMATCH
      DOMULTIST = .FALSE.
C
C     --- DEFINE NDSR
C
      NDSR = NSTAT
C
C     --- DEFINE NOCC & NVIR
C
      LX = NQMT
      NOCC = NA
      NVIR = LX-NOCC
      L1 = NUM
      L2 = (L1*(L1+1))/2
      L3 = L1*L1
      L7 = NOCC*NVIR
C
C     --- DEFINE ITF
C
      ITF = 0
      IF(TAMMD) ITF = 1
C
      IF(TDPRP) THEN
         CALL VALFM(LOADFM)
         LINXOV  = LOADFM + 1
         LVROALL = LINXOV + L7 * 2
         LVLOALL = LVROALL + L7*NDSR
         LTDEN   = LVLOALL + L7*NDSR
         LDEN    = LTDEN + L1*L1
         LZAN    = LDEN + L2
         LV      = LZAN + NAT
         LSCR    = LV + L1*L1
         LAST    = LSCR + L1*LX
         NEED    = LAST - LOADFM - 1
         CALL GETFM(NEED)
C
         CALL VCLR(X(LINXOV), 1,L7*2)
         CALL VCLR(X(LVROALL),1,L7*NDSR)
         CALL VCLR(X(LVLOALL),1,L7*NDSR)
         CALL VCLR(X(LTDEN),  1,L1*L1)
         CALL VCLR(X(LDEN),   1,L2)
         CALL VCLR(X(LZAN),   1,NAT)
         CALL VCLR(X(LV),     1,L1*L1)
         CALL VCLR(X(LSCR),   1,L1*LX)
C
         IF(EXETYP.NE.CHECK.AND.NPRINT.NE.-5) THEN
            CALL DAREAD(IDAF,IODA,X(LVROALL),L7*NDSR,471,0)
            CALL DAREAD(IDAF,IODA,X(LVLOALL),L7*NDSR,479,0)
            CALL DAREAD(IDAF,IODA,X(LV),L1*LX,15,0)
C
            CALL DAREAD(IDAF,IODA,X(LDEN),L2,16,0)
            NTHSTSV = NTHST
            CALL DCOPY(NAT,ZAN,1,X(LZAN),1)
C
C           --- UNRELAXED EXCITED STATE DENSITY
C
            DO IST=1,NDSR
               IF(IST.EQ.NTHSTSV.OR.ALTRN) THEN
                  CALL TCONST2(X(LTDEN),X(LVROALL),X(LVLOALL),
     *                         LX,L7,NOCC,ITF,IST,IST)
                  CALL DGEMM('N','N',L1,LX,LX,ONE,X(LV),L1,X(LTDEN),LX,
     *                        ZERO,X(LSCR),L1)
                  CALL DGEMM('N','T',L1,L1,LX,ONE,X(LSCR),L1,X(LV),L1,
     *                       ZERO,X(LTDEN),L1)
                  CALL SQ2TRI(L1,L1,X(LTDEN),X(LSCR),HALF)
                  CALL DAXPY(L2,ONE,X(LDEN),1,X(LSCR),1)
                  CALL DAWRIT(IDAF,IODA,X(LSCR),L2,16,0)
                  NTHST = IST
                  CALL PROPTY('TDD3')
                  NTHST = NTHSTSV
                  CALL DAWRIT(IDAF,IODA,X(LDEN),L2,16,0)
               END IF
            END DO
C
C           --- GROUND TO EXCITED STATE TRANSITION DENSITY
C
            IF(TRNSD) THEN
               CALL EXCnstLab(X(LINXOV),NOCC,LX)
               DO IST=1,NDSR
                  IF(IST.EQ.NTHSTSV.OR.ALTRN) THEN
                     LVRO = LVROALL+(IST-1)*L7
                     LVLO = LVLOALL+(IST-1)*L7
                     CALL TRAD(X(LTDEN),X(LVRO),X(LVLO),LX,L7,
     *                         X(LINXOV),1,1)
C                    possibly better (check!):
C                    CALL TRAD(X(LTDEN),X(LVROALL),X(LVLOALL),LX,L7,
C    *                         X(LINXOV),NDSR,IST)
                     CALL DSCAL(L1*L1,SQRT_TWO_HALF,X(LTDEN),1)
                     CALL DGEMM('N','N',L1,LX,LX,ONE,X(LV),L1,
     *                          X(LTDEN),LX,ZERO,X(LSCR),L1)
                     CALL DGEMM('N','T',L1,L1,LX,ONE,X(LSCR),L1,
     *                          X(LV),L1,ZERO,X(LTDEN),L1)
                     CALL SQ2TRI(L1,L1,X(LTDEN),X(LSCR),HALF)
                     CALL DAWRIT(IDAF,IODA,X(LSCR),L2,16,0)
                     CALL VCLR(ZAN,1,NAT)
                     NTHST = IST
                     CALL PROPTY('TDD2')
                     NTHST = NTHSTSV
                     CALL DCOPY(NAT,X(LZAN),1,ZAN,1)
                     CALL DAWRIT(IDAF,IODA,X(LDEN),L2,16,0)
                  END IF
               END DO
            END IF
C
C           --- UNRELAXED EXCITED TO EXCITED STATE TRANSITION DENSITY
C
            IF(ALTRN) THEN
               DO IST=1,NDSR
                  DO JST=IST+1,NDSR
                     CALL TCONST2(X(LTDEN),X(LVROALL),X(LVLOALL),
     *                          LX,L7,NOCC,ITF,IST,JST)
                     CALL DGEMM('N','N',L1,LX,LX,ONE,X(LV),
     *                          L1,X(LTDEN),LX,ZERO,X(LSCR),L1)
                     CALL DGEMM('N','T',L1,L1,LX,ONE,X(LSCR),
     *                          L1,X(LV),L1,ZERO,X(LTDEN),L1)
                     CALL SQ2TRI(L1,L1,X(LTDEN),X(LSCR),HALF)
                     CALL DAWRIT(IDAF,IODA,X(LSCR),L2,16,0)
                     CALL VCLR(ZAN,1,NAT)
                     NTHST = IST
                     MTHST = JST
                     CALL PROPTY('TDD4')
                     NTHST = NTHSTSV
                     CALL DCOPY(NAT,X(LZAN),1,ZAN,1)
                     CALL DAWRIT(IDAF,IODA,X(LDEN),L2,16,0)
                  END DO
               END DO
            END IF
         END IF
C
         CALL RETFM(NEED)
      END IF
C
      IF(FED) THEN
C
C        --- ALLOCATE MEMORY FOR FED/FCD RESULTS ---
C
         CALL VALFM(LOADFM)
         LX12D   = LOADFM +1
         LX12A   = LX12D + NDSR*NDSR
         LEE     = LX12A + NDSR*NDSR
         LQ12D   = LEE + NDSR
         LQ12A   = LQ12D + NDSR*NDSR
         LQ02D   = LQ12A + NDSR*NDSR
         LQ02A   = LQ02D + NDSR
         LTXYZ   = LQ02A + NDSR
         LTXYZ2  = LTXYZ + 3*NDSR
         LAST    = LTXYZ2 + 3*NDSR*NDSR
         NEED    = LAST - LOADFM - 1
         CALL GETFM(NEED)
C
C        --- ALLOCATE WORKING MEMORY FOR FED/FCD ---
C
         CALL VALFM(LOADFM)
         LINXOV  = LOADFM + 1
         LVROALL = LINXOV + L7 * 2
         LVLOALL = LVROALL + L7*NDSR
         LTDEN   = LVLOALL + L7*NDSR
         LDEN    = LTDEN + L1*L1
         LV      = LDEN + L2
         LSCR    = LV + L1*L1
         LOVLP   = LSCR + L1*LX
         LMULPOP = LOVLP + L2
         LMULCH  = LMULPOP + L2
         LWRK1   = LMULCH + NAT
         LWRK2   = LWRK1 + L2
         LSHMAT  = LWRK2 + L1
         LOS     = LSHMAT + L3
         LSCR2   = LOS + NDSR
         LAX     = LSCR2 + L1*L1
         LAY     = LAX + L2
         LAZ     = LAY + L2
         LAST    = LAZ + L2
         NEED1   = LAST - LOADFM - 1
         CALL GETFM(NEED1)
C
         CALL VCLR(X(LINXOV), 1,L7*2)
         CALL VCLR(X(LVROALL),1,L7*NDSR)
         CALL VCLR(X(LVLOALL),1,L7*NDSR)
         CALL VCLR(X(LTDEN),  1,L1*L1)
         CALL VCLR(X(LDEN),   1,L2)
         CALL VCLR(X(LV),     1,L1*L1)
         CALL VCLR(X(LSCR),   1,L1*LX)
         CALL VCLR(X(LOVLP),  1,L2)
         CALL VCLR(X(LMULPOP),1,L2)
         CALL VCLR(X(LMULCH), 1,NAT)
         CALL VCLR(X(LWRK1),  1,L2)
         CALL VCLR(X(LWRK2),  1,L1)
         CALL VCLR(X(LSHMAT), 1,L3)
         CALL VCLR(X(LOS),    1,NDSR)
         CALL VCLR(X(LSCR2),  1,L1*L1)
         CALL VCLR(X(LAX),    1,L2)
         CALL VCLR(X(LAY),    1,L2)
         CALL VCLR(X(LAZ),    1,L2)
         CALL VCLR(X(LX12D),  1,NDSR*NDSR)
         CALL VCLR(X(LX12A),  1,NDSR*NDSR)
         CALL VCLR(X(LEE),    1,NDSR)
         CALL VCLR(X(LQ12D),  1,NDSR*NDSR)
         CALL VCLR(X(LQ12A),  1,NDSR*NDSR)
         CALL VCLR(X(LQ02D),  1,NDSR)
         CALL VCLR(X(LQ02A),  1,NDSR)
         CALL VCLR(X(LTXYZ),  1,3*NDSR)
         CALL VCLR(X(LTXYZ2), 1,3*NDSR*NDSR)
C
         IF(EXETYP.NE.CHECK) THEN
C
C           ----- READ IN GROUND-STATE DENSITY,MO COEFFICIENTS -----
C           ----- TDDFT-RESPONSE VECTORS, TDDFT-EXCITATION ENERGIES -----
C
            CALL DAREAD(IDAF,IODA,X(LDEN),L2,16,0)
            CALL DAREAD(IDAF,IODA,X(LV),L1*LX,15,0)
            CALL DAREAD(IDAF,IODA,X(LVROALL),L7*NDSR,471,0)
            CALL DAREAD(IDAF,IODA,X(LVLOALL),L7*NDSR,479,0)
            CALL DAREAD(IDAF,IODA,X(LEE),NDSR,480,0)
C
            CALL AOLIM
C
C           ----- READ IN OVERLAP MATRIX -----
C           ----- IF NEEDED, TRANSFORM TO SPHERICAL BASIS -----
C
            IF (ISPHER.GE.0.AND..NOT.DFTBFL) THEN
                IFINV = 1
                IFTSP = 0
                CALL SPHMAT(X(LSHMAT),L1,IFINV,IFTSP)
                CALL DAREAD(IDAF,IODA,X(LWRK1),L2,12,0)
                CALL TFTRI(X(LOVLP),X(LWRK1),X(LSHMAT),
     *                     X(LWRK2),L1,L1,L1)
                IFINV = -1
                IFTSP = 1
                CALL SPHMAT(X(LSHMAT),L1,IFINV,IFTSP)
            ELSE
                CALL DAREAD(IDAF,IODA,X(LOVLP),L2,12,0)
            END IF
C
C           --- GROUND STATE DENSITY
C
            CALL DCOPY(L2,X(LDEN),1,X(LSCR),1)
C
C           ----- TRANSFORM TO SPHERICAL HARMONICS IF NEEDED -----
C
            IF (ISPHER.GE.0.AND..NOT.DFTBFL) THEN
               CALL TFTRI(X(LWRK1),X(LSCR),X(LSHMAT),
     *                    X(LWRK2),L1,L1,L1)
            ELSE
               CALL DCOPY(L2,X(LSCR),1,X(LWRK1),1)
            END IF
C
C           ------ CALCULATE MULLIKEN CHARGES (FCD) -----
C
            CALL OVLPOP(X(LWRK1),X(LOVLP),L2)
            CALL ATPOP(X(LWRK1),IA,X(LMULPOP),NAT)
            CALL GROSSC(X(LMULPOP),X(LMULCH),IA,NAT)
C
C           ----- SUM UP RESPECTIVE MULLIKEN CHARGES (FCD) -----
C           ----- THIS ALSO INCLUDES NUCLEAR CHARGES -----
C
            IF(NFG.EQ.0) THEN
               Q0D = ZERO
               DO I=IFEDAT(1),IFEDAT(2)
                  Q0D = Q0D + ZAN(I) - X(LMULCH+I-1)
               END DO
               Q0A = ZERO
               DO I=IFEDAT(3),IFEDAT(4)
                  Q0A = Q0A + ZAN(I) - X(LMULCH+I-1)
               END DO
            ELSE
               CALL FEDFMOSUM(0,X(LMULCH),Q0D,Q0A)
               CALL FEDFMOSUM(1,ZAN,Q0D,Q0A)
            END IF
C
C           --- UNRELAXED EXCITED STATE DENSITY
C
            DO IST=1,NDSR
C
C              ----- CREATE UNRELAXED EXCITED STATE DENSITY (FCD) -----
C
               CALL TCONST2(X(LTDEN),X(LVROALL),X(LVLOALL),
     *                     LX,L7,NOCC,ITF,IST,IST)
               CALL DGEMM('N','N',L1,LX,LX,ONE,X(LV),L1,X(LTDEN),LX,
     *                    ZERO,X(LSCR),L1)
               CALL DGEMM('N','T',L1,L1,LX,ONE,X(LSCR),L1,X(LV),L1,
     *                    ZERO,X(LTDEN),L1)
               CALL SQ2TRI(L1,L1,X(LTDEN),X(LSCR),HALF)
               CALL DAXPY(L2,ONE,X(LDEN),1,X(LSCR),1)
C
C              ----- TRANSFORM TO SPHERICAL HARMONICS IF NEEDED -----
C
               IF (ISPHER.GE.0.AND..NOT.DFTBFL) THEN
                   CALL TFTRI(X(LWRK1),X(LSCR),X(LSHMAT),
     *                        X(LWRK2),L1,L1,L1)
               ELSE
                  CALL DCOPY(L2,X(LSCR),1,X(LWRK1),1)
               END IF
C
C              ------ CALCULATE MULLIKEN CHARGES (FCD) -----
C
               CALL OVLPOP(X(LWRK1),X(LOVLP),L2)
               CALL ATPOP(X(LWRK1),IA,X(LMULPOP),NAT)
               CALL GROSSC(X(LMULPOP),X(LMULCH),IA,NAT)
C
C              ----- SUM UP RESPECTIVE MULLIKEN CHARGES (FCD) -----
C              ----- THIS ALSO INCLUDES NUCLEAR CHARGES -----
C
               IF(NFG.EQ.0) THEN
                  Q12D = ZERO
                  DO I=IFEDAT(1),IFEDAT(2)
                     Q12D = Q12D + ZAN(I) - X(LMULCH+I-1)
                  END DO
                  Q12A = ZERO
                  DO I=IFEDAT(3),IFEDAT(4)
                     Q12A = Q12A + ZAN(I) - X(LMULCH+I-1)
                  END DO
               ELSE
                  CALL FEDFMOSUM(0,X(LMULCH),Q12D,Q12A)
                  CALL FEDFMOSUM(1,ZAN,Q12D,Q12A)
               END IF
               X(LQ12D + (IST-1)*NDSR + (IST-1)) = Q12D
               X(LQ12A + (IST-1)*NDSR + (IST-1)) = Q12A
C
C              ----- CREATE UNRELAXED EXCITED STATE DENSITY (FED) -----
C
               CALL TCONST2(X(LTDEN),X(LVROALL),X(LVLOALL),
     *                     LX,L7,NOCC,ITF,IST,IST)
               DO I=1,NOCC
                  DO J=1,NOCC
                     ITMP = LTDEN+(I-1)*LX+(J-1)
                     X(ITMP) = -X(ITMP)
                  END DO
               END DO
               CALL DGEMM('N','N',L1,LX,LX,ONE,X(LV),L1,X(LTDEN),LX,
     *                    ZERO,X(LSCR),L1)
               CALL DGEMM('N','T',L1,L1,LX,ONE,X(LSCR),L1,X(LV),L1,
     *                    ZERO,X(LTDEN),L1)
               CALL SQ2TRI(L1,L1,X(LTDEN),X(LSCR),HALF)
C
C              ----- TRANSFORM TO SPHERICAL HARMONICS IF NEEDED -----
C
               IF (ISPHER.GE.0.AND..NOT.DFTBFL) THEN
                  CALL TFTRI(X(LWRK1),X(LSCR),X(LSHMAT),
     *                       X(LWRK2),L1,L1,L1)
               ELSE
                  CALL DCOPY(L2,X(LSCR),1,X(LWRK1),1)
               END IF
C
C              ------ CALCULATE MULLIKEN CHARGES (FED) -----
C
               CALL OVLPOP(X(LWRK1),X(LOVLP),L2)
               CALL ATPOP(X(LWRK1),IA,X(LMULPOP),NAT)
               CALL GROSSC(X(LMULPOP),X(LMULCH),IA,NAT)
C
C              ----- SUM UP RESPECTIVE MULLIKEN CHARGES (FED) -----
C
               IF(NFG.EQ.0) THEN
                  X12D = ZERO
                  DO I=IFEDAT(1),IFEDAT(2)
                     X12D = X12D - X(LMULCH+I-1)
                  END DO
                  X12A = ZERO
                  DO I=IFEDAT(3),IFEDAT(4)
                     X12A = X12A - X(LMULCH+I-1)
                  END DO
               ELSE
                  CALL FEDFMOSUM(0,X(LMULCH),X12D,X12A)
               END IF
               X(LX12D + (IST-1)*NDSR + (IST-1)) = X12D
               X(LX12A + (IST-1)*NDSR + (IST-1)) = X12A
            END DO
C
C           --- GROUND TO EXCITED STATE TRANSITION DENSITY
C
            CALL EXCnstLab(X(LINXOV),NOCC,LX)
            DO IST=1,NDSR
C
C              ----- CREATE TRANSITION DENSITY (FCD) -----
C
               LVRO = LVROALL+(IST-1)*L7
               LVLO = LVLOALL+(IST-1)*L7
               CALL TRAD(X(LTDEN),X(LVRO),X(LVLO),LX,L7,
     *                   X(LINXOV),1,1)
C              possibly better (check!):
C              CALL TRAD(X(LTDEN),X(LVROALL),X(LVLOALL),LX,L7,
C    *                   X(LINXOV),NDSR,IST)
               CALL DSCAL(L1*L1,SQRT_TWO_HALF,X(LTDEN),1)
               CALL DGEMM('N','N',L1,LX,LX,ONE,X(LV),L1,X(LTDEN),LX,
     *                    ZERO,X(LSCR),L1)
               CALL DGEMM('N','T',L1,L1,LX,ONE,X(LSCR),L1,X(LV),L1,
     *                    ZERO,X(LTDEN),L1)
               CALL SQ2TRI(L1,L1,X(LTDEN),X(LSCR),HALF)
C
C              ----- TRANSFORM TO SPHERICAL HARMONICS IF NEEDED -----
C
               IF (ISPHER.GE.0.AND..NOT.DFTBFL) THEN
                   CALL TFTRI(X(LWRK1),X(LSCR),X(LSHMAT),
     *                        X(LWRK2),L1,L1,L1)
               ELSE
                  CALL DCOPY(L2,X(LSCR),1,X(LWRK1),1)
               END IF
C
C              ------ CALCULATE MULLIKEN CHARGES (FCD) -----
C
               CALL OVLPOP(X(LWRK1),X(LOVLP),L2)
               CALL ATPOP(X(LWRK1),IA,X(LMULPOP),NAT)
               CALL GROSSC(X(LMULPOP),X(LMULCH),IA,NAT)
C
C              ----- SUM UP RESPECTIVE MULLIKEN CHARGES (FCD)  -----
C
               IF(NFG.EQ.0) THEN
                  Q02D = ZERO
                  DO I=IFEDAT(1),IFEDAT(2)
                     Q02D = Q02D - X(LMULCH+I-1)
                  END DO
                  Q02A = ZERO
                  DO I=IFEDAT(3),IFEDAT(4)
                      Q02A = Q02A - X(LMULCH+I-1)
                  END DO
               ELSE
                  CALL FEDFMOSUM(0,X(LMULCH),Q02D,Q02A)
               END IF
               X(LQ02D + (IST-1)) = Q02D
               X(LQ02A + (IST-1)) = Q02A
            END DO
C
C           --- CALCULATE TRANSITION DIPOLE MOMENTS ---
C
            CALL TDOSCALC(X(LOS),X(LTXYZ),X(LTXYZ2),X(LTDEN),X(LSCR2),
     *                    X(LVROALL),X(LVLOALL),X(LEE),X(LV),
     *                    X(LAX),X(LAY),X(LAZ),L1,LX,L2,L7,
     *                    NOCC,NDSR,X(LINXOV),ITF)
C
C           --- UNRELAXED EXCITED TO EXCITED STATE TRANSITION DENSITY
C
            DO IST=1,NDSR
               DO JST=IST+1,NDSR
C
C                 ----- CREATE TRANSITION DENSITY (FCD) -----
C
                  CALL TCONST2(X(LTDEN),X(LVROALL),X(LVLOALL),
     *                         LX,L7,NOCC,ITF,IST,JST)
                  CALL DGEMM('N','N',L1,LX,LX,ONE,X(LV),L1,X(LTDEN),LX,
     *                       ZERO,X(LSCR),L1)
                  CALL DGEMM('N','T',L1,L1,LX,ONE,X(LSCR),L1,X(LV),L1,
     *                       ZERO,X(LTDEN),L1)
                  CALL SQ2TRI(L1,L1,X(LTDEN),X(LSCR),HALF)
C
C                 ----- TRANSFORM TO SPHERICAL HARMONICS IF NEEDED -----
C
                  IF (ISPHER.GE.0.AND..NOT.DFTBFL) THEN
                     CALL TFTRI(X(LWRK1),X(LSCR),X(LSHMAT),
     *                          X(LWRK2),L1,L1,L1)
                  ELSE
                     CALL DCOPY(L2,X(LSCR),1,X(LWRK1),1)
                  END IF
C
C                 ------ CALCULATE MULLIKEN CHARGES (FCD) -----
C
                  CALL OVLPOP(X(LWRK1),X(LOVLP),L2)
                  CALL ATPOP(X(LWRK1),IA,X(LMULPOP),NAT)
                  CALL GROSSC(X(LMULPOP),X(LMULCH),IA,NAT)
C
C                 ----- SUM UP RESPECTIVE MULLIKEN CHARGES (FCD)  -----
C
                  IF(NFG.EQ.0) THEN
                     Q12D = ZERO
                     DO I=IFEDAT(1),IFEDAT(2)
                        Q12D = Q12D - X(LMULCH+I-1)
                     END DO
                     Q12A = ZERO
                     DO I=IFEDAT(3),IFEDAT(4)
                         Q12A = Q12A - X(LMULCH+I-1)
                     END DO
                  ELSE
                     CALL FEDFMOSUM(0,X(LMULCH),Q12D,Q12A)
                  END IF
                  X(LQ12D + (IST-1)*NDSR + (JST-1)) = Q12D
                  X(LQ12D + (JST-1)*NDSR + (IST-1)) = Q12D
                  X(LQ12A + (IST-1)*NDSR + (JST-1)) = Q12A
                  X(LQ12A + (JST-1)*NDSR + (IST-1)) = Q12A
C
C                 ----- CREATE TRANSITION DENSITY (FED) -----
C
                  CALL TCONST2(X(LTDEN),X(LVROALL),X(LVLOALL),
     *                        LX,L7,NOCC,ITF,IST,JST)
                  DO I=1,NOCC
                     DO J=1,NOCC
                        ITMP = LTDEN+(I-1)*LX+(J-1)
                        X(ITMP) = -X(ITMP)
                     END DO
                  END DO
                  CALL DGEMM('N','N',L1,LX,LX,ONE,X(LV),L1,X(LTDEN),LX,
     *                       ZERO,X(LSCR),L1)
                  CALL DGEMM('N','T',L1,L1,LX,ONE,X(LSCR),L1,X(LV),L1,
     *                       ZERO,X(LTDEN),L1)
                  CALL SQ2TRI(L1,L1,X(LTDEN),X(LSCR),HALF)
C
C                 ----- TRANSFORM TO SPHERICAL HARMONICS IF NEEDED -----
C
                  IF (ISPHER.GE.0.AND..NOT.DFTBFL) THEN
                     CALL TFTRI(X(LWRK1),X(LSCR),X(LSHMAT),
     *                          X(LWRK2),L1,L1,L1)
                  ELSE
                     CALL DCOPY(L2,X(LSCR),1,X(LWRK1),1)
                  END IF
C
C                 ------ CALCULATE MULLIKEN CHARGES (FED) -----
C
                  CALL OVLPOP(X(LWRK1),X(LOVLP),L2)
                  CALL ATPOP(X(LWRK1),IA,X(LMULPOP),NAT)
                  CALL GROSSC(X(LMULPOP),X(LMULCH),IA,NAT)
C
C                 ----- SUM UP RESPECTIVE MULLIKEN CHARGES (FED) -----
C
                  IF(NFG.EQ.0) THEN
                     X12D = ZERO
                     DO I=IFEDAT(1),IFEDAT(2)
                        X12D = X12D - X(LMULCH+I-1)
                     END DO
                     X12A = ZERO
                     DO I=IFEDAT(3),IFEDAT(4)
                        X12A = X12A - X(LMULCH+I-1)
                     END DO
                  ELSE
                     CALL FEDFMOSUM(0,X(LMULCH),X12D,X12A)
                  END IF
                  X(LX12D + (IST-1)*NDSR + (JST-1)) = X12D
                  X(LX12D + (JST-1)*NDSR + (IST-1)) = X12D
                  X(LX12A + (IST-1)*NDSR + (JST-1)) = X12A
                  X(LX12A + (JST-1)*NDSR + (IST-1)) = X12A
               END DO
            END DO
         END IF
C
C        --- FREE WORKING MEMORY FOR FED/FCD ---
C
         CALL RETFM(NEED1)
C
C        ----- MULTI-FED-FCD ---
C
C        The general idea is to just take all NSTATE
C        states for the diabatic Hamiltonian.
C
C        See: https://doi.org/10.1007/s11120-018-0492-1
C
C        --- ALLOCATE MEMORY FOR MULTI-FED-FCD RESULTS ---
C
         CALL VALFM(LOADFM)
         LHFIN    = LOADFM + 1
         LTXYZMF  = LHFIN + NDSR*NDSR
         LTXYZ2MF = LTXYZMF + 3*NDSR
         LBCT1D   = LTXYZ2MF + 3*NDSR*NDSR
         LX12DL   = LBCT1D + NDSR*NDSR
         LX12AL   = LX12DL + NDSR
         LQ12DL   = LX12AL + NDSR
         LQ12AL   = LQ12DL + NDSR
         LAST     = LQ12AL + NDSR
         NEED1    = LAST - LOADFM - 1
         CALL GETFM(NEED1)
C
C        --- ALLOCATE WORKING MEMORY FOR MULTI-FED-FCD ---
C
         CALL VALFM(LOADFM)
         LXMAT   = LOADFM + 1
         LQMAT   = LXMAT + NDSR*NDSR
         LDMAT   = LQMAT + NDSR*NDSR
         LW      = LDMAT + NDSR*NDSR
         LWRK1   = LW + NDSR*NDSR
         LWRK2   = LWRK1 + 3*NDSR-1
         LWRK3   = LWRK2 + NDSR*NDSR
         LBLE1   = LWRK3 + NDSR*NDSR
         LBLE2   = LBLE1 + NDSR*NDSR
         LBCT1   = LBLE2 + NDSR*NDSR
         LBCT2   = LBCT1 + NDSR*NDSR
         LHLE1   = LBCT2 + NDSR*NDSR
         LHLE2   = LHLE1 + NDSR*NDSR
         LHCT1   = LHLE2 + NDSR*NDSR
         LHCT2   = LHCT1 + NDSR*NDSR
         LBLE1D  = LHCT2 + NDSR*NDSR
         LBLE2D  = LBLE1D + NDSR*NDSR
         LBCT2D  = LBLE2D + NDSR*NDSR
         LDONI   = LBCT2D + NDSR*NDSR
         LACCI   = LDONI + NDSR
         LAST    = LACCI + NDSR
         NEED2   = LAST - LOADFM - 1
         CALL GETFM(NEED2)
C
         CALL VCLR(X(LHFIN),   1, NDSR*NDSR)
         CALL VCLR(X(LTXYZMF), 1, 3*NDSR)
         CALL VCLR(X(LTXYZ2MF),1, 3*NDSR*NDSR)
         CALL VCLR(X(LBCT1D),  1, NDSR*NDSR)
         CALL VCLR(X(LX12DL),  1, NDSR)
         CALL VCLR(X(LX12DL),  1, NDSR)
         CALL VCLR(X(LX12DL),  1, NDSR)
         CALL VCLR(X(LX12DL),  1, NDSR)
         CALL VCLR(X(LXMAT),   1, NDSR*NDSR)
         CALL VCLR(X(LQMAT),   1, NDSR*NDSR)
         CALL VCLR(X(LDMAT),   1, NDSR*NDSR)
         CALL VCLR(X(LW),      1, NDSR*NDSR)
         CALL VCLR(X(LWRK1),   1, 3*NDSR-1)
         CALL VCLR(X(LWRK2),   1, NDSR*NDSR)
         CALL VCLR(X(LWRK3),   1, NDSR*NDSR)
         CALL VCLR(X(LBLE1),   1, NDSR*NDSR)
         CALL VCLR(X(LBLE2),   1, NDSR*NDSR)
         CALL VCLR(X(LBCT1),   1, NDSR*NDSR)
         CALL VCLR(X(LBCT2),   1, NDSR*NDSR)
         CALL VCLR(X(LHLE1),   1, NDSR*NDSR)
         CALL VCLR(X(LHLE2),   1, NDSR*NDSR)
         CALL VCLR(X(LHCT1),   1, NDSR*NDSR)
         CALL VCLR(X(LHCT2),   1, NDSR*NDSR)
         CALL VCLR(X(LBLE1D),  1, NDSR*NDSR)
         CALL VCLR(X(LBLE2D),  1, NDSR*NDSR)
         CALL VCLR(X(LBCT2D),  1, NDSR*NDSR)
         CALL VICLR(X(LDONI),  1, NDSR)
         CALL VICLR(X(LACCI),  1, NDSR)
C
         IF(EXETYP.NE.CHECK) THEN
C
C           1)  We can easily generate the D matrix.
C               It is written as: D = Q Q - X X
C               This is \delta x_mn (eq. 4).
C
            CALL DCOPY(NDSR*NDSR,X(LX12D),1,X(LXMAT),1)
            CALL DAXPY(NDSR*NDSR,-ONE,X(LX12A),1,X(LXMAT),1)
            CALL DSCAL(NDSR*NDSR,HALF,X(LXMAT),1)
            CALL DCOPY(NDSR*NDSR,X(LQ12D),1,X(LQMAT),1)
            CALL DAXPY(NDSR*NDSR,-ONE,X(LQ12A),1,X(LQMAT),1)
            CALL DSCAL(NDSR*NDSR,HALF,X(LQMAT),1)
            CALL DGEMM('N','N',NDSR,NDSR,NDSR,
     *                 ONE,
     *                 X(LQMAT),NDSR,
     *                 X(LQMAT),NDSR,
     *                 ONE,
     *                 X(LDMAT),NDSR)
            CALL DGEMM('N','N',NDSR,NDSR,NDSR,
     *                 -ONE,
     *                 X(LXMAT),NDSR,
     *                 X(LXMAT),NDSR,
     *                 ONE,
     *                 X(LDMAT),NDSR)
C
C           2)  The eigenvectors of the D matrix
C               belonging to eigenvalue +1 are the LE
C               states while the eigenvectors belonging
C               to eigenvalue -1 are the CT states (eq. 7).
C
            CALL DSYEV('V','U',NDSR,X(LDMAT),NDSR,X(LW),
     *                 X(LWRK1),3*NDSR-1,INFO)
            DO ILECT=1,NDSR
               IF(X(LW+ILECT-1).GT.ZERO) EXIT
            END DO
C
C           3a) We calculate the projection Xdash of X into
C               the subspace of LE states. Then we diagonalize
C               Xdash (eq. 8). The eigenvectors of Xdash belonging
C               to eigenvalue +1 are LE1 states while the
C               eigenvectors of Xdash belonging to eigenvalue -1
C               are LE2 states. We have to transform back in order
C               to get eigenvectors of the D matrix.
C
            IF(ILECT.GT.1) THEN
               CALL DGEMM('N','N',NDSR,ILECT-1,NDSR,
     *                    ONE,
     *                    X(LXMAT),NDSR,
     *                    X(LDMAT),NDSR,
     *                    ZERO,
     *                    X(LWRK2),NDSR)
               CALL DGEMM('T','N',ILECT-1,ILECT-1,NDSR,
     *                    ONE,
     *                    X(LDMAT),NDSR,
     *                    X(LWRK2),NDSR,
     *                    ZERO,
     *                    X(LXMAT),ILECT-1)
               CALL DSYEV('V','U',ILECT-1,X(LXMAT),ILECT-1,X(LW),
     *                    X(LWRK1),3*NDSR-1,INFO)
               DO ILE2LE1=1,ILECT-1
                  IF(X(LW+ILE2LE1-1).GT.ZERO) EXIT
               END DO
               CALL DGEMM('N','N',NDSR,ILE2LE1-1,ILECT-1,
     *                    ONE,
     *                    X(LDMAT),NDSR,
     *                    X(LXMAT),ILECT-1,
     *                    ZERO,X(LBLE2),NDSR)
               CALL DGEMM('N','N',NDSR,ILECT-ILE2LE1,ILECT-1,
     *                    ONE,
     *                    X(LDMAT),NDSR,
     *                    X(LXMAT+(ILE2LE1-1)*(ILECT-1)),ILECT-1,
     *                    ZERO,X(LBLE1),NDSR)
            ELSE
               ILE2LE1=1
            END IF
C
C           3b) We calculate the projection Qdash of Q into
C               the subspace of CT states. Then we diagonalize
C               Qdash (eq. 8). The eigenvectors of Qdash belonging
C               to eigenvalue +1 are CT1 states while the
C               eigenvectors of Qdash belonging to eigenvalue -1
C               are CT2 states. We transform back in order to
C               get eigenvectors of the D matrix.
C
            IF(ILECT.LE.NDSR) THEN
               CALL DGEMM('N','N',NDSR,NDSR-ILECT+1,NDSR,
     *                    ONE,
     *                    X(LQMAT),NDSR,
     *                    X(LDMAT+(ILECT-1)*NDSR),NDSR,
     *                    ZERO,
     *                    X(LWRK2),NDSR)
               CALL DGEMM('T','N',NDSR-ILECT+1,NDSR-ILECT+1,NDSR,
     *                    ONE,
     *                    X(LDMAT+(ILECT-1)*NDSR),NDSR,
     *                    X(LWRK2),NDSR,
     *                    ZERO,
     *                    X(LQMAT),NDSR-ILECT+1)
               CALL DSYEV('V','U',NDSR-ILECT+1,
     *                    X(LQMAT),NDSR-ILECT+1,X(LW),
     *                    X(LWRK1),3*NDSR-1,INFO)
               DO ICT2CT1=1,NDSR-ILECT+1
                  IF(X(LW+ICT2CT1-1).GT.ZERO) EXIT
               END DO
               CALL DGEMM('N','N',NDSR,ICT2CT1-1,NDSR-ILECT+1,
     *                    ONE,
     *                    X(LDMAT+(ILECT-1)*NDSR),NDSR,
     *                    X(LQMAT),NDSR-ILECT+1,
     *                    ZERO,
     *                    X(LBCT2),NDSR)
               CALL DGEMM('N','N',
     *                    NDSR,NDSR-ILECT-ICT2CT1+2,NDSR-ILECT+1,
     *                    ONE,
     *                    X(LDMAT+(ILECT-1)*NDSR),NDSR,
     *                    X(LQMAT+(ICT2CT1-1)*(NDSR-ILECT+1)),
     *                    NDSR-ILECT+1,
     *                    ZERO,
     *                    X(LBCT1),NDSR)
            ELSE
               ICT2CT1 = 1
            END IF
C
C           4)  Now each of our localized states has one of the
C               four labels CT1,CT2,LE1,LE2 but the Hamiltonian
C               in each of the subspaces HCT1, HCT2, HLE1, HLE2
C               may not be diagonal. So we need to diagonalize
C               in the subspaces and then transform back.
C
            IF(ILECT.GT.1.AND.ILE2LE1.GT.1) THEN
               CALL VCLR(X(LWRK2),1,NDSR*NDSR)
               DO I=1,NDSR
                  X(LWRK2+(I-1)*NDSR+I-1)=X(LEE+I-1)
               END DO
               CALL DGEMM('N','N',NDSR,ILE2LE1-1,NDSR,
     *                    ONE,
     *                    X(LWRK2),NDSR,
     *                    X(LBLE2),NDSR,
     *                    ZERO,
     *                    X(LWRK3),NDSR)
               CALL DGEMM('T','N',ILE2LE1-1,ILE2LE1-1,NDSR,
     *                    ONE,
     *                    X(LBLE2),NDSR,
     *                    X(LWRK3),NDSR,
     *                    ZERO,
     *                    X(LHLE2),ILE2LE1-1)
               CALL DSYEV('V','U',ILE2LE1-1,
     *                    X(LHLE2),ILE2LE1-1,
     *                    X(LW),X(LWRK1),3*NDSR-1,INFO)
               CALL DGEMM('N','N',NDSR,ILE2LE1-1,ILE2LE1-1,
     *                    ONE,
     *                    X(LBLE2),NDSR,
     *                    X(LHLE2),ILE2LE1-1,
     *                    ZERO,
     *                    X(LBLE2D),NDSR)
            END IF
            IF(ILECT.GT.1.AND.ILE2LE1.LT.ILECT) THEN
               CALL VCLR(X(LWRK2),1,NDSR*NDSR)
               DO I=1,NDSR
                  X(LWRK2+(I-1)*NDSR+I-1)=X(LEE+I-1)
               END DO
               CALL DGEMM('N','N',NDSR,ILECT-ILE2LE1,NDSR,
     *                    ONE,
     *                    X(LWRK2),NDSR,
     *                    X(LBLE1),NDSR,
     *                    ZERO,
     *                    X(LWRK3),NDSR)
               CALL DGEMM('T','N',ILECT-ILE2LE1,ILECT-ILE2LE1,NDSR,
     *                    ONE,
     *                    X(LBLE1),NDSR,
     *                    X(LWRK3),NDSR,
     *                    ZERO,
     *                    X(LHLE1),ILECT-ILE2LE1)
               CALL DSYEV('V','U',ILECT-ILE2LE1,
     *                    X(LHLE1),ILECT-ILE2LE1,
     *                    X(LW),X(LWRK1),3*NDSR-1,INFO)
               CALL DGEMM('N','N',NDSR,ILECT-ILE2LE1,ILECT-ILE2LE1,
     *                    ONE,
     *                    X(LBLE1),NDSR,
     *                    X(LHLE1),ILECT-ILE2LE1,
     *                    ZERO,
     *                    X(LBLE1D),NDSR)
            END IF
            IF(ILECT.LE.NDSR.AND.ICT2CT1.GT.1) THEN
               CALL VCLR(X(LWRK2),1,NDSR*NDSR)
               DO I=1,NDSR
                  X(LWRK2+(I-1)*NDSR+I-1)=X(LEE+I-1)
               END DO
               CALL DGEMM('N','N',NDSR,ICT2CT1-1,NDSR,
     *                    ONE,
     *                    X(LWRK2),NDSR,
     *                    X(LBCT2),NDSR,
     *                    ZERO,
     *                    X(LWRK3),NDSR)
               CALL DGEMM('T','N',ICT2CT1-1,ICT2CT1-1,NDSR,
     *                    ONE,
     *                    X(LBCT2),NDSR,
     *                    X(LWRK3),NDSR,
     *                    ZERO,
     *                    X(LHCT2),ICT2CT1-1)
               CALL DSYEV('V','U',ICT2CT1-1,
     *                    X(LHCT2),ICT2CT1-1,
     *                    X(LW),X(LWRK1),3*NDSR-1,INFO)
               CALL DGEMM('N','N',NDSR,ICT2CT1-1,ICT2CT1-1,
     *                    ONE,
     *                    X(LBCT2),NDSR,
     *                    X(LHCT2),ICT2CT1-1,
     *                    ZERO,
     *                    X(LBCT2D),NDSR)
            END IF
            IF(ILECT.LE.NDSR.AND.ICT2CT1.LE.NDSR-ILECT+1) THEN
               CALL VCLR(X(LWRK2),1,NDSR*NDSR)
               DO I=1,NDSR
                  X(LWRK2+(I-1)*NDSR+I-1)=X(LEE+I-1)
               END DO
               CALL DGEMM('N','N',NDSR,NDSR-ILECT-ICT2CT1+2,NDSR,
     *                    ONE,
     *                    X(LWRK2),NDSR,
     *                    X(LBCT1),NDSR,
     *                    ZERO,
     *                    X(LWRK3),NDSR)
               CALL DGEMM('T','N',
     *                    NDSR-ILECT-ICT2CT1+2,
     *                    NDSR-ILECT-ICT2CT1+2,NDSR,
     *                    ONE,
     *                    X(LBCT1),NDSR,
     *                    X(LWRK3),NDSR,
     *                    ZERO,
     *                    X(LHCT1),NDSR-ILECT-ICT2CT1+2)
               CALL DSYEV('V','U',NDSR-ILECT-ICT2CT1+2,
     *                    X(LHCT1),NDSR-ILECT-ICT2CT1+2,
     *                    X(LW),X(LWRK1),3*NDSR-1,INFO)
               CALL DGEMM('N','N',
     *                    NDSR,NDSR-ILECT-ICT2CT1+2,
     *                    NDSR-ILECT-ICT2CT1+2,
     *                    ONE,
     *                    X(LBCT1),NDSR,
     *                    X(LHCT1),NDSR-ILECT-ICT2CT1+2,
     *                    ZERO,
     *                    X(LBCT1D),NDSR)
            END IF
C
C           5)  We assemble the basis vectors of the subspaces
C               into one big basis. With this we can finally
C               determine H''' (eq. 9).
C
            IF(ILECT.LE.NDSR.AND.ICT2CT1.GT.1) THEN
                CALL DCOPY((ICT2CT1-1)*NDSR,
     *                     X(LBCT2D),1,
     *                     X(LBCT1D+(NDSR-ILECT-ICT2CT1+2)*NDSR),1)
            END IF
            IF(ILECT.GT.1.AND.ILE2LE1.LT.ILECT) THEN
                CALL DCOPY((ILECT-ILE2LE1)*NDSR,
     *                     X(LBLE1D),1,
     *                     X(LBCT1D+(NDSR-ILECT+1)*NDSR),1)
            END IF
            IF(ILECT.GT.1.AND.ILE2LE1.GT.1) THEN
                CALL DCOPY((ILE2LE1-1)*NDSR,
     *                     X(LBLE2D),1,
     *                     X(LBCT1D+(NDSR-ILE2LE1+1)*NDSR),1)
            END IF
            CALL VCLR(X(LWRK2),1,NDSR*NDSR)
            DO I=1,NDSR
               X(LWRK2+(I-1)*NDSR+I-1)=X(LEE+I-1)
            END DO
            CALL DGEMM('N','N',NDSR,NDSR,NDSR,
     *                 ONE,
     *                 X(LWRK2),NDSR,
     *                 X(LBCT1D),NDSR,
     *                 ZERO,
     *                 X(LWRK3),NDSR)
            CALL DGEMM('T','N',NDSR,NDSR,NDSR,
     *                 ONE,
     *                 X(LBCT1D),NDSR,
     *                 X(LWRK3),NDSR,
     *                 ZERO,
     *                 X(LHFIN),NDSR)
            CALL DSCAL(NDSR*NDSR,TOEV,X(LHFIN),1)
            CALL SQ2TRI(NDSR,NDSR,X(LHFIN),X(LWRK2),HALF)
            CALL DCOPY((NDSR*(NDSR+1))/2,X(LWRK2),1,X(LHFIN),1)
C
C           6)  Transform TXYZ and TXYZ2 using the final basis.
C
            CALL DGEMM('N','N',3,NDSR,NDSR,
     *                 ONE,
     *                 X(LTXYZ),3,
     *                 X(LBCT1D),NDSR,
     *                 ZERO,
     *                 X(LTXYZMF),3)
C           TODO: Transform TXYZ2
C
C           7) Transform X and Q matrices using the final basis.
C
            CALL DGEMM('N','N',NDSR,NDSR,NDSR,
     *                 ONE,
     *                 X(LX12D),NDSR,
     *                 X(LBCT1D),NDSR,
     *                 ZERO,
     *                 X(LWRK2),NDSR)
            CALL DGEMM('T','N',NDSR,NDSR,NDSR,
     *                 ONE,
     *                 X(LBCT1D),NDSR,
     *                 X(LWRK2),NDSR,
     *                 ZERO,
     *                 X(LWRK3),NDSR)
            DO I=1,NDSR
               X(LX12DL+I-1)=X(LWRK3+(I-1)*NDSR+I-1)
            END DO
            CALL DGEMM('N','N',NDSR,NDSR,NDSR,
     *                 ONE,
     *                 X(LX12A),NDSR,
     *                 X(LBCT1D),NDSR,
     *                 ZERO,
     *                 X(LWRK2),NDSR)
            CALL DGEMM('T','N',NDSR,NDSR,NDSR,
     *                 ONE,
     *                 X(LBCT1D),NDSR,
     *                 X(LWRK2),NDSR,
     *                 ZERO,
     *                 X(LWRK3),NDSR)
            DO I=1,NDSR
               X(LX12AL+I-1)=X(LWRK3+(I-1)*NDSR+I-1)
            END DO
            CALL DGEMM('N','N',NDSR,NDSR,NDSR,
     *                 ONE,
     *                 X(LQ12D),NDSR,
     *                 X(LBCT1D),NDSR,
     *                 ZERO,
     *                 X(LWRK2),NDSR)
            CALL DGEMM('T','N',NDSR,NDSR,NDSR,
     *                 ONE,
     *                 X(LBCT1D),NDSR,
     *                 X(LWRK2),NDSR,
     *                 ZERO,
     *                 X(LWRK3),NDSR)
            DO I=1,NDSR
               X(LQ12DL+I-1)=X(LWRK3+(I-1)*NDSR+I-1)
            END DO
            CALL DGEMM('N','N',NDSR,NDSR,NDSR,
     *                 ONE,
     *                 X(LQ12A),NDSR,
     *                 X(LBCT1D),NDSR,
     *                 ZERO,
     *                 X(LWRK2),NDSR)
            CALL DGEMM('T','N',NDSR,NDSR,NDSR,
     *                 ONE,
     *                 X(LBCT1D),NDSR,
     *                 X(LWRK2),NDSR,
     *                 ZERO,
     *                 X(LWRK3),NDSR)
            DO I=1,NDSR
               X(LQ12AL+I-1)=X(LWRK3+(I-1)*NDSR+I-1)
            END DO
C
C           8) For FMO: Match localized LE states to
C              states from FMO1-TDDFT. Write energies
C              and coupling.
C
            IF(NFG.NE.0) THEN
               CALL FEDFMOMATCH(
     *                NDSR,X(LDONI),X(LACCI),X(LHFIN),
     *                X(LX12DL),X(LX12AL),
     *                ILECT,ILE2LE1,ICURFGST,JCURFGST,
     *                DOMULTIST,NDON,NACC)
            END IF
         END IF
C
C         --- FREE WORKING MEMORY FOR MULTI-FED-FCD ---
C
         CALL RETFM(NEED2)
         IF(EXETYP.NE.CHECK) THEN
C
C           ----- PRINT OUT RESULTS -----
C
            IF(MASWRK) THEN
               WRITE(IW,"(/10X,'-----------------------------------')")
               WRITE(IW,"(10X,'FED ELECTRONIC-COUPLING CALCULATION')")
               IF(NFG.NE.0) THEN
                  WRITE(IW,"(10X,'IN FMO, "//
     *                     "DONOR=',I5,' ACCEPTOR=',I5)") jcurfg,icurfg
               END IF
               WRITE(IW,"(10X,'-----------------------------------')")
               WRITE(IW,"(/7X,'FRAGMENT EXCITATIONS OF "//
     *                  "SINGLET EXCITED STATE')")
               WRITE(IW,"(/3X,'STATE',9X,'X(D)',8X,'X(A)',10X,'DX')")
            END IF
            DO IST=1,NDSR
               X12D = X(LX12D + (IST-1)*NDSR + (IST-1))
               X12A = X(LX12A + (IST-1)*NDSR + (IST-1))
               DX12 = X12D - X12A
               IF(MASWRK) THEN
                  WRITE(IW,"(I8,F13.6,F12.6,F12.6)") IST,X12D,X12A,DX12
               END IF
            END DO
            IF(MASWRK) THEN
               WRITE(IW,"(/15X,'FRAGMENT EXCITATIONS OF "//
     *                  "TRANSITION DENSITIES AND')")
               WRITE(IW,"(12X,'FED COUPLINGS BETWEEN "//
     *                  "SINGLET EXCITED STATES')")
               WRITE(IW,"(/4X,'STATES',5X,'X12(D)',5X,'X12(A)',7X,"//
     *                  "'DX12',3X,'COUPLING(EV)',4X,"//
     *                  "'FED-DX1',4X,'FED-DX2',12X,'EL1',12X,'EL2')")
            END IF
            DO IST=1,NDSR
               DO JST=1,NDSR
                  IF(IST.GE.JST) CYCLE
                  X12D = X(LX12D + (IST-1)*NDSR + (JST-1))
                  X12A = X(LX12A + (IST-1)*NDSR + (JST-1))
                  DX12 = X12D - X12A
                  DX1  =       X(LX12D + (IST-1)*NDSR + (IST-1))
                  DX1  = DX1 - X(LX12A + (IST-1)*NDSR + (IST-1))
                  DX2  =       X(LX12D + (JST-1)*NDSR + (JST-1))
                  DX2  = DX2 - X(LX12A + (JST-1)*NDSR + (JST-1))
                  E1   = X(LEE + IST-1)
                  E2   = X(LEE + JST-1)
                  IF(ABS(DX1-DX2).GT.TINY(DX1)) THEN
                     COUPL = (E1-E2)*ABS(DX12)/
     *                       SQRT((DX1-DX2)**2+4*DX12**2)
                     TAN2TH  = 2*DX12 / (DX2-DX1)
                     TANTH   = TAN2TH / (ONE + SQRT(1+TAN2TH**2))
                     COSSQTH = ONE / (ONE + TANTH**2)
                     SINSQTH = TANTH**2 * COSSQTH
                     EL1 = E1*COSSQTH + E2*SINSQTH
                     EL2 = E1*SINSQTH + E2*COSSQTH
C                    SINCOSTH = TANTH * COSSQTH
C                    COUPL = (E1-E2)*SINCOSTH
                  ELSE
                     COUPL = ZERO
                     EL1 = E1
                     EL2 = E2
                  ENDIF
                  IF(MASWRK) THEN
                     WRITE(IW,"(I5,I5,"//
     *                        "F11.6,F11.6,F11.6,"//
     *                        "F15.7,"//
     *                        "F11.6,F11.6,"//
     *                        "F15.7,F15.7)")
     *                        IST,JST,X12D,X12A,DX12,COUPL*TOEV,DX1,DX2,
     *                        EL1*TOEV,EL2*TOEV
                  END IF
               END DO
            END DO
            IF(MASWRK) THEN
               WRITE(IW,"(/1X,'...... END OF FED CALCULATION ......'/)")
            END IF
            IST=0
            DQ0 = Q0D - Q0A
            IF(MASWRK) THEN
               WRITE(IW,"(10X,'-----------------------------------')")
               WRITE(IW,"(10X,'FCD ELECTRONIC-COUPLING CALCULATION')")
               IF(NFG.NE.0) THEN
                  WRITE(IW,"(10X,'IN FMO, "//
     *                     "DONOR=',I5,' ACCEPTOR=',I5)") jcurfg,icurfg
               END IF
               WRITE(IW,"(10X,'-----------------------------------')")
               WRITE(IW,"(/7X,'FRAGMENT CHARGES OF "//
     *                  "GROUND STATE "//
     *                  "WITH NUCLEAR CHARGES')")
               WRITE(IW,"(/3X,'STATE',9X,'Q(D)',8X,'Q(A)',10X,'DQ')")
               WRITE(IW,"(I8,F13.6,F12.6,F12.6)") IST,Q0D,Q0A,DQ0
            END IF
            IF(MASWRK) THEN
               WRITE(IW,"(/7X,'FRAGMENT CHARGES OF "//
     *                  "SINGLET EXCITED STATE "//
     *                  "WITH NUCLEAR CHARGES')")
               WRITE(IW,"(/3X,'STATE',9X,'Q(D)',8X,'Q(A)',10X,'DQ')")
            END IF
            DO IST=1,NDSR
               Q12D = X(LQ12D + (IST-1)*NDSR + (IST-1))
               Q12A = X(LQ12A + (IST-1)*NDSR + (IST-1))
               DQ12 = Q12D - Q12A
               IF(MASWRK) THEN
                  WRITE(IW,"(I8,F13.6,F12.6,F12.6)") IST,Q12D,Q12A,DQ12
               END IF
            END DO
            IF(MASWRK) THEN
               WRITE(IW,"(/15X,'FRAGMENT CHARGES OF "//
     *                  "TRANSITION DENSITIES AND')")
               WRITE(IW,"(12X,'FCD COUPLINGS BETWEEN "//
     *                  "GROUND AND SINGLET EXCITED STATES')")
               WRITE(IW,"(/4X,'STATES',5X,'Q12(D)',5X,'Q12(A)',7X,"//
     *                  "'DQ12',3X,'COUPLING(EV)',4X,"//
     *                  "'FCD-DQ1',4X,'FCD-DQ2',12X,'EL1',12X,'EL2')")
            END IF
            IST=0
            DO JST=1,NDSR
               Q02D = X(LQ02D + (JST-1))
               Q02A = X(LQ02A + (JST-1))
               DQ02 = Q02D - Q02A
               DQ2  =       X(LQ12D + (JST-1)*NDSR + (JST-1))
               DQ2  = DQ2 - X(LQ12A + (JST-1)*NDSR + (JST-1))
               E2   = X(LEE + JST-1)
               IF(ABS(DQ0-DQ2).GT.TINY(DQ0)) THEN
                  COUPL=E2*ABS(DQ02)/
     *                  SQRT((DQ0-DQ2)**2+4*DQ02**2)
                  TAN2TH  = 2*DQ02 / (DQ0-DQ2)
                  TANTH   = TAN2TH / (ONE + SQRT(1+TAN2TH**2))
                  COSSQTH = ONE / (ONE + TANTH**2)
                  SINSQTH = TANTH**2 * COSSQTH
                  EL1 = E2*SINSQTH
                  EL2 = E2*COSSQTH
C                 SINCOSTH = TANTH * COSSQTH
C                 COUPL = E2*SINCOSTH
               ELSE
                  COUPL = ZERO
                  EL1 = ZERO
                  EL2 = E2
               ENDIF
               IF(MASWRK) THEN
                  WRITE(IW,"(I5,I5,"//
     *                     "F11.6,F11.6,F11.6,"//
     *                     "F15.7,"//
     *                     "F11.6,F11.6,"//
     *                     "F15.7,F15.7)")
     *                     IST,JST,Q02D,Q02A,DQ02,COUPL*TOEV,DQ0,DQ2,
     *                     EL1*TOEV,EL2*TOEV
               END IF
            END DO
            IF(MASWRK) THEN
               WRITE(IW,"(/15X,'FRAGMENT CHARGES OF "//
     *                  "TRANSITION DENSITIES AND')")
               WRITE(IW,"(12X,'FCD COUPLINGS BETWEEN "//
     *                  "SINGLET EXCITED STATES')")
               WRITE(IW,"(/4X,'STATES',5X,'Q12(D)',5X,'Q12(A)',7X,"//
     *                  "'DQ12',3X,'COUPLING(EV)',4X,"//
     *                  "'FCD-DQ1',4X,'FCD-DQ2',12X,'EL1',12X,'EL2')")
            END IF
            DO IST=1,NDSR
               DO JST=1,NDSR
                  IF(IST.GE.JST) CYCLE
                  Q12D = X(LQ12D + (IST-1)*NDSR + (JST-1))
                  Q12A = X(LQ12A + (IST-1)*NDSR + (JST-1))
                  DQ12 = Q12D - Q12A
                  DQ1  =       X(LQ12D + (IST-1)*NDSR + (IST-1))
                  DQ1  = DQ1 - X(LQ12A + (IST-1)*NDSR + (IST-1))
                  DQ2  =       X(LQ12D + (JST-1)*NDSR + (JST-1))
                  DQ2  = DQ2 - X(LQ12A + (JST-1)*NDSR + (JST-1))
                  E1   = X(LEE + IST-1)
                  E2   = X(LEE + JST-1)
                  IF(ABS(DQ1-DQ2).GT.TINY(DQ1)) THEN
                     COUPL=(E2-E1)*ABS(DQ12)/
     *                     SQRT((DQ1-DQ2)**2+4*DQ12**2)
                     TAN2TH  = 2*DQ12 / (DQ2-DQ1)
                     TANTH   = TAN2TH / (ONE + SQRT(1+TAN2TH**2))
                     COSSQTH = ONE / (ONE + TANTH**2)
                     SINSQTH = TANTH**2 * COSSQTH
                     EL1 = E1*COSSQTH + E2*SINSQTH
                     EL2 = E1*SINSQTH + E2*COSSQTH
C                    SINCOSTH = TANTH * COSSQTH
C                    COUPL = (E1-E2)*SINCOSTH
                  ELSE
                     COUPL=ZERO
                     EL1 = E1
                     EL2 = E2
                  ENDIF
                  IF(MASWRK) THEN
                     WRITE(IW,"(I5,I5,"//
     *                        "F11.6,F11.6,F11.6,"//
     *                        "F15.7,"//
     *                        "F11.6,F11.6,"//
     *                        "F15.7,F15.7)")
     *                        IST,JST,Q12D,Q12A,DQ12,COUPL*TOEV,DQ1,DQ2,
     *                        EL1*TOEV,EL2*TOEV
                  END IF
               END DO
            END DO
            IF(MASWRK) THEN
               WRITE(IW,"(/1X,'...... END OF FCD CALCULATION ......'/)")
            END IF
            IF(MASWRK) THEN
               WRITE(IW,"(10X,'----------------------------------"//
     *                  "-----------')")
               WRITE(IW,"(10X,'MULTI-FED-FCD ELECTRONIC-COUPLING "//
     *                  "CALCULATION')")
               IF(NFG.NE.0) THEN
                  WRITE(IW,"(10X,'IN FMO, "//
     *                     "DONOR=',I5,' ACCEPTOR=',I5)") jcurfg,icurfg
               END IF
               WRITE(IW,"(10X,'----------------------------------"//
     *                  "-----------')")
               WRITE(IW,"(/6X,'THE MULTI-FED-FCD MATRIX IS (EV)')")
               CALL PRTRI(X(LHFIN),NDSR)
               WRITE(IW,"(/1X,'CT1 STATES:',I5)") NDSR-ILECT-ICT2CT1+2
               WRITE(IW,"(1X,'CT2 STATES:',I5)") ICT2CT1-1
               WRITE(IW,"(1X,'LE1 STATES:',I5)") ILECT-ILE2LE1
               WRITE(IW,"(1X,'LE2 STATES:',I5)") ILE2LE1-1
               WRITE(IW,
     *               "(/6X,"//
     *               "'SUMMARY OF MULTI-FED-FCD LOCALIZED STATES')")
               WRITE(IW,
     *               "(/2X,'STATE',2X,'EXCITATION',"//
     *               "2X,'TRANSITION DIPOLE, A.U.',2X,'OSCILLATOR'/"//
     *               "13X,'eV',10X,'X',7X,'Y',7X,'Z',5X,'STRENGTH')")
               DO IST=1,NDSR
                  E1 = X(LHFIN+IST*(IST+1)/2-1) / TOEV
                  DX = X(LTXYZMF+3*(IST-1))
                  DY = X(LTXYZMF+3*(IST-1)+1)
                  DZ = X(LTXYZMF+3*(IST-1)+2)
                  OS =      DX**2 * E1 * TWO/THREE
                  OS = OS + DY**2 * E1 * TWO/THREE
                  OS = OS + DZ**2 * E1 * TWO/THREE
                  WRITE(IW,'(2X,I3,4X,F8.3,3X,3F8.4,1X,F8.3)')
     *                     IST,E1*TOEV,DX,DY,DZ,OS
               END DO
               WRITE(IW,"(/6X,'MIXING COEFFICIENTS FOR MULTI-FED-FCD "//
     *                  "LOCALIZED STATES')")
               CALL PRSQ(X(LBCT1D),NDSR,NDSR,NDSR)
               WRITE(IW,"(/6X,'ANALYSIS FOR MULTI-FED-FCD "//
     *                  "LOCALIZED STATES')")
               WRITE(IW,
     *               "(/2X,'STATE',2X,'EXCITATION',5X,"//
     *               "'X(D)',4X,'X(A)',4X,'Q(D)',4X,'Q(A)')")
               DO IST=1,NDSR
                  WRITE(IW,'(2X,I3,4X,F8.3,3X,4F8.4)')
     *                     IST,
     *                     X(LHFIN+IST*(IST+1)/2-1),
     *                     X(LX12DL+IST-1),
     *                     X(LX12AL+IST-1),
     *                     X(LQ12DL+IST-1),
     *                     X(LQ12AL+IST-1)
               END DO
               IF(NFG.NE.0) THEN
                  IF(DOMULTIST) THEN
                     WRITE(IW,"(/1X,'FED search found',I3,"//
     *                        "' donor (need',I3,') and',I3,"//
     *                        "' acceptor (need',I3,') states.')")
     *                        NDON,JCURFGST,NACC,ICURFGST
                     WRITE(IW,"(1x,A8,' states:',100I3)")
     *                        'Donor   ',(IXFTCH(X(LDONI),I),I=1,NDON)
                     WRITE(IW,"(1x,A8,' states:',100I3)")
     *                        'Acceptor',(IXFTCH(X(LACCI),I),I=1,NACC)
                  ELSE
                     IF(ICURFGST.NE.0) THEN
                        WRITE(IW,"(/1X,'ASSIGNING STATE',I5,"//
     *                           "' TO FMO FRAGMENT',I5)")
     *                           ICURFGST,ICURFG
                     ELSE
                        WRITE(IW,"(/1X,'NO MATCHING STATE FOR "//
     *                            "FMO FRAGMENT',I5)") ICURFG
                     END IF
                     IF(JCURFGST.NE.0) THEN
                        WRITE(IW,"(1X,'ASSIGNING STATE',I5,"//
     *                           "' TO FMO FRAGMENT',I5)")
     *                           JCURFGST,JCURFG
                     ELSE
                        WRITE(IW,"(1X,'NO MATCHING STATE FOR "//
     *                            "FMO FRAGMENT',I5)") JCURFG
                     END IF
                     IF(ICURFGST.EQ.0.OR.JCURFGST.EQ.0) THEN
                        CALL ABRTX("STATE REQUESTED IN IACTFG(1) "//
     *                             "CAN NOT BE MATCHED.")
                     END IF
                  END IF
               END IF
               WRITE(IW,"(/1X,"//
     *                  "'... END OF MULTI-FED-FCD CALCULATION ...'/)")
            END IF
         END IF
C
C        --- FREE MEMORY FOR MULTI-FED-FCD RESULTS ---
C
         CALL RETFM(NEED1)
C
C        --- FREE MEMORY FOR FED/FCD RESULTS ---
C
         CALL RETFM(NEED)
      END IF
      RETURN
      END
C*MODULE TDUPRP *DECK FEDFMOSUM
C>
C>    @brief   Split properties between two fragments.
C>
C>    @details Store FED results.
C>
C>    @author Dmitri Fedorov
C>
      SUBROUTINE FEDFMOSUM(IOPT,PROPV,PROPD,PROPA)
      USE MX_LIMITS, ONLY: MXATM
      IMPLICIT NONE
      DOUBLE PRECISION PROPV(*),PROPD,PROPA,ZAN,C,espscf,e0scf,emp2s,X,A
      INTEGER IOPT,I,iatfg,nfg,nlayer,natfmo,nbdfg,naotyp,
     *        nbody,NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,IAN,IDAFMO,icurfg,
     *        jcurfg,kcurfg,icurlay,icurunt,nat1e,ncursh,ngau,icurpop,
     *        ifmostp,moncor,needr,modrst,norbproj,nunesp,iskipesp,
     *        IESDPPC,idoprop,mp2run,icurit,idmfmo,iddfmo,iddcur,nddleft
     *       ,ivmfmo,nzmtfmo,ifmobas,itmfmo,ixftch,lichfg,lmulfg,lidmrec
     *       ,lfrgnam,llayfrg,lindat,lnCBS,lfmozan,lfmoc,lfmomas,lizbas,
     *        liaglob,libdgh,liabdfg,ljabdfg,lnCAO,lidxCAO,liaprjo,
     *        ljaprjo,lCoreAO,lOccCor,lshiftb,liodfmo,lfmoda,lfmodb,
     *        lfmoespa,lfmoespb,llocfmo,lscffrg,lfmoscf,lrij,
     *        lpopmul,lpopmat,lialoc,lindbd,liatfrg,lindfrg,
     *        lindgfrg,lnatfrg,lnat0frg,lianfrg,lzanfrg,lcfrg,
     *        llibish,llibnsh,llibng,lindatg,lfmobuf,lfmode,
     *        lnumfrg,lloctat,liaoglob,lloadm,lfmoge,ldgrid,
     *        liodcfmo,ljob2grp,lfmopg,lemocdr,luntxyz,luntrot,
     *        lstonep,lmapsu,lfrgmul,lclmo,lialmo,lindlmo,
     *        latclmo,llmobdf,lfgflmo,lnfglmo,llfglmo,lpfglmo,
     *        LPOPDMAT,lidmpnt,liddpnt,livmpnt,liactfg,lcrfrg,
     *        lzlmfrgv,lYlmfrgv,lndtfrg,lf_mm,lg_mm,lmaxl30,
     *        libuffg,lindatp,lsmon,lsdim,ledimfed,lexcit3d,lifgfret,
     *        iat,nsegm
      COMMON /FMCOM / X(1)
      Common /fmoinf/ nfg,nlayer,natfmo,nbdfg,naotyp,nbody,nsegm
      common /fmopnt/ lichfg,lmulfg,lidmrec,lfrgnam,llayfrg,lindat,
     *                lnCBS,lfmozan,lfmoc,lfmomas,lizbas,liaglob,libdgh,
     *                liabdfg,ljabdfg,lnCAO,lidxCAO,liaprjo,ljaprjo,
     *                lCoreAO,lOccCor,lshiftb,liodfmo,lfmoda,lfmodb,
     *                lfmoespa,lfmoespb,llocfmo,lscffrg,lfmoscf,lrij,
     *                lpopmul,lpopmat,lialoc,lindbd,liatfrg,lindfrg,
     *                lindgfrg,lnatfrg,lnat0frg,lianfrg,lzanfrg,lcfrg,
     *                llibish,llibnsh,llibng,lindatg,lfmobuf(3),lfmode,
     *                lnumfrg,lloctat,liaoglob,lloadm,lfmoge,ldgrid,
     *                liodcfmo,ljob2grp,lfmopg,lemocdr,luntxyz,luntrot,
     *                lstonep,lmapsu,lfrgmul,lclmo,lialmo,lindlmo,
     *                latclmo,llmobdf,lfgflmo,lnfglmo,llfglmo,lpfglmo,
     *                LPOPDMAT,lidmpnt,liddpnt,livmpnt,liactfg,lcrfrg,
     *                lzlmfrgv,lYlmfrgv,lndtfrg,lf_mm,lg_mm,lmaxl30,
     *                libuffg,lindatp,lsmon,lsdim,ledimfed,lexcit3d,
     *                lifgfret
      common /fmorun/ espscf,e0scf(2),emp2s,IDAFMO,icurfg,jcurfg,kcurfg,
     *                icurlay,icurunt,nat1e,ncursh,ngau,icurpop,ifmostp,
     *                moncor,needr,modrst,norbproj,nunesp,iskipesp,
     *                IESDPPC,idoprop,mp2run,icurit,idmfmo,iddfmo,
     *                iddcur,nddleft,ivmfmo,nzmtfmo,ifmobas,itmfmo(2)
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
C
c     IOPT=0 initialize and add (put)
c     IOPT=1 add
C
      A=1.0D+00
      IF(IOPT.eq.0) THEN
         PROPD=0
         PROPA=0
         A=-1.0D+00
      ENDIF
      DO I=1,NAT
         iat=ixftch(x(liaglob),i)
         iatfg=ixftch(x(lindat),iat)
c        if(iatfg.eq.icurfg.and.iatfg.eq.jcurfg) THEN
c           PROP1 = PROP1 + PROPV(I)/2
c           PROP2 = PROP2 + PROPV(I)/2
c           Equal split for now: it is better to split as IAN/ZAN.
c        ELSE
            if(iatfg.eq.icurfg) PROPA = PROPA + A*PROPV(I)
            if(iatfg.eq.jcurfg) PROPD = PROPD + A*PROPV(I)
c        The case of shared atoms (BDA) is treated by assigning BDA
c        values in a connected dimer involving BDA to the BDA fragment.
c        This may not be the best but it is the simplest.
c        BDA should contribute little to an excitation; if it is not so,
c        the bond should not be detached.
c        ENDIF
      END DO
      RETURN
      END
C*MODULE TDUPRP *DECK FEDFMOMATCH
C>
C>    @brief   Match FMO1-TDDFT states to Multi-FED-FCD states.
C>
C>    @author Christian Friedl
C>
C>    @date   Dec, 2022 - Christian Friedl
C>    - Added subroutine
C>    @date   Feb, 2023 - Dmitri Fedorov
C>    - Added multi-state handling
C>
      SUBROUTINE FEDFMOMATCH(NDSR,KDON,KACC,HFIN,X12DL,X12AL,
     *              ILECT,ILE2LE1,ICURFGST,JCURFGST,
     *              domultist,ndon,nacc)
C
      USE MX_LIMITS, ONLY: MXRT
      IMPLICIT NONE
C     Declarations of arguments
      INTEGER, INTENT(IN) :: NDSR
      INTEGER, INTENT(INOUT), DIMENSION(NDSR) :: KDON(NDSR)
      INTEGER, INTENT(INOUT), DIMENSION(NDSR) :: KACC(NDSR)
      DOUBLE PRECISION, INTENT(IN), DIMENSION(NDSR*NDSR) :: HFIN
      DOUBLE PRECISION, INTENT(IN), DIMENSION(NDSR) :: X12DL
      DOUBLE PRECISION, INTENT(IN), DIMENSION(NDSR) :: X12AL
      INTEGER, INTENT(IN) :: ILECT,ILE2LE1
      INTEGER, INTENT(OUT) :: ICURFGST,JCURFGST
      LOGICAL, INTENT(OUT) :: domultist
      INTEGER, INTENT(OUT) :: ndon,nacc
C     Declarations of common blocks
      DOUBLE PRECISION ENUCR,EELCT,ETOT,SZ,SZZ,ECORE,ESCF,EERD,E1,E2
      DOUBLE PRECISION VEN,VEE,EPOT,EKIN,ESTATE,STATN,EDFT,EDISP
      DOUBLE PRECISION X
      DOUBLE PRECISION espsca,RESPAP,rESPPC,rESDIM
      DOUBLE PRECISION restri,rcorsd,respct,convfg,cnvdmp,coroff,rflmo
      DOUBLE PRECISION orshft,orshft2,cnvafo,ascreen
      INTEGER IXESP,mxitfg
      INTEGER nguess,NBSSE,modorb,modpar,irststp,irstlay,nprfmo
      INTEGER nfmopal,modprp,maxl1c,ipieda,modgrd,modesp,ivmul
      INTEGER modlmo,nopden,mofock,modfd,modfmm,ncentm,ndualb
      INTEGER ngab,modpan
      INTEGER lichfg,lmulfg,lidmrec,lfrgnam,llayfrg,lindat
      INTEGER lnCBS,lfmozan,lfmoc,lfmomas,lizbas,liaglob,libdgh
      INTEGER liabdfg,ljabdfg,lnCAO,lidxCAO,liaprjo,ljaprjo
      INTEGER lCoreAO,lOccCor,lshiftb,liodfmo,lfmoda,lfmodb
      INTEGER lfmoespa,lfmoespb,llocfmo,lscffrg,lfmoscf,lrij
      INTEGER lpopmul,lpopmat,lialoc,lindbd,liatfrg,lindfrg
      INTEGER lindgfrg,lnatfrg,lnat0frg,lianfrg,lzanfrg,lcfrg
      INTEGER llibish,llibnsh,llibng,lindatg,lfmobuf,lfmode
      INTEGER lnumfrg,lloctat,liaoglob,lloadm,lfmoge,ldgrid
      INTEGER liodcfmo,ljob2grp,lfmopg,lemocdr,luntxyz,luntrot
      INTEGER lstonep,lmapsu,lfrgmul,lclmo,lialmo,lindlmo
      INTEGER latclmo,llmobdf,lfgflmo,lnfglmo,llfglmo,lpfglmo
      INTEGER LPOPDMAT,lidmpnt,liddpnt,livmpnt,liactfg,lcrfrg
      INTEGER lzlmfrgv,lYlmfrgv,lndtfrg,lf_mm,lg_mm,lmaxl30
      INTEGER libuffg,lindatp,lsmon,lsdim,ledimfed,lexcit3d,
     *        lifgfret
      DOUBLE PRECISION espscf,e0scf,emp2s
      INTEGER IDAFMO,icurfg,jcurfg,kcurfg
      INTEGER icurlay,icurunt,nat1e,ncursh,ngau,icurpop,ifmostp
      INTEGER moncor,needr,modrst,norbproj,nunesp,iskipesp
      INTEGER IESDPPC,idoprop,mp2run,icurit,idmfmo,iddfmo
      INTEGER iddcur,nddleft,ivmfmo,nzmtfmo,ifmobas,itmfmo
C     Other declarations
      INTEGER IDONMX,IACCMX,IDON,IACC,ILE1,ILE2,IST,JST,KST
      INTEGER ILE1B,ILE1E,ILE2B,ILE2E
      INTEGER i,j,ii,jj,ih,jh,ii0,jj0,loop,looph
      DOUBLE PRECISION, PARAMETER :: TOEV = 27.21138386D+00
      INTEGER, EXTERNAL :: IXFTCH
C
      COMMON /ENRGYS/ ENUCR,EELCT,ETOT,SZ,SZZ,ECORE,ESCF,EERD,E1,E2,
     *                VEN,VEE,EPOT,EKIN,ESTATE(MXRT),STATN,EDFT(3),EDISP
      COMMON /FMCOM / X(1)
      common /fmoopt/ espsca(9),RESPAP(2),rESPPC(2),rESDIM,restri(4),
     *                rcorsd,respct,convfg,cnvdmp,coroff,rflmo(4),
     *                orshft,orshft2,cnvafo,ascreen(4),IXESP,mxitfg,
     *                nguess,NBSSE,modorb,modpar,irststp,irstlay,nprfmo,
     *                nfmopal,modprp,maxl1c,ipieda,modgrd,modesp,ivmul,
     *                modlmo,nopden,mofock,modfd,modfmm,ncentm,ndualb,
     *                ngab,modpan
      common /fmopnt/ lichfg,lmulfg,lidmrec,lfrgnam,llayfrg,lindat,
     *                lnCBS,lfmozan,lfmoc,lfmomas,lizbas,liaglob,libdgh,
     *                liabdfg,ljabdfg,lnCAO,lidxCAO,liaprjo,ljaprjo,
     *                lCoreAO,lOccCor,lshiftb,liodfmo,lfmoda,lfmodb,
     *                lfmoespa,lfmoespb,llocfmo,lscffrg,lfmoscf,lrij,
     *                lpopmul,lpopmat,lialoc,lindbd,liatfrg,lindfrg,
     *                lindgfrg,lnatfrg,lnat0frg,lianfrg,lzanfrg,lcfrg,
     *                llibish,llibnsh,llibng,lindatg,lfmobuf(3),lfmode,
     *                lnumfrg,lloctat,liaoglob,lloadm,lfmoge,ldgrid,
     *                liodcfmo,ljob2grp,lfmopg,lemocdr,luntxyz,luntrot,
     *                lstonep,lmapsu,lfrgmul,lclmo,lialmo,lindlmo,
     *                latclmo,llmobdf,lfgflmo,lnfglmo,llfglmo,lpfglmo,
     *                LPOPDMAT,lidmpnt,liddpnt,livmpnt,liactfg,lcrfrg,
     *                lzlmfrgv,lYlmfrgv,lndtfrg,lf_mm,lg_mm,lmaxl30,
     *                libuffg,lindatp,lsmon,lsdim,ledimfed,lexcit3d,
     *                lifgfret
      common /fmorun/ espscf,e0scf(2),emp2s,IDAFMO,icurfg,jcurfg,kcurfg,
     *                icurlay,icurunt,nat1e,ncursh,ngau,icurpop,ifmostp,
     *                moncor,needr,modrst,norbproj,nunesp,iskipesp,
     *                IESDPPC,idoprop,mp2run,icurit,idmfmo,iddfmo,
     *                iddcur,nddleft,ivmfmo,nzmtfmo,ifmobas,itmfmo(2)
C
C     Create two lists of state indices - each with a
C     maximum of NDSR elements. One list KDON is for
C     states located on the donor, another one KACC for
C     states located on the acceptor. Go through both LE1 and
C     LE2 states, check whether they are located on either
C     donor or acceptor (via the expectation value  of X(D) and
C     X(A)) and add them sorted by energy to one of the two
C     lists.
C
      IF(ILECT.LE.1) THEN
         ICURFGST=0
         JCURFGST=0
         RETURN
      END IF
      domultist=iand(modprp,262144).ne.0
      if(domultist) then
c        IST (JST) is the number of states in IFG (JFG).
c        IFG is acceptor, JFG is donor
         IST=IXFTCH(X(LIACTFG),ICURFG)
         JST=IXFTCH(X(LIACTFG),JCURFG)
c        copy diagonal of H to a work array
         do i=1,ndsr
           x(ledimfed+i-1)=HFIN((i*i+i)/2)
         enddo
         call rindsort(0,ndsr,x(ledimfed),kdon)
c        save a copy of the index array
         CALL iCOPY(NDSR,kdon,1,x(ledimfed),1)
         nacc=0
         ndon=0
c        loop back because sort is in the decreasing order
         do i=ndsr,1,-1
c           j is the state number in FED
            j=IXFTCH(X(ledimfed),i)
            IF(ABS(X12DL(j)).GE.ABS(X12AL(j))) THEN
c              found a donor state
               ndon=ndon+1
               kdon(ndon)=j
            else
c              found an acceptor state
               nacc=nacc+1
               kacc(nacc)=j
            endif
         enddo
         ICURFGST=IST
         JCURFGST=JST
         if(ndon.lt.JST.or.nacc.lt.IST) then
            call abrtx("Failure in FED search")
         endif
c        i=ist+jst
c        call vclr(x(ledimfed),1,(i*i+i)/2)
c        I-I block
         do i=1,ist
            ii=i+ixftch(x(lifgfret),icurfg)-1
            do j=1,i
               jj=j+ixftch(x(lifgfret),icurfg)-1
               ih=max(kacc(i),kacc(j))
               jh=min(kacc(i),kacc(j))
               loop=(ii*ii-ii)/2+jj
               looph=(ih*ih-ih)/2+jh
               x(lexcit3d+loop-1)=x(lexcit3d+loop-1)+HFIN(looph)/TOEV
c              x(ledimfed+loop-1)=HFIN(looph)
            enddo
         enddo
c        J-J block
         do i=1,jst
            ii=i+ixftch(x(lifgfret),jcurfg)-1
            do j=1,i
               jj=j+ixftch(x(lifgfret),jcurfg)-1
               ih=max(kdon(i),kdon(j))
               jh=min(kdon(i),kdon(j))
               loop=(ii*ii-ii)/2+jj
               looph=(ih*ih-ih)/2+jh
               x(lexcit3d+loop-1)=x(lexcit3d+loop-1)+HFIN(looph)/TOEV
c              x(ledimfed+loop-1)=HFIN(looph)
            enddo
         enddo
c        J-I block
         loop=0
         do i=1,jst
            ii0=i+ixftch(x(lifgfret),jcurfg)-1
            do j=1,ist
               jj0=j+ixftch(x(lifgfret),icurfg)-1
               ii=max(ii0,jj0)
               jj=min(ii0,jj0)
               ih=max(kdon(i),kacc(j))
               jh=min(kdon(i),kacc(j))
               loop=(ii*ii-ii)/2+jj
               looph=(ih*ih-ih)/2+jh
               x(lexcit3d+loop-1)=x(lexcit3d+loop-1)+HFIN(looph)/TOEV
c              x(ledimfed+loop-1)=HFIN(looph)
            enddo
         enddo
c        CALL PRTRI(X(ledimfed),ist+jst)
      else
         IDON=1
         IACC=1
         ILE1B=NDSR-ILECT+2
         ILE1E=NDSR-ILE2LE1+1
         ILE2B=NDSR-ILE2LE1+2
         ILE2E=NDSR
         ILE1=ILE1B
         ILE2=ILE2B
         DO WHILE(ILE1.LE.ILE1E.OR.ILE2.LE.ILE2E)
C
C          Pick the lower energy state
C          The lower energy state is saved in IST,
C          the higher energy state is saved in JST.
C
           IF(ILE1.GT.ILE1E) THEN
              IST=ILE2
              JST=0
              ILE2=ILE2+1
           ELSE IF(ILE2.GT.ILE2E) THEN
              IST=ILE1
              JST=0
              ILE1=ILE1+1
           ELSE
             IF(HFIN(ILE1*(ILE1+1)/2).LE.HFIN(ILE2*(ILE2+1)/2)) THEN
                IST=ILE1
                JST=ILE2
             ELSE
                IST=ILE2
                JST=ILE1
             END IF
             ILE1=ILE1+1
             ILE2=ILE2+1
           END IF
C
C          Find out whether state IST
C          belongs to donor or acceptor and
C          append it to the appropriate list
C
           IF(ABS(X12DL(IST)).GE.ABS(X12AL(IST))) THEN
              KDON(IDON)=IST
              IDON=IDON+1
           ELSE
              KACC(IACC)=IST
              IACC=IACC+1
           END IF
C
C          Find out whether state JST
C          belongs to donor or acceptor and
C          append it to the appropriate list
C
           IF(JST.NE.0) THEN
             IF(ABS(X12DL(JST)).GE.ABS(X12AL(JST))) THEN
                KDON(IDON)=JST
                IDON=IDON+1
              ELSE
                KACC(IACC)=JST
                IACC=IACC+1
              END IF
           END IF
         END DO
C
C        Match the states to the ones requested IACTFG(1)
C
         IST=IXFTCH(X(LIACTFG),ICURFG)
         IF(IST.LT.IACC) THEN
            IST=KACC(IST)
         ELSE
            IST=0
         END IF
         JST=IXFTCH(X(LIACTFG),JCURFG)
         IF(JST.LT.IDON) THEN
            JST=KDON(JST)
         ELSE
            JST=0
         END IF
         ICURFGST=IST
         JCURFGST=JST
         nacc=1
         ndon=1
         IF(ICURFGST.EQ.0) nacc=0
         IF(JCURFGST.EQ.0) ndon=0
         IF(ICURFGST.EQ.0.OR.JCURFGST.EQ.0) RETURN
C
C        Write states energies into EDIMFED(1), EDIMFED(2)
C        and their coupling into EDIMFED(3)
C
         x(ledimfed)=HFIN((IST*(IST+1))/2)/TOEV
         x(ledimfed+1)=HFIN((JST*(JST+1))/2)/TOEV
         KST=MAX(IST,JST)
         KST=(KST*(KST-1))/2+MIN(IST,JST)
         x(ledimfed+2)=HFIN(KST)/TOEV
      endif
      RETURN
      END
