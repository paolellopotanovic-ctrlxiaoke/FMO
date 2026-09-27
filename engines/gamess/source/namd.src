C*MODULE NAMD   *DECK VVERTULLY
C>
C>     @brief Velocity verlet propagator
C>            with Tully's fewest switches surface hopping
C>
C>     @author Seunghoon Lee
C>
      SUBROUTINE VVERTULLY(ISTEP,TDECOE,RANDOM,NDNINT,INTMD,NPROP)
      USE comm_NONAD, ONLY: NAMD,NDSWCH,NDRST,NDTLF,COLD,THRSHE
      USE mx_limits, ONLY: mxatm, mxfrg
      USE prec, ONLY: dp
      IMPLICIT NONE
C
        REAL(KIND=dp), DIMENSION(MXFRG) :: OX, OY,  OZ,  QW, QW1, QW2,  &
     &                                     QX, QX1, QX2, QY, QY1, QY2,  &
     &                                     QZ, QZ1, QZ2, VX, VY,  VZ
        REAL(KIND=dp), DIMENSION(MXATM) :: VXQM, VYQM, VZQM
      COMMON /ATMDAT/ VX, VY, VZ, QW, QX, QY, QZ, OX, OY, OZ, QW1, QX1, &
     &                QY1, QZ1, QW2, QX2, QY2, QZ2, VXQM, VYQM, VZQM
        REAL(KIND=dp), DIMENSION(1) :: X
      COMMON /FMCOM / X
        REAL(KIND=dp) :: E
        REAL(KIND=dp), DIMENSION(3*MXATM) :: EG
      COMMON /FUNCT / E, EG
        REAL(KIND=dp), DIMENSION(3,MXATM) :: C
        INTEGER, DIMENSION(MXATM) :: IAN
        INTEGER :: ICH, MUL, NA, NAT, NB, NE, NQMT, NUM
        REAL(KIND=dp), DIMENSION(MXATM) :: ZAN
      COMMON /INFOA / NAT, ICH, MUL, NUM, NQMT, NE, NA, NB, ZAN, C, IAN
        LOGICAL :: ALPHKWD, BETAKWD, MRDEA, MREKT, SG1T, TAMMD, TPA,    &
     &             TRIPLET
        REAL(KIND=dp) :: CNVTOL
        INTEGER :: IRECTD, ITDFG, ITDPRP, JANST, MAXVEC, MODTD, MTHST,  &
     &             MULTD, NLEBT, NONEQR, NPHIT, NRADT, NSTAT, NTHET,    &
     &             NTHST, NTRIAL
        INTEGER, DIMENSION(4) :: IFEDAT
        REAL(KIND=dp), DIMENSION(2) :: PFREQ
        REAL(KIND=dp), DIMENSION(3) :: SPCP
      COMMON /INFOTD/ CNVTOL, PFREQ, MODTD, JANST, NRADT, NTHET, NPHIT, &
     &                NLEBT, NSTAT, NTRIAL, MAXVEC, NTHST, IRECTD,      &
     &                ITDFG, ITDPRP, TRIPLET, SG1T, NONEQR, TAMMD, TPA, &
     &                ALPHKWD, BETAKWD, SPCP, MULTD, MREKT, MRDEA,      &
     &                MTHST, IFEDAT
        INTEGER :: IDAF, IJKO, IJKT, IP, IR, IW, NAV
        INTEGER, DIMENSION(950) :: IODA
      COMMON /IOFILE/ IR, IW, IP, IJKO, IJKT, IDAF, NAV, IODA
        REAL(KIND=dp), DIMENSION(MXATM) :: ZMASS
      COMMON /MASSES/ ZMASS
        INTEGER :: NSTEPS
      COMMON /MDSIM4/ NSTEPS
        LOGICAL :: DSKWRK, GOPARR, MASWRK
        INTEGER :: IBTYP, IPTIM, MASTER, ME, NPROC
      COMMON /PAR   / ME, MASTER, NPROC, IBTYP, IPTIM, GOPARR, DSKWRK,  &
     &                MASWRK
C
      INTEGER :: INTMD, ISTEP, NDNINT, NPROP
      REAL(KIND=dp), DIMENSION(NSTEPS+1) :: RANDOM
      REAL(KIND=dp), DIMENSION(2,NSTAT) :: TDECOE

      INTEGER, DIMENSION(NAT, 2) :: ICOBONDS, ICNBONDS, ICCBONDS
      REAL(KIND=dp), DIMENSION(3,nat,3) :: RATPBC
C
      INTEGER :: I, IEGAO, IEGBO, ITMP, IVAOD, IVBOD, IOX, INIT, ICAR,  &
     &           LAST, LOADFM, NATCC, NEED3, NTHSTOLD
      LOGICAL :: LMRSF, MRSFS, MRSFT, NEWRND, DOPBC
      REAL(KIND=dp), DIMENSION(NSTEPS+1) :: RSRND
      REAL(KIND=dp), DIMENSION(2,NSTAT) :: TDECOEO
C
      LOGICAL, SAVE :: FIRST
      DATA FIRST/.TRUE./
C
C   SAVE REQUIRED OLD VALUES
c   SAVE OLD NTHST FOR DEBUGGING
      NTHSTOLD=NTHST
C   SAVE OLD COORDINATE FOR RESTARTING
      DO I=1,NSTAT
         TDECOEO(1,I) = TDECOE(1,I)
         TDECOEO(2,I) = TDECOE(2,I)
      ENDDO
C   SAVE OLD COORDINATES.
      CALL DCOPY(3*NAT,C,1,COLD,1)
C   SAVE OLD ALPHA BETA MO COEFFICIENTS, ORB ENERGIES
      CALL VALFM(LOADFM)
      IVAOD = LOADFM + 1
      IVBOD = IVAOD + NUM*NQMT
      IEGAO = IVBOD + NUM*NQMT
      IEGBO = IEGAO + NQMT
      LAST  = IEGBO + NQMT
      NEED3 = LAST  - LOADFM-1
      CALL GETFM(NEED3)
C
      CALL DAREAD(IDAF,IODA,X(IVAOD),NUM*NQMT,15 ,0)
      CALL DAREAD(IDAF,IODA,X(IVBOD),NUM*NQMT,19 ,0)
      CALL DAREAD(IDAF,IODA,X(IEGAO),NQMT    ,17 ,0)
      CALL DAREAD(IDAF,IODA,X(IEGBO),NQMT    ,21 ,0)
C
      CALL DAWRIT(IDAF,IODA,X(IVAOD),NUM*NQMT,700,0)
      CALL DAWRIT(IDAF,IODA,X(IVBOD),NUM*NQMT,701,0)
      CALL DAWRIT(IDAF,IODA,X(IEGAO),NQMT    ,702,0)
      CALL DAWRIT(IDAF,IODA,X(IEGBO),NQMT    ,703,0)
C
      CALL RETFM(NEED3)
C
C VELOCITY-VERLET PROPAGATION
      CALL VVERMOVE(INTMD,ISTEP,NPROP,VXQM,VYQM,VZQM,
     *              ZMASS,E,EG,C,IAN,NAT,NATCC,ZAN,COLD,ICOBONDS,
     *              ICNBONDS,ICCBONDS,IOX,INIT,ICAR,DOPBC,RATPBC)

      IF(MASWRK) WRITE(IW,1000)
      CALL TIMIT(1)
C
C SURFACE HOPPING ALGORITHM FOR SINGLET MRSF-TDDFT
      IF(.NOT.FIRST) CALL MRSFHP(ISTEP,TDECOE,TDECOEO,RANDOM,NDNINT,
     *                           NTHSTOLD)
      FIRST=.FALSE.
      IF(MASWRK) WRITE(IW,1100)
      CALL TIMIT(1)
C
C SAVE RESTART INFO
      DO I=1,NSTEPS+1
         RSRND(I)=0.D+00
      ENDDO
      DO I=1,NSTEPS+1-ISTEP
         ITMP=I+ISTEP
         RSRND(I)=RANDOM(ITMP)
      ENDDO
C
 1000 FORMAT( 1X,'.... DONE WITH VELOCITY-VERLET MOVE ....')
 1100 FORMAT(/1X,'.... DONE WITH TULLY''S SURFACE HOPPING ....')
      END SUBROUTINE VVERTULLY
C*MODULE NAMD   *DECK MRSFHP
C>
C>     @brief Tully's fewest switches surface hopping algorithm
C>            for MRSF-TDDFT
C>
C>     @author Seunghoon Lee
C>
      SUBROUTINE MRSFHP(ISTEP,TDECOE,TDECOEO,RANDOM,NDNINT,NTHSTOLD)
C
      USE comm_NONAD, ONLY: NAMD,NDSWCH,NDRST,NDTLF,COLD,THRSHE
      USE mx_limits, ONLY: mxatm, mxrt
      USE prec, ONLY: dp
      IMPLICIT NONE
C
        REAL(KIND=dp) :: AU2AMU, AU2ANG, AU2KCAL, AU2SEC, EK2BAR,       &
     &                   EK2KCAL
      COMMON /CONVMD/ AU2KCAL, EK2KCAL, AU2ANG, AU2SEC, AU2AMU, EK2BAR
        REAL(KIND=dp) :: E1, E2, ECORE, EDISP, EELCT, EERD, EKIN, ENUCR,&
     &                   EPOT, ESCF, ETOT, STATN, SZ, SZZ, VEE, VEN
        REAL(KIND=dp), DIMENSION(2) :: EDFT
        REAL(KIND=dp), DIMENSION(MXRT) :: ESTATE
      COMMON /ENRGYS/ ENUCR, EELCT, ETOT, SZ, SZZ, ECORE, ESCF, EERD,   &
     &                E1, E2, VEN, VEE, EPOT, EKIN, ESTATE, STATN, EDFT,&
     &                EDISP
        REAL(KIND=dp), DIMENSION(1) :: X
      COMMON /FMCOM / X
        REAL(KIND=dp), DIMENSION(3,MXATM) :: C
        INTEGER, DIMENSION(MXATM) :: IAN
        INTEGER :: ICH, MUL, NA, NAT, NB, NE, NQMT, NUM
        REAL(KIND=dp), DIMENSION(MXATM) :: ZAN
      COMMON /INFOA / NAT, ICH, MUL, NUM, NQMT, NE, NA, NB, ZAN, C, IAN
        LOGICAL :: ALPHKWD, BETAKWD, MRDEA, MREKT, SG1T, TAMMD, TPA,    &
     &             TRIPLET
        REAL(KIND=dp) :: CNVTOL
        INTEGER :: IRECTD, ITDFG, ITDPRP, JANST, MAXVEC, MODTD, MTHST,  &
     &             MULTD, NLEBT, NONEQR, NPHIT, NRADT, NSTAT, NTHET,    &
     &             NTHST, NTRIAL
        INTEGER, DIMENSION(4) :: IFEDAT
        REAL(KIND=dp), DIMENSION(2) :: PFREQ
        REAL(KIND=dp), DIMENSION(3) :: SPCP
      COMMON /INFOTD/ CNVTOL, PFREQ, MODTD, JANST, NRADT, NTHET, NPHIT, &
     &                NLEBT, NSTAT, NTRIAL, MAXVEC, NTHST, IRECTD,      &
     &                ITDFG, ITDPRP, TRIPLET, SG1T, NONEQR, TAMMD, TPA, &
     &                ALPHKWD, BETAKWD, SPCP, MULTD, MREKT, MRDEA,      &
     &                MTHST, IFEDAT
        INTEGER :: IDAF, IJKO, IJKT, IP, IR, IW, NAV
        INTEGER, DIMENSION(950) :: IODA
      COMMON /IOFILE/ IR, IW, IP, IJKO, IJKT, IDAF, NAV, IODA
        REAL(KIND=dp) :: BATHT, DT, EKINR, EKINT
        INTEGER :: MDTYP, PRESSMD
      COMMON /MDSIM / DT, BATHT, EKINT, EKINR, MDTYP, PRESSMD
        INTEGER :: NSTEPS
      COMMON /MDSIM4/ NSTEPS
        LOGICAL :: DSKWRK, GOPARR, MASWRK
        INTEGER :: IBTYP, IPTIM, MASTER, ME, NPROC
      COMMON /PAR   / ME, MASTER, NPROC, IBTYP, IPTIM, GOPARR, DSKWRK,  &
     &                MASWRK
C
      INTEGER :: ISTEP, NDNINT, NTHSTOLD
      REAL(KIND=dp), DIMENSION(NSTEPS+1) :: RANDOM
      REAL(KIND=dp), DIMENSION(2,NSTAT) :: TDECOE, TDECOEO
C
      REAL(KIND=dp) :: DCHECK, DDNINT, DTAU, RNDND
      INTEGER :: I, ICMHP, IHPPR, IINXN, IINXO, INCME, ISTAS, ISTEG,    &
     &           IXPYN, IXPYNT, IXPYO, IXPYOT, JSTATE, L1, L2, L3, L4,  &
     &           L4LR, L5, L6, L7, L7LR, L7MAX, LAST, LOADFM, LP, LPLR, &
     &           LX, MXVEC, NDSR, NEED, NOCA, NOCB, NTST, NVIRA, NVIRB
      LOGICAL :: LMRSF, MRSFS, MRSFT, MRSFQ, NEWRND
      REAL(KIND=dp), DIMENSION(NSTAT) :: TMP
C
      IF(MULTD.EQ.1) THEN
          MRSFS = .TRUE.
      ELSE IF(MULTD.EQ.3) THEN
          MRSFT = .TRUE.
      ELSE IF(MULTD.EQ.5) THEN
          MRSFQ = .TRUE.
      ELSE
          IF(MASWRK) WRITE (IW,*)
     *    'ERROR:: CANNOT FIND PROPER MRSF STATES FOR NAMD'
          CALL ABRT
      ENDIF

      IF(MASWRK) WRITE (IW,2000)
C
      LX    = NQMT
      NOCA  = NA
      NOCB  = NB
      NVIRA = LX-NOCA
      NVIRB = LX-NOCB
      L1    = NUM
      L3    = L1*L1
      L2    = (L3+L1)/2
      L7    = NOCA*NVIRB
      L7LR  = NOCA*NVIRB-1
      L7MAX = NOCA*(NUM-NOCB)
      NDSR  = NSTAT
      IF(NDSR.GT.L7) NDSR=L7
      NTST  = NDSR
      L5    = NTST*NTST
      MXVEC = MAXVEC*NDSR
      IF(MXVEC.GT.L7) MXVEC=L7
      L4    = L7*NDSR
      L4LR  = L7LR*NDSR
      DTAU  = DT/AU2SEC
C
      CALL DAREAD(IDAF,IODA,TMP,NDSR,706,0)
C   ALLOCATE STATE ENERGIES
      DO I = 1, NDSR
         TMP(I) = ESCF + TMP(I)
      ENDDO
C
      NDNINT = INT( ABS( DTAU*TMP(NDSR)/0.1 ) )
      DDNINT = REAL(NDNINT)
      L6 = NTST*NTST*NDNINT
C
C     NOCA : THE NUMBER OF OCCUPIED ALPHA SPATIAL ORBITALS (NA) : N+1
C     NOCB : THE NUMBER OF OCCUPIED BETA SPATIAL ORBITALS (NB) : N-1
C     NVIRA : THE NUMBER OF VIRTUAL ALPHA SPATIAL ORBITALS : V-1
C     NVIRB : THE NUMBER OF VIRTUAL BETA SPATIAL ORBITALS : V+1
C     L1 : THE NUMBER OF BASIS SET FUNCTIONS (NUM)
C     L2 : THE NUMBER OF SYMMETRIC KS OVERLAP INTEGRAL MATRIX
C     L3 : THE NUMBER OF UN-SYMMETRIC KS OVERLAP INTEGRAL MATRIX
C     L4 : THE NUMBER OF ELEMENTS OF X+Y WITH ALL STATES
C     L5 : THE NUMBER OF UN-SYMMETRIC STATE OVERLAP INTEGRAL MATRIX
C     L6 : THE NUMBER OF HOPPING PROBABILITY ARRAY AT EACH DTAU/NDNINT STEP.
C     L7 : THE NUMBER OF OCCUPIED ALPHA * VIRTUAL BETA SPATIAL ORBITALS
C     MXVEC : THE NUMVER OF EXCITATION ENERGIES
C     NDSR : THE NUMBER OF CONSIDERED EXCITED STATES IN SFDFT
C     NTST : THE NUMBER OF CONSIDERED TOTAL STATES IN SFDFT
C     DTAU : DT in atomic unit
C
C MEMORY ALLOCATE
      CALL VALFM(LOADFM)
      IXPYO = LOADFM + 1
      IXPYN = IXPYO + L4
      IXPYOT= IXPYN + L4
      IXPYNT= IXPYOT+ L4LR
      ISTAS = IXPYNT+ L4LR
      INCME = ISTAS + L5
      ISTEG = INCME + L5
      IHPPR = ISTEG + NTST
      ICMHP = IHPPR + L6
      IINXO = ICMHP + L5
      IINXN = IINXO + 2*L7
      LAST  = IINXN + 2*L7
      NEED  = LAST - LOADFM-1
      CALL GETFM(NEED)
C READ OLD AND NEW CSF COEFFICIENTS AND TDDFT EXCITATION ENERGIES
      CALL DAREAD(IDAF,IODA,X(IXPYOT),L4LR,704,0)
      CALL DAREAD(IDAF,IODA,X(IXPYNT),L4LR,705,0)
C GET DIMENSIONAL TRANSFORMED X AMPLITUDES OF MRSF-TDDFT
      LPLR = IXPYOT - L7LR
      LP   = IXPYO  - L7
      DO JSTATE=1,NDSR
         LPLR  = LPLR + L7LR
         LP    = LP   + L7
         call mrsfxvec(x(lplr),x(lp),noca,nocb,lx,mrsfs,mrsft)
      END DO
      LPLR = IXPYNT - L7LR
      LP   = IXPYN  - L7
      DO JSTATE=1,NDSR
         LPLR  = LPLR + L7LR
         LP    = LP   + L7
         call mrsfxvec(x(lplr),x(lp),noca,nocb,lx,mrsfs,mrsft)
      END DO
C
      CALL DCOPY(NTST,TMP,1,X(ISTEG),1)
C
C COMPUTE NUMERICAL NACVS BY FINITE DIFFERENCE METHOD
C     OVERLAP
      CALL MRSFOV(X(ISTAS),X(IXPYO),X(IXPYN),
     *            LX,NOCA,NOCB,L1,L2,L3,NDSR,L7)
      CALL NACVFD(X(INCME),X(ISTAS),DTAU,NDSR)
C
      IF(MASWRK) WRITE(IW,1000)
      CALL TIMIT(1)
C
C PROPAGATE TIME-DEPENDENT EXPANSION COEFFICIENTS
C AND COMPUTE HOPPING PROBABILITY
C IN ELECTRONIC SUB-TIME SCALE (DTAU/DDNINT)
C
      DO I=1,NDNINT
         CALL PPTDECOE(TDECOE,X(INCME),X(ISTEG),NDSR,DTAU/DDNINT)
         CALL FSSHPRST(X(IHPPR),I,TDECOE,X(INCME),NDSR,DTAU/DDNINT,
     *               NDNINT)
      ENDDO
C
      IF(MASWRK) WRITE(IW,1100)
      CALL TIMIT(1)
C
C COMPUTE HOPPING PROBABILITY IN NUCLEUS TIME SCALE
      CALL FSSHPR(X(ICMHP),X(IHPPR),NDSR,NDNINT)
C
C PRINT RESULTS
      DCHECK=DTAU/DDNINT*X(ISTEG+NDSR-1)
      RNDND=RANDOM(ISTEP+1)
      IF(MASWRK) CALL NAMDPRNT(X(ISTAS),X(INCME),X(ICMHP),DCHECK,
     *                         TDECOE,TDECOEO,NDTLF,
     *                         NDNINT,NDSR,RNDND,NTHST)
C
C DETERMINE SURFACE HOPPING BASED ON
C FEWEST SWITCHES SURFACE HOPPING ALGORITHM
      CALL FSSH(X(ICMHP),X(ISTEG),RNDND,NDSR)
C
C      MEMORY DEALLOCATE
      CALL RETFM(NEED)
C
 1000 FORMAT( 1X,'.... DONE WITH NONADIABATIC COUPLING TERMS ....')
 1100 FORMAT( 1X,'.... DONE WITH HOPPING PROBABILITIES ....')
 2000 FORMAT(/1X,71(1H-)/
     *   7X,'TULLY''S SURFACE HOPPING NONADIABATIC MOLECULAR DYNAMICS'/
     *               19X,'FOR MRSF-TDDFT CALCULATION WITH SCFTYP=ROHF'/
     *                                    40X,'CODED BY SEUNGHOON LEE'/
     *       1X,71(1H-))
      END SUBROUTINE MRSFHP
C*MODULE NAMD   *DECK MRSFOV
C>
C>     @brief Compute overlap integrals between MRSF response states
C>            at different MD time steps
C>
C>     @details Fast overlap evaluations using the TLF approximation
C>              introduced in JCTC 15 882 (2019)
C>
C>     @author Seunghoon Lee
C>
      SUBROUTINE MRSFOV(DSTAS,DXPYO,DXPYN,LX,NOCA,NOCB,L1,L2,L3,NDSR,L7)
C
      USE prec, ONLY: dp
      USE constants, ONLY: two
      IMPLICIT NONE
C
        REAL(KIND=dp), DIMENSION(1) :: X
      COMMON /FMCOM / X
        INTEGER :: IDAF, IJKO, IJKT, IP, IR, IW, NAV
        INTEGER, DIMENSION(950) :: IODA
      COMMON /IOFILE/ IR, IW, IP, IJKO, IJKT, IDAF, NAV, IODA
        LOGICAL :: DSKWRK, GOPARR, MASWRK
        INTEGER :: IBTYP, IPTIM, MASTER, ME, NPROC
      COMMON /PAR   / ME, MASTER, NPROC, IBTYP, IPTIM, GOPARR, DSKWRK,  &
     &                MASWRK
C
      INTEGER :: L1, L2, L3, L7, LX, NDSR, NOCA, NOCB
      REAL(KIND=dp), DIMENSION(NDSR,NDSR) :: DSTAS
      REAL(KIND=dp), DIMENSION(L7,NDSR) :: DXPYN, DXPYO
C
      REAL(KIND=dp), DIMENSION(NDSR,NOCA,LX-NOCB) :: ALPHAM, BETAM
      REAL(KIND=dp), DIMENSION(NDSR,NOCA,NOCA) :: DELTAM, GAMMAM
      LOGICAL :: DLB, SLB
      INTEGER :: I, IIN, IIO, IIP, IIQ, IIR, IIS, IJS, IJX, IKSMSA,     &
     &           IKSMSB, IOC, IOC1, IOC2, IPCOUNT, IPQ, IPR, IPS, IQS,  &
     &           IRQ, IRS, ISTSAB, ISTSIA, ISTSIJ, IVIR, J, JOC, JVIR,  &
     &           LAST, LOADFM, MINE, NEED, NEXT, NOC, NVIRB
      REAL(KIND=dp) :: STSACM
C
C     ----- INITIALIZATION FOR PARALLEL WORK -----
C     BOTH STATIC AND DYNAMIC LOAD BALANCING ARE IMPLEMENTED BELOW
C
      SLB = GOPARR  .AND.  IBTYP.EQ.0
      DLB = GOPARR  .AND.  IBTYP.EQ.1
      NEXT = -1
      MINE = -1
      IPCOUNT = ME - 1
C
C INITIALIZE
      DO I=1,NDSR
      DO J=1,NDSR
         DSTAS(I,J)=0.0D0
      ENDDO
      ENDDO
C
      NOC      = NOCA-1
      NVIRB    = LX-NOCB
C
C MEMORY ALLOCATION
      CALL VALFM(LOADFM)
      IKSMSA= LOADFM + 1
      IKSMSB= IKSMSA + LX*LX
      ISTSIJ= IKSMSB + LX*LX
      ISTSAB= ISTSIJ + NOCA*NOCA
      ISTSIA= ISTSAB + NVIRB*NVIRB
      LAST  = ISTSIA + NOCA*NVIRB
      NEED  = LAST - LOADFM-1
      CALL GETFM(NEED)
C
C     < MO_M (T-DT) | MO_N (T) >
      CALL OVKSMO(X(IKSMSA),X(IKSMSB),L1,L2,L3,LX)
C
C     S_ij, S_ab, S_ia
      CALL MRSFTLF(X(ISTSIJ),X(ISTSAB),X(ISTSIA),X(IKSMSA),X(IKSMSB),
     *             NOCA,NOCB,NVIRB,LX)
C
C CONTRACTION
      DO IIO=1,NDSR
         DO IOC=1,NOCA
         DO IVIR=1,NVIRB
            ALPHAM(IIO,IOC,IVIR)=0.D+00
            DO JVIR=1,NVIRB
               IF((IOC.GT.NOCB).AND.(JVIR.LE.2)) CYCLE
               IJX=(JVIR-1)*NOCA +IOC
               IJS=(IVIR-1)*NVIRB+JVIR
               ALPHAM(IIO,IOC,IVIR)=ALPHAM(IIO,IOC,IVIR)+
     *         DXPYO(IJX,IIO)*X(ISTSAB+IJS-1)
            ENDDO
         ENDDO
         ENDDO
      ENDDO
C
      DO IIN=1,NDSR
         DO IOC=1,NOCA
         DO IVIR=1,NVIRB
            BETAM(IIN,IOC,IVIR)=0.D+00
            DO JOC=1,NOCA
              IF((JOC.GT.NOCB).AND.(IVIR.LE.2)) CYCLE
              IJS=(JOC-1)*NOCA+IOC
              IJX=(IVIR-1)*NOCA+JOC
              BETAM(IIN,IOC,IVIR)=BETAM(IIN,IOC,IVIR)+
     *        X(ISTSIJ+IJS-1)*DXPYN(IJX,IIN)
            ENDDO
         ENDDO
         ENDDO
      ENDDO
C
      DO IIO=1,NDSR
         DO IOC1=1,NOCA
         DO IOC2=1,NOCA
            GAMMAM(IIO,IOC1,IOC2)=0.D+00
            DO JVIR=1,NVIRB
              IF((IOC1.GT.NOCB).AND.(JVIR.LE.2)) CYCLE
              IJS=(JVIR-1)*NOCA+IOC1
              IJX=(JVIR-1)*NOCA+IOC2
              GAMMAM(IIO,IOC1,IOC2)=GAMMAM(IIO,IOC1,IOC2)+
     *        DXPYO(IJS,IIO)*X(ISTSIA+IJX-1)
            ENDDO
         ENDDO
         ENDDO
      ENDDO
C
      DO IIN=1,NDSR
         DO IOC1=1,NOCA
         DO IOC2=1,NOCA
            DELTAM(IIN,IOC1,IOC2)=0.D+00
            DO JVIR=1,NVIRB
              IF((IOC2.GT.NOCB).AND.(JVIR.LE.2)) CYCLE
              IJS=(JVIR-1)*NOCA+IOC1
              IJX=(JVIR-1)*NOCA+IOC2
              DELTAM(IIN,IOC1,IOC2)=DELTAM(IIN,IOC1,IOC2)+
     *        X(ISTSIA+IJS-1)*DXPYN(IJX,IIN)
            ENDDO
         ENDDO
         ENDDO
      ENDDO
C
C 4INDEX SUMMATION
      STSACM = 0.0D+00
      DO IIO=1,NDSR
      DO IIN=1,NDSR
C
C        IF ( IIO .EQ. IIN ) CYCLE
C
         STSACM=0.D+00
         DO IIP=NOCB+1,NOCA
         DO IIQ=NOCB+1,NOCA
            IPQ=(IIQ-NOCB-1)*NOCA+IIP
            DO IIR=NOCB+1,NOCA
               IPR=(IIR-1)*NOCA+IIP
            DO IIS=NOCB+1,NOCA
               IRS=(IIS-NOCB-1)*NOCA+IIR
               IQS=(IIS-NOCB-1)*NVIRB+IIQ-NOCB
C
               STSACM = STSACM +
     *         DXPYO(IPQ,IIO)*X(ISTSIJ+IPR-1)
     *        *X(ISTSAB+IQS-1)*DXPYN(IRS,IIN)
            ENDDO
            ENDDO
         ENDDO
         ENDDO
C
         DO IIP=NOCB+1,NOCA
         DO IIQ=NOCB+1,NOCA
            IPQ=(IIQ-NOCB-1)*NOCA+IIP
            DO IIR=1,NOCA
               IPR=(IIR-1)*NOCA+IIP
               IRQ=(IIQ-NOCB-1)*NOCA+IIR
            DO IIS=NOCB+1,LX
               IF((IIR.GE.NOCB+1).AND.(IIS.LE.NOCA)) CYCLE
               IPS=(IIS-NOCB-1)*NOCA+IIP
               IRS=(IIS-NOCB-1)*NOCA+IIR
               IQS=(IIS-NOCB-1)*NVIRB+IIQ-NOCB
C
               STSACM = STSACM +
     *         DXPYO(IPQ,IIO)*DXPYN(IRS,IIN)
     *         *(X(ISTSIJ+IPR-1)*X(ISTSAB+IQS-1)
     *          +X(ISTSIA+IPS-1)*X(ISTSIA+IRQ-1))/sqrt(two)
            ENDDO
            ENDDO
         ENDDO
         ENDDO
C
         DO IIP=1,NOCA
         DO IIQ=NOCB+1,LX
            IF((IIP.GE.NOCB+1).AND.(IIQ.LE.NOCA)) CYCLE
C
            IPQ=(IIQ-NOCB-1)*NOCA+IIP
            DO IIR=NOCB+1,NOCA
               IPR=(IIR-1)*NOCA+IIP
               IRQ=(IIQ-NOCB-1)*NOCA+IIR
            DO IIS=NOCB+1,NOCA
               IPS=(IIS-NOCB-1)*NOCA+IIP
               IRS=(IIS-NOCB-1)*NOCA+IIR
               IQS=(IIS-NOCB-1)*NVIRB+IIQ-NOCB
C
               STSACM = STSACM +
     *         DXPYO(IPQ,IIO)*DXPYN(IRS,IIN)
     *         *(X(ISTSIJ+IPR-1)*X(ISTSAB+IQS-1)
     *          +X(ISTSIA+IPS-1)*X(ISTSIA+IRQ-1))/sqrt(two)
            ENDDO
            ENDDO
         ENDDO
         ENDDO
C
         DO IOC =1,NOCA
         DO IVIR=1,NVIRB
            STSACM = STSACM +
     *      ALPHAM(IIO,IOC,IVIR)*BETAM(IIN,IOC,IVIR)
         ENDDO
         ENDDO
C
         DO IOC1=1,NOCA
         DO IOC2=1,NOCA
            STSACM = STSACM +
     *      GAMMAM(IIO,IOC1,IOC2)*DELTAM(IIN,IOC1,IOC2)
         ENDDO
         ENDDO
         DSTAS(IIO,IIN)=STSACM
C
      ENDDO
      ENDDO
C
C MEMORY DEALLOCATION
      CALL RETFM(NEED)
      END SUBROUTINE MRSFOV
C*MODULE SFDFT   *DECK MRSFTLF
C>
C>     @brief Compute overlap integrals between CSFs
C>            of MRSF-TDDFT using TLF approximation
C>
C>     @details TLF approximation for SF- and LR-TDDFT
C>              is introduced in JCTC 15 882 (2019)
C>
C>     @author Seunghoon Lee
C>
      SUBROUTINE MRSFTLF(DSTSIJ,DSTSAB,DSTSIA,DKSMSA,DKSMSB,NOCA,NOCB,  &
     &                   NVIRB,LX)
      USE comm_NONAD, ONLY: NAMD,NDSWCH,NDRST,NDTLF,COLD,THRSHE
      USE mx_limits, ONLY: mxatm
      USE prec, ONLY: dp
      IMPLICIT NONE
C
        INTEGER :: IDAF, IJKO, IJKT, IP, IR, IW, NAV
        INTEGER, DIMENSION(950) :: IODA
      COMMON /IOFILE/ IR, IW, IP, IJKO, IJKT, IDAF, NAV, IODA
        LOGICAL :: DSKWRK, GOPARR, MASWRK
        INTEGER :: IBTYP, IPTIM, MASTER, ME, NPROC
      COMMON /PAR   / ME, MASTER, NPROC, IBTYP, IPTIM, GOPARR, DSKWRK,  &
     &                MASWRK
C
      INTEGER :: LX, NOCA, NOCB, NVIRB
      REAL(KIND=dp), DIMENSION(LX,LX) :: DKSMSA, DKSMSB
      REAL(KIND=dp), DIMENSION(NVIRB,NVIRB) :: DSTSAB
      REAL(KIND=dp), DIMENSION(NOCA,NVIRB) :: DSTSIA
      REAL(KIND=dp), DIMENSION(NOCA,NOCA) :: DSTSIJ
C
      LOGICAL :: DLB, NEWRND, SLB
      INTEGER :: I, I1, I2, IA1, IA2, IPCOUNT, J1, J2, MINE, NEXT, NOC
      REAL(KIND=dp) :: PRECOMP, TEMP1, TEMP2, TMP, TMP1, TMP2
C
C     ----- INITIALIZATION FOR PARALLEL WORK -----
C     BOTH STATIC AND DYNAMIC LOAD BALANCING ARE IMPLEMENTED BELOW
C
      SLB = GOPARR  .AND.  IBTYP.EQ.0
      DLB = GOPARR  .AND.  IBTYP.EQ.1
      NEXT = -1
      MINE = -1
      IPCOUNT = ME - 1
      NOC=NOCA-1

c     ndtlf=0 : no tlf
c     ndtlf=1 : tlf(1)
c     ndtlf=2 : tlf(2)

      CALL VCLR(DSTSIJ,1,NOCA*NOCA)
      CALL VCLR(DSTSAB,1,NVIRB*NVIRB)
      CALL VCLR(DSTSIA,1,NOCA*NVIRB)

      IF(NDTLF.EQ.1) THEN

         PRECOMP=1.D+00
         DO I=1,NOCA
            PRECOMP=PRECOMP*DKSMSA(I,I)
         ENDDO

C        ALPHA dETERMINANT
         DO I1 = 1, NOCA
            DO I2 = 1, NOCA
C
C           ----- GO PARALLEL! -----
C
            IF (DLB) THEN
               MINE = MINE + 1
               IF (MINE.GT.NEXT) CALL DDI_DLBNEXT(NEXT)
               IF (NEXT.NE.MINE) CYCLE
            END IF
            IF (SLB) THEN
               IPCOUNT = IPCOUNT + 1
               IF (MOD(IPCOUNT,NPROC).NE.0) CYCLE
            END IF

            CALL TLFEXP(TMP,11,I1,I2,DKSMSA,PRECOMP,NOCA,LX)
            IF ( I1. NE. I2 ) THEN
               TMP= -1.D+00 * TMP
            ENDIF

            DSTSIJ(I1,I2)=TMP

            END DO
         END DO

         IF(DLB) CALL DDI_DLBRESET
         IF(GOPARR) CALL DDI_GSUMF(2361,DSTSIJ,NOCA*NOCA)
         NEXT = -1
         MINE = -1
         IPCOUNT = ME - 1

      ELSE IF(NDTLF .EQ. 2) THEN

         PRECOMP=1.D+00
         DO I=1,NOCA
            PRECOMP=PRECOMP*DKSMSA(I,I)
         ENDDO

C        alpha Determinant
         DO I1 = 1, NOCA
            DO I2 = 1, NOCA
C
C           ----- GO PARALLEL! -----
C
            IF (DLB) THEN
               MINE = MINE + 1
               IF (MINE.GT.NEXT) CALL DDI_DLBNEXT(NEXT)
               IF (NEXT.NE.MINE) CYCLE
            END IF
            IF (SLB) THEN
               IPCOUNT = IPCOUNT + 1
               IF (MOD(IPCOUNT,NPROC).NE.0) CYCLE
            END IF

            CALL TLFEXP(TMP1,11,I1,I2,DKSMSA,PRECOMP,NOCA,LX)
            CALL TLFEXP(TMP2,12,I1,I2,DKSMSA,PRECOMP,NOCA,LX)
            TMP = TMP1+TMP2
            IF ( I1. NE. I2 ) THEN
               TMP= -1.D+00 * TMP
            ENDIF

            DSTSIJ(I1,I2)=TMP

            END DO
         END DO

         IF(DLB) CALL DDI_DLBRESET
         IF(GOPARR) CALL DDI_GSUMF(2361,DSTSIJ,NOCA*NOCA)
         NEXT = -1
         MINE = -1
         IPCOUNT = ME - 1

      ELSE

C        alpha Determinant
         DO I1 = 1, NOCA
            DO I2 = 1, NOCA
C
C           ----- GO PARALLEL! -----
C
            IF (DLB) THEN
               MINE = MINE + 1
               IF (MINE.GT.NEXT) CALL DDI_DLBNEXT(NEXT)
               IF (NEXT.NE.MINE) CYCLE
            END IF
            IF (SLB) THEN
               IPCOUNT = IPCOUNT + 1
               IF (MOD(IPCOUNT,NPROC).NE.0) CYCLE
            END IF

            CALL OVEXACT(TEMP1,I1,I2,IA1,IA2,DKSMSA,LX,1,NOC,1)

            DSTSIJ(I1,I2)=TEMP1

            END DO
         END DO

         IF(DLB) CALL DDI_DLBRESET
         IF(GOPARR) CALL DDI_GSUMF(2361,DSTSIJ,NOCA*NOCA)
         NEXT = -1
         MINE = -1
         IPCOUNT = ME - 1

      ENDIF

      IF(NDTLF .EQ. 1) THEN

        PRECOMP=1.D+00
        DO I=1,NOCA-2
           PRECOMP=PRECOMP*DKSMSB(I,I)
        ENDDO

C       1-2 DET
         DO J1 = 1, NVIRB
            IA1 = NOCB + J1
            DO J2 = 1, NVIRB
               IA2 = NOCB + J2
C
C                ----- GO PARALLEL! -----
C
               IF (DLB) THEN
                  MINE = MINE + 1
                  IF (MINE.GT.NEXT) CALL DDI_DLBNEXT(NEXT)
                  IF (NEXT.NE.MINE) CYCLE
               END IF
               IF (SLB) THEN
                  IPCOUNT = IPCOUNT + 1
                  IF (MOD(IPCOUNT,NPROC).NE.0) CYCLE
               END IF

               CALL TLFEXP(TMP,21,IA1,IA2,DKSMSB,PRECOMP,NOCA,LX)

               DSTSAB(J1,J2)=TMP

            END DO
         END DO

        IF(DLB) CALL DDI_DLBRESET
        IF(GOPARR) CALL DDI_GSUMF(2362,DSTSAB,NVIRB*NVIRB)
        NEXT = -1
        MINE = -1
        IPCOUNT = ME - 1

      ELSE IF(NDTLF .EQ. 2) THEN

        PRECOMP=1.D+00
        DO I=1,NOCA-2
           PRECOMP=PRECOMP*DKSMSB(I,I)
        ENDDO

C       1-2 DET
         DO J1 = 1, NVIRB
            IA1 = NOCB + J1
            DO J2 = 1, NVIRB
               IA2 = NOCB + J2
C
C          ----- GO PARALLEL! -----
C
               IF (DLB) THEN
                  MINE = MINE + 1
                  IF (MINE.GT.NEXT) CALL DDI_DLBNEXT(NEXT)
                  IF (NEXT.NE.MINE) CYCLE
               END IF
               IF (SLB) THEN
                  IPCOUNT = IPCOUNT + 1
                  IF (MOD(IPCOUNT,NPROC).NE.0) CYCLE
               END IF

               CALL TLFEXP(TMP1,21,IA1,IA2,DKSMSB,PRECOMP,NOCA,LX)
               CALL TLFEXP(TMP2,22,IA1,IA2,DKSMSB,PRECOMP,NOCA,LX)

               DSTSAB(J1,J2)=TMP1+TMP2

            END DO
         END DO

        IF(DLB) CALL DDI_DLBRESET
        IF(GOPARR) CALL DDI_GSUMF(2362,DSTSAB,NVIRB*NVIRB)
        NEXT = -1
        MINE = -1
        IPCOUNT = ME - 1

      ELSE

C        1-2 DET
         DO J1 = 1, NVIRB
            IA1 = NOCB + J1
            DO J2 = 1, NVIRB
               IA2 = NOCB + J2
C
C        ----- GO PARALLEL! -----
C
               IF (DLB) THEN
                  MINE = MINE + 1
                  IF (MINE.GT.NEXT) CALL DDI_DLBNEXT(NEXT)
                  IF (NEXT.NE.MINE) CYCLE
               END IF
               IF (SLB) THEN
                  IPCOUNT = IPCOUNT + 1
                  IF (MOD(IPCOUNT,NPROC).NE.0) CYCLE
               END IF

               CALL OVEXACT(TEMP2,I1,I2,IA1,IA2,DKSMSB,LX,1,NOC,2)

               DSTSAB(J1,J2)=TEMP2

            END DO
         END DO

         IF(DLB) CALL DDI_DLBRESET
         IF(GOPARR) CALL DDI_GSUMF(2362,DSTSAB,NVIRB*NVIRB)
         NEXT = -1
         MINE = -1
         IPCOUNT = ME - 1

      ENDIF

C        ALPHA DETERMINANT
      DO I1 = 1, NOCA
         DO J1 = 1, NVIRB
            IA1 = NOCB + J1
C
C           ----- GO PARALLEL! -----
C
            IF (DLB) THEN
               MINE = MINE + 1
               IF (MINE.GT.NEXT) CALL DDI_DLBNEXT(NEXT)
               IF (NEXT.NE.MINE) CYCLE
            END IF
            IF (SLB) THEN
               IPCOUNT = IPCOUNT + 1
               IF (MOD(IPCOUNT,NPROC).NE.0) CYCLE
            END IF

            CALL OVEXACT(TEMP1,I1,I2,IA1,IA2,DKSMSA,LX,1,NOC,3)

            DSTSIA(I1,J1)=TEMP1

         END DO
      END DO

      IF(DLB) CALL DDI_DLBRESET
      IF(GOPARR) CALL DDI_GSUMF(2363,DSTSIA,NOCA*NVIRB)
      NEXT = -1
      MINE = -1
      IPCOUNT = ME - 1

      END SUBROUTINE MRSFTLF
C*MODULE SFDFT   *DECK TLFEXP
C>
C>     @brief Compute each order of expansion terms
C>            in Leibniz formula for MRSF-TDDFT
C>
C>     @author Seunghoon Lee
C>
      SUBROUTINE TLFEXP(DPTOV,ITYP,I1,I2,DKSMS,PRECOMP,NOCASUB,L1)
C
      USE prec, ONLY: dp
      IMPLICIT NONE
C
        INTEGER :: IDAF, IJKO, IJKT, IP, IR, IW, NAV
        INTEGER, DIMENSION(950) :: IODA
      COMMON /IOFILE/ IR, IW, IP, IJKO, IJKT, IDAF, NAV, IODA
        LOGICAL :: DSKWRK, GOPARR, MASWRK
        INTEGER :: IBTYP, IPTIM, MASTER, ME, NPROC
      COMMON /PAR   / ME, MASTER, NPROC, IBTYP, IPTIM, GOPARR, DSKWRK,  &
     &                MASWRK
C
      REAL(KIND=dp) :: DPTOV, PRECOMP
      INTEGER :: I1, I2, ITYP, L1, NOCASUB
      REAL(KIND=dp), DIMENSION(L1,L1) :: DKSMS
C
      REAL(KIND=dp) :: DPTOV1, DPTOV2
      INTEGER :: IA1, IA2, L, LP

C     ityp=11 : Alpha 1st order
C     ityp=21 : Beta  1st order
C     ityp=12 : Alpha 2nd order
C     ityp=22 : Beta  2nd order

      IF ( ITYP .EQ. 11 ) THEN

         DPTOV = PRECOMP*DKSMS(I2,I1)/(DKSMS(I1,I1)*DKSMS(I2,I2))
         RETURN

      ELSE IF ( ITYP .EQ. 21 ) THEN

         IA1=I1
         IA2=I2
         DPTOV=PRECOMP*DKSMS(IA1,IA2)
         RETURN

      ELSE IF ( ITYP .EQ. 12 ) THEN

         IF ( I1 .NE. I2 ) THEN

            DPTOV=0.D+00

            DO L=1,NOCASUB
              IF (L.NE.I1 .AND. L.NE.I2) THEN
                 DPTOV = DPTOV + DKSMS(I2,L)*DKSMS(L,I1)/DKSMS(L,L)
              ENDIF
            ENDDO
            DPTOV =
     *      (-1.D+0*PRECOMP)*DPTOV/(DKSMS(I1,I1)*DKSMS(I2,I2))

            RETURN

         ELSE

            DPTOV=0.D+00
            DO L=1,NOCASUB-1
            IF (L.NE.I1) THEN
            DO LP=L+1,NOCASUB
            IF (LP.NE.I1) THEN
               DPTOV = DPTOV +
     *         DKSMS(L,LP)*DKSMS(LP,L)/(DKSMS(L,L)*DKSMS(LP,LP))
            ENDIF
            ENDDO
            ENDIF
            ENDDO

            DPTOV = DPTOV * (-1.D+0*PRECOMP)/DKSMS(I1,I1)
            RETURN

         ENDIF

      ELSE IF ( ITYP .EQ. 22 ) THEN

         IA1=I1
         IA2=I2

         IF ( IA1 .NE. IA2 ) THEN

            DPTOV=0.D+00
            DO L=1,NOCASUB-2
               DPTOV = DPTOV +
     *         DKSMS(IA1,L)*DKSMS(L,IA2)/DKSMS(L,L)
            ENDDO
            DPTOV = DPTOV * (-1.D+0*PRECOMP)
            RETURN

         ELSE

            DPTOV1=0.D+00
            DO L=1,NOCASUB-2
               DPTOV1 = DPTOV1 +
     *         DKSMS(IA1,L)*DKSMS(L,IA2)/DKSMS(L,L)
            ENDDO
            DPTOV2=0.D+00
            DO L=1,NOCASUB-3
            DO LP=L+1,NOCASUB-2
               DPTOV2 = DPTOV2 +
     *         DKSMS(L,LP)*DKSMS(LP,L)/(DKSMS(L,L)*DKSMS(LP,LP))
            ENDDO
            ENDDO
            DPTOV2 = DPTOV2 * DKSMS(IA1,IA1)

            DPTOV = (-1.D+0*PRECOMP)*(DPTOV1+DPTOV2)
            RETURN

         ENDIF

      ENDIF

      END SUBROUTINE TLFEXP
C*MODULE SFDFT   *DECK OVEXACT
      SUBROUTINE OVEXACT(TEMP1,I1,I2,IA1,IA2,DKSMS,L1,ILOW,NOCSUB,ITYP)
C
      USE prec, ONLY: dp
      IMPLICIT NONE
C
        INTEGER :: IDAF, IJKO, IJKT, IP, IR, IW, NAV
        INTEGER, DIMENSION(950) :: IODA
      COMMON /IOFILE/ IR, IW, IP, IJKO, IJKT, IDAF, NAV, IODA
        LOGICAL :: DSKWRK, GOPARR, MASWRK
        INTEGER :: IBTYP, IPTIM, MASTER, ME, NPROC
      COMMON /PAR   / ME, MASTER, NPROC, IBTYP, IPTIM, GOPARR, DSKWRK,  &
     &                MASWRK
C
      INTEGER :: I1, I2, IA1, IA2, ILOW, ITYP, L1, NOCSUB
      REAL(KIND=dp) :: TEMP1
      REAL(KIND=dp), DIMENSION(L1,L1) :: DKSMS
C
      REAL(KIND=dp), DIMENSION(NOCSUB*NOCSUB) :: DDET
      INTEGER :: I, IIPP, IMAX, IMIN, IPP
C
      IF (ITYP.EQ.1) THEN
         IMIN=MIN(I1,I2)
         IMAX=MAX(I1,I2)
C     (1,1) BLOCK
         DO I = 1, IMIN-1
         DO IPP = 1, IMIN-1
           IIPP = (IPP-1)*NOCSUB+I
           DDET(IIPP) = DKSMS(I+ILOW-1,IPP+ILOW-1)
         END DO
         END DO
C     (1,2) BLOCK
         DO I = 1, IMIN-1
         DO IPP = IMIN, IMAX-2
           IIPP = (IPP-1)*NOCSUB+I
           DDET(IIPP) = DKSMS(I+ILOW-1,IPP+1+ILOW-1)
         END DO
         END DO
C     (1,3) BLOCK
         DO I = 1, IMIN-1
         DO IPP = IMAX-1, NOCSUB-1
           IIPP = (IPP-1)*NOCSUB+I
           DDET(IIPP) = DKSMS(I+ILOW-1,IPP+2+ILOW-1)
         END DO
         END DO
C     (2,1) BLOCK
         DO I = IMIN, IMAX-2
         DO IPP = 1, IMIN-1
            IIPP = (IPP-1)*NOCSUB+I
            DDET(IIPP) = DKSMS(I+1+ILOW-1,IPP+ILOW-1)
         END DO
         END DO
C     (2,2) BLOCK
         DO I = IMIN, IMAX-2
         DO IPP = IMIN, IMAX-2
           IIPP = (IPP-1)*NOCSUB+I
           DDET(IIPP) = DKSMS(I+1+ILOW-1,IPP+1+ILOW-1)
         END DO
         END DO
C     (2,3) BLOCK
         DO I = IMIN, IMAX-2
         DO IPP = IMAX-1, NOCSUB-1
           IIPP = (IPP-1)*NOCSUB+I
           DDET(IIPP) = DKSMS(I+1+ILOW-1,IPP+2+ILOW-1)
         END DO
         END DO
C     (3,1) BLOCK
         DO I = IMAX-1, NOCSUB-1
         DO IPP = 1, IMIN-1
            IIPP = (IPP-1)*NOCSUB+I
            DDET(IIPP) = DKSMS(I+2+ILOW-1,IPP+ILOW-1)
         END DO
         END DO
C     (3,2) BLOCK
         DO I = IMAX-1, NOCSUB-1
         DO IPP = IMIN, IMAX-2
           IIPP = (IPP-1)*NOCSUB+I
           DDET(IIPP) = DKSMS(I+2+ILOW-1,IPP+1+ILOW-1)
         END DO
         END DO
C     (3,3) BLOCK
         DO I = IMAX-1, NOCSUB-1
         DO IPP = IMAX-1, NOCSUB-1
           IIPP = (IPP-1)*NOCSUB+I
           DDET(IIPP) = DKSMS(I+2+ILOW-1,IPP+2+ILOW-1)
         END DO
         END DO

C     (1,4) BLOCK
         DO I = 1, IMIN-1
         DO IPP = NOCSUB, NOCSUB
           IIPP = (IPP-1)*NOCSUB+I
           DDET(IIPP) = DKSMS(I+ILOW-1,I1+ILOW-1)
         END DO
         END DO
C     (2,4) BLOCK
         DO I = IMIN, IMAX-2
         DO IPP = NOCSUB, NOCSUB
           IIPP = (IPP-1)*NOCSUB+I
           DDET(IIPP) = DKSMS(I+1+ILOW-1,I1+ILOW-1)
         END DO
         END DO
C     (3,4) BLOCK
         DO I = IMAX-1, NOCSUB-1
         DO IPP = NOCSUB, NOCSUB
           IIPP = (IPP-1)*NOCSUB+I
           DDET(IIPP) = DKSMS(I+2+ILOW-1,I1+ILOW-1)
         END DO
         END DO

C     (4,1) BLOCK
         DO I = NOCSUB, NOCSUB
         DO IPP = 1, IMIN-1
           IIPP = (IPP-1)*NOCSUB+I
           DDET(IIPP) = DKSMS(I2+ILOW-1,IPP+ILOW-1)
         END DO
         END DO
C     (4,2) BLOCK
         DO I = NOCSUB, NOCSUB
         DO IPP = IMIN, IMAX-2
           IIPP = (IPP-1)*NOCSUB+I
           DDET(IIPP) = DKSMS(I2+ILOW-1,IPP+1+ILOW-1)
         END DO
         END DO
C     (4,3) BLOCK
         DO I = NOCSUB, NOCSUB
         DO IPP = IMAX-1, NOCSUB-1
           IIPP = (IPP-1)*NOCSUB+I
           DDET(IIPP) = DKSMS(I2+ILOW-1,IPP+2+ILOW-1)
         END DO
         END DO

C     (4,4) BLOCK
         DO I = NOCSUB, NOCSUB
         DO IPP = NOCSUB, NOCSUB
           IIPP = (IPP-1)*NOCSUB+I
           DDET(IIPP) = DKSMS(I2+ILOW-1,I1+ILOW-1)
         END DO
         END DO

C     CALCULATE ALPHA dETERMINANT
         CALL POLDET(TEMP1,DDET,NOCSUB,NOCSUB)

         IF ( I1. EQ. I2 ) THEN
            RETURN
         ELSE IF ( I1. NE. I2 ) THEN
            TEMP1 = -1.D+00 * TEMP1
            RETURN
         ENDIF

      ELSEIF (ITYP.EQ.2) THEN

C     (1,1) BLOCK
         DO I = 1, NOCSUB-1
         DO IPP = 1, NOCSUB-1
            IIPP = (IPP-1)*NOCSUB+I
            DDET(IIPP) = DKSMS(I+ILOW-1,IPP+ILOW-1)
         END DO
         END DO
C     (1,2) BLOCK
         IPP = NOCSUB
         DO I = 1, NOCSUB-1
            IIPP = (IPP-1)*NOCSUB+I
            DDET(IIPP) = DKSMS(I+ILOW-1,IA2)
         END DO
C     (2,1) BLOCK
         I = NOCSUB
         DO IPP = 1, NOCSUB-1
            IIPP = (IPP-1)*NOCSUB+I
            DDET(IIPP) = DKSMS(IA1,IPP+ILOW-1)
         END DO
C     (2,2) BLOCK
         I = NOCSUB
         IPP = NOCSUB
         IIPP = (IPP-1)*NOCSUB+I
         DDET(IIPP) = DKSMS(IA1,IA2)
C     CALCULATE 2 DET
         CALL POLDET(TEMP1,DDET,NOCSUB,NOCSUB)

         RETURN

      ELSE IF (ITYP.EQ.3) THEN
C     (1,1) BLOCK
         DO I = 1, I1-1
         DO IPP = 1, I1-1
           IIPP = (IPP-1)*NOCSUB+I
           DDET(IIPP) = DKSMS(I,IPP)
         END DO
         END DO
C     (1,2) BLOCK
         DO I = 1, I1-1
         DO IPP = I1, NOCSUB-2
           IIPP = (IPP-1)*NOCSUB+I
           DDET(IIPP) = DKSMS(I,IPP+1)
         END DO
         END DO
C     (1,3) BLOCK
         DO I = 1, I1-1
           IPP = NOCSUB-1
           IIPP = (IPP-1)*NOCSUB+I
           DDET(IIPP) = DKSMS(I,I1)
         END DO
C     (1,4) BLOCK
         DO I = 1, I1-1
           IPP = NOCSUB
           IIPP = (IPP-1)*NOCSUB+I
           DDET(IIPP) = DKSMS(I,IA1)
         END DO
C     (2,1) BLOCK
         DO I = I1, NOCSUB-2
         DO IPP = 1, I1-1
            IIPP = (IPP-1)*NOCSUB+I
            DDET(IIPP) = DKSMS(I+1,IPP)
         END DO
         END DO
C     (2,2) BLOCK
         DO I = I1, NOCSUB-2
         DO IPP = I1, NOCSUB-2
           IIPP = (IPP-1)*NOCSUB+I
           DDET(IIPP) = DKSMS(I+1,IPP+1)
         END DO
         END DO
C     (2,3) BLOCK
         DO I = I1, NOCSUB-2
           IPP = NOCSUB-1
           IIPP = (IPP-1)*NOCSUB+I
           DDET(IIPP) = DKSMS(I+1,I1)
         END DO
C     (2,4) BLOCK
         DO I = I1, NOCSUB-2
           IPP = NOCSUB
           IIPP = (IPP-1)*NOCSUB+I
           DDET(IIPP) = DKSMS(I+1,IA1)
         END DO
C     (3,1) BLOCK
         I = NOCSUB-1
         DO IPP = 1, I1-1
           IIPP = (IPP-1)*NOCSUB+I
           DDET(IIPP) = DKSMS(I+1,IPP)
         END DO
C     (3,2) BLOCK
         I = NOCSUB-1
         DO IPP = I1, NOCSUB-2
           IIPP = (IPP-1)*NOCSUB+I
           DDET(IIPP) = DKSMS(I+1,IPP+1)
         END DO
C     (3,3) BLOCK
         I   = NOCSUB-1
         IPP = NOCSUB-1
         IIPP = (IPP-1)*NOCSUB+I
         DDET(IIPP) = DKSMS(I+1,I1)
C     (3,4) BLOCK
         I   = NOCSUB-1
         IPP = NOCSUB
         IIPP = (IPP-1)*NOCSUB+I
         DDET(IIPP) = DKSMS(I+1,IA1)
C     (4,1) BLOCK
         I = NOCSUB
         DO IPP = 1, I1-1
            IIPP = (IPP-1)*NOCSUB+I
            DDET(IIPP) = DKSMS(I+1,IPP)
         END DO
C     (4,2) BLOCK
         I = NOCSUB
         DO IPP = I1, NOCSUB-2
           IIPP = (IPP-1)*NOCSUB+I
           DDET(IIPP) = DKSMS(I+1,IPP+1)
         END DO
C     (4,3) BLOCK
         I   = NOCSUB
         IPP = NOCSUB-1
         IIPP = (IPP-1)*NOCSUB+I
         DDET(IIPP) = DKSMS(I+1,I1)
C     (4,4) BLOCK
         I   = NOCSUB
         IPP = NOCSUB
         IIPP = (IPP-1)*NOCSUB+I
         DDET(IIPP) = DKSMS(I+1,IA1)

C     CALCULATE ALPHA dETERMINANT
         CALL POLDET(TEMP1,DDET,NOCSUB,NOCSUB)

      ENDIF
      END SUBROUTINE OVEXACT
C*MODULE NAMD   *DECK NACVFD
C>
C>     @brief Compute nonadiabatic coupling vectors
C>            using the finite difference method
C>
C>     @author Seunghoon Lee
C>
      SUBROUTINE NACVFD(DNCME,DSTAS,DTAU,NTST)
C
      USE prec, ONLY: dp
      IMPLICIT NONE
C
      REAL(KIND=dp) :: DTAU
      INTEGER :: NTST
      REAL(KIND=dp), DIMENSION(NTST,NTST) :: DNCME, DSTAS
C
      INTEGER :: I, J
C
      DO I = 1, NTST
         DO J = 1, NTST
            DNCME(I,J) = (DSTAS(I,J)-DSTAS(J,I))/(2*DTAU)
         END DO
      END DO
C
      END SUBROUTINE NACVFD
C*MODULE NAMD   *DECK PPTDECOE
C>
C>     @brief Propagate time-dependent expansion coefficients
C>            using Runge-Kutta 4th order method (RK4)
C>
C>     @author Seunghoon Lee
C>
      SUBROUTINE PPTDECOE(DTCOE,DNCME,DSTEG,NTST,DTAUNINT)
C
      USE prec, ONLY: dp
      IMPLICIT NONE
C
        LOGICAL :: ALPHKWD, BETAKWD, MRDEA, MREKT, SG1T, TAMMD, TPA,    &
     &             TRIPLET
        REAL(KIND=dp) :: CNVTOL
        INTEGER :: IRECTD, ITDFG, ITDPRP, JANST, MAXVEC, MODTD, MTHST,  &
     &             MULTD, NLEBT, NONEQR, NPHIT, NRADT, NSTAT, NTHET,    &
     &             NTHST, NTRIAL
        INTEGER, DIMENSION(4) :: IFEDAT
        REAL(KIND=dp), DIMENSION(2) :: PFREQ
        REAL(KIND=dp), DIMENSION(3) :: SPCP
      COMMON /INFOTD/ CNVTOL, PFREQ, MODTD, JANST, NRADT, NTHET, NPHIT, &
     &                NLEBT, NSTAT, NTRIAL, MAXVEC, NTHST, IRECTD,      &
     &                ITDFG, ITDPRP, TRIPLET, SG1T, NONEQR, TAMMD, TPA, &
     &                ALPHKWD, BETAKWD, SPCP, MULTD, MREKT, MRDEA,      &
     &                MTHST, IFEDAT
C
      REAL(KIND=dp) :: DTAUNINT
      INTEGER :: NTST
      REAL(KIND=dp), DIMENSION(NTST,NTST) :: DNCME
      REAL(KIND=dp), DIMENSION(NTST) :: DSTEG
      REAL(KIND=dp), DIMENSION(2,NSTAT+1) :: DTCOE
C
      REAL(KIND=dp) :: DNORM
      REAL(KIND=dp), DIMENSION(2,NTST) :: DTMP1, DTMP2, DTMP3, DTMP4
      REAL(KIND=dp), DIMENSION(2,NSTAT+1) :: DTMPK1, DTMPK2, DTMPK3
      INTEGER :: I, J
C
C  RUNGE-KUTTA METHOD : FORT77 NUMERICAL RECIPE 16.1.3
C  K1 = DTMP1
      DO I=1,NTST
         CALL NDDTCR(DTMP1(1,I),I,DTCOE,DNCME,DSTEG,NTST,DTAUNINT)
         CALL NDDTCC(DTMP1(2,I),I,DTCOE,DNCME,DSTEG,NTST,DTAUNINT)
      ENDDO
C  Yn+K1/2 = DTMPK1
      DO I=1,2
      DO J=1,NTST
         DTMPK1(I,J)=DTCOE(I,J)+DTMP1(I,J)/2
      ENDDO
      ENDDO
C
C  K2 = DTMP2
      DO I=1,NTST
         CALL NDDTCR(DTMP2(1,I),I,DTMPK1,DNCME,DSTEG,NTST,DTAUNINT)
         CALL NDDTCC(DTMP2(2,I),I,DTMPK1,DNCME,DSTEG,NTST,DTAUNINT)
      ENDDO
C  Yn+K2/2 = DTMPK2
      DO I=1,2
      DO J=1,NTST
         DTMPK2(I,J)=DTCOE(I,J)+DTMP2(I,J)/2
      ENDDO
      ENDDO
C
C  K3 = DTMP3
      DO I=1,NTST
         CALL NDDTCR(DTMP3(1,I),I,DTMPK2,DNCME,DSTEG,NTST,DTAUNINT)
         CALL NDDTCC(DTMP3(2,I),I,DTMPK2,DNCME,DSTEG,NTST,DTAUNINT)
      ENDDO
C  Yn+K3 = DTMPK3
      DO I=1,2
      DO J=1,NTST
         DTMPK3(I,J)=DTCOE(I,J)+DTMP3(I,J)
      ENDDO
      ENDDO
C
C  K4 = DTMP4
      DO I=1,NTST
         CALL NDDTCR(DTMP4(1,I),I,DTMPK3,DNCME,DSTEG,NTST,DTAUNINT)
         CALL NDDTCC(DTMP4(2,I),I,DTMPK3,DNCME,DSTEG,NTST,DTAUNINT)
      ENDDO
C
      DO I=1,2
      DO J=1,NTST
         DTCOE(I,J) = DTCOE(I,J)+(DTMP1(I,J)+2*DTMP2(I,J)+
     *                2*DTMP3(I,J)+DTMP4(I,J))/6
      ENDDO
      ENDDO
C
C     NORMALIZED PROCEDURE
      DNORM=0.0D+00
      DO I=1,2
      DO J=1,NTST
         DNORM=DNORM+DTCOE(I,J)*DTCOE(I,J)
      ENDDO
      ENDDO
      DNORM = SQRT(DNORM)
      DO I=1,2
      DO J=1,NTST
         DTCOE(I,J) = DTCOE(I,J)/DNORM
      ENDDO
      ENDDO
      END SUBROUTINE PPTDECOE
C*MODULE NAMD   *DECK NDDTCR
C>
C>     @brief Real part propagation
C>
C>     @author Seunghoon Lee
C>
      SUBROUTINE NDDTCR(DTCR,NDCOEN,DTMPK,DNCME,DSTEG,NTST,DTAUNINT)
C
      USE prec, ONLY: dp
      IMPLICIT NONE
C
        LOGICAL :: ALPHKWD, BETAKWD, MRDEA, MREKT, SG1T, TAMMD, TPA,    &
     &             TRIPLET
        REAL(KIND=dp) :: CNVTOL
        INTEGER :: IRECTD, ITDFG, ITDPRP, JANST, MAXVEC, MODTD, MTHST,  &
     &             MULTD, NLEBT, NONEQR, NPHIT, NRADT, NSTAT, NTHET,    &
     &             NTHST, NTRIAL
        INTEGER, DIMENSION(4) :: IFEDAT
        REAL(KIND=dp), DIMENSION(2) :: PFREQ
        REAL(KIND=dp), DIMENSION(3) :: SPCP
      COMMON /INFOTD/ CNVTOL, PFREQ, MODTD, JANST, NRADT, NTHET, NPHIT, &
     &                NLEBT, NSTAT, NTRIAL, MAXVEC, NTHST, IRECTD,      &
     &                ITDFG, ITDPRP, TRIPLET, SG1T, NONEQR, TAMMD, TPA, &
     &                ALPHKWD, BETAKWD, SPCP, MULTD, MREKT, MRDEA,      &
     &                MTHST, IFEDAT
C
      REAL(KIND=dp) :: DTAUNINT, DTCR
      INTEGER :: NDCOEN, NTST
      REAL(KIND=dp), DIMENSION(NTST,NTST) :: DNCME
      REAL(KIND=dp), DIMENSION(NTST) :: DSTEG
      REAL(KIND=dp), DIMENSION(2,NSTAT+1) :: DTMPK
C
      INTEGER :: I
C
      DTCR = 0.0D+00
C
      DO I=1, NTST
         DTCR = DTCR - DNCME(NDCOEN,I)*DTMPK(1,I)
      ENDDO
C
      DTCR = DTCR + DSTEG(NDCOEN)*DTMPK(2,NDCOEN)
C
      DTCR = DTCR * DTAUNINT
C
      END SUBROUTINE NDDTCR
C*MODULE NAMD   *DECK NDDTCC
C>
C>     @brief Imaginary part propagation
C>
C>     @author Seunghoon Lee
C>
      SUBROUTINE NDDTCC(DTCC,NDCOEN,DTMPK,DNCME,DSTEG,NTST,DTAUNINT)
C
      USE prec, ONLY: dp
      IMPLICIT NONE
C
        LOGICAL :: ALPHKWD, BETAKWD, MRDEA, MREKT, SG1T, TAMMD, TPA,    &
     &             TRIPLET
        REAL(KIND=dp) :: CNVTOL
        INTEGER :: IRECTD, ITDFG, ITDPRP, JANST, MAXVEC, MODTD, MTHST,  &
     &             MULTD, NLEBT, NONEQR, NPHIT, NRADT, NSTAT, NTHET,    &
     &             NTHST, NTRIAL
        INTEGER, DIMENSION(4) :: IFEDAT
        REAL(KIND=dp), DIMENSION(2) :: PFREQ
        REAL(KIND=dp), DIMENSION(3) :: SPCP
      COMMON /INFOTD/ CNVTOL, PFREQ, MODTD, JANST, NRADT, NTHET, NPHIT, &
     &                NLEBT, NSTAT, NTRIAL, MAXVEC, NTHST, IRECTD,      &
     &                ITDFG, ITDPRP, TRIPLET, SG1T, NONEQR, TAMMD, TPA, &
     &                ALPHKWD, BETAKWD, SPCP, MULTD, MREKT, MRDEA,      &
     &                MTHST, IFEDAT
C
      REAL(KIND=dp) :: DTAUNINT, DTCC
      INTEGER :: NDCOEN, NTST
      REAL(KIND=dp), DIMENSION(NTST,NTST) :: DNCME
      REAL(KIND=dp), DIMENSION(NTST) :: DSTEG
      REAL(KIND=dp), DIMENSION(2,NSTAT+1) :: DTMPK
C
      INTEGER :: I
C
      DTCC = 0.0D+00
C
      DO I=1, NTST
         DTCC = DTCC - DNCME(NDCOEN,I)*DTMPK(2,I)
      ENDDO
C
      DTCC = DTCC - DSTEG(NDCOEN)*DTMPK(1,NDCOEN)
C
      DTCC = DTCC * DTAUNINT
C
      END SUBROUTINE NDDTCC
C*MODULE NAMD   *DECK FSSHPRST
C>
C>     @brief Compute surface hopping probability in sub-time scale
C>
C>     @author Seunghoon Lee
C>
      SUBROUTINE FSSHPRST(DHPPR,ITIME,DTCOE,DNCME,NTST,DTAUNINT,NDNINT)
C
      USE prec, ONLY: dp
      IMPLICIT NONE
C
        LOGICAL :: ALPHKWD, BETAKWD, MRDEA, MREKT, SG1T, TAMMD, TPA,    &
     &             TRIPLET
        REAL(KIND=dp) :: CNVTOL
        INTEGER :: IRECTD, ITDFG, ITDPRP, JANST, MAXVEC, MODTD, MTHST,  &
     &             MULTD, NLEBT, NONEQR, NPHIT, NRADT, NSTAT, NTHET,    &
     &             NTHST, NTRIAL
        INTEGER, DIMENSION(4) :: IFEDAT
        REAL(KIND=dp), DIMENSION(2) :: PFREQ
        REAL(KIND=dp), DIMENSION(3) :: SPCP
      COMMON /INFOTD/ CNVTOL, PFREQ, MODTD, JANST, NRADT, NTHET, NPHIT, &
     &                NLEBT, NSTAT, NTRIAL, MAXVEC, NTHST, IRECTD,      &
     &                ITDFG, ITDPRP, TRIPLET, SG1T, NONEQR, TAMMD, TPA, &
     &                ALPHKWD, BETAKWD, SPCP, MULTD, MREKT, MRDEA,      &
     &                MTHST, IFEDAT
C
      REAL(KIND=dp) :: DTAUNINT
      INTEGER :: ITIME, NDNINT, NTST
      REAL(KIND=dp), DIMENSION(NTST,NTST,NDNINT) :: DHPPR
      REAL(KIND=dp), DIMENSION(NTST,NTST) :: DNCME
      REAL(KIND=dp), DIMENSION(2,NSTAT+1) :: DTCOE
C
      REAL(KIND=dp) :: DTMP, TMP
      INTEGER :: I, J
C
C     DHPPR(I,J,ITIME) denotes hopping probability from I state to J state
C     at (Nuclear_time_step + ITIME*Nuclear_time_step/NDNINT)
      DO I = 1, NTST
      DO J = 1, NTST
         DTMP=DTCOE(1,I)*DTCOE(1,I)+DTCOE(2,I)*DTCOE(2,I)
         TMP=
     *   2*(DNCME(I,J)*(DTCOE(1,I)*DTCOE(1,J)+DTCOE(2,I)*DTCOE(2,J)))*
     *   DTAUNINT/DTMP
         DHPPR(I,J,ITIME)=TMP
      ENDDO
      ENDDO
C
      END SUBROUTINE FSSHPRST
C*MODULE NAMD   *DECK FSSHPR
C>
C>     @brief Obtain hopping probabilities
C>
C>     @author Seunghoon Lee
C>
      SUBROUTINE FSSHPR(DCMHP,DHPPR,NTST,NDNINT)
C
      USE prec, ONLY: dp
      IMPLICIT NONE
C
      INTEGER :: NDNINT, NTST
      REAL(KIND=dp), DIMENSION(NTST,NTST) :: DCMHP
      REAL(KIND=dp), DIMENSION(NTST,NTST,NDNINT) :: DHPPR
C
      REAL(KIND=dp) :: DNORM
      INTEGER :: I, ITIME, J
C
      DO I=1,NTST
      DO J=1,NTST
        DCMHP(I,J) = 0.0D+00
      ENDDO
      ENDDO
C
      DO I=1,NTST
      DO J=1,NTST
      DO ITIME=1,NDNINT
        DCMHP(I,J)=DCMHP(I,J)+DHPPR(I,J,ITIME)
      ENDDO
      IF(DCMHP(I,J).LT.0.0D+00) THEN
        DCMHP(I,J)=0.0D+00
      ENDIF
      ENDDO
      ENDDO
C
      DNORM = 0.0D+00
      DO I=1,NTST
        DO J=1,NTST
          DNORM = DNORM + DCMHP(I,J)
        ENDDO
        IF (DNORM .GT. 1.0D+00) THEN
          DO J=1,NTST
            DCMHP(I,J) = DCMHP(I,J)/DNORM
          ENDDO
        ENDIF
        DNORM = 0.0D+00
      ENDDO
C
      END SUBROUTINE FSSHPR
C*MODULE NAMD   *DECK FSSH
C>
C>     @brief Tully's fewest switches surface hopping algorithm
C>
C>     @author Seunghoon Lee
C>
      SUBROUTINE FSSH(DCMHP,DSTEG,RANDOM,NTST)
C
      USE comm_NONAD, ONLY: NAMD,NDSWCH,NDRST,NDTLF,COLD,THRSHE
      USE mx_limits, ONLY: mxatm
      USE constants, ONLY: zero
      USE prec, ONLY: dp
      IMPLICIT NONE
C
        REAL(KIND=dp) :: AU2AMU, AU2ANG, AU2KCAL, AU2SEC, EK2BAR,       &
     &                   EK2KCAL
      COMMON /CONVMD/ AU2KCAL, EK2KCAL, AU2ANG, AU2SEC, AU2AMU, EK2BAR
        REAL(KIND=dp), DIMENSION(1) :: X
      COMMON /FMCOM / X
        LOGICAL :: ALPHKWD, BETAKWD, MRDEA, MREKT, SG1T, TAMMD, TPA,    &
     &             TRIPLET
        REAL(KIND=dp) :: CNVTOL
        INTEGER :: IRECTD, ITDFG, ITDPRP, JANST, MAXVEC, MODTD, MTHST,  &
     &             MULTD, NLEBT, NONEQR, NPHIT, NRADT, NSTAT, NTHET,    &
     &             NTHST, NTRIAL
        INTEGER, DIMENSION(4) :: IFEDAT
        REAL(KIND=dp), DIMENSION(2) :: PFREQ
        REAL(KIND=dp), DIMENSION(3) :: SPCP
      COMMON /INFOTD/ CNVTOL, PFREQ, MODTD, JANST, NRADT, NTHET, NPHIT, &
     &                NLEBT, NSTAT, NTRIAL, MAXVEC, NTHST, IRECTD,      &
     &                ITDFG, ITDPRP, TRIPLET, SG1T, NONEQR, TAMMD, TPA, &
     &                ALPHKWD, BETAKWD, SPCP, MULTD, MREKT, MRDEA,      &
     &                MTHST, IFEDAT
        LOGICAL :: DSKWRK, GOPARR, MASWRK
        INTEGER :: IBTYP, IPTIM, MASTER, ME, NPROC
      COMMON /PAR   / ME, MASTER, NPROC, IBTYP, IPTIM, GOPARR, DSKWRK,  &
     &                MASWRK
      integer IR,IW,IP,IS,IPK,IDAF,NAV,IODA
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
C
      INTEGER :: NTST
      REAL(KIND=dp) :: RANDOM
      REAL(KIND=dp), DIMENSION(NTST,NTST) :: DCMHP
      REAL(KIND=dp), DIMENSION(NTST) :: DSTEG
C
      REAL(KIND=dp) :: CMPR, DKIN, DSTAE, FLAG1, FLAG2
      INTEGER :: I, ICMPR, LAST, LOADFM, NCRST, NEED
      LOGICAL :: LMRSF, MRSFS, MRSFT, NEWRND
C
      REAL(KIND=dp), PARAMETER :: EXCITE = transfer('EXCITE  ',1.0d0)
      REAL(KIND=dp), PARAMETER :: RNONE = transfer('NONE    ',1.0d0)
      REAL(KIND=dp), PARAMETER :: SPNFLP = transfer('SPNFLP  ',1.0d0)
C
C  DEFINE CURRENT STATE NUMBER
      NCRST = NTHST
C  COMPUTE CUMULATE PROBABILITY TO LINE UP
      CALL VALFM(LOADFM)
      ICMPR = LOADFM + 1
      LAST = ICMPR + NTST
      NEED = LAST - LOADFM-1
      CALL GETFM(NEED)
C
      DO I=1,NTST
        CALL CMHPPR(CMPR,NCRST,I,DCMHP,NTST)
        X(ICMPR+I-1)=CMPR
      ENDDO
C  COMPARE RANDOM NUMBER TO LINED UP PROBABILITY
      DO I=1,NTST
C  SETUP LINED UP PROBABILITY
         IF(NCRST.NE.1) THEN
            IF(I.EQ.1) THEN
               FLAG1=0.0D+00
               FLAG2=X(ICMPR+I-1)
            ELSEIF(I.NE.1.AND.I.NE.NCRST) THEN
               FLAG1=X(ICMPR+I-2)
               FLAG2=X(ICMPR+I-1)
            ELSEIF(I.EQ.NCRST) THEN
               CYCLE
            ENDIF
         ELSE
            IF(I.NE.NCRST) THEN
               FLAG1=X(ICMPR+I-2)
               FLAG2=X(ICMPR+I-1)
            ELSEIF(I.EQ.NCRST) THEN
               CYCLE
            ENDIF
         ENDIF
C  IF HOPPING OCCUR VELOCITY RESCALING
         IF(RANDOM.GT.FLAG1.AND.RANDOM.LT.FLAG2) THEN
            CALL MDQKIN(DKIN)
            DSTAE = DSTEG(NCRST)-DSTEG(I)
c     CHECK KINETIC ENERGY > EXCITED ENERGY
            IF(DSTAE.LT.ZERO.AND.DKIN.LT.ABS(DSTAE)) CYCLE
c     CHECK THRESHOLD OF ENERGY FOR HOPPING > ENERGY DIFFERENCE
            IF(THRSHE/AU2KCAL.LT.ABS(DSTAE)) THEN
               IF(MASWRK) WRITE(IW,1100)
     *                    NCRST,I,THRSHE/AU2KCAL,ABS(DSTAE)
               CYCLE
            ENDIF
C     HOPPING AND QUANTUM VELOCITY RESCALING
            NTHST=I
            CALL RESCALV(DKIN,DSTAE)
            IF (MASWRK) THEN
               WRITE(IW,1000) NCRST,NTHST
            ENDIF
         ENDIF
      ENDDO
C
      CALL RETFM(NEED)
 1000 format(/5x,'SURFACE HOPPING OCCURS FROM',I3,'TH STATE TO',
     *                                         I3,'TH STATE')
 1100 format(/5x,'SURFACE HOPPING FROM',I3,'TH STATE TO',
     *            I3,'TH STATE',1X,
     *       /5x,'IS BLOCKED BY THRSHE (',F10.5,') < |DELE| (',
     *            F10.5,')')
      END SUBROUTINE FSSH
C*MODULE NAMD   *DECK RESCALV
C>
C>     @brief Rescales velocities of nuclei after surface hopping
C>
C>     @author Seunghoon Lee
C>
      SUBROUTINE RESCALV(DKIN,DSTAE)
C
      USE mx_limits, ONLY: mxfrg, mxatm
      USE prec, ONLY: dp
      IMPLICIT NONE
C
        REAL(KIND=dp), DIMENSION(MXFRG) :: OX, OY,  OZ,  QW, QW1, QW2,  &
     &                                     QX, QX1, QX2, QY, QY1, QY2,  &
     &                                     QZ, QZ1, QZ2, VX, VY,  VZ
        REAL(KIND=dp), DIMENSION(MXATM) :: VXQM, VYQM, VZQM
      COMMON /ATMDAT/ VX, VY, VZ, QW, QX, QY, QZ, OX, OY, OZ, QW1, QX1, &
     &                QY1, QZ1, QW2, QX2, QY2, QZ2, VXQM, VYQM, VZQM
        REAL(KIND=dp), DIMENSION(3,MXATM) :: C
        INTEGER, DIMENSION(MXATM) :: IAN
        INTEGER :: ICH, MUL, NA, NAT, NB, NE, NQMT, NUM
        REAL(KIND=dp), DIMENSION(MXATM) :: ZAN
      COMMON /INFOA / NAT, ICH, MUL, NUM, NQMT, NE, NA, NB, ZAN, C, IAN
C
      REAL(KIND=dp) :: DKIN, DSTAE
C
      INTEGER :: I
      REAL(KIND=dp) :: SCALE1
C
      SCALE1 = 1+DSTAE/DKIN
      SCALE1 = SQRT(SCALE1)
C
      DO I = 1, NAT
         VXQM(I) = SCALE1*VXQM(I)
         VYQM(I) = SCALE1*VYQM(I)
         VZQM(I) = SCALE1*VZQM(I)
      ENDDO
      END SUBROUTINE RESCALV
C*MODULE NAMD   *DECK MDQKIN
C>
C>     @brief Computes the kinetic energy
C>
C>     @author Seunghoon Lee
C>
      SUBROUTINE MDQKIN(DQKIN)
C
      USE mx_limits, ONLY: mxfrg, mxatm
      USE prec, ONLY: dp
      IMPLICIT NONE
C
        REAL(KIND=dp), DIMENSION(MXFRG) :: OX, OY,  OZ,  QW, QW1, QW2,  &
     &                                     QX, QX1, QX2, QY, QY1, QY2,  &
     &                                     QZ, QZ1, QZ2, VX, VY,  VZ
        REAL(KIND=dp), DIMENSION(MXATM) :: VXQM, VYQM, VZQM
      COMMON /ATMDAT/ VX, VY, VZ, QW, QX, QY, QZ, OX, OY, OZ, QW1, QX1, &
     &                QY1, QZ1, QW2, QX2, QY2, QZ2, VXQM, VYQM, VZQM
        REAL(KIND=dp) :: AU2AMU, AU2ANG, AU2KCAL, AU2SEC, EK2BAR,       &
     &                   EK2KCAL
      COMMON /CONVMD/ AU2KCAL, EK2KCAL, AU2ANG, AU2SEC, AU2AMU, EK2BAR
        REAL(KIND=dp), DIMENSION(3,MXATM) :: C
        INTEGER, DIMENSION(MXATM) :: IAN
        INTEGER :: ICH, MUL, NA, NAT, NB, NE, NQMT, NUM
        REAL(KIND=dp), DIMENSION(MXATM) :: ZAN
      COMMON /INFOA / NAT, ICH, MUL, NUM, NQMT, NE, NA, NB, ZAN, C, IAN
        REAL(KIND=dp), DIMENSION(MXATM) :: ZMASS
      COMMON /MASSES/ ZMASS
C
      REAL(KIND=dp) :: DQKIN
C
      REAL(KIND=dp) :: DQKINX, DQKINY, DQKINZ, ZMAS
      INTEGER :: I
C
      DQKIN  = 0.0D+00
      DQKINX = 0.0D+00
      DQKINY = 0.0D+00
      DQKINZ = 0.0D+00
c  TRANSLATION ENERGY. (COMPUTE IN ATOMIC UNIT)
      DO I = 1, NAT
         ZMAS  = ZMASS(I)
         DQKINX = DQKINX + ZMAS*VXQM(I)**2
         DQKINY = DQKINY + ZMAS*VYQM(I)**2
         DQKINZ = DQKINZ + ZMAS*VZQM(I)**2
      ENDDO
      DQKIN  = 0.5D+00*(DQKINX+DQKINY+DQKINZ)

c  CONVERT TO AU
      DQKIN = DQKIN * EK2KCAL / AU2KCAL
C
      END SUBROUTINE MDQKIN
C*MODULE NAMD   *DECK CMHPPR
C>
C>     @brief Computes cumulative hopping probabilities
C>
C>     @author Seunghoon Lee
C>
      SUBROUTINE CMHPPR(CMPR,NCRST,K,DCMHP,NTST)
C
      USE prec, ONLY: dp
      IMPLICIT NONE
C
      REAL(KIND=dp) :: CMPR
      INTEGER :: K, NCRST, NTST
      REAL(KIND=dp), DIMENSION(NTST,NTST) :: DCMHP
C
      INTEGER :: I
C
      CMPR = 0.0D+00
      DO I=1,K
         CMPR=CMPR+DCMHP(NCRST,I)
      ENDDO
C
      END SUBROUTINE CMHPPR
C*MODULE SFDFT   *DECK OVKSMO
C>
C>     @brief Compute overlap integrals between KS MOs
C>            at different time steps
C>
C>     @author Seunghoon Lee
C>
      SUBROUTINE OVKSMO(DKSMSA,DKSMSB,L1,L2,L3,LX)
C
      USE comm_NONAD, ONLY: NAMD,NDSWCH,NDRST,NDTLF,COLD,THRSHE
      USE mx_limits, ONLY: mxatm, mxgtot, mxsh
      USE prec, ONLY: dp
      IMPLICIT NONE
C
        REAL(KIND=dp), DIMENSION(1) :: X
      COMMON /FMCOM / X
        REAL(KIND=dp), DIMENSION(3,MXATM) :: C
        INTEGER, DIMENSION(MXATM) :: IAN
        INTEGER :: ICH, MUL, NA, NAT, NB, NE, NQMT, NUM
        REAL(KIND=dp), DIMENSION(MXATM) :: ZAN
      COMMON /INFOA / NAT, ICH, MUL, NUM, NQMT, NE, NA, NB, ZAN, C, IAN
        INTEGER :: IDAF, IJKO, IJKT, IP, IR, IW, NAV
        INTEGER, DIMENSION(950) :: IODA
      COMMON /IOFILE/ IR, IW, IP, IJKO, IJKT, IDAF, NAV, IODA
        REAL(KIND=dp), DIMENSION(MXGTOT) :: CD, CF, CG, CH, CI, CP, CS, &
     &                                      EX
        INTEGER, DIMENSION(MXSH) :: KATOM, KLOC, KMAX, KMIN, KNG,       &
     &                              KSTART, KTYPE
        INTEGER :: NSHELL
      COMMON /NSHEL / EX, CS, CP, CD, CF, CG, CH, CI, KSTART, KATOM,    &
     &                KTYPE, KNG, KLOC, KMIN, KMAX, NSHELL
C
      INTEGER :: L1, L2, L3, LX
      REAL(KIND=dp), DIMENSION(LX*LX) :: DKSMSA, DKSMSB
C
      REAL(KIND=dp) :: DUMMY, SIGNA, SIGNB
      INTEGER :: I, IJ, INDAS, ITMP1, IVECNA, IVECNB, IVECOA, IVECOB, J,&
     &           L3T, LAST, LOADFM, NEED
      LOGICAL :: NEWRND
c
C        MEMORY ALLOCATION (AO 1ELEC INTEG, NEW AND OLD COEFFI)
C
      L3T=LX*L1
      CALL VALFM(LOADFM)
      IVECOA = LOADFM + 1
      IVECOB = IVECOA + L3T
      IVECNA = IVECOB + L3T
      IVECNB = IVECNA + L3T
      INDAS = IVECNB + L3T
      ITMP1 = INDAS + L3
      LAST = ITMP1 + L3
      NEED = LAST - LOADFM - 1
      CALL GETFM(NEED)

      CALL VCLR(X(IVECOA),1,NEED)
c
C        GET AO 1ELECTRON INTEGRAL MATRIX NDAS (Non-Adiabatic Atomic S matrix)
C          <OLD(a)|NEW(b')>

C      CALL VCLR(X(INDAS),1,LX*LX)    ! compute overlap of two different GTO basis sets in GUESS
      CALL COOVLP(1,DUMMY,X(INDAS),L1,L1,L2,NAT,
     *        MXGTOT,NSHELL,
     *        EX,CS,CP,CD,CF,CG,CH,CI,
     *        KSTART,KATOM,KTYPE,KNG,KLOC,KMIN,KMAX,
     *        MXGTOT,NSHELL,
     *        EX,CS,CP,CD,CF,CG,CH,CI,
     *        KSTART,KATOM,KTYPE,KNG,KLOC,KMIN,KMAX,
     *        C,COLD)
Cvap
C      CALL VCLR(X(IVECNA),1,L3T)
C      CALL VCLR(X(IVECNB),1,L3T)
C      CALL VCLR(X(IVECOA),1,L3T)
C      CALL VCLR(X(IVECOB),1,L3T)
Cvap
      CALL DAREAD(IDAF,IODA,X(IVECOA),L3T,700,0)
      CALL DAREAD(IDAF,IODA,X(IVECOB),L3T,701,0)
      CALL DAREAD(IDAF,IODA,X(IVECNA),L3T,15,0)
      CALL DAREAD(IDAF,IODA,X(IVECNB),L3T,19,0)
C
C        EIGENVECTOR SIGN CORRECTION OF DIFFERENT TIME
      DO I=1,LX
         SIGNA=0.D+00
         SIGNB=0.D+00
         DO J=1,L1
           IJ=L1*(I-1)+J
           SIGNA=SIGNA+X(IVECOA+IJ-1)*X(IVECNA+IJ-1)
           SIGNB=SIGNB+X(IVECOB+IJ-1)*X(IVECNB+IJ-1)
         ENDDO
C
         IF(SIGNA.LE.0.D+00) THEN
           SIGNA=-1.D+00
           DO J=1,L1
             IJ=L1*(I-1)+J
             X(IVECNA+IJ-1) = -1.D+00 * X(IVECNA+IJ-1)
           ENDDO
         ENDIF
         IF(SIGNB.LE.0.D+00) THEN
           SIGNB=-1.D+00
           DO J=1,L1
             IJ=L1*(I-1)+J
             X(IVECNB+IJ-1) = -1.D+00 * X(IVECNB+IJ-1) ! BUG IVECNA->IVECNB ?
           ENDDO
         ENDIF
      ENDDO
C
      CALL MRTRBR(X(IVECOA),L1,L1,LX,X(INDAS),L1,L1,X(ITMP1),LX) !a C = A^ * B
      CALL MRARBR(X(ITMP1),LX,LX,L1,X(IVECNA),L1,LX,DKSMSA,LX)   !a C = A * B
      CALL MRTRBR(X(IVECOB),L1,L1,LX,X(INDAS),L1,L1,X(ITMP1),LX) !b C = A^ * B
      CALL MRARBR(X(ITMP1),LX,LX,L1,X(IVECNB),L1,LX,DKSMSB,LX)   !b C = A * B
C
      CALL RETFM(NEED)

      END SUBROUTINE OVKSMO
C*MODULE NAMD   *DECK NDRSINI
C>
C>     @brief Read restarting information of NAMD
C>
C>     @author Seunghoon Lee
C>
      SUBROUTINE NDRSINI(TDCOE,DRND,NSTEPRND,NEWRND)
C
      USE comm_NONAD, ONLY: NAMD,NDSWCH,NDRST,NDTLF,COLD,THRSHE
      USE mx_limits, ONLY: mxatm
      USE prec, ONLY: dp
      IMPLICIT NONE
C
        REAL(KIND=dp), DIMENSION(3,MXATM) :: C
        INTEGER, DIMENSION(MXATM) :: IAN
        INTEGER :: ICH, MUL, NA, NAT, NB, NE, NQMT, NUM
        REAL(KIND=dp), DIMENSION(MXATM) :: ZAN
      COMMON /INFOA / NAT, ICH, MUL, NUM, NQMT, NE, NA, NB, ZAN, C, IAN
        LOGICAL :: ALPHKWD, BETAKWD, MRDEA, MREKT, SG1T, TAMMD, TPA,    &
     &             TRIPLET
        REAL(KIND=dp) :: CNVTOL
        INTEGER :: IRECTD, ITDFG, ITDPRP, JANST, MAXVEC, MODTD, MTHST,  &
     &             MULTD, NLEBT, NONEQR, NPHIT, NRADT, NSTAT, NTHET,    &
     &             NTHST, NTRIAL
        INTEGER, DIMENSION(4) :: IFEDAT
        REAL(KIND=dp), DIMENSION(2) :: PFREQ
        REAL(KIND=dp), DIMENSION(3) :: SPCP
      COMMON /INFOTD/ CNVTOL, PFREQ, MODTD, JANST, NRADT, NTHET, NPHIT, &
     &                NLEBT, NSTAT, NTRIAL, MAXVEC, NTHST, IRECTD,      &
     &                ITDFG, ITDPRP, TRIPLET, SG1T, NONEQR, TAMMD, TPA, &
     &                ALPHKWD, BETAKWD, SPCP, MULTD, MREKT, MRDEA,      &
     &                MTHST, IFEDAT
        INTEGER :: IDAF, IJKO, IJKT, IP, IR, IW, NAV
        INTEGER, DIMENSION(950) :: IODA
      COMMON /IOFILE/ IR, IW, IP, IJKO, IJKT, IDAF, NAV, IODA
        LOGICAL :: DSKWRK, GOPARR, MASWRK
        INTEGER :: IBTYP, IPTIM, MASTER, ME, NPROC
      COMMON /PAR   / ME, MASTER, NPROC, IBTYP, IPTIM, GOPARR, DSKWRK,  &
     &                MASWRK
C
      LOGICAL :: NEWRND
      INTEGER :: NSTEPRND
      REAL(KIND=dp), DIMENSION(NSTEPRND) :: DRND
      REAL(KIND=dp), DIMENSION(2,NSTAT) :: TDCOE
C
      INTEGER :: INIT, L1, LX
      LOGICAL :: LMRSF, MRSFS, MRSFT
C
      REAL(KIND=dp), PARAMETER :: EXCITE = transfer('EXCITE  ',1.0d0)
      REAL(KIND=dp), PARAMETER :: RNONE = transfer('NONE    ',1.0d0)
      REAL(KIND=dp), PARAMETER :: SPNFLP = transfer('SPNFLP  ',1.0d0)
C
      LX = NQMT
      L1 = NUM
C
      IF(MASWRK) THEN
         CALL RREADND(IR, IW,' $TDC',2*NSTAT,TDCOE)
         IF(NEWRND) THEN
            INIT=0
            CALL RNGEN(DRND,NSTEPRND,3,INIT)
            INIT=1
         ELSE
            CALL RREADND(IR, IW,' $RND',NSTEPRND,DRND)
         END IF
         CALL RREADND(IR, IW,' $CDO',3*NAT,COLD)
      ENDIF
C
      IF (MASWRK) THEN
        CALL DDI_BCAST(2602,'F',TDCOE,2*NSTAT,MASTER)
        CALL DDI_BCAST(2608,'F',DRND,NSTEPRND,MASTER)
        CALL DDI_BCAST(2610,'F',COLD,3*NAT,MASTER)
      ELSE
        CALL DDI_BCAST(2602,'F',TDCOE,2*NSTAT,MASTER)
        CALL DDI_BCAST(2608,'F',DRND,NSTEPRND,MASTER)
        CALL DDI_BCAST(2610,'F',COLD,3*NAT,MASTER)
      ENDIF
C
      CALL DAWRIT(IDAF,IODA,TDCOE,2*NSTAT,707,0)
      CALL DAWRIT(IDAF,IODA,NTHST,1      ,708,1)
      CALL DAWRIT(IDAF,IODA,DRND,NSTEPRND,709,0)

      IF(MASWRK) WRITE(IW,9001)
C
 9000 FORMAT(1X,'MEMORY REQUIRED FOR RESTART INITIALIZATION IS=',I10,
     *          ' WORDS.'/)
 9001 FORMAT(1X,'..... DONE WITH RESTART INITIALIZATION .....')
      END SUBROUTINE NDRSINI
C*MODULE NAMD   *DECK SFPRTCORVEL
C>
C>     @brief clone of MDPRTCORVEL
C>            write output data for NAMD
C>
C>     @author Seunghoon Lee
C>
      SUBROUTINE SFPRTCORVEL(ISTEP,TTOTAL,INTMD,iptraj,NTHSTO,TDC,RND,  &
     &                       L1,LX,NEWRND)
C
      USE comm_NONAD, ONLY: NAMD,NDSWCH,NDRST,NDTLF,COLD,THRSHE
      USE mx_limits, ONLY: mxatm, mxao, mxfrg, mxfgpt
      USE constants, ONLY: zero
      USE comm_EFMULT
      USE comm_FRGMSS
      USE comm_EFPFRC
      USE comm_FGRAD
      USE comm_FRGINF
      USE prec, ONLY: dp
      IMPLICIT NONE
C
      INTEGER, PARAMETER :: NREPCA = 100
      INTEGER, PARAMETER :: NBIN = 1000
      INTEGER, PARAMETER :: MAXRDF = 100
C
        REAL(KIND=dp), DIMENSION(MXFRG) :: OX, OY,  OZ,  QW, QW1, QW2,  &
     &                                     QX, QX1, QX2, QY, QY1, QY2,  &
     &                                     QZ, QZ1, QZ2, VX, VY,  VZ
        REAL(KIND=dp), DIMENSION(MXATM) :: VXQM, VYQM, VZQM
      COMMON /ATMDAT/ VX, VY, VZ, QW, QX, QY, QZ, OX, OY, OZ, QW1, QX1, &
     &                QY1, QZ1, QW2, QX2, QY2, QZ2, VXQM, VYQM, VZQM
        REAL(KIND=dp) :: AU2AMU, AU2ANG, AU2KCAL, AU2SEC, EK2BAR,       &
     &                   EK2KCAL
      COMMON /CONVMD/ AU2KCAL, EK2KCAL, AU2ANG, AU2SEC, AU2AMU, EK2BAR
        REAL(KIND=dp), DIMENSION(nrEPCA) :: ABATHT
        LOGICAL :: CCMS, PROD, RSRAND, RSTEMP, SSBP, USAMP
        REAL(KIND=dp) :: CFORCE, DTEMP, MDINT, SFORCE, TMIN, UV
        REAL(KIND=dp), DIMENSION(2) :: DROFF
        INTEGER, DIMENSION(50) :: INDEXOH, IPAIR
        INTEGER, DIMENSION(200) :: ISOLUT, MOVEMM
        INTEGER :: IUSTYP, JEVERY, KEVERY, LEVERY, MREMD, NCST, NEFPMV, &
     &             NRAND, NVTOFF
        REAL(KIND=dp), DIMENSION(20) :: PRVEC
        REAL(KIND=dp), DIMENSION(10) :: RZERO, UFORCE
        REAL(KIND=dp), DIMENSION(3) :: TMASS
      COMMON /CSTPOT/ DROFF, SFORCE, UFORCE, RZERO, DTEMP, TMIN, UV,    &
     &                ABATHT, CFORCE, TMASS, PRVEC, MDINT, IPAIR,       &
     &                IUSTYP, NCST, NRAND, NVTOFF, JEVERY, KEVERY,      &
     &                LEVERY, MREMD, INDEXOH, MOVEMM, NEFPMV, RSTEMP,   &
     &                RSRAND, PROD, SSBP, CCMS, USAMP, ISOLUT
        REAL(KIND=dp), DIMENSION(1) :: X
      COMMON /FMCOM / X
        CHARACTER(6), DIMENSION(MXFRG) :: FRGNAM
        CHARACTER(8), DIMENSION(50,MXFRG) :: PTNAM
        REAL(KIND=dp), DIMENSION(50,MXFRG) :: XCRD, YCRD, ZCRD
      COMMON /FRAGMT/ XCRD, YCRD, ZCRD, PTNAM, FRGNAM
        REAL(KIND=dp) :: E
        REAL(KIND=dp), DIMENSION(3*MXATM) :: EG
      COMMON /FUNCT / E, EG
        LOGICAL :: INITGDDI, ISGDDI, MLGDDI, PAROUT, WASGDDI
        INTEGER :: ISCOPE, JBTYP, MEGLOB, MEUNIV, MYGROUP, MYWORLD,     &
     &             NGROUPS, NNGLOB, NPGLOB, NPUNIV, NSUBGR, NUMDLB,     &
     &             NWORLDS, MOGDDI
      COMMON /GDDI  / ISCOPE, NGROUPS, MYGROUP, MEGLOB, NPGLOB, NNGLOB, &
     &                JBTYP, ISGDDI, PAROUT, INITGDDI, WASGDDI, MLGDDI, &
     &                NSUBGR, MEUNIV, NPUNIV, NUMDLB, MYWORLD, NWORLDS, &
     &                MOGDDI 
        REAL(KIND=dp), DIMENSION(NBIN) :: CCRDF
        REAL(KIND=dp) :: DELR
        REAL(KIND=dp), DIMENSION(NBIN,2) :: GEN2FR
        REAL(KIND=dp), DIMENSION(NBIN,MAXRDF) :: GENFR
        REAL(KIND=dp), DIMENSION(NBIN,3) :: GOFR
        INTEGER :: NUMG, NUMSAM1, NUMSAM2, NUMSAM3
      COMMON /GOFR  / GOFR, GENFR, CCRDF, GEN2FR, DELR, NUMG, NUMSAM1,  &
     &                NUMSAM2, NUMSAM3
        REAL(KIND=dp), DIMENSION(3,MXATM) :: C
        INTEGER, DIMENSION(MXATM) :: IAN
        INTEGER :: ICH, MUL, NA, NAT, NB, NE, NQMT, NUM
        REAL(KIND=dp), DIMENSION(MXATM) :: ZAN
      COMMON /INFOA / NAT, ICH, MUL, NUM, NQMT, NE, NA, NB, ZAN, C, IAN
        LOGICAL :: ALPHKWD, BETAKWD, MRDEA, MREKT, SG1T, TAMMD, TPA,    &
     &             TRIPLET
        REAL(KIND=dp) :: CNVTOL
        INTEGER :: IRECTD, ITDFG, ITDPRP, JANST, MAXVEC, MODTD, MTHST,  &
     &             MULTD, NLEBT, NONEQR, NPHIT, NRADT, NSTAT, NTHET,    &
     &             NTHST, NTRIAL
        INTEGER, DIMENSION(4) :: IFEDAT
        REAL(KIND=dp), DIMENSION(2) :: PFREQ
        REAL(KIND=dp), DIMENSION(3) :: SPCP
      COMMON /INFOTD/ CNVTOL, PFREQ, MODTD, JANST, NRADT, NTHET, NPHIT, &
     &                NLEBT, NSTAT, NTRIAL, MAXVEC, NTHST, IRECTD,      &
     &                ITDFG, ITDPRP, TRIPLET, SG1T, NONEQR, TAMMD, TPA, &
     &                ALPHKWD, BETAKWD, SPCP, MULTD, MREKT, MRDEA,      &
     &                MTHST, IFEDAT
        INTEGER :: IDAF, IJKO, IJKT, IP, IR, IW, NAV
        INTEGER, DIMENSION(950) :: IODA
      COMMON /IOFILE/ IR, IW, IP, IJKO, IJKT, IDAF, NAV, IODA
        REAL(KIND=dp) :: BATHT, DT, EKINR, EKINT
        INTEGER :: MDTYP, PRESSMD
      COMMON /MDSIM / DT, BATHT, EKINT, EKINR, MDTYP, PRESSMD
        INTEGER :: NSTEPS
      COMMON /MDSIM4/ NSTEPS
        REAL(KIND=dp), DIMENSION(3,MXFRG) :: TQI
        REAL(KIND=dp), DIMENSION(MXFRG) :: WX, WY, WZ
      COMMON /MDVEL / WX, WY, WZ, TQI
        REAL(KIND=dp) :: EKT0, GFREE0
        REAL(KIND=dp), DIMENSION(2) :: GNH, QNH, VNH, XNH
        INTEGER :: NVTNH
      COMMON /NOSEHO/ GNH, VNH, XNH, QNH, EKT0, GFREE0, NVTNH
        LOGICAL :: DSKWRK, GOPARR, MASWRK
        INTEGER :: IBTYP, IPTIM, MASTER, ME, NPROC
      COMMON /PAR   / ME, MASTER, NPROC, IBTYP, IPTIM, GOPARR, DSKWRK,  &
     &                MASWRK
        LOGICAL :: PBCOUT, PBCWRI
      COMMON /PBC   / PBCOUT, PBCWRI
        REAL(KIND=dp), DIMENSION(3,MXFGPT) :: CREPPB, EFDPPB, EFPPB,    &
     &                                        FRGCPB
      COMMON /PBCCO / FRGCPB, EFPPB, EFDPPB, CREPPB
        REAL(KIND=dp), DIMENSION(MXATM) :: ANAM, BNAM
        REAL(KIND=dp), DIMENSION(MXAO) :: BFLAB
        REAL(KIND=dp), DIMENSION(10) :: TITLE
      COMMON /RUNLAB/ TITLE, ANAM, BNAM, BFLAB
C
      INTEGER :: INTMD, IPTRAJ, ISTEP, L1, LX, NTHSTO
      LOGICAL :: NEWRND
      REAL(KIND=dp) :: TTOTAL
      REAL(KIND=dp), DIMENSION(NSTEPS+1) :: RND
      REAL(KIND=dp), DIMENSION(2*NSTAT) :: TDC
C
      CHARACTER(3) :: CHROUT, CHROUTND
      REAL(KIND=dp) :: CONVFAC, EKIN, EPOT, ETOT, F2F, TTOTALRS, VFACT, &
     &                 XF, Y, YF, Z, ZF, ZNUC
      INTEGER :: I, IAT, III, IN, IPNDRS, IPREMD, J, K, KMASS, N, NQMMM
      LOGICAL :: LMRSF, MRSFQ, MRSFS, MRSFT
C
      REAL(KIND=dp), PARAMETER :: EXCITE = transfer('EXCITE  ',1.0d0)
      REAL(KIND=dp), PARAMETER :: RNONE = transfer('NONE    ',1.0d0)
      REAL(KIND=dp), PARAMETER :: SPNFLP = transfer('SPNFLP  ',1.0d0)
c
      IF(NAMD) THEN
        IPNDRS  =45
        IF(IPTRAJ.EQ.IPNDRS) THEN
          TTOTALRS=TTOTAL*1.0D+15
          IF(.NOT.(NDRST .OR. TTOTALRS.GT.0.0D+00)) GO TO 200
        END IF
      END IF
c
c        OUTPUT PARTICLE COORDINATES AND VELOCITY INFORMATION
c                 USED FOR qm, qm/efp, OR efp RUNS
c
      if(.not.maswrk) return
      ipndrs =45
      ipremd =44
      IF ((IPTRAJ.NE.IPREMD).AND.(IPTRAJ.NE.IPNDRS)) THEN
c
c        FIRST SECTION IS PRINTING PARTICLE COORDINATES TO LOG FILE
c
      IF(NAT.GT.0) WRITE(IW,1000)
      DO IAT=1,NAT
         ZNUC = IAN(IAT)
         X = AU2ANG * C(1,IAT)
         Y = AU2ANG * C(2,IAT)
         Z = AU2ANG * C(3,IAT)
         WRITE(IW,1005) ANAM(IAT),BNAM(IAT),ZNUC,X,Y,Z
      ENDDO
C
      IF(NFRG.GT.0) WRITE(IW,1009)
      IAT=0
      KMASS=0
      DO I = 1,NFRG
         WRITE(IW,1010) FRGNAM(I), I
         DO III=1,NMPTS(I)
            IAT=IAT+1
            IF(FMASS(III+KMASS).GT.ZERO) THEN
               XF = EFC(1,IAT)*AU2ANG
               YF = EFC(2,IAT)*AU2ANG
               ZF = EFC(3,IAT)*AU2ANG
               WRITE(IW,1020) FRGNME(IAT),XF,YF,ZF
            END IF
         END DO
         KMASS=KMASS+NMPTS(I)
      ENDDO
c
c         OPTION TO REPORT efp PARTICLES SHIFTED INSIDE A SINGLE BOX
c
      IF (ISTEP.EQ.NSTEPS  .AND.  PBCWRI) THEN
         IF(NFRG.GT.0) WRITE(IW,1011)
         IAT=0
         KMASS=0
         DO I= 1,NFRG
            WRITE(IW,1010) FRGNAM(I), I
            DO III=1,NMPTS(I)
               IAT=IAT+1
               IF(FMASS(III+KMASS).GT.ZERO) THEN
                  XF =FRGCPB(1,IAT)*AU2ANG
                  YF =FRGCPB(2,IAT)*AU2ANG
                  ZF =FRGCPB(3,IAT)*AU2ANG
                  WRITE(IW,1020) FRGNME(IAT),XF,YF,ZF
               END IF
            END DO
            KMASS=KMASS+NMPTS(I)
         ENDDO
      END IF
C
      ENDIF
c
c        SECOND SECTION IS COORDINATE OUTPUT, TO THE TRAJECTORY FILE
c
      EPOT = E*AU2KCAL
      EKIN = EKINT+EKINR
      ETOT = EPOT+EKIN
      NQMMM=0
      WRITE(IPTRAJ,2000) NAT,NFRG,NQMMM,1.0D+15*TTOTAL,ETOT,
     *                  EPOT,batht,EKIN,EKINT,EKINR
C
      IF(NAT.GT.0) WRITE(IPTRAJ,2001)
      DO I=1,NAT
         ZNUC = IAN(I)
         WRITE(IPTRAJ,2010) ANAM(I),BNAM(I),ZNUC,(AU2ANG*C(J,I),J=1,3)
      ENDDO
C
      IAT=0
      KMASS=0
      IF(NFRG.GT.0) WRITE(IPTRAJ,2002)
      write(IPTRAJ,2012)
      DO I=1,NFRG
         WRITE(IPTRAJ,2020) FRGNAM(I),I
         DO III=1,NMPTS(I)
            IAT=IAT+1
            IF(FMASS(III+KMASS).GT.ZERO) THEN
               WRITE(IPTRAJ,2030) FRGNME(IAT),(AU2ANG*EFC(J,IAT),J=1,3)
            END IF
         END DO
         KMASS=KMASS+NMPTS(I)
      ENDDO
c
      WRITE(IPTRAJ,2013)
c
c         THIRD SECTION IS GRADIENT DATA, TO THE TRAJECTORY FILE
c     efp GRADIENT OUTPUT IS NEEDED FOR efp COARSE GRAINING.
c
      F2F = -AU2KCAL / EK2KCAL
      IF(NFRG.GT.0) WRITE(IPTRAJ,2050)
      N=0
      DO I = 1, NFRG
         WRITE(IPTRAJ,2051) I, FRGNAM(I)
         WRITE(IPTRAJ,2052) EFCENT(1,I),EFCENT(2,I),EFCENT(3,I)
         WRITE(IPTRAJ,2053) DEFT(1,I)/F2F,DEFT(2,I)/F2F,DEFT(3,I)/F2F
         WRITE(IPTRAJ,2054) TORQ(1,I)/F2F,TORQ(2,I)/F2F,TORQ(3,I)/F2F
         IF (IFRCPNT) THEN
           DO K=1, NFRCPNT(I)
             N=N+1
             WRITE(IPTRAJ,2055)I,K,FRCCRD(1,N),FRCCRD(2,N),FRCCRD(3,N),
     *                (FRCTRQ(IN,N),IN=1,6)
           END DO
          END IF
      ENDDO
 2050 FORMAT(5X,'GRADIENT DATA (NOT USED BY RESTARTS)...')
 2051 FORMAT('FRAGMENT #',I6,'  ',A6)
 2052 FORMAT('    EFCENT',3F16.10)
 2053 FORMAT('     FORCE',3F16.10)
 2054 FORMAT('      TORQ',3F16.10)
 2055 FORMAT(1X,2I4,9F16.10)
c
c         4TH SECTION IS VELOCITY RESTART DATA, TO THE TRAJECTORY FILE
c
      WRITE(IPTRAJ,3000) TTOTAL
      WRITE(IPTRAJ,3001) 'VVERLET',DT,NVTNH,NSTEPS-ISTEP
      write(chrout,'(A3)') '.F.'
      if (rstemp) write(chrout,'(A3)') '.T.'
      WRITE(IPTRAJ,3002) chrout,dtemp,levery
      write(chrout,'(A3)') '.F.'
      if (rsrand) write(chrout,'(A3)') '.T.'
      WRITE(IPTRAJ,3003) chrout,nrand,nvtoff,jevery
      write(chrout,'(A3)') '.F.'
      if (prod) write(chrout,'(A3)') '.T.'
      WRITE(IPTRAJ,3004) chrout,kevery,delr
      write(chroutnd,'(A3)') '.f.'
      if (NAMD) write(chroutnd,'(A3)') '.t.'
      write(chrout,'(A3)') '.f.'
      if (newrnd) write(chrout,'(A3)') '.t.'
      WRITE(IPTRAJ,3006) chroutnd,ndtlf,'.t.',chrout
c
      IF (ISGDDI) THEN
         WRITE(IPTRAJ,'(a,$)') " BathT(1)="
         DO I=1,NGROUPS
            IF (MOD(I,5).NE.0) THEN
               IF (I.EQ.NGROUPS) THEN
                  write(iptraj,"(1X,f6.2)") ABATHT(I)
               ELSE
                  write(iptraj,"(1X,f6.2,$)") ABATHT(I)
               ENDIF
            ELSE
               write(iptraj,"(1X,f6.2)") ABATHT(I)
            ENDIF
         ENDDO
      ELSE
         write(iptraj,3005) batht
      ENDIF
c
      CONVFAC=AU2KCAL/AU2ANG**2
      IF(SSBP) WRITE(iptraj,3010) SFORCE*CONVFAC,DROFF(1)*AU2ANG,
     *                            DROFF(2)*AU2ANG
      IF(CCMS) WRITE(iptraj,3019) CFORCE*CONVFAC
      IF(USAMP) THEN
         WRITE(iptraj,3011) UFORCE(1)*CONVFAC,UFORCE(2)*CONVFAC,
     *    UFORCE(3)*CONVFAC,RZERO(1)*AU2ANG,RZERO(2)*AU2ANG,
     *    RZERO(3)*AU2ANG
         WRITE(iptraj,3018) IUSTYP
          IF (IUSTYP.EQ.0) THEN
             IF (NCST.EQ.1) THEN
                write(iptraj,3012) IPAIR(1),IPAIR(2)
             ELSEIF (NCST.EQ.2) THEN
                write(iptraj,3013) IPAIR(1),IPAIR(2),IPAIR(3),IPAIR(4)
             ENDIF
          ELSEIF (IUSTYP.EQ.1) THEN
             IF (NCST.EQ.1) THEN
                write(iptraj,3014) IPAIR(1),IPAIR(2),IPAIR(3)
             ELSEIF (NCST.EQ.2) THEN
                write(iptraj,3015) IPAIR(1),IPAIR(2),IPAIR(3),IPAIR(4),
     *             IPAIR(5),IPAIR(6)
             ELSEIF (NCST.EQ.3) THEN
                write(iptraj,3016) IPAIR(1),IPAIR(2),IPAIR(3),IPAIR(4),
     *             IPAIR(5),IPAIR(6),IPAIR(7),IPAIR(8),IPAIR(9)
             ELSEIF (NCST.EQ.4) THEN
                write(iptraj,3017) IPAIR(1),IPAIR(2),IPAIR(3),IPAIR(4),
     *             IPAIR(5),IPAIR(6),IPAIR(7),IPAIR(8),IPAIR(9),
     *             IPAIR(10),IPAIR(11),IPAIR(12)
             ENDIF
         ENDIF
      ENDIF
C
      vfact = 1.0d-12
      if(nat.gt.0) then
         write(iptraj,3020)
         do i=1,nat
            write(iptraj,3030) vxqm(i)*vfact,vyqm(i)*vfact,vzqm(i)*vfact
         enddo
      end if
      if(nfrg.gt.0) then
         write(iptraj,3050)
         do i=1,nfrg
            write(iptraj,3030) vx(i)*vfact, vy(i)*vfact, vz(i)*vfact
         enddo
         write(iptraj,3060)
         do i=1,nfrg
            write(iptraj,3040) qw(i),qx(i),qy(i),qz(i)
         enddo
c
c           dISTINGUISH BETWEEN LEAPFROG AND VELOCITY VERLET...
c
         if(intmd.eq.1) then
            write(iptraj,3070)
            do i=1,nfrg
               write(iptraj,3040) ox(i)*vfact, oy(i)*vfact, oz(i)*vfact
            enddo
         else
            write(iptraj,3075)
            do i=1,nfrg
               write(iptraj,3030) wx(i)*vfact, wy(i)*vfact, wz(i)*vfact
            enddo
            write(iptraj,3080)
            do i=1,nfrg
               write(iptraj,3040) qw1(i),qx1(i),qy1(i),qz1(i)
            enddo
            write(iptraj,3090)
            do i=1,nfrg
               write(iptraj,3040) qw2(i),qx2(i),qy2(i),qz2(i)
            enddo
         endif
      end if
      write(iptraj,3095)
c
c         5TH SECTION IS EIGEN-VECTORS AND x+y RESTART DATA, TO THE RST FILE
c
      if(NAMD) then
        IF(IPTRAJ.EQ.IPNDRS) THEN
          WRITE(IPTRAJ,3110) NSTAT, NTHSTO, MULTD, '.T.'
          WRITE(IPTRAJ,3100) ' $TDC   '
          CALL NDPUSQLF(IPTRAJ,TDC,2*NSTAT,1)
          WRITE(IPTRAJ,3100) ' $END   '
          WRITE(IPTRAJ,3100) ' $RND   '
          CALL NDPUSQLF(IPTRAJ,RND,NSTEPS+1,1)
          WRITE(IPTRAJ,3100) ' $END   '
          WRITE(IPTRAJ,3100) ' $CDO   '
          CALL NDPUSQLF(IPTRAJ,COLD,3*NAT,1)
          WRITE(IPTRAJ,3100) ' $END   '
        END IF
      END IF
 200  CONTINUE
c
c        ASK UNIX TO FLUSH ITS BUFFERS TO DISK
c
      RETURN
c
 1000 FORMAT(1X,'QM ATOM COORDINATES (ANG)')
 1005 FORMAT(1X,A8,A2,F5.1,3F15.10)
 1009 FORMAT(11X,'CARTESIAN COORDINATES (ANG)')
 1010 FORMAT(1X,'FRAGNAME=',A6,'  !',I4)
 1011 FORMAT(11X,'PBC CARTESIAN COORDINATES (ANG)')
 1020 FORMAT(1X,A8,7X,3F18.12)
C
 2000 FORMAT('===== MD DATA PACKET ====='/
     *       'NAT=',I8,' NFRG=',I8,' NQMMM=',I8/
     *       'TTOTAL=',F12.2,' FS    TOT. E=',F20.6,' KCAL/MOL'/
     *       'POT. E= ',F20.6,' KCAL/MOL  BATHT= ',F20.6/
     *       'KIN. E= ',F20.6,'  TRANS KE=',F11.6,
     *       '  ROT KE=',F11.6,' KCAL/MOL')
 2001 FORMAT('----- QM PARTICLE COORDINATES FOR $DATA GROUP -----')
 2002 FORMAT('----- EFP PARTICLE COORDINATES FOR $EFRAG GROUP -----')
 2010 FORMAT(A8,A2,F5.1,3F20.10)
 2012 FORMAT(' $EFRAG'/
     *       'COORD=CART ','POSITION=OPTIMIZE')
 2013 FORMAT(' $END')
 2020 FORMAT('FRAGNAME=',A6,'  !',I4)
 2030 FORMAT(A8,7X,3F20.10)
C
 3000 FORMAT('----- RESTART VELOCITIES FOR $MD GROUP -----'/
     *       ' $MD READ=.TRUE. MBT=.FALSE. MBR=.FALSE.',
     *       ' TTOTAL=',1P,E12.2,0P)
 3001 FORMAT(' MDINT= ',A8,4x,'DT= ',E8.2,' NVTNH= ',I1,2x,'NSTEPS=',I8)
 3002 FORMAT(' RSTEMP=',A3,1x,'DTEMP= ',F8.2, ' LEVERY=',I6)
 3003 FORMAT(' RSRAND=',A3,1x,'NRAND=',I6,1x,'NVTOFF=',I3,1x,
     *       'JEVERY=',I6)
 3004 FORMAT(' PROD=',A3,3x,'KEVERY=',I6,1x,'DELR=',F8.3)
 3005 FORMAT('Batht(1)=',F6.2)
 3006 FORMAT(' NAMD=',A3,1X,'NDTLF=',I3,1X,'NDRST=',A3,1X,
     *        'NEWRND=',A3)
 3010 FORMAT(' SSBP=.T.   SFORCE= ',F5.1,' DROFF=',2F5.1)
 3011 FORMAT(' USAMP=.T.  UFORCE(1)=',3F5.1/
     *       10x,'  RZERO(1)= ',3F5.1)
 3012 FORMAT(' ipair(1)=',2I3)
 3013 FORMAT(' ipair(1)=',4I3)
 3014 FORMAT(' ipair(1)=',3I3)
 3015 FORMAT(' ipair(1)=',6I3)
 3016 FORMAT(' ipair(1)=',9I3)
 3017 FORMAT(' ipair(1)=',12I3)
 3018 FORMAT(' IUSTYP=',I3)
 3019 FORMAT(' CCMS=.T.   CFORCE= ',F5.1)
 3020 FORMAT('TVELQM(1)=',5X,'! QM ATOM TRANS. VELOCITIES (BOHR/PS) !')
 3030 FORMAT(2X,1P,E16.9,2X,E16.9,2X,E16.9)
 3040 FORMAT(2X,1P,E16.9,2X,E16.9,2X,E16.9,2X,E16.9)
 3050 FORMAT('TVEL(1)=',5X,'! EFP TRANSLATIONAL VELOCITIES (BOHR/PS) !')
 3060 FORMAT('QUAT(1)=',5X,'! EFP QUATERNIONS !')
 3070 FORMAT('RMOM(1)=',5X,'! EFP ANGULAR MOMENTA (RAD/PS) !')
 3075 FORMAT('RVEL(1)=',5X,'! EFP ANGULAR VELOCITY (RAD/PS) !')
 3080 FORMAT('QUAT1D(1)=',5X,'! EFP QUATERNION 1ST DERIV. !')
 3090 FORMAT('QUAT2D(1)=',5X,'! EFP QUATERNION 2ND DERIV. !')
 3095 FORMAT(' $END')
 3100 FORMAT(A8)
 3110 FORMAT(' $TDDFT'/
     *       ' NSTATE=',I2,' IROOT=',I2,' MULT=',I2,' TAMMD=',A3/
     *       ' $END')
      END SUBROUTINE SFPRTCORVEL
c*MODULE NAMD   *DECK NDPUSQLF
C>
C>     @brief Punch vectors in restarting format
C>
C>     @author Seunghoon Lee
C>
      SUBROUTINE NDPUSQLF(LUFILE,V,M,N)
C
      USE prec, ONLY: dp
      IMPLICIT NONE
C
      INTEGER :: LUFILE, M, N
      REAL(KIND=dp), DIMENSION(M*N) :: V
C
      INTEGER :: IC, IJ, IJ2, MAX, MIN, MODI, MODJ
C
C     ----- PUNCH A RECTANGULAR MATRIX WITH ORDERING LABELS -----
C     -V- IS -M- ROWS BY -N- COLUMNS
C
      MAX = 0
      IC = 0
      DO IJ = 1, M*N, 5
         MIN = MAX+1
         MAX = MAX+5
         IC = IC+1
         IF (MAX .GT. M*N) MAX = M*N
C
         MODI = MOD(MOD(IC,INT((M-1)/5)+1) ,100 )
         MODJ = MOD(INT((IJ-1)/M)+1  ,1000)
C
         WRITE (LUFILE,9008) MODJ,MODI,(V(IJ2),IJ2 = MIN,MAX)
      ENDDO
C
      RETURN
C
 9008 FORMAT(I2,I3,1P,5E15.8)
      END SUBROUTINE NDPUSQLF
C*MODULE NAMD   *DECK RREADND
C>
C>     @brief read vectors of restarting info
C>
C>     @author Seunghoon Lee
C>
      SUBROUTINE RREADND(IR,IW,KEY,N,A)
C
      USE prec, ONLY: dp
      USE comm_PAR, ONLY: MASWRK
      IMPLICIT NONE
C
      INTEGER :: IR, IW, N
      CHARACTER(*) :: KEY
      REAL(KIND=dp), DIMENSION(N) :: A
C
      INTEGER :: I, ITRY, IYES, K, L
      CHARACTER(80) :: LINE
C
      CALL VCLR(A,1,N)
C
      IYES=0
      ITRY=0
      L=LEN(KEY)
 50   DO WHILE(IYES.EQ.0)
         READ(IR,'(A80)',END=100,ERR=100) LINE
         IF (INDEX(LINE,KEY(1:L)).NE.0) THEN
            IYES=1
            EXIT
         ENDIF
      END DO
c
 100  ITRY=ITRY+1
      IF (IYES.EQ.0) THEN
         IF (ITRY.EQ.1) THEN
            REWIND(IR)
            GOTO 50
         ELSE
            CALL ABRTX('NO '//KEY(1:L)//' INPUT FOUND - STOP')
         ENDIF
      ELSE
         READ(IR,'(5X,1P,5E15.8)') (A(K),K=1,N)
         READ(IR,*,END=200)
      ENDIF
C
 200  RETURN
      END SUBROUTINE RREADND
c*MODULE NAMD   *deck rearrangeIDX
C>
C>     @brief rearrange indices
C>
C>     @author Yong Su Baek
C>
      SUBROUTINE rearrangeIDX(ipr,idxatom)
C
      USE mx_limits, only: mxatm
      USE prec, ONLY: dp
      IMPLICIT NONE
C
      INTEGER, PARAMETER :: NREPCA = 100
C
        REAL(KIND=dp), DIMENSION(nrEPCA) :: BATHT
        LOGICAL :: CCMS, PROD, RSRAND, RSTEMP, SSBP, USAMP
        REAL(KIND=dp) :: CFORCE, DTEMP, MDINT, SFORCE, TMIN, UV
        REAL(KIND=dp), DIMENSION(2) :: DROFF
        INTEGER, DIMENSION(50) :: INDEXOH, IPAIR
        INTEGER, DIMENSION(200) :: ISOLUT, MOVEMM
        INTEGER :: IUSTYP, JEVERY, KEVERY, LEVERY, MREMD, NCST, NEFPMV, &
     &             NRAND, NVTOFF
        REAL(KIND=dp), DIMENSION(20) :: PRVEC
        REAL(KIND=dp), DIMENSION(10) :: RZERO, UFORCE
        REAL(KIND=dp), DIMENSION(3) :: TMASS
      COMMON /CSTPOT/ DROFF, SFORCE, UFORCE, RZERO, DTEMP, TMIN, UV,    &
     &                BATHT, CFORCE, TMASS, PRVEC, MDINT, IPAIR,        &
     &                IUSTYP, NCST, NRAND, NVTOFF, JEVERY, KEVERY,      &
     &                LEVERY, MREMD, INDEXOH, MOVEMM, NEFPMV, RSTEMP,   &
     &                RSRAND, PROD, SSBP, CCMS, USAMP, ISOLUT
        REAL(KIND=dp), DIMENSION(3,MXATM) :: C
        INTEGER, DIMENSION(MXATM) :: IAN
        INTEGER :: ICH, MUL, NA, NAT, NB, NE, NQMT, NUM
        REAL(KIND=dp), DIMENSION(MXATM) :: ZAN
      COMMON /INFOA / NAT, ICH, MUL, NUM, NQMT, NE, NA, NB, ZAN, C, IAN
C
      INTEGER, DIMENSION(3) :: IDXATOM
      INTEGER, DIMENSION(*) :: IPR
C
      REAL(KIND=dp), DIMENSION(10) :: COMPDIS, DIS
      INTEGER :: I, IH, IO,  NUMH, NUMO
      INTEGER, DIMENSION(9) :: INDEXH
      INTEGER, DIMENSION(6) :: INDEXO, ISAVE
      REAL(KIND=dp), DIMENSION(3) :: UVEC1, UVEC2
      REAL(KIND=dp) :: UVECN1, VECN1, VECN_SAVE
c
      do i=1,6
        ISave(i)=0.0D+00
      enddo
      if (ipr(1) .eq. 3) then
      NumO =3
      NumH =4
      IndexO(1)=ipr(2)
      IndexO(2)=ipr(3)
      IndexO(3)=ipr(4)
      IndexH(1)=ipr(5)
      IndexH(2)=ipr(6)
      IndexH(3)=ipr(7)
      IndexH(4)=ipr(8)
      vecn_save = 1000000.0d+00
      do iO=1, NumO-1
        do iH=1, NumH
          uvecn1=0
          do i=1,3
             uvec1(i)=c(i,IndexO(iO))-c(i,IndexH(iH))
             uvecn1=uvecn1+uvec1(i)**2
          enddo
c
          vecn1=sqrt(uvecn1)
          if ( vecn1 .lt. vecn_save) then
            vecn_save = vecn1
            ISave(1)=IndexO(iO)
            ISave(2)=IndexH(iH)
            ISave(3)=IndexO(3)
          endif
        enddo
      enddo
c     rearrange ipair
         idxatom(1)=ISave(2)
         idxatom(2)=ISave(1)
         idxatom(3)=ISave(3)
      else if (ipr(1) .eq. 4) then
      NumO =2
      NumH =1
      IndexO(1)=ipr(2)
      IndexO(2)=ipr(3)
      IndexH(1)=ipr(4)
      vecn_save = 1000000.0d+00
      do iO=1, NumO
        do iH=1, NumH
          uvecn1=0
          do i=1,3
             uvec1(i)=c(i,IndexO(iO))-c(i,IndexH(iH))
             uvecn1=uvecn1+uvec1(i)**2
          enddo
c
          vecn1=sqrt(uvecn1)
          if ( vecn1 .lt. vecn_save) then
            vecn_save = vecn1
            ISave(1)=IndexO(iO)
            ISave(2)=IndexH(1)
          endif
        enddo
      enddo
c     rearrange ipair
         idxatom(1)=ISave(2)
         idxatom(2)=ISave(1)
      end if
c
      RETURN
      END SUBROUTINE rearrangeIDX
C*MODULE NAMD   *DECK ODP
C>
C>     @brief One-Dimensional Projection
C>            of Collective Variables technique
C>
C>     @author Yong Su Baek
C>
      SUBROUTINE ODP(ncst,prvec,f2f,ipair,uforce,rzero,uV)
C
      USE mx_limits, ONLY: mxatm
      USE constants, ONLY: half, pi
      USE prec, ONLY: dp
      IMPLICIT NONE
c
c     One-Dimensional Projection of Collective Variables technique.
c     IPAIR 1 X X X : Asymmetric coordinate requires three values in ipair
c     IPAIR 2 X X : Normal constraining between 2 atoms
c     IPAIR 3 X X X X X X X : Flexible asymmtric coordinate
c     IPAIR 4 X X X : Normal constraining for upto 2 bonds
c                     based on the minimal energy atomic pair(OH)
c     IPAIR 5 X X X : Angle constrain.
c
      REAL(KIND=dp), PARAMETER :: ANG2AU = 1.8897259877
C
        REAL(KIND=dp) :: AU2AMU, AU2ANG, AU2KCAL, AU2SEC, EK2BAR,       &
     &                   EK2KCAL
      COMMON /CONVMD/ AU2KCAL, EK2KCAL, AU2ANG, AU2SEC, AU2AMU, EK2BAR
        REAL(KIND=dp) :: E
        REAL(KIND=dp), DIMENSION(3*MXATM) :: EG
      COMMON /FUNCT / E, EG
        REAL(KIND=dp), DIMENSION(3,MXATM) :: C
        INTEGER, DIMENSION(MXATM) :: IAN
        INTEGER :: ICH, MUL, NA, NAT, NB, NE, NQMT, NUM
        REAL(KIND=dp), DIMENSION(MXATM) :: ZAN
      COMMON /INFOA / NAT, ICH, MUL, NUM, NQMT, NE, NA, NB, ZAN, C, IAN
        INTEGER :: IDAF, IJK, IPK, IP, IR, IW, NAV
        INTEGER, DIMENSION(950) :: IODA
      COMMON /IOFILE/ IR, IW, IP, IJK, IPK, IDAF, NAV, IODA
        LOGICAL :: DSKWRK, GOPARR, MASWRK
        INTEGER :: IBTYP, IPTIM, MASTER, ME, NPROC
      COMMON /PAR   / ME, MASTER, NPROC, IBTYP, IPTIM, GOPARR, DSKWRK,  &
     &                MASWRK
C
      REAL(KIND=dp) :: F2F, RZERO, UV
      INTEGER :: NCST
      INTEGER, DIMENSION(*) :: IPAIR
      REAL(KIND=dp), DIMENSION(*) :: PRVEC, UFORCE
C
      REAL(KIND=dp), DIMENSION(15) :: A
      REAL(KIND=dp) :: CBFACTOR, CBOFF, RDISP, TXDIS, UFACTOR, VECNORM, &
     &                XPVEC, XRDIS
      REAL(KIND=dp), DIMENSION(10) :: CTHETA, DOT, DPp, DS, TANVEC,      &
     &                                THETA, VECINT, XVEC
      INTEGER :: I, J, NPOS
      INTEGER, DIMENSION(30) :: IA
      INTEGER, DIMENSION(4) :: IDXATOM
      REAL(KIND=dp), DIMENSION(30) :: RDIS
      REAL(KIND=dp), DIMENSION(3,30) :: RIJ

c     initializing the variables
      vecnorm=0.0d+00
      XRdis=0.0d+00
      TXdis=0.0d+0
      do i=1,30
         do j=1,3
         rdis(i)=0.0d+00
         ia(i)=0.0d+00
         rij(j,i)=0.0d+00
         enddo
      enddo
      do i=1,4
         idxatom(i)=0.0d+00
      enddo
      do i=1,10
         tanvec(i)=0.0d+00
         Xvec(i)=0.0d+00
         DS(i)=0.0d+00
         Dpp(i)=0.0d+00
         dot(i)=0.0d+00
         ctheta(i)=0.0d+00
      enddo
      XPvec=0.0d+00
      uV=0.0d+00
c     tangent vector
      do i=1,ncst
         tanvec(i)=prvec(i)-prvec(i+ncst)
         vecnorm=vecnorm+tanvec(i)**2
      enddo
      vecnorm=sqrt(vecnorm)
      do i=1,ncst
         tanvec(i)=tanvec(i)/vecnorm
      enddo
c     Cartesian --> Internal Coordinates
c     Get atomic index for internal coordinates of Xvec
      npos=1
      do i=1,ncst
c     in case of asymmetric coordinates
         if (ipair(npos).eq.1) then
            idxatom(1)=ipair(npos+1)
            idxatom(2)=ipair(npos+2)
            idxatom(3)=ipair(npos+3)
            do j=1,3
               rij(j,npos)=c(j,idxatom(2))-c(j,idxatom(1))
               rij(j,npos+1)=c(j,idxatom(3))-c(j,idxatom(1))
               rdis(npos)=rdis(npos)+rij(j,npos)**2
               rdis(npos+1)=rdis(npos+1)+rij(j,npos+1)**2
            enddo
            rdis(npos)=sqrt(rdis(npos))
            rdis(npos+1)=sqrt(rdis(npos+1))
            Xvec(i)=(rdis(npos+1)-rdis(npos))*sqrt(half)
c     save the atom index
            do j=1,3
               rij(j,npos)=rij(j,npos)/rdis(npos)
               rij(j,npos+1)=rij(j,npos+1)/rdis(npos+1)
            enddo
            ia(npos+1)=3*(idxatom(1)-1)
            ia(npos+2)=3*(idxatom(2)-1)
            ia(npos+3)=3*(idxatom(3)-1)
            npos=npos+4
c     in case of bond coordinates
         else if (ipair(npos).eq.2) then
            idxatom(1)=ipair(npos+1)
            idxatom(2)=ipair(npos+2)
            do j=1,3
               rij(j,npos)=c(j,idxatom(2))-c(j,idxatom(1))
               rdis(npos)=rdis(npos)+rij(j,npos)**2
            enddo
            rdis(npos)=sqrt(rdis(npos))
            Xvec(i)=rdis(npos)
c     save the atom index
            do j=1,3
               rij(j,npos)=rij(j,npos)/rdis(npos)
            enddo
            ia(npos+1)=3*(idxatom(1)-1)
            ia(npos+2)=3*(idxatom(2)-1)
            npos=npos+3
c     in case of flexible asymmetric coordinates_4H
         else if (ipair(npos).eq.3) then
            call rearrangeIDX (ipair(npos), idxatom)
            do j=1,3
               rij(j,npos)=c(j,idxatom(2))-c(j,idxatom(1))
               rij(j,npos+1)=c(j,idxatom(3))-c(j,idxatom(1))
               rdis(npos)=rdis(npos)+rij(j,npos)**2
               rdis(npos+1)=rdis(npos+1)+rij(j,npos+1)**2
            enddo
            rdis(npos)=sqrt(rdis(npos))
            rdis(npos+1)=sqrt(rdis(npos+1))
            rdis(npos)=rdis(npos)*sqrt(half)
            rdis(npos+1)=rdis(npos+1)*sqrt(half)
            Xvec(i)=rdis(npos+1)-rdis(npos)
c     save the atom index
            do j=1,3
               rij(j,npos)=rij(j,npos)/rdis(npos)
               rij(j,npos+1)=rij(j,npos+1)/rdis(npos+1)
            enddo
            ia(npos+1)=3*(idxatom(1)-1)
            ia(npos+2)=3*(idxatom(2)-1)
            ia(npos+3)=3*(idxatom(3)-1)
            npos=npos+8
c     in case of flexible bond coordinates
         else if (ipair(npos).eq.4) then
            call rearrangeIDX (ipair(npos), idxatom)
            do j=1,3
               rij(j,npos)=c(j,idxatom(2))-c(j,idxatom(1))
               rdis(npos)=rdis(npos)+rij(j,npos)**2
            enddo
            rdis(npos)=sqrt(rdis(npos))
            Xvec(i)=rdis(npos)
c     save the atom index
            do j=1,3
               rij(j,npos)=rij(j,npos)/rdis(npos)
            enddo
            ia(npos+1)=3*(idxatom(1)-1)
            ia(npos+2)=3*(idxatom(2)-1)
            npos=npos+4
c      in case of angle
         else if (ipair(npos).eq.5) then
            idxatom(1)=ipair(npos+1)
            idxatom(2)=ipair(npos+2)
            idxatom(3)=ipair(npos+3)
            do j=1,3
               rij(j,npos)=c(j,idxatom(2))-c(j,idxatom(1))
               rij(j,npos+1)=c(j,idxatom(3))-c(j,idxatom(1))
               rdis(npos)=rdis(npos)+rij(j,npos)**2
               rdis(npos+1)=rdis(npos+1)+rij(j,npos+1)**2
               dot(i)=dot(i)+rij(j,npos)*rij(j,npos+1)
            enddo
            rdis(npos)=sqrt(rdis(npos))
            rdis(npos+1)=sqrt(rdis(npos+1))
            ctheta(i)=dot(i)/(rdis(npos)*rdis(npos+1))
            Xvec(i)=dacos(ctheta(i))
c     save the atom index
            do j=1,3
               rij(j,npos)=rij(j,npos)/rdis(npos)
               rij(j,npos+1)=rij(j,npos+1)/rdis(npos+1)
            enddo
            DS(i)=1.0d+00/dsin(Xvec(i))
            Dpp(i)=1.0d+00/(rdis(npos)*rdis(npos+1))
            ia(npos+1)=3*(idxatom(1)-1)
            ia(npos+2)=3*(idxatom(2)-1)
            ia(npos+3)=3*(idxatom(3)-1)
            npos=npos+4
         endif
      enddo
c     Get XR vector
      do i=1,ncst
         vecint(i)=Xvec(i)-prvec(i+ncst)
         XRdis=XRdis+vecint(i)**2
      enddo
      XRdis=sqrt(XRdis)
c     constraint to XR vector
c     Get XP vector Projected onto tanvec
      do i=1,ncst
        XPvec=XPvec+(vecint(i)*tanvec(i))
      enddo
c     Get TX vector
      TXdis=TXdis+sqrt((XRdis**2)-(XPvec**2))
c     applying Cylinder Boundary Force
      CBoff=0.0d+00
      if (TXdis.GT.CBoff) then
               A(1)=eg(ia(2)+1)/f2f
               A(2)=eg(ia(2)+2)/f2f
               A(3)=eg(ia(2)+3)/f2f
               A(4)=eg(ia(6)+1)/f2f
               A(5)=eg(ia(6)+2)/f2f
               A(6)=eg(ia(6)+3)/f2f
               A(7)=eg(ia(3)+1)/f2f
               A(8)=eg(ia(3)+2)/f2f
               A(9)=eg(ia(3)+3)/f2f
               A(10)=eg(ia(8)+1)/f2f
               A(11)=eg(ia(8)+2)/f2f
               A(12)=eg(ia(8)+3)/f2f
               A(13)=eg(ia(10)+1)/f2f
               A(14)=eg(ia(10)+2)/f2f
               A(15)=eg(ia(10)+3)/f2f
         uV=uV+uforce(2)*TXdis**2*half
         CBfactor=uforce(2)*f2f*TXdis
         npos=1
         do i=1,ncst
            if (ipair(npos).eq.1) then
               do j=1,3
               eg(ia(npos+1)+j)=eg(ia(npos+1)+j)-CBfactor
     *                          * ((vecint(i)-XPvec*tanvec(i))/TXdis)
     *                          * (rij(j,npos+1)-rij(j,npos))
               eg(ia(npos+2)+j)=eg(ia(npos+2)+j)-CBfactor
     *                          * ((vecint(i)-XPvec*tanvec(i))/TXdis)
     *                          * rij(j,npos)
               eg(ia(npos+3)+j)=eg(ia(npos+3)+j)+CBfactor
     *                          * ((vecint(i)-XPvec*tanvec(i))/TXdis)
     *                          * rij(j,npos+1)
               enddo
               npos=npos+4
            else if (ipair(npos).eq.2) then
               do j=1,3
               eg(ia(npos+1)+j)=eg(ia(npos+1)+j)-CBfactor
     *                          * ((vecint(i)-XPvec*tanvec(i))/TXdis)
     *                          * rij(j,npos)
               eg(ia(npos+2)+j)=eg(ia(npos+2)+j)+CBfactor
     *                          * ((vecint(i)-XPvec*tanvec(i))/TXdis)
     *                          * rij(j,npos)
               enddo
               npos=npos+3
            else if (ipair(npos).eq.3) then
               do j=1,3
               eg(ia(npos+1)+j)=eg(ia(npos+1)+j)-CBfactor
     *                          * ((vecint(i)-XPvec*tanvec(i))/TXdis)
     *                          * (rij(j,npos+1)-rij(j,npos))
               eg(ia(npos+2)+j)=eg(ia(npos+2)+j)-CBfactor
     *                          * ((vecint(i)-XPvec*tanvec(i))/TXdis)
     *                          * rij(j,npos)
               eg(ia(npos+3)+j)=eg(ia(npos+3)+j)-CBfactor
     *                          * ((vecint(i)-XPvec*tanvec(i))/TXdis)
     *                          * rij(j,npos+1)
               enddo
               npos=npos+8
            else if (ipair(npos).eq.4) then
               do j=1,3
               eg(ia(npos+1)+j)=eg(ia(npos+1)+j)-CBfactor
     *                          * ((vecint(i)-XPvec*tanvec(i))/TXdis)
     *                          * rij(j,npos)
               eg(ia(npos+2)+j)=eg(ia(npos+2)+j)+CBfactor
     *                          * ((vecint(i)-XPvec*tanvec(i))/TXdis)
     *                          * rij(j,npos)
               enddo
               npos=npos+4
            else if (ipair(npos).eq.5) then
               do j=1,3
               eg(ia(npos+1)+j)=eg(ia(npos+1)+j)+CBfactor
     *                         * ((vecint(i)-XPvec*tanvec(i))/TXdis)
     *                         *DS(i)*Dpp(i)*(rij(j,npos+1)+rij(j,npos))
               eg(ia(npos+2)+j)=eg(ia(npos+2)+j)-CBfactor
     *                         * ((vecint(i)-XPvec*tanvec(i))/TXdis)
     *                         *DS(i)*Dpp(i)*rij(j,npos+1)
               eg(ia(npos+3)+j)=eg(ia(npos+3)+j)-CBfactor
     *                         * ((vecint(i)-XPvec*tanvec(i))/TXdis)
     *                         *DS(i)*Dpp(i)*rij(j,npos)
               enddo
               npos=npos+4
            endif
         enddo
               A(1)=eg(ia(2)+1)/f2f
               A(2)=eg(ia(2)+2)/f2f
               A(3)=eg(ia(2)+3)/f2f
               A(4)=eg(ia(6)+1)/f2f
               A(5)=eg(ia(6)+2)/f2f
               A(6)=eg(ia(6)+3)/f2f
               A(7)=eg(ia(3)+1)/f2f
               A(8)=eg(ia(3)+2)/f2f
               A(9)=eg(ia(3)+3)/f2f
               A(10)=eg(ia(8)+1)/f2f
               A(11)=eg(ia(8)+2)/f2f
               A(12)=eg(ia(8)+3)/f2f
               A(13)=eg(ia(10)+1)/f2f
               A(14)=eg(ia(10)+2)/f2f
               A(15)=eg(ia(10)+3)/f2f
      endif
c     Divide the umbrella by rzero
      rdisp=XPvec-(rzero*vecnorm)
c     Calculate the uforce
      uV=uV+uforce(1)*rdisp**2*half
c     applying constraints
      ufactor=uforce(1)*f2f*rdisp
      npos=1
      do i=1,ncst
         if (ipair(npos).eq.1) then
      if (maswrk) write(iw,2005) ipair(npos+1),
     *               ipair(npos+2),ipair(npos+3)
      if (maswrk) write(iw,2002) Xvec(i)*au2ang
            npos=npos+4
         else if (ipair(npos).eq.2) then
      if (maswrk) write(iw,2006) ipair(npos+1),ipair(npos+2)
      if (maswrk) write(iw,2002) rdis(npos)*au2ang
            npos=npos+3
         else if (ipair(npos).eq.3) then
            call rearrangeIDX (ipair(npos), idxatom)
      if (maswrk) write(iw,2005) idxatom(1),
     *               idxatom(2),idxatom(3)
      if (maswrk) write(iw,2002) Xvec(i)*au2ang
            npos=npos+8
         else if (ipair(npos).eq.4) then
            call rearrangeIDX (ipair(npos), idxatom)
      if (maswrk) write(iw,2006) idxatom(1),idxatom(2)
      if (maswrk) write(iw,2002) rdis(npos)*au2ang
            npos=npos+4
         else if (ipair(npos).eq.5) then
      if (maswrk) write(iw,2012) ipair(npos+1),
     *               ipair(npos+2),ipair(npos+3)
      if (maswrk) write(iw,2001) Xvec(i)*180.0d+00/pi
            npos=npos+4
         endif
      enddo
      A(1)=eg(ia(2)+1)/f2f
      A(2)=eg(ia(2)+2)/f2f
      A(3)=eg(ia(2)+3)/f2f
      if (maswrk) write(iw,2003) A(1),A(2),A(3)
      npos=1
      do i=1,ncst
         if (ipair(npos).eq.1) then
            do j=1,3
            eg(ia(npos+1)+j)=eg(ia(npos+1)+j)-ufactor*tanvec(i)
     *                       * (rij(j,npos+1)-rij(j,npos))
            eg(ia(npos+2)+j)=eg(ia(npos+2)+j)
     *                       - ufactor*tanvec(i)*rij(j,npos)
            eg(ia(npos+3)+j)=eg(ia(npos+3)+j)
     *                       + ufactor*tanvec(i)*rij(j,npos+1)
            enddo
            npos=npos+4
         else if (ipair(npos).eq.2) then
            do j=1,3
            eg(ia(npos+1)+j)=eg(ia(npos+1)+j)
     *         -ufactor*tanvec(i)*rij(j,npos)
            eg(ia(npos+2)+j)=eg(ia(npos+2)+j)
     *         +ufactor*tanvec(i)*rij(j,npos)
            enddo
            npos=npos+3
         else if (ipair(npos).eq.3) then
            do j=1,3
            eg(ia(npos+1)+j)=eg(ia(npos+1)+j)-ufactor*tanvec(i)
     *                       * (rij(j,npos+1)-rij(j,npos))
            eg(ia(npos+2)+j)=eg(ia(npos+2)+j)
     *                       - ufactor*tanvec(i)*rij(j,npos)
            eg(ia(npos+3)+j)=eg(ia(npos+3)+j)
     *                       + ufactor*tanvec(i)*rij(j,npos+1)
            enddo
            npos=npos+8
         else if (ipair(npos).eq.4) then
            do j=1,3
            eg(ia(npos+1)+j)=eg(ia(npos+1)+j)
     *         -ufactor*tanvec(i)*rij(j,npos)
            eg(ia(npos+2)+j)=eg(ia(npos+2)+j)
     *         +ufactor*tanvec(i)*rij(j,npos)
            enddo
            npos=npos+4
         else if (ipair(npos).eq.5) then
            do j=1,3
            eg(ia(npos+1)+j)=eg(ia(npos+1)+j)+ufactor*tanvec(i)
     *                       *DS(i)*Dpp(i)*(rij(j,npos+1)+rij(j,npos))
            eg(ia(npos+2)+j)=eg(ia(npos+2)+j)-ufactor*tanvec(i)
     *                       *DS(i)*Dpp(i)*rij(j,npos+1)
            eg(ia(npos+3)+j)=eg(ia(npos+3)+j)-ufactor*tanvec(i)
     *                       *DS(i)*Dpp(i)*rij(j,npos)
            enddo
            npos=npos+4
         endif
      enddo
      if (maswrk) write(iw,2004) eg(ia(2)+1)/f2f,eg(ia(2)+2)/f2f,
     *   eg(ia(2)+3)/f2f
      if (maswrk) write(iw,2013) vecnorm
      if (maswrk) write(iw,2007) XPvec/vecnorm
      if (maswrk) write(iw,2008) TXdis
c
      RETURN
 2001 format(5x,'Current Angle (Deg.) = ',F15.10,/)
 2002 format(5x,'Current Distance (Ang.) = ',F15.10,/)
 2003 format(5x,'Gradient on Moving Atom = '/
     *      ,5x,3ES20.10)
 2004 format(5x,'Gradient on Moving Atom with Uforce = '/
     *      ,5x,3ES20.10,/)
 2005 format(5x,'Atoms (asymmetric constrain)',3i5)
 2012 format(5x,'Atoms (angle)',3i5)
 2006 format(5x,'Atoms (bond constrain)',2i5)
 2007 format(5x,'Current XPvec/prvec = ',F15.10)
 2008 format(5x,'Current TXdistance = ',F15.10)
 2009 format(5x,'Before CB force '/
     *      ,5x,3ES20.10)
 2010 format(5x,'After CB force '/
     *      ,5x,3ES20.10)
 2011 format(5x,'Current XPvec(Ang.) = ',F15.10)
 2013 format(5x,'vecnorm = ',F15.10)
      END SUBROUTINE ODP
C*MODULE NAMD   *DECK NAMDPRNT
C>
C>     @brief prints output for NAMD
C>
C>     @author Seunghoon Lee
C>
      SUBROUTINE NAMDPRNT(STAS,NACT,HOPPROB,DCHECK,TDECOE,TDECOEO,NDTLF,&
     &                    NDNINT,NDSR,RANDOM,NTHST)
C
      USE prec, ONLY: dp
      IMPLICIT NONE
C
      INTEGER :: IDAF, IJK, IJKT, IP, IR, IW, NAV
      INTEGER, DIMENSION(950) :: IODA
      COMMON /IOFILE/ IR, IW, IP, IJK, IJKT, IDAF, NAV, IODA
C
      REAL(KIND=dp) :: DCHECK, RANDOM
      INTEGER :: NDNINT, NDSR, NDTLF, NTHST
      REAL(KIND=dp), DIMENSION(NDSR,NDSR) :: HOPPROB, NACT, STAS
      REAL(KIND=dp), DIMENSION(2,NDSR) :: TDECOE, TDECOEO
C
      INTEGER :: I, IMAX, IMIN, J, MAX
C
C     OVERLAP INTEGRAL
      write (iw,2100) ndtlf
      MAX = 10
      IMAX = 0
c
      DO
         IMIN = IMAX+1
         IMAX = IMAX+MAX
         IF (IMAX .GT. ndsr) IMAX = ndsr
         WRITE (IW,2120) (I,       I=IMIN,IMAX)

         DO J = 1,ndsr
            WRITE (IW,2130) J,(stas(J,I),I = IMIN,IMAX)
         ENDDO
         IF (IMAX .LT. ndsr) THEN
            CYCLE
         ELSE
            EXIT
         END IF
      END DO

C     nonadiabatic coupling term
      write (iw,2200)
      MAX = 10
      IMAX = 0
C
      DO
         IMIN = IMAX+1
         IMAX = IMAX+MAX
         IF (IMAX .GT. ndsr) IMAX = ndsr
         WRITE (IW,2220) (I,       I=IMIN,IMAX)

         DO J = 1,ndsr
            WRITE (IW,2230) J,(nact(J,I),I = IMIN,IMAX)
         ENDDO
         IF (IMAX .LT. ndsr) THEN
            CYCLE
         ELSE
            EXIT
         END IF
      END DO

C     time-dependent coefficient
      write (iw,2300)
      write (iw,2310) ndnint
      write (iw,2320) abs(dcheck)
      write (iw,2330) abs(dcheck)
      write (iw,2340)
      do i = 1, ndsr
         write (iw,2350) i, tdecoeo(1,i)**2+tdecoeo(2,i)**2,
     *                      tdecoe(1,i)**2 +tdecoe(2,i)**2
      enddo

c     hopping probability
      write (iw,2400)
      write (iw,2410) random
      write (iw,2420)
      do i = 1, ndsr
         if (i.eq.nthst) cycle
         write (iw,2430) nthst, i, hopprob(nthst,i)
      enddo

 2100 FORMAT(/5X,40(1H-)/
     *        5X,'STATE OVERLAP INTEGRAL BETWEEN DIFFERENT'/
     *        5X,'   TIME STEPS BY USING TLF(',I1,') APPROX'/
     *        5X,'     (<PHI^{I}(T-DT)|PHI^{J}(T)>)'/
     *        5X,40(1H-))
 2120 FORMAT(5X,10(4X,I4,3X))
 2130 FORMAT(I5,10F11.6)
 2200 FORMAT(/5X,33(1H-)/
     *        5X,'NONADIABATIC COUPLING TERM (A.U.)'/
     *        5X,'BY USING FINITE DIFFERENCE APPROX'/
     *        5X,'    (<PHI^{I}|D/DT|PHI^{J}>)'/
     *        5X,33(1H-))
 2220 FORMAT(5X,10(4X,I4,3X))
 2230 FORMAT(I5,10F11.6)
 2300 FORMAT(/5X,45(1H-)/
     *        5X,'RESULT OF TIME-DEPENDENT SCHRODINGER EQUATION'/
     *        5X,' BY USING RUNGE-KUTTA 4TH ORDER METHOD (RK4)'/
     *        5X,45(1H-))
 2310 FORMAT( 5X,'THE NUMBER OF SUB-TIMESTEPS =',I10)
 2320 FORMAT( 5X,'LOCAL TRUNCATION ERROR OF RK4 = (',F7.3,')^5')
 2330 FORMAT( 5X,'TOTAL ACCUMULATE ERROR OF RK4 = (',F7.3,')^4')
 2340 FORMAT(/5X,'STATE #',5X,'|C(T-DT)|^2',4X,'|C(T)|^2  ')
 2350 FORMAT( 5X,I4,7X,F9.3,3X,'->',F9.3)
 2400 FORMAT(/5X,34(1H-)/
     *        5X,'RESULT OF SURFACE HOPPING BASED ON'/
     *        5X,'     FEWEST SWICHES ALGORITHM'/
     *        5X,34(1H-))
 2410 FORMAT( 5X,'RANDOM NUMBER AT THIS TIME STEP =',F10.5)
 2420 FORMAT(/5X,'STATE #', 2X ,'STATE #',2X,'HOPPING PROB')
 2430 FORMAT( 5X,    I4,3X,'->',    I4,3X,2X, F10.5)
      END SUBROUTINE NAMDPRNT
