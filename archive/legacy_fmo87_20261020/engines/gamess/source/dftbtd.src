C  6 Jun 18 - YN -  updates for FMO 5.3
C 18 Apr 16 - YN  - Implement TD-DFTB
C
C*MODULE DFTBTD    *DECK DFTB_TD_K
C>
C>    @brief Calculate K (coupling) matrix for TD-DFTB
C>
C>    @details Calculate K (coupling) matrix for TD-DFTB.
C>             Details of implementation related to TD-DFTB in GAMESS
C>             - TD-DFTB2 and TD-DFTB3
C>               Nishimoto, Y. J. Chem. Phys. 2015, 143, 091408.
C>             - TD-DFTB/PCM
C>               Nishimoto, Y. J. Phys. Chem. A 2016, 120, 771-784.
C>             - TD-LC-DFTB/PCM
C>               Nishimoto, Y. J. Phys. Chem. A 2019, 123, 5649-5659.
C>
C>    @author Yoshio Nishimoto
C>    - Feb, 2016- Subroutine written
C>    @data   Yoshio Nishimoto
C>    - Oct, 2019- Added TD-LC-DFTB
C>
C>           --- INPUT ---
C>    @param PA      BEC_AO in AO basis
C>    @param IMF     Calculate (A-B) matrix or not
C>    @param IPF     Calculate (A+B) matrix or not
C>    @param ITF     Calculate A matrix (ITF=1) or not
C>    @param NV      Number of vectors to be calculated
C>    @param L1      Number of basis functions
C>    @param TRIPLET Flag for singlet-triplet excitation
C>    @param TDER3   Flag for third-derivative of energy
C>           --- OUTPUT ---
C>    @param F       (A+B)*PA matrix in AO basis
C>    @param FM2     (A-B)*PA matrix in AO basis
C>    @param DERIV3  Third-derivative of energy, used in TD-DFTB3 grad
C>
C
      SUBROUTINE DFTB_TD_K(PA,F,FM2,DERIV3,IMF,IPF,ITF,NV,L1,TRIPLET,
     *                     TDER3)
      USE lrcdft, ONLY: LCFLAG, EMU, EMU2, LRFILE
      use mx_limits, only: mxatm,mxgtot,mxsh
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER (MXSPE=10)
C
      LOGICAL DFTBFL,SCC,SRSCC,DFTB3,DAMPXH,LCDFTB
C
      COMMON /DFTB  / DFTBFL,SCC,SRSCC,DFTB3,DAMPXH,LCDFTB
      COMMON /DFTBPR/ ETEMP,DFTBDP(MXSPE*14),DAMPXHE,HUBDER(MXSPE),
     *                ZREF(MXATM),SPNCST(6,MXSPE),SPE(MXATM),NSPE,
     *                MAXANG(MXATM),ISPE(MXATM),IND(MXATM+1),IDFTBD,
     *                PARAMDIR
      COMMON /DFTBSK/ SKDIM(MXSPE,MXSPE),SKSPIN(MXSPE),QREFL(3,MXSPE),
     *                HUBBL(3,MXSPE),QREF(MXSPE),HUBB(MXSPE),SKCUT2,
     *                NEEDSK,LSKHTAB(MXSPE,MXSPE),LSKSTAB(MXSPE,MXSPE),
     *                LSKGRID(MXSPE,MXSPE),LSKSELF(MXSPE)
      COMMON /FMCOM / X(1)
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON /NSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),
     *                CF(MXGTOT),CG(MXGTOT),CH(MXGTOT),CI(MXGTOT),
     *                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),
     *                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
C
      DIMENSION PA(L1,L1,NV)       ! X(K10)=BVEC_AO(L1,L1,NMAX) in AO basis
      DIMENSION F(L1*(L1+1)/2,NV)  ! X(K20)=B1PROA(L2,NMAX)    = A+B
      DIMENSION FM2(L1,L1,NV)      ! X(K30)=B2PROA(L1,L1,NMAX) = A-B
      DIMENSION DERIV3(L1*(L1+1)/2)! Third-order derivative of energy (DFTB3)
      LOGICAL   TRIPLET,TDER3
C
      L2 = L1*(L1+1)/2
      L3 = L1*L1
      IF (SRSCC) THEN !! .OR.SCFTYP.EQ.UHF) THEN
        NDIMGAM = NSHELL
      ELSE
        NDIMGAM = NAT
      END IF
C
C     ----- ALLOCATE MEMORY -----
C
      CALL VALFM(LOADFM)
      CALL GOTFM(NGOTMX)
      LWRK     = LOADFM  + 1
      LVEC     = LWRK    + MAX(L2,NAT*2)
      LS       = LVEC    + L3
      LGAMMA   = LS      + L2
      LDIST    = LGAMMA  + NDIMGAM*NDIMGAM
      LWRK1    = LDIST   + NAT*(NAT+1)/2
      LSHIFT   = LWRK1   + L3+L1
      LSHIFTS  = LSHIFT  + NAT
      LCHAMUL  = LSHIFTS + NSHELL
      LCHAMULS = LCHAMUL + NAT
      LAST     = LCHAMULS+ NSHELL
      IF (DFTB3) THEN
        LGAM3  = LAST
        LDQ    = LGAM3   + NDIMGAM*NDIMGAM
        LAST   = LDQ     + NAT
        IF (SRSCC) THEN
          LDQS = LAST
          LAST = LDQS    + NSHELL
        END IF
        IF (TDER3) THEN
          LSHIFT3 = LAST
          LAST    = LSHIFT3 + NDIMGAM
        END IF
      END IF
      IF (LCDFTB) THEN
        LGAMLC = LAST
        LAST   = LGAMLC  + NAT*NAT
      END IF
      NEED    = LAST - LOADFM - 1
C
C     ----- GET MEMORY -----
C
      CALL GETFM(NEED)
C
      CALL DAREAD(IDAF,IODA,X(LS),L2,12,0)
      CALL DAREAD(IDAF,IODA,X(LVEC),L3,15,0)
      IF (DFTB3) THEN
        !! GET MULLIKEN CHARGES (DQ=Q-Q0)
        CALL DAREAD(IDAF,IODA,X(LDQ),NAT,556,0)
        DO NI = 1, NAT
          X(LDQ+NI-1) = X(LDQ+NI-1) - ZREF(NI)
        END DO
        IF (SRSCC) THEN
          CALL DAREAD(IDAF,IODA,X(LDQS),NSHELL,557,0)
          ISH = 0
          DO I = 1, NAT
            DO J = 1, MAXANG(ISPE(I))
              ISH = ISH + 1
              X(LDQS+ISH-1) = X(LDQS+ISH-1) - QREFL(J,ISPE(I))
            END DO
          END DO
        END IF
      END IF
C
C     ----- CALCULATE DIPOLE INTEGRALS FOR OSCILLATOR STRENGTH ---
C     MAYBE DIFFERENT DEFINITION IS BETTER ?
C     USE <psi^i|r|\psi^j> = \sum_A R_A q_A^{ij} ?
C
C     IF(IST.EQ.1) THEN
C       CALL CALCOM(XP,YP,ZP)
C       CALL DFTB_DIPINT(IND,XP,YP,ZP,X(LS),.FALSE.)
C     END IF
C
C     ----- CALCULATE INTER-ATOMIC DISTANCES -----
C
      CALL VCLR(X(LDIST),1,NAT*(NAT+1)/2)
      NSEQ = 0
      DO I = 1, NAT
        DO J = 1, I
          NSEQ = NSEQ + 1
          DIST = SQRT( (C(1,I)-C(1,J))**2
     &               + (C(2,I)-C(2,J))**2
     &               + (C(3,I)-C(3,J))**2 )
          X(LDIST+NSEQ-1) = DIST
        END DO
      END DO
C
C     ----- CALCULATE GAMMA OF DFTB2 AND DFTB3 -----
C
      CALL DFTB_GAMMA(DAMPXHE,HUBBL,X(LWRK),X(LDIST),NAT,NSPE,
     *  ISPE,MAXANG,SPE,DAMPXH,SRSCC)
      CALL CPYTSQ(X(LWRK),X(LGAMMA),NDIMGAM,1)
      IF (DFTB3) THEN
        CALL VCLR(X(LGAM3),1,NDIMGAM*NDIMGAM)
        CALL DFTB_GAMMA3(NAT,HUBBL,X(LGAM3),X(LDIST))
      END IF
      IF (LCDFTB) THEN
        CALL DFTB_LCGAMMA(EMU,HUBB,NAT,ISPE,X(LDIST),X(LWRK))
        CALL CPYTSQ(X(LWRK),X(LGAMLC),NAT,1)
      END IF
C
C     ----- CALCULATE K (COUPLING MATRIX) -----
C
      CALL DFTB_TD_CPMAT(PA,X(LVEC),X(LS),X(LGAMMA),X(LGAM3),SPNCST,F,
     *                   FM2,DERIV3,X(LSHIFT),X(LSHIFTS),X(LCHAMUL),
     *                   X(LCHAMULS),X(LSHIFT3),X(LDQ),X(LDQS),
     *                   IND,L1,L2,NAT,NSHELL,NDIMGAM,NV,NSPE,ISPE,
     *                   MAXANG,IMF,IPF,ITF,TRIPLET,DFTB3,TDER3,
     *                   X(LWRK1),X(LGAMLC),LCDFTB,SRSCC)
C
C     ----- RELEASE ALLOCATED MEMORY -----
C
      CALL RETFM(NEED)
C
      RETURN
C
      END SUBROUTINE DFTB_TD_K
C
C-----------------------------------------------------------------------
C*MODULE DFTBTD    *DECK DFTB_TD_CPMAT
C>
C>    @brief Calculate K (coupling) matrix for TD-DFTB
C>
C>    @details Calculate K matrix in TD-DFTB calculations.
C>
C>    @author Yoshio Nishimoto
C>    - Feb, 2016- Subroutine written
C>    @data   Yoshio Nishimoto
C>    - Oct, 2019- Added TD-LC-DFTB
C>
C>           --- INPUT ---
C>    @param PA      BEC_AO in AO basis
C>    @param VEC     Temporary matrix to store vectors
C>    @param S       Overlap matrix
C>    @param GAMMA   Matrix of gamma used in DFTB2 and DFTB3
C>    @param GAMMA3  Matrix of gamma used in DFTB3
C>    @param SPNCST  Spin constants used for singlet-triplet excitation
C>    @param SHIFT   Tmpoerary matrix to store shift
C>    @param CHAMUL  Mulliken populations with PA charges
C>    @param CHAMULS Mulliken populations with PA charges
C>    @param SHIFT3  Temporary matrix to store shift for third deriv.
C>    @param DQ      Mulliken charges
C>    @param IND     AO index of each atom
C>    @param L1      Number of basis functions
C>    @param L2      L1*(L1+1)/2
C>    @param NAT     Number of atoms
C>    @param NSHELL  Number of shells
C>    @param NDIMGAM Number of atoms (RHF w/o SRSCC) or shells (else)
C>    @param NVEC    Number of vectors to be calculated
C>    @param NSPE    Number of species
C>    @param ISPE    Index of species
C>    @parma MAXANG  Maximum angular momentum of each species
C>    @param IMF     Calculate (A-B) matrix or not
C>    @param IPF     Calculate (A+B) matrix or not
C>    @param ITF     Calculate A matrix (ITF=1) or not
C>    @param TRIPLET Flag for singlet-triplet excitation
C>    @param DFTB3   Whether DFTB3 or not
C>    @param TDER3   Flag for third-derivative of energy
C>    @param WRK     Working array
C>    @param GAMLC   Gamma for long-range part
C>    @param LCDFTB  Enable long-range corrections
C>    @param SRSCC   Flag for shell-resolved SCC
C>           --- OUTPUT ---
C>    @param F       (A+B)*PA matrix in AO basis
C>    @param FM2     (A-B)*PA matrix in AO basis
C>    @param DERIV3  Third-derivative of energy, used in TD-DFTB3 grad.
C>
C
      SUBROUTINE DFTB_TD_CPMAT(PA,VEC,S,GAMMA,GAMMA3,SPNCST,F,FM2,
     *                         DERIV3,SHIFT,SHIFTS,CHAMUL,CHAMULS,
     *                         SHIFT3,DQ,DQS,
     *                         IND,L1,L2,NAT,NSHELL,NDIMGAM,NVEC,NSPE,
     *                         ISPE,MAXANG,IMF,IPF,ITF,TRIPLET,DFTB3,
     *                         TDER3,WRK,GAMLC,LCDFTB,SRSCC)
C
      IMPLICIT NONE
C
      DOUBLE PRECISION, INTENT(IN) :: PA(L1*L1,NVEC),S(*),
     *  GAMMA(NDIMGAM,NDIMGAM),GAMMA3(NDIMGAM,NDIMGAM),SPNCST(6,NSPE),
     *  DQ(NAT),DQS(NSHELL),GAMLC(*)
      DOUBLE PRECISION, INTENT(INOUT) :: VEC(L1,L1),F(L2,NVEC),
     *  FM2(L1*L1,NVEC),DERIV3(L2),SHIFT(NAT),SHIFTS(NSHELL),
     *  CHAMUL(NAT),CHAMULS(NSHELL),SHIFT3(NAT),WRK(L1,L1+1)
      INTEGER, INTENT(IN) :: IND(*),L1,L2,NAT,NSHELL,NDIMGAM,NVEC,NSPE,
     *  ISPE(*),MAXANG(*),IMF,IPF,ITF
      LOGICAL, INTENT(IN) :: TRIPLET,DFTB3,TDER3,LCDFTB,SRSCC
C
      INTEGER :: IVEC,NI,NJ,NC,MM,MU,NU,II,JJ,KK,ISH,JSH,KSH,
     *           ISH0,JSH0,KSH0
      DOUBLE PRECISION :: TMP,VAL,DDOT,XDOT
      DOUBLE PRECISION, PARAMETER :: ZERO=0.0D+00,TWO=2.0D+00,
     *  ONE_THIRD=1.0D+00/3.0D+00,TWO_THIRD=TWO*ONE_THIRD
C
      DO IVEC = 1, NVEC
C
C       ----- CALCULATE q_C^{jb} FIRST -----
C
        CALL SYMSQT(L1,PA(1,IVEC),L1,VEC,1)
        CALL DFTB_MULLIKEN(L1,L2,VEC,S,WRK)
        CALL DFTB_MULSA(WRK,CHAMULS,CHAMUL,L1,NSHELL,NAT)
        CALL DSCAL(NAT,2.0D+00,CHAMUL,1)
        IF (SRSCC) CALL DSCAL(NSHELL,2.0D+00,CHAMULS,1)
C
C       ----- CALCULATE SUM OF SHIFT CONTRIBUTION ON EACH ATOM -----
C       IMF = 1 : CALCULATE (A-B) MATRIX
C       IPF = 1 : CALCULATE (A+B) MATRIX
C       ITF = 1 : ONLY CALCULATE A MATRIX
C
C       IF NOT LC-DFTB,
C       A-B MATRIX IS ALWAYS DIAGONAL, MEANING THAT A-B IS ZERO MATRIX
C       WHEN REFERENCE ELECTRONIC STRUCTURE IS A CLOSED SHELL SINGLET.
C
C       CALCULATE A+B MATRIX
        IF (TRIPLET) THEN
          IF (SRSCC) THEN
            CALL VCLR(SHIFTS,1,NSHELL)
            CALL DFTB_SPIN_SHIFT(SHIFTS,CHAMULS,SPNCST,NAT,NSHELL,
     *                           NSPE,ISPE,MAXANG)
          ELSE
            CALL VCLR(SHIFT,1,NAT)
            DO NI = 1, NAT
              VAL = CHAMUL(NI)
              SHIFT(NI) = SHIFT(NI) + VAL*SPNCST(1,ISPE(NI))
            END DO
          END IF
        ELSE
          IF (SRSCC) THEN
            CALL DGEMV('T',NSHELL,NSHELL,
     *                 1.0D+00,GAMMA,NSHELL,CHAMULS,1,
     *                 0.0D+00,SHIFTS,1)
            IF (DFTB3) THEN
              CALL DFTB_3RD_SHIFT2(SHIFTS,DQ,DQS,CHAMUL,CHAMULS,
     *                             GAMMA3,NAT,NSHELL,NDIMGAM,MAXANG,
     *                             ISPE,SRSCC)
            END IF
          ELSE
            CALL DGEMV('T',NAT,NAT,
     *                 1.0D+00,GAMMA,NAT,CHAMUL,1,
     *                 0.0D+00,SHIFT,1)
            IF (DFTB3) THEN
              CALL DFTB_3RD_SHIFT2(SHIFT,DQ,DQS,CHAMUL,CHAMULS,
     *                             GAMMA3,NAT,NSHELL,NDIMGAM,MAXANG,
     *                             ISPE,SRSCC)
C             DO NI = 1, NAT
C               TMP = ZERO
C               DO NC = 1, NAT
C                 TMP = TMP + GAMMA3(NI,NC)*CHAMUL(NC)*DQ(NI)
C    *                      + GAMMA3(NI,NC)*CHAMUL(NI)*DQ(NC)
C    *                      + GAMMA3(NC,NI)*CHAMUL(NC)*DQ(NC)
C               END DO
C               SHIFT(NI) = SHIFT(NI) + TMP*TWO_THIRD
C             END DO
            END IF
          END IF
        END IF
        IF (ITF.NE.1) THEN
          IF (SRSCC) THEN
            CALL DSCAL(NSHELL,2.0D+00,SHIFTS,1) !! total coeff.=4
          ELSE
            CALL DSCAL(NAT,2.0D+00,SHIFT,1) !! total coeff.=4
          END IF
        END IF
        IF (TDER3) THEN
          CALL VCLR(SHIFT3,1,NDIMGAM)
          CALL DFTB_3RD_SHIFT2(SHIFT3,CHAMUL,CHAMULS,CHAMUL,CHAMULS,
     *                         GAMMA3,NAT,NSHELL,NDIMGAM,MAXANG,
     *                         ISPE,SRSCC)
C         DO NI = 1, NAT
C           TMP = ZERO
C           DO NC = 1, NAT
C             TMP = TMP + GAMMA3(NI,NC)*CHAMUL(NC)*CHAMUL(NI)
C    *                  + GAMMA3(NI,NC)*CHAMUL(NI)*CHAMUL(NC)
C    *                  + GAMMA3(NC,NI)*CHAMUL(NC)*CHAMUL(NC)
C           END DO
C           SHIFT3(NI) = SHIFT3(NI) + TMP*TWO_THIRD !! total coeff=8/3
C           !! note that DERIV3 will be multiplied by two outside
C           !! the subroutine
C         END DO
        END IF
C
C       ----- ADD SHIFT CONTRIBUTION TO F MATRIX -----
C
        CALL VCLR(WRK,1,L2)
        IF (SRSCC) THEN
          CALL DFTB_SHIFT_FOCK_SRSCC(WRK,SHIFTS,S,NAT,L1,L2,NSHELL,NSPE,
     *                         ISPE,MAXANG)
        ELSE
          CALL DFTB_SHIFT_FOCK(SHIFT,WRK,S,NAT,L2,NSPE,ISPE,IND,
     *                         MAXANG)
        END IF
        IF (ITF.EQ.1) CALL CPYTSQ(WRK,FM2(1,IVEC),L1,1)
        IF (ITF.EQ.0 .AND. IPF.EQ.1) CALL DCOPY(L2,WRK,1,F(1,IVEC),1)
C
C       ----- ADD EXCHANGE CONTRIBUTIONS -----
C
        IF (LCDFTB) THEN
          !! Compute exchange contributions with asymmetric pseudo-DM
          CALL VCLR(WRK,1,L1*L1)
          CALL DFTB_LCSHIFT3(NAT,L1,ISPE,MAXANG,IND,PA(1,IVEC),
     *                       S,GAMLC,WRK)
          IF (ITF.EQ.0) THEN
            CALL DSCAL(L1*L1,4.0D+00,WRK,1)
            IF (IPF.EQ.1) THEN
              !! For (A+B), symmetrize
              CALL SYMSQT(L1,WRK,L1,VEC,1)
              CALL DAXPY(L2,1.0D+00,VEC,1,F(1,IVEC),1)
            END IF
C
            IF (IMF.EQ.1) THEN
              !! For (A-B), anti-symmetrize
              CALL DCOPY(L1*L1,WRK,1,FM2(1,IVEC),1)
              CALL ASYMTRZE(FM2(1,IVEC),L1,L1,TMP)
            END IF
          ELSE IF (ITF.EQ.1) THEN
            !! This is not necessarily symmetric
            CALL DAXPY(L1*L1,2.0D+00,WRK,1,FM2(1,IVEC),1)
          END IF
        END IF
      END DO
C
      IF (TDER3) THEN
        CALL VCLR(DERIV3,1,L2)
        IF (.NOT.TRIPLET) THEN
          IF (SRSCC) THEN
            CALL DFTB_SHIFT_FOCK_SRSCC(DERIV3,SHIFT3,S,NAT,L1,L2,NSHELL,
     *                           NSPE,ISPE,MAXANG)
          ELSE
            CALL DFTB_SHIFT_FOCK(SHIFT3,DERIV3,S,NAT,L2,NSPE,ISPE,IND,
     *                           MAXANG)
          END IF
        END IF
      END IF
C
      RETURN
C
      END SUBROUTINE DFTB_TD_CPMAT
C
C-----------------------------------------------------------------------
C*MODULE DFTBTD    *DECK DFTB_TD_GRAD_SHIFT
C>
C>    @brief Gradient for TD-DFTB
C>
C>    @details Calculate DQES and DQXY (Mulliken charges with excited
C>             density matrix (P=T+Z) and (X+Y) matrix in AO basis) and
C>             construct shift matrices with those charges.
C>
C>    @author Yoshio Nishimoto
C>    - Feb, 2016- Subroutine written
C>
C>           --- INPUT ---
C>    @param WRK     Working array
C>    @param S       Overlap matrix
C>    @param DAMPXHE Exponent of X-H damping
C>    @param DGS     Ground state density matrix
C>    @param DQ      Ground state Mulliken charges
C>    @param DIST    Distance matrix
C>    @param HUBBL   Shell-resolved Hubbar values
C>    @param SPE     Name of species
C>    @param ICPM    PCM or not
C>    @param NSPE    Number of species
C>    @param NAT     Number of atoms
C>    @param L1      Number of basis functions
C>    @param L2      L1*(L1+1)/2
C>    @param IDAF    10
C>    @param IODA    Length of unit
C>    @param ISPE    Index of species
C>    @parma MAXANG  Maximum angular momentum of each species
C>    @param IND     AO index of each atom
C>    @param DFTB3   Whether DFTB3 or not
C>    @param DAMPXH  Whether X-H damping is used or not
C>    @parma SRSCC   Whether Shell-resolved SCC or not
C>    @parma TRIPLET Whether triplet or not
C>           --- OUTPUT ---
C>    @param XYAO    (X+Y)_{\mu \nu}
C>    @param DQES    Mulliken population with excited P (=T+Z)
C>    @param DQXY    Mulliken population with (X+Y)
C>    @param SHIFTES Shift matrix with DQES
C>    @param SHIFTXY Shift matrix with DQXY
C>
C
      SUBROUTINE DFTB_TD_GRAD_SHIFT(WRK,S,DAMPXHE,DGS,DQ,DQS,DIST,HUBBL,
     *                              SPE,SPNCST,XYAO,DQES,DQESS,DQXY,
     *                              DQXYS,SHIFTES,SHIFTXY,GAMMA,GAMMA3,
     *                              IPCM,NSPE,NAT,L1,L2,IDAF,IODA,ISPE,
     *                              MAXANG,IND,NDIMGAM,NSHELL,
     *                              DFTB3,DAMPXH,SRSCC,TRIPLET)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DOUBLE PRECISION, INTENT(INOUT) :: WRK(L1,L1),S(L2)
      DOUBLE PRECISION, INTENT(IN)    :: DAMPXHE,DGS(L2),DQ(NAT),
     *                                   DQS(NSHELL),DIST(*),HUBBL(*),
     *                                   SPE(*),SPNCST(6,NSPE)
      DOUBLE PRECISION, INTENT(OUT)   :: XYAO(L2),DQES(NAT),
     *                                   DQESS(NSHELL),DQXY(NAT),
     *                                   DQXYS(NSHELL),SHIFTES(NAT),
     *                                   SHIFTXY(NAT),GAMMA(*),GAMMA3(*)
      INTEGER, INTENT(IN) :: IPCM,NSPE,NAT,L1,L2,IDAF,IODA(950),
     *                       ISPE(NAT),MAXANG(NSPE),IND(NAT+1)
      LOGICAL, INTENT(IN) :: DFTB3,DAMPXH,SRSCC,TRIPLET
C
      DOUBLE PRECISION, PARAMETER :: ZERO=0.0D+00
C
      COMMON /FFPARM/ NFFAT,NBOND,NANGL,NDIHR,NDIHB,NCMAP,NWAGG,
     *                N1213J,N14J,NLKQMM,IDOCHG,IDOPOL,IDOLJ,IDOCMAP
C
      L3 = L1*L1
C
C     ----- MAKE q^ES -----
C
      CALL DAREAD(IDAF,IODA,WRK,L3,476,0) !! PRISTINE P MATRIX
      CALL VCLR(DQES,1,NAT)
      IF (SRSCC) CALL VCLR(DQESS,1,NSHELL)
      ISH = 0
      DO NI = 1, NAT
        TMP = ZERO
        DO MM = 1, IND(NI+1)-IND(NI)
          MU = IND(NI)+MM
          DO NU = 1, L1
            CALL DFTB_CNVSQ(MU,NU,NSEQ)
            TMP = TMP + (WRK(MU,NU)+WRK(NU,MU))*S(NSEQ)
          END DO
        END DO
        DQES(NI) = TMP
C
        IF (SRSCC) THEN
          DO MM = 1, IND(NI+1)-IND(NI)
            IF (MM.EQ.1.OR.MM.EQ.2.OR.MM.EQ.5) THEN
              TMP = 0.0D+00
              ISH = ISH + 1
            END IF
            MU = IND(NI)+MM
            DO NU = 1, L1
              CALL DFTB_CNVSQ(MU,NU,NSEQ)
              TMP = TMP + (WRK(MU,NU)+WRK(NU,MU))*S(NSEQ)
            END DO
            IF (MM.EQ.1.OR.MM.EQ.4.OR.MM.EQ.9) DQESS(ISH) = TMP
          END DO
        END IF
      END DO
C
C     ----- MAKE q^XY -----
C
      CALL DAREAD(IDAF,IODA,WRK,L3,474,0) !! (X+Y) in AO, Forget LP
      CALL VCLR(DQXY,1,NAT)
      IF (SRSCC) CALL VCLR(DQXYS,1,NSHELL)
      ISH = 0
      DO NI = 1, NAT
        TMP = ZERO
        DO MM = 1, IND(NI+1)-IND(NI)
          MU = IND(NI)+MM
          DO NU = 1, L1
            CALL DFTB_CNVSQ(MU,NU,NSEQ)
            TMP = TMP + (WRK(MU,NU)+WRK(NU,MU))*S(NSEQ)
          END DO
        END DO
        DQXY(NI) = TMP
        IF (SRSCC) THEN
          DO MM = 1, IND(NI+1)-IND(NI)
            IF (MM.EQ.1.OR.MM.EQ.2.OR.MM.EQ.5) THEN
              TMP = 0.0D+00
              ISH = ISH + 1
            END IF
            MU = IND(NI)+MM
            DO NU = 1, L1
              CALL DFTB_CNVSQ(MU,NU,NSEQ)
              TMP = TMP + (WRK(MU,NU)+WRK(NU,MU))*S(NSEQ)
            END DO
            IF (MM.EQ.1.OR.MM.EQ.4.OR.MM.EQ.9) DQXYS(ISH) = TMP
          END DO
        END IF
      END DO
C
C     ----- CONSTRUCT SYMMETRIC (X+Y) IN AO BASIS -----
C
      CALL SYMTRZE(WRK,L1,L1) !! multiplied with 0.5 in the subroutine
      NSEQ = 0
      DO I = 1, L1
        DO J = 1, I
          NSEQ = NSEQ + 1
          XYAO(NSEQ) = WRK(J,I)
        END DO
      END DO
C
C     ----- DO SHIFT WITH XY AND ES CHARGES -----
C
      CALL DFTB_GAMMA(DAMPXHE,HUBBL,GAMMA,DIST,NAT,NSPE,
     *                ISPE,MAXANG,SPE,DAMPXH,SRSCC)
      IF (DFTB3) CALL DFTB_GAMMA3(NAT,HUBBL,GAMMA3,DIST)
      !! Shift by excited state charge
      CALL VCLR(SHIFTES,1,NDIMGAM)
      IF (SRSCC) THEN
        CALL DFTB_2ND_SHIFT(SHIFTES,DQESS,GAMMA,NDIMGAM)
      ELSE
        CALL DFTB_2ND_SHIFT(SHIFTES,DQES,GAMMA,NDIMGAM)
      END IF
      !! Shift by (X+Y) vector
      CALL VCLR(SHIFTXY,1,NDIMGAM)
      IF (TRIPLET) THEN
        IF (SRSCC) THEN
          CALL DFTB_SPIN_SHIFT(SHIFTXY,DQXYS,SPNCST,NAT,NSHELL,
     *                         NSPE,ISPE,MAXANG)
        ELSE
          DO NI = 1, NAT
            SHIFTXY(NI) = SHIFTXY(NI) + DQXY(NI)*SPNCST(1,ISPE(NI))
          END DO
        END IF
      ELSE
        IF (SRSCC) THEN
          CALL DFTB_2ND_SHIFT(SHIFTXY,DQXYS,GAMMA,NDIMGAM)
        ELSE
          CALL DFTB_2ND_SHIFT(SHIFTXY,DQXY,GAMMA,NDIMGAM)
        END IF
      END IF
      !! DFTB3 contributions
      IF (DFTB3) THEN
        IF (SRSCC) THEN
          CALL DFTB_3RD_SHIFT2(SHIFTES,DQ,DQS,DQES,DQESS,GAMMA3,
     *                         NAT,NSHELL,NDIMGAM,MAXANG,ISPE,SRSCC)
          IF (.NOT.TRIPLET) THEN
            CALL DFTB_3RD_SHIFT2(SHIFTXY,DQ,DQS,DQXY,DQXYS,GAMMA3,
     *                           NAT,NSHELL,NDIMGAM,MAXANG,ISPE,SRSCC)
            CALL DFTB_3RD_SHIFT2(SHIFTES,DQXY,DQXYS,DQXY,DQXYS,GAMMA3,
     *                           NAT,NSHELL,NDIMGAM,MAXANG,ISPE,SRSCC)
          END IF
        ELSE
          CALL DFTB_TD_3RD_SHIFT(SHIFTES,SHIFTXY,DQES,DQXY,DQ,
     *                                    GAMMA3,NAT,TRIPLET)
        END IF
      END IF
C
C     ----- PCM CORRECTION FOR SHIFT MATRICES -----
C
      IF (IPCM.EQ.1) THEN
        CALL VCLR(S,1,L2)
        CALL DENTD(S,WRK,L1) !! S = SYMMETRIC 2P MATRIX
        CALL DFTB_TD_GRAD_PCM(DGS,S,XYAO,SHIFTES,SHIFTXY,
     *                        MAXANG,ISPE,IND,L1,L2)
      END IF
      IF(NFFAT.GT.0) THEN
        CALL VCLR(S,1,L2)
        CALL DENTD(S,WRK,L1) !! S = SYMMETRIC 2P MATRIX
        CALL DFTB_TD_GRAD_POL(S,XYAO,SHIFTES,SHIFTXY,L1,L2)
      END IF
C
      RETURN
C
      END SUBROUTINE DFTB_TD_GRAD_SHIFT
C
C-----------------------------------------------------------------------
C*MODULE DFTBTD    *DECK DFTB_TD_3RD_SHIFT
C>
C>    @brief Gradient for TD-DFTB
C>
C>    @details Add shift contribution of DQES and DQXY for TD-DFTB3
C>
C>    @author Yoshio Nishimoto
C>    - Feb, 2016- Subroutine written
C>
C>           --- INPUT ---
C>    @param DQES    Mulliken population with excited P (=T+Z)
C>    @param DQXY    Mulliken population with (X+Y)
C>    @param DQ      Ground state Mulliken charges
C>    @param GAMMA   Matrix of gamma used in DFTB3
C>    @param NDIM    Dimension of gamma (NAT)
C>           --- OUTPUT ---
C>    @param SHIFTES Shift matrix with DQES
C>    @param SHIFTXY Shift matrix with DQXY
C>
C
      SUBROUTINE DFTB_TD_3RD_SHIFT(SHIFTES,SHIFTXY,DQES,DQXY,DQ,GAMMA,
     *                             NDIM,TRIPLET)
C
      IMPLICIT NONE
C
      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(IN) :: DQES(NDIM),DQXY(NDIM),DQ(NDIM),
     *                                GAMMA(NDIM,NDIM)
      DOUBLE PRECISION, INTENT(INOUT) :: SHIFTES(NDIM),SHIFTXY(NDIM)
      LOGICAL, INTENT(IN) :: TRIPLET
C
      INTEGER :: I,J
      DOUBLE PRECISION :: GAMI,GAMJ
      DOUBLE PRECISION, PARAMETER :: TWO_THIRD=2.0D+00/3.0D+00
C
C     DQES terms come from the derivative of Fock matrix
C     DQXY terms come from the derivative of coupling matrix (4K)
C     DQES and DQXY are already twiced
C
      DO I = 1, NDIM !! a
        DO J = 1, NDIM !! c
          GAMI = GAMMA(I,J) !! G_ac
          GAMJ = GAMMA(J,I) !! G_ca
          !! a) mulliken derivative term (d)
          !! b) Dmn term (b) (same to mulliken derivative)
          !! c) mulliken derivative for (X+Y)4K(X+Y)
          SHIFTES(I) = SHIFTES(I)
     *      + TWO_THIRD*(DQ(I)*DQES(J)*GAMI
     *      +            DQ(J)*DQES(I)*GAMI
     *      +            DQ(J)*DQES(J)*GAMJ)*0.5D+00 !! a)
     *      + TWO_THIRD*(DQ(J)*DQES(J)*GAMJ
     *      +            DQ(I)*DQES(J)*GAMI
     *      +            DQ(J)*DQES(I)*GAMI)*0.5D+00 !! b)
          IF (.NOT.TRIPLET) THEN
            SHIFTES(I) = SHIFTES(I)
     *        + TWO_THIRD*(DQXY(I)*DQXY(J)*GAMI
     *        +            DQXY(I)*DQXY(J)*GAMI
     *        +            DQXY(J)*DQXY(J)*GAMJ) !! c)
            SHIFTXY(I) = SHIFTXY(I)
     *        + TWO_THIRD*(DQ(I)*DQXY(J)*GAMI
     *        +            DQ(J)*DQXY(I)*GAMI
     *        +            DQ(J)*DQXY(J)*GAMJ)
          END IF
        END DO
      END DO
C
      RETURN
C
      END SUBROUTINE DFTB_TD_3RD_SHIFT
C
C-----------------------------------------------------------------------
C*MODULE DFTBTD    *DECK DFTB_TD_GRAD_PCM
C>
C>    @brief Gradient for TD-DFTB
C>
C>    @details Add PCM contribution to shift matrices
C>
C>    @author Yoshio Nishimoto
C>    - Feb, 2016- Subroutine written
C>
C>           --- INPUT ---
C>    @param DM      Ground state density matrix
C>    @param DMES    P(=T+Z)
C>    @param DMXY    (X+Y) in AO basis
C>    @parma MAXANG  Maximum angular momentum of each species
C>    @param ISPE    Index of species
C>    @param IND     AO index of each atom
C>    @param L1      Number of basis functions
C>    @param L2      L1*(L1+1)/2
C>           --- OUTPUT ---
C>    @param SHIFTES Shift matrix with DQES
C>    @param SHIFTXY Shift matrix with DQXY
C>
C
      SUBROUTINE DFTB_TD_GRAD_PCM(DM,DMES,DMXY,SHIFTES,SHIFTXY,MAXANG,
     *                            ISPE,IND,L1,L2)
      use mx_limits, only: mxatm
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      LOGICAL DFTBFL,SCC,SRSCC,DFTB3,DAMPXH,LCDFTB
C
      PARAMETER (HALF=0.5D+00)
C
      DOUBLE PRECISION, INTENT(IN) :: DM(L2),DMES(*),DMXY(*)
      DOUBLE PRECISION, INTENT(INOUT) :: SHIFTES(NAT),SHIFTXY(NAT)
      INTEGER, INTENT(IN) :: L1,L2,MAXANG(*),ISPE(*),IND(*)
      LOGICAL TRIPLET,SG1T,NONEQR,TAMMD,TPA,ALPHKWD,BETAKWD
      LOGICAL MREKT,MRDEA
C
      COMMON /DFTB  / DFTBFL,SCC,SRSCC,DFTB3,DAMPXH,LCDFTB
      COMMON /DFTBS / LS
      COMMON /FMCOM / X(1)
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /INFOTD/ CNVTOL,PFREQ(2),MODTD,
     *                JANST,NRADT,NTHET,NPHIT,NLEBT,
     *                NSTAT,NTRIAL,MAXVEC,NTHST,IRECTD,ITDFG,ITDPRP,
     *                TRIPLET,SG1T,NONEQR,TAMMD,TPA,ALPHKWD,BETAKWD,
     *                SPCP(3),MULTD,MREKT,MRDEA,MTHST,IFEDAT(4)
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON /PCMDIM/ MXSP,MXTS,MEMPCM1,MEMPCM2,NTS
      COMMON /PCMPAR/ IPCM,NFT26,NFT27,IRPPCM,IEF,IP_F,NFMOPCM,IHET
      COMMON /PCMPNT/ LXYZRE,LSSFE,LLIST,LALPCM,LRINPCM,LINAPCM,LINFPCM,
     *                LINIPCM,LMEPCM,LAXYZCT,LXYZCT2,LNVERT,LQSN,LQSE,
     *                LQFS,LQIND,LISPHE,LVAD,LQOR,LVECMUL,LAIPRJ,LFIPRJ,
     *                LPEL,LPCMCDR,LDAI,LIDDAI,LQSND,LQSED,LEPSPCM
C
      LOGICAL IPCFP
C
      ITDFG = 1
      L3 = L1*L1
C
      CALL VALFM(LOADFM)
      LFAO = LOADFM + 1
      LS   = LFAO   + L2
      LAST = LS     + L2
C
C     -- PCM MEMORY
C
      ISD=1
      IPCFP=.FALSE.
      CALL PCMMEM3(IPCFP,L2,LAST,NTSPAR,LSOL,LCHG,LDMATM,LQPOT,LQFLD,
     *             LVPOT,LSE,LDE,LQ0,LQ1,LQ2,LQ3,LD0,LQA,LDIMAT,LQREP,
     *             LTMP0,LTMP1,LTMP2,LIPVT,LPOTTMP,LRMUL,LCHG2,
     *             LASCCRD,LASCCHG,LASCDIP,LASCQAD,LDISV,LDIS1,LDIS2,
     *             LCQEF,LQEFF,LELD,LXDINT,LYDINT,LZDINT,LEFLD,MADD,
     *             LABFLD,LVECMP,LVEC_2,LVEC_S,LAPROJ,LFPROJ)
C     IF(IPCM.EQ.1 .AND. IEFP.EQ.0) THEN
         LWRK1 = LAST
         LWRK3 = LWRK1 + L3
         LSCR  = LWRK3 + L3
         LQSET = LSCR  + L1
         LQSNT = LQSET + NTS
         LQSEDT= LQSNT + NTS
         LAST  = LQSEDT+ NTS
C
         LZAN  = LAST
         LAST  = LZAN  + NAT
         NEED  = LAST  - LOADFM - 1
C
C     ----- GET MEMORY -----
C
      CALL GETFM(NEED)
C
      CALL DFTB_GET_ZAN(0,NAT,X(LZAN),ZAN)
      CALL DAREAD(IDAF,IODA,X(LS),L2,12,0)
C
C     ----- UPDATE SHIFTES -----
C
      CALL VCLR(X(LQ0),1,NTS)
      CALL PCMFLD(X(LFAO),DUMMY,DMES,X(LSOL),
     *           X(LCHG),X(LWRK1),X(LDISV),
     *           X(LDIS1),X(LDIS2),X(LWRK3),X(LSCR),X(LDMATM),
     *           X(LCQEF),X(LQEFF),X(LELD),X(LQPOT),X(LQFLD),
     *           X(LVPOT),X(LSE),X(LDE),
     *           X(LQ0),X(LQ1),X(LQ2),X(LQ3),
     *           X(LD0),X(LQA),X(LDIMAT),X(LQREP),
     *           X(LTMP0),X(LTMP1),X(LTMP2),X(LIPVT),
     *           X(LPOTTMP),X(LRMUL),
     *           TCH,QET,QETN,QESC,ISD,NFT27,IPCFP,L1,L2,NTSPAR,
     *           X(LAXYZCT),X(LAXYZCT+MXTS),X(LAXYZCT+MXTS*2),
     *           X(LXYZCT2),X(LXYZCT2+MXTS),X(LXYZCT2+MXTS*2),
     *           X(LAXYZCT+MXTS*3),X(LQSNT),X(LQSET),
     *           X(LQSEDT),X(LQFS),
     *           X(LVECMUL),X(LCHG2),X(LQIND+MXTS),X(LPEL),
     *           X(LASCCRD),X(LASCCHG),X(LASCDIP),X(LASCQAD),
     *           X(LXYZRE),X(LXYZRE+MXSP),X(LXYZRE+MXSP*2),
     *           X(LXYZRE+MXSP*3),X(LISPHE),X(LLIST),X(LEPSPCM))
      !! P IS ALREADY DOUBLED IN DENTD SUBOUROUTINE.
      !! HOWEVER, X(LSOL) IS DOUBLED IN PCMFLD (DUE TO SPIN?),
      !! SO X(LSOL) HAS TO BE SCALED BY HALF.
      NSEQ = 0
      ISH = 0
      DO IAT = 1, NAT
        DO K = 1, MAXANG(ISPE(IAT))*MAXANG(ISPE(IAT))
          MU = IND(IAT) + K
          NSEQ = NSEQ + MU
          !! Note that DMES is already multiplied with two
          IF (SRSCC) THEN
            IF (K.EQ.1.OR.K.EQ.2.OR.K.EQ.5) THEN
              ISH = ISH + 1
              SHIFTES(ISH) = SHIFTES(ISH) + X(LSOL+NSEQ-1)*HALF
            END IF
          ELSE
           IF (K.EQ.1) SHIFTES(IAT) = SHIFTES(IAT) + X(LSOL+NSEQ-1)*HALF
          END IF
        END DO
      END DO
C
C     ----- UPDATE SHIFTXY -----
C
      CALL VCLR(X(LSOL),1,L2)
      CALL VCLR(X(LQ0),1,NTS)
      CALL PCMFLD(X(LFAO),DUMMY,DMXY,X(LSOL),
     *           X(LCHG),X(LWRK1),X(LDISV),
     *           X(LDIS1),X(LDIS2),X(LWRK3),X(LSCR),X(LDMATM),
     *           X(LCQEF),X(LQEFF),X(LELD),X(LQPOT),X(LQFLD),
     *           X(LVPOT),X(LSE),X(LDE),
     *           X(LQ0),X(LQ1),X(LQ2),X(LQ3),
     *           X(LD0),X(LQA),X(LDIMAT),X(LQREP),
     *           X(LTMP0),X(LTMP1),X(LTMP2),X(LIPVT),
     *           X(LPOTTMP),X(LRMUL),
     *           TCH,QET,QETN,QESC,ISD,NFT27,IPCFP,L1,L2,NTSPAR,
     *           X(LAXYZCT),X(LAXYZCT+MXTS),X(LAXYZCT+MXTS*2),
     *           X(LXYZCT2),X(LXYZCT2+MXTS),X(LXYZCT2+MXTS*2),
     *           X(LAXYZCT+MXTS*3),X(LQSNT),X(LQSET),
     *           X(LQSEDT),X(LQFS),
     *           X(LVECMUL),X(LCHG2),X(LQIND+MXTS),X(LPEL),
     *           X(LASCCRD),X(LASCCHG),X(LASCDIP),X(LASCQAD),
     *           X(LXYZRE),X(LXYZRE+MXSP),X(LXYZRE+MXSP*2),
     *           X(LXYZRE+MXSP*3),X(LISPHE),X(LLIST),X(LEPSPCM))
      !! (X+Y) IS NOT DOUBLED, AND X(LSOL) IS DOUBLED DUE TO SPIN (?),
      !! SO NO NEED TO SCALE.
      NSEQ = 0
      ISH = 0
      DO IAT = 1, NAT
        DO K = 1, MAXANG(ISPE(IAT))*MAXANG(ISPE(IAT))
          MU = IND(IAT) + K
          NSEQ = NSEQ + MU
          IF (SRSCC) THEN
            IF (K.EQ.1.OR.K.EQ.2.OR.K.EQ.5) THEN
              ISH = ISH + 1
              SHIFTXY(ISH) = SHIFTXY(ISH) + X(LSOL+NSEQ-1)
            END IF
          ELSE
            IF (K.EQ.1) SHIFTXY(IAT) = SHIFTXY(IAT) + X(LSOL+NSEQ-1)
          END IF
        END DO
      END DO
C
      CALL DFTB_GET_ZAN(1,NAT,X(LZAN),ZAN)
C
      CALL RETFM(NEED)
C
      ITDFG = 0
C
      if(l2.lt.0) write(6,*) dm(1)
      RETURN
C
      END SUBROUTINE DFTB_TD_GRAD_PCM
C
C-----------------------------------------------------------------------
C*MODULE DFTBTD    *DECK DFTB_LAMBDA
C>
C>    @brief   TD-DFTB
C>
C>    @details Calculate lambda (degree of charge transfer)
C>             See J. Chem. Phy. 2015, 143, 134120.
C>
C>    @author Yoshio Nishimoto
C>    - Oct, 2019- Subroutine written
C>
C>           --- INPUT ---
C>    @param NAT     Number of atoms
C>    @param L1      Number of basis functions
C>    @param I       Index of molecular orbital
C>    @param J       Index of molecular orbital
C>    @param C       Coordinates
C>    @parma VEC     MO coefficients
C>    @param S       Overlap matrix
C>           --- OUTPUT ---
C>    @param DUMO    lambda?
C>
C
      SUBROUTINE DFTB_LAMBDA(DUMO,L1,NAT,I,J,C,VEC,S)
      USE MX_LIMITS,ONLY:MXATM
C
      IMPLICIT NONE
C     Declarations of arguments
      DOUBLE PRECISION DUMO,C(3,NAT),VEC(L1,*),S(*)
      INTEGER L1,NAT,I,J
C     Declarations of common blocks
      DOUBLE PRECISION ETEMP,DFTBDP,DAMPXHE,HUBDER
      DOUBLE PRECISION ZREF,SPNCST,SPE
      DOUBLE PRECISION PARAMDIR
      INTEGER NSPE
      INTEGER MAXANG,ISPE,IND,IDFTBD
      DOUBLE PRECISION SKDIM,SKSPIN,QREFL
      DOUBLE PRECISION HUBBL,QREF,HUBB,SKCUT2
      INTEGER NEEDSK,LSKHTAB,LSKSTAB
      INTEGER LSKGRID,LSKSELF
C     Other declarations
      INTEGER, PARAMETER :: MXSPE=10
C
      COMMON /DFTBPR/ ETEMP,DFTBDP(MXSPE*14),DAMPXHE,HUBDER(MXSPE),
     *                ZREF(MXATM),SPNCST(6,MXSPE),SPE(MXATM),NSPE,
     *                MAXANG(MXATM),ISPE(MXATM),IND(MXATM+1),IDFTBD,
     *                PARAMDIR
      COMMON /DFTBSK/ SKDIM(MXSPE,MXSPE),SKSPIN(MXSPE),QREFL(3,MXSPE),
     *                HUBBL(3,MXSPE),QREF(MXSPE),HUBB(MXSPE),SKCUT2,
     *                NEEDSK,LSKHTAB(MXSPE,MXSPE),LSKSTAB(MXSPE,MXSPE),
     *                LSKGRID(MXSPE,MXSPE),LSKSELF(MXSPE)
C
      DOUBLE PRECISION QII(NAT),QJJ(NAT)
      DOUBLE PRECISION aa,bb,dist
      INTEGER iat,jat,mu,nseq,nu
      DOUBLE PRECISION oii,oij,ojj,pi,pisq,sigi,sigj,ss,tmpii,tmpjj,val
C
      PI = DACOS(-1.0D+00)
      PISQ = SQRT(PI)
C
      !! eq. (45): 1/2 is omitted
      CALL VCLR(QII,1,NAT)
      CALL VCLR(QJJ,1,NAT)
      DO IAT = 1, NAT
        TMPII = 0.0D+00
        TMPJJ = 0.0D+00
        DO MU = IND(IAT)+1, IND(IAT+1)
          DO NU = 1, L1
            CALL DFTB_CNVSQ(MU,NU,NSEQ)
            SS = S(NSEQ)
            TMPII = TMPII + (VEC(MU,I)*VEC(NU,I)+VEC(NU,I)*VEC(MU,I))*SS
            TMPJJ = TMPJJ + (VEC(MU,J)*VEC(NU,J)+VEC(NU,J)*VEC(MU,J))*SS
          END DO
        END DO
        QII(IAT) = TMPII
        QJJ(IAT) = TMPJJ
      END DO
C
      !! eq 52
      OII = 0.0D+00
      OIJ = 0.0D+00
      OJJ = 0.0D+00
      DO IAT = 1, NAT
        ! the scaling with 0.3125 is incorrect
        SIGI = 1.0D+00/(PISQ*HUBB(ISPE(IAT)))*0.3125D+00
        DO JAT = 1, NAT
          SIGJ = 1.0D+00/(PISQ*HUBB(ISPE(JAT)))*0.3125D+00
          DIST = (C(1,IAT)-C(1,JAT))**2
     *         + (C(2,IAT)-C(2,JAT))**2
     *         + (C(3,IAT)-C(3,JAT))**2
          AA = (2.0D+00*PI*(SIGI*SIGI+SIGJ*SIGJ))**1.5D+00
          BB = -0.5D+00*DIST/(SIGI*SIGI+SIGJ*SIGJ)
          VAL = 1.0D+00/AA * EXP(BB)
          OII = OII + QII(IAT)*QII(JAT)*VAL
          OIJ = OIJ + QII(IAT)*QJJ(JAT)*VAL
          OJJ = OJJ + QJJ(IAT)*QJJ(JAT)*VAL
        END DO
      END DO
C
      !! the last coefficient of eq (48)
      DUMO = OIJ/SQRT(OII*OJJ)
C
      RETURN
C
      END SUBROUTINE DFTB_LAMBDA
C
C-----------------------------------------------------------------------
C*MODULE DFTBTD    *DECK DFTB_LCSHIFT3
C>
C>    @brief   TD-LC-DFTB(/PCM) gradient
C>
C>    @details Calculate the coupling matrix of long-range correction
C>             contributions (exchange-type term)
C>
C>    @author Yoshio Nishimoto
C>    - Oct, 2019- Subroutine written
C>
C>           --- INPUT ---
C>    @param NAT     Number of atoms
C>    @param L1      Number of basis functions
C>    @param ISPE    Index of species
C>    @parma MAXANG  Maximum angular momentum of each species
C>    @param IND     AO index of each atom
C>    @param DD      pseudo-density matrix
C>    @param S       Overlap matrix
C>    @param GAMLCS  LC-gamma in squared matrix
C>           --- OUTPUT ---
C>    @param FAO     coupling matrix of the exchange-type term
C>
C
      SUBROUTINE DFTB_LCSHIFT3(NAT,L1,ISPE,MAXANG,IND,DD,S,GAMLCS,FAO)
      USE constants,ONLY:ZERO,ONE
      IMPLICIT NONE
C     Declarations of arguments
      DOUBLE PRECISION DD(*),S(*),GAMLCS(NAT,NAT),FAO(*)
      INTEGER NAT,L1,ISPE(*),MAXANG(*),IND(*)
C     Declarations of common blocks
      DOUBLE PRECISION X
C     Other declarations
      DOUBLE PRECISION gamtmp
      INTEGER i,iat,jat,l2,l3,last,ldsq,loadfm,lss,lwrk,lwrk1,lwrk2,mu
      INTEGER need,ngotmx,nu
C
      COMMON /FMCOM / X(1)
C
      L2=L1*(L1+1)/2
      L3=L1*L1
      CALL VALFM(LOADFM)
      CALL GOTFM(NGOTMX)
      LSS     = LOADFM  + 1
      LDSQ    = LSS     + L3
      LWRK    = LDSQ    + L3
      LWRK1   = LWRK    + L3
      LWRK2   = LWRK1   + L3
      LAST    = LWRK2   + L3
      NEED    = LAST - LOADFM - 1
      CALL GETFM(NEED)
C
      DO IAT = 1, NAT
        DO JAT = 1, NAT
          GAMTMP = GAMLCS(IAT,JAT)
          DO MU = IND(IAT)+1, IND(IAT)+MAXANG(ISPE(IAT))**2
            DO NU = IND(JAT)+1, IND(JAT)+MAXANG(ISPE(JAT))**2
              X(LWRK2+MU-1+L1*(NU-1)) = GAMTMP
            END DO
          END DO
        END DO
      END DO
C
      CALL CPYTSQ(S,X(LSS),L1,1)
      CALL DCOPY(L3,DD,1,X(LDSQ),1)
C
      !! gamma_{AB}
      CALL DGEMM('N','N',L1,L1,L1,ONE,X(LSS),L1,X(LDSQ),L1,
     *           ZERO,X(LWRK1),L1)
      CALL DGEMM('N','N',L1,L1,L1,ONE,X(LWRK1),L1,X(LSS),L1,
     *           ZERO,X(LWRK),L1)
      DO I = 1, L3
        X(LWRK+I-1) = X(LWRK+I-1)*X(LWRK2+I-1)
      END DO
      !! gamma_{AD}
      DO I = 1, L3
        X(LWRK1+I-1) = X(LWRK1+I-1)*X(LWRK2+I-1)
      END DO
       
      !! Prepare for gamma_{CD}
      DO I = 1, L3
        X(LDSQ+I-1) = X(LDSQ+I-1)*X(LWRK2+I-1)
      END DO
      CALL DGEMM('N','N',L1,L1,L1,ONE,X(LSS),L1,X(LDSQ),L1,
     *           ONE,X(LWRK1),L1)
C
      !! gamma_{AD} + gamma_{CD}
      CALL DGEMM('N','N',L1,L1,L1,ONE,X(LWRK1),L1,X(LSS),L1,
     *           ONE,X(LWRK),L1)
C
      !! gamma_{BC}
      CALL DCOPY(L3,DD,1,X(LDSQ),1)
      CALL DGEMM('N','N',L1,L1,L1,ONE,X(LDSQ),L1,X(LSS),L1,
     *           ZERO,X(LWRK1),L1)
      DO I = 1, L3
        X(LWRK1+I-1) = X(LWRK1+I-1)*X(LWRK2+I-1)
      END DO
      CALL DGEMM('N','N',L1,L1,L1,ONE,X(LSS),L1,X(LWRK1),L1,
     *           ONE,X(LWRK),L1)
C
      CALL DAXPY(L3,-0.125D+00,X(LWRK),1,FAO,1)
C
      CALL RETFM(NEED)
C
      RETURN
C
      END SUBROUTINE DFTB_LCSHIFT3
