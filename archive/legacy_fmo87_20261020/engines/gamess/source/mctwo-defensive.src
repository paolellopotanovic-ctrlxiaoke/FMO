C*MODULE MCTWO-DEFENSIVE   *DECK AHPDDI
      SUBROUTINE AHPDDI(NMOS,NCOR,NOCC,NROT
     *,                 OPDM,TPDM,FCOR,FVAL
     *,                 LAGN,IROT,DIAH,BUFF,TVEC,PROD)
      use mx_limits, only: mxao
C
C -----------------------------------------------------------------
C
C  NEWTON-RAPHSON MCSCF. CALLED FROM SUBROUTINE NTNDVD.
C  COMPUTE THE PRODUCT OF THE AUGMENTED HESSIAN WITH THE TRIAL VECTOR.
C  SEE MOTECC-90 (ESCOM), PAGE 293, EQNS (B.71-75).
C  SEE YARKONY, CHEM. PHYS. LETT, VOLUME 77, PAGE 634.
C  TO MAKE USE OF PARALLEL TRANSFORMATION (TRANDDI), CORE AND
C  ACTIVE MO INDICES ARE SUBSETS OF THE OCCUPIED MO LIST.
C  TRANDDI IS CALLED FROM SUBROUTINE TRFMCX (IN MODULE TRANS.SRC).
C
C  SYMBOLS:
C     NMOS = TOTAL NUMBER OF MOS
C     NCOR = NUMBER OF CORE MOS
C     NOCC = NUMBER OF OCCUPIED MOS
C     NROT = NUMBER OF ROTATIONS
C     OPDM = 1-EL DENSITY SPANNING ACTIVE INDICES
C     TPDM = 2-EL DENSITY SPANNING ACTIVE INDICES
C     FCOR = CORE FOCK OPERATOR
C     FVAL = VALENCE FOCK OPERATOR
C     LAGN = LAGRANGIAN MATRIX
C     IROT = INDEX OF ROTATIONS
C     DIAH = DIAGONAL ELEMENTS OF AUGMENTED HESSIAN
C     TVEC = TRIAL VECTOR OF DAVIDSON SOLVER
C     BUFF = MESSAGE BUFFER
C     PROD = PRODUCT VECTOR OF DAVIDSON SOLVER [OUTPUT]
C
C -----------------------------------------------------------------
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (HALF=0.5D+00)
      PARAMETER (TWO=2.0D+00)
      PARAMETER (FOUR=4.0D+00)
      PARAMETER (EIGHT=8.0D+00)
      LOGICAL INJ,INK,JNK,KEL,LEJ,ANB, BENCH
      INTEGER IROT(NOCC,*), DDI_NP,DDI_ME, A,AR,B,BR,AB,BA
      DOUBLE PRECISION  OPDM(*),TPDM(*),FCOR(*),FVAL(*)
      DOUBLE PRECISION  LAGN(NOCC,*),DIAH(*),BUFF(*),TVEC(*),PROD(*)
      COMMON /IJPAIR/ IA(MXAO)
C
      INTEGER         D_OOOO,D_VOOO,D_VVOO,D_VOVO,D_VVVO,D_VVVV,
     *                D_OOOOAB,D_OOOOBB,D_VOOOAB,D_VOOOBA,D_VOOOBB,
     *                D_VVOOAB,D_VVOOBA,D_VVOOBB,D_VOVOAB,D_VOVOBB,
     *                D_U,D_UB,D_E,D_EB
      LOGICAL         NDOOOO,NDVOOO,NDVVOO,NDVOVO,NDVVVO,NDVVVV,NDCORE,
     *                NDVVOOBA,NDVVOOAB,NDVVOOBB,NDVOVOAB,NDVOVOBB,
     *                NDVOOOBA,NDVOOOAB,NDVOOOBB,NDOOOOAB,NDOOOOBB
      COMMON /TRFDMS/ D_OOOO,D_VOOO,D_VVOO,D_VOVO,D_VVVO,D_VVVV,
     *                D_OOOOAB,D_OOOOBB,D_VOOOAB,D_VOOOBA,D_VOOOBB,
     *                D_VVOOAB,D_VVOOBA,D_VVOOBB,D_VOVOAB,D_VOVOBB,
     *                D_U,D_UB,D_E,D_EB,
     *                NDOOOO,NDVOOO,NDVVOO,NDVOVO,NDVVVO,NDVVVV,NDCORE,
     *                NDVVOOBA,NDVVOOAB,NDVVOOBB,NDVOVOAB,NDVOVOBB,
     *                NDVOOOBA,NDVOOOAB,NDVOOOBB,NDOOOOAB,NDOOOOBB
C
      NCP1 = NCOR + 1
      NVIR = NMOS - NOCC
      NOTR = (NOCC*NOCC+NOCC)/2
      NVTR = (NVIR*NVIR+NVIR)/2
      NVSQ = NVIR*NVIR
      CALL DDI_NPROC(DDI_NP,DDI_ME)
C
      BENCH = DDI_ME.EQ.0   !  SWITCH FOR BENCHMARK TIMING
      BENCH = .FALSE.
      ICHANL = 6            !  CAN BE UNIQUE TO A PROCESS
C
C  1-EL CONTRIBUTIONS TO PRODUCT
C
      CALL DCOPY(NROT+1,0.0D+00,0,PROD,1)
      DO M = 1, NMOS
        IF (MOD(M,DDI_NP).EQ.DDI_ME) THEN
C
C  1. HESS(IM|IN) = 2*( FCOR(MN) + FVAL(MN) ), I CORE, M,N NOT CORE
C
          IF (M.GT.NCOR) THEN
            DO N = NCP1, M
              MN = IA(M) + N
              FAC = TWO*( FCOR(MN) + FVAL(MN) )
              DO I = 1, NCOR
                IMR = IROT(I,M)
                INR = IROT(I,N)
                IF (IMR.NE.0.AND.INR.NE.0) THEN
                  IF (IMR.NE.INR) THEN
                    IX = IMR + 1
                    JX = INR + 1
                    PROD(IX) = PROD(IX) + FAC*TVEC(JX)
                    PROD(JX) = PROD(JX) + FAC*TVEC(IX)
                  END IF
                END IF
              END DO
            END DO
          END IF
C
C  2. H(IM|JN) = OPDM(IJ)*FCOR(MN), I,J VALENCE, M,N GENERAL
C
          DO N = 1, M
            MN = IA(M) + N
            DO I = NCP1, NOCC
              IMR = IROT(I,M)
              IF (IMR.NE.0) THEN
                DO J = NCP1, I
                  JNR = IROT(J,N)
                  IF (JNR.NE.0) THEN
                    IF (IMR.NE.JNR) THEN
                      ID = I - NCOR
                      JD = J - NCOR
                      IJD = IA(ID) + JD
                      DIJ = OPDM(IJD)
                      FAC = DIJ*FCOR(MN)
                      IF (M.EQ.N) FAC = FAC*HALF
                      IF (I.EQ.J) FAC = FAC*HALF
                      IF (M.LT.I) FAC = -FAC
                      IF (N.LT.J) FAC = -FAC
                      IX = IMR + 1
                      JX = JNR + 1
                      PROD(IX) = PROD(IX) + FAC*TVEC(JX)
                      PROD(JX) = PROD(JX) + FAC*TVEC(IX)
                    END IF
                  END IF
                END DO
              END IF
            END DO
            DO I = NCP1, NOCC
              INR = IROT(I,N)
              IF (INR.NE.0) THEN
                DO J = NCP1, I
                  JMR = IROT(J,M)
                  IF (JMR.NE.0) THEN
                    IF (INR.NE.JMR) THEN
                      ID = I - NCOR
                      JD = J - NCOR
                      IJD = IA(ID) + JD
                      DIJ = OPDM(IJD)
                      FAC = DIJ*FCOR(MN)
                      IF (M.EQ.N) FAC = FAC*HALF
                      IF (I.EQ.J) FAC = FAC*HALF
                      IF (N.LT.I) FAC = -FAC
                      IF (M.LT.J) FAC = -FAC
                      IX = INR + 1
                      JX = JMR + 1
                      PROD(IX) = PROD(IX) + FAC*TVEC(JX)
                      PROD(JX) = PROD(JX) + FAC*TVEC(IX)
                    END IF
                  END IF
                END DO
              END IF
            END DO
          END DO
C
C  LAGRANGIAN CONTRIBUTIONS
C
          DO J = 1, NOCC
C
C  3. ROW-COLUMN OF AUGMENTED HESSIAN
C
            JMR = IROT(J,M)
            IF (JMR.NE.0) THEN
              FAC = LAGN(J,M)
              IF (J.GT.M) FAC = -FAC
              IX = 1
              JX = JMR + 1
              PROD(IX) = PROD(IX) + FAC*TVEC(JX)
              PROD(JX) = PROD(JX) + FAC*TVEC(IX)
            END IF
C
C  4. (UNDOCUMENTED)
C
            DO K = 1, NMOS
              IF (M.LE.NOCC.OR.K.LE.NOCC) THEN
                IF (M.LE.NOCC) MKR = IROT(M,K)
                IF (M.GT.NOCC) MKR = IROT(K,M)
                IF (MKR.NE.0) THEN
                  JKR = IROT(J,K)
                  IF (JKR.NE.0) THEN
                    IF (MKR.NE.JKR) THEN
                      FAC = LAGN(J,M)*0.5D+00
                      IF (K.GT.M) FAC = -FAC
                      IF (K.LT.J) FAC = -FAC
                      IX = MKR + 1
                      JX = JKR + 1
                      PROD(IX) = PROD(IX) + FAC*TVEC(JX)
                      PROD(JX) = PROD(JX) + FAC*TVEC(IX)
                    END IF
                  END IF
                END IF
              END IF
            END DO
          END DO
        END IF   !  PARALLEL
      END DO   ! M
C
      IF (BENCH) THEN
        WRITE(ICHANL,9000) 'ONE-ELEC'
        CALL TIMIT(1)
      END IF
C
C  2-EL CONTRIBUTIONS TO PRODUCT
C
      CALL DDI_DLBRESET()
      CALL DDI_DLBNEXT(MYTASK)
      LOCTSK = -1
C
C ----------- CONTRIBUTIONS FROM (OO|OO) TYPE INTEGRALS -----------
C
      CALL DDI_DISTRIB(D_OOOO,DDI_ME,ILO,IHI,JLO,JHI)
C
C  (ACTIVE-CORE|**) TYPES
C
      DO I = NCP1, NOCC
        DO J = 1, NCOR
          IJ = IA(MAX(I,J)) + MIN(I,J)
          LOCTSK = LOCTSK + 1
          IF (LOCTSK.EQ.MYTASK) THEN
            CALL DDI_GET(D_OOOO,1,NOTR,IJ,IJ,BUFF)
C
C  1. H(IJ|KL) <- 8*(IJ|KL), J,L CORE
C
            IJR = IROT(I,J)
            IF (IJR.NE.0) THEN
              DO K = NCP1, I
                MX = J
                IF (I.NE.K) MX = NCOR
                DO L = 1, MX
                  KLR = IROT(K,L)
                  IF (KLR.NE.0) THEN
                    IF (IJR.NE.KLR) THEN
                      KL = IA(MAX(K,L)) + MIN(K,L)
                      ERI = BUFF(KL)
                      FAC = EIGHT*ERI
                      IX = IJR + 1
                      JX = KLR + 1
                      PROD(IX) = PROD(IX) + FAC*TVEC(JX)
                      PROD(JX) = PROD(JX) + FAC*TVEC(IX)
                    END IF
                  END IF
                END DO
              END DO
            END IF
C
C  2. H(IL|JK) <- -2*(IJ|KL), J,L CORE
C
            DO K = I, NOCC
              MX = J
              IF (I.NE.K) MX = NCOR
              JKR = IROT(J,K)
              IF (JKR.NE.0) THEN
                DO L = 1, MX
                  ILR = IROT(I,L)
                  IF (ILR.NE.0) THEN
                    IF (JKR.NE.ILR) THEN
                      KL = IA(MAX(K,L)) + MIN(K,L)
                      ERI = BUFF(KL)
                      FAC = -TWO*ERI
                      IX = JKR + 1
                      JX = ILR + 1
                      PROD(IX) = PROD(IX) + FAC*TVEC(JX)
                      PROD(JX) = PROD(JX) + FAC*TVEC(IX)
                    END IF
                  END IF
                END DO
              END IF
            END DO
C
C  3. H(JK|LM) <- -4*D(IK)*(IJ|LM), J,M CORE
C
            ID = I - NCOR
            DO K = NCP1, NOCC
              KD = K - NCOR
              IKD = IA(MAX(ID,KD)) + MIN(ID,KD)
              JKR = IROT(K,J)
              IF (JKR.NE.0) THEN
                DIK = FOUR*OPDM(IKD)
                DO L = NCP1, K
                  MX = J
                  IF (K.NE.L) MX = NCOR
                  DO M = 1, MX
                    LMR = IROT(L,M)
                    IF (LMR.NE.0) THEN
                      IF (JKR.NE.LMR) THEN
                        LM  = IA(MAX(L,M)) + MIN(L,M)
                        ERI = BUFF(LM)
                        FAC = -DIK*ERI
                        IX = JKR + 1
                        JX = LMR + 1
                        PROD(IX) = PROD(IX) + FAC*TVEC(JX)
                        PROD(JX) = PROD(JX) + FAC*TVEC(IX)
                      END IF
                    END IF
                  END DO
                END DO
              END IF
            END DO
C
C  4. H(IL|JK) <- D(KM)*(IJ|LM), J,L CORE
C
            DO K = NCP1, NOCC
              IF (I.LE.K) THEN
                KD = K - NCOR
                JKR = IROT(K,J)
                IF (JKR.NE.0) THEN
                  MX = J
                  IF (I.NE.K) MX = NCOR
                  DO L = 1, MX
                    ILR = IROT(I,L)
                    IF (ILR.NE.0) THEN
                      IF (JKR.NE.ILR) THEN
                        DO M = NCP1, NOCC
                          LM = IA(MAX(L,M)) + MIN(L,M)
                          ERI = BUFF(LM)
                          MD = M - NCOR
                          KMD = IA(MAX(KD,MD)) + MIN(KD,MD)
                          DKM = OPDM(KMD)
                          FAC = DKM*ERI
                          IX = JKR + 1
                          JX = ILR + 1
                          PROD(IX) = PROD(IX) + FAC*TVEC(JX)
                          PROD(JX) = PROD(JX) + FAC*TVEC(IX)
                        END DO
                      END IF
                    END IF
                  END DO
                END IF
              END IF
            END DO
C
C  5. H(IJ|KL) <- -4*D(KM)*(IJ|LM), J,L CORE
C
            IJR = IROT(I,J)
            IF (IJR.NE.0) THEN
              DO K = NCP1, I
                KD = K - NCOR
                MX = J
                IF (I.NE.K) MX = NCOR
                DO M = NCP1, NOCC
                  MD = M - NCOR
                  KMD = IA(MAX(KD,MD)) + MIN(KD,MD)
                  DKM = OPDM(KMD)*FOUR
                  DO L = 1, MX
                    KLR = IROT(K,L)
                    IF (KLR.NE.0) THEN
                      IF (IJR.NE.KLR) THEN
                        LM  = IA(MAX(L,M)) + MIN(L,M)
                        ERI = BUFF(LM)
                        FAC = -DKM*ERI
                        IX = IJR + 1
                        JX = KLR + 1
                        PROD(IX) = PROD(IX) + FAC*TVEC(JX)
                        PROD(JX) = PROD(JX) + FAC*TVEC(IX)
                      END IF
                    END IF
                  END DO
                END DO
              END DO
            END IF
C
C  6. H(JK|LM) <- D(IL)*(IJ|KM), J,M CORE
C
            ID = I - NCOR
            DO K = NCP1, NOCC
              JKR = IROT(K,J)
              IF (JKR.NE.0) THEN
                DO L = NCP1, K
                  LD = L - NCOR
                  ILD = IA(MAX(ID,LD)) + MIN(ID,LD)
                  DIL = OPDM(ILD)
                  MX = J
                  IF (K.NE.L) MX = NCOR
                  DO M = 1, MX
                    LMR = IROT(L,M)
                    IF (LMR.NE.0) THEN
                      IF (JKR.NE.LMR) THEN
                        KM = IA(MAX(K,M)) + MIN(K,M)
                        ERI = BUFF(KM)
                        FAC = DIL*ERI
                        IX = JKR + 1
                        JX = LMR + 1
                        PROD(IX) = PROD(IX) + FAC*TVEC(JX)
                        PROD(JX) = PROD(JX) + FAC*TVEC(IX)
                      END IF
                    END IF
                  END DO
                END DO
              END IF
            END DO
C
C  7. H(JK|LN) <- 2*D(IK|LM)*(IJ|MN), J,N CORE
C
            ID = I - NCOR
            DO K = NCP1, NOCC
              KD = K - NCOR
              IKD = IA(MAX(ID,KD)) + MIN(ID,KD)
              JKR = IROT(K,J)
              IF (JKR.NE.0) THEN
                DO L = NCP1, K
                  LD = L - NCOR
                  MX = J
                  IF (K.NE.L) MX = NCOR
                  DO M = NCP1, NOCC
                    MD = M - NCOR
                    LMD = IA(MAX(LD,MD)) + MIN(LD,MD)
                    IKLMD = IA(MAX(IKD,LMD)) + MIN(IKD,LMD)
                    DIKLM = TPDM(IKLMD)*TWO
                    DO N = 1, MX
                      LNR = IROT(L,N)
                      IF (LNR.NE.0) THEN
                        IF (JKR.NE.LNR) THEN
                          MN = IA(MAX(M,N)) + MIN(M,N)
                          ERI = BUFF(MN)
                          FAC = DIKLM*ERI
                          IX = JKR + 1
                          JX = LNR + 1
                          PROD(IX) = PROD(IX) + FAC*TVEC(JX)
                          PROD(JX) = PROD(JX) + FAC*TVEC(IX)
                        END IF
                      END IF
                    END DO
                  END DO
                END DO
              END IF
            END DO
C
C  8. H(IJ|KL) <- 4*D(LM)*(IJ|KM), J CORE
C
            ID = I - NCOR
            IJR = IROT(I,J)
            IF (IJR.NE.0) THEN
              DO K = NCP1, NOCC
                DO L = NCP1, K-1
                  KLR = IROT(K,L)
                  IF (KLR.NE.0) THEN
                    IF (IJR.NE.KLR) THEN
                      LD = L - NCOR
                      DO M = NCP1, NOCC
                        KM = IA(MAX(K,M)) + MIN(K,M)
                        ERI = BUFF(KM)
                        MD = M - NCOR
                        LMD = IA(MAX(LD,MD)) + MIN(LD,MD)
                        DLM = OPDM(LMD)*FOUR
                        FAC = DLM*ERI
                        IX = IJR + 1
                        JX = KLR + 1
                        PROD(IX) = PROD(IX) + FAC*TVEC(JX)
                        PROD(JX) = PROD(JX) + FAC*TVEC(IX)
                      END DO
                    END IF
                  END IF
                END DO
              END DO
            END IF
C
C  9. H(JM|KL) <- -D(IL)*(IJ|KM), J CORE
C
            ID = I - NCOR
            DO K = NCP1, NOCC
              DO L = NCP1, K-1
                KLR = IROT(K,L)
                IF (KLR.NE.0) THEN
                  LD = L - NCOR
                  ILD = IA(MAX(ID,LD)) + MIN(ID,LD)
                  DIL = OPDM(ILD)
                  DO M = NCP1, NOCC
                    JMR = IROT(J,M)
                    IF (JMR.NE.0) THEN
                      IF (KLR.NE.JMR) THEN
                        KM = IA(MAX(K,M)) + MIN(K,M)
                        ERI = BUFF(KM)
                        FAC = -DIL*ERI
                        IX = KLR + 1
                        JX = JMR + 1
                        PROD(IX) = PROD(IX) + FAC*TVEC(JX)
                        PROD(JX) = PROD(JX) + FAC*TVEC(IX)
                      END IF
                    END IF
                  END DO
                END IF
              END DO
            END DO
C
C  10. H(IL|JK) <- -D(LM)*(IJ|KM), J CORE
C
            DO K = NCP1, NOCC
              JKR = IROT(K,J)
              IF (JKR.NE.0) THEN
                DO L = NCP1, I-1
                  LD = L - NCOR
                  ILR = IROT(I,L)
                  IF (ILR.NE.0) THEN
                    IF (JKR.NE.ILR) THEN
                      DO M = NCP1, NOCC
                        KM  = IA(MAX(K,M)) + MIN(K,M)
                        ERI = BUFF(KM)
                        MD = M - NCOR
                        LMD = IA(MAX(LD,MD)) + MIN(LD,MD)
                        DLM = OPDM(LMD)
                        FAC = -DLM*ERI
                        IX = JKR + 1
                        JX = ILR + 1
                        PROD(IX) = PROD(IX) + FAC*TVEC(JX)
                        PROD(JX) = PROD(JX) + FAC*TVEC(IX)
                      END DO
                    END IF
                  END IF
                END DO
              END IF
            END DO
C
C  11. H(KL|IJ) <- -4*D(KM)*(IJ|LM),  J CORE
C
            ID = I - NCOR
            IJR = IROT(I,J)
            IF (IJR.NE.0) THEN
              DO K = NCP1, NOCC
                KD = K - NCOR
                DO L = NCP1, K-1
                  KLR = IROT(K,L)
                  IF (KLR.NE.0) THEN
                    IF (IJR.NE.KLR) THEN
                      DO M = NCP1, NOCC
                        LM = IA(MAX(L,M)) + MIN(L,M)
                        ERI = BUFF(LM)
                        MD = M - NCOR
                        KMD = IA(MAX(KD,MD)) + MIN(KD,MD)
                        DKM = OPDM(KMD)*FOUR
                        FAC = -DKM*ERI
                        IX = IJR + 1
                        JX = KLR + 1
                        PROD(IX) = PROD(IX) + FAC*TVEC(JX)
                        PROD(JX) = PROD(JX) + FAC*TVEC(IX)
                      END DO
                    END IF
                  END IF
                END DO
              END DO
            END IF
C
C  12. H(JM|KL) <- D(IK)*(IJ|LM), J CORE
C
            ID = I - NCOR
            DO K = NCP1, NOCC
              KD = K - NCOR
              IKD = IA(MAX(ID,KD)) + MIN(ID,KD)
              DIK = OPDM(IKD)
              DO L = NCP1, K-1
                KLR = IROT(K,L)
                IF (KLR.NE.0) THEN
                  DO M = NCP1, NOCC
                    JMR = IROT(M,J)
                    IF (JMR.NE.0) THEN
                      IF (KLR.NE.JMR) THEN
                        LM = IA(MAX(L,M)) + MIN(L,M)
                        ERI = BUFF(LM)
                        FAC = DIK*ERI
                        IX = KLR + 1
                        JX = JMR + 1
                        PROD(IX) = PROD(IX) + FAC*TVEC(JX)
                        PROD(JX) = PROD(JX) + FAC*TVEC(IX)
                      END IF
                    END IF
                  END DO
                END IF
              END DO
            END DO
C
C  13. H(IM|JN) <- -D(KL|MN)*(IJ|KL), J CORE
C
            DO K = NCP1, NOCC
              KD = K - NCOR
              DO L = NCP1, NOCC
                LD = L - NCOR
                KLD = IA(MAX(KD,LD)) + MIN(KD,LD)
                KL = IA(MAX(K,L)) + MIN(K,L)
                ERI = BUFF(KL)
                DO M = NCP1, I-1
                  MD = M - NCOR
                  IMR = IROT(I,M)
                  IF (IMR.NE.0) THEN
                    DO N = NCP1, NOCC
                      JNR = IROT(N,J)
                      IF (JNR.NE.0) THEN
                        IF (IMR.NE.JNR) THEN
                          ND = N - NCOR
                          MND = IA(MAX(MD,ND)) + MIN(MD,ND)
                          KLMND = IA(MAX(KLD,MND)) + MIN(KLD,MND)
                          DKLMN = TPDM(KLMND)
                          FAC = -DKLMN*ERI
                          IX = IMR + 1
                          JX = JNR + 1
                          PROD(IX) = PROD(IX) + FAC*TVEC(JX)
                          PROD(JX) = PROD(JX) + FAC*TVEC(IX)
                        END IF
                      END IF
                    END DO
                  END IF
                END DO
              END DO
            END DO
C
C  14. H(KL|JN) <- -2*D(IN|LM)*(IJ|KM), J CORE
C
            ID = I - NCOR
            DO K = NCP1, NOCC
              DO L = NCP1, K-1
                LD = L - NCOR
                KLR = IROT(K,L)
                IF (KLR.NE.0) THEN
                  DO M = NCP1, NOCC
                    MD = M - NCOR
                    LMD = IA(MAX(LD,MD)) + MIN(LD,MD)
                    KM = IA(MAX(K,M)) + MIN(K,M)
                    ERI = BUFF(KM)*TWO
                    DO N = NCP1, NOCC
                      JNR = IROT(N,J)
                      IF (JNR.NE.0) THEN
                        IF (KLR.NE.JNR) THEN
                          ND = N - NCOR
                          IND = IA(MAX(ID,ND)) + MIN(ID,ND)
                          INLMD = IA(MAX(IND,LMD)) + MIN(IND,LMD)
                          DINLM = TPDM(INLMD)
                          FAC = -DINLM*ERI
                          IX = KLR + 1
                          JX = JNR + 1
                          PROD(IX) = PROD(IX) + FAC*TVEC(JX)
                          PROD(JX) = PROD(JX) + FAC*TVEC(IX)
                        END IF
                      END IF
                    END DO
                  END DO
                END IF
              END DO
            END DO
C
C  15. H(KL|JN) <- 2*D(KM|IN)*(IJ|LM), J CORE
C
            ID = I - NCOR
            DO K = NCP1, NOCC
              KD = K - NCOR
              DO L = NCP1, K-1
                KLR = IROT(K,L)
                IF (KLR.NE.0) THEN
                  DO M = NCP1, NOCC
                    MD = M - NCOR
                    KMD = IA(MAX(KD,MD)) + MIN(KD,MD)
                    LM = IA(MAX(L,M)) + MIN(L,M)
                    ERI = BUFF(LM)*TWO
                    DO N = NCP1, NOCC
                      JNR = IROT(N,J)
                      IF (JNR.NE.0) THEN
                        IF (KLR.NE.JNR) THEN
                          ND = N - NCOR
                          IND = IA(MAX(ID,ND)) + MIN(ID,ND)
                          KMIND = IA(MAX(KMD,IND)) + MIN(KMD,IND)
                          DKMIN = TPDM(KMIND)
                          FAC = DKMIN*ERI
                          IX = KLR + 1
                          JX = JNR + 1
                          PROD(IX) = PROD(IX) + FAC*TVEC(JX)
                          PROD(JX) = PROD(JX) + FAC*TVEC(IX)
                        END IF
                      END IF
                    END DO
                  END DO
                END IF
              END DO
            END DO
C
C  END DISTRIBUTION LOOPS
C
            CALL DDI_DLBNEXT(MYTASK)
          END IF  ! DLB
        END DO  ! J
      END DO  ! I
C
C  (ACTIVE-ACTIVE|**) TYPES
C
      DO I = NCP1, NOCC
        DO J = NCP1, I
          IJ = IA(MAX(I,J)) + MIN(I,J)
          LOCTSK = LOCTSK + 1
          IF (LOCTSK.EQ.MYTASK) THEN
            CALL DDI_GET(D_OOOO,1,NOTR,IJ,IJ,BUFF)
C
C  16. H(IK|JL) <- -2*(IJ|KL), K,L CORE
C
            INJ = I.NE.J
            DO K = 1, NCOR
              IKR = IROT(I,K)
              IF (IKR.NE.0) THEN
                MX = K
                IF (INJ) MX = NCOR
                DO L = 1, MX
                  JLR = IROT(J,L)
                  IF (JLR.NE.0) THEN
                    IF (IKR.NE.JLR) THEN
                      KL = IA(MAX(K,L)) + MIN(K,L)
                      ERI = BUFF(KL)
                      FAC = -TWO*ERI
                      IX = IKR + 1
                      JX = JLR + 1
                      PROD(IX) = PROD(IX) + FAC*TVEC(JX)
                      PROD(JX) = PROD(JX) + FAC*TVEC(IX)
                    END IF
                  END IF
                END DO
              END IF
            END DO
C
C  17. H(KL|JM) <- D(IK)*(IJ|LM), L,M CORE
C
            ID = I - NCOR
            JD = J - NCOR
            INJ = I.NE.J
            DO K = NCP1, NOCC
              KD = K - NCOR
              INK = I.NE.K
              JNK = J.NE.K
              IF (J.LE.K) THEN
                IKD = IA(MAX(ID,KD)) + MIN(ID,KD)
                DIK = OPDM(IKD)
                DO L = 1, NCOR
                  KLR = IROT(K,L)
                  IF (KLR.NE.0) THEN
                    MX = L
                    IF (JNK) MX = NCOR
                    DO M = 1, MX
                      JMR = IROT(J,M)
                      IF (JMR.NE.0) THEN
                        IF (KLR.NE.JMR) THEN
                          LM = IA(MAX(L,M)) + MIN(L,M)
                          ERI = BUFF(LM)
                          FAC = DIK*ERI
                          IX = KLR + 1
                          JX = JMR + 1
                          PROD(IX) = PROD(IX) + FAC*TVEC(JX)
                          PROD(JX) = PROD(JX) + FAC*TVEC(IX)
                        END IF
                      END IF
                    END DO
                  END IF
                END DO
              END IF
              IF (INJ.AND.(I.LE.K)) THEN
                JKD = IA(MAX(JD,KD)) + MIN(JD,KD)
                DJK = OPDM(JKD)
                DO L = 1, NCOR
                  KLR = IROT(K,L)
                  IF (KLR.NE.0) THEN
                    MX = L
                    IF (INK) MX = NCOR
                    DO M = 1, MX
                      IMR = IROT(I,M)
                      IF (IMR.NE.0) THEN
                        IF (KLR.NE.IMR) THEN
                          LM = IA(MAX(L,M)) + MIN(L,M)
                          ERI = BUFF(LM)
                          FAC = DJK*ERI
                          IX = KLR + 1
                          JX = IMR + 1
                          PROD(IX) = PROD(IX) + FAC*TVEC(JX)
                          PROD(JX) = PROD(JX) + FAC*TVEC(IX)
                        END IF
                      END IF
                    END DO
                  END IF
                END DO
              END IF
            END DO   !  K
C
C  19. H(IK|LM) <- D(JL)*(IJ|KM), K,M CORE
C
            ID = I - NCOR
            JD = J - NCOR
            INJ = I.NE.J
            DO K = 1, NCOR
              IKR = IROT(I,K)
              IF (IKR.NE.0) THEN
                DO L = NCP1, I
                  LD = L - NCOR
                  JLD = IA(MAX(JD,LD)) + MIN(JD,LD)
                  DJL = OPDM(JLD)
                  MX = K
                  IF (I.NE.L) MX = NCOR
                  DO M = 1, MX
                    LMR = IROT(L,M)
                    IF (LMR.NE.0) THEN
                      IF (IKR.NE.LMR) THEN
                        KM = IA(MAX(K,M)) + MIN(K,M)
                        ERI = BUFF(KM)
                        FAC = DJL*ERI
                        IX = IKR + 1
                        JX = LMR + 1
                        PROD(IX) = PROD(IX) + FAC*TVEC(JX)
                        PROD(JX) = PROD(JX) + FAC*TVEC(IX)
                      END IF
                    END IF
                  END DO
                END DO
              END IF
              IF (INJ) THEN
                JKR = IROT(J,K)
                IF (JKR.NE.0) THEN
                  DO L = NCP1, J
                    LD = L - NCOR
                    ILD = IA(MAX(ID,LD)) + MIN(ID,LD)
                    DIL = OPDM(ILD)
                    MX = K
                    IF (J.NE.L) MX = NCOR
                    DO M = 1, MX
                      LMR = IROT(L,M)
                      IF (LMR.NE.0) THEN
                        IF (JKR.NE.LMR) THEN
                          KM = IA(MAX(K,M)) + MIN(K,M)
                          ERI = BUFF(KM)
                          FAC = DIL*ERI
                          IX = JKR + 1
                          JX = LMR + 1
                          PROD(IX) = PROD(IX) + FAC*TVEC(JX)
                          PROD(JX) = PROD(JX) + FAC*TVEC(IX)
                        END IF
                      END IF
                    END DO
                  END DO
                END IF
              END IF
            END DO   !  K
C
C  20. H(KL|MN) <- D(IJ|KM)*(IJ|LN), L,N CORE
C
            ID = I - NCOR
            JD = J - NCOR
            IJD = IA(MAX(ID,JD)) + MIN(ID,JD)
            INJ = I.NE.J
            DO K = NCP1, NOCC
              KD = K - NCOR
              DO L = 1, NCOR
                KLR = IROT(K,L)
                IF (KLR.NE.0) THEN
                  DO M = NCP1, K
                    MD = M - NCOR
                    KMD = IA(MAX(KD,MD)) + MIN(KD,MD)
                    IJKMD = IA(MAX(IJD,KMD)) + MIN(IJD,KMD)
                    DIJKM = TPDM(IJKMD)
                    IF (INJ) DIJKM = DIJKM*TWO
                    MX = L
                    IF (K.NE.M) MX = NCOR
                    DO N = 1, MX
                      MNR = IROT(M,N)
                      IF (MNR.NE.0) THEN
                        IF (KLR.NE.MNR) THEN
                          LN = IA(MAX(L,N)) + MIN(L,N)
                          ERI = BUFF(LN)
                          FAC = DIJKM*ERI
                          IX = KLR + 1
                          JX = MNR + 1
                          PROD(IX) = PROD(IX) + FAC*TVEC(JX)
                          PROD(JX) = PROD(JX) + FAC*TVEC(IX)
                        END IF
                      END IF
                    END DO
                  END DO
                END IF
              END DO
            END DO
C
C  21. H(IK|LM) <- D(JL)*(IJ|KM), K CORE
C
            ID = I - NCOR
            JD = J - NCOR
            INJ = I.NE.J
            DO K = 1, NCOR
              IKR = IROT(I,K)
              IF (IKR.NE.0) THEN
                DO L = NCP1, NOCC
                  LD = L - NCOR
                  JLD = IA(MAX(JD,LD)) + MIN(JD,LD)
                  DJL = OPDM(JLD)
                  DO M = NCP1, L-1
                    LMR = IROT(L,M)
                    IF (LMR.NE.0) THEN
                      IF (IKR.NE.LMR) THEN
                        KM = IA(MAX(K,M)) + MIN(K,M)
                        ERI = BUFF(KM)
                        FAC = DJL*ERI
                        IX = IKR + 1
                        JX = LMR + 1
                        PROD(IX) = PROD(IX) + FAC*TVEC(JX)
                        PROD(JX) = PROD(JX) + FAC*TVEC(IX)
                      END IF
                    END IF
                  END DO
                END DO   !  L
              END IF
              IF (INJ) THEN
                JKR = IROT(J,K)
                IF (JKR.NE.0) THEN
                  DO L = NCP1, NOCC
                    LD = L - NCOR
                    ILD = IA(MAX(ID,LD)) + MIN(ID,LD)
                    DIL = OPDM(ILD)
                    DO M = NCP1, L-1
                      LMR = IROT(L,M)
                      IF (LMR.NE.0) THEN
                        IF (JKR.NE.LMR) THEN
                          KM = IA(MAX(K,M)) + MIN(K,M)
                          ERI = BUFF(KM)
                          FAC = DIL*ERI
                          IX = JKR + 1
                          JX = LMR + 1
                          PROD(IX) = PROD(IX) + FAC*TVEC(JX)
                          PROD(JX) = PROD(JX) + FAC*TVEC(IX)
                        END IF
                      END IF
                    END DO
                  END DO
                END IF
              END IF
            END DO   !  K
C
C  22. H(KL|MN) <- D(IJ|KM)*(IJ|LN), N CORE
C
            ID = I - NCOR
            JD = J - NCOR
            IJD = IA(MAX(ID,JD)) + MIN(ID,JD)
            INJ = I.NE.J
            DO K = NCP1, NOCC
              KD = K - NCOR
              DO L = NCP1, K-1
                KLR = IROT(K,L)
                IF (KLR.NE.0) THEN
                  DO M = NCP1, NOCC
                    MD = M - NCOR
                    KMD = IA(MAX(KD,MD)) + MIN(KD,MD)
                    IJKMD = IA(MAX(IJD,KMD)) + MIN(IJD,KMD)
                    DIJKM = TPDM(IJKMD)
                    IF (INJ) DIJKM = DIJKM*TWO
                    DO N = 1, NCOR
                      MNR = IROT(M,N)
                      IF (MNR.NE.0) THEN
                        IF (KLR.NE.MNR) THEN
                          LN = IA(MAX(L,N)) + MIN(L,N)
                          ERI = BUFF(LN)
                          FAC = DIJKM*ERI
                          IX = KLR + 1
                          JX = MNR + 1
                          PROD(IX) = PROD(IX) + FAC*TVEC(JX)
                          PROD(JX) = PROD(JX) + FAC*TVEC(IX)
                        END IF
                      END IF
                    END DO
                  END DO
                END IF
              END DO
            END DO
C
C  23. H(IM|LN) <- 2*D(JM|KN)*(IJ|KL)
C
            ID = I - NCOR
            JD = J - NCOR
            INJ = I.NE.J
            DO K = NCP1, NOCC
              KD = K - NCOR
              DO L = NCP1, I
                KL = IA(MAX(K,L)) + MIN(K,L)
                ERI = BUFF(KL)*TWO
                DO M = NCP1, I-1
                  IMR = IROT(I,M)
                  IF (IMR.NE.0) THEN
                    MD = M - NCOR
                    JMD = IA(MAX(JD,MD)) + MIN(JD,MD)
                    MX = L-1
                    IF (L.EQ.I) MX = M
                    DO N = NCP1, MX
                      LNR = IROT(L,N)
                      IF (LNR.NE.0) THEN
                        IF (IMR.NE.LNR) THEN
                          ND = N - NCOR
                          KND = IA(MAX(KD,ND)) + MIN(KD,ND)
                          JMKND = IA(MAX(JMD,KND)) + MIN(JMD,KND)
                          DJMKN = TPDM(JMKND)
                          FAC = DJMKN*ERI
                          IX = IMR + 1
                          JX = LNR + 1
                          PROD(IX) = PROD(IX) + FAC*TVEC(JX)
                          PROD(JX) = PROD(JX) + FAC*TVEC(IX)
                        END IF
                      END IF
                    END DO
                  END IF
                END DO
              END DO
            END DO
            IF (INJ) THEN
              DO K = NCP1, NOCC
                KD = K - NCOR
                DO L = NCP1, J
                  KL = IA(MAX(K,L)) + MIN(K,L)
                  ERI = BUFF(KL)*TWO
                  LEJ = L.EQ.J
                  DO M = NCP1, J-1
                    JMR = IROT(J,M)
                    IF (JMR.NE.0) THEN
                      MD  = M - NCOR
                      IMD = IA(MAX(ID,MD)) + MIN(ID,MD)
                      MX = L-1
                      IF (LEJ) MX = M
                      DO N = NCP1, MX
                        LNR = IROT(L,N)
                        IF (LNR.NE.0) THEN
                          IF (JMR.NE.LNR) THEN
                            ND = N - NCOR
                            KND = IA(MAX(KD,ND)) + MIN(KD,ND)
                            IMKND = IA(MAX(IMD,KND)) + MIN(IMD,KND)
                            DIMKN = TPDM(IMKND)
                            FAC = DIMKN*ERI
                            IX = JMR + 1
                            JX = LNR + 1
                            PROD(IX) = PROD(IX) + FAC*TVEC(JX)
                            PROD(JX) = PROD(JX) + FAC*TVEC(IX)
                          END IF
                        END IF
                      END DO
                    END IF
                  END DO
                END DO
              END DO
            END IF
C
C  24. H(KM|LN) <- D(IJ|MN)*(IJ|KL)
C
            ID = I - NCOR
            JD = J - NCOR
            IJD = IA(MAX(ID,JD)) + MIN(ID,JD)
            INJ = I.NE.J
            DO K = NCP1, NOCC
              DO L = NCP1, K
                KL = IA(MAX(K,L)) + MIN(K,L)
                ERI = BUFF(KL)
                IF (INJ) ERI = ERI*TWO
                KEL = K.EQ.L
                DO M = NCP1, K-1
                  KMR = IROT(K,M)
                  IF (KMR.NE.0) THEN
                    MD = M - NCOR
                    MX = L-1
                    IF (KEL) MX = M
                    DO N = NCP1, MX
                      LNR = IROT(L,N)
                      IF (LNR.NE.0) THEN
                        IF (KMR.NE.LNR) THEN
                          ND = N - NCOR
                          MND = IA(MAX(MD,ND)) + MIN(MD,ND)
                          IJMND = IA(MAX(IJD,MND)) + MIN(IJD,MND)
                          DIJMN = TPDM(IJMND)
                          FAC = DIJMN*ERI
                          IX = KMR + 1
                          JX = LNR + 1
                          PROD(IX) = PROD(IX) + FAC*TVEC(JX)
                          PROD(JX) = PROD(JX) + FAC*TVEC(IX)
                        END IF
                      END IF
                    END DO
                  END IF
                END DO
              END DO
            END DO
C
C  25. H(KL|MN) <- D(IJ|KM)*(IJ|LN)
C
            ID = I - NCOR
            JD = J - NCOR
            IJD = IA(MAX(ID,JD)) + MIN(ID,JD)
            INJ = I.NE.J
            DO K = NCP1, NOCC
              KD = K - NCOR
              DO L = NCP1, K-1
                KLR = IROT(K,L)
                IF (KLR.NE.0) THEN
                  DO M = NCP1, K
                    MD = M - NCOR
                    KMD = IA(MAX(KD,MD)) + MIN(KD,MD)
                    IJKMD = IA(MAX(IJD,KMD)) + MIN(IJD,KMD)
                    DIJKM = TPDM(IJKMD)
                    IF (INJ) DIJKM = DIJKM*TWO
                    MX = M-1
                    IF (K.EQ.M) MX = L
                    DO N = NCP1, MX
                      MNR = IROT(M,N)
                      IF (MNR.NE.0) THEN
                        IF (KLR.NE.MNR) THEN
                          LN = IA(MAX(L,N)) + MIN(L,N)
                          ERI = BUFF(LN)
                          FAC = DIJKM*ERI
                          IX = KLR + 1
                          JX = MNR + 1
                          PROD(IX) = PROD(IX) + FAC*TVEC(JX)
                          PROD(JX) = PROD(JX) + FAC*TVEC(IX)
                        END IF
                      END IF
                    END DO
                  END DO
                END IF
              END DO
            END DO
C
C  26. H(JK|LN) <- 2*D(IK|LM)*(IJ|MN)
C
            ID = I - NCOR
            JD = J - NCOR
            INJ = I.NE.J
            DO K = NCP1, NOCC
              KD = K - NCOR
              IF (J.LT.K) THEN
                JKR = IROT(K,J)
                IF (JKR.NE.0) THEN
                  IKD = IA(MAX(ID,KD)) + MIN(ID,KD)
                  DO L = NCP1, K
                    LD = L - NCOR
                    MX = L-1
                    IF (L.EQ.K) MX = J
                    DO M = NCP1, NOCC
                      MD = M - NCOR
                      LMD = IA(MAX(LD,MD)) + MIN(LD,MD)
                      IKLMD = IA(MAX(IKD,LMD)) + MIN(IKD,LMD)
                      DIKLM = TPDM(IKLMD)*TWO
                      DO N = NCP1, MX
                        LNR = IROT(L,N)
                        IF (LNR.NE.0) THEN
                          IF (JKR.NE.LNR) THEN
                            MN = IA(MAX(M,N)) + MIN(M,N)
                            ERI = BUFF(MN)
                            FAC = DIKLM*ERI
                            IX = JKR + 1
                            JX = LNR + 1
                            PROD(IX) = PROD(IX) + FAC*TVEC(JX)
                            PROD(JX) = PROD(JX) + FAC*TVEC(IX)
                          END IF
                        END IF
                      END DO
                    END DO
                  END DO
                END IF
              END IF
              IF (INJ.AND.(I.LT.K)) THEN
                IKR = IROT(I,K)
                IF (IKR.NE.0) THEN
                  JKD = IA(MAX(JD,KD)) + MIN(JD,KD)
                  DO L = NCP1, K
                    LD = L - NCOR
                    MX = L-1
                    IF (L.EQ.K) MX = I
                    DO M = NCP1, NOCC
                      MD = M - NCOR
                      LMD = IA(MAX(LD,MD)) + MIN(LD,MD)
                      JKLMD = IA(MAX(JKD,LMD)) + MIN(JKD,LMD)
                      DJKLM = TPDM(JKLMD)*TWO
                      DO N = NCP1, MX
                        LNR = IROT(L,N)
                        IF (LNR.NE.0) THEN
                          IF (IKR.NE.LNR) THEN
                            MN = IA(MAX(M,N)) + MIN(M,N)
                            ERI = BUFF(MN)
                            FAC = DJKLM*ERI
                            IX = IKR + 1
                            JX = LNR + 1
                            PROD(IX) = PROD(IX) + FAC*TVEC(JX)
                            PROD(JX) = PROD(JX) + FAC*TVEC(IX)
                          END IF
                        END IF
                      END DO
                    END DO
                  END DO
                END IF
              END IF
            END DO   !  K
C
C  27. H(KL|MN) <- -D(IJ|KN)*(IJ|LM)
C
            ID = I - NCOR
            JD = J - NCOR
            IJD = IA(MAX(ID,JD)) + MIN(ID,JD)
            INJ = I.NE.J
            DO K = NCP1, NOCC
              KD = K - NCOR
              DO L = NCP1, K-1
                KLR = IROT(K,L)
                IF (KLR.NE.0) THEN
                  DO M = NCP1, K
                    LM = IA(MAX(L,M)) + MIN(L,M)
                    ERI = BUFF(LM)
                    IF (INJ) ERI = ERI*TWO
                    MX = M-1
                    IF (K.EQ.M) MX = L
                    DO N = NCP1, MX
                      MNR = IROT(M,N)
                      IF (MNR.NE.0) THEN
                        IF (KLR.NE.MNR) THEN
                          ND = N - NCOR
                          KND = IA(MAX(KD,ND)) + MIN(KD,ND)
                          IJKND = IA(MAX(IJD,KND)) + MIN(IJD,KND)
                          DIJKN = TPDM(IJKND)
                          FAC = -DIJKN*ERI
                          IX = KLR + 1
                          JX = MNR + 1
                          PROD(IX) = PROD(IX) + FAC*TVEC(JX)
                          PROD(JX) = PROD(JX) + FAC*TVEC(IX)
                        END IF
                      END IF
                    END DO
                  END DO
                END IF
              END DO
            END DO
C
C  28.  H(JK|LN) <- -2*D(IK|MN)*(IJ|LM)
C
            ID = I - NCOR
            JD = J - NCOR
            INJ = I.NE.J
            DO K = NCP1, NOCC
              KD = K - NCOR
              IF (J.LT.K) THEN
                JKR = IROT(K,J)
                IF (JKR.NE.0) THEN
                  IKD = IA(MAX(ID,KD)) + MIN(ID,KD)
                  DO L = NCP1, K
                    MX = L-1
                    IF (K.EQ.L) MX = J
                    DO M = NCP1, NOCC
                      MD = M - NCOR
                      LM = IA(MAX(L,M)) + MIN(L,M)
                      ERI = BUFF(LM)*TWO
                      DO N = NCP1, MX
                        LNR = IROT(L,N)
                        IF (LNR.NE.0) THEN
                          IF (JKR.NE.LNR) THEN
                            ND = N - NCOR
                            MND = IA(MAX(MD,ND)) + MIN(MD,ND)
                            IKMND = IA(MAX(IKD,MND)) + MIN(IKD,MND)
                            DIKMN = TPDM(IKMND)
                            FAC = -DIKMN*ERI
                            IX = JKR + 1
                            JX = LNR + 1
                            PROD(IX) = PROD(IX) + FAC*TVEC(JX)
                            PROD(JX) = PROD(JX) + FAC*TVEC(IX)
                          END IF
                        END IF
                      END DO
                    END DO
                  END DO
                END IF
              END IF
              IF (INJ.AND.(I.LT.K)) THEN
                IKR = IROT(K,I)
                IF (IKR.NE.0) THEN
                  JKD = IA(MAX(JD,KD)) + MIN(JD,KD)
                  DO L = NCP1, K
                    MX = L-1
                    IF (K.EQ.L) MX = I
                    DO M = NCP1, NOCC
                      MD = M - NCOR
                      LM = IA(MAX(L,M)) + MIN(L,M)
                      ERI = BUFF(LM)*TWO
                      DO N = NCP1, MX
                        LNR = IROT(L,N)
                        IF (LNR.NE.0) THEN
                          IF (IKR.NE.LNR) THEN
                            ND = N - NCOR
                            MND = IA(MAX(MD,ND)) + MIN(MD,ND)
                            JKMND = IA(MAX(JKD,MND)) + MIN(JKD,MND)
                            DJKMN = TPDM(JKMND)
                            FAC = -DJKMN*ERI
                            IX = IKR + 1
                            JX = LNR + 1
                            PROD(IX) = PROD(IX) + FAC*TVEC(JX)
                            PROD(JX) = PROD(JX) + FAC*TVEC(IX)
                          END IF
                        END IF
                      END DO
                    END DO
                  END DO
                END IF
              END IF
            END DO   !  K
C
C  29. H(KL|MN) <- -D(IJ|LM)*(IJ|KN)
C
            ID = I - NCOR
            JD = J - NCOR
            IJD = IA(MAX(ID,JD)) + MIN(ID,JD)
            INJ = I.NE.J
            DO K = NCP1, NOCC
              DO L = NCP1, K-1
                LD = L - NCOR
                KLR = IROT(K,L)
                IF (KLR.NE.0) THEN
                  DO M = NCP1, K
                    MD = M - NCOR
                    LMD = IA(MAX(LD,MD)) + MIN(LD,MD)
                    IJLMD = IA(MAX(IJD,LMD)) + MIN(IJD,LMD)
                    DIJLM = TPDM(IJLMD)
                    IF (INJ) DIJLM = DIJLM*TWO
                    MX = M-1
                    IF (K.EQ.M) MX = L
                    DO N = NCP1, MX
                      MNR = IROT(M,N)
                      IF (MNR.NE.0) THEN
                        IF (KLR.NE.MNR) THEN
                          KN = IA(MAX(K,N)) + MIN(K,N)
                          ERI = BUFF(KN)
                          FAC = -DIJLM*ERI
                          IX = KLR + 1
                          JX = MNR + 1
                          PROD(IX) = PROD(IX) + FAC*TVEC(JX)
                          PROD(JX) = PROD(JX) + FAC*TVEC(IX)
                        END IF
                      END IF
                    END DO
                  END DO
                END IF
              END DO
            END DO
C
C  30. H(IK|LN) <- -2*D(JK|LM)*(IJ|MN)
C
            ID = I - NCOR
            JD = J - NCOR
            INJ = I.NE.J
            DO K = NCP1, I-1
              IKR = IROT(I,K)
              IF (IKR.NE.0) THEN
                KD = K - NCOR
                JKD = IA(MAX(JD,KD)) + MIN(JD,KD)
                DO L = NCP1, I
                  LD = L - NCOR
                  MX = L-1
                  IF (L.EQ.I) MX = K
                  DO M = NCP1, NOCC
                    MD = M - NCOR
                    LMD = IA(MAX(LD,MD)) + MIN(LD,MD)
                    JKLMD = IA(MAX(JKD,LMD)) + MIN(JKD,LMD)
                    DJKLM = TPDM(JKLMD)*TWO
                    DO N = NCP1, MX
                      LNR = IROT(L,N)
                      IF (LNR.NE.0) THEN
                        IF (IKR.NE.LNR) THEN
                          MN = IA(MAX(M,N)) + MIN(M,N)
                          ERI = BUFF(MN)
                          FAC = -DJKLM*ERI
                          IX = IKR + 1
                          JX = LNR + 1
                          PROD(IX) = PROD(IX) + FAC*TVEC(JX)
                          PROD(JX) = PROD(JX) + FAC*TVEC(IX)
                        END IF
                      END IF
                    END DO
                  END DO
                END DO
              END IF
            END DO
            IF (INJ) THEN
              DO K = NCP1, J-1
                JKR = IROT(J,K)
                IF (JKR.NE.0) THEN
                  KD = K - NCOR
                  IKD = IA(MAX(ID,KD)) + MIN(ID,KD)
                  DO L = NCP1, J
                    LD = L - NCOR
                    MX = L-1
                    IF (L.EQ.J) MX = K
                    DO M = NCP1, NOCC
                      MD = M - NCOR
                      LMD = IA(MAX(LD,MD)) + MIN(LD,MD)
                      IKLMD = IA(MAX(IKD,LMD)) + MIN(IKD,LMD)
                      DIKLM = TPDM(IKLMD)*TWO
                      DO N = NCP1, MX
                        LNR = IROT(L,N)
                        IF (LNR.NE.0) THEN
                          IF (JKR.NE.LNR) THEN
                            MN = IA(MAX(M,N)) + MIN(M,N)
                            ERI = BUFF(MN)
                            FAC = -DIKLM*ERI
                            IX = JKR + 1
                            JX = LNR + 1
                            PROD(IX) = PROD(IX) + FAC*TVEC(JX)
                            PROD(JX) = PROD(JX) + FAC*TVEC(IX)
                          END IF
                        END IF
                      END DO
                    END DO
                  END DO
                END IF
              END DO
            END IF
C
C  END DISTRIBUTION LOOPS
C
            CALL DDI_DLBNEXT(MYTASK)
          END IF  ! DLB
        END DO  ! J
      END DO  ! I
C
      IF (BENCH) THEN
        WRITE(ICHANL,9000) '(OO|OO)'
        CALL TIMIT(1)
      END IF
C
C ----------- CONTRIBUTIONS FROM (VO|OO) TYPE INTEGRALS -----------
C
C  NOTE: ALL THESE HESSIAN ELEMENTS ARE OFF-DIAGONAL
C
      CALL DDI_DISTRIB(D_VOOO,DDI_ME,ILO,IHI,JLO,JHI)
C
C  (ACTIVE-CORE|CORE-VIRTUAL)
C
      DO I = NCP1, NOCC
        DO J = 1, NCOR
          IJ = IA(MAX(I,J)) + MIN(I,J)
          DO K = 1, NCOR
            IJK = (IJ-1)*NOCC + K
            LOCTSK = LOCTSK + 1
            IF (LOCTSK.EQ.MYTASK) THEN
              CALL DDI_GET(D_VOOO,1,NVIR,IJK,IJK,BUFF)
C
C  31. H(AK|IJ) <- 8*(AK|IJ), J,K CORE
C
              IJR = IROT(I,J)
              IF (IJR.NE.0) THEN
                DO A = 1, NVIR
                  AR = A + NOCC
                  KAR = IROT(K,AR)
                  IF (KAR.NE.0) THEN
                    ERI = BUFF(A)
                    FAC = EIGHT*ERI
                    IX = IJR + 1
                    JX = KAR + 1
                    PROD(IX) = PROD(IX) + FAC*TVEC(JX)
                    PROD(JX) = PROD(JX) + FAC*TVEC(IX)
                  END IF
                END DO
              END IF
C
C  32. H(AJ|IK) <- -2*(AJ|IK), J,K CORE
C
              IKR = IROT(I,K)
              IF (IKR.NE.0) THEN
                DO A = 1, NVIR
                  AR = A + NOCC
                  JAR = IROT(J,AR)
                  IF (JAR.NE.0) THEN
                    ERI = BUFF(A)
                    FAC = -TWO*ERI
                    IX = IKR + 1
                    JX = JAR + 1
                    PROD(IX) = PROD(IX) + FAC*TVEC(JX)
                    PROD(JX) = PROD(JX) + FAC*TVEC(IX)
                  END IF
                END DO
              END IF
C
C  33. H(AK|LJ) <- -4*D(IL)*(AK|IJ), J,K CORE
C
              ID = I - NCOR
              DO A = 1, NVIR
                ERI = BUFF(A)
                AR = A + NOCC
                KAR = IROT(K,AR)
                IF (KAR.NE.0) THEN
                  DO L = NCP1, NOCC
                    JLR = IROT(L,J)
                    IF (JLR.NE.0) THEN
                      LD = L - NCOR
                      ILD = IA(MAX(ID,LD)) + MIN(ID,LD)
                      DIL = FOUR*OPDM(ILD)
                      FAC = -DIL*ERI
                      IX = KAR + 1
                      JX = JLR + 1
                      PROD(IX) = PROD(IX) + FAC*TVEC(JX)
                      PROD(JX) = PROD(JX) + FAC*TVEC(IX)
                    END IF
                  END DO
                END IF
C
C  34. H(AJ|LK) <- D(IL)*(AK|IJ), J,K CORE
C
                JAR = IROT(J,AR)
                IF (JAR.NE.0) THEN
                  DO L = NCP1, NOCC
                    KLR = IROT(L,K)
                    IF (KLR.NE.0) THEN
                      LD = L - NCOR
                      ILD = IA(MAX(ID,LD)) + MIN(ID,LD)
                      DIL = OPDM(ILD)
                      FAC = DIL*ERI
                      IX = JAR + 1
                      JX = KLR + 1
                      PROD(IX) = PROD(IX) + FAC*TVEC(JX)
                      PROD(JX) = PROD(JX) + FAC*TVEC(IX)
                    END IF
                  END DO
                END IF
              END DO   !  A
C
C  END DISTRIBUTION LOOPS
C
              CALL DDI_DLBNEXT(MYTASK)
            END IF  ! DLB
          END DO  ! K
        END DO  ! J
      END DO  ! I
C
C  (CORE-CORE|ACTIVE-VIRTUAL)
C
      DO I = 1, NCOR
        DO J = 1, I
          IJ = IA(MAX(I,J)) + MIN(I,J)
          DO K = NCP1, NOCC
            IJK = (IJ-1)*NOCC + K
            LOCTSK = LOCTSK + 1
            IF (LOCTSK.EQ.MYTASK) THEN
              CALL DDI_GET(D_VOOO,1,NVIR,IJK,IJK,BUFF)
              INJ = I.NE.J
C
C  35. H(AI|JK) <- -2*(AK|IJ), I,J CORE
C
              JKR = IROT(K,J)
              IF (JKR.NE.0) THEN
                DO A = 1, NVIR
                  AR = A + NOCC
                  IAR = IROT(I,AR)
                  IF (IAR.NE.0) THEN
                    ERI = BUFF(A)
                    FAC = -TWO*ERI
                    IX = JKR + 1
                    JX = IAR + 1
                    PROD(IX) = PROD(IX) + FAC*TVEC(JX)
                    PROD(JX) = PROD(JX) + FAC*TVEC(IX)
                  END IF
                END DO
              END IF
              IF (INJ) THEN
                IKR = IROT(K,I)
                IF (IKR.NE.0) THEN
                  DO A = 1, NVIR
                    AR = A + NOCC
                    JAR = IROT(J,AR)
                    IF (JAR.NE.0) THEN
                      ERI = BUFF(A)
                      FAC = -TWO*ERI
                      IX = IKR + 1
                      JX = JAR + 1
                      PROD(IX) = PROD(IX) + FAC*TVEC(JX)
                      PROD(JX) = PROD(JX) + FAC*TVEC(IX)
                    END IF
                  END DO
                END IF
              END IF
C
C  36. H(AI|LJ) <- D(KL)*(AK|IJ), I,J CORE
C
              KD = K - NCOR
              DO A = 1, NVIR
                ERI = BUFF(A)
                AR = A + NOCC
                IAR = IROT(I,AR)
                IF (IAR.NE.0) THEN
                  DO L = NCP1, NOCC
                    LJR = IROT(L,J)
                    IF (LJR.NE.0) THEN
                      LD = L - NCOR
                      KLD = IA(MAX(KD,LD)) + MIN(KD,LD)
                      DKL = OPDM(KLD)
                      FAC = DKL*ERI
                      IX = IAR + 1
                      JX = LJR + 1
                      PROD(IX) = PROD(IX) + FAC*TVEC(JX)
                      PROD(JX) = PROD(JX) + FAC*TVEC(IX)
                    END IF
                  END DO
                END IF
                IF (INJ) THEN
                  JAR = IROT(J,AR)
                  IF (JAR.NE.0) THEN
                    DO L = NCP1, NOCC
                      LIR = IROT(L,I)
                      IF (LIR.NE.0) THEN
                        LD = L - NCOR
                        KLD = IA(MAX(KD,LD)) + MIN(KD,LD)
                        DKL = OPDM(KLD)
                        FAC = DKL*ERI
                        IX = JAR + 1
                        JX = LIR + 1
                        PROD(IX) = PROD(IX) + FAC*TVEC(JX)
                        PROD(JX) = PROD(JX) + FAC*TVEC(IX)
                      END IF
                    END DO
                  END IF
                END IF
              END DO   !  A
C
C  END DISTRIBUTION LOOPS
C
              CALL DDI_DLBNEXT(MYTASK)
            END IF  ! DLB
          END DO  ! K
        END DO  ! J
      END DO  ! I
C
C  (ACTIVE-ACTIVE|CORE-VIRTUAL)
C
      DO I = NCP1, NOCC
        DO J = NCP1, I
          IJ = IA(MAX(I,J)) + MIN(I,J)
          DO K = 1, NCOR
            IJK = (IJ-1)*NOCC + K
            LOCTSK = LOCTSK + 1
            IF (LOCTSK.EQ.MYTASK) THEN
              CALL DDI_GET(D_VOOO,1,NVIR,IJK,IJK,BUFF)
C
C  37. H(AK|IL) <- 4*D(JL)*(AK|IJ), K CORE
C
              INJ = I.NE.J
              ID = I - NCOR
              JD = J - NCOR
              DO A = 1, NVIR
                ERI = BUFF(A)*FOUR
                AR = A + NOCC
                KAR = IROT(K,AR)
                IF (KAR.NE.0) THEN
                  DO L = NCP1, I-1
                    ILR = IROT(I,L)
                    IF (ILR.NE.0) THEN
                      LD = L - NCOR
                      JLD = IA(MAX(JD,LD)) + MIN(JD,LD)
                      DJL = OPDM(JLD)
                      FAC = DJL*ERI
                      IX = KAR + 1
                      JX = ILR + 1
                      PROD(IX) = PROD(IX) + FAC*TVEC(JX)
                      PROD(JX) = PROD(JX) + FAC*TVEC(IX)
                    END IF
                  END DO
                  IF (INJ) THEN
                    DO L = NCP1, J-1
                      JLR = IROT(J,L)
                      IF (JLR.NE.0) THEN
                        LD = L - NCOR
                        ILD = IA(MAX(ID,LD)) + MIN(ID,LD)
                        DIL = OPDM(ILD)
                        FAC = DIL*ERI
                        IX = KAR + 1
                        JX = JLR + 1
                        PROD(IX) = PROD(IX) + FAC*TVEC(JX)
                        PROD(JX) = PROD(JX) + FAC*TVEC(IX)
                      END IF
                    END DO
                  END IF
                END IF
              END DO   !  A
C
C  38. H(AK|JL) <- -4*D(IL)*(AK|IJ), K CORE
C
              INJ = I.NE.J
              ID = I - NCOR
              JD = J - NCOR
              DO A = 1, NVIR
                ERI = BUFF(A)*FOUR
                AR = A + NOCC
                KAR = IROT(K,AR)
                IF (KAR.NE.0) THEN
                  DO L = J+1, NOCC
                    JLR = IROT(L,J)
                    IF (JLR.NE.0) THEN
                      LD = L - NCOR
                      ILD = IA(MAX(ID,LD)) + MIN(ID,LD)
                      DIL = OPDM(ILD)
                      FAC = -DIL*ERI
                      IX = KAR + 1
                      JX = JLR + 1
                      PROD(IX) = PROD(IX) + FAC*TVEC(JX)
                      PROD(JX) = PROD(JX) + FAC*TVEC(IX)
                    END IF
                  END DO
                  IF (INJ) THEN
                    DO L = I+1, NOCC
                      ILR = IROT(L,I)
                      IF (ILR.NE.0) THEN
                        LD = L - NCOR
                        JLD = IA(MAX(JD,LD)) + MIN(JD,LD)
                        DJL = OPDM(JLD)
                        FAC = -DJL*ERI
                        IX = KAR + 1
                        JX = ILR + 1
                        PROD(IX) = PROD(IX) + FAC*TVEC(JX)
                        PROD(JX) = PROD(JX) + FAC*TVEC(IX)
                      END IF
                    END DO
                  END IF
                END IF
              END DO   !  A
C
C  39. H(AL|IK) <- -D(JL)*(AK|IJ), K CORE
C
              INJ = I.NE.J
              ID = I - NCOR
              JD = J - NCOR
              IKR = IROT(I,K)
              IF (IKR.NE.0) THEN
                DO A = 1, NVIR
                  ERI = BUFF(A)
                  AR = A + NOCC
                  DO L = NCP1, NOCC
                    LAR = IROT(L,AR)
                    IF (LAR.NE.0) THEN
                      LD = L - NCOR
                      JLD = IA(MAX(JD,LD)) + MIN(JD,LD)
                      DJL = OPDM(JLD)
                      FAC = -DJL*ERI
                      IX = IKR + 1
                      JX = LAR + 1
                      PROD(IX) = PROD(IX) + FAC*TVEC(JX)
                      PROD(JX) = PROD(JX) + FAC*TVEC(IX)
                    END IF
                  END DO
                END DO
              END IF
              IF (INJ) THEN
                JKR = IROT(J,K)
                IF (JKR.NE.0) THEN
                  DO A = 1, NVIR
                    ERI = BUFF(A)
                    AR = A + NOCC
                    DO L = NCP1, NOCC
                      LAR = IROT(L,AR)
                      IF (LAR.NE.0) THEN
                        LD = L - NCOR
                        ILD = IA(MAX(ID,LD)) + MIN(ID,LD)
                        DIL = OPDM(ILD)
                        FAC = -DIL*ERI
                        IX = JKR + 1
                        JX = LAR + 1
                        PROD(IX) = PROD(IX) + FAC*TVEC(JX)
                        PROD(JX) = PROD(JX) + FAC*TVEC(IX)
                      END IF
                    END DO
                  END DO
                END IF
              END IF
C
C  40. H(AL|KM) <- -D(IJ|LM)*(AK|IJ), K CORE
C
              INJ = I.NE.J
              ID = I - NCOR
              JD = J - NCOR
              IJD = IA(MAX(ID,JD)) + MIN(ID,JD)
              DO A = 1, NVIR
                ERI = BUFF(A)
                IF (INJ) ERI = ERI*TWO
                AR = A + NOCC
                DO L = NCP1, NOCC
                  LD = L - NCOR
                  LAR = IROT(L,AR)
                  IF (LAR.NE.0) THEN
                    DO M = NCP1, NOCC
                      KMR = IROT(M,K)
                      IF (KMR.NE.0) THEN
                        MD = M - NCOR
                        LMD = IA(MAX(LD,MD)) + MIN(LD,MD)
                        IJLMD = IA(MAX(IJD,LMD)) + MIN(IJD,LMD)
                        DIJLM = TPDM(IJLMD)
                        FAC = -DIJLM*ERI
                        IX = LAR + 1
                        JX = KMR + 1
                        PROD(IX) = PROD(IX) + FAC*TVEC(JX)
                        PROD(JX) = PROD(JX) + FAC*TVEC(IX)
                      END IF
                    END DO
                  END IF
                END DO
              END DO
C
C  END DISTRIBUTION LOOPS
C
              CALL DDI_DLBNEXT(MYTASK)
            END IF  ! DLB
          END DO  ! K
        END DO  ! J
      END DO  ! I
C
C  (ACTIVE-CORE|ACTIVE-VIRTUAL)
C
      DO I = NCP1, NOCC
        DO J = 1, NCOR
          IJ = IA(MAX(I,J)) + MIN(I,J)
          DO K = NCP1, NOCC
            IJK = (IJ-1)*NOCC + K
            LOCTSK = LOCTSK + 1
            IF (LOCTSK.EQ.MYTASK) THEN
              CALL DDI_GET(D_VOOO,1,NVIR,IJK,IJK,BUFF)
C
C  41. H(AJ|IL) <- -D(KL)*(AK|IJ), J CORE
C
              ID = I - NCOR
              KD = K - NCOR
              DO A = 1, NVIR
                ERI = BUFF(A)
                AR = A + NOCC
                JAR = IROT(J,AR)
                IF (JAR.NE.0) THEN
                  DO L = NCP1, I-1
                    ILR = IROT(I,L)
                    IF (ILR.NE.0) THEN
                      LD = L - NCOR
                      KLD = IA(MAX(KD,LD)) + MIN(KD,LD)
                      DKL = OPDM(KLD)
                      FAC = -DKL*ERI
                      IX = JAR + 1
                      JX = ILR + 1
                      PROD(IX) = PROD(IX) + FAC*TVEC(JX)
                      PROD(JX) = PROD(JX) + FAC*TVEC(IX)
                    END IF
                  END DO
C
C  42. H(AJ|KL) <- -D(IL)*(AK|IJ), J CORE
C
                  DO L = NCP1, K-1
                    KLR = IROT(K,L)
                    IF (KLR.NE.0) THEN
                      LD = L - NCOR
                      ILD = IA(MAX(ID,LD)) + MIN(ID,LD)
                      DIL = OPDM(ILD)
                      FAC = -DIL*ERI
                      IX = JAR + 1
                      JX = KLR + 1
                      PROD(IX) = PROD(IX) + FAC*TVEC(JX)
                      PROD(JX) = PROD(JX) + FAC*TVEC(IX)
                    END IF
                  END DO
                END IF
              END DO   !  A
C
C  43. H(AJ|IL) <- D(KL)*(AK|IJ), J CORE
C
              ID = I - NCOR
              KD = K - NCOR
              DO A = 1, NVIR
                ERI = BUFF(A)
                AR = A + NOCC
                JAR = IROT(J,AR)
                IF (JAR.NE.0) THEN
                  DO L = I+1, NOCC
                    ILR = IROT(L,I)
                    IF (ILR.NE.0) THEN
                      LD = L - NCOR
                      KLD = IA(MAX(KD,LD)) + MIN(KD,LD)
                      DKL = OPDM(KLD)
                      FAC = DKL*ERI
                      IX = JAR + 1
                      JX = ILR + 1
                      PROD(IX) = PROD(IX) + FAC*TVEC(JX)
                      PROD(JX) = PROD(JX) + FAC*TVEC(IX)
                    END IF
                  END DO
C
C  44. H(AJ|KL) <- D(IL)*(AK|IJ), J CORE
C
                  DO L = K+1, NOCC
                    KLR = IROT(L,K)
                    IF (KLR.NE.0) THEN
                      LD = L - NCOR
                      ILD = IA(MAX(ID,LD)) + MIN(ID,LD)
                      DIL = OPDM(ILD)
                      FAC = DIL*ERI
                      IX = JAR + 1
                      JX = KLR + 1
                      PROD(IX) = PROD(IX) + FAC*TVEC(JX)
                      PROD(JX) = PROD(JX) + FAC*TVEC(IX)
                    END IF
                  END DO
                END IF
              END DO   !  A
C
C  45. H(AL|IJ) <- 4*D(KL)*(AK|IJ), J CORE
C
              KD = K - NCOR
              IJR = IROT(I,J)
              IF (IJR.NE.0) THEN
                DO A = 1, NVIR
                  ERI = BUFF(A)*FOUR
                  AR = A + NOCC
                  DO L = NCP1, NOCC
                    LAR = IROT(L,AR)
                    IF (LAR.NE.0) THEN
                      LD = L - NCOR
                      KLD = IA(MAX(KD,LD)) + MIN(KD,LD)
                      DKL = OPDM(KLD)
                      FAC = DKL*ERI
                      IX = IJR + 1
                      JX = LAR + 1
                      PROD(IX) = PROD(IX) + FAC*TVEC(JX)
                      PROD(JX) = PROD(JX) + FAC*TVEC(IX)
                    END IF
                  END DO
                END DO
              END IF
C
C  46. H(AL|KJ) <- -D(IL)*(AK|IJ), J CORE
C
              ID = I - NCOR
              KJR = IROT(K,J)
              IF (KJR.NE.0) THEN
                DO A = 1, NVIR
                  ERI = BUFF(A)
                  AR = A + NOCC
                  DO L = NCP1, NOCC
                    LAR = IROT(L,AR)
                    IF (LAR.NE.0) THEN
                      LD = L - NCOR
                      ILD = IA(MAX(ID,LD)) + MIN(ID,LD)
                      DIL = OPDM(ILD)
                      FAC = -DIL*ERI
                      IX = KJR + 1
                      JX = LAR + 1
                      PROD(IX) = PROD(IX) + FAC*TVEC(JX)
                      PROD(JX) = PROD(JX) + FAC*TVEC(IX)
                    END IF
                  END DO
                END DO
              END IF
C
C  47. H(AL|JM) <- -2*D(KL|IM)*(AK|IJ), J CORE
C
              ID = I - NCOR
              KD = K - NCOR
              DO A = 1, NVIR
                ERI = BUFF(A)*TWO
                AR = A + NOCC
                DO L = NCP1, NOCC
                  LAR = IROT(L,AR)
                  IF (LAR.NE.0) THEN
                    LD = L - NCOR
                    KLD = IA(MAX(KD,LD)) + MIN(KD,LD)
                    DO M = NCP1, NOCC
                      JMR = IROT(M,J)
                      IF (JMR.NE.0) THEN
                        MD = M - NCOR
                        IMD = IA(MAX(ID,MD)) + MIN(ID,MD)
                        KLIMD = IA(MAX(KLD,IMD)) + MIN(KLD,IMD)
                        DKLIM = TPDM(KLIMD)
                        FAC = -DKLIM*ERI
                        IX = LAR + 1
                        JX = JMR + 1
                        PROD(IX) = PROD(IX) + FAC*TVEC(JX)
                        PROD(JX) = PROD(JX) + FAC*TVEC(IX)
                      END IF
                    END DO
                  END IF
                END DO
              END DO
C
C  END DISTRIBUTION LOOPS
C
              CALL DDI_DLBNEXT(MYTASK)
            END IF  ! DLB
          END DO  ! K
        END DO  ! J
      END DO  ! I
C
C  (ACTIVE-ACTIVE|ACTIVE-VIRTUAL)
C
      DO I = NCP1, NOCC
        DO J = NCP1, I
          IJ = IA(MAX(I,J)) + MIN(I,J)
          DO K = NCP1, NOCC
            IJK = (IJ-1)*NOCC + K
            LOCTSK = LOCTSK + 1
            IF (LOCTSK.EQ.MYTASK) THEN
              CALL DDI_GET(D_VOOO,1,NVIR,IJK,IJK,BUFF)
C
C  48. H(AL|IM) <- 2*D(KL|JM)*(AK|IJ)
C
              INJ = I.NE.J
              ID = I - NCOR
              JD = J - NCOR
              KD = K - NCOR
              DO A = 1, NVIR
                ERI = BUFF(A)*TWO
                AR = A + NOCC
                DO L = NCP1, NOCC
                  LAR = IROT(L,AR)
                  IF (LAR.NE.0) THEN
                    LD = L - NCOR
                    KLD = IA(MAX(KD,LD)) + MIN(KD,LD)
                    DO M = NCP1, I-1
                      IMR = IROT(I,M)
                      IF (IMR.NE.0) THEN
                        MD = M - NCOR
                        JMD = IA(MAX(JD,MD)) + MIN(JD,MD)
                        KLJMD = IA(MAX(KLD,JMD)) + MIN(KLD,JMD)
                        DKLJM = TPDM(KLJMD)
                        FAC = DKLJM*ERI
                        IX = LAR + 1
                        JX = IMR + 1
                        PROD(IX) = PROD(IX) + FAC*TVEC(JX)
                        PROD(JX) = PROD(JX) + FAC*TVEC(IX)
                      END IF
                    END DO
C
C  49. H(AL|JM) <- -2*D(KL|IM)*(AK|IJ)
C
                    DO M = J+1, NOCC
                      JMR = IROT(J,M)
                      IF (JMR.NE.0) THEN
                        MD = M - NCOR
                        IMD = IA(MAX(ID,MD)) + MIN(ID,MD)
                        KLIMD = IA(MAX(KLD,IMD)) + MIN(KLD,IMD)
                        DKLIM = TPDM(KLIMD)
                        FAC = -DKLIM*ERI
                        IX = LAR + 1
                        JX = JMR + 1
                        PROD(IX) = PROD(IX) + FAC*TVEC(JX)
                        PROD(JX) = PROD(JX) + FAC*TVEC(IX)
                      END IF
                    END DO
C
C  50. H(AL|JM) <- 2*D(KL|IM)*(AK|IJ)
C
                    IF (INJ) THEN
                      DO M = NCP1, J-1
                        JMR = IROT(J,M)
                        IF (JMR.NE.0) THEN
                          MD = M - NCOR
                          IMD = IA(MAX(ID,MD)) + MIN(ID,MD)
                          KLIMD = IA(MAX(KLD,IMD)) + MIN(KLD,IMD)
                          DKLIM = TPDM(KLIMD)
                          FAC = DKLIM*ERI
                          IX = LAR + 1
                          JX = JMR + 1
                          PROD(IX) = PROD(IX) + FAC*TVEC(JX)
                          PROD(JX) = PROD(JX) + FAC*TVEC(IX)
                        END IF
                      END DO
C
C  51. H(AL|MI) <- -2*D(KL|JM)*(AK|IJ)
C
                      DO M = I+1, NOCC
                        IMR = IROT(I,M)
                        IF (IMR.NE.0) THEN
                          MD = M - NCOR
                          JMD = IA(MAX(JD,MD)) + MIN(JD,MD)
                          KLJMD = IA(MAX(KLD,JMD)) + MIN(KLD,JMD)
                          DKLJM = TPDM(KLJMD)
                          FAC = -DKLJM*ERI
                          IX = LAR + 1
                          JX = IMR + 1
                          PROD(IX) = PROD(IX) + FAC*TVEC(JX)
                          PROD(JX) = PROD(JX) + FAC*TVEC(IX)
                        END IF
                      END DO
                    END IF
                  END IF
                END DO   !  L
              END DO   !  A
C
C  52. H(AL|KM) <- D(LM|IJ)*(AK|IJ)
C
              INJ = I.NE.J
              ID = I - NCOR
              JD = J - NCOR
              IJD = IA(MAX(ID,JD)) + MIN(ID,JD)
              DO A = 1, NVIR
                ERI = BUFF(A)
                IF (INJ) ERI = ERI*TWO
                AR = A + NOCC
                DO L = NCP1, NOCC
                  LAR = IROT(L,AR)
                  IF (LAR.NE.0) THEN
                    LD = L - NCOR
                    DO M = NCP1, K-1
                      KMR = IROT(K,M)
                      IF (KMR.NE.0) THEN
                        MD = M - NCOR
                        LMD = IA(MAX(LD,MD)) + MIN(LD,MD)
                        IJLMD = IA(MAX(IJD,LMD)) + MIN(IJD,LMD)
                        DIJLM = TPDM(IJLMD)
                        FAC = DIJLM*ERI
                        IX = LAR + 1
                        JX = KMR + 1
                        PROD(IX) = PROD(IX) + FAC*TVEC(JX)
                        PROD(JX) = PROD(JX) + FAC*TVEC(IX)
                      END IF
                    END DO
C
C  53. H(AL|MK) <- -D(LM|IJ)*(AK|IJ)
C
                    DO M = K+1, NOCC
                      KMR = IROT(K,M)
                      IF (KMR.NE.0) THEN
                        MD = M - NCOR
                        LMD = IA(MAX(LD,MD)) + MIN(LD,MD)
                        IJLMD = IA(MAX(IJD,LMD)) + MIN(IJD,LMD)
                        DIJLM = TPDM(IJLMD)
                        FAC = -DIJLM*ERI
                        IX = LAR + 1
                        JX = KMR + 1
                        PROD(IX) = PROD(IX) + FAC*TVEC(JX)
                        PROD(JX) = PROD(JX) + FAC*TVEC(IX)
                      END IF
                    END DO
                  END IF
                END DO    !  L
              END DO    !  A
C
C  END DISTRIBUTION LOOPS
C
              CALL DDI_DLBNEXT(MYTASK)
            END IF  ! DLB
          END DO  ! K
        END DO  ! J
      END DO  ! I
C
      IF (BENCH) THEN
        WRITE(ICHANL,9000) '(VO|OO)'
        CALL TIMIT(1)
      END IF
C
C ----------- CONTRIBUTIONS FROM (VV|OO) TYPE INTEGRALS -----------
C
      CALL DDI_DISTRIB(D_VVOO,DDI_ME,ILO,IHI,JLO,JHI)
C
C  (VIRTUAL-VIRTUAL|CORE-CORE)
C
      DO I = 1, NCOR
        DO J = 1, I
          IJ = IA(MAX(I,J)) + MIN(I,J)
          LOCTSK = LOCTSK + 1
          IF (LOCTSK.EQ.MYTASK) THEN
            CALL DDI_GET(D_VVOO,1,NVTR,IJ,IJ,BUFF)
C
C  54. H(AI|BJ) <- -2*(AB|IJ), I,J CORE
C
            DO A = 1, NVIR
              AR = A + NOCC
              IAR = IROT(I,AR)
              IF (IAR.NE.0) THEN
                DO B = 1, A
                  BR = B + NOCC
                  JBR = IROT(J,BR)
                  IF (JBR.NE.0) THEN
                    IF (IAR.NE.JBR) THEN
                      AB = IA(MAX(A,B)) + MIN(A,B)
                      ERI = BUFF(AB)
                      FAC = -TWO*ERI
                      IX = IAR + 1
                      JX = JBR + 1
                      PROD(IX) = PROD(IX) + FAC*TVEC(JX)
                      PROD(JX) = PROD(JX) + FAC*TVEC(IX)
                    END IF
                  END IF
                END DO
              END IF
            END DO
            IF (I.NE.J) THEN
              DO A = 1, NVIR
                AR = A + NOCC
                JAR = IROT(J,AR)
                IF (JAR.NE.0) THEN
                  DO B = 1, A-1
                    BR = B + NOCC
                    IBR = IROT(I,BR)
                    IF (IBR.NE.0) THEN
                      IF (JAR.NE.IBR) THEN
                        AB = IA(MAX(A,B)) + MIN(A,B)
                        ERI = BUFF(AB)
                        FAC = -TWO*ERI
                        IX = JAR + 1
                        JX = IBR + 1
                        PROD(IX) = PROD(IX) + FAC*TVEC(JX)
                        PROD(JX) = PROD(JX) + FAC*TVEC(IX)
                      END IF
                    END IF
                  END DO
                END IF
              END DO
            END IF
C
C  END DISTRIBUTION LOOPS
C
            CALL DDI_DLBNEXT(MYTASK)
          END IF  ! DLB
        END DO  ! J
      END DO  ! I
C
C  (VIRTUAL-VIRTUAL|ACTIVE-CORE)
C
      DO I = NCP1, NOCC
        DO J = 1, NCOR
          IJ = IA(MAX(I,J)) + MIN(I,J)
          LOCTSK = LOCTSK + 1
          IF (LOCTSK.EQ.MYTASK) THEN
            CALL DDI_GET(D_VVOO,1,NVTR,IJ,IJ,BUFF)
C
C  55. H(AK|BJ) <- -D(IK)*(AB|IJ), J CORE
C
            ID = I - NCOR
            DO A = 1, NVIR
              AR = A + NOCC
              DO B = 1, A
                BR = B + NOCC
                JBR = IROT(J,BR)
                IF (JBR.NE.0) THEN
                  AB = IA(MAX(A,B)) + MIN(A,B)
                  ERI = BUFF(AB)
                  DO K = NCP1, NOCC
                    KAR = IROT(K,AR)
                    IF (KAR.NE.0) THEN
                      IF (JBR.NE.KAR) THEN
                        KD = K - NCOR
                        IKD = IA(MAX(ID,KD)) + MIN(ID,KD)
                        DIK = OPDM(IKD)
                        FAC = -DIK*ERI
                        IX = JBR + 1
                        JX = KAR + 1
                        PROD(IX) = PROD(IX) + FAC*TVEC(JX)
                        PROD(JX) = PROD(JX) + FAC*TVEC(IX)
                      END IF
                    END IF
                  END DO
                END IF
              END DO
            END DO
            DO A = 1, NVIR
              AR = A + NOCC
              JAR = IROT(J,AR)
              IF (JAR.NE.0) THEN
                DO B = 1, A-1
                  BR = B + NOCC
                  AB = IA(MAX(A,B)) + MIN(A,B)
                  ERI = BUFF(AB)
                  DO K = NCP1, NOCC
                    KBR = IROT(K,BR)
                    IF (KBR.NE.0) THEN
                      IF (JAR.NE.KBR) THEN
                        KD = K - NCOR
                        IKD = IA(MAX(ID,KD)) + MIN(ID,KD)
                        DIK = OPDM(IKD)
                        FAC = -DIK*ERI
                        IX = JAR + 1
                        JX = KBR + 1
                        PROD(IX) = PROD(IX) + FAC*TVEC(JX)
                        PROD(JX) = PROD(JX) + FAC*TVEC(IX)
                      END IF
                    END IF
                  END DO
                END DO
              END IF
            END DO
C
C  END DISTRIBUTION LOOPS
C
            CALL DDI_DLBNEXT(MYTASK)
          END IF  ! DLB
        END DO  ! J
      END DO  ! I
C
C  (VIRTUAL-VIRTUAL|ACTIVE-ACTIVE)
C
      DO I = NCP1, NOCC
        DO J = NCP1, I
          IJ = IA(MAX(I,J)) + MIN(I,J)
          LOCTSK = LOCTSK + 1
          IF (LOCTSK.EQ.MYTASK) THEN
            CALL DDI_GET(D_VVOO,1,NVTR,IJ,IJ,BUFF)
C
C  56. H(AK|BL) <- D(IJ|KL)*(AB|IJ)
C
            INJ = I.NE.J
            ID = I - NCOR
            JD = J - NCOR
            IJD = IA(MAX(ID,JD)) + MIN(ID,JD)
            DO A = 1, NVIR
              AR = A + NOCC
              DO B = 1, A
                BR = B + NOCC
                AB = IA(MAX(A,B)) + MIN(A,B)
                ERI = BUFF(AB)
                IF (INJ) ERI = ERI*TWO
                ANB = A.NE.B
                DO K = NCP1, NOCC
                  KAR = IROT(K,AR)
                  IF (KAR.NE.0) THEN
                    KD = K - NCOR
                    MX = K
                    IF (ANB) MX = NOCC
                    DO L = NCP1, MX
                      LBR = IROT(L,BR)
                      IF (LBR.NE.0) THEN
                        IF (KAR.NE.LBR) THEN
                          LD = L - NCOR
                          KLD = IA(MAX(KD,LD)) + MIN(KD,LD)
                          IJKLD = IA(MAX(IJD,KLD)) + MIN(IJD,KLD)
                          DIJKL = TPDM(IJKLD)
                          FAC = DIJKL*ERI
                          IX = KAR + 1
                          JX = LBR + 1
                          PROD(IX) = PROD(IX) + FAC*TVEC(JX)
                          PROD(JX) = PROD(JX) + FAC*TVEC(IX)
                        END IF
                      END IF
                    END DO
                  END IF
                END DO
              END DO
            END DO
C
C  END DISTRIBUTION LOOPS
C
            CALL DDI_DLBNEXT(MYTASK)
          END IF  ! DLB
        END DO  ! J
      END DO  ! I
C
      IF (BENCH) THEN
        WRITE(ICHANL,9000) '(VV|OO)'
        CALL TIMIT(1)
      END IF
C
C ----------- CONTRIBUTIONS FROM (VO|VO) TYPE INTEGRALS -----------
C
      CALL DDI_DISTRIB(D_VOVO,DDI_ME,ILO,IHI,JLO,JHI)
C
C  (VIRTUAL-CORE|VIRTUAL-CORE)
C
      DO I = 1, NCOR
        DO J = 1, I
          IJ = IA(MAX(I,J)) + MIN(I,J)
          LOCTSK = LOCTSK + 1
          IF (LOCTSK.EQ.MYTASK) THEN
            CALL DDI_GET(D_VOVO,1,NVSQ,IJ,IJ,BUFF)
C
C  57. H(AI|BJ) <- 8*(AI|BJ) -2*(AJ|BI), I,J CORE
C
            INJ = I.NE.J
            DO A = 1, NVIR
              AR = A + NOCC
              IAR = IROT(I,AR)
              IF (IAR.NE.0) THEN
                MX = A
                IF (INJ) MX = NVIR
                DO B = 1, MX
                  BR = B + NOCC
                  JBR = IROT(J,BR)
                  IF (JBR.NE.0) THEN
                    IF (IAR.NE.JBR) THEN
                      AB = (B-1)*NVIR + A
                      BA = (A-1)*NVIR + B
                      ERI1 = BUFF(AB)*EIGHT
                      ERI2 = BUFF(BA)*TWO
                      FAC = ERI1 - ERI2
                      IX = IAR + 1
                      JX = JBR + 1
                      PROD(IX) = PROD(IX) + FAC*TVEC(JX)
                      PROD(JX) = PROD(JX) + FAC*TVEC(IX)
                    END IF
                  END IF
                END DO
              END IF
            END DO
C
C  END DISTRIBUTION LOOPS
C
            CALL DDI_DLBNEXT(MYTASK)
          END IF  ! DLB
        END DO  ! J
      END DO  ! I
C
C  (VIRTUAL-ACTIVE|VIRTUAL-CORE)
C
      DO I = NCP1, NOCC
        DO J = 1, NCOR
          IJ = IA(MAX(I,J)) + MIN(I,J)
          LOCTSK = LOCTSK + 1
          IF (LOCTSK.EQ.MYTASK) THEN
            CALL DDI_GET(D_VOVO,1,NVSQ,IJ,IJ,BUFF)
C
C  58. H(AK|BJ) <- 4*D(IK)*(AI|BJ) -D(IK)*(AJ|BI), J CORE
C
            ID = I - NCOR
            DO A = 1, NVIR
              AR = A + NOCC
              DO B = 1, NVIR
                BR = B + NOCC
                JBR = IROT(J,BR)
                IF (JBR.NE.0) THEN
                  AB = (B-1)*NVIR + A
                  BA = (A-1)*NVIR + B
                  ERI1 = BUFF(AB)*FOUR
                  ERI2 = BUFF(BA)
                  DO K = NCP1, NOCC
                    KAR = IROT(K,AR)
                    IF (KAR.NE.0) THEN
                      IF (JBR.NE.KAR) THEN
                        KD = K - NCOR
                        IKD = IA(MAX(ID,KD)) + MIN(ID,KD)
                        DIK = OPDM(IKD)
                        FAC = DIK*(ERI1-ERI2)
                        IX = JBR + 1
                        JX = KAR + 1
                        PROD(IX) = PROD(IX) + FAC*TVEC(JX)
                        PROD(JX) = PROD(JX) + FAC*TVEC(IX)
                      END IF
                    END IF
                  END DO
                END IF
              END DO
            END DO
C
C  END DISTRIBUTION LOOPS
C
            CALL DDI_DLBNEXT(MYTASK)
          END IF  ! DLB
        END DO  ! J
      END DO  ! I
C
C  (VIRTUAL-ACTIVE|VIRTUAL-ACTIVE)
C
      DO I = NCP1, NOCC
        DO J = NCP1, I
          IJ = IA(MAX(I,J)) + MIN(I,J)
          LOCTSK = LOCTSK + 1
          IF (LOCTSK.EQ.MYTASK) THEN
            CALL DDI_GET(D_VOVO,1,NVSQ,IJ,IJ,BUFF)
C
C  59. H(AK|BL) <- 2*D(IK|JL)*(AI|BJ)
C
            ID = I - NCOR
            JD = J - NCOR
            DO A = 1, NVIR
              AR = A + NOCC
              DO B = 1, A
                BR = B + NOCC
                AB = (B-1)*NVIR + A
                ERI = BUFF(AB)*TWO
                ANB = A.NE.B
                DO K = NCP1, NOCC
                  KAR = IROT(K,AR)
                  IF (KAR.NE.0) THEN
                    KD = K - NCOR
                    IKD = IA(MAX(ID,KD)) + MIN(ID,KD)
                    MX = K
                    IF (ANB) MX = NOCC
                    DO L = NCP1, MX
                      LBR = IROT(L,BR)
                      IF (LBR.NE.0) THEN
                        IF (KAR.NE.LBR) THEN
                          LD = L - NCOR
                          JLD = IA(MAX(JD,LD)) + MIN(JD,LD)
                          IKJLD = IA(MAX(IKD,JLD)) + MIN(IKD,JLD)
                          DIKJL = TPDM(IKJLD)
                          FAC = DIKJL*ERI
                          IX = KAR + 1
                          JX = LBR + 1
                          PROD(IX) = PROD(IX) + FAC*TVEC(JX)
                          PROD(JX) = PROD(JX) + FAC*TVEC(IX)
                        END IF
                      END IF
                    END DO
                  END IF
                END DO
              END DO
            END DO
            IF (I.NE.J) THEN
              DO A = 1, NVIR
                AR = A + NOCC
                DO B = 1, A
                  BR = B + NOCC
                  BA = (A-1)*NVIR + B
                  ERI = BUFF(BA)*TWO
                  ANB = A.NE.B
                  DO K = NCP1, NOCC
                    KAR = IROT(K,AR)
                    IF (KAR.NE.0) THEN
                      KD = K - NCOR
                      JKD = IA(MAX(JD,KD)) + MIN(JD,KD)
                      MX = K
                      IF (ANB) MX = NOCC
                      DO L = NCP1, MX
                        LBR = IROT(L,BR)
                        IF (LBR.NE.0) THEN
                          IF (KAR.NE.LBR) THEN
                            LD = L - NCOR
                            ILD = IA(MAX(ID,LD)) + MIN(ID,LD)
                            ILJKD = IA(MAX(ILD,JKD)) + MIN(ILD,JKD)
                            DILJK = TPDM(ILJKD)
                            FAC = DILJK*ERI
                            IX = KAR + 1
                            JX = LBR + 1
                            PROD(IX) = PROD(IX) + FAC*TVEC(JX)
                            PROD(JX) = PROD(JX) + FAC*TVEC(IX)
                          END IF
                        END IF
                      END DO
                    END IF
                  END DO
                END DO
              END DO
            END IF
C
C  END DISTRIBUTION LOOPS
C
            CALL DDI_DLBNEXT(MYTASK)
          END IF  ! DLB
        END DO  ! J
      END DO  ! I
C
      IF (BENCH) THEN
        WRITE(ICHANL,9000) '(VO|VO)'
        CALL TIMIT(1)
      END IF
C
C  GLOBALLY SUM THE PRODUCT VECTOR
C
      CALL DDI_GSUMF(2294,PROD,NROT+1)
C
C  DOUBLE PRODUCT VECTOR
C
      CALL DSCAL(NROT+1,TWO,PROD,1)
C
C  CONTRIBUTION FROM DIAGONAL ELEMENTS OF HESSIAN
C
      DO I = 1, NROT+1
        PROD(I) = PROD(I) + DIAH(I)*TVEC(I)
      END DO
      CALL DDI_DLBRESET()
      RETURN
9000  FORMAT(' ..... DONE WITH ',A8,' CONTRIBUTIONS .....')
      END