C  8 MAY 06 - JMR - NTNDVD: PRINT ITERATION ERROR ONLY ON MASTER
C 22 FEB 06 - GDF - AHPDDI: FIX ACTIVE/CORE,ACTIVE/ACTIVE TERM 21
C 17 JAN 06 - MWS - USE NEW DDI TRANSF COMMON TRFDMS
C  5 JUL 05 - MWS - SELECT NEW ATOM,BASIS,EFP,PCM,DAF DIMENSIONS
C  9 JUN 04 - MWS - AHPDDI: RESET DLB COUNTER BEFORE LEAVING
C 16 JAN 04 - GDF - PARALLEL RUNS FORM HESSIAN*TRIAL VEC ON THE FLY
C 12 DEC 02 - MWS - PNTNHES: ENSURE ENTIRE INDEX ARRAY IS FILLED
C 22 MAY 02 - GDF - FORM LAGRANGIAN AND HESSIAN FROM DDI INTEGRALS
C 12 NOV 98 - GDF - CHANGE BIT PACKING TO ISHIFT, DROP ALL M3 CODE
C  8 JAN 97 - GMC - MC1COU,M1EGH,M1EXC,M1LGR: CHANGES TO DROP MCC ORBS
C 13 JUN 96 - MWS - REMOVE FTNCHEK WARNINGS
C 10 NOV 94 - MWS - REMOVE FTNCHEK WARNINGS
C 10 AUG 94 - MWS - INCREASE NUMBER OF DAF RECORDS
C 13 DEC 93 - TLW - M1EGHR: DON'T ASSUME TRANS. INTS ARE IN ANY ORDER
C 11 JAN 92 - TLW - MAKE READ PARALLEL
C 10 JAN 92 - TLW - CHANGE REWINDS TO CALL SEQREW;CALL ABRT BEFORE STOP
C  7 JAN 92 - TLW - MAKE WRITES PARALLEL; ADD COMMON PAR
C  7 OCT 91 - MWS - CHANGE UNIX BYTE UNPACKING
C  3 JUN 90 - JAM - CHANGE VAX BYTE PACKING TO FULL 8 BITS
C 21 DEC 89 - STE - ALPHABETIZE THE ROUTINES
C  3 JAN 89 - MWS - FIX VAX,CEL,CRY,FPS CODE IN M3SRT (IEH, NOT NEH)
C 12 OCT 88 - MWS - IMPLEMENT MICHEL'S HONDO7 MCSCF PROGRAM IN GAMESS,
C                   THIS FILE FORMS THE 2E- LAGRANGIAN AND HESSIAN
C
C*MODULE MCTWO   *DECK M1COU
      SUBROUTINE M1COU(NEH,EH,X,NFTG,G,GIJ,IX,
     *                 NIA,IA,NORB,IC,NINTMX,MI,MA,MAB,NCORBS,DROPC)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      LOGICAL IANDJ,KANDL,SAME,DROPC
C
      DIMENSION EH(NEH)
      DIMENSION IA(NIA),IC(NORB,*)
      DIMENSION X(MAB,*),G(NINTMX),IX(NINTMX),GIJ(*)
C
      COMMON /PCKLAB/ LABSIZ
C
      PARAMETER (GTOL=1.0D-08)
      PARAMETER (ZERO=0.0D+00)
C
      DO 10 I=1,NIA
   10 IA(I)=(I*(I-1))/2
C
      IF (NCORBS.EQ.0 .OR. .NOT. DROPC) GO TO 80
C
C     CONTRIBUTION TO H(RI,SJ) FROM I J CORE
C
C     EH(RI,SJ) = EH(RI,SJ)  - 2 <RS|IJ>
C
      DO 35 M = 1,MA
C
       DO 30 I=1,NCORBS
       MIROT=IC(M+MI,I)
       IF(MIROT.EQ.0) GO TO 30
C
        DO 20 N = 1,MA
        MN=IA(MAX0(M,N))+MIN0(M,N)
C
         DO 15 J=1,NCORBS
         NJROT=IC(N+MI,J)
         IF(NJROT.EQ.0) GO TO 15
C
         IJ=IA(MAX0(I,J)) + MIN0(I,J)
         XC = X(MN,IJ)
         IF(XC.EQ.ZERO) GO TO 15
         XC = XC + XC
C
         NRX = IA(MAX0(MIROT,NJROT))+MIN0(MIROT,NJROT)
         EH(NRX) = EH(NRX) - XC
C
  15     CONTINUE
  20    CONTINUE
  30   CONTINUE
  35  CONTINUE
C
C
C     CONTRIBUTION TO H(RI,SJ) FROM I CORE J VALENCE
C
C     EH(RI,SJ) = EH(RI,SJ) - SUM GJL * ( <RS|IL> )
C                             (OVER L - VALENCE)
C
      DO 60 M = 1,MA
C
       DO 55 I = 1,NCORBS
       MIROT=IC(M+MI,I)
       IF(MIROT.EQ.0) GO TO 55
C
       DO 50 N = 1,MA
       MN=IA(MAX0(M,N))+MIN0(M,N)
C
        DO 45 L = NCORBS+1,MI
        IL=IA(L)+I
        XC = X(MN,IL)
        IF(XC.EQ.ZERO) GO TO 45
C
          DO 40 J = NCORBS+1,MI
          NJROT=IC(N+MI,J)
          IF(NJROT.EQ.0) GO TO 40
          J1=MAX0(J,L)-NCORBS
          L1=MIN0(J,L)-NCORBS
          NGJL = IA(J1)+L1
          XGJL = GIJ(NGJL)
          IF(XGJL.EQ.ZERO) GO TO 40
         FF = XGJL * (XC + XC)
         NRX = IA(MAX0(MIROT,NJROT))+MIN0(MIROT,NJROT)
         EH(NRX) = EH(NRX) - FF
C
  40   CONTINUE
  45   CONTINUE
  50  CONTINUE
  55  CONTINUE
  60  CONTINUE
C
  80  CONTINUE
C
C     ----- READ IN -DM2- FROM -NFTG- -----
C
      CALL SEQREW(NFTG)
  100 CALL PREAD(NFTG,G,IX,NG,NINTMX)
      IF(NG.EQ.0) GO TO 4100
      MGIJKL=IABS(NG)
      IF(MGIJKL.GT.NINTMX) THEN
         CALL ABRT
         STOP
      END IF
      DO 4000 MG=1,MGIJKL
      GGIJKL=G(MG)
      IF( ABS(GGIJKL).LT.GTOL) GO TO 4000
C
                       NPACK = MG
                       IF (LABSIZ .EQ. 2) THEN
                         LABEL = IX(NPACK)
                         IPACK = ISHFT( LABEL, -48 )
                         JPACK = IAND( ISHFT( LABEL, -32 ), 65535 )
                         KPACK = IAND( ISHFT( LABEL, -16 ), 65535 )
                         LPACK = IAND( LABEL, 65535 )
                       ELSE IF (LABSIZ .EQ. 1) THEN
                         IF ( MOD(NPACK,2) .EQ. 0 ) THEN
                           LABEL = IX( NPACK/2 )
                           IPACK = IAND( ISHFT( LABEL, -24 ), 255 )
                           JPACK = IAND( ISHFT( LABEL, -16 ), 255 )
                           KPACK = IAND( ISHFT( LABEL,  -8 ), 255 )
                           LPACK = IAND( LABEL, 255 )
                         ELSE
                           LABEL = IX( (NPACK/2)+1 )
                           IPACK = ISHFT( LABEL, -56 )
                           JPACK = IAND( ISHFT( LABEL, -48 ), 255 )
                           KPACK = IAND( ISHFT( LABEL, -40 ), 255 )
                           LPACK = IAND( ISHFT( LABEL, -32 ), 255 )
                         END IF
                       END IF
                       IGIJKL = IPACK
                       JGIJKL = JPACK
                       KGIJKL = KPACK
                       LGIJKL = LPACK
C
C     -----  I AND J VALENCE -------
C
      IANDJ=IGIJKL.EQ.JGIJKL
      KANDL=KGIJKL.EQ.LGIJKL
      SAME =IGIJKL.EQ.KGIJKL.AND.JGIJKL.EQ.LGIJKL
C
      IF (DROPC) THEN
         IGIJKL = IGIJKL + NCORBS
         JGIJKL = JGIJKL + NCORBS
         KGIJKL = KGIJKL + NCORBS
         LGIJKL = LGIJKL + NCORBS
      END IF
C
C     ----- COULOMB CONTRIBUTION TO HESSIAN -----
C
      GAB=GGIJKL
      IF(.NOT.IANDJ) GAB=GAB+GAB
      NAB=IA(IGIJKL)+JGIJKL
      GCD=GGIJKL
      IF(.NOT.KANDL) GCD=GCD+GCD
      NCD=IA(KGIJKL)+LGIJKL
C
      IG=IGIJKL
      JG=JGIJKL
      KL=NCD
      GMINJ=GCD
      NPASS=1
 2000 CONTINUE
C
      DO 2100 M=1,MA
      GM=GMINJ
      MIROT=IC(M+MI,IG)
      IF(MIROT.EQ.0) GO TO 2100
      MAXN=M
      IF(JG.GT.IG) MAXN=M-1
      IF(MAXN.EQ.0) GO TO 2100
      DO 2050 N=1,MAXN
      GIJKL=GM
      NJROT=IC(N+MI,JG)
      IF(NJROT.EQ.0) GO TO 2050
      MN=IA(M)+N
      XMNKL=X(MN,KL)
      IF(XMNKL.EQ.ZERO) GO TO 2050
      IF(N.NE.M.OR.JG.NE.IG) GIJKL=GIJKL+GIJKL
      MIH=MAX(MIROT,NJROT)
      NJH=MIN(MIROT,NJROT)
      MINJH=IA(MIH)+NJH
      EH(MINJH)=EH(MINJH)+ GIJKL*XMNKL
 2050 CONTINUE
 2100 CONTINUE
C
 2200 NPASS=NPASS+1
      GO TO (2300,2210,2220,2230,2300),NPASS
 2210 IF(IANDJ) GO TO 2200
      IG=JGIJKL
      JG=IGIJKL
      GO TO 2000
 2220 IF(SAME) GO TO 2300
      IG=KGIJKL
      JG=LGIJKL
      KL=NAB
      GMINJ=GAB
      GO TO 2000
 2230 IF(KANDL) GO TO 2300
      IG=LGIJKL
      JG=KGIJKL
      GO TO 2000
 2300 CONTINUE
C
 4000 CONTINUE
      IF(NG.GT.0) GO TO 100
 4100 CONTINUE
      RETURN
      END
C*MODULE MCTWO   *DECK M1COUR
      SUBROUTINE M1COUR(NFTI,OX,X,IX,IA,NINTMX,MI,MA,MIJ,MAB,NOSQUR)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      LOGICAL SKIP
C
      DIMENSION OX(MAB,*),X(NINTMX),IX(NINTMX),IA(*)
C
      COMMON /PCKLAB/ LABSIZ
C
C     ----- READ (VIR,VIR//OCC,OCC) INTEGRALS FOR M1COU -----
C
      NIA=MAX(MI,MA)
      DO 10 I=1,NIA
   10 IA(I)=(I*(I-1))/2
      CALL VCLR(OX,1,MAB*MIJ)
C
C     ----- READ IN INTEGRALS FROM -NFTI- -----
C
  110 CALL PREAD(NFTI,X,IX,NX,NINTMX)
      IF(NX.EQ.0) GO TO 200
      MX=IABS(NX)
      IF(MX.GT.NINTMX) THEN
         CALL ABRT
         STOP
      END IF
      DO 120 M=1,MX
C
                       NPACK = M
                       IF (LABSIZ .EQ. 2) THEN
                         LABEL = IX(NPACK)
                         IPACK = ISHFT( LABEL, -48 )
                         JPACK = IAND( ISHFT( LABEL, -32 ), 65535 )
                         KPACK = IAND( ISHFT( LABEL, -16 ), 65535 )
                         LPACK = IAND( LABEL, 65535 )
                       ELSE IF (LABSIZ .EQ. 1) THEN
                         IF ( MOD(NPACK,2) .EQ. 0 ) THEN
                           LABEL = IX( NPACK/2 )
                           IPACK = IAND( ISHFT( LABEL, -24 ), 255 )
                           JPACK = IAND( ISHFT( LABEL, -16 ), 255 )
                           KPACK = IAND( ISHFT( LABEL,  -8 ), 255 )
                           LPACK = IAND( LABEL, 255 )
                         ELSE
                           LABEL = IX( (NPACK/2)+1 )
                           IPACK = ISHFT( LABEL, -56 )
                           JPACK = IAND( ISHFT( LABEL, -48 ), 255 )
                           KPACK = IAND( ISHFT( LABEL, -40 ), 255 )
                           LPACK = IAND( ISHFT( LABEL, -32 ), 255 )
                         END IF
                       END IF
                       K = IPACK
                       L = JPACK
                       I = KPACK
                       J = LPACK
C
         IF(NOSQUR.EQ.1) THEN
            IDUM=I
            JDUM=J
            I=K
            J=L
            K=IDUM
            L=JDUM
         END IF
         SKIP=(I.LE.MI).OR.(J.LE.MI).OR.(K.GT.MI).OR.(L.GT.MI)
         IF(SKIP) GO TO 120
         IJ=IA(I-MI)+(J-MI)
         KL=IA(K)+L
         OX(IJ,KL)=X(M)
  120 CONTINUE
      IF(NX.GT.0) GO TO 110
  200 CONTINUE
      RETURN
      END
C*MODULE MCTWO   *DECK M1EGH
      SUBROUTINE M1EGH(NEG,NEH,EG,EH,X,NFTG,G,GIJ,IX,
     *                 NIA,IA,NORB,IB,IC,NINTMX,MI,NCORBS,DROPC)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      LOGICAL IANDJ,KANDL,SAME,DROPC
C
      DIMENSION EG(NEG),EH(NEH)
      DIMENSION IA(NIA),IB(NORB),IC(NORB,*)
      DIMENSION X(*),G(NINTMX),IX(NINTMX),GIJ(*)
C
      COMMON /PCKLAB/ LABSIZ
C
      PARAMETER (GTOL=1.0D-08)
      PARAMETER (ZERO=0.0D+00)
C
      DO 10 I=1,NIA
   10 IA(I)=(I*(I-1))/2
C
      IF (.NOT. DROPC) THEN
         CALL VCLR(EG,1,NEG)
         CALL VCLR(EH,1,NEH)
      END IF
C
      IF (NCORBS.EQ.0 .OR. .NOT. DROPC) GO TO 80
C
C     CONTRIBUTION TO H(RI,SJ) FROM I J CORE
C
C     EH(RI,SJ) = EH(RI,SJ) + 8 <RI|SJ> - 2 <RJ|SI> - 2 <RS|IJ>
C
      DO 35 M = NCORBS+1,MI
      M1 = IA(M)
C
       DO 30 I=1,NCORBS
       MIROT=IC(M,I)
       IF(MIROT.EQ.0) GO TO 30
       MII= M1 + I
C
        DO 25 N = NCORBS+1,M
        N1 = IA(N)
        MN=IA(M) + N
        NII= N1 + I
C
          DO 20 J=1,NCORBS
          NJROT=IC(N,J)
          IF(NJROT.EQ.0) GO TO 20
          IJ=IA(MAX0(I,J)) + MIN0(I,J)
C
          NXC = IA(MN) + IJ
          XNXC = X(NXC)
          NJJ= N1 + J
          NEC1=IA(MAX0(MII,NJJ))+ MIN0(NJJ,MII)
          XNEC1=X(NEC1)
          MJJ= M1 + J
          NEC2=IA(MAX0(MJJ,NII)) + MIN0(MJJ,NII)
          XNEC2=X(NEC2)
C
          XC  = (XNEC1 + XNEC1)+(XNEC1 + XNEC1) - XNEC2 -XNXC
          IF (XC.EQ.ZERO) GO TO 20
          XC = XC + XC
          IF(M.NE.N) XC = XC+XC
          NRX = IA(MAX0(MIROT,NJROT))+MIN0(MIROT,NJROT)
          EH(NRX) = EH(NRX) + XC
C
 20      CONTINUE
 25     CONTINUE
 30    CONTINUE
 35   CONTINUE
C
C     CONTRIBUTION TO H(RI,SJ) FROM I CORE J VALENCE
C
C     EH(RI,SJ) = EH(RI,SJ) + SUM GJL * ( 4 <RI|SL> - <RL|SI> - <RS|IL>)
C                             (OVER L - VALENCE)
C
      DO 65 M = 1,MI
      MINM=1
      IF(M.LE.NCORBS) MINM=NCORBS+1
C
       DO 60 I=1,NCORBS
       MIROT=IC(M,I)
       IF(MIROT.EQ.0) GO TO 60
C
        DO 50 N = MINM,MI
        MN=IA(MAX0(M,N))+MIN0(M,N)
C
        DO 40 L = NCORBS+1,MI
        IL=IA(L)+I
        NXC = IA(MAX0(MN,IL)) + MIN0(MN,IL)
        XNXC=X(NXC)
        MMM = IA(MAX0(M,I))+ MIN0(M,I)
        NNN = IA(MAX0(N,L)) + MIN0(N,L)
        NEC1=IA(MAX0(MMM,NNN)) + MIN0(MMM,NNN)
        XNEC1=X(NEC1)
        MMM = IA(MAX0(M,L))+ MIN0(M,L)
        NNN = IA(MAX0(N,I))+ MIN0(N,I)
        NEC2=IA(MAX0(MMM,NNN)) + MIN0(MMM,NNN)
        XNEC2=X(NEC2)
         XC  = (XNEC1+XNEC1)+(XNEC1+XNEC1) -XNEC2 - XNXC
C
         DO 40 J = NCORBS+1,MI
          NJROT=IC(N,J)
          IF(NJROT.EQ.0) GO TO 40
          J1=MAX0(J,L)-NCORBS
          L1=MIN0(J,L)-NCORBS
          NGJL = IA(J1) + L1
          XGJL = GIJ(NGJL)
         IF(XGJL.EQ.ZERO) GO TO 40
C
          FF =  XGJL * (XC + XC )
          IF (M.LT.I) FF = -FF
          IF (N.LT.J) FF = -FF
         NRX = IA(MAX0(MIROT,NJROT)) + MIN0(MIROT,NJROT)
         EH(NRX) = EH(NRX) + FF
  40    CONTINUE
  50   CONTINUE
  60  CONTINUE
  65  CONTINUE
C
  80  CONTINUE
C
C     ----- READ IN -DM2- FROM -NFTG- -----
C
      CALL SEQREW(NFTG)
  100 CALL PREAD(NFTG,G,IX,NG,NINTMX)
      IF(NG.EQ.0) GO TO 4100
      MGIJKL=IABS(NG)
      IF(MGIJKL.GT.NINTMX) THEN
         CALL ABRT
         STOP
      END IF
      DO 4000 MG=1,MGIJKL
      GGIJKL=G(MG)
      IF( ABS(GGIJKL).LT.GTOL) GO TO 4000
C
                       NPACK = MG
                       IF (LABSIZ .EQ. 2) THEN
                         LABEL = IX(NPACK)
                         IPACK = ISHFT( LABEL, -48 )
                         JPACK = IAND( ISHFT( LABEL, -32 ), 65535 )
                         KPACK = IAND( ISHFT( LABEL, -16 ), 65535 )
                         LPACK = IAND( LABEL, 65535 )
                       ELSE IF (LABSIZ .EQ. 1) THEN
                         IF ( MOD(NPACK,2) .EQ. 0 ) THEN
                           LABEL = IX( NPACK/2 )
                           IPACK = IAND( ISHFT( LABEL, -24 ), 255 )
                           JPACK = IAND( ISHFT( LABEL, -16 ), 255 )
                           KPACK = IAND( ISHFT( LABEL,  -8 ), 255 )
                           LPACK = IAND( LABEL, 255 )
                         ELSE
                           LABEL = IX( (NPACK/2)+1 )
                           IPACK = ISHFT( LABEL, -56 )
                           JPACK = IAND( ISHFT( LABEL, -48 ), 255 )
                           KPACK = IAND( ISHFT( LABEL, -40 ), 255 )
                           LPACK = IAND( ISHFT( LABEL, -32 ), 255 )
                         END IF
                       END IF
                       IGIJKL = IPACK
                       JGIJKL = JPACK
                       KGIJKL = KPACK
                       LGIJKL = LPACK
C
C
C     ----- I AND J VALENCE -------
C
      IANDJ=IGIJKL.EQ.JGIJKL
      KANDL=KGIJKL.EQ.LGIJKL
      SAME =IGIJKL.EQ.KGIJKL.AND.JGIJKL.EQ.LGIJKL
C
      IF (DROPC) THEN
         IGIJKL = IGIJKL + NCORBS
         JGIJKL = JGIJKL + NCORBS
         KGIJKL = KGIJKL + NCORBS
         LGIJKL = LGIJKL + NCORBS
      END IF
C
C     ----- LAGRANGIAN CONTRIBUTION -----
C
      GAB=GGIJKL
      IF(.NOT.IANDJ) GAB=GAB+GAB
      NAB=IA(IGIJKL)+JGIJKL
      GCD=GGIJKL
      IF(.NOT.KANDL) GCD=GCD+GCD
      NCD=IA(KGIJKL)+LGIJKL
C
      IG=IGIJKL
      JG=JGIJKL
      KG=KGIJKL
      LG=LGIJKL
      KL=NCD
      GIJKL=GCD
      NPASS=1
 1000 CONTINUE
C
C     ----- COMPUTE CONTRIBUTION TO THE LAGRANGIAN -----
C
      DO 1100 M=1,MI
         MM=MAX(M,JG)
         JJ=MIN(M,JG)
         MJ=IA(MM)+JJ
         MJX=MAX(MJ,KL)
         KLX=MIN(MJ,KL)
         MJKLX=IA(MJX)+KLX
         XMJKL=X(MJKLX)
         IF(XMJKL.EQ.ZERO) GO TO 1100
         MILGR=IB(M)+IG
         EG(MILGR)=EG(MILGR)+ GIJKL*XMJKL
 1100 CONTINUE
C
C     ----- CONTRIBUTIONS TO ORBITAL HESSIAN -----
C
 1200 NPASS=NPASS+1
      GO TO (1300,1210,1220,1230,1300),NPASS
 1210 IF(IANDJ) GO TO 1200
      IG=JGIJKL
      JG=IGIJKL
      GO TO 1000
 1220 IF(SAME) GO TO 1300
      IG=KGIJKL
      JG=LGIJKL
      KG=IGIJKL
      LG=JGIJKL
      KL=NAB
      GIJKL=GAB
      GO TO 1000
 1230 IF(KANDL) GO TO 1300
      IG=LGIJKL
      JG=KGIJKL
      GO TO 1000
 1300 CONTINUE
C
C     ----- COULOMB CONTRIBUTION TO HESSIAN -----
C
      IG=IGIJKL
      JG=JGIJKL
      KG=KGIJKL
      LG=LGIJKL
      KL=NCD
      GMINJ=GCD
      NPASS=1
 2000 CONTINUE
C
      DO 2100 M=1,MI
      GM=GMINJ
      IF(M -IG) 2010,2100,2020
 2010 GM=-GM
 2020 MIROT=IC(M,IG)
      IF(MIROT.EQ.0) GO TO 2100
      MAXN=M
      IF(JG.GT.IG) MAXN=M-1
      IF(MAXN.EQ.0) GO TO 2100
      DO 2050 N=1,MAXN
      GIJKL=GM
      IF(N -JG) 2030,2050,2040
 2030 GIJKL=-GIJKL
 2040 NJROT=IC(N,JG)
      IF(NJROT.EQ.0) GO TO 2050
      MN=IA(M)+N
      MNX=MAX(MN,KL)
      KLX=MIN(MN,KL)
      MNKLX=IA(MNX)+KLX
      XMNKL=X(MNKLX)
      IF(XMNKL.EQ.ZERO) GO TO 2050
      IF(N.NE.M.OR.JG.NE.IG) GIJKL=GIJKL+GIJKL
      MIH=MAX(MIROT,NJROT)
      NJH=MIN(MIROT,NJROT)
      MINJH=IA(MIH)+NJH
      EH(MINJH)=EH(MINJH)+ GIJKL*XMNKL
 2050 CONTINUE
 2100 CONTINUE
C
 2200 NPASS=NPASS+1
      GO TO (2300,2210,2220,2230,2300),NPASS
 2210 IF(IANDJ) GO TO 2200
      IG=JGIJKL
      JG=IGIJKL
      GO TO 2000
 2220 IF(SAME) GO TO 2300
      IG=KGIJKL
      JG=LGIJKL
      KG=IGIJKL
      LG=JGIJKL
      KL=NAB
      GMINJ=GAB
      GO TO 2000
 2230 IF(KANDL) GO TO 2300
      IG=LGIJKL
      JG=KGIJKL
      GO TO 2000
 2300 CONTINUE
C
C     ----- EXCHANGE CONTRIBUTION TO HESSIAN -----
C
      GDM2=GGIJKL+GGIJKL
      IG=IGIJKL
      JG=JGIJKL
      KG=KGIJKL
      LG=LGIJKL
      NPASS=1
 3000 CONTINUE
C
      DO 3100 M=1,MI
      GM=GDM2
      IF(M -IG) 3010,3100,3020
 3010 GM=-GM
 3020 MIROT=IC(M,IG)
      IF(MIROT.EQ.0) GO TO 3100
      MM=MAX(M,JG)
      JJ=MIN(M,JG)
      MJ=IA(MM)+JJ
      MAXN=M
      IF(KG.GT.IG) MAXN=M-1
      IF(MAXN.EQ.0) GO TO 3100
      DO 3050 N=1,MAXN
      GIJKL=GM
      IF(N -KG) 3030,3050,3040
 3030 GIJKL=-GIJKL
 3040 NKROT=IC(N,KG)
      IF(NKROT.EQ.0) GO TO 3050
      NN=MAX(N,LG)
      LL=MIN(N,LG)
      NL=IA(NN)+LL
      MJX=MAX(MJ,NL)
      NLX=MIN(MJ,NL)
      MJNLX=IA(MJX)+NLX
      XMJNL=X(MJNLX)
      IF(XMJNL.EQ.ZERO) GO TO 3050
      IF(N.NE.M.OR.KG.NE.IG) GIJKL=GIJKL+GIJKL
      MIH=MAX(MIROT,NKROT)
      NKH=MIN(MIROT,NKROT)
      MINKH=IA(MIH)+NKH
      EH(MINKH)=EH(MINKH)+ GIJKL*XMJNL
 3050 CONTINUE
 3100 CONTINUE
C
 3200 NPASS=NPASS+1
      GO TO (3300,3210,3220,3230,3240,3250,3260,3270,3300),NPASS
 3210 IF(IANDJ) GO TO 3200
      IG=JGIJKL
      JG=IGIJKL
      GO TO 3000
 3220 IF(KANDL) GO TO 3200
      IG=IGIJKL
      JG=JGIJKL
      KG=LGIJKL
      LG=KGIJKL
      GO TO 3000
 3230 IF(IANDJ.OR.KANDL) GO TO 3200
      IG=JGIJKL
      JG=IGIJKL
      GO TO 3000
 3240 IF(SAME) GO TO 3300
      IG=KGIJKL
      JG=LGIJKL
      KG=IGIJKL
      LG=JGIJKL
      GO TO 3000
 3250 IF(IANDJ) GO TO 3200
      KG=JGIJKL
      LG=IGIJKL
      GO TO 3000
 3260 IF(KANDL) GO TO 3200
      IG=LGIJKL
      JG=KGIJKL
      KG=IGIJKL
      LG=JGIJKL
      GO TO 3000
 3270 IF(IANDJ.OR.KANDL) GO TO 3300
      KG=JGIJKL
      LG=IGIJKL
      GO TO 3000
 3300 CONTINUE
C
 4000 CONTINUE
      IF(NG.GT.0) GO TO 100
 4100 CONTINUE
      RETURN
      END
C*MODULE MCTWO   *DECK M1EGHR
      SUBROUTINE M1EGHR(NFTI,OX,X,IX,IA,NINTMX,MI,NOSQUR)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      DIMENSION OX(*),X(NINTMX),IX(NINTMX),IA(*)
C
      COMMON /PCKLAB/ LABSIZ
C
C     ----- READ (OCC,OCC//OCC,OCC) INTEGRALS FOR M1EGH -----
C     INTEGRALS FROM THE FULL TRANSFORMATION ARE EXPECTED TO
C     BE IN REVERSE CANONICAL ORDER, INTEGRALS FROM THE PARTIAL
C     TRANSFORMATION IN NORMAL CANONICAL ORDER.  THE FORM OF
C     THE SORTED AO INTEGRAL FILE DISTINQUISHES THESE TWO CASES.
C
      M2=(MI*(MI+1))/2
      M4=(M2*(M2+1))/2
      DO 10 I=1,M2
   10 IA(I)=(I*(I-1))/2
      CALL VCLR(OX,1,M4)
C
C     ----- READ IN INTEGRALS FROM -NFTI- -----
C
  110 CALL PREAD(NFTI,X,IX,NX,NINTMX)
      IF(NX.EQ.0) GO TO 200
      MX=IABS(NX)
      IF(MX.GT.NINTMX) CALL ABRT
      DO 120 M=1,MX
C
                       NPACK = M
                       IF (LABSIZ .EQ. 2) THEN
                         LABEL = IX(NPACK)
                         IPACK = ISHFT( LABEL, -48 )
                         JPACK = IAND( ISHFT( LABEL, -32 ), 65535 )
                         KPACK = IAND( ISHFT( LABEL, -16 ), 65535 )
                         LPACK = IAND( LABEL, 65535 )
                       ELSE IF (LABSIZ .EQ. 1) THEN
                         IF ( MOD(NPACK,2) .EQ. 0 ) THEN
                           LABEL = IX( NPACK/2 )
                           IPACK = IAND( ISHFT( LABEL, -24 ), 255 )
                           JPACK = IAND( ISHFT( LABEL, -16 ), 255 )
                           KPACK = IAND( ISHFT( LABEL,  -8 ), 255 )
                           LPACK = IAND( LABEL, 255 )
                         ELSE
                           LABEL = IX( (NPACK/2)+1 )
                           IPACK = ISHFT( LABEL, -56 )
                           JPACK = IAND( ISHFT( LABEL, -48 ), 255 )
                           KPACK = IAND( ISHFT( LABEL, -40 ), 255 )
                           LPACK = IAND( ISHFT( LABEL, -32 ), 255 )
                         END IF
                       END IF
                       K = IPACK
                       L = JPACK
                       I = KPACK
                       J = LPACK
C
         IF(NOSQUR.EQ.1) THEN
            IDUM=I
            JDUM=J
            I=K
            J=L
            K=IDUM
            L=JDUM
            IF(I.GT.MI) GO TO 200
         ELSE
            IF(K.GT.MI) GO TO 120
            IF(I.GT.MI) GO TO 120
         END IF
         IJ=IA(I)+J
         KL=IA(K)+L
         IJKL=IA(IJ)+KL
         OX(IJKL)=X(M)
  120 CONTINUE
      IF(NX.GT.0) GO TO 110
C
  200 CONTINUE
      RETURN
      END
C*MODULE MCTWO   *DECK M1EXC
      SUBROUTINE M1EXC(NEH,EH,X,NFTG,G,GIJ,IX,
     *                 NIA,IA,NORB,IC,NINTMX,MI,MA,NCORBS,DROPC)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      LOGICAL IANDJ,KANDL,SAME,DROPC
C
      DIMENSION EH(NEH),IA(NIA),IC(NORB,*),X(*)
      DIMENSION G(NINTMX),IX(NINTMX),GIJ(*)
C
      COMMON /PCKLAB/ LABSIZ
C
      PARAMETER (GTOL=1.0D-08)
      PARAMETER (ZERO=0.0D+00)
C
      DO 10 I=1,NIA
   10 IA(I)=(I*I-I)/2
C
      IF (NCORBS.EQ.0 .OR. .NOT. DROPC) GO TO 80
C
C     CONTRIBUTION TO H(RI,SJ) FROM I J CORE
C
C     EH(RI,SJ) = EH(RI,SJ) + 8 <RI|SJ> - 2 <RJ|SI>
C
      DO 35 M = 1,MA
      M1 = (M-1)*MI
C
      DO 30 I=1,NCORBS
       MIROT=IC(M+MI,I)
       IF(MIROT.EQ.0) GO TO 30
       MMI = M1 + I
C
       DO 20 N = 1,MA
       N1 = (N-1)*MI
       NNI = N1 + I
C
       DO 15 J=1,NCORBS
       NJROT=IC(N+MI,J)
       IF(NJROT.EQ.0) GO TO 15
C
       NNJ = N1 + J
       NEC1=IA(MAX0(MMI,NNJ)) + MIN0(MMI,NNJ)
       MMJ = M1 + J
       NEC2=IA(MAX0(MMJ,NNI))+MIN0(MMJ,NNI)
C
C      ----- EXCH -----
C
       XC=(X(NEC1) + X(NEC1))+(X(NEC1)+X(NEC1))-X(NEC2)
       IF(XC.EQ.ZERO) GO TO 15
       XC = XC + XC
       NRX = IA(MAX0(MIROT,NJROT))+MIN0(MIROT,NJROT)
       EH(NRX) = EH(NRX) + XC
C
 15   CONTINUE
 20   CONTINUE
 30   CONTINUE
 35   CONTINUE
C
C     ----- CONTRIBUTION TO H(RI,SJ) FROM I CORE J VALENCE -----
C
      DO 65 M = 1,MA
      M1 = (M-1)*MI
C
       DO 60 I=1,NCORBS
       MIROT=IC(M+MI,I)
       IF(MIROT.EQ.0) GO TO 60
       MMI = M1 + I
C
        DO 50 N = 1,MA
C       MN=IA(MAX0(M,N))+MIN0(M,N)
        N1 = (N-1)*MI
        NNI = N1 + I
C
         DO 45 L = NCORBS+1,MI
         NNL = N1 + L
         NEC1=IA(MAX0(MMI,NNL))+MIN0(MMI,NNL)
         MML = M1 + L
         NEC2=IA(MAX0(MML,NNI))+MIN0(MML,NNI)
         XC1 = (X(NEC1)+X(NEC1))+(X(NEC1)+X(NEC1)) - X(NEC2)
         IF(XC1.EQ.ZERO) GO TO 45
C
          DO 40 J = NCORBS+1,MI
          NJROT=IC(N+MI,J)
          IF(NJROT.EQ.0) GO TO 40
          J1=MAX0(J,L)-NCORBS
          L1=MIN0(J,L)-NCORBS
          NGJL = IA(J1)+L1
          XGJL = GIJ(NGJL)
          IF(XGJL.EQ.ZERO) GO TO 40
C
C    EH(RI,SJ) = EH(RI,SJ) + SUM GJL * ( 4 <RI|SL> - <RL|SI> )
C                             (OVER L - VALENCE)
C
           FF = XGJL * (XC1+XC1)
           NRX = IA(MAX0(MIROT,NJROT))+MIN0(MIROT,NJROT)
           EH(NRX) = EH(NRX) + FF
C
  40      CONTINUE
  45     CONTINUE
  50    CONTINUE
  60   CONTINUE
  65  CONTINUE
C
  80  CONTINUE
C
C     ----- READ IN -DM2- FROM -NFTG- -----
C
      CALL SEQREW(NFTG)
  100 CALL PREAD(NFTG,G,IX,NG,NINTMX)
      IF(NG.EQ.0) GO TO 4100
      MGIJKL=IABS(NG)
      IF(MGIJKL.GT.NINTMX) THEN
         CALL ABRT
         STOP
      END IF
      DO 4000 MG=1,MGIJKL
      GGIJKL=G(MG)
      IF( ABS(GGIJKL).LT.GTOL) GO TO 4000
C
                       NPACK = MG
                       IF (LABSIZ .EQ. 2) THEN
                         LABEL = IX(NPACK)
                         IPACK = ISHFT( LABEL, -48 )
                         JPACK = IAND( ISHFT( LABEL, -32 ), 65535 )
                         KPACK = IAND( ISHFT( LABEL, -16 ), 65535 )
                         LPACK = IAND( LABEL, 65535 )
                       ELSE IF (LABSIZ .EQ. 1) THEN
                         IF ( MOD(NPACK,2) .EQ. 0 ) THEN
                           LABEL = IX( NPACK/2 )
                           IPACK = IAND( ISHFT( LABEL, -24 ), 255 )
                           JPACK = IAND( ISHFT( LABEL, -16 ), 255 )
                           KPACK = IAND( ISHFT( LABEL,  -8 ), 255 )
                           LPACK = IAND( LABEL, 255 )
                         ELSE
                           LABEL = IX( (NPACK/2)+1 )
                           IPACK = ISHFT( LABEL, -56 )
                           JPACK = IAND( ISHFT( LABEL, -48 ), 255 )
                           KPACK = IAND( ISHFT( LABEL, -40 ), 255 )
                           LPACK = IAND( ISHFT( LABEL, -32 ), 255 )
                         END IF
                       END IF
                       IGIJKL = IPACK
                       JGIJKL = JPACK
                       KGIJKL = KPACK
                       LGIJKL = LPACK
C
C     ----- I AND J VALENCE -------
C
      IANDJ=IGIJKL.EQ.JGIJKL
      KANDL=KGIJKL.EQ.LGIJKL
      SAME =IGIJKL.EQ.KGIJKL.AND.JGIJKL.EQ.LGIJKL
C
      IF (DROPC) THEN
         IGIJKL = IGIJKL + NCORBS
         JGIJKL = JGIJKL + NCORBS
         KGIJKL = KGIJKL + NCORBS
         LGIJKL = LGIJKL + NCORBS
      END IF
C
C     ----- EXCHANGE CONTRIBUTION TO HESSIAN -----
C
      GDM2=GGIJKL+GGIJKL
      IG=IGIJKL
      JG=JGIJKL
      KG=KGIJKL
      LG=LGIJKL
      NPASS=1
 3000 CONTINUE
C
      DO 3100 M=1,MA
      GM=GDM2
      MIROT=IC(M+MI,IG)
      IF(MIROT.EQ.0) GO TO 3100
      MJ=(M-1)*MI+JG
      MAXN=M
      IF(KG.GT.IG) MAXN=M-1
      IF(MAXN.EQ.0) GO TO 3100
      DO 3050 N=1,MAXN
      GIJKL=GM
      NKROT=IC(N+MI,KG)
      IF(NKROT.EQ.0) GO TO 3050
      NL=(N-1)*MI+LG
      MJX=MAX(MJ,NL)
      NLX=MIN(MJ,NL)
      MJNLX=IA(MJX)+NLX
      XMJNL=X(MJNLX)
      IF(XMJNL.EQ.ZERO) GO TO 3050
      IF(N.NE.M.OR.KG.NE.IG) GIJKL=GIJKL+GIJKL
      MIH=MAX(MIROT,NKROT)
      NKH=MIN(MIROT,NKROT)
      MINKH=IA(MIH)+NKH
      EH(MINKH)=EH(MINKH)+ GIJKL*XMJNL
 3050 CONTINUE
 3100 CONTINUE
C
 3200 NPASS=NPASS+1
      GO TO (3300,3210,3220,3230,3240,3250,3260,3270,3300),NPASS
 3210 IF(IANDJ) GO TO 3200
      IG=JGIJKL
      JG=IGIJKL
      GO TO 3000
 3220 IF(KANDL) GO TO 3200
      IG=IGIJKL
      JG=JGIJKL
      KG=LGIJKL
      LG=KGIJKL
      GO TO 3000
 3230 IF(IANDJ.OR.KANDL) GO TO 3200
      IG=JGIJKL
      JG=IGIJKL
      GO TO 3000
 3240 IF(SAME) GO TO 3300
      IG=KGIJKL
      JG=LGIJKL
      KG=IGIJKL
      LG=JGIJKL
      GO TO 3000
 3250 IF(IANDJ) GO TO 3200
      KG=JGIJKL
      LG=IGIJKL
      GO TO 3000
 3260 IF(KANDL) GO TO 3200
      IG=LGIJKL
      JG=KGIJKL
      KG=IGIJKL
      LG=JGIJKL
      GO TO 3000
 3270 IF(IANDJ.OR.KANDL) GO TO 3300
      KG=JGIJKL
      LG=IGIJKL
      GO TO 3000
 3300 CONTINUE
C
 4000 CONTINUE
      IF(NG.GT.0) GO TO 100
 4100 CONTINUE
      RETURN
      END
C*MODULE MCTWO   *DECK M1EXCR
      SUBROUTINE M1EXCR(NFTI,OX,X,IX,IA,NINTMX,MI,MAJ,NOSQUR)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      LOGICAL SKIP
C
      DIMENSION OX(*),X(NINTMX),IX(NINTMX),IA(*)
C
      COMMON /PCKLAB/ LABSIZ
C
C     ----- READ (VIR,OCC/VIR,OCC) INTEGRALS FOR M1EXC -----
C
      M2=MAJ
      M4=(M2*(M2+1))/2
      DO 10 I=1,M2
   10 IA(I)=(I*(I-1))/2
      CALL VCLR(OX,1,M4)
C
C     ----- READ IN INTEGRALS FROM -NFTI- -----
C     INTEGRALS ASSUMED IN REVERSE CANONICAL ORDER
C
  110 CALL PREAD(NFTI,X,IX,NX,NINTMX)
      IF(NX.EQ.0) GO TO 200
      MX=IABS(NX)
      IF(MX.GT.NINTMX) THEN
         CALL ABRT
         STOP
      END IF
      DO 120 M=1,MX
C
                       NPACK = M
                       IF (LABSIZ .EQ. 2) THEN
                         LABEL = IX(NPACK)
                         IPACK = ISHFT( LABEL, -48 )
                         JPACK = IAND( ISHFT( LABEL, -32 ), 65535 )
                         KPACK = IAND( ISHFT( LABEL, -16 ), 65535 )
                         LPACK = IAND( LABEL, 65535 )
                       ELSE IF (LABSIZ .EQ. 1) THEN
                         IF ( MOD(NPACK,2) .EQ. 0 ) THEN
                           LABEL = IX( NPACK/2 )
                           IPACK = IAND( ISHFT( LABEL, -24 ), 255 )
                           JPACK = IAND( ISHFT( LABEL, -16 ), 255 )
                           KPACK = IAND( ISHFT( LABEL,  -8 ), 255 )
                           LPACK = IAND( LABEL, 255 )
                         ELSE
                           LABEL = IX( (NPACK/2)+1 )
                           IPACK = ISHFT( LABEL, -56 )
                           JPACK = IAND( ISHFT( LABEL, -48 ), 255 )
                           KPACK = IAND( ISHFT( LABEL, -40 ), 255 )
                           LPACK = IAND( ISHFT( LABEL, -32 ), 255 )
                         END IF
                       END IF
                       K = IPACK
                       L = JPACK
                       I = KPACK
                       J = LPACK
C
         IF(NOSQUR.EQ.1) THEN
            IDUM=I
            JDUM=J
            I=K
            J=L
            K=IDUM
            L=JDUM
         END IF
         SKIP=(I.LE.MI).OR.(J.GT.MI).OR.(K.LE.MI).OR.(L.GT.MI)
         IF(SKIP) GO TO 120
         IJ=(I-MI-1)*MI+J
         KL=(K-MI-1)*MI+L
         IJKL=IA(IJ)+KL
         OX(IJKL)=X(M)
  120 CONTINUE
      IF(NX.GT.0) GO TO 110
  200 CONTINUE
      RETURN
      END
C*MODULE MCTWO   *DECK M1LGR
      SUBROUTINE M1LGR(NEG,NEH,EG,EH,X,NFTG,G,GIJ,IX,
     *                 NIA,IA,NORB,IB,IC,NINTMX,MI,MA,MIJ,NCORBS,DROPC)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      LOGICAL IANDJ,KANDL,SAME,DROPC
C
      DIMENSION EG(NEG),EH(NEH)
      DIMENSION IA(NIA),IB(NORB),IC(NORB,*)
      DIMENSION X(MA,*),G(NINTMX),IX(NINTMX),GIJ(*)
C
      COMMON /PCKLAB/ LABSIZ
C
      PARAMETER (GTOL=1.0D-08)
      PARAMETER (ZERO=0.0D+00)
C
      DO 10 I=1,NIA
   10 IA(I)=(I*(I-1))/2
C
      IF (NCORBS.EQ.0 .OR. .NOT. DROPC) GO TO 80
C
C     ---- COULOMB CONTRIBUTION FROM I AND J CORE, I CORE J VAL,
C                                                  I VAL J CORE
C
      DO 25 M = 1,MA
C
      DO 20 I=1,NCORBS
       MIROT=IC(M+MI,I)
       IF(MIROT.EQ.0) GO TO 20
C
      DO 15 N = NCORBS+1,MI
C
        DO 5 J=1,NCORBS
         NJROT=IC(N,J)
         IF(NJROT.EQ.0) GO TO 5
C
         IJ=IA(MAX0(I,J))+MIN0(I,J)
         NXC =(N-1)*MIJ+IJ
         NN1 = IA(N)+J
         NEC1=(I-1)*MIJ +NN1
         NNN = IA(N)+ I
         NEC2=(J-1)*MIJ+NNN
C
      XC=(X(M,NEC1)+X(M,NEC1))+(X(M,NEC1)+X(M,NEC1))-X(M,NEC2)-X(M,NXC)
         IF(XC.EQ.ZERO) GO TO 5
      FF = (XC+XC) + (XC+XC)
         NRX = IA(MAX0(MIROT,NJROT))+ MIN0(MIROT,NJROT)
         EH(NRX) = EH(NRX) + FF
C
  5   CONTINUE
 15   CONTINUE
 20   CONTINUE
 25   CONTINUE
C
C     ----- I CORE J VALENCE OR I VALENCE AND J CORE -----
C
      DO 65 M = 1,MA
C
       DO 60 I=1,NCORBS
       MIROT=IC(M+MI,I)
C
        DO 50 N = 1,MI
        NIROT=IC(N,I)
        IF(MIROT.EQ.0.AND.NIROT.EQ.0) GO TO 50
C
         DO 45 L = NCORBS+1,MI
         IL=IA(L)+I
C
       NXC =(N-1)*MIJ+IL
       XNXC = X(M,NXC)
       NNN = IA(MAX0(N,L))+MIN0(N,L)
       NEC1=(I-1)*MIJ +NNN
       XNEC1 = X(M,NEC1)
       NNN = IA(MAX0(N,I))+MIN0(N,I)
       NEC2=(L-1)*MIJ+NNN
       XNEC2 = X(M,NEC2)
       XC1 = (XNEC1+XNEC1)+(XNEC1+XNEC1) - XNEC2 - XNXC
       XC2 = (XNEC2+XNEC2)+(XNEC2+XNEC2) - XNEC1 - XNXC
       IF(XC1.EQ.ZERO.AND.XC2.EQ.ZERO) GO TO 45
C
C      ----- EXCH -----
C
        DO 40 J = NCORBS+1,MI
         J1=MAX0(J,L)-NCORBS
         L1=MIN0(J,L)-NCORBS
         NGJL = IA(J1)+L1
         XGJL = GIJ(NGJL)
         IF(XGJL.EQ.ZERO) GO TO 40
C
         NJROT=IC(N,J)
          IF(NJROT.EQ.0.OR.MIROT.EQ.0) GO TO 35
         FF = XGJL * (XC1+XC1)
         IF (N.LT.J) FF = -FF
         NRX = IA(MAX0(MIROT,NJROT))+MIN0(MIROT,NJROT)
         EH(NRX) = EH(NRX) + FF
C
  35     MJROT=IC(M+MI,J)
          IF(MJROT.EQ.0.OR.NIROT.EQ.0) GO TO 40
         FF = XGJL * (XC2+XC2)
          IF (N.LT.I) FF = -FF
         NRX = IA(MAX0(MJROT,NIROT))+MIN0(MJROT,NIROT)
         EH(NRX) = EH(NRX) + FF
C
  40   CONTINUE
  45   CONTINUE
  50  CONTINUE
  60  CONTINUE
  65  CONTINUE
C
  80  CONTINUE
C
C     ----- READ IN -DM2- FROM -NFTG- -----
C
      CALL SEQREW(NFTG)
  100 CALL PREAD(NFTG,G,IX,NG,NINTMX)
      IF(NG.EQ.0) GO TO 4100
      MGIJKL=IABS(NG)
      IF(MGIJKL.GT.NINTMX) THEN
         CALL ABRT
         STOP
      END IF
      DO 4000 MG=1,MGIJKL
      GGIJKL=G(MG)
      IF( ABS(GGIJKL).LT.GTOL) GO TO 4000
C
                       NPACK = MG
                       IF (LABSIZ .EQ. 2) THEN
                         LABEL = IX(NPACK)
                         IPACK = ISHFT( LABEL, -48 )
                         JPACK = IAND( ISHFT( LABEL, -32 ), 65535 )
                         KPACK = IAND( ISHFT( LABEL, -16 ), 65535 )
                         LPACK = IAND( LABEL, 65535 )
                       ELSE IF (LABSIZ .EQ. 1) THEN
                         IF ( MOD(NPACK,2) .EQ. 0 ) THEN
                           LABEL = IX( NPACK/2 )
                           IPACK = IAND( ISHFT( LABEL, -24 ), 255 )
                           JPACK = IAND( ISHFT( LABEL, -16 ), 255 )
                           KPACK = IAND( ISHFT( LABEL,  -8 ), 255 )
                           LPACK = IAND( LABEL, 255 )
                         ELSE
                           LABEL = IX( (NPACK/2)+1 )
                           IPACK = ISHFT( LABEL, -56 )
                           JPACK = IAND( ISHFT( LABEL, -48 ), 255 )
                           KPACK = IAND( ISHFT( LABEL, -40 ), 255 )
                           LPACK = IAND( ISHFT( LABEL, -32 ), 255 )
                         END IF
                       END IF
                       IGIJKL = IPACK
                       JGIJKL = JPACK
                       KGIJKL = KPACK
                       LGIJKL = LPACK
C
C     ----- I AND J VALENCE -------
C
      IANDJ=IGIJKL.EQ.JGIJKL
      KANDL=KGIJKL.EQ.LGIJKL
      SAME =IGIJKL.EQ.KGIJKL.AND.JGIJKL.EQ.LGIJKL
C
      IF (DROPC) THEN
         IGIJKL = IGIJKL + NCORBS
         JGIJKL = JGIJKL + NCORBS
         KGIJKL = KGIJKL + NCORBS
         LGIJKL = LGIJKL + NCORBS
      END IF
C
C     ----- LAGRANGIAN CONTRIBUTION -----
C
      GAB=GGIJKL
      IF(.NOT.IANDJ) GAB=GAB+GAB
      NAB=IA(IGIJKL)+JGIJKL
      GCD=GGIJKL
      IF(.NOT.KANDL) GCD=GCD+GCD
      NCD=IA(KGIJKL)+LGIJKL
C
      IG=IGIJKL
      JG=JGIJKL
      KG=KGIJKL
      LG=LGIJKL
      KL=NCD
      GIJKL=GCD
      NPASS=1
 1000 CONTINUE
C
      DO 1100 M=1,MA
         JKLX=(JG-1)*MIJ+KL
         XMJKL=X(M,JKLX)
         IF(XMJKL.EQ.ZERO) GO TO 1100
         MILGR=IB(M+MI)+IG
         EG(MILGR)=EG(MILGR)+ GIJKL*XMJKL
 1100 CONTINUE
C
C     ----- CONTRIBUTIONS TO ORBITAL HESSIAN -----
C
 1200 NPASS=NPASS+1
      GO TO (1300,1210,1220,1230,1300),NPASS
 1210 IF(IANDJ) GO TO 1200
      IG=JGIJKL
      JG=IGIJKL
      GO TO 1000
 1220 IF(SAME) GO TO 1300
      IG=KGIJKL
      JG=LGIJKL
      KG=IGIJKL
      LG=JGIJKL
      KL=NAB
      GIJKL=GAB
      GO TO 1000
 1230 IF(KANDL) GO TO 1300
      IG=LGIJKL
      JG=KGIJKL
      GO TO 1000
 1300 CONTINUE
C
C     ----- COULOMB CONTRIBUTION TO HESSIAN -----
C
      IG=IGIJKL
      JG=JGIJKL
      KG=KGIJKL
      LG=LGIJKL
      KL=NCD
      GMINJ=GCD
      NPASS=1
 2000 CONTINUE
C
      DO 2100 M=1,MA
      GM=GMINJ
      MIROT=IC(M+MI,IG)
      IF(MIROT.EQ.0) GO TO 2100
      DO 2050 N=1,MI
      GIJKL=GM
      IF(N -JG) 2030,2050,2040
 2030 GIJKL=-GIJKL
 2040 NJROT=IC(N,JG)
      IF(NJROT.EQ.0) GO TO 2050
      NKLX=(N-1)*MIJ+KL
      XMNKL=X(M,NKLX)
      IF(XMNKL.EQ.ZERO) GO TO 2050
      GIJKL=GIJKL+GIJKL
      MIH=MAX(MIROT,NJROT)
      NJH=MIN(MIROT,NJROT)
      MINJH=IA(MIH)+NJH
      EH(MINJH)=EH(MINJH)+ GIJKL*XMNKL
 2050 CONTINUE
 2100 CONTINUE
C
 2200 NPASS=NPASS+1
      GO TO (2300,2210,2220,2230,2300),NPASS
 2210 IF(IANDJ) GO TO 2200
      IG=JGIJKL
      JG=IGIJKL
      GO TO 2000
 2220 IF(SAME) GO TO 2300
      IG=KGIJKL
      JG=LGIJKL
      KG=IGIJKL
      LG=JGIJKL
      KL=NAB
      GMINJ=GAB
      GO TO 2000
 2230 IF(KANDL) GO TO 2300
      IG=LGIJKL
      JG=KGIJKL
      GO TO 2000
 2300 CONTINUE
C
C     ----- EXCHANGE CONTRIBUTION TO HESSIAN -----
C
      GDM2=GGIJKL+GGIJKL
      IG=IGIJKL
      JG=JGIJKL
      KG=KGIJKL
      LG=LGIJKL
      NPASS=1
 3000 CONTINUE
C
      DO 3100 M=1,MA
      GM=GDM2
      MIROT=IC(M+MI,IG)
      IF(MIROT.EQ.0) GO TO 3100
      MJ=(JG-1)*MIJ
      DO 3050 N=1,MI
      GIJKL=GM
      IF(N -KG) 3030,3050,3040
 3030 GIJKL=-GIJKL
 3040 NKROT=IC(N,KG)
      IF(NKROT.EQ.0) GO TO 3050
      NN=MAX(N,LG)
      LL=MIN(N,LG)
      NL=IA(NN)+LL
      JNLX=MJ+NL
      XMJNL=X(M,JNLX)
      IF(XMJNL.EQ.ZERO) GO TO 3050
      GIJKL=GIJKL+GIJKL
      MIH=MAX(MIROT,NKROT)
      NKH=MIN(MIROT,NKROT)
      MINKH=IA(MIH)+NKH
      EH(MINKH)=EH(MINKH)+ GIJKL*XMJNL
 3050 CONTINUE
 3100 CONTINUE
C
 3200 NPASS=NPASS+1
      GO TO (3300,3210,3220,3230,3240,3250,3260,3270,3300),NPASS
 3210 IF(IANDJ) GO TO 3200
      IG=JGIJKL
      JG=IGIJKL
      GO TO 3000
 3220 IF(KANDL) GO TO 3200
      IG=IGIJKL
      JG=JGIJKL
      KG=LGIJKL
      LG=KGIJKL
      GO TO 3000
 3230 IF(IANDJ.OR.KANDL) GO TO 3200
      IG=JGIJKL
      JG=IGIJKL
      GO TO 3000
 3240 IF(SAME) GO TO 3300
      IG=KGIJKL
      JG=LGIJKL
      KG=IGIJKL
      LG=JGIJKL
      GO TO 3000
 3250 IF(IANDJ) GO TO 3200
      KG=JGIJKL
      LG=IGIJKL
      GO TO 3000
 3260 IF(KANDL) GO TO 3200
      IG=LGIJKL
      JG=KGIJKL
      KG=IGIJKL
      LG=JGIJKL
      GO TO 3000
 3270 IF(IANDJ.OR.KANDL) GO TO 3300
      KG=JGIJKL
      LG=IGIJKL
      GO TO 3000
 3300 CONTINUE
C
 4000 CONTINUE
      IF(NG.GT.0) GO TO 100
 4100 CONTINUE
      RETURN
      END
C*MODULE MCTWO   *DECK M1LGRR
      SUBROUTINE M1LGRR(NFTI,OX,X,IX,IA,NINTMX,MI,MA,MIJ,NOSQUR)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      LOGICAL SKIP
C
      DIMENSION OX(MA,*),X(NINTMX),IX(NINTMX),IA(*)
C
      COMMON /PCKLAB/ LABSIZ
C
C     ----- READ (VIR,OCC//OCC,OCC) INTEGRALS FOR M1LGR -----
C
      DO 10 I=1,MI
   10 IA(I)=(I*I-I)/2
      MOX=MI*MIJ
      CALL VCLR(OX,1,MA*MOX)
C
C     ----- READ IN INTEGRALS FROM -NFTI- -----
C
  110 CALL PREAD(NFTI,X,IX,NX,NINTMX)
      IF(NX.EQ.0) GO TO 200
      MX=IABS(NX)
      IF(MX.GT.NINTMX) CALL ABRT
      DO 120 M=1,MX
C
                       NPACK = M
                       IF (LABSIZ .EQ. 2) THEN
                         LABEL = IX(NPACK)
                         IPACK = ISHFT( LABEL, -48 )
                         JPACK = IAND( ISHFT( LABEL, -32 ), 65535 )
                         KPACK = IAND( ISHFT( LABEL, -16 ), 65535 )
                         LPACK = IAND( LABEL, 65535 )
                       ELSE IF (LABSIZ .EQ. 1) THEN
                         IF ( MOD(NPACK,2) .EQ. 0 ) THEN
                           LABEL = IX( NPACK/2 )
                           IPACK = IAND( ISHFT( LABEL, -24 ), 255 )
                           JPACK = IAND( ISHFT( LABEL, -16 ), 255 )
                           KPACK = IAND( ISHFT( LABEL,  -8 ), 255 )
                           LPACK = IAND( LABEL, 255 )
                         ELSE
                           LABEL = IX( (NPACK/2)+1 )
                           IPACK = ISHFT( LABEL, -56 )
                           JPACK = IAND( ISHFT( LABEL, -48 ), 255 )
                           KPACK = IAND( ISHFT( LABEL, -40 ), 255 )
                           LPACK = IAND( ISHFT( LABEL, -32 ), 255 )
                         END IF
                       END IF
                       K = IPACK
                       L = JPACK
                       I = KPACK
                       J = LPACK
C
         IF(NOSQUR.EQ.1) THEN
            IDUM=I
            JDUM=J
            I=K
            J=L
            K=IDUM
            L=JDUM
         END IF
         SKIP=(I.LE.MI).OR.(J.GT.MI).OR.(K.GT.MI).OR.(L.GT.MI)
         IF(SKIP) GO TO 120
         KL=IA(K)+L
         JKL=(J-1)*MIJ+KL
         OX(I-MI,JKL)=X(M)
  120 CONTINUE
      IF(NX.GT.0) GO TO 110
  200 CONTINUE
      RETURN
      END
C*MODULE MCTWO   *DECK M2DM2R
      SUBROUTINE M2DM2R(NFTG,NORBS,OG,G,IX,NINTMX,IA)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      DIMENSION OG(*),G(NINTMX),IX(NINTMX),IA(*)
C
      COMMON /PCKLAB/ LABSIZ
C
C     ----- READ 2 PARTICLE DENSITY MATRIX FOR USE BY M2TEI -----
C
      NORB2=(NORBS*(NORBS+1))/2
      NORB4=(NORB2*(NORB2+1))/2
      DO 10 I=1,NORB2
   10 IA(I)=(I*(I-1))/2
      CALL VCLR(OG,1,NORB4)
C
C     ----- READ IN -DM2- FROM -NFTG- -----
C
      CALL SEQREW(NFTG)
  110 CALL PREAD(NFTG,G,IX,NG,NINTMX)
      IF(NG.EQ.0) GO TO 200
      MG=IABS(NG)
      IF(MG.GT.NINTMX) THEN
         CALL ABRT
         STOP
      END IF
      DO 120 M=1,MG
C
                       NPACK = M
                       IF (LABSIZ .EQ. 2) THEN
                         LABEL = IX(NPACK)
                         IPACK = ISHFT( LABEL, -48 )
                         JPACK = IAND( ISHFT( LABEL, -32 ), 65535 )
                         KPACK = IAND( ISHFT( LABEL, -16 ), 65535 )
                         LPACK = IAND( LABEL, 65535 )
                       ELSE IF (LABSIZ .EQ. 1) THEN
                         IF ( MOD(NPACK,2) .EQ. 0 ) THEN
                           LABEL = IX( NPACK/2 )
                           IPACK = IAND( ISHFT( LABEL, -24 ), 255 )
                           JPACK = IAND( ISHFT( LABEL, -16 ), 255 )
                           KPACK = IAND( ISHFT( LABEL,  -8 ), 255 )
                           LPACK = IAND( LABEL, 255 )
                         ELSE
                           LABEL = IX( (NPACK/2)+1 )
                           IPACK = ISHFT( LABEL, -56 )
                           JPACK = IAND( ISHFT( LABEL, -48 ), 255 )
                           KPACK = IAND( ISHFT( LABEL, -40 ), 255 )
                           LPACK = IAND( ISHFT( LABEL, -32 ), 255 )
                         END IF
                       END IF
                       I = IPACK
                       J = JPACK
                       K = KPACK
                       L = LPACK
C
         IF(I.GT.NORBS) GO TO 200
         IJ=IA(I)+J
         KL=IA(K)+L
         IJKL=IA(IJ)+KL
         OG(IJKL)=G(M)
  120 CONTINUE
      IF(NG.GT.0) GO TO 110
  200 CONTINUE
      RETURN
      END
C*MODULE MCTWO   *DECK M2TEI
      SUBROUTINE M2TEI(NE,NH2,EH,G,NFTX,X,IX,
     *                 IA,IB,IC,NINTMX,NIA,NORB,NORBS)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      LOGICAL IANDJ,KANDL,SAME
C
      DIMENSION EH(*),IA(NIA),IB(NORB),IC(NORB,NORBS)
      DIMENSION G(*),X(NINTMX),IX(NINTMX)
C
      COMMON /PCKLAB/ LABSIZ
C
      PARAMETER (TOL=1.0D-09)
C
      DO 10 I=1,NIA
   10 IA(I)=(I*I-I)/2
      NEH=NE+NH2
      CALL VCLR(EH,1,NEH)
C
C     ----- READ IN -INT2- FROM -NFTX- -----
C
  100 CALL PREAD(NFTX,X,IX,NX,NINTMX)
      IF(NX.EQ.0) GO TO 4100
      MXIJKL=IABS(NX)
      IF(MXIJKL.GT.NINTMX) THEN
         CALL ABRT
         STOP
      END IF
      DO 4000 MX=1,MXIJKL
         XMX=X(MX)
C
                       NPACK = MX
                       IF (LABSIZ .EQ. 2) THEN
                         LABEL = IX(NPACK)
                         IPACK = ISHFT( LABEL, -48 )
                         JPACK = IAND( ISHFT( LABEL, -32 ), 65535 )
                         KPACK = IAND( ISHFT( LABEL, -16 ), 65535 )
                         LPACK = IAND( LABEL, 65535 )
                       ELSE IF (LABSIZ .EQ. 1) THEN
                         IF ( MOD(NPACK,2) .EQ. 0 ) THEN
                           LABEL = IX( NPACK/2 )
                           IPACK = IAND( ISHFT( LABEL, -24 ), 255 )
                           JPACK = IAND( ISHFT( LABEL, -16 ), 255 )
                           KPACK = IAND( ISHFT( LABEL,  -8 ), 255 )
                           LPACK = IAND( LABEL, 255 )
                         ELSE
                           LABEL = IX( (NPACK/2)+1 )
                           IPACK = ISHFT( LABEL, -56 )
                           JPACK = IAND( ISHFT( LABEL, -48 ), 255 )
                           KPACK = IAND( ISHFT( LABEL, -40 ), 255 )
                           LPACK = IAND( ISHFT( LABEL, -32 ), 255 )
                         END IF
                       END IF
                       IXIJKL = IPACK
                       JXIJKL = JPACK
                       KXIJKL = KPACK
                       LXIJKL = LPACK
C
      IANDJ=IXIJKL.EQ.JXIJKL
      KANDL=KXIJKL.EQ.LXIJKL
      SAME =IXIJKL.EQ.KXIJKL.AND.JXIJKL.EQ.LXIJKL
C
      XAB=XMX
      IF(.NOT.IANDJ) XAB=XAB+XAB
      NAB=IA(IXIJKL)+JXIJKL
      XCD=XMX
      IF(.NOT.KANDL) XCD=XCD+XCD
      NCD=IA(KXIJKL)+LXIJKL
C
C     ----- LAGRANGIAN CONTRIBUTION -----
C
      IG=IXIJKL
      JG=JXIJKL
      KG=KXIJKL
      LG=LXIJKL
      KL=NCD
      XIJKL=XCD
      NPASS=1
 1000 CONTINUE
      IF(JG.GT.NORBS.OR.KG.GT.NORBS.OR.LG.GT.NORBS) GO TO 1200
C
      DO 1100 M=1,NORBS
         MM=MAX(M,JG)
         JJ=MIN(M,JG)
         MJ=IA(MM)+JJ
         MJG=MAX(MJ,KL)
         KLG=MIN(MJ,KL)
         MJKLG=IA(MJG)+KLG
         GMJKL=G(MJKLG)
         IF( ABS(GMJKL).LT.TOL) GO TO 1100
         MILGR=IB(IG)+M
         EH(MILGR)=EH(MILGR)+XIJKL*GMJKL
 1100 CONTINUE
C
 1200 NPASS=NPASS+1
      GO TO (1300,1210,1220,1230,1300),NPASS
 1210 IF(IANDJ) GO TO 1200
      IG=JXIJKL
      JG=IXIJKL
      GO TO 1000
 1220 IF(SAME) GO TO 1300
      IG=KXIJKL
      JG=LXIJKL
      KG=IXIJKL
      LG=JXIJKL
      KL=NAB
      XIJKL=XAB
      GO TO 1000
 1230 IF(KANDL) GO TO 1300
      IG=LXIJKL
      JG=KXIJKL
      GO TO 1000
 1300 CONTINUE
C
C     ----- COULOMB CONTRIBUTION TO HESSIAN -----
C
      IG=IXIJKL
      JG=JXIJKL
      KG=KXIJKL
      LG=LXIJKL
      KL=NCD
      XMINJ=XCD
      NPASS=1
 2000 CONTINUE
      IF(KG.GT.NORBS.OR.LG.GT.NORBS) GO TO 2200
C
      DO 2100 M=1,NORBS
      XM=XMINJ
      IF(M -IG) 2020,2100,2010
 2010 XM=-XM
 2020 MIROT=IC(IG,M)
      IF(MIROT.EQ.0) GO TO 2100
      MAXN=M
      IF(JG.GT.IG) MAXN=M-1
      IF(MAXN.EQ.0) GO TO 2100
      DO 2050 N=1,MAXN
      XIJKL=XM
      IF(N -JG) 2040,2050,2030
 2030 XIJKL=-XIJKL
 2040 NJROT=IC(JG,N)
      IF(NJROT.EQ.0) GO TO 2050
      IF(N.NE.M.OR.JG.NE.IG) XIJKL=XIJKL+XIJKL
      MN=IA(M)+N
      MNG=MAX(MN,KL)
      KLG=MIN(MN,KL)
      MNKLG=IA(MNG)+KLG
      GMNKL=G(MNKLG)
      IF( ABS(GMNKL).LT.TOL) GO TO 2050
      MIH=MAX(MIROT,NJROT)
      NJH=MIN(MIROT,NJROT)
      MINJH=IA(MIH)+NJH
      MINJH=MINJH+NE
      EH(MINJH)=EH(MINJH)+XIJKL*GMNKL
 2050 CONTINUE
 2100 CONTINUE
C
 2200 NPASS=NPASS+1
      GO TO (2300,2210,2220,2230,2300),NPASS
 2210 IF(IANDJ) GO TO 2200
      IG=JXIJKL
      JG=IXIJKL
      GO TO 2000
 2220 IF(SAME) GO TO 2300
      IG=KXIJKL
      JG=LXIJKL
      KG=IXIJKL
      LG=JXIJKL
      KL=NAB
      XMINJ=XAB
      GO TO 2000
 2230 IF(KANDL) GO TO 2300
      IG=LXIJKL
      JG=KXIJKL
      GO TO 2000
 2300 CONTINUE
C
C     ----- EXCHANGE CONTRIBUTION TO HESSIAN -----
C
      XINT2=XMX+XMX
      IG=IXIJKL
      JG=JXIJKL
      KG=KXIJKL
      LG=LXIJKL
      NPASS=1
 3000 CONTINUE
      IF(JG.GT.NORBS.OR.LG.GT.NORBS) GO TO 3200
C
      DO 3100 M=1,NORBS
      XM=XINT2
      IF(M -IG) 3020,3100,3010
 3010 XM=-XM
 3020 MIROT=IC(IG,M)
      IF(MIROT.EQ.0) GO TO 3100
      MM=MAX(M,JG)
      JJ=MIN(M,JG)
      MJ=IA(MM)+JJ
      MAXN=M
      IF(KG.GT.IG) MAXN=M-1
      IF(MAXN.EQ.0) GO TO 3100
      DO 3050 N=1,MAXN
      XIJKL=XM
      IF(N -KG) 3040,3050,3030
 3030 XIJKL=-XIJKL
 3040 NKROT=IC(KG,N)
      IF(NKROT.EQ.0) GO TO 3050
      IF(N.NE.M.OR.KG.NE.IG) XIJKL=XIJKL+XIJKL
      NN=MAX(N,LG)
      LL=MIN(N,LG)
      NL=IA(NN)+LL
      MJG=MAX(MJ,NL)
      NLG=MIN(MJ,NL)
      MJNLG=IA(MJG)+NLG
      GMJNL=G(MJNLG)
      IF( ABS(GMJNL).LT.TOL) GO TO 3050
      MIH=MAX(MIROT,NKROT)
      NKH=MIN(MIROT,NKROT)
      MINKH=IA(MIH)+NKH
      MINKH=MINKH+NE
      EH(MINKH)=EH(MINKH)+XIJKL*GMJNL
 3050 CONTINUE
 3100 CONTINUE
C
 3200 NPASS=NPASS+1
      GO TO (3300,3210,3220,3230,3240,3250,3260,3270,3300),NPASS
 3210 IF(IANDJ) GO TO 3200
      IG=JXIJKL
      JG=IXIJKL
      GO TO 3000
 3220 IF(KANDL) GO TO 3200
      IG=IXIJKL
      JG=JXIJKL
      KG=LXIJKL
      LG=KXIJKL
      GO TO 3000
 3230 IF(IANDJ.OR.KANDL) GO TO 3200
      IG=JXIJKL
      JG=IXIJKL
      GO TO 3000
 3240 IF(SAME) GO TO 3300
      IG=KXIJKL
      JG=LXIJKL
      KG=IXIJKL
      LG=JXIJKL
      GO TO 3000
 3250 IF(IANDJ) GO TO 3200
      KG=JXIJKL
      LG=IXIJKL
      GO TO 3000
 3260 IF(KANDL) GO TO 3200
      IG=LXIJKL
      JG=KXIJKL
      KG=IXIJKL
      LG=JXIJKL
      GO TO 3000
 3270 IF(IANDJ.OR.KANDL) GO TO 3300
      KG=JXIJKL
      LG=IXIJKL
      GO TO 3000
 3300 CONTINUE
C
 4000 CONTINUE
      IF(NX.GT.0) GO TO 100
 4100 CONTINUE
      RETURN
      END
C
C*MODULE MCTWO   *DECK FVADDI
      SUBROUTINE FVADDI(NMOS,NCOR,NACT,DVAL,FVAL,BUFF)
      use mx_limits, only: mxao
C
C ----------------------------------------------------------------------
C
C  NEWTON-RAPHSON MCSCF. CALLED FROM NEWTON (MCSCF.SRC).
C  CONSTRUCT THE VALENCE FOCK OPERATOR DIRECTLY IN THE MO BASIS
C  FROM TRANSFORMED INTEGRALS STORED IN THE DDI ARRAYS.
C
C  SYMBOLS:
C           NMOS = TOTAL NUMBER OF MOS
C           NCOR = NUMBER OF CORE MOS
C           NACT = NUMBER OF ACTIVE MOS
C           DVAL = 1-EL DENSITY SPANNING ACTIVE INDICES
C           FVAL = VALENCE FOCK OPERATOR IN THE MO BASIS
C           BUFF = A MESSAGE BUFFER
C
C  INTEGRALS NEED AT LEAST TWO ACTIVE INDICES IN THE FORMULA
C  FVAL(IJ) = SUM_OVER_KL{DVAL(KL)*[(IJ|KL)-1/2(IK|JL)]}, K,L ACTIVE,
C  I,J GENERAL. FVAL,DVAL STORED IN TRIANGULAR FORM.
C
C ----------------------------------------------------------------------
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (ZERO=0.0D+00)
      PARAMETER (HALF=0.5D+00)
      LOGICAL INJ
      DOUBLE PRECISION  FVAL(*),DVAL(*),BUFF(*)
      INTEGER DDI_NP, DDI_ME, A,B,AB,BA,AI
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
      NOCC = NCOR + NACT
      NVIR = NMOS - NOCC
      NMTR = (NMOS*NMOS+NMOS)/2
      NOTR = (NOCC*NOCC+NOCC)/2
      NVTR = (NVIR*NVIR+NVIR)/2
      NVSQ = NVIR*NVIR
      CALL DCOPY(NMTR,ZERO,0,FVAL,1)
      CALL DDI_NPROC(DDI_NP, DDI_ME)
C
C  (OO|OO) CLASS
C
      CALL DDI_DISTRIB(D_OOOO,DDI_ME,ILO,IHI,JLO,JHI)
      DO I = 1, NOCC
        DO J = 1, I
          IJCOL = IA(I) + J
          IF ((IJCOL.GE.JLO).AND.(IJCOL.LE.JHI)) THEN
            CALL DDI_GET(D_OOOO,1,NOTR,IJCOL,IJCOL,BUFF)
            IJ = IJCOL
C
C  COULOMB TERMS
C
            DO K = 1, NACT
              KN = K + NCOR
              DO L = 1, NACT
                LN = L + NCOR
                KL   = IA(MAX(K,L))   + MIN(K,L)
                KLN  = IA(MAX(KN,LN)) + MIN(KN,LN)
                FVAL(IJ) = FVAL(IJ)   + DVAL(KL)*BUFF(KLN)
              END DO
            END DO
C
C  EXCHANGE TERMS
C
            IF (I.GT.NCOR) THEN
              IN = I - NCOR
              DO K = NCP1, NOCC
                KN = K - NCOR
                IKN = IA(MAX(IN,KN)) + MIN(IN,KN)
                DIK = HALF*DVAL(IKN)
                DO L = 1, J
                  JL = IA(J) + L
                  KL = IA(MAX(K,L)) + MIN(K,L)
                  FVAL(JL) = FVAL(JL) - DIK*BUFF(KL)
                END DO
              END DO
              IF ((J.GT.NCOR).AND.(J.NE.I)) THEN
                JN = J - NCOR
                DO K = NCP1, NOCC
                  KN = K - NCOR
                  JKN = IA(MAX(JN,KN)) + MIN(JN,KN)
                  DJK = HALF*DVAL(JKN)
                  DO L = 1, I
                    IL = IA(I) + L
                    KL = IA(MAX(K,L)) + MIN(K,L)
                    FVAL(IL) = FVAL(IL) - DJK*BUFF(KL)
                  END DO
                END DO
              END IF
            END IF
          END IF  ! LOCAL
        END DO  ! J
      END DO  ! I
C
C  (VO|OO) CLASS
C
      CALL DDI_DISTRIB(D_VOOO,DDI_ME,ILO,IHI,JLO,JHI)
C
C  COULOMB TERMS (VO|AA)
C
      DO I = NCP1, NOCC
        IN = I - NCOR
        DO J = NCP1, I
          JN = J - NCOR
          IJN = IA(IN) + JN
          DIJ = DVAL(IJN)
          IF (I.NE.J) DIJ = DIJ + DIJ
          IJCOL = IA(I) + J
          DO K = 1, NOCC
            IJKCOL = (IJCOL-1)*NOCC + K
            IF ((IJKCOL.GE.JLO).AND.(IJKCOL.LE.JHI)) THEN
              CALL DDI_GET(D_VOOO,1,NVIR,IJKCOL,IJKCOL,BUFF)
              DO A = 1, NVIR
                NA = A + NOCC
                KA = IA(NA) + K
                FVAL(KA) = FVAL(KA) + DIJ*BUFF(A)
              END DO
            END IF  ! LOCAL
          END DO  ! K
        END DO  ! J
      END DO  ! I
C
C  EXCHANGE TERMS (VA|AO)
C
      DO I = 1, NOCC
        DO J = NCP1, NOCC
          JN = J - NCOR
          IJCOL = IA(MAX(I,J)) + MIN(I,J)
          DO K = NCP1, NOCC
            KN = K - NCOR
            IJKCOL = (IJCOL-1)*NOCC + K
            IF ((IJKCOL.GE.JLO).AND.(IJKCOL.LE.JHI)) THEN
              CALL DDI_GET(D_VOOO,1,NVIR,IJKCOL,IJKCOL,BUFF)
              JKN = IA(MAX(JN,KN)) + MIN(JN,KN)
              DJK = HALF*DVAL(JKN)
              DO A = 1, NVIR
                NA = A + NOCC
                AI = IA(NA) + I
                FVAL(AI) = FVAL(AI) - DJK*BUFF(A)
              END DO
            END IF  ! LOCAL
          END DO  ! K
        END DO  ! J
      END DO  ! I
C
C  (VV|AA) CLASS - COULOMB TERMS
C
      CALL DDI_DISTRIB(D_VVOO,DDI_ME,ILO,IHI,JLO,JHI)
      DO I = 1, NACT
        IN = I + NCOR
        DO J = 1, I
          IJ = IA(I) + J
          JN = J + NCOR
          IJCOL = IA(IN) + JN
          IF ((IJCOL.GE.JLO).AND.(IJCOL.LE.JHI)) THEN
            CALL DDI_GET(D_VVOO,1,NVTR,IJCOL,IJCOL,BUFF)
            DIJ = DVAL(IJ)
            IF (I.NE.J) DIJ = DIJ + DIJ
            DO A = 1, NVIR
              NA = A + NOCC
              DO B = 1, A
                AB   = IA(A) + B
                NB   = B + NOCC
                NAB  = IA(NA) + NB
                FVAL(NAB) = FVAL(NAB) + DIJ*BUFF(AB)
              END DO
            END DO
          END IF  ! LOCAL
        END DO  ! J
      END DO  ! I
C
C  (VA|VA) CLASS - EXCHANGE TERMS
C
      CALL DDI_DISTRIB(D_VOVO,DDI_ME,ILO,IHI,JLO,JHI)
      DO I = 1, NACT
        IN = I + NCOR
        DO J = 1, I
          IJ = IA(I) + J
          JN = J + NCOR
          IJCOL = IA(IN) + JN
          INJ = I.NE.J
          IF ((IJCOL.GE.JLO).AND.(IJCOL.LE.JHI)) THEN
            CALL DDI_GET(D_VOVO,1,NVSQ,IJCOL,IJCOL,BUFF)
            DIJ = HALF*DVAL(IJ)
            DO A = 1, NVIR
              NA = A + NOCC
              DO B = 1, A
                AB   = (A-1)*NVIR + B
                NB   = B + NOCC
                NAB  = IA(NA) + NB
                FVAL(NAB) = FVAL(NAB) - DIJ*BUFF(AB)
              END DO
            END DO
            IF (INJ) THEN
              DO A = 1, NVIR
                NA = A + NOCC
                DO B = 1, A
                  BA   = (B-1)*NVIR + A
                  NB   = B + NOCC
                  NAB  = IA(NA) + NB
                  FVAL(NAB) = FVAL(NAB) - DIJ*BUFF(BA)
                END DO
              END DO
            END IF
          END IF  ! LOCAL
        END DO  ! J
      END DO  ! I
C
C  GLOBALLY SUM THE VALENCE FOCK MATRIX
C
      CALL DDI_GSUMF( 2520, FVAL, NMTR )
      RETURN
      END
C*MODULE MCTWO   *DECK RD2PDM
      SUBROUTINE RD2PDM(INFILE,TPDM,GBUF,LABS,NINTMX)
      use mx_limits, only: mxao
C
C -----------------------------------------------------------------
C
C  READ 2-EL DENSITY SPANNING ACTIVE INDICES. WRITTEN BY WTDM12
C  (ALDECI.SRC). NOTE LABELS START AT 1 FOR FIRST ACTIVE MO.
C
C  SYMBOLS:
C           INFILE = DISC FILE IDENTIFIER
C           TPDM   = 2-EL DENSITY SPANNING ACTIVE INDICES
C           GBUF   = DENSITY MATRIX ELEMENT I/O BUFFER
C           LABS   = I/O BUFFER FOR 2-PDM ACTIVE LABELS
C           NINTMX = I/O BUFFER LENGTH IN WORDS
C
C  NOTE: USE OF IA() WITH MXAO=8192 ASSUMES NUMBER OF
C        ACTIVE ORBITALS IS LESS THAN APPROX. 63.
C
C -----------------------------------------------------------------
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      LOGICAL READMORE
      INTEGER LABS(*)
      DOUBLE PRECISION  TPDM(*),GBUF(*)
      COMMON /PCKLAB/ LABSIZ
      COMMON /IJPAIR/ IA(MXAO)
C
      CALL SEQREW(INFILE)
      READMORE = .TRUE.
      DO WHILE (READMORE)
        NG = 0
        CALL PREAD(INFILE,GBUF,LABS,NG,NINTMX)
        IF (NG.LE.0) READMORE = .FALSE.
        MG = IABS(NG)
        DO M = 1, MG
C
C  UNPACK LABELS
C
                       NPACK = M
                       IF (LABSIZ .EQ. 2) THEN
                         LABEL = LABS(NPACK)
                         IPACK = ISHFT( LABEL, -48 )
                         JPACK = IAND( ISHFT( LABEL, -32 ), 65535 )
                         KPACK = IAND( ISHFT( LABEL, -16 ), 65535 )
                         LPACK = IAND( LABEL, 65535 )
                       ELSE IF (LABSIZ .EQ. 1) THEN
                         IF ( MOD(NPACK,2) .EQ. 0 ) THEN
                           LABEL = LABS( NPACK/2 )
                           IPACK = IAND( ISHFT( LABEL, -24 ), 255 )
                           JPACK = IAND( ISHFT( LABEL, -16 ), 255 )
                           KPACK = IAND( ISHFT( LABEL,  -8 ), 255 )
                           LPACK = IAND( LABEL, 255 )
                         ELSE
                           LABEL = LABS( (NPACK/2)+1 )
                           IPACK = ISHFT( LABEL, -56 )
                           JPACK = IAND( ISHFT( LABEL, -48 ), 255 )
                           KPACK = IAND( ISHFT( LABEL, -40 ), 255 )
                           LPACK = IAND( ISHFT( LABEL, -32 ), 255 )
                         END IF
                       END IF
                       I = IPACK
                       J = JPACK
                       K = KPACK
                       L = LPACK
          IJ = IA(MAX(I,J)) + MIN(I,J)
          KL = IA(MAX(K,L)) + MIN(K,L)
          IJKL = IA(MAX(IJ,KL)) + MIN(IJ,KL)
          TPDM(IJKL) = GBUF(M)
        END DO  ! CURRENT BATCH
      END DO  ! NEXT BATCH
      RETURN
      END
C*MODULE MCTWO   *DECK CHNRGY
      SUBROUTINE CHNRGY(NCOR,NACT,HAMO,OPDM,ECOR)
      use mx_limits, only: mxao
C
C -----------------------------------------------------------------
C
C  COMPUTE THE 1-EL (OR CORE-HAMILTONIAN) ENERGY
C
C  SYMBOLS:
C           NCOR = NUMBER OF CORE MOS
C           NACT = NUMBER OF ACTIVE MOS
C           HAMO = 1-EL CORE HAMILTONIAN INTEGRALS
C           OPDM = 1-EL DENSITY SPANNING ACTIVE INDICES
C           ECOR = 1-EL 'CORE' ENERGY
C
C -----------------------------------------------------------------
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (ZERO=0.0D+00)
      PARAMETER (TWO=2.0D+00)
      DOUBLE PRECISION HAMO(*),OPDM(*)
      COMMON /IJPAIR/ IA(MXAO)
C
      ECOR = ZERO
      II = 0
      DO I = 1, NCOR
        II = II + I
        ECOR = ECOR + TWO*HAMO(II)
      END DO
C
      DO I = 1, NACT
        IN = I + NCOR
        DO J = 1, NACT
          JN   = J + NCOR
          IJ   = IA(MAX(I,J))   + MIN(I,J)
          IJN  = IA(MAX(IN,JN)) + MIN(IN,JN)
          ECOR = ECOR + OPDM(IJ)*HAMO(IJN)
        END DO
      END DO
      RETURN
      END
C*MODULE MCTWO   *DECK LAGDDI
      SUBROUTINE LAGDDI(NMOS,NCOR,NOCC
     *,                 OPDM,TPDM,FCOR,FVAL
     *,                 BUFF,LAGN,ENERGY)
      use mx_limits, only: mxao
C
C -----------------------------------------------------------------
C
C  NEWTON-RAPHSON MCSCF. CALLED FROM SUBROUTINE NEWTON.
C  COMPUTE THE LAGRANGIAN MATRIX AND THE MCSCF ENERGY.
C  SEE MOTECC-90 (ESCOM), PAGE 293, EQNS (B.71-75).
C  TO MAKE USE OF PARALLEL TRANSFORMATION (TRANDDI), CORE AND
C  ACTIVE MO INDICES ARE SUBSETS OF THE OCCUPIED MO LIST.
C  TRANDDI IS CALLED FROM SUBROUTINE TRFMCX (IN MODULE TRANS.SRC).
C
C  SYMBOLS:
C     NMOS = TOTAL NUMBER OF MOS
C     NCOR = NUMBER OF CORE MOS
C     NOCC = NUMBER OF OCCUPIED MOS
C     OPDM = 1-EL DENSITY SPANNING ACTIVE INDICES
C     TPDM = 2-EL DENSITY SPANNING ACTIVE INDICES
C     FCOR = CORE FOCK OPERATOR
C     FVAL = VALENCE FOCK OPERATOR
C     BUFF = MESSAGE BUFFER
C     LAGN = LAGRANGIAN MATRIX  [OUTPUT]
C     ENERGY  = MCSCF ENERGY [OUTPUT]
C
C -----------------------------------------------------------------
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (ZERO=0.0D+00)
      PARAMETER (HALF=0.5D+00)
      PARAMETER (TWO=2.0D+00)
      LOGICAL INJ
      INTEGER DDI_NP,DDI_ME, A,AN
      DOUBLE PRECISION  OPDM(*),TPDM(*),FCOR(*),FVAL(*)
      DOUBLE PRECISION  LAGN(NOCC,*),BUFF(*)
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
      CALL DDI_NPROC(DDI_NP,DDI_ME)
C
C  1-EL TERMS OF THE LAGRANGIAN
C
      CALL DCOPY(NMOS*NOCC,ZERO,0,LAGN,1)
      DO M = 1, NMOS
        IF (MOD(M,DDI_NP).EQ.DDI_ME) THEN
C
C  1. L(MI) = 2*(FC(MI) + FV(MI)),  I CORE, M GENERAL
C
          DO I = 1, NCOR
            IM = IA(MAX(I,M)) + MIN(I,M)
            FF = FCOR(IM) + FVAL(IM)
            LAGN(I,M) = LAGN(I,M) + TWO*FF
          END DO
C
C  2. L(MI) = SUM_J  D(IJ)*FC(MJ),  I,J VALENCE, M GENERAL
C
          DO I = NCP1, NOCC
            ID = I - NCOR
            DO J = NCP1, NOCC
              JD = J - NCOR
              IJD = IA(MAX(ID,JD)) + MIN(ID,JD)
              MJ = IA(MAX(M,J)) + MIN(M,J)
              LAGN(I,M) = LAGN(I,M) + FCOR(MJ)*OPDM(IJD)
            END DO
          END DO
C
        END IF  !  PARALLEL
      END DO
C
C  2-EL TERMS OF THE LAGRANGIAN
C
C ----------- CONTRIBUTIONS FROM (OO|OO) TYPE INTEGRALS -----------
C
      CALL DDI_DISTRIB(D_OOOO,DDI_ME,ILO,IHI,JLO,JHI)
C
C  (ACTIVE-CORE|**) TERMS
C
      DO I = NCP1, NOCC
        DO J = 1, NCOR
          IJ = IA(MAX(I,J)) + MIN(I,J)
          IF ((IJ.GE.JLO).AND.(IJ.LE.JHI)) THEN
            CALL DDI_GET(D_OOOO,1,NOTR,IJ,IJ,BUFF)
C
C  1. L(JM) <- D(IM|KL)*(IJ|KL), J CORE
C
            ID = I - NCOR
            DO K = NCP1, NOCC
              KD = K - NCOR
              DO L = NCP1, NOCC
                LD = L - NCOR
                KL  = IA(MAX(K ,L )) + MIN(K ,L )
                KLD = IA(MAX(KD,LD)) + MIN(KD,LD)
                ERI = BUFF(KL)
                DO M = NCP1, NOCC
                  MD = M - NCOR
                  IMD = IA(MAX(ID,MD)) + MIN(ID,MD)
                  IMKLD = IA(MAX(IMD,KLD)) + MIN(IMD,KLD)
                  LAGN(M,J) = LAGN(M,J) + TPDM(IMKLD)*ERI
                END DO
              END DO
            END DO
          END IF  ! LOCAL
        END DO  ! J
      END DO  ! I
C
C  (ACTIVE-ACTIVE|**) TERMS
C
      DO I = NCP1, NOCC
        DO J = NCP1, I
          IJ = IA(MAX(I,J)) + MIN(I,J)
          IF ((IJ.GE.JLO).AND.(IJ.LE.JHI)) THEN
            CALL DDI_GET(D_OOOO,1,NOTR,IJ,IJ,BUFF)
C
C  2. L(KM) <- D(IJ|LM)*(IJ|KL)
C
            INJ = I.NE.J
            ID = I - NCOR
            JD = J - NCOR
            IJD = IA(MAX(ID,JD)) + MIN(ID,JD)
            DO K = NCP1, NOCC
              DO L = NCP1, NOCC
                LD = L - NCOR
                KL = IA(MAX(K,L)) + MIN(K,L)
                ERI = BUFF(KL)
                IF (INJ) ERI = ERI*TWO
                DO M = NCP1, NOCC
                  MD = M - NCOR
                  LMD = IA(MAX(LD,MD)) + MIN(LD,MD)
                  IJLMD = IA(MAX(IJD,LMD)) + MIN(IJD,LMD)
                  LAGN(M,K) = LAGN(M,K) + TPDM(IJLMD)*ERI
                END DO
              END DO
            END DO
          END IF  ! LOCAL
        END DO  ! J
      END DO  ! I
C
C ----------- CONTRIBUTIONS FROM (VO|OO) TYPE INTEGRALS -----------
C
      CALL DDI_DISTRIB(D_VOOO,DDI_ME,ILO,IHI,JLO,JHI)
C
C  (ACTIVE-ACTIVE|ACTIVE-VIRTUAL)
C
      DO I = NCP1, NOCC
        DO J = NCP1, I
          IJ = IA(MAX(I,J)) + MIN(I,J)
          DO K = NCP1, NOCC
            IJK = (IJ-1)*NOCC + K
            IF ((IJK.GE.JLO).AND.(IJK.LE.JHI)) THEN
              CALL DDI_GET(D_VOOO,1,NVIR,IJK,IJK,BUFF)
C
C  3. L(AL) <- D(IJ|KL)*(AK|IJ)
C
              INJ = I.NE.J
              ID = I - NCOR
              JD = J - NCOR
              KD = K - NCOR
              IJD = IA(MAX(ID,JD)) + MIN(ID,JD)
              DO A = 1, NVIR
                ERI = BUFF(A)
                IF (INJ) ERI = ERI*TWO
                AN = A + NOCC
                DO L = NCP1, NOCC
                  LD = L - NCOR
                  KLD = IA(MAX(KD,LD)) + MIN(KD,LD)
                  IJKLD = IA(MAX(IJD,KLD)) + MIN(IJD,KLD)
                  LAGN(L,AN) = LAGN(L,AN) + TPDM(IJKLD)*ERI
                END DO
              END DO
            END IF  ! LOCAL
          END DO  ! K
        END DO  ! J
      END DO  ! I
C
C  GLOBALLY SUM LAGRANGIAN
C
      CALL DDI_GSUMF(2290,LAGN,NMOS*NOCC)
C
C  COMPUTE MCSCF ENERGY
C
      ENERGY = ZERO
      DO I = 1, NOCC
        ENERGY = ENERGY + LAGN(I,I)
      END DO
      ENERGY = ENERGY*HALF
      RETURN
      END
C*MODULE MCTWO   *DECK DIADDI
      SUBROUTINE DIADDI(NMOS,NCOR,NOCC,NROT
     *,                 OPDM,TPDM,FCOR,FVAL
     *,                 LAGN,IROT,BUFF,DIAH)
      use mx_limits, only: mxao
C
C -----------------------------------------------------------------
C
C  NEWTON-RAPHSON MCSCF. CALLED FROM SUBROUTINE NEWTON.
C  COMPUTE THE DIAGONAL ELEMENTS OF THE AUGMENTED HESSIAN MATRIX
C  SO THAT WE CAN APPLY LEVEL SHIFTING TO THEM.
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
C     NROT = NUMBER OF ROTATIONS TO BE OPTIMIZED
C     OPDM = 1-EL DENSITY SPANNING ACTIVE INDICES
C     TPDM = 2-EL DENSITY SPANNING ACTIVE INDICES
C     FCOR = CORE FOCK OPERATOR
C     FVAL = VALENCE FOCK OPERATOR
C     LAGN = LAGRANGIAN MATRIX
C     IROT = INDEX OF ROTATIONS TO BE OPTIMIZED
C     BUFF = MESSAGE BUFFER
C     DIAH = DIAGONAL ELEMENTS OF THE AUGMENTED HESSIAN [OUTPUT]
C
C -----------------------------------------------------------------
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (ZERO=0.0D+00)
      PARAMETER (TWO=2.0D+00)
      PARAMETER (FOUR=4.0D+00)
      PARAMETER (SIX=6.0D+00)
      PARAMETER (HSHIFT=0.5D+00)
      LOGICAL INJ
      INTEGER IROT(NOCC,*), DDI_NP,DDI_ME, A,AA,AR
      DOUBLE PRECISION  OPDM(*),TPDM(*),FCOR(*),FVAL(*)
      DOUBLE PRECISION  LAGN(NOCC,*),BUFF(*),DIAH(*)
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
      CALL DCOPY(NROT,ZERO,0,DIAH,1)
C
C  1-EL TERMS OF THE AUGMENTED HESSIAN
C
      DO I = 1, NMOS
        IF (MOD(I,DDI_NP).EQ.DDI_ME) THEN
          II = IA(I) + I
C
C  1. H(IJ|IJ) <- 2*( FCOR(II) + FVAL(II) ), J CORE, I NOT CORE
C
          IF (I.GT.NCOR) THEN
            FAC = TWO*( FCOR(II) + FVAL(II) )
            DO J = 1, NCOR
              JIR = IROT(J,I)
              IF (JIR.NE.0) DIAH(JIR) = DIAH(JIR) + FAC
            END DO
          END IF
C
C  2. H(IJ|IJ) <- D(JJ)*FCOR(II), J VALENCE, I GENERAL
C
          DO J = NCP1, NOCC
            JIR = IROT(J,I)
            IF (JIR.NE.0) THEN
              JD = J - NCOR
              JJD = IA(JD) + JD
              DJJ = OPDM(JJD)
              DIAH(JIR) = DIAH(JIR) + DJJ*FCOR(II)
            END IF
          END DO
          IF (I.GT.NCOR.AND.I.LE.NOCC) THEN
            ID = I - NCOR
            DO J = NCP1, NOCC
              JIR = IROT(J,I)
              IF (JIR.NE.0) THEN
                JD = J - NCOR
                IJD = IA(MAX(ID,JD)) + MIN(ID,JD)
                DIJ = OPDM(IJD)
                IJ = IA(MAX(I,J)) + MIN(I,J)
                DIAH(JIR) = DIAH(JIR) - DIJ*FCOR(IJ)
              END IF
            END DO
          END IF
C
C  3. LAGRANGIAN TERM (UNDOCUMENTED?)
C
            DO J = 1, NOCC
              JIR = IROT(J,I)
              FAC = LAGN(J,J)
              IF (I.NE.J) FAC = -FAC
              IF (JIR.NE.0) DIAH(JIR) = DIAH(JIR) + FAC
            END DO
        END IF   !  PARALLEL
      END DO
C
C  2-EL TERMS OF THE AUGMENTED HESSIAN
C
C ------------ CONTRIBUTIONS FROM (OO|OO) INTEGRALS -----------
C
      CALL DDI_DISTRIB(D_OOOO,DDI_ME,ILO,IHI,JLO,JHI)
C
C  (ACTIVE-CORE|ACTIVE-CORE) TERMS
C
      DO I = NCP1, NOCC
        DO J = 1, NCOR
          IJ = IA(I) + J
          IF ((IJ.GE.JLO).AND.(IJ.LE.JHI)) THEN
            CALL DDI_GET(D_OOOO,1,NOTR,IJ,IJ,BUFF)
C
C  1. H(IJ|IJ) <- 6(IJ|IJ), J CORE
C
            IJR = IROT(I,J)
            IF (IJR.NE.0)
     *      DIAH(IJR) = DIAH(IJR) + SIX*BUFF(IJ)
C
C  2. H(IJ|IJ) <- -6*D(IK)*(IJ|JK), J CORE
C
            IJR = IROT(I,J)
            IF (IJR.NE.0) THEN
              ID = I - NCOR
              DO K = NCP1, NOCC
                JK = IA(K) + J
                ERI = BUFF(JK)
                KD = K - NCOR
                IKD = IA(MAX(ID,KD)) + MIN(ID,KD)
                DIK = OPDM(IKD)
                FAC = -DIK*ERI
                DIAH(IJR) = DIAH(IJR) + SIX*FAC
              END DO
            END IF
C
C  3. H(JK|JK) <- 2*D(IK|KL)*(IJ|JL), J CORE
C
            ID = I - NCOR
            DO K = NCP1, NOCC
              JKR = IROT(K,J)
              IF (JKR.NE.0) THEN
                KD = K - NCOR
                IKD = IA(MAX(ID,KD)) + MIN(ID,KD)
                DO L = NCP1, NOCC
                  JL = IA(L) + J
                  ERI = BUFF(JL)
                  LD = L - NCOR
                  KLD = IA(MAX(KD,LD)) + MIN(KD,LD)
                  IKKLD = IA(MAX(IKD,KLD)) + MIN(IKD,KLD)
                  DIKKL = TPDM(IKKLD)
                  FAC = DIKKL*ERI
                  DIAH(JKR) = DIAH(JKR) + TWO*FAC
                END DO
              END IF
            END DO
C
C  END OF DISTRIBUTION LOOPS
C
          END IF   ! LOCAL
        END DO  ! J
      END DO  ! I
C
C  (ACTIVE-ACTIVE|**) TERMS
C
      DO I = NCP1, NOCC
        DO J = NCP1, I
          IJ = IA(I) + J
          IF ((IJ.GE.JLO).AND.(IJ.LE.JHI)) THEN
            CALL DDI_GET(D_OOOO,1,NOTR,IJ,IJ,BUFF)
C
C  4. H(IK|IK) <- -2*(II|KK), K CORE
C
            IF (I.EQ.J) THEN
              DO K = 1, NCOR
                IKR = IROT(I,K)
                IF (IKR.NE.0) THEN
                  KK = IA(K) + K
                  DIAH(IKR) = DIAH(IKR) - TWO*BUFF(KK)
                END IF
              END DO
            END IF
C
C  5. H(IK|IK) <- 2*D(IJ)*(IJ|KK), K CORE
C
            INJ = I.NE.J
            ID = I - NCOR
            JD = J - NCOR
            IJD = IA(MAX(ID,JD)) + MIN(ID,JD)
            DO K = 1, NCOR
              KK = IA(K) + K
              ERI = BUFF(KK)*TWO
              DIJ = OPDM(IJD)
              FAC = DIJ*ERI
              IKR = IROT(I,K)
              IF (IKR.NE.0) DIAH(IKR) = DIAH(IKR) + FAC
              IF (INJ) THEN
                JKR = IROT(J,K)
                IF (JKR.NE.0) DIAH(JKR) = DIAH(JKR) + FAC
              END IF
            END DO
C
C  6. H(LK|LK) <- D(LL|IJ)*(KK|IJ), K CORE
C
            ID = I - NCOR
            JD = J - NCOR
            IJD = IA(ID) + JD
            INJ = I.NE.J
            DO K = 1, NCOR
              KK = IA(K) + K
              ERI = BUFF(KK)
              IF (INJ) ERI = ERI*TWO
              DO L = NCP1, NOCC
                KLR = IROT(L,K)
                IF (KLR.NE.0) THEN
                  LD = L - NCOR
                  LLD = IA(LD) + LD
                  IJLLD = IA(MAX(IJD,LLD)) + MIN(IJD,LLD)
                  DIJLL = TPDM(IJLLD)
                  DIAH(KLR) = DIAH(KLR) + DIJLL*ERI
                END IF
              END DO
            END DO
C
C  7. H(LK|LK) <- D(LL|IJ)*(KK|IJ)
C
            ID = I - NCOR
            JD = J - NCOR
            IJD = IA(ID) + JD
            INJ = I.NE.J
            DO K = NCP1, NOCC
              KK = IA(K) + K
              ERI = BUFF(KK)
              IF (INJ) ERI = ERI*TWO
              DO L = NCP1, NOCC
                KLR = IROT(K,L)
                IF (KLR.NE.0) THEN
                  LD = L - NCOR
                  LLD = IA(LD) + LD
                  IJLLD = IA(MAX(IJD,LLD)) + MIN(IJD,LLD)
                  DIJLL = TPDM(IJLLD)
                  FAC = DIJLL*ERI
                  DIAH(KLR) = DIAH(KLR) + FAC
                END IF
              END DO
            END DO
C
C  8. H(IK|IK) <- 2*D(JK|KL)*(IJ|IL)
C
            ID = I - NCOR
            JD = J - NCOR
            INJ = I.NE.J
            DO K = NCP1, NOCC
              IKR = IROT(I,K)
              IF (IKR.NE.0) THEN
                KD = K - NCOR
                JKD = IA(MAX(JD,KD)) + MIN(JD,KD)
                DO L = NCP1, NOCC
                  IL = IA(MAX(I,L)) + MIN(I,L)
                  ERI = BUFF(IL)
                  LD = L - NCOR
                  KLD = IA(MAX(KD,LD)) + MIN(KD,LD)
                  JKKLD = IA(MAX(JKD,KLD)) + MIN(JKD,KLD)
                  DJKKL = TPDM(JKKLD)
                  FAC = DJKKL*ERI
                  DIAH(IKR) = DIAH(IKR) + FAC*TWO
                END DO
              END IF
            END DO
            IF (INJ) THEN
              DO K = NCP1, NOCC
                JKR = IROT(J,K)
                IF (JKR.NE.0) THEN
                  KD = K - NCOR
                  IKD = IA(MAX(ID,KD)) + MIN(ID,KD)
                  DO L = NCP1, NOCC
                    JL = IA(MAX(J,L)) + MIN(J,L)
                    ERI = BUFF(JL)
                    LD = L - NCOR
                    KLD = IA(MAX(KD,LD)) + MIN(KD,LD)
                    IKKLD = IA(MAX(IKD,KLD)) + MIN(IKD,KLD)
                    DIKKL = TPDM(IKKLD)
                    FAC = DIKKL*ERI
                    DIAH(JKR) = DIAH(JKR) + FAC*TWO
                  END DO
                END IF
              END DO
            END IF
C
C  9. H(IJ|IJ) <- -2*D(KL|IJ)*(KL|IJ)
C
            IJR = IROT(I,J)
            IF (IJR.NE.0) THEN
              ID = I - NCOR
              JD = J - NCOR
              IJD = IA(ID) + JD
              DO K = NCP1, NOCC
                KD = K - NCOR
                DO L = NCP1, K
                  KL = IA(MAX(K,L)) + MIN(K,L)
                  ERI = BUFF(KL)
                  IF (K.NE.L) ERI = ERI*TWO
                  LD = L - NCOR
                  KLD = IA(KD) + LD
                  IJKLD = IA(MAX(IJD,KLD)) + MIN(IJD,KLD)
                  DIJKL = TPDM(IJKLD)
                  FAC = -DIJKL*ERI
                  DIAH(IJR) = DIAH(IJR) + FAC*TWO
                END DO
              END DO
            END IF
C
C  10. H(IK|IK) <- -4*D(IL|JK)*(KL|IJ)
C
            ID = I - NCOR
            JD = J - NCOR
            DO K = NCP1, I-1
              IKR = IROT(I,K)
              IF (IKR.NE.0) THEN
                KD = K - NCOR
                JKD = IA(MAX(JD,KD)) + MIN(JD,KD)
                DO L = NCP1, NOCC
                  KL = IA(MAX(K,L)) + MIN(K,L)
                  ERI = BUFF(KL)
                  LD = L - NCOR
                  ILD = IA(MAX(ID,LD)) + MIN(ID,LD)
                  ILJKD = IA(MAX(ILD,JKD)) + MIN(ILD,JKD)
                  DILJK = TPDM(ILJKD)
                  FAC = -DILJK*ERI
                  DIAH(IKR) = DIAH(IKR) + FAC*FOUR
                END DO
              END IF
            END DO
            IF (I.NE.J) THEN
              DO K = NCP1, J-1
                JKR = IROT(J,K)
                IF (JKR.NE.0) THEN
                  KD = K - NCOR
                  IKD = IA(MAX(ID,KD)) + MIN(ID,KD)
                  DO L = NCP1, NOCC
                    KL = IA(MAX(K,L)) + MIN(K,L)
                    ERI = BUFF(KL)
                    LD = L - NCOR
                    JLD = IA(MAX(JD,LD)) + MIN(JD,LD)
                    IKJLD = IA(MAX(IKD,JLD)) + MIN(IKD,JLD)
                    DIKJL = TPDM(IKJLD)
                    FAC = -DIKJL*ERI
                    DIAH(JKR) = DIAH(JKR) + FAC*FOUR
                  END DO
                END IF
              END DO
            END IF
C
C  END OF DISTRIBUTION LOOPS
C
          END IF   ! LOCAL
        END DO  ! J
      END DO  ! I
C
C ------------ CONTRIBUTIONS FROM (VV|OO) INTEGRALS -----------
C
      CALL DDI_DISTRIB(D_VVOO,DDI_ME,ILO,IHI,JLO,JHI)
C
C  (VIRTUAL-VIRTUAL|CORE-CORE) TERM
C
      DO I = 1, NCOR
        II = IA(I) + I
        IF ((II.GE.JLO).AND.(II.LE.JHI)) THEN
          CALL DDI_GET(D_VVOO,1,NVTR,II,II,BUFF)
C
C  11. H(AI|AI) <- -2*(AA|II), I CORE
C
          DO A = 1, NVIR
            AA = IA(A) + A
            AR = A + NOCC
            IAR = IROT(I,AR)
            IF (IAR.NE.0)
     *      DIAH(IAR) = DIAH(IAR) - TWO*BUFF(AA)
          END DO
        END IF   ! LOCAL
      END DO   !  I
C
C  (VIRTUAL-VIRTUAL|ACTIVE-ACTIVE) TERM
C
      DO I = NCP1, NOCC
        DO J = NCP1, I
          IJ = IA(I) + J
          IF ((IJ.GE.JLO).AND.(IJ.LE.JHI)) THEN
            CALL DDI_GET(D_VVOO,1,NVTR,IJ,IJ,BUFF)
C
C  12. H(AK|AK) <- D(KK|IJ)*(AA|IJ)
C
            INJ = I.NE.J
            ID = I - NCOR
            JD = J - NCOR
            IJD = IA(ID) + JD
            DO A = 1, NVIR
              AA = IA(A) + A
              AR = A + NOCC
              ERI = BUFF(AA)
              IF (INJ) ERI = ERI*TWO
              DO K = NCP1, NOCC
                KAR = IROT(K,AR)
                IF (KAR.NE.0) THEN
                  KD = K - NCOR
                  KKD = IA(KD) + KD
                  IJKKD = IA(MAX(IJD,KKD)) + MIN(IJD,KKD)
                  DIJKK = TPDM(IJKKD)
                  DIAH(KAR) = DIAH(KAR) + DIJKK*ERI
                END IF
              END DO
            END DO
C
C  END OF DISTRIBUTION LOOPS
C
          END IF   ! LOCAL
        END DO   !  J
      END DO   !  I
C
C ------------ CONTRIBUTIONS FROM (VO|VO) INTEGRALS -----------
C
      CALL DDI_DISTRIB(D_VOVO,DDI_ME,ILO,IHI,JLO,JHI)
C
C  (VIRTUAL-CORE|VIRTUAL-CORE) TERM
C
      DO I = 1, NCOR
        II = IA(I) + I
        IF ((II.GE.JLO).AND.(II.LE.JHI)) THEN
          CALL DDI_GET(D_VOVO,1,NVSQ,II,II,BUFF)
C
C  13. H(AI|AI) <- 6*(AI|AI), I CORE
C
          DO A = 1, NVIR
            AA = (A-1)*NVIR + A
            AR = A + NOCC
            IAR = IROT(I,AR)
            IF (IAR.NE.0)
     *      DIAH(IAR) = DIAH(IAR) + SIX*BUFF(AA)
          END DO
        END IF   ! LOCAL
      END DO   !  I
C
C  (VIRTUAL-ACTIVE|VIRTUAL-ACTIVE) TERM
C
      DO I = NCP1, NOCC
        DO J = NCP1, I
          IJ = IA(I) + J
          IF ((IJ.GE.JLO).AND.(IJ.LE.JHI)) THEN
            CALL DDI_GET(D_VOVO,1,NVSQ,IJ,IJ,BUFF)
C
C  14. H(AK|AK) <- 2*D(IK|JK)*(AI|AJ)
C
            INJ = I.NE.J
            ID = I - NCOR
            JD = J - NCOR
            DO A = 1, NVIR
              AA = (A-1)*NVIR + A
              AR = A + NOCC
              ERI = BUFF(AA)*TWO
              IF (INJ) ERI = ERI*TWO
              DO K = NCP1, NOCC
                KAR = IROT(K,AR)
                IF (KAR.NE.0) THEN
                  KD = K - NCOR
                  IKD = IA(MAX(ID,KD)) + MIN(ID,KD)
                  JKD = IA(MAX(JD,KD)) + MIN(JD,KD)
                  IKJKD = IA(MAX(IKD,JKD)) + MIN(IKD,JKD)
                  DIKJK = TPDM(IKJKD)
                  DIAH(KAR) = DIAH(KAR) + DIKJK*ERI
                END IF
              END DO
            END DO
C
C  END OF DISTRIBUTION LOOPS
C
          END IF   ! LOCAL
        END DO   !  J
      END DO   !  I
C
C  GLOBALLY SUM DIAGONAL ELEMENTS
C
      CALL DDI_GSUMF(2289,DIAH,NROT)
C
C  DOUBLE ELEMENTS OF DIAGONAL
C
      CALL DSCAL(NROT,TWO,DIAH,1)
C
C  B. LENGSFIELD'S LEVEL SHIFTING
C
      DO I = 1, NROT
        IF (DIAH(I).LE.ZERO) DIAH(I) = HSHIFT
      END DO
C
C  MOVE ELEMENTS FORWARD ONE PLACE
C
      DO I = NROT, 1, -1
        DIAH(I+1) = DIAH(I)
      END DO
C
C  FIRST ELEMENT IS ARBITRARY BUT SHOULD BE ZEROED
C
      DIAH(1) = ZERO
      RETURN
      END
C*MODULE MCTWO   *DECK NTNDVD
      SUBROUTINE NTNDVD(H,HD,V,W,T
     *,                 AA,A,VEC,EIG,WRK,IWRK, IA
     *,                 NWKS,MXXPAN,CVGTOL,MAXIT
     *,                 OUT,DBUG,PRTTOL
     *,                 NMOS,NCOR,NOCC,NROT
     *,                 OPDM,TPDM,FCOR,FVAL,LAGN,IROT,BUFF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C  DATA FOR PARALLEL CALL COMPUTING THE AUGMENTED HESSIAN ON THE FLY IS
C    NMOS,NCOR,NOCC,NROT,OPDM,TPDM,FCOR,FVAL,LAGN,IROT,BUFF
C  THESE ARGUMENTS ARE NOT USED BY THE SERIAL VERSION.
C
      LOGICAL CVGED,SOME,OUT,DBUG,PART1,GOPARR,DSKWRK,MASWRK
C
      COMMON /IOFILE/ IR,IW,IP,IJKO,IJKT,IDAF,NAV,IODA(950)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
C
      DIMENSION H(1),HD(1), V(NWKS,1),W(NWKS,1),T(NWKS)
      DIMENSION AA(1),A(1),VEC(MXXPAN,1),EIG(1),WRK(*),IWRK(*)
      DIMENSION IA(1), KCOEF(6),COEFK(6)
C
      INTEGER IROT(NOCC,*)
      DOUBLE PRECISION  OPDM(*),TPDM(*),FCOR(*),FVAL(*),LAGN(*),BUFF(*)
C
      DATA ZERO,ONE /0.0D+00,1.0D+00/
      DATA TOL      /1.0D-09/
C
C     ----- DAVIDSON'S METHOD MODIFIED TO DIAGONALIZE -----
C                 THE AUGMENTED HESSIAN MATRIX
C
C        ONLY ONE ROOT AND A THRESOLD VERY STRONG ( 1.0D-9 )
C
C        MINIMUM OF EXPANSIONS IS ONE ET MAXIMUM DEFAULT 50
C
C        MAXIMUM NUMBER OF ITERATIONS 200 (MAXDIA IN NUTON )
C
C     ----- GENERATE TRIAL EXPANSION VECTOR SET -----
C
      SOME =.FALSE.  .AND.  MASWRK
      THRES=TOL*TOL
C
      KSTAT = 1
      V(1   ,KSTAT)=ONE
      DO 10 IWKS=2,NWKS
   10 V(IWKS,KSTAT)=ZERO
C
      IF(DBUG) WRITE(IW,9996) KSTAT,(IWKS,V(IWKS,1),IWKS=1,NWKS)
C
      DO 20 I=1,MXXPAN
      DO 20 J=1,I
      IJ=IA(I)+J
   20 AA(IJ)=ZERO
C
C     ----- START ITERATION -----
C
      IF(SOME) WRITE(IW,9998)
      ITER=0
      NXPAN0=1
  100 CONTINUE
      ITER=ITER+1
      IF(ITER.GT.MAXIT) GO TO 9000
      IF(OUT) WRITE(IW,9997) ITER,NXPAN0
      NVEC=NXPAN0-1
      IF(NXPAN0.GT.MXXPAN) GO TO 700
  110 CONTINUE
C
C     ----- CALCULATE W = H * V FOR NEW EXPANSION VECTORS -----
C
      IF (GOPARR) THEN
C
C  PARALLEL W = H * V, WITHOUT STORING AUGMENTED HESSIAN (H)
C
        CALL AHPDDI(NMOS,NCOR,NOCC,NROT
     *,             OPDM,TPDM,FCOR,FVAL
     *,             LAGN,IROT,HD,BUFF
     *,             V(1,NXPAN0),W(1,NXPAN0))
      ELSE
        CALL DGHNTN(NWKS,V(1,NXPAN0),W(1,NXPAN0),H,IA)
      END IF
C
C     ----- CALCULATE NEW TRIANGULAR PART OF INTERACTION MATRIX -----
C
      DUMIJ=ZERO
      DO 120 IWKS=1,NWKS
  120 DUMIJ=DUMIJ+V(IWKS,NXPAN0)*W(IWKS,NXPAN0)
      IJ=IA(NXPAN0)+NXPAN0
      AA(IJ)=DUMIJ
C
C     ----- CALCULATE NEW BAND OF INTERACTION MATRIX -----
C
      IF(NVEC.EQ.0) GO TO 200
C
      DO 170 J=1,NVEC
      DUMIJ=ZERO
      DO 160 IWKS=1,NWKS
  160 DUMIJ=DUMIJ+V(IWKS,NXPAN0)*W(IWKS,J)
      IJ=IA(NXPAN0)+J
      AA(IJ)=DUMIJ
  170 CONTINUE
C
C     ----- SOLVE (NXPAN0*NXPAN0) EIGENVALUE PROBLEM -----
C
  200 CONTINUE
      NDUM=(NXPAN0*(NXPAN0+1))/2
      DO 210 IDUM=1,NDUM
  210 A(IDUM)=AA(IDUM)
      IF(DBUG) CALL PRTRI(A,NXPAN0)
C
C              PERHAPS SET KDIAG=3 FOR JACOBI?
C
      CALL GLDIAG(NXPAN0,NXPAN0,NXPAN0,A,WRK,EIG,VEC,IERR,IWRK)
C
      IF(OUT) WRITE(IW,9988) EIG(1)
C
C     ----- FORM CORRECTION VECTORS -----
C
      DO 310 IWKS=1,NWKS
  310 T(IWKS)=ZERO
      PART1=.TRUE.
  320 CONTINUE
C
      EIGK=EIG(1)
      DO 360 I=1,NXPAN0
      AIK=VEC(I,1)
      IF(PART1) AIK=-AIK*EIGK
      DO 350 IWKS=1,NWKS
  350 T(IWKS)=T(IWKS)+AIK*V(IWKS,I)
  360 CONTINUE
C
      DO 380 I=1,NXPAN0
      DO 380 IWKS=1,NWKS
      DUM=V(IWKS,I)
      V(IWKS,I)=W(IWKS,I)
  380 W(IWKS,I)=DUM
      IF(.NOT.PART1) GO TO 390
      PART1=.FALSE.
      GO TO 320
  390 CONTINUE
C
C     ----- CHECK CONVERGENCE -----
C
      CVG=ZERO
      CVGED=.TRUE.
      DUM=ZERO
      DO 410 IWKS=1,NWKS
  410 DUM=DUM+T(IWKS)*T(IWKS)
      DUM= SQRT(DUM)
      IF(OUT) WRITE(IW,9989) DUM
      IF(DUM.GT.CVG)  CVG=DUM
      CVGED=CVGED.AND.(DUM.LT.CVGTOL)
      IF(SOME) WRITE(IW,9995) ITER,CVG,EIG(1)
      NVEC=NXPAN0
      IF(CVGED) GO TO 700
C
      IF(ITER.EQ.1) GO TO 460
      EIGK=EIG(1)
      DUM=ZERO
      DO 430 IWKS=1,NWKS
      DENOM=EIGK-HD(IWKS)
      IF( ABS(DENOM).LT.TOL) DENOM=TOL
      T(IWKS)=T(IWKS)/DENOM
  430 DUM=DUM+T(IWKS)*T(IWKS)
      DUM=ONE/ SQRT(DUM)
      DO 440 IWKS=1,NWKS
  440 T(IWKS)=T(IWKS)*DUM
C
  460 CONTINUE
C
C     ----- ORTHOGONALIZE CORRECTION VECTORS AND EXPANSION VECTORS -----
C     ----- UPDATE SET OF EXPANSION VECTORS                        -----
C
      DO 550 I=1,NXPAN0
      DUMIK=ZERO
      DO 530 IWKS=1,NWKS
  530 DUMIK=DUMIK+T(IWKS)*V(IWKS,I)
      DO 540 IWKS=1,NWKS
  540 T(IWKS)=T(IWKS)-DUMIK*V(IWKS,I)
  550 CONTINUE
C
      IVEC = 0
      DUM=ZERO
      DO 610 IWKS=1,NWKS
  610 DUM=DUM+T(IWKS)*T(IWKS)
      IF(DUM.LT.THRES) GO TO 670
      DUM=ONE/ SQRT(DUM)
      DO 620 IWKS=1,NWKS
  620 T(IWKS)=T(IWKS)*DUM
C
      IVEC = 1
C
  670 CONTINUE
      IF(OUT) THEN
         WRITE(IW,9991) IVEC
         IF(IVEC.EQ.0) WRITE(IW,9990)
      END IF
C
C     ----- END OF CYCLE -----
C
      NXPAN0 = NXPAN0+1
      DO 680 IWKS=1,NWKS
  680 V(IWKS,NXPAN0)=T(IWKS)
      IF(DBUG) WRITE(IW,9996) NXPAN0,(IWKS,V(IWKS,NXPAN0),IWKS=1,NWKS)
      GO TO 100
C
C     ----- RE-ORTHONORMALIZE EXPANSION COEFFICIENT MATRIX -----
C
  700 CONTINUE
      DUMIJ=ZERO
      DO 710 K=1,NVEC
  710 DUMIJ=DUMIJ+VEC(K,1)*VEC(K,1)
      DUMIJ=ONE/ SQRT(DUMIJ)
      DO 740 K=1,NVEC
  740 VEC(K,1)=VEC(K,1)*DUMIJ
C
C     ----- GET APPROXIMATE OR CONVERGED -CI- VECTORS -----
C
      PART1=.TRUE.
  800 CONTINUE
      DO 810 IWKS=1,NWKS
  810 T(IWKS)=ZERO
C
      DO 830 IVEC=1,NVEC
      AIK=VEC(IVEC,1)
      DO 820 IWKS=1,NWKS
  820 T(IWKS)=T(IWKS)+AIK*V(IWKS,IVEC)
  830 CONTINUE
C
      IF(CVGED) GO TO 900
C
C     ----- IF .NOT.CVGED, USE VECTORS AS NEW EXPANSION VECTORS -----
C
      IF(.NOT.PART1) GO TO 870
      DO 850 IVEC=1,NVEC
      DO 850 IWKS=1,NWKS
  850 V(IWKS,IVEC)=W(IWKS,IVEC)
C
      IF(DBUG) WRITE(IW,9996) KSTAT,(IWKS,T(IWKS),IWKS=1,NWKS)
      DO 860 IWKS=1,NWKS
  860 W(IWKS,1)=T(IWKS)
      PART1=.FALSE.
      GO TO 800
C
  870 CONTINUE
      DO 880 IWKS=1,NWKS
      V(IWKS,1)=W(IWKS,1)
  880 W(IWKS,1)=T(IWKS)
C
      AA(1)=EIG(1)
      NXPAN0=1
      IVEC=1
      IF(OUT) WRITE(IW,9987)
      GO TO 110
C
C     ----- PRINT FINAL VECTORS ----
C
  900 CONTINUE
C     IF(SOME) WRITE(IW,9993)   EIG(1)
C9993 FORMAT(/,1X,'FIRST ROOT EIGENVALUE  = ',F18.9,/,1X,10(1H-),
C    1 8X,6(1H-))
      NCOEF=0
      DO 910 IWKS=1,NWKS
      V(IWKS,1)=T(IWKS)
      DUM=V(IWKS,1)
      IF( ABS(DUM).LT.PRTTOL) GO TO 910
      NCOEF=NCOEF+1
      KCOEF(NCOEF)=IWKS
      COEFK(NCOEF)=DUM
      IF(NCOEF.LT.6) GO TO 910
      IF(SOME) WRITE(IW,9992) (KCOEF(K),COEFK(K),K=1,NCOEF)
      NCOEF=0
  910 CONTINUE
      IF(NCOEF.EQ.0) GO TO 920
      IF(SOME) WRITE(IW,9992) (KCOEF(K),COEFK(K),K=1,NCOEF)
      NCOEF=0
  920 CONTINUE
C
      RETURN
 9000 CONTINUE
      IF(MASWRK) WRITE(IW,9994) CVG, CVGTOL
      RETURN
 9998 FORMAT(/,' ITER.    MAX.DEV.    STATE ENERGIES ',/,
     1       '        STATE NORM(D)                ')
 9997 FORMAT(' ITERATION ',I3,' WITH ',I3,' EXPANSION VECTORS.')
 9996 FORMAT(' EXPANSION VECTORS ',I3,/,(6(I5,F15.8)) )
 9995 FORMAT(I4  ,F12.8, F17.9)
 9994 FORMAT(' EXCESSIVE NUMBER OF ITERATIONS IN -NTNDVD-',
     1 ' DURING AUGMENTED HESSIAN MATRIX DIAGONALIZATION. STOP',/,
     2 ' CVG = ',E12.4,' CVGTOL = ',E12.4)
 9992 FORMAT(1X,6(I7,F14.7))
 9991 FORMAT(I3,' VECTOR(S) ADDED TO EXPANSION SET.')
 9990 FORMAT(' NO VECTORS ADDED TO THE EXPANSION SET IN -NTNDVD- ',/,
     1 ' EVEN THOUGH DIAGONALIZATION NOT CONVERGED TO REQUESTED',
     2 ' THRESHOLD.',/,' INCREASE CONVERGENCE THRESHOLD, IN',/,
     3 ' NAMELIST $NEWTON - STOP')
 9989 FORMAT(' CONVERGENCE CHECK FOR STATE  IS = ',F15.8)
 9988 FORMAT(' FIRST ROOT  EIGENVALUE = ',F17.9)
 9987 FORMAT(' .... EXPANSION VECTORS RESET .... ')
      END
C*MODULE MCTWO   *DECK DGHNTN
      SUBROUTINE DGHNTN(MAX,V,W,H,IA)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION V(1),W(1),H(1),IA(1)
C
C     ----- CALCULATE W = H * V -----
C
      DO 30 I=1,MAX
      IJ0=IA(I)
      W(I)=     H(IJ0+I)*V(I)
      DO 10 J=1,I-1
      W(I)=W(I)+H(IJ0+J)*V(J)
  10  CONTINUE
      DO 20 J=1,I-1
      W(J)=W(J)+H(IJ0+J)*V(I)
  20  CONTINUE
  30  CONTINUE
      RETURN
      END
C*MODULE MCTWO   *DECK NTNLEF
      SUBROUTINE NTNLEF(V,B,E,IROT,NROT,NORB,NORBS,OUT,SQCDF)
C
C     ----- B. LENGSFIELD'S EMPIRICAL FUDGING -----
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL OUT
      DIMENSION V(*),E(NORBS,*),B(*),IROT(NORBS,*)
      COMMON /IOFILE/ IR,IW,IP,IJKO,IJKT,IDAF,NAV,IODA(950)
      PARAMETER (ZERO=0.0D+00)
      PARAMETER (PT5=0.5D+00)
      PARAMETER (TEST=0.5D+00)
      PARAMETER (TEST1=0.1D+00)
C
      TESTV =TEST
      TESTV1=TEST1
      TESTV2= SQRT(TEST)
      C1=V(1)
      IF( ABS(C1).GT.TESTV2) GO TO 340
      C1SQ=C1**2
C
  300 CONTINUE
      C2SQ=ZERO
      DO 310 I=1,NROT
         C2=V(I+1)
         IF( ABS(C2).LT.TESTV1) GO TO 310
         C2SQ=C2SQ+C2**2
  310 CONTINUE
      IF(OUT) WRITE(IW,9994) C1SQ,C2SQ,TESTV1
      IF(C1SQ+C2SQ.GT.TESTV) GO TO 320
      TESTV1=TESTV1*PT5
      IF(OUT) WRITE(IW,9993)
      GO TO 300
C
  320 CONTINUE
      V(1)= SIGN(TESTV2,V(1))
      SCALE= SQRT((C1SQ+C2SQ-TESTV)/C2SQ)
      DO 330 I=1,NROT
         C2=V(I+1)
         IF( ABS(C2).LT.TESTV1) GO TO 330
         V(I+1)=C2*SCALE
  330 CONTINUE
C
  340 CONTINUE
      IF(V(1).GT.ZERO) GO TO 360
      DO 350 I=1,NROT+1
         V(I)=-V(I)
  350 CONTINUE
C
  360 CONTINUE
      SQCDF=ZERO
      DO 370 I=1,NROT
         DUM=V(I+1)
         SQCDF=SQCDF+DUM*DUM
         B(I)=DUM
  370 CONTINUE
      IF(OUT) WRITE(IW,9995) SQCDF
C
C     ----- CONSTRUCT ROTATION PARAMETERS IN -E- -----
C
      DO 420 IX=1,NORB
         DO 410 IG=1,NORBS
            E(IG,IX)=ZERO
            IB=IROT(IG,IX)
            IF(IB.EQ.0) GO TO 410
            VAL=B(IB)
            IF(IG.GT.IX) VAL=-VAL
            E(IG,IX)=VAL
  410    CONTINUE
  420 CONTINUE
C
      RETURN
 9993 FORMAT(' TESTV1 REDUCED BY FACTOR OF 2. ')
 9994 FORMAT(' VECTOR FIXUP USED. C1SQ = ',F10.8,' C2SQ = ',F10.8,
     1 ' TESTV1 = ',F7.3)
 9995 FORMAT(' SQCDF = ',E20.4)
      END
C*MODULE MCTWO   *DECK PNTNWT
      SUBROUTINE PNTNWT(NFT15,GIJ,LEN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION GIJ(LEN)
      LOGICAL GOPARR,DSKWRK,MASWRK
      COMMON /IOFILE/ IR,IW,IP,IJKO,IJKT,IDAF,NAV,IODA(950)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
C
      IF(MASWRK) READ(NFT15) GIJ
      CALL DAWRIT(IDAF,IODA,GIJ,LEN,68,0)
      RETURN
      END
