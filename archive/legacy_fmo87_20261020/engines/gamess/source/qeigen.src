C 21 Jun 15 - YN  - add QSANDL,QSANDTLUT for q.p. LUT-IOTC runs
C  7 Sep 12 - MWS - CHANGE EPSLON TO DEPSLON
C 22 NOV 11 - MWS - DQDIAG: USE 16 BIT ALIGNMENT FOR Q.P. DATA
C 14 OCT 09 - MWS - PAD GAUSS-HERMITE QUADRATURE COMMON
C  6 NOV 06 - MWS - ADJUST WAVEFUNCTION COMMON BLOCK
C 19 SEP 05 - MWS - ADD TRUE NUCLEAR CHARGE ARRAY TO INFOA COMMON
C  5 JUL 05 - MWS - SELECT NEW ATOM,BASIS,EFP,PCM,DAF DIMENSIONS
C 13 FEB 05 - MWS - PAD COMMON BLOCK HERMIT, WERMIT, NSHEL, ROOT
C  9 DEC 03 - MWS - SYNCH COMMON BLOCK RUNOPT
C  4 NOV 03 - DGF - ADD QUADRUPLE PRECISION DIAGONALIZATION, S+T INTS
C
C     THIS CODE IS USED FOR 128 BIT PRECISION DURING THE RESOLUTION
C     OF THE IDENTITY APPROXIMATION TO RELATIVISTIC PVP TRANSFORMATIONS.
C
C     IF THE REAL*16 DATA TYPE IS NOT AVAILABLE ON A COMPUTER,
C     THE COMPILING SCRIPT MUST CONVERT THIS SOURCE TO DOUBLE PRECISION,
C     AFTER WHICH THE INPUT OPTION $RELWFN MODEQR=8 WILL REFUSE TO RUN.
C
C*MODULE QEIGEN  *DECK QJACDG
      SUBROUTINE QJACDG(A,VEC,EIG,JBIG,BIG,LDVEC,N)
      USE prec, ONLY: qp
C
      IMPLICIT NONE
C
      INTEGER, INTENT(IN) :: LDVEC, N, JBIG(N)
      REAL(kind=qp) :: A(*),VEC(LDVEC,N),EIG(N),BIG(N)
C
      REAL(kind=qp), PARAMETER :: ONE = 1.0_qp
C
      INTEGER :: I, NB1, NB2, NMIN, NMAX
C
C     ----- JACOBI DIAGONALIZATION OF SYMMETRIC MATRIX -----
C     SYMMETRIC MATRIX -A- OF DIMENSION -N- IS DESTROYED ON EXIT.
C     ALL EIGENVECTORS ARE FOUND, SO -VEC- MUST BE SQUARE,
C     UNLESS SOMEONE TAKES THE TROUBLE TO LOOK AT -NMAX- BELOW.
C     -BIG- AND -JBIG- ARE SCRATCH WORK ARRAYS.
C
      CALL QCLR(VEC,1,LDVEC*N)
      DO I = 1,N
        VEC(I,I) = ONE
      END DO
C
      NB1 = N
      NB2 = (NB1*NB1+NB1)/2
      NMIN = 1
      NMAX = NB1
C
      CALL QJACDIA(A,VEC,NB1,NB2,LDVEC,NMIN,NMAX,BIG,JBIG)
C
      DO I=1,N
        EIG(I) = A((I*I+I)/2)
      END DO
C
      CALL QJACORD(VEC,EIG,NB1,LDVEC)
      RETURN
      END
C*MODULE EIGEN   *DECK QJACDIA
      SUBROUTINE QJACDIA(F,VEC,NB1,NB2,LDVEC,NMIN,NMAX,BIG,JBIG)
      USE prec, ONLY: qp
      IMPLICIT NONE
      INTEGER :: NB1,NB2,LDVEC,NMIN,NMAX,JBIG(NB1)
      REAL(kind=qp) :: F(NB2),VEC(LDVEC,NB1),BIG(NB1)
C
      LOGICAL GOPARR,DSKWRK,MASWRK
      INTEGER ME,MASTER,NPROC,IBTYP,IPTIM
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
C
C     PARAMETER (ROOT2=0.707106781186548Q+00 )
      REAL(kind=qp), PARAMETER :: ZERO=0.0_qp
      REAL(kind=qp), PARAMETER :: ONE=1.0_qp
      REAL(kind=qp), PARAMETER :: D1050=1.05_qp
      REAL(kind=qp), PARAMETER :: D1500=1.5_qp
      REAL(kind=qp), PARAMETER :: D3875=3.875_qp
      REAL(kind=qp), PARAMETER :: D0500=0.5_qp
      REAL(kind=qp), PARAMETER :: D1375=1.375_qp
      REAL(kind=qp), PARAMETER :: D0250=0.25_qp
      REAL(kind=qp), PARAMETER :: C2=1.0E-18_qp
      REAL(kind=qp), PARAMETER :: C3=4.0E-25_qp
      REAL(kind=qp), PARAMETER :: C4=2.0E-25_qp
      REAL(kind=qp), PARAMETER :: C5=8.0E-16_qp
      REAL(kind=qp), PARAMETER :: C6=3.0E-12_qp
      REAL(kind=qp), PARAMETER :: ROOT2=SQRT(D0500)
C
      INTEGER :: IEAA, IEAB, IEAR, IEBR
      INTEGER :: I, J, K, II, JJ, IR, IR1, ITER
      INTEGER :: I1, IT, IA, IAA, IB, IBB
      INTEGER :: KQ, MAXIT
      REAL(kind=qp) :: TT, T, T1, T2, T2X2, T2X25
      REAL(kind=qp) :: EPS, CX, SX, SD, TEST, DIF
      REAL(kind=qp) QEPSLON
C
C     NOTE THAT THE DOUBLE PRECISION CODE USES
C     PARAMETER (C2=1.0D-12, C3=4.0D-16,
C    *           C4=2.0D-16, C5=8.0D-09, C6=3.0D-06 )
C     WARNING: IT IS NOT CLEAR WHAT THESE NUMBERS SHOULD BE IN QUADRUPLE
C              PRECISION, AND IN PARTICULAR IF THE NUMBER OF SIGNIFICANT
C              FIGURES IN QUADRUPLE PRECISION IS MACHINE-INDEPENDENT.
C              ONE THING IS CLEAR:
C         *** DO NOT CALL THIS ROUTINE ON MACHINES WHERE 128-BIT
C         PRECISION IS NOT AVAILABLE!! GARBAGE WILL BE PRODUCED DUE TO
C         TRYING TO MAKE SENSE OUT OF ZEROES. THE THRESHOLDS ASSUME
C         AT LEAST 25 SIGNIFICANT FIGURES (AIX APPEARS TO HAVE 30).
C         25 MAY NOT BE ENOUGH FOR TRUE QUADRUPLE PRECISION MATRICES!
C         (HOPEFULLY ENOUGH FOR D.P. DIAGONALISED IN QUADRUPLE)
C
C      F IS THE MATRIX TO BE DIAGONALIZED, F IS STORED TRIANGULAR
C      VEC IS THE ARRAY OF EIGENVECTORS, DIMENSION NB1*NB1
C      BIG AND JBIG ARE TEMPORARY SCRATCH AREAS OF DIMENSION NB1
C      THE ROTATIONS AMONG THE FIRST NMIN BASIS FUNCTIONS ARE NOT
C      ACCOUNTED FOR.
C      THE ROTATIONS AMONG THE LAST NB1-NMAX BASIS FUNCTIONS ARE NOT
C      ACCOUNTED FOR.
C
      IEAA=0
      IEAB=0
      TT=ZERO
      EPS = 64.0_qp*QEPSLON(ONE)
C
C      LOOP OVER COLUMNS (K) OF TRIANGULAR MATRIX TO DETERMINE
C      LARGEST OFF-DIAGONAL ELEMENTS IN ROW(I).
C
      DO I=1,NB1
         BIG(I)=ZERO
         JBIG(I)=0
         IF(I.LT.NMIN  .OR.  I.EQ.1) CYCLE
         II = (I*I-I)/2
         J=MIN(I-1,NMAX)
         DO K=1,J
            IF(ABS(BIG(I)).GE.ABS(F(II+K))) CYCLE
            BIG(I)=F(II+K)
            JBIG(I)=K
         END DO
C        WRITE(6,*) I,'-TH LARGEST ',BIG(I)
      END DO
C
C     ----- 2X2 JACOBI ITERATIONS BEGIN HERE -----
C
      MAXIT=MAX(NB2*20,500)
      ITER=0
      DO
      ITER=ITER+1
C
C      FIND SMALLEST DIAGONAL ELEMENT
C
      SD=D1050
      JJ=0
      DO J=1,NB1
         JJ=JJ+J
         SD= MIN(SD,ABS(F(JJ)))
      END DO
      TEST = MAX(EPS, C2*MAX(SD,C6))
C
C      FIND LARGEST OFF-DIAGONAL ELEMENT
C
      T=ZERO
      I1=MAX(2,NMIN)
      IB = I1
      DO I=I1,NB1
         IF(T.GE.ABS(BIG(I))) CYCLE
         T= ABS(BIG(I))
         IB=I
      END DO
C
C      TEST FOR CONVERGENCE, THEN DETERMINE ROTATION.
C
      IF(T.LT.TEST) RETURN
C                   ******
C
      IF(ITER.GT.MAXIT) THEN
         IF (MASWRK) THEN
            WRITE(6,*) 'JACOBI DIAGONALIZATION FAILS, DIMENSION=',NB1
            WRITE(6,9020) ITER,T,TEST,SD
         ENDIF
         CALL ABRT
         STOP
      END IF
C
      IA=JBIG(IB)
      IAA=IA*(IA-1)/2
      IBB=IB*(IB-1)/2
      DIF=F(IAA+IA)-F(IBB+IB)
      IF(ABS(DIF).GT.C3*T) GO TO 70
      SX=ROOT2
      CX=ROOT2
      GO TO 110
   70 T2X2=BIG(IB)/DIF
      T2X25=T2X2*T2X2
      IF(T2X25 . GT . C4) GO TO 80
      CX=ONE
      SX=T2X2
      GO TO 110
   80 IF(T2X25 . GT . C5) GO TO 90
      SX=T2X2*(ONE-D1500*T2X25)
      CX=ONE-D0500*T2X25
      GO TO 110
   90 IF(T2X25 . GT . C6) GO TO 100
      CX=ONE+T2X25*(T2X25*D1375 - D0500)
      SX= T2X2*(ONE + T2X25*(T2X25*D3875 - D1500))
      GO TO 110
  100 T=D0250  / SQRT(D0250   + T2X25)
      CX= SQRT(D0500   + T)
      SX= SIGN( SQRT(D0500   - T),T2X2)
  110 IEAR=IAA+1
      IEBR=IBB+1
C
      DO IR=1,NB1
         T=F(IEAR)*SX
         F(IEAR)=F(IEAR)*CX+F(IEBR)*SX
         F(IEBR)=T-F(IEBR)*CX
         IF(IR-IA) 220,120,130
  120    TT=F(IEBR)
         IEAA=IEAR
         IEAB=IEBR
         F(IEBR)=BIG(IB)
         IEAR=IEAR+IR-1
         IF(JBIG(IR)) 200,220,200
  130    T=F(IEAR)
         IT=IA
         IEAR=IEAR+IR-1
         IF(IR-IB) 180,150,160
  150    F(IEAA)=F(IEAA)*CX+F(IEAB)*SX
         F(IEAB)=TT*CX+F(IEBR)*SX
         F(IEBR)=TT*SX-F(IEBR)*CX
         IEBR=IEBR+IR-1
         GO TO 200
  160    IF(  ABS(T) . GE .  ABS(F(IEBR))) GO TO 170
         IF(IB.GT.NMAX) GO TO 170
         T=F(IEBR)
         IT=IB
  170    IEBR=IEBR+IR-1
  180    IF(  ABS(T) . LT .  ABS(BIG(IR))) GO TO 190
         BIG(IR) = T
         JBIG(IR) = IT
         GO TO 220
  190    IF(IA . NE . JBIG(IR) . AND . IB . NE . JBIG(IR))  GO TO 220
  200    KQ=IEAR-IR-IA+1
         BIG(IR)=ZERO
         IR1=MIN(IR-1,NMAX)
         DO I=1,IR1
            K=KQ+I
            IF(ABS(BIG(IR)) . GE . ABS(F(K))) CYCLE
            BIG(IR) = F(K)
            JBIG(IR)=I
         END DO
  220    IEAR=IEAR+1
         IEBR=IEBR+1
      END DO
C
      DO I=1,NB1
         T1=VEC(I,IA)*CX + VEC(I,IB)*SX
         T2=VEC(I,IA)*SX - VEC(I,IB)*CX
         VEC(I,IA)=T1
         VEC(I,IB)=T2
      END DO
      END DO
C
 9020 FORMAT(1X,'ITER=',I6,' T,TEST,SD=',1P,3E20.10)
      END
C*MODULE EIGEN   *DECK QJACORD
      SUBROUTINE QJACORD(VEC,EIG,N,LDVEC)
      USE prec, ONLY: qp
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: N, LDVEC
      REAL(kind=qp) :: VEC(LDVEC,N),EIG(N)
C
      INTEGER :: I, J, JJ
      REAL(kind=qp) :: T
C
C     ---- SORT EIGENDATA INTO ASCENDING ORDER -----
C
      DO I = 1, N
         JJ = I
         DO J = I, N
            IF (EIG(J) .LT. EIG(JJ)) JJ = J
         END DO
         IF (JJ .EQ. I) CYCLE
         T = EIG(JJ)
         EIG(JJ) = EIG(I)
         EIG(I) = T
         DO J = 1, N
            T = VEC(J,JJ)
            VEC(J,JJ) = VEC(J,I)
            VEC(J,I) = T
         END DO
      END DO
      RETURN
      END
C*MODULE QEIGEN  *DECK QEPSLON
      FUNCTION QEPSLON (X)
      USE prec, ONLY: qp
      IMPLICIT NONE
      REAL(kind=qp) QEPSLON
C*
C*    AUTHORS -
C*       THIS ROUTINE WAS TAKEN FROM EISPACK EDITION 3 DATED 4/6/83
C*       THIS VERSION IS BY S. T. ELBERT, AMES LABORATORY-USDOE NOV 1986
C*
C*    PURPOSE -
C*       ESTIMATE UNIT ROUNDOFF IN QUANTITIES OF SIZE X.
C*
C*    ON ENTRY -
C*       X      - WORKING PRECISION REAL
C*                VALUES TO FIND QEPSLON FOR
C*
C*    ON EXIT -
C*       QEPSLON - WORKING PRECISION REAL
C*                SMALLEST POSITIVE VALUE SUCH THAT X+QEPSLON.NE.ZERO
C*
C*    QUALIFICATIONS -
C*       THIS ROUTINE SHOULD PERFORM PROPERLY ON ALL SYSTEMS
C*       SATISFYING THE FOLLOWING TWO ASSUMPTIONS,
C*          1.  THE BASE USED IN REPRESENTING FLOATING POINT
C*              NUMBERS IS NOT A POWER OF THREE.
C*          2.  THE QUANTITY  A  IN STATEMENT 10 IS REPRESENTED TO
C*              THE ACCURACY USED IN FLOATING POINT VARIABLES
C*              THAT ARE STORED IN MEMORY.
C*       THE STATEMENT NUMBER 10 AND THE GO TO 10 ARE INTENDED TO
C*       FORCE OPTIMIZING COMPILERS TO GENERATE CODE SATISFYING
C*       ASSUMPTION 2.
C*       UNDER THESE ASSUMPTIONS, IT SHOULD BE TRUE THAT,
C*              A  IS NOT EXACTLY EQUAL TO FOUR-THIRDS,
C*              B  HAS A ZERO FOR ITS LAST BIT OR DIGIT,
C*              C  IS NOT EXACTLY EQUAL TO ONE,
C*              EPS  MEASURES THE SEPARATION OF 1.0 FROM
C*                   THE NEXT LARGER FLOATING POINT NUMBER.
C*       THE DEVELOPERS OF EISPACK WOULD APPRECIATE BEING INFORMED
C*       ABOUT ANY SYSTEMS WHERE THESE ASSUMPTIONS DO NOT HOLD.
C*
C*    DIFFERENCES FROM EISPACK 3 -
C*       USE IS MADE OF PARAMETER STATEMENTS AND INTRINSIC FUNCTIONS
C*       --NO EXECUTEABLE CODE CHANGES--
C*
C*    NOTE -
C*       QUESTIONS AND COMMENTS CONCERNING EISPACK SHOULD BE DIRECTED TO
C*       B. S. GARBOW, APPLIED MATH. DIVISION, ARGONNE NATIONAL LAB.
C
      REAL(kind=qp) A,B,C,EPS,X,ABS
      REAL(kind=qp) ZERO, ONE, THREE, FOUR
C
      PARAMETER (ZERO=0.0_qp)
      PARAMETER (ONE=1.0_qp)
      PARAMETER (THREE=3.0_qp)
      PARAMETER (FOUR=4.0_qp)
C
C-----------------------------------------------------------------------
C
      A = FOUR/THREE
      DO
         B = A - ONE
         C = B + B + B
         EPS = ABS(C - ONE)
         IF (EPS .NE. ZERO) EXIT
      END DO
      QEPSLON = EPS*ABS(X)
      RETURN
      END
C*MODULE QEIGEN  *DECK QCLR
      SUBROUTINE QCLR(A,INCA,N)
      USE prec, ONLY: qp
C
      IMPLICIT NONE
C
      INTEGER, INTENT(IN) :: INCA, N
      REAL(kind=qp) :: A(*)
C
      REAL(kind=qp), PARAMETER :: ZERO=0.0_qp
C
      INTEGER :: L, LA
C
C     ----- ZERO OUT VECTOR -A-, USING INCREMENT -INCA- -----
C
      IF (INCA .NE. 1) GO TO 200
      DO L=1,N
         A(L) = ZERO
      END DO
      RETURN
C
  200 CONTINUE
      LA=1-INCA
      DO L=1,N
         LA=LA+INCA
         A(LA) = ZERO
      END DO
      RETURN
      END
C*MODULE QEIGEN  *DECK DQDIAG
      SUBROUTINE DQDIAG(A,VEC,EIG,LDVEC,N,IPREC)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NDQ=2)
C
      DIMENSION A(*),VEC(LDVEC,N),EIG(N)
      COMMON /FMCOM / X(1)
C
C     DIAGONALISE A DOUBLE PRECISION MATRIX WITHIN QUADRUPLE PRECISION
C     THE MATRIX A OF SIZE N IS STORED IN TRIANGULAR MODE.
C     IPREC=0: A IS STORED AS DOUBLE PRECISION
C     IPREC=1: A IS STORED AS QUADRUPLE PRECISION
C
      N1=(N*N+N)/2
C
C     INITIAL OFFSET ENSURES 16-BYTE ALIGNMENT OF ALL Q.P. ENTRIES!
C     NOTE: GAMESS' FAST MEMORY POOL IS ALLOCATED ON 8-BYTE BOUNDARY,
C     BUT THE OFFSET -LOADFM- FROM /FMCOM/ MUST SOMEHOW INVOLVE AN
C     ODD NUMBER OF WORDS.  EXPERIMENTATION WITH GFORTRAN SHOWED
C     THAT THE LOGIC BELOW GETS TO 16-BYTE ALIGNMENT FOR STORAGE
C     USED AS Q.P. INSIDE THE CALLS BELOW.
C
      CALL VALFM(LOADFM)
      IF(MOD(LOADFM,2).EQ.0) THEN
         IOFF=1
      ELSE
         IOFF=2
      ENDIF
      LQA   = LOADFM+ IOFF
      LQVEC = LQA   + NDQ*N1
      LQEIG = LQVEC + NDQ*N*N
      LBIG  = LQEIG + NDQ*N
      LJBIG = LBIG  + NDQ*N
      LAST  = LJBIG +     N
      NEED=LAST-LOADFM
      CALL GETFM(NEED)
C
C        JACOBI-TYPE DIAGONALIZATION IN Q.P. PRECISION
C
      IF(IPREC.EQ.0) THEN
        CALL DQCOPY(N1,A,1,X(LQA),1)
        CALL QJACDG(X(LQA),X(LQVEC),X(LQEIG),X(LJBIG),X(LBIG),N,N)
      ELSE
        CALL QJACDG(A,X(LQVEC),X(LQEIG),X(LJBIG),X(LBIG),N,N)
      ENDIF
C
C        COPY EIGENDATA TO D.P. STORAGE
C
      CALL QDCOPY(N,X(LQEIG),1,EIG,1)
      DO I=1,N
        CALL QDCOPY(N,X(LQVEC+(I-1)*N*NDQ),1,VEC(1,I),1)
      ENDDO
C
      CALL RETFM(NEED)
      RETURN
      END
C*MODULE QEIGEN  *DECK DQCOPY
      SUBROUTINE DQCOPY(N,DX,INCX,DY,INCY)
      USE prec, ONLY: qp
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION DX(*)
      REAL(kind=qp) DY(*)
C
C     CLONE OF DCOPY, TO COPY FROM DOUBLE TO QUADRUPLE PRECISION.
C           DY(I) <== DX(I)
C
      IF(N.LE.0) RETURN
      IF(INCX.EQ.1.AND.INCY.EQ.1) GO TO 20
C
C        CODE FOR UNEQUAL INCREMENTS OR EQUAL NON-UNIT INCREMENTS
C
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DY(IY) = DX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
   20 CONTINUE
      DO I=1,N
        DY(I) = DX(I)
      ENDDO
      RETURN
      END
C*MODULE QEIGEN  *DECK QDCOPY
      SUBROUTINE QDCOPY(N,DX,INCX,DY,INCY)
      USE prec, ONLY: qp
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION DY(*)
      REAL(kind=qp) DX(*)
C
C     CLONE OF DCOPY, TO COPY FROM QUADRUPLE TO DOUBLE PRECISION.
C           DY(I) <== DX(I)
C
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20
C
C        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
C          NOT EQUAL TO 1
C
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DY(IY) = DX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
   20 CONTINUE
      DO I=1,N
        DY(I) = DX(I)
      ENDDO
      RETURN
      END
C*MODULE QEIGEN  *DECK IQTEST
      SUBROUTINE IQTEST(SOME)
      USE prec, ONLY: qp, dp
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      REAL(kind=qp) QEPSLON,QE
      LOGICAL SOME,GOPARR,DSKWRK,MASWRK
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
C
      E =DEPSLON(1.0_dp)
      QE=QEPSLON(1.0_qp)
      IF(MASWRK  .AND.  SOME) WRITE(IW,9100) E,QE
      IF(dp.EQ.qp) THEN
         IF(MASWRK) WRITE(IW,9000)
         CALL ABRT
      ENDIF
      RETURN
 9000 FORMAT(/1X,'128 BIT PRECISION APPEARS TO BE UNAVAILABLE!')
 9100 FORMAT(1X,'64 BIT EPSILON:',E9.2,' AND 128 BIT EPSILON:',E9.2/)
      END
C*MODULE INT1    *DECK QSANDT
C>    
C> @author - Dmitri G. Fedorov
C>    
C> @date 8-2-23 George Schoendorff
C>  - Added support for h and i functions
C>
      SUBROUTINE QSANDT(SS,TT,S,T,LL2,SOME)
      USE prec, ONLY: qp
      use mx_limits, only: mxsh,mxgtot,mxatm,mxang,mxang2,mxang3,mxdim
      use constants, only: ix,iy,iz,jx,jy,jz
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION MOROKM
C
      LOGICAL SOME,IANDJ,NORM,DOUBLE,GOPARR,DSKWRK,MASWRK
      REAL(kind=qp) S(LL2),T(LL2),SBLK(MXDIM),TBLK(MXDIM),
     *        FT(MXDIM),DIJ(MXDIM),
     *        XIN(MXANG3),YIN(MXANG3),ZIN(MXANG3),
     *        XINT,YINT,ZINT,TAA,X0,Y0,Z0,XI,YI,ZI,XJ,YJ,ZJ,
     *        CSI,CPI,CDI,CFI,CGI,CHI,CII,CSJ,CPJ,CDJ,CFJ,CGJ,CHJ,CIJ,
     *        QTOL,RR,ARRI,AXI,AYI,AZI,AI,AJ,AA,AA1,DUM,DUM1,DUM2,T1,T2,
     *        YZ,AX,AY,AZ,FAC,ZERO,PT5,ONE,TWO,THREE,FIVE,SEVEN,NINE,
     *        ELEVEN,THIRTEEN,FIFTEEN,SQRT3,SQRT5,SQRT7,SQRT9,SQRT11,
     *        RLN10
C
      DIMENSION IJX(MXDIM),IJY(MXDIM),IJZ(MXDIM),SS(*),TT(*)
C
      PARAMETER (NDQ=2)
C
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON /NSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),
     *                CF(MXGTOT),CG(MXGTOT),CH(MXGTOT),CI(MXGTOT),
     *                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),
     *                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
      COMMON /OUTPUT/ NPRINT,ITOL,ICUT,NORMF,NORMP,NOPK
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /ROOT  / XX,U(13),W(13),NROOTS
      COMMON /RUNOPT/ RUNTYP,EXETYP,NEVALS,NGLEVL,NHLEVL
      COMMON /SYMIND/ TOL,II,JJ,LIT,LJT,MINI,MINJ,MAXI,MAXJ,IANDJ
      COMMON /WFNOPT/ SCFTYP,VBTYP,DFTYPE,TDDFTYP,CITYP,CCTYP,
     *                MPLEVL,MPCTYP
C
      PARAMETER (ZERO=0.0_qp)
      PARAMETER (PT5=0.5_qp)
      PARAMETER (ONE=1.0_qp)
      PARAMETER (TWO=2.0_qp)
      PARAMETER (THREE=3.0_qp)
      PARAMETER (FIVE=5.0_qp)
      PARAMETER (SEVEN=7.0_qp)
      PARAMETER (NINE=9.0_qp)
      PARAMETER (ELEVEN=11.0_qp)
      PARAMETER (THIRTEEN=13.0_qp)
      PARAMETER (FIFTEEN=15.0_qp)
C
      DATA MOROKM/8HMOROKUMA/,NONE/4HNONE/
C
C     ----- COMPUTE CONVENTIONAL S, AND T INTEGRALS -----
C             ----- IN QUADRUPLE PRECISION -----
C     HERMIT COEFFICIENTS ARE ROUNDED TO DOUBLE PRECISION, HOWEVER,
C     THE LOSS OF ACCURACY TO LARGE EXTENT COMES FROM MULTIPLYING
C     LARGE EXPONENTS IN THIS ROUTINE AND WE GET FULL ACCURACY HERE.
C
C     PI212=SQRT(ONE/ACOS(ZERO))
      SQRT3=SQRT(THREE)
      SQRT5=SQRT(FIVE)
      SQRT7=SQRT(SEVEN)
      SQRT9=THREE
      SQRT11=SQRT(ELEVEN)
      RLN10=LOG(NINE+ONE)
C     UNCON=RMETHOD.NE.ANONE.AND.MOD(MODQR,2).EQ.1
      QTOL = RLN10*ITOL
      NORM = NORMF .NE. 1 .OR. NORMP .NE. 1
C
C     ----- MOPAC INTEGRALS ARE DONE ELSEWHERE -----
C
      IF(MPCTYP.NE.NONE) THEN
C        CALL MPCINT
         CALL ABRT
         RETURN
      END IF
C
      IF (RUNTYP.EQ.MOROKM) THEN
         CALL ABRT
C        CALL STINT1(ISTART,IEND,JSTART,LOCIJ,NATST,NATED,ISAVE,L1,L2)
      ELSE
         ISTART = 1
         IEND   = NSHELL
         JSTART = 1
         LOCIJ  = 0
C        NATST  = 1
C        NATED  = NAT+NCHMAT
C        IF(ISEPS) NATED = NAT+NPS
C        ISAVE  = 0
C        L1 = NUM
C        IF(UNCON) L1=NUMU
C        L2 = (L1*(L1+1))/2
      END IF
C
      CALL QCLR(S,1,LL2)
      CALL QCLR(T,1,LL2)
C
C     ----- INTIALIZE PARALLEL -----
C
C     IPCOUNT = ME - 1
C
C     ----- I SHELL -----
C
      DO 720 II = ISTART,IEND
         I = KATOM(II)
         XI = C(1,I)
         YI = C(2,I)
         ZI = C(3,I)
         I1 = KSTART(II)
         I2 = I1+KNG(II)-1
         LIT = KTYPE(II)
         MINI = KMIN(II)
         MAXI = KMAX(II)
         LOCI = KLOC(II)-MINI-LOCIJ
C
C     ----- J SHELL -----
C
         DO 700 JJ = JSTART,II
C
C     ----- GO PARALLEL! -----
C     PENDING GLOBAL SUM FOR QUADRUPLE PRECISION
C
C           IF (GOPARR) THEN
C              IPCOUNT = IPCOUNT + 1
C              IF (MOD(IPCOUNT,NPROC).NE.0) GO TO 700
C           END IF
            J = KATOM(JJ)
            XJ = C(1,J)
            YJ = C(2,J)
            ZJ = C(3,J)
            J1 = KSTART(JJ)
            J2 = J1+KNG(JJ)-1
            LJT = KTYPE(JJ)
            MINJ = KMIN(JJ)
            MAXJ = KMAX(JJ)
            LOCJ = KLOC(JJ)-MINJ-LOCIJ
            NROOTS = (LIT+LJT-2)/2+1
            RR = (XI-XJ)**2+(YI-YJ)**2+(ZI-ZJ)**2
            IANDJ = II .EQ. JJ
C
C     ----- PREPARE INDICES FOR PAIRS OF (I,J) FUNCTIONS
C
            IJ = 0
            MAX = MAXJ
            DO 160 I = MINI,MAXI
               NX = IX(I)
               NY = IY(I)
               NZ = IZ(I)
               IF (IANDJ) MAX = I
               DO 140 J = MINJ,MAX
                  IJ = IJ+1
                  IJX(IJ) = NX+JX(J)
                  IJY(IJ) = NY+JY(J)
                  IJZ(IJ) = NZ+JZ(J)
                  IF (J.LE.1) FT(IJ) = THREE
                  IF ((J.GT.1).AND.(J.LE.4)) FT(IJ) = FIVE
                  IF ((J.GT.4).AND.(J.LE.10)) FT(IJ) = SEVEN
                  IF ((J.GT.10).AND.(J.LE.20)) FT(IJ) = NINE
                  IF ((J.GT.20).AND.(J.LE.35)) FT(IJ) = ELEVEN
                  IF ((J.GT.35).AND.(J.LE.56)) FT(IJ) = THIRTEEN
                  IF (J.GT.56) FT(IJ) = FIFTEEN
  140          CONTINUE
  160       CONTINUE
C
            CALL QCLR(SBLK,1,IJ)
            CALL QCLR(TBLK,1,IJ)
C
C     ----- I PRIMITIVE
C
            JGMAX = J2
            DO 520 IG = I1,I2
               AI = EX(IG)
               ARRI = AI*RR
               AXI = AI*XI
               AYI = AI*YI
               AZI = AI*ZI
               CSI = CS(IG)
               CPI = CP(IG)
               CDI = CD(IG)
               CFI = CF(IG)
               CGI = CG(IG)
               CHI = CH(IG)
               CII = CI(IG)
C
C     ----- J PRIMITIVE
C
               IF (IANDJ) JGMAX = IG
               DO 500 JG = J1,JGMAX
                  AJ = EX(JG)
                  AA = AI+AJ
                  AA1 = ONE/AA
                  DUM = AJ*ARRI*AA1
                  IF (DUM .GT. QTOL) GO TO 500
                  FAC = EXP(-DUM)
                  CSJ = CS(JG)
                  CPJ = CP(JG)
                  CDJ = CD(JG)
                  CFJ = CF(JG)
                  CGJ = CG(JG)
                  CHJ = CH(JG)
                  CIJ = CI(JG)
                  AX = (AXI+AJ*XJ)*AA1
                  AY = (AYI+AJ*YJ)*AA1
                  AZ = (AZI+AJ*ZJ)*AA1
C
C     ----- DENSITY FACTOR
C
                  DOUBLE=IANDJ.AND.IG.NE.JG
                  MAX = MAXJ
                  NN = 0
                  DUM1 = ZERO
                  DUM2 = ZERO
                  DO 220 I = MINI,MAXI
                     IF (I.EQ.1) DUM1=CSI*FAC
                     IF (I.EQ.2) DUM1=CPI*FAC
                     IF (I.EQ.5) DUM1=CDI*FAC
                     IF ((I.EQ. 8).AND.NORM) DUM1=DUM1*SQRT3
                     IF (I.EQ.11) DUM1=CFI*FAC
                     IF ((I.EQ.14).AND.NORM) DUM1=DUM1*SQRT5
                     IF ((I.EQ.20).AND.NORM) DUM1=DUM1*SQRT3
                     IF (I.EQ.21) DUM1=CGI*FAC
                     IF ((I.EQ.24).AND.NORM) DUM1=DUM1*SQRT7
                     IF ((I.EQ.30).AND.NORM) DUM1=DUM1*SQRT5/SQRT3
                     IF ((I.EQ.33).AND.NORM) DUM1=DUM1*SQRT3
                     IF ( I.EQ.36)           DUM1 = CHI*FAC
                     IF ((I.EQ.39).AND.NORM) DUM1 = DUM1*SQRT9
                     IF ((I.EQ.45).AND.NORM) DUM1 = DUM1*SQRT7/SQRT3
                     IF ((I.EQ.51).AND.NORM) DUM1 = DUM1*SQRT3
                     IF ((I.EQ.54).AND.NORM) DUM1 = DUM1*SQRT5/SQRT3
                     IF ( I.EQ.57)           DUM1 = CII*FAC
                     IF ((I.EQ.60).AND.NORM) DUM1 = DUM1*SQRT11
                     IF ((I.EQ.66).AND.NORM) DUM1 = DUM1*SQRT9/SQRT3
                     IF ((I.EQ.72).AND.NORM) DUM1 = DUM1*SQRT3
                     IF ((I.EQ.75).AND.NORM) THEN
                       DUM1 = DUM1*SQRT7/(SQRT3*SQRT5)
                     END IF
                     IF ((I.EQ.78).AND.NORM) DUM1 = DUM1*SQRT5
                     IF ((I.EQ.84).AND.NORM) DUM1 = DUM1*SQRT5/SQRT3
                     IF (IANDJ) MAX = I
                     DO 200 J = MINJ,MAX
                        IF (J.EQ.1) THEN
                           DUM2=DUM1*CSJ
                           IF (DOUBLE) THEN
                              IF (I.LE.1) THEN
                                 DUM2=DUM2+DUM2
                              ELSE
                                 DUM2=DUM2+CSI*CPJ*FAC
                              END IF
                           END IF
                        ELSE IF (J.EQ.2) THEN
                           DUM2=DUM1*CPJ
                           IF (DOUBLE) DUM2=DUM2+DUM2
                        ELSE IF (J.EQ.5) THEN
                           DUM2=DUM1*CDJ
                           IF (DOUBLE) DUM2=DUM2+DUM2
                        ELSE IF ((J.EQ.8).AND.NORM) THEN
                           DUM2=DUM2*SQRT3
                        ELSE IF (J.EQ.11) THEN
                           DUM2=DUM1*CFJ
                           IF (DOUBLE) DUM2=DUM2+DUM2
                        ELSE IF ((J.EQ.14).AND.NORM) THEN
                           DUM2=DUM2*SQRT5
                        ELSE IF ((J.EQ.20).AND.NORM) THEN
                           DUM2=DUM2*SQRT3
                        ELSE IF (J.EQ.21) THEN
                           DUM2=DUM1*CGJ
                           IF (DOUBLE) DUM2=DUM2+DUM2
                        ELSE IF ((J.EQ.24).AND.NORM) THEN
                           DUM2=DUM2*SQRT7
                        ELSE IF ((J.EQ.30).AND.NORM) THEN
                           DUM2=DUM2*SQRT5/SQRT3
                        ELSE IF ((J.EQ.33).AND.NORM) THEN
                           DUM2=DUM2*SQRT3
                        ELSE IF ( J.EQ.36) THEN
                           DUM2 = DUM1*CHJ
                           IF (DOUBLE) DUM2=DUM2+DUM2
                        ELSE IF((J.EQ.39).AND.NORM) THEN
                           DUM2 = DUM2*SQRT9
                        ELSE IF ((J.EQ.45).AND.NORM) THEN
                           DUM2 = DUM2*SQRT7/SQRT3
                        ELSE IF ((J.EQ.51).AND.NORM) THEN
                           DUM2 = DUM2*SQRT3
                        ELSE IF ((J.EQ.54).AND.NORM) THEN
                           DUM2 = DUM2*SQRT5/SQRT3
                        ELSE IF ( J.EQ.57) THEN
                           DUM2 = DUM1*CIJ
                           IF (DOUBLE) DUM2=DUM2+DUM2
                        ELSE IF ((J.EQ.60).AND.NORM) THEN
                           DUM2 = DUM2*SQRT11
                        ELSE IF ((J.EQ.66).AND.NORM) THEN
                           DUM2 = DUM2*SQRT9/SQRT3
                        ELSE IF ((J.EQ.72).AND.NORM) THEN
                           DUM2 = DUM2*SQRT3
                        ELSE IF ((J.EQ.75).AND.NORM) THEN
                           DUM2 = DUM2*SQRT7/(SQRT3*SQRT5)
                        ELSE IF ((J.EQ.78).AND.NORM) THEN
                           DUM2 = DUM2*SQRT5
                        ELSE IF ((J.EQ.84).AND.NORM) THEN
                           DUM2 = DUM2*SQRT5/SQRT3
                        END IF
                        NN = NN+1
                        DIJ(NN) = DUM2
  200                CONTINUE
  220             CONTINUE
C
C     ----- OVERLAP AND KINETIC ENERGY
C
                  TAA = SQRT(AA1)
                  T1 = -TWO*AJ*AJ*TAA
                  T2 = -PT5*TAA
                  X0 = AX
                  Y0 = AY
                  Z0 = AZ
                  IN = -MXANG
                  DO 320 I = 1,LIT
                     IN = IN+MXANG
                     NI = I
                     DO 300 J = 1,LJT
                        JN = IN+J
                        NJ = J
                        CALL QSTVINT(XINT,YINT,ZINT,
     *                          TAA,X0,Y0,Z0,XI,YI,ZI,XJ,YJ,ZJ,NI,NJ)
                        XIN(JN) = XINT*TAA
                        YIN(JN) = YINT*TAA
                        ZIN(JN) = ZINT*TAA
                        NJ = J+2
                        CALL QSTVINT(XINT,YINT,ZINT,
     *                          TAA,X0,Y0,Z0,XI,YI,ZI,XJ,YJ,ZJ,NI,NJ)
                        XIN(JN+MXANG2) = XINT*T1
                        YIN(JN+MXANG2) = YINT*T1
                        ZIN(JN+MXANG2) = ZINT*T1
                        NJ = J-2
                        IF (NJ .GT. 0) THEN
                           CALL QSTVINT(XINT,YINT,ZINT,
     *                          TAA,X0,Y0,Z0,XI,YI,ZI,XJ,YJ,ZJ,NI,NJ)
                        ELSE
                           XINT = ZERO
                           YINT = ZERO
                           ZINT = ZERO
                        END IF
                        N = (J-1)*(J-2)
                        DUM = N * T2
                        XIN(JN+MXANG2*2) = XINT*DUM
                        YIN(JN+MXANG2*2) = YINT*DUM
                        ZIN(JN+MXANG2*2) = ZINT*DUM
  300                CONTINUE
  320             CONTINUE
                  DO 340 I = 1,IJ
                     NX = IJX(I)
                     NY = IJY(I)
                     NZ = IJZ(I)
                     YZ = YIN(NY)*ZIN(NZ)
                     DUM = YZ*XIN(NX)
                     DUM1 = (XIN(NX+MXANG2)+XIN(NX+MXANG2*2))*YZ+
     *                      (YIN(NY+MXANG2)+YIN(NY+MXANG2*2))
     *                      *XIN(NX)*ZIN(NZ)+
     *                      (ZIN(NZ+MXANG2)+ZIN(NZ+MXANG2*2))
     *                      *XIN(NX)*YIN(NY)
                     SBLK(I) = SBLK(I) + DIJ(I)*DUM
                     TBLK(I) = TBLK(I) + DIJ(I)*(DUM*AJ*FT(I)+DUM1)
  340             CONTINUE
C
C     ----- END OF PRIMITIVE LOOPS -----
C
  500          CONTINUE
  520       CONTINUE
C
C     ----- COPY BLOCK OVERLAP AND KINETIC ENERGY MATRICES
C
            MAX = MAXJ
            NN = 0
            DO 620 I = MINI,MAXI
               LI = LOCI+I
               IN = (LI*(LI-1))/2
               IF (IANDJ) MAX = I
               DO 600 J = MINJ,MAX
                  LJ = LOCJ+J
                  JN = LJ+IN
                  NN = NN+1
                  S(JN) = SBLK(NN)
                  T(JN) = TBLK(NN)
  600          CONTINUE
  620       CONTINUE
C
C     ----- END OF SHELL LOOPS -----
C
  700    CONTINUE
  720 CONTINUE
C
C     ----- SUM UP PARTIAL CONTRIBUTIONS IF PARALLEL -----
C
      IF (GOPARR) THEN
C        CALL DDI_GSUMF(911,S,L2)
C        CALL DDI_GSUMF(912,T,L2)
      END IF
C
C     ----- SAVE S, AND T MATRICES ON THE DAF -----
C
      CALL DAWRIT(IDAF,IODA,S,LL2*NDQ,394,0)
      CALL DAWRIT(IDAF,IODA,T,LL2*NDQ,395,0)
C
C     REWRITE THE OVERLAP AND KINETIC INTEGRALS
C
      CALL QDCOPY(LL2,S,1,SS,1)
      CALL QDCOPY(LL2,T,1,TT,1)
      CALL DAWRIT(IDAF,IODA,SS,LL2,12,0)
      CALL DAWRIT(IDAF,IODA,TT,LL2,13,0)
C     IF (ISAVE.EQ.1) THEN
C        CALL DAWRIT(IDAF,IODA,H,LL2,311,0)
C        CALL DAWRIT(IDAF,IODA,S,LL2,312,0)
C     END IF
      IF(SOME) WRITE(6,*)
     *     'QUADRUPLE PRECISION RELATIVISTIC INTEGRALS WILL BE USED.'
      RETURN
      END
C
C*MODULE INT1    *DECK QSANDTLUT
C>
C>    @brief   Non-relativistic one-electron integral evaluation
C>             for the LUT scheme (quadruple precision)
C>
C>    @details This is a modified version of QSANDT to obtain
C>             one-electron integrals for atomic transformation.
C>
C>    @author  Yuya Nakajima, June, 2015
C>
C>    @date 8-2-23 George Schoendorff
C>     - Added support for h and i functions
C>                      
C>    @param   SS    d.p. storage for overlap integrals
C>    @param   TT    d.p. storage for kinetic energy integrals
C>    @param   S     q.p. storage for overlap integrals
C>    @param   T     q.p. storage for kinetic energy integrals
C>    @param   LL2   size of the four integral arrays
C>    @param   SOME  set .TRUE. to see debug printing
C>
      SUBROUTINE QSANDTLUT(SS,TT,S,T,LL2,SOME)
      USE prec, ONLY: qp
      use mx_limits, only: mxsh,mxgtot,mxatm,mxang,mxang2,mxang3,mxdim
      use constants, only: ix,iy,iz,jx,jy,jz
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION MOROKM
C
      LOGICAL SOME,IANDJ,NORM,DOUBLE,GOPARR,DSKWRK,MASWRK
      REAL(kind=qp) S(LL2),T(LL2),SBLK(MXDIM),TBLK(MXDIM),
     *        FT(MXDIM),DIJ(MXDIM),
     *        XIN(MXANG3),YIN(MXANG3),ZIN(MXANG3),
     *        XINT,YINT,ZINT,TAA,X0,Y0,Z0,XI,YI,ZI,XJ,YJ,ZJ,
     *        CSI,CPI,CDI,CFI,CGI,CHI,CII,CSJ,CPJ,CDJ,CFJ,CGJ,CHJ,CIJ,
     *        QTOL,RR,ARRI,AXI,AYI,AZI,AI,AJ,AA,AA1,DUM,DUM1,DUM2,T1,T2,
     *        YZ,AX,AY,AZ,FAC,ZERO,PT5,ONE,TWO,THREE,FIVE,SEVEN,NINE,
     *        ELEVEN,THIRTEEN,FIFTEEN,SQRT3,SQRT5,SQRT7,SQRT9,SQRT11,
     *        RLN10
C
      DIMENSION IJX(MXDIM),IJY(MXDIM),IJZ(MXDIM),SS(*),TT(*)
C
C
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON /NSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),
     *                CF(MXGTOT),CG(MXGTOT),CH(MXGTOT),CI(MXGTOT),
     *                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),
     *                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
      COMMON /OUTPUT/ NPRINT,ITOL,ICUT,NORMF,NORMP,NOPK
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /ROOT  / XX,U(13),W(13),NROOTS
      COMMON /RUNOPT/ RUNTYP,EXETYP,NEVALS,NGLEVL,NHLEVL
      COMMON /SYMIND/ TOL,II,JJ,LIT,LJT,MINI,MINJ,MAXI,MAXJ,IANDJ
      COMMON /WFNOPT/ SCFTYP,VBTYP,DFTYPE,TDDFTYP,CITYP,CCTYP,
     *                MPLEVL,MPCTYP
C
      PARAMETER (ZERO=0.0_qp)
      PARAMETER (PT5=0.5_qp)
      PARAMETER (ONE=1.0_qp)
      PARAMETER (TWO=2.0_qp)
      PARAMETER (THREE=3.0_qp)
      PARAMETER (FIVE=5.0_qp)
      PARAMETER (SEVEN=7.0_qp)
      PARAMETER (NINE=9.0_qp)
      PARAMETER (ELEVEN=11.0_qp)
      PARAMETER (THIRTEEN=13.0_qp)
      PARAMETER (FIFTEEN=15.0_qp)
C
      DATA MOROKM/8HMOROKUMA/,NONE/4HNONE/
C
C     ----- COMPUTE CONVENTIONAL S, AND T INTEGRALS -----
C             ----- IN QUADRUPLE PRECISION -----
C     HERMIT COEFFICIENTS ARE ROUNDED TO DOUBLE PRECISION, HOWEVER,
C     THE LOSS OF ACCURACY TO LARGE EXTENT COMES FROM MULTIPLYING
C     LARGE EXPONENTS IN THIS ROUTINE AND WE GET FULL ACCURACY HERE.
C
C     PI212=SQRT(ONE/ACOS(ZERO))
      SQRT3=SQRT(THREE)
      SQRT5=SQRT(FIVE)
      SQRT7=SQRT(SEVEN)
      SQRT9=THREE
      SQRT11=SQRT(ELEVEN)
      RLN10=LOG(NINE+ONE)
C     UNCON=RMETHOD.NE.ANONE.AND.MOD(MODQR,2).EQ.1
      QTOL = RLN10*ITOL
      NORM = NORMF .NE. 1 .OR. NORMP .NE. 1
C
C     ----- MOPAC INTEGRALS ARE DONE ELSEWHERE -----
C
      IF(MPCTYP.NE.NONE) THEN
C        CALL MPCINT
         CALL ABRT
         RETURN
      END IF
C
      IF (RUNTYP.EQ.MOROKM) THEN
         CALL ABRT
C        CALL STINT1(ISTART,IEND,JSTART,LOCIJ,NATST,NATED,ISAVE,L1,L2)
      ELSE
         ISTART = 1
         IEND   = NSHELL
         JSTART = 1
         LOCIJ  = 0
C        NATST  = 1
C        NATED  = NAT+NCHMAT
C        IF(ISEPS) NATED = NAT+NPS
C        ISAVE  = 0
C        L1 = NUM
C        IF(UNCON) L1=NUMU
C        L2 = (L1*(L1+1))/2
      END IF
C
      CALL QCLR(S,1,LL2)
      CALL QCLR(T,1,LL2)
C
C     ----- INTIALIZE PARALLEL -----
C
C     IPCOUNT = ME - 1
C
C     ----- I SHELL -----
C
      DO 720 II = ISTART,IEND
         I = KATOM(II)
         XI = C(1,I)
         YI = C(2,I)
         ZI = C(3,I)
         I1 = KSTART(II)
         I2 = I1+KNG(II)-1
         LIT = KTYPE(II)
         MINI = KMIN(II)
         MAXI = KMAX(II)
         LOCI = KLOC(II)-MINI-LOCIJ
C
C     ----- J SHELL -----
C
         DO 700 JJ = JSTART,II
C
C     ----- GO PARALLEL! -----
C     PENDING GLOBAL SUM FOR QUADRUPLE PRECISION
C
C           IF (GOPARR) THEN
C              IPCOUNT = IPCOUNT + 1
C              IF (MOD(IPCOUNT,NPROC).NE.0) GO TO 700
C           END IF
            J = KATOM(JJ)
            XJ = C(1,J)
            YJ = C(2,J)
            ZJ = C(3,J)
            J1 = KSTART(JJ)
            J2 = J1+KNG(JJ)-1
            LJT = KTYPE(JJ)
            MINJ = KMIN(JJ)
            MAXJ = KMAX(JJ)
            LOCJ = KLOC(JJ)-MINJ-LOCIJ
            NROOTS = (LIT+LJT-2)/2+1
            RR = (XI-XJ)**2+(YI-YJ)**2+(ZI-ZJ)**2
            IANDJ = II .EQ. JJ
C
C     ----- PREPARE INDICES FOR PAIRS OF (I,J) FUNCTIONS
C
            IJ = 0
            MAX = MAXJ
            DO 160 I = MINI,MAXI
               NX = IX(I)
               NY = IY(I)
               NZ = IZ(I)
               IF (IANDJ) MAX = I
               DO 140 J = MINJ,MAX
                  IJ = IJ+1
                  IJX(IJ) = NX+JX(J)
                  IJY(IJ) = NY+JY(J)
                  IJZ(IJ) = NZ+JZ(J)
                  IF (J.LE.1) FT(IJ) = THREE
                  IF ((J.GT.1).AND.(J.LE.4)) FT(IJ) = FIVE
                  IF ((J.GT.4).AND.(J.LE.10)) FT(IJ) = SEVEN
                  IF ((J.GT.10).AND.(J.LE.20)) FT(IJ) = NINE
                  IF ((J.GT.20).AND.(J.LE.35)) FT(IJ) = ELEVEN
                  IF ((J.GT.35).AND.(J.LE.56)) FT(IJ) = THIRTEEN
                  IF (J.GT.56) FT(IJ) = FIFTEEN
  140          CONTINUE
  160       CONTINUE
C
            CALL QCLR(SBLK,1,IJ)
            CALL QCLR(TBLK,1,IJ)
C
C     ----- I PRIMITIVE
C
            JGMAX = J2
            DO 520 IG = I1,I2
               AI = EX(IG)
               ARRI = AI*RR
               AXI = AI*XI
               AYI = AI*YI
               AZI = AI*ZI
               CSI = CS(IG)
               CPI = CP(IG)
               CDI = CD(IG)
               CFI = CF(IG)
               CGI = CG(IG)
               CHI = CH(IG)
               CII = CI(IG)
C
C     ----- J PRIMITIVE
C
               IF (IANDJ) JGMAX = IG
               DO 500 JG = J1,JGMAX
                  AJ = EX(JG)
                  AA = AI+AJ
                  AA1 = ONE/AA
                  DUM = AJ*ARRI*AA1
                  IF (DUM .GT. QTOL) GO TO 500
                  FAC = EXP(-DUM)
                  CSJ = CS(JG)
                  CPJ = CP(JG)
                  CDJ = CD(JG)
                  CFJ = CF(JG)
                  CGJ = CG(JG)
                  CHJ = CH(JG)
                  CIJ = CI(JG)
                  AX = (AXI+AJ*XJ)*AA1
                  AY = (AYI+AJ*YJ)*AA1
                  AZ = (AZI+AJ*ZJ)*AA1
C
C     ----- DENSITY FACTOR
C
                  DOUBLE=IANDJ.AND.IG.NE.JG
                  MAX = MAXJ
                  NN = 0
                  DUM1 = ZERO
                  DUM2 = ZERO
                  DO 220 I = MINI,MAXI
                     IF (I.EQ.1) DUM1=CSI*FAC
                     IF (I.EQ.2) DUM1=CPI*FAC
                     IF (I.EQ.5) DUM1=CDI*FAC
                     IF ((I.EQ. 8).AND.NORM) DUM1=DUM1*SQRT3
                     IF (I.EQ.11) DUM1=CFI*FAC
                     IF ((I.EQ.14).AND.NORM) DUM1=DUM1*SQRT5
                     IF ((I.EQ.20).AND.NORM) DUM1=DUM1*SQRT3
                     IF (I.EQ.21) DUM1=CGI*FAC
                     IF ((I.EQ.24).AND.NORM) DUM1=DUM1*SQRT7
                     IF ((I.EQ.30).AND.NORM) DUM1=DUM1*SQRT5/SQRT3
                     IF ((I.EQ.33).AND.NORM) DUM1=DUM1*SQRT3
                     IF ( I.EQ.36)           DUM1 = CHI*FAC
                     IF ((I.EQ.39).AND.NORM) DUM1 = DUM1*SQRT9
                     IF ((I.EQ.45).AND.NORM) DUM1 = DUM1*SQRT7/SQRT3
                     IF ((I.EQ.51).AND.NORM) DUM1 = DUM1*SQRT3
                     IF ((I.EQ.54).AND.NORM) DUM1 = DUM1*SQRT5/SQRT3
                     IF ( I.EQ.57)           DUM1 = CII*FAC
                     IF ((I.EQ.60).AND.NORM) DUM1 = DUM1*SQRT11
                     IF ((I.EQ.66).AND.NORM) DUM1 = DUM1*SQRT9/SQRT3
                     IF ((I.EQ.72).AND.NORM) DUM1 = DUM1*SQRT3
                     IF ((I.EQ.75).AND.NORM) THEN
                       DUM1 = DUM1*SQRT7/(SQRT3*SQRT5)
                     END IF
                     IF ((I.EQ.78).AND.NORM) DUM1 = DUM1*SQRT5
                     IF ((I.EQ.84).AND.NORM) DUM1 = DUM1*SQRT5/SQRT3
                     IF (IANDJ) MAX = I
                     DO 200 J = MINJ,MAX
                        IF (J.EQ.1) THEN
                           DUM2=DUM1*CSJ
                           IF (DOUBLE) THEN
                              IF (I.LE.1) THEN
                                 DUM2=DUM2+DUM2
                              ELSE
                                 DUM2=DUM2+CSI*CPJ*FAC
                              END IF
                           END IF
                        ELSE IF (J.EQ.2) THEN
                           DUM2=DUM1*CPJ
                           IF (DOUBLE) DUM2=DUM2+DUM2
                        ELSE IF (J.EQ.5) THEN
                           DUM2=DUM1*CDJ
                           IF (DOUBLE) DUM2=DUM2+DUM2
                        ELSE IF ((J.EQ.8).AND.NORM) THEN
                           DUM2=DUM2*SQRT3
                        ELSE IF (J.EQ.11) THEN
                           DUM2=DUM1*CFJ
                           IF (DOUBLE) DUM2=DUM2+DUM2
                        ELSE IF ((J.EQ.14).AND.NORM) THEN
                           DUM2=DUM2*SQRT5
                        ELSE IF ((J.EQ.20).AND.NORM) THEN
                           DUM2=DUM2*SQRT3
                        ELSE IF (J.EQ.21) THEN
                           DUM2=DUM1*CGJ
                           IF (DOUBLE) DUM2=DUM2+DUM2
                        ELSE IF ((J.EQ.24).AND.NORM) THEN
                           DUM2=DUM2*SQRT7
                        ELSE IF ((J.EQ.30).AND.NORM) THEN
                           DUM2=DUM2*SQRT5/SQRT3
                        ELSE IF ((J.EQ.33).AND.NORM) THEN
                           DUM2=DUM2*SQRT3
                        ELSE IF ( J.EQ.36) THEN
                           DUM2 = DUM1*CHJ
                           IF (DOUBLE) DUM2=DUM2+DUM2
                        ELSE IF ((J.EQ.39).AND.NORM) THEN
                           DUM2 = DUM2*SQRT9
                        ELSE IF ((J.EQ.45).AND.NORM) THEN
                           DUM2 = DUM2*SQRT7/SQRT3
                        ELSE IF ((J.EQ.51).AND.NORM) THEN
                           DUM2 = DUM2*SQRT3
                        ELSE IF ((J.EQ.54).AND.NORM) THEN
                           DUM2 = DUM2*SQRT5/SQRT3
                        ELSE IF ( J.EQ.57) THEN
                           DUM2 = DUM1*CIJ
                           IF (DOUBLE) DUM2=DUM2+DUM2
                        ELSE IF ((J.EQ.60).AND.NORM) THEN
                           DUM2 = DUM2*SQRT11
                        ELSE IF ((J.EQ.66).AND.NORM) THEN
                           DUM2 = DUM2*SQRT9/SQRT3
                        ELSE IF ((J.EQ.72).AND.NORM) THEN
                           DUM2 = DUM2*SQRT3
                        ELSE IF ((J.EQ.75).AND.NORM) THEN
                           DUM2 = DUM2*SQRT7/(SQRT3*SQRT5)
                        ELSE IF ((J.EQ.78).AND.NORM) THEN
                           DUM2 = DUM2*SQRT5
                        ELSE IF ((J.EQ.84).AND.NORM) THEN
                           DUM2 = DUM2*SQRT5/SQRT3
                        END IF
                        NN = NN+1
                        DIJ(NN) = DUM2
  200                CONTINUE
  220             CONTINUE
C
C     ----- OVERLAP AND KINETIC ENERGY
C
                  TAA = SQRT(AA1)
                  T1 = -TWO*AJ*AJ*TAA
                  T2 = -PT5*TAA
                  X0 = AX
                  Y0 = AY
                  Z0 = AZ
                  IN = -MXANG
                  DO 320 I = 1,LIT
                     IN = IN+MXANG
                     NI = I
                     DO 300 J = 1,LJT
                        JN = IN+J
                        NJ = J
                        CALL QSTVINT(XINT,YINT,ZINT,
     *                          TAA,X0,Y0,Z0,XI,YI,ZI,XJ,YJ,ZJ,NI,NJ)
                        XIN(JN) = XINT*TAA
                        YIN(JN) = YINT*TAA
                        ZIN(JN) = ZINT*TAA
                        NJ = J+2
                        CALL QSTVINT(XINT,YINT,ZINT,
     *                          TAA,X0,Y0,Z0,XI,YI,ZI,XJ,YJ,ZJ,NI,NJ)
                        XIN(JN+MXANG2) = XINT*T1
                        YIN(JN+MXANG2) = YINT*T1
                        ZIN(JN+MXANG2) = ZINT*T1
                        NJ = J-2
                        IF (NJ .GT. 0) THEN
                           CALL QSTVINT(XINT,YINT,ZINT,
     *                          TAA,X0,Y0,Z0,XI,YI,ZI,XJ,YJ,ZJ,NI,NJ)
                        ELSE
                           XINT = ZERO
                           YINT = ZERO
                           ZINT = ZERO
                        END IF
                        N = (J-1)*(J-2)
                        DUM = N * T2
                        XIN(JN+MXANG2*2) = XINT*DUM
                        YIN(JN+MXANG2*2) = YINT*DUM
                        ZIN(JN+MXANG2*2) = ZINT*DUM
  300                CONTINUE
  320             CONTINUE
                  DO 340 I = 1,IJ
                     NX = IJX(I)
                     NY = IJY(I)
                     NZ = IJZ(I)
                     YZ = YIN(NY)*ZIN(NZ)
                     DUM = YZ*XIN(NX)
                     DUM1 = (XIN(NX+MXANG2)+XIN(NX+MXANG2*2))*YZ+
     *                      (YIN(NY+MXANG2)+YIN(NY+MXANG2*2))
     *                      *XIN(NX)*ZIN(NZ)+
     *                      (ZIN(NZ+MXANG2)+ZIN(NZ+MXANG2*2))
     *                      *XIN(NX)*YIN(NY)
                     SBLK(I) = SBLK(I) + DIJ(I)*DUM
                     TBLK(I) = TBLK(I) + DIJ(I)*(DUM*AJ*FT(I)+DUM1)
  340             CONTINUE
C
C     ----- END OF PRIMITIVE LOOPS -----
C
  500          CONTINUE
  520       CONTINUE
C
C     ----- COPY BLOCK OVERLAP AND KINETIC ENERGY MATRICES
C
            MAX = MAXJ
            NN = 0
            DO 620 I = MINI,MAXI
               LI = LOCI+I
               IN = (LI*(LI-1))/2
               IF (IANDJ) MAX = I
               DO 600 J = MINJ,MAX
                  LJ = LOCJ+J
                  JN = LJ+IN
                  NN = NN+1
                  S(JN) = SBLK(NN)
                  T(JN) = TBLK(NN)
  600          CONTINUE
  620       CONTINUE
C
C     ----- END OF SHELL LOOPS -----
C
  700    CONTINUE
  720 CONTINUE
C
C     ----- SUM UP PARTIAL CONTRIBUTIONS IF PARALLEL -----
C
      IF (GOPARR) THEN
C        CALL DDI_GSUMF(911,S,L2)
C        CALL DDI_GSUMF(912,T,L2)
      END IF
C
C     ----- SAVE S, AND T MATRICES ON THE DAF -----
C
C     NDQ=2
C     CALL DAWRIT(IDAF,IODA,S,LL2*NDQ,394,0)
C     CALL DAWRIT(IDAF,IODA,T,LL2*NDQ,395,0)
C
C     REWRITE THE OVERLAP AND KINETIC INTEGRALS
C
      CALL QDCOPY(LL2,S,1,SS,1)
      CALL QDCOPY(LL2,T,1,TT,1)
      CALL DAWRIT(IDAF,IODA,SS,LL2,12,0)
      CALL DAWRIT(IDAF,IODA,TT,LL2,13,0)
C     IF (ISAVE.EQ.1) THEN
C        CALL DAWRIT(IDAF,IODA,H,LL2,311,0)
C        CALL DAWRIT(IDAF,IODA,S,LL2,312,0)
C     END IF
      IF(SOME) WRITE(6,*)
     *     'QUADRUPLE PRECISION RELATIVISTIC INTEGRALS WILL BE USED.'
      RETURN
      END
C
C*MODULE INT1    *DECK QSANDTL
C>
C>    @brief   Non-relativistic one-electron integral evaluation
C>             for atomic IOTC transformation (quadruple precision)
C>
C>    @details Thie is a modified version of HSANDT to obtain
C>             one-electron integrals foratomic transformation.
C>
C>    @author  Yuya Nakajima, June, 2015
C>
C>    @date 8-2-23 George Schoendorff
C>     - Added support for h and i functions
C>
C>    @param   SS     d.p. storage for overlap integrals
C>    @param   TT     d.p. storage for kinetic energy integrals
C>    @param   S      q.p. storage for overlap integrals
C>    @param   T      q.p. storage for kinetic energy integrals
C>    @param   ISUB   indexes the atom being worked on.
C>    @param   LL1    size of atomic basis
C>    @param   LL2    size of the four integral arrays
C>    @param   MAXLL2 largest size of the 4 integral arrays
C>    @param   SOME   set .TRUE. to see debug printing
C>
      SUBROUTINE QSANDTL(SS,TT,S,T,ISUB,LL1,LL2,MAXLL2,SOME)
      USE prec, ONLY: qp
      use mx_limits, only: mxsh,mxgtot,mxatm,mxang,mxang2,mxang3,mxdim
      use constants, only: ix,iy,iz,jx,jy,jz
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION MOROKM
C
      LOGICAL SOME,IANDJ,NORM,DOUBLE,GOPARR,DSKWRK,MASWRK
      REAL(kind=qp) S(MAXLL2),T(MAXLL2),SBLK(MXDIM),TBLK(MXDIM),
     *        FT(MXDIM),DIJ(MXDIM),
     *        XIN(MXANG3),YIN(MXANG3),ZIN(MXANG3),
     *        XINT,YINT,ZINT,TAA,X0,Y0,Z0,XI,YI,ZI,XJ,YJ,ZJ,
     *        CSI,CPI,CDI,CFI,CGI,CHI,CII,CSJ,CPJ,CDJ,CFJ,CGJ,CHJ,CIJ,
     *        QTOL,RR,ARRI,AXI,AYI,AZI,AI,AJ,AA,AA1,DUM,DUM1,DUM2,T1,T2,
     *        YZ,AX,AY,AZ,FAC,ZERO,PT5,ONE,TWO,THREE,FIVE,SEVEN,NINE,
     *        ELEVEN,THIRTEEN,FIFTEEN,SQRT3,SQRT5,SQRT7,SQRT9,SQRT11,
     *        RLN10
C
      DIMENSION IJX(MXDIM),IJY(MXDIM),IJZ(MXDIM),SS(*),TT(*)
C
      PARAMETER (NDQ=2)
C
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON /NSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),
     *                CF(MXGTOT),CG(MXGTOT),CH(MXGTOT),CI(MXGTOT),
     *                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),
     *                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
      COMMON /OUTPUT/ NPRINT,ITOL,ICUT,NORMF,NORMP,NOPK
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /ROOT  / XX,U(13),W(13),NROOTS
      COMMON /RUNOPT/ RUNTYP,EXETYP,NEVALS,NGLEVL,NHLEVL
      COMMON /SYMIND/ TOL,II,JJ,LIT,LJT,MINI,MINJ,MAXI,MAXJ,IANDJ
      COMMON /WFNOPT/ SCFTYP,VBTYP,DFTYPE,TDDFTYP,CITYP,CCTYP,
     *                MPLEVL,MPCTYP
      COMMON /SETATM/ NSHL(MXATM),NBSF(MXATM),NNUM(MXATM),NLOC(MXATM),
     *                NUNV,NUNPV,ISUBX
C
      PARAMETER (ZERO=0.0_qp)
      PARAMETER (PT5=0.5_qp)
      PARAMETER (ONE=1.0_qp)
      PARAMETER (TWO=2.0_qp)
      PARAMETER (THREE=3.0_qp)
      PARAMETER (FIVE=5.0_qp)
      PARAMETER (SEVEN=7.0_qp)
      PARAMETER (NINE=9.0_qp)
      PARAMETER (ELEVEN=11.0_qp)
      PARAMETER (THIRTEEN=13.0_qp)
      PARAMETER (FIFTEEN=15.0_qp)
C
      DATA MOROKM/8HMOROKUMA/,NONE/4HNONE/
C
C     ----- COMPUTE CONVENTIONAL S, AND T INTEGRALS -----
C             ----- IN QUADRUPLE PRECISION -----
C     HERMIT COEFFICIENTS ARE ROUNDED TO DOUBLE PRECISION, HOWEVER,
C     THE LOSS OF ACCURACY TO LARGE EXTENT COMES FROM MULTIPLYING
C     LARGE EXPONENTS IN THIS ROUTINE AND WE GET FULL ACCURACY HERE.
C
C     PI212=SQRT(ONE/ACOS(ZERO))
      SQRT3=SQRT(THREE)
      SQRT5=SQRT(FIVE)
      SQRT7=SQRT(SEVEN)
      SQRT9=THREE
      SQRT11=SQRT(ELEVEN)
      RLN10=LOG(NINE+ONE)
C     UNCON=RMETHOD.NE.ANONE.AND.MOD(MODQR,2).EQ.1
      QTOL = RLN10*ITOL
      NORM = NORMF .NE. 1 .OR. NORMP .NE. 1
C
C     ----- MOPAC INTEGRALS ARE DONE ELSEWHERE -----
C
      IF(MPCTYP.NE.NONE) THEN
C        CALL MPCINT
         CALL ABRT
         RETURN
      END IF
C
      IF (RUNTYP.EQ.MOROKM) THEN
         CALL ABRT
C        CALL STINT1(ISTART,IEND,JSTART,LOCIJ,NATST,NATED,ISAVE,L1,L2)
      ELSE
         ISTART = 1
         IEND   = NSHELL
         JSTART = 1
         LOCIJ  = 0
C        NATST  = 1
C        NATED  = NAT+NCHMAT
C        IF(ISEPS) NATED = NAT+NPS
C        ISAVE  = 0
C        L1 = NUM
C        IF(UNCON) L1=NUMU
C        L2 = (L1*(L1+1))/2
      END IF
C
      CALL STLUT(ISTART,IEND,JSTART,JEND,NATST,NATED,ISUB,ISUB)
C
      CALL QCLR(S,1,MAXLL2)
      CALL QCLR(T,1,MAXLL2)
C
C     ----- INTIALIZE PARALLEL -----
C
C     IPCOUNT = ME - 1
C
C     ----- I SHELL -----
C
      DO 720 II = ISTART,IEND
         I = KATOM(II)
         XI = C(1,I)
         YI = C(2,I)
         ZI = C(3,I)
         I1 = KSTART(II)
         I2 = I1+KNG(II)-1
         LIT = KTYPE(II)
         MINI = KMIN(II)
         MAXI = KMAX(II)
         LOCI = KLOC(II)-MINI-LOCIJ-NLOC(ISUB)+1
C
C     ----- J SHELL -----
C
         DO 700 JJ = JSTART,II
C
C     ----- GO PARALLEL! -----
C     PENDING GLOBAL SUM FOR QUADRUPLE PRECISION
C
C           IF (GOPARR) THEN
C              IPCOUNT = IPCOUNT + 1
C              IF (MOD(IPCOUNT,NPROC).NE.0) GO TO 700
C           END IF
            J = KATOM(JJ)
            XJ = C(1,J)
            YJ = C(2,J)
            ZJ = C(3,J)
            J1 = KSTART(JJ)
            J2 = J1+KNG(JJ)-1
            LJT = KTYPE(JJ)
            MINJ = KMIN(JJ)
            MAXJ = KMAX(JJ)
            LOCJ = KLOC(JJ)-MINJ-LOCIJ-NLOC(ISUB)+1
            NROOTS = (LIT+LJT-2)/2+1
            RR = (XI-XJ)**2+(YI-YJ)**2+(ZI-ZJ)**2
            IANDJ = II .EQ. JJ
C
C     ----- PREPARE INDICES FOR PAIRS OF (I,J) FUNCTIONS
C
            IJ = 0
            MAX = MAXJ
            DO 160 I = MINI,MAXI
               NX = IX(I)
               NY = IY(I)
               NZ = IZ(I)
               IF (IANDJ) MAX = I
               DO 140 J = MINJ,MAX
                  IJ = IJ+1
                  IJX(IJ) = NX+JX(J)
                  IJY(IJ) = NY+JY(J)
                  IJZ(IJ) = NZ+JZ(J)
                  IF (J.LE.1) FT(IJ) = THREE
                  IF ((J.GT.1).AND.(J.LE.4)) FT(IJ) = FIVE
                  IF ((J.GT.4).AND.(J.LE.10)) FT(IJ) = SEVEN
                  IF ((J.GT.10).AND.(J.LE.20)) FT(IJ) = NINE
                  IF ((J.GT.20).AND.(J.LE.35)) FT(IJ) = ELEVEN
                  IF ((J.GT.35).AND.(J.LE.56)) FT(IJ) = THIRTEEN
                  IF (J.GT.56) FT(IJ) = FIFTEEN
  140          CONTINUE
  160       CONTINUE
C
            CALL QCLR(SBLK,1,IJ)
            CALL QCLR(TBLK,1,IJ)
C
C     ----- I PRIMITIVE
C
            JGMAX = J2
            DO 520 IG = I1,I2
               AI = EX(IG)
               ARRI = AI*RR
               AXI = AI*XI
               AYI = AI*YI
               AZI = AI*ZI
               CSI = CS(IG)
               CPI = CP(IG)
               CDI = CD(IG)
               CFI = CF(IG)
               CGI = CG(IG)
               CHI = CH(IG)
               CII = CI(IG)
C
C     ----- J PRIMITIVE
C
               IF (IANDJ) JGMAX = IG
               DO 500 JG = J1,JGMAX
                  AJ = EX(JG)
                  AA = AI+AJ
                  AA1 = ONE/AA
                  DUM = AJ*ARRI*AA1
                  IF (DUM .GT. QTOL) GO TO 500
                  FAC = EXP(-DUM)
                  CSJ = CS(JG)
                  CPJ = CP(JG)
                  CDJ = CD(JG)
                  CFJ = CF(JG)
                  CGJ = CG(JG)
                  CHJ = CH(JG)
                  CIJ = CI(JG)
                  AX = (AXI+AJ*XJ)*AA1
                  AY = (AYI+AJ*YJ)*AA1
                  AZ = (AZI+AJ*ZJ)*AA1
C
C     ----- DENSITY FACTOR
C
                  DOUBLE=IANDJ.AND.IG.NE.JG
                  MAX = MAXJ
                  NN = 0
                  DUM1 = ZERO
                  DUM2 = ZERO
                  DO 220 I = MINI,MAXI
                     IF (I.EQ.1) DUM1=CSI*FAC
                     IF (I.EQ.2) DUM1=CPI*FAC
                     IF (I.EQ.5) DUM1=CDI*FAC
                     IF ((I.EQ. 8).AND.NORM) DUM1=DUM1*SQRT3
                     IF (I.EQ.11) DUM1=CFI*FAC
                     IF ((I.EQ.14).AND.NORM) DUM1=DUM1*SQRT5
                     IF ((I.EQ.20).AND.NORM) DUM1=DUM1*SQRT3
                     IF (I.EQ.21) DUM1=CGI*FAC
                     IF ((I.EQ.24).AND.NORM) DUM1=DUM1*SQRT7
                     IF ((I.EQ.30).AND.NORM) DUM1=DUM1*SQRT5/SQRT3
                     IF ((I.EQ.33).AND.NORM) DUM1=DUM1*SQRT3
                     IF ( I.EQ.36)           DUM1 = CHI*FAC
                     IF ((I.EQ.39).AND.NORM) DUM1 = DUM1*SQRT9
                     IF ((I.EQ.45).AND.NORM) DUM1 = DUM1*SQRT7/SQRT3
                     IF ((I.EQ.51).AND.NORM) DUM1 = DUM1*SQRT3
                     IF ((I.EQ.54).AND.NORM) DUM1 = DUM1*SQRT5/SQRT3
                     IF ( I.EQ.57)           DUM1 = CII*FAC
                     IF ((I.EQ.60).AND.NORM) DUM1 = DUM1*SQRT11
                     IF ((I.EQ.66).AND.NORM) DUM1 = DUM1*SQRT9/SQRT3
                     IF ((I.EQ.72).AND.NORM) DUM1 = DUM1*SQRT3
                     IF ((I.EQ.75).AND.NORM) THEN
                       DUM1 = DUM1*SQRT7/(SQRT3*SQRT5)
                     END IF
                     IF ((I.EQ.78).AND.NORM) DUM1 = DUM1*SQRT5
                     IF ((I.EQ.84).AND.NORM) DUM1 = DUM1*SQRT5/SQRT3
                     IF (IANDJ) MAX = I
                     DO 200 J = MINJ,MAX
                        IF (J.EQ.1) THEN
                           DUM2=DUM1*CSJ
                           IF (DOUBLE) THEN
                              IF (I.LE.1) THEN
                                 DUM2=DUM2+DUM2
                              ELSE
                                 DUM2=DUM2+CSI*CPJ*FAC
                              END IF
                           END IF
                        ELSE IF (J.EQ.2) THEN
                           DUM2=DUM1*CPJ
                           IF (DOUBLE) DUM2=DUM2+DUM2
                        ELSE IF (J.EQ.5) THEN
                           DUM2=DUM1*CDJ
                           IF (DOUBLE) DUM2=DUM2+DUM2
                        ELSE IF ((J.EQ.8).AND.NORM) THEN
                           DUM2=DUM2*SQRT3
                        ELSE IF (J.EQ.11) THEN
                           DUM2=DUM1*CFJ
                           IF (DOUBLE) DUM2=DUM2+DUM2
                        ELSE IF ((J.EQ.14).AND.NORM) THEN
                           DUM2=DUM2*SQRT5
                        ELSE IF ((J.EQ.20).AND.NORM) THEN
                           DUM2=DUM2*SQRT3
                        ELSE IF (J.EQ.21) THEN
                           DUM2=DUM1*CGJ
                           IF (DOUBLE) DUM2=DUM2+DUM2
                        ELSE IF ((J.EQ.24).AND.NORM) THEN
                           DUM2=DUM2*SQRT7
                        ELSE IF ((J.EQ.30).AND.NORM) THEN
                           DUM2=DUM2*SQRT5/SQRT3
                        ELSE IF ((J.EQ.33).AND.NORM) THEN
                           DUM2=DUM2*SQRT3
                        ELSE IF ( J.EQ.36) THEN
                           DUM2 = DUM1*CHJ
                           IF (DOUBLE) DUM2=DUM2+DUM2
                        ELSE IF ((J.EQ.39).AND.NORM) THEN
                           DUM2 = DUM2*SQRT9
                        ELSE IF ((J.EQ.45).AND.NORM) THEN
                           DUM2 = DUM2*SQRT7/SQRT3
                        ELSE IF ((J.EQ.51).AND.NORM) THEN
                           DUM2 = DUM2*SQRT3
                        ELSE IF ((J.EQ.54).AND.NORM) THEN
                           DUM2 = DUM2*SQRT5/SQRT3
                        ELSE IF ( J.EQ.57) THEN
                           DUM2 = DUM1*CIJ
                           IF (DOUBLE) DUM2=DUM2+DUM2
                        ELSE IF ((J.EQ.60).AND.NORM) THEN
                           DUM2 = DUM2*SQRT11
                        ELSE IF ((J.EQ.66).AND.NORM) THEN
                           DUM2 = DUM2*SQRT9/SQRT3
                        ELSE IF ((J.EQ.72).AND.NORM) THEN
                           DUM2 = DUM2*SQRT3
                        ELSE IF ((J.EQ.75).AND.NORM) THEN
                           DUM2 = DUM2*SQRT7/(SQRT3*SQRT5)
                        ELSE IF ((J.EQ.78).AND.NORM) THEN
                           DUM2 = DUM2*SQRT5
                        ELSE IF ((J.EQ.84).AND.NORM) THEN
                           DUM2 = DUM2*SQRT5/SQRT3
                        END IF
                        NN = NN+1
                        DIJ(NN) = DUM2
  200                CONTINUE
  220             CONTINUE
C
C     ----- OVERLAP AND KINETIC ENERGY
C
                  TAA = SQRT(AA1)
                  T1 = -TWO*AJ*AJ*TAA
                  T2 = -PT5*TAA
                  X0 = AX
                  Y0 = AY
                  Z0 = AZ
                  IN = -MXANG
                  DO 320 I = 1,LIT
                     IN = IN+MXANG
                     NI = I
                     DO 300 J = 1,LJT
                        JN = IN+J
                        NJ = J
                        CALL QSTVINT(XINT,YINT,ZINT,
     *                          TAA,X0,Y0,Z0,XI,YI,ZI,XJ,YJ,ZJ,NI,NJ)
                        XIN(JN) = XINT*TAA
                        YIN(JN) = YINT*TAA
                        ZIN(JN) = ZINT*TAA
                        NJ = J+2
                        CALL QSTVINT(XINT,YINT,ZINT,
     *                          TAA,X0,Y0,Z0,XI,YI,ZI,XJ,YJ,ZJ,NI,NJ)
                        XIN(JN+MXANG2) = XINT*T1
                        YIN(JN+MXANG2) = YINT*T1
                        ZIN(JN+MXANG2) = ZINT*T1
                        NJ = J-2
                        IF (NJ .GT. 0) THEN
                           CALL QSTVINT(XINT,YINT,ZINT,
     *                          TAA,X0,Y0,Z0,XI,YI,ZI,XJ,YJ,ZJ,NI,NJ)
                        ELSE
                           XINT = ZERO
                           YINT = ZERO
                           ZINT = ZERO
                        END IF
                        N = (J-1)*(J-2)
                        DUM = N * T2
                        XIN(JN+MXANG2*2) = XINT*DUM
                        YIN(JN+MXANG2*2) = YINT*DUM
                        ZIN(JN+MXANG2*2) = ZINT*DUM
  300                CONTINUE
  320             CONTINUE
                  DO 340 I = 1,IJ
                     NX = IJX(I)
                     NY = IJY(I)
                     NZ = IJZ(I)
                     YZ = YIN(NY)*ZIN(NZ)
                     DUM = YZ*XIN(NX)
                     DUM1 = (XIN(NX+MXANG2)+XIN(NX+MXANG2*2))*YZ+
     *                      (YIN(NY+MXANG2)+YIN(NY+MXANG2*2))
     *                      *XIN(NX)*ZIN(NZ)+
     *                      (ZIN(NZ+MXANG2)+ZIN(NZ+MXANG2*2))
     *                      *XIN(NX)*YIN(NY)
                     SBLK(I) = SBLK(I) + DIJ(I)*DUM
                     TBLK(I) = TBLK(I) + DIJ(I)*(DUM*AJ*FT(I)+DUM1)
  340             CONTINUE
C
C     ----- END OF PRIMITIVE LOOPS -----
C
  500          CONTINUE
  520       CONTINUE
C
C     ----- COPY BLOCK OVERLAP AND KINETIC ENERGY MATRICES
C
            MAX = MAXJ
            NN = 0
            DO 620 I = MINI,MAXI
               LI = LOCI+I
               IN = (LI*(LI-1))/2
               IF (IANDJ) MAX = I
               DO 600 J = MINJ,MAX
                  LJ = LOCJ+J
                  JN = LJ+IN
                  NN = NN+1
                  S(JN) = SBLK(NN)
                  T(JN) = TBLK(NN)
  600          CONTINUE
  620       CONTINUE
C
C     ----- END OF SHELL LOOPS -----
C
  700    CONTINUE
  720 CONTINUE
C
C     ----- SUM UP PARTIAL CONTRIBUTIONS IF PARALLEL -----
C
      IF (GOPARR) THEN
C        CALL DDI_GSUMF(911,S,L2)
C        CALL DDI_GSUMF(912,T,L2)
      END IF
C
C     ----- SAVE S, AND T MATRICES ON THE DAF -----
C
      CALL DAWRIT(IDAF,IODA,S,MAXLL2*NDQ,394,0)
      CALL DAWRIT(IDAF,IODA,T,MAXLL2*NDQ,395,0)
C
C     REWRITE THE OVERLAP AND KINETIC INTEGRALS
C
      CALL QDCOPY(LL2,S,1,SS,1)
      CALL QDCOPY(LL2,T,1,TT,1)
C
C     NO NEEDS TO SAVE S & T TO DAF
C     CALL DAWRIT(IDAF,IODA,SS,LL2,12,0)
C     CALL DAWRIT(IDAF,IODA,TT,LL2,13,0)
C
      IF(SOME) THEN
         CALL PRTRIEXL(SS,LL1,'OVL.MAT(Q)')
         CALL PRTRIEXL(TT,LL1,'KIN.MAT(Q)')
      END IF
C
C     IF (ISAVE.EQ.1) THEN
C        CALL DAWRIT(IDAF,IODA,H,LL2,311,0)
C        CALL DAWRIT(IDAF,IODA,S,LL2,312,0)
C     END IF
      IF(SOME) WRITE(IW,*)
     *     'QUADRUPLE PRECISION RELATIVISTIC INTEGRALS WILL BE USED.'
      RETURN
      END
C
C*MODULE INT1    *DECK QSTVINT
C>
C> @date 8-1-23 George Schoendorff
C>  - Added support for h and i functions
C>
      SUBROUTINE QSTVINT(XINT,YINT,ZINT,
     *                            T,X0,Y0,Z0,XI,YI,ZI,XJ,YJ,ZJ,NI,NJ)
      USE gausshermite, ONLY: H, W
      USE prec, ONLY: qp
C
      IMPLICIT NONE
      REAL(kind=qp), INTENT(OUT) :: XINT,YINT,ZINT
      REAL(kind=qp), INTENT(IN)  :: T,X0,Y0,Z0,XI,YI,ZI,XJ,YJ,ZJ
      INTEGER,       INTENT(IN)  :: NI,NJ
C
      INTEGER, PARAMETER :: MIN(10) = (/1,2,4,7,11,16,22,29,37,46/)
      INTEGER, PARAMETER :: MAX(10) = (/1,3,6,10,15,21,28,36,45,55/)
C
      REAL(kind=qp), PARAMETER :: ZERO=0.0_qp
C
      INTEGER :: I,NPTS,IMIN,IMAX
      REAL(kind=qp) :: DUM, PTX, PTY, PTZ
      REAL(kind=qp) :: PX, PY, PZ, AX, AY, AZ, BX, BY, BZ
C
C     ----- GAUSS-HERMITE QUADRATURE USING MINIMUM POINT FORMULA -----
C     HERMIT COEFFICIENTS ARE IN DOUBLE PRECISION!
C
      XINT = ZERO
      YINT = ZERO
      ZINT = ZERO
      NPTS = (NI+NJ-2)/2+1
      IMIN = MIN(NPTS)
      IMAX = MAX(NPTS)
      DO I = IMIN,IMAX
         DUM = W(I)
         PX = DUM
         PY = DUM
         PZ = DUM
         DUM = H(I)*T
         PTX = DUM+X0
         PTY = DUM+Y0
         PTZ = DUM+Z0
         AX = PTX-XI
         AY = PTY-YI
         AZ = PTZ-ZI
         BX = PTX-XJ
         BY = PTY-YJ
         BZ = PTZ-ZJ
         GO TO (170,160,150,140,130,120,110,100),NI
  100       PX = PX*AX
            PY = PY*AY
            PZ = PZ*AZ
  110       PX = PX*AX
            PY = PY*AY
            PZ = PZ*AZ
  120       PX = PX*AX
            PY = PY*AY
            PZ = PZ*AZ
  130       PX = PX*AX
            PY = PY*AY
            PZ = PZ*AZ
  140       PX = PX*AX
            PY = PY*AY
            PZ = PZ*AZ
  150       PX = PX*AX
            PY = PY*AY
            PZ = PZ*AZ
  160       PX = PX*AX
            PY = PY*AY
            PZ = PZ*AZ
  170       GO TO (290,280,270,260,250,240,230,220,210,200),NJ
  200          PX = PX*BX
               PY = PY*BY
               PZ = PZ*BZ
  210          PX = PX*BX
               PY = PY*BY
               PZ = PZ*BZ
  220          PX = PX*BX
               PY = PY*BY
               PZ = PZ*BZ
  230          PX = PX*BX
               PY = PY*BY
               PZ = PZ*BZ
  240          PX = PX*BX
               PY = PY*BY
               PZ = PZ*BZ
  250          PX = PX*BX
               PY = PY*BY
               PZ = PZ*BZ
  260          PX = PX*BX
               PY = PY*BY
               PZ = PZ*BZ
  270          PX = PX*BX
               PY = PY*BY
               PZ = PZ*BZ
  280          PX = PX*BX
               PY = PY*BY
               PZ = PZ*BZ
  290          XINT = XINT+PX
               YINT = YINT+PY
               ZINT = ZINT+PZ
      END DO
      RETURN
      END
C*MODULE MTHLIB  *DECK QTFTRI
      SUBROUTINE QTFTRI(H,F,T,WRK,M,N,LDT)
      USE prec, ONLY: qp
C
      IMPLICIT NONE
C
      INTEGER, INTENT(IN) :: N, M, LDT
      REAL(kind=qp) :: H(*),F(*),T(LDT,M),WRK(N)
C
      INTEGER, PARAMETER :: MXROWS=5
      REAL(kind=qp), PARAMETER :: ZERO=0.0_qp
      REAL(kind=qp), PARAMETER :: SMALL=1.0E-25_qp
C
      INTEGER :: I, J, K, JJ, IK, IJ
      INTEGER :: JJMAX, IM1
      REAL(kind=qp) :: DUM, HIJ, TDUM
C
      REAL(kind=qp) :: QDOT
C
C     ----- TRANSFORM THE TRIANGULAR MATRIX F USING VECTORS T -----
C                      H = T-DAGGER * F * T
C     THE ORDER OF THE TRIANGULAR MATRICES H AND F ARE M AND N.
C
C
      IJ = 0
      DO J = 1,M,MXROWS
         JJMAX = MIN(M,J+MXROWS-1)
C
C             FIRST CALCULATE T-DAGGER TIMES -F-, A ROW AT A TIME
C
         DO JJ=J,JJMAX
            IK = 0
            DO I = 1,N
               IM1 = I-1
               DUM = ZERO
               TDUM = T(I,JJ)
               IF (IM1.GT.0) THEN
                  DO K = 1,IM1
                     IK = IK+1
                     WRK(K) = WRK(K)+F(IK)*TDUM
                     DUM = DUM+F(IK)*T(K,JJ)
                  END DO
               END IF
               IK = IK+1
               WRK(I) = DUM+F(IK)*TDUM
            END DO
C
C             THEN TAKE THAT ROW TIMES EVERY COLUMN IN -T-
C
            DO I = 1,JJ
               IJ = IJ+1
               HIJ = QDOT(N,T(1,I),1,WRK,1)
               IF(ABS(HIJ).LT.SMALL) HIJ=ZERO
               H(IJ)=HIJ
            END DO
         END DO
      END DO
C
      RETURN
      END
C*MODULE BLAS1   *DECK QDOT
      FUNCTION QDOT(N,DX,INCX,DY,INCY)
      USE prec, ONLY: qp
      IMPLICIT NONE
      REAL(kind=qp) QDOT
      INTEGER,       INTENT(IN) :: N, INCX, INCY
      REAL(kind=qp), INTENT(IN) :: DX(*),DY(*)
C
      REAL(kind=qp) :: DTEMP
      INTEGER :: I, IX, IY, M, MP1
C
C     FORMS THE DOT PRODUCT OF TWO VECTORS.
C           DOT = DX(I) * DY(I)
C     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      QDOT = 0.0_qp
      DTEMP = 0.0_qp
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20
C
C        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
C          NOT EQUAL TO 1
C
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO I = 1,N
        DTEMP = DTEMP + DX(IX)*DY(IY)
        IX = IX + INCX
        IY = IY + INCY
      END DO
      QDOT = DTEMP
      RETURN
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
C        CLEAN-UP LOOP
C
   20 M = MOD(N,5)
      IF( M .EQ. 0 ) GO TO 40
      DO I = 1,M
        DTEMP = DTEMP + DX(I)*DY(I)
      END DO
      IF( N .LT. 5 ) GO TO 60
   40 MP1 = M + 1
      DO I = MP1,N,5
        DTEMP = DTEMP + DX(I)*DY(I) + DX(I + 1)*DY(I + 1) +
     *   DX(I + 2)*DY(I + 2) + DX(I + 3)*DY(I + 3) + DX(I + 4)*DY(I + 4)
      END DO
   60 QDOT = DTEMP
      RETURN
      END
