C  6 Jun 18 - YN -  updates for FMO 5.3
C 18 Apr 16 - YN  - updates for FMO 5.2
C 29 Sep 14 - YN  - Implement Density-Functional Tight-Binding Method
C
C*MODULE DFTBSK    *DECK DFTB_SKTAB
C>
C>    @brief Read tabulated Slater-Koster transformation files
C>
C>    @details Read tabulated Slater-Koster transformation files
C>
C>    @author Yoshio Nishimoto
C>    - Sep, 2014- Subroutine written
C>    @date  Feb, 2016 - Yoshio Nishimoto
C>    - Changes due to the DFTBSK common block
C>
C>           --- INPUT ---
C>    @param I,J      Index of species of atom A and B
C>    @param SKPATH   Absolute path of Slater-Koster file for the
C>           pair of atom A-B
C>           --- OUTPUT ---
C>    @param SKHTAB   Slater-Koster parameters for H_{mn}^0
C>    @param SKSTAB   Slater-Koster parameters for S_{mn}
C>    @param SKSTAB   Slater-Koster parameters for diagonal H_{mn}^0
C>    @param SKGRID   Grid size of SK parameters
C>    @param REPCOEFF Coefficients for repulsive energy
C>    @param REPINTV  Intervals of discretized repulsive energy
C>    @param REPSHORT Coefficients for short-range repulsive energy
C>
C
      SUBROUTINE DFTB_SKTAB(I,J,SKPATH,SKHTAB,SKSTAB,SKSELF,SKGRID,
     *                      REPCOEFF,REPINTV,REPSHORT)
      use mx_limits, only: mxatm
C
      IMPLICIT NONE
C
      INTEGER,PARAMETER :: MXSPE=10
      INTEGER,PARAMETER :: MXTAB=1500, MXINT=150
C
      INTEGER, INTENT(IN) :: I, J
      CHARACTER(*), INTENT(IN) :: SKPATH
      DOUBLE PRECISION SKHTAB(10,MXTAB),SKSTAB(10,MXTAB),SKSELF(3),
     *                 SKGRID(2),REPCOEFF(6,MXINT),REPINTV(2,MXINT),
     *                 REPSHORT(3)
C
      COMMON /DFTBPR/ ETEMP,DFTBDP(MXSPE*14),DAMPXHE,HUBDER(MXSPE),
     *                ZREF(MXATM),SPNCST(6,MXSPE),SPE(MXATM),NSPE,
     *                MAXANG(MXATM),ISPE(MXATM),IND(MXATM+1),IDFTBD,
     *                PARAMDIR
      DOUBLE PRECISION :: ETEMP,DFTBDP,DAMPXHE,HUBDER,ZREF,SPNCST,SPE
      INTEGER :: NSPE,MAXANG,ISPE,IND,IDFTBD
C
      COMMON /DFTBRE/ REPCUT(MXSPE,MXSPE),NUMREP(MXSPE,MXSPE),
     *                LREPINTV(MXSPE,MXSPE),LREPSHORT(MXSPE,MXSPE),
     *                LREPCOEFF(MXSPE,MXSPE)
      DOUBLE PRECISION REPCUT,PARAMDIR
      INTEGER LREPCOEFF,LREPINTV,LREPSHORT,NUMREP
C
      COMMON /DFTBSK/ SKDIM(MXSPE,MXSPE),SKSPIN(MXSPE),QREFL(3,MXSPE),
     *                HUBBL(3,MXSPE),QREF(MXSPE),HUBB(MXSPE),SKCUT2,
     *                NEEDSK,LSKHTAB(MXSPE,MXSPE),LSKSTAB(MXSPE,MXSPE),
     *                LSKGRID(MXSPE,MXSPE),LSKSELF(MXSPE)
      DOUBLE PRECISION SKDIM,SKSPIN,QREFL,HUBBL,QREF,HUBB,SKCUT2
      INTEGER NEEDSK,LSKHTAB,LSKSTAB,LSKGRID,LSKSELF
C
      COMMON /IOFILE/ IR,IW,IP,IJK,IJKT,IDAF,NAV,IODA(950)
      INTEGER :: IR,IW,IP,IJK,IJKT,IDAF,NAV,IODA
C
      DOUBLE PRECISION :: SLKCUTOFF,QREFTMP

      INTEGER :: K, L, NFTDFTB, NSKIP

C     CHARACTER(4) :: DUM4
      CHARACTER(128) CHDUMMY
C
      nftdftb=21
      WRITE (IW,'(I3,I3," (",A8," - ",A8,") = ",A)')
     &         I,J,SPE(I),SPE(J),TRIM(ADJUSTL(SKPATH))
      OPEN (NFTDFTB,FILE=ADJUSTL(TRIM(SKPATH)),STATUS='OLD',
     *      FORM='FORMATTED',ACTION='READ')
      REWIND NFTDFTB
C       --- READ HEADERS ---
C     FIRST DISTINGUISH HOMO- OR HETERO-SPECIES
      IF (I.NE.J) THEN
        READ (nftdftb,*) SKGRID(1), SKDIM(I,J)
        SKGRID(2) = 1.0D+00/SKGRID(1)
      ELSE
        READ (nftdftb,*) SKGRID(1), SKDIM(I,J)
        SKGRID(2) = 1.0D+00/SKGRID(1)
        READ (nftdftb,*) (SKSELF(L),L= 1,3), SKSPIN(I),
     &    (HUBBL(4-L,I),L=1,3), (QREFL(4-L,I),L=1,3)
        HUBB(I) = HUBBL(1,I)
        QREFTMP = QREF(I)
        QREF(I) = 0.0D+00
        DO K = 1, 3
          QREF(I) = QREF(I) + QREFL(K,I)
          SKSPIN(I) = SKSPIN(I) + SKSELF(K)*QREFL(4-K,I)
        END DO
        IF (QREFTMP.NE.-1.0D+00) QREF(I) = QREFTMP
      END IF
      IF (MXTAB.LT.SKDIM(I,J)) WRITE (IW,*) "WARNING: MXTAB < SKDIM"
C       ----- READ HAMILTONIAN AND OVERLAP MATRICES -----
      DO K = 1, INT(SKDIM(I,J)+1.0D-01)
        READ (nftdftb,*) (SKHTAB(L,K), L = 1, 10),
     &                   (SKSTAB(L,K), L = 1, 10)
      END DO
  100 CONTINUE
      NSKIP = 0
      READ (nftdftb,'(A)') CHDUMMY
      CALL UPRCAS(CHDUMMY,6)
      IF (CHDUMMY.EQ.'SPLINE') THEN
        GOTO 110
      END IF
      NSKIP = NSKIP + 1
      GOTO 100
C
C       --- TRY TO READ REPULSIVE POTENTIALS ---
C
  110 CONTINUE
      IF (NSKIP.NE.0) WRITE (IW,*) NSKIP," LINES ARE SKIPPED"
      READ (nftdftb,*) NUMREP(I,J), REPCUT(I,J)
      IF (NUMREP(I,J).GT.MXINT) THEN
        WRITE(IW,*) "TOO MANY INTERVALS FOR PAIR", I, J
c       EXIT
      END IF
C       SHORT DISTANCE REGION
      READ  (nftdftb,*) (REPSHORT(K), K = 1, 3)
c       MEDIUM DISTANCE REGION
      DO K = 1, NUMREP(I,J)
        IF (K.LT.NUMREP(I,J)) THEN
          READ (nftdftb,*) REPINTV(1,K), REPINTV(2,K),
     &      REPCOEFF(1,K),REPCOEFF(2,K),REPCOEFF(3,K),REPCOEFF(4,K)
        ELSE
C           LONG DISTANCE REGION
          READ (nftdftb,*) REPINTV(1,K), REPINTV(2,K),
     &      REPCOEFF(1,K), REPCOEFF(2,K), REPCOEFF(3,K),
     &      REPCOEFF(4,K), REPCOEFF(5,K), REPCOEFF(6,K)
        END IF
      END DO
      IF (REPINTV(2,NUMREP(I,J)).NE.REPCUT(I,J)) THEN
        WRITE (IW,*) "ERROR IN DATAFILE"
        RETURN
      END IF
C     PRINT *, '(SKFILE FOR PAIR ',I,J,' :', SKFILE
      CLOSE (NFTDFTB)

C     Check the cutoff for Slater-Koster file

      SLKCUTOFF = SKGRID(1)*SKDIM(I,J) + 0.3D+00
      IF (SLKCUTOFF**2.GT.SKCUT2) SKCUT2 = SLKCUTOFF**2
C
C     K = 1
C     CALL DFTB_SETLAB(K,DUM4)
C
      RETURN
C
      END SUBROUTINE DFTB_SKTAB
C
C-----------------------------------------------------------------------
C*MODULE DFTBSK    *DECK DFTB_SETLAB
C>
C>    @brief Set label of AOs for DFTB
C>
C>    @details Set label (BFLAB = s, px, py, ...) of AOs for DFTB.
C>             This is a clone of SETLAB subroutine
C>             Also, check ranges of each shell (KMAX and KMIN)
C>
C>    @author Yoshio Nishimoto
C>    - Sep, 2014- Subroutine written
C>
C>           --- INPUT ---
C>    @param IARG Choose mode, currently only 1 is used
C>           --- OUTPUT ---
C>    @param ATMSYM Currently not used, but left as it is for future?
C
      SUBROUTINE DFTB_SETLAB(IARG,ATMSYM)
      use mx_limits, only: mxatm,mxsh,mxgtot,mxao
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      CHARACTER*4 ATMSYM(*)
      CHARACTER*4 BFNAM1(35)
      CHARACTER*4 ATMLAB(106),BONDF,LABEL
      CHARACTER*2 ATMNUM(100)
C
      PARAMETER (MXSPE=10)
C
      COMMON /ECP2  / CLP(400),ZLP(400),NLP(400),KFIRST(MXATM,6),
     *                KLAST(MXATM,6),LMAX(MXATM),LPSKIP(MXATM),
     *                IZCORE(MXATM)
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      INTEGER :: NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,IAN
      DOUBLE PRECISION :: ZAN,C
      COMMON /NSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),
     *                CF(MXGTOT),CG(MXGTOT),CH(MXGTOT),CI(MXGTOT),
     *                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),
     *                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
      COMMON /RUNLAB/ TITLE(10),ANAM(MXATM),BNAM(MXATM),BFLAB(MXAO)
C
      DATA BFNAM1/'  S ','  X ','  Y ','  Z ',
     *            ' XY ',' YZ ',' ZX ','XXYY',' ZZ ',' YZ ',
     *            ' XXX',' YYY',' ZZZ',' XXY',' XXZ',
     *            ' YYX',' YYZ',' ZZX',' ZZY',' XYZ',
     *            'XXXX','YYYY','ZZZZ','XXXY','XXXZ',
     *            'YYYX','YYYZ','ZZZX','ZZZY','XXYY',
     *            'XXZZ','YYZZ','XXYZ','YYXZ','ZZXY'/
      DATA BONDF/' BF '/
      DATA ATMNUM/' 1',' 2',' 3',' 4',' 5',' 6',' 7',' 8',' 9','10',
     *            '11','12','13','14','15','16','17','18','19','20',
     *            '21','22','23','24','25','26','27','28','29','30',
     *            '31','32','33','34','35','36','37','38','39','40',
     *            '41','42','43','44','45','46','47','48','49','50',
     *            '51','52','53','54','55','56','57','58','59','60',
     *            '61','62','63','64','65','66','67','68','69','70',
     *            '71','72','73','74','75','76','77','78','79','80',
     *            '81','82','83','84','85','86','87','88','89','90',
     *            '91','92','93','94','95','96','97','98','99','00'/
C
C                           104 TRUE ELEMENTS, AND 2 DUMMIES
C
      DATA ATMLAB/'H   ','HE  ','LI  ','BE  ','B   ','C   ',
     *            'N   ','O   ','F   ','NE  ','NA  ','MG  ',
     *            'AL  ','SI  ','P   ','S   ','CL  ','AR  ',
     *            'K   ','CA  ','SC  ','TI  ','V   ','CR  ',
     *            'MN  ','FE  ','CO  ','NI  ','CU  ','ZN  ',
     *            'GA  ','GE  ','AS  ','SE  ','BR  ','KR  ',
     *            'RB  ','SR  ','Y   ','ZR  ','NB  ','MO  ',
     *            'TC  ','RU  ','RH  ','PD  ','AG  ','CD  ',
     *            'IN  ','SN  ','SB  ','TE  ','I   ','XE  ',
     *            'CS  ','BA  ','LA  ','CE  ','PR  ','ND  ',
     *            'PM  ','SM  ','EU  ','GD  ','TB  ','DY  ',
     *            'HO  ','ER  ','TM  ','YB  ','LU  ','HF  ',
     *            'TA  ','W   ','RE  ','OS  ','IR  ','PT  ',
     *            'AU  ','HG  ','TL  ','PB  ','BI  ','PO  ',
     *            'AT  ','RN  ','FR  ','RA  ','AC  ','TH  ',
     *            'PA  ','U   ','NP  ','PU  ','AM  ','CM  ',
     *            'BK  ','CF  ','ES  ','FM  ','MD  ','NO  ',
     *            'LR  ','RF  ','X   ','BQ  '/
C
      COMMON /DFTBPR/ ETEMP,DFTBDP(MXSPE*14),DAMPXHE,HUBDER(MXSPE),
     *                ZREF(MXATM),SPNCST(6,MXSPE),SPE(MXATM),NSPE,
     *                MAXANG(MXATM),ISPE(MXATM),IND(MXATM+1),IDFTBD,
     *                PARAMDIR
      DOUBLE PRECISION :: ETEMP,DFTBDP,DAMPXHE,HUBDER,ZREF,SPNCST,SPE
      INTEGER :: NSPE,MAXANG,ISPE,IND,IDFTBD
C
      CHARACTER(8) :: TEMP
c
c     write(6,*) 'called labels',nat
      IF(IARG.EQ.1) GO TO 100
      IF(IARG.EQ.2) GO TO 200
      IF(IARG.EQ.3) GO TO 300
      CALL ABRT
C
C         ----- ZERO FILE FOR EXPONENTS AND COEFFICIENTS -----
C     NOW IN NSHEL COMMON BLOCK, KSTART AND KNG ARE NOT FILLED
C
      CALL VCLR(EX,1,MXGTOT)
      CALL VCLR(CS,1,MXGTOT)
      CALL VCLR(CP,1,MXGTOT)
      CALL VCLR(CD,1,MXGTOT)
      CALL VCLR(CF,1,MXGTOT)
      CALL VCLR(CG,1,MXGTOT)
      CALL VCLR(CH,1,MXGTOT)
      CALL VCLR(CI,1,MXGTOT)
C
C         ----- FILL BASIS FUNCTION SYMBOL TABLE -----
C
  100 IJ = 0
      KLOC(1) = 1
      DO I = 1, NAT
        DO J = 1, MAXANG(ISPE(I))
          IJ = IJ + 1
          KATOM(IJ) = I
          KTYPE(IJ) = J
          IF (J.EQ.1) THEN
            KMIN(IJ) = 1
            KMAX(IJ) = 1
          ELSE IF (J.EQ.2) THEN
            KMIN(IJ) = 2
            KMAX(IJ) = 4
          ELSE IF (J.EQ.3) THEN
            KMIN(IJ) = 5
            KMAX(IJ) = 9
          END IF
          IF (IJ.NE.1)
     &        KLOC(IJ) = KLOC(IJ-1) + KMAX(IJ-1) - KMIN(IJ-1) + 1
        END DO
      END DO
C
      N = 0
      DO II = 1, NSHELL
        IAT = KATOM(II)
        J = INT(ZAN(IAT))
        IF (J .GT. 106) J = 106
        IF (J .GT. 0) LABEL=ATMLAB(J)
        IF (J .LE. 0) LABEL=BONDF
        MINI = KMIN(II)
        MAXI = KMAX(II)
        DO I = MINI,MAXI
          N = N+1
          TEMP(1:2) = LABEL(1:2)
          K = MOD(IAT,100)
          IF (K.EQ.0) THEN
            TEMP(3:4) = '00'
          ELSE
            TEMP(3:4) = ATMNUM(K)
          END IF
          TEMP(5:8) = BFNAM1(I)
          CALL PUTCHARA(TEMP,BFLAB(N))
        END DO
      END DO
      RETURN
C
C        ----- FILL ATOM SYMBOL TABLE -----
C
  200 CONTINUE
      DO 210 IAT=1,NAT
         NUCZ = INT(ZAN(IAT)) + IZCORE(IAT)
         ATMSYM(IAT) = ATMLAB(105)
         IF(NUCZ.GE.1  .AND.  NUCZ.LE.104) ATMSYM(IAT) = ATMLAB(NUCZ)
  210 CONTINUE
      RETURN
C
C        ----- RETURN THE PERIODIC TABLE -----
C
  300 CONTINUE
      DO 310 IAT=1,106
         ATMSYM(IAT) = ATMLAB(IAT)
  310 CONTINUE
      RETURN
      END
C
C-----------------------------------------------------------------------
C*MODULE DFTBSK    *DECK DFTB_PUTCHARA
C>
C>    @brief CHARA2 = CHARA1
C>
C>    @details CHARA2 = CHARA1
C>
C>    @author Yoshio Nishimoto
C>    - Sep, 2014- Subroutine written
C>
C>           --- INPUT ---
C>    @param CHARA1 Input character (length = 8)
C>           --- OUTPUT ---
C>    @param CHARA2 Output character (length = 8)
C
      SUBROUTINE PUTCHARA(CHARA1,CHARA2)

      IMPLICIT NONE

      CHARACTER(8), INTENT(IN ) :: CHARA1
      CHARACTER(8), INTENT(OUT) :: CHARA2

      CHARA2 = CHARA1

      END SUBROUTINE PUTCHARA
C
C-----------------------------------------------------------------------
C*MODULE DFTBSK    *DECK DFTB_SKODE
C>
C>    @brief Get Hamiltonian and overlap values in AO
C>
C>    @details Calculate Hamiltonian and overlap for the given distance
C>             This subroutine is used for numerical derivative and now
C>             only for zeroth order value
C>
C>    @author Yoshio Nishimoto
C>    - Sep, 2014- Subroutine written
C>
C>           --- INPUT ---
C>    @param I,J Index of species of atom A and B
C>    @param MAXANGI,MAXANGJ Max angular momentum of atom A and B
C>    @param VECIN Three dimensional vector
C>    @param SKGRIDIJ,SKGRIDJI Grid length in AU, specified in
C>           so-called as Slater-Koster parameter files
C>    @param SKDIMIJ,SKDIMJI Number of grids for combinations of A-B
C>           and B-A, specified in parameter file
C>    @param SKHTABIJ,SKHTABJI Array of non-perturbed Hamiltonian for
C>           combinations of A-B and B-A
C>    @param SKSTABIJ,SKSTABJI Array of non-perturbed overlap element
C>           for combinations of A-B and B-A
C>    @param SKSELF On-site interaction potential
C>    @param SKCUT2 Cut-off value: longer than this distance results in
C>           zero contribution for both Hamiltonian and overlap
C>    @param THAM,TOVER Whether Hamiltonian and/or overlap elements are
C>           calculated or not
C>           --- OUTPUT ---
C>    @param AU,BU Block Hamiltonian and overlap matrix elements
C
      SUBROUTINE DFTB_SKODE(I,J,MAXANGI,MAXANGJ,VECIN,SKGRIDIJ,SKGRIDJI,
     & SKDIMIJ,SKDIMJI,SKHTABIJ,SKHTABJI,SKSTABIJ,SKSTABJI,SKSELF,
     & SKCUT2,AU,BU,THAM,TOVER)

      IMPLICIT NONE
C
      INTEGER, PARAMETER :: MXTAB=1500
C
      INTEGER         , INTENT(IN)    :: I, J ! SPECIE NUMBER
      INTEGER         , INTENT(IN)    :: MAXANGI, MAXANGJ ! MAX ANGULAR MOMENTUM OF I AND J
      DOUBLE PRECISION, INTENT(IN)    :: VECIN(3)
      DOUBLE PRECISION, INTENT(IN)    :: SKGRIDIJ
      DOUBLE PRECISION, INTENT(IN)    :: SKGRIDJI
      DOUBLE PRECISION, INTENT(IN)    :: SKDIMIJ
      DOUBLE PRECISION, INTENT(IN)    :: SKDIMJI
      DOUBLE PRECISION, INTENT(IN)    :: SKHTABIJ(10,MXTAB)
      DOUBLE PRECISION, INTENT(IN)    :: SKHTABJI(10,MXTAB)
      DOUBLE PRECISION, INTENT(IN)    :: SKSTABIJ(10,MXTAB)
      DOUBLE PRECISION, INTENT(IN)    :: SKSTABJI(10,MXTAB)
      DOUBLE PRECISION, INTENT(IN)    :: SKSELF(3)
      DOUBLE PRECISION, INTENT(IN)    :: SKCUT2
      DOUBLE PRECISION, INTENT(OUT)   :: AU(9,9), BU(9,9) ! AU -> HAMILTONIAN, BU -> OVERLAP
C
      DOUBLE PRECISION :: VEC(6), VEC2(6)
      DOUBLE PRECISION :: DUMMY(9,9)
      DOUBLE PRECISION :: R2, R2I, RI
      DOUBLE PRECISION :: HIJ(10), HJI(10), SIJ(10), SJI(10)
      INTEGER  :: MINMAX, MAXMAX, L, INU
      LOGICAL  :: THAM, TOVER
C
C
      MINMAX = MIN(MAXANGI,MAXANGJ)
      MAXMAX = MAX(MAXANGI,MAXANGJ)
C
      R2 = 0.0D+00 ! : distance
      DO L = 1, 3
        VEC(L) = VECIN(L)
        VEC2(L) = VEC(L)*VEC(L)
        R2 = R2 + VEC2(L)
      END DO
C
      IF (R2.GE.SKCUT2) RETURN
C
      IF (R2.GE.1.0D-08) THEN
        R2I = 1.0D+00/R2
        RI = SQRT(R2I)
C       NORMALIZE
        DO L = 1, 3
          VEC(L) = VEC(L) * RI
          VEC(L+3) = VEC(L)
          VEC2(L) = VEC2(L) * R2I
          VEC2(L+3) = VEC2(L)
        END DO
C
        IF (MAXMAX.LE.1) THEN
          INU = 10 ! ss interaction
        ELSE
          IF (MAXMAX.LE.2) THEN
            IF (MINMAX.LE.1) THEN
              INU = 9 ! ss and sp interaction
            ELSE
              INU = 6 ! ss, sp, pp interaction
            END IF
          ELSE
            IF (MINMAX.LE.1) THEN
              INU = 8 ! ss, sp, sd interaction
            ELSE
              IF (MINMAX.LE.2) THEN
                INU = 4 ! ss, sp, pp, sd, pd interaction
              ELSE
                INU = 1 ! ss, sp, pp, sd, pd, dd interaction
              END IF
            END IF
          END IF
        END IF
C
        IF (THAM)  CALL DFTB_SKPAR(R2,INU,HIJ,SKGRIDIJ,SKDIMIJ,SKHTABIJ)
        IF (TOVER) CALL DFTB_SKPAR(R2,INU,SIJ,SKGRIDIJ,SKDIMIJ,SKSTABIJ)
C
C s Interaction
C
        IF (THAM)  CALL DFTB_SKSS(HIJ(10),AU(1,1))
        IF (TOVER) CALL DFTB_SKSS(SIJ(10),BU(1,1))
        IF (MAXMAX.LE.1) RETURN
C
        IF (THAM)  CALL DFTB_SKPAR(R2,INU,HJI,SKGRIDJI,SKDIMJI,SKHTABJI)
        IF (TOVER) CALL DFTB_SKPAR(R2,INU,SJI,SKGRIDJI,SKDIMJI,SKSTABJI)
C
C p Interaction
C
        IF (MINMAX.GE.2) THEN ! (I,J) = (p,p)
          IF (THAM) THEN
            CALL DFTB_SKPP(VEC,VEC2,HIJ(6:7),AU(2:4,2:4))
            CALL DFTB_SKSP(VEC,HIJ(9),AU(1,2:4),AU(2:4,1))
          END IF
          IF (TOVER) THEN
            CALL DFTB_SKPP(VEC,VEC2,SIJ(6:7),BU(2:4,2:4))
            CALL DFTB_SKSP(VEC,SIJ(9),BU(1,2:4),BU(2:4,1))
          END IF
          IF (I.NE.J) THEN
            IF (THAM)  CALL DFTB_SKSP(VEC,HJI(9),DUMMY(1,2:4),
     *        AU(2:4,1))
            IF (TOVER) CALL DFTB_SKSP(VEC,SJI(9),DUMMY(1,2:4),
     *        BU(2:4,1))
          END IF
        ELSE ! (I,J) = (s,p) or (p,s)
          IF (MAXANGJ.GE.2) THEN ! (I,J) = (s,p)
            IF (THAM)  CALL DFTB_SKSP(VEC,HIJ(9),AU(1,2:4),
     :        AU(2:4,1))
            IF (TOVER) CALL DFTB_SKSP(VEC,SIJ(9),BU(1,2:4),
     *        BU(2:4,1))
          ELSE ! (I,J) = (p,s)
            IF (THAM)  CALL DFTB_SKSP(VEC,HJI(9),DUMMY(1,2:4),
     *        AU(2:4,1))
            IF (TOVER) CALL DFTB_SKSP(VEC,SJI(9),DUMMY(1,2:4),
     *        BU(2:4,1))
          END IF
        END IF
C
        IF (MAXMAX.LE.2) RETURN
C
C d Interaction
C
        IF (MINMAX.EQ.3) THEN
          IF (THAM) THEN
            CALL DFTB_SKDD(VEC,VEC2,HIJ(1:3),AU(5:9,5:9))
            CALL DFTB_SKSD(VEC,VEC2,HIJ(8),AU(1,5:9),AU(5:9,1))
            CALL DFTB_SKPD(VEC,VEC2,HIJ(4:5),AU(2:4,5:9),AU(5:9,2:4))
          END IF
          IF (TOVER) THEN
            CALL DFTB_SKDD(VEC,VEC2,SIJ(1:3),BU(5:9,5:9))
            CALL DFTB_SKSD(VEC,VEC2,SIJ(8),BU(1,5:9),BU(5:9,1))
            CALL DFTB_SKPD(VEC,VEC2,SIJ(4:5),BU(2:4,5:9),BU(5:9,2:4))
          END IF
C
          IF (I.NE.J) THEN
            IF (THAM) THEN
              CALL DFTB_SKSD(VEC,VEC2,HJI(8),DUMMY(1,5:9),AU(5:9,1))
              CALL DFTB_SKPD(VEC,VEC2,HJI(4:5),DUMMY(2:4,5:9),
     *          AU(5:9,2:4))
            END IF
            IF (TOVER) THEN
              CALL DFTB_SKSD(VEC,VEC2,SJI(8),DUMMY(1,5:9),BU(5:9,1))
              CALL DFTB_SKPD(VEC,VEC2,SJI(4:5),DUMMY(2:4,5:9),
     *          BU(5:9,2:4))
            END IF
          END IF
        ELSE
          IF (MAXANGI.EQ.1) THEN
            IF (THAM)  CALL DFTB_SKSD(VEC,VEC2,HIJ(8),AU(1,5:9),
     *        AU(5:9,1))
            IF (TOVER) CALL DFTB_SKSD(VEC,VEC2,SIJ(8),BU(1,5:9),
     *        BU(5:9,1))
          ELSE
            IF (MAXANGI.EQ.2) THEN
              IF (THAM) THEN
                CALL DFTB_SKPD(VEC,VEC2,HIJ(4:5),AU(2:4,5:9),
     *            AU(5:9,2:4))
                CALL DFTB_SKSD(VEC,VEC2,HIJ(8),AU(1,5:9),AU(5:9,1))
              END IF
              IF (TOVER) THEN
                CALL DFTB_SKPD(VEC,VEC2,SIJ(4:5),BU(2:4,5:9),
     *            BU(5:9,2:4))
                CALL DFTB_SKSD(VEC,VEC2,SIJ(8),BU(1,5:9),BU(5:9,1))
              END IF
            ELSE
              IF (MAXANGJ.EQ.1) THEN
                IF (THAM)
     *            CALL DFTB_SKSD(VEC,VEC2,HJI(8),DUMMY(1,5:9),AU(5:9,1))
                IF (TOVER)
     *            CALL DFTB_SKSD(VEC,VEC2,SJI(8),DUMMY(1,5:9),BU(5:9,1))
              ELSE
                IF (THAM) THEN
                  CALL DFTB_SKSD(VEC,VEC2,HJI(8),DUMMY(1,5:9),AU(5:9,1))
                  CALL DFTB_SKPD(VEC,VEC2,HJI(4:5),DUMMY(2:4,5:9),
     *              AU(5:9,2:4))
                END IF
                IF (TOVER) THEN
                  CALL DFTB_SKSD(VEC,VEC2,SJI(8),DUMMY(1,5:9),BU(5:9,1))
                  CALL DFTB_SKPD(VEC,VEC2,SJI(4:5),DUMMY(2:4,5:9),
     *              BU(5:9,2:4))
                END IF
              END IF
            END IF
          END IF
        END IF
      ELSE IF (I.NE.J) THEN
        IF (I.NE.J) RETURN
      ELSE
C
C zero distance
C
C ON-SITE S HAMILTONIAN AND OVERLAP (ALWAYS 1.0)
        IF (THAM)  AU(1,1) = SKSELF(3)
        IF (TOVER) BU(1,1) = 1.0D+00
        IF (MAXANGI.LE.1) RETURN
C ON-SITE P HAMILTONIAN AND OVERLAP (ALWAYS 1.0)
         IF (THAM) THEN
           AU(2,2) = SKSELF(2)
           AU(3,3) = SKSELF(2)
           AU(4,4) = SKSELF(2)
         END IF
         IF (TOVER) THEN
           BU(2,2) = 1.0D+00
           BU(3,3) = 1.0D+00
           BU(4,4) = 1.0D+00
         END IF
        IF (MAXANGI.LE.2) RETURN
C ON-SITE D HAMILTONIAN AND OVERLAP (ALWAYS 1.0)
        IF (THAM) THEN
          AU(5,5) = SKSELF(1)
          AU(6,6) = SKSELF(1)
          AU(7,7) = SKSELF(1)
          AU(8,8) = SKSELF(1)
          AU(9,9) = SKSELF(1)
        END IF
        IF (TOVER) THEN
          BU(5,5) = 1.0D+00
          BU(6,6) = 1.0D+00
          BU(7,7) = 1.0D+00
          BU(8,8) = 1.0D+00
          BU(9,9) = 1.0D+00
        END IF
      END IF

      END SUBROUTINE DFTB_SKODE
C
C-----------------------------------------------------------------------
C*MODULE DFTBSK    *DECK DFTB_SKODEA
C>
C>    @brief Get Hamiltonian and overlap values in AO
C>
C>    @details Calculate Hamiltonian and overlap for the given distance
C>    This subroutine is used for analytical derivatives
C>
C>    @author Yoshio Nishimoto
C>    - Sep, 2014- Subroutine written
C>
C>           --- INPUT ---
C>    @param Same as DFTB_SKODE
C>    @param NDER Order of derivative to be calculated
C
      SUBROUTINE DFTB_SKODEA(SAME,I,J,NDER,KI,KJ,KK,KL,MAXANGI,MAXANGJ,
     &                       VECIN,SKGRIDIJ,SKGRIDJI,SKDIMIJ,SKDIMJI,
     *                       SKHTABIJ,SKHTABJI,SKSTABIJ,SKSTABJI,SKSELF,
     *                       SKCUT2,AU,BU,THAM,TOVER,SKP)

      IMPLICIT NONE

      INTEGER,PARAMETER :: MXTAB=1500

      LOGICAL         , INTENT(IN)    :: SAME
      INTEGER         , INTENT(IN)    :: I, J
      INTEGER         , INTENT(IN)    :: KI,KJ,KK,KL,NDER
      INTEGER         , INTENT(IN)    :: MAXANGI, MAXANGJ
      DOUBLE PRECISION, INTENT(IN)    :: VECIN(3)
      DOUBLE PRECISION, INTENT(IN)    :: SKGRIDIJ(2)
      DOUBLE PRECISION, INTENT(IN)    :: SKGRIDJI(2)
      DOUBLE PRECISION, INTENT(IN)    :: SKDIMIJ
      DOUBLE PRECISION, INTENT(IN)    :: SKDIMJI
      DOUBLE PRECISION, INTENT(IN)    :: SKHTABIJ(10,MXTAB)
      DOUBLE PRECISION, INTENT(IN)    :: SKHTABJI(10,MXTAB)
      DOUBLE PRECISION, INTENT(IN)    :: SKSTABIJ(10,MXTAB)
      DOUBLE PRECISION, INTENT(IN)    :: SKSTABJI(10,MXTAB)
      DOUBLE PRECISION, INTENT(IN)    :: SKSELF(3)
      DOUBLE PRECISION, INTENT(IN)    :: SKCUT2
      DOUBLE PRECISION, INTENT(OUT)   :: AU(9,9), BU(9,9)
      LOGICAL         , INTENT(IN)    :: THAM, TOVER
      LOGICAL         , INTENT(OUT)   :: SKP

      DOUBLE PRECISION :: VEC(6), VEC2(6)
      DOUBLE PRECISION :: DUMMY(9,9)
      DOUBLE PRECISION :: R2, R2I, RI, R
      DOUBLE PRECISION :: HIJ(10),HJI(10),SIJ(10),SJI(10),HIJ1(10),
     *  HJI1(10),SIJ1(10),SJI1(10),HIJ2(10),HJI2(10),SIJ2(10),SJI2(10),
     *  HIJ3(10),HJI3(10),SIJ3(10),SJI3(10),HIJ4(10),HJI4(10),SIJ4(10),
     *  SJI4(10)
      INTEGER  :: MINMAX, MAXMAX, INU, MAXIND, MAXPAR
C
      SKP = .FALSE.
      IF (SAME) THEN
C
C zero distance
C
C ON-SITE S HAMILTONIAN AND OVERLAP (ALWAYS 1.0)
        IF (THAM)  AU(1,1) = SKSELF(3)
        IF (TOVER) BU(1,1) = 1.0D+00
        IF (MAXANGI.LE.1) RETURN
C ON-SITE P HAMILTONIAN AND OVERLAP (ALWAYS 1.0)
         IF (THAM) THEN
           AU(2,2) = SKSELF(2)
           AU(3,3) = SKSELF(2)
           AU(4,4) = SKSELF(2)
         END IF
         IF (TOVER) THEN
           BU(2,2) = 1.0D+00
           BU(3,3) = 1.0D+00
           BU(4,4) = 1.0D+00
         END IF
        IF (MAXANGI.LE.2) RETURN
C ON-SITE D HAMILTONIAN AND OVERLAP (ALWAYS 1.0)
        IF (THAM) THEN
          AU(5,5) = SKSELF(1)
          AU(6,6) = SKSELF(1)
          AU(7,7) = SKSELF(1)
          AU(8,8) = SKSELF(1)
          AU(9,9) = SKSELF(1)
        END IF
        IF (TOVER) THEN
          BU(5,5) = 1.0D+00
          BU(6,6) = 1.0D+00
          BU(7,7) = 1.0D+00
          BU(8,8) = 1.0D+00
          BU(9,9) = 1.0D+00
        END IF
      ELSE
        VEC(1)  = VECIN(1)
        VEC(2)  = VECIN(2)
        VEC(3)  = VECIN(3)
        VEC2(1) = VEC(1)*VEC(1)
        VEC2(2) = VEC(2)*VEC(2)
        VEC2(3) = VEC(3)*VEC(3)
        R2      = VEC2(1)+VEC2(2)+VEC2(3)
C
        IF (R2.GE.SKCUT2) THEN
          SKP = .TRUE.
          RETURN
        END IF
C
C       NORMALIZE
        R       = SQRT(R2)
        RI      = 1.0D+00/R
        R2I     = RI*RI
C       R2I     = 1.0D+00/R2
C       RI      = SQRT(R2I)
        VEC(1)  = VEC(1)*RI
        VEC(2)  = VEC(2)*RI
        VEC(3)  = VEC(3)*RI
        VEC(4)  = VEC(1)
        VEC(5)  = VEC(2)
        VEC(6)  = VEC(3)
        VEC2(1) = VEC2(1)*R2I
        VEC2(2) = VEC2(2)*R2I
        VEC2(3) = VEC2(3)*R2I
        VEC2(4) = VEC2(1)
        VEC2(5) = VEC2(2)
        VEC2(6) = VEC2(3)
C
        MINMAX = MIN(MAXANGI,MAXANGJ)
        MAXMAX = MAX(MAXANGI,MAXANGJ)
C
        IF (MAXMAX.LE.1) THEN
          INU = 10 ! ss interaction
        ELSE
          IF (MAXMAX.LE.2) THEN
            IF (MINMAX.LE.1) THEN
              INU = 9 ! ss and sp interaction
            ELSE
              INU = 6 ! ss, sp, pp interaction
            END IF
          ELSE
            IF (MINMAX.LE.1) THEN
              INU = 8 ! ss, sp, sd interaction
            ELSE
              IF (MINMAX.LE.2) THEN
                INU = 4 ! ss, sp, pp, sd, pd interaction
              ELSE
                INU = 1 ! ss, sp, pp, sd, pd, dd interaction
              END IF
            END IF
          END IF
        END IF

        MAXIND = INT(SKDIMIJ + (0.3D+00*SKGRIDIJ(2) - 1.0D+00))
        MAXPAR = INT(R*SKGRIDIJ(2) + 1.0D+00)
        IF (THAM)  CALL DFTB_SKPARA(R,INU,HIJ,HIJ1,HIJ2,HIJ3,HIJ4,
     *                  SKGRIDIJ,SKDIMIJ,SKHTABIJ,NDER,MAXIND,MAXPAR)
        IF (TOVER) CALL DFTB_SKPARA(R,INU,SIJ,SIJ1,SIJ2,SIJ3,SIJ4,
     *                  SKGRIDIJ,SKDIMIJ,SKSTABIJ,NDER,MAXIND,MAXPAR)
C
C s Interaction
C
        IF (THAM)  CALL DFTB_SKSSA(NDER,KI,KJ,KK,KL,RI,VEC,HIJ(10),
     *                        HIJ1(10),HIJ2(10),HIJ3(10),HIJ4(10),
     *                        AU(1,1))
        IF (TOVER) CALL DFTB_SKSSA(NDER,KI,KJ,KK,KL,RI,VEC,SIJ(10),
     *                        SIJ1(10),SIJ2(10),SIJ3(10),SIJ4(10),
     *                        BU(1,1))
        IF (MAXMAX.LE.1) RETURN

C       this can be omitted?
C       IF ( (MAXANGI.GT.MAXANGJ) .OR.
C    &       (MAXANGI.EQ.MAXANGJ .AND. I.NE.J) ) THEN
        MAXIND = INT(SKDIMJI + (0.3D+00*SKGRIDJI(2) - 1.0D+00))
        MAXPAR = INT(R*SKGRIDJI(2) + 1.0D+00)
        IF (THAM)  CALL DFTB_SKPARA(R,INU,HJI,HJI1,HJI2,HJI3,HJI4,
     *                  SKGRIDJI,SKDIMJI,SKHTABJI,NDER,MAXIND,MAXPAR)
        IF (TOVER) CALL DFTB_SKPARA(R,INU,SJI,SJI1,SJI2,SJI3,SJI4,
     *                  SKGRIDJI,SKDIMJI,SKSTABJI,NDER,MAXIND,MAXPAR)
C       END IF
C
C p Interaction
C
        IF (MINMAX.GE.2) THEN ! (I,J) = (p,p)
          IF (THAM) THEN
            CALL DFTB_SKPPA(NDER,KI,KJ,KK,KL,RI,VEC,HIJ(6:7),HIJ1(6:7),
     *                 HIJ2(6:7),HIJ3(6:7),HIJ4(6:7),AU(2:4,2:4))
            CALL DFTB_SKSPA(NDER,KI,KJ,KK,KL,RI,VEC,HIJ(9),HIJ1(9),
     *                 HIJ2(9),HIJ3(9),HIJ4(9),AU(1,2:4),AU(2:4,1))
          END IF
          IF (TOVER) THEN
            CALL DFTB_SKPPA(NDER,KI,KJ,KK,KL,RI,VEC,SIJ(6:7),SIJ1(6:7),
     *                 SIJ2(6:7),SIJ3(6:7),SIJ4(6:7),BU(2:4,2:4))
            CALL DFTB_SKSPA(NDER,KI,KJ,KK,KL,RI,VEC,SIJ(9),SIJ1(9),
     *                 SIJ2(9),SIJ3(9),SIJ4(9),BU(1,2:4),BU(2:4,1))
          END IF
          IF (I.NE.J) THEN
            IF (THAM)  THEN
              CALL VCLR(DUMMY,1,36)
              CALL DFTB_SKSPA(NDER,KI,KJ,KK,KL,RI,VEC,HJI(9),HJI1(9),
     *                        HJI2(9),HJI3(9),HJI4(9),DUMMY(1,2:4),
     *                        AU(2:4,1))
            END IF
            IF (TOVER) THEN
              CALL VCLR(DUMMY,1,36)
              CALL DFTB_SKSPA(NDER,KI,KJ,KK,KL,RI,VEC,SJI(9),SJI1(9),
     *                        SJI2(9),SJI3(9),SJI4(9),DUMMY(1,2:4),
     *                        BU(2:4,1))
            END IF
          END IF
        ELSE ! (I,J) = (s,p) or (p,s)
          IF (MAXANGJ.GE.2) THEN ! (I,J) = (s,p)
            IF (THAM)  CALL DFTB_SKSPA(NDER,KI,KJ,KK,KL,RI,VEC,HIJ(9),
     *                            HIJ1(9),HIJ2(9),HIJ3(9),HIJ4(9),
     *                            AU(1,2:4),AU(2:4,1))
            IF (TOVER) CALL DFTB_SKSPA(NDER,KI,KJ,KK,KL,RI,VEC,SIJ(9),
     *                            SIJ1(9),SIJ2(9),SIJ3(9),SIJ3(9),
     *                            BU(1,2:4),BU(2:4,1))
          ELSE ! (I,J) = (p,s)
            IF (THAM)  THEN
              CALL VCLR(DUMMY,1,36)
              CALL DFTB_SKSPA(NDER,KI,KJ,KK,KL,RI,VEC,HJI(9),HJI1(9),
     *                        HJI2(9),HJI3(9),HJI4(9),DUMMY(1,2:4),
     *                        AU(2:4,1))
            END IF
            IF (TOVER) THEN
              CALL VCLR(DUMMY,1,36)
              CALL DFTB_SKSPA(NDER,KI,KJ,KK,KL,RI,VEC,SJI(9),SJI1(9),
     *                        SJI2(9),SJI3(9),SJI4(9),DUMMY(1,2:4),
     *                        BU(2:4,1))
            END IF
          END IF
        END IF

        IF (MAXMAX.LE.2) RETURN
C
C d Interaction
C
        IF (MINMAX.EQ.3) THEN
          IF (THAM) THEN
            CALL DFTB_SKDDA(NDER,KI,KJ,KK,KL,RI,VEC,HIJ(1:3),HIJ1(1:3),
     *           HIJ2(1:3),HIJ3(1:3),HIJ4(1:3),AU(5:9,5:9))
            CALL DFTB_SKPDA(NDER,KI,KJ,KK,KL,RI,VEC,HIJ(4:5),HIJ1(4:5),
     *           HIJ2(4:5),HIJ3(4:5),HIJ4(4:5),AU(2:4,5:9),AU(5:9,2:4))
            CALL DFTB_SKSDA(NDER,KI,KJ,KK,KL,RI,VEC,HIJ(8),HIJ1(8),
     *           HIJ2(8),HIJ3(8),HIJ4(8),AU(1,5:9),AU(5:9,1))
          END IF
          IF (TOVER) THEN
            CALL DFTB_SKDDA(NDER,KI,KJ,KK,KL,RI,VEC,SIJ(1:3),SIJ1(1:3),
     *           SIJ2(1:3),SIJ3(1:3),SIJ4(1:3),BU(5:9,5:9))
            CALL DFTB_SKPDA(NDER,KI,KJ,KK,KL,RI,VEC,SIJ(4:5),SIJ1(4:5),
     *           SIJ2(4:5),SIJ3(4:5),SIJ4(4:5),BU(2:4,5:9),BU(5:9,2:4))
            CALL DFTB_SKSDA(NDER,KI,KJ,KK,KL,RI,VEC,SIJ(8),SIJ1(8),
     *           SIJ2(8),SIJ3(8),SIJ4(8),BU(1,5:9),BU(5:9,1))
          END IF
C
          IF (I.NE.J) THEN
            IF (THAM) THEN
              CALL VCLR(DUMMY(1,5),1,45)
              CALL DFTB_SKSDA(NDER,KI,KJ,KK,KL,RI,VEC,HJI(8),HJI1(8),
     *             HJI2(8),HJI3(8),HJI4(8),DUMMY(1,5:9),AU(5:9,1))
              CALL DFTB_SKPDA(NDER,KI,KJ,KK,KL,RI,VEC,HJI(4:5),
     *             HJI1(4:5),HJI2(4:5),HJI3(4:5),HJI4(4:5),
     *             DUMMY(2:4,5:9),AU(5:9,2:4))
            END IF
            IF (TOVER) THEN
              CALL VCLR(DUMMY(1,5),1,45)
              CALL DFTB_SKSDA(NDER,KI,KJ,KK,KL,RI,VEC,SJI(8),SJI1(8),
     *             SJI2(8),SJI3(8),SJI4(8),DUMMY(1,5:9),BU(5:9,1))
              CALL DFTB_SKPDA(NDER,KI,KJ,KK,KL,RI,VEC,SJI(4:5),
     *             SJI1(4:5),SJI2(4:5),SJI3(4:5),SJI4(4:5),
     *             DUMMY(2:4,5:9),BU(5:9,2:4))
            END IF
          END IF
        ELSE
          IF (MAXANGI.EQ.1) THEN
            IF (THAM)  CALL DFTB_SKSDA(NDER,KI,KJ,KK,KL,RI,VEC,HIJ(8),
     *                      HIJ1(8),HIJ2(8),HIJ3(8),HIJ4(8),AU(1,5:9),
     *                      AU(5:9,1))
            IF (TOVER) CALL DFTB_SKSDA(NDER,KI,KJ,KK,KL,RI,VEC,SIJ(8),
     *                      SIJ1(8),SIJ2(8),SIJ3(8),SIJ4(8),BU(1,5:9),
     *                      BU(5:9,1))
          ELSE
            IF (MAXANGI.EQ.2) THEN
              IF (THAM) THEN
                CALL DFTB_SKSDA(NDER,KI,KJ,KK,KL,RI,VEC,HIJ(8),HIJ1(8),
     *               HIJ2(8),HIJ3(8),HIJ4(8),AU(1,5:9),AU(5:9,1))
                CALL DFTB_SKPDA(NDER,KI,KJ,KK,KL,RI,VEC,HIJ(4:5),
     *               HIJ1(4:5),HIJ2(4:5),HIJ3(4:5),HIJ4(4:5),
     *               AU(2:4,5:9),AU(5:9,2:4))
              END IF
              IF (TOVER) THEN
                CALL DFTB_SKSDA(NDER,KI,KJ,KK,KL,RI,VEC,SIJ(8),SIJ1(8),
     *               SIJ2(8),SIJ3(8),SIJ4(8),BU(1,5:9),BU(5:9,1))
                CALL DFTB_SKPDA(NDER,KI,KJ,KK,KL,RI,VEC,SIJ(4:5),
     *               SIJ1(4:5),SIJ2(4:5),SIJ3(4:5),SIJ4(4:5),
     *               BU(2:4,5:9),BU(5:9,2:4))
              END IF
            ELSE
              IF (MAXANGJ.EQ.1) THEN
                IF (THAM)  THEN
                  CALL VCLR(DUMMY(1,5),1,45)
                  CALL DFTB_SKSDA(NDER,KI,KJ,KK,KL,RI,VEC,HJI(8),
     *                 HJI1(8),HJI2(8),HJI3(8),HJI4(8),DUMMY(1,5:9),
     *                 AU(5:9,1))
                END IF
                IF (TOVER) THEN
                  CALL VCLR(DUMMY(1,5),1,45)
                  CALL DFTB_SKSDA(NDER,KI,KJ,KK,KL,RI,VEC,SJI(8),
     *                 SJI1(8),SJI2(8),SJI3(8),SJI4(8),DUMMY(1,5:9),
     *                 BU(5:9,1))
                END IF
              ELSE
                IF (THAM) THEN
                  CALL VCLR(DUMMY(1,5),1,45)
                  CALL DFTB_SKSDA(NDER,KI,KJ,KK,KL,RI,VEC,HJI(8),
     *                 HJI1(8),HJI2(8),HJI3(8),HJI4(8),DUMMY(1,5:9),
     *                 AU(5:9,1))
                  CALL DFTB_SKPDA(NDER,KI,KJ,KK,KL,RI,VEC,HJI(4:5),
     *                 HJI1(4:5),HJI2(4:5),HJI3(4:5),HJI4(4:5),
     *                 DUMMY(2:4,5:9),AU(5:9,2:4))
                END IF
                IF (TOVER) THEN
                  CALL VCLR(DUMMY(1,5),1,45)
                  CALL DFTB_SKSDA(NDER,KI,KJ,KK,KL,RI,VEC,SJI(8),
     *                 SJI1(8),SJI2(8),SJI3(8),SJI4(8),DUMMY(1,5:9),
     *                 BU(5:9,1))
                  CALL DFTB_SKPDA(NDER,KI,KJ,KK,KL,RI,VEC,SJI(4:5),
     *                 SJI1(4:5),SJI2(4:5),SJI3(4:5),SJI4(4:5),
     *                 DUMMY(2:4,5:9),BU(5:9,2:4))
                END IF
              END IF
            END IF
          END IF
        END IF
C     ELSE IF (I.NE.J) THEN
C       RETURN
      END IF

      END SUBROUTINE DFTB_SKODEA
C
C-----------------------------------------------------------------------
C*MODULE DFTBSK    *DECK DFTB_SKPAR
C>
C>    @brief Interpolate and get Hamiltonian and overlap
C>
C>    @details Looking for Hamiltonian/overlap matrix parametes using
C>             Newton's appproximation
C>             This can calculate only zeroth order derivatives
C>
C>    @author Yoshio Nishimoto
C>    - Sep, 2014- Subroutine written
C>
C>           --- INPUT ---
C>    @param R2 Square distance
C>    @param INU Minimum index to be calculated
C>    @param SKINTV Size of grid
C>    @param SKDIM Number of grids
C>    @param SKTAB Tabulated Hamiltonian/overlap matrix elements
C>           --- OUTPUT ---
C>    @param DD Hamiltonian/overlap values in distance space
C
      SUBROUTINE DFTB_SKPAR(R2,INU,DD,SKINTV,SKDIM,SKTAB)
C
      IMPLICIT NONE
C
      INTEGER, PARAMETER :: MXTAB=1500
C
      DOUBLE PRECISION,INTENT(IN) :: R2
      INTEGER, INTENT(IN) :: INU
      DOUBLE PRECISION,INTENT(INOUT) :: DD(10)
      DOUBLE PRECISION,INTENT(IN) :: SKINTV
      DOUBLE PRECISION,INTENT(IN) :: SKDIM
      DOUBLE PRECISION,INTENT(IN) :: SKTAB(10,MXTAB)
C
      INTEGER :: MAXIND, MAXPAR, IN
      DOUBLE PRECISION :: R,X0,X1,X2,F0,F1,F2,XH,HL,GRDR,SKINTVI
      DOUBLE PRECISION :: DFTB_CUBICSPLINE, DFTB_SPLINE5TH
C
      INTEGER, PARAMETER :: IW = 6
C
C mxind = Maximaler Index bis zu dem der Spline weitergefuehrt wird
C MAXIND  = MAXIMUM INDEX FOR SPLINE ... ?
C MAXPAR = MAX NUMBER OF H AND S PARAMETERS
C
C SKINTV = around 0.02
C
      SKINTVI = 1.0D+00/SKINTV
      MAXIND = INT(SKDIM + (0.3D+00*SKINTVI - 1.0D+00))
      R = SQRT(R2)
      MAXPAR = INT(R*SKINTVI + 1.0D+00)
C
C If distance is almost zero
C
      DD = 0.0D+00
      IF (R2.LE.1.0D-08) THEN
        WRITE (IW,*) " TOO CLOSE DISTANCE ALTHOUGH NOT SAME ATOM"
        CALL ABRT
      ELSE
        IF (MAXPAR+2.GT.SKDIM) THEN
          IF (MAXPAR.LT.SKDIM) THEN
C
C free cubic spline
C
            X0 = (SKDIM-3.0D+00) * SKINTV
            X1 = X0 + SKINTV
            X2 = X1 + SKINTV
            XH = R - X1
            HL = X2 - X1
            DO IN = INU, 10
              F0 = SKTAB(IN,INT(SKDIM)-2)
              F1 = SKTAB(IN,INT(SKDIM)-1)
              F2 = SKTAB(IN,INT(SKDIM)  )
              DD(IN) = DFTB_CUBICSPLINE(F0,F1,F2,F0,F1,XH,HL,SKINTV)
            END DO
          ELSE
            IF ((MAXPAR.GE.SKDIM).AND.(MAXPAR.LT.MAXIND)) THEN
C
C 5th degree spline
C
              X0 = (SKDIM-3.0D+00) * SKINTV
              X1 = X0 + SKINTV
              X2 = X1 + SKINTV
              XH = R - (MAXIND-1.0D+00) * SKINTV
              DO IN = INU, 10
                F0 = SKTAB(IN,INT(SKDIM)-2)
                F1 = SKTAB(IN,INT(SKDIM)-1)
                F2 = SKTAB(IN,INT(SKDIM)  )
                DD(IN) = DFTB_SPLINE5TH(F0,F1,F2,X0,X1,X2,XH,SKINTV,
     *            MAXIND)
              END DO
            ELSE
C
C zero
C
              DO IN = INU, 10
                DD(IN) = 0.0D+00
              END DO
            END IF
          END IF
        ELSE
C
C  MAXPAR+2.GT.SKDIM :
C
          GRDR = (R-(MAXPAR-1.0D+00)*SKINTV)*SKINTVI
          DO IN = INU, 10
            F0 = SKTAB(IN,MAXPAR)
            F1 = SKTAB(IN,MAXPAR+1)
            F2 = SKTAB(IN,MAXPAR+2)
            DD(IN) = F0 + (F1-F0)*GRDR +
     &               (F2+F0-2.0D+00*F1)*GRDR*(GRDR-1.0D+00)/2.0D+00
          END DO
        END IF
      END IF

      RETURN

      END SUBROUTINE DFTB_SKPAR
C
C-----------------------------------------------------------------------
C*MODULE DFTBSK    *DECK DFTB_SKPARA
C>
C>    @brief Interpolate and get Hamiltonian and overlap
C>
C>    @details Looking for Hamiltonian/overlap matrix parametes using
C>             Newton's appproximation
C>             This is used for analytical Slater-Koster transformation,
C>             and this subroutine gives SK values differentiated with
C>             respect to distance
C>             This is a clone of DFTB_SKPAR
C>
C>    @author Yoshio Nishimoto
C>    - Sep, 2014- Subroutine written
C>    @date   Jun, 2018 - Yoshio Nishimoto
C>    - Third-order derivatives
C>
C>           --- INPUT ---
C>    @param R2     Square distance
C>    @param INU    Minimum index to be calculated
C>    @param SKINTV Size of grid
C>    @param SKDIM  Number of grids
C>    @param SKTAB  Tabulated Hamiltonian/overlap matrix elements
C>    @param NDER   Order of derivative to be done
C>           --- OUTPUT ---
C>    @param DD     Hamiltonian/overlap values in distance space
C>    @param DD1    First-order  der of H and S
C>    @param DD2    Second-order der of H and S
C>    @param DD3    Third-order  der of H and S
C>    @param DD4    Fourth-order der of H and S
C
      SUBROUTINE DFTB_SKPARA(R,INU,DD,DD1,DD2,DD3,DD4,SKGRID,SKDIM,
     *                       SKTAB,NDER,MAXIND,MAXPAR)

      IMPLICIT NONE

      INTEGER,PARAMETER :: MXTAB=1500

C     INTEGER,INTENT(IN) :: I, J
      INTEGER,INTENT(IN) :: NDER
      INTEGER,INTENT(IN) :: INU
      INTEGER,INTENT(IN) :: MAXIND, MAXPAR
      DOUBLE PRECISION,INTENT(INOUT) :: DD(10),DD1(10),DD2(10),DD3(10),
     *  DD4(10)
      DOUBLE PRECISION,INTENT(IN) :: R
      DOUBLE PRECISION,INTENT(IN) :: SKGRID(2)
      DOUBLE PRECISION,INTENT(IN) :: SKDIM
      DOUBLE PRECISION,INTENT(IN) :: SKTAB(10,MXTAB)
C
      INTEGER :: IN
      DOUBLE PRECISION :: X0,X1,X2,F0,F1,F2,XH,HL,GRDR,SKINTV,SKINTVI,
     *                    SKINTVI2,XH2,XH3,HLI,HLI2,HLI3,HL2,HLMX,
     *                    HLMXI,HLMXI2,HLMXI3,HLMXI4,HLMXI5
C
      SKINTV  = SKGRID(1)
      SKINTVI = SKGRID(2)
C
C If distance is almost zero
C
      DD = 0.0D+00
      IF (NDER.GE.1) DD1 = 0.0D+00
      IF (NDER.GE.2) DD2 = 0.0D+00
      IF (NDER.GE.3) DD3 = 0.0D+00
      IF (NDER.GE.4) DD4 = 0.0D+00
      IF (MAXPAR+2.GT.SKDIM) THEN
        SKINTVI2 = SKINTVI*SKINTVI
        IF (MAXPAR.LT.SKDIM) THEN
C
C free cubic spline
C
          X0   = (SKDIM-3.0D+00) * SKINTV
          X1   = X0 + SKINTV
          X2   = X1 + SKINTV
          XH   = R - X1
          XH2  = XH*XH
          XH3  = XH*XH2
          HL   = X2 - X1
          HLI  = 1.0D+00/HL
          HLI2 = HLI*HLI
          HLI3 = HLI*HLI2
          DO IN = INU, 10
            F0 = SKTAB(IN,INT(SKDIM)-2)
            F1 = SKTAB(IN,INT(SKDIM)-1)
            F2 = SKTAB(IN,INT(SKDIM)  )
            CALL DFTB_CSPL(NDER,DD(IN),DD1(IN),DD2(IN),DD3(IN),
     *                     DD4(IN),F0,F1,F2,X0,X1,XH,XH2,XH3,HLI,
     *                     HLI2,HLI3,SKINTVI,SKINTVI2)
C           write (*,*) in, dd(in),dd1(in)
          END DO
        ELSE
          IF ((MAXPAR.GE.SKDIM).AND.(MAXPAR.LT.MAXIND)) THEN
C
C 5th degree spline
C
            X0     = (SKDIM-3.0D+00) * SKINTV
            X1     = X0 + SKINTV
            X2     = X1 + SKINTV
            XH     = R - (MAXIND-1.0D+00) * SKINTV
            XH2    = XH*XH
            XH3    = XH*XH2
            HL     = X2-X1
            HL2    = HL*HL
            HLI    = 1.0D+00/HL
            HLI2   = HLI*HLI
            HLI3   = HLI*HLI2
            HLMX   = X2-DBLE(MAXIND-1)*SKINTV
            HLMXI  = 1.0D+00/HLMX
            HLMXI2 = HLMXI*HLMXI
            HLMXI3 = HLMXI*HLMXI2
            HLMXI4 = HLMXI*HLMXI3
            HLMXI5 = HLMXI*HLMXI4
            DO IN = INU, 10
              F0 = SKTAB(IN,INT(SKDIM)-2)
              F1 = SKTAB(IN,INT(SKDIM)-1)
              F2 = SKTAB(IN,INT(SKDIM)  )
              CALL DFTB_SPL5(NDER,DD(IN),DD1(IN),DD2(IN),DD3(IN),
     *                       DD4(IN),F0,F1,F2,X0,X1,XH,XH2,XH3,HL,HL2,
     *                       HLI,HLI2,HLI3,HLMXI,HLMXI2,HLMXI3,HLMXI4,
     *                       HLMXI5,SKINTVI,SKINTVI2)
C           write (*,*) in, dd(in),dd1(in)
            END DO
          ELSE
C
C zero
C
C           DO IN = INU, 10
C             DD(IN) = 0.0D+00
C           END DO
          END IF
        END IF
      ELSE
C
C  MAXPAR+2.GT.SKDIM :
C
C       GRDR = (R-(DBLE(MAXPAR)-1.0D+00)*SKINTV)*SKINTVI
        GRDR = (R-(MAXPAR-1.0D+00)*SKINTV)*SKINTVI
        DO IN = INU, 10
          F0 = SKTAB(IN,MAXPAR)
          F1 = SKTAB(IN,MAXPAR+1)
          F2 = SKTAB(IN,MAXPAR+2)
          IF (NDER.EQ.0) THEN
            DD(IN)  = F0 + (F1-F0)*GRDR +
     &               (F2+F0-2.0D+00*F1)*GRDR*(GRDR-1.0D+00)*0.5D+00
          ELSE IF (NDER.EQ.1) THEN
            DD(IN)  = F0 + (F1-F0)*GRDR +
     &               (F2+F0-2.0D+00*F1)*GRDR*(GRDR-1.0D+00)*0.5D+00
            DD1(IN) = -(F2-4.0D+00*F1+3.0D+00*F0)*0.5D+00
     *              + (F2+F0-2.0D+00*F1)*GRDR
            DD1(IN) = DD1(IN) * SKINTVI
          ELSE IF (NDER.EQ.2) THEN
            DD(IN)  = F0 + (F1-F0)*GRDR +
     &               (F2+F0-2.0D+00*F1)*GRDR*(GRDR-1.0D+00)*0.5D+00
            DD1(IN) = -(F2-4.0D+00*F1+3.0D+00*F0)*0.5D+00
     *              + (F2+F0-2.0D+00*F1)*GRDR
            DD1(IN) = DD1(IN) * SKINTVI
            DD2(IN) = (F2+F0-2.0D+00*F1)*SKINTVI*SKINTVI
          ELSE IF (NDER.GE.3) THEN
            DD(IN)  = F0 + (F1-F0)*GRDR +
     &               (F2+F0-2.0D+00*F1)*GRDR*(GRDR-1.0D+00)*0.5D+00
            DD1(IN) = -(F2-4.0D+00*F1+3.0D+00*F0)*0.5D+00
     *              + (F2+F0-2.0D+00*F1)*GRDR
            DD1(IN) = DD1(IN) * SKINTVI
            DD2(IN) = (F2+F0-2.0D+00*F1)*SKINTVI*SKINTVI
            DD3(IN) = 0.0D+00
            DD4(IN) = 0.0D+00
          END IF
C           write (*,*) in, dd(in),dd1(in)
        END DO
      END IF

      RETURN

      END SUBROUTINE DFTB_SKPARA
C
C-----------------------------------------------------------------------
C*MODULE DFTBSK    *DECK DFTB_CSPL
C>
C>    @brief Interpolation with cubicspline
C>
C>    @details Interpolate Hamiltonian/overlap elements using three
C>             adjuacent points with cubic spline, giving 0th order
C>             derivatives
C>
C>    @author Yoshio Nishimoto
C>    - Sep, 2014- Subroutine written
C>    @date   Jun, 2018 - Yoshio Nishimoto
C>    - Third- and fourth-order derivatives
C>
C>           --- INPUT ---
C>    @param F0,F1,F2 adjacent three points for interpolation
C>    @param X0       R of first point
C>    @param X1       R of second point
C>    @param XH       R-X1
C>    @param XH2      XH*XH
C>    @param XH3      XH*XH*XH
C>    @param HL       Size of grid
C>    @param HLI      H / LI
C>    @param HLI2     H / LI**2
C>    @param HLI3     H / LI**3
C>    @param DRI      Size of grid
C>    @param DRI2     Size of grid square
C>           --- OUTPUT ---
C>    @param DD0      0th-order value
C>    @param DD1      1st-order value
C>    @param DD2      2nd-order value
C>    @param DD3      3rd-order value
C>    @param DD4      4th-order value
C
      SUBROUTINE DFTB_CSPL(NDER,DD0,DD1,DD2,DD3,DD4,F0,F1,F2,X0,X1,XH,
     *  XH2,XH3,HLI,HLI2,HLI3,DRI,DRI2)
C
      IMPLICIT NONE
C
      INTEGER NDER
      DOUBLE PRECISION DD0,DD1,DD2,DD3,DD4,F0,F1,F2,X0,X1,XH,XH2,XH3,
     *                 HLI,HLI2,HLI3,DRI,DRI2
      DOUBLE PRECISION F1ABL,F2ABL,A,B,C,D
C
      F2ABL = (F2+F0-2.0D+00*F1)*DRI2
      F1ABL = (F1-F0)*DRI + 0.5D+00*F2ABL*(X1-X0)
      A = F1
      B = F1ABL
      C = F2ABL*0.5D+00
      D = (F2-A)*HLI3 - B*HLI2 - C*HLI
      DD0 = A + B*XH + C*XH2+ D*XH3
      IF (NDER.GE.1) DD1 = B + 2.0D+00*C*XH + 3.0D+00*D*XH2
      IF (NDER.GE.2) DD2 = 2.0D+00*C + 6.0D+00*D*XH
      IF (NDER.GE.3) DD3 = 6.0D+00*D
      IF (NDER.GE.4) DD4 = 0.0D+00
C
      END SUBROUTINE DFTB_CSPL
C
C-----------------------------------------------------------------------
C*MODULE DFTBSK    *DECK DFTB_CUBICSPLINE
C>
C>    @brief Interpolation with cubicspline
C>
C>    @details Interpolate Hamiltonian/overlap elements using three
C>             adjuacent points with cubic spline, giving 0th order
C>             derivatives
C>
C>    @author Yoshio Nishimoto
C>    - Sep, 2014- Subroutine written
C>
C>           --- INPUT ---
C>    @param F0,F1,F2 adjacent three points for interpolation
C>    @param XH (R-x0)?
C>    @param HL,DR Size of grid
C>           --- OUTPUT ---
C>    @param CUBICSPLINE Interpolated zeroth derivative value
C
      DOUBLE PRECISION FUNCTION DFTB_CUBICSPLINE(F0,F1,F2,X0,X1,XH,HL,
     *  DR)
C
      IMPLICIT NONE
C
      DOUBLE PRECISION F0,F1,F2,X0,X1,XH,HL,DR
      DOUBLE PRECISION F1ABL,F2ABL,A,B,C,D
C
      F2ABL = (F2+F0-2.0D+00*F1)/(DR*DR)
      F1ABL = (F1-F0)/DR + 0.5D+00*F2ABL*(X1-X0)
      A = F1
      B = F1ABL
      C = F2ABL/2.0D+00
      D = (F2-A)/(HL*HL*HL) - B/(HL*HL) - C/HL
      DFTB_CUBICSPLINE = A + B*XH + C*XH*XH + D*XH*XH*XH
C
      END FUNCTION DFTB_CUBICSPLINE
C
C-----------------------------------------------------------------------
C*MODULE DFTBSK    *DECK DFTB_SPL5
C>
C>    @brief Interpolation with 5th-order spline
C>
C>    @details Interpolate Hamiltonian/overlap elements using three
C>             adjuacent points with 5th-order spline, giving 0th--4th
C>             order derivatives
C>
C>    @author Yoshio Nishimoto
C>    - Sep, 2014- Subroutine written
C>    @date   Jun, 2016 - Yoshio Nishimoto
C>    - Bug fix
C>    @date   Feb, 2018 - Yoshio Nishimoto
C>    - Third- and fourth-order derivatives
C>
C>           --- INPUT ---
C>    @param NDER     Order of derivatives
C>    @param F0,F1,F2 Adjacent three points for interpolation
C>    @param X0       R of first point
C>    @param X1       R of second point
C>    @param XH       Distance between R and the last point
C>    @param XH2      XH*XH
C>    @param XH3      XH*XH*XH
C>    @param HL       Size of grid
C>    @param HL2      HL*HL
C>    @param HLI      H / LI
C>    @param HLI2     H / LI**2
C>    @param HLI3     H / LI**3
C>    @param HLMXI    Inverse distance between X2 and the last point
C>    @param HLMXI2   HLMXI*HLMXI
C>    @param HLMXI3   HLMXI*HLMXI*HLMXI
C>    @param HLMXI4   HLMXI*HLMXI*HLMXI*HLMXI
C>    @param HLMXI5   HLMXI*HLMXI*HLMXI*HLMXI*HLMXI
C>    @param DRI      Size of grid
C>    @param DRI2     Size of grid square
C>           --- OUTPUT ---
C>    @param DD0      0th-order derivatives
C>    @param DD1      1st-order derivatives
C>    @param DD2      2nd-order derivatives
C>    @param DD3      3rd-order derivatives
C>    @param DD4      4th-order derivatives
C
      SUBROUTINE DFTB_SPL5(NDER,DD0,DD1,DD2,DD3,DD4,F0,F1,F2,X0,X1,XH,
     *                     XH2,XH3,HL,HL2,HLI,HLI2,HLI3,HLMXI,HLMXI2,
     *                     HLMXI3,HLMXI4,HLMXI5,DRI,DRI2)
C
      IMPLICIT NONE
C
      INTEGER NDER
      DOUBLE PRECISION DD0,DD1,DD2,DD3,DD4,F0,F1,F2,X0,X1,XH,XH2,XH3,
     *                 HL,HL2,HLI,HLI2,HLI3,HLMXI,HLMXI2,HLMXI3,HLMXI4,
     *                 HLMXI5,DRI,DRI2
      DOUBLE PRECISION F1ABL,F2ABL,A,B,C,D,HSP,ISP,JSP
C
      F2ABL = (F2+F0-2.0D+00*F1)*DRI2
      F1ABL = (F1-F0)*DRI + 0.5D+00*F2ABL*(X1-X0)
      A = F1
      B = F1ABL
      C = F2ABL*0.5D+00
      D = (F2-A)*HLI3 - B*HLI2 - C*HLI
      F1ABL = B + 2.0D+00*C*HL + 3.0D+00*D*HL2
      F2ABL = 2.0D+00*C + 6.0D+00*D*HL
C
      HSP =  1.0D+01*F2   *HLMXI3
     *    -  4.0D+00*F1ABL*HLMXI2
     *    +  0.5D+00*F2ABL*HLMXI
      ISP = -1.5D+01*F2   *HLMXI4
     *    +  7.0D+00*F1ABL*HLMXI3
     *    -          F2ABL*HLMXI2
      JSP =  6.0D+00*F2   *HLMXI5
     *    -  3.0D+00*F1ABL*HLMXI4
     *    +  0.5D+00*F2ABL*HLMXI3
C
      DD0 = (HSP + ISP*XH + JSP*XH2) * XH3
      IF (NDER.GE.1)
     *   DD1 = (3.0D+00*HSP+4.0D+00*ISP*XH+5.0D+00*JSP*XH2)*XH2
      IF (NDER.GE.2)
     *   DD2 = (6.0D+00*HSP+1.2D+01*ISP*XH+2.0D+01*JSP*XH2)*XH
      IF (NDER.GE.3)
     *   DD3 = (6.0D+00*HSP+2.4D+01*ISP*XH+6.0D+01*JSP*XH2)
      IF (NDER.GE.4)
     *   DD4 = (            2.4D+01*ISP   +1.2D+02*JSP*XH )
C
      END SUBROUTINE DFTB_SPL5
C
C-----------------------------------------------------------------------
C*MODULE DFTBSK    *DECK DFTB_SPLINE5TH
C>
C>    @brief Interpolation with 5th-order spline
C>
C>    @details Interpolate Hamiltonian/overlap elements using three
C>             adjuacent points with 5th-order spline, giving 0th order
C>             derivatives
C>
C>    @author Yoshio Nishimoto
C>    - Sep, 2014- Subroutine written
C>
C>           --- INPUT ---
C>    @param F0,F1,F2 adjacent three points for interpolation
C>    @param X0,X1,X2 x of F0,F1,F2
C>    @param XH (R-x0)?
C>    @param DR Size of grid
C>    @param MXIND Upper limit of Slater-Koster files
C>           --- OUTPUT ---
C>    @param SPLINE5TH Interpolated zeroth derivative value
C
      DOUBLE PRECISION FUNCTION DFTB_SPLINE5TH(F0,F1,F2,X0,X1,X2,XH,DR,
     *  MXIND)
C
      IMPLICIT NONE
C
      DOUBLE PRECISION F0,F1,F2,X0,X1,X2,XH,DR
      INTEGER MXIND
      DOUBLE PRECISION HL,F1ABL,F2ABL,A,B,C,D,HSP,ISP,JSP
C
      F2ABL = (F2+F0-2.0D+00*F1)/(DR*DR)
      F1ABL = (F1-F0)/DR + 0.5D+00*F2ABL*(X1-X0)
      A = F1
      B = F1ABL
      C = F2ABL/2.0D+00
      HL = X2-X1
      D = (F2-A)/(HL*HL*HL) - B/(HL*HL) - C/HL
      F1ABL = B + 2.0D+00*C*HL + 3.0D+00*D*HL*HL
      F2ABL = 2.0D+00*C + 6.0D+00*D*HL
      HL = X2 - (MXIND-1)*DR

      HSP = 1.0D+01*F2/(HL*HL*HL) - 4.0D+00*F1ABL/(HL*HL)
     &      + F2ABL/(2.0D+00*HL)
      ISP = -1.5D+01*F2/(HL*HL*HL*HL) + 7.0D+00*F1ABL/(HL*HL*HL)
     &      - F2ABL/(HL*HL)
      JSP = 6.0D+00*F2/(HL*HL*HL*HL*HL) - 3.0D+00*F1ABL/(HL*HL*HL*HL)
     &      + F2ABL/(2.0D+00*HL*HL*HL)
      HL = XH*XH*XH
      DFTB_SPLINE5TH = (HSP + ISP*XH + JSP*XH*XH) * HL
C
      END FUNCTION DFTB_SPLINE5TH
C
C-----------------------------------------------------------------------
C*MODULE DFTBSK    *DECK DFTB_SKSS
C>
C>    @brief Slater-Koster transformation of s-s interaction
C>
C>    @details zeroth order Slater-Koster transformation of s-s
C>             interaction
C>
C>    @author Yoshio Nishimoto
C>    - Sep, 2014- Subroutine written
C>
C>           --- INPUT ---
C>    @param X Vector of A-B
C>    @param X2 Squared vector of A-B
C>    @param DD Distance-based Hamiltonian/overlap elements
C>           --- OUTPUT ---
C>    @param MAT Block Hamiltonian/overlap elements in Cartesian
C
      SUBROUTINE DFTB_SKSS(DD,MAT)
C
      IMPLICIT NONE
C
      DOUBLE PRECISION,INTENT(IN) :: DD(1)
      DOUBLE PRECISION,INTENT(INOUT) :: MAT(1,1)
C
      MAT(1,1) = DD(1)
C
      RETURN
C
      END SUBROUTINE DFTB_SKSS
C
C-----------------------------------------------------------------------
C*MODULE DFTBSK    *DECK DFTB_SKSP
C>
C>    @brief Slater-Koster transformation of s-p interaction
C>
C>    @details zeroth order Slater-Koster transformation of s-p
C>             interaction
C>
C>    @author Yoshio Nishimoto
C>    - Sep, 2014- Subroutine written
C>
C>           --- INPUT ---
C>    @param X Vector of A-B
C>    @param X2 Squared vector of A-B
C>    @param DD Distance-based Hamiltonian/overlap elements
C>           --- OUTPUT ---
C>    @param MAT1,MAT2 Block Hamiltonian/overlap elements in Cartesian
C
      SUBROUTINE DFTB_SKSP(X,DD,MAT1,MAT2)
C
      IMPLICIT NONE
C
      DOUBLE PRECISION,INTENT(IN) :: X(6)
      DOUBLE PRECISION,INTENT(IN) :: DD(1)
      DOUBLE PRECISION,INTENT(INOUT) :: MAT1(1,3), MAT2(3,1)
C
      INTEGER :: I
C
C     INDEX ::
C       1 -> X
C       2 -> Y
C       3 -> Z
C
C     THIS ORDER IS DIFFERENT FROM THAT OF DFTB+
C     IN THE CASE OF DFTB+, THE ORDER IS 1 -> Y, 2 -> Z, 3 -> X
C
      DO I = 1, 3
        MAT1(1,I) =  X(I) * DD(1)
        MAT2(I,1) = -MAT1(1,I)
      END DO
C
      RETURN
C
      END SUBROUTINE DFTB_SKSP
C
C-----------------------------------------------------------------------
C*MODULE DFTBSK    *DECK DFTB_SKSD
C>
C>    @brief Slater-Koster transformation of s-d interaction
C>
C>    @details zeroth order Slater-Koster transformation of s-d
C>             interaction
C>
C>    @author Yoshio Nishimoto
C>    - Sep, 2014- Subroutine written
C>
C>           --- INPUT ---
C>    @param X Vector of A-B
C>    @param X2 Squared vector of A-B
C>    @param DD Distance-based Hamiltonian/overlap elements
C>           --- OUTPUT ---
C>    @param MAT1,MAT2 Block Hamiltonian/overlap elements in Cartesian
C
      SUBROUTINE DFTB_SKSD(X,X2,DD,MAT1,MAT2)
C
      IMPLICIT NONE
C
      DOUBLE PRECISION,INTENT(IN) :: X(6),X2(6)
      DOUBLE PRECISION,INTENT(IN) :: DD(1)
      DOUBLE PRECISION,INTENT(INOUT) :: MAT1(1,5), MAT2(5,1)
C
      INTEGER :: I
      DOUBLE PRECISION :: R3, ALPHA, BETA
C
C     R3    : 1.732050808...
C     ALPHA : X^2 + Y^2
C     BETA  : X^2 - Y^2
C
      R3 = SQRT(3.0D+00)
      ALPHA = X2(1) + X2(2)
      BETA  = X2(1) - X2(2)
C
C     S-D(XY) = SQRT(3)*XY*sigma
C     S-D(YZ) = SQRT(3)*YZ*sigma
C     S-D(ZX) = SQRT(3)*ZX*sigma
C     S-D(X^2-Y^2) = 0.5*SQRT(3)*(X^2-Y^2) = 0.5*SQRT(3)*BETA
C     S-D(3Z^2-R^2) = -0.5X^2-0.5Y^2+Z^2 = Z^2-0.5*ALPHA
C
C     INDEX::
C       1 -> XY
C       2 -> YZ
C       3 -> ZX
C       4 -> X^2 - Y^2
C       5 -> 3*Z^2 - R^2
C
      DO I = 1, 3
        MAT1(1,I) = R3*X(I)*X(I+1)
      END DO
      MAT1(1,4) = 0.5D+00*R3*BETA
      MAT1(1,5) = X2(3) - 0.5D+00*ALPHA
C
      DO I = 1, 5
        MAT1(1,I) = MAT1(1,I) * DD(1)
        MAT2(I,1) = MAT1(1,I)
      END DO
C
      RETURN
C
      END SUBROUTINE DFTB_SKSD
C
C-----------------------------------------------------------------------
C*MODULE DFTBSK    *DECK DFTB_SKPP
C>
C>    @brief Slater-Koster transformation of p-p interaction
C>
C>    @details zeroth order Slater-Koster transformation of p-p
C>             interaction
C>
C>    @author Yoshio Nishimoto
C>    - Sep, 2014- Subroutine written
C>
C>           --- INPUT ---
C>    @param X Vector of A-B
C>    @param X2 Squared vector of A-B
C>    @param DD Distance-based Hamiltonian/overlap elements
C>           --- OUTPUT ---
C>    @param MAT1 Block Hamiltonian/overlap elements in Cartesian
C
      SUBROUTINE DFTB_SKPP(X,X2,DD,MAT)
C
      IMPLICIT NONE
C
      DOUBLE PRECISION,INTENT(IN) :: X(6), X2(6)
      DOUBLE PRECISION,INTENT(IN) :: DD(2)
      DOUBLE PRECISION,INTENT(INOUT) :: MAT(3,3)
C
C TM : TRANSFORMATION ...
C PROD : PRODUCT
C        PROD(1) = X * X
C        PROD(2) = Y * Y
C        PROD(3) = Z * Z
C        PROD(4) = X * Y
C        PROD(5) = Y * Z
C        PROD(6) = Z * X
C
      DOUBLE PRECISION :: TM(6), PROD(6)
      INTEGER :: I, J, K, L
C
      DO I = 1, 3
        PROD(I) = X2(I)
        PROD(I+3) = X(I) * X(I+1)
      END DO
C
C DD(1) : SIGMA INTERACTION
C DD(2) : PI    INTERACTION
C
      DO I = 1, 3
        TM(I)   = PROD(I) * DD(1) + (1.0D+00-PROD(I)) * DD(2)
        TM(I+3) = PROD(I+3) * (DD(1)-DD(2))
      END DO
C
C PX-PX = X^2*sigma + (1-X^2)*pi -> TM(1)
C PY-PY = Y^2*sigma + (1-Y^2)*pi -> TM(2)
C PZ-PZ = Z^2*sigma + (1-Y^2)*pi -> TM(3)
C PX-PY = X*Y*(sigma-pi) -> TM(4) ... and so on
C
      DO I = 1, 3
        DO J = 1, I
          K = I - J
          L = 3*K - (K*(K-1))/2 + J
          MAT(I,J) = TM(L)
          MAT(J,I) = TM(L)
        END DO
      END DO
C
      RETURN
C
      END SUBROUTINE DFTB_SKPP
C
C-----------------------------------------------------------------------
C*MODULE DFTBSK    *DECK DFTB_SKPD
C>
C>    @brief Slater-Koster transformation of p-d interaction
C>
C>    @details zeroth order Slater-Koster transformation of p-d
C>             interaction
C>
C>    @author Yoshio Nishimoto
C>    - Sep, 2014- Subroutine written
C>
C>           --- INPUT ---
C>    @param X Vector of A-B
C>    @param X2 Squared vector of A-B
C>    @param DD Distance-based Hamiltonian/overlap elements
C>           --- OUTPUT ---
C>    @param MAT1,MAT2 Block Hamiltonian/overlap elements in Cartesian
C
      SUBROUTINE DFTB_SKPD(X,X2,DD,MAT1,MAT2)
C
      IMPLICIT NONE
C
      DOUBLE PRECISION,INTENT(IN) :: X(6),X2(6)
      DOUBLE PRECISION,INTENT(IN) :: DD(2)
      DOUBLE PRECISION,INTENT(INOUT) :: MAT1(3,5), MAT2(5,3)
C
      INTEGER :: I, J, N
      DOUBLE PRECISION :: R3, ALPHA, BETA, XYZ, Z_ALP, BETA2
      DOUBLE PRECISION :: TM(13,2)
C
C     XYZ   = X * Y *Z
C     Z_ALP = X^2 - ALPHA
C     BETA2 = SQRT(3)*0.5*BETA
C
      R3    = SQRT(3.0D+00)
      ALPHA = X2(1) + X2(2)
      BETA  = X2(1) - X2(2)
      XYZ   = X(1) * X(2) * X(3)
      Z_ALP = X2(3) - 0.5D+00*ALPHA
      BETA2 = R3 * 0.5D+00 * BETA
C
C     TM(1 ,1) :: P(X)-D(XY) SIGMA INTERACTION
C     TM(1 ,2) :: P(X)-D(XY) PI    INTERACTION
C     TM(2 ,1) :: P(Y)-D(YZ) SIGMA INTERACTION
C                   .
C                   .
C                   .
C     TM(4 ,1) :: P(X)-D(YZ)
C                 P(Y)-D(ZX)
C                 P(Z)-D(XY) SIGMA INTERACTION
C     TM(4 ,2) :: P(X)-D(YZ)
C                 P(Y)-D(ZX)
C                 P(Z)-D(XY) PI    INTERACTION
C     TM(5 ,1) :: P(X)-D(ZX) SIGMA INTERACTION
C     TM(5 ,2) :: P(X)-D(ZX) PI    INTERACTION
C     TM(6 ,1) :: P(Y)-D(XY) SIGMA INTERACTION
C                   .
C                   .
C                   .
C     TM(8 ,1) :: P(X)-D(X^2-Y^2) SIGMA INTERACTION
C     TM(8 ,2) :: P(X)-D(X^2-Y^2) PI    INTERACTION
C     TM(9 ,1) :: P(Y)-D(X^2-Y^2) SIGMA INTERACTION
C                   .
C                   .
C                   .
C     TM(11,1) :: P(X)-D(3Z^2-R^2) SIGMA INTERACTION
C     TM(11,2) :: P(X)-D(3Z^2-R^2) PI    INTERACTION
C
      DO I = 1, 3
        TM(I   ,1) = R3 * X2(I) * X(I+1)
        TM(I+4 ,1) = R3 * X2(I) * X(I+2)
        TM(I   ,2) = X(I+1) * (1.0D+00-2.0D+00*X2(I))
        TM(I+4 ,2) = X(I+2) * (1.0D+00-2.0D+00*X2(I))
        TM(I+7 ,1) = BETA2 * X(I)
        TM(I+10,1) = X(I) * Z_ALP
      END DO
      TM( 4,1) =  R3 * XYZ
      TM( 4,2) = -2.0D+00 * XYZ
      TM( 8,2) =  X(1) * (1.0D+00-BETA)
      TM( 9,2) = -X(2) * (1.0D+00+BETA)
      TM(10,2) = -X(3) * BETA
      TM(11,2) = -R3 * X(1) * X2(3)
      TM(12,2) = -R3 * X(2) * X2(3)
      TM(13,2) =  R3 * X(3) * ALPHA
C
C     MAT1(1,1) : P(X)-D(XY)
C     MAT1(1,2) : P(X)-D(YZ)
C     MAT1(1,3) : P(X)-D(ZX)
C     MAT1(1,4) : P(X)-D(X^2-Y^2)
C     MAT1(1,5) : P(X)-D(3Z^2-R^2)
C
      MAT1(1,1) = TM( 1,1)*DD(1) + TM( 1,2)*DD(2)
      MAT1(2,1) = TM( 6,1)*DD(1) + TM( 6,2)*DD(2)
      MAT1(3,1) = TM( 4,1)*DD(1) + TM( 4,2)*DD(2)
      MAT1(1,2) = MAT1(3,1)
      MAT1(2,2) = TM( 2,1)*DD(1) + TM( 2,2)*DD(2)
      MAT1(3,2) = TM( 7,1)*DD(1) + TM( 7,2)*DD(2)
      MAT1(1,3) = TM( 5,1)*DD(1) + TM( 5,2)*DD(2)
      MAT1(2,3) = MAT1(1,2)
      MAT1(3,3) = TM( 3,1)*DD(1) + TM( 3,2)*DD(2)
      DO I = 1, 2
        DO J = 1, 3
          N = J+7 + (I-1)*3
          MAT1(J,I+3) = TM(N,1)*DD(1) + TM(N,2)*DD(2)
        END DO
      END DO
C
      DO I = 1, 5
        DO J = 1, 3
          MAT2(I,J) = -MAT1(J,I)
        END DO
      END DO
C
      RETURN
C
      END SUBROUTINE DFTB_SKPD
C
C-----------------------------------------------------------------------
C*MODULE DFTBSK    *DECK DFTB_SKDD
C>
C>    @brief Slater-Koster transformation of d-d interaction
C>
C>    @details zeroth order Slater-Koster transformation of d-d
C>             interaction
C>
C>    @author Yoshio Nishimoto
C>    - Sep, 2014- Subroutine written
C>    @data Nov, 2014- Yoshio Nishimoto
C>    - Bug fix of transformation
C>
C>           --- INPUT ---
C>    @param X Vector of A-B
C>    @param X2 Squared vector of A-B
C>    @param DD Distance-based Hamiltonian/overlap elements
C>           --- OUTPUT ---
C>    @param MAT1,MAT2 Block Hamiltonian/overlap elements in Cartesian
C
      SUBROUTINE DFTB_SKDD(X,X2,DD,MAT)
C
      IMPLICIT NONE
C
      DOUBLE PRECISION,INTENT(IN) :: X(6),X2(6)
      DOUBLE PRECISION,INTENT(IN) :: DD(3)
      DOUBLE PRECISION,INTENT(INOUT) :: MAT(5,5)
C
      INTEGER :: I, J
      DOUBLE PRECISION :: ALPHA, BETA, XYZ, BETA2,
     &                    VAL, VAL1, VAL2, VAL3, VAL4
      DOUBLE PRECISION :: TM(15,3)
      DOUBLE PRECISION R3
C
      R3    = SQRT(3.0D+00)
      XYZ   = X(1) * X(2) * X(3)
C
C     TM(1,1) :: D(XY)-D(XY) SIGMA INTERACTION
C     TM(1,2) :: D(XY)-D(XY) PI    INTERACTION
C     TM(1,3) :: D(XY)-D(XY) DELTA INTERACTION
C     TM(2,1) :: D(YZ)-D(YZ) SIGMA INTERACTION
C                   .
C                   .
C                   .
C     TM(4,1) :: D(XY)-D(YZ) SIGMA INTERACTION
C     TM(4,2) :: D(XY)-D(YZ) PI    INTERACTION
C     TM(4,3) :: D(XY)-D(YZ) DELTA INTERACTION
C     TM(5,1) :: D(YZ)-D(ZX) SIGMA INTERACTION
C                   .
C                   .
C                   .
C     TM(7,1) :: D(XY)-D(ZX) SIGMA INTERACTION
C     TM(7,2) :: D(XY)-D(ZX) PI    INTERACTION
C     TM(7,3) :: D(XY)-D(ZX) DELTA INTERACTION
C     TM(8,1) :: D(YZ)-D(XY) SIGMA INTERACTION
C
      DO I = 1, 3
        ALPHA = X2(I) + X2(I+1)
        VAL = X2(I) * X2(I+1)
        TM(I,1) = 3.0D+00 * VAL
        TM(I,2) = ALPHA - 4.0D+00*VAL
        TM(I,3) = X2(I+2) + VAL
        VAL1 = XYZ * X(I+1)  !! X * Y^2 * Z
        VAL2 = X(I) * X(I+2) !! X * Z
        TM(I+3,1) = 3.0D+00 * VAL1
        TM(I+3,2) = VAL2 - 4.0D+00*VAL1
        TM(I+3,3) = VAL1 - VAL2
      END DO
      ALPHA = X2(1) + X2(2)
      BETA  = X2(1) - X2(2)
      BETA2 = R3 * 0.5D+00 * BETA
C     TM(7,1) = D(XY)-D(X^2-Y^2) SIGMA INTERACTION
C     TM(7,2) = D(XY)-D(X^2-Y^2) PI    INTERACTION
C     TM(7,3) = D(XY)-D(X^2-Y^2) DELTA INTERACTION
C     TM(8,1) = D(YZ)-D(X^2-Y^2) SIGMA INTERACTION
C     TM(9,1) = D(ZX)-D(X^2-Y^2) SIGMA INTERACTION
C
C     TM(10,1) = D(XY)-D(3Z^2-R^2) SIGMA INTERACTION
C     TM(11,1) = D(YZ)-D(3Z^2-R^2) SIGMA INTERACTION
C     TM(12,1) = D(ZX)-D(3Z^2-R^2) SIGMA INTERACTION
C
      VAL1 = X(1) * X(2)   !! X * Y
      VAL2 = VAL1 * BETA   !! X * Y * BETA
      VAL3 = VAL1 * X2(3)  !! X * Y * Z^2
      VAL4 = VAL1 * ALPHA  !! X * Y * ALPHA
      TM(7,1) =  1.5D+00 * VAL2
      TM(7,2) = -2.0D+00 * VAL2
      TM(7,3) =  0.5D+00 * VAL2
      TM(10,1) =  R3 * (VAL3 - 0.5D+00 * VAL4)
      TM(10,2) = -R3 * 2.0D+00 * VAL3
      TM(10,3) =  R3 * 0.5D+00 * (VAL1 + VAL3)
      VAL1 = X(2) * X(3)
      VAL2 = VAL1 * BETA
      VAL3 = VAL1 * X2(3)
      VAL4 = VAL1 * ALPHA
      TM(8,1) =  1.5D+00 * VAL2
      TM(8,2) = -VAL1 - 2.0D+00 * VAL2
      TM(8,3) =  VAL1 + 0.5D+00 * VAL2
      TM(11,1) =  R3 * (VAL3 - 0.5D+00 * VAL4)
      TM(11,2) =  R3 * (VAL4 - VAL3)
      TM(11,3) = -R3 * 0.5D+00 * VAL4
      VAL1 = X(3) * X(4)
      VAL2 = VAL1 * BETA
      VAL3 = VAL1 * X2(3)
      VAL4 = VAL1 * ALPHA
      TM(9,1) =  1.5D+00 * VAL2
      TM(9,2) =  VAL1 - 2.0D+00 * VAL2
      TM(9,3) = -VAL1 + 0.5D+00 * VAL2
      TM(12,1) =  R3 * (VAL3 - 0.5D+00 * VAL4)
      TM(12,2) =  R3 * (VAL4 - VAL3)
      TM(12,3) = -R3 * 0.5D+00 * VAL4
C
      BETA2 = BETA * BETA
C     TM(13,:) :: D(X^2-Y^2) -D(X^2-Y^2)  INTERACTION
C     TM(14,:) :: D(X^2-Y^2) -D(3Z^2-R^2) INTERACTION
C     TM(15,:) :: D(3Z^2-R^2)-D(3Z^2-R^2) INTERACTION
      TM(13,1) = 0.75D+00 * BETA2
      TM(13,2) = ALPHA - BETA2
      TM(13,3) = X2(3)+ 0.25D+00 * BETA2
      VAL = BETA * X2(3)
      TM(14,1) =  R3 * 0.5D+00 * (VAL - 0.5D+00 * ALPHA * BETA)
      TM(14,2) = -R3 * VAL
      TM(14,3) =  R3 * 0.25D+00 * (BETA + VAL)
      TM(15,1) = (X2(3)-0.5D+00 * ALPHA)**2
      TM(15,2) = 3.0D+00 * X2(3) * ALPHA
      TM(15,3) = 0.75D+00 * ALPHA * ALPHA
C
C     D(XY), D(YZ), D(ZX), D(X^2-Y^2), D(3Z^2-R^2)
C
C     D(XY)-D(XY), D(YZ)-D(YZ), D(ZX)-D(ZX) INTERACTION
      DO I = 1, 3
        MAT(I,I) = TM(I,1)*DD(1) + TM(I,2)*DD(2) + TM(I,3)*DD(3)
      END DO
C     TM(4,:) :: D(XY)-D(YZ) -> MAT(2,1)
C     TM(5,:) :: D(YZ)-D(ZX) -> MAT(3,2)
C     TM(6,:) :: D(ZX)-D(XY) -> MAT(3,1)
      MAT(2,1) = TM(4,1)*DD(1) + TM(4,2)*DD(2) + TM(4,3)*DD(3)
      MAT(3,1) = TM(6,1)*DD(1) + TM(6,2)*DD(2) + TM(6,3)*DD(3)
      MAT(3,2) = TM(5,1)*DD(1) + TM(5,2)*DD(2) + TM(5,3)*DD(3)
      DO I = 1, 3
C     D(X^2-Y^2) -[D(XY),D(YZ),D(ZX)] INTERACTION (7, 8, 9)
        MAT(4,I) = TM(I+6,1)*DD(1) + TM(I+6,2)*DD(2) + TM(I+6,3)*DD(3)
C     D(3Z^2-R^2)-[D(XY),D(YZ),D(ZX)] INTERACTION (10, 11, 12)
        MAT(5,I) = TM(I+9,1)*DD(1) + TM(I+9,2)*DD(2) + TM(I+9,3)*DD(3)
      END DO
C     D(X^2-Y^2) - D(X^2-Y^2) INTERACTION
      MAT(4,4) = TM(13,1)*DD(1) + TM(13,2)*DD(2) + TM(13,3)*DD(3)
C     D(X^2-Y^2) - D(3Z^2-R^2) INTERACTION
      MAT(5,4) = TM(14,1)*DD(1) + TM(14,2)*DD(2) + TM(14,3)*DD(3)
C     D(3Z^2-R^2) - D(3Z^2-R^2) INTERACTION
      MAT(5,5) = TM(15,1)*DD(1) + TM(15,2)*DD(2) + TM(15,3)*DD(3)
      DO I = 1, 5
        DO J = 1, I-1
          MAT(J,I) = MAT(I,J)
        END DO
      END DO
C
      RETURN
C
      END SUBROUTINE DFTB_SKDD
C
C-----------------------------------------------------------------------
C*MODULE DFTBSK    *DECK DFTB_AD_PREP
C>
C>    @brief Prepare automatic differentiation (forward accumulation)
C>
C>    @details Prepare automatic-differentiation (AD)-like Slater-
C>             Koster transformation derivative calculation.
C>
C>    @author Yoshio Nishimoto
C>    - Sep, 2014- Subroutine written
C>    @date  Feb, 2018 - Yoshio Nishimoto
C>    - Third-order derivatives
C>
C>           --- INPUT ---
C>    @param NDER  Order of derivativesd
C>    @param KI    First axis to be differentiated
C>    @param KJ    Second axis to be differentiated
C>    @param KK    Third axis to be differentiated
C>    @param KL    Fourth axis to be differentiated
C>    @param DISTI Inverse of distance
C>    @param VEC   Vector of A-B
C>    @param F0    0th-order derivative of H or S wrt distance
C>    @param F1    1st-order derivative of H or S wrt distance
C>    @param F2    2nd-order derivative of H or S wrt distance
C>    @param F3    3rd-order derivative of H or S wrt distance
C>    @param F4    4th-order derivative of H or S wrt distance
C>           --- OUTPUT ---
C>    @param MAT1  Derivative seeds of Cartesian
C>    @param MAT2  Derivative seeds of Hamiltonian/overlap elements
C
      SUBROUTINE DFTB_AD_PREP(NDER,KI,KJ,KK,KL,DISTI,VEC,F0,F1,F2,F3,F4,
     *                        MAT1,MAT2)
C
      IMPLICIT NONE
C
      INTEGER, INTENT(IN) :: NDER,KI,KJ,KK,KL
      DOUBLE PRECISION, INTENT(IN) :: VEC(3),DISTI,F0,F1,F2,F3,F4
      DOUBLE PRECISION, INTENT(OUT) :: MAT1(2**NDER,3),MAT2(*)
C
      INTEGER :: I
      DOUBLE PRECISION :: DK_DL,DK_DM,DL_DM,DX_DK,DX_DL,DX_DM
      DOUBLE PRECISION, PARAMETER :: ZERO=0.0D+00,ONE=1.0D+00,
     *  THREE=3.0D+00
C
C     PREPARE AUTOMATIC-DFFIERENTIATION(AD)-LIKE SLATER-KOSTER
C     TRANSFORMATION DERIVATIVE CALCULATION.  THESE TRANSFORMATIONS
C     USE VERY BASIC AD.
C
C     MAT1( 1,1) =               X/R
C     MAT1( 2,1) = (d /dk      ) X/R
C     MAT1( 3,1) = (d /  dl    ) X/R
C     MAT1( 4,1) = (d2/dkdl    ) X/R
C     MAT1( 5,1) = (d2/    dm  ) X/R
C     MAT1( 6,1) = (d2/dk  dm  ) X/R
C     MAT1( 7,1) = (d2/  dldm  ) X/R
C     MAT1( 8,1) = (d3/dkdldm  ) X/R
C     MAT1( 9,1) = (d /      dn) X/R
C     MAT1(10,1) = (d2/dk    dn) X/R
C     MAT1(11,1) = (d2/  dl  dn) X/R
C     MAT1(12,1) = (d2/    dmdn) X/R
C     MAT1(13,1) = (d3/dkdl  dn) X/R
C     MAT1(14,1) = (d3/dk  dmdn) X/R
C     MAT1(15,1) = (d3/  dldmdn) X/R
C     MAT1(16,1) = (d4/dkdldmdn) X/R
C     MAT1( 1,2) =               Y/R
C     MAT1( 2,2) = (d /dk      ) Y/R
C     MAT1( 3,2) = (d /  dl    ) Y/R
C       ......
C
C     F0 =          C
C     F1 = (d /dR ) C
C     F2 = (d2/dR2) C
C     F3 = (d3/dR3) C
C     F4 = (d4/dR4) C
C     MAT2( 1)   =               C
C     MAT2( 2)   = (d /dk      ) C
C     MAT2( 3)   = (d /  dl    ) C
C     MAT2( 4)   = (d2/dkdl    ) C
C     MAT2( 5)   = (d /    dm  ) C
C     MAT2( 6)   = (d2/dk  dm  ) C
C     MAT2( 7)   = (d2/  dldm  ) C
C     MAT2( 8)   = (d3/dkdldm  ) C
C     MAT2( 9)   = (d /      dn) C
C     MAT2(10)   = (d2/dk    dn) C
C     MAT2(11)   = (d2/  dl  dn) C
C     MAT2(12)   = (d2/    dmdn) C
C     MAT2(13)   = (d3/dkdl  dn) C
C     MAT2(14)   = (d3/dk  dmdn) C
C     MAT2(15)   = (d3/  dldmdn) C
C     MAT2(16)   = (d4/dkdldmdn) C
C
      DK_DL = ZERO
      DK_DM = ZERO
      DL_DM = ZERO
      IF (KI.EQ.KJ) DK_DL = ONE
      IF (KI.EQ.KK) DK_DM = ONE
      IF (KJ.EQ.KK) DL_DM = ONE
C
      DO I = 1, 3
        MAT1(1,I) = VEC(I)
        IF (NDER.EQ.0) CYCLE
C
        DX_DK = ZERO
        IF (I.EQ.KI) DX_DK = ONE
        MAT1(2,I) = DISTI*(DX_DK-VEC(I)*VEC(KI))
        IF (NDER.EQ.1) CYCLE
C
        DX_DL = ZERO
        IF (I.EQ.KJ) DX_DL = ONE
        MAT1(3,I) = DISTI*(DX_DL-VEC(I)*VEC(KJ))
        MAT1(4,I) = -DISTI*DISTI*(VEC(KI)*DX_DL+VEC(KJ)*DX_DK
     *    +VEC(I)*DK_DL-THREE*VEC(I)*VEC(KI)*VEC(KJ))
        IF (NDER.EQ.2) CYCLE
C
        DX_DM = ZERO
        IF (I.EQ.KK) DX_DM = ONE
        MAT1(5,I) = DISTI*(DX_DM-VEC(I)*VEC(KK))
        MAT1(6,I) = -DISTI*DISTI*(VEC(KI)*DX_DM+VEC(KK)*DX_DK
     *    +VEC(I)*DK_DM-THREE*VEC(I)*VEC(KI)*VEC(KK))
        MAT1(7,I) = -DISTI*DISTI*(VEC(KJ)*DX_DM+VEC(KK)*DX_DL
     *    +VEC(I)*DL_DM-THREE*VEC(I)*VEC(KJ)*VEC(KK))
        MAT1(8,I) = -DISTI*DISTI*DISTI*(
     *    (DX_DK*DL_DM + DX_DL*DK_DM + DX_DM*DK_DL)
     *    -THREE*(VEC(I) *VEC(KI)*DL_DM
     *          + VEC(I) *VEC(KJ)*DK_DM
     *          + VEC(I) *VEC(KK)*DK_DL
     *          + VEC(KI)*VEC(KJ)*DX_DM
     *          + VEC(KJ)*VEC(KK)*DX_DK
     *          + VEC(KK)*VEC(KI)*DX_DL)
     *    + 1.5D+01*VEC(I)*VEC(KI)*VEC(KJ)*VEC(KK))
      END DO
C
      MAT2(1) = F0
      IF (NDER.EQ.0) RETURN
C
      MAT2(2) = VEC(KI)*F1
      IF (NDER.EQ.1) RETURN
C
      MAT2(3) = VEC(KJ)*F1
      MAT2(4) = DISTI*(DK_DL-VEC(KI)*VEC(KJ))*F1
     *  +VEC(KI)*VEC(KJ)*F2
      IF (NDER.EQ.2) RETURN
C
      MAT2(5) = VEC(KK)*F1
      MAT2(6) = DISTI*(DK_DM-VEC(KI)*VEC(KK))*F1
     *  +VEC(KI)*VEC(KK)*F2
      MAT2(7) = DISTI*(DL_DM-VEC(KJ)*VEC(KK))*F1
     *  +VEC(KJ)*VEC(KK)*F2
      MAT2(8) =
     * (VEC(KI)*DL_DM + VEC(KJ)*DK_DM +
     *  VEC(KK)*DK_DL)*(-DISTI*F1+F2)*DISTI
     * - THREE*VEC(KI)*VEC(KJ)*VEC(KK)*(-DISTI*F1+F2)*DISTI
     * + VEC(KI)*VEC(KJ)*VEC(KK)*F3
      IF (NDER.EQ.3) RETURN
C
C
      END SUBROUTINE DFTB_AD_PREP
C
C-----------------------------------------------------------------------
C*MODULE DFTBSK    *DECK DFTB_AD_ADD
C>
C>    @brief Calculate 1st and 2nd derivative with AD
C>
C>    @details Calculate 1st and 2nd derivative of sum with automatic
C>             differentiation, note that this does only simple sum.
C>
C>    @author Yoshio Nishimoto
C>    - Sep, 2014- Subroutine written
C>    @date  Feb, 2018 - Yoshio Nishimoto
C>    - Arbitrary-order derivatives
C>
C>           --- INPUT ---
C>    @param NDER  Order of derivative
C>    @param A     Coefficient for the first matrix
C>    @param B     Coefficient for the second matrix
C>    @param MAT2  Second matrix
C>           --- OUTPUT ---
C>    @param MAT1  Sirst matrix
C
      SUBROUTINE DFTB_AD_ADD(NDER,MAT1,A,MAT2,B)
C
      IMPLICIT NONE
C
      INTEGER, INTENT(IN) :: NDER
      DOUBLE PRECISION, INTENT(IN) :: A,B
      DOUBLE PRECISION, INTENT(INOUT) :: MAT1(*),MAT2(*)
C
C     ALL TERMS ARE SIMPLE SUM
C
C     TMP1 :: 0TH ORDER DERIVATIVE
C     TMP2 :: 1ST ORDER DERIVATIVE FOR KI
C     TMP3 :: 1ST ORDER DERIVATIVE FOR KJ
C     TMP4 :: 2ND ORDER DERIVATIVE FOR KI-KJ
C
      
      IF (A.NE.1.0D+00) CALL DSCAL(2**NDER,A,MAT1,1)
      CALL DAXPY(2**NDER,B,MAT2,1,MAT1,1)
C     MAT1(1) = A*MAT1(1) + B*MAT2(1)
C     MAT1(2) = A*MAT1(2) + B*MAT2(2)
C     MAT1(3) = A*MAT1(3) + B*MAT2(3)
C     MAT1(4) = A*MAT1(4) + B*MAT2(4)
C
      END SUBROUTINE DFTB_AD_ADD
C
C-----------------------------------------------------------------------
C*MODULE DFTBSK    *DECK DFTB_AD_MUL
C>
C>    @brief Calculate 1st and 2nd derivative with AD
C>
C>    @details Calculate 1st and 2nd derivative of multiplication with
C>             automatic differentiation.
C>
C>    @author Yoshio Nishimoto
C>    - Sep, 2014- Subroutine written
C>    @date  Feb, 2018 - Yoshio Nishimoto
C>    - Third-order derivatives
C>
C>           --- INPUT ---
C>    @param NDER Order of derivative
C>    @param MAT2 Second matrix
C>           --- OUTPUT ---
C>    @param MAT1 First matrix
C
      SUBROUTINE DFTB_AD_MUL(NDER,MAT1,MAT2)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DOUBLE PRECISION MAT1(*),MAT2(*)
C
C     TMP1 :: 0TH ORDER DERIVATIVE
C       F = F1 * F2
C     TMP2 :: 1ST ORDER DERIVATIVE FOR KI
C       dF/dK = F1*(dF2/dK) + (dF1/dK)*F2
C     TMP3 :: 1ST ORDER DERIVATIVE FOR KJ
C       dF/dL = F1*(dF2/dL) + (dF1/dL)*F2
C     TMP4 :: 2ND ORDER DERIVATIVE FOR KI-KJ
C       d2F/dKdL = F1*(d2F2/dKdL) + (dF1/dK)*(dF2/dL)
C                + (dF1/dL)*(dF2/dK) + (d2F1/dKdL)*F2
C     TMP8 :: 3RD ORDER DERIVATIVE FOR KI-KJ-KK
C       d3F/dKdLdM
C          = F1*F2^kml
C          + F1^k*F2^ml + F1^l*F2^mk + F1^m*F2^kl
C          + F1^kl*F2^m + F1^lm*F2^k + F1^mk*F2^l
C          + F1^klm*F2
C
      TMP1 = MAT1(1)*MAT2(1)
      IF (NDER.EQ.0) GO TO 10
C
      TMP2 = MAT1(1)*MAT2(2) + MAT1(2)*MAT2(1)
      IF (NDER.EQ.1) GO TO 10
C
      TMP3 = MAT1(1)*MAT2(3) + MAT1(3)*MAT2(1)
      TMP4 = MAT1(4)*MAT2(1) + MAT1(3)*MAT2(2)
     *     + MAT1(2)*MAT2(3) + MAT1(1)*MAT2(4)
      IF (NDER.EQ.2) GO TO 10
C
      TMP5 = MAT1(1)*MAT2(5) + MAT1(5)*MAT2(1)
      TMP6 = MAT1(1)*MAT2(6) + MAT1(2)*MAT2(5)
     *     + MAT1(5)*MAT2(2) + MAT1(6)*MAT2(1)
      TMP7 = MAT1(1)*MAT2(7) + MAT1(3)*MAT2(5)
     *     + MAT1(5)*MAT2(3) + MAT1(7)*MAT2(1)
      TMP8 = MAT1(1)*MAT2(8)
     *     + MAT1(2)*MAT2(7) + MAT1(3)*MAT2(6) + MAT1(5)*MAT2(4)
     *     + MAT1(4)*MAT2(5) + MAT1(7)*MAT2(2) + MAT1(6)*MAT2(3)
     *     + MAT1(8)*MAT2(1)
C
   10 CONTINUE
C
      MAT1(1) = TMP1
      IF (NDER.EQ.0) RETURN
C
      MAT1(2) = TMP2
      IF (NDER.EQ.1) RETURN
C
      MAT1(3) = TMP3
      MAT1(4) = TMP4
      IF (NDER.EQ.2) RETURN
C
      MAT1(5) = TMP5
      MAT1(6) = TMP6
      MAT1(7) = TMP7
      MAT1(8) = TMP8
      IF (NDER.EQ.3) RETURN
C
      END SUBROUTINE DFTB_AD_MUL
C
C-----------------------------------------------------------------------
C*MODULE DFTBSK    *DECK DFTB_AD_VAL
C>
C>    @brief Calculate 1st and 2nd derivative with AD
C>
C>    @details Calculate 1st and 2nd derivative of multiplication with
C>             automatic differentiation.
C>
C>    @author Yoshio Nishimoto
C>    - Sep, 2014- Subroutine written
C>    @date  Feb, 2018 - Yoshio Nishimoto
C>    - Arbitrary-order derivatives
C>
C>           --- INPUT ---
C>    @param NDER Order of derivative to be calculated
C>    @param TMP Array of 0th, 1st, and 2nd derivatives
C>    @param ALPHA,BETA Coefficients of first and second matrix
C>           --- OUTPUT ---
C>    @param VAL Returned value
C
      SUBROUTINE DFTB_AD_VAL(NDER,VAL,TMP,ALPHA,BETA)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DOUBLE PRECISION :: TMP(*)
C
      IF (NDER.EQ.0) THEN
        VAL = VAL*ALPHA + TMP(2**NDER)*BETA
      ELSE
        VAL = VAL*ALPHA - TMP(2**NDER)*BETA
      END IF
C
      END SUBROUTINE DFTB_AD_VAL
C
C-----------------------------------------------------------------------
C*MODULE DFTBSK    *DECK DFTB_SKSSA
C>
C>    @brief Slater-Koster transformation of s-s interaction
C>
C>    @details first and second order Slater-Koster transformation of
C>             s-s interaction analytically
C>             Symbolic differentiation for 1st- and 2nd-order, and
C>             automatic differentiation for higher-order derivatives
C>
C>    @author Yoshio Nishimoto
C>    - Sep, 2014- Subroutine written
C>    @date  Feb, 2018 - Yoshio Nishimoto
C>    - Arbitrary-order derivatives
C>
C>           --- INPUT ---
C>    @param NDER  Order of derivative to be calculated
C>    @param KI,KJ,KK,KL Axis to be differentiated
C>    @param DISTI Inverse of distance
C>    @param X     Vector of A-B
C>    @param DD,DD1,DD2,DD3,DD4 0th-, 1st-, 2nd-, 3rd-, and 4th-order
C>           derivative matrices with respect to disance
C>           --- OUTPUT ---
C>    @param MAT   Block Hamiltonian/overlap elements in Cartesian
C
      SUBROUTINE DFTB_SKSSA(NDER,KI,KJ,KK,KL,DISTI,X,DD,DD1,DD2,DD3,DD4,
     *                      MAT)
C
      IMPLICIT NONE
C
      INTEGER         ,INTENT(IN) :: NDER,KI,KJ,KK,KL
      DOUBLE PRECISION,INTENT(IN) :: DISTI,X(6)
      DOUBLE PRECISION,INTENT(IN) :: DD(1),DD1(1),DD2(1),DD3(1),DD4(1)
      DOUBLE PRECISION,INTENT(INOUT) :: MAT(1,1)
      DOUBLE PRECISION :: TMP1(2**NDER,3),TMP2(2**NDER)
C
      IF (NDER.EQ.0) THEN
        MAT(1,1) = DD(1)
      ELSE IF (NDER.EQ.1) THEN
        MAT(1,1) = -DD1(1)*X(KI)
      ELSE IF (NDER.EQ.2) THEN
        IF (KI.EQ.KJ) THEN
          MAT(1,1) = -DISTI*(1.0D+00-X(KI)*X(KJ))*DD1(1)
     *      -X(KI)*X(KJ)*DD2(1)
        ELSE
          MAT(1,1) = +DISTI*X(KI)*X(KJ)*DD1(1)-X(KI)*X(KJ)*DD2(1)
        END IF
      ELSE
        CALL DFTB_AD_PREP(NDER,KI,KJ,KK,KL,DISTI,X,DD,DD1,DD2,DD3,DD4,
     *                    TMP1,TMP2)
        CALL DFTB_AD_VAL(NDER,MAT(1,1),TMP2,0.0D+00,1.0D+00)
      END IF
C
C     CALL DFTB_AD_PREP(KI,KJ,DISTI,X,DD,DD1,DD2,TMP1,TMP2)
C
      RETURN
C
      END SUBROUTINE DFTB_SKSSA
C
C-----------------------------------------------------------------------
C*MODULE DFTBSK    *DECK DFTB_SKSPA
C>
C>    @brief Slater-Koster transformation of s-p interaction
C>
C>    @details first and second order Slater-Koster transformation of
C>             s-p interaction analytically
C>             Symbolic differentiation for 1st- and 2nd-order, and
C>             automatic differentiation for higher-order derivatives
C>
C>    @author Yoshio Nishimoto
C>    - Sep, 2014- Subroutine written
C>    @date  Feb, 2018 - Yoshio Nishimoto
C>    - Arbitrary-order derivatives
C>
C>           --- INPUT ---
C>    @param NDER  Order of derivative to be calculated
C>    @param KI,KJ,KK,KL Axis to be differentiated
C>    @param DISTI Inverse of distance
C>    @param X     Vector of A-B
C>    @param DD,DD1,DD2,DD3,DD4 0th-, 1st-, 2nd-, 3rd-, and 4th-order
C>           derivative matrices with respect to disance
C>           --- OUTPUT ---
C>    @param MAT1,MAT2 Block Hamiltonian/overlap elements in Cartesian
C
      SUBROUTINE DFTB_SKSPA(NDER,KI,KJ,KK,KL,DISTI,X,DD,DD1,DD2,DD3,DD4,
     *                      MAT1,MAT2)
C
      IMPLICIT NONE
C
      INTEGER         ,INTENT(IN) :: NDER,KI,KJ,KK,KL
      DOUBLE PRECISION,INTENT(IN) :: DISTI,X(3)
      DOUBLE PRECISION,INTENT(IN) :: DD(1),DD1(1),DD2(1),DD3(1),DD4(1)
      DOUBLE PRECISION,INTENT(INOUT) :: MAT1(1,3), MAT2(3,1)
C
      INTEGER :: I,NNN
      DOUBLE PRECISION :: TMP
      DOUBLE PRECISION :: TMP1(2**NDER,3),TMP2(2**NDER),TTMP(2**NDER)
C
C     INDEX ::
C       1 -> X
C       2 -> Y
C       3 -> Z
C
C     THIS ORDER IS DIFFERENT FROM THAT OF DFTB+
C     IN THE CASE OF DFTB+, THE ORDER IS 1 -> Y, 2 -> Z, 3 -> X
C
      IF (NDER.EQ.0) THEN
        DO I = 1, 3
          MAT1(1,I) =  X(I) * DD(1)
          MAT2(I,1) = -MAT1(1,I)
        END DO
      ELSE IF (NDER.EQ.1) THEN
        DO I = 1, 3
          IF (I.EQ.KI) THEN
            MAT1(1,I) = DISTI*(1.0D+00-X(I)*X(I))*DD(1)
     *                + X(I)*X(I)*DD1(1)
            MAT1(1,I) = -MAT1(1,I)
            MAT2(I,1) = -MAT1(1,I)
          ELSE
            MAT1(1,I) = X(I)*X(KI)*(DD1(1)-DISTI*DD(1))
            MAT1(1,I) = -MAT1(1,I)
            MAT2(I,1) = -MAT1(1,I)
          END IF
        END DO
      ELSE IF (NDER.EQ.2) THEN
        DO I = 1, 3
          TMP = 0.0D+00
          IF (KI.EQ.KJ) TMP = TMP + X(I)
          IF (I.EQ.KI) TMP = TMP + X(KJ)
          IF (I.EQ.KJ) TMP = TMP + X(KI)
            MAT1(1,I) = +DISTI*(DD1(1)-DISTI*DD(1))*
     *        (TMP-3.0D+00*X(I)*X(KI)*X(KJ)) + X(I)*X(KI)*X(KJ)*DD2(1)
            MAT1(1,I) = -MAT1(1,I)
            MAT2(I,1) = -MAT1(1,I)
        END DO
      ELSE
        CALL DFTB_AD_PREP(NDER,KI,KJ,KK,KL,DISTI,X,DD(1),DD1(1),DD2(1),
     *                    DD3(1),DD4(1),TMP1,TMP2)
        NNN = 2**NDER
        DO I = 1, 3
          CALL DCOPY(NNN,TMP1(1,I),1,TTMP,1)
          CALL DFTB_AD_MUL(NDER,TTMP,TMP2)
          CALL DFTB_AD_VAL(NDER,MAT1(1,I),TTMP,0.0D+00,1.0D+00)
          MAT2(I,1) = -MAT1(1,I)
        END DO
      END IF
C
C     write (*,*) "skspa nder = ", nder
C     CALL DFTB_AD_PREP(KI,KJ,DISTI,X,DD,DD1,DD2,TMP1,TMP2)
C     write (*,*) "result"
C     DO I = 1, 3
C       CALL VCLR(TTMP,1,4)
C       CALL DCOPY(4,TMP1(1,i),1,TTMP,1)
C       CALL DFTB_AD_MUL(TTMP,TMP2)
C       if (nder.eq.0) then
C         write (*,*) mat1(1,I),ttmp(1)
C       else if (nder.eq.1) then
C         write (*,*) mat1(1,I),ttmp(2)
C       else
C         write (*,*) mat1(1,I),ttmp(4)
C       end if
C     END DO
      RETURN
C
      END SUBROUTINE DFTB_SKSPA
C
C-----------------------------------------------------------------------
C*MODULE DFTBSK    *DECK DFTB_SKSDA
C>
C>    @brief Slater-Koster transformation of s-d interaction
C>
C>    @details first and second order Slater-Koster transformation of
C>             s-d interaction analytically
C>             This is automatic differentiation.
C>
C>    @author Yoshio Nishimoto
C>    - Sep, 2014- Subroutine written
C>    @date  Feb, 2018 - Yoshio Nishimoto
C>    - Arbitrary-order derivatives
C>
C>           --- INPUT ---
C>    @param NDER  Order of derivative to be calculated
C>    @param KI,KJ,KK,KL Axis to be differentiated
C>    @param DISTI Inverse of distance
C>    @param X     Vector of A-B
C>    @param DD,DD1,DD2,DD3,DD4 0th-, 1st-, 2nd-, 3rd-, and 4th-order
C>           derivative matrices with respect to disance
C>           --- OUTPUT ---
C>    @param MAT1,MAT2 Block Hamiltonian/overlap elements in Cartesian
C
      SUBROUTINE DFTB_SKSDA(NDER,KI,KJ,KK,KL,DISTI,X,DD,DD1,DD2,DD3,DD4,
     *                      MAT1,MAT2)
C
      IMPLICIT NONE
C
      INTEGER         ,INTENT(IN) :: NDER,KI,KJ,KK,KL
      DOUBLE PRECISION,INTENT(IN) :: DISTI,X(3)
      DOUBLE PRECISION,INTENT(IN) :: DD(1),DD1(1),DD2(1),DD3(1),DD4(1)
      DOUBLE PRECISION,INTENT(INOUT) :: MAT1(1,5), MAT2(5,1)
C
      INTEGER :: I,NNN
      DOUBLE PRECISION, PARAMETER :: ONE=1.0D+00, HALF=0.5D+00
      DOUBLE PRECISION :: TMP1(2**NDER,3),TMP2(2**NDER),TTMP(2**NDER),
     *                    TTMP2(2**NDER)
      DOUBLE PRECISION R3
C
      R3 = SQRT(3.0D+00)
      NNN =2**NDER
C
C     INDEX::
C       1 -> XY
C       2 -> YZ
C       3 -> ZX
C       4 -> X^2 - Y^2
C       5 -> 3*Z^2 - R^2
C
C     write (*,*) "sksda nder = ", nder
      CALL DFTB_AD_PREP(NDER,KI,KJ,KK,KL,DISTI,X,DD(1),DD1(1),DD2(1),
     *                  DD3(1),DD4(1),TMP1,TMP2)
C
      DO I = 1, 3
        CALL DCOPY(NNN,TMP1(1,I),1,TTMP,1)
        IF (I.EQ.1.OR.I.EQ.2) THEN
          CALL DFTB_AD_MUL(NDER,TTMP,TMP1(1,I+1))
        ELSE
          CALL DFTB_AD_MUL(NDER,TTMP,TMP1(1,1))
        END IF
        CALL DFTB_AD_MUL(NDER,TTMP,TMP2)
        CALL DFTB_AD_VAL(NDER,MAT1(1,I),TTMP,0.0D+00,R3)
      END DO
C
      CALL DCOPY(NNN,TMP1(1,1),1,TTMP,1)
      CALL DFTB_AD_MUL(NDER,TTMP,TMP1(1,1))
      CALL DCOPY(NNN,TMP1(1,2),1,TTMP2,1)
      CALL DFTB_AD_MUL(NDER,TTMP2,TMP1(1,2))
      CALL DFTB_AD_ADD(NDER,TTMP,ONE,TTMP2,-ONE) !! TTMP = BETA DERIVATIVE
      CALL DFTB_AD_MUL(NDER,TTMP,TMP2)
      CALL DFTB_AD_VAL(NDER,MAT1(1,4),TTMP,0.0D+00,R3*HALF)
C
      CALL DCOPY(NNN,TMP1(1,1),1,TTMP,1)
      CALL DFTB_AD_MUL(NDER,TTMP,TMP1(1,1))
      CALL DCOPY(NNN,TMP1(1,2),1,TTMP2,1)
      CALL DFTB_AD_MUL(NDER,TTMP2,TMP1(1,2))
      CALL DFTB_AD_ADD(NDER,TTMP,ONE,TTMP2,ONE) !! TTMP = ALPHA DERIVATIVE
      CALL DCOPY(NNN,TMP1(1,3),1,TTMP2,1)
      CALL DFTB_AD_MUL(NDER,TTMP2,TMP1(1,3)) !! TTMP2 = Z^2
      CALL DFTB_AD_ADD(NDER,TTMP,-HALF,TTMP2,ONE)
      CALL DFTB_AD_MUL(NDER,TTMP,TMP2)
      CALL DFTB_AD_VAL(NDER,MAT1(1,5),TTMP,0.0D+00,ONE)
C
      DO I = 1, 5
        MAT2(I,1) = MAT1(1,I)
      END DO
C
      RETURN
C
      END SUBROUTINE DFTB_SKSDA
C
C-----------------------------------------------------------------------
C*MODULE DFTBSK    *DECK DFTB_SKPPA
C>
C>    @brief Slater-Koster transformation of p-p interaction
C>
C>    @details first and second order Slater-Koster transformation of
C>             p-p interaction analytically
C>             Symbolic differentiation for 1st- and 2nd-order, and
C>             automatic differentiation for higher-order derivatives
C>
C>    @author Yoshio Nishimoto
C>    - Sep, 2014- Subroutine written
C>    @date  Feb, 2018 - Yoshio Nishimoto
C>    - Arbitrary-order derivatives
C>
C>           --- INPUT ---
C>    @param NDER  Order of derivative to be calculated
C>    @param KI,KJ,KK,KL Axis to be differentiated
C>    @param DISTI Inverse of distance
C>    @param X     Vector of A-B
C>    @param DD,DD1,DD2,DD3,DD4 0th-, 1st-, 2nd-, 3rd-, and 4th-order
C>           derivative matrices with respect to disance
C>           --- OUTPUT ---
C>    @param MAT1 Block Hamiltonian/overlap elements in Cartesian
C
      SUBROUTINE DFTB_SKPPA(NDER,KI,KJ,KK,KL,DISTI,X,DD,DD1,DD2,DD3,
     *                      DD4,MAT)
C
      IMPLICIT NONE
C
      INTEGER         ,INTENT(IN) :: NDER,KI,KJ,KK,KL
      DOUBLE PRECISION,INTENT(IN) :: DISTI,X(6)
      DOUBLE PRECISION,INTENT(IN) :: DD(2),DD1(2),DD2(2),DD3(2),DD4(2)
      DOUBLE PRECISION,INTENT(INOUT) :: MAT(3,3)
C
C TM : TRANSFORMATION ...
C PROD : PRODUCT
C        PROD(1) = X * X
C        PROD(2) = Y * Y
C        PROD(3) = Z * Z
C        PROD(4) = X * Y
C        PROD(5) = Y * Z
C        PROD(6) = Z * X
C
      DOUBLE PRECISION :: TM(6), PROD(6),DX_DL,DX_DK,DK_DL,
     * DY_DK,DY_DL,TMP
      DOUBLE PRECISION TMP1(2**NDER,3),TMP2(2**NDER),TTMP(2**NDER),
     *  TTMP2(2**NDER)
      DOUBLE PRECISION, PARAMETER :: ZERO=0.0D+00,ONE=1.0D+00,
     * TWO=2.0D+00
      INTEGER :: I, J, K, L, NNN
C
      IF (NDER.EQ.0) THEN
        DO I = 1, 3
          PROD(I) = X(I)*X(I)
          PROD(I+3) = X(I) * X(I+1)
        END DO
C
C DD(1) : SIGMA INTERACTION
C DD(2) : PI    INTERACTION
C
        DO I = 1, 3
          TM(I)   = PROD(I) * DD(1) + (1.0D+00-PROD(I)) * DD(2)
          TM(I+3) = PROD(I+3) * (DD(1)-DD(2))
        END DO
C
C PX-PX = X^2*sigma + (1-X^2)*pi -> TM(1)
C PY-PY = Y^2*sigma + (1-Y^2)*pi -> TM(2)
C PZ-PZ = Z^2*sigma + (1-Y^2)*pi -> TM(3)
C PX-PY = X*Y*(sigma-pi) -> TM(4) ... and so on
C
        DO I = 1, 3
          DO J = 1, I
            K = I - J
            L = 3*K - K*(K-1)/2 + J
            MAT(I,J) = TM(L)
            MAT(J,I) = TM(L)
          END DO
        END DO
      ELSE IF (NDER.EQ.1) THEN
C       write (*,*) "k = ", k
        DO I = 1, 3
          DO J = 1, 3
            IF (I.EQ.J) THEN
              IF (I.EQ.KI) THEN
                MAT(I,I) =
     *            2.0D+00*X(I)*DISTI*(1.0D+00-X(I)*X(I))*(DD(1)-DD(2))
     *              + X(I)*X(I)*X(I)*(DD1(1)-DD1(2)) + X(I)*DD1(2)
C               write (*,*) i,i,mat(i,i)
              ELSE
                MAT(I,I) =
     *            X(I)*X(I)*X(KI)
     *            *((DD1(1)-DD1(2))-2.0D+00*DISTI*(DD(1)-DD(2)))
     *            +X(KI)*DD1(2)
C               write (*,*) i,i,mat(i,i)
C               write (*,*) x(i),x(k),disti
              END IF
            ELSE
              IF (I.NE.KI.AND.J.NE.KI) THEN
C               MAT(I,J) = X(I)*X(J)*X(K)*(DD1(1)-DD1(2))
                MAT(I,J) = X(I)*X(J)*X(KI)*
     *            ((DD1(1)-DD1(2)) - 2.0D+00*DISTI*(DD(1)-DD(2)))
              ELSE
                IF (I.EQ.KI) THEN
                  MAT(I,J) =
     *              DISTI*X(J)*(1.0D+00-2.0D+00*X(I)*X(I))*(DD(1)-DD(2))
     *              + X(I)*X(J)*X(KI)*(DD1(1)-DD1(2))
                ELSE IF (J.EQ.KI) THEN
                  MAT(I,J) =
     *              DISTI*X(I)*(1.0D+00-2.0D+00*X(J)*X(J))*(DD(1)-DD(2))
     *              + X(I)*X(J)*X(KI)*(DD1(1)-DD1(2))
                END IF
              END IF
            END IF
          END DO
        END DO
        CALL DSCAL(9,-1.0D+00,MAT,1)
      ELSE IF (NDER.EQ.2) THEN
        DO I = 1, 3 !! x,y,z for p-orbital
          DO J = 1, 3 !! x,y,z for p-orbital
            IF (I.EQ.J) THEN
              DX_DK = ZERO
              DX_DL = ZERO
              DK_DL = ZERO
              IF (I.EQ.KI) DX_DK = ONE
              IF (I.EQ.KJ) DX_DL = ONE
              IF (KI.EQ.KJ) DK_DL = ONE
              MAT(I,J) = TWO*DISTI*DISTI*(DX_DL*DX_DK
     *          -TWO*X(I)*(X(KI)*DX_DL+X(KJ)*DX_DK)
     *          -X(I)*X(I)*DK_DL
     *          +4.0D+00*X(I)*X(I)*X(KI)*X(KJ))*(DD(1)-DD(2))
     *          +DISTI*X(I)*(TWO*(X(KI)*DX_DL+X(KJ)*DX_DK)
     *          +X(I)*DK_DL-5.0D+00*X(I)*X(KI)*X(KJ))
     *          *(DD1(1)-DD1(2))
     *          +X(I)*X(I)*X(KI)*X(KJ)*(DD2(1)-DD2(2))
     *          +DISTI*(DK_DL-X(KI)*X(KJ))*DD1(2)
     *          +X(KI)*X(KJ)*DD2(2)
            ELSE
              DX_DK = ZERO
              DX_DL = ZERO
              DY_DK = ZERO
              DY_DL = ZERO
              DK_DL = ZERO
              IF (I.EQ.KI)  DX_DK = ONE
              IF (I.EQ.KJ)  DX_DL = ONE
              IF (J.EQ.KI)  DY_DK = ONE
              IF (J.EQ.KJ)  DY_DL = ONE
              IF (KI.EQ.KJ) DK_DL = ONE
              TMP = X(J)*X(KJ)*DX_DK+X(I)*X(KJ)*DY_DK+X(J)*X(KI)*DX_DL
     *            + X(I)*X(KI)*DY_DL+X(I)*X(J)*DK_DL
              MAT(I,J) = DISTI*DISTI*(DX_DK*DY_DL+DX_DL*DY_DK
     *          -TWO*TMP+8.0D+00*X(I)*X(J)*X(KI)*X(KJ))*(DD(1)-DD(2))
     *          +DISTI*(TMP-5.0D+00*X(I)*X(J)*X(KI)*X(KJ))
     *            *(DD1(1)-DD1(2))
     *          +X(I)*X(J)*X(KI)*X(KJ)*(DD2(1)-DD2(2))
            END IF
          END DO
        END DO
        CALL DSCAL(9,-1.0D+00,MAT,1)
      ELSE
        NNN = 2**NDER
        DO I = 1, 3
          DO J = 1, 3
            MAT(I,J) = 0.0D+00
            IF (I.EQ.J) THEN
              CALL DFTB_AD_PREP(NDER,KI,KJ,KK,KL,DISTI,X,DD(1),DD1(1),
     *             DD2(1),DD3(1),DD4(1),TMP1,TMP2)
              CALL DCOPY(NNN,TMP1(1,I),1,TTMP,1)
              CALL DFTB_AD_MUL(NDER,TTMP,TMP1(1,J))
              CALL DFTB_AD_MUL(NDER,TTMP,TMP2)
C
              CALL DFTB_AD_PREP(NDER,KI,KJ,KK,KL,DISTI,X,DD(2),DD1(2),
     *             DD2(2),DD3(2),DD4(2),TMP1,TMP2)
              CALL DCOPY(NNN,TMP1(1,I),1,TTMP2,1)
              CALL DFTB_AD_MUL(NDER,TTMP2,TMP1(1,I))
              CALL DSCAL(NNN,-ONE,TTMP2,1)
              TTMP2(1) = TTMP2(1) + ONE
              CALL DFTB_AD_MUL(NDER,TTMP2,TMP2)
              CALL DFTB_AD_ADD(NDER,TTMP,ONE,TTMP2,ONE)
            ELSE
              CALL DFTB_AD_PREP(NDER,KI,KJ,KK,KL,DISTI,X,DD(1)-DD(2),
     *          DD1(1)-DD1(2),DD2(1)-DD2(2),DD3(1)-DD3(2),DD4(1)-DD4(2),
     *          TMP1,TMP2)
              CALL DCOPY(NNN,TMP1(1,I),1,TTMP,1)
              CALL DFTB_AD_MUL(NDER,TTMP,TMP1(1,J))
              CALL DFTB_AD_MUL(NDER,TTMP,TMP2)
            END IF
            CALL DFTB_AD_VAL(NDER,MAT(I,J),TTMP,ZERO,ONE)
          END DO
        END DO
      END IF
C
      return
C
      END SUBROUTINE DFTB_SKPPA
C
C-----------------------------------------------------------------------
C*MODULE DFTBSK    *DECK DFTB_SKPDA
C>
C>    @brief Slater-Koster transformation of p-d interaction
C>
C>    @details first and second order Slater-Koster transformation of
C>             p-d interaction analytically
C>             This is automatic differentiation.
C>
C>    @author Yoshio Nishimoto
C>    - Sep, 2014- Subroutine written
C>    @date  Feb, 2018 - Yoshio Nishimoto
C>    - Arbitrary-order derivatives
C>
C>           --- INPUT ---
C>    @param NDER  Order of derivative to be calculated
C>    @param KI,KJ,KK,KL Axis to be differentiated
C>    @param DISTI Inverse of distance
C>    @param X     Vector of A-B
C>    @param DD,DD1,DD2,DD3,DD4 0th-, 1st-, 2nd-, 3rd-, and 4th-order
C>           derivative matrices with respect to disance
C>           --- OUTPUT ---
C>    @param MAT1,MAT2 Block Hamiltonian/overlap elements in Cartesian
C
      SUBROUTINE DFTB_SKPDA(NDER,KI,KJ,KK,KL,DISTI,X,DD,DD1,DD2,DD3,DD4,
     *                      MAT1,MAT2)
C
      IMPLICIT NONE
C
      INTEGER         ,INTENT(IN) :: NDER,KI,KJ,KK,KL
      DOUBLE PRECISION,INTENT(IN) :: DISTI,X(3)
      DOUBLE PRECISION,INTENT(IN) :: DD(2),DD1(2),DD2(2),DD3(2),DD4(2)
      DOUBLE PRECISION,INTENT(INOUT) :: MAT1(3,5), MAT2(5,3)
C
      INTEGER :: I,J,K,NNN
      DOUBLE PRECISION, PARAMETER :: ZERO=0.0D+00,ONE=1.0D+00,
     *  HALF=0.5D+00,TWO=2.0D+00
      DOUBLE PRECISION :: TMP1(2**NDER,3),TMPS(2**NDER),TMPP(2**NDER),
     * TTMP(2**NDER),TTMP2(2**NDER),TMPA(2**NDER),TMPB(2**NDER),
     * TMPXYZ(2**NDER,4),TMPSQ(2**NDER,3),TTMP3(2**NDER)
      DOUBLE PRECISION R3
C
      R3 = SQRT(3.0D+00)
      NNN = 2**NDER
C
C     MAT1(1,1) : P(X)-D(XY)
C     MAT1(1,2) : P(X)-D(YZ)
C     MAT1(1,3) : P(X)-D(ZX)
C     MAT1(1,4) : P(X)-D(X^2-Y^2)
C     MAT1(1,5) : P(X)-D(3Z^2-R^2)
C
      CALL DFTB_AD_PREP(NDER,KI,KJ,KK,KL,DISTI,X,DD(1),DD1(1),DD2(1),
     *                  DD3(1),DD4(1),TMP1,TMPS)
      CALL DFTB_AD_PREP(NDER,KI,KJ,KK,KL,DISTI,X,DD(2),DD1(2),DD2(2),
     *                  DD3(1),DD4(1),TMP1,TMPP)
C
      CALL DCOPY(3*NNN,TMP1(1,1),1,TMPXYZ(1,1),1)
      CALL DFTB_AD_MUL(NDER,TMPXYZ(1,1),TMP1(1,2)) !! TMPXYZ(1,1) = XY
      CALL DCOPY(NNN,TMPXYZ(1,1),1,TMPXYZ(1,4),1)
      CALL DFTB_AD_MUL(NDER,TMPXYZ(1,4),TMP1(1,3)) !! TMPXYZ(1,4) = XYZ
      CALL DFTB_AD_MUL(NDER,TMPXYZ(1,2),TMP1(1,3)) !! TMPXYZ(1,2) = YZ
      CALL DFTB_AD_MUL(NDER,TMPXYZ(1,3),TMP1(1,1)) !! TMPXYZ(1,3) = ZX
C
      CALL DCOPY(3*NNN,TMP1,1,TMPSQ,1)
      DO I = 1, 3
        CALL DFTB_AD_MUL(NDER,TMPSQ(1,I),TMP1(1,I)) !! X^2,Y^2,Z^2
      END DO
C
      CALL DCOPY(NNN,TMPSQ(1,1),1,TMPA,1)
      CALL DCOPY(NNN,TMPSQ(1,1),1,TMPB,1)
      CALL DFTB_AD_ADD(NDER,TMPA,ONE,TMPSQ(1,2),+ONE) !! ALPHA
      CALL DFTB_AD_ADD(NDER,TMPB,ONE,TMPSQ(1,2),-ONE) !! BETA
C
      CALL DCOPY(NNN,TMPS,1,TTMP3,1)
      CALL DFTB_AD_ADD(NDER,TTMP3,R3,TMPP,-TWO) !! SQRT(3)C(s)-2C(p)
      DO I = 1, 3 !! x,y,z
        DO J = 1, 3 !! xy,yz,zx
          CALL DCOPY(NNN,TTMP3,1,TTMP,1)
          IF ((I.EQ.1.AND.J.EQ.2).OR.(I.EQ.2.AND.J.EQ.3)
     *        .OR.(I.EQ.3.AND.J.EQ.1)) THEN
            CALL DFTB_AD_MUL(NDER,TTMP,TMPXYZ(1,4))
            CALL DFTB_AD_VAL(NDER,MAT1(I,J),TTMP,ZERO,ONE)
          ELSE
            IF (I.EQ.1) THEN
              IF (J.EQ.1) K = 2
              IF (J.EQ.3) K = 3
            ELSE IF (I.EQ.2) THEN
              IF (J.EQ.1) K = 1
              IF (J.EQ.2) K = 3
            ELSE IF (I.EQ.3) THEN
              IF (J.EQ.2) K = 2
              IF (J.EQ.3) K = 1
            END IF
            CALL DCOPY(NNN,TMP1(1,K),1,TTMP2,1)
            CALL DFTB_AD_MUL(NDER,TTMP2,TMPP)
            CALL DFTB_AD_VAL(NDER,MAT1(I,J),TTMP2,ZERO,ONE)
            CALL DCOPY(NNN,TTMP3,1,TTMP,1)
            CALL DFTB_AD_MUL(NDER,TTMP,TMPSQ(1,I))
            CALL DFTB_AD_MUL(NDER,TTMP,TMP1(1,K))
            CALL DFTB_AD_VAL(NDER,MAT1(I,J),TTMP,ONE,ONE)
          END IF
        END DO
      END DO
C
      DO I = 1, 3
        CALL DCOPY(NNN,TMPB,1,TTMP,1)
        CALL DFTB_AD_MUL(NDER,TTMP,TMP1(1,I))
        CALL DFTB_AD_MUL(NDER,TTMP,TMPS)
        CALL DFTB_AD_VAL(NDER,MAT1(I,4),TTMP,ZERO,HALF*R3)
      END DO
C     MAT1(1,4)
      CALL DCOPY(NNN,TMPB,1,TTMP,1)
      CALL DSCAL(NNN,-ONE,TTMP,1)
      TTMP(1) = TTMP(1) + ONE
      CALL DFTB_AD_MUL(NDER,TTMP,TMP1(1,1))
      CALL DFTB_AD_MUL(NDER,TTMP,TMPP)
      CALL DFTB_AD_VAL(NDER,MAT1(1,4),TTMP,ONE,ONE)
C     MAT1(2,4)
      CALL DCOPY(NNN,TMPB,1,TTMP,1)
      TTMP(1) = TTMP(1) + ONE
      CALL DFTB_AD_MUL(NDER,TTMP,TMP1(1,2))
      CALL DFTB_AD_MUL(NDER,TTMP,TMPP)
      CALL DFTB_AD_VAL(NDER,MAT1(2,4),TTMP,ONE,-ONE)
C     MAT1(3,4)
      CALL DCOPY(NNN,TMPB,1,TTMP,1)
      CALL DFTB_AD_MUL(NDER,TTMP,TMP1(1,3))
      CALL DFTB_AD_MUL(NDER,TTMP,TMPP)
      CALL DFTB_AD_VAL(NDER,MAT1(3,4),TTMP,ONE,-ONE)
C
      CALL DCOPY(NNN,TMPSQ(1,3),1,TTMP,1)
      CALL DFTB_AD_ADD(NDER,TTMP,ONE,TMPA,-HALF)
      DO I = 1, 3
        CALL DCOPY(NNN,TTMP,1,TTMP2,1)
        CALL DFTB_AD_MUL(NDER,TTMP2,TMP1(1,I))
        CALL DFTB_AD_MUL(NDER,TTMP2,TMPS)
        CALL DFTB_AD_VAL(NDER,MAT1(I,5),TTMP2,ZERO,ONE)
      END DO
C     MAT(1,5)
      CALL DCOPY(NNN,TMPSQ(1,3),1,TTMP,1)
      CALL DFTB_AD_MUL(NDER,TTMP,TMP1(1,1))
      CALL DFTB_AD_MUL(NDER,TTMP,TMPP)
      CALL DFTB_AD_VAL(NDER,MAT1(1,5),TTMP,ONE,-R3)
C     MAT(2,5)
      CALL DCOPY(NNN,TMPSQ(1,3),1,TTMP,1)
      CALL DFTB_AD_MUL(NDER,TTMP,TMP1(1,2))
      CALL DFTB_AD_MUL(NDER,TTMP,TMPP)
      CALL DFTB_AD_VAL(NDER,MAT1(2,5),TTMP,ONE,-R3)
C     MAT(3,5)
      CALL DCOPY(NNN,TMPA,1,TTMP,1)
      CALL DFTB_AD_MUL(NDER,TTMP,TMP1(1,3))
      CALL DFTB_AD_MUL(NDER,TTMP,TMPP)
      CALL DFTB_AD_VAL(NDER,MAT1(3,5),TTMP,ONE,R3)
C
      DO I = 1, 5
        DO J = 1, 3
          MAT2(I,J) = -MAT1(J,I)
        END DO
      END DO
C
      RETURN
C
      END SUBROUTINE DFTB_SKPDA
C
C-----------------------------------------------------------------------
C*MODULE DFTBSK    *DECK DFTB_SKDDA
C>
C>    @brief Slater-Koster transformation of d-d interaction
C>
C>    @details first and second order Slater-Koster transformation of
C>             d-d interaction analytically
C>             This is automatic differentiation.
C>
C>    @author Yoshio Nishimoto
C>    - Sep, 2014- Subroutine written
C>    @data Nov, 2014- Yoshio Nishimoto
C>    - Bug fix of transformation
C>    @date  Feb, 2018 - Yoshio Nishimoto
C>    - Arbitrary-order derivatives
C>
C>           --- INPUT ---
C>    @param NDER  Order of derivative to be calculated
C>    @param KI,KJ,KK,KL Axis to be differentiated
C>    @param DISTI Inverse of distance
C>    @param X     Vector of A-B
C>    @param DD,DD1,DD2,DD3,DD4 0th-, 1st-, 2nd-, 3rd-, and 4th-order
C>           derivative matrices with respect to disance
C>           --- OUTPUT ---
C>    @param MAT Block Hamiltonian/overlap elements in Cartesian
C
      SUBROUTINE DFTB_SKDDA(NDER,KI,KJ,KK,KL,DISTI,X,DD,DD1,DD2,DD3,DD4,
     *                      MAT)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      INTEGER         ,INTENT(IN) :: NDER,KI,KJ
      DOUBLE PRECISION,INTENT(IN) :: DISTI,X(3)
      DOUBLE PRECISION,INTENT(IN) :: DD(3),DD1(3),DD2(3),DD3(3),DD4(3)
      DOUBLE PRECISION,INTENT(INOUT) :: MAT(5,5)
C
      INTEGER :: I,J,K,NNN
      DOUBLE PRECISION, PARAMETER :: ZERO=0.0D+00,ONE=1.0D+00,
     *  HALF=0.5D+00,TWO=2.0D+00,THREE=3.0D+00,FOUR=4.0D+00
      DOUBLE PRECISION :: TMP1(2**NDER,3),TMPS(2**NDER),TMPP(2**NDER),
     * TTMP(2**NDER),TMPD(2**NDER),TTMP2(2**NDER),TMPA(2**NDER),
     * TMPB(2**NDER),TMPXYZ(2**NDER,4),TMPSQ(2**NDER,3),TTMP3(2**NDER)
C
      R3 = SQRT(3.0D+00)
      NNN = 2**NDER
C
      CALL DFTB_AD_PREP(NDER,KI,KJ,KK,KL,DISTI,X,DD(1),DD1(1),DD2(1),
     *                  DD3(1),DD4(1),TMP1,TMPS)
      CALL DFTB_AD_PREP(NDER,KI,KJ,KK,KL,DISTI,X,DD(2),DD1(2),DD2(2),
     *                  DD3(2),DD4(2),TMP1,TMPP)
      CALL DFTB_AD_PREP(NDER,KI,KJ,KK,KL,DISTI,X,DD(3),DD1(3),DD2(3),
     *                  DD3(3),DD4(3),TMP1,TMPD)
C
      CALL DCOPY(3*NNN,TMP1(1,1),1,TMPXYZ(1,1),1)
      CALL DFTB_AD_MUL(NDER,TMPXYZ(1,1),TMP1(1,2)) !! TMPXYZ(1,1) = XY
      CALL DCOPY(NNN,TMPXYZ(1,1),1,TMPXYZ(1,4),1)
      CALL DFTB_AD_MUL(NDER,TMPXYZ(1,4),TMP1(1,3)) !! TMPXYZ(1,4) = XYZ
      CALL DFTB_AD_MUL(NDER,TMPXYZ(1,2),TMP1(1,3)) !! TMPXYZ(1,2) = YZ
      CALL DFTB_AD_MUL(NDER,TMPXYZ(1,3),TMP1(1,1)) !! TMPXYZ(1,3) = ZX
C
      CALL DCOPY(3*NNN,TMP1,1,TMPSQ,1)
      DO I = 1, 3
        CALL DFTB_AD_MUL(NDER,TMPSQ(1,I),TMP1(1,I)) !! X^2,Y^2,Z^2
      END DO
C
C     MAT(1,1) ~ MAT(3,3)
      CALL DCOPY(NNN,TMPS,1,TTMP2,1)
      CALL DFTB_AD_ADD(NDER,TTMP2,THREE,TMPP,-FOUR)
      CALL DFTB_AD_ADD(NDER,TTMP2,ONE,TMPD,ONE) !! TMP2=3C(s)-4C(p)+C(d)
      CALL DCOPY(NNN,TMPP,1,TTMP3,1)
      CALL DFTB_AD_ADD(NDER,TTMP3,ONE,TMPD,-ONE) !! TMP3=C(p)-C(d)
      ! MAT(1,1) is wrong?
      DO I = 1, 3 !! x,y,z
        DO J = 1, 3 !! xy,yz,zx
          CALL DCOPY(NNN,TMPXYZ(1,I),1,TTMP,1)
          CALL DFTB_AD_MUL(NDER,TTMP,TMPXYZ(1,J))
          CALL DFTB_AD_MUL(NDER,TTMP,TTMP2)
          CALL DFTB_AD_VAL(NDER,MAT(I,J),TTMP,ZERO,ONE)
          IF (I.EQ.J) THEN
            CALL DCOPY(NNN,TMPSQ(1,I),1,TTMP,1)
            IF (I.EQ.1) K = 2
            IF (I.EQ.2) K = 3
            IF (I.EQ.3) K = 1
            CALL DFTB_AD_ADD(NDER,TTMP,ONE,TMPSQ(1,K),ONE)
            CALL DFTB_AD_MUL(NDER,TTMP,TMPP)
            CALL DFTB_AD_VAL(NDER,MAT(I,J),TTMP,ONE,ONE)
            IF (I.EQ.1) K = 3
            IF (I.EQ.2) K = 1
            IF (I.EQ.3) K = 2
            CALL DCOPY(NNN,TMPSQ(1,K),1,TTMP,1)
            CALL DFTB_AD_MUL(NDER,TTMP,TMPD)
            CALL DFTB_AD_VAL(NDER,MAT(I,J),TTMP,ONE,ONE)
          ELSE
C     (I,J) = (1,2) -> 3
C           = (1,3) -> 2
C           = (2,1) -> 3
C           = (2,3) -> 1
C           = (3,1) -> 2
C           = (3,2) -> 1
            IF (I.EQ.1.AND.J.EQ.2) K = 3
            IF (I.EQ.1.AND.J.EQ.3) K = 2
            IF (I.EQ.2.AND.J.EQ.1) K = 3
            IF (I.EQ.2.AND.J.EQ.3) K = 1
            IF (I.EQ.3.AND.J.EQ.1) K = 2
            IF (I.EQ.3.AND.J.EQ.2) K = 1
            CALL DCOPY(NNN,TMPXYZ(1,K),1,TTMP,1)
            CALL DFTB_AD_MUL(NDER,TTMP,TTMP3)
            CALL DFTB_AD_VAL(NDER,MAT(I,J),TTMP,ONE,ONE)
          END IF
        END DO
      END DO
C
      CALL DCOPY(NNN,TMPSQ(1,1),1,TMPA,1)
      CALL DCOPY(NNN,TMPSQ(1,1),1,TMPB,1)
      CALL DFTB_AD_ADD(NDER,TMPA,ONE,TMPSQ(1,2),+ONE) !! ALPHA
      CALL DFTB_AD_ADD(NDER,TMPB,ONE,TMPSQ(1,2),-ONE) !! BETA
C
      CALL DCOPY(NNN,TMPS,1,TTMP2,1)
      CALL DFTB_AD_ADD(NDER,TTMP2,1.5D+00,TMPP,-TWO)
      CALL DFTB_AD_ADD(NDER,TTMP2,ONE,TMPD,HALF)
      CALL DFTB_AD_MUL(NDER,TTMP2,TMPB) !! TMP2=B(1.5C(s)-2C(p)+0.5C(d))
      CALL DCOPY(NNN,TMPP,1,TTMP3,1)
      CALL DFTB_AD_ADD(NDER,TTMP3,ONE,TMPD,-ONE) !! TMP3=C(p)-C(d)
C     MAT(4,1) -- d(xy)-d(x^2-y^2)
      CALL DCOPY(NNN,TTMP2,1,TTMP,1)
      CALL DFTB_AD_MUL(NDER,TTMP,TMPXYZ(1,1))
      CALL DFTB_AD_VAL(NDER,MAT(4,1),TTMP,ZERO,ONE)
C     MAT(4,2) -- d(yz)-d(x^2-y^2)
      CALL DCOPY(NNN,TTMP2,1,TTMP,1)
      CALL DFTB_AD_MUL(NDER,TTMP,TMPXYZ(1,2))
      CALL DFTB_AD_VAL(NDER,MAT(4,2),TTMP,ZERO,ONE)
      CALL DCOPY(NNN,TTMP3,1,TTMP,1)
      CALL DFTB_AD_MUL(NDER,TTMP,TMPXYZ(1,2))
      CALL DFTB_AD_VAL(NDER,MAT(4,2),TTMP,ONE,-ONE)
C     MAT(4,3) -- d(zx)-d(x^2-y^2)
      CALL DCOPY(NNN,TTMP2,1,TTMP,1)
      CALL DFTB_AD_MUL(NDER,TTMP,TMPXYZ(1,3))
      CALL DFTB_AD_VAL(NDER,MAT(4,3),TTMP,ZERO,ONE)
      CALL DCOPY(NNN,TTMP3,1,TTMP,1)
      CALL DFTB_AD_MUL(NDER,TTMP,TMPXYZ(1,3))
      CALL DFTB_AD_VAL(NDER,MAT(4,3),TTMP,ONE,ONE)
C
      CALL DCOPY(NNN,TMPSQ(1,3),1,TTMP2,1)
      CALL DFTB_AD_ADD(NDER,TTMP2,ONE,TMPA,-HALF) !! TMP2=z^2-0.5*alpha
      CALL DFTB_AD_MUL(NDER,TTMP2,TMPS)
      DO I = 1, 3
        CALL DCOPY(NNN,TTMP2,1,TTMP,1)
        CALL DFTB_AD_MUL(NDER,TTMP,TMPXYZ(1,I))
        CALL DFTB_AD_VAL(NDER,MAT(5,I),TTMP,ZERO,R3)
      END DO
C     MAT(5,1) -- d(xy)-d(3z^2-r^2)
      CALL DCOPY(NNN,TMPP,1,TTMP,1)
      CALL DFTB_AD_ADD(NDER,TTMP,-TWO*R3,TMPD,HALF*R3)
      CALL DFTB_AD_MUL(NDER,TTMP,TMPXYZ(1,4))
      CALL DFTB_AD_MUL(NDER,TTMP,TMP1(1,3))
      CALL DFTB_AD_VAL(NDER,MAT(5,1),TTMP,ONE,ONE)
      CALL DCOPY(NNN,TMPXYZ(1,1),1,TTMP,1)
      CALL DFTB_AD_MUL(NDER,TTMP,TMPD)
      CALL DFTB_AD_VAL(NDER,MAT(5,1),TTMP,ONE,HALF*R3)
C     MAT(5,2) -- d(yz)-d(3z^2-r^2)
C     MAT(5,3) -- d(zx)-d(3z^2-r^2)
      CALL DCOPY(NNN,TMPA,1,TTMP2,1)
      CALL DFTB_AD_ADD(NDER,TTMP2,ONE,TMPSQ(1,3),-ONE) !! TMP2=alpha-z^2
      CALL DFTB_AD_MUL(NDER,TTMP2,TMPP)
      CALL DCOPY(NNN,TTMP2,1,TTMP,1)
      CALL DFTB_AD_MUL(NDER,TTMP,TMPXYZ(1,2))
      CALL DFTB_AD_VAL(NDER,MAT(5,2),TTMP,ONE,R3)
      CALL DCOPY(NNN,TTMP2,1,TTMP,1)
      CALL DFTB_AD_MUL(NDER,TTMP,TMPXYZ(1,3))
      CALL DFTB_AD_VAL(NDER,MAT(5,3),TTMP,ONE,R3)
C     CONTINUE
      CALL DCOPY(NNN,TMPA,1,TTMP,1)
      CALL DFTB_AD_MUL(NDER,TTMP,TMPXYZ(1,2))
      CALL DFTB_AD_MUL(NDER,TTMP,TMPD)
      CALL DFTB_AD_VAL(NDER,MAT(5,2),TTMP,ONE,-HALF*R3)
      CALL DCOPY(NNN,TMPA,1,TTMP,1)
      CALL DFTB_AD_MUL(NDER,TTMP,TMPXYZ(1,3))
      CALL DFTB_AD_MUL(NDER,TTMP,TMPD)
      CALL DFTB_AD_VAL(NDER,MAT(5,3),TTMP,ONE,-HALF*R3)
C
C     MAT(4,4) -- d(x^2-y^2)-d(x^2-y^2)
      CALL DCOPY(NNN,TMPS,1,TTMP,1)
      CALL DFTB_AD_ADD(NDER,TTMP,7.5D-01,TMPP,-ONE)
      CALL DFTB_AD_ADD(NDER,TTMP,ONE,TMPD,2.5D-01)
      CALL DFTB_AD_MUL(NDER,TTMP,TMPB)
      CALL DFTB_AD_MUL(NDER,TTMP,TMPB)
      CALL DFTB_AD_VAL(NDER,MAT(4,4),TTMP,ZERO,ONE)
      CALL DCOPY(NNN,TMPA,1,TTMP,1)
      CALL DFTB_AD_MUL(NDER,TTMP,TMPP)
      CALL DFTB_AD_VAL(NDER,MAT(4,4),TTMP,ONE,ONE)
      CALL DCOPY(NNN,TMPSQ(1,3),1,TTMP,1)
      CALL DFTB_AD_MUL(NDER,TTMP,TMPD)
      CALL DFTB_AD_VAL(NDER,MAT(4,4),TTMP,ONE,ONE)
C
C     MAT(5,4) -- d(x^2-y^2)-d(3z^2-r^2)
      CALL DCOPY(NNN,TMPS,1,TTMP,1)
      CALL DFTB_AD_ADD(NDER,TTMP,HALF,TMPP,-ONE)
      CALL DFTB_AD_ADD(NDER,TTMP,ONE,TMPD,2.5D-01)
      CALL DFTB_AD_MUL(NDER,TTMP,TMPB)
      CALL DFTB_AD_MUL(NDER,TTMP,TMPSQ(1,3))
      CALL DFTB_AD_VAL(NDER,MAT(5,4),TTMP,ZERO,R3)
      CALL DCOPY(NNN,TMPA,1,TTMP,1)
      CALL DFTB_AD_MUL(NDER,TTMP,TMPB)
      CALL DFTB_AD_MUL(NDER,TTMP,TMPS)
      CALL DFTB_AD_VAL(NDER,MAT(5,4),TTMP,ONE,-R3*2.5D-01)
      CALL DCOPY(NNN,TMPB,1,TTMP,1)
      CALL DFTB_AD_MUL(NDER,TTMP,TMPD)
      CALL DFTB_AD_VAL(NDER,MAT(5,4),TTMP,ONE,R3*2.5D-01)
C
C     MAT(5,5) -- d(3z^2-r^2)-d(3z^2-r^2)
      CALL DCOPY(NNN,TMPA,1,TTMP,1)
      CALL DFTB_AD_ADD(NDER,TTMP,-HALF,TMPSQ(1,3),ONE)
      CALL DFTB_AD_MUL(NDER,TTMP,TTMP) !! TMP=(z^2-0.5*alpha^2)^2
      CALL DFTB_AD_MUL(NDER,TTMP,TMPS)
      CALL DFTB_AD_VAL(NDER,MAT(5,5),TTMP,ZERO,ONE)
      CALL DCOPY(NNN,TMPA,1,TTMP,1)
      CALL DFTB_AD_MUL(NDER,TTMP,TMPSQ(1,3))
      CALL DFTB_AD_MUL(NDER,TTMP,TMPP)
      CALL DFTB_AD_VAL(NDER,MAT(5,5),TTMP,ONE,THREE)
      CALL DCOPY(NNN,TMPA,1,TTMP,1)
      CALL DFTB_AD_MUL(NDER,TTMP,TTMP)
      CALL DFTB_AD_MUL(NDER,TTMP,TMPD)
      CALL DFTB_AD_VAL(NDER,MAT(5,5),TTMP,ONE,7.5D-01)
C
      DO I = 1, 5
        DO J = 1, I-1
          MAT(J,I) = MAT(I,J)
        END DO
      END DO
C
      RETURN
C
      END SUBROUTINE DFTB_SKDDA
C
C-----------------------------------------------------------------------
C*MODULE DFTBSK    *DECK DFTB_SKMEM
C>
C>    @brief   Manage memory to store Slater-Koster parameters
C>
C>    @details Allocate and deallocate the memory above
C>             The first call with MODE=0 should be done from
C>             DFTB_INPUT in the group scope (somehow). The second call
C>             with MODE=1 from the main program is done in the world
C>             scope.
C>
C>    @author Yoshio Nishimoto
C>    - Feb, 2016- Subroutine written
C>    @date   Yoshio Nishimoto
C>    - Oct, 2019- Added DFTB/PBC
C>
C>           --- INPUT ---
C>    @param MODE
C>           0: Allocate memory
C>           1: 0+Read parameters from the unformatted file
C>              It should have been generated in DFTB_INPUT subroutine.
C>           2: Deallocate memory
C>
C
      SUBROUTINE DFTB_SKMEM(MODE)
      use mx_limits, only: mxatm
      USE DFTBPB_MOD,ONLY: PERIOD,MAXINT,MXKPTS,PME
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER (MXSPE=10,MXTAB=1500,MXINT=150)
C
      LOGICAL DFTBFL,SCC,SRSCC,DFTB3,DAMPXH,LRDFTB
C
      COMMON /DFTB  / DFTBFL,SCC,SRSCC,DFTB3,DAMPXH,LRDFTB
      COMMON /DFTBPR/ ETEMP,DFTBDP(MXSPE*14),DAMPXHE,HUBDER(MXSPE),
     *                ZREF(MXATM),SPNCST(6,MXSPE),SPE(MXATM),NSPE,
     *                MAXANG(MXATM),ISPE(MXATM),IND(MXATM+1),IDFTBD,
     *                PARAMDIR
      COMMON /DFTBRE/ REPCUT(MXSPE,MXSPE),NUMREP(MXSPE,MXSPE),
     *                LREPINTV(MXSPE,MXSPE),LREPSHORT(MXSPE,MXSPE),
     *                LREPCOEFF(MXSPE,MXSPE)
      COMMON /DFTBSK/ SKDIM(MXSPE,MXSPE),SKSPIN(MXSPE),QREFL(3,MXSPE),
     *                HUBBL(3,MXSPE),QREF(MXSPE),HUBB(MXSPE),SKCUT2,
     *                NEEDSK,LSKHTAB(MXSPE,MXSPE),LSKSTAB(MXSPE,MXSPE),
     *                LSKGRID(MXSPE,MXSPE),LSKSELF(MXSPE)
      COMMON /FMCOM / X(1)
      COMMON /IOFILE/ IR,IW,IP,IJK,IJKT,IDAF,NAV,IODA(950)
C
C     RETURN IMMEDIATELY IF NOT DFTB CALCULATION
      IF (.NOT.DFTBFL) RETURN
C
      IF (MODE.EQ.2) THEN
        CALL RETFM(NEEDSK)
        RETURN
      END IF
C
      CALL VALFM(LOADFM)
      CALL GOTFM(NGOTMX)
      !! FOR SLATER-KOSTER PARAMETERS
      LSKHTAB(1,1)   = LOADFM         + 1
      LSKSTAB(1,1)   = LSKHTAB(1,1)   + 10*NSPE*NSPE*MXTAB
      LSKGRID(1,1)   = LSKSTAB(1,1)   + 10*NSPE*NSPE*MXTAB
      LSKSELF(1)     = LSKGRID(1,1)   +  2*NSPE*NSPE
      LAST           = LSKSELF(1)     +  3*NSPE
C
      !! FOR REPULSIVE PARAMETERS
      LREPINTV(1,1)  = LAST
      LREPSHORT(1,1) = LREPINTV(1,1)  + 2*NSPE*NSPE*MXINT
      LREPCOEFF(1,1) = LREPSHORT(1,1) + 3*NSPE*NSPE
      LAST           = LREPCOEFF(1,1) + 6*NSPE*NSPE*MXINT
C
      !! FOR BONDING LIST
      LINDBND        = LAST
      IF (MAXINT.NE.0) THEN
        LAST         = LINDBND        + 3*(MAXINT+1)*MXATM
      END IF
C
      !! FOR PBC
      IF (PERIOD) THEN
        LVKPTS       = LAST
        LVKWGHT      = LVKPTS         + 3*MXKPTS
        LAST         = LVKWGHT        +   MXKPTS
        LBSQ         = LAST
        IF (PME) THEN
          LAST       = LBSQ           +   256 + 1
        END IF
      END IF
C
      !! Make pointer for INDBND and PBC
      CALL DFTB_PBCPTR1(X(LINDBND),X(LVKPTS),X(LVKWGHT),X(LBSQ),MXATM)
C
      NEEDSK = LAST - LOADFM - 1
      CALL GETFM(NEEDSK)
C
C     ----- SET MEMORY LOCATION FOR EACH (PAIR-)SPECIES -----
C
      DO I = 1, NSPE
        DO J = 1, NSPE
          LSKHTAB(J,I)   = LSKHTAB(1,1)   + 10*MXTAB*((J-1)+NSPE*(I-1))
          LSKSTAB(J,I)   = LSKSTAB(1,1)   + 10*MXTAB*((J-1)+NSPE*(I-1))
          LSKGRID(J,I)   = LSKGRID(1,1)   +  2*      ((J-1)+NSPE*(I-1))
          LREPINTV(J,I)  = LREPINTV(1,1)  +  2*MXINT*((J-1)+NSPE*(I-1))
          LREPSHORT(J,I) = LREPSHORT(1,1) +  3*      ((J-1)+NSPE*(I-1))
          LREPCOEFF(J,I) = LREPCOEFF(1,1) +  6*MXINT*((J-1)+NSPE*(I-1))
        END DO
        LSKSELF(I) = LSKSELF(1) + 3*(I-1)
      END DO
C
      IF (MODE.EQ.0) THEN
        CALL VCLR(X(LSKHTAB(1,1)),1,NEEDSK)
      ELSE IF (MODE.EQ.1) THEN
C       RECOVER SK PARAMETERS FROM THE DICTIONARY FILE
        CALL DAREAD(IDAF,IODA,X(LSKHTAB(1,1)),NEEDSK,568,0)
      END IF
C
      RETURN
C
      END SUBROUTINE DFTB_SKMEM
C
C-----------------------------------------------------------------------
C*MODULE DFTBSK    *DECK DFTB_PBCPTR1
C>
C>    @brief   DFTB/PBC
C>
C>    @details Call DFTB_PBCPTR
C>
C>    @author Yoshio Nishimoto
C>    - Oct, 2019- Subroutine written
C>
      SUBROUTINE DFTB_PBCPTR1(INDBND,VKPTS,VKWGHT,BSQ,MXATM)
      USE DFTBPB_MOD,ONLY: DFTB_PBCPTR
      IMPLICIT NONE
      INTEGER INDBND(*),MXATM
      DOUBLE PRECISION VKPTS(*),VKWGHT(*),BSQ(*)
      CALL DFTB_PBCPTR(INDBND,VKPTS,VKWGHT,BSQ,MXATM)
      END SUBROUTINE DFTB_PBCPTR1
C
C-----------------------------------------------------------------------
C*MODULE DFTBSK    *DECK DFTB_TRMAT
C>
C>    @brief   FMO-DFTB analytic gradient
C>
C>    @details Rotate orbitals in 5D (spherical coordinate).
C>             Not applicable to F and G functions. The rotations are
C>             determined by the solutions of the Schrodinger equation
C>             for a hydrogen-like atom:
C>             http://en.wikipedia.org/wiki/Table_of_spherical_harmonics
C>
C>    @author Yoshio Nishimoto
C>    - Feb, 2018- Subroutine written
C>
C
      SUBROUTINE DFTB_TRMAT
      use mx_limits, only: mxsh,mxatm
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      LOGICAL GOPARR,DSKWRK,MASWRK
C
      DIMENSION PNAME(3),DNAME(6),FNAME(10)
C
C
      COMMON /IOFILE/ IR,IW,IP,IJK,IJKT,IDAF,NAV,IODA(950)
      COMMON /OUTPUT/ NPRINT,ITOL,ICUT,NORMF,NORMP,NOPK
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /SYMSPD/ PTR(3,144),DTR(6,288),FTR(10,480),GTR(15,720)
      COMMON /SYMTRY/ MAPSHL(MXSH,48),MAPCTR(MXATM,48),
     *                TT(432),INVT(48),NT
      COMMON /TRANSF/ XSMAL,YSMAL,ZSMAL,XNEW,YNEW,ZNEW,XP,YP,ZP
      COMMON /FRAME / U1,U2,U3,V1,V2,V3,W1,W2,W3,X0,Y0,Z0
C
      PARAMETER (ONE=1.0D+00,SQRT3=1.73205080756888D+00,
     *           SQRT5=2.23606797749979D+00,PI=3.14159265358979312D+00)
C
      DATA PNAME /4HX   ,4HY   ,4HZ   /
      DATA DNAME /4HXX  ,4HYY  ,4HZZ  ,4HXY  ,4HXZ  ,4HYZ  /
      DATA FNAME /4HXXX ,4HYYY ,4HZZZ ,4HXXY ,4HXXZ ,4HYYX ,
     *            4HYYZ ,4HZZX ,4HZZY ,4HXYZ /
C
C     ----- CALCULATE TRANSFORMS OF P, D, F, AND G FUNCTIONS
C           FOR ALL SYMETRY OPERATIONS.
C
      XS = ONE
      YS = ONE
      ZS = ONE
      X = X0+ONE
      Y = Y0
      Z = Z0
      CALL LOCAL(X,Y,Z,XS,YS,ZS)
      XSMAL = XS
      YSMAL = YS
      ZSMAL = ZS
      DO 100 IT = 1,NT
      NN = 9*(IT-1)
      CALL TRANS(NN)
      CALL ROT
      N = 3*(IT-1)
      PTR(1,N+1) = XP-X0
      PTR(2,N+1) = YP-Y0
      PTR(3,N+1) = ZP-Z0
  100 CONTINUE
      X = X0
      Y = Y0+ONE
      Z = Z0
      CALL LOCAL(X,Y,Z,XS,YS,ZS)
      XSMAL = XS
      YSMAL = YS
      ZSMAL = ZS
      DO 120 IT = 1,NT
      NN = 9*(IT-1)
      CALL TRANS(NN)
      CALL ROT
      N = 3*(IT-1)
      PTR(1,N+2) = XP-X0
      PTR(2,N+2) = YP-Y0
      PTR(3,N+2) = ZP-Z0
  120 CONTINUE
      X = X0
      Y = Y0
      Z = Z0+ONE
      CALL LOCAL(X,Y,Z,XS,YS,ZS)
      XSMAL = XS
      YSMAL = YS
      ZSMAL = ZS
      DO 140 IT = 1,NT
      NN = 9*(IT-1)
      CALL TRANS(NN)
      CALL ROT
      N = 3*(IT-1)
      PTR(1,N+3) = XP-X0
      PTR(2,N+3) = YP-Y0
      PTR(3,N+3) = ZP-Z0
  140 CONTINUE
      DO 640 IT = 1,NT
      NP = 3*(IT-1)
      ND = 6*(IT-1)
      NF = 10*(IT-1)
      DO 300 I = 1,5
      GO TO (160,180,200,220,240,260) ,I
  160 J = 1
      K = 2
      GO TO 280
C
  180 J = 2
      K = 3
      GO TO 280
C
  200 J = 3
      K = 1
      GO TO 280
C
  220 J = 1
      K = 2
      GO TO 280
C
  240 J = 3
      K = 3
      GO TO 280
C
  260 J = 2
      K = 3
C 280 DTR(1,ND+I) = PTR(1,NP+J)*PTR(1,NP+K)
C     DTR(2,ND+I) = PTR(2,NP+J)*PTR(2,NP+K)
C     DTR(3,ND+I) = PTR(3,NP+J)*PTR(3,NP+K)
C     DTR(4,ND+I) = PTR(1,NP+J)*PTR(2,NP+K) +PTR(2,NP+J)*PTR(1,NP+K)
C     DTR(5,ND+I) = PTR(1,NP+J)*PTR(3,NP+K) +PTR(3,NP+J)*PTR(1,NP+K)
C     DTR(6,ND+I) = PTR(2,NP+J)*PTR(3,NP+K) +PTR(3,NP+J)*PTR(2,NP+K)
  280 DTR(1,ND+I) = PTR(1,NP+J)*PTR(1,NP+K)*SQRT(3.0D+00/(2.0D+01*PI))
      DTR(2,ND+I) = PTR(2,NP+J)*PTR(2,NP+K)*SQRT(3.0D+00/(2.0D+01*PI))
      DTR(3,ND+I) = PTR(3,NP+J)*PTR(3,NP+K)*SQRT(3.0D+00/(2.0D+01*PI))
      DTR(4,ND+I) = PTR(1,NP+J)*PTR(1,NP+K) -PTR(2,NP+J)*PTR(2,NP+K)
      DTR(4,ND+I) = DTR(4,ND+I)*SQRT(3.0D+00/(5.0D+00*PI))
      DTR(5,ND+I) = 2.0D+00*PTR(3,NP+J)*PTR(3,NP+K)
     *            - PTR(1,NP+J)*PTR(1,NP+K)
     *            - PTR(2,NP+J)*PTR(2,NP+K)
      DTR(5,ND+I) = DTR(5,ND+I)*SQRT(9.0D+00/(5.0D+00*PI))
  300 CONTINUE
      DO 440 I=1,10
      GO TO (320,330,340,350,360,370,380,390,400,410),I
  320 J=1
      K=1
      GO TO 420
  330 J=2
      K=2
      GO TO 420
  340 J=3
      K=3
      GO TO 420
  350 J=1
      K=2
      GO TO 420
  360 J=1
      K=3
      GO TO 420
  370 J=2
      K=1
      GO TO 420
  380 J=2
      K=3
      GO TO 420
  390 J=3
      K=1
      GO TO 420
  400 J=3
      K=2
      GO TO 420
  410 J=4
      K=3
  420 FTR(1,NF+I)=DTR(1,ND+J)*PTR(1,NP+K)
      FTR(2,NF+I)=DTR(2,ND+J)*PTR(2,NP+K)
      FTR(3,NF+I)=DTR(3,ND+J)*PTR(3,NP+K)
      FTR(4,NF+I)=DTR(1,ND+J)*PTR(2,NP+K)
     1           +DTR(4,ND+J)*PTR(1,NP+K)
      FTR(5,NF+I)=DTR(1,ND+J)*PTR(3,NP+K)
     1           +DTR(5,ND+J)*PTR(1,NP+K)
      FTR(6,NF+I)=DTR(2,ND+J)*PTR(1,NP+K)
     1           +DTR(4,ND+J)*PTR(2,NP+K)
      FTR(7,NF+I)=DTR(2,ND+J)*PTR(3,NP+K)
     1           +DTR(6,ND+J)*PTR(2,NP+K)
      FTR(8,NF+I)=DTR(3,ND+J)*PTR(1,NP+K)
     1           +DTR(5,ND+J)*PTR(3,NP+K)
      FTR(9,NF+I)=DTR(3,ND+J)*PTR(2,NP+K)
     1           +DTR(6,ND+J)*PTR(3,NP+K)
      FTR(10,NF+I)=DTR(4,ND+J)*PTR(3,NP+K)
     1           + DTR(5,ND+J)*PTR(2,NP+K)
     2           + DTR(6,ND+J)*PTR(1,NP+K)
  440 CONTINUE
  640 CONTINUE
      IF (NORMF .EQ. 1 .AND. NORMP .EQ. 1) GO TO 750
      DO 740 IT = 1,NT
      ND = 6*(IT-1)
      NF = 10*(IT-1)
      DO 660 I = 1,6
      IF (I .GT. 3) GO TO 650
C     DTR(4,ND+I) = DTR(4,ND+I)/SQRT3
C     DTR(5,ND+I) = DTR(5,ND+I)/SQRT3
C     DTR(6,ND+I) = DTR(6,ND+I)/SQRT3
      GO TO 660
C
  650 CONTINUE
C 650 DTR(1,ND+I) = DTR(1,ND+I)*SQRT3
C     DTR(2,ND+I) = DTR(2,ND+I)*SQRT3
C     DTR(3,ND+I) = DTR(3,ND+I)*SQRT3
  660 CONTINUE
      DO 690 I=1,10
      IF(I.GT.3) GO TO 670
      FTR(4,NF+I)=FTR(4,NF+I)/SQRT5
      FTR(5,NF+I)=FTR(5,NF+I)/SQRT5
      FTR(6,NF+I)=FTR(6,NF+I)/SQRT5
      FTR(7,NF+I)=FTR(7,NF+I)/SQRT5
      FTR(8,NF+I)=FTR(8,NF+I)/SQRT5
      FTR(9,NF+I)=FTR(9,NF+I)/SQRT5
      FTR(10,NF+I)=FTR(10,NF+I)/(SQRT5*SQRT3)
      GO TO 690
  670 IF(I.GT.9) GO TO 680
      FTR(1,NF+I)=FTR(1,NF+I)*SQRT5
      FTR(2,NF+I)=FTR(2,NF+I)*SQRT5
      FTR(3,NF+I)=FTR(3,NF+I)*SQRT5
      FTR(10,NF+I)=FTR(10,NF+I)/SQRT3
      GO TO 690
  680 FTR(1,NF+I)=FTR(1,NF+I)*SQRT5*SQRT3
      FTR(2,NF+I)=FTR(2,NF+I)*SQRT5*SQRT3
      FTR(3,NF+I)=FTR(3,NF+I)*SQRT5*SQRT3
      FTR(4,NF+I)=FTR(4,NF+I)*SQRT3
      FTR(5,NF+I)=FTR(5,NF+I)*SQRT3
      FTR(6,NF+I)=FTR(6,NF+I)*SQRT3
      FTR(7,NF+I)=FTR(7,NF+I)*SQRT3
      FTR(8,NF+I)=FTR(8,NF+I)*SQRT3
      FTR(9,NF+I)=FTR(9,NF+I)*SQRT3
  690 CONTINUE
  740 CONTINUE
  750 CONTINUE
C
C     ----- PRINT MATRICES IF NPRINT.EQ.5 -----
C
      IF (NPRINT .NE. 5) GO TO 860
      IF (MASWRK) THEN
      WRITE (IW,9028)
      DO 840 IT = 1,NT
      WRITE (IW,9088)
      WRITE (IW,9108) IT
      NP = 3*(IT-1)
      WRITE (IW,9048) (PNAME(J),J = 1,3)
      WRITE (IW,9128)
      DO 760 I = 1,3
  760 WRITE (IW,9068) PNAME(I),(PTR(I,NP+J),J = 1,3)
      WRITE (IW,9008)
      ND = 6*(IT-1)
      WRITE (IW,9048) (DNAME(J),J = 1,6)
      WRITE (IW,9128)
      DO 780 I = 1,6
  780 WRITE (IW,9068) DNAME(I),(DTR(I,ND+J),J = 1,6)
      WRITE (IW,9008)
      NF = 10*(IT-1)
      WRITE (IW,9048) (FNAME(J),J = 1,10)
      WRITE (IW,9128)
      DO 800 I = 1,10
  800 WRITE (IW,9068) FNAME(I),(FTR(I,ND+J),J = 1,10)
      WRITE (IW,9008)
  840 CONTINUE
      END IF
  860 CONTINUE
      RETURN
C
 9008 FORMAT(//)
 9028 FORMAT(/,38H TRANSFORMATION OF THE BASIS FUNCTIONS,/)
 9048 FORMAT(7X,15(2X,A4,2X))
 9068 FORMAT(2X,A4,1X,15F8.4)
 9088 FORMAT(1H1)
 9108 FORMAT(/,21X,21HTRANSFORMATION NUMBER,I4,/)
 9128 FORMAT(/)
      END
