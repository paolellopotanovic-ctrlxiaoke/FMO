C 19 Oct 12 - MWS - synchronize FRGINF common
C 21 JAN 11 - DK  - EFP2EI: ADJUST OVERLY LONG PARAMETER STATEMENT
C 14 OCT 09 - DMK - 3-CENTER QM, 1- CENTER EFP 2E- INTS FOR EFP2 EXREP 
C
C*MODULE EFTEI   *DECK EFP2EI
C>    @brief    QM-EFP 2EI driver 
C>
C>    @author   Dan Kemp
C>            - 2010?
C>    @details  calculate QM-EFP two-electron integrals 
C>    @date     Jul 2015   Peng Xu
C>            - modified to allow multiple EFP fragments
C>    @date     Oct 2018   Peng Xu and Tosaporn Sattasathuchana
C>            - create two modes for two types of QM/EFP TEI
C>              and remove overlap screening
C>    @param SMATL: overlap between QM and all EFP fragments
C>    @param NUMEFT: total number of EFP AOs
C>    @param XINTB: largest exchange integrals of each shell between QM and EFP 
C>    @param XINTA: largest exchange integrals of each shell (QM only)
C>    @param NSHELT: total number of EFP shells
C>    @param MODE: different type of TEI involve in QM/EFP EXREP
C>                 MODE=1 ALL 4 center are in QM region.
C>                 MODE=2 two-center on QM, two-center on EFP
C>                 MODE=3 three-center on QM, One-center on EFP
C>    @see   subroutine TWOEI in int2a.src
C>    
      SUBROUTINE EFP2EI(MODE,TYPSCF,DIRSCF,DIRNLO,DIRTRF,DIRCIS,
     *                  SCHWRZ,NINT,NSCHWZ,L1,L2,
     *                  BUFP,BUFK,IX,NINTMX,
     *                  GHONDO,MAXG,
     *                  IA,DA,FA,DB,FB,DNLO,FNLO,NFLMAT,
     *                  SMATL,NUMEFT,XINTB,XINTA,NSHELT)
C
C     DRIVER FOR AB INITIO/EFP2 2E- INTEGRALS.  BASED ON TWOEI
C
      USE MX_LIMITS,ONLY:MXGEFP,MXDFG,MXFGPT,MXFRG,MXDPPT,mxsh,mxgtot
      USE comm_FRGTYP
      USE comm_EFPBAS
      USE comm_FRGINF
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      LOGICAL DIRSCF,DIRNLO,DIRTRF,DIRCIS,CMBDIR
      LOGICAL SCHWRZ,GOPARR,DSKWRK,MASWRK,DLB,SLB,SCHSKP
      LOGICAL PK,PANDK,NOTPK,BLOCK,GPSAVE,OUT
C
      DIMENSION BUFP(NINTMX),BUFK(NINTMX),IX(*),
     *          GHONDO(MAXG),IA(L1),DA(L2),FA(L2),DB(L2),FB(L2),
     *          DNLO(L1,L1),FNLO(L1,L1)
      DIMENSION SMATL(L1,NUMEFT),XINTB(NSHELT,NSHELL),
     *          XINTA((NSHELL*NSHELL+NSHELL)/2),XINTASQ(NSHELL,NSHELL)
C
      PARAMETER (HUND=100.0D+00)
      PARAMETER (ZERO=0.0D+00)
C
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON /NSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),
     *                CF(MXGTOT),CG(MXGTOT),CH(MXGTOT),CI(MXGTOT),
     *                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),
     *                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
      COMMON /OUTPUT/ NPRINT,ITOL,ICUT,NORMF,NORMP,NOPK
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /PKFIL / PK,PANDK,BLOCK
      COMMON /RESTAR/ TIMLIM,IREST,NREC,INTLOC,IST,JST,KST,LST
      COMMON /SHLEXC/ NORGSH(3),NORGSP(3),IEXCH,NANGM,NGTH(4)
      COMMON /SHLG70/ ISH,JSH,KSH,LSH,IJKLXX(4)
      COMMON /SHLNOS/ QQ4,LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,
     *                MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL,
     *                NIJ,IJ,KL,IJKL
      COMMON /SHLT  / TOL,CUTOFF,ICOUNT,OUT
C
      STOL=CUTOFF/HUND
      ICNTSV = ICOUNT
      NRECSV = NREC 
      INTLSV = INTLOC
      ICOUNT = 1
      NREC   = 0
      INTLOC = 1
      IF(MODE.EQ.2) THEN
C OPEN FILE 28 FOR (IJ|IJ) EFP-QM: EFAI2IJ
      NFT=28
      CALL SEQOPN(NFT,'MLTPL','UNKNOWN',.FALSE.,'UNFORMATTED')
      ELSE IF(MODE.EQ.3) THEN
C OPEN FILE 39 FOR (IJ|KL) EFP_I-QM_JKL: QMMMTEI
      NFT=39
      CALL SEQOPN(NFT,'QMMMTEI','UNKNOWN',.FALSE.,'UNFORMATTED')
      ENDIF
      TIM = ZERO
      CALL TSECND(TIM)
C
C   P. XU :read in QM exchange integral XINTA from disk
C           for Schwartz inequality screening
C           transform XINTA to a square matrix XINTASQ

      IF(SCHWRZ) THEN
      NSH2 = (NSHELL*NSHELL+NSHELL)/2
      CALL DAREAD(IDAF,IODA,XINTA,NSH2,54,0)
      CALL CPYTSQ(XINTA,XINTASQ,NSHELL,1)
      END IF
C
      CMBDIR= DIRSCF .OR. DIRNLO .OR. DIRTRF .OR. DIRCIS
C
C        THE OLD FASHIONED PARALLEL INTEGRAL TRANSFORMATIONS DO NOT
C        ALLOW THE AO INTEGRAL WORK TO BE RUN IN PARALLEL.  THE MODERN
C        DISTRIBUTED MEMORY TRANSFORMATIONS DO NOT PASS THROUGH HERE.
C
      ICONT=  0
      GPSAVE = GOPARR
      IF(DIRTRF) GOPARR=.FALSE.
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
      NT1 = 1
C
      LMAX=0
      DO MJ = 1,NFRG
      JM = ISET(MJ)
      CALL BASEFP(LMAX2,JM)
      IF (LMAX2 .GT. LMAX) LMAX=LMAX2
      END DO
                    NANGM =  4
      IF(LMAX.EQ.2) NANGM =  6
      IF(LMAX.EQ.3) NANGM = 10
      IF(LMAX.EQ.4) NANGM = 15
      IF(LMAX.EQ.5) NANGM = 21
      IF(LMAX.EQ.6) NANGM = 28
      NGTH(4) = 1
      NGTH(3) = NGTH(4) * NANGM
      NGTH(2) = NGTH(3) * NANGM
      NGTH(1) = NGTH(2) * NANGM
      IF(NOPK.EQ.0) THEN
         NORGSH(1) = 0
         NORGSH(2) = NORGSH(1) + NANGM**4
         NORGSH(3) = NORGSH(2) + NANGM**4
         NORGSP(1) = 0
         NORGSP(2) = 256
         NORGSP(3) = 512
      ELSE
         DO I=1,3
            NORGSH(I) = 0
            NORGSP(I) = 0
         ENDDO
      END IF
C
      NOTPK = .NOT.PK
      NINT  = 0
      NSCHWZ= 0
      SCHSKP=.FALSE.
C
C        NOW WE ARE READY TO LOOP OVER ALL NSHELL**4 SHELL QUARTETS
C
C     ----- I SHELL -----
C
C
      JAO  = 0
      JNAT = 0
      ISHLEF = 0

      DO 1000 MJ = 1,NFRG
      JM = ISET(MJ)
      DO 920 II = 1,NSHLEF(JM)
C
C     ----- CHECK CPU TIME -----
C
      CALL TSECND(TIM)
      IF(TIM.GE.TIMLIM) THEN
C        NOTHING CAN BE DONE FOR IN-CORE INTEGRALS: JUST FORGET THEM
         IF(.NOT.CMBDIR)
     *      CALL FINAL(0,II,1,1,1,PANDK,BUFP,BUFK,IX,NINTMX)
         IF(MASWRK) WRITE(IW,9030)
         CALL ABRT
      END IF
C
C
C     ----- J SHELL -----
C
C     J0 = JST
      DO 900 JJ = 1,NSHELL
C
C     ----- GO PARALLEL! -----
C
      IF (DLB) THEN
         MINE = MINE + 1
         IF (MINE.GT.NEXT) CALL DDI_DLBNEXT(NEXT)
         IF (NEXT.NE.MINE) GO TO 900
      END IF
C
C     ----- K SHELL -----
C
C MODE 2 FOR (ij|ij) EFP-QM
C MODE 3 FOR (ij|kl) EFP-QM, k is qm
      IF(MODE.EQ.2) THEN
          KNSHELL=NSHLEF(JM)
      ELSE
        KNSHELL=NSHELL
      ENDIF
      DO 880 KK = 1,KNSHELL
C
C     ----- GO PARALLEL! -----
C
      IF(SLB) THEN
         IPCOUNT = IPCOUNT + 1
         IF (MOD(IPCOUNT,NPROC).NE.0) GO TO 880
      END IF
C
C     ----- L SHELL ----
C
      IF(MODE.EQ.2) THEN
          LNSHELL=NSHELL
      ELSE
        LNSHELL=KK
      ENDIF
      DO 860 LL = 1,LNSHELL
      N4=1
C
C         THE LOOP STRUCTURE IN THIS ROUTINE IS DESIGNED TO FACILITATE
C         SUPERMATRIX CONSTRUCTION BY HAVING UP TO THREE "EXCHANGED"
C         QUARTETS AVAILABLE AT ONCE.  THE LOOP STRUCTURE TO GENERATE
C         A MORE NORMAL CANONICAL ORDERING OF THE QUARTETS HITS THE
C         SAME QUARTETS IN A SLIGHTLY DIFFERENT ORDER, BUT BOTH LOOPS
C         WILL DO EXACTLY THE SAME QUARTETS.
C
C             CANONICAL                      SUPERMATRIX
C         DO ISH=1,NSHELL                 DO II=1,NSHELL
C           DO JSH=1,ISH                    DO JJ=1,II
C             IJSH = IA(ISH)+JSH
C             DO KSH=1,ISH                    DO KK=1,JJ
C               DO LSH=1,KSH                    DO LL=1,KK
C                 KLSH=IA(KSH)+LSH
C                 IF(IJSH.LT.KLSH),               [II JJ|KK LL],
C                    CYCLE KSH LOOP               [II KK|JJ LL],
C                 [ISH JSH|KSH LSH]               [II LL|JJ KK]
C               ENDDO                           ENDDO
C             ENDDO                           ENDDO
C           ENDDO                           ENDDO
C         ENDDO                           ENDDO
C
      Q4 = NT1
      Q4 = Q4 / N4
C
C     ----- (II,JJ//KK,LL) -----
C
      IEXCH = 1
      ISH = II
      JSH = JJ
      KSH = KK
      LSH = LL
      QQ4 = Q4
      GO TO 780
C
C        ----- COMPUTE TWO-ELECTRON INTEGRALS ----
C
  780 CONTINUE
      IF(PK .AND. IEXCH.EQ.1) CALL ZPKOUT(ISH,JSH,KSH,LSH,GHONDO,
     *                                    .FALSE.,.FALSE.,
     *                                    .FALSE.,.FALSE.)
C
C     APPLY THE SCHWARZ INEQUALITY TO SCREEN OUT SMALL INTEGRALS,
C      (II,JJ//KK,LL) .LE.  SQRT( (II,JJ//II,JJ)*(KK,LL//KK,LL) )
C     SEE, FOR EXAMPLE, J.L.WHITTEN, J.CHEM.PHYS. 58,4496-4501(1973)
C
      SCHSKP = .FALSE.
      IF(SCHWRZ) THEN
C     P. XU: Schwarz inequality screening
       IF(MODE.EQ.2) THEN
         TEST = XINTB(II+ISHLEF,JJ)*XINTB(KK+ISHLEF,LL)
       ELSE IF(MODE.EQ.3) THEN
         TEST = XINTB(II+ISHLEF,JJ)*XINTASQ(KK,LL)
       ENDIF
         SCHSKP = TEST .LT. CUTOFF
      END IF
      IF(SCHSKP) GO TO 820
C
      CALL SHELLQUART_EFP(MODE,ISH,JSH,KSH,LSH,GHONDO,JM,JNAT)
C
C        USE THE INTEGRALS JUST FORMED.  AT MOST, 1 OF THESE IS CALLED
C
           IF(DIRSCF) THEN
         CALL DIRFCK(TYPSCF,IA,DA,FA,DB,FB,GHONDO,L2,NINT,NFLMAT)
      ELSE IF(DIRCIS) THEN
         CALL DRFCIS(DA,FA,DB,FB,GHONDO,L1,NINT,NFLMAT)
      ELSE IF(DIRNLO) THEN
         CALL DFCKNS(DNLO,FNLO,GHONDO,L1,NINT,NFLMAT)
      ELSE IF(DIRTRF) THEN
         CALL DIRTRN(BUFP,IX,NINTMX,GHONDO,ICONT,NINT)
      ELSE
         IF(NOTPK) CALL QOUTEFP(MODE,BUFP,IX,NINTMX,GHONDO,JM)
      END IF
C
  820 CONTINUE
C
  860 CONTINUE
  880 CONTINUE
  900 CONTINUE
  920 CONTINUE

C----P. XU: call FINAL_EFP inside fragment loop
C     ----- OUTPUT THE LAST BITS OF INTEGRALS -----
C
      IF(DIRTRF) CALL ONEIDX(BUFP,IX,ICONT)
      IF(.NOT.CMBDIR) THEN
         CALL FINAL_EFP(1,II,II,II,II,PANDK,BUFP,BUFK,IX,NINTMX,MODE)
      END IF
C
C RESET ICOUNT=1 AND NREC=0 FOR NEW FRAGMENT
      ICOUNT= 1
      NREC=   0

      JAO = JAO + NUMEF(JM)
      JNAT = JNAT + NATEF(JM)
      ISHLEF = ISHLEF + NSHLEF(JM)

 1000 CONTINUE      
C
      IF(DLB) CALL DDI_DLBRESET
C
      GOPARR = GPSAVE
      ICOUNT = ICNTSV
      NREC   = NRECSV
      INTLOC = INTLSV
      RETURN
C
 9030 FORMAT(//1X,'*** THIS JOB HAS EXHAUSTED ITS CPU TIME ***'/
     *         1X,'     (WHILE COMPUTING 2E- INTEGRALS)'///)
      END
C
C*MODULE EFTEI   *DECK BASEFP
C>    @brief  highest angular momentum number in the basis set 
C>
C>    @author Dan Kemp
C>            - 2010?
C>    @details determine the highest angular momentum present 
C>             in the basis of both QM and EFPs
C>             KTYPE=1,2,3,4,5 means S,P(L),D,F,G function  
C>    @date    Jul 2015  Peng Xu
C>             -bug fix 
C>    @param   LMAX: max value of angular momentum
C>    @param   JM: index for the type of EFP potentials
C>    @see     BASCHK in int2a.src
C
      SUBROUTINE BASEFP(LMAX,JM)
C
      USE MX_LIMITS,ONLY:MXDFG,MXFGPT,MXSHEF,MXGEFP,mxsh,mxgtot
      USE comm_EFPBAS
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C
      COMMON /NSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),
     *                CF(MXGTOT),CG(MXGTOT),CH(MXGTOT),CI(MXGTOT),
     *                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),
     *                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
C
C
C     RETURN THE HIGHEST ANGULAR MOMENTUM PRESENT IN THE BASIS.
C     NOTE THAT KTYPE=1,2,3,4,5 MEANS S, P(L), D, F, G FUNCTION.
C     SCAN EFP2 BASIS AS WELL
C
      KANG = 0
      DO 100 N=1,NSHELL
          IF(KTYPE(N).GT.KANG) KANG = KTYPE(N)
  100 CONTINUE
      DO 200 N=1,NSHLEF(JM)
          IF(KTYPEF(N,JM).GT.KANG) KANG = KTYPEF(N,JM)
  200 CONTINUE
      LMAX = KANG-1
      RETURN
      END
C
C*MODULE EFTEI   *DECK FINAL_EFP
C>    @brief  output the last bit of 2EI
C>
C>    @author Dan Kemp
C>            - 2010?
C>    @details  output the last bit of QM-EFP 2EI to disk
C>    @date   Jul 2015  Peng Xu
C>            -bug fix 
C>    @date   Oct 2018  Peng Xu and Tosaporn Sattasathuchana
C>            - adding MODE for the differenty types of QM/EFP TEI
C>            - the initialization of ICOUNT and NREC is move to EFP2EI
C>              subroutine
C>    @param MODE: different type of TEI involve in QM/EFP EXREP
C>                 MODE=1 ALL 4 center are in QM region.
C>                 MODE=2 two-center on QM, two-center on EFP
C>                 MODE=3 three-center on QM, One-center on EFP
C>    @see    subroutine FINAL in int2a.src

      SUBROUTINE FINAL_EFP(INDEX,II,JJ,KK,LL,PANDK,BUFP,BUFK,IX,NINTMX,
     &                     MODE)
      use mx_limits, only: mxsh,mxgtot
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     LOGICAL PANDK,OUT,GOPARR,DSKWRK,MASWRK,LTRMST
      LOGICAL PANDK,OUT,GOPARR,DSKWRK,MASWRK
C     LOGICAL LCFLAG,LRINT
C
      DIMENSION BUFP(NINTMX),BUFK(NINTMX),IX(*)
C
C
C     COMMON /ELGPMT/ NELONG,NATM,NASPIN,NCT,NBNDAB,NTMLB,IPRI,LDOS
C     COMMON /ELGTRM/ LTRMST,NFLTRM,NRCTRM,NPSTRM,NHTSHL
      COMMON /FMOINF/ NFG,NLAYER,NATFMO,NBDFG,NAOTYP,NBODY,NSEGM
      COMMON /INT2IC/ NINTIC,ININTIC,NXXIC,LBUFPIC,LIXIC,LABSIX,NINTIX
      COMMON /INTPR / Q(2),V(2),JC,N1(2),J1(2),J2(2),J3(2),J4(2)
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON /NSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),
     *                CF(MXGTOT),CG(MXGTOT),CH(MXGTOT),CI(MXGTOT),
     *                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),
     *                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
      COMMON /OUTPUT/ NPRINT,ITOL,ICUT,NORMF,NORMP,NOPK
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /RESTAR/ TIMLIM,IREST,NREC,INTLOC,IST,JST,KST,LST
      COMMON /SHLT  / TOL,CUTOFF,ICOUNT,OUT
C
      IF(MODE.EQ.2) IT=28
      IF(MODE.EQ.3) IT=39
C     ISSAVE=IS
C     IF(LRINT)IS=LRFILE
C
      IF(OUT .AND. JC.GT.0)
     * WRITE(IW,9088) (J1(M),J2(M),J3(M),J4(M),Q(M),V(M),M=1,JC)
C
      IF (INDEX.EQ.1) GO TO 140
C
C         RAN OUT OF TIME, FLUSH PARTIAL BUFFER, PRINT RESTART DATA
C
      IREST = 1
      IST = II
      JST = JJ
      KST = KK
      LST = LL+1
      IF(LST.LE.KK) GO TO 120
      LST = 1
      KST = KK+1
      IF(KST.LE.JJ) GO TO 120
      KST = 1
      JST = JJ+1
      IF(JST.LE.II) GO TO 120
      JST = 1
      IST = II+1
      IF(IST.GT.NSHELL) GO TO 140
C
  120 CONTINUE
CNB   DO SOMETHING FOR IN CORE
      NXX = ICOUNT-1
      IF(.NOT.PANDK) CALL PWRIT(IT,BUFP,IX,NXX,NINTMX)
      IF(     PANDK) CALL PKWRIT(IT,BUFP,BUFK,IX,NXX,NINTMX)
      NINT = NINTMX*(NREC-1)+ICOUNT-1
C     PRINT *, 'NINT IN FINAL', NINT
      IF (MASWRK) THEN
         WRITE(IW,9010) NINT,NREC,IT
         WRITE(IW,9020) NREC,ICOUNT,IST,JST,KST,LST
      END IF
C     IF(LRINT)IS=ISSAVE
      RETURN
C
C        DONE WITH INTEGRALS, WRITE LAST BUFFER, PRINT STATISTICS
C
  140 CONTINUE
      IREST = 0
      IST = 1
      JST = 1
      KST = 1
      LST = 1
      NXX = ICOUNT-1
      IF(NXX.GE.NINTIC) THEN
C
C     GE INSTEAD OF GT ENFORCES WRITING A RECORD WITH 0 INTEGRALS
C     (IF NINTIC IS EQUAL TO NXX) TO COMPLY WITH THE INTEGRAL FILE
C     STRUCTURE
C
         NXX=NXX-NINTIC
         NXX = -NXX
C        WRITE(6,*) 'SAVING',NXX
         IF(.NOT.PANDK) CALL PWRIT(IT,BUFP(NINTIC+1),IX(ININTIC+1),
     *                             NXX,NINTMX)
         NREC = NREC + 1
         NINT = NINTMX*(NREC-1)+ICOUNT-1
C        NOTE THAT ICOUNT-1 INCLUDES NINTIC
      ELSE
         NINT = NXX
         PRINT *, 'NINT IN FINAL', NINT
      ENDIF
C     THE # OF INTEGRALS MUST BE STORED ON EACH NODE BEFORE GSUMI.
      IF(NINTIC.NE.0) NXXIC=MIN(NINT,NINTIC)
      NINTMY=NINT
C
C          ELONGATION METHOD INTEGRAL FILES
C
C     IF(NELONG.GT.1.AND.(.NOT.LTRMST)) THEN
C        NRCTRM = NREC
C        NPSTRM = ICOUNT
C     END IF
C
      IF (GOPARR) THEN
         CALL DDI_GSUMI(1056,NINT,1)
         CALL DDI_GSUMI(1057,NREC,1)
      END IF
      IF (MASWRK.AND.(NFG.EQ.0.OR.NPRINT.NE.-5)) THEN
         IF(NINTIC.EQ.0) THEN
C           WRITE(IW,9010) NINT,NREC,IT
         ELSE IF(NINTMY.LE.NINTIC) THEN
C           WRITE(IW,9015) NINT,(NINTMY*1.0D+02)/NINTIC
         ELSE
C           WRITE(IW,9017) NINT,NINTIC,(NINTIC*1.0D+02)/NINTMY,
C    *                     NREC,IT
         ENDIF
      ENDIF
C     IF(LRINT)IS=ISSAVE
      RETURN
C
 9010 FORMAT(1X,'TOTAL NUMBER OF NONZERO TWO-ELECTRON INTEGRALS =',I20/
     *       1X,I10,' INTEGRAL RECORDS WERE STORED ON DISK ',
     *       'FILE',I3,'EFTEI.')
C9015 FORMAT(1X,'TOTAL NUMBER OF NONZERO TWO-ELECTRON INTEGRALS =',I20/
C    *       1X,'ALL INTEGRALS FITTED INTO MEMORY (',F5.1,
C    *          '% OF BUFFER USED)!')
C9017 FORMAT(1X,'TOTAL NUMBER OF NONZERO TWO-ELECTRON INTEGRALS =',I20
C    *      /1X,'ON NODE 0',I13,' INTEGRALS STORED IN MEMORY(',F4.1,'%)'
C    *     ,/I10,' INTEGRAL RECORDS WERE STORED ON DISK FILE',I3,'.')
 9020 FORMAT(/1X,'...... WARNING   .......   WARNING   .......'/
     *        1X,'TIME LIMIT HAS EXPIRED.  THIS JOB MUST BE RESTARTED.'/
     *        1X,'IF YOU SAVED THE INTEGRALS, RESTART WITH IREST=1,'/
     *        1X,'NREC=',I8,' INTLOC=',I6,' IST,JST,KST,LST=',4I6)
 9088 FORMAT(2(4I4,F5.1,F17.9,1X))
      END
C
C*MODULE EFTEI   *DECK IJPRIMEFP
C>    @brief  prepare primitive functions for EFP-AI 2EI calculation 
C>
C>    @author Dan Kemp
C>            - 2010?
C>    @details
C>    @see    subroutine IJPRIM in int2a.src

      SUBROUTINE IJPRIMEFP(DDIJ)
      use mx_limits, only: mxgsh,mxg2
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     LOGICAL IANDJ,KANDL,SAME,OUT,NORM
      LOGICAL OUT,NORM
C
C
      DIMENSION DDIJ(49*MXG2)
C
      COMMON /IJGNRL/ A(MXG2),R(MXG2),X1(MXG2),Y1(MXG2),Z1(MXG2),
     *                IJD(784)
C     COMMON /MISC  / IANDJ,KANDL,SAME
      COMMON /OUTPUT/ NPRINT,ITOL,ICUT,NORMF,NORMP,NOPK
      COMMON /SHLINF/  AG(MXGSH),CSA(MXGSH),CPA(MXGSH),CDA(MXGSH),
     *                CFA(MXGSH),CGA(MXGSH),CHA(MXGSH),CIA(MXGSH),
     *                 BG(MXGSH),CSB(MXGSH),CPB(MXGSH),CDB(MXGSH),
     *                CFB(MXGSH),CGB(MXGSH),CHB(MXGSH),CIB(MXGSH),
     *                 CG(MXGSH),CSC(MXGSH),CPC(MXGSH),CDC(MXGSH),
     *                CFC(MXGSH),CGC(MXGSH),CHC(MXGSH),CIC(MXGSH),
     *                 DG(MXGSH),CSD(MXGSH),CPD(MXGSH),CDD(MXGSH),
     *                CFD(MXGSH),CGD(MXGSH),CHD(MXGSH),CID(MXGSH),
     *                XI,YI,ZI,XJ,YJ,ZJ,RRI,XK,YK,ZK,XL,YL,ZL,RRK,
     *                NGA,NGB,NGC,NGD
      COMMON /SHLNOS/ QQ4,LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,
     *                MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL,
     *                NIJ,IJ,KL,IJKL
      COMMON /SHLT  / TOL,CUTOFF,ICOUNT,OUT
C
      PARAMETER (SQRT3=1.73205080756888D+00)
      PARAMETER (SQRT5=2.23606797749979D+00)
      PARAMETER (SQRT7=2.64575131106459D+00)
      PARAMETER (SQRT9=3.0D+00)
      PARAMETER (SQRT11=3.3166247903553998D+00)
      PARAMETER (ZERO=0.0D+00)
      PARAMETER (ONE=1.0D+00)
C
      NORM = NORMF .NE. 1 .OR. NORMP .NE. 1
      MAX = MAXJ
      N = 0
      NN = 0
      NM = -2**20
      DO 180 I = MINI,MAXI
         GO TO (100,100,120,120,100,120,120,100,120,120,
     1          100,120,120,100,120,120,120,120,120,100,
     1          100,120,120,100,120,120,120,120,120,100,
     1          120,120,100,120,120,
     1          100,120,120,100,120,120,120,120,120,100,
     1          120,120,120,120,120,100,120,120,100,120,
     1          120,
     1          100,120,120,100,120,120,120,120,120,100,
     1          120,120,120,120,120,100,120,120,100,120,
     1          120,100,120,120,120,120,120,100),I
  100    NM = NN
  120    NN = NM
C        END IF
         DO 170 J = MINJ,MAX
            GO TO (140,140,160,160,140,160,160,140,160,160,
     1             140,160,160,140,160,160,160,160,160,140,
     1             140,160,160,140,160,160,160,160,160,140,
     1             160,160,140,160,160,
     1             140,160,160,140,160,160,160,160,160,140,
     1             160,160,160,160,160,140,160,160,140,160,
     1             160,
     1             140,160,160,140,160,160,160,160,160,140,
     1             160,160,160,160,160,140,160,160,140,160,
     1             160,140,160,160,160,160,160,140),J
  140       NN = NN+1
  160       N = N+1
            IJD(N) = NN
  170    CONTINUE
  180 CONTINUE
C
C     ----- I PRIMITIVE
C
      NIJ = 0
      JBMAX = NGB
      DO 540 IA = 1,NGA
         AI = AG(IA)
         ARRI = AI*RRI
         AXI = AI*XI
         AYI = AI*YI
         AZI = AI*ZI
         CSI = CSA(IA)
         CPI = CPA(IA)
         CDI = CDA(IA)
         CFI = CFA(IA)
         CGI = CGA(IA)
         CHI = CHA(IA)
         CII = CIA(IA)
C
C        ----- J PRIMITIVE
C
         DO 520 JB = 1,JBMAX
            AJ = BG(JB)
            AA = AI+AJ
            AAINV = ONE/AA
            DUM = AJ*ARRI*AAINV
            IF (DUM .GT. TOL) GO TO 520
            CSJ = CSB(JB)
            CPJ = CPB(JB)
            CDJ = CDB(JB)
            CFJ = CFB(JB)
            CGJ = CGB(JB)
            CHJ = CHB(JB)
            CIJ = CIB(JB)
            NM = 49*NIJ
            NN = NM
            NIJ = NIJ+1
            R(NIJ) = DUM
            A(NIJ) = AA
            X1(NIJ) = (AXI+AJ*XJ)*AAINV
            Y1(NIJ) = (AYI+AJ*YJ)*AAINV
            Z1(NIJ) = (AZI+AJ*ZJ)*AAINV
C
C           ----- DENSITY FACTOR
C
            DUM1 = ZERO
            DUM2 = ZERO
            DO 420 I = MINI,MAXI
               GO TO (200,220,420,420,240,420,420,260,420,420,
     1                261,420,420,262,420,420,420,420,420,263,
     1                264,420,420,265,420,420,420,420,420,266,
     1                420,420,267,420,420,
     1                268,420,420,269,420,420,420,420,420,270,
     1                420,420,420,420,420,271,420,420,272,420,
     1                420,
     1                273,420,420,274,420,420,420,420,420,275,
     1                420,420,420,420,420,276,420,420,277,420,
     1                420,278,420,420,420,420,420,279),I
  200          DUM1 = CSI*AAINV
               GO TO 280
  220          DUM1 = CPI*AAINV
               GO TO 280
  240          DUM1 = CDI*AAINV
               GO TO 280
  260          IF (NORM) DUM1 = DUM1*SQRT3
               GO TO 280
  261          DUM1 = CFI*AAINV
               GO TO 280
  262          IF (NORM) DUM1 = DUM1*SQRT5
               GO TO 280
  263          IF (NORM) DUM1 = DUM1*SQRT3
               GO TO 280
  264          DUM1 = CGI*AAINV
               GO TO 280
  265          IF (NORM) DUM1 = DUM1*SQRT7
               GO TO 280
  266          IF (NORM) DUM1 = DUM1*SQRT5/SQRT3
               GO TO 280
  267          IF (NORM) DUM1 = DUM1*SQRT3
               GO TO 280
  268          DUM1 = CHI*AAINV
               GO TO 280
  269          IF (NORM) DUM1 = DUM1*SQRT9
               GO TO 280
  270          IF (NORM) DUM1 = DUM1*SQRT7/SQRT3
               GO TO 280
  271          IF (NORM) DUM1 = DUM1*SQRT3
               GO TO 280
  272          IF (NORM) DUM1 = DUM1*SQRT5/SQRT3
               GO TO 280
  273          DUM1 = CII*AAINV
               GO TO 280
  274          IF (NORM) DUM1 = DUM1*SQRT11
               GO TO 280
  275          IF (NORM) DUM1 = DUM1*SQRT3
               GO TO 280
  276          IF (NORM) DUM1 = DUM1*SQRT3
               GO TO 280
  277          IF (NORM) DUM1 = DUM1*SQRT7/(SQRT5*SQRT3)
               GO TO 280
  278          IF (NORM) DUM1 = DUM1*SQRT5
               GO TO 280
  279          IF (NORM) DUM1 = DUM1*SQRT5/SQRT3
C
  280          CONTINUE
C              END IF
               DO 400 J = MINJ,MAX
                  GO TO (300,320,400,400,340,400,400,360,400,400,
     1                   361,400,400,362,400,400,400,400,400,363,
     1                   364,400,400,365,400,400,400,400,400,366,
     1                   400,400,367,400,400,
     1                   368,400,400,369,400,400,400,400,400,370,
     1                   400,400,400,400,400,371,400,400,372,400,
     1                   400,
     1                   373,400,400,374,400,400,400,400,400,375,
     1                   400,400,400,400,400,376,400,400,377,400,
     1                   400,378,400,400,400,400,400,379),J
  300             DUM2 = DUM1*CSJ
                  GO TO 380
  320             DUM2 = DUM1*CPJ
                  GO TO 380
  340             DUM2 = DUM1*CDJ
                  GO TO 380
  360             IF (NORM) DUM2 = DUM2*SQRT3
                  GO TO 380
  361             DUM2 = DUM1*CFJ
                  GO TO 380
  362             IF (NORM) DUM2 = DUM2*SQRT5
                  GO TO 380
  363             IF (NORM) DUM2 = DUM2*SQRT3
                  GO TO 380
  364             DUM2 = DUM1*CGJ
                  GO TO 380
  365             IF (NORM) DUM2 = DUM2*SQRT7
                  GO TO 380
  366             IF (NORM) DUM2 = DUM2*SQRT5/SQRT3
                  GO TO 380
  367             IF (NORM) DUM2 = DUM2*SQRT3
                  GO TO 380
  368             DUM2 = DUM1*CHJ
                  GO TO 380
  369             IF (NORM) DUM2 = DUM2*SQRT9
                  GO TO 380
  370             IF (NORM) DUM2 = DUM2*SQRT7/SQRT3
                  GO TO 380
  371             IF (NORM) DUM2 = DUM2*SQRT3
                  GO TO 380
  372             IF (NORM) DUM2 = DUM2*SQRT5/SQRT3
                  GO TO 380
  373             DUM2 = DUM1*CIJ
                  GO TO 380
  374             IF (NORM) DUM2 = DUM2*SQRT11
                  GO TO 380
  375             IF (NORM) DUM2 = DUM2*SQRT3
                  GO TO 380
  376             IF (NORM) DUM2 = DUM2*SQRT3
                  GO TO 380
  377             IF (NORM) DUM2 = DUM2*SQRT7/(SQRT5*SQRT3)
                  GO TO 380
  378             IF (NORM) DUM2 = DUM2*SQRT5
                  GO TO 380
  379             IF (NORM) DUM2 = DUM2*SQRT5/SQRT3
C
  380             NN = NN+1
                  DDIJ(NN) = DUM2
  400          CONTINUE
  420       CONTINUE
  520    CONTINUE
  540 CONTINUE
      RETURN
      END
C
C*MODULE EFTEI   *DECK QOUTEFP
C>    @brief  ouput QM-EFP 2EI to disk
C>
C>    @author Dan Kemp
C>            - 2010?
C>    @date   Oct 2018  Peng Xu and Tosaporn Sattasathuchana
C>            - adding MODE for the differenty types of QM/EFP TEI
C>    @details pack the 4 indecis into one
C>    @param MODE: different type of TEI involve in QM/EFP EXREP
C>                 MODE=1 ALL 4 center are in QM region.
C>                 MODE=2 two-center on QM, two-center on EFP
C>                 MODE=3 three-center on QM, One-center on EFP
C>    @see     subroutine QOUT in int2a.src
C> 
      SUBROUTINE QOUTEFP(MODE,BUFP,IX,NINTMX,GHONDO,JM)
C
      USE lrcdft, ONLY: LRFILE
      USE MX_LIMITS,ONLY:MXDFG,MXFGPT,MXSHEF,MXGEFP,mxsh,mxgtot
      USE comm_EFPBAS
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      LOGICAL IANDJ,KANDL,SAME
      LOGICAL OUT,LRINT
C
C
      DIMENSION BUFP(NINTMX),IX(*),GHONDO(*)
C
C DDK
C
C     COMMON /ELGPMT/ NELONG,NATM,NASPIN,NCT,NBNDAB,NTMLB,IPRI,LDOS
      COMMON /ERIOUT/ ISH,JSH,KSH,LSH,LSTRI,LSTRJ,LSTRK,LSTRL
C     COMMON /INT2IC/ NINTIC,ININTIC,NXXIC,LBUFPIC,LIXIC,LABSIX,NINTIX
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON /MISC  / IANDJ,KANDL,SAME
      COMMON /NLRCF / LRINT
      COMMON /NSHEL / EXX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),
     *                CF(MXGTOT),CG(MXGTOT),CH(MXGTOT),CI(MXGTOT),
     *                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),
     *                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
      COMMON /PCKLAB/ LABSIZ
      COMMON /RESTAR/ TIMLIM,IREST,NREC,INTLOC,IST,JST,KST,LST
      COMMON /SHLNOS/ QQ4,LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,
     *                MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL,
     *                NIJ,IJ,KL,IJKL
      COMMON /SHLT  / TOL,CUTOFF,ICOUNT,OUT
C
C-NEXT STATEMENT IS FOR VARIOUS IBM XLF 3.X AND 5.X COMPILERS-
C
      SAVE IJN,KLN
C
      DATA HALF /0.5D+00/
      ISSAVE=IS
      IF(LRINT)IS=LRFILE
C
C     ----- PACK THE 4 INDICES OF INTEGRAL INTO ONE WORD
C     ----- WRITE LABEL + INTEGRAL ON TAPE (IS)
C
       IF(MODE.EQ.2) THEN
         KANDL=.FALSE.
       ELSE
         KANDL = KSH .EQ. LSH
       ENDIF
C
      MINI = KMINEF(ISH,JM)
      MAXI = KMAXEF(ISH,JM)
      LOCI = KLOCEF(ISH,JM)-MINI
C
      MINJ = KMIN(JSH)
      MINL = KMIN(LSH)
      MAXJ = KMAX(JSH)
      MAXL = KMAX(LSH)
      LOCJ = KLOC(JSH)-MINJ
      LOCL = KLOC(LSH)-MINL
C
      IF(MODE.EQ.2) THEN
      MINK = KMINEF(KSH,JM)
      MAXK = KMAXEF(KSH,JM)
      LOCK = KLOCEF(KSH,JM)-MINK
      ELSE
      MINK = KMIN(KSH)
      MAXK = KMAX(KSH)
      LOCK = KLOC(KSH)-MINK
      ENDIF
C
      IJN = 0
      JMAX = MAXJ
      DO 260 I = MINI,MAXI
         I_INDEX = (I-MINI)*LSTRI + 1
C        IF (IANDJ) JMAX = I
         DO 240 J = MINJ,JMAX
            IJ_INDEX = (J-MINJ)*LSTRJ + I_INDEX
            IJN = IJN+1
            LMAX = MAXL
            KLN = 0
            DO 220 K =  MINK,MAXK
               IJK_INDEX = (K-MINK)*LSTRK + IJ_INDEX
               IF (KANDL) LMAX = K
               DO 200 L = MINL,LMAX
                  KLN = KLN+1
C                 IF(SAME  .AND.  KLN.GT.IJN) GO TO 240
                  IJKL_INDEX = (L-MINL)*LSTRL + IJK_INDEX
C
                  VAL = GHONDO( IJKL_INDEX )
                  IF(ABS(VAL).LT.CUTOFF) GO TO 200
C
                  I1 = LOCI+I
                  I2 = LOCJ+J
                  I3 = LOCK+K
                  I4 = LOCL+L
C                 IF (I1 .GE. I2) GO TO 100
C                 N = I1
C                 I1 = I2
C                 I2 = N
C                 IF (I3 .GE. I4) GO TO 180
C                 N = I3
C                 I3 = I4
C                 I4 = N
C 120             IF (I1-I3) 140,160,180
C 140             N = I1
C                 I1 = I3
C                 I3 = N
C                 N = I2
C                 I2 = I4
C                 I4 = N
C                 GO TO 180
C 160             IF (I2 .LT. I4) GO TO 140
C
                  IF (OUT) CALL INTOUT(I1,I2,I3,I4,QQ4,IJKL_INDEX,VAL)
C
                  IF ((I3 .EQ. I4).AND.(MODE.NE.2)) VAL = VAL*HALF
C
                  NPACK = ICOUNT
                  IPACK = I1
                  JPACK = I2
                  KPACK = I3
                  LPACK = I4
                  IF(LABSIZ .EQ. 2) THEN
                     LABEL = ISHFT( IPACK, 48 ) + ISHFT( JPACK, 32 ) +
     *                       ISHFT( KPACK, 16 ) + LPACK
                     IX(NPACK) = LABEL
                  ELSE IF (LABSIZ .EQ. 1) THEN
                     IF ( MOD(NPACK,2) .EQ. 0 ) THEN
                       LABEL = ISHFT( IPACK, 24 ) + ISHFT( JPACK, 16 ) +
     *                         ISHFT( KPACK,  8 ) + LPACK
                       IX( NPACK/2 ) = IX( NPACK/2 ) + LABEL
                     ELSE
                       LABEL = ISHFT( IPACK, 56 ) + ISHFT( JPACK, 48 ) +
     *                         ISHFT( KPACK, 40 ) + ISHFT( LPACK, 32 )
                       IX( (NPACK/2)+1 ) = LABEL
                     END IF
                  END IF
C
                  BUFP(ICOUNT) = VAL
                  ICOUNT = ICOUNT+1
                  IF (ICOUNT .GT. NINTMX) THEN
                     NXX = NINTMX
                     IF(MODE.EQ.2) NFT=28
                     IF(MODE.EQ.3) NFT=39
                     CALL PWRIT(NFT,BUFP,IX,NXX,NINTMX)
                     ICOUNT = 1
                     NREC = NREC+1
                  END IF
  200          CONTINUE
  220       CONTINUE
  240    CONTINUE
  260 CONTINUE
      IS=ISSAVE
      RETURN
      END
C*MODULE EFTEI   *DECK SHELLSEFP
C>    @brief  prepare shell information for QM-EFP 2EI integration
C>
C>    @author Dan Kemp
C>            - 2010?
C>    @date   Jul 2015  Peng Xu
C>            -bug fix
C>    @date   Oct 2018  Peng Xu and Tosaporn Sattasathuchana
C>            - adding MODE for the differenty types of QM/EFP TEI
C>    @param   JM: index for the type of EFP potentials
C>    @param   JNAT: counter for number of EFP atoms
C>    @param   MODE: different type of TEI involve in QM/EFP EXREP
C>                 MODE=1 ALL 4 center are in QM region.
C>                 MODE=2 two-center on QM, two-center on EFP
C>                 MODE=3 three-center on QM, One-center on EFP
C>    @see     subroutine SHELLS in int2a.src

      SUBROUTINE SHELLSEFP(MODE,NELEC,ISH,JSH,KSH,LSH,JM,JNAT)
C
      USE MX_LIMITS,ONLY:MXDFG,MXFGPT,MXSHEF,MXGEFP,mxsh,mxgsh,mxgtot,
     *     mxatm
      USE comm_EFPBAS
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     LOGICAL FLIP
      LOGICAL IANDJ,KANDL,SAME
C
      DIMENSION IX(84),IY(84),IZ(84),
     *          JX(84),JY(84),JZ(84),
     *          KX(84),KY(84),KZ(84),
     *          LX(84),LY(84),LZ(84)
C
C
      COMMON /ERIOUT/ INU,JNU,KNU,LNU,NGTI,NGTJ,NGTK,NGTL
C
      COMMON /INTDEX/ IJGT(784),IJX(784),IJY(784),IJZ(784),IK(784),
     *                KLGT(784),KLX(784),KLY(784),KLZ(784)
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /MISC  / IANDJ,KANDL,SAME
      COMMON /NSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),
     *                CF(MXGTOT),CG(MXGTOT),CH(MXGTOT),CI(MXGTOT),
     *                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),
     *                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
      COMMON /ROOT  / XX,U(13),W(13),NROOTS
      COMMON /SHLEXC/ NORGSH(3),NORGSP(3),IEXCH,NANGM,NGTH(4)
      COMMON /SHLINF/  GA(MXGSH),CSA(MXGSH),CPA(MXGSH),CDA(MXGSH),
     *                CFA(MXGSH),CGA(MXGSH),CHA(MXGSH),CIA(MXGSH),
     *                 GB(MXGSH),CSB(MXGSH),CPB(MXGSH),CDB(MXGSH),
     *                CFB(MXGSH),CGB(MXGSH),CHB(MXGSH),CIB(MXGSH),
     *                 GC(MXGSH),CSC(MXGSH),CPC(MXGSH),CDC(MXGSH),
     *                CFC(MXGSH),CGC(MXGSH),CHC(MXGSH),CIC(MXGSH),
     *                 GD(MXGSH),CSD(MXGSH),CPD(MXGSH),CDD(MXGSH),
     *                CFD(MXGSH),CGD(MXGSH),CHD(MXGSH),CID(MXGSH),
     *                AX,AY,AZ,BX,BY,BZ,RAB,CX,CY,CZ,DX,DY,DZ,RCD,
     *                NGA,NGB,NGC,NGD
      COMMON /SHLNOS/ QQ4,LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,
     +                MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL,
     +                NIJ,IJ,KL,IJKL
C
      DATA LX /   0,   1,   0,   0,   2,   0,   0,   1,   1,   0,
     *            3,   0,   0,   2,   2,   1,   0,   1,   0,   1,
     *            4,   0,   0,   3,   3,   1,   0,   1,   0,   2,
     *            2,   0,   2,   1,   1,
     *            5,   0,   0,   4,   4,   1,   0,   1,   0,   3,
     *            3,   2,   0,   2,   0,   3,   1,   1,   2,   2,
     *            1,
     *            6,   0,   0,   5,   5,   1,   0,   1,   0,   4,
     *            4,   2,   0,   2,   0,   4,   1,   1,   3,   3,
     *            0,   3,   3,   2,   1,   2,   1,   2/
      DATA KX /   0,   7,   0,   0,  14,   0,   0,   7,   7,   0,
     *           21,   0,   0,  14,  14,   7,   0,   7,   0,   7,
     *           28,   0,   0,  21,  21,   7,   0,   7,   0,  14,
     *           14,   0,  14,   7,   7,
     *           35,   0,   0,  28,  28,   7,   0,   7,   0,  21,
     *           21,  14,   0,  14,   0,  21,   7,   7,  14,  14,
     *            7,
     *           42,   0,   0,  35,  35,   7,   0,   7,   0,  28,
     *           28,  14,   0,  14,   0,  28,   7,   7,  21,  21,
     *            0,  21,  21,  14,   7,  14,   7,  14/
      DATA JX /   0,  49,   0,   0,  98,   0,   0,  49,  49,   0,
     *          147,   0,   0,  98,  98,  49,   0,  49,   0,  49,
     *          196,   0,   0, 147, 147,  49,   0,  49,   0,  98,
     *           98,   0,  98,  49,  49,
     *          245,   0,   0, 196, 196,  49,   0,  49,   0, 147,
     *          147,  98,   0,  98,   0, 147,  49,  49,  98,  98,
     *           49,
     *          294,   0,   0, 245, 245,  49,   0,  49,   0, 196,
     *          196,  98,   0,  98,   0, 196,  49,  49, 147, 147,
     *            0, 147, 147,  98,  49,  98,  49,  98/
      DATA IX /   1, 344,   1,   1, 687,   1,   1, 344, 344,   1,
     *         1030,   1,   1, 687, 687, 344,   1, 344,   1, 344,
     *         1373,   1,   1,1030,1030, 344,   1, 344,   1, 687,
     *          687,   1, 687, 344, 344,
     *         1716,   1,   1,1373,1373, 344,   1, 344,   1,1030,
     *         1030, 687,   1, 687,   1,1030, 344, 344, 687, 687,
     *          344,
     *         2059,   1,   1,1716,1716, 344,   1, 344,   1,1373,
     *         1373, 687,   1, 687,   1,1373, 344, 344,1030,1030,
     *            1,1030,1030, 687, 344, 687, 344, 687/
      DATA LY /   0,   0,   1,   0,   0,   2,   0,   1,   0,   1,
     *            0,   3,   0,   1,   0,   2,   2,   0,   1,   1,
     *            0,   4,   0,   1,   0,   3,   3,   0,   1,   2,
     *            0,   2,   1,   2,   1,
     *            0,   5,   0,   1,   0,   4,   4,   0,   1,   2,
     *            0,   3,   3,   0,   2,   1,   3,   1,   2,   1,
     *            2,
     *            0,   6,   0,   1,   0,   5,   5,   0,   1,   2,
     *            0,   4,   4,   0,   2,   1,   4,   1,   3,   0,
     *            3,   2,   1,   3,   3,   1,   2,   2/
      DATA KY /   0,   0,   7,   0,   0,  14,   0,   7,   0,   7,
     *            0,  21,   0,   7,   0,  14,  14,   0,   7,   7,
     *            0,  28,   0,   7,   0,  21,  21,   0,   7,  14,
     *            0,  14,   7,  14,   7,
     *            0,  35,   0,   7,   0,  28,  28,   0,   7,  14,
     *            0,  21,  21,   0,  14,   7,  21,   7,  14,   7,
     *           14,
     *            0,  42,   0,   7,   0,  35,  35,   0,   7,  14,
     *            0,  28,  28,   0,  14,   7,  28,   7,  21,   0,
     *           21,  14,   7,  21,  21,   7,  14,  14/
      DATA JY /   0,   0,  49,   0,   0,  98,   0,  49,   0,  49,
     *            0, 147,   0,  49,   0,  98,  98,   0,  49,  49,
     *            0, 196,   0,  49,   0, 147, 147,   0,  49,  98,
     *            0,  98,  49,  98,  49,
     *            0, 245,   0,  49,   0, 196, 196,   0,  49,  98,
     *            0, 147, 147,   0,  98,  49, 147,  49,  98,  49,
     *           98,
     *            0, 294,   0,  49,   0, 245, 245,   0,  49,  98,
     *            0, 196, 196,   0,  98,  49, 196,  49, 147,   0,
     *          147,  98,  49, 147, 147,  49,  98,  98/
      DATA IY /   1,   1, 344,   1,   1, 687,   1, 344,   1, 344,
     *            1,1030,   1, 344,   1, 687, 687,   1, 344, 344,
     *            1,1373,   1, 344,   1,1030,1030,   1, 344, 687,
     *            1, 687, 344, 687, 344,
     *            1,1716,   1, 344,   1,1373,1373,   1, 344, 687,
     *            1,1030,1030,   1, 687, 344,1030, 344, 687, 344,
     *          687,
     *            1,2059,   1, 344,   1,1716,1716,   1, 344, 687,
     *            1,1373,1373,   1, 687, 344,1373, 344,1030,   1,
     *         1030, 687, 344,1030,1030, 344, 687, 687/
      DATA LZ /   0,   0,   0,   1,   0,   0,   2,   0,   1,   1,
     *            0,   0,   3,   0,   1,   0,   1,   2,   2,   1,
     *            0,   0,   4,   0,   1,   0,   1,   3,   3,   0,
     *            2,   2,   1,   1,   2,
     *            0,   0,   5,   0,   1,   0,   1,   4,   4,   0,
     *            2,   0,   2,   3,   3,   1,   1,   3,   1,   2,
     *            2,
     *            0,   0,   6,   0,   1,   0,   1,   5,   5,   0,
     *            2,   0,   2,   4,   4,   1,   1,   4,   0,   3,
     *            3,   1,   2,   1,   2,   3,   3,   2/
      DATA KZ /   0,   0,   0,   7,   0,   0,  14,   0,   7,   7,
     *            0,   0,  21,   0,   7,   0,   7,  14,  14,   7,
     *            0,   0,  28,   0,   7,   0,   7,  21,  21,   0,
     *           14,  14,   7,   7,  14,
     *            0,   0,  35,   0,   7,   0,   7,  28,  28,   0,
     *           14,   0,  14,  21,  21,   7,   7,  21,   7,  14,
     *           14,
     *            0,   0,  42,   0,   7,   0,   7,  35,  35,   0,
     *           14,   0,  14,  28,  28,   7,   7,  28,   0,  21,
     *           21,   7,  14,   7,  14,  21,  21,  14/
      DATA JZ /   0,   0,   0,  49,   0,   0,  98,   0,  49,  49,
     *            0,   0, 147,   0,  49,   0,  49,  98,  98,  49,
     *            0,   0, 196,   0,  49,   0,  49, 147, 147,   0,
     *           98,  98,  49,  49,  98,
     *            0,   0, 245,   0,  49,   0,  49, 196, 196,   0,
     *           98,   0,  98, 147, 147,  49,  49, 147,  49,  98,
     *           98,
     *            0,   0, 294,   0,  49,   0,  49, 245, 245,   0,
     *           98,   0,  98, 196, 196,  49,  49, 196,   0, 147,
     *          147,  49,  98,  49,  98, 147, 147,  98/
      DATA IZ /   1,   1,   1, 344,   1,   1, 687,   1, 344, 344,
     *            1,   1,1030,   1, 344,   1, 344, 687, 687, 344,
     *            1,   1,1373,   1, 344,   1, 344,1030,1030,   1,
     *          687, 687, 344, 344, 687,
     *            1,   1,1716,   1, 344,   1, 344,1373,1373,   1,
     *          687,   1, 687,1030,1030, 344, 344,1030, 344, 687,
     *          687,
     *            1,   1,2059,   1, 344,   1, 344,1716,1716,   1,
     *          687,   1, 687,1373,1373, 344, 344,1373,   1,1030,
     *         1030, 344, 687, 344, 687,1030,1030, 687/
C
C     PREPARE SHELL INFORMATION/FOR HONDO INTEGRATION
C
      IF(NELEC.EQ.2) GO TO 200
C
C     ----- PERMUTE ISH AND JSH SHELLS, FOR THEIR TYPE
C     THIS IS DONE FOR SPEED REASONS.  THE CODE GETS THE RIGHT ANSWER
C     WITHOUT THE ANGULAR MOMENTUM FLIPPING, AND THEREFORE A CALLING
C     ARGUMENT ALLOWS ONE DO EXACTLY THE INTEGRAL BLOCK AS SPECIFIED,
C     SHOULD THAT BE DESIRED.
C
C
C     KTYPII = KTYPEF(ISH,JM)
C     KMINII = KMINEF(ISH,JM)
C     KMAXII = KMAXEF(ISH,JM)
C     KNGII  = KNGEF(ISH,JM)
C
C     KTYPJJ = KTYPE(JSH)
C     KMINJJ = KMIN(JSH)
C     KMAXJJ = KMAX(JSH)
C     KNGJJ  = KNG(JSH)
C
      INU = ISH
      JNU = JSH
      NGTI = NGTH(1)
      NGTJ = NGTH(2)
C
C     END IF
C
C     ----- ISHELL
C
C     IF (.NOT. TWOEFP) THEN
C     I = KATOM(INU)
C     AX = C(1,I)
C     AY = C(2,I)
C     AZ = C(3,I)
C     I1 = KSTART(INU)
C     I2 = I1+KNG(INU)-1
C     LIT = KTYPE(INU)
C     MINI = KMIN(INU)
C     MAXI = KMAX(INU)
C     LOCI = KLOC(INU)-MINI
C     NGA = 0
C     DO 140 I = I1,I2
C        NGA = NGA+1
C        GA(NGA) = EX(I)
C        CSA(NGA) = CS(I)
C        CPA(NGA) = CP(I)
C        CDA(NGA) = CD(I)
C        CFA(NGA) = CF(I)
C        CGA(NGA) = CG(I)
C        CHA(NGA) = CH(I)
C        CIA(NGA) = CI(I)
C 140 CONTINUE
C
C     ELSE
      I = KATMEF(INU,JM)
      AX = PRCORD(1,I+JNAT)
      AY = PRCORD(2,I+JNAT)
      AZ = PRCORD(3,I+JNAT)
      I1 = KSTREF(INU,JM)
      I2 = I1+KNGEF(INU,JM)-1
      LIT = KTYPEF(INU,JM)
      MINI = KMINEF(INU,JM)
      MAXI = KMAXEF(INU,JM)
      LOCI = KLOCEF(INU,JM)-MINI
      NGA = 0
      DO 141 I = I1,I2
         NGA = NGA+1                                        
         GA(NGA) = EXEF(I,JM)
         CSA(NGA) = CSEF(I,JM)
         CPA(NGA) = CPEF(I,JM)
         CDA(NGA) = CDEF(I,JM)
         CFA(NGA) = CFEF(I,JM)
         CGA(NGA) = CGEF(I,JM)
  141 CONTINUE
C     END IF
C
C     ----- JSHELL
C
      J = KATOM(JNU)
      BX = C(1,J)
      BY = C(2,J)
      BZ = C(3,J)
      J1 = KSTART(JNU)
      J2 = J1+KNG(JNU)-1
      LJT = KTYPE(JNU)
      MINJ = KMIN(JNU)
      MAXJ = KMAX(JNU)
      LOCJ = KLOC(JNU)-MINJ
      NGB = 0
      DO 160 J = J1,J2
         NGB = NGB+1
         GB(NGB) = EX(J)
         CSB(NGB) = CS(J)
         CPB(NGB) = CP(J)
         CDB(NGB) = CD(J)
         CFB(NGB) = CF(J)
         CGB(NGB) = CG(J)
         CHB(NGB) = CH(J)
         CIB(NGB) = CI(J)
  160 CONTINUE
      RAB = ((AX-BX)*(AX-BX) + (AY-BY)*(AY-BY) + (AZ-BZ)*(AZ-BZ))
C 
C     ----- PREPARE INDICES FOR PAIRS OF (I,J) FUNCTIONS
C
      IJ = 0
      JMAX = MAXJ
      DO 190 I = MINI,MAXI
         NX = IX(I)
         NY = IY(I)
         NZ = IZ(I)
C        IF (IANDJ .AND. TWOEFP .EQ. .FALSE.) JMAX = I
         DO 180 J = MINJ,JMAX
            IJ = IJ+1
            IJX(IJ) = NX+JX(J)
            IJY(IJ) = NY+JY(J)
            IJZ(IJ) = NZ+JZ(J)
            IJGT(IJ) = NGTI*(I-MINI)+NGTJ*(J-MINJ)+1
  180    CONTINUE
  190 CONTINUE
      RETURN
C     ******
C
C        K AND L SHELL
C
  200 CONTINUE
      IF(MODE.EQ.2) THEN
      KANDL=.FALSE.
      ELSE
      KANDL = KSH .EQ. LSH
      ENDIF
C
C     IF (TWOEFP .EQ. .FALSE.) THEN
C     SAME = ISH .EQ. KSH .AND. JSH .EQ. LSH
C     END IF
C
C     ----- PERMUTE KSH AND LSH SHELLS, FOR THEIR TYPE
C
C     KTYPKK = KTYPE(KSH)
C     KMINKK = KMIN(KSH)
C     KMAXKK = KMAX(KSH)
C     KNGKK  = KNG(KSH)
C
C     KTYPLL = KTYPE(LSH)
C     KMINLL = KMIN(LSH)
C     KMAXLL = KMAX(LSH)
C     KNGLL  = KNG(LSH)
C
C     IF (KTYPKK .LT. KTYPLL  .AND.  FLIP) THEN
C        KNU = LSH
C        LNU = KSH
C        NGTK = NGTH(4)
C        NGTL = NGTH(3)
      IF (MODE .EQ. 2) THEN
         KNU = KSH
         LNU = LSH
         NGTK = NGTH(3)
         NGTL = NGTH(4)
      ELSE
       IF (KTYPE(KSH) .LT. KTYPE(LSH)) THEN
          KNU = LSH
          LNU = KSH
          NGTK = NGTH(4)
          NGTL = NGTH(3)
       ELSE
          KNU = KSH
          LNU = LSH
          NGTK = NGTH(3)
          NGTL = NGTH(4)
       END IF
      END IF
C
C     ----- K SHELL
C
      IF(MODE.EQ.2) THEN
      K = KATMEF(KNU,JM)
      CX = PRCORD(1,K+JNAT)
      CY = PRCORD(2,K+JNAT)
      CZ = PRCORD(3,K+JNAT)
      K1 = KSTREF(KNU,JM)
      K2 = K1+KNGEF(KNU,JM)-1
      LKT = KTYPEF(KNU,JM)
      MINK = KMINEF(KNU,JM)
      MAXK = KMAXEF(KNU,JM)
      LOCK = KLOCEF(KNU,JM)-MINK
      NGC = 0
      DO 261 K = K1,K2
         NGC = NGC+1
         GC(NGC) = EXEF(K,JM)
         CSC(NGC) = CSEF(K,JM)
         CPC(NGC) = CPEF(K,JM)
         CDC(NGC) = CDEF(K,JM)
         CFC(NGC) = CFEF(K,JM)
         CGC(NGC) = CGEF(K,JM)
  261 CONTINUE
      ELSE
      K = KATOM(KNU)
      CX = C(1,K)
      CY = C(2,K)
      CZ = C(3,K)
      K1 = KSTART(KNU)
      K2 = K1+KNG(KNU)-1
      LKT = KTYPE(KNU)
      MINK = KMIN(KNU)
      MAXK = KMAX(KNU)
      LOCK = KLOC(KNU)-MINK
      NGC = 0
      DO 260 K = K1,K2
         NGC = NGC+1
         GC(NGC) = EX(K)
         CSC(NGC) = CS(K)
         CPC(NGC) = CP(K)
         CDC(NGC) = CD(K)
         CFC(NGC) = CF(K)
         CGC(NGC) = CG(K)
         CHC(NGC) = CH(K)
         CIC(NGC) = CI(K)
  260 CONTINUE
      ENDIF
C
C     ----- LSHELL
C
      L = KATOM(LNU)
      DX = C(1,L)
      DY = C(2,L)
      DZ = C(3,L)
      L1 = KSTART(LNU)
      L2 = L1+KNG(LNU)-1
      LLT = KTYPE(LNU)
      MINL = KMIN(LNU)
      MAXL = KMAX(LNU)
      LOCL = KLOC(LNU)-MINL
      NGD = 0
      DO 280 L = L1,L2
         NGD = NGD+1
         GD(NGD) = EX(L)
         CSD(NGD) = CS(L)
         CPD(NGD) = CP(L)
         CDD(NGD) = CD(L)
         CFD(NGD) = CF(L)
         CGD(NGD) = CG(L)
         CHD(NGD) = CH(L)
         CID(NGD) = CI(L)
  280 CONTINUE
C
C
      NROOTS = (LIT+LJT+LKT+LLT-2)/2
      RCD = ((CX-DX)*(CX-DX) + (CY-DY)*(CY-DY) + (CZ-DZ)*(CZ-DZ))
C
C     ----- PREPARE INDICES FOR PAIRS OF (K,L) FUNCTIONS
C
      KL = 0
      LMAX = MAXL
      IF(MODE.EQ.2) KANDL=.FALSE.
      DO 310 K = MINK,MAXK
         NX = KX(K)
         NY = KY(K)
         NZ = KZ(K)
         IF (KANDL) LMAX = K
         DO 300 L = MINL,LMAX
            KL = KL+1
            KLX(KL) = NX+LX(L)
            KLY(KL) = NY+LY(L)
            KLZ(KL) = NZ+LZ(L)
            KLGT(KL) = NGTK*(K-MINK)+NGTL*(L-MINL)
  300    CONTINUE
  310 CONTINUE
      MAX = KL
      DO 320 I = 1,IJ
C     IF (SAME) MAX = I
  320 IK(I) = MAX
      IJKL = IJ*KL
C     IF (SAME) IJKL = IJ*(IJ+1)/2
      RETURN
      END
C
C*MODULE INT2C   *DECK SHELLQUART_EFP
C>    @brief  compute QM-EFP 2EI
C>
C>    @author Dan Kemp
C>            - 2010?
C>    @details compute QM-EFP 2EI
C>    @date    Jul 2015   Peng XU
C>             -bug fix 
C>    @date   Oct 2018  Peng Xu and Tosaporn Sattasathuchana
C>            - adding MODE for the differenty types of QM/EFP TEI
C>
C>    @date   Jan 2020  Peng Xu and Tosaporn Sattasathuchana
C>            - add calls for rotated axis and ERIC schemes for QMEFP
C>
C>    @param   JM: index for the type of EFP potentials 
C>    @param   JNAT: counter for number of EFP atoms 
C>    @param   MODE: different type of TEI involve in QM/EFP EXREP
C>                 MODE=1 ALL 4 center are in QM region.
C>                 MODE=2 two-center on QM, two-center on EFP
C>                 MODE=3 three-center on QM, One-center on EFP
C>    @see     subroutine SHELLQUART in int2a.src

      SUBROUTINE SHELLQUART_EFP(MODE,ISH,JSH,KSH,LSH,GHONDO,JM,JNAT)
C
      USE MX_LIMITS,ONLY:MXDFG,MXFGPT,MXSHEF,MXGEFP,MXFRG,MXDPPT,mxsh,
     *     mxgtot,mxgsh,mxg2
      USE comm_EFPBAS
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      DIMENSION GHONDO(*)
      DIMENSION GROTSPD(1296)
C
      LOGICAL PACK2E,IANDJ,KANDL,SAME
C
C
      COMMON /FLIPS / IB(4,3)
C
      COMMON /INTFIL/ NINTMX,NHEX,NTUPL,PACK2E,INTTYP,IGRDTYP
      COMMON /NSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),
     *                CF(MXGTOT),CG(MXGTOT),CH(MXGTOT),CI(MXGTOT),
     *                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),
     *                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
      COMMON /OUTPUT/ NPRINT,ITOL,ICUT,NORMF,NORMP,NOPK
      COMMON /SHLEXC/ NORGSH(3),NORGSP(3),IEXCH,NANGM,NGTH(4)
C
C  FOR RYS QUADRATURE CODE
C
      LOGICAL LSHEL
      COMMON /FMCOM /XX(1)
C
C  FOR ERIC CODE, BUT ALL CODES MUST SET ERIOUT FOR OUTPUT ROUTINES
C
      LOGICAL ERICQT
      COMMON /ERIDAT/ LEN1,LEN2,LEN3,LEN4
      COMMON /ERIOUT/ INW,JNW,KNW,LNW,LSTRI,LSTRJ,LSTRK,LSTRL
C
C  FOR THE TWO ROTATED AXIS CODES
C
      INTEGER IDPOP(4,10)
      LOGICAL SP,SPD,ROTSP,ROTSPD,SPDFG
      COMMON /GOUT  / GPOPLE(768),NORGP
      COMMON /POPOUT/ LPOPI,LPOPJ,LPOPK,LPOPL
      COMMON /SHLG70/ IPL,JPL,KPL,LPL,INEW,JNEW,KNEW,LNEW
      DATA IDPOP/0,0,0,0,216,36,6,1,432,72,12,2,648,108,18,3,
     *           0,0,0,0,216,36,6,1,432,72,12,2,648,108,18,3,
     *           864,144,24,4,1080,180,30,5/
C
C     ----- SELECT THE INTEGRAL CODE FOR THIS SHELL QUARTET -----
C     THE USER INPUT SELECTION -INTTYP- HAS THE FOLLOWING MEANING:
C      INTTYP=0    BEST TIMING:
C                  USE ROTATED AXIS CODES FOR ANY S,P,D,L SHELLS,
C                  OTHERWISE PICK ERIC CODE WHENEVER POSSIBLE, BUT
C                  USE RYS CODE FOR L SHELL OR IF TOO MUCH ANG.MOM.
C      INTTYP=1    USE THE S,P,L OR S,P,D,L ROTATED AXIS CODE WHENEVER
C                  POSSIBLE, OTHERWISE RYS QUADRATURE (NO ERIC).
C      INTTYP=2    USE ERIC CODE AS MUCH AS POSSIBLE, OTHERWISE
C                  USE THE RYS QUADRATURE (NO ROTATED AXIS).
C      INTTYP=3    USE RYS POLYNOMIAL QUADRATURE FOR EVERYTHING.
C
C
      KTYPII = KTYPEF(ISH,JM)
      KMINII = KMINEF(ISH,JM)
      KMAXII = KMAXEF(ISH,JM)
      KNGII  = KNGEF(ISH,JM)
C
      KTYPJJ = KTYPE(JSH)
      KMINJJ = KMIN(JSH)
      KMAXJJ = KMAX(JSH)
      KNGJJ  = KNG(JSH)
C
      IF(MODE.EQ.2) THEN
      KTYPKK = KTYPEF(KSH,JM)
      KMINKK = KMINEF(KSH,JM)
      KMAXKK = KMAXEF(KSH,JM)
      KNGKK  = KNGEF(KSH,JM)
      ELSE
      KTYPKK = KTYPE(KSH)
      KMINKK = KMIN(KSH)
      KMAXKK = KMAX(KSH)
      KNGKK  = KNG(KSH)
      ENDIF
C
      KTYPLL = KTYPE(LSH)
      KMINLL = KMIN(LSH)
      KMAXLL = KMAX(LSH)
      KNGLL  = KNG(LSH)
C
      SP    = KTYPII.LE.2.AND.
     *        KTYPJJ.LE.2.AND.
     *        KTYPKK.LE.2.AND.
     *        KTYPLL.LE.2
      SPD   = KTYPII.LE.3.AND.
     *        KTYPJJ.LE.3.AND.
     *        KTYPKK.LE.3.AND.
     *        KTYPLL.LE.3
      IF(SP) SPD=.FALSE.
      LSHEL = (KMAXII-KMINII+1).EQ.4.OR.
     *        (KMAXJJ-KMINJJ+1).EQ.4.OR.
     *        (KMAXKK-KMINKK+1).EQ.4.OR.
     *        (KMAXLL-KMINLL+1).EQ.4
      LQSUM = KTYPII + KTYPJJ + KTYPKK + KTYPLL - 4
C
      ROTSP  = SP
      ROTSPD = SPD
      SPDFG = KTYPII.LE.5.AND.
     *        KTYPJJ.LE.5.AND.
     *        KTYPKK.LE.5.AND.
     *        KTYPLL.LE.5
      ERICQT = .NOT.LSHEL .AND. LQSUM.LE.5 .AND. SPDFG
C
C        RYS QUADRATURE IS A BIT FASTER AT UNCONTRACTED QUARTETS
C
      KQCON = KNGII * KNGJJ * KNGKK * KNGLL + INTTYP
      IF(KQCON.EQ.1  .AND.  .NOT.SP) THEN
         ROTSP  = .FALSE.
         ROTSPD = .FALSE.
         ERICQT = .FALSE.
      END IF
C
C        INPUT OVERRIDES
C
      IF(INTTYP.EQ.1) ERICQT = .FALSE.
      IF(INTTYP.GE.2) ROTSP  = .FALSE.
      IF(INTTYP.GE.2) ROTSPD = .FALSE.
      IF(INTTYP.EQ.3) ERICQT = .FALSE.
C-- Keep this for testing individual scheme 
c      ROTSP  = .FALSE.
c      ROTSPD = .FALSE.
c      ERICQT = .FALSE.
C
C        THE VARIOUS PACKAGES ARE CALLED BELOW IN THE ORDER OF
C        FIRST ROTATED AXIS, THEN ERIC, FINALLY RYS QUADRATURE.
C        SINCE EACH CODE RETURNS, A QUARTET IS NEVER DONE TWICE.
C
C  ROTATED AXIS CODE FOR PURE SP SHELL QUARTET
C
      IF (ROTSP) THEN
        IPL = ISH
        JPL = JSH
        KPL = KSH
        LPL = LSH
        INW = ISH
        JNW = JSH
        KNW = KSH
        LNW = LSH
        NORGP = NORGSP(IEXCH)
        NORGH = NORGSH(IEXCH)
C
        CALL GENR70_QMEFP(JM,JNAT,MODE)
C
C  SAVE TO OUTPUT ARRAY WITH HONDO INDEXING
C
        MINI = KMINII
        MAXI = KMAXII
        MINJ = KMINJJ
        MAXJ = KMAXJJ
        MINK = KMINKK
        MAXK = KMAXKK
        MINL = KMINLL
        MAXL = KMAXLL
C
        II = 1
        DO I = MINI, MAXI
          IP = (I-1)*LPOPI + 1
          IJ  = II
          DO J = MINJ, MAXJ
            IJP = (J-1)*LPOPJ + IP
            IJK  = IJ
            DO K = MINK, MAXK
              IJKP = (K-1)*LPOPK + IJP
              IJKL  = IJK
              DO L = MINL, MAXL
                IJKLP = (L-1)*LPOPL + IJKP
                GHONDO(IJKL+NORGH) = GPOPLE(IJKLP+NORGP)
                IJKL = IJKL  + LEN1
              END DO
              IJK  = IJK  + LEN2
            END DO
            IJ  = IJ  + LEN3
          END DO
          II = II + LEN4
        END DO
        LSTRI = LEN4
        LSTRJ = LEN3
        LSTRK = LEN2
        LSTRL = LEN1
        RETURN
C       ******
C
C  ROTATED AXIS CODE FOR QUARTET CONTAINING AT LEAST ONE D FUNCTION
C
      ELSE IF (ROTSPD) THEN
        IPL = ISH
        JPL = JSH
        KPL = KSH
        LPL = LSH
        INW = ISH
        JNW = JSH
        KNW = KSH
        LNW = LSH
        LSTRI = NGTH(1)                
        LSTRJ = NGTH(2)
        LSTRK = NGTH(3)
        LSTRL = NGTH(4)

        IF((MODE.EQ.3).AND.(KTYPKK.LT.KTYPLL)) THEN
          KPL = LSH
          LPL = KSH
          KNW = LSH
          LNW = KSH

          LSTRK = NGTH(4)
          LSTRL = NGTH(3)
        ENDIF
C
        CALL GENR03_QMEFP(GROTSPD,JM,JNAT,MODE)
C
C  SAVE TO OUTPUT ARRAY WITH HONDO INDEXING
C
        NORGH = NORGSH(IEXCH)
c        IANDJ = ISH.EQ.JSH
        KANDL = .FALSE. 
c        SAME  = ISH.EQ.KSH  .AND.  JSH.EQ.LSH
        IF(MODE.EQ.3) KANDL = KSH.EQ.LSH 
        IF(NOPK.EQ.0) SAME=.FALSE.
C
        IEX=1
        IBB = IB(1,IEX)
        JBB = IB(2,IEX)
        KBB = IB(3,IEX)
        LBB = IB(4,IEX)
C
        MINI = KMINII
        MAXI = KMAXII
        MINJ = KMINJJ
        MAXJ = KMAXJJ
        MINK = KMINKK
        MAXK = KMAXKK
        MINL = KMINLL
        MAXL = KMAXLL
        IF((MODE.EQ.3).AND.(KTYPKK.LT.KTYPLL)) THEN                                   
        MINK = KMINLL                                                                 
        MAXK = KMAXLL                                                                 
        MINL = KMINKK                                                                 
        MAXL = KMAXKK                                                                 
        ENDIF                                                                         
C
        IJN = 0
        JMAX = MAXJ
        DO I = MINI, MAXI
          IHONDO = (I-MINI)*LSTRI + 1
          IROTAX = IDPOP(IBB,I)  + 1
          DO 340 J = MINJ, JMAX
            IJHONDO = (J-MINJ)*LSTRJ + IHONDO
            IJROTAX = IDPOP(JBB,J)  + IROTAX
            IJN = IJN+1
            LMAX=MAXL
            KLN=0
            DO K = MINK, MAXK
              IJKHONDO = (K-MINK)*LSTRK + IJHONDO
              IJKROTAX = IDPOP(KBB,K)  + IJROTAX
              IF(KANDL) LMAX=K
              DO L = MINL, LMAX
                KLN = KLN+1
                IJKLHONDO = (L-MINL)*LSTRL + IJKHONDO
                IJKLROTAX = IDPOP(LBB,L)  + IJKROTAX
                GHONDO(IJKLHONDO+NORGH) = GROTSPD(IJKLROTAX)
              END DO
            END DO
  340     CONTINUE
        END DO
c        LSTRI = LEN4
c        LSTRJ = LEN3
c        LSTRK = LEN2
c        LSTRL = LEN1
        RETURN
C       ******
C
C  USE ERIC FAST CODES, REQUIRES THAT LQSYM.LE.5 AND NO L-SHELLS
C  NOTE THAT THERE IS CODE COPYING ERIC BUFFERS INTO THE HONDO
C  FORMAT BUFFER AT THE END OF THE -ERIC- ROUTINE.
C
      ELSE IF (ERICQT) THEN
        NORGH = NORGSH(IEXCH)
        CALL ERIC_QMEFP(ISH,JSH,KSH,LSH,GHONDO(1+NORGH),JM,JNAT,MODE)
        RETURN
C       ******
C
C  GENERAL CASE = HONDO/RYS QUADRATURE: ANY S,P,D,F,G, OR L SHELLS
C
      ELSE
        CALL VALFM(LOADFM)
        IDDIJ = LOADFM + 1
        NEED  = 49*MXG2
        CALL GETFM(NEED)
        CALL SHELLSEFP(MODE,1,ISH,JSH,KSH,LSH,JM,JNAT)
        CALL SHELLSEFP(MODE,2,ISH,JSH,KSH,LSH,JM,JNAT)
        CALL IJPRIMEFP(XX(IDDIJ))
        NORGH = NORGSH(IEXCH)
        IF(NOPK.NE.0) THEN 
        CALL ZQOUTEFP(MODE,GHONDO)
        END IF
        IF(LQSUM.EQ.0) THEN
           CALL S0000(GHONDO(1+NORGH),XX(IDDIJ))
        ELSE
           CALL GENRAL(GHONDO(1+NORGH),XX(IDDIJ))
        END IF
        CALL RETFM(NEED)
        RETURN
C       ******
C
      END IF
      END
C
C*MODULE EFTEI   *DECK ZQOUTEFP
C>    @brief  zeroing the integrals
C>
C>    @author Dan Kemp
C>            - 2010?
C>    @date   Oct 2018  Peng Xu and Tosaporn Sattasathuchana
C>            - adding MODE for the differenty types of QM/EFP TEI
C>    @details zeroing QM-EFP 2EI
C>    @param MODE: different type of TEI involve in QM/EFP EXREP
C>                 MODE=1 ALL 4 center are in QM region.
C>                 MODE=2 two-center on QM, two-center on EFP
C>                 MODE=3 three-center on QM, One-center on EFP
C>    @see    subroutine ZQOUT in int2a.src
C>
      SUBROUTINE ZQOUTEFP(MODE,GHONDO)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      DIMENSION GHONDO(*)
C
      LOGICAL IANDJ,KANDL,SAME
C
      COMMON /INTDEX/ IJGT(784),IJX(784),IJY(784),IJZ(784),IK(784),
     *                KLGT(784),KLX(784),KLY(784),KLZ(784)
      COMMON /MISC  / IANDJ,KANDL,SAME
      COMMON /SHLNOS/ QQ4,LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,
     *                MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL,
     *                NIJX,IJ,KL,IJKL
C
      PARAMETER (ZERO=0.0D+00)
C
C     ----- ZERO HONDO CONVENTIONAL INTEGRAL OUTPUT REGION -----
C
      IJN = 0
      JMAX = MAXJ
      IF(MODE.EQ.2) KANDL=.FALSE.
      DO 260 I = MINI,MAXI
C        IF (IANDJ) JMAX = I
         DO 240 J = MINJ,JMAX
            IJN = IJN+1
            N1 = IJGT(IJN)
            LMAX = MAXL
            KLN = 0
            DO 220 K =  MINK,MAXK
              IF (KANDL) LMAX = K
               DO 200 L = MINL,LMAX
                  KLN = KLN+1
                  NN = N1+KLGT(KLN)
                  GHONDO(NN) = ZERO
  200          CONTINUE
  220       CONTINUE
  240    CONTINUE
  260 CONTINUE
      RETURN
      END
C*MODULE EFTEI   *DECK EXCHNG_EF
C>    @brief Schwarz inequality screening for QM/EFP EXREP
C>    
C>    @author   Peng Xu and Tosaporn Sattasathuchana
C>            - 2018
C>    @details  Calculate largest exchange integrals of each shell
C>              between QM and EFP
C>    @date     Sep 2019   Peng Xu and Tosaporn Sattasathuchana
C>            - adding static parallel 
C>            - adding OpenMP parallel
C>    @date     Jan 2020  Peng Xu and Tosaporn Sattasathuchana
C>            - add rotated axis scheme 
C>    @param    NSHELT : total number of EFP shells
C>    @see    subroutine EXCHNG in int2a.src
C>
      SUBROUTINE EXCHNG_EF(XINTS,GHONDO,DDIJ,MAXG,NSHELT)
C
C$    USE params, ONLY: intomp
C$    USE efteiomp
      USE MX_LIMITS,ONLY:MXSHEF,MXGEFP,MXDFG,MXFGPT,MXFRG,MXDPPT,mxsh,
     *     mxgtot
      use comm_FRGTYP
      USE comm_EFPBAS
      USE comm_FRGINF
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      DIMENSION XINTS(NSHELT,NSHELL),GHONDO(MAXG),DDIJ(*)
C
      LOGICAL POPLE,OUT,SOME,IANDJ,GOPARR,DSKWRK,MASWRK,SAVEINT,NXT,
     *        DOESP,LCUT

      LOGICAL SLB_EFAI
      PARAMETER (HUND=100.0D+00)
      PARAMETER (ZERO=0.0D+00)
      PARAMETER (ONE=1.0D+00)
C
      COMMON /INTAC2/ EI1,EI2,CUX
      COMMON /FLIPS / IB(4,3)
      COMMON /FMOINF/ NFG,NLAYER,NATFMO,NBDFG,NAOTYP,NBODY,NSEGM
      COMMON /FMORUN/ ESPSCF,E0SCF(2),EMP2S,IDAFMO,ICURFG,JCURFG,KCURFG,
     *                ICURLAY,ICURUNT,NAT1E,NCURSH,NGAU,ICURPOP,IFMOSTP,
     *                MONCOR,NEEDR,MODRST,NORBPROJ,NUNESP,ISKIPESP,
     *                IESDPPC,IDOPROP,MP2RUN,ICURIT,IDMFMO,IDDFMO,
     *                IDDCUR,NDDLEFT,IVMFMO,nzmtfmo,ifmobas,itmfmo(2)
      COMMON /GOUT  / GPOPLE(768),NORGP
      COMMON /INTDEX/ IJGT(784),IJX(784),IJY(784),IJZ(784),IK(784),
     *                KLGT(784),KLX(784),KLY(784),KLZ(784)
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON /ELGIDX/ LCUT
      COMMON /NSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),
     *                CF(MXGTOT),CG(MXGTOT),CH(MXGTOT),CI(MXGTOT),
     *                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),
     *                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
      COMMON /OUTPUT/ NPRINT,ITOL,ICUT,NORMF,NORMP,NOPK
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /SHLEXC/ NORGSH(3),NORGSP(3),IEXCH,NANGM,NGTH(4)
      COMMON /SHLG70/ ISH,JSH,KSH,LSH,IJKLXX(4)
      COMMON /SHLNOS/ QQ4,LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,
     *                MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL,
     *                NIJ,IJ,KL,IJKL
      COMMON /SHLT  / TOL,CUTOFF,ICOUNT,OUT
      COMMON /INTFIL/ NINTMX,NHEX,NTUPL,PACK2E,INTTYP,IGRDTYP
      COMMON /POPOUT/ LPOPI,LPOPJ,LPOPK,LPOPL
C

C     ----- COMPUTE ALL EXCHANGE INTEGRALS [II,JJ/II,JJ] -----
C     WHERE II is EFP and JJ is QM
C
      MODE = 2
      NSH2= NSHELT*NSHELL
      NXT=.FALSE.
      IPCOUNT = ME - 1
      IF(NXT) IPCOUNT = - 1
      NEXT = -1
      DOESP=NFG.NE.0.AND.NCURSH.NE.0
      IF (GOPARR.OR.DOESP) CALL VCLR(XINTS,1,NSH2)
      SLB_EFAI= GOPARR.AND.NFRG.GT.10
C
      SOME = NPRINT.NE.-5 .AND. MASWRK
      SAVEINT=NFG.EQ.0.OR.NCURSH.EQ.0.OR.IFMOSTP.EQ.6

C     CALL THREADED VERSION IF POSSIBLE
C
C       The code below is only active
C       with -fopenmp/-qopenmp compiler option
C
C$    IF (intomp.ne.0) THEN
C$      CALL efp2ei_omp_exchng(xints,maxg,nshelt,inttyp)
C$      RETURN
C$    END IF
C
C
C      IF(SOME) THEN
C         TIM = ZERO
C         CALL TSECND(TIM)
C         TIM0 = TIM
C      ENDIF
C
C      CALL BASCHK(LMAX)
      LMAX=0
      DO MJ = 1,NFRG
      JM = ISET(MJ)
      CALL BASEFP(LMAX2,JM)
      IF (LMAX2 .GT. LMAX) LMAX=LMAX2
      END DO
                    NANGM =  4
      IF(LMAX.EQ.2) NANGM =  6
      IF(LMAX.EQ.3) NANGM = 10
      IF(LMAX.EQ.4) NANGM = 15
      IF(LMAX.EQ.5) NANGM = 21
      IF(LMAX.EQ.6) NANGM = 28
      NGTH(4) = 1
      NGTH(3) = NGTH(4) * NANGM
      NGTH(2) = NGTH(3) * NANGM
      NGTH(1) = NGTH(2) * NANGM
      IF(NOPK.EQ.0) THEN
         NORGSH(1) = 0
         NORGSH(2) = NORGSH(1) + NANGM**4
         NORGSH(3) = NORGSH(2) + NANGM**4
         NORGSP(1) = 0
         NORGSP(2) = 256
         NORGSP(3) = 512
      ELSE
         DO I=1,3
            NORGSH(I) = 0
            NORGSP(I) = 0
         ENDDO
      END IF
C
C        THE IDEA IS TO DO EVEN SMALL INTEGRALS, BELOW THE USUAL
C        CUTOFF THRESHHOLDS, BY RESETTING TOLERANCES TIGHTLY.
C
      TOLSV = TOL
      TOL = 75.0D+00
C
      EI1SV = EI1
      EI2SV = EI2
      CUXSV = CUX
      EI1 = 1.0D-17
      EI2 = 1.0D-17
      CUX = 50.0D+00
C
      IEXCH = 1
      NORGP = 0
      QQ4   = ONE
      NINT  = 0
C
C     ----- LOOP OVER ALL SHELL BLOCKS -----
C
C      IJIJ = 0
C     ----- I SHELL -----
      JAO  = 0
      JNAT = 0
      ISHLEF = 0

      DO 1000 MJ = 1,NFRG
      JM = ISET(MJ)
C
C           ----- GO PARALLEL! -----
C
      IF (SLB_EFAI) THEN
         IPCOUNT = IPCOUNT + 1
         IF (MOD(IPCOUNT,NPROC).NE.0) GO TO 700
      END IF
      DO 600 ISH = 1,NSHLEF(JM)
C ----- TURNING OFF NOW FOR TESTING ---------------------
C         IF (GOPARR.AND.NXT) THEN
C            IPCOUNT = IPCOUNT + 1
C            IF (IPCOUNT.GT.NEXT) CALL DDI_DLBNEXT(NEXT)
C            IF (NEXT.NE.IPCOUNT) THEN
C               IJIJ = IJIJ+ISH
C               GO TO 600
C            ENDIF
C         END IF
C -------------------------------------------------------
         DO 500 JSH = 1,NSHELL
C            IJIJ = IJIJ+1
C
C           SKIP UNNEEDED OFF-DIAGONAL BLOCKS FOR FMO ESP SCREENING.
C
            IF(DOESP.AND.ISH.GT.NCURSH.AND.JSH.LE.NCURSH) GO TO 500
C
C           ----- GO PARALLEL! -----
C
C            IF (GOPARR.AND..NOT.NXT) THEN
            IF (GOPARR.AND..NOT.SLB_EFAI) THEN
               IPCOUNT = IPCOUNT + 1
               IF (MOD(IPCOUNT,NPROC).NE.0) GO TO 500
            END IF
C
C     USE POPLE CODE FOR ANY PURE SP INTEGRAL BLOCKS,
C     USE HONDO RYS POLYNOMIAL CODE FOR OTHER BLOCKS
C
         KTYPII = KTYPEF(ISH,JM)
         KTYPJJ = KTYPE(JSH)
         LQSUM  = KTYPII + KTYPJJ + KTYPII + KTYPJJ - 4
c
         POPLE=.TRUE.
         IF(INTTYP.GE.2)     POPLE=.FALSE.
         IF(KTYPII.GT.2)     POPLE=.FALSE.
         IF(KTYPJJ.GT.2)     POPLE=.FALSE.

         IF(POPLE) THEN
               KSH=ISH
               LSH=JSH
               CALL GENR70_QMEFP(JM,JNAT,MODE)
            ELSE
               CALL SHELLSEFP(MODE,1,ISH,JSH,ISH,JSH,JM,JNAT)
               CALL IJPRIMEFP(DDIJ)
               CALL SHELLSEFP(MODE,2,ISH,JSH,ISH,JSH,JM,JNAT)
               CALL ZQOUTEFP(MODE,GHONDO)

            IF(LQSUM.EQ.0) THEN
               CALL S0000(GHONDO,DDIJ)
            ELSE 
               CALL GENRAL(GHONDO,DDIJ)
            END IF
         END IF
C
C     ----- PICK OUT LARGEST EXCHANGE INTEGRAL FOR THIS BLOCK -----
C
         VMAX = ZERO
         MINI = KMINEF(ISH,JM)
         MAXI = KMAXEF(ISH,JM)

         MINJ = KMIN(JSH)
         JMAX = KMAX(JSH)
         IANDJ=.FALSE.
         IJN = 0

         DO 300 I=MINI,MAXI
            IF(IANDJ) JMAX = I
            DO 200 J=MINJ,JMAX
               IF(POPLE) THEN
                  NN = (I-1)*LPOPI + (J-1)*LPOPJ
     *                +(I-1)*LPOPK + (J-1)*LPOPL + 1
                  VAL = GPOPLE(NN)
               ELSE
                  IJN = IJN+1
                  NN = IJGT(IJN) + KLGT(IJN)
                  VAL = GHONDO(NN)
               END IF
               IF(VAL.GT.ZERO) NINT=NINT+1
               IF(VAL.GT.VMAX) VMAX=VAL
  200       CONTINUE
  300    CONTINUE
C         XINTS(ISH+ISHEF,JSH)=SQRT(VMAX)
         XINTS(ISH+ISHLEF,JSH)=SQRT(VMAX)
  500    CONTINUE
  600 CONTINUE
  700 CONTINUE
      JAO = JAO + NUMEF(JM)
      JNAT = JNAT + NATEF(JM)
      ISHLEF = ISHLEF + NSHLEF(JM)
 1000 CONTINUE
C
C     ----- SUM UP PARTIAL CONTRIBUTIONS IF PARALLEL -----
C
      IF (GOPARR) THEN
         CALL DDI_GSUMF(1050,XINTS,NSH2)
         CALL DDI_GSUMI(1051,NINT ,1)
         !IF(NXT) CALL DDI_DLBRESET
      END IF
C
C      IF(OUT) THEN
C         WRITE(IW,*) 'MAX EXCHANGE INTEGRAL IN SHELL'
C         CALL PRTRI(XINTS,NSHELL)
C      END IF
C
C      IF(SOME) THEN
C         CALL TSECND(TIM)
C         TEXCH = TIM-TIM0
C         WRITE(IW,9000) NINT,TEXCH
C      ENDIF
C
      TOL = TOLSV
      EI1 = EI1SV
      EI2 = EI2SV
      CUX = CUXSV
C
C     DURING FMO ESP RUNS, EXCHANGE INTEGRALS HAVE DIFFERENT SIZE SO
C     ONE CANNOT WRITE THEM TO THE SAME RECORD. THE ONLY EXCEPTION IS
C     THE SEPARATED DIMER ENERGIES WHERE THERE IS JUST ONE SET OF 2E
C     INTEGRALS.
C     ELONGATION METHOD ALSO MUST DECIDE ON THIS.
C
C      IF(SAVEINT.AND.(.NOT.LCUT)) CALL DAWRIT(IDAF,IODA,XINTS,NSH2,54,0)
      RETURN
 9000 FORMAT(1X,'SCHWARZ INEQUALITY OVERHEAD:',I10,' INTEGRALS, T=',
     *       F12.2)
      END
C*MODULE EFTEI  *DECK EFP2EI_MODE1
C>    @brief    form fock matrix for QM/EFP EXCHANG REPULSION from TEI 
C>              MODE 1
C>    @author   Peng Xu and Tosaporn Sattasathuchana
C>            - 2018
C>    @details  Calculate fock matrix for QM/EFP EXCHANG REPULSION 
C>              on the fly where all 4 center are in QM region
C>    @date     Sep 2019   Peng Xu and Tosaporn Sattasathuchana
C>            - adding OpenMP parallel
C>    @param VEC: MO coefficient of QM
C>    @param SAB: OVERLAP B/W QM and EFP in MO basis
C>    @param SMAT: OVERLAP B/W QM and EFP in AO basis
C>    @param TMAT: Kinetic energy integrals in AO basis
C>    @param COREA: number of QM core orbitals
C>    @param MXBF: maximum basis functions of EFP
C>    @param PROVEC: EFP LMO coefficients
C>    @param MXMO: maximum MO for EFP fragment
C>    @param NFRG: number of EFP fragment
C>    @param VMO: FOCK contribution of QM/EFP EXREP in MO basis
C>    @param SKN: square of SAB
C>    @param SPRIME: is SKN backtransformed in AO basis
C>    @param S:  triagular version of SPRIME
C>    @param VAO: FOCK contribution of QM/EFP EXREP in AO basis
C>    @param VAOSQ: is VAO in square matrix
C>    @see SUBROUTINE TWOEI in int2a.src

      SUBROUTINE EFP2EI_MODE1(TYPSCF,DIRSCF,DIRNLO,DIRTRF,DIRCIS,
     *                 INTTYP,SCHWRZ,NINT,NSCHWZ,L1,L2,BUFP,BUFK,IX,
     *                 NINTMX,XINTS,NSH2,GHONDO,MAXG,DDIJ,DSH,DNLO,FNLO,
     *                 NFLMAT,VEC,SAB,SMAT,TMAT,COREA,MXBF,PROVEC,WRK,
     *                 MXMO,NFRG,SPRIME,SKN,S,VAO,VAOSQ,WRK2,VMO,VMI)
C
C       MANY ARGUMENTS ARE OPTIONAL, YOU MUST ALLOCATE STORAGE FOR
C              ALL CALLS: GHONDO, DDIJ, XINTS
C           CONVENTIONAL: BUFP, IX, AND POSSIBLY BUFK
C             DIRECT SCF: IA, DSH, DA, FA, AND POSSIBLY DB, FB
C             DIRECT CIS: IA, DSH, DA, FA, AND POSSIBLY DB, FB
C             DIRECT NLO: DNLO, FNLO, NFLMAT
C  DIRECT TRANSFORMATION: BUFP, IX
C     RESPONSE EQUATIONS: MUST DEFINE NFLMAT.NE.1
C  NOTE THAT THE TYPE OF FOCK MATRIX BUILT (TYPSCF) DOES NOT
C  NECESSARILY HAVE TO MATCH THIS RUNS SCFTYP VALUE.
C
C$    USE params, ONLY: intomp
C$    USE efteiomp
C
      USE MX_LIMITS,ONLY:MXDFG,MXFRG,MXFGPT,MXSHEF,MXGEFP,mxatm,mxsh,
     *     mxgtot,mxao
      USE comm_PAULMO
      use comm_FRGTYP
      USE comm_EFPBAS
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      LOGICAL DIRSCF,DIRNLO,DIRTRF,DIRCIS,CMBDIR
      LOGICAL OUT,SCHWRZ,SCHSKP,GOPARR,DSKWRK,MASWRK,DLB,SLB,C1GRP
      LOGICAL SKIPA,SKIPB,SKIPC,NPSYM
      LOGICAL PK,PANDK,NOTPK,BLOCK,GPSAVE,SCREEN
      LOGICAL LTRMST
      LOGICAL DLB_EFAI
C
      CHARACTER*8 INAOFL
C
      DIMENSION BUFP(NINTMX),BUFK(NINTMX),IX(*),XINTS(NSH2),
     *          GHONDO(MAXG),FA(L2),DB(L2),FB(L2),
     *          DSH(NSH2),DDIJ(*),DNLO(L1,L1),FNLO(L1,L1)
      DIMENSION MI(48),MJ(48),MK(48),M0(48)
C
      INTEGER COREA
      DIMENSION VXR(L1,L1),
     *          SAB(NA,MXMO),VEC(L1,L1),WRK(MXBF),
     *          PROVEC(MXBF,NTMO),WRK2(L1)
      DIMENSION SMAT(L1,MXBF),TMAT(L1,MXBF)
      DIMENSION SPRIME(L1,L1), SKN(NA,NA), VAO(L2),VMO(L1,L1)
      DIMENSION VAOSQ(L1,L1), S(L2), VMI(L1,L1)
C
      PARAMETER (FOUR=4.0D+00)
      PARAMETER (HALF=0.5D+00)
      PARAMETER (TWO=2.0D+00)
      PARAMETER (ZERO=0.0D+00)
C
      COMMON /ELGFIL/ NFILE,INAOFL
      COMMON /ELGTRM/ LTRMST,NFLTRM,NRCTRM,NPSTRM,NHTSHL
      COMMON /INT2IC/ NINTIC,ININTIC,NXXIC,LBUFPIC,LIXIC,LABSIX,NINTIX
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON /NSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),
     *                CF(MXGTOT),CG(MXGTOT),CH(MXGTOT),CI(MXGTOT),
     *                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),
     *                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
      COMMON /OUTPUT/ NPRINT,ITOL,ICUT,NORMF,NORMP,NOPK
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /PKFIL / PK,PANDK,BLOCK
      COMMON /RESTAR/ TIMLIM,IREST,NREC,INTLOC,IST,JST,KST,LST
      COMMON /SCINP / VLAMB,SCREEN
      COMMON /SHLEXC/ NORGSH(3),NORGSP(3),IEXCH,NANGM,NGTH(4)
      COMMON /SHLG70/ ISH,JSH,KSH,LSH,IJKLXX(4)
      COMMON /SHLNOS/ QQ4,LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,
     *                MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL,
     *                NIJ,IJ,KL,IJKL
      COMMON /SHLT  / TOL,CUTOFF,ICOUNT,OUT
      COMMON /SYMTRY/ MAPSHL(MXSH,48),MAPCTR(MXATM,48),
     *                T(432),INVT(48),NT
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /IJPAIR/ IA(MXAO)
C
C          ----- TWO-ELECTRON INTEGRALS -----
C     THIS VERSION CAN HANDLE S,P,D,F,G,H,I AND L SHELLS
C
      TIM = ZERO
      CALL TSECND(TIM)
C
C
      CMBDIR= DIRSCF .OR. DIRNLO .OR. DIRTRF .OR. DIRCIS
C
C        THE OLD FASHIONED PARALLEL INTEGRAL TRANSFORMATIONS DO NOT
C        ALLOW THE AO INTEGRAL WORK TO BE RUN IN PARALLEL.  THE MODERN
C        DISTRIBUTED MEMORY TRANSFORMATIONS DO NOT PASS THROUGH HERE.
C
      ICONT=  0
      GPSAVE = GOPARR
      IF(DIRTRF) GOPARR=.FALSE.
      MODE=1
      IDUMMY=0
      DUMMY= 0.D0
      ICNTSV = ICOUNT
      NRECSV = NREC
      INTLSV = INTLOC
      ICOUNT = 1
      NREC   = 0
      INTLOC = 1
C
C     ----- INITIALIZATION FOR PARALLEL WORK -----
C     BOTH STATIC AND DYNAMIC LOAD BALANCING ARE IMPLEMENTED BELOW
C
      SLB = GOPARR  .AND.  IBTYP.EQ.0
      DLB = GOPARR  .AND.  IBTYP.EQ.1
      DLB_EFAI= DLB .AND. NFRG.GT.10 .AND. L1.LE.100
C$    IF (DIRSCF.and.(intomp.ne.0)) DLB_EFAI=.false.
      NEXT = -1
      MINE = -1
      IPCOUNT = ME - 1
C
      C1GRP = NT.EQ.1
C
      CALL BASCHK(LMAX)
                    NANGM =  4
      IF(LMAX.EQ.2) NANGM =  6
      IF(LMAX.EQ.3) NANGM = 10
      IF(LMAX.EQ.4) NANGM = 15
      IF(LMAX.EQ.5) NANGM = 21
      IF(LMAX.EQ.6) NANGM = 28
      NGTH(4) = 1
      NGTH(3) = NGTH(4) * NANGM
      NGTH(2) = NGTH(3) * NANGM
      NGTH(1) = NGTH(2) * NANGM
      IF(NOPK.EQ.0) THEN
         NORGSH(1) = 0
         NORGSH(2) = NORGSH(1) + NANGM**4
         NORGSH(3) = NORGSH(2) + NANGM**4
         NORGSP(1) = 0
         NORGSP(2) = 256
         NORGSP(3) = 512
      ELSE
         DO I=1,3
            NORGSH(I) = 0
            NORGSP(I) = 0
         ENDDO
      END IF
C COMPUTE QM/EFP OVERLAP MATRIX
      CALL VCLR(VAO,1,L2)
      CALL VCLR(VAOSQ,1,L1*L1)
      JAO=0
      JMO=0
      JNAT=0
      DO 1000 NJ=1,NFRG
      JM = ISET(NJ)

C
C     ----- GO PARALLEL! -----
C
      IF (DLB_EFAI) THEN
         MINE = MINE + 1
         IF (MINE.GT.NEXT) CALL DDI_DLBNEXT(NEXT)
         IF (NEXT.NE.MINE) GO TO 950
      END IF
      CALL VCLR(SMAT,1,L1*MXBF)
      CALL VCLR(SAB,1,NA*MXMO)
      CALL VCLR(SPRIME,1,L1*L1)
      CALL VCLR(SKN,1,NA*NA)

      CALL STINT(SMAT,TMAT,MXBF,L1,NAT,NATEF(JM),NGSSEF(JM),NSHLEF(JM),
     *           EXEF(1,JM),CSEF(1,JM),CPEF(1,JM),CDEF(1,JM),CFEF(1,JM),
     *           CGEF(1,JM),KSTREF(1,JM),KATMEF(1,JM),KTYPEF(1,JM),
     *           KNGEF(1,JM),KLOCEF(1,JM),KMINEF(1,JM),KMAXEF(1,JM),
     *           L1,NSHELL,EX,CS,CP,CD,CF,CG,KSTART,KATOM,KTYPE,KNG,
     *           KLOC,KMIN,KMAX,PRCORD(1,JNAT+1),C)

C          SMAT ADDRESSING CHANGED FEB 2010
      CALL TFSQP(SAB,SMAT,VEC,PROVEC(1,JMO+1),WRK,
     *           L1,NUMEF(JM),NA,NORBEF(NJ),NA,NORBEF(NJ),.FALSE.)

      DO K=COREA,NA
       DO N=COREA,NA
        DO J=1,NORBEF(NJ)
        SKN(K,N)=SKN(K,N)+(SAB(K,J)*SAB(N,J))
        END DO
       END DO
      END DO
      DO NALPHA=1,L1
       DO LAMBDA=1,L1
        TEMP2=ZERO
        DO N=COREA,NA
         TEMP1=ZERO
         DO K=COREA,NA
         TEMP1=TEMP1+SKN(K,N)*VEC(LAMBDA,K)
         END DO
        TEMP2=TEMP2+TEMP1*VEC(NALPHA,N)
        END DO
       SPRIME(LAMBDA,NALPHA)=SPRIME(LAMBDA,NALPHA)+TEMP2
       END DO
      END DO
C --- P. XU: S is triangular version of SPRIME ---
      DO I = 1,L1
       DO J = 1,I
          IJ = IA(I) + J
          S(IJ) = SPRIME(I,J)
       END DO
      END DO
      CALL VCLR(DSH,1,NSH2)
      IF(SCHWRZ) THEN
        DUMMY=0.0D+00
        CALL SHLDEN(TYPSCF,S,S,DUMMY,DSH,IA,L1,L2,NSH2,1)
      ENDIF

C
C     CALL THREADED VERSION IF POSSIBLE
C
C       The code below is only active
C       with openmp enabled by compiler
C
C$    IF (.not.(DIRNLO.OR.DIRTRF).and.DIRSCF.and.(intomp.ne.0)) THEN
C$      IF (intomp.EQ.1) THEN
C$              CALL efp2ei_omp_jk_mode1(typscf,schwrz,nint,nschwz,
C$   +                               l1,l2,1,xints,
C$   +                               nsh2,maxg,
C$   +                               ia,S,VAO,db,fb,dsh,nflmat,
C$   +                               nshell,ktype,cutoff,.true.)
C$
C$      ELSEIF (intomp.EQ.2) THEN
C$
C$              CALL efp2ei_omp_kl_mode1(typscf,schwrz,nint,nschwz,
C$   +                                   l1,l2,1,xints,
C$   +                                   nsh2,maxg,
C$   +                                   ia,S,VAO,db,fb,dsh,nflmat,
C$   +                                   nshell,ktype,cutoff,out)
C$      ELSE
C$!       Should never get here, just in case
C$        IF(MASWRK) WRITE(IW,'(" INCONSISTENT OPENMP PARAMETERS")')
C$        IF(MASWRK) WRITE(IW,'(" ABORTING CALCULATION")')
C$        CALL ABRT
C$      END IF
C$
C$      GOTO 950 
C$    END IF


C END COMPUTING QM/EFP OVERLAP
      NOTPK = .NOT.PK
      NINT  = 0
      NSCHWZ= 0
      SCHSKP=.FALSE.
      DENMAX = ZERO
C
C        NOW WE ARE READY TO LOOP OVER ALL NSHELL**4 SHELL QUARTETS
C
C     ----- I SHELL -----
C
      DO 920 II = IST,NSHELL
C
C     ----- CHECK CPU TIME -----
C
      CALL TSECND(TIM)
      IF(TIM.GE.TIMLIM) THEN
C        NOTHING CAN BE DONE FOR IN-CORE INTEGRALS: JUST FORGET THEM
c         IF(.NOT.CMBDIR)
c     *      CALL FINAL(0,II,1,1,1,PANDK,BUFP,BUFK,IX,NINTMX)
         IF(MASWRK) WRITE(IW,9030)
         CALL ABRT
      END IF
C
C         ELONGATION METHOD
      IF(LTRMST) THEN
         IF(II.EQ.NHTSHL) THEN
            NFLTRM = NFILE
            NRCTRM = NREC
            NPSTRM = ICOUNT
         ENDIF
      ENDIF
C
C     ----- PRINT INTERMEDIATE RESTART DATA -----
C
      IF(NPRINT.NE.-5  .AND.  .NOT.CMBDIR .AND. MASWRK) THEN
         IF(ICOUNT.LE.NINTIC) THEN
            WRITE(IW,9015) II,JST,KST,LST,ICOUNT
         ELSE
            WRITE(IW,9010) II,JST,KST,LST,NREC,ICOUNT-NINTIC
         ENDIF
      ENDIF
C
C     ----- SKIP I SHELL IF NOT SYMMETRY UNIQUE -----
C     THIS, AND THE SIMILAR BRANCHINGS FOR THE J, K, AND L LOOPS IS
C     WHAT GENERATES THE "PETITE" RATHER THAN "GRANDE" INTEGRAL LIST.
C
      IF(C1GRP) THEN
         MI(1)=II
      ELSE
         DO 120 IT = 1,NT
            ID = MAPSHL(II,IT)
            IF (ID .GT. II) GO TO 920
            MI(IT) = ID
  120    CONTINUE
      END IF
C
C     ----- GO PARALLEL! -----
C
      IF(SLB) THEN
         IPCOUNT = IPCOUNT + 1
         IF (MOD(IPCOUNT,NPROC).NE.0) GO TO 920
      END IF
C     ----- J SHELL -----
C
      J0 = JST
      DO 900 JJ = J0,II
      JST = 1
C
      IF(C1GRP) THEN
         MJ(1)=JJ
      ELSE
         DO 200 IT = 1,NT
            ID = MI(IT)
            JD = MAPSHL(JJ,IT)
            MJ(IT) = JD
            IF (ID .GE. JD) GO TO 160
            ND = ID
            ID = JD
            JD = ND
  160       IF (ID-II) 200,180,900
  180       IF (JD-JJ) 200,200,900
  200    CONTINUE
      END IF
C
C     ----- GO PARALLEL! -----
C
      IF (DLB.AND.(.NOT.DLB_EFAI)) THEN
         MINE = MINE + 1
         IF (MINE.GT.NEXT) CALL DDI_DLBNEXT(NEXT)
         IF (NEXT.NE.MINE) GO TO 900
      END IF
C
C     ----- K SHELL -----
C
      K0 = KST
      DO 880 KK = K0,JJ
      KST = 1
C
      IF(C1GRP) THEN
         MK(1)=KK
      ELSE
         DO 340 IT = 1,NT
            ID = MI(IT)
            JD = MJ(IT)
            KD = MAPSHL(KK,IT)
            MK(IT) = KD
  240       IF (ID .GE. JD) GO TO 260
            ND = ID
            ID = JD
            JD = ND
  260       IF (JD .GE. KD) GO TO 280
            ND = JD
            JD = KD
            KD = ND
            GO TO 240
  280       IF (ID-II) 340,300,880
  300       IF (JD-JJ) 340,320,880
  320       IF (KD-KK) 340,340,880
  340    CONTINUE
      END IF
C
C     ----- GO PARALLEL! -----
C
!      IF(SLB) THEN
!         IPCOUNT = IPCOUNT + 1
!         IF (MOD(IPCOUNT,NPROC).NE.0) GO TO 880
!      END IF
C
C     ----- L SHELL ----
C
      L0 = LST
      DO 860 LL = L0,KK
      LST = 1
C
      IF(C1GRP) THEN
         M0(1)=1
         N4=1
      ELSE
         N4 = 0
         DO 540 IT = 1,NT
            ID = MI(IT)
            JD = MJ(IT)
            KD = MK(IT)
            LD = MAPSHL(LL,IT)
  380       IF (ID .GE. JD) GO TO 400
            ND = ID
            ID = JD
            JD = ND
  400       IF (JD .GE. KD) GO TO 420
            ND = JD
            JD = KD
            KD = ND
            GO TO 380
  420       IF (KD .GE. LD) GO TO 440
            ND = KD
            KD = LD
            LD = ND
            GO TO 400
  440       IF (ID-II) 540,460,860
  460       IF (JD-JJ) 540,480,860
  480       IF (KD-KK) 540,500,860
  500       IF (LD-LL) 540,520,860
  520       N4 = N4+1
            M0(N4) = IT
  540    CONTINUE
      END IF
C
C         THE LOOP STRUCTURE IN THIS ROUTINE IS DESIGNED TO FACILITATE
C         SUPERMATRIX CONSTRUCTION BY HAVING UP TO THREE "EXCHANGED"
C         QUARTETS AVAILABLE AT ONCE.  THE LOOP STRUCTURE TO GENERATE
C         A MORE NORMAL CANONICAL ORDERING OF THE QUARTETS HITS THE
C         SAME QUARTETS IN A SLIGHTLY DIFFERENT ORDER, BUT BOTH LOOPS
C         WILL DO EXACTLY THE SAME QUARTETS.
C
C             CANONICAL                      SUPERMATRIX
C         DO ISH=1,NSHELL                 DO II=1,NSHELL
C           DO JSH=1,ISH                    DO JJ=1,II
C             IJSH = IA(ISH)+JSH
C             DO KSH=1,ISH                    DO KK=1,JJ
C               DO LSH=1,KSH                    DO LL=1,KK
C                 KLSH=IA(KSH)+LSH
C                 IF(IJSH.LT.KLSH),               [II JJ|KK LL],
C                    CYCLE KSH LOOP               [II KK|JJ LL],
C                 [ISH JSH|KSH LSH]               [II LL|JJ KK]
C               ENDDO                           ENDDO
C             ENDDO                           ENDDO
C           ENDDO                           ENDDO
C         ENDDO                           ENDDO
C
C     ----- CHECK FOR REDUNDANIES BETWEEN THE 3 COMBINATIONS -----
C            (II,JJ//KK,LL), (II,KK//JJ,LL), (II,LL//JJ,KK)
C
      SKIPA =  JJ.EQ.KK
      SKIPB = (II.EQ.KK) .OR. (JJ.EQ.LL)
      SKIPC = (II.EQ.JJ) .OR. (KK.EQ.LL)
      NPSYM = .FALSE.
      IF (SKIPA .OR. SKIPB .OR. SKIPC) GO TO 720
      NPSYM = .TRUE.
      DO 640 M = 1,N4
         IT = M0(M)
         IH = MI(IT)
         JH = MJ(IT)
         IF(JH.LE.IH) THEN
            ID = IH
            JD = JH
         ELSE
            ID = JH
            JD = IH
         END IF
         IF(.NOT.SKIPA) SKIPA = (ID.EQ.II .AND. JD.EQ.KK) .OR.
     *                          (ID.EQ.JJ .AND. JD.EQ.LL)
         IF(.NOT.SKIPB) SKIPB = (ID.EQ.II .AND. JD.EQ.LL) .OR.
     *                          (ID.EQ.JJ .AND. JD.EQ.KK)
         IF (SKIPA .AND. SKIPB) GO TO 660
         KH = MK(IT)
         IF(KH.LE.IH) THEN
            ID = IH
            KD = KH
         ELSE
            ID = KH
            KD = IH
         END IF
         IF(.NOT.SKIPC) SKIPC = (ID.EQ.II .AND. KD.EQ.LL) .OR.
     *                          (ID.EQ.JJ .AND. KD.EQ.KK)
         IF(SKIPA .AND. SKIPC) GO TO 680
         IF(SKIPB .AND. SKIPC) GO TO 700
  640 CONTINUE
      GO TO 720
C
  660 SKIPC = .TRUE.
      GO TO 720
  680 SKIPB = .TRUE.
      GO TO 720
  700 SKIPA = .TRUE.
C
C        GENERATE SYMMETRY FACTOR -Q4- FOR THIS QUARTET IN PETITE LIST
C
  720 CONTINUE
      Q4 = NT
      Q4 = Q4 / N4
C
C     ----- (II,JJ//KK,LL) -----
C
      IEXCH = 1
      ISH = II
      JSH = JJ
      KSH = KK
      LSH = LL
      QQ4 = Q4
      IF(SKIPA .AND. NPSYM) QQ4 = QQ4+Q4
      IF(SKIPB .AND. NPSYM) QQ4 = QQ4+Q4
      GO TO 780
C
C     ----- (II,KK//JJ,LL) -----
C
  740 IF (SKIPA) GO TO 760
      IEXCH = 2
      ISH = II
      JSH = KK
      KSH = JJ
      LSH = LL
      QQ4 = Q4
      IF (SKIPC .AND. NPSYM) QQ4 = QQ4+Q4
      GO TO 780
C
C     ----- (II,LL//JJ,KK) -----
C
  760 IF (SKIPB .OR. SKIPC) GO TO 840
      IEXCH = 3
      ISH = II
      JSH = LL
      KSH = JJ
      LSH = KK
      QQ4 = Q4
C
C        ----- COMPUTE TWO-ELECTRON INTEGRALS ----
C
  780 CONTINUE
      IF(PK .AND. IEXCH.EQ.1) CALL ZPKOUT(ISH,JSH,KSH,LSH,GHONDO,
     *                                    SKIPA,SKIPB,SKIPC,NPSYM)
C
C     APPLY THE SCHWARZ INEQUALITY TO SCREEN OUT SMALL INTEGRALS,
C      (II,JJ//KK,LL) .LE.  SQRT( (II,JJ//II,JJ)*(KK,LL//KK,LL) )
C     SEE, FOR EXAMPLE, J.L.WHITTEN, J.CHEM.PHYS. 58,4496-4501(1973)
C
      IF(SCHWRZ) THEN
         IJIJ = (ISH*ISH-ISH)/2 + JSH
         KLKL = (KSH*KSH-KSH)/2 + LSH
         TEST = QQ4*XINTS(IJIJ)*XINTS(KLKL)
         IF(DIRSCF) THEN
            DENMAX = SCHWDN(DSH,ISH,JSH,KSH,LSH,IA)
            TEST = TEST*DENMAX
         END IF
         SCHSKP = TEST.LT.CUTOFF
         IF(SCHSKP) NSCHWZ = NSCHWZ + 1
      END IF
      IF(SCHSKP) GO TO 820
C
C        ----- ELECTRON REPULSION INTEGRAL CALCULATION -----
C     THIS MAY USE ROTATED AXIS, ERIC, OR RYS QUADRATURE METHODS
C
      CALL SHELLQUART(ISH,JSH,KSH,LSH,GHONDO)
C
C        USE THE INTEGRALS JUST FORMED.  AT MOST, 1 OF THESE IS CALLED
C
C MODE=1 IS QM ONLY, NO EFP CONTRIBUTION
      IF(DIRSCF) THEN
         CALL DIRFCK_EFP(TYPSCF,IA,S,VAO,GHONDO,L1,NINT,
     *        DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,IDUMMY,IDUMMY,MODE)

      ELSE IF(DIRCIS) THEN
         CALL DRFCIS(DA,FA,DB,FB,GHONDO,L1,NINT,NFLMAT)
      ELSE IF(DIRNLO) THEN
         CALL DFCKNS(DNLO,FNLO,GHONDO,L1,NINT,NFLMAT)
      ELSE IF(DIRTRF) THEN
         CALL DIRTRN(BUFP,IX,NINTMX,GHONDO,ICONT,NINT)
      ELSE
         IF(NOTPK) CALL QOUT(BUFP,IX,NINTMX,GHONDO)
      END IF
C
  820 CONTINUE
      GO TO (740,760,840),IEXCH
C
C     ----- WRITE THE P (OR PK) SUPERMATRIX TO DISK FILE -----
C
  840 CONTINUE
C
C      IF(PK) CALL PKFILE(II,JJ,KK,LL,SKIPA,SKIPB,SKIPC,NPSYM,
C     *                   BUFP,BUFK,IX,NINTMX,GHONDO)
C
C     ----- END OF SHELL LOOPS -----
C
  860 CONTINUE
  880 CONTINUE
  900 CONTINUE
  920 CONTINUE
C
  950 CONTINUE !OMP

      JMO = JMO + NORBEF(NJ)
      JAO = JAO + NUMEF(JM)
      JNAT = JNAT + NATEF(JM)

 1000 CONTINUE
      IF (GOPARR) CALL DDI_GSUMF(1500,VAO,L2)
      IF(DLB) CALL DDI_DLBRESET
C
C     ----- OUTPUT THE LAST BITS OF INTEGRALS -----
C
      IF(DIRTRF) CALL ONEIDX(BUFP,IX,ICONT)
c      IF(.NOT.CMBDIR) THEN
c         IF(SCHWRZ) THEN
c            IF(GOPARR) CALL DDI_GSUMI(1055,NSCHWZ,1)
c            IF(NPRINT.NE.-5 .AND. MASWRK) WRITE(IW,9020) NSCHWZ
c         END IF
c         CALL FINAL(1,II,II,II,II,PANDK,BUFP,BUFK,IX,NINTMX)
c      END IF
C
C END OF PARALLEL GET VAO
       CALL DSCAL(L2,HALF,VAO,1)
       II = 0
       DO I = 1,L1
          II = II + I
          VAO(II) = VAO(II) + VAO(II)
       END DO

C P. XU: triangular VAO --> square VAOSQ
       CALL CPYTSQ(VAO,VAOSQ,L1,1)
C TRANSFORM VAO TO MO BASIS
      CALL TFSQU(VMO,VAOSQ,VEC,WRK2,L1,L1)

      DO I=COREA,NA
       DO M=COREA,NA
C P. XU: need to multiply by TWO to agree with the formula
        VMI(M,I)=TWO*VMO(M,I)
       END DO
      END DO

      GOPARR = GPSAVE
      ICOUNT =  ICNTSV
      NREC   =  NRECSV
      INTLOC =  INTLSV
      RETURN
C
 9010 FORMAT(1X,'II,JST,KST,LST =',4I3,' NREC =',I10,' INTLOC =',I5)
 9015 FORMAT(1X,'II,JST,KST,LST =',4I3,' IN CORE, INTLOC =',I12)
 9020 FORMAT(1X,'SCHWARZ INEQUALITY TEST SKIPPED',I12,
     *        ' INTEGRAL BLOCKS.')
 9030 FORMAT(//1X,'*** THIS JOB HAS EXHAUSTED ITS CPU TIME ***'/
     *         1X,'     (WHILE COMPUTING 2E- INTEGRALS)'///)
      END
C*MODULE EFTEI   *DECK EFP2EI_MODE2
C>    @brief    form fock matrix for QM/EFP EXCHANG REPULSION from TEI 
C>              MODE 2
C>    @author   Peng Xu and Tosaporn Sattasathuchana
C>            - 2018
C>    @details  Calculate fock matrix for QM/EFP EXCHANG REPULSION 
C>              on the fly where 2 center are in QM region
C>    @date     Sep 2019   Peng Xu and Tosaporn Sattasathuchana
C>            - adding OpenMP parallel
C>    @param VEC: MO coefficient of QM
C>    @param COREA: number of QM core orbitals
C>    @param MXBF: maximum basis functions of EFP
C>    @param PROVEC: EFP LMO coefficients
C>    @param VMI: FOCK contribution of QM/EFP EXREP in MO basis
C>    @param VAO2: TEI FOCK CONTRIBUTION USED FOR (ij|kl) in AO
C>                    basis where, i and k are EFP centers.
C>    @param DPRIME: density of EFP LMO
C>    @param ABAB: diagonal terms of VAO2
C>    @param EIJIJ: energy vector after diagonalization  
C>    @see SUBROUTINE EFP2EI in eftei.src
C>    
      SUBROUTINE EFP2EI_MODE2(MODE,TYPSCF,DIRSCF,DIRNLO,DIRTRF,DIRCIS,
     *                  SCHWRZ,NINT,NSCHWZ,L1,L2,BUFP,BUFK,IX,NINTMX,
     *                  GHONDO,MAXG,DA,FA,DB,FB,DNLO,FNLO,NFLMAT,
     *                  XINTB,NSHELT,PROVEC,MXBF,VEC,NA,VIJ,WRK,
     *                  VMI,COREA,EXR1,VAO2,DPRIME,ABAB,EIJIJ)
C
C$    USE params,ONLY: intomp
C$    USE efteiomp

      USE MX_LIMITS,ONLY:MXSHEF,MXGEFP,MXDFG,MXFGPT,MXFRG,MXDPPT,mxsh,
     *     mxgtot,mxao
      USE comm_PAULMO
      use comm_FRGTYP
      USE comm_EFPBAS
      USE comm_FRGINF
      USE EFAIEXREP, only : EFAITOL
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      LOGICAL DIRSCF,DIRNLO,DIRTRF,DIRCIS,CMBDIR
      LOGICAL SCHWRZ,GOPARR,DSKWRK,MASWRK,DLB,SLB,SCHSKP

      LOGICAL DLB_EFAI
      LOGICAL PK,PANDK,NOTPK,BLOCK,GPSAVE,OUT
      LOGICAL SDSKWRK
C
      DIMENSION DPRIME(MXBF,MXBF)
      DIMENSION PROVEC(MXBF,NTMO), VEC(L1,L1)
      DIMENSION WRK(L1),VMO(L2),VIJ(NA,NA),VAO2(L1,L1)
      DIMENSION VMI(L1,L1),ABAB(L2),EIJIJ(NA)
      INTEGER COREA
C
      DIMENSION BUFP(NINTMX),BUFK(NINTMX),IX(*),
     *          GHONDO(MAXG),DA(L2),FA(L2),DB(L2),FB(L2),
     *          DNLO(L1,L1),FNLO(L1,L1)
      DIMENSION XINTB(NSHELT,NSHELL)
C

      PARAMETER (HUND=100.0D+00)
      PARAMETER (ZERO=0.0D+00)
C
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON /NSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),
     *                CF(MXGTOT),CG(MXGTOT),CH(MXGTOT),CI(MXGTOT),
     *                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),
     *                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
      COMMON /OUTPUT/ NPRINT,ITOL,ICUT,NORMF,NORMP,NOPK
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /PKFIL / PK,PANDK,BLOCK
      COMMON /RESTAR/ TIMLIM,IREST,NREC,INTLOC,IST,JST,KST,LST
      COMMON /SHLEXC/ NORGSH(3),NORGSP(3),IEXCH,NANGM,NGTH(4)
      COMMON /SHLG70/ ISH,JSH,KSH,LSH,IJKLXX(4)
      COMMON /SHLNOS/ QQ4,LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,
     *                MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL,
     *                NIJ,IJ,KL,IJKL
      COMMON /SHLT  / TOL,CUTOFF,ICOUNT,OUT
      COMMON /IJPAIR/ IA(MXAO)
C
      STOL=CUTOFF/HUND
      ICNTSV = ICOUNT
      NRECSV = NREC 
      INTLSV = INTLOC
      ICOUNT = 1
      NREC   = 0
      INTLOC = 1
      IDUMMY=0
      DUMMY=ZERO

      TIM = ZERO
      CALL TSECND(TIM)
C
      CMBDIR= DIRSCF .OR. DIRNLO .OR. DIRTRF .OR. DIRCIS
C
C        THE OLD FASHIONED PARALLEL INTEGRAL TRANSFORMATIONS DO NOT
C        ALLOW THE AO INTEGRAL WORK TO BE RUN IN PARALLEL.  THE MODERN
C        DISTRIBUTED MEMORY TRANSFORMATIONS DO NOT PASS THROUGH HERE.
C
      ICONT=  0
      GPSAVE = GOPARR
      IF(DIRTRF) GOPARR=.FALSE.
C
C     ----- INITIALIZATION FOR PARALLEL WORK -----
C     BOTH STATIC AND DYNAMIC LOAD BALANCING ARE IMPLEMENTED BELOW
C
      SLB = GOPARR  .AND.  IBTYP.EQ.0
      DLB = GOPARR  .AND.  IBTYP.EQ.1
      DLB_EFAI= DLB .AND. NFRG.GT.10 .AND. L1.LE.100
C$    IF (DIRSCF.and.(intomp.ne.0)) DLB_EFAI=.false.
      NEXT = -1
      MINE = -1
      IPCOUNT = ME - 1
C
      NT1 = 1
C
      LMAX=0
      DO MJ = 1,NFRG
      JM = ISET(MJ)
      CALL BASEFP(LMAX2,JM)
      IF (LMAX2 .GT. LMAX) LMAX=LMAX2
      END DO
                    NANGM =  4
      IF(LMAX.EQ.2) NANGM =  6
      IF(LMAX.EQ.3) NANGM = 10
      IF(LMAX.EQ.4) NANGM = 15
      IF(LMAX.EQ.5) NANGM = 21
      IF(LMAX.EQ.6) NANGM = 28
      NGTH(4) = 1
      NGTH(3) = NGTH(4) * NANGM
      NGTH(2) = NGTH(3) * NANGM
      NGTH(1) = NGTH(2) * NANGM
      IF(NOPK.EQ.0) THEN
         NORGSH(1) = 0
         NORGSH(2) = NORGSH(1) + NANGM**4
         NORGSH(3) = NORGSH(2) + NANGM**4
         NORGSP(1) = 0
         NORGSP(2) = 256
         NORGSP(3) = 512
      ELSE
         DO I=1,3
            NORGSH(I) = 0
            NORGSP(I) = 0
         ENDDO
      END IF
C
      NOTPK = .NOT.PK
      NINT  = 0
      NSCHWZ= 0
      SCHSKP=.FALSE.
C
C        NOW WE ARE READY TO LOOP OVER ALL NSHELL**4 SHELL QUARTETS
C
C     ----- I SHELL -----
C
C
      JAO  = 0
      JMO  = 0
      JNAT = 0
      ISHLEF = 0
      CALL VCLR(VAO2,1,L1*L1)
      CALL VCLR(VIJ,1,NA*NA)
      DO 1000 MJ = 1,NFRG
      JM = ISET(MJ)
C
C     ----- GO PARALLEL! -----
C
      IF (DLB_EFAI) THEN
         MINE = MINE + 1
         IF (MINE.GT.NEXT) CALL DDI_DLBNEXT(NEXT)
         IF (NEXT.NE.MINE) GO TO 950
      END IF
      CALL VCLR(DPRIME,1,MXBF*MXBF)
       DO KSIGMA=1,NUMEF(JM)
        DO LAMBDA=1,NUMEF(JM)
        DO JEF=1,NORBEF(MJ)
          DPRIME(KSIGMA,LAMBDA)=DPRIME(KSIGMA,LAMBDA)+
     $                  (PROVEC(KSIGMA,JEF+JMO)*PROVEC(LAMBDA,JEF+JMO))
        END DO
       END DO
      END DO
C
C     CALL THREADED VERSION IF POSSIBLE
C
C       The code below is only active
C       with openmp enabled by compiler
C
C$    IF (.not.(DIRNLO.OR.DIRTRF).and.DIRSCF.and.(intomp.ne.0)) THEN
C$              CALL efp2ei_omp_jk_mode2(typscf,schwrz,nint,nschwz,
C$   +                               l1,l2,1,
C$   +                               nsh2,maxg,
C$   +                               ia,da,fa,db,fb,nflmat,
C$   +                               nshell,ktype,cutoff,.true.,
C$   +                               jm,ishlef,jnat,mxbf,nshelt,
C$   +                               dprime,vao2,xintb)
C$
C$    GOTO 950 
C$    END IF
C$
      DO 920 II = 1,NSHLEF(JM)
C
C     ----- GO PARALLEL! -----
C
      IF (DLB.AND.(.NOT.DLB_EFAI)) THEN
         MINE = MINE + 1
         IF (MINE.GT.NEXT) CALL DDI_DLBNEXT(NEXT)
         IF (NEXT.NE.MINE) GO TO 920
      END IF
      IF(SLB) THEN
         IPCOUNT = IPCOUNT + 1
         IF (MOD(IPCOUNT,NPROC).NE.0) GO TO 920
      END IF
C
C
C     ----- J SHELL -----
C
      DO 900 JJ = 1,NSHELL
C
C     ----- GO PARALLEL! -----
C
!      IF (DLB) THEN
!         MINE = MINE + 1
!         IF (MINE.GT.NEXT) CALL DDI_DLBNEXT(NEXT)
!         IF (NEXT.NE.MINE) GO TO 900
!      END IF
C
C     ----- K SHELL -----
C
C MODE 2 FOR (ij|kl) EFP-QM, k is efp
          KNSHELL=NSHLEF(JM)
      DO 880 KK = 1,KNSHELL
C
C     ----- GO PARALLEL! -----
C
!      IF(SLB) THEN
!         IPCOUNT = IPCOUNT + 1
!         IF (MOD(IPCOUNT,NPROC).NE.0) GO TO 880
!      END IF
C
C     ----- L SHELL ----
C
          LNSHELL=NSHELL
      DO 860 LL = 1,LNSHELL
      N4=1
C
      Q4 = NT1
      Q4 = Q4 / N4
C
C     ----- (II,JJ//KK,LL) -----
C
      IEXCH = 1
      ISH = II
      JSH = JJ
      KSH = KK
      LSH = LL
      QQ4 = Q4
C      GO TO 780
C
C        ----- COMPUTE TWO-ELECTRON INTEGRALS ----
C
C  780 CONTINUE
C      IF(PK .AND. IEXCH.EQ.1) CALL ZPKOUT(ISH,JSH,KSH,LSH,GHONDO,
C     *                                    .FALSE.,.FALSE.,
C     *                                    .FALSE.,.FALSE.)
C
C     APPLY THE SCHWARZ INEQUALITY TO SCREEN OUT SMALL INTEGRALS,
C      (II,JJ//KK,LL) .LE.  SQRT( (II,JJ//II,JJ)*(KK,LL//KK,LL) )
C     SEE, FOR EXAMPLE, J.L.WHITTEN, J.CHEM.PHYS. 58,4496-4501(1973)
C
      SCHSKP = .FALSE.
      IF(SCHWRZ) THEN
C     P. XU: Schwarz inequality screening
         TEST = XINTB(II+ISHLEF,JJ)*XINTB(KK+ISHLEF,LL)
C
         SCHSKP = TEST .LT. CUTOFF
      END IF
      IF(SCHSKP) GO TO 820
C
      CALL SHELLQUART_EFP(MODE,ISH,JSH,KSH,LSH,GHONDO,JM,JNAT)
C
C        USE THE INTEGRALS JUST FORMED.  AT MOST, 1 OF THESE IS CALLED
C
         CALL DIRFCK_EFP(TYPSCF,IA,DA,FA,GHONDO,L1,NINT,
     *                   DUMMY,DPRIME,DUMMY,VAO2,DUMMY,MXBF,JM,MODE)
C
  820 CONTINUE
C
  860 CONTINUE
  880 CONTINUE
  900 CONTINUE
  920 CONTINUE
C
  950 CONTINUE ! OMP
C
C----P. XU: call FINAL_EFP inside fragment loop
C RESET ICOUNT=1 AND NREC=0 FOR NEW FRAGMENT
      ICOUNT= 1
      NREC=   0

      JAO = JAO + NUMEF(JM)
      JMO = JMO + NORBEF(MJ)
      JNAT = JNAT + NATEF(JM)
      ISHLEF = ISHLEF + NSHLEF(JM)

 1000 CONTINUE      
C
      IF (GOPARR) CALL DDI_GSUMF(1501,VAO2,L1*L1)
      IF(DLB) CALL DDI_DLBRESET
C
C-------- GETTING EXR1 ENERGY
      CALL VCLR(ABAB,1,L2)
      DO 550 MU=1,L1
        DO 500 NU=1,MU
        K= (MU*MU-MU)/2 + NU
        ABAB(K)= ABAB(K) + VAO2(MU,NU)
  500   CONTINUE
  550 CONTINUE
       CALL TFDIAG(ABAB,VEC,EIJIJ,WRK,L1,L2,L1,NA)
       EXR1= ZERO
       DO L= COREA,NA
       EXR1 = EXR1 +EIJIJ(L)
       ENDDO
C-------- END EXR1
C-------- GETTING VMI CONTRIBUTION
      CALL TFSQU(VIJ,VAO2,VEC,WRK,L1,NA)
      DO M = COREA, NA
        DO I = COREA, NA
        VMI(M,I)=VMI(M,I) + VIJ(M,I)
        ENDDO
      ENDDO
C--------- END VMI
      GOPARR = GPSAVE
      ICOUNT = ICNTSV
      NREC   = NRECSV
      INTLOC = INTLSV

      RETURN
C
 9030 FORMAT(//1X,'*** THIS JOB HAS EXHAUSTED ITS CPU TIME ***'/
     *         1X,'     (WHILE COMPUTING 2E- INTEGRALS)'///)
      END
C*MODULE EFTEI   *DECK EFP2EI_MODE3
C>    @brief    form fock matrix for QM/EFP EXCHANG REPULSION from TEI 
C>              MODE 3
C>    @author   Peng Xu and Tosaporn Sattasathuchana
C>            - 2018
C>    @details  Calculate fock matrix for QM/EFP EXCHANG REPULSION 
C>              on the fly where 3 center are in QM region
C>    @date     Sep 2019   Peng Xu and Tosaporn Sattasathuchana
C>            - adding OpenMP parallel
C>
C>    @param XINTB: largest exchange integrals of each shell between QM
C                   and EFP 
C>    @param XINTA: largest exchange integrals of each shell (QM only)
C>    @param NSHELT: total number of EFP shells
C>    @param PROVEC: EFP LMO coefficients
C>    @param MXBF: maximum basis functions of EFP
C>    @param VEC: MO coefficient of QM
C>    @param MXMO: maximum number of MO of EFP
C>    @param DQM: density of QM
C>    @param SMAT: OVERLAP B/W QM and EFP in AO basis
C>    @param TMAT: Kinetic energy integrals in AO basis
C>    @param SAB: OVERLAP B/W QM and EFP in MO basis
C>    @param SSTAR: backtransformation of SAB in AO basis
C>    @param SPRIME: backtransformation of SAB where dimension of
C>                   SPRIME is MO basis of QM x AO basis of EFP
C>    @param VAO: TEI FOCK CONTRIBUTION, details below
C>    @param VAO2: TEI FOCK CONTRIBUTION, details below
C>    @param VXR: FOCK contribution of QM/EFP EXREP in MO basis
C>    @param ETERM34:ENERGY OF 2(G=2J-K)
C>    @param COREA: number of QM core orbitals
C>    @see SUBROUTINE EFP2EI in eftei.src
C>       
      SUBROUTINE EFP2EI_MODE3(MODE,TYPSCF,DIRSCF,DIRNLO,DIRTRF,DIRCIS,
     *                  SCHWRZ,NINT,NSCHWZ,L1,L2,BUFP,BUFK,IX,NINTMX,
     *                  GHONDO,MAXG,DA,FA,DB,FB,DNLO,FNLO,NFLMAT,
     *                  XINTB,XINTA,NSHELT,PROVEC,MXBF,
     *                  VEC,MXMO,DQM,SSTAR,VAO,VAO2,WRK,VXR,COREA,
     *                  SPRIME,SAB,SMAT,TMAT,ETERM34,XINTASQ,VMITMP)
C
C     DRIVER FOR AB INITIO/EFP2 2E- INTEGRALS.  BASED ON TWOEI

C$    USE params, ONLY: intomp 
C$    USE efteiomp 

      USE MX_LIMITS,ONLY:MXSHEF,MXGEFP,MXDFG,MXFGPT,MXFRG,MXDPPT,mxsh,
     *     mxgtot,mxao,mxatm
      USE comm_PAULMO
      use comm_FRGTYP
      USE comm_EFPBAS
      USE comm_FRGINF
      USE EFAIEXREP, only : EFAITOL
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      LOGICAL DIRSCF,DIRNLO,DIRTRF,DIRCIS,CMBDIR
      LOGICAL SCHWRZ,GOPARR,DSKWRK,MASWRK,DLB,SLB,SCHSKP
      LOGICAL PK,PANDK,NOTPK,BLOCK,GPSAVE,OUT
      LOGICAL DLB_EFAI
C
      DIMENSION BUFP(NINTMX),BUFK(NINTMX),IX(*),
     *          GHONDO(MAXG),DA(*),FA(L2),DB(L2),FB(L2),
     *          DNLO(L1,L1),FNLO(L1,L1)
      DIMENSION XINTB(NSHELT,NSHELL),
     *          XINTA((NSHELL*NSHELL+NSHELL)/2),XINTASQ(NSHELL,NSHELL)
C
      INTEGER COREA
      DIMENSION DQM(L1,L1),VAO(MXBF,L1),SPRIME(NA,MXBF),
     *          SSTAR(L1,MXBF),VAO2(L1,L1)
      DIMENSION SMAT(L1,MXBF), TMAT(L1,MXBF)
      DIMENSION WRK(MXBF),
     *          SAB(NA,MXMO),VEC(L1,L1),
     *          PROVEC(MXBF,NTMO),VXR(L1,L1), VMITMP(L1,L1)
C
      PARAMETER (HUND=100.0D+00)
      PARAMETER (ZERO=0.0D+00)
      PARAMETER (TWO=2.0D+00)
C
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON /NSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),
     *                CF(MXGTOT),CG(MXGTOT),CH(MXGTOT),CI(MXGTOT),
     *                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),
     *                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
      COMMON /OUTPUT/ NPRINT,ITOL,ICUT,NORMF,NORMP,NOPK
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /PKFIL / PK,PANDK,BLOCK
      COMMON /RESTAR/ TIMLIM,IREST,NREC,INTLOC,IST,JST,KST,LST
      COMMON /SHLEXC/ NORGSH(3),NORGSP(3),IEXCH,NANGM,NGTH(4)
      COMMON /SHLG70/ ISH,JSH,KSH,LSH,IJKLXX(4)
      COMMON /SHLNOS/ QQ4,LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,
     *                MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL,
     *                NIJ,IJ,KL,IJKL
      COMMON /SHLT  / TOL,CUTOFF,ICOUNT,OUT
C
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /IJPAIR/ IA(MXAO)
C
      STOL=CUTOFF/HUND
      ICNTSV = ICOUNT
      NRECSV = NREC 
      INTLSV = INTLOC
      ICOUNT = 1
      NREC   = 0
      INTLOC = 1
      ETERM34 = ZERO
      DUMMY=0.D0
      IDUMMY= 0
      TIM = ZERO
      CALL VCLR(VMITMP,1,L1*L1)
      CALL TSECND(TIM)
C
C           transform XINTA to a square matrix XINTASQ

      IF(SCHWRZ) THEN
      CALL CPYTSQ(XINTA,XINTASQ,NSHELL,1)
      END IF
C
      CMBDIR= DIRSCF .OR. DIRNLO .OR. DIRTRF .OR. DIRCIS
C
C        THE OLD FASHIONED PARALLEL INTEGRAL TRANSFORMATIONS DO NOT
C        ALLOW THE AO INTEGRAL WORK TO BE RUN IN PARALLEL.  THE MODERN
C        DISTRIBUTED MEMORY TRANSFORMATIONS DO NOT PASS THROUGH HERE.
C
      ICONT=  0
      GPSAVE = GOPARR
      IF(DIRTRF) GOPARR=.FALSE.
C
C     ----- INITIALIZATION FOR PARALLEL WORK -----
C     BOTH STATIC AND DYNAMIC LOAD BALANCING ARE IMPLEMENTED BELOW
C
      SLB = GOPARR  .AND.  IBTYP.EQ.0
      DLB = GOPARR  .AND.  IBTYP.EQ.1
      DLB_EFAI= DLB .AND. NFRG.GT.10 .AND. L1.LE.100
C$    IF (DIRSCF.and.(intomp.ne.0)) DLB_EFAI=.false.
      NEXT = -1
      MINE = -1
      IPCOUNT = ME - 1
C
      NT1 = 1
C
      LMAX=0
      DO MJ = 1,NFRG
      JM = ISET(MJ)
      CALL BASEFP(LMAX2,JM)
      IF (LMAX2 .GT. LMAX) LMAX=LMAX2
      END DO
                    NANGM =  4
      IF(LMAX.EQ.2) NANGM =  6
      IF(LMAX.EQ.3) NANGM = 10
      IF(LMAX.EQ.4) NANGM = 15
      IF(LMAX.EQ.5) NANGM = 21
      IF(LMAX.EQ.6) NANGM = 28
      NGTH(4) = 1
      NGTH(3) = NGTH(4) * NANGM
      NGTH(2) = NGTH(3) * NANGM
      NGTH(1) = NGTH(2) * NANGM
      IF(NOPK.EQ.0) THEN
         NORGSH(1) = 0
         NORGSH(2) = NORGSH(1) + NANGM**4
         NORGSH(3) = NORGSH(2) + NANGM**4
         NORGSP(1) = 0
         NORGSP(2) = 256
         NORGSP(3) = 512
      ELSE
         DO I=1,3
            NORGSH(I) = 0
            NORGSP(I) = 0
         ENDDO
      END IF
C
      NOTPK = .NOT.PK
      NINT  = 0
      NSCHWZ= 0
      SCHSKP=.FALSE.

      CALL VCLR(DQM,1,L1*L1)
      CALL VCLR(VAO2,1,L1*L1)

      DO LAMBDA=1,L1
       DO KSIGMA=1,L1
        DO K=COREA,NA
          DQM(LAMBDA,KSIGMA)=DQM(LAMBDA,KSIGMA)+
     $                    (VEC(LAMBDA,K)*VEC(KSIGMA,K))
        END DO
       END DO
      END DO
C
C        NOW WE ARE READY TO LOOP OVER ALL NSHELL**4 SHELL QUARTETS
C
C     ----- I SHELL -----
C
C
      JAO  = 0
      JMO  = 0
      JNAT = 0
      ISHLEF = 0
      DO 1000 MJ = 1,NFRG
      JM = ISET(MJ)
      CALL VCLR(VAO,1,MXBF*L1)
      CALL VCLR(SMAT,1,L1*MXBF)
      CALL VCLR(SAB,1,NA*MXMO)
      CALl VCLR(SSTAR,1,L1*MXBF)
      CALL VCLR(SPRIME,1,NA*MXBF)
C
C     ----- GO PARALLEL! -----
C
      IF (DLB_EFAI) THEN
         MINE = MINE + 1
         IF (MINE.GT.NEXT) CALL DDI_DLBNEXT(NEXT)
         IF (NEXT.NE.MINE) GO TO 960
      ENDIF

      CALL STINT(SMAT,TMAT,MXBF,L1,NAT,NATEF(JM),
     *           NGSSEF(JM),NSHLEF(JM),EXEF(1,JM),CSEF(1,JM),CPEF(1,JM),
     *           CDEF(1,JM),CFEF(1,JM),CGEF(1,JM),KSTREF(1,JM),
     *           KATMEF(1,JM),
     *           KTYPEF(1,JM),KNGEF(1,JM),KLOCEF(1,JM),KMINEF(1,JM),
     *           KMAXEF(1,JM),
     *           L1,NSHELL,EX,CS,
     *           CP,CD,CF,CG,KSTART,
     *           KATOM,KTYPE,KNG,KLOC,
     *           KMIN,KMAX,PRCORD(1,JNAT+1),C)
C          SMAT ADDRESSING CHANGED FEB 2010
      CALL TFSQP(SAB,SMAT,VEC,PROVEC(1,JMO+1),WRK,
     *           L1,NUMEF(JM),NA,NORBEF(MJ),NA,NORBEF(MJ),.FALSE.)

      DO JZETA=1,NUMEF(JM)
       DO M=COREA,NA
        DO J=1,NORBEF(MJ)
         SPRIME(M,JZETA)=SPRIME(M,JZETA)+
     *            SAB(M,J)*PROVEC(JZETA,J+JMO)
        END DO
       END DO
      END DO
      DO JZETA=1,NUMEF(JM)
       DO LAMBDA=1,L1
        DO J=1,NORBEF(MJ)
         DO K=COREA,NA
         SSTAR(LAMBDA,JZETA)=SSTAR(LAMBDA,JZETA)+
     *        VEC(LAMBDA,K)*(SAB(K,J)*PROVEC(JZETA,J+JMO))
         END DO
        END DO
       END DO
      END DO
C
C     CALL THREADED VERSION IF POSSIBLE
C
C       The code below is only active
C       with openmp enabled by compiler
C
C$    IF (.not.(DIRNLO.OR.DIRTRF).and.DIRSCF.and.(intomp.ne.0)) THEN
C$              CALL efp2ei_omp_jk_mode3(typscf,schwrz,nint,nschwz,
C$   +                               l1,l2,1,xinta,
C$   +                               nsh2,maxg,
C$   +                               ia,da,fa,db,fb,nflmat,
C$   +                               nshell,ktype,cutoff,.true.,
C$   +                               jm,ishlef,jnat,mxbf,nshelt,
C$   +                               dqm,vao,vao2,sstar,xintasq,xintb)
C$
C$    GOTO 950
C$    END IF
C$
      DO 920 II = 1,NSHLEF(JM)
C
C
C     ----- GO PARALLEL! -----
C
      IF (DLB.AND.(.NOT.DLB_EFAI)) THEN
         MINE = MINE + 1
         IF (MINE.GT.NEXT) CALL DDI_DLBNEXT(NEXT)
         IF (NEXT.NE.MINE) GO TO 920
      END IF
      IF(SLB) THEN
         IPCOUNT = IPCOUNT + 1
         IF (MOD(IPCOUNT,NPROC).NE.0) GO TO 920
      END IF
C
C     ----- J SHELL -----
C
C     J0 = JST
      DO 900 JJ = 1,NSHELL
C
C     ----- GO PARALLEL! -----
C
!      IF (DLB) THEN
!         MINE = MINE + 1
!         IF (MINE.GT.NEXT) CALL DDI_DLBNEXT(NEXT)
!         IF (NEXT.NE.MINE) GO TO 900
!      END IF
C
C     ----- K SHELL -----
C
C MODE 3 FOR (ij|kl) EFP-QM, only i is efp.
      DO 880 KK = 1,NSHELL
C
C     ----- GO PARALLEL! -----
C
!      IF(SLB) THEN
!         IPCOUNT = IPCOUNT + 1
!         IF (MOD(IPCOUNT,NPROC).NE.0) GO TO 880
!      END IF
C
C     ----- L SHELL ----
C
      DO 860 LL = 1,KK
      N4=1
C
C         THE LOOP STRUCTURE IN THIS ROUTINE IS DESIGNED TO FACILITATE
C         SUPERMATRIX CONSTRUCTION BY HAVING UP TO THREE "EXCHANGED"
C         QUARTETS AVAILABLE AT ONCE.  THE LOOP STRUCTURE TO GENERATE
C         A MORE NORMAL CANONICAL ORDERING OF THE QUARTETS HITS THE
C         SAME QUARTETS IN A SLIGHTLY DIFFERENT ORDER, BUT BOTH LOOPS
C         WILL DO EXACTLY THE SAME QUARTETS.
C
C             CANONICAL                      SUPERMATRIX
C         DO ISH=1,NSHELL                 DO II=1,NSHELL
C           DO JSH=1,ISH                    DO JJ=1,II
C             IJSH = IA(ISH)+JSH
C             DO KSH=1,ISH                    DO KK=1,JJ
C               DO LSH=1,KSH                    DO LL=1,KK
C                 KLSH=IA(KSH)+LSH
C                 IF(IJSH.LT.KLSH),               [II JJ|KK LL],
C                    CYCLE KSH LOOP               [II KK|JJ LL],
C                 [ISH JSH|KSH LSH]               [II LL|JJ KK]
C               ENDDO                           ENDDO
C             ENDDO                           ENDDO
C           ENDDO                           ENDDO
C         ENDDO                           ENDDO
C
      Q4 = NT1
      Q4 = Q4 / N4
C
C     ----- (II,JJ//KK,LL) -----
C
      IEXCH = 1
      ISH = II
      JSH = JJ
      KSH = KK
      LSH = LL
      QQ4 = Q4
      GO TO 780
C
C        ----- COMPUTE TWO-ELECTRON INTEGRALS ----
C
  780 CONTINUE
C
C     APPLY THE SCHWARZ INEQUALITY TO SCREEN OUT SMALL INTEGRALS,
C      (II,JJ//KK,LL) .LE.  SQRT( (II,JJ//II,JJ)*(KK,LL//KK,LL) )
C     SEE, FOR EXAMPLE, J.L.WHITTEN, J.CHEM.PHYS. 58,4496-4501(1973)
C
      SCHSKP = .FALSE.
      IF(SCHWRZ) THEN
C     P. XU: Schwarz inequality screening
         TEST = XINTB(II+ISHLEF,JJ)*XINTASQ(KK,LL)
         SCHSKP = TEST .LT. CUTOFF
      END IF
      IF(SCHSKP) GO TO 820
C
      CALL SHELLQUART_EFP(MODE,ISH,JSH,KSH,LSH,GHONDO,JM,JNAT)
C
C        USE THE INTEGRALS JUST FORMED.  AT MOST, 1 OF THESE IS CALLED
C
C FOR VAO: (-1/2)*S(m,j)*2*G(i,j) and (-1/2)*S(i,j)*2*G(m,j)
C FOR VAO2: -S(k,j)*[4(kj|mi)-(km|ji)-(ki|jm)] contribution
          CALL DIRFCK_EFP(TYPSCF,IA,DA,FA,GHONDO,L1,NINT,
     *                    DQM,DUMMY,VAO,VAO2,SSTAR,MXBF,JM,MODE)
C
  820 CONTINUE
C
  860 CONTINUE
  880 CONTINUE
  900 CONTINUE
  920 CONTINUE
C
  950 CONTINUE !OMP
  970 CONTINUE !sabtest
C
C----P. XU: call FINAL_EFP inside fragment loop
C
C RESET ICOUNT=1 AND NREC=0 FOR NEW FRAGMENT
      ICOUNT= 1
      NREC=   0
      !IF (GOPARR) CALL DDI_GSUMF(1502,VAO,MXBF*L1)
C ------------------- CALCULATE ENERGY OF 2(G=2J-K)
      DO JZETA=1,NUMEF(JM)
       DO LAMBDA=1,L1
         ETERM34= ETERM34 + SSTAR(LAMBDA,JZETA)*VAO(JZETA,LAMBDA)
       END DO
      END DO
C
C ------ COMPUTE VMITMP ------------------------------
       DO M=COREA,NA
        DO I=COREA,NA
        DO MU=1,L1
         DO JZETA=1,NUMEF(JM)
          VMITMP(M,I)=VMITMP(M,I)+(VAO(JZETA,MU)*
     *            SPRIME(M,JZETA)*VEC(MU,I))
         END DO
        END DO
        END DO
       END DO
       DO M=COREA,NA
        DO I=COREA,NA
        DO MALPHA=1,L1
         DO JZETA=1,NUMEF(JM)
          VMITMP(M,I)=VMITMP(M,I)+(VAO(JZETA,MALPHA)*
     *            SPRIME(I,JZETA)*VEC(MALPHA,M))
         END DO
        END DO
        END DO
       END DO
 960  CONTINUE      
      JAO = JAO + NUMEF(JM)
      JMO = JMO + NORBEF(MJ)
      JNAT = JNAT + NATEF(JM)
      ISHLEF = ISHLEF + NSHLEF(JM)

 1000 CONTINUE      
C
      !IF (GOPARR) CALL DDI_GSUMF(1503,VAO2,L1*L1)
      !IF(DLB) CALL DDI_DLBRESET
C
      !ETERM34= TWO*ETERM34
C
      DO M=COREA,NA
       DO I=COREA,NA
        DO MALPHA=1,L1
         DO MU=1,L1
         VMITMP(M,I)=VMITMP(M,I)+VAO2(MALPHA,MU)*
     *         VEC(MALPHA,M)*VEC(MU,I)
         END DO
        END DO
       END DO
      END DO
C
      IF (GOPARR) CALL DDI_GSUMF(1502,ETERM34,1)
      IF (GOPARR) CALL DDI_GSUMF(1503,VMITMP,L1*L1)
      ETERM34= TWO*ETERM34

      DO M = COREA,NA
      DO I = COREA,NA
      VXR(M,I) = VXR(M,I) + VMItmp(M,I)
      ENDDO
      ENDDO

      GOPARR = GPSAVE
      ICOUNT = ICNTSV
      NREC   = NRECSV
      INTLOC = INTLSV
      IF(DLB) CALL DDI_DLBRESET
      RETURN
C
 9030 FORMAT(//1X,'*** THIS JOB HAS EXHAUSTED ITS CPU TIME ***'/
     *         1X,'     (WHILE COMPUTING 2E- INTEGRALS)'///)
      END
C*MODULE EFTEI *DECK DIRFCK_EFP
C>    @brief    form fock matrix for QM/EFP EXCHANG REPULSION from TEI 
C>    
C>    @author   Peng Xu and Tosaporn Sattasathuchana
C>            - 2018
C>    @details  Calculate fock matrix for QM/EFP EXCHANG REPULSION 
C>              on the fly
C>    @param DQM    : density of QM 
C>    @param DPRIME : density of EFP LMO
C>    @param VAO    : TEI FOCK CONTRIBUTION USED IN MODE3, DETAILS BELOW
C>    @param VAO2   : TEI FOCK CONTRIBUTION USED IN MODE 2 and 3,
C>                    DETAILS BELOW
C>    @param SSTAR  : QM/EFP OVERLAP MATRIX IN AO BASIS 
C>    @param MXBF   : maximum basis functions of EFP
C>    @param JM     : COUNTER OF EFP TYPE 
C>    @param MODE   : MODE FOR QM/EFP TEI EXREP
C>    @see SUBROUTINE DIRFCK in scflib.src
C>
      SUBROUTINE DIRFCK_EFP(SCFTYP,IA,DA,FA,GHONDO,L1,NINT,
     *           DQM,DPRIME,VAO,VAO2,SSTAR,MXBF,JM,MODE)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      DIMENSION IA(*),DA(*),FA(*),GHONDO(*)
      DIMENSION DQM(L1,L1),VAO(MXBF,L1),SSTAR(L1,MXBF),VAO2(L1,L1)
      DIMENSION DPRIME(MXBF,MXBF)
C
      LOGICAL RHF,CPHF,UROHF,GVB,CPUHF
      LOGICAL OUT
C
      COMMON /DFTPAR/ DFTTYP(20),EXENA,EXENB,EXENC,
     *                IDFT34,NAUXFUN,NAUXSHL
      COMMON /ERIOUT/ ISH,JSH,KSH,LSH,LSTRI,LSTRJ,LSTRK,LSTRL
      COMMON /SHLT  / TOL,CUTOFFAO,ICOUNT,OUT
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
C
      DATA DHRHF,DHUHF,DHROHF /8HRHF     ,8HUHF     ,8HROHF     /
      DATA DHGVB /8HGVB     /
C
C     ----- FORM FOCK OPERATOR DIRECTLY FROM INTEGRALS -----
C     THIS ROUTINE WAS MODIFIED BASED ON DIRFCK IN SCFLIB.SRC
C
C     THE DFT CASES ARE NOT IMPLEMENT FOR QM/EFP YET.
C
C     DIRFCK_EF CALLS ONE OF THE FOLLOWING:
C         DIREFAI_MODE1 FOR 4-QM-CENTER TEI IN QM/EFP EXREP
C         DIREFAI_MODE2 FOR 2-QM-CENTER TEI IN QM/EFP EXREP
C         DIREFAI_MODE3 FOR 3-QM-CENTER TEI IN QM/EFP EXREP
C
C        SCFTYP    [IN] SCF TYPE 8 BYTE HOLLERITH CONSTANT.
C        IA(*)     [IN] TRIANGULAR INDEX ADDRESSES
C        DA(*)     [IN] RHF OR ROHF/UHF ALPHA DENSITY MATRIX.
C        FA(*)     [IN,OUT] RHF FOCK MATRIX.
C        GHONDO(*) [IN] HONDO INDEXED INTEGRALS.
C        L2        [IN] SIZE OF SYMMETRIC PACKED DENSITY OR FOCK MATRIX.
C        NINT      [IN,OUT] NUMBER OF NON-ZERO INTEGRALS.
C
      HFSCAL=DFTTYP(3)
      CSCALT=1.0D+00
C
C SEE Theor Chem Acc (2010) 125:481-491 for reference
C
      IF (MODE.EQ.1) THEN
C DA here is actually S^2 (overlap between QM and EFP), in AO basis
         CALL DIRFCK_RHF(IA,DA,FA,ISH,JSH,KSH,LSH,
     &                   GHONDO,LSTRI,LSTRJ,LSTRK,LSTRL,
     &                   HFSCAL,CSCALT,CUTOFFAO,NINT)
      ELSE IF(MODE.EQ.2) THEN
C SGO approximation is removed, so we compute the following integrals in
C this mode2
C  VAO2: (ij|ij) where i is efp and j is qm.

         CALL DIREFAI_MODE2(IA,DPRIME,VAO2,ISH,JSH,KSH,LSH,
     &                   GHONDO,LSTRI,LSTRJ,LSTRK,LSTRL,
     &                   HFSCAL,CSCALT,CUTOFFAO,NINT,JM,MXBF,L1)
      ELSE IF(MODE.EQ.3) THEN
C FOR VAO: (-1/2)*S(m,j)*2*G(i,j) and (-1/2)*S(i,j)*2*G(m,j)
C contribution
C FOR VAO2: -S(k,j)*[4(kj|mi)-(km|ji)-(ki|jm)] contribution
         CALL DIREFAI_MODE3(IA,DQM,VAO,SSTAR,VAO2,ISH,JSH,KSH,LSH,
     &                   GHONDO,LSTRI,LSTRJ,LSTRK,LSTRL,
     &                   HFSCAL,CSCALT,CUTOFFAO,NINT,JM,MXBF,L1)
      ENDIF
C
      RETURN
      END
C*MODULE EFTEI *DECK DIREFAI_MODE1
C>    @brief    Fork contributions of QM/EFP Exchange
C>              repulsion for 3-QM-center TEI  
C>    
C>    @author   Peng Xu and Tosaporn Sattasathuchana
C>            - 2018
C>    @details  Calculate fock matrix on the fly. This subroutine is
C>              modified from SUBROUTINE DIRFCK_RHF in scflib.src

C June 20, 2018
C THIS IS THE CLONE OF DIRFCK_RHF SUBROUTINE.
      SUBROUTINE DIREFAI_MODE1(IA,DA,FA,ISH,JSH,KSH,LSH,
     &                      GHONDO,ISTRIDE,JSTRIDE,KSTRIDE,LSTRIDE,
     &                      HFSCAL,CSCALT,CUTOFF,NINT,LL1)
      use mx_limits, only: mxsh,mxgtot
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      DIMENSION IA(*),DA(*),FA(*),GHONDO(*)
C
      LOGICAL IANDJ,KANDL,SAME
C
C
      COMMON /NSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),
     *                CF(MXGTOT),CG(MXGTOT),CH(MXGTOT),CI(MXGTOT),
     *                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),
     *                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
C
      PARAMETER (HALF=0.5D+00)
      PARAMETER (FOUR=4.0D+00)
C
C     COMPUTES FOCK MATRIX ELEMENTS OF RHF OR RHF/DFT WAVEFUNCTION.
C     SHOULD NOT BE CALLED DIRECTLY, USE THE -DIRFCK- INTERFACE.
C        IA(*)                  [IN] TRIANGULAR INDEX ARRAY
C        DA(*)                  [IN] DENSITY MATRIX.
C        FA(*)                  [IN,OUT] FOCK MATRIX.
C        ISH,JSH,KSH,LSH        [IN] SHELL INDICES
C        GHONDO(*)              [IN] HONDO INDEXED INTEGRALS.
C        ISTRIDE,JSTRIDE,KSTRIDE,LSTRIDE -
C             [IN] HONDO INDEXED INTEGRALS SHELL STRIDES.
C        HFSCAL,CSCALT          [IN] DFT SCALING FACTORS.
C        CUTOFF                 [IN] NON-ZERO INTEGRAL CUTOFF
C        NINT                   [IN,OUT] NUMBER OF NONZERO INTEGRALS
C
      SAME  = ISH.EQ.KSH.AND.JSH.EQ.LSH
      IANDJ = ISH.EQ.JSH
      KANDL = KSH.EQ.LSH
C
      MINI = KMIN(ISH)
      MINJ = KMIN(JSH)
      MINK = KMIN(KSH)
      MINL = KMIN(LSH)
      MAXI = KMAX(ISH)
      MAXJ = KMAX(JSH)
      MAXK = KMAX(KSH)
      MAXL = KMAX(LSH)
      LOCI = KLOC(ISH)-MINI
      LOCJ = KLOC(JSH)-MINJ
      LOCK = KLOC(KSH)-MINK
      LOCL = KLOC(LSH)-MINL
C
      XVAL1 = HFSCAL
      XVAL4 = FOUR*CSCALT
C
      NIJ = 0
      MAXJ2 = MAXJ
      I_INDEX = 1
      DO 360 I = MINI,MAXI
         IF (IANDJ) MAXJ2 = I
C
         I1 = I+LOCI
C
         IJ_INDEX = I_INDEX
         I_INDEX = I_INDEX + ISTRIDE
C
         DO 340 J = MINJ,MAXJ2
            NIJ = NIJ+1
            MAXL2 = MAXL
C
            J1 = J+LOCJ
            I2 = I1
            J2 = J1
            IF (I1.LT.J1) THEN ! SORT <IJ|
               I2 = J1
               J2 = I1
            ENDIF
C
            IJK_INDEX = IJ_INDEX
            IJ_INDEX = IJ_INDEX + JSTRIDE
C
            NKL = NIJ
C
            DO 320 K =  MINK,MAXK
               IF (KANDL) MAXL2 = K
C
               K1 = K + LOCK
C
               IF(SAME) THEN ! ACCOUNT FOR NON-UNIQUE PERMUTATIONS
                  ITMP = MIN(MAXL2-MINL+1,NKL)
                  IF (ITMP.EQ.0) GOTO 340
                  MAXL2 = MINL + ITMP - 1
                  NKL = NKL - ITMP
               ENDIF
C
               IJKL_INDEX = IJK_INDEX
               IJK_INDEX = IJK_INDEX + KSTRIDE
C
               DO 300 L=MINL,MAXL2
C
                  VAL = GHONDO( IJKL_INDEX )
                  IJKL_INDEX = IJKL_INDEX + LSTRIDE
                  IF(ABS(VAL).LT.CUTOFF) GOTO 300
                  NINT = NINT + 1
C
                  L1 = L + LOCL
                  K2 = K1
                  L2 = L1
C
                  IF (K2.LT.L2) THEN ! SORT |KL>
                     K2 = L1
                     L2 = K1
                  ENDIF
C
                  II = I2
                  JJ = J2
                  KK = K2
                  LL = L2
C
                  IF (II.LT.KK) THEN ! SORT <IJ|KL>
                     II = K2
                     JJ = L2
                     KK = I2
                     LL = J2
                  ELSE IF (II.EQ.KK.AND.JJ.LT.LL) THEN ! SORT <IJ|IL>
                     JJ = L2
                     LL = J2
                  ENDIF
C
                  II2 = IA(II)
                  JJ2 = IA(JJ)
                  KK2 = IA(KK)
C
                  IJ = II2 + JJ
                  IK = II2 + KK
                  IL = II2 + LL
                  JK = JJ2 + KK
                  JL = JJ2 + LL
                  KL = KK2 + LL
                  IF (JJ.LT.KK) JK = KK2 + JJ
                  IF (JJ.LT.LL) JL = IA(LL) + JJ
C
C       ACCOUNT FOR IDENTICAL PERMUTATIONS.
C
                  IF(II.EQ.JJ) VAL = VAL*HALF
                  IF(KK.EQ.LL) VAL = VAL*HALF
                  IF(II.EQ.KK.AND.JJ.EQ.LL) VAL = VAL*HALF
                  VAL1 = VAL*XVAL1
                  VAL4 = VAL*XVAL4
C
                  FA(IJ) = FA(IJ) + VAL4*DA(KL)
                  FA(KL) = FA(KL) + VAL4*DA(IJ)
                  FA(IK) = FA(IK) - VAL1*DA(JL)
                  FA(JL) = FA(JL) - VAL1*DA(IK)
                  FA(IL) = FA(IL) - VAL1*DA(JK)
                  FA(JK) = FA(JK) - VAL1*DA(IL)
 300           ENDDO
 320        ENDDO
 340     ENDDO
 360  ENDDO
C
      RETURN
      END
C END OF MODE=1
C*MODULE EFTEI  *DECK DIREFAI_MODE2
C>    @brief    Fork contributions of QM/EFP Exchange
C>              repulsion for 2-QM-center TEI  
C>    
C>    @author   Peng Xu and Tosaporn Sattasathuchana
C>            - 2018
C>    @details  Calculate fock matrix on the fly. This subroutine is
C>              modified from SUBROUTINE DIRFCK_RHF in scflib.src
C>              This term was previously approximated by SGO. This
C>              subroutine calculates exact integral and provides better
C>              EXREP results.
C>    @param DA   : density of EFP LMO
C>    @param FA   : TEI FOCK CONTRIBUTION USED IN MODE 2
C>    @param MXBF : maximum basis functions of EFP
C>    @param NINT : NUMBER OF NONZERO INTEGRALS
C>    @param JM   : COUNTER OF EFP TYPE 
C>    @param HFSCAL,CSCALT: SCALING FACTORS.
C>    @param CUTOFF: NON-ZERO INTEGRAL CUTOFF
C>    @param ISH,JSH,KSH,LSH: SHELL INDICES
C>    @param GHONDO: HONDO INDEXED INTEGRALS.
C>    @param ISTRIDE,JSTRIDE,KSTRIDE,LSTRIDE: HONDO INDEXED INTEGRALS 
C>                                            SHELL STRIDES.
C>    @see SUBROUTINE DIRFCK_RHF in scflib.src
C>
      SUBROUTINE DIREFAI_MODE2(IA,DA,FA,ISH,JSH,KSH,LSH,
     &                      GHONDO,ISTRIDE,JSTRIDE,KSTRIDE,LSTRIDE,
     &                      HFSCAL,CSCALT,CUTOFF,NINT,JM,MXBF,LL1)
C
      USE MX_LIMITS,ONLY:MXDFG,MXFGPT,MXSHEF,MXGEFP,mxsh,mxgtot
      USE comm_EFPBAS
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      DIMENSION IA(*),DA(MXBF,MXBF),FA(LL1,LL1),GHONDO(*)
C
      LOGICAL IANDJ,KANDL,SAME
C
      PARAMETER (HALF=0.5D+00)
      PARAMETER (FOUR=4.0D+00)
C
      COMMON /NSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),
     *                CF(MXGTOT),CG(MXGTOT),CH(MXGTOT),CI(MXGTOT),
     *                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),
     *                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
C
C----------------------------------------------------------------------
C THERE IS NO SYMMETRY FOR (IJ|IJ) 
C I IS EFP CENTER, AND J IS QM CENTER.
C      SAME  = ISH.EQ.KSH.AND.JSH.EQ.LSH
C      IANDJ = ISH.EQ.JSH
C      KANDL = KSH.EQ.LSH
      SAME=.FALSE.
      IANDJ=.FALSE.
      KANDL=.FALSE.
C
      MINI = KMINEF(ISH,JM)
      MINJ = KMIN(JSH)
      MINK = KMINEF(KSH,JM)
      MINL = KMIN(LSH)
      MAXI = KMAXEF(ISH,JM)
      MAXJ = KMAX(JSH)
      MAXK = KMAXEF(KSH,JM)
      MAXL = KMAX(LSH)
      LOCI = KLOCEF(ISH,JM)-MINI
      LOCJ = KLOC(JSH)-MINJ
      LOCK = KLOCEF(KSH,JM)-MINK
      LOCL = KLOC(LSH)-MINL
C
      XVAL1 = HFSCAL
      XVAL4 = FOUR*CSCALT
C
      NIJ = 0
      MAXJ2 = MAXJ
      I_INDEX = 1
      DO 360 I = MINI,MAXI
C         IF (IANDJ) MAXJ2 = I
C
         I_INDEX = (I-MINI)*ISTRIDE + 1
C
         DO 340 J = MINJ,MAXJ2
            NIJ = NIJ+1
            MAXL2 = MAXL
C
            IJ_INDEX = (J-MINJ)*JSTRIDE + I_INDEX
C
             NKL = 0
C
            DO 320 K =  MINK,MAXK
               IJK_INDEX = (K-MINK)*KSTRIDE + IJ_INDEX
C               IF (KANDL) MAXL2 = K
C
               DO 300 L=MINL,MAXL2
                  NKL = NKL+1
                  IJKL_INDEX = (L-MINL)*LSTRIDE + IJK_INDEX
C
                  VAL = GHONDO( IJKL_INDEX )
                  IF(ABS(VAL).LT.CUTOFF) GOTO 300
                  NINT = NINT + 1
C
                  I1 = I + LOCI
                  J1 = J + LOCJ
                  K1 = K + LOCK
                  L1 = L + LOCL
C
                  VAL1 = VAL*XVAL1
                  VAL4 = VAL*XVAL4
C
                  FA(J1,L1) =FA(J1,L1) - VAL*DA(I1,K1)

 300           ENDDO
 320        ENDDO
 340     ENDDO
 360  ENDDO
C
      RETURN
      END    
C*MODULE EFTEI  *DECK DIREFAI_MODE3
C>    @brief    Fork contributions of QM/EFP Exchange
C>              repulsion for 3-QM-center TEI  
C>    
C>    @author   Peng Xu and Tosaporn Sattasathuchana
C>            - 2018
C>    @details  Calculate fock matrix on the fly. This subroutine is
C>              modified from SUBROUTINE DIRFCK_RHF in scflib.src
C>    @param DA   : density of QM
C>    @param FA   : TEI FOCK CONTRIBUTION, details below
C>    @param DB   : density of EFP LMO
C>    @param FB   : TEI FOCK CONTRIBUTION, details below
C>    @param MXBF : maximum basis functions of EFP
C>    @param NINT : NUMBER OF NONZERO INTEGRALS
C>    @param JM   : COUNTER OF EFP TYPE 
C>    @param HFSCAL,CSCALT: SCALING FACTORS.
C>    @param CUTOFF: NON-ZERO INTEGRAL CUTOFF
C>    @param ISH,JSH,KSH,LSH: SHELL INDICES
C>    @param GHONDO: HONDO INDEXED INTEGRALS.
C>    @param ISTRIDE,JSTRIDE,KSTRIDE,LSTRIDE: HONDO INDEXED INTEGRALS 
C>                                            SHELL STRIDES.
C>    @see SUBROUTINE DIRFCK_RHF in scflib.src
C>
      SUBROUTINE DIREFAI_MODE3(IA,DA,FA,DB,FB,ISH,JSH,KSH,LSH,
     &                      GHONDO,ISTRIDE,JSTRIDE,KSTRIDE,LSTRIDE,
     &                      HFSCAL,CSCALT,CUTOFF,NINT,JM,MXBF,LL1)
C
      USE MX_LIMITS,ONLY:MXDFG,MXFGPT,MXSHEF,MXGEFP,mxsh,mxgtot
      USE comm_EFPBAS
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      DIMENSION IA(*),DA(LL1,LL1),FA(MXBF,LL1),GHONDO(*)
     &          ,DB(LL1,MXBF),FB(LL1,LL1)
C
      LOGICAL IANDJ,KANDL,SAME
C
      PARAMETER (HALF=0.5D+00)
      PARAMETER (FOUR=4.0D+00)
C
      COMMON /NSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),
     *                CF(MXGTOT),CG(MXGTOT),CH(MXGTOT),CI(MXGTOT),
     *                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),
     *                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL

C
C
C I IS EFP CENTER, THE RESTS ARE QM CENTERS.
C      SAME  = ISH.EQ.KSH.AND.JSH.EQ.LSH
C      IANDJ = ISH.EQ.JSH
      KANDL = KSH.EQ.LSH
      SAME=.FALSE.
      IANDJ=.FALSE.
C
      MINI = KMINEF(ISH,JM)
      MINJ = KMIN(JSH)
      MINK = KMIN(KSH)
      MINL = KMIN(LSH)
      MAXI = KMAXEF(ISH,JM)
      MAXJ = KMAX(JSH)
      MAXK = KMAX(KSH)
      MAXL = KMAX(LSH)
      LOCI = KLOCEF(ISH,JM)-MINI
      LOCJ = KLOC(JSH)-MINJ
      LOCK = KLOC(KSH)-MINK
      LOCL = KLOC(LSH)-MINL
C
      XVAL1 = HFSCAL
      XVAL4 = FOUR*CSCALT
C
      NIJ = 0
      MAXJ2 = MAXJ
      I_INDEX = 1
      DO 360 I = MINI,MAXI
         IF (IANDJ) MAXJ2 = I
C
         I1 = I+LOCI
C
         IJ_INDEX = I_INDEX
         I_INDEX = I_INDEX + ISTRIDE
C
         DO 340 J = MINJ,MAXJ2
            NIJ = NIJ+1
            MAXL2 = MAXL
C
            J1 = J+LOCJ
            I2 = I1
            J2 = J1
C
            IJK_INDEX = IJ_INDEX
            IJ_INDEX = IJ_INDEX + JSTRIDE
C
            NKL = NIJ
C
            DO 320 K =  MINK,MAXK
               IF (KANDL) MAXL2 = K
C
               K1 = K + LOCK
C
C               IF(SAME) THEN ! ACCOUNT FOR NON-UNIQUE PERMUTATIONS
C                  ITMP = MIN(MAXL2-MINL+1,NKL)
C                  IF (ITMP.EQ.0) GOTO 340
C                  MAXL2 = MINL + ITMP - 1
C                  NKL = NKL - ITMP
C               ENDIF
C
               IJKL_INDEX = IJK_INDEX
               IJK_INDEX = IJK_INDEX + KSTRIDE
C
               DO 300 L=MINL,MAXL2
C
                  VAL = GHONDO( IJKL_INDEX )
                  IJKL_INDEX = IJKL_INDEX + LSTRIDE
                  IF(ABS(VAL).LT.CUTOFF) GOTO 300
                  NINT = NINT + 1
C
                  L1 = L + LOCL
                  K2 = K1
                  L2 = L1
C
                  II = I2
                  JJ = J2
                  KK = K2
                  LL = L2
C
C       ACCOUNT FOR IDENTICAL PERMUTATIONS.
C
                  IF(KK.EQ.LL) VAL = VAL*HALF

                  VAL1 = VAL*XVAL1
                  VAL4 = VAL*XVAL4
C FOR FA: (-1/2)*S(m,j)*2*G(i,j) and (-1/2)*S(i,j)*2*G(m,j)
         FA(II,JJ)=FA(II,JJ)-DA(LL,KK)*VAL4
C
         FA(II,LL)=FA(II,LL)+DA(JJ,KK)*VAL1
C
         FA(II,KK)=FA(II,KK)+DA(JJ,LL)*VAL1
C FOR FB: -S(k,j)*[4(kj|mi)-(km|ji)-(ki|jm)] contribution
         FB(LL,KK)=FB(LL,KK)-DB(JJ,II)*VAL4
         FB(KK,LL)=FB(KK,LL)-DB(JJ,II)*VAL4

         FB(LL,JJ)=FB(LL,JJ)+DB(KK,II)*VAL1
         FB(KK,JJ)=FB(KK,JJ)+DB(LL,II)*VAL1

         FB(JJ,KK)=FB(JJ,KK)+DB(LL,II)*VAL1
         FB(JJ,LL)=FB(JJ,LL)+DB(KK,II)*VAL1
                  
 300           ENDDO
 320        ENDDO
 340     ENDDO
 360  ENDDO
C
      RETURN
      END

