C*MODULE QUANPOC  *DECK FFMD1
!>
!> @brief    pure MM MD
!>
!> @author   Nandun Thellamurege, Hui Li
!>           - Jan 2011
!>
!> @details  MD simulation for MM atoms
!>
      SUBROUTINE FFMD1(ATMNAM,CORD,CORDSV,ZANF,
     *                 ZMAS,ONEMAS,QMZMAS,QM1MAS,
     *                 CHARG,POL,POLSV,DIP,
     *                 FIELD1,FIELD2,FIELD3,
     *                 SIG,EPS,SIG2,EPS2,
     *                 BOND0,FCBOND,FCSTBD,
     *                 ANGL0,FCANGL,FCWAGG,
     *                 DIHB0,FCDIHB,FCDIHR,
     *                 VROT,NNN,GAMA,IPAIR,
     *                 KLIST,LLIST,KBLST,MLIST,NLIST,
     *                 L1213J,L14J,
     *                 VEL,QMVEL,FCLJTP,NTYPE,
     *                 FFGRD0,FFGRD1,FFGRD2,
     *                 XTS,YTS,ZTS,CMAT1,
     *                 POT1,POT2,QRXN1,QRXN2,NTS,
     *                 LISTQM,NONLS1,NONLSTQ,MAPLST,CMAPCO,
     *                 LSTCELL,NONLS2,CORDSV2,CORDSVQ,
     *                 MVFASTS2,MVFASTS3,MVFASTS4,
     *                 MVFASTL2,MVFASTL3,MVFASTL4,
     *                 AFIX,QFIX,
     *                 RFIX,IDATOM,DAI,IDDAI,
     *                 VFIX1,VFIX2,XTSFIX,YTSFIX,ZTSFIX,
     *                 NONLSA,NONLSB,L1213A,L1213B,L14A,L14B,
     *                 CHARGB,
     *                 SIGB,EPSB,SIG2B,EPS2B,
     *                 OLDCORD,LSTRAT,DSTRAT,NONLSPMA,
     *                 L1213PMA,L14PMA,CORDB,
     *                 LSBONDPMA,LSANGLPMA,
     *                 LSDIHRPMA,LSDIHBPMA,
     *                 LSWAGGPMA,LSCMAPPMA,
     *                 NONLSPMB,L1213PMB,L14PMB,
     *                 LSBONDPMB,LSANGLPMB,
     *                 LSDIHRPMB,LSDIHBPMB,
     *                 LSWAGGPMB,LSCMAPPMB,
     *                 UMBHIS,UM2HIS,VELSV,DFSC0,
     *                 CORDG,CORDGSV,CORDGSV2,CORDGSVQ,
     *                 FFGRDG0,FFGRDG1,FFGRDG2)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      LOGICAL GOPARR,DSKWRK,MASWRK
C
      PARAMETER (TOANGS=0.52917724924D+00)
      PARAMETER (TOBOHR=1.0D+00/TOANGS)
      PARAMETER (TOKCAL=627.509469D+00)
      PARAMETER (TOHART=1.0D+00/TOKCAL)
      PARAMETER (PT5=0.5D+00)
      PARAMETER (TWOTHIRD=2.0D+00/3.0D+00)
      PARAMETER (ONETHIRD=1.0D+00/3.0D+00)
      PARAMETER (ONESIX=1.0D+00/6.0D+00)
      PARAMETER (FIVESIX=5.0D+00/6.0D+00)
C
      CHARACTER*10  ATMNAM
C
      DIMENSION ATMNAM(NFFAT),CORD(3,NFFAT),CORDSV(3,NFFAT),ZANF(NFFAT),
     *          ZMAS(NFFAT),ONEMAS(NFFAT),QMZMAS(*),QM1MAS(*),
     *          CHARG(NFFAT),POL(NFFAT),DIP(3,NFFAT),
     *          FIELD1(3,NFFAT),FIELD2(3,NFFAT),FIELD3(3,NFFAT),
     *          SIG(NFFAT),EPS(NFFAT),SIG2(NFFAT),EPS2(NFFAT),
     *          BOND0(NBOND),FCBOND(NBOND),
     *          ANGL0(NANGL),FCANGL(NANGL),
     *          FCWAGG(NWAGG),DIHB0(NDIHB),FCDIHB(NDIHB),FCDIHR(3,*),
     *          VROT(NDIHR),NNN(NDIHR),GAMA(NDIHR),IPAIR(2,NBOND),
     *          KLIST(3,NANGL),
     *          LLIST(4,NDIHR),MLIST(4,NWAGG),
     *          L1213J(2,*),L14J(2,NDIHR),
     *          VEL(3,NFFAT),QMVEL(3,*),NTYPE(*),
     *          FFGRD0(3,NFFAT),FFGRD1(3,NFFAT),FFGRD2(3,NFFAT),
     *          XTS(NTS),YTS(NTS),ZTS(NTS),CMAT1(NTS,NTS),
     *          POT1(NTS),POT2(NTS),QRXN1(NTS),QRXN2(NTS),NONLS1(2,*),
     *          NONLSTQ(*),MAPLST(6,*),CMAPCO(4,4,24,24,3),
     *          CORDB(3,*),OLDCORD(3,*),DFSC0(3,*),CORDSV2(3,*)
      DIMENSION TIMSTR(3)
      DIMENSION ENALL(100)
      DIMENSION VELSV(3,*)
      DIMENSION IMDCATM(20020)
C
      COMMON /FFDFS / TIMDFS,QDION,AMION,TEFF,NDFS,NATMGAS,
     *                LFFDFSC,
     *                LFFDFSC0,LFFDFSA,LFFDFSN,LFFDFCOM,KDFS,LFFDFSCAV
      COMMON /FFENGY/ EN12,EN123,EN123R4,EN123B4,EN234W1,ENCHAR,ENLJR,
     *                ENLJD,ENPOL,XENPOL,ENRXN,XENRXN,ENRXNPOL,ENRXNR,
     *                EN12312,ENQUANP(30),
     *                ENBIAS,ENCENT,ENUCCH,ENCMAP,ENPOT,ENKIN,ENTOT
      COMMON /FFEWLD/ SPLIT,IEWALD,KEWALD,NKVEC,
     *                L1KVEC,L2KVEC,LFFRKEXPEL,LFFRKVEC,
     *                LFFKVEC,LFFTCHCH,LFFCOSCH,LFFSINCH
      COMMON /FFFIXS/ ENFIXSO,FIXEPS,FIXTOL,FIXA,FIXQ,RALLMM,RALLQM,
     *                RADMM(200),RADQM(200),NRADMM,NRADQM,IFIXSOL,
     *                LFFDAI,LFFDAIT,LFFIDDAI,LFFIDTMP,LFFTMPTS,
     *                LFFAFIX,LFFIDATOM,LFFRFIX,LFFQFIX,NTSATM,
     *                LFFQFIXMP,LFFQFIXTA,LFFQFIXXY,
     *                LFFXTSFIX,LFFYTSFIX,LFFZTSFIX,
     *                LFFVFIX1,LFFVFIX2,NCYCLE,MXFFTS,NFFTS
      COMMON /FFFREE/ SOL1CH,SOL1LJ,SOL1IM,SOLFRE1MM,
     *                SOL2CH,SOL2LJ,SOL2IM,SOLFRE2MM,
     *                ASOL1CH,ASOL1LJ,ASOL1IM,ASOL1MM,
     *                ASOL2CH,ASOL2LJ,ASOL2IM,ASOL2MM,
     *                PMF1BD,PMF1AG,PMF1DR,PMF1DB,PMF1WG,PMF1CM,
     *                PMF1CH,PMF1PO,PMF1LJ,PMF1IM,PMFFRE1MM,
     *                EN12A,EN12B,EN123A,EN123B,EN123R4A,EN123R4B,
     *                EN234W1A,EN234W1B,EN123B4A,EN123B4B,
     *                EN12312A,EN12312B,
     *                ENCHARA,ENCHARB,ENPOLA,ENPOLB,
     *                ENRXNA,ENRXNB,ENRXNPOLA,ENRXNPOLB,
     *                ENRXNRA,ENRXNRB,ENLJRA,ENLJRB,ENLJDA,ENLJDB,
     *                ENCMAPA,ENCMAPB,ENFIXSOA,ENFIXSOB,
     *                ENCENTA,ENCENTB,WSIMUL,WPERT1,WPERT2,
     *                IFEPTYP
      COMMON /FFFRE2/ N1FFAT,N1BOND,N1ANGL,N1DIHR,N1DIHB,N1CMAP,N1WAGG,
     *                N2FFAT,N2BOND,N2ANGL,N2DIHR,N2DIHB,N2CMAP,N2WAGG,
     *                LFF2ATMNAM,LFF2CORD,LFF2ZANF,LFF2ZMAS,
     *                LFF2ONEMAS,LFF2CHARG,LFF2POL,
     *                LFF2SIG,LFF2EPS,LFF2SIG2,LFF2EPS2,
     *                LFF2BOND0,LFF2FCBOND,LFF2ANGL0,LFF2FCANGL,
     *                LFF2FCWAGG,LFF2DIHB0,LFF2FCDIHB,
     *                LFF2VROT,LFF2NNN,LFF2GAMA,LFF2IPAIR,
     *                LFF2KLIST,LFF2LLIST,LFF2MLIST,LFF2NLIST,
     *                LFF2VEL,LFF2QMVEL,LFF2CLPR,LFF2ZLPR,
     *                LFF2NLPR,LFF2MAPLST,
     *                LFFLISTB2A,NTODOA,LFFNONLSA,NTODOB,LFFNONLSB,
     *                N1213A,LFFL1213A,N1213B,LFFL1213B,
     *                N14A,LFFL14A,N14B,LFFL14B,
     *                NTODOPMA,LFFNONLSPMA,NBONDPMA,LFFLSBONDPMA,
     *                NANGLPMA,LFFLSANGLPMA,NDIHRPMA,LFFLSDIHRPMA,
     *                NDIHBPMA,LFFLSDIHBPMA,NWAGGPMA,LFFLSWAGGPMA,
     *                NCMAPPMA,LFFLSCMAPPMA,
     *                N1213PMA,LFFL1213PMA,
     *                N14PMA,LFFL14PMA,
     *                NTODOPMB,LFFNONLSPMB,NBONDPMB,LFFLSBONDPMB,
     *                NANGLPMB,LFFLSANGLPMB,NDIHRPMB,LFFLSDIHRPMB,
     *                NDIHBPMB,LFFLSDIHBPMB,NWAGGPMB,LFFLSWAGGPMB,
     *                NCMAPPMB,LFFLSCMAPPMB,
     *                N1213PMB,LFFL1213PMB,
     *                N14PMB,LFFL14PMB
      COMMON /FFPARM/ NFFAT,NBOND,NANGL,NDIHR,NDIHB,NCMAP,NWAGG,
     *                N1213J,N14J,NLKQMM,IDOCHG,IDOPOL,IDOLJ,IDOCMAP
      COMMON /FFMDPA/ DT,DT2,TEMP0,PRES0,POLTOL,VIR(3),PMEAN,VOLAV,
     *                ENPAV,ENKAV,TEMPAV,BERENDT,BERENDP,VELMAX,
     *                PMEANX,PMEANY,PMEANZ,NSTEP,KMASTER,KOUTACT(2),
     *                IHESS,INTALG,ITSTAT,IPSTAT,JOUT,KOUT,LOUT
      COMMON /FFMDOP/ LFFCORDG,LFFCORDGSV,LFFCORDGSV2,LFFCORDGSVQ,
     *                LFFFFGRDG0,LFFFFGRDG1,LFFFFGRDG2,MDOPT
      COMMON /FFMPT3/ NACTMM,LACTMM(2020),LFFDIPOLD,JUMBUP,
     *                NACTQM,LACTQM(2020),LFFOLDC,LFFQMVELSV,MMHESS,
     *                LFFQMCHG,LFFQMCHGB,ISWAP,R2SWAP,DFTBMM
      COMMON /FFPBSW/ XBOX,YBOX,ZBOX,SWF,SWFDX,SWFDY,SWFDZ,
     *                SWRA,ONESWRA,SWRA2,ONESWRA2,
     *                SWRB,ONESWRB,SWRB2,ONESWRB2,ONESWRB4,
     *                SWFDUM3,SWFDUM4,SWFDUM5,
     *                SWRAQ,ONESWRAQ,SWRAQ2,ONESWRAQ2,
     *                SWRBQ,ONESWRBQ,SWRBQ2,ONESWRBQ2,ONESWRBQ4,
     *                SWFDUM3Q,SWFDUM4Q,SWFDUM5Q,
     *                QMSIZE,QMCX,QMCY,QMCZ,QMCXSV,QMCYSV,QMCZSV,
     *                CENTX,CENTY,CENTZ,BUFWID1,BUFWID2,RDAMP,
     *                EFIELDX,EFIELDY,EFIELDZ,QMCXSV2,QMCYSV2,QMCZSV2,
     *                EPS1RB,EPS1RB3,ONEXBOX,ONEYBOX,ONEZBOX,
     *                LQMCT,MXLIST1,NTODO,NTODOSV,NTODOQ,
     *                ISWITCH,ISHIFT,IPOLSHF,
     *                LFFLSTCELL,LFFCORDSV,
     *                LFFPOLSV,LFFCORDSV2,LFFNONLS2,LFFCORDSVQ,
     *                LFFMVFASTS2,LFFMVFASTS3,LFFMVFASTS4,
     *                LFFMVFASTL2,LFFMVFASTL3,LFFMVFASTL4,
     *                MXCHECK,MXLIST2,NTODO2,NTODO2SV
      COMMON /FFRATT/ RATOLC,RATOLV,SCALRAT,VIRRAT(3),IRATTLE,JRATTLE,
     *                NRATTLE,MXRATT,LFFOLDCORD,LFFLSTRAT,LFFDSTRAT,
     *                LFFVELSV,IRATQM
      COMMON /FFRMDF/ RDIST1,RDIST2,RINTRV,NRATM1,NRATM2,INTRSTP
      COMMON /FFRMSD/ DIPT(8),TIMGYRA,TIMRALL,
     *                NATPDB,NGYRA,NDIEL,NRALL,
     *                LFFRALL0,NRMSD,LFFRMSD0,KFREEAB(201),
     *                NRIJMM,IJRMM(2,100),NRIJQM,IJRQM(2,100),
     *                NAIJKMM,IJKMM(3,100),NAIJKQM,IJKQM(3,100),
     *                NFIXMM,IFIXMM(200),NFIXQM,IFIXQM(200)
      COMMON /FFRXN / RXNEPS,RSPHSOL,ISPHSOL
      COMMON /FFSPH / SPHRAD,SPHSIG,SPHEPS,IADDWAT
      COMMON /FFUMBR/ UMBFC,UMBR0,UMBSIZE,
     *                NUMBBIN,NUMBATM(6),NUMBTYP,LFFUMBHIS,
     *                UM2FC,UM2R0,UM2SIZE,
     *                NUM2BIN,NUM2ATM(6),NUM2TYP,LFFUM2HIS
      COMMON /IOFILE/ IR,IW,IP,IJK,IJKT,IDAF,NAV,IODA(950)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
C
      INTEGER, PARAMETER :: K15 = SELECTED_INT_KIND(15)
      INTEGER, PARAMETER :: MAGIC2 = TRANSFER(20170303193144_K15,1)
C
C     NANDUN THELLAMUREGE, HUI LI, JAN 2011, LINCOLN
C     HUI LI, MAR 2012
C     FENGCHAO CUI, HUI LI, JUN 2012 (ADD FREE ENERGY)
C
      IF(NSTEP.LT.0) RETURN
C
      CALL SEQOPN(36,'TRAJECT','NEW',.FALSE.,'FORMATTED')
      CALL SEQOPN(38,'MDDIP'  ,'NEW',.FALSE.,'FORMATTED')
      ISTEP = 0
      IF(IADDWAT.EQ.2 .AND. ISPHSOL.GT.0) RSPHSOL=RSPHSOL*1.0D+40
C
C     -- READY TO TAKE OFF ?
C
      MDSTEP = -1
      DO ISTEP = 0, NSTEP
         MDSTEP = MDSTEP + 1
         IF(ISTEP.EQ.100.AND.IADDWAT.EQ.2 .AND. ISPHSOL.GT.0)
     *   RSPHSOL = RSPHSOL*1.0D-40
         IF(MDSTEP.LE.1) CALL VCLR(UMBHIS,1,NUMBBIN)
         IF(MDSTEP.LE.1) CALL VCLR(UM2HIS,1,NUM2BIN*NUMBBIN)
C
         IF(MDSTEP.EQ.0) THEN
            IF(MDOPT.GT.0) THEN
               IFATM  = 0
               ICHOSE = 0
               ION    = 0
               DO IFFAT = 1,NFFAT
                  NUCZ = INT(ZANF(IFFAT)+0.001D+00)
                  IF(NUCZ.EQ.  3) THEN
                     ION  = IFFAT
                     RCH0 = 1.8D+00*TOBOHR
                  END IF
                  IF(NUCZ.EQ. 11) THEN
                     ION  = IFFAT
                     RCH0 = 2.3D+00*TOBOHR
                  END IF
                  IF(NUCZ.EQ. 19) THEN
                     ION  = IFFAT
                     RCH0 = 2.7D+00*TOBOHR
                  END IF
                  IF(NUCZ.EQ. 37) THEN
                     ION  = IFFAT
                     RCH0 = 2.8D+00*TOBOHR
                  END IF
                  IF(NUCZ.EQ. 55) THEN
                     ION  = IFFAT
                     RCH0 = 3.0D+00*TOBOHR
                  END IF
                  IF(NUCZ.EQ.7.OR.NUCZ.EQ.8.OR.NUCZ.EQ.16) THEN
                     IFATM= IFATM+1
                     IF(IFATM.GE.20000) THEN
                        IF(MASWRK) WRITE(IW,*)
     *                  'ERROR: IMDCATM(IFATM) EXCEEDED 20000.'
                        CALL ABRT
                     END IF
                     IMDCATM(IFATM)=IFFAT
                  END IF
               ENDDO
               IFSTEP=NSTEP/MAX(1,IFATM)
            END IF
            GOTO 50
         END IF
C
         CALL DCOPY(3*NFFAT,CORD,1,OLDCORD,1)
C
C        - CALCULATE X(T+DT)
         IF(INTALG.EQ.1) THEN
            IF(MDSTEP.EQ.1) CALL DCOPY(3*NFFAT,FFGRD1,1,FFGRD0,1)
            DO IFFAT = 1, NFFAT
               DUM  = ONEMAS(IFFAT)*DT2
               DO III = 1, 3
                  CHANGE = +VEL(III,IFFAT)*DT
     *                     -TWOTHIRD*FFGRD1(III,IFFAT)*DUM
     *                     +ONESIX*FFGRD0(III,IFFAT)*DUM
                  CORD(III,IFFAT)=CORD(III,IFFAT) + CHANGE
               ENDDO
            ENDDO
         END IF
         IF(INTALG.EQ.2) THEN
            DO IFFAT = 1, NFFAT
               DUM  = PT5*ONEMAS(IFFAT)*DT2
               DO III = 1, 3
                  CHANGE = +VEL(III,IFFAT)*DT
     *                     -FFGRD1(III,IFFAT)*DUM
                  CORD(III,IFFAT)=CORD(III,IFFAT) + CHANGE
               ENDDO
            ENDDO
         END IF
         IF(NACTMM.GT.0) THEN
            DO KOPT = 1, NACTMM
               IFFAT = LACTMM(KOPT)
               IF(IFFAT.GT.0) THEN
                  CORD(1,IFFAT)=(CORD(1,IFFAT)+3.1416D+00)*1.0D-40
                  CORD(2,IFFAT)=(CORD(2,IFFAT)+2.7183D+00)*1.0D-40
                  CORD(3,IFFAT)=(CORD(3,IFFAT)+1.4427D+00)*1.0D-40
               END IF
            ENDDO
            DO IFFAT = 1, NFFAT
               IF(ABS(CORD(1,IFFAT))+ABS(CORD(2,IFFAT))
     *           +ABS(CORD(3,IFFAT)).LT.1.0D-20) THEN
                  CORD(1,IFFAT) = CORD(1,IFFAT)*1.0D+40 - 3.1416D+00
                  CORD(2,IFFAT) = CORD(2,IFFAT)*1.0D+40 - 2.7183D+00
                  CORD(3,IFFAT) = CORD(3,IFFAT)*1.0D+40 - 1.4427D+00
               ELSE
                  DO III = 1, 3
                     CORD(III,IFFAT) = OLDCORD(III,IFFAT)
                  ENDDO
               END IF
            ENDDO
         END IF
         DO KFIX=1,NFIXMM
            DO III = 1, 3
               CORD(III,IFIXMM(KFIX)) = OLDCORD(III,IFIXMM(KFIX))
            ENDDO
         ENDDO
C
         CALL VCLR(VIRRAT,1,3)
         IF(NRATTLE.GT.0) THEN
C           - FIXED ATOMS DO NOT MOVE, SO MAKE THE MASS INFINITE
            IF(NACTMM.GT.0) THEN
               DO IFFAT = 1, NFFAT
                  ONEMAS(IFFAT) = 0.0D+00
               ENDDO
               DO KOPT = 1, NACTMM
                  IFFAT = LACTMM(KOPT)
                  IF(IFFAT.GT.0) ONEMAS(IFFAT) = 1.0D+00/ZMAS(IFFAT)
               ENDDO
            END IF
            DO KFIX = 1, NFIXMM
               KFFAT = IFIXMM(KFIX)
               IF(KFFAT.GT.0) ONEMAS(KFFAT) = 0.0D+00
            ENDDO
            CALL RATTLE1(CORD,OLDCORD,VEL,DSTRAT,ONEMAS,LSTRAT,
     *                   MDSTEP)
            DO KFIX = 1, NFIXMM
               KFFAT = IFIXMM(KFIX)
               IF(KFFAT.GT.0) ONEMAS(KFFAT) = 1.0D+00/ZMAS(KFFAT)
            ENDDO
            IF(NACTMM.GT.0) THEN
               DO IFFAT = 1, NFFAT
                  ONEMAS(IFFAT) = 1.0D+00/ZMAS(IFFAT)
               ENDDO
            END IF
         END IF
C
         IF(NRATM1.GT.0) THEN
            CALL IRMDF1(CORD,OLDCORD,VEL,ONEMAS,MDSTEP,RMDFCE)
         END IF
C
  50     CONTINUE
C
C        -- CALCULATE ENERGY AND GRADIENT AT X(T+DT)
         CALL NONBOND(MDSTEP,CORD,CORDSV,CORDSV2,CORDSVQ,
     *                NONLS1,NONLS2,
     *                NONLSTQ,LSTCELL,
     *                MVFASTS2,MVFASTS3,MVFASTS4,
     *                MVFASTL2,MVFASTL3,MVFASTL4,
     *                NONLSA,NONLSB,
     *                NONLSPMA,NONLSPMB)
         CALL VCLR(FFGRD2,1,3*NFFAT)
         CALL VCLR(VIR,1,3)
         CALL E00012(CORD,FFGRD2,BOND0,FCBOND,IPAIR,CORDB,
     *               LSBONDPMA,LSBONDPMB)
         CALL E00123(CORD,FFGRD2,ANGL0,FCANGL,KLIST,CORDB,
     *               LSANGLPMA,LSANGLPMB)
         CALL E12312(CORD,FFGRD2,ANGL0,KLIST,BOND0,FCSTBD,
     *               KBLST,CORDB,LSANGLPMA,LSANGLPMB)
         CALL E123B4(CORD,FFGRD2,DIHB0,FCDIHB,NLIST,CORDB,
     *               LSDIHBPMA,LSDIHBPMB)
         CALL E234W1(CORD,FFGRD2,      FCWAGG,MLIST,CORDB,
     *               LSWAGGPMA,LSWAGGPMB)
         CALL E123R4(CORD,FFGRD2,VROT,GAMA,NNN,LLIST,CORDB,
     *               LSDIHRPMA,LSDIHRPMB,FCDIHR)
         CALL ECMAP (CORD,FFGRD2,MAPLST,CMAPCO,CORDB,
     *               LSCMAPPMA,LSCMAPPMB)
         CALL ELJ126(CORD,FFGRD2,SIG,EPS,SIG2,EPS2,L14J,NONLS1,
     *               L1213J,SIGB,EPSB,SIG2B,EPS2B,
     *               NONLSA,NONLSB,L1213A,L1213B,L14A,L14B,
     *               NONLSPMA,L1213PMA,L14PMA,CORDB,
     *               NONLSPMB,L1213PMB,L14PMB,FCLJTP,NTYPE)
         CALL ESPHER(CORD,FFGRD2)
         CALL ECHARG(CORD,FFGRD2,CHARG,NONLS1,L1213J,L14J,
     *               CHARGB,NONLSA,NONLSB,L1213A,L1213B,L14A,L14B,
     *               NONLSPMA,L1213PMA,L14PMA,CORDB,
     *               NONLSPMB,L1213PMB,L14PMB)
C
C        -- MDOPT: APPLY FORCE BETWEEN METAL ION AND A CHOSEN ATOM
C
         IF(MDOPT.GT.0) THEN
         IF(MOD(ISTEP-1,IFSTEP).EQ.0) ICHOSE=ICHOSE+1
         IF(MDSTEP.GT.0.AND.ION.GT.0.AND.IFATM.GT.1
     *      .AND.ICHOSE.LE.IFATM) THEN
            UMBFC = 5.0D+02*TOANGS*TOANGS*TOHART
            IMD1  = ION
            IMD2  = IMDCATM(ICHOSE)
            X     = CORD(1,IMD1) - CORD(1,IMD2)
            Y     = CORD(2,IMD1) - CORD(2,IMD2)
            Z     = CORD(3,IMD1) - CORD(3,IMD2)
            R2    = X*X + Y*Y + Z*Z
            R     = SQRT(R2)
            DUM   = UMBFC*(R-RCH0)/R
            DEX   = DUM*X
            DEY   = DUM*Y
            DEZ   = DUM*Z
            FFGRD2(1,IMD1)=FFGRD2(1,IMD1) + DEX
            FFGRD2(2,IMD1)=FFGRD2(2,IMD1) + DEY
            FFGRD2(3,IMD1)=FFGRD2(3,IMD1) + DEZ
            FFGRD2(1,IMD2)=FFGRD2(1,IMD2) - DEX
            FFGRD2(2,IMD2)=FFGRD2(2,IMD2) - DEY
            FFGRD2(3,IMD2)=FFGRD2(3,IMD2) - DEZ
            IF(MASWRK.AND.MOD(ISTEP,MDOPT).EQ.0) THEN
               WRITE(IW,'(A,X,I6,A,I3)')
     *      ' BIAS FORCE FOR ATOMS ',IMD1,' AND ',IMD2
            END IF
         END IF
         END IF
C
         CALL UMBRELLA(CORD,FFGRD2,UMBHIS,UM2HIS)
         IF(IFIXSOL.EQ.0) THEN
            CALL CHGRXN(CORD,FFGRD2,CHARG,XTS,YTS,ZTS,
     *                  CMAT1,POT1,QRXN1,NTS)
            CALL POLRXN(CORD,FFGRD2,CHARG,POL,POLSV,DIP,
     *                  FIELD1,FIELD2,FIELD3,
     *                  XTS,YTS,ZTS,CMAT1,POT1,POT2,QRXN1,QRXN2,NTS,
     *                  NONLS1,L1213J)
         END IF
         IF(IFIXSOL.EQ.1) THEN
            CALL FIXSOL(CORD,FFGRD2,CHARG,ZANF,AFIX,QFIX,
     *                  VFIX1,VFIX2,XTSFIX,YTSFIX,ZTSFIX,
     *                  RFIX,IDATOM,DAI,IDDAI,
     *                  POL,POLSV,DIP,FIELD1,FIELD2,FIELD3,
     *                  NONLS1,L1213J)
         END IF
         IF(GOPARR) THEN
            CALL VCLR(ENALL,1,100)
            ENALL( 1) = EN12
            ENALL( 2) = EN123
            ENALL( 3) = EN123R4
            ENALL( 4) = EN123B4
            ENALL( 5) = EN234W1
            ENALL( 6) = ENCHAR
            ENALL( 7) = ENLJR
            ENALL( 8) = ENLJD
            ENALL( 9) = VIR(1)
            ENALL(10) = VIR(2)
            ENALL(11) = VIR(3)
C           ENALL(12) = ENUCCH
            ENALL(13) = ENRXN
            ENALL(14) = ENRXNR
C           ENALL(15) = ENCENT
            ENALL(16) = ENCMAP
            ENALL(18) = SOL1CH
            ENALL(19) = SOL1LJ
            ENALL(21) = SOL2CH
            ENALL(22) = SOL2LJ
            ENALL(33) = PMF1BD
            ENALL(34) = PMF1AG
            ENALL(35) = PMF1DR
            ENALL(36) = PMF1DB
            ENALL(37) = PMF1WG
            ENALL(38) = PMF1CM
            ENALL(39) = PMF1CH
            ENALL(40) = PMF1LJ
            ENALL(41) = ENBIAS
            ENALL(42) = EN12312
            CALL DDI_GSUMF(2410,ENALL  ,42)
            CALL DDI_GSUMF(2411,FFGRD2,3*NFFAT)
            EN12      = ENALL( 1)
            EN123     = ENALL( 2)
            EN123R4   = ENALL( 3)
            EN123B4   = ENALL( 4)
            EN234W1   = ENALL( 5)
            ENCHAR    = ENALL( 6)
            ENLJR     = ENALL( 7)
            ENLJD     = ENALL( 8)
            VIR(1)    = ENALL( 9)
            VIR(2)    = ENALL(10)
            VIR(3)    = ENALL(11)
C           ENUCCH    = ENALL(12)
            ENRXN     = ENALL(13)
            ENRXNR    = ENALL(14)
C           ENCENT    = ENALL(15)
            ENCMAP    = ENALL(16)
            SOL1CH    = ENALL(18)
            SOL1LJ    = ENALL(19)
            SOL2CH    = ENALL(21)
            SOL2LJ    = ENALL(22)
            SOL2IM    = ENALL(23)
            PMF1BD    = ENALL(33)
            PMF1AG    = ENALL(34)
            PMF1DR    = ENALL(35)
            PMF1DB    = ENALL(36)
            PMF1WG    = ENALL(37)
            PMF1CM    = ENALL(38)
            PMF1CH    = ENALL(39)
            PMF1LJ    = ENALL(40)
            ENBIAS    = ENALL(41)
            EN12312   = ENALL(42)
         END IF
C
         IF(IFEPTYP.GT.0.AND.(IDOPOL.GT.0.OR.IEWALD.GT.0)) THEN
            CALL SAVEABPROP(1)
            CALL SAVEFFDATA
            CALL SETFFDATAB(2)
            CALL VCLR(FFGRD2,1,3*NFFAT)
            CALL VCLR(VIR,1,3)
            CALL E00012(CORD,FFGRD2,BOND0,FCBOND,IPAIR,CORDB,
     *                  LSBONDPMA,LSBONDPMB)
            CALL E00123(CORD,FFGRD2,ANGL0,FCANGL,KLIST,CORDB,
     *                  LSANGLPMA,LSANGLPMB)
            CALL E12312(CORD,FFGRD2,ANGL0,KLIST,BOND0,FCSTBD,
     *                  KBLST,CORDB,LSANGLPMA,LSANGLPMB)
            CALL E123B4(CORD,FFGRD2,DIHB0,FCDIHB,NLIST,CORDB,
     *                  LSDIHBPMA,LSDIHBPMB)
            CALL E234W1(CORD,FFGRD2,      FCWAGG,MLIST,CORDB,
     *                  LSWAGGPMA,LSWAGGPMB)
            CALL E123R4(CORD,FFGRD2,VROT,GAMA,NNN,LLIST,CORDB,
     *                  LSDIHRPMA,LSDIHRPMB,FCDIHR)
            CALL ECMAP (CORD,FFGRD2,MAPLST,CMAPCO,CORDB,
     *                  LSCMAPPMA,LSCMAPPMB)
            CALL ELJ126(CORD,FFGRD2,SIG,EPS,SIG2,EPS2,L14J,NONLS1,
     *                  L1213J,SIGB,EPSB,SIG2B,EPS2B,
     *                  NONLSA,NONLSB,L1213A,L1213B,L14A,L14B,
     *                  NONLSPMA,L1213PMA,L14PMA,CORDB,
     *                  NONLSPMB,L1213PMB,L14PMB,FCLJTP,NTYPE)
            CALL ESPHER(CORD,FFGRD2)
            CALL ECHARG(CORD,FFGRD2,CHARG,NONLS1,L1213J,L14J,
     *                  CHARGB,NONLSA,NONLSB,L1213A,L1213B,L14A,L14B,
     *                  NONLSPMA,L1213PMA,L14PMA,CORDB,
     *                  NONLSPMB,L1213PMB,L14PMB)
            CALL UMBRELLA(CORD,FFGRD2,UMBHIS,UM2HIS)
            IF(IFIXSOL.EQ.0) THEN
               CALL CHGRXN(CORD,FFGRD2,CHARG,XTS,YTS,ZTS,
     *                     CMAT1,POT1,QRXN1,NTS)
               CALL POLRXN(CORD,FFGRD2,CHARG,POL,POLSV,DIP,
     *                     FIELD1,FIELD2,FIELD3,
     *                     XTS,YTS,ZTS,CMAT1,POT1,POT2,QRXN1,QRXN2,NTS,
     *                     NONLS1,L1213J)
            END IF
            IF(IFIXSOL.EQ.1) THEN
               CALL FIXSOL(CORD,FFGRD2,CHARG,ZANF,AFIX,QFIX,
     *                     VFIX1,VFIX2,XTSFIX,YTSFIX,ZTSFIX,
     *                     RFIX,IDATOM,DAI,IDDAI,
     *                     POL,POLSV,DIP,FIELD1,FIELD2,FIELD3,
     *                     NONLS1,L1213J)
            END IF
            IF(GOPARR) THEN
               CALL VCLR(ENALL,1,100)
               ENALL( 1) = EN12
               ENALL( 2) = EN123
               ENALL( 3) = EN123R4
               ENALL( 4) = EN123B4
               ENALL( 5) = EN234W1
               ENALL( 6) = ENCHAR
               ENALL( 7) = ENLJR
               ENALL( 8) = ENLJD
               ENALL( 9) = VIR(1)
               ENALL(10) = VIR(2)
               ENALL(11) = VIR(3)
C              ENALL(12) = ENUCCH
               ENALL(13) = ENRXN
               ENALL(14) = ENRXNR
C              ENALL(15) = ENCENT
               ENALL(16) = ENCMAP
               ENALL(18) = SOL1CH
               ENALL(19) = SOL1LJ
               ENALL(21) = SOL2CH
               ENALL(22) = SOL2LJ
               ENALL(33) = PMF1BD
               ENALL(34) = PMF1AG
               ENALL(35) = PMF1DR
               ENALL(36) = PMF1DB
               ENALL(37) = PMF1WG
               ENALL(38) = PMF1CM
               ENALL(39) = PMF1CH
               ENALL(40) = PMF1LJ
               ENALL(41) = ENBIAS
               ENALL(42) = EN12312
               CALL DDI_GSUMF(2410,ENALL  ,42)
               CALL DDI_GSUMF(2411,FFGRD2,3*NFFAT)
               EN12      = ENALL( 1)
               EN123     = ENALL( 2)
               EN123R4   = ENALL( 3)
               EN123B4   = ENALL( 4)
               EN234W1   = ENALL( 5)
               ENCHAR    = ENALL( 6)
               ENLJR     = ENALL( 7)
               ENLJD     = ENALL( 8)
               VIR(1)    = ENALL( 9)
               VIR(2)    = ENALL(10)
               VIR(3)    = ENALL(11)
C              ENUCCH    = ENALL(12)
               ENRXN     = ENALL(13)
               ENRXNR    = ENALL(14)
C              ENCENT    = ENALL(15)
               ENCMAP    = ENALL(16)
               SOL1CH    = ENALL(18)
               SOL1LJ    = ENALL(19)
               SOL2CH    = ENALL(21)
               SOL2LJ    = ENALL(22)
               PMF1BD    = ENALL(33)
               PMF1AG    = ENALL(34)
               PMF1DR    = ENALL(35)
               PMF1DB    = ENALL(36)
               PMF1WG    = ENALL(37)
               PMF1CM    = ENALL(38)
               PMF1CH    = ENALL(39)
               PMF1LJ    = ENALL(40)
               ENBIAS    = ENALL(41)
               EN12312   = ENALL(42)
            END IF
            CALL SAVEABPROP(2)
            CALL SETFFDATAB(1)
            CALL MIXABPROP
         END IF
C
         IF(MDSTEP.EQ.0) THEN
            CALL DCOPY(3*NFFAT,FFGRD2,1,FFGRD1,1)
            GOTO 51
         END IF
C
C        - FIXED ATOMS DO NOT MOVE, SO MAKE THE MASS INFINITE
         IF(NACTMM.GT.0) THEN
            DO IFFAT = 1, NFFAT
               ONEMAS(IFFAT) = 0.0D+00
            ENDDO
            DO KOPT = 1, NACTMM
               IFFAT = LACTMM(KOPT)
               IF(IFFAT.GT.0) ONEMAS(IFFAT) = 1.0D+00/ZMAS(IFFAT)
            ENDDO
         END IF
         DO KFIX = 1, NFIXMM
            KFFAT = IFIXMM(KFIX)
            IF(KFFAT.GT.0) ONEMAS(KFFAT) = 0.0D+00
         ENDDO
C
C        -- CALCULATE V(T+DT)
         IF(INTALG.EQ.1) THEN
            IF(MDSTEP.EQ.1)THEN
               DO IFFAT = 1, NFFAT
                  DUM  = ONEMAS(IFFAT)*DT
                  DO III = 1, 3
                     VEL(III,IFFAT)=VEL(III,IFFAT)
     *                            -PT5*(FFGRD1(III,IFFAT)
     *                                 +FFGRD2(III,IFFAT))*DUM
                     FFGRD1(III,IFFAT)=FFGRD2(III,IFFAT)
                  ENDDO
               ENDDO
            ELSE
               DO IFFAT = 1, NFFAT
                  DUM  = ONEMAS(IFFAT)*DT
                  DO III = 1, 3
                     VEL(III,IFFAT)=VEL(III,IFFAT)
     *                            -(ONETHIRD*FFGRD2(III,IFFAT)
     *                              +FIVESIX*FFGRD1(III,IFFAT)
     *                               -ONESIX*FFGRD0(III,IFFAT))*DUM
                     FFGRD0(III,IFFAT)=FFGRD1(III,IFFAT)
                     FFGRD1(III,IFFAT)=FFGRD2(III,IFFAT)
                  ENDDO
               ENDDO
            END IF
         END IF
         IF(INTALG.EQ.2) THEN
            DO IFFAT = 1, NFFAT
               DUM  = ONEMAS(IFFAT)*DT
               DO III = 1, 3
                  VEL(III,IFFAT)=VEL(III,IFFAT)
     *                         -PT5*(FFGRD1(III,IFFAT)
     *                              +FFGRD2(III,IFFAT))*DUM
                  FFGRD1(III,IFFAT)=FFGRD2(III,IFFAT)
               ENDDO
            ENDDO
         END IF
C
         IF(NRATTLE.GT.0) THEN
            CALL RATTLE2A(CORD,VEL,DSTRAT,ONEMAS,LSTRAT,MDSTEP)
         END IF
C
         IF(NRATM1.GT.0) THEN
            CALL IRMDF2(CORD,VEL,ONEMAS,MDSTEP)
         END IF
C
         DO KFIX = 1, NFIXMM
            KFFAT = IFIXMM(KFIX)
            IF(KFFAT.GT.0) ONEMAS(KFFAT) = 1.0D+00/ZMAS(KFFAT)
         ENDDO
         IF(NACTMM.GT.0) THEN
            DO IFFAT = 1, NFFAT
               ONEMAS(IFFAT) = 1.0D+00/ZMAS(IFFAT)
            ENDDO
         END IF
C
C        -- SYNCHRONIZE CORD AND VEL EVERY 200 STEPS --
         IF(GOPARR.AND.MOD(ISTEP,200).EQ.0) THEN
            CALL DDI_BCAST(462,'F',CORD,3*NFFAT,MASTER)
            CALL DDI_BCAST(464,'F',VEL,3*NFFAT,MASTER)
         END IF
C
  51     CONTINUE
C
C        -- CALCULATE PROPERTIES
         CALL MDPROP(CORD,VEL,QMVEL,ZMAS,QMZMAS,MDSTEP,
     *               LISTQM,TEMP,PRES,UMBHIS,UM2HIS,
     *               PRESX,PRESY,PRESZ)
C
C
         IF(NDFS.EQ.MAGIC2)THEN
C
C        -- RELOCATE BUFFER GAS MOLECULES TO MASTER BOX --
C           MUST BE DONE AFTER MDPROP AND BEFORE PRINT OUT
C           (EVERY 100 MD STEPS)
C
         IF(MIN(XBOX,YBOX,ZBOX).LT.1.0D+30.AND.
     *     (MOD(ISTEP,100).EQ.0 .OR.
     *      MOD(ISTEP,KOUT).EQ.0 .OR. ISTEP.EQ.NSTEP)) THEN
C           - COM OF NATPDB -
            XCOM = 0.0D+00
            YCOM = 0.0D+00
            ZCOM = 0.0D+00
            TMAS = 0.0D+00
            DO IFFAT = 1,NATPDB
               XCOM = XCOM + CORD(1,IFFAT)*ZMAS(IFFAT)
               YCOM = YCOM + CORD(2,IFFAT)*ZMAS(IFFAT)
               ZCOM = ZCOM + CORD(3,IFFAT)*ZMAS(IFFAT)
               TMAS = TMAS + ZMAS(IFFAT)
            ENDDO
            ONETMAS = 1.0D+00/TMAS
            XCOM    = XCOM*ONETMAS
            YCOM    = YCOM*ONETMAS
            ZCOM    = ZCOM*ONETMAS
            DO IFFAT=NATPDB+1,NFFAT
               NUCZ  = INT(ZANF(IFFAT)+0.01D+00)
               Q     = CHARG(IFFAT)
C              -- NON WATER BUFFER GAS --
               IF(NUCZ.NE.8) THEN
               JYES1 = 0
               CX    = CORD(1,IFFAT) - XCOM
               CY    = CORD(2,IFFAT) - YCOM
               CZ    = CORD(3,IFFAT) - ZCOM
               PBCX  = XBOX * ANINT(CX*ONEXBOX)
               PBCY  = YBOX * ANINT(CY*ONEYBOX)
               PBCZ  = ZBOX * ANINT(CZ*ONEZBOX)
               JFFAT=IFFAT+1
               CXH = CORD(1,IFFAT)-CORD(1,JFFAT)
               CYH = CORD(2,IFFAT)-CORD(2,JFFAT)
               CZH = CORD(3,IFFAT)-CORD(3,JFFAT)
               R2H = CXH*CXH+CYH*CYH+CZH*CZH
C              --  2.0 A  --
               IF(R2H.LE.14.28D+00) JYES1 = JFFAT
               JFFAT=IFFAT-1
               CXH = CORD(1,IFFAT)-CORD(1,JFFAT)
               CYH = CORD(2,IFFAT)-CORD(2,JFFAT)
               CZH = CORD(3,IFFAT)-CORD(3,JFFAT)
               R2H = CXH*CXH+CYH*CYH+CZH*CZH
C              --  2.0 A  --
               IF(R2H.LE.14.28D+00) JYES1 = -1
               IF(JYES1.LT.0) THEN
               ELSE IF(JYES1.EQ.0) THEN
                  CORD(1,IFFAT) = CORD(1,IFFAT) - PBCX
                  CORD(2,IFFAT) = CORD(2,IFFAT) - PBCY
                  CORD(3,IFFAT) = CORD(3,IFFAT) - PBCZ
                  CORDSV(1,IFFAT) = CORDSV(1,IFFAT) - PBCX
                  CORDSV(2,IFFAT) = CORDSV(2,IFFAT) - PBCY
                  CORDSV(3,IFFAT) = CORDSV(3,IFFAT) - PBCZ
                  CORDSV2(1,IFFAT) = CORDSV2(1,IFFAT) - PBCX
                  CORDSV2(2,IFFAT) = CORDSV2(2,IFFAT) - PBCY
                  CORDSV2(3,IFFAT) = CORDSV2(3,IFFAT) - PBCZ
               ELSE IF(JYES1.GT.0) THEN
                  CORD(1,IFFAT) = CORD(1,IFFAT) - PBCX
                  CORD(2,IFFAT) = CORD(2,IFFAT) - PBCY
                  CORD(3,IFFAT) = CORD(3,IFFAT) - PBCZ
                  CORDSV(1,IFFAT) = CORDSV(1,IFFAT) - PBCX
                  CORDSV(2,IFFAT) = CORDSV(2,IFFAT) - PBCY
                  CORDSV(3,IFFAT) = CORDSV(3,IFFAT) - PBCZ
                  CORDSV2(1,IFFAT) = CORDSV2(1,IFFAT) - PBCX
                  CORDSV2(2,IFFAT) = CORDSV2(2,IFFAT) - PBCY
                  CORDSV2(3,IFFAT) = CORDSV2(3,IFFAT) - PBCZ
                  CORD(1,JYES1) = CORD(1,JYES1) - PBCX
                  CORD(2,JYES1) = CORD(2,JYES1) - PBCY
                  CORD(3,JYES1) = CORD(3,JYES1) - PBCZ
                  CORDSV(1,JYES1) = CORDSV(1,JYES1) - PBCX
                  CORDSV(2,JYES1) = CORDSV(2,JYES1) - PBCY
                  CORDSV(3,JYES1) = CORDSV(3,JYES1) - PBCZ
                  CORDSV2(1,JYES1) = CORDSV2(1,JYES1) - PBCX
                  CORDSV2(2,JYES1) = CORDSV2(2,JYES1) - PBCY
                  CORDSV2(3,JYES1) = CORDSV2(3,JYES1) - PBCZ
               END IF
               END IF
C
C              -- WATER --
               IF(NUCZ.EQ.8) THEN
                  JYES1 = 0
                  JYES2 = 0
                  JYES3 = 0
                  JYES4 = 0
                  CX    = CORD(1,IFFAT) - CENTX
                  CY    = CORD(2,IFFAT) - CENTY
                  CZ    = CORD(3,IFFAT) - CENTZ
                  PBCX  = XBOX * ANINT(CX*ONEXBOX)
                  PBCY  = YBOX * ANINT(CY*ONEYBOX)
                  PBCZ  = ZBOX * ANINT(CZ*ONEZBOX)
                  IF((ABS(PBCX)+ABS(PBCY)+ABS(PBCZ)).GT.0.0D+00)THEN
                     DO JFFAT=IFFAT+1,IFFAT+4
                        IF(INT(ZANF(JFFAT)+0.01D+00).EQ.1) THEN
                           CXH = CORD(1,IFFAT)-CORD(1,JFFAT)
                           CYH = CORD(2,IFFAT)-CORD(2,JFFAT)
                           CZH = CORD(3,IFFAT)-CORD(3,JFFAT)
                           R2H = CXH*CXH+CYH*CYH+CZH*CZH
C                          --  1.4 A  --
                           IF(R2H.LE.7.0D+00) THEN
                             IF(JYES3.GT.0.AND.JYES4.EQ.0) JYES4=JFFAT
                             IF(JYES2.GT.0.AND.JYES3.EQ.0) JYES3=JFFAT
                             IF(JYES1.GT.0.AND.JYES2.EQ.0) JYES2=JFFAT
                             IF(JYES1.EQ.0)                JYES1=JFFAT
                           END IF
                        END IF
                     ENDDO
                     IF(JYES1.GT.0.AND.JYES2.GT.0) THEN
                        CORD(1,IFFAT) = CORD(1,IFFAT) - PBCX
                        CORD(2,IFFAT) = CORD(2,IFFAT) - PBCY
                        CORD(3,IFFAT) = CORD(3,IFFAT) - PBCZ
                        CORD(1,JYES1) = CORD(1,JYES1) - PBCX
                        CORD(2,JYES1) = CORD(2,JYES1) - PBCY
                        CORD(3,JYES1) = CORD(3,JYES1) - PBCZ
                        CORD(1,JYES2) = CORD(1,JYES2) - PBCX
                        CORD(2,JYES2) = CORD(2,JYES2) - PBCY
                        CORD(3,JYES2) = CORD(3,JYES2) - PBCZ
                     END IF
                     IF(JYES3.GT.0.AND.JYES4.GT.0) THEN
                        CORD(1,JYES3) = CORD(1,JYES3) - PBCX
                        CORD(2,JYES3) = CORD(2,JYES3) - PBCY
                        CORD(3,JYES3) = CORD(3,JYES3) - PBCZ
                        CORD(1,JYES4) = CORD(1,JYES4) - PBCX
                        CORD(2,JYES4) = CORD(2,JYES4) - PBCY
                        CORD(3,JYES4) = CORD(3,JYES4) - PBCZ
                     END IF
                  END IF
               END IF
C
            ENDDO
         END IF
C
         ELSE
C
C        -- RELOCATE WATER MOLECULES AND IONS TO THE MASTER BOX --
C           MUST BE DONE AFTER MDPROP AND BEFORE PRINT OUT
C
         IF(KMASTER.GT.0.AND.MIN(XBOX,YBOX,ZBOX).LT.1.0D+30.AND.
     *      MOD(ISTEP,KMASTER).EQ.0) THEN
         DO IFFAT=1,NFFAT
         JYES1 = 0
         JYES2 = 0
         JYES3 = 0
         JYES4 = 0
         NUCZ  = INT(ZANF(IFFAT)+0.01D+00)
         Q     = CHARG(IFFAT)
C        - IONS AND MONOATOMIC GAS -
         IF((NUCZ.EQ. 2.AND.Q.EQ. 0.0D+00).OR.
     *      (NUCZ.EQ. 3.AND.Q.EQ.+1.0D+00).OR.
     *      (NUCZ.EQ. 9.AND.Q.EQ.-1.0D+00).OR.
     *      (NUCZ.EQ.10.AND.Q.EQ. 0.0D+00).OR.
     *      (NUCZ.EQ.11.AND.Q.EQ.+1.0D+00).OR.
     *      (NUCZ.EQ.12.AND.Q.EQ.+2.0D+00).OR.
     *      (NUCZ.EQ.13.AND.Q.EQ.+3.0D+00).OR.
     *      (NUCZ.EQ.17.AND.Q.EQ.-1.0D+00).OR.
     *      (NUCZ.EQ.18.AND.Q.EQ. 0.0D+00).OR.
     *      (NUCZ.EQ.19.AND.Q.EQ.+1.0D+00).OR.
     *      (NUCZ.EQ.20.AND.Q.EQ.+2.0D+00).OR.
     *      (NUCZ.GE.21.AND.NUCZ.LE.30.AND.Q.EQ.+1.0D+00).OR.
     *      (NUCZ.GE.21.AND.NUCZ.LE.30.AND.Q.EQ.+2.0D+00).OR.
     *      (NUCZ.GE.21.AND.NUCZ.LE.30.AND.Q.EQ.+3.0D+00).OR.
     *      (NUCZ.GE.21.AND.NUCZ.LE.30.AND.Q.EQ.+4.0D+00).OR.
     *      (NUCZ.EQ.35.AND.Q.EQ.-1.0D+00).OR.
     *      (NUCZ.EQ.36.AND.Q.EQ. 0.0D+00).OR.
     *      (NUCZ.EQ.37.AND.Q.EQ.+1.0D+00).OR.
     *      (NUCZ.EQ.53.AND.Q.EQ.-1.0D+00).OR.
     *      (NUCZ.EQ.54.AND.Q.EQ. 0.0D+00).OR.
     *      (NUCZ.EQ.55.AND.Q.EQ.+1.0D+00)     ) THEN
            CX    = CORD(1,IFFAT) - CENTX
            CY    = CORD(2,IFFAT) - CENTY
            CZ    = CORD(3,IFFAT) - CENTZ
            PBCX  = XBOX * ANINT(CX*ONEXBOX)
            PBCY  = YBOX * ANINT(CY*ONEYBOX)
            PBCZ  = ZBOX * ANINT(CZ*ONEZBOX)
            CORD(1,IFFAT) = CORD(1,IFFAT) - PBCX
            CORD(2,IFFAT) = CORD(2,IFFAT) - PBCY
            CORD(3,IFFAT) = CORD(3,IFFAT) - PBCZ
C           - ONLY DIFFUSION NEEDS THIS CORRECTION -
            IF(NDFS.GT.0) THEN
               DFSC0(1,IFFAT) = DFSC0(1,IFFAT) - PBCX
               DFSC0(2,IFFAT) = DFSC0(2,IFFAT) - PBCY
               DFSC0(3,IFFAT) = DFSC0(3,IFFAT) - PBCZ
            END IF
         END IF
C        - WATER -
         IF(NUCZ.EQ.8) THEN
            CX    = CORD(1,IFFAT) - CENTX
            CY    = CORD(2,IFFAT) - CENTY
            CZ    = CORD(3,IFFAT) - CENTZ
            PBCX  = XBOX * ANINT(CX*ONEXBOX)
            PBCY  = YBOX * ANINT(CY*ONEYBOX)
            PBCZ  = ZBOX * ANINT(CZ*ONEZBOX)
            IF((ABS(PBCX)+ABS(PBCY)+ABS(PBCZ)).GT.0.0D+00)THEN
               NHWAT = 0
               DO 78 JFFAT=IFFAT+1,IFFAT+4
                  IF(JFFAT.LE.NFFAT) THEN
                  IF(INT(ZANF(JFFAT)+0.01D+00).EQ.1) NHWAT = NHWAT + 1
                  IF(INT(ZANF(JFFAT)+0.01D+00).EQ.8) GOTO 78
                  END IF
 78            CONTINUE
               DO JFFAT=IFFAT+1,IFFAT+NHWAT
                  IF(INT(ZANF(JFFAT)+0.01D+00).EQ.1) THEN
                     CXH = CORD(1,IFFAT)-CORD(1,JFFAT)
                     CYH = CORD(2,IFFAT)-CORD(2,JFFAT)
                     CZH = CORD(3,IFFAT)-CORD(3,JFFAT)
                     R2H = CXH*CXH+CYH*CYH+CZH*CZH
C                    --  1.4 A  --
                     IF(R2H.LE.7.0D+00) THEN
                        IF(JYES3.GT.0.AND.JYES4.EQ.0) JYES4 = JFFAT
                        IF(JYES2.GT.0.AND.JYES3.EQ.0) JYES3 = JFFAT
                        IF(JYES1.GT.0.AND.JYES2.EQ.0) JYES2 = JFFAT
                        IF(JYES1.EQ.0)                JYES1 = JFFAT
                     END IF
                  END IF
               ENDDO
               IF(JYES1.GT.0.AND.JYES2.GT.0) THEN
                  CORD(1,IFFAT) = CORD(1,IFFAT) - PBCX
                  CORD(2,IFFAT) = CORD(2,IFFAT) - PBCY
                  CORD(3,IFFAT) = CORD(3,IFFAT) - PBCZ
                  CORD(1,JYES1) = CORD(1,JYES1) - PBCX
                  CORD(2,JYES1) = CORD(2,JYES1) - PBCY
                  CORD(3,JYES1) = CORD(3,JYES1) - PBCZ
                  CORD(1,JYES2) = CORD(1,JYES2) - PBCX
                  CORD(2,JYES2) = CORD(2,JYES2) - PBCY
                  CORD(3,JYES2) = CORD(3,JYES2) - PBCZ
                  IF(NDFS.GT.0) THEN
                     DFSC0(1,IFFAT) = DFSC0(1,IFFAT) - PBCX
                     DFSC0(2,IFFAT) = DFSC0(2,IFFAT) - PBCY
                     DFSC0(3,IFFAT) = DFSC0(3,IFFAT) - PBCZ
                     DFSC0(1,JYES1) = DFSC0(1,JYES1) - PBCX
                     DFSC0(2,JYES1) = DFSC0(2,JYES1) - PBCY
                     DFSC0(3,JYES1) = DFSC0(3,JYES1) - PBCZ
                     DFSC0(1,JYES2) = DFSC0(1,JYES2) - PBCX
                     DFSC0(2,JYES2) = DFSC0(2,JYES2) - PBCY
                     DFSC0(3,JYES2) = DFSC0(3,JYES2) - PBCZ
                  END IF
               END IF
               IF(JYES3.GT.0.AND.JYES4.GT.0) THEN
                  CORD(1,JYES3) = CORD(1,JYES3) - PBCX
                  CORD(2,JYES3) = CORD(2,JYES3) - PBCY
                  CORD(3,JYES3) = CORD(3,JYES3) - PBCZ
                  CORD(1,JYES4) = CORD(1,JYES4) - PBCX
                  CORD(2,JYES4) = CORD(2,JYES4) - PBCY
                  CORD(3,JYES4) = CORD(3,JYES4) - PBCZ
                  IF(NDFS.GT.0) THEN
                     DFSC0(1,JYES3) = DFSC0(1,JYES3) - PBCX
                     DFSC0(2,JYES3) = DFSC0(2,JYES3) - PBCY
                     DFSC0(3,JYES3) = DFSC0(3,JYES3) - PBCZ
                     DFSC0(1,JYES4) = DFSC0(1,JYES4) - PBCX
                     DFSC0(2,JYES4) = DFSC0(2,JYES4) - PBCY
                     DFSC0(3,JYES4) = DFSC0(3,JYES4) - PBCZ
                  END IF
               END IF
            END IF
         END IF
C        - N2 AND N2 GAS -
         IF((NUCZ.EQ.7.OR.NUCZ.EQ.8).AND.Q.EQ.0.0D+00) THEN
            JYES1 = 0
            CX    = CORD(1,IFFAT) - CENTX
            CY    = CORD(2,IFFAT) - CENTY
            CZ    = CORD(3,IFFAT) - CENTZ
            PBCX  = XBOX * ANINT(CX*ONEXBOX)
            PBCY  = YBOX * ANINT(CY*ONEYBOX)
            PBCZ  = ZBOX * ANINT(CZ*ONEZBOX)
            JFFAT=IFFAT+1
            IF(JFFAT.LE.NFFAT) THEN
            NUCZ1 = INT(ZANF(JFFAT)+0.01D+00)
            Q1    = CHARG(JFFAT)
            IF((NUCZ1.EQ.NUCZ).AND.Q1.EQ.0.0D+00) THEN
            CXH = CORD(1,IFFAT)-CORD(1,JFFAT)
            CYH = CORD(2,IFFAT)-CORD(2,JFFAT)
            CZH = CORD(3,IFFAT)-CORD(3,JFFAT)
            R2H = CXH*CXH+CYH*CYH+CZH*CZH
            IF(R2H.LE.14.28D+00) JYES1 = JFFAT   !  2.0 A
            END IF
            END IF
            JFFAT=IFFAT-1
            NUCZ1 = INT(ZANF(JFFAT)+0.01D+00)
            Q1    = CHARG(JFFAT)
            IF((NUCZ1.EQ.NUCZ).AND.Q1.EQ.0.0D+00) THEN
            CXH = CORD(1,IFFAT)-CORD(1,JFFAT)
            CYH = CORD(2,IFFAT)-CORD(2,JFFAT)
            CZH = CORD(3,IFFAT)-CORD(3,JFFAT)
            R2H = CXH*CXH+CYH*CYH+CZH*CZH
            IF(R2H.LE.14.28D+00) JYES1 = 0       !  2.0 A
            END IF
            IF(JYES1.GT.0) THEN
               CORD(1,IFFAT) = CORD(1,IFFAT) - PBCX
               CORD(2,IFFAT) = CORD(2,IFFAT) - PBCY
               CORD(3,IFFAT) = CORD(3,IFFAT) - PBCZ
               CORD(1,JYES1) = CORD(1,JYES1) - PBCX
               CORD(2,JYES1) = CORD(2,JYES1) - PBCY
               CORD(3,JYES1) = CORD(3,JYES1) - PBCZ
            END IF
         END IF
C
         ENDDO
         END IF
C
         END IF
C
C
         CALL TMDATE(TIMSTR)
         IF(MOD(ISTEP,KOUT).EQ.0) CALL TIMIT(1)
         IF(MOD(ISTEP,KOUT).NE.0 .AND. ISTEP.NE.NSTEP) GOTO 200
C
         IF(MASWRK)THEN
C
C           -- PRINT KOUTACT ATOMS --
C
            IF(KOUTACT(1).GT.0) THEN
            WRITE(IW,'(/A,I10)')'!KOUTACT ATOMS AROUND $FFDATA ATOM ',
     *      KOUTACT(1)
            WRITE(IW,'(1X,A,I10,10X,A,F19.2,A)')
     *      '$FFDATA ! KOUTACT     MD STEP',ISTEP,
     *      'TIME=',ISTEP*DT*2.418884326505D-02,' FS'
            ACTX = CORD(1,KOUTACT(1))
            ACTY = CORD(2,KOUTACT(1))
            ACTZ = CORD(3,KOUTACT(1))
            ACTR = DBLE(KOUTACT(2))*1.0D-08*TOBOHR
            ACTR2= ACTR**2
            DO IFFAT=1,NFFAT
               CX    = CORD(1,IFFAT) - ACTX
               CY    = CORD(2,IFFAT) - ACTY
               CZ    = CORD(3,IFFAT) - ACTZ
               PBCX  = XBOX*ANINT(CX*ONEXBOX)
               PBCY  = YBOX*ANINT(CY*ONEYBOX)
               PBCZ  = ZBOX*ANINT(CZ*ONEZBOX)
               CX    = CX - PBCX
               CY    = CY - PBCY
               CZ    = CZ - PBCZ
               RCXYZ = CX**2 + CY**2 + CZ**2
               IF(RCXYZ.LE.ACTR2) THEN
                  WRITE(IW,1000)ATMNAM(IFFAT),ZANF(IFFAT),
     *            (CX+ACTX)*TOANGS,(CY+ACTY)*TOANGS,(CZ+ACTZ)*TOANGS
               END IF
            ENDDO
            WRITE(IW,*)'$END ! KOUTACT'
            END IF
C
C           -- PRINT OUT RESTART COORDINATES + VEL --
C
            WRITE(36,'(1X,A,I10,1X,A,F15.2,A)')
     *      'RESTART COORD + VEL FOR QUANPOL AT MD STEP',ISTEP,
     *      'TIME=',ISTEP*DT*2.418884326505D-02,' FS'
            IF(MIN(XBOX,YBOX,ZBOX).LT.1.0D+30) THEN
              WRITE(36,'(A,3(A,F15.10,1X),/9X,3(A,F15.10,1X),A)')
     *        ' $QUANPO ',
     *        'CENTX=',CENTX*TOANGS,
     *        'CENTY=',CENTY*TOANGS,
     *        'CENTZ=',CENTZ*TOANGS,
     *        'XBOX=',XBOX*TOANGS,
     *        'YBOX=',YBOX*TOANGS,
     *        'ZBOX=',ZBOX*TOANGS,
     *        '$END'
            ELSE IF(SPHRAD.LT.1.0D+30) THEN
              WRITE(36,'(A,3(A,F15.10,1X),/9X,2(A,F15.10,1X),A)')
     *        ' $QUANPO ',
     *        'CENTX=',CENTX*TOANGS,
     *        'CENTY=',CENTY*TOANGS,
     *        'CENTZ=',CENTZ*TOANGS,
     *        'SPHRAD=',SPHRAD*TOANGS,
     *        'RSPHSOL=',RSPHSOL*TOANGS,
     *        '$END'
            ELSE
              WRITE(36,'(A,3(A,F15.10,1X),A)')
     *        ' $QUANPO ',
     *        'CENTX=',CENTX*TOANGS,
     *        'CENTY=',CENTY*TOANGS,
     *        'CENTZ=',CENTZ*TOANGS,
     *        '$END'
            END IF
C
            WRITE(36,'(1X,A,I10,10X,A,F19.2,A)')
     *      '$FFDATA          !    MD STEP',ISTEP,
     *      'TIME=',ISTEP*DT*2.418884326505D-02,' FS'
            IF(IFEPTYP.EQ.1) THEN
               WRITE(36,'(1X,A,I10,10X,A,F19.2,A)')
     *         '$FFDATB          !    MD STEP',ISTEP,
     *         'TIME=',ISTEP*DT*2.418884326505D-02,' FS'
            END IF
            WRITE(36,*)'COORDINATES  NUC                   X',
     *          '                   Y                   Z'
            IF(NDFS.EQ.MAGIC2)THEN
            DO IFFAT = 1, N1FFAT
               WRITE(36,1002)ATMNAM(IFFAT),ZANF(IFFAT),
     *         CORD(1,IFFAT)*TOANGS,
     *         CORD(2,IFFAT)*TOANGS,
     *         CORD(3,IFFAT)*TOANGS
            ENDDO
            ELSE
            DO IFFAT = 1, N1FFAT
               WRITE(36,1000)ATMNAM(IFFAT),ZANF(IFFAT),
     *         CORD(1,IFFAT)*TOANGS,
     *         CORD(2,IFFAT)*TOANGS,
     *         CORD(3,IFFAT)*TOANGS
            ENDDO
            END IF
            WRITE(36,*)'STOP'
            WRITE(36,*)'MMVELOCITY                    VX',
     *          '                    VY                    VZ'
            DO IFFAT = 1, N1FFAT
               WRITE(36,1001) ATMNAM(IFFAT),
     *         VEL(1,IFFAT),VEL(2,IFFAT),VEL(3,IFFAT)
            ENDDO
            WRITE(36,*)'STOP'
C
            IF(IFEPTYP.EQ.2) THEN
               WRITE(36,'(1X,A,I10,10X,A,F19.2,A)')
     *         '$FFDATB          !    MD STEP',ISTEP,
     *         'TIME=',ISTEP*DT*2.418884326505D-02,' FS'
               WRITE(36,*)'COORDINATES  NUC                   X',
     *             '                   Y                   Z'
               DO IFFAT = 1, N2FFAT
                  WRITE(36,1000)ATMNAM(IFFAT),ZANF(IFFAT),
     *            CORDB(1,IFFAT)*TOANGS,
     *            CORDB(2,IFFAT)*TOANGS,
     *            CORDB(3,IFFAT)*TOANGS
               ENDDO
               WRITE(36,*)'STOP'
               WRITE(36,*)'MMVELOCITY                    VX',
     *             '                    VY                    VZ'
               DO IFFAT = 1, N2FFAT
                  WRITE(36,1001) ATMNAM(IFFAT),
     *            VEL(1,IFFAT),VEL(2,IFFAT),VEL(3,IFFAT)
               ENDDO
               WRITE(36,*)'STOP'
            END IF
C
         END IF
C
         CALL FLSHBF(IW)
         CALL FLSHBF(36)
         CALL TIMIT(1)
 200     CONTINUE
C
C        -- SCALE VELOCITY AND VOLUME --
         CALL DCOPY(3*NFFAT,CORD,1,OLDCORD,1)
         CALL DCOPY(3*NFFAT,VEL,1,VELSV,1)
         CALL TPSTAT(CORD,VEL,QMVEL,ZMAS,QMZMAS,ONEMAS,QM1MAS,
     *               MDSTEP,LISTQM,TEMP,PRES,
     *               PRESX,PRESY,PRESZ,NONLSTQ)
         IF(NACTMM.GT.0) THEN
            DO KOPT = 1, NACTMM
               IFFAT = LACTMM(KOPT)
               IF(IFFAT.GT.0) THEN
                  CORD(1,IFFAT)=(CORD(1,IFFAT)+3.1416D+00)*1.0D-40
                  CORD(2,IFFAT)=(CORD(2,IFFAT)+2.7183D+00)*1.0D-40
                  CORD(3,IFFAT)=(CORD(3,IFFAT)+1.4427D+00)*1.0D-40
                  VEL (1,IFFAT)=VEL(1,IFFAT) + 1.0D+00
               END IF
            ENDDO
            DO IFFAT = 1, NFFAT
               IF(ABS(CORD(1,IFFAT))+ABS(CORD(2,IFFAT))
     *           +ABS(CORD(3,IFFAT)).LT.1.0D-20) THEN
                  CORD(1,IFFAT) = CORD(1,IFFAT)*1.0D+40 - 3.1416D+00
                  CORD(2,IFFAT) = CORD(2,IFFAT)*1.0D+40 - 2.7183D+00
                  CORD(3,IFFAT) = CORD(3,IFFAT)*1.0D+40 - 1.4427D+00
               ELSE
                  DO III = 1, 3
                     CORD(III,IFFAT) = OLDCORD(III,IFFAT)
                  ENDDO
               END IF
               IF(VEL(1,IFFAT).GT.0.5D+00) THEN
                  VEL(1,IFFAT) = VEL(1,IFFAT) - 1.0D+00
               ELSE
                  DO III = 1, 3
                     VEL(III,IFFAT) = VELSV(III,IFFAT)
                  ENDDO
               END IF
            ENDDO
         END IF
         DO KFIX=1,NFIXMM
            DO III = 1, 3
               CORD(III,IFIXMM(KFIX)) = OLDCORD(III,IFIXMM(KFIX))
               VEL (III,IFIXMM(KFIX)) = VELSV(III,IFIXMM(KFIX))
            ENDDO
         ENDDO
         IF(MDOPT.GT.0.AND.ISTEP.GT.0.AND.MOD(ISTEP,MDOPT).EQ.0) THEN
            CALL DCOPY(3*NFFAT,CORD,1,CORDG,1)
            CALL DCOPY(3*NFFAT,CORDSV,1,CORDGSV,1)
            CALL DCOPY(3*NFFAT,CORDSV2,1,CORDGSV2,1)
            CALL DCOPY(3*NFFAT,CORDSVQ,1,CORDGSVQ,1)
            CALL DCOPY(3*NFFAT,FFGRD0,1,FFGRDG0,1)
            CALL DCOPY(3*NFFAT,FFGRD1,1,FFGRDG1,1)
            CALL DCOPY(3*NFFAT,FFGRD2,1,FFGRDG2,1)
            NSTEPSV=NSTEP
            NSTEP  =20000
            JOUTSV =JOUT
            JOUT   =20000
            CALL FFOPT1(ATMNAM,CORD,CORDSV,ZANF,CHARG,POL,POLSV,DIP,
     *                  FIELD1,FIELD2,FIELD3,SIG,EPS,SIG2,EPS2,
     *                  BOND0,FCBOND,FCSTBD,ANGL0,FCANGL,FCWAGG,
     *                  DIHB0,FCDIHB,FCDIHR,VROT,NNN,GAMA,IPAIR,
     *                  KLIST,LLIST,KBLST,MLIST,NLIST,L1213J,L14J,
     *                  FFGRD0,FFGRD2,FFGRD1,FCLJTP,NTYPE,
     *                  XTS,YTS,ZTS,CMAT1,POT1,POT2,QRXN1,QRXN2,NTS,
     *                  NONLS1,NONLSTQ,MAPLST,CMAPCO,LSTCELL,NONLS2,
     *                  CORDSV2,CORDSVQ,MVFASTS2,MVFASTS3,MVFASTS4,
     *                  MVFASTL2,MVFASTL3,MVFASTL4,AFIX,QFIX,
     *                  RFIX,IDATOM,DAI,IDDAI,
     *                  VFIX1,VFIX2,XTSFIX,YTSFIX,ZTSFIX,
     *                  NONLSA,NONLSB,L1213A,L1213B,L14A,L14B,
     *                  CHARGB,
     *                  SIGB,EPSB,SIG2B,EPS2B,NONLSPMA,
     *                  L1213PMA,L14PMA,CORDB,LSBONDPMA,LSANGLPMA,
     *                  LSDIHRPMA,LSDIHBPMA,LSWAGGPMA,LSCMAPPMA,
     *                  NONLSPMB,L1213PMB,L14PMB,
     *                  LSBONDPMB,LSANGLPMB,LSDIHRPMB,LSDIHBPMB,
     *                  LSWAGGPMB,LSCMAPPMB,UMBHIS,UM2HIS)
            NSTEP  = NSTEPSV
            JOUT   = JOUTSV
            CALL DCOPY(3*NFFAT,CORDG,1,CORD,1)
            CALL DCOPY(3*NFFAT,CORDGSV,1,CORDSV,1)
            CALL DCOPY(3*NFFAT,CORDGSV2,1,CORDSV2,1)
            CALL DCOPY(3*NFFAT,CORDGSVQ,1,CORDSVQ,1)
            CALL DCOPY(3*NFFAT,FFGRDG0,1,FFGRD0,1)
            CALL DCOPY(3*NFFAT,FFGRDG1,1,FFGRD1,1)
            CALL DCOPY(3*NFFAT,FFGRDG2,1,FFGRD2,1)
         END IF
      ENDDO
C
 1000 FORMAT(1X,A10,1X,F5.1,1X,F19.13,1X,F19.13,1X,F19.13)
 1002 FORMAT(1X,A10,1X,F5.1,1X,F19.9,1X,F19.9,1X,F19.9)
 1001 FORMAT(1X,A10,1X,F21.18,1X,F21.18,1X,F21.18)
C
      CALL SEQCLO(36,'KEEP')
      CALL SEQCLO(38,'KEEP')
C
      RETURN
      END
C*MODULE QUANPOC  *DECK FFMD2
!>
!> @brief    QM/MM MD driver
!>
!> @author   Nandun Thellamurege, Hui Li
!>           - Jan 2011
!>
!> @details  MD simulation for QM/MM atoms
!>
      SUBROUTINE FFMD2(ATMNAM,CORD,CORDSV,ZANF,
     *                 ZMAS,ONEMAS,QMZMAS,QM1MAS,
     *                 CHARG,VEL,QMVEL,
     *                 FFGRD0,FFGRD1,FFGRD2,
     *                 QMGRD0,QMGRD1,QMGRD2,
     *                 LISTQM,NONLS1,NONLSTQ,
     *                 LSTCELL,NONLS2,CORDSV2,CORDSVQ,
     *                 MVFASTS2,MVFASTS3,MVFASTS4,
     *                 MVFASTL2,MVFASTL3,MVFASTL4,
     *                 NONLSA,NONLSB,
     *                 OLDCORD,LSTRAT,DSTRAT,NONLSPMA,NONLSPMB,
     *                 CORDB,UMBHIS,UM2HIS,VELSV,
     *                 OLDC,QMVELSV,DFSC0)
      use mx_limits, only: mxatm,mxao,mxrt
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      LOGICAL TRIPLET,SG1T,TAMMD,TPA,ALPHKWD,BETAKWD
      LOGICAL MREKT,MRDEA
      LOGICAL GOPARR,DSKWRK,MASWRK
      LOGICAL MPTEST
C
      PARAMETER (TOANGS=0.52917724924D+00)
      PARAMETER (TOBOHR=1.0D+00/TOANGS)
      PARAMETER (ZERO=0.0D+00)
      PARAMETER (PT5=0.5D+00)
      PARAMETER (TWOTHIRD=2.0D+00/3.0D+00)
      PARAMETER (ONETHIRD=1.0D+00/3.0D+00)
      PARAMETER (ONESIX=1.0D+00/6.0D+00)
      PARAMETER (FIVESIX=5.0D+00/6.0D+00)
      PARAMETER (MAXML=1024)
C
      CHARACTER*10  ATMNAM
C
      DIMENSION ATMNAM(NFFAT),CORD(3,NFFAT),CORDSV(3,NFFAT),ZANF(NFFAT),
     *          ZMAS(NFFAT),ONEMAS(NFFAT),QMZMAS(NAT),QM1MAS(NAT),
     *          CHARG(NFFAT),VEL(3,NFFAT),QMVEL(3,NAT),
     *          FFGRD0(3,NFFAT),FFGRD1(3,NFFAT),FFGRD2(3,NFFAT),
     *          QMGRD0(3,NAT),QMGRD1(3,NAT),QMGRD2(3,NAT),LISTQM(*),
     *          NONLSTQ(*),CORDB(3,*),OLDCORD(3,*),
     *          DFSC0(3,*)
      DIMENSION TIMSTR(3)
      DIMENSION FIXQM(3,400)
      DIMENSION VELSV(3,*),
     *          OLDC(3,*),QMVELSV(3,*)
C
      COMMON /ENRGMP/ EMP2,EMP3,EMP4,EMP2A
      COMMON /ENRGYS/ ENUCR,EELCT,ETOT,SZ,SZZ,ECORE,ESCF,EERD,E1,E2,
     *                VEN,VEE,EPOT,EKIN,ESTATE(MXRT),STATN,EDFT(3),EDISP
      COMMON /FFDFS / TIMDFS,QDION,AMION,TEFF,NDFS,NATMGAS,
     *                LFFDFSC,
     *                LFFDFSC0,LFFDFSA,LFFDFSN,LFFDFCOM,KDFS,LFFDFSCAV
      COMMON /FFDIMR/ IDIMER,IBREAK(81),N1213JMM,NMOLE,MATOM(MAXML),
     *                MCHARG(MAXML),MMULT(MAXML),MELEC(MAXML),RDIMER
      COMMON /FFFREE/ SOL1CH,SOL1LJ,SOL1IM,SOLFRE1MM,
     *                SOL2CH,SOL2LJ,SOL2IM,SOLFRE2MM,
     *                ASOL1CH,ASOL1LJ,ASOL1IM,ASOL1MM,
     *                ASOL2CH,ASOL2LJ,ASOL2IM,ASOL2MM,
     *                PMF1BD,PMF1AG,PMF1DR,PMF1DB,PMF1WG,PMF1CM,
     *                PMF1CH,PMF1PO,PMF1LJ,PMF1IM,PMFFRE1MM,
     *                EN12A,EN12B,EN123A,EN123B,EN123R4A,EN123R4B,
     *                EN234W1A,EN234W1B,EN123B4A,EN123B4B,
     *                EN12312A,EN12312B,
     *                ENCHARA,ENCHARB,ENPOLA,ENPOLB,
     *                ENRXNA,ENRXNB,ENRXNPOLA,ENRXNPOLB,
     *                ENRXNRA,ENRXNRB,ENLJRA,ENLJRB,ENLJDA,ENLJDB,
     *                ENCMAPA,ENCMAPB,ENFIXSOA,ENFIXSOB,
     *                ENCENTA,ENCENTB,WSIMUL,WPERT1,WPERT2,
     *                IFEPTYP
      COMMON /FFFRE2/ N1FFAT,N1BOND,N1ANGL,N1DIHR,N1DIHB,N1CMAP,N1WAGG,
     *                N2FFAT,N2BOND,N2ANGL,N2DIHR,N2DIHB,N2CMAP,N2WAGG,
     *                LFF2ATMNAM,LFF2CORD,LFF2ZANF,LFF2ZMAS,
     *                LFF2ONEMAS,LFF2CHARG,LFF2POL,
     *                LFF2SIG,LFF2EPS,LFF2SIG2,LFF2EPS2,
     *                LFF2BOND0,LFF2FCBOND,LFF2ANGL0,LFF2FCANGL,
     *                LFF2FCWAGG,LFF2DIHB0,LFF2FCDIHB,
     *                LFF2VROT,LFF2NNN,LFF2GAMA,LFF2IPAIR,
     *                LFF2KLIST,LFF2LLIST,LFF2MLIST,LFF2NLIST,
     *                LFF2VEL,LFF2QMVEL,LFF2CLPR,LFF2ZLPR,
     *                LFF2NLPR,LFF2MAPLST,
     *                LFFLISTB2A,NTODOA,LFFNONLSA,NTODOB,LFFNONLSB,
     *                N1213A,LFFL1213A,N1213B,LFFL1213B,
     *                N14A,LFFL14A,N14B,LFFL14B,
     *                NTODOPMA,LFFNONLSPMA,NBONDPMA,LFFLSBONDPMA,
     *                NANGLPMA,LFFLSANGLPMA,NDIHRPMA,LFFLSDIHRPMA,
     *                NDIHBPMA,LFFLSDIHBPMA,NWAGGPMA,LFFLSWAGGPMA,
     *                NCMAPPMA,LFFLSCMAPPMA,
     *                N1213PMA,LFFL1213PMA,
     *                N14PMA,LFFL14PMA,
     *                NTODOPMB,LFFNONLSPMB,NBONDPMB,LFFLSBONDPMB,
     *                NANGLPMB,LFFLSANGLPMB,NDIHRPMB,LFFLSDIHRPMB,
     *                NDIHBPMB,LFFLSDIHBPMB,NWAGGPMB,LFFLSWAGGPMB,
     *                NCMAPPMB,LFFLSCMAPPMB,
     *                N1213PMB,LFFL1213PMB,
     *                N14PMB,LFFL14PMB
      COMMON /FFMDPA/ DT,DT2,TEMP0,PRES0,POLTOL,VIR(3),PMEAN,VOLAV,
     *                ENPAV,ENKAV,TEMPAV,BERENDT,BERENDP,VELMAX,
     *                PMEANX,PMEANY,PMEANZ,NSTEP,KMASTER,KOUTACT(2),
     *                IHESS,INTALG,ITSTAT,IPSTAT,JOUT,KOUT,LOUT
      COMMON /FFMPT3/ NACTMM,LACTMM(2020),LFFDIPOLD,JUMBUP,
     *                NACTQM,LACTQM(2020),LFFOLDC,LFFQMVELSV,MMHESS,
     *                LFFQMCHG,LFFQMCHGB,ISWAP,R2SWAP,DFTBMM
      COMMON /FFPARM/ NFFAT,NBOND,NANGL,NDIHR,NDIHB,NCMAP,NWAGG,
     *                N1213J,N14J,NLKQMM,IDOCHG,IDOPOL,IDOLJ,IDOCMAP
      COMMON /FFPBSW/ XBOX,YBOX,ZBOX,SWF,SWFDX,SWFDY,SWFDZ,
     *                SWRA,ONESWRA,SWRA2,ONESWRA2,
     *                SWRB,ONESWRB,SWRB2,ONESWRB2,ONESWRB4,
     *                SWFDUM3,SWFDUM4,SWFDUM5,
     *                SWRAQ,ONESWRAQ,SWRAQ2,ONESWRAQ2,
     *                SWRBQ,ONESWRBQ,SWRBQ2,ONESWRBQ2,ONESWRBQ4,
     *                SWFDUM3Q,SWFDUM4Q,SWFDUM5Q,
     *                QMSIZE,QMCX,QMCY,QMCZ,QMCXSV,QMCYSV,QMCZSV,
     *                CENTX,CENTY,CENTZ,BUFWID1,BUFWID2,RDAMP,
     *                EFIELDX,EFIELDY,EFIELDZ,QMCXSV2,QMCYSV2,QMCZSV2,
     *                EPS1RB,EPS1RB3,ONEXBOX,ONEYBOX,ONEZBOX,
     *                LQMCT,MXLIST1,NTODO,NTODOSV,NTODOQ,
     *                ISWITCH,ISHIFT,IPOLSHF,
     *                LFFLSTCELL,LFFCORDSV,
     *                LFFPOLSV,LFFCORDSV2,LFFNONLS2,LFFCORDSVQ,
     *                LFFMVFASTS2,LFFMVFASTS3,LFFMVFASTS4,
     *                LFFMVFASTL2,LFFMVFASTL3,LFFMVFASTL4,
     *                MXCHECK,MXLIST2,NTODO2,NTODO2SV
      COMMON /FFQMFE/ ETOTQA,EMP2QA,ESTATEQA,ETOTQB,EMP2QB,ESTATEQB,
     *                FRE1QMMM,FRE2QMMM,
     *                VIRA(3),VIRB(3),
     *                LFFFFGRDA,LFFFFGRDB,LFFQMGRDA,LFFQMGRDB,
     *                LFFCORDX,LFFCHARGX,LFFPOLX,LFFSIGX,LFFEPSX,
     *                LFFSIG2X,LFFEPS2X,LFFCLPRX,LFFZLPRX,LFFCHGIMX,
     *                LFFCRDIMX,LFFLISTAQM,LFFLISTBQM
      COMMON /FFQMPA/ ENFFQM2,ENPAV2,
     *                SCFTYP2,TDDFT2,MPLEVL2,CITYP2,
     *                ICHARG2,MULT2,IDOQM2,IREDOX,IQMPKA,IQMRXN,
     *                MATOMA,MCHARGA,MULTA,MELEA,
     *                MATOMB,MCHARGB,MULTB,MELEB,
     *                IECPX,NSHELLX,IMP,JMP,ICORSH,IGTF,
     *                LFFZANX,LFFCLPX,LFFZLPX,LFFNLPX,LFFKFRSTX,
     *                LFFKLASTX,LFFLMAXX,LFFLPSKIPX,LFFIZCOREX,
     *                LFFCX,LFFIANX,LFFEXX,LFFCSX,LFFCPX,LFFCDX,
     *                LFFCFX,LFFCGX,LFFCHX,LFFCIX,LFFKSTARTX,
     *                LFFKATOMX,LFFKTYPEX,LFFKNGX,LFFKLOCX,
     *                LFFMINX,LFFMAXX,LFFMPTYPX,LFFAN0X,
     *                LFFALPN0X,LFFAN1X,LFFALPN1X,LFFMPSKPX,
     *                LFFNOAN0X,LFFNOAN1X,LFFBPARX,LFFEXPMPX,
     *                LFFCSMPX,LFFCPMPX,LFFCDMPX,LFFCFMPX,LFFMPSKIPX,
     *                LFFNOCOSHX,LFFMPKSTAX,LFFMPKNGX,LFFMPKTYPX,
     *                LFFMPKMINX,LFFMPKMAXX,LFFMPKLOCX,LFFANAMX
      COMMON /FFRATT/ RATOLC,RATOLV,SCALRAT,VIRRAT(3),IRATTLE,JRATTLE,
     *                NRATTLE,MXRATT,LFFOLDCORD,LFFLSTRAT,LFFDSTRAT,
     *                LFFVELSV,IRATQM
      COMMON /FFRMDF/ RDIST1,RDIST2,RINTRV,NRATM1,NRATM2,INTRSTP
      COMMON /FFRMSD/ DIPT(8),TIMGYRA,TIMRALL,
     *                NATPDB,NGYRA,NDIEL,NRALL,
     *                LFFRALL0,NRMSD,LFFRMSD0,KFREEAB(201),
     *                NRIJMM,IJRMM(2,100),NRIJQM,IJRQM(2,100),
     *                NAIJKMM,IJKMM(3,100),NAIJKQM,IJKQM(3,100),
     *                NFIXMM,IFIXMM(200),NFIXQM,IFIXQM(200)
      COMMON /FFRXN / RXNEPS,RSPHSOL,ISPHSOL
      COMMON /FFSPH / SPHRAD,SPHSIG,SPHEPS,IADDWAT
      COMMON /FFUMBR/ UMBFC,UMBR0,UMBSIZE,
     *                NUMBBIN,NUMBATM(6),NUMBTYP,LFFUMBHIS,
     *                UM2FC,UM2R0,UM2SIZE,
     *                NUM2BIN,NUM2ATM(6),NUM2TYP,LFFUM2HIS
      COMMON /FMCOM / XX(1)
      COMMON /FUNCT / E,EG(3,MXATM)
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /INFOTD/ CNVTOL,PFREQ(2),MODTD,
     *                JANST,NRADT,NTHET,NPHIT,NLEBT,
     *                NSTAT,NTRIAL,MAXVEC,NTHST,IRECTD,ITDFG,ITDPRP,
     *                TRIPLET,SG1T,NONEQR,TAMMD,TPA,ALPHKWD,BETAKWD,
     *                SPCP(3),MULTD,MREKT,MRDEA,MTHST,IFEDAT(4)
      COMMON /IOFILE/ IR,IW,IP,IJK,IJKT,IDAF,NAV,IODA(950)
      COMMON /OUTPUT/ NPRINT,ITOL,ICUT,NORMF,NORMP,NOPK
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /RUNLAB/ TITLE(10),ANAM(MXATM),BNAM(MXATM),BFLAB(MXAO)
      COMMON /WFNOPT/ SCFTYP,VBTYP,DFTYPE,TDDFTYP,CITYP,CCTYP,
     *                MPLEVL,MPCTYP
C
      DATA RNONE       /8HNONE    /
C
C     NANDUN THELLAMUREGE, HUI LI, JAN 2011, LINCOLN
C     FENGCHAO CUI, HUI LI, JUN 2012 (ADD FREE ENERGY)
C
      IF(NSTEP.LT.0) RETURN
C
      CALL SEQOPN(36,'TRAJECT','NEW',.FALSE.,'FORMATTED')
      CALL SEQOPN(37,'QMWAVE' ,'NEW',.FALSE.,'FORMATTED')
      CALL SEQOPN(38,'MDDIP'  ,'NEW',.FALSE.,'FORMATTED')
      ISTEP  =  0
      NPRINT = -5
      IF(IADDWAT.EQ.2 .AND. ISPHSOL.GT.0) RSPHSOL=RSPHSOL*1.0D+40
C
C     -- SAVE NFIXQM COORDINATE AND VELOCITY
C
      DO KFIX=1,NFIXQM
         DO III = 1, 3
            FIXQM(III,KFIX)     =     C(III,IFIXQM(KFIX))
            FIXQM(III,KFIX+200) = QMVEL(III,IFIXQM(KFIX))
         ENDDO
      ENDDO
C
C     -- READY TO TAKE OFF ?
C
      MDSTEP = -1
      DO ISTEP = 0, NSTEP
         MDSTEP = MDSTEP + 1
         IF(ISTEP.EQ.100.AND.IADDWAT.EQ.2 .AND. ISPHSOL.GT.0)
     *   RSPHSOL = RSPHSOL*1.0D-40
         IF(MDSTEP.LE.1) CALL VCLR(UMBHIS,1,NUMBBIN)
         IF(MDSTEP.LE.1) CALL VCLR(UM2HIS,1,NUM2BIN*NUMBBIN)
C
         IF(MDSTEP.EQ.0) GOTO 50
C
         CALL DCOPY(3*NAT,C,1,OLDC,1)
         CALL DCOPY(3*NFFAT,CORD,1,OLDCORD,1)
C
C        - CALCULATE X(T+DT)
         IF(INTALG.EQ.1) THEN
            IF(MDSTEP.EQ.1) CALL DCOPY(3*NAT,QMGRD1,1,QMGRD0,1)
            DO IAT = 1, NAT
               DUM  = QM1MAS(IAT)*DT2
               DO III = 1, 3
                  C(III,IAT)=C(III,IAT)
     *                      +QMVEL(III,IAT)*DT
     *                      -TWOTHIRD*QMGRD1(III,IAT)*DUM
     *                      +ONESIX*QMGRD0(III,IAT)*DUM
               ENDDO
            ENDDO
            IF(MDSTEP.EQ.1) CALL DCOPY(3*NFFAT,FFGRD1,1,FFGRD0,1)
            DO IFFAT = 1, NFFAT
               DUM  = ONEMAS(IFFAT)*DT2
               DO III = 1, 3
                  CHANGE = +VEL(III,IFFAT)*DT
     *                     -TWOTHIRD*FFGRD1(III,IFFAT)*DUM
     *                     +ONESIX*FFGRD0(III,IFFAT)*DUM
                  CORD(III,IFFAT)=CORD(III,IFFAT) + CHANGE
               ENDDO
            ENDDO
         END IF
         IF(INTALG.EQ.2) THEN
            DO IAT = 1, NAT
               DUM  = QM1MAS(IAT)*DT2
               DO III = 1, 3
                  C(III,IAT)=C(III,IAT)
     *                      +QMVEL(III,IAT)*DT
     *                      -PT5*QMGRD1(III,IAT)*DUM
               ENDDO
            ENDDO
            DO IFFAT = 1, NFFAT
               DUM  = ONEMAS(IFFAT)*DT2
               DO III = 1, 3
                  CHANGE = +VEL(III,IFFAT)*DT
     *                     -PT5*FFGRD1(III,IFFAT)*DUM
                  CORD(III,IFFAT)=CORD(III,IFFAT) + CHANGE
               ENDDO
            ENDDO
         END IF
         IF(NACTQM.GT.0) THEN
            DO KOPT = 1, NACTQM
               IAT = LACTQM(KOPT)
               IF(IAT.GT.0) THEN
                  C(1,IAT)=(C(1,IAT)+3.1416D+00)*1.0D-40
                  C(2,IAT)=(C(2,IAT)+2.7183D+00)*1.0D-40
                  C(3,IAT)=(C(3,IAT)+1.4427D+00)*1.0D-40
               END IF
            ENDDO
            DO IAT = 1, NAT
               IF(ABS(C(1,IAT))+ABS(C(2,IAT))
     *           +ABS(C(3,IAT)).LT.1.0D-20) THEN
                  C(1,IAT) = C(1,IAT)*1.0D+40 - 3.1416D+00
                  C(2,IAT) = C(2,IAT)*1.0D+40 - 2.7183D+00
                  C(3,IAT) = C(3,IAT)*1.0D+40 - 1.4427D+00
               ELSE
                  DO III = 1, 3
                     C(III,IAT) = OLDC(III,IAT)
                  ENDDO
               END IF
            ENDDO
         END IF
         DO KFIX=1,NFIXQM
            DO III = 1, 3
               C(III,IFIXQM(KFIX)) = FIXQM(III,KFIX)
            ENDDO
         ENDDO
         IF(NACTMM.GT.0) THEN
            DO KOPT = 1, NACTMM
               IFFAT = LACTMM(KOPT)
               IF(IFFAT.GT.0) THEN
                  CORD(1,IFFAT)=(CORD(1,IFFAT)+3.1416D+00)*1.0D-40
                  CORD(2,IFFAT)=(CORD(2,IFFAT)+2.7183D+00)*1.0D-40
                  CORD(3,IFFAT)=(CORD(3,IFFAT)+1.4427D+00)*1.0D-40
               END IF
            ENDDO
            DO IFFAT = 1, NFFAT
               IF(ABS(CORD(1,IFFAT))+ABS(CORD(2,IFFAT))
     *           +ABS(CORD(3,IFFAT)).LT.1.0D-20) THEN
                  CORD(1,IFFAT) = CORD(1,IFFAT)*1.0D+40 - 3.1416D+00
                  CORD(2,IFFAT) = CORD(2,IFFAT)*1.0D+40 - 2.7183D+00
                  CORD(3,IFFAT) = CORD(3,IFFAT)*1.0D+40 - 1.4427D+00
               ELSE
                  DO III = 1, 3
                     CORD(III,IFFAT) = OLDCORD(III,IFFAT)
                  ENDDO
               END IF
            ENDDO
         END IF
         DO KFIX=1,NFIXMM
            DO III = 1, 3
               CORD(III,IFIXMM(KFIX)) = OLDCORD(III,IFIXMM(KFIX))
            ENDDO
         ENDDO
C
         CALL VCLR(VIRRAT,1,3)
         IF(NRATTLE.GT.0) THEN
C           - FIXED ATOMS DO NOT MOVE, SO MAKE THE MASS INFINITE
            IF(NACTMM.GT.0) THEN
               DO IFFAT = 1, NFFAT
                  ONEMAS(IFFAT) = 0.0D+00
               ENDDO
               DO KOPT = 1, NACTMM
                  IFFAT = LACTMM(KOPT)
                  IF(IFFAT.GT.0) ONEMAS(IFFAT) = 1.0D+00/ZMAS(IFFAT)
               ENDDO
            END IF
            DO KFIX = 1, NFIXMM
               KFFAT = IFIXMM(KFIX)
               IF(KFFAT.GT.0) ONEMAS(KFFAT) = 0.0D+00
            ENDDO
            CALL RATTLE1(CORD,OLDCORD,VEL,DSTRAT,ONEMAS,LSTRAT,ISTEP)
            DO KFIX = 1, NFIXMM
               KFFAT = IFIXMM(KFIX)
               IF(KFFAT.GT.0) ONEMAS(KFFAT) = 1.0D+00/ZMAS(KFFAT)
            ENDDO
            IF(NACTMM.GT.0) THEN
               DO IFFAT = 1, NFFAT
                  ONEMAS(IFFAT) = 1.0D+00/ZMAS(IFFAT)
               ENDDO
            END IF
C           - RATTLE CODE WORKS DIRECTLY ON MM ATOMS, BUT NOW PASS TO QM
            IF(LISTQM(NFFAT+NAT+1).GT.0)THEN
               DO IAT = 1, NAT
                  KFFAT = LISTQM(NFFAT+IAT)
                  IF(KFFAT.GT.0) THEN
                     C(1,IAT)     = CORD(1,KFFAT)
                     C(2,IAT)     = CORD(2,KFFAT)
                     C(3,IAT)     = CORD(3,KFFAT)
                     QMVEL(1,IAT) = VEL(1,KFFAT)
                     QMVEL(2,IAT) = VEL(2,KFFAT)
                     QMVEL(3,IAT) = VEL(3,KFFAT)
                  END IF
               ENDDO
            END IF
         END IF
C
         IF(NRATM1.GT.0) THEN
            CALL IRMDF1(CORD,OLDCORD,VEL,ONEMAS,ISTEP,RMDFCE)
C           - IRMDF CODE WORKS DIRECTLY ON MM ATOMS, BUT NOW PASS TO QM
            IF(LISTQM(NFFAT+NAT+1).GT.0)THEN
               DO IAT = 1, NAT
                  KFFAT = LISTQM(NFFAT+IAT)
                  IF(KFFAT.GT.0) THEN
                     C(1,IAT)     = CORD(1,KFFAT)
                     C(2,IAT)     = CORD(2,KFFAT)
                     C(3,IAT)     = CORD(3,KFFAT)
                     QMVEL(1,IAT) = VEL(1,KFFAT)
                     QMVEL(2,IAT) = VEL(2,KFFAT)
                     QMVEL(3,IAT) = VEL(3,KFFAT)
                  END IF
               ENDDO
            END IF
         END IF
C
  50     CONTINUE
C
         CALL NONBOND(MDSTEP,CORD,CORDSV,CORDSV2,CORDSVQ,
     *                NONLS1,NONLS2,
     *                NONLSTQ,LSTCELL,
     *                MVFASTS2,MVFASTS3,MVFASTS4,
     *                MVFASTL2,MVFASTL3,MVFASTL4,
     *                NONLSA,NONLSB,
     *                NONLSPMA,NONLSPMB)
C
         IWSAVE = IW
         IF(ISTEP.GT.0) IW = 37  ! PRINT TO QMWAVE
         IPSAVE = IP
         IP = 37  ! PUNCH TO QMWAVE
C
C        -- RUN A DIFFERENT QM CALCULATION (ONLY FOR ENERGY)
         IF(IDOQM2.EQ.1) THEN
            SCFTYPSV   = SCFTYP
            TDDFTSV    = TDDFTYP
            CITYPSV    = CITYP
            MPLEVLSV   = MPLEVL
            ICHARGSV   = ICH
            MULTSV     = MUL
            NESV       = NE
            NASV       = NA
            NBSV       = NB
            SCFTYP     = SCFTYP2
            TDDFTYP    = TDDFT2
            CITYP      = CITYP2
            MPLEVL     = MPLEVL2
            ICH        = ICHARG2
            MUL        = MULT2
            NE         = NESV - ICH + ICHARGSV
            NA         = (NE-1+MUL)/2
            NB         = (NE+1-MUL)/2
            IF(MPLEVL2.EQ.2) CALL MP2INP(MPTEST)
            IF(TDDFT2.NE.RNONE) CALL TDDINP
            CALL VCLR(QMGRD2,1,3*NAT  )
            CALL VCLR(FFGRD2,1,3*NFFAT)
            CALL VCLR(VIR   ,1,3)
            IF(IDIMER.EQ.0)THEN
               CALL GRADX
            ELSE
               CALL DIMERX
            END IF
C           -- ESCF WAS ZERO IF SCF NOT CONVERGED --
            IF(ABS(ESCF).LT.1.0D-12) THEN
               IF(MASWRK) WRITE(IWSAVE,'(/A/)')
     *            ' ERROR: SCF WAS NOT CONVERGED. QUANPOL MUST STOP.'
               CALL ABRT
            END IF
                                 ENFFQM2 = ETOT
            IF(MPLEVL.EQ.2)      ENFFQM2 = EMP2
            IF(TDDFTYP.NE.RNONE) ENFFQM2 = ESTATE(NTHST)
            SCFTYP     = SCFTYPSV
            TDDFTYP    = TDDFTSV
            CITYP      = CITYPSV
            MPLEVL     = MPLEVLSV
            ICH        = ICHARGSV
            MUL        = MULTSV
            NE         = NESV
            NA         = NASV
            NB         = NBSV
         END IF
C
C        -- CALCULATE ENERGY AND GRADIENT AT X(T+DT)
         CALL VCLR(QMGRD2,1,3*NAT  )
         CALL VCLR(FFGRD2,1,3*NFFAT)
         CALL VCLR(VIR   ,1,3)
         IF(IDIMER.EQ.0)THEN
            CALL GRADX
         ELSE
            CALL DIMERX
         END IF
C        -- ESCF WAS ZERO IF SCF NOT CONVERGED --
         IF(ABS(ESCF).LT.1.0D-12) THEN
            IF(MASWRK) WRITE(IWSAVE,'(/A/)')
     *         ' ERROR: SCF WAS NOT CONVERGED. QUANPOL MUST STOP.'
            CALL ABRT
         END IF
         CALL DCOPY(3*NAT,EG,1,QMGRD2,1)
         IF(IFEPTYP.GT.0.AND.MATOMB.GT.0) THEN
            CALL SAVEABPROP(1)
            CALL DCOPY(3*MATOMA,C,1,XX(LFFCX),1)
            CALL SAVEFFDATA
            CALL SETQMAB(2)
            CALL SETFFDATAB(2)
            CALL VCLR(QMGRD2,1,3*NAT  )
            CALL VCLR(FFGRD2,1,3*NFFAT)
            CALL VCLR(VIR   ,1,3)
            CALL GRADX
C           -- ESCF WAS ZERO IF SCF NOT CONVERGED --
            IF(ABS(ESCF).LT.1.0D-12) THEN
               IF(MASWRK) WRITE(IWSAVE,'(/A/)')
     *            ' ERROR: SCF WAS NOT CONVERGED. QUANPOL MUST STOP.'
               CALL ABRT
            END IF
            CALL DCOPY(3*NAT,EG,1,QMGRD2,1)
            CALL SAVEABPROP(2)
            CALL SETQMAB(1)
            CALL SETFFDATAB(1)
            CALL MIXABPROP
         END IF
C
         IW = IWSAVE
         IP = IPSAVE
C
C        - COMBINE QM AND MM GRADIENTS -
C
         IF(LISTQM(NFFAT+NAT+1).GT.0)THEN
            DO IAT = 1, NAT
               KFFAT = LISTQM(NFFAT+IAT)
               IF(KFFAT.GT.0) THEN
                  QMGRD2(1,IAT)   = QMGRD2(1,IAT) + FFGRD2(1,KFFAT)
                  QMGRD2(2,IAT)   = QMGRD2(2,IAT) + FFGRD2(2,KFFAT)
                  QMGRD2(3,IAT)   = QMGRD2(3,IAT) + FFGRD2(3,KFFAT)
                  FFGRD2(1,KFFAT) = QMGRD2(1,IAT)
                  FFGRD2(2,KFFAT) = QMGRD2(2,IAT)
                  FFGRD2(3,KFFAT) = QMGRD2(3,IAT)
               END IF
            ENDDO
         END IF
C
         IF(MDSTEP.EQ.0) THEN
            CALL DCOPY(3*NAT  ,QMGRD2,1,QMGRD1,1)
            CALL DCOPY(3*NFFAT,FFGRD2,1,FFGRD1,1)
            GOTO 51
         END IF
C
C        - FIXED ATOMS DO NOT MOVE, SO MAKE THE MASS INFINITE
         IF(NACTQM.GT.0) THEN
            DO IAT = 1, NAT
               QM1MAS(IAT) = 0.0D+00
            ENDDO
            DO KOPT = 1, NACTQM
               IAT = LACTQM(KOPT)
               IF(IAT.GT.0) QM1MAS(IAT) = 1.0D+00/QMZMAS(IAT)
            ENDDO
         END IF
         DO KFIX = 1, NFIXQM
            KAT = IFIXQM(KFIX)
            IF(KAT.GT.0) QM1MAS(KAT) = 0.0D+00
         ENDDO
         IF(NACTMM.GT.0) THEN
            DO IFFAT = 1, NFFAT
               ONEMAS(IFFAT) = 0.0D+00
            ENDDO
            DO KOPT = 1, NACTMM
               IFFAT = LACTMM(KOPT)
               IF(IFFAT.GT.0) ONEMAS(IFFAT) = 1.0D+00/ZMAS(IFFAT)
            ENDDO
         END IF
         DO KFIX = 1, NFIXMM
            KFFAT = IFIXMM(KFIX)
            IF(KFFAT.GT.0) ONEMAS(KFFAT) = 0.0D+00
         ENDDO
C
C        -- CALCULATE V(T+DT)
         IF(INTALG.EQ.1) THEN
            IF(MDSTEP.EQ.1) THEN
               DO IAT = 1, NAT
                  DUM  = QM1MAS(IAT)*DT
                  DO III = 1, 3
                     QMVEL(III,IAT)=QMVEL(III,IAT)
     *                            -PT5*(QMGRD1(III,IAT)
     *                                 +QMGRD2(III,IAT))*DUM
                     QMGRD1(III,IAT)=QMGRD2(III,IAT)
                  ENDDO
               ENDDO
               DO IFFAT = 1, NFFAT
                  DUM  = ONEMAS(IFFAT)*DT
                  DO III = 1, 3
                     VEL(III,IFFAT)=VEL(III,IFFAT)
     *                            -PT5*(FFGRD1(III,IFFAT)
     *                                 +FFGRD2(III,IFFAT))*DUM
                     FFGRD1(III,IFFAT)=FFGRD2(III,IFFAT)
                  ENDDO
               ENDDO
            ELSE
               DO IAT = 1, NAT
                  DUM  = QM1MAS(IAT)*DT
                  DO III = 1, 3
                     QMVEL(III,IAT)=QMVEL(III,IAT)
     *                           -(ONETHIRD*QMGRD2(III,IAT)
     *                             +FIVESIX*QMGRD1(III,IAT)
     *                              -ONESIX*QMGRD0(III,IAT))*DUM
                     QMGRD0(III,IAT)=QMGRD1(III,IAT)
                     QMGRD1(III,IAT)=QMGRD2(III,IAT)
                  ENDDO
               ENDDO
               DO IFFAT = 1, NFFAT
                  DUM  = ONEMAS(IFFAT)*DT
                  DO III = 1, 3
                        VEL(III,IFFAT)=VEL(III,IFFAT)
     *                               -(ONETHIRD*FFGRD2(III,IFFAT)
     *                                 +FIVESIX*FFGRD1(III,IFFAT)
     *                                  -ONESIX*FFGRD0(III,IFFAT))*DUM
                        FFGRD0(III,IFFAT)=FFGRD1(III,IFFAT)
                        FFGRD1(III,IFFAT)=FFGRD2(III,IFFAT)
                  ENDDO
               ENDDO
            END IF
         END IF
         IF(INTALG.EQ.2) THEN
            DO IAT = 1, NAT
               DUM  = QM1MAS(IAT)*DT
               DO III = 1, 3
                  QMVEL(III,IAT)=QMVEL(III,IAT)
     *                         -PT5*(QMGRD1(III,IAT)
     *                              +QMGRD2(III,IAT))*DUM
                  QMGRD1(III,IAT)=QMGRD2(III,IAT)
               ENDDO
            ENDDO
            DO IFFAT = 1, NFFAT
               DUM  = ONEMAS(IFFAT)*DT
               DO III = 1, 3
                  VEL(III,IFFAT)=VEL(III,IFFAT)
     *                         -PT5*(FFGRD1(III,IFFAT)
     *                              +FFGRD2(III,IFFAT))*DUM
                  FFGRD1(III,IFFAT)=FFGRD2(III,IFFAT)
               ENDDO
            ENDDO
         END IF
C
         IF(NRATTLE.GT.0) THEN
            CALL RATTLE2A(CORD,VEL,DSTRAT,ONEMAS,LSTRAT,ISTEP)
C           - RATTLE CODE WORKS DIRECTLY ON MM ATOMS, BUT NOW PASS TO QM
            IF(LISTQM(NFFAT+NAT+1).GT.0)THEN
               DO IAT = 1, NAT
                  KFFAT = LISTQM(NFFAT+IAT)
                  IF(KFFAT.GT.0) THEN
                     QMVEL(1,IAT) = VEL(1,KFFAT)
                     QMVEL(2,IAT) = VEL(2,KFFAT)
                     QMVEL(3,IAT) = VEL(3,KFFAT)
                  END IF
               ENDDO
            END IF
         END IF
C
         IF(NRATM1.GT.0) THEN
            CALL IRMDF2(CORD,VEL,ONEMAS,ISTEP)
C           - IRMDF CODE WORKS DIRECTLY ON MM ATOMS, BUT NOW PASS TO QM
            IF(LISTQM(NFFAT+NAT+1).GT.0)THEN
               DO IAT = 1, NAT
                  KFFAT = LISTQM(NFFAT+IAT)
                  IF(KFFAT.GT.0) THEN
                     QMVEL(1,IAT) = VEL(1,KFFAT)
                     QMVEL(2,IAT) = VEL(2,KFFAT)
                     QMVEL(3,IAT) = VEL(3,KFFAT)
                  END IF
               ENDDO
            END IF
         END IF
C
         DO KFIX = 1, NFIXQM
            KAT = IFIXQM(KFIX)
            IF(KAT.GT.0) QM1MAS(KAT) = 1.0D+00/QMZMAS(KAT)
         ENDDO
         DO KFIX = 1, NFIXMM
            KFFAT = IFIXMM(KFIX)
            IF(KFFAT.GT.0) ONEMAS(KFFAT) = 1.0D+00/ZMAS(KFFAT)
         ENDDO
         IF(NACTQM.GT.0) THEN
            DO IAT = 1, NAT
               QM1MAS(IAT) = 1.0D+00/QMZMAS(IAT)
            ENDDO
         END IF
         IF(NACTMM.GT.0) THEN
            DO IFFAT = 1, NFFAT
               ONEMAS(IFFAT) = 1.0D+00/ZMAS(IFFAT)
            ENDDO
         END IF
C
C        -- SYNCHRONIZE CORD AND VEL EVERY 200 STEPS --
         IF(MOD(ISTEP,200).EQ.0) THEN
            IF(LISTQM(NFFAT+NAT+1).GT.0)THEN
               DO IAT = 1, NAT
                  KFFAT = LISTQM(NFFAT+IAT)
                  IF(KFFAT.GT.0) THEN
                     C(1,IAT)     = CORD(1,KFFAT)
                     C(2,IAT)     = CORD(2,KFFAT)
                     C(3,IAT)     = CORD(3,KFFAT)
                     QMVEL(1,IAT) = VEL(1,KFFAT)
                     QMVEL(2,IAT) = VEL(2,KFFAT)
                     QMVEL(3,IAT) = VEL(3,KFFAT)
                  END IF
               ENDDO
            END IF
            IF(GOPARR) THEN
               CALL DDI_BCAST(461,'F',C,3*NAT,MASTER)
               CALL DDI_BCAST(462,'F',CORD,3*NFFAT,MASTER)
               CALL DDI_BCAST(463,'F',QMVEL,3*NAT,MASTER)
               CALL DDI_BCAST(464,'F',VEL,3*NFFAT,MASTER)
            END IF
         END IF
C
  51     CONTINUE
C
C        -- CALCULATE PROPERTIES
         CALL MDPROP(CORD,VEL,QMVEL,ZMAS,QMZMAS,MDSTEP,
     *               LISTQM,TEMP,PRES,UMBHIS,UM2HIS,
     *               PRESX,PRESY,PRESZ)
C
C        -- RELOCATE WATER MOLECULES AND IONS TO THE MASTER BOX --
C           MUST BE DONE AFTER MDPROP AND BEFORE PRINT OUT
C
         IF(KMASTER.GT.0.AND.MIN(XBOX,YBOX,ZBOX).LT.1.0D+30.AND.
     *      MOD(ISTEP,KMASTER).EQ.0) THEN
         DO IFFAT=1,NFFAT
         JYES1 = 0
         JYES2 = 0
         JYES3 = 0
         JYES4 = 0
         NUCZ  = INT(ZANF(IFFAT)+0.01D+00)
         Q     = CHARG(IFFAT)
C        - IONS AND MONOATOMIC GAS -
         IF((NUCZ.EQ. 2.AND.Q.EQ. 0.0D+00).OR.
     *      (NUCZ.EQ. 3.AND.Q.EQ.+1.0D+00).OR.
     *      (NUCZ.EQ. 9.AND.Q.EQ.-1.0D+00).OR.
     *      (NUCZ.EQ.10.AND.Q.EQ. 0.0D+00).OR.
     *      (NUCZ.EQ.11.AND.Q.EQ.+1.0D+00).OR.
     *      (NUCZ.EQ.12.AND.Q.EQ.+2.0D+00).OR.
     *      (NUCZ.EQ.13.AND.Q.EQ.+3.0D+00).OR.
     *      (NUCZ.EQ.17.AND.Q.EQ.-1.0D+00).OR.
     *      (NUCZ.EQ.18.AND.Q.EQ. 0.0D+00).OR.
     *      (NUCZ.EQ.19.AND.Q.EQ.+1.0D+00).OR.
     *      (NUCZ.EQ.20.AND.Q.EQ.+2.0D+00).OR.
     *      (NUCZ.GE.21.AND.NUCZ.LE.30.AND.Q.EQ.+1.0D+00).OR.
     *      (NUCZ.GE.21.AND.NUCZ.LE.30.AND.Q.EQ.+2.0D+00).OR.
     *      (NUCZ.GE.21.AND.NUCZ.LE.30.AND.Q.EQ.+3.0D+00).OR.
     *      (NUCZ.GE.21.AND.NUCZ.LE.30.AND.Q.EQ.+4.0D+00).OR.
     *      (NUCZ.EQ.35.AND.Q.EQ.-1.0D+00).OR.
     *      (NUCZ.EQ.36.AND.Q.EQ. 0.0D+00).OR.
     *      (NUCZ.EQ.37.AND.Q.EQ.+1.0D+00).OR.
     *      (NUCZ.EQ.53.AND.Q.EQ.-1.0D+00).OR.
     *      (NUCZ.EQ.54.AND.Q.EQ. 0.0D+00).OR.
     *      (NUCZ.EQ.55.AND.Q.EQ.+1.0D+00)     ) THEN
            CX    = CORD(1,IFFAT) - CENTX
            CY    = CORD(2,IFFAT) - CENTY
            CZ    = CORD(3,IFFAT) - CENTZ
            PBCX  = XBOX * ANINT(CX*ONEXBOX)
            PBCY  = YBOX * ANINT(CY*ONEYBOX)
            PBCZ  = ZBOX * ANINT(CZ*ONEZBOX)
            CORD(1,IFFAT) = CORD(1,IFFAT) - PBCX
            CORD(2,IFFAT) = CORD(2,IFFAT) - PBCY
            CORD(3,IFFAT) = CORD(3,IFFAT) - PBCZ
C           - ONLY DIFFUSION NEEDS THIS CORRECTION -
            IF(NDFS.GT.0) THEN
               DFSC0(1,IFFAT) = DFSC0(1,IFFAT) - PBCX
               DFSC0(2,IFFAT) = DFSC0(2,IFFAT) - PBCY
               DFSC0(3,IFFAT) = DFSC0(3,IFFAT) - PBCZ
            END IF
         END IF
         IF(NUCZ.EQ.8) THEN
            CX    = CORD(1,IFFAT) - CENTX
            CY    = CORD(2,IFFAT) - CENTY
            CZ    = CORD(3,IFFAT) - CENTZ
            PBCX  = XBOX * ANINT(CX*ONEXBOX)
            PBCY  = YBOX * ANINT(CY*ONEYBOX)
            PBCZ  = ZBOX * ANINT(CZ*ONEZBOX)
            IF((ABS(PBCX)+ABS(PBCY)+ABS(PBCZ)).GT.0.0D+00)THEN
               DO JFFAT=IFFAT+1,IFFAT+4
                  IF(INT(ZANF(JFFAT)+0.01D+00).EQ.1) THEN
                     CXH = CORD(1,IFFAT)-CORD(1,JFFAT)
                     CYH = CORD(2,IFFAT)-CORD(2,JFFAT)
                     CZH = CORD(3,IFFAT)-CORD(3,JFFAT)
                     R2H = CXH*CXH+CYH*CYH+CZH*CZH
C                    --  1.4 A  --
                     IF(R2H.LE.7.0D+00) THEN
                        IF(JYES3.GT.0.AND.JYES4.EQ.0) JYES4 = JFFAT
                        IF(JYES2.GT.0.AND.JYES3.EQ.0) JYES3 = JFFAT
                        IF(JYES1.GT.0.AND.JYES2.EQ.0) JYES2 = JFFAT
                        IF(JYES1.EQ.0)                JYES1 = JFFAT
                     END IF
                  END IF
               ENDDO
               IF(JYES1.GT.0.AND.JYES2.GT.0) THEN
                  CORD(1,IFFAT) = CORD(1,IFFAT) - PBCX
                  CORD(2,IFFAT) = CORD(2,IFFAT) - PBCY
                  CORD(3,IFFAT) = CORD(3,IFFAT) - PBCZ
                  CORD(1,JYES1) = CORD(1,JYES1) - PBCX
                  CORD(2,JYES1) = CORD(2,JYES1) - PBCY
                  CORD(3,JYES1) = CORD(3,JYES1) - PBCZ
                  CORD(1,JYES2) = CORD(1,JYES2) - PBCX
                  CORD(2,JYES2) = CORD(2,JYES2) - PBCY
                  CORD(3,JYES2) = CORD(3,JYES2) - PBCZ
                  IF(NDFS.GT.0) THEN
                     DFSC0(1,IFFAT) = DFSC0(1,IFFAT) - PBCX
                     DFSC0(2,IFFAT) = DFSC0(2,IFFAT) - PBCY
                     DFSC0(3,IFFAT) = DFSC0(3,IFFAT) - PBCZ
                     DFSC0(1,JYES1) = DFSC0(1,JYES1) - PBCX
                     DFSC0(2,JYES1) = DFSC0(2,JYES1) - PBCY
                     DFSC0(3,JYES1) = DFSC0(3,JYES1) - PBCZ
                     DFSC0(1,JYES2) = DFSC0(1,JYES2) - PBCX
                     DFSC0(2,JYES2) = DFSC0(2,JYES2) - PBCY
                     DFSC0(3,JYES2) = DFSC0(3,JYES2) - PBCZ
                  END IF
               END IF
               IF(JYES3.GT.0.AND.JYES4.GT.0) THEN
                  CORD(1,JYES3) = CORD(1,JYES3) - PBCX
                  CORD(2,JYES3) = CORD(2,JYES3) - PBCY
                  CORD(3,JYES3) = CORD(3,JYES3) - PBCZ
                  CORD(1,JYES4) = CORD(1,JYES4) - PBCX
                  CORD(2,JYES4) = CORD(2,JYES4) - PBCY
                  CORD(3,JYES4) = CORD(3,JYES4) - PBCZ
                  IF(NDFS.GT.0) THEN
                     DFSC0(1,JYES3) = DFSC0(1,JYES3) - PBCX
                     DFSC0(2,JYES3) = DFSC0(2,JYES3) - PBCY
                     DFSC0(3,JYES3) = DFSC0(3,JYES3) - PBCZ
                     DFSC0(1,JYES4) = DFSC0(1,JYES4) - PBCX
                     DFSC0(2,JYES4) = DFSC0(2,JYES4) - PBCY
                     DFSC0(3,JYES4) = DFSC0(3,JYES4) - PBCZ
                  END IF
               END IF
            END IF
         END IF
         ENDDO
C        - FOR QM/MM, PASS TO QM
         IF(LISTQM(NFFAT+NAT+1).GT.0)THEN
            DO IAT = 1, NAT
               KFFAT = LISTQM(NFFAT+IAT)
               IF(KFFAT.GT.0) THEN
                  C(1,IAT)     = CORD(1,KFFAT)
                  C(2,IAT)     = CORD(2,KFFAT)
                  C(3,IAT)     = CORD(3,KFFAT)
               END IF
            ENDDO
C           -- ADJUST QMCX,Y,Z, ALSO UPDATE LQMCT
            IF(NAT.GT.0 .AND. NFFAT.GT.0 .AND. SWRB2.LT.1.0D+08) THEN
               XMAX = -1.0D+30
               YMAX = -1.0D+30
               ZMAX = -1.0D+30
               XMIN =  1.0D+30
               YMIN =  1.0D+30
               ZMIN =  1.0D+30
               DO IAT = 1,NAT
                  XMAX = MAX(XMAX,C(1,IAT))
                  YMAX = MAX(YMAX,C(2,IAT))
                  ZMAX = MAX(ZMAX,C(3,IAT))
                  XMIN = MIN(XMIN,C(1,IAT))
                  YMIN = MIN(YMIN,C(2,IAT))
                  ZMIN = MIN(ZMIN,C(3,IAT))
               ENDDO
               QMSIZE = ZERO
               QMSIZE = MAX(QMSIZE,XMAX-XMIN)
               QMSIZE = MAX(QMSIZE,YMAX-YMIN)
               QMSIZE = MAX(QMSIZE,ZMAX-ZMIN)
               QMSIZE = QMSIZE*1.732D+00
               QMCX   = (XMAX+XMIN)*PT5
               QMCY   = (YMAX+YMIN)*PT5
               QMCZ   = (ZMAX+ZMIN)*PT5
               R2NEAR = 100.0D+00
               NEAR   = 0
               DO IAT = 1,NAT
                  XI = C(1,IAT) - QMCX
                  YI = C(2,IAT) - QMCY
                  ZI = C(3,IAT) - QMCZ
                  R2 = XI**2 + YI**2 + ZI**2
                  IF(R2.LT.R2NEAR) THEN
                     R2NEAR = R2
                     NEAR   = IAT
                  END IF
               ENDDO
               LQMCT = NEAR
            END IF
         END IF
         END IF
C
         CALL TMDATE(TIMSTR)
         IF(MOD(ISTEP,KOUT).EQ.0) CALL TIMIT(1)
         IF(MOD(ISTEP,KOUT).NE.0 .AND. ISTEP.NE.NSTEP) GOTO 200
C
         IF(MASWRK)THEN
C
C           -- PRINT KOUTACT ATOMS --
C              HERE ONLY $FFDATA IS CONSIDERED, $FFDATB IS NOT.
C
            IF(KOUTACT(1).GT.0) THEN
            WRITE(IW,'(/A,I10)')'!KOUTACT ATOMS AROUND $FFDATA ATOM ',
     *      KOUTACT(1)
            WRITE(IW,'(1X,A,I10,10X,A,F19.2,A)')
     *      '$FFDATA ! KOUTACT     MD STEP',ISTEP,
     *      'TIME=',ISTEP*DT*2.418884326505D-02,' FS'
            ACTX = CORD(1,KOUTACT(1))
            ACTY = CORD(2,KOUTACT(1))
            ACTZ = CORD(3,KOUTACT(1))
            ACTR = DBLE(KOUTACT(2))*1.0D-08*TOBOHR
            ACTR2= ACTR**2
            DO IFFAT=1,NFFAT
               CX    = CORD(1,IFFAT) - ACTX
               CY    = CORD(2,IFFAT) - ACTY
               CZ    = CORD(3,IFFAT) - ACTZ
               PBCX  = XBOX*ANINT(CX*ONEXBOX)
               PBCY  = YBOX*ANINT(CY*ONEYBOX)
               PBCZ  = ZBOX*ANINT(CZ*ONEZBOX)
               CX    = CX - PBCX
               CY    = CY - PBCY
               CZ    = CZ - PBCZ
               RCXYZ = CX**2 + CY**2 + CZ**2
               IF(RCXYZ.LE.ACTR2) THEN
                  WRITE(IW,1000)ATMNAM(IFFAT),ZANF(IFFAT),
     *            (CX+ACTX)*TOANGS,(CY+ACTY)*TOANGS,(CZ+ACTZ)*TOANGS
               END IF
            ENDDO
            WRITE(IW,*)'$END ! KOUTACT'
            END IF
C
            WRITE(36,'(1X,A,I10,1X,A,F15.2,A)')
     *      'RESTART COORD + VEL FOR QUANPOL AT MD STEP',ISTEP,
     *      'TIME=',ISTEP*DT*2.418884326505D-02,' FS'
            IF(MIN(XBOX,YBOX,ZBOX).LT.1.0D+30) THEN
              WRITE(36,'(A,3(A,F15.10,1X),/9X,3(A,F15.10,1X),A)')
     *        ' $QUANPO ',
     *        'CENTX=',CENTX*TOANGS,
     *        'CENTY=',CENTY*TOANGS,
     *        'CENTZ=',CENTZ*TOANGS,
     *        'XBOX=',XBOX*TOANGS,
     *        'YBOX=',YBOX*TOANGS,
     *        'ZBOX=',ZBOX*TOANGS,
     *        '$END'
            ELSE IF(SPHRAD.LT.1.0D+30) THEN
              WRITE(36,'(A,3(A,F15.10,1X),/9X,2(A,F15.10,1X),A)')
     *        ' $QUANPO ',
     *        'CENTX=',CENTX*TOANGS,
     *        'CENTY=',CENTY*TOANGS,
     *        'CENTZ=',CENTZ*TOANGS,
     *        'SPHRAD=',SPHRAD*TOANGS,
     *        'RSPHSOL=',RSPHSOL*TOANGS,
     *        '$END'
            ELSE
              WRITE(36,'(A,3(A,F15.10,1X),A)')
     *        ' $QUANPO ',
     *        'CENTX=',CENTX*TOANGS,
     *        'CENTY=',CENTY*TOANGS,
     *        'CENTZ=',CENTZ*TOANGS,
     *        '$END'
            END IF
C
            WRITE(36,*)'$DATA'
            WRITE(36,'(A,3(F8.3,1X),A,I3,A,I10)')
     *      ' QM CENTER = ',QMCX*TOANGS,QMCY*TOANGS,QMCZ*TOANGS,
     *      'NEAR ATOM',LQMCT,'     AT MD STEP',ISTEP
            WRITE(36,*)'C1'
            DO IAT = 1, NAT
               WRITE(36,999)ANAM(IAT),ZAN(IAT),
     *         C(1,IAT)*TOANGS,
     *         C(2,IAT)*TOANGS,
     *         C(3,IAT)*TOANGS
            ENDDO
            IF(IFEPTYP.GT.0)THEN
               DO IATB = 1, MATOMB
                  WRITE(36,999)ANAM(MATOMA+IATB),  ! ANAM IS THE ORIGINAL
     *            XX(LFFZANX +MATOMA+IATB-1),
     *            XX(LFFCX+3*(MATOMA+IATB-1)  )*TOANGS,
     *            XX(LFFCX+3*(MATOMA+IATB-1)+1)*TOANGS,
     *            XX(LFFCX+3*(MATOMA+IATB-1)+2)*TOANGS
               ENDDO
            END IF
            WRITE(36,*)'$END'
            WRITE(36,'(1X,A,I10,10X,A,F19.2,A)')
     *      '$FFDATA          !    MD STEP',ISTEP,
     *      'TIME=',ISTEP*DT*2.418884326505D-02,' FS'
            IF(IFEPTYP.EQ.1) THEN
               WRITE(36,'(1X,A,I10,10X,A,F19.2,A)')
     *         '$FFDATB          !    MD STEP',ISTEP,
     *         'TIME=',ISTEP*DT*2.418884326505D-02,' FS'
            END IF
            WRITE(36,*)'COORDINATES  NUC                   X',
     *          '                   Y                   Z'
            DO IFFAT = 1, N1FFAT
               WRITE(36,1000)ATMNAM(IFFAT),ZANF(IFFAT),
     *         CORD(1,IFFAT)*TOANGS,
     *         CORD(2,IFFAT)*TOANGS,
     *         CORD(3,IFFAT)*TOANGS
            ENDDO
            WRITE(36,*)'STOP'
            WRITE(36,*)'QMVELOCITY                    VX',
     *          '                    VY                    VZ'
            DO IAT = 1, NAT
               WRITE(36,1002) ANAM(IAT),
     *         QMVEL(1,IAT),QMVEL(2,IAT),QMVEL(3,IAT)
            ENDDO
            IF(IFEPTYP.GT.0)THEN
               DO IATB = 1, MATOMB
                  CALL GETIFFAT(XX(LFFLISTBQM),IATB,IATA)
                  WRITE(36,1002) ANAM(MATOMA+IATB),
     *            QMVEL(1,IATA),QMVEL(2,IATA),QMVEL(3,IATA)
               ENDDO
            END IF
            WRITE(36,*)'STOP'
            WRITE(36,*)'MMVELOCITY                    VX',
     *          '                    VY                    VZ'
            DO IFFAT = 1, N1FFAT
               WRITE(36,1001) ATMNAM(IFFAT),
     *         VEL(1,IFFAT),VEL(2,IFFAT),VEL(3,IFFAT)
            ENDDO
            WRITE(36,*)'STOP'
C
            IF(IFEPTYP.EQ.2) THEN
               WRITE(36,'(1X,A,I10,10X,A,F19.2,A)')
     *         '$FFDATB          !    MD STEP',ISTEP,
     *         'TIME=',ISTEP*DT*2.418884326505D-02,' FS'
               WRITE(36,*)'COORDINATES  NUC                   X',
     *             '                   Y                   Z'
               DO IFFAT = 1, N2FFAT
                  WRITE(36,1000)ATMNAM(IFFAT),ZANF(IFFAT),
     *            CORDB(1,IFFAT)*TOANGS,
     *            CORDB(2,IFFAT)*TOANGS,
     *            CORDB(3,IFFAT)*TOANGS
               ENDDO
               WRITE(36,*)'STOP'
               WRITE(36,*)'QMVELOCITY                    VX',
     *             '                    VY                    VZ'
               DO IAT = 1, NAT
                  WRITE(36,1002) ANAM(IAT),
     *            QMVEL(1,IAT),QMVEL(2,IAT),QMVEL(3,IAT)
               ENDDO
C              IF(IFEPTYP.GT.0)THEN  ! ALWAYS DO THIS
                  DO IATB = 1, MATOMB
                     CALL GETIFFAT(XX(LFFLISTBQM),IATB,IATA)
                     WRITE(36,1002) ANAM(MATOMA+IATB),
     *               QMVEL(1,IATA),QMVEL(2,IATA),QMVEL(3,IATA)
                  ENDDO
C              END IF
               WRITE(36,*)'STOP'
               WRITE(36,*)'MMVELOCITY                    VX',
     *             '                    VY                    VZ'
               DO IFFAT = 1, N2FFAT
                  WRITE(36,1001) ATMNAM(IFFAT),
     *            VEL(1,IFFAT),VEL(2,IFFAT),VEL(3,IFFAT)
               ENDDO
               WRITE(36,*)'STOP'
            END IF
C
         END IF
C
         CALL FLSHBF(IW)
         CALL FLSHBF(36)
         CALL FLSHBF(37)
         CALL TIMIT(1)
 200     CONTINUE
C
C        -- SCALE VELOCITY AND VOLUME --
         CALL DCOPY(3*NAT,C,1,OLDC,1)
         CALL DCOPY(3*NAT,QMVEL,1,QMVELSV,1)
         CALL DCOPY(3*NFFAT,CORD,1,OLDCORD,1)
         CALL DCOPY(3*NFFAT,VEL,1,VELSV,1)
         CALL TPSTAT(CORD,VEL,QMVEL,ZMAS,QMZMAS,ONEMAS,QM1MAS,
     *               MDSTEP,LISTQM,TEMP,PRES,
     *               PRESX,PRESY,PRESZ,NONLSTQ)
         IF(NACTQM.GT.0) THEN
            DO KOPT = 1, NACTQM
               IAT = LACTQM(KOPT)
               IF(IAT.GT.0) THEN
                  C(1,IAT)=(C(1,IAT)+3.1416D+00)*1.0D-40
                  C(2,IAT)=(C(2,IAT)+2.7183D+00)*1.0D-40
                  C(3,IAT)=(C(3,IAT)+1.4427D+00)*1.0D-40
                  QMVEL(1,IAT)=QMVEL(1,IAT) + 1.0D+00
               END IF
            ENDDO
            DO IAT = 1, NAT
               IF(ABS(C(1,IAT))+ABS(C(2,IAT))
     *           +ABS(C(3,IAT)).LT.1.0D-20) THEN
                  C(1,IAT) = C(1,IAT)*1.0D+40 - 3.1416D+00
                  C(2,IAT) = C(2,IAT)*1.0D+40 - 2.7183D+00
                  C(3,IAT) = C(3,IAT)*1.0D+40 - 1.4427D+00
               ELSE
                  DO III = 1, 3
                     C(III,IAT) = OLDC(III,IAT)
                  ENDDO
               END IF
               IF(QMVEL(1,IAT).GT.0.5D+00) THEN
                  QMVEL(1,IAT) = QMVEL(1,IAT) - 1.0D+00
               ELSE
                  DO III = 1, 3
                     QMVEL(III,IAT) = QMVELSV(III,IAT)
                  ENDDO
               END IF
            ENDDO
         END IF
         DO KFIX=1,NFIXQM
            DO III = 1, 3
               C(III,IFIXQM(KFIX))     = FIXQM(III,KFIX)
               QMVEL(III,IFIXQM(KFIX)) = FIXQM(III,KFIX+200)
            ENDDO
         ENDDO
         IF(NACTMM.GT.0) THEN
            DO KOPT = 1, NACTMM
               IFFAT = LACTMM(KOPT)
               IF(IFFAT.GT.0) THEN
                  CORD(1,IFFAT)=(CORD(1,IFFAT)+3.1416D+00)*1.0D-40
                  CORD(2,IFFAT)=(CORD(2,IFFAT)+2.7183D+00)*1.0D-40
                  CORD(3,IFFAT)=(CORD(3,IFFAT)+1.4427D+00)*1.0D-40
                  VEL (1,IFFAT)=VEL(1,IFFAT) + 1.0D+00
               END IF
            ENDDO
            DO IFFAT = 1, NFFAT
               IF(ABS(CORD(1,IFFAT))+ABS(CORD(2,IFFAT))
     *           +ABS(CORD(3,IFFAT)).LT.1.0D-20) THEN
                  CORD(1,IFFAT) = CORD(1,IFFAT)*1.0D+40 - 3.1416D+00
                  CORD(2,IFFAT) = CORD(2,IFFAT)*1.0D+40 - 2.7183D+00
                  CORD(3,IFFAT) = CORD(3,IFFAT)*1.0D+40 - 1.4427D+00
               ELSE
                  DO III = 1, 3
                     CORD(III,IFFAT) = OLDCORD(III,IFFAT)
                  ENDDO
               END IF
               IF(VEL(1,IFFAT).GT.0.5D+00) THEN
                  VEL(1,IFFAT) = VEL(1,IFFAT) - 1.0D+00
               ELSE
                  DO III = 1, 3
                     VEL(III,IFFAT) = VELSV(III,IFFAT)
                  ENDDO
               END IF
            ENDDO
         END IF
         DO KFIX=1,NFIXMM
            DO III = 1, 3
               CORD(III,IFIXMM(KFIX)) = OLDCORD(III,IFIXMM(KFIX))
               VEL (III,IFIXMM(KFIX)) = VELSV(III,IFIXMM(KFIX))
            ENDDO
         ENDDO
C        -- ADJUST QMCX,Y,Z WHEN T OR V ARE ADJUSTED --
C           ALSO UPDATE LQMCT
C           (DO THIS ONLY IF QM ATOMS ARE SCALED)
         IF(ITSTAT.GT.0 .OR. IPSTAT.GT.0) THEN
            IF(NAT.GT.0 .AND. NFFAT.GT.0 .AND. SWRB2.LT.1.0D+08) THEN
               XMAX = -1.0D+30
               YMAX = -1.0D+30
               ZMAX = -1.0D+30
               XMIN =  1.0D+30
               YMIN =  1.0D+30
               ZMIN =  1.0D+30
               DO IAT = 1,NAT
                  XMAX = MAX(XMAX,C(1,IAT))
                  YMAX = MAX(YMAX,C(2,IAT))
                  ZMAX = MAX(ZMAX,C(3,IAT))
                  XMIN = MIN(XMIN,C(1,IAT))
                  YMIN = MIN(YMIN,C(2,IAT))
                  ZMIN = MIN(ZMIN,C(3,IAT))
               ENDDO
               QMSIZE = ZERO
               QMSIZE = MAX(QMSIZE,XMAX-XMIN)
               QMSIZE = MAX(QMSIZE,YMAX-YMIN)
               QMSIZE = MAX(QMSIZE,ZMAX-ZMIN)
               QMSIZE = QMSIZE*1.732D+00
               QMCX   = (XMAX+XMIN)*PT5
               QMCY   = (YMAX+YMIN)*PT5
               QMCZ   = (ZMAX+ZMIN)*PT5
               R2NEAR = 100.0D+00
               NEAR   = 0
               DO IAT = 1,NAT
                  XI = C(1,IAT) - QMCX
                  YI = C(2,IAT) - QMCY
                  ZI = C(3,IAT) - QMCZ
                  R2 = XI**2 + YI**2 + ZI**2
                  IF(R2.LT.R2NEAR) THEN
                     R2NEAR = R2
                     NEAR   = IAT
                  END IF
               ENDDO
               LQMCT = NEAR
            END IF
         END IF
C        -- SWAP A QM WATER MOLECULE WITH A MM WATER MOLECULE --
         IF(ISWAP.GT.0) THEN
            IF(MOD(ISTEP,ISWAP).EQ.0) THEN
               CALL SWAPWATER(CORD,ZANF,VEL,QMVEL,LISTQM,ATMNAM)
            END IF
         END IF
      ENDDO
C
  999 FORMAT(1X,A8,3X,F5.1,1X,F19.13,1X,F19.13,1X,F19.13)
 1000 FORMAT(1X,A10,1X,F5.1,1X,F19.13,1X,F19.13,1X,F19.13)
 1001 FORMAT(1X,A10,1X,F21.18,1X,F21.18,1X,F21.18)
 1002 FORMAT(1X,A8,3X,F21.18,1X,F21.18,1X,F21.18)
C
      CALL SEQCLO(36,'KEEP')
      CALL SEQCLO(37,'KEEP')
      CALL SEQCLO(38,'KEEP')
C
      RETURN
      END
C*MODULE QUANPOC  *DECK MDPROP
!>
!> @brief    MD properties
!>
!> @author   Nandun Thellamurege, Hui Li
!>           - Mar 2011
!>
!> @details  calculate and average MD properties
!>           pressure, volume, temperature ...
!>
      SUBROUTINE MDPROP(CORD,VEL,QMVEL,ZMAS,QMZMAS,ISTEP,
     *                  LISTQM,TEMP,PRES,UMBHIS,UM2HIS,
     *                  PRESX,PRESY,PRESZ)
      use mx_limits, only: mxatm,mxrt
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      LOGICAL TRIPLET,SG1T,TAMMD,TPA,ALPHKWD,BETAKWD
      LOGICAL MREKT,MRDEA
      LOGICAL GOPARR,DSKWRK,MASWRK
C
      PARAMETER (PI=3.14159265358979323846264338D+00)
      PARAMETER (TOKCAL=627.509469D+00)
      PARAMETER (DEGREE=57.2957795130823D+00)
      PARAMETER (TORAD=1.0D+00/DEGREE)
      PARAMETER (TOANGS=0.52917724924D+00)
      PARAMETER (TOBOHR=1.0D+00/TOANGS)
      PARAMETER (TOKELVIN=3.15774646D+05)
      PARAMETER (BOLTZK=1.0D+00/TOKELVIN)
      PARAMETER (ZERO=0.0D+00)
      PARAMETER (PT5=0.50D+00)
      PARAMETER (ONE=1.0D+00)
      PARAMETER (FOURTHIRD=4.0D+00/3.0D+00)
      PARAMETER (ONETHIRD=1.0D+00/3.0D+00)
      PARAMETER (TOBAR=2.942191219D+08)
      PARAMETER (TOEV=27.21138386D+00)
C
      DIMENSION CORD(3,*),VEL(3,*),QMVEL(3,*),ZMAS(*),QMZMAS(*),
     *          LISTQM(*),UMBHIS(*),UM2HIS(NUM2BIN,*)
C
      COMMON /ENRGMP/ EMP2,EMP3,EMP4,EMP2A
      COMMON /ENRGYS/ ENUCR,EELCT,ETOT,SZ,SZZ,ECORE,ESCF,EERD,E1,E2,
     *                VEN,VEE,EPOT,EKIN,ESTATE(MXRT),STATN,EDFT(3),EDISP
      COMMON /FFDFS / TIMDFS,QDION,AMION,TEFF,NDFS,NATMGAS,
     *                LFFDFSC,
     *                LFFDFSC0,LFFDFSA,LFFDFSN,LFFDFCOM,KDFS,LFFDFSCAV
      COMMON /FFENGY/ EN12,EN123,EN123R4,EN123B4,EN234W1,ENCHAR,ENLJR,
     *                ENLJD,ENPOL,XENPOL,ENRXN,XENRXN,ENRXNPOL,ENRXNR,
     *                EN12312,ENQUANP(30),
     *                ENBIAS,ENCENT,ENUCCH,ENCMAP,ENPOT,ENKIN,ENTOT
      COMMON /FFEWLD/ SPLIT,IEWALD,KEWALD,NKVEC,
     *                L1KVEC,L2KVEC,LFFRKEXPEL,LFFRKVEC,
     *                LFFKVEC,LFFTCHCH,LFFCOSCH,LFFSINCH
      COMMON /FFFIXS/ ENFIXSO,FIXEPS,FIXTOL,FIXA,FIXQ,RALLMM,RALLQM,
     *                RADMM(200),RADQM(200),NRADMM,NRADQM,IFIXSOL,
     *                LFFDAI,LFFDAIT,LFFIDDAI,LFFIDTMP,LFFTMPTS,
     *                LFFAFIX,LFFIDATOM,LFFRFIX,LFFQFIX,NTSATM,
     *                LFFQFIXMP,LFFQFIXTA,LFFQFIXXY,
     *                LFFXTSFIX,LFFYTSFIX,LFFZTSFIX,
     *                LFFVFIX1,LFFVFIX2,NCYCLE,MXFFTS,NFFTS
      COMMON /FFFREE/ SOL1CH,SOL1LJ,SOL1IM,SOLFRE1MM,
     *                SOL2CH,SOL2LJ,SOL2IM,SOLFRE2MM,
     *                ASOL1CH,ASOL1LJ,ASOL1IM,ASOL1MM,
     *                ASOL2CH,ASOL2LJ,ASOL2IM,ASOL2MM,
     *                PMF1BD,PMF1AG,PMF1DR,PMF1DB,PMF1WG,PMF1CM,
     *                PMF1CH,PMF1PO,PMF1LJ,PMF1IM,PMFFRE1MM,
     *                EN12A,EN12B,EN123A,EN123B,EN123R4A,EN123R4B,
     *                EN234W1A,EN234W1B,EN123B4A,EN123B4B,
     *                EN12312A,EN12312B,
     *                ENCHARA,ENCHARB,ENPOLA,ENPOLB,
     *                ENRXNA,ENRXNB,ENRXNPOLA,ENRXNPOLB,
     *                ENRXNRA,ENRXNRB,ENLJRA,ENLJRB,ENLJDA,ENLJDB,
     *                ENCMAPA,ENCMAPB,ENFIXSOA,ENFIXSOB,
     *                ENCENTA,ENCENTB,WSIMUL,WPERT1,WPERT2,
     *                IFEPTYP
      COMMON /FFPARM/ NFFAT,NBOND,NANGL,NDIHR,NDIHB,NCMAP,NWAGG,
     *                N1213J,N14J,NLKQMM,IDOCHG,IDOPOL,IDOLJ,IDOCMAP
      COMMON /FFMDPA/ DT,DT2,TEMP0,PRES0,POLTOL,VIR(3),PMEAN,VOLAV,
     *                ENPAV,ENKAV,TEMPAV,BERENDT,BERENDP,VELMAX,
     *                PMEANX,PMEANY,PMEANZ,NSTEP,KMASTER,KOUTACT(2),
     *                IHESS,INTALG,ITSTAT,IPSTAT,JOUT,KOUT,LOUT
      COMMON /FFMPNT/ LFFATMNAM,LFFCORD,LFFZANF,
     *                LFFZMAS,LFFONEMAS,LFFQMZMAS,LFFQM1MAS,
     *                LFFCHARG,LFFPOL,LFFDIP,
     *                LFFFIELD1,LFFFIELD2,LFFFIELD3,
     *                LFFSIG,LFFEPS,LFFSIG2,LFFEPS2,
     *                LFFBOND0,LFFFCBOND,
     *                LFFANGL0,LFFFCANGL,LFFFCWAGG,
     *                LFFDIHB0,LFFFCDIHB,
     *                LFFVROT,LFFNNN,LFFGAMA,LFFIPAIR,
     *                LFFKLIST,LFFLLIST,LFFL1213J,LFFL14J,
     *                LFFMLIST,LFFNLIST,LFFLKQMMM,
     *                LFFVEL,LFFQMVEL,
     *                LFFFFGRD0,LFFFFGRD1,LFFFFGRD2,
     *                LFFQMGRD0,LFFQMGRD1,LFFQMGRD2,LFFDETMP,
     *                LFFCLPR,LFFZLPR,LFFNLPR,
     *                LFFXTS,LFFYTS,LFFZTS,LFFCMAT1,
     *                LFFQRXN1,LFFQRXN2,LFFPOT1,LFFPOT2,LFFQRXNMP,
     *                LFFQRXNTA,LFFQRXNXY,LFFNONLSTQ,
     *                LFFDIPMP,LFFDIPTA,LFFDIPXY,LFFLISTQM,LFFNONLS1,
     *                LFFMAPLST,LFFCMAPCO
      COMMON /FFPBSW/ XBOX,YBOX,ZBOX,SWF,SWFDX,SWFDY,SWFDZ,
     *                SWRA,ONESWRA,SWRA2,ONESWRA2,
     *                SWRB,ONESWRB,SWRB2,ONESWRB2,ONESWRB4,
     *                SWFDUM3,SWFDUM4,SWFDUM5,
     *                SWRAQ,ONESWRAQ,SWRAQ2,ONESWRAQ2,
     *                SWRBQ,ONESWRBQ,SWRBQ2,ONESWRBQ2,ONESWRBQ4,
     *                SWFDUM3Q,SWFDUM4Q,SWFDUM5Q,
     *                QMSIZE,QMCX,QMCY,QMCZ,QMCXSV,QMCYSV,QMCZSV,
     *                CENTX,CENTY,CENTZ,BUFWID1,BUFWID2,RDAMP,
     *                EFIELDX,EFIELDY,EFIELDZ,QMCXSV2,QMCYSV2,QMCZSV2,
     *                EPS1RB,EPS1RB3,ONEXBOX,ONEYBOX,ONEZBOX,
     *                LQMCT,MXLIST1,NTODO,NTODOSV,NTODOQ,
     *                ISWITCH,ISHIFT,IPOLSHF,
     *                LFFLSTCELL,LFFCORDSV,
     *                LFFPOLSV,LFFCORDSV2,LFFNONLS2,LFFCORDSVQ,
     *                LFFMVFASTS2,LFFMVFASTS3,LFFMVFASTS4,
     *                LFFMVFASTL2,LFFMVFASTL3,LFFMVFASTL4,
     *                MXCHECK,MXLIST2,NTODO2,NTODO2SV
      COMMON /FFQMFE/ ETOTQA,EMP2QA,ESTATEQA,ETOTQB,EMP2QB,ESTATEQB,
     *                FRE1QMMM,FRE2QMMM,
     *                VIRA(3),VIRB(3),
     *                LFFFFGRDA,LFFFFGRDB,LFFQMGRDA,LFFQMGRDB,
     *                LFFCORDX,LFFCHARGX,LFFPOLX,LFFSIGX,LFFEPSX,
     *                LFFSIG2X,LFFEPS2X,LFFCLPRX,LFFZLPRX,LFFCHGIMX,
     *                LFFCRDIMX,LFFLISTAQM,LFFLISTBQM
      COMMON /FFQMPA/ ENFFQM2,ENPAV2,
     *                SCFTYP2,TDDFT2,MPLEVL2,CITYP2,
     *                ICHARG2,MULT2,IDOQM2,IREDOX,IQMPKA,IQMRXN,
     *                MATOMA,MCHARGA,MULTA,MELEA,
     *                MATOMB,MCHARGB,MULTB,MELEB,
     *                IECPX,NSHELLX,IMP,JMP,ICORSH,IGTF,
     *                LFFZANX,LFFCLPX,LFFZLPX,LFFNLPX,LFFKFRSTX,
     *                LFFKLASTX,LFFLMAXX,LFFLPSKIPX,LFFIZCOREX,
     *                LFFCX,LFFIANX,LFFEXX,LFFCSX,LFFCPX,LFFCDX,
     *                LFFCFX,LFFCGX,LFFCHX,LFFCIX,LFFKSTARTX,
     *                LFFKATOMX,LFFKTYPEX,LFFKNGX,LFFKLOCX,
     *                LFFMINX,LFFMAXX,LFFMPTYPX,LFFAN0X,
     *                LFFALPN0X,LFFAN1X,LFFALPN1X,LFFMPSKPX,
     *                LFFNOAN0X,LFFNOAN1X,LFFBPARX,LFFEXPMPX,
     *                LFFCSMPX,LFFCPMPX,LFFCDMPX,LFFCFMPX,LFFMPSKIPX,
     *                LFFNOCOSHX,LFFMPKSTAX,LFFMPKNGX,LFFMPKTYPX,
     *                LFFMPKMINX,LFFMPKMAXX,LFFMPKLOCX,LFFANAMX
      COMMON /FFRATT/ RATOLC,RATOLV,SCALRAT,VIRRAT(3),IRATTLE,JRATTLE,
     *                NRATTLE,MXRATT,LFFOLDCORD,LFFLSTRAT,LFFDSTRAT,
     *                LFFVELSV,IRATQM
      COMMON /FFRDF / DELRDF,NUMGRD,NUMSUM,NRDF,LFFGOFR,
     *                LFFNFRAG1,LFFNFRAG2,LFFFRAG1,LFFFRAG2,
     *                NRDEN,NBINRDEN,LFFPRO,LFFNRDPRATM,
     *                LFFRDPRATM,LFFDIESTEP,LFFDI1STEP
      COMMON /FFRMSD/ DIPT(8),TIMGYRA,TIMRALL,
     *                NATPDB,NGYRA,NDIEL,NRALL,
     *                LFFRALL0,NRMSD,LFFRMSD0,KFREEAB(201),
     *                NRIJMM,IJRMM(2,100),NRIJQM,IJRQM(2,100),
     *                NAIJKMM,IJKMM(3,100),NAIJKQM,IJKQM(3,100),
     *                NFIXMM,IFIXMM(200),NFIXQM,IFIXQM(200)
      COMMON /FFSPH / SPHRAD,SPHSIG,SPHEPS,IADDWAT
      COMMON /FFUMBR/ UMBFC,UMBR0,UMBSIZE,
     *                NUMBBIN,NUMBATM(6),NUMBTYP,LFFUMBHIS,
     *                UM2FC,UM2R0,UM2SIZE,
     *                NUM2BIN,NUM2ATM(6),NUM2TYP,LFFUM2HIS
      COMMON /FFVIB / JVIBOUT,NVIBMM,LFFDIPSTEP,LFFVELSTEP,LFFIVIBMM,
     *                LFFDQMSTEP,LFFVQMSTEP,LFFDM1STEP,LFFVM1STEP,
     *                LFFDMMSTEP,LFFVMMSTEP,LFFQMVSTEP,LFFMMVSTEP
      COMMON /FMCOM / X(1)
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /INFOTD/ CNVTOL,PFREQ(2),MODTD,
     *                JANST,NRADT,NTHET,NPHIT,NLEBT,
     *                NSTAT,NTRIAL,MAXVEC,NTHST,IRECTD,ITDFG,ITDPRP,
     *                TRIPLET,SG1T,NONEQR,TAMMD,TPA,ALPHKWD,BETAKWD,
     *                SPCP(3),MULTD,MREKT,MRDEA,MTHST,IFEDAT(4)
      COMMON /IOFILE/ IR,IW,IP,IJK,IJKT,IDAF,NAV,IODA(950)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /WFNOPT/ SCFTYP,VBTYP,DFTYPE,TDDFTYP,CITYP,CCTYP,
     *                MPLEVL,MPCTYP
C
      DATA RNONE/8HNONE    /
C
C     NANDUN THELLAMUREGE, HUI LI, MAR 2011, LINCOLN
C     FENGCHAO CUI, HUI LI, JUN 2012 (ADD FREE ENERGY)
C
      NNNATM = NFFAT + NAT - LISTQM(NFFAT+NAT+1)
      IF(NNNATM.GT.2)     NDOF = 3*NNNATM - 6
      IF(NNNATM.EQ.2)     NDOF = 3*NNNATM - 5
      IF(NRATTLE.GT.0)    NDOF = NDOF - NRATTLE
      IF(XBOX.LT.1.0D+30) NDOF = NDOF + 1
      IF(YBOX.LT.1.0D+30) NDOF = NDOF + 1
      IF(ZBOX.LT.1.0D+30) NDOF = NDOF + 1
      IF(NNNATM.EQ.1)     NDOF = 3
C
      IF(ISTEP.EQ.0.AND.MASWRK) THEN
         WRITE(IW,*)' '
         WRITE(IW,*)'======================= QUANPOL MD SIMULATION',
     *              ' INITIATED ======================='
         WRITE(IW,*)' '
      END IF
C
      IF(IFIXSOL.GT.0) THEN
         IF(MASWRK.AND.
     *      (MOD(ISTEP,JOUT).EQ.0.OR.ISTEP.EQ.NSTEP))
     *      WRITE(IW,'(A,F20.10,A,A,I8)')
     *      ' FIXSOL TOTAL SURFACE AREA =',FIXA,' A**2,',
     *      ' NFFTS=',NFFTS
         IF(MASWRK.AND.NCYCLE.LT.200.AND.
     *      (MOD(ISTEP,JOUT).EQ.0.OR.ISTEP.EQ.NSTEP))
     *      WRITE(IW,'(A,I3,A,F11.6)')
     *      ' FIXSOL CONVERGED IN ',NCYCLE,
     *      ' ITERATIONS, TOTAL SURFACE CHARGE=',FIXQ
         IF(MASWRK.AND.NCYCLE.EQ.200)
     *       WRITE(IW,'(A,I3,A,F10.6,A,F12.10)')
     *      ' FIXSOL NOT CONVERGED IN ',NCYCLE,
     *      ' ITERATIONS.  TOTAL SURFACE CHARGE=',FIXQ
      END IF
C
      AMU2QM = ZERO
      DO IAT =1, NAT
         DUMY = ZERO
         DO III = 1, 3
            DUMY = DUMY + QMVEL(III,IAT)*QMVEL(III,IAT)
         ENDDO
         AMU2QM = AMU2QM + DUMY*QMZMAS(IAT)
      ENDDO
C
      AMU2MM = ZERO
      DO IFFAT = 1, NFFAT
         IF(LISTQM(IFFAT).EQ.0.OR.LISTQM(IFFAT).GT.NAT) THEN
            DUMY = ZERO
            DO III = 1, 3
               DUMY = DUMY + VEL(III,IFFAT)*VEL(III,IFFAT)
            ENDDO
            AMU2MM = AMU2MM + DUMY*ZMAS(IFFAT)
         END IF
      ENDDO
C
      AMU2 = AMU2QM + AMU2MM
C
      VOL   = MIN(XBOX*YBOX*ZBOX, FOURTHIRD*PI*SPHRAD**3)
      ENKIN = PT5*AMU2
      TEMP  = TOKELVIN*AMU2/NDOF
      IF(VOL.LT.1.0D+30) THEN
         PRESX = (AMU2*ONETHIRD - (VIR(1)+VIRRAT(1)))/VOL*TOBAR
         PRESY = (AMU2*ONETHIRD - (VIR(2)+VIRRAT(2)))/VOL*TOBAR
         PRESZ = (AMU2*ONETHIRD - (VIR(3)+VIRRAT(3)))/VOL*TOBAR
         PRES  = (PRESX+PRESY+PRESZ)*ONETHIRD
      ELSE
         PRESX = ZERO
         PRESY = ZERO
         PRESZ = ZERO
         PRES  = ZERO
      END IF
      IF(GOPARR) THEN
         CALL DDI_BCAST(465,'F',TEMP,1,MASTER)
         CALL DDI_BCAST(466,'F',PRES,1,MASTER)
         CALL DDI_BCAST(467,'F',PRESX,1,MASTER)
         CALL DDI_BCAST(468,'F',PRESY,1,MASTER)
         CALL DDI_BCAST(469,'F',PRESZ,1,MASTER)
      END IF
C
C     -- NOT EASY TO CORRECT THE PRESSURE FOR THE INITIAL STEP --
C        MAKING IT ZERO IS OK BECAUSE IT IS NOT AVERAGED
      IF(ISTEP.EQ.0 .AND. NRATTLE.GT.0) THEN
         PRESX= ZERO
         PRESY= ZERO
         PRESZ= ZERO
         PRES = ZERO
      END IF
C
C     -- FORM FREE ENERGY ITEMS --
      SOL1MM = ZERO
      SOL2MM = ZERO
      PMF1MM = ZERO
      POT1QMMM = ZERO
      POT2QMMM = ZERO
      IF(IFEPTYP.GT.0) THEN
         IF(IDOPOL.EQ.0.AND.IEWALD.EQ.0) THEN
C           -SOL1CH,SOL1LJ,PMF1CH,PMF1LJ,PMF1BD... ARE READY
            SOL1MM = SOL1CH+SOL1LJ
            SOL2MM = SOL2CH+SOL2LJ
            PMF1MM = PMF1CH+PMF1LJ+
     *               PMF1BD+PMF1AG+PMF1DR+PMF1DB+PMF1WG+PMF1CM+
     *               SOL1CH+SOL1LJ
         END IF
         IF(IDOPOL.GT.0.OR.IEWALD.GT.0) THEN
            SOL1CH  = (WPERT1-WSIMUL)*(ENCHARB -ENCHARA)
            SOL2CH  = (WPERT2-WSIMUL)*(ENCHARB -ENCHARA)
            SOL1LJ  = (WPERT1-WSIMUL)*(ENLJRB  +ENLJDB -ENLJRA -ENLJDA)
            SOL2LJ  = (WPERT2-WSIMUL)*(ENLJRB  +ENLJDB -ENLJRA -ENLJDA)
            SOL1MM  = SOL1CH+SOL1LJ
            SOL2MM  = SOL2CH+SOL2LJ
            ENPOTA  = EN12A   + EN123A + EN123R4A + EN234W1A  + EN123B4A
     *              + ENCHARA + ENPOLA + ENRXNA   + ENRXNPOLA + ENRXNRA
     *              + ENLJRA  + ENLJDA + ENCMAPA  + ENFIXSOA
     *              + ENBIAS  + EN12312A
            ENPOTB  = EN12B   + EN123B + EN123R4B + EN234W1B  + EN123B4B
     *              + ENCHARB + ENPOLB + ENRXNB   + ENRXNPOLB + ENRXNRB
     *              + ENLJRB  + ENLJDB + ENCMAPB  + ENFIXSOB
     *              + ENBIAS  + EN12312B
            PMF1PO  = (WPERT1-WSIMUL)*(ENPOLB  -ENPOLA)
            PMF1MM  = (WPERT1-WSIMUL)*(ENPOTB - ENPOTA)
         END IF
         IF(NAT.GT.0) THEN
            POT1QMMM  = (WPERT1-WSIMUL)*(ETOTQB - ETOTQA)
            POT2QMMM  = (WPERT2-WSIMUL)*(ETOTQB - ETOTQA)
            IF(MPLEVL.EQ.2)THEN
               POT1QMMM  = (WPERT1-WSIMUL)*(EMP2QB - EMP2QA)
               POT2QMMM  = (WPERT2-WSIMUL)*(EMP2QB - EMP2QA)
            END IF
            IF(TDDFTYP.NE.RNONE) THEN
               POT1QMMM  = (WPERT1-WSIMUL)*(ESTATEQB - ESTATEQA)
               POT2QMMM  = (WPERT2-WSIMUL)*(ESTATEQB - ESTATEQA)
            END IF
         END IF
      END IF
C
      IF(ISTEP.GT.0) THEN
         PMEANX     = (PRESX     +PMEANX*(ISTEP-1))/ISTEP
         PMEANY     = (PRESY     +PMEANY*(ISTEP-1))/ISTEP
         PMEANZ     = (PRESZ     +PMEANZ*(ISTEP-1))/ISTEP
         PMEAN      = (PRES      +PMEAN *(ISTEP-1))/ISTEP
         VOLAV      = (VOL       +VOLAV *(ISTEP-1))/ISTEP
         ENKAV      = (ENKIN     +ENKAV *(ISTEP-1))/ISTEP
         TEMPAV     = (TEMP      +TEMPAV*(ISTEP-1))/ISTEP
         ASOL1CH    = (SOL1CH    +ASOL1CH*(ISTEP-1))/ISTEP
         ASOL2CH    = (SOL2CH    +ASOL2CH*(ISTEP-1))/ISTEP
         ASOL1LJ    = (SOL1LJ    +ASOL1LJ*(ISTEP-1))/ISTEP
         ASOL2LJ    = (SOL2LJ    +ASOL2LJ*(ISTEP-1))/ISTEP
         ASOL1MM    = (SOL1MM    +ASOL1MM*(ISTEP-1))/ISTEP
         ASOL2MM    = (SOL2MM    +ASOL2MM*(ISTEP-1))/ISTEP
         SOLFRE1MMX = (EXP(-SOL1MM/(BOLTZK*TEMP0)) +
     *                 EXP(-SOLFRE1MM/(BOLTZK*TEMP0))*(ISTEP-1))/ISTEP
         SOLFRE1MM  = -LOG(SOLFRE1MMX)*BOLTZK*TEMP0
         SOLFRE2MMX = (EXP(-SOL2MM/(BOLTZK*TEMP0)) +
     *                 EXP(-SOLFRE2MM/(BOLTZK*TEMP0))*(ISTEP-1))/ISTEP
         SOLFRE2MM  = -LOG(SOLFRE2MMX)*BOLTZK*TEMP0
         PMFFRE1MMX = (EXP(-PMF1MM/(BOLTZK*TEMP0)) +
     *                 EXP(-PMFFRE1MM/(BOLTZK*TEMP0))*(ISTEP-1))/ISTEP
         PMFFRE1MM  = -LOG(PMFFRE1MMX)*BOLTZK*TEMP0
         IF(NAT.GT.0) THEN
         XFRE1QMMM  = (EXP(-POT1QMMM/(BOLTZK*TEMP0))
     *                +EXP(-FRE1QMMM/(BOLTZK*TEMP0))*(ISTEP-1))/ISTEP
          FRE1QMMM  = -LOG(XFRE1QMMM)*BOLTZK*TEMP0
         XFRE2QMMM  = (EXP(-POT2QMMM/(BOLTZK*TEMP0))
     *                +EXP(-FRE2QMMM/(BOLTZK*TEMP0))*(ISTEP-1))/ISTEP
          FRE2QMMM  = -LOG(XFRE2QMMM)*BOLTZK*TEMP0
         END IF
      END IF
C
C     -- FFMD1 --
      IF(NAT.LE.0) THEN
         ENPOT = EN12    + EN123  + EN123R4 + EN234W1  + EN123B4
     *         + ENCHAR  + ENPOL  + ENRXN   + ENRXNPOL + ENRXNR
     *         + ENLJR   + ENLJD  + ENCMAP  + ENFIXSO
     *         + ENBIAS  + EN12312
         ENTOT = ENPOT  + ENKIN
         IF(ISTEP.GT.0) THEN
            ENPAV = (ENPOT+(ENPAV*(ISTEP-1)))/ISTEP
         END IF
         ENRXN = ENRXN + ENRXNPOL
C
         IF(MASWRK.AND.(MOD(ISTEP,JOUT).EQ.0.OR.ISTEP.EQ.NSTEP)) THEN
            WRITE(IW,*)' '
            WRITE(IW,'(1X,A,I10,33X,A,F19.2,A)')'MD STEP ',ISTEP,
     *        'TIME= ',ISTEP*DT*2.418884326505D-02,' FS'
            WRITE(IW,9000)
     *      'BOND STRETCHING               ENERGY = ', EN12   *TOKCAL
            WRITE(IW,9000)
     *      'BOND ANGLE BENDING            ENERGY = ', EN123  *TOKCAL
            WRITE(IW,9000)
     *      'STRETCHING BENDING            ENERGY = ', EN12312*TOKCAL
            WRITE(IW,9000)
     *      'DIHEDRAL ROTATION             ENERGY = ', EN123R4*TOKCAL
            WRITE(IW,9000)
     *      'DIHEDRAL BENDING              ENERGY = ', EN123B4*TOKCAL
            WRITE(IW,9000)
     *      'CMAP                          ENERGY = ', ENCMAP *TOKCAL
            WRITE(IW,9000)
     *      'WAGGING                       ENERGY = ', EN234W1*TOKCAL
            IF(NUMBTYP.GT.0)
     *      WRITE(IW,9000)
     *      'UMBRELLA SAMPLING BIAS        ENERGY = ', ENBIAS *TOKCAL
            WRITE(IW,9000)
     *      'LJ REPULSION                  ENERGY = ', ENLJR  *TOKCAL
            WRITE(IW,9000)
     *      'LJ DISPERSION                 ENERGY = ', ENLJD  *TOKCAL
            WRITE(IW,9000)
     *      'CHARGE                        ENERGY = ', ENCHAR *TOKCAL
            WRITE(IW,9000)
     *      'INDUCED DIPOLE                ENERGY = ', ENPOL  *TOKCAL
            WRITE(IW,9000)
     *      'SPHSOL                        ENERGY = ', ENRXN  *TOKCAL
            WRITE(IW,9000)
     *      'FIXSOL                        ENERGY = ', ENFIXSO*TOKCAL
            WRITE(IW,9000)
     *      'QM CENTER                     ENERGY = ', ENCENT *TOKCAL
            WRITE(IW,9000)
     *      'SPHERE                        ENERGY = ', ENRXNR *TOKCAL
            WRITE(IW,9000)
     *      'POTENTIAL                     ENERGY = ', ENPOT  *TOKCAL
            WRITE(IW,9000)
     *      'KINETIC                       ENERGY = ', ENKIN  *TOKCAL
            WRITE(IW,9000)
     *      'TOTAL                         ENERGY = ', ENTOT  *TOKCAL
            WRITE(IW,9001)
     *      'TEMPERATURE                          = ', TEMP
            WRITE(IW,9002)
     *      'PRESSURE                             = ', PRES
            IF(IPSTAT.EQ.3) THEN
            WRITE(IW,9002)
     *      'PRESSUREX                            = ', PRESX
            WRITE(IW,9002)
     *      'PRESSUREY                            = ', PRESY
            WRITE(IW,9002)
     *      'PRESSUREZ                            = ', PRESZ
            END IF
            IF(VOL.LE.1.0D+30)
     *      WRITE(IW,9003)
     *      'VOLUME                               = ', VOL*TOANGS**3
            IF(VOL.GT.1.0D+30)
     *      WRITE(IW,9004)
     *      'VOLUME                               = '
C
            IF(ISTEP.GT.0) THEN
              WRITE(IW,*)' '
              WRITE(IW,9000)
     *        'AVERAGE POTENTIAL             ENERGY = ', ENPAV*TOKCAL
              WRITE(IW,9000)
     *        'AVERAGE KINETIC               ENERGY = ', ENKAV*TOKCAL
              WRITE(IW,9000)
     *        'AVERAGE TOTAL                 ENERGY = ',(ENPAV+ENKAV)
     *                                                        *TOKCAL
              WRITE(IW,9001)
     *        'AVERAGE TEMPERATURE                  = ', TEMPAV
              WRITE(IW,9002)
     *        'AVERAGE PRESSURE                     = ', PMEAN
              IF(IPSTAT.EQ.3) THEN
              WRITE(IW,9002)
     *        'AVERAGE PRESSUREX                    = ', PMEANX
              WRITE(IW,9002)
     *        'AVERAGE PRESSUREY                    = ', PMEANY
              WRITE(IW,9002)
     *        'AVERAGE PRESSUREZ                    = ', PMEANZ
              END IF
              IF(VOL.LE.1.0D+30)
     *        WRITE(IW,9003)
     *        'AVERAGE VOLUME                       = ',VOLAV*TOANGS**3
            END IF
C
C           -- MM SOL1 IN FFMD1 --
            IF(IFEPTYP.EQ.1) THEN
              WRITE(IW,*)' '
              WRITE(IW,'(A,F9.6,A,F9.6)')
     *        ' RELATIVE SOLVATION ENERGY: FROM    WSIMUL= ',WSIMUL,
     *        '    TO    WPERT1= ',WPERT1
              WRITE(IW,9000)
     *        'RELATIVE CHARGE SOLVATION     ENERGY = ', SOL1CH*TOKCAL
              WRITE(IW,9000)
     *        'RELATIVE LJ SOLVATION         ENERGY = ', SOL1LJ*TOKCAL
              WRITE(IW,9000)
     *        'RELATIVE MM SOLVATION         ENERGY = ', SOL1MM*TOKCAL
              IF(ISTEP.GT.0) THEN
              WRITE(IW,9000)
     *        'AVERAGE REL CHARGE SOLVATION  ENERGY = ', ASOL1CH*TOKCAL
              WRITE(IW,9000)
     *        'AVERAGE REL LJ SOLVATION      ENERGY = ', ASOL1LJ*TOKCAL
              WRITE(IW,9000)
     *        'AVERAGE REL MM SOLVATION      ENERGY = ', ASOL1MM*TOKCAL
              WRITE(IW,9000)
     *        'RELATIVE SOLVATION FREE       ENERGY = ',
     *                                                 SOLFRE1MM*TOKCAL
              END IF
            END IF
C
C           -- MM SOL2 IN FFMD1 --
            IF(IFEPTYP.EQ.1.AND.WPERT2.NE.WPERT1) THEN
              WRITE(IW,*)' '
              WRITE(IW,'(A,F9.6,A,F9.6)')
     *        ' RELATIVE SOLVATION ENERGY: FROM    WSIMUL= ',WSIMUL,
     *             '    TO    WPERT2= ',WPERT2
              WRITE(IW,9000)
     *        'RELATIVE CHARGE SOLVATION     ENERGY = ', SOL2CH*TOKCAL
              WRITE(IW,9000)
     *        'RELATIVE LJ SOLVATION         ENERGY = ', SOL2LJ*TOKCAL
              WRITE(IW,9000)
     *        'RELATIVE MM SOLVATION         ENERGY = ', SOL2MM*TOKCAL
              IF(ISTEP.GT.0) THEN
              WRITE(IW,9000)
     *        'AVERAGE REL CHARGE SOLVATION  ENERGY = ', ASOL2CH*TOKCAL
              WRITE(IW,9000)
     *        'AVERAGE REL LJ SOLVATION      ENERGY = ', ASOL2LJ*TOKCAL
              WRITE(IW,9000)
     *        'AVERAGE REL MM SOLVATION      ENERGY = ', ASOL2MM*TOKCAL
              WRITE(IW,9000)
     *        'RELATIVE SOLVATION FREE       ENERGY = ',
     *                                                 SOLFRE2MM*TOKCAL
              END IF
            END IF
C
C           -- MM PMF IN FFMD1 --
            IF(IFEPTYP.EQ.2) THEN
              WRITE(IW,*)' '
              WRITE(IW,*)'RELATIVE ENERGY: FROM $FFDATA TO $FFDATB',
     *                   '              (POTENTIAL OF MEAN FORCE)'
              IF(IDOPOL.EQ.0.AND.IEWALD.EQ.0) THEN
              WRITE(IW,9000)
     *        'RELATIVE BOND STRETCHING  PMF ENERGY = ', PMF1BD*TOKCAL
              WRITE(IW,9000)
     *        'RELATIVE BOND ANGLE BEND  PMF ENERGY = ', PMF1AG*TOKCAL
              WRITE(IW,9000)
     *        'RELATIVE DIHEDRAL ROTAT   PMF ENERGY = ', PMF1DR*TOKCAL
              WRITE(IW,9000)
     *        'RELATIVE DIHEDRAL BENDING PMF ENERGY = ', PMF1DB*TOKCAL
              WRITE(IW,9000)
     *        'RELATIVE CMAP             PMF ENERGY = ', PMF1CM*TOKCAL
              WRITE(IW,9000)
     *        'RELATIVE WAGGING          PMF ENERGY = ', PMF1WG*TOKCAL
              WRITE(IW,9000)
     *        'RELATIVE CHARGE           PMF ENERGY = ', PMF1CH*TOKCAL
              WRITE(IW,9000)
     *        'RELATIVE LJ               PMF ENERGY = ', PMF1LJ*TOKCAL
              WRITE(IW,9000)
     *        'RELATIVE CHARGE           SOL ENERGY = ', SOL1CH*TOKCAL
              WRITE(IW,9000)
     *        'RELATIVE LJ               SOL ENERGY = ', SOL1LJ*TOKCAL
              END IF
              WRITE(IW,9000)
     *        'RELATIVE POTENTIAL            ENERGY = ', PMF1MM*TOKCAL
              IF(ISTEP.GT.0)
     *        WRITE(IW,9000)
     *        'RELATIVE FREE                 ENERGY = ',
     *                                                PMFFRE1MM*TOKCAL
            END IF
            WRITE(IW,*)' '
            CALL TIMIT(1)
            WRITE(IW,*)' '
         END IF
C
C     -- FFMD2 --
      ELSE IF(NAT.GT.0.AND.NFFAT.GT.0) THEN
                              ENPOT = ETOT
         IF(MPLEVL.EQ.2)      ENPOT = EMP2
         IF(TDDFTYP.NE.RNONE) ENPOT = ESTATE(NTHST)
         ENPOT2= ENFFQM2
         ENTOT = ENPOT  + ENKIN
         IF(ISTEP.GT.0) THEN
            ENPAV  = (ENPOT +(ENPAV *(ISTEP-1)))/ISTEP
            ENPAV2 = (ENPOT2+(ENPAV2*(ISTEP-1)))/ISTEP
         ELSE
            ENPAV  = ENPOT
            ENPAV2 = ENPOT2
         END IF
         ENRXN = ENRXN + ENRXNPOL
C
         IF(MASWRK.AND.(MOD(ISTEP,JOUT).EQ.0.OR.ISTEP.EQ.NSTEP)) THEN
            WRITE(IW,*)' '
            WRITE(IW,'(1X,A,I10,33X,A,F19.2,A)')'MD STEP ',ISTEP,
     *        'TIME= ',ISTEP*DT*2.418884326505D-02,' FS'
            WRITE(IW,9000)
     *      'BOND STRETCHING               ENERGY = ', EN12   *TOKCAL
            WRITE(IW,9000)
     *      'BOND ANGLE BENDING            ENERGY = ', EN123  *TOKCAL
            WRITE(IW,9000)
     *      'STRETCHING BENDING            ENERGY = ', EN12312*TOKCAL
            WRITE(IW,9000)
     *      'DIHEDRAL ROTATION             ENERGY = ', EN123R4*TOKCAL
            WRITE(IW,9000)
     *      'DIHEDRAL BENDING              ENERGY = ', EN123B4*TOKCAL
            WRITE(IW,9000)
     *      'CMAP                          ENERGY = ', ENCMAP *TOKCAL
            WRITE(IW,9000)
     *      'WAGGING                       ENERGY = ', EN234W1*TOKCAL
            IF(NUMBTYP.GT.0)
     *      WRITE(IW,9000)
     *      'UMBRELLA SAMPLING BIAS        ENERGY = ', ENBIAS *TOKCAL
            WRITE(IW,9000)
     *      'LJ REPULSION                  ENERGY = ', ENLJR  *TOKCAL
            WRITE(IW,9000)
     *      'LJ DISPERSION                 ENERGY = ', ENLJD  *TOKCAL
            WRITE(IW,9000)
     *      'CHARGE                        ENERGY = ', ENCHAR *TOKCAL
            WRITE(IW,9000)
     *      'INDUCED DIPOLE                ENERGY = ', ENPOL  *TOKCAL
            WRITE(IW,9000)
     *      'SPHSOL                        ENERGY = ', ENRXN  *TOKCAL
            WRITE(IW,9000)
     *      'FIXSOL                        ENERGY = ', ENFIXSO*TOKCAL
            WRITE(IW,9000)
     *      'QM CENTER                     ENERGY = ', ENCENT *TOKCAL
            WRITE(IW,9000)
     *      'SPHERE                        ENERGY = ', ENRXNR *TOKCAL
            WRITE(IW,9000)
     *      'POTENTIAL                     ENERGY = ', ENPOT  *TOKCAL
            IF(IDOQM2.EQ.1)
     *      WRITE(IW,9000)
     *      '2ND POTENTIAL                 ENERGY = ', ENPOT2 *TOKCAL
            WRITE(IW,9000)
     *      'KINETIC                       ENERGY = ', ENKIN  *TOKCAL
            WRITE(IW,9000)
     *      'TOTAL                         ENERGY = ', ENTOT  *TOKCAL
            WRITE(IW,9001)
     *      'TEMPERATURE                          = ', TEMP
            WRITE(IW,9002)
     *      'PRESSURE                             = ', PRES
            IF(IPSTAT.EQ.3) THEN
            WRITE(IW,9002)
     *      'PRESSUREX                            = ', PRESX
            WRITE(IW,9002)
     *      'PRESSUREY                            = ', PRESY
            WRITE(IW,9002)
     *      'PRESSUREZ                            = ', PRESZ
            END IF
            IF(VOL.LE.1.0D+30)
     *      WRITE(IW,9003)
     *      'VOLUME                               = ', VOL*TOANGS**3
            IF(VOL.GT.1.0D+30)
     *      WRITE(IW,9004)
     *      'VOLUME                               = '
C
            IF(ISTEP.GT.0) THEN
              WRITE(IW,*)' '
              WRITE(IW,9000)
     *        'AVERAGE POTENTIAL             ENERGY = ', ENPAV *TOKCAL
              IF(IDOQM2.EQ.1)
     *        WRITE(IW,9000)
     *        'AVERAGE 2ND POTENTIAL         ENERGY = ', ENPAV2*TOKCAL
              WRITE(IW,9000)
     *        'AVERAGE KINETIC               ENERGY = ', ENKAV *TOKCAL
              WRITE(IW,9000)
     *        'AVERAGE TOTAL                 ENERGY = ' ,(ENPAV+ENKAV)
     *                                                         *TOKCAL
              WRITE(IW,9001)
     *        'AVERAGE TEMPERATURE                  = ', TEMPAV
              WRITE(IW,9002)
     *        'AVERAGE PRESSURE                     = ', PMEAN
              IF(IPSTAT.EQ.3) THEN
              WRITE(IW,9002)
     *        'AVERAGE PRESSUREX                    = ', PMEANX
              WRITE(IW,9002)
     *        'AVERAGE PRESSUREY                    = ', PMEANY
              WRITE(IW,9002)
     *        'AVERAGE PRESSUREZ                    = ', PMEANZ
              END IF
              IF(VOL.LE.1.0D+30)
     *        WRITE(IW,9003)
     *        'AVERAGE VOLUME                       = ',VOLAV*TOANGS**3
            END IF
C
            IF(IFEPTYP.EQ.1) THEN
              WRITE(IW,*)' '
              WRITE(IW,'(A,F9.6,A,F9.6)')
     *        ' RELATIVE FREE ENERGY:      FROM    WSIMUL= ',WSIMUL,
     *             '    TO    WPERT1= ',WPERT1
              WRITE(IW,9000)
     *        'RELATIVE QMMM POTENTIAL       ENERGY = ',POT1QMMM*TOKCAL
              IF(ISTEP.GT.0)
     *        WRITE(IW,9000)
     *        'RELATIVE QMMM FREE            ENERGY = ',
     *                                                  FRE1QMMM*TOKCAL
            END IF
            IF(IFEPTYP.EQ.1.AND.WPERT2.NE.WPERT1) THEN
              WRITE(IW,*)' '
              WRITE(IW,'(A,F9.6,A,F9.6)')
     *        ' RELATIVE FREE ENERGY:      FROM    WSIMUL= ',WSIMUL,
     *             '    TO    WPERT2= ',WPERT2
              WRITE(IW,9000)
     *        'RELATIVE QMMM POTENTIAL       ENERGY = ',POT2QMMM*TOKCAL
              IF(ISTEP.GT.0)
     *        WRITE(IW,9000)
     *        'RELATIVE QMMM FREE            ENERGY = ',
     *                                                  FRE2QMMM*TOKCAL
            END IF
C
            IF(IFEPTYP.EQ.2) THEN
              WRITE(IW,*)' '
              WRITE(IW,*)'RELATIVE ENERGY: FROM $FFDATA TO $FFDATB',
     *                   '              (POTENTIAL OF MEAN FORCE)'
              WRITE(IW,9000)
     *        'RELATIVE QMMM POTENTIAL       ENERGY = ',POT1QMMM*TOKCAL
              IF(ISTEP.GT.0)
     *        WRITE(IW,9000)
     *        'RELATIVE QMMM FREE            ENERGY = ',
     *                                                  FRE1QMMM*TOKCAL
            END IF
C
            WRITE(IW,*)' '
            CALL TIMIT(1)
            WRITE(IW,*)' '
         END IF
C        - ALWAYS PRINT OUT TDDFT RESULTS -
         IF(MASWRK.AND.TDDFTYP.NE.RNONE)THEN
            WRITE(IW,'(1X,A,9(1X,I6))')
     *      'TDDFT STATES:',
     *      (I,I=0,MIN(8,NSTAT))
            WRITE(IW,'(1X,A,9(1X,F6.3))')
     *      'ENERGY IN EV:',
     *      0.0, ((ESTATE(I)-ESCF)*TOEV,I=1,MIN(8,NSTAT))
            WRITE(IW,*)' '
         END IF
      END IF
C
C     -- ALWAYS PRINT OUT DIPOLE AND VEL --
C
      IF(ISTEP.LE.100000) THEN
         CALL VELDIP(CORD,ZMAS,X(LFFCHARG),X(LFFDIP),ISTEP,
     *            X(LFFDIPSTEP),X(LFFVELSTEP),VEL,X(LFFIVIBMM),
     *            X(LFFDQMSTEP),X(LFFDMMSTEP),X(LFFDM1STEP),
     *            X(LFFVQMSTEP),X(LFFVMMSTEP),X(LFFVM1STEP),
     *            QMVEL,X(LFFQMVSTEP),X(LFFMMVSTEP),X(LFFLISTQM))
      END IF
C
C     -- CALCULATE DIELECTRIC CONSTANT --
C
      IF(NDIEL.GT.0) THEN
         CALL DIELECT(CORD,ZMAS,X(LFFCHARG),X(LFFDIP),ISTEP,
     *                X(LFFDIESTEP),X(LFFDI1STEP))
      END IF
C
C     -- CALCULATE RDF EVERY STEP --
C
      IF(NRDF.GT.0) THEN
         IF(ISTEP.GE.0) THEN
            IF(SPHRAD.GE.1.0D+30) THEN
               CALL RDF(CORD,ISTEP,X(LFFATMNAM),X(LFFFRAG1),
     *                  X(LFFFRAG2),X(LFFNFRAG1),X(LFFNFRAG2),
     *                  X(LFFGOFR))
            ELSE
               CALL RDFSPH(CORD,ISTEP,X(LFFATMNAM),X(LFFFRAG1),
     *                     X(LFFFRAG2),X(LFFNFRAG1),X(LFFNFRAG2),
     *                     X(LFFGOFR))
            END IF
         END IF
      END IF
C
C     -- CALCULATE RDEN EVERY STEP --
C
      IF(NRDEN.GT.0) THEN
         IF(ISTEP.GE.0) THEN
            CALL RDEN(CORD,ISTEP,X(LFFATMNAM),X(LFFRDPRATM),
     *                X(LFFNRDPRATM),X(LFFPRO))
         END IF
      END IF
C
C     -- CALCULATE DIFFUSION COEFFICIENT --
C
      IF(NDFS.GT.0) THEN
         IF(MOD(ISTEP,INT(TIMDFS/DT+0.1D+00)).EQ.0) THEN
            CALL DFS(CORD,ISTEP,X(LFFATMNAM),X(LFFDFSC0),
     *               X(LFFDFSA),X(LFFDFSC),X(LFFDFSCAV),X(LFFDFSN))
         END IF
      END IF
C
C     -- CALCULATE RMSD EVERY JOUT STEPS --
C
      IF(NRMSD.EQ.1) THEN
         IF(MOD(ISTEP,JOUT).EQ.0) THEN
            CALL RMSD(CORD,ISTEP,X(LFFZANF),X(LFFRMSD0))
         END IF
      END IF
C
C     -- CALCULATE RADIUS OF GYRATION --
C
      IF(NGYRA.GT.0) THEN
         IF(MOD(ISTEP,INT(TIMGYRA/DT+0.1D+00)).EQ.0) THEN
            CALL GYRA(CORD,X(LFFZANF),ZMAS)
         END IF
      END IF
C
C     -- CALCULATE RALL --
C
      IF(NRALL.EQ.1) THEN
         IF(MOD(ISTEP,INT(TIMRALL/DT+0.1D+00)).EQ.0) THEN
            CALL RALL(X(LFFRALL0),CORD,ISTEP)
         END IF
      END IF
C
C     -- CALCULATE SELECT DISTANCES EVERY JOUT STEPS --
C
      IF(((NRIJMM+NRIJQM+NAIJKMM+NAIJKQM).GT.0) .AND.
     *    MOD(ISTEP,JOUT).EQ.0) CALL DISIJ(CORD)
C
C     -- PRINT OUT UMBRELLA SAMPLING HISTOGRAM AND PMF --
C
      IF(MASWRK.AND.NUMBTYP.GT.0.AND.NUM2TYP.EQ.0) THEN
      IF((ISTEP.GT.0.AND.MOD(ISTEP,JOUT).EQ.0).OR.ISTEP.EQ.NSTEP) THEN
         UMBRLOW = UMBR0 - (DBLE(NUMBBIN-1)/2)*UMBSIZE
         MIDPTI  = (NUMBBIN+1)/2
         IF(NUMBTYP.EQ.12.OR.NUMBTYP.EQ.1212) THEN
            WRITE(IW,'(/1X,A,I10)')
     *      '1D UMBRELLA SAMPLING HISTOGRAM AT MD STEP=',ISTEP
            WRITE(IW,'(1X,A,F10.4,1X,A)')
     *      'K=',UMBFC*TOBOHR*TOBOHR*TOKCAL,'KCAL MOL-1 A-2'
            WRITE(IW,'(1X,A,F10.4,1X,A)') 'T=',TEMP0,'KELVIN'
            WRITE(IW,*)
     *      '  IBIN       R(A)    HISTOGRAM       PMF KCAL/MOL'
            RUNITI= TOANGS
         ELSE IF(NUMBTYP.EQ.123.OR.NUMBTYP.EQ.1234) THEN
            WRITE(IW,'(/1X,A,I10)')
     *      '1D UMBRELLA SAMPLING HISTOGRAM AT MD STEP=',ISTEP
            WRITE(IW,'(1X,A,F10.4,1X,A)')
     *      'K=',UMBFC*TORAD*TORAD*TOKCAL,'KCAL MOL-1 DEG-2'
            WRITE(IW,'(1X,A,F10.4,1X,A)') 'T=',TEMP0,'KELVIN'
            WRITE(IW,*)
     *      '  IBIN     R(DEG)    HISTOGRAM       PMF KCAL/MOL'
            RUNITI= DEGREE
         END IF
C
         FACTOR  = ONE/UMBHIS(MIDPTI)
         DO IBIN=1,NUMBBIN
            COUNTS = UMBHIS(IBIN)
            PMF    = ZERO
            IF(COUNTS.GT.ZERO.AND.UMBHIS(MIDPTI).GT.ZERO)
     *      PMF    = -BOLTZK*TEMP0*LOG(COUNTS*FACTOR)
     *               -PT5*UMBFC*((IBIN-MIDPTI)*UMBSIZE)
     *                         *((IBIN-MIDPTI)*UMBSIZE)
            WRITE(IW,'(1X,I6,1X,F10.4,1X,I12,1X,F18.10)')
     *      IBIN-MIDPTI,(UMBRLOW+UMBSIZE*(IBIN-1))*RUNITI,
     *      NINT(COUNTS),PMF*TOKCAL
         ENDDO
         WRITE(IW,*)' '
      END IF
      END IF
C
      IF(MASWRK.AND.NUMBTYP.GT.0.AND.NUM2TYP.GT.0) THEN
      IF((ISTEP.GT.0.AND.MOD(ISTEP,JOUT).EQ.0).OR.ISTEP.EQ.NSTEP) THEN
         WRITE(IW,'(/1X,A,A/)')
     *   '2D UMBRELLA SAMPLING HISTOGRAM WRITTEN TO THE ',
     *   'TRAJECTORY FILE EVERY KOUT STEPS.'
      END IF
      IF((ISTEP.GT.0.AND.MOD(ISTEP,KOUT).EQ.0).OR.ISTEP.EQ.NSTEP) THEN
         UMBRLOW = UMBR0 - (DBLE(NUMBBIN-1)/2)*UMBSIZE
         MIDPTI  = (NUMBBIN+1)/2
         IF(NUMBTYP.EQ.12.OR.NUMBTYP.EQ.1212) THEN
            WRITE(36,'(/1X,A,I10)')
     *      '2D UMBRELLA SAMPLING HISTOGRAM AT MD STEP=',ISTEP
            WRITE(36,'(1X,A,F10.4,1X,A)')
     *      'K=',UMBFC*TOBOHR*TOBOHR*TOKCAL,'KCAL MOL-1 A-2'
            RUNITI= TOANGS
         ELSE IF(NUMBTYP.EQ.123.OR.NUMBTYP.EQ.1234) THEN
            WRITE(36,'(/1X,A,I10)')
     *      '2D UMBRELLA SAMPLING HISTOGRAM AT MD STEP=',ISTEP
            WRITE(36,'(1X,A,F10.4,1X,A)')
     *      'K=',UMBFC*TORAD*TORAD*TOKCAL,'KCAL MOL-1 DEG-2'
            RUNITI= DEGREE
         END IF
         UM2RLOW = UM2R0 - (DBLE(NUM2BIN-1)/2)*UM2SIZE
         MIDPTJ  = (NUM2BIN+1)/2
         IF(NUM2TYP.EQ.12.OR.NUM2TYP.EQ.1212) THEN
            WRITE(36,'(1X,A,F10.4,1X,A)')
     *      'K=',UM2FC*TOBOHR*TOBOHR*TOKCAL,'KCAL MOL-1 A-2'
            RUNITJ= TOANGS
         ELSE IF(NUM2TYP.EQ.123.OR.NUM2TYP.EQ.1234) THEN
            WRITE(36,'(1X,A,F10.4,1X,A)')
     *      'K=',UM2FC*TORAD*TORAD*TOKCAL,'KCAL MOL-1 DEG-2'
            RUNITJ= DEGREE
         END IF
         WRITE(36,'(1X,A,F10.4,1X,A)') 'T=',TEMP0,'KELVIN'
C
         FACTOR = ONE/UM2HIS(MIDPTJ,MIDPTI)
         DO IBIN=1,NUMBBIN
         IF(NUMBTYP.EQ.12.OR.NUMBTYP.EQ.1212) THEN
            WRITE(36,*)
     *      '  IBIN       R(A)    HISTOGRAM       PMF KCAL/MOL'
         ELSE IF(NUMBTYP.EQ.123.OR.NUMBTYP.EQ.1234) THEN
            WRITE(36,*)
     *      '  IBIN     R(DEG)    HISTOGRAM       PMF KCAL/MOL'
         END IF
         WRITE(36,'(1X,I6,1X,F10.4)')
     *   IBIN-MIDPTI,(UMBRLOW+UMBSIZE*(IBIN-1))*RUNITI
         IF(NUM2TYP.EQ.12.OR.NUM2TYP.EQ.1212) THEN
            WRITE(36,*)
     *      '  JBIN       R(A)'
            RUNITJ= TOANGS
         ELSE IF(NUM2TYP.EQ.123.OR.NUM2TYP.EQ.1234) THEN
            WRITE(36,*)
     *      '  JBIN     R(DEG)'
            RUNITJ= DEGREE
         END IF
         DO JBIN=1,NUM2BIN
            COUNTS = UM2HIS(JBIN,IBIN)
            PMF    = ZERO
            IF(COUNTS.GT.ZERO.AND.UM2HIS(MIDPTJ,MIDPTI).GT.ZERO)
     *      PMF    = -BOLTZK*TEMP0*LOG(COUNTS*FACTOR)
     *               -PT5*UMBFC*((IBIN-MIDPTI)*UMBSIZE)
     *                         *((IBIN-MIDPTI)*UMBSIZE)
     *               -PT5*UM2FC*((JBIN-MIDPTJ)*UM2SIZE)
     *                         *((JBIN-MIDPTJ)*UM2SIZE)
            WRITE(36,'(1X,I6,1X,F10.4,1X,I12,1X,F18.10)')
     *      JBIN-MIDPTJ,(UM2RLOW+UM2SIZE*(JBIN-1))*RUNITJ,
     *      NINT(COUNTS),PMF*TOKCAL
         ENDDO
         ENDDO
         WRITE(36,*)' '
      END IF
      END IF
C
 9000 FORMAT(1X,A,F30.10,2X,'KCAL/MOL')
 9001 FORMAT(1X,A,F30.10,2X,'K')
 9002 FORMAT(1X,A,F30.10,2X,'BAR')
 9003 FORMAT(1X,A,F30.10,2X,'A**3')
 9004 FORMAT(1X,A,30X,   2X,'OPEN SYSTEM')
C
      CALL FLSHBF(IW)
C
      RETURN
      END
C*MODULE QUANPOC  *DECK TPSTAT
!>
!> @brief    T and P scaling
!>
!> @author   Nandun Thellamurege, Hui Li
!>           - Mar 2011
!>
!> @details  scale T and P in MD simulation
!>
      SUBROUTINE TPSTAT(CORD,VEL,QMVEL,ZMAS,QMZMAS,ONEMAS,QM1MAS,
     *                  ISTEP,LISTQM,TEMP,PRES,
     *                  PRESX,PRESY,PRESZ,LSTPOL)
      use mx_limits, only: mxatm
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      LOGICAL GOPARR,DSKWRK,MASWRK
C
      PARAMETER (PI=3.14159265358979323846264338D+00)
      PARAMETER (TOKELVIN=3.15774646D+05)
      PARAMETER (BOLTZK=1.0D+00/TOKELVIN)
      PARAMETER (ZERO=0.0D+00)
      PARAMETER (PT5=0.50D+00)
      PARAMETER (TWO=2.0D+00)
      PARAMETER (ONE=1.0D+00)
      PARAMETER (FOURTHIRD=4.0D+00/3.0D+00)
      PARAMETER (ONETHIRD=1.0D+00/3.0D+00)
      PARAMETER (TOANGS=0.52917724924D+00)
      PARAMETER (TOBOHR=1.0D+00/TOANGS)
C
      DIMENSION CORD(3,*),VEL(3,*),QMVEL(3,*),ZMAS(*),QMZMAS(*),
     *          ONEMAS(*),QM1MAS(*),LISTQM(*),TIMAT(3,3)
      DIMENSION INEAR(16384),LSTPOL(*)
C
      COMMON /FFDFS / TIMDFS,QDION,AMION,TEFF,NDFS,NATMGAS,
     *                LFFDFSC,
     *                LFFDFSC0,LFFDFSA,LFFDFSN,LFFDFCOM,KDFS,LFFDFSCAV
      COMMON /FFENGY/ EN12,EN123,EN123R4,EN123B4,EN234W1,ENCHAR,ENLJR,
     *                ENLJD,ENPOL,XENPOL,ENRXN,XENRXN,ENRXNPOL,ENRXNR,
     *                EN12312,ENQUANP(30),
     *                ENBIAS,ENCENT,ENUCCH,ENCMAP,ENPOT,ENKIN,ENTOT
      COMMON /FFPARM/ NFFAT,NBOND,NANGL,NDIHR,NDIHB,NCMAP,NWAGG,
     *                N1213J,N14J,NLKQMM,IDOCHG,IDOPOL,IDOLJ,IDOCMAP
      COMMON /FFMDPA/ DT,DT2,TEMP0,PRES0,POLTOL,VIR(3),PMEAN,VOLAV,
     *                ENPAV,ENKAV,TEMPAV,BERENDT,BERENDP,VELMAX,
     *                PMEANX,PMEANY,PMEANZ,NSTEP,KMASTER,KOUTACT(2),
     *                IHESS,INTALG,ITSTAT,IPSTAT,JOUT,KOUT,LOUT
      COMMON /FFPBSW/ XBOX,YBOX,ZBOX,SWF,SWFDX,SWFDY,SWFDZ,
     *                SWRA,ONESWRA,SWRA2,ONESWRA2,
     *                SWRB,ONESWRB,SWRB2,ONESWRB2,ONESWRB4,
     *                SWFDUM3,SWFDUM4,SWFDUM5,
     *                SWRAQ,ONESWRAQ,SWRAQ2,ONESWRAQ2,
     *                SWRBQ,ONESWRBQ,SWRBQ2,ONESWRBQ2,ONESWRBQ4,
     *                SWFDUM3Q,SWFDUM4Q,SWFDUM5Q,
     *                QMSIZE,QMCX,QMCY,QMCZ,QMCXSV,QMCYSV,QMCZSV,
     *                CENTX,CENTY,CENTZ,BUFWID1,BUFWID2,RDAMP,
     *                EFIELDX,EFIELDY,EFIELDZ,QMCXSV2,QMCYSV2,QMCZSV2,
     *                EPS1RB,EPS1RB3,ONEXBOX,ONEYBOX,ONEZBOX,
     *                LQMCT,MXLIST1,NTODO,NTODOSV,NTODOQ,
     *                ISWITCH,ISHIFT,IPOLSHF,
     *                LFFLSTCELL,LFFCORDSV,
     *                LFFPOLSV,LFFCORDSV2,LFFNONLS2,LFFCORDSVQ,
     *                LFFMVFASTS2,LFFMVFASTS3,LFFMVFASTS4,
     *                LFFMVFASTL2,LFFMVFASTL3,LFFMVFASTL4,
     *                MXCHECK,MXLIST2,NTODO2,NTODO2SV
      COMMON /FFRATT/ RATOLC,RATOLV,SCALRAT,VIRRAT(3),IRATTLE,JRATTLE,
     *                NRATTLE,MXRATT,LFFOLDCORD,LFFLSTRAT,LFFDSTRAT,
     *                LFFVELSV,IRATQM
      COMMON /FFRMSD/ DIPT(8),TIMGYRA,TIMRALL,
     *                NATPDB,NGYRA,NDIEL,NRALL,
     *                LFFRALL0,NRMSD,LFFRMSD0,KFREEAB(201),
     *                NRIJMM,IJRMM(2,100),NRIJQM,IJRQM(2,100),
     *                NAIJKMM,IJKMM(3,100),NAIJKQM,IJKQM(3,100),
     *                NFIXMM,IFIXMM(200),NFIXQM,IFIXQM(200)
      COMMON /FFSPH / SPHRAD,SPHSIG,SPHEPS,IADDWAT
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /IOFILE/ IR,IW,IP,IJK,IJKT,IDAF,NAV,IODA(950)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
C
C     NANDUN THELLAMUREGE, HUI LI, MAR 2011, LINCOLN
C
      NNNATM = NFFAT + NAT - LISTQM(NFFAT+NAT+1)
      IF(NNNATM.GT.2)     NDOF = 3*NNNATM - 6
      IF(NNNATM.EQ.2)     NDOF = 3*NNNATM - 5
      IF(NRATTLE.GT.0)    NDOF = NDOF - NRATTLE
      IF(XBOX.LT.1.0D+30) NDOF = NDOF + 1
      IF(YBOX.LT.1.0D+30) NDOF = NDOF + 1
      IF(ZBOX.LT.1.0D+30) NDOF = NDOF + 1
      IF(NNNATM.EQ.1)     NDOF = 3
C
C     -- VELOCITY SCALING --
C
      IF(ITSTAT.GT.0) THEN
      IF(ABS(EFIELDX)+ABS(EFIELDY)+ABS(EFIELDZ).EQ.ZERO) THEN
C
C        -- (1) BERENDSEN THERMOSTAT --
         IF (ITSTAT.EQ.1) THEN
            DUM1  = ZERO
            IF(TEMP.GT.ZERO) DUM1  = TEMP0/TEMP
            DUM2  = DUM1 - 1.0D+00
            ONETT = 2.418884326505D-17/BERENDT
            IF((TEMP-TEMP0).GT.100.0D+00)
     *      ONETT = 0.1D+00/DT
            IF((TEMP-TEMP0).GT.200.0D+00)
     *      ONETT = 1.0D+00/DT
            TFACT = SQRT(1.0D+00 + (DT*ONETT)*DUM2)
            DO IAT = 1, NAT
               QMVEL(1,IAT) = QMVEL(1,IAT)*TFACT
               QMVEL(2,IAT) = QMVEL(2,IAT)*TFACT
               QMVEL(3,IAT) = QMVEL(3,IAT)*TFACT
               IF(ABS(QMVEL(1,IAT)).GT.VELMAX)
     *                QMVEL(1,IAT) = SIGN(ONE,QMVEL(1,IAT))*VELMAX
               IF(ABS(QMVEL(2,IAT)).GT.VELMAX)
     *                QMVEL(2,IAT) = SIGN(ONE,QMVEL(2,IAT))*VELMAX
               IF(ABS(QMVEL(3,IAT)).GT.VELMAX)
     *                QMVEL(3,IAT) = SIGN(ONE,QMVEL(3,IAT))*VELMAX
               IF(LISTQM(NFFAT+IAT).GT.0)THEN
                  IFFAT=LISTQM(NFFAT+IAT)
                  VEL(1,IFFAT)=QMVEL(1,IAT)
                  VEL(2,IFFAT)=QMVEL(2,IAT)
                  VEL(3,IFFAT)=QMVEL(3,IAT)
               END IF
            ENDDO
            DO IFFAT = 1, NFFAT
               IF(LISTQM(IFFAT).EQ.0.OR.LISTQM(IFFAT).GT.NAT) THEN
                  VEL(1,IFFAT) = VEL(1,IFFAT)*TFACT
                  VEL(2,IFFAT) = VEL(2,IFFAT)*TFACT
                  VEL(3,IFFAT) = VEL(3,IFFAT)*TFACT
                  IF(ABS(VEL(1,IFFAT)).GT.VELMAX)
     *                   VEL(1,IFFAT) = SIGN(ONE,VEL(1,IFFAT))*VELMAX
                  IF(ABS(VEL(2,IFFAT)).GT.VELMAX)
     *                   VEL(2,IFFAT) = SIGN(ONE,VEL(2,IFFAT))*VELMAX
                  IF(ABS(VEL(3,IFFAT)).GT.VELMAX)
     *                   VEL(3,IFFAT) = SIGN(ONE,VEL(3,IFFAT))*VELMAX
               END IF
            ENDDO
         END IF
C
C        -- (2) ANDERSEN THERMOSTAT --
C
         IF (ITSTAT.EQ.2) THEN
            FREQ   = 0.2D+00
            FACT   = SQRT(DBLE(NDOF)/DBLE(3*NNNATM))
            DO IAT = 1,NAT
               CALL FFRAND(TEST)
               IF(TEST.LT.FREQ) THEN
                  SIGMA=FACT*SQRT(BOLTZK*TEMP0*QM1MAS(IAT))
                  DO IDIM = 1,3
                     CALL FFRAND(U1)
                     CALL FFRAND(U2)
                     SET=SQRT(-2.0D+00*LOG(U1))*COS(2.0D+00*PI*U2)
                     QMVEL(IDIM,IAT)=SIGMA*SET
                     IF(ABS(QMVEL(IDIM,IAT)).GT.VELMAX)
     *               QMVEL(IDIM,IAT)=SIGN(ONE,QMVEL(IDIM,IAT))*VELMAX
                  ENDDO
               END IF
               IF(LISTQM(NFFAT+IAT).GT.0)THEN
                  IFFAT=LISTQM(NFFAT+IAT)
                  VEL(1,IFFAT)=QMVEL(1,IAT)
                  VEL(2,IFFAT)=QMVEL(2,IAT)
                  VEL(3,IFFAT)=QMVEL(3,IAT)
               END IF
            ENDDO
            DO IFFAT = 1,NFFAT
               IF(LISTQM(IFFAT).EQ.0.OR.LISTQM(IFFAT).GT.NAT) THEN
                  CALL FFRAND(TEST)
                  IF(TEST.LT.FREQ) THEN
                     SIGMA=FACT*SQRT(BOLTZK*TEMP0*ONEMAS(IFFAT))
                     DO IDIM = 1,3
                        CALL FFRAND(U1)
                        CALL FFRAND(U2)
                        SET=SQRT(-2.0D+00*LOG(U1))*COS(2.0D+00*PI*U2)
                        VEL(IDIM,IFFAT)=SIGMA*SET
                        IF(ABS(VEL(IDIM,IFFAT)).GT.VELMAX)
     *                  VEL(IDIM,IFFAT)=SIGN(ONE,VEL(IDIM,IFFAT))*
     *                  VELMAX
                     ENDDO
                  END IF
               END IF
            ENDDO
         END IF
C
C        -- REMOVE TRANSLATIONAL MOMENTUM --
C           BEFORE CALCULATING T
C
         XCOM = ZERO
         YCOM = ZERO
         ZCOM = ZERO
         XMNT = ZERO
         YMNT = ZERO
         ZMNT = ZERO
         TMAS = ZERO
         DO IAT = 1,NAT
            XCOM = XCOM + C(1,IAT)*QMZMAS(IAT)
            YCOM = YCOM + C(2,IAT)*QMZMAS(IAT)
            ZCOM = ZCOM + C(3,IAT)*QMZMAS(IAT)
            XMNT = XMNT + QMVEL(1,IAT)*QMZMAS(IAT)
            YMNT = YMNT + QMVEL(2,IAT)*QMZMAS(IAT)
            ZMNT = ZMNT + QMVEL(3,IAT)*QMZMAS(IAT)
            TMAS = TMAS + QMZMAS(IAT)
         ENDDO
         DO IFFAT = 1,NFFAT
         IF(LISTQM(IFFAT).EQ.0.OR.LISTQM(IFFAT).GT.NAT) THEN
            XCOM = XCOM + CORD(1,IFFAT)*ZMAS(IFFAT)
            YCOM = YCOM + CORD(2,IFFAT)*ZMAS(IFFAT)
            ZCOM = ZCOM + CORD(3,IFFAT)*ZMAS(IFFAT)
            XMNT = XMNT + VEL(1,IFFAT)*ZMAS(IFFAT)
            YMNT = YMNT + VEL(2,IFFAT)*ZMAS(IFFAT)
            ZMNT = ZMNT + VEL(3,IFFAT)*ZMAS(IFFAT)
            TMAS = TMAS + ZMAS(IFFAT)
         END IF
         ENDDO
         ONETMAS = ONE/TMAS
         XCOM = XCOM*ONETMAS
         YCOM = YCOM*ONETMAS
         ZCOM = ZCOM*ONETMAS
         XMNT = XMNT*ONETMAS
         YMNT = YMNT*ONETMAS
         ZMNT = ZMNT*ONETMAS
C        -- MUST CORRECT THE TRANSLATION BEFORE ROTATION
         DO IAT = 1,NAT
            QMVEL(1,IAT) = QMVEL(1,IAT) - XMNT
            QMVEL(2,IAT) = QMVEL(2,IAT) - YMNT
            QMVEL(3,IAT) = QMVEL(3,IAT) - ZMNT
         ENDDO
         DO IFFAT = 1,NFFAT
            VEL(1,IFFAT) = VEL(1,IFFAT) - XMNT
            VEL(2,IFFAT) = VEL(2,IFFAT) - YMNT
            VEL(3,IFFAT) = VEL(3,IFFAT) - ZMNT
         ENDDO
C
C        -- ANGULAR MOMENTUM AROUND COM
C           MOMENT OF INTERTIA AROUND COM
C
         SUMAX= ZERO
         SUMAY= ZERO
         SUMAZ= ZERO
         TXX  = ZERO
         TYY  = ZERO
         TZZ  = ZERO
         TXY  = ZERO
         TXZ  = ZERO
         TYZ  = ZERO
         DO IAT=1,NAT
         SUMAX=SUMAX+((C(2,IAT)-YCOM)*QMVEL(3,IAT)
     *               -(C(3,IAT)-ZCOM)*QMVEL(2,IAT))*QMZMAS(IAT)
         SUMAY=SUMAY+((C(3,IAT)-ZCOM)*QMVEL(1,IAT)
     *               -(C(1,IAT)-XCOM)*QMVEL(3,IAT))*QMZMAS(IAT)
         SUMAZ=SUMAZ+((C(1,IAT)-XCOM)*QMVEL(2,IAT)
     *               -(C(2,IAT)-YCOM)*QMVEL(1,IAT))*QMZMAS(IAT)
         TXX = TXX + QMZMAS(IAT)*((C(2,IAT)-YCOM)*(C(2,IAT)-YCOM)
     *                           +(C(3,IAT)-ZCOM)*(C(3,IAT)-ZCOM))
         TYY = TYY + QMZMAS(IAT)*((C(1,IAT)-XCOM)*(C(1,IAT)-XCOM)
     *                           +(C(3,IAT)-ZCOM)*(C(3,IAT)-ZCOM))
         TZZ = TZZ + QMZMAS(IAT)*((C(1,IAT)-XCOM)*(C(1,IAT)-XCOM)
     *                           +(C(2,IAT)-YCOM)*(C(2,IAT)-YCOM))
         TXY = TXY - QMZMAS(IAT)* (C(1,IAT)-XCOM)*(C(2,IAT)-YCOM)
         TXZ = TXZ - QMZMAS(IAT)* (C(1,IAT)-XCOM)*(C(3,IAT)-ZCOM)
         TYZ = TYZ - QMZMAS(IAT)* (C(2,IAT)-YCOM)*(C(3,IAT)-ZCOM)
         ENDDO
         DO IFFAT=1,NFFAT
         IF(LISTQM(IFFAT).EQ.0.OR.LISTQM(IFFAT).GT.NAT) THEN
         SUMAX=SUMAX+((CORD(2,IFFAT)-YCOM)*VEL(3,IFFAT)
     *               -(CORD(3,IFFAT)-ZCOM)*VEL(2,IFFAT))*ZMAS(IFFAT)
         SUMAY=SUMAY+((CORD(3,IFFAT)-ZCOM)*VEL(1,IFFAT)
     *               -(CORD(1,IFFAT)-XCOM)*VEL(3,IFFAT))*ZMAS(IFFAT)
         SUMAZ=SUMAZ+((CORD(1,IFFAT)-XCOM)*VEL(2,IFFAT)
     *               -(CORD(2,IFFAT)-YCOM)*VEL(1,IFFAT))*ZMAS(IFFAT)
         TXX = TXX + ZMAS(IFFAT)*
     *              ((CORD(2,IFFAT)-YCOM)*(CORD(2,IFFAT)-YCOM)
     *              +(CORD(3,IFFAT)-ZCOM)*(CORD(3,IFFAT)-ZCOM))
         TYY = TYY + ZMAS(IFFAT)*
     *              ((CORD(1,IFFAT)-XCOM)*(CORD(1,IFFAT)-XCOM)
     *              +(CORD(3,IFFAT)-ZCOM)*(CORD(3,IFFAT)-ZCOM))
         TZZ = TZZ + ZMAS(IFFAT)*
     *              ((CORD(1,IFFAT)-XCOM)*(CORD(1,IFFAT)-XCOM)
     *              +(CORD(2,IFFAT)-YCOM)*(CORD(2,IFFAT)-YCOM))
         TXY = TXY - ZMAS(IFFAT)*
     *               (CORD(1,IFFAT)-XCOM)*(CORD(2,IFFAT)-YCOM)
         TXZ = TXZ - ZMAS(IFFAT)*
     *               (CORD(1,IFFAT)-XCOM)*(CORD(3,IFFAT)-ZCOM)
         TYZ = TYZ - ZMAS(IFFAT)*
     *               (CORD(2,IFFAT)-YCOM)*(CORD(3,IFFAT)-ZCOM)
         END IF
         ENDDO
C
C        -- CALCULATE INVERSE INERTIA TENSOR
         IF(NFFAT.GE.3) THEN
         CALL TINV(TXX,TYY,TZZ,TXY,TXZ,TYZ,TIMAT)
C        -- MULTIPLY INVERSE INERTIA TENSOR WITH ANGULAR MOMENTUM
         PRODX = SUMAX*TIMAT(1,1) + SUMAY*TIMAT(1,2) + SUMAZ*TIMAT(1,3)
         PRODY = SUMAX*TIMAT(2,1) + SUMAY*TIMAT(2,2) + SUMAZ*TIMAT(2,3)
         PRODZ = SUMAX*TIMAT(3,1) + SUMAY*TIMAT(3,2) + SUMAZ*TIMAT(3,3)
         ELSE
         PRODX = ZERO
         PRODY = ZERO
         PRODZ = ZERO
         END IF
         DO IAT =1 ,NAT
C        -- GET THE CROSS PRODUCT WITH THE DISTANCES --
         CROSX = PRODY*(C(3,IAT)-ZCOM) - PRODZ*(C(2,IAT)-YCOM)
         CROSY = PRODZ*(C(1,IAT)-XCOM) - PRODX*(C(3,IAT)-ZCOM)
         CROSZ = PRODX*(C(2,IAT)-YCOM) - PRODY*(C(1,IAT)-XCOM)
         QMVEL(1,IAT)=QMVEL(1,IAT)-CROSX
         QMVEL(2,IAT)=QMVEL(2,IAT)-CROSY
         QMVEL(3,IAT)=QMVEL(3,IAT)-CROSZ
         ENDDO
         DO IFFAT =1,NFFAT
C        -- GET THE CROSS PRODUCT WITH THE DISTANCES --
         CROSX = PRODY*(CORD(3,IFFAT)-ZCOM)-PRODZ*(CORD(2,IFFAT)-YCOM)
         CROSY = PRODZ*(CORD(1,IFFAT)-XCOM)-PRODX*(CORD(3,IFFAT)-ZCOM)
         CROSZ = PRODX*(CORD(2,IFFAT)-YCOM)-PRODY*(CORD(1,IFFAT)-XCOM)
         VEL(1,IFFAT)=VEL(1,IFFAT)-CROSX
         VEL(2,IFFAT)=VEL(2,IFFAT)-CROSY
         VEL(3,IFFAT)=VEL(3,IFFAT)-CROSZ
         ENDDO
C
C
         ENKIN = ZERO
         DO IAT =1, NAT
            DUMY = ZERO
            DO III = 1, 3
               DUMY = DUMY + QMVEL(III,IAT)*QMVEL(III,IAT)
            ENDDO
            ENKIN = ENKIN + DUMY*QMZMAS(IAT)
         ENDDO
         DO IFFAT = 1, NFFAT
            IF(LISTQM(IFFAT).EQ.0.OR.LISTQM(IFFAT).GT.NAT) THEN
               DUMY = ZERO
               DO III = 1,3
                  DUMY = DUMY + VEL(III,IFFAT)*VEL(III,IFFAT)
               ENDDO
               ENKIN = ENKIN + DUMY*ZMAS(IFFAT)
            END IF
         ENDDO
         ENKIN = PT5*ENKIN
         TEMP  = TOKELVIN*TWO*ENKIN/NDOF
C
         IF(MASWRK.AND.(MOD(ISTEP,JOUT).EQ.0.OR.ISTEP.EQ.NSTEP)) THEN
            WRITE(IW,9001) 'RESCALED TEMPERATURE    =', TEMP
            WRITE(IW,*)' '
         END IF
C
      ELSE
C
C        -- (1) BERENDSEN THERMOSTAT --
         IF (ITSTAT.EQ.1) THEN
C
C           - FOR DRIFT ION IN MD-CCS -
C
C           - DETERMINE TRANSLATIONAL MOMENTUM
            XCOM = ZERO
            YCOM = ZERO
            ZCOM = ZERO
            XMNT = ZERO
            YMNT = ZERO
            ZMNT = ZERO
            TMAS = ZERO
            DO IFFAT = 1, NATPDB
               XCOM = XCOM + CORD(1,IFFAT)*ZMAS(IFFAT)
               YCOM = YCOM + CORD(2,IFFAT)*ZMAS(IFFAT)
               ZCOM = ZCOM + CORD(3,IFFAT)*ZMAS(IFFAT)
               XMNT = XMNT + VEL(1,IFFAT)*ZMAS(IFFAT)
               YMNT = YMNT + VEL(2,IFFAT)*ZMAS(IFFAT)
               ZMNT = ZMNT + VEL(3,IFFAT)*ZMAS(IFFAT)
               TMAS = TMAS + ZMAS(IFFAT)
            ENDDO
            ONETMAS = ONE/TMAS
            XCOM = XCOM*ONETMAS
            YCOM = YCOM*ONETMAS
            ZCOM = ZCOM*ONETMAS
            XMNT = XMNT*ONETMAS
            YMNT = YMNT*ONETMAS
            ZMNT = ZMNT*ONETMAS
C
C           - REMOVE TRANSLATION
            DO IFFAT = 1, NATPDB
               VEL(1,IFFAT) = VEL(1,IFFAT) - XMNT
               VEL(2,IFFAT) = VEL(2,IFFAT) - YMNT
               VEL(3,IFFAT) = VEL(3,IFFAT) - ZMNT
            ENDDO
C
C           - NOW KINETIC ENERGY (NOT COUNTING TRANSLATION)
            ENKINM = ZERO
            DO IFFAT = 1, NATPDB
                  DUMY = ZERO
                  DO III = 1,3
                     DUMY = DUMY + VEL(III,IFFAT)*VEL(III,IFFAT)
                  ENDDO
                  ENKINM = ENKINM + DUMY*ZMAS(IFFAT)
            ENDDO
            ENKINM = PT5*ENKINM
            TEMPM  = TOKELVIN*TWO*ENKINM/(NATPDB*3-3)  !  CAN ROTATE
            IF(NATPDB.EQ.1) TEMPM = TEFF
C
C           - TEMPERATURE SCALING
C
C           DUM1  = ZERO
C           IF(TEMPM.GT.ZERO) DUM1  = TEFF/TEMPM       !  USE TEFF
C           DUM2  = DUM1 - 1.0D+00
C           ONETT = 0.0001D+00/DT                      !  USE 0.0001
C            ! IF NEED 0.0001D+00/DT, USE BERENDT=10000*DT
C            ! IF NEED 0.0010D+00/DT, USE BERENDT=1000*DT
C           ONETT = 2.418884326505D-17/BERENDT
C           TFACT = SQRT(1.0D+00 + (DT*ONETT)*DUM2)
C           DO IFFAT = 1, NATPDB
C                 VEL(1,IFFAT) = VEL(1,IFFAT)*TFACT
C                 VEL(2,IFFAT) = VEL(2,IFFAT)*TFACT
C                 VEL(3,IFFAT) = VEL(3,IFFAT)*TFACT
C                 IF(ABS(VEL(1,IFFAT)).GT.VELMAX)
C    *                   VEL(1,IFFAT) = SIGN(ONE,VEL(1,IFFAT))*VELMAX
C                 IF(ABS(VEL(2,IFFAT)).GT.VELMAX)
C    *                   VEL(2,IFFAT) = SIGN(ONE,VEL(2,IFFAT))*VELMAX
C                 IF(ABS(VEL(3,IFFAT)).GT.VELMAX)
C    *                   VEL(3,IFFAT) = SIGN(ONE,VEL(3,IFFAT))*VELMAX
C           ENDDO
C
C           - ADD TRANSLATION BACK
C             BUT PURPOSELY ADD ZERO FOR ISTEP=0
            IF(ISTEP.GT.0) THEN
            DO IFFAT = 1, NATPDB
               VEL(1,IFFAT) = VEL(1,IFFAT) + XMNT
               VEL(2,IFFAT) = VEL(2,IFFAT) + YMNT
               VEL(3,IFFAT) = VEL(3,IFFAT) + ZMNT
            ENDDO
            END IF
C
C
C           - FOR BUFFER GAS MOLECULES -
C
C           - DETERMINE IF GAS MOLECULE IS NEAR ION
            CALL VICLR(INEAR,1,NFFAT)
            RCUT2 = (12.0D+00*TOBOHR)**2
            DO III = 1,LSTPOL(NFFAT)
               IFFAT = LSTPOL(III)
               IF(MOD(IFFAT+NATMGAS-1-NATPDB,NATMGAS).EQ.0) THEN
               INEAR(IFFAT) = 0
               IYES=0
               DO JFFAT=1,NATPDB
                  CX    = CORD(1,IFFAT) - CORD(1,JFFAT)
                  CY    = CORD(2,IFFAT) - CORD(2,JFFAT)
                  CZ    = CORD(3,IFFAT) - CORD(3,JFFAT)
                  PBCX  = XBOX*ANINT(CX*ONEXBOX)
                  PBCY  = YBOX*ANINT(CY*ONEYBOX)
                  PBCZ  = ZBOX*ANINT(CZ*ONEZBOX)
                  CX    = CX - PBCX
                  CY    = CY - PBCY
                  CZ    = CZ - PBCZ
                  R2    = CX*CX + CY*CY + CZ*CZ
                  IF(R2.LT.RCUT2) IYES=1
                  IF(NATMGAS.EQ.2)THEN
                  CX    = CORD(1,IFFAT+1) - CORD(1,JFFAT)
                  CY    = CORD(2,IFFAT+1) - CORD(2,JFFAT)
                  CZ    = CORD(3,IFFAT+1) - CORD(3,JFFAT)
                  PBCX  = XBOX*ANINT(CX*ONEXBOX)
                  PBCY  = YBOX*ANINT(CY*ONEYBOX)
                  PBCZ  = ZBOX*ANINT(CZ*ONEZBOX)
                  CX    = CX - PBCX
                  CY    = CY - PBCY
                  CZ    = CZ - PBCZ
                  R2    = CX*CX + CY*CY + CZ*CZ
                  IF(R2.LT.RCUT2) IYES=1
                  END IF
                  IF(NATMGAS.EQ.3)THEN
                  CX    = CORD(1,IFFAT+2) - CORD(1,JFFAT)
                  CY    = CORD(2,IFFAT+2) - CORD(2,JFFAT)
                  CZ    = CORD(3,IFFAT+2) - CORD(3,JFFAT)
                  PBCX  = XBOX*ANINT(CX*ONEXBOX)
                  PBCY  = YBOX*ANINT(CY*ONEYBOX)
                  PBCZ  = ZBOX*ANINT(CZ*ONEZBOX)
                  CX    = CX - PBCX
                  CY    = CY - PBCY
                  CZ    = CZ - PBCZ
                  R2    = CX*CX + CY*CY + CZ*CZ
                  IF(R2.LT.RCUT2) IYES=1
                  END IF
               ENDDO
               INEAR(IFFAT)=IYES
               END IF
            ENDDO
C
C           - DETERMINE TRANSLATIONAL MOMENTUM
            XCOM = ZERO
            YCOM = ZERO
            ZCOM = ZERO
            XMNT = ZERO
            YMNT = ZERO
            ZMNT = ZERO
            TMAS = ZERO
            DO IFFAT = NATPDB+1,NFFAT-NATMGAS+1,NATMGAS
               IYES=INEAR(IFFAT)
               IF(IYES.EQ.0) THEN
               XCOM = XCOM + CORD(1,IFFAT)*ZMAS(IFFAT)
               YCOM = YCOM + CORD(2,IFFAT)*ZMAS(IFFAT)
               ZCOM = ZCOM + CORD(3,IFFAT)*ZMAS(IFFAT)
               XMNT = XMNT + VEL(1,IFFAT)*ZMAS(IFFAT)
               YMNT = YMNT + VEL(2,IFFAT)*ZMAS(IFFAT)
               ZMNT = ZMNT + VEL(3,IFFAT)*ZMAS(IFFAT)
               TMAS = TMAS + ZMAS(IFFAT)
               IF(NATMGAS.EQ.2)THEN
               XCOM = XCOM + CORD(1,IFFAT+1)*ZMAS(IFFAT+1)
               YCOM = YCOM + CORD(2,IFFAT+1)*ZMAS(IFFAT+1)
               ZCOM = ZCOM + CORD(3,IFFAT+1)*ZMAS(IFFAT+1)
               XMNT = XMNT + VEL(1,IFFAT+1)*ZMAS(IFFAT+1)
               YMNT = YMNT + VEL(2,IFFAT+1)*ZMAS(IFFAT+1)
               ZMNT = ZMNT + VEL(3,IFFAT+1)*ZMAS(IFFAT+1)
               TMAS = TMAS + ZMAS(IFFAT+1)
               END IF
               IF(NATMGAS.EQ.3)THEN
               XCOM = XCOM + CORD(1,IFFAT+2)*ZMAS(IFFAT+2)
               YCOM = YCOM + CORD(2,IFFAT+2)*ZMAS(IFFAT+2)
               ZCOM = ZCOM + CORD(3,IFFAT+2)*ZMAS(IFFAT+2)
               XMNT = XMNT + VEL(1,IFFAT+2)*ZMAS(IFFAT+2)
               YMNT = YMNT + VEL(2,IFFAT+2)*ZMAS(IFFAT+2)
               ZMNT = ZMNT + VEL(3,IFFAT+2)*ZMAS(IFFAT+2)
               TMAS = TMAS + ZMAS(IFFAT+2)
               END IF
               END IF
            ENDDO
            ONETMAS = ONE/TMAS
            XCOM = XCOM*ONETMAS
            YCOM = YCOM*ONETMAS
            ZCOM = ZCOM*ONETMAS
            XMNT = XMNT*ONETMAS
            YMNT = YMNT*ONETMAS
            ZMNT = ZMNT*ONETMAS
C           -- CORRECT TRANSLATION
            DO IFFAT = NATPDB+1,NFFAT-NATMGAS+1,NATMGAS
               IYES=INEAR(IFFAT)
               IF(IYES.EQ.0) THEN
               VEL(1,IFFAT) = VEL(1,IFFAT) - XMNT
               VEL(2,IFFAT) = VEL(2,IFFAT) - YMNT
               VEL(3,IFFAT) = VEL(3,IFFAT) - ZMNT
               IF(NATMGAS.EQ.2)THEN
               VEL(1,IFFAT+1) = VEL(1,IFFAT+1) - XMNT
               VEL(2,IFFAT+1) = VEL(2,IFFAT+1) - YMNT
               VEL(3,IFFAT+1) = VEL(3,IFFAT+1) - ZMNT
               END IF
               IF(NATMGAS.EQ.3)THEN
               VEL(1,IFFAT+2) = VEL(1,IFFAT+2) - XMNT
               VEL(2,IFFAT+2) = VEL(2,IFFAT+2) - YMNT
               VEL(3,IFFAT+2) = VEL(3,IFFAT+2) - ZMNT
               END IF
               END IF
            ENDDO
C
C           -- ANGULAR MOMENTUM AND MOMENT OF INERTIA
C
            SUMAX= ZERO
            SUMAY= ZERO
            SUMAZ= ZERO
            TXX  = ZERO
            TYY  = ZERO
            TZZ  = ZERO
            TXY  = ZERO
            TXZ  = ZERO
            TYZ  = ZERO
            DO IFFAT=NATPDB+1,NFFAT-NATMGAS+1,NATMGAS
               IYES=INEAR(IFFAT)
               IF(IYES.EQ.0) THEN
               SUMAX=SUMAX+((CORD(2,IFFAT)-YCOM)*VEL(3,IFFAT)
     *                  -(CORD(3,IFFAT)-ZCOM)*VEL(2,IFFAT))*ZMAS(IFFAT)
               SUMAY=SUMAY+((CORD(3,IFFAT)-ZCOM)*VEL(1,IFFAT)
     *                  -(CORD(1,IFFAT)-XCOM)*VEL(3,IFFAT))*ZMAS(IFFAT)
               SUMAZ=SUMAZ+((CORD(1,IFFAT)-XCOM)*VEL(2,IFFAT)
     *                  -(CORD(2,IFFAT)-YCOM)*VEL(1,IFFAT))*ZMAS(IFFAT)
               TXX = TXX + ZMAS(IFFAT)*
     *                 ((CORD(2,IFFAT)-YCOM)*(CORD(2,IFFAT)-YCOM)
     *                 +(CORD(3,IFFAT)-ZCOM)*(CORD(3,IFFAT)-ZCOM))
               TYY = TYY + ZMAS(IFFAT)*
     *                 ((CORD(1,IFFAT)-XCOM)*(CORD(1,IFFAT)-XCOM)
     *                 +(CORD(3,IFFAT)-ZCOM)*(CORD(3,IFFAT)-ZCOM))
               TZZ = TZZ + ZMAS(IFFAT)*
     *                 ((CORD(1,IFFAT)-XCOM)*(CORD(1,IFFAT)-XCOM)
     *                 +(CORD(2,IFFAT)-YCOM)*(CORD(2,IFFAT)-YCOM))
               TXY = TXY - ZMAS(IFFAT)*
     *                  (CORD(1,IFFAT)-XCOM)*(CORD(2,IFFAT)-YCOM)
               TXZ = TXZ - ZMAS(IFFAT)*
     *                  (CORD(1,IFFAT)-XCOM)*(CORD(3,IFFAT)-ZCOM)
               TYZ = TYZ - ZMAS(IFFAT)*
     *                  (CORD(2,IFFAT)-YCOM)*(CORD(3,IFFAT)-ZCOM)
               IF(NATMGAS.EQ.2)THEN
               SUMAX=SUMAX+((CORD(2,IFFAT+1)-YCOM)*VEL(3,IFFAT+1)
     *         -(CORD(3,IFFAT+1)-ZCOM)*VEL(2,IFFAT+1))*ZMAS(IFFAT+1)
               SUMAY=SUMAY+((CORD(3,IFFAT+1)-ZCOM)*VEL(1,IFFAT+1)
     *         -(CORD(1,IFFAT+1)-XCOM)*VEL(3,IFFAT+1))*ZMAS(IFFAT+1)
               SUMAZ=SUMAZ+((CORD(1,IFFAT+1)-XCOM)*VEL(2,IFFAT+1)
     *         -(CORD(2,IFFAT+1)-YCOM)*VEL(1,IFFAT+1))*ZMAS(IFFAT+1)
               TXX = TXX + ZMAS(IFFAT+1)*
     *            ((CORD(2,IFFAT+1)-YCOM)*(CORD(2,IFFAT+1)-YCOM)
     *            +(CORD(3,IFFAT+1)-ZCOM)*(CORD(3,IFFAT+1)-ZCOM))
               TYY = TYY + ZMAS(IFFAT+1)*
     *            ((CORD(1,IFFAT+1)-XCOM)*(CORD(1,IFFAT+1)-XCOM)
     *            +(CORD(3,IFFAT+1)-ZCOM)*(CORD(3,IFFAT+1)-ZCOM))
               TZZ = TZZ + ZMAS(IFFAT+1)*
     *            ((CORD(1,IFFAT+1)-XCOM)*(CORD(1,IFFAT+1)-XCOM)
     *            +(CORD(2,IFFAT+1)-YCOM)*(CORD(2,IFFAT+1)-YCOM))
               TXY = TXY - ZMAS(IFFAT+1)*
     *             (CORD(1,IFFAT+1)-XCOM)*(CORD(2,IFFAT+1)-YCOM)
               TXZ = TXZ - ZMAS(IFFAT+1)*
     *             (CORD(1,IFFAT+1)-XCOM)*(CORD(3,IFFAT+1)-ZCOM)
               TYZ = TYZ - ZMAS(IFFAT+1)*
     *             (CORD(2,IFFAT+1)-YCOM)*(CORD(3,IFFAT+1)-ZCOM)
               END IF
               IF(NATMGAS.EQ.3)THEN
               SUMAX=SUMAX+((CORD(2,IFFAT+2)-YCOM)*VEL(3,IFFAT+2)
     *         -(CORD(3,IFFAT+2)-ZCOM)*VEL(2,IFFAT+2))*ZMAS(IFFAT+2)
               SUMAY=SUMAY+((CORD(3,IFFAT+2)-ZCOM)*VEL(1,IFFAT+2)
     *         -(CORD(1,IFFAT+2)-XCOM)*VEL(3,IFFAT+2))*ZMAS(IFFAT+2)
               SUMAZ=SUMAZ+((CORD(1,IFFAT+2)-XCOM)*VEL(2,IFFAT+2)
     *         -(CORD(2,IFFAT+2)-YCOM)*VEL(1,IFFAT+2))*ZMAS(IFFAT+2)
               TXX = TXX + ZMAS(IFFAT+2)*
     *            ((CORD(2,IFFAT+2)-YCOM)*(CORD(2,IFFAT+2)-YCOM)
     *            +(CORD(3,IFFAT+2)-ZCOM)*(CORD(3,IFFAT+2)-ZCOM))
               TYY = TYY + ZMAS(IFFAT+2)*
     *            ((CORD(1,IFFAT+2)-XCOM)*(CORD(1,IFFAT+2)-XCOM)
     *            +(CORD(3,IFFAT+2)-ZCOM)*(CORD(3,IFFAT+2)-ZCOM))
               TZZ = TZZ + ZMAS(IFFAT+2)*
     *            ((CORD(1,IFFAT+2)-XCOM)*(CORD(1,IFFAT+2)-XCOM)
     *            +(CORD(2,IFFAT+2)-YCOM)*(CORD(2,IFFAT+2)-YCOM))
               TXY = TXY - ZMAS(IFFAT+2)*
     *             (CORD(1,IFFAT+2)-XCOM)*(CORD(2,IFFAT+2)-YCOM)
               TXZ = TXZ - ZMAS(IFFAT+2)*
     *             (CORD(1,IFFAT+2)-XCOM)*(CORD(3,IFFAT+2)-ZCOM)
               TYZ = TYZ - ZMAS(IFFAT+2)*
     *             (CORD(2,IFFAT+2)-YCOM)*(CORD(3,IFFAT+2)-ZCOM)
               END IF
               END IF
            ENDDO
C           -- CALCULATE INVERSE INERTIA TENSOR
            IF(NFFAT-NATPDB.GE.3) THEN
            CALL TINV(TXX,TYY,TZZ,TXY,TXZ,TYZ,TIMAT)
C           -- MULTIPLY INVERSE INERTIA TENSOR WITH ANGULAR MOMENTUM
            PRODX =SUMAX*TIMAT(1,1)+SUMAY*TIMAT(1,2)+SUMAZ*TIMAT(1,3)
            PRODY =SUMAX*TIMAT(2,1)+SUMAY*TIMAT(2,2)+SUMAZ*TIMAT(2,3)
            PRODZ =SUMAX*TIMAT(3,1)+SUMAY*TIMAT(3,2)+SUMAZ*TIMAT(3,3)
            ELSE
            PRODX = ZERO
            PRODY = ZERO
            PRODZ = ZERO
            END IF
            DO IFFAT =NATPDB+1,NFFAT-NATMGAS+1,NATMGAS
               IYES=INEAR(IFFAT)
               IF(IYES.EQ.0) THEN
C              -- GET THE CROSS PRODUCT WITH THE DISTANCES --
               CROSX = PRODY*(CORD(3,IFFAT)-ZCOM)-
     *                 PRODZ*(CORD(2,IFFAT)-YCOM)
               CROSY = PRODZ*(CORD(1,IFFAT)-XCOM)-
     *                 PRODX*(CORD(3,IFFAT)-ZCOM)
               CROSZ = PRODX*(CORD(2,IFFAT)-YCOM)-
     *                 PRODY*(CORD(1,IFFAT)-XCOM)
               VEL(1,IFFAT)=VEL(1,IFFAT)-CROSX
               VEL(2,IFFAT)=VEL(2,IFFAT)-CROSY
               VEL(3,IFFAT)=VEL(3,IFFAT)-CROSZ
               IF(NATMGAS.EQ.2)THEN
               CROSX = PRODY*(CORD(3,IFFAT+1)-ZCOM)-
     *                 PRODZ*(CORD(2,IFFAT+1)-YCOM)
               CROSY = PRODZ*(CORD(1,IFFAT+1)-XCOM)-
     *                 PRODX*(CORD(3,IFFAT+1)-ZCOM)
               CROSZ = PRODX*(CORD(2,IFFAT+1)-YCOM)-
     *                 PRODY*(CORD(1,IFFAT+1)-XCOM)
               VEL(1,IFFAT+1)=VEL(1,IFFAT+1)-CROSX
               VEL(2,IFFAT+1)=VEL(2,IFFAT+1)-CROSY
               VEL(3,IFFAT+1)=VEL(3,IFFAT+1)-CROSZ
               END IF
               IF(NATMGAS.EQ.3)THEN
               CROSX = PRODY*(CORD(3,IFFAT+2)-ZCOM)-
     *                 PRODZ*(CORD(2,IFFAT+2)-YCOM)
               CROSY = PRODZ*(CORD(1,IFFAT+2)-XCOM)-
     *                 PRODX*(CORD(3,IFFAT+2)-ZCOM)
               CROSZ = PRODX*(CORD(2,IFFAT+2)-YCOM)-
     *                 PRODY*(CORD(1,IFFAT+2)-XCOM)
               VEL(1,IFFAT+2)=VEL(1,IFFAT+2)-CROSX
               VEL(2,IFFAT+2)=VEL(2,IFFAT+2)-CROSY
               VEL(3,IFFAT+2)=VEL(3,IFFAT+2)-CROSZ
               END IF
               END IF
            ENDDO
C
            ENKING = ZERO
            NCOUNT = 0
            DO IFFAT = NATPDB+1, NFFAT-NATMGAS+1,NATMGAS
               IYES=INEAR(IFFAT)
               IF(IYES.EQ.0) THEN
                  NCOUNT = NCOUNT + 1
                  DUMY = ZERO
                  DO III = 1,3
                     DUMY = DUMY + VEL(III,IFFAT)*VEL(III,IFFAT)
                  ENDDO
                  ENKING = ENKING + DUMY*ZMAS(IFFAT)
                  IF(NATMGAS.EQ.2)THEN
                  NCOUNT = NCOUNT + 1
                  DUMY = ZERO
                  DO III = 1,3
                     DUMY = DUMY + VEL(III,IFFAT+1)*VEL(III,IFFAT+1)
                  ENDDO
                  ENKING = ENKING + DUMY*ZMAS(IFFAT+1)
                  END IF
                  IF(NATMGAS.EQ.3)THEN
                  NCOUNT = NCOUNT + 1
                  DUMY = ZERO
                  DO III = 1,3
                     DUMY = DUMY + VEL(III,IFFAT+2)*VEL(III,IFFAT+2)
                  ENDDO
                  ENKING = ENKING + DUMY*ZMAS(IFFAT+2)
                  END IF
               END IF
            ENDDO
            ENKING = PT5*ENKING
            TEMPG  = TOKELVIN*TWO*ENKING/(NCOUNT*3-6)   !  CANNOT ROTATE
C
            DUM1  = ZERO
            IF(TEMPG.GT.ZERO) DUM1  = TEMP0/TEMPG
            DUM2  = DUM1 - 1.0D+00
            ONETT = 1.0D+00/DT             !  USE 1.0
            TFACT = SQRT(1.0D+00 + (DT*ONETT)*DUM2)
            DO IFFAT = NATPDB+1, NFFAT-NATMGAS+1,NATMGAS
               IYES=INEAR(IFFAT)
               IF(IYES.EQ.0) THEN
                  VEL(1,IFFAT) = VEL(1,IFFAT)*TFACT
                  VEL(2,IFFAT) = VEL(2,IFFAT)*TFACT
                  VEL(3,IFFAT) = VEL(3,IFFAT)*TFACT
                  IF(ABS(VEL(1,IFFAT)).GT.VELMAX)
     *                   VEL(1,IFFAT) = SIGN(ONE,VEL(1,IFFAT))*VELMAX
                  IF(ABS(VEL(2,IFFAT)).GT.VELMAX)
     *                   VEL(2,IFFAT) = SIGN(ONE,VEL(2,IFFAT))*VELMAX
                  IF(ABS(VEL(3,IFFAT)).GT.VELMAX)
     *                   VEL(3,IFFAT) = SIGN(ONE,VEL(3,IFFAT))*VELMAX
                  IF(NATMGAS.EQ.2)THEN
                  VEL(1,IFFAT+1) = VEL(1,IFFAT+1)*TFACT
                  VEL(2,IFFAT+1) = VEL(2,IFFAT+1)*TFACT
                  VEL(3,IFFAT+1) = VEL(3,IFFAT+1)*TFACT
                  IF(ABS(VEL(1,IFFAT+1)).GT.VELMAX)
     *               VEL(1,IFFAT+1) = SIGN(ONE,VEL(1,IFFAT+1))*VELMAX
                  IF(ABS(VEL(2,IFFAT+1)).GT.VELMAX)
     *               VEL(2,IFFAT+1) = SIGN(ONE,VEL(2,IFFAT+1))*VELMAX
                  IF(ABS(VEL(3,IFFAT+1)).GT.VELMAX)
     *               VEL(3,IFFAT+1) = SIGN(ONE,VEL(3,IFFAT+1))*VELMAX
                  END IF
                  IF(NATMGAS.EQ.3)THEN
                  VEL(1,IFFAT+2) = VEL(1,IFFAT+2)*TFACT
                  VEL(2,IFFAT+2) = VEL(2,IFFAT+2)*TFACT
                  VEL(3,IFFAT+2) = VEL(3,IFFAT+2)*TFACT
                  IF(ABS(VEL(1,IFFAT+2)).GT.VELMAX)
     *               VEL(1,IFFAT+2) = SIGN(ONE,VEL(1,IFFAT+2))*VELMAX
                  IF(ABS(VEL(2,IFFAT+2)).GT.VELMAX)
     *               VEL(2,IFFAT+2) = SIGN(ONE,VEL(2,IFFAT+2))*VELMAX
                  IF(ABS(VEL(3,IFFAT+2)).GT.VELMAX)
     *               VEL(3,IFFAT+2) = SIGN(ONE,VEL(3,IFFAT+2))*VELMAX
                  END IF
               END IF
            ENDDO
         END IF
C
C
         ENKIN = ZERO
         DO IFFAT = 1, NFFAT
            DUMY = ZERO
            DO III = 1,3
               DUMY = DUMY + VEL(III,IFFAT)*VEL(III,IFFAT)
            ENDDO
            ENKIN = ENKIN + DUMY*ZMAS(IFFAT)
         ENDDO
         ENKIN = PT5*ENKIN
         TEMP  = TOKELVIN*TWO*ENKIN/NDOF
C
         IF(MASWRK.AND.(MOD(ISTEP,JOUT).EQ.0.OR.ISTEP.EQ.NSTEP)) THEN
            WRITE(IW,*)'-- MD SIMULATION OF ION MOBILITY PROGRAM',
     *                 ' WRITTEN BY HUI LI --'
            WRITE(IW,9001) 'RESCALED TEMPERATURE     =', TEMP
            WRITE(IW,9001) 'TEMPERATURE OF DRIFT ION =', TEMPM
            WRITE(IW,9001) 'TEMPERATURE OF BUFFER GAS=', TEMPG
            WRITE(IW,*)' '
         END IF
      END IF
      END IF
C
C     -- VOLUME SCALING --
C
      VOL   = MIN(XBOX*YBOX*ZBOX, FOURTHIRD*PI*SPHRAD**3)
      IF(IPSTAT.GT.0 .AND. ABS(PRES).GT.ZERO .AND.
     *   VOL.LE.1.0D+30.AND.MIN(XBOX,YBOX,ZBOX).LT.1.0D+30) THEN
         IF(IPSTAT.EQ.1.OR.IPSTAT.EQ.3)THEN
            BETA = 4.9D-05  ! ISOTHERMAL COMPRESSIBILITY OF WATER, BAR-1
            ONETP= 2.418884326505D-17/BERENDP
            DUM1 = MAX(-ONE,(PRES - PRES0)*BETA*(DT*ONETP))
            DUM1X= MAX(-ONE,(PRESX- PRES0)*BETA*(DT*ONETP))
            DUM1Y= MAX(-ONE,(PRESY- PRES0)*BETA*(DT*ONETP))
            DUM1Z= MAX(-ONE,(PRESZ- PRES0)*BETA*(DT*ONETP))
            SFAC = (1.0D+00+DUM1 )**ONETHIRD
            SFACX= (1.0D+00+DUM1X)**ONETHIRD
            SFACY= (1.0D+00+DUM1Y)**ONETHIRD
            SFACZ= (1.0D+00+DUM1Z)**ONETHIRD
            IF(SFAC .GT.1.0001D+00) SFAC =1.0001D+00
            IF(SFACX.GT.1.0001D+00) SFACX=1.0001D+00
            IF(SFACY.GT.1.0001D+00) SFACY=1.0001D+00
            IF(SFACZ.GT.1.0001D+00) SFACZ=1.0001D+00
            IF(SFAC .LT.0.9999D+00) SFAC =0.9999D+00
            IF(SFACX.LT.0.9999D+00) SFACX=0.9999D+00
            IF(SFACY.LT.0.9999D+00) SFACY=0.9999D+00
            IF(SFACZ.LT.0.9999D+00) SFACZ=0.9999D+00
            IF(IPSTAT.EQ.1) THEN
               SFACX = SFAC
               SFACY = SFAC
               SFACZ = SFAC
            END IF
            IF(NAT.GT.0) THEN
               CALL DSCAL(3*NAT  ,SFACX,C(1,1),3)
               CALL DSCAL(3*NAT  ,SFACY,C(2,1),3)
               CALL DSCAL(3*NAT  ,SFACZ,C(3,1),3)
            END IF
            IF(NFFAT.GT.0) THEN
               CALL DSCAL(3*NFFAT,SFACX,CORD(1,1),3)
               CALL DSCAL(3*NFFAT,SFACY,CORD(2,1),3)
               CALL DSCAL(3*NFFAT,SFACZ,CORD(3,1),3)
            END IF
            DO IAT = 1, NAT
               IF(LISTQM(NFFAT+IAT).GT.0)THEN
                  IFFAT=LISTQM(NFFAT+IAT)
                  CORD(1,IFFAT)=C(1,IAT)
                  CORD(2,IFFAT)=C(2,IAT)
                  CORD(3,IFFAT)=C(3,IAT)
               END IF
            ENDDO
            XBOX   = XBOX*SFACX
            YBOX   = YBOX*SFACY
            ZBOX   = ZBOX*SFACZ
            ONEXBOX= ONE/XBOX
            ONEYBOX= ONE/YBOX
            ONEZBOX= ONE/ZBOX
         END IF
         IF(GOPARR) THEN
            CALL DDI_BCAST(467,'F',   XBOX,1,MASTER)
            CALL DDI_BCAST(468,'F',   YBOX,1,MASTER)
            CALL DDI_BCAST(469,'F',   ZBOX,1,MASTER)
            CALL DDI_BCAST(470,'F',ONEXBOX,1,MASTER)
            CALL DDI_BCAST(471,'F',ONEYBOX,1,MASTER)
            CALL DDI_BCAST(472,'F',ONEZBOX,1,MASTER)
         END IF
      END IF
C
 9001 FORMAT(1X,A,F30.10,4X,'K')
C
      CALL FLSHBF(IW)
C
      RETURN
      END
C*MODULE QUANPOC  *DECK FFOPTX
!>
!> @brief    QuanPol main driver for geometry optimization
!>
!> @author   Nandun Thellamurege, Hui Li
!>           - May 2011
!>
!> @details  Main driver for geometry optimization and
!>           saddle point search
!>
      SUBROUTINE FFOPTX
      use mx_limits, only: mxatm
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      LOGICAL GOPARR,DSKWRK,MASWRK
      LOGICAL LINEAR
C
      COMMON /FFFIXS/ ENFIXSO,FIXEPS,FIXTOL,FIXA,FIXQ,RALLMM,RALLQM,
     *                RADMM(200),RADQM(200),NRADMM,NRADQM,IFIXSOL,
     *                LFFDAI,LFFDAIT,LFFIDDAI,LFFIDTMP,LFFTMPTS,
     *                LFFAFIX,LFFIDATOM,LFFRFIX,LFFQFIX,NTSATM,
     *                LFFQFIXMP,LFFQFIXTA,LFFQFIXXY,
     *                LFFXTSFIX,LFFYTSFIX,LFFZTSFIX,
     *                LFFVFIX1,LFFVFIX2,NCYCLE,MXFFTS,NFFTS
      COMMON /FFFRE2/ N1FFAT,N1BOND,N1ANGL,N1DIHR,N1DIHB,N1CMAP,N1WAGG,
     *                N2FFAT,N2BOND,N2ANGL,N2DIHR,N2DIHB,N2CMAP,N2WAGG,
     *                LFF2ATMNAM,LFF2CORD,LFF2ZANF,LFF2ZMAS,
     *                LFF2ONEMAS,LFF2CHARG,LFF2POL,
     *                LFF2SIG,LFF2EPS,LFF2SIG2,LFF2EPS2,
     *                LFF2BOND0,LFF2FCBOND,LFF2ANGL0,LFF2FCANGL,
     *                LFF2FCWAGG,LFF2DIHB0,LFF2FCDIHB,
     *                LFF2VROT,LFF2NNN,LFF2GAMA,LFF2IPAIR,
     *                LFF2KLIST,LFF2LLIST,LFF2MLIST,LFF2NLIST,
     *                LFF2VEL,LFF2QMVEL,LFF2CLPR,LFF2ZLPR,
     *                LFF2NLPR,LFF2MAPLST,
     *                LFFLISTB2A,NTODOA,LFFNONLSA,NTODOB,LFFNONLSB,
     *                N1213A,LFFL1213A,N1213B,LFFL1213B,
     *                N14A,LFFL14A,N14B,LFFL14B,
     *                NTODOPMA,LFFNONLSPMA,NBONDPMA,LFFLSBONDPMA,
     *                NANGLPMA,LFFLSANGLPMA,NDIHRPMA,LFFLSDIHRPMA,
     *                NDIHBPMA,LFFLSDIHBPMA,NWAGGPMA,LFFLSWAGGPMA,
     *                NCMAPPMA,LFFLSCMAPPMA,
     *                N1213PMA,LFFL1213PMA,
     *                N14PMA,LFFL14PMA,
     *                NTODOPMB,LFFNONLSPMB,NBONDPMB,LFFLSBONDPMB,
     *                NANGLPMB,LFFLSANGLPMB,NDIHRPMB,LFFLSDIHRPMB,
     *                NDIHBPMB,LFFLSDIHBPMB,NWAGGPMB,LFFLSWAGGPMB,
     *                NCMAPPMB,LFFLSCMAPPMB,
     *                N1213PMB,LFFL1213PMB,
     *                N14PMB,LFFL14PMB
      COMMON /FFPARM/ NFFAT,NBOND,NANGL,NDIHR,NDIHB,NCMAP,NWAGG,
     *                N1213J,N14J,NLKQMM,IDOCHG,IDOPOL,IDOLJ,IDOCMAP
      COMMON /FFMDPA/ DT,DT2,TEMP0,PRES0,POLTOL,VIR(3),PMEAN,VOLAV,
     *                ENPAV,ENKAV,TEMPAV,BERENDT,BERENDP,VELMAX,
     *                PMEANX,PMEANY,PMEANZ,NSTEP,KMASTER,KOUTACT(2),
     *                IHESS,INTALG,ITSTAT,IPSTAT,JOUT,KOUT,LOUT
      COMMON /FFMPNT/ LFFATMNAM,LFFCORD,LFFZANF,
     *                LFFZMAS,LFFONEMAS,LFFQMZMAS,LFFQM1MAS,
     *                LFFCHARG,LFFPOL,LFFDIP,
     *                LFFFIELD1,LFFFIELD2,LFFFIELD3,
     *                LFFSIG,LFFEPS,LFFSIG2,LFFEPS2,
     *                LFFBOND0,LFFFCBOND,
     *                LFFANGL0,LFFFCANGL,LFFFCWAGG,
     *                LFFDIHB0,LFFFCDIHB,
     *                LFFVROT,LFFNNN,LFFGAMA,LFFIPAIR,
     *                LFFKLIST,LFFLLIST,LFFL1213J,LFFL14J,
     *                LFFMLIST,LFFNLIST,LFFLKQMMM,
     *                LFFVEL,LFFQMVEL,
     *                LFFFFGRD0,LFFFFGRD1,LFFFFGRD2,
     *                LFFQMGRD0,LFFQMGRD1,LFFQMGRD2,LFFDETMP,
     *                LFFCLPR,LFFZLPR,LFFNLPR,
     *                LFFXTS,LFFYTS,LFFZTS,LFFCMAT1,
     *                LFFQRXN1,LFFQRXN2,LFFPOT1,LFFPOT2,LFFQRXNMP,
     *                LFFQRXNTA,LFFQRXNXY,LFFNONLSTQ,
     *                LFFDIPMP,LFFDIPTA,LFFDIPXY,LFFLISTQM,LFFNONLS1,
     *                LFFMAPLST,LFFCMAPCO
      COMMON /FFMPT2/ MXMMTP,LFFKBLST,LFFFCSTBD,LFFFCDIHR,
     *                LFFFCLJTP,LFFNTYPE,
     *                LFF2KBLST,LFF2FCSTBD,LFF2FCDIHR,
     *                LFF2FCLJTP,LFF2NTYPE
      COMMON /FFMPT3/ NACTMM,LACTMM(2020),LFFDIPOLD,JUMBUP,
     *                NACTQM,LACTQM(2020),LFFOLDC,LFFQMVELSV,MMHESS,
     *                LFFQMCHG,LFFQMCHGB,ISWAP,R2SWAP,DFTBMM
      COMMON /FFPBSW/ XBOX,YBOX,ZBOX,SWF,SWFDX,SWFDY,SWFDZ,
     *                SWRA,ONESWRA,SWRA2,ONESWRA2,
     *                SWRB,ONESWRB,SWRB2,ONESWRB2,ONESWRB4,
     *                SWFDUM3,SWFDUM4,SWFDUM5,
     *                SWRAQ,ONESWRAQ,SWRAQ2,ONESWRAQ2,
     *                SWRBQ,ONESWRBQ,SWRBQ2,ONESWRBQ2,ONESWRBQ4,
     *                SWFDUM3Q,SWFDUM4Q,SWFDUM5Q,
     *                QMSIZE,QMCX,QMCY,QMCZ,QMCXSV,QMCYSV,QMCZSV,
     *                CENTX,CENTY,CENTZ,BUFWID1,BUFWID2,RDAMP,
     *                EFIELDX,EFIELDY,EFIELDZ,QMCXSV2,QMCYSV2,QMCZSV2,
     *                EPS1RB,EPS1RB3,ONEXBOX,ONEYBOX,ONEZBOX,
     *                LQMCT,MXLIST1,NTODO,NTODOSV,NTODOQ,
     *                ISWITCH,ISHIFT,IPOLSHF,
     *                LFFLSTCELL,LFFCORDSV,
     *                LFFPOLSV,LFFCORDSV2,LFFNONLS2,LFFCORDSVQ,
     *                LFFMVFASTS2,LFFMVFASTS3,LFFMVFASTS4,
     *                LFFMVFASTL2,LFFMVFASTL3,LFFMVFASTL4,
     *                MXCHECK,MXLIST2,NTODO2,NTODO2SV
      COMMON /FFRXN / RXNEPS,RSPHSOL,ISPHSOL
      COMMON /FFUMBR/ UMBFC,UMBR0,UMBSIZE,
     *                NUMBBIN,NUMBATM(6),NUMBTYP,LFFUMBHIS,
     *                UM2FC,UM2R0,UM2SIZE,
     *                NUM2BIN,NUM2ATM(6),NUM2TYP,LFFUM2HIS
      COMMON /FMCOM / X(1)
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /IOFILE/ IR,IW,IP,IJK,IJKT,IDAF,NAV,IODA(950)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /RUNOPT/ RUNTYP,EXETYP,NEVALS,NGLEVL,NHLEVL
      COMMON /ZMAT  / NZMAT,NZVAR,NVAR,NSYMC,LINEAR
C
      DATA HSSIAN   /8HHESSIAN /
      DATA OPTMIZ   /8HOPTIMIZE/
      DATA SADPOINT /8HSADPOINT/
C
C     NANDUN THELLAMUREGE, HUI LI, MAY 2011, LINCOLN
C
C     -- QUANPOL OPT CANNOT USE INTERNAL COORDINATES
      IF(NZVAR.GT.0) THEN
         IF(MASWRK)WRITE(IW,'(1X,/1X,A,/)')
     *   'ERROR: QUANPOL CANNOT USE INTERNAL COORDINATES.'
         IF(MASWRK)WRITE(IW,*)' '
         CALL ABRT
      END IF
C
      NTS = ISPHSOL
C
      IF(RUNTYP.EQ.OPTMIZ.OR.RUNTYP.EQ.SADPOINT) THEN
      IF(NAT.EQ.0 .AND. NFFAT.GT.0 .AND. IHESS .EQ.0) THEN
         CALL FFOPT1(X(LFFATMNAM),X(LFFCORD),X(LFFCORDSV),X(LFFZANF),
     *             X(LFFCHARG),X(LFFPOL),X(LFFPOLSV),X(LFFDIP),
     *             X(LFFFIELD1),X(LFFFIELD2),X(LFFFIELD3),
     *             X(LFFSIG),X(LFFEPS),X(LFFSIG2),X(LFFEPS2),
     *             X(LFFBOND0),X(LFFFCBOND),X(LFFFCSTBD),
     *             X(LFFANGL0),X(LFFFCANGL),X(LFFFCWAGG),
     *             X(LFFDIHB0),X(LFFFCDIHB),X(LFFFCDIHR),
     *             X(LFFVROT),X(LFFNNN),X(LFFGAMA),X(LFFIPAIR),
     *             X(LFFKLIST),X(LFFLLIST),X(LFFKBLST),
     *             X(LFFMLIST),X(LFFNLIST),
     *             X(LFFL1213J),X(LFFL14J),X(LFFFFGRD0),
     *             X(LFFFFGRD2),X(LFFFFGRD1),X(LFFFCLJTP),X(LFFNTYPE),
     *             X(LFFXTS),X(LFFYTS),X(LFFZTS),X(LFFCMAT1),
     *             X(LFFPOT1),X(LFFPOT2),X(LFFQRXN1),X(LFFQRXN2),NTS,
     *             X(LFFNONLS1),X(LFFNONLSTQ),
     *             X(LFFMAPLST),X(LFFCMAPCO),
     *             X(LFFLSTCELL),
     *             X(LFFNONLS2),X(LFFCORDSV2),X(LFFCORDSVQ),
     *             X(LFFMVFASTS2),X(LFFMVFASTS3),X(LFFMVFASTS4),
     *             X(LFFMVFASTL2),X(LFFMVFASTL3),X(LFFMVFASTL4),
     *             X(LFFAFIX),X(LFFQFIX),
     *             X(LFFRFIX),X(LFFIDATOM),X(LFFDAI),X(LFFIDDAI),
     *             X(LFFVFIX1),X(LFFVFIX2),
     *             X(LFFXTSFIX),X(LFFYTSFIX),X(LFFZTSFIX),
     *             X(LFFNONLSA),X(LFFNONLSB),X(LFFL1213A),
     *             X(LFFL1213B),X(LFFL14A),X(LFFL14B),
     *             X(LFF2CHARG),
     *             X(LFF2SIG),X(LFF2EPS),X(LFF2SIG2),X(LFF2EPS2),
     *             X(LFFNONLSPMA),X(LFFL1213PMA),X(LFFL14PMA),
     *             X(LFF2CORD),X(LFFLSBONDPMA),X(LFFLSANGLPMA),
     *             X(LFFLSDIHRPMA),X(LFFLSDIHBPMA),
     *             X(LFFLSWAGGPMA),X(LFFLSCMAPPMA),
     *             X(LFFNONLSPMB),X(LFFL1213PMB),X(LFFL14PMB),
     *             X(LFFLSBONDPMB),X(LFFLSANGLPMB),
     *             X(LFFLSDIHRPMB),X(LFFLSDIHBPMB),
     *             X(LFFLSWAGGPMB),X(LFFLSCMAPPMB),
     *             X(LFFUMBHIS),X(LFFUM2HIS))
      ELSE IF (NAT.GT.0 .AND. NFFAT.GT.0 .AND. IHESS .EQ.0) THEN
         CALL FFOPT2(X(LFFATMNAM),X(LFFCORD),X(LFFCORDSV),X(LFFZANF),
     *             X(LFFFFGRD0),
     *             X(LFFFFGRD2),X(LFFFFGRD1),X(LFFQMGRD2),X(LFFQMGRD1),
     *             X(LFFQMGRD0),
     *             X(LFFLISTQM),X(LFFNONLS1),X(LFFNONLSTQ),
     *             X(LFFLSTCELL),X(LFFNONLS2),
     *             X(LFFCORDSV2),X(LFFCORDSVQ),
     *             X(LFFMVFASTS2),X(LFFMVFASTS3),X(LFFMVFASTS4),
     *             X(LFFMVFASTL2),X(LFFMVFASTL3),X(LFFMVFASTL4),
     *             X(LFFNONLSA),X(LFFNONLSB),
     *             X(LFFNONLSPMA),X(LFFNONLSPMB))
      ELSE IF (NAT.EQ.0 .AND. NFFAT.GT.0 .AND. IHESS .GT.0) THEN
         CALL VALFM(LOADFM)
         LXOLD   = LOADFM + 1
         LXNEW   = LXOLD  + 3*NFFAT
         LAST    = LXNEW  + 3*NFFAT
         NEED    = LAST   - LOADFM -1
         CALL GETFM(NEED)
         CALL FFOPT3(X(LFFATMNAM),X(LFFCORD),X(LFFCORDSV),X(LFFZANF),
     *             X(LFFZMAS),
     *             X(LFFCHARG),X(LFFPOL),X(LFFPOLSV),X(LFFDIP),
     *             X(LFFFIELD1),X(LFFFIELD2),X(LFFFIELD3),
     *             X(LFFSIG),X(LFFEPS),X(LFFSIG2),X(LFFEPS2),
     *             X(LFFBOND0),X(LFFFCBOND),X(LFFFCSTBD),
     *             X(LFFANGL0),X(LFFFCANGL),X(LFFFCWAGG),
     *             X(LFFDIHB0),X(LFFFCDIHB),X(LFFFCDIHR),
     *             X(LFFVROT),X(LFFNNN),X(LFFGAMA),X(LFFIPAIR),
     *             X(LFFKLIST),X(LFFLLIST),X(LFFKBLST),
     *             X(LFFMLIST),X(LFFNLIST),
     *             X(LFFL1213J),X(LFFL14J),
     *             X(LFFFFGRD2),X(LFFFCLJTP),X(LFFNTYPE),
     *             X(LFFXTS),X(LFFYTS),X(LFFZTS),X(LFFCMAT1),
     *             X(LFFPOT1),X(LFFPOT2),X(LFFQRXN1),X(LFFQRXN2),NTS,
     *             X(LFFNONLS1),X(LFFNONLSTQ),
     *             X(LFFMAPLST),X(LFFCMAPCO),
     *             X(LXNEW),X(LXOLD),
     *             X(LFFLSTCELL),
     *             X(LFFNONLS2),X(LFFCORDSV2),X(LFFCORDSVQ),
     *             X(LFFMVFASTS2),X(LFFMVFASTS3),X(LFFMVFASTS4),
     *             X(LFFMVFASTL2),X(LFFMVFASTL3),X(LFFMVFASTL4),
     *             X(LFFAFIX),X(LFFQFIX),
     *             X(LFFRFIX),X(LFFIDATOM),X(LFFDAI),X(LFFIDDAI),
     *             X(LFFVFIX1),X(LFFVFIX2),
     *             X(LFFXTSFIX),X(LFFYTSFIX),X(LFFZTSFIX),
     *             X(LFFNONLSA),X(LFFNONLSB),X(LFFL1213A),
     *             X(LFFL1213B),X(LFFL14A),X(LFFL14B),
     *             X(LFF2CHARG),
     *             X(LFF2SIG),X(LFF2EPS),X(LFF2SIG2),X(LFF2EPS2),
     *             X(LFFNONLSPMA),X(LFFL1213PMA),X(LFFL14PMA),
     *             X(LFF2CORD),X(LFFLSBONDPMA),X(LFFLSANGLPMA),
     *             X(LFFLSDIHRPMA),X(LFFLSDIHBPMA),
     *             X(LFFLSWAGGPMA),X(LFFLSCMAPPMA),
     *             X(LFFNONLSPMB),X(LFFL1213PMB),X(LFFL14PMB),
     *             X(LFFLSBONDPMB),X(LFFLSANGLPMB),
     *             X(LFFLSDIHRPMB),X(LFFLSDIHBPMB),
     *             X(LFFLSWAGGPMB),X(LFFLSCMAPPMB),
     *             X(LFFUMBHIS),X(LFFUM2HIS))
         CALL RETFM(NEED)
      ELSE IF (NAT.GT.0 .AND. NFFAT.GT.0 .AND. IHESS .GT.0) THEN
         CALL VALFM(LOADFM)
         LXOLD   = LOADFM + 1
         LXNEW   = LXOLD  + 3*(NFFAT+NAT)
         LAST    = LXNEW  + 3*(NFFAT+NAT)
         NEED    = LAST   - LOADFM -1
         CALL GETFM(NEED)
         CALL FFOPT4(X(LFFATMNAM),X(LFFCORD),X(LFFCORDSV),X(LFFZANF),
     *             X(LFFZMAS),
     *             X(LFFFFGRD2),X(LFFQMGRD2),
     *             X(LFFLISTQM),X(LFFNONLS1),X(LFFNONLSTQ),
     *             X(LXNEW),X(LXOLD),
     *             X(LFFLSTCELL),X(LFFNONLS2),
     *             X(LFFCORDSV2),X(LFFCORDSVQ),
     *             X(LFFMVFASTS2),X(LFFMVFASTS3),X(LFFMVFASTS4),
     *             X(LFFMVFASTL2),X(LFFMVFASTL3),X(LFFMVFASTL4),
     *             X(LFFNONLSA),X(LFFNONLSB),
     *             X(LFFNONLSPMA),X(LFFNONLSPMB))
         CALL RETFM(NEED)
      END IF
      END IF
C
      IF(RUNTYP.EQ.HSSIAN) THEN
      IF(NUMBTYP.GT.0.OR.NUM2TYP.GT.0) THEN
         WRITE(IW,'(/1X,A,A/)')'ERROR: JUMBPOT OR JUM2POT ',
     *   'SHOULD NOT BE USED FOR RUNTYP=HESSIAN.'
         CALL ABRT
      END IF
      IF(NAT.EQ.0) THEN
         NATVIB  = NFFAT
         IF(NACTMM.GT.0) NATVIB = NACTMM
         CALL VALFM(LOADFM)
         LVIBGRD = LOADFM + 1
         LVIBHSS = LVIBGRD+ 3*NATVIB*3*NATVIB*2
         LVIBDIP = LVIBHSS+ 3*NATVIB*3*NATVIB
         LVIBDDM = LVIBDIP+ 3*3*NATVIB*2
         LAST    = LVIBDDM+ 3*3*NATVIB
         NEED    = LAST   - LOADFM -1
         CALL GETFM(NEED)
         CALL HESSMM(X(LFFATMNAM),X(LFFCORD),X(LFFCORDSV),X(LFFZANF),
     *             X(LFFZMAS),
     *             X(LFFCHARG),X(LFFPOL),X(LFFPOLSV),X(LFFDIP),
     *             X(LFFFIELD1),X(LFFFIELD2),X(LFFFIELD3),
     *             X(LFFSIG),X(LFFEPS),X(LFFSIG2),X(LFFEPS2),
     *             X(LFFBOND0),X(LFFFCBOND),X(LFFFCSTBD),
     *             X(LFFANGL0),X(LFFFCANGL),X(LFFFCWAGG),
     *             X(LFFDIHB0),X(LFFFCDIHB),X(LFFFCDIHR),
     *             X(LFFVROT),X(LFFNNN),X(LFFGAMA),X(LFFIPAIR),
     *             X(LFFKLIST),X(LFFLLIST),X(LFFKBLST),
     *             X(LFFMLIST),X(LFFNLIST),
     *             X(LFFL1213J),X(LFFL14J),
     *             X(LFFFFGRD2),X(LFFFCLJTP),X(LFFNTYPE),
     *             X(LFFXTS),X(LFFYTS),X(LFFZTS),X(LFFCMAT1),
     *             X(LFFPOT1),X(LFFPOT2),X(LFFQRXN1),X(LFFQRXN2),NTS,
     *             X(LFFNONLS1),X(LFFNONLSTQ),
     *             X(LFFMAPLST),X(LFFCMAPCO),
     *             X(LFFLSTCELL),
     *             X(LFFNONLS2),X(LFFCORDSV2),X(LFFCORDSVQ),
     *             X(LFFMVFASTS2),X(LFFMVFASTS3),X(LFFMVFASTS4),
     *             X(LFFMVFASTL2),X(LFFMVFASTL3),X(LFFMVFASTL4),
     *             X(LFFAFIX),X(LFFQFIX),
     *             X(LFFRFIX),X(LFFIDATOM),X(LFFDAI),X(LFFIDDAI),
     *             X(LFFVFIX1),X(LFFVFIX2),
     *             X(LFFXTSFIX),X(LFFYTSFIX),X(LFFZTSFIX),
     *             X(LFFNONLSA),X(LFFNONLSB),X(LFFL1213A),
     *             X(LFFL1213B),X(LFFL14A),X(LFFL14B),
     *             X(LFF2CHARG),
     *             X(LFF2SIG),X(LFF2EPS),X(LFF2SIG2),X(LFF2EPS2),
     *             X(LFFNONLSPMA),X(LFFL1213PMA),X(LFFL14PMA),
     *             X(LFF2CORD),X(LFFLSBONDPMA),X(LFFLSANGLPMA),
     *             X(LFFLSDIHRPMA),X(LFFLSDIHBPMA),
     *             X(LFFLSWAGGPMA),X(LFFLSCMAPPMA),
     *             X(LFFNONLSPMB),X(LFFL1213PMB),X(LFFL14PMB),
     *             X(LFFLSBONDPMB),X(LFFLSANGLPMB),
     *             X(LFFLSDIHRPMB),X(LFFLSDIHBPMB),
     *             X(LFFLSWAGGPMB),X(LFFLSCMAPPMB),
     *             X(LFFUMBHIS),X(LFFUM2HIS),
     *             X(LVIBGRD),X(LVIBDIP),X(LVIBHSS),X(LVIBDDM),
     *             NATVIB)
         CALL RETFM(NEED)
      ELSE IF(NAT.GT.0.AND.NACTMM.EQ.0) THEN
         CALL VALFM(LOADFM)
         LVIBGRD = LOADFM + 1
         LVIBHSS = LVIBGRD+ 3*NAT*3*NAT*2
         LVIBDIP = LVIBHSS+ 3*NAT*3*NAT
         LVIBE   = LVIBDIP+ 3*3*NAT*2
         LVIBDDM = LVIBE  + 3*NAT*2
         LAST    = LVIBDDM+ 3*3*NAT
         NEED    = LAST   - LOADFM -1
         CALL GETFM(NEED)
         CALL HESSQM(X(LFFCORD),X(LFFCORDSV),X(LFFDIP),
     *             X(LFFFFGRD2),X(LFFQMGRD2),
     *             X(LFFLISTQM),X(LFFNONLS1),X(LFFNONLSTQ),
     *             X(LFFLSTCELL),X(LFFNONLS2),
     *             X(LFFCORDSV2),X(LFFCORDSVQ),
     *             X(LFFMVFASTS2),X(LFFMVFASTS3),X(LFFMVFASTS4),
     *             X(LFFMVFASTL2),X(LFFMVFASTL3),X(LFFMVFASTL4),
     *             X(LFFNONLSA),X(LFFNONLSB),
     *             X(LFFNONLSPMA),X(LFFNONLSPMB),
     *             X(LVIBGRD),X(LVIBDIP),X(LVIBE),X(LVIBHSS),
     *             X(LVIBDDM))
         CALL RETFM(NEED)
      ELSE IF(NAT.GT.0.AND.NACTMM.GT.2) THEN
         CALL VALFM(LOADFM)
         LVIBGRD = LOADFM   + 1
         LVIBE   = LVIBGRD  + 3*NACTMM*3*NACTMM*2
         LVIBGRD0= LVIBE    + 3*MXATM*2
         LAST    = LVIBGRD0 + 3*NACTMM
         NEED    = LAST     - LOADFM -1
         CALL GETFM(NEED)
         CALL HESSQMMM(X(LFFCORD),X(LFFCORDSV),
     *             X(LFFFFGRD2),X(LFFQMGRD2),
     *             X(LFFLISTQM),X(LFFNONLS1),X(LFFNONLSTQ),
     *             X(LFFLSTCELL),X(LFFNONLS2),
     *             X(LFFCORDSV2),X(LFFCORDSVQ),
     *             X(LFFMVFASTS2),X(LFFMVFASTS3),X(LFFMVFASTS4),
     *             X(LFFMVFASTL2),X(LFFMVFASTL3),X(LFFMVFASTL4),
     *             X(LFFNONLSA),X(LFFNONLSB),
     *             X(LFFNONLSPMA),X(LFFNONLSPMB),
     *             X(LVIBGRD),X(LVIBE),X(LVIBGRD0))
         CALL RETFM(NEED)
      END IF
      END IF
C
      CALL TIMIT(1)
C
      RETURN
      END
C*MODULE QUANPOC  *DECK FFOPT1
!>
!> @brief    pure MM optimization
!>
!> @author   Hui Li
!>           - Apr 2011
!>
!> @details  relaxed steepest descent method
!>
      SUBROUTINE FFOPT1(ATMNAM,CORD,CORDSV,ZANF,
     *                  CHARG,POL,POLSV,DIP,
     *                  FIELD1,FIELD2,FIELD3,
     *                  SIG,EPS,SIG2,EPS2,
     *                  BOND0,FCBOND,FCSTBD,
     *                  ANGL0,FCANGL,FCWAGG,
     *                  DIHB0,FCDIHB,FCDIHR,
     *                  VROT,NNN,GAMA,IPAIR,
     *                  KLIST,LLIST,KBLST,MLIST,NLIST,
     *                  L1213J,L14J,FFGRD0,
     *                  FFGRD2,FFGRD1,FCLJTP,NTYPE,
     *                  XTS,YTS,ZTS,CMAT1,
     *                  POT1,POT2,QRXN1,QRXN2,NTS,
     *                  NONLS1,NONLSTQ,MAPLST,CMAPCO,
     *                  LSTCELL,NONLS2,
     *                  CORDSV2,CORDSVQ,
     *                  MVFASTS2,MVFASTS3,MVFASTS4,
     *                  MVFASTL2,MVFASTL3,MVFASTL4,
     *                  AFIX,QFIX,
     *                  RFIX,IDATOM,DAI,IDDAI,
     *                  VFIX1,VFIX2,XTSFIX,YTSFIX,ZTSFIX,
     *                  NONLSA,NONLSB,L1213A,L1213B,L14A,L14B,
     *                  CHARGB,
     *                  SIGB,EPSB,SIG2B,EPS2B,
     *                  NONLSPMA,
     *                  L1213PMA,L14PMA,CORDB,
     *                  LSBONDPMA,LSANGLPMA,
     *                  LSDIHRPMA,LSDIHBPMA,
     *                  LSWAGGPMA,LSCMAPPMA,
     *                  NONLSPMB,L1213PMB,L14PMB,
     *                  LSBONDPMB,LSANGLPMB,
     *                  LSDIHRPMB,LSDIHBPMB,
     *                  LSWAGGPMB,LSCMAPPMB,UMBHIS,UM2HIS)
      use mx_limits, only: mxatm
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      LOGICAL GOPARR,DSKWRK,MASWRK
      LOGICAL LINEAR,OUT,STPT,PROJCT,GOTEG
C
      PARAMETER (TOKCAL=627.509469D+00)
      PARAMETER (TOANGS=0.52917724924D+00)
      PARAMETER (TOBOHR=1.0D+00/TOANGS)
      PARAMETER (ZERO=0.0D+00)
      PARAMETER (THREE=3.0D+00)
C
      CHARACTER*10  ATMNAM
C
      DOUBLE PRECISION METHOD
C
      DIMENSION ATMNAM(NFFAT),CORD(3,NFFAT),CORDSV(3,NFFAT),ZANF(NFFAT),
     *          CHARG(NFFAT),
     *          POL(NFFAT),DIP(3,NFFAT),
     *          FIELD1(3,NFFAT),FIELD2(3,NFFAT),FIELD3(3,NFFAT),
     *          SIG(NFFAT),EPS(NFFAT),SIG2(NFFAT),EPS2(NFFAT),
     *          BOND0(NBOND),FCBOND(NBOND),
     *          ANGL0(NANGL),FCANGL(NANGL),
     *          FCWAGG(NWAGG),FCDIHB(NDIHB),FCDIHR(3,*),
     *          VROT(NDIHR),NNN(NDIHR),GAMA(NDIHR),IPAIR(2,NBOND),
     *          KLIST(3,NANGL),
     *          LLIST(4,NDIHR),MLIST(4,NWAGG),
     *          L1213J(2,*),L14J(2,NDIHR),FFGRD0(3,NFFAT),
     *          FFGRD2(3,NFFAT),FFGRD1(3,NFFAT),NTYPE(*),
     *          XTS(NTS),YTS(NTS),ZTS(NTS),CMAT1(NTS,NTS),
     *          POT1(NTS),POT2(NTS),QRXN1(NTS),QRXN2(NTS),NONLS1(2,*),
     *          NONLSTQ(*),MAPLST(6,*),CMAPCO(4,4,24,24,3)
      DIMENSION TIMSTR(3)
      DIMENSION ENALL(100)
C
      COMMON /FFENGY/ EN12,EN123,EN123R4,EN123B4,EN234W1,ENCHAR,ENLJR,
     *                ENLJD,ENPOL,XENPOL,ENRXN,XENRXN,ENRXNPOL,ENRXNR,
     *                EN12312,ENQUANP(30),
     *                ENBIAS,ENCENT,ENUCCH,ENCMAP,ENPOT,ENKIN,ENTOT
      COMMON /FFFIXS/ ENFIXSO,FIXEPS,FIXTOL,FIXA,FIXQ,RALLMM,RALLQM,
     *                RADMM(200),RADQM(200),NRADMM,NRADQM,IFIXSOL,
     *                LFFDAI,LFFDAIT,LFFIDDAI,LFFIDTMP,LFFTMPTS,
     *                LFFAFIX,LFFIDATOM,LFFRFIX,LFFQFIX,NTSATM,
     *                LFFQFIXMP,LFFQFIXTA,LFFQFIXXY,
     *                LFFXTSFIX,LFFYTSFIX,LFFZTSFIX,
     *                LFFVFIX1,LFFVFIX2,NCYCLE,MXFFTS,NFFTS
      COMMON /FFFREE/ SOL1CH,SOL1LJ,SOL1IM,SOLFRE1MM,
     *                SOL2CH,SOL2LJ,SOL2IM,SOLFRE2MM,
     *                ASOL1CH,ASOL1LJ,ASOL1IM,ASOL1MM,
     *                ASOL2CH,ASOL2LJ,ASOL2IM,ASOL2MM,
     *                PMF1BD,PMF1AG,PMF1DR,PMF1DB,PMF1WG,PMF1CM,
     *                PMF1CH,PMF1PO,PMF1LJ,PMF1IM,PMFFRE1MM,
     *                EN12A,EN12B,EN123A,EN123B,EN123R4A,EN123R4B,
     *                EN234W1A,EN234W1B,EN123B4A,EN123B4B,
     *                EN12312A,EN12312B,
     *                ENCHARA,ENCHARB,ENPOLA,ENPOLB,
     *                ENRXNA,ENRXNB,ENRXNPOLA,ENRXNPOLB,
     *                ENRXNRA,ENRXNRB,ENLJRA,ENLJRB,ENLJDA,ENLJDB,
     *                ENCMAPA,ENCMAPB,ENFIXSOA,ENFIXSOB,
     *                ENCENTA,ENCENTB,WSIMUL,WPERT1,WPERT2,
     *                IFEPTYP
      COMMON /FFPARM/ NFFAT,NBOND,NANGL,NDIHR,NDIHB,NCMAP,NWAGG,
     *                N1213J,N14J,NLKQMM,IDOCHG,IDOPOL,IDOLJ,IDOCMAP
      COMMON /FFMDOP/ LFFCORDG,LFFCORDGSV,LFFCORDGSV2,LFFCORDGSVQ,
     *                LFFFFGRDG0,LFFFFGRDG1,LFFFFGRDG2,MDOPT
      COMMON /FFMDPA/ DT,DT2,TEMP0,PRES0,POLTOL,VIR(3),PMEAN,VOLAV,
     *                ENPAV,ENKAV,TEMPAV,BERENDT,BERENDP,VELMAX,
     *                PMEANX,PMEANY,PMEANZ,NSTEP,KMASTER,KOUTACT(2),
     *                IHESS,INTALG,ITSTAT,IPSTAT,JOUT,KOUT,LOUT
      COMMON /FFMPT3/ NACTMM,LACTMM(2020),LFFDIPOLD,JUMBUP,
     *                NACTQM,LACTQM(2020),LFFOLDC,LFFQMVELSV,MMHESS,
     *                LFFQMCHG,LFFQMCHGB,ISWAP,R2SWAP,DFTBMM
      COMMON /FFPBSW/ XBOX,YBOX,ZBOX,SWF,SWFDX,SWFDY,SWFDZ,
     *                SWRA,ONESWRA,SWRA2,ONESWRA2,
     *                SWRB,ONESWRB,SWRB2,ONESWRB2,ONESWRB4,
     *                SWFDUM3,SWFDUM4,SWFDUM5,
     *                SWRAQ,ONESWRAQ,SWRAQ2,ONESWRAQ2,
     *                SWRBQ,ONESWRBQ,SWRBQ2,ONESWRBQ2,ONESWRBQ4,
     *                SWFDUM3Q,SWFDUM4Q,SWFDUM5Q,
     *                QMSIZE,QMCX,QMCY,QMCZ,QMCXSV,QMCYSV,QMCZSV,
     *                CENTX,CENTY,CENTZ,BUFWID1,BUFWID2,RDAMP,
     *                EFIELDX,EFIELDY,EFIELDZ,QMCXSV2,QMCYSV2,QMCZSV2,
     *                EPS1RB,EPS1RB3,ONEXBOX,ONEYBOX,ONEZBOX,
     *                LQMCT,MXLIST1,NTODO,NTODOSV,NTODOQ,
     *                ISWITCH,ISHIFT,IPOLSHF,
     *                LFFLSTCELL,LFFCORDSV,
     *                LFFPOLSV,LFFCORDSV2,LFFNONLS2,LFFCORDSVQ,
     *                LFFMVFASTS2,LFFMVFASTS3,LFFMVFASTS4,
     *                LFFMVFASTL2,LFFMVFASTL3,LFFMVFASTL4,
     *                MXCHECK,MXLIST2,NTODO2,NTODO2SV
      COMMON /FFRMSD/ DIPT(8),TIMGYRA,TIMRALL,
     *                NATPDB,NGYRA,NDIEL,NRALL,
     *                LFFRALL0,NRMSD,LFFRMSD0,KFREEAB(201),
     *                NRIJMM,IJRMM(2,100),NRIJQM,IJRQM(2,100),
     *                NAIJKMM,IJKMM(3,100),NAIJKQM,IJKQM(3,100),
     *                NFIXMM,IFIXMM(200),NFIXQM,IFIXQM(200)
      COMMON /FFRXN / RXNEPS,RSPHSOL,ISPHSOL
      COMMON /FFSPH / SPHRAD,SPHSIG,SPHEPS,IADDWAT
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /IOFILE/ IR,IW,IP,IJK,IJKT,IDAF,NAV,IODA(950)
      COMMON /OPTGRD/ XX(3*MXATM),ES,FE(20),
     *                CONVF,FMAXT,DXMAXT,RMAX,RMIN,RLIM,
     *                EIGMAX,EIGMIN,GRDERR,FRMS,FMAX,TRMAX,TRMIN,
     *                IC(20),MSTEP,NSERCH,NPMAX,NP,IFOLOW,
     *                NNEG,IUPHSS,IEXIT,ITRUPD,IPAD,KDIAGH,NPRICO
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /RUNOPT/ RUNTYP,EXETYP,NEVALS,NGLEVL,NHLEVL
      COMMON /ZMAT  / NZMAT,NZVAR,NVAR,NSYMC,LINEAR
C
C     HUI LI, APR 2011, LINCOLN
C
      IF(NSTEP.LT.0) RETURN
C
      CALL SEQOPN(36,'TRAJECT','NEW',.FALSE.,'FORMATTED')
      IF(IADDWAT.EQ.2 .AND. ISPHSOL.GT.0) RSPHSOL=RSPHSOL*1.0D+40
C
      CALL VCLR(FFGRD2,1,3*NFFAT)
      ICONV   = 0
C
C     -- SET SOME CONTROL VARIABLES --
      NAT    = NFFAT
      IF(NFFAT.GE.10) NAT = 10
      NCOORD = 3*NAT
      NCVAL  = NCOORD
      MSTEP  = NSTEP
      NZMAT  = 0
      MODE   = 0
      OUT    = .FALSE.
      MAX2   = NCVAL
C
      IF(MDOPT.EQ.0) THEN
C     -- THIS IS MAINLY TO GET CONVF FROM THE INPUT FILE
C        TO AVOID ZERO HESSIAN, USE SOME FAKE QM COORDINATES
      CALL DCOPY(3*NAT,CORD,1,C,1)
      CALL SIGINI(MODE,RUNTYP,NCVAL,NCOORD,METHOD,OUT,
     *            GOTEG,NPRT,NPUN,ITBMAT,STPT,STSTEP,
     *            PROJCT,DUMMY,DUMMY,MAX2,NPRTHS)
      ELSE
      CONVF=1.0D-05
      END IF
      NAT    = 0
C
      DO ISTEP = 0, NSTEP
         IF(ISTEP.EQ.100.AND.IADDWAT.EQ.2 .AND. ISPHSOL.GT.0)
     *   RSPHSOL = RSPHSOL*1.0D-40
C
C        -- CALCULATE ENERGY AND GRADIENT
         CALL NONBOND(ISTEP,CORD,CORDSV,CORDSV2,CORDSVQ,
     *                NONLS1,NONLS2,
     *                NONLSTQ,LSTCELL,
     *                MVFASTS2,MVFASTS3,MVFASTS4,
     *                MVFASTL2,MVFASTL3,MVFASTL4,
     *                NONLSA,NONLSB,
     *                NONLSPMA,NONLSPMB)
         CALL VCLR(FFGRD2,1,3*NFFAT)
         CALL VCLR(VIR,1,3)    !  NEED TO DO THIS - HUI LI
         CALL E00012(CORD,FFGRD2,BOND0,FCBOND,IPAIR,CORDB,
     *               LSBONDPMA,LSBONDPMB)
         CALL E00123(CORD,FFGRD2,ANGL0,FCANGL,KLIST,CORDB,
     *               LSANGLPMA,LSANGLPMB)
         CALL E12312(CORD,FFGRD2,ANGL0,KLIST,BOND0,FCSTBD,
     *               KBLST,CORDB,LSANGLPMA,LSANGLPMB)
         CALL E123B4(CORD,FFGRD2,DIHB0,FCDIHB,NLIST,CORDB,
     *               LSDIHBPMA,LSDIHBPMB)
         CALL E234W1(CORD,FFGRD2,      FCWAGG,MLIST,CORDB,
     *               LSWAGGPMA,LSWAGGPMB)
         CALL E123R4(CORD,FFGRD2,VROT,GAMA,NNN,LLIST,CORDB,
     *               LSDIHRPMA,LSDIHRPMB,FCDIHR)
         CALL ECMAP (CORD,FFGRD2,MAPLST,CMAPCO,CORDB,
     *               LSCMAPPMA,LSCMAPPMB)
         CALL ELJ126(CORD,FFGRD2,SIG,EPS,SIG2,EPS2,L14J,NONLS1,
     *               L1213J,SIGB,EPSB,SIG2B,EPS2B,
     *               NONLSA,NONLSB,L1213A,L1213B,L14A,L14B,
     *               NONLSPMA,L1213PMA,L14PMA,CORDB,
     *               NONLSPMB,L1213PMB,L14PMB,FCLJTP,NTYPE)
         CALL ESPHER(CORD,FFGRD2)
         CALL ECHARG(CORD,FFGRD2,CHARG,NONLS1,L1213J,L14J,
     *               CHARGB,NONLSA,NONLSB,L1213A,L1213B,L14A,L14B,
     *               NONLSPMA,L1213PMA,L14PMA,CORDB,
     *               NONLSPMB,L1213PMB,L14PMB)
         CALL UMBRELLA(CORD,FFGRD2,UMBHIS,UM2HIS)
         IF(IFIXSOL.EQ.0) THEN
            CALL CHGRXN(CORD,FFGRD2,CHARG,XTS,YTS,ZTS,
     *                  CMAT1,POT1,QRXN1,NTS)
            CALL POLRXN(CORD,FFGRD2,CHARG,POL,POLSV,DIP,
     *                  FIELD1,FIELD2,FIELD3,
     *                  XTS,YTS,ZTS,CMAT1,POT1,POT2,QRXN1,QRXN2,NTS,
     *                  NONLS1,L1213J)
         END IF
         IF(IFIXSOL.EQ.1) THEN
            CALL FIXSOL(CORD,FFGRD2,CHARG,ZANF,AFIX,QFIX,
     *                  VFIX1,VFIX2,XTSFIX,YTSFIX,ZTSFIX,
     *                  RFIX,IDATOM,DAI,IDDAI,
     *                  POL,POLSV,DIP,FIELD1,FIELD2,FIELD3,
     *                  NONLS1,L1213J)
         END IF
         IF(GOPARR) THEN
            CALL VCLR(ENALL,1,100)
            ENALL( 1) = EN12
            ENALL( 2) = EN123
            ENALL( 3) = EN123R4
            ENALL( 4) = EN123B4
            ENALL( 5) = EN234W1
            ENALL( 6) = ENCHAR
            ENALL( 7) = ENLJR
            ENALL( 8) = ENLJD
C           ENALL(12) = ENUCCH
            ENALL(13) = ENRXN
            ENALL(14) = ENRXNR
C           ENALL(15) = ENCENT
            ENALL(16) = ENCMAP
            ENALL(18) = SOL1CH
            ENALL(19) = SOL1LJ
            ENALL(21) = SOL2CH
            ENALL(22) = SOL2LJ
            ENALL(33) = PMF1BD
            ENALL(34) = PMF1AG
            ENALL(35) = PMF1DR
            ENALL(36) = PMF1DB
            ENALL(37) = PMF1WG
            ENALL(38) = PMF1CM
            ENALL(39) = PMF1CH
            ENALL(40) = PMF1LJ
            ENALL(41) = ENBIAS
            ENALL(42) = EN12312
            CALL DDI_GSUMF(2410,ENALL  ,42)
            CALL DDI_GSUMF(2411,FFGRD2,3*NFFAT)
            EN12      = ENALL( 1)
            EN123     = ENALL( 2)
            EN123R4   = ENALL( 3)
            EN123B4   = ENALL( 4)
            EN234W1   = ENALL( 5)
            ENCHAR    = ENALL( 6)
            ENLJR     = ENALL( 7)
            ENLJD     = ENALL( 8)
C           ENUCCH    = ENALL(12)
            ENRXN     = ENALL(13)
            ENRXNR    = ENALL(14)
C           ENCENT    = ENALL(15)
            ENCMAP    = ENALL(16)
            SOL1CH    = ENALL(18)
            SOL1LJ    = ENALL(19)
            SOL2CH    = ENALL(21)
            SOL2LJ    = ENALL(22)
            PMF1BD    = ENALL(33)
            PMF1AG    = ENALL(34)
            PMF1DR    = ENALL(35)
            PMF1DB    = ENALL(36)
            PMF1WG    = ENALL(37)
            PMF1CM    = ENALL(38)
            PMF1CH    = ENALL(39)
            PMF1LJ    = ENALL(40)
            ENBIAS    = ENALL(41)
            EN12312   = ENALL(42)
         END IF
C
C        - ZERO OFF SOME FORCES -
C
         IF(NACTMM.GT.0) THEN
            DO KOPT = 1, NACTMM
               IFFAT = LACTMM(KOPT)
               IF(IFFAT.GT.0)FFGRD2(1,IFFAT)=FFGRD2(1,IFFAT)+1.0D+03
            ENDDO
            DO IFFAT = 1, NFFAT
               IF(FFGRD2(1,IFFAT).GT.0.5D+03) THEN
                  FFGRD2(1,IFFAT) = FFGRD2(1,IFFAT)-1.0D+03
               ELSE
                  DO III = 1, 3
                     FFGRD2(III,IFFAT) = ZERO
                  ENDDO
               END IF
            ENDDO
         END IF
         DO KFIX=1,NFIXMM
            DO III = 1, 3
               FFGRD2(III,IFIXMM(KFIX)) = ZERO
            ENDDO
         ENDDO
C
C        - TEST CONVERGENCE
C
         GRDMAX = ZERO
         GRDRMS = ZERO
         DO IFFAT=1, NFFAT
            GRDMAX=MAX(GRDMAX,ABS(FFGRD2(1,IFFAT)))
            GRDMAX=MAX(GRDMAX,ABS(FFGRD2(2,IFFAT)))
            GRDMAX=MAX(GRDMAX,ABS(FFGRD2(3,IFFAT)))
            GRDRMS=GRDRMS + FFGRD2(1,IFFAT)**2
     *                    + FFGRD2(2,IFFAT)**2
     *                    + FFGRD2(3,IFFAT)**2
         ENDDO
         IF(NACTMM.EQ.0) GRDRMS = SQRT(GRDRMS/(3*(NFFAT-NFIXMM)))
         IF(NACTMM.GT.0) GRDRMS = SQRT(GRDRMS/(3*(NACTMM-NFIXMM)))
         IF(GRDMAX.LT.CONVF.AND.GRDRMS.LT.(CONVF/THREE)) ICONV = 1
C
C        -- CALCULATE PROPERTIES
         IIISTEP=ISTEP
         CALL OPTPROP(IIISTEP,ICONV)
         IF(MDOPT.EQ.0) THEN
         IF(MASWRK) WRITE(IW,'(A,I6,A,F12.9,A,F12.9,A,F24.10,A)')
     *              ' OPT STEP=',ISTEP,
     *              '  MAX GRAD=',GRDMAX,
     *              '  RMS GRAD=',GRDRMS,
     *              '  E=',ENTOT*TOKCAL,' KCAL/MOL'
         END IF
C
         IF(MDOPT.GT.0.AND.(ISTEP.EQ.NSTEP.OR.ICONV.EQ.1)) THEN
         IF(MASWRK) WRITE(IW,'(A,I6,A,F24.10,A)')
     *              ' MDOPT STEP=',ISTEP,
     *              '  E=',ENTOT*TOKCAL,' KCAL/MOL'
         END IF
         CALL TMDATE(TIMSTR)
         IF(MOD(ISTEP,KOUT).EQ.0) CALL TIMIT(1)
         IF(MASWRK .AND. (MOD(ISTEP,KOUT).EQ.0.OR.ISTEP.EQ.NSTEP.OR.
     *      ICONV.EQ.1).AND.MDOPT.EQ.0)THEN
C
            WRITE(36,'(1X,A,I10)')
     *      'RESTART COORDINATES FOR QUANPOL AT OPT STEP',ISTEP
            WRITE(36,'(A,3(A,F15.10,1X),A)')
     *      ' $QUANPO ',
     *      'CENTX=',CENTX*TOANGS,
     *      'CENTY=',CENTY*TOANGS,
     *      'CENTZ=',CENTZ*TOANGS,
     *      '$END'
C
            WRITE(36,'(1X,A,I10,10X,A,F19.2,A)')
     *      '$FFDATA          !   OPT STEP',ISTEP
            WRITE(36,*)'COORDINATES  NUC                   X',
     *          '                   Y                   Z'
            DO IFFAT = 1, NFFAT
               WRITE(36,1000)ATMNAM(IFFAT),ZANF(IFFAT),
     *         CORD(1,IFFAT)*TOANGS,
     *         CORD(2,IFFAT)*TOANGS,
     *         CORD(3,IFFAT)*TOANGS
            ENDDO
            WRITE(36,*)'STOP'
C
            IF(NFFAT.LE.100) THEN
            WRITE(IW,'(1X,A,I10,10X,A,F19.2,A)')
     *      '$FFDATA          !   OPT STEP',ISTEP
            WRITE(IW,*)'COORDINATES  NUC                   X',
     *          '                   Y                   Z'
            DO IFFAT = 1, NFFAT
               WRITE(IW,1000)ATMNAM(IFFAT),ZANF(IFFAT),
     *         CORD(1,IFFAT)*TOANGS,
     *         CORD(2,IFFAT)*TOANGS,
     *         CORD(3,IFFAT)*TOANGS
            ENDDO
            WRITE(IW,*)'STOP'
            END IF
C
C           -- PRINT KOUTACT ATOMS --
C              HERE ONLY $FFDATA IS CONSIDERED, $FFDATB IS NOT.
C
            IF(KOUTACT(1).GT.0) THEN
            WRITE(IW,'(/A,I10)')'!KOUTACT ATOMS AROUND $FFDATA ATOM ',
     *      KOUTACT(1)
            WRITE(IW,'(1X,A,I10)')
     *      '$FFDATA ! KOUTACT     OPT STEP',ISTEP
            ACTX = CORD(1,KOUTACT(1))
            ACTY = CORD(2,KOUTACT(1))
            ACTZ = CORD(3,KOUTACT(1))
            ACTR = DBLE(KOUTACT(2))*1.0D-08*TOBOHR
            ACTR2= ACTR**2
            DO IFFAT=1,NFFAT
               CX    = CORD(1,IFFAT) - ACTX
               CY    = CORD(2,IFFAT) - ACTY
               CZ    = CORD(3,IFFAT) - ACTZ
               PBCX  = XBOX*ANINT(CX*ONEXBOX)
               PBCY  = YBOX*ANINT(CY*ONEYBOX)
               PBCZ  = ZBOX*ANINT(CZ*ONEZBOX)
               CX    = CX - PBCX
               CY    = CY - PBCY
               CZ    = CZ - PBCZ
               RCXYZ = CX**2 + CY**2 + CZ**2
               IF(RCXYZ.LE.ACTR2) THEN
                  WRITE(IW,1000)ATMNAM(IFFAT),ZANF(IFFAT),
     *            (CX+ACTX)*TOANGS,(CY+ACTY)*TOANGS,(CZ+ACTZ)*TOANGS
               END IF
            ENDDO
            WRITE(IW,*)'$END ! KOUTACT'
            END IF
C
            CALL TIMIT(1)
            CALL FLSHBF(IW)
            CALL FLSHBF(36)
         END IF
C
C        -- PRINT FOR MDOPT
C
         IF(MASWRK .AND. MDOPT.GT.0.AND.(ISTEP.EQ.NSTEP.OR.
     *      ICONV.EQ.1))THEN
C
            WRITE(IW,'(1X,A,I10,10X,A,F19.2,A)')
     *      '$FFDATA          !   MDOPT   STEP',ISTEP
            WRITE(IW,*)'COORDINATES  NUC                   X',
     *          '                   Y                   Z'
            DO IFFAT = 1, NFFAT
               WRITE(IW,1000)ATMNAM(IFFAT),ZANF(IFFAT),
     *         CORD(1,IFFAT)*TOANGS,
     *         CORD(2,IFFAT)*TOANGS,
     *         CORD(3,IFFAT)*TOANGS
            ENDDO
            WRITE(IW,*)'STOP'
C
            CALL TIMIT(1)
            CALL FLSHBF(IW)
         END IF
C
         IF(ICONV.EQ.1) GOTO 500
C
C        - CALCULATE NEW COORDINATES
         IF(NSTEP.GT.0)THEN
         DO 200 IFFAT = 1, NFFAT
            IF((FFGRD2(1,IFFAT)+FFGRD2(2,IFFAT)+FFGRD2(3,IFFAT))
     *         .EQ.ZERO) GOTO 200
            DO II = 1, 3
               DISPLC =  -(0.58D+00*FFGRD2(II,IFFAT)+
     *                     0.28D+00*FFGRD1(II,IFFAT)+
     *                     0.14D+00*FFGRD0(II,IFFAT) )*
     *                     ABS(LOG(ABS(CONVF)+ABS(FFGRD2(II,IFFAT))))
     *                     *0.23D+00
               IF(DISPLC.GT. 0.20D+00) DISPLC = 0.20D+00
               IF(DISPLC.LT.-0.20D+00) DISPLC =-0.20D+00
               CORD(II,IFFAT) = CORD(II,IFFAT) + DISPLC
            ENDDO
  200    CONTINUE
         END IF
         CALL DCOPY(3*NFFAT,FFGRD1,1,FFGRD0,1)
         CALL DCOPY(3*NFFAT,FFGRD2,1,FFGRD1,1)
C        -- SYNCHRONIZE CORD EVERY 200 STEPS --
         IF(MOD(ISTEP,200).EQ.0) THEN
            IF(GOPARR) THEN
               CALL DDI_BCAST(462,'F',CORD,3*NFFAT,MASTER)
            END IF
         END IF
C
      ENDDO
C
C     -- NOT LOCATED --
C
      IF(MASWRK) THEN
         WRITE(IW,*)' '
         WRITE(IW,*)'QUANPOL OPTIMIZATION FAILED. ',
     *              'PLEASE RESTART THE JOB USING THE LAST GEOMETRY.'
         WRITE(IW,*)' '
      END IF
      CALL FLSHBF(IW)
      CALL FLSHBF(36)
      IF(MDOPT.EQ.0) CALL SEQCLO(36,'KEEP')
      RETURN
C
C     -- LOCATED --
C
  500 CONTINUE
      IF(MASWRK) THEN
         WRITE(IW,*)' '
         WRITE(IW,*)'================= QUANPOL OPTIMIZATION',
     *              ' SUCCESSFULLY COMPLETED ================='
         WRITE(IW,*)' '
      END IF
C
 1000 FORMAT(1X,A10,1X,F5.1,1X,F19.13,1X,F19.13,1X,F19.13)
C
      CALL FLSHBF(IW)
      CALL FLSHBF(36)
      IF(MDOPT.EQ.0) CALL SEQCLO(36,'KEEP')
C
      RETURN
      END
C*MODULE QUANPOC  *DECK FFOPT2
!>
!> @brief    QM/MM optimization
!>
!> @author   Hui Li
!>           - Apr 2011
!>
!> @details  relaxed steepest descent method
!>
      SUBROUTINE FFOPT2(ATMNAM,CORD,CORDSV,ZANF,FFGRD0,
     *                  FFGRD2,FFGRD1,QMGRD2,QMGRD1,QMGRD0,
     *                  LISTQM,NONLS1,NONLSTQ,
     *                  LSTCELL,NONLS2,CORDSV2,CORDSVQ,
     *                  MVFASTS2,MVFASTS3,MVFASTS4,
     *                  MVFASTL2,MVFASTL3,MVFASTL4,
     *                  NONLSA,NONLSB,
     *                  NONLSPMA,NONLSPMB)
      use mx_limits, only: mxatm,mxao,mxrt
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      LOGICAL GOPARR,DSKWRK,MASWRK
      LOGICAL LINEAR,OUT,STPT,PROJCT,GOTEG
C
      PARAMETER (TOKCAL=627.509469D+00)
      PARAMETER (TOANGS=0.52917724924D+00)
      PARAMETER (TOBOHR=1.0D+00/TOANGS)
      PARAMETER (ZERO=0.0D+00)
      PARAMETER (PT5=0.5D+00)
      PARAMETER (THREE=3.0D+00)
      PARAMETER (MAXML=1024)
C
      CHARACTER*10  ATMNAM
C
      DOUBLE PRECISION METHOD
C
      DIMENSION ATMNAM(NFFAT),CORD(3,NFFAT),CORDSV(3,NFFAT),ZANF(NFFAT),
     *          FFGRD2(3,NFFAT),FFGRD1(3,NFFAT),QMGRD1(3,NAT),
     *          QMGRD2(3,NAT),LISTQM(*),QMGRD0(3,NAT),FFGRD0(3,NFFAT),
     *          NONLSTQ(*)
      DIMENSION TIMSTR(3)
C
      COMMON /ENRGYS/ ENUCR,EELCT,ETOT,SZ,SZZ,ECORE,ESCF,EERD,E1,E2,
     *                VEN,VEE,EPOT,EKIN,ESTATE(MXRT),STATN,EDFT(3),EDISP
      COMMON /FFDIMR/ IDIMER,IBREAK(81),N1213JMM,NMOLE,MATOM(MAXML),
     *                MCHARG(MAXML),MMULT(MAXML),MELEC(MAXML),RDIMER
      COMMON /FFENGY/ EN12,EN123,EN123R4,EN123B4,EN234W1,ENCHAR,ENLJR,
     *                ENLJD,ENPOL,XENPOL,ENRXN,XENRXN,ENRXNPOL,ENRXNR,
     *                EN12312,ENQUANP(30),
     *                ENBIAS,ENCENT,ENUCCH,ENCMAP,ENPOT,ENKIN,ENTOT
      COMMON /FFFREE/ SOL1CH,SOL1LJ,SOL1IM,SOLFRE1MM,
     *                SOL2CH,SOL2LJ,SOL2IM,SOLFRE2MM,
     *                ASOL1CH,ASOL1LJ,ASOL1IM,ASOL1MM,
     *                ASOL2CH,ASOL2LJ,ASOL2IM,ASOL2MM,
     *                PMF1BD,PMF1AG,PMF1DR,PMF1DB,PMF1WG,PMF1CM,
     *                PMF1CH,PMF1PO,PMF1LJ,PMF1IM,PMFFRE1MM,
     *                EN12A,EN12B,EN123A,EN123B,EN123R4A,EN123R4B,
     *                EN234W1A,EN234W1B,EN123B4A,EN123B4B,
     *                EN12312A,EN12312B,
     *                ENCHARA,ENCHARB,ENPOLA,ENPOLB,
     *                ENRXNA,ENRXNB,ENRXNPOLA,ENRXNPOLB,
     *                ENRXNRA,ENRXNRB,ENLJRA,ENLJRB,ENLJDA,ENLJDB,
     *                ENCMAPA,ENCMAPB,ENFIXSOA,ENFIXSOB,
     *                ENCENTA,ENCENTB,WSIMUL,WPERT1,WPERT2,
     *                IFEPTYP
      COMMON /FFPARM/ NFFAT,NBOND,NANGL,NDIHR,NDIHB,NCMAP,NWAGG,
     *                N1213J,N14J,NLKQMM,IDOCHG,IDOPOL,IDOLJ,IDOCMAP
      COMMON /FFMDPA/ DT,DT2,TEMP0,PRES0,POLTOL,VIR(3),PMEAN,VOLAV,
     *                ENPAV,ENKAV,TEMPAV,BERENDT,BERENDP,VELMAX,
     *                PMEANX,PMEANY,PMEANZ,NSTEP,KMASTER,KOUTACT(2),
     *                IHESS,INTALG,ITSTAT,IPSTAT,JOUT,KOUT,LOUT
      COMMON /FFMPT3/ NACTMM,LACTMM(2020),LFFDIPOLD,JUMBUP,
     *                NACTQM,LACTQM(2020),LFFOLDC,LFFQMVELSV,MMHESS,
     *                LFFQMCHG,LFFQMCHGB,ISWAP,R2SWAP,DFTBMM
      COMMON /FFPBSW/ XBOX,YBOX,ZBOX,SWF,SWFDX,SWFDY,SWFDZ,
     *                SWRA,ONESWRA,SWRA2,ONESWRA2,
     *                SWRB,ONESWRB,SWRB2,ONESWRB2,ONESWRB4,
     *                SWFDUM3,SWFDUM4,SWFDUM5,
     *                SWRAQ,ONESWRAQ,SWRAQ2,ONESWRAQ2,
     *                SWRBQ,ONESWRBQ,SWRBQ2,ONESWRBQ2,ONESWRBQ4,
     *                SWFDUM3Q,SWFDUM4Q,SWFDUM5Q,
     *                QMSIZE,QMCX,QMCY,QMCZ,QMCXSV,QMCYSV,QMCZSV,
     *                CENTX,CENTY,CENTZ,BUFWID1,BUFWID2,RDAMP,
     *                EFIELDX,EFIELDY,EFIELDZ,QMCXSV2,QMCYSV2,QMCZSV2,
     *                EPS1RB,EPS1RB3,ONEXBOX,ONEYBOX,ONEZBOX,
     *                LQMCT,MXLIST1,NTODO,NTODOSV,NTODOQ,
     *                ISWITCH,ISHIFT,IPOLSHF,
     *                LFFLSTCELL,LFFCORDSV,
     *                LFFPOLSV,LFFCORDSV2,LFFNONLS2,LFFCORDSVQ,
     *                LFFMVFASTS2,LFFMVFASTS3,LFFMVFASTS4,
     *                LFFMVFASTL2,LFFMVFASTL3,LFFMVFASTL4,
     *                MXCHECK,MXLIST2,NTODO2,NTODO2SV
      COMMON /FFQMPA/ ENFFQM2,ENPAV2,
     *                SCFTYP2,TDDFT2,MPLEVL2,CITYP2,
     *                ICHARG2,MULT2,IDOQM2,IREDOX,IQMPKA,IQMRXN,
     *                MATOMA,MCHARGA,MULTA,MELEA,
     *                MATOMB,MCHARGB,MULTB,MELEB,
     *                IECPX,NSHELLX,IMP,JMP,ICORSH,IGTF,
     *                LFFZANX,LFFCLPX,LFFZLPX,LFFNLPX,LFFKFRSTX,
     *                LFFKLASTX,LFFLMAXX,LFFLPSKIPX,LFFIZCOREX,
     *                LFFCX,LFFIANX,LFFEXX,LFFCSX,LFFCPX,LFFCDX,
     *                LFFCFX,LFFCGX,LFFCHX,LFFCIX,LFFKSTARTX,
     *                LFFKATOMX,LFFKTYPEX,LFFKNGX,LFFKLOCX,
     *                LFFMINX,LFFMAXX,LFFMPTYPX,LFFAN0X,
     *                LFFALPN0X,LFFAN1X,LFFALPN1X,LFFMPSKPX,
     *                LFFNOAN0X,LFFNOAN1X,LFFBPARX,LFFEXPMPX,
     *                LFFCSMPX,LFFCPMPX,LFFCDMPX,LFFCFMPX,LFFMPSKIPX,
     *                LFFNOCOSHX,LFFMPKSTAX,LFFMPKNGX,LFFMPKTYPX,
     *                LFFMPKMINX,LFFMPKMAXX,LFFMPKLOCX,LFFANAMX
      COMMON /FFRMSD/ DIPT(8),TIMGYRA,TIMRALL,
     *                NATPDB,NGYRA,NDIEL,NRALL,
     *                LFFRALL0,NRMSD,LFFRMSD0,KFREEAB(201),
     *                NRIJMM,IJRMM(2,100),NRIJQM,IJRQM(2,100),
     *                NAIJKMM,IJKMM(3,100),NAIJKQM,IJKQM(3,100),
     *                NFIXMM,IFIXMM(200),NFIXQM,IFIXQM(200)
      COMMON /FFRXN / RXNEPS,RSPHSOL,ISPHSOL
      COMMON /FFSPH / SPHRAD,SPHSIG,SPHEPS,IADDWAT
      COMMON /FFUMBR/ UMBFC,UMBR0,UMBSIZE,
     *                NUMBBIN,NUMBATM(6),NUMBTYP,LFFUMBHIS,
     *                UM2FC,UM2R0,UM2SIZE,
     *                NUM2BIN,NUM2ATM(6),NUM2TYP,LFFUM2HIS
      COMMON /FMCOM / X(1)
      COMMON /FUNCT / E,EG(3,MXATM)
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /IOFILE/ IR,IW,IP,IJK,IJKT,IDAF,NAV,IODA(950)
      COMMON /OPTGRD/ XX(3*MXATM),ES,FE(20),
     *                CONVF,FMAXT,DXMAXT,RMAX,RMIN,RLIM,
     *                EIGMAX,EIGMIN,GRDERR,FRMS,FMAX,TRMAX,TRMIN,
     *                IC(20),MSTEP,NSERCH,NPMAX,NP,IFOLOW,
     *                NNEG,IUPHSS,IEXIT,ITRUPD,IPAD,KDIAGH,NPRICO
      COMMON /OUTPUT/ NPRINT,ITOL,ICUT,NORMF,NORMP,NOPK
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /RUNLAB/ TITLE(10),ANAM(MXATM),BNAM(MXATM),BFLAB(MXAO)
      COMMON /RUNOPT/ RUNTYP,EXETYP,NEVALS,NGLEVL,NHLEVL
      COMMON /ZMAT  / NZMAT,NZVAR,NVAR,NSYMC,LINEAR
C
C     HUI LI, APR 2011, LINCOLN
C
      IF(NSTEP.LT.0) RETURN
C
      CALL SEQOPN(36,'TRAJECT','NEW',.FALSE.,'FORMATTED')
      IF(IADDWAT.EQ.2 .AND. ISPHSOL.GT.0) RSPHSOL=RSPHSOL*1.0D+40
      NPRINT = -5
C
      CALL VCLR(QMGRD2,1,3*NAT  )
      CALL VCLR(FFGRD2,1,3*NFFAT)
      ICONV   = 0
C
C     -- SET SOME CONTROL VARIABLES --
      NCOORD = 3*NAT
      NCVAL  = NCOORD
      MSTEP  = NSTEP
      NZMAT  = 0
      MODE   = 0
      OUT    = .FALSE.
      MAX2   = NCVAL
C
      CALL SIGINI(MODE,RUNTYP,NCVAL,NCOORD,METHOD,OUT,
     *            GOTEG,NPRT,NPUN,ITBMAT,STPT,STSTEP,
     *            PROJCT,DUMMY,DUMMY,MAX2,NPRTHS)
C
C     -- READY TO TAKE OFF ?
C
      DO ISTEP = 0, NSTEP
         IF(ISTEP.EQ.100.AND.IADDWAT.EQ.2 .AND. ISPHSOL.GT.0)
     *   RSPHSOL = RSPHSOL*1.0D-40
C
C        - CALCULATE ENERGY AND GRADIENT AT NEW COORDINATES
         CALL NONBOND(ISTEP,CORD,CORDSV,CORDSV2,CORDSVQ,
     *                NONLS1,NONLS2,
     *                NONLSTQ,LSTCELL,
     *                MVFASTS2,MVFASTS3,MVFASTS4,
     *                MVFASTL2,MVFASTL3,MVFASTL4,
     *                NONLSA,NONLSB,
     *                NONLSPMA,NONLSPMB)
         CALL VCLR(QMGRD2,1,3*NAT  )
         CALL VCLR(FFGRD2,1,3*NFFAT)
         CALL VCLR(VIR,1,3)    !  NEED TO DO THIS - HUI LI
         IF(IDIMER.EQ.0)THEN
            CALL GRADX
         ELSE
            CALL DIMERX
         END IF
C        -- ESCF WAS ZERO IF SCF NOT CONVERGED --
         IF(ABS(ESCF).LT.1.0D-12) THEN
            IF(MASWRK) WRITE(IW,'(/A/)')
     *         ' ERROR: SCF WAS NOT CONVERGED. QUANPOL MUST STOP.'
            CALL ABRT
         END IF
         CALL DCOPY(3*NAT  ,EG,1,QMGRD2,1)
         IF(IFEPTYP.GT.0.AND.MATOMB.GT.0) THEN
            CALL SAVEABPROP(1)
            CALL DCOPY(3*MATOMA,C,1,X(LFFCX),1)
            CALL SAVEFFDATA
            CALL SETQMAB(2)
            CALL SETFFDATAB(2)
            CALL VCLR(QMGRD2,1,3*NAT  )
            CALL VCLR(FFGRD2,1,3*NFFAT)
            CALL VCLR(VIR   ,1,3)
            CALL GRADX
C           -- ESCF WAS ZERO IF SCF NOT CONVERGED --
            IF(ABS(ESCF).LT.1.0D-12) THEN
               IF(MASWRK) WRITE(IW,'(/A/)')
     *            ' ERROR: SCF WAS NOT CONVERGED. QUANPOL MUST STOP.'
               CALL ABRT
            END IF
            CALL DCOPY(3*NAT,EG,1,QMGRD2,1)
            CALL SAVEABPROP(2)
            CALL SETQMAB(1)
            CALL SETFFDATAB(1)
            CALL MIXABPROP
         END IF
         E = ETOT
C
C        - COMBINE QM AND MM GRADIENTS -
C
         IF(LISTQM(NFFAT+NAT+1).GT.0)THEN
            DO IAT = 1, NAT
               KFFAT = LISTQM(NFFAT+IAT)
               IF(KFFAT.GT.0) THEN
                  QMGRD2(1,IAT)   = QMGRD2(1,IAT) + FFGRD2(1,KFFAT)
                  QMGRD2(2,IAT)   = QMGRD2(2,IAT) + FFGRD2(2,KFFAT)
                  QMGRD2(3,IAT)   = QMGRD2(3,IAT) + FFGRD2(3,KFFAT)
                  FFGRD2(1,KFFAT) = QMGRD2(1,IAT)
                  FFGRD2(2,KFFAT) = QMGRD2(2,IAT)
                  FFGRD2(3,KFFAT) = QMGRD2(3,IAT)
               END IF
            ENDDO
         END IF
C
C        - ZERO OFF SOME FORCES -
C
         IF(NACTQM.GT.0) THEN
            DO KOPT = 1, NACTQM
               IAT = LACTQM(KOPT)
               IF(IAT.GT.0)QMGRD2(1,IAT)=QMGRD2(1,IAT)+1.0D+03
            ENDDO
            DO IAT = 1, NAT
               IF(QMGRD2(1,IAT).GT.0.5D+03) THEN
                  QMGRD2(1,IAT) = QMGRD2(1,IAT)-1.0D+03
               ELSE
                  DO III = 1, 3
                     QMGRD2(III,IAT) = ZERO
                     EG(III,IAT)     = ZERO
                  ENDDO
               END IF
            ENDDO
         END IF
         DO KFIX=1,NFIXQM
            DO III = 1, 3
               QMGRD2(III,IFIXQM(KFIX)) = ZERO
            ENDDO
         ENDDO
         IF(NACTMM.GT.0) THEN
            DO KOPT = 1, NACTMM
               IFFAT = LACTMM(KOPT)
               IF(IFFAT.GT.0)FFGRD2(1,IFFAT)=FFGRD2(1,IFFAT)+1.0D+03
            ENDDO
            DO IFFAT = 1, NFFAT
               IF(FFGRD2(1,IFFAT).GT.0.5D+03) THEN
                  FFGRD2(1,IFFAT) = FFGRD2(1,IFFAT)-1.0D+03
               ELSE
                  DO III = 1, 3
                     FFGRD2(III,IFFAT) = ZERO
                  ENDDO
               END IF
            ENDDO
         END IF
         DO KFIX=1,NFIXMM
            DO III = 1, 3
               FFGRD2(III,IFIXMM(KFIX)) = ZERO
            ENDDO
         ENDDO
C
C        - TEST CONVERGENCE
C
         GRDMAX = ZERO
         GRDRMS = ZERO
         DO IAT=1, NAT
            GRDMAX=MAX(GRDMAX,ABS(QMGRD2(1,IAT)))
            GRDMAX=MAX(GRDMAX,ABS(QMGRD2(2,IAT)))
            GRDMAX=MAX(GRDMAX,ABS(QMGRD2(3,IAT)))
            GRDRMS=GRDRMS + QMGRD2(1,IAT)**2
     *                    + QMGRD2(2,IAT)**2
     *                    + QMGRD2(3,IAT)**2
         ENDDO
         DO IFFAT=1, NFFAT
            IF(LISTQM(IFFAT).EQ.0) THEN
            GRDMAX=MAX(GRDMAX,ABS(FFGRD2(1,IFFAT)))
            GRDMAX=MAX(GRDMAX,ABS(FFGRD2(2,IFFAT)))
            GRDMAX=MAX(GRDMAX,ABS(FFGRD2(3,IFFAT)))
            GRDRMS=GRDRMS + FFGRD2(1,IFFAT)**2
     *                    + FFGRD2(2,IFFAT)**2
     *                    + FFGRD2(3,IFFAT)**2
            END IF
         ENDDO
         NQMMM = LISTQM(NFFAT+NAT+1)
         MQMMM = 0
         DO KFIX=1,NFIXQM
            IAT = IFIXQM(KFIX)
            IF(LISTQM(NFFAT+IAT).GT.0) MQMMM = MQMMM + 1
         ENDDO
         IF((NACTMM+NACTQM).EQ.0)
     *   GRDRMS = SQRT(GRDRMS/(3*(NAT+NFFAT-NQMMM
     *                           -NFIXQM-NFIXMM+MQMMM)))
         IF((NACTMM+NACTQM).GT.0)
     *   GRDRMS = SQRT(GRDRMS/(3*(NACTMM+NACTQM-NQMMM
     *                           -NFIXMM-NFIXQM+MQMMM)))
         IF(GRDMAX.LT.CONVF.AND.GRDRMS.LT.(CONVF/THREE)) ICONV = 1
C
C        - RESET JUMBPOT=12 R0 -
         IF(NUMBTYP.EQ.12.AND.JUMBUP.NE.0) THEN
            IF(ISTEP.GT.0.AND.MOD(ISTEP,100).EQ.0) THEN
               N1    = NUMBATM(1)
               N2    = NUMBATM(2)
               CX    = CORD(1,N1) - CORD(1,N2)
               CY    = CORD(2,N1) - CORD(2,N2)
               CZ    = CORD(3,N1) - CORD(3,N2)
               R2    = CX*CX + CY*CY + CZ*CZ
               R     = SQRT(R2)
               IF(JUMBUP.GT.0)UMBR0=UMBR0+(UMBR0 - R)
               IF(JUMBUP.LT.0)UMBR0=R    +(R - UMBR0)
               IF(MASWRK) THEN
                  WRITE(IW,'(/1X,A,F10.6/)')
     *            'JUMBPOT=12 R0 IS RESET TO BE ',UMBR0*TOANGS
               END IF
               IF(UMBR0.LT.0.7D+00*TOBOHR.OR.
     *            UMBR0.GT.3.0D+00*TOBOHR) THEN
                  IF(MASWRK) THEN
                  WRITE(IW,'(/1X,A,A/)')
     *            'ERROR: JUMBPOT=12 R0 IS TOO SHORT OR TOO LONG. ',
     *            'THERE IS UNLIKELY A SADDLE POINT.'
                  END IF
                  RETURN
               END IF
            END IF
            IF(ENBIAS.GT.1.6D-06) ICONV = 0
         END IF
C
C        -- CALCULATE PROPERTIES
         IIISTEP=ISTEP
         CALL OPTPROP(IIISTEP,ICONV)
         IF(MASWRK) WRITE(IW,'(A,I6,A,F24.10,A,F10.7,A,F10.7)')
     *              ' NSERCH=',ISTEP,
     *              '  E=',ENTOT,
     *              '  GRAD. MAX=',GRDMAX,
     *              '  R.M.S.=',GRDRMS
         IF(MASWRK) WRITE(IW,'(A,I6,A,F12.9,A,F12.9,A,F24.10,A)')
     *              ' OPT STEP=',ISTEP,
     *              '  MAX GRAD=',GRDMAX,
     *              '  RMS GRAD=',GRDRMS,
     *              '  E=',ENTOT*TOKCAL,' KCAL/MOL'
C
         CALL TMDATE(TIMSTR)
         IF(MOD(ISTEP,KOUT).EQ.0) CALL TIMIT(1)
         IF(MASWRK .AND. (MOD(ISTEP,KOUT).EQ.0.OR.ISTEP.EQ.NSTEP.OR.
     *      ICONV.EQ.1))THEN
C
            WRITE(36,'(1X,A,I10)')
     *      'RESTART COORDINATES FOR QUANPOL AT OPT STEP',ISTEP
            WRITE(36,'(A,3(A,F15.10,1X),A)')
     *      ' $QUANPO ',
     *      'CENTX=',CENTX*TOANGS,
     *      'CENTY=',CENTY*TOANGS,
     *      'CENTZ=',CENTZ*TOANGS,
     *      '$END'
C
            WRITE(36,*)'$DATA'
            WRITE(36,'(A,3(F8.3,1X),A,I3,A,I10)')
     *      ' QM CENTER = ',QMCX*TOANGS,QMCY*TOANGS,QMCZ*TOANGS,
     *      'NEAR ATOM',LQMCT,'     AT OPT STEP',ISTEP
            WRITE(36,*)'C1'
            DO IAT = 1, NAT
               WRITE(36,999)ANAM(IAT),ZAN(IAT),
     *         C(1,IAT)*TOANGS,
     *         C(2,IAT)*TOANGS,
     *         C(3,IAT)*TOANGS
            ENDDO
            IF(IFEPTYP.GT.0)THEN
               DO IATB = 1, MATOMB
                  WRITE(36,999)ANAM(MATOMA+IATB),  ! ANAM IS THE ORIGINAL
     *            X(LFFZANX +MATOMA+IATB-1),
     *            X(LFFCX+3*(MATOMA+IATB-1)  )*TOANGS,
     *            X(LFFCX+3*(MATOMA+IATB-1)+1)*TOANGS,
     *            X(LFFCX+3*(MATOMA+IATB-1)+2)*TOANGS
               ENDDO
            END IF
            WRITE(36,*)'$END'
            WRITE(36,'(1X,A,I10,10X,A,F19.2,A)')
     *      '$FFDATA          !   OPT STEP',ISTEP
            IF(IFEPTYP.EQ.1) THEN
               WRITE(36,'(1X,A,I10,10X,A,F19.2,A)')
     *         '$FFDATB          !   OPT STEP',ISTEP
            END IF
            WRITE(36,*)'COORDINATES  NUC                   X',
     *          '                   Y                   Z'
            DO IFFAT = 1, NFFAT
               WRITE(36,1000)ATMNAM(IFFAT),ZANF(IFFAT),
     *         CORD(1,IFFAT)*TOANGS,
     *         CORD(2,IFFAT)*TOANGS,
     *         CORD(3,IFFAT)*TOANGS
            ENDDO
            WRITE(36,*)'STOP'
C
            IF(NFFAT.LE.100) THEN
            WRITE(IW,*)'$DATA'
            WRITE(IW,'(A,3(F8.3,1X),A,I3,A,I10)')
     *      ' QM CENTER = ',QMCX*TOANGS,QMCY*TOANGS,QMCZ*TOANGS,
     *      'NEAR ATOM',LQMCT,'     AT OPT STEP',ISTEP
            WRITE(IW,*)'C1'
            DO IAT = 1, NAT
               WRITE(IW,999)ANAM(IAT),ZAN(IAT),
     *         C(1,IAT)*TOANGS,
     *         C(2,IAT)*TOANGS,
     *         C(3,IAT)*TOANGS
            ENDDO
            IF(IFEPTYP.GT.0)THEN
               DO IATB = 1, MATOMB
                  WRITE(IW,999)ANAM(MATOMA+IATB),  ! ANAM IS THE ORIGINAL
     *            X(LFFZANX +MATOMA+IATB-1),
     *            X(LFFCX+3*(MATOMA+IATB-1)  )*TOANGS,
     *            X(LFFCX+3*(MATOMA+IATB-1)+1)*TOANGS,
     *            X(LFFCX+3*(MATOMA+IATB-1)+2)*TOANGS
               ENDDO
            END IF
            WRITE(IW,*)'$END'
            WRITE(IW,'(1X,A,I10,10X,A,F19.2,A)')
     *      '$FFDATA          !   OPT STEP',ISTEP
            IF(IFEPTYP.EQ.1) THEN
               WRITE(IW,'(1X,A,I10,10X,A,F19.2,A)')
     *         '$FFDATB          !   OPT STEP',ISTEP
            END IF
            WRITE(IW,*)'COORDINATES  NUC                   X',
     *          '                   Y                   Z'
            DO IFFAT = 1, NFFAT
               WRITE(IW,1000)ATMNAM(IFFAT),ZANF(IFFAT),
     *         CORD(1,IFFAT)*TOANGS,
     *         CORD(2,IFFAT)*TOANGS,
     *         CORD(3,IFFAT)*TOANGS
            ENDDO
            WRITE(IW,*)'STOP'
            END IF
            IF(NFFAT.GT.100.AND.NACTMM.GT.0) THEN
            WRITE(IW,*)'$DATA         ! THESE ARE LACTMM ATOMS'
            WRITE(IW,'(A,3(F8.3,1X),A,I3,A,I10)')
     *      ' QM CENTER = ',QMCX*TOANGS,QMCY*TOANGS,QMCZ*TOANGS,
     *      'NEAR ATOM',LQMCT,'     AT OPT STEP',ISTEP
            WRITE(IW,*)'C1'
            DO IAT = 1, NAT
               ZCHG = ZAN(IAT)
               IF(LISTQM(NFFAT+IAT).GT.0)THEN
                  ZCHG = ZANF(LISTQM(NFFAT+IAT))
               END IF
               CX    = C(1,IAT) - QMCX
               CY    = C(2,IAT) - QMCY
               CZ    = C(3,IAT) - QMCZ
               PBCX  = XBOX*ANINT(CX*ONEXBOX)
               PBCY  = YBOX*ANINT(CY*ONEYBOX)
               PBCZ  = ZBOX*ANINT(CZ*ONEZBOX)
               CX    = CX - PBCX + QMCX
               CY    = CY - PBCY + QMCY
               CZ    = CZ - PBCZ + QMCZ
               WRITE(IW,999)ANAM(IAT),ZCHG,
     *         CX*TOANGS,
     *         CY*TOANGS,
     *         CZ*TOANGS
            ENDDO
            DO I=1,NACTMM
               IFFAT  = LACTMM(I)
               IF(LISTQM(IFFAT).EQ.0.OR.LISTQM(IFFAT).GT.NAT) THEN
                  CX    = CORD(1,IFFAT) - QMCX
                  CY    = CORD(2,IFFAT) - QMCY
                  CZ    = CORD(3,IFFAT) - QMCZ
                  PBCX  = XBOX*ANINT(CX*ONEXBOX)
                  PBCY  = YBOX*ANINT(CY*ONEYBOX)
                  PBCZ  = ZBOX*ANINT(CZ*ONEZBOX)
                  CX    = CORD(1,IFFAT) - PBCX
                  CY    = CORD(2,IFFAT) - PBCY
                  CZ    = CORD(3,IFFAT) - PBCZ
                  WRITE(IW,999)ATMNAM(IFFAT)(1:8),ZANF(IFFAT),
     *            CX*TOANGS,
     *            CY*TOANGS,
     *            CZ*TOANGS
               END IF
            ENDDO
            WRITE(IW,*)'$END'
            END IF
C
C           -- PRINT KOUTACT ATOMS --
C              HERE ONLY $FFDATA IS CONSIDERED, $FFDATB IS NOT.
C
            IF(KOUTACT(1).GT.0) THEN
            WRITE(IW,'(/A,I10)')'!KOUTACT ATOMS AROUND $FFDATA ATOM ',
     *      KOUTACT(1)
            WRITE(IW,'(1X,A,I10)')
     *      '$FFDATA ! KOUTACT     OPT STEP',ISTEP
            ACTX = CORD(1,KOUTACT(1))
            ACTY = CORD(2,KOUTACT(1))
            ACTZ = CORD(3,KOUTACT(1))
            ACTR = DBLE(KOUTACT(2))*1.0D-08*TOBOHR
            ACTR2= ACTR**2
            DO IFFAT=1,NFFAT
               CX    = CORD(1,IFFAT) - ACTX
               CY    = CORD(2,IFFAT) - ACTY
               CZ    = CORD(3,IFFAT) - ACTZ
               PBCX  = XBOX*ANINT(CX*ONEXBOX)
               PBCY  = YBOX*ANINT(CY*ONEYBOX)
               PBCZ  = ZBOX*ANINT(CZ*ONEZBOX)
               CX    = CX - PBCX
               CY    = CY - PBCY
               CZ    = CZ - PBCZ
               RCXYZ = CX**2 + CY**2 + CZ**2
               IF(RCXYZ.LE.ACTR2) THEN
                  WRITE(IW,1000)ATMNAM(IFFAT),ZANF(IFFAT),
     *            (CX+ACTX)*TOANGS,(CY+ACTY)*TOANGS,(CZ+ACTZ)*TOANGS
               END IF
            ENDDO
            WRITE(IW,*)'$END ! KOUTACT'
            END IF
C
            CALL TIMIT(1)
            CALL FLSHBF(IW)
            CALL FLSHBF(36)
         END IF
C
         IF(ICONV.EQ.1) GOTO 500
C
C        - CALCULATE NEW COORDINATES
         IF(NSTEP.GT.0) THEN
         DO 201 IAT = 1, NAT
            IF((QMGRD2(1,IAT)+QMGRD2(2,IAT)+QMGRD2(3,IAT))
     *         .EQ.ZERO) GOTO 201
            DO II = 1, 3
               DISPLC =  -(0.58D+00*QMGRD2(II,IAT)+
     *                     0.28D+00*QMGRD1(II,IAT)+
     *                     0.14D+00*QMGRD0(II,IAT) )*
     *                     ABS(LOG(ABS(CONVF)+ABS(QMGRD2(II,IAT))))
     *                     *0.23D+00
               IF(DISPLC.GT. 0.20D+00) DISPLC = 0.20D+00
               IF(DISPLC.LT.-0.20D+00) DISPLC =-0.20D+00
               C(II,IAT) = C(II,IAT) + DISPLC
            ENDDO
  201    CONTINUE
         END IF
         CALL DCOPY(3*NAT,QMGRD1,1,QMGRD0,1)
         CALL DCOPY(3*NAT,QMGRD2,1,QMGRD1,1)
         DO 200 IFFAT = 1, NFFAT
            IF((FFGRD2(1,IFFAT)+FFGRD2(2,IFFAT)+FFGRD2(3,IFFAT))
     *         .EQ.ZERO) GOTO 200
            DO II = 1, 3
               DISPLC =  -(0.58D+00*FFGRD2(II,IFFAT)+
     *                     0.28D+00*FFGRD1(II,IFFAT)+
     *                     0.14D+00*FFGRD0(II,IFFAT) )*
     *                     ABS(LOG(ABS(CONVF)+ABS(FFGRD2(II,IFFAT))))
     *                     *0.23D+00
               IF(DISPLC.GT. 0.20D+00) DISPLC = 0.20D+00
               IF(DISPLC.LT.-0.20D+00) DISPLC =-0.20D+00
               CORD(II,IFFAT) = CORD(II,IFFAT) + DISPLC
            ENDDO
  200    CONTINUE
         CALL DCOPY(3*NFFAT,FFGRD1,1,FFGRD0,1)
         CALL DCOPY(3*NFFAT,FFGRD2,1,FFGRD1,1)
C        -- SYNCHRONIZE CORD EVERY STEP --
         IF(LISTQM(NFFAT+NAT+1).GT.0)THEN
            DO IAT = 1, NAT
               KFFAT = LISTQM(NFFAT+IAT)
               IF(KFFAT.GT.0) THEN
                  C(1,IAT) = CORD(1,KFFAT)
                  C(2,IAT) = CORD(2,KFFAT)
                  C(3,IAT) = CORD(3,KFFAT)
               END IF
            ENDDO
         END IF
         IF(GOPARR) THEN
            CALL DDI_BCAST(461,'F',C,3*NAT,MASTER)
            CALL DDI_BCAST(462,'F',CORD,3*NFFAT,MASTER)
         END IF
C
C        -- ADJUST QMCX,Y,Z ON THE FLY --
C           ALSO UPDATE LQMCT
C           (DO THIS ONLY WHEN PBC AND SWRB2 ARE USED)
         IF(SWRB2.LT.1.0D+08) THEN
            XMAX = -1.0D+30
            YMAX = -1.0D+30
            ZMAX = -1.0D+30
            XMIN =  1.0D+30
            YMIN =  1.0D+30
            ZMIN =  1.0D+30
            DO IAT = 1,NAT
               XMAX = MAX(XMAX,C(1,IAT))
               YMAX = MAX(YMAX,C(2,IAT))
               ZMAX = MAX(ZMAX,C(3,IAT))
               XMIN = MIN(XMIN,C(1,IAT))
               YMIN = MIN(YMIN,C(2,IAT))
               ZMIN = MIN(ZMIN,C(3,IAT))
            ENDDO
            QMSIZE = ZERO
            QMSIZE = MAX(QMSIZE,XMAX-XMIN)
            QMSIZE = MAX(QMSIZE,YMAX-YMIN)
            QMSIZE = MAX(QMSIZE,ZMAX-ZMIN)
            QMSIZE = QMSIZE*1.732D+00
            QMCX   = (XMAX+XMIN)*PT5
            QMCY   = (YMAX+YMIN)*PT5
            QMCZ   = (ZMAX+ZMIN)*PT5
            R2NEAR = 100.0D+00
            NEAR   = 0
            DO IAT = 1,NAT
               XI = C(1,IAT) - QMCX
               YI = C(2,IAT) - QMCY
               ZI = C(3,IAT) - QMCZ
               R2 = XI**2 + YI**2 + ZI**2
               IF(R2.LT.R2NEAR) THEN
                  R2NEAR = R2
                  NEAR   = IAT
               END IF
            ENDDO
            LQMCT = NEAR
            IF(MASWRK) WRITE(IW,'(A,3(F16.10,1X),A,I4,A,I7)')
     *      ' QM CENTER = ',QMCX*TOANGS,QMCY*TOANGS,QMCZ*TOANGS,
     *      'NEAR ATOM',LQMCT,' AT STEP ',ISTEP+1
         END IF
      ENDDO
C
C     -- NOT LOCATED --
C
      IF(MASWRK) THEN
         WRITE(IW,*)' '
         WRITE(IW,*)'QUANPOL OPTIMIZATION FAILED. ',
     *              'PLEASE RESTART THE JOB USING THE LAST GEOMETRY.'
         WRITE(IW,*)' '
      END IF
      CALL FLSHBF(IW)
      CALL FLSHBF(36)
      CALL SEQCLO(36,'KEEP')
      RETURN
C
C     -- LOCATED --
C
  500 CONTINUE
      IF(MASWRK) THEN
         WRITE(IW,*)' '
         WRITE(IW,*)'================= QUANPOL OPTIMIZATION',
     *              ' SUCCESSFULLY COMPLETED ================='
         WRITE(IW,*)' '
      END IF
C
  999 FORMAT(1X,A8,3X,F5.1,1X,F19.13,1X,F19.13,1X,F19.13)
 1000 FORMAT(1X,A10,1X,F5.1,1X,F19.13,1X,F19.13,1X,F19.13)
C
      CALL FLSHBF(IW)
      CALL FLSHBF(36)
      CALL SEQCLO(36,'KEEP')
C
      RETURN
      END
C*MODULE QUANPOC  *DECK FFOPT3
!>
!> @brief    pure MM optimization
!>
!> @author   Nandun Thellamurege
!>           - May 2011
!>
!> @details  BFGS Hessian method
!>
      SUBROUTINE FFOPT3(ATMNAM,CORD,CORDSV,ZANF,
     *                  ZMAS,CHARG,POL,POLSV,DIP,
     *                  FIELD1,FIELD2,FIELD3,
     *                  SIG,EPS,SIG2,EPS2,
     *                  BOND0,FCBOND,FCSTBD,
     *                  ANGL0,FCANGL,FCWAGG,
     *                  DIHB0,FCDIHB,FCDIHR,
     *                  VROT,NNN,GAMA,IPAIR,
     *                  KLIST,LLIST,KBLST,MLIST,NLIST,
     *                  L1213J,L14J,
     *                  FFGRD2,FCLJTP,NTYPE,
     *                  XTS,YTS,ZTS,CMAT1,
     *                  POT1,POT2,QRXN1,QRXN2,NTS,
     *                  NONLS1,NONLSTQ,MAPLST,CMAPCO,
     *                  CRDNEW,CRDOLD,
     *                  LSTCELL,NONLS2,CORDSV2,CORDSVQ,
     *                  MVFASTS2,MVFASTS3,MVFASTS4,
     *                  MVFASTL2,MVFASTL3,MVFASTL4,
     *                  AFIX,QFIX,
     *                  RFIX,IDATOM,DAI,IDDAI,
     *                  VFIX1,VFIX2,XTSFIX,YTSFIX,ZTSFIX,
     *                  NONLSA,NONLSB,L1213A,L1213B,L14A,L14B,
     *                  CHARGB,
     *                  SIGB,EPSB,SIG2B,EPS2B,
     *                  NONLSPMA,
     *                  L1213PMA,L14PMA,CORDB,
     *                  LSBONDPMA,LSANGLPMA,
     *                  LSDIHRPMA,LSDIHBPMA,
     *                  LSWAGGPMA,LSCMAPPMA,
     *                  NONLSPMB,L1213PMB,L14PMB,
     *                  LSBONDPMB,LSANGLPMB,
     *                  LSDIHRPMB,LSDIHBPMB,
     *                  LSWAGGPMB,LSCMAPPMB,UMBHIS,UM2HIS)
      USE MX_LIMITS,ONLY:MXFRG,mxatm,mxao
C
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      LOGICAL GOPARR,DSKWRK,MASWRK
      LOGICAL CVGED,OUT,STPT,PROJCT,GOTEG
C
      PARAMETER (TOANGS=0.52917724924D+00)
      PARAMETER (TOBOHR=1.0D+00/TOANGS)
      PARAMETER (ZERO=0.0D+00)
C
      CHARACTER*10  ATMNAM
C
      DOUBLE PRECISION METHOD
      LOGICAL LINEAR
C
      DIMENSION ATMNAM(NFFAT),CORD(3,NFFAT),CORDSV(3,NFFAT),ZANF(NFFAT),
     *          ZMAS(NFFAT),
     *          CHARG(NFFAT),POL(NFFAT),DIP(3,NFFAT),
     *          FIELD1(3,NFFAT),FIELD2(3,NFFAT),FIELD3(3,NFFAT),
     *          SIG(NFFAT),EPS(NFFAT),SIG2(NFFAT),EPS2(NFFAT),
     *          BOND0(NBOND),FCBOND(NBOND),
     *          ANGL0(NANGL),FCANGL(NANGL),
     *          FCWAGG(NWAGG),DIHB0(NDIHB),FCDIHB(NDIHB),FCDIHR(3,*),
     *          VROT(NDIHR),NNN(NDIHR),GAMA(NDIHR),IPAIR(2,NBOND),
     *          KLIST(3,NANGL),
     *          LLIST(4,NDIHR),MLIST(4,NWAGG),
     *          L1213J(2,*),L14J(2,NDIHR),
     *          FFGRD2(3,NFFAT),NTYPE(*),
     *          XTS(NTS),YTS(NTS),ZTS(NTS),CMAT1(NTS,NTS),
     *          POT1(NTS),POT2(NTS),QRXN1(NTS),QRXN2(NTS),NONLS1(2,*),
     *          NONLSTQ(*),MAPLST(6,*),CMAPCO(4,4,24,24,3),
     *          CRDNEW(*)
      DIMENSION TIMSTR(3)
      DIMENSION ENALL(100)
C
      COMMON /DLCFRZ/ FVALUE(50),ITABLE(50),IFTYPE(50),NCONST
      COMMON /FFENGY/ EN12,EN123,EN123R4,EN123B4,EN234W1,ENCHAR,ENLJR,
     *                ENLJD,ENPOL,XENPOL,ENRXN,XENRXN,ENRXNPOL,ENRXNR,
     *                EN12312,ENQUANP(30),
     *                ENBIAS,ENCENT,ENUCCH,ENCMAP,ENPOT,ENKIN,ENTOT
      COMMON /FFFIXS/ ENFIXSO,FIXEPS,FIXTOL,FIXA,FIXQ,RALLMM,RALLQM,
     *                RADMM(200),RADQM(200),NRADMM,NRADQM,IFIXSOL,
     *                LFFDAI,LFFDAIT,LFFIDDAI,LFFIDTMP,LFFTMPTS,
     *                LFFAFIX,LFFIDATOM,LFFRFIX,LFFQFIX,NTSATM,
     *                LFFQFIXMP,LFFQFIXTA,LFFQFIXXY,
     *                LFFXTSFIX,LFFYTSFIX,LFFZTSFIX,
     *                LFFVFIX1,LFFVFIX2,NCYCLE,MXFFTS,NFFTS
      COMMON /FFFREE/ SOL1CH,SOL1LJ,SOL1IM,SOLFRE1MM,
     *                SOL2CH,SOL2LJ,SOL2IM,SOLFRE2MM,
     *                ASOL1CH,ASOL1LJ,ASOL1IM,ASOL1MM,
     *                ASOL2CH,ASOL2LJ,ASOL2IM,ASOL2MM,
     *                PMF1BD,PMF1AG,PMF1DR,PMF1DB,PMF1WG,PMF1CM,
     *                PMF1CH,PMF1PO,PMF1LJ,PMF1IM,PMFFRE1MM,
     *                EN12A,EN12B,EN123A,EN123B,EN123R4A,EN123R4B,
     *                EN234W1A,EN234W1B,EN123B4A,EN123B4B,
     *                EN12312A,EN12312B,
     *                ENCHARA,ENCHARB,ENPOLA,ENPOLB,
     *                ENRXNA,ENRXNB,ENRXNPOLA,ENRXNPOLB,
     *                ENRXNRA,ENRXNRB,ENLJRA,ENLJRB,ENLJDA,ENLJDB,
     *                ENCMAPA,ENCMAPB,ENFIXSOA,ENFIXSOB,
     *                ENCENTA,ENCENTB,WSIMUL,WPERT1,WPERT2,
     *                IFEPTYP
      COMMON /FFMDPA/ DT,DT2,TEMP0,PRES0,POLTOL,VIR(3),PMEAN,VOLAV,
     *                ENPAV,ENKAV,TEMPAV,BERENDT,BERENDP,VELMAX,
     *                PMEANX,PMEANY,PMEANZ,NSTEP,KMASTER,KOUTACT(2),
     *                IHESS,INTALG,ITSTAT,IPSTAT,JOUT,KOUT,LOUT
      COMMON /FFMPT3/ NACTMM,LACTMM(2020),LFFDIPOLD,JUMBUP,
     *                NACTQM,LACTQM(2020),LFFOLDC,LFFQMVELSV,MMHESS,
     *                LFFQMCHG,LFFQMCHGB,ISWAP,R2SWAP,DFTBMM
      COMMON /FFPARM/ NFFAT,NBOND,NANGL,NDIHR,NDIHB,NCMAP,NWAGG,
     *                N1213J,N14J,NLKQMM,IDOCHG,IDOPOL,IDOLJ,IDOCMAP
      COMMON /FFPBSW/ XBOX,YBOX,ZBOX,SWF,SWFDX,SWFDY,SWFDZ,
     *                SWRA,ONESWRA,SWRA2,ONESWRA2,
     *                SWRB,ONESWRB,SWRB2,ONESWRB2,ONESWRB4,
     *                SWFDUM3,SWFDUM4,SWFDUM5,
     *                SWRAQ,ONESWRAQ,SWRAQ2,ONESWRAQ2,
     *                SWRBQ,ONESWRBQ,SWRBQ2,ONESWRBQ2,ONESWRBQ4,
     *                SWFDUM3Q,SWFDUM4Q,SWFDUM5Q,
     *                QMSIZE,QMCX,QMCY,QMCZ,QMCXSV,QMCYSV,QMCZSV,
     *                CENTX,CENTY,CENTZ,BUFWID1,BUFWID2,RDAMP,
     *                EFIELDX,EFIELDY,EFIELDZ,QMCXSV2,QMCYSV2,QMCZSV2,
     *                EPS1RB,EPS1RB3,ONEXBOX,ONEYBOX,ONEZBOX,
     *                LQMCT,MXLIST1,NTODO,NTODOSV,NTODOQ,
     *                ISWITCH,ISHIFT,IPOLSHF,
     *                LFFLSTCELL,LFFCORDSV,
     *                LFFPOLSV,LFFCORDSV2,LFFNONLS2,LFFCORDSVQ,
     *                LFFMVFASTS2,LFFMVFASTS3,LFFMVFASTS4,
     *                LFFMVFASTL2,LFFMVFASTL3,LFFMVFASTL4,
     *                MXCHECK,MXLIST2,NTODO2,NTODO2SV
      COMMON /FFRMSD/ DIPT(8),TIMGYRA,TIMRALL,
     *                NATPDB,NGYRA,NDIEL,NRALL,
     *                LFFRALL0,NRMSD,LFFRMSD0,KFREEAB(201),
     *                NRIJMM,IJRMM(2,100),NRIJQM,IJRQM(2,100),
     *                NAIJKMM,IJKMM(3,100),NAIJKQM,IJKQM(3,100),
     *                NFIXMM,IFIXMM(200),NFIXQM,IFIXQM(200)
      COMMON /FFRXN / RXNEPS,RSPHSOL,ISPHSOL
      COMMON /FFSPH / SPHRAD,SPHSIG,SPHEPS,IADDWAT
      COMMON /FMCOM / X(1)
      COMMON /FUNCT / E,EG(3,MXATM)
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON /MACHSW/ KDIAG,ICORFL,IXDR,modio,mem10,lpnt10,mem10m
      COMMON /MASSES/ ZMASS(MXATM)
      COMMON /RUNOPT/ RUNTYP,EXETYP,NEVALS,NGLEVL,NHLEVL
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /OPTGRD/ XX(3,MXATM),ES,FE(20),
     *                CONVF,FMAXT,DXMAXT,RMAX,RMIN,RLIM,
     *                EIGMAX,EIGMIN,GRDERR,FRMS,FMAX,TRMAX,TRMIN,
     *                IC(20),MSTEP,NSERCH,NPMAX,NP,IFOLOW,
     *                NNEG,IUPHSS,IEXIT,ITRUPD,IPAD,KDIAGH,NPRICO
      COMMON /OPTEF / D(3*MXATM+6*MXFRG),OLDF(3*MXATM+6*MXFRG),
     *                VMODE(3*MXATM+6*MXFRG),RADIUS
      COMMON /ZMAT  / NZMAT,NZVAR,NVAR,NSYMC,LINEAR
      COMMON /ZMTALT/ NZMAT2,NZVAR2,NVAR2,NZMTRD,ICOORD
      COMMON /RUNLAB/ TITLE(10),ANAM(MXATM),BNAM(MXATM),BFLAB(MXAO)
      COMMON /FROZNC/ IFREEZ(3*MXATM),NFRZ
C
      DATA BLANK/8H        /
C
C     NANDUN THELLAMUREGE, MAY 2011, LINCOLN
C     HUI LI, MAY 6, 2014, USE LACTMM AND PRJGRD
C
      IF(NSTEP.LT.0) RETURN
C
      CALL SEQOPN(36,'TRAJECT','NEW',.FALSE.,'FORMATTED')
      CALL SEQOPN(39,'OPTHES1','NEW',.FALSE.,'FORMATTED')
      CALL SEQOPN(40,'OPTHES2','NEW',.FALSE.,'FORMATTED')
      IF(IADDWAT.EQ.2 .AND. ISPHSOL.GT.0) RSPHSOL=RSPHSOL*1.0D+40
C
C     -- DEFINE INFOA VARIABLES --
      IF(NACTMM.EQ.0) THEN
         NCOORD = 3*NFFAT
         NCVAL  = NCOORD
         CALL DCOPY(3*NFFAT,CORD,1,C,1)
         CALL DCOPY(NFFAT,ZANF,1,ZAN,1)
         DO IFFAT =1,NFFAT
            READ(UNIT=ATMNAM(IFFAT)(1:8),FMT='(A8)') ANAM(IFFAT)
            BNAM(IFFAT)=BLANK
            ZMASS(IFFAT)=ZMAS(IFFAT)/1822.88850204D+00
         ENDDO
      ELSE
         NCOORD = 3*NACTMM
         NCVAL  = NCOORD
         DO I=1,NACTMM
            IFFAT  = LACTMM(I)
            C(1,I) = CORD(1,IFFAT)
            C(2,I) = CORD(2,IFFAT)
            C(3,I) = CORD(3,IFFAT)
            ZAN(I) = ZANF(IFFAT)
            READ(UNIT=ATMNAM(IFFAT)(1:8),FMT='(A8)') ANAM(I)
            BNAM(I)=BLANK
            ZMASS(I)=ZMAS(IFFAT)/1822.88850204D+00
         ENDDO
      END IF
C
C     -- SET SOME CONTROL VARIABLES --
      MSTEP  = NSTEP
      NZMAT  = 0
      NFRZ   = 0
      NCONST = 0
      ICOORD = 8
      CVGED  = .FALSE.
      MODE   = 0
      OUT    = .FALSE.
      MAX2   = NCVAL
C
      NAT    = NFFAT
      IF(NACTMM.GT.0) NAT = NACTMM
      CALL SIGINI(MODE,RUNTYP,NCVAL,NCOORD,METHOD,OUT,
     *            GOTEG,NPRT,NPUN,ITBMAT,STPT,STSTEP,
     *            PROJCT,DUMMY,DUMMY,MAX2,NPRTHS)
      NAT    = 0
C
      DO ISTEP = 0, NSTEP
         NSERCH= ISTEP
         IF(ISTEP.EQ.100.AND.IADDWAT.EQ.2 .AND. ISPHSOL.GT.0)
     *   RSPHSOL = RSPHSOL*1.0D-40
C
         NPUN = -1
         IF(MOD(ISTEP,KOUT).EQ.0 .OR. ISTEP.EQ.NSTEP) NPUN = -1
C
C        -- CALCULATE ENERGY AND GRADIENT
         CALL NONBOND(ISTEP,CORD,CORDSV,CORDSV2,CORDSVQ,
     *                NONLS1,NONLS2,
     *                NONLSTQ,LSTCELL,
     *                MVFASTS2,MVFASTS3,MVFASTS4,
     *                MVFASTL2,MVFASTL3,MVFASTL4,
     *                NONLSA,NONLSB,
     *                NONLSPMA,NONLSPMB)
         CALL VCLR(FFGRD2,1,3*NFFAT)
         CALL VCLR(VIR,1,3)    !  NEED TO DO THIS - HUI LI
         CALL E00012(CORD,FFGRD2,BOND0,FCBOND,IPAIR,CORDB,
     *               LSBONDPMA,LSBONDPMB)
         CALL E00123(CORD,FFGRD2,ANGL0,FCANGL,KLIST,CORDB,
     *               LSANGLPMA,LSANGLPMB)
         CALL E12312(CORD,FFGRD2,ANGL0,KLIST,BOND0,FCSTBD,
     *               KBLST,CORDB,LSANGLPMA,LSANGLPMB)
         CALL E123B4(CORD,FFGRD2,DIHB0,FCDIHB,NLIST,CORDB,
     *               LSDIHBPMA,LSDIHBPMB)
         CALL E234W1(CORD,FFGRD2,      FCWAGG,MLIST,CORDB,
     *               LSWAGGPMA,LSWAGGPMB)
         CALL E123R4(CORD,FFGRD2,VROT,GAMA,NNN,LLIST,CORDB,
     *               LSDIHRPMA,LSDIHRPMB,FCDIHR)
         CALL ECMAP (CORD,FFGRD2,MAPLST,CMAPCO,CORDB,
     *               LSCMAPPMA,LSCMAPPMB)
         CALL ELJ126(CORD,FFGRD2,SIG,EPS,SIG2,EPS2,L14J,NONLS1,
     *               L1213J,SIGB,EPSB,SIG2B,EPS2B,
     *               NONLSA,NONLSB,L1213A,L1213B,L14A,L14B,
     *               NONLSPMA,L1213PMA,L14PMA,CORDB,
     *               NONLSPMB,L1213PMB,L14PMB,FCLJTP,NTYPE)
         CALL ESPHER(CORD,FFGRD2)
         CALL ECHARG(CORD,FFGRD2,CHARG,NONLS1,L1213J,L14J,
     *               CHARGB,NONLSA,NONLSB,L1213A,L1213B,L14A,L14B,
     *               NONLSPMA,L1213PMA,L14PMA,CORDB,
     *               NONLSPMB,L1213PMB,L14PMB)
         CALL UMBRELLA(CORD,FFGRD2,UMBHIS,UM2HIS)
         IF(IFIXSOL.EQ.0) THEN
            CALL CHGRXN(CORD,FFGRD2,CHARG,XTS,YTS,ZTS,
     *                  CMAT1,POT1,QRXN1,NTS)
            CALL POLRXN(CORD,FFGRD2,CHARG,POL,POLSV,DIP,
     *                  FIELD1,FIELD2,FIELD3,
     *                  XTS,YTS,ZTS,CMAT1,POT1,POT2,QRXN1,QRXN2,NTS,
     *                  NONLS1,L1213J)
         END IF
         IF(IFIXSOL.EQ.1) THEN
            NATSV = NAT   ! NEED THIS BECAUSE FIXSOL USES NFFAT+NAT
            NAT   = 0
            CALL FIXSOL(CORD,FFGRD2,CHARG,ZANF,AFIX,QFIX,
     *                  VFIX1,VFIX2,XTSFIX,YTSFIX,ZTSFIX,
     *                  RFIX,IDATOM,DAI,IDDAI,
     *                  POL,POLSV,DIP,FIELD1,FIELD2,FIELD3,
     *                  NONLS1,L1213J)
            NAT   = NATSV
         END IF
         IF(GOPARR) THEN
            CALL VCLR(ENALL,1,100)
            ENALL( 1) = EN12
            ENALL( 2) = EN123
            ENALL( 3) = EN123R4
            ENALL( 4) = EN123B4
            ENALL( 5) = EN234W1
            ENALL( 6) = ENCHAR
            ENALL( 7) = ENLJR
            ENALL( 8) = ENLJD
C           ENALL(12) = ENUCCH
            ENALL(13) = ENRXN
            ENALL(14) = ENRXNR
C           ENALL(15) = ENCENT
            ENALL(16) = ENCMAP
            ENALL(18) = SOL1CH
            ENALL(19) = SOL1LJ
            ENALL(21) = SOL2CH
            ENALL(22) = SOL2LJ
            ENALL(33) = PMF1BD
            ENALL(34) = PMF1AG
            ENALL(35) = PMF1DR
            ENALL(36) = PMF1DB
            ENALL(37) = PMF1WG
            ENALL(38) = PMF1CM
            ENALL(39) = PMF1CH
            ENALL(40) = PMF1LJ
            ENALL(41) = ENBIAS
            ENALL(42) = EN12312
            CALL DDI_GSUMF(2410,ENALL  ,42)
            CALL DDI_GSUMF(2411,FFGRD2,3*NFFAT)
            EN12      = ENALL( 1)
            EN123     = ENALL( 2)
            EN123R4   = ENALL( 3)
            EN123B4   = ENALL( 4)
            EN234W1   = ENALL( 5)
            ENCHAR    = ENALL( 6)
            ENLJR     = ENALL( 7)
            ENLJD     = ENALL( 8)
C           ENUCCH    = ENALL(12)
            ENRXN     = ENALL(13)
            ENRXNR    = ENALL(14)
C           ENCENT    = ENALL(15)
            ENCMAP    = ENALL(16)
            SOL1CH    = ENALL(18)
            SOL1LJ    = ENALL(19)
            SOL2CH    = ENALL(21)
            SOL2LJ    = ENALL(22)
            PMF1BD    = ENALL(33)
            PMF1AG    = ENALL(34)
            PMF1DR    = ENALL(35)
            PMF1DB    = ENALL(36)
            PMF1WG    = ENALL(37)
            PMF1CM    = ENALL(38)
            PMF1CH    = ENALL(39)
            PMF1LJ    = ENALL(40)
            ENBIAS    = ENALL(41)
            EN12312   = ENALL(42)
         END IF
C
C        - ZERO OFF SOME FORCES -
C
         IF(NACTMM.GT.0) THEN
            DO KOPT = 1, NACTMM
               IFFAT = LACTMM(KOPT)
               IF(IFFAT.GT.0)FFGRD2(1,IFFAT)=FFGRD2(1,IFFAT)+1.0D+03
            ENDDO
            DO IFFAT = 1, NFFAT
               IF(FFGRD2(1,IFFAT).GT.0.5D+03) THEN
                  FFGRD2(1,IFFAT) = FFGRD2(1,IFFAT)-1.0D+03
               ELSE
                  DO III = 1, 3
                     FFGRD2(III,IFFAT) = ZERO
                  ENDDO
               END IF
            ENDDO
         END IF
         DO KFIX=1,NFIXMM
            DO III = 1, 3
               FFGRD2(III,IFIXMM(KFIX)) = ZERO
            ENDDO
         ENDDO
C
C        -- PRINT COORDINATES --
  499    CONTINUE
         CALL TMDATE(TIMSTR)
         IF(MASWRK .AND. (MOD(ISTEP,KOUT).EQ.0.OR.ISTEP.EQ.NSTEP.OR.
     *                    CVGED))THEN
C
            WRITE(36,'(1X,A,I10)')
     *      'RESTART COORDINATES FOR QUANPOL AT OPT STEP',ISTEP
            WRITE(36,'(A,3(A,F15.10,1X))')
     *      ' $QUANPO ',
     *      'CENTX=',CENTX*TOANGS,
     *      'CENTY=',CENTY*TOANGS,
     *      'CENTZ=',CENTZ*TOANGS
            WRITE(36,*)'$END'
C
            WRITE(36,'(1X,A,I10,10X,A,F19.2,A)')
     *      '$FFDATA          !   OPT STEP',ISTEP
            WRITE(36,*)'COORDINATES  NUC                   X',
     *          '                   Y                   Z'
            DO IFFAT = 1, NFFAT
               WRITE(36,1000)ATMNAM(IFFAT),ZANF(IFFAT),
     *         CORD(1,IFFAT)*TOANGS,
     *         CORD(2,IFFAT)*TOANGS,
     *         CORD(3,IFFAT)*TOANGS
            ENDDO
            WRITE(36,*)'STOP'
C
            IF(NFFAT.LE.100) THEN
            WRITE(IW,'(1X,A,I10,10X,A,F19.2,A)')
     *      '$FFDATA          !   OPT STEP',ISTEP
            WRITE(IW,*)'COORDINATES  NUC                   X',
     *          '                   Y                   Z'
            DO IFFAT = 1, NFFAT
               WRITE(IW,1000)ATMNAM(IFFAT),ZANF(IFFAT),
     *         CORD(1,IFFAT)*TOANGS,
     *         CORD(2,IFFAT)*TOANGS,
     *         CORD(3,IFFAT)*TOANGS
            ENDDO
            WRITE(IW,*)'STOP'
            END IF
C
C           -- PRINT KOUTACT ATOMS --
C              HERE ONLY $FFDATA IS CONSIDERED, $FFDATB IS NOT.
C
            IF(KOUTACT(1).GT.0) THEN
            WRITE(IW,'(/A,I10)')'!KOUTACT ATOMS AROUND $FFDATA ATOM ',
     *      KOUTACT(1)
            WRITE(IW,'(1X,A,I10)')
     *      '$FFDATA ! KOUTACT     OPT STEP',ISTEP
            ACTX = CORD(1,KOUTACT(1))
            ACTY = CORD(2,KOUTACT(1))
            ACTZ = CORD(3,KOUTACT(1))
            ACTR = DBLE(KOUTACT(2))*1.0D-08*TOBOHR
            ACTR2= ACTR**2
            DO IFFAT=1,NFFAT
               CX    = CORD(1,IFFAT) - ACTX
               CY    = CORD(2,IFFAT) - ACTY
               CZ    = CORD(3,IFFAT) - ACTZ
               PBCX  = XBOX*ANINT(CX*ONEXBOX)
               PBCY  = YBOX*ANINT(CY*ONEYBOX)
               PBCZ  = ZBOX*ANINT(CZ*ONEZBOX)
               CX    = CX - PBCX
               CY    = CY - PBCY
               CZ    = CZ - PBCZ
               RCXYZ = CX**2 + CY**2 + CZ**2
               IF(RCXYZ.LE.ACTR2) THEN
                  WRITE(IW,1000)ATMNAM(IFFAT),ZANF(IFFAT),
     *            (CX+ACTX)*TOANGS,(CY+ACTY)*TOANGS,(CZ+ACTZ)*TOANGS
               END IF
            ENDDO
            WRITE(IW,*)'$END ! KOUTACT'
            END IF
C
            CALL TIMIT(1)
            CALL FLSHBF(IW)
            CALL FLSHBF(36)
         END IF
         IF(CVGED) GOTO 500
C
C        -- CALCULATE PROPERTIES
         NAT = 0
         IIISTEP=ISTEP
         CALL OPTPROP(IIISTEP,1)
         NAT = NFFAT
         IF(NACTMM.GT.0) NAT = NACTMM
         E= ENPOT
C
         IF(NACTMM.EQ.0) THEN
            CALL DCOPY(3*NFFAT,CORD,1,C,1)
            CALL DCOPY(3*NFFAT,FFGRD2,1,EG,1)
            CALL DCOPY(3*NFFAT,C,1,XX,1)
         ELSE
            DO I=1,NACTMM
               IFFAT = LACTMM(I)
               C(1,I)  = CORD(1,IFFAT)
               C(2,I)  = CORD(2,IFFAT)
               C(3,I)  = CORD(3,IFFAT)
               EG(1,I) = FFGRD2(1,IFFAT)
               EG(2,I) = FFGRD2(2,IFFAT)
               EG(3,I) = FFGRD2(3,IFFAT)
               XX(1,I) = C(1,I)
               XX(2,I) = C(2,I)
               XX(3,I) = C(3,I)
            ENDDO
         END IF
         PROJCT=.FALSE.
         KDIAGS= KDIAG
         KDIAG = 1    ! USE EVVRSP, WHICH IS PARALLEL
         CALL DISPLC(CVGED,CRDOLD,CRDNEW,NCVAL,NCOORD,
     *               NPUN,OUT,METHOD,ITBMAT,STPT,STSTEP,PROJCT,RUNTYP)
         KDIAG = KDIAGS
C
         CALL VALFM(LOADFM)
         LHESS   = LOADFM + 1
         LAST    = LHESS  + NCOORD*NCOORD
         NEED    = LAST   - LOADFM -1
         CALL GETFM(NEED)
         IPSV2=IP
         CALL DAREAD(IDAF,IODA,X(LHESS),NCOORD*NCOORD,4,0)
         IF(MOD(ISTEP,2).EQ.0) THEN
            CALL SEQCLO(39,'DELETE')
            CALL SEQOPN(39,'OPTHES1','NEW',.FALSE.,'FORMATTED')
            IP=39
         END IF
         IF(MOD(ISTEP,2).EQ.1) THEN
            CALL SEQCLO(40,'DELETE')
            CALL SEQOPN(40,'OPTHES2','NEW',.FALSE.,'FORMATTED')
            IP=40
         END IF
         IF(MASWRK) WRITE(IP,*) 'THE APPROXIMATE HESSIAN IS'
         CALL FCMPUN(X(LHESS),NCOORD)
         CALL FLSHBF(IP)
         IP=IPSV2
         CALL RETFM(NEED)
C
         IF(CVGED) THEN
            IF(MOD(ISTEP,KOUT).EQ.0.OR.ISTEP.EQ.NSTEP) THEN
               GOTO 500
            ELSE
               GOTO 499  ! MAKE IT UP
            END IF
         END IF
C
C        -- UPDATE THE NUCLEAR COORDINATES --
C
         IF(NSTEP.GT.0) THEN
         CALL DCOPY(NCVAL,EG,1,OLDF,1)
         CALL SYMDR(CRDNEW)
         IF(NACTMM.EQ.0) THEN
            DO 200 IFFAT =1,NFFAT
               IF((FFGRD2(1,IFFAT)+FFGRD2(2,IFFAT)+FFGRD2(3,IFFAT))
     *            .EQ.ZERO) GOTO 200
               C(1,IFFAT)=C(1,IFFAT) + CRDNEW(IFFAT*3-2)
               C(2,IFFAT)=C(2,IFFAT) + CRDNEW(IFFAT*3-1)
               C(3,IFFAT)=C(3,IFFAT) + CRDNEW(IFFAT*3  )
               CORD(1,IFFAT)=C(1,IFFAT)
               CORD(2,IFFAT)=C(2,IFFAT)
               CORD(3,IFFAT)=C(3,IFFAT)
  200       CONTINUE
         ELSE
            DO 202 I=1,NACTMM
               IFFAT = LACTMM(I)
               IF((FFGRD2(1,IFFAT)+FFGRD2(2,IFFAT)+FFGRD2(3,IFFAT))
     *            .EQ.ZERO) GOTO 202
               C(1,I)  = C(1,I) + CRDNEW(I*3-2)
               C(2,I)  = C(2,I) + CRDNEW(I*3-1)
               C(3,I)  = C(3,I) + CRDNEW(I*3  )
               CORD(1,IFFAT) = C(1,I)
               CORD(2,IFFAT) = C(2,I)
               CORD(3,IFFAT) = C(3,I)
  202       CONTINUE
         END IF
         END IF
C
C        -- SYNCHRONIZE CORD EVERY 200 STEPS --
         IF(MOD(ISTEP,200).EQ.0) THEN
            IF(GOPARR) THEN
               CALL DDI_BCAST(462,'F',CORD,3*NFFAT,MASTER)
            END IF
         END IF
C
      ENDDO
C
C     -- NOT LOCATED --
C
      IF(MASWRK) THEN
         WRITE(IW,*)' '
         WRITE(IW,*)'QUANPOL OPTIMIZATION FAILED. ',
     *              'PLEASE RESTART THE JOB USING THE BEST GEOMETRY.'
         WRITE(IW,*)' '
      END IF
      CALL FLSHBF(IW)
      CALL FLSHBF(36)
      CALL FLSHBF(39)
      CALL FLSHBF(40)
      CALL SEQCLO(36,'KEEP')
      CALL SEQCLO(39,'KEEP')
      CALL SEQCLO(40,'KEEP')
      RETURN
C
C     -- LOCATED --
C
  500 CONTINUE
      IF(MASWRK) THEN
         WRITE(IW,*)' '
         WRITE(IW,*)'================= QUANPOL OPTIMIZATION',
     *              ' SUCCESSFULLY COMPLETED ================='
         WRITE(IW,*)' '
      END IF
      CALL FLSHBF(IW)
      CALL FLSHBF(36)
      CALL FLSHBF(39)
      CALL FLSHBF(40)
      CALL SEQCLO(36,'KEEP')
      CALL SEQCLO(39,'KEEP')
      CALL SEQCLO(40,'KEEP')
      RETURN
 1000 FORMAT(1X,A10,1X,F5.1,1X,F19.13,1X,F19.13,1X,F19.13)
      END
C*MODULE QUANPOC  *DECK FFOPT4
!>
!> @brief    QM/MM optimization
!>
!> @author   Nandun Thellamurege
!>           - May 2011
!>
!> @details  BFGS Hessian method
!>
      SUBROUTINE FFOPT4(ATMNAM,CORD,CORDSV,ZANF,
     *                  ZMAS,
     *                  FFGRD2,QMGRD2,
     *                  LISTQM,NONLS1,NONLSTQ,
     *                  CRDNEW,CRDOLD,
     *                  LSTCELL,NONLS2,CORDSV2,CORDSVQ,
     *                  MVFASTS2,MVFASTS3,MVFASTS4,
     *                  MVFASTL2,MVFASTL3,MVFASTL4,
     *                  NONLSA,NONLSB,
     *                  NONLSPMA,NONLSPMB)
      USE MX_LIMITS,ONLY:MXFRG,mxatm,mxao,mxrt
C
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      LOGICAL GOPARR,DSKWRK,MASWRK
      LOGICAL CVGED,OUT,STPT,PROJCT,GOTEG,GOTEH
C
      PARAMETER (TOANGS=0.52917724924D+00)
      PARAMETER (TOBOHR=1.0D+00/TOANGS)
      PARAMETER (ZERO=0.0D+00)
      PARAMETER (PT5=0.50D+00)
      PARAMETER (MAXML=1024)
C
      CHARACTER*10  ATMNAM,TXT1,TXT2,TXT3,TXT4
C
      DOUBLE PRECISION METHOD
      LOGICAL LINEAR
C
      DIMENSION ATMNAM(NFFAT),CORD(3,NFFAT),CORDSV(3,NFFAT),ZANF(NFFAT),
     *          ZMAS(NFFAT),
     *          FFGRD2(3,NFFAT),QMGRD2(3,NAT),
     *          LISTQM(*),NONLSTQ(*),CRDNEW(*)
      DIMENSION TIMSTR(3)
C
      COMMON /DLCFRZ/ FVALUE(50),ITABLE(50),IFTYPE(50),NCONST
      COMMON /ENRGYS/ ENUCR,EELCT,ETOT,SZ,SZZ,ECORE,ESCF,EERD,E1,E2,
     *                VEN,VEE,EPOT,EKIN,ESTATE(MXRT),STATN,EDFT(3),EDISP
      COMMON /FFDIMR/ IDIMER,IBREAK(81),N1213JMM,NMOLE,MATOM(MAXML),
     *                MCHARG(MAXML),MMULT(MAXML),MELEC(MAXML),RDIMER
      COMMON /FFENGY/ EN12,EN123,EN123R4,EN123B4,EN234W1,ENCHAR,ENLJR,
     *                ENLJD,ENPOL,XENPOL,ENRXN,XENRXN,ENRXNPOL,ENRXNR,
     *                EN12312,ENQUANP(30),
     *                ENBIAS,ENCENT,ENUCCH,ENCMAP,ENPOT,ENKIN,ENTOT
      COMMON /FFFREE/ SOL1CH,SOL1LJ,SOL1IM,SOLFRE1MM,
     *                SOL2CH,SOL2LJ,SOL2IM,SOLFRE2MM,
     *                ASOL1CH,ASOL1LJ,ASOL1IM,ASOL1MM,
     *                ASOL2CH,ASOL2LJ,ASOL2IM,ASOL2MM,
     *                PMF1BD,PMF1AG,PMF1DR,PMF1DB,PMF1WG,PMF1CM,
     *                PMF1CH,PMF1PO,PMF1LJ,PMF1IM,PMFFRE1MM,
     *                EN12A,EN12B,EN123A,EN123B,EN123R4A,EN123R4B,
     *                EN234W1A,EN234W1B,EN123B4A,EN123B4B,
     *                EN12312A,EN12312B,
     *                ENCHARA,ENCHARB,ENPOLA,ENPOLB,
     *                ENRXNA,ENRXNB,ENRXNPOLA,ENRXNPOLB,
     *                ENRXNRA,ENRXNRB,ENLJRA,ENLJRB,ENLJDA,ENLJDB,
     *                ENCMAPA,ENCMAPB,ENFIXSOA,ENFIXSOB,
     *                ENCENTA,ENCENTB,WSIMUL,WPERT1,WPERT2,
     *                IFEPTYP
      COMMON /FFPARM/ NFFAT,NBOND,NANGL,NDIHR,NDIHB,NCMAP,NWAGG,
     *                N1213J,N14J,NLKQMM,IDOCHG,IDOPOL,IDOLJ,IDOCMAP
      COMMON /FFMDPA/ DT,DT2,TEMP0,PRES0,POLTOL,VIR(3),PMEAN,VOLAV,
     *                ENPAV,ENKAV,TEMPAV,BERENDT,BERENDP,VELMAX,
     *                PMEANX,PMEANY,PMEANZ,NSTEP,KMASTER,KOUTACT(2),
     *                IHESS,INTALG,ITSTAT,IPSTAT,JOUT,KOUT,LOUT
      COMMON /FFMPT3/ NACTMM,LACTMM(2020),LFFDIPOLD,JUMBUP,
     *                NACTQM,LACTQM(2020),LFFOLDC,LFFQMVELSV,MMHESS,
     *                LFFQMCHG,LFFQMCHGB,ISWAP,R2SWAP,DFTBMM
      COMMON /FFPBSW/ XBOX,YBOX,ZBOX,SWF,SWFDX,SWFDY,SWFDZ,
     *                SWRA,ONESWRA,SWRA2,ONESWRA2,
     *                SWRB,ONESWRB,SWRB2,ONESWRB2,ONESWRB4,
     *                SWFDUM3,SWFDUM4,SWFDUM5,
     *                SWRAQ,ONESWRAQ,SWRAQ2,ONESWRAQ2,
     *                SWRBQ,ONESWRBQ,SWRBQ2,ONESWRBQ2,ONESWRBQ4,
     *                SWFDUM3Q,SWFDUM4Q,SWFDUM5Q,
     *                QMSIZE,QMCX,QMCY,QMCZ,QMCXSV,QMCYSV,QMCZSV,
     *                CENTX,CENTY,CENTZ,BUFWID1,BUFWID2,RDAMP,
     *                EFIELDX,EFIELDY,EFIELDZ,QMCXSV2,QMCYSV2,QMCZSV2,
     *                EPS1RB,EPS1RB3,ONEXBOX,ONEYBOX,ONEZBOX,
     *                LQMCT,MXLIST1,NTODO,NTODOSV,NTODOQ,
     *                ISWITCH,ISHIFT,IPOLSHF,
     *                LFFLSTCELL,LFFCORDSV,
     *                LFFPOLSV,LFFCORDSV2,LFFNONLS2,LFFCORDSVQ,
     *                LFFMVFASTS2,LFFMVFASTS3,LFFMVFASTS4,
     *                LFFMVFASTL2,LFFMVFASTL3,LFFMVFASTL4,
     *                MXCHECK,MXLIST2,NTODO2,NTODO2SV
      COMMON /FFQMPA/ ENFFQM2,ENPAV2,
     *                SCFTYP2,TDDFT2,MPLEVL2,CITYP2,
     *                ICHARG2,MULT2,IDOQM2,IREDOX,IQMPKA,IQMRXN,
     *                MATOMA,MCHARGA,MULTA,MELEA,
     *                MATOMB,MCHARGB,MULTB,MELEB,
     *                IECPX,NSHELLX,IMP,JMP,ICORSH,IGTF,
     *                LFFZANX,LFFCLPX,LFFZLPX,LFFNLPX,LFFKFRSTX,
     *                LFFKLASTX,LFFLMAXX,LFFLPSKIPX,LFFIZCOREX,
     *                LFFCX,LFFIANX,LFFEXX,LFFCSX,LFFCPX,LFFCDX,
     *                LFFCFX,LFFCGX,LFFCHX,LFFCIX,LFFKSTARTX,
     *                LFFKATOMX,LFFKTYPEX,LFFKNGX,LFFKLOCX,
     *                LFFMINX,LFFMAXX,LFFMPTYPX,LFFAN0X,
     *                LFFALPN0X,LFFAN1X,LFFALPN1X,LFFMPSKPX,
     *                LFFNOAN0X,LFFNOAN1X,LFFBPARX,LFFEXPMPX,
     *                LFFCSMPX,LFFCPMPX,LFFCDMPX,LFFCFMPX,LFFMPSKIPX,
     *                LFFNOCOSHX,LFFMPKSTAX,LFFMPKNGX,LFFMPKTYPX,
     *                LFFMPKMINX,LFFMPKMAXX,LFFMPKLOCX,LFFANAMX
      COMMON /FFRMSD/ DIPT(8),TIMGYRA,TIMRALL,
     *                NATPDB,NGYRA,NDIEL,NRALL,
     *                LFFRALL0,NRMSD,LFFRMSD0,KFREEAB(201),
     *                NRIJMM,IJRMM(2,100),NRIJQM,IJRQM(2,100),
     *                NAIJKMM,IJKMM(3,100),NAIJKQM,IJKQM(3,100),
     *                NFIXMM,IFIXMM(200),NFIXQM,IFIXQM(200)
      COMMON /FFRXN / RXNEPS,RSPHSOL,ISPHSOL
      COMMON /FFSPH / SPHRAD,SPHSIG,SPHEPS,IADDWAT
      COMMON /FFUMBR/ UMBFC,UMBR0,UMBSIZE,
     *                NUMBBIN,NUMBATM(6),NUMBTYP,LFFUMBHIS,
     *                UM2FC,UM2R0,UM2SIZE,
     *                NUM2BIN,NUM2ATM(6),NUM2TYP,LFFUM2HIS
      COMMON /FMCOM / X(1)
      COMMON /FUNCT / E,EG(3,MXATM)
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /IOFILE/ IR,IW,IP,IJK,IJKT,IDAF,NAV,IODA(950)
      COMMON /MACHSW/ KDIAG,ICORFL,IXDR,modio,mem10,lpnt10,mem10m
      COMMON /MASSES/ ZMASS(MXATM)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /RUNLAB/ TITLE(10),ANAM(MXATM),BNAM(MXATM),BFLAB(MXAO)
      COMMON /RUNOPT/ RUNTYP,EXETYP,NEVALS,NGLEVL,NHLEVL
      COMMON /OPTEF / D(3*MXATM+6*MXFRG),OLDF(3*MXATM+6*MXFRG),
     *                VMODE(3*MXATM+6*MXFRG),RADIUS
      COMMON /OPTGRD/ XX(3,MXATM),ES,FE(20),
     *                CONVF,FMAXT,DXMAXT,RMAX,RMIN,RLIM,
     *                EIGMAX,EIGMIN,GRDERR,FRMS,FMAX,TRMAX,TRMIN,
     *                IC(20),MSTEP,NSERCH,NPMAX,NP,IFOLOW,
     *                NNEG,IUPHSS,IEXIT,ITRUPD,IPAD,KDIAGH,NPRICO
      COMMON /SCFOPT/ CONVHF,MAXIT,MCONV,NPUNCH,NPREO(4),FSHIFT
      COMMON /ZMAT  / NZMAT,NZVAR,NVAR,NSYMC,LINEAR
      COMMON /ZMTALT/ NZMAT2,NZVAR2,NVAR2,NZMTRD,ICOORD
      COMMON /FROZNC/ IFREEZ(3*MXATM),NFRZ
      COMMON /OUTPUT/ NPRINT,ITOL,ICUT,NORMF,NORMP,NOPK
C
      DATA BLANK/8H        /
      DATA OPTMIZ   /8HOPTIMIZE/
C
C     NANDUN THELLAMUREGE, MAY 2011, LINCOLN
C     HUI LI, MAY 6, 2014, USE LACTMM AND PRJGRD
C
      IF(NSTEP.LT.0) RETURN
C
      CALL SEQOPN(36,'TRAJECT','NEW',.FALSE.,'FORMATTED')
      CALL SEQOPN(39,'OPTHES1','NEW',.FALSE.,'FORMATTED')
      CALL SEQOPN(40,'OPTHES2','NEW',.FALSE.,'FORMATTED')
      IF(IADDWAT.EQ.2 .AND. ISPHSOL.GT.0) RSPHSOL=RSPHSOL*1.0D+40
C
      NPRINT  = -5
      NATSV   = NAT
      IF(NACTMM.EQ.0) THEN
         DO IFFAT =1,NFFAT
            IF(LISTQM(IFFAT).EQ.0.OR.LISTQM(IFFAT).GT.NATSV) THEN
               NAT = NAT + 1
               C(1,NAT)=CORD(1,IFFAT)
               C(2,NAT)=CORD(2,IFFAT)
               C(3,NAT)=CORD(3,IFFAT)
               ZAN(NAT)=ZANF(IFFAT)
               READ(UNIT=ATMNAM(IFFAT)(1:8),FMT='(A8)') ANAM(NAT)
               BNAM(NAT)=BLANK
               ZMASS(NAT)=ZMAS(IFFAT)/1822.88850204D+00
            END IF
         ENDDO
      ELSE
         DO I=1,NACTMM
            IFFAT  = LACTMM(I)
            IF(LISTQM(IFFAT).EQ.0.OR.LISTQM(IFFAT).GT.NATSV) THEN
               NAT = NAT + 1
               C(1,NAT) = CORD(1,IFFAT)
               C(2,NAT) = CORD(2,IFFAT)
               C(3,NAT) = CORD(3,IFFAT)
               ZAN(NAT) = ZANF(IFFAT)
               READ(UNIT=ATMNAM(IFFAT)(1:8),FMT='(A8)') ANAM(NAT)
               BNAM(NAT)=BLANK
               ZMASS(NAT)=ZMAS(IFFAT)/1822.88850204D+00
            END IF
         ENDDO
      END IF
      IF(NAT.GT.2000) THEN
         IF(MASWRK)WRITE(IW,*)
     *   'ERROR: NUMBER OF OPTIMIZE ATOMS EXCEEDED 2000'
         IF(MASWRK)WRITE(IW,*)' '
         CALL ABRT
      END IF
C
      NATTOT=NAT
      NCOORD=NATTOT*3
      NCVAL =NCOORD
C
C     -- SET SOME CONTROL VARIABLES --
      MSTEP  = NSTEP
      NZMAT  = 0
      NFRZ   = 0
      NCONST = 0
      ICOORD = 8
      CVGED  = .FALSE.
      MODE   = 0
      OUT    = .FALSE.
      MAX2   = NCVAL
C
      DO IAT =1,NATSV
         IF(LISTQM(NFFAT+IAT).GT.0) THEN
            IFFAT       = LISTQM(NFFAT+IAT)
            ZDUM        = ZAN(IAT)
            ZAN(IAT)    = ZANF(IFFAT)
            ZANF(IFFAT) = ZDUM
         END IF
      ENDDO
      CALL SIGINI(MODE,RUNTYP,NCVAL,NCOORD,METHOD,OUT,
     *            GOTEG,NPRT,NPUN,ITBMAT,STPT,STSTEP,
     *            PROJCT,DUMMY,DUMMY,MAX2,NPRTHS)
      DO IAT =1,NATSV
         IF(LISTQM(NFFAT+IAT).GT.0) THEN
            IFFAT       = LISTQM(NFFAT+IAT)
            ZDUM        = ZAN(IAT)
            ZAN(IAT)    = ZANF(IFFAT)
            ZANF(IFFAT) = ZDUM
         END IF
      ENDDO
C
C     -- CHECK THE $HESS E(NUC) --
C
      CALL SEQREW(IR)
      CALL FNDGRP(IR,' $HESS  ',IEOF)
      IF (IEOF.EQ.1) THEN
         GOTEH=.FALSE.
      ELSE
         GOTEH=.TRUE.
         IF (MASWRK) THEN
            READ(IR,*)TXT1,TXT2,AEN1,TXT3,TXT4,AEN2
            IF(TXT1.EQ.'ENERGY'.AND.TXT2.EQ.'IS'.AND.
     *         TXT3.EQ.'E(NUC)'.AND.TXT4.EQ.'IS') THEN
            IF(AEN2.EQ.ZERO.AND.MMHESS.EQ.0) THEN
               WRITE(IW,'(/1X,A,/)')'ERROR: MMHESS=1 SHOULD BE USED.'
               CALL ABRT
            END IF
            IF(AEN1.LT.ZERO.AND.AEN2.GT.ZERO.AND.MMHESS.EQ.1) THEN
               WRITE(IW,'(/1X,A,/)')'ERROR: MMHESS=0 SHOULD BE USED.'
               CALL ABRT
            END IF
            END IF
         END IF
      END IF
C
C     -- RE-ORDER THE READ IN HESSIAN IN CASE OF MM TO QM/MM --
C
      IF(MMHESS.EQ.1) THEN
         IF(LISTQM(NFFAT+NATSV+1).LT.NATSV) THEN
            IF(MASWRK) WRITE(IW,'(/1X,A,A/)')
     *      'ERROR: WHEN MM $HESS IS USED FOR QM/MM, ALL QM ATOMS ',
     *      'MUST BE IN $FFDATA.'
            CALL ABRT
         END IF
         CALL VALFM(LOADFM)
         LHESB   = LOADFM + 1
         LHESS   = LHESB  + NCOORD*NCOORD
         LAST    = LHESS  + NCOORD*NCOORD
         NEED    = LAST   - LOADFM -1
         CALL GETFM(NEED)
C
         GOTEH = .FALSE.
         CALL FCMIN(X(LHESS),NCOORD,GOTEH)
         IF(.NOT.GOTEH) THEN
            IF(MASWRK) WRITE(IW,'(/1X,A/)')
     *      'ERROR: MMHESS=1 MUST READ IN $HESS FROM INPUT FILE.'
            CALL ABRT
         END IF
C
         MFFAT1 = NATSV
         IF(NACTMM.EQ.0) THEN
            DO IFFAT =1,NFFAT
               IF(LISTQM(IFFAT).EQ.0.OR.LISTQM(IFFAT).GT.NATSV) THEN
                  MFFAT1 = MFFAT1 + 1
                  MFFAT2 = NATSV
                  DO JFFAT =1,NFFAT
                     IF(LISTQM(JFFAT).EQ.0.OR.LISTQM(JFFAT).GT.NATSV)
     *                                                           THEN
                        MFFAT2 = MFFAT2 + 1
                        DO III = 1, 3
                        DO JJJ = 1, 3
              X(LHESB+((MFFAT1-1)*3+III-1)*NCOORD+(MFFAT2-1)*3+JJJ-1)
     *      = X(LHESS+((IFFAT -1)*3+III-1)*NCOORD+(JFFAT -1)*3+JJJ-1)
                        ENDDO
                        ENDDO
                     ELSE
                        MAT2 = LISTQM(JFFAT)
                        DO III = 1, 3
                        DO JJJ = 1, 3
              X(LHESB+((MFFAT1-1)*3+III-1)*NCOORD+(  MAT2-1)*3+JJJ-1)
     *      = X(LHESS+((IFFAT -1)*3+III-1)*NCOORD+( JFFAT-1)*3+JJJ-1)
                        ENDDO
                        ENDDO
                     END IF
                  ENDDO
               ELSE
                  MAT1 = LISTQM(IFFAT)
                  MFFAT2 = NATSV
                  DO JFFAT =1,NFFAT
                     IF(LISTQM(JFFAT).EQ.0.OR.LISTQM(JFFAT).GT.NATSV)
     *                                                           THEN
                        MFFAT2 = MFFAT2 + 1
                        DO III = 1, 3
                        DO JJJ = 1, 3
              X(LHESB+((  MAT1-1)*3+III-1)*NCOORD+(MFFAT2-1)*3+JJJ-1)
     *      = X(LHESS+((IFFAT -1)*3+III-1)*NCOORD+( JFFAT-1)*3+JJJ-1)
                        ENDDO
                        ENDDO
                     ELSE
                        MAT2 = LISTQM(JFFAT)
                        DO III = 1, 3
                        DO JJJ = 1, 3
              X(LHESB+((  MAT1-1)*3+III-1)*NCOORD+( MAT2-1)*3+JJJ-1)
     *      = X(LHESS+((IFFAT -1)*3+III-1)*NCOORD+(JFFAT-1)*3+JJJ-1)
                        ENDDO
                        ENDDO
                     END IF
                  ENDDO
               END IF
            ENDDO
         END IF
         IF(NACTMM.GT.0) THEN
            DO I=1,NACTMM
               IFFAT  = LACTMM(I)
               IF(LISTQM(IFFAT).EQ.0.OR.LISTQM(IFFAT).GT.NATSV) THEN
                  MFFAT1 = MFFAT1 + 1
                  MFFAT2 = NATSV
                  DO J=1,NACTMM
                     JFFAT  = LACTMM(J)
                     IF(LISTQM(JFFAT).EQ.0.OR.LISTQM(JFFAT).GT.NATSV)
     *                                                           THEN
                        MFFAT2 = MFFAT2 + 1
                        DO III = 1, 3
                        DO JJJ = 1, 3
              X(LHESB+((MFFAT1-1)*3+III-1)*NCOORD+(MFFAT2-1)*3+JJJ-1)
     *      = X(LHESS+((    I -1)*3+III-1)*NCOORD+(     J-1)*3+JJJ-1)
                        ENDDO
                        ENDDO
                     ELSE
                        MAT2 = LISTQM(JFFAT)
                        DO III = 1, 3
                        DO JJJ = 1, 3
              X(LHESB+((MFFAT1-1)*3+III-1)*NCOORD+(  MAT2-1)*3+JJJ-1)
     *      = X(LHESS+((    I -1)*3+III-1)*NCOORD+(     J-1)*3+JJJ-1)
                        ENDDO
                        ENDDO
                     END IF
                  ENDDO
               ELSE
                  MAT1   = LISTQM(IFFAT)
                  MFFAT2 = NATSV
                  DO J=1,NACTMM
                     JFFAT  = LACTMM(J)
                     IF(LISTQM(JFFAT).EQ.0.OR.LISTQM(JFFAT).GT.NATSV)
     *                                                           THEN
                        MFFAT2 = MFFAT2 + 1
                        DO III = 1, 3
                        DO JJJ = 1, 3
              X(LHESB+((  MAT1-1)*3+III-1)*NCOORD+(MFFAT2-1)*3+JJJ-1)
     *      = X(LHESS+((    I -1)*3+III-1)*NCOORD+(     J-1)*3+JJJ-1)
                        ENDDO
                        ENDDO
                     ELSE
                        MAT2 = LISTQM(JFFAT)
                        DO III = 1, 3
                        DO JJJ = 1, 3
              X(LHESB+((  MAT1-1)*3+III-1)*NCOORD+( MAT2-1)*3+JJJ-1)
     *      = X(LHESS+((    I -1)*3+III-1)*NCOORD+(    J-1)*3+JJJ-1)
                        ENDDO
                        ENDDO
                     END IF
                  ENDDO
               END IF
            ENDDO
         END IF
         CALL DAWRIT(IDAF,IODA,X(LHESB),NCOORD*NCOORD,4,0)
         CALL RETFM(NEED)
      END IF
C
C     -- READY TO TAKE OFF ?
C
      DO ISTEP = 0, NSTEP
         NSERCH = ISTEP
         IF(ISTEP.EQ.50.AND.IADDWAT.EQ.2 .AND. ISPHSOL.GT.0)
     *   RSPHSOL = RSPHSOL*1.0D-40
C
         NPUN = -1
         NPUNCH = 0
         IF(MOD(ISTEP,KOUT).EQ.0 .OR. ISTEP.EQ.NSTEP) NPUNCH = 0
         IF(ISTEP.EQ.0) NPUNCH = 2
C
C        -- CALCULATE ENERGY AND GRADIENT -
         NAT    = NATSV    !  BECAUSE NONBOND USES NAT
         CALL NONBOND(ISTEP,CORD,CORDSV,CORDSV2,CORDSVQ,
     *                NONLS1,NONLS2,
     *                NONLSTQ,LSTCELL,
     *                MVFASTS2,MVFASTS3,MVFASTS4,
     *                MVFASTL2,MVFASTL3,MVFASTL4,
     *                NONLSA,NONLSB,
     *                NONLSPMA,NONLSPMB)
         CALL VCLR(QMGRD2,1,3*NATSV)
         CALL VCLR(FFGRD2,1,3*NFFAT)
         CALL VCLR(VIR,1,3)    !  NEED TO DO THIS - HUI LI
         IF(IDIMER.EQ.0)THEN
            CALL GRADX
         ELSE
            CALL DIMERX
         END IF
C        -- ESCF WAS ZERO IF SCF NOT CONVERGED --
         IF(ABS(ESCF).LT.1.0D-12) THEN
            IF(MASWRK) WRITE(IW,'(/A/)')
     *         ' ERROR: SCF WAS NOT CONVERGED. QUANPOL MUST STOP.'
            CALL ABRT
         END IF
         CALL DCOPY(3*NATSV,EG,1,QMGRD2,1)
         IF(IFEPTYP.GT.0.AND.MATOMB.GT.0) THEN
            CALL SAVEABPROP(1)
            CALL DCOPY(3*MATOMA,C,1,X(LFFCX),1)
            CALL SAVEFFDATA
            CALL SETQMAB(2)
            CALL SETFFDATAB(2)
            CALL VCLR(QMGRD2,1,3*NAT  )
            CALL VCLR(FFGRD2,1,3*NFFAT)
            CALL VCLR(VIR   ,1,3)
            CALL GRADX
C           -- ESCF WAS ZERO IF SCF NOT CONVERGED --
            IF(ABS(ESCF).LT.1.0D-12) THEN
               IF(MASWRK) WRITE(IW,'(/A/)')
     *            ' ERROR: SCF WAS NOT CONVERGED. QUANPOL MUST STOP.'
               CALL ABRT
            END IF
            CALL DCOPY(3*NAT,EG,1,QMGRD2,1)
            CALL SAVEABPROP(2)
            CALL SETQMAB(1)
            CALL SETFFDATAB(1)
            CALL MIXABPROP
         END IF
         E      = ETOT
         NAT    = NATTOT
C
C        -- PRINT COORDINATES --
C
  499    CONTINUE
         CALL TMDATE(TIMSTR)
         IF(MASWRK .AND. (MOD(ISTEP,KOUT).EQ.0.OR.ISTEP.EQ.NSTEP.OR.
     *                    CVGED))THEN
C
            WRITE(36,'(1X,A,I10)')
     *      'RESTART COORDINATES FOR QUANPOL AT OPT STEP',ISTEP
            WRITE(36,'(A,3(A,F15.10,1X),A)')
     *      ' $QUANPO ',
     *      'CENTX=',CENTX*TOANGS,
     *      'CENTY=',CENTY*TOANGS,
     *      'CENTZ=',CENTZ*TOANGS,
     *      '$END'
C
            WRITE(36,*)'$DATA'
            WRITE(36,'(A,3(F8.3,1X),A,I3,A,I10)')
     *      ' QM CENTER = ',QMCX*TOANGS,QMCY*TOANGS,QMCZ*TOANGS,
     *      'NEAR ATOM',LQMCT,'     AT OPT STEP',ISTEP
            WRITE(36,*)'C1'
            DO IAT = 1, NATSV
               WRITE(36,999)ANAM(IAT),ZAN(IAT),
     *         C(1,IAT)*TOANGS,
     *         C(2,IAT)*TOANGS,
     *         C(3,IAT)*TOANGS
            ENDDO
            IF(IFEPTYP.GT.0)THEN
               DO IATB = 1, MATOMB
                  WRITE(36,999)ANAM(MATOMA+IATB),  ! ANAM IS THE ORIGINAL
     *            X(LFFZANX +MATOMA+IATB-1),
     *            X(LFFCX+3*(MATOMA+IATB-1)  )*TOANGS,
     *            X(LFFCX+3*(MATOMA+IATB-1)+1)*TOANGS,
     *            X(LFFCX+3*(MATOMA+IATB-1)+2)*TOANGS
               ENDDO
            END IF
            WRITE(36,*)'$END'
            WRITE(36,'(1X,A,I10,10X,A,F19.2,A)')
     *      '$FFDATA          !   OPT STEP',ISTEP
            IF(IFEPTYP.EQ.1) THEN
               WRITE(36,'(1X,A,I10,10X,A,F19.2,A)')
     *         '$FFDATB          !   OPT STEP',ISTEP
            END IF
            WRITE(36,*)'COORDINATES  NUC                   X',
     *          '                   Y                   Z'
            DO IFFAT = 1, NFFAT
               WRITE(36,1000)ATMNAM(IFFAT),ZANF(IFFAT),
     *         CORD(1,IFFAT)*TOANGS,
     *         CORD(2,IFFAT)*TOANGS,
     *         CORD(3,IFFAT)*TOANGS
            ENDDO
            WRITE(36,*)'STOP'
            IF(NFFAT.LE.100) THEN
            WRITE(IW,*)'$DATA'
            WRITE(IW,'(A,3(F8.3,1X),A,I3,A,I10)')
     *      ' QM CENTER = ',QMCX*TOANGS,QMCY*TOANGS,QMCZ*TOANGS,
     *      'NEAR ATOM',LQMCT,'     AT OPT STEP',ISTEP
            WRITE(IW,*)'C1'
            DO IAT = 1, NATSV
               WRITE(IW,999)ANAM(IAT),ZAN(IAT),
     *         C(1,IAT)*TOANGS,
     *         C(2,IAT)*TOANGS,
     *         C(3,IAT)*TOANGS
            ENDDO
            IF(IFEPTYP.GT.0)THEN
               DO IATB = 1, MATOMB
                  WRITE(IW,999)ANAM(MATOMA+IATB),  ! ANAM IS THE ORIGINAL
     *            X(LFFZANX +MATOMA+IATB-1),
     *            X(LFFCX+3*(MATOMA+IATB-1)  )*TOANGS,
     *            X(LFFCX+3*(MATOMA+IATB-1)+1)*TOANGS,
     *            X(LFFCX+3*(MATOMA+IATB-1)+2)*TOANGS
               ENDDO
            END IF
            WRITE(IW,*)'$END'
            WRITE(IW,'(1X,A,I10,10X,A,F19.2,A)')
     *      '$FFDATA          !   OPT STEP',ISTEP
            IF(IFEPTYP.EQ.1) THEN
               WRITE(IW,'(1X,A,I10,10X,A,F19.2,A)')
     *         '$FFDATB          !   OPT STEP',ISTEP
            END IF
            WRITE(IW,*)'COORDINATES  NUC                   X',
     *          '                   Y                   Z'
            DO IFFAT = 1, NFFAT
               WRITE(IW,1000)ATMNAM(IFFAT),ZANF(IFFAT),
     *         CORD(1,IFFAT)*TOANGS,
     *         CORD(2,IFFAT)*TOANGS,
     *         CORD(3,IFFAT)*TOANGS
            ENDDO
            WRITE(IW,*)'STOP'
            END IF
C
            IF(NFFAT.GT.100.AND.NACTMM.GT.0) THEN
            WRITE(IW,*)'$DATA         ! THESE ARE LACTMM ATOMS'
            WRITE(IW,'(A,3(F8.3,1X),A,I3,A,I10)')
     *      ' QM CENTER = ',QMCX*TOANGS,QMCY*TOANGS,QMCZ*TOANGS,
     *      'NEAR ATOM',LQMCT,'     AT OPT STEP',ISTEP
            WRITE(IW,*)'C1'
            DO IAT = 1, NAT
               ZCHG = ZAN(IAT)
               IF(LISTQM(NFFAT+IAT).GT.0)THEN
                  ZCHG = ZANF(LISTQM(NFFAT+IAT))
               END IF
               CX    = C(1,IAT) - QMCX
               CY    = C(2,IAT) - QMCY
               CZ    = C(3,IAT) - QMCZ
               PBCX  = XBOX*ANINT(CX*ONEXBOX)
               PBCY  = YBOX*ANINT(CY*ONEYBOX)
               PBCZ  = ZBOX*ANINT(CZ*ONEZBOX)
               CX    = CX - PBCX + QMCX
               CY    = CY - PBCY + QMCY
               CZ    = CZ - PBCZ + QMCZ
               WRITE(IW,999)ANAM(IAT),ZCHG,
     *         CX*TOANGS,
     *         CY*TOANGS,
     *         CZ*TOANGS
            ENDDO
            WRITE(IW,*)'$END'
            END IF
C
C           -- PRINT KOUTACT ATOMS --
C              HERE ONLY $FFDATA IS CONSIDERED, $FFDATB IS NOT.
C
            IF(KOUTACT(1).GT.0) THEN
            WRITE(IW,'(/A,I10)')'!KOUTACT ATOMS AROUND $FFDATA ATOM ',
     *      KOUTACT(1)
            WRITE(IW,'(1X,A,I10)')
     *      '$FFDATA ! KOUTACT     OPT STEP',ISTEP
            ACTX = CORD(1,KOUTACT(1))
            ACTY = CORD(2,KOUTACT(1))
            ACTZ = CORD(3,KOUTACT(1))
            ACTR = DBLE(KOUTACT(2))*1.0D-08*TOBOHR
            ACTR2= ACTR**2
            DO IFFAT=1,NFFAT
               CX    = CORD(1,IFFAT) - ACTX
               CY    = CORD(2,IFFAT) - ACTY
               CZ    = CORD(3,IFFAT) - ACTZ
               PBCX  = XBOX*ANINT(CX*ONEXBOX)
               PBCY  = YBOX*ANINT(CY*ONEYBOX)
               PBCZ  = ZBOX*ANINT(CZ*ONEZBOX)
               CX    = CX - PBCX
               CY    = CY - PBCY
               CZ    = CZ - PBCZ
               RCXYZ = CX**2 + CY**2 + CZ**2
               IF(RCXYZ.LE.ACTR2) THEN
                  WRITE(IW,1000)ATMNAM(IFFAT),ZANF(IFFAT),
     *            (CX+ACTX)*TOANGS,(CY+ACTY)*TOANGS,(CZ+ACTZ)*TOANGS
               END IF
            ENDDO
            WRITE(IW,*)'$END ! KOUTACT'
            END IF
C
            CALL TIMIT(1)
            CALL FLSHBF(IW)
            CALL FLSHBF(36)
         END IF
         IF(CVGED) GOTO 500
C
C        - COMBINE QM AND MM GRADIENTS -
C
         IF(LISTQM(NFFAT+NATSV+1).GT.0)THEN
            DO IAT = 1, NATSV
               KFFAT = LISTQM(NFFAT+IAT)
               IF(KFFAT.GT.0) THEN
                  QMGRD2(1,IAT)   = QMGRD2(1,IAT) + FFGRD2(1,KFFAT)
                  QMGRD2(2,IAT)   = QMGRD2(2,IAT) + FFGRD2(2,KFFAT)
                  QMGRD2(3,IAT)   = QMGRD2(3,IAT) + FFGRD2(3,KFFAT)
                  EG(1,IAT)       = QMGRD2(1,IAT)
                  EG(2,IAT)       = QMGRD2(2,IAT)
                  EG(3,IAT)       = QMGRD2(3,IAT)
                  FFGRD2(1,KFFAT) = QMGRD2(1,IAT)
                  FFGRD2(2,KFFAT) = QMGRD2(2,IAT)
                  FFGRD2(3,KFFAT) = QMGRD2(3,IAT)
               END IF
            ENDDO
         END IF
C
C        - ZERO OFF SOME FORCES -
C
         IF(NACTQM.GT.0) THEN
            DO KOPT = 1, NACTQM
               IAT = LACTQM(KOPT)
               IF(IAT.GT.0)QMGRD2(1,IAT)=QMGRD2(1,IAT)+1.0D+03
            ENDDO
            DO IAT = 1, NATSV
               IF(QMGRD2(1,IAT).GT.0.5D+03) THEN
                  QMGRD2(1,IAT) = QMGRD2(1,IAT)-1.0D+03
               ELSE
                  DO III = 1, 3
                     QMGRD2(III,IAT) = ZERO
                     EG(III,IAT)     = ZERO
                  ENDDO
               END IF
            ENDDO
         END IF
         DO KFIX=1,NFIXQM
            DO III = 1, 3
               QMGRD2(III,IFIXQM(KFIX)) = ZERO
               EG(III,IFIXQM(KFIX))     = ZERO
            ENDDO
         ENDDO
         IF(NACTMM.GT.0) THEN
            DO KOPT = 1, NACTMM
               IFFAT = LACTMM(KOPT)
               IF(IFFAT.GT.0)FFGRD2(1,IFFAT)=FFGRD2(1,IFFAT)+1.0D+03
            ENDDO
            DO IFFAT = 1, NFFAT
               IF(FFGRD2(1,IFFAT).GT.0.5D+03) THEN
                  FFGRD2(1,IFFAT) = FFGRD2(1,IFFAT)-1.0D+03
               ELSE
                  DO III = 1, 3
                     FFGRD2(III,IFFAT) = ZERO
                  ENDDO
               END IF
            ENDDO
         END IF
         DO KFIX=1,NFIXMM
            DO III = 1, 3
               FFGRD2(III,IFIXMM(KFIX)) = ZERO
            ENDDO
         ENDDO
C
C        - EXPAND EG(3,*) TO INCLUDE MM ATOMS -
         IAT = NATSV
         IF(NACTMM.EQ.0) THEN
            DO IFFAT = 1, NFFAT
               IF(LISTQM(IFFAT).EQ.0.OR.LISTQM(IFFAT).GT.NATSV)THEN
                  IAT = IAT + 1
                  EG(1,IAT) = FFGRD2(1,IFFAT)
                  EG(2,IAT) = FFGRD2(2,IFFAT)
                  EG(3,IAT) = FFGRD2(3,IFFAT)
               END IF
            ENDDO
         ELSE
            DO I=1,NACTMM
               IFFAT  = LACTMM(I)
C              - LACTMM DOES NOT AFFECT PURE QM EG()
C              - LACTMM DOES NOT AFFECT QM(MM)  EG()
               IF(LISTQM(IFFAT).EQ.0.OR.LISTQM(IFFAT).GT.NATSV) THEN
                  IAT = IAT + 1
                  EG(1,IAT) = FFGRD2(1,IFFAT)
                  EG(2,IAT) = FFGRD2(2,IFFAT)
                  EG(3,IAT) = FFGRD2(3,IFFAT)
               END IF
            ENDDO
         END IF
C
C        - CALCULATE PROPERTIES -
         IIISTEP=ISTEP
         CALL OPTPROP(IIISTEP,1)
C
C        - CALCULATE DISPLACEMENT -
C
         CALL DCOPY(3*NATTOT,C,1,XX,1)
         PROJCT=.FALSE.
         KDIAGS= KDIAG
         KDIAG = 1    ! USE EVVRSP, WHICH IS PARALLEL
         CALL DISPLC(CVGED,CRDOLD,CRDNEW,NCVAL,NCOORD,
     *               NPUN,OUT,METHOD,ITBMAT,STPT,STSTEP,PROJCT,RUNTYP)
         KDIAG = KDIAGS
C
         CALL VALFM(LOADFM)
         LHESS   = LOADFM + 1
         LAST    = LHESS  + NCOORD*NCOORD
         NEED    = LAST   - LOADFM -1
         CALL GETFM(NEED)
         IPSV2=IP
         CALL DAREAD(IDAF,IODA,X(LHESS),NCOORD*NCOORD,4,0)
         IF(MOD(ISTEP,2).EQ.0) THEN
            CALL SEQCLO(39,'DELETE')
            CALL SEQOPN(39,'OPTHES1','NEW',.FALSE.,'FORMATTED')
            IP=39
         END IF
         IF(MOD(ISTEP,2).EQ.1) THEN
            CALL SEQCLO(40,'DELETE')
            CALL SEQOPN(40,'OPTHES2','NEW',.FALSE.,'FORMATTED')
            IP=40
         END IF
         IF(MASWRK) WRITE(IP,*) 'THE APPROXIMATE HESSIAN IS'
         CALL FCMPUN(X(LHESS),NCOORD)
         CALL FLSHBF(IP)
         IP=IPSV2
         CALL RETFM(NEED)
C
C        - RESET JUMBPOT=12 R0 -
         IF(RUNTYP.EQ.OPTMIZ.AND.NUMBTYP.EQ.12.AND.JUMBUP.NE.0) THEN
            IF(ISTEP.GT.0.AND.MOD(ISTEP,8).EQ.0) THEN
               N1    = NUMBATM(1)
               N2    = NUMBATM(2)
               CX    = CORD(1,N1) - CORD(1,N2)
               CY    = CORD(2,N1) - CORD(2,N2)
               CZ    = CORD(3,N1) - CORD(3,N2)
               R2    = CX*CX + CY*CY + CZ*CZ
               R     = SQRT(R2)
               IF(JUMBUP.GT.0)UMBR0=UMBR0+(UMBR0 - R)
               IF(JUMBUP.LT.0)UMBR0=R    +(R - UMBR0)
               IF(MASWRK) THEN
                  WRITE(IW,'(/1X,A,F10.6/)')
     *            'JUMBPOT=12 R0 IS RESET TO BE ',UMBR0*TOANGS
               END IF
               IF(UMBR0.LT.0.7D+00*TOBOHR.OR.
     *            UMBR0.GT.3.0D+00*TOBOHR) THEN
                  IF(MASWRK) THEN
                  WRITE(IW,'(/1X,A,A/)')
     *            'ERROR: JUMBPOT=12 R0 IS TOO SHORT OR TOO LONG. ',
     *            'THERE IS UNLIKELY A SADDLE POINT.'
                  END IF
                  RETURN
               END IF
            END IF
            IF(ENBIAS.GT.1.6D-06) CVGED = .FALSE.
         END IF
C
         IF(CVGED) THEN
            IF(MOD(ISTEP,KOUT).EQ.0.OR.ISTEP.EQ.NSTEP) THEN
               GOTO 500
            ELSE
               GOTO 499  ! MAKE IT UP
            END IF
         END IF
C
         IF(NSTEP.GT.0) THEN
         CALL DCOPY(NCVAL,EG,1,OLDF,1)
         CALL SYMDR(CRDNEW)
         DO 201 IAT =1,NATTOT
            IF((EG(1,IAT)+EG(2,IAT)+EG(3,IAT)).EQ.ZERO) GOTO 201
            IF(MMHESS.EQ.1.AND.ISTEP.LE.3) THEN
               C(1,IAT)=C(1,IAT) - EG(1,IAT)*0.50D+00
               C(2,IAT)=C(2,IAT) - EG(2,IAT)*0.50D+00
               C(3,IAT)=C(3,IAT) - EG(3,IAT)*0.50D+00
            ELSE
               C(1,IAT)=C(1,IAT) + CRDNEW(IAT*3-2)
               C(2,IAT)=C(2,IAT) + CRDNEW(IAT*3-1)
               C(3,IAT)=C(3,IAT) + CRDNEW(IAT*3  )
            END IF
  201    CONTINUE
C
C        - COPY SOME C(3,*) BACK TO BE CORD(3,*) -
         IAT = NATSV
         IF(NACTMM.EQ.0) THEN
            DO 200 IFFAT =1,NFFAT
               IF(LISTQM(IFFAT).EQ.0.OR.LISTQM(IFFAT).GT.NATSV)THEN
                  IAT=IAT+1
                  CORD(1,IFFAT) = C(1,IAT)
                  CORD(2,IFFAT) = C(2,IAT)
                  CORD(3,IFFAT) = C(3,IAT)
               ELSE IF(LISTQM(IFFAT).GT.0)THEN
                  CORD(1,IFFAT) = C(1,LISTQM(IFFAT))
                  CORD(2,IFFAT) = C(2,LISTQM(IFFAT))
                  CORD(3,IFFAT) = C(3,LISTQM(IFFAT))
               END IF
  200       CONTINUE
         ELSE
            DO 202 I=1,NACTMM
               IFFAT = LACTMM(I)
               IF(LISTQM(IFFAT).EQ.0.OR.LISTQM(IFFAT).GT.NATSV)THEN
                  IAT=IAT+1
                  CORD(1,IFFAT) = C(1,IAT)
                  CORD(2,IFFAT) = C(2,IAT)
                  CORD(3,IFFAT) = C(3,IAT)
               ELSE IF(LISTQM(IFFAT).GT.0)THEN
                  CORD(1,IFFAT) = C(1,LISTQM(IFFAT))
                  CORD(2,IFFAT) = C(2,LISTQM(IFFAT))
                  CORD(3,IFFAT) = C(3,LISTQM(IFFAT))
               END IF
  202       CONTINUE
         END IF
         END IF
C
         CALL DAWRIT(IDAF,IODA,C,NATSV*3,1,0)
C
C        -- SYNCHRONIZE CORD EVERY STEP --
         IF(LISTQM(NFFAT+NATSV+1).GT.0)THEN
            DO IAT = 1, NATSV
               KFFAT = LISTQM(NFFAT+IAT)
               IF(KFFAT.GT.0) THEN
                  C(1,IAT) = CORD(1,KFFAT)
                  C(2,IAT) = CORD(2,KFFAT)
                  C(3,IAT) = CORD(3,KFFAT)
               END IF
            ENDDO
         END IF
         IF(GOPARR) THEN
            CALL DDI_BCAST(461,'F',C,3*NATSV,MASTER)
            CALL DDI_BCAST(462,'F',CORD,3*NFFAT,MASTER)
         END IF
C
C        -- ADJUST QMCX,Y,Z ON THE FLY --
C           ALSO UPDATE LQMCT
C           (DO THIS ONLY WHEN PBC AND SWRB2 ARE USED)
         IF(SWRB2.LT.1.0D+08) THEN
            XMAX = -1.0D+30
            YMAX = -1.0D+30
            ZMAX = -1.0D+30
            XMIN =  1.0D+30
            YMIN =  1.0D+30
            ZMIN =  1.0D+30
            DO IAT = 1,NATSV
               XMAX = MAX(XMAX,C(1,IAT))
               YMAX = MAX(YMAX,C(2,IAT))
               ZMAX = MAX(ZMAX,C(3,IAT))
               XMIN = MIN(XMIN,C(1,IAT))
               YMIN = MIN(YMIN,C(2,IAT))
               ZMIN = MIN(ZMIN,C(3,IAT))
            ENDDO
            QMSIZE = ZERO
            QMSIZE = MAX(QMSIZE,XMAX-XMIN)
            QMSIZE = MAX(QMSIZE,YMAX-YMIN)
            QMSIZE = MAX(QMSIZE,ZMAX-ZMIN)
            QMSIZE = QMSIZE*1.732D+00
            QMCX   = (XMAX+XMIN)*PT5
            QMCY   = (YMAX+YMIN)*PT5
            QMCZ   = (ZMAX+ZMIN)*PT5
            R2NEAR = 100.0D+00
            NEAR   = 0
            DO IAT = 1,NATSV
               XI = C(1,IAT) - QMCX
               YI = C(2,IAT) - QMCY
               ZI = C(3,IAT) - QMCZ
               R2 = XI**2 + YI**2 + ZI**2
               IF(R2.LT.R2NEAR) THEN
                  R2NEAR = R2
                  NEAR   = IAT
               END IF
            ENDDO
            LQMCT = NEAR
            IF(MASWRK) WRITE(IW,'(A,3(F16.10,1X),A,I4,A,I7)')
     *      ' QM CENTER = ',QMCX*TOANGS,QMCY*TOANGS,QMCZ*TOANGS,
     *      'NEAR ATOM',LQMCT,' AT STEP ',ISTEP+1
         END IF
      ENDDO
C
C     -- NOT LOCATED --
C
      IF(MASWRK) THEN
         WRITE(IW,*)' '
         WRITE(IW,*)'QUANPOL OPTIMIZATION FAILED. ',
     *              'PLEASE RESTART THE JOB USING THE BEST GEOMETRY.'
         WRITE(IW,*)' '
      END IF
      CALL FLSHBF(IW)
      CALL FLSHBF(36)
      CALL FLSHBF(39)
      CALL FLSHBF(40)
      CALL SEQCLO(36,'KEEP')
      CALL SEQCLO(39,'KEEP')
      CALL SEQCLO(40,'KEEP')
      RETURN
C
C     -- LOCATED --
C
  500 CONTINUE
      IF(MASWRK) THEN
         WRITE(IW,*)' '
         WRITE(IW,*)'================= QUANPOL OPTIMIZATION',
     *              ' SUCCESSFULLY COMPLETED ================='
         WRITE(IW,*)' '
      END IF
      CALL FLSHBF(IW)
      CALL FLSHBF(36)
      CALL FLSHBF(39)
      CALL FLSHBF(40)
      CALL SEQCLO(36,'KEEP')
      CALL SEQCLO(39,'KEEP')
      CALL SEQCLO(40,'KEEP')
      RETURN
C
  999 FORMAT(1X,A8,3X,F5.1,1X,F19.13,1X,F19.13,1X,F19.13)
 1000 FORMAT(1X,A10,1X,F5.1,1X,F19.13,1X,F19.13,1X,F19.13)
C
      END
C*MODULE QUANPOC  *DECK HESSMM
!>
!> @brief    pure MM Hessian
!>
!> @author   Hui Li
!>           - Dec 2014
!>
!> @details  double displacement
!>
      SUBROUTINE HESSMM(ATMNAM,CORD,CORDSV,ZANF,
     *                  ZMAS,CHARG,POL,POLSV,DIP,
     *                  FIELD1,FIELD2,FIELD3,
     *                  SIG,EPS,SIG2,EPS2,
     *                  BOND0,FCBOND,FCSTBD,
     *                  ANGL0,FCANGL,FCWAGG,
     *                  DIHB0,FCDIHB,FCDIHR,
     *                  VROT,NNN,GAMA,IPAIR,
     *                  KLIST,LLIST,KBLST,MLIST,NLIST,
     *                  L1213J,L14J,
     *                  FFGRD2,FCLJTP,NTYPE,
     *                  XTS,YTS,ZTS,CMAT1,
     *                  POT1,POT2,QRXN1,QRXN2,NTS,
     *                  NONLS1,NONLSTQ,MAPLST,CMAPCO,
     *                  LSTCELL,NONLS2,CORDSV2,CORDSVQ,
     *                  MVFASTS2,MVFASTS3,MVFASTS4,
     *                  MVFASTL2,MVFASTL3,MVFASTL4,
     *                  AFIX,QFIX,
     *                  RFIX,IDATOM,DAI,IDDAI,
     *                  VFIX1,VFIX2,XTSFIX,YTSFIX,ZTSFIX,
     *                  NONLSA,NONLSB,L1213A,L1213B,L14A,L14B,
     *                  CHARGB,
     *                  SIGB,EPSB,SIG2B,EPS2B,
     *                  NONLSPMA,
     *                  L1213PMA,L14PMA,CORDB,
     *                  LSBONDPMA,LSANGLPMA,
     *                  LSDIHRPMA,LSDIHBPMA,
     *                  LSWAGGPMA,LSCMAPPMA,
     *                  NONLSPMB,L1213PMB,L14PMB,
     *                  LSBONDPMB,LSANGLPMB,
     *                  LSDIHRPMB,LSDIHBPMB,
     *                  LSWAGGPMB,LSCMAPPMB,UMBHIS,UM2HIS,
     *                  VIBGRD,VIBDIP,VIBHSS,VIBDDM,
     *                  NATVIB)
      use mx_limits, only: mxatm,mxao
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      LOGICAL GOPARR,DSKWRK,MASWRK
      LOGICAL CVGED,PROJCT,GOTDDM,GOTADM,STATPT,SCFOK,
     *        PRTSCN,GOTFRQ,GOTEH,EOF,RSTART
C
      PARAMETER (MXIRR=14)
      PARAMETER (ZERO=0.0D+00)
      PARAMETER (THREE=3.0D+00)
      PARAMETER (DB2AU=1.0D+00/2.541766D+00)
C
      CHARACTER*10  ATMNAM
C
      LOGICAL EFLDL
C
      DIMENSION ATMNAM(NFFAT),CORD(3,NFFAT),CORDSV(3,NFFAT),ZANF(NFFAT),
     *          ZMAS(NFFAT),
     *          CHARG(NFFAT),POL(NFFAT),DIP(3,NFFAT),
     *          FIELD1(3,NFFAT),FIELD2(3,NFFAT),FIELD3(3,NFFAT),
     *          SIG(NFFAT),EPS(NFFAT),SIG2(NFFAT),EPS2(NFFAT),
     *          BOND0(NBOND),FCBOND(NBOND),
     *          ANGL0(NANGL),FCANGL(NANGL),
     *          FCWAGG(NWAGG),DIHB0(NDIHB),FCDIHB(NDIHB),FCDIHR(3,*),
     *          VROT(NDIHR),NNN(NDIHR),GAMA(NDIHR),IPAIR(2,NBOND),
     *          KLIST(3,NANGL),
     *          LLIST(4,NDIHR),MLIST(4,NWAGG),
     *          L1213J(2,*),L14J(2,NDIHR),
     *          FFGRD2(3,NFFAT),NTYPE(*),
     *          XTS(NTS),YTS(NTS),ZTS(NTS),CMAT1(NTS,NTS),
     *          POT1(NTS),POT2(NTS),QRXN1(NTS),QRXN2(NTS),NONLS1(2,*),
     *          NONLSTQ(*),MAPLST(6,*),CMAPCO(4,4,24,24,3)
      DIMENSION ENALL(100)
      DIMENSION VIBGRD(3*NATVIB,3,NATVIB,2),
     *          VIBDIP(3,3,NATVIB,2),
     *          VIBHSS(3*NATVIB,3*NATVIB),
     *          VIBDDM(3,3*NATVIB)
      DIMENSION DIPOLE(3)
C
      COMMON /EFLDC / EVEC(3),EFLDL
      COMMON /FFENGY/ EN12,EN123,EN123R4,EN123B4,EN234W1,ENCHAR,ENLJR,
     *                ENLJD,ENPOL,XENPOL,ENRXN,XENRXN,ENRXNPOL,ENRXNR,
     *                EN12312,ENQUANP(30),
     *                ENBIAS,ENCENT,ENUCCH,ENCMAP,ENPOT,ENKIN,ENTOT
      COMMON /FFFIXS/ ENFIXSO,FIXEPS,FIXTOL,FIXA,FIXQ,RALLMM,RALLQM,
     *                RADMM(200),RADQM(200),NRADMM,NRADQM,IFIXSOL,
     *                LFFDAI,LFFDAIT,LFFIDDAI,LFFIDTMP,LFFTMPTS,
     *                LFFAFIX,LFFIDATOM,LFFRFIX,LFFQFIX,NTSATM,
     *                LFFQFIXMP,LFFQFIXTA,LFFQFIXXY,
     *                LFFXTSFIX,LFFYTSFIX,LFFZTSFIX,
     *                LFFVFIX1,LFFVFIX2,NCYCLE,MXFFTS,NFFTS
      COMMON /FFFREE/ SOL1CH,SOL1LJ,SOL1IM,SOLFRE1MM,
     *                SOL2CH,SOL2LJ,SOL2IM,SOLFRE2MM,
     *                ASOL1CH,ASOL1LJ,ASOL1IM,ASOL1MM,
     *                ASOL2CH,ASOL2LJ,ASOL2IM,ASOL2MM,
     *                PMF1BD,PMF1AG,PMF1DR,PMF1DB,PMF1WG,PMF1CM,
     *                PMF1CH,PMF1PO,PMF1LJ,PMF1IM,PMFFRE1MM,
     *                EN12A,EN12B,EN123A,EN123B,EN123R4A,EN123R4B,
     *                EN234W1A,EN234W1B,EN123B4A,EN123B4B,
     *                EN12312A,EN12312B,
     *                ENCHARA,ENCHARB,ENPOLA,ENPOLB,
     *                ENRXNA,ENRXNB,ENRXNPOLA,ENRXNPOLB,
     *                ENRXNRA,ENRXNRB,ENLJRA,ENLJRB,ENLJDA,ENLJDB,
     *                ENCMAPA,ENCMAPB,ENFIXSOA,ENFIXSOB,
     *                ENCENTA,ENCENTB,WSIMUL,WPERT1,WPERT2,
     *                IFEPTYP
      COMMON /FFMDPA/ DT,DT2,TEMP0,PRES0,POLTOL,VIR(3),PMEAN,VOLAV,
     *                ENPAV,ENKAV,TEMPAV,BERENDT,BERENDP,VELMAX,
     *                PMEANX,PMEANY,PMEANZ,NSTEP,KMASTER,KOUTACT(2),
     *                IHESS,INTALG,ITSTAT,IPSTAT,JOUT,KOUT,LOUT
      COMMON /FFMPT3/ NACTMM,LACTMM(2020),LFFDIPOLD,JUMBUP,
     *                NACTQM,LACTQM(2020),LFFOLDC,LFFQMVELSV,MMHESS,
     *                LFFQMCHG,LFFQMCHGB,ISWAP,R2SWAP,DFTBMM
      COMMON /FFPARM/ NFFAT,NBOND,NANGL,NDIHR,NDIHB,NCMAP,NWAGG,
     *                N1213J,N14J,NLKQMM,IDOCHG,IDOPOL,IDOLJ,IDOCMAP
      COMMON /FFRMSD/ DIPT(8),TIMGYRA,TIMRALL,
     *                NATPDB,NGYRA,NDIEL,NRALL,
     *                LFFRALL0,NRMSD,LFFRMSD0,KFREEAB(201),
     *                NRIJMM,IJRMM(2,100),NRIJQM,IJRQM(2,100),
     *                NAIJKMM,IJKMM(3,100),NAIJKQM,IJKQM(3,100),
     *                NFIXMM,IFIXMM(200),NFIXQM,IFIXQM(200)
      COMMON /FMCOM / X(1)
      COMMON /FUNCT / E,EG(3,MXATM)
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON /MASSES/ ZMASS(MXATM)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /SYMBLK/ NIRRED,NSALC,NSALC2,NSALC3,NSAFMO
      COMMON /SYMREP/ IRPNAM(MXIRR),IPA(MXIRR),LAMBDA(MXIRR),
     *                LAMBD0(MXIRR),IADDR1(MXIRR),IADDR2(MXIRR),
     *                IADDR3(MXIRR)
      COMMON /THERMD/ FREQ(3*MXATM),TEMP(10),SCLFAC,NTEMP,PRTSCN,GOTFRQ
      COMMON /RUNLAB/ TITLE(10),ANAM(MXATM),BNAM(MXATM),BFLAB(MXAO)
C
      DATA BLANK/8H        /
      DATA IRPA /4HA   /
C
C     HUI LI, DEC 25, 2014
C
      IF(NSTEP.LT.0) RETURN
C
C     -- DEFINE INFOA VARIABLES --
      IF(NACTMM.EQ.0) THEN
         NCOORD = 3*NFFAT
         CALL DCOPY(3*NFFAT,CORD,1,C,1)
         CALL DCOPY(NFFAT,ZANF,1,ZAN,1)
         DO IFFAT =1,NFFAT
            READ(UNIT=ATMNAM(IFFAT)(1:8),FMT='(A8)') ANAM(IFFAT)
            BNAM(IFFAT)=BLANK
            ZMASS(IFFAT)=ZMAS(IFFAT)/1822.88850204D+00
         ENDDO
      ELSE
         NCOORD = 3*NACTMM
         DO I=1,NACTMM
            IFFAT  = LACTMM(I)
            C(1,I) = CORD(1,IFFAT)
            C(2,I) = CORD(2,IFFAT)
            C(3,I) = CORD(3,IFFAT)
            ZAN(I) = ZANF(IFFAT)
            READ(UNIT=ATMNAM(IFFAT)(1:8),FMT='(A8)') ANAM(I)
            BNAM(I)=BLANK
            ZMASS(I)=ZMAS(IFFAT)/1822.88850204D+00
         ENDDO
      END IF
C
C     -- PRINT ATOMIC COORDINATES FOR MACMOLPLT --
C
      IF (MASWRK) THEN
         WRITE (IW,50)
         NAT = NFFAT
         IF(NACTMM.GT.0) NAT = NACTMM
         DO IAT = 1,NAT
            WRITE (IW,60) ANAM(IAT),BNAM(IAT),ZAN(IAT),
     *                      C(1,IAT),C(2,IAT),C(3,IAT)
         ENDDO
         NAT = 0
      END IF
   50 FORMAT(/1X,'ATOM',6X,'ATOMIC',22X,'COORDINATES (BOHR)'/
     *         11X,'CHARGE',9X,'X',19X,'Y',19X,'Z')
   60 FORMAT(1X,A8,A2,F5.1,F17.10,2F20.10)
C
         CVGED  = .FALSE.
C
C        -- CALCULATE ENERGY AND GRADIENT
         CALL NONBOND(0,CORD,CORDSV,CORDSV2,CORDSVQ,
     *                NONLS1,NONLS2,
     *                NONLSTQ,LSTCELL,
     *                MVFASTS2,MVFASTS3,MVFASTS4,
     *                MVFASTL2,MVFASTL3,MVFASTL4,
     *                NONLSA,NONLSB,
     *                NONLSPMA,NONLSPMB)
         CALL VCLR(FFGRD2,1,3*NFFAT)
         CALL E00012(CORD,FFGRD2,BOND0,FCBOND,IPAIR,CORDB,
     *               LSBONDPMA,LSBONDPMB)
         CALL E00123(CORD,FFGRD2,ANGL0,FCANGL,KLIST,CORDB,
     *               LSANGLPMA,LSANGLPMB)
         CALL E12312(CORD,FFGRD2,ANGL0,KLIST,BOND0,FCSTBD,
     *               KBLST,CORDB,LSANGLPMA,LSANGLPMB)
         CALL E123B4(CORD,FFGRD2,DIHB0,FCDIHB,NLIST,CORDB,
     *               LSDIHBPMA,LSDIHBPMB)
         CALL E234W1(CORD,FFGRD2,      FCWAGG,MLIST,CORDB,
     *               LSWAGGPMA,LSWAGGPMB)
         CALL E123R4(CORD,FFGRD2,VROT,GAMA,NNN,LLIST,CORDB,
     *               LSDIHRPMA,LSDIHRPMB,FCDIHR)
         CALL ECMAP (CORD,FFGRD2,MAPLST,CMAPCO,CORDB,
     *               LSCMAPPMA,LSCMAPPMB)
         CALL ELJ126(CORD,FFGRD2,SIG,EPS,SIG2,EPS2,L14J,NONLS1,
     *               L1213J,SIGB,EPSB,SIG2B,EPS2B,
     *               NONLSA,NONLSB,L1213A,L1213B,L14A,L14B,
     *               NONLSPMA,L1213PMA,L14PMA,CORDB,
     *               NONLSPMB,L1213PMB,L14PMB,FCLJTP,NTYPE)
         CALL ESPHER(CORD,FFGRD2)
         CALL ECHARG(CORD,FFGRD2,CHARG,NONLS1,L1213J,L14J,
     *               CHARGB,NONLSA,NONLSB,L1213A,L1213B,L14A,L14B,
     *               NONLSPMA,L1213PMA,L14PMA,CORDB,
     *               NONLSPMB,L1213PMB,L14PMB)
         CALL UMBRELLA(CORD,FFGRD2,UMBHIS,UM2HIS)
         IF(IFIXSOL.EQ.0) THEN
            CALL CHGRXN(CORD,FFGRD2,CHARG,XTS,YTS,ZTS,
     *                  CMAT1,POT1,QRXN1,NTS)
            CALL POLRXN(CORD,FFGRD2,CHARG,POL,POLSV,DIP,
     *                  FIELD1,FIELD2,FIELD3,
     *                  XTS,YTS,ZTS,CMAT1,POT1,POT2,QRXN1,QRXN2,NTS,
     *                  NONLS1,L1213J)
         END IF
         IF(IFIXSOL.EQ.1) THEN
            NATSV = NAT   ! NEED THIS BECAUSE FIXSOL USES NFFAT+NAT
            NAT   = 0
            CALL FIXSOL(CORD,FFGRD2,CHARG,ZANF,AFIX,QFIX,
     *                  VFIX1,VFIX2,XTSFIX,YTSFIX,ZTSFIX,
     *                  RFIX,IDATOM,DAI,IDDAI,
     *                  POL,POLSV,DIP,FIELD1,FIELD2,FIELD3,
     *                  NONLS1,L1213J)
            NAT   = NATSV
         END IF
         IF(GOPARR) THEN
            CALL VCLR(ENALL,1,100)
            ENALL( 1) = EN12
            ENALL( 2) = EN123
            ENALL( 3) = EN123R4
            ENALL( 4) = EN123B4
            ENALL( 5) = EN234W1
            ENALL( 6) = ENCHAR
            ENALL( 7) = ENLJR
            ENALL( 8) = ENLJD
C           ENALL(12) = ENUCCH
            ENALL(13) = ENRXN
            ENALL(14) = ENRXNR
C           ENALL(15) = ENCENT
            ENALL(16) = ENCMAP
            ENALL(18) = SOL1CH
            ENALL(19) = SOL1LJ
            ENALL(21) = SOL2CH
            ENALL(22) = SOL2LJ
            ENALL(33) = PMF1BD
            ENALL(34) = PMF1AG
            ENALL(35) = PMF1DR
            ENALL(36) = PMF1DB
            ENALL(37) = PMF1WG
            ENALL(38) = PMF1CM
            ENALL(39) = PMF1CH
            ENALL(40) = PMF1LJ
            ENALL(41) = ENBIAS
            ENALL(42) = EN12312
            CALL DDI_GSUMF(2410,ENALL  ,42)
            CALL DDI_GSUMF(2411,FFGRD2,3*NFFAT)
            EN12      = ENALL( 1)
            EN123     = ENALL( 2)
            EN123R4   = ENALL( 3)
            EN123B4   = ENALL( 4)
            EN234W1   = ENALL( 5)
            ENCHAR    = ENALL( 6)
            ENLJR     = ENALL( 7)
            ENLJD     = ENALL( 8)
C           ENUCCH    = ENALL(12)
            ENRXN     = ENALL(13)
            ENRXNR    = ENALL(14)
C           ENCENT    = ENALL(15)
            ENCMAP    = ENALL(16)
            SOL1CH    = ENALL(18)
            SOL1LJ    = ENALL(19)
            SOL2CH    = ENALL(21)
            SOL2LJ    = ENALL(22)
            PMF1BD    = ENALL(33)
            PMF1AG    = ENALL(34)
            PMF1DR    = ENALL(35)
            PMF1DB    = ENALL(36)
            PMF1WG    = ENALL(37)
            PMF1CM    = ENALL(38)
            PMF1CH    = ENALL(39)
            PMF1LJ    = ENALL(40)
            ENBIAS    = ENALL(41)
            EN12312   = ENALL(42)
         END IF
         DMMX= ZERO
         DMMY= ZERO
         DMMZ= ZERO
         AMX = ZERO
         AMY = ZERO
         AMZ = ZERO
         AMT = ZERO
         DO JFFAT = 1, NFFAT
            AMX = AMX  + CORD(1,JFFAT)*ZMAS(JFFAT)
            AMY = AMY  + CORD(2,JFFAT)*ZMAS(JFFAT)
            AMZ = AMZ  + CORD(3,JFFAT)*ZMAS(JFFAT)
            AMT = AMT  + ZMAS(JFFAT)
         ENDDO
         AMX = AMX/AMT
         AMY = AMY/AMT
         AMZ = AMZ/AMT
         DO JFFAT = 1, NFFAT
            DMMX=DMMX + CHARG(JFFAT)*(CORD(1,JFFAT)-AMX)
            DMMY=DMMY + CHARG(JFFAT)*(CORD(2,JFFAT)-AMY)
            DMMZ=DMMZ + CHARG(JFFAT)*(CORD(3,JFFAT)-AMZ)
         ENDDO
         IF(IDOPOL.GT.0) THEN
            DO JFFAT = 1, NFFAT
               DMMX=DMMX + DIP(1,JFFAT)
               DMMY=DMMY + DIP(2,JFFAT)
               DMMZ=DMMZ + DIP(3,JFFAT)
            ENDDO
         END IF
         DIPOLE(1) = DMMX/DB2AU
         DIPOLE(2) = DMMY/DB2AU
         DIPOLE(3) = DMMZ/DB2AU
C
C        - ZERO OFF SOME FORCES -
C
         IF(NACTMM.GT.0) THEN
            DO KOPT = 1, NACTMM
               IFFAT = LACTMM(KOPT)
               IF(IFFAT.GT.0)FFGRD2(1,IFFAT)=FFGRD2(1,IFFAT)+1.0D+03
            ENDDO
            DO IFFAT = 1, NFFAT
               IF(FFGRD2(1,IFFAT).GT.0.5D+03) THEN
                  FFGRD2(1,IFFAT) = FFGRD2(1,IFFAT)-1.0D+03
               ELSE
                  DO III = 1, 3
                     FFGRD2(III,IFFAT) = ZERO
                  ENDDO
               END IF
            ENDDO
         END IF
         DO KFIX=1,NFIXMM
            DO III = 1, 3
               FFGRD2(III,IFIXMM(KFIX)) = ZERO
            ENDDO
         ENDDO
C
C        - TEST CONVERGENCE
C
         CONVF  = 1.0D-05
         GRDMAX = ZERO
         GRDRMS = ZERO
         DO IFFAT=1, NFFAT
            GRDMAX=MAX(GRDMAX,ABS(FFGRD2(1,IFFAT)))
            GRDMAX=MAX(GRDMAX,ABS(FFGRD2(2,IFFAT)))
            GRDMAX=MAX(GRDMAX,ABS(FFGRD2(3,IFFAT)))
            GRDRMS=GRDRMS + FFGRD2(1,IFFAT)**2
     *                    + FFGRD2(2,IFFAT)**2
     *                    + FFGRD2(3,IFFAT)**2
         ENDDO
         IF(NACTMM.EQ.0) GRDRMS = SQRT(GRDRMS/(3*(NFFAT-NFIXMM)))
         IF(NACTMM.GT.0) GRDRMS = SQRT(GRDRMS/(3*(NACTMM-NFIXMM)))
         IF(GRDMAX.LT.CONVF.AND.GRDRMS.LT.(CONVF/THREE)) CVGED=.TRUE.
C
C        -- CALCULATE PROPERTIES
         NAT = 0
         CALL OPTPROP(0,1)
         NAT = NFFAT
         IF(NACTMM.GT.0) NAT = NACTMM
         E= ENPOT
C
         IF(NACTMM.EQ.0) THEN
            CALL DCOPY(3*NFFAT,CORD,1,C,1)
            CALL DCOPY(3*NFFAT,FFGRD2,1,EG,1)
         END IF
         IF(NACTMM.GT.0) THEN
            DO I=1,NACTMM
               IFFAT = LACTMM(I)
               C(1,I)  = CORD(1,IFFAT)
               C(2,I)  = CORD(2,IFFAT)
               C(3,I)  = CORD(3,IFFAT)
               EG(1,I) = FFGRD2(1,IFFAT)
               EG(2,I) = FFGRD2(2,IFFAT)
               EG(3,I) = FFGRD2(3,IFFAT)
            ENDDO
         END IF
C        -- MUST SAVE THE VIB 0 GRADIENT
         CALL DAWRIT(IDAF,IODA,EG,NCOORD, 3,0)
C
      CALL TIMIT(1)
C
      IF(MASWRK) WRITE(IW,'(//1X,A/1X,A/1X,A//)')
     *'***********************************',
     *'QUANPOL PURE MM HESSIAN CALCULATION',
     *'***********************************'
C
      CALL VCLR(VIBHSS,3*NAT*3*NAT,1)
      CALL VCLR(VIBDDM,3*3*NAT,1)
C
C     -- TRY TO READ IN $HESS --
      CALL FCMIN(VIBHSS,NCOORD,GOTEH)
      IF(GOTEH) GOTO 700
      IF(.NOT.GOTEH) THEN
         IF(MASWRK) WRITE(IW,'(/1X,A/)')
     *   'NO $HESS IS FOUND IN INPUT FILE. MUST COMPUTE IT.'
      END IF
C
C     -- TRY TO READ IN $VIB FROM THE INPUT FILE
C
      CALL SEQREW(IR)
      CALL FNDGRP(IR,' $VIB   ',IEOF)
      IF(IEOF.EQ.1) THEN
         RSTART = .FALSE.
      ELSE
         RSTART = .TRUE.
         IF(MASWRK) WRITE(IW,9120)
      END IF
 9120 FORMAT(/5X,'A $VIB GROUP WAS FOUND IN YOUR INPUT.'/
     *        5X,'THEREFORE, THIS IS A NUMERICAL HESSIAN RESTART.')
C
      IFCM=35
      IF (MASWRK) CALL SEQOPN(IFCM,'RESTART','NEW',.FALSE.,'FORMATTED')
      IF (MASWRK) WRITE(IFCM,8010) TITLE
 8010 FORMAT('ENERGY/GRADIENT/DIPOLE RESTART DATA FOR RUNTYP=HESSIAN'/
     *         10A8)
      IF (MASWRK) WRITE(IW,9130) BLANK,0,0,0,E
C     CALL PUVIB(IFCM,IW,.FALSE.,NCOORD,0,0,0,
C    *           E,EG,DIPOLE)
      IF(RSTART) THEN
         IF(MASWRK) WRITE(IW,9050)
 9050    FORMAT(/5X,'READING $VIB RESTART DATA FROM $VIB CARDS...')
         ECOMPUT = E
         IREAD   = 0
         IVIB0   = 0
         IAT0    = 0
         ICOORD0 = 0
         CALL RDVIB(IR,IW,IREAD,NCOORD,IVIB0,IAT0,ICOORD0,
     *              E,EG,DIPOLE,EOF,DUMY,DUMY)
         IF(EOF) THEN
            IF(MASWRK) WRITE(IW,*) 'PROBLEM READING $VIB GROUP.'
         END IF
         IF(ABS(ECOMPUT-E).GT.1.0D-05) THEN
            IF(MASWRK) WRITE(IW,'(/1X,A,A/)')
     *      'ERROR: $VIB IVIB=0 E= DOES NOT MATCH THE COMPUTED ',
     *      'ENERGY. QUANPOL WILL QUIT.'
            CALL ABRT
         END IF
      END IF
C
      DISPL=0.01D+00
      DO IVIB = 1, 2
      IF(IVIB.EQ.1) SIGN=+1.0D+00
      IF(IVIB.EQ.2) SIGN=-1.0D+00
      DO IAT = 1, NAT
      IFFAT = IAT
      IF(NACTMM.GT.0) IFFAT = LACTMM(IAT)
      DO III = 1, 3
C
         IF(RSTART) THEN
            KVIB    = IVIB
            KAT     = IAT
            KCOORD  = III
            CALL RDVIB(IR,IW,IREAD,NCOORD,KVIB,KAT,KCOORD,
     *                 E,EG,DIPOLE,EOF,DUMY,DUMY)
            IF(EOF) RSTART = .FALSE.
         END IF
         IF(RSTART) GOTO 630
C
         CSAVE=C(III,IAT)
         C(III,IAT)=C(III,IAT)+SIGN*DISPL
         CORD(III,IFFAT)=C(III,IAT)
C
C        -- CALCULATE ENERGY AND GRADIENT
C           NO NEED TO UPDATE NONBOND LIST
         CALL VCLR(FFGRD2,1,3*NFFAT)
         CALL E00012(CORD,FFGRD2,BOND0,FCBOND,IPAIR,CORDB,
     *               LSBONDPMA,LSBONDPMB)
         CALL E00123(CORD,FFGRD2,ANGL0,FCANGL,KLIST,CORDB,
     *               LSANGLPMA,LSANGLPMB)
         CALL E12312(CORD,FFGRD2,ANGL0,KLIST,BOND0,FCSTBD,
     *               KBLST,CORDB,LSANGLPMA,LSANGLPMB)
         CALL E123B4(CORD,FFGRD2,DIHB0,FCDIHB,NLIST,CORDB,
     *               LSDIHBPMA,LSDIHBPMB)
         CALL E234W1(CORD,FFGRD2,      FCWAGG,MLIST,CORDB,
     *               LSWAGGPMA,LSWAGGPMB)
         CALL E123R4(CORD,FFGRD2,VROT,GAMA,NNN,LLIST,CORDB,
     *               LSDIHRPMA,LSDIHRPMB,FCDIHR)
         CALL ECMAP (CORD,FFGRD2,MAPLST,CMAPCO,CORDB,
     *               LSCMAPPMA,LSCMAPPMB)
         CALL ELJ126(CORD,FFGRD2,SIG,EPS,SIG2,EPS2,L14J,NONLS1,
     *               L1213J,SIGB,EPSB,SIG2B,EPS2B,
     *               NONLSA,NONLSB,L1213A,L1213B,L14A,L14B,
     *               NONLSPMA,L1213PMA,L14PMA,CORDB,
     *               NONLSPMB,L1213PMB,L14PMB,FCLJTP,NTYPE)
         CALL ESPHER(CORD,FFGRD2)
         CALL ECHARG(CORD,FFGRD2,CHARG,NONLS1,L1213J,L14J,
     *               CHARGB,NONLSA,NONLSB,L1213A,L1213B,L14A,L14B,
     *               NONLSPMA,L1213PMA,L14PMA,CORDB,
     *               NONLSPMB,L1213PMB,L14PMB)
         CALL UMBRELLA(CORD,FFGRD2,UMBHIS,UM2HIS)
         IF(IFIXSOL.EQ.0) THEN
            CALL CHGRXN(CORD,FFGRD2,CHARG,XTS,YTS,ZTS,
     *                  CMAT1,POT1,QRXN1,NTS)
            CALL POLRXN(CORD,FFGRD2,CHARG,POL,POLSV,DIP,
     *                  FIELD1,FIELD2,FIELD3,
     *                  XTS,YTS,ZTS,CMAT1,POT1,POT2,QRXN1,QRXN2,NTS,
     *                  NONLS1,L1213J)
         END IF
         IF(IFIXSOL.EQ.1) THEN
            NATSV = NAT   ! NEED THIS BECAUSE FIXSOL USES NFFAT+NAT
            NAT   = 0
            CALL FIXSOL(CORD,FFGRD2,CHARG,ZANF,AFIX,QFIX,
     *                  VFIX1,VFIX2,XTSFIX,YTSFIX,ZTSFIX,
     *                  RFIX,IDATOM,DAI,IDDAI,
     *                  POL,POLSV,DIP,FIELD1,FIELD2,FIELD3,
     *                  NONLS1,L1213J)
            NAT   = NATSV
         END IF
         IF(GOPARR) THEN
            CALL VCLR(ENALL,1,100)
            ENALL( 1) = EN12
            ENALL( 2) = EN123
            ENALL( 3) = EN123R4
            ENALL( 4) = EN123B4
            ENALL( 5) = EN234W1
            ENALL( 6) = ENCHAR
            ENALL( 7) = ENLJR
            ENALL( 8) = ENLJD
            ENALL(12) = ENUCCH
            ENALL(13) = ENRXN
            ENALL(14) = ENRXNR
C           ENALL(15) = ENCENT
            ENALL(16) = ENCMAP
            ENALL(18) = SOL1CH
            ENALL(19) = SOL1LJ
            ENALL(21) = SOL2CH
            ENALL(22) = SOL2LJ
            ENALL(33) = PMF1BD
            ENALL(34) = PMF1AG
            ENALL(35) = PMF1DR
            ENALL(36) = PMF1DB
            ENALL(37) = PMF1WG
            ENALL(38) = PMF1CM
            ENALL(39) = PMF1CH
            ENALL(40) = PMF1LJ
            ENALL(41) = ENBIAS
            ENALL(42) = EN12312
            CALL DDI_GSUMF(2410,ENALL  ,42)
            CALL DDI_GSUMF(2411,FFGRD2,3*NFFAT)
            EN12      = ENALL( 1)
            EN123     = ENALL( 2)
            EN123R4   = ENALL( 3)
            EN123B4   = ENALL( 4)
            EN234W1   = ENALL( 5)
            ENCHAR    = ENALL( 6)
            ENLJR     = ENALL( 7)
            ENLJD     = ENALL( 8)
            ENUCCH    = ENALL(12)
            ENRXN     = ENALL(13)
            ENRXNR    = ENALL(14)
C           ENCENT    = ENALL(15)
            ENCMAP    = ENALL(16)
            SOL1CH    = ENALL(18)
            SOL1LJ    = ENALL(19)
            SOL2CH    = ENALL(21)
            SOL2LJ    = ENALL(22)
            PMF1BD    = ENALL(33)
            PMF1AG    = ENALL(34)
            PMF1DR    = ENALL(35)
            PMF1DB    = ENALL(36)
            PMF1WG    = ENALL(37)
            PMF1CM    = ENALL(38)
            PMF1CH    = ENALL(39)
            PMF1LJ    = ENALL(40)
            ENBIAS    = ENALL(41)
            EN12312   = ENALL(42)
         END IF
         ENPOT = EN12 + EN123 + EN123R4 + EN234W1 + EN123B4 + ENCHAR
     *         + ENPOL  + ENRXN  + ENRXNPOL+ENRXNR + ENLJR + ENLJD
     *         + ENCMAP + ENFIXSO+ EN12312
         E= ENPOT
         IF(NACTMM.EQ.0) THEN
            DO JAT = 1, NFFAT
               JFFAT = JAT
               DO JJJ = 1, 3
                  EG(JJJ,JAT) =FFGRD2(JJJ,JFFAT)
               ENDDO
            ENDDO
         END IF
         IF(NACTMM.GT.0) THEN
            DO JAT = 1, NACTMM
               JFFAT = LACTMM(JAT)
               DO JJJ = 1, 3
                  EG(JJJ,JAT) =FFGRD2(JJJ,JFFAT)
               ENDDO
            ENDDO
         END IF
C        -- CALCULATE DIPOLE MOMENT
         DMMX= ZERO
         DMMY= ZERO
         DMMZ= ZERO
         AMX = ZERO
         AMY = ZERO
         AMZ = ZERO
         AMT = ZERO
         DO JFFAT = 1, NFFAT
            AMX = AMX  + CORD(1,JFFAT)*ZMAS(JFFAT)
            AMY = AMY  + CORD(2,JFFAT)*ZMAS(JFFAT)
            AMZ = AMZ  + CORD(3,JFFAT)*ZMAS(JFFAT)
            AMT = AMT  + ZMAS(JFFAT)
         ENDDO
         AMX = AMX/AMT
         AMY = AMY/AMT
         AMZ = AMZ/AMT
         DO JFFAT = 1, NFFAT
            DMMX=DMMX + CHARG(JFFAT)*(CORD(1,JFFAT)-AMX)
            DMMY=DMMY + CHARG(JFFAT)*(CORD(2,JFFAT)-AMY)
            DMMZ=DMMZ + CHARG(JFFAT)*(CORD(3,JFFAT)-AMZ)
         ENDDO
         IF(IDOPOL.GT.0) THEN
            DO JFFAT = 1, NFFAT
               DMMX=DMMX + DIP(1,JFFAT)
               DMMY=DMMY + DIP(2,JFFAT)
               DMMZ=DMMZ + DIP(3,JFFAT)
            ENDDO
         END IF
         DIPOLE(1) = DMMX/DB2AU
         DIPOLE(2) = DMMY/DB2AU
         DIPOLE(3) = DMMZ/DB2AU
C
         C(III,IAT)=CSAVE
         CORD(III,IFFAT)=C(III,IAT)
C
  630    CONTINUE
         IF (MASWRK) WRITE(IW,9130) BLANK,IVIB,IAT,III,E
 9130 FORMAT(1X,A8,' IVIB=',I4,' IATOM=',I4,' ICOORD=',I4,' E=',F20.10)
C        CALL PUVIB(IFCM,IW,RSTART,NCOORD,IVIB,IAT,III,
C    *              E,EG,DIPOLE)
         CALL DCOPY(3*NAT,EG,1,VIBGRD(1,III,IAT,IVIB),1)
         VIBDIP(1,III,IAT,IVIB)=DIPOLE(1)
         VIBDIP(2,III,IAT,IVIB)=DIPOLE(2)
         VIBDIP(3,III,IAT,IVIB)=DIPOLE(3)
      ENDDO
      ENDDO
      ENDDO
C
      DO IAT = 1, NAT
         DO III = 1, 3
            DO JAT = 1, NAT
               DO JJJ = 1, 3
                  VIBHSS((JAT-1)*3+JJJ,(IAT-1)*3+III)
     *            =(VIBGRD((JAT-1)*3+JJJ,III,IAT,1)
     *             -VIBGRD((JAT-1)*3+JJJ,III,IAT,2))/(DISPL+DISPL)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      DO IAT = 1, NAT
         DO III = 1, 3
            DO JJJ = 1, 3
               VIBDDM(III,(IAT-1)*3+JJJ)
     *         =(VIBDIP(III,JJJ,IAT,1)
     *          -VIBDIP(III,JJJ,IAT,2))/(DISPL+DISPL)
            ENDDO
         ENDDO
      ENDDO
C
  700 CONTINUE
      IF(MASWRK) WRITE(IP,'(A)')
     *'----- HESSIAN MATRIX (NOT MASS WEIGHTED) -----'
      IF(MASWRK) WRITE(IP,'(A,A/A/1X,A/1X,A)')
     *'MAKE SURE TO USE MMHESS=1 WHEN MM $HESS IS USED FOR QM/MM ',
     *'GEOMETRY SEARCH',
     *'MAKE SURE TO USE MMHESS=0 WHEN QM/MM $HESS IS USED FOR QM/MM ',
     *'$QUANPO MMHESS=1 $END ! IF E(NUC) IS     0.0000000000 BELOW',
     *'$QUANPO MMHESS=0 $END ! IF E(NUC) IS NOT 0.0000000000 BELOW'
      CALL FCMPUN(VIBHSS,3*NAT)
      CALL DAWRIT(IDAF,IODA,VIBHSS,3*NAT*3*NAT,4,0)
      CALL DAWRIT(IDAF,IODA,VIBDDM,3*3*NAT,34,0)
C
C     - VIBANL -
      PROJCT =.TRUE.
      GOTDDM =.TRUE.
      GOTADM =.FALSE.
      IF((.NOT.CVGED).AND.MASWRK) WRITE(IW,9010)
 9010 FORMAT(/
     *   5X,'*******************************************************'/
     *   5X,'* THIS IS NOT A STATIONARY POINT ON THE MOLECULAR PES *'/
     *   5X,'*     THE VIBRATIONAL ANALYSIS IS NOT VALID !!!       *'/
     *   5X,'*******************************************************')
      STATPT =.TRUE.
      SCFOK  =.TRUE.
      NTEMP  = 1
      TEMP(1)= TEMP0
      SCLFAC = 1.0D+00
      NIRRED = MXIRR
      DO I =1, MXIRR
         IRPNAM(I)=IRPA
      ENDDO
      EFLDL  = .FALSE.  ! THIS HELPS PROJECT THE ROTATIONS
C
      NPART  = NAT
      NC1    = 3*NAT
      NC2    = (NC1**2+NC1)/2
      NC3    = NC1**2
      CALL VALFM(LOADFM)
      LVEC   = LOADFM + 1
      LFCM   = LVEC   + NC3
      LE     = LFCM   + NC2
      LSCR   = LE     + NC1
      LIA    = LSCR   + NC1*8
      LRM    = LIA    + NC1
      LSVT   = LRM    + NCOORD
      LSVR   = LSVT   + NC1*3
      LSVTT  = LSVR   + NC1*3
      LSVRT  = LSVTT  + NC1
      LCC    = LSVRT  + NC1
      LCOM   = LCC    + 3*NPART
      LZMS   = LCOM   + 3*NPART
      LBUF1  = LZMS   +   NPART
      LBUF2  = LBUF1  + NC3
      LDDM   = LBUF2  + NC3
      LADM   = LDDM   + NC1*3
      LAST   = LADM   + NC1*6
      NEED   = LAST   - LOADFM
      CALL GETFM(NEED)
      IF(MASWRK) WRITE(IP,'(A,A/A/1X,A/1X,A)')
     *'MAKE SURE TO USE MMHESS=1 WHEN MM $HESS IS USED FOR QM/MM ',
     *'GEOMETRY SEARCH',
     *'MAKE SURE TO USE MMHESS=0 WHEN QM/MM $HESS IS USED FOR QM/MM ',
     *'$QUANPO MMHESS=1 $END ! IF E(NUC) IS     0.0000000000 BELOW',
     *'$QUANPO MMHESS=0 $END ! IF E(NUC) IS NOT 0.0000000000 BELOW'
      CALL FGMTRX(X(LVEC),X(LFCM),X(LE),X(LSCR),X(LIA),X(LRM),
     *            X(LSVT),X(LSVR),X(LSVTT),X(LSVRT),
     *            freq,X(LCC),X(LCOM),X(LZMS),X(LBUF1),X(LBUF2),
     *            NC1,NC2,NPART,X(LDDM),GOTDDM,X(LADM),GOTADM,
     *            PROJCT,STATPT,SCFOK,0,0)
      CALL RETFM(NEED)
      IF((.NOT.CVGED).AND.MASWRK) WRITE(IW,9010)
      IF (MASWRK) WRITE(IW,9020)
      CALL TIMIT(1)
      NAT = 0
C
 9020 FORMAT(1X,'......END OF NORMAL COORDINATE ANALYSIS......')
C
      RETURN
      END
C*MODULE QUANPOC  *DECK HESSQM
!>
!> @brief    QM Hessian in QM/MM system
!>
!> @author   Rui Lai and Hui Li
!>           - Dec 2014
!>
!> @details  double displacement of selected QM atoms
!>
      SUBROUTINE HESSQM(CORD,CORDSV,DIP,
     *                  FFGRD2,QMGRD2,
     *                  LISTQM,NONLS1,NONLSTQ,
     *                  LSTCELL,NONLS2,CORDSV2,CORDSVQ,
     *                  MVFASTS2,MVFASTS3,MVFASTS4,
     *                  MVFASTL2,MVFASTL3,MVFASTL4,
     *                  NONLSA,NONLSB,
     *                  NONLSPMA,NONLSPMB,
     *                  VIBGRD,VIBDIP,VIBE,VIBHSS,VIBDDM)
      use mx_limits, only: mxatm,mxao,mxrt
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      LOGICAL GOPARR,DSKWRK,MASWRK
      LOGICAL CVGED,PROJCT,GOTDDM,GOTADM,STATPT,SCFOK,
     *        PRTSCN,GOTFRQ,GOTEH,EOF,RSTART
C
      PARAMETER (MXIRR=14)
      PARAMETER (ZERO=0.0D+00)
      PARAMETER (THREE=3.0D+00)
      PARAMETER (DB2AU=1.0D+00/2.541766D+00)
      PARAMETER (MAXML=1024)
C
      LOGICAL EFLDL
C
      DIMENSION CORD(3,NFFAT),CORDSV(3,NFFAT),
     *          DIP(3,NFFAT),
     *          FFGRD2(3,NFFAT),QMGRD2(3,NAT),
     *          LISTQM(*),NONLSTQ(*)
      DIMENSION VIBGRD(3*NAT,3,NAT,2),
     *          VIBDIP(3,3,NAT,2),VIBE(3,NAT,2),
     *          VIBHSS(3*NAT,3*NAT),
     *          VIBDDM(3,3*NAT)
      DIMENSION DIPOLE(3)
C
      COMMON /EFLDC / EVEC(3),EFLDL
      COMMON /ENRGYS/ ENUCR,EELCT,ETOT,SZ,SZZ,ECORE,ESCF,EERD,E1,E2,
     *                VEN,VEE,EPOT,EKIN,ESTATE(MXRT),STATN,EDFT(3),EDISP
      COMMON /FFDIMR/ IDIMER,IBREAK(81),N1213JMM,NMOLE,MATOM(MAXML),
     *                MCHARG(MAXML),MMULT(MAXML),MELEC(MAXML),RDIMER
      COMMON /FFFREE/ SOL1CH,SOL1LJ,SOL1IM,SOLFRE1MM,
     *                SOL2CH,SOL2LJ,SOL2IM,SOLFRE2MM,
     *                ASOL1CH,ASOL1LJ,ASOL1IM,ASOL1MM,
     *                ASOL2CH,ASOL2LJ,ASOL2IM,ASOL2MM,
     *                PMF1BD,PMF1AG,PMF1DR,PMF1DB,PMF1WG,PMF1CM,
     *                PMF1CH,PMF1PO,PMF1LJ,PMF1IM,PMFFRE1MM,
     *                EN12A,EN12B,EN123A,EN123B,EN123R4A,EN123R4B,
     *                EN234W1A,EN234W1B,EN123B4A,EN123B4B,
     *                EN12312A,EN12312B,
     *                ENCHARA,ENCHARB,ENPOLA,ENPOLB,
     *                ENRXNA,ENRXNB,ENRXNPOLA,ENRXNPOLB,
     *                ENRXNRA,ENRXNRB,ENLJRA,ENLJRB,ENLJDA,ENLJDB,
     *                ENCMAPA,ENCMAPB,ENFIXSOA,ENFIXSOB,
     *                ENCENTA,ENCENTB,WSIMUL,WPERT1,WPERT2,
     *                IFEPTYP
      COMMON /FFPARM/ NFFAT,NBOND,NANGL,NDIHR,NDIHB,NCMAP,NWAGG,
     *                N1213J,N14J,NLKQMM,IDOCHG,IDOPOL,IDOLJ,IDOCMAP
      COMMON /FFMDPA/ DT,DT2,TEMP0,PRES0,POLTOL,VIR(3),PMEAN,VOLAV,
     *                ENPAV,ENKAV,TEMPAV,BERENDT,BERENDP,VELMAX,
     *                PMEANX,PMEANY,PMEANZ,NSTEP,KMASTER,KOUTACT(2),
     *                IHESS,INTALG,ITSTAT,IPSTAT,JOUT,KOUT,LOUT
      COMMON /FFMPT3/ NACTMM,LACTMM(2020),LFFDIPOLD,JUMBUP,
     *                NACTQM,LACTQM(2020),LFFOLDC,LFFQMVELSV,MMHESS,
     *                LFFQMCHG,LFFQMCHGB,ISWAP,R2SWAP,DFTBMM
      COMMON /FFQMPA/ ENFFQM2,ENPAV2,
     *                SCFTYP2,TDDFT2,MPLEVL2,CITYP2,
     *                ICHARG2,MULT2,IDOQM2,IREDOX,IQMPKA,IQMRXN,
     *                MATOMA,MCHARGA,MULTA,MELEA,
     *                MATOMB,MCHARGB,MULTB,MELEB,
     *                IECPX,NSHELLX,IMP,JMP,ICORSH,IGTF,
     *                LFFZANX,LFFCLPX,LFFZLPX,LFFNLPX,LFFKFRSTX,
     *                LFFKLASTX,LFFLMAXX,LFFLPSKIPX,LFFIZCOREX,
     *                LFFCX,LFFIANX,LFFEXX,LFFCSX,LFFCPX,LFFCDX,
     *                LFFCFX,LFFCGX,LFFCHX,LFFCIX,LFFKSTARTX,
     *                LFFKATOMX,LFFKTYPEX,LFFKNGX,LFFKLOCX,
     *                LFFMINX,LFFMAXX,LFFMPTYPX,LFFAN0X,
     *                LFFALPN0X,LFFAN1X,LFFALPN1X,LFFMPSKPX,
     *                LFFNOAN0X,LFFNOAN1X,LFFBPARX,LFFEXPMPX,
     *                LFFCSMPX,LFFCPMPX,LFFCDMPX,LFFCFMPX,LFFMPSKIPX,
     *                LFFNOCOSHX,LFFMPKSTAX,LFFMPKNGX,LFFMPKTYPX,
     *                LFFMPKMINX,LFFMPKMAXX,LFFMPKLOCX,LFFANAMX
      COMMON /FFRMSD/ DIPT(8),TIMGYRA,TIMRALL,
     *                NATPDB,NGYRA,NDIEL,NRALL,
     *                LFFRALL0,NRMSD,LFFRMSD0,KFREEAB(201),
     *                NRIJMM,IJRMM(2,100),NRIJQM,IJRQM(2,100),
     *                NAIJKMM,IJKMM(3,100),NAIJKQM,IJKQM(3,100),
     *                NFIXMM,IFIXMM(200),NFIXQM,IFIXQM(200)
      COMMON /FMCOM / X(1)
      COMMON /FUNCT / E,EG(3,MXATM)
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /IOFILE/ IR,IW,IP,IJK,IJKT,IDAF,NAV,IODA(950)
      COMMON /MASSES/ ZMASS(MXATM)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /RUNLAB/ TITLE(10),ANAM(MXATM),BNAM(MXATM),BFLAB(MXAO)
      COMMON /OPTGRD/ XX(3,MXATM),ES,FE(20),
     *                CONVF,FMAXT,DXMAXT,RMAX,RMIN,RLIM,
     *                EIGMAX,EIGMIN,GRDERR,FRMS,FMAX,TRMAX,TRMIN,
     *                IC(20),MSTEP,NSERCH,NPMAX,NP,IFOLOW,
     *                NNEG,IUPHSS,IEXIT,ITRUPD,IPAD,KDIAGH,NPRICO
      COMMON /SCFOPT/ CONVHF,MAXIT,MCONV,NPUNCH,NPREO(4),FSHIFT
      COMMON /SYMBLK/ NIRRED,NSALC,NSALC2,NSALC3,NSAFMO
      COMMON /SYMREP/ IRPNAM(MXIRR),IPA(MXIRR),LAMBDA(MXIRR),
     *                LAMBD0(MXIRR),IADDR1(MXIRR),IADDR2(MXIRR),
     *                IADDR3(MXIRR)
      COMMON /THERMD/ FREQ(3*MXATM),TEMP(10),SCLFAC,NTEMP,PRTSCN,GOTFRQ
      COMMON /XYZPRP/ XP,YP,ZP,
     *                DMX,DMY,DMZ,
     *                QXX,QYY,QZZ,QXY,QXZ,QYZ,
     *                QMXX,QMYY,QMZZ,QMXY,QMXZ,QMYZ,
     *                OXXX,OXXY,OXXZ,OXYY,OYYY,OYYZ,
     *                OXZZ,OYZZ,OZZZ,OXYZ,
     *                OMXXX,OMXXY,OMXXZ,OMXYY,OMYYY,
     *                OMYYZ,OMXZZ,OMYZZ,OMZZZ,OMXYZ
      COMMON /OUTPUT/ NPRINT,ITOL,ICUT,NORMF,NORMP,NOPK
C
      DATA IRPA /4HA   /
C
C     -- CALCULATE VIBRATIONAL FREQUENCIES
C        RUI LAI AND HUI LI, DEC 22, 2014, LINCOLN
C
      IF(NSTEP.LT.0) RETURN
C
      NPRINT  = -5
C
      NCOORD=NAT*3
      CALL DAWRIT(IDAF,IODA,EG,NCOORD, 3,0)
      CVGED =.FALSE.
C
C     -- TRY TO READ IN $HESS --
      GOTEH = .FALSE.
      CALL VCLR(VIBGRD,1,3*NAT*3*NAT*2)
      CALL VCLR(VIBDIP,1,3*3*NAT*2)
      CALL VCLR(VIBE  ,1,3*NAT*2)
      CALL VCLR(VIBHSS,1,3*NAT*3*NAT)
      CALL VCLR(VIBDDM,1,3*3*NAT)
      CALL FCMIN(VIBHSS,NCOORD,GOTEH)
      IF(GOTEH) THEN
         CVGED =.TRUE.
         GOTO 700
      END IF
C
C        -- VIB 0 --
C
         NPUNCH = 2
C
C        -- CALCULATE ENERGY AND GRADIENT -
         CALL NONBOND(0,CORD,CORDSV,CORDSV2,CORDSVQ,
     *                NONLS1,NONLS2,
     *                NONLSTQ,LSTCELL,
     *                MVFASTS2,MVFASTS3,MVFASTS4,
     *                MVFASTL2,MVFASTL3,MVFASTL4,
     *                NONLSA,NONLSB,
     *                NONLSPMA,NONLSPMB)
         CALL VCLR(QMGRD2,1,3*NAT  )
         CALL VCLR(FFGRD2,1,3*NFFAT)
         CALL VCLR(VIR   ,1,3)
         IF(IDIMER.EQ.0)THEN
            CALL GRADX
         ELSE
            CALL DIMERX
         END IF
C        -- ESCF WAS ZERO IF SCF NOT CONVERGED --
         IF(ABS(ESCF).LT.1.0D-12) THEN
            IF(MASWRK) WRITE(IW,'(/A/)')
     *         ' ERROR: SCF WAS NOT CONVERGED. QUANPOL MUST STOP.'
            CALL ABRT
         END IF
         CALL DCOPY(3*NAT,EG,1,QMGRD2,1)
         CALL ELMOMC
         DMXA = DMX
         DMYA = DMY
         DMZA = DMZ
         IF(IFEPTYP.GT.0.AND.MATOMB.GT.0) THEN
            CALL SAVEABPROP(1)
            CALL DCOPY(3*MATOMA,C,1,X(LFFCX),1)
            CALL SAVEFFDATA
            CALL SETQMAB(2)
            CALL SETFFDATAB(2)
            CALL VCLR(QMGRD2,1,3*NAT  )
            CALL VCLR(FFGRD2,1,3*NFFAT)
            CALL VCLR(VIR   ,1,3)
            CALL GRADX
C           -- ESCF WAS ZERO IF SCF NOT CONVERGED --
            IF(ABS(ESCF).LT.1.0D-12) THEN
               IF(MASWRK) WRITE(IW,'(/A/)')
     *            ' ERROR: SCF WAS NOT CONVERGED. QUANPOL MUST STOP.'
               CALL ABRT
            END IF
            CALL DCOPY(3*NAT,EG,1,QMGRD2,1)
            CALL ELMOMC
            DMXB = DMX
            DMYB = DMY
            DMZB = DMZ
            CALL SAVEABPROP(2)
            CALL SETQMAB(1)
            CALL SETFFDATAB(1)
            CALL MIXABPROP
            DMX = (1-WSIMUL)*DMXA + WSIMUL*DMXB
            DMY = (1-WSIMUL)*DMYA + WSIMUL*DMYB
            DMZ = (1-WSIMUL)*DMZA + WSIMUL*DMZB
         END IF
         E      = ETOT
C        -- CALCULATE DIPOLE MOMENT
         DIPOLE(1)= ZERO
         DIPOLE(2)= ZERO
         DIPOLE(3)= ZERO
         IF(IDOPOL.GT.0) THEN
            DO JFFAT = 1, NFFAT
               DIPOLE(1)=DIPOLE(1) + DIP(1,JFFAT)
               DIPOLE(2)=DIPOLE(2) + DIP(2,JFFAT)
               DIPOLE(3)=DIPOLE(3) + DIP(3,JFFAT)
            ENDDO
         END IF
         DIPOLE(1)=DIPOLE(1)/DB2AU + DMX
         DIPOLE(2)=DIPOLE(2)/DB2AU + DMY
         DIPOLE(3)=DIPOLE(3)/DB2AU + DMZ
C
C        - COMBINE QM AND MM GRADIENTS -
C
         IF(LISTQM(NFFAT+NAT+1).GT.0)THEN
            DO IAT = 1, NAT
               KFFAT = LISTQM(NFFAT+IAT)
               IF(KFFAT.GT.0) THEN
                  QMGRD2(1,IAT)   = QMGRD2(1,IAT) + FFGRD2(1,KFFAT)
                  QMGRD2(2,IAT)   = QMGRD2(2,IAT) + FFGRD2(2,KFFAT)
                  QMGRD2(3,IAT)   = QMGRD2(3,IAT) + FFGRD2(3,KFFAT)
                  EG(1,IAT)       = QMGRD2(1,IAT)
                  EG(2,IAT)       = QMGRD2(2,IAT)
                  EG(3,IAT)       = QMGRD2(3,IAT)
                  FFGRD2(1,KFFAT) = QMGRD2(1,IAT)
                  FFGRD2(2,KFFAT) = QMGRD2(2,IAT)
                  FFGRD2(3,KFFAT) = QMGRD2(3,IAT)
               END IF
            ENDDO
         END IF
C
C        - ZERO OFF SOME FORCES -
C
         IF(NACTQM.GT.0) THEN
            DO KOPT = 1, NACTQM
               IAT = LACTQM(KOPT)
               IF(IAT.GT.0)QMGRD2(1,IAT)=QMGRD2(1,IAT)+1.0D+03
            ENDDO
            DO IAT = 1, NAT
               IF(QMGRD2(1,IAT).GT.0.5D+03) THEN
                  QMGRD2(1,IAT) = QMGRD2(1,IAT)-1.0D+03
               ELSE
                  DO III = 1, 3
                     QMGRD2(III,IAT) = ZERO
                     EG(III,IAT)     = ZERO
                  ENDDO
               END IF
            ENDDO
         END IF
         DO KFIX=1,NFIXQM
            DO III = 1, 3
               QMGRD2(III,IFIXQM(KFIX)) = ZERO
               EG(III,IFIXQM(KFIX))     = ZERO
            ENDDO
         ENDDO
C
C        - TEST CONVERGENCE
C
         CONVF  = 1.0D-04
         GRDMAX = ZERO
         GRDRMS = ZERO
         DO IAT=1, NAT
            GRDMAX=MAX(GRDMAX,ABS(QMGRD2(1,IAT)))
            GRDMAX=MAX(GRDMAX,ABS(QMGRD2(2,IAT)))
            GRDMAX=MAX(GRDMAX,ABS(QMGRD2(3,IAT)))
            GRDRMS=GRDRMS + QMGRD2(1,IAT)**2
     *                    + QMGRD2(2,IAT)**2
     *                    + QMGRD2(3,IAT)**2
         ENDDO
         IF(NACTQM.EQ.0)
     *   GRDRMS = SQRT(GRDRMS/(3*(NAT-NFIXQM)))
         IF(NACTQM.GT.0)
     *   GRDRMS = SQRT(GRDRMS/(3*(NACTQM-NFIXQM)))
         IF(GRDMAX.LT.CONVF.AND.GRDRMS.LT.(CONVF/THREE)) CVGED=.TRUE.
C
C        - CALCULATE PROPERTIES -
         CALL OPTPROP(0,1)
C        -- MUST SAVE THE VIB 0 GRADIENT
         CALL DAWRIT(IDAF,IODA,EG,NCOORD, 3,0)
C
      CALL TIMIT(1)
C
      IF(MASWRK) WRITE(IW,'(//1X,A/1X,A/1X,A//)')
     *'***********************************************',
     *'QUANPOL QM/MM HESSIAN CALCULATION: QM FREQUENCY',
     *'***********************************************'
      IF(GOTEH) GOTO 700
      IF(MASWRK) WRITE(IW,'(/1X,A/)')
     *'NO $HESS IS FOUND IN INPUT FILE. MUST COMPUTE IT.'
C
C     -- TRY TO READ IN $VIB FROM THE INPUT FILE
C
      CALL SEQREW(IR)
      CALL FNDGRP(IR,' $VIB   ',IEOF)
      IF(IEOF.EQ.1) THEN
         RSTART = .FALSE.
      ELSE
         RSTART = .TRUE.
         IF(MASWRK) WRITE(IW,9120)
      END IF
 9120 FORMAT(/5X,'A $VIB GROUP WAS FOUND IN YOUR INPUT.'/
     *        5X,'THEREFORE, THIS IS A NUMERICAL HESSIAN RESTART.')
C
      IFCM=35
      IF (MASWRK) CALL SEQOPN(IFCM,'RESTART','NEW',.FALSE.,'FORMATTED')
      IF (MASWRK) WRITE(IFCM,8010) TITLE
 8010 FORMAT('ENERGY/GRADIENT/DIPOLE RESTART DATA FOR RUNTYP=HESSIAN'/
     *         10A8)
      CALL PUVIB(IFCM,IW,.FALSE.,NCOORD,0,0,0,
     *           E,EG,DIPOLE)
      IF(RSTART) THEN
         IF(MASWRK) WRITE(IW,9050)
 9050    FORMAT(/5X,'READING $VIB RESTART DATA FROM $VIB CARDS...')
         ECOMPUT = E
         IREAD   = 0
         IVIB0   = 0
         IAT0    = 0
         ICOORD0 = 0
         CALL RDVIB(IR,IW,IREAD,NCOORD,IVIB0,IAT0,ICOORD0,
     *              E,EG,DIPOLE,EOF,DUMY,DUMY)
         IF(EOF) THEN
            IF(MASWRK) WRITE(IW,*) 'PROBLEM READING $VIB GROUP.'
         END IF
         IF(ABS(ECOMPUT-E).GT.1.0D-05) THEN
            IF(MASWRK) WRITE(IW,'(/1X,A,A/)')
     *      'ERROR: $VIB IVIB=0 E= DOES NOT MATCH THE COMPUTED ',
     *      'ENERGY. QUANPOL WILL QUIT.'
            CALL ABRT
         END IF
C
         DO I = 1, 2*3*NAT
            KVIB   = 0
            KAT    = 0
            KCOORD = 0
            CALL RDVIB(IR,IW,IREAD,NCOORD,KVIB,KAT,KCOORD,
     *                 E,EG,DIPOLE,EOF,DUMY,DUMY)
            IF(EOF) GOTO 500
C           - ONLY MASWRK CHANGES KVIB, KAT, KCOORD IN RDVIB -
            IF(GOPARR) THEN
               CALL DDI_BCAST(476,'I',KVIB  ,1,MASTER)
               CALL DDI_BCAST(477,'I',KAT   ,1,MASTER)
               CALL DDI_BCAST(478,'I',KCOORD,1,MASTER)
            END IF
            IVIB    = KVIB
            IAT     = KAT
            III     = KCOORD
            VIBE(III,IAT,IVIB) = E
            CALL DCOPY(3*NAT,EG,1,VIBGRD(1,III,IAT,IVIB),1)
            VIBDIP(1,III,IAT,IVIB)=DIPOLE(1)
            VIBDIP(2,III,IAT,IVIB)=DIPOLE(2)
            VIBDIP(3,III,IAT,IVIB)=DIPOLE(3)
         ENDDO
  500    CONTINUE
      END IF
C
      DISPL= 0.01D+00
      DO IVIB = 1, 2
      IF(IVIB.EQ.1) SIGN=+1.0D+00
      IF(IVIB.EQ.2) SIGN=-1.0D+00
      DO 610 IAT = 1, NAT
      IYES = 1
      IF(NACTQM.GT.0) THEN
         IYES = 0
         DO LLL = 1, NACTQM
            IF(IAT.EQ.LACTQM(LLL)) IYES = 1
         ENDDO
      END IF
      IF(IYES.EQ.0) GOTO 610
      DO III = 1, 3
C
         SUM = ZERO
         DO J=1, 3*NAT
            SUM=SUM+ABS(VIBGRD(J,III,IAT,IVIB))
         ENDDO
         IF(SUM.GT.1.0D-12) THEN
            CALL PUVIB(IFCM,IW,.FALSE.,NCOORD,IVIB,IAT,III,
     *                 VIBE(III,IAT,IVIB),VIBGRD(1,III,IAT,IVIB),
     *                 VIBDIP(1,III,IAT,IVIB))
            GOTO 630
         END IF
C
         CSAVE=C(III,IAT)
         C(III,IAT)=C(III,IAT)+SIGN*DISPL
         IF(LISTQM(NFFAT+IAT).GT.0)THEN
            CORD(III,LISTQM(NFFAT+IAT)) = C(III,IAT)
         END IF
C
C        -- CALCULATE ENERGY AND GRADIENT
C           NO NEED TO UPDATE NONBOND LIST
         CALL VCLR(QMGRD2,1,3*NAT)
         CALL VCLR(FFGRD2,1,3*NFFAT)
         NAT    = NAT
         NPUNCH = 0
         CALL VCLR(QMGRD2,1,3*NAT  )
         CALL VCLR(FFGRD2,1,3*NFFAT)
         CALL VCLR(VIR   ,1,3)
         IF(IDIMER.EQ.0)THEN
            CALL GRADX
         ELSE
            CALL DIMERX
         END IF
C        -- ESCF WAS ZERO IF SCF NOT CONVERGED --
         IF(ABS(ESCF).LT.1.0D-12) THEN
            IF(MASWRK) WRITE(IW,'(/A/)')
     *         ' ERROR: SCF WAS NOT CONVERGED. QUANPOL MUST STOP.'
            CALL ABRT
         END IF
         CALL DCOPY(3*NAT,EG,1,QMGRD2,1)
         CALL ELMOMC
         DMXA = DMX
         DMYA = DMY
         DMZA = DMZ
         IF(IFEPTYP.GT.0.AND.MATOMB.GT.0) THEN
            CALL SAVEABPROP(1)
            CALL DCOPY(3*MATOMA,C,1,X(LFFCX),1)
            CALL SAVEFFDATA
            CALL SETQMAB(2)
            CALL SETFFDATAB(2)
            CALL VCLR(QMGRD2,1,3*NAT  )
            CALL VCLR(FFGRD2,1,3*NFFAT)
            CALL VCLR(VIR   ,1,3)
            CALL GRADX
C           -- ESCF WAS ZERO IF SCF NOT CONVERGED --
            IF(ABS(ESCF).LT.1.0D-12) THEN
               IF(MASWRK) WRITE(IW,'(/A/)')
     *            ' ERROR: SCF WAS NOT CONVERGED. QUANPOL MUST STOP.'
               CALL ABRT
            END IF
            CALL DCOPY(3*NAT,EG,1,QMGRD2,1)
            CALL ELMOMC
            DMXB = DMX
            DMYB = DMY
            DMZB = DMZ
            CALL SAVEABPROP(2)
            CALL SETQMAB(1)
            CALL SETFFDATAB(1)
            CALL MIXABPROP
            DMX = (1-WSIMUL)*DMXA + WSIMUL*DMXB
            DMY = (1-WSIMUL)*DMYA + WSIMUL*DMYB
            DMZ = (1-WSIMUL)*DMZA + WSIMUL*DMZB
         END IF
         E      = ETOT
C
C        - COMBINE QM AND MM GRADIENTS -
C
         IF(LISTQM(NFFAT+NAT+1).GT.0)THEN
            DO KAT = 1, NAT
               KFFAT = LISTQM(NFFAT+KAT)
               IF(KFFAT.GT.0) THEN
                  QMGRD2(1,KAT)   = QMGRD2(1,KAT) + FFGRD2(1,KFFAT)
                  QMGRD2(2,KAT)   = QMGRD2(2,KAT) + FFGRD2(2,KFFAT)
                  QMGRD2(3,KAT)   = QMGRD2(3,KAT) + FFGRD2(3,KFFAT)
                  EG(1,KAT)       = QMGRD2(1,KAT)
                  EG(2,KAT)       = QMGRD2(2,KAT)
                  EG(3,KAT)       = QMGRD2(3,KAT)
                  FFGRD2(1,KFFAT) = QMGRD2(1,KAT)
                  FFGRD2(2,KFFAT) = QMGRD2(2,KAT)
                  FFGRD2(3,KFFAT) = QMGRD2(3,KAT)
               END IF
            ENDDO
         END IF
C        -- CALCULATE DIPOLE MOMENT
         DIPOLE(1)= ZERO
         DIPOLE(2)= ZERO
         DIPOLE(3)= ZERO
         IF(IDOPOL.GT.0) THEN
            DO JFFAT = 1, NFFAT
               DIPOLE(1)=DIPOLE(1) + DIP(1,JFFAT)
               DIPOLE(2)=DIPOLE(2) + DIP(2,JFFAT)
               DIPOLE(3)=DIPOLE(3) + DIP(3,JFFAT)
            ENDDO
         END IF
         DIPOLE(1)=DIPOLE(1)/DB2AU + DMX
         DIPOLE(2)=DIPOLE(2)/DB2AU + DMY
         DIPOLE(3)=DIPOLE(3)/DB2AU + DMZ
C
         C(III,IAT)=CSAVE
         IF(LISTQM(NFFAT+IAT).GT.0)THEN
            CORD(III,LISTQM(NFFAT+IAT)) = C(III,IAT)
         END IF
C
         CALL PUVIB(IFCM,IW,.FALSE.,NCOORD,IVIB,IAT,III,
     *              E,EG,DIPOLE)
         CALL DCOPY(3*NAT,EG,1,VIBGRD(1,III,IAT,IVIB),1)
         VIBDIP(1,III,IAT,IVIB)=DIPOLE(1)
         VIBDIP(2,III,IAT,IVIB)=DIPOLE(2)
         VIBDIP(3,III,IAT,IVIB)=DIPOLE(3)
  630    CONTINUE
      ENDDO
  610 CONTINUE
      ENDDO
C
      DO IAT = 1, NAT
         IYES = 1
         IF(NACTQM.GT.0) THEN
            IYES = 0
            DO LLL = 1, NACTQM
               IF(IAT.EQ.LACTQM(LLL)) IYES = 1
            ENDDO
         END IF
         DO III = 1, 3
            DO JAT = 1, NAT
               JYES = 1
               IF(NACTQM.GT.0) THEN
                  JYES = 0
                  DO LLL = 1, NACTQM
                     IF(JAT.EQ.LACTQM(LLL)) JYES = 1
                  ENDDO
               END IF
               DO JJJ = 1, 3
                  VIBHSS((JAT-1)*3+JJJ,(IAT-1)*3+III)
     *            =(VIBGRD((JAT-1)*3+JJJ,III,IAT,1)
     *             -VIBGRD((JAT-1)*3+JJJ,III,IAT,2))/(2*DISPL)
                  IF(IYES+JYES.LT.2)
     *            VIBHSS((JAT-1)*3+JJJ,(IAT-1)*3+III) = ZERO
                  IF(IYES+JYES.EQ.0.AND.IAT.EQ.JAT.AND.III.EQ.JJJ)
     *            VIBHSS((JAT-1)*3+JJJ,(IAT-1)*3+III) =
     *            -3.7843877795D-08*ZMASS(IAT)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      DO IAT = 1, NAT
         IYES = 1
         IF(NACTQM.GT.0) THEN
            IYES = 0
            DO LLL = 1, NACTQM
               IF(IAT.EQ.LACTQM(LLL)) IYES = 1
            ENDDO
         END IF
         DO III = 1, 3
            DO JJJ = 1, 3
               VIBDDM(III,(IAT-1)*3+JJJ)
     *         =(VIBDIP(III,JJJ,IAT,1)
     *          -VIBDIP(III,JJJ,IAT,2))/(2*DISPL)
               IF(IYES.EQ.0) VIBDDM(III,(IAT-1)*3+JJJ) = ZERO
            ENDDO
         ENDDO
      ENDDO
C
  700 CONTINUE
      IF(MASWRK) WRITE(IP,'(A)')
     *'----- HESSIAN MATRIX (NOT MASS WEIGHTED) -----'
      CALL FCMPUN(VIBHSS,NCOORD)
      CALL DAWRIT(IDAF,IODA,VIBHSS,NCOORD*NCOORD,4,0)
      CALL DAWRIT(IDAF,IODA,VIBDDM,3*NCOORD,34,0)
C
C     - VIBANL -
      PROJCT =.TRUE.
      GOTDDM =.TRUE.
      GOTADM =.FALSE.
      IF((.NOT.CVGED).AND.MASWRK) WRITE(IW,9010)
 9010 FORMAT(/
     *   5X,'*******************************************************'/
     *   5X,'* THIS IS NOT A STATIONARY POINT ON THE MOLECULAR PES *'/
     *   5X,'*     THE VIBRATIONAL ANALYSIS IS NOT VALID !!!       *'/
     *   5X,'*******************************************************')
      STATPT =.TRUE.
      SCFOK  =.TRUE.
      NTEMP  = 1
      TEMP(1)= TEMP0
      SCLFAC = 1.0D+00
      NIRRED = MXIRR
      DO I =1, MXIRR
         IRPNAM(I)=IRPA
      ENDDO
      EFLDL  = .FALSE.  ! THIS HELPS PROJECT THE ROTATIONS
C
      NPART  = NAT
      NC1    = 3*NAT
      NC2    = (NC1**2+NC1)/2
      NC3    = NC1**2
      CALL VALFM(LOADFM)
      LVEC   = LOADFM + 1
      LFCM   = LVEC   + NC3
      LE     = LFCM   + NC2
      LSCR   = LE     + NC1
      LIA    = LSCR   + NC1*8
      LRM    = LIA    + NC1
      LSVT   = LRM    + NCOORD
      LSVR   = LSVT   + NC1*3
      LSVTT  = LSVR   + NC1*3
      LSVRT  = LSVTT  + NC1
      LCC    = LSVRT  + NC1
      LCOM   = LCC    + 3*NPART
      LZMS   = LCOM   + 3*NPART
      LBUF1  = LZMS   +   NPART
      LBUF2  = LBUF1  + NC3
      LDDM   = LBUF2  + NC3
      LADM   = LDDM   + NC1*3
      LAST   = LADM   + NC1*6
      NEED   = LAST   - LOADFM
      CALL GETFM(NEED)
      CALL FGMTRX(X(LVEC),X(LFCM),X(LE),X(LSCR),X(LIA),X(LRM),
     *            X(LSVT),X(LSVR),X(LSVTT),X(LSVRT),
     *            freq,X(LCC),X(LCOM),X(LZMS),X(LBUF1),X(LBUF2),
     *            NC1,NC2,NPART,X(LDDM),GOTDDM,X(LADM),GOTADM,
     *            PROJCT,STATPT,SCFOK,0,0)
      CALL RETFM(NEED)
      IF((.NOT.CVGED).AND.MASWRK) WRITE(IW,9010)
      IF (MASWRK) WRITE(IW,9020)
      CALL TIMIT(1)
C
 9020 FORMAT(1X,'......END OF NORMAL COORDINATE ANALYSIS......')
C
      RETURN
      END
C*MODULE QUANPOC  *DECK HESSQMMM
!>
!> @brief    QM/MM Hessian
!>
!> @author   Hui Li
!>           - Oct 2016
!>
!> @details  single displacement of selected QM atoms
!>           if NSTEP=2, double displacement
!>
      SUBROUTINE HESSQMMM(CORD,CORDSV,
     *                  FFGRD2,QMGRD2,
     *                  LISTQM,NONLS1,NONLSTQ,
     *                  LSTCELL,NONLS2,CORDSV2,CORDSVQ,
     *                  MVFASTS2,MVFASTS3,MVFASTS4,
     *                  MVFASTL2,MVFASTL3,MVFASTL4,
     *                  NONLSA,NONLSB,
     *                  NONLSPMA,NONLSPMB,
     *                  VIBGRD,VIBE,VIBGRD0)
      use mx_limits, only: mxatm,mxao,mxrt
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      LOGICAL GOPARR,DSKWRK,MASWRK
C     LOGICAL GOTEH,EOF,RSTART
      LOGICAL PROJCT,GOTDDM,GOTADM,STATPT,SCFOK,
     *        PRTSCN,GOTFRQ,GOTEH,EOF,RSTART
      LOGICAL EFLDL
C
      PARAMETER (MXIRR=14)
      PARAMETER (ZERO=0.0D+00)
C
      CHARACTER*10  TXT1,TXT2,TXT3,TXT4
C
      DIMENSION CORD(3,NFFAT),CORDSV(3,NFFAT),
     *          FFGRD2(3,NFFAT),QMGRD2(3,MXATM),
     *          LISTQM(*),NONLSTQ(*)
      DIMENSION VIBGRD(3*NACTMM,3,NACTMM,2),VIBGRD0(3*NACTMM)
      DIMENSION DIPOLE(3),VIBE(3,MXATM,2)
C
      COMMON /EFLDC / EVEC(3),EFLDL
      COMMON /ENRGYS/ ENUCR,EELCT,ETOT,SZ,SZZ,ECORE,ESCF,EERD,E1,E2,
     *                VEN,VEE,EPOT,EKIN,ESTATE(MXRT),STATN,EDFT(3),EDISP
      COMMON /FFFREE/ SOL1CH,SOL1LJ,SOL1IM,SOLFRE1MM,
     *                SOL2CH,SOL2LJ,SOL2IM,SOLFRE2MM,
     *                ASOL1CH,ASOL1LJ,ASOL1IM,ASOL1MM,
     *                ASOL2CH,ASOL2LJ,ASOL2IM,ASOL2MM,
     *                PMF1BD,PMF1AG,PMF1DR,PMF1DB,PMF1WG,PMF1CM,
     *                PMF1CH,PMF1PO,PMF1LJ,PMF1IM,PMFFRE1MM,
     *                EN12A,EN12B,EN123A,EN123B,EN123R4A,EN123R4B,
     *                EN234W1A,EN234W1B,EN123B4A,EN123B4B,
     *                EN12312A,EN12312B,
     *                ENCHARA,ENCHARB,ENPOLA,ENPOLB,
     *                ENRXNA,ENRXNB,ENRXNPOLA,ENRXNPOLB,
     *                ENRXNRA,ENRXNRB,ENLJRA,ENLJRB,ENLJDA,ENLJDB,
     *                ENCMAPA,ENCMAPB,ENFIXSOA,ENFIXSOB,
     *                ENCENTA,ENCENTB,WSIMUL,WPERT1,WPERT2,
     *                IFEPTYP
      COMMON /FFPARM/ NFFAT,NBOND,NANGL,NDIHR,NDIHB,NCMAP,NWAGG,
     *                N1213J,N14J,NLKQMM,IDOCHG,IDOPOL,IDOLJ,IDOCMAP
      COMMON /FFMDPA/ DT,DT2,TEMP0,PRES0,POLTOL,VIR(3),PMEAN,VOLAV,
     *                ENPAV,ENKAV,TEMPAV,BERENDT,BERENDP,VELMAX,
     *                PMEANX,PMEANY,PMEANZ,NSTEP,KMASTER,KOUTACT(2),
     *                IHESS,INTALG,ITSTAT,IPSTAT,JOUT,KOUT,LOUT
      COMMON /FFMPT3/ NACTMM,LACTMM(2020),LFFDIPOLD,JUMBUP,
     *                NACTQM,LACTQM(2020),LFFOLDC,LFFQMVELSV,MMHESS,
     *                LFFQMCHG,LFFQMCHGB,ISWAP,R2SWAP,DFTBMM
      COMMON /FFQMPA/ ENFFQM2,ENPAV2,
     *                SCFTYP2,TDDFT2,MPLEVL2,CITYP2,
     *                ICHARG2,MULT2,IDOQM2,IREDOX,IQMPKA,IQMRXN,
     *                MATOMA,MCHARGA,MULTA,MELEA,
     *                MATOMB,MCHARGB,MULTB,MELEB,
     *                IECPX,NSHELLX,IMP,JMP,ICORSH,IGTF,
     *                LFFZANX,LFFCLPX,LFFZLPX,LFFNLPX,LFFKFRSTX,
     *                LFFKLASTX,LFFLMAXX,LFFLPSKIPX,LFFIZCOREX,
     *                LFFCX,LFFIANX,LFFEXX,LFFCSX,LFFCPX,LFFCDX,
     *                LFFCFX,LFFCGX,LFFCHX,LFFCIX,LFFKSTARTX,
     *                LFFKATOMX,LFFKTYPEX,LFFKNGX,LFFKLOCX,
     *                LFFMINX,LFFMAXX,LFFMPTYPX,LFFAN0X,
     *                LFFALPN0X,LFFAN1X,LFFALPN1X,LFFMPSKPX,
     *                LFFNOAN0X,LFFNOAN1X,LFFBPARX,LFFEXPMPX,
     *                LFFCSMPX,LFFCPMPX,LFFCDMPX,LFFCFMPX,LFFMPSKIPX,
     *                LFFNOCOSHX,LFFMPKSTAX,LFFMPKNGX,LFFMPKTYPX,
     *                LFFMPKMINX,LFFMPKMAXX,LFFMPKLOCX,LFFANAMX
      COMMON /FMCOM / X(1)
      COMMON /FUNCT / E,EG(3,MXATM)
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /IOFILE/ IR,IW,IP,IJK,IJKT,IDAF,NAV,IODA(950)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /RUNLAB/ TITLE(10),ANAM(MXATM),BNAM(MXATM),BFLAB(MXAO)
      COMMON /SCFOPT/ CONVHF,MAXIT,MCONV,NPUNCH,NPREO(4),FSHIFT
      COMMON /SYMBLK/ NIRRED,NSALC,NSALC2,NSALC3,NSAFMO
      COMMON /SYMREP/ IRPNAM(MXIRR),IPA(MXIRR),LAMBDA(MXIRR),
     *                LAMBD0(MXIRR),IADDR1(MXIRR),IADDR2(MXIRR),
     *                IADDR3(MXIRR)
      COMMON /THERMD/ FREQ(3*MXATM),TEMP(10),SCLFAC,NTEMP,PRTSCN,GOTFRQ
      COMMON /OUTPUT/ NPRINT,ITOL,ICUT,NORMF,NORMP,NOPK
C
      DATA IRPA /4HA   /
C
C     -- READ IN A QM/MM HESSIAN, UPDATE THE ELEMENTS FOR SELECTED QM ATOMS
C        SO THE NEW QM/MM HESSIAN CAN BE USED FOR TS SEARCH.
C        MM HESSIAN CANNOT BE USED, ONLY QMMM HESSIAN CAN BE USED.
C        HUI LI, OCT 04, 2016, LINCOLN
C
      IF(NSTEP.LT.0) RETURN
C
      CALL VCLR(VIBGRD,1,3*NACTMM*3*NACTMM*2)
      CALL VCLR(VIBE  ,1,3*MXATM*2)
C
C     -- DETECT $HESS --
C
      CALL SEQREW(IR)
      CALL FNDGRP(IR,' $HESS  ',IEOF)
      IF (IEOF.EQ.1) THEN
         GOTEH=.FALSE.
         IF(MASWRK) WRITE(IW,'(/1X,A/A/)')
     *   'ERROR: NO $HESS IS FOUND IN INPUT FILE FOR UPDATING.',
     *   '       A REGULAR HESSIAN JOB SHOULD NOT HAVE LACTMM.'
         CALL ABRT
      ELSE
         GOTEH=.TRUE.
         IF (MASWRK) THEN
            READ(IR,*)TXT1,TXT2,AEN1,TXT3,TXT4,AEN2
            IF(TXT1.EQ.'ENERGY'.AND.TXT2.EQ.'IS'.AND.
     *         TXT3.EQ.'E(NUC)'.AND.TXT4.EQ.'IS') THEN
            IF(AEN2.EQ.ZERO) THEN
               WRITE(IW,'(/1X,A,A,/)')'ERROR: MM $HESS CANNOT BE USED.',
     *         ' ONLY QMMM $HESS CAN BE USED.'
               CALL ABRT
            END IF
            IF(AEN1.LT.ZERO.AND.AEN2.GT.ZERO.AND.MMHESS.EQ.1) THEN
               WRITE(IW,'(/1X,A,/)')'A GOOD QMMM $HESS IS IN INPUT.'
            END IF
            END IF
         END IF
      END IF
C
      DIPOLE(1)= ZERO
      DIPOLE(2)= ZERO
      DIPOLE(3)= ZERO
C
      NPRINT  = -5
      NATSV   = NAT
      IF(NACTMM.EQ.0) THEN
         DO IFFAT =1,NFFAT
            IF(LISTQM(IFFAT).EQ.0.OR.LISTQM(IFFAT).GT.NATSV) THEN
               NAT = NAT + 1
               C(1,NAT)=CORD(1,IFFAT)
               C(2,NAT)=CORD(2,IFFAT)
               C(3,NAT)=CORD(3,IFFAT)
            END IF
         ENDDO
      ELSE
         DO I=1,NACTMM
            IFFAT  = LACTMM(I)
            IF(LISTQM(IFFAT).EQ.0.OR.LISTQM(IFFAT).GT.NATSV) THEN
               NAT = NAT + 1
               C(1,NAT) = CORD(1,IFFAT)
               C(2,NAT) = CORD(2,IFFAT)
               C(3,NAT) = CORD(3,IFFAT)
            END IF
         ENDDO
      END IF
      IF(NAT.GT.2000) THEN
         IF(MASWRK)WRITE(IW,*)
     *   'ERROR: NUMBER OF OPTIMIZE ATOMS EXCEEDED 2000.'
         IF(MASWRK)WRITE(IW,*)' '
         CALL ABRT
      END IF
C
      NATTOT=NAT
      NCOORD=NATTOT*3
      IF(NATTOT.NE.NACTMM)THEN
         IF(MASWRK)WRITE(IW,*)
     *   'ERROR: LACTMM LIST SHOULD INCLUDE ALL QM ATOMS.'
         IF(MASWRK)WRITE(IW,*)' '
         CALL ABRT
      END IF
C     -- MUST SAVE THE VIB 0 GRADIENT
C        HERE TO OCCUPY THE DAF SPACE
      CALL DAWRIT(IDAF,IODA,EG,NCOORD, 3,0)
C
C     -- VIB 0 --
C
      NPUNCH = 2
C
      IF(NSTEP.EQ.0) GOTO 700
C
C     -- CALCULATE ENERGY AND GRADIENT -
      NAT    = NATSV    !  BECAUSE NONBOND USES NAT
      CALL NONBOND(0,CORD,CORDSV,CORDSV2,CORDSVQ,
     *             NONLS1,NONLS2,
     *             NONLSTQ,LSTCELL,
     *             MVFASTS2,MVFASTS3,MVFASTS4,
     *             MVFASTL2,MVFASTL3,MVFASTL4,
     *             NONLSA,NONLSB,
     *             NONLSPMA,NONLSPMB)
      CALL VCLR(QMGRD2,1,3*NAT)
      CALL VCLR(FFGRD2,1,3*NFFAT)
      CALL VCLR(VIR   ,1,3)
      CALL GRADX
C     -- ESCF WAS ZERO IF SCF NOT CONVERGED --
      IF(ABS(ESCF).LT.1.0D-12) THEN
         IF(MASWRK) WRITE(IW,'(/A/)')
     *      ' ERROR: SCF WAS NOT CONVERGED. QUANPOL MUST STOP.'
         CALL ABRT
      END IF
      CALL DCOPY(3*NAT,EG,1,QMGRD2,1)
      IF(IFEPTYP.GT.0.AND.MATOMB.GT.0) THEN
         CALL SAVEABPROP(1)
         CALL DCOPY(3*MATOMA,C,1,X(LFFCX),1)
         CALL SAVEFFDATA
         CALL SETQMAB(2)
         CALL SETFFDATAB(2)
         CALL VCLR(QMGRD2,1,3*NAT  )
         CALL VCLR(FFGRD2,1,3*NFFAT)
         CALL VCLR(VIR   ,1,3)
         CALL GRADX
C        -- ESCF WAS ZERO IF SCF NOT CONVERGED --
         IF(ABS(ESCF).LT.1.0D-12) THEN
            IF(MASWRK) WRITE(IW,'(/A/)')
     *         ' ERROR: SCF WAS NOT CONVERGED. QUANPOL MUST STOP.'
            CALL ABRT
         END IF
         CALL DCOPY(3*NAT,EG,1,QMGRD2,1)
         CALL SAVEABPROP(2)
         CALL SETQMAB(1)
         CALL SETFFDATAB(1)
         CALL MIXABPROP
      END IF
      E      = ETOT
      NAT    = NATTOT
C
C     - COMBINE QM AND MM GRADIENTS -
C
      IF(LISTQM(NFFAT+NATSV+1).GT.0)THEN
         DO IAT = 1, NATSV
            KFFAT = LISTQM(NFFAT+IAT)
            IF(KFFAT.GT.0) THEN
               QMGRD2(1,IAT)   = QMGRD2(1,IAT) + FFGRD2(1,KFFAT)
               QMGRD2(2,IAT)   = QMGRD2(2,IAT) + FFGRD2(2,KFFAT)
               QMGRD2(3,IAT)   = QMGRD2(3,IAT) + FFGRD2(3,KFFAT)
               EG(1,IAT)       = QMGRD2(1,IAT)
               EG(2,IAT)       = QMGRD2(2,IAT)
               EG(3,IAT)       = QMGRD2(3,IAT)
               FFGRD2(1,KFFAT) = QMGRD2(1,IAT)
               FFGRD2(2,KFFAT) = QMGRD2(2,IAT)
               FFGRD2(3,KFFAT) = QMGRD2(3,IAT)
            END IF
         ENDDO
      END IF
C     - EXPAND EG(3,*) TO INCLUDE MM ATOMS -
      IAT = NATSV
      IF(NACTMM.EQ.0) THEN
         DO IFFAT = 1, NFFAT
            IF(LISTQM(IFFAT).EQ.0.OR.LISTQM(IFFAT).GT.NATSV)THEN
               IAT = IAT + 1
               EG(1,IAT) = FFGRD2(1,IFFAT)
               EG(2,IAT) = FFGRD2(2,IFFAT)
               EG(3,IAT) = FFGRD2(3,IFFAT)
            END IF
         ENDDO
      ELSE
         DO I=1,NACTMM
            IFFAT  = LACTMM(I)
C           - LACTMM DOES NOT AFFECT PURE QM EG()
C           - LACTMM DOES NOT AFFECT QM(MM)  EG()
            IF(LISTQM(IFFAT).EQ.0.OR.LISTQM(IFFAT).GT.NATSV) THEN
               IAT = IAT + 1
               EG(1,IAT) = FFGRD2(1,IFFAT)
               EG(2,IAT) = FFGRD2(2,IFFAT)
               EG(3,IAT) = FFGRD2(3,IFFAT)
            END IF
         ENDDO
      END IF
C
C     -- MUST SAVE THE VIB 0 GRADIENT
      CALL DAWRIT(IDAF,IODA,EG,NCOORD, 3,0)
C
      CALL TIMIT(1)
C
      IF(MASWRK) WRITE(IW,'(//1X,A/1X,A/1X,A//)')
     *'********************************************',
     *'QUANPOL QM/MM HESSIAN CALCULATION: QM UPDATE',
     *'********************************************'
C
C     -- TRY TO READ IN $VIB FROM THE INPUT FILE
C
      CALL SEQREW(IR)
      CALL FNDGRP(IR,' $VIB   ',IEOF)
      IF(IEOF.EQ.1) THEN
         RSTART = .FALSE.
      ELSE
         RSTART = .TRUE.
         IF(MASWRK) WRITE(IW,9120)
      END IF
 9120 FORMAT(/5X,'A $VIB GROUP WAS FOUND IN YOUR INPUT.'/
     *        5X,'THEREFORE, THIS IS A NUMERICAL HESSIAN RESTART.')
C
      IFCM=35
      IF (MASWRK) CALL SEQOPN(IFCM,'RESTART','NEW',.FALSE.,'FORMATTED')
      IF (MASWRK) WRITE(IFCM,8010) TITLE
 8010 FORMAT('ENERGY/GRADIENT/DIPOLE RESTART DATA FOR RUNTYP=HESSIAN'/
     *         10A8)
      CALL PUVIB(IFCM,IW,.FALSE.,NCOORD,0,0,0,
     *           E,EG,DIPOLE)
      CALL DCOPY(3*NAT,EG,1,VIBGRD0,1)
      IF(RSTART) THEN
         IF(MASWRK) WRITE(IW,9050)
 9050    FORMAT(/5X,'READING $VIB RESTART DATA FROM $VIB CARDS...')
         ECOMPUT = E
         IREAD   = 0
         IVIB0   = 0
         IAT0    = 0
         ICOORD0 = 0
         CALL RDVIB(IR,IW,IREAD,NCOORD,IVIB0,IAT0,ICOORD0,
     *              E,EG,DIPOLE,EOF,DUMY,DUMY)
         IF(EOF) THEN
            IF(MASWRK) WRITE(IW,*) 'PROBLEM READING $VIB GROUP.'
         END IF
         IF(ABS(ECOMPUT-E).GT.1.0D-05) THEN
            IF(MASWRK) WRITE(IW,'(/1X,A,A/)')
     *      'ERROR: $VIB IVIB=0 E= DOES NOT MATCH THE COMPUTED ',
     *      'ENERGY. QUANPOL WILL QUIT.'
            CALL ABRT
         END IF
C
         DO I = 1, 2*3*NATSV
            KVIB   = 0
            KAT    = 0
            KCOORD = 0
            CALL RDVIB(IR,IW,IREAD,NCOORD,KVIB,KAT,KCOORD,
     *                 E,EG,DIPOLE,EOF,DUMY,DUMY)
            IF(EOF) GOTO 500
C           - ONLY MASWRK CHANGES KVIB, KAT, KCOORD IN RDVIB -
            IF(GOPARR) THEN
               CALL DDI_BCAST(476,'I',KVIB  ,1,MASTER)
               CALL DDI_BCAST(477,'I',KAT   ,1,MASTER)
               CALL DDI_BCAST(478,'I',KCOORD,1,MASTER)
            END IF
            IVIB    = KVIB
            IAT     = KAT
            III     = KCOORD
            VIBE(III,IAT,IVIB) = E
            CALL DCOPY(3*NAT,EG,1,VIBGRD(1,III,IAT,IVIB),1)
         ENDDO
  500    CONTINUE
      END IF
C
      NVIB = 1
      IF(NSTEP.EQ.2) NVIB = 2
      DISPL= 0.01D+00
C
      DO IVIB = 1, NVIB
      IF(IVIB.EQ.1) SIGN=+1.0D+00
      IF(IVIB.EQ.2) SIGN=-1.0D+00
      DO 610 IAT = 1, NATSV
      IYES = 1
      IF(NACTQM.GT.0) THEN
         IYES = 0
         DO LLL = 1, NACTQM
            IF(IAT.EQ.LACTQM(LLL)) IYES = 1
         ENDDO
      END IF
      IF(IYES.EQ.0) GOTO 610
      DO III = 1, 3
C
         SUM = ZERO
         DO J=1, 3*NAT
            SUM=SUM+ABS(VIBGRD(J,III,IAT,IVIB))
         ENDDO
         IF(SUM.GT.1.0D-12) THEN
            CALL PUVIB(IFCM,IW,.FALSE.,NCOORD,IVIB,IAT,III,
     *                 VIBE(III,IAT,IVIB),VIBGRD(1,III,IAT,IVIB),
     *                 DIPOLE)
            GOTO 630
         END IF
C
         CSAVE=C(III,IAT)
         C(III,IAT)=C(III,IAT)+SIGN*DISPL
         IF(LISTQM(NFFAT+IAT).GT.0)THEN
            CORD(III,LISTQM(NFFAT+IAT)) = C(III,IAT)
         END IF
C
C        -- CALCULATE ENERGY AND GRADIENT
C           NO NEED TO UPDATE NONBOND LIST
         NAT    = NATSV
         NPUNCH = 0
         CALL VCLR(QMGRD2,1,3*NAT  )
         CALL VCLR(FFGRD2,1,3*NFFAT)
         CALL VCLR(VIR   ,1,3)
         CALL GRADX
C        -- ESCF WAS ZERO IF SCF NOT CONVERGED --
         IF(ABS(ESCF).LT.1.0D-12) THEN
            IF(MASWRK) WRITE(IW,'(/A/)')
     *         ' ERROR: SCF WAS NOT CONVERGED. QUANPOL MUST STOP.'
            CALL ABRT
         END IF
         CALL DCOPY(3*NAT,EG,1,QMGRD2,1)
         IF(IFEPTYP.GT.0.AND.MATOMB.GT.0) THEN
            CALL SAVEABPROP(1)
            CALL DCOPY(3*MATOMA,C,1,X(LFFCX),1)
            CALL SAVEFFDATA
            CALL SETQMAB(2)
            CALL SETFFDATAB(2)
            CALL VCLR(QMGRD2,1,3*NAT  )
            CALL VCLR(FFGRD2,1,3*NFFAT)
            CALL VCLR(VIR   ,1,3)
            CALL GRADX
C           -- ESCF WAS ZERO IF SCF NOT CONVERGED --
            IF(ABS(ESCF).LT.1.0D-12) THEN
               IF(MASWRK) WRITE(IW,'(/A/)')
     *            ' ERROR: SCF WAS NOT CONVERGED. QUANPOL MUST STOP.'
               CALL ABRT
            END IF
            CALL DCOPY(3*NAT,EG,1,QMGRD2,1)
            CALL SAVEABPROP(2)
            CALL SETQMAB(1)
            CALL SETFFDATAB(1)
            CALL MIXABPROP
         END IF
         E      = ETOT
         NAT    = NATTOT
C
C        - COMBINE QM AND MM GRADIENTS -
C
         IF(LISTQM(NFFAT+NATSV+1).GT.0)THEN
            DO KAT = 1, NATSV
               KFFAT = LISTQM(NFFAT+KAT)
               IF(KFFAT.GT.0) THEN
                  QMGRD2(1,KAT)   = QMGRD2(1,KAT) + FFGRD2(1,KFFAT)
                  QMGRD2(2,KAT)   = QMGRD2(2,KAT) + FFGRD2(2,KFFAT)
                  QMGRD2(3,KAT)   = QMGRD2(3,KAT) + FFGRD2(3,KFFAT)
                  EG(1,KAT)       = QMGRD2(1,KAT)
                  EG(2,KAT)       = QMGRD2(2,KAT)
                  EG(3,KAT)       = QMGRD2(3,KAT)
                  FFGRD2(1,KFFAT) = QMGRD2(1,KAT)
                  FFGRD2(2,KFFAT) = QMGRD2(2,KAT)
                  FFGRD2(3,KFFAT) = QMGRD2(3,KAT)
               END IF
            ENDDO
         END IF
C        - EXPAND EG(3,*) TO INCLUDE MM ATOMS -
         KAT = NATSV
         IF(NACTMM.EQ.0) THEN
            DO IFFAT = 1, NFFAT
               IF(LISTQM(IFFAT).EQ.0.OR.LISTQM(IFFAT).GT.NATSV)THEN
                  KAT = KAT + 1
                  EG(1,KAT) = FFGRD2(1,IFFAT)
                  EG(2,KAT) = FFGRD2(2,IFFAT)
                  EG(3,KAT) = FFGRD2(3,IFFAT)
               END IF
            ENDDO
         ELSE
            DO I=1,NACTMM
               IFFAT  = LACTMM(I)
C              - LACTMM DOES NOT AFFECT PURE QM EG()
C              - LACTMM DOES NOT AFFECT QM(MM)  EG()
               IF(LISTQM(IFFAT).EQ.0.OR.LISTQM(IFFAT).GT.NATSV) THEN
                  KAT = KAT + 1
                  EG(1,KAT) = FFGRD2(1,IFFAT)
                  EG(2,KAT) = FFGRD2(2,IFFAT)
                  EG(3,KAT) = FFGRD2(3,IFFAT)
               END IF
            ENDDO
         END IF
C
         C(III,IAT)=CSAVE
         IF(LISTQM(NFFAT+IAT).GT.0)THEN
            CORD(III,LISTQM(NFFAT+IAT)) = C(III,IAT)
         END IF
C
         CALL PUVIB(IFCM,IW,.FALSE.,NCOORD,IVIB,IAT,III,
     *              E,EG,DIPOLE)
         CALL DCOPY(3*NAT,EG,1,VIBGRD(1,III,IAT,IVIB),1)
  630    CONTINUE
      ENDDO
  610 CONTINUE
      ENDDO
C
  700 CONTINUE
C
C     -- READ IN $HESS --
      GOTEH = .FALSE.
      CALL VALFM(LOADFM)
      LHESS   = LOADFM + 1
      LAST    = LHESS  + NCOORD*NCOORD
      NEED    = LAST   - LOADFM -1
      CALL GETFM(NEED)
      CALL FCMIN(X(LHESS),NCOORD,GOTEH)
      IF(GOTEH)THEN
      ELSE
         IF(MASWRK) WRITE(IW,'(/1X,A/)')
     *   'ERROR: NO $HESS IS FOUND IN INPUT FILE FOR UPDATING.'
         CALL ABRT
      END IF
C
      IF(NSTEP.EQ.0) GOTO 800
C
      DO IAT = 1, NAT
         IYES = 1
         IF(NACTQM.GT.0) THEN
            IYES = 0
            DO LLL = 1, NACTQM
               IF(IAT.EQ.LACTQM(LLL)) IYES = 1
            ENDDO
         END IF
         IF(IYES.EQ.1) THEN
         DO III = 1, 3
            DO JAT = 1, NAT    !  ALL NATVIB QM+MM ATOMS
               DO JJJ = 1, 3
                 IF(NVIB.EQ.1) THEN
                  X(LHESS+((IAT-1)*3+III-1)*NCOORD+(JAT-1)*3+JJJ-1)
     *            =(VIBGRD((JAT-1)*3+JJJ,III,IAT,1)
     *             -VIBGRD0((JAT-1)*3+JJJ))/DISPL
                  X(LHESS+((JAT-1)*3+JJJ-1)*NCOORD+(IAT-1)*3+III-1)
     *           =X(LHESS+((IAT-1)*3+III-1)*NCOORD+(JAT-1)*3+JJJ-1)
                 END IF
                 IF(NVIB.EQ.2) THEN
                  X(LHESS+((IAT-1)*3+III-1)*NCOORD+(JAT-1)*3+JJJ-1)
     *            =(VIBGRD((JAT-1)*3+JJJ,III,IAT,1)
     *             -VIBGRD((JAT-1)*3+JJJ,III,IAT,2))/(2*DISPL)
                  X(LHESS+((JAT-1)*3+JJJ-1)*NCOORD+(IAT-1)*3+III-1)
     *           =X(LHESS+((IAT-1)*3+III-1)*NCOORD+(JAT-1)*3+JJJ-1)
                 END IF
               ENDDO
            ENDDO
         ENDDO
         END IF
      ENDDO
C
      CALL SEQOPN(39,'OPTHES1','NEW',.FALSE.,'FORMATTED')
      IF(MASWRK) WRITE(39,*)
     *'----- UPDATED HESSIAN MATRIX (NOT MASS WEIGHTED) -----'
      IPSV = IP
      IP   = 39
      CALL FCMPUN(X(LHESS),NCOORD)
      IP   = IPSV
      CALL FLSHBF(39)
      CALL SEQCLO(39,'KEEP')
C
      IF(MASWRK) WRITE(IW,'(/1X,A/)')
     * 'THE UPDATED HESSIAN IS IN THE *.HS1 FILE.'
C
  800 CONTINUE
C
C     - SHRINK THE QM/MM FCM TO BE QM FCM -
      DO IAT = 1, NATSV
         DO III = 1, 3
            DO JAT = 1, NATSV
               DO JJJ = 1, 3
                  X(LHESS+((IAT-1)*3+III-1)*NATSV*3+(JAT-1)*3+JJJ-1)
     *           =X(LHESS+((IAT-1)*3+III-1)*NCOORD+(JAT-1)*3+JJJ-1)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      CALL DAWRIT(IDAF,IODA,X(LHESS),NATSV*NATSV*3*3,4,0)
      CALL RETFM(NEED)
      NAT    = NATSV
      NCOORD = NAT*3
C
C     - VIBANL -
      PROJCT =.TRUE.
      GOTDDM =.FALSE.
      GOTADM =.FALSE.
      IF(MASWRK) WRITE(IW,9010)
 9010 FORMAT(/
     *   5X,'*******************************************************'/
     *   5X,'* THIS IS NOT A STATIONARY POINT ON THE MOLECULAR PES *'/
     *   5X,'*     THE VIBRATIONAL ANALYSIS IS NOT VALID !!!       *'/
     *   5X,'*******************************************************')
      STATPT =.TRUE.
      SCFOK  =.TRUE.
      NTEMP  = 1
      TEMP(1)= TEMP0
      SCLFAC = 1.0D+00
      NIRRED = MXIRR
      DO I =1, MXIRR
         IRPNAM(I)=IRPA
      ENDDO
      EFLDL  = .FALSE.  ! THIS HELPS PROJECT THE ROTATIONS
C
      NPART  = NAT
      NC1    = 3*NAT
      NC2    = (NC1**2+NC1)/2
      NC3    = NC1**2
      CALL VALFM(LOADFM)
      LVEC   = LOADFM + 1
      LFCM   = LVEC   + NC3
      LE     = LFCM   + NC2
      LSCR   = LE     + NC1
      LIA    = LSCR   + NC1*8
      LRM    = LIA    + NC1
      LSVT   = LRM    + NCOORD
      LSVR   = LSVT   + NC1*3
      LSVTT  = LSVR   + NC1*3
      LSVRT  = LSVTT  + NC1
      LCC    = LSVRT  + NC1
      LCOM   = LCC    + 3*NPART
      LZMS   = LCOM   + 3*NPART
      LBUF1  = LZMS   +   NPART
      LBUF2  = LBUF1  + NC3
      LDDM   = LBUF2  + NC3
      LADM   = LDDM   + NC1*3
      LAST   = LADM   + NC1*6
      NEED   = LAST   - LOADFM
      CALL GETFM(NEED)
      CALL FGMTRX(X(LVEC),X(LFCM),X(LE),X(LSCR),X(LIA),X(LRM),
     *            X(LSVT),X(LSVR),X(LSVTT),X(LSVRT),
     *            freq,X(LCC),X(LCOM),X(LZMS),X(LBUF1),X(LBUF2),
     *            NC1,NC2,NPART,X(LDDM),GOTDDM,X(LADM),GOTADM,
     *            PROJCT,STATPT,SCFOK,0,0)
      NAT    = NATTOT
      NCOORD = NAT*3
      CALL RETFM(NEED)
      IF(MASWRK) WRITE(IW,9010)
      IF(MASWRK) WRITE(IW,9020)
      IF(MASWRK) WRITE(IW,'(/1X,A/1X,A,A/)')
     * 'THE UPDATED HESSIAN IS IN THE *.HS1 FILE.',
     * 'PLEASE CHECK THE FREQUENCIES.  FOR TS SEARCH, ONE ',
     * 'IMAGINARY FREQUENCY IS REQUIRED.'
      CALL TIMIT(1)
C
 9020 FORMAT(1X,'......END OF NORMAL COORDINATE ANALYSIS......')
C
      RETURN
      END
C*MODULE QUANPOC  *DECK OPTPROP
!>
!> @brief    properties in optimization process
!>
!> @author   Hui Li
!>           - Apr 2011
!>
!> @details  calculate properties in optimization
!>
      SUBROUTINE OPTPROP(ISTEP,ICONV)
      use mx_limits, only: mxatm,mxrt
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      LOGICAL TRIPLET,SG1T,TAMMD,TPA,ALPHKWD,BETAKWD
      LOGICAL MREKT,MRDEA
      LOGICAL GOPARR,DSKWRK,MASWRK
C
      PARAMETER (PI=3.14159265358979323846264338D+00)
      PARAMETER (TOKCAL=627.509469D+00)
      PARAMETER (TOANGS=0.52917724924D+00)
      PARAMETER (FOURTHIRD=4.0D+00/3.0D+00)
C
      COMMON /ENRGMP/ EMP2,EMP3,EMP4,EMP2A
      COMMON /ENRGYS/ ENUCR,EELCT,ETOT,SZ,SZZ,ECORE,ESCF,EERD,E1,E2,
     *                VEN,VEE,EPOT,EKIN,ESTATE(MXRT),STATN,EDFT(3),EDISP
      COMMON /FFENGY/ EN12,EN123,EN123R4,EN123B4,EN234W1,ENCHAR,ENLJR,
     *                ENLJD,ENPOL,XENPOL,ENRXN,XENRXN,ENRXNPOL,ENRXNR,
     *                EN12312,ENQUANP(30),
     *                ENBIAS,ENCENT,ENUCCH,ENCMAP,ENPOT,ENKIN,ENTOT
      COMMON /FFPARM/ NFFAT,NBOND,NANGL,NDIHR,NDIHB,NCMAP,NWAGG,
     *                N1213J,N14J,NLKQMM,IDOCHG,IDOPOL,IDOLJ,IDOCMAP
      COMMON /FFMDPA/ DT,DT2,TEMP0,PRES0,POLTOL,VIR(3),PMEAN,VOLAV,
     *                ENPAV,ENKAV,TEMPAV,BERENDT,BERENDP,VELMAX,
     *                PMEANX,PMEANY,PMEANZ,NSTEP,KMASTER,KOUTACT(2),
     *                IHESS,INTALG,ITSTAT,IPSTAT,JOUT,KOUT,LOUT
      COMMON /FFPBSW/ XBOX,YBOX,ZBOX,SWF,SWFDX,SWFDY,SWFDZ,
     *                SWRA,ONESWRA,SWRA2,ONESWRA2,
     *                SWRB,ONESWRB,SWRB2,ONESWRB2,ONESWRB4,
     *                SWFDUM3,SWFDUM4,SWFDUM5,
     *                SWRAQ,ONESWRAQ,SWRAQ2,ONESWRAQ2,
     *                SWRBQ,ONESWRBQ,SWRBQ2,ONESWRBQ2,ONESWRBQ4,
     *                SWFDUM3Q,SWFDUM4Q,SWFDUM5Q,
     *                QMSIZE,QMCX,QMCY,QMCZ,QMCXSV,QMCYSV,QMCZSV,
     *                CENTX,CENTY,CENTZ,BUFWID1,BUFWID2,RDAMP,
     *                EFIELDX,EFIELDY,EFIELDZ,QMCXSV2,QMCYSV2,QMCZSV2,
     *                EPS1RB,EPS1RB3,ONEXBOX,ONEYBOX,ONEZBOX,
     *                LQMCT,MXLIST1,NTODO,NTODOSV,NTODOQ,
     *                ISWITCH,ISHIFT,IPOLSHF,
     *                LFFLSTCELL,LFFCORDSV,
     *                LFFPOLSV,LFFCORDSV2,LFFNONLS2,LFFCORDSVQ,
     *                LFFMVFASTS2,LFFMVFASTS3,LFFMVFASTS4,
     *                LFFMVFASTL2,LFFMVFASTL3,LFFMVFASTL4,
     *                MXCHECK,MXLIST2,NTODO2,NTODO2SV
      COMMON /FFFIXS/ ENFIXSO,FIXEPS,FIXTOL,FIXA,FIXQ,RALLMM,RALLQM,
     *                RADMM(200),RADQM(200),NRADMM,NRADQM,IFIXSOL,
     *                LFFDAI,LFFDAIT,LFFIDDAI,LFFIDTMP,LFFTMPTS,
     *                LFFAFIX,LFFIDATOM,LFFRFIX,LFFQFIX,NTSATM,
     *                LFFQFIXMP,LFFQFIXTA,LFFQFIXXY,
     *                LFFXTSFIX,LFFYTSFIX,LFFZTSFIX,
     *                LFFVFIX1,LFFVFIX2,NCYCLE,MXFFTS,NFFTS
      COMMON /FFMPNT/ LFFATMNAM,LFFCORD,LFFZANF,
     *                LFFZMAS,LFFONEMAS,LFFQMZMAS,LFFQM1MAS,
     *                LFFCHARG,LFFPOL,LFFDIP,
     *                LFFFIELD1,LFFFIELD2,LFFFIELD3,
     *                LFFSIG,LFFEPS,LFFSIG2,LFFEPS2,
     *                LFFBOND0,LFFFCBOND,
     *                LFFANGL0,LFFFCANGL,LFFFCWAGG,
     *                LFFDIHB0,LFFFCDIHB,
     *                LFFVROT,LFFNNN,LFFGAMA,LFFIPAIR,
     *                LFFKLIST,LFFLLIST,LFFL1213J,LFFL14J,
     *                LFFMLIST,LFFNLIST,LFFLKQMMM,
     *                LFFVEL,LFFQMVEL,
     *                LFFFFGRD0,LFFFFGRD1,LFFFFGRD2,
     *                LFFQMGRD0,LFFQMGRD1,LFFQMGRD2,LFFDETMP,
     *                LFFCLPR,LFFZLPR,LFFNLPR,
     *                LFFXTS,LFFYTS,LFFZTS,LFFCMAT1,
     *                LFFQRXN1,LFFQRXN2,LFFPOT1,LFFPOT2,LFFQRXNMP,
     *                LFFQRXNTA,LFFQRXNXY,LFFNONLSTQ,
     *                LFFDIPMP,LFFDIPTA,LFFDIPXY,LFFLISTQM,LFFNONLS1,
     *                LFFMAPLST,LFFCMAPCO
      COMMON /FFRMSD/ DIPT(8),TIMGYRA,TIMRALL,
     *                NATPDB,NGYRA,NDIEL,NRALL,
     *                LFFRALL0,NRMSD,LFFRMSD0,KFREEAB(201),
     *                NRIJMM,IJRMM(2,100),NRIJQM,IJRQM(2,100),
     *                NAIJKMM,IJKMM(3,100),NAIJKQM,IJKQM(3,100),
     *                NFIXMM,IFIXMM(200),NFIXQM,IFIXQM(200)
      COMMON /FFSPH / SPHRAD,SPHSIG,SPHEPS,IADDWAT
      COMMON /FFUMBR/ UMBFC,UMBR0,UMBSIZE,
     *                NUMBBIN,NUMBATM(6),NUMBTYP,LFFUMBHIS,
     *                UM2FC,UM2R0,UM2SIZE,
     *                NUM2BIN,NUM2ATM(6),NUM2TYP,LFFUM2HIS
      COMMON /FMCOM / X(1)
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /INFOTD/ CNVTOL,PFREQ(2),MODTD,
     *                JANST,NRADT,NTHET,NPHIT,NLEBT,
     *                NSTAT,NTRIAL,MAXVEC,NTHST,IRECTD,ITDFG,ITDPRP,
     *                TRIPLET,SG1T,NONEQR,TAMMD,TPA,ALPHKWD,BETAKWD,
     *                SPCP(3),MULTD,MREKT,MRDEA,MTHST,IFEDAT(4)
      COMMON /IOFILE/ IR,IW,IP,IJK,IJKT,IDAF,NAV,IODA(950)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /WFNOPT/ SCFTYP,VBTYP,DFTYPE,TDDFTYP,CITYP,CCTYP,
     *                MPLEVL,MPCTYP
C
      DATA RNONE/8HNONE    /
C
C     HUI LI, APR 2011, LINCOLN
C
      VOL   = MIN(XBOX*YBOX*ZBOX, FOURTHIRD*PI*SPHRAD**3)
C
      IF(ISTEP.EQ.0.AND.MASWRK) THEN
         WRITE(IW,*)' '
         WRITE(IW,*)'======================= QUANPOL OPTIMIZATION',
     *              ' INITIATED ========================'
         WRITE(IW,*)' '
      END IF
C
      IF(IFIXSOL.GT.0) THEN
         IF(MASWRK.AND.
     *      (MOD(ISTEP,JOUT).EQ.0.OR.ISTEP.EQ.NSTEP))
     *      WRITE(IW,'(A,F20.10,A,A,I8)')
     *      ' FIXSOL TOTAL SURFACE AREA =',FIXA,' A**2,',
     *      ' NFFTS=',NFFTS
         IF(MASWRK.AND.NCYCLE.LT.200.AND.
     *      (MOD(ISTEP,JOUT).EQ.0.OR.ISTEP.EQ.NSTEP))
     *      WRITE(IW,'(A,I3,A,F11.6)')
     *      ' FIXSOL CONVERGED IN ',NCYCLE,
     *      ' ITERATIONS, TOTAL SURFACE CHARGE=',FIXQ
         IF(MASWRK.AND.NCYCLE.EQ.200)
     *       WRITE(IW,'(A,I3,A,F10.6,A,F12.10)')
     *      ' FIXSOL NOT CONVERGED IN ',NCYCLE,
     *      ' ITERATIONS.  TOTAL SURFACE CHARGE=',FIXQ
      END IF
C
      IF(NAT.LE.0) THEN
         ENPOT = EN12 + EN123 + EN123R4 + EN234W1 + EN123B4 + ENCHAR
     *         + ENPOL  + ENRXN  + ENRXNPOL+ENRXNR + ENLJR + ENLJD
     *         + ENCMAP + ENFIXSO+ EN12312
         ENTOT = ENPOT
         ENRXN = ENRXN + ENRXNPOL
         IF(MASWRK.AND.(MOD(ISTEP,JOUT).EQ.0.OR.ISTEP.EQ.NSTEP.OR.
     *      ICONV.EQ.1)) THEN
            WRITE(IW,*)' '
            WRITE(IW,'(1X,A,I10,33X,A,F19.2,A)')'OPT STEP',ISTEP
            WRITE(IW,9000)
     *      'BOND STRETCHING               ENERGY = ', EN12   *TOKCAL
            WRITE(IW,9000)
     *      'BOND ANGLE BENDING            ENERGY = ', EN123  *TOKCAL
            WRITE(IW,9000)
     *      'STRETCHING BENDING            ENERGY = ', EN12312*TOKCAL
            WRITE(IW,9000)
     *      'DIHEDRAL ROTATION             ENERGY = ', EN123R4*TOKCAL
            WRITE(IW,9000)
     *      'DIHEDRAL BENDING              ENERGY = ', EN123B4*TOKCAL
            WRITE(IW,9000)
     *      'CMAP                          ENERGY = ', ENCMAP *TOKCAL
            WRITE(IW,9000)
     *      'WAGGING                       ENERGY = ', EN234W1*TOKCAL
            IF(NUMBTYP.GT.0)
     *      WRITE(IW,9000)
     *      'UMBRELLA SAMPLING BIAS        ENERGY = ', ENBIAS *TOKCAL
            WRITE(IW,9000)
     *      'LJ REPULSION                  ENERGY = ', ENLJR  *TOKCAL
            WRITE(IW,9000)
     *      'LJ DISPERSION                 ENERGY = ', ENLJD  *TOKCAL
            WRITE(IW,9000)
     *      'CHARGE                        ENERGY = ', ENCHAR *TOKCAL
            WRITE(IW,9000)
     *      'INDUCED DIPOLE                ENERGY = ', ENPOL  *TOKCAL
            WRITE(IW,9000)
     *      'SPHSOL                        ENERGY = ', ENRXN  *TOKCAL
            WRITE(IW,9000)
     *      'FIXSOL                        ENERGY = ', ENFIXSO*TOKCAL
            WRITE(IW,9000)
     *      'QM CENTER                     ENERGY = ', ENCENT *TOKCAL
            WRITE(IW,9000)
     *      'SPHERE                        ENERGY = ', ENRXNR *TOKCAL
            WRITE(IW,9000)
     *      'POTENTIAL                     ENERGY = ', ENPOT  *TOKCAL
            WRITE(IW,9000)
     *      'TOTAL                         ENERGY = ', ENTOT  *TOKCAL
            IF(VOL.LE.1.0D+30)
     *      WRITE(IW,9003)
     *      'VOLUME                               = ',  VOL*TOANGS**3
            IF(VOL.GT.1.0D+30)
     *      WRITE(IW,9004)
     *      'VOLUME                               = '
            WRITE(IW,*)' '
            CALL TIMIT(1)
            WRITE(IW,*)' '
         END IF
C
      ELSE IF(NAT.GT.0.AND.NFFAT.GT.0) THEN
                              ENPOT = ETOT
         IF(MPLEVL.EQ.2)      ENPOT = EMP2
         IF(TDDFTYP.NE.RNONE) ENPOT = ESTATE(NTHST)
         ENTOT = ENPOT
         ENRXN = ENRXN + ENRXNPOL
         IF(MASWRK.AND.(MOD(ISTEP,JOUT).EQ.0.OR.ISTEP.EQ.NSTEP.OR.
     *      ICONV.EQ.1)) THEN
            WRITE(IW,*)' '
            WRITE(IW,'(1X,A,I10,33X,A,F19.2,A)')'OPT STEP',ISTEP
            WRITE(IW,9000)
     *      'BOND STRETCHING               ENERGY = ', EN12   *TOKCAL
            WRITE(IW,9000)
     *      'BOND ANGLE BENDING            ENERGY = ', EN123  *TOKCAL
            WRITE(IW,9000)
     *      'STRETCHING BENDING            ENERGY = ', EN12312*TOKCAL
            WRITE(IW,9000)
     *      'DIHEDRAL ROTATION             ENERGY = ', EN123R4*TOKCAL
            WRITE(IW,9000)
     *      'DIHEDRAL BENDING              ENERGY = ', EN123B4*TOKCAL
            WRITE(IW,9000)
     *      'CMAP                          ENERGY = ', ENCMAP *TOKCAL
            WRITE(IW,9000)
     *      'WAGGING                       ENERGY = ', EN234W1*TOKCAL
            IF(NUMBTYP.GT.0)
     *      WRITE(IW,9000)
     *      'UMBRELLA SAMPLING BIAS        ENERGY = ', ENBIAS *TOKCAL
            WRITE(IW,9000)
     *      'LJ REPULSION                  ENERGY = ', ENLJR  *TOKCAL
            WRITE(IW,9000)
     *      'LJ DISPERSION                 ENERGY = ', ENLJD  *TOKCAL
            WRITE(IW,9000)
     *      'CHARGE                        ENERGY = ', ENCHAR *TOKCAL
            WRITE(IW,9000)
     *      'INDUCED DIPOLE                ENERGY = ', ENPOL  *TOKCAL
            WRITE(IW,9000)
     *      'SPHSOL                        ENERGY = ', ENRXN  *TOKCAL
            WRITE(IW,9000)
     *      'FIXSOL                        ENERGY = ', ENFIXSO*TOKCAL
            WRITE(IW,9000)
     *      'QM CENTER                     ENERGY = ', ENCENT *TOKCAL
            WRITE(IW,9000)
     *      'SPHERE                        ENERGY = ', ENRXNR *TOKCAL
            WRITE(IW,9000)
     *      'POTENTIAL                     ENERGY = ', ENPOT  *TOKCAL
            WRITE(IW,9000)
     *      'TOTAL                         ENERGY = ', ENTOT  *TOKCAL
            IF(VOL.LE.1.0D+30)
     *      WRITE(IW,9003)
     *      'VOLUME                               = ',  VOL*TOANGS**3
            IF(VOL.GT.1.0D+30)
     *      WRITE(IW,9004)
     *      'VOLUME                               = '
            WRITE(IW,*)' '
            CALL TIMIT(1)
            WRITE(IW,*)' '
         END IF
      END IF
C
C     -- CALCULATE RMSD EVERY JOUT STEPS --
C
      IF(NRMSD.EQ.1) THEN
         IF(MOD(ISTEP,JOUT).EQ.0) THEN
            CALL RMSD(X(LFFCORD),ISTEP,X(LFFZANF),X(LFFRMSD0))
         END IF
      END IF
C
C     -- CALCULATE RADIUS OF GYRATION --
C
      IF(NGYRA.GT.0) THEN
         IF(MOD(ISTEP,JOUT).EQ.0) THEN
            CALL GYRA(X(LFFCORD),X(LFFZANF),X(LFFZMAS))
         END IF
      END IF
C
C     -- CALCULATE RALL --
C
      IF(NRALL.EQ.1) THEN
         IF(MOD(ISTEP,JOUT).EQ.0) THEN
            CALL RALL(X(LFFRALL0),X(LFFCORD),ISTEP)
         END IF
      END IF
C
C     -- CALCULATE SELECT DISTANCES EVERY JOUT STEPS --
C
      IF(((NRIJMM+NRIJQM+NAIJKMM+NAIJKQM).GT.0) .AND.
     *    MOD(ISTEP,JOUT).EQ.0) CALL DISIJ(X(LFFCORD))
C
 9000 FORMAT(1X,A,F30.10,2X,'KCAL/MOL')
 9003 FORMAT(1X,A,F30.10,2X,'A**3')
 9004 FORMAT(1X,A,30X,   2X,'OPEN SYSTEM')
C
      CALL FLSHBF(IW)
      RETURN
      END
C*MODULE QUANPOC  *DECK SHIFT
!>
!> @brief    shifting function
!>
!> @author   Hui Li
!>           - Feb 2012
!>
!> @details  several shifting functions
!>
      SUBROUTINE SHIFT(R2,R,ONER,CX,CY,CZ)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      PARAMETER (ZERO=0.0D+00)
      PARAMETER (ONE=1.0D+00)
C
      COMMON /FFPBSW/ XBOX,YBOX,ZBOX,SWF,SWFDX,SWFDY,SWFDZ,
     *                SWRA,ONESWRA,SWRA2,ONESWRA2,
     *                SWRB,ONESWRB,SWRB2,ONESWRB2,ONESWRB4,
     *                SWFDUM3,SWFDUM4,SWFDUM5,
     *                SWRAQ,ONESWRAQ,SWRAQ2,ONESWRAQ2,
     *                SWRBQ,ONESWRBQ,SWRBQ2,ONESWRBQ2,ONESWRBQ4,
     *                SWFDUM3Q,SWFDUM4Q,SWFDUM5Q,
     *                QMSIZE,QMCX,QMCY,QMCZ,QMCXSV,QMCYSV,QMCZSV,
     *                CENTX,CENTY,CENTZ,BUFWID1,BUFWID2,RDAMP,
     *                EFIELDX,EFIELDY,EFIELDZ,QMCXSV2,QMCYSV2,QMCZSV2,
     *                EPS1RB,EPS1RB3,ONEXBOX,ONEYBOX,ONEZBOX,
     *                LQMCT,MXLIST1,NTODO,NTODOSV,NTODOQ,
     *                ISWITCH,ISHIFT,IPOLSHF,
     *                LFFLSTCELL,LFFCORDSV,
     *                LFFPOLSV,LFFCORDSV2,LFFNONLS2,LFFCORDSVQ,
     *                LFFMVFASTS2,LFFMVFASTS3,LFFMVFASTS4,
     *                LFFMVFASTL2,LFFMVFASTL3,LFFMVFASTL4,
     *                MXCHECK,MXLIST2,NTODO2,NTODO2SV
C
C     HUI LI, FEB 27, 2012
C
      IF(ISHIFT.EQ.0) THEN
         SWF   = ONE
         SWFDX = ZERO
         SWFDY = ZERO
         SWFDZ = ZERO
         RETURN
      END IF
C
      IF(R2.LE.SWRB2) THEN
         IF(     ISHIFT.EQ.1) THEN
            SWF   = ONE - R*ONESWRB
            SWF   = SWF*SWF
            DSWF  = ONESWRB2 - ONESWRB*ONER
            DSWF  = DSWF + DSWF
            SWFDX = DSWF*CX
            SWFDY = DSWF*CY
            SWFDZ = DSWF*CZ
         ELSE IF(ISHIFT.EQ.2) THEN
            SWF   = ONE - R*EPS1RB + R2*R*EPS1RB3
            DSWF  = R*EPS1RB3
            DSWF  = DSWF + DSWF + DSWF
            DSWF  = DSWF - ONER*EPS1RB
            SWFDX = DSWF*CX
            SWFDY = DSWF*CY
            SWFDZ = DSWF*CZ
         ELSE IF(ISHIFT.EQ.3) THEN
            SWF   = ONE - R*ONESWRB
            DSWF  = -ONER*ONESWRB
            SWFDX = DSWF*CX
            SWFDY = DSWF*CY
            SWFDZ = DSWF*CZ
         ELSE IF(ISHIFT.EQ.4) THEN
            SWF   = ONE - R2*ONESWRB2
            SWF   = SWF*SWF
            DSWF  = R2*ONESWRB4 - ONESWRB2
            DSWF  = DSWF + DSWF
            DSWF  = DSWF + DSWF
            SWFDX = DSWF*CX
            SWFDY = DSWF*CY
            SWFDZ = DSWF*CZ
         ELSE
            SWF   = ONE
            SWFDX = ZERO
            SWFDY = ZERO
            SWFDZ = ZERO
         END IF
      ELSE
         SWF   = ZERO
         SWFDX = ZERO
         SWFDY = ZERO
         SWFDZ = ZERO
      END IF
C
      RETURN
      END
C*MODULE QUANPOC  *DECK SWFUNC
!>
!> @brief    switching function for MM atoms
!>
      SUBROUTINE SWFUNC(R2,CX,CY,CZ)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      PARAMETER (ZERO=0.0D+00)
      PARAMETER (ONE=1.0D+00)
C
      COMMON /FFPBSW/ XBOX,YBOX,ZBOX,SWF,SWFDX,SWFDY,SWFDZ,
     *                SWRA,ONESWRA,SWRA2,ONESWRA2,
     *                SWRB,ONESWRB,SWRB2,ONESWRB2,ONESWRB4,
     *                SWFDUM3,SWFDUM4,SWFDUM5,
     *                SWRAQ,ONESWRAQ,SWRAQ2,ONESWRAQ2,
     *                SWRBQ,ONESWRBQ,SWRBQ2,ONESWRBQ2,ONESWRBQ4,
     *                SWFDUM3Q,SWFDUM4Q,SWFDUM5Q,
     *                QMSIZE,QMCX,QMCY,QMCZ,QMCXSV,QMCYSV,QMCZSV,
     *                CENTX,CENTY,CENTZ,BUFWID1,BUFWID2,RDAMP,
     *                EFIELDX,EFIELDY,EFIELDZ,QMCXSV2,QMCYSV2,QMCZSV2,
     *                EPS1RB,EPS1RB3,ONEXBOX,ONEYBOX,ONEZBOX,
     *                LQMCT,MXLIST1,NTODO,NTODOSV,NTODOQ,
     *                ISWITCH,ISHIFT,IPOLSHF,
     *                LFFLSTCELL,LFFCORDSV,
     *                LFFPOLSV,LFFCORDSV2,LFFNONLS2,LFFCORDSVQ,
     *                LFFMVFASTS2,LFFMVFASTS3,LFFMVFASTS4,
     *                LFFMVFASTL2,LFFMVFASTL3,LFFMVFASTL4,
     *                MXCHECK,MXLIST2,NTODO2,NTODO2SV
C
C     HUI LI, JAN 2011, LINCOLN
C     HUI LI, FEB 2012
C
      IF(ISWITCH.EQ.0) THEN
         SWF   = ONE
         SWFDX = ZERO
         SWFDY = ZERO
         SWFDZ = ZERO
         RETURN
      END IF
C
      IF(R2.LE.SWRA2) THEN
         SWF   = ONE
         SWFDX = ZERO
         SWFDY = ZERO
         SWFDZ = ZERO
      ELSE IF(R2.LE.SWRB2) THEN
         VUMY1 = R2 - SWRA2
         VUMY2 = VUMY1*VUMY1
         VUMY3 = VUMY2*VUMY1
         VUMY4 = VUMY2*VUMY2
         VUMY5 = VUMY2*VUMY3
         SWF   = ONE - 10.0D+00*SWFDUM3*VUMY3
     *               + 15.0D+00*SWFDUM4*VUMY4
     *               -  6.0D+00*SWFDUM5*VUMY5
         DSWF  = -30.0D+00*SWFDUM3*VUMY2
     *           +60.0D+00*SWFDUM4*VUMY3
     *           -30.0D+00*SWFDUM5*VUMY4
         DSWF  = DSWF + DSWF
         SWFDX = DSWF*CX
         SWFDY = DSWF*CY
         SWFDZ = DSWF*CZ
      ELSE
         SWF   = ZERO
         SWFDX = ZERO
         SWFDY = ZERO
         SWFDZ = ZERO
      END IF
C
      RETURN
      END
C*MODULE QUANPOC  *DECK SWFUNCQ
!>
!> @brief    switching function for QM atoms
!>
      SUBROUTINE SWFUNCQ(R2,CX,CY,CZ)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      PARAMETER (ZERO=0.0D+00)
      PARAMETER (ONE=1.0D+00)
C
      COMMON /FFPBSW/ XBOX,YBOX,ZBOX,SWF,SWFDX,SWFDY,SWFDZ,
     *                SWRA,ONESWRA,SWRA2,ONESWRA2,
     *                SWRB,ONESWRB,SWRB2,ONESWRB2,ONESWRB4,
     *                SWFDUM3,SWFDUM4,SWFDUM5,
     *                SWRAQ,ONESWRAQ,SWRAQ2,ONESWRAQ2,
     *                SWRBQ,ONESWRBQ,SWRBQ2,ONESWRBQ2,ONESWRBQ4,
     *                SWFDUM3Q,SWFDUM4Q,SWFDUM5Q,
     *                QMSIZE,QMCX,QMCY,QMCZ,QMCXSV,QMCYSV,QMCZSV,
     *                CENTX,CENTY,CENTZ,BUFWID1,BUFWID2,RDAMP,
     *                EFIELDX,EFIELDY,EFIELDZ,QMCXSV2,QMCYSV2,QMCZSV2,
     *                EPS1RB,EPS1RB3,ONEXBOX,ONEYBOX,ONEZBOX,
     *                LQMCT,MXLIST1,NTODO,NTODOSV,NTODOQ,
     *                ISWITCH,ISHIFT,IPOLSHF,
     *                LFFLSTCELL,LFFCORDSV,
     *                LFFPOLSV,LFFCORDSV2,LFFNONLS2,LFFCORDSVQ,
     *                LFFMVFASTS2,LFFMVFASTS3,LFFMVFASTS4,
     *                LFFMVFASTL2,LFFMVFASTL3,LFFMVFASTL4,
     *                MXCHECK,MXLIST2,NTODO2,NTODO2SV
C
C     HUI LI, JAN 2011, LINCOLN
C     HUI LI, FEB 2012
C
      IF(ISWITCH.EQ.0) THEN
         SWF   = ONE
         SWFDX = ZERO
         SWFDY = ZERO
         SWFDZ = ZERO
         RETURN
      END IF
C
      IF(R2.LE.SWRAQ2) THEN
         SWF   = ONE
         SWFDX = ZERO
         SWFDY = ZERO
         SWFDZ = ZERO
      ELSE IF(R2.LE.SWRBQ2) THEN
         VUMY1 = R2 - SWRAQ2
         VUMY2 = VUMY1*VUMY1
         VUMY3 = VUMY2*VUMY1
         VUMY4 = VUMY2*VUMY2
         VUMY5 = VUMY2*VUMY3
         DSWF  = -30.0D+00*SWFDUM3Q*VUMY2
     *           +60.0D+00*SWFDUM4Q*VUMY3
     *           -30.0D+00*SWFDUM5Q*VUMY4
         DSWF  = DSWF + DSWF
         SWF   = ONE - 10.0D+00*SWFDUM3Q*VUMY3
     *               + 15.0D+00*SWFDUM4Q*VUMY4
     *               -  6.0D+00*SWFDUM5Q*VUMY5
         SWFDX = DSWF*CX
         SWFDY = DSWF*CY
         SWFDZ = DSWF*CZ
      ELSE
         SWF   = ZERO
         SWFDX = ZERO
         SWFDY = ZERO
         SWFDZ = ZERO
      END IF
C
      RETURN
      END
C*MODULE QUANPOC  *DECK QMMM1EINT
!>
!> @brief    1-e integrals of MM charges
!>
!> @author   Hui Li
!>           - Apr 2011
!>
!> @details  1-e integrals of MM charges
!>           MM force field energies
!>
      SUBROUTINE QMMM1EINT(WRK,H,L2,CORD,CHARG,
     *                     NONLSTQ,
     *                     CLPR,ZLPR,NLPR,POL,POLSV,
     *                     FIELD1,NONLS1,
     *                     XTS,YTS,ZTS,POT1,
     *                     FFGRD2,BOND0,FCBOND,
     *                     IPAIR,FCSTBD,KBLST,
     *                     ANGL0,FCANGL,KLIST,
     *                     DIHB0,FCDIHB,FCDIHR,NLIST,
     *                     FCWAGG,MLIST,
     *                     VROT,GAMA,NNN,LLIST,
     *                     MAPLST,CMAPCO,
     *                     SIG,EPS,
     *                     SIG2,EPS2,L14J,L1213J,
     *                     DETMP,LISTQM,FCLJTP,NTYPE,
     *                     NONLSA,NONLSB,L1213A,
     *                     L1213B,L14A,L14B,CHARGB,
     *                     SIGB,EPSB,SIG2B,EPS2B,
     *                     NONLSPMA,
     *                     L1213PMA,L14PMA,CORDB,
     *                     LSBONDPMA,LSANGLPMA,
     *                     LSDIHRPMA,LSDIHBPMA,
     *                     LSWAGGPMA,LSCMAPPMA,
     *                     NONLSPMB,L1213PMB,L14PMB,
     *                     LSBONDPMB,LSANGLPMB,
     *                     LSDIHRPMB,LSDIHBPMB,
     *                     LSWAGGPMB,LSCMAPPMB,UMBHIS,UM2HIS)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      LOGICAL GOPARR,DSKWRK,MASWRK
      LOGICAL DFTBFL,SCC,SRSCC,DFTB3,DAMPXH,LRDFTB
C
      DIMENSION ENALL(100)
C
      COMMON /DFTB  / DFTBFL,SCC,SRSCC,DFTB3,DAMPXH,LRDFTB
      COMMON /FFENGY/ EN12,EN123,EN123R4,EN123B4,EN234W1,ENCHAR,ENLJR,
     *                ENLJD,ENPOL,XENPOL,ENRXN,XENRXN,ENRXNPOL,ENRXNR,
     *                EN12312,ENQUANP(30),
     *                ENBIAS,ENCENT,ENUCCH,ENCMAP,ENPOT,ENKIN,ENTOT
      COMMON /FFFREE/ SOL1CH,SOL1LJ,SOL1IM,SOLFRE1MM,
     *                SOL2CH,SOL2LJ,SOL2IM,SOLFRE2MM,
     *                ASOL1CH,ASOL1LJ,ASOL1IM,ASOL1MM,
     *                ASOL2CH,ASOL2LJ,ASOL2IM,ASOL2MM,
     *                PMF1BD,PMF1AG,PMF1DR,PMF1DB,PMF1WG,PMF1CM,
     *                PMF1CH,PMF1PO,PMF1LJ,PMF1IM,PMFFRE1MM,
     *                EN12A,EN12B,EN123A,EN123B,EN123R4A,EN123R4B,
     *                EN234W1A,EN234W1B,EN123B4A,EN123B4B,
     *                EN12312A,EN12312B,
     *                ENCHARA,ENCHARB,ENPOLA,ENPOLB,
     *                ENRXNA,ENRXNB,ENRXNPOLA,ENRXNPOLB,
     *                ENRXNRA,ENRXNRB,ENLJRA,ENLJRB,ENLJDA,ENLJDB,
     *                ENCMAPA,ENCMAPB,ENFIXSOA,ENFIXSOB,
     *                ENCENTA,ENCENTB,WSIMUL,WPERT1,WPERT2,
     *                IFEPTYP
      COMMON /FFPARM/ NFFAT,NBOND,NANGL,NDIHR,NDIHB,NCMAP,NWAGG,
     *                N1213J,N14J,NLKQMM,IDOCHG,IDOPOL,IDOLJ,IDOCMAP
      COMMON /FFRXN / RXNEPS,RSPHSOL,ISPHSOL
      COMMON /IOFILE/ IR,IW,IP,IJK,IJKT,IDAF,NAV,IODA(950)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /WFNOPT/ SCFTYP,VBTYP,DFTYPE,TDDFTYP,CITYP,CCTYP,
     *                MPLEVL,MPCTYP
C
      DATA NONE/4HNONE/
C
C     HUI LI, APR 2011, LINCOLN
C     HUI LI, MAY 2012
C
      CALL VCLR(FFGRD2,1,3*NFFAT)
C
C     QUANPOL CANNOT USE G, H, I TYPE FUNCTIONS
C
      IF(NFFAT.GT.0) THEN
         CALL BASCHK(LMAX)
         IF(LMAX.GE.4) THEN
            IF(MASWRK)WRITE(IW,'(1X,/1X,A,/)')
     *      'ERROR: QUANPOL CANNOT HANDLE G, H, I TYPE FUNCTIONS.'
            IF(MASWRK)WRITE(IW,*)' '
            CALL ABRT
         END IF
      END IF
C
C     - MOPAC -
      IF(MPCTYP.NE.NONE) GOTO 150
C     - DFTB -
      IF(DFTBFL) THEN
         CALL DFTB_QMMMCHGINT(WRK,CORD,CHARG,NONLSTQ)
         GOTO 150
      END IF
C
      CALL VCLR(WRK,1,L2)
      CALL VCLR(H,1,L2)
      CALL DAWRIT(IDAF,IODA,H,L2,89,0)
      CALL QMMMCHGINT(WRK,CORD,CHARG,NONLSTQ)
      CALL QMMMREPINT(WRK,CORD,CLPR,ZLPR,NLPR,NONLSTQ)
      IF (GOPARR) CALL DDI_GSUMF(920,WRK,L2)
C
      CALL DAREAD(IDAF,IODA,H,L2,11,0)
      CALL VADD(H,1,WRK,1,H,1,L2)
      CALL DAWRIT(IDAF,IODA,H,L2,11,0)
C
      CALL DAREAD(IDAF,IODA,H,L2,89,0)
      CALL VADD(H,1,WRK,1,H,1,L2)
      CALL DAWRIT(IDAF,IODA,H,L2,89,0)
C
 150  CONTINUE
C
C     - PREPARE FOR MM POLARIZATION FIELD -
C
      CALL POLSCALE(CORD,POL,POLSV)
      CALL QMMMPOLFLDNUCMM(CORD,CHARG,POL,FIELD1,
     *                     NONLS1,NONLSTQ,L1213J)
      CALL QMMMRXNPOTNUCMM(CORD,CHARG,XTS,YTS,ZTS,POT1,ISPHSOL)
C
C     - FORCE FIELD TERMS -
C         MM GRADIENT CAN BE DONE HERE, BUT NO DDI_GSUMF
C         QM GRADIENT MUST WAIT FOR STVDER
C
      CALL E00012(CORD,FFGRD2,BOND0,FCBOND,IPAIR,CORDB,
     *            LSBONDPMA,LSBONDPMB)
      CALL E00123(CORD,FFGRD2,ANGL0,FCANGL,KLIST,CORDB,
     *            LSANGLPMA,LSANGLPMB)
      CALL E12312(CORD,FFGRD2,ANGL0,KLIST,BOND0,FCSTBD,
     *            KBLST,CORDB,LSANGLPMA,LSANGLPMB)
      CALL E123B4(CORD,FFGRD2,DIHB0,FCDIHB,NLIST,CORDB,
     *            LSDIHBPMA,LSDIHBPMB)
      CALL E234W1(CORD,FFGRD2,      FCWAGG,MLIST,CORDB,
     *            LSWAGGPMA,LSWAGGPMB)
      CALL E123R4(CORD,FFGRD2,VROT,GAMA,NNN,LLIST,CORDB,
     *            LSDIHRPMA,LSDIHRPMB,FCDIHR)
      CALL ECMAP (CORD,FFGRD2,MAPLST,CMAPCO,CORDB,
     *            LSCMAPPMA,LSCMAPPMB)
      CALL ELJ126(CORD,FFGRD2,SIG,EPS,SIG2,EPS2,L14J,NONLS1,
     *            L1213J,SIGB,EPSB,SIG2B,EPS2B,
     *            NONLSA,NONLSB,L1213A,L1213B,L14A,L14B,
     *            NONLSPMA,L1213PMA,L14PMA,CORDB,
     *            NONLSPMB,L1213PMB,L14PMB,FCLJTP,NTYPE)
      CALL ELJ126QM(CORD,FFGRD2,SIG,EPS,
     *              FCLJTP,NTYPE,LISTQM)
      CALL ESPHER(CORD,FFGRD2)
      CALL ECHARG(CORD,FFGRD2,CHARG,NONLS1,L1213J,L14J,
     *            CHARGB,NONLSA,NONLSB,L1213A,L1213B,L14A,L14B,
     *            NONLSPMA,L1213PMA,L14PMA,CORDB,
     *            NONLSPMB,L1213PMB,L14PMB)
      CALL UMBRELLA(CORD,FFGRD2,UMBHIS,UM2HIS)
      IF(MPCTYP.EQ.NONE)
     *CALL QMMMCHGNUC(CORD,FFGRD2,DETMP,CHARG,NONLSTQ)
      IF(MPCTYP.NE.NONE)
     *CALL QMMMCHGCHG(CORD,FFGRD2,DETMP,CHARG,NONLSTQ)
      CALL ESPHQM(DETMP,LISTQM)
C
      IF(GOPARR) THEN
         CALL VCLR(ENALL,1,100)
         ENALL( 1) = EN12
         ENALL( 2) = EN123
         ENALL( 3) = EN123R4
         ENALL( 4) = EN123B4
         ENALL( 5) = EN234W1
         ENALL( 6) = ENCHAR
         ENALL( 7) = ENLJR
         ENALL( 8) = ENLJD
C        ENALL( 9) = VIR(1)    ! NOT HERE
C        ENALL(10) = VIR(2)
C        ENALL(11) = VIR(3)
         ENALL(12) = ENUCCH    ! ONLY HERE
         ENALL(13) = ENRXN
         ENALL(14) = ENRXNR
         ENALL(15) = ENCENT
         ENALL(16) = ENCMAP
         ENALL(18) = SOL1CH
         ENALL(19) = SOL1LJ
         ENALL(21) = SOL2CH
         ENALL(22) = SOL2LJ
         ENALL(33) = PMF1BD
         ENALL(34) = PMF1AG
         ENALL(35) = PMF1DR
         ENALL(36) = PMF1DB
         ENALL(37) = PMF1WG
         ENALL(38) = PMF1CM
         ENALL(39) = PMF1CH
         ENALL(40) = PMF1LJ
         ENALL(41) = ENBIAS
         ENALL(42) = EN12312
         CALL DDI_GSUMF(2410,ENALL  ,42)
C        CALL DDI_GSUMF(2411,FFGRD2,3*NFFAT)   ! NOT HERE
         EN12      = ENALL( 1)
         EN123     = ENALL( 2)
         EN123R4   = ENALL( 3)
         EN123B4   = ENALL( 4)
         EN234W1   = ENALL( 5)
         ENCHAR    = ENALL( 6)
         ENLJR     = ENALL( 7)
         ENLJD     = ENALL( 8)
C        VIR(1)    = ENALL( 9)
C        VIR(2)    = ENALL(10)
C        VIR(3)    = ENALL(11)
         ENUCCH    = ENALL(12)   ! ONLY HERE
         ENRXN     = ENALL(13)
         ENRXNR    = ENALL(14)
         ENCENT    = ENALL(15)
         ENCMAP    = ENALL(16)
         SOL1CH    = ENALL(18)
         SOL1LJ    = ENALL(19)
         SOL2CH    = ENALL(21)
         SOL2LJ    = ENALL(22)
         PMF1BD    = ENALL(33)
         PMF1AG    = ENALL(34)
         PMF1DR    = ENALL(35)
         PMF1DB    = ENALL(36)
         PMF1WG    = ENALL(37)
         PMF1CM    = ENALL(38)
         PMF1CH    = ENALL(39)
         PMF1LJ    = ENALL(40)
         ENBIAS    = ENALL(41)
         EN12312   = ENALL(42)
      END IF
      RETURN
C
      END
C*MODULE QUANPOC  *DECK QMMMGRD
!>
!> @brief    QM/MM gradient
!>
!> @author   Hui Li, Dejun Si
!>           - Mar 2011
!>
!> @details  calculate forces on QM and MM atoms
!>
      SUBROUTINE QMMMGRD(DTOT)
      use mx_limits, only: mxatm,mxsh,mxgtot
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      LOGICAL TRIPLET,SG1T,TAMMD,TPA,ALPHKWD,BETAKWD
      LOGICAL MREKT,MRDEA
      LOGICAL DFTBFL,SCC,SRSCC,DFTB3,DAMPXH,LRDFTB
      LOGICAL GOPARR,MASWRK,DSKWRK
C
      PARAMETER (MXSPE=10)
      PARAMETER (ONE=1.0D+00)
C
      DIMENSION DTOT(*)
C
      COMMON /DFTB  / DFTBFL,SCC,SRSCC,DFTB3,DAMPXH,LRDFTB
      COMMON /DFTBPR/ ETEMP,DFTBDP(MXSPE*14),DAMPXHE,HUBDER(MXSPE),
     *                ZREF(MXATM),SPNCST(6,MXSPE),SPE(MXATM),NSPE,
     *                MAXANG(MXATM),ISPE(MXATM),IND(MXATM+1),IDFTBD,
     *                PARAMDIR
      COMMON /FFFIXS/ ENFIXSO,FIXEPS,FIXTOL,FIXA,FIXQ,RALLMM,RALLQM,
     *                RADMM(200),RADQM(200),NRADMM,NRADQM,IFIXSOL,
     *                LFFDAI,LFFDAIT,LFFIDDAI,LFFIDTMP,LFFTMPTS,
     *                LFFAFIX,LFFIDATOM,LFFRFIX,LFFQFIX,NTSATM,
     *                LFFQFIXMP,LFFQFIXTA,LFFQFIXXY,
     *                LFFXTSFIX,LFFYTSFIX,LFFZTSFIX,
     *                LFFVFIX1,LFFVFIX2,NCYCLE,MXFFTS,NFFTS
      COMMON /FFMDPA/ DT,DT2,TEMP0,PRES0,POLTOL,VIR(3),PMEAN,VOLAV,
     *                ENPAV,ENKAV,TEMPAV,BERENDT,BERENDP,VELMAX,
     *                PMEANX,PMEANY,PMEANZ,NSTEP,KMASTER,KOUTACT(2),
     *                IHESS,INTALG,ITSTAT,IPSTAT,JOUT,KOUT,LOUT
      COMMON /FFMPNT/ LFFATMNAM,LFFCORD,LFFZANF,
     *                LFFZMAS,LFFONEMAS,LFFQMZMAS,LFFQM1MAS,
     *                LFFCHARG,LFFPOL,LFFDIP,
     *                LFFFIELD1,LFFFIELD2,LFFFIELD3,
     *                LFFSIG,LFFEPS,LFFSIG2,LFFEPS2,
     *                LFFBOND0,LFFFCBOND,
     *                LFFANGL0,LFFFCANGL,LFFFCWAGG,
     *                LFFDIHB0,LFFFCDIHB,
     *                LFFVROT,LFFNNN,LFFGAMA,LFFIPAIR,
     *                LFFKLIST,LFFLLIST,LFFL1213J,LFFL14J,
     *                LFFMLIST,LFFNLIST,LFFLKQMMM,
     *                LFFVEL,LFFQMVEL,
     *                LFFFFGRD0,LFFFFGRD1,LFFFFGRD2,
     *                LFFQMGRD0,LFFQMGRD1,LFFQMGRD2,LFFDETMP,
     *                LFFCLPR,LFFZLPR,LFFNLPR,
     *                LFFXTS,LFFYTS,LFFZTS,LFFCMAT1,
     *                LFFQRXN1,LFFQRXN2,LFFPOT1,LFFPOT2,LFFQRXNMP,
     *                LFFQRXNTA,LFFQRXNXY,LFFNONLSTQ,
     *                LFFDIPMP,LFFDIPTA,LFFDIPXY,LFFLISTQM,LFFNONLS1,
     *                LFFMAPLST,LFFCMAPCO
      COMMON /FFPARM/ NFFAT,NBOND,NANGL,NDIHR,NDIHB,NCMAP,NWAGG,
     *                N1213J,N14J,NLKQMM,IDOCHG,IDOPOL,IDOLJ,IDOCMAP
      COMMON /FFPBSW/ XBOX,YBOX,ZBOX,SWF,SWFDX,SWFDY,SWFDZ,
     *                SWRA,ONESWRA,SWRA2,ONESWRA2,
     *                SWRB,ONESWRB,SWRB2,ONESWRB2,ONESWRB4,
     *                SWFDUM3,SWFDUM4,SWFDUM5,
     *                SWRAQ,ONESWRAQ,SWRAQ2,ONESWRAQ2,
     *                SWRBQ,ONESWRBQ,SWRBQ2,ONESWRBQ2,ONESWRBQ4,
     *                SWFDUM3Q,SWFDUM4Q,SWFDUM5Q,
     *                QMSIZE,QMCX,QMCY,QMCZ,QMCXSV,QMCYSV,QMCZSV,
     *                CENTX,CENTY,CENTZ,BUFWID1,BUFWID2,RDAMP,
     *                EFIELDX,EFIELDY,EFIELDZ,QMCXSV2,QMCYSV2,QMCZSV2,
     *                EPS1RB,EPS1RB3,ONEXBOX,ONEYBOX,ONEZBOX,
     *                LQMCT,MXLIST1,NTODO,NTODOSV,NTODOQ,
     *                ISWITCH,ISHIFT,IPOLSHF,
     *                LFFLSTCELL,LFFCORDSV,
     *                LFFPOLSV,LFFCORDSV2,LFFNONLS2,LFFCORDSVQ,
     *                LFFMVFASTS2,LFFMVFASTS3,LFFMVFASTS4,
     *                LFFMVFASTL2,LFFMVFASTL3,LFFMVFASTL4,
     *                MXCHECK,MXLIST2,NTODO2,NTODO2SV
      COMMON /FFRXN / RXNEPS,RSPHSOL,ISPHSOL
      COMMON /FMCOM / X(1)
      COMMON /FUNCT / E,EG(3,MXATM)
      COMMON /GRAD  / DE(3,MXATM)
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /INFOTD/ CNVTOL,PFREQ(2),MODTD,
     *                JANST,NRADT,NTHET,NPHIT,NLEBT,
     *                NSTAT,NTRIAL,MAXVEC,NTHST,IRECTD,ITDFG,ITDPRP,
     *                TRIPLET,SG1T,NONEQR,TAMMD,TPA,ALPHKWD,BETAKWD,
     *                SPCP(3),MULTD,MREKT,MRDEA,MTHST,IFEDAT(4)
      COMMON /IOFILE/ IR,IW,IP,IJK,IJKT,IDAF,NAV,IODA(950)
      COMMON /NSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),
     *                CF(MXGTOT),CG(MXGTOT),CH(MXGTOT),CI(MXGTOT),
     *                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),
     *                KLOC(MXSH),MIN(MXSH),MAX(MXSH),NSHELL
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /RUNOPT/ RUNTYP,EXETYP,NEVALS,NGLEVL,NHLEVL
      COMMON /WFNOPT/ SCFTYP,VBTYP,DFTYPE,TDDFTYP,CITYP,CCTYP,
     *                MPLEVL,MPCTYP
C
      DATA CHECK/8HCHECK   /
      DATA RHF,UHF,ROHF/8HRHF     ,8HUHF     ,8HROHF    /
      DATA NONE,RNONE/4HNONE,8HNONE    /
C
C     HUI LI, DEJUN SI, MAR 2011, LINCOLN
C     HUI LI, DEC 22, 2011
C     HUI LI, NOV 23, 2016, ADD DFTB
C     HUI LI, JUL 10, 2020, DETECT UNCONVERGED SCF
C
C     -- ALL GRADIENTS SHOULD BE ZERO IF SCF WAS NOT CONVERGED --
      IF(ABS(E).LT.1.0D-12) THEN
         CALL VCLR(DE,1,3*NAT)
         CALL VCLR(EG,1,3*NAT)
         CALL VCLR(X(LFFQMGRD2),1,3*NAT)
         CALL VCLR(X(LFFFFGRD2),1,3*NFFAT)
         RETURN
      END IF
C
C     -- MOPAC --
      IF(MPCTYP.NE.NONE) THEN
C        ----- GRADIENT WITH RESPECT TO QM COORDINATES -----
C              FORCES ON QM CHARGE BY MM CHARGE ALREADY DONE IN DETMP
         CALL DDI_GSUMF(2417,X(LFFDETMP),3*NAT)
         CALL VADD(X(LFFDETMP),1,EG,1,EG,1,3*NAT)
C        ----- GRADIENT WITH RESPECT TO MM COORDINATES -----
C              SOME FORCES ON MM ATOMS HAVE BEEN DONE IN ROUTINE QMMM1EINT
C        ----- GLOBAL SUM OF FORCES ON MM ATOMS ----
C              -VIR- AND -FFGRD2- WERE ZEROED BEFORE CALLING GRADX
         CALL DDI_GSUMF(2418,VIR   ,3)
         CALL DDI_GSUMF(2419,X(LFFFFGRD2),3*NFFAT)
         RETURN
      END IF
C
C     ----- GRADIENT FOR QMMM RUNS -----
C           NOTE: FOR MP2,   DTOT = DHF  + P(2)
C                 FOR TDDFT, DTOT = DDFT + DTA
C
      L1 = NUM
      L2 = (NUM*NUM+NUM)/2
      L3 = NUM*NUM
      CALL VALFM(LOADFM)
      LCH    = LOADFM   + 1
      LCH2   = LCH      + 1200*10  ! GIVE 10 TIMES TO OVER KILL
      LCH3   = LCH2     + 300*10
      LDHFA  = LCH3     + 300*10
      LDHFB  = LDHFA    + L2
      LDXY   = LDHFB    + L2
      LDTEMP = LDXY     + L2
      LAST   = LDTEMP   + L3
      NEED   = LAST     - LOADFM -1
      CALL GETFM(NEED)
C     - GET -DHF-
      IF(MPLEVL.EQ.2) THEN
         IF(SCFTYP.EQ.RHF .OR. SCFTYP.EQ.ROHF) THEN
            CALL DAREAD(IDAF,IODA,X(LDHFA),L2,308,0)
         ELSE IF (SCFTYP.EQ.UHF) THEN
            CALL DAREAD(IDAF,IODA,X(LDHFA),L2,418,0)
            CALL DAREAD(IDAF,IODA,X(LDHFB),L2,428,0)
            CALL VADD(X(LDHFA),1,X(LDHFB),1,X(LDHFA),1,L2)
         END IF
      END IF
      IF(TDDFTYP.NE.RNONE) THEN
         IF(SCFTYP.EQ.RHF) THEN
            CALL DAREAD(IDAF,IODA,X(LDHFA),L2,308,0)
            CALL DAREAD(IDAF,IODA,X(LDTEMP),L3,IRECTD+2,0)
            CALL VCLR(X(LDXY),1,L2)
            CALL TDPCMDEN(X(LDXY),X(LDTEMP),L1)
            CALL DSCAL(L2,2.0D+00,X(LDXY),1)
         END IF
      END IF
C
      IF(EXETYP.EQ.CHECK) GOTO 100
C
      IF(DFTBFL) GOTO 50
C
C     ----- GRADIENT WITH RESPECT TO QM COORDINATES -----
C           FORCES ON QM ELECTRONS BY MM CHARGE, REP, DIP, ASC
C           FORCES ON QM NUCLEI BY MM CHARGE ALREADY DONE IN DETMP
C
      CALL QMMMREPFQM(DTOT,X(LFFCORD),X(LFFCLPR),
     *                X(LFFZLPR),X(LFFNLPR),X(LFFNONLSTQ))
      CALL QMMMCHGFQM(DTOT,X(LFFCORD),X(LFFCHARG),X(LFFNONLSTQ))
      CALL QMMMPOLFQM(DTOT,X(LDHFA),X(LDXY),X(LFFDIP),
     *                X(LFFDIPMP),X(LFFDIPTA),X(LFFDIPXY),
     *                X(LFFCORD),X(LFFNONLSTQ))
      IF(ISPHSOL.GE.60 .AND. RSPHSOL.LT.1.0D+30) THEN
         CALL QMMMRXNFQM(DTOT,X(LDHFA),X(LDXY),X(LFFQRXN1),
     *                   X(LFFQRXN2),X(LFFQRXNMP),X(LFFQRXNTA),
     *                   X(LFFQRXNXY),ISPHSOL)
      END IF
      IF(IFIXSOL.EQ.1) THEN
         CALL QMMMRXNFQM(DTOT,X(LDHFA),X(LDXY),X(LFFQFIX),
     *                   X(LFFQRXN2),X(LFFQFIXMP),X(LFFQFIXTA),
     *                   X(LFFQFIXXY),NFFTS)
      END IF
      CALL VADD(X(LFFDETMP),1,DE,1,DE,1,3*NAT)
C
C     ----- GRADIENT WITH RESPECT TO MM COORDINATES -----
C           FORCES ON MM ATOMS(CHARGE, REP, DIP, ASC) BY QM ELECTRONS
C           SOME FORCES ON MM ATOMS HAVE BEEN DONE IN ROUTINE QMMM1EINT
C
      CALL QMMMCHGFMM(DTOT,X(LCH),X(LFFCORD),X(LFFCHARG),
     *                X(LFFFFGRD2),X(LFFNONLSTQ))
      CALL QMMMREPFMM(DTOT,X(LCH),X(LFFCORD),X(LFFCLPR),
     *                X(LFFZLPR),X(LFFNLPR),
     *                X(LFFFFGRD2),X(LFFNONLSTQ))
      CALL QMMMPOLFMM(DTOT,X(LDHFA),X(LDXY),X(LCH),X(LCH2),X(LCH3),
     *                X(LFFCORD),X(LFFFFGRD2),X(LFFNONLSTQ))
      CALL QMMMRXNFMM(DTOT,X(LDHFA),X(LDXY),X(LCH),
     *                X(LFFQFIX),X(LFFQFIXMP),X(LFFQFIXTA),
     *                X(LFFQFIXXY),X(LFFFFGRD2),X(LFFIDATOM),
     *                X(LFFXTSFIX),X(LFFYTSFIX),X(LFFZTSFIX))
C
 50   CONTINUE
C
      IF(DFTBFL) THEN
C        -- GRADIENT WITH RESPECT TO QM AND MM COORDINATES --
C           FORCES ON QM ELECTRONS BY MM CHARGE, DIP, ASC
C           FORCES ON MM ATOMS (CHARGE, DIP, ASC) BY QM ELECTRONS
C           IT SEEMS DFTB*.SRC USE ONLY -EG-, NOT -DE-
C           BUT -DE- IS USED IN THIS SUBROUTINE FOR POL AND ASC FORCES
         CALL DCOPY(NAT,ZAN,1,X(LDTEMP),1)
         CALL DCOPY(NAT,ZREF,1,ZAN,1)
         CALL VALFM(LOADFM)
         LQMQTOT  = LOADFM   + 1
         LQMQHF   = LQMQTOT  + NAT
         LQMQXY   = LQMQHF   + NAT
         LMUL1    = LQMQXY   + NAT
         LMUL1S   = LMUL1    + L1
         LS       = LMUL1S   + NSHELL
         LDA      = LS       + L2
         LDB      = LDA      + L2
         LAST     = LDB      + L2
         NEEDX    = LAST     - LOADFM -1
         CALL GETFM(NEEDX)
         CALL DENDD1(X(LDA),X(LDB),L2)
         CALL DFTB_QMMMFQMFMM(X(LDA),X(LDHFA),X(LDXY),X(LFFQFIX),
     *                   X(LFFQFIXTA),X(LFFQFIXXY),
     *                   X(LFFDIP),X(LFFDIPTA),X(LFFDIPXY),
     *                   X(LFFCORD),X(LFFFFGRD2),
     *                   X(LFFDETMP),X(LFFCHARG),X(LFFNONLSTQ),
     *                   X(LFFIDATOM),X(LQMQTOT),X(LQMQHF),X(LQMQXY),
     *                   X(LMUL1),X(LMUL1S),X(LS),L1,L2)
         CALL RETFM(NEEDX)
         CALL DDI_GSUMF(2417,X(LFFDETMP),3*NAT)
         CALL VADD(X(LFFDETMP),1,EG,1,EG,1,3*NAT)
         CALL DCOPY(NAT*3,EG,1,DE,1)
         IF(GOPARR) CALL DSCAL(NAT*3,ONE/NPROC,DE,1)
      END IF
C
C     ----- FORCE BETWEEN MM POL AND QM NUC
C                 BETWEEN MM POL AND MM CHARGE
C                 BETWEEN MM POL AND MM POL
C
      CALL QMMMPOLFMMNUCMM(X(LFFCORD),X(LFFFFGRD2),X(LFFCHARG),
     *                     X(LFFPOL),X(LFFPOLSV),X(LFFDIP),
     *                     X(LFFNONLS1),X(LFFNONLSTQ),
     *                     X(LFFDIPMP),X(LFFDIPTA),X(LFFDIPXY),
     *                     X(LFFL1213J))
      CALL POLSCALEFORCE(X(LFFCORD),X(LFFFFGRD2),X(LFFPOL),
     *                   X(LFFDIP))
C
C     ----- FORCE BETWEEN SPHSOL ASC AND QM NUC
C                 BETWEEN SPHSOL ASC AND MM CHARGE
C                 BETWEEN SPHSOL ASC AND MM POL
C        NO FORCE BETWEEN SPHSOL ASC AND ASC
C
      IF(ISPHSOL.GE.60 .AND. RSPHSOL.LT.1.0D+30)
     *CALL QMMMRXNFMMNUCMM(X(LFFCORD),X(LFFFFGRD2),X(LFFCHARG),
     *                     X(LFFDIP),X(LFFQRXN1),X(LFFQRXN2),ISPHSOL,
     *                     X(LFFIDATOM),X(LFFIDDAI),X(LFFDAI),
     *                     X(LFFAFIX),X(LFFQRXNMP),
     *                     X(LFFQRXNTA),X(LFFQRXNXY),
     *                     X(LFFDIPMP),X(LFFDIPTA),X(LFFDIPXY))
C
C     ----- FORCE BETWEEN FIXSOL ASC AND ASC
C                 BETWEEN FIXSOL ASC AND QM NUC
C                 BETWEEN FIXSOL ASC AND MM CHARGE
C                 BETWEEN FIXSOL ASC AND MM POL
C
      IF(IFIXSOL.EQ.1)
     *CALL QMMMRXNFMMNUCMM(X(LFFCORD),X(LFFFFGRD2),X(LFFCHARG),
     *                     X(LFFDIP),X(LFFQFIX),X(LFFQRXN2),NFFTS,
     *                     X(LFFIDATOM),X(LFFIDDAI),X(LFFDAI),
     *                     X(LFFAFIX),X(LFFQFIXMP),
     *                     X(LFFQFIXTA),X(LFFQFIXXY),
     *                     X(LFFDIPMP),X(LFFDIPTA),X(LFFDIPXY))
C
C
      IF(.NOT.DFTBFL) THEN
C     ----- GLOBAL SUM OF FORCES ON MM ATOMS ----
C           -VIR- AND -FFGRD2- WERE ZEROED BEFORE CALLING GRADX
C           QM FORCES ARE IN -DE-, TO BE SUMMED IN STVDER
      CALL DDI_GSUMF(2418,VIR   ,3)
      CALL DDI_GSUMF(2419,X(LFFFFGRD2),3*NFFAT)
      END IF
C
      IF(DFTBFL) THEN
C        -- DFTB USES -EG- INSTEAD OF -DE-
C           MUST SAVE -EG-
         IF(GOPARR) CALL DDI_GSUMF(2417,DE,NAT*3)
         CALL DCOPY(NAT*3,DE,1,EG,1)
         CALL DAWRIT(IDAF,IODA,EG,3*NAT,3,0)
C        ----- GLOBAL SUM OF FORCES ON MM ATOMS ----
C              -VIR- AND -FFGRD2- WERE ZEROED BEFORE CALLING GRADX
         CALL DDI_GSUMF(2418,VIR   ,3)
         CALL DDI_GSUMF(2419,X(LFFFFGRD2),3*NFFAT)
         CALL DCOPY(NAT,X(LDTEMP),1,ZAN,1)
      END IF
C
  100 CONTINUE
      CALL RETFM(NEED)
C
      RETURN
      END
C*MODULE QUANPOC  *DECK QMMMCHGINT
!>
!> @brief    1-e integral of MM point charges
!>
      SUBROUTINE QMMMCHGINT(QQ,CORD,CHARG,NONLSTQ)
      use mx_limits, only: mxsh,mxgtot,mxatm
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      LOGICAL IANDJ,NORM,DOUBLE
      LOGICAL GOPARR,MASWRK,DSKWRK
C
      DIMENSION QQ(*),CHARG(*),CORD(3,*),NONLSTQ(*)
C
      DIMENSION DIJ(225),XIN(125),YIN(125),ZIN(125),
     *          IX(35),IY(35),IZ(35),JX(35),JY(35),JZ(35),
     *          IJX(225),IJY(225),IJZ(225)
      DIMENSION CHCINT(225)
      DIMENSION FIJ(225)
C
C
      COMMON /FFPARM/ NFFAT,NBOND,NANGL,NDIHR,NDIHB,NCMAP,NWAGG,
     *                N1213J,N14J,NLKQMM,IDOCHG,IDOPOL,IDOLJ,IDOCMAP
      COMMON /FFPBSW/ XBOX,YBOX,ZBOX,SWF,SWFDX,SWFDY,SWFDZ,
     *                SWRA,ONESWRA,SWRA2,ONESWRA2,
     *                SWRB,ONESWRB,SWRB2,ONESWRB2,ONESWRB4,
     *                SWFDUM3,SWFDUM4,SWFDUM5,
     *                SWRAQ,ONESWRAQ,SWRAQ2,ONESWRAQ2,
     *                SWRBQ,ONESWRBQ,SWRBQ2,ONESWRBQ2,ONESWRBQ4,
     *                SWFDUM3Q,SWFDUM4Q,SWFDUM5Q,
     *                QMSIZE,QMCX,QMCY,QMCZ,QMCXSV,QMCYSV,QMCZSV,
     *                CENTX,CENTY,CENTZ,BUFWID1,BUFWID2,RDAMP,
     *                EFIELDX,EFIELDY,EFIELDZ,QMCXSV2,QMCYSV2,QMCZSV2,
     *                EPS1RB,EPS1RB3,ONEXBOX,ONEYBOX,ONEZBOX,
     *                LQMCT,MXLIST1,NTODO,NTODOSV,NTODOQ,
     *                ISWITCH,ISHIFT,IPOLSHF,
     *                LFFLSTCELL,LFFCORDSV,
     *                LFFPOLSV,LFFCORDSV2,LFFNONLS2,LFFCORDSVQ,
     *                LFFMVFASTS2,LFFMVFASTS3,LFFMVFASTS4,
     *                LFFMVFASTL2,LFFMVFASTL3,LFFMVFASTL4,
     *                MXCHECK,MXLIST2,NTODO2,NTODO2SV
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /NSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),
     *                CF(MXGTOT),CG(MXGTOT),CH(MXGTOT),CI(MXGTOT),
     *                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),
     *                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
      COMMON /OUTPUT/ NPRINT,ITOL,ICUT,NORMF,NORMP,NOPK
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /ROOT  / XX,U(13),W(13),NROOTS
      COMMON /STV   / XINT,YINT,ZINT,T,X0,Y0,Z0,XI,YI,ZI,XJ,YJ,ZJ,NI,NJ
      COMMON /SYMIND/ TOL,II,JJ,LIT,LJT,MINI,MINJ,MAXI,MAXJ,IANDJ
C
      PARAMETER (ZERO=0.0D+00)
      PARAMETER (ONE=1.0D+00)
      PARAMETER (PI212=1.1283791670955D+00)
      PARAMETER (SQRT3=1.73205080756888D+00)
      PARAMETER (SQRT5=2.23606797749979D+00)
      PARAMETER (SQRT7=2.64575131106459D+00)
      PARAMETER (RLN10=2.30258D+00)
C
      DATA JX / 0, 1, 0, 0, 2, 0, 0, 1, 1, 0,
     *          3, 0, 0, 2, 2, 1, 0, 1, 0, 1,
     *          4, 0, 0, 3, 3, 1, 0, 1, 0, 2,
     *          2, 0, 2, 1, 1/
      DATA IX / 1, 6, 1, 1,11, 1, 1, 6, 6, 1,
     *         16, 1, 1,11,11, 6, 1, 6, 1, 6,
     *         21, 1, 1,16,16, 6, 1, 6, 1,11,
     *         11, 1,11, 6, 6/
      DATA JY / 0, 0, 1, 0, 0, 2, 0, 1, 0, 1,
     *          0, 3, 0, 1, 0, 2, 2, 0, 1, 1,
     *          0, 4, 0, 1, 0, 3, 3, 0, 1, 2,
     *          0, 2, 1, 2, 1/
      DATA IY / 1, 1, 6, 1, 1,11, 1, 6, 1, 6,
     *          1,16, 1, 6, 1,11,11, 1, 6, 6,
     *          1,21, 1, 6, 1,16,16, 1, 6,11,
     *          1,11, 6,11, 6/
      DATA JZ / 0, 0, 0, 1, 0, 0, 2, 0, 1, 1,
     *          0, 0, 3, 0, 1, 0, 1, 2, 2, 1,
     *          0, 0, 4, 0, 1, 0, 1, 3, 3, 0,
     *          2, 2, 1, 1, 2/
      DATA IZ / 1, 1, 1, 6, 1, 1,11, 1, 6, 6,
     *          1, 1,16, 1, 6, 1, 6,11,11, 6,
     *          1, 1,21, 1, 6, 1, 6,16,16, 1,
     *         11,11, 6, 6,11/
C
C     HUI LI, JAN 2011, LINCOLN
C
C     CALCULATE MM CHARGE CONTRIBUTION TO QM 1-E INTEGRALS
C
      IF(IDOCHG.EQ.0) RETURN
C
      IPCOUNT = ME - 1
C
      TOL = RLN10*ITOL
      NORM = NORMF .NE. 1 .OR. NORMP .NE. 1
C
C       HERE (AND SIMILAR PLACES) INTENTIONALLY DEACTIVATES SOME LOOPS,
C       WHILE ALLOWING THEM TO REMAIN PRESENT FOR FUTURE USE.
C
      L2 = (NUM*NUM+NUM)/2
      CALL VCLR(QQ,1,L2)
C
C     -- ISHELL
C
      DO 720 II = 1,NSHELL
         I = KATOM(II)
         XI = C(1,I)
         YI = C(2,I)
         ZI = C(3,I)
         I1 = KSTART(II)
         I2 = I1+KNG(II)-1
         LIT = KTYPE(II)
         MINI = KMIN(II)
         MAXI = KMAX(II)
         LOCI = KLOC(II)-MINI
C
C        -- JSHELL
C
         DO 700 JJ = 1,II
            IF (GOPARR) THEN
               IPCOUNT = IPCOUNT + 1
               IF (MOD(IPCOUNT,NPROC).NE.0) GOTO 700
            END IF
            J = KATOM(JJ)
            XJ = C(1,J)
            YJ = C(2,J)
            ZJ = C(3,J)
            J1 = KSTART(JJ)
            J2 = J1+KNG(JJ)-1
            LJT = KTYPE(JJ)
            MINJ = KMIN(JJ)
            MAXJ = KMAX(JJ)
            LOCJ = KLOC(JJ)-MINJ
            NROOTS = (LIT+LJT-2)/2+1
            RR = (XI-XJ)**2+(YI-YJ)**2+(ZI-ZJ)**2
            IANDJ = II .EQ. JJ
C
C           -- PREPARE INDICES FOR PAIRS OF (I,J) FUNCTIONS
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
  140          CONTINUE
  160       CONTINUE
            DO 180 I = 1,IJ
               CHCINT(I) = ZERO
  180       CONTINUE
C
C           -- I PRIMITIVE
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
C
C              -- J PRIMITIVE
C
               IF (IANDJ) JGMAX = IG
               DO 500 JG = J1,JGMAX
                  AJ = EX(JG)
                  AA = AI+AJ
                  AA1 = ONE/AA
                  DUM = AJ*ARRI*AA1
                  IF (DUM .GT. TOL) GOTO 500
                  FAC = EXP(-DUM)
                  CSJ = CS(JG)
                  CPJ = CP(JG)
                  CDJ = CD(JG)
                  CFJ = CF(JG)
                  CGJ = CG(JG)
                  AX = (AXI+AJ*XJ)*AA1
                  AY = (AYI+AJ*YJ)*AA1
                  AZ = (AZI+AJ*ZJ)*AA1
C
C                 -- DENSITY FACTOR
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
                     IF ((I.EQ.8).AND.NORM) DUM1=DUM1*SQRT3
                     IF (I.EQ.11) DUM1=CFI*FAC
                     IF ((I.EQ.14).AND.NORM) DUM1=DUM1*SQRT5
                     IF ((I.EQ.20).AND.NORM) DUM1=DUM1*SQRT3
                     IF (I.EQ.21) DUM1=CGI*FAC
                     IF ((I.EQ.24).AND.NORM) DUM1=DUM1*SQRT7
                     IF ((I.EQ.30).AND.NORM) DUM1=DUM1*SQRT5/SQRT3
                     IF ((I.EQ.33).AND.NORM) DUM1=DUM1*SQRT3
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
                        END IF
                        NN = NN+1
                        DIJ(NN) = DUM2
  200                CONTINUE
  220             CONTINUE
C
C                 -- CHARGES INTEGRALS.
C
                  DUM = PI212*AA1
                  DO 400 I=1,IJ
                     FIJ(I) = DIJ(I)*DUM
  400             CONTINUE
                  AAX = AA*AX
                  AAY = AA*AY
                  AAZ = AA*AZ
C
                  DO 495 IIQ = 1, NTODOQ
                     IFFAT = NONLSTQ(IIQ)
                     IF(ISWITCH.LE.1) THEN
                     CXGRP = CORD(1,IFFAT) - QMCX
                     CYGRP = CORD(2,IFFAT) - QMCY
                     CZGRP = CORD(3,IFFAT) - QMCZ
                     PBCX  = XBOX*ANINT(CXGRP*ONEXBOX)
                     PBCY  = YBOX*ANINT(CYGRP*ONEYBOX)
                     PBCZ  = ZBOX*ANINT(CZGRP*ONEZBOX)
                     CXGRP = CXGRP - PBCX
                     CYGRP = CYGRP - PBCY
                     CZGRP = CZGRP - PBCZ
                     R2    = CXGRP*CXGRP+CYGRP*CYGRP+CZGRP*CZGRP
                     IF(R2.GT.SWRBQ2) GOTO 495
                     CALL SWFUNCQ(R2,CXGRP,CYGRP,CZGRP)
                     ELSE
                     IAT   = KATOM(II)
                     JAT   = KATOM(JJ)
                     XIJ   = 0.5D+00*(C(1,IAT) + C(1,JAT))
                     YIJ   = 0.5D+00*(C(2,IAT) + C(2,JAT))
                     ZIJ   = 0.5D+00*(C(3,IAT) + C(3,JAT))
                     X     = CORD(1,IFFAT) - XIJ
                     Y     = CORD(2,IFFAT) - YIJ
                     Z     = CORD(3,IFFAT) - ZIJ
                     PBCX  = XBOX*ANINT(X*ONEXBOX)
                     PBCY  = YBOX*ANINT(Y*ONEYBOX)
                     PBCZ  = ZBOX*ANINT(Z*ONEZBOX)
                     X     = X - PBCX
                     Y     = Y - PBCY
                     Z     = Z - PBCZ
                     R2    = X*X+Y*Y+Z*Z
                     IF(R2.GT.SWRB2) GOTO 495
                     IF(R2.LT.1.0D-10) GOTO 495
                     R     = SQRT(R2)
                     ONER  = ONE/R
                     CALL SHIFT(R2,R,ONER,X,Y,Z)
                     END IF
C
                     CX    = CORD(1,IFFAT) - PBCX
                     CY    = CORD(2,IFFAT) - PBCY
                     CZ    = CORD(3,IFFAT) - PBCZ
                     ZNUC = -CHARG(IFFAT)
                     IF(ZNUC.EQ.ZERO) GOTO 495
                     XX = AA*((AX-CX)**2+(AY-CY)**2+(AZ-CZ)**2)
                     IF (NROOTS.LE.3) CALL RT123
                     IF (NROOTS.EQ.4) CALL ROOT4
                     IF (NROOTS.EQ.5) CALL ROOT5
                     MM = 0
                     DO 477 K = 1,NROOTS
                        UU = AA*U(K)
                        WW = W(K)*ZNUC
                        TT = ONE/(AA+UU)
                        T = SQRT(TT)
                        X0 = (AAX+UU*CX)*TT
                        Y0 = (AAY+UU*CY)*TT
                        Z0 = (AAZ+UU*CZ)*TT
                        IN = -5+MM
                        DO 476 I = 1,LIT
                           IN = IN+5
                           NI = I
                           DO 475 J = 1,LJT
                              JN = IN+J
                              NJ = J
                              CALL STVINT
                              XIN(JN) = XINT
                              YIN(JN) = YINT
                              ZIN(JN) = ZINT*WW
  475                      CONTINUE
  476                   CONTINUE
                        MM = MM+25
  477                CONTINUE
                     DO 481 I = 1,IJ
                        NX = IJX(I)
                        NY = IJY(I)
                        NZ = IJZ(I)
                        DUM = ZERO
                        MM = 0
                        DO 479 K = 1,NROOTS
                           DUM = DUM+XIN(NX+MM)*YIN(NY+MM)*ZIN(NZ+MM)
                           MM = MM+25
  479                   CONTINUE
                        CHCINT(I) = CHCINT(I)+DUM*FIJ(I)*SWF
  481                CONTINUE
  495             CONTINUE
  500          CONTINUE
  520       CONTINUE
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
                  QQ(JN) = QQ(JN) + CHCINT(NN)
  600          CONTINUE
  620       CONTINUE
  700    CONTINUE
  720 CONTINUE
C
      RETURN
      END
C*MODULE QUANPOC  *DECK QMMMCHGFMMSWF
!>
!> @brief    QM electron -MM charge interaction energy
!>
!> @author   Hui Li
!>           - Jan 2011
!>
!> @details  This energy is scaled by SWF to produce force
!>
      SUBROUTINE QMMMCHGFMMSWF(DM,CHARG,CORD,FFGRD,NONLSTQ)
      use mx_limits, only: mxsh,mxgtot,mxatm
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      LOGICAL IANDJ,NORM
      LOGICAL GOPARR,MASWRK,DSKWRK
C
      DIMENSION DM(*),CHARG(*),CORD(3,*),FFGRD(3,*),NONLSTQ(*)
C
      DIMENSION DIJ(225),XIN(125),YIN(125),ZIN(125),
     *          IX(35),IY(35),IZ(35),JX(35),JY(35),JZ(35),
     *          IJX(225),IJY(225),IJZ(225)
      DIMENSION FIJ(225)
C
C
      COMMON /FFPARM/ NFFAT,NBOND,NANGL,NDIHR,NDIHB,NCMAP,NWAGG,
     *                N1213J,N14J,NLKQMM,IDOCHG,IDOPOL,IDOLJ,IDOCMAP
      COMMON /FFMDPA/ DT,DT2,TEMP0,PRES0,POLTOL,VIR(3),PMEAN,VOLAV,
     *                ENPAV,ENKAV,TEMPAV,BERENDT,BERENDP,VELMAX,
     *                PMEANX,PMEANY,PMEANZ,NSTEP,KMASTER,KOUTACT(2),
     *                IHESS,INTALG,ITSTAT,IPSTAT,JOUT,KOUT,LOUT
      COMMON /FFPBSW/ XBOX,YBOX,ZBOX,SWF,SWFDX,SWFDY,SWFDZ,
     *                SWRA,ONESWRA,SWRA2,ONESWRA2,
     *                SWRB,ONESWRB,SWRB2,ONESWRB2,ONESWRB4,
     *                SWFDUM3,SWFDUM4,SWFDUM5,
     *                SWRAQ,ONESWRAQ,SWRAQ2,ONESWRAQ2,
     *                SWRBQ,ONESWRBQ,SWRBQ2,ONESWRBQ2,ONESWRBQ4,
     *                SWFDUM3Q,SWFDUM4Q,SWFDUM5Q,
     *                QMSIZE,QMCX,QMCY,QMCZ,QMCXSV,QMCYSV,QMCZSV,
     *                CENTX,CENTY,CENTZ,BUFWID1,BUFWID2,RDAMP,
     *                EFIELDX,EFIELDY,EFIELDZ,QMCXSV2,QMCYSV2,QMCZSV2,
     *                EPS1RB,EPS1RB3,ONEXBOX,ONEYBOX,ONEZBOX,
     *                LQMCT,MXLIST1,NTODO,NTODOSV,NTODOQ,
     *                ISWITCH,ISHIFT,IPOLSHF,
     *                LFFLSTCELL,LFFCORDSV,
     *                LFFPOLSV,LFFCORDSV2,LFFNONLS2,LFFCORDSVQ,
     *                LFFMVFASTS2,LFFMVFASTS3,LFFMVFASTS4,
     *                LFFMVFASTL2,LFFMVFASTL3,LFFMVFASTL4,
     *                MXCHECK,MXLIST2,NTODO2,NTODO2SV
      COMMON /FFRMSD/ DIPT(8),TIMGYRA,TIMRALL,
     *                NATPDB,NGYRA,NDIEL,NRALL,
     *                LFFRALL0,NRMSD,LFFRMSD0,KFREEAB(201),
     *                NRIJMM,IJRMM(2,100),NRIJQM,IJRQM(2,100),
     *                NAIJKMM,IJKMM(3,100),NAIJKQM,IJKQM(3,100),
     *                NFIXMM,IFIXMM(200),NFIXQM,IFIXQM(200)
      COMMON /FMCOM / X(1)
      COMMON /GRAD  / DE(3,MXATM)
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /NSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),
     *                CF(MXGTOT),CG(MXGTOT),CH(MXGTOT),CI(MXGTOT),
     *                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),
     *                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
      COMMON /OUTPUT/ NPRINT,ITOL,ICUT,NORMF,NORMP,NOPK
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /ROOT  / XX,U(13),W(13),NROOTS
      COMMON /STV   / XINT,YINT,ZINT,T,X0,Y0,Z0,XI,YI,ZI,XJ,YJ,ZJ,NI,NJ
      COMMON /SYMIND/ TOL,II,JJ,LIT,LJT,MINI,MINJ,MAXI,MAXJ,IANDJ
C
      PARAMETER (ZERO=0.0D+00)
      PARAMETER (PT5=0.5D+00)
      PARAMETER (ONE=1.0D+00)
      PARAMETER (PI212=1.1283791670955D+00)
      PARAMETER (SQRT3=1.73205080756888D+00)
      PARAMETER (SQRT5=2.23606797749979D+00)
      PARAMETER (SQRT7=2.64575131106459D+00)
      PARAMETER (RLN10=2.30258D+00)
C
      DATA JX / 0, 1, 0, 0, 2, 0, 0, 1, 1, 0,
     *          3, 0, 0, 2, 2, 1, 0, 1, 0, 1,
     *          4, 0, 0, 3, 3, 1, 0, 1, 0, 2,
     *          2, 0, 2, 1, 1/
      DATA IX / 1, 6, 1, 1,11, 1, 1, 6, 6, 1,
     *         16, 1, 1,11,11, 6, 1, 6, 1, 6,
     *         21, 1, 1,16,16, 6, 1, 6, 1,11,
     *         11, 1,11, 6, 6/
      DATA JY / 0, 0, 1, 0, 0, 2, 0, 1, 0, 1,
     *          0, 3, 0, 1, 0, 2, 2, 0, 1, 1,
     *          0, 4, 0, 1, 0, 3, 3, 0, 1, 2,
     *          0, 2, 1, 2, 1/
      DATA IY / 1, 1, 6, 1, 1,11, 1, 6, 1, 6,
     *          1,16, 1, 6, 1,11,11, 1, 6, 6,
     *          1,21, 1, 6, 1,16,16, 1, 6,11,
     *          1,11, 6,11, 6/
      DATA JZ / 0, 0, 0, 1, 0, 0, 2, 0, 1, 1,
     *          0, 0, 3, 0, 1, 0, 1, 2, 2, 1,
     *          0, 0, 4, 0, 1, 0, 1, 3, 3, 0,
     *          2, 2, 1, 1, 2/
      DATA IZ / 1, 1, 1, 6, 1, 1,11, 1, 6, 6,
     *          1, 1,16, 1, 6, 1, 6,11,11, 6,
     *          1, 1,21, 1, 6, 1, 6,16,16, 1,
     *         11,11, 6, 6,11/
C
C     HUI LI, JAN 2011, LINCOLN
C
C     MM CHARGE - QM ELECTRON INTERACTION ENERGY
C     THIS ENERGY IS ONLY USED WITH SWFDX TO PRODUCE GRADIENT
C
      IF(IDOCHG.EQ.0) RETURN
C
      NNQ = NTODOQ
C
      CALL VALFM(LOADFM)
      LSWF    = LOADFM + 1
      LAST    = LSWF   + 3*225*NNQ
      NEED    = LAST   - LOADFM
      CALL GETFM(NEED)
      CALL VCLR(X(LSWF),1,3*225*NNQ)
C
      IPCOUNT = ME - 1
      TOL = RLN10*ITOL
      NORM = NORMF .NE. 1 .OR. NORMP .NE. 1
C
C     -- ISHELL
C
      DO 720 II = 1,NSHELL
         I = KATOM(II)
         XI = C(1,I)
         YI = C(2,I)
         ZI = C(3,I)
         I1 = KSTART(II)
         I2 = I1+KNG(II)-1
         LIT = KTYPE(II)
         MINI = KMIN(II)
         MAXI = KMAX(II)
         LOCI = KLOC(II)-MINI
C
C        -- JSHELL
C
         DO 700 JJ = 1,NSHELL
            IF (GOPARR) THEN
               IPCOUNT = IPCOUNT + 1
               IF (MOD(IPCOUNT,NPROC).NE.0) GOTO 700
            END IF
            J = KATOM(JJ)
            XJ = C(1,J)
            YJ = C(2,J)
            ZJ = C(3,J)
            J1 = KSTART(JJ)
            J2 = J1+KNG(JJ)-1
            LJT = KTYPE(JJ)
            MINJ = KMIN(JJ)
            MAXJ = KMAX(JJ)
            LOCJ = KLOC(JJ)-MINJ
            NROOTS = (LIT+LJT-2)/2+1
            RR = (XI-XJ)**2+(YI-YJ)**2+(ZI-ZJ)**2
C
C           -- PREPARE INDICES FOR PAIRS OF (I,J) FUNCTIONS
C
            IJ = 0
            MAX = MAXJ
            DO 160 I = MINI,MAXI
               NX = IX(I)
               NY = IY(I)
               NZ = IZ(I)
               DO 140 J = MINJ,MAX
                  IJ = IJ+1
                  IJX(IJ) = NX+JX(J)
                  IJY(IJ) = NY+JY(J)
                  IJZ(IJ) = NZ+JZ(J)
  140          CONTINUE
  160       CONTINUE
            CALL VCLR(X(LSWF),1,3*225*NNQ)
C
C           -- I PRIMITIVE
C
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
C
C              -- J PRIMITIVE
C
               DO 500 JG = J1,J2
                  AJ = EX(JG)
                  AA = AI+AJ
                  AA1 = ONE/AA
                  DUM = AJ*ARRI*AA1
                  IF (DUM .GT. TOL) GOTO 500
                  FAC = EXP(-DUM)
                  CSJ = CS(JG)
                  CPJ = CP(JG)
                  CDJ = CD(JG)
                  CFJ = CF(JG)
                  CGJ = CG(JG)
                  AX = (AXI+AJ*XJ)*AA1
                  AY = (AYI+AJ*YJ)*AA1
                  AZ = (AZI+AJ*ZJ)*AA1
C
C                 -- DENSITY FACTOR
C
                  MAX = MAXJ
                  NN = 0
                  DUM1 = ZERO
                  DUM2 = ZERO
                  DO 220 I = MINI,MAXI
                     IF (I.EQ.1) DUM1=CSI*FAC
                     IF (I.EQ.2) DUM1=CPI*FAC
                     IF (I.EQ.5) DUM1=CDI*FAC
                     IF ((I.EQ.8).AND.NORM) DUM1=DUM1*SQRT3
                     IF (I.EQ.11) DUM1=CFI*FAC
                     IF ((I.EQ.14).AND.NORM) DUM1=DUM1*SQRT5
                     IF ((I.EQ.20).AND.NORM) DUM1=DUM1*SQRT3
                     IF (I.EQ.21) DUM1=CGI*FAC
                     IF ((I.EQ.24).AND.NORM) DUM1=DUM1*SQRT7
                     IF ((I.EQ.30).AND.NORM) DUM1=DUM1*SQRT5/SQRT3
                     IF ((I.EQ.33).AND.NORM) DUM1=DUM1*SQRT3
                     DO 200 J = MINJ,MAX
                        IF (J.EQ.1) THEN
                           DUM2=DUM1*CSJ
                        ELSE IF (J.EQ.2) THEN
                           DUM2=DUM1*CPJ
                        ELSE IF (J.EQ.5) THEN
                           DUM2=DUM1*CDJ
                        ELSE IF ((J.EQ.8).AND.NORM) THEN
                           DUM2=DUM2*SQRT3
                        ELSE IF (J.EQ.11) THEN
                           DUM2=DUM1*CFJ
                        ELSE IF ((J.EQ.14).AND.NORM) THEN
                           DUM2=DUM2*SQRT5
                        ELSE IF ((J.EQ.20).AND.NORM) THEN
                           DUM2=DUM2*SQRT3
                        ELSE IF (J.EQ.21) THEN
                           DUM2=DUM1*CGJ
                        ELSE IF ((J.EQ.24).AND.NORM) THEN
                           DUM2=DUM2*SQRT7
                        ELSE IF ((J.EQ.30).AND.NORM) THEN
                           DUM2=DUM2*SQRT5/SQRT3
                        ELSE IF ((J.EQ.33).AND.NORM) THEN
                           DUM2=DUM2*SQRT3
                        END IF
                        NN = NN+1
                        DIJ(NN) = DUM2
  200                CONTINUE
  220             CONTINUE
C
C                 -- CHARGES INTEGRALS.
C
                  DUM = PI212*AA1
                  DO 400 I=1,IJ
                     FIJ(I) = DIJ(I)*DUM
  400             CONTINUE
                  AAX = AA*AX
                  AAY = AA*AY
                  AAZ = AA*AZ
C
                  DO 495 IIQ = 1, NNQ
                     IFFAT = NONLSTQ(IIQ)
                     IF(ISWITCH.LE.1) THEN
                     CXGRP = CORD(1,IFFAT) - QMCX
                     CYGRP = CORD(2,IFFAT) - QMCY
                     CZGRP = CORD(3,IFFAT) - QMCZ
                     PBCX  = XBOX*ANINT(CXGRP*ONEXBOX)
                     PBCY  = YBOX*ANINT(CYGRP*ONEYBOX)
                     PBCZ  = ZBOX*ANINT(CZGRP*ONEZBOX)
                     CXGRP = CXGRP - PBCX
                     CYGRP = CYGRP - PBCY
                     CZGRP = CZGRP - PBCZ
                     R2    = CXGRP*CXGRP+CYGRP*CYGRP+CZGRP*CZGRP
                     IF(R2.GT.SWRBQ2) GOTO 495
                     CALL SWFUNCQ(R2,CXGRP,CYGRP,CZGRP)
                     ELSE
                     IAT   = KATOM(II)
                     JAT   = KATOM(JJ)
                     XIJ   = 0.5D+00*(C(1,IAT) + C(1,JAT))
                     YIJ   = 0.5D+00*(C(2,IAT) + C(2,JAT))
                     ZIJ   = 0.5D+00*(C(3,IAT) + C(3,JAT))
                     DX    = CORD(1,IFFAT) - XIJ
                     DY    = CORD(2,IFFAT) - YIJ
                     DZ    = CORD(3,IFFAT) - ZIJ
                     PBCX  = XBOX*ANINT(DX*ONEXBOX)
                     PBCY  = YBOX*ANINT(DY*ONEYBOX)
                     PBCZ  = ZBOX*ANINT(DZ*ONEZBOX)
                     DX    = DX - PBCX
                     DY    = DY - PBCY
                     DZ    = DZ - PBCZ
                     R2    = DX*DX+DY*DY+DZ*DZ
                     IF(R2.GT.SWRB2) GOTO 495
                     IF(R2.LT.1.0D-10) GOTO 495
                     R     = SQRT(R2)
                     ONER  = ONE/R
                     CALL SHIFT(R2,R,ONER,DX,DY,DZ)
                     END IF
C
                     CX    = CORD(1,IFFAT) - PBCX
                     CY    = CORD(2,IFFAT) - PBCY
                     CZ    = CORD(3,IFFAT) - PBCZ
                     ZNUC = -CHARG(IFFAT)
                     XX = AA*((AX-CX)**2+(AY-CY)**2+(AZ-CZ)**2)
                     IF (NROOTS.LE.3) CALL RT123
                     IF (NROOTS.EQ.4) CALL ROOT4
                     IF (NROOTS.EQ.5) CALL ROOT5
                     MM = 0
                     DO 477 K = 1,NROOTS
                        UU = AA*U(K)
                        WW = W(K)*ZNUC
                        TT = ONE/(AA+UU)
                        T = SQRT(TT)
                        X0 = (AAX+UU*CX)*TT
                        Y0 = (AAY+UU*CY)*TT
                        Z0 = (AAZ+UU*CZ)*TT
                        IN = -5+MM
                        DO 476 I = 1,LIT
                           IN = IN+5
                           NI = I
                           DO 475 J = 1,LJT
                              JN = IN+J
                              NJ = J
                              CALL STVINT
                              XIN(JN) = XINT
                              YIN(JN) = YINT
                              ZIN(JN) = ZINT*WW
  475                      CONTINUE
  476                   CONTINUE
                        MM = MM+25
  477                CONTINUE
                     DO 481 I = 1,IJ
                        NX = IJX(I)
                        NY = IJY(I)
                        NZ = IJZ(I)
                        DUM = ZERO
                        MM = 0
                        DO 479 K = 1,NROOTS
                           DUM = DUM+XIN(NX+MM)*YIN(NY+MM)*ZIN(NZ+MM)
                           MM = MM+25
  479                   CONTINUE
                        DUMFIJ = DUM*FIJ(I)
                        X(LSWF+(IIQ-1)*3*225+(I-1)*3  )=
     *                  X(LSWF+(IIQ-1)*3*225+(I-1)*3  )+DUMFIJ*SWFDX
                        X(LSWF+(IIQ-1)*3*225+(I-1)*3+1)=
     *                  X(LSWF+(IIQ-1)*3*225+(I-1)*3+1)+DUMFIJ*SWFDY
                        X(LSWF+(IIQ-1)*3*225+(I-1)*3+2)=
     *                  X(LSWF+(IIQ-1)*3*225+(I-1)*3+2)+DUMFIJ*SWFDZ
  481                CONTINUE
  495             CONTINUE
  500          CONTINUE
  520       CONTINUE
C
C
            MAX = MAXJ
            KK = 0
            DO 620 I = MINI,MAXI
               LI = LOCI+I
               DO 600 J = MINJ,MAX
                  KK = KK + 1
                  LJ = LOCJ+J
                  IF(LI-LJ) 920,940,940
  920             ID = LJ
                  JD = LI
                  GOTO 960
  940             ID = LI
                  JD = LJ
  960             NN = (ID*(ID-1))/2+JD
                  DUM = DM(NN)
                  DO IIQ = 1, NNQ
                     IFFAT = NONLSTQ(IIQ)
                     IF(ISWITCH.LE.1) THEN
                     CXGRP = CORD(1,IFFAT) - QMCX
                     CYGRP = CORD(2,IFFAT) - QMCY
                     CZGRP = CORD(3,IFFAT) - QMCZ
                     ELSE
                     IAT   = KATOM(II)
                     JAT   = KATOM(JJ)
                     XIJ   = 0.5D+00*(C(1,IAT) + C(1,JAT))
                     YIJ   = 0.5D+00*(C(2,IAT) + C(2,JAT))
                     ZIJ   = 0.5D+00*(C(3,IAT) + C(3,JAT))
                     CXGRP = CORD(1,IFFAT) - XIJ
                     CYGRP = CORD(2,IFFAT) - YIJ
                     CZGRP = CORD(3,IFFAT) - ZIJ
                     END IF
                     PBCX  = XBOX*ANINT(CXGRP*ONEXBOX)
                     PBCY  = YBOX*ANINT(CYGRP*ONEYBOX)
                     PBCZ  = ZBOX*ANINT(CZGRP*ONEZBOX)
                     CXGRP = CXGRP - PBCX
                     CYGRP = CYGRP - PBCY
                     CZGRP = CZGRP - PBCZ
                     DEX   = DUM*X(LSWF+(IIQ-1)*3*225+(KK-1)*3  )
                     DEY   = DUM*X(LSWF+(IIQ-1)*3*225+(KK-1)*3+1)
                     DEZ   = DUM*X(LSWF+(IIQ-1)*3*225+(KK-1)*3+2)
                     FFGRD(1,IFFAT)=FFGRD(1,IFFAT)+DEX
                     FFGRD(2,IFFAT)=FFGRD(2,IFFAT)+DEY
                     FFGRD(3,IFFAT)=FFGRD(3,IFFAT)+DEZ
                     IF(ISWITCH.EQ.2) THEN    ! ONLY ATOM-ATOM SHIFTING NEEDS THIS
                     DE(1,IAT) = DE(1,IAT) - PT5*DEX
                     DE(2,IAT) = DE(2,IAT) - PT5*DEY
                     DE(3,IAT) = DE(3,IAT) - PT5*DEZ
                     DE(1,JAT) = DE(1,JAT) - PT5*DEX
                     DE(2,JAT) = DE(2,JAT) - PT5*DEY
                     DE(3,JAT) = DE(3,JAT) - PT5*DEZ
                     END IF
                     VIR(1)=VIR(1)+DEX*CXGRP
                     VIR(2)=VIR(2)+DEY*CYGRP
                     VIR(3)=VIR(3)+DEZ*CZGRP
                     IYES = 0
                     DO KFIX=1,NFIXMM
                        IF(IFFAT.EQ.IFIXMM(KFIX)) IYES = 1
                     ENDDO
                     IF(IYES.EQ.1.AND.NFIXQM.GT.0) THEN
                        VIR(1)=VIR(1)-DEX*CXGRP
                        VIR(2)=VIR(2)-DEY*CYGRP
                        VIR(3)=VIR(3)-DEZ*CZGRP
                     END IF
                  ENDDO
  600          CONTINUE
  620       CONTINUE
  700    CONTINUE
  720 CONTINUE
C
      CALL RETFM(NEED)
      RETURN
      END
C*MODULE QUANPOC  *DECK QMMMCHGFMM
!>
!> @brief    QM electron -MM charge force
!>
!> @author   Hui Li
!>           - Jan 2011
!>
!> @details  The forces act on MM charges by QM electrons
!>
      SUBROUTINE QMMMCHGFMM(DM,CHGINT,CORD,CHARG,FFGRD,NONLSTQ)
      use mx_limits, only: mxsh,mxgtot,mxatm
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      LOGICAL NORM,GOPARR,DSKWRK,MASWRK,NXT
C
C
      DIMENSION DM(*),CHGINT(*),CORD(3,*),CHARG(*),FFGRD(3,*),
     *          NONLSTQ(*)
C
      DIMENSION DIJ(100),XIN(432),YIN(432),ZIN(432)
      DIMENSION FIJ(100)
      DIMENSION IX(20),IY(20),IZ(20),JX(20),JY(20),JZ(20)
      DIMENSION IJX(100),IJY(100),IJZ(100)
C
      COMMON /CSSTV / CX,CY,CZ
      COMMON /FFPARM/ NFFAT,NBOND,NANGL,NDIHR,NDIHB,NCMAP,NWAGG,
     *                N1213J,N14J,NLKQMM,IDOCHG,IDOPOL,IDOLJ,IDOCMAP
      COMMON /FFMDPA/ DT,DT2,TEMP0,PRES0,POLTOL,VIR(3),PMEAN,VOLAV,
     *                ENPAV,ENKAV,TEMPAV,BERENDT,BERENDP,VELMAX,
     *                PMEANX,PMEANY,PMEANZ,NSTEP,KMASTER,KOUTACT(2),
     *                IHESS,INTALG,ITSTAT,IPSTAT,JOUT,KOUT,LOUT
      COMMON /FFMPT3/ NACTMM,LACTMM(2020),LFFDIPOLD,JUMBUP,
     *                NACTQM,LACTQM(2020),LFFOLDC,LFFQMVELSV,MMHESS,
     *                LFFQMCHG,LFFQMCHGB,ISWAP,R2SWAP,DFTBMM
      COMMON /FFPBSW/ XBOX,YBOX,ZBOX,SWF,SWFDX,SWFDY,SWFDZ,
     *                SWRA,ONESWRA,SWRA2,ONESWRA2,
     *                SWRB,ONESWRB,SWRB2,ONESWRB2,ONESWRB4,
     *                SWFDUM3,SWFDUM4,SWFDUM5,
     *                SWRAQ,ONESWRAQ,SWRAQ2,ONESWRAQ2,
     *                SWRBQ,ONESWRBQ,SWRBQ2,ONESWRBQ2,ONESWRBQ4,
     *                SWFDUM3Q,SWFDUM4Q,SWFDUM5Q,
     *                QMSIZE,QMCX,QMCY,QMCZ,QMCXSV,QMCYSV,QMCZSV,
     *                CENTX,CENTY,CENTZ,BUFWID1,BUFWID2,RDAMP,
     *                EFIELDX,EFIELDY,EFIELDZ,QMCXSV2,QMCYSV2,QMCZSV2,
     *                EPS1RB,EPS1RB3,ONEXBOX,ONEYBOX,ONEZBOX,
     *                LQMCT,MXLIST1,NTODO,NTODOSV,NTODOQ,
     *                ISWITCH,ISHIFT,IPOLSHF,
     *                LFFLSTCELL,LFFCORDSV,
     *                LFFPOLSV,LFFCORDSV2,LFFNONLS2,LFFCORDSVQ,
     *                LFFMVFASTS2,LFFMVFASTS3,LFFMVFASTS4,
     *                LFFMVFASTL2,LFFMVFASTL3,LFFMVFASTL4,
     *                MXCHECK,MXLIST2,NTODO2,NTODO2SV
      COMMON /FFRMSD/ DIPT(8),TIMGYRA,TIMRALL,
     *                NATPDB,NGYRA,NDIEL,NRALL,
     *                LFFRALL0,NRMSD,LFFRMSD0,KFREEAB(201),
     *                NRIJMM,IJRMM(2,100),NRIJQM,IJRQM(2,100),
     *                NAIJKMM,IJKMM(3,100),NAIJKQM,IJKQM(3,100),
     *                NFIXMM,IFIXMM(200),NFIXQM,IFIXQM(200)
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /IOFILE/ IR,IW,IP,IJK,IJKT,IDAF,NAV,IODA(950)
      COMMON /NSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),
     *                CF(MXGTOT),CG(MXGTOT),CH(MXGTOT),CI(MXGTOT),
     *                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),
     *                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
      COMMON /OUTPUT/ NPRINT,ITOL,ICUT,NORMF,NORMP,NOPK
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /ROOT  / XX,U(13),W(13),NROOTS
      COMMON /STV   / XINT,YINT,ZINT,T,X0,Y0,Z0,XI,YI,ZI,XJ,YJ,ZJ,NI,NJ
C
      DATA JX / 0, 1, 0, 0, 2, 0, 0, 1, 1, 0,
     1          3, 0, 0, 2, 2, 1, 0, 1, 0, 1/
      DATA JY / 0, 0, 1, 0, 0, 2, 0, 1, 0, 1,
     1          0, 3, 0, 1, 0, 2, 2, 0, 1, 1/
      DATA JZ / 0, 0, 0, 1, 0, 0, 2, 0, 1, 1,
     1          0, 0, 3, 0, 1, 0, 1, 2, 2, 1/
      DATA IX / 1, 5, 1, 1, 9, 1, 1, 5, 5, 1,
     1         13, 1, 1, 9, 9, 5, 1, 5, 1, 5/
      DATA IY / 1, 1, 5, 1, 1, 9, 1, 5, 1, 5,
     1          1,13, 1, 5, 1, 9, 9, 1, 5, 5/
      DATA IZ / 1, 1, 1, 5, 1, 1, 9, 1, 5, 5,
     1          1, 1,13, 1, 5, 1, 5, 9, 9, 5/
C
      DATA PI212 /1.1283791670955D+00/
      DATA SQRT3 /1.73205080756888D+00/
      DATA SQRT5 /2.23606797749979D+00/
      DATA ZERO,ONE/0.0D+00,1.0D+00/
      DATA RLN10 /2.30258D+00/
C
C     HUI LI, JAN 2011, LINCOLN
C
C     FORCES ON MM CHARGES BY QM ELECTRONS
C
      IF(IDOCHG.EQ.0) RETURN
C
      CALL QMMMCHGFMMSWF(DM,CHARG,CORD,FFGRD,NONLSTQ)
C
C
      NXT = IBTYP.EQ.1
      IPCOUNT = ME - 1
      NEXT = -1
      LCNT = -1
      TOL = RLN10*ITOL
      NORM = NORMF .NE. 1 .OR. NORMP .NE. 1
C
      NNQ = NTODOQ
      IF(NACTMM.GT.0) NNQ = NACTMM
      DO 100 IIQ = 1, NNQ
      IFFAT = NONLSTQ(IIQ)
      IF(NACTMM.GT.0) IFFAT = LACTMM(IIQ)
      ZNUC = CHARG(IFFAT)
      IF(ZNUC.EQ.ZERO) GOTO 100
      IF(ISWITCH.LE.1) THEN
      CXGRP = CORD(1,IFFAT) - QMCX
      CYGRP = CORD(2,IFFAT) - QMCY
      CZGRP = CORD(3,IFFAT) - QMCZ
      PBCX  = XBOX*ANINT(CXGRP*ONEXBOX)
      PBCY  = YBOX*ANINT(CYGRP*ONEYBOX)
      PBCZ  = ZBOX*ANINT(CZGRP*ONEZBOX)
      CXGRP = CXGRP - PBCX
      CYGRP = CYGRP - PBCY
      CZGRP = CZGRP - PBCZ
      R2    = CXGRP*CXGRP+CYGRP*CYGRP+CZGRP*CZGRP
      IF(R2.GT.SWRBQ2) GOTO 100
      CALL SWFUNCQ(R2,CXGRP,CYGRP,CZGRP)
      CX    = CORD(1,IFFAT) - PBCX
      CY    = CORD(2,IFFAT) - PBCY
      CZ    = CORD(3,IFFAT) - PBCZ
      END IF
C
C     -- ISHELL
C
      DO 600 II = 1,NSHELL
C
      IF(NXT .AND. GOPARR) THEN
         LCNT = LCNT + 1
         IF(LCNT.GT.NEXT) CALL DDI_DLBNEXT(NEXT)
         IF(NEXT.NE.LCNT) GOTO 600
      END IF
      IAT= KATOM(II)
      XI = C(1,IAT)
      YI = C(2,IAT)
      ZI = C(3,IAT)
      I1 = KSTART(II)
      I2 = I1+KNG(II)-1
      LIT = KTYPE(II)
      MINI = KMIN(II)
      MAXI = KMAX(II)
      LOCI = KLOC(II)-MINI
C
C     -- JSHELL
C
      DO 580 JJ = 1,NSHELL
      IF((.NOT.NXT) .AND. GOPARR) THEN
         IPCOUNT = IPCOUNT + 1
         IF(MOD(IPCOUNT,NPROC).NE.0) GOTO 580
      END IF
      JAT= KATOM(JJ)
      XJ = C(1,JAT)
      YJ = C(2,JAT)
      ZJ = C(3,JAT)
      J1 = KSTART(JJ)
      J2 = J1+KNG(JJ)-1
      LJT = KTYPE(JJ)
      MINJ = KMIN(JJ)
      MAXJ = KMAX(JJ)
      LOCJ = KLOC(JJ)-MINJ
      NROOTS = (LIT+LJT+1-2)/2 + 1
      RR = (XI-XJ)**2+(YI-YJ)**2+(ZI-ZJ)**2
C
C     -- PREPARE INDICES FOR PAIRS OF (I,J) FUNCTIONS
C
      IJ = 0
      MAX = MAXJ
      DO 50 I = MINI,MAXI
      NX = IX(I)
      NY = IY(I)
      NZ = IZ(I)
      DO 50 J = MINJ,MAX
      IJ = IJ+1
      IJX(IJ) = NX+JX(J)
      IJY(IJ) = NY+JY(J)
      IJZ(IJ) = NZ+JZ(J)
  50  CONTINUE
      CALL VCLR(CHGINT,1,3*IJ)
C
C     -- I PRIMITIVE
C
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
C
C     -- J PRIMITIVE
C
      DO 500 JG = J1,J2
      AJ = EX(JG)
      AA = AI+AJ
      AA1 = ONE/AA
      DUM = AJ*ARRI*AA1
      IF (DUM .GT. TOL) GOTO 500
      FAC =  EXP(-DUM)
      CSJ = CS(JG)
      CPJ = CP(JG)
      CDJ = CD(JG)
      CFJ = CF(JG)
      AX = (AXI+AJ*XJ)*AA1
      AY = (AYI+AJ*YJ)*AA1
      AZ = (AZI+AJ*ZJ)*AA1
C
C     -- DENSITY FACTOR
C
      MAX = MAXJ
      NN = 0
      DUM1=ZERO
      DUM2=DUM1
      DO 170 I=MINI,MAXI
      GOTO (70,80,110,110,90,110,110,95,110,110,
     1       102,110,110,104,110,110,110,110,110,108),I
  70  DUM1=CSI*FAC
      GOTO 110
  80  DUM1=CPI*FAC
      GOTO 110
  90  DUM1=CDI*FAC
      GOTO 110
  95  IF(NORM) DUM1=DUM1*SQRT3
      GOTO 110
 102  DUM1=CFI*FAC
      GOTO 110
 104  DUM1 = DUM1 *SQRT5
      GOTO 110
 108  DUM1 = DUM1 * SQRT3
 110  CONTINUE
      DO 170 J=MINJ,MAX
      GOTO (125,130,160,160,140,160,160,150,160,160,
     1       152,160,160,154,160,160,160,160,160,156),J
  125 DUM2=DUM1*CSJ
      GOTO 160
  130 DUM2=DUM1*CPJ
      GOTO 160
  140 DUM2=DUM1*CDJ
      GOTO 160
  150 IF(NORM) DUM2=DUM2*SQRT3
      GOTO 160
  152 DUM2 = DUM1 * CFJ
      GOTO 160
  154 DUM2 = DUM2 *SQRT5
      GOTO 160
  156 DUM2 = DUM2 * SQRT3
  160 NN=NN+1
  170 DIJ(NN)=DUM2
C
C     ..... HELLMANN-FEYNMAN TERM .....
C
      DUM = PI212*AA1
      DUM=DUM+DUM
      DO 380 I = 1,IJ
      FIJ(I) = DIJ(I)*DUM
  380 CONTINUE
      AAX = AA*AX
      AAY = AA*AY
      AAZ = AA*AZ
C
      ICC=1
C
      IF(ISWITCH.EQ.2) THEN
      IAT   = KATOM(II)
      JAT   = KATOM(JJ)
      XIJ   = 0.5D+00*(C(1,IAT) + C(1,JAT))
      YIJ   = 0.5D+00*(C(2,IAT) + C(2,JAT))
      ZIJ   = 0.5D+00*(C(3,IAT) + C(3,JAT))
      X     = CORD(1,IFFAT) - XIJ
      Y     = CORD(2,IFFAT) - YIJ
      Z     = CORD(3,IFFAT) - ZIJ
      PBCX  = XBOX*ANINT(X*ONEXBOX)
      PBCY  = YBOX*ANINT(Y*ONEYBOX)
      PBCZ  = ZBOX*ANINT(Z*ONEZBOX)
      X     = X - PBCX
      Y     = Y - PBCY
      Z     = Z - PBCZ
      R2    = X*X+Y*Y+Z*Z
      IF(R2.GT.SWRB2) GOTO 500
      IF(R2.LT.1.0D-10) GOTO 500
      R     = SQRT(R2)
      ONER  = ONE/R
      CALL SHIFT(R2,R,ONER,X,Y,Z)
      CX    = CORD(1,IFFAT) - PBCX
      CY    = CORD(2,IFFAT) - PBCY
      CZ    = CORD(3,IFFAT) - PBCZ
      END IF
C
      XX   = AA*((AX-CX)**2+(AY-CY)**2+(AZ-CZ)**2)
      IF(NROOTS.LE.3) CALL RT123
      IF(NROOTS.EQ.4) CALL ROOT4
      IF(NROOTS.EQ.5) CALL ROOT5
      IF(NROOTS.EQ.6) CALL ROOT6
      IF(NROOTS.GE.7) THEN
         WRITE(IW,9008)
         CALL ABRT
      END IF
      MM = 0
      DO 401 K = 1,NROOTS
      UU = AA*U(K)
      WW = W(K)*ZNUC
      WW = WW*UU
      TT = ONE/(AA+UU)
      T  = SQRT(TT)
      X0 = (AAX+UU*CX)*TT
      Y0 = (AAY+UU*CY)*TT
      Z0 = (AAZ+UU*CZ)*TT
      IN = -4+MM
      DO 400 I = 1,LIT
      IN = IN+4
      NI = I
      DO 400 J = 1,LJT
      JN = IN+J
      NJ = J
      CALL STVINT
      XIN(JN   ) = XINT
      YIN(JN   ) = YINT
      ZIN(JN   ) = ZINT*WW
      CALL POLXYZ
      XIN(JN+125) = XINT
      YIN(JN+125) = YINT
      ZIN(JN+125) = ZINT*WW
  400 CONTINUE
  401 MM = MM+16
      DO 403 I = 1,IJ
      NX    = IJX(I)
      NY    = IJY(I)
      NZ    = IJZ(I)
      DUMX  = ZERO
      DUMY  = ZERO
      DUMZ  = ZERO
      MM    = 0
      DO 402 K = 1,NROOTS
      DUMX= DUMX+XIN(NX+MM+125)*YIN(NY+MM    )*ZIN(NZ+MM    )
      DUMY= DUMY+XIN(NX+MM    )*YIN(NY+MM+125)*ZIN(NZ+MM    )
      DUMZ= DUMZ+XIN(NX+MM    )*YIN(NY+MM    )*ZIN(NZ+MM+125)
  402 MM  = MM+16
      DUM = FIJ(I)*SWF
      CHGINT(ICC  )=CHGINT(ICC  )-DUM*DUMX
      CHGINT(ICC+1)=CHGINT(ICC+1)-DUM*DUMY
      CHGINT(ICC+2)=CHGINT(ICC+2)-DUM*DUMZ
      ICC=ICC+3
 403  CONTINUE
 500  CONTINUE
 520  CONTINUE
C
C     ----- END OF *PRIMITIVE* LOOPS -----
C
      MAX=MAXJ
      ICC=1
      DO 550 I=MINI,MAXI
      LI=LOCI+I
      DO 550 J=MINJ,MAX
      LJ=LOCJ+J
      IF (LI-LJ) 920,940,940
  920 ID = LJ
      JD = LI
      GOTO 960
  940 ID = LI
      JD = LJ
  960 NN = (ID*(ID-1))/2+JD
      DUM = DM(NN)
      FFGRD(1,IFFAT)=FFGRD(1,IFFAT) + DUM*CHGINT(ICC  )
      FFGRD(2,IFFAT)=FFGRD(2,IFFAT) + DUM*CHGINT(ICC+1)
      FFGRD(3,IFFAT)=FFGRD(3,IFFAT) + DUM*CHGINT(ICC+2)
C
C     -- HUI LI: BEFORE DEC 10, 2016, VIRIAL IS CALCULATED WITH CX,
C        WHICH IS MM CORD INSTEAD OF QM-MM DISTANCE, BECAUSE
C        WE DO NOT KNOW WHICH ELECTRON (AND WHERE IT IS) GIVES
C        THE FROCE TO THE MM CHARGE POINT (AND IT SHOULD NOT
C        BE DONE WITH QMCX BECAUSE THE FORCE COMES FROM ELECTRONS
C        INSTEAD OF QMCX SWITCHING FUNCTION)
C     VIR(1)        =VIR(1)         + DUM*CHGINT(ICC  )*CX
C     VIR(2)        =VIR(2)         + DUM*CHGINT(ICC+1)*CY
C     VIR(3)        =VIR(3)         + DUM*CHGINT(ICC+2)*CZ
C
C     -- HUI LI: AFTER DEC 10, 2016, VIRIAL IS CALCULATED WITH CXGRP,
C        WHICH IS QM-MM DISTANCE, AND THE
C        QM ELECTRON IS APPROXIMATELY AT XIJ, WHERE IT GIVES THE FORCE TO
C        THE MM CHARGE POINT (AND IT SHOULD NOT BE DONE WITH
C        QMCX BECAUSE THE FORCE COMES FROM ELECTRONS
C        INSTEAD OF QMCX SWITCHING FUNCTION).
      IAT   = KATOM(II)
      JAT   = KATOM(JJ)
      XIJ   = 0.5D+00*(C(1,IAT) + C(1,JAT))
      YIJ   = 0.5D+00*(C(2,IAT) + C(2,JAT))
      ZIJ   = 0.5D+00*(C(3,IAT) + C(3,JAT))
      CXGRP = CORD(1,IFFAT) - XIJ
      CYGRP = CORD(2,IFFAT) - YIJ
      CZGRP = CORD(3,IFFAT) - ZIJ
      PBCX  = XBOX*ANINT(CXGRP*ONEXBOX)
      PBCY  = YBOX*ANINT(CYGRP*ONEYBOX)
      PBCZ  = ZBOX*ANINT(CZGRP*ONEZBOX)
      CXGRP = CXGRP - PBCX
      CYGRP = CYGRP - PBCY
      CZGRP = CZGRP - PBCZ
      VIR(1)        =VIR(1)         + DUM*CHGINT(ICC  )*CXGRP
      VIR(2)        =VIR(2)         + DUM*CHGINT(ICC+1)*CYGRP
      VIR(3)        =VIR(3)         + DUM*CHGINT(ICC+2)*CZGRP
      IYES = 0
      DO KFIX=1,NFIXMM
         IF(IFFAT.EQ.IFIXMM(KFIX)) IYES = 1
      ENDDO
      IF(IYES.EQ.1.AND.NFIXQM.GT.0) THEN
         VIR(1)     =VIR(1)         - DUM*CHGINT(ICC  )*CXGRP
         VIR(2)     =VIR(2)         - DUM*CHGINT(ICC+1)*CYGRP
         VIR(3)     =VIR(3)         - DUM*CHGINT(ICC+2)*CZGRP
      END IF
      ICC=ICC+3
  550 CONTINUE
  580 CONTINUE
  600 CONTINUE
  100 CONTINUE
C
      IF(GOPARR .AND. NXT) CALL DDI_DLBRESET
      RETURN
 9008 FORMAT(/' NUMBER OF POLYNOMIAL ROOTS NEEDED (NROOTS) IS GREATER',
     *        ' THAN 6 IN EFDINT.  CALL A PROGRAMMER/QUANTUM CHEMIST.')
      END
C*MODULE QUANPOC  *DECK QMMMCHGFQM
      SUBROUTINE QMMMCHGFQM(DM,CORD,CHARG,NONLSTQ)
      use mx_limits, only: mxsh,mxgtot,mxatm
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      LOGICAL ISKIP
      LOGICAL NORM,GOPARR,DSKWRK,MASWRK,NXT
C
      DIMENSION DM(*),CORD(3,*),CHARG(*),NONLSTQ(*)
C
      DIMENSION ISKIP(35),INDX(70)
      DIMENSION IX(35),IY(35),IZ(35),JX(20),JY(20),JZ(20)
      DIMENSION IJG(210),IJX(210),IJY(210),IJZ(210)
      DIMENSION XIN(200),YIN(200),ZIN(200),G(210),DIJ(210)
      DIMENSION FIJ(210)
C
C
      COMMON /DSTV  / XINT,YINT,ZINT,T,X0,Y0,Z0,XI,YI,ZI,
     *                XJ,YJ,ZJ,NI,NJ,CX,CY,CZ
      COMMON /FFPARM/ NFFAT,NBOND,NANGL,NDIHR,NDIHB,NCMAP,NWAGG,
     *                N1213J,N14J,NLKQMM,IDOCHG,IDOPOL,IDOLJ,IDOCMAP
      COMMON /FFPBSW/ XBOX,YBOX,ZBOX,SWF,SWFDX,SWFDY,SWFDZ,
     *                SWRA,ONESWRA,SWRA2,ONESWRA2,
     *                SWRB,ONESWRB,SWRB2,ONESWRB2,ONESWRB4,
     *                SWFDUM3,SWFDUM4,SWFDUM5,
     *                SWRAQ,ONESWRAQ,SWRAQ2,ONESWRAQ2,
     *                SWRBQ,ONESWRBQ,SWRBQ2,ONESWRBQ2,ONESWRBQ4,
     *                SWFDUM3Q,SWFDUM4Q,SWFDUM5Q,
     *                QMSIZE,QMCX,QMCY,QMCZ,QMCXSV,QMCYSV,QMCZSV,
     *                CENTX,CENTY,CENTZ,BUFWID1,BUFWID2,RDAMP,
     *                EFIELDX,EFIELDY,EFIELDZ,QMCXSV2,QMCYSV2,QMCZSV2,
     *                EPS1RB,EPS1RB3,ONEXBOX,ONEYBOX,ONEZBOX,
     *                LQMCT,MXLIST1,NTODO,NTODOSV,NTODOQ,
     *                ISWITCH,ISHIFT,IPOLSHF,
     *                LFFLSTCELL,LFFCORDSV,
     *                LFFPOLSV,LFFCORDSV2,LFFNONLS2,LFFCORDSVQ,
     *                LFFMVFASTS2,LFFMVFASTS3,LFFMVFASTS4,
     *                LFFMVFASTL2,LFFMVFASTL3,LFFMVFASTL4,
     *                MXCHECK,MXLIST2,NTODO2,NTODO2SV
      COMMON /GRAD  / DE(3,MXATM)
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /NSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),
     *                CF(MXGTOT),CG(MXGTOT),CH(MXGTOT),CI(MXGTOT),
     *                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),
     *                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
      COMMON /OUTPUT/ NPRINT,ITOL,ICUT,NORMF,NORMP,NOPK
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /RUNOPT/ RUNTYP,EXETYP,NEVALS,NGLEVL,NHLEVL
      COMMON /ROOT  / XX,U(13),W(13),NROOTS
C
      PARAMETER (ZERO=0.0D+00)
      PARAMETER (ONE=1.0D+00)
      PARAMETER (RLN10=2.30258D+00)
      PARAMETER (PI212=1.1283791670955D+00)
      PARAMETER (SQRT3=1.73205080756888D+00)
      PARAMETER (SQRT5=2.23606797749979D+00)
C
      DATA CHECK/8HCHECK   /
      DATA INDX / 1, 2, 3, 4, 5, 6, 7, 8, 9,10,
     1           11,12,13,14,15,16,17,18,19,20,
     2           -0,-0,-0,-0,-0,-0,-0,-0,-0,-0,
     3           -0,-0,-0,-0,-0,
     4           -0, 1, 2, 3, 1, 2, 3, 4, 5, 6,
     5            4, 5, 6, 7, 8, 9,10,11,12,13,
     6            7, 8, 9,10,11,12,13,14,15,16,
     7           17,18,19,20,21/
      DATA JX / 0, 1, 0, 0, 2, 0, 0, 1, 1, 0,
     1          3, 0, 0, 2, 2, 1, 0, 1, 0, 1/
      DATA IX / 1, 5, 1, 1, 9, 1, 1, 5, 5, 1,
     1         13, 1, 1, 9, 9, 5, 1, 5, 1, 5,
     2         17, 1, 1,13,13, 5, 1, 5, 1, 9,
     3          9, 1, 9, 5, 5/
      DATA JY / 0, 0, 1, 0, 0, 2, 0, 1, 0, 1,
     1          0, 3, 0, 1, 0, 2, 2, 0, 1, 1/
      DATA IY / 1, 1, 5, 1, 1, 9, 1, 5, 1, 5,
     1          1,13, 1, 5, 1, 9, 9, 1, 5, 5,
     2          1,17, 1, 5, 1,13,13, 1, 5, 9,
     3          1, 9, 5, 9, 5/
      DATA JZ / 0, 0, 0, 1, 0, 0, 2, 0, 1, 1,
     1          0, 0, 3, 0, 1, 0, 1, 2, 2, 1/
      DATA IZ / 1, 1, 1, 5, 1, 1, 9, 1, 5, 5,
     1          1, 1,13, 1, 5, 1, 5, 9, 9, 5,
     2          1, 1,17, 1, 5, 1, 5,13,13, 1,
     3          9, 9, 5, 5, 9/
C
C     HUI LI, JAN 2011, LINCOLN
C
C     FORCES ON QM ELECTRONS BY MM CHARGES
C
      IF(IDOCHG.EQ.0) RETURN
C
      NXT = IBTYP.EQ.1
      IPCOUNT = ME - 1
      NEXT = -1
      LCNT = -1
      TOL = RLN10*ITOL
      NORM = NORMF .NE. 1 .OR. NORMP .NE. 1
C
      IF (EXETYP .EQ. CHECK) GOTO 1120
C
C     -- ISHELL
C
      DO 1060 II = 1,NSHELL
C
      IF(NXT .AND. GOPARR) THEN
         LCNT = LCNT + 1
         IF(LCNT.GT.NEXT) CALL DDI_DLBNEXT(NEXT)
         IF(NEXT.NE.LCNT) GOTO 1060
      END IF
      IAT = KATOM(II)
      XI = C(1,IAT)
      YI = C(2,IAT)
      ZI = C(3,IAT)
      I1 = KSTART(II)
      I2 = I1+KNG(II)-1
      LIT = KTYPE(II)+1
      MINI = KMIN(II)
      MAXI = KMAX(II)
      LOCI = KLOC(II)-MINI
      DO 100 I=1,35
 100  ISKIP(I)=.TRUE.
      DO 260 I=MINI,MAXI
      GOTO (110,140,260,260,180,260,260,260,260,260,
     *       220,260,260,260,260,260,260,260,260,260),I
 110  DO 120 K=2,4
 120  ISKIP(K)=.FALSE.
      GOTO 260
 140  DO 160 K=5,10
 160  ISKIP(K)=.FALSE.
      ISKIP(1)=.FALSE.
      GOTO 260
 180  DO 200 K=2,4
 200  ISKIP(K)=.FALSE.
      DO 210 K=11,20
 210  ISKIP(K)=.FALSE.
      GOTO 260
 220  DO 240 K=5,10
 240  ISKIP(K)=.FALSE.
      DO 250 K=21,35
 250  ISKIP(K)=.FALSE.
 260  CONTINUE
C
C     -- JSHELL
C
      DO 1040 JJ = 1,NSHELL
        IF((.NOT.NXT) .AND. GOPARR) THEN
           IPCOUNT = IPCOUNT + 1
           IF(MOD(IPCOUNT,NPROC).NE.0) GOTO 1040
        END IF
C
      JAT = KATOM(JJ)
      XJ = C(1,JAT)
      YJ = C(2,JAT)
      ZJ = C(3,JAT)
      J1 = KSTART(JJ)
      J2 = J1+KNG(JJ)-1
      LJT = KTYPE(JJ)
      MINJ = KMIN(JJ)
      MAXJ = KMAX(JJ)
      LOCJ = KLOC(JJ)-MINJ
      NROOTS = (LIT+LJT+1-2)/2 + 1
      RR = (XI-XJ)**2+(YI-YJ)**2+(ZI-ZJ)**2
C
C     -- PREPARE INDICES FOR PAIRS OF (I,J) FUNCTIONS
C
      N0 = 0
      IF (LIT .GE. 4) N0 = 35
      IJ = 0
      DO 340 I = 1,35
      IF (ISKIP(I)) GOTO 340
      IN = INDX(I+N0)
      NX = IX(I)
      NY = IY(I)
      NZ = IZ(I)
      DO 320 J = MINJ,MAXJ
      IJ = IJ+1
      IJX(IJ) = NX+JX(J)
      IJY(IJ) = NY+JY(J)
      IJZ(IJ) = NZ+JZ(J)
      IJG(IJ) = IN+21*(J-MINJ)
  320 CONTINUE
  340 CONTINUE
      DO 360 I = 1,IJ
      N = IJG(I)
  360 G(N) = ZERO
C
C     -- I PRIMITIVE
C
      DO 840 IG = I1,I2
      AI = EX(IG)
      ARRI = AI*RR
      AXI = AI*XI
      AYI = AI*YI
      AZI = AI*ZI
      DUM = AI+AI
      CSI = CP(IG)
      CPI = CS(IG)*DUM
      IF (LIT .EQ. 4) CPI = CD(IG)
      CDI = CP(IG)*DUM
      IF (LIT .EQ. 5) CDI = CF(IG)
      CFI = CD(IG)*DUM
      CGI = CF(IG)*DUM
C
C     -- J PRIMITIVE
C
      DO 820 JG = J1,J2
      AJ = EX(JG)
      AA = AI+AJ
      AA1 = ONE/AA
      DUM = AJ*ARRI*AA1
      IF (DUM .GT. TOL) GOTO 820
      FAC = EXP(-DUM)
      CSJ = CS(JG)
      CPJ = CP(JG)
      CDJ = CD(JG)
      CFJ = CF(JG)
      AX = (AXI+AJ*XJ)*AA1
      AY = (AYI+AJ*YJ)*AA1
      AZ = (AZI+AJ*ZJ)*AA1
C
C     -- DENSITY FACTOR
C
      NN = 0
      DUM1 = ZERO
      DUM2 = DUM1
      DO 600 I=1,35
      IF(ISKIP(I)) GOTO 600
      GOTO (370,380,420,420,390,420,420,420,420,420,
     *       400,420,420,420,420,420,420,420,420,420,
     *       410,420,420,420,420,420,420,420,420,420,
     *       420,420,420,420,420),I
 370  DUM1=CSI*FAC
      GOTO 420
 380  DUM1=CPI*FAC
      GOTO 420
 390  DUM1=CDI*FAC
      GOTO 420
 400  DUM1=CFI*FAC
      GOTO 420
 410  DUM1 = CGI*FAC
 420  CONTINUE
      DO 580 J=MINJ,MAXJ
      GOTO (430,440,560,560,460,560,560,480,560,560,
     *       500,560,560,520,560,560,560,560,560,540),J
  430 DUM2=DUM1*CSJ
      GOTO 560
  440 DUM2=DUM1*CPJ
      GOTO 560
  460 DUM2=DUM1*CDJ
      GOTO 560
  480 IF(NORM) DUM2=DUM2*SQRT3
      GOTO 560
  500 DUM2 = DUM1 * CFJ
      GOTO 560
  520 DUM2 = DUM2 *SQRT5
      GOTO 560
  540 DUM2 = DUM2 * SQRT3
  560 NN=NN+1
  580 DIJ(NN)=DUM2
  600 CONTINUE
C
C     ....BASIC CHARGE TERM
C
      DUM = PI212*AA1
      DO 700 I=1,IJ
 700  FIJ(I)=DIJ(I)*DUM
      AAX=AA*AX
      AAY=AA*AY
      AAZ=AA*AZ
C
      DO 710 IIQ = 1, NTODOQ
      IFFAT = NONLSTQ(IIQ)
      IF(ISWITCH.LE.1) THEN
      CXGRP = CORD(1,IFFAT) - QMCX
      CYGRP = CORD(2,IFFAT) - QMCY
      CZGRP = CORD(3,IFFAT) - QMCZ
      PBCX  = XBOX*ANINT(CXGRP*ONEXBOX)
      PBCY  = YBOX*ANINT(CYGRP*ONEYBOX)
      PBCZ  = ZBOX*ANINT(CZGRP*ONEZBOX)
      CXGRP = CXGRP - PBCX
      CYGRP = CYGRP - PBCY
      CZGRP = CZGRP - PBCZ
      R2    = CXGRP*CXGRP+CYGRP*CYGRP+CZGRP*CZGRP
      IF(R2.GT.SWRBQ2) GOTO 710
      CALL SWFUNCQ(R2,CXGRP,CYGRP,CZGRP)
      ELSE
      IAT   = KATOM(II)
      JAT   = KATOM(JJ)
      XIJ   = 0.5D+00*(C(1,IAT) + C(1,JAT))
      YIJ   = 0.5D+00*(C(2,IAT) + C(2,JAT))
      ZIJ   = 0.5D+00*(C(3,IAT) + C(3,JAT))
      X     = CORD(1,IFFAT) - XIJ
      Y     = CORD(2,IFFAT) - YIJ
      Z     = CORD(3,IFFAT) - ZIJ
      PBCX  = XBOX*ANINT(X*ONEXBOX)
      PBCY  = YBOX*ANINT(Y*ONEYBOX)
      PBCZ  = ZBOX*ANINT(Z*ONEZBOX)
      X     = X - PBCX
      Y     = Y - PBCY
      Z     = Z - PBCZ
      R2    = X*X+Y*Y+Z*Z
      IF(R2.GT.SWRB2) GOTO 710
      IF(R2.LT.1.0D-10) GOTO 710
      R     = SQRT(R2)
      ONER  = ONE/R
      CALL SHIFT(R2,R,ONER,X,Y,Z)
      END IF
C
      CX    = CORD(1,IFFAT) - PBCX
      CY    = CORD(2,IFFAT) - PBCY
      CZ    = CORD(3,IFFAT) - PBCZ
C
      ZNUC = -CHARG(IFFAT)
      IF(ZNUC.EQ.ZERO) GOTO 710
      XX=AA*((AX-CX)**2+(AY-CY)**2+(AZ-CZ)**2)
      IF(NROOTS.LE.3)CALL RT123
      IF(NROOTS.EQ.4)CALL ROOT4
      IF(NROOTS.EQ.5)CALL ROOT5
      MM=0
      DO 806 K=1,NROOTS
      UU=AA*U(K)
      WW=W(K)*ZNUC
      TT=ONE/(AA+UU)
      T =SQRT(TT)
      X0=(AAX+UU*CX)*TT
      Y0=(AAY+UU*CY)*TT
      Z0=(AAZ+UU*CZ)*TT
      IN=-4+MM
      DO 804 I=1,LIT
      IN=IN+4
      NI=I
      DO 804 J=1,LJT
      JN=IN+J
      NJ=J
      CALL VINT
      XIN(JN)=XINT
      YIN(JN)=YINT
      ZIN(JN)=ZINT*WW
 804  CONTINUE
 806  MM=MM+20
      DO 810 I=1,IJ
      N=IJG(I)
      NX=IJX(I)
      NY=IJY(I)
      NZ=IJZ(I)
      DUM=ZERO
      MM=0
      DO 808 K=1,NROOTS
      DUM=DUM+XIN(NX+MM)*YIN(NY+MM)*ZIN(NZ+MM)
 808  MM=MM+20
 810  G(N)=G(N)+DUM*FIJ(I)*SWF
 710  CONTINUE
C
C     -- END OF PRIMITIVE LOOPS -----
C
  820 CONTINUE
  840 CONTINUE
C
C     -- FORM INTEGRALS OVER DERIVATIVES -----
C
      NN = 0
      N = 1
      DO 900 J = MINJ,MAXJ
      IF (MINI .GT. 1) GOTO 860
      NN = NN+1
      XIN(NN) = G(N+ 1)
      YIN(NN) = G(N+ 2)
      ZIN(NN) = G(N+ 3)
      IF (MAXI .EQ. 1) GOTO 900
  860 IF (MINI .GT. 2) GOTO 880
      NN = NN+1
      XIN(NN) = (G(N+ 4)-G(N ))
      YIN(NN) = G(N+ 7)
      ZIN(NN) = G(N+ 8)
      NN = NN+1
      XIN(NN) = G(N+ 7)
      YIN(NN) = (G(N+ 5)-G(N ))
      ZIN(NN) = G(N+ 9)
      NN = NN+1
      XIN(NN) = G(N+ 8)
      YIN(NN) = G(N+ 9)
      ZIN(NN) = (G(N+ 6)-G(N ))
      IF (MAXI .EQ. 4) GOTO 900
  880 CONTINUE
      IF (MINI .GT. 5) GOTO 890
      NN = NN+1
      XIN(NN) = (G(N+ 3)-G(N )-G(N ))
      YIN(NN) = G(N+ 6)
      ZIN(NN) = G(N+ 7)
      NN = NN+1
      XIN(NN) = G(N+ 8)
      YIN(NN) = (G(N+ 4)-G(N+ 1)-G(N+ 1))
      ZIN(NN) = G(N+ 9)
      NN = NN+1
      XIN(NN) = G(N+10)
      YIN(NN) = G(N+11)
      ZIN(NN) = (G(N+ 5)-G(N+ 2)-G(N+ 2))
      NN = NN+1
      DUM = ONE
      IF (NORM) DUM = SQRT3
      XIN(NN) = DUM*(G(N+ 6)-G(N+ 1))
      YIN(NN) = DUM*(G(N+ 8)-G(N ))
      ZIN(NN) = DUM* G(N+12)
      NN = NN+1
      XIN(NN) = DUM*(G(N+ 7)-G(N+ 2))
      YIN(NN) = DUM* G(N+12)
      ZIN(NN) = DUM*(G(N+10)-G(N ))
      NN = NN+1
      XIN(NN) = DUM* G(N+12)
      YIN(NN) = DUM*(G(N+ 9)-G(N+ 2))
      ZIN(NN) = DUM*(G(N+11)-G(N+ 1))
      IF(MAXI.EQ.10) GOTO 900
 890  CONTINUE
      NN=NN+1
      XIN(NN)=(G(N+ 6)-G(N   )-G(N   )-G(N   ))
      YIN(NN)= G(N+ 9)
      ZIN(NN)= G(N+10)
      NN=NN+1
      XIN(NN)= G(N+11)
      YIN(NN)=(G(N+ 7)-G(N+ 1)-G(N+ 1)-G(N+ 1))
      ZIN(NN)= G(N+12)
      NN=NN+1
      XIN(NN)= G(N+13)
      YIN(NN)= G(N+14)
      ZIN(NN)=(G(N+ 8)-G(N+ 2)-G(N+ 2)-G(N+ 2))
      NN=NN+1
      DUM=ONE
      IF(NORM) DUM=SQRT5
      XIN(NN)=DUM* (G(N+ 9)-G(N+ 3)-G(N+ 3))
      YIN(NN)=DUM* (G(N+15)-G(N  ))
      ZIN(NN)=DUM* G(N+18)
      NN=NN+1
      XIN(NN)=DUM* (G(N+10)-G(N+ 4)-G(N+ 4))
      YIN(NN)=DUM* G(N+18)
      ZIN(NN)=DUM* (G(N+16)-G(N  ))
      NN=NN+1
      XIN(NN)=DUM* (G(N+15)-G(N+ 1))
      YIN(NN)=DUM* (G(N+11)-G(N+ 3)-G(N+ 3))
      ZIN(NN)=DUM*  G(N+19)
      NN=NN+1
      XIN(NN)=DUM* G(N+19)
      YIN(NN)=DUM* (G(N+12)-G(N+ 5)-G(N+ 5))
      ZIN(NN)=DUM* (G(N+17)-G(N+ 1))
      NN=NN+1
      XIN(NN)=DUM* (G(N+16)-G(N+ 2))
      YIN(NN)=DUM*  G(N+20)
      ZIN(NN)=DUM* (G(N+13)-G(N+ 4)-G(N+ 4))
      NN=NN+1
      XIN(NN)=DUM*  G(N+20)
      YIN(NN)=DUM* (G(N+17)-G(N+ 2))
      ZIN(NN)=DUM* (G(N+14)-G(N+ 5)-G(N+ 5))
      NN=NN+1
      IF(NORM) DUM=DUM*SQRT3
      XIN(NN)=DUM* (G(N+18)-G(N+ 5))
      YIN(NN)=DUM* (G(N+19)-G(N+ 4))
      ZIN(NN)=DUM* (G(N+20)-G(N+ 3))
  900 N = N+21
C
C     -- CALCULATE CONTRIBUTION TO GRADIENT -----
C
      N = 0
      DO 980 J = MINJ,MAXJ
      JN = LOCJ+J
      DO 980 I = MINI,MAXI
      N = N+1
      IN = LOCI+I
      IF (IN-JN) 920,940,940
  920 ID = JN
      JD = IN
      GOTO 960
  940 ID = IN
      JD = JN
  960 NN = (ID*(ID-1))/2+JD
      DUM = DM(NN)
      DUM = DUM+DUM
      DE(1,IAT) = DE(1,IAT)+DUM*XIN(N)
      DE(2,IAT) = DE(2,IAT)+DUM*YIN(N)
      DE(3,IAT) = DE(3,IAT)+DUM*ZIN(N)
  980 CONTINUE
 1040 CONTINUE
 1060 CONTINUE
C
C     ----- END OF SHELL LOOPS -----
C
 1120 CONTINUE
      IF(GOPARR .AND. NXT) CALL DDI_DLBRESET
      RETURN
C
      END
C*MODULE QUANPOC  *DECK QMMMCHGNUC
!>
!> @brief    QM nuc -MM charge force
!>
!> @author   Hui Li
!>           - Jan 2011
!>
!> @details  mutual forces are calculated
!>
      SUBROUTINE QMMMCHGNUC(CORD,FFGRD,DETMP,CHARG,NONLSTQ)
      use mx_limits, only: mxatm
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      LOGICAL GOPARR,DSKWRK,MASWRK
C
      PARAMETER (ZERO=0.0D+00)
      PARAMETER (ONE=1.0D+00)
C
      DIMENSION CORD(3,*),FFGRD(3,*),CHARG(*),DETMP(3,*),
     *          NONLSTQ(*)
C
      COMMON /FFENGY/ EN12,EN123,EN123R4,EN123B4,EN234W1,ENCHAR,ENLJR,
     *                ENLJD,ENPOL,XENPOL,ENRXN,XENRXN,ENRXNPOL,ENRXNR,
     *                EN12312,ENQUANP(30),
     *                ENBIAS,ENCENT,ENUCCH,ENCMAP,ENPOT,ENKIN,ENTOT
      COMMON /FFMPT3/ NACTMM,LACTMM(2020),LFFDIPOLD,JUMBUP,
     *                NACTQM,LACTQM(2020),LFFOLDC,LFFQMVELSV,MMHESS,
     *                LFFQMCHG,LFFQMCHGB,ISWAP,R2SWAP,DFTBMM
      COMMON /FFPARM/ NFFAT,NBOND,NANGL,NDIHR,NDIHB,NCMAP,NWAGG,
     *                N1213J,N14J,NLKQMM,IDOCHG,IDOPOL,IDOLJ,IDOCMAP
      COMMON /FFMDPA/ DT,DT2,TEMP0,PRES0,POLTOL,VIR(3),PMEAN,VOLAV,
     *                ENPAV,ENKAV,TEMPAV,BERENDT,BERENDP,VELMAX,
     *                PMEANX,PMEANY,PMEANZ,NSTEP,KMASTER,KOUTACT(2),
     *                IHESS,INTALG,ITSTAT,IPSTAT,JOUT,KOUT,LOUT
      COMMON /FFPBSW/ XBOX,YBOX,ZBOX,SWF,SWFDX,SWFDY,SWFDZ,
     *                SWRA,ONESWRA,SWRA2,ONESWRA2,
     *                SWRB,ONESWRB,SWRB2,ONESWRB2,ONESWRB4,
     *                SWFDUM3,SWFDUM4,SWFDUM5,
     *                SWRAQ,ONESWRAQ,SWRAQ2,ONESWRAQ2,
     *                SWRBQ,ONESWRBQ,SWRBQ2,ONESWRBQ2,ONESWRBQ4,
     *                SWFDUM3Q,SWFDUM4Q,SWFDUM5Q,
     *                QMSIZE,QMCX,QMCY,QMCZ,QMCXSV,QMCYSV,QMCZSV,
     *                CENTX,CENTY,CENTZ,BUFWID1,BUFWID2,RDAMP,
     *                EFIELDX,EFIELDY,EFIELDZ,QMCXSV2,QMCYSV2,QMCZSV2,
     *                EPS1RB,EPS1RB3,ONEXBOX,ONEYBOX,ONEZBOX,
     *                LQMCT,MXLIST1,NTODO,NTODOSV,NTODOQ,
     *                ISWITCH,ISHIFT,IPOLSHF,
     *                LFFLSTCELL,LFFCORDSV,
     *                LFFPOLSV,LFFCORDSV2,LFFNONLS2,LFFCORDSVQ,
     *                LFFMVFASTS2,LFFMVFASTS3,LFFMVFASTS4,
     *                LFFMVFASTL2,LFFMVFASTL3,LFFMVFASTL4,
     *                MXCHECK,MXLIST2,NTODO2,NTODO2SV
      COMMON /FFRMSD/ DIPT(8),TIMGYRA,TIMRALL,
     *                NATPDB,NGYRA,NDIEL,NRALL,
     *                LFFRALL0,NRMSD,LFFRMSD0,KFREEAB(201),
     *                NRIJMM,IJRMM(2,100),NRIJQM,IJRQM(2,100),
     *                NAIJKMM,IJKMM(3,100),NAIJKQM,IJKQM(3,100),
     *                NFIXMM,IFIXMM(200),NFIXQM,IFIXQM(200)
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
C
C     HUI LI, JAN 2011, LINCOLN
C
C     MM CHARGE - QM NUC INTERACTION
C
      CALL VCLR(DETMP,1,3*MXATM)
      ENUCCH =ZERO
      IF(IDOCHG.EQ.0) RETURN
C
      IPCOUNT = ME - 1
      DO 200 IIQ = 1, NTODOQ
         IFFAT = NONLSTQ(IIQ)
         QI = CHARG(IFFAT)
         IF(QI.EQ.ZERO) GOTO 200
         IF(GOPARR) THEN
            IPCOUNT = IPCOUNT + 1
            IF (MOD(IPCOUNT,NPROC).NE.0) GOTO 200
         END IF
         IF(ISWITCH.LE.1) THEN
         CXGRP = CORD(1,IFFAT) - QMCX
         CYGRP = CORD(2,IFFAT) - QMCY
         CZGRP = CORD(3,IFFAT) - QMCZ
         PBCX  = XBOX*ANINT(CXGRP*ONEXBOX)
         PBCY  = YBOX*ANINT(CYGRP*ONEYBOX)
         PBCZ  = ZBOX*ANINT(CZGRP*ONEZBOX)
         CXGRP = CXGRP - PBCX
         CYGRP = CYGRP - PBCY
         CZGRP = CZGRP - PBCZ
         R2    = CXGRP*CXGRP+CYGRP*CYGRP+CZGRP*CZGRP
         IF(R2.GT.SWRBQ2) GOTO 200
         CALL SWFUNCQ(R2,CXGRP,CYGRP,CZGRP)
         CX    = CORD(1,IFFAT) - PBCX
         CY    = CORD(2,IFFAT) - PBCY
         CZ    = CORD(3,IFFAT) - PBCZ
         END IF
         DO 210 JAT=1,NAT
            IF(ZAN(JAT).EQ.ZERO) GOTO 210
            IF(ISWITCH.LE.1) THEN
            X     = CX - C(1,JAT)
            Y     = CY - C(2,JAT)
            Z     = CZ - C(3,JAT)
            R2    = X*X+Y*Y+Z*Z+DFTBMM   ! FOR BOTH QM AND DFTB
            IF(R2.LT.1.0D-10) GOTO 210
            ONER2 = ONE/R2
            ONER  = SQRT(ONER2)
            ELSE
            X     = CORD(1,IFFAT) - C(1,JAT)
            Y     = CORD(2,IFFAT) - C(2,JAT)
            Z     = CORD(3,IFFAT) - C(3,JAT)
            PBCX  = XBOX*ANINT(X*ONEXBOX)
            PBCY  = YBOX*ANINT(Y*ONEYBOX)
            PBCZ  = ZBOX*ANINT(Z*ONEZBOX)
            X     = X - PBCX
            Y     = Y - PBCY
            Z     = Z - PBCZ
            CXGRP = X
            CYGRP = Y
            CZGRP = Z
            R2    = X*X+Y*Y+Z*Z+DFTBMM   ! FOR BOTH QM AND DFTB
            IF(R2.GT.SWRB2) GOTO 210
            IF(R2.LT.1.0D-10) GOTO 210
            R     = SQRT(R2)
            ONER  = ONE/R
            ONER2 = ONER*ONER
            CALL SHIFT(R2,R,ONER,X,Y,Z)
            END IF
C
            EPAIR = QI*ZAN(JAT)*ONER
            ENUCCH= ENUCCH + EPAIR*SWF
            DUM   = -EPAIR*ONER2*SWF
            DEX   = DUM*X
            DEY   = DUM*Y
            DEZ   = DUM*Z
            FFGRD(1,IFFAT) = FFGRD(1,IFFAT) + DEX
            FFGRD(2,IFFAT) = FFGRD(2,IFFAT) + DEY
            FFGRD(3,IFFAT) = FFGRD(3,IFFAT) + DEZ
            FFGRD(1,IFFAT) = FFGRD(1,IFFAT) + EPAIR*SWFDX
            FFGRD(2,IFFAT) = FFGRD(2,IFFAT) + EPAIR*SWFDY
            FFGRD(3,IFFAT) = FFGRD(3,IFFAT) + EPAIR*SWFDZ
            DETMP(1,JAT )  = DETMP(1,JAT )  - DEX
            DETMP(2,JAT )  = DETMP(2,JAT )  - DEY
            DETMP(3,JAT )  = DETMP(3,JAT )  - DEZ
            IF(ISWITCH.EQ.2)THEN    ! ONLY ATOM-ATOM SHIFTING NEEDS THIS
            DETMP(1,JAT )  = DETMP(1,JAT )  - EPAIR*SWFDX
            DETMP(2,JAT )  = DETMP(2,JAT )  - EPAIR*SWFDY
            DETMP(3,JAT )  = DETMP(3,JAT )  - EPAIR*SWFDZ
            END IF
            VIR(1)         = VIR(1) + DEX*X + EPAIR*SWFDX*CXGRP
            VIR(2)         = VIR(2) + DEY*Y + EPAIR*SWFDY*CYGRP
            VIR(3)         = VIR(3) + DEZ*Z + EPAIR*SWFDZ*CZGRP
            IYES = 0
            DO KFIX=1,NFIXMM
               IF(IFFAT.EQ.IFIXMM(KFIX)) IYES = 1
            ENDDO
            IF(IYES.EQ.1.AND.NFIXQM.GT.0) THEN
               VIR(1)      = VIR(1) - DEX*X - EPAIR*SWFDX*CXGRP
               VIR(2)      = VIR(2) - DEY*Y - EPAIR*SWFDY*CYGRP
               VIR(3)      = VIR(3) - DEZ*Z - EPAIR*SWFDZ*CZGRP
            END IF
  210    CONTINUE
  200 CONTINUE
C
      RETURN
      END
C*MODULE QUANPOC  *DECK QMMMCHGCHG
!>
!> @brief    QM charge -MM charge force
!>
!> @author   Hui Li
!>           - Oct 2016
!>
!> @details  mutual forces are calculated
!>
      SUBROUTINE QMMMCHGCHG(CORD,FFGRD,DETMP,CHARG,NONLSTQ)
      use mx_limits, only: mxatm
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      LOGICAL GOPARR,DSKWRK,MASWRK
C
      PARAMETER (ZERO=0.0D+00)
      PARAMETER (ONE=1.0D+00)
C
      DIMENSION CORD(3,*),FFGRD(3,*),CHARG(*),DETMP(3,*),
     *          NONLSTQ(*)
C
      COMMON /FFENGY/ EN12,EN123,EN123R4,EN123B4,EN234W1,ENCHAR,ENLJR,
     *                ENLJD,ENPOL,XENPOL,ENRXN,XENRXN,ENRXNPOL,ENRXNR,
     *                EN12312,ENQUANP(30),
     *                ENBIAS,ENCENT,ENUCCH,ENCMAP,ENPOT,ENKIN,ENTOT
      COMMON /FFPARM/ NFFAT,NBOND,NANGL,NDIHR,NDIHB,NCMAP,NWAGG,
     *                N1213J,N14J,NLKQMM,IDOCHG,IDOPOL,IDOLJ,IDOCMAP
      COMMON /FFMDPA/ DT,DT2,TEMP0,PRES0,POLTOL,VIR(3),PMEAN,VOLAV,
     *                ENPAV,ENKAV,TEMPAV,BERENDT,BERENDP,VELMAX,
     *                PMEANX,PMEANY,PMEANZ,NSTEP,KMASTER,KOUTACT(2),
     *                IHESS,INTALG,ITSTAT,IPSTAT,JOUT,KOUT,LOUT
      COMMON /FFMPT3/ NACTMM,LACTMM(2020),LFFDIPOLD,JUMBUP,
     *                NACTQM,LACTQM(2020),LFFOLDC,LFFQMVELSV,MMHESS,
     *                LFFQMCHG,LFFQMCHGB,ISWAP,R2SWAP,DFTBMM
      COMMON /FFPBSW/ XBOX,YBOX,ZBOX,SWF,SWFDX,SWFDY,SWFDZ,
     *                SWRA,ONESWRA,SWRA2,ONESWRA2,
     *                SWRB,ONESWRB,SWRB2,ONESWRB2,ONESWRB4,
     *                SWFDUM3,SWFDUM4,SWFDUM5,
     *                SWRAQ,ONESWRAQ,SWRAQ2,ONESWRAQ2,
     *                SWRBQ,ONESWRBQ,SWRBQ2,ONESWRBQ2,ONESWRBQ4,
     *                SWFDUM3Q,SWFDUM4Q,SWFDUM5Q,
     *                QMSIZE,QMCX,QMCY,QMCZ,QMCXSV,QMCYSV,QMCZSV,
     *                CENTX,CENTY,CENTZ,BUFWID1,BUFWID2,RDAMP,
     *                EFIELDX,EFIELDY,EFIELDZ,QMCXSV2,QMCYSV2,QMCZSV2,
     *                EPS1RB,EPS1RB3,ONEXBOX,ONEYBOX,ONEZBOX,
     *                LQMCT,MXLIST1,NTODO,NTODOSV,NTODOQ,
     *                ISWITCH,ISHIFT,IPOLSHF,
     *                LFFLSTCELL,LFFCORDSV,
     *                LFFPOLSV,LFFCORDSV2,LFFNONLS2,LFFCORDSVQ,
     *                LFFMVFASTS2,LFFMVFASTS3,LFFMVFASTS4,
     *                LFFMVFASTL2,LFFMVFASTL3,LFFMVFASTL4,
     *                MXCHECK,MXLIST2,NTODO2,NTODO2SV
      COMMON /FFRMSD/ DIPT(8),TIMGYRA,TIMRALL,
     *                NATPDB,NGYRA,NDIEL,NRALL,
     *                LFFRALL0,NRMSD,LFFRMSD0,KFREEAB(201),
     *                NRIJMM,IJRMM(2,100),NRIJQM,IJRQM(2,100),
     *                NAIJKMM,IJKMM(3,100),NAIJKQM,IJKQM(3,100),
     *                NFIXMM,IFIXMM(200),NFIXQM,IFIXQM(200)
      COMMON /FMCOM / XX(1)
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
C
C     HUI LI, OCT 2016, LINCOLN
C
C     MM CHARGE - QM CHARGE INTERACTION
C
      CALL VCLR(DETMP,1,3*MXATM)
      ENUCCH =ZERO
      IF(IDOCHG.EQ.0) RETURN
C
      IPCOUNT = ME - 1
      DO 200 IIQ = 1, NTODOQ
         IFFAT = NONLSTQ(IIQ)
         QI = CHARG(IFFAT)
         IF(QI.EQ.ZERO) GOTO 200
         IF(GOPARR) THEN
            IPCOUNT = IPCOUNT + 1
            IF (MOD(IPCOUNT,NPROC).NE.0) GOTO 200
         END IF
         DO 210 JAT=1,NAT
            ZANJAT=XX(LFFQMCHG+JAT-1)
            IF(ZANJAT.EQ.ZERO) GOTO 210
C           - ALWAYS USE SHIFTING FUNCTION
            X     = CORD(1,IFFAT) - C(1,JAT)
            Y     = CORD(2,IFFAT) - C(2,JAT)
            Z     = CORD(3,IFFAT) - C(3,JAT)
            PBCX  = XBOX*ANINT(X*ONEXBOX)
            PBCY  = YBOX*ANINT(Y*ONEYBOX)
            PBCZ  = ZBOX*ANINT(Z*ONEZBOX)
            X     = X - PBCX
            Y     = Y - PBCY
            Z     = Z - PBCZ
            R2    = X*X+Y*Y+Z*Z
            IF(R2.GT.SWRB2) GOTO 210
            IF(R2.LT.1.0D-10) GOTO 210
            R     = SQRT(R2)
            ONER  = ONE/R
            ONER2 = ONER*ONER
            CALL SHIFT(R2,R,ONER,X,Y,Z)
C
            EPAIR = QI*ZANJAT*ONER
            ENUCCH= ENUCCH + EPAIR*SWF
            DUM   = -EPAIR*ONER2*SWF
            DEX   = DUM*X
            DEY   = DUM*Y
            DEZ   = DUM*Z
            FFGRD(1,IFFAT) = FFGRD(1,IFFAT) + DEX
            FFGRD(2,IFFAT) = FFGRD(2,IFFAT) + DEY
            FFGRD(3,IFFAT) = FFGRD(3,IFFAT) + DEZ
            FFGRD(1,IFFAT) = FFGRD(1,IFFAT) + EPAIR*SWFDX
            FFGRD(2,IFFAT) = FFGRD(2,IFFAT) + EPAIR*SWFDY
            FFGRD(3,IFFAT) = FFGRD(3,IFFAT) + EPAIR*SWFDZ
            DETMP(1,JAT )  = DETMP(1,JAT )  - DEX
            DETMP(2,JAT )  = DETMP(2,JAT )  - DEY
            DETMP(3,JAT )  = DETMP(3,JAT )  - DEZ
C           - QM ATOMS NEED SHIFTING SWFDX TERM
            DETMP(1,JAT )  = DETMP(1,JAT )  - EPAIR*SWFDX
            DETMP(2,JAT )  = DETMP(2,JAT )  - EPAIR*SWFDY
            DETMP(3,JAT )  = DETMP(3,JAT )  - EPAIR*SWFDZ
            VIR(1)         = VIR(1) + DEX*X + EPAIR*SWFDX*X
            VIR(2)         = VIR(2) + DEY*Y + EPAIR*SWFDY*Y
            VIR(3)         = VIR(3) + DEZ*Z + EPAIR*SWFDZ*Z
            IYES = 0
            DO KFIX=1,NFIXMM
               IF(IFFAT.EQ.IFIXMM(KFIX)) IYES = 1
            ENDDO
            IF(IYES.EQ.1.AND.NFIXQM.GT.0) THEN
               VIR(1)      = VIR(1) - DEX*X - EPAIR*SWFDX*X
               VIR(2)      = VIR(2) - DEY*Y - EPAIR*SWFDY*Y
               VIR(3)      = VIR(3) - DEZ*Z - EPAIR*SWFDZ*Z
            END IF
  210    CONTINUE
  200 CONTINUE
C
      RETURN
      END
C*MODULE QUANPOC  *DECK QMMMREPINT
      SUBROUTINE QMMMREPINT(QQ,CORD,CLPR,ZLPR,NLPR,NONLSTQ)
      use mx_limits, only: mxsh,mxgtot,mxatm
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      LOGICAL IANDJ,NORM,DOUBLE
      LOGICAL GOPARR,MASWRK,DSKWRK
C
      DIMENSION QQ(*),CORD(3,*),CLPR(4,*),ZLPR(4,*),NLPR(*),
     *          NONLSTQ(*)
C
      DIMENSION RINT(100)
      DIMENSION DIJ(100),XIN(432),YIN(432),ZIN(432)
      DIMENSION GIJ(100)
      DIMENSION IX(20),IY(20),IZ(20),JX(20),JY(20),JZ(20)
      DIMENSION IJX(100),IJY(100),IJZ(100)
C
C
      COMMON /CSSTV / CX,CY,CZ
      COMMON /FFPBSW/ XBOX,YBOX,ZBOX,SWF,SWFDX,SWFDY,SWFDZ,
     *                SWRA,ONESWRA,SWRA2,ONESWRA2,
     *                SWRB,ONESWRB,SWRB2,ONESWRB2,ONESWRB4,
     *                SWFDUM3,SWFDUM4,SWFDUM5,
     *                SWRAQ,ONESWRAQ,SWRAQ2,ONESWRAQ2,
     *                SWRBQ,ONESWRBQ,SWRBQ2,ONESWRBQ2,ONESWRBQ4,
     *                SWFDUM3Q,SWFDUM4Q,SWFDUM5Q,
     *                QMSIZE,QMCX,QMCY,QMCZ,QMCXSV,QMCYSV,QMCZSV,
     *                CENTX,CENTY,CENTZ,BUFWID1,BUFWID2,RDAMP,
     *                EFIELDX,EFIELDY,EFIELDZ,QMCXSV2,QMCYSV2,QMCZSV2,
     *                EPS1RB,EPS1RB3,ONEXBOX,ONEYBOX,ONEZBOX,
     *                LQMCT,MXLIST1,NTODO,NTODOSV,NTODOQ,
     *                ISWITCH,ISHIFT,IPOLSHF,
     *                LFFLSTCELL,LFFCORDSV,
     *                LFFPOLSV,LFFCORDSV2,LFFNONLS2,LFFCORDSVQ,
     *                LFFMVFASTS2,LFFMVFASTS3,LFFMVFASTS4,
     *                LFFMVFASTL2,LFFMVFASTL3,LFFMVFASTL4,
     *                MXCHECK,MXLIST2,NTODO2,NTODO2SV
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /IOFILE/ IR,IW,IP,IJK,IJKT,IDAF,NAV,IODA(950)
      COMMON /NSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),
     *                CF(MXGTOT),CG(MXGTOT),CH(MXGTOT),CI(MXGTOT),
     *                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),
     *                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
      COMMON /OUTPUT/ NPRINT,ITOL,ICUT,NORMF,NORMP,NOPK
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /ROOT  / XX,U(13),W(13),NROOTS
      COMMON /STV   / XINT,YINT,ZINT,T,X0,Y0,Z0,XI,YI,ZI,XJ,YJ,ZJ,NI,NJ
C
      DATA JX / 0, 1, 0, 0, 2, 0, 0, 1, 1, 0,
     1          3, 0, 0, 2, 2, 1, 0, 1, 0, 1/
      DATA JY / 0, 0, 1, 0, 0, 2, 0, 1, 0, 1,
     1          0, 3, 0, 1, 0, 2, 2, 0, 1, 1/
      DATA JZ / 0, 0, 0, 1, 0, 0, 2, 0, 1, 1,
     1          0, 0, 3, 0, 1, 0, 1, 2, 2, 1/
      DATA IX / 1, 5, 1, 1, 9, 1, 1, 5, 5, 1,
     1         13, 1, 1, 9, 9, 5, 1, 5, 1, 5/
      DATA IY / 1, 1, 5, 1, 1, 9, 1, 5, 1, 5,
     1          1,13, 1, 5, 1, 9, 9, 1, 5, 5/
      DATA IZ / 1, 1, 1, 5, 1, 1, 9, 1, 5, 5,
     1          1, 1,13, 1, 5, 1, 5, 9, 9, 5/
C
      DATA PI212 /1.1283791670955D+00/
      DATA SQRT3 /1.73205080756888D+00/
      DATA SQRT5 /2.23606797749979D+00/
      DATA ZERO,ONE/0.0D+00,1.0D+00/
      DATA RLN10 /2.30258D+00/
C
C     HUI LI, JAN 2011, LINCOLN
C
C     --- CHARGE-REPULSIVE POTENTIAL INTEGRALS FOR POWERS 0 AND -1 OF R.
C
      TOL = RLN10*ITOL
      NORM = NORMF .NE. 1 .OR. NORMP .NE. 1
C
      IPCOUNT = ME - 1
C
C     ----- ISHELL
C
      DO 600 II = 1,NSHELL
      IAT= KATOM(II)
      XI = C(1,IAT)
      YI = C(2,IAT)
      ZI = C(3,IAT)
      I1 = KSTART(II)
      I2 = I1+KNG(II)-1
      LIT = KTYPE(II)
      MINI = KMIN(II)
      MAXI = KMAX(II)
      LOCI = KLOC(II)-MINI
C
C     ----- JSHELL
C
      DO 580 JJ = 1,II
      IF (GOPARR) THEN
         IPCOUNT = IPCOUNT + 1
         IF (MOD(IPCOUNT,NPROC).NE.0) GOTO 580
      END IF
      JAT= KATOM(JJ)
      XJ = C(1,JAT)
      YJ = C(2,JAT)
      ZJ = C(3,JAT)
      J1 = KSTART(JJ)
      J2 = J1+KNG(JJ)-1
      LJT = KTYPE(JJ)
      MINJ = KMIN(JJ)
      MAXJ = KMAX(JJ)
      LOCJ = KLOC(JJ)-MINJ
      NROOTS = (LIT+LJT-2)/2 + 1
      RR = (XI-XJ)**2+(YI-YJ)**2+(ZI-ZJ)**2
      IANDJ = II .EQ. JJ
C
C     ----- PREPARE INDICES FOR PAIRS OF (I,J) FUNCTIONS
C
      IJ = 0
      MAX = MAXJ
      DO 50 I = MINI,MAXI
      NX = IX(I)
      NY = IY(I)
      NZ = IZ(I)
      IF (IANDJ) MAX = I
      DO 50 J = MINJ,MAX
      IJ = IJ+1
      IJX(IJ) = NX+JX(J)
      IJY(IJ) = NY+JY(J)
      IJZ(IJ) = NZ+JZ(J)
  50  CONTINUE
      DO 60 I=1,IJ
  60  RINT(I) = ZERO
C
C     ----- I PRIMITIVE
C
      JGMAX = J2
      DO 520 IG = I1,I2
      AI  = EX(IG)
      ARRI= AI*RR
      AXI = AI*XI
      AYI = AI*YI
      AZI = AI*ZI
      CSI = CS(IG)
      CPI = CP(IG)
      CDI = CD(IG)
      CFI = CF(IG)
C
C     ----- J PRIMITIVE
C
      IF (IANDJ) JGMAX = IG
      DO 500 JG = J1,JGMAX
      AJ = EX(JG)
      AA = AI+AJ
      AA1 = ONE/AA
      DUM = AJ*ARRI*AA1
      IF (DUM .GT. TOL) GOTO 500
      FAC =  EXP(-DUM)
      CSJ = CS(JG)
      CPJ = CP(JG)
      CDJ = CD(JG)
      CFJ = CF(JG)
      AX = (AXI+AJ*XJ)*AA1
      AY = (AYI+AJ*YJ)*AA1
      AZ = (AZI+AJ*ZJ)*AA1
C
C     ----- DENSITY FACTOR
C
      DOUBLE=IANDJ.AND.IG.NE.JG
      MAX = MAXJ
      NN = 0
      DO 170 I=MINI,MAXI
      GOTO (70,80,110,110,90,110,110,100,110,110,
     *       102,110,110,104,110,110,110,110,110,106),I
   70 DUM1=CSI*FAC
      GOTO 110
   80 DUM1=CPI*FAC
      GOTO 110
   90 DUM1=CDI*FAC
      GOTO 110
  100 IF(NORM) DUM1=DUM1*SQRT3
      GOTO 110
  102 DUM1=CFI*FAC
      GOTO 110
  104 IF(NORM) DUM1=DUM1*SQRT5
      GOTO 110
  106 IF(NORM) DUM1=DUM1*SQRT3
  110 IF(IANDJ) MAX=I
      DO 170 J=MINJ,MAX
      GOTO (125,130,160,160,140,160,160,150,160,160,
     1       152,160,160,154,160,160,160,160,160,156),J
  125 DUM2=DUM1*CSJ
      IF(.NOT.DOUBLE) GOTO 160
      IF(I.GT.1) GOTO 126
      DUM2=DUM2+DUM2
      GOTO 160
  126 DUM2=DUM2+CSI*CPJ*FAC
      GOTO 160
  130 DUM2=DUM1*CPJ
      IF(DOUBLE) DUM2=DUM2+DUM2
      GOTO 160
  140 DUM2=DUM1*CDJ
      IF(DOUBLE) DUM2=DUM2+DUM2
      GOTO 160
  150 IF(NORM) DUM2=DUM2*SQRT3
      GOTO 160
  152 DUM2=DUM1*CFJ
      IF(DOUBLE) DUM2=DUM2+DUM2
      GOTO 160
  154 IF(NORM) DUM2=DUM2*SQRT5
      GOTO 160
  156 IF(NORM) DUM2=DUM2*SQRT3
  160 NN=NN+1
  170 DIJ(NN)=DUM2
      AAX = AA*AX
      AAY = AA*AY
      AAZ = AA*AZ
C
      DO 481 IIQ = 1, NTODOQ
      IFFAT = NONLSTQ(IIQ)
      IF((CLPR(1,IFFAT)+CLPR(2,IFFAT)+
     *    CLPR(3,IFFAT)+CLPR(4,IFFAT)).LT.1.0D-06) GOTO 481
      X     = CORD(1,IFFAT) - AX
      Y     = CORD(2,IFFAT) - AY
      Z     = CORD(3,IFFAT) - AZ
      PBCX  = XBOX*ANINT(X*ONEXBOX)
      PBCY  = YBOX*ANINT(Y*ONEYBOX)
      PBCZ  = ZBOX*ANINT(Z*ONEZBOX)
      CX    = CORD(1,IFFAT) - PBCX
      CY    = CORD(2,IFFAT) - PBCY
      CZ    = CORD(3,IFFAT) - PBCZ
      PCSQ = (AX-CX)**2+(AY-CY)**2+(AZ-CZ)**2
      IF(PCSQ.GT.SWRB2) GOTO 481     ! CUTOFF, NO SWF/SHIFT
      DO 479 LTERM=1,4
      ALFA = ZLPR(LTERM,IFFAT)
      BETA = CLPR(LTERM,IFFAT)
      IF(BETA.EQ.ZERO) GOTO 479
      PREI = EXP(-AA*ALFA*PCSQ/(AA+ALFA))
      IF(NLPR(IFFAT).EQ.2) THEN
       TT = ONE/(AA+ALFA)
       T  =  SQRT(TT)
       X0 = (AAX+ALFA*CX)*TT
       Y0 = (AAY+ALFA*CY)*TT
       Z0 = (AAZ+ALFA*CZ)*TT
       IN = -4
       DO 206 I = 1,LIT
       IN = IN+4
       NI = I
       DO 206 J = 1,LJT
       JN = IN+J
       NJ = J
       CALL STVINT
       XIN(JN   ) = XINT*T
       YIN(JN   ) = YINT*T
       ZIN(JN   ) = ZINT*T
  206  CONTINUE
       DO 266 I = 1,IJ
       NX    = IJX(I)
       NY    = IJY(I)
       NZ    = IJZ(I)
       RINT(I)=RINT(I)+DIJ(I)*PREI*BETA*XIN(NX)*YIN(NY)*ZIN(NZ)
 266   CONTINUE
      ELSE
C      ONE OVER R.
       DUM = PI212/(AA+ALFA)
       DO 385 I = 1,IJ
       GIJ(I) = DIJ(I)*DUM
  385  CONTINUE
       XX=AA*AA*PCSQ/(AA+ALFA)
       IF(NROOTS .LE. 3) CALL RT123
       IF(NROOTS .EQ. 4) CALL ROOT4
       IF(NROOTS.EQ.5) CALL ROOT5
       IF(NROOTS.EQ.6) CALL ROOT6
       IF(NROOTS.GE.7) THEN
          IF (MASWRK) WRITE(IW,9008)
          CALL ABRT
       END IF
       MM = 0
       DO 425 K = 1,NROOTS
       UU = (ALFA+AA)*U(K)
       WW = W(K)
       TT = ONE/(AA+UU+ALFA)
       T  =  SQRT(TT)
       X0 = (AAX+(UU+ALFA)*CX)*TT
       Y0 = (AAY+(UU+ALFA)*CY)*TT
       Z0 = (AAZ+(UU+ALFA)*CZ)*TT
       IN = -4+MM
       DO 405 I = 1,LIT
       IN = IN+4
       NI = I
       DO 405 J = 1,LJT
       JN = IN+J
       NJ = J
       CALL STVINT
       XIN(JN   ) = XINT
       YIN(JN   ) = YINT
       ZIN(JN   ) = ZINT*WW
  405  CONTINUE
  425  MM = MM+16
       DO 465 I = 1,IJ
       NX    = IJX(I)
       NY    = IJY(I)
       NZ    = IJZ(I)
       DUM = ZERO
       MM    = 0
       DO 445 K = 1,NROOTS
       DUM= DUM+XIN(NX+MM)*YIN(NY+MM)*ZIN(NZ+MM)
  445  MM  = MM+16
       RINT(I)=RINT(I)+GIJ(I)*PREI*BETA*DUM
 465   CONTINUE
      END IF
 479  CONTINUE
 481  CONTINUE
 500  CONTINUE
 520  CONTINUE
C
C     ----- END OF *PRIMITIVE* LOOPS -----
C     ----- SET QQ MATRIX
C
      MAX=MAXJ
      NN=0
      DO 550 I=MINI,MAXI
      LI=LOCI+I
      IN = (LI*(LI-1))/2
      IF(IANDJ) MAX=I
      DO 550 J=MINJ,MAX
      LJ=LOCJ+J
      JN=LJ+IN
      NN=NN+1
      QQ(JN)=QQ(JN) + RINT(NN)
  550 CONTINUE
  580 CONTINUE
  600 CONTINUE
C
      RETURN
C
 9008 FORMAT(/' NUMBER OF POLYNOMIAL ROOTS NEEDED (NROOTS) IS GREATER',
     *        ' THAN 6 IN REPINT.  CALL A PROGRAMMER/QUANTUM CHEMIST.')
      END
C*MODULE QUANPOC  *DECK QMMMREPFMM
!>
!> @brief    QM electron - MM ex-rep force
!>
!> @author   Hui Li
!>           - Jan 2011
!>
!> @details  forces on MM ex-rep points due to QM electrons
!>
      SUBROUTINE QMMMREPFMM(DM,CHGINT,CORD,CLPR,ZLPR,NLPR,FFGRD,
     *                      NONLSTQ)
      use mx_limits, only: mxsh,mxgtot,mxatm
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C
      LOGICAL NORM,GOPARR,DSKWRK,MASWRK,NXT
C
      DIMENSION DM(*),CHGINT(*),CORD(3,*),CLPR(4,*),ZLPR(4,*),
     *          NLPR(*),FFGRD(3,*),NONLSTQ(*)
C
      DIMENSION DIJ(100),XIN(432),YIN(432),ZIN(432)
      DIMENSION GIJ(100)
      DIMENSION IX(20),IY(20),IZ(20),JX(20),JY(20),JZ(20)
      DIMENSION IJX(100),IJY(100),IJZ(100)
C
      COMMON /CSSTV / CX,CY,CZ
      COMMON /FFMDPA/ DT,DT2,TEMP0,PRES0,POLTOL,VIR(3),PMEAN,VOLAV,
     *                ENPAV,ENKAV,TEMPAV,BERENDT,BERENDP,VELMAX,
     *                PMEANX,PMEANY,PMEANZ,NSTEP,KMASTER,KOUTACT(2),
     *                IHESS,INTALG,ITSTAT,IPSTAT,JOUT,KOUT,LOUT
      COMMON /FFPBSW/ XBOX,YBOX,ZBOX,SWF,SWFDX,SWFDY,SWFDZ,
     *                SWRA,ONESWRA,SWRA2,ONESWRA2,
     *                SWRB,ONESWRB,SWRB2,ONESWRB2,ONESWRB4,
     *                SWFDUM3,SWFDUM4,SWFDUM5,
     *                SWRAQ,ONESWRAQ,SWRAQ2,ONESWRAQ2,
     *                SWRBQ,ONESWRBQ,SWRBQ2,ONESWRBQ2,ONESWRBQ4,
     *                SWFDUM3Q,SWFDUM4Q,SWFDUM5Q,
     *                QMSIZE,QMCX,QMCY,QMCZ,QMCXSV,QMCYSV,QMCZSV,
     *                CENTX,CENTY,CENTZ,BUFWID1,BUFWID2,RDAMP,
     *                EFIELDX,EFIELDY,EFIELDZ,QMCXSV2,QMCYSV2,QMCZSV2,
     *                EPS1RB,EPS1RB3,ONEXBOX,ONEYBOX,ONEZBOX,
     *                LQMCT,MXLIST1,NTODO,NTODOSV,NTODOQ,
     *                ISWITCH,ISHIFT,IPOLSHF,
     *                LFFLSTCELL,LFFCORDSV,
     *                LFFPOLSV,LFFCORDSV2,LFFNONLS2,LFFCORDSVQ,
     *                LFFMVFASTS2,LFFMVFASTS3,LFFMVFASTS4,
     *                LFFMVFASTL2,LFFMVFASTL3,LFFMVFASTL4,
     *                MXCHECK,MXLIST2,NTODO2,NTODO2SV
      COMMON /FFRMSD/ DIPT(8),TIMGYRA,TIMRALL,
     *                NATPDB,NGYRA,NDIEL,NRALL,
     *                LFFRALL0,NRMSD,LFFRMSD0,KFREEAB(201),
     *                NRIJMM,IJRMM(2,100),NRIJQM,IJRQM(2,100),
     *                NAIJKMM,IJKMM(3,100),NAIJKQM,IJKQM(3,100),
     *                NFIXMM,IFIXMM(200),NFIXQM,IFIXQM(200)
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /IOFILE/ IR,IW,IP,IJK,IJKT,IDAF,NAV,IODA(950)
      COMMON /NSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),
     *                CF(MXGTOT),CG(MXGTOT),CH(MXGTOT),CI(MXGTOT),
     *                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),
     *                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
      COMMON /OUTPUT/ NPRINT,ITOL,ICUT,NORMF,NORMP,NOPK
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /ROOT  / XX,U(13),W(13),NROOTS
      COMMON /STV   / XINT,YINT,ZINT,T,X0,Y0,Z0,XI,YI,ZI,XJ,YJ,ZJ,NI,NJ
C
      DATA JX / 0, 1, 0, 0, 2, 0, 0, 1, 1, 0,
     1          3, 0, 0, 2, 2, 1, 0, 1, 0, 1/
      DATA JY / 0, 0, 1, 0, 0, 2, 0, 1, 0, 1,
     1          0, 3, 0, 1, 0, 2, 2, 0, 1, 1/
      DATA JZ / 0, 0, 0, 1, 0, 0, 2, 0, 1, 1,
     1          0, 0, 3, 0, 1, 0, 1, 2, 2, 1/
      DATA IX / 1, 5, 1, 1, 9, 1, 1, 5, 5, 1,
     1         13, 1, 1, 9, 9, 5, 1, 5, 1, 5/
      DATA IY / 1, 1, 5, 1, 1, 9, 1, 5, 1, 5,
     1          1,13, 1, 5, 1, 9, 9, 1, 5, 5/
      DATA IZ / 1, 1, 1, 5, 1, 1, 9, 1, 5, 5,
     1          1, 1,13, 1, 5, 1, 5, 9, 9, 5/
C
      DATA ZERO,ONE,TWO/0.0D+00,1.0D+00,2.0D+00/
      DATA RLN10 /2.30258D+00/
      DATA PI212 /1.1283791670955D+00/
      DATA SQRT3 /1.73205080756888D+00/
      DATA SQRT5 /2.23606797749979D+00/
C
C     HUI LI, JAN 2011, LINCOLN
C
      NXT = IBTYP.EQ.1
      IPCOUNT = ME - 1
      NEXT = -1
      LCNT = -1
      TOL = RLN10*ITOL
      NORM = NORMF .NE. 1 .OR. NORMP .NE. 1
C
      DO 100 IIQ = 1, NTODOQ
      IFFAT = NONLSTQ(IIQ)
      IF((CLPR(1,IFFAT)+CLPR(2,IFFAT)+
     *    CLPR(3,IFFAT)+CLPR(4,IFFAT)).LT.1.0D-06) GOTO 100
C
C     ----- ISHELL
C
      DO 600 II = 1,NSHELL
      IF(NXT .AND. GOPARR) THEN
         LCNT = LCNT + 1
         IF(LCNT.GT.NEXT) CALL DDI_DLBNEXT(NEXT)
         IF(NEXT.NE.LCNT) GOTO 600
      END IF
      IAT= KATOM(II)
      XI = C(1,IAT)
      YI = C(2,IAT)
      ZI = C(3,IAT)
      I1 = KSTART(II)
      I2 = I1+KNG(II)-1
      LIT = KTYPE(II)
      MINI = KMIN(II)
      MAXI = KMAX(II)
      LOCI = KLOC(II)-MINI
C
C     ----- JSHELL
C
      DO 580 JJ = 1,NSHELL
      IF((.NOT.NXT) .AND. GOPARR) THEN
         IPCOUNT = IPCOUNT + 1
         IF(MOD(IPCOUNT,NPROC).NE.0) GOTO 580
      END IF
      JAT= KATOM(JJ)
      XJ = C(1,JAT)
      YJ = C(2,JAT)
      ZJ = C(3,JAT)
      J1 = KSTART(JJ)
      J2 = J1+KNG(JJ)-1
      LJT = KTYPE(JJ)
      MINJ = KMIN(JJ)
      MAXJ = KMAX(JJ)
      LOCJ = KLOC(JJ)-MINJ
      NROOTS = (LIT+LJT-2)/2 + 1
      RR = (XI-XJ)**2+(YI-YJ)**2+(ZI-ZJ)**2
C
C     ----- PREPARE INDICES FOR PAIRS OF (I,J) FUNCTIONS
C
      IJ = 0
      MAX = MAXJ
      DO 50 I = MINI,MAXI
      NX = IX(I)
      NY = IY(I)
      NZ = IZ(I)
      DO 50 J = MINJ,MAX
      IJ = IJ+1
      IJX(IJ) = NX+JX(J)
      IJY(IJ) = NY+JY(J)
      IJZ(IJ) = NZ+JZ(J)
  50  CONTINUE
      CALL VCLR(CHGINT,1,3*IJ*4)
C
C     ----- I PRIMITIVE
C
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
C
C     ----- J PRIMITIVE
C
      DO 500 JG = J1,J2
      AJ = EX(JG)
      AA = AI+AJ
      AA1 = ONE/AA
      DUM = AJ*ARRI*AA1
      IF (DUM .GT. TOL) GOTO 500
      FAC =  EXP(-DUM)
      CSJ = CS(JG)
      CPJ = CP(JG)
      CDJ = CD(JG)
      CFJ = CF(JG)
      AX = (AXI+AJ*XJ)*AA1
      AY = (AYI+AJ*YJ)*AA1
      AZ = (AZI+AJ*ZJ)*AA1
C
C     ----- DENSITY FACTOR
C
      MAX = MAXJ
      NN = 0
      DUM1=ZERO
      DUM2=DUM1
      DO 170 I=MINI,MAXI
      GOTO (70,80,110,110,90,110,110,95,110,110,
     1       102,110,110,104,110,110,110,110,110,108),I
  70  DUM1=CSI*FAC
      GOTO 110
  80  DUM1=CPI*FAC
      GOTO 110
  90  DUM1=CDI*FAC
      GOTO 110
  95  IF(NORM) DUM1=DUM1*SQRT3
      GOTO 110
 102  DUM1=CFI*FAC
      GOTO 110
 104  DUM1 = DUM1 *SQRT5
      GOTO 110
 108  DUM1 = DUM1 * SQRT3
 110  CONTINUE
      DO 170 J=MINJ,MAX
      GOTO (125,130,160,160,140,160,160,150,160,160,
     1       152,160,160,154,160,160,160,160,160,156),J
  125 DUM2=DUM1*CSJ
      GOTO 160
  130 DUM2=DUM1*CPJ
      GOTO 160
  140 DUM2=DUM1*CDJ
      GOTO 160
  150 IF(NORM) DUM2=DUM2*SQRT3
      GOTO 160
  152 DUM2 = DUM1 * CFJ
      GOTO 160
  154 DUM2 = DUM2 *SQRT5
      GOTO 160
  156 DUM2 = DUM2 * SQRT3
  160 NN=NN+1
  170 DIJ(NN)=DUM2
C
      AAX = AA*AX
      AAY = AA*AY
      AAZ = AA*AZ
C
      ICC=1
C
      X     = CORD(1,IFFAT) - AX
      Y     = CORD(2,IFFAT) - AY
      Z     = CORD(3,IFFAT) - AZ
      PBCX  = XBOX*ANINT(X*ONEXBOX)
      PBCY  = YBOX*ANINT(Y*ONEYBOX)
      PBCZ  = ZBOX*ANINT(Z*ONEZBOX)
      CX    = CORD(1,IFFAT) - PBCX
      CY    = CORD(2,IFFAT) - PBCY
      CZ    = CORD(3,IFFAT) - PBCZ
      PCSQ = (AX-CX)**2+(AY-CY)**2+(AZ-CZ)**2
      IF(PCSQ.GT.SWRB2) GOTO 500    ! CUTOFF, NO SWF/SHIFT
      DO 480 LTERM=1,4
      ALFA = ZLPR(LTERM,IFFAT)
      BETA = CLPR(LTERM,IFFAT)
      IF(BETA.EQ.ZERO) GOTO 480
      PREI = EXP(-AA*ALFA*PCSQ/(AA+ALFA))
      IF(NLPR(IFFAT).EQ.2) THEN
C      R TO THE ZERO POWER.
       TT = ONE/(AA+ALFA)
       T  =  SQRT(TT)
       X0 = (AAX+ALFA*CX)*TT
       Y0 = (AAY+ALFA*CY)*TT
       Z0 = (AAZ+ALFA*CZ)*TT
       IN = -4
       DO 206 I = 1,LIT
       IN = IN+4
       NI = I
       DO 206 J = 1,LJT
       JN = IN+J
       NJ = J
       CALL STVINT
       XIN(JN   ) = XINT*T
       YIN(JN   ) = YINT*T
       ZIN(JN   ) = ZINT*T
      CALL POLXYZ
      XIN(JN+125) = XINT*T
      YIN(JN+125) = YINT*T
      ZIN(JN+125) = ZINT*T
  206  CONTINUE
       DO 266 I = 1,IJ
       NX    = IJX(I)
       NY    = IJY(I)
       NZ    = IJZ(I)
      DUMX= XIN(NX+125)*YIN(NY    )*ZIN(NZ    )
      DUMY= XIN(NX    )*YIN(NY+125)*ZIN(NZ    )
      DUMZ= XIN(NX    )*YIN(NY    )*ZIN(NZ+125)
      DUM = TWO*DIJ(I)*ALFA*PREI*BETA
      CHGINT(ICC  )=CHGINT(ICC  )+DUM*DUMX
      CHGINT(ICC+1)=CHGINT(ICC+1)+DUM*DUMY
      CHGINT(ICC+2)=CHGINT(ICC+2)+DUM*DUMZ
      ICC=ICC+3
 266   CONTINUE
      ELSE
C      ONE OVER R.
       DUM = PI212/(AA+ALFA)
       DO 385 I = 1,IJ
       GIJ(I) = DIJ(I)*DUM*TWO
  385  CONTINUE
       XX=AA*AA*PCSQ/(AA+ALFA)
       IF(NROOTS.LE.3) CALL RT123
       IF(NROOTS.EQ.4) CALL ROOT4
       IF(NROOTS.EQ.5) CALL ROOT5
       IF(NROOTS.EQ.6) CALL ROOT6
       IF(NROOTS.GE.7) THEN
          IF (MASWRK) WRITE(IW,9008)
          CALL ABRT
       END IF
       MM = 0
       DO 425 K = 1,NROOTS
       UU = (ALFA+AA)*U(K)
       WW = W(K)
       WW = WW*(UU+ALFA)
       TT = ONE/(AA+UU+ALFA)
       T  =  SQRT(TT)
       X0 = (AAX+(UU+ALFA)*CX)*TT
       Y0 = (AAY+(UU+ALFA)*CY)*TT
       Z0 = (AAZ+(UU+ALFA)*CZ)*TT
       IN = -4+MM
       DO 405 I = 1,LIT
       IN = IN+4
       NI = I
       DO 405 J = 1,LJT
       JN = IN+J
       NJ = J
       CALL STVINT
       XIN(JN   ) = XINT
       YIN(JN   ) = YINT
       ZIN(JN   ) = ZINT*WW
      CALL POLXYZ
      XIN(JN+125) = XINT
      YIN(JN+125) = YINT
      ZIN(JN+125) = ZINT*WW
  405  CONTINUE
  425  MM = MM+16
       DO 403 I = 1,IJ
       NX    = IJX(I)
       NY    = IJY(I)
       NZ    = IJZ(I)
       MM    = 0
      DUMX = ZERO
      DUMY = ZERO
      DUMZ = ZERO
      DO 402 K = 1,NROOTS
      DUMX= DUMX+XIN(NX+MM+125)*YIN(NY+MM    )*ZIN(NZ+MM    )
      DUMY= DUMY+XIN(NX+MM    )*YIN(NY+MM+125)*ZIN(NZ+MM    )
      DUMZ= DUMZ+XIN(NX+MM    )*YIN(NY+MM    )*ZIN(NZ+MM+125)
  402 MM  = MM+16
      DUM = GIJ(I)*PREI*BETA
      CHGINT(ICC  )=CHGINT(ICC  )+DUM*DUMX
      CHGINT(ICC+1)=CHGINT(ICC+1)+DUM*DUMY
      CHGINT(ICC+2)=CHGINT(ICC+2)+DUM*DUMZ
      ICC=ICC+3
 403  CONTINUE
      END IF
C
 480  CONTINUE
 500  CONTINUE
 520  CONTINUE
C
C     ----- SET GRADIENT MATRIX
C
      MAX=MAXJ
      NN=0
      ICC=1
C
      DO 551 LTERM =1,4
      BETA = CLPR(LTERM,IFFAT)
      IF(BETA.EQ.ZERO) GOTO 551
      DO 550 I=MINI,MAXI
      LI=LOCI+I
      DO 550 J=MINJ,MAX
      LJ=LOCJ+J
      IF (LI-LJ) 920,940,940
  920 ID = LJ
      JD = LI
      GOTO 960
  940 ID = LI
      JD = LJ
  960 NN = (ID*(ID-1))/2+JD
      DUM = DM(NN)
      FFGRD(1,IFFAT)=FFGRD(1,IFFAT) + DUM*CHGINT(ICC  )
      FFGRD(2,IFFAT)=FFGRD(2,IFFAT) + DUM*CHGINT(ICC+1)
      FFGRD(3,IFFAT)=FFGRD(3,IFFAT) + DUM*CHGINT(ICC+2)
C
C     -- HUI LI: BEFORE DEC 10, 2016, VIRIAL IS CALCULATED WITH CX,
C        WHICH IS MM CORD INSTEAD OF QM-MM DISTANCE, BECAUSE
C        WE DO NOT KNOW WHICH ELECTRON (AND WHERE IT IS) GIVES
C        THE FROCE TO THE MM REP POINT (AND IT SHOULD NOT
C        BE DONE WITH QMCX BECAUSE THE FORCE COMES FROM ELECTRONS
C        INSTEAD OF QMCX SWITCHING FUNCTION)
C     VIR(1)        =VIR(1)         + DUM*CHGINT(ICC  )*CX
C     VIR(2)        =VIR(2)         + DUM*CHGINT(ICC+1)*CY
C     VIR(3)        =VIR(3)         + DUM*CHGINT(ICC+2)*CZ
C
C     -- HUI LI: AFTER DEC 10, 2016, VIRIAL IS CALCULATED WITH CXGRP,
C        WHICH IS QM-MM DISTANCE, AND THE
C        QM ELECTRON IS APPROXIMATELY AT XIJ, WHERE IT GIVES THE FORCE TO
C        THE MM REP POINT (AND IT SHOULD NOT BE DONE WITH
C        QMCX BECAUSE THE FORCE COMES FROM ELECTRONS
C        INSTEAD OF QMCX SWITCHING FUNCTION).
      IAT   = KATOM(II)
      JAT   = KATOM(JJ)
      XIJ   = 0.5D+00*(C(1,IAT) + C(1,JAT))
      YIJ   = 0.5D+00*(C(2,IAT) + C(2,JAT))
      ZIJ   = 0.5D+00*(C(3,IAT) + C(3,JAT))
      CXGRP = CORD(1,IFFAT) - XIJ
      CYGRP = CORD(2,IFFAT) - YIJ
      CZGRP = CORD(3,IFFAT) - ZIJ
      PBCX  = XBOX*ANINT(CXGRP*ONEXBOX)
      PBCY  = YBOX*ANINT(CYGRP*ONEYBOX)
      PBCZ  = ZBOX*ANINT(CZGRP*ONEZBOX)
      CXGRP = CXGRP - PBCX
      CYGRP = CYGRP - PBCY
      CZGRP = CZGRP - PBCZ
      VIR(1)        =VIR(1)         + DUM*CHGINT(ICC  )*CXGRP
      VIR(2)        =VIR(2)         + DUM*CHGINT(ICC+1)*CYGRP
      VIR(3)        =VIR(3)         + DUM*CHGINT(ICC+2)*CZGRP
      IYES = 0
      DO KFIX=1,NFIXMM
         IF(IFFAT.EQ.IFIXMM(KFIX)) IYES = 1
      ENDDO
      IF(IYES.EQ.1.AND.NFIXQM.GT.0) THEN
         VIR(1)     =VIR(1)         - DUM*CHGINT(ICC  )*CXGRP
         VIR(2)     =VIR(2)         - DUM*CHGINT(ICC+1)*CYGRP
         VIR(3)     =VIR(3)         - DUM*CHGINT(ICC+2)*CZGRP
      END IF
      ICC=ICC+3
  550 CONTINUE
  551 CONTINUE
  580 CONTINUE
  600 CONTINUE
  100 CONTINUE
C
      IF(GOPARR .AND. NXT) CALL DDI_DLBRESET
      RETURN
 9008 FORMAT(/' NUMBER OF POLYNOMIAL ROOTS NEEDED (NROOTS) IS GREATER',
     *        ' THAN 6 IN REPINT.  CALL A PROGRAMMER/QUANTUM CHEMIST.')
      END
C*MODULE QUANPOC  *DECK QMMMREPFQM
      SUBROUTINE QMMMREPFQM(DM,CORD,CLPR,ZLPR,NLPR,NONLSTQ)
      use mx_limits, only: mxsh,mxgtot,mxatm
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      LOGICAL ISKIP
      LOGICAL NORM,GOPARR,DSKWRK,MASWRK,NXT
C
C
      DIMENSION CORD(3,*),CLPR(4,*),ZLPR(4,*),NLPR(*),
     *          NONLSTQ(*)
C
      COMMON /DSTV  / XINT,YINT,ZINT,T,X0,Y0,Z0,XI,YI,ZI,
     *                XJ,YJ,ZJ,NI,NJ,CX,CY,CZ
      COMMON /FFPBSW/ XBOX,YBOX,ZBOX,SWF,SWFDX,SWFDY,SWFDZ,
     *                SWRA,ONESWRA,SWRA2,ONESWRA2,
     *                SWRB,ONESWRB,SWRB2,ONESWRB2,ONESWRB4,
     *                SWFDUM3,SWFDUM4,SWFDUM5,
     *                SWRAQ,ONESWRAQ,SWRAQ2,ONESWRAQ2,
     *                SWRBQ,ONESWRBQ,SWRBQ2,ONESWRBQ2,ONESWRBQ4,
     *                SWFDUM3Q,SWFDUM4Q,SWFDUM5Q,
     *                QMSIZE,QMCX,QMCY,QMCZ,QMCXSV,QMCYSV,QMCZSV,
     *                CENTX,CENTY,CENTZ,BUFWID1,BUFWID2,RDAMP,
     *                EFIELDX,EFIELDY,EFIELDZ,QMCXSV2,QMCYSV2,QMCZSV2,
     *                EPS1RB,EPS1RB3,ONEXBOX,ONEYBOX,ONEZBOX,
     *                LQMCT,MXLIST1,NTODO,NTODOSV,NTODOQ,
     *                ISWITCH,ISHIFT,IPOLSHF,
     *                LFFLSTCELL,LFFCORDSV,
     *                LFFPOLSV,LFFCORDSV2,LFFNONLS2,LFFCORDSVQ,
     *                LFFMVFASTS2,LFFMVFASTS3,LFFMVFASTS4,
     *                LFFMVFASTL2,LFFMVFASTL3,LFFMVFASTL4,
     *                MXCHECK,MXLIST2,NTODO2,NTODO2SV
      COMMON /GRAD  / DE(3,MXATM)
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /IOFILE/ IR,IW,IP,IJK,IJKT,IDAF,NAV,IODA(950)
      COMMON /NSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),
     *                CF(MXGTOT),CG(MXGTOT),CH(MXGTOT),CI(MXGTOT),
     *                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),
     *                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
      COMMON /OUTPUT/ NPRINT,ITOL,ICUT,NORMF,NORMP,NOPK
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /RUNOPT/ RUNTYP,EXETYP,NEVALS,NGLEVL,NHLEVL
      COMMON /ROOT  / XX,U(13),W(13),NROOTS
C
      PARAMETER (ZERO=0.0D+00)
      PARAMETER (ONE=1.0D+00)
      PARAMETER (RLN10=2.30258D+00)
      PARAMETER (PI212=1.1283791670955D+00)
      PARAMETER (SQRT3=1.73205080756888D+00)
      PARAMETER (SQRT5=2.23606797749979D+00)
C
      DIMENSION DM(*)
      DIMENSION ISKIP(35),INDX(70)
      DIMENSION IX(35),IY(35),IZ(35),JX(20),JY(20),JZ(20)
      DIMENSION IJG(210),IJX(210),IJY(210),IJZ(210)
      DIMENSION XIN(200),YIN(200),ZIN(200),G(210),DIJ(210)
      DIMENSION GIJ(210)
C
      DATA CHECK/8HCHECK   /
      DATA INDX / 1, 2, 3, 4, 5, 6, 7, 8, 9,10,
     1           11,12,13,14,15,16,17,18,19,20,
     2           -0,-0,-0,-0,-0,-0,-0,-0,-0,-0,
     3           -0,-0,-0,-0,-0,
     4           -0, 1, 2, 3, 1, 2, 3, 4, 5, 6,
     5            4, 5, 6, 7, 8, 9,10,11,12,13,
     6            7, 8, 9,10,11,12,13,14,15,16,
     7           17,18,19,20,21/
      DATA JX / 0, 1, 0, 0, 2, 0, 0, 1, 1, 0,
     1          3, 0, 0, 2, 2, 1, 0, 1, 0, 1/
      DATA IX / 1, 5, 1, 1, 9, 1, 1, 5, 5, 1,
     1         13, 1, 1, 9, 9, 5, 1, 5, 1, 5,
     2         17, 1, 1,13,13, 5, 1, 5, 1, 9,
     3          9, 1, 9, 5, 5/
      DATA JY / 0, 0, 1, 0, 0, 2, 0, 1, 0, 1,
     1          0, 3, 0, 1, 0, 2, 2, 0, 1, 1/
      DATA IY / 1, 1, 5, 1, 1, 9, 1, 5, 1, 5,
     1          1,13, 1, 5, 1, 9, 9, 1, 5, 5,
     2          1,17, 1, 5, 1,13,13, 1, 5, 9,
     3          1, 9, 5, 9, 5/
      DATA JZ / 0, 0, 0, 1, 0, 0, 2, 0, 1, 1,
     1          0, 0, 3, 0, 1, 0, 1, 2, 2, 1/
      DATA IZ / 1, 1, 1, 5, 1, 1, 9, 1, 5, 5,
     1          1, 1,13, 1, 5, 1, 5, 9, 9, 5,
     2          1, 1,17, 1, 5, 1, 5,13,13, 1,
     3          9, 9, 5, 5, 9/
C
C     HUI LI, JAN 2011, LINCOLN
C
C
C     INITIALIZATION FOR PARALLEL
C
      NXT = IBTYP.EQ.1
      IPCOUNT = ME - 1
      NEXT = -1
      LCNT = -1
      TOL = RLN10*ITOL
      NORM = NORMF .NE. 1 .OR. NORMP .NE. 1
C
      IF (EXETYP .EQ. CHECK) GOTO 1120
C
C     GRADIENT CONTRIBUTION FROM REPULSIVE POTENTIALS FOR POWERS 0
C     AND -1 OF R.
C
C     ----- ISHELL
C
      DO 1060 II = 1,NSHELL
      IF(NXT .AND. GOPARR) THEN
         LCNT = LCNT + 1
         IF(LCNT.GT.NEXT) CALL DDI_DLBNEXT(NEXT)
         IF(NEXT.NE.LCNT) GOTO 1060
      END IF
      IAT = KATOM(II)
      XI = C(1,IAT)
      YI = C(2,IAT)
      ZI = C(3,IAT)
      I1 = KSTART(II)
      I2 = I1+KNG(II)-1
      LIT = KTYPE(II)+1
      MINI = KMIN(II)
      MAXI = KMAX(II)
      LOCI = KLOC(II)-MINI
      DO 100 I=1,35
 100  ISKIP(I)=.TRUE.
      DO 260 I=MINI,MAXI
      GOTO (110,140,260,260,180,260,260,260,260,260,
     1       220,260,260,260,260,260,260,260,260,260),I
 110  DO 120 K=2,4
 120  ISKIP(K)=.FALSE.
      GOTO 260
 140  DO 160 K=5,10
 160  ISKIP(K)=.FALSE.
      ISKIP(1)=.FALSE.
      GOTO 260
 180  DO 200 K=2,4
 200  ISKIP(K)=.FALSE.
      DO 210 K=11,20
 210  ISKIP(K)=.FALSE.
      GOTO 260
 220  DO 240 K=5,10
 240  ISKIP(K)=.FALSE.
      DO 250 K=21,35
 250  ISKIP(K)=.FALSE.
 260  CONTINUE
C
C     ----- JSHELL
C
      DO 1040 JJ = 1,NSHELL
        IF((.NOT.NXT) .AND. GOPARR) THEN
           IPCOUNT = IPCOUNT + 1
           IF(MOD(IPCOUNT,NPROC).NE.0) GOTO 1040
        END IF
C
      JAT = KATOM(JJ)
      XJ = C(1,JAT)
      YJ = C(2,JAT)
      ZJ = C(3,JAT)
      J1 = KSTART(JJ)
      J2 = J1+KNG(JJ)-1
      LJT = KTYPE(JJ)
      MINJ = KMIN(JJ)
      MAXJ = KMAX(JJ)
      LOCJ = KLOC(JJ)-MINJ
      NROOTS = (LIT+LJT-2)/2 + 1
      RR = (XI-XJ)**2+(YI-YJ)**2+(ZI-ZJ)**2
C
C     ----- PREPARE INDICES FOR PAIRS OF (I,J) FUNCTIONS
C
      N0 = 0
      IF (LIT .GE. 4) N0 = 35
      IJ = 0
      DO 340 I = 1,35
      IF (ISKIP(I)) GOTO 340
      IN = INDX(I+N0)
      NX = IX(I)
      NY = IY(I)
      NZ = IZ(I)
      DO 320 J = MINJ,MAXJ
      IJ = IJ+1
      IJX(IJ) = NX+JX(J)
      IJY(IJ) = NY+JY(J)
      IJZ(IJ) = NZ+JZ(J)
      IJG(IJ) = IN+21*(J-MINJ)
  320 CONTINUE
  340 CONTINUE
      DO 360 I = 1,IJ
      N = IJG(I)
  360 G(N) = ZERO
C
C     ----- I PRIMITIVE
C
      DO 840 IG = I1,I2
      AI = EX(IG)
      ARRI = AI*RR
      AXI = AI*XI
      AYI = AI*YI
      AZI = AI*ZI
      DUM = AI+AI
      CSI = CP(IG)
      CPI = CS(IG)*DUM
      IF (LIT .EQ. 4) CPI = CD(IG)
      CDI = CP(IG)*DUM
      IF (LIT .EQ. 5) CDI = CF(IG)
      CFI = CD(IG)*DUM
      CGI = CF(IG)*DUM
C
C     ----- J PRIMITIVE
C
      DO 820 JG = J1,J2
      AJ = EX(JG)
      AA = AI+AJ
      AA1 = ONE/AA
      DUM = AJ*ARRI*AA1
      IF (DUM .GT. TOL) GOTO 820
      FAC = EXP(-DUM)
      CSJ = CS(JG)
      CPJ = CP(JG)
      CDJ = CD(JG)
      CFJ = CF(JG)
      AX = (AXI+AJ*XJ)*AA1
      AY = (AYI+AJ*YJ)*AA1
      AZ = (AZI+AJ*ZJ)*AA1
C
C     ----- DENSITY FACTOR
C
      NN = 0
      DUM1 = ZERO
      DUM2 = DUM1
      DO 600 I=1,35
      IF(ISKIP(I)) GOTO 600
      GOTO (370,380,420,420,390,420,420,420,420,420,
     1       400,420,420,420,420,420,420,420,420,420,
     2       410,420,420,420,420,420,420,420,420,420,
     3       420,420,420,420,420),I
 370  DUM1=CSI*FAC
      GOTO 420
 380  DUM1=CPI*FAC
      GOTO 420
 390  DUM1=CDI*FAC
      GOTO 420
 400  DUM1=CFI*FAC
      GOTO 420
 410  DUM1 = CGI*FAC
 420  CONTINUE
      DO 580 J=MINJ,MAXJ
      GOTO (430,440,560,560,460,560,560,480,560,560,
     1       500,560,560,520,560,560,560,560,560,540),J
  430 DUM2=DUM1*CSJ
      GOTO 560
  440 DUM2=DUM1*CPJ
      GOTO 560
  460 DUM2=DUM1*CDJ
      GOTO 560
  480 IF(NORM) DUM2=DUM2*SQRT3
      GOTO 560
  500 DUM2 = DUM1 * CFJ
      GOTO 560
  520 DUM2 = DUM2 *SQRT5
      GOTO 560
  540 DUM2 = DUM2 * SQRT3
  560 NN=NN+1
  580 DIJ(NN)=DUM2
  600 CONTINUE
      AAX = AA*AX
      AAY = AA*AY
      AAZ = AA*AZ
C
      DO 781 IIQ = 1, NTODOQ
      IFFAT = NONLSTQ(IIQ)
      IF((CLPR(1,IFFAT)+CLPR(2,IFFAT)+
     *    CLPR(3,IFFAT)+CLPR(4,IFFAT)).LT.1.0D-06) GOTO 781
      X     = CORD(1,IFFAT) - AX
      Y     = CORD(2,IFFAT) - AY
      Z     = CORD(3,IFFAT) - AZ
      PBCX  = XBOX*ANINT(X*ONEXBOX)
      PBCY  = YBOX*ANINT(Y*ONEYBOX)
      PBCZ  = ZBOX*ANINT(Z*ONEZBOX)
      CX    = CORD(1,IFFAT) - PBCX
      CY    = CORD(2,IFFAT) - PBCY
      CZ    = CORD(3,IFFAT) - PBCZ
      PCSQ = (AX-CX)**2+(AY-CY)**2+(AZ-CZ)**2
      IF(PCSQ.GT.SWRB2) GOTO 781    ! CUTOFF, NO SWF/SHIFT
      DO 780 LTERM=1,4
      ALFA = ZLPR(LTERM,IFFAT)
      BETA = CLPR(LTERM,IFFAT)
      IF(BETA.EQ.ZERO) GOTO 780
      PREI = EXP(-AA*ALFA*PCSQ/(AA+ALFA))
      IF(NLPR(IFFAT).EQ.2)THEN
C      R TO THE ZERO POWER.
       TT = ONE/(AA+ALFA)
       T  =  SQRT(TT)
       X0 = (AAX+ALFA*CX)*TT
       Y0 = (AAY+ALFA*CY)*TT
       Z0 = (AAZ+ALFA*CZ)*TT
       IN = -4
       DO 606 I = 1,LIT
       IN = IN+4
       NI = I
       DO 606 J = 1,LJT
       JN = IN+J
       NJ = J
       CALL VINT
       XIN(JN   ) = XINT*T
       YIN(JN   ) = YINT*T
       ZIN(JN   ) = ZINT*T
  606  CONTINUE
       DO 666 I = 1,IJ
       N = IJG(I)
       NX    = IJX(I)
       NY    = IJY(I)
       NZ    = IJZ(I)
       G(N)=G(N)+DIJ(I)*PREI*BETA*XIN(NX)*YIN(NY)*ZIN(NZ)
 666   CONTINUE
      ELSE
C      ONE OVER R.
       DUM = PI212/(AA+ALFA)
       DO 685 I = 1,IJ
       GIJ(I) = DIJ(I)*DUM
  685  CONTINUE
       XX=AA*AA*PCSQ/(AA+ALFA)
       IF(NROOTS .LE. 3) CALL RT123
       IF(NROOTS .EQ. 4) CALL ROOT4
       IF(NROOTS.EQ.5) CALL ROOT5
       IF(NROOTS.EQ.6) CALL ROOT6
       IF(NROOTS.GE.7)THEN
          IF (MASWRK) WRITE(IW,9008)
          CALL ABRT
       END IF
       MM = 0
       DO 725 K = 1,NROOTS
       UU = (ALFA+AA)*U(K)
       WW = W(K)
       TT = ONE/(AA+UU+ALFA)
       T  =  SQRT(TT)
       X0 = (AAX+(UU+ALFA)*CX)*TT
       Y0 = (AAY+(UU+ALFA)*CY)*TT
       Z0 = (AAZ+(UU+ALFA)*CZ)*TT
       IN = -4+MM
       DO 705 I = 1,LIT
       IN = IN+4
       NI = I
       DO 705 J = 1,LJT
       JN = IN+J
       NJ = J
       CALL VINT
       XIN(JN   ) = XINT
       YIN(JN   ) = YINT
       ZIN(JN   ) = ZINT*WW
  705  CONTINUE
  725  MM = MM+16
       DO 765 I = 1,IJ
       N = IJG(I)
       NX    = IJX(I)
       NY    = IJY(I)
       NZ    = IJZ(I)
       DUM = ZERO
       MM    = 0
       DO 745 K = 1,NROOTS
       DUM= DUM+XIN(NX+MM)*YIN(NY+MM)*ZIN(NZ+MM)
  745  MM  = MM+16
       G(N)=G(N)+GIJ(I)*PREI*BETA*DUM
 765   CONTINUE
      END IF
 780  CONTINUE
 781  CONTINUE
C
C     ----- END OF PRIMITIVE LOOPS -----
C
  820 CONTINUE
  840 CONTINUE
C
C     ----- FORM INTEGRALS OVER DERIVATIVES -----
C
      NN = 0
      N = 1
      DO 900 J = MINJ,MAXJ
      IF (MINI .GT. 1) GOTO 860
      NN = NN+1
      XIN(NN) = G(N+ 1)
      YIN(NN) = G(N+ 2)
      ZIN(NN) = G(N+ 3)
      IF (MAXI .EQ. 1) GOTO 900
  860 IF (MINI .GT. 2) GOTO 880
      NN = NN+1
      XIN(NN) = (G(N+ 4)-G(N ))
      YIN(NN) = G(N+ 7)
      ZIN(NN) = G(N+ 8)
      NN = NN+1
      XIN(NN) = G(N+ 7)
      YIN(NN) = (G(N+ 5)-G(N ))
      ZIN(NN) = G(N+ 9)
      NN = NN+1
      XIN(NN) = G(N+ 8)
      YIN(NN) = G(N+ 9)
      ZIN(NN) = (G(N+ 6)-G(N ))
      IF (MAXI .EQ. 4) GOTO 900
  880 CONTINUE
      IF (MINI .GT. 5) GOTO 890
      NN = NN+1
      XIN(NN) = (G(N+ 3)-G(N )-G(N ))
      YIN(NN) = G(N+ 6)
      ZIN(NN) = G(N+ 7)
      NN = NN+1
      XIN(NN) = G(N+ 8)
      YIN(NN) = (G(N+ 4)-G(N+ 1)-G(N+ 1))
      ZIN(NN) = G(N+ 9)
      NN = NN+1
      XIN(NN) = G(N+10)
      YIN(NN) = G(N+11)
      ZIN(NN) = (G(N+ 5)-G(N+ 2)-G(N+ 2))
      NN = NN+1
      DUM = ONE
      IF (NORM) DUM = SQRT3
      XIN(NN) = DUM*(G(N+ 6)-G(N+ 1))
      YIN(NN) = DUM*(G(N+ 8)-G(N ))
      ZIN(NN) = DUM* G(N+12)
      NN = NN+1
      XIN(NN) = DUM*(G(N+ 7)-G(N+ 2))
      YIN(NN) = DUM* G(N+12)
      ZIN(NN) = DUM*(G(N+10)-G(N ))
      NN = NN+1
      XIN(NN) = DUM* G(N+12)
      YIN(NN) = DUM*(G(N+ 9)-G(N+ 2))
      ZIN(NN) = DUM*(G(N+11)-G(N+ 1))
      IF(MAXI.EQ.10) GOTO 900
 890  CONTINUE
      NN=NN+1
      XIN(NN)=(G(N+ 6)-G(N   )-G(N   )-G(N   ))
      YIN(NN)= G(N+ 9)
      ZIN(NN)= G(N+10)
      NN=NN+1
      XIN(NN)= G(N+11)
      YIN(NN)=(G(N+ 7)-G(N+ 1)-G(N+ 1)-G(N+ 1))
      ZIN(NN)= G(N+12)
      NN=NN+1
      XIN(NN)= G(N+13)
      YIN(NN)= G(N+14)
      ZIN(NN)=(G(N+ 8)-G(N+ 2)-G(N+ 2)-G(N+ 2))
      NN=NN+1
      DUM=ONE
      IF(NORM) DUM=SQRT5
      XIN(NN)=DUM* (G(N+ 9)-G(N+ 3)-G(N+ 3))
      YIN(NN)=DUM* (G(N+15)-G(N  ))
      ZIN(NN)=DUM* G(N+18)
      NN=NN+1
      XIN(NN)=DUM* (G(N+10)-G(N+ 4)-G(N+ 4))
      YIN(NN)=DUM* G(N+18)
      ZIN(NN)=DUM* (G(N+16)-G(N  ))
      NN=NN+1
      XIN(NN)=DUM* (G(N+15)-G(N+ 1))
      YIN(NN)=DUM* (G(N+11)-G(N+ 3)-G(N+ 3))
      ZIN(NN)=DUM*  G(N+19)
      NN=NN+1
      XIN(NN)=DUM* G(N+19)
      YIN(NN)=DUM* (G(N+12)-G(N+ 5)-G(N+ 5))
      ZIN(NN)=DUM* (G(N+17)-G(N+ 1))
      NN=NN+1
      XIN(NN)=DUM* (G(N+16)-G(N+ 2))
      YIN(NN)=DUM*  G(N+20)
      ZIN(NN)=DUM* (G(N+13)-G(N+ 4)-G(N+ 4))
      NN=NN+1
      XIN(NN)=DUM*  G(N+20)
      YIN(NN)=DUM* (G(N+17)-G(N+ 2))
      ZIN(NN)=DUM* (G(N+14)-G(N+ 5)-G(N+ 5))
      NN=NN+1
      IF(NORM) DUM=DUM*SQRT3
      XIN(NN)=DUM* (G(N+18)-G(N+ 5))
      YIN(NN)=DUM* (G(N+19)-G(N+ 4))
      ZIN(NN)=DUM* (G(N+20)-G(N+ 3))
  900 N = N+21
C
C     ----- CALCULATE CONTRIBUTION TO GRADIENT -----
C
      N = 0
      DO 980 J = MINJ,MAXJ
      JN = LOCJ+J
      DO 980 I = MINI,MAXI
      N = N+1
      IN = LOCI+I
      IF (IN-JN) 920,940,940
  920 ID = JN
      JD = IN
      GOTO 960
  940 ID = IN
      JD = JN
  960 NN = (ID*(ID-1))/2+JD
      DUM = DM(NN)
      DUM = DUM+DUM
      DE(1,IAT) = DE(1,IAT)+DUM*XIN(N)
      DE(2,IAT) = DE(2,IAT)+DUM*YIN(N)
      DE(3,IAT) = DE(3,IAT)+DUM*ZIN(N)
  980 CONTINUE
 1040 CONTINUE
 1060 CONTINUE
C
C
 1120 CONTINUE
      IF(GOPARR .AND. NXT) CALL DDI_DLBRESET
      RETURN
 9008 FORMAT(/' NUMBER OF POLYNOMIAL ROOTS NEEDED (NROOTS) IS GREATER',
     *        ' THAN 6 IN REPINT.  CALL A PROGRAMMER/QUANTUM CHEMIST.')
C
      END
C*MODULE QUANPOC  *DECK QMMMPOLFLDINT
      SUBROUTINE QMMMPOLFLDINT(IFFAT,CORD,PEX,PEY,PEZ,
     *                         NUM2)
      use mx_limits, only: mxsh,mxgtot,mxatm
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      LOGICAL IANDJ,NORM,DOUBLE
C
C
      DIMENSION CORD(3,*),PEX(NUM2),PEY(NUM2),PEZ(NUM2)
C
      DIMENSION XIN(128),YIN(128),ZIN(128),PLX(100),PLY(100),PLZ(100)
      DIMENSION DIJ(100),FIJ(100)
      DIMENSION IX(20),IY(20),IZ(20),JX(20),JY(20),JZ(20)
      DIMENSION IJX(100),IJY(100),IJZ(100)
C
      COMMON /CSSTV / CX,CY,CZ
      COMMON /FFPARM/ NFFAT,NBOND,NANGL,NDIHR,NDIHB,NCMAP,NWAGG,
     *                N1213J,N14J,NLKQMM,IDOCHG,IDOPOL,IDOLJ,IDOCMAP
      COMMON /FFPBSW/ XBOX,YBOX,ZBOX,SWF,SWFDX,SWFDY,SWFDZ,
     *                SWRA,ONESWRA,SWRA2,ONESWRA2,
     *                SWRB,ONESWRB,SWRB2,ONESWRB2,ONESWRB4,
     *                SWFDUM3,SWFDUM4,SWFDUM5,
     *                SWRAQ,ONESWRAQ,SWRAQ2,ONESWRAQ2,
     *                SWRBQ,ONESWRBQ,SWRBQ2,ONESWRBQ2,ONESWRBQ4,
     *                SWFDUM3Q,SWFDUM4Q,SWFDUM5Q,
     *                QMSIZE,QMCX,QMCY,QMCZ,QMCXSV,QMCYSV,QMCZSV,
     *                CENTX,CENTY,CENTZ,BUFWID1,BUFWID2,RDAMP,
     *                EFIELDX,EFIELDY,EFIELDZ,QMCXSV2,QMCYSV2,QMCZSV2,
     *                EPS1RB,EPS1RB3,ONEXBOX,ONEYBOX,ONEZBOX,
     *                LQMCT,MXLIST1,NTODO,NTODOSV,NTODOQ,
     *                ISWITCH,ISHIFT,IPOLSHF,
     *                LFFLSTCELL,LFFCORDSV,
     *                LFFPOLSV,LFFCORDSV2,LFFNONLS2,LFFCORDSVQ,
     *                LFFMVFASTS2,LFFMVFASTS3,LFFMVFASTS4,
     *                LFFMVFASTL2,LFFMVFASTL3,LFFMVFASTL4,
     *                MXCHECK,MXLIST2,NTODO2,NTODO2SV
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /NSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),
     *                CF(MXGTOT),CG(MXGTOT),CH(MXGTOT),CI(MXGTOT),
     *                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),
     *                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
      COMMON /OUTPUT/ NPRINT,ITOL,ICUT,NORMF,NORMP,NOPK
      COMMON /ROOT  / XX,U(13),W(13),NROOTS
      COMMON /STV   / XINT,YINT,ZINT,T,X0,Y0,Z0,XI,YI,ZI,XJ,YJ,ZJ,NI,NJ
C
      DATA JX / 0, 1, 0, 0, 2, 0, 0, 1, 1, 0,
     1          3, 0, 0, 2, 2, 1, 0, 1, 0, 1/
      DATA IX / 1, 5, 1, 1, 9, 1, 1, 5, 5, 1,
     1         13, 1, 1, 9, 9, 5, 1, 5, 1, 5/
      DATA JY / 0, 0, 1, 0, 0, 2, 0, 1, 0, 1,
     1          0, 3, 0, 1, 0, 2, 2, 0, 1, 1/
      DATA IY / 1, 1, 5, 1, 1, 9, 1, 5, 1, 5,
     1          1,13, 1, 5, 1, 9, 9, 1, 5, 5/
      DATA JZ / 0, 0, 0, 1, 0, 0, 2, 0, 1, 1,
     1          0, 0, 3, 0, 1, 0, 1, 2, 2, 1/
      DATA IZ / 1, 1, 1, 5, 1, 1, 9, 1, 5, 5,
     1          1, 1,13, 1, 5, 1, 5, 9, 9, 5/
C
      DATA ZERO,ONE/0.0D+00,1.0D+00/
      DATA PI212 /1.1283791670955D+00/
      DATA SQRT3 /1.73205080756888D+00/
      DATA SQRT5 /2.23606797749979D+00/
      DATA RLN10 /2.30258D+00/
C
C     HUI LI, MAR 2011, LINCOLN
C
      IF(IDOPOL.EQ.0) RETURN
C
      TOL=RLN10*ITOL
      NORM=NORMF.NE.1.OR.NORMP.NE.1
C
C     -- FIELD INTEGRALS AT MM POL DUE TO QM BASIS FUNC
C
      CALL VCLR(PEX,1,NUM2)
      CALL VCLR(PEY,1,NUM2)
      CALL VCLR(PEZ,1,NUM2)
C
C     ----- ISHELL
C
      DO 9000 II=1,NSHELL
      I=KATOM(II)
      XI=C(1,I)
      YI=C(2,I)
      ZI=C(3,I)
      I1=KSTART(II)
      I2=I1+KNG(II)-1
      LIT=KTYPE(II)
      MINI=KMIN(II)
      MAXI=KMAX(II)
      LOCI=KLOC(II)-MINI
C
C     ----- JSHELL
C
      DO 8000 JJ=1,II
C
      J=KATOM(JJ)
      XJ=C(1,J)
      YJ=C(2,J)
      ZJ=C(3,J)
      J1=KSTART(JJ)
      J2=J1+KNG(JJ)-1
      LJT=KTYPE(JJ)
      MINJ=KMIN(JJ)
      MAXJ=KMAX(JJ)
      LOCJ=KLOC(JJ)-MINJ
      NROOTS=(LIT+LJT+1-2)/2 + 1
      RR=(XI-XJ)**2+(YI-YJ)**2+(ZI-ZJ)**2
      IANDJ=II.EQ.JJ
C
C     ----- PREPARE INDICES FOR PAIRS OF (I,J) FUNCTIONS
C
      IJ=0
      MAX=MAXJ
      DO 50 I=MINI,MAXI
      NX=IX(I)
      NY=IY(I)
      NZ=IZ(I)
      IF(IANDJ) MAX=I
      DO 50 J=MINJ,MAX
      IJ=IJ+1
      IJX(IJ)=NX+JX(J)
      IJY(IJ)=NY+JY(J)
      IJZ(IJ)=NZ+JZ(J)
   50 CONTINUE
      DO 60 I=1,IJ
      PLX(I) = ZERO
      PLY(I) = ZERO
  60  PLZ(I) = ZERO
C
C     ----- I PRIMITIVE
C
      JGMAX=J2
      DO 7000 IG=I1,I2
      AI=EX(IG)
      ARRI=AI*RR
      AXI=AI*XI
      AYI=AI*YI
      AZI=AI*ZI
      CSI=CS(IG)
      CPI=CP(IG)
      CDI=CD(IG)
      CFI=CF(IG)
C
C     ----- J PRIMITIVE
C
      IF(IANDJ) JGMAX=IG
      DO 6000 JG=J1,JGMAX
      AJ=EX(JG)
      AA=AI+AJ
      AA1=ONE/AA
      DUM=AJ*ARRI*AA1
      IF(DUM.GT.TOL) GOTO 6000
      FAC=EXP(-DUM)
      CSJ=CS(JG)
      CPJ=CP(JG)
      CDJ=CD(JG)
      CFJ=CF(JG)
      AX=(AXI+AJ*XJ)*AA1
      AY=(AYI+AJ*YJ)*AA1
      AZ=(AZI+AJ*ZJ)*AA1
C
C     ----- DENSITY FACTOR
C
      DOUBLE=IANDJ.AND.IG.NE.JG
      MAX=MAXJ
      NN=0
      DO 310 I=MINI,MAXI
      GOTO ( 70, 80,180,180, 90,180,180,100,180,180,
     *       110,180,180,120,180,180,180,180,180,130),I
   70 DUM1=CSI*FAC
      GOTO 180
   80 DUM1=CPI*FAC
      GOTO 180
   90 DUM1=CDI*FAC
      GOTO 180
  100 IF(NORM) DUM1=DUM1*SQRT3
      GOTO 180
  110 DUM1=CFI*FAC
      GOTO 180
  120 IF(NORM) DUM1=DUM1*SQRT5
      GOTO 180
  130 IF(NORM) DUM1=DUM1*SQRT3
  180 IF(IANDJ) MAX=I
      DO 310 J=MINJ,MAX
      GOTO (190,200,300,300,210,300,300,220,300,300,
     *       230,300,300,240,300,300,300,300,300,250),J
  190 DUM2=DUM1*CSJ
      IF(.NOT.DOUBLE) GOTO 300
      IF(I.GT.1) GOTO 195
      DUM2=DUM2+DUM2
      GOTO 300
  195 DUM2=DUM2+CSI*CPJ*FAC
      GOTO 300
  200 DUM2=DUM1*CPJ
      IF(DOUBLE) DUM2=DUM2+DUM2
      GOTO 300
  210 DUM2=DUM1*CDJ
      IF(DOUBLE) DUM2=DUM2+DUM2
      GOTO 300
  220 IF(NORM) DUM2=DUM2*SQRT3
      GOTO 300
  230 DUM2=DUM1*CFJ
      IF(DOUBLE) DUM2=DUM2+DUM2
      GOTO 300
  240 IF(NORM) DUM2=DUM2*SQRT5
      GOTO 300
  250 IF(NORM) DUM2=DUM2*SQRT3
  300 NN=NN+1
  310 DIJ(NN)=DUM2
C
      AAX=AA*AX
      AAY=AA*AY
      AAZ=AA*AZ
      DUM=PI212*AA1
      DUM=DUM+DUM
      DO 800 I=1,IJ
  800 FIJ(I)=DIJ(I)*DUM
C
      IAT   = KATOM(II)
      JAT   = KATOM(JJ)
      XIJ   = 0.5D+00*(C(1,IAT) + C(1,JAT))
      YIJ   = 0.5D+00*(C(2,IAT) + C(2,JAT))
      ZIJ   = 0.5D+00*(C(3,IAT) + C(3,JAT))
      X     = CORD(1,IFFAT) - XIJ
      Y     = CORD(2,IFFAT) - YIJ
      Z     = CORD(3,IFFAT) - ZIJ
      PBCX  = XBOX*ANINT(X*ONEXBOX)
      PBCY  = YBOX*ANINT(Y*ONEYBOX)
      PBCZ  = ZBOX*ANINT(Z*ONEZBOX)
      X     = X - PBCX
      Y     = Y - PBCY
      Z     = Z - PBCZ
      R2    = X*X+Y*Y+Z*Z
      IF(R2.GT.SWRB2) GOTO 6000
      IF(R2.LT.1.0D-10) THEN    !  R COULD BE ZERO
      SWF   = ONE
      SWFDX = ZERO
      SWFDY = ZERO
      SWFDZ = ZERO
      ELSE
      R     = SQRT(R2)
      ONER  = ONE/R
      IF(IPOLSHF.EQ.1) CALL SHIFT(R2,R,ONER,X,Y,Z)
      IF(IPOLSHF.EQ.0) CALL SWFUNC(R2,X,Y,Z)
      END IF
      CX    = CORD(1,IFFAT) - PBCX
      CY    = CORD(2,IFFAT) - PBCY
      CZ    = CORD(3,IFFAT) - PBCZ
C
      XX=AA*((AX-CX)**2+(AY-CY)**2+(AZ-CZ)**2)
      IF(NROOTS.LE.3) CALL RT123
      IF(NROOTS.EQ.4) CALL ROOT4
      MM=0
      DO 830 K=1,NROOTS
      UU=AA*U(K)
      WW=W(K)
      WW=WW*UU
      TT=ONE/(AA+UU)
      T=SQRT(TT)
      X0=(AAX+UU*CX)*TT
      Y0=(AAY+UU*CY)*TT
      Z0=(AAZ+UU*CZ)*TT
      IN=-4+MM
      DO 820 I=1,LIT
      IN=IN+4
      NI=I
      DO 820 J=1,LJT
      JN=IN+J
      NJ=J
      CALL STVINT
      XIN(JN)=XINT
      YIN(JN)=YINT
      ZIN(JN)=ZINT*WW
      CALL POLXYZ
      XIN(JN+64)=XINT
      YIN(JN+64)=YINT
      ZIN(JN+64)=ZINT*WW
  820 CONTINUE
  830 MM=MM+16
      DO 850 I=1,IJ
      NX=IJX(I)
      NY=IJY(I)
      NZ=IJZ(I)
      DUMX=ZERO
      DUMY=ZERO
      DUMZ=ZERO
      MM=0
      DO 840 K=1,NROOTS
      DUMX=DUMX+XIN(NX+MM+64)*YIN(NY+MM   )*ZIN(NZ+MM   )
      DUMY=DUMY+XIN(NX+MM   )*YIN(NY+MM+64)*ZIN(NZ+MM   )
      DUMZ=DUMZ+XIN(NX+MM   )*YIN(NY+MM   )*ZIN(NZ+MM+64)
  840 MM=MM+16
      DUM=FIJ(I)
      PLX(I) = PLX(I) + DUM*DUMX*SWF
      PLY(I) = PLY(I) + DUM*DUMY*SWF
      PLZ(I) = PLZ(I) + DUM*DUMZ*SWF
  850 CONTINUE
C
 6000 CONTINUE
 7000 CONTINUE
C
C
      MAX=MAXJ
      NN=0
      DO 7500 I=MINI,MAXI
      LI=LOCI+I
      IN = (LI*(LI-1))/2
      IF(IANDJ) MAX=I
      DO 7500 J=MINJ,MAX
      LJ=LOCJ+J
      JN=LJ+IN
      NN=NN+1
      PEX(JN)=PLX(NN)
      PEY(JN)=PLY(NN)
      PEZ(JN)=PLZ(NN)
 7500 CONTINUE
 8000 CONTINUE
 9000 CONTINUE
C
      RETURN
      END
C*MODULE QUANPOC  *DECK QMMMPOLFLDNUCMM
!>
!> @brief    QM nuclear and MM charge electric field at MM pol sites
!>
!> @author   Hui Li
!>           - Jan 2011
!>
!> @details  includes only QM nucear charges and MM charges
!>
      SUBROUTINE QMMMPOLFLDNUCMM(CORD,CHARG,POL,FIELD1,
     *                           NONLS1,NONLSTQ,L1213J)
      use mx_limits, only: mxatm
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      LOGICAL GOPARR,DSKWRK,MASWRK
C
      PARAMETER (ONE=1.0D+00)
      PARAMETER (ZERO=0.0D+00)
C
      DIMENSION CORD(3,*),CHARG(*),POL(*),FIELD1(3,*),
     *          NONLS1(2,*),NONLSTQ(*),L1213J(2,*)
C
      COMMON /FFNODE/ L1BOND,L2BOND,L1ANGL,L2ANGL,L1DIHR,L2DIHR,
     *                L1DIHB,L2DIHB,L1CMAP,L2CMAP,L1WAGG,L2WAGG,
     *                L11213,L21213,L1N14J,L2N14J,
     *                L11213A,L21213A,L1N14A,L2N14A,
     *                L11213B,L21213B,L1N14B,L2N14B,
     *                L1BONDPMA,L2BONDPMA,L1ANGLPMA,L2ANGLPMA,
     *                L1DIHRPMA,L2DIHRPMA,L1DIHBPMA,L2DIHBPMA,
     *                L1WAGGPMA,L2WAGGPMA,L1CMAPPMA,L2CMAPPMA,
     *                L11213PMA,L21213PMA,L1N14PMA,L2N14PMA,
     *                L1BONDPMB,L2BONDPMB,L1ANGLPMB,L2ANGLPMB,
     *                L1DIHRPMB,L2DIHRPMB,L1DIHBPMB,L2DIHBPMB,
     *                L1WAGGPMB,L2WAGGPMB,L1CMAPPMB,L2CMAPPMB,
     *                L11213PMB,L21213PMB,L1N14PMB,L2N14PMB,
     *                L1FFAT,L2FFAT
      COMMON /FFPARM/ NFFAT,NBOND,NANGL,NDIHR,NDIHB,NCMAP,NWAGG,
     *                N1213J,N14J,NLKQMM,IDOCHG,IDOPOL,IDOLJ,IDOCMAP
      COMMON /FFPBSW/ XBOX,YBOX,ZBOX,SWF,SWFDX,SWFDY,SWFDZ,
     *                SWRA,ONESWRA,SWRA2,ONESWRA2,
     *                SWRB,ONESWRB,SWRB2,ONESWRB2,ONESWRB4,
     *                SWFDUM3,SWFDUM4,SWFDUM5,
     *                SWRAQ,ONESWRAQ,SWRAQ2,ONESWRAQ2,
     *                SWRBQ,ONESWRBQ,SWRBQ2,ONESWRBQ2,ONESWRBQ4,
     *                SWFDUM3Q,SWFDUM4Q,SWFDUM5Q,
     *                QMSIZE,QMCX,QMCY,QMCZ,QMCXSV,QMCYSV,QMCZSV,
     *                CENTX,CENTY,CENTZ,BUFWID1,BUFWID2,RDAMP,
     *                EFIELDX,EFIELDY,EFIELDZ,QMCXSV2,QMCYSV2,QMCZSV2,
     *                EPS1RB,EPS1RB3,ONEXBOX,ONEYBOX,ONEZBOX,
     *                LQMCT,MXLIST1,NTODO,NTODOSV,NTODOQ,
     *                ISWITCH,ISHIFT,IPOLSHF,
     *                LFFLSTCELL,LFFCORDSV,
     *                LFFPOLSV,LFFCORDSV2,LFFNONLS2,LFFCORDSVQ,
     *                LFFMVFASTS2,LFFMVFASTS3,LFFMVFASTS4,
     *                LFFMVFASTL2,LFFMVFASTL3,LFFMVFASTL4,
     *                MXCHECK,MXLIST2,NTODO2,NTODO2SV
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
C
C     HUI LI, JAN 2011, LINCOLN
C
      IF(IDOPOL.EQ.0) RETURN
C
C     -- FIELD1: FIELD AT POL DUE TO MM CHARGE --
C
      CALL VCLR(FIELD1,1,3*NFFAT)
C
      DO LL=1,2
      IF(LL.EQ.1) THEN
         NN1  = 1
         NN2  = NTODO
         SIGN = 1.0D+00
      END IF
      IF(LL.EQ.2) THEN
         NN1  = L11213
         NN2  = L21213
         SIGN = -1.0D+00
      END IF
      DO 100 III=NN1, NN2
         IF(LL.EQ.1) THEN
            IFFAT = NONLS1(1,III)
            JFFAT = NONLS1(2,III)
         END IF
         IF(LL.EQ.2) THEN
            IFFAT = L1213J(1,III)
            JFFAT = L1213J(2,III)
         END IF
         IF(IFFAT.EQ.0.OR.JFFAT.EQ.0) GOTO 100
         IF(CHARG(IFFAT).EQ.ZERO.AND.  POL(IFFAT).EQ.ZERO) GOTO 100
         IF(CHARG(JFFAT).EQ.ZERO.AND.  POL(JFFAT).EQ.ZERO) GOTO 100
         IF(  POL(IFFAT).EQ.ZERO.AND.  POL(JFFAT).EQ.ZERO) GOTO 100
         IF(CHARG(IFFAT).EQ.ZERO.AND.CHARG(JFFAT).EQ.ZERO) GOTO 100
         QI    = CHARG(IFFAT)
         QJ    = CHARG(JFFAT)
C
         X     = CORD(1,IFFAT) - CORD(1,JFFAT)
         Y     = CORD(2,IFFAT) - CORD(2,JFFAT)
         Z     = CORD(3,IFFAT) - CORD(3,JFFAT)
         PBCX  = XBOX * ANINT(X*ONEXBOX)
         PBCY  = YBOX * ANINT(Y*ONEYBOX)
         PBCZ  = ZBOX * ANINT(Z*ONEZBOX)
         X     = X - PBCX
         Y     = Y - PBCY
         Z     = Z - PBCZ
         R2    = X*X+Y*Y+Z*Z
         IF(R2.GT.SWRB2) GOTO 100
         IF(R2.LT.0.01D+00) GOTO 100
         R     = SQRT(R2)
         ONER  = ONE/R
         IF(IPOLSHF.EQ.1) CALL SHIFT(R2,R,ONER,X,Y,Z)
         IF(IPOLSHF.EQ.0) CALL SWFUNC(R2,X,Y,Z)
C
         ONER2 = ONER*ONER
         ONER3 = ONER2*ONER
         SWF   = SWF*SIGN
         DUMI  = QJ*ONER3*SWF
         DUMJ  = QI*ONER3*SWF
         FIELD1(1,IFFAT)=FIELD1(1,IFFAT)+DUMI*X
         FIELD1(2,IFFAT)=FIELD1(2,IFFAT)+DUMI*Y
         FIELD1(3,IFFAT)=FIELD1(3,IFFAT)+DUMI*Z
         FIELD1(1,JFFAT)=FIELD1(1,JFFAT)-DUMJ*X
         FIELD1(2,JFFAT)=FIELD1(2,JFFAT)-DUMJ*Y
         FIELD1(3,JFFAT)=FIELD1(3,JFFAT)-DUMJ*Z
 100  CONTINUE
      ENDDO
C
C     -- FIELD1: FIELD AT POL DUE TO QM NUC --
C
      IF(NAT.LE.0) GOTO 330
      IPCOUNT = ME - 1
      DO 300 IIQ = 1, NTODOQ
         IFFAT = NONLSTQ(IIQ)
         IF(GOPARR) THEN
            IPCOUNT = IPCOUNT + 1
            IF(MOD(IPCOUNT,NPROC).NE.0) GOTO 300
         END IF
         IF(POL(IFFAT).EQ.ZERO) GOTO 300
         DO 310 JAT=1,NAT
            QJ    = ZAN(JAT)
C
            X     = CORD(1,IFFAT) - C(1,JAT)
            Y     = CORD(2,IFFAT) - C(2,JAT)
            Z     = CORD(3,IFFAT) - C(3,JAT)
            PBCX  = XBOX*ANINT(X*ONEXBOX)
            PBCY  = YBOX*ANINT(Y*ONEYBOX)
            PBCZ  = ZBOX*ANINT(Z*ONEZBOX)
            X     = X - PBCX
            Y     = Y - PBCY
            Z     = Z - PBCZ
            R2    = X*X+Y*Y+Z*Z
            IF(R2.GT.SWRB2) GOTO 310
            IF(R2.LT.1.0D-10) GOTO 310
            R     = SQRT(R2)
            ONER  = ONE/R
            ONER2 = ONER*ONER
            ONER3 = ONER2*ONER
            IF(IPOLSHF.EQ.1) CALL SHIFT(R2,R,ONER,X,Y,Z)
            IF(IPOLSHF.EQ.0) CALL SWFUNC(R2,X,Y,Z)
C
            DUMI  = QJ*ONER3*SWF
            FIELD1(1,IFFAT)=FIELD1(1,IFFAT)+DUMI*X
            FIELD1(2,IFFAT)=FIELD1(2,IFFAT)+DUMI*Y
            FIELD1(3,IFFAT)=FIELD1(3,IFFAT)+DUMI*Z
 310     CONTINUE
 300  CONTINUE
 330  CONTINUE
C
      IF(GOPARR) THEN
         CALL DDI_GSUMF(2405,FIELD1,3*NFFAT)
      END IF
C
      RETURN
      END
