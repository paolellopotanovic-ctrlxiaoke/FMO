C*MODULE QUANPOA_1  *DECK FFDATA
!>
!> @brief    read $FFDATA
!>
!> @author   Hui Li group
!>           - Nov 2015
!>
!> @details  read in all $FFDATA sections
!>
      SUBROUTINE FFDATA(ATMNAM,CORD,ZANF,
     *                  ZMAS,ONEMAS,
     *                  CHARG,POL,
     *                  SIG,EPS,
     *                  SIG2,EPS2,
     *                  BOND0,FCBOND,
     *                  ANGL0,FCANGL,FCWAGG,
     *                  DIHB0,FCDIHB,
     *                  VROT,NNN,GAMA,IPAIR,
     *                  KLIST,LLIST,
     *                  MLIST,NLIST,
     *                  VEL,QMVEL,
     *                  CLPR,ZLPR,NLPR,
     *                  MXFFAT,MXBOND,MXANGL,MXDIHR,MXCMAP,
     *                  MXWAGG,MXDIHB,NFOLD,
     *                  ITYPWAT,MAPLST,
     *                  INPQMV,INPMMV,IFFDATA,IFFPDB,
     *                  IDELETE,NORDER,ISCOOP,
     *                  KBLST,FCSTBD,FCDIHR,
     *                  FCLJTP,NTYPE)
      use mx_limits, only: mxatm
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      LOGICAL GOPARR,DSKWRK,MASWRK
C
      PARAMETER (TOANGS=0.52917724924D+00)
      PARAMETER (TOBOHR=1.0D+00/TOANGS)
      PARAMETER (PT5=0.5D+00)
      PARAMETER (TWO=2.0D+00)
      PARAMETER (HUGE=1.0D+30*TOBOHR)
C
      CHARACTER*10 WORD,ATMNAM
C
      DIMENSION ATMNAM(MXFFAT),CORD(3,MXFFAT),
     *          ZANF(MXFFAT),ZMAS(MXFFAT),ONEMAS(MXFFAT),
     *          CHARG(MXFFAT),POL(MXFFAT),SIG(MXFFAT),EPS(MXFFAT),
     *          SIG2(MXFFAT),EPS2(MXFFAT),
     *          BOND0(MXBOND),FCBOND(MXBOND),
     *          ANGL0(MXANGL),FCANGL(MXANGL),
     *          FCWAGG(MXWAGG),
     *          DIHB0(MXDIHB),FCDIHB(MXDIHB),
     *          VROT(MXDIHR),NNN(MXDIHR),GAMA(MXDIHR),IPAIR(2,MXBOND),
     *          KLIST(3,MXANGL),
     *          LLIST(4,MXDIHR),MLIST(4,MXWAGG),
     *          NLIST(4,MXDIHB),
     *          VEL(3,MXFFAT),QMVEL(3,MXATM),
     *          CLPR(4,MXFFAT),ZLPR(4,MXFFAT),NLPR(MXFFAT),
     *          MAPLST(6,*),NORDER(MXFFAT),
     *          KBLST(2,*),FCSTBD(2,*),FCDIHR(3,*),
     *          FCLJTP(2,MXMMTP,*),NTYPE(*)
C
      COMMON /FFMPT2/ MXMMTP,LFFKBLST,LFFFCSTBD,LFFFCDIHR,
     *                LFFFCLJTP,LFFNTYPE,
     *                LFF2KBLST,LFF2FCSTBD,LFF2FCDIHR,
     *                LFF2FCLJTP,LFF2NTYPE
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
      COMMON /FFSPH / SPHRAD,SPHSIG,SPHEPS,IADDWAT
      COMMON /FMCOM / X(1)
      COMMON /IOFILE/ IR,IW,IP,IJK,IJKT,IDAF,NAV,IODA(950)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
C
C     NANDUN THELLAMUREGE, DEJUN SI, HUI LI, MAR 2011, LINCOLN
C     HONGBO ZHU, HUI LI, OCT 30, 2012, LINCOLN
C
C     -- READ FROM INPUT FILE --
C
      CALL OPNCRD(IR,-IW)
      IEOF = 0
      IF(IFFDATA.EQ.0.AND.IFFPDB.GT.0) GOTO 290
C
  110 CONTINUE
      CALL RDCARD('$FFDATA ',IEOF)
      IF(IEOF.EQ.1) THEN
         IF(IFFDATA.EQ.1)
     *      WRITE(IW,*)'ERROR: END OF FILE READING $FFDATA'
         IF(IFFDATA.EQ.2)
     *      WRITE(IW,*)'ERROR: END OF FILE READING $FFDATB'
         CALL ABRT
      END IF
      WORD ='          '
      KSIZE = -10
      CALL GSTRNG(WORD,KSIZE)
C
      IF(WORD.NE.'COORDINATE')GOTO 120
      CALL RDXYZ(ATMNAM,CORD,ZANF,ZMAS,ONEMAS,MXFFAT,NFOLD)
      GOTO 110
C
  120 CONTINUE
      IF(WORD.NE.'MMVELOCITY')GOTO 125
      INPMMV = 1
      CALL RDMMVEL(VEL)
      GOTO 110
C
  125 CONTINUE
      IF(WORD.NE.'QMVELOCITY')GOTO 126
      INPQMV = 1
      CALL RDQMVEL(QMVEL)
      GOTO 110
C
  126 CONTINUE
      IF(WORD.NE.'PARAMETERS')GOTO 130
      CALL RDPARA(ATMNAM,ZMAS,ONEMAS,CHARG,POL,SIG,EPS,SIG2,EPS2)
      GOTO 110
C
  130 CONTINUE
      IF(WORD.NE.'QMMMREP   ')GOTO 135
      CALL RDQMMM(CLPR,ZLPR,NLPR)
      GOTO 110
C
  135 CONTINUE
      IF(WORD.NE.'BOND      ')GOTO 140
      CALL RDBOND(BOND0,FCBOND,IPAIR,MXBOND)
      GOTO 110
C
  140 CONTINUE
      IF(WORD.NE.'ANGLE     ')GOTO 150
      CALL RDANGL(ANGL0,FCANGL,KLIST,MXANGL)
      GOTO 110
C
  150 CONTINUE
      IF(WORD.NE.'DIHROT    ')GOTO 160
      CALL RDDIHR(LLIST,GAMA,NNN,VROT,MXDIHR)
      GOTO 110
C
  160 CONTINUE
      IF(WORD.NE.'DIHBND    ')GOTO 170
      CALL RDDIHB(NLIST,DIHB0,FCDIHB,MXDIHB)
      GOTO 110
C
  170 CONTINUE
      IF(WORD.NE.'WAGGING   ')GOTO 180
      CALL RDWAGG(MLIST,FCWAGG,MXWAGG)
      GOTO 110
C
  180 CONTINUE
      IF(WORD.NE.'CMAP      ')GOTO 190
      CALL RDCMAP(MAPLST,MXCMAP)
      GOTO 110
C
  190 CONTINUE
      IF(WORD.NE.'STRBEND   ')GOTO 200
      CALL RDSTRBEND(FCSTBD,KBLST)
      GOTO 110
C
  200 CONTINUE
      IF(WORD.NE.'DIHR3V    ')GOTO 210
      CALL RDDIHR3V(LLIST,FCDIHR)
      GOTO 110
C
  210 CONTINUE
      IF(WORD.NE.'MMFFLJ    ')GOTO 220
      CALL RDMMFFLJ(FCLJTP)
      GOTO 110
C
  220 CONTINUE
      IF(WORD.NE.'MMTYPE    ')GOTO 230
      CALL RDMMTYPE(NTYPE)
      GOTO 110
C
  230 CONTINUE
      IF(WORD.EQ.'$END      ') THEN
         GOTO 300
      ELSE
         IF(IFFDATA.EQ.1)
     *   WRITE(IW,*)'ERROR: UNRECOGNIZED KEYWORD IN $FFDATA'
         IF(IFFDATA.EQ.2)
     *   WRITE(IW,*)'ERROR: UNRECOGNIZED KEYWORD IN $FFDATB'
         CALL ABRT
      END IF
C
  290 CONTINUE
C
C     --- READ $FFPDB AND LOAD PROTEIN PARAMETERS ---
C
      CALL VALFM(LOADFM)
      LSSBDID  = LOADFM   + 1
      LSSBOND  = LSSBDID  + MXFFAT/14
      LPDBNAM  = LSSBOND  + MXFFAT/14
      LRESNAM  = LPDBNAM  + MXFFAT
      LNATAAA  = LRESNAM  + MXFFAT
      LNNNCCC  = LNATAAA  + MXFFAT
      LAMBNAM  = LNNNCCC  + MXFFAT/7
      LLLIST1  = LAMBNAM  + MXFFAT
      LTEXTA   = LLLIST1  + MXDIHR*4
      LTEXTB   = LTEXTA   + 20000*100*2
      LTEXTC   = LTEXTB   + 20000*100*2
      LTEXTD   = LTEXTC   + 20000*100*2
      LAST     = LTEXTD   + 20000*100*2
      NEED     = LAST     - LOADFM -1
      CALL GETFM(NEED)
      CALL FFPDB(ATMNAM,CORD,ZANF,
     *           ZMAS,ONEMAS,
     *           CHARG,POL,SIG,EPS,
     *           SIG2,EPS2,
     *           BOND0,FCBOND,
     *           ANGL0,FCANGL,FCWAGG,
     *           DIHB0,FCDIHB,
     *           VROT,NNN,GAMA,IPAIR,
     *           KLIST,LLIST,MLIST,NLIST,
     *           CLPR,ZLPR,NLPR,
     *           X(LSSBDID),X(LSSBOND),
     *           X(LPDBNAM),X(LRESNAM),X(LNATAAA),
     *           X(LNNNCCC),
     *           MXFFAT,MXBOND,MXANGL,MXDIHR,
     *           MXWAGG,MXDIHB,
     *           ITYPWAT,MAPLST,X(LAMBNAM),X(LLLIST1),
     *           X(LTEXTA),X(LTEXTB),X(LTEXTC),X(LTEXTD))
      CALL RETFM(NEED)
      INPMMV  = 0
      CALL TIMIT(1)
 300  CONTINUE
C
C     -- APPLY IDELETE OR ISCOOP FUNCTION --
C
      IF(IDELETE.GT.0 .OR. ISCOOP.GT.0) THEN
         CALL VICLR(NORDER,1,NFFAT)
C
         IF(IDELETE.GT.0) THEN
            DO JFFAT = 1, IDELETE
               DO IFFAT = IDELETE+1, NFFAT
                  DX = CORD(1,IFFAT) - CORD(1,JFFAT)
                  DY = CORD(2,IFFAT) - CORD(2,JFFAT)
                  DZ = CORD(3,IFFAT) - CORD(3,JFFAT)
                  R2 = DX*DX+DY*DY+DZ*DZ
                  IF(R2.LT.3.571D+00) ZMAS(IFFAT) = -1.0D+00
               ENDDO
            ENDDO
         END IF
C
         IF(ISCOOP.GT.0) THEN
         IF(CENTX.EQ.HUGE.OR.CENTY.EQ.HUGE.OR.CENTZ.EQ.HUGE) THEN
            XMAX = -1.0D+30
            YMAX = -1.0D+30
            ZMAX = -1.0D+30
            XMIN =  1.0D+30
            YMIN =  1.0D+30
            ZMIN =  1.0D+30
            DO IFFAT = 1,NFFAT
               XMAX = MAX(XMAX,CORD(1,IFFAT))
               YMAX = MAX(YMAX,CORD(2,IFFAT))
               ZMAX = MAX(ZMAX,CORD(3,IFFAT))
               XMIN = MIN(XMIN,CORD(1,IFFAT))
               YMIN = MIN(YMIN,CORD(2,IFFAT))
               ZMIN = MIN(ZMIN,CORD(3,IFFAT))
            ENDDO
            CENTX= (XMAX+XMIN)*PT5
            CENTY= (YMAX+YMIN)*PT5
            CENTZ= (ZMAX+ZMIN)*PT5
         END IF
         END IF
C
         IF(ISCOOP.EQ.1) THEN
            DO IFFAT = 1, NFFAT
               DX = TWO*ABS(CORD(1,IFFAT) - CENTX)
               DY = TWO*ABS(CORD(2,IFFAT) - CENTY)
               DZ = TWO*ABS(CORD(3,IFFAT) - CENTZ)
               IF(DX.GT.XBOX.OR.DY.GT.YBOX.OR.DZ.GT.ZBOX)
     *         ZMAS(IFFAT) = -1.0D+00
            ENDDO
         END IF
C
         IF(ISCOOP.EQ.2) THEN
            TEST = SPHRAD**2
            DO IFFAT = 1, NFFAT
               DX = CORD(1,IFFAT) - CENTX
               DY = CORD(2,IFFAT) - CENTY
               DZ = CORD(3,IFFAT) - CENTZ
               R2 = DX*DX+DY*DY+DZ*DZ
               IF(R2.GT.TEST)
     *         ZMAS(IFFAT) = -1.0D+00
            ENDDO
         END IF
C
 310     CONTINUE
         ICOUNT = 0
         DO IBOND=1, NBOND
            IFFAT = IPAIR(1,IBOND)
            JFFAT = IPAIR(2,IBOND)
            IF(ZMAS(IFFAT).EQ.-1.0D+00.AND.ZMAS(JFFAT).NE.-1.0D+00) THEN
               ZMAS(JFFAT) = -1.0D+00
               ICOUNT = ICOUNT + 1
            END IF
            IF(ZMAS(JFFAT).EQ.-1.0D+00.AND.ZMAS(IFFAT).NE.-1.0D+00) THEN
               ZMAS(IFFAT) = -1.0D+00
               ICOUNT = ICOUNT + 1
            END IF
         ENDDO
         IF(ICOUNT.NE.0) GOTO 310
C
         KFFAT = 0
         DO IFFAT=1, NFFAT
            IF(ZMAS(IFFAT).GT.-1.0D+00) THEN
               KFFAT           = KFFAT + 1
               NORDER(IFFAT)   = KFFAT
               ATMNAM(KFFAT)   = ATMNAM(IFFAT)
               ZANF  (KFFAT)   = ZANF  (IFFAT)
               ZMAS  (KFFAT)   = ZMAS  (IFFAT)
               ONEMAS(KFFAT)   = ONEMAS(IFFAT)
               CHARG (KFFAT)   = CHARG (IFFAT)
               POL   (KFFAT)   = POL   (IFFAT)
               SIG   (KFFAT)   = SIG   (IFFAT)
               EPS   (KFFAT)   = EPS   (IFFAT)
               SIG2  (KFFAT)   = SIG2  (IFFAT)
               EPS2  (KFFAT)   = EPS2  (IFFAT)
               CORD  (1,KFFAT) = CORD  (1,IFFAT)
               CORD  (2,KFFAT) = CORD  (2,IFFAT)
               CORD  (3,KFFAT) = CORD  (3,IFFAT)
               VEL   (1,KFFAT) = VEL   (1,IFFAT)
               VEL   (2,KFFAT) = VEL   (2,IFFAT)
               VEL   (3,KFFAT) = VEL   (3,IFFAT)
               CLPR  (1,KFFAT) = CLPR  (1,IFFAT)
               CLPR  (2,KFFAT) = CLPR  (2,IFFAT)
               CLPR  (3,KFFAT) = CLPR  (3,IFFAT)
               CLPR  (4,KFFAT) = CLPR  (4,IFFAT)
               ZLPR  (1,KFFAT) = ZLPR  (1,IFFAT)
               ZLPR  (2,KFFAT) = ZLPR  (2,IFFAT)
               ZLPR  (3,KFFAT) = ZLPR  (3,IFFAT)
               ZLPR  (4,KFFAT) = ZLPR  (4,IFFAT)
               NLPR  (KFFAT)   = NLPR  (IFFAT)
               NTYPE (KFFAT)   = NTYPE (IFFAT)
            END IF
         ENDDO
         NFFAT = KFFAT
C
         KBOND = 0
         DO IBOND=1, NBOND
            IFFAT = IPAIR(1,IBOND)
            JFFAT = IPAIR(2,IBOND)
            IF(NORDER(IFFAT).GT.0) THEN
               KBOND          = KBOND + 1
               IPAIR(1,KBOND) = NORDER(IFFAT)
               IPAIR(2,KBOND) = NORDER(JFFAT)
               FCBOND (KBOND) = FCBOND(IBOND)
               BOND0  (KBOND) = BOND0 (IBOND)
            END IF
         ENDDO
         NBOND = KBOND
C
         KANGL = 0
         DO IANGL=1, NANGL
            IFFAT = KLIST(1,IANGL)
            JFFAT = KLIST(2,IANGL)
            KFFAT = KLIST(3,IANGL)
            IF(NORDER(IFFAT).GT.0) THEN
               KANGL          = KANGL + 1
               KLIST(1,KANGL) = NORDER(IFFAT)
               KLIST(2,KANGL) = NORDER(JFFAT)
               KLIST(3,KANGL) = NORDER(KFFAT)
               FCANGL (KANGL) = FCANGL(IANGL)
               ANGL0  (KANGL) = ANGL0 (IANGL)
               FCSTBD(1,KANGL)= FCSTBD(1,IANGL)
               FCSTBD(2,KANGL)= FCSTBD(2,IANGL)
            END IF
         ENDDO
         NANGL = KANGL
C        - MUST REDO KBLST -
         DO IANGL=1,NANGL
            K1 = KLIST(1,IANGL)
            K2 = KLIST(2,IANGL)
            K3 = KLIST(3,IANGL)
            DO IBOND=1,NBOND
               KK1 = IPAIR(1,IBOND)
               KK2 = IPAIR(2,IBOND)
               IF((KK1.EQ.K1.AND.KK2.EQ.K2).OR.
     *            (KK1.EQ.K2.AND.KK2.EQ.K1)) THEN
                  KBLST(1,IANGL) = IBOND
               END IF
               IF((KK1.EQ.K2.AND.KK2.EQ.K3).OR.
     *            (KK1.EQ.K3.AND.KK2.EQ.K2)) THEN
                  KBLST(2,IANGL) = IBOND
               END IF
            ENDDO
         ENDDO
C
         KDIHR = 0
         DO IDIHR=1, NDIHR
            IFFAT = LLIST(1,IDIHR)
            JFFAT = LLIST(2,IDIHR)
            KFFAT = LLIST(3,IDIHR)
            LFFAT = LLIST(4,IDIHR)
            IF(NORDER(IFFAT).GT.0) THEN
               KDIHR          = KDIHR + 1
               LLIST(1,KDIHR) = NORDER(IFFAT)
               LLIST(2,KDIHR) = NORDER(JFFAT)
               LLIST(3,KDIHR) = NORDER(KFFAT)
               LLIST(4,KDIHR) = NORDER(LFFAT)
               NNN    (KDIHR) = NNN (IDIHR)
               GAMA   (KDIHR) = GAMA(IDIHR)
               VROT   (KDIHR) = VROT(IDIHR)
               FCDIHR(1,KDIHR)= FCDIHR(1,IDIHR)
               FCDIHR(2,KDIHR)= FCDIHR(2,IDIHR)
               FCDIHR(3,KDIHR)= FCDIHR(3,IDIHR)
            END IF
         ENDDO
         NDIHR = KDIHR
C
         KDIHB = 0
         DO IDIHB=1, NDIHB
            IFFAT = NLIST(1,IDIHB)
            JFFAT = NLIST(2,IDIHB)
            KFFAT = NLIST(3,IDIHB)
            LFFAT = NLIST(4,IDIHB)
            IF(NORDER(IFFAT).GT.0) THEN
               KDIHB          = KDIHB + 1
               NLIST(1,KDIHB) = NORDER(IFFAT)
               NLIST(2,KDIHB) = NORDER(JFFAT)
               NLIST(3,KDIHB) = NORDER(KFFAT)
               NLIST(4,KDIHB) = NORDER(LFFAT)
               DIHB0 (KDIHB)  = DIHB0 (IDIHB)
               FCDIHB(KDIHB)  = FCDIHB(IDIHB)
            END IF
         ENDDO
         NDIHB = KDIHB
C
         KWAGG = 0
         DO IWAGG=1, NWAGG
            IFFAT = MLIST(1,IWAGG)
            JFFAT = MLIST(2,IWAGG)
            KFFAT = MLIST(3,IWAGG)
            LFFAT = MLIST(4,IWAGG)
            IF(NORDER(IFFAT).GT.0) THEN
               KWAGG          = KWAGG + 1
               MLIST(1,KWAGG) = NORDER(IFFAT)
               MLIST(2,KWAGG) = NORDER(JFFAT)
               MLIST(3,KWAGG) = NORDER(KFFAT)
               MLIST(4,KWAGG) = NORDER(LFFAT)
               FCWAGG (KWAGG) = FCWAGG(IWAGG)
            END IF
         ENDDO
         NWAGG = KWAGG
C
         KCMAP = 0
         DO ICMAP=1, NCMAP
            II1 = MAPLST(1,ICMAP)
            II2 = MAPLST(2,ICMAP)
            II3 = MAPLST(3,ICMAP)
            II4 = MAPLST(4,ICMAP)
            II5 = MAPLST(5,ICMAP)
            II6 = MAPLST(6,ICMAP)
            IF(NORDER(II1).GT.0) THEN
               KCMAP = KCMAP + 1
               MAPLST(1,KCMAP) = II1
               MAPLST(2,KCMAP) = II2
               MAPLST(3,KCMAP) = II3
               MAPLST(4,KCMAP) = II4
               MAPLST(5,KCMAP) = II5
               MAPLST(6,KCMAP) = II6
            END IF
         ENDDO
         NCMAP = KCMAP
      END IF
C
C     -- FINE TUNE THE INPUT BOND, ANGLE LISTS --
C        FORCE THEM TO BE FROM LOWER TO HIGHER NUMBER.
C        WATER ADDED LATER WILL BE DONE SEPARATELY
C        IN THE SAME WAY.
C
      DO IBOND = 1, NBOND
         IF(IPAIR(1,IBOND).LE.0 .OR.
     *      IPAIR(2,IBOND).LE.0     ) THEN
            IF(MASWRK)WRITE(IW,'(2I7,/1X,A/)')
     *         IPAIR(1,IBOND),IPAIR(2,IBOND),
     *         'ERROR: BOND ATM1 OR ATM2 MUST BE > 0'
            CALL ABRT
         END IF
         IF(IPAIR(2,IBOND).EQ.IPAIR(1,IBOND)) THEN
            IF(MASWRK)WRITE(IW,'(2I7,/1X,A/)')
     *         IPAIR(1,IBOND),IPAIR(2,IBOND),
     *         'ERROR: BOND ATM1, ATM2 MUST BE DIFFERENT'
            CALL ABRT
         END IF
         IF(IPAIR(2,IBOND).LT.IPAIR(1,IBOND)) THEN
            IZ5G9          = IPAIR(1,IBOND)
            IPAIR(1,IBOND) = IPAIR(2,IBOND)
            IPAIR(2,IBOND) = IZ5G9
          END IF
      ENDDO
C
      DO IANGL = 1, NANGL
         IF(KLIST(1,IANGL).LE.0.OR.
     *      KLIST(2,IANGL).LE.0.OR.
     *      KLIST(3,IANGL).LE.0    ) THEN
            IF(MASWRK)WRITE(IW,'(3I7,/1X,A/)')
     *         KLIST(1,IANGL),KLIST(2,IANGL),KLIST(3,IANGL),
     *         'ERROR: ANGLE ATM1, ATM2, ATM3 MUST BE > 0'
               CALL ABRT
         END IF
         IF(KLIST(1,IANGL).EQ.KLIST(2,IANGL).OR.
     *      KLIST(1,IANGL).EQ.KLIST(3,IANGL).OR.
     *      KLIST(2,IANGL).EQ.KLIST(3,IANGL)) THEN
            IF(MASWRK)WRITE(IW,'(3I7,/1X,A/)')
     *         KLIST(1,IANGL),KLIST(2,IANGL),KLIST(3,IANGL),
     *         'ERROR: ANGLE ATM1, ATM2, ATM3 MUST BE DIFFERENT'
               CALL ABRT
         END IF
         IF(KLIST(1,IANGL).GT.KLIST(3,IANGL)) THEN
            IZ5G9          = KLIST(1,IANGL)
            KLIST(1,IANGL) = KLIST(3,IANGL)
            KLIST(3,IANGL) = IZ5G9
            IZ5G9          = KBLST(1,IANGL)
            KBLST(1,IANGL) = KBLST(2,IANGL)
            KBLST(2,IANGL) = IZ5G9
            AZ5G9          = FCSTBD(1,IANGL)
            FCSTBD(1,IANGL)= FCSTBD(2,IANGL)
            FCSTBD(2,IANGL)= AZ5G9
         END IF
      ENDDO
C
      DO IDIHR = 1, NDIHR
         IF(LLIST(1,IDIHR).LE.0.OR.
     *      LLIST(2,IDIHR).LE.0.OR.
     *      LLIST(3,IDIHR).LE.0.OR.
     *      LLIST(4,IDIHR).LE.0    ) THEN
            IF(MASWRK)WRITE(IW,'(4I7,/1X,A/)')
     *         LLIST(1,IDIHR),LLIST(2,IDIHR),
     *         LLIST(3,IDIHR),LLIST(4,IDIHR),
     *        'ERROR: DIHR ATM1, ATM2, ATM3, ATM4 MUST BE > 0'
               CALL ABRT
         END IF
         IF(LLIST(1,IDIHR).EQ.LLIST(2,IDIHR).OR.
     *      LLIST(1,IDIHR).EQ.LLIST(3,IDIHR).OR.
     *      LLIST(1,IDIHR).EQ.LLIST(4,IDIHR).OR.
     *      LLIST(2,IDIHR).EQ.LLIST(3,IDIHR).OR.
     *      LLIST(2,IDIHR).EQ.LLIST(4,IDIHR).OR.
     *      LLIST(3,IDIHR).EQ.LLIST(4,IDIHR)) THEN
            IF(MASWRK)WRITE(IW,'(4I7,/1X,A/)')
     *         LLIST(1,IDIHR),LLIST(2,IDIHR),
     *         LLIST(3,IDIHR),LLIST(4,IDIHR),
     *        'ERROR: DIHR ATM1, ATM2, ATM3, ATM4 MUST BE DIFFERENT'
               CALL ABRT
         END IF
         IF(LLIST(1,IDIHR).GT.LLIST(4,IDIHR)) THEN
            IZ5G9          = LLIST(1,IDIHR)
            LLIST(1,IDIHR) = LLIST(4,IDIHR)
            LLIST(4,IDIHR) = IZ5G9
            IZ5G9          = LLIST(2,IDIHR)
            LLIST(2,IDIHR) = LLIST(3,IDIHR)
            LLIST(3,IDIHR) = IZ5G9
         END IF
      ENDDO
C
      DO IDIHB = 1, NDIHB
         IF(NLIST(1,IDIHB).LE.0.OR.
     *      NLIST(2,IDIHB).LE.0.OR.
     *      NLIST(3,IDIHB).LE.0.OR.
     *      NLIST(4,IDIHB).LE.0    ) THEN
            IF(MASWRK)WRITE(IW,'(4I7,/1X,A/)')
     *         NLIST(1,IDIHB),NLIST(2,IDIHB),
     *         NLIST(3,IDIHB),NLIST(4,IDIHB),
     *        'ERROR: DIHB ATM1, ATM2, ATM3, ATM4 MUST BE > 0'
               CALL ABRT
         END IF
         IF(NLIST(1,IDIHB).EQ.NLIST(2,IDIHB).OR.
     *      NLIST(1,IDIHB).EQ.NLIST(3,IDIHB).OR.
     *      NLIST(1,IDIHB).EQ.NLIST(4,IDIHB).OR.
     *      NLIST(2,IDIHB).EQ.NLIST(3,IDIHB).OR.
     *      NLIST(2,IDIHB).EQ.NLIST(4,IDIHB).OR.
     *      NLIST(3,IDIHB).EQ.NLIST(4,IDIHB)) THEN
            IF(MASWRK)WRITE(IW,'(4I7,/1X,A/)')
     *         NLIST(1,IDIHB),NLIST(2,IDIHB),
     *         NLIST(3,IDIHB),NLIST(4,IDIHB),
     *        'ERROR: DIHB ATM1, ATM2, ATM3, ATM4 MUST BE DIFFERENT'
               CALL ABRT
         END IF
         IF(NLIST(1,IDIHB).GT.NLIST(4,IDIHB)) THEN
            IZ5G9          = NLIST(1,IDIHB)
            NLIST(1,IDIHB) = NLIST(4,IDIHB)
            NLIST(4,IDIHB) = IZ5G9
            IZ5G9          = NLIST(2,IDIHB)
            NLIST(2,IDIHB) = NLIST(3,IDIHB)
            NLIST(3,IDIHB) = IZ5G9
         END IF
      ENDDO
C
      RETURN
      END