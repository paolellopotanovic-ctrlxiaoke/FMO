C*MODULE QUANPOB  *DECK JADDIONS
!>
!> @brief    ADD NA+ OR K+ IONS TO DNA/RNA PO4 SITES
!>
!> @author   Hui Li
!>           - Dec 2011
!>
!> @details  ADD NA+ OR K+ IONS TO DNA/RNA PO4 SITES
!>
      SUBROUTINE JADDIONS(ATMNAM,CORD,ZANF,ZMAS,ONEMAS,CHARG,POL,
     *                    SIG,EPS,SIG2,EPS2,CLPR,ZLPR,NLPR,
     *                    JADDNA1,JADDK1)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      CHARACTER*8  RNAME
      CHARACTER*10 ATMNAM
C
      PARAMETER (PT5=0.5D+00)
      PARAMETER (FOUR=4.0D+00)
C
      DIMENSION ATMNAM(*),CORD(3,*),ZANF(*),ZMAS(*),
     *          ONEMAS(*),CHARG(*),POL(*),
     *          SIG(*),EPS(*),SIG2(*),EPS2(*),CLPR(4,*),ZLPR(4,*),
     *          NLPR(*)
C
      COMMON /FFPARM/ NFFAT,NBOND,NANGL,NDIHR,NDIHB,NCMAP,NWAGG,
     *                N1213J,N14J,NLKQMM,IDOCHG,IDOPOL,IDOLJ,IDOCMAP
      COMMON /FFTYPE/ WT14LJ,WT14CH,C3BOND,C4BOND,C3ANGL,
     *                NFFTYP,NFFFILE,LJQMMM,LJQM,INTCHG,
     *                LJSIGMA,JTOPFILE(90),JPARFILE(90),
     *                JTOPAMIA(90),JTOPNTER(90),JTOPCTER(90),
     *                JTOPNUCA(90),JPARFIL2(90),JPARFIL3(90)
C
C     ADD NA+ OR K+ IONS TO DNA/RNA PO4 SITES.
C
C     HUI LI, DEC 6, 2011, LINCOLN
C
      IF(JADDNA1.EQ.0 .AND. JADDK1.EQ.0) RETURN
C
      MADDNA1= 0
      MADDK1 = 0
      DO IFFAT = 1, NFFAT
         IF(ATMNAM(IFFAT  ).EQ.'P         '.AND.
     *      ATMNAM(IFFAT+1).EQ.'O         '.AND.
     *      ATMNAM(IFFAT+2).EQ.'O         ') THEN
            DIS1 = (CORD(1,IFFAT)-CORD(1,IFFAT+1))**2 +
     *             (CORD(2,IFFAT)-CORD(2,IFFAT+1))**2 +
     *             (CORD(3,IFFAT)-CORD(3,IFFAT+1))**2
            DIS2 = (CORD(1,IFFAT)-CORD(1,IFFAT+2))**2 +
     *             (CORD(2,IFFAT)-CORD(2,IFFAT+2))**2 +
     *             (CORD(3,IFFAT)-CORD(3,IFFAT+2))**2
            IF(DIS1.LT.11.6D+00 .AND. DIS2.LT.11.6D+00) THEN
               XXX = CORD(1,IFFAT)
     *             + ((CORD(1,IFFAT+1)+CORD(1,IFFAT+2))*PT5
     *                -CORD(1,IFFAT))*FOUR
               YYY = CORD(2,IFFAT)
     *             + ((CORD(2,IFFAT+1)+CORD(2,IFFAT+2))*PT5
     *                -CORD(2,IFFAT))*FOUR
               ZZZ = CORD(3,IFFAT)
     *             + ((CORD(3,IFFAT+1)+CORD(3,IFFAT+2))*PT5
     *                -CORD(3,IFFAT))*FOUR
               IF(JADDNA1.EQ.1) THEN
                  MADDNA1         = MADDNA1 + 1
                  NFFAT           = NFFAT + 1
                  CORD(1,NFFAT)   = XXX
                  CORD(2,NFFAT)   = YYY
                  CORD(3,NFFAT)   = ZZZ
                  ZANF  (NFFAT)   = 11.0D+00
                  ATMNAM(NFFAT)   = 'NA'
                  IF(NFFTYP/10000.LE.9)WRITE(RNAME,'(A5,I1,A2)')
     *                           ' $NA1',NFFTYP/10000,'  '
                  IF(NFFTYP/10000.GT.9)WRITE(RNAME,'(A5,I2,A1)')
     *                           ' $NA1',NFFTYP/10000,' '
                  CALL RDIONS(RNAME,ZMAS,ONEMAS,
     *                        CHARG,POL,SIG,EPS,
     *                        SIG2,EPS2,
     *                        CLPR,ZLPR,NLPR,NFFAT)
               ELSE IF(JADDK1.EQ.1) THEN
                  MADDK1          = MADDK1 + 1
                  NFFAT           = NFFAT + 1
                  CORD(1,NFFAT)   = XXX
                  CORD(2,NFFAT)   = YYY
                  CORD(3,NFFAT)   = ZZZ
                  ZANF  (NFFAT)   = 19.0D+00
                  ATMNAM(NFFAT)   = 'K'
                  IF(NFFTYP/10000.LE.9)WRITE(RNAME,'(A4,I1,A3)')
     *                           ' $K1',NFFTYP/10000,'   '
                  IF(NFFTYP/10000.GT.9)WRITE(RNAME,'(A4,I2,A2)')
     *                           ' $K1',NFFTYP/10000,'  '
                  CALL RDIONS(RNAME,ZMAS,ONEMAS,
     *                        CHARG,POL,SIG,EPS,
     *                        SIG2,EPS2,
     *                        CLPR,ZLPR,NLPR,NFFAT)
               END IF
            END IF
         END IF
      ENDDO
      JADDNA1 = MADDNA1
      JADDK1  = MADDK1
C
      RETURN
      END
C*MODULE QUANPOB  *DECK ADDPBCWAT
!>
!> @brief    Add water molecules to a PBC box
!>
!> @author   Hui Li
!>           - Mar 2011
!>
!> @details  add ions first
!>           water parameters from external file
!>
      SUBROUTINE ADDPBCWAT(ATMNAM,CORD,ZANF,ZMAS,ONEMAS,CHARG,POL,
     *                     SIG,EPS,SIG2,EPS2,CLPR,ZLPR,NLPR,
     *                     IPAIR,FCBOND,BOND0,
     *                     KLIST,FCANGL,ANGL0,L1213J,
     *                     MXFFAT,MXBOND,MXANGL,
     *                     WATO1,WATH2,WATH3,NWATER,ITYPWAT,
     *                     IADDNA1,IADDK1,IADDCA2,IADDMG2,IADDCL1,
     *                     JADDNA1,JADDK1,LSTRAT,DSTRAT)
      use mx_limits, only: mxatm
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      LOGICAL GOPARR,MASWRK,DSKWRK
C
      CHARACTER*8  RNAME
      CHARACTER*10 ATMNAM
      CHARACTER*256 QPFILE
C
      PARAMETER (ZERO=0.0D+00)
      PARAMETER (PT5=0.5D+00)
      PARAMETER (FOUR=4.0D+00)
      PARAMETER (TOANGS=0.52917724924D+00)
      PARAMETER (TOBOHR=1.0D+00/TOANGS)
C
      DIMENSION ATMNAM(*),CORD(3,*),ZANF(*),ZMAS(*),
     *          ONEMAS(*),CHARG(*),POL(*),
     *          SIG(*),EPS(*),SIG2(*),EPS2(*),CLPR(4,*),ZLPR(4,*),
     *          NLPR(*),IPAIR(2,*),FCBOND(*),BOND0(*),
     *          KLIST(3,*),FCANGL(*),ANGL0(*),L1213J(2,*),
     *          WATO1(3,*),WATH2(3,*),WATH3(3,*),
     *          O1(3),H2(3),H3(3),LSTRAT(2,*),DSTRAT(*)
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
      COMMON /FFRATT/ RATOLC,RATOLV,SCALRAT,VIRRAT(3),IRATTLE,JRATTLE,
     *                NRATTLE,MXRATT,LFFOLDCORD,LFFLSTRAT,LFFDSTRAT,
     *                LFFVELSV,IRATQM
      COMMON /FFTYPE/ WT14LJ,WT14CH,C3BOND,C4BOND,C3ANGL,
     *                NFFTYP,NFFFILE,LJQMMM,LJQM,INTCHG,
     *                LJSIGMA,JTOPFILE(90),JPARFILE(90),
     *                JTOPAMIA(90),JTOPNTER(90),JTOPCTER(90),
     *                JTOPNUCA(90),JPARFIL2(90),JPARFIL3(90)
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /IOFILE/ IR,IW,IP,IJK,IJKT,IDAF,NAV,IODA(950)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
C
C     ADD WATER MOLECULES TO THE SYSTEM TO FILL A RECTANGULAR PBC BOX.
C
C     HUI LI, MAR 2011, LINCOLN
C
      NFFAT0 = NFFAT
C
      XMAXMOL = -1.0D+30
      YMAXMOL = -1.0D+30
      ZMAXMOL = -1.0D+30
      XMINMOL =  1.0D+30
      YMINMOL =  1.0D+30
      ZMINMOL =  1.0D+30
      DO IAT = 1,NAT
         XMAXMOL = MAX(XMAXMOL,C(1,IAT))
         YMAXMOL = MAX(YMAXMOL,C(2,IAT))
         ZMAXMOL = MAX(ZMAXMOL,C(3,IAT))
         XMINMOL = MIN(XMINMOL,C(1,IAT))
         YMINMOL = MIN(YMINMOL,C(2,IAT))
         ZMINMOL = MIN(ZMINMOL,C(3,IAT))
      ENDDO
      DO IFFAT = 1,NFFAT0
         XMAXMOL = MAX(XMAXMOL,CORD(1,IFFAT))
         YMAXMOL = MAX(YMAXMOL,CORD(2,IFFAT))
         ZMAXMOL = MAX(ZMAXMOL,CORD(3,IFFAT))
         XMINMOL = MIN(XMINMOL,CORD(1,IFFAT))
         YMINMOL = MIN(YMINMOL,CORD(2,IFFAT))
         ZMINMOL = MIN(ZMINMOL,CORD(3,IFFAT))
      ENDDO
      IF(XMAXMOL.EQ.-1.0D+30) XMAXMOL = ZERO
      IF(YMAXMOL.EQ.-1.0D+30) YMAXMOL = ZERO
      IF(ZMAXMOL.EQ.-1.0D+30) ZMAXMOL = ZERO
      IF(XMINMOL.EQ. 1.0D+30) XMINMOL = ZERO
      IF(YMINMOL.EQ. 1.0D+30) YMINMOL = ZERO
      IF(ZMINMOL.EQ. 1.0D+30) ZMINMOL = ZERO
      IF(MASWRK) THEN
         WRITE(IW,*)'THE DIMENSION OF THE SOLUTE MOLECULE IS (A):'
         WRITE(IW,'(6(1X,A,F8.2))')
     *     'XMIN=',XMINMOL*TOANGS,'XMAX=',XMAXMOL*TOANGS,
     *     'YMIN=',YMINMOL*TOANGS,'YMAX=',YMAXMOL*TOANGS,
     *     'ZMIN=',ZMINMOL*TOANGS,'ZMAX=',ZMAXMOL*TOANGS
      END IF
C     - THE AVERAGE DISTANCE BETWEEN AN ATOM
C       AND WATER OXYGEN ATOM SHOULD BE 2.3 ANGSTROM
C       THIS VALUE DOES NOT AFFECT THE NUMBER OF WATER
C       BECAUSE IT IS USED ONLY FOR PRE-SCREENING
      VDW     = 2.3D+00*TOBOHR
      XMAXMOL = XMAXMOL + VDW
      YMAXMOL = YMAXMOL + VDW
      ZMAXMOL = ZMAXMOL + VDW
      XMINMOL = XMINMOL - VDW
      YMINMOL = YMINMOL - VDW
      ZMINMOL = ZMINMOL - VDW
      IF((XBOX-(XMAXMOL-XMINMOL)).LT.ZERO)THEN
         IF(MASWRK) WRITE(IW,*)'ERROR: XBOX MUST BE LARGER THAN ',
     *      (XMAXMOL-XMINMOL)*TOANGS,' ANGSTROM.'
         IF(MASWRK) WRITE(IW,*) XMAXMOL*TOANGS,XMINMOL*TOANGS
         IF(MASWRK)WRITE(IW,*)' '
         CALL ABRT
      END IF
      IF((YBOX-(YMAXMOL-YMINMOL)).LT.ZERO)THEN
         IF(MASWRK) WRITE(IW,*)'ERROR: YBOX MUST BE LARGER THAN ',
     *      (YMAXMOL-YMINMOL)*TOANGS,' ANGSTROM.'
         IF(MASWRK)WRITE(IW,*)' '
         CALL ABRT
      END IF
      IF((ZBOX-(ZMAXMOL-ZMINMOL)).LT.ZERO)THEN
         IF(MASWRK) WRITE(IW,*)'ERROR: ZBOX MUST BE LARGER THAN ',
     *      (ZMAXMOL-ZMINMOL)*TOANGS,' ANGSTROM.'
         IF(MASWRK)WRITE(IW,*)' '
         CALL ABRT
      END IF
C
C     -- SET UP A WATER BOX CENTERED AT THE DEFINED PBC CENTER
      XMAX    = CENTX + PT5*XBOX
      YMAX    = CENTY + PT5*YBOX
      ZMAX    = CENTZ + PT5*ZBOX
      XMIN    = CENTX - PT5*XBOX
      YMIN    = CENTY - PT5*YBOX
      ZMIN    = CENTZ - PT5*ZBOX
C
C     - DETERMINE THE MAXIMUM NUMBER OF WATER
      VOLM    = (MAX(XBOX,YBOX,ZBOX))**3
C     - VOLUME OF ONE WATER =  29.998696 A**3 AT 298.15 K, 1 BAR
C                           = 202.441191 BOHR**3
      MAXWAT  = INT(VOLM/202.441191D+00) + 10
      NFOLD   = 0
      DO K = 0, 30, 3
         IF(MAXWAT.GT.4096*(2**K)) NFOLD = K+3
      ENDDO
      KKKWAT = 4096*(2**NFOLD)
      IF(KKKWAT.GT.8*MXFFAT) THEN
         IF(MASWRK) WRITE(IW,*)
     *      'ERROR: TOO MANY WATER IN ADDPBCWAT. INCREASE MXFFAT'
         IF(MASWRK)WRITE(IW,*)' '
         CALL ABRT
      END IF
C
C     -- READ IN 4096 WATER MOLECULES --
C        THEY ARE IN A CUBE WITH SIDE LENGTH = 49.7150 A.
C        CREATE MORE WATER IF NECESSARY
C
      IF(MASWRK) THEN
      CALL GENQPFILE('WATER4096.DAT ',QPFILE,LENQP)
      OPEN(12,FILE=QPFILE(1:LENQP),STATUS='OLD')
      DO IWAT=1,4096
         READ(12,*) WATO1(1,IWAT),WATO1(2,IWAT),WATO1(3,IWAT),
     *              WATH2(1,IWAT),WATH2(2,IWAT),WATH2(3,IWAT),
     *              WATH3(1,IWAT),WATH3(2,IWAT),WATH3(3,IWAT)
      ENDDO
      CLOSE(12)
      END IF
      IF(GOPARR) CALL DDI_BCAST(457,'F',WATO1,3*4096,MASTER)
      IF(GOPARR) CALL DDI_BCAST(458,'F',WATH2,3*4096,MASTER)
      IF(GOPARR) CALL DDI_BCAST(459,'F',WATH3,3*4096,MASTER)
      SIZE  = 49.7150D+00*TOBOHR*PT5
      DO IFOLD = 1, NFOLD
         IXYZ1 = MOD(IFOLD,3)
         IF(IXYZ1.EQ.1) IXYZ2 = 2
         IF(IXYZ1.EQ.1) IXYZ3 = 3
         IF(IXYZ1.EQ.2) IXYZ2 = 1
         IF(IXYZ1.EQ.2) IXYZ3 = 3
         IF(IXYZ1.EQ.0) IXYZ2 = 1
         IF(IXYZ1.EQ.0) IXYZ3 = 2
         IF(IXYZ1.EQ.0) IXYZ1 = 3
         IF(IXYZ1.EQ.1) SIZE  = SIZE*2.0D+00
         NDONE  = 2**(IFOLD-1)
         LENGTH = NDONE*4096
         DO III = 1, LENGTH
            WATO1(IXYZ1,LENGTH+III) = WATO1(IXYZ1,III) + SIZE
            WATO1(IXYZ2,LENGTH+III) = WATO1(IXYZ2,III)
            WATO1(IXYZ3,LENGTH+III) = WATO1(IXYZ3,III)
            WATH2(IXYZ1,LENGTH+III) = WATH2(IXYZ1,III) + SIZE
            WATH2(IXYZ2,LENGTH+III) = WATH2(IXYZ2,III)
            WATH2(IXYZ3,LENGTH+III) = WATH2(IXYZ3,III)
            WATH3(IXYZ1,LENGTH+III) = WATH3(IXYZ1,III) + SIZE
            WATH3(IXYZ2,LENGTH+III) = WATH3(IXYZ2,III)
            WATH3(IXYZ3,LENGTH+III) = WATH3(IXYZ3,III)
         ENDDO
      ENDDO
C
C     -- ADD IONS BEFORE WATER --
C        (1) ADD NA+ OR K+ TO DNA/RNA, IF REQUESTED.
C        (2) ADD OTHER IONS
C
      MADDNA1= 0
      MADDK1 = 0
      IF(JADDNA1.EQ.1 .OR. JADDK1.EQ.1) THEN
      DO IFFAT = 1, NFFAT
         IF(ATMNAM(IFFAT  ).EQ.'P         '.AND.
     *      ATMNAM(IFFAT+1).EQ.'O         '.AND.
     *      ATMNAM(IFFAT+2).EQ.'O         ') THEN
            DIS1 = (CORD(1,IFFAT)-CORD(1,IFFAT+1))**2 +
     *             (CORD(2,IFFAT)-CORD(2,IFFAT+1))**2 +
     *             (CORD(3,IFFAT)-CORD(3,IFFAT+1))**2
            DIS2 = (CORD(1,IFFAT)-CORD(1,IFFAT+2))**2 +
     *             (CORD(2,IFFAT)-CORD(2,IFFAT+2))**2 +
     *             (CORD(3,IFFAT)-CORD(3,IFFAT+2))**2
            IF(DIS1.LT.11.6D+00 .AND. DIS2.LT.11.6D+00) THEN
               XXX = CORD(1,IFFAT)
     *             + ((CORD(1,IFFAT+1)+CORD(1,IFFAT+2))*PT5
     *                -CORD(1,IFFAT))*FOUR
               YYY = CORD(2,IFFAT)
     *             + ((CORD(2,IFFAT+1)+CORD(2,IFFAT+2))*PT5
     *                -CORD(2,IFFAT))*FOUR
               ZZZ = CORD(3,IFFAT)
     *             + ((CORD(3,IFFAT+1)+CORD(3,IFFAT+2))*PT5
     *                -CORD(3,IFFAT))*FOUR
               IF(JADDNA1.EQ.1) THEN
                  MADDNA1         = MADDNA1 + 1
                  NFFAT           = NFFAT + 1
                  CORD(1,NFFAT)   = XXX
                  CORD(2,NFFAT)   = YYY
                  CORD(3,NFFAT)   = ZZZ
                  ZANF  (NFFAT)   = 11.0D+00
                  ATMNAM(NFFAT)   = 'NA'
                  IF(NFFTYP/10000.LE.9)WRITE(RNAME,'(A5,I1,A2)')
     *                           ' $NA1',NFFTYP/10000,'  '
                  IF(NFFTYP/10000.GT.9)WRITE(RNAME,'(A5,I2,A1)')
     *                           ' $NA1',NFFTYP/10000,' '
                  CALL RDIONS(RNAME,ZMAS,ONEMAS,
     *                        CHARG,POL,SIG,EPS,
     *                        SIG2,EPS2,
     *                        CLPR,ZLPR,NLPR,NFFAT)
               ELSE IF(JADDK1.EQ.1) THEN
                  MADDK1          = MADDK1 + 1
                  NFFAT           = NFFAT + 1
                  CORD(1,NFFAT)   = XXX
                  CORD(2,NFFAT)   = YYY
                  CORD(3,NFFAT)   = ZZZ
                  ZANF  (NFFAT)   = 19.0D+00
                  ATMNAM(NFFAT)   = 'K'
                  IF(NFFTYP/10000.LE.9)WRITE(RNAME,'(A4,I1,A3)')
     *                           ' $K1',NFFTYP/10000,'   '
                  IF(NFFTYP/10000.GT.9)WRITE(RNAME,'(A4,I2,A2)')
     *                           ' $K1',NFFTYP/10000,'  '
                  CALL RDIONS(RNAME,ZMAS,ONEMAS,
     *                        CHARG,POL,SIG,EPS,
     *                        SIG2,EPS2,
     *                        CLPR,ZLPR,NLPR,NFFAT)
               END IF
            END IF
         END IF
      ENDDO
      JADDNA1 = MADDNA1
      JADDK1  = MADDK1
      END IF
C
      NADDNA1= 0
      NADDK1 = 0
      NADDMG2= 0
      NADDCA2= 0
      NADDCL1= 0
C
 100  CONTINUE
      IF(IADDNA1.EQ.0) GOTO 101
      CALL FFRAND(XXX)
      CALL FFRAND(YYY)
      CALL FFRAND(ZZZ)
      XXX   = CENTX + XBOX*(XXX-PT5)
      YYY   = CENTY + YBOX*(YYY-PT5)
      ZZZ   = CENTZ + ZBOX*(ZZZ-PT5)
      DO IFFAT = 1, NFFAT
         CX    = XXX - CORD(1,IFFAT)
         CY    = YYY - CORD(2,IFFAT)
         CZ    = ZZZ - CORD(3,IFFAT)
         PBCX  = XBOX*ANINT(CX*ONEXBOX)
         PBCY  = YBOX*ANINT(CY*ONEYBOX)
         PBCZ  = ZBOX*ANINT(CZ*ONEZBOX)
         CX    = CX - PBCX
         CY    = CY - PBCY
         CZ    = CZ - PBCZ
         R2    = CX*CX+CY*CY+CZ*CZ
C        -- CANNOT BE SMALLER THAN 4.0 ANGSTROM
         IF(R2.LT.57.14D+00) GOTO 100
      ENDDO
      NADDNA1         = NADDNA1 + 1
      NFFAT           = NFFAT + 1
      CORD(1,NFFAT)   = XXX
      CORD(2,NFFAT)   = YYY
      CORD(3,NFFAT)   = ZZZ
      ZANF  (NFFAT)   = 11.0D+00
      ATMNAM(NFFAT)   = 'NA'
      IF(NFFTYP/10000.LE.9)WRITE(RNAME,'(A5,I1,A2)')
     *                     ' $NA1',NFFTYP/10000,'  '
      IF(NFFTYP/10000.GT.9)WRITE(RNAME,'(A5,I2,A1)')
     *                     ' $NA1',NFFTYP/10000,' '
      CALL RDIONS(RNAME,ZMAS,ONEMAS,
     *            CHARG,POL,SIG,EPS,
     *            SIG2,EPS2,
     *            CLPR,ZLPR,NLPR,NFFAT)
      IF(NADDNA1.LT.IADDNA1) GOTO 100
C
 101  CONTINUE
      IF(IADDCL1.EQ.0) GOTO 102
      CALL FFRAND(XXX)
      CALL FFRAND(YYY)
      CALL FFRAND(ZZZ)
      XXX   = CENTX + XBOX*(XXX-PT5)
      YYY   = CENTY + YBOX*(YYY-PT5)
      ZZZ   = CENTZ + ZBOX*(ZZZ-PT5)
      DO IFFAT = 1, NFFAT
         CX    = XXX - CORD(1,IFFAT)
         CY    = YYY - CORD(2,IFFAT)
         CZ    = ZZZ - CORD(3,IFFAT)
         PBCX  = XBOX*ANINT(CX*ONEXBOX)
         PBCY  = YBOX*ANINT(CY*ONEYBOX)
         PBCZ  = ZBOX*ANINT(CZ*ONEZBOX)
         CX    = CX - PBCX
         CY    = CY - PBCY
         CZ    = CZ - PBCZ
         R2    = CX*CX+CY*CY+CZ*CZ
C        -- CANNOT BE SMALLER THAN 4.0 ANGSTROM
         IF(R2.LT.57.14D+00) GOTO 101
      ENDDO
      NADDCL1         = NADDCL1 + 1
      NFFAT           = NFFAT + 1
      CORD(1,NFFAT)   = XXX
      CORD(2,NFFAT)   = YYY
      CORD(3,NFFAT)   = ZZZ
      ZANF  (NFFAT)   = 17.0D+00
      ATMNAM(NFFAT)   = 'CL'
      IF(NFFTYP/10000.LE.9)WRITE(RNAME,'(A5,I1,A2)')
     *                     ' $CL1',NFFTYP/10000,'  '
      IF(NFFTYP/10000.GT.9)WRITE(RNAME,'(A5,I2,A1)')
     *                     ' $CL1',NFFTYP/10000,' '
      CALL RDIONS(RNAME,ZMAS,ONEMAS,
     *            CHARG,POL,SIG,EPS,
     *            SIG2,EPS2,
     *            CLPR,ZLPR,NLPR,NFFAT)
      IF(NADDCL1.LT.IADDCL1) GOTO 101
C
 102  CONTINUE
      IF(IADDK1.EQ.0) GOTO 103
      CALL FFRAND(XXX)
      CALL FFRAND(YYY)
      CALL FFRAND(ZZZ)
      XXX   = CENTX + XBOX*(XXX-PT5)
      YYY   = CENTY + YBOX*(YYY-PT5)
      ZZZ   = CENTZ + ZBOX*(ZZZ-PT5)
      DO IFFAT = 1, NFFAT
         CX    = XXX - CORD(1,IFFAT)
         CY    = YYY - CORD(2,IFFAT)
         CZ    = ZZZ - CORD(3,IFFAT)
         PBCX  = XBOX*ANINT(CX*ONEXBOX)
         PBCY  = YBOX*ANINT(CY*ONEYBOX)
         PBCZ  = ZBOX*ANINT(CZ*ONEZBOX)
         CX    = CX - PBCX
         CY    = CY - PBCY
         CZ    = CZ - PBCZ
         R2    = CX*CX+CY*CY+CZ*CZ
C        -- CANNOT BE SMALLER THAN 4.0 ANGSTROM
         IF(R2.LT.57.14D+00) GOTO 102
      ENDDO
      NADDK1          = NADDK1 + 1
      NFFAT           = NFFAT + 1
      CORD(1,NFFAT)   = XXX
      CORD(2,NFFAT)   = YYY
      CORD(3,NFFAT)   = ZZZ
      ZANF  (NFFAT)   = 19.0D+00
      ATMNAM(NFFAT)   = 'K'
      IF(NFFTYP/10000.LE.9)WRITE(RNAME,'(A4,I1,A3)')
     *                     ' $K1',NFFTYP/10000,'   '
      IF(NFFTYP/10000.GT.9)WRITE(RNAME,'(A4,I2,A2)')
     *                     ' $K1',NFFTYP/10000,'  '
      CALL RDIONS(RNAME,ZMAS,ONEMAS,
     *            CHARG,POL,SIG,EPS,
     *            SIG2,EPS2,
     *            CLPR,ZLPR,NLPR,NFFAT)
      IF(NADDK1.LT.IADDK1) GOTO 102
C
 103  CONTINUE
      IF(IADDCA2.EQ.0) GOTO 104
      CALL FFRAND(XXX)
      CALL FFRAND(YYY)
      CALL FFRAND(ZZZ)
      XXX   = CENTX + XBOX*(XXX-PT5)
      YYY   = CENTY + YBOX*(YYY-PT5)
      ZZZ   = CENTZ + ZBOX*(ZZZ-PT5)
      DO IFFAT = 1, NFFAT
         CX    = XXX - CORD(1,IFFAT)
         CY    = YYY - CORD(2,IFFAT)
         CZ    = ZZZ - CORD(3,IFFAT)
         PBCX  = XBOX*ANINT(CX*ONEXBOX)
         PBCY  = YBOX*ANINT(CY*ONEYBOX)
         PBCZ  = ZBOX*ANINT(CZ*ONEZBOX)
         CX    = CX - PBCX
         CY    = CY - PBCY
         CZ    = CZ - PBCZ
         R2    = CX*CX+CY*CY+CZ*CZ
C        -- CANNOT BE SMALLER THAN 4.0 ANGSTROM
         IF(R2.LT.57.14D+00) GOTO 103
      ENDDO
      NADDCA2         = NADDCA2+ 1
      NFFAT           = NFFAT + 1
      CORD(1,NFFAT)   = XXX
      CORD(2,NFFAT)   = YYY
      CORD(3,NFFAT)   = ZZZ
      ZANF  (NFFAT)   = 20.0D+00
      ATMNAM(NFFAT)   = 'CA'
      IF(NFFTYP/10000.LE.9)WRITE(RNAME,'(A5,I1,A2)')
     *                     ' $CA2',NFFTYP/10000,'  '
      IF(NFFTYP/10000.GT.9)WRITE(RNAME,'(A5,I2,A1)')
     *                     ' $CA2',NFFTYP/10000,' '
      CALL RDIONS(RNAME,ZMAS,ONEMAS,
     *            CHARG,POL,SIG,EPS,
     *            SIG2,EPS2,
     *            CLPR,ZLPR,NLPR,NFFAT)
      IF(NADDCA2.LT.IADDCA2) GOTO 103
C
 104  CONTINUE
      IF(IADDMG2.EQ.0) GOTO 105
      CALL FFRAND(XXX)
      CALL FFRAND(YYY)
      CALL FFRAND(ZZZ)
      XXX   = CENTX + XBOX*(XXX-PT5)
      YYY   = CENTY + YBOX*(YYY-PT5)
      ZZZ   = CENTZ + ZBOX*(ZZZ-PT5)
      DO IFFAT = 1, NFFAT
         CX    = XXX - CORD(1,IFFAT)
         CY    = YYY - CORD(2,IFFAT)
         CZ    = ZZZ - CORD(3,IFFAT)
         PBCX  = XBOX*ANINT(CX*ONEXBOX)
         PBCY  = YBOX*ANINT(CY*ONEYBOX)
         PBCZ  = ZBOX*ANINT(CZ*ONEZBOX)
         CX    = CX - PBCX
         CY    = CY - PBCY
         CZ    = CZ - PBCZ
         R2    = CX*CX+CY*CY+CZ*CZ
C        -- CANNOT BE SMALLER THAN 4.0 ANGSTROM
         IF(R2.LT.57.14D+00) GOTO 104
      ENDDO
      NADDMG2         = NADDMG2+ 1
      NFFAT           = NFFAT + 1
      CORD(1,NFFAT)   = XXX
      CORD(2,NFFAT)   = YYY
      CORD(3,NFFAT)   = ZZZ
      ZANF  (NFFAT)   = 12.0D+00
      ATMNAM(NFFAT)   = 'MG'
      IF(NFFTYP/10000.LE.9)WRITE(RNAME,'(A5,I1,A2)')
     *                     ' $MG2',NFFTYP/10000,'  '
      IF(NFFTYP/10000.GT.9)WRITE(RNAME,'(A5,I2,A1)')
     *                     ' $MG2',NFFTYP/10000,' '
      CALL RDIONS(RNAME,ZMAS,ONEMAS,
     *            CHARG,POL,SIG,EPS,
     *            SIG2,EPS2,
     *            CLPR,ZLPR,NLPR,NFFAT)
      IF(NADDMG2.LT.IADDMG2) GOTO 104
C
 105  CONTINUE
C
C     -- MUST UPDATE NFFAT0
      NFFAT0 = NFFAT
C
C     -- CHECK EACH WATER MOLECULE --
C
      LPTS = ITYPWAT/100
      NWATER = 0
      DO 200 IWAT = 1, KKKWAT
         O1(1) = WATO1(1,IWAT) + XMIN
         O1(2) = WATO1(2,IWAT) + YMIN
         O1(3) = WATO1(3,IWAT) + ZMIN
         H2(1) = WATH2(1,IWAT) + XMIN
         H2(2) = WATH2(2,IWAT) + YMIN
         H2(3) = WATH2(3,IWAT) + ZMIN
         H3(1) = WATH3(1,IWAT) + XMIN
         H3(2) = WATH3(2,IWAT) + YMIN
         H3(3) = WATH3(3,IWAT) + ZMIN
C        - IF A BIG BOX, USE SMALLER BOX
         IF(KKKWAT.GT.4096*8) THEN
            IF(O1(1).GT.(XMAX-3.0D+00).OR.
     *         O1(2).GT.(YMAX-3.0D+00).OR.
     *         O1(3).GT.(ZMAX-3.0D+00)    )THEN
               GOTO 200
            ELSE
               IF(NFFAT0.EQ.0.AND.NAT.EQ.0) THEN
                  NEAR = 0
                  GOTO 221
               END IF
            END IF
         ELSE
            IF(O1(1).GT.(XMAX+9.0D+00).OR.
     *         O1(2).GT.(YMAX+9.0D+00).OR.
     *         O1(3).GT.(ZMAX+9.0D+00))THEN
               GOTO 200
            END IF
         END IF
C        - WATER TOO CLOSE TO PROTEIN ATOMS AND ALREADY ADDED WATERS
C          IS EXCLUDED, SO USE NFFAT
         NEAR = 0
         DO 210 IFFAT = 1, NFFAT
            IF(ZANF(IFFAT).GT.1.0001D+00) THEN
               CX    = O1(1) - CORD(1,IFFAT)
               CY    = O1(2) - CORD(2,IFFAT)
               CZ    = O1(3) - CORD(3,IFFAT)
               PBCX  = XBOX*ANINT(CX*ONEXBOX)
               PBCY  = YBOX*ANINT(CY*ONEYBOX)
               PBCZ  = ZBOX*ANINT(CZ*ONEZBOX)
               CX    = CX - PBCX
               CY    = CY - PBCY
               CZ    = CZ - PBCZ
               IF(ABS(CX).GT.5.68D+00) GOTO 210
               IF(ABS(CY).GT.5.68D+00) GOTO 210
               IF(ABS(CZ).GT.5.68D+00) GOTO 210
               R2    = CX*CX+CY*CY+CZ*CZ
               DUM   = 32.14D+00   !  3.0 A
               IF(IFFAT.LE.NFFAT0.AND.R2.LT.DUM) NEAR=NEAR+1
C              -- WATER-WATER CANNOT BE SMALLER THAN 2.4 ANGSTROM
               IF(IFFAT.GT.NFFAT0.AND.R2.LT.20.57D+00) NEAR=NEAR+1
            END IF
 210     CONTINUE
         DO 220 IAT = 1, NAT
            IF(ZAN(IAT).NE.1.0D+00) THEN
               CX    = O1(1) - C(1,IAT)
               CY    = O1(2) - C(2,IAT)
               CZ    = O1(3) - C(3,IAT)
               PBCX  = XBOX*ANINT(CX*ONEXBOX)
               PBCY  = YBOX*ANINT(CY*ONEYBOX)
               PBCZ  = ZBOX*ANINT(CZ*ONEZBOX)
               CX    = CX - PBCX
               CY    = CY - PBCY
               CZ    = CZ - PBCZ
               IF(ABS(CX).GT.5.68D+00) GOTO 220
               IF(ABS(CY).GT.5.68D+00) GOTO 220
               IF(ABS(CZ).GT.5.68D+00) GOTO 220
               R2    = CX*CX+CY*CY+CZ*CZ
C              -- CANNOT BE SMALLER THAN 3.0 ANGSTROM
               IF(R2.LT.32.14D+00) NEAR=NEAR+1
            END IF
 220     CONTINUE
 221     CONTINUE
         IF(NEAR.GT.0) GOTO 200
C        - WATER NOT CLOSE TO PROTEIN ATOMS MUST STAY
         NWATER= NWATER+ 1
         NFFAT = NFFAT + 1
         CORD(1,NFFAT) = O1(1)
         CORD(2,NFFAT) = O1(2)
         CORD(3,NFFAT) = O1(3)
         ZANF  (NFFAT) = 8.0D+00
         NFFAT = NFFAT + 1
         CORD(1,NFFAT) = H2(1)
         CORD(2,NFFAT) = H2(2)
         CORD(3,NFFAT) = H2(3)
         ZANF  (NFFAT) = 1.0D+00
         NFFAT = NFFAT + 1
         CORD(1,NFFAT) = H3(1)
         CORD(2,NFFAT) = H3(2)
         CORD(3,NFFAT) = H3(3)
         ZANF  (NFFAT) = 1.0D+00
         IF(LPTS.EQ.5) THEN
            X12   = O1(1) - H2(1)
            Y12   = O1(2) - H2(2)
            Z12   = O1(3) - H2(3)
            X13   = O1(1) - H3(1)
            Y13   = O1(2) - H3(2)
            Z13   = O1(3) - H3(3)
            X14   = Y12*Z13 - Z12*Y13
            Y14   = Z12*X13 - X12*Z13
            Z14   = X12*Y13 - Y12*X13
            R14   = SQRT(X14*X14 + Y14*Y14 + Z14*Z14)
            ONER14= 1.0D+00/R14
            X14   = X14*ONER14
            Y14   = Y14*ONER14
            Z14   = Z14*ONER14
            X1T   = X12+X13
            Y1T   = Y12+Y13
            Z1T   = Z12+Z13
            R1T   = SQRT(X1T*X1T + Y1T*Y1T + Z1T*Z1T)
            ONER1T= 1.0D+00/R1T
            X1T   = X1T*ONER1T
            Y1T   = Y1T*ONER1T
            Z1T   = Z1T*ONER1T
            NFFAT = NFFAT + 1
            CORD(1,NFFAT) = O1(1) + X14 + X1T
            CORD(2,NFFAT) = O1(2) + Y14 + Y1T
            CORD(3,NFFAT) = O1(3) + Z14 + Z1T
            ZANF  (NFFAT) = 1.0D+00
            NFFAT = NFFAT + 1
            CORD(1,NFFAT) = O1(1) - X14 + X1T
            CORD(2,NFFAT) = O1(2) - Y14 + Y1T
            CORD(3,NFFAT) = O1(3) - Z14 + Z1T
            ZANF  (NFFAT) = 1.0D+00
         END IF
 200  CONTINUE
C
      IF(NFFAT.GT.MXFFAT-10) THEN
         IF(MASWRK) WRITE(IW,*)
     *   'ERROR: TOO MANY NFFAT ATOMS IN ADDPBCWAT. INCREASE MXFFAT'
         IF(MASWRK)WRITE(IW,*)' '
         CALL ABRT
      END IF
C
      LPTS = ITYPWAT/100
      IF(LPTS.EQ.3) LREAL = 2
      IF(LPTS.EQ.5) LREAL = 4
      IF(LPTS.EQ.3) LBON = 3
      IF(LPTS.EQ.5) LBON = 9
      IF(LPTS.EQ.3) LANG = 1
      IF(LPTS.EQ.5) LANG = 6
      IF(LPTS.EQ.3) L123 = 3
      IF(LPTS.EQ.5) L123 =10
C
      DO IFFAT = NFFAT0+1, NFFAT-LPTS+1, LPTS
         IF(IFFAT.EQ.NFFAT0+1) THEN
            CALL RDHOH(ATMNAM,ZMAS,ONEMAS,CHARG,POL,SIG,EPS,
     *                 SIG2,EPS2,BOND0,FCBOND,
     *                 ANGL0,FCANGL,IPAIR,KLIST,
     *                 CLPR,ZLPR,NLPR,NFFAT0+1,L1213J,ITYPWAT,
     *                 LSTRAT,DSTRAT)
         ELSE
            DO KK = 0, LPTS-1
               ATMNAM(IFFAT+KK) = ATMNAM(IFFAT+KK-LPTS)
               ZMAS  (IFFAT+KK) = ZMAS(IFFAT+KK-LPTS)
               ONEMAS(IFFAT+KK) = ONEMAS(IFFAT+KK-LPTS)
               CHARG (IFFAT+KK) = CHARG(IFFAT+KK-LPTS)
               POL   (IFFAT+KK) = POL(IFFAT+KK-LPTS)
               SIG   (IFFAT+KK) = SIG(IFFAT+KK-LPTS)
               EPS   (IFFAT+KK) = EPS(IFFAT+KK-LPTS)
               SIG2  (IFFAT+KK) = SIG2(IFFAT+KK-LPTS)
               EPS2  (IFFAT+KK) = EPS2(IFFAT+KK-LPTS)
               CLPR(1,IFFAT+KK) = CLPR(1,IFFAT+KK-LPTS)
               ZLPR(1,IFFAT+KK) = ZLPR(1,IFFAT+KK-LPTS)
               CLPR(2,IFFAT+KK) = CLPR(2,IFFAT+KK-LPTS)
               ZLPR(2,IFFAT+KK) = ZLPR(2,IFFAT+KK-LPTS)
               CLPR(3,IFFAT+KK) = CLPR(3,IFFAT+KK-LPTS)
               ZLPR(3,IFFAT+KK) = ZLPR(3,IFFAT+KK-LPTS)
               CLPR(4,IFFAT+KK) = CLPR(4,IFFAT+KK-LPTS)
               ZLPR(4,IFFAT+KK) = ZLPR(4,IFFAT+KK-LPTS)
               NLPR  (IFFAT+KK) = 2
            ENDDO
C
            DO KK=1, LBON
               NBOND = NBOND + 1
               IPAIR(1,NBOND)  = IPAIR(1,NBOND-LBON)+LPTS
               IPAIR(2,NBOND)  = IPAIR(2,NBOND-LBON)+LPTS
               FCBOND(NBOND)   = FCBOND(NBOND-LBON)
               BOND0(NBOND)    = BOND0(NBOND-LBON)
               IF((IRATTLE.GT.0.AND.KK.LE.LREAL)  .OR.
     *           ((IRATTLE.EQ.10.OR.IRATTLE.EQ.20).AND.
     *             KK.GT.LREAL))THEN
                  NRATTLE = NRATTLE + 1
                  LSTRAT(1,NRATTLE)=IPAIR(1,NBOND)
                  LSTRAT(2,NRATTLE)=IPAIR(2,NBOND)
                  DSTRAT(NRATTLE)  =BOND0(NBOND)*BOND0(NBOND)
               END IF
            ENDDO
C
            DO KK=1, LANG
               NANGL = NANGL + 1
               KLIST(1,NANGL)  = KLIST(1,NANGL-LANG)+LPTS
               KLIST(2,NANGL)  = KLIST(2,NANGL-LANG)+LPTS
               KLIST(3,NANGL)  = KLIST(3,NANGL-LANG)+LPTS
               FCANGL (NANGL)  = FCANGL (NANGL-LANG)
               ANGL0  (NANGL)  = ANGL0  (NANGL-LANG)
            ENDDO
C
            DO KK=1, L123
               N1213J = N1213J + 1
               L1213J(1,N1213J) = L1213J(1,N1213J-L123)+LPTS
               L1213J(2,N1213J) = L1213J(2,N1213J-L123)+LPTS
            ENDDO
         END IF
      ENDDO
      IF(NBOND.GT.MXBOND) THEN
         IF(MASWRK) WRITE(IW,*)
     *   'ERROR: TOO MANY BONDS IN ADDPBCWAT. INCREASE MXBOND'
         IF(MASWRK)WRITE(IW,*)' '
         CALL ABRT
      END IF
      IF(NANGL.GT.MXANGL) THEN
         IF(MASWRK) WRITE(IW,*)
     *   'ERROR: TOO MANY ANGLES IN ADDPBCWAT. INCREASE MXANGL'
         IF(MASWRK)WRITE(IW,*)' '
         CALL ABRT
      END IF
C
      RETURN
      END
C*MODULE QUANPOB  *DECK ADDSPHWAT
!>
!> @brief    Add water molecules to a sphere
!>
!> @author   Hui Li
!>           - Mar 2011
!>
!> @details  add ions first
!>           water parameters from external file
!>
      SUBROUTINE ADDSPHWAT(ATMNAM,CORD,ZANF,ZMAS,ONEMAS,CHARG,POL,
     *                     SIG,EPS,SIG2,EPS2,CLPR,ZLPR,NLPR,
     *                     IPAIR,FCBOND,BOND0,
     *                     KLIST,FCANGL,ANGL0,L1213J,
     *                     MXFFAT,MXBOND,MXANGL,
     *                     WATO1,WATH2,WATH3,NWATER,ITYPWAT,
     *                     IADDNA1,IADDK1,IADDCA2,IADDMG2,IADDCL1,
     *                     JADDNA1,JADDK1,LSTRAT,DSTRAT)
      use mx_limits, only: mxatm
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      LOGICAL GOPARR,MASWRK,DSKWRK
C
      CHARACTER*8  RNAME
      CHARACTER*10 ATMNAM
      CHARACTER*256 QPFILE
C
      PARAMETER (ZERO=0.0D+00)
      PARAMETER (PT5=0.5D+00)
      PARAMETER (TWO=2.0D+00)
      PARAMETER (FOUR=4.0D+00)
      PARAMETER (TOANGS=0.52917724924D+00)
      PARAMETER (TOBOHR=1.0D+00/TOANGS)
C
      DIMENSION ATMNAM(*),CORD(3,*),ZANF(*),ZMAS(*),
     *          ONEMAS(*),CHARG(*),POL(*),
     *          SIG(*),EPS(*),SIG2(*),EPS2(*),CLPR(4,*),ZLPR(4,*),
     *          NLPR(*),IPAIR(2,*),FCBOND(*),BOND0(*),
     *          KLIST(3,*),FCANGL(*),ANGL0(*),L1213J(2,*),
     *          WATO1(3,*),WATH2(3,*),WATH3(3,*),
     *          O1(3),H2(3),H3(3),LSTRAT(2,*),DSTRAT(*)
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
      COMMON /FFRATT/ RATOLC,RATOLV,SCALRAT,VIRRAT(3),IRATTLE,JRATTLE,
     *                NRATTLE,MXRATT,LFFOLDCORD,LFFLSTRAT,LFFDSTRAT,
     *                LFFVELSV,IRATQM
      COMMON /FFSPH / SPHRAD,SPHSIG,SPHEPS,IADDWAT
      COMMON /FFTYPE/ WT14LJ,WT14CH,C3BOND,C4BOND,C3ANGL,
     *                NFFTYP,NFFFILE,LJQMMM,LJQM,INTCHG,
     *                LJSIGMA,JTOPFILE(90),JPARFILE(90),
     *                JTOPAMIA(90),JTOPNTER(90),JTOPCTER(90),
     *                JTOPNUCA(90),JPARFIL2(90),JPARFIL3(90)
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /IOFILE/ IR,IW,IP,IJK,IJKT,IDAF,NAV,IODA(950)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
C
C     ADD WATER MOLECULES TO THE SYSTEM TO FILL A SPHERICAL SPACE
C     WITH R=SPHRAD, AND CENTERTED AT (CENTX,CENTY,CENTZ).
C     NOTE: SPHERICAL SYSTEMS HAVE 10E+6 BAR PRESSURE (SURFACE TENSION).
C           ISPHSOL MAY REDUCE THE PRESSURE.
C
C     HUI LI, MAR 2011, LINCOLN
C
      NFFAT0 = NFFAT
C
      RRMOL   = ZERO
      DO IAT = 1, NAT
         DX      = C(1,IAT) - CENTX
         DY      = C(2,IAT) - CENTY
         DZ      = C(3,IAT) - CENTZ
         RR      = DX*DX+DY*DY+DZ*DZ
         RRMOL   = MAX(RRMOL,RR)
      ENDDO
      DO IFFAT = 1,NFFAT0
         DX      = CORD(1,IFFAT) - CENTX
         DY      = CORD(2,IFFAT) - CENTY
         DZ      = CORD(3,IFFAT) - CENTZ
         RR      = DX*DX+DY*DY+DZ*DZ
         RRMOL   = MAX(RRMOL,RR)
      ENDDO
      RMOL = SQRT(RRMOL)
      IF(MASWRK) THEN
         WRITE(IW,'(1X,A,F8.2," ANGSTROM")')
     *     'RMOL=',RMOL*TOANGS
      END IF
C     - THE AVERAGE DISTANCE BETWEEN AN ATOM
C       AND WATER OXYGEN ATOM SHOULD BE 2.3 ANGSTROM
C       THIS VALUE DOES NOT AFFECT THE NUMBER OF WATER
C       BECAUSE IT IS USED FOR PRE-SCREENING
      VDW     = 2.3D+00*TOBOHR
      RMOL    = RMOL + VDW
      IF((SPHRAD-RMOL).LT.ZERO)THEN
         IF(MASWRK) WRITE(IW,*)'ERROR: SPHRAD MUST BE LARGER THAN ',
     *      RMOL*TOANGS,' ANGSTROM.'
         IF(MASWRK)WRITE(IW,*)' '
         CALL ABRT
      END IF
C     - DETERMINE THE BOUNDARY OF THE WATER BOX
      XMAX    = CENTX + SPHRAD
      YMAX    = CENTY + SPHRAD
      ZMAX    = CENTZ + SPHRAD
      XMIN    = CENTX - SPHRAD
      YMIN    = CENTY - SPHRAD
      ZMIN    = CENTZ - SPHRAD
C     - DETERMINE THE MAXIMUM NUMBER OF WATER
      VOLM    = (TWO*SPHRAD)**3
C     - VOLUME OF ONE WATER =  29.998696 A**3 AT 298.15 K, 1 BAR
C                           = 202.441191 BOHR**3
      MAXWAT  = INT(VOLM/202.441191D+00)
      NFOLD   = 0
      DO K = 0, 30, 3
         IF(MAXWAT.GT.4096*(2**K)) NFOLD = K+3
      ENDDO
      KKKWAT  = 4096*(2**NFOLD)
      IF(KKKWAT.GT.8*MXFFAT) THEN
         IF(MASWRK) WRITE(IW,*)
     *      'ERROR: TOO MANY WATER IN ADDSPHWAT. INCREASE MXFFAT'
         IF(MASWRK)WRITE(IW,*)' '
         CALL ABRT
      END IF
C
C     -- READ IN 4096 WATER MOLECULES --
C        THEY ARE IN A CUBE WITH SIDE LENGTH = 49.7150 A.
C        CREATE MORE WATER IF NECESSARY
C
      IF(MASWRK) THEN
      CALL GENQPFILE('WATER4096.DAT ',QPFILE,LENQP)
      OPEN(12,FILE=QPFILE(1:LENQP),STATUS='OLD')
      DO IWAT=1,4096
         READ(12,*) WATO1(1,IWAT),WATO1(2,IWAT),WATO1(3,IWAT),
     *              WATH2(1,IWAT),WATH2(2,IWAT),WATH2(3,IWAT),
     *              WATH3(1,IWAT),WATH3(2,IWAT),WATH3(3,IWAT)
      ENDDO
      CLOSE(12)
      END IF
      IF(GOPARR) CALL DDI_BCAST(457,'F',WATO1,3*4096,MASTER)
      IF(GOPARR) CALL DDI_BCAST(458,'F',WATH2,3*4096,MASTER)
      IF(GOPARR) CALL DDI_BCAST(459,'F',WATH3,3*4096,MASTER)
      SIZE  = 49.7150D+00*TOBOHR*PT5
      DO IFOLD = 1, NFOLD
         IXYZ1 = MOD(IFOLD,3)
         IF(IXYZ1.EQ.1) IXYZ2 = 2
         IF(IXYZ1.EQ.1) IXYZ3 = 3
         IF(IXYZ1.EQ.2) IXYZ2 = 1
         IF(IXYZ1.EQ.2) IXYZ3 = 3
         IF(IXYZ1.EQ.0) IXYZ2 = 1
         IF(IXYZ1.EQ.0) IXYZ3 = 2
         IF(IXYZ1.EQ.0) IXYZ1 = 3
         IF(IXYZ1.EQ.1) SIZE  = SIZE*2.0D+00
         NDONE  = 2**(IFOLD-1)
         LENGTH = NDONE*4096
         DO III = 1, LENGTH
            WATO1(IXYZ1,LENGTH+III) = WATO1(IXYZ1,III) + SIZE
            WATO1(IXYZ2,LENGTH+III) = WATO1(IXYZ2,III)
            WATO1(IXYZ3,LENGTH+III) = WATO1(IXYZ3,III)
            WATH2(IXYZ1,LENGTH+III) = WATH2(IXYZ1,III) + SIZE
            WATH2(IXYZ2,LENGTH+III) = WATH2(IXYZ2,III)
            WATH2(IXYZ3,LENGTH+III) = WATH2(IXYZ3,III)
            WATH3(IXYZ1,LENGTH+III) = WATH3(IXYZ1,III) + SIZE
            WATH3(IXYZ2,LENGTH+III) = WATH3(IXYZ2,III)
            WATH3(IXYZ3,LENGTH+III) = WATH3(IXYZ3,III)
         ENDDO
      ENDDO
C
C     -- ADD IONS BEFORE WATER --
C        (1) ADD NA+ OR K+ TO DNA/RNA, IF REQUESTED.
C        (2) ADD OTHER IONS
C
      MADDNA1= 0
      MADDK1 = 0
      IF(JADDNA1.EQ.1 .OR. JADDK1.EQ.1) THEN
      DO IFFAT = 1, NFFAT
         IF(ATMNAM(IFFAT  ).EQ.'P         '.AND.
     *      ATMNAM(IFFAT+1).EQ.'O         '.AND.
     *      ATMNAM(IFFAT+2).EQ.'O         ') THEN
            DIS1 = (CORD(1,IFFAT)-CORD(1,IFFAT+1))**2 +
     *             (CORD(2,IFFAT)-CORD(2,IFFAT+1))**2 +
     *             (CORD(3,IFFAT)-CORD(3,IFFAT+1))**2
            DIS2 = (CORD(1,IFFAT)-CORD(1,IFFAT+2))**2 +
     *             (CORD(2,IFFAT)-CORD(2,IFFAT+2))**2 +
     *             (CORD(3,IFFAT)-CORD(3,IFFAT+2))**2
            IF(DIS1.LT.11.6D+00 .AND. DIS2.LT.11.6D+00) THEN
               XXX = CORD(1,IFFAT)
     *             + ((CORD(1,IFFAT+1)+CORD(1,IFFAT+2))*PT5
     *                -CORD(1,IFFAT))*FOUR
               YYY = CORD(2,IFFAT)
     *             + ((CORD(2,IFFAT+1)+CORD(2,IFFAT+2))*PT5
     *                -CORD(2,IFFAT))*FOUR
               ZZZ = CORD(3,IFFAT)
     *             + ((CORD(3,IFFAT+1)+CORD(3,IFFAT+2))*PT5
     *                -CORD(3,IFFAT))*FOUR
               IF(JADDNA1.EQ.1) THEN
               MADDNA1         = MADDNA1 + 1
               NFFAT           = NFFAT + 1
               CORD(1,NFFAT)   = XXX
               CORD(2,NFFAT)   = YYY
               CORD(3,NFFAT)   = ZZZ
               ZANF  (NFFAT)   = 11.0D+00
               ATMNAM(NFFAT)   = 'NA'
               IF(NFFTYP/10000.LE.9)WRITE(RNAME,'(A5,I1,A2)')
     *                        ' $NA1',NFFTYP/10000,'  '
               IF(NFFTYP/10000.GT.9)WRITE(RNAME,'(A5,I2,A1)')
     *                        ' $NA1',NFFTYP/10000,' '
               CALL RDIONS(RNAME,ZMAS,ONEMAS,
     *                     CHARG,POL,SIG,EPS,
     *                     SIG2,EPS2,
     *                     CLPR,ZLPR,NLPR,NFFAT)
               ELSE IF(JADDK1.EQ.1) THEN
               MADDK1          = MADDK1 + 1
               NFFAT           = NFFAT + 1
               CORD(1,NFFAT)   = XXX
               CORD(2,NFFAT)   = YYY
               CORD(3,NFFAT)   = ZZZ
               ZANF  (NFFAT)   = 19.0D+00
               ATMNAM(NFFAT)   = 'K'
               IF(NFFTYP/10000.LE.9)WRITE(RNAME,'(A4,I1,A3)')
     *                        ' $K1',NFFTYP/10000,'   '
               IF(NFFTYP/10000.GT.9)WRITE(RNAME,'(A4,I2,A2)')
     *                        ' $K1',NFFTYP/10000,'  '
               CALL RDIONS(RNAME,ZMAS,ONEMAS,
     *                     CHARG,POL,SIG,EPS,
     *                     SIG2,EPS2,
     *                     CLPR,ZLPR,NLPR,NFFAT)
               END IF
            END IF
         END IF
      ENDDO
      JADDNA1 = MADDNA1
      JADDK1  = MADDK1
      END IF
C
      NADDNA1= 0
      NADDK1 = 0
      NADDMG2= 0
      NADDCA2= 0
      NADDCL1= 0
C
 100  CONTINUE
      IF(IADDNA1.EQ.0) GOTO 101
      CALL FFRAND(XXX)
      CALL FFRAND(YYY)
      CALL FFRAND(ZZZ)
      XXX   = TWO*SPHRAD*(XXX-PT5)
      YYY   = TWO*SPHRAD*(YYY-PT5)
      ZZZ   = TWO*SPHRAD*(ZZZ-PT5)
C     - IONS BEYOND THE SPHRAD-RWAT LIMIT IS EXCLUDED
      RR = XXX**2 + YYY**2 + ZZZ**2
      R  = SQRT(RR)
      IF(R.GT.(SPHRAD-2.936D+00)) GOTO 100
      XXX   = CENTX + XXX
      YYY   = CENTY + YYY
      ZZZ   = CENTZ + ZZZ
      DO IFFAT = 1, NFFAT
         CX    = XXX - CORD(1,IFFAT)
         CY    = YYY - CORD(2,IFFAT)
         CZ    = ZZZ - CORD(3,IFFAT)
         R2    = CX*CX+CY*CY+CZ*CZ
C        -- CANNOT BE SMALLER THAN 4.0 ANGSTROM
         IF(R2.LT.57.14D+00) GOTO 100
      ENDDO
      NADDNA1         = NADDNA1 + 1
      NFFAT           = NFFAT + 1
      CORD(1,NFFAT)   = XXX
      CORD(2,NFFAT)   = YYY
      CORD(3,NFFAT)   = ZZZ
      ZANF  (NFFAT)   = 11.0D+00
      ATMNAM(NFFAT)   = 'NA'
      IF(NFFTYP/10000.LE.9)WRITE(RNAME,'(A5,I1,A2)')
     *                     ' $NA1',NFFTYP/10000,'  '
      IF(NFFTYP/10000.GT.9)WRITE(RNAME,'(A5,I2,A1)')
     *                     ' $NA1',NFFTYP/10000,' '
      CALL RDIONS(RNAME,ZMAS,ONEMAS,
     *            CHARG,POL,SIG,EPS,
     *            SIG2,EPS2,
     *            CLPR,ZLPR,NLPR,NFFAT)
      IF(NADDNA1.LT.IADDNA1) GOTO 100
C
 101  CONTINUE
      IF(IADDCL1.EQ.0) GOTO 102
      CALL FFRAND(XXX)
      CALL FFRAND(YYY)
      CALL FFRAND(ZZZ)
      XXX   = TWO*SPHRAD*(XXX-PT5)
      YYY   = TWO*SPHRAD*(YYY-PT5)
      ZZZ   = TWO*SPHRAD*(ZZZ-PT5)
C     - IONS BEYOND THE SPHRAD-RWAT LIMIT IS EXCLUDED
      RR = XXX**2 + YYY**2 + ZZZ**2
      R  = SQRT(RR)
      IF(R.GT.(SPHRAD-2.936D+00)) GOTO 101
      XXX   = CENTX + XXX
      YYY   = CENTY + YYY
      ZZZ   = CENTZ + ZZZ
      DO IFFAT = 1, NFFAT
         CX    = XXX - CORD(1,IFFAT)
         CY    = YYY - CORD(2,IFFAT)
         CZ    = ZZZ - CORD(3,IFFAT)
         R2    = CX*CX+CY*CY+CZ*CZ
C        -- CANNOT BE SMALLER THAN 4.0 ANGSTROM
         IF(R2.LT.57.14D+00) GOTO 101
      ENDDO
      NADDCL1         = NADDCL1 + 1
      NFFAT           = NFFAT + 1
      CORD(1,NFFAT)   = XXX
      CORD(2,NFFAT)   = YYY
      CORD(3,NFFAT)   = ZZZ
      ZANF  (NFFAT)   = 17.0D+00
      ATMNAM(NFFAT)   = 'CL'
      IF(NFFTYP/10000.LE.9)WRITE(RNAME,'(A5,I1,A2)')
     *                     ' $CL1',NFFTYP/10000,'  '
      IF(NFFTYP/10000.GT.9)WRITE(RNAME,'(A5,I2,A1)')
     *                     ' $CL1',NFFTYP/10000,' '
      CALL RDIONS(RNAME,ZMAS,ONEMAS,
     *            CHARG,POL,SIG,EPS,
     *            SIG2,EPS2,
     *            CLPR,ZLPR,NLPR,NFFAT)
      IF(NADDCL1.LT.IADDCL1) GOTO 101
C
 102  CONTINUE
      IF(IADDK1.EQ.0) GOTO 103
      CALL FFRAND(XXX)
      CALL FFRAND(YYY)
      CALL FFRAND(ZZZ)
      XXX   = TWO*SPHRAD*(XXX-PT5)
      YYY   = TWO*SPHRAD*(YYY-PT5)
      ZZZ   = TWO*SPHRAD*(ZZZ-PT5)
C     - IONS BEYOND THE SPHRAD-RWAT LIMIT IS EXCLUDED
      RR = XXX**2 + YYY**2 + ZZZ**2
      R  = SQRT(RR)
      IF(R.GT.(SPHRAD-2.936D+00)) GOTO 102
      XXX   = CENTX + XXX
      YYY   = CENTY + YYY
      ZZZ   = CENTZ + ZZZ
      DO IFFAT = 1, NFFAT
         CX    = XXX - CORD(1,IFFAT)
         CY    = YYY - CORD(2,IFFAT)
         CZ    = ZZZ - CORD(3,IFFAT)
         R2    = CX*CX+CY*CY+CZ*CZ
C        -- CANNOT BE SMALLER THAN 4.0 ANGSTROM
         IF(R2.LT.57.14D+00) GOTO 102
      ENDDO
      NADDK1          = NADDK1 + 1
      NFFAT           = NFFAT + 1
      CORD(1,NFFAT)   = XXX
      CORD(2,NFFAT)   = YYY
      CORD(3,NFFAT)   = ZZZ
      ZANF  (NFFAT)   = 19.0D+00
      ATMNAM(NFFAT)   = 'K'
      IF(NFFTYP/10000.LE.9)WRITE(RNAME,'(A4,I1,A3)')
     *                     ' $K1',NFFTYP/10000,'   '
      IF(NFFTYP/10000.GT.9)WRITE(RNAME,'(A4,I2,A2)')
     *                     ' $K1',NFFTYP/10000,'  '
      CALL RDIONS(RNAME,ZMAS,ONEMAS,
     *            CHARG,POL,SIG,EPS,
     *            SIG2,EPS2,
     *            CLPR,ZLPR,NLPR,NFFAT)
      IF(NADDK1.LT.IADDK1) GOTO 102
C
 103  CONTINUE
      IF(IADDCA2.EQ.0) GOTO 104
      CALL FFRAND(XXX)
      CALL FFRAND(YYY)
      CALL FFRAND(ZZZ)
      XXX   = TWO*SPHRAD*(XXX-PT5)
      YYY   = TWO*SPHRAD*(YYY-PT5)
      ZZZ   = TWO*SPHRAD*(ZZZ-PT5)
C     - IONS BEYOND THE SPHRAD-RWAT LIMIT IS EXCLUDED
      RR = XXX**2 + YYY**2 + ZZZ**2
      R  = SQRT(RR)
      IF(R.GT.(SPHRAD-2.936D+00)) GOTO 103
      XXX   = CENTX + XXX
      YYY   = CENTY + YYY
      ZZZ   = CENTZ + ZZZ
      DO IFFAT = 1, NFFAT
         CX    = XXX - CORD(1,IFFAT)
         CY    = YYY - CORD(2,IFFAT)
         CZ    = ZZZ - CORD(3,IFFAT)
         PBCX  = XBOX*ANINT(CX*ONEXBOX)
         PBCY  = YBOX*ANINT(CY*ONEYBOX)
         PBCZ  = ZBOX*ANINT(CZ*ONEZBOX)
         CX    = CX - PBCX
         CY    = CY - PBCY
         CZ    = CZ - PBCZ
         R2    = CX*CX+CY*CY+CZ*CZ
C        -- CANNOT BE SMALLER THAN 4.0 ANGSTROM
         IF(R2.LT.57.14D+00) GOTO 103
      ENDDO
      NADDCA2         = NADDCA2+ 1
      NFFAT           = NFFAT + 1
      CORD(1,NFFAT)   = XXX
      CORD(2,NFFAT)   = YYY
      CORD(3,NFFAT)   = ZZZ
      ZANF  (NFFAT)   = 20.0D+00
      ATMNAM(NFFAT)   = 'CA'
      IF(NFFTYP/10000.LE.9)WRITE(RNAME,'(A5,I1,A2)')
     *                     ' $CA2',NFFTYP/10000,'  '
      IF(NFFTYP/10000.GT.9)WRITE(RNAME,'(A5,I2,A1)')
     *                     ' $CA2',NFFTYP/10000,' '
      CALL RDIONS(RNAME,ZMAS,ONEMAS,
     *            CHARG,POL,SIG,EPS,
     *            SIG2,EPS2,
     *            CLPR,ZLPR,NLPR,NFFAT)
      IF(NADDCA2.LT.IADDCA2) GOTO 103
C
 104  CONTINUE
      IF(IADDMG2.EQ.0) GOTO 105
      CALL FFRAND(XXX)
      CALL FFRAND(YYY)
      CALL FFRAND(ZZZ)
      XXX   = TWO*SPHRAD*(XXX-PT5)
      YYY   = TWO*SPHRAD*(YYY-PT5)
      ZZZ   = TWO*SPHRAD*(ZZZ-PT5)
C     - IONS BEYOND THE SPHRAD-RWAT LIMIT IS EXCLUDED
      RR = XXX**2 + YYY**2 + ZZZ**2
      R  = SQRT(RR)
      IF(R.GT.(SPHRAD-2.936D+00)) GOTO 104
      XXX   = CENTX + XXX
      YYY   = CENTY + YYY
      ZZZ   = CENTZ + ZZZ
      DO IFFAT = 1, NFFAT
         CX    = XXX - CORD(1,IFFAT)
         CY    = YYY - CORD(2,IFFAT)
         CZ    = ZZZ - CORD(3,IFFAT)
         PBCX  = XBOX*ANINT(CX*ONEXBOX)
         PBCY  = YBOX*ANINT(CY*ONEYBOX)
         PBCZ  = ZBOX*ANINT(CZ*ONEZBOX)
         CX    = CX - PBCX
         CY    = CY - PBCY
         CZ    = CZ - PBCZ
         R2    = CX*CX+CY*CY+CZ*CZ
C        -- CANNOT BE SMALLER THAN 4.0 ANGSTROM
         IF(R2.LT.57.14D+00) GOTO 104
      ENDDO
      NADDMG2         = NADDMG2+ 1
      NFFAT           = NFFAT + 1
      CORD(1,NFFAT)   = XXX
      CORD(2,NFFAT)   = YYY
      CORD(3,NFFAT)   = ZZZ
      ZANF  (NFFAT)   = 12.0D+00
      ATMNAM(NFFAT)   = 'MG'
      IF(NFFTYP/10000.LE.9)WRITE(RNAME,'(A5,I1,A2)')
     *                     ' $MG2',NFFTYP/10000,'  '
      IF(NFFTYP/10000.GT.9)WRITE(RNAME,'(A5,I2,A1)')
     *                     ' $MG2',NFFTYP/10000,' '
      CALL RDIONS(RNAME,ZMAS,ONEMAS,
     *            CHARG,POL,SIG,EPS,
     *            SIG2,EPS2,
     *            CLPR,ZLPR,NLPR,NFFAT)
      IF(NADDMG2.LT.IADDMG2) GOTO 104
C
 105  CONTINUE
C
C     -- MUST UPDATE NFFAT0
      NFFAT0 = NFFAT
C
C     -- CHECK EACH WATER MOLECULE --
C
      LPTS = ITYPWAT/100
      NWATER = 0
      DO 200 IWAT = 1, KKKWAT
         O1(1) = WATO1(1,IWAT) + XMIN
         O1(2) = WATO1(2,IWAT) + YMIN
         O1(3) = WATO1(3,IWAT) + ZMIN
         H2(1) = WATH2(1,IWAT) + XMIN
         H2(2) = WATH2(2,IWAT) + YMIN
         H2(3) = WATH2(3,IWAT) + ZMIN
         H3(1) = WATH3(1,IWAT) + XMIN
         H3(2) = WATH3(2,IWAT) + YMIN
         H3(3) = WATH3(3,IWAT) + ZMIN
C        - WATER BEYOND THE LIMIT IS EXCLUDED
         IF(O1(1).GT.XMAX.OR.O1(2).GT.YMAX.OR.O1(3).GT.ZMAX)THEN
            GOTO 200
         END IF
C        - WATER BEYOND THE SPHRAD LIMIT IS EXCLUDED
         RR = (O1(1)-CENTX)**2
     *      + (O1(2)-CENTY)**2
     *      + (O1(3)-CENTZ)**2
         R  = SQRT(RR)
         IF(R.GT.SPHRAD) GOTO 200
         IF(NFFAT0.EQ.0.AND.NAT.EQ.0) THEN
            NEAR = 0
            GOTO 221
         END IF
C        - WATER TOO CLOSE TO PROTEIN ATOMS IS EXCLUDED
C          NO PBC HERE, SO USE 'NFFAT0'.
         NEAR = 0
         DO 210 IFFAT = 1, NFFAT0
            IF(ZANF(IFFAT).GT.1.0001D+00) THEN
               CX    = O1(1) - CORD(1,IFFAT)
               CY    = O1(2) - CORD(2,IFFAT)
               CZ    = O1(3) - CORD(3,IFFAT)
               IF(ABS(CX).GT.5.68D+00) GOTO 210
               IF(ABS(CY).GT.5.68D+00) GOTO 210
               IF(ABS(CZ).GT.5.68D+00) GOTO 210
               R2    = CX*CX+CY*CY+CZ*CZ
C              -- CANNOT BE SMALLER THAN 3.0 ANGSTROM
               DUM   = 32.14D+00   !  3.0 A
               IF(R2.LT.DUM) NEAR=NEAR+1
            END IF
 210     CONTINUE
         DO 220 IAT = 1, NAT
            IF(ZAN(IAT).NE.1.0D+00) THEN
               CX    = O1(1) - C(1,IAT)
               CY    = O1(2) - C(2,IAT)
               CZ    = O1(3) - C(3,IAT)
               IF(ABS(CX).GT.5.68D+00) GOTO 220
               IF(ABS(CY).GT.5.68D+00) GOTO 220
               IF(ABS(CZ).GT.5.68D+00) GOTO 220
               R2    = CX*CX+CY*CY+CZ*CZ
C              -- CANNOT BE SMALLER THAN 3.0 ANGSTROM
               IF(R2.LT.32.14D+00) NEAR=NEAR+1
            END IF
 220     CONTINUE
 221     CONTINUE
         IF(NEAR.GT.0) GOTO 200
C        - WATER NOT CLOSE TO PROTEIN ATOMS MUST STAY
         NWATER= NWATER+ 1
         NFFAT = NFFAT + 1
         CORD(1,NFFAT) = O1(1)
         CORD(2,NFFAT) = O1(2)
         CORD(3,NFFAT) = O1(3)
         ZANF  (NFFAT) = 8.0D+00
         NFFAT = NFFAT + 1
         CORD(1,NFFAT) = H2(1)
         CORD(2,NFFAT) = H2(2)
         CORD(3,NFFAT) = H2(3)
         ZANF  (NFFAT) = 1.0D+00
         NFFAT = NFFAT + 1
         CORD(1,NFFAT) = H3(1)
         CORD(2,NFFAT) = H3(2)
         CORD(3,NFFAT) = H3(3)
         ZANF  (NFFAT) = 1.0D+00
         IF(LPTS.EQ.5) THEN
            X12   = O1(1) - H2(1)
            Y12   = O1(2) - H2(2)
            Z12   = O1(3) - H2(3)
            X13   = O1(1) - H3(1)
            Y13   = O1(2) - H3(2)
            Z13   = O1(3) - H3(3)
            X14   = Y12*Z13 - Z12*Y13
            Y14   = Z12*X13 - X12*Z13
            Z14   = X12*Y13 - Y12*X13
            R14   = SQRT(X14*X14 + Y14*Y14 + Z14*Z14)
            ONER14= 1.0D+00/R14
            X14   = X14*ONER14
            Y14   = Y14*ONER14
            Z14   = Z14*ONER14
            X1T   = X12+X13
            Y1T   = Y12+Y13
            Z1T   = Z12+Z13
            R1T   = SQRT(X1T*X1T + Y1T*Y1T + Z1T*Z1T)
            ONER1T= 1.0D+00/R1T
            X1T   = X1T*ONER1T
            Y1T   = Y1T*ONER1T
            Z1T   = Z1T*ONER1T
            NFFAT = NFFAT + 1
            CORD(1,NFFAT) = O1(1) + X14 + X1T
            CORD(2,NFFAT) = O1(2) + Y14 + Y1T
            CORD(3,NFFAT) = O1(3) + Z14 + Z1T
            ZANF  (NFFAT) = 1.0D+00
            NFFAT = NFFAT + 1
            CORD(1,NFFAT) = O1(1) - X14 + X1T
            CORD(2,NFFAT) = O1(2) - Y14 + Y1T
            CORD(3,NFFAT) = O1(3) - Z14 + Z1T
            ZANF  (NFFAT) = 1.0D+00
         END IF
 200  CONTINUE
C
      IF(NFFAT.GT.MXFFAT-10) THEN
         IF(MASWRK) WRITE(IW,*)
     *   'ERROR: TOO MANY NFFAT ATOMS IN ADDSPHWAT. INCREASE MXFFAT'
         IF(MASWRK)WRITE(IW,*)' '
         CALL ABRT
      END IF
C
      LPTS = ITYPWAT/100
      IF(LPTS.EQ.3) LREAL = 2
      IF(LPTS.EQ.5) LREAL = 4
      IF(LPTS.EQ.3) LBON = 3
      IF(LPTS.EQ.5) LBON = 9
      IF(LPTS.EQ.3) LANG = 1
      IF(LPTS.EQ.5) LANG = 6
      IF(LPTS.EQ.3) L123 = 3
      IF(LPTS.EQ.5) L123 =10
C
      DO IFFAT = NFFAT0+1, NFFAT-LPTS+1, LPTS
         IF(IFFAT.EQ.NFFAT0+1) THEN
            CALL RDHOH(ATMNAM,ZMAS,ONEMAS,CHARG,POL,SIG,EPS,
     *                 SIG2,EPS2,BOND0,FCBOND,
     *                 ANGL0,FCANGL,IPAIR,KLIST,
     *                 CLPR,ZLPR,NLPR,NFFAT0+1,L1213J,ITYPWAT,
     *                 LSTRAT,DSTRAT)
         ELSE
            DO KK = 0, LPTS-1
               ATMNAM(IFFAT+KK) = ATMNAM(IFFAT+KK-LPTS)
               ZMAS  (IFFAT+KK) = ZMAS(IFFAT+KK-LPTS)
               ONEMAS(IFFAT+KK) = ONEMAS(IFFAT+KK-LPTS)
               CHARG (IFFAT+KK) = CHARG(IFFAT+KK-LPTS)
               POL   (IFFAT+KK) = POL(IFFAT+KK-LPTS)
               SIG   (IFFAT+KK) = SIG(IFFAT+KK-LPTS)
               EPS   (IFFAT+KK) = EPS(IFFAT+KK-LPTS)
               SIG2  (IFFAT+KK) = SIG2(IFFAT+KK-LPTS)
               EPS2  (IFFAT+KK) = EPS2(IFFAT+KK-LPTS)
               CLPR(1,IFFAT+KK) = CLPR(1,IFFAT+KK-LPTS)
               ZLPR(1,IFFAT+KK) = ZLPR(1,IFFAT+KK-LPTS)
               CLPR(2,IFFAT+KK) = CLPR(2,IFFAT+KK-LPTS)
               ZLPR(2,IFFAT+KK) = ZLPR(2,IFFAT+KK-LPTS)
               CLPR(3,IFFAT+KK) = CLPR(3,IFFAT+KK-LPTS)
               ZLPR(3,IFFAT+KK) = ZLPR(3,IFFAT+KK-LPTS)
               CLPR(4,IFFAT+KK) = CLPR(4,IFFAT+KK-LPTS)
               ZLPR(4,IFFAT+KK) = ZLPR(4,IFFAT+KK-LPTS)
               NLPR  (IFFAT+KK) = 2
            ENDDO
C
            DO KK=1, LBON
               NBOND = NBOND + 1
               IPAIR(1,NBOND)  = IPAIR(1,NBOND-LBON)+LPTS
               IPAIR(2,NBOND)  = IPAIR(2,NBOND-LBON)+LPTS
               FCBOND (NBOND)  = FCBOND (NBOND-LBON)
               BOND0  (NBOND)  = BOND0  (NBOND-LBON)
               IF((IRATTLE.GT.0.AND.KK.LE.LREAL)  .OR.
     *           ((IRATTLE.EQ.10.OR.IRATTLE.EQ.20).AND.
     *             KK.GT.LREAL))THEN
                  NRATTLE = NRATTLE + 1
                  LSTRAT(1,NRATTLE)=IPAIR(1,NBOND)
                  LSTRAT(2,NRATTLE)=IPAIR(2,NBOND)
                  DSTRAT(NRATTLE)  =BOND0(NBOND)*BOND0(NBOND)
               END IF
            ENDDO
C
            DO KK=1, LANG
               NANGL = NANGL + 1
               KLIST(1,NANGL)  = KLIST(1,NANGL-LANG)+LPTS
               KLIST(2,NANGL)  = KLIST(2,NANGL-LANG)+LPTS
               KLIST(3,NANGL)  = KLIST(3,NANGL-LANG)+LPTS
               FCANGL (NANGL)  = FCANGL (NANGL-LANG)
               ANGL0  (NANGL)  = ANGL0  (NANGL-LANG)
            ENDDO
C
            DO KK=1, L123
               N1213J = N1213J + 1
               L1213J(1,N1213J) = L1213J(1,N1213J-L123)+LPTS
               L1213J(2,N1213J) = L1213J(2,N1213J-L123)+LPTS
            ENDDO
         END IF
      ENDDO
      IF(NBOND.GT.MXBOND) THEN
         IF(MASWRK) WRITE(IW,*)
     *   'ERROR: TOO MANY BONDS IN ADDSPHWAT. INCREASE MXBOND'
         IF(MASWRK)WRITE(IW,*)' '
         CALL ABRT
      END IF
      IF(NANGL.GT.MXANGL) THEN
         IF(MASWRK) WRITE(IW,*)
     *   'ERROR: TOO MANY ANGLES IN ADDSPHWAT. INCREASE MXANGL'
         IF(MASWRK)WRITE(IW,*)' '
         CALL ABRT
      END IF
C
      RETURN
      END
C*MODULE QUANPOB  *DECK RDIONS
      SUBROUTINE RDIONS(RNAME,ZMAS,ONEMAS,
     *                  CHARG,POL,SIG,EPS,
     *                  SIG2,EPS2,
     *                  CLPR,ZLPR,NLPR,IFFAT)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      LOGICAL GOPARR,DSKWRK,MASWRK
C
      PARAMETER (TOANGS=0.52917724924D+00)
      PARAMETER (TOBOHR=1.0D+00/TOANGS)
      PARAMETER (TOKCAL=627.509469D+00)
      PARAMETER (TOHART=1.0D+00/TOKCAL)
C
      CHARACTER*8 RNAME
      CHARACTER*10 WORD,FFNAME
      CHARACTER*256 QPFILE
C
      DIMENSION ZMAS(*),ONEMAS(*),CHARG(*),POL(*),
     *          SIG(*),EPS(*),
     *          SIG2(*),EPS2(*),
     *          CLPR(4,*),ZLPR(4,*),NLPR(*)
C
      COMMON /IOFILE/ IR,IW,IP,IJK,IJKT,IDAF,NAV,IODA(950)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
C
C     HUI LI, MAY 2011, LINCOLN
C
C     -- READ IN IONS PARAMETERS FROM WATERIONS.LIB --
C
      IF(MASWRK) THEN
         CALL GENQPFILE('WATERIONS.LIB ',QPFILE,LENQP)
         OPEN(12,FILE=QPFILE(1:LENQP),STATUS='OLD')
      END IF
      IEOF=0
      CALL SEQREW(12)
      CALL FNDGRP(12,RNAME,IEOF)
      IF(IEOF.EQ.1) THEN
         IF(MASWRK) WRITE(IW,*)
     *   'ERROR: ',RNAME,' NOT IN WATERIONS.LIB'
         IF(MASWRK)WRITE(IW,*)' '
         CALL ABRT
      END IF
      CALL OPNCRD(12,-IW)
C
  100 CONTINUE
      CALL RDCARD('WATERION',IEOF)
      IF(IEOF.EQ.1) THEN
         IF(MASWRK) WRITE(IW,*)
     *   'ERROR: END OF FILE READING WATERIONS.LIB'
         IF(MASWRK)WRITE(IW,*)' '
         CALL ABRT
      END IF
      WORD ='          '
      KSIZE = -10
      CALL GSTRNG(WORD,KSIZE)
C
      IF(WORD.EQ.'COORDINATE'.OR.
     *   WORD.EQ.'BOND      '.OR.
     *   WORD.EQ.'ANGLE     '.OR.
     *   WORD.EQ.'DIHROT    '.OR.
     *   WORD.EQ.'DIHBND    '.OR.
     *   WORD.EQ.'WAGGING   ') THEN
         CALL FNDGRP(12,' STOP   ',IEOF)
         GOTO 100
      END IF
C
      IF(WORD.EQ.'PARAMETERS') THEN
         IERR = 0
         IEOF = 0
         CALL RDCARD('RDPARAMT',IEOF)
         FFNAME = '          '
         LGSTR=-10
         CALL GSTRNG(FFNAME,LGSTR)
         ZMAS(IFFAT)  = RFIND('FFMASS  ',IERR)*1822.88850204D+00
         ONEMAS(IFFAT)= 1.0D+00/ZMAS(IFFAT)
         CHARG(IFFAT) = RFIND('FFCHARGE',IERR)
         POL(IFFAT)   = RFIND('FFPOL   ',IERR)*TOBOHR**3
         SIG(IFFAT)   = RFIND('FFSIGMA ',IERR)*TOBOHR
         EPS(IFFAT)   = RFIND('FFEPSILN',IERR)*TOHART
         SIG2(IFFAT)  = RFIND('FFSIGMA ',IERR)*TOBOHR
         EPS2(IFFAT)  = RFIND('FFEPSILN',IERR)*TOHART
C        - WATERIONS.LIB USES RMIN/2 INSTAED OF SIGMA
         SIG(IFFAT)  = SIG(IFFAT) *1.781797436280679D+00
         SIG2(IFFAT) = SIG2(IFFAT)*1.781797436280679D+00
         CALL FNDGRP(12,' STOP   ',IEOF)
         GOTO 100
      END IF
C
      IF(WORD.EQ.'QMMMREP   ') THEN
         IERR = 0
         IEOF = 0
         CALL RDCARD('RDQMMMRP',IEOF)
         FFNAME = '          '
         LGSTR=-10
         CALL GSTRNG(FFNAME,LGSTR)
         NTERMS = IFIND('NTERM   ',IERR)
         DO II = 1, NTERMS
            CLPR(II,IFFAT) = RFIND('CLPR    ',IERR)
            ZLPR(II,IFFAT) = RFIND('ZLPR    ',IERR)
         ENDDO
         NLPR(IFFAT) = 2
         CALL FNDGRP(12,' STOP   ',IEOF)
         GOTO 100
      END IF
C
      IF(WORD.EQ.'$END      ') THEN
         IF(MASWRK) CLOSE(12)
         RETURN
      END IF
C
      END
C*MODULE QUANPOB  *DECK RDHOH
!>
!> @brief    read water parameters
!>
!> @author   Hui Li
!>           - May 2011
!>
!> @details  read water parameters from WATERIONS.LIB
!>
      SUBROUTINE RDHOH(ATMNAM,ZMAS,ONEMAS,CHARG,POL,SIG,EPS,
     *                 SIG2,EPS2,BOND0,FCBOND,ANGL0,FCANGL,
     *                 IPAIR,KLIST,CLPR,ZLPR,NLPR,IFFAT,
     *                 L1213J,ITYPWAT,LSTRAT,DSTRAT)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      LOGICAL GOPARR,DSKWRK,MASWRK
C
      PARAMETER (TOANGS=0.52917724924D+00)
      PARAMETER (TOBOHR=1.0D+00/TOANGS)
      PARAMETER (DEGREE=57.2957795130823D+00)
      PARAMETER (TORAD=1.0D+00/DEGREE)
      PARAMETER (TOKCAL=627.509469D+00)
      PARAMETER (TOHART=1.0D+00/TOKCAL)
C
      CHARACTER*8 RNAME
      CHARACTER*10 WORD,ATMNAM,FFNAME,HEAD
      CHARACTER*256 QPFILE
C
      DIMENSION ATMNAM(*),ZMAS(*),ONEMAS(*),CHARG(*),POL(*),
     *          SIG(*),EPS(*),SIG2(*),EPS2(*),BOND0(*),FCBOND(*),
     *          ANGL0(*),FCANGL(*),IPAIR(2,*),KLIST(3,*),CLPR(4,*),
     *          ZLPR(4,*),NLPR(*),L1213J(2,*),LSTRAT(2,*),DSTRAT(*)
C
      COMMON /FFPARM/ NFFAT,NBOND,NANGL,NDIHR,NDIHB,NCMAP,NWAGG,
     *                N1213J,N14J,NLKQMM,IDOCHG,IDOPOL,IDOLJ,IDOCMAP
      COMMON /FFRATT/ RATOLC,RATOLV,SCALRAT,VIRRAT(3),IRATTLE,JRATTLE,
     *                NRATTLE,MXRATT,LFFOLDCORD,LFFLSTRAT,LFFDSTRAT,
     *                LFFVELSV,IRATQM
      COMMON /IOFILE/ IR,IW,IP,IJK,IJKT,IDAF,NAV,IODA(950)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
C
C     HUI LI, MAY 2011, LINCOLN
C
C     -- READ IN WATER PARAMETERS FROM WATERIONS.LIB --
C
      WRITE(RNAME,'(A5,I3)') ' $HOH', ITYPWAT
      IF(MASWRK) THEN
         CALL GENQPFILE('WATERIONS.LIB ',QPFILE,LENQP)
         OPEN(12,FILE=QPFILE(1:LENQP),STATUS='OLD')
      END IF
      IEOF=0
      CALL SEQREW(12)
      CALL FNDGRP(12,RNAME,IEOF)
      IF(IEOF.EQ.1) THEN
         IF(MASWRK) WRITE(IW,*)
     *   'ERROR: ',RNAME,' NOT IN WATERIONS.LIB'
         IF(MASWRK)WRITE(IW,*)' '
         CALL ABRT
      END IF
      CALL OPNCRD(12,-IW)
C
  100 CONTINUE
      CALL RDCARD('WATERION',IEOF)
      IF(IEOF.EQ.1) THEN
         IF(MASWRK) WRITE(IW,*)
     *   'ERROR: END OF FILE READING WATERIONS.LIB'
         IF(MASWRK)WRITE(IW,*)' '
         CALL ABRT
      END IF
      WORD ='          '
      KSIZE = -10
      CALL GSTRNG(WORD,KSIZE)
C
      IF(WORD.EQ.'COORDINATE'.OR.
     *   WORD.EQ.'DIHROT    '.OR.
     *   WORD.EQ.'DIHBND    '.OR.
     *   WORD.EQ.'WAGGING   ') THEN
         CALL FNDGRP(12,' STOP   ',IEOF)
         GOTO 100
      END IF
C
      IF(WORD.EQ.'PARAMETERS') THEN
         IERR = 0
         IEOF = 0
         DO K=0,ITYPWAT/100-1
            CALL RDCARD('RDPARAMT',IEOF)
            FFNAME = '          '
            LGSTR=-10
            CALL GSTRNG(FFNAME,LGSTR)
            ATMNAM(IFFAT+K) = FFNAME
            ZMAS(IFFAT+K)  = RFIND('FFMASS  ',IERR)*1822.88850204D+00
            ONEMAS(IFFAT+K)= 1.0D+00/ZMAS(IFFAT+K)
            CHARG(IFFAT+K) = RFIND('FFCHARGE',IERR)
            POL(IFFAT+K)   = RFIND('FFPOL   ',IERR)*TOBOHR**3
            SIG(IFFAT+K)   = RFIND('FFSIGMA ',IERR)*TOBOHR
            EPS(IFFAT+K)   = RFIND('FFEPSILN',IERR)*TOHART
            SIG2(IFFAT+K)  = RFIND('FFSIGMA ',IERR)*TOBOHR
            EPS2(IFFAT+K)  = RFIND('FFEPSILN',IERR)*TOHART
C           - WATERIONS.LIB USES RMIN/2 INSTAED OF SIGMA
            SIG(IFFAT+K)  = SIG(IFFAT+K) *1.781797436280679D+00
            SIG2(IFFAT+K) = SIG2(IFFAT+K)*1.781797436280679D+00
         ENDDO
         CALL FNDGRP(12,' STOP   ',IEOF)
         GOTO 100
      END IF
C
      IF(WORD.EQ.'QMMMREP   ') THEN
         IERR = 0
         IEOF = 0
         DO K=0,ITYPWAT/100-1
            CALL RDCARD('RDQMMMRP',IEOF)
            FFNAME = '          '
            LGSTR=-10
            CALL GSTRNG(FFNAME,LGSTR)
            NTERMS = IFIND('NTERM   ',IERR)
            DO II = 1, NTERMS
               CLPR(II,IFFAT+K) = RFIND('CLPR    ',IERR)
               ZLPR(II,IFFAT+K) = RFIND('ZLPR    ',IERR)
            ENDDO
            NLPR(IFFAT+K) = 2
         ENDDO
         CALL FNDGRP(12,' STOP   ',IEOF)
         GOTO 100
      END IF
C
      IF(WORD.EQ.'BOND      ') THEN
         IERR = 0
         IEOF = 0
C        - REAL BONDS FOR ENERGY -
         IF(ITYPWAT/100.EQ.3) NREAL = 2
         IF(ITYPWAT/100.EQ.5) NREAL = 4
         DO K=1,NREAL
            CALL RDCARD('RDBOND  ',IEOF)
            HEAD='          '
            LGSTR=-10
            CALL GSTRNG(HEAD,LGSTR)
            NBOND = NBOND + 1
            IPAIR(1,NBOND) = IFIND('IPAIR  ',IERR) + IFFAT - 1
            IPAIR(2,NBOND) = IFIND('IPAIR  ',IERR) + IFFAT - 1
            FCBOND(NBOND)  = RFIND('FCBOND  ',IERR)*TOANGS**2*TOHART
            BOND0(NBOND)   = RFIND('R0      ',IERR)*TOBOHR
            N1213J = N1213J + 1
            L1213J(1,N1213J) = IPAIR(1,NBOND)
            L1213J(2,N1213J) = IPAIR(2,NBOND)
            IF(IRATTLE.GT.0) THEN
               NRATTLE = NRATTLE + 1
               LSTRAT(1,NRATTLE)=IPAIR(1,NBOND)
               LSTRAT(2,NRATTLE)=IPAIR(2,NBOND)
               DSTRAT(NRATTLE)  =BOND0(NBOND)*BOND0(NBOND)
            END IF
         ENDDO
C        - FAKE BONDS FOR CONSTRAINT -
         IF(ITYPWAT/100.EQ.3) NFAKE = 1
         IF(ITYPWAT/100.EQ.5) NFAKE = 5
         DO K=1,NFAKE
            CALL RDCARD('RDBOND  ',IEOF)
            HEAD='          '
            LGSTR=-10
            CALL GSTRNG(HEAD,LGSTR)
            NBOND = NBOND + 1
            IPAIR(1,NBOND) = IFIND('IPAIR  ',IERR) + IFFAT - 1
            IPAIR(2,NBOND) = IFIND('IPAIR  ',IERR) + IFFAT - 1
            FCBOND(NBOND)  = RFIND('FCBOND  ',IERR)*TOANGS**2*TOHART
            BOND0(NBOND)   = RFIND('R0      ',IERR)*TOBOHR
            IF(IRATTLE.EQ.10.OR.IRATTLE.EQ.20) THEN
               NRATTLE = NRATTLE + 1
               LSTRAT(1,NRATTLE)=IPAIR(1,NBOND)
               LSTRAT(2,NRATTLE)=IPAIR(2,NBOND)
               DSTRAT(NRATTLE)  =BOND0(NBOND)*BOND0(NBOND)
            END IF
         ENDDO
         CALL FNDGRP(12,' STOP   ',IEOF)
         GOTO 100
      END IF
C
      IF(WORD.EQ.'ANGLE     ') THEN
         IERR = 0
         IEOF = 0
         IF(ITYPWAT/100.EQ.3) NANG = 1
         IF(ITYPWAT/100.EQ.5) NANG = 6
         DO K=1,NANG
            CALL RDCARD('RDANGL  ',IEOF)
            HEAD='          '
            LGSTR=-10
            CALL GSTRNG(HEAD,LGSTR)
            NANGL = NANGL + 1
            KLIST(1,NANGL) = IFIND('KLIST1  ',IERR) + IFFAT - 1
            KLIST(2,NANGL) = IFIND('KLIST2  ',IERR) + IFFAT - 1
            KLIST(3,NANGL) = IFIND('KLIST3  ',IERR) + IFFAT - 1
            FCANGL(NANGL)  = RFIND('FCANGL  ',IERR)*TOHART
            ANGL0(NANGL)   = RFIND('ANGL0   ',IERR)*TORAD
            N1213J = N1213J + 1
            L1213J(1,N1213J) = KLIST(1,NANGL)
            L1213J(2,N1213J) = KLIST(3,NANGL)
         ENDDO
         CALL FNDGRP(12,' STOP   ',IEOF)
         GOTO 100
      END IF
C
      IF(WORD.EQ.'$END      ') THEN
         IF(MASWRK) CLOSE(12)
         RETURN
      END IF
C
      END
C*MODULE QUANPOB  *DECK FFBOND
      SUBROUTINE FFBOND(KFFAT,IPAIR,ZANF,CORD,KBOND)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      LOGICAL GOPARR,DSKWRK,MASWRK
      DIMENSION RCOV(86),ZANF(*),CORD(3,*)
C
      PARAMETER (TOANGS=0.52917724924D+00)
      PARAMETER (TOBOHR=1.0D+00/TOANGS)
C
      DIMENSION IPAIR(2,*)
      COMMON /IOFILE/ IR,IW,IP,IJK,IJKT,IDAF,NAV,IODA(950)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
C
C      COVALENT RADII FROM J.EMSLEY, "THE ELEMENTS", 2ND EDITION, 1991
C      EXCEPT VAN DER WAALS RADII FOR HE,NE,AR,KR (SAME SOURCE),
C      AND GUESSES FOR NA,V,CR,RB,TC,PM,EU,YB,AT,RN
C
      DATA (RCOV(NUCZ),NUCZ=1,2)/0.30D+00,1.22D+00/
      DATA (RCOV(NUCZ),NUCZ=3,10)
     *  /1.23D+00,0.89D+00,0.88D+00,0.77D+00,
     *   0.70D+00,0.66D+00,0.58D+00,1.60D+00/
      DATA (RCOV(NUCZ),NUCZ=11,18)
     *  /1.66D+00,1.36D+00,1.25D+00,1.09D+00,
     *   1.10D+00,1.04D+00,1.02D+00,1.91D+00/
      DATA (RCOV(NUCZ),NUCZ=19,36)
     *  /2.03D+00,1.74D+00,
     *   1.44D+00,1.32D+00,1.22D+00,1.19D+00,1.17D+00,
     *   1.165D+00,1.16D+00,1.15D+00,1.17D+00,1.25D+00,
     *   1.25D+00,1.22D+00,1.21D+00,1.17D+00,1.14D+00,1.98D+00/
      DATA (RCOV(NUCZ),NUCZ=37,54)
     *  /2.22D+00,1.92D+00,
     *   1.62D+00,1.45D+00,1.34D+00,1.29D+00,1.27D+00,
     *   1.24D+00,1.25D+00,1.28D+00,1.34D+00,1.41D+00,
     *   1.50D+00,1.40D+00,1.41D+00,1.37D+00,1.33D+00,2.09D+00/
      DATA (RCOV(NUCZ),NUCZ=55,86)
     *  /2.35D+00,1.98D+00,
     *   1.69D+00,1.65D+00,1.65D+00,1.64D+00,1.65D+00,1.66D+00,1.65D+00,
     *   1.61D+00,1.59D+00,1.59D+00,1.58D+00,1.57D+00,1.56D+00,1.56D+00,
     *   1.56D+00,1.44D+00,1.34D+00,1.30D+00,1.28D+00,
     *   1.26D+00,1.26D+00,1.29D+00,1.34D+00,1.44D+00,
     *   1.55D+00,1.54D+00,1.52D+00,1.53D+00,1.50D+00,2.20D+00/
C
C     ADAPTED BY NANDUN THELLAMUREGE, JAN 2011, LINCOLN
C
      KBOND=0
      DO 130 I=1,KFFAT
         INUCZ = INT(ZANF(I))
         RADI = 1.6D+00
         IF(INUCZ.EQ.1)                 RADI =         RCOV(1)
         IF(INUCZ.GT.1.AND.INUCZ.LE.86) RADI = 1.2D+00*RCOV(INUCZ)
         DO 120 J=I+1,KFFAT
            JNUCZ = INT(ZANF(J))
            RADJ = 1.6D+00
            IF(JNUCZ.EQ.1)                 RADJ =         RCOV(1)
            IF(JNUCZ.GT.1.AND.JNUCZ.LE.86) RADJ = 1.2D+00*RCOV(JNUCZ)
            DIST = SQRT((CORD(1,I)*TOANGS-CORD(1,J)*TOANGS)**2
     *                + (CORD(2,I)*TOANGS-CORD(2,J)*TOANGS)**2
     *                + (CORD(3,I)*TOANGS-CORD(3,J)*TOANGS)**2)
            BOND = RADI + RADJ
            IF((INUCZ.EQ. 2.OR.INUCZ.EQ.10.OR.INUCZ.EQ.18.OR.
     *          INUCZ.EQ.36.OR.INUCZ.EQ.54).OR.
     *         (JNUCZ.EQ. 2.OR.JNUCZ.EQ.10.OR.JNUCZ.EQ.18.OR.
     *          JNUCZ.EQ.36.OR.JNUCZ.EQ.54)) BOND = -1.0D+00
            IF(DIST.LE.BOND) THEN
               KBOND = KBOND + 1
               IPAIR(1,KBOND) = MIN(I,J)
               IPAIR(2,KBOND) = MAX(I,J)
            END IF
  120    CONTINUE
  130 CONTINUE
C
C     -- SOMETIMES H ATOMS ARE TOO FAR AWAY FROM HEAVY ATOMS --
C     -- QM/MM LINK H ATOMS: TYPICALLY 1.5-2.5 A --
C        HUI LI, DEC 17, 2019
C
      KKBOND = KBOND
      DO 200 I=1,KFFAT
         INUCZ = INT(ZANF(I))
         IF(INUCZ.EQ.1) THEN
            IYES = 0
            DO IBOND=1, KKBOND
               IF(I.EQ.IPAIR(1,IBOND)) IYES = 1
               IF(I.EQ.IPAIR(2,IBOND)) IYES = 1
            ENDDO
            DO IIIDO = 1, 16
             IF(IYES.EQ.0) THEN
               R2CUT = ((1.0D+00 + 0.1D+00*IIIDO)*TOBOHR)**2
               DO J=1,KFFAT
                JNUCZ = INT(ZANF(J))
                IF(JNUCZ.GT.1) THEN
                  DIST2 = (CORD(1,I)-CORD(1,J))**2
     *                  + (CORD(2,I)-CORD(2,J))**2
     *                  + (CORD(3,I)-CORD(3,J))**2
                  IF(DIST2.LE.R2CUT.AND.DIST2.GE.0.10D+00) THEN
                     KBOND = KBOND + 1
                     IPAIR(1,KBOND) = MIN(I,J)
                     IPAIR(2,KBOND) = MAX(I,J)
                     IYES = 1
                  END IF
                END IF
               ENDDO
             END IF
            ENDDO
            IF(IYES.EQ.0) THEN
               IF(MASWRK)WRITE(IW,'(/1X,A,I8,A/)')
     *         'ERROR: H ATOM ',I,' HAS BONDING PROBLEM.'
               CALL ABRT
            END IF
         END IF
  200 CONTINUE
C
      RETURN
      END
C*MODULE QUANPOB  *DECK FFANGL
      SUBROUTINE FFANGL(IPAIR,KBOND,KLIST,KANGL)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      INTEGER P1,P2,P3,P4
C
      DIMENSION IPAIR(2,*),KLIST(3,*)
C
C     HUI LI, MAR 2011, LINCOLN
C
C     ---- GET ANGLE LIST ----
C
      KANGL=0
      DO IBOND = 1, KBOND
         P1 = IPAIR(1,IBOND)
         P2 = IPAIR(2,IBOND)
         DO JBOND = IBOND+1, KBOND
            P3 = IPAIR(1,JBOND)
            P4 = IPAIR(2,JBOND)
            IF(P3.EQ.P2) THEN
               KANGL = KANGL + 1
               KLIST(1,KANGL) = MIN(P1,P4)
               KLIST(2,KANGL) = P2
               KLIST(3,KANGL) = MAX(P1,P4)
            END IF
         ENDDO
      ENDDO
      DO IBOND = 1, KBOND
         P1 = IPAIR(1,IBOND)
         P2 = IPAIR(2,IBOND)
         DO JBOND = IBOND+1, KBOND
            P3 = IPAIR(1,JBOND)
            P4 = IPAIR(2,JBOND)
            IF(P3.EQ.P1) THEN
               KANGL = KANGL + 1
               KLIST(1,KANGL) = MIN(P2,P4)
               KLIST(2,KANGL) = P1
               KLIST(3,KANGL) = MAX(P2,P4)
            END IF
         ENDDO
      ENDDO
      DO IBOND = 1, KBOND
         P1 = IPAIR(1,IBOND)
         P2 = IPAIR(2,IBOND)
         DO JBOND = IBOND+1, KBOND
            P3 = IPAIR(1,JBOND)
            P4 = IPAIR(2,JBOND)
            IF(P2.EQ.P4) THEN
               KANGL = KANGL + 1
               KLIST(1,KANGL) = MIN(P1,P3)
               KLIST(2,KANGL) = P2
               KLIST(3,KANGL) = MAX(P1,P3)
            END IF
         ENDDO
      ENDDO
C
      RETURN
      END
C*MODULE QUANPOB  *DECK FFDIHR
      SUBROUTINE FFDIHR(KLIST,KANGL,LLIST,KDIHR)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      DIMENSION KLIST(3,*),LLIST(4,*)
C
C     HUI LI, APR 2011, LINCOLN
C
C     ---- GET DIHEDRAL ROTATION ANGLE LIST ----
C          USE ONLY THE ANGLE LIST
C
      KDIHR=0
      DO I=1, KANGL
         DO J=I+1,KANGL
            IF(KLIST(2,I).EQ.KLIST(1,J).AND.
     *         KLIST(3,I).EQ.KLIST(2,J).AND.
     *         KLIST(1,I).LT.KLIST(3,J)) THEN
                  KDIHR=KDIHR+1
                  LLIST(1,KDIHR)=KLIST(1,I)
                  LLIST(2,KDIHR)=KLIST(2,I)
                  LLIST(3,KDIHR)=KLIST(3,I)
                  LLIST(4,KDIHR)=KLIST(3,J)
            END IF
            IF(KLIST(2,I).EQ.KLIST(1,J).AND.
     *         KLIST(3,I).EQ.KLIST(2,J).AND.
     *         KLIST(1,I).GT.KLIST(3,J)) THEN
                  KDIHR=KDIHR+1
                  LLIST(1,KDIHR)=KLIST(3,J)
                  LLIST(2,KDIHR)=KLIST(3,I)
                  LLIST(3,KDIHR)=KLIST(2,I)
                  LLIST(4,KDIHR)=KLIST(1,I)
            END IF
            IF(KLIST(2,I).EQ.KLIST(3,J).AND.
     *         KLIST(3,I).EQ.KLIST(2,J).AND.
     *         KLIST(1,I).LT.KLIST(1,J)) THEN
                  KDIHR=KDIHR+1
                  LLIST(1,KDIHR)=KLIST(1,I)
                  LLIST(2,KDIHR)=KLIST(2,I)
                  LLIST(3,KDIHR)=KLIST(3,I)
                  LLIST(4,KDIHR)=KLIST(1,J)
            END IF
            IF(KLIST(2,I).EQ.KLIST(3,J).AND.
     *         KLIST(3,I).EQ.KLIST(2,J).AND.
     *         KLIST(1,I).GT.KLIST(1,J)) THEN
                  KDIHR=KDIHR+1
                  LLIST(1,KDIHR)=KLIST(1,J)
                  LLIST(2,KDIHR)=KLIST(3,I)
                  LLIST(3,KDIHR)=KLIST(2,I)
                  LLIST(4,KDIHR)=KLIST(1,I)
            END IF
            IF(KLIST(1,I).EQ.KLIST(2,J).AND.
     *         KLIST(2,I).EQ.KLIST(1,J).AND.
     *         KLIST(3,I).LT.KLIST(3,J)) THEN
                  KDIHR=KDIHR+1
                  LLIST(1,KDIHR)=KLIST(3,I)
                  LLIST(2,KDIHR)=KLIST(2,I)
                  LLIST(3,KDIHR)=KLIST(1,I)
                  LLIST(4,KDIHR)=KLIST(3,J)
            END IF
            IF(KLIST(1,I).EQ.KLIST(2,J).AND.
     *         KLIST(2,I).EQ.KLIST(1,J).AND.
     *         KLIST(3,I).GT.KLIST(3,J)) THEN
                  KDIHR=KDIHR+1
                  LLIST(1,KDIHR)=KLIST(3,J)
                  LLIST(2,KDIHR)=KLIST(1,I)
                  LLIST(3,KDIHR)=KLIST(2,I)
                  LLIST(4,KDIHR)=KLIST(3,I)
            END IF
            IF(KLIST(1,I).EQ.KLIST(2,J).AND.
     *         KLIST(2,I).EQ.KLIST(3,J).AND.
     *         KLIST(3,I).LT.KLIST(1,J)) THEN
                  KDIHR=KDIHR+1
                  LLIST(1,KDIHR)=KLIST(3,I)
                  LLIST(2,KDIHR)=KLIST(2,I)
                  LLIST(3,KDIHR)=KLIST(1,I)
                  LLIST(4,KDIHR)=KLIST(1,J)
            END IF
            IF(KLIST(1,I).EQ.KLIST(2,J).AND.
     *         KLIST(2,I).EQ.KLIST(3,J).AND.
     *         KLIST(3,I).GT.KLIST(1,J)) THEN
                  KDIHR=KDIHR+1
                  LLIST(1,KDIHR)=KLIST(1,J)
                  LLIST(2,KDIHR)=KLIST(1,I)
                  LLIST(3,KDIHR)=KLIST(2,I)
                  LLIST(4,KDIHR)=KLIST(3,I)
            END IF
         ENDDO
      ENDDO
C
      RETURN
      END
C*MODULE QUANPOB  *DECK FFDIHB
!>
!> @brief    determine dihedral bending cases
!>
!> @author   Hui Li
!>           - Mar 2011
!>
!> @details  check all possible cases but save those
!>           close to 0 or 180 degree
!>
      SUBROUTINE FFDIHB(CORD,IPAIR,KLIST,NLIST,AAA)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      PARAMETER (DEGREE=57.2957795130823D+00)
      PARAMETER (TORAD=1.0D+00/DEGREE)
      PARAMETER (PT5=0.5D+00)
      PARAMETER (ONE=1.0D+00)
C
      INTEGER P1, P2, P3, P4
C
      DIMENSION CORD(3,*),IPAIR(2,*),KLIST(3,*),NLIST(4,*)
C
      COMMON /FFPARM/ NFFAT,NBOND,NANGL,NDIHR,NDIHB,NCMAP,NWAGG,
     *                N1213J,N14J,NLKQMM,IDOCHG,IDOPOL,IDOLJ,IDOCMAP
C
C     HUI LI, MAR 2011, LINCOLN
C
C     ---- GET DIHEDRAL BENDING ANGLE LIST ----
C          ONLY THOSE CLOSE TO 0.0 OR 180.0 DEGREE
C          WILL BE SAVED.
C
      COSAAA = COS(AAA*TORAD)
C
      NDIHB=0
      DO I=1, NANGL
         DO J=1,NBOND
            IF((KLIST(2,I).EQ.IPAIR(1,J)).AND.
     *         (KLIST(1,I).NE.IPAIR(2,J)).AND.
     *         (KLIST(3,I).NE.IPAIR(2,J)))THEN
                  P1=KLIST(2,I)
                  P2=KLIST(1,I)
                  P3=KLIST(3,I)
                  P4=IPAIR(2,J)
                  GOTO 100
            END IF
            IF((KLIST(2,I).EQ.IPAIR(2,J)).AND.
     *         (KLIST(1,I).NE.IPAIR(1,J)).AND.
     *         (KLIST(3,I).NE.IPAIR(1,J)))THEN
                  P1=KLIST(2,I)
                  P2=KLIST(1,I)
                  P3=KLIST(3,I)
                  P4=IPAIR(1,J)
                  GOTO 100
            END IF
            GOTO 200
C
 100        CONTINUE
            X12=CORD(1,P1)-CORD(1,P2)
            Y12=CORD(2,P1)-CORD(2,P2)
            Z12=CORD(3,P1)-CORD(3,P2)
            X13=CORD(1,P1)-CORD(1,P3)
            Y13=CORD(2,P1)-CORD(2,P3)
            Z13=CORD(3,P1)-CORD(3,P3)
            X14=CORD(1,P1)-CORD(1,P4)
            Y14=CORD(2,P1)-CORD(2,P4)
            Z14=CORD(3,P1)-CORD(3,P4)
            X23=CORD(1,P2)-CORD(1,P3)
            Y23=CORD(2,P2)-CORD(2,P3)
            Z23=CORD(3,P2)-CORD(3,P3)
            X34=CORD(1,P3)-CORD(1,P4)
            Y34=CORD(2,P3)-CORD(2,P4)
            Z34=CORD(3,P3)-CORD(3,P4)
            X24=CORD(1,P2)-CORD(1,P4)
            Y24=CORD(2,P2)-CORD(2,P4)
            Z24=CORD(3,P2)-CORD(3,P4)
C
            R14=SQRT(X14*X14+Y14*Y14+Z14*Z14)
            R13=SQRT(X13*X13+Y13*Y13+Z13*Z13)
            R12=SQRT(X12*X12+Y12*Y12+Z12*Z12)
            R23=SQRT(X23*X23+Y23*Y23+Z23*Z23)
            R34=SQRT(X34*X34+Y34*Y34+Z34*Z34)
            R24=SQRT(X24*X24+Y24*Y24+Z24*Z24)
C
C           HERE WE CHECK THE THREE ANGLES: 213, 214, 314
C
            ONEBC  =ONE/(R12*R13)
            COSA   =(R12*R12 + R13*R13 - R23*R23)*ONEBC*PT5
            IF(COSA.GT. ONE) COSA = ONE
            IF(COSA.LT.-ONE) COSA =-ONE
            ALPHA213 =ACOS(COSA)
            ONEBC  =ONE/(R12*R14)
            COSA   =(R12*R12 + R14*R14 - R24*R24)*ONEBC*PT5
            IF(COSA.GT. ONE) COSA = ONE
            IF(COSA.LT.-ONE) COSA =-ONE
            ALPHA214 =ACOS(COSA)
            ONEBC  =ONE/(R13*R14)
            COSA   =(R13*R13 + R14*R14 - R34*R34)*ONEBC*PT5
            IF(COSA.GT. ONE) COSA = ONE
            IF(COSA.LT.-ONE) COSA =-ONE
            ALPHA314 =ACOS(COSA)
            IF(ALPHA213*DEGREE.GE.160.0D+00) GOTO 200
            IF(ALPHA214*DEGREE.GE.160.0D+00) GOTO 200
            IF(ALPHA314*DEGREE.GE.160.0D+00) GOTO 200
C
            COS123=(-(X12*X23)-(Y12*Y23)-(Z12*Z23))/(R12*R23)
            COS234=(-(X23*X34)-(Y23*Y34)-(Z23*Z34))/(R23*R34)
            SIN2123= 1.0D+00-COS123*COS123
            SIN2234= 1.0D+00-COS234*COS234
            SIN123 = SQRT(ABS(SIN2123))
            SIN234 = SQRT(ABS(SIN2234))
            IF(ABS(SIN123).LT.1.0D-06) GOTO 200
            IF(ABS(SIN234).LT.1.0D-06) GOTO 200
            ONESIN = 1.0D+00/(SIN123*SIN234)
C
            COSTOR = ONESIN*(COS123*COS234-
     *               ((+X12*X34+Y12*Y34+Z12*Z34)/(R12*R34)))
            IF(ABS(COSTOR).GT.COSAAA) THEN   !   5 OR 30 DEGREE
               NDIHB=NDIHB+1
               NLIST(1,NDIHB)=P1
               NLIST(2,NDIHB)=P2
               NLIST(3,NDIHB)=P3
               NLIST(4,NDIHB)=P4
            END IF
 200        CONTINUE
         ENDDO
      ENDDO
C
      RETURN
      END
C*MODULE QUANPOB  *DECK FFWAGG
!>
!> @brief    MMFF94 wagging terms
!>
!> @author   Hongbo Zhu and Hui Li
!>           - Oct 2012
!>
!> @details  determine wagging terms for MMFF
!>
      SUBROUTINE FFWAGG(MLIST,NFFAT,NBOND,NWAGG,IPAIR,CORD)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      LOGICAL GOPARR,DSKWRK,MASWRK
      INTEGER P1, P2, P3, P4
      PARAMETER (PI=3.14159265358979323846264338D+00)
C
      DIMENSION MLIST(4,*),IPAIR(2,*),CORD(3,*)
C
      COMMON /FFMAX / MXFFAT,MXBOND,MXANGL,MXDIHR,MXDIHB,MXCMAP,
     *                MXWAGG
      COMMON /IOFILE/ IR,IW,IP,IJK,IJKT,IDAF,NAV,IODA(950)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
C
C     HONGBO ZHU, HUI LI, OCT 29, 2012, LINCOLN
C
C     -- ATOM 4 IS THE CENTRAL (LINKING) ATOM
C
      CALL VICLR(MLIST,1,4*MXWAGG)
      IWAGG=0
      DO IFFAT=1,NFFAT
         P4 = IFFAT
         LCOUNT1=0
         DO IBOND=1,NBOND
            J1=IPAIR(1,IBOND)
            J2=IPAIR(2,IBOND)
            IF(IFFAT.EQ.J1.OR.IFFAT.EQ.J2) THEN
               LCOUNT1=LCOUNT1+1
               IF(LCOUNT1.EQ.1) THEN
                  P1 = J1+J2-IFFAT
               END IF
               IF(LCOUNT1.EQ.2) THEN
                  P2 = J1+J2-IFFAT
               END IF
               IF(LCOUNT1.EQ.3) THEN
                  P3 = J1+J2-IFFAT
               END IF
            END IF
         ENDDO
         IF(LCOUNT1.EQ.3)THEN
            X41    = CORD(1,P4)-CORD(1,P1)
            Y41    = CORD(2,P4)-CORD(2,P1)
            Z41    = CORD(3,P4)-CORD(3,P1)
            X42    = CORD(1,P4)-CORD(1,P2)
            Y42    = CORD(2,P4)-CORD(2,P2)
            Z42    = CORD(3,P4)-CORD(3,P2)
            X43    = CORD(1,P4)-CORD(1,P3)
            Y43    = CORD(2,P4)-CORD(2,P3)
            Z43    = CORD(3,P4)-CORD(3,P3)
            X23    = CORD(1,P2)-CORD(1,P3)
            Y23    = CORD(2,P2)-CORD(2,P3)
            Z23    = CORD(3,P2)-CORD(3,P3)
C
            R41SQ  = X41*X41 + Y41*Y41 + Z41*Z41
            R41    = SQRT(R41SQ)
            ONER41 = 1.0D+00/R41
            R42SQ  = X42*X42 +Y42*Y42 + Z42*Z42
            R42    = SQRT(R42SQ)
            R43SQ  = X43*X43 + Y43*Y43 + Z43*Z43
            R43    = SQRT(R43SQ)
            R23SQ  = X23*X23 + Y23*Y23 + Z23*Z23
            R23    = SQRT(R23SQ)
            P4243X = Y42*Z43 - Z42*Y43
            P4243Y = Z42*X43 - X42*Z43
            P4243Z = X42*Y43 - Y42*X43
C
C           --- ANGLE 243 CAN BE 90 TO 120 ---
            COS243 = (R42*R42 + R43*R43 - R23*R23)/(2.0D+00*R42*R43)
            SIN243 = SQRT(ABS(1.0D+00 - COS243*COS243))
            DUM    = 1.0D+00/(SIN243*R42*R43)
            AX     = P4243X*DUM
            AY     = P4243Y*DUM
            AZ     = P4243Z*DUM
C           --- WWW IS LIKELY -90 TO +90 ---
            SINW   = -(AX*X41 + AY*Y41 + AZ*Z41)*ONER41
            WWW    = ASIN(SINW)
            AL     = PI/2.0D+00
            IF(ABS(WWW).LT.AL)THEN
               IWAGG=IWAGG+1
               MLIST(1,IWAGG) = P1
               MLIST(2,IWAGG) = P2
               MLIST(3,IWAGG) = P3
               MLIST(4,IWAGG) = P4
               IWAGG=IWAGG+1
               MLIST(1,IWAGG) = P2
               MLIST(2,IWAGG) = P3
               MLIST(3,IWAGG) = P1
               MLIST(4,IWAGG) = P4
               IWAGG=IWAGG+1
               MLIST(1,IWAGG) = P3
               MLIST(2,IWAGG) = P1
               MLIST(3,IWAGG) = P2
               MLIST(4,IWAGG) = P4
            END IF
         END IF
      ENDDO
      NWAGG=IWAGG
      IF(NWAGG.GT.MXWAGG) THEN
         IF(MASWRK)WRITE(IW,*)
     *   'ERROR: TOO MANY NWAGG. INCREASE MXWAGG.'
         IF(MASWRK)WRITE(IW,*)' '
         CALL ABRT
      END IF
C
      RETURN
      END
C*MODULE QUANPOB  *DECK NONBOND
!>
!> @brief    cell-list
!>
!> @author   Hui Li
!>           - Mar 2011
!>
!> @details  use the fast-list method
!>
      SUBROUTINE NONBOND(ISTEP,CORD,CORDSV,CORDSV2,CORDSVQ,
     *                   NONLS1,NONLS2,
     *                   NONLSTQ,LSTCELL,
     *                   MVFASTS2,MVFASTS3,MVFASTS4,
     *                   MVFASTL2,MVFASTL3,MVFASTL4,
     *                   NONLSA,NONLSB,
     *                   NONLSPMA,NONLSPMB)
      use mx_limits, only: mxatm
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      LOGICAL GOPARR, DSKWRK, MASWRK
C
      PARAMETER (TOANGS=0.52917724924D+00)
      PARAMETER (TOBOHR=1.0D+00/TOANGS)
      PARAMETER (PT5=0.50D+00)
      PARAMETER (HUGE=1.0D+30*TOBOHR)
C
      DIMENSION CORD(3,*),CORDSV(3,*),NONLS1(2,*),NONLS2(2,*),
     *          NONLSTQ(*),
     *          LSTCELL(*),CORDSV2(3,*),CORDSVQ(3,*)
      DIMENSION MVFASTS2(*),MVFASTS3(*),MVFASTS4(*),
     *          MVFASTS5(300),MVFASTS6(100),MVFASTS7(100),
     *          MVFASTL2(*),MVFASTL3(*),MVFASTL4(*),
     *          MVFASTL5(300),MVFASTL6(100),MVFASTL7(100)
      DIMENSION NFASTS(8),NFASTL(8)
      DIMENSION NONLSA(2,*),NONLSB(2,*),
     *          NONLSPMA(2,*),NONLSPMB(2,*)
C
      COMMON /FFDFS / TIMDFS,QDION,AMION,TEFF,NDFS,NATMGAS,
     *                LFFDFSC,
     *                LFFDFSC0,LFFDFSA,LFFDFSN,LFFDFCOM,KDFS,LFFDFSCAV
      COMMON /FFFAST/ KLARGE(6),KSMALL(6),KQMMM(6)
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
      COMMON /FFMAX / MXFFAT,MXBOND,MXANGL,MXDIHR,MXDIHB,MXCMAP,
     *                MXWAGG
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
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /IOFILE/ IR,IW,IP,IJK,IJKT,IDAF,NAV,IODA(950)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
C
      INTEGER, PARAMETER :: K15 = SELECTED_INT_KIND(15)
      INTEGER, PARAMETER :: MAGIC1 = TRANSFER(201603081127_K15,1)
      INTEGER, PARAMETER :: MAGIC2 = TRANSFER(20170303193144_K15,1)
C
C     HUI LI, MAR 2011, LINCOLN
C     HUI LI AND FENGCHAO CUI, OCT 2011, LINCOLN,
C     ADD CELL-LIST AND FAST-LIST METHODS, AND LARGE-TO-SMALL SCHEME
C
C     -- FOR SMALL AND SPARSE SYSTEMS JUST DO ALL ATOMS --
C
      IF(NFFAT.LE. 128) THEN
         IF(MXLIST1.EQ.MAGIC1) GOTO 251
         NTODO = 0
         IPCOUNT  = ME - 1
         DO 100 IFFAT=1, NFFAT-1
            DO 110 JFFAT=IFFAT+1,NFFAT
               IF(GOPARR) THEN
                  IPCOUNT = IPCOUNT + 1
                  IF(MOD(IPCOUNT,NPROC).NE.0) GOTO 110
               END IF
               NTODO = NTODO + 1
               NONLS1(1,NTODO) = IFFAT
               NONLS1(2,NTODO) = JFFAT
 110        CONTINUE
 100     CONTINUE
         IF(NDFS.EQ.MAGIC2) THEN
          NONLSTQ(NFFAT) = NFFAT-NATPDB      !  BORROW IT
          DO IFFAT = NATPDB+1,NFFAT
            NONLSTQ(IFFAT-NATPDB) = IFFAT
          ENDDO
         END IF
         NTODOQ  = 0
         IF(NAT.GT.0) THEN
            DO IFFAT=1, NFFAT
              NTODOQ = NTODOQ + 1
              NONLSTQ(NTODOQ) = IFFAT
            ENDDO
         END IF
         MXLIST1 = MAGIC1
         RETURN
      END IF
C
C     -- DO NOT UPDATE NEIGHBOR LISTS IF SPHSOL OR FIXSOL IS USED --
C
      IF(ISTEP.EQ.0) THEN
         CALL VICLR(KLARGE,1,6)
         CALL VICLR(KSMALL,1,6)
         CALL VICLR(KQMMM,1,6)
      END IF
      IF(ISTEP.GT.0 .AND. (ISTEP.LT.100.AND.IADDWAT.EQ.2)) RETURN
      IF(ISTEP.GT.0 .AND. (ISPHSOL.GT.0.OR. IFIXSOL.GT.0)) RETURN
C
C     -- SOMETIMES IT IS SAFE TO SKIP CHECKING THE SMALL AND QMMM LISTS
      IF((ISTEP-KSMALL(5)).LT.KSMALL(6) .AND. NAT.LE.0) RETURN
      IF((ISTEP-KSMALL(5)).LT.KSMALL(6) .AND. NAT.GT.0) THEN
         IF((ISTEP-KQMMM(5)).LT.KQMMM(6)) THEN
            RETURN
         ELSE
            GOTO 299
         END IF
      END IF
C
      NTODO   = NTODOSV      ! NTODOSV  IS DIFFERENT AT DIFFERENT CPUS
      TSMALL2 = (0.20D+00*BUFWID1)**2
      TSMALL3 = (0.30D+00*BUFWID1)**2
      TSMALL4 = (0.40D+00*BUFWID1)**2
      TSMALL5 = (0.50D+00*BUFWID1)**2
      TSMALL6 = (0.60D+00*BUFWID1)**2
      TSMALL7 = (0.70D+00*BUFWID1)**2
      TSMALL8 = (0.80D+00*BUFWID1)**2
      RSMALL  = SWRB + BUFWID1
      R2SMAL  = RSMALL**2
C
      NTODO2  = NTODO2SV     ! NTODO2SV IS DIFFERENT AT DIFFERENT CPUS
      TLARGE2 = (0.20D+00*(BUFWID2-BUFWID1))**2
      TLARGE3 = (0.30D+00*(BUFWID2-BUFWID1))**2
      TLARGE4 = (0.40D+00*(BUFWID2-BUFWID1))**2
      TLARGE5 = (0.50D+00*(BUFWID2-BUFWID1))**2
      TLARGE6 = (0.60D+00*(BUFWID2-BUFWID1))**2
      TLARGE7 = (0.70D+00*(BUFWID2-BUFWID1))**2
      TLARGE8 = (0.80D+00*(BUFWID2-BUFWID1))**2
      RLARGE  = SWRB + BUFWID2
      R2LARG  = RLARGE**2
C
      CALL VICLR(NFASTS,1,8)
      DO IFFAT=1, NFFAT
         DX = CORD(1,IFFAT)-CORDSV(1,IFFAT)
         DY = CORD(2,IFFAT)-CORDSV(2,IFFAT)
         DZ = CORD(3,IFFAT)-CORDSV(3,IFFAT)
         DIS= DX*DX+DY*DY+DZ*DZ
         IF     (DIS.GT.TSMALL8.AND.NFASTS(8).LT.100)THEN
            NFASTS(8)           = NFASTS(8) + 1
         ELSE IF(DIS.GT.TSMALL7.AND.NFASTS(7).LT.100)THEN
            NFASTS(7)           = NFASTS(7) + 1
            MVFASTS7(NFASTS(7)) = IFFAT
         ELSE IF(DIS.GT.TSMALL6.AND.NFASTS(6).LT.100)THEN
            NFASTS(6)           = NFASTS(6) + 1
            MVFASTS6(NFASTS(6)) = IFFAT
         ELSE IF(DIS.GT.TSMALL5.AND.NFASTS(5).LT.300)THEN
            NFASTS(5)           = NFASTS(5) + 1
            MVFASTS5(NFASTS(5)) = IFFAT
         ELSE IF(DIS.GT.TSMALL4.AND.NFASTS(4).LT.(NFFAT/4 +10))THEN
            NFASTS(4)           = NFASTS(4) + 1
            MVFASTS4(NFASTS(4)) = IFFAT
         ELSE IF(DIS.GT.TSMALL3.AND.NFASTS(3).LT.(NFFAT/2 +10))THEN
            NFASTS(3)           = NFASTS(3) + 1
            MVFASTS3(NFASTS(3)) = IFFAT
         ELSE IF(DIS.GT.TSMALL2.AND.NFASTS(2).LT.(NFFAT/2 +10))THEN
            NFASTS(2)           = NFASTS(2) + 1
            MVFASTS2(NFASTS(2)) = IFFAT
         END IF
      ENDDO
      ISMALL = 0
      IF(NFASTS(2).GE.(NFFAT/2 +10).OR.
     *   NFASTS(3).GE.(NFFAT/2 +10).OR.
     *   NFASTS(4).GE.(NFFAT/4 +10).OR.
     *   NFASTS(5).GE.MXCHECK.OR.
     *   NFASTS(6).GE.100.OR.
     *   NFASTS(7).GE.100.OR.
     *   NFASTS(8).GE.1) THEN
         ISMALL = 1
      ELSE IF((NFASTS(5)+NFASTS(6)+NFASTS(7)).EQ.0)THEN
         ISMALL = 0
      ELSE
         DO IBIN=5,7
         DO IFAST=1,NFASTS(IBIN)
            IF(IBIN.EQ.5) I=MVFASTS5(IFAST)
            IF(IBIN.EQ.6) I=MVFASTS6(IFAST)
            IF(IBIN.EQ.7) I=MVFASTS7(IFAST)
            XISV=CORDSV(1,I)
            YISV=CORDSV(2,I)
            ZISV=CORDSV(3,I)
            XI  =CORD(1,I)
            YI  =CORD(2,I)
            ZI  =CORD(3,I)
            DO 51 JBIN=2,7
               IF(JBIN.EQ.2.AND.IBIN.LT.7) GOTO 51
               IF(JBIN.EQ.3.AND.IBIN.LT.6) GOTO 51
               IF(JBIN.EQ.5.AND.IBIN.GT.5) GOTO 51
               IF(JBIN.EQ.6.AND.IBIN.GT.6) GOTO 51
               DO 52 JFAST=1,NFASTS(JBIN)
               IF(JBIN.EQ.2) J=MVFASTS2(JFAST)
               IF(JBIN.EQ.3) J=MVFASTS3(JFAST)
               IF(JBIN.EQ.4) J=MVFASTS4(JFAST)
               IF(JBIN.EQ.5) J=MVFASTS5(JFAST)
               IF(JBIN.EQ.6) J=MVFASTS6(JFAST)
               IF(JBIN.EQ.7) J=MVFASTS7(JFAST)
               IF(IBIN.EQ.JBIN .AND. J.LE.I) GOTO 52
               XJSV=CORDSV(1,J)
               YJSV=CORDSV(2,J)
               ZJSV=CORDSV(3,J)
               XJ  =CORD(1,J)
               YJ  =CORD(2,J)
               ZJ  =CORD(3,J)
               DX  =XI-XJ
               DY  =YI-YJ
               DZ  =ZI-ZJ
               PBCX=XBOX*ANINT(DX*ONEXBOX)
               PBCY=YBOX*ANINT(DY*ONEYBOX)
               PBCZ=ZBOX*ANINT(DZ*ONEZBOX)
               DX  =DX - PBCX
               DY  =DY - PBCY
               DZ  =DZ - PBCZ
               IF(ABS(DX).GT.SWRB) GOTO 52
               IF(ABS(DY).GT.SWRB) GOTO 52
               IF(ABS(DZ).GT.SWRB) GOTO 52
               DIS=DX*DX+DY*DY+DZ*DZ
               IF(DIS.LE.SWRB2) THEN
                  DX  =XISV-XJSV
                  DY  =YISV-YJSV
                  DZ  =ZISV-ZJSV
                  PBCX=XBOX*ANINT(DX*ONEXBOX)
                  PBCY=YBOX*ANINT(DY*ONEYBOX)
                  PBCZ=ZBOX*ANINT(DZ*ONEZBOX)
                  DX  =DX - PBCX
                  DY  =DY - PBCY
                  DZ  =DZ - PBCZ
                  DISSV=DX*DX+DY*DY+DZ*DZ
                  IF(DISSV.GT.R2SMAL) THEN
                     IF(MASWRK) THEN
                        NTODO           = NTODO + 1
                        NONLS1(1,NTODO) = I
                        NONLS1(2,NTODO) = J
                        IF(NTODO.GE.MXLIST1-100) THEN
                           WRITE(IW,*)
     *                     'ERROR: SMALL NEIGHBOR LIST IS FULL.',
     *                     ' INCREASE MXLIST1.'
                           WRITE(IW,*)' '
                           CALL ABRT
                        END IF
                     END IF
                  END IF
               END IF
  52           CONTINUE
  51        CONTINUE
         ENDDO
         ENDDO
      END IF
C
C     -- NO NEED TO UPDATE THE LARGE LIST IF THE SMALL LIST
C        DOES NOT NEED TO UPDATE
      IF(ISMALL.EQ.0.AND.ISTEP.GT.0.AND.NAT.LE.0) RETURN
      IF(ISMALL.EQ.0.AND.ISTEP.GT.0.AND.NAT.GT.0) GOTO 299
C
      IF(ISMALL.GT.0 .AND. BUFWID2.EQ.BUFWID1) THEN
         ILARGE = 1
         CALL ICOPY(8,NFASTS,1,NFASTL,1)
         GOTO 69
      END IF
C
C     -- SOMETIMES IT IS SAFE TO SKIP CHECKING THE LARGE LIST
C        AND GO DIRECTLY TO UPDATE THE SMALL LIST --
      IF((ISTEP-KLARGE(5)).LT.KLARGE(6)) GOTO 230
C
      CALL VICLR(NFASTL,1,8)
      DO IFFAT=1, NFFAT
         DX = CORD(1,IFFAT)-CORDSV2(1,IFFAT)
         DY = CORD(2,IFFAT)-CORDSV2(2,IFFAT)
         DZ = CORD(3,IFFAT)-CORDSV2(3,IFFAT)
         DIS= DX*DX+DY*DY+DZ*DZ
         IF     (DIS.GT.TLARGE8.AND.NFASTL(8).LT.100)THEN
            NFASTL(8)           = NFASTL(8) + 1
         ELSE IF(DIS.GT.TLARGE7.AND.NFASTL(7).LT.100)THEN
            NFASTL(7)           = NFASTL(7) + 1
            MVFASTL7(NFASTL(7)) = IFFAT
         ELSE IF(DIS.GT.TLARGE6.AND.NFASTL(6).LT.100)THEN
            NFASTL(6)           = NFASTL(6) + 1
            MVFASTL6(NFASTL(6)) = IFFAT
         ELSE IF(DIS.GT.TLARGE5.AND.NFASTL(5).LT.300)THEN
            NFASTL(5)           = NFASTL(5) + 1
            MVFASTL5(NFASTL(5)) = IFFAT
         ELSE IF(DIS.GT.TLARGE4.AND.NFASTL(4).LT.(NFFAT/4 -1))THEN
            NFASTL(4)           = NFASTL(4) + 1
            MVFASTL4(NFASTL(4)) = IFFAT
         ELSE IF(DIS.GT.TLARGE3.AND.NFASTL(3).LT.(NFFAT/2 -1))THEN
            NFASTL(3)           = NFASTL(3) + 1
            MVFASTL3(NFASTL(3)) = IFFAT
         ELSE IF(DIS.GT.TLARGE2.AND.NFASTL(2).LT.(NFFAT/2 -1))THEN
            NFASTL(2)           = NFASTL(2) + 1
            MVFASTL2(NFASTL(2)) = IFFAT
         END IF
      ENDDO
      ILARGE = 0
      IF(NFASTL(2).GE.(NFFAT/2-1).OR.
     *   NFASTL(3).GE.(NFFAT/2-1).OR.
     *   NFASTL(4).GE.(NFFAT/4-1).OR.
     *   NFASTL(5).GE.MXCHECK.OR.
     *   NFASTL(6).GE.100.OR.
     *   NFASTL(7).GE.100.OR.
     *   NFASTL(8).GE.1) THEN
         ILARGE = 1
      ELSE IF((NFASTL(5)+NFASTL(6)+NFASTL(7)).EQ.0)THEN
         ILARGE = 0
      ELSE
         DO IBIN=5,7
         DO IFAST=1,NFASTL(IBIN)
            IF(IBIN.EQ.5) I=MVFASTL5(IFAST)
            IF(IBIN.EQ.6) I=MVFASTL6(IFAST)
            IF(IBIN.EQ.7) I=MVFASTL7(IFAST)
            XISV=CORDSV2(1,I)
            YISV=CORDSV2(2,I)
            ZISV=CORDSV2(3,I)
            XI  =CORD(1,I)
            YI  =CORD(2,I)
            ZI  =CORD(3,I)
            DO 61 JBIN=2,7
               IF(JBIN.EQ.2.AND.IBIN.LT.7) GOTO 61
               IF(JBIN.EQ.3.AND.IBIN.LT.6) GOTO 61
               IF(JBIN.EQ.5.AND.IBIN.GT.5) GOTO 61
               IF(JBIN.EQ.6.AND.IBIN.GT.6) GOTO 61
               DO 62 JFAST=1,NFASTL(JBIN)
               IF(JBIN.EQ.2) J=MVFASTL2(JFAST)
               IF(JBIN.EQ.3) J=MVFASTL3(JFAST)
               IF(JBIN.EQ.4) J=MVFASTL4(JFAST)
               IF(JBIN.EQ.5) J=MVFASTL5(JFAST)
               IF(JBIN.EQ.6) J=MVFASTL6(JFAST)
               IF(JBIN.EQ.7) J=MVFASTL7(JFAST)
               IF(IBIN.EQ.JBIN .AND. J.LE.I) GOTO 62
               XJSV=CORDSV2(1,J)
               YJSV=CORDSV2(2,J)
               ZJSV=CORDSV2(3,J)
               XJ  =CORD(1,J)
               YJ  =CORD(2,J)
               ZJ  =CORD(3,J)
               DX  =XI-XJ
               DY  =YI-YJ
               DZ  =ZI-ZJ
               PBCX=XBOX*ANINT(DX*ONEXBOX)
               PBCY=YBOX*ANINT(DY*ONEYBOX)
               PBCZ=ZBOX*ANINT(DZ*ONEZBOX)
               DX  =DX - PBCX
               DY  =DY - PBCY
               DZ  =DZ - PBCZ
               IF(ABS(DX).GT.RSMALL) GOTO 62
               IF(ABS(DY).GT.RSMALL) GOTO 62
               IF(ABS(DZ).GT.RSMALL) GOTO 62
               DIS=DX*DX+DY*DY+DZ*DZ
               IF(DIS.LT.R2SMAL) THEN
                  DX  =XISV-XJSV
                  DY  =YISV-YJSV
                  DZ  =ZISV-ZJSV
                  PBCX=XBOX*ANINT(DX*ONEXBOX)
                  PBCY=YBOX*ANINT(DY*ONEYBOX)
                  PBCZ=ZBOX*ANINT(DZ*ONEZBOX)
                  DX  =DX - PBCX
                  DY  =DY - PBCY
                  DZ  =DZ - PBCZ
                  DISSV=DX*DX+DY*DY+DZ*DZ
                  IF(DISSV.GT.R2LARG) THEN
                     IF(MASWRK) THEN
                        NTODO2           = NTODO2 + 1
                        NONLS2(1,NTODO2) = I
                        NONLS2(2,NTODO2) = J
                        IF(NTODO2.GE.MXLIST2-100) THEN
                           WRITE(IW,*)
     *                     'ERROR: LARGE NEIGHBOR LIST IS FULL.',
     *                     ' INCREASE MXLIST2 (OR MXFFAT).'
                           WRITE(IW,*)' '
                           CALL ABRT
                        END IF
                     END IF
                  END IF
               END IF
  62           CONTINUE
  61        CONTINUE
         ENDDO
         ENDDO
      END IF
  69  CONTINUE
C
C     -- UPDATE ONLY THE SMALL NEIGHBOR LIST --
      IF(ILARGE.EQ.0 .AND. ISTEP.GT.0) GOTO 230
C
C
C     -- UPDATE BOTH LARGE AND SMALL NEIGHBOR LISTS --
C
      IF(MASWRK.AND.ISTEP.EQ.0) WRITE(IW,*)' '
      IF(ISTEP.EQ.0) CALL TIMIT(1)
      IF(MASWRK.AND.ISTEP.EQ.0) THEN
         WRITE(IW,*)' '
         WRITE(IW,'(1X,A,I10)')
     *           'UPDATING LARGE NEIGHBOR LIST AT STEP ',ISTEP
         WRITE(IW,*)' '
      END IF
      IF(MASWRK.AND.ISTEP.LE.10000) WRITE(IP,'(1X,A,I10)')
     *           'UPDATING LARGE NEIGHBOR LIST AT STEP ',ISTEP
      KLARGE(1) = KLARGE(2)
      KLARGE(2) = KLARGE(3)
      KLARGE(3) = KLARGE(4)
      KLARGE(4) = KLARGE(5)
      KLARGE(5) = ISTEP
      KLARGE(6) = (KLARGE(5) - KLARGE(1))/4
      KLARGE(6) = KLARGE(6) - NINT(SQRT(DBLE(KLARGE(6))))
      IF(MASWRK.AND.ISTEP.GT.0.AND.ISTEP.LE.10000)THEN
         WRITE(IP,'(A,I10,6I9)')' LARGE FAST-LIST:',
     *                     (IBIN,IBIN=2,8)
         WRITE(IP,'(A,I10,6I9)')' NUMBER OF ATOMS:',
     *                     (NFASTL(IBIN),IBIN=2,8)
      END IF
      CALL DCOPY(3*NFFAT,CORD,1,CORDSV2,1)
C
C     -- CREATE A CELL LIST --
      IF(XBOX.EQ.HUGE.OR.YBOX.EQ.HUGE.OR.ZBOX.EQ.HUGE) THEN
         NCELX    = 1
         NCELY    = 1
         NCELZ    = 1
         MXXX     = MXFFAT
      ELSE
         NCELX    = MAX(1,INT(XBOX/RLARGE))
         NCELY    = MAX(1,INT(YBOX/RLARGE))
         NCELZ    = MAX(1,INT(ZBOX/RLARGE))
         MXXX     = INT(XBOX*YBOX*ZBOX*0.02D+00)
      END IF
      MXATCEL  = MIN(MXFFAT,2*MXXX/(NCELX*NCELY*NCELZ))
      CALL VICLR(LSTCELL,1,MXATCEL*NCELZ*NCELY*NCELX)
      DO IFFAT=1, NFFAT
         XGRP  = CORD(1,IFFAT)
         YGRP  = CORD(2,IFFAT)
         ZGRP  = CORD(3,IFFAT)
         PBCX  = XBOX * ANINT(XGRP*ONEXBOX)
         PBCY  = YBOX * ANINT(YGRP*ONEYBOX)
         PBCZ  = ZBOX * ANINT(ZGRP*ONEZBOX)
         XGRP  = XGRP - PBCX + PT5*XBOX
         YGRP  = YGRP - PBCY + PT5*YBOX
         ZGRP  = ZGRP - PBCZ + PT5*ZBOX
         IIX   = INT(NCELX*XGRP*ONEXBOX) + 1
         IIY   = INT(NCELY*YGRP*ONEYBOX) + 1
         IIZ   = INT(NCELZ*ZGRP*ONEZBOX) + 1
         INUM  = ((IIX-1)*NCELY+IIY-1)*NCELZ+IIZ
         LSTCELL(INUM*MXATCEL)=LSTCELL(INUM*MXATCEL)+1
         IKK   = LSTCELL(INUM*MXATCEL)
         LSTCELL((INUM-1)*MXATCEL+IKK) =IFFAT
      ENDDO
C
      KFULL    = 0
      NTODO2   = 0    !  NTODO2 IS DIFFERENT FOR DIFFERENT CPUS.
      IPCOUNT  = ME - 1
      DO 201 ICELX = 1, NCELX
      DO 202 ICELY = 1, NCELY
      DO 203 ICELZ = 1, NCELZ
         INUM= ((ICELX-1)*NCELY+ICELY-1)*NCELZ+ICELZ
         DO 211 JCELX = 1, NCELX
         DO 212 JCELY = 1, NCELY
         DO 213 JCELZ = 1, NCELZ
            JNUM= ((JCELX-1)*NCELY+JCELY-1)*NCELZ+JCELZ
            IF(JNUM.LT.INUM)                   GOTO 213
            IDIFX = ABS(JCELX-ICELX)
            IF(JCELX.EQ.NCELX.AND.ICELX.EQ.1) IDIFX = 1
            IF(ICELX.EQ.NCELX.AND.JCELX.EQ.1) IDIFX = 1
            IF(JCELX.EQ.ICELX               ) IDIFX = 0
            IF(IDIFX.GE.2)                     GOTO 213
            IDIFY = ABS(JCELY-ICELY)
            IF(JCELY.EQ.NCELY.AND.ICELY.EQ.1) IDIFY = 1
            IF(ICELY.EQ.NCELY.AND.JCELY.EQ.1) IDIFY = 1
            IF(JCELY.EQ.ICELY               ) IDIFY = 0
            IF(IDIFY.GE.2)                     GOTO 213
            IDIFZ = ABS(JCELZ-ICELZ)
            IF(JCELZ.EQ.NCELZ.AND.ICELZ.EQ.1) IDIFZ = 1
            IF(ICELZ.EQ.NCELZ.AND.JCELZ.EQ.1) IDIFZ = 1
            IF(JCELZ.EQ.ICELZ               ) IDIFZ = 0
            IF(IDIFZ.GE.2)                     GOTO 213
            DO 200 IKK = 1, LSTCELL(INUM*MXATCEL)
               IFFAT = LSTCELL((INUM-1)*MXATCEL+IKK)
               IF(IFFAT.EQ.0) GOTO 200
C              -- HUI LI: THE PARALLEL DISTRIBUTION MUST BE HERE
C                 TO ACHIEVE EVEN-LOADING AND SCALABILITY:
C                    N_TASK = N_ATOM * (N_CELL + 1)/2
               IF(GOPARR) THEN
                  IPCOUNT = IPCOUNT + 1
                  IF(MOD(IPCOUNT,NPROC).NE.0) GOTO 200
               END IF
               DO 210 JKK = 1, LSTCELL(JNUM*MXATCEL)
                  JFFAT = LSTCELL((JNUM-1)*MXATCEL+JKK)
                  IF(JFFAT.EQ.0) GOTO 210
                  IF(INUM.EQ.JNUM.AND.JFFAT.LE.IFFAT) GOTO 210
                  XGRP=CORD(1,IFFAT)-CORD(1,JFFAT)
                  YGRP=CORD(2,IFFAT)-CORD(2,JFFAT)
                  ZGRP=CORD(3,IFFAT)-CORD(3,JFFAT)
                  PBCX=XBOX * ANINT(XGRP*ONEXBOX)
                  PBCY=YBOX * ANINT(YGRP*ONEYBOX)
                  PBCZ=ZBOX * ANINT(ZGRP*ONEZBOX)
                  XGRP=XGRP - PBCX
                  YGRP=YGRP - PBCY
                  ZGRP=ZGRP - PBCZ
                  IF(ABS(XGRP).GT.RLARGE) GOTO 210
                  IF(ABS(YGRP).GT.RLARGE) GOTO 210
                  IF(ABS(ZGRP).GT.RLARGE) GOTO 210
                  R2    = XGRP*XGRP+YGRP*YGRP+ZGRP*ZGRP
                  IF(R2.GT.R2LARG) GOTO 210
                  NTODO2 = NTODO2 + 1
                  NONLS2(1,NTODO2) = IFFAT
                  NONLS2(2,NTODO2) = JFFAT
 210           CONTINUE
               IF(NTODO2.GT.MXLIST2-MXATCEL) THEN
                  KFULL = 1
                  GOTO 220
               END IF
 200        CONTINUE
 213     CONTINUE
 212     CONTINUE
 211     CONTINUE
 203  CONTINUE
 202  CONTINUE
 201  CONTINUE
C
 220  CONTINUE
      NTODO2SV = NTODO2
      IF(GOPARR) CALL DDI_GSUMI(2408,KFULL,1)
      IF(MASWRK.AND.KFULL.GT.0) THEN
         WRITE(IW,*)'ERROR: LARGE NEIGHBOR LIST IS FULL (MXATCEL).',
     *              ' INCREASE MXLIST2 (OR MXFFAT).'
         WRITE(IW,*)' '
         CALL ABRT
      END IF
C
C     -- UPDATE ION-GAS POLARIZATION LIST --
C
      IF(NDFS.EQ.MAGIC2) THEN
      NONLSTQ(NFFAT) = 0        !  BORROW IT
      DO IFFAT = NATPDB+1,NFFAT
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
            IF(R2.LT.R2LARG) IYES=1
         ENDDO
         IF(IYES.EQ.1) THEN
            NONLSTQ(NFFAT) = NONLSTQ(NFFAT) + 1
            NONLSTQ(NONLSTQ(NFFAT)) = IFFAT
         END IF
      ENDDO
      END IF
C
 230  CONTINUE
C
C     -- UPDATE THE SMALL MM-MM NEIGHBOR LIST --
C
      IF(BUFWID2.GT.BUFWID1) THEN
         IF(ISTEP.EQ.0) CALL TIMIT(1)
         IF(MASWRK.AND.ISTEP.EQ.0) THEN
            WRITE(IW,*)' '
            WRITE(IW,'(1X,A,I10)')
     *           'UPDATING SMALL NEIGHBOR LIST AT STEP ',ISTEP
            WRITE(IW,*)' '
         END IF
         IF(MASWRK.AND.ISTEP.LE.10000) WRITE(IP,'(1X,A,I10)')
     *           'UPDATING SMALL NEIGHBOR LIST AT STEP ',ISTEP
         KSMALL(1) = KSMALL(2)
         KSMALL(2) = KSMALL(3)
         KSMALL(3) = KSMALL(4)
         KSMALL(4) = KSMALL(5)
         KSMALL(5) = ISTEP
         KSMALL(6) = (KSMALL(5) - KSMALL(1))/4
         KSMALL(6) = KSMALL(6) - NINT(SQRT(DBLE(KSMALL(6))))
         IF(MASWRK.AND.ISTEP.GT.0.AND.ISTEP.LE.10000)THEN
            WRITE(IP,'(A,I10,6I9)')' SMALL FAST-LIST:',
     *                        (IBIN,IBIN=2,8)
            WRITE(IP,'(A,I10,6I9)')' NUMBER OF ATOMS:',
     *                        (NFASTS(IBIN),IBIN=2,8)
         END IF
      END IF
      CALL DCOPY(3*NFFAT,CORD,1,CORDSV,1)
C
C     - A SPECIAL CASE IS BUFWID2 = BUFWID1 -
C       NOTE IN THIS CASE, LFFNONLS2 = LFFNONLS1
C       SO SIMPLY DEFINING NTODO=NTODO2 IS ENOUGH
      IF(BUFWID2.EQ.BUFWID1) THEN
         KFULL = 0
         NTODO = NTODO2
         GOTO 250
      END IF
C
      KFULL   = 0
      NTODO   = 0    !  NTODO IS DIFFERENT FOR DIFFERENT PROCESSOR.
      DO 240 III=1,NTODO2    !  NTODO2 IS DIFFERENT FOR DIFFERENT PROCESSOR.
         I=NONLS2(1,III)
         J=NONLS2(2,III)
         XGRP=CORD(1,I)-CORD(1,J)
         YGRP=CORD(2,I)-CORD(2,J)
         ZGRP=CORD(3,I)-CORD(3,J)
         PBCX=XBOX * ANINT(XGRP*ONEXBOX)
         PBCY=YBOX * ANINT(YGRP*ONEYBOX)
         PBCZ=ZBOX * ANINT(ZGRP*ONEZBOX)
         XGRP=XGRP - PBCX
         YGRP=YGRP - PBCY
         ZGRP=ZGRP - PBCZ
         IF(ABS(XGRP).GT.RSMALL) GOTO 240
         IF(ABS(YGRP).GT.RSMALL) GOTO 240
         IF(ABS(ZGRP).GT.RSMALL) GOTO 240
         R2    = XGRP*XGRP+YGRP*YGRP+ZGRP*ZGRP
         IF(R2.GT.R2SMAL) GOTO 240
         NTODO = NTODO + 1
         NONLS1(1,NTODO) = I
         NONLS1(2,NTODO) = J
         IF(NTODO.GT.MXLIST1-NFFAT)THEN
            KFULL = 1
            GOTO 250
         END IF
 240  CONTINUE
C
 250  CONTINUE
      NTODOSV = NTODO
      IF(GOPARR) CALL DDI_GSUMI(2408,KFULL,1)
      IF(MASWRK.AND.KFULL.GT.0) THEN
         WRITE(IW,*)'ERROR: SMALL NEIGHBOR LIST IS FULL.',
     *              ' INCREASE MXLIST1.'
         WRITE(IW,*)' '
         CALL ABRT
      END IF
 251  CONTINUE
C
C     - IDENTIFY PAIRS ASSOCIATED WITH SOL FREE ENERGY CALCULATION -
C       FOR IFEPTYP=1, MUST USE KFREEAB, CANNOT USE NFIXMM
      NTODOA = 0
      NTODOB = 0
      IF(IFEPTYP.GT.0) THEN
         DO I=1,NTODO
            I1    = NONLS1(1,I)
            I2    = NONLS1(2,I)
            I1YES = 0
            I2YES = 0
            DO III=1,KFREEAB(201)
               IA=KFREEAB(III)
               IF(I1.EQ.IA) I1YES = 1
               IF(I2.EQ.IA) I2YES = 1
            ENDDO
            IF((I1YES+I2YES).EQ.1) THEN
               NTODOA           = NTODOA + 1
               NONLSA(1,NTODOA) = (I1+I2)  - I1*I1YES - I2*I2YES
               NONLSA(2,NTODOA) =            I1*I1YES + I2*I2YES
            END IF
         ENDDO
         NTODOB = NTODOA
         CALL ICOPY(2*NTODOB,NONLSA,1,NONLSB,1)
      END IF
C
C     - IDENTIFY PAIRS ASSOCIATED WITH PMF FREE ENERGY CALCULATION -
      NTODOPMA = 0
      NTODOPMB = 0
      IF(IFEPTYP.EQ.2) THEN  !  NFIXMM=KFREEAB FOR IFEPTYP=2
         DO I=1,NTODO
            I1    = NONLS1(1,I)
            I2    = NONLS1(2,I)
            IYESA  = 0
            JYESA  = 0
            DO III=1,NFIXMM
               IA = IFIXMM(III)
               IF(I1.EQ.IA) IYESA = 1
               IF(I2.EQ.IA) JYESA = 1
            ENDDO
            IF(IYESA.EQ.1.AND.JYESA.EQ.1) THEN
               NTODOPMA             = NTODOPMA + 1
               NONLSPMA(1,NTODOPMA) = I1
               NONLSPMA(2,NTODOPMA) = I2
            END IF
         ENDDO
         NTODOPMB = NTODOPMA
         CALL ICOPY(2*NTODOPMA,NONLSPMA,1,NONLSPMB,1)
      END IF
C
      IF(ISTEP.EQ.0) CALL TIMIT(1)
      IF(MASWRK.AND.ISTEP.EQ.0) WRITE(IW,*)' '
C
 299  CONTINUE
C
C     -- CHECK IF A NEW QM-MM LIST IS NEEDED --
C
      IF(MXLIST1.EQ.MAGIC1) RETURN
C
      IF(NAT.GT.0) THEN
      DISMAX = -1.0D+30
      DO IFFAT=1, NFFAT
         DX = CORD(1,IFFAT)-CORDSVQ(1,IFFAT)
         DY = CORD(2,IFFAT)-CORDSVQ(2,IFFAT)
         DZ = CORD(3,IFFAT)-CORDSVQ(3,IFFAT)
         DIS= DX*DX+DY*DY+DZ*DZ
         DISMAX = MAX(DISMAX,DIS)
      ENDDO
      QDX = QMCX-QMCXSV
      QDY = QMCY-QMCYSV
      QDZ = QMCZ-QMCZSV
      DIS = QDX*QDX+QDY*QDY+QDZ*QDZ
      DISMAX = MAX(DISMAX,DIS)
      IF(ISTEP.GT.0.AND.DISMAX.LT.(BUFWID2/2.0D+00)**2) GOTO 999
C
      IF(MASWRK.AND.ISTEP.EQ.0) THEN
         WRITE(IW,'(1X,A,I10)')
     *           'UPDATING QM/MM NEIGHBOR LIST AT STEP ',ISTEP
         WRITE(IW,*) ' '
      END IF
      IF(MASWRK.AND.ISTEP.LE.10000) WRITE(IP,'(1X,A,I10)')
     *           'UPDATING QM/MM NEIGHBOR LIST AT STEP ',ISTEP
      KQMMM(1) = KQMMM(2)
      KQMMM(2) = KQMMM(3)
      KQMMM(3) = KQMMM(4)
      KQMMM(4) = KQMMM(5)
      KQMMM(5) = ISTEP
      KQMMM(6) = (KQMMM(5) - KQMMM(1))/4
      KQMMM(6) = KQMMM(6) - NINT(SQRT(DBLE(KQMMM(6))))
      CALL DCOPY(3*NFFAT,CORD,1,CORDSVQ,1)
      QMCXSV = QMCX
      QMCYSV = QMCY
      QMCZSV = QMCZ
      END IF
C
      IF(ISWITCH.LE.1) THEN
C     -- GENERATE QM-MM LIST: MM ATOMS CLOSE TO QMCX,Y,Z
C        THIS IS DONE ON EVERY CPU BECAUSE
C        WE NEED THE CORRECT AND TOTAL NTODOQ EVERYWHERE.
C        PARALLEL IS AT THE QM LEVEL, NOT MM.
C
      RSMALLQ  = SWRBQ + BUFWID2
      RSMALLQ2 = RSMALLQ*RSMALLQ
C
      NTODOQ  = 0
      IF(NAT.GT.0) THEN
      DO 300 IFFAT=1, NFFAT
         XGRP  = CORD(1,IFFAT) - QMCX
         YGRP  = CORD(2,IFFAT) - QMCY
         ZGRP  = CORD(3,IFFAT) - QMCZ
         PBCX  = XBOX * ANINT(XGRP*ONEXBOX)
         PBCY  = YBOX * ANINT(YGRP*ONEYBOX)
         PBCZ  = ZBOX * ANINT(ZGRP*ONEZBOX)
         XGRP  = XGRP - PBCX
         YGRP  = YGRP - PBCY
         ZGRP  = ZGRP - PBCZ
         IF(ABS(XGRP).GT.RSMALLQ) GOTO 300
         IF(ABS(YGRP).GT.RSMALLQ) GOTO 300
         IF(ABS(ZGRP).GT.RSMALLQ) GOTO 300
         R2    = XGRP*XGRP+YGRP*YGRP+ZGRP*ZGRP
         IF(R2.GT.RSMALLQ2) GOTO 300
         NTODOQ = NTODOQ + 1
         NONLSTQ(NTODOQ) = IFFAT
 300  CONTINUE
      END IF
      END IF
C
      IF(ISWITCH.EQ.2) THEN
C     -- GENERATE QM-MM LIST: MM ATOMS CLOSE TO QM ATOMS
C
      RSMALLQ  = SWRB + BUFWID2
      RSMALLQ2 = RSMALLQ*RSMALLQ
C
      NTODOQ  = 0
      IF(NAT.GT.0) THEN
      DO 310 IFFAT=1, NFFAT
         RMIN2 = 1.0D+30
         DO IAT=1,NAT
            XGRP  = CORD(1,IFFAT)-C(1,IAT)
            YGRP  = CORD(2,IFFAT)-C(2,IAT)
            ZGRP  = CORD(3,IFFAT)-C(3,IAT)
            PBCX  = XBOX * ANINT(XGRP*ONEXBOX)
            PBCY  = YBOX * ANINT(YGRP*ONEYBOX)
            PBCZ  = ZBOX * ANINT(ZGRP*ONEZBOX)
            XGRP  = XGRP - PBCX
            YGRP  = YGRP - PBCY
            ZGRP  = ZGRP - PBCZ
            R2    = XGRP*XGRP+YGRP*YGRP+ZGRP*ZGRP
            RMIN2 = MIN(RMIN2,R2)
         ENDDO
         IF(RMIN2.GT.RSMALLQ2) GOTO 310
         NTODOQ = NTODOQ + 1
         NONLSTQ(NTODOQ) = IFFAT
 310  CONTINUE
      END IF
      END IF
C
 999  CONTINUE
C
      CALL FLSHBF(IW)
      RETURN
      END
C*MODULE QUANPOB  *DECK E00012
!>
!> @brief    bond stretching energy
!>
!> @author   Nandun Thellamurege
!>           - Jan 2011
!>
!> @details  force field bond stretching energy
!>
      SUBROUTINE E00012(CORD,FFGRD,BOND0,FCBOND,IPAIR,CORDB,
     *                  LSBONDPMA,LSBONDPMB)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      INTEGER P1, P2
C
      LOGICAL GOPARR,DSKWRK,MASWRK
C
      PARAMETER (TWO=2.0D+00)
      PARAMETER (THREE=3.0D+00)
      PARAMETER (FOUR=4.0D+00)
      PARAMETER (ONE=1.0D+00)
C
      DIMENSION CORD(3,*),FFGRD(3,*),IPAIR(2,*),BOND0(*),FCBOND(*),
     *          CORDB(3,*),LSBONDPMA(*),LSBONDPMB(*)
C
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
      COMMON /FFMDPA/ DT,DT2,TEMP0,PRES0,POLTOL,VIR(3),PMEAN,VOLAV,
     *                ENPAV,ENKAV,TEMPAV,BERENDT,BERENDP,VELMAX,
     *                PMEANX,PMEANY,PMEANZ,NSTEP,KMASTER,KOUTACT(2),
     *                IHESS,INTALG,ITSTAT,IPSTAT,JOUT,KOUT,LOUT
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
      COMMON /FFRMSD/ DIPT(8),TIMGYRA,TIMRALL,
     *                NATPDB,NGYRA,NDIEL,NRALL,
     *                LFFRALL0,NRMSD,LFFRMSD0,KFREEAB(201),
     *                NRIJMM,IJRMM(2,100),NRIJQM,IJRQM(2,100),
     *                NAIJKMM,IJKMM(3,100),NAIJKQM,IJKQM(3,100),
     *                NFIXMM,IFIXMM(200),NFIXQM,IFIXQM(200)
      COMMON /FFTYPE/ WT14LJ,WT14CH,C3BOND,C4BOND,C3ANGL,
     *                NFFTYP,NFFFILE,LJQMMM,LJQM,INTCHG,
     *                LJSIGMA,JTOPFILE(90),JPARFILE(90),
     *                JTOPAMIA(90),JTOPNTER(90),JTOPCTER(90),
     *                JTOPNUCA(90),JPARFIL2(90),JPARFIL3(90)
      COMMON /IOFILE/ IR,IW,IP,IJK,IJKT,IDAF,NAV,IODA(950)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
C
      INTEGER, PARAMETER :: K15 = SELECTED_INT_KIND(15)
      INTEGER, PARAMETER :: MAGIC3 =
     *   TRANSFER(-31415926535897932_K15,1)
C
C     NANDUN THELLAMUREGE, JAN 2011, LINCOLN
C     FENGCHAO CUI, HUI LI, MAY 2012
C     HONGBO ZHU, HUI LI, NOV 1, 2012
C
C     FORMULA:  E = K*(R-R0)**2
C               E = K*(R-R0)**2*[1+C3*(R-R0)+C4*(R-R0)**2]
C
      NERROR = 0
C
      PMF1BD = 0.0D+00
      EN12   = 0.0D+00
      DO 100 IBOND=L1BOND,L2BOND
         FCONST  = FCBOND(IBOND)
         IF(FCONST.EQ.0.0D+00) GOTO 100
         P1      = IPAIR(1,IBOND)
         P2      = IPAIR(2,IBOND)
         R0      = BOND0(IBOND)
         X       = CORD(1,P1) - CORD(1,P2)
         Y       = CORD(2,P1) - CORD(2,P2)
         Z       = CORD(3,P1) - CORD(3,P2)
         R2      = X*X + Y*Y + Z*Z
         IF(R2.LE.1.0D-10) GOTO 100
         IF(R2.GT.100.0D+00) THEN
            WRITE(*,'(1X,A/1X,A,1X,I8,1X,A,1X,I8,1X,A,1X,F10.4,1X,A/)')
     *      'ERROR: TOO LONG BOND LENGTH.',
     *      'ATOM1=',P1,'ATOM2=',P2,'R=',SQRT(R2),' BOHR.'
            NERROR = NERROR + 1
         END IF
         R       = SQRT(R2)
         DR      = R-R0
         DR2     = DR*DR
         DR3     = DR*DR2
         EN12    = EN12+FCONST*DR2*
     *             (ONE+C3BOND*DR+C4BOND*DR2)
         DUM     = FCONST*(TWO*DR+THREE*C3BOND*DR2+
     *                     FOUR*C4BOND*DR3)/R
         DEX     = DUM*X
         DEY     = DUM*Y
         DEZ     = DUM*Z
         FFGRD(1,P1)=FFGRD(1,P1) + DEX
         FFGRD(2,P1)=FFGRD(2,P1) + DEY
         FFGRD(3,P1)=FFGRD(3,P1) + DEZ
         FFGRD(1,P2)=FFGRD(1,P2) - DEX
         FFGRD(2,P2)=FFGRD(2,P2) - DEY
         FFGRD(3,P2)=FFGRD(3,P2) - DEZ
         VIR(1)     =VIR(1)      + DEX*X
         VIR(2)     =VIR(2)      + DEY*Y
         VIR(3)     =VIR(3)      + DEZ*Z
         IYES = 0
         JYES = 0
         DO KFIX=1,NFIXMM
            IF(P1.EQ.IFIXMM(KFIX)) IYES = 1
            IF(P2.EQ.IFIXMM(KFIX)) JYES = 1
         ENDDO
         IF(IYES.EQ.1.AND.JYES.EQ.1) THEN
            VIR(1)    =VIR(1)     - DEX*X
            VIR(2)    =VIR(2)     - DEY*Y
            VIR(3)    =VIR(3)     - DEZ*Z
         END IF
 100  CONTINUE
C
      DO 200 III=L1BONDPMA,L2BONDPMA
         IBOND   = LSBONDPMA(III)
         FCONST  = FCBOND(IBOND)
         IF(FCONST.EQ.0.0D+00) GOTO 200
         P1      = IPAIR(1,IBOND)
         P2      = IPAIR(2,IBOND)
         R0      = BOND0(IBOND)
         X       = CORD(1,P1) - CORD(1,P2)
         Y       = CORD(2,P1) - CORD(2,P2)
         Z       = CORD(3,P1) - CORD(3,P2)
         R2      = X*X + Y*Y + Z*Z
         IF(R2.LE.1.0D-10) GOTO 200
         R       = SQRT(R2)
         DR      = R-R0
         DR2     = DR*DR
         PMF1BD  = PMF1BD - FCONST*DR2*
     *             (ONE+C3BOND*DR+C4BOND*DR2)
 200  CONTINUE
C
      DO 210 III=L1BONDPMB,L2BONDPMB
         IBOND   = LSBONDPMB(III)
         FCONST  = FCBOND(IBOND)
         IF(FCONST.EQ.0.0D+00) GOTO 210
         P1      = IPAIR(1,IBOND)
         P2      = IPAIR(2,IBOND)
         R0      = BOND0(IBOND)
         X       = CORD(1,P1) - CORD(1,P2)
         Y       = CORD(2,P1) - CORD(2,P2)
         Z       = CORD(3,P1) - CORD(3,P2)
         IF(IFEPTYP.EQ.2) THEN  !  NFIXMM=KFREEAB FOR IFEPTYP=2
            DO KFIX=1,NFIXMM
               IF(P1.EQ.IFIXMM(KFIX)) THEN
                  X = X - CORD(1,P1) + CORDB(1,P1)
                  Y = Y - CORD(2,P1) + CORDB(2,P1)
                  Z = Z - CORD(3,P1) + CORDB(3,P1)
               END IF
               IF(P2.EQ.IFIXMM(KFIX)) THEN
                  X = X + CORD(1,P2) - CORDB(1,P2)
                  Y = Y + CORD(2,P2) - CORDB(2,P2)
                  Z = Z + CORD(3,P2) - CORDB(3,P2)
               END IF
            ENDDO
         END IF
         R2      = X*X + Y*Y + Z*Z
         IF(R2.LE.1.0D-10) GOTO 210
         R       = SQRT(R2)
         R       = SQRT(R2)
         DR      = R-R0
         DR2     = DR*DR
         PMF1BD  = PMF1BD + FCONST*DR2*
     *             (ONE+C3BOND*DR+C4BOND*DR2)
 210  CONTINUE
C
      IF(GOPARR) CALL DDI_GSUMI(2422,NERROR,1)
      IF(NERROR.GT.0) THEN
         IF(LOUT.EQ.MAGIC3) CALL ABRT
         LOUT = MAGIC3
         JOUT = 1
         KOUT = 1
      END IF
C
      RETURN
      END
C*MODULE QUANPOB  *DECK E00123
!>
!> @brief    angle bending energy
!>
!> @author   Nandun Thellamurege
!>           - Jan 2011
!>
!> @details  force field angle bending energy
!>
      SUBROUTINE E00123(CORD,FFGRD,ANGL0,FCANGL,KLIST,CORDB,
     *                  LSANGLPMA,LSANGLPMB)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      INTEGER P1, P2, P3
C
      PARAMETER (PI=3.14159265358979323846264338D+00)
      PARAMETER (PT5=0.5D+00)
      PARAMETER (ZERO=0.0D+00)
      PARAMETER (ONE=1.0D+00)
      PARAMETER (TWO=2.0D+00)
      PARAMETER (THREE=3.0D+00)
C
      DIMENSION CORD(3,*),FFGRD(3,*),KLIST(3,*),ANGL0(*),FCANGL(*),
     *          CORDB(3,*),LSANGLPMA(*),LSANGLPMB(*)
C
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
      COMMON /FFMDPA/ DT,DT2,TEMP0,PRES0,POLTOL,VIR(3),PMEAN,VOLAV,
     *                ENPAV,ENKAV,TEMPAV,BERENDT,BERENDP,VELMAX,
     *                PMEANX,PMEANY,PMEANZ,NSTEP,KMASTER,KOUTACT(2),
     *                IHESS,INTALG,ITSTAT,IPSTAT,JOUT,KOUT,LOUT
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
      COMMON /FFRMSD/ DIPT(8),TIMGYRA,TIMRALL,
     *                NATPDB,NGYRA,NDIEL,NRALL,
     *                LFFRALL0,NRMSD,LFFRMSD0,KFREEAB(201),
     *                NRIJMM,IJRMM(2,100),NRIJQM,IJRQM(2,100),
     *                NAIJKMM,IJKMM(3,100),NAIJKQM,IJKQM(3,100),
     *                NFIXMM,IFIXMM(200),NFIXQM,IFIXQM(200)
      COMMON /FFTYPE/ WT14LJ,WT14CH,C3BOND,C4BOND,C3ANGL,
     *                NFFTYP,NFFFILE,LJQMMM,LJQM,INTCHG,
     *                LJSIGMA,JTOPFILE(90),JPARFILE(90),
     *                JTOPAMIA(90),JTOPNTER(90),JTOPCTER(90),
     *                JTOPNUCA(90),JPARFIL2(90),JPARFIL3(90)
C
C     NANDUN THELLAMUREGE, JAN 2011, LINCOLN
C     FENGCHAO CUI, HUI LI, MAY 2012
C     HONGBO ZHU, HUI LI, NOV 1, 2012
C
C     FORMULA:  E = K*[(A-A0)**2]*[1 + C3*(A-A0)]
C                 + 4.0*K*(A-A0-20)**2  (TOO LARGE POSITIVE DEVIATION)
C                 + 4.0*K*(A-A0+20)**2  (TOO LARGE NEGATIVE DEVIATION)
C     FORMULA:  E = 2*K*(1+COS(A))       (FOR A0 = 180 DEGREE)
C
      PMF1AG = ZERO
      EN123  = ZERO
      DO 100 IANGL=L1ANGL,L2ANGL
         FCONST =FCANGL(IANGL)
         IF(FCONST.EQ.ZERO) GOTO 100
         P1     =KLIST(1,IANGL)
         P2     =KLIST(2,IANGL)
         P3     =KLIST(3,IANGL)
         X13    =CORD(1,P1)-CORD(1,P3)
         Y13    =CORD(2,P1)-CORD(2,P3)
         Z13    =CORD(3,P1)-CORD(3,P3)
         X12    =CORD(1,P1)-CORD(1,P2)
         Y12    =CORD(2,P1)-CORD(2,P2)
         Z12    =CORD(3,P1)-CORD(3,P2)
         X23    =CORD(1,P2)-CORD(1,P3)
         Y23    =CORD(2,P2)-CORD(2,P3)
         Z23    =CORD(3,P2)-CORD(3,P3)
C
         R13R13 =X13*X13+Y13*Y13+Z13*Z13
         R12R12 =X12*X12+Y12*Y12+Z12*Z12
         R23R23 =X23*X23+Y23*Y23+Z23*Z23
         R12    =SQRT(R12R12)
         R23    =SQRT(R23R23)
C
         ONEBC  =ONE/(R12*R23)
         COSA   =(R12R12 + R23R23 - R13R13)*ONEBC*PT5
C
         IF(ABS(ANGL0(IANGL)-PI).GT.0.053D+00) THEN
            IF(COSA.GT. ONE) COSA = ONE
            IF(COSA.LT.-ONE) COSA =-ONE
            ALPHA  =ACOS(COSA)
            SINA   =SQRT(ABS(ONE - COSA*COSA))
            DELT   =ALPHA - ANGL0(IANGL)
            DELT2  =DELT*DELT
            EN123  =EN123 + FCONST*DELT2*(ONE+C3ANGL*DELT)
            DUMY   =-FCONST*(TWO*DELT+THREE*C3ANGL*DELT2)/SINA
         ELSE
            EN123 = EN123 + TWO*FCONST*(ONE+COSA)
            DUMY   = TWO*FCONST
         END IF
C
C        - CALCULATE BENDING ENERGY GRADIENTS
C
         DR23X2 =-X23*ONEBC
         DR23Y2 =-Y23*ONEBC
         DR23Z2 =-Z23*ONEBC
         DR12X1A=-COSA*X12/R12R12
         DR12Y1A=-COSA*Y12/R12R12
         DR12Z1A=-COSA*Z12/R12R12
C
         DR23X2A=COSA*X23/R23R23
         DR23Y2A=COSA*Y23/R23R23
         DR23Z2A=COSA*Z23/R23R23
         DR12X1 =X12*ONEBC
         DR12Y1 =Y12*ONEBC
         DR12Z1 =Z12*ONEBC
C
         DEX1   = DUMY*(DR12X1A+DR23X2)
         DEY1   = DUMY*(DR12Y1A+DR23Y2)
         DEZ1   = DUMY*(DR12Z1A+DR23Z2)
         FFGRD(1,P1)=FFGRD(1,P1) + DEX1
         FFGRD(2,P1)=FFGRD(2,P1) + DEY1
         FFGRD(3,P1)=FFGRD(3,P1) + DEZ1
C
         DEX3   = DUMY*(DR23X2A+DR12X1)
         DEY3   = DUMY*(DR23Y2A+DR12Y1)
         DEZ3   = DUMY*(DR23Z2A+DR12Z1)
         FFGRD(1,P3)=FFGRD(1,P3)+DEX3
         FFGRD(2,P3)=FFGRD(2,P3)+DEY3
         FFGRD(3,P3)=FFGRD(3,P3)+DEZ3
C
         FFGRD(1,P2)=FFGRD(1,P2)-DEX1-DEX3
         FFGRD(2,P2)=FFGRD(2,P2)-DEY1-DEY3
         FFGRD(3,P2)=FFGRD(3,P2)-DEZ1-DEZ3
C
         VIR(1)     =VIR(1)  + DEX1*X12 - DEX3*X23
         VIR(2)     =VIR(2)  + DEY1*Y12 - DEY3*Y23
         VIR(3)     =VIR(3)  + DEZ1*Z12 - DEZ3*Z23
         IYES = 0
         JYES = 0
         KYES = 0
         DO KFIX=1,NFIXMM
            IF(P1.EQ.IFIXMM(KFIX)) IYES = 1
            IF(P2.EQ.IFIXMM(KFIX)) JYES = 1
            IF(P3.EQ.IFIXMM(KFIX)) KYES = 1
         ENDDO
         IF(IYES.EQ.1.AND.JYES.EQ.1) THEN
            VIR(1)    =VIR(1)     - DEX1*X12
            VIR(2)    =VIR(2)     - DEY1*Y12
            VIR(3)    =VIR(3)     - DEZ1*Z12
         END IF
         IF(JYES.EQ.1.AND.KYES.EQ.1) THEN
            VIR(1)    =VIR(1)     + DEX3*X23
            VIR(2)    =VIR(2)     + DEY3*Y23
            VIR(3)    =VIR(3)     + DEZ3*Z23
         END IF
 100  CONTINUE
C
      DO 200 III=L1ANGLPMA,L2ANGLPMA
         IANGL  =LSANGLPMA(III)
         FCONST =FCANGL(IANGL)
         IF(FCONST.EQ.ZERO) GOTO 200
         P1     =KLIST(1,IANGL)
         P2     =KLIST(2,IANGL)
         P3     =KLIST(3,IANGL)
         X13    =CORD(1,P1)-CORD(1,P3)
         Y13    =CORD(2,P1)-CORD(2,P3)
         Z13    =CORD(3,P1)-CORD(3,P3)
         X12    =CORD(1,P1)-CORD(1,P2)
         Y12    =CORD(2,P1)-CORD(2,P2)
         Z12    =CORD(3,P1)-CORD(3,P2)
         X23    =CORD(1,P2)-CORD(1,P3)
         Y23    =CORD(2,P2)-CORD(2,P3)
         Z23    =CORD(3,P2)-CORD(3,P3)
C
         R13R13 =X13*X13+Y13*Y13+Z13*Z13
         R12R12 =X12*X12+Y12*Y12+Z12*Z12
         R23R23 =X23*X23+Y23*Y23+Z23*Z23
         R12    =SQRT(R12R12)
         R23    =SQRT(R23R23)
C
         ONEBC  =ONE/(R12*R23)
         COSA   =(R12R12 + R23R23 - R13R13)*ONEBC*PT5
C
         IF(ABS(ANGL0(IANGL)-PI).GT.0.053D+00) THEN
            IF(COSA.GT. ONE) COSA = ONE
            IF(COSA.LT.-ONE) COSA =-ONE
            ALPHA  =ACOS(COSA)
            DELT   =ALPHA - ANGL0(IANGL)
            DELT2  =DELT*DELT
            PMF1AG =PMF1AG - FCONST*DELT2*(ONE+C3ANGL*DELT)
         ELSE
            PMF1AG = PMF1AG - TWO*FCONST*(ONE+COSA)
         END IF
 200  CONTINUE
C
      DO 210 III=L1ANGLPMB,L2ANGLPMB
         IANGL  =LSANGLPMB(III)
         FCONST =FCANGL(IANGL)
         IF(FCONST.EQ.ZERO) GOTO 210
         P1     =KLIST(1,IANGL)
         P2     =KLIST(2,IANGL)
         P3     =KLIST(3,IANGL)
         X13    =CORD(1,P1)-CORD(1,P3)
         Y13    =CORD(2,P1)-CORD(2,P3)
         Z13    =CORD(3,P1)-CORD(3,P3)
         X12    =CORD(1,P1)-CORD(1,P2)
         Y12    =CORD(2,P1)-CORD(2,P2)
         Z12    =CORD(3,P1)-CORD(3,P2)
         X23    =CORD(1,P2)-CORD(1,P3)
         Y23    =CORD(2,P2)-CORD(2,P3)
         Z23    =CORD(3,P2)-CORD(3,P3)
C
         IF(IFEPTYP.EQ.2) THEN  !  NFIXMM=KFREEAB FOR IFEPTYP=2
            DO KFIX=1,NFIXMM
               IF(P1.EQ.IFIXMM(KFIX)) THEN
                  X13 = X13 - CORD(1,P1) + CORDB(1,P1)
                  Y13 = Y13 - CORD(2,P1) + CORDB(2,P1)
                  Z13 = Z13 - CORD(3,P1) + CORDB(3,P1)
                  X12 = X12 - CORD(1,P1) + CORDB(1,P1)
                  Y12 = Y12 - CORD(2,P1) + CORDB(2,P1)
                  Z12 = Z12 - CORD(3,P1) + CORDB(3,P1)
               END IF
               IF(P2.EQ.IFIXMM(KFIX)) THEN
                  X12 = X12 + CORD(1,P2) - CORDB(1,P2)
                  Y12 = Y12 + CORD(2,P2) - CORDB(2,P2)
                  Z12 = Z12 + CORD(3,P2) - CORDB(3,P2)
                  X23 = X23 - CORD(1,P2) + CORDB(1,P2)
                  Y23 = Y23 - CORD(2,P2) + CORDB(2,P2)
                  Z23 = Z23 - CORD(3,P2) + CORDB(3,P2)
               END IF
               IF(P3.EQ.IFIXMM(KFIX)) THEN
                  X13 = X13 + CORD(1,P3) - CORDB(1,P3)
                  Y13 = Y13 + CORD(2,P3) - CORDB(2,P3)
                  Z13 = Z13 + CORD(3,P3) - CORDB(3,P3)
                  X23 = X23 + CORD(1,P3) - CORDB(1,P3)
                  Y23 = Y23 + CORD(2,P3) - CORDB(2,P3)
                  Z23 = Z23 + CORD(3,P3) - CORDB(3,P3)
               END IF
            ENDDO
         END IF
C
         R13R13 =X13*X13+Y13*Y13+Z13*Z13
         R12R12 =X12*X12+Y12*Y12+Z12*Z12
         R23R23 =X23*X23+Y23*Y23+Z23*Z23
         R12    =SQRT(R12R12)
         R23    =SQRT(R23R23)
C
         ONEBC  =ONE/(R12*R23)
         COSA   =(R12R12 + R23R23 - R13R13)*ONEBC*PT5
C
         IF(ABS(ANGL0(IANGL)-PI).GT.0.053D+00) THEN
            IF(COSA.GT. ONE) COSA = ONE
            IF(COSA.LT.-ONE) COSA =-ONE
            ALPHA  =ACOS(COSA)
            DELT   =ALPHA - ANGL0(IANGL)
            DELT2  =DELT*DELT
            PMF1AG =PMF1AG + FCONST*DELT2*(ONE+C3ANGL*DELT)
         ELSE
            PMF1AG = PMF1AG + TWO*FCONST*(ONE+COSA)
         END IF
 210  CONTINUE
C
      RETURN
      END
C*MODULE QUANPOB  *DECK E12312
!>
!> @brief    MMFF94 bond-angle energy
!>
!> @author   Hongbo Zhu
!>           - Nov 2012
!>
!> @details  force field bond-angle coupled term
!>
      SUBROUTINE E12312(CORD,FFGRD,ANGL0,KLIST,BOND0,FCSTBD,
     *                  KBLST,CORDB,LSANGLPMA,LSANGLPMB)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      PARAMETER (PT5=0.5D+00)
      PARAMETER (ZERO=0.0D+00)
      PARAMETER (ONE=1.0D+00)
C
      DIMENSION CORD(3,*),KLIST(3,*),FCSTBD(2,*),ANGL0(*),
     *          BOND0(*),FFGRD(3,*),KBLST(2,*),CORDB(3,*),
     *          LSANGLPMA(*),LSANGLPMB(*)
C
      INTEGER P1, P2, P3
C
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
      COMMON /FFMDPA/ DT,DT2,TEMP0,PRES0,POLTOL,VIR(3),PMEAN,VOLAV,
     *                ENPAV,ENKAV,TEMPAV,BERENDT,BERENDP,VELMAX,
     *                PMEANX,PMEANY,PMEANZ,NSTEP,KMASTER,KOUTACT(2),
     *                IHESS,INTALG,ITSTAT,IPSTAT,JOUT,KOUT,LOUT
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
      COMMON /FFRMSD/ DIPT(8),TIMGYRA,TIMRALL,
     *                NATPDB,NGYRA,NDIEL,NRALL,
     *                LFFRALL0,NRMSD,LFFRMSD0,KFREEAB(201),
     *                NRIJMM,IJRMM(2,100),NRIJQM,IJRQM(2,100),
     *                NAIJKMM,IJKMM(3,100),NAIJKQM,IJKQM(3,100),
     *                NFIXMM,IFIXMM(200),NFIXQM,IFIXQM(200)
C
C     HONGBO ZHU, HUI LI, NOV 5, 2012, LINCOLN
C
      EN12312=ZERO
C
C     IF(NFFTYP/10000.NE.5) RETURN
C
      DO 100 IANGL=L1ANGL,L2ANGL
         IF(FCSTBD(1,IANGL).EQ.ZERO.AND.FCSTBD(2,IANGL).EQ.ZERO)
     *   GOTO 100
         A0     =ANGL0(IANGL)
         P1     =KLIST(1,IANGL)
         P2     =KLIST(2,IANGL)
         P3     =KLIST(3,IANGL)
         X13    =CORD(1,P1)-CORD(1,P3)
         Y13    =CORD(2,P1)-CORD(2,P3)
         Z13    =CORD(3,P1)-CORD(3,P3)
         X12    =CORD(1,P1)-CORD(1,P2)
         Y12    =CORD(2,P1)-CORD(2,P2)
         Z12    =CORD(3,P1)-CORD(3,P2)
         X23    =CORD(1,P2)-CORD(1,P3)
         Y23    =CORD(2,P2)-CORD(2,P3)
         Z23    =CORD(3,P2)-CORD(3,P3)
C
         R13R13 =X13*X13+Y13*Y13+Z13*Z13
         R12R12 =X12*X12+Y12*Y12+Z12*Z12
         R23R23 =X23*X23+Y23*Y23+Z23*Z23
         R12    =SQRT(R12R12)
         R23    =SQRT(R23R23)
C
         ONEBC  =ONE/(R12*R23)
         COSA   =(R12R12 + R23R23 - R13R13)*ONEBC*PT5
         IF(COSA.GT. ONE) COSA = ONE
         IF(COSA.LT.-ONE) COSA =-ONE
C
         ALPHA  =ACOS(COSA)
         SINA   =SQRT(ABS(ONE - COSA*COSA))
         DELT   =ALPHA - A0
         IBOND  =KBLST(1,IANGL)
         JBOND  =KBLST(2,IANGL)
         R0IJ   =BOND0(IBOND)
         DRIJ   =R12-R0IJ
         R0JK   =BOND0(JBOND)
         DRJK   =R23-R0JK
         CCC1   =FCSTBD(1,IANGL)*DRIJ+FCSTBD(2,IANGL)*DRJK
         EN12312=EN12312+CCC1*DELT
         DUMY1  =DELT*FCSTBD(1,IANGL)/R12
         DUMY2  =-CCC1/SINA
         DUMY3  =DELT*FCSTBD(2,IANGL)/R23
C
C        - CALCULATE BENDING ENERGY GRADIENTS
C
         DR23X2 =-X23*ONEBC
         DR23Y2 =-Y23*ONEBC
         DR23Z2 =-Z23*ONEBC
         DR12X1A=-COSA*X12/R12R12
         DR12Y1A=-COSA*Y12/R12R12
         DR12Z1A=-COSA*Z12/R12R12
C
         DR23X2A=COSA*X23/R23R23
         DR23Y2A=COSA*Y23/R23R23
         DR23Z2A=COSA*Z23/R23R23
         DR12X1 =X12*ONEBC
         DR12Y1 =Y12*ONEBC
         DR12Z1 =Z12*ONEBC
C
         DEX1   = DUMY1*X12
         DEY1   = DUMY1*Y12
         DEZ1   = DUMY1*Z12
         FFGRD(1,P1)=FFGRD(1,P1) + DEX1
         FFGRD(2,P1)=FFGRD(2,P1) + DEY1
         FFGRD(3,P1)=FFGRD(3,P1) + DEZ1
         DEX3   = DUMY3*X23
         DEY3   = DUMY3*Y23
         DEZ3   = DUMY3*Z23
         FFGRD(1,P3)=FFGRD(1,P3) - DEX3
         FFGRD(2,P3)=FFGRD(2,P3) - DEY3
         FFGRD(3,P3)=FFGRD(3,P3) - DEZ3
C
         FFGRD(1,P2)=FFGRD(1,P2)-DEX1+DEX3
         FFGRD(2,P2)=FFGRD(2,P2)-DEY1+DEY3
         FFGRD(3,P2)=FFGRD(3,P2)-DEZ1+DEZ3
C
         VIR(1)     =VIR(1)  + DEX1*X12 + DEX3*X23
         VIR(2)     =VIR(2)  + DEY1*Y12 + DEY3*Y23
         VIR(3)     =VIR(3)  + DEZ1*Z12 + DEZ3*Z23
         IYES = 0
         JYES = 0
         KYES = 0
         DO KFIX=1,NFIXMM
            IF(P1.EQ.IFIXMM(KFIX)) IYES = 1
            IF(P2.EQ.IFIXMM(KFIX)) JYES = 1
            IF(P3.EQ.IFIXMM(KFIX)) KYES = 1
         ENDDO
         IF(IYES.EQ.1.AND.JYES.EQ.1) THEN
            VIR(1)    =VIR(1)     - DEX1*X12
            VIR(2)    =VIR(2)     - DEY1*Y12
            VIR(3)    =VIR(3)     - DEZ1*Z12
         END IF
         IF(JYES.EQ.1.AND.KYES.EQ.1) THEN
            VIR(1)    =VIR(1)     - DEX3*X23
            VIR(2)    =VIR(2)     - DEY3*Y23
            VIR(3)    =VIR(3)     - DEZ3*Z23
         END IF
C
         DEX1   = DUMY2*(DR12X1A+DR23X2)
         DEY1   = DUMY2*(DR12Y1A+DR23Y2)
         DEZ1   = DUMY2*(DR12Z1A+DR23Z2)
         FFGRD(1,P1)=FFGRD(1,P1) + DEX1
         FFGRD(2,P1)=FFGRD(2,P1) + DEY1
         FFGRD(3,P1)=FFGRD(3,P1) + DEZ1
C
         DEX3   = DUMY2*(DR23X2A+DR12X1)
         DEY3   = DUMY2*(DR23Y2A+DR12Y1)
         DEZ3   = DUMY2*(DR23Z2A+DR12Z1)
         FFGRD(1,P3)=FFGRD(1,P3)+DEX3
         FFGRD(2,P3)=FFGRD(2,P3)+DEY3
         FFGRD(3,P3)=FFGRD(3,P3)+DEZ3
C
         FFGRD(1,P2)=FFGRD(1,P2)-DEX1-DEX3
         FFGRD(2,P2)=FFGRD(2,P2)-DEY1-DEY3
         FFGRD(3,P2)=FFGRD(3,P2)-DEZ1-DEZ3
C
         VIR(1)     =VIR(1)  + DEX1*X12 - DEX3*X23
         VIR(2)     =VIR(2)  + DEY1*Y12 - DEY3*Y23
         VIR(3)     =VIR(3)  + DEZ1*Z12 - DEZ3*Z23
         IYES = 0
         JYES = 0
         KYES = 0
         DO KFIX=1,NFIXMM
            IF(P1.EQ.IFIXMM(KFIX)) IYES = 1
            IF(P2.EQ.IFIXMM(KFIX)) JYES = 1
            IF(P3.EQ.IFIXMM(KFIX)) KYES = 1
         ENDDO
         IF(IYES.EQ.1.AND.JYES.EQ.1) THEN
            VIR(1)    =VIR(1)     - DEX1*X12
            VIR(2)    =VIR(2)     - DEY1*Y12
            VIR(3)    =VIR(3)     - DEZ1*Z12
         END IF
         IF(JYES.EQ.1.AND.KYES.EQ.1) THEN
            VIR(1)    =VIR(1)     + DEX3*X23
            VIR(2)    =VIR(2)     + DEY3*Y23
            VIR(3)    =VIR(3)     + DEZ3*Z23
         END IF
 100  CONTINUE
C
      DO 200 III=L1ANGLPMA,L2ANGLPMA
         IANGL  =LSANGLPMA(III)
         IF(FCSTBD(1,IANGL).EQ.ZERO.AND.FCSTBD(2,IANGL).EQ.ZERO)
     *   GOTO 200
         P1     =KLIST(1,IANGL)
         P2     =KLIST(2,IANGL)
         P3     =KLIST(3,IANGL)
         X13    =CORD(1,P1)-CORD(1,P3)
         Y13    =CORD(2,P1)-CORD(2,P3)
         Z13    =CORD(3,P1)-CORD(3,P3)
         X12    =CORD(1,P1)-CORD(1,P2)
         Y12    =CORD(2,P1)-CORD(2,P2)
         Z12    =CORD(3,P1)-CORD(3,P2)
         X23    =CORD(1,P2)-CORD(1,P3)
         Y23    =CORD(2,P2)-CORD(2,P3)
         Z23    =CORD(3,P2)-CORD(3,P3)
C
         R13R13 =X13*X13+Y13*Y13+Z13*Z13
         R12R12 =X12*X12+Y12*Y12+Z12*Z12
         R23R23 =X23*X23+Y23*Y23+Z23*Z23
         R12    =SQRT(R12R12)
         R23    =SQRT(R23R23)
C
         ONEBC  =ONE/(R12*R23)
         COSA   =(R12R12 + R23R23 - R13R13)*ONEBC*PT5
         IF(COSA.GT. ONE) COSA = ONE
         IF(COSA.LT.-ONE) COSA =-ONE
C
         ALPHA  =ACOS(COSA)
         SINA   =SQRT(ABS(ONE - COSA*COSA))
         DELT   =ALPHA - A0
         IBOND  =KBLST(1,IANGL)
         JBOND  =KBLST(2,IANGL)
         R0IJ   =BOND0(IBOND)
         DRIJ   =R12-R0IJ
         R0JK   =BOND0(JBOND)
         DRJK   =R23-R0JK
         CCC1   =FCSTBD(1,IANGL)*DRIJ+FCSTBD(2,IANGL)*DRJK
         PMF1AG = PMF1AG-CCC1*DELT
 200  CONTINUE
C
      DO 210 III=L1ANGLPMB,L2ANGLPMB
         IANGL  =LSANGLPMB(III)
         IF(FCSTBD(1,IANGL).EQ.ZERO.AND.FCSTBD(2,IANGL).EQ.ZERO)
     *   GOTO 210
         P1     =KLIST(1,IANGL)
         P2     =KLIST(2,IANGL)
         P3     =KLIST(3,IANGL)
         X13    =CORD(1,P1)-CORD(1,P3)
         Y13    =CORD(2,P1)-CORD(2,P3)
         Z13    =CORD(3,P1)-CORD(3,P3)
         X12    =CORD(1,P1)-CORD(1,P2)
         Y12    =CORD(2,P1)-CORD(2,P2)
         Z12    =CORD(3,P1)-CORD(3,P2)
         X23    =CORD(1,P2)-CORD(1,P3)
         Y23    =CORD(2,P2)-CORD(2,P3)
         Z23    =CORD(3,P2)-CORD(3,P3)
C
         IF(IFEPTYP.EQ.2) THEN  !  NFIXMM=KFREEAB FOR IFEPTYP=2
            DO KFIX=1,NFIXMM
               IF(P1.EQ.IFIXMM(KFIX)) THEN
                  X13 = X13 - CORD(1,P1) + CORDB(1,P1)
                  Y13 = Y13 - CORD(2,P1) + CORDB(2,P1)
                  Z13 = Z13 - CORD(3,P1) + CORDB(3,P1)
                  X12 = X12 - CORD(1,P1) + CORDB(1,P1)
                  Y12 = Y12 - CORD(2,P1) + CORDB(2,P1)
                  Z12 = Z12 - CORD(3,P1) + CORDB(3,P1)
               END IF
               IF(P2.EQ.IFIXMM(KFIX)) THEN
                  X12 = X12 + CORD(1,P2) - CORDB(1,P2)
                  Y12 = Y12 + CORD(2,P2) - CORDB(2,P2)
                  Z12 = Z12 + CORD(3,P2) - CORDB(3,P2)
                  X23 = X23 - CORD(1,P2) + CORDB(1,P2)
                  Y23 = Y23 - CORD(2,P2) + CORDB(2,P2)
                  Z23 = Z23 - CORD(3,P2) + CORDB(3,P2)
               END IF
               IF(P3.EQ.IFIXMM(KFIX)) THEN
                  X13 = X13 + CORD(1,P3) - CORDB(1,P3)
                  Y13 = Y13 + CORD(2,P3) - CORDB(2,P3)
                  Z13 = Z13 + CORD(3,P3) - CORDB(3,P3)
                  X23 = X23 + CORD(1,P3) - CORDB(1,P3)
                  Y23 = Y23 + CORD(2,P3) - CORDB(2,P3)
                  Z23 = Z23 + CORD(3,P3) - CORDB(3,P3)
               END IF
            ENDDO
         END IF
C
         R13R13 =X13*X13+Y13*Y13+Z13*Z13
         R12R12 =X12*X12+Y12*Y12+Z12*Z12
         R23R23 =X23*X23+Y23*Y23+Z23*Z23
         R12    =SQRT(R12R12)
         R23    =SQRT(R23R23)
C
         ONEBC  =ONE/(R12*R23)
         COSA   =(R12R12 + R23R23 - R13R13)*ONEBC*PT5
         IF(COSA.GT. ONE) COSA = ONE
         IF(COSA.LT.-ONE) COSA =-ONE
C
         ALPHA  =ACOS(COSA)
         SINA   =SQRT(ABS(ONE - COSA*COSA))
         DELT   =ALPHA - A0
         IBOND  =KBLST(1,IANGL)
         JBOND  =KBLST(2,IANGL)
         R0IJ   =BOND0(IBOND)
         DRIJ   =R12-R0IJ
         R0JK   =BOND0(JBOND)
         DRJK   =R23-R0JK
         CCC1   =FCSTBD(1,IANGL)*DRIJ+FCSTBD(2,IANGL)*DRJK
         PMF1AG = PMF1AG+CCC1*DELT
 210  CONTINUE
C
      RETURN
      END
C*MODULE QUANPOB  *DECK E123R4
!>
!> @brief    dihedral rotation
!>
!> @author   Nandun Thellamurege
!>           - Jan 2011
!>
!> @details  force field dihedral rotation energy
!>
      SUBROUTINE E123R4(CORD,FFGRD,VROT,GAMA,NNN,LLIST,CORDB,
     *                  LSDIHRPMA,LSDIHRPMB,FCDIHR)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      INTEGER P1, P2, P3, P4
C
      PARAMETER (PI=3.14159265358979323846264338D+00)
      PARAMETER (ONE=1.0D+00)
      PARAMETER (PT5=0.5D+00)
      PARAMETER (TWO=2.0D+00)
      PARAMETER (FOUR=4.0D+00)
      PARAMETER (THREE=3.0D+00)
C
      DIMENSION CORD(3,*),FFGRD(3,*),VROT(*),
     *          GAMA(*),NNN(*),LLIST(4,*),CORDB(3,*),LSDIHRPMA(*),
     *          LSDIHRPMB(*),FCDIHR(3,*)
C
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
      COMMON /FFMDPA/ DT,DT2,TEMP0,PRES0,POLTOL,VIR(3),PMEAN,VOLAV,
     *                ENPAV,ENKAV,TEMPAV,BERENDT,BERENDP,VELMAX,
     *                PMEANX,PMEANY,PMEANZ,NSTEP,KMASTER,KOUTACT(2),
     *                IHESS,INTALG,ITSTAT,IPSTAT,JOUT,KOUT,LOUT
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
      COMMON /FFRMSD/ DIPT(8),TIMGYRA,TIMRALL,
     *                NATPDB,NGYRA,NDIEL,NRALL,
     *                LFFRALL0,NRMSD,LFFRMSD0,KFREEAB(201),
     *                NRIJMM,IJRMM(2,100),NRIJQM,IJRQM(2,100),
     *                NAIJKMM,IJKMM(3,100),NAIJKQM,IJKQM(3,100),
     *                NFIXMM,IFIXMM(200),NFIXQM,IFIXQM(200)
C
C     NANDUN THELLAMUREGE, JAN 2011, LINCOLN
C     FENGCHAO CUI, HUI LI, MAY 2012
C
C     FORMULA:  E = VROT*[1 + COS(NNN*TOR - GAMA)]
C     FORMULA:  E = 0.5*(VROT1*(1+COS(1TOR))
C                       +VROT2*(1-COS(2TOR))
C                       +VROT3*(1+COS(3TOR)))
C
      PMF1DR   =0.0D+00
      EN123R4  =0.0D+00
C
      DO 100 I=L1DIHR,L2DIHR
C        - CONSIDER ALL NFFTYP -
         VR      = VROT(I)
         REALN   = DBLE(NNN(I))
         GAMMA   = GAMA(I)
         VROT1   = FCDIHR(1,I)
         VROT2   = FCDIHR(2,I)
         VROT3   = FCDIHR(3,I)
         IF(ABS(VR)+ABS(VROT1)+ABS(VROT2)+ABS(VROT3).EQ.0.0D+00)
     *   GOTO 100
         P1=LLIST(1,I)
         P2=LLIST(2,I)
         P3=LLIST(3,I)
         P4=LLIST(4,I)
C
         X12=CORD(1,P1)-CORD(1,P2)
         Y12=CORD(2,P1)-CORD(2,P2)
         Z12=CORD(3,P1)-CORD(3,P2)
         X13=CORD(1,P1)-CORD(1,P3)
         Y13=CORD(2,P1)-CORD(2,P3)
         Z13=CORD(3,P1)-CORD(3,P3)
         X23=CORD(1,P2)-CORD(1,P3)
         Y23=CORD(2,P2)-CORD(2,P3)
         Z23=CORD(3,P2)-CORD(3,P3)
         X34=CORD(1,P3)-CORD(1,P4)
         Y34=CORD(2,P3)-CORD(2,P4)
         Z34=CORD(3,P3)-CORD(3,P4)
C
         R12=SQRT(X12*X12+Y12*Y12+Z12*Z12)
         R23=SQRT(X23*X23+Y23*Y23+Z23*Z23)
         R34=SQRT(X34*X34+Y34*Y34+Z34*Z34)
C
         COS123=(-(X12*X23)-(Y12*Y23)-(Z12*Z23))/(R12*R23)
         COS234=(-(X23*X34)-(Y23*Y34)-(Z23*Z34))/(R23*R34)
         SIN2123= 1.0D+00-COS123*COS123
         SIN2234= 1.0D+00-COS234*COS234
         SIN123 = SQRT(ABS(SIN2123))
         SIN234 = SQRT(ABS(SIN2234))
         IF(ABS(SIN123).LT.1.0D-06) GOTO 100
         IF(ABS(SIN234).LT.1.0D-06) GOTO 100
         ONESIN = 1.0D+00/(SIN123*SIN234)
C
         COSTOR = ONESIN*(COS123*COS234-
     *                    ((+X12*X34+Y12*Y34+Z12*Z34)/(R12*R34)))
         IF(COSTOR.GT. ONE) COSTOR= ONE
         IF(COSTOR.LT.-ONE) COSTOR=-ONE
         TOR    = ACOS(COSTOR)
C        -- DIHEDRAL ANGLE IS 0 - 360 DEGREES
         XNORM  = -Y23*Z34 + Z23*Y34
         YNORM  = -Z23*X34 + X23*Z34
         ZNORM  = -X23*Y34 + Y23*X34
         DOTN12 = X12*XNORM + Y12*YNORM + Z12*ZNORM
         IF(DOTN12.LT.0.0D+00) TOR = 2.0D+00*PI - TOR
C
         IF(ABS(VR).GT.0.0D+00) THEN
            ROTA    = REALN*TOR-GAMMA
            EN123R4 = EN123R4 + VR*(ONE + COS(ROTA))
            FACT    = VR*REALN*(-SIN(ROTA))
         END IF
         IF(ABS(VROT1)+ABS(VROT2)+ABS(VROT3).GT.0.0D+00) THEN
            COSTOR2 = COSTOR*COSTOR
            SINTOR  = SIN(TOR)
            EN123R4 = EN123R4+PT5*
     *                (VROT1*(ONE+COSTOR)+
     *                 VROT2*(TWO-TWO*COSTOR2)+
     *                 VROT3*(ONE+COSTOR*(FOUR*COSTOR2-THREE)))
            FACT    = PT5*(-SINTOR)*(+VROT1
     *                               -VROT2*FOUR*COSTOR
     *                               +VROT3*(12.0D+00*COSTOR2-THREE))
         END IF
C
C        - CALCULATING DERIVATIVES
C
         A123=1.0D+00/(R12*SIN2123)
         A432=1.0D+00/(R34*SIN2234)
         B123=R12*COS123/R23
         B432=R34*COS234/R23
         C123=B123-1.0D+00
C
         DUMY       =  FACT*A123/(R12*R23)
         DEX1       =  DUMY*(-Y12*Z23+Z12*Y23)
         DEY1       =  DUMY*(-Z12*X23+X12*Z23)
         DEZ1       =  DUMY*(-X12*Y23+Y12*X23)
         FFGRD(1,P1)=FFGRD(1,P1)+DEX1
         FFGRD(2,P1)=FFGRD(2,P1)+DEY1
         FFGRD(3,P1)=FFGRD(3,P1)+DEZ1
C
         DUMY       =  FACT*(-A432)/(R23*R34)
         DEX4       =  DUMY*(-Y23*Z34+Z23*Y34)
         DEY4       =  DUMY*(-Z23*X34+X23*Z34)
         DEZ4       =  DUMY*(-X23*Y34+Y23*X34)
         FFGRD(1,P4)=FFGRD(1,P4)+DEX4
         FFGRD(2,P4)=FFGRD(2,P4)+DEY4
         FFGRD(3,P4)=FFGRD(3,P4)+DEZ4
C
         DEX2       =  C123*DEX1 - B432*DEX4
         DEY2       =  C123*DEY1 - B432*DEY4
         DEZ2       =  C123*DEZ1 - B432*DEZ4
         FFGRD(1,P2)=FFGRD(1,P2)+DEX2
         FFGRD(2,P2)=FFGRD(2,P2)+DEY2
         FFGRD(3,P2)=FFGRD(3,P2)+DEZ2
C
         FFGRD(1,P3)=FFGRD(1,P3)-DEX1-DEX2-DEX4
         FFGRD(2,P3)=FFGRD(2,P3)-DEY1-DEY2-DEY4
         FFGRD(3,P3)=FFGRD(3,P3)-DEZ1-DEZ2-DEZ4
C
         VIR(1)     =VIR(1)  + DEX1*X13 + DEX2*X23 - DEX4*X34
         VIR(2)     =VIR(2)  + DEY1*Y13 + DEY2*Y23 - DEY4*Y34
         VIR(3)     =VIR(3)  + DEZ1*Z13 + DEZ2*Z23 - DEZ4*Z34
         IYES = 0
         JYES = 0
         KYES = 0
         LYES = 0
         DO KFIX=1,NFIXMM
            IF(P1.EQ.IFIXMM(KFIX)) IYES = 1
            IF(P2.EQ.IFIXMM(KFIX)) JYES = 1
            IF(P3.EQ.IFIXMM(KFIX)) KYES = 1
            IF(P4.EQ.IFIXMM(KFIX)) LYES = 1
         ENDDO
         IF(IYES.EQ.1.AND.KYES.EQ.1) THEN
            VIR(1)    =VIR(1)     - DEX1*X13
            VIR(2)    =VIR(2)     - DEY1*Y13
            VIR(3)    =VIR(3)     - DEZ1*Z13
         END IF
         IF(JYES.EQ.1.AND.KYES.EQ.1) THEN
            VIR(1)    =VIR(1)     - DEX2*X23
            VIR(2)    =VIR(2)     - DEY2*Y23
            VIR(3)    =VIR(3)     - DEZ2*Z23
         END IF
         IF(KYES.EQ.1.AND.LYES.EQ.1) THEN
            VIR(1)    =VIR(1)     + DEX4*X34
            VIR(2)    =VIR(2)     + DEY4*Y34
            VIR(3)    =VIR(3)     + DEZ4*Z34
         END IF
 100  CONTINUE
C
      DO 200 III=L1DIHRPMA,L2DIHRPMA
         IDIHR = LSDIHRPMA(III)
         IDIHR = LSDIHRPMA(III)
         VR      = VROT(IDIHR)
         REALN   = DBLE(NNN(IDIHR))
         GAMMA   = GAMA(IDIHR)
         VROT1   = FCDIHR(1,IDIHR)
         VROT2   = FCDIHR(2,IDIHR)
         VROT3   = FCDIHR(3,IDIHR)
         IF(ABS(VR)+ABS(VROT1)+ABS(VROT2)+ABS(VROT3).EQ.0.0D+00)
     *   GOTO 200
         P1=LLIST(1,IDIHR)
         P2=LLIST(2,IDIHR)
         P3=LLIST(3,IDIHR)
         P4=LLIST(4,IDIHR)
C
         X12=CORD(1,P1)-CORD(1,P2)
         Y12=CORD(2,P1)-CORD(2,P2)
         Z12=CORD(3,P1)-CORD(3,P2)
         X13=CORD(1,P1)-CORD(1,P3)
         Y13=CORD(2,P1)-CORD(2,P3)
         Z13=CORD(3,P1)-CORD(3,P3)
         X23=CORD(1,P2)-CORD(1,P3)
         Y23=CORD(2,P2)-CORD(2,P3)
         Z23=CORD(3,P2)-CORD(3,P3)
         X34=CORD(1,P3)-CORD(1,P4)
         Y34=CORD(2,P3)-CORD(2,P4)
         Z34=CORD(3,P3)-CORD(3,P4)
C
         R12=SQRT(X12*X12+Y12*Y12+Z12*Z12)
         R23=SQRT(X23*X23+Y23*Y23+Z23*Z23)
         R34=SQRT(X34*X34+Y34*Y34+Z34*Z34)
C
         COS123=(-(X12*X23)-(Y12*Y23)-(Z12*Z23))/(R12*R23)
         COS234=(-(X23*X34)-(Y23*Y34)-(Z23*Z34))/(R23*R34)
         SIN2123= 1.0D+00-COS123*COS123
         SIN2234= 1.0D+00-COS234*COS234
         SIN123 = SQRT(ABS(SIN2123))
         SIN234 = SQRT(ABS(SIN2234))
         IF(ABS(SIN123).LT.1.0D-06) GOTO 200
         IF(ABS(SIN234).LT.1.0D-06) GOTO 200
         ONESIN = 1.0D+00/(SIN123*SIN234)
C
         COSTOR = ONESIN*(COS123*COS234-
     *                    ((+X12*X34+Y12*Y34+Z12*Z34)/(R12*R34)))
         IF(COSTOR.GT. ONE) COSTOR= ONE
         IF(COSTOR.LT.-ONE) COSTOR=-ONE
         TOR    = ACOS(COSTOR)
C        -- DIHEDRAL ANGLE IS 0 - 360 DEGREES
         XNORM  = -Y23*Z34 + Z23*Y34
         YNORM  = -Z23*X34 + X23*Z34
         ZNORM  = -X23*Y34 + Y23*X34
         DOTN12 = X12*XNORM + Y12*YNORM + Z12*ZNORM
         IF(DOTN12.LT.0.0D+00) TOR = 2.0D+00*PI - TOR
C
         IF(ABS(VR).GT.0.0D+00) THEN
            ROTA    = REALN*TOR-GAMMA
            PMF1DR  = PMF1DR - VR*(ONE + COS(ROTA))
         END IF
         IF(ABS(VROT1)+ABS(VROT2)+ABS(VROT3).GT.0.0D+00) THEN
            COSTOR2 = COSTOR*COSTOR
            PMF1DR  = PMF1DR - PT5*
     *                (VROT1*(ONE+COSTOR)+
     *                 VROT2*(TWO-TWO*COSTOR2)+
     *                 VROT3*(ONE+COSTOR*(FOUR*COSTOR2-THREE)))
         END IF
 200  CONTINUE
C
      DO 210 III=L1DIHRPMB,L2DIHRPMB
         IDIHR = LSDIHRPMB(III)
         VR      = VROT(IDIHR)
         REALN   = DBLE(NNN(IDIHR))
         GAMMA   = GAMA(IDIHR)
         VROT1   = FCDIHR(1,IDIHR)
         VROT2   = FCDIHR(2,IDIHR)
         VROT3   = FCDIHR(3,IDIHR)
         IF(ABS(VR)+ABS(VROT1)+ABS(VROT2)+ABS(VROT3).EQ.0.0D+00)
     *   GOTO 210
         P1=LLIST(1,IDIHR)
         P2=LLIST(2,IDIHR)
         P3=LLIST(3,IDIHR)
         P4=LLIST(4,IDIHR)
C
         X12=CORD(1,P1)-CORD(1,P2)
         Y12=CORD(2,P1)-CORD(2,P2)
         Z12=CORD(3,P1)-CORD(3,P2)
         X13=CORD(1,P1)-CORD(1,P3)
         Y13=CORD(2,P1)-CORD(2,P3)
         Z13=CORD(3,P1)-CORD(3,P3)
         X23=CORD(1,P2)-CORD(1,P3)
         Y23=CORD(2,P2)-CORD(2,P3)
         Z23=CORD(3,P2)-CORD(3,P3)
         X34=CORD(1,P3)-CORD(1,P4)
         Y34=CORD(2,P3)-CORD(2,P4)
         Z34=CORD(3,P3)-CORD(3,P4)
C
         IF(IFEPTYP.EQ.2) THEN  !  NFIXMM=KFREEAB FOR IFEPTYP=2
            DO KFIX=1,NFIXMM
               IF(P1.EQ.IFIXMM(KFIX)) THEN
                  X12 = X12 - CORD(1,P1) + CORDB(1,P1)
                  Y12 = Y12 - CORD(2,P1) + CORDB(2,P1)
                  Z12 = Z12 - CORD(3,P1) + CORDB(3,P1)
                  X13 = X13 - CORD(1,P1) + CORDB(1,P1)
                  Y13 = Y13 - CORD(2,P1) + CORDB(2,P1)
                  Z13 = Z13 - CORD(3,P1) + CORDB(3,P1)
               END IF
               IF(P2.EQ.IFIXMM(KFIX)) THEN
                  X12 = X12 + CORD(1,P2) - CORDB(1,P2)
                  Y12 = Y12 + CORD(2,P2) - CORDB(2,P2)
                  Z12 = Z12 + CORD(3,P2) - CORDB(3,P2)
                  X23 = X23 - CORD(1,P2) + CORDB(1,P2)
                  Y23 = Y23 - CORD(2,P2) + CORDB(2,P2)
                  Z23 = Z23 - CORD(3,P2) + CORDB(3,P2)
               END IF
               IF(P3.EQ.IFIXMM(KFIX)) THEN
                  X13 = X13 + CORD(1,P3) - CORDB(1,P3)
                  Y13 = Y13 + CORD(2,P3) - CORDB(2,P3)
                  Z13 = Z13 + CORD(3,P3) - CORDB(3,P3)
                  X23 = X23 + CORD(1,P3) - CORDB(1,P3)
                  Y23 = Y23 + CORD(2,P3) - CORDB(2,P3)
                  Z23 = Z23 + CORD(3,P3) - CORDB(3,P3)
                  X34 = X34 - CORD(1,P3) + CORDB(1,P3)
                  Y34 = Y34 - CORD(2,P3) + CORDB(2,P3)
                  Z34 = Z34 - CORD(3,P3) + CORDB(3,P3)
               END IF
               IF(P4.EQ.IFIXMM(KFIX)) THEN
                  X34 = X34 + CORD(1,P4) - CORDB(1,P4)
                  Y34 = Y34 + CORD(2,P4) - CORDB(2,P4)
                  Z34 = Z34 + CORD(3,P4) - CORDB(3,P4)
               END IF
            ENDDO
         END IF
C
         R12=SQRT(X12*X12+Y12*Y12+Z12*Z12)
         R23=SQRT(X23*X23+Y23*Y23+Z23*Z23)
         R34=SQRT(X34*X34+Y34*Y34+Z34*Z34)
C
         COS123=(-(X12*X23)-(Y12*Y23)-(Z12*Z23))/(R12*R23)
         COS234=(-(X23*X34)-(Y23*Y34)-(Z23*Z34))/(R23*R34)
         SIN2123= 1.0D+00-COS123*COS123
         SIN2234= 1.0D+00-COS234*COS234
         SIN123 = SQRT(ABS(SIN2123))
         SIN234 = SQRT(ABS(SIN2234))
         IF(ABS(SIN123).LT.1.0D-06) GOTO 210
         IF(ABS(SIN234).LT.1.0D-06) GOTO 210
         ONESIN = 1.0D+00/(SIN123*SIN234)
C
         COSTOR = ONESIN*(COS123*COS234-
     *                    ((+X12*X34+Y12*Y34+Z12*Z34)/(R12*R34)))
         IF(COSTOR.GT. ONE) COSTOR= ONE
         IF(COSTOR.LT.-ONE) COSTOR=-ONE
         TOR    = ACOS(COSTOR)
C        -- DIHEDRAL ANGLE IS 0 - 360 DEGREES
         XNORM  = -Y23*Z34 + Z23*Y34
         YNORM  = -Z23*X34 + X23*Z34
         ZNORM  = -X23*Y34 + Y23*X34
         DOTN12 = X12*XNORM + Y12*YNORM + Z12*ZNORM
         IF(DOTN12.LT.0.0D+00) TOR = 2.0D+00*PI - TOR
C
         IF(ABS(VR).GT.0.0D+00) THEN
            ROTA    = REALN*TOR-GAMMA
            PMF1DR  = PMF1DR + VR*(ONE + COS(ROTA))
         END IF
         IF(ABS(VROT1)+ABS(VROT2)+ABS(VROT3).GT.0.0D+00) THEN
            COSTOR2 = COSTOR*COSTOR
            PMF1DR  = PMF1DR + PT5*
     *                (VROT1*(ONE+COSTOR)+
     *                 VROT2*(TWO-TWO*COSTOR2)+
     *                 VROT3*(ONE+COSTOR*(FOUR*COSTOR2-THREE)))
         END IF
 210  CONTINUE
C
      RETURN
      END
C*MODULE QUANPOB  *DECK E123B4
!>
!> @brief    dihedral bending
!>
!> @author   Nandun Thellamurege
!>           - Feb 2011
!>
!> @details  force field dihedral bending energy
!>           (out-of-plane energy)
!>
      SUBROUTINE E123B4(CORD,FFGRD,DIHB0,FCDIHB,NLIST,CORDB,
     *                  LSDIHBPMA,LSDIHBPMB)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      INTEGER P1, P2, P3, P4
C
      PARAMETER (PI=3.14159265358979323846264338D+00)
      PARAMETER (TWO=2.0D+00)
      PARAMETER (ONE=1.0D+00)
C
      DIMENSION CORD(3,*),FFGRD(3,*),DIHB0(*),FCDIHB(*),NLIST(4,*),
     *          CORDB(3,*),LSDIHBPMA(*),LSDIHBPMB(*)
C
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
      COMMON /FFMDPA/ DT,DT2,TEMP0,PRES0,POLTOL,VIR(3),PMEAN,VOLAV,
     *                ENPAV,ENKAV,TEMPAV,BERENDT,BERENDP,VELMAX,
     *                PMEANX,PMEANY,PMEANZ,NSTEP,KMASTER,KOUTACT(2),
     *                IHESS,INTALG,ITSTAT,IPSTAT,JOUT,KOUT,LOUT
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
      COMMON /FFRMSD/ DIPT(8),TIMGYRA,TIMRALL,
     *                NATPDB,NGYRA,NDIEL,NRALL,
     *                LFFRALL0,NRMSD,LFFRMSD0,KFREEAB(201),
     *                NRIJMM,IJRMM(2,100),NRIJQM,IJRQM(2,100),
     *                NAIJKMM,IJKMM(3,100),NAIJKQM,IJKQM(3,100),
     *                NFIXMM,IFIXMM(200),NFIXQM,IFIXQM(200)
C
C     NANDUN THELLAMUREGE, FEB 2011, LINCOLN
C     FENGCHAO CUI, HUI LI, MAY 2012
C
C     FORMULA:  E = K*(B-B0)**2
C
      PMF1DB =0.0D+00
      EN123B4=0.0D+00
      DO 100 I=L1DIHB,L2DIHB
         IF(FCDIHB(I).EQ.0.0D+00) GOTO 100
         FC=FCDIHB(I)
         A0=DIHB0(I)
         P1=NLIST(1,I)
         P2=NLIST(2,I)
         P3=NLIST(3,I)
         P4=NLIST(4,I)
C
         X12=CORD(1,P1)-CORD(1,P2)
         Y12=CORD(2,P1)-CORD(2,P2)
         Z12=CORD(3,P1)-CORD(3,P2)
         X13=CORD(1,P1)-CORD(1,P3)
         Y13=CORD(2,P1)-CORD(2,P3)
         Z13=CORD(3,P1)-CORD(3,P3)
         X23=CORD(1,P2)-CORD(1,P3)
         Y23=CORD(2,P2)-CORD(2,P3)
         Z23=CORD(3,P2)-CORD(3,P3)
         X34=CORD(1,P3)-CORD(1,P4)
         Y34=CORD(2,P3)-CORD(2,P4)
         Z34=CORD(3,P3)-CORD(3,P4)
C
         R12=SQRT(X12*X12+Y12*Y12+Z12*Z12)
         R23=SQRT(X23*X23+Y23*Y23+Z23*Z23)
         R34=SQRT(X34*X34+Y34*Y34+Z34*Z34)
C
         COS123=(-(X12*X23)-(Y12*Y23)-(Z12*Z23))/(R12*R23)
         COS234=(-(X23*X34)-(Y23*Y34)-(Z23*Z34))/(R23*R34)
         SIN2123= 1.0D+00-COS123*COS123
         SIN2234= 1.0D+00-COS234*COS234
         SIN123 = SQRT(ABS(SIN2123))
         SIN234 = SQRT(ABS(SIN2234))
         IF(ABS(SIN123).LT.1.0D-06) GOTO 100
         IF(ABS(SIN234).LT.1.0D-06) GOTO 100
         ONESIN = 1.0D+00/(SIN123*SIN234)
C
         COSTOR = ONESIN*(COS123*COS234-
     *                    ((+X12*X34+Y12*Y34+Z12*Z34)/(R12*R34)))
         IF(COSTOR.GT. ONE) COSTOR= ONE
         IF(COSTOR.LT.-ONE) COSTOR=-ONE
         TOR    = ACOS(COSTOR)
         XNORM  = -Y23*Z34 + Z23*Y34
         YNORM  = -Z23*X34 + X23*Z34
         ZNORM  = -X23*Y34 + Y23*X34
         DOTN12 = X12*XNORM + Y12*YNORM + Z12*ZNORM
C
C        -- TOR FROM ACOS() IS ALWAYS 0 - 180 DEGREE
C           IF A0 =   0, TOR SHOULD BE  -20 TO  +20 DEGREE
C           IF A0 = 180, TOR SHOULD BE +160 TO +200 DEGREE
C
         IF(ABS(A0).LE.0.2D+00.AND.DOTN12.LT.0.0D+00)
     *      TOR= -TOR
         IF(ABS(A0-PI).LE.0.2D+00.AND.DOTN12.LT.0.0D+00)
     *      TOR = 2.0D+00*PI - TOR
C
         EN123B4 = EN123B4 + FC*(TOR-A0)*(TOR-A0)
C
C        -- DERIVATIVES --
C
         A123=1.0D+00/(R12*SIN2123)
         A432=1.0D+00/(R34*SIN2234)
         B123=R12*COS123/R23
         B432=R34*COS234/R23
         C123=B123-1.0D+00
         FACT=TWO*FC*(TOR-A0)
C
         DUMY       =  FACT*A123/(R12*R23)
         DEX1       =  DUMY*(-Y12*Z23+Z12*Y23)
         DEY1       =  DUMY*(-Z12*X23+X12*Z23)
         DEZ1       =  DUMY*(-X12*Y23+Y12*X23)
         FFGRD(1,P1)=FFGRD(1,P1)+DEX1
         FFGRD(2,P1)=FFGRD(2,P1)+DEY1
         FFGRD(3,P1)=FFGRD(3,P1)+DEZ1
C
         DUMY       =  FACT*(-A432)/(R23*R34)
         DEX4       =  DUMY*(-Y23*Z34+Z23*Y34)
         DEY4       =  DUMY*(-Z23*X34+X23*Z34)
         DEZ4       =  DUMY*(-X23*Y34+Y23*X34)
         FFGRD(1,P4)=FFGRD(1,P4)+DEX4
         FFGRD(2,P4)=FFGRD(2,P4)+DEY4
         FFGRD(3,P4)=FFGRD(3,P4)+DEZ4
C
         DEX2       =  C123*DEX1 - B432*DEX4
         DEY2       =  C123*DEY1 - B432*DEY4
         DEZ2       =  C123*DEZ1 - B432*DEZ4
         FFGRD(1,P2)=FFGRD(1,P2)+DEX2
         FFGRD(2,P2)=FFGRD(2,P2)+DEY2
         FFGRD(3,P2)=FFGRD(3,P2)+DEZ2
C
         FFGRD(1,P3)=FFGRD(1,P3)-DEX1-DEX2-DEX4
         FFGRD(2,P3)=FFGRD(2,P3)-DEY1-DEY2-DEY4
         FFGRD(3,P3)=FFGRD(3,P3)-DEZ1-DEZ2-DEZ4
C
         VIR(1)     =VIR(1)  + DEX1*X13 + DEX2*X23 - DEX4*X34
         VIR(2)     =VIR(2)  + DEY1*Y13 + DEY2*Y23 - DEY4*Y34
         VIR(3)     =VIR(3)  + DEZ1*Z13 + DEZ2*Z23 - DEZ4*Z34
         IYES = 0
         JYES = 0
         KYES = 0
         LYES = 0
         DO KFIX=1,NFIXMM
            IF(P1.EQ.IFIXMM(KFIX)) IYES = 1
            IF(P2.EQ.IFIXMM(KFIX)) JYES = 1
            IF(P3.EQ.IFIXMM(KFIX)) KYES = 1
            IF(P4.EQ.IFIXMM(KFIX)) LYES = 1
         ENDDO
         IF(IYES.EQ.1.AND.KYES.EQ.1) THEN
            VIR(1)    =VIR(1)     - DEX1*X13
            VIR(2)    =VIR(2)     - DEY1*Y13
            VIR(3)    =VIR(3)     - DEZ1*Z13
         END IF
         IF(JYES.EQ.1.AND.KYES.EQ.1) THEN
            VIR(1)    =VIR(1)     - DEX2*X23
            VIR(2)    =VIR(2)     - DEY2*Y23
            VIR(3)    =VIR(3)     - DEZ2*Z23
         END IF
         IF(KYES.EQ.1.AND.LYES.EQ.1) THEN
            VIR(1)    =VIR(1)     + DEX4*X34
            VIR(2)    =VIR(2)     + DEY4*Y34
            VIR(3)    =VIR(3)     + DEZ4*Z34
         END IF
 100  CONTINUE
C
C
      DO 200 III=L1DIHBPMA,L2DIHBPMA
         IDIHB = LSDIHBPMA(III)
         IF(FCDIHB(IDIHB).EQ.0.0D+00) GOTO 200
         FC=FCDIHB(IDIHB)
         A0=DIHB0(IDIHB)
         P1=NLIST(1,IDIHB)
         P2=NLIST(2,IDIHB)
         P3=NLIST(3,IDIHB)
         P4=NLIST(4,IDIHB)
C
         X12=CORD(1,P1)-CORD(1,P2)
         Y12=CORD(2,P1)-CORD(2,P2)
         Z12=CORD(3,P1)-CORD(3,P2)
         X13=CORD(1,P1)-CORD(1,P3)
         Y13=CORD(2,P1)-CORD(2,P3)
         Z13=CORD(3,P1)-CORD(3,P3)
         X23=CORD(1,P2)-CORD(1,P3)
         Y23=CORD(2,P2)-CORD(2,P3)
         Z23=CORD(3,P2)-CORD(3,P3)
         X34=CORD(1,P3)-CORD(1,P4)
         Y34=CORD(2,P3)-CORD(2,P4)
         Z34=CORD(3,P3)-CORD(3,P4)
C
         R12=SQRT(X12*X12+Y12*Y12+Z12*Z12)
         R23=SQRT(X23*X23+Y23*Y23+Z23*Z23)
         R34=SQRT(X34*X34+Y34*Y34+Z34*Z34)
C
         COS123=(-(X12*X23)-(Y12*Y23)-(Z12*Z23))/(R12*R23)
         COS234=(-(X23*X34)-(Y23*Y34)-(Z23*Z34))/(R23*R34)
         SIN2123= 1.0D+00-COS123*COS123
         SIN2234= 1.0D+00-COS234*COS234
         SIN123 = SQRT(ABS(SIN2123))
         SIN234 = SQRT(ABS(SIN2234))
         IF(ABS(SIN123).LT.1.0D-06) GOTO 200
         IF(ABS(SIN234).LT.1.0D-06) GOTO 200
         ONESIN = 1.0D+00/(SIN123*SIN234)
C
         COSTOR = ONESIN*(COS123*COS234-
     *                    ((+X12*X34+Y12*Y34+Z12*Z34)/(R12*R34)))
         IF(COSTOR.GT. ONE) COSTOR= ONE
         IF(COSTOR.LT.-ONE) COSTOR=-ONE
         TOR    = ACOS(COSTOR)
         XNORM  = -Y23*Z34 + Z23*Y34
         YNORM  = -Z23*X34 + X23*Z34
         ZNORM  = -X23*Y34 + Y23*X34
         DOTN12 = X12*XNORM + Y12*YNORM + Z12*ZNORM
C
C        -- TOR FROM ACOS() IS ALWAYS 0 - 180 DEGREE
C           IF A0 =   0, TOR SHOULD BE  -20 TO  +20 DEGREE
C           IF A0 = 180, TOR SHOULD BE +160 TO +200 DEGREE
C
         IF(ABS(A0).LE.0.2D+00.AND.DOTN12.LT.0.0D+00)
     *      TOR  = -TOR
         IF(ABS(A0-PI).LE.0.2D+00.AND.DOTN12.LT.0.0D+00)
     *      TOR  = 2.0D+00*PI - TOR
C
         PMF1DB  = PMF1DB - FC*(TOR-A0)*(TOR-A0)
 200  CONTINUE
C
      DO 210 III=L1DIHBPMB,L2DIHBPMB
         IDIHB = LSDIHBPMB(III)
         IF(FCDIHB(IDIHB).EQ.0.0D+00) GOTO 210
         FC=FCDIHB(IDIHB)
         A0=DIHB0(IDIHB)
         P1=NLIST(1,IDIHB)
         P2=NLIST(2,IDIHB)
         P3=NLIST(3,IDIHB)
         P4=NLIST(4,IDIHB)
C
         X12=CORD(1,P1)-CORD(1,P2)
         Y12=CORD(2,P1)-CORD(2,P2)
         Z12=CORD(3,P1)-CORD(3,P2)
         X13=CORD(1,P1)-CORD(1,P3)
         Y13=CORD(2,P1)-CORD(2,P3)
         Z13=CORD(3,P1)-CORD(3,P3)
         X23=CORD(1,P2)-CORD(1,P3)
         Y23=CORD(2,P2)-CORD(2,P3)
         Z23=CORD(3,P2)-CORD(3,P3)
         X34=CORD(1,P3)-CORD(1,P4)
         Y34=CORD(2,P3)-CORD(2,P4)
         Z34=CORD(3,P3)-CORD(3,P4)
C
         IF(IFEPTYP.EQ.2) THEN  !  NFIXMM=KFREEAB FOR IFEPTYP=2
            DO KFIX=1,NFIXMM
               IF(P1.EQ.IFIXMM(KFIX)) THEN
                  X12 = X12 - CORD(1,P1) + CORDB(1,P1)
                  Y12 = Y12 - CORD(2,P1) + CORDB(2,P1)
                  Z12 = Z12 - CORD(3,P1) + CORDB(3,P1)
                  X13 = X13 - CORD(1,P1) + CORDB(1,P1)
                  Y13 = Y13 - CORD(2,P1) + CORDB(2,P1)
                  Z13 = Z13 - CORD(3,P1) + CORDB(3,P1)
               END IF
               IF(P2.EQ.IFIXMM(KFIX)) THEN
                  X12 = X12 + CORD(1,P2) - CORDB(1,P2)
                  Y12 = Y12 + CORD(2,P2) - CORDB(2,P2)
                  Z12 = Z12 + CORD(3,P2) - CORDB(3,P2)
                  X23 = X23 - CORD(1,P2) + CORDB(1,P2)
                  Y23 = Y23 - CORD(2,P2) + CORDB(2,P2)
                  Z23 = Z23 - CORD(3,P2) + CORDB(3,P2)
               END IF
               IF(P3.EQ.IFIXMM(KFIX)) THEN
                  X13 = X13 + CORD(1,P3) - CORDB(1,P3)
                  Y13 = Y13 + CORD(2,P3) - CORDB(2,P3)
                  Z13 = Z13 + CORD(3,P3) - CORDB(3,P3)
                  X23 = X23 + CORD(1,P3) - CORDB(1,P3)
                  Y23 = Y23 + CORD(2,P3) - CORDB(2,P3)
                  Z23 = Z23 + CORD(3,P3) - CORDB(3,P3)
                  X34 = X34 - CORD(1,P3) + CORDB(1,P3)
                  Y34 = Y34 - CORD(2,P3) + CORDB(2,P3)
                  Z34 = Z34 - CORD(3,P3) + CORDB(3,P3)
               END IF
               IF(P4.EQ.IFIXMM(KFIX)) THEN
                  X34 = X34 + CORD(1,P4) - CORDB(1,P4)
                  Y34 = Y34 + CORD(2,P4) - CORDB(2,P4)
                  Z34 = Z34 + CORD(3,P4) - CORDB(3,P4)
               END IF
            ENDDO
         END IF
C
         R12=SQRT(X12*X12+Y12*Y12+Z12*Z12)
         R23=SQRT(X23*X23+Y23*Y23+Z23*Z23)
         R34=SQRT(X34*X34+Y34*Y34+Z34*Z34)
C
         COS123=(-(X12*X23)-(Y12*Y23)-(Z12*Z23))/(R12*R23)
         COS234=(-(X23*X34)-(Y23*Y34)-(Z23*Z34))/(R23*R34)
         SIN2123= 1.0D+00-COS123*COS123
         SIN2234= 1.0D+00-COS234*COS234
         SIN123 = SQRT(ABS(SIN2123))
         SIN234 = SQRT(ABS(SIN2234))
         IF(ABS(SIN123).LT.1.0D-06) GOTO 210
         IF(ABS(SIN234).LT.1.0D-06) GOTO 210
         ONESIN = 1.0D+00/(SIN123*SIN234)
C
         COSTOR = ONESIN*(COS123*COS234-
     *                    ((+X12*X34+Y12*Y34+Z12*Z34)/(R12*R34)))
         IF(COSTOR.GT. ONE) COSTOR= ONE
         IF(COSTOR.LT.-ONE) COSTOR=-ONE
         TOR    = ACOS(COSTOR)
         XNORM  = -Y23*Z34 + Z23*Y34
         YNORM  = -Z23*X34 + X23*Z34
         ZNORM  = -X23*Y34 + Y23*X34
         DOTN12 = X12*XNORM + Y12*YNORM + Z12*ZNORM
C
C        -- TOR FROM ACOS() IS ALWAYS 0 - 180 DEGREE
C           IF A0 =   0, TOR SHOULD BE  -20 TO  +20 DEGREE
C           IF A0 = 180, TOR SHOULD BE +160 TO +200 DEGREE
C
         IF(ABS(A0).LE.0.2D+00.AND.DOTN12.LT.0.0D+00)
     *      TOR  = -TOR
         IF(ABS(A0-PI).LE.0.2D+00.AND.DOTN12.LT.0.0D+00)
     *      TOR  = 2.0D+00*PI - TOR
C
         PMF1DB  = PMF1DB + FC*(TOR-A0)*(TOR-A0)
 210  CONTINUE
C
      RETURN
      END
C*MODULE QUANPOB  *DECK E234W1
!>
!> @brief    wagging energy
!>
!> @author   Nandun Thellamurege
!>           - Jan 2011
!>
!> @details  force field wagging energy
!>
      SUBROUTINE E234W1(CORD,FFGRD,FCWAGG,MLIST,CORDB,LSWAGGPMA,
     *                  LSWAGGPMB)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      INTEGER P1, P2, P3, P4
C
      DIMENSION CORD(3,*),FFGRD(3,*),MLIST(4,*),FCWAGG(*),
     *          CORDB(3,*),LSWAGGPMA(*),LSWAGGPMB(*)
C
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
      COMMON /FFMDPA/ DT,DT2,TEMP0,PRES0,POLTOL,VIR(3),PMEAN,VOLAV,
     *                ENPAV,ENKAV,TEMPAV,BERENDT,BERENDP,VELMAX,
     *                PMEANX,PMEANY,PMEANZ,NSTEP,KMASTER,KOUTACT(2),
     *                IHESS,INTALG,ITSTAT,IPSTAT,JOUT,KOUT,LOUT
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
      COMMON /FFRMSD/ DIPT(8),TIMGYRA,TIMRALL,
     *                NATPDB,NGYRA,NDIEL,NRALL,
     *                LFFRALL0,NRMSD,LFFRMSD0,KFREEAB(201),
     *                NRIJMM,IJRMM(2,100),NRIJQM,IJRQM(2,100),
     *                NAIJKMM,IJKMM(3,100),NAIJKQM,IJKQM(3,100),
     *                NFIXMM,IFIXMM(200),NFIXQM,IFIXQM(200)
C
C     NANDUN THELLAMUREGE, JAN 2011, LINCOLN
C     FENGCHAO CUI, HUI LI, MAY 2012
C
C     FORMULA:  E = K*W**2
C
      PMF1WG  = 0.0D+00
      EN234W1 = 0.0D+00
      DO 100 IWAGG = L1WAGG, L2WAGG
         FCONST = FCWAGG(IWAGG)
         IF(FCONST.EQ.0.0D+00) GOTO 100
C
C        --- P1 IS THE WAGGING ATOM, P4 IS THE CENTER LINKING ATOM ---
C
         P1     = MLIST(1,IWAGG)
         P2     = MLIST(2,IWAGG)
         P3     = MLIST(3,IWAGG)
         P4     = MLIST(4,IWAGG)
C
         X41    = CORD(1,P4)-CORD(1,P1)
         Y41    = CORD(2,P4)-CORD(2,P1)
         Z41    = CORD(3,P4)-CORD(3,P1)
         X42    = CORD(1,P4)-CORD(1,P2)
         Y42    = CORD(2,P4)-CORD(2,P2)
         Z42    = CORD(3,P4)-CORD(3,P2)
         X43    = CORD(1,P4)-CORD(1,P3)
         Y43    = CORD(2,P4)-CORD(2,P3)
         Z43    = CORD(3,P4)-CORD(3,P3)
         X23    = CORD(1,P2)-CORD(1,P3)
         Y23    = CORD(2,P2)-CORD(2,P3)
         Z23    = CORD(3,P2)-CORD(3,P3)
C
         R41SQ  = X41*X41 + Y41*Y41 + Z41*Z41
         R41    = SQRT(R41SQ)
         ONER41 = 1.0D+00/R41
         R42SQ  = X42*X42 +Y42*Y42 + Z42*Z42
         R42    = SQRT(R42SQ)
         ONER42 = 1.0D+00/R42
         R43SQ  = X43*X43 + Y43*Y43 + Z43*Z43
         R43    = SQRT(R43SQ)
         ONER43 = 1.0D+00/R43
         R23SQ  = X23*X23 + Y23*Y23 + Z23*Z23
         R23    = SQRT(R23SQ)
         P4243X = Y42*Z43 - Z42*Y43
         P4243Y = Z42*X43 - X42*Z43
         P4243Z = X42*Y43 - Y42*X43
C
C        --- ANGLE 243 CAN BE 90 -120 ---
         COS243 = (R42*R42 + R43*R43 - R23*R23)/(2.0D+00*R42*R43)
         SIN243 = SQRT(ABS(1.0D+00 - COS243*COS243))
         DUM    = 1.0D+00/(SIN243*R42*R43)
         AX     = P4243X*DUM
         AY     = P4243Y*DUM
         AZ     = P4243Z*DUM
C        --- WWW IS LIKELY -30 TO +30 ---
         SINW   = -(AX*X41 + AY*Y41 + AZ*Z41)*ONER41
         COSW   = SQRT(ABS(1.0D+00 - SINW*SINW))
         ONECOSW= 1.0D+00/COSW
         WWW    = ASIN(SINW)
         EN234W1= EN234W1 + FCONST*WWW*WWW
C
         C243   = 1.0D+00/(R42*SIN243)
         C342   = 1.0D+00/(R43*SIN243)
         B2     = R42*SINW
         B3     = R43*SINW
C
C        --- CALCULATING DERIVATIVES ---
C
         FACT   = 2.0D+00*FCONST*WWW*ONECOSW
C
         DEX1   = FACT*ONER41*(AX - SINW*(-X41)*ONER41)
         DEY1   = FACT*ONER41*(AY - SINW*(-Y41)*ONER41)
         DEZ1   = FACT*ONER41*(AZ - SINW*(-Z41)*ONER41)
C
         D2X    = -B3*ONER43*(X43*COS243*ONER43 - X42*ONER42)/SIN243
         D2Y    = -B3*ONER43*(Y43*COS243*ONER43 - Y42*ONER42)/SIN243
         D2Z    = -B3*ONER43*(Z43*COS243*ONER43 - Z42*ONER42)/SIN243
         P4341X = (Y43*Z41 - Z43*Y41)/(R43*R41)
         P4341Y = (Z43*X41 - X43*Z41)/(R43*R41)
         P4341Z = (X43*Y41 - Y43*X41)/(R43*R41)
         DEX2   = FACT*C243*(D2X + P4341X)
         DEY2   = FACT*C243*(D2Y + P4341Y)
         DEZ2   = FACT*C243*(D2Z + P4341Z)
C
         D3X    = -B2*ONER42*(X42*COS243*ONER42 - X43*ONER43)/SIN243
         D3Y    = -B2*ONER42*(Y42*COS243*ONER42 - Y43*ONER43)/SIN243
         D3Z    = -B2*ONER42*(Z42*COS243*ONER42 - Z43*ONER43)/SIN243
         P4142X = (Y41*Z42 - Z41*Y42)/(R41*R42)
         P4142Y = (Z41*X42 - X41*Z42)/(R41*R42)
         P4142Z = (X41*Y42 - Y41*X42)/(R41*R42)
         DEX3   = FACT*C342*(D3X + P4142X)
         DEY3   = FACT*C342*(D3Y + P4142Y)
         DEZ3   = FACT*C342*(D3Z + P4142Z)
C
         FFGRD(1,P1)=FFGRD(1,P1)+DEX1
         FFGRD(2,P1)=FFGRD(2,P1)+DEY1
         FFGRD(3,P1)=FFGRD(3,P1)+DEZ1
         FFGRD(1,P2)=FFGRD(1,P2)+DEX2
         FFGRD(2,P2)=FFGRD(2,P2)+DEY2
         FFGRD(3,P2)=FFGRD(3,P2)+DEZ2
         FFGRD(1,P3)=FFGRD(1,P3)+DEX3
         FFGRD(2,P3)=FFGRD(2,P3)+DEY3
         FFGRD(3,P3)=FFGRD(3,P3)+DEZ3
         FFGRD(1,P4)=FFGRD(1,P4)-DEX1-DEX2-DEX3
         FFGRD(2,P4)=FFGRD(2,P4)-DEY1-DEY2-DEY3
         FFGRD(3,P4)=FFGRD(3,P4)-DEZ1-DEZ2-DEZ3
C
         VIR(1)     =VIR(1)  - DEX1*X41 - DEX2*X42 - DEX3*X43
         VIR(2)     =VIR(2)  - DEY1*Y41 - DEY2*Y42 - DEY3*Y43
         VIR(3)     =VIR(3)  - DEZ1*Z41 - DEZ2*Z42 - DEZ3*Z43
         IYES = 0
         JYES = 0
         KYES = 0
         LYES = 0
         DO KFIX=1,NFIXMM
            IF(P1.EQ.IFIXMM(KFIX)) IYES = 1
            IF(P2.EQ.IFIXMM(KFIX)) JYES = 1
            IF(P3.EQ.IFIXMM(KFIX)) KYES = 1
            IF(P4.EQ.IFIXMM(KFIX)) LYES = 1
         ENDDO
         IF(IYES.EQ.1.AND.LYES.EQ.1) THEN
            VIR(1)    =VIR(1)     + DEX1*X41
            VIR(2)    =VIR(2)     + DEY1*Y41
            VIR(3)    =VIR(3)     + DEZ1*Z41
         END IF
         IF(JYES.EQ.1.AND.LYES.EQ.1) THEN
            VIR(1)    =VIR(1)     + DEX2*X42
            VIR(2)    =VIR(2)     + DEY2*Y42
            VIR(3)    =VIR(3)     + DEZ2*Z42
         END IF
         IF(KYES.EQ.1.AND.LYES.EQ.1) THEN
            VIR(1)    =VIR(1)     + DEX3*X43
            VIR(2)    =VIR(2)     + DEY3*Y43
            VIR(3)    =VIR(3)     + DEZ3*Z43
         END IF
 100  CONTINUE
C
      DO 200 III=L1WAGGPMA,L2WAGGPMA
         IWAGG  = LSWAGGPMA(III)
         FCONST = FCWAGG(IWAGG)
         IF(FCONST.EQ.0.0D+00) GOTO 200
C
C        --- P1 IS THE WAGGING ATOM, P4 IS THE LINKING ATOM ---
C
         P1     = MLIST(1,IWAGG)
         P2     = MLIST(2,IWAGG)
         P3     = MLIST(3,IWAGG)
         P4     = MLIST(4,IWAGG)
C
         X41    = CORD(1,P4)-CORD(1,P1)
         Y41    = CORD(2,P4)-CORD(2,P1)
         Z41    = CORD(3,P4)-CORD(3,P1)
         X42    = CORD(1,P4)-CORD(1,P2)
         Y42    = CORD(2,P4)-CORD(2,P2)
         Z42    = CORD(3,P4)-CORD(3,P2)
         X43    = CORD(1,P4)-CORD(1,P3)
         Y43    = CORD(2,P4)-CORD(2,P3)
         Z43    = CORD(3,P4)-CORD(3,P3)
         X23    = CORD(1,P2)-CORD(1,P3)
         Y23    = CORD(2,P2)-CORD(2,P3)
         Z23    = CORD(3,P2)-CORD(3,P3)
C
         R41SQ  = X41*X41 + Y41*Y41 + Z41*Z41
         R41    = SQRT(R41SQ)
         ONER41 = 1.0D+00/R41
         R42SQ  = X42*X42 +Y42*Y42 + Z42*Z42
         R42    = SQRT(R42SQ)
         ONER42 = 1.0D+00/R42
         R43SQ  = X43*X43 + Y43*Y43 + Z43*Z43
         R43    = SQRT(R43SQ)
         ONER43 = 1.0D+00/R43
         R23SQ  = X23*X23 + Y23*Y23 + Z23*Z23
         R23    = SQRT(R23SQ)
         P4243X = Y42*Z43 - Z42*Y43
         P4243Y = Z42*X43 - X42*Z43
         P4243Z = X42*Y43 - Y42*X43
C
C        --- ANGLE 243 CAN BE 90 -120 ---
         COS243 = (R42*R42 + R43*R43 - R23*R23)/(2.0D+00*R42*R43)
         SIN243 = SQRT(ABS(1.0D+00 - COS243*COS243))
         DUM    = 1.0D+00/(SIN243*R42*R43)
         AX     = P4243X*DUM
         AY     = P4243Y*DUM
         AZ     = P4243Z*DUM
C        --- WWW IS LIKELY -30 TO +30 ---
         SINW   = -(AX*X41 + AY*Y41 + AZ*Z41)*ONER41
         COSW   = SQRT(ABS(1.0D+00 - SINW*SINW))
         ONECOSW= 1.0D+00/COSW
         WWW     = ASIN(SINW)
C
         PMF1WG = PMF1WG - FCONST*WWW*WWW
 200  CONTINUE
C
      DO 210 III=L1WAGGPMB,L2WAGGPMB
         IWAGG  = LSWAGGPMB(III)
         FCONST = FCWAGG(IWAGG)
         IF(FCONST.EQ.0.0D+00) GOTO 210
C
C        --- P1 IS THE WAGGING ATOM, P4 IS THE LINKING ATOM ---
C
         P1     = MLIST(1,IWAGG)
         P2     = MLIST(2,IWAGG)
         P3     = MLIST(3,IWAGG)
         P4     = MLIST(4,IWAGG)
C
         X41    = CORD(1,P4)-CORD(1,P1)
         Y41    = CORD(2,P4)-CORD(2,P1)
         Z41    = CORD(3,P4)-CORD(3,P1)
         X42    = CORD(1,P4)-CORD(1,P2)
         Y42    = CORD(2,P4)-CORD(2,P2)
         Z42    = CORD(3,P4)-CORD(3,P2)
         X43    = CORD(1,P4)-CORD(1,P3)
         Y43    = CORD(2,P4)-CORD(2,P3)
         Z43    = CORD(3,P4)-CORD(3,P3)
         X23    = CORD(1,P2)-CORD(1,P3)
         Y23    = CORD(2,P2)-CORD(2,P3)
         Z23    = CORD(3,P2)-CORD(3,P3)
C
         IF(IFEPTYP.EQ.2) THEN   !  NFIXMM=KFREEAB FOR IFEPTYP=2
            DO KFIX=1,NFIXMM
               IF(P1.EQ.IFIXMM(KFIX)) THEN
                  X41 = X41 + CORD(1,P1) - CORDB(1,P1)
                  Y41 = Y41 + CORD(2,P1) - CORDB(2,P1)
                  Z41 = Z41 + CORD(3,P1) - CORDB(3,P1)
               END IF
               IF(P2.EQ.IFIXMM(KFIX)) THEN
                  X42 = X42 + CORD(1,P2) - CORDB(1,P2)
                  Y42 = Y42 + CORD(2,P2) - CORDB(2,P2)
                  Z42 = Z42 + CORD(3,P2) - CORDB(3,P2)
                  X23 = X23 - CORD(1,P2) + CORDB(1,P2)
                  Y23 = Y23 - CORD(2,P2) + CORDB(2,P2)
                  Z23 = Z23 - CORD(3,P2) + CORDB(3,P2)
               END IF
               IF(P3.EQ.IFIXMM(KFIX)) THEN
                  X43 = X43 + CORD(1,P3) - CORDB(1,P3)
                  Y43 = Y43 + CORD(2,P3) - CORDB(2,P3)
                  Z43 = Z43 + CORD(3,P3) - CORDB(3,P3)
                  X23 = X23 + CORD(1,P3) - CORDB(1,P3)
                  Y23 = Y23 + CORD(2,P3) - CORDB(2,P3)
                  Z23 = Z23 + CORD(3,P3) - CORDB(3,P3)
               END IF
               IF(P4.EQ.IFIXMM(KFIX)) THEN
                  X41 = X41 - CORD(1,P4) + CORDB(1,P4)
                  Y41 = Y41 - CORD(2,P4) + CORDB(2,P4)
                  Z41 = Z41 - CORD(3,P4) + CORDB(3,P4)
                  X42 = X42 - CORD(1,P4) + CORDB(1,P4)
                  Y42 = Y42 - CORD(2,P4) + CORDB(2,P4)
                  Z42 = Z42 - CORD(3,P4) + CORDB(3,P4)
                  X43 = X43 - CORD(1,P4) + CORDB(1,P4)
                  Y43 = Y43 - CORD(2,P4) + CORDB(2,P4)
                  Z43 = Z43 - CORD(3,P4) + CORDB(3,P4)
               END IF
            ENDDO
         END IF
C
         R41SQ  = X41*X41 + Y41*Y41 + Z41*Z41
         R41    = SQRT(R41SQ)
         ONER41 = 1.0D+00/R41
         R42SQ  = X42*X42 +Y42*Y42 + Z42*Z42
         R42    = SQRT(R42SQ)
         ONER42 = 1.0D+00/R42
         R43SQ  = X43*X43 + Y43*Y43 + Z43*Z43
         R43    = SQRT(R43SQ)
         ONER43 = 1.0D+00/R43
         R23SQ  = X23*X23 + Y23*Y23 + Z23*Z23
         R23    = SQRT(R23SQ)
         P4243X = Y42*Z43 - Z42*Y43
         P4243Y = Z42*X43 - X42*Z43
         P4243Z = X42*Y43 - Y42*X43
C
C        --- ANGLE 243 CAN BE 90 -120 ---
         COS243 = (R42*R42 + R43*R43 - R23*R23)/(2.0D+00*R42*R43)
         SIN243 = SQRT(ABS(1.0D+00 - COS243*COS243))
         DUM    = 1.0D+00/(SIN243*R42*R43)
         AX     = P4243X*DUM
         AY     = P4243Y*DUM
         AZ     = P4243Z*DUM
C        --- WWW IS LIKELY -30 TO +30 ---
         SINW   = -(AX*X41 + AY*Y41 + AZ*Z41)*ONER41
         COSW   = SQRT(ABS(1.0D+00 - SINW*SINW))
         ONECOSW= 1.0D+00/COSW
         WWW     = ASIN(SINW)
C
         PMF1WG = PMF1WG + FCONST*WWW*WWW
 210  CONTINUE
C
      RETURN
      END
C*MODULE QUANPOB  *DECK ECMAP
!>
!> @brief    charmm correction map
!>
!> @author   Nandun Thellamurege
!>           - Apr 2011
!>
!> @details  this is slightly different from charmm
!>
      SUBROUTINE ECMAP(CORD,FFGRD,MAPLST,CMAPCO,CORDB,LSCMAPPMA,
     *                 LSCMAPPMB)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      PARAMETER (DEGREE=57.2957795130823D+00)
      PARAMETER (TOKCAL=627.509469D+00)
      PARAMETER (TOHART=1.0D+00/TOKCAL)
      PARAMETER (ONEFIF=1.0D+00/15.0D+00)
      PARAMETER (ONE=1.0D+00)
C
      INTEGER P1PSI, P2PSI, P3PSI, P4PSI, P1PHI, P2PHI, P3PHI, P4PHI
C
      DIMENSION CORD(3,*),FFGRD(3,*),MAPLST(6,*),CMAPCO(4,4,24,24,3),
     *          CORDB(3,*),LSCMAPPMA(*),LSCMAPPMB(*)
C
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
      COMMON /FFMDPA/ DT,DT2,TEMP0,PRES0,POLTOL,VIR(3),PMEAN,VOLAV,
     *                ENPAV,ENKAV,TEMPAV,BERENDT,BERENDP,VELMAX,
     *                PMEANX,PMEANY,PMEANZ,NSTEP,KMASTER,KOUTACT(2),
     *                IHESS,INTALG,ITSTAT,IPSTAT,JOUT,KOUT,LOUT
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
      COMMON /FFRMSD/ DIPT(8),TIMGYRA,TIMRALL,
     *                NATPDB,NGYRA,NDIEL,NRALL,
     *                LFFRALL0,NRMSD,LFFRMSD0,KFREEAB(201),
     *                NRIJMM,IJRMM(2,100),NRIJQM,IJRQM(2,100),
     *                NAIJKMM,IJKMM(3,100),NAIJKQM,IJKQM(3,100),
     *                NFIXMM,IFIXMM(200),NFIXQM,IFIXQM(200)
      COMMON /FFPARM/ NFFAT,NBOND,NANGL,NDIHR,NDIHB,NCMAP,NWAGG,
     *                N1213J,N14J,NLKQMM,IDOCHG,IDOPOL,IDOLJ,IDOCMAP
      COMMON /FFTYPE/ WT14LJ,WT14CH,C3BOND,C4BOND,C3ANGL,
     *                NFFTYP,NFFFILE,LJQMMM,LJQM,INTCHG,
     *                LJSIGMA,JTOPFILE(90),JPARFILE(90),
     *                JTOPAMIA(90),JTOPNTER(90),JTOPCTER(90),
     *                JTOPNUCA(90),JPARFIL2(90),JPARFIL3(90)
C
C     NANDUN THELLAMUREGE, APR 2011, LINCOLN
C     FENGCHAO CUI, HUI LI, MAY 2012
C
      PMF1CM=0.0D+00
      ENCMAP=0.0D+00
      IF(IDOCMAP.EQ.0 .OR. NFFTYP/10000.NE.2) RETURN
C
      DO 100 ICMAP=L1CMAP,L2CMAP
C        -- CALCULATE PHI AND PSI --
         ITYPE=MAPLST(6,ICMAP)
         P1PHI=MAPLST(1,ICMAP)
         P2PHI=MAPLST(2,ICMAP)
         P3PHI=MAPLST(3,ICMAP)
         P4PHI=MAPLST(4,ICMAP)
         X12PHI=CORD(1,P1PHI)-CORD(1,P2PHI)
         Y12PHI=CORD(2,P1PHI)-CORD(2,P2PHI)
         Z12PHI=CORD(3,P1PHI)-CORD(3,P2PHI)
         X13PHI=CORD(1,P1PHI)-CORD(1,P3PHI)
         Y13PHI=CORD(2,P1PHI)-CORD(2,P3PHI)
         Z13PHI=CORD(3,P1PHI)-CORD(3,P3PHI)
         X23PHI=CORD(1,P2PHI)-CORD(1,P3PHI)
         Y23PHI=CORD(2,P2PHI)-CORD(2,P3PHI)
         Z23PHI=CORD(3,P2PHI)-CORD(3,P3PHI)
         X34PHI=CORD(1,P3PHI)-CORD(1,P4PHI)
         Y34PHI=CORD(2,P3PHI)-CORD(2,P4PHI)
         Z34PHI=CORD(3,P3PHI)-CORD(3,P4PHI)
         R12PHI=SQRT(X12PHI*X12PHI+Y12PHI*Y12PHI+Z12PHI*Z12PHI)
         R23PHI=SQRT(X23PHI*X23PHI+Y23PHI*Y23PHI+Z23PHI*Z23PHI)
         R34PHI=SQRT(X34PHI*X34PHI+Y34PHI*Y34PHI+Z34PHI*Z34PHI)
         COS123PHI=(-(X12PHI*X23PHI)-(Y12PHI*Y23PHI)-(Z12PHI*Z23PHI))/
     *               (R12PHI*R23PHI)
         COS234PHI=(-(X23PHI*X34PHI)-(Y23PHI*Y34PHI)-(Z23PHI*Z34PHI))/
     *               (R23PHI*R34PHI)
         SIN2123PHI= 1.0D+00-COS123PHI*COS123PHI
         SIN2234PHI= 1.0D+00-COS234PHI*COS234PHI
         SIN123PHI = SQRT(ABS(SIN2123PHI))
         SIN234PHI = SQRT(ABS(SIN2234PHI))
         IF(ABS(SIN123PHI).LT.1.0D-06) GOTO 100
         IF(ABS(SIN234PHI).LT.1.0D-06) GOTO 100
         ONESINPHI = 1.0D+00/(SIN123PHI*SIN234PHI)
         COSPHI = ONESINPHI*(COS123PHI*COS234PHI-
     *             ((+X12PHI*X34PHI+Y12PHI*Y34PHI+Z12PHI*Z34PHI)/
     *              (R12PHI*R34PHI)))
         IF(COSPHI.GT. ONE) COSPHI= ONE
         IF(COSPHI.LT.-ONE) COSPHI=-ONE
         PHI   = ACOS(COSPHI)
C        -- DIHEDRAL ANGLE IS -180 TO +180 DEGREE
         XNORMPHI  = -Y23PHI*Z34PHI + Z23PHI*Y34PHI
         YNORMPHI  = -Z23PHI*X34PHI + X23PHI*Z34PHI
         ZNORMPHI  = -X23PHI*Y34PHI + Y23PHI*X34PHI
         DOTN12PHI =  X12PHI*XNORMPHI + Y12PHI*YNORMPHI +
     *                Z12PHI*ZNORMPHI
         IF(DOTN12PHI.LT.0.0D+00) PHI = -PHI       !  -180 TO 0 DEGREE
         PHI = PHI*DEGREE
C
C
         P1PSI=MAPLST(2,ICMAP)
         P2PSI=MAPLST(3,ICMAP)
         P3PSI=MAPLST(4,ICMAP)
         P4PSI=MAPLST(5,ICMAP)
         X12PSI=CORD(1,P1PSI)-CORD(1,P2PSI)
         Y12PSI=CORD(2,P1PSI)-CORD(2,P2PSI)
         Z12PSI=CORD(3,P1PSI)-CORD(3,P2PSI)
         X13PSI=CORD(1,P1PSI)-CORD(1,P3PSI)
         Y13PSI=CORD(2,P1PSI)-CORD(2,P3PSI)
         Z13PSI=CORD(3,P1PSI)-CORD(3,P3PSI)
         X23PSI=CORD(1,P2PSI)-CORD(1,P3PSI)
         Y23PSI=CORD(2,P2PSI)-CORD(2,P3PSI)
         Z23PSI=CORD(3,P2PSI)-CORD(3,P3PSI)
         X34PSI=CORD(1,P3PSI)-CORD(1,P4PSI)
         Y34PSI=CORD(2,P3PSI)-CORD(2,P4PSI)
         Z34PSI=CORD(3,P3PSI)-CORD(3,P4PSI)
         R12PSI=SQRT(X12PSI*X12PSI+Y12PSI*Y12PSI+Z12PSI*Z12PSI)
         R23PSI=SQRT(X23PSI*X23PSI+Y23PSI*Y23PSI+Z23PSI*Z23PSI)
         R34PSI=SQRT(X34PSI*X34PSI+Y34PSI*Y34PSI+Z34PSI*Z34PSI)
         COS123PSI=(-(X12PSI*X23PSI)-(Y12PSI*Y23PSI)-(Z12PSI*Z23PSI))/
     *               (R12PSI*R23PSI)
         COS234PSI=(-(X23PSI*X34PSI)-(Y23PSI*Y34PSI)-(Z23PSI*Z34PSI))/
     *               (R23PSI*R34PSI)
         SIN2123PSI= 1.0D+00-COS123PSI*COS123PSI
         SIN2234PSI= 1.0D+00-COS234PSI*COS234PSI
         SIN123PSI = SQRT(ABS(SIN2123PSI))
         SIN234PSI = SQRT(ABS(SIN2234PSI))
         IF(ABS(SIN123PSI).LT.1.0D-06) GOTO 100
         IF(ABS(SIN234PSI).LT.1.0D-06) GOTO 100
         ONESINPSI = 1.0D+00/(SIN123PSI*SIN234PSI)
         COSPSI = ONESINPSI*(COS123PSI*COS234PSI-
     *             ((+X12PSI*X34PSI+Y12PSI*Y34PSI+Z12PSI*Z34PSI)/
     *              (R12PSI*R34PSI)))
         IF(COSPSI.GT. ONE) COSPHI= ONE
         IF(COSPSI.LT.-ONE) COSPHI=-ONE
         PSI   = ACOS(COSPSI)
C        -- DIHEDRAL ANGLE IS -180 TO +180 DEGREE
         XNORMPSI  = -Y23PSI*Z34PSI + Z23PSI*Y34PSI
         YNORMPSI  = -Z23PSI*X34PSI + X23PSI*Z34PSI
         ZNORMPSI  = -X23PSI*Y34PSI + Y23PSI*X34PSI
         DOTN12PSI =  X12PSI*XNORMPSI + Y12PSI*YNORMPSI +
     *                Z12PSI*ZNORMPSI
         IF(DOTN12PSI.LT.0.0D+00) PSI = -PSI       !  -180 TO 0 DEGREE
         PSI = PSI*DEGREE
C
C
         IPHI  = INT((PHI+180.0D+00)*ONEFIF) + 1  ! FROM 1 TO 24
         IF(IPHI.EQ.25) IPHI = 1
         JPSI  = INT((PSI+180.0D+00)*ONEFIF) + 1  ! FROM 1 TO 24
         IF(JPSI.EQ.25) JPSI = 1
         DPHI  = (PHI+180.0D+00)*ONEFIF - IPHI + 1
         DPSI  = (PSI+180.0D+00)*ONEFIF - JPSI + 1
C
C
         DO ICPHI =1,4
            CCCPHI = 1.0D+00
            IF(ICPHI.GT.1) CCCPHI = DPHI**(ICPHI-1)
            DO JCPSI=1,4
               CCCPSI = 1.0D+00
               IF(JCPSI.GT.1) CCCPSI = DPSI**(JCPSI-1)
               CIJ = CMAPCO(JCPSI,ICPHI,JPSI,IPHI,ITYPE)
               ENCMAP=ENCMAP + CIJ*CCCPHI*CCCPSI
            ENDDO
         ENDDO
C
C        - CALCULATING DERIVATIVES
C
         PHIFAC = 0.0D+00
         DO ICPHI =1,4
            DO JCPSI=1,4
               CIJ = CMAPCO(JCPSI,ICPHI,JPSI,IPHI,ITYPE)
               PHIFAC=PHIFAC +
     *                (ICPHI-1)*CIJ*DPHI**(ICPHI-2)*DPSI**(JCPSI-1)
            ENDDO
         ENDDO
         PHIFAC = PHIFAC*TOHART*ONEFIF*DEGREE
C
         A123PHI = 1.0D+00/(R12PHI*SIN2123PHI)
         A432PHI = 1.0D+00/(R34PHI*SIN2234PHI)
         B123PHI = R12PHI*COS123PHI/R23PHI
         B432PHI = R34PHI*COS234PHI/R23PHI
         C123PHI = B123PHI-1.0D+00
C
         DUMY        =  PHIFAC*A123PHI/(R12PHI*R23PHI)
         DEX1        =  DUMY*(-Y12PHI*Z23PHI+Z12PHI*Y23PHI)
         DEY1        =  DUMY*(-Z12PHI*X23PHI+X12PHI*Z23PHI)
         DEZ1        =  DUMY*(-X12PHI*Y23PHI+Y12PHI*X23PHI)
         FFGRD(1,P1PHI)=FFGRD(1,P1PHI)+DEX1
         FFGRD(2,P1PHI)=FFGRD(2,P1PHI)+DEY1
         FFGRD(3,P1PHI)=FFGRD(3,P1PHI)+DEZ1
C
         DUMY        =  PHIFAC*(-A432PHI)/(R23PHI*R34PHI)
         DEX4        =  DUMY*(-Y23PHI*Z34PHI+Z23PHI*Y34PHI)
         DEY4        =  DUMY*(-Z23PHI*X34PHI+X23PHI*Z34PHI)
         DEZ4        =  DUMY*(-X23PHI*Y34PHI+Y23PHI*X34PHI)
         FFGRD(1,P4PHI)=FFGRD(1,P4PHI)+DEX4
         FFGRD(2,P4PHI)=FFGRD(2,P4PHI)+DEY4
         FFGRD(3,P4PHI)=FFGRD(3,P4PHI)+DEZ4
C
         DEX2        =  C123PHI*DEX1 - B432PHI*DEX4
         DEY2        =  C123PHI*DEY1 - B432PHI*DEY4
         DEZ2        =  C123PHI*DEZ1 - B432PHI*DEZ4
         FFGRD(1,P2PHI)=FFGRD(1,P2PHI)+DEX2
         FFGRD(2,P2PHI)=FFGRD(2,P2PHI)+DEY2
         FFGRD(3,P2PHI)=FFGRD(3,P2PHI)+DEZ2
C
         FFGRD(1,P3PHI)=FFGRD(1,P3PHI)-DEX1-DEX2-DEX4
         FFGRD(2,P3PHI)=FFGRD(2,P3PHI)-DEY1-DEY2-DEY4
         FFGRD(3,P3PHI)=FFGRD(3,P3PHI)-DEZ1-DEZ2-DEZ4
C
         VIR(1)  =VIR(1) + DEX1*X13PHI + DEX2*X23PHI - DEX4*X34PHI
         VIR(2)  =VIR(2) + DEY1*Y13PHI + DEY2*Y23PHI - DEY4*Y34PHI
         VIR(3)  =VIR(3) + DEZ1*Z13PHI + DEZ2*Z23PHI - DEZ4*Z34PHI
         IYES = 0
         JYES = 0
         KYES = 0
         LYES = 0
         DO KFIX=1,NFIXMM
            IF(P1PHI.EQ.IFIXMM(KFIX)) IYES = 1
            IF(P2PHI.EQ.IFIXMM(KFIX)) JYES = 1
            IF(P3PHI.EQ.IFIXMM(KFIX)) KYES = 1
            IF(P4PHI.EQ.IFIXMM(KFIX)) LYES = 1
         ENDDO
         IF(IYES.EQ.1.AND.KYES.EQ.1) THEN
            VIR(1)    =VIR(1)     - DEX1*X13PHI
            VIR(2)    =VIR(2)     - DEY1*Y13PHI
            VIR(3)    =VIR(3)     - DEZ1*Z13PHI
         END IF
         IF(JYES.EQ.1.AND.KYES.EQ.1) THEN
            VIR(1)    =VIR(1)     - DEX2*X23PHI
            VIR(2)    =VIR(2)     - DEY2*Y23PHI
            VIR(3)    =VIR(3)     - DEZ2*Z23PHI
         END IF
         IF(KYES.EQ.1.AND.LYES.EQ.1) THEN
            VIR(1)    =VIR(1)     + DEX4*X34PHI
            VIR(2)    =VIR(2)     + DEY4*Y34PHI
            VIR(3)    =VIR(3)     + DEZ4*Z34PHI
         END IF
C
C
         PSIFAC=0.0D+00
         DO ICPHI =1,4
            DO JCPSI=1,4
               CIJ = CMAPCO(JCPSI,ICPHI,JPSI,IPHI,ITYPE)
               PSIFAC=PSIFAC+
     *                (JCPSI-1)*CIJ*DPHI**(ICPHI-1)*DPSI**(JCPSI-2)
            ENDDO
         ENDDO
         PSIFAC = PSIFAC*TOHART*ONEFIF*DEGREE
C
         A123PSI = 1.0D+00/(R12PSI*SIN2123PSI)
         A432PSI = 1.0D+00/(R34PSI*SIN2234PSI)
         B123PSI = R12PSI*COS123PSI/R23PSI
         B432PSI = R34PSI*COS234PSI/R23PSI
         C123PSI = B123PSI-1.0D+00
C
         DUMY         = PSIFAC*A123PSI/(R12PSI*R23PSI)
         DEX1         = DUMY*(-Y12PSI*Z23PSI+Z12PSI*Y23PSI)
         DEY1         = DUMY*(-Z12PSI*X23PSI+X12PSI*Z23PSI)
         DEZ1         = DUMY*(-X12PSI*Y23PSI+Y12PSI*X23PSI)
         FFGRD(1,P1PSI) = FFGRD(1,P1PSI)+DEX1
         FFGRD(2,P1PSI) = FFGRD(2,P1PSI)+DEY1
         FFGRD(3,P1PSI) = FFGRD(3,P1PSI)+DEZ1
C
         DUMY         = PSIFAC*(-A432PSI)/(R23PSI*R34PSI)
         DEX4         = DUMY*(-Y23PSI*Z34PSI+Z23PSI*Y34PSI)
         DEY4         = DUMY*(-Z23PSI*X34PSI+X23PSI*Z34PSI)
         DEZ4         = DUMY*(-X23PSI*Y34PSI+Y23PSI*X34PSI)
         FFGRD(1,P4PSI) = FFGRD(1,P4PSI)+DEX4
         FFGRD(2,P4PSI) = FFGRD(2,P4PSI)+DEY4
         FFGRD(3,P4PSI) = FFGRD(3,P4PSI)+DEZ4
C
         DEX2         = C123PSI*DEX1 - B432PSI*DEX4
         DEY2         = C123PSI*DEY1 - B432PSI*DEY4
         DEZ2         = C123PSI*DEZ1 - B432PSI*DEZ4
         FFGRD(1,P2PSI) = FFGRD(1,P2PSI)+DEX2
         FFGRD(2,P2PSI) = FFGRD(2,P2PSI)+DEY2
         FFGRD(3,P2PSI) = FFGRD(3,P2PSI)+DEZ2
C
         FFGRD(1,P3PSI) = FFGRD(1,P3PSI)-DEX1-DEX2-DEX4
         FFGRD(2,P3PSI) = FFGRD(2,P3PSI)-DEY1-DEY2-DEY4
         FFGRD(3,P3PSI) = FFGRD(3,P3PSI)-DEZ1-DEZ2-DEZ4
C
         VIR(1)  =VIR(1) + DEX1*X13PSI + DEX2*X23PSI - DEX4*X34PSI
         VIR(2)  =VIR(2) + DEY1*Y13PSI + DEY2*Y23PSI - DEY4*Y34PSI
         VIR(3)  =VIR(3) + DEZ1*Z13PSI + DEZ2*Z23PSI - DEZ4*Z34PSI
         IYES = 0
         JYES = 0
         KYES = 0
         LYES = 0
         DO KFIX=1,NFIXMM
            IF(P1PSI.EQ.IFIXMM(KFIX)) IYES = 1
            IF(P2PSI.EQ.IFIXMM(KFIX)) JYES = 1
            IF(P3PSI.EQ.IFIXMM(KFIX)) KYES = 1
            IF(P4PSI.EQ.IFIXMM(KFIX)) LYES = 1
         ENDDO
         IF(IYES.EQ.1.AND.KYES.EQ.1) THEN
            VIR(1)    =VIR(1)     - DEX1*X13PSI
            VIR(2)    =VIR(2)     - DEY1*Y13PSI
            VIR(3)    =VIR(3)     - DEZ1*Z13PSI
         END IF
         IF(JYES.EQ.1.AND.KYES.EQ.1) THEN
            VIR(1)    =VIR(1)     - DEX2*X23PSI
            VIR(2)    =VIR(2)     - DEY2*Y23PSI
            VIR(3)    =VIR(3)     - DEZ2*Z23PSI
         END IF
         IF(KYES.EQ.1.AND.LYES.EQ.1) THEN
            VIR(1)    =VIR(1)     + DEX4*X34PSI
            VIR(2)    =VIR(2)     + DEY4*Y34PSI
            VIR(3)    =VIR(3)     + DEZ4*Z34PSI
         END IF
 100  CONTINUE
      ENCMAP = ENCMAP*TOHART
C
C
      DO 200 III=L1CMAPPMA,L2CMAPPMA
C        -- CALCULATE PHI AND PSI --
         ICMAP=LSCMAPPMA(III)
         ITYPE=MAPLST(6,ICMAP)
         P1PHI=MAPLST(1,ICMAP)
         P2PHI=MAPLST(2,ICMAP)
         P3PHI=MAPLST(3,ICMAP)
         P4PHI=MAPLST(4,ICMAP)
         P1PSI=MAPLST(2,ICMAP)
         P2PSI=MAPLST(3,ICMAP)
         P3PSI=MAPLST(4,ICMAP)
         P4PSI=MAPLST(5,ICMAP)
         X12PHI=CORD(1,P1PHI)-CORD(1,P2PHI)
         Y12PHI=CORD(2,P1PHI)-CORD(2,P2PHI)
         Z12PHI=CORD(3,P1PHI)-CORD(3,P2PHI)
         X13PHI=CORD(1,P1PHI)-CORD(1,P3PHI)
         Y13PHI=CORD(2,P1PHI)-CORD(2,P3PHI)
         Z13PHI=CORD(3,P1PHI)-CORD(3,P3PHI)
         X23PHI=CORD(1,P2PHI)-CORD(1,P3PHI)
         Y23PHI=CORD(2,P2PHI)-CORD(2,P3PHI)
         Z23PHI=CORD(3,P2PHI)-CORD(3,P3PHI)
         X34PHI=CORD(1,P3PHI)-CORD(1,P4PHI)
         Y34PHI=CORD(2,P3PHI)-CORD(2,P4PHI)
         Z34PHI=CORD(3,P3PHI)-CORD(3,P4PHI)
         R12PHI=SQRT(X12PHI*X12PHI+Y12PHI*Y12PHI+Z12PHI*Z12PHI)
         R23PHI=SQRT(X23PHI*X23PHI+Y23PHI*Y23PHI+Z23PHI*Z23PHI)
         R34PHI=SQRT(X34PHI*X34PHI+Y34PHI*Y34PHI+Z34PHI*Z34PHI)
         COS123PHI=(-(X12PHI*X23PHI)-(Y12PHI*Y23PHI)-(Z12PHI*Z23PHI))/
     *               (R12PHI*R23PHI)
         COS234PHI=(-(X23PHI*X34PHI)-(Y23PHI*Y34PHI)-(Z23PHI*Z34PHI))/
     *               (R23PHI*R34PHI)
         SIN2123PHI= 1.0D+00-COS123PHI*COS123PHI
         SIN2234PHI= 1.0D+00-COS234PHI*COS234PHI
         SIN123PHI = SQRT(ABS(SIN2123PHI))
         SIN234PHI = SQRT(ABS(SIN2234PHI))
         IF(ABS(SIN123PHI).LT.1.0D-06) GOTO 200
         IF(ABS(SIN234PHI).LT.1.0D-06) GOTO 200
         ONESINPHI = 1.0D+00/(SIN123PHI*SIN234PHI)
         COSPHI = ONESINPHI*(COS123PHI*COS234PHI-
     *             ((+X12PHI*X34PHI+Y12PHI*Y34PHI+Z12PHI*Z34PHI)/
     *              (R12PHI*R34PHI)))
         IF(COSPHI.GT. ONE) COSPHI= ONE
         IF(COSPHI.LT.-ONE) COSPHI=-ONE
         PHI   = ACOS(COSPHI)
C        -- DIHEDRAL ANGLE IS -180 TO +180 DEGREE
         XNORMPHI  = -Y23PHI*Z34PHI + Z23PHI*Y34PHI
         YNORMPHI  = -Z23PHI*X34PHI + X23PHI*Z34PHI
         ZNORMPHI  = -X23PHI*Y34PHI + Y23PHI*X34PHI
         DOTN12PHI =  X12PHI*XNORMPHI + Y12PHI*YNORMPHI +
     *                Z12PHI*ZNORMPHI
         IF(DOTN12PHI.LT.0.0D+00) PHI = -PHI       !  -180 TO 0 DEGREE
         PHI  = PHI*DEGREE
C
C
         X12PSI=CORD(1,P1PSI)-CORD(1,P2PSI)
         Y12PSI=CORD(2,P1PSI)-CORD(2,P2PSI)
         Z12PSI=CORD(3,P1PSI)-CORD(3,P2PSI)
         X13PSI=CORD(1,P1PSI)-CORD(1,P3PSI)
         Y13PSI=CORD(2,P1PSI)-CORD(2,P3PSI)
         Z13PSI=CORD(3,P1PSI)-CORD(3,P3PSI)
         X23PSI=CORD(1,P2PSI)-CORD(1,P3PSI)
         Y23PSI=CORD(2,P2PSI)-CORD(2,P3PSI)
         Z23PSI=CORD(3,P2PSI)-CORD(3,P3PSI)
         X34PSI=CORD(1,P3PSI)-CORD(1,P4PSI)
         Y34PSI=CORD(2,P3PSI)-CORD(2,P4PSI)
         Z34PSI=CORD(3,P3PSI)-CORD(3,P4PSI)
         R12PSI=SQRT(X12PSI*X12PSI+Y12PSI*Y12PSI+Z12PSI*Z12PSI)
         R23PSI=SQRT(X23PSI*X23PSI+Y23PSI*Y23PSI+Z23PSI*Z23PSI)
         R34PSI=SQRT(X34PSI*X34PSI+Y34PSI*Y34PSI+Z34PSI*Z34PSI)
         COS123PSI=(-(X12PSI*X23PSI)-(Y12PSI*Y23PSI)-(Z12PSI*Z23PSI))/
     *               (R12PSI*R23PSI)
         COS234PSI=(-(X23PSI*X34PSI)-(Y23PSI*Y34PSI)-(Z23PSI*Z34PSI))/
     *               (R23PSI*R34PSI)
         SIN2123PSI= 1.0D+00-COS123PSI*COS123PSI
         SIN2234PSI= 1.0D+00-COS234PSI*COS234PSI
         SIN123PSI = SQRT(ABS(SIN2123PSI))
         SIN234PSI = SQRT(ABS(SIN2234PSI))
         IF(ABS(SIN123PSI).LT.1.0D-06) GOTO 200
         IF(ABS(SIN234PSI).LT.1.0D-06) GOTO 200
         ONESINPSI = 1.0D+00/(SIN123PSI*SIN234PSI)
         COSPSI = ONESINPSI*(COS123PSI*COS234PSI-
     *             ((+X12PSI*X34PSI+Y12PSI*Y34PSI+Z12PSI*Z34PSI)/
     *              (R12PSI*R34PSI)))
         IF(COSPSI.GT. ONE) COSPHI= ONE
         IF(COSPSI.LT.-ONE) COSPHI=-ONE
         PSI   = ACOS(COSPSI)
C        -- DIHEDRAL ANGLE IS -180 TO +180 DEGREE
         XNORMPSI  = -Y23PSI*Z34PSI + Z23PSI*Y34PSI
         YNORMPSI  = -Z23PSI*X34PSI + X23PSI*Z34PSI
         ZNORMPSI  = -X23PSI*Y34PSI + Y23PSI*X34PSI
         DOTN12PSI =  X12PSI*XNORMPSI + Y12PSI*YNORMPSI +
     *                Z12PSI*ZNORMPSI
         IF(DOTN12PSI.LT.0.0D+00) PSI = -PSI       !  -180 TO 0 DEGREE
         PSI  = PSI*DEGREE
C
C
         IPHI  = INT((PHI+180.0D+00)*ONEFIF) + 1  ! FROM 1 TO 24
         IF(IPHI.EQ.25) IPHI = 1
         JPSI  = INT((PSI+180.0D+00)*ONEFIF) + 1  ! FROM 1 TO 24
         IF(JPSI.EQ.25) JPSI = 1
         DPHI  = (PHI+180.0D+00)*ONEFIF - IPHI + 1
         DPSI  = (PSI+180.0D+00)*ONEFIF - JPSI + 1
C
C
         DO ICPHI =1,4
            CCCPHI = 1.0D+00
            IF(ICPHI.GT.1) CCCPHI = DPHI**(ICPHI-1)
            DO JCPSI=1,4
               CCCPSI = 1.0D+00
               IF(JCPSI.GT.1) CCCPSI = DPSI**(JCPSI-1)
               CIJ = CMAPCO(JCPSI,ICPHI,JPSI,IPHI,ITYPE)
               PMF1CM = PMF1CM - CIJ*CCCPHI*CCCPSI
            ENDDO
         ENDDO
C
 200  CONTINUE
C
C
      DO 210 III=L1CMAPPMB,L2CMAPPMB
C        -- CALCULATE PHI AND PSI --
         ICMAP=LSCMAPPMB(III)
         ITYPE=MAPLST(6,ICMAP)
         P1PHI=MAPLST(1,ICMAP)
         P2PHI=MAPLST(2,ICMAP)
         P3PHI=MAPLST(3,ICMAP)
         P4PHI=MAPLST(4,ICMAP)
         P1PSI=MAPLST(2,ICMAP)
         P2PSI=MAPLST(3,ICMAP)
         P3PSI=MAPLST(4,ICMAP)
         P4PSI=MAPLST(5,ICMAP)
         X12PHI=CORD(1,P1PHI)-CORD(1,P2PHI)
         Y12PHI=CORD(2,P1PHI)-CORD(2,P2PHI)
         Z12PHI=CORD(3,P1PHI)-CORD(3,P2PHI)
         X13PHI=CORD(1,P1PHI)-CORD(1,P3PHI)
         Y13PHI=CORD(2,P1PHI)-CORD(2,P3PHI)
         Z13PHI=CORD(3,P1PHI)-CORD(3,P3PHI)
         X23PHI=CORD(1,P2PHI)-CORD(1,P3PHI)
         Y23PHI=CORD(2,P2PHI)-CORD(2,P3PHI)
         Z23PHI=CORD(3,P2PHI)-CORD(3,P3PHI)
         X34PHI=CORD(1,P3PHI)-CORD(1,P4PHI)
         Y34PHI=CORD(2,P3PHI)-CORD(2,P4PHI)
         Z34PHI=CORD(3,P3PHI)-CORD(3,P4PHI)
C
         IF(IFEPTYP.EQ.2) THEN   !   NFIXMM=KFREEAB FOR IFEPTYP=2
            DO KFIX=1,NFIXMM
               IF(P1PHI.EQ.IFIXMM(KFIX)) THEN
                  X12PHI = X12PHI - CORD(1,P1PHI) + CORDB(1,P1PHI)
                  Y12PHI = Y12PHI - CORD(2,P1PHI) + CORDB(2,P1PHI)
                  Z12PHI = Z12PHI - CORD(3,P1PHI) + CORDB(3,P1PHI)
                  X13PHI = X13PHI - CORD(1,P1PHI) + CORDB(1,P1PHI)
                  Y13PHI = Y13PHI - CORD(2,P1PHI) + CORDB(2,P1PHI)
                  Z13PHI = Z13PHI - CORD(3,P1PHI) + CORDB(3,P1PHI)
               END IF
               IF(P2PHI.EQ.IFIXMM(KFIX)) THEN
                  X12PHI = X12PHI + CORD(1,P2PHI) - CORDB(1,P2PHI)
                  Y12PHI = Y12PHI + CORD(2,P2PHI) - CORDB(2,P2PHI)
                  Z12PHI = Z12PHI + CORD(3,P2PHI) - CORDB(3,P2PHI)
                  X23PHI = X23PHI - CORD(1,P2PHI) + CORDB(1,P2PHI)
                  Y23PHI = Y23PHI - CORD(2,P2PHI) + CORDB(2,P2PHI)
                  Z23PHI = Z23PHI - CORD(3,P2PHI) + CORDB(3,P2PHI)
               END IF
               IF(P3PHI.EQ.IFIXMM(KFIX)) THEN
                  X13PHI = X13PHI + CORD(1,P3PHI) - CORDB(1,P3PHI)
                  Y13PHI = Y13PHI + CORD(2,P3PHI) - CORDB(2,P3PHI)
                  Z13PHI = Z13PHI + CORD(3,P3PHI) - CORDB(3,P3PHI)
                  X23PHI = X23PHI + CORD(1,P3PHI) - CORDB(1,P3PHI)
                  Y23PHI = Y23PHI + CORD(2,P3PHI) - CORDB(2,P3PHI)
                  Z23PHI = Z23PHI + CORD(3,P3PHI) - CORDB(3,P3PHI)
                  X34PHI = X34PHI - CORD(1,P3PHI) + CORDB(1,P3PHI)
                  Y34PHI = Y34PHI - CORD(2,P3PHI) + CORDB(2,P3PHI)
                  Z34PHI = Z34PHI - CORD(3,P3PHI) + CORDB(3,P3PHI)
               END IF
               IF(P4PHI.EQ.IFIXMM(KFIX)) THEN
                  X34PHI = X34PHI + CORD(1,P4PHI) - CORDB(1,P4PHI)
                  Y34PHI = Y34PHI + CORD(2,P4PHI) - CORDB(2,P4PHI)
                  Z34PHI = Z34PHI + CORD(3,P4PHI) - CORDB(3,P4PHI)
               END IF
            ENDDO
         END IF
C
         R12PHI=SQRT(X12PHI*X12PHI+Y12PHI*Y12PHI+Z12PHI*Z12PHI)
         R23PHI=SQRT(X23PHI*X23PHI+Y23PHI*Y23PHI+Z23PHI*Z23PHI)
         R34PHI=SQRT(X34PHI*X34PHI+Y34PHI*Y34PHI+Z34PHI*Z34PHI)
         COS123PHI=(-(X12PHI*X23PHI)-(Y12PHI*Y23PHI)-(Z12PHI*Z23PHI))/
     *               (R12PHI*R23PHI)
         COS234PHI=(-(X23PHI*X34PHI)-(Y23PHI*Y34PHI)-(Z23PHI*Z34PHI))/
     *               (R23PHI*R34PHI)
         SIN2123PHI= 1.0D+00-COS123PHI*COS123PHI
         SIN2234PHI= 1.0D+00-COS234PHI*COS234PHI
         SIN123PHI = SQRT(ABS(SIN2123PHI))
         SIN234PHI = SQRT(ABS(SIN2234PHI))
         IF(ABS(SIN123PHI).LT.1.0D-06) GOTO 210
         IF(ABS(SIN234PHI).LT.1.0D-06) GOTO 210
         ONESINPHI = 1.0D+00/(SIN123PHI*SIN234PHI)
         COSPHI = ONESINPHI*(COS123PHI*COS234PHI-
     *             ((+X12PHI*X34PHI+Y12PHI*Y34PHI+Z12PHI*Z34PHI)/
     *              (R12PHI*R34PHI)))
         IF(COSPHI.GT. ONE) COSPHI= ONE
         IF(COSPHI.LT.-ONE) COSPHI=-ONE
         PHI   = ACOS(COSPHI)
C        -- DIHEDRAL ANGLE IS -180 TO +180 DEGREE
         XNORMPHI  = -Y23PHI*Z34PHI + Z23PHI*Y34PHI
         YNORMPHI  = -Z23PHI*X34PHI + X23PHI*Z34PHI
         ZNORMPHI  = -X23PHI*Y34PHI + Y23PHI*X34PHI
         DOTN12PHI =  X12PHI*XNORMPHI + Y12PHI*YNORMPHI +
     *                Z12PHI*ZNORMPHI
         IF(DOTN12PHI.LT.0.0D+00) PHI = -PHI       !  -180 TO 0 DEGREE
         PHI  = PHI*DEGREE
C
C
         X12PSI=CORD(1,P1PSI)-CORD(1,P2PSI)
         Y12PSI=CORD(2,P1PSI)-CORD(2,P2PSI)
         Z12PSI=CORD(3,P1PSI)-CORD(3,P2PSI)
         X13PSI=CORD(1,P1PSI)-CORD(1,P3PSI)
         Y13PSI=CORD(2,P1PSI)-CORD(2,P3PSI)
         Z13PSI=CORD(3,P1PSI)-CORD(3,P3PSI)
         X23PSI=CORD(1,P2PSI)-CORD(1,P3PSI)
         Y23PSI=CORD(2,P2PSI)-CORD(2,P3PSI)
         Z23PSI=CORD(3,P2PSI)-CORD(3,P3PSI)
         X34PSI=CORD(1,P3PSI)-CORD(1,P4PSI)
         Y34PSI=CORD(2,P3PSI)-CORD(2,P4PSI)
         Z34PSI=CORD(3,P3PSI)-CORD(3,P4PSI)
C
         IF(IFEPTYP.EQ.2) THEN  !  NFIXMM=KFREEAB FOR IFEPTYP=2
            DO KFIX=1,NFIXMM
               IF(P1PSI.EQ.IFIXMM(KFIX)) THEN
                  X12PSI = X12PSI - CORD(1,P1PSI) + CORDB(1,P1PSI)
                  Y12PSI = Y12PSI - CORD(2,P1PSI) + CORDB(2,P1PSI)
                  Z12PSI = Z12PSI - CORD(3,P1PSI) + CORDB(3,P1PSI)
                  X13PSI = X13PSI - CORD(1,P1PSI) + CORDB(1,P1PSI)
                  Y13PSI = Y13PSI - CORD(2,P1PSI) + CORDB(2,P1PSI)
                  Z13PSI = Z13PSI - CORD(3,P1PSI) + CORDB(3,P1PSI)
               END IF
               IF(P2PSI.EQ.IFIXMM(KFIX)) THEN
                  X12PSI = X12PSI + CORD(1,P2PSI) - CORDB(1,P2PSI)
                  Y12PSI = Y12PSI + CORD(2,P2PSI) - CORDB(2,P2PSI)
                  Z12PSI = Z12PSI + CORD(3,P2PSI) - CORDB(3,P2PSI)
                  X23PSI = X23PSI - CORD(1,P2PSI) + CORDB(1,P2PSI)
                  Y23PSI = Y23PSI - CORD(2,P2PSI) + CORDB(2,P2PSI)
                  Z23PSI = Z23PSI - CORD(3,P2PSI) + CORDB(3,P2PSI)
               END IF
               IF(P3PSI.EQ.IFIXMM(KFIX)) THEN
                  X13PSI = X13PSI + CORD(1,P3PSI) - CORDB(1,P3PSI)
                  Y13PSI = Y13PSI + CORD(2,P3PSI) - CORDB(2,P3PSI)
                  Z13PSI = Z13PSI + CORD(3,P3PSI) - CORDB(3,P3PSI)
                  X23PSI = X23PSI + CORD(1,P3PSI) - CORDB(1,P3PSI)
                  Y23PSI = Y23PSI + CORD(2,P3PSI) - CORDB(2,P3PSI)
                  Z23PSI = Z23PSI + CORD(3,P3PSI) - CORDB(3,P3PSI)
                  X34PSI = X34PSI - CORD(1,P3PSI) + CORDB(1,P3PSI)
                  Y34PSI = Y34PSI - CORD(2,P3PSI) + CORDB(2,P3PSI)
                  Z34PSI = Z34PSI - CORD(3,P3PSI) + CORDB(3,P3PSI)
               END IF
               IF(P4PSI.EQ.IFIXMM(KFIX)) THEN
                  X34PSI = X34PSI + CORD(1,P4PSI) - CORDB(1,P4PSI)
                  Y34PSI = Y34PSI + CORD(2,P4PSI) - CORDB(2,P4PSI)
                  Z34PSI = Z34PSI + CORD(3,P4PSI) - CORDB(3,P4PSI)
               END IF
            ENDDO
         END IF
C
         R12PSI=SQRT(X12PSI*X12PSI+Y12PSI*Y12PSI+Z12PSI*Z12PSI)
         R23PSI=SQRT(X23PSI*X23PSI+Y23PSI*Y23PSI+Z23PSI*Z23PSI)
         R34PSI=SQRT(X34PSI*X34PSI+Y34PSI*Y34PSI+Z34PSI*Z34PSI)
         COS123PSI=(-(X12PSI*X23PSI)-(Y12PSI*Y23PSI)-(Z12PSI*Z23PSI))/
     *               (R12PSI*R23PSI)
         COS234PSI=(-(X23PSI*X34PSI)-(Y23PSI*Y34PSI)-(Z23PSI*Z34PSI))/
     *               (R23PSI*R34PSI)
         SIN2123PSI= 1.0D+00-COS123PSI*COS123PSI
         SIN2234PSI= 1.0D+00-COS234PSI*COS234PSI
         SIN123PSI = SQRT(ABS(SIN2123PSI))
         SIN234PSI = SQRT(ABS(SIN2234PSI))
         IF(ABS(SIN123PSI).LT.1.0D-06) GOTO 210
         IF(ABS(SIN234PSI).LT.1.0D-06) GOTO 210
         ONESINPSI = 1.0D+00/(SIN123PSI*SIN234PSI)
         COSPSI = ONESINPSI*(COS123PSI*COS234PSI-
     *             ((+X12PSI*X34PSI+Y12PSI*Y34PSI+Z12PSI*Z34PSI)/
     *              (R12PSI*R34PSI)))
         IF(COSPSI.GT. ONE) COSPHI= ONE
         IF(COSPSI.LT.-ONE) COSPHI=-ONE
         PSI   = ACOS(COSPSI)
C        -- DIHEDRAL ANGLE IS -180 TO +180 DEGREE
         XNORMPSI  = -Y23PSI*Z34PSI + Z23PSI*Y34PSI
         YNORMPSI  = -Z23PSI*X34PSI + X23PSI*Z34PSI
         ZNORMPSI  = -X23PSI*Y34PSI + Y23PSI*X34PSI
         DOTN12PSI =  X12PSI*XNORMPSI + Y12PSI*YNORMPSI +
     *                Z12PSI*ZNORMPSI
         IF(DOTN12PSI.LT.0.0D+00) PSI = -PSI       !  -180 TO 0 DEGREE
         PSI  = PSI*DEGREE
C
C
         IPHI  = INT((PHI+180.0D+00)*ONEFIF) + 1  ! FROM 1 TO 24
         IF(IPHI.EQ.25) IPHI = 1
         JPSI  = INT((PSI+180.0D+00)*ONEFIF) + 1  ! FROM 1 TO 24
         IF(JPSI.EQ.25) JPSI = 1
         DPHI  = (PHI+180.0D+00)*ONEFIF - IPHI + 1
         DPSI  = (PSI+180.0D+00)*ONEFIF - JPSI + 1
C
C
         DO ICPHI =1,4
            CCCPHI = 1.0D+00
            IF(ICPHI.GT.1) CCCPHI = DPHI**(ICPHI-1)
            DO JCPSI=1,4
               CCCPSI = 1.0D+00
               IF(JCPSI.GT.1) CCCPSI = DPSI**(JCPSI-1)
               CIJ = CMAPCO(JCPSI,ICPHI,JPSI,IPHI,ITYPE)
               PMF1CM = PMF1CM + CIJ*CCCPHI*CCCPSI
            ENDDO
         ENDDO
C
 210  CONTINUE
C
      PMF1CM = PMF1CM*TOHART
C
      RETURN
      END
C*MODULE QUANPOB  *DECK ELJ126
!>
!> @brief    LJ energy
!>
!> @author   Nandun Thellamurege, Hui Li
!>           - Jan 2011
!>
!> @details  force field LJ term
!>
      SUBROUTINE ELJ126(CORD,FFGRD,SIG,EPS,SIG2,EPS2,L14J,
     *                  NONLS1,L1213J,
     *                  SIGB,EPSB,SIG2B,EPS2B,
     *                  NONLSA,NONLSB,L1213A,L1213B,L14A,L14B,
     *                  NONLSPMA,L1213PMA,L14PMA,CORDB,
     *                  NONLSPMB,L1213PMB,L14PMB,FCLJTP,
     *                  NTYPE)
      use mx_limits, only: mxatm
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      LOGICAL GOPARR,MASWRK,DSKWRK
C
      PARAMETER (ZERO=0.0D+00)
      PARAMETER (ONE=1.0D+00)
      PARAMETER (SEVEN=7.0D+00)
      PARAMETER (TWO=2.0D+00)
C
      DIMENSION CORD(3,*),FFGRD(3,*),SIG(*),EPS(*),SIG2(*),EPS2(*),
     *          L14J(2,*),NONLS1(2,*),
     *          L1213J(2,*),SIGB(*),EPSB(*),SIG2B(*),EPS2B(*),
     *          NONLSA(2,*),NONLSB(2,*),L1213A(2,*),L1213B(2,*),
     *          L14A(2,*),L14B(2,*),NONLSPMA(2,*),L1213PMA(2,*),
     *          L14PMA(2,*),CORDB(3,*),
     *          NONLSPMB(2,*),L1213PMB(2,*),L14PMB(2,*),
     *          FCLJTP(2,MXMMTP,*),NTYPE(*)
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
      COMMON /FFMPT2/ MXMMTP,LFFKBLST,LFFFCSTBD,LFFFCDIHR,
     *                LFFFCLJTP,LFFNTYPE,
     *                LFF2KBLST,LFF2FCSTBD,LFF2FCDIHR,
     *                LFF2FCLJTP,LFF2NTYPE
      COMMON /FFPARM/ NFFAT,NBOND,NANGL,NDIHR,NDIHB,NCMAP,NWAGG,
     *                N1213J,N14J,NLKQMM,IDOCHG,IDOPOL,IDOLJ,IDOCMAP
      COMMON /FFMDPA/ DT,DT2,TEMP0,PRES0,POLTOL,VIR(3),PMEAN,VOLAV,
     *                ENPAV,ENKAV,TEMPAV,BERENDT,BERENDP,VELMAX,
     *                PMEANX,PMEANY,PMEANZ,NSTEP,KMASTER,KOUTACT(2),
     *                IHESS,INTALG,ITSTAT,IPSTAT,JOUT,KOUT,LOUT
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
      COMMON /FFTYPE/ WT14LJ,WT14CH,C3BOND,C4BOND,C3ANGL,
     *                NFFTYP,NFFFILE,LJQMMM,LJQM,INTCHG,
     *                LJSIGMA,JTOPFILE(90),JPARFILE(90),
     *                JTOPAMIA(90),JTOPNTER(90),JTOPCTER(90),
     *                JTOPNUCA(90),JPARFIL2(90),JPARFIL3(90)
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /IOFILE/ IR,IW,IP,IJK,IJKT,IDAF,NAV,IODA(950)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
C
      INTEGER, PARAMETER :: K15 = SELECTED_INT_KIND(15)
      INTEGER, PARAMETER :: MAGIC2 = TRANSFER(20170303193144_K15,1)
C
C     NANDUN THELLAMUREGE, HUI LI, JAN 2011, LINCOLN
C     FENGCHAO CUI, HUI LI, MAY 2012, LINCOLN
C     HONGBO ZHU, HUI LI, NOV 3, LINCOLN
C
C     FORMULA:  E = 4*EPSILON*[(SIGMA/R)**12 - (SIGMA/R)**6)]
C
C               SIGMA(I,J)   = 0.5*(SIGMA_I   + SIGMA_J  )
C                         OR = SQRT(SIGMA_I*SIGMA_J)
C               EPSILON(I,J) = SQRT(EPSILON_I * EPSILON_J)
C
      ENLJR  =ZERO
      ENLJD  =ZERO
      SOL1LJ =ZERO
      SOL2LJ =ZERO
      PMF1LJ =ZERO
      IF(IDOLJ.EQ.0) RETURN
C
      IF(NDFS.EQ.MAGIC2)THEN
         CALL ELJ126CCS(CORD,FFGRD,SIG,EPS,
     *                  NONLS1,L1213J,FCLJTP,NTYPE)
         RETURN
      END IF
C
      LLTODO = 2
      IF(IFEPTYP.GT.0.AND.IDOPOL.EQ.0.AND.IEWALD.EQ.0) THEN
         IF(    NAT.LE.0) THEN       ! FFMD1
            IF(IFEPTYP.EQ.1)LLTODO= 8
            IF(IFEPTYP.EQ.2)LLTODO=14 ! NO (9,10)
         END IF
      END IF
      DO 110 LL=1,LLTODO
      IF(IFEPTYP.EQ.2.AND.NAT.LE.0)THEN
         IF(LL.GE.9.AND.LL.LE.10) GOTO 110
      END IF
      IF(LL.EQ.1) THEN
         NN1  = 1
         NN2  = NTODO
         SIGN = 1.0D+00
         DWT1 = SIGN
         DWT2 = SIGN
      END IF
      IF(LL.EQ.2) THEN
         NN1  = L11213
         NN2  = L21213
         SIGN = -1.0D+00
         DWT1 = SIGN
         DWT2 = SIGN
      END IF
      IF(LL.EQ.3) THEN
         NN1  = 1
         NN2  = NTODOA            ! REMOVE A SOL
         SIGN = -1.0D+00
         DWT1 = SIGN
         DWT2 = SIGN
      END IF
      IF(LL.EQ.4) THEN
         NN1  = L11213A
         NN2  = L21213A
         SIGN = 1.0D+00
         DWT1 = SIGN
         DWT2 = SIGN
      END IF
      IF(LL.EQ.5) THEN
         NN1  = 1
         NN2  = NTODOA
         SIGN = 1 - WSIMUL        ! ADD (1-WSIMUL) A SOL
         DWT1 = 1 - WPERT1
         DWT2 = 1 - WPERT2
      END IF
      IF(LL.EQ.6) THEN
         NN1  = L11213A
         NN2  = L21213A
         SIGN = -1 + WSIMUL
         DWT1 = -1 + WPERT1
         DWT2 = -1 + WPERT2
      END IF
      IF(LL.EQ.7) THEN
         NN1  = 1
         NN2  = NTODOB
         SIGN = WSIMUL            ! ADD WSIMUL B SOL
         DWT1 = WPERT1
         DWT2 = WPERT2
      END IF
      IF(LL.EQ.8) THEN
         NN1  = L11213B
         NN2  = L21213B
         SIGN = -WSIMUL
         DWT1 = -WPERT1
         DWT2 = -WPERT2
      END IF
      IF(LL.EQ.9) THEN
         NN1  = 1
         NN2  = NTODOB            ! REMOVE B SOL
         SIGN = -1.0D+00
         DWT1 = SIGN
         DWT2 = SIGN
      END IF
      IF(LL.EQ.10) THEN
         NN1  = L11213B
         NN2  = L21213B
         SIGN = 1.0D+00
         DWT1 = SIGN
         DWT2 = SIGN
      END IF
      IF(LL.EQ.11) THEN
         NN1  = 1
         NN2  = NTODOPMA
         SIGN = -WSIMUL
         DWT1 = -WPERT1     !  NEGATIVE A
         DWT2 = -WPERT2
      END IF
      IF(LL.EQ.12) THEN
         NN1  = L11213PMA
         NN2  = L21213PMA
         SIGN = +WSIMUL
         DWT1 = +WPERT1
         DWT2 = +WPERT2
      END IF
      IF(LL.EQ.13) THEN
         NN1  = 1
         NN2  = NTODOPMB
         SIGN = +WSIMUL
         DWT1 = +WPERT1     !  POSITIVE B
         DWT2 = +WPERT2
      END IF
      IF(LL.EQ.14) THEN
         NN1  = L11213PMB
         NN2  = L21213PMB
         SIGN = -WSIMUL
         DWT1 = -WPERT1
         DWT2 = -WPERT2
      END IF
      DO 100 III=NN1, NN2
         IF(LL.EQ.1) THEN
            I = NONLS1(1,III)
            J = NONLS1(2,III)
         END IF
         IF(LL.EQ.2) THEN
            I = L1213J(1,III)
            J = L1213J(2,III)
         END IF
         IF(LL.EQ.3) THEN
            I = NONLSA(1,III)
            J = NONLSA(2,III)
         END IF
         IF(LL.EQ.4) THEN
            I = L1213A(1,III)
            J = L1213A(2,III)
         END IF
         IF(LL.EQ.5) THEN
            I = NONLSA(1,III)
            J = NONLSA(2,III)
         END IF
         IF(LL.EQ.6) THEN
            I = L1213A(1,III)
            J = L1213A(2,III)
         END IF
         IF(LL.EQ.7) THEN
            I = NONLSB(1,III)
            J = NONLSB(2,III)
         END IF
         IF(LL.EQ.8) THEN
            I = L1213B(1,III)
            J = L1213B(2,III)
         END IF
         IF(LL.EQ.9) THEN
            I = NONLSB(1,III)
            J = NONLSB(2,III)
         END IF
         IF(LL.EQ.10) THEN
            I = L1213B(1,III)
            J = L1213B(2,III)
         END IF
         IF(LL.EQ.11) THEN
            I = NONLSPMA(1,III)
            J = NONLSPMA(2,III)
         END IF
         IF(LL.EQ.12) THEN
            I = L1213PMA(1,III)
            J = L1213PMA(2,III)
         END IF
         IF(LL.EQ.13) THEN
            I = NONLSPMB(1,III)
            J = NONLSPMB(2,III)
         END IF
         IF(LL.EQ.14) THEN
            I = L1213PMB(1,III)
            J = L1213PMB(2,III)
         END IF
         IF(I.EQ.0.OR.J.EQ.0) GOTO 100
C
         IF(NFFTYP/10000.NE.5) THEN
            SIGI  = SIG(I)
            EPSI  = EPS(I)
            SIGJ  = SIG(J)
            EPSJ  = EPS(J)
            IF(IFEPTYP.GT.0.AND.(LL.EQ. 7.OR.LL.EQ. 8.OR.
     *                           LL.EQ.13.OR.LL.EQ.14)) THEN
               SIGI  = SIGB(I)
               EPSI  = EPSB(I)
               SIGJ  = SIGB(J)
               EPSJ  = EPSB(J)
            END IF
            IF(SIGI.EQ.ZERO .OR. EPSI.EQ.ZERO) GOTO 100
            IF(SIGJ.EQ.ZERO .OR. EPSJ.EQ.ZERO) GOTO 100
         END IF
         IF(NFFTYP/10000.EQ.5) THEN
            II    = NTYPE(I)
            JJ    = NTYPE(J)
            RRIJ  = FCLJTP(1,JJ,II)
            EPSNO = FCLJTP(2,JJ,II)
            IF(RRIJ.EQ.ZERO.OR.EPSNO.EQ.ZERO) GOTO 100
            RRIJ2 = RRIJ*RRIJ
            RRIJ6 = RRIJ2*RRIJ2*RRIJ2
            RRIJ7 = RRIJ6*RRIJ
         END IF
C
         X     = CORD(1,I) - CORD(1,J)
         Y     = CORD(2,I) - CORD(2,J)
         Z     = CORD(3,I) - CORD(3,J)
         IF(IFEPTYP.GT.0.AND.(LL.EQ. 7.OR.LL.EQ. 8.OR.
     *                        LL.EQ.13.OR.LL.EQ.14)) THEN
            DO IFIXB = 1, NFIXMM
               IF(I.EQ.IFIXMM(IFIXB)) THEN
                  X = X + CORDB(1,I) - CORD(1,I)
                  Y = Y + CORDB(2,I) - CORD(2,I)
                  Z = Z + CORDB(3,I) - CORD(3,I)
               END IF
               IF(J.EQ.IFIXMM(IFIXB)) THEN
                  X = X - CORDB(1,J) + CORD(1,J)
                  Y = Y - CORDB(2,J) + CORD(2,J)
                  Z = Z - CORDB(3,J) + CORD(3,J)
               END IF
            ENDDO
         END IF
         PBCX  = XBOX * ANINT(X*ONEXBOX)
         PBCY  = YBOX * ANINT(Y*ONEYBOX)
         PBCZ  = ZBOX * ANINT(Z*ONEZBOX)
         X     = X - PBCX
         Y     = Y - PBCY
         Z     = Z - PBCZ
         R2    = X*X+Y*Y+Z*Z
         IF(R2.GT.SWRB2) GOTO 100
         IF(R2.LT.0.01D+00) GOTO 100
         CALL SWFUNC(R2,X,Y,Z)
         ONER2 = ONE/R2
C
         IF(NFFTYP/10000.NE.5) THEN
            SIG1  = 0.5D+00*(SIGI+SIGJ)
            FOUREP= EPSI*EPSJ
            SIGSQ = SIG1*SIG1*ONER2
            IF(NFFTYP/10000.EQ.4) SIGSQ  = SIGI*SIGJ*ONER2
            DISP  = -SIGSQ*SIGSQ*SIGSQ
            REP   = DISP*DISP
               EPAIRR  = FOUREP*REP
               EPAIRD  = FOUREP*DISP
               EPAIRR  = EPAIRR*SIGN
               EPAIRD  = EPAIRD*SIGN
               DUM     =-SWF*6.0D+00*FOUREP*(2.0D+00*REP+DISP)*
     *                   ONER2*SIGN
               IF(LL.GT.4.AND.LL.NE.9.AND.LL.NE.10) THEN
                  EPAIRRW = FOUREP*REP
                  EPAIRDW = FOUREP*DISP
                  EPAIRRW1= EPAIRRW*DWT1
                  EPAIRDW1= EPAIRDW*DWT1
                  EPAIRRW2= EPAIRRW*DWT2
                  EPAIRDW2= EPAIRDW*DWT2
               END IF
         END IF
         IF(NFFTYP/10000.EQ.5) THEN
            R       = SQRT(R2)
            ONER    = ONE/R
            R6      = R2*R2*R2
            R7      = R6*R
            C1LJ    = ONE/(R+0.07D+00*RRIJ)
            C2LJ    = ONE/(R7+0.12D+00*RRIJ7)
            FM1VDW  = 1.07D+00*RRIJ*C1LJ
            FM1VDW2 = FM1VDW*FM1VDW
            FM1VDW6 = FM1VDW2*FM1VDW2*FM1VDW2
            FM1VDW7 = FM1VDW6*FM1VDW
            FM2VDW  = 1.12D+00*RRIJ7*C2LJ
C
            EPAIRR   = EPSNO*FM1VDW7*FM2VDW
            EPAIRD   =-EPSNO*FM1VDW7*TWO
            EPAIRR   = EPAIRR*SIGN
            EPAIRD   = EPAIRD*SIGN
            DUM      = -SWF*EPSNO*SEVEN*FM1VDW7*
     *                 ( (FM2VDW-TWO)*C1LJ + FM2VDW*C2LJ*R6 )
     *                 *ONER*SIGN
            IF(LL.GT.4.AND.LL.NE.9.AND.LL.NE.10) THEN
               EPAIRRW = EPSNO*FM1VDW7*FM2VDW
               EPAIRDW =-EPSNO*FM1VDW7*TWO
               EPAIRRW1= EPAIRRW*DWT1
               EPAIRDW1= EPAIRDW*DWT1
               EPAIRRW2= EPAIRRW*DWT2
               EPAIRDW2= EPAIRDW*DWT2
            END IF
         END IF
C
         EPAIR = EPAIRR+EPAIRD
         ENLJR = ENLJR + EPAIRR*SWF
         ENLJD = ENLJD + EPAIRD*SWF
         IF(LL.GE. 5.AND.LL.LE. 8) THEN
            SOL1LJ= SOL1LJ + (EPAIRRW1+EPAIRDW1-EPAIRR-EPAIRD)*SWF
            SOL2LJ= SOL2LJ + (EPAIRRW2+EPAIRDW2-EPAIRR-EPAIRD)*SWF
         END IF
         IF(LL.GE.11.AND.LL.LE.14) THEN
            PMF1LJ= PMF1LJ + (EPAIRRW1+EPAIRDW1-EPAIRR-EPAIRD)*SWF
         END IF
C
         DEX   = DUM*X + EPAIR*SWFDX
         DEY   = DUM*Y + EPAIR*SWFDY
         DEZ   = DUM*Z + EPAIR*SWFDZ
         FFGRD(1,I)=FFGRD(1,I) + DEX
         FFGRD(2,I)=FFGRD(2,I) + DEY
         FFGRD(3,I)=FFGRD(3,I) + DEZ
         FFGRD(1,J)=FFGRD(1,J) - DEX
         FFGRD(2,J)=FFGRD(2,J) - DEY
         FFGRD(3,J)=FFGRD(3,J) - DEZ
         VIR(1)    =VIR(1)     + DEX*X
         VIR(2)    =VIR(2)     + DEY*Y
         VIR(3)    =VIR(3)     + DEZ*Z
         IYES = 0
         JYES = 0
         DO KFIX=1,NFIXMM
            IF(I.EQ.IFIXMM(KFIX)) IYES = 1
            IF(J.EQ.IFIXMM(KFIX)) JYES = 1
         ENDDO
         IF(IYES.EQ.1.AND.JYES.EQ.1) THEN
            VIR(1)    =VIR(1)     - DEX*X
            VIR(2)    =VIR(2)     - DEY*Y
            VIR(3)    =VIR(3)     - DEZ*Z
         END IF
 100  CONTINUE
 110  CONTINUE
C
C
C     -- REMOVE ALL 1-4 PAIRS, THEN ADD THEM BACK.
C        WHEN ADD BACK, SELECT THE CORRECT
C        L-J PARAMETERS. -HL.
C
      IF(WT14LJ.EQ.1.0D+00.AND.NFFTYP/10000.NE.5) THEN
      LLTODO = 1
      IF(IFEPTYP.GT.0.AND.IDOPOL.EQ.0) THEN
         IF(    NAT.LE.0) THEN       ! FFMD1
            IF(IFEPTYP.EQ.1)LLTODO= 4
            IF(IFEPTYP.EQ.2)LLTODO= 7
         END IF
      END IF
      DO 390 LL=1,LLTODO
      IF(IFEPTYP.EQ.2.AND.NAT.LE.0)THEN
         IF(LL.EQ.5) GOTO 390
      END IF
      IF(LL.EQ.1) THEN
         NN1  = L1N14J
         NN2  = L2N14J
         SIGN = 1.0D+00
         DWT1 = ZERO
         DWT2 = ZERO
      END IF
      IF(LL.EQ.2) THEN
         NN1  = L1N14A
         NN2  = L2N14A
         SIGN = -1.0D+00
         DWT1 = ZERO
         DWT2 = ZERO
      END IF
      IF(LL.EQ.3) THEN
         NN1  = L1N14A
         NN2  = L2N14A
         SIGN = 1 - WSIMUL
         DWT1 = 1 - WPERT1
         DWT2 = 1 - WPERT2
      END IF
      IF(LL.EQ.4) THEN
         NN1  = L1N14B
         NN2  = L2N14B
         SIGN = WSIMUL
         DWT1 = WPERT1
         DWT2 = WPERT2
      END IF
      IF(LL.EQ.5) THEN
         NN1  = L1N14B
         NN2  = L2N14B
         SIGN = -1.0D+00
         DWT1 = ZERO
         DWT2 = ZERO
      END IF
      IF(LL.EQ.6) THEN
         NN1  = L1N14PMA
         NN2  = L2N14PMA
         SIGN = -WSIMUL
         DWT1 = -WPERT1
         DWT2 = -WPERT2
      END IF
      IF(LL.EQ.7) THEN
         NN1  = L1N14PMB
         NN2  = L2N14PMB
         SIGN = +WSIMUL
         DWT1 = +WPERT1
         DWT2 = +WPERT2
      END IF
      DO 400 III=NN1, NN2
         IF(LL.EQ.1) THEN
            I = L14J(1,III)
            J = L14J(2,III)
         END IF
         IF(LL.EQ.2) THEN
            I = L14A(1,III)
            J = L14A(2,III)
         END IF
         IF(LL.EQ.3) THEN
            I = L14A(1,III)
            J = L14A(2,III)
         END IF
         IF(LL.EQ.4) THEN
            I = L14B(1,III)
            J = L14B(2,III)
         END IF
         IF(LL.EQ.5) THEN
            I = L14B(1,III)
            J = L14B(2,III)
         END IF
         IF(LL.EQ.6) THEN
            I = L14PMA(1,III)
            J = L14PMA(2,III)
         END IF
         IF(LL.EQ.7) THEN
            I = L14PMB(1,III)
            J = L14PMB(2,III)
         END IF
         IF(I.LE.0.OR.J.LE.0) GOTO 400
         SIG2I  = SIG2(I)
         EPS2I  = EPS2(I)
         SIG2J  = SIG2(J)
         EPS2J  = EPS2(J)
         IF(IFEPTYP.GT.0.AND.(LL.EQ.4.OR.LL.EQ.7))THEN
            SIG2I = SIG2B(I)
            EPS2I = EPS2B(I)
            SIG2J = SIG2B(J)
            EPS2J = EPS2B(J)
         END IF
         IF(SIG2I.EQ.ZERO.AND.EPS2I.EQ.ZERO .AND.
     *      SIG2J.EQ.ZERO.AND.EPS2J.EQ.ZERO) GOTO 400
         DO 410 LLL = 1, 2
            IF(LLL.EQ.1) THEN
               FACT = -1.0D+00     ! REMOVE
               SIGI  = SIG(I)
               EPSI  = EPS(I)
               SIGJ  = SIG(J)
               EPSJ  = EPS(J)
               IF(IFEPTYP.GT.0.AND.(LL.EQ.4.OR.LL.EQ.7))THEN
                  SIGI  = SIGB(I)
                  EPSI  = EPSB(I)
                  SIGJ  = SIGB(J)
                  EPSJ  = EPSB(J)
               END IF
            END IF
            IF(LLL.EQ.2) THEN
               FACT =  1.0D+00     ! ADD
               SIGI  = SIG(I)
               EPSI  = EPS(I)
               SIGJ  = SIG(J)
               EPSJ  = EPS(J)
               IF(SIG2I.GT.ZERO.AND.EPS2I.GT.ZERO) THEN
                  SIGI = SIG2I
                  EPSI = EPS2I
               END IF
               IF(SIG2J.GT.ZERO.AND.EPS2J.GT.ZERO) THEN
                  SIGJ = SIG2J
                  EPSJ = EPS2J
               END IF
            END IF
            IF(SIGI.EQ.ZERO .OR. EPSI.EQ.ZERO) GOTO 410
            IF(SIGJ.EQ.ZERO .OR. EPSJ.EQ.ZERO) GOTO 410
            X     = CORD(1,I) - CORD(1,J)
            Y     = CORD(2,I) - CORD(2,J)
            Z     = CORD(3,I) - CORD(3,J)
            IF(IFEPTYP.GT.0.AND.(LL.EQ.4.OR.LL.EQ.7))THEN
               DO IFIXB = 1, NFIXMM
                  IF(I.EQ.IFIXMM(IFIXB)) THEN
                     X = X + CORDB(1,I) - CORD(1,I)
                     Y = Y + CORDB(2,I) - CORD(2,I)
                     Z = Z + CORDB(3,I) - CORD(3,I)
                  END IF
                  IF(J.EQ.IFIXMM(IFIXB)) THEN
                     X = X - CORDB(1,J) + CORD(1,J)
                     Y = Y - CORDB(2,J) + CORD(2,J)
                     Z = Z - CORDB(3,J) + CORD(3,J)
                  END IF
               ENDDO
            END IF
            R2    = X*X+Y*Y+Z*Z
            CALL SWFUNC(R2,X,Y,Z)
            ONER2 = ONE/R2
C
            SIG1  = 0.5D+00*(SIGI+SIGJ)
            FOUREP= EPSI*EPSJ
            SIGSQ = SIG1*SIG1*ONER2
            IF(NFFTYP/10000.EQ.4) SIGSQ  = SIGI*SIGJ*ONER2
            DISP  = -SIGSQ*SIGSQ*SIGSQ
            REP   = DISP*DISP
               EPAIRR  = FOUREP*REP
               EPAIRD  = FOUREP*DISP
               EPAIRR  = EPAIRR*FACT*SIGN
               EPAIRD  = EPAIRD*FACT*SIGN
               DUM     =-SWF*6.0D+00*FOUREP*(2.0D+00*REP
     *                                      +DISP)*ONER2*FACT*SIGN
               IF(LL.GT.2.AND.LL.NE.5) THEN
                  EPAIRRW = FOUREP*REP*FACT
                  EPAIRDW = FOUREP*DISP*FACT
                  EPAIRRW1= EPAIRRW*DWT1
                  EPAIRDW1= EPAIRDW*DWT1
                  EPAIRRW2= EPAIRRW*DWT2
                  EPAIRDW2= EPAIRDW*DWT2
               END IF
C
            EPAIR = EPAIRR+EPAIRD
            ENLJR = ENLJR + EPAIRR*SWF
            ENLJD = ENLJD + EPAIRD*SWF
            IF(LL.GE. 3.AND.LL.LE. 4) THEN
               SOL1LJ= SOL1LJ + (EPAIRRW1+EPAIRDW1-EPAIRR-EPAIRD)*SWF
               SOL2LJ= SOL2LJ + (EPAIRRW2+EPAIRDW2-EPAIRR-EPAIRD)*SWF
            END IF
            IF(LL.GE. 6.AND.LL.LE. 7) THEN
               PMF1LJ= PMF1LJ + (EPAIRRW1+EPAIRDW1-EPAIRR-EPAIRD)*SWF
            END IF
            DEX   = DUM*X + EPAIR*SWFDX
            DEY   = DUM*Y + EPAIR*SWFDY
            DEZ   = DUM*Z + EPAIR*SWFDZ
            FFGRD(1,I)=FFGRD(1,I) + DEX
            FFGRD(2,I)=FFGRD(2,I) + DEY
            FFGRD(3,I)=FFGRD(3,I) + DEZ
            FFGRD(1,J)=FFGRD(1,J) - DEX
            FFGRD(2,J)=FFGRD(2,J) - DEY
            FFGRD(3,J)=FFGRD(3,J) - DEZ
            VIR(1)    =VIR(1)     + DEX*X
            VIR(2)    =VIR(2)     + DEY*Y
            VIR(3)    =VIR(3)     + DEZ*Z
            IYES = 0
            JYES = 0
            DO KFIX=1,NFIXMM
               IF(I.EQ.IFIXMM(KFIX)) IYES = 1
               IF(J.EQ.IFIXMM(KFIX)) JYES = 1
            ENDDO
            IF(IYES.EQ.1.AND.JYES.EQ.1) THEN
               VIR(1)    =VIR(1)     - DEX*X
               VIR(2)    =VIR(2)     - DEY*Y
               VIR(3)    =VIR(3)     - DEZ*Z
            END IF
 410     CONTINUE
 400  CONTINUE
 390  CONTINUE
      END IF
C
C     -- IF WT14LJ IS NOT 1.0, ADD (WT14LJ-1)
C
      IF(WT14LJ.NE.1.0D+00.AND.NFFTYP/10000.NE.5) THEN
         FACT   =  WT14LJ - 1.0D+00
         LLTODO = 1
         IF(IFEPTYP.GT.0.AND.IDOPOL.EQ.0) THEN
            IF(    NAT.LE.0) THEN       ! FFMD1
               IF(IFEPTYP.EQ.1)LLTODO= 4
               IF(IFEPTYP.EQ.2)LLTODO= 7
            END IF
         END IF
         DO 490 LL=1,LLTODO
         IF(IFEPTYP.EQ.2.AND.NAT.LE.0)THEN
            IF(LL.EQ.5) GOTO 490
         END IF
         IF(LL.EQ.1) THEN
            NN1  = L1N14J
            NN2  = L2N14J
            SIGN = 1.0D+00
            DWT1 = ZERO
            DWT2 = ZERO
         END IF
         IF(LL.EQ.2) THEN
            NN1  = L1N14A
            NN2  = L2N14A
            SIGN = -1.0D+00
            DWT1 = ZERO
            DWT2 = ZERO
         END IF
         IF(LL.EQ.3) THEN
            NN1  = L1N14A
            NN2  = L2N14A
            SIGN = 1 - WSIMUL
            DWT1 = 1 - WPERT1
            DWT2 = 1 - WPERT2
         END IF
         IF(LL.EQ.4) THEN
            NN1  = L1N14B
            NN2  = L2N14B
            SIGN = WSIMUL
            DWT1 = WPERT1
            DWT2 = WPERT2
         END IF
         IF(LL.EQ.5) THEN
            NN1  = L1N14B
            NN2  = L2N14B
            SIGN = -1.0D+00
            DWT1 = ZERO
            DWT2 = ZERO
         END IF
         IF(LL.EQ.6) THEN
            NN1  = L1N14PMA
            NN2  = L2N14PMA
            SIGN = -WSIMUL
            DWT1 = -WPERT1
            DWT2 = -WPERT2
         END IF
         IF(LL.EQ.7) THEN
            NN1  = L1N14PMB
            NN2  = L2N14PMB
            SIGN = +WSIMUL
            DWT1 = +WPERT1
            DWT2 = +WPERT2
         END IF
         DO 500 III=NN1, NN2
            IF(LL.EQ.1) THEN
               I = L14J(1,III)
               J = L14J(2,III)
            END IF
            IF(LL.EQ.2) THEN
               I = L14A(1,III)
               J = L14A(2,III)
            END IF
            IF(LL.EQ.3) THEN
               I = L14A(1,III)
               J = L14A(2,III)
            END IF
            IF(LL.EQ.4) THEN
               I = L14B(1,III)
               J = L14B(2,III)
            END IF
            IF(LL.EQ.5) THEN
               I = L14B(1,III)
               J = L14B(2,III)
            END IF
            IF(LL.EQ.6) THEN
               I = L14PMA(1,III)
               J = L14PMA(2,III)
            END IF
            IF(LL.EQ.7) THEN
               I = L14PMB(1,III)
               J = L14PMB(2,III)
            END IF
            SIGI  = SIG(I)
            EPSI  = EPS(I)
            SIGJ  = SIG(J)
            EPSJ  = EPS(J)
            IF(IFEPTYP.GT.0.AND.(LL.EQ.4.OR.LL.EQ.7))THEN
               SIGI = SIGB(I)
               EPSI = EPSB(I)
               SIGJ = SIGB(J)
               EPSJ = EPSB(J)
            END IF
            IF(SIGI.EQ.ZERO .OR. EPSI.EQ.ZERO) GOTO 500
            IF(SIGJ.EQ.ZERO .OR. EPSJ.EQ.ZERO) GOTO 500
            X     = CORD(1,I) - CORD(1,J)
            Y     = CORD(2,I) - CORD(2,J)
            Z     = CORD(3,I) - CORD(3,J)
            IF(IFEPTYP.GT.0.AND.(LL.EQ.4.OR.LL.EQ.7))THEN
               DO IFIXB = 1, NFIXMM
                  IF(I.EQ.IFIXMM(IFIXB)) THEN
                     X = X + CORDB(1,I) - CORD(1,I)
                     Y = Y + CORDB(2,I) - CORD(2,I)
                     Z = Z + CORDB(3,I) - CORD(3,I)
                  END IF
                  IF(J.EQ.IFIXMM(IFIXB)) THEN
                     X = X - CORDB(1,J) + CORD(1,J)
                     Y = Y - CORDB(2,J) + CORD(2,J)
                     Z = Z - CORDB(3,J) + CORD(3,J)
                  END IF
               ENDDO
            END IF
            R2    = X*X+Y*Y+Z*Z
            CALL SWFUNC(R2,X,Y,Z)
            ONER2 = ONE/R2
C
            SIG1  = 0.5D+00*(SIGI+SIGJ)
            FOUREP= EPSI*EPSJ
            SIGSQ = SIG1*SIG1*ONER2
            IF(NFFTYP/10000.EQ.4) SIGSQ  = SIGI*SIGJ*ONER2
            DISP  = -SIGSQ*SIGSQ*SIGSQ
            REP   = DISP*DISP
C
               EPAIRR  = FOUREP*REP
               EPAIRD  = FOUREP*DISP
               EPAIRR  = EPAIRR*FACT*SIGN
               EPAIRD  = EPAIRD*FACT*SIGN
               DUM     =-SWF*6.0D+00*FOUREP*(2.0D+00*REP
     *                                      +DISP)*ONER2*FACT*SIGN
               IF(LL.GT.2.AND.LL.NE.5) THEN
                  EPAIRRW = FOUREP*REP *FACT
                  EPAIRDW = FOUREP*DISP*FACT
                  EPAIRRW1= EPAIRRW*DWT1
                  EPAIRDW1= EPAIRDW*DWT1
                  EPAIRRW2= EPAIRRW*DWT2
                  EPAIRDW2= EPAIRDW*DWT2
               END IF
C
            EPAIR = EPAIRR+EPAIRD
            ENLJR = ENLJR + EPAIRR*SWF
            ENLJD = ENLJD + EPAIRD*SWF
            IF(LL.GE. 3.AND.LL.LE. 4) THEN
               SOL1LJ= SOL1LJ + (EPAIRRW1+EPAIRDW1-EPAIRR-EPAIRD)*SWF
               SOL2LJ= SOL2LJ + (EPAIRRW2+EPAIRDW2-EPAIRR-EPAIRD)*SWF
            END IF
            IF(LL.GE. 6.AND.LL.LE. 7) THEN
               PMF1LJ= PMF1LJ + (EPAIRRW1+EPAIRDW1-EPAIRR-EPAIRD)*SWF
            END IF
            DEX   = DUM*X + EPAIR*SWFDX
            DEY   = DUM*Y + EPAIR*SWFDY
            DEZ   = DUM*Z + EPAIR*SWFDZ
            FFGRD(1,I)=FFGRD(1,I) + DEX
            FFGRD(2,I)=FFGRD(2,I) + DEY
            FFGRD(3,I)=FFGRD(3,I) + DEZ
            FFGRD(1,J)=FFGRD(1,J) - DEX
            FFGRD(2,J)=FFGRD(2,J) - DEY
            FFGRD(3,J)=FFGRD(3,J) - DEZ
            VIR(1)    =VIR(1)     + DEX*X
            VIR(2)    =VIR(2)     + DEY*Y
            VIR(3)    =VIR(3)     + DEZ*Z
            IYES = 0
            JYES = 0
            DO KFIX=1,NFIXMM
               IF(I.EQ.IFIXMM(KFIX)) IYES = 1
               IF(J.EQ.IFIXMM(KFIX)) JYES = 1
            ENDDO
            IF(IYES.EQ.1.AND.JYES.EQ.1) THEN
               VIR(1)    =VIR(1)     - DEX*X
               VIR(2)    =VIR(2)     - DEY*Y
               VIR(3)    =VIR(3)     - DEZ*Z
            END IF
 500     CONTINUE
 490     CONTINUE
      END IF
C
      RETURN
      END
C*MODULE QUANPOB  *DECK ESPHER
!>
!> @brief    calculate LJ energy for MM atoms due to sphere
!>
      SUBROUTINE ESPHER(CORD,FFGRD)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      LOGICAL GOPARR,DSKWRK,MASWRK
C
      PARAMETER (ZERO=0.0D+00)
      PARAMETER (ONE=1.0D+00)
      PARAMETER (TWO=2.0D+00)
      PARAMETER (FOUR=4.0D+00)
      PARAMETER (SIX=6.0D+00)
C
      DIMENSION CORD(3,*),FFGRD(3,*)
C
      COMMON /FFENGY/ EN12,EN123,EN123R4,EN123B4,EN234W1,ENCHAR,ENLJR,
     *                ENLJD,ENPOL,XENPOL,ENRXN,XENRXN,ENRXNPOL,ENRXNR,
     *                EN12312,ENQUANP(30),
     *                ENBIAS,ENCENT,ENUCCH,ENCMAP,ENPOT,ENKIN,ENTOT
      COMMON /FFMDPA/ DT,DT2,TEMP0,PRES0,POLTOL,VIR(3),PMEAN,VOLAV,
     *                ENPAV,ENKAV,TEMPAV,BERENDT,BERENDP,VELMAX,
     *                PMEANX,PMEANY,PMEANZ,NSTEP,KMASTER,KOUTACT(2),
     *                IHESS,INTALG,ITSTAT,IPSTAT,JOUT,KOUT,LOUT
      COMMON /FFPARM/ NFFAT,NBOND,NANGL,NDIHR,NDIHB,NCMAP,NWAGG,
     *                N1213J,N14J,NLKQMM,IDOCHG,IDOPOL,IDOLJ,IDOCMAP
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
      COMMON /FFRMSD/ DIPT(8),TIMGYRA,TIMRALL,
     *                NATPDB,NGYRA,NDIEL,NRALL,
     *                LFFRALL0,NRMSD,LFFRMSD0,KFREEAB(201),
     *                NRIJMM,IJRMM(2,100),NRIJQM,IJRQM(2,100),
     *                NAIJKMM,IJKMM(3,100),NAIJKQM,IJKQM(3,100),
     *                NFIXMM,IFIXMM(200),NFIXQM,IFIXQM(200)
      COMMON /FFSPH / SPHRAD,SPHSIG,SPHEPS,IADDWAT
      COMMON /FMCOM / XX(1)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
C
C     HUI LI, FEB 2011, LINCOLN
C
C     --  ENRXNR: LENNARD-JONES POTENTIAL ENERGY DUE TO SPHERE --
C                 (REQUIRED IF SPHRAD IS APPLIED)
C
      ENRXNR =ZERO
      IF(SPHRAD.GE.1.0D+30.OR.SPHEPS.EQ.ZERO) RETURN
C
      R     = SPHRAD + (TWO**(ONE/SIX)-ONE)*SPHSIG
      FOUREP= FOUR*SPHEPS
      SIG6  = SPHSIG**6
      SIG12 = SIG6*SIG6
C
      IPCOUNT = ME - 1
      DO 100 IFFAT = 1, NFFAT
         IF(GOPARR) THEN
            IPCOUNT = IPCOUNT + 1
            IF(MOD(IPCOUNT,NPROC).NE.0) GOTO 100
         END IF
         IF(XX(LFFZANF+IFFAT-1).EQ.1.0D+00) GOTO 100
         XI    = CORD(1,IFFAT) - CENTX
         YI    = CORD(2,IFFAT) - CENTY
         ZI    = CORD(3,IFFAT) - CENTZ
         RI2   = XI*XI+YI*YI+ZI*ZI
         RI    = SQRT(RI2)
         R2    = (RI-R)**2
         ONER2 = ONE/R2
         ONER6 = ONER2*ONER2*ONER2
         ONER12= ONER6*ONER6
         REP   = SIG12*ONER12
         DISP  =-SIG6 *ONER6
         ENRXNR= ENRXNR + FOUREP*(REP+DISP)
         IF(R2.LT.(0.50D+00*SPHSIG)**2)ENRXNR=ENRXNR-FOUREP*(REP+DISP)
         IF(RI.EQ.ZERO) THEN
            DUM = ZERO
         ELSE
            DUM   = -SIX*FOUREP*(TWO*REP+DISP)/((RI-R)*RI)
         END IF
         IF(RI.GT.R) DUM = -DUM
         IF(DUM.GT. 0.01D+00) DUM= 0.01D+00
         IF(DUM.LT.-0.01D+00) DUM=-0.01D+00
         FFGRD(1,IFFAT)=FFGRD(1,IFFAT) + DUM*XI
         FFGRD(2,IFFAT)=FFGRD(2,IFFAT) + DUM*YI
         FFGRD(3,IFFAT)=FFGRD(3,IFFAT) + DUM*ZI
         VIR(1)        =VIR(1)         + DUM*XI*CORD(1,IFFAT)
         VIR(2)        =VIR(2)         + DUM*YI*CORD(2,IFFAT)
         VIR(3)        =VIR(3)         + DUM*ZI*CORD(3,IFFAT)
         IYES = 0
         DO KFIX=1,NFIXMM
            IF(IFFAT.EQ.IFIXMM(KFIX)) IYES = 1
         ENDDO
         IF(IYES.EQ.1) THEN
            VIR(1)    =VIR(1)     - DUM*XI*CORD(1,IFFAT)
            VIR(2)    =VIR(2)     - DUM*YI*CORD(2,IFFAT)
            VIR(3)    =VIR(3)     - DUM*ZI*CORD(3,IFFAT)
         END IF
 100  CONTINUE
C
      RETURN
      END
C*MODULE QUANPOB  *DECK ESPHQM
!>
!> @brief    calculate center point and sphere LJ energy for QM
!>
      SUBROUTINE ESPHQM(DETMP,LISTQM)
      use mx_limits, only: mxatm
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      LOGICAL GOPARR,DSKWRK,MASWRK
C
      PARAMETER (ZERO=0.0D+00)
      PARAMETER (ONE=1.0D+00)
      PARAMETER (TWO=2.0D+00)
      PARAMETER (FOUR=4.0D+00)
      PARAMETER (SIX=6.0D+00)
      PARAMETER (TOANGS=0.52917724924D+00)
      PARAMETER (TOBOHR=1.0D+00/TOANGS)
C
      DIMENSION DETMP(3,*),LISTQM(*)
C
      COMMON /FFENGY/ EN12,EN123,EN123R4,EN123B4,EN234W1,ENCHAR,ENLJR,
     *                ENLJD,ENPOL,XENPOL,ENRXN,XENRXN,ENRXNPOL,ENRXNR,
     *                EN12312,ENQUANP(30),
     *                ENBIAS,ENCENT,ENUCCH,ENCMAP,ENPOT,ENKIN,ENTOT
      COMMON /FFMDPA/ DT,DT2,TEMP0,PRES0,POLTOL,VIR(3),PMEAN,VOLAV,
     *                ENPAV,ENKAV,TEMPAV,BERENDT,BERENDP,VELMAX,
     *                PMEANX,PMEANY,PMEANZ,NSTEP,KMASTER,KOUTACT(2),
     *                IHESS,INTALG,ITSTAT,IPSTAT,JOUT,KOUT,LOUT
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
      COMMON /FFSPH / SPHRAD,SPHSIG,SPHEPS,IADDWAT
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
C
C     HUI LI, MAR 2011, LINCOLN
C
C     -- ENCENT: POTENTIAL ENERGY DUE TO A BOND-LIKE POTENTIAL --
C                (REQUIRED IN QM/MM SWITCHING FUNCTION METHOD)
C     -- ENRXNR: POTENTIAL ENERGY DUE TO SPHERE --
C                (REQUIRED IN REACTION FIELD METHOD)
C
      ENCENT = ZERO
      IF(MASWRK) THEN   ! ONLY ONE CASE, LET THE MASWRK DO IT
         IF(NAT.GT.0 .AND. NFFAT.GT.0 .AND. SWRB2.LT.1.0D+08) THEN
            XI  = C(1,LQMCT)-QMCX
            YI  = C(2,LQMCT)-QMCY
            ZI  = C(3,LQMCT)-QMCZ
            RI2 = XI*XI+YI*YI+ZI*ZI
            RI  = SQRT(RI2)
            R0  = 2.0D+00*TOBOHR         !   R0=2.0 ANGSTROM
            IF(RI.GT.R0) THEN
               ENCENT = ENCENT + (RI-R0)**2  ! K=1.0 HARTREE/BOHR**2
               DUM    = 2.0D+00*(RI-R0)/RI
               DETMP(1,LQMCT) = DETMP(1,LQMCT) + DUM*XI
               DETMP(2,LQMCT) = DETMP(2,LQMCT) + DUM*YI
               DETMP(3,LQMCT) = DETMP(3,LQMCT) + DUM*ZI
            END IF
         END IF
      END IF
C
      IF(SPHRAD.GE.1.0D+30.OR.SPHEPS.EQ.ZERO) RETURN
C
      R     = SPHRAD + (TWO**(ONE/SIX)-ONE)*SPHSIG
      FOUREP= FOUR*SPHEPS
      SIG6  = SPHSIG**6
      SIG12 = SIG6*SIG6
C
      IPCOUNT = ME - 1
      DO 200 IAT=1,NAT
       IF(LISTQM(NFFAT+IAT).EQ.0)THEN  ! DO NOT DO IT TWICE
         IF(GOPARR) THEN
            IPCOUNT = IPCOUNT + 1
            IF(MOD(IPCOUNT,NPROC).NE.0) GOTO 200
         END IF
         IF(IAN(IAT).EQ.1) GOTO 200
         XI    = C(1,IAT) - CENTX
         YI    = C(2,IAT) - CENTY
         ZI    = C(3,IAT) - CENTZ
         RI2   = XI*XI+YI*YI+ZI*ZI
         RI    = SQRT(RI2)
         R2    = (RI-R)**2
         ONER2 = ONE/R2
         ONER6 = ONER2*ONER2*ONER2
         ONER12= ONER6*ONER6
         REP   = SIG12*ONER12
         DISP  =-SIG6 *ONER6
         ENRXNR= ENRXNR + FOUREP*(REP+DISP)
         DUM   = -SIX*FOUREP*(TWO*REP+DISP)/((RI-R)*RI)
         IF(RI.GT.R) DUM = -DUM
         IF(DUM.GT. 0.01D+00) DUM= 0.01D+00
         IF(DUM.LT.-0.01D+00) DUM=-0.01D+00
         DETMP(1,IAT) = DETMP(1,IAT) + DUM*XI
         DETMP(2,IAT) = DETMP(2,IAT) + DUM*YI
         DETMP(3,IAT) = DETMP(3,IAT) + DUM*ZI
         VIR(1)       = VIR(1)       + DUM*XI*C(1,IAT)
         VIR(2)       = VIR(2)       + DUM*YI*C(2,IAT)
         VIR(3)       = VIR(3)       + DUM*ZI*C(3,IAT)
         IYES = 0
         DO KFIX=1,NFIXQM
            IF(IAT.EQ.IFIXQM(KFIX)) IYES = 1
         ENDDO
         IF(IYES.EQ.1) THEN
            VIR(1)    =VIR(1)     - DUM*XI*C(1,IAT)
            VIR(2)    =VIR(2)     - DUM*YI*C(2,IAT)
            VIR(3)    =VIR(3)     - DUM*ZI*C(3,IAT)
         END IF
       END IF
 200  CONTINUE
C
      RETURN
      END
C*MODULE QUANPOB  *DECK ECHARG
!>
!> @brief    charge-charge interaction
!>
!> @author   Hui Li
!>           - Jan 2011
!>
!> @details  force field charge-charge term
!>
      SUBROUTINE ECHARG(CORD,FFGRD,CHARG,NONLS1,L1213J,L14J,
     *                  CHARGB,NONLSA,NONLSB,L1213A,L1213B,L14A,L14B,
     *                  NONLSPMA,L1213PMA,L14PMA,CORDB,
     *                  NONLSPMB,L1213PMB,L14PMB)
      use mx_limits, only: mxatm
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      LOGICAL GOPARR,MASWRK,DSKWRK
C
      PARAMETER (ZERO=0.0D+00)
      PARAMETER (ONE=1.0D+00)
      PARAMETER (TWO=2.0D+00)
      PARAMETER (PI=3.14159265358979323846264338D+00)
      PARAMETER (TOANGS=0.52917724924D+00)
      PARAMETER (TOBOHR=1.0D+00/TOANGS)
C
      DIMENSION CORD(3,*),FFGRD(3,*),CHARG(*),NONLS1(2,*),
     *          L1213J(2,*),L14J(2,*),CHARGB(*),
     *          NONLSA(2,*),NONLSB(2,*),L1213A(2,*),L1213B(2,*),
     *          L14A(2,*),L14B(2,*),NONLSPMA(2,*),L1213PMA(2,*),
     *          L14PMA(2,*),CORDB(3,*),
     *          NONLSPMB(2,*),L1213PMB(2,*),L14PMB(2,*)
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
      COMMON /FFRMSD/ DIPT(8),TIMGYRA,TIMRALL,
     *                NATPDB,NGYRA,NDIEL,NRALL,
     *                LFFRALL0,NRMSD,LFFRMSD0,KFREEAB(201),
     *                NRIJMM,IJRMM(2,100),NRIJQM,IJRQM(2,100),
     *                NAIJKMM,IJKMM(3,100),NAIJKQM,IJKQM(3,100),
     *                NFIXMM,IFIXMM(200),NFIXQM,IFIXQM(200)
      COMMON /FFTYPE/ WT14LJ,WT14CH,C3BOND,C4BOND,C3ANGL,
     *                NFFTYP,NFFFILE,LJQMMM,LJQM,INTCHG,
     *                LJSIGMA,JTOPFILE(90),JPARFILE(90),
     *                JTOPAMIA(90),JTOPNTER(90),JTOPCTER(90),
     *                JTOPNUCA(90),JPARFIL2(90),JPARFIL3(90)
      COMMON /FMCOM / XX(1)
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /IOFILE/ IR,IW,IP,IJK,IJKT,IDAF,NAV,IODA(950)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
C
      INTEGER, PARAMETER :: K15 = SELECTED_INT_KIND(15)
      INTEGER, PARAMETER :: MAGIC2 = TRANSFER(20170303193144_K15,1)
C
C     HUI LI, JAN 2011, LINCOLN
C     FENGCHAO CUI, HUI LI, MAY 2012, LINCOLN
C     HUI LI, JUNE 19, 2012, ADD EWALD
C
C     FORMULA:  E = QI*QJ/RIJ  (IN VACUUM)
C          OR:  E = SHF*(QI*QJ/RIJ)
C                   SHF = SHIFTING  FUNCTION IN RANGE 0 - SWRB
C
CMMFF
CMMFF           E = QI*QJ/[D*(RIJ+0.05)]  !  0.05 ANGSTROM
CMMFF               D = DIELECTRIC CONSTANT, USUALLY 1.0
CMMFF
      ENCHAR    = ZERO
      SOL1CH    = ZERO
      SOL2CH    = ZERO
      PMF1CH    = ZERO
      IF(IDOCHG.EQ.0) RETURN
C
C     - THIS IS FOR MD-CCS SIMULATION -
      IF(ISHIFT.EQ.0.AND.NDFS.EQ.MAGIC2) THEN
         CALL ECHARGCCS(CORD,FFGRD,CHARG,L1213J,L14J)
         RETURN
      END IF
C
      SQRTPI=SQRT(PI)
C
C     -- MMFF94 BUFFER DISTANCE IS NOT USED FOR
C        SHIFTING FUNCTION, EWALD, IFIXSOL
C
      RBUF      = ZERO
      IF(NFFTYP/10000.EQ.5) RBUF = 0.05D+00*TOBOHR
      IF(ISHIFT .GT.0) RBUF = ZERO
      IF(IEWALD .GT.0) RBUF = ZERO
      IF(IFIXSOL.GT.0) RBUF = ZERO
C
      TWOROFF   = ZERO
      ONEROFF2  = ZERO
      ONEROFF   = ZERO
      ONEROFF3  = ZERO
      TWOROFF3  = ZERO
      TWOROFF2  = ZERO
      ONEROFF4  = ZERO
      THREEROFF4= ZERO
      IF(ISHIFT.EQ.1) THEN
         TWOROFF = ONESWRB + ONESWRB
         ONEROFF2= ONESWRB2
      END IF
      IF(ISHIFT.EQ.2) THEN
         ONEROFF3  = EPS1RB3
         TWOROFF3  = ONEROFF3 + ONEROFF3
         ONEROFF   = EPS1RB
      END IF
      IF(ISHIFT.EQ.3) THEN
         ONEROFF = ONESWRB
      END IF
      IF(ISHIFT.EQ.4) THEN
         TWOROFF2  = ONESWRB2 + ONESWRB2
         ONEROFF4  = ONESWRB4
         THREEROFF4= ONEROFF4 + ONEROFF4 + ONEROFF4
      END IF
      IF(ISHIFT.EQ.0.OR.IEWALD.GT.0) THEN
         TWOROFF   = ZERO
         ONEROFF2  = ZERO
         ONEROFF   = ZERO
         ONEROFF3  = ZERO
         TWOROFF3  = ZERO
         TWOROFF2  = ZERO
         ONEROFF4  = ZERO
         THREEROFF4= ZERO
      END IF
C
      LLTODO = 2
      IF(IFEPTYP.GT.0.AND.IDOPOL.EQ.0.AND.IEWALD.EQ.0) THEN
         IF(    NAT.LE.0) THEN       ! FFMD1
            IF(IFEPTYP.EQ.1)LLTODO= 8
            IF(IFEPTYP.EQ.2)LLTODO=14 ! NO (9,10)
         END IF
      END IF
      DO 110 LL=1,LLTODO
      IF(IFEPTYP.EQ.2.AND.NAT.LE.0)THEN
         IF(LL.GE.9.AND.LL.LE.10) GOTO 110
      END IF
      IF(LL.EQ.1) THEN
         NN1  = 1
         NN2  = NTODO
         SIGN = 1.0D+00           !  ALL I-J PAIRS
         DWT1 = ZERO
         DWT2 = ZERO
      END IF
      IF(LL.EQ.2) THEN
         NN1  = L11213
         NN2  = L21213
         SIGN = -1.0D+00          !  REMOVE 1-2, 1-3 PAIRS
         DWT1 = ZERO
         DWT2 = ZERO
      END IF
      IF(LL.EQ.3) THEN
         NN1  = 1
         NN2  = NTODOA
         SIGN = -1.0D+00          !  REMOVE A ALCHEMICAL ATOMS
         DWT1 = ZERO
         DWT2 = ZERO
      END IF
      IF(LL.EQ.4) THEN
         NN1  = L11213A
         NN2  = L21213A
         SIGN = 1.0D+00
         DWT1 = ZERO
         DWT2 = ZERO
      END IF
      IF(LL.EQ.5) THEN
         NN1  = 1
         NN2  = NTODOA
         SIGN = 1 - WSIMUL        !  ADD (1-WSIMUL) A
         DWT1 = 1 - WPERT1
         DWT2 = 1 - WPERT2
      END IF
      IF(LL.EQ.6) THEN
         NN1  = L11213A
         NN2  = L21213A
         SIGN = -1 + WSIMUL
         DWT1 = -1 + WPERT1
         DWT2 = -1 + WPERT2
      END IF
      IF(LL.EQ.7) THEN
         NN1  = 1
         NN2  = NTODOB
         SIGN = WSIMUL            ! ADD WSIMUL B
         DWT1 = WPERT1
         DWT2 = WPERT2
      END IF
      IF(LL.EQ.8) THEN
         NN1  = L11213B
         NN2  = L21213B
         SIGN = -WSIMUL
         DWT1 = -WPERT1
         DWT2 = -WPERT2
      END IF
      IF(LL.EQ.9) THEN
         NN1  = 1
         NN2  = NTODOB            !  REMOVE B ALCHEMICAL ATOMS
         SIGN = -1.0D+00
         DWT1 = ZERO
         DWT2 = ZERO
      END IF
      IF(LL.EQ.10) THEN
         NN1  = L11213B
         NN2  = L21213B
         SIGN = 1.0D+00
         DWT1 = ZERO
         DWT2 = ZERO
      END IF
      IF(LL.EQ.11) THEN
         NN1  = 1
         NN2  = NTODOPMA
         SIGN = -WSIMUL
         DWT1 = -WPERT1     !  NEGATIVE A
         DWT2 = -WPERT2
      END IF
      IF(LL.EQ.12) THEN
         NN1  = L11213PMA
         NN2  = L21213PMA
         SIGN = +WSIMUL
         DWT1 = +WPERT1
         DWT2 = +WPERT2
      END IF
      IF(LL.EQ.13) THEN
         NN1  = 1
         NN2  = NTODOPMB
         SIGN = +WSIMUL
         DWT1 = +WPERT1     !  POSITIVE B
         DWT2 = +WPERT2
      END IF
      IF(LL.EQ.14) THEN
         NN1  = L11213PMB
         NN2  = L21213PMB
         SIGN = -WSIMUL
         DWT1 = -WPERT1
         DWT2 = -WPERT2
      END IF
      DO 100 III=NN1, NN2
         IF(LL.EQ.1) THEN
            I = NONLS1(1,III)
            J = NONLS1(2,III)
         END IF
         IF(LL.EQ.2) THEN
            I = L1213J(1,III)
            J = L1213J(2,III)
         END IF
         IF(LL.EQ.3) THEN
            I = NONLSA(1,III)
            J = NONLSA(2,III)
         END IF
         IF(LL.EQ.4) THEN
            I = L1213A(1,III)
            J = L1213A(2,III)
         END IF
         IF(LL.EQ.5) THEN
            I = NONLSA(1,III)
            J = NONLSA(2,III)
         END IF
         IF(LL.EQ.6) THEN
            I = L1213A(1,III)
            J = L1213A(2,III)
         END IF
         IF(LL.EQ.7) THEN
            I = NONLSB(1,III)
            J = NONLSB(2,III)
         END IF
         IF(LL.EQ.8) THEN
            I = L1213B(1,III)
            J = L1213B(2,III)
         END IF
         IF(LL.EQ.9) THEN
            I = NONLSB(1,III)
            J = NONLSB(2,III)
         END IF
         IF(LL.EQ.10) THEN
            I = L1213B(1,III)
            J = L1213B(2,III)
         END IF
         IF(LL.EQ.11) THEN
            I = NONLSPMA(1,III)
            J = NONLSPMA(2,III)
         END IF
         IF(LL.EQ.12) THEN
            I = L1213PMA(1,III)
            J = L1213PMA(2,III)
         END IF
         IF(LL.EQ.13) THEN
            I = NONLSPMB(1,III)
            J = NONLSPMB(2,III)
         END IF
         IF(LL.EQ.14) THEN
            I = L1213PMB(1,III)
            J = L1213PMB(2,III)
         END IF
         IF(I.EQ.0.OR.J.EQ.0) GOTO 100
         QI  = CHARG(I)
         QJ  = CHARG(J)
         IF(IFEPTYP.GT.0.AND.(LL.EQ. 7.OR.LL.EQ. 8.OR.
     *                        LL.EQ.13.OR.LL.EQ.14)) THEN
            QI = CHARGB(I)
            QJ = CHARGB(J)
         END IF
         IF(QI.EQ.ZERO.OR.QJ.EQ.ZERO) GOTO 100
         QIQJ  = QI*QJ
C
         X     = CORD(1,I) - CORD(1,J)
         Y     = CORD(2,I) - CORD(2,J)
         Z     = CORD(3,I) - CORD(3,J)
         IF(IFEPTYP.GT.0.AND.(LL.EQ. 7.OR.LL.EQ. 8.OR.
     *                        LL.EQ.13.OR.LL.EQ.14)) THEN
            DO IFIXB = 1, NFIXMM
               IF(I.EQ.IFIXMM(IFIXB)) THEN
                  X = X + CORDB(1,I) - CORD(1,I)
                  Y = Y + CORDB(2,I) - CORD(2,I)
                  Z = Z + CORDB(3,I) - CORD(3,I)
               END IF
               IF(J.EQ.IFIXMM(IFIXB)) THEN
                  X = X - CORDB(1,J) + CORD(1,J)
                  Y = Y - CORDB(2,J) + CORD(2,J)
                  Z = Z - CORDB(3,J) + CORD(3,J)
               END IF
            ENDDO
         END IF
         PBCX  = XBOX * ANINT(X*ONEXBOX)
         PBCY  = YBOX * ANINT(Y*ONEYBOX)
         PBCZ  = ZBOX * ANINT(Z*ONEZBOX)
         X     = X - PBCX
         Y     = Y - PBCY
         Z     = Z - PBCZ
         R2    = X*X+Y*Y+Z*Z
         IF(R2.GT.SWRB2) GOTO 100
         IF(R2.LT.0.01D+00) GOTO 100
C
         R     = SQRT(R2)
         ONERX = ONE/R
         R     = R + RBUF
         R2    = R*R
         R3    = R*R2
         ONER  = ONE/R
         ONER2 = ONER*ONER
         ONER3 = ONER2*ONERX
C
            EDUM  = QIQJ*(ONER-TWOROFF+R*ONEROFF2
     *                           +R2*ONEROFF3-ONEROFF
     *                        -R*TWOROFF2+R3*ONEROFF4)*SIGN
            DUM   = -QIQJ*(ONER3-ONEROFF2*ONER-TWOROFF3
     *                +TWOROFF2*ONER-THREEROFF4*R)*SIGN
C           - EWALD CORRECTION APPLIES TO ALL PAIRS
C             (HUI LI: NO EWALD CORRECTION TO 1213 PAIRS)
            IF(IEWALD.GT.0.AND.LL.EQ.1) THEN
               F0      = GMSERFC(SPLIT*R)
               EDUM    = EDUM*F0
               EXPBR   = EXP(-SPLIT*SPLIT*R2)
               F0_GRAD = TWO*SPLIT*R*EXPBR/SQRTPI + F0
               DUM     = DUM*F0_GRAD
            END IF
            IF(LL.GT.4.AND.LL.NE.9.AND.LL.NE.10) THEN
               EDUMW  = QIQJ*(ONER-TWOROFF+R*ONEROFF2
     *                           +R2*ONEROFF3-ONEROFF
     *                        -R*TWOROFF2+R3*ONEROFF4)
               EDUMW1 = EDUMW*DWT1
               EDUMW2 = EDUMW*DWT2
            END IF
         EPAIR = EDUM
         ENCHAR= ENCHAR  + EPAIR
         IF(LL.GE. 5.AND.LL.LE. 8) THEN
            SOL1CH= SOL1CH + (EDUMW1-EPAIR)
            SOL2CH= SOL2CH + (EDUMW2-EPAIR)
         END IF
         IF(LL.GE.11.AND.LL.LE.14) THEN
            PMF1CH= PMF1CH + (EDUMW1-EPAIR)
         END IF
         DEX   = DUM*X
         DEY   = DUM*Y
         DEZ   = DUM*Z
         FFGRD(1,I)=FFGRD(1,I) + DEX
         FFGRD(2,I)=FFGRD(2,I) + DEY
         FFGRD(3,I)=FFGRD(3,I) + DEZ
         FFGRD(1,J)=FFGRD(1,J) - DEX
         FFGRD(2,J)=FFGRD(2,J) - DEY
         FFGRD(3,J)=FFGRD(3,J) - DEZ
         VIR(1)    =VIR(1)     + DEX*X
         VIR(2)    =VIR(2)     + DEY*Y
         VIR(3)    =VIR(3)     + DEZ*Z
         IYES = 0
         JYES = 0
         DO KFIX=1,NFIXMM
            IF(I.EQ.IFIXMM(KFIX)) IYES = 1
            IF(J.EQ.IFIXMM(KFIX)) JYES = 1
         ENDDO
         IF(IYES.EQ.1.AND.JYES.EQ.1) THEN
            VIR(1)    =VIR(1)     - DEX*X
            VIR(2)    =VIR(2)     - DEY*Y
            VIR(3)    =VIR(3)     - DEZ*Z
         END IF
 100  CONTINUE
 110  CONTINUE
C
C     -- EWALD SELF AND RECIPROCAL TERM --
      IF(IEWALD.GT.0) THEN
         CALL EWALDRECIPR(ESELF,ERECIPR,CORD,CHARG,FFGRD,
     *                    XX(LFFRKEXPEL),XX(LFFRKVEC),XX(LFFTCHCH),
     *                    XX(LFFCOSCH),XX(LFFSINCH),XX(LFFKVEC))
         ENCHAR =ENCHAR + ESELF + ERECIPR
      END IF
C
C
C     -- IF WT14CH IS NOT 1.0, ADD (WT14CH-1)
C        (HUI LI: NO EWALD CORRECTION TO THIS PART)
C
      IF(WT14CH.NE.1.0D+00) THEN
      FACT   =  WT14CH - 1.0D+00
      LLTODO = 1
      IF(IFEPTYP.GT.0.AND.IDOPOL.EQ.0) THEN
         IF(    NAT.LE.0) THEN       ! FFMD1
            IF(IFEPTYP.EQ.1)LLTODO= 4
            IF(IFEPTYP.EQ.2)LLTODO= 7
         END IF
      END IF
      DO 290 LL=1,LLTODO
      IF(IFEPTYP.EQ.2.AND.NAT.LE.0)THEN
         IF(LL.EQ.5) GOTO 290
      END IF
      IF(LL.EQ.1) THEN
         NN1  = L1N14J
         NN2  = L2N14J
         SIGN = 1.0D+00
         DWT1 = ZERO
         DWT2 = ZERO
      END IF
      IF(LL.EQ.2) THEN
         NN1  = L1N14A
         NN2  = L2N14A
         SIGN = -1.0D+00
         DWT1 = ZERO
         DWT2 = ZERO
      END IF
      IF(LL.EQ.3) THEN
         NN1  = L1N14A
         NN2  = L2N14A
         SIGN = 1 - WSIMUL
         DWT1 = 1 - WPERT1
         DWT2 = 1 - WPERT2
      END IF
      IF(LL.EQ.4) THEN
         NN1  = L1N14B
         NN2  = L2N14B
         SIGN = WSIMUL
         DWT1 = WPERT1
         DWT2 = WPERT2
      END IF
      IF(LL.EQ.5) THEN
         NN1  = L1N14B
         NN2  = L2N14B
         SIGN = -1.0D+00
         DWT1 = ZERO
         DWT2 = ZERO
      END IF
      IF(LL.EQ.6) THEN
         NN1  = L1N14PMA
         NN2  = L2N14PMA
         SIGN = -WSIMUL
         DWT1 = -WPERT1
         DWT2 = -WPERT2
      END IF
      IF(LL.EQ.7) THEN
         NN1  = L1N14PMB
         NN2  = L2N14PMB
         SIGN = +WSIMUL
         DWT1 = +WPERT1
         DWT2 = +WPERT2
      END IF
      DO 300 III=NN1, NN2
         IF(LL.EQ.1) THEN
            I = L14J(1,III)
            J = L14J(2,III)
         END IF
         IF(LL.EQ.2) THEN
            I = L14A(1,III)
            J = L14A(2,III)
         END IF
         IF(LL.EQ.3) THEN
            I = L14A(1,III)
            J = L14A(2,III)
         END IF
         IF(LL.EQ.4) THEN
            I = L14B(1,III)
            J = L14B(2,III)
         END IF
         IF(LL.EQ.5) THEN
            I = L14B(1,III)
            J = L14B(2,III)
         END IF
         IF(LL.EQ.6) THEN
            I = L14PMA(1,III)
            J = L14PMA(2,III)
         END IF
         IF(LL.EQ.7) THEN
            I = L14PMB(1,III)
            J = L14PMB(2,III)
         END IF
         IF(I.EQ.0.OR.J.EQ.0) GOTO 300
         QI  = CHARG(I)
         QJ  = CHARG(J)
         IF(IFEPTYP.GT.0.AND.(LL.EQ.4.OR.LL.EQ.7))THEN
            QI = CHARGB(I)
            QJ = CHARGB(J)
         END IF
         IF(QI.EQ.ZERO.OR.QJ.EQ.ZERO) GOTO 300
         QIQJ  = QI*QJ
C
         X     = CORD(1,I) - CORD(1,J)
         Y     = CORD(2,I) - CORD(2,J)
         Z     = CORD(3,I) - CORD(3,J)
         IF(IFEPTYP.GT.0.AND.(LL.EQ.4.OR.LL.EQ.7))THEN
            DO IFIXB = 1, NFIXMM
               IF(I.EQ.IFIXMM(IFIXB)) THEN
                  X = X + CORDB(1,I) - CORD(1,I)
                  Y = Y + CORDB(2,I) - CORD(2,I)
                  Z = Z + CORDB(3,I) - CORD(3,I)
               END IF
               IF(J.EQ.IFIXMM(IFIXB)) THEN
                  X = X - CORDB(1,J) + CORD(1,J)
                  Y = Y - CORDB(2,J) + CORD(2,J)
                  Z = Z - CORDB(3,J) + CORD(3,J)
               END IF
            ENDDO
         END IF
         PBCX  = XBOX * ANINT(X*ONEXBOX)
         PBCY  = YBOX * ANINT(Y*ONEYBOX)
         PBCZ  = ZBOX * ANINT(Z*ONEZBOX)
         X     = X - PBCX
         Y     = Y - PBCY
         Z     = Z - PBCZ
         R2    = X*X+Y*Y+Z*Z
         IF(R2.GT.SWRB2) GOTO 300
         IF(R2.LT.0.01D+00) GOTO 300
C
         R     = SQRT(R2)
         ONERX = ONE/R
         R     = R + RBUF
         R2    = R*R
         R3    = R*R2
         ONER  = ONE/R
         ONER2 = ONER*ONER
         ONER3 = ONER2*ONERX
C
            EDUM  = QIQJ*(ONER-TWOROFF+R*ONEROFF2
     *                       +R2*ONEROFF3-ONEROFF
     *                    -R*TWOROFF2+R3*ONEROFF4)*FACT*SIGN
            DUM   = -QIQJ*(ONER3-ONEROFF2*ONER-TWOROFF3
     *                      +TWOROFF2*ONER-THREEROFF4*R)*FACT*SIGN
            IF(LL.GT.2.AND.LL.NE.5) THEN
               EDUMW  = QIQJ*(ONER-TWOROFF+R*ONEROFF2
     *                           +R2*ONEROFF3-ONEROFF
     *                        -R*TWOROFF2+R3*ONEROFF4)*FACT
               EDUMW1 = EDUMW*DWT1
               EDUMW2 = EDUMW*DWT2
            END IF
         EPAIR = EDUM
         ENCHAR= ENCHAR  + EPAIR
         IF(LL.GE. 3.AND.LL.LE. 4) THEN
            SOL1CH= SOL1CH + (EDUMW1-EPAIR)
            SOL2CH= SOL2CH + (EDUMW2-EPAIR)
         END IF
         IF(LL.GE. 6.AND.LL.LE. 7) THEN
            PMF1CH= PMF1CH + (EDUMW1-EPAIR)
         END IF
         DEX   = DUM*X
         DEY   = DUM*Y
         DEZ   = DUM*Z
         FFGRD(1,I)=FFGRD(1,I) + DEX
         FFGRD(2,I)=FFGRD(2,I) + DEY
         FFGRD(3,I)=FFGRD(3,I) + DEZ
         FFGRD(1,J)=FFGRD(1,J) - DEX
         FFGRD(2,J)=FFGRD(2,J) - DEY
         FFGRD(3,J)=FFGRD(3,J) - DEZ
         VIR(1)    =VIR(1)     + DEX*X
         VIR(2)    =VIR(2)     + DEY*Y
         VIR(3)    =VIR(3)     + DEZ*Z
         IYES = 0
         JYES = 0
         DO KFIX=1,NFIXMM
            IF(I.EQ.IFIXMM(KFIX)) IYES = 1
            IF(J.EQ.IFIXMM(KFIX)) JYES = 1
         ENDDO
         IF(IYES.EQ.1.AND.JYES.EQ.1) THEN
            VIR(1)    =VIR(1)     - DEX*X
            VIR(2)    =VIR(2)     - DEY*Y
            VIR(3)    =VIR(3)     - DEZ*Z
         END IF
 300  CONTINUE
 290  CONTINUE
      END IF
C
C     -- EXTERNAL ELECTRIC FIELD --
      CALL EEFIELD(CORD,CHARG,FFGRD,XX(LFFZMAS))
C
      RETURN
      END
C*MODULE QUANPOB  *DECK EEFIELD
!>
!> @brief    calculate charge-field energy for MM atoms
!>
      SUBROUTINE EEFIELD(CORD,CHARG,FFGRD,ZMAS)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      LOGICAL GOPARR,DSKWRK,MASWRK
C
      PARAMETER (ZERO=0.0D+00)
C
      DIMENSION CORD(3,*),CHARG(*),FFGRD(3,*),ZMAS(*)
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
      COMMON /IOFILE/ IR,IW,IP,IJK,IJKT,IDAF,NAV,IODA(950)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
C
      INTEGER, PARAMETER :: K15 = SELECTED_INT_KIND(15)
      INTEGER, PARAMETER :: MAGIC2 = TRANSFER(20170303193144_K15,1)
C
C     HUI LI, FEB 26, 2017, LINCOLN
C
C     --  ENFIELD: MM CHARGE-ELECTRIC FIELD ENERGY
C
      ENFIELD =ZERO
      IF(ABS(EFIELDX)+ABS(EFIELDY)+ABS(EFIELDZ).EQ.ZERO) RETURN
C
      IF(NDFS.EQ.MAGIC2) THEN
C
      EFORCEX = EFIELDX*QDION/AMION
      EFORCEY = EFIELDY*QDION/AMION
      EFORCEZ = EFIELDZ*QDION/AMION
C
      IPCOUNT = ME - 1
      DO 100 IFFAT = 1, NATPDB
         IF(GOPARR) THEN
            IPCOUNT = IPCOUNT + 1
            IF(MOD(IPCOUNT,NPROC).NE.0) GOTO 100
         END IF
         FFGRD(1,IFFAT)=FFGRD(1,IFFAT) + EFORCEX*ZMAS(IFFAT)
         FFGRD(2,IFFAT)=FFGRD(2,IFFAT) + EFORCEY*ZMAS(IFFAT)
         FFGRD(3,IFFAT)=FFGRD(3,IFFAT) + EFORCEZ*ZMAS(IFFAT)
 100  CONTINUE
C
      ELSE
C
      IPCOUNT = ME - 1
      DO 200 IFFAT = 1, NFFAT
         IF(GOPARR) THEN
            IPCOUNT = IPCOUNT + 1
            IF(MOD(IPCOUNT,NPROC).NE.0) GOTO 200
         END IF
         IF(CHARG(IFFAT).EQ.ZERO) GOTO 200
         EFORCEX = EFIELDX*CHARG(IFFAT)
         EFORCEY = EFIELDY*CHARG(IFFAT)
         EFORCEZ = EFIELDZ*CHARG(IFFAT)
         FFGRD(1,IFFAT)=FFGRD(1,IFFAT) + EFORCEX
         FFGRD(2,IFFAT)=FFGRD(2,IFFAT) + EFORCEY
         FFGRD(3,IFFAT)=FFGRD(3,IFFAT) + EFORCEZ
         ENFIELD= ENFIELD + EFORCEX*CORD(1,IFFAT)
     *                    + EFORCEY*CORD(2,IFFAT)
     *                    + EFORCEZ*CORD(3,IFFAT)
C        - NOT SURE HOW TO DO IT RIGHT NOW -
C        VIR(1)        =VIR(1)         + DUM*XI*CORD(1,IFFAT)
C        VIR(2)        =VIR(2)         + DUM*YI*CORD(2,IFFAT)
C        VIR(3)        =VIR(3)         + DUM*ZI*CORD(3,IFFAT)
C
 200  CONTINUE
C
C     -- ADD THE ENERGY TO ENCHAR --
      ENCHAR = ENCHAR + ENFIELD
C
      END IF
C
      RETURN
      END
C*MODULE QUANPOB  *DECK CHGRXN
!>
!> @brief    reaction field
!>
!> @author   Dejun Si, Hui Li
!>           - Feb 2011
!>
!> @details  calculate the reaction field of MM charges
!>
      SUBROUTINE CHGRXN(CORD,FFGRD,CHARG,XTS,YTS,ZTS,
     *                  CMAT1,POT1,QRXN1,NTS)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      LOGICAL GOPARR,DSKWRK,MASWRK
C
      PARAMETER (ZERO=0.0D+00)
      PARAMETER (PT5=0.5D+00)
      PARAMETER (ONE=1.0D+00)
      PARAMETER (TWO=2.0D+00)
C
      DIMENSION CORD(3,*),FFGRD(3,*),CHARG(*),
     *          XTS(NTS),YTS(NTS),ZTS(NTS),POT1(NTS),QRXN1(NTS),
     *          CMAT1(NTS,NTS)
C
      COMMON /FFENGY/ EN12,EN123,EN123R4,EN123B4,EN234W1,ENCHAR,ENLJR,
     *                ENLJD,ENPOL,XENPOL,ENRXN,XENRXN,ENRXNPOL,ENRXNR,
     *                EN12312,ENQUANP(30),
     *                ENBIAS,ENCENT,ENUCCH,ENCMAP,ENPOT,ENKIN,ENTOT
      COMMON /FFMDPA/ DT,DT2,TEMP0,PRES0,POLTOL,VIR(3),PMEAN,VOLAV,
     *                ENPAV,ENKAV,TEMPAV,BERENDT,BERENDP,VELMAX,
     *                PMEANX,PMEANY,PMEANZ,NSTEP,KMASTER,KOUTACT(2),
     *                IHESS,INTALG,ITSTAT,IPSTAT,JOUT,KOUT,LOUT
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
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
C
C     DEJUN SI, HUI LI, FEB 2011, LINCOLN
C
C     --  CONTINUUM SOLVATION ENERGY --
C
      ENRXN =ZERO
      IF(IDOCHG.EQ.0) RETURN
      IF(RSPHSOL.GE.1.0D+30) RETURN
C
      R     = RSPHSOL
      R2    = R*R
      IF(ISPHSOL.EQ.1) THEN
         SCALE = (RXNEPS-ONE)/RXNEPS
         IPCOUNT = ME - 1
         DO 200 IFFAT = 1, NFFAT
            IF(CHARG(IFFAT).EQ.ZERO) GOTO 200
            QI    = CHARG(IFFAT)
            XI    = CORD(1,IFFAT) - CENTX
            YI    = CORD(2,IFFAT) - CENTY
            ZI    = CORD(3,IFFAT) - CENTZ
            RI2   = XI*XI+YI*YI+ZI*ZI
            IF(RI2.GT.(R-1.0D+00)**2) GOTO 200
            IF(GOPARR) THEN
               IPCOUNT = IPCOUNT + 1
               IF(MOD(IPCOUNT,NPROC).NE.0) GOTO 200
            END IF
            FAC   = R - RI2/R
            FAC2  = FAC*FAC
            RFAC2 = R*FAC2
            DUM   = -SCALE*QI*QI/RFAC2
            ENRXN = ENRXN - PT5*QI*QI/FAC
            FFGRD(1,IFFAT)=FFGRD(1,IFFAT) + DUM*XI
            FFGRD(2,IFFAT)=FFGRD(2,IFFAT) + DUM*YI
            FFGRD(3,IFFAT)=FFGRD(3,IFFAT) + DUM*ZI
            VIR(1)        =VIR(1)         + DUM*XI*CORD(1,IFFAT)
            VIR(2)        =VIR(2)         + DUM*YI*CORD(2,IFFAT)
            VIR(3)        =VIR(3)         + DUM*ZI*CORD(3,IFFAT)
            IYES = 0
            DO KFIX=1,NFIXMM
               IF(IFFAT.EQ.IFIXMM(KFIX)) IYES = 1
            ENDDO
            IF(IYES.EQ.1) THEN
               VIR(1)    =VIR(1)     - DUM*XI*CORD(1,IFFAT)
               VIR(2)    =VIR(2)     - DUM*YI*CORD(2,IFFAT)
               VIR(3)    =VIR(3)     - DUM*ZI*CORD(3,IFFAT)
            END IF
            DO 210 JFFAT = IFFAT+1, NFFAT
               IF(CHARG(JFFAT).EQ.ZERO) GOTO 210
               QJ    = CHARG(JFFAT)
               XJ    = CORD(1,JFFAT) - CENTX
               YJ    = CORD(2,JFFAT) - CENTY
               ZJ    = CORD(3,JFFAT) - CENTZ
               RJ2   = XJ*XJ+YJ*YJ+ZJ*ZJ
               IF(RJ2.GT.(R-1.0D+00)**2) GOTO 210
               RIRJ  = XI*XJ+YI*YJ+ZI*ZJ
               FAC   = RI2*RJ2/R2 + R2 - TWO*RIRJ
               FACRT = SQRT(FAC)
               DUM   = SCALE*QI*QJ/(FAC*FACRT)
               ENRXN = ENRXN - QI*QJ/FACRT
               FFGRD(1,IFFAT)=FFGRD(1,IFFAT) + DUM*(XI*RJ2/R2-XJ)
               FFGRD(2,IFFAT)=FFGRD(2,IFFAT) + DUM*(YI*RJ2/R2-YJ)
               FFGRD(3,IFFAT)=FFGRD(3,IFFAT) + DUM*(ZI*RJ2/R2-ZJ)
               FFGRD(1,JFFAT)=FFGRD(1,JFFAT) + DUM*(XJ*RI2/R2-XI)
               FFGRD(2,JFFAT)=FFGRD(2,JFFAT) + DUM*(YJ*RI2/R2-YI)
               FFGRD(3,JFFAT)=FFGRD(3,JFFAT) + DUM*(ZJ*RI2/R2-ZI)
               VIR(1) =VIR(1) + DUM*(XI*RJ2/R2-XJ)*CORD(1,IFFAT)
               VIR(2) =VIR(2) + DUM*(YI*RJ2/R2-YJ)*CORD(2,IFFAT)
               VIR(3) =VIR(3) + DUM*(ZI*RJ2/R2-ZJ)*CORD(3,IFFAT)
               VIR(1) =VIR(1) + DUM*(XJ*RI2/R2-XI)*CORD(1,JFFAT)
               VIR(2) =VIR(2) + DUM*(YJ*RI2/R2-YI)*CORD(2,JFFAT)
               VIR(3) =VIR(3) + DUM*(ZJ*RI2/R2-ZI)*CORD(3,JFFAT)
               IYES = 0
               JYES = 0
               DO KFIX=1,NFIXMM
                  IF(IFFAT.EQ.IFIXMM(KFIX)) IYES = 1
                  IF(JFFAT.EQ.IFIXMM(KFIX)) JYES = 1
               ENDDO
               IF(IYES.EQ.1.AND.JYES.EQ.1) THEN
                  VIR(1) =VIR(1) - DUM*(XI*RJ2/R2-XJ)*CORD(1,IFFAT)
                  VIR(2) =VIR(2) - DUM*(YI*RJ2/R2-YJ)*CORD(2,IFFAT)
                  VIR(3) =VIR(3) - DUM*(ZI*RJ2/R2-ZJ)*CORD(3,IFFAT)
                  VIR(1) =VIR(1) - DUM*(XJ*RI2/R2-XI)*CORD(1,JFFAT)
                  VIR(2) =VIR(2) - DUM*(YJ*RI2/R2-YI)*CORD(2,JFFAT)
                  VIR(3) =VIR(3) - DUM*(ZJ*RI2/R2-ZI)*CORD(3,JFFAT)
               END IF
 210        CONTINUE
 200     CONTINUE
         ENRXN = ENRXN*SCALE
      END IF
C
      IF(ISPHSOL.GE.60) THEN
         SCALE = (RXNEPS-ONE)/RXNEPS
         IPCOUNT = ME - 1
         DO 300 ITS=1, NTS
            POT1(ITS) = ZERO
            DO 310 JFFAT=1,NFFAT
               IF(GOPARR) THEN
                  IPCOUNT = IPCOUNT + 1
                  IF(MOD(IPCOUNT,NPROC).NE.0) GOTO 310
               END IF
               R2 =  (XTS(ITS)-CORD(1,JFFAT))**2
     *              +(YTS(ITS)-CORD(2,JFFAT))**2
     *              +(ZTS(ITS)-CORD(3,JFFAT))**2
               R  = SQRT(R2)
               POT1(ITS) = POT1(ITS)+CHARG(JFFAT)/R
 310        CONTINUE
 300     CONTINUE
         IF(GOPARR) CALL DDI_GSUMF(2405,POT1,NTS)
         IPCOUNT = ME - 1
         DO 320 ITS = 1, NTS
            QRXN1(ITS) = ZERO
            IF(GOPARR) THEN
               IPCOUNT = IPCOUNT + 1
               IF(MOD(IPCOUNT,NPROC).NE.0) GOTO 320
            END IF
            DO JTS = 1, NTS
               QRXN1(ITS) = QRXN1(ITS)-CMAT1(ITS,JTS)*POT1(JTS)
            ENDDO
            QRXN1(ITS) = QRXN1(ITS)*SCALE
            ENRXN      = ENRXN + QRXN1(ITS)*POT1(ITS)
 320     CONTINUE
         IF(GOPARR) CALL DDI_GSUMF(2406,QRXN1,NTS)
         ENRXN = PT5*ENRXN
      END IF
      IF(ISPHSOL.GE.60 .AND. IDOPOL.EQ.0) THEN
         IPCOUNT = ME - 1
         DO IFFAT=1,NFFAT
            QIFFAT = CHARG(IFFAT)
            XI     = CORD(1,IFFAT)
            YI     = CORD(2,IFFAT)
            ZI     = CORD(3,IFFAT)
            RI2    = (XI-CENTX)**2 + (YI-CENTY)**2 + (ZI-CENTZ)**2
            DO 340 JTS=1, NTS
               IF(GOPARR) THEN
                  IPCOUNT = IPCOUNT + 1
                  IF(MOD(IPCOUNT,NPROC).NE.0) GOTO 340
               END IF
               XJ   = XTS(JTS)
               YJ   = YTS(JTS)
               ZJ   = ZTS(JTS)
               QJTS = QRXN1(JTS)
               X    = XI - XJ
               Y    = YI - YJ
               Z    = ZI - ZJ
               R2   = X*X + Y*Y + Z*Z
               R    = SQRT(R2)
               R3   = R*R2
               DUM  = -QIFFAT*QJTS/R3
               FFGRD(1,IFFAT) = FFGRD(1,IFFAT) + DUM*X
               FFGRD(2,IFFAT) = FFGRD(2,IFFAT) + DUM*Y
               FFGRD(3,IFFAT) = FFGRD(3,IFFAT) + DUM*Z
               VIR(1)         = VIR(1)         + DUM*X*CORD(1,IFFAT)
               VIR(2)         = VIR(2)         + DUM*Y*CORD(2,IFFAT)
               VIR(3)         = VIR(3)         + DUM*Z*CORD(3,IFFAT)
               IYES = 0
               DO KFIX=1,NFIXMM
                  IF(IFFAT.EQ.IFIXMM(KFIX)) IYES = 1
               ENDDO
               IF(IYES.EQ.1) THEN
                  VIR(1)      = VIR(1)         - DUM*X*CORD(1,IFFAT)
                  VIR(2)      = VIR(2)         - DUM*Y*CORD(2,IFFAT)
                  VIR(3)      = VIR(3)         - DUM*Z*CORD(3,IFFAT)
               END IF
 340        CONTINUE
         ENDDO
      END IF
C
      RETURN
      END
C*MODULE QUANPOB  *DECK POLRXN
!>
!> @brief    reaction field and induced dipole
!>
!> @author   Hui Li, Dejun Si
!>           - Jan 2011
!>
!> @details  calculate the reaction field and induced dipole
!>           due to MM charges
!>
      SUBROUTINE POLRXN(CORD,FFGRD,CHARG,POL,POLSV,DIP,
     *                  FIELD1,FIELD2,FIELD3,
     *                  XTS,YTS,ZTS,CMAT1,POT1,POT2,QRXN1,QRXN2,NTS,
     *                  NONLS1,L1213J)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      LOGICAL GOPARR,DSKWRK,MASWRK
C
      PARAMETER (ONE=1.0D+00)
      PARAMETER (ZERO=0.0D+00)
      PARAMETER (PT5=0.5D+00)
      PARAMETER (TWO=2.0D+00)
      PARAMETER (THREE=3.0D+00)
      PARAMETER (FOUR=4.0D+00)
      PARAMETER (ONESIX=1.0D+00/6.0D+00)
C
      DIMENSION CORD(3,*),FFGRD(3,*),CHARG(*),
     *          POL(*),POLSV(*),DIP(3,*),FIELD1(3,*),FIELD2(3,*),
     *          FIELD3(3,*),
     *          XTS(NTS),YTS(NTS),ZTS(NTS),QRXN1(NTS),QRXN2(NTS),
     *          POT1(NTS),POT2(NTS),CMAT1(NTS,NTS),NONLS1(2,*),
     *          L1213J(2,*)
C
      COMMON /FFDAMP/ IPODAMP,IPO1213,SCRFAC,RSCRFAC,SCRF2,SCRF3,SCRF4
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
      COMMON /FFRMSD/ DIPT(8),TIMGYRA,TIMRALL,
     *                NATPDB,NGYRA,NDIEL,NRALL,
     *                LFFRALL0,NRMSD,LFFRMSD0,KFREEAB(201),
     *                NRIJMM,IJRMM(2,100),NRIJQM,IJRQM(2,100),
     *                NAIJKMM,IJKMM(3,100),NAIJKQM,IJKQM(3,100),
     *                NFIXMM,IFIXMM(200),NFIXQM,IFIXQM(200)
      COMMON /FFRXN / RXNEPS,RSPHSOL,ISPHSOL
      COMMON /FMCOM / XXX(1)
      COMMON /IOFILE/ IR,IW,IP,IJK,IJKT,IDAF,NAV,IODA(950)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
C
      INTEGER, PARAMETER :: K15 = SELECTED_INT_KIND(15)
      INTEGER, PARAMETER :: MAGIC2 = TRANSFER(20170303193144_K15,1)
C
C     HUI LI, DEJUN SI, JAN 2011, LINCOLN
C
      ENPOL    = ZERO
      ENRXNPOL = ZERO
      IF(IDOCHG.EQ.0) RETURN
      IF(IDOPOL.EQ.0) RETURN
C
      IF(NDFS.EQ.MAGIC2) THEN
         CALL POLCCS(CORD,FFGRD,CHARG,POL,DIP,FIELD1,
     *               XXX(LFFNONLSTQ))
         RETURN
      END IF
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
         QI    = CHARG(IFFAT)
         QJ    = CHARG(JFFAT)
         POLI  = POL(IFFAT)
         POLJ  = POL(JFFAT)
         IF(QI  .EQ.ZERO.AND.POLI.EQ.ZERO) GOTO 100
         IF(QJ  .EQ.ZERO.AND.POLJ.EQ.ZERO) GOTO 100
         IF(POLI.EQ.ZERO.AND.POLJ.EQ.ZERO) GOTO 100
         IF(QI  .EQ.ZERO.AND.  QJ.EQ.ZERO) GOTO 100
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
      IF(GOPARR) CALL DDI_GSUMF(2405,FIELD1,3*NFFAT)
C
C     -- FIELD2: FIELD AT POL DUE TO REACTION FIELD OF MM CHARGE --
C
      CALL VCLR(FIELD2,1,3*NFFAT)
C
      IF(ISPHSOL.EQ.1 .AND. RSPHSOL.LT.1.0D+30) THEN
      SCALE = (RXNEPS-ONE)/RXNEPS
      R     = RSPHSOL
      R2    = R*R
      R3    = R*R2
      R4    = R*R3
      IPCOUNT = ME - 1
      DO 200 IFFAT = 1, NFFAT
         IF(CHARG(IFFAT).EQ.ZERO.AND.POL(IFFAT).EQ.ZERO) GOTO 200
         QI    = CHARG(IFFAT)
         XI    = CORD(1,IFFAT) - CENTX
         YI    = CORD(2,IFFAT) - CENTY
         ZI    = CORD(3,IFFAT) - CENTZ
         RI2   = XI*XI+YI*YI+ZI*ZI
         IF(RI2.GT.(R-1.0D+00)**2) GOTO 200
         IF(GOPARR) THEN
            IPCOUNT = IPCOUNT + 1
            IF(MOD(IPCOUNT,NPROC).NE.0) GOTO 200
         END IF
         FAC   = R - RI2/R
         FAC2  = FAC*FAC
         RFAC2 = R*FAC2
         DUM   = SCALE*QI/RFAC2
         FIELD2(1,IFFAT) = FIELD2(1,IFFAT) + DUM*XI
         FIELD2(2,IFFAT) = FIELD2(2,IFFAT) + DUM*YI
         FIELD2(3,IFFAT) = FIELD2(3,IFFAT) + DUM*ZI
         DO 210 JFFAT = IFFAT+1, NFFAT
            IF(CHARG(JFFAT).EQ.ZERO.AND.POL(JFFAT).EQ.ZERO) GOTO 210
            QJ    = CHARG(JFFAT)
            XJ    = CORD(1,JFFAT) - CENTX
            YJ    = CORD(2,JFFAT) - CENTY
            ZJ    = CORD(3,JFFAT) - CENTZ
            RJ2   = XJ*XJ+YJ*YJ+ZJ*ZJ
            IF(RJ2.GT.(R-1.0D+00)**2) GOTO 210
            RIRJ  = XI*XJ+YI*YJ+ZI*ZJ
            FAC   = RI2*RJ2/R2 + R2 - TWO*RIRJ
            FACRT = SQRT(FAC)
            DUM   = -SCALE/(FAC*FACRT)
            DUMI  = DUM*QJ
            DUMJ  = DUM*QI
            FIELD2(1,IFFAT)=FIELD2(1,IFFAT) + DUMI*(XI*RJ2/R2-XJ)
            FIELD2(2,IFFAT)=FIELD2(2,IFFAT) + DUMI*(YI*RJ2/R2-YJ)
            FIELD2(3,IFFAT)=FIELD2(3,IFFAT) + DUMI*(ZI*RJ2/R2-ZJ)
            FIELD2(1,JFFAT)=FIELD2(1,JFFAT) + DUMJ*(XJ*RI2/R2-XI)
            FIELD2(2,JFFAT)=FIELD2(2,JFFAT) + DUMJ*(YJ*RI2/R2-YI)
            FIELD2(3,JFFAT)=FIELD2(3,JFFAT) + DUMJ*(ZJ*RI2/R2-ZI)
 210     CONTINUE
 200  CONTINUE
      IF(GOPARR) CALL DDI_GSUMF(2406,FIELD2,3*NFFAT)
      END IF
C
C     -- SOLVE INDUCED DIPOLES --
      CALL SOLVEDIP(CORD,POL,POLSV,DIP,ENPOL,
     *              FIELD1,FIELD2,FIELD3,
     *              NONLS1,L1213J,
     *              XTS,YTS,ZTS,CMAT1,POT2,QRXN1,QRXN2,NTS)
C
      ENRXNPOL = ZERO
      IF(ISPHSOL.EQ.1 .AND. RSPHSOL.LT.1.0D+30) THEN
         DO IFFAT=1,NFFAT
            ENRXNPOL = ENRXNPOL - FIELD2(1,IFFAT)*DIP(1,IFFAT)
     *                          - FIELD2(2,IFFAT)*DIP(2,IFFAT)
     *                          - FIELD2(3,IFFAT)*DIP(3,IFFAT)
         ENDDO
      ELSE IF(ISPHSOL.GE.60 .AND. RSPHSOL.LT.1.0D+30) THEN
         DO ITS = 1, NTS
            ENRXNPOL = ENRXNPOL + POT1(ITS)*QRXN2(ITS)
         ENDDO
      END IF
      ENRXNPOL = PT5*ENRXNPOL
C
C     -- COMPUTE GRADIENTS --
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
      DO 400 III=NN1, NN2
         IF(LL.EQ.1) THEN
            IFFAT = NONLS1(1,III)
            JFFAT = NONLS1(2,III)
         END IF
         IF(LL.EQ.2) THEN
            IFFAT = L1213J(1,III)
            JFFAT = L1213J(2,III)
         END IF
         IF(IFFAT.EQ.0.OR.JFFAT.EQ.0) GOTO 400
         IF(POL(IFFAT).EQ.ZERO.AND.CHARG(IFFAT).EQ.ZERO) GOTO 400
         QI    = CHARG(IFFAT)
         DIPIX = DIP(1,IFFAT)
         DIPIY = DIP(2,IFFAT)
         DIPIZ = DIP(3,IFFAT)
         IF(POL(JFFAT).EQ.ZERO.AND.CHARG(JFFAT).EQ.ZERO) GOTO 400
         QJ    = CHARG(JFFAT)
         DIPJX = DIP(1,JFFAT)
         DIPJY = DIP(2,JFFAT)
         DIPJZ = DIP(3,JFFAT)
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
         IF(R2.GT.SWRB2) GOTO 400
         IF(R2.LT.0.01D+00) GOTO 400
         R     = SQRT(R2)
C
         FE    = ONE
         FT    = ONE
         FEGRD = ZERO
         FTGRD = ZERO
         IF(IPODAMP.GT.0) THEN
            POLIJ = POLSV(IFFAT)*POLSV(JFFAT) + 1.0D-60
            POL6  = POLIJ**ONESIX
            RPOL6 = ONE/POL6
            RSFAC = RSCRFAC*RPOL6
            VFAC  = R*RSFAC
            VFAC2 = VFAC*VFAC
            VFAC3 = VFAC2*VFAC
            VFAC4 = VFAC3*VFAC
            UFAC  = R*RPOL6
            UFAC2 = UFAC*UFAC
            UFAC3 = UFAC2*UFAC
            IF(IPODAMP.EQ.1.AND.VFAC.LT.ONE) THEN
C              -- LINEAR THOLE MODEL
               FE    = FOUR*VFAC3-THREE*VFAC4
               FT    = VFAC4
               FEGRD = 12.0D+00*(VFAC2-VFAC3)*RSFAC
               FTGRD = FOUR*VFAC3*RSFAC
            ELSE IF(IPODAMP.EQ.2) THEN
C              -- EXPONENTIAL THOLE MODEL
               VALEXP= EXP(-SCRFAC*UFAC)
               FE    = ONE-(PT5*SCRF2*UFAC2+SCRFAC*UFAC+ONE)*VALEXP
               FT    = ONE-(ONESIX*SCRF3*UFAC3+PT5*SCRF2*UFAC2+
     *                      SCRFAC*UFAC+ONE)*VALEXP
               FEGRD = PT5*SCRF3*UFAC2*RPOL6*VALEXP
               FTGRD = ONESIX*SCRF4*UFAC3*RPOL6*VALEXP
            ELSE IF(IPODAMP.EQ.3) THEN
C              -- THOLE-TINKER EXPONENTIAL MODEL
               VALEXP= EXP(-SCRFAC*UFAC3)
               FE    = ONE-VALEXP
               FT    = ONE-(ONE+SCRFAC*UFAC3)*VALEXP
               FEGRD = THREE*SCRFAC*UFAC2*RPOL6*VALEXP
               FTGRD = FEGRD*SCRFAC*UFAC3
            END IF
         END IF
C
         ONER  = ONE/R
         IF(IPOLSHF.EQ.1) CALL SHIFT(R2,R,ONER,X,Y,Z)
         IF(IPOLSHF.EQ.0) CALL SWFUNC(R2,X,Y,Z)
C
         XX    = X*X
         YY    = Y*Y
         ZZ    = Z*Z
         XY    = X*Y
         XZ    = X*Z
         YZ    = Y*Z
         ONER2 = ONER*ONER
         ONER3 = ONER2*ONER
         ONER4 = ONER2*ONER2
         ONER5 = ONER2*ONER3
         ONER6 = ONER3*ONER3
         ONER7 = ONER2*ONER5
C
C        - FORCES BETWEEN CHARGE AND DIPOLE
C
         QIDOT = THREE*ONER5*(DIPJX*X+DIPJY*Y+DIPJZ*Z)*QI
         QJDOT = THREE*ONER5*(DIPIX*X+DIPIY*Y+DIPIZ*Z)*QJ
         QIONER3= QI*ONER3
         QJONER3= QJ*ONER3
C              NEGATIVE FORCE ON QI BY DIPJ
         DXI   = DIPJX*QIONER3 - QIDOT*X
         DYI   = DIPJY*QIONER3 - QIDOT*Y
         DZI   = DIPJZ*QIONER3 - QIDOT*Z
C              NEGATIVE FORCE ON QJ BY DIPI
C                       FORCE ON DIPI BY QJ
         DXJ   = DIPIX*QJONER3 - QJDOT*X
         DYJ   = DIPIY*QJONER3 - QJDOT*Y
         DZJ   = DIPIZ*QJONER3 - QJDOT*Z
C              NEGATIVE FORCE ON IFFAT
         DX    = DXI - DXJ
         DY    = DYI - DYJ
         DZ    = DZI - DZJ
C
C        - FORCES BETWEEN DIPOLE AND DIPOLE -
C
         IF(IPO1213.EQ.1.AND.LL.EQ.2) GOTO 330
         IF(IDOPOL.GT.1) THEN
         DOTM  = DIPJX*DIPIX + DIPJY*DIPIY + DIPJZ*DIPIZ
         DOTXY = DIPJY*DIPIX + DIPIY*DIPJX
         DOTXZ = DIPJZ*DIPIX + DIPIZ*DIPJX
         DOTYZ = DIPJY*DIPIZ + DIPIY*DIPJZ
         DOTX  = TWO*DIPJX*DIPIX*X + DOTXY*Y + DOTXZ*Z
         DOTY  = TWO*DIPJY*DIPIY*Y + DOTXY*X + DOTYZ*Z
         DOTZ  = TWO*DIPJZ*DIPIZ*Z + DOTXZ*X + DOTYZ*Y
         DUM   = THREE*ONER5
         TEMP  = DIPIX*DIPJX*XX
     *          +DIPIY*DIPJY*YY
     *          +DIPIZ*DIPJZ*ZZ
     *          +(DIPIX*DIPJY+DIPIY*DIPJX)*XY
     *          +(DIPIX*DIPJZ+DIPIZ*DIPJX)*XZ
     *          +(DIPIY*DIPJZ+DIPIZ*DIPJY)*YZ
         DUM7  = 15.0D+00*TEMP*ONER7*FT
         FEDOTM= FE*DOTM
         FTDOTX= FT*DOTX
         FTDOTY= FT*DOTY
         FTDOTZ= FT*DOTZ
         DUM4  = FEGRD*DOTM*ONER4
         DUM6  = THREE*FTGRD*TEMP*ONER6
C        - NEGATIVE FORCE ON IFFAT
         DX  = DX - DUM*(FEDOTM*X+FTDOTX)+(DUM4-DUM6+DUM7)*X
         DY  = DY - DUM*(FEDOTM*Y+FTDOTY)+(DUM4-DUM6+DUM7)*Y
         DZ  = DZ - DUM*(FEDOTM*Z+FTDOTZ)+(DUM4-DUM6+DUM7)*Z
         END IF
C
 330     CONTINUE
C
C        - ENERGY-SWFDX TERMS
C
         FLDIX = QJONER3*X      ! FIELD AT DIPI DUE TO QJ
         FLDIY = QJONER3*Y
         FLDIZ = QJONER3*Z
         FLDJX = QIONER3*(-X)   ! FIELD AT DIPJ DUE TO QI
         FLDJY = QIONER3*(-Y)
         FLDJZ = QIONER3*(-Z)
         EDIQJ = -(FLDIX*DIPIX+FLDIY*DIPIY+FLDIZ*DIPIZ) ! E = -F*D
         EDJQI = -(FLDJX*DIPJX+FLDJY*DIPJY+FLDJZ*DIPJZ)
C
         IF(IDOPOL.EQ.1) THEN
            EDIDJ = ZERO
         ELSE
            DOTJ  = THREE*ONER5*FT*(DIPJX*X+DIPJY*Y+DIPJZ*Z)
            FLDIX = -DIPJX*ONER3*FE + DOTJ*X  ! FIELD AT DIPI DUE TO DIPJ
            FLDIY = -DIPJY*ONER3*FE + DOTJ*Y
            FLDIZ = -DIPJZ*ONER3*FE + DOTJ*Z
            EDIDJ = -(DIPIX*FLDIX + DIPIY*FLDIY + DIPIZ*FLDIZ) ! E=-F*D
         END IF
         IF(IPO1213.EQ.1.AND.LL.EQ.2) EDIDJ=ZERO
         EPAIR = (EDIQJ + EDJQI + EDIDJ)*SIGN
C
         DEX   = DX*SWF*SIGN
         DEY   = DY*SWF*SIGN
         DEZ   = DZ*SWF*SIGN
         FFGRD(1,IFFAT)=FFGRD(1,IFFAT) + DEX
         FFGRD(2,IFFAT)=FFGRD(2,IFFAT) + DEY
         FFGRD(3,IFFAT)=FFGRD(3,IFFAT) + DEZ
         FFGRD(1,IFFAT)=FFGRD(1,IFFAT)+EPAIR*SWFDX
         FFGRD(2,IFFAT)=FFGRD(2,IFFAT)+EPAIR*SWFDY
         FFGRD(3,IFFAT)=FFGRD(3,IFFAT)+EPAIR*SWFDZ
         FFGRD(1,JFFAT)=FFGRD(1,JFFAT) - DEX
         FFGRD(2,JFFAT)=FFGRD(2,JFFAT) - DEY
         FFGRD(3,JFFAT)=FFGRD(3,JFFAT) - DEZ
         FFGRD(1,JFFAT)=FFGRD(1,JFFAT)-EPAIR*SWFDX
         FFGRD(2,JFFAT)=FFGRD(2,JFFAT)-EPAIR*SWFDY
         FFGRD(3,JFFAT)=FFGRD(3,JFFAT)-EPAIR*SWFDZ
         VIR(1)    =VIR(1)     + DEX*X + EPAIR*SWFDX*X
         VIR(2)    =VIR(2)     + DEY*Y + EPAIR*SWFDY*Y
         VIR(3)    =VIR(3)     + DEZ*Z + EPAIR*SWFDZ*Z
         IYES = 0
         JYES = 0
         DO KFIX=1,NFIXMM
            IF(IFFAT.EQ.IFIXMM(KFIX)) IYES = 1
            IF(JFFAT.EQ.IFIXMM(KFIX)) JYES = 1
         ENDDO
         IF(IYES.EQ.1.AND.JYES.EQ.1) THEN
            VIR(1)    =VIR(1)     - DEX*X - EPAIR*SWFDX*X
            VIR(2)    =VIR(2)     - DEY*Y - EPAIR*SWFDY*Y
            VIR(3)    =VIR(3)     - DEZ*Z - EPAIR*SWFDZ*Z
         END IF
 400  CONTINUE
      ENDDO
C
C     -- FORCE CORRECTION DUE TO REACTION FIELD DIRECT METHOD --
C        (IMAGE CHARGE METHOD)
C
      IF(ISPHSOL.EQ.1 .AND. RSPHSOL.LT.1.0D+30) THEN
      R    = RSPHSOL
      R2   = R*R
      R3   = R*R2
      R4   = R2*R2
      IPCOUNT = ME - 1
      DO 500 IFFAT=1,NFFAT
         IF(POL(IFFAT).EQ.ZERO.AND.CHARG(IFFAT).EQ.ZERO) GOTO 500
         QI     = CHARG(IFFAT)
         XI     = CORD(1,IFFAT) - CENTX
         YI     = CORD(2,IFFAT) - CENTY
         ZI     = CORD(3,IFFAT) - CENTZ
         RI2    = XI*XI+YI*YI+ZI*ZI
         IF(RI2.GT.(R-1.0D+00)**2) GOTO 500
         IF(GOPARR) THEN
            IPCOUNT = IPCOUNT + 1
            IF(MOD(IPCOUNT,NPROC).NE.0) GOTO 500
         END IF
         DIX    = DIP(1,IFFAT)
         DIY    = DIP(2,IFFAT)
         DIZ    = DIP(3,IFFAT)
C
         DEDX1 = ZERO
         DEDY1 = ZERO
         DEDZ1 = ZERO
         DEDX2 = ZERO
         DEDY2 = ZERO
         DEDZ2 = ZERO
         DEDX3 = ZERO
         DEDY3 = ZERO
         DEDZ3 = ZERO
         DO 510 JFFAT=1,NFFAT
            QJ     = CHARG(JFFAT)
            XJ     = CORD(1,JFFAT) - CENTX
            YJ     = CORD(2,JFFAT) - CENTY
            ZJ     = CORD(3,JFFAT) - CENTZ
            RJ2    = XJ*XJ+YJ*YJ+ZJ*ZJ
            IF(RJ2.GT.(R-1.0D+00)**2) GOTO 510
            DJX    = DIP(1,JFFAT)
            DJY    = DIP(2,JFFAT)
            DJZ    = DIP(3,JFFAT)
            DJRJ   = DJX*XJ + DJY*YJ + DJZ*ZJ
            DJRI   = DJX*XI + DJY*YI + DJZ*ZI
            RIRJ   = XI*XJ+YI*YJ+ZI*ZJ
            FAC2   = ONE/(RI2*RJ2 + R4 -TWO*R2*RIRJ)
            FAC    = SQRT(FAC2)
            FAC3   = FAC*FAC2
            FAC5   = FAC3*FAC2
            FAC7   = FAC5*FAC2
C
C         (A) FORCE ON CHARGE(IFFAT) BY REACTION FIELD OF DIPOLE(JFFAT)
C
            FDIMGX = - (TWO*R*DJRJ*XI - R3*DJX)*FAC3
     *               + (R*RI2*DJRJ - R3*DJRI)*FAC5
     *                      *THREE*(RJ2*XI - R2*XJ)
            FDIMGY = - (TWO*R*DJRJ*YI - R3*DJY)*FAC3
     *               + (R*RI2*DJRJ - R3*DJRI)*FAC5
     *                      *THREE*(RJ2*YI - R2*YJ)
            FDIMGZ = - (TWO*R*DJRJ*ZI - R3*DJZ)*FAC3
     *               + (R*RI2*DJRJ - R3*DJRI)*FAC5
     *                      *THREE*(RJ2*ZI - R2*ZJ)
            DEDX1 = DEDX1 - SCALE*QI*FDIMGX
            DEDY1 = DEDY1 - SCALE*QI*FDIMGY
            DEDZ1 = DEDZ1 - SCALE*QI*FDIMGZ
C
C         (B) FORCE ON DIPOLE(IFFAT) BY REACTION FIELD OF CHARGE(JFFAT)
C
            DUM2   = FAC2*R2
            DUM    = FAC*R
            DUM3   = DUM2*DUM
            DUM5   = DUM2*DUM3
            FDXDX  = DUM5*QJ*THREE*(RJ2/R2*XI-XJ)*(RJ2/R2*XI-XJ)
     *              -DUM3*QJ*(RJ2/R2)
            FDXDY  = DUM5*QJ*THREE*(RJ2/R2*XI-XJ)*(RJ2/R2*YI-YJ)
            FDXDZ  = DUM5*QJ*THREE*(RJ2/R2*XI-XJ)*(RJ2/R2*ZI-ZJ)
            FDYDX  = DUM5*QJ*THREE*(RJ2/R2*YI-YJ)*(RJ2/R2*XI-XJ)
            FDYDY  = DUM5*QJ*THREE*(RJ2/R2*YI-YJ)*(RJ2/R2*YI-YJ)
     *              -DUM3*QJ*(RJ2/R2)
            FDYDZ  = DUM5*QJ*THREE*(RJ2/R2*YI-YJ)*(RJ2/R2*ZI-ZJ)
            FDZDX  = DUM5*QJ*THREE*(RJ2/R2*ZI-ZJ)*(RJ2/R2*XI-XJ)
            FDZDY  = DUM5*QJ*THREE*(RJ2/R2*ZI-ZJ)*(RJ2/R2*YI-YJ)
            FDZDZ  = DUM5*QJ*THREE*(RJ2/R2*ZI-ZJ)*(RJ2/R2*ZI-ZJ)
     *              -DUM3*QJ*(RJ2/R2)
C
            DEDX2 = DEDX2 - SCALE*(DIX*FDXDX+DIY*FDYDX+DIZ*FDZDX)
            DEDY2 = DEDY2 - SCALE*(DIX*FDXDY+DIY*FDYDY+DIZ*FDZDY)
            DEDZ2 = DEDZ2 - SCALE*(DIX*FDXDZ+DIY*FDYDZ+DIZ*FDZDZ)
C
C         (C) FORCE ON DIPOLE(IFFAT) BY REACTION FIELD OF DIPOLE(JFFAT)
C
            FDXDX  =-FAC3*TWO*R*DJRJ
     *              +FAC5*(TWO*R*DJRJ*XI-R3*DJX)*THREE*(RJ2*XI-R2*XJ)
     *              +FAC5*(TWO*R*DJRJ*XI-R3*DJX)*THREE*(RJ2*XI-R2*XJ)
     *              +FAC5*(R*RI2*DJRJ-R3*DJRI)*THREE*RJ2
     *              -FAC7*(R*RI2*DJRJ-R3*DJRI)*15.0D+00
     *                   *(RJ2*XI-R2*XJ)*(RJ2*XI-R2*XJ)
            FDXDY  =+FAC5*(TWO*R*DJRJ*XI-R3*DJX)*THREE*(RJ2*YI-R2*YJ)
     *              +FAC5*(TWO*R*DJRJ*YI-R3*DJY)*THREE*(RJ2*XI-R2*XJ)
     *              -FAC7*(R*RI2*DJRJ-R3*DJRI)*15.0D+00
     *                   *(RJ2*XI-R2*XJ)*(RJ2*YI-R2*YJ)
            FDXDZ  =+FAC5*(TWO*R*DJRJ*XI-R3*DJX)*THREE*(RJ2*ZI-R2*ZJ)
     *              +FAC5*(TWO*R*DJRJ*ZI-R3*DJZ)*THREE*(RJ2*XI-R2*XJ)
     *              -FAC7*(R*RI2*DJRJ-R3*DJRI)*15.0D+00
     *                   *(RJ2*XI-R2*XJ)*(RJ2*ZI-R2*ZJ)
            FDYDX  =+FAC5*(TWO*R*DJRJ*YI-R3*DJY)*THREE*(RJ2*XI-R2*XJ)
     *              +FAC5*(TWO*R*DJRJ*XI-R3*DJX)*THREE*(RJ2*YI-R2*YJ)
     *              -FAC7*(R*RI2*DJRJ-R3*DJRI)*15.0D+00
     *                   *(RJ2*YI-R2*YJ)*(RJ2*XI-R2*XJ)
            FDYDY  =-FAC3*TWO*R*DJRJ
     *              +FAC5*(TWO*R*DJRJ*YI-R3*DJY)*THREE*(RJ2*YI-R2*YJ)
     *              +FAC5*(TWO*R*DJRJ*YI-R3*DJY)*THREE*(RJ2*YI-R2*YJ)
     *              +FAC5*(R*RI2*DJRJ-R3*DJRI)*THREE*RJ2
     *              -FAC7*(R*RI2*DJRJ-R3*DJRI)*15.0D+00
     *                   *(RJ2*YI-R2*YJ)*(RJ2*YI-R2*YJ)
            FDYDZ  =+FAC5*(TWO*R*DJRJ*YI-R3*DJY)*THREE*(RJ2*ZI-R2*ZJ)
     *              +FAC5*(TWO*R*DJRJ*ZI-R3*DJZ)*THREE*(RJ2*YI-R2*YJ)
     *              -FAC7*(R*RI2*DJRJ-R3*DJRI)*15.0D+00
     *                   *(RJ2*YI-R2*YJ)*(RJ2*ZI-R2*ZJ)
            FDZDX  =+FAC5*(TWO*R*DJRJ*ZI-R3*DJZ)*THREE*(RJ2*XI-R2*XJ)
     *              +FAC5*(TWO*R*DJRJ*XI-R3*DJX)*THREE*(RJ2*ZI-R2*ZJ)
     *              -FAC7*(R*RI2*DJRJ-R3*DJRI)*15.0D+00
     *                   *(RJ2*ZI-R2*ZJ)*(RJ2*XI-R2*XJ)
            FDZDY  =+FAC5*(TWO*R*DJRJ*ZI-R3*DJZ)*THREE*(RJ2*YI-R2*YJ)
     *              +FAC5*(TWO*R*DJRJ*YI-R3*DJY)*THREE*(RJ2*ZI-R2*ZJ)
     *              -FAC7*(R*RI2*DJRJ-R3*DJRI)*15.0D+00
     *                   *(RJ2*ZI-R2*ZJ)*(RJ2*YI-R2*YJ)
            FDZDZ  =-FAC3*TWO*R*DJRJ
     *              +FAC5*(TWO*R*DJRJ*ZI-R3*DJZ)*THREE*(RJ2*ZI-R2*ZJ)
     *              +FAC5*(TWO*R*DJRJ*ZI-R3*DJZ)*THREE*(RJ2*ZI-R2*ZJ)
     *              +FAC5*(R*RI2*DJRJ-R3*DJRI)*THREE*RJ2
     *              -FAC7*(R*RI2*DJRJ-R3*DJRI)*15.0D+00
     *                   *(RJ2*ZI-R2*ZJ)*(RJ2*ZI-R2*ZJ)
C
            DEDX3 = DEDX3 - SCALE*(DIX*FDXDX+DIY*FDYDX+DIZ*FDZDX)
            DEDY3 = DEDY3 - SCALE*(DIX*FDXDY+DIY*FDYDY+DIZ*FDZDY)
            DEDZ3 = DEDZ3 - SCALE*(DIX*FDXDZ+DIY*FDYDZ+DIZ*FDZDZ)
C
 510     CONTINUE
         FFGRD(1,IFFAT)=FFGRD(1,IFFAT) + DEDX1 + DEDX2 + DEDX3
         FFGRD(2,IFFAT)=FFGRD(2,IFFAT) + DEDY1 + DEDY2 + DEDY3
         FFGRD(3,IFFAT)=FFGRD(3,IFFAT) + DEDZ1 + DEDZ2 + DEDZ3
         VIR(1) =VIR(1) + (DEDX1 + DEDX2 + DEDX3)*CORD(1,IFFAT)
         VIR(2) =VIR(2) + (DEDY1 + DEDY2 + DEDY3)*CORD(2,IFFAT)
         VIR(3) =VIR(3) + (DEDZ1 + DEDZ2 + DEDZ3)*CORD(3,IFFAT)
         IYES = 0
         DO KFIX=1,NFIXMM
            IF(IFFAT.EQ.IFIXMM(KFIX)) IYES = 1
         ENDDO
         IF(IYES.EQ.1) THEN
            VIR(1) =VIR(1) - (DEDX1 + DEDX2 + DEDX3)*CORD(1,IFFAT)
            VIR(2) =VIR(2) - (DEDY1 + DEDY2 + DEDY3)*CORD(2,IFFAT)
            VIR(3) =VIR(3) - (DEDZ1 + DEDZ2 + DEDZ3)*CORD(3,IFFAT)
         END IF
 500  CONTINUE
      END IF
C
C
C     -- FORCE DUE TO REACTION FIELD SURFACE CHARGE --
C        (CPCM OR COSMO STYLE)
C
      IF(ISPHSOL.GE.60 .AND. RSPHSOL.LT.1.0D+30) THEN
         IPCOUNT = ME - 1
         DO 600 IFFAT=1,NFFAT
            QIFFAT= CHARG(IFFAT)
            DIX   = DIP(1,IFFAT)
            DIY   = DIP(2,IFFAT)
            DIZ   = DIP(3,IFFAT)
            XI    = CORD(1,IFFAT)
            YI    = CORD(2,IFFAT)
            ZI    = CORD(3,IFFAT)
            RI2   = (XI-CENTX)**2 + (YI-CENTY)**2 + (ZI-CENTZ)**2
            DO 610 JTS = 1, NTS
               IF(GOPARR) THEN
                  IPCOUNT = IPCOUNT + 1
                  IF(MOD(IPCOUNT,NPROC).NE.0) GOTO 610
               END IF
               QJ    = QRXN1(JTS) + QRXN2(JTS)
               XJ    = XTS(JTS)
               YJ    = YTS(JTS)
               ZJ    = ZTS(JTS)
               X     = XI - XJ
               Y     = YI - YJ
               Z     = ZI - ZJ
               R2    = X*X + Y*Y + Z*Z
               ONER2 = ONE/R2
               ONER  = SQRT(ONER2)
               ONER3 = ONER2*ONER
               ONER5 = ONER3*ONER2
               DUM   = -QIFFAT*QJ*ONER3
               DXA   = DUM*X
               DYA   = DUM*Y
               DZA   = DUM*Z
               QJONER3= QJ*ONER3
               QJDOT = THREE*ONER5*(DIX*X+DIY*Y+DIZ*Z)*QJ
               DXB   = -DIX*QJONER3 + QJDOT*X
               DYB   = -DIY*QJONER3 + QJDOT*Y
               DZB   = -DIZ*QJONER3 + QJDOT*Z
               FFGRD(1,IFFAT)=FFGRD(1,IFFAT) + DXA + DXB
               FFGRD(2,IFFAT)=FFGRD(2,IFFAT) + DYA + DYB
               FFGRD(3,IFFAT)=FFGRD(3,IFFAT) + DZA + DZB
               VIR(1) =VIR(1) + (DXA + DXB)*CORD(1,IFFAT)
               VIR(2) =VIR(2) + (DYA + DYB)*CORD(2,IFFAT)
               VIR(3) =VIR(3) + (DZA + DZB)*CORD(3,IFFAT)
               IYES = 0
               DO KFIX=1,NFIXMM
                  IF(IFFAT.EQ.IFIXMM(KFIX)) IYES = 1
               ENDDO
               IF(IYES.EQ.1) THEN
                  VIR(1) =VIR(1) - (DXA + DXB)*CORD(1,IFFAT)
                  VIR(2) =VIR(2) - (DYA + DYB)*CORD(2,IFFAT)
                  VIR(3) =VIR(3) - (DZA + DZB)*CORD(3,IFFAT)
               END IF
 610        CONTINUE
 600     CONTINUE
      END IF
C
      RETURN
      END
C*MODULE QUANPOB  *DECK SOLVEDIP
!>
!> @brief    solve for the induced dipole
!>
!> @author   Hui Li, Dejun Si
!>           - Jan 2011
!>
!> @details  iteratively solve for the induced dipoles
!>
      SUBROUTINE SOLVEDIP(CORD,POL,POLSV,DIP,ENPOL,FIELD1,FIELD2,
     *                    FIELD3,NONLS1,L1213J,
     *                    XTS,YTS,ZTS,CMAT1,POT2,QRXN1,QRXN2,NTS)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      LOGICAL GOPARR,DSKWRK,MASWRK
C
      PARAMETER (ONE=1.0D+00)
      PARAMETER (ZERO=0.0D+00)
      PARAMETER (PT5=0.5D+00)
      PARAMETER (TWO=2.0D+00)
      PARAMETER (THREE=3.0D+00)
      PARAMETER (FOUR=4.0D+00)
      PARAMETER (ONESIX=1.0D+00/6.0D+00)
C
      DIMENSION CORD(3,*),POL(*),DIP(3,*),FIELD1(3,*),FIELD2(3,*),
     *          FIELD3(3,*),NONLS1(2,*),L1213J(2,*),POLSV(*),
     *          XTS(NTS),YTS(NTS),ZTS(NTS),POT2(NTS),
     *          QRXN1(NTS),QRXN2(NTS),CMAT1(NTS,NTS)
C
      COMMON /FMCOM / XX(1)
      COMMON /FFDAMP/ IPODAMP,IPO1213,SCRFAC,RSCRFAC,SCRF2,SCRF3,SCRF4
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
      COMMON /FFMPT3/ NACTMM,LACTMM(2020),LFFDIPOLD,JUMBUP,
     *                NACTQM,LACTQM(2020),LFFOLDC,LFFQMVELSV,MMHESS,
     *                LFFQMCHG,LFFQMCHGB,ISWAP,R2SWAP,DFTBMM
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
      COMMON /IOFILE/ IR,IW,IP,IJK,IJKT,IDAF,NAV,IODA(950)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
C
C     HUI LI, DEJUN SI, JAN 2011, LINCOLN
C     HUI LI, MAR 2, 2014, ADD DIPOLD
C
C     -- ITERATIVE METHODS ARE EQUIVALENT TO MATRIX INVERSION
C
C
      IF(IDOPOL.EQ.0) THEN
         ENPOL = ZERO
         CALL VCLR(DIP,1,3*NFFAT)
         CALL VCLR(QRXN2,1,NTS)
         RETURN
      END IF
C
      SUMOLD=ZERO
      DO 305 ITER = 1, IDOPOL
      SUM = ZERO
      CALL VCLR(FIELD3,1,3*NFFAT)
      IF(IDOPOL.EQ.1) GOTO 301
C
C     -- FIELD3: FIELD AT POL DUE TO INDUCED DIPOLES --
C
      NTIMES=2
      IF(IPO1213.EQ.1) NTIMES=1
      DO LL=1,NTIMES
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
      DO 300 III=NN1, NN2
         IF(LL.EQ.1) THEN
            IFFAT = NONLS1(1,III)
            JFFAT = NONLS1(2,III)
         END IF
         IF(LL.EQ.2) THEN
            IFFAT = L1213J(1,III)
            JFFAT = L1213J(2,III)
         END IF
         IF(IFFAT.EQ.0.OR.JFFAT.EQ.0) GOTO 300
         IF(POL(IFFAT).EQ.ZERO) GOTO 300
         DIPIX = DIP(1,IFFAT)
         DIPIY = DIP(2,IFFAT)
         DIPIZ = DIP(3,IFFAT)
         IF(POL(JFFAT).EQ.ZERO) GOTO 300
         DIPJX = DIP(1,JFFAT)
         DIPJY = DIP(2,JFFAT)
         DIPJZ = DIP(3,JFFAT)
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
         IF(R2.GT.SWRB2) GOTO 300
         IF(R2.LT.0.01D+00) GOTO 300
         R     = SQRT(R2)
C
         FE    = ONE
         FT    = ONE
         IF(IPODAMP.GT.0) THEN
            POL6  = (POLSV(IFFAT)*POLSV(JFFAT))**ONESIX
            RPOL6 = ONE/POL6
            RSFAC = RSCRFAC*RPOL6
            VFAC  = R*RSFAC
            VFAC3 = VFAC*VFAC*VFAC
            VFAC4 = VFAC3*VFAC
            UFAC  = R*RPOL6
            UFAC2 = UFAC*UFAC
            UFAC3 = UFAC2*UFAC
            IF(IPODAMP.EQ.1.AND.VFAC.LT.ONE) THEN
C              -- LINEAR THOLE MODEL
               FE    = FOUR*VFAC3-THREE*VFAC4
               FT    = VFAC4
            ELSE IF(IPODAMP.EQ.2) THEN
C              -- EXPONENTIAL THOLE MODEL
               VALEXP= EXP(-SCRFAC*UFAC)
               FE    = ONE-(PT5*SCRF2*UFAC2+SCRFAC*UFAC+ONE)*VALEXP
               FT    = ONE-(ONESIX*SCRF3*UFAC3+PT5*SCRF2*UFAC2+
     *                      SCRFAC*UFAC+ONE)*VALEXP
            ELSE IF(IPODAMP.EQ.3) THEN
C              -- THOLE-TINKER EXPONENTIAL MODEL
               VALEXP= EXP(-SCRFAC*UFAC3)
               FE    = ONE-VALEXP
               FT    = ONE-(ONE+SCRFAC*UFAC3)*VALEXP
            END IF
         END IF
C
         ONER  = ONE/R
         IF(IPOLSHF.EQ.1) CALL SHIFT(R2,R,ONER,X,Y,Z)
         IF(IPOLSHF.EQ.0) CALL SWFUNC(R2,X,Y,Z)
C
         ONER2 = ONER*ONER
         ONER3 = ONER2*ONER
         ONER5 = ONER2*ONER3
         FER3  = ONER3*FE
         FTR5  = ONER5*FT
         DOTJ  = THREE*FTR5*(DIPJX*X+DIPJY*Y+DIPJZ*Z)
         DOTI  = THREE*FTR5*(DIPIX*X+DIPIY*Y+DIPIZ*Z)
         SWF   = SWF*SIGN
         FIELD3(1,IFFAT)=FIELD3(1,IFFAT)+(-DIPJX*FER3+DOTJ*X)*SWF
         FIELD3(2,IFFAT)=FIELD3(2,IFFAT)+(-DIPJY*FER3+DOTJ*Y)*SWF
         FIELD3(3,IFFAT)=FIELD3(3,IFFAT)+(-DIPJZ*FER3+DOTJ*Z)*SWF
         FIELD3(1,JFFAT)=FIELD3(1,JFFAT)+(-DIPIX*FER3+DOTI*X)*SWF
         FIELD3(2,JFFAT)=FIELD3(2,JFFAT)+(-DIPIY*FER3+DOTI*Y)*SWF
         FIELD3(3,JFFAT)=FIELD3(3,JFFAT)+(-DIPIZ*FER3+DOTI*Z)*SWF
 300  CONTINUE
      ENDDO
 301  CONTINUE
C
C
C     -- FIELD3: FIELD AT POL DUE TO REACTION FIELD OF INDUCED DIPOLES
C                DIRECT METHOD (IMAGE CHARGE)
C
      IF(ISPHSOL.EQ.1 .AND. RSPHSOL.LT.1.0D+30) THEN
         SCALE = (RXNEPS-ONE)/RXNEPS
         R     = RSPHSOL
         R2    = R*R
         R3    = R*R2
         R4    = R2*R2
         IPCOUNT = ME - 1
         DO 400 IFFAT = 1, NFFAT
            FDIMGX = ZERO
            FDIMGY = ZERO
            FDIMGZ = ZERO
            XI     = CORD(1,IFFAT) - CENTX
            YI     = CORD(2,IFFAT) - CENTY
            ZI     = CORD(3,IFFAT) - CENTZ
            RI2    = XI*XI+YI*YI+ZI*ZI
            IF(RI2.GT.(R-1.0D+00)**2) GOTO 400
            IF(GOPARR) THEN
               IPCOUNT = IPCOUNT + 1
               IF(MOD(IPCOUNT,NPROC).NE.0) GOTO 400
            END IF
            DIPIX  = DIP(1,IFFAT)
            DIPIY  = DIP(2,IFFAT)
            DIPIZ  = DIP(3,IFFAT)
            DO 410 JFFAT = 1,NFFAT
               XJ     = CORD(1,JFFAT) - CENTX
               YJ     = CORD(2,JFFAT) - CENTY
               ZJ     = CORD(3,JFFAT) - CENTZ
               RJ2    = XJ*XJ+YJ*YJ+ZJ*ZJ
               IF(RJ2.GT.(R-1.0D+00)**2) GOTO 410
               DIPJX  = DIP(1,JFFAT)
               DIPJY  = DIP(2,JFFAT)
               DIPJZ  = DIP(3,JFFAT)
               DJRJ   = DIPJX*XJ+DIPJY*YJ+DIPJZ*ZJ
               DJRI   = DIPJX*XI+DIPJY*YI+DIPJZ*ZI
               RIRJ   = XI*XJ+YI*YJ+ZI*ZJ
               FAC2   = ONE/(RI2*RJ2 + R4 - TWO*R2*RIRJ)
               FAC    = SQRT(FAC2)
               FAC3   = FAC*FAC2
               FAC5   = FAC3*FAC2
               FDIMGX = FDIMGX-SCALE*(TWO*R*DJRJ*XI-R3*DIPJX)*FAC3
     *                        +SCALE*(R*RI2*DJRJ - R3*DJRI)*FAC5
     *                                     *THREE*(RJ2*XI - R2*XJ)
               FDIMGY = FDIMGY-SCALE*(TWO*R*DJRJ*YI-R3*DIPJY)*FAC3
     *                        +SCALE*(R*RI2*DJRJ - R3*DJRI)*FAC5
     *                                     *THREE*(RJ2*YI - R2*YJ)
               FDIMGZ = FDIMGZ-SCALE*(TWO*R*DJRJ*ZI-R3*DIPJZ)*FAC3
     *                        +SCALE*(R*RI2*DJRJ - R3*DJRI)*FAC5
     *                                     *THREE*(RJ2*ZI - R2*ZJ)
 410        CONTINUE
            FIELD3(1,IFFAT) = FIELD3(1,IFFAT) + FDIMGX
            FIELD3(2,IFFAT) = FIELD3(2,IFFAT) + FDIMGY
            FIELD3(3,IFFAT) = FIELD3(3,IFFAT) + FDIMGZ
 400     CONTINUE
      END IF
C
C     -- FIELD3: FIELD AT POL DUE TO REACTION FIELD OF INDUCED DIPOLES
C                (SPHERICAL BOUNDARY SURFACE CHARGE METHOD)
C        POT2  : POTENTIAL AT SURFACE DUE TO INDUCED DIPOLES
C        QRXN2 : SURFACE CHARGE DUE TO INDUCED DIPOLES
C
      IF(ISPHSOL.GE.60 .AND. RSPHSOL.LT.1.0D+30) THEN
         SCALE = (RXNEPS-ONE)/RXNEPS
         IPCOUNT = ME - 1
         DO 500 ITS = 1, NTS
            POT2(ITS) = ZERO
            XI     = XTS(ITS)
            YI     = YTS(ITS)
            ZI     = ZTS(ITS)
            DO 510 JFFAT = 1, NFFAT
               IF(GOPARR) THEN
                  IPCOUNT = IPCOUNT + 1
                  IF(MOD(IPCOUNT,NPROC).NE.0) GOTO 510
               END IF
               XJ    = CORD(1,JFFAT)
               YJ    = CORD(2,JFFAT)
               ZJ    = CORD(3,JFFAT)
               QJ    = XX(LFFCHARG+JFFAT-1)
               DIPJX = DIP(1,JFFAT)
               DIPJY = DIP(2,JFFAT)
               DIPJZ = DIP(3,JFFAT)
               X     = XI - XJ
               Y     = YI - YJ
               Z     = ZI - ZJ
               R2    = X*X + Y*Y + Z*Z
               ONER2 = ONE/R2
               ONER  = SQRT(ONER2)
               ONER3 = ONER2*ONER
               POT2(ITS)=POT2(ITS)+(DIPJX*X+DIPJY*Y+DIPJZ*Z)*ONER3
 510        CONTINUE
 500     CONTINUE
         IF(GOPARR) CALL DDI_GSUMF(2405,POT2,NTS)
         IPCOUNT = ME - 1
         DO 520 ITS = 1, NTS
            QRXN2(ITS) = ZERO
            IF(GOPARR) THEN
               IPCOUNT = IPCOUNT + 1
               IF(MOD(IPCOUNT,NPROC).NE.0) GOTO 520
            END IF
            DO JTS = 1, NTS
               QRXN2(ITS) = QRXN2(ITS) - CMAT1(ITS,JTS)*POT2(JTS)
            ENDDO
            QRXN2(ITS) = QRXN2(ITS)*SCALE
 520     CONTINUE
         IF(GOPARR) CALL DDI_GSUMF(2404,QRXN2,NTS)
         IPCOUNT = ME - 1
         DO 600 IFFAT=1,NFFAT
            XI = CORD(1,IFFAT)
            YI = CORD(2,IFFAT)
            ZI = CORD(3,IFFAT)
            DO 610 JTS=1, NTS
               IF(GOPARR) THEN
                  IPCOUNT = IPCOUNT + 1
                  IF(MOD(IPCOUNT,NPROC).NE.0) GOTO 610
               END IF
               XJ   = XTS(JTS)
               YJ   = YTS(JTS)
               ZJ   = ZTS(JTS)
               QJ   = QRXN1(JTS) + QRXN2(JTS)
               X    = XI - XJ
               Y    = YI - YJ
               Z    = ZI - ZJ
               R2   = X*X + Y*Y + Z*Z
               ONER2= ONE/R2
               ONER = SQRT(ONER2)
               ONER3= ONER2*ONER
               DUM  = QJ*ONER3
               FIELD3(1,IFFAT)=FIELD3(1,IFFAT)+DUM*X
               FIELD3(2,IFFAT)=FIELD3(2,IFFAT)+DUM*Y
               FIELD3(3,IFFAT)=FIELD3(3,IFFAT)+DUM*Z
 610        CONTINUE
 600     CONTINUE
      END IF
C
      IF(GOPARR) CALL DDI_GSUMF(2406,FIELD3,3*NFFAT)
C
C     - UPDATE DIPOLES AND COMPUTE POLARIZATION ENERGY -
      ENPOL = ZERO
      DO IFFAT=1,NFFAT
         DIP(1,IFFAT)=POL(IFFAT)*(FIELD1(1,IFFAT)
     *                           +FIELD2(1,IFFAT)
     *                           +FIELD3(1,IFFAT))
         DIP(2,IFFAT)=POL(IFFAT)*(FIELD1(2,IFFAT)
     *                           +FIELD2(2,IFFAT)
     *                           +FIELD3(2,IFFAT))
         DIP(3,IFFAT)=POL(IFFAT)*(FIELD1(3,IFFAT)
     *                           +FIELD2(3,IFFAT)
     *                           +FIELD3(3,IFFAT))
         SUM=SUM+ABS(DIP(1,IFFAT))+ABS(DIP(2,IFFAT))+ABS(DIP(3,IFFAT))
C        -- ONLY EXTRENAL FIELD SHOULD BE USED
         ENPOL = ENPOL - FIELD1(1,IFFAT)*DIP(1,IFFAT)
     *                 - FIELD1(2,IFFAT)*DIP(2,IFFAT)
     *                 - FIELD1(3,IFFAT)*DIP(3,IFFAT)
      ENDDO
      IF(ITER.EQ.1) CALL DCOPY(3*NFFAT,DIP,1,XX(LFFDIPOLD),1)
      DO IFFAT=1,NFFAT
         DIP(1,IFFAT)=0.75D+00*DIP(1,IFFAT)+
     *                0.25D+00*XX(LFFDIPOLD+(IFFAT-1)*3  )
         DIP(2,IFFAT)=0.75D+00*DIP(2,IFFAT)+
     *                0.25D+00*XX(LFFDIPOLD+(IFFAT-1)*3+1)
         DIP(3,IFFAT)=0.75D+00*DIP(3,IFFAT)+
     *                0.25D+00*XX(LFFDIPOLD+(IFFAT-1)*3+2)
      ENDDO
      CALL DCOPY(3*NFFAT,DIP,1,XX(LFFDIPOLD),1)
      ENPOL    = PT5*ENPOL
C
      IF(ITER.GT.2) THEN
         DIFF = ABS(SUM - SUMOLD)
         IF(DIFF.LT.POLTOL)THEN
C           IF(MASWRK) WRITE(IW,'(A,F16.14,A,I4,A)')
C    *      'INDUCED DIPOLE CONVERGED TO ',POLTOL,' IN ',ITER,' STEPS.'
            GOTO 306
         END IF
      END IF
      SUMOLD=SUM
      IF(IDOPOL.EQ.1) GOTO 306
 305  CONTINUE
      IF(MASWRK) WRITE(IW,'(A,F16.14,A,I4,A)')
     *'INDUCED DIPOLE NOT CONVERGED TO ',POLTOL,' IN ',IDOPOL,' STEPS.'
      IF(MASWRK) WRITE(IW,'(1X,A,F20.16)')'DIFF=',DIFF
      IF(MASWRK) WRITE(IW,*)'THIS IS A PROBLEM ONLY WHEN THIS MESSAGE',
     *                      ' APPEARS FREQUENTLY'
 306  CONTINUE
C
      RETURN
      END
C*MODULE QUANPOB  *DECK INIVEL
!>
!> @brief    guess initial velocity
!>
!> @author   Nandun Thellamurege, Hui Li
!>           - Jan 2011
!>
!> @details  assign initial velocity to all QM and MM atoms
!>
      SUBROUTINE INIVEL(VEL,QMVEL,ONEMAS,QM1MAS,CORD,ZMAS,QMZMAS,
     *                  SET,DX,DY,DZ,LISTQM)
      use mx_limits, only: mxatm
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      DIMENSION VEL(3,*),QMVEL(3,*),ONEMAS(*),QM1MAS(*),SET(*),
     *          DX(*),DY(*),DZ(*),ZMAS(*),QMZMAS(*),CORD(3,*),
     *          TIMAT(3,3),LISTQM(*)
C
      PARAMETER (TOKELVIN=3.15774646D+05)
      PARAMETER (BOLTZK=1.0D+00/TOKELVIN)
      PARAMETER (PT5=0.50D+00)
      PARAMETER (TWO=2.0D+00)
      PARAMETER (ZERO=0.0D+00)
C
      COMMON /FFMDPA/ DT,DT2,TEMP0,PRES0,POLTOL,VIR(3),PMEAN,VOLAV,
     *                ENPAV,ENKAV,TEMPAV,BERENDT,BERENDP,VELMAX,
     *                PMEANX,PMEANY,PMEANZ,NSTEP,KMASTER,KOUTACT(2),
     *                IHESS,INTALG,ITSTAT,IPSTAT,JOUT,KOUT,LOUT
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
      COMMON /FFRATT/ RATOLC,RATOLV,SCALRAT,VIRRAT(3),IRATTLE,JRATTLE,
     *                NRATTLE,MXRATT,LFFOLDCORD,LFFLSTRAT,LFFDSTRAT,
     *                LFFVELSV,IRATQM
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
C
C     NANDUN THELLAMUREGE, HUI LI, JAN 2011, LINCOLN
C
C     -- FOR U1,U2 UNIFORMLY DISTRIB ON (0,1)
C        X = SQRT(-2*LN(U1)) * COS(2*PI*U2)
C        GIVES NORMALLY DISTRIB X ON (-INF,+INF)
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
      PI     = 3.14159265358979323846264338D+00
      DO IDIM = 1, (NFFAT+NAT)*3
         CALL FFRAND(U1)
         CALL FFRAND(U2)
         SET(IDIM) = SQRT(-2.0D+00*LOG(U1)) * COS(2.0D+00*PI*U2)
      ENDDO
C
C     -- CALCULATE COM --
C
      COMX =ZERO
      COMY =ZERO
      COMZ =ZERO
      SMAS =ZERO
      DO IAT=1,NAT
         COMX =COMX +C(1,IAT)*QMZMAS(IAT)
         COMY =COMY +C(2,IAT)*QMZMAS(IAT)
         COMZ =COMZ +C(3,IAT)*QMZMAS(IAT)
         SMAS =SMAS +QMZMAS(IAT)
      ENDDO
      DO IFFAT=1,NFFAT
         IF(LISTQM(IFFAT).EQ.0.OR.LISTQM(IFFAT).GT.NAT)THEN
            COMX =COMX +CORD(1,IFFAT)*ZMAS(IFFAT)
            COMY =COMY +CORD(2,IFFAT)*ZMAS(IFFAT)
            COMZ =COMZ +CORD(3,IFFAT)*ZMAS(IFFAT)
            SMAS =SMAS +ZMAS(IFFAT)
         END IF
      ENDDO
      ONESMAS=1.0D+00/SMAS
      COMX = COMX*ONESMAS
      COMY = COMY*ONESMAS
      COMZ = COMZ*ONESMAS
C
C     -- FIRST ASSIGNMENT --
C
      SUMTX = ZERO
      SUMTY = ZERO
      SUMTZ = ZERO
      SUMAX = ZERO
      SUMAY = ZERO
      SUMAZ = ZERO
      J=0
      DO IAT=1,NAT
         SIGMA=SQRT(BOLTZK*TEMP0*QM1MAS(IAT))
         QMVEL(1,IAT)=SIGMA*SET(J+1)
         QMVEL(2,IAT)=SIGMA*SET(J+2)
         QMVEL(3,IAT)=SIGMA*SET(J+3)
         J=J+3
         DX(IAT) = C(1,IAT)-COMX
         DY(IAT) = C(2,IAT)-COMY
         DZ(IAT) = C(3,IAT)-COMZ
         SUMTX=SUMTX+QMVEL(1,IAT)*QMZMAS(IAT)
         SUMTY=SUMTY+QMVEL(2,IAT)*QMZMAS(IAT)
         SUMTZ=SUMTZ+QMVEL(3,IAT)*QMZMAS(IAT)
      ENDDO
      DO IFFAT=1,NFFAT
         IF(LISTQM(IFFAT).EQ.0.OR.LISTQM(IFFAT).GT.NAT)THEN
            SIGMA=SQRT(BOLTZK*TEMP0*ONEMAS(IFFAT))
            VEL(1,IFFAT)=SIGMA*SET(J+1)
            VEL(2,IFFAT)=SIGMA*SET(J+2)
            VEL(3,IFFAT)=SIGMA*SET(J+3)
            J=J+3
            DX(IFFAT+NAT) = CORD(1,IFFAT)-COMX
            DY(IFFAT+NAT) = CORD(2,IFFAT)-COMY
            DZ(IFFAT+NAT) = CORD(3,IFFAT)-COMZ
            SUMTX=SUMTX+VEL(1,IFFAT)*ZMAS(IFFAT)
            SUMTY=SUMTY+VEL(2,IFFAT)*ZMAS(IFFAT)
            SUMTZ=SUMTZ+VEL(3,IFFAT)*ZMAS(IFFAT)
         END IF
      ENDDO
C
C     -- MAKE TRANS AND ROT MOMENTA = ZERO --
C        (ZHOU ET AL BIOPHYSICAL JOURNAL,79,2902-2908)
C
C     -- VELOCITY OF COM --
      SUMTX=SUMTX*ONESMAS
      SUMTY=SUMTY*ONESMAS
      SUMTZ=SUMTZ*ONESMAS
C
      SUMAX = ZERO
      SUMAY = ZERO
      SUMAZ = ZERO
      TXX   = ZERO
      TYY   = ZERO
      TZZ   = ZERO
      TXY   = ZERO
      TXZ   = ZERO
      TYZ   = ZERO
      DO IAT=1,NAT
         QMVEL(1,IAT)=QMVEL(1,IAT)-SUMTX
         QMVEL(2,IAT)=QMVEL(2,IAT)-SUMTY
         QMVEL(3,IAT)=QMVEL(3,IAT)-SUMTZ
C        -- FIND ANGULAR MOMENTUM AROUND COM
         ANGVX=DY(IAT)*QMVEL(3,IAT)-DZ(IAT)*QMVEL(2,IAT)
         ANGVY=DZ(IAT)*QMVEL(1,IAT)-DX(IAT)*QMVEL(3,IAT)
         ANGVZ=DX(IAT)*QMVEL(2,IAT)-DY(IAT)*QMVEL(1,IAT)
         SUMAX=SUMAX+ANGVX*QMZMAS(IAT)
         SUMAY=SUMAY+ANGVY*QMZMAS(IAT)
         SUMAZ=SUMAZ+ANGVZ*QMZMAS(IAT)
C        -- CALCULATE MOMENT OF INTERTIA AROUND COM
         TXX = TXX + QMZMAS(IAT)*(DY(IAT)*DY(IAT)+DZ(IAT)*DZ(IAT))
         TYY = TYY + QMZMAS(IAT)*(DX(IAT)*DX(IAT)+DZ(IAT)*DZ(IAT))
         TZZ = TZZ + QMZMAS(IAT)*(DX(IAT)*DX(IAT)+DY(IAT)*DY(IAT))
         TXY = TXY - QMZMAS(IAT)* DX(IAT)*DY(IAT)
         TXZ = TXZ - QMZMAS(IAT)* DX(IAT)*DZ(IAT)
         TYZ = TYZ - QMZMAS(IAT)* DY(IAT)*DZ(IAT)
      ENDDO
      DO IFFAT=1,NFFAT
         IF(LISTQM(IFFAT).EQ.0.OR.LISTQM(IFFAT).GT.NAT)THEN
            VEL(1,IFFAT)=VEL(1,IFFAT)-SUMTX
            VEL(2,IFFAT)=VEL(2,IFFAT)-SUMTY
            VEL(3,IFFAT)=VEL(3,IFFAT)-SUMTZ
C           --- FIND ANGULAR MOMENTUM AROUND COM ---
            ANGVX=DY(IFFAT+NAT)*VEL(3,IFFAT)-DZ(IFFAT+NAT)*VEL(2,IFFAT)
            ANGVY=DZ(IFFAT+NAT)*VEL(1,IFFAT)-DX(IFFAT+NAT)*VEL(3,IFFAT)
            ANGVZ=DX(IFFAT+NAT)*VEL(2,IFFAT)-DY(IFFAT+NAT)*VEL(1,IFFAT)
            SUMAX=SUMAX+ANGVX*ZMAS(IFFAT)
            SUMAY=SUMAY+ANGVY*ZMAS(IFFAT)
            SUMAZ=SUMAZ+ANGVZ*ZMAS(IFFAT)
C           --- CALCULATE MOMENT OF INTERTIA AROUND COM ---
            TXX = TXX + ZMAS(IFFAT)*(DY(IFFAT+NAT)*DY(IFFAT+NAT)
     *                              +DZ(IFFAT+NAT)*DZ(IFFAT+NAT))
            TYY = TYY + ZMAS(IFFAT)*(DX(IFFAT+NAT)*DX(IFFAT+NAT)
     *                              +DZ(IFFAT+NAT)*DZ(IFFAT+NAT))
            TZZ = TZZ + ZMAS(IFFAT)*(DX(IFFAT+NAT)*DX(IFFAT+NAT)
     *                              +DY(IFFAT+NAT)*DY(IFFAT+NAT))
            TXY = TXY - ZMAS(IFFAT)* DX(IFFAT+NAT)*DY(IFFAT+NAT)
            TXZ = TXZ - ZMAS(IFFAT)* DX(IFFAT+NAT)*DZ(IFFAT+NAT)
            TYZ = TYZ - ZMAS(IFFAT)* DY(IFFAT+NAT)*DZ(IFFAT+NAT)
         END IF
      ENDDO
C
C     --- CALCULATE INVERSE INERTIA TENSOR ---
      IF(NFFAT.GE.3) THEN
      CALL TINV(TXX,TYY,TZZ,TXY,TXZ,TYZ,TIMAT)
C     --- MULTIPLY INVERSE INERTIA TENSOR WITH ANGULAR MOMENTUM ---
      PRODX = SUMAX*TIMAT(1,1) + SUMAY*TIMAT(1,2) + SUMAZ*TIMAT(1,3)
      PRODY = SUMAX*TIMAT(2,1) + SUMAY*TIMAT(2,2) + SUMAZ*TIMAT(2,3)
      PRODZ = SUMAX*TIMAT(3,1) + SUMAY*TIMAT(3,2) + SUMAZ*TIMAT(3,3)
      ELSE
      PRODX = ZERO
      PRODY = ZERO
      PRODZ = ZERO
      END IF
C
C     -- SECOND ASSIGNMENT --
C
      ENKIN = ZERO
      DO IAT =1 ,NAT
C        -- GET THE CROSS PRODUCT WITH THE DISTANCES --
         CROSX = PRODY*DZ(IAT) - PRODZ*DY(IAT)
         CROSY = PRODZ*DX(IAT) - PRODX*DZ(IAT)
         CROSZ = PRODX*DY(IAT) - PRODY*DX(IAT)
C        -- MAKE IT (3N-6) --
         QMVEL(1,IAT)=QMVEL(1,IAT)-CROSX
         QMVEL(2,IAT)=QMVEL(2,IAT)-CROSY
         QMVEL(3,IAT)=QMVEL(3,IAT)-CROSZ
C        -- CALCULATE CURRENT TEMPERATURE --
         DUMY = ZERO
         DO III = 1, 3
            DUMY = DUMY + QMVEL(III,IAT)*QMVEL(III,IAT)
         ENDDO
         ENKIN = ENKIN + DUMY*QMZMAS(IAT)
      ENDDO
      DO IFFAT =1 ,NFFAT
         IF(LISTQM(IFFAT).EQ.0.OR.LISTQM(IFFAT).GT.NAT)THEN
C           -- GET THE CROSS PRODUCT WITH THE DISTANCES --
            CROSX = PRODY*DZ(IFFAT+NAT) - PRODZ*DY(IFFAT+NAT)
            CROSY = PRODZ*DX(IFFAT+NAT) - PRODX*DZ(IFFAT+NAT)
            CROSZ = PRODX*DY(IFFAT+NAT) - PRODY*DX(IFFAT+NAT)
C           -- MAKE IT (3N-6) --
            VEL(1,IFFAT)=VEL(1,IFFAT)-CROSX
            VEL(2,IFFAT)=VEL(2,IFFAT)-CROSY
            VEL(3,IFFAT)=VEL(3,IFFAT)-CROSZ
C           -- CALCULATE CURRENT TEMPERATURE --
            DUMY = ZERO
            DO III = 1, 3
               DUMY = DUMY + VEL(III,IFFAT)*VEL(III,IFFAT)
            ENDDO
            ENKIN = ENKIN + DUMY*ZMAS(IFFAT)
         END IF
      ENDDO
      ENKIN = PT5*ENKIN
      TEMP  = TOKELVIN*TWO*ENKIN/NDOF
C
C     -- THIRD ASSIGNMENT --
C        RESCALING VELOCITIES TO GIVE TEMPERATURE AS TEMP0
C
      TFACT =SQRT(TEMP0/TEMP)
      TMOMX = ZERO
      TMOMY = ZERO
      TMOMZ = ZERO
      AMOMX = ZERO
      AMOMY = ZERO
      AMOMZ = ZERO
      ENKIN = ZERO
      DO IAT=1,NAT
         QMVEL(1,IAT)=QMVEL(1,IAT)*TFACT
         QMVEL(2,IAT)=QMVEL(2,IAT)*TFACT
         QMVEL(3,IAT)=QMVEL(3,IAT)*TFACT
         TMOMX= TMOMX+ QMVEL(1,IAT)*QMZMAS(IAT)
         TMOMY= TMOMY+ QMVEL(2,IAT)*QMZMAS(IAT)
         TMOMZ= TMOMZ+ QMVEL(3,IAT)*QMZMAS(IAT)
         AMOMX= AMOMX+ QMVEL(3,IAT)*DY(IAT)*QMZMAS(IAT)
     *               - QMVEL(2,IAT)*DZ(IAT)*QMZMAS(IAT)
         AMOMY= AMOMY+ QMVEL(1,IAT)*DZ(IAT)*QMZMAS(IAT)
     *               - QMVEL(3,IAT)*DX(IAT)*QMZMAS(IAT)
         AMOMZ= AMOMZ+ QMVEL(2,IAT)*DX(IAT)*QMZMAS(IAT)
     *               - QMVEL(1,IAT)*DY(IAT)*QMZMAS(IAT)
         DUMY = ZERO
         DO III = 1, 3
            DUMY = DUMY + QMVEL(III,IAT)*QMVEL(III,IAT)
         ENDDO
         ENKIN = ENKIN + DUMY*QMZMAS(IAT)
      ENDDO
      DO IFFAT=1,NFFAT
         IF(LISTQM(IFFAT).EQ.0.OR.LISTQM(IFFAT).GT.NAT)THEN
            VEL(1,IFFAT)=VEL(1,IFFAT)*TFACT
            VEL(2,IFFAT)=VEL(2,IFFAT)*TFACT
            VEL(3,IFFAT)=VEL(3,IFFAT)*TFACT
            TMOMX= TMOMX+ VEL(1,IFFAT)*ZMAS(IFFAT)
            TMOMY= TMOMY+ VEL(2,IFFAT)*ZMAS(IFFAT)
            TMOMZ= TMOMZ+ VEL(3,IFFAT)*ZMAS(IFFAT)
            AMOMX= AMOMX+ VEL(3,IFFAT)*DY(IFFAT+NAT)*ZMAS(IFFAT)
     *                  - VEL(2,IFFAT)*DZ(IFFAT+NAT)*ZMAS(IFFAT)
            AMOMY= AMOMY+ VEL(1,IFFAT)*DZ(IFFAT+NAT)*ZMAS(IFFAT)
     *                  - VEL(3,IFFAT)*DX(IFFAT+NAT)*ZMAS(IFFAT)
            AMOMZ= AMOMZ+ VEL(2,IFFAT)*DX(IFFAT+NAT)*ZMAS(IFFAT)
     *                  - VEL(1,IFFAT)*DY(IFFAT+NAT)*ZMAS(IFFAT)
            DUMY = ZERO
            DO III = 1, 3
               DUMY = DUMY + VEL(III,IFFAT)*VEL(III,IFFAT)
            ENDDO
            ENKIN = ENKIN + DUMY*ZMAS(IFFAT)
         END IF
      ENDDO
      ENKIN = PT5*ENKIN
      TEMP  = TOKELVIN*TWO*ENKIN/NDOF
C     WRITE(IW,*)'TMOMX=',TMOMX,TMOMY,TMOMZ
C     WRITE(IW,*)'AMOMX=',AMOMX,AMOMY,AMOMZ
C
C     -- SOME ATOMS HAVE NOT BEEN GIVEN VELOCITY --
      IF(LISTQM(NFFAT+NAT+1).GT.0)THEN
         DO IAT = 1, NAT
            IFFAT = LISTQM(NFFAT+IAT)
            IF(IFFAT.GT.0) THEN
               VEL(1,IFFAT)   = QMVEL(1,IAT)
               VEL(2,IFFAT)   = QMVEL(2,IAT)
               VEL(3,IFFAT)   = QMVEL(3,IAT)
            END IF
         ENDDO
      END IF
C
      RETURN
      END
C*MODULE QUANPOB  *DECK TINV
      SUBROUTINE TINV(TXX,TYY,TZZ,TXY,TXZ,TYZ,TIMAT)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      DIMENSION TEMP(3,3),AUGMAT(3,6),TIMAT(3,3)
C
C     NANDUN THELLAMUREGE, JAN 2011, LINCOLN
C
C     --- STORE INPUT MATRIX ---
      TEMP(1,1) = TXX
      TEMP(2,2) = TYY
      TEMP(3,3) = TZZ
      TEMP(1,2) = TXY
      TEMP(1,3) = TXZ
      TEMP(2,3) = TYZ
      TEMP(2,1) = TEMP(1,2)
      TEMP(3,1) = TEMP(1,3)
      TEMP(3,2) = TEMP(2,3)
C     --- AUGMENT INPUT MATRIX WITH AN IDENTITY MATRIX ---
      DO I = 1,3
         DO J = 1,6
            IF (J.LE.3) THEN
               AUGMAT(I,J) = TEMP(I,J)
            ELSE
               IF (I+3.EQ.J) THEN
                  AUGMAT(I,J) = 1.0D+00
               ELSE
                  AUGMAT(I,J) = 0.0D+00
               END IF
            END IF
         ENDDO
      ENDDO
C     --- REDUCE AUGMENTED MATRIX TO UPPER TRIANGULAR FORM ---
      DO K=1,2
         IF(AUGMAT(K,K).EQ.0.0D+00) THEN
            DO I = K+1,3
               IF(AUGMAT(I,K).NE.0.0D+00) THEN
                  DO J=1,6
                     AUGMAT(K,J)=AUGMAT(K,J)+AUGMAT(I,J)
                  ENDDO
               END IF
             ENDDO
         END IF
         DO J=K+1,3
            VAR = AUGMAT(J,K)/AUGMAT(K,K)
            DO I=K,6
               AUGMAT(J,I)=AUGMAT(J,I) - VAR*AUGMAT(K,I)
            ENDDO
         ENDDO
      ENDDO
C     --- MAKE DIAGONAL ELEMENTS 1.0D+00 ---
      DO I=1,3
         VAR = AUGMAT(I,I)
         DO J= I,6
            AUGMAT(I,J)=AUGMAT(I,J)/VAR
         ENDDO
      ENDDO
C     -- REDUCE RIGHT HALF OF THE AUGMENTED MATRIX TO IDENTITY MATRIX --
      DO K = -1,2
      DO I = 1,K
         VAR = AUGMAT(I,K+1)
         DO J = K,6
            AUGMAT(I,J) = AUGMAT(I,J) - AUGMAT(K+1,J)*VAR
         ENDDO
      ENDDO
      ENDDO
C     --- STORE THE INVERTED MATRIX ---
      DO I = 1,3
         DO J=1,3
            TIMAT(I,J) = AUGMAT(I,J+3)
         ENDDO
      ENDDO
C
      RETURN
      END
C*MODULE QUANPOB  *DECK FFMDX
!>
!> @brief    main driver for QuanPol MD
!>
!> @author   Hui Li
!>           - Jan 2011
!>
!> @details  main driver for QuanPol MD
!>
      SUBROUTINE FFMDX
      use mx_limits, only: mxatm
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      LOGICAL GOPARR,DSKWRK,MASWRK
C
      COMMON /FFDFS / TIMDFS,QDION,AMION,TEFF,NDFS,NATMGAS,
     *                LFFDFSC,
     *                LFFDFSC0,LFFDFSA,LFFDFSN,LFFDFCOM,KDFS,LFFDFSCAV
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
      COMMON /FFMDOP/ LFFCORDG,LFFCORDGSV,LFFCORDGSV2,LFFCORDGSVQ,
     *                LFFFFGRDG0,LFFFFGRDG1,LFFFFGRDG2,MDOPT
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
      COMMON /FFRATT/ RATOLC,RATOLV,SCALRAT,VIRRAT(3),IRATTLE,JRATTLE,
     *                NRATTLE,MXRATT,LFFOLDCORD,LFFLSTRAT,LFFDSTRAT,
     *                LFFVELSV,IRATQM
      COMMON /FFRDF / DELRDF,NUMGRD,NUMSUM,NRDF,LFFGOFR,
     *                LFFNFRAG1,LFFNFRAG2,LFFFRAG1,LFFFRAG2,
     *                NRDEN,NBINRDEN,LFFPRO,LFFNRDPRATM,
     *                LFFRDPRATM,LFFDIESTEP,LFFDI1STEP
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
C
C     HUI LI, JAN 2011, LINCOLN
C
      IF(NRDF.GT.0) CALL INIRDF(X(LFFATMNAM),X(LFFFRAG1),X(LFFFRAG2),
     *                   X(LFFNFRAG1),X(LFFNFRAG2),X(LFFGOFR))
      IF(NDFS.GT.0) CALL INIDFS(X(LFFATMNAM),X(LFFDFSA),
     *                   X(LFFDFSC),X(LFFDFSN),X(LFFCHARG),X(LFFZMAS))
      NTS = ISPHSOL
      IF(NAT.EQ.0 .AND. NFFAT.GT.0) THEN
        CALL FFMD1(X(LFFATMNAM),X(LFFCORD),X(LFFCORDSV),X(LFFZANF),
     *             X(LFFZMAS),X(LFFONEMAS),X(LFFQMZMAS),X(LFFQM1MAS),
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
     *             X(LFFVEL),X(LFFQMVEL),X(LFFFCLJTP),X(LFFNTYPE),
     *             X(LFFFFGRD0),X(LFFFFGRD1),X(LFFFFGRD2),
     *             X(LFFXTS),X(LFFYTS),X(LFFZTS),X(LFFCMAT1),
     *             X(LFFPOT1),X(LFFPOT2),X(LFFQRXN1),X(LFFQRXN2),NTS,
     *             X(LFFLISTQM),X(LFFNONLS1),
     *             X(LFFNONLSTQ),
     *             X(LFFMAPLST),X(LFFCMAPCO),
     *             X(LFFLSTCELL),X(LFFNONLS2),
     *             X(LFFCORDSV2),X(LFFCORDSVQ),
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
     *             X(LFFOLDCORD),X(LFFLSTRAT),X(LFFDSTRAT),
     *             X(LFFNONLSPMA),X(LFFL1213PMA),X(LFFL14PMA),
     *             X(LFF2CORD),X(LFFLSBONDPMA),X(LFFLSANGLPMA),
     *             X(LFFLSDIHRPMA),X(LFFLSDIHBPMA),
     *             X(LFFLSWAGGPMA),X(LFFLSCMAPPMA),
     *             X(LFFNONLSPMB),X(LFFL1213PMB),X(LFFL14PMB),
     *             X(LFFLSBONDPMB),X(LFFLSANGLPMB),
     *             X(LFFLSDIHRPMB),X(LFFLSDIHBPMB),
     *             X(LFFLSWAGGPMB),X(LFFLSCMAPPMB),
     *             X(LFFUMBHIS),X(LFFUM2HIS),
     *             X(LFFVELSV),X(LFFDFSC0),
     *             X(LFFCORDG),X(LFFCORDGSV),X(LFFCORDGSV2),
     *             X(LFFCORDGSVQ),X(LFFFFGRDG0),X(LFFFFGRDG1),
     *             X(LFFFFGRDG2))
      ELSE IF (NAT.GT.0 .AND. NFFAT.GT.0) THEN
        CALL FFMD2(X(LFFATMNAM),X(LFFCORD),X(LFFCORDSV),X(LFFZANF),
     *             X(LFFZMAS),X(LFFONEMAS),X(LFFQMZMAS),X(LFFQM1MAS),
     *             X(LFFCHARG),X(LFFVEL),X(LFFQMVEL),
     *             X(LFFFFGRD0),X(LFFFFGRD1),X(LFFFFGRD2),
     *             X(LFFQMGRD0),X(LFFQMGRD1),X(LFFQMGRD2),
     *             X(LFFLISTQM),X(LFFNONLS1),
     *             X(LFFNONLSTQ),
     *             X(LFFLSTCELL),X(LFFNONLS2),
     *             X(LFFCORDSV2),X(LFFCORDSVQ),
     *             X(LFFMVFASTS2),X(LFFMVFASTS3),X(LFFMVFASTS4),
     *             X(LFFMVFASTL2),X(LFFMVFASTL3),X(LFFMVFASTL4),
     *             X(LFFNONLSA),X(LFFNONLSB),
     *             X(LFFOLDCORD),X(LFFLSTRAT),X(LFFDSTRAT),
     *             X(LFFNONLSPMA),X(LFFNONLSPMB),X(LFF2CORD),
     *             X(LFFUMBHIS),X(LFFUM2HIS),
     *             X(LFFVELSV),
     *             X(LFFOLDC),X(LFFQMVELSV),X(LFFDFSC0))
      END IF
C
      IF(MASWRK.AND.NSTEP.GE.0) THEN
         WRITE(IW,*)' '
         WRITE(IW,*)'================ QUANPOL MD SIMULATION',
     *              ' SUCCESSFULLY COMPLETED ================='
         WRITE(IW,*)'========= RESTART $DATA, $FFDATA, $FFDATB',
     *              ' ARE IN .DAT AND .TRJ FILES =========='
         WRITE(IW,*)' '
      END IF
      CALL TIMIT(1)
C
      RETURN
      END
