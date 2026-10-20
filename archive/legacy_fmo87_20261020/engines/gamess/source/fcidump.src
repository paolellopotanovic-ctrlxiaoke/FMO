C*MODULE FCIDUMP  *DECK DRIVE_FCIDUMP
C     ------------------------------
C> @brief      This routine generates a FCIDUMP file containing one- and
C>             two-body molecular integrals and nuclear repulsion
C>             energies. For the specification see:
C>             https://doi.org/10.1016/0010-4655(89)90033-7
C>
C> @author     J. Emiliano Deustua
C>             -2020
C>
C> @date August 17, 2020-J. Emiliano Deustua
C> - Based on the ALDECI routine

      SUBROUTINE DRIVE_FCIDUMP()
C     ------------------------------
      USE MX_LIMITS,ONLY:mxatm
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      LOGICAL GOPARR,DSKWRK,MASWRK,DOEXCH,
     *        DDITRF,DOCORE,
     *        DOOOOO,DOVOOO,DOVVOO,DOVOVO,DOVVVO,DOVVVV
C
      COMMON /IOFILE/ IR,IW,IP,IS,IJKT,IDAF,NAV,IODA(950)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,MA,MB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
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
C
C        driver for generation of FCIDUMP files
C
C        ----- integral transformation -----
C
      DDITRF=GOPARR
      IF(DDITRF) THEN
        DOOOOO=.TRUE.
        DOVOOO=.FALSE.
        DOVVOO=.FALSE.
        DOVOVO=.FALSE.
        DOVVVO=.FALSE.
        DOVVVV=.FALSE.
        DOCORE=.FALSE.
        DOEXCH=.FALSE.
      ELSE
        DOOOOO=.TRUE.
        DOVOOO=.FALSE.
        DOVVOO=.FALSE.
        DOVOVO=.FALSE.
        DOVVVO=.FALSE.
        DOVVVV=.FALSE.
        DOCORE=.FALSE.
        DOEXCH=.FALSE.
      ENDIF

      NCOR = 0
      NACT = NQMT
      NORB = NQMT

      ! [TODO] Find out what these flags are. Might be only
      ! applicable to DDI (see gamess.src)
C      IF (DDITRF.AND.MASWRK) THEN
C        WRITE(IW, '(/2X,A/)')
C     * 'ERROR: DDI INTEGRAL TRANSFORMATION IS NOT AVAILABLE FOR FCIDUMP'
C        CALL ABRT()
C      END IF

      CALL TRFMCX(0,NCOR,NORB,NORB,.FALSE.,DOEXCH,
     *            DDITRF,DOOOOO,DOVOOO,DOVVOO,DOVOVO,
     *            DOVVVO,DOVVVV,DOCORE)
C
      CALL PRINT_FCIDUMP(2, DDITRF,0 )
      IF(DDITRF) CALL DDI_DESTROY(D_OOOO)
C
      RETURN
C
      END

C*MODULE FCIDUMP  *DECK PRINT_FCIDUMP
C     ------------------------------
C> @brief      This routine generates a FCIDUMP file containing one- and
C>             two-body molecular integrals.
C>
C> @author     J. Emiliano Deustua
C>             -2020
C>
C> @date August 17, 2020-J. Emiliano Deustua
C> - Based on the ALDECI routine

      SUBROUTINE PRINT_FCIDUMP(NPRINT,DDITRF,ICIMALMQ)
C
      USE MX_LIMITS,ONLY:mxrt,mxatm,mxfrg,mxdppt
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      LOGICAL PACK2E,GOPARR,DSKWRK,MASWRK,DDITRF
C
      COMMON /FRGINF/ NMPTS(MXFRG),NMTTPT,IEFC,IEFD,IEFQ,IEFO,
     *                NPPTS(MXFRG),NPTTPT,IEFP,
     *                NRPTS(MXFRG),NRTTPT,IREP,ICHGP,NFRG,
     *                NDPPTS(MXDPPT),NDPTTPT,IEFDP,LSTMPTS(MXFRG),
     *                NBSFN(MXFRG),NMXMO(MXFRG)
C
      COMMON /CIFILS/ NFT11,NFT12,NFT13,NFT14,NFT15,NFT16,IDAF20,NEMEMX
      COMMON /DETWFN/ WSTATE(MXRT),SPINS(MXRT),CRIT,PRTTOL,S,SZ,
     *                GRPDET,STSYM,GLIST,DWPARM,
     *                NFLGDM(MXRT),IWTS(MXRT),NCORSV,NCOR,NACT,NORB,
     *                NA,NB,K,KST,IROOT,IPURES,MAXW1,NITER,MAXP,NCI,
     *                IGPDET,KSTSYM,NFTGCI,IDWEIGH,
     *                fstate(mxrt),ifts(mxrt)
      COMMON /ENRGYS/ ENUCR,EELCT,ETOT,STOT,SZZ,ECORE,ESCF,EERD,E1,E2,
     *                VEN,VEE,EPOT,EKIN,ESTATE(MXRT),STATN,EDFT(2),EDISP
      COMMON /FMCOM / X(1)
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,MA,MB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /INTFIL/ NINTMX,NHEX,NTUPL,PACK2E,INTTYP,IGRDTYP
      COMMON /IOFILE/ IR,IW,IP,IS,IJKT,IDAF,NAV,IODA(950)
      COMMON /MACHIN/ NWDVAR,MAXFM,MAXSM,LIMFM,LIMSM
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
C
C     ----- driver for FCIDUMP integral printer -----
      IF(MASWRK) WRITE(IW, 9000)

C     Define orbital sets and sizes
      ! [TODO] make sure that NQMT is the right variable
      ECONST = ECORE + ENUCR
      NCORE = 0
      NTCO = 0
      NTOT = NQMT
      NORB = NQMT
      NOCC = MA
      NUNOCC = NORB - NOCC
      NCOR = 0
      NSYM = 2**IGPDET
      M1 = NQMT
      M2 = (M1*M1+M1)/2
      M4 = (M2*M2+M2)/2
      MXRTMQ = MXRT
C
C        integral buffers for distributed/disk file transformed ints
C
      IF(DDITRF) THEN
        NOCC  = NACT + NCORSV
        NOTR  = (NOCC*NOCC+NOCC)/2
        LENXX = NOTR
        LENIXX= 0
      ELSE
        LENXX = NINTMX
        LENIXX= NINTMX
      END IF

C
C        allocate memory for for loading integrals
C
      CALL VALFM(LOADFM)
      LSINT1 = LOADFM + 1
      LSINT2 = LSINT1 + M2
      LIA    = LSINT2 + M4
      LXX    = LIA    + M2/NWDVAR + 1
      LIXX   = LXX    + LENXX
      LAST   = LIXX   + MXRTMQ
      NEED2  = LAST - LOADFM - 1
      CALL GETFM(NEED2)
C
C     -- print 1 and 2 e- transformed integrals over active orbitals --
C
      CALL PRINT_INTEGRALS(DDITRF,IJKT,X(LSINT1),X(LSINT2),
     *            NCORE,M1,M2,M4,
     *            X(LIA),X(LXX),X(LIXX),NINTMX)

      CALL RETFM(NEED2)
C
      RETURN
C
 9000 FORMAT(/5X,50(1H-)/
     *       5X,'     FCIDUMP INTEGRAL FILE GENERATION'/
     *       5X,'     WRITTEN BY J. EMILIANO DEUSTUA'/
     *       5X,50(1H-))
 9110 FORMAT(/1X,'THE NUMBER OF DETERMINANTS HAVING SPACE SYMMETRY ',A3/
     *        1X,'IN POINT GROUP ',A4,' WITH SZ=',F5.1,' IS',I15)
 9120 FORMAT(1X,'WHICH INCLUDES',I15,' CSFS WITH S=',F5.1)
 9130 FORMAT(1X,'THE DETERMINANT FULL CI REQUIRES',I16,' WORDS')
 9020 FORMAT(1X,'SOLVATION MODEL REQUIRES ',7X,I12,
     *    ' WORDS')
      END

C*MODULE ALDECI  *DECK RDCI12
      SUBROUTINE PRINT_INTEGRALS(DDITRF,NFT,X1,X2,NCORE,M1,M2,M4,
     *                           IA,XX,IX,NINTMX)
C
      USE MX_LIMITS,ONLY:mxrt,mxatm
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      LOGICAL GOPARR,DSKWRK,MASWRK,DDITRF
C          non-dditrf needx xx(nintmx), and ix(nintmx d.p.)
C              dditrf needs xx(m2), and no ix array
      DIMENSION X1(M2),X2(M4),IA(M2),XX(*),IX(*)

       PARAMETER (ZERO=0.0D+00, ONE=1.0D+00)
C
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /PCKLAB/ LABSIZ
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,MA,MB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)

C EDS modification
      COMMON /ENRGYS/ ENUCR,EELCT,ETOT,STOT,SZZ,ECORE,ESCF,EERD,E1,E2,
     *                VEN,VEE,EPOT,EKIN,ESTATE(MXRT),STATN,EDFT(2),EDISP
C
C  DDI ARRAY HANDLES
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
      INTEGER         I, J, NOCC
C
C     -- read 1 and 2 e- transformed integrals into replicated memory --
C     Only integrals in the active space, between NCORE and NCORE+M1
C     are returned in X1 and X2 arrays.  The 2e- integrals might be
C     in distributed memory, depending on DDITRF flag.
C
      IROW = 0
      DO 110 I=1,M2
         IA(I) = IROW
         IROW = IROW+I
  110 CONTINUE
C
      CALL VCLR(X2,1,M4)
      CALL SEQOPN(297, 'FCIDUMP', 'NEW', .FALSE., 'FORMATTED')


C
      NOCC = MA
      NORB = NQMT
      NUNOCC = NORB - NOCC
      NCORE = 0

C     Write FCIDUMP header
      WRITE(297, 9200) NORB, NOCC*2, 0
      CALL PRINT_SYM(NOCC, NUNOCC, NUM)
      WRITE(297, 9201)
C
C
      IF (DDITRF) THEN
C
C         integrals are to be obtained from distributed memory
C         obtain the one electron integrals, always read from disk.
C
        CALL SEQREW(NFT)
        CALL SQREAD(NFT,X1,M2)
        CALL SEQREW(NFT)

        IF (MASWRK) THEN
          WRITE(IW,'(4X, A/)')
     *      'WRITING INTEGRALS FROM DDI TRANSFORMATION'
          NACT = M1
          NOCC = NACT + NCORE
          NOTR = (NOCC*NOCC+NOCC)/2
          DO I = 1, NACT
            IN = I + NCORE
            DO J = 1, I
              JN = J + NCORE
              IJ = (I*I-I)/2 + J
              IJN = (IN*IN-IN)/2 + JN
                CALL DDI_GET(D_OOOO,1,NOTR,IJN,IJN,XX)
                DO K = 1, NACT
                  KN = K + NCORE
                  DO L = 1, K
                    LN = L + NCORE
                    KL = (K*K-K)/2 + L
                    IF (IJ.GE.KL) THEN
                     KLN = (KN*KN-KN)/2 + LN
                     IJKL = (IJ*IJ-IJ)/2 + KL
                     VAL = XX(KLN)
                     IF (VAL.NE.ZERO) THEN
                       WRITE(297, 9202) VAL, I, J, K, L
                     END IF
                    END IF
                  END DO
                END DO
            END DO ! J
          END DO ! I

C         EDS modification to write integrals down
          DO I=1, NORB
            DO J=1, I
              IF (X1(I*(I-1)/2+J).NE.ZERO) THEN
                WRITE(297, 9202) X1(I*(I-1)/2+J), I, J, 0, 0
              ENDIF
            END DO
          ENDDO
          WRITE(297, 9202) ENUCR, 0, 0, 0, 0
        ENDIF

C
      ELSE

C
C         obtain the one electron integrals, always read from disk.
C         only the master has the 1e- integrals, but if the 2e-
C         integrals are on disk, all nodes must process them.
C
         CALL SEQREW(NFT)
         IF(MASWRK) THEN
           CALL SQREAD(NFT,X1,M2)
           WRITE(IW,'(4X, A/)')
     *       'WRITING INTEGRALS FROM SERIAL TRANSFORMATION'
         END IF
C
C         Read transformed 2e- integral file in reverse canonical order.
C
  200    CONTINUE
         LABSIZ=1
         CALL PREAD(NFT,XX,IX,NX,NINTMX)
         IF (NX.EQ.0) GO TO 240
         MX = ABS(NX)
         IF (MX.GT.NINTMX) THEN
           IF(MASWRK) WRITE(IW,*)
     *       'INTEGRAL CONFUSION IN FCIDUMP PRINT_INTEGRALS'
           CALL ABRT
         END IF

         DO 220 M = 1,MX
            VAL = XX(M)
            NPACK = M
            IF(LABSIZ .EQ. 2) THEN
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
C
C        note index reversals to convert from reverse canonical order
C
            K = IPACK
            L = JPACK
            I = KPACK
            J = LPACK

C
            IF(I.LE.0  .OR.  I.GT.M1) GO TO 220
            IF(J.LE.0  .OR.  J.GT.M1) GO TO 220
            IF(K.LE.0  .OR.  K.GT.M1) GO TO 220
            IF(L.LE.0  .OR.  L.GT.M1) GO TO 220

            IF (MASWRK.AND.VAL.NE.ZERO) THEN
              WRITE(297, 9202) VAL, I, J, K, L
            END IF
C
            IJ = IA(I)+J
            KL = IA(K)+L
            IJKL = IA(IJ) + KL
            X2(IJKL) = VAL
  220    CONTINUE


C        EDS modification to write integrals down
         IF (MASWRK) THEN
           DO I=1, NORB
             DO J=1, I
               IF (X1(I*(I-1)/2+J).NE.ZERO) THEN
                 WRITE(297, 9202) X1(I*(I-1)/2+J), I, J, 0, 0
               ENDIF
             END DO
           ENDDO
         ENDIF

         WRITE(297, 9202) ENUCR, 0, 0, 0, 0
         IF(NX.GT.0) GO TO 200
C
  240    CONTINUE
         CALL SEQREW(NFT)
      END IF
C
C  GLOBAL SUM ALSO ACTS AS A SYNC
C
      CALL DDI_GSUMF(2500,X2,M4)
      RETURN

 9200 FORMAT(1X,'&FCI NORB=', I0, ',NELEC=', I0, ',MS2=', I0)
 9201 FORMAT(2X,'ISYM=1', 1X, 'UHF=.FALSE.',/,
     *       1X,'&END')
 9202 FORMAT(D24.16, 4I4)

      END

C*MODULE FCIDUMP  *DECK PRINT_SYM
C     ------------------------------
C> @brief      This routine prints the orbital symmetries toe
C>             the FCIDUMP file.
C>
C> @author     J. Emiliano Deustua
C>             -2020
C>
C> @date August 17, 2020-J. Emiliano Deustua
      SUBROUTINE PRINT_SYM(NO,NU,L1)
C
      USE MX_LIMITS,ONLY:mxatm
      IMPLICIT DOUBLE PRECISION  (A-H,O-Z)
      LOGICAL PACK2E
C
      COMMON /FMCOM / X(1)
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /INTFIL/ NINTMX,NHEX,NTUPL,PACK2E,INTTYP,IGRDTYP
      COMMON /IOFILE/ IR,IW,IP,IJK,IJKT,IDAF,NAV,IODA(950)
      COMMON /MACHIN/ NWDVAR,MAXFM,MAXSM,LIMFM,LIMSM
      COMMON /RUNOPT/ RUNTYP,EXETYP,NEVALS

      LOGICAL ABELPT
      DIMENSION GRPS(8)
      DIMENSION LABMO(L1),ISORB(NO+NU)
      CHARACTER*3 GR,REP
      CHARACTER*4 LBTEMP
      DIMENSION NREX(4),MCX(8,8,4)
      DIMENSION LABREP(8)
      COMMON /EOMSYMXX/ IG,NRE(4),MC(8,8,4),REP(8)

      DATA GRPS/8HC1      ,8HCI      ,8HCS      ,8HC2      ,
     *          8HD2      ,8HC2V     ,8HC2H     ,8HD2H     /
      COMMON /SYMMOL/ GROUP,COMPLEX,IGROUP,NAXIS,ILABMO,ABEL

      DATA NREX/1,2,4,8/
      DATA MCx/
     &1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     &0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     &
     & 1,1,0,0,0,0,0,0, 1,-1,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,
     & 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,
     & 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,
     &
     & 1,1,1,1,0,0,0,0,      1,1,-1,-1,0,0,0,0, 1,-1,1,-1,0,0,0,0,
     & 1,-1,-1,1,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,
     & 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,
     &
     & 1,1,1,1,1,1,1,1,      1,1,1,1,-1,-1,-1,-1,  1,1,-1,-1,1,1,-1,-1,
     & 1,1,-1,-1,-1,-1,1,1,  1,-1,1,-1,1,-1,1,-1,  1,-1,1,-1,-1,1,-1,1,
     & 1,-1,-1,1,1,-1,-1,1,  1,-1,-1,1,-1,1,1,-1/
C
C
C
      GROUP_MOL  = GRPS(1)
      IF(IGROUP.EQ.1)                GROUP_MOL = GRPS(1)
      IF(IGROUP.EQ.3)                GROUP_MOL = GRPS(2)
      IF(IGROUP.EQ.2)                GROUP_MOL = GRPS(3)
      IF(IGROUP.EQ.4.AND.NAXIS.EQ.2) GROUP_MOL = GRPS(4)
      IF(IGROUP.EQ.8.AND.NAXIS.EQ.2) GROUP_MOL = GRPS(5)
      IF(IGROUP.EQ.7.AND.NAXIS.EQ.2) GROUP_MOL = GRPS(6)
      IF(IGROUP.EQ.6.AND.NAXIS.EQ.2) GROUP_MOL = GRPS(7)
      IF(IGROUP.EQ.9.AND.NAXIS.EQ.2) GROUP_MOL = GRPS(8)
C
C         OBTAIN ORBITAL SYMMETRY INFORMATION
C
      WRITE(UNIT=GR,FMT='(A3)') GROUP_MOL
      IF(.NOT.ABELPT()) GR='C1 '
C
      CALL REPFIX(IW,GR,REP,IG)
C
      NIRREP = NREX(IG)
      DO I=1,NIRREP
         LBTEMP(1:3) = REP(I)
         LBTEMP(4:4) = ' '
         READ(UNIT=LBTEMP,FMT='(A4)') LABREP(I)
      ENDDO
C
C        4 BYTE MO SYMMETRY LABELS ARE AVAILABLE ON DISK FROM THE SCF
C        USE THESE TO FILL -ISORB- WITH CORRECT INTEGER VALUES
C
      IF(GR.EQ.'C1 ') THEN
         DO I=1,NO+NU
            ISORB(I) = 1
         ENDDO
      ELSE
         CALL DAREAD(IDAF,IODA,LABMO,L1,255,1)
         NERR=0
         DO I=1,NO+NU
            LTEMP = LABMO(I)
            MATCH=0
            DO J=1,NIRREP
               IF(LTEMP.EQ.LABREP(J)) MATCH=J
            ENDDO
            MODI = MATCH
            ! [TODO] Implement conversions to the other
            ! point-group symmetries
            IF (GR.EQ.'D2H') THEN
              IF (MATCH.EQ.2) MODI=5
              IF (MATCH.EQ.3) MODI=2
              IF (MATCH.EQ.4) MODI=6
              IF (MATCH.EQ.5) MODI=3
              IF (MATCH.EQ.6) MODI=7
              IF (MATCH.EQ.7) MODI=4
            END IF
            ISORB(I) = MODI

            IF(MATCH.EQ.0) THEN
               WRITE(IW,9020) I,LTEMP
               NERR=NERR+1
            END IF
         ENDDO

         IF(NERR.GT.0) THEN
            WRITE(IW,9030) NERR
            CALL ABRT
            STOP
         END IF
      END IF

C     Write orbital symmetries to the FCIDUMP header.
      WRITE(297, '(A)', ADVANCE='NO') '  ORBSYM='
      DO IJK=1,NO+NU
          !WRITE(503,933) ISORB(IJK),REP(ISORB(IJK)),ISORB_JAKAL(IJK)
          WRITE(297,934,ADVANCE='NO') ISORB(IJK)
      END DO
      WRITE(297, *)
      CALL FLSHBF(297)

 9020 FORMAT(1X,'ORBITAL',I5,' HAS UNKNOWN SYMMETRY LABEL ',A4)
 9030 FORMAT(1X,'UNABLE TO ASSIGN SYMMETRY TO',I5,' ORBITAL(S).'/
     *       1X,'PLEASE PROVIDE CRISPLY CONVERGED SCF ORBITALS!')
 933  FORMAT(I4,2X,A3,2X,I4)
 934  FORMAT(I0,',')

      END
