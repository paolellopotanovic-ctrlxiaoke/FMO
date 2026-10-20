C 10 Feb 21 - HL  - Add QuanPol
C  6 Jun 18 - YN -  updates for FMO 5.3
C 18 Apr 16 - YN  - updates for FMO 5.2
C 29 Sep 14 - YN  - Implement Density-Functional Tight-Binding Method
C
C*MODULE DFTBGR    *DECK DFTB_GRAD
C>
C>    @brief Gradient of DFTB
C>
C>    @details Calculate gradient of DFTB
C>             If FMO-DFTB, subtract ESP contribution by modifying
C>             Lagrangian.
C>
C>    @author Yoshio Nishimoto
C>    - Sep, 2014- Subroutine written
C>    @date Feb, 2016- Yoshio Nishimoto
C>    - Added time-dependent and PCM
C>    @date Jan, 2018- Van Quan Vuong
C>    - Added LC-DFTB
C>    @date Feb, 2018- Yoshio Nishimoto
C>    - Added FMO-DFTB/AFO and modified derivative calculations
C>
C
      SUBROUTINE DFTB_GRAD
      USE lrcdft, ONLY: EMU
      use mx_limits, only: mxsh,mxgtot,mxatm,mxao
      USE DFTBPB_MOD,ONLY: PERIOD,L2PBC,NKPTS,INDBND,MAXINT,PBCBOX,
     *                     STRESS,DFTB_GAMMA_PBC,DFTB_PBCSTR,
     *                     LATOPT,NDIMPBC,
     * dftb_gamder_pbc,dftb_gdpbc_prep,dftb_disp_pbc
C
      USE comm_EFPFMO
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER (MXSPE=10)
      PARAMETER (ZERO=0.0D+00,HALF=0.5D+00,TWO=2.0D+00,ONE=1.0D+00)
C
      LOGICAL DFTBFL,SCC,SRSCC,DFTB3,DAMPXH,LCDFTB,DC,MFRZ
      LOGICAL TRIPLET,SG1T,TAMMD,TPA,ALPHKWD,BETAKWD
      LOGICAL MREKT,MRDEA
      LOGICAL GOPARR,DSKWRK,MASWRK,mdout
      DOUBLE PRECISION MD
C
      COMMON /DFTB  / DFTBFL,SCC,SRSCC,DFTB3,DAMPXH,LCDFTB
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
      COMMON /EFLDC / EVEC(3),EFLDL
      LOGICAL :: EFLDL
C     COMMON /ENRGYS/ ENUCR,EELCT,ETOT1,SZ,SZZ,ECORE,ESCF,EERD,E1,E2,
C    *                VEN,VEE,EPOT,EKIN,ESTATE(MXRT),STATN,EDFT(3),EDISP
C      COMMON /EFPFMO/ RUNEFP,REFFMO,REPNUCEFP,IEFPFMO,ISWENUC,EFPEFP
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
      COMMON /FMCOM / X(1)
      COMMON /FMODTB/ LGAMMAGLOB,LZREF,LISPE,LSHIFTG,LSHIFTCT,
     *                LGAMMAD,MODESD,MODGAMMA,LDCTMP,LGAMMA3,LGAMMA3D,
     *                LNBOND,LCTMUL,LESPDFTB,LESPASC,NDFTBRST,LDFTBRST
      COMMON /FMOINF/ NFG,NLAYER,NATFMO,NBDFG,NAOTYP,NBODY,NSEGM
      COMMON /FMOOPT/ ESPSCA(9),RESPAP(2),RESPPC(2),RESDIM,RESTRI(4),
     *                RCORSD,RESPCT,CONVFG,CNVDMP,COROFF,RFLMO(4),
     *                ORSHFT,ORSHFT2,CNVAFO,ASCREEN(4),IXESP,MXITFG,
     *                NGUESS,NBSSE,MODORB,MODPAR,IRSTSTP,IRSTLAY,NPRFMO,
     *                NFMOPAL,MODPRP,MAXL1C,IPIEDA,MODGRD,MODESP,IVMUL,
     *                MODLMO,NOPDEN,MOFOCK,MODFD,MODFMM,NCENTM,ndualb,
     *                ngab,modpan
      COMMON /FMOPNT/ LICHFG,LMULFG,LIDMREC,LFRGNAM,LLAYFRG,LINDAT,
     *                LNCBS,LFMOZAN,LFMOC,LFMOMAS,LIZBAS,LIAGLOB,LIBDGH,
     *                LIABDFG,LJABDFG,LNCAO,LIDXCAO,LIAPRJO,LJAPRJO,
     *                LCOREAO,LOCCCOR,LSHIFTB,LIODFMO,LFMODA,LFMODB,
     *                LFMOESPA,LFMOESPB,LLOCFMO,LSCFFRG,LFMOSCF,LRIJ,
     *                LPOPMUL,LPOPMAT,LIALOC,LINDBD,LIATFRG,LINDFRG,
     *                LINDGFRG,LNATFRG,LNAT0FRG,LIANFRG,LZANFRG,LCFRG,
     *                LLIBISH,LLIBNSH,LLIBNG,LINDATG,LFMOBUF(3),LFMODE,
     *                LNUMFRG,LLOCTAT,LIAOGLOB,LLOADM,LFMOGE,LDGRID,
     *                LIODCFMO,LJOB2GRP,LFMOPG,LEMOCDR,LUNTXYZ,LUNTROT,
     *                LSTONEP,LMAPSU,LFRGMUL,LCLMO,LIALMO,LINDLMO,
     *                LATCLMO,LLMOBDF,LFGFLMO,LNFGLMO,LLFGLMO,LPFGLMO,
     *                LPOPDMAT,LIDMPNT,LIDDPNT,LIVMPNT,LIACTFG,lcrfrg,
     *                lzlmfrgv,lYlmfrgv,lndtfrg,lf_mm,lg_mm,lmaxl30,
     *                libuffg,lindatp,lsmon,lsdim,ledimfed,lexcit3d,
     *                lifgfret
      COMMON /FMORUN/ ESPSCF,E0SCF(2),EMP2S,IDAFMO,ICURFG,JCURFG,KCURFG,
     *                ICURLAY,ICURUNT,NAT1E,NCURSH,NGAU,ICURPOP,IFMOSTP,
     *                MONCOR,NEEDR,MODRST,NORBPROJ,NUNESP,ISKIPESP,
     *                IESDPPC,IDOPROP,MP2RUN,ICURIT,IDMFMO,IDDFMO,
     *                IDDCUR,NDDLEFT,IVMFMO,nzmtfmo,ifmobas,itmfmo(2)
      COMMON /FUNCT / ENERGY,EGRAD(3,MXATM)
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /INFOTD/ CNVTOL,PFREQ(2),MODTD,
     *                JANST,NRADT,NTHET,NPHIT,NLEBT,
     *                NSTAT,NTRIAL,MAXVEC,NTHST,IRECTD,ITDFG,ITDPRP,
     *                TRIPLET,SG1T,NONEQR,TAMMD,TPA,ALPHKWD,BETAKWD,
     *                SPCP(3),MULTD,MREKT,MRDEA,MTHST,IFEDAT(4)
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON /MACHSW/ KDIAGG,ICORFL,IXDR,MODIO,mem10,lpnt10,mem10m
      COMMON /MFRPAR/ MFRZ,NUMFRZ,NORFRZ,IFRZ(MXAO)
      COMMON /NSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),
     *                CF(MXGTOT),CG(MXGTOT),CH(MXGTOT),CI(MXGTOT),
     *                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),
     *                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
      COMMON /OUTPUT/ NPRINT,ITOL,ICUT,NORMF,NORMP,NOPK
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /PCMPAR/ IPCM,NFT26,NFT27,IKREP,IEF,IP_F,NFMOPCM,IHET
C     COMMON /RESTAR/ TIMLIM,IREST,NREC,INTLOC,IST,JST,KST,LST
C     COMMON /RUNLAB/ TITLE(10),ANAM(MXATM),BNAM(MXATM),BFLAB(MXAO)
      COMMON /RUNOPT/ RUNTYP,EXETYP,NEVALS,NGLEVL,NHLEVL
      COMMON /WFNOPT/ SCFTYP,VBTYP,DFTYPE,TDDFTYP,CITYP,CCTYP,
     *                MPLEVL,MPCTYP
C
      DOUBLE PRECISION :: SHIFT(NAT),DQ(NAT)
      DOUBLE PRECISION :: AUM(9,9),BUM(9,9),VEC(3)  !! VECTRA(3)
C
      DATA CHECK/8HCHECK   /,DEBUG/8HDEBUG   /,DBUGME/8HRHFCL   /,
     *     RHF  /8HRHF     /,ROHF /8HROHF    /
      DATA UHF  /8HUHF     /,MD   /8HMD      /,RNONE /8HNONE    /
C
      INTEGER :: O,MN0(3),IJ0(3)
      CHARACTER(3) :: CHARA
      LOGICAL :: SOME,DFTD,SKP
c     DOUBLE PRECISION :: MAXDD, SUMS
      dimension gam2tmp(3),gam3itmp(3),gam3jtmp(3),
     *          gam2st(9),gam3ist(9),gam3jst(9)
      logical dampxh1
      data h/8hH       /
C
      SOME = .FALSE.
      IF(MASWRK) THEN
        SOME   = NPRINT.NE.-5 .OR. EXETYP.EQ.DEBUG .OR. EXETYP.EQ.DBUGME
      END IF
      mdout=nfg.eq.0.or.(runefp.ne.md.or.iand(modio,16).eq.0)
     *      .and.iand(nprfmo,3).lt.2
      if (iand(nprfmo,3).eq.3) mdout=.false.
C
      IF (MASWRK.AND.SOME) THEN
        WRITE (IW,*)
        WRITE (IW,'(5X,"-------------------------",5X)')
        WRITE (IW,'(5X,"FIRST DERIVATIVE FOR DFTB",5X)')
        WRITE (IW,'(5X,"-------------------------",5X)')
        WRITE (IW,*)
      END IF
C
C     ----- START THE CLOCKS -----
C
      TIM0 = ZERO
C     TIM1 = ZERO
      CALL TSECND(TIM0)
C     TLEFTS = TIMLIM - TIM0
C
      DFTD = .FALSE.
      IF (IDFTBD.NE.0) THEN
        DC = .TRUE.
        IF (IDFTBD.EQ.3) DFTD = .TRUE.
      ELSE
        DC = .FALSE.
      END IF
C
      L1 = NUM
      L2 = (L1*L1+L1)/2
      L3 = L1*L1
      NORBPRJ = NORFRZ - NUMFRZ
      IF (NORBPRJ.LE.0) NORBPRJ=0
      L2PBC = L2
      IF (PERIOD.AND.NKPTS.GT.0) L2PBC = IXFTCH(X(LINDBND),3)
      ISCAL = 1
      IF(PERIOD.AND.NKPTS.GT.0) ISCAL=2
      NSPIN = 1
      IF(SCFTYP.EQ.UHF.OR.SCFTYP.EQ.ROHF) NSPIN = 2
      IF (SRSCC.OR.SCFTYP.EQ.UHF) THEN
        NDIMGAM = NSHELL
      ELSE
        NDIMGAM = NAT
      END IF
      IF (SRSCC) THEN
        NNSPE = 0
        DO I = 1, NSPE
          NNSPE = NNSPE + MAXANG(I)
        END DO
      ELSE
        NNSPE = NSPE
      END IF
C
C     ----- ALLOCATE MEMORY -----
C
      CALL VALFM(LOADFM)
      CALL GOTFM(NGOTMX)
      LSCR    = LOADFM  + 1
      LIWRK   = LSCR    + L1*8
      LWRK    = LIWRK   + MAX(L1,NAT+1)
      LWRK1   = LWRK    + L1
      LVEC    = LWRK1   + L3*ISCAL
      LD      = LVEC    + L3*(1+NKPTS)*NSPIN*ISCAL
      LD2     = LD      + L2PBC*NSPIN
      LEIG    = LD2     + L2PBC*NSPIN
      LOCC    = LEIG    + L1*(1+NKPTS)*NSPIN
      LMULA   = LOCC    + L1*(1+NKPTS)*NSPIN
      LMULS   = LMULA   + NAT
      LGAMMA  = LMULS   + NSHELL
      LDIST   = LGAMMA  + NDIMGAM*(NDIMGAM+1)/2
      LAST    = LDIST   + (NAT*NAT+NAT)/2
      IF (DFTB3) THEN
        LGAM3 = LAST
        LAST  = LGAM3 + NDIMGAM*NDIMGAM
      ELSE
        LGAM3 = LAST
      END IF
      IF (SCFTYP.EQ.UHF) THEN
        LVA   = LVEC
        LDA   = LD
        LDA2  = LD2
        LEA   = LEIG
        LAOC  = LOCC
C
        LDB   = LDA     + L2PBC
        LDB2  = LDA2    + L2PBC
        LVB   = LVA     + L3*(1+NKPTS)*ISCAL
        LEB   = LEA     + L1*(1+NKPTS)
        LBOC  = LAOC    + L1*(1+NKPTS)
      END IF
      IF (SRSCC.OR.SCFTYP.EQ.UHF) THEN
        LSHIFTS  = LAST
        LSHIFTSH = LSHIFTS + NSHELL
        LISHELL  = LSHIFTSH+ NSHELL
        LAST     = LISHELL + NAT
      END IF
      IF (SCC .AND. PERIOD) THEN
        LTMPGAM  = LAST
        LAST     = LTMPGAM + 11*NNSPE*NNSPE

        ltmpd    = last
        ltmp3d   = ltmpd   + 10*nspe*nspe
        last     = ltmp3d  + 11*nspe*nspe
      END IF
      IF (MFRZ.AND.JCURFG.NE.0) THEN
        TDDFTYP = CHECK
        TRIPLET = .TRUE.
      END IF
      IF (TDDFTYP.NE.RNONE.OR.MFRZ) THEN
        IF (NFG.NE.0.AND..NOT.MFRZ) THEN
          WRITE (IW,*) "NOT YET CONSIDERED"
          CALL ABRT
        END IF
        LDQES    = LAST
        LDQXY    = LDQES   + NAT
        LSHIFTES = LDQXY   + NAT
        LSHIFTXY = LSHIFTES+ NDIMGAM
        LDGS     = LSHIFTXY+ NDIMGAM
        LXYAO    = LDGS    + L2
        LAST     = LXYAO   + L2
        IF (SRSCC) THEN
          LDQESS = LAST
          LDQXYS = LDQESS  + NSHELL
          LAST   = LDQXYS  + NSHELL
        END IF
      END IF
C
      IF (LCDFTB) THEN
        LGAMLC   = LAST
        LGAMDERLC= LGAMLC   + NAT*(NAT+1)/2
        LS       = LGAMDERLC+ NAT*(NAT+1)/2
        LDDS     = LS       + L2
        LAST     = LDDS     + L3
      END IF
C
      NEED    = LAST - LOADFM - 1
C
C     ----- GET MEMORY -----
C
      CALL GETFM(NEED)
C
      IF (EXETYP.EQ.CHECK) GO TO 100
C
C     ----- INITIALIZE GRADIENT -----
C
      CALL VCLR(EGRAD,1,3*NAT)
C
C     ----- RESTORE SEVERAL QUANTITIES -----
C
      !! EIGENVECTORS
      CALL DAREAD(IDAF,IODA,X(LVEC),L3*(1+NKPTS)*ISCAL,15,0)
      !! DENSITY MATRIX
      CALL DAREAD(IDAF,IODA,X(LD),L2PBC,16,0)
      IF (SCFTYP.EQ.ROHF) THEN
        CALL DAREAD(IDAF,IODA,X(LWRK1),L2PBC,20,0)
        CALL VADD(X(LD),1,X(LWRK1),1,X(LD),1,L2PBC)
      END IF
      !! EIGENVALUES
      CALL DAREAD(IDAF,IODA,X(LEIG),L1*(1+NKPTS),17,0)
      !! OCCUPATION NUMBER (ALPHA)
      CALL DAREAD(IDAF,IODA,X(LOCC),L1*(1+NKPTS),562,0)
      IF (SCFTYP.EQ.UHF) THEN
        CALL DAREAD(IDAF,IODA,X(LVB),L3*(1+NKPTS)*ISCAL,19,0)
        CALL DAREAD(IDAF,IODA,X(LDB),L2PBC,20,0)
        CALL DAREAD(IDAF,IODA,X(LEB),L1*(1+NKPTS),21,0)
        !! SPIN SHIFT MATRIX
        IF (SCC) CALL DAREAD(IDAF,IODA,X(LSHIFTS),NSHELL,561,0)
        !! OCCUPATION NUMBER (BETA)
        CALL DAREAD(IDAF,IODA,X(LBOC),L1*(1+NKPTS),563,0)
      END IF
      IF (SCC) THEN
        !! ATOMIC MULLIKEN POPULATION
        CALL DAREAD(IDAF,IODA,X(LMULA),NAT,556,0)
        !! SHELL-RESOLVED MULLIKEN CHARGE (NOT POPULATION)
        CALL DAREAD(IDAF,IODA,X(LMULS),NSHELL,557,0)
        !! SHIFT MATRIX
        CALL DAREAD(IDAF,IODA,SHIFT,NAT,559,0)
        DO I = 1, NAT
          DQ(I) = X(LMULA+I-1) - ZREF(I) !! \Delta q_A
        END DO
        IF (SRSCC) THEN
          CALL DAREAD(IDAF,IODA,X(LSHIFTSH),NSHELL,560,0)
          IF (SCFTYP.NE.UHF) CALL VCLR(X(LSHIFTS),1,NSHELL)
          ISH = 0
          DO I = 1, NAT
            DO J = 1, MAXANG(ISPE(I))
              ISH = ISH + 1
              X(LMULS+ISH-1) = X(LMULS+ISH-1) - QREFL(J,ISPE(I))
            END DO
          END DO
        END IF
      END IF
C
      IF (SRSCC.OR.SCFTYP.EQ.UHF) THEN
        CALL IXSTOR(X(LISHELL),1,0)
        DO IAT = 2, NAT
          CALL IXSTOR(X(LISHELL),IAT,
     *                IXFTCH(X(LISHELL),IAT-1)+MAXANG(ISPE(IAT-1)))
        END DO
      END IF
C
C     ----- CALCULATE DERIVATIVE OF EREP -----
C
      IF (PERIOD) CALL VCLR(STRESS,1,9)
      CALL VCLR(X(LDIST),1,(NAT*NAT+NAT)/2)
      CALL DFTB_REPGRAD(NAT,C,ISPE,NSPE,NUMREP,REPCUT,X(LREPINTV(1,1)),
     *       X(LREPSHORT(1,1)),X(LREPCOEFF(1,1)),EGRAD,X(LDIST),ZREF,
     *       QREF)
C
C     ----- CALCULATE EIGENVALUE WEIGHTED DENSITY MATRIX (X(LD2)) -----
C     SUBROUTINE FROM DMTX
C
      IF (MFRZ) THEN
        CALL DAREAD(IDAF,IODA,X(LD2),L2,36,0)
        CALL DSCAL(L2,-1.0D+00,X(LD2),1)
      ELSE
C       IF (NKPTS.GT.0) THEN
C         CALL DFTB_DMTEIG2(X(LD2),X(LVEC),X(LOCC),X(LEIG),X(LWRK1),L1,
C    *                      L2PBC,NKPTS,NAT,C,PBCBOX,X(LINDBND),MAXINT,
C    *                      X(LVKPT),X(LVKWGHT),IND,NSPIN)
C       ELSE
          CALL DFTB_DMTEIG(X(LD2),X(LVEC),X(LOCC),X(LEIG),NA,L1,L1)
          IF (SCFTYP.EQ.UHF)
     *      CALL DFTB_DMTEIG(X(LDB2),X(LVB),X(LBOC),X(LEB),NB,L1,L1)
C       END IF
      END IF
C     if(nfg.ne.0 .and. scc) call fmolag(x(ld2),x(lwrk1),x(lscr),L1,l2)
C
C     ----- TIME-DEPENDENT DFTB -----
C
      IF(TDDFTYP.NE.RNONE) THEN
        IF (SCFTYP.EQ.UHF) THEN
          WRITE (*,*) "NOT YET"
          CALL ABRT
        END IF
        CALL DAREAD(IDAF,IODA,X(LDGS),L2,308,0) !! GS DENSITY MATRIX
        IF (MFRZ) THEN
          CALL DAREAD(IDAF,IODA,X(LWRK1),L3,476,0)
          CALL SYMSQT(L1,X(LWRK1),L1,X(LXYAO),1)
          CALL DAXPY(L2,ONE,X(LXYAO),1,X(LD),1)
          IF (NFMOPCM.NE.0) CALL DAWRIT(IDAF,IODA,X(LD),L2,16,0)
        END IF
        CALL DAREAD(IDAF,IODA,X(LWRK1),L2,12,0) !! OVERLAP MATRIX
C
C       ----- MAKE q^ES AND q^XY (destroy LVEC) -----
C
        CALL DFTB_TD_GRAD_SHIFT(X(LVEC),X(LWRK1),DAMPXHE,X(LDGS),DQ,
     *                          X(LMULS),X(LDIST),HUBBL,SPE,SPNCST,
     *                          X(LXYAO),X(LDQES),X(LDQESS),X(LDQXY),
     *                          X(LDQXYS),X(LSHIFTES),X(LSHIFTXY),
     *                          X(LGAMMA),X(LGAM3),
     *                          IPCM,NSPE,NAT,L1,L2,IDAF,IODA,ISPE,
     *                          MAXANG,IND,NDIMGAM,NSHELL,
     *                          DFTB3,DAMPXH,SRSCC,TRIPLET)
C
C       ----- UPDATE ENERGY WEIGHTED DENSITY MATRIX -----
C
        IF (MFRZ) THEN
          CALL DAREAD(IDAF,IODA,X(LWRK1),L3,477,0)
          CALL SYMSQT(L1,X(LWRK1),L1,X(LXYAO),1)
          CALL DAXPY(L2,ONE,X(LXYAO),1,X(LD2),1)
          CALL VCLR(X(LXYAO),1,L2)
          CALL DSCAL(NAT,HALF,X(LDQES),1)
          CALL DSCAL(NAT,HALF,X(LSHIFTES),1)
C
          IF (NFG.NE.0.AND.DFTB3) THEN
            !! add ESP contribution
            CALL DAREAD(IDAF,IODA,X(LWRK1),NAT*2,567,0)
            DO NI = 1, NAT
              X(LSHIFTES+NI-1) = X(LSHIFTES+NI-1)
     *          + TWO*X(LWRK1+NI-1)*X(LDQES+NI-1)
            END DO
          END IF
        ELSE
          CALL DSCAL(L2,-1.0D+00,X(LD2),1)
          CALL TDDENS(X(LD2),X(LVEC),L1)
          CALL DSCAL(L2,-1.0D+00,X(LD2),1)
        END IF
      END IF
C
      IF (MFRZ.AND.(NFG.EQ.0.OR.
     *    (IAND(MODGRD,32).EQ.0.AND.JCURFG.EQ.0))) THEN
        TDDFTYP = CHECK
        TRIPLET=.TRUE.
        CALL DCOPY(L2,X(LD),1,X(LDGS),1)
        CALL DAREAD(IDAF,IODA,X(LXYAO),L2,12,0)
        CALL DAREAD(IDAF,IODA,X(LVEC),L3,15,0)
        !! WORKING STORAGE: X(LOCC),X(LWRK1),X(LXYAO),X(LVEC)
        CALL DFTB_MFRZ(X(LD),X(LD2),X(LDQES),X(LSHIFTES),
     *                 X(LOCC),X(LWRK1),X(LEIG),X(LVEC),X(LDIST),
     *                 X(LGAMMA),X(LGAM3),X(LXYAO),DQ)
        CALL VCLR(X(LSHIFTXY),1,NAT)
        CALL VCLR(X(LXYAO),1,L2)
      END IF
C
C     --- ADD ELECTRIC FIELD ---
C
      IF (EFLDL) THEN !! USE EVEC(1:3)
        DO IEVEC = 1, 3
          IF (EVEC(IEVEC).EQ.ZERO) CYCLE
          TEVEC = EVEC(IEVEC)
          DO I = 1, NAT
            EGRAD(IEVEC,I) = EGRAD(IEVEC,I) + TEVEC*DQ(I)
            IF (TDDFTYP.NE.RNONE) THEN
              EGRAD(IEVEC,I) = EGRAD(IEVEC,I) + TEVEC*X(LDQES+I-1)
            END IF
          END DO
        END DO
      END IF
C
C     --- FOR LC-DFTB ---
C
      IF (LCDFTB) THEN
        SCAL = ONE/DBLE(NSPIN)
        IF (TDDFTYP.EQ.RNONE) LDGS = LD
        CALL VCLR(X(LWRK1),1,L2)
        CALL DFTB_MAKEREFD(NAT,X(LWRK1))
        CALL DCOPY(L2,X(LDGS),1,X(LVEC),1)
        CALL DAXPY(L2,-SCAL,X(LWRK1),1,X(LVEC),1)
        CALL CPYTSQ(X(LVEC),X(LDDS),L1,1)
C
        CALL DFTB_LCGAMMADER(EMU,HUBB,NAT,ISPE,X(LDIST),
     *                       X(LGAMLC),X(LGAMDERLC))
        CALL DAREAD(IDAF,IODA,X(LS),L2,12,0)
        IF (TDDFTYP.EQ.RNONE) THEN
C
C         --- CONSTRUCT \widetilde{F}_mn (\Delta D, \Delta D) ---
C
          CALL DFTB_LCGRAD1(0,NAT,L1,ISPE,MAXANG,IND,X(LDDS),X(LDDS),
     *                      X(LS),X(LGAMLC),X(LWRK1))
          CALL DAXPY(L2,-DBLE(NSPIN),X(LWRK1),1,X(LD2),1)
C
C         --- ADD GAMMA DERIVATIVE CONTRIBUTIONS ---
C
          CALL CPYTSQ(X(LS),X(LWRK1),L1,1)
          CALL DFTB_LCGRAD2(0,DBLE(NSPIN),NAT,L1,ISPE,MAXANG,IND,
     *                      X(LDDS),X(LDDS),X(LWRK1),X(LWRK1),X(LDIST),
     *                      X(LGAMDERLC),C,EGRAD,X(LAST))
          IF (SCFTYP.EQ.UHF) THEN
            CALL VCLR(X(LWRK1),1,L2)
            CALL DFTB_MAKEREFD(NAT,X(LWRK1))
            CALL DCOPY(L2,X(LDB),1,X(LVEC),1)
            CALL DAXPY(L2,-SCAL,X(LWRK1),1,X(LVEC),1)
            CALL CPYTSQ(X(LVEC),X(LDDS),L1,1)
            !
            CALL DFTB_LCGRAD1(0,NAT,L1,ISPE,MAXANG,IND,X(LDDS),X(LDDS),
     *                        X(LS),X(LGAMLC),X(LWRK1))
            CALL DAXPY(L2,-DBLE(NSPIN),X(LWRK1),1,X(LDB2),1)
            !
            CALL CPYTSQ(X(LS),X(LWRK1),L1,1)
            CALL DFTB_LCGRAD2(0,DBLE(NSPIN),NAT,L1,ISPE,MAXANG,IND,
     *                        X(LDDS),X(LDDS),X(LWRK1),X(LWRK1),
     *                        X(LDIST),X(LGAMDERLC),C,EGRAD,X(LAST))
          END IF
        ELSE IF (TDDFTYP.NE.RNONE) THEN
          !! Add dF/da*P*D terms
          !! First, S derivative contributions
          CALL DAREAD(IDAF,IODA,X(LVEC),L3,476,0) !! P in AO
          CALL DSCAL(L3,2.0D+00,X(LVEC),1)
          CALL SYMTRZE(X(LVEC),L1,L1)
          CALL DAXPY(L3,0.5D+00,X(LDDS),1,X(LVEC),1)
          !! X(LVEC) = P + (D-D0)/2
          !! X(LDDS) = (D-D0)
          !! Here computes (only exchange)
          !! R^a = 2*(P+(D-D0)/2) * dF(D-D0)/da
          !!     = 2*P * dH(D-D0)/da + (D-D0) * dF(D-D0)/da
          !! The former comes from the solution of Z-vector
          !! The latter comes from the ground state contributions
          CALL DFTB_LCGRAD1(1,NAT,L1,ISPE,MAXANG,IND,X(LVEC),X(LDDS),
     *                      X(LS),X(LGAMLC),X(LWRK1))
          CALL DAXPY(L2,-2.0d+00,X(LWRK1),1,X(LD2),1)
          !! Then, gamma derivative contributions
          CALL CPYTSQ(X(LS),X(LWRK1),L1,1)
          CALL DFTB_LCGRAD2(1,TWO,NAT,L1,ISPE,MAXANG,IND,X(LDDS),
     *                      X(LVEC),X(LWRK1),X(LWRK1),X(LDIST),
     *                      X(LGAMDERLC),C,EGRAD,X(LAST))
C
          !! d(A+B)/da*(X+Y)*(X+Y) (overlap derivative)
          CALL DAREAD(IDAF,IODA,X(LDDS),L3,474,0) !! (X+Y) in AO
          CALL DSCAL(L3,2.0D+00,X(LDDS),1)
          CALL SYMTRZE(X(LDDS),L1,L1)
          CALL DFTB_LCGRAD1(0,NAT,L1,ISPE,MAXANG,IND,X(LDDS),X(LDDS),
     *                      X(LS),X(LGAMLC),X(LWRK1))
          CALL DAXPY(L2,-2.0D+00,X(LWRK1),1,X(LD2),1)
          !! Gamma derivative contributions
          CALL CPYTSQ(X(LS),X(LWRK1),L1,1)
          CALL DFTB_LCGRAD2(0,TWO,NAT,L1,ISPE,MAXANG,IND,X(LDDS),
     *                      X(LDDS),X(LWRK1),X(LWRK1),X(LDIST),
     *                      X(LGAMDERLC),C,EGRAD,X(LAST))
C
          !! d(A+B)/da*(X-Y)*(X-Y) (overlap derivative)
          CALL DAREAD(IDAF,IODA,X(LDDS),L3,475,0) !! (X-Y) in AO
          CALL DSCAL(L3,2.0D+00,X(LDDS),1)
          CALL ASYMTRZE(X(LDDS),L1,L1,TMP)
          CALL DCOPY(L3,X(LDDS),1,X(LVEC),1)
          CALL TRPOSQ(X(LVEC),L1)
          CALL DFTB_LCGRAD1(1,NAT,L1,ISPE,MAXANG,IND,X(LDDS),X(LVEC),
     *                      X(LS),X(LGAMLC),X(LWRK1))
          CALL DAXPY(L2,+2.0D+00,X(LWRK1),1,X(LD2),1)
          !! Gamma derivative contributions
          CALL CPYTSQ(X(LS),X(LWRK1),L1,1)
          CALL DFTB_LCGRAD2(1,-TWO,NAT,L1,ISPE,MAXANG,IND,X(LDDS),
     *                      X(LVEC),X(LWRK1),X(LWRK1),X(LDIST),
     *                      X(LGAMDERLC),C,EGRAD,X(LAST))
        END IF
      END IF
C
C     --- CALCULATE DH/DR AND DS/DR ANALYTICALLY ---
C
      NNN = 4210
      SHIFTIJ   = ZERO
      SHIFTESIJ = ZERO
      SHIFTXYIJ = ZERO
      CALL VCLR(AUM,1,81)
      CALL VCLR(BUM,1,81)
      IF ((SCFTYP.EQ.RHF.OR.SCFTYP.EQ.ROHF).AND..NOT.SRSCC) THEN
      DO I = 1, NAT
        SHIFTI = SHIFT(I)
        IF(TDDFTYP.NE.RNONE) THEN
          SHIFTESI = X(LSHIFTES+I-1)
          SHIFTXYI = X(LSHIFTXY+I-1)
        END IF
        IF (MAXINT.EQ.0) THEN
          NINTPAIR = NAT
        ELSE
          NINTPAIR = INDBND(1,1,I)
        END IF
        DO JJ = 1, NINTPAIR
          IF (MAXINT.EQ.0) THEN
            J = JJ
          ELSE
            J = INDBND(1,JJ+1,I)
          END IF
          IF (PERIOD) THEN
            NNN = INDBND(2,JJ+1,I)
            IXTRV =  NNN/400-10
            IYTRV = (NNN-(IXTRV+10)*400)/20-10
            IZTRV = (NNN-(IXTRV+10)*400-(IYTRV+10)*20)-10
          END IF
          IF (I.EQ.J .AND. NNN.EQ.4210) CYCLE
          FACTORH = 1.0D+00
          FACTORS = 1.0D+00
          IF(NFFAT.GT.0.AND.NLKQMM.GT.0) THEN
            DO ILK=1,NLKQMM
              LK1 = INT(X(LFFLKQMMM+2*(ILK-1)  )+1.0D-03)
              LK2 = INT(X(LFFLKQMMM+2*(ILK-1)+1)+1.0D-03)
              IF((LK1.EQ.I.AND.LK2.EQ.J).OR.
     *           (LK1.EQ.J.AND.LK2.EQ.I)    ) THEN
                 FACTORH = 1.8D+00
                 FACTORS = 1.9D+00
              END IF
            ENDDO
          END IF
          SHIFTJ = SHIFT(J)
          IF (SCC) SHIFTIJ = HALF*(SHIFTI+SHIFTJ)
          IF(TDDFTYP.NE.RNONE) THEN
            SHIFTESJ = X(LSHIFTES+J-1)
            SHIFTXYJ = X(LSHIFTXY+J-1)
            SHIFTESIJ = HALF*(SHIFTESI+SHIFTESJ)
            !! SCALE BY TWO (BUT NOT ONE ACCORDING TO FORMUALTION)
            !! BECAUSE OF SPIN (X(LXYAO) IS NOT DOUBLED)
            SHIFTXYIJ = TWO*(SHIFTXYI+SHIFTXYJ)
          END IF
          VEC(1) = C(1,J) - C(1,I)
          VEC(2) = C(2,J) - C(2,I)
          VEC(3) = C(3,J) - C(3,I)
          IF (PERIOD) THEN
            VEC(1) = VEC(1) + PBCBOX(1,1)*IXTRV
     *             + PBCBOX(2,1)*IYTRV + PBCBOX(3,1)*IZTRV
            VEC(2) = VEC(2) + PBCBOX(1,2)*IXTRV
     *             + PBCBOX(2,2)*IYTRV + PBCBOX(3,2)*IZTRV
            VEC(3) = VEC(3) + PBCBOX(1,3)*IXTRV
     *             + PBCBOX(2,3)*IYTRV + PBCBOX(3,3)*IZTRV
          END IF
          DO K = 1, 3
            CALL DFTB_SKODEA(.FALSE.,ISPE(I),ISPE(J),1,K,0,0,0,
     *        MAXANG(ISPE(I)),MAXANG(ISPE(J)),VEC,
     *        X(LSKGRID(ISPE(I),ISPE(J))),X(LSKGRID(ISPE(J),ISPE(I))),
     *        SKDIM(ISPE(I),ISPE(J)),SKDIM(ISPE(J),ISPE(I)),
     *        X(LSKHTAB(ISPE(I),ISPE(J))),X(LSKHTAB(ISPE(J),ISPE(I))),
     *        X(LSKSTAB(ISPE(I),ISPE(J))),X(LSKSTAB(ISPE(J),ISPE(I))),
     *        X(LSKSELF(ISPE(I))),SKCUT2,AUM,BUM,
     *        .TRUE.,.TRUE.,SKP)
            IF (SKP) EXIT !! CYCLE
            IF (NKPTS.GT.0) NSEQPBC = L2 + INDBND(3,JJ+1,I)
            STRTMP = 0.0D+00
            DO L = 1, IND(J+1) - IND(J)
              O = IND(J)+L
              DO M = 1, IND(I+1) - IND(I)
                N = IND(I)+M
                CALL DFTB_CNVSQ(N,O,NSEQ)
C
C               ----- CALCULATE 2*DH/DR AND 2*(epsilon)*DS/DR -----
C
                DHDR = AUM(M,L)*TWO*FACTORH
                DSDR = BUM(M,L)*TWO*FACTORS
C
C               ----- ADD SCC CONTRIBUTIONS -----
C
                IF (NKPTS.GT.0) THEN
                  NSEQPBC = NSEQPBC + 1
                  TMP = (DHDR+DSDR*SHIFTIJ)*X(LD+NSEQPBC-1)
     *                - DSDR*X(LD2+NSEQPBC-1)
                ELSE IF (SCC) THEN
                  TMP = (DHDR+DSDR*SHIFTIJ)*X(LD+NSEQ-1)
     *                - DSDR*X(LD2+NSEQ-1)
                  IF (TDDFTYP.NE.RNONE) THEN
                    TMP = TMP + DSDR*(SHIFTESIJ*X(LDGS+NSEQ-1)
     *                              + SHIFTXYIJ*X(LXYAO+NSEQ-1))
                  END IF
                ELSE
                  TMP = DHDR*X(LD+NSEQ-1) - DSDR*X(LD2+NSEQ-1)
                END IF
                STRTMP = STRTMP + TMP
              END DO
            END DO
            EGRAD(K,I) = EGRAD(K,I) + STRTMP
            IF (PERIOD) THEN
              DO K2 = 1, 3
                STRESS(K,K2) = STRESS(K,K2) - STRTMP*VEC(K2)
              END DO
            END IF
          END DO
        END DO
      END DO
C
      ELSE IF (SCFTYP.EQ.UHF.OR.SRSCC) THEN
        IJ0(1) = 1
        IJ0(2) = 3
        IJ0(3) = 5
        MN0(1) = 0
        MN0(2) = 1
        MN0(3) = 4
        IF (.NOT.SRSCC) CALL DFTB_SHIFT_ATOSH(SHIFT,X(LSHIFTSH),NAT,
     *    NSHELL,NSPE,ISPE,MAXANG)
        ISH0 = 0
        DO I = 1, NAT
          ISH0 = IXFTCH(X(LISHELL),I)
          IF (MAXINT.EQ.0) THEN
            NINTPAIR = NAT
          ELSE
            NINTPAIR = INDBND(1,1,I)
          END IF
          DO JJ = 1, NINTPAIR
            IF (MAXINT.EQ.0) THEN
              J = JJ
            ELSE
              J = INDBND(1,JJ+1,I)
            END IF
            IF (PERIOD) THEN
              NNN = INDBND(2,JJ+1,I)
              IXTRV =  NNN/400-10
              IYTRV = (NNN-(IXTRV+10)*400)/20-10
              IZTRV = (NNN-(IXTRV+10)*400-(IYTRV+10)*20)-10
            END IF
            IF (I.EQ.J .AND. NNN.EQ.4210) CYCLE
            FACTORH = 1.0D+00
            FACTORS = 1.0D+00
            IF(NFFAT.GT.0.AND.NLKQMM.GT.0) THEN
            DO ILK=1,NLKQMM
              LK1 = INT(X(LFFLKQMMM+2*(ILK-1)  )+1.0D-03)
              LK2 = INT(X(LFFLKQMMM+2*(ILK-1)+1)+1.0D-03)
              IF((LK1.EQ.I.AND.LK2.EQ.J).OR.
     *           (LK1.EQ.J.AND.LK2.EQ.I)    ) THEN
                 FACTORH = 1.8D+00
                 FACTORS = 1.9D+00
              END IF
            ENDDO
            END IF
            JSH0 = IXFTCH(X(LISHELL),J)
            DO K = 1, 3
              VEC(1) = C(1,J) - C(1,I)
              VEC(2) = C(2,J) - C(2,I)
              VEC(3) = C(3,J) - C(3,I)
              IF (PERIOD) THEN
                VEC(1) = VEC(1) + PBCBOX(1,1)*IXTRV
     *                 + PBCBOX(2,1)*IYTRV + PBCBOX(3,1)*IZTRV
                VEC(2) = VEC(2) + PBCBOX(1,2)*IXTRV
     *                 + PBCBOX(2,2)*IYTRV + PBCBOX(3,2)*IZTRV
                VEC(3) = VEC(3) + PBCBOX(1,3)*IXTRV
     *                 + PBCBOX(2,3)*IYTRV + PBCBOX(3,3)*IZTRV
              END IF
              CALL DFTB_SKODEA(I.EQ.J,ISPE(I),ISPE(J),1,K,0,0,0,
     *          MAXANG(ISPE(I)),MAXANG(ISPE(J)),VEC,
     *          X(LSKGRID(ISPE(I),ISPE(J))),X(LSKGRID(ISPE(J),ISPE(I))),
     *          SKDIM(ISPE(I),ISPE(J)),SKDIM(ISPE(J),ISPE(I)),
     *          X(LSKHTAB(ISPE(I),ISPE(J))),X(LSKHTAB(ISPE(J),ISPE(I))),
     *          X(LSKSTAB(ISPE(I),ISPE(J))),X(LSKSTAB(ISPE(J),ISPE(I))),
     *          X(LSKSELF(ISPE(I))),SKCUT2,AUM,BUM,
     *          .TRUE.,.TRUE.,SKP)
              IF (SKP) CYCLE
              IF (NKPTS.GT.0) THEN
                NSEQPBC
     *            = IXFTCH(X(LINDBND),3*(MAXINT+1)*(I-1)+3*(JJ+1))+L2
              END IF
              NU0 = IND(J)
              DO JANG = 1, MAXANG(ISPE(J))
                JSH = JSH0 + JANG
                SHIFTJ = X(LSHIFTSH+JSH-1)
                SHIFTJS = X(LSHIFTS+JSH-1)
                IF(TDDFTYP.NE.RNONE) THEN
                  SHIFTESJ = X(LSHIFTES+JSH-1)
                  SHIFTXYJ = X(LSHIFTXY+JSH-1)
                END IF
                DO N = 1, IJ0(JANG)
                  NU = NU0 + N
                  NN = MN0(JANG) + N
                  MU0 = IND(I)
                  DO IANG = 1, MAXANG(ISPE(I))
                    ISH = ISH0 + IANG
                    SHIFTI = X(LSHIFTSH+ISH-1)
                    SHIFTIS = X(LSHIFTS+ISH-1)
                    SHIFTIJ = HALF*(SHIFTI+SHIFTJ)
                    SHIFTIJS = HALF*(SHIFTIS+SHIFTJS)
                    IF(TDDFTYP.NE.RNONE) THEN
                      SHIFTESI = X(LSHIFTES+ISH-1)
                      SHIFTXYI = X(LSHIFTXY+ISH-1)
                      SHIFTESIJ = HALF*(SHIFTESI+SHIFTESJ)
                      SHIFTXYIJ = TWO*(SHIFTXYI+SHIFTXYJ)
                    END IF
                    DO M = 1, IJ0(IANG)
                      MU = MU0 + M
                      MM = MN0(IANG) + M
                      CALL DFTB_CNVSQ(MU,NU,NSEQ)
C                     if (k.eq.1) write (*,'(9i3)')
C    *                  iang,jang,ish,jsh,mu,nu,m,n,nseq
C
C                     ----- CALCULATE 2*DH/DR AND 2*(epsilon)*DS/DR -----
C
                      DHDR = AUM(MM,NN)*TWO*FACTORH
                      DSDR = BUM(MM,NN)*TWO*FACTORS
C
C                     ----- ADD SCC CONTRIBUTIONS -----
C
                      IF (SCFTYP.EQ.UHF) THEN
                        DGRA = DSDR*(SHIFTIJ+SHIFTIJS)
                        DGRB = DSDR*(SHIFTIJ-SHIFTIJS)
                        IF (NKPTS.GT.0) THEN
                          NSEQPBC = NSEQPBC + 1
                          TMP = (DHDR+DGRA)*X(LDA+NSEQPBC-1)
     *                        + (DHDR+DGRB)*X(LDB+NSEQPBC-1)
     &                      - DSDR*(X(LDA2+NSEQPBC-1)+X(LDB2+NSEQPBC-1))
                        ELSE
                          TMP = (DHDR+DGRA)*X(LDA+NSEQ-1)
     *                        + (DHDR+DGRB)*X(LDB+NSEQ-1)
     &                        - DSDR*(X(LDA2+NSEQ-1)+X(LDB2+NSEQ-1))
                        END IF
                      ELSE
                        IF (NKPTS.GT.0) THEN
                          NSEQPBC = NSEQPBC + 1
                          TMP = (DHDR+DSDR*SHIFTIJ)*X(LD+NSEQPBC-1)
     *                        - DSDR*X(LD2+NSEQPBC-1)
                        ELSE
                          TMP = (DHDR+DSDR*SHIFTIJ)*X(LD+NSEQ-1)
     *                        - DSDR*X(LD2+NSEQ-1)
                        END IF
                        IF (TDDFTYP.NE.RNONE) THEN
                          TMP = TMP + DSDR*(SHIFTESIJ*X(LDGS+NSEQ-1)
     *                                    + SHIFTXYIJ*X(LXYAO+NSEQ-1))
                        END IF
                      END IF
                      EGRAD(K,I) = EGRAD(K,I) + TMP
                      IF (PERIOD) THEN
                        DO K2 = 1, 3
                          STRESS(K,K2) = STRESS(K,K2) - TMP*VEC(K2)
                        END DO
                      END IF
                    END DO
                    MU0 = MU0 + IJ0(IANG)
                  END DO
                END DO
                NU0 = NU0 + IJ0(JANG)
              END DO
            END DO
          END DO
        END DO
      ELSE
        WRITE (IW,*) "SOMETHING IS STRANGE IN GRADIENT CALCULATION"
        CALL ABRT
      END IF
C
C     ---- ADD SCC CONTRIBUTION (2) -----
C
      IF (SCC) THEN
        IF (PERIOD) THEN
          IF (SRSCC) THEN
            NNN = 0
            DO I = 1, NSPE
              NNN = NNN + MAXANG(I)
            END DO
          ELSE
            NNN = NSPE
          END IF
          IF (TDDFTYP.NE.RNONE .AND. TRIPLET)
     *      CALL VCLR(X(LDQXY),1,NAT) !! NO CONTRIBUTION
          IF (.NOT.SRSCC) CALL DCOPY(NAT,DQ,1,X(LMULS),1)
C     call vclr(egrad,1,3*nat)
C     call vclr(stress,1,9)
          CALL DFTB_GAMMA_PBC(1,0,0,NAT,NSPE,NSHELL,NNN,ISPE,MAXANG,
     *                        C,SPE,HUBBL,DAMPXHE,HUBDER,DAMPXH,SRSCC,
     *                        DFTB3,X(LGAMMA),X(LGAM3),X(LTMPGAM),EGRAD,
     *                        X(LAST),X(LMULS),X(LDQES),X(LDQXY),
     *                        TDDFTYP.NE.RNONE,DUM,X(LAST))
!     WRITE (IW,'(X,"STRESS TENSOR IN AU")')
!     DO I = 1, 3
!       WRITE (IW,'(3F20.10)') (STRESS(I,J),J=1,3)
!     END DO
!     do i = 1, nat
!       write (*,'(i3,3f20.10)') i,(egrad(j,i),j=1,3)
!     end do
!     do i = 1, 3
!       write (*,'(i3,3f20.10)') i,(stress(i,j),j=1,3)
!     end do
!     call dftb_gdpbc_prep(x(ltmpd),x(ltmp3d),nspe,hubbl,dftb3,dampxh)
!     call vclr(egrad,1,3*nat)
!     call vclr(stress,1,9)
!     do i = 1, nat
!       isp = ispe(i)
!       dqi = dq(i)
!       do j = 1, i !! -1
!         jsp = ispe(j)
!         dqj = dq(j)
!         ish = isp-1+nspe*(jsp-1)
!         jsh = jsp-1+nspe*(isp-1)
!         CHARA2 = SPE(ISPE(JAT))
!         DAMPXH1 = DAMPXH.AND.(SPE(ISP).EQ.H.OR.SPE(JSP).EQ.H)
!         CALL DFTB_DMPCHK(ISP,JSP,DAMPXH1)
!         call dftb_gamder_pbc(i,j,c(1,i),c(1,j),gam2tmp,gam3itmp,
!    *                         gam3jtmp,gam2st,gam3ist,gam3jst,
!    *         x(ltmpd+10*ish),x(ltmp3d+11*ish),x(ltmp3d+11*jsh),
!    *         dampxh1,isp,jsp)
!     write (*,*) "i,j = ", i,j
!     write (*,*) gam2tmp
!          call prsq(gam2st,3,3,3)
!         egrad(1,i) = egrad(1,i) + gam2tmp(1)*dqi*dqj
!         egrad(2,i) = egrad(2,i) + gam2tmp(2)*dqi*dqj
!         egrad(3,i) = egrad(3,i) + gam2tmp(3)*dqi*dqj
!         egrad(1,j) = egrad(1,j) - gam2tmp(1)*dqi*dqj
!         egrad(2,j) = egrad(2,j) - gam2tmp(2)*dqi*dqj
!         egrad(3,j) = egrad(3,j) - gam2tmp(3)*dqi*dqj
!         call daxpy(9,dqi*dqj,gam2st,1,stress,1)
!         if (dftb3) then
!     write (*,*) "i,j = ", i,j
!     write (*,*) gam3itmp
!     write (*,*) gam3jtmp
!         egrad(1,i) = egrad(1,i)
!    *      +(gam3itmp(1)*dqi+gam3jtmp(1)*dqj)*dqi*dqj/3.0d+00
!         egrad(2,i) = egrad(2,i)
!    *      +(gam3itmp(2)*dqi+gam3jtmp(2)*dqj)*dqi*dqj/3.0d+00
!         egrad(3,i) = egrad(3,i)
!    *      +(gam3itmp(3)*dqi+gam3jtmp(3)*dqj)*dqi*dqj/3.0d+00
!         egrad(1,j) = egrad(1,j)
!    *      -(gam3itmp(1)*dqi+gam3jtmp(1)*dqj)*dqi*dqj/3.0d+00
!         egrad(2,j) = egrad(2,j)
!    *      -(gam3itmp(2)*dqi+gam3jtmp(2)*dqj)*dqi*dqj/3.0d+00
!         egrad(3,j) = egrad(3,j)
!    *      -(gam3itmp(3)*dqi+gam3jtmp(3)*dqj)*dqi*dqj/3.0d+00
!         call daxpy(9,dqi*dqi*dqj/3.0d+00,gam3ist,1,stress,1)
!         call daxpy(9,dqi*dqj*dqj/3.0d+00,gam3jst,1,stress,1)
!         end if
!       end do
!     end do
!     do i = 1, nat
!       write (*,'(i3,3f20.10)') i,(egrad(j,i),j=1,3)
!     end do
!     do i = 1, 3
!       write (*,'(i3,3f20.10)') i,(stress(i,j),j=1,3)
!     end do
        ELSE IF (TDDFTYP.NE.RNONE) THEN
C
C         ----- CALCULATE K TERMS AT THE SAME TIME -----
C
          CALL DFTB_GAMDER(NAT,NSPE,HUBBL,ISPE,X(LGAMMA),DAMPXH,
     *      DAMPXHE,SPE,X(LDIST),DFTB3,X(LGAM3),HUBDER,SRSCC,MAXANG,
     *      NDIMGAM,NNSPE)
          IF (TRIPLET) THEN !! NO CONTRIBUTION
            CALL VCLR(X(LDQXY),1,NAT)
            IF (SRSCC) CALL VCLR(X(LDQXYS),1,NSHELL)
          END IF
          ISH0 = 0
          DO I = 1, NAT
            JSH0 = 0
            DO J = 1, I-1
              CALL DFTB_CNVSQ(I,J,NSEQ)
              VEC(1) = C(1,J) - C(1,I)
              VEC(2) = C(2,J) - C(2,I)
              VEC(3) = C(3,J) - C(3,I)
              DIST = X(LDIST+NSEQ-1)
              IF (SRSCC) THEN
                VAL = 0.0D+00
                DO II = 1, MAXANG(ISPE(I))
                  ISH = ISH0 + II
                  DQI = X(LMULS+ISH-1)
                  DO JJ = 1, MAXANG(ISPE(J))
                    JSH = JSH0 + JJ
                    DQJ = X(LMULS+JSH-1)
                    CALL DFTB_CNVSQ(ISH,JSH,NSEQ)
                    GAMMAV = X(LGAMMA+NSEQ-1)
                    VAL = VAL + GAMMAV*(DQI*DQJ
     *                  + X(LDQESS+ISH-1)*DQJ+DQI*X(LDQESS+JSH-1)
     *                  + 2.0D+00*X(LDQXYS+ISH-1)*X(LDQXYS+JSH-1))
                    IF (DFTB3) THEN
                      GAMI = X(LGAM3+ISH-1+(JSH-1)*NSHELL)
                      GAMJ = X(LGAM3+JSH-1+(ISH-1)*NSHELL)
                      VAL = VAL
     *                  + DQI*DQJ*(DQ(I)*GAMI+DQ(J)*GAMJ)
     *                  + (X(LDQESS+ISH-1)*DQJ+DQI*X(LDQESS+JSH-1))
     *                     *(GAMI*DQ(I)+GAMJ*DQ(J))
     *                  + (GAMI*X(LDQES+I-1)+GAMJ*X(LDQES+J-1))*DQI*DQJ
     *                  + 2.0D+00*(
     *                    X(LDQXYS+ISH-1)*X(LDQXYS+JSH-1)
     *                    *(GAMI*DQ(I)+GAMJ*DQ(J))
     *                    + (X(LDQXYS+ISH-1)*DQJ+DQI*X(LDQXYS+JSH-1))
     *                      *(GAMI*X(LDQXY+I-1)+GAMJ*X(LDQXY+J-1)))
                    END IF
                  END DO
                END DO
              ELSE
                GAMMAV = X(LGAMMA+NSEQ-1)
                VAL = DQ(I)*DQ(J)*GAMMAV !! GS
     *              + X(LDQES+I-1)*DQ(J)*GAMMAV !! DQES
     *              + DQ(I)*X(LDQES+J-1)*GAMMAV !! DQES
     *              + X(LDQXY+I-1)*X(LDQXY+J-1)*GAMMAV !! DQXY
     *              + X(LDQXY+J-1)*X(LDQXY+I-1)*GAMMAV !! DQXY
                IF (DFTB3) THEN
                  GAMI = X(LGAM3+I-1+(J-1)*NAT) !! G_ac
                  GAMJ = X(LGAM3+J-1+(I-1)*NAT) !! G_ca
                  VAL = VAL
     *              + DQ(I)*DQ(J)*(DQ(I)*GAMI+DQ(J)*GAMJ) !! GS
     *              + GAMI*X(LDQES+I-1)*DQ(I)*DQ(J)*2.0D+00 !! DQES
     *              + GAMJ*X(LDQES+I-1)*DQ(J)*DQ(J)         !! DQES
     *              + GAMJ*X(LDQES+J-1)*DQ(J)*DQ(I)*2.0D+00 !! DQES
     *              + GAMI*X(LDQES+J-1)*DQ(I)*DQ(I)         !! DQES
     *              +(GAMI*DQ(I)*X(LDQXY+I-1)*X(LDQXY+J-1)*2.0D+00
     *              + GAMJ*DQ(J)*X(LDQXY+I-1)*X(LDQXY+J-1)*2.0D+00
     *              + GAMI*DQ(J)*X(LDQXY+I-1)*X(LDQXY+I-1)
     *              + GAMJ*DQ(I)*X(LDQXY+J-1)*X(LDQXY+J-1))*2.0D+00!!DQXY
                END IF
              END IF
              VAL = VAL / DIST
              DO K = 1, 3
                EGRAD(K,I) = EGRAD(K,I) - VEC(K)*VAL
                EGRAD(K,J) = EGRAD(K,J) + VEC(K)*VAL
              END DO
              JSH0 = JSH0 + MAXANG(ISPE(J))
            END DO
            ISH0 = ISH0 + MAXANG(ISPE(I))
          END DO
        ELSE IF (SRSCC) THEN
          NN = 0
          DO I = 1, NSPE
            NN = NN + MAXANG(I)
          END DO
          CALL DFTB_GAMDER(NAT,NSPE,HUBBL,ISPE,X(LGAMMA),DAMPXH,
     *      DAMPXHE,SPE,X(LDIST),DFTB3,X(LGAM3),HUBDER,SRSCC,MAXANG,
     *      NSHELL,NN)
          ISH0 = 0
          DO I = 1, NAT
            JSH0 = 0
            DO J = 1, I-1
              CALL DFTB_CNVSQ(I,J,NSEQ)
              VEC(1) = C(1,J) - C(1,I)
              VEC(2) = C(2,J) - C(2,I)
              VEC(3) = C(3,J) - C(3,I)
              DIST   = X(LDIST+NSEQ-1)
              VAL    = ZERO
              VALTMP1= ZERO
              VALTMP2= ZERO
              DO II = 1, MAXANG(ISPE(I))
                ISH = ISH0 + II
                DQI = X(LMULS+ISH-1)
                DO JJ = 1, MAXANG(ISPE(J))
                  JSH = JSH0 + JJ
                  DQJ = X(LMULS+JSH-1)
                  CALL DFTB_CNVSQ(ISH,JSH,NSEQ)
                  GAMMAV = X(LGAMMA+NSEQ-1)
                  VAL = VAL + DQI*DQJ*GAMMAV
                  IF (DFTB3) THEN
                    GAMI = X(LGAM3+ISH-1+(JSH-1)*NSHELL) !! G_ac
                    GAMJ = X(LGAM3+JSH-1+(ISH-1)*NSHELL) !! G_ca
C                   VALTMP = VALTMP + DQI*DQJ*GAMI + DQJ*DQI*GAMJ
                    VALTMP1 = VALTMP1 + GAMI*DQI*DQJ
                    VALTMP2 = VALTMP2 + GAMJ*DQJ*DQI
                  END IF
                END DO
              END DO
              IF (DFTB3) VAL = VAL + VALTMP1*DQ(I) + VALTMP2*DQ(J)
              VAL = VAL / DIST
              DO K = 1, 3
                EGRAD(K,I) = EGRAD(K,I) - VEC(K)*VAL
                EGRAD(K,J) = EGRAD(K,J) + VEC(K)*VAL
              END DO
              JSH0 = JSH0 + MAXANG(ISPE(J))
            END DO
            ISH0 = ISH0 + MAXANG(ISPE(I))
          END DO
        ELSE
          CALL DFTB_GAMDER(NAT,NSPE,HUBBL,ISPE,X(LGAMMA),DAMPXH,
     *      DAMPXHE,SPE,X(LDIST),DFTB3,X(LGAM3),HUBDER,SRSCC,MAXANG,
     *      NAT,NSPE)
          DO I = 1, NAT
            DO J = 1, I-1
              CALL DFTB_CNVSQ(I,J,NSEQ)
              VEC(1) = C(1,J) - C(1,I)
              VEC(2) = C(2,J) - C(2,I)
              VEC(3) = C(3,J) - C(3,I)
              DIST = X(LDIST+NSEQ-1)
              GAMMAV = X(LGAMMA+NSEQ-1)
              VAL = DQ(I)*DQ(J)*GAMMAV
              IF (DFTB3) THEN
                GAMI = X(LGAM3+I-1+(J-1)*NAT) !! G_ac
                GAMJ = X(LGAM3+J-1+(I-1)*NAT) !! G_ca
                VAL = VAL + DQ(I)*
     *            (DQ(I)*DQ(J)*GAMI + DQ(J)*DQ(J)*GAMJ)
              END IF
              VAL = VAL / DIST
              DO K = 1, 3
                EGRAD(K,I) = EGRAD(K,I) - VEC(K)*VAL
                EGRAD(K,J) = EGRAD(K,J) + VEC(K)*VAL
              END DO
            END DO
          END DO
        END IF
      END IF
C
C     ---- ADD DISPERSION CONTRIBUTION -----
C
      IF (DC) THEN
        IF (DFTD) THEN
          CALL DFTD3(2,DUMMY,EGRAD) !! EXPERIMENTAL
        ELSE
          X(LIWRK+NAT) = ZERO
          IF (NFG.NE.0.AND.IDFTBD.EQ.4) THEN
            CALL IXSTOR(X(LIWRK+NAT),1,NFG)
            CALL DFTB_FMONBOND(NAT,X(LIAGLOB),X(LNBOND),X(LIWRK))
          END IF
          IF (PERIOD) THEN
            CALL DFTB_DISP_PBC(1,EDISP,EGRAD,X(LAST))
          ELSE
            CALL DFTB_DISP_GRAD(QREF,X(LDIST),EGRAD,X(LIWRK))
          END IF
        END IF
      END IF
C
      IF (MFRZ) THEN
        TDDFTYP = RNONE
        TRIPLET = .FALSE.
      END IF
C
      IF (PERIOD.AND.NFG.EQ.0) THEN
        CALL DFTB_PBCSTR(X(LTMPGAM),NDIMPBC,SOME)
        CALL DCOPY(3*NDIMPBC,X(LTMPGAM),1,EGRAD(1,NAT+1),1)
      END IF
C
  100 CONTINUE
c     write(6,*) 'wwwhoi',PERIOD,LATOPT,NAT,NDIMPBC
      IF (PERIOD.AND.LATOPT.NE.0) THEN
        CALL DAWRIT(IDAF,IODA,EGRAD,3*(NAT+NDIMPBC),3,0)
      ELSE
        CALL DAWRIT(IDAF,IODA,EGRAD,3*NAT,3,0)
      END IF
C
      IF (MASWRK.AND.SOME) THEN
        CHARA(1:3) = "   "
        IF (.NOT.SCC) CHARA(1:1) = "1"
        IF (SCC.AND..NOT.DFTB3) CHARA(1:1) = "2"
        IF (SCC.AND.DFTB3) CHARA(1:1) = "3"
        IF (DC) CHARA(2:3) = "-D"
        WRITE(IW,'(" ...... END OF DFTB",A3," GRADIENT ......")') CHARA
      END IF
      if(mdout) CALL TIMIT(1)
C
      CALL RETFM(NEED)
C
      RETURN
C
      END SUBROUTINE DFTB_GRAD
C
C-----------------------------------------------------------------------
C*MODULE DFTBGR    *DECK DFTB_PREP_GAMMA_GRAD
C>
C>    @brief Gradient of DFTB
C>
C>    @details Calculate TMP matrix used for gamma values
C>
C>    @author Yoshio Nishimoto
C>    - Sep, 2014- Subroutine written
C>
C>           --- INPUT ---
C>    @param HUBA,HUBB Hubbard values of atom A and B
C>    @param DFTB3 Whether DFTB3 or not
C>    @param DAMPXH Whether X-H damping or not
C>           --- OUTPUT ---
C>    @param TMP Constants for gamma calculation
C>
C
      SUBROUTINE DFTB_PREP_GAMMA_GRAD(HUBA,HUBB,TMP,DFTB3,DAMPXH)
C
      IMPLICIT NONE
C
      DOUBLE PRECISION, INTENT(IN) :: HUBA,HUBB
      LOGICAL, INTENT(IN) :: DFTB3,DAMPXH
      DOUBLE PRECISION, INTENT(OUT) :: TMP(11)
C
      DOUBLE PRECISION :: TAUA,TAUA2,TAUA3,TAUA4,TAUA6,TAUB,TAUB2,TAUB3,
     *  TAUB4,TAUB6,TAUAB,TAUABI,TAUABI2,TAUABI3,TAUABI4
      DOUBLE PRECISION,PARAMETER :: HUBTOL=1.0D-04,
     *  HALF=0.5D+00,ONE=1.0D+00,TWO=2.0D+00,THREE=3.0D+00,
     *  ONE_24=ONE/2.4D+01,TWELVE=1.2D+01,ONE_15=ONE/1.5D+01
C
      TAUA = 3.2D+00*HUBA
      TAUA2 = TAUA**2
      TAUA3 = TAUA*TAUA2
      TAUA4 = TAUA2**2
      TAUA6 = TAUA2*TAUA4
      TAUB = 3.2D+00*HUBB
      TAUB2 = TAUB**2
      TAUB3 = TAUB*TAUB2
      TAUB4 = TAUB2**2
      TAUB6 = TAUB2*TAUB4
      IF (ABS(HUBA-HUBB).LE.HUBTOL) THEN
        TMP(1) = 1.875D-01*TAUA2 !! 3/16 a^2
        TMP(2) = ONE_24*TAUA3    !! 1/24 a^3
        TMP(3) = 6.875D-01*TAUA  !! 11/16 a
        TMP(4) = HALF*TMP(2) !! 1/48 a^3
        IF (DFTB3) THEN
          TMP(5) = -6.0D-01*TAUA
          TMP(6) = 2.0D-01*TAUA2
          TMP(7) = 2.0D-01*TAUA3
          TMP(8) = TAUA4*ONE_15
          IF (DAMPXH) THEN
            TMP(9)  = -TAUA
            TMP(10) = -4.0D-01*TAUA2
            TMP(11) = -TAUA3*ONE_15
          END IF
        END IF
      ELSE
        TAUAB = TAUA2 - TAUB2
        TAUABI = ONE/TAUAB
        TAUABI2 = TAUABI**2
        TAUABI3 = TAUABI**3
        TMP(1) =  (TAUB6-THREE*TAUA2*TAUB4)*TAUABI3
        TMP(2) = HALF*TAUA*TAUB4*TAUABI2
        TMP(3) = -(TAUA6-THREE*TAUB2*TAUA4)*TAUABI3
        TMP(4) = HALF*TAUB*TAUA4*TAUABI2
        IF (DFTB3) THEN
          TAUABI4 = TAUABI**4
          TMP(5) = -HALF*(TAUB6+THREE*TAUA2*TAUB4)*TAUABI3
          TMP(6) = -TWELVE*TAUA3*TAUB4*TAUABI4
          TMP(7) = -TWO*TAUB3*TAUA3*TAUABI3
          TMP(8) =  TWELVE*TAUB4*TAUA3*TAUABI4
        END IF
      END IF
C
      RETURN
C
      END SUBROUTINE DFTB_PREP_GAMMA_GRAD
C
C-----------------------------------------------------------------------
C*MODULE DFTBGR    *DECK DFTB_GAMDERF
C>
C>    @brief Gradient of DFTB
C>
C>    @details Function which calculates gamma 1st-order derivative
C>
C>    @author Yoshio Nishimoto
C>    - Sep, 2014- Subroutine written
C>
C>           --- INPUT ---
C>    @param HUBA,HUBB Hubbard values of atom A and B
C>    @param HUBDER Hubbard derivative (DFTB3)
C>    @param DIST Distance between atom A and B
C>    @param DAMPXHE Exponent of X-H damping
C>    @param TMP Constants of the pair atom A and B
C>    @param DFTB3 Whether DFTB3 or not
C>    @param HDAMP Whether X-H damping or not, decided before calling
C>           this function
C>           --- OUTPUT ---
C>    @param GAMMA2 Gamma used in DFTB2
C>    @param GAMMA3 Gamma used in DFTB3
C>
C
      SUBROUTINE DFTB_GAMDERF(GAMMA2,GAMMA3,HUBA,HUBB,HUBDER,DIST,
     *  DAMPXHE,TMP,DFTB3,HDAMP)
C
      IMPLICIT NONE
C
      DOUBLE PRECISION, INTENT(IN) :: HUBA,HUBB,DIST,HUBDER,TMP(11),
     *  DAMPXHE
      DOUBLE PRECISION, INTENT(OUT) :: GAMMA2,GAMMA3
      LOGICAL, INTENT(IN) :: DFTB3,HDAMP
C
      DOUBLE PRECISION :: TAUA,TAUB,EXPA,EXPB,DIST2,DISTI,DISTI2,
     * HUBMEAN,HUBMEANEXP,EXPH,GAMTMP,GAMTMP2,GA,DGA_DR,FAB,FBA,
     * DFAB_DR,DFBA_DR,GAMTMP3,DIST3,DFAB_DA,DFBA_DA,D2FAB_DADR,S,
     * D2FBA_DADR,GAMTMP4,DS_DR,DS_DU,DH_DR,DH_DU,D2H_DUDR,HUBMEANEXP1
      DOUBLE PRECISION,PARAMETER :: DISTTOL=1.0D-08,HUBTOL=1.0D-04,
     *  ZERO=0.0D+00,HALF=0.5D+00,ONE=1.0D+00,TWO=2.0D+00,THREE=3.0D+00,
     *  ONE_THIRD=ONE/THREE
C
      TAUA = 3.2D+00*HUBA
      TAUB = 3.2D+00*HUBB
      GAMTMP2 = ZERO
      GAMTMP4 = ZERO
      IF (DIST.LE.DISTTOL) THEN
        GAMMA2 = ZERO
        IF (DFTB3) GAMMA3 = ZERO
      ELSE
        DIST2 = DIST**2
        DISTI = ONE/DIST
        DISTI2 = DISTI**2
        EXPA = EXP(-TAUA*DIST)
        IF (ABS(HUBA-HUBB).LE.HUBTOL) THEN
          !! CALCULATE HUBA EQ HUBB
          GA = DISTI+TMP(3)+TMP(1)*DIST
     *       + TMP(4)*DIST2
          DGA_DR = -DISTI2+TMP(1)+TMP(2)*DIST
          DS_DR = EXPA*DGA_DR - TAUA*EXPA*GA
          GAMTMP = -DS_DR
          IF (HDAMP) THEN
            HUBMEAN = HUBA
            HUBMEANEXP = HUBMEAN**DAMPXHE
            EXPH = EXP(-HUBMEANEXP*DIST2)
            GAMTMP = GAMTMP * EXPH !! -h*dSg/dR
            GAMTMP2 = EXPA*GA
            GAMTMP2 = GAMTMP2*HUBMEANEXP*TWO*DIST*EXPH
          END IF
          IF (DFTB3) THEN !! 3rd
            DIST3 = DIST2*DIST
            GAMTMP3 = EXPA*(TMP(6)*DIST
     *        + TMP(7)*DIST2 + TMP(8)*DIST3)
            GAMTMP3 = -GAMTMP3
            IF (HDAMP) THEN !! damp correction with 3rd
              HUBMEANEXP1 = HUBMEAN**(DAMPXHE-ONE)
              DS_DU = EXPA*(-ONE+TMP(9)*DIST
     *          + TMP(10)*DIST2+TMP(11)*DIST3)
              DH_DR = -TWO*DIST*HUBMEANEXP*EXPH
              DH_DU = -HALF*DIST2*DAMPXHE*HUBMEANEXP1*EXPH
              D2H_DUDR = DAMPXHE*DIST*HUBMEANEXP1
     *                   *(DIST2*HUBMEANEXP-ONE)*EXPH
              S = EXPA * GA
              GAMTMP4 = DS_DU*DH_DR + DS_DR*DH_DU + S*D2H_DUDR
              GAMTMP3 = GAMTMP3*EXPH - GAMTMP4
            END IF
            GAMMA3 = GAMTMP3 * HUBDER * ONE_THIRD
          END IF
        ELSE
          !! CALCULATE HUBA NE HUBB
          EXPB = EXP(-TAUB*DIST)
          FAB = TMP(2)-TMP(1)*DISTI
          FBA = TMP(4)-TMP(3)*DISTI
          DFAB_DR = TMP(1)*DISTI2
          DFBA_DR = TMP(3)*DISTI2
          !! CALCULATE -dSf/dR
          DS_DR = EXPA*DFAB_DR - TAUA*EXPA*FAB
     *          + EXPB*DFBA_DR - TAUB*EXPB*FBA
          GAMTMP = -DS_DR
          IF (HDAMP) THEN
            HUBMEAN = HALF*(HUBA+HUBB)
            HUBMEANEXP = HUBMEAN**DAMPXHE
            EXPH = EXP(-HUBMEANEXP*DIST2)
            GAMTMP = GAMTMP * EXPH !! -h*dSf/dR
            GAMTMP2 = EXPA*FAB + EXPB*FBA
            GAMTMP2 = GAMTMP2*HUBMEANEXP*TWO*DIST*EXPH
          END IF
          IF (DFTB3) THEN !! 3rd
            DFAB_DA = TMP(5) + TMP(6)*DISTI
            DFBA_DA = TMP(7) + TMP(8)*DISTI
            D2FAB_DADR = -TMP(6)*DISTI2
            D2FBA_DADR = -D2FAB_DADR !! CHECK
C           D2FBA_DADR = D2FAB_DADR
            GAMTMP3 = EXPA*((TAUA*DIST-ONE)*FAB-TAUA*DFAB_DA
     *        + D2FAB_DADR - DIST*DFAB_DR)
     *        + EXPB*(D2FBA_DADR - TAUB*DFBA_DA)
            GAMTMP3 = -3.2D+00*GAMTMP3
            IF (HDAMP) THEN  !! damp correction with 3rd
              HUBMEANEXP1 = HUBMEAN**(DAMPXHE-ONE)
              DS_DU = EXPA*(DFAB_DA-DIST*FAB)+EXPB*DFBA_DA
              DS_DU = 3.2D+00*DS_DU
              DH_DR = -TWO*DIST*HUBMEANEXP*EXPH
              DH_DU = -HALF*DIST2*DAMPXHE*HUBMEANEXP1*EXPH
              D2H_DUDR = DAMPXHE*DIST*HUBMEANEXP1
     *                   *(DIST2*HUBMEANEXP-ONE)*EXPH
              S = EXPA*FAB + EXPB*FBA
              GAMTMP4 = DS_DU*DH_DR + DS_DR*DH_DU + S*D2H_DUDR
              GAMTMP3 = GAMTMP3*EXPH - GAMTMP4
            END IF
            GAMMA3 = GAMTMP3 * HUBDER * ONE_THIRD
          END IF
        END IF
        GAMMA2 = -DISTI2 + GAMTMP + GAMTMP2
C       WRITE (*,*) GAMMA(NSEQ),GAM3(I,J)
      END IF
C
      RETURN
C
      END SUBROUTINE DFTB_GAMDERF
C
C-----------------------------------------------------------------------
C*MODULE DFTBGR    *DECK DFTB_GAMDER
C>
C>    @brief Gradient of DFTB
C>
C>    @details Calculates gamma 1st-order derivative
C>             First, calculate constants which do not depend on
C>             inter-atomic distance, but only depends on species of
C>             atom A and B, which is done by DFTB_PREP_GAMMA_GRAD
C>             subroutine.  Then, calculate actual gamma derivatives
C>             with DFTB_GAMDERF function.
C>
C>    @author Yoshio Nishimoto
C>    - Sep, 2014- Subroutine written
C>    @date Feb, 2016- Yoshio Nishimoto
C>    - Bug fix
C>
C>           --- INPUT ---
C>    @param NAT Number of atoms
C>    @param NSPE Number of species
C>    @param HUBBL Shell resolved Hubbard derivatives
C>    @param ISPE Index of species
C>    @param DAMPXH Whether H-X damping or not
C>    @param DAMPXHE Exponent of X-H damping
C>    @param SPE Name of each species
C>    @param DISTMAT Matrix of inter-atomic distances
C>    @param DFTB3 Whether DFTB3 or not
C>    @param HUBDER Hubbard derivative (DFTB3)
C>    @param SRSCC Shell-resolved SCC or not
C>    @param MAXANG Maximum angular momentum of each species
C>    @param NDIMGAM Dimension of gamma matrix, if SRSCC=.T., NDIMGAM=
C>           NSHELL, if not, NDIMGAM=NAT
C>    @param NN Number of species to be calculated.  If SRSCC=.T.,
C>           this is the sum of the number of shell of each atom, and
C>           if not, NN=NSPE
C>           --- OUTPUT ---
C>    @param GAMMA Gamma used in DFTB2
C>    @param GAMMA3 Gamma used in DFTB3
C>
C
      SUBROUTINE DFTB_GAMDER(NAT,NSPE,HUBBL,ISPE,GAMMA,DAMPXH,
     *  DAMPXHE,SPE,DISTMAT,DFTB3,GAMMA3,HUBDER,SRSCC,MAXANG,NDIMGAM,NN)
C
      IMPLICIT NONE
C
      INTEGER,INTENT(IN) :: NAT,NSPE,ISPE(NAT),NDIMGAM,MAXANG(NSPE),NN
      DOUBLE PRECISION,INTENT(IN) :: HUBBL(3,NSPE),DISTMAT(*),DAMPXHE,
     *  HUBDER(*),SPE(*)
      LOGICAL,INTENT(IN) :: DAMPXH,DFTB3,SRSCC
      DOUBLE PRECISION,INTENT(OUT) :: GAMMA(*),GAMMA3(NDIMGAM,NDIMGAM)
C
      INTEGER :: I,II,ISH,J,JJ,JSH,NSP1,NSP2,NSEQ,ISHELL(NDIMGAM),
     *  IND(NSPE),ISH0,JSH0
      DOUBLE PRECISION :: HUBA,HUBB,DIST,HUBDERTMP,CHARA1,CHARA2
      DOUBLE PRECISION :: TMP(11,NN,NN)
      DOUBLE PRECISION H/8HH       /
      LOGICAL :: HDAMP
C
      CALL VCLR(TMP,1,11*NN*NN)
      IF (SRSCC) THEN
        IND(1) = 0
        DO I = 2, NSPE
          IND(I) = IND(I-1) + MAXANG(I-1) !MAXANG(ISPE(I-1))
        END DO
        ISH = 0
        JSH = 0
        DO I = 1, NAT
          JSH = IND(ISPE(I))
          DO J = 1, MAXANG(ISPE(I))
            ISH = ISH + 1
            ISHELL(ISH) = JSH + J
          END DO
        END DO
        ISH = 0
        DO I = 1, NSPE
          DO II = 1, MAXANG(I)
            ISH = ISH + 1
            HUBA = HUBBL(II,I)
            JSH = 0
            DO J = 1, NSPE
              DO JJ = 1, MAXANG(J)
                JSH = JSH + 1
                HUBB = HUBBL(JJ,J)
C               write (*,*) ish,jsh
                CALL DFTB_PREP_GAMMA_GRAD(HUBA,HUBB,TMP(1,ISH,JSH),
     *            DFTB3,DAMPXH)
              END DO
            END DO
          END DO
        END DO
C
        ISH0 = 0
        DO I = 1, NAT
          NSP1 = ISPE(I)
          CHARA1 = SPE(ISPE(I))
          JSH0 = 0
          DO J = 1, NAT
C           NSP2 = ISPE(J)
            CHARA2 = SPE(ISPE(J))
            CALL DFTB_CNVSQ(I,J,NSEQ)
            DIST = DISTMAT(NSEQ)
            IF (DAMPXH.AND.(CHARA1.EQ.H.OR.CHARA2.EQ.H)) THEN
              HDAMP = .TRUE.
              CALL DFTB_DMPCHK(ISPE(I),ISPE(J),HDAMP)
            ELSE
              HDAMP = .FALSE.
            END IF
            DO II = 1, MAXANG(ISPE(I))
              ISH = ISH0 + II
              HUBA = HUBBL(II,ISPE(I))
              NSP1 = ISHELL(ISH)
              HUBDERTMP = HUBDER(NSP1)
              DO JJ = 1, MAXANG(ISPE(J))
                JSH = JSH0 + JJ
                HUBB = HUBBL(JJ,ISPE(J))
                NSP2 = ISHELL(JSH)
C               WRITE (*,*) I, J
                CALL DFTB_CNVSQ(ISH,JSH,NSEQ)
                CALL DFTB_GAMDERF(GAMMA(NSEQ),GAMMA3(ISH,JSH),HUBA,HUBB,
     *            HUBDERTMP,DIST,DAMPXHE,TMP(1,NSP1,NSP2),DFTB3,HDAMP)
C      write (*,'(6i3,f20.10)') i,j,ish,jsh,nsp1,nsp2,gamma(nseq)
              END DO
            END DO
            JSH0 = JSH0 + MAXANG(ISPE(J))
          END DO
          ISH0 = ISH0 + MAXANG(ISPE(I))
        END DO
      ELSE
        DO I = 1, NSPE
          HUBA = HUBBL(1,I)
          DO J = 1, NSPE
            HUBB = HUBBL(1,J)
            CALL DFTB_PREP_GAMMA_GRAD(HUBA,HUBB,TMP(1,I,J),DFTB3,DAMPXH)
          END DO
        END DO
C
        DO I = 1, NAT
          NSP1 = ISPE(I)
          HUBA = HUBBL(1,NSP1)
          CHARA1 = SPE(ISPE(I))
          DO J = 1, NAT
            NSP2 = ISPE(J)
            HUBB = HUBBL(1,NSP2)
            CHARA2 = SPE(ISPE(J))
            CALL DFTB_CNVSQ(I,J,NSEQ)
            DIST = DISTMAT(NSEQ)
            IF (DAMPXH.AND.(CHARA1.EQ.H.OR.CHARA2.EQ.H)) THEN
              HDAMP = .TRUE.
              CALL DFTB_DMPCHK(NSP1,NSP2,HDAMP)
            ELSE
              HDAMP = .FALSE.
            END IF
            CALL DFTB_GAMDERF(GAMMA(NSEQ),GAMMA3(I,J),HUBA,HUBB,
     *        HUBDER(NSP1),DIST,DAMPXHE,TMP(1,NSP1,NSP2),DFTB3,HDAMP)
C      write (*,'(6i3,f20.10)') i,j,ish,jsh,nsp1,nsp2,gamma(nseq)
          END DO
        END DO
      END IF
C
      RETURN
C
      END SUBROUTINE DFTB_GAMDER
C
C-----------------------------------------------------------------------
C*MODULE DFTBGR    *DECK DFTB_EXTSDER
C>
C>    @brief Gradient of DFTB
C>
C>    @details Calculates 1st-order derivative of overlap matrix in AO
C>             basis.  Maybe called only from HOPDER?
C>
C>    @author Yoshio Nishimoto
C>    - Sep, 2014- Subroutine written
C>    @date Feb, 2016- Yoshio Nishimoto
C>    - Analytic derivatives are now used.
C>
C>           --- INPUT ---
C>    @param L1 Number of basis functions
C>    @param L2 L2=L1*(L1+1)/2
C>    @param NSTART Lower limit of atom index to be calculated
C>    @param NLAST Upper limit of atom index to be calculated
C>           --- OUTPUT ---
C>    @param DSINT Overlap derivative
C>
C
      SUBROUTINE DFTB_EXTSDER(DSINT,L2,NSTART,NLAST)
      use mx_limits, only: mxatm
      USE DFTBPB_MOD,ONLY: PERIOD,INDBND,PBCBOX,DFTB_INTPAIR
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER (MXSPE=10)
C
      COMMON /DFTBPR/ ETEMP,DFTBDP(MXSPE*14),DAMPXHE,HUBDER(MXSPE),
     *                ZREF(MXATM),SPNCST(6,MXSPE),SPE(MXATM),NSPE,
     *                MAXANG(MXATM),ISPE(MXATM),IND(MXATM+1),IDFTBD,
     *                PARAMDIR
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
C
      DOUBLE PRECISION :: AUP(9,9),BUP(9,9),VEC(3),VECTMP(3)
      DOUBLE PRECISION :: DSINT(L2,3,NLAST-NSTART+1)
      dimension tmp(l2,3)
C
      INTEGER :: NSEQ
      LOGICAL :: SKP
C
C     ----- CALCULATE DS/DR NUMERICALLY -----
C     DS/DR = (BUP - BUM)/DELTA * 0.5
C
      CALL VCLR(DSINT,1,L2*3*NAT)
C
      IF (PERIOD) THEN
        CALL DFTB_INTPAIR(C,NAT,ISPE,MAXANG,NAOPBC)
        NNN = 4210
        CALL VCLR(BUP,1,81)
        DO NI = 1, NAT
          NINTPAIR = INDBND(1,1,NI)
          DO JJ = 2, NINTPAIR !! skip JJ=1, because the same atom
            NJ = INDBND(1,JJ+1,NI)
            IF (NI.EQ.NJ) CYCLE
            NNN = INDBND(2,JJ+1,NI)
            IXTRV =  NNN/400-10
            IYTRV = (NNN-(IXTRV+10)*400)/20-10
            IZTRV = (NNN-(IXTRV+10)*400-(IYTRV+10)*20)-10
            VEC(1) = C(1,NJ) - C(1,NI) + PBCBOX(1,1)*IXTRV
     *             + PBCBOX(2,1)*IYTRV + PBCBOX(3,1)*IZTRV
            VEC(2) = C(2,NJ) - C(2,NI) + PBCBOX(1,2)*IXTRV
     *             + PBCBOX(2,2)*IYTRV + PBCBOX(3,2)*IZTRV
            VEC(3) = C(3,NJ) - C(3,NI) + PBCBOX(1,3)*IXTRV
     *             + PBCBOX(2,3)*IYTRV + PBCBOX(3,3)*IZTRV
            DO KI = 1, 3
              CALL DFTB_SKODE_HELPA(1,NI,NJ,KI,0,0,0,VEC,AUP,BUP,
     *                              .FALSE.,.TRUE.,SKP)
              IF (SKP) CYCLE
              DO M = 1, IND(NJ+1)-IND(NJ)
                MU = IND(NJ)+M
                DO N = 1, IND(NI+1)-IND(NI)
                  NU = IND(NI)+N
                  CALL DFTB_CNVSQ(MU,NU,NSEQ)
                  DSINT(NSEQ,KI,NI) = DSINT(NSEQ,KI,NI) + BUP(N,M)
                END DO
              END DO
            END DO
          END DO
        END DO
      ELSE
        CALL VCLR(BUP,1,81)
        DO NI = 1, NAT
          DO KI = 1, 3
            DO NJ = 1, NAT
              IF (NI.EQ.NJ) CYCLE
              VEC(1) = C(1,NJ) - C(1,NI)
              VEC(2) = C(2,NJ) - C(2,NI)
              VEC(3) = C(3,NJ) - C(3,NI)
              CALL DFTB_SKODE_HELPA(1,NI,NJ,KI,0,0,0,VEC,AUP,BUP,
     *                              .FALSE.,.TRUE.,SKP)
              IF (SKP) CYCLE
              DO M = 1, IND(NJ+1)-IND(NJ)
                MU = IND(NJ)+M
                DO N = 1, IND(NI+1)-IND(NI)
                  NU = IND(NI)+N
                  CALL DFTB_CNVSQ(MU,NU,NSEQ)
                  DSINT(NSEQ,KI,NI) = BUP(N,M) !! SDERAO
                END DO
              END DO
            END DO
          END DO
        END DO
      END IF
C
      RETURN
C
      END SUBROUTINE DFTB_EXTSDER
C
C-----------------------------------------------------------------------
C*MODULE DFTBGR    *DECK DFTB_GAMD3
C>
C>    @brief Gradient of DFTB
C>
C>    @details Calculates gamma 1st-order derivative of DFTB3.
C>             It seems that it is only used in FMO-DFTB
C>
C>    @author Yoshio Nishimoto
C>    - Sep, 2014- Subroutine written
C>    @date Feb, 2016- Yoshio Nishimoto
C>    - Bug fix
C>
C>           --- INPUT ---
C>    @param DIST Inter-atomic distances of atom A and B
C>    @param HUBA,HUBB Hubbard values of atom A and B
C>    @param HUBDER1,HUBDER2 Hubbard derivatives of atom A and B
C>    @param VALUE1,VALUE2 Constants, they are prepared in DFTB_GAMFMO
C>    @param DAMPXH Whether H-X damping or not
C>    @param DAMPXHE Exponent of X-H damping
C>    @param SPE1,SPE2 Name of atom A and B
C>           --- OUTPUT ---
C>    @param VAL1,VAL2 Derivatives of DFTB3 gamma
C>
C
      SUBROUTINE DFTB_GAMD3(DIST,DISTI,VAL1,VAL2,HUBA,HUBB,HUBDER1,
     * HUBDER2,VALUE1,VALUE2,DAMPXH,DAMPXHE,SPE1,SPE2)
C
      IMPLICIT NONE
C
      DOUBLE PRECISION, INTENT(IN) :: DIST,DISTI,HUBA,HUBB,VALUE1(11),
     *  VALUE2(11),DAMPXHE,HUBDER1,HUBDER2,SPE1,SPE2
      DOUBLE PRECISION, INTENT(OUT) :: VAL1,VAL2
      LOGICAL, INTENT(IN) :: DAMPXH
C
      DOUBLE PRECISION, PARAMETER :: ZERO=0.0D+00,ONE=1.0D+00, 
     * HALF=5.0D-01,TWO=2.0D+00!! ,ONE_THIRD=ONE/3.0D+00
      DOUBLE PRECISION, PARAMETER :: DISTTOL=1.0D-08,HUBTOL=1.0D-04
C
      DOUBLE PRECISION :: TAUA,TAUB,EXPA,EXPB,DIST2,DIST3,
     * DFAB_DA,DFBA_DA,DFAB_DB,DFBA_DB,FAB,FBA,HUBMEAN,
     * HUBMEANEXP,HUBMEANEXP1,EXPH,DH_DU,D2FAB_DBDR,D2FBA_DBDR,
     * DS_DR,GAMTMP3BA,DISTI2,GA,DGA_DR,GAMTMP3,DS_DU,DH_DR,D2H_DUDR,S,
     * DFAB_DR,DFBA_DR,D2FAB_DADR,D2FBA_DADR,GAMTMP4,GAMTMP5,SBA,
     * DS_DUBA,DS_DRBA
      DOUBLE PRECISION H/8HH       /
C
      TAUA = 3.2D+00*HUBA
      TAUB = 3.2D+00*HUBB
      IF (DIST.LE.DISTTOL) THEN
        VAL1 = ZERO
        VAL2 = ZERO
      ELSE
        DIST2 = DIST**2
C       DISTI = ONE/DIST
        DISTI2 = DISTI**2
        EXPA = EXP(-TAUA*DIST)
        IF (ABS(HUBA-HUBB).LE.HUBTOL) THEN
          !! CALCULATE HUBA EQ HUBB
          GA = DISTI+VALUE1(3)+VALUE1(1)*DIST+VALUE1(4)*DIST2
          DGA_DR = -DISTI2+VALUE1(1)+VALUE1(2)*DIST
          DS_DR = EXPA*DGA_DR-TAUA*EXPA*GA
          DIST3 = DIST2*DIST
          GAMTMP3 = EXPA*(VALUE1(6)*DIST + VALUE1(7)*DIST2
     *            + VALUE1(8)*DIST3)
          GAMTMP3 = -GAMTMP3
          IF (DAMPXH) THEN !! damp correction with 3rd
            IF (SPE1.EQ.H.OR.SPE2.EQ.H) THEN
              HUBMEAN = HUBA
              HUBMEANEXP = HUBMEAN**DAMPXHE
              HUBMEANEXP1 = HUBMEAN**(DAMPXHE-ONE)
              EXPH = EXP(-HUBMEANEXP*DIST2)
              DS_DU = EXPA*(-ONE+VALUE1(9)*DIST
     *          + VALUE1(10)*DIST2+VALUE1(11)*DIST3)
              DH_DR = -TWO*DIST*HUBMEANEXP*EXPH
              DH_DU = -HALF*DIST2*DAMPXHE*HUBMEANEXP1*EXPH
              D2H_DUDR = DAMPXHE*DIST*HUBMEANEXP1
     *                   *(DIST2*HUBMEANEXP-ONE)*EXPH
              S = EXPA * GA
              GAMTMP4 = DS_DU*DH_DR + DS_DR*DH_DU + S*D2H_DUDR
              GAMTMP3 = GAMTMP3*EXPH - GAMTMP4
            END IF
          END IF
          VAL1 = GAMTMP3 * HUBDER1! * ONE_THIRD
          VAL2 = VAL1
        ELSE
          !! CALCULATE HUBA NE HUBB
          EXPB = EXP(-TAUB*DIST)
          FAB = VALUE1(2)-VALUE1(1)*DISTI
          FBA = VALUE1(4)-VALUE1(3)*DISTI
          DFAB_DR = VALUE1(1)*DISTI2
          DFBA_DR = VALUE1(3)*DISTI2
          !! CALCULATE -dSf/dR
          DS_DR   = EXPA*DFAB_DR - TAUA*EXPA*FAB
     *            + EXPB*DFBA_DR - TAUB*EXPB*FBA
          DS_DRBA = EXPB*DFBA_DR - TAUB*EXPB*FBA
     *            + EXPA*DFAB_DR - TAUA*EXPA*FAB
          DFAB_DA = VALUE1(5) + VALUE1(6)*DISTI
          DFBA_DA = VALUE1(7) + VALUE1(8)*DISTI
          DFAB_DB = VALUE2(7) + VALUE2(8)*DISTI
          DFBA_DB = VALUE2(5) + VALUE2(6)*DISTI
          D2FAB_DADR = -VALUE1(6)*DISTI2
          D2FBA_DADR = -D2FAB_DADR
          D2FBA_DBDR = -VALUE2(6)*DISTI2
          D2FAB_DBDR = -D2FBA_DBDR
          GAMTMP3 = EXPA*((TAUA*DIST-ONE)*FAB-TAUA*DFAB_DA
     *      + D2FAB_DADR - DIST*DFAB_DR)
     *      + EXPB*(D2FBA_DADR - TAUB*DFBA_DA)
          GAMTMP3 = -3.2D+00*GAMTMP3
          GAMTMP3BA = EXPB*((TAUB*DIST-ONE)*FBA-TAUB*DFBA_DB
     *      + D2FBA_DBDR - DIST*DFBA_DR)
     *      + EXPA*(D2FAB_DBDR - TAUA*DFAB_DB)
          GAMTMP3BA = -3.2D+00*GAMTMP3BA
          IF (DAMPXH) THEN  !! damp correction with 3rd
            IF (SPE1.EQ.H.OR.SPE2.EQ.H) THEN
              HUBMEAN = HALF*(HUBA+HUBB)
              HUBMEANEXP = HUBMEAN**DAMPXHE
              EXPH = EXP(-HUBMEANEXP*DIST2)
              HUBMEANEXP1 = HUBMEAN**(DAMPXHE-ONE)
              DS_DU   = EXPA*(DFAB_DA-DIST*FAB)+EXPB*DFBA_DA
              DS_DUBA = EXPB*(DFBA_DB-DIST*FBA)+EXPA*DFAB_DB
              DS_DU   = 3.2D+00*DS_DU
              DS_DUBA = 3.2D+00*DS_DUBA
              DH_DR = -TWO*DIST*HUBMEANEXP*EXPH
              DH_DU = -HALF*DIST2*DAMPXHE*HUBMEANEXP1*EXPH
              D2H_DUDR = DAMPXHE*DIST*HUBMEANEXP1
     *                   *(DIST2*HUBMEANEXP-ONE)*EXPH
              S   = EXPA*FAB + EXPB*FBA
              SBA = EXPB*FBA + EXPA*FAB
              GAMTMP4   = DS_DU  *DH_DR + DS_DR  *DH_DU + S  *D2H_DUDR
              GAMTMP5   = DS_DUBA*DH_DR + DS_DRBA*DH_DU + SBA*D2H_DUDR
              GAMTMP3   = GAMTMP3*EXPH - GAMTMP4
              GAMTMP3BA = GAMTMP3BA*EXPH - GAMTMP5
            END IF
          END IF
          VAL1 = GAMTMP3 * HUBDER1! * ONE_THIRD
          VAL2 = GAMTMP3BA * HUBDER2! * ONE_THIRD
        END IF
      END IF
C
C
      END SUBROUTINE DFTB_GAMD3
C
C-----------------------------------------------------------------------
C*MODULE DFTBGR    *DECK DFTB_DISP_GRAD
C>
C>    @brief Gradient of DFTB
C>
C>    @details Calculates 1st-order derivative of dispersion correction
C>
C>    @author Yoshio Nishimoto
C>    - Sep, 2014- Subroutine written
C>    @date Feb, 2016- Yoshio Nishimoto
C>    - Added Slater-Kirkwoord (hybrid-polarized)
C>
C>           --- INPUT ---
C>    @param ZREF    Nuclear charges of a given atom. In FMO-DFTB, ZREF
C>                   of BDA or BAA is differet from QREF
C>    @param QREF    Nuclear charges of a given species. This is taken
C>                   from Slater-Koster files
C>    @param DISTMAT Matrix of inter-atomic distances
C>    @param IWRK    Integer working array
C>           --- OUTPUT ---
C>    @param EGRAD   Gradient
C>
C
      SUBROUTINE DFTB_DISP_GRAD(QREF,DISTMAT,EGRAD,IWRK)
      use mx_limits, only: mxatm

      IMPLICIT NONE
C
      INTEGER,PARAMETER :: MXSPE=10
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     & ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      INTEGER :: NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,IAN
      DOUBLE PRECISION :: ZAN,C,PARAMDIR
C
      COMMON /DFTBPR/ ETEMP,DFTBDP(MXSPE*14),DAMPXHE,HUBDER(MXSPE),
     *                ZREF(MXATM),SPNCST(6,MXSPE),SPE(MXATM),NSPE,
     *                MAXANG(MXATM),ISPE(MXATM),IND(MXATM+1),IDFTBD,
     *                PARAMDIR
      DOUBLE PRECISION :: ETEMP,DFTBDP,DAMPXHE,HUBDER,ZREF,SPNCST,SPE
      INTEGER :: NSPE,MAXANG,ISPE,IND,IDFTBD
C
      INTEGER :: I,J,K,NSP1,NSP2,NSEQ
      DOUBLE PRECISION :: R5,R6,R10,R12,C1,C2,CR0,RI,RJ,DI,DJ,DIST,
     * RIJ,DIJ,R0,TMP,VAL,RDMAT(3,NSPE,NSPE),CI,CJ,C6A,C6B,C6AB,R0A,R0B,
     * R0AB,NEA,PA,NEB,PB,R0AB7I,DISTI,DIST6I,DIST7I,DAMP1,DAMP3,
     * VEXP,DIST7,R7R07
      DOUBLE PRECISION, INTENT(IN) :: QREF(*),DISTMAT(*)
      DOUBLE PRECISION, INTENT(INOUT) :: EGRAD(3,*)
      INTEGER, INTENT(INOUT) :: IWRK(*)
C
      NSP1 = 0
      NSP2 = 0
      IF (IDFTBD.EQ.1) THEN !! "UFF"
C       C0 = 1.584D+01 !! 396/25
        C1 = 4.78947150872246D+01 !! 2^(5/6)*672/25
        C2 = 3.50498152274578D+01 !! 2^(2/3)*552/25
        CR0 = 8.908987181D-01
        !! FIRST CALCULATE SOME CONSTANTS FOR ALL PAIRS OF SPECIES
        DO I = 1, NSPE
          RI = DFTBDP(I*2-1)
          DI = DFTBDP(I*2  )
          DO J = I, NSPE
            RJ = DFTBDP(J*2-1)
            DJ = DFTBDP(J*2  )
            RIJ = SQRT(RI*RJ) !! van der Waals distance
            DIJ = SQRT(DI*DJ) !! well depth
            R0 = CR0*RIJ
            RDMAT(1,I,J) = R0
            RDMAT(2,I,J) = RIJ
            RDMAT(3,I,J) = DIJ
            RDMAT(1,J,I) = R0
            RDMAT(2,J,I) = RIJ
            RDMAT(3,J,I) = DIJ
          END DO
        END DO
        !! CALCULATE DISPERSION CORRECTION
        DO I = 1, NAT
          NSP1 = ISPE(I)
          CI = ZREF(I)/QREF(NSP1)
          DO J = I+1, NAT
            NSP2 = ISPE(J)
            CJ = ZREF(J)/QREF(NSP2)
            R0  = RDMAT(1,NSP1,NSP2)
            RIJ = RDMAT(2,NSP1,NSP2)
            DIJ = RDMAT(3,NSP1,NSP2)
            CALL DFTB_CNVSQ(I,J,NSEQ)
            DIST = DISTMAT(NSEQ)
            DIJ = DIJ * CI * CJ
            IF (DIST.GE.R0) THEN !! ATTRACTIVE POTENTIAL
              TMP = RIJ/DIST
              R6 = TMP**6
              R12 = R6**2
              VAL = 1.2D+01*DIJ*(R6-R12)/(DIST**2)
              DO K = 1, 3
                EGRAD(K,I) = EGRAD(K,I) + (C(K,I)-C(K,J))*VAL
                EGRAD(K,J) = EGRAD(K,J) - (C(K,I)-C(K,J))*VAL
              END DO
            ELSE !! REPULSIVE POTENTIAL
              TMP = DIST/RIJ
              R5 = TMP**5
              R10 = R5**2
              VAL = DIJ*(-C1*R5*5.0D+00+C2*R10*1.0D+01)/(DIST**2)
              DO K = 1, 3
                EGRAD(K,I) = EGRAD(K,I) + (C(K,I)-C(K,J))*VAL
                EGRAD(K,J) = EGRAD(K,J) - (C(K,I)-C(K,J))*VAL
              END DO
            END IF
          END DO
        END DO
      ELSE IF (IDFTBD.EQ.2) THEN !! "SK"
        !! FIRST CALCULATE SOME CONSTANTS FOR ALL PAIRS OF SPECIES
        DO I = 1, NSPE
          PA  = DFTBDP(I*3-2)
          R0A = DFTBDP(I*3-1)
          NEA = DFTBDP(I*3  )
          C6A = 7.5D-01*SQRT(NEA*PA*PA*PA)
          DO J = I, NSPE
            PB  = DFTBDP(J*3-2)
            R0B = DFTBDP(J*3-1)
            NEB = DFTBDP(J*3  )
            C6B = 7.5D-01*SQRT(NEB*PB*PB*PB)
C
            C6AB = (2.0D+00*C6A*C6B*PA*PB) / (PA*PA*C6B+PB*PB*C6A)
            R0AB = (R0A**3 + R0B**3) / (R0A**2 + R0B**2)
            R0AB7I = 1.0D+00 / R0AB**7
            RDMAT(1,I,J) = C6AB
            RDMAT(2,I,J) = R0AB7I
            RDMAT(1,J,I) = C6AB
            RDMAT(2,J,I) = R0AB7I
          END DO
        END DO
        DO I = 1, NAT
          NSP1 = ISPE(I)
          CI = ZREF(I)/QREF(NSP1)
          DO J = I+1, NAT
            NSP2 = ISPE(J)
            CJ = ZREF(J)/QREF(NSP2)
            C6AB   = RDMAT(1,NSP1,NSP2)
            R0AB7I = RDMAT(2,NSP1,NSP2)
            CALL DFTB_CNVSQ(I,J,NSEQ)
            DIST = DISTMAT(NSEQ)
            C6AB = C6AB * CI * CJ
            DISTI = 1.0D+00 / DIST
            DIST6I = DISTI**6
            DIST7  = DIST**7
            DIST7I = DIST6I*DISTI
            R7R07 = DIST7*R0AB7I
            VEXP = EXP(-3.0D+00*R7R07)
            DAMP1 = 1.0D+00 - VEXP
            DAMP3 = DAMP1**3
            VAL = 8.4D+01*VEXP*R7R07-6.0D+00*DAMP1
            VAL = -VAL*C6AB*DAMP3*DIST7I*DISTI
            DO K = 1, 3
              EGRAD(K,I) = EGRAD(K,I) + (C(K,I)-C(K,J))*VAL
              EGRAD(K,J) = EGRAD(K,J) - (C(K,I)-C(K,J))*VAL
            END DO
          END DO
        END DO
      ELSE IF (IDFTBD.EQ.4) THEN !! "SKHP"
        IF (IWRK(NAT+1).EQ.0) THEN
          !! COUNT NUMBER OF BONDS FOR ALL ATOMS
          CALL VICLR(IWRK,1,NAT)
          NSEQ = 0
          DO I = 1, NAT
            NSP1 = ISPE(I)
            R0A = DFTBDP(NSP1*14-13) !! COVALENT RADIUS
            DO J = 1, I-1
              NSP2 = ISPE(I)
              R0B = DFTBDP(NSP2*14-13)
              NSEQ = NSEQ + 1
              DIST = DISTMAT(NSEQ)
              IF (DIST.LT.(R0A+R0B)) THEN
                IWRK(I) = IWRK(I) + 1
                IWRK(J) = IWRK(J) + 1
              END IF
            END DO
            NSEQ = NSEQ + 1 !! SKIP THE INDEX FOR SAME ATOM
          END DO
          DO I = 1, NAT
            IF (IWRK(I).GT.5) IWRK(I) = 5 !! MAX NEIGHBORS SHOULD BE 5
          END DO
        END IF
        !! CALCULATE DISPERSION CORRECTION
        NSEQ = 0
        DO I = 1, NAT
          NSP1 = ISPE(I)
          PA   = DFTBDP(NSP1*14-12+IWRK(I))
          R0A  = DFTBDP(NSP1*14- 6+IWRK(I))
          NEA  = DFTBDP(NSP1*14)
          C6A  = 7.5D-01*SQRT(NEA*PA*PA*PA)
          CI   = ZREF(I)/QREF(NSP1)
          DO J = 1, I-1
            NSP2 = ISPE(J)
            PB   = DFTBDP(NSP2*14-12+IWRK(J))
            R0B  = DFTBDP(NSP2*14- 6+IWRK(J))
            NEB  = DFTBDP(NSP2*14)
            C6B  = 7.5D-01*SQRT(NEB*PB*PB*PB)
            CJ   = ZREF(J)/QREF(NSP2)
            C6AB = (2.0D+00*C6A*C6B*PA*PB) / (PA*PA*C6B+PB*PB*C6A)
            R0AB = (R0A**3 + R0B**3) / (R0A**2 + R0B**2)
            R0AB7I = 1.0D+00 / R0AB**7
            NSEQ = NSEQ + 1
            DIST = DISTMAT(NSEQ)
            C6AB = C6AB * CI * CJ
            DISTI = 1.0D+00 / DIST
            DIST6I = DISTI**6
            DIST7  = DIST**7
            DIST7I = DIST6I*DISTI
            R7R07 = DIST7*R0AB7I
            VEXP = EXP(-3.0D+00*R7R07)
            DAMP1 = 1.0D+00 - VEXP
            DAMP3 = DAMP1**3
            VAL = 8.4D+01*VEXP*R7R07-6.0D+00*DAMP1
            VAL = -VAL*C6AB*DAMP3*DIST7I*DISTI
            DO K = 1, 3
              EGRAD(K,I) = EGRAD(K,I) + (C(K,I)-C(K,J))*VAL
              EGRAD(K,J) = EGRAD(K,J) - (C(K,I)-C(K,J))*VAL
            END DO
          END DO
          NSEQ = NSEQ + 1
        END DO
      END IF
C
      END SUBROUTINE DFTB_DISP_GRAD
C
C-----------------------------------------------------------------------
C*MODULE DFTBGR    *DECK DFTB_REPGRAD
C>
C>    @brief Calculate gradient for DFTB
C>
C>    @details Calculate gradient for DFTB
C>
C>    @author Yoshio Nishimoto
C>    - Sep, 2014- Subroutine written
C>
C>           --- INPUT ---
C>    @param NAT      Number of atoms
C>    @param C        Coordinate
C>    @param ISPE     Index of species
C>    @param NSPE     Number of species
C>    @param NUMREP   Number of grids for specific pairs
C>    @param REPCUT   Cutoff distance for specific pairs
C>    @param REPINTV  Lower and upper bound of each grid of repulsive
C>                    potential for specific pairs
C>    @param REPSHORT Coefficients (and exponent) for short range
C>                    respulsive potential
C>    @param REPSHORT Coefficients for medium (4 terms) and long (6
C>                    terms) range respulsive potential
C>           --- OUTPUT ---
C>    @param GRAD     Gradient of repulsive potential
C>    @param DISTMAT  Matrix of distance between each atom (up. triang.)
C>
C
      SUBROUTINE DFTB_REPGRAD(NAT,C,ISPE,NSPE,NUMREP,REPCUT,
     & REPINTV,REPSHORT,REPCOEFF,GRAD,DISTMAT,ZREF,QREF)
C
      USE DFTBPB_MOD,ONLY: PBCBOX,INDBND,MAXINT,PERIOD,STRESS
C
      IMPLICIT NONE
C
      INTEGER,PARAMETER :: MXINT=150,MXSPE=10
C
      INTEGER,INTENT(IN) :: NAT
      DOUBLE PRECISION,INTENT(IN) :: C(3,NAT)
      INTEGER,INTENT(IN) :: ISPE(NAT)
      INTEGER,INTENT(IN) :: NSPE
      INTEGER,INTENT(IN) :: NUMREP(MXSPE,MXSPE)
      DOUBLE PRECISION,INTENT(IN) :: REPCUT(MXSPE,MXSPE)
      DOUBLE PRECISION,INTENT(IN) :: REPINTV(2,MXINT,NSPE,NSPE)
      DOUBLE PRECISION,INTENT(IN) :: REPSHORT(3,NSPE,NSPE)
      DOUBLE PRECISION,INTENT(IN) :: REPCOEFF(6,MXINT,NSPE,NSPE)
      DOUBLE PRECISION,INTENT(IN) :: ZREF(NAT),QREF(*)
C
      DOUBLE PRECISION,INTENT(INOUT) :: GRAD(3,NAT)
      DOUBLE PRECISION,INTENT(INOUT) :: DISTMAT(*)
C
      INTEGER :: I,J,K,L,SPE1,SPE2,NSEQ,NNN,JJ,IXTRV,IYTRV,IZTRV,
     *           NINTPAIR
      DOUBLE PRECISION :: VAL,X,X1,DIST,CI,CJ,TQ1,TQ2,VEC1(3),VEC2(3),
     *                    SCAL
C
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
      INTEGER :: LFFATMNAM,LFFCORD,LFFZANF
      INTEGER :: LFFZMAS,LFFONEMAS,LFFQMZMAS,LFFQM1MAS
      INTEGER :: LFFCHARG,LFFPOL,LFFDIP
      INTEGER :: LFFFIELD1,LFFFIELD2,LFFFIELD3
      INTEGER :: LFFSIG,LFFEPS,LFFSIG2,LFFEPS2
      INTEGER :: LFFBOND0,LFFFCBOND
      INTEGER :: LFFANGL0,LFFFCANGL,LFFFCWAGG
      INTEGER :: LFFDIHB0,LFFFCDIHB
      INTEGER :: LFFVROT,LFFNNN,LFFGAMA,LFFIPAIR
      INTEGER :: LFFKLIST,LFFLLIST,LFFL1213J,LFFL14J
      INTEGER :: LFFMLIST,LFFNLIST,LFFLKQMMM
      INTEGER :: LFFVEL,LFFQMVEL
      INTEGER :: LFFFFGRD0,LFFFFGRD1,LFFFFGRD2
      INTEGER :: LFFQMGRD0,LFFQMGRD1,LFFQMGRD2,LFFDETMP
      INTEGER :: LFFCLPR,LFFZLPR,LFFNLPR
      INTEGER :: LFFXTS,LFFYTS,LFFZTS,LFFCMAT1
      INTEGER :: LFFQRXN1,LFFQRXN2,LFFPOT1,LFFPOT2,LFFQRXNMP
      INTEGER :: LFFQRXNTA,LFFQRXNXY,LFFNONLSTQ
      INTEGER :: LFFDIPMP,LFFDIPTA,LFFDIPXY,LFFLISTQM,LFFNONLS1
      INTEGER :: LFFMAPLST,LFFCMAPCO
      COMMON /FFPARM/ NFFAT,NBOND,NANGL,NDIHR,NDIHB,NCMAP,NWAGG,
     *                N1213J,N14J,NLKQMM,IDOCHG,IDOPOL,IDOLJ,IDOCMAP
      INTEGER :: NFFAT,NBOND,NANGL,NDIHR,NDIHB,NCMAP,NWAGG
      INTEGER :: N1213J,N14J,NLKQMM,IDOCHG,IDOPOL,IDOLJ,IDOCMAP
      INTEGER :: ILK, LK1, LK2
      COMMON /FMCOM / XX(1)
      DOUBLE PRECISION :: SHRINK, XX
C
      NNN = 4210
      SCAL = 1.0D+00
      DO I = 1, NAT
        SPE1 = ISPE(I)
        TQ1 = 1.0D+00 / QREF(SPE1)
        CI = ZREF(I)*TQ1
        CALL DCOPY(3,C(1,I),1,VEC1,1)
        IF (MAXINT.EQ.0) THEN
          NINTPAIR = I-1
        ELSE
          NINTPAIR = INDBND(2,1,I)
        END IF
        DO JJ = 1, NINTPAIR
          IF (MAXINT.EQ.0) THEN
            J = JJ
          ELSE
            J = INDBND(1,JJ+1,I)
            IF (I.LT.J) CYCLE
          END IF
          CALL DCOPY(3,C(1,J),1,VEC2,1)
          SPE2 = ISPE(J)
          TQ2 = 1.0D+00 / QREF(SPE2)
          CJ = ZREF(J)*TQ2
          IF (PERIOD) THEN
            NNN = INDBND(2,JJ+1,I)
            IXTRV =  NNN/400-10
            IYTRV = (NNN-(IXTRV+10)*400)/20-10
            IZTRV = (NNN-(IXTRV+10)*400-(IYTRV+10)*20)-10
            VEC2(1) = VEC2(1) + PBCBOX(1,1)*IXTRV
     *              + PBCBOX(2,1)*IYTRV + PBCBOX(3,1)*IZTRV
            VEC2(2) = VEC2(2) + PBCBOX(1,2)*IXTRV
     *              + PBCBOX(2,2)*IYTRV + PBCBOX(3,2)*IZTRV
            VEC2(3) = VEC2(3) + PBCBOX(1,3)*IXTRV
     *              + PBCBOX(2,3)*IYTRV + PBCBOX(3,3)*IZTRV
            SCAL = 2.0D+00
            IF (I.EQ.J) SCAL = 1.0D+00
          END IF
          DIST = SQRT( (VEC1(1)-VEC2(1))**2
     *               + (VEC1(2)-VEC2(2))**2
     *               + (VEC1(3)-VEC2(3))**2 )
          IF (NNN.EQ.4210) THEN
            CALL DFTB_CNVSQ(I,J,NSEQ)
            DISTMAT(NSEQ) = DIST
          END IF
          SHRINK = 0.0D+00
          IF(NFFAT.GT.0.AND.NLKQMM.GT.0) THEN
            DO ILK=1,NLKQMM
              LK1 = INT(XX(LFFLKQMMM+2*(ILK-1)  )+1.0D-03)
              LK2 = INT(XX(LFFLKQMMM+2*(ILK-1)+1)+1.0D-03)
              IF((LK1.EQ.I.AND.LK2.EQ.J).OR.
     *           (LK1.EQ.J.AND.LK2.EQ.I)    ) THEN
                 SHRINK = 0.83D+00
              END IF
            ENDDO
          END IF
          DIST = DIST - SHRINK
          IF (DIST.LT.1.0D-02) THEN
C         ----- VAL = 0.0 IF DISTANCE IS LESS THAN 0.01 AU
            VAL = 0.0D+00
          ELSE IF (DIST.GT.REPCUT(SPE1,SPE2)) THEN
C         ----- VAL = 0.0 IF DISTANCE IS MORE THAN CUTOFF
            VAL = 0.0D+00
          ELSE
            IF (DIST.LT.REPINTV(1,1,SPE1,SPE2)) THEN
C            --- SHORT RANGE REPULSIVE PART ---
              VAL = -REPSHORT(1,SPE1,SPE2)
     &            * EXP(-REPSHORT(1,SPE1,SPE2)*DIST
     &                 + REPSHORT(2,SPE1,SPE2))
            ELSE
              DO K = 1, NUMREP(SPE1,SPE2)
                IF (DIST.GE.REPINTV(1,K,SPE1,SPE2) 
     &              .AND. DIST.LE.REPINTV(2,K,SPE1,SPE2)) EXIT
              END DO
              VAL = 0.0D+00 ! REPCOEFF(1,K,SPE1,SPE2)
              X1  = DIST - REPINTV(1,K,SPE1,SPE2)
              X   = 1.0D+00
              IF (K.LT.NUMREP(SPE1,SPE2)) THEN
C            --- MEDIUM RANGE REPULSIVE PART ---
                DO L = 2, 4
                  VAL = VAL + (L-1)*REPCOEFF(L,K,SPE1,SPE2) * X
                  X   = X * X1
                END DO
              ELSE
C            --- LONG RANGE REPULSIVE PART ---
                DO L = 2, 6
                  VAL = VAL + (L-1)*REPCOEFF(L,K,SPE1,SPE2) * X
                  X   = X * X1
                END DO
              END IF
            END IF
            VAL = VAL / (DIST+SHRINK) *CI*CJ
          END IF
C         WRITE (6,'(3I3,X,F14.7,X,F8.5)') I, J, K, VAL, DIST
          DO K = 1, 3
            GRAD(K,I) = GRAD(K,I) + (VEC1(K)-VEC2(K))*VAL
            GRAD(K,J) = GRAD(K,J) - (VEC1(K)-VEC2(K))*VAL
          END DO
          !! (dV/dR)*(dR/dH_{ax})*H_{bx}
          IF (PERIOD) THEN
            VAL = VAL*SCAL
            DO K = 1, 3
              DO L = 1, 3
                STRESS(K,L) = STRESS(K,L)
     *            + VAL*(VEC1(K)-VEC2(K))*(VEC1(L)-VEC2(L))
              END DO
            END DO
          END IF
        END DO
      END DO
C
      !! COMPLETE DISTMAT
C     IF (.NOT.PERIOD) THEN
      IF (MAXINT.NE.0) THEN
        NSEQ = 0
        DO I = 1, NAT
          DO J = 1, I
            NSEQ = NSEQ + 1
            IF (DISTMAT(NSEQ).EQ.0.0D+00) THEN
              DISTMAT(NSEQ) = SQRT( (C(1,I)-C(1,J))**2
     *                            + (C(2,I)-C(2,J))**2
     *                            + (C(3,I)-C(3,J))**2 )
            END IF
          END DO
        END DO
      END IF
C      
      END SUBROUTINE DFTB_REPGRAD
C
C-----------------------------------------------------------------------
C*MODULE DFTBGR    *DECK REPGRADFMO
C>
C>    @brief Calculate respulsive gradient
C>
C>    @details Calculate gradient repulsive terms in ES-DIM calculations
C>             for FMO-DFTB.  This function may possibly be required
C>             for ES-DIM calculation when repulsive term has values
C>             at long-range region.  We have assumed that inter-fragm.
C>             contribution in ES-DIM is (close to) zero, since the
C>             repulsive terms are much short-ranged.
C>
C>    @author Yoshio Nishimoto
C>    - Sep, 2014- Subroutine written
C>
C>           --- INPUT ---
C>    @param DIST Distance between atom A and B
C>    @param NUMREP Number of grids for specific pairs
C>    @param REPCUT Cutoff distance for specific pairs
C>    @param ISPE Index of species
C>    @param REPINTV Lower and upper bound of each grid of repulsive
C>                   potential for specific pairs
C>    @param REPSHORT Coefficients (and exponent) for short range
C>                    respulsive potential
C>    @param REPCOEFF Coefficients for medium (4 terms) and long (6
C>                    terms) range respulsive potential
C>           --- OUTPUT ---
C>    @param REPGRADFMO Gradient
C>
C
      DOUBLE PRECISION FUNCTION REPGRADFMO(DIST,NUMREP,
     & REPCUT,REPINTV,REPSHORT,REPCOEFF)
C
      IMPLICIT NONE
C
      INTEGER,PARAMETER :: MXINT=150
C
      DOUBLE PRECISION,INTENT(IN) :: DIST
      INTEGER,INTENT(IN) :: NUMREP
      DOUBLE PRECISION,INTENT(IN) :: REPCUT
      DOUBLE PRECISION,INTENT(IN) :: REPINTV(2,MXINT)
      DOUBLE PRECISION,INTENT(IN) :: REPSHORT(3)
      DOUBLE PRECISION,INTENT(IN) :: REPCOEFF(6,MXINT)
C
      INTEGER :: K, L
      DOUBLE PRECISION :: X, X1
C
      IF (DIST.LT.1.0D-02) THEN
C     ----- REPGRADFMO = 0.0 IF DISTANCE IS LESS THAN 0.01 AU
        REPGRADFMO = 0.0D+00
        RETURN
      ELSE IF (DIST.GT.REPCUT) THEN
C     ----- REPGRADFMO = 0.0 IF DISTANCE IS MORE THAN CUTOFF
        REPGRADFMO = 0.0D+00
        RETURN
      ELSE
        IF (DIST.LT.REPINTV(1,1)) THEN
C        --- SHORT RANGE REPULSIVE PART ---
          REPGRADFMO = -REPSHORT(1)
     &        * EXP(-REPSHORT(1)*DIST+ REPSHORT(2))
        ELSE
          DO K = 1, NUMREP
            IF (DIST.GE.REPINTV(1,K) .AND. DIST.LE.REPINTV(2,K)) EXIT
          END DO
          REPGRADFMO = 0.0D+00 ! REPCOEFF(1,K,SPE1,SPE2)
          X1  = DIST - REPINTV(1,K)
          X   = 1.0D+00
          IF (K.LT.NUMREP) THEN
C        --- MEDIUM RANGE REPULSIVE PART ---
            DO L = 2, 4
              REPGRADFMO = REPGRADFMO + (L-1)*REPCOEFF(L,K) * X
              X   = X * X1
            END DO
          ELSE
C        --- LONG RANGE REPULSIVE PART ---
            DO L = 2, 6
              REPGRADFMO = REPGRADFMO + (L-1)*REPCOEFF(L,K) * X
              X   = X * X1
            END DO
          END IF
        END IF
      END IF
C     WRITE (6,'(3I3,X,F14.7,X,F8.5)') I, J, K, REPGRADFMO, DIST
      REPGRADFMO = REPGRADFMO / DIST
C
      END FUNCTION REPGRADFMO
C
C-----------------------------------------------------------------------
C*MODULE DFTBGR    *DECK DFTB_PRJGRD
C>
C>    @brief   Projection of gradients
C>
C>    @details Projection of gradients for DFTB.  This subroutine should
C>             be called for DFTB/PCM calculations (gradient).
C>             Since DFTB does not call JKDER, the PRJGRD subroutine is
C>             never called.  For this reason, it has to be called later
C>             after the addition of PCM contributions.  This is copied
C>             from the JKDER subroutine.
C>
C>    @author Yoshio Nishimoto
C>    - Feb, 2016- Subroutine written
C>
C>           --- INPUT ---
C>    @param C      Coordinates
C>    @param RUNTYP Runtype.
C>    @param NAT    Number of aotms
C>           --- IN/OUTPUT ---
C>    @param EGRAD  Gradient on entry, and gradient whose rotational
C>                  vectors have been projected out on exit
C>
      SUBROUTINE DFTB_PRJGRD(EGRAD,C,RUNTYP,NAT)
      USE MX_LIMITS,ONLY:MXFRG,MXDFG,MXDPPT,mxatm
      USE comm_FRGINF
      USE DFTBPB_MOD,ONLY:PERIOD
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C
      LOGICAL PROJGRAD,LINEAR
C
      COMMON /DFTBHS/ DEGTHR(2),PFREQ(2),CPCONV,OCCTHR,MODHSS,
     *                MODPRJ,MXCPIT,ARAMAN
      COMMON /FMCOM / X(1)
      COMMON /FMOINF/ NFG,NLAYER,NATFMO,NBDFG,NAOTYP,NBODY,NSEGM
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON /MASSES/ ZMASS(MXATM)
      COMMON /RUNOPT/ RUNTYP0,EXETYP,NEVALS,NGLEVL,NHLEVL
      COMMON /ZMAT  / NZMAT,NZVAR,NVAR,NSYMC,LINEAR
C
      DATA OPT,SADPT,MD/8HOPTIMIZE,8HSADPOINT,8HMD      /
      DATA HESS/8HHESSIAN /
C
      DIMENSION EGRAD(3,NAT)
C
C            PROJECT ROTATIONAL CONTAMINANT FROM GRADIENTS
C       WHICH MAY HAVE TORQUES DUE TO FINITE GRIDS IN DFT OR PCM
C     GEOMETRY SEARCHES CAN SKIP IF THEY TRANSFORM INTO INTERNAL
C     COORDINATES, A HEALTHIER WAY TO REMOVE THE PROBLEM.
C     WHAT IS THE THINKING IS BEHIND EXCLUDING MD?
C     THIS WILL BE DONE LATER FOR THE FULL MOLECULE IN CASE OF FMO.
C
      PROJGRAD = .TRUE.
      IF(period) PROJGRAD=.FALSE.
C
      IF(RUNTYP.EQ.OPT   .AND. NZVAR.GT.0) PROJGRAD=.FALSE.
      IF(RUNTYP.EQ.SADPT .AND. NZVAR.GT.0) PROJGRAD=.FALSE.
      IF(RUNTYP.EQ.MD)                     PROJGRAD=.FALSE.
      IF(NFG.GT.0)                         PROJGRAD=.FALSE.
      IF(RUNTYP.EQ.HESS  .AND. NHLEVL.EQ.1)PROJGRAD=.FALSE.
c     IF(period.and.IAND(MODPRJ,8).NE.0)   PROJGRAD=.TRUE.
      IF(IAND(MODPRJ,4).NE.0) PROJGRAD=.NOT.PROJGRAD
C
      IF(PROJGRAD) THEN
c        write(6,*) 'Gradient is projected out.'
         NCCF = 3*NAT + 6*NFRG
         CALL VALFM(LOADFM)
         LGTOT  = LOADFM + 1
         LRM    = LGTOT  + 3*NAT +  6*NFRG
         LX     = LRM    + 3*NAT + 21*NFRG
         LP     = LX     + 3*NAT +  6*NFRG
         LAST   = LP     + 3*NAT +  6*NFRG 
         NEEDP  = LAST-LOADFM-1
         CALL GETFM(NEEDP)
         CALL PRJGRD(EGRAD,C,X(LX),X(LGTOT),X(LRM),X(LP),ZMASS,
     *               NAT,3*NAT,NCCF,.FALSE.)
         CALL RETFM(NEEDP)
         CALL DAWRIT(IDAF,IODA,EGRAD,3*NAT,3,0)
      END IF
C
      RETURN
C
      END SUBROUTINE DFTB_PRJGRD
C
C-----------------------------------------------------------------------
C*MODULE DFTBGR    *DECK DFTB_ZVEC
C>
C>    @brief   FMO-DFTB/AFO analytic derivative
C>
C>    @details Solve three types of Z-vector equations. Types are
C>             specified by MODE, and find the details below.
C>
C>    @author  Yoshio Nishimoto
C>    - Feb, 2018- Subroutine written
C>
C>           --- INPUT ---
C>    @param MODE     See below
C>    @param L1       Number of AOs
C>    @param L2       L1*(L1+1)/2
C>    @param NNUM     Number of independent elements
C>    @param NA       Number of alpha electrons
C>    @param NAT      Number of atoms
C>    @param NUMFRZ   Number of frozen occupied orbitals
C>    @param R0       Initial residues
C>    @param Z0       Preconditioned residues
C>    @param P0       Temporary Z
C>    @param AP       Work space for A*P0 (this should be P0*A)
C>    @param PREC     Preconditioner
C>    @param WRK      Work space
C>    @param WRK1     Work space
C>    @param WRK2     Work space
C>    @param WRK3     Work space
C>    @param VEC      MO coefficient
C>    @param FAO      Fock matrix in AO
C>    @param S        Overlap matrix in AO
C>    @param EIG      Orbital energy
C>    @param GAMMA    gamma for DFTB2
C>    @param GAMMA3   Gamma for DFTB3
C>    @param DQ       Mulliken charge (not population)
C>    @param NSPE     Number of species
C>    @param ISPE     Index of species
C>    @param IND      AO index of each atom
C>    @param IW       6
C>    @param MASWRK   Master node?
C>    @param VECL     LMO coefficient
C>    @param MDOUT    Controls the amount of output
C>    @param DFTB3    Whether DFTB3 or not
C>    @param SHIFT3   Interfragment contribution for DFTB3
C>    @param MODGAMMA How gamma/Gamma are stored
C>           --- OUTPUT ---
C>    @param Z1       Z after solving Z-vector (in MO)
C>    @param ZTMP     Z1 in AO
C>
C
      SUBROUTINE DFTB_ZVEC(MODE,L1,L2,NNUM,NA,NAT,NUMFRZ,
     *                     R0,Z0,P0,AP,Z1,PREC,ZTMP,
     *                     WRK,WRK1,WRK2,WRK3,
     *                     VEC,FAO,S,EIG,GAMMA,GAMMA3,DQ,
     *                     NSPE,ISPE,IND,MAXANG,IW,MASWRK,VECL,NORBPRJ,
     *                     MDOUT,DFTB3,SHIFT3,MODGAMMA)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER (ZERO=0.0D+00,ONE=1.0D+00,TWO=2.0D+00)
C
      DIMENSION R0(*),Z0(*),P0(L1-NA,NA-NUMFRZ),AP(L1-NA,NA-NUMFRZ),
     *          Z1(*),PREC(*),ZTMP(*),
     *          WRK(L1),WRK1(L1*L1),WRK2(L1*L1),WRK3(L1,L1),
     *          VEC(L1,L1),FAO(*),S(*),EIG(*),GAMMA(*),GAMMA3(NAT,NAT),
     *          DQ(*),VECL(L1,L1),SHIFT3(*)
      DIMENSION ISPE(*),IND(*),MAXANG(*)
      LOGICAL   MASWRK,MDOUT,DFTB3
C
C     MODE 1: CPHF TYPE Z-VECTOR WITH PROJECTION (Eq. (67))
C     MODE 2: CPHF TYPE STANDARD Z-VECTOR (Eq. (92))
C     MODE 3: CPL  TYPE Z-VECTOR (Eq. (87))
C
      CALL DCOPY(NNUM,Z0,1,P0,1)
      CALL VCLR(Z1,1,NNUM)
      IF (MODE.EQ.3) NLOC = NUMFRZ
C
      THRES=1.0D-08
C     THRES=1.0D-10
      DO ITER = 1, 50
        ! COMPUTE AP = A*P0
        IF (MODE.EQ.1.OR.MODE.EQ.2) THEN
          !! first, construct Z_{\mu \nu}
          CALL DGEMM('N','N',L1,NA-NUMFRZ,L1-NA,ONE,VEC(1,NA+1),L1,
     *               P0,L1-NA,ZERO,WRK2,L1)
          CALL DGEMM('N','T',L1,L1,NA-NUMFRZ,ONE,WRK2,L1,
     *               VEC(1,NUMFRZ+1),L1,ZERO,ZTMP,L1)
          IF (MODE.EQ.1 .AND. NORBPRJ.GT.0) THEN
          !! second, construct \widetilde{L}_{\mu \nu}^\prime for term 3
          !! This is only for the case with projection
            CALL DGEMM('N','N',L1,NORBPRJ,L1,ONE,ZTMP,L1,
     *                 WRK3(1,L1-NORBPRJ+1),L1,ZERO,WRK2,L1)
            CALL DCOPY(L1*NORBPRJ,VECL(1,L1-NORBPRJ+1),1,WRK1,1)
            DO KO = 1, NORBPRJ
              TMP = DDOT(L1,WRK2(1+L1*(KO-1)),1,WRK3(1,L1-NORBPRJ+KO),1)
              CALL DSCAL(L1,TMP,WRK1(1+L1*(KO-1)),1)
            END DO
            CALL DGEMM('N','T',L1,L1,NORBPRJ,-ONE,WRK1,L1,
     *                 VECL(1,L1-NORBPRJ+1),L1,ONE,ZTMP,L1)
          END IF
C
          !! term 2 + 3
          !! USE FAO AS A TEMPORARY STORAGE
          CALL SYMSQT(L1,ZTMP,L1,FAO,1)
          CALL DFTB_MULLIKEN(L1,L2,S,FAO,WRK2)
          DO NC = 1, NAT
            VAL = ZERO
            DO L = 1, IND(NC+1) - IND(NC)
              MU = IND(NC)+L
              VAL = VAL + WRK2(MU)
            END DO
            WRK1(NC) = VAL !! Q_C^Z - Q_C^{\prime Z}
          END DO
          CALL VCLR(WRK,1,NAT)
          CALL DFTB_2ND_SHIFT(WRK,WRK1,GAMMA,NAT)
          IF (DFTB3) THEN
            IF (MODE.EQ.1.AND.IAND(MODGAMMA,4096).EQ.0) THEN
              !! ADD INTER-FRAGMENT
              DO NI = 1, NAT
                WRK(NI) = WRK(NI) + TWO*SHIFT3(NI)*WRK1(NI)
              END DO
            END IF
            CALL DFTB_TD_3RD_SHIFT(WRK,WRK,WRK1,WRK1,DQ,GAMMA3,
     *           NAT,.TRUE.)
          END IF
          CALL VCLR(FAO,1,L2)
          CALL DFTB_SHIFT_FOCK(WRK,FAO,S,NAT,L2,NSPE,ISPE,
     *      IND,MAXANG)
          CALL CPYTSQ(FAO,WRK1,L1,1)
          CALL DGEMM('T','N',L1-NA,L1,L1,ONE,VEC(1,NA+1),L1,
     *               WRK1,L1,ZERO,WRK2,L1-NA)
          CALL DGEMM('N','N',L1-NA,NA-NUMFRZ,L1,ONE,WRK2,L1-NA,
     *               VEC(1,NUMFRZ+1),L1,ZERO,AP,L1-NA)
C
          CALL DSCAL(NNUM,-4.0D+00,AP,1)
          !! Term 1
          DO IO = NA+1, L1
            DO JO = NUMFRZ+1, NA
              AP(IO-NA,JO-NUMFRZ) = AP(IO-NA,JO-NUMFRZ)
     *          + (EIG(JO)-EIG(IO))*P0(IO-NA,JO-NUMFRZ)
            END DO
          END DO
        ELSE IF (MODE.EQ.3) THEN
          CALL ASYMTRZE(P0,NLOC,NLOC,ASYM)
          CALL DFTB_CCC(NLOC,WRK3,P0,AP)
        END IF
C
        DENO  = DDOT(NNUM,P0,1,AP,1)
        RDNO  = DDOT(NNUM,R0,1,Z0,1)
        IF (MASWRK) THEN
          RNORM = DDOT(NNUM,R0,1,R0,1)
          RNORM = SQRT(RNORM)
        END IF
        IF (MASWRK.AND.MDOUT) WRITE(IW,1000) ITER, RNORM
 1000   FORMAT(1X,I4,'-TH LOCAL CONJUGATE GRADIENT RMSD:',1F12.8)
        IF (RNORM.LT.THRES) THEN
          IF (MASWRK.AND.MDOUT) WRITE(IW,*) 'CG CONVERGES AT ',ITER
          GO TO 100
        END IF
C       IF (MODE.EQ.0) RETURN
        ALPHA1 = RDNO/DENO
        ! NEW X = ZVEC1
        CALL DAXPY(NNUM,ALPHA1,P0,1,Z1,1)
        ! NEW R
        CALL DAXPY(NNUM,-ALPHA1,AP,1,R0,1)
        ! NEW Z (SCALE BY PRECONDITIONER)
        DO I = 1, NNUM
          Z0(I) = R0(I)*PREC(I)
        END DO
C
        BETA1 = DDOT(NNUM,R0,1,Z0,1)/RDNO
        ! NEW P
        CALL DSCAL(NNUM,BETA1,P0,1)
        CALL DAXPY(NNUM,ONE,Z0,1,P0,1)
      END DO
  100 CONTINUE
C
C     Z-VECTOR ITERATION MUST BE CONVERGED
C
      !! CONSTUCT X(LZTMP) = Z_{\mu \nu}
      IF (MODE.EQ.1.OR.MODE.EQ.2) THEN
        CALL DGEMM('N','N',L1,NA-NUMFRZ,L1-NA,ONE,VEC(1,NA+1),L1,
     *             Z1,L1-NA,ZERO,WRK2,L1)
        CALL DGEMM('N','T',L1,L1,NA-NUMFRZ,ONE,WRK2,L1,
     *             VEC(1,NUMFRZ+1),L1,ZERO,ZTMP,L1)
      END IF
C
      RETURN
C
      END SUBROUTINE DFTB_ZVEC
C
C-----------------------------------------------------------------------
C*MODULE DFTBGR    *DECK DFTB_MFRZ
C>
C>    @brief   FMO-DFTB/AFO analytic derivative
C>
C>    @details Solve the standard or projected Z-vector equations
C>             and compute Z*B and Z*G.
C>
C>    @author  Yoshio Nishimoto
C>    - Feb, 2018- Subroutine written
C>
C>           --- INPUT ---
C>    @param DQES     Charges induced by relaxed density (Z)
C>    @param SHIFTES  Shift contributions from DQES
C>    @param WRK      Work space
C>    @param WRK1     Work space
C>    @param EIG      Orbital energy
C>    @param VEC      MO coefficient
C>    @param DIST     Distance between atom A and B
C>    @param GAMMA    gamma for DFTB2
C>    @param GAMMA3   Gamma for DFTB3
C>    @param S        Overlap matrix in AO
C>    @param DQ       Mulliken charge (not population)
C>           --- INPUT ---
C>    @param D        Density matrix in AO (on entry)
C>                    Density + relaxed density matrix in AO (on exit)
C>    @param W        Energy-weighted density matrix in AO (on entry)
C>
C
      SUBROUTINE DFTB_MFRZ(D,W,DQES,SHIFTES,
     *                     WRK,WRK1,EIG,VEC,DIST,GAMMA,GAMMA3,S,DQ)
      use mx_limits, only: mxatm,mxao

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      LOGICAL DFTBFL,SCC,SRSCC,DFTB3,DAMPXH,LCDFTB
      LOGICAL MFRZ,GOPARR,DSKWRK,MASWRK
C
      PARAMETER (ZERO=0.0D+00,ONE=1.0D+00,MXSPE=10)
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
      COMMON /FMOAFO/ LINDFRZ,LGMK,LGMKSAV
      COMMON /FMODTB/ LGAMMAGLOB,LZREF,LISPE,LSHIFTG,LSHIFTCT,
     *                LGAMMAD,MODESD,MODGAMMA,LDCTMP,LGAMMA3,LGAMMA3D,
     *                LNBOND,LCTMUL,LESPDFTB,LESPASC,NDFTBRST,LDFTBRST
      COMMON /FMOINF/ NFG,NLAYER,NATFMO,NBDFG,NAOTYP,NBODY,NSEGM
      COMMON /FMOOPT/ ESPSCA(9),RESPAP(2),RESPPC(2),RESDIM,RESTRI(4),
     *                RCORSD,RESPCT,CONVFG,CNVDMP,COROFF,RFLMO(4),
     *                ORSHFT,ORSHFT2,CNVAFO,ASCREEN(4),IXESP,MXITFG,
     *                NGUESS,NBSSE,MODORB,MODPAR,IRSTSTP,IRSTLAY,NPRFMO,
     *                NFMOPAL,MODPRP,MAXL1C,IPIEDA,MODGRD,MODESP,IVMUL,
     *                MODLMO,NOPDEN,MOFOCK,MODFD,MODFMM,NCENTM,ndualb,
     *                ngab,modpan
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON /MFRPAR/ MFRZ,NUMFRZ,NORFRZ,IFRZ(MXAO)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /PRPOPT/ ETOLLZ,ILOCAL,IAHARD
C
      DIMENSION D(*),W(*),DQES(*),SHIFTES(*),
     *          WRK(NUM),WRK1(NUM,*),EIG(*),VEC(NUM,NUM),DIST(*),
     *          GAMMA(*),GAMMA3(NAT,NAT),S(*),DQ(*)
      LOGICAL DEPLMO
C
      DEPLMO = .TRUE.
      IF (RFLMO(1).EQ.0) DEPLMO = .FALSE.
C
      L1 = NUM
      L2 = L1*(L1+1)/2
      L3 = L1*L1
      NORBPRJ = NORFRZ - NUMFRZ
      IF (NORBPRJ.LT.0) NORBPRJ=0
      NNUM = (L1-NA)*(NA-NUMFRZ)
      NFRZ = MAX(NUMFRZ,NORFRZ)
      CALL VCLR(DQES,1,NAT)
      CALL VCLR(SHIFTES,1,NAT)
C
      IF (NFRZ.EQ.0) THEN
        CALL VCLR(WRK1,1,L3)
        CALL DAWRIT(IDAF,IODA,D,L2,308,0)
        CALL DAWRIT(IDAF,IODA,WRK1,L3,474,0)
        CALL DAWRIT(IDAF,IODA,WRK1,L3,476,0)
        CALL DAWRIT(IDAF,IODA,WRK1,L3,477,0)
        IF (NFG.NE.0) CALL VCLR(X(LGMK),1,L1*NFRZ)
        RETURN
      END IF
C
      CALL VALFM(LOADFM)
      LEIGL = LOADFM + 1
      LFAO  = LEIGL  + L2
      LWRK1 = LFAO   + L2
      LWRK2 = LWRK1  + L3
      LWRK3 = LWRK2  + L3
      LWRK4 = LWRK3  + L3
      LZTMP = LWRK4  + L3
      LVECL = LZTMP  + L3
      LAST  = LVECL  + L3
C
      LR0   = LAST
      LZ0   = LR0    + MAX(NA*NA,NNUM)
      LP0   = LZ0    + MAX(NA*NA,NNUM)
      LAP   = LP0    + MAX(NA*NA,NNUM)
      LZ1   = LAP    + MAX(NA*NA,NNUM)
      LPREC = LZ1    + MAX(NA*NA,NNUM)
      LAST  = LPREC  + MAX(NA*NA,NNUM)
C
      IF (DEPLMO) THEN
        LRROT     = LAST
        LUROT     = LRROT     + NFRZ*NFRZ
        LLSC      = LUROT     + NFRZ*NFRZ
        LAST      = LLSC      + L3
C
        !! For non-AFO calculations.
        !! For AFO, it has been allocated somewhere
        IF (RFLMO(1).EQ.0) THEN
          LGMK    = LAST
          LAST    = LGMK      + L3
        END IF
      END IF
C
      IF (DFTB3) THEN
        LSHIFT3   = LAST
        LAST      = LSHIFT3   + NAT*2
      END IF
C
      NEED  = LAST - LOADFM - 1
      CALL GETFM(NEED)
C
      CALL DAREAD(IDAF,IODA,EIG,L1,17,0)
      CALL DAREAD(IDAF,IODA,S,L2,12,0)
      CALL DAREAD(IDAF,IODA,X(LVECL),L3,318,0)
C
      CALL DAWRIT(IDAF,IODA,D,L2,308,0)
C
      CALL DFTB_GAMMA(DAMPXHE,HUBBL,GAMMA,DIST,NAT,NSPE,
     *  ISPE,MAXANG,SPE,DAMPXH,SRSCC)
      IF (DFTB3) THEN
        CALL VCLR(GAMMA3,1,NAT*NAT)
        CALL DFTB_GAMMA3(NAT,HUBBL,GAMMA3,DIST)
        IF (NFG.NE.0) THEN
          CALL DAREAD(IDAF,IODA,X(LSHIFT3),NAT*2,567,0)
        ELSE
          CALL VCLR(X(LSHIFT3),1,NAT*2)
        END IF
      END IF
C
      IF (DEPLMO) THEN
        CALL VCLR(X(LGMK),1,L1*NFRZ)
        CALL DAREAD(IDAF,IODA,X(LUROT),NFRZ*NFRZ,569,0)
        CALL CPYTSQ(S,X(LWRK2),L1,1)
        CALL DGEMM('T','N',L1,L1,L1,ONE,X(LVECL),L1,
     *             X(LWRK2),L1,ZERO,WRK1,L1)
        CALL DAREAD(IDAF,IODA,X(LWRK1),L3,15,0)
        IF (NORBPRJ.GT.0) THEN
          CALL DCOPY(L1*(L1-NUMFRZ),X(LWRK1+L1*NUMFRZ),1,X(LWRK2),1)
          CALL DCOPY(L1*NORBPRJ,X(LWRK1+L1*(L1-NORBPRJ)),1,
     *               X(LWRK1+L1*NUMFRZ),1)
          CALL DCOPY(L1*(L1-NUMFRZ),X(LWRK2),1,X(LWRK1+L1*NORFRZ),1)
        END IF
        CALL DGEMM('N','N',L1,L1,L1,ONE,WRK1,L1,
     *             X(LWRK1),L1,ZERO,X(LLSC),L1)
        CALL DGEMM('N','N',NFRZ,NFRZ,NFRZ,ONE,X(LUROT),NFRZ,
     *             X(LLSC),L1,ZERO,X(LRROT),NFRZ)
      END IF
      IF (NORBPRJ.GT.0) THEN
        CALL DCOPY(L1*NORBPRJ,X(LVECL+L1*NUMFRZ),1,
     *             X(LVECL+L1*(L1-NORBPRJ)),1)
      END IF
C
C
C
      IF (NORBPRJ.EQ.0) THEN
        CALL DAREAD(IDAF,IODA,X(LFAO),L2,317,0)
        GO TO 1000
      END IF
C
C
C
      !! CONSTRUCT SL
      CALL CPYTSQ(S,X(LWRK1),L1,1)
      CALL DGEMM('N','N',L1,L1,L1,ONE,X(LWRK1),L1,X(LVECL),L1,ZERO,
     *           X(LWRK2),L1)
      !! X(LWRK2) CONTAINS SL; IT SHOULD NOT BE OVERWRITTEN
      !! UNTIL Z_{\mu \nu}^\prime IS OBTAINED
C
C     ----- SOLVE ZA=B -----
C
C     (A) CONSTRUCT INITIAL RESIDUE R0
      CALL DAREAD(IDAF,IODA,X(LFAO),L2,317,0)
      CALL DAREAD(IDAF,IODA,VEC,L3,15,0)
      CALL CPYTSQ(X(LFAO),X(LWRK1),L1,1)
      CALL DGEMM('T','N',L1-NA,L1,L1,ONE,VEC(1,NA+1),L1,X(LWRK1),L1,
     *           ZERO,WRK1,L1-NA)
      CALL DGEMM('N','N',L1-NA,NA-NUMFRZ,L1,ONE,WRK1,L1-NA,
     *           VEC(1,NUMFRZ+1),L1,ZERO,X(LR0),L1-NA)
      CALL DSCAL(NNUM,4.0D+00,X(LR0),1)
C     CALL PRSQ(X(LR0),NA-NUMFRZ,L1-NA,L1-NA)
C
C     (B) CONSTRUCT Z0 FROM R0 AND MAKE P0=Z0
C     (C) PREPARE PRECONDITIONER
      NSEQ = 0
      DO JO = NUMFRZ+1, NA
        DO IO = NA+1, L1
          NSEQ = NSEQ + 1
          X(LPREC+NSEQ-1) = ONE/(EIG(JO)-EIG(IO))
          X(LZ0+NSEQ-1) = X(LR0+NSEQ-1)*X(LPREC+NSEQ-1)
        END DO
      END DO
C     CALL PRSQ(X(LZ0),NA-NUMFRZ,L1-NA,L1-NA)
C
C     ----- START ITERATION -----
C
      IF (IAND(MODGRD,128).EQ.0) THEN
        CALL DFTB_ZVEC(1,L1,L2,NNUM,NA,NAT,NUMFRZ,X(LR0),
     *                 X(LZ0),X(LP0),X(LAP),X(LZ1),X(LPREC),X(LZTMP),
     *                 WRK,WRK1,X(LWRK1),X(LWRK2),
     *                 VEC,X(LFAO),S,EIG,GAMMA,GAMMA3,DQ,
     *                 NSPE,ISPE,IND,MAXANG,IW,MASWRK,X(LVECL),NORBPRJ,
     *                 IAND(NPRFMO,3).NE.3,DFTB3,X(LSHIFT3),MODGAMMA)
      ELSE
        CALL DCOPY(NNUM,X(LZ0),1,X(LZ1),1)
        CALL DGEMM('N','N',L1,NA-NUMFRZ,L1-NA,ONE,VEC(1,NA+1),L1,
     *             X(LZ1),L1-NA,ZERO,X(LWRK2),L1)
        CALL DGEMM('N','T',L1,L1,NA-NUMFRZ,ONE,X(LWRK2),L1,
     *             VEC(1,NUMFRZ+1),L1,ZERO,X(LZTMP),L1)
      END IF
C
C     ----- COMPUTE B*Z -----
C
      ! D appears twice, but no problem
      ! D and W are updated
      CALL DFTB_ZVEC_ZB(0,L1,NNUM,NFRZ,X(LZ1),D,W,DQ,DQES,SHIFTES,
     *                  X(LSHIFT3),EIG,S,D,X(LFAO),X(LEIGL),VEC,
     *                  X(LVECL),X(LUROT),X(LRROT),X(LGMK),X(LZTMP),
     *                  X(LWRK2),GAMMA,GAMMA3,X(LGAMMA3),WRK,WRK1,
     *                  X(LWRK1),X(LWRK3),X(LWRK4),DEPLMO,DFTB3,
     *                  MODGAMMA,X(LAST),X(LAST))
C     ! On exit, X(LFAO) should have the pristine Fock matrix, read from
C     ! 317, in the current dimer/trimer
C
 1000 CONTINUE
      IF (.NOT.DEPLMO) GO TO 200
C
C     INCLUDE GEOMETRY-DEPENDENCE (RESPONSE IN MODEL SYSTEMS)
C
      !! compute G_{\mu k} from non-CP contributions
      IF (NORBPRJ.EQ.0) CALL DAREAD(IDAF,IODA,VEC,L3,15,0)
      IF (NUMFRZ.GT.0) THEN
        CALL DFTB_CNSTGMK(L1,NA,NUMFRZ,NFRZ,ONE,
     *                   X(LGMK),S,X(LFAO),VEC,X(LRROT),X(LWRK1),
     *                   X(LWRK2),WRK1)
      END IF
C     WRITE (*,*) "G_{\MU K}"
C     CALL PRSQ(X(LGMK),NFRZ,L1,L1)
C
      CALL VCLR(X(LWRK1),1,L3)
      CALL DAWRIT(IDAF,IODA,X(LWRK1),L3,474,0) !! (X+Y) is always zero
      CALL CPYTSQ(W,WRK1,L1,1)
      CALL DAWRIT(IDAF,IODA,WRK1,L3,477,0)
      IF (RFLMO(1).NE.0) THEN
        CALL RETFM(NEED)
        RETURN
      END IF
C
C     ----- NON-AFO CALCULATIONS -----
C
      !! some projection of G_{\mu k} is needed
      CALL DAREAD(IDAF,IODA,X(LWRK2),NFRZ*NFRZ,569,0)
      CALL DTRTRI('U','N',NFRZ,X(LWRK2),NFRZ,INFO)
      CALL DAREAD(IDAF,IODA,X(LZTMP),L3,318,0)
      CALL DGEMM('N','N',L1,NFRZ,NFRZ,ONE,X(LZTMP),L1,
     *           X(LWRK2),NFRZ,ZERO,WRK1,L1)
      DO KO = 1, NFRZ
        DO MU = 1, L1
          IF (ABS(WRK1(MU,KO)).LE.1.0D-12) X(LGMK+MU-1+L1*(KO-1))=ZERO
        END DO
      END DO
      !! here is manual manupulation
      call dcopy(l3,x(LGMK),1,x(lwrk2),1)
      call vclr(x(LGMK),1,l3)
      call dcopy(l1*3,x(lwrk2),1,x(lwrk1),1)
      call dscal(l1,-one,x(lwrk1+l1),1)
C     write (*,*) "G_{\mu k} for model system"
C     call prsq(x(lwrk1),l1,l1,l1)
      CALL DCOPY(L1*NA,X(LWRK1),1,X(LGMK),1)
      CALL DCOPY(L3,WRK1,1,X(LWRK2),1)
C
      CALL DAREAD(IDAF,IODA,S,L2,12,0)
      CALL DAREAD(IDAF,IODA,VEC,L3,15,0)
      CALL DAREAD(IDAF,IODA,EIG,L1,17,0)
      CALL DFTB_LOCZVEC(0,X(LGMK),.TRUE.)
      CALL DAWRIT(IDAF,IODA,VEC,L3,15,0)
C
C     ---- FINALIZING... -----
C
      CALL DAREAD(IDAF,IODA,X(LWRK1),L3,318,0)
      CALL DAWRIT(IDAF,IODA,X(LWRK1),L3,15,0)
      !!!!!
      CALL DFTBX(DUM,DUM)
      !!!!!
  200 CONTINUE
      CALL DAREAD(IDAF,IODA,S,L2,12,0)
      CALL DAREAD(IDAF,IODA,VEC,L3,15,0)
      CALL DAREAD(IDAF,IODA,EIG,L1,17,0)
C
      ILOCAL = 0
C
      CALL RETFM(NEED)
C
      END SUBROUTINE DFTB_MFRZ
C
C-----------------------------------------------------------------------
C*MODULE DFTBGR    *DECK  LCGAMMADER
C>
C>    @brief Gradient of DFTB
C>
C>    @details Calculate lc-gamma 1st-order derivative
C>
C>    @author Van Quan Vuong 
C>    - Jan, 2018- Subroutine written
C>
C>           --- INPUT ---
C>    @param EMU range-separation parameter for LC-DFTB
C>    @param HUBB Hubbard values of atoms
C>    @param NAT Number of atoms
C>    @param ISPE Index of species
C>    @param DISTMAT Matrix of inter-atomic distances
C>
C>           --- OUTPUT ---
C>    @param GAMMA LC-Gamma used in LC-DFTB
C>    @param GAMMADER LC-Gamma 1st-order derivative used in LC-DFTB
C>
C
      SUBROUTINE DFTB_LCGAMMADER(EMU,HUBB,NAT,ISPE,DISTMAT,
     *                           GAMMA,GAMMADER)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      INTEGER :: NAT,ISPE(NAT)
      DOUBLE PRECISION :: EMU,HUBB(*),DISTMAT(*)
      DOUBLE PRECISION :: GAMMA(*),GAMMADER(*)

      DOUBLE PRECISION irij,ita2EMU2,itb2EMU2,ita2tb2,itb2ta2,
     *                 i48tau5,irij2,it2EMU24,i48

      DOUBLE PRECISION, PARAMETER  :: TOLSAMEDIST = 1.0D-05
      DOUBLE PRECISION, PARAMETER  :: TOLSAMEHUBB = 0.3125D-05

C     WRITE (*,*) "CALCULATE GAMMA VALUES FOR LC-DFTB"
C     WRITE (*,*) "OMEGA = ", EMU 

      nseq = 0
      DO I = 1, NAT
C       WRITE(*,*)  HUBB(ISPE(I))
        DO J = 1, I
          ta = HUBB(ISPE(I)) * 3.2D+00
          tb = HUBB(ISPE(J)) * 3.2D+00
          nseq = nseq + 1
          rij = DISTMAT(nseq)
          EMU2 = EMU *  EMU
          EMU4 = EMU2 * EMU2
          EMU6 = EMU4 * EMU2

          IF (rij .lt. TOLSAMEDIST) THEN
            IF (ABS(ta - tb) .lt. TOLSAMEHUBB) THEN
              tau = 0.5D+00 * (ta + tb)
              tau2= tau * tau
              tau4= tau2 * tau2
              tmp = 5.0D+00 * tau4 *tau2 + 15.0D+00 * tau4 * EMU2 
     *              - 5.0D+00 * tau2 * EMU4  + EMU6
              tmp = tmp * 0.0625D+00 / (tau4*tau) - EMU
              tmp = tmp * tau4*tau4 / (tau2 - EMU2)**4
              GAMMA(nseq) = tau * 0.3125D+00 - tmp
              GAMMADER(nseq) = 0.0D+00
            ELSE
              write(6,'(a)') "Error(RangeSep): R = 0, Ua != Ub"
              call abrt ! or simply gamma(nseq)=zero?
            END IF
          ELSE
            irij = 1.0D+00 / rij
            irij2= irij * irij
            IF (ABS(ta - tb) .lt. TOLSAMEHUBB ) THEN
              i48 = 1.0D+00 / 48.0D+00
              tau = 0.5D+00 * (ta + tb)
              tau2 = tau * tau
              tau4 = tau2 * tau2
              tau6 = tau4 * tau2
              tau8 = tau4 * tau4
              taurij = tau * rij
              it2EMU24 = 1.0D+00 / (tau2 - EMU2)**4
               
              expmtaurij = EXP(-taurij)
              expmWrij = EXP(-EMU * rij)
              rij2 = rij * rij

              i48tau5 = i48 / (tau4*tau)
              tcons1  = (3.0D+00 * tau4 * EMU4 - 3.0D+00
     *                  * tau6 * EMU2 - tau2 * EMU6)*i48tau5
              tcons2  = (15.0D+00 * tau2*tau * EMU4 - 21.0D+00
     *                  * tau4*tau * EMU2 - 3.0D+00 * tau*EMU6)
     *                  *i48tau5
              tcons3  = (15.0D+00 * tau2 * EMU4 - 45.0D+00 * tau4
     *                  * EMU2 - 3.0D+00 * EMU6)*i48tau5 

              tmp0 = (taurij**3 * i48 + 0.1875D+00 * taurij**2
     *               + 0.6875D+00 * taurij + 1.0D+00) * expmtaurij
     *               * irij

              tmp = -tau8 * it2EMU24 * (tmp0 + expmtaurij 
     *              * (rij2 * tcons1 + rij * tcons2 + tcons3))        

              GAMMA(nseq) = irij - tmp0 - (tau8 * it2EMU24
     *                      * expmWrij * irij + tmp)


              td0     = tau2*taurij * 2.0D+00 * i48 + 0.1875D+00 * tau2 
     *                  - irij2

              GAMMADER(nseq) = - irij2 + tau * tmp0 - expmtaurij * td0
     *                         - (tau8 * it2EMU24
     *                         * expmWrij * irij * (-irij-EMU) 
     *                         - tau * tmp - tau8 * it2EMU24 
     *                         * expmtaurij * (td0 + 2.0D+00 * rij 
     *                         * tcons1 + tcons2))

            ELSE

              ta2 = ta * ta
              tb2 = tb * tb
              ta4 = ta2 * ta2
              tb4 = tb2 * tb2
  
              ita2tb2 = 1.0D+00 / (ta2 - tb2)
              itb2ta2 = -ita2tb2

              ita2EMU2 = 1.0D+00 / (ta2 - EMU2)
              itb2EMU2 = 1.0D+00 / (tb2 - EMU2)
             
              expmtarij=EXP(-ta * rij)
              expmtbrij=EXP(-tb * rij)
              expmWrij = EXP(-EMU * rij)            

              tconsa1 = ta2 * ita2EMU2 * ta*tb4 * itb2ta2**2 * 0.5D+00
              tconsa2 = ta4 * ita2EMU2**2 * (tb4 * tb2 - 3.0D+00 
     *                  * ta2 * tb4 + 2.0D+00 * EMU2 * tb4) 
     *                  * ita2tb2**3
              tconsb1 = tb2 * itb2EMU2 * tb*ta4 * ita2tb2**2 * 0.5D+00
              tconsb2 = tb4 * itb2EMU2**2 * (ta4 * ta2 - 3.0D+00 
     *                  * tb2 * ta4 + 2.0D+00 * EMU2 * ta4) 
     *                  * itb2ta2**3

 
              tmp = expmtarij * (tconsa1 - tconsa2 * irij)
     *              + expmtbrij * (tconsb1 - tconsb2 * irij)


              tconsa01 = ta*tb4 * itb2ta2**2 * 0.5D+00 
              tconsa02 = (tb4 * tb2 - 3.0D+00 * ta2 * tb4) * ita2tb2**3
              tconsb01 = tb*ta4 * ita2tb2**2 * 0.5D+00
              tconsb02 = (ta4 * ta2 - 3.0D+00 * tb2 * ta4) * itb2ta2**3

              tmp0     = expmtarij * (tconsa01 -  tconsa02 * irij)
     *                   + expmtbrij * (tconsb01 - tconsb02 * irij)
              
              GAMMA(nseq) = irij - tmp0 - (ta4 * tb4 * ita2EMU2**2 
     *                      * itb2EMU2**2 * expmWrij * irij
     *                      - tmp) 

              tmp0    = expmtarij * (-ta * tconsa01
     *                  + (ta + irij) * tconsa02 * irij)
     *                  + expmtbrij * (-tb * tconsb01 
     *                  + (tb + irij) * tconsb02 * irij)


              tmp     = expmtarij * (-ta * tconsa1
     *                  + (ta + irij) * tconsa2 * irij)
     *                  + expmtbrij * (-tb * tconsb1 
     *                  + (tb + irij) * tconsb2 * irij)


              GAMMADER(nseq) = - irij2 - tmp0  - (ta4 * tb4 
     *                         * ita2EMU2**2 * itb2EMU2**2 
     *                         * expmWrij * irij2 * ( -1.0D+00 
     *                         - EMU * rij ) - tmp )
            END IF
          END IF 
        END DO
      END DO
      END SUBROUTINE DFTB_LCGAMMADER
C
C-----------------------------------------------------------------------
C*MODULE DFTBGR    *DECK  DFTB_LCDER
C>
C>    @brief Gradient of DFTB
C>
C>    @details Calculate long-range contribution of gradient for LC-DFTB
C>
C>    @author Van Quan Vuong 
C>    - Jan, 2018- Subroutine written
C>
C>           --- INPUT ---
C>    @param NAT Number of atoms
C>    @param L1 Number of basis functions
C>    @param XDIST Matrix of inter-atomic distances
C>    @param XDSDR Matrix of AO overlap 1st-order derivative dS/dR 
C>    @param XDGDR Matrix of gamma 1st-order derivative dgamma/dR 
C>    @param XSS AO overlap matrix
C>    @param XDDS \Delta Density matrix
C>    @param XLCGAMMAS Matrix of LC-gamma
C>    @param XSUMSA Matrix of sum of AO overlap for atoms
C>    @param XSUMSS Matrix of sum of AO overlap for atomic pairs
C>    @param XMAXDDA Maximum of \Delta density matrix for atoms
C>    @param XMAXDDS Maximum of \Delta density matrix for atomic pairs
C>    @param XSUMDSDR Sum of dS/dR for atoms
C>    @param SUMS Sum of AO overlap
C>    @param MAXDD Maximum of \Delta density matrix 
C>    @param TTHRESLCDFTB Switch for threshold algorithm of LC-DFTB
C>    @param THRESLCDFTB Cut-off for threshold algorithm of LC-DFTB 
C>    @param SKCUT2 Cut-off distance both Hamiltonian and overlap
C>
C>           --- OUTPUT ---
C>    @param EGRAD   Gradient
C>
C
      SUBROUTINE DFTB_LCDER(NAT,L1,XDIST,XDSDR,XDGDR,XSS,XDDS,XLCGAMMAS,
     *           XSUMSA,XSUMSS,XMAXDDA,XMAXDDS,XSUMDSDR,EGRAD,
     *           SUMS,MAXDD,TTHRESLCDFTB,ICUT,SKCUT2)
      use mx_limits, only: mxatm

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      INTEGER,PARAMETER :: MXSPE=10
      DOUBLE PRECISION,PARAMETER::ZERO=0.0D+00

      COMMON /DFTBPR/ ETEMP,DFTBDP(MXSPE*14),DAMPXHE,HUBDER(MXSPE),
     *                ZREF(MXATM),SPNCST(6,MXSPE),SPE(MXATM),NSPE,
     *                MAXANG(MXATM),ISPE(MXATM),IND(MXATM+1),IDFTBD,
     *                PARAMDIR
      DOUBLE PRECISION :: ETEMP,DFTBDP,DAMPXHE,HUBDER,ZREF,SPNCST,SPE
      INTEGER :: NSPE,MAXANG,ISPE,IND,IDFTBD

      DOUBLE PRECISION :: XDIST(*),XDSDR(*),XDGDR(*),XSS(L1,L1),       
     *                    XDDS(L1,L1),XLCGAMMAS(NAT,NAT),XSUMSA(*),
     *                    XSUMSS(NAT,NAT),XMAXDDA(*),XMAXDDS(NAT,NAT),
     *                    XSUMDSDR(*),EGRAD(3,MXATM)

      LOGICAL TTHRESLCDFTB
      INTEGER :: NAT,L1,ICUT
      DOUBLE PRECISION :: SUMS,MAXDD,THRESLCDFTB,SKCUT2

C     CALL CPU_TIME(T1)

      NAT2=NAT*(NAT+1)/2
      NAT3=NAT*NAT
      L3 = L1*L1
      IF (TTHRESLCDFTB) THEN
        THRESLCDFTB  = 1.0D+01**(-ICUT)
      ELSE
        THRESLCDFTB = ZERO
      END IF

      DO IA = 1, NAT                                               !! atom A
        DO IB = IA + 1, NAT                                             !! atom B
          IF (IA.EQ.IB) CYCLE
          DO K = 1, 3
            CALL DFTB_CNVSQ(IA,IB,NSEQAB)
            DIST = XDIST(NSEQAB)
            IF (DIST*DIST.GE.SKCUT2) THEN
              CYCLE
            END IF

            DGDR = XDGDR(NAT3*(K-1)+(IB-1)*NAT+IA)
            NAOA = MAXANG(ISPE(IA))**2
            NAOB = MAXANG(ISPE(IB))**2

            TMPS = ZERO 
            TMPG1= ZERO 
            TMPG2= ZERO 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C                 DSDR CONTRIBUTION

            IF(TTHRESLCDFTB) THEN
              TDELTA1 = XSUMDSDR(NAT2*(K-1)+NSEQAB)
              IF ((XMAXDDA(IA)*XMAXDDA(IB)
     *            * TDELTA1*SUMS).LT.THRESLCDFTB) GO TO 111 
            END IF
            DO IC = 1, NAT
              IF(TTHRESLCDFTB) THEN
                TDELTA2 = TDELTA1 * XMAXDDS(IA,IC)
                IF ((XMAXDDA(IB)
     *              *TDELTA2*XSUMSA(IC)).LT.THRESLCDFTB) CYCLE
              END IF
              NAOC = MAXANG(ISPE(IC))**2
              GAMACBC = XLCGAMMAS(IA,IC) + XLCGAMMAS(IB,IC)
              DO ID = 1, NAT
                IF(TTHRESLCDFTB) THEN
                  TDELTA3 = TDELTA2 * XMAXDDS(IB,ID)
                  IF((TDELTA3*XSUMSS(IC,ID)).LT.THRESLCDFTB)CYCLE
                END IF
                NAOD = MAXANG(ISPE(ID))**2
                TMP0 = ZERO 
                DO IAT = 1, NAOA           
                  IOA = IND(IA) + IAT                        !! AO atom A
                  DO IBT = 1, NAOB           
                    IOB = IND(IB) + IBT                        !! AO atom B
                    TMP1 = ZERO 
                    DO ICT = 1, NAOC           
                      IOC = IND(IC) + ICT                        !! AO atom C
                      DDOAOC = XDDS(IOA,IOC)
                      DO IDT = 1, NAOD 
                        IOD = IND(ID) + IDT                      !! AO atom D
                        TMP1 = TMP1 + DDOAOC*XSS(IOC,IOD)*XDDS(IOD,IOB)
                      END DO 
                    END DO 
                    TMP0 = TMP0 - TMP1 * XDSDR(L3*(K-1)+(IOB-1)*L1+IOA)
                  END DO
                END DO
                TMPS = TMPS + TMP0 *(GAMACBC
     *                             +XLCGAMMAS(IA,ID) + XLCGAMMAS(IB,ID))
              END DO
            END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C                 DGDR1 CONTRIBUTION

  111       IF(TTHRESLCDFTB) THEN
              TDELTA1 = ABS(DGDR)
              IF ((TDELTA1*XMAXDDA(IA)*XMAXDDA(IB)
     *             *XSUMSA(IA)*XSUMSA(IB)).LT.THRESLCDFTB) 
     *             GO TO 222 
            END IF
            DO IC = 1, NAT
              IF(TTHRESLCDFTB) THEN
                TDELTA2 = TDELTA1 * XMAXDDS(IA,IC)
     *                    *XSUMSS(IB,IC)
                IF ((TDELTA2*XSUMSA(IA)*XMAXDDA(IB))
     *               .LT.THRESLCDFTB) CYCLE 
              END IF
              NAOC = MAXANG(ISPE(IC))**2
              DO ID = 1, NAT
                IF(TTHRESLCDFTB) THEN
                  TDELTA3 = TDELTA2*XMAXDDS(IB,ID) 
                  IF (TDELTA3 
     *             *XSUMSS(IA,ID).LT.THRESLCDFTB) CYCLE 
                END IF
                NAOD = MAXANG(ISPE(ID))**2
                DO IAT = 1, NAOA           
                  IOA = IND(IA) + IAT                        !! AO atom A
                  DO IBT = 1, NAOB           
                    IOB = IND(IB) + IBT                        !! AO atom B
                    DO ICT = 1, NAOC           
                      TMP0 = ZERO 
                      IOC = IND(IC) + ICT                        !! AO atom C
                      DO IDT = 1, NAOD           
                        IOD = IND(ID) + IDT                      !! AO atom D
                        TMP0 = TMP0  + XSS(IOA,IOD)*XDDS(IOB,IOD)
                      END DO 
                      TMPG1 = TMPG1 + TMP0*XSS(IOB,IOC)*XDDS(IOA,IOC)
                    END DO 
                  END DO
                END DO
              END DO
            END DO
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C                 DGDR2 CONTRIBUTION
 
  222       IF(TTHRESLCDFTB) THEN
              TDELTA1 = ABS(DGDR)*XMAXDDS(IA,IB)
              IF ((TDELTA1*MAXDD*XSUMSA(IA)
     *           *XSUMSA(IB)).LT.THRESLCDFTB) GO TO 333 
            END IF
            DO IC = 1, NAT
              IF(TTHRESLCDFTB) THEN
                TDELTA2 =  TDELTA1*XSUMSS(IB,IC)
                IF ((TDELTA2*XMAXDDA(IC)
     *              *XSUMSA(IA)).LT.THRESLCDFTB) CYCLE 
              END IF
              NAOC = MAXANG(ISPE(IC))**2
              DO ID = 1, NAT
                IF(TTHRESLCDFTB) THEN
                  TDELTA3 =  TDELTA2 * XMAXDDS(IC,ID)
                  IF ((TDELTA3
     *             *XSUMSS(IA,ID)).LT.THRESLCDFTB) CYCLE
                END IF
                NAOD = MAXANG(ISPE(ID))**2
                DO IAT = 1, NAOA           
                  IOA = IND(IA) + IAT                        !! AO atom A
                  DO IBT = 1, NAOB           
                    IOB = IND(IB) + IBT                        !! AO atom B
                    TMP0 = ZERO 
                    DO ICT = 1, NAOC 
                      IOC = IND(IC) + ICT                        !! AO atom C
                      TMP1 = ZERO 
                      DO IDT = 1, NAOD         
                        IOD = IND(ID) + IDT                      !! AO atom D
                        TMP1 = TMP1 +  XDDS(IOC,IOD) * XSS(IOA,IOD)
                      END DO 
                      TMP0 = TMP0 + TMP1 * XSS(IOB,IOC)
                    END DO 
                    TMPG2= TMPG2 + TMP0 * XDDS(IOA,IOB) 
                  END DO
                END DO
              END DO
            END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  333       EGRAD(K,IA) = EGRAD(K,IA) - TMPS - (TMPG1+TMPG2)*DGDR 
            EGRAD(K,IB) = EGRAD(K,IB) + TMPS + (TMPG1+TMPG2)*DGDR 
          END DO
        END DO
      END DO

C     CALL CPU_TIME(T2)
C     WRITE (*,'(" CPU TIME IN LCDER ",F8.4)') T2-T1

      END SUBROUTINE DFTB_LCDER
C
C-----------------------------------------------------------------------
C*MODULE DFTBGR    *DECK  DFTB_LCGRAD1
C>
C>    @brief   Gradient of LC-DFTB
C>
C>    @details Calculate pre-long-range contributions of gradient which
C>             comes from the derivative of overlap matrix. Actual
C>             contributions are computed elsewhere, when overlap
C>             derivatives are computed. After standardSCF, DD1 and DD2
C>             should be identical. After solving response equations,
C>             DD1 and DD2 are usually different.
C>             This subroutine does not use the screening threshold.
C>
C>    @author Yoshio Nishimoto
C>    - May, 2018- Subroutine written
C>
C>           --- INPUT ---
C>    @param MODE    MODE=0: DD1 = DD2; MODE=1: DD1 \= DD2
C>    @param NAT     Number of atoms
C>    @param L1      Number of AOs
C>    @param ISPE    Index of species
C>    @param MAXANG  Maximum angular momentum of species
C>    @param IND     AO index of each atom
C>    @param DD1     Input density matrix 1
C>    @param DD2     Input density matrix 2
C>    @param S       Overlap matrix in AO
C>    @param GAMMA   Long-range gamma
C>           --- OUTPUT ---
C>    @param FMOLC   Pre contributions
C>
      SUBROUTINE DFTB_LCGRAD1(MODE,NAT,L1,ISPE,MAXANG,IND,DD1,DD2,S,
     *                        GAMMA,FAOLC)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER (ZERO=0.0D+00,ONE=1.0D+00)
C
      COMMON /FMCOM / X(1)
C
      DIMENSION ISPE(*),MAXANG(*),IND(*),DD1(*),DD2(*),S(*),GAMMA(*),
     *          FAOLC(L1*L1)
C
      L3=L1*L1
      CALL VALFM(LOADFM)
      CALL GOTFM(NGOTMX)
      LSS     = LOADFM  + 1
      LWRK1   = LSS     + L3
      LWRK2   = LWRK1   + L3
      LAST    = LWRK2   + L3
      NEED    = LAST - LOADFM - 1
      CALL GETFM(NEED)
C
      CALL VCLR(FAOLC,1,L3)
      DO IMODE = 1, 1+MODE
        CALL CPYTSQ(S,X(LSS),L1,1)
        IF (IMODE.EQ.1) THEN
          CALL DGEMM('N','N',L1,L1,L1,ONE,DD1,L1,X(LSS),L1,
     *               ZERO,X(LWRK1),L1)
        ELSE
          CALL DGEMM('N','T',L1,L1,L1,ONE,X(LSS),L1,DD2,L1,
     *               ZERO,X(LWRK1),L1)
        END IF
        CALL DCOPY(L3,X(LWRK1),1,X(LWRK2),1)
C
        DO IAT = 1, NAT
          DO JAT = 1, NAT
            CALL DFTB_CNVSQ(IAT,JAT,NSEQ)
            GAMTMP = GAMMA(NSEQ)
            DO MU = IND(IAT)+1, IND(IAT)+MAXANG(ISPE(IAT))**2
              DO NU = IND(JAT)+1, IND(JAT)+MAXANG(ISPE(JAT))**2
                X(LSS+MU-1+L1*(NU-1)) = GAMTMP
              END DO
            END DO
          END DO
        END DO
C
        DO I = 1, L3
          X(LWRK1+I-1) = X(LWRK1+I-1)*X(LSS+I-1)
        END DO
        IF (IMODE.EQ.1) THEN
          CALL DGEMM('N','T',L1,L1,L1,ONE,X(LWRK1),L1,DD2,L1,
     *               ONE,FAOLC,L1)
        ELSE
          CALL DGEMM('N','N',L1,L1,L1,ONE,DD1,L1,X(LWRK1),L1,
     *               ONE,FAOLC,L1)
        END IF
C
        IF (IMODE.EQ.1) THEN
          CALL DCOPY(L3,DD2,1,X(LWRK1),1)
        ELSE
          CALL DCOPY(L3,DD1,1,X(LWRK1),1)
        END IF
        DO I = 1, L3
          X(LWRK1+I-1) = X(LWRK1+I-1)*X(LSS+I-1)
        END DO
        IF (IMODE.EQ.1) THEN
          CALL DGEMM('N','T',L1,L1,L1,ONE,X(LWRK2),L1,X(LWRK1),L1,
     *               ONE,FAOLC,L1)
        ELSE
          CALL DGEMM('N','N',L1,L1,L1,ONE,X(LWRK1),L1,X(LWRK2),L1,
     *               ONE,FAOLC,L1)
        END IF
      END DO
C
      CALL DCOPY(L3,FAOLC,1,X(LWRK1),1)
      CALL SYMSQT(L1,X(LWRK1),L1,FAOLC,1)
      IF (MODE.EQ.0) CALL DSCAL(L1*(L1+1)/2,-0.25D+00,FAOLC,1)
      IF (MODE.EQ.1) CALL DSCAL(L1*(L1+1)/2,-0.125D+00,FAOLC,1)
C
      CALL RETFM(NEED)
C
      RETURN
C
      END SUBROUTINE DFTB_LCGRAD1
C
C-----------------------------------------------------------------------
C*MODULE DFTBGR    *DECK  DFTB_LCGRAD2
C>
C>    @brief   Gradient of LC-DFTB
C>
C>    @details Calculate long-range contributions which come from the
C>             derivative of long-range gamma functions. 
C>             This subroutine does not use the screening threshold.
C>
C>    @author Yoshio Nishimoto
C>    - May, 2018- Subroutine written
C>
C>           --- INPUT ---
C>    @param MODE     MODE=0: DD1 = DD2; MODE=1: DD1 \= DD2
C>    @param NAT      Number of atoms
C>    @param L1       Number of AOs
C>    @param ISPE     Index of species
C>    @param MAXANG   Maximum angular momentum of species
C>    @param IND      AO index of each atom
C>    @param DD1      Input density matrix 1
C>    @param DD2      Input density matrix 2
C>    @param SS       Overlap matrix in AO
C>    @param DISTMAT  Inter-atomic distances
C>    @param GAMDERLC Derivatives of LC-gamma wrt distance
C>    @param C        Coordinates
C>           --- IN/OUTPUT ---
C>    @param EGRAD    Gradient
C>
      SUBROUTINE DFTB_LCGRAD2(MODE,SCAL0,NAT,L1,ISPE,MAXANG,IND,DD1,DD2,
     *                        SS1,SS2,DISTMAT,GAMDERLC,C,EGRAD,HESSIAN)
C 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER (ZERO=0.0D+00,ONE=1.0D+00)
C
      COMMON /FMCOM / X(1)
C
      DIMENSION ISPE(*),MAXANG(*),IND(*),DD1(L1,L1),DD2(L1,L1),SS1(*),
     *          SS2(*),DISTMAT(*),GAMDERLC(*),C(3,*),EGRAD(3,*),
     *          HESSIAN(3*NAT,3*NAT)
C
      DOUBLE PRECISION VEC(3)
C
      L3 = L1*L1
      CALL VALFM(LOADFM)
      LWRK1   = LOADFM  + 1
      LWRK2   = LWRK1   + L3
      LWRK3   = LWRK2   + L3
      LWRK4   = LWRK3   + L3
      LAST    = LWRK4   + L3
      IF (MODE.EQ.1) THEN
        LWRK5 = LAST
        LWRK6 = LWRK5   + L3
        LAST  = LWRK6   + L3
      ELSE
        LWRK5 = LWRK3
        LWRK6 = LWRK4
      END IF
      NEED    = LAST - LOADFM - 1
      CALL GETFM(NEED)
C
      IF (MODE.LE.2) THEN
        !! For \gamma_{AB} and \gamma_{CD}
        CALL DGEMM('N','N',L1,L1,L1,ONE,SS1,L1,DD1,L1,
     *             ZERO,X(LWRK3),L1)
        CALL DGEMM('N','N',L1,L1,L1,ONE,X(LWRK3),L1,SS1,L1,
     *             ZERO,X(LWRK1),L1)
        IF (MODE.EQ.0.OR.MODE.EQ.2) THEN
          CALL DCOPY(L3,X(LWRK1),1,X(LWRK2),1)
        ELSE
          CALL DGEMM('N','N',L1,L1,L1,ONE,SS1,L1,DD2,L1,
     *               ZERO,X(LWRK3),L1)
          CALL DGEMM('N','N',L1,L1,L1,ONE,X(LWRK3),L1,SS1,L1,
     *               ZERO,X(LWRK2),L1)
        END IF
C
        !! For \gamma_{AD} and \gamma_{BC}
        CALL DGEMM('N','N',L1,L1,L1,ONE,DD1,L1,SS1,L1,
     *             ZERO,X(LWRK3),L1)
        IF (MODE.EQ.0.OR.MODE.EQ.2) THEN
          CALL DCOPY(L3,X(LWRK3),1,X(LWRK4),1)
          CALL TRPOSQ(X(LWRK4),L1)
        ELSE
          CALL DGEMM('N','N',L1,L1,L1,ONE,SS1,L1,DD2,L1,
     *               ZERO,X(LWRK4),L1)
        END IF
C
        IF (MODE.EQ.1) THEN
          CALL DGEMM('N','N',L1,L1,L1,ONE,DD2,L1,SS1,L1,
     *               ZERO,X(LWRK5),L1)
          CALL DGEMM('N','N',L1,L1,L1,ONE,SS1,L1,DD1,L1,
     *               ZERO,X(LWRK6),L1)
        END IF
      ELSE IF (MODE.EQ.3) THEN
        !! DD1 = DD2 and SS1 \= SS2
        CALL DGEMM('N','N',L1,L1,L1,ONE,SS1,L1,DD1,L1,
     *             ZERO,X(LWRK3),L1)
        CALL DGEMM('N','N',L1,L1,L1,ONE,X(LWRK3),L1,SS2,L1,
     *             ZERO,X(LWRK1),L1)
        CALL DCOPY(L3,X(LWRK1),1,X(LWRK2),1)
C
        CALL DGEMM('N','N',L1,L1,L1,ONE,DD1,L1,SS2,L1,
     *             ZERO,X(LWRK3),L1)
        CALL DGEMM('N','N',L1,L1,L1,ONE,SS1,L1,DD2,L1,
     *             ZERO,X(LWRK4),L1)
      END IF
C
      SCAL= -0.125D+00*SCAL0
      DO IAT = 1, NAT
        DO JAT = 1, IAT-1
          VEC(1) = C(1,JAT) - C(1,IAT)
          VEC(2) = C(2,JAT) - C(2,IAT)
          VEC(3) = C(3,JAT) - C(3,IAT)
          TMP = ZERO
          DO NU = IND(JAT)+1, IND(JAT)+MAXANG(ISPE(JAT))**2
            DO MU = IND(IAT)+1, IND(IAT)+MAXANG(ISPE(IAT))**2
              TMP = TMP
     *            + X(LWRK3+MU-1+L1*(NU-1))*X(LWRK4+MU-1+L1*(NU-1))
     *            + X(LWRK5+MU-1+L1*(NU-1))*X(LWRK6+MU-1+L1*(NU-1))
     *            + X(LWRK1+MU-1+L1*(NU-1))*DD2(MU,NU)
     *            + X(LWRK2+MU-1+L1*(NU-1))*DD1(MU,NU)
            END DO
          END DO
          CALL DFTB_CNVSQ(IAT,JAT,NSEQ)
          DIST = DISTMAT(NSEQ)
          DISTI= ONE/DIST
          IF (MODE.EQ.0.OR.MODE.EQ.1.OR.MODE.EQ.3) THEN
            GAMMAV = GAMDERLC(NSEQ)
            VAL = TMP*GAMMAV
            VAL = VAL*DISTI*SCAL
            DO K = 1, 3
              EGRAD(K,IAT) = EGRAD(K,IAT) - VEC(K)*VAL
              EGRAD(K,JAT) = EGRAD(K,JAT) + VEC(K)*VAL
            END DO
          ELSE IF (MODE.EQ.2) THEN
            CALL DFTB_LCHESS(IAT,JAT,DIST,GRAD,HESS)
            IGZ = 3*(IAT-1)
            JGZ = 3*(JAT-1)
            DO KI = 1, 3
              DO KJ = 1, 3
                VAL = VEC(KI)*VEC(KJ)*(HESS-DISTI*GRAD)*DISTI*DISTI
                IF (KI.EQ.KJ) VAL = VAL + DISTI*GRAD
                VAL = VAL*SCAL*TMP
                HESSIAN(IGZ+KI,IGZ+KJ) = HESSIAN(IGZ+KI,IGZ+KJ) + VAL
                HESSIAN(JGZ+KI,JGZ+KJ) = HESSIAN(JGZ+KI,JGZ+KJ) + VAL
                HESSIAN(IGZ+KI,JGZ+KJ) = HESSIAN(IGZ+KI,JGZ+KJ) - VAL
                HESSIAN(JGZ+KI,IGZ+KJ) = HESSIAN(JGZ+KI,IGZ+KJ) - VAL
              END DO
            END DO
          END IF
        END DO
      END DO
C
      CALL RETFM(NEED)
C
      RETURN
C
      END SUBROUTINE DFTB_LCGRAD2
C-----------------------------------------------------------------------
C*MODULE DFTBGR    *DECK DFTB_RMNUM
C>
C>    @brief Remove numbers, period, and comma from INC
C>
C>    @details Remove numbers, period, and comma from INC and put OUTC
C>
C>    @author Yoshio Nishimoto
C>    - Sep, 2014- Subroutine written
C>
C>           --- INPUT ---
C>    @param INC Input characters
C>           --- OUTPUT ---
C>    @param OUTC Output characters, which should not contain numbers,
C>           period, and comma
C
      SUBROUTINE DFTB_RMNUM(INC,OUTC)
C
      IMPLICIT NONE
C
      CHARACTER(8), INTENT(IN) :: INC
      CHARACTER(8), INTENT(OUT) :: OUTC
C
      INTEGER :: I,II
      CHARACTER(8) :: TMP
      CHARACTER(12) :: CHARA,BLANK
      DATA CHARA /'.,1234567890'/
      DATA BLANK /'            '/
C
      TMP = INC
      DO I = 1, 8
        II = INDEX(CHARA,TMP(I:I))
        IF (II.GT.0) TMP(I:I) = BLANK(II:II)
      END DO
      II = LEN_TRIM(TMP)
      OUTC = INC(1:II)
C
      RETURN
C
      END SUBROUTINE DFTB_RMNUM
C
