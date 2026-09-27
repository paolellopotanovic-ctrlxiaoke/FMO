c 17 Oct 19 - DGF - tweaks for FMO 5.4
C  6 Jun 18 - DGF - tweaks for FMO 5.3
C 18 Apr 16 - HN  - NEW MODULE FOR PCM HESSIANS
C
C*MODULE PCMHSS  *DECK HSSCCM
C>    @brief Second order derivatiev for PCM Hessian 
C>
C>    @details  Main driver
C>
C>    @author Hiroya Nakata
C>    - Apr, 2013- Subroutine writte
C>
      SUBROUTINE HSSCCM(imode)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /FMCOM / X(1)
      COMMON /PCMDIM/ MXSP,MXTS,MEMPCM1,MEMPCM2,NTS
      COMMON /PCMPNT/ LXYZRE,LSSFE,LLIST,LALPCM,LRINPCM,LINAPCM,LINFPCM,
     *                LINIPCM,LMEPCM,LAXYZCT,LXYZCT2,LNVERT,LQSN,LQSE,
     *                LQFS,LQIND,LISPHE,LVAD,LQOR,LVECMUL,LAIPRJ,LFIPRJ,
     *                LPEL,LPCMCDR,LDAI,LIDDAI,LQSND,LQSED,LEPSPCM
C
      CALL HSSCCMX(X(LAXYZCT),X(LAXYZCT+MXTS),X(LAXYZCT+MXTS*2),
     *             X(LAXYZCT+MXTS*3),X(LQSN),X(LQSE),
     *             X(LQSND),X(LQSED),X(LQFS),X(LQIND),
     *             X(LQIND+MXTS),X(LXYZRE),X(LXYZRE+MXSP),
     *             X(LXYZRE+MXSP*2),X(LXYZRE+MXSP*3),X(LISPHE),
     *             X(LPCMCDR),imode)
      RETURN
      END
C>    @brief Second order derivatiev for PCM Hessian 
C>
C>    @details  Main driver
C>
C>    @author Hiroya Nakata
C>    - Apr, 2013- Subroutine writte
C>   @param   XCTS    x coordinate for teserae
C>   @param   YCTS    y coordinate for teserae
C>   @param   ZCTS    z coordinate for teserae
C>   @param   AS      Areae surface charge
C>   @param   QSN     nuclear charge
C>   @param   QSE     electronic charge
C>   @param   QSND    ASC for het
C>   @param   QSED    ASC for het
C>   @param   Q_FS    ASC
C>   @param   Q_IND   ASC for induction
C>   @param   Q2      ASC
C>   @param   XE      Solute geometry informtion  
C>   @param   YE      Solute geometry informtion
C>   @param   ZE      Solute geometry informtion
C>   @param   RE      Solute geometry informtion
C>   @param   ISPHE   Indice for atom
C>   @param   PCMCDR  Derivative of coordinate
C>
C*MODULE PCMHSS  *DECK DERCCMX
      SUBROUTINE HSSCCMX(XCTS,YCTS,ZCTS,AS,QSN,QSE,QSND,QSED,
     *                   Q_FS,Q_IND,Q2,
     *                   XE,YE,ZE,RE,ISPHE,PCMCDR,imode)
      use mx_limits, only: mxatm
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C
      LOGICAL GOPARR, DSKWRK, MASWRK
c    , TDDFT, MP2PCM
      LOGICAL DFTBFL,SCC,SRSCC,DFTB3,DAMPXH,LCDFTB
C
      DIMENSION XCTS(*),YCTS(*),ZCTS(*),AS(*),QSN(*),QSE(*),
     *          QSND(*),QSED(*),Q_FS(*),
     *          Q_IND(*),Q2(*),XE(*),YE(*),ZE(*),RE(*),ISPHE(*),
     *          PCMCDR(3,*)
C
      COMMON /DFTB  / DFTBFL,SCC,SRSCC,DFTB3,DAMPXH,LCDFTB
      COMMON /FMCOM / XX(1)
      COMMON /FMOINF/ NFG,NLAYER,NATFMO,NBDFG,NAOTYP,NBODY,NSEGM
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
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON /MACHIN/ NWDVAR,MAXFM,MAXSM,LIMFM,LIMSM
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /PCMDIM/ MXSP,MXTS,MEMPCM1,MEMPCM2,NTS
      COMMON /PCMOPT/ RABI,RASC,REFPOL,THRSLS,DENSLS,WB,WA,ETA2,GD,EVAC,
     *                RHOW,PM,AREATL,AREAKP,BONDRY,OMEGA,RET,FRO,EPSINF,
     *                EPS,DR,RSOLV,VMOL,TCE,STEN,DSTEN,CMF,TABS,IDIRCT,
     *                IPCDER,IDP,ICOMP,IFIELD,ICAV,IDISP,IPRINT,IRETCAV,
     *                ICENT,IFAST,NEVAL,IEFPOL,KEEPSM,IMGABI,IMGASC,NADD
      COMMON /PCMPAR/ IPCM,NFT26,NFT27,IKREP,IEF,IP_F,NFMOPCM,IHET
      COMMON /PCMPNT/ LXYZRE,LSSFE,LLIST,LALPCM,LRINPCM,LINAPCM,LINFPCM,
     *                LINIPCM,LMEPCM,LAXYZCT,LXYZCT2,LNVERT,LQSN,LQSE,
     *                LQFS,LQIND,LISPHE,LVAD,LQOR,LVECMUL,LAIPRJ,LFIPRJ,
     *                LPEL,LPCMCDR,LDAI,LIDDAI,LQSND,LQSED,LEPSPCM
      COMMON /RUNOPT/ RUNTYP,EXETYP,NEVALS,NGLEVL,NHLEVL
C
      DATA CHECK/8HCHECK   /
c     data excite,spnflp/8HEXCITE  ,8HSPNFLP  /
C
C     ----- C-PCM AND IEF-PCM GRADIENT CORRECTIONS -----
C
c     TDDFT=.FALSE.
c     MP2PCM=.false.
c     IF(TDDFTYP.NE.RNONE) TDDFT=.TRUE.
c     MP2PCM=MPLEVL.EQ.2
c     IF (ISWNEW1.EQ.1) MP2PCM = .FALSE.
c     L1 = NUM
      L2 = (NUM*NUM+NUM)/2
C
C
      LFIXPV = 0
      IF(IPCDER.EQ.3) LFIXPV = 3*20*MXTS
      LFIXPV2 = 0
      IF(IPCDER.EQ.3) LFIXPV2= 21*MXTS/NWDVAR + 1
C
      
      if(nfg.ne.0) then
        CALL VALFM(LOADFM)
        liwrk  = LOADFM + 1
        Lbglb  = liwrk  + natfmo
        LAST   = lbglb  + nat
        NEEDFMO= LAST   - LOADFM -1
        CALL GETFM(NEEDFMO)
C
        call mapfmohess(natfmo,xx(LIAGLOB),xx(LINDAT),xx(LLAYFRG),
     *     xx(liwrk),natb,4)
        call mapfmohess(nat,xx(LIAGLOB),xx(LINDAT),xx(LLAYFRG),
     *     xx(lbglb),natc,0)
CC
c       do iii = 1,natfmo
c         write(6,'(I4)') IXFTCH(XX(liwrk),iii)
c       end do
c       do iii = 1,nat
c         write(6,'(I4)') IXFTCH(XX(lbglb),iii)
c       end do
C
c       write(6,*) "check natb=",nat,natb,natc,natfmo
      end if
C
C
      CALL VALFM(LOADFM)
      LD     = LOADFM + 1
      LSCR   = LD     + L2
      LFLD   = LSCR   + L2
      LFLW   = LFLD   + 3*L2
      LSOLQM = LFLW   + 3*225
c     LSOLQM = LSOLEF + 6*NFRG
C   For gradient/hess/FD
      LEH    = LSOLQM + 3 * NAT
      LFEQ   = LEH    + 9 * (NAT*NAT+NAT)/2
      LFEQWRK= LFEQ   + 3 * NAT*NTS
      LFD    = LFEQWRK+ 9 * (NAT*NAT+NAT)/2
      LEHQ   = LFD    + 3 * NAT*L2
      if(nfg.ne.0) LEHQ   = LFD    + 3 * natb * L2
      LFDQ   = LEHQ   + 9 * (NAT*NAT+NAT)/2
      IMAP   = LFDQ   + 3 * NAT*L2
      LQT    = IMAP   + NAT + 1
      if(nfg.ne.0) then
         LQT = IMAP   + NATFMO + 1
      endif
C
      NAT2   = (NAT*NAT+NAT)/2
C      
C
      LVECTMP= LQT    + NTS
      if(imode.eq.1) LVECTMP= LQT    + NTS * 2
      LAST   = LVECTMP+ NTS
c     No efp currently
c     LAST   = LFFRGQ + 6*NFRG*NTS
C
c     LTA    = LAST
c     LV1    = LAST
      LQSEXY = LAST
      LQSEMN = LAST
c     LPMP   = LAST
      LQSEMP = LAST
      IF(NFMOPCM.EQ.0) THEN
C       ELSE ALREADY ALLOCATED IN PCMMEM2
        LDAI=LAST
        LIDDAI = LDAI   + LFIXPV
        LAST   = LIDDAI + LFIXPV2
      ENDIF
      IF(DFTBFL) THEN
         LZAN  = LAST
         LAST  = LZAN   + NAT
      END IF
      NEED     = LAST   - LOADFM - 1
      CALL GETFM(NEED)
      IF(EXETYP.EQ.CHECK) GO TO 600
C

      if(nfg.eq.0) natM = nat 
      if(nfg.ne.0) natM = natfmo 
      
      DO IAT=1,NATM
        ICK = 0
        DO ITS=1,NTS
         IF(IAT.EQ.ISPHE(ITS).and.ICK.EQ.0) THEN 
           call ixstor(xx(IMAP),IAT,ITS)
           ICK  = 1
         END IF
        END DO
        IF(ICK.EQ.0) call ixstor(xx(IMAP),IAT,0)
      END DO
C
      call ixstor(xx(IMAP),NATM+1,NTS+1)
C
      do iat = natm,1,-1
        its =  ixftch(xx(imap),iat)
        if(its.eq.0) then
          jts  = ixftch(xx(imap),iat+1)
          call ixstor(xx(IMAP),IAT,jts)
        end if
      end do
C
C
C
C
C
      IF(IPCDER.EQ.2.OR.IPCDER.EQ.3) THEN
c       LEXT1=LAST
        LEXT2=LAST
C       LFMOPG AND LIAGLOB ARE ONLY DEFINED FOR FMO!
        IF(NFG.NE.0) THEN
c        LEXT1=LFMOPG
         LEXT2=LIAGLOB
        ENDIF
C
        IF(IPCDER.EQ.3.AND.NFMOPCM.EQ.0) THEN
         CALL DAREAD(IDAF,IODA,XX(LDAI),LFIXPV,334,0)
         CALL DAREAD(IDAF,IODA,XX(LIDDAI),LFIXPV2,335,1)
        ENDIF
c
C
        IF(DFTBFL) CALL DFTB_GET_ZAN(0,NAT,XX(LZAN),ZAN)
        if(nfg.eq.0) then 
          CALL CCMBEMH(XX(LD),XX(LSOLQM),XX(LFEQ),XX(LQT),XX(LVECTMP),
     *              L2,NAT,XCTS,YCTS,ZCTS,AS,QSN,QSE,QSND,QSED,Q_FS,
     *              Q_IND,Q2,XE,YE,ZE,RE,ISPHE,
     *              PCMCDR,XX(LEXT2),XX(LDAI),XX(LIDDAI),
     *              XX(LQSEXY),XX(LQSEMN),XX(LFEQWRK),
     *              XX(LQSEMP),XX(LEH),XX(LFD),XX(LEHQ),
     *              XX(LFDQ),nat2,xx(IMAP),XX(LZAN))
        else if(nfg.ne.0) then
          CALL CCMBEMFMOH(XX(LD),XX(LSOLQM),XX(LFEQ),XX(LQT),
     *    XX(LVECTMP),L2,NAT,XCTS,YCTS,ZCTS,AS,QSN,QSE,QSND,QSED,Q_FS,
     *              Q_IND,Q2,XE,YE,ZE,RE,ISPHE,
     *              PCMCDR,XX(LEXT2),XX(LDAI),XX(LIDDAI),
     *              XX(LQSEXY),XX(LQSEMN),XX(LFEQWRK),
     *              XX(LQSEMP),XX(LEH),XX(LFD),XX(LEHQ),
     *              XX(LFDQ),nat2,xx(IMAP),imode,xx(liwrk),natb,
     *              XX(lbglb))
        end if
        IF(DFTBFL) CALL DFTB_GET_ZAN(1,NAT,XX(LZAN),ZAN)
      END IF
C
  600 CONTINUE
      CALL RETFM(NEED)
C
      if(nfg.ne.0)  CALL RETFM(NEEDFMO)
C
      IF(MASWRK) WRITE(IW,*)
     *      '.... DONE WITH PCM CONTRIBUTION TO HESSIAN ....'
C
      CALL TIMIT(1)
      RETURN
      END
C
C*MODULE PCMHSS  *DECK HSSCCM
C>    @brief Second order derivatiev for PCM Hessian 
C>
C>    @details  Main solver for Hessian
C>
C>    @author Hiroya Nakata
C>    - Apr, 2013- Subroutine writte
C>
C
C>   @param   XCTS    x coordinate for teserae
C>   @param   YCTS    y coordinate for teserae
C>   @param   ZCTS    z coordinate for teserae
C>   @param   AS      Areae surface charge
C>   @param   QSN     nuclear charge
C>   @param   QSE     electronic charge
C>   @param   QSND    ASC for het
C>   @param   QSED    ASC for het
C>   @param   Q_FS    ASC
C>   @param   Q_IND   ASC for induction
C>   @param   Q2      ASC
C>   @param   XE      Solute geometry informtion  
C>   @param   YE      Solute geometry informtion
C>   @param   ZE      Solute geometry informtion
C>   @param   RE      Solute geometry informtion
C>   @param   ISPHE   Indice for atom
C>   @param   PCMCDR 

C>   @param   DEN     Density matrix
C>   @param   DSOLQM  Derivative of coordinate
C>   @param   FEQ     Grdient elements for coordinate
C>   @param   QT      total charge
C>   @param   VECTMP  potential
C>   @param   L2      number of triangular marix
C>   @param   NATM    number of atom
C>   @param   XCTS    x coordinate for teserae
C>   @param   YCTS    y coordinate for teserae
C>   @param   ZCTS    z coordinate for teserae
C>   @param   AS      Surface charge
C>   @param   QSN     nuclear charge
C>   @param   QSE     electronic charge
C>   @param   QSND    ASC for het
C>   @param   QSED    ASC for het
C>   @param   Q_FS    ASC
C>   @param   Q_IND   ASC for induction
C>   @param   Q2      ASC
C>   @param   XE      Solute geometry informtion  
C>   @param   YE      Solute geometry informtion
C>   @param   ZE      Solute geometry informtion
C>   @param   RE      Solute geometry informtion
C>   @param   ISPHE   Indice for atom
C>   @param   PCMCDR  Grdient
C>   @param   FMOPG   FMO gradient
C>   @param   IAGLOB  global indice
C>   @param   DAI     PCM parameter
C>   @param   IDDAI   PCM parameter indice
C>   @param   QSEXY   FMO charge
C>   @param   QSEMN   FMO charge
C>   @param   FEQWRK  Work for derivative
C>   @param   QSEMP   FMO  charge
C>   @param   EH      Hessian
C>   @param   FD      Fock derivative
C>   @param   EHQ     work for hessian
C>   @param   FDQ     work for derivative
C>   @param   nat2    triangular matrix
C>   @param   MAP     map indice
C>   @param   ZANTMP  temp
      SUBROUTINE CCMBEMH(DEN,DSOLQM,FEQ,QT,VECTMP,
     *           L2,NATM,XCTS,YCTS,ZCTS,AS,QSN,QSE,QSND,QSED,Q_FS,
     *           Q_IND,Q2,XE,YE,ZE,RE,ISPHE,
     *           PCMCDR,IAGLOB,DAI,IDDAI,
     *           QSEXY,QSEMN,FEQWRK,QSEMP,
     *           EH,FD,EHQ,FDQ,nat2,MAP,ZANTMP)
      use mx_limits, only: mxatm
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      PARAMETER (MXSPE=10) 
c     MXFRG=1050, MXDFG=5, MXDPPT=MXFRG*MXDFG*12)
      PARAMETER (ONE=1.0D+00)
C
      LOGICAL GOPARR, DSKWRK, MASWRK
      LOGICAL ISGDDI,PAROUT,INITGDDI,wasgddi,MLGDDI
c     LOGICAL SAVGOP
c          MOPAC
      LOGICAL LFMOPCM1
      LOGICAL DFTBFL,SCC,SRSCC,DFTB3,DAMPXH,LCDFTB
C
      DIMENSION DEN(L2),DSOLQM(3,NATM),
     *          FEQ(3,NATM,*),QT(NTS),VECTMP(*),
     *          XCTS(*),YCTS(*),ZCTS(*),AS(*),QSN(*),QSE(*),
     *          QSND(*),QSED(*),Q_FS(*),
     *          Q_IND(*),Q2(*),XE(*),YE(*),ZE(*),RE(*),ISPHE(*),
     *          PCMCDR(3,*),IAGLOB(*),
     *          QSEXY(*),QSEMN(*),FEQWRK(9,*),QSEMP(*),
     *          EH(9,*),FD(L2,3,natm),
     *          EHQ(9,nat2),FDQ(L2,3,natm,*),MAP(*),ZANTMP(*)
      DIMENSION DAI(3,20,MXTS),IDDAI(21,MXTS)
C
      COMMON /DFTB  / DFTBFL,SCC,SRSCC,DFTB3,DAMPXH,LCDFTB
      COMMON /DFTBPR/ ETEMP,DFTBDP(MXSPE*14),DAMPXHE,HUBDER(MXSPE),
     *                ZREF(MXATM),SPNCST(6,MXSPE),SPE(MXATM),NSPE,
     *                MAXANG(MXATM),ISPE(MXATM),IND(MXATM+1),IDFTBD,
     *                PARAMDIR
      COMMON /FMCOM / XX(1)
      COMMON /GDDI/   ISCOPE,NGROUPS,MYGROUP,MEGLOB,NPGLOB,NNGLOB,JBTYP,
     *                ISGDDI,PAROUT,INITGDDI,wasgddi,MLGDDI,NSUBGR,
     *                MeUniv,NPUniv,numdlb,myworld,nworlds,mogddi
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON /MLTFMO/ Q_MUL,IMLTFMO,ISWPFD,ISWNEW1,lfmoicm
      COMMON /MLTPNT/ LDELEPOT,LDELQSE,LQSETMP,LFZVEC,LMAPCMFD,LQBL2FD,
     *                LFRZQ
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /PCMDIM/ MXSP,MXTS,MEMPCM1,MEMPCM2,NTS
      COMMON /PCMOPT/ RABI,RASC,REFPOL,THRSLS,DENSLS,WB,WA,ETA2,GD,EVAC,
     *                RHOW,PM,AREATL,AREAKP,BONDRY,OMEGA,RET,FRO,EPSINF,
     *                EPS,DR,RSOLV,VMOL,TCE,STEN,DSTEN,CMF,TABS,IDIRCT,
     *                IPCDER,IDP,ICOMP,IFIELD,ICAV,IDISP,IPRINT,IRETCAV,
     *                ICENT,IFAST,NEVAL,IEFPOL,KEEPSM,IMGABI,IMGASC,NADD
      COMMON /PCMPAR/ IPCM,NFT26,NFT27,IKREP,IEF,IP_F,NFMOPCM,IHET
C
      DATA ZERO/0.0D+00/
c     DATA RHF,UHF,ROHF,GVB/8HRHF     ,8HUHF     ,8HROHF    ,8HGVB     /
c     DATA NONE,RNONE/4HNONE,8HNONE    /
c     data excite,spnflp/8HEXCITE  ,8HSPNFLP  /
C---------------------------------------------------------------------
C     COMPUTE THE C-PCM AND IEF-PCM Hessian 
C---------------------------------------------------------------------
C
C     IN AN QM/EFP/PCM RUN, AN EFP(EFPA) AND ITS
C     ASC(ASCA) EXPERIENCE THE FOLLOWING FORCES/TORQUES:
C        EFPA <-- OTHER ASC
C        ASCA <-- OTHER ASC, OTHER EFP, QM ATOMS
C     FOR AN QM ATOM(ATOMB) AND ITS ASC(ASCB):
C        ATOMB <-- OTHER ASC
C        ASCB  <-- OTHER ASC, OTHER QM ATOMS, EFP
C     ***********************************************
C
C
C     ISWPFD=1: FOR FMO/PCM<1> COMPLETELY ANALYTIC GRADIENT
C               NOTE THAT ONLY FMO/CPCM[1] COMPLETELY ANALYTIC GRADIENT
C                         IS CURRENTLY AVAILABLE.
      if(nat.eq.0) write(6,*) "PCMCDR=",PCMCDR(1,1)
C
      LFMOPCM1 = ISWPFD.EQ.1
      LFMOPCM1 = .false.
c     IF (ISWNEW1.EQ.1) MP2PCM = .FALSE.
c     IF (LFMOPCM1) LTMPTMP = LQSETMP - 1
C
c     IF(NFMOPCM.NE.0) CALL VCLR(FMOPG,1,3*NATFMO)
C     -- READ GRAD. HESS. and Fock^a 
c     CALL DAREAD(IDAF,IODA,DSOLQM,NEGH,67,0)
c     write(6,*) "Check val0"
c     CALL EGOUT(DSOLQM,NAT)
      L1   = NUM
      L3   = L1*L1
      NEGH = 3*NAT + 9*(NAT*NAT+NAT)/2
      CALL VCLR(DSOLQM,1,NEGH)
      CALL VCLR(FD,1,3*NAT*L2)
C
C
C     -- READ DENSITY MATRIX
C
      IF(NAT.GT.0) THEN
         CALL DAREAD(IDAF,IODA,DEN,L2,16,0)
c        IF(SCFTYP.EQ.ROHF .OR. SCFTYP.EQ.UHF
c    *                     .OR. SCFTYP.EQ.GVB) THEN
c           CALL VALFM(LOADFM)
c           LDENB = LOADFM + 1
c           LAST  = LDENB  + L2
c           NEED  = LAST - LOADFM - 1
c           CALL GETFM(NEED)
c           CALL DENDD1(DEN,XX(LDENB),L2)
c           CALL RETFM(NEED)
c        END IF
      END IF
C
c     -- COMBINE THE SURFACE CHARGES
c     write(6,'("Check charge")')
c     write(6,'("VAL="4F12.7)') 
c    *     (QSE(III),QSN(III),Q_FS(III),Q_IND(III),iii=1,nts)
C
      DO ITS=1,NTS
         QT(ITS)=QSE(ITS)+QSN(ITS)+Q_FS(ITS)+Q_IND(ITS)
      ENDDO
C
C     -- ZERO OUT Q2 FOR IEF=10
C
      IF(IEF.EQ.10) THEN
         DO ITS=1,NTS
            Q2(ITS)=ZERO
         END DO
      END IF
C
c     DO ITS=1,nts
c       write(6,*) "ISP=",ITS,ISPHE(its)
c     END DO
C
      IF (DFTBFL) THEN
        CALL VALFM(LOADFM)
        LWRK = LOADFM + 1
        LAST = LWRK   + NAT
        NEED = LAST   - LOADFM -1
        CALL GETFM(NEED)
        CALL DAREAD(IDAF,IODA,XX(LWRK),NAT,556,0)
        CALL DAXPY(NAT,-1.0D+00,XX(LWRK),1,ZAN,1)
      END IF
C
      DO 110 IAT=1,NAT
       NINITS=MAP(IAT)
       NENDTS=MAP(IAT+1)
       NLEN  =NENDTS-NINITS
       IF(NLEN.EQ.0) GO TO 110
c      write(6,'("CHK DIM=",3I4)') NINITS,NENDTS,NLEN
C
       CALL VCLR(FEQ,1,3*NAT)
       CALL VCLR(EHQ,1,9*nat2)
       CALL VCLR(FEQWRK,1,9*nat2)
       CALL VCLR(FDQ,1,3*L2*NAT)
C      TVD! One-Electron interaction
       IF (.NOT.DFTBFL)
     * CALL PCMTVDSPD(FEQ,EHQ,DEN,FDQ,XCTS(NINITS),YCTS(NINITS),
     *      ZCTS(NINITS),NLEN,QT(NINITS))
C      NUCLEAR REPULSION
       CALL PCMNUCHSS(FEQ,EHQ,XCTS(NINITS),YCTS(NINITS),
     *      ZCTS(NINITS),NLEN,QT(NINITS),IAT)
C      Translation invraiance is used
       CALL TRANINVHSS(IAT,NAT,FEQ,FDQ,EHQ,DSOLQM,FD,EH,FEQWRK,L2,nat2)
C
c      WRITE(IW,*) "wwwcheck hessian",iat
c      CALL HSSPRT(NAT,DSOLQM,EH)
C
 110  CONTINUE
C
c     WRITE(IW,*) "wwwcheck hessian"
c     CALL HSSPRT(NAT,DSOLQM,EH)
C
      IF (DFTBFL) THEN
        CALL DAXPY(NAT,1.0D+00,XX(LWRK),1,ZAN,1)
        CALL RETFM(NEED)
      END IF
C     check
c     CALL VCLR(FD,1,3*NAT*L2)
c     do i = 1,nat
c      do j=1,3
c       val = TRACEP(DEN,FD(1,j,i),num) 
c       write(*,'("gradient ="2I4,F10.7)') i,j,val
c      end do
c     end do
c     write(6,*) "check grd."
c     CALL EGOUT(DSOLQM,NAT)
C
C
      CALL VALFM(LOADFM)
      NDASC = 0
      LFPNT = LOADFM
      NTSPAR=(NTS-1)/npglob + 1
      CALL INDQPROD(0,NTSPAR,LFPNT,NDASC,DUM,DUM,DUM,0,0)
      LQTMP   = LFPNT   + NDASC  + 1
      LWRK    = LQTMP   + nts    
      LFCMTMP = LWRK    + L2
      LAST    = LFCMTMP + nat * nat * 9
      NEEDW   = LAST    - LOADFM    - 1
      CALL GETFM(NEEDW)
C
C     V^a * C^{-1} * V^b ||  V(ij) * C^{-1} V^a (Derivative of W)
c     CALL VCLR(FEQ,1,3*NAT*nts)
      IF (LFMOPCM1) THEN
        CALL FELECQ(DEN,XX(LQSETMP),FEQ,L2,XCTS,YCTS,ZCTS)
      ELSE
        call dacopy(nts,one,VECTMP,1)
        CALL FELECQ(DEN,VECTMP,FEQ,L2,XCTS,YCTS,ZCTS)
        CALL FNUCNT(FEQ,XCTS,YCTS,ZCTS,NTS,ISPHE)
      END IF
C
c     SAVGOP=GOPARR
c     GOPARR=.false.
C
      IF(GOPARR) THEN
        CALL DSCAL(NAT*3,ONE/NPROC,DSOLQM,1)
        CALL DSCAL(NAT*3*NTS,ONE/NPROC,FEQ,1)
      END IF
C     C^ab
C
      IF (DFTBFL) CALL DFTB_GET_ZAN(1,NAT,ZANTMP,ZAN)
      CALL QCHESS(NATM,0,FEQWRK,QT,VECTMP,XCTS,YCTS,ZCTS,
     *            AS,QSN,QSE,QSND,QSED,Q_FS,Q_IND,Q2,XE,YE,ZE,
     *            RE,ISPHE,DAI,IDDAI,XX(LQSETMP),QSEMP,QSEXY,QSEMN,
     *            DSOLQM,LFMOPCM1,IAGLOB,XX(LFCMTMP),FEQ,DUM,0)
      IF (DFTBFL) CALL DCOPY(NAT,ZREF,1,ZAN,1)
c
C     - SUM UP -
      IF(GOPARR) THEN 
        CALL DDI_GSUMF(2461,DSOLQM,3*NAT)
        CALL DDI_GSUMF(2462,XX(LFCMTMP),9*NAT*NAT)
        CALL DDI_GSUMF(2463,FEQ,NAT*3*NTS)
      END IF
C
c     GOPARR=SAVGOP
C     - C^{-1} * V^a
c     call pcmaddhss(EH,XX(LFCMTMP),NAT)
c     CALL VCLR(xx(lfcmtmp),1,9*nat*nat)
C
      if (dftbfl) go to 100
      do iat=1,nat
       do ixyz=1,3
C        FEQ --> VECTMP
         DO ITS=1,NTS
          VECTMP(ITS)=FEQ(ixyz,iat,its)
         END DO
         ITER0 = 1
c        write(6,*) "www get q^a"
         CALL INDQPROD(ITER0,NTSPAR,LFPNT,NDASC,VECTMP,
     *                 XX(LQTMP),XX(LQTMP),1,0)
C        Fock derivative terms
         CALL PCMPOT(XX(lwrk),XX(lwrk),NTS,XX(LQTMP),XCTS,L2,1,.false.)
         CALL DAXPY(L2,1.0D+00,xx(lwrk),1,FD(1,ixyz,iat),1)
c        Complete HESS cont
         CALL PCMCMPHSS(XX(LFCMTMP),NTS,XX(LQTMP),FEQ,NAT,iat,ixyz)
       end do
      end do
c     call prsq(xx(lfcmtmp),nat*3,nat*3,nat*3)
  100 continue
      call pcmaddhss(EH,XX(LFCMTMP),NAT)
C
      CALL RETFM(NEEDW)
C
      !! Nishimoto added on 2016/10/05. It is correct?
      IF(ICAV.EQ.1.OR.IDISP.EQ.1)CALL DAXPY(3*NAT,ONE,PCMCDR,1,DSOLQM,1)
C
C     Save gradient and Hessian
      IF (.NOT.DFTBFL) CALL DAREAD(IDAF,IODA,FEQ,NEGH,67,0)
      IF (.NOT.DFTBFL) CALL DAXPY(NEGH,1.0D+00,FEQ,1,DSOLQM,1)
      CALL DAWRIT(IDAF,IODA,DSOLQM,NEGH,67,0)
      IF (DFTBFL)  RETURN
C     Save Fock derivative
      CALL VALFM(LOADFM)
      LVEC  = LOADFM + 1
      LWRK1 = LVEC   + L3
      LWRK2 = LWRK1  + L3
      LWRK3 = LWRK2  + L3
      LAST  = LWRK3  + L3
      NEED  = LAST   - LOADFM - 1
      CALL GETFM(NEED)
C
      CALL DAREAD(IDAF,IODA,XX(LVEC),L3,15,0)
C
      NFT18 = 18
c     NXYZ  = NAT * 3
      CALL SEQREW(NFT18)
      DO IAT=1,NAT
       DO IXYZ=1,3
C
c       write(6,*) "www fd check 0 =",iat,ixyz
c       CALL PRTRI(FD(1,IXYZ,IAT),L1)
C
        CALL TFTRI(XX(LWRK1),FD(1,IXYZ,IAT),XX(LVEC),XX(LWRK2),L1,L1,L1)
        CALL SQREAD(NFT18,XX(LWRK3),L2)
        CALL VADD(XX(LWRK3),1,XX(LWRK1),1,FD(1,IXYZ,IAT),1,L2)
       END DO
      END DO
      CALL SEQREW(NFT18)
      DO IAT=1,NAT
       DO IXYZ=1,3
C
c        write(6,*) "www fd check =",iat,ixyz
c        CALL PRTRI(FD(1,IXYZ,IAT),L1)
C
         CALL SQWRIT(NFT18,FD(1,IXYZ,IAT),L2)
       END DO
      END DO
C
      CALL RETFM(NEED)
c     write(6,*) "Check final grad."
c     CALL EGOUT(DSOLQM,NAT)
C
      RETURN
      END
C
C
C*MODULE PCMHSS  *DECK PCMTVDSPD
C>    @brief Second order derivatiev for PCM 
C>
C>    @details  Derivative of Wuv
C>
C>    @author Hiroya Nakata
C>    - Apr, 2013- Subroutine writte
C>
C>   @param   EG      gradient
C>   @param   EH      Hessian
C>   @param   DAB     Density matrix
C>   @param   FD      Fock derivative
C>   @param   XCTS    x coordinate for teserae
C>   @param   YCTS    y coordinate for teserae
C>   @param   ZCTS    z coordinate for teserae
C>   @param   NTS     number of tesserae
C>   @param   QT      total charge
      SUBROUTINE PCMTVDSPD(EG,EH,DAB,FD,XCTS,YCTS,ZCTS,NTS,QT)
      use mx_limits, only: mxgtot,mxsh,mxgsh,mxatm
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION FOE(36,9)
C
      LOGICAL DBUG,OUT,GOPARR,DSKWRK,MASWRK,NXT
      LOGICAL FIRST,SECND,CPHF,BOTH,MFIRST,MSECND,MCPHF,BFDD
C
      DIMENSION EG(3,*),EH(9,*),DAB(*),FD(*)
      DIMENSION DIJ10( 78),DIJ11(169),DIJ20(132),
     1          FFT10( 78),FFT11(169),FFT20(132),
     2          GIJ10( 78),GIJ11(169),GIJ20(132),
     3          XIN(   80),YIN(   80),ZIN(   80)
      DIMENSION DF(3),DS(9,2),DOE(108)
      DIMENSION XCTS(*),YCTS(*),ZCTS(*),QT(*)
C
C
      COMMON /HSSPAR/ FIRST,SECND,CPHF,BOTH,MFIRST,MSECND,MCPHF
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /IOFILE/ IR,IW,IP,IJK,IPK,IDAF,NAV,IODA(950)
      COMMON /NSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),
     *                CF(MXGTOT),CG(MXGTOT),CH(MXGTOT),CI(MXGTOT),
     *                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),
     *                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
      COMMON /OUTPUT/ NPRINT,ITOL,ICUT,NORMF,NORMP,NOPK
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /RESTAR/ TIMLIM,NREST,NREC,INTLOC,IST,JST,KST,LST
      COMMON /ROOT  / XX,U(13),W(13),NROOTS
      COMMON /SHLOED/ LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,
     *                MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL,
     *                NIJ,IJD,KLD,IJ,KL
      COMMON /XYZDER/ XINT,YINT,ZINT,T,X0,Y0,Z0,XI,YI,ZI,XJ,YJ,ZJ,
     *                NI,NJ,CX,CY,CZ
      COMMON /INFOED/ GA(MXGSH),CSA(MXGSH),CPA(MXGSH),CDA(MXGSH),
     *                GB(MXGSH),CSB(MXGSH),CPB(MXGSH),CDB(MXGSH),
     *                GC(MXGSH),CSC(MXGSH),CPC(MXGSH),CDC(MXGSH),
     *                GD(MXGSH),CSD(MXGSH),CPD(MXGSH),CDD(MXGSH),
     *                XA,YA,ZA,XB,YB,ZB,RR,XC,YC,ZC,XD,YD,ZD,RRC,
     *                NGA,NGB,NGC,NGD
      COMMON /DIJOED/ NADIJ1( 78),NADIJ2(169),NADIJ3(132)
      COMMON /IJXOED/ NPIJX1( 78),NPIJX2(169),NPIJX3(132)
      COMMON /IJYOED/ NPIJY1( 78),NPIJY2(169),NPIJY3(132)
      COMMON /IJZOED/ NPIJZ1( 78),NPIJZ2(169),NPIJZ3(132)
C
C     COSMO INFORMATION
C
c     PARAMETER (MAXDEN=25*MXATM, MXABC=6000)
c     LOGICAL ISEPS,USEPS
c     COMMON /ISEPS / ISEPS,USEPS
c     LOGICAL COSBUG,COSWRT,DCOSMO,PRFCND,IOUTCH
c     COMMON /COSDAT/ SE2,SECORR,QVCOSMO,ELAST,EMP2COS,EMP2LAST,
c    *                COSVOL,COSSAR,EDIEL,EOC1,DEOC_RS,SUMQSC,
c    *                SUMQSCOLD,ZSUM,ZSUM2,ZSUM3,FEPSI,RDS,DISEX2,
c    *                EPSI,COSRAD,DISEX,OUTCHG,EDIEL_SAVE,
c    *                MAXNPS,ICORR,ITRIP,NQS,MP2TRIP,MP2ITER,
c    *                ICFREQ,NSPA,NSPH,NPSD,NPS,NPS2,NDEN,NPSPHER,
c    *                COSBUG,COSWRT,DCOSMO,PRFCND,IOUTCH
c     COMMON /COSVCE/ QSCNET(MXABC),CORZAN(3,MXABC),
c    *                QDEN(MAXDEN),QSCNET_SAVE(MXABC),IATSP(MXABC+1)
C
      DATA PI212 /1.1283791670955D+00/
      DATA ZERO,ONE,TWO /0.0D+00,1.0D+00,2.0D+00/
      DATA RLN10 /2.30258D+00/
      DATA NFDIJ1,NFDIJ2,NFDIJ3,NW / 78,169,132, 80/
C
      OUT =.FALSE.
      BFDD=.FALSE.
      DBUG=.FALSE. .AND. MASWRK
      IF (OUT  .AND.  MASWRK) WRITE (IW,9008)
      TOL = RLN10*ITOL
      NRD=0
      IF(MFIRST) NRD = 1
      IF(MSECND) NRD = 2
      L2=(NUM*(NUM+1))/2
C
c     write(6,*) "Check=",MFIRST,MSECND,MCPHF
C
      CALL STVJDD
C
C INITIALIZE PARALLEL
C
      NXT = IBTYP.EQ.1
      IPCOUNT = ME - 1
      NEXT = -1
      MINE = -1
C
C     ----- I SHELL -----
C
      DO 1060 II = 1,NSHELL
C
C     ----- GO PARALLEL! -----
C
      IF (NXT .AND. GOPARR) THEN
         MINE = MINE + 1
         IF (MINE.GT.NEXT) CALL DDI_DLBNEXT(NEXT)
         IF (NEXT.NE.MINE) GO TO 1060
      END IF
C
      IAT = KATOM(II)
      IFDX  = 3*L2*(IAT-1)
      IFDY  = IFDX + L2
      IFDZ  = IFDY + L2
      XI = C(1,IAT)
      YI = C(2,IAT)
      ZI = C(3,IAT)
      XA = XI
      YA = YI
      ZA = ZI
      I1 = KSTART(II)
      I2 = I1+KNG(II)-1
      LIT = KTYPE(II) + NRD
      MINI = KMIN(II)
      MAXI = KMAX(II)
      LOCI = KLOC(II)-MINI
      NGA = 0
      DO 10 I = I1,I2
      NGA = NGA + 1
      GA( NGA) = EX(I)
      CSA(NGA) = CS(I)
      CPA(NGA) = CP(I)
      CDA(NGA) = CD(I)
   10 CONTINUE
C
C     ----- J SHELL -----
C
      DO 1040 JJ = 1,NSHELL
C
C           GO PARALLEL!
C
      IF ((.NOT.NXT) .AND. GOPARR) THEN
         IPCOUNT = IPCOUNT + 1
         IF (MOD(IPCOUNT,NPROC).NE.0) GO TO 1040
      END IF
C
      JAT = KATOM(JJ)
      XJ = C(1,JAT)
      YJ = C(2,JAT)
      ZJ = C(3,JAT)
      XB = XJ
      YB = YJ
      ZB = ZJ
      J1 = KSTART(JJ)
      J2 = J1+KNG(JJ)-1
      LJT = KTYPE(JJ)
      IF(MSECND) LJT = LJT + 1
      MINJ = KMIN(JJ)
      MAXJ = KMAX(JJ)
      LOCJ = KLOC(JJ)-MINJ
      NGB = 0
      DO 20 J = J1,J2
      NGB = NGB + 1
      GB( NGB) = EX(J)
      CSB(NGB) = CS(J)
      CPB(NGB) = CP(J)
      CDB(NGB) = CD(J)
   20 CONTINUE
      IF(MFIRST) NROOTS = (LIT+LJT-2  )/2 + 1
      IF(MSECND) NROOTS = (LIT+LJT-2  )/2 + 1
      RR = (XI-XJ)**2+(YI-YJ)**2+(ZI-ZJ)**2
C
C     ----- PREPARE INDICES FOR PAIRS OF (I,J) FUNCTIONS
C
      CALL STVIDX(0)
C
      CALL VCLR(GIJ10,1,NFDIJ1)
      CALL VCLR(GIJ11,1,NFDIJ2)
      CALL VCLR(GIJ20,1,NFDIJ3)
C
C     ----- I PRIMITIVE
C
      IIG = 0
      DO 840 IG = I1,I2
      IIG = IIG + 1
      AI = EX(IG)
      ARRI = AI*RR
      AXI = AI*XI
      AYI = AI*YI
      AZI = AI*ZI
      DUM = AI+AI
C
C     ----- J PRIMITIVE
C
      JJG = 0
      DO 820 JG = J1,J2
      JJG = JJG + 1
      AJ = EX(JG)
      AA = AI+AJ
      AA1 = ONE/AA
      DUM = AJ*ARRI*AA1
      IF (DUM .GT. TOL) GO TO 820
C
      AX = (AXI+AJ*XJ)*AA1
      AY = (AYI+AJ*YJ)*AA1
      AZ = (AZI+AJ*ZJ)*AA1
C
C     ----- DENSITY FACTOR
C
      CALL STVDIJ(0,IIG,JJG,NIJ10,NIJ11,NIJ20,
     1            DIJ10,DIJ11,DIJ20,FFT10,FFT11,FFT20,
     2            NFDIJ1,NFDIJ2,NFDIJ3)
C
C
C     ----- NUCLEAR ATTRACTION -----
C
      DUM = PI212*AA1
      IF(MFIRST) THEN
         DO 50 I = 1,NIJ10
   50    DIJ10(I) = DIJ10(I) * DUM
      END IF
      IF(MSECND) THEN
         DO 51 I = 1,NIJ11
   51    DIJ11(I) = DIJ11(I) * DUM
         DO 52 I = 1,NIJ20
   52    DIJ20(I) = DIJ20(I) * DUM
      END IF
C
      AAX = AA*AX
      AAY = AA*AY
      AAZ = AA*AZ
C
c     ICMAX = NAT
c     IF(ISEPS) ICMAX=NAT+NPS
C
      DO 800 IC = 1,NTS
C     COSMO PARTIAL CHARGES
      ZNUC = -QT(IC)
      CX   = XCTS(IC)
      CY   = YCTS(IC)
      CZ   = ZCTS(IC)
C
      XX = AA*((AX-CX)**2+(AY-CY)**2+(AZ-CZ)**2)
      IF(NROOTS .LE. 3) CALL RT123
      IF(NROOTS .EQ. 4) CALL ROOT4
      MM = 0
      DO 740 K = 1,NROOTS
      UU = AA*U(K)
      WW = W(K)*ZNUC
      TT = ONE/(AA+UU)
      T =  SQRT(TT)
      X0 = (AAX+UU*CX)*TT
      Y0 = (AAY+UU*CY)*TT
      Z0 = (AAZ+UU*CZ)*TT
      IN = -4+MM
      DO 720 I = 1,LIT
      IN = IN+4
      NI = I
      DO 720 J = 1,LJT
      JN = IN+J
      NJ = J
      CALL DERXYZ(0)
      XIN(JN) = XINT
      YIN(JN) = YINT
      ZIN(JN) = ZINT*WW
  720 CONTINUE
  740 MM = MM+20
C
      IF(MFIRST) THEN
         DO 60 I = 1,NIJ10
            N  = NADIJ1(I)
            NX = NPIJX1(I)
            NY = NPIJY1(I)
            NZ = NPIJZ1(I)
            DUM = ZERO
            MM  = 0
            DO 61 K = 1,NROOTS
               DUM = DUM+XIN(NX+MM)*YIN(NY+MM)*ZIN(NZ+MM)
               MM = MM+20
   61       CONTINUE
            GIJ10(N) = GIJ10(N) + DUM*DIJ10(I)
   60    CONTINUE
      END IF
C
      IF(MSECND) THEN
         IF(JJ .GT. II ) GO TO 790
         DO 62 I =1,NIJ11
            N  = NADIJ2(I)
            NX = NPIJX2(I)
            NY = NPIJY2(I)
            NZ = NPIJZ2(I)
            DUM = ZERO
            MM  = 0
            DO K = 1,NROOTS
               DUM = DUM+XIN(NX+MM)*YIN(NY+MM)*ZIN(NZ+MM)
               MM = MM+20
            ENDDO
            GIJ11(N) = GIJ11(N) + DUM*DIJ11(I)
   62    CONTINUE
  790    CONTINUE
         DO 64 I =1,NIJ20
            N  = NADIJ3(I)
            NX = NPIJX3(I)
            NY = NPIJY3(I)
            NZ = NPIJZ3(I)
            DUM = ZERO
            MM  = 0
            DO 65 K = 1,NROOTS
               DUM = DUM+XIN(NX+MM)*YIN(NY+MM)*ZIN(NZ+MM)
               MM = MM+20
   65       CONTINUE
            GIJ20(N) = GIJ20(N) + DUM*DIJ20(I)
   64    CONTINUE
      END IF
  800 CONTINUE
C
C     ----- END OF PRIMITIVE LOOPS -----
C
  820 CONTINUE
  840 CONTINUE
C
C     ----- FORM INTEGRALS OVER DERIVATIVES -----
C
      DO 70 I = 1,9
   70 DS(I,2) = ZERO
      CALL TVFSIJ(GIJ10,GIJ20,DF,DS(1,1),DOE,NFDIJ1,NFDIJ3,DAB,FOE)
      IF(.NOT.MCPHF) GO TO 846
      N = -3
      DO 75 J = MINJ,MAXJ
      NJ = LOCJ + J
      DO 75 I = MINI,MAXI
      NI = LOCI + I
      IF(NI .GE. NJ) GO TO 842
      NIJ = NJ * (NJ-1)/2 + NI
      GO TO 844
  842 NIJ = NI * (NI-1)/2 + NJ
  844 N = N + 3
      FAC = ONE
      IF(NI .EQ. NJ) FAC = TWO
      FD(IFDX+NIJ) = FD(IFDX+NIJ) + FAC * DOE(N+1)
      FD(IFDY+NIJ) = FD(IFDY+NIJ) + FAC * DOE(N+2)
      FD(IFDZ+NIJ) = FD(IFDZ+NIJ) + FAC * DOE(N+3)
   75 CONTINUE
  846 CONTINUE
      IAT = KATOM(II)
      JAT = KATOM(JJ)
      FAC = ONE
      IF(IAT .EQ. JAT .AND. II .NE. JJ) FAC = TWO
      IF(BFDD) CALL ADDFOE(FDD,IAT,IAT,FOE,ONE,TWO)
C
      IF(.NOT.MSECND .OR. JJ .GT. II) GO TO 850
         CALL TVFIFJ(GIJ11,DS(1,2),XIN,YIN,ZIN,NW,NFDIJ2,DAB,FOE)
         IF(BFDD) CALL ADDFOE(FDD,IAT,JAT,FOE,FAC,TWO)
  850 CONTINUE
C
C     ----- CALCULATE CONTRIBUTION TO GRADIENT AND FORCE CONSTANTS -----
C
      IAT = KATOM(II)
      JAT = KATOM(JJ)
      FAC = ONE
      IF(IAT .EQ. JAT .AND. II .NE. JJ) FAC = TWO
      IF(JAT .GT. IAT) GO TO 860
      IIAT = IAT*(IAT-1)/2 + IAT
      IJAT = IAT*(IAT-1)/2 + JAT
      GO TO 865
  860 IIAT = IAT*(IAT-1)/2 + IAT
      IJAT = JAT*(JAT-1)/2 + IAT
  865 CONTINUE
      DO 80 I = 1,3
   80 EG(I, IAT) = EG(I, IAT) + DF(I)
      IF(MSECND) THEN
         DO 81 I = 1,9
            EH(I,IIAT) = EH(I,IIAT) + DS(I,1)
            EH(I,IJAT) = EH(I,IJAT) + DS(I,2) * FAC
   81    CONTINUE
         IF(DBUG) WRITE(IW,9028) II,JJ
         IF(DBUG) CALL HSSPRT(NAT,EG,EH)
      END IF
 1040 CONTINUE
 1060 CONTINUE
C
C     ----- END OF SHELL LOOPS -----
C
      IF (GOPARR) THEN
        IF(NXT) CALL DDI_DLBRESET
        CALL DDI_GSUMF(1605,FD,3*NAT*L2)
        CALL DDI_GSUMF(1624,EG,NAT*3)
        CALL DDI_GSUMF(1625,EH,9*(NAT*(NAT+1))/2)
      END IF
C
      IF(MSECND) THEN
         DO 1080 IAT = 1,NAT
            IIAT = ( IAT * (IAT-1) ) / 2 + IAT
            FAC  = ( EH(2,IIAT) + EH(4,IIAT) ) / TWO
            EH(2,IIAT) = FAC
            EH(4,IIAT) = FAC
            FAC  = ( EH(3,IIAT) + EH(7,IIAT) ) / TWO
            EH(3,IIAT) = FAC
            EH(7,IIAT) = FAC
            FAC  = ( EH(6,IIAT) + EH(8,IIAT) ) / TWO
            EH(6,IIAT) = FAC
            EH(8,IIAT) = FAC
 1080    CONTINUE
      END IF
C
C     ----- PRINTING SECTION -----
C
      IF(OUT) THEN
         NEGH = 3*NAT+9*(NAT*(NAT+1))/2
         IF(GOPARR) CALL DDI_GSUMF(1624,EG,NEGH)
         IF(MASWRK) THEN
            WRITE(IW,9058)
            CALL HSSPRT(NAT,EG,EH)
         END IF
         IF(GOPARR) CALL DSCAL(NEGH,ONE/NPROC,EG,1)
C
         IF(MCPHF) THEN
            IF(MASWRK) WRITE(IW,9068)
            DO 1100 IC=1,NAT
               IF (MASWRK) WRITE(IW,9048) IC
               IFDX=3*L2*(IC-1)+1
               IFDY=IFDX+L2
               IFDZ=IFDY+L2
               CALL PRTRI(FD(IFDX),NUM)
               CALL PRTRI(FD(IFDY),NUM)
               CALL PRTRI(FD(IFDZ),NUM)
 1100       CONTINUE
         END IF
      END IF
C
C     ----- SET UP RESTART -----
C
      NREST = 4
      IST = 1
      JST = 1
      KST = 1
      LST = 1
      RETURN
C
 9008 FORMAT(/5X,50(1H-)/
     *        5X,'1E- INTEGRAL CONTRIBUTIONS TO GRADIENT AND HESSIAN'/
     *        5X,50(1H-))
 9028 FORMAT(' SHELLS II,JJ = ',2I5)
 9048 FORMAT(1X,'CONTRIBUTION TO FOCK DERIVATIVE MATRIX, ATOM= ',I5)
 9058 FORMAT(/1X,'GRADIENT AND HESSIAN AFTER -TVDSPD-')
 9068 FORMAT(/1X,'FOCK DERIVATIVES AFTER -TVDSPD-')
      END
C
C*MODULE PCMHSS  *DECK TRANINVHSS
C>    @brief Second order derivatiev for PCM 
C>
C>    @details  Translation invariance for hessian
C>
C>    @author Hiroya Nakata
C>    - Apr, 2013- Subroutine writte
C>   @param   NCNTR   center of atom
C>   @param   NAT     number of atom
C>   @param   FEQ     derivatives
C>   @param   FDQ     Work for Fock derivative elements
C>   @param   EHQ     Work for Hessian
C>   @param   DSOLQM  gradient
C>   @param   FD      Fock derivative 
C>   @param   EH      Hessian
C>   @param   FEQWRK  Derivative elements
C>   @param   L2      triangular matrix for orbital
C>   @param   nat2    triangular matrix for atom
      SUBROUTINE TRANINVHSS(NCNTR,NAT,FEQ,FDQ,EHQ,DSOLQM,FD,EH,
     *           FEQWRK,L2,nat2)

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
c     COMMON /FMOINF/ NFG,NLAYER,NATFMO,NBDFG,NAOTYP,NBODY,NSEGM
C
      DIMENSION FEQ(3,NAT),FDQ(L2,3,nat),EHQ(9,*)
      DIMENSION DSOLQM(3,NAT),FD(L2,3,nat),EH(9,*)
      DIMENSION FEQWRK(9,*)
c     DIMENSION FCMWRK(natfmo*3,*)
      PARAMETER (ONE=1.0D+00)
C
      DO JAT=1,NAT
       DO IXYZ=1,3
        IF(NCNTR.NE.JAT) THEN
C       for gradient
         DSOLQM(IXYZ,JAT)  =DSOLQM(IXYZ,JAT)   + FEQ(IXYZ,JAT)
         if(NCNTR.ne.0)
     *      DSOLQM(IXYZ,NCNTR)=DSOLQM(IXYZ,NCNTR) - FEQ(IXYZ,JAT)
C       for hessian
         CALL DAXPY(L2, ONE,FDQ(1,IXYZ,JAT),1,FD(1,IXYZ,JAT),1)
         if(NCNTR.ne.0)
     *      CALL DAXPY(L2,-ONE,FDQ(1,IXYZ,JAT),1,FD(1,IXYZ,NCNTR),1)
        END IF
       END DO
      END DO

C     W^{ab} 
C       1. Normal terms
      DO 700 JATM=1,NAT
        IF(JATM.EQ.NCNTR) GOTO 700
        DO 710 KATM=1,JATM
         IF(KATM.EQ.NCNTR) GOTO 710
         JKAT=(JATM*JATM-JATM)/2 +  KATM
         CALL DAXPY(9,ONE,EHQ(1,JKAT),1,EH(1,JKAT),1)
 710    CONTINUE
 700  CONTINUE

      if(NCNTR.eq.0) return
c     return
 
C     2. Translation invariance for term 5
      DO 720 JATM=1,NAT
        IF (JATM.EQ.NCNTR) GOTO 720
        DO 730  KATM=1,NAT
         IF (KATM.EQ.NCNTR) GOTO 730
         IF(JATM.GT.NCNTR) THEN
           NCNT2 = (JATM *JATM  - JATM )/ 2 + NCNTR
           IF(JATM.GT.KATM) THEN
             JKAT=(JATM*JATM-JATM)/2 +  KATM
             FEQWRK(1,NCNT2)=FEQWRK(1,NCNT2)-EHQ(1,JKAT)
             FEQWRK(2,NCNT2)=FEQWRK(2,NCNT2)-EHQ(2,JKAT)
             FEQWRK(3,NCNT2)=FEQWRK(3,NCNT2)-EHQ(3,JKAT)
             FEQWRK(4,NCNT2)=FEQWRK(4,NCNT2)-EHQ(4,JKAT)
             FEQWRK(5,NCNT2)=FEQWRK(5,NCNT2)-EHQ(5,JKAT)
             FEQWRK(6,NCNT2)=FEQWRK(6,NCNT2)-EHQ(6,JKAT)
             FEQWRK(7,NCNT2)=FEQWRK(7,NCNT2)-EHQ(7,JKAT)
             FEQWRK(8,NCNT2)=FEQWRK(8,NCNT2)-EHQ(8,JKAT)
             FEQWRK(9,NCNT2)=FEQWRK(9,NCNT2)-EHQ(9,JKAT)
           END IF
           IF(KATM.GE.JATM) THEN
             JKAT=(KATM*KATM-KATM)/2 + JATM
             FEQWRK(1,NCNT2)=FEQWRK(1,NCNT2)-EHQ(1,JKAT)
             FEQWRK(2,NCNT2)=FEQWRK(2,NCNT2)-EHQ(4,JKAT)
             FEQWRK(3,NCNT2)=FEQWRK(3,NCNT2)-EHQ(7,JKAT)
             FEQWRK(4,NCNT2)=FEQWRK(4,NCNT2)-EHQ(2,JKAT)
             FEQWRK(5,NCNT2)=FEQWRK(5,NCNT2)-EHQ(5,JKAT)
             FEQWRK(6,NCNT2)=FEQWRK(6,NCNT2)-EHQ(8,JKAT)
             FEQWRK(7,NCNT2)=FEQWRK(7,NCNT2)-EHQ(3,JKAT)
             FEQWRK(8,NCNT2)=FEQWRK(8,NCNT2)-EHQ(6,JKAT)
             FEQWRK(9,NCNT2)=FEQWRK(9,NCNT2)-EHQ(9,JKAT)
           END IF
         END IF
         IF(JATM.LT.NCNTR) THEN
           NCNT2 = (NCNTR*NCNTR - NCNTR)/ 2 + JATM 
           IF(KATM.GE.JATM) THEN
             JKAT=(KATM*KATM-KATM)/2 + JATM
             FEQWRK(1,NCNT2)=FEQWRK(1,NCNT2)-EHQ(1,JKAT)
             FEQWRK(2,NCNT2)=FEQWRK(2,NCNT2)-EHQ(2,JKAT)
             FEQWRK(3,NCNT2)=FEQWRK(3,NCNT2)-EHQ(3,JKAT)
             FEQWRK(4,NCNT2)=FEQWRK(4,NCNT2)-EHQ(4,JKAT)
             FEQWRK(5,NCNT2)=FEQWRK(5,NCNT2)-EHQ(5,JKAT)
             FEQWRK(6,NCNT2)=FEQWRK(6,NCNT2)-EHQ(6,JKAT)
             FEQWRK(7,NCNT2)=FEQWRK(7,NCNT2)-EHQ(7,JKAT)
             FEQWRK(8,NCNT2)=FEQWRK(8,NCNT2)-EHQ(8,JKAT)
             FEQWRK(9,NCNT2)=FEQWRK(9,NCNT2)-EHQ(9,JKAT)
           END IF
           IF(JATM.GT.KATM) THEN
             JKAT=(JATM*JATM-JATM)/2 + KATM
             FEQWRK(1,NCNT2)=FEQWRK(1,NCNT2)-EHQ(1,JKAT)
             FEQWRK(2,NCNT2)=FEQWRK(2,NCNT2)-EHQ(4,JKAT)
             FEQWRK(3,NCNT2)=FEQWRK(3,NCNT2)-EHQ(7,JKAT)
             FEQWRK(4,NCNT2)=FEQWRK(4,NCNT2)-EHQ(2,JKAT)
             FEQWRK(5,NCNT2)=FEQWRK(5,NCNT2)-EHQ(5,JKAT)
             FEQWRK(6,NCNT2)=FEQWRK(6,NCNT2)-EHQ(8,JKAT)
             FEQWRK(7,NCNT2)=FEQWRK(7,NCNT2)-EHQ(3,JKAT)
             FEQWRK(8,NCNT2)=FEQWRK(8,NCNT2)-EHQ(6,JKAT)
             FEQWRK(9,NCNT2)=FEQWRK(9,NCNT2)-EHQ(9,JKAT)
           END IF
         END IF
 730    CONTINUE
 720  CONTINUE
      NCNT2 = (NCNTR *NCNTR -NCNTR)/2 + NCNTR
      DO JATM=1,NAT
        IF (JATM.NE.NCNTR) THEN
          IF(JATM.GT.NCNTR) THEN
           JNAT2 = (JATM*JATM-JATM)/2   +NCNTR
           FEQWRK(1,NCNT2)=FEQWRK(1,NCNT2)-FEQWRK(1,JNAT2)
           FEQWRK(2,NCNT2)=FEQWRK(2,NCNT2)-FEQWRK(4,JNAT2)
           FEQWRK(3,NCNT2)=FEQWRK(3,NCNT2)-FEQWRK(7,JNAT2)
           FEQWRK(4,NCNT2)=FEQWRK(4,NCNT2)-FEQWRK(2,JNAT2)
           FEQWRK(5,NCNT2)=FEQWRK(5,NCNT2)-FEQWRK(5,JNAT2)
           FEQWRK(6,NCNT2)=FEQWRK(6,NCNT2)-FEQWRK(8,JNAT2)
           FEQWRK(7,NCNT2)=FEQWRK(7,NCNT2)-FEQWRK(3,JNAT2)
           FEQWRK(8,NCNT2)=FEQWRK(8,NCNT2)-FEQWRK(6,JNAT2)
           FEQWRK(9,NCNT2)=FEQWRK(9,NCNT2)-FEQWRK(9,JNAT2)
          ELSE IF(NCNTR.GT.JATM) THEN
           JNAT2 = (NCNTR*NCNTR-NCNTR)/2+JATM 
           FEQWRK(1,NCNT2)=FEQWRK(1,NCNT2)-FEQWRK(1,JNAT2)
           FEQWRK(2,NCNT2)=FEQWRK(2,NCNT2)-FEQWRK(2,JNAT2)
           FEQWRK(3,NCNT2)=FEQWRK(3,NCNT2)-FEQWRK(3,JNAT2)
           FEQWRK(4,NCNT2)=FEQWRK(4,NCNT2)-FEQWRK(4,JNAT2)
           FEQWRK(5,NCNT2)=FEQWRK(5,NCNT2)-FEQWRK(5,JNAT2)
           FEQWRK(6,NCNT2)=FEQWRK(6,NCNT2)-FEQWRK(6,JNAT2)
           FEQWRK(7,NCNT2)=FEQWRK(7,NCNT2)-FEQWRK(7,JNAT2)
           FEQWRK(8,NCNT2)=FEQWRK(8,NCNT2)-FEQWRK(8,JNAT2)
           FEQWRK(9,NCNT2)=FEQWRK(9,NCNT2)-FEQWRK(9,JNAT2)
          END IF
        END IF 
      END DO
      CALL DAXPY(NAT2*9,ONE,FEQWRK,1,EH,1)
C
      RETURN
      END
C
C*MODULE PCMHSS  *DECK PCMNUCHSS
C>    @brief Second order derivatiev for PCM 
C>
C>    @details  Nuclear interaction derivative 
C>
C>    @author Hiroya Nakata
C>    - Apr, 2013- Subroutine writte
C>   @param   FEQ    derivative element
C>   @param   EHQ    Work for Hessian 
C>   @param   XCTS    x coordinate for teserae
C>   @param   YCTS    y coordinate for teserae
C>   @param   ZCTS    z coordinate for teserae
C>   @param   NTS     number of teserae
C>   @param   QT      total charge
C>   @param   IAT     number of atom
      SUBROUTINE PCMNUCHSS(FEQ,EHQ,XCTS,YCTS,ZCTS,NTS,QT,IAT)
      use mx_limits, only: mxatm
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      DIMENSION FEQ(3,*),EHQ(9,*),QT(*)
      DIMENSION XCTS(*),YCTS(*),ZCTS(*)
C
C
c     COMMON /FMOINF/ NFG,NLAYER,NATFMO,NBDFG,NAOTYP,NBODY,NSEGM
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)

      DO ITS = 1, NTS
       DO JAT = 1, NAT
        IF(JAT.NE.IAT) THEN
C         - OTHER NUCLEI -
          XI = XCTS(ITS)
          YI = YCTS(ITS)
          ZI = ZCTS(ITS)
          XN = C(1,JAT)
          YN = C(2,JAT)
          ZN = C(3,JAT)
          R2 = (XN-XI)**2+(YN-YI)**2+(ZN-ZI)**2
          R  = SQRT(R2)
          R3 = R*R2
          Z  = ZAN(JAT)
C
          DUM = Z*QT(ITS)/R3
c         gradient
          FEQ(1,JAT)=FEQ(1,JAT) - DUM*(XN-XI) 
          FEQ(2,JAT)=FEQ(2,JAT) - DUM*(YN-YI)
          FEQ(3,JAT)=FEQ(3,JAT) - DUM*(ZN-ZI)
c         hessian
          SCAL1=3.0D+00
          SCAL2=1.0D+00
C
          DUMXX=SCAL1*DUM*(XN-XI)*(XN-XI)/R2 - SCAL2*DUM 
          DUMXY=SCAL1*DUM*(XN-XI)*(YN-YI)/R2 
          DUMXZ=SCAL1*DUM*(XN-XI)*(ZN-ZI)/R2
          DUMYX=SCAL1*DUM*(YN-YI)*(XN-XI)/R2 
          DUMYY=SCAL1*DUM*(YN-YI)*(YN-YI)/R2 - SCAL2*DUM 
          DUMYZ=SCAL1*DUM*(YN-YI)*(ZN-ZI)/R2
          DUMZX=SCAL1*DUM*(ZN-ZI)*(XN-XI)/R2 
          DUMZY=SCAL1*DUM*(ZN-ZI)*(YN-YI)/R2
          DUMZZ=SCAL1*DUM*(ZN-ZI)*(ZN-ZI)/R2 - SCAL2*DUM 
C
          JAT2 = (JAT *JAT - JAT)/2 + JAT
          EHQ(1,JAT2) = EHQ(1,JAT2) + DUMXX
          EHQ(2,JAT2) = EHQ(2,JAT2) + DUMXY
          EHQ(3,JAT2) = EHQ(3,JAT2) + DUMXZ
          EHQ(4,JAT2) = EHQ(4,JAT2) + DUMYX
          EHQ(5,JAT2) = EHQ(5,JAT2) + DUMYY
          EHQ(6,JAT2) = EHQ(6,JAT2) + DUMYZ
          EHQ(7,JAT2) = EHQ(7,JAT2) + DUMZX
          EHQ(8,JAT2) = EHQ(8,JAT2) + DUMZY
          EHQ(9,JAT2) = EHQ(9,JAT2) + DUMZZ
          
        END IF
       END DO
      END DO
      RETURN
      END
C*MODULE PCMHSS  *DECK FNUCNT
C>    @brief Second order derivatiev for PCM 
C>
C>    @details  Nuclear interaction derivative for each tessere
C>
C>    @author Hiroya Nakata
C>    - Apr, 2013- Subroutine writte
C>   @param   FEQ   derivative element
C>   @param   XCTS  x coordinate for teserae
C>   @param   YCTS  y coordinate for teserae
C>   @param   ZCTS  z coordinate for teserae
C>   @param   NTS   number of tesserae
C>   @param   ISPHE indice for atom
C>
      SUBROUTINE FNUCNT(FEQ,XCTS,YCTS,ZCTS,NTS,ISPHE)
      use mx_limits, only: mxatm
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      DIMENSION FEQ(3,nat,*)
      DIMENSION XCTS(*),YCTS(*),ZCTS(*)
      DIMENSION ISPHE(*)
C
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
C
C
      DO ITS = 1, NTS
       IAT = ISPHE(ITS)
       DO JAT = 1, NAT
        IF(JAT.NE.IAT) THEN
C         - OTHER NUCLEI -
          XI = XCTS(ITS)
          YI = YCTS(ITS)
          ZI = ZCTS(ITS)
          XN = C(1,JAT)
          YN = C(2,JAT)
          ZN = C(3,JAT)
          R2 = (XN-XI)**2+(YN-YI)**2+(ZN-ZI)**2
          R  = SQRT(R2)
          R3 = R*R2
          Z  = ZAN(JAT)
C
          DUM = Z/R3
c         gradient
          FEQ(1,JAT,ITS)=FEQ(1,JAT,ITS) - DUM*(XN-XI) 
          FEQ(2,JAT,ITS)=FEQ(2,JAT,ITS) - DUM*(YN-YI)
          FEQ(3,JAT,ITS)=FEQ(3,JAT,ITS) - DUM*(ZN-ZI)
          
        END IF
       END DO
C      Translation invraiance
       CALL VCLR(FEQ(1,IAT,ITS),1,3)
       DO JAT=1,NAT
        IF(JAT.NE.IAT) THEN
         FEQ(1,IAT,ITS)=FEQ(1,IAT,ITS)-FEQ(1,JAT,ITS)
         FEQ(2,IAT,ITS)=FEQ(2,IAT,ITS)-FEQ(2,JAT,ITS)
         FEQ(3,IAT,ITS)=FEQ(3,IAT,ITS)-FEQ(3,JAT,ITS)
        END IF   
       END DO
      END DO
      RETURN
      END
C
C*MODULE PCMHSS  *DECK PCMCMPHSS
C>    @brief Second order derivatiev for PCM 
C>
C>    @details  Sum up hessian contribution
C>
C>    @author Hiroya Nakata
C>    - Apr, 2013- Subroutine writte
C>   @param   FCM   force constant
C>   @param   NTS   number of tesserae
C>   @param   QA    charge
C>   @param   FEQ   derivative element
C>   @param   NAT   number of atom
C>   @param   ICURATM indice for atom
C>   @param   ICURXYZ indice for xyz
C>
      SUBROUTINE PCMCMPHSS(FCM,NTS,QA,FEQ,NAT,ICURATM,ICURXYZ)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION FCM(nat*3,nat*3),QA(*),FEQ(3,nat,*)
C
      jxyz   = (ICURATM-1) * 3 + ICURXYZ
C
      DO IATM=1,NAT
        DO IXYZ=1,3
          iixyz= (IATM-1) * 3 + IXYZ
          val  = 0.0D+00
          DO ITS=1,NTS
           FCM(iixyz,jxyz)=FCM(iixyz,jxyz)+ QA(ITS)*FEQ(ixyz,IATM,ITS) 
           val = val + QA(ITS)*FEQ(ixyz,IATM,ITS)
          END DO
c         write(6,*) "www pcmcmphss =",jxyz,iixyz,val
        END DO
      END DO
C
C
      RETURN
      END
C*MODULE PCMHSS  *DECK pcmaddhss
C
C>    @brief Second order derivatiev for PCM 
C>
C>    @details  add hessian to fc matrix 
C>
C>    @author Hiroya Nakata
C>    - Apr, 2013- Subroutine writte
C>   @param EH   Hessian
C>   @param FCM   Force constant matrix
C>   @param NAT   Number of atoms
      subroutine pcmaddhss(EH,FCM,NAT)
      use mx_limits, only: mxatm
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION EH(9,*),FCM(nat*3,nat*3)
C
      DO IAT=1,NAT
        iix=(IAT-1)*3+1
        iiy=(IAT-1)*3+2
        iiz=(IAT-1)*3+3
        DO JAT=1,IAT
         jjx=(JAT-1)*3+1
         jjy=(JAT-1)*3+2
         jjz=(JAT-1)*3+3
C
         IIAT = IAT*(IAT-1)/2 + JAT
         EH(1,IIAT) = EH(1,IIAT) + FCM(iix,jjx)
         EH(2,IIAT) = EH(2,IIAT) + FCM(iix,jjy)
         EH(3,IIAT) = EH(3,IIAT) + FCM(iix,jjz)
         EH(4,IIAT) = EH(4,IIAT) + FCM(iiy,jjx)
         EH(5,IIAT) = EH(5,IIAT) + FCM(iiy,jjy)
         EH(6,IIAT) = EH(6,IIAT) + FCM(iiy,jjz)
         EH(7,IIAT) = EH(7,IIAT) + FCM(iiz,jjx)
         EH(8,IIAT) = EH(8,IIAT) + FCM(iiz,jjy)
         EH(9,IIAT) = EH(9,IIAT) + FCM(iiz,jjz)
        END DO
      END DO
C
      RETURN
      END
C
C
c     SUBROUTINE CHKGRD(DEN,VPOT,NTS,QA,XCTS,YCTS,ZCTS)
c     IMPLICIT DOUBLE PRECISION(A-H,O-Z)
c     DIMENSION DEN(*),VPOT(*),QA(*),XCTS(*),YCTS(*),ZCTS(*)
c     COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
c    *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
c
c     val=TRACEP(den,vpot,num)
c     DO ITS = 1, NTS
c      DO JAT = 1, NAT
c         XI  = XCTS(ITS)
c         YI  = YCTS(ITS)
c         ZI  = ZCTS(ITS)
c         XN  = C(1,JAT)
c         YN  = C(2,JAT)
c         ZN  = C(3,JAT)
c         R2  = (XN-XI)**2+(YN-YI)**2+(ZN-ZI)**2
c         R   = SQRT(R2)
c         Z   = ZAN(JAT)
C
c         val = val + Z*QA(ITS)/R
c      END DO
c     END DO
c     write(6,'("val=",F10.7)') val
c
c     RETURN
c     END
C
C
C*MODULE PCMHSS  *DECK AOORBHSSCP
C>    @brief Second order derivatiev for PCM 
C>
C>    @details  add  orbital hessian to CPHF
C>
C>    @author Hiroya Nakata
C>    - Apr, 2013- Subroutine writte
C>   @param PMN   Density matrix
C>   @param FCM   Fock derivative matrix
C>   @param  L1   Number of orbital
C>   @param NAT   Number of triangular matrix
C>   @param NXYZ  Number of xyz elements
      SUBROUTINE AOORBHSSCP(PMN,FMN,L1,L2,NXYZ) 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      COMMON /FMCOM / X(1)
      COMMON /PCMDIM/ MXSP,MXTS,MEMPCM1,MEMPCM2,NTS
      COMMON /PCMPNT/ LXYZRE,LSSFE,LLIST,LALPCM,LRINPCM,LINAPCM,LINFPCM,
     *                LINIPCM,LMEPCM,LAXYZCT,LXYZCT2,LNVERT,LQSN,LQSE,
     *                LQFS,LQIND,LISPHE,LVAD,LQOR,LVECMUL,LAIPRJ,LFIPRJ,
     *                LPEL,LPCMCDR,LDAI,LIDDAI,LQSND,LQSED,LEPSPCM
C
      LOGICAL GOPARR, DSKWRK, MASWRK
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
C
      DIMENSION PMN(NXYZ,*),FMN(NXYZ,*)
C
      if(L1.EQ.0) write(6,*) "wwwww"
C
      CALL VALFM(LOADFM)
      NDASC = 0
      LFPNT = LOADFM
      NTSPAR=(NTS-1)/nproc + 1
      CALL INDQPROD(0,NTSPAR,LFPNT,NDASC,DUM,DUM,DUM,0,0)
      LQTMP   = LFPNT   + NDASC  + 1
      LEPOT   = LQTMP   + nts    
      LWRK    = LEPOT   + nts    
      LAST    = LWRK    + L2 
      NEEDW   = LAST    - LOADFM    - 1
      CALL GETFM(NEEDW)
C
      CALL ORBHSSCPX(nts,X(LAXYZCT),X(LAXYZCT+MXTS),X(LAXYZCT+MXTS*2),
     *     X(LQTMP),X(LEPOT),X(LWRK),NTSPAR,LFPNT,NDASC,L2,NXYZ,
     *     PMN,FMN)
C
      CALL RETFM(NEEDW)
C
      RETURN
      END
C
C
C*MODULE PCMHSS  *DECK ORBHSSCPX
C>    @brief Second order derivatiev for PCM 
C>
C>    @details  construct orbital hessian to CPHF
C>
C>    @author Hiroya Nakata
C>    - Apr, 2013- Subroutine writte
C>   @param nts    number of tesserae
C>   @param XCTS    X coordinate 
C>   @param YCTS    Y coordinate
C>   @param ZCTS    Z coordinate
C>   @param QTMP    temporary charge
C>   @param EPOT    electronic potential
C>   @param WRK     work
C>   @param NTSPAR  ntesserae for parallel
C>   @param LFPNT   start point
C>   @param NDASC   size for PCM memory 
C>   @param L2      number of triangular matrix
C>   @param NXYZ    number of xyz coordinate
C>   @param PMN     Density like  matrix
C>   @param FMN     Fock    like  matrix
      SUBROUTINE ORBHSSCPX(nts,XCTS,YCTS,ZCTS,QTMP,EPOT,WRK,
     *           NTSPAR,LFPNT,NDASC,L2,NXYZ,PMN,FMN)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION XCTS(*),YCTS(*),ZCTS(*),QTMP(*),EPOT(*),WRK(*)
      DIMENSION PMN(NXYZ,*),FMN(NXYZ,*)

      LOGICAL GOPARR,DSKWRK,MASWRK
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
C
      PARAMETER (ONE=1.0D+00)

c     call dacopy(nts,one,QTMP,1)
      do ixyz=1,nxyz
       do iloop=1,L2
          WRK(iloop)=PMN(ixyz,iloop)
       end do
       CALL VCLR(EPOT,1,nts)
       call EPOTEN(-one,XCTS,YCTS,ZCTS,EPOT,wrk,NTS,L2)
       if(GOPARR) call ddi_gsumf(2418,EPOT,NTS)
c
       ITER0 = 1
       CALL INDQPROD(ITER0,NTSPAR,LFPNT,NDASC,EPOT,QTMP,QTMP,1,0)
C      Fock derivative terms
       CALL PCMPOT(wrk,wrk,NTS,QTMP,XCTS,L2,1,.false.)
c
       do iloop=1,L2
          FMN(ixyz,iloop)=FMN(ixyz,iloop) + WRK(iloop) 
       end do
      end do
      RETURN
      END
C
C*MODULE PCMHSS  *DECK QCHESS
C>    @brief Second order derivatiev for PCM 
C>
C>    @details  Derivative of  C
C>
C>    @author Hiroya Nakata
C>
C>   @pram  NATM     number of atom
C>   @pram  IMODE    mode
C>   @pram  FEQWRK   derivative element
C>   @pram  QT       total charge
C>   @pram  VECTMP   potential
C>   @pram  XCTS     x coordinate for teserae
C>   @pram  YCTS     y coordinate for teserae
C>   @pram  ZCTS     z coordinate for teserae
C>   @pram  AS       area surface charge
C>   @pram  QSN      nuclear charge
C>   @pram  QSE      electronic charge
C>   @pram  QSND     ASC for het
C>   @pram  QSED     ASC for het
C>   @pram  Q_FS     ASC
C>   @pram  Q_IND    ASC for induction
C>   @pram  Q2       ASC
C>   @pram  XE       Solute geometry informtion  
C>   @pram  YE       Solute geometry informtion
C>   @pram  ZE       Solute geometry informtion
C>   @pram  RE       Solute geometry informtion
C>   @pram  ISPHE    indice for sphere
C>   @pram  DAI      parameter for PCM
C>   @pram  IDDAI    indice for DAI
C>   @pram  QSETMP   FMO  charge
C>   @pram  QSEMP    MP charge
C>   @pram  QSEXY    PCM charge
C>   @pram  QSEMN    PCM charge
C>   @pram  DSOLQM   PCM charge
C>   @pram  LFMOPCM1 logic for PCM
C>   @pram  IAGLOB   global atom indice
C>   @pram  FCM      force constnt matrix
C>   @pram  FEQ      derivative element
      SUBROUTINE QCHESS(NATM,IMODE,FEQWRK,QT,VECTMP,
     *     XCTS,YCTS,ZCTS,AS,QSN,QSE,QSND,QSED,Q_FS,Q_IND,Q2,
     *     XE,YE,ZE,RE,ISPHE,DAI,IDDAI,QSETMP,QSEMP,QSEXY,QSEMN,
     *     DSOLQM,LFMOPCM1,IAGLOB,FCM,FEQ,FMOC,jmode)
      use mx_limits, only: mxatm

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C
      LOGICAL   LFMOPCM1
C
      COMMON /FMCOM / X(1)
      COMMON /FMOINF/ NFG,NLAYER,NATFMO,NBDFG,NAOTYP,NBODY,NSEGM
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON /MACHIN/ NWDVAR,MAXFM,MAXSM,LIMFM,LIMSM
      COMMON /PCMDIM/ MXSP,MXTS,MEMPCM1,MEMPCM2,NTS
      COMMON /PCMOPT/ RABI,RASC,REFPOL,THRSLS,DENSLS,WB,WA,ETA2,GD,EVAC,
     *                RHOW,PM,AREATL,AREAKP,BONDRY,OMEGA,RET,FRO,EPSINF,
     *                EPS,DR,RSOLV,VMOL,TCE,STEN,DSTEN,CMF,TABS,IDIRCT,
     *                IPCDER,IDP,ICOMP,IFIELD,ICAV,IDISP,IPRINT,IRETCAV,
     *                ICENT,IFAST,NEVAL,IEFPOL,KEEPSM,IMGABI,IMGASC,NADD
      COMMON /PCMPAR/ IPCM,NFT26,NFT27,IKREP,IEF,IP_F,NFMOPCM,IHET
      COMMON /PCMPNT/ LXYZRE,LSSFE,LLIST,LALPCM,LRINPCM,LINAPCM,LINFPCM,
     *                LINIPCM,LMEPCM,LAXYZCT,LXYZCT2,LNVERT,LQSN,LQSE,
     *                LQFS,LQIND,LISPHE,LVAD,LQOR,LVECMUL,LAIPRJ,LFIPRJ,
     *                LPEL,LPCMCDR,LDAI,LIDDAI,LQSND,LQSED,LEPSPCM
      COMMON /OUTPUT/ NPRINT,ITOL,ICUT,NORMF,NORMP,NOPK
C
      DIMENSION FEQ(3,NATM,*),QT(*),VECTMP(*),XCTS(*),YCTS(*),ZCTS(*)
      DIMENSION AS(*),QSN(*),QSE(*),QSND(*),QSED(*) 
      DIMENSION Q_FS(*),Q_IND(*),Q2(*),XE(*),YE(*),ZE(*),RE(*),ISPHE(*)
      DIMENSION DAI(3,20,MXTS),IDDAI(21,MXTS),QSETMP(*),QSEMP(*)
      DIMENSION QSEXY(*),QSEMN(*),DSOLQM(3,*)
      DIMENSION IAGLOB(*)
      DIMENSION D(2),FCM(natm*3,natm*3)
      DIMENSION FEQWRK(3,NATM)
      DIMENSION FMOC(3,*)
C
      if(imode.ne.0) write(6,*) "PCM HESS"
      NPRSAV = NPRINT
      NPRINT = 817
C
      if(nfg.ne.0) then
        if(jmode.eq.1) LFMOPCM1=.true. 
      end if
C
      NVIB   =  2 
      VIBSIZ =  0.000001D+00
c     VIBSIZ =  0.00001D+00
      D(1)   =  VIBSIZ
      D(2)   = -VIBSIZ
      LFIXPV = 0
      IF(IPCDER.EQ.3) LFIXPV = 3*20*MXTS
      LFIXPV2 = 0
      IF(IPCDER.EQ.3) LFIXPV2= 21*MXTS/NWDVAR + 1
C     Numerical derivatives of dC/da  :  (d^2 C/da db) 
c     write(6,*) "wwww begin numerical"
      DO 390 IVIB = 1,NVIB
        DO 380 IAT = 1,NAT
          DO 370 IXYZ = 1,3
C
           JXYZ       = IAT * 3 - 3 + IXYZ
           if(nfg.ne.0) then
            izat = iaglob(iat)
            fmoc(ixyz,izat) = fmoc(ixyz,izat) + D(IVIB)
           else
            C(IXYZ,IAT)= C(IXYZ,IAT) + D(IVIB)
           end if
C
           CALL MAKCVM(0)
           CALL TESIN(NATM,XE,YE,ZE,X(LALPCM),X(LRINPCM),X(LINIPCM),
     *                X(LMEPCM))
c          write(6,*) "enter pedram"
           CALL PEDRAM
           IF(IPCDER.EQ.3.AND.NFMOPCM.EQ.0) THEN
            CALL DAREAD(IDAF,IODA,DAI,LFIXPV,334,0)
            CALL DAREAD(IDAF,IODA,IDDAI,LFIXPV2,335,1)
           ENDIF
C
           CALL VCLR(FEQWRK,1,3*NAT)
           DUM=0
c          It seems that the call below will not work with EFP.
           CALL QCDERQ(NAT,NATM,NTS,0,FEQWRK,DUM,QT,VECTMP,XCTS,YCTS,
     *          ZCTS,AS,QSN,QSE,QSND,QSED,Q_FS,Q_IND,Q2,XE,YE,ZE,
     *          RE,ISPHE,DAI,IDDAI,QSETMP,QSEMP,QSEXY,QSEMN,
     *          FEQWRK,DUM,.false.,.false.,LFMOPCM1,IAGLOB,X(LEPSPCM),
     *          DUM,natfmo,nbody,.false.)
C
c          write(6,*) "end qcderq"
C

           IF(IVIB.EQ.1) CALL DCOPY(nat*3,FEQWRK,1,FCM(1,JXYZ),1)
           IF(IVIB.EQ.2) 
     *       CALL DAXPY(3*NAT,-1.0D+00,FEQWRK,1,FCM(1,JXYZ),1)

           if(nfg.ne.0) then
            izat = iaglob(iat)
            fmoc(ixyz,izat) = fmoc(ixyz,izat) - D(IVIB)
           else
            C(IXYZ,IAT)= C(IXYZ,IAT) - D(IVIB)
           end if
C
 370      CONTINUE
 380    CONTINUE
 390  CONTINUE
C
      DUM = 1.0D+00/(VIBSIZ*NVIB)
      DO 140 IAT = 1,NAT*3
         DO 130 JAT = 1,NAT*3
            FCM(IAT,JAT) = DUM*FCM(IAT,JAT)
  130    CONTINUE
  140 CONTINUE
      DO 590 IAT = 2,NAT*3
         DO 580 JAT = 1,IAT-1
            AVE = (FCM(IAT,JAT)+FCM(JAT,IAT))*0.5D+00
            FCM(JAT,IAT) = AVE
            FCM(IAT,JAT) = AVE
  580    CONTINUE
  590 CONTINUE
C
c     write(6,*) "wwww end numerical"
C
      CALL MAKCVM(0)
      CALL TESIN(NATM,XE,YE,ZE,X(LALPCM),X(LRINPCM),X(LINIPCM),
     *           X(LMEPCM))
      CALL PEDRAM
      IF(IPCDER.EQ.3.AND.NFMOPCM.EQ.0) THEN
       CALL DAREAD(IDAF,IODA,DAI,LFIXPV,334,0)
       CALL DAREAD(IDAF,IODA,IDDAI,LFIXPV2,335,1)
      ENDIF
c
c     write(6,*) "begin qcdfck"
c
c     if(nfg.eq.0)
       CALL QCDFCK(NAT,NATM,NTS,FEQ,QT,XCTS,YCTS,ZCTS,
     *           AS,QSN,QSE,QSND,QSED,Q_FS,Q_IND,
     *           ISPHE,DAI,IDDAI,QSETMP,QSEMP,QSEXY,QSEMN,
     *           DSOLQM,.false.,.false.,LFMOPCM1,IAGLOB,X(LEPSPCM))
C
C
c     write(6,*) "end qcdfck"
c     call prsq(fcm,nat*3,nat*3,nat*3)
      if(nfg.ne.0) then
        LFMOPCM1=.false. 
      end if
C
      NPRINT = NPRSAV
c
      RETURN
      END
C
C*MODULE PCMHSS  *DECK QCDFCK
C>    @brief Second order derivatiev for PCM 
C>
C>    @details  Derivative of  C
C>
C>    @author Hiroya Nakata
C>
C>   @param  NATA    number of atom
C>   @param  NATM    number of atom
C>   @param  NTS     number of tesserae
C>   @param  FEQ     derivative of elements
C>   @param  QT      total charge
C>   @param  XCTS    x coordinate for teserae 
C>   @param  YCTS    y coordinate for teserae
C>   @param  ZCTS    z coordinate for teserae
C>   @param  AS      area charge
C>   @param  QSN     nuclear charge
C>   @param  QSE     electronic charge
C>   @param  QSND    ASC for het
C>   @param  QSED    ASC for het
C>   @param  Q_FS    ASC
C>   @param  Q_IND   ASC for induction
C>   @param  ISPHE   ASC indice
C>   @param  DAI     PCM parameter
C>   @param  IDDAI   PCM parameter indice
C>   @param  QSETMP  ASC charge
C>   @param  QSEMP   ASC charge
C>   @param  QSEXY   ASC charge
C>   @param  QSEMN   ASC charge
C>   @param  DSOLQM  derivative
C>   @param  MP2PCM  MP2 logical
C>   @param  TDDFT   TDDFT logical
C>   @param  LFMOPCM1 FMO logical
C>   @param  IAGLOB   indice for global atm
C>   @param  EPSHET   epsilon
C>
      SUBROUTINE QCDFCK(NATA,NATM,NTS,FEQ,QT,XCTS,YCTS,ZCTS,
     *                  AS,QSN,QSE,QSND,QSED,Q_FS,Q_IND,
     *                  ISPHE,DAI,IDDAI,QSETMP,QSEMP,QSEXY,QSEMN,
     *                  DSOLQM,MP2PCM,TDDFT,LFMOPCM1,IAGLOB,EPSHET)
      use mx_limits, only: mxatm,mxfrg,mxfgpt
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C
      LOGICAL GOPARR, DSKWRK, MASWRK, TDDFT, MP2PCM
      LOGICAL LFMOPCM1
C
      DIMENSION FEQ(3,NATM,NTS)
      DIMENSION XCTS(1),YCTS(1),ZCTS(1),AS(1)
      DIMENSION QT(1),QSN(1),QSE(1),QSND(1),QSED(1),QSEMP(1),
     *          QSETMP(1),Q_FS(1),Q_IND(1),QSEXY(1),QSEMN(1)
      DIMENSION ISPHE(1),DSOLQM(3,*),DAI(3,20,*),IDDAI(21,*)
      DIMENSION IAGLOB(1)
      DIMENSION EPSHET(*)

      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /PCMOPT/ RABI,RASC,REFPOL,THRSLS,DENSLS,WB,WA,ETA2,GD,EVAC,
     *                RHOW,PM,AREATL,AREAKP,BONDRY,OMEGA,RET,FRO,EPSINF,
     *                EPS,DR,RSOLV,VMOL,TCE,STEN,DSTEN,CMF,TABS,IDIRCT,
     *                IPCDER,IDP,ICOMP,IFIELD,ICAV,IDISP,IPRINT,IRETCAV,
     *                ICENT,IFAST,NEVAL,IEFPOL,KEEPSM,IMGABI,IMGASC,NADD
      COMMON /PCMPAR/ IPCM,NFT26,NFT27,IKREP,IEF,IP_F,NFMOPCM,IHET
      COMMON /PCMPLY/ STOT,VOL,CCX,CCY,CCZ,RDIF,NESF,NESFP,NESFF,I_NESF,
     *                L_AST
C
      DATA ZERO, ONE /0.0D+00, 1.0D+00/
      DATA PT5/0.5D+00/
c     DATA FPI/12.56637062D+00/
C
C     CALCULATE -1/2*Q*DC/DA*Q
C     -- COMPUTE THE VECTOR A^(-1)*S*A^(-1)*Q2 STORE IN VECTMP
c     write(6,*) "IEF=",IEF,IPCDER
C
C     - INITIALIZE PARALLEL -
C
      IPCOUNT = ME - 1
C
      DO 100 IAT = 1, NATA
C        GO PARALLEL!
         IF(GOPARR) THEN
           IPCOUNT = IPCOUNT + 1
           IF(MOD(IPCOUNT,NPROC).NE.0) GO TO 100
         END IF
C
         IATGLB = IAT
         IF(NFMOPCM.GT.0) IATGLB = IAGLOB(IAT)
         IF (IAN(IAT).NE.1.AND.ABS(ZAN(IAT)-ONE).LT.1.0D-08) GOTO 100
C
C        TERM 2.A FORCE IMPOSED ON ASC(IAT) BY OTHER ASC
C
         FQQX  = ZERO
         FQQY  = ZERO
         FQQZ  = ZERO
         IF(NESFP.EQ.0)THEN
         DO ITS = 1, NTS
            IF(ISPHE(ITS).EQ.IATGLB) THEN
               XI = XCTS(ITS)
               YI = YCTS(ITS)
               ZI = ZCTS(ITS)
               DO JTS = 1, NTS
                  IF(ISPHE(JTS).NE.IATGLB) THEN
                     XJ = XCTS(JTS)
                     YJ = YCTS(JTS)
                     ZJ = ZCTS(JTS)
                     R2 = (XJ-XI)**2+(YJ-YI)**2+(ZJ-ZI)**2
                     R  = SQRT(R2)
                     R3 = R*R2
                     IF(IEF.EQ.10) THEN
                        EI=EPSHET(ISPHE(ITS))/(EPSHET(ISPHE(ITS))-ONE)
                        EJ=EPSHET(ISPHE(JTS))/(EPSHET(ISPHE(JTS))-ONE)
                        QTIJTS=QSE(JTS)+QSE(JTS)-QSED(JTS)
     *                        +QSN(JTS)+QSN(JTS)-QSND(JTS)
     *                                            +Q_FS(JTS)+Q_IND(JTS)
                        QTIITS=QSE(ITS)+QSE(ITS)-QSED(ITS)
     *                        +QSN(ITS)+QSN(ITS)-QSND(ITS)
     *                                            +Q_FS(ITS)+Q_IND(ITS)
                        IF (LFMOPCM1) THEN
                         QTDJTS=QSETMP(JTS)
                         QTDITS=QSETMP(ITS)
                         FT = (QTDITS*EJ*QTIJTS+QTIITS*EI*QTDJTS)/R3
                        ELSE
                         QTDJTS=QSED(JTS)+QSND(JTS)+Q_FS(JTS)+Q_IND(JTS)
                         QTDITS=QSED(ITS)+QSND(ITS)+Q_FS(ITS)+Q_IND(ITS)
                         FT = PT5*(QTDITS*EJ*QTIJTS+QTIITS*EI*QTDJTS)/R3
c                        write(6,'("QI=",2F12.8)') QTDITS,QTIITS
c                        write(6,'("EPS=",2F12.8)') EI,EJ
                         QITS  = EI*QTIJTS/R3
                         QJTS  = EJ*QTIITS/R3
C
                        END IF
C
                        IF(TDDFT) THEN
                           SCALE = EPS/(EPS-ONE)
                           FT=FT+SCALE*(2.0D+00*QSEXY(ITS)*QSEXY(JTS)/R3
     *                                + QTIJTS*QSEMN(ITS)/R3
     *                                + QTIITS*QSEMN(JTS)/R3)
                        END IF
                        IF(MP2PCM)THEN
                           SCALE = EPS/(EPS-ONE)
                           FT=FT+SCALE*(QTIJTS*QSEMP(ITS)/R3
     *                                + QTIITS*QSEMP(JTS)/R3)
                        END IF
C
C                           ALL IEF OTHER THAN CPCM:
                     ELSE
                        IF (LFMOPCM1) THEN
                          ! MAYBE THIS IS WRONG
                          FT = QT(ITS)*QSETMP(JTS)/R3
                        ELSE
                          FT = QT(ITS)*QT(JTS)/R3
                        END IF
                     END IF
C
                     FQQX = FQQX - (XI-XJ)*FT
                     FQQY = FQQY - (YI-YJ)*FT
                     FQQZ = FQQZ - (ZI-ZJ)*FT
C
                     FEQ(1,IAT,ITS)= FEQ(1,IAT,ITS) - (XI-XJ)*QITS
                     FEQ(2,IAT,ITS)= FEQ(2,IAT,ITS) - (YI-YJ)*QITS
                     FEQ(3,IAT,ITS)= FEQ(3,IAT,ITS) - (ZI-ZJ)*QITS
                     FEQ(1,IAT,JTS)= FEQ(1,IAT,JTS) - (XI-XJ)*QJTS
                     FEQ(2,IAT,JTS)= FEQ(2,IAT,JTS) - (YI-YJ)*QJTS
                     FEQ(3,IAT,JTS)= FEQ(3,IAT,JTS) - (ZI-ZJ)*QJTS
C
                  END IF
               END DO
            END IF
         END DO
         END IF
C
C
C        TERM 5. FORCE DUE TO TESSARA AREA CHANGE (FIXPVA)
C
         FAX = ZERO
         FAY = ZERO
         FAZ = ZERO
         IF(IPCDER.EQ.3 .AND. NESFP.EQ.0) THEN
         SQRTPI = SQRT(3.1415926535897932D+00)
         DO ITS = 1, NTS
            DO III = 1, IDDAI(21,ITS)
               IF(IATGLB.EQ.IDDAI(III,ITS)) THEN
                  IF(IEF.EQ.10) THEN
                     QTIITS=QSE(ITS)+QSE(ITS)-QSED(ITS)
     *                     +QSN(ITS)+QSN(ITS)-QSND(ITS)
     *                                         +Q_FS(ITS)+Q_IND(ITS)
                     IF (LFMOPCM1) THEN
                       ! THIS PART SHOULD BE FIXED?
                       QTDITS=QSETMP(ITS)
                       DUM=-1.07D+00*SQRTPI/SQRT(AS(ITS)**3)*
     *                      QTDITS*QTIITS
                     ELSE
                       QTDITS=QSED(ITS)+QSND(ITS)+Q_FS(ITS)+Q_IND(ITS)
                       DUM =-1.07D+00*SQRTPI/SQRT(AS(ITS)**3)*
     *                      0.5D+00*QTDITS*QTIITS
                       DUM2=-1.07D+00*SQRTPI/SQRT(AS(ITS)**3)*QTDITS
                     END IF
C
                     IF(TDDFT) THEN
                        DUM = DUM -1.07D+00*SQRTPI/SQRT(AS(ITS)**3)*
     *                     (QSEXY(ITS)*QSEXY(ITS)+QTIITS*QSEMN(ITS))
                     END IF
                     IF(MP2PCM)THEN
                        DUM = DUM -1.07D+00*SQRTPI/SQRT(AS(ITS)**3)*
     *                     (QTIITS*QSEMP(ITS))
                     END IF
C
                     SCALE=EPSHET(ISPHE(ITS))/(EPSHET(ISPHE(ITS))-ONE)
                     FAX=FAX+DUM*DAI(1,III,ITS)*SCALE
                     FAY=FAY+DUM*DAI(2,III,ITS)*SCALE
                     FAZ=FAZ+DUM*DAI(3,III,ITS)*SCALE
C
                     FX = DUM2*DAI(1,III,ITS)*SCALE 
                     FY = DUM2*DAI(2,III,ITS)*SCALE
                     FZ = DUM2*DAI(3,III,ITS)*SCALE
                     FEQ(1,IAT,ITS)= FEQ(1,IAT,ITS) + FX
                     FEQ(2,IAT,ITS)= FEQ(2,IAT,ITS) + FY
                     FEQ(3,IAT,ITS)= FEQ(3,IAT,ITS) + FZ
C
                  ELSE
                     IF (LFMOPCM1) THEN
                       DUM=-1.07D+00*SQRTPI/SQRT(AS(ITS)**3)*
     *                      QSETMP(ITS)*QT(ITS)
                     ELSE
                       DUM=-1.07D+00*SQRTPI/SQRT(AS(ITS)**3)*
     *                      0.5D+00*QT(ITS)*QT(ITS)
                     END IF
                     FAX=FAX+DUM*DAI(1,III,ITS)
                     FAY=FAY+DUM*DAI(2,III,ITS)
                     FAZ=FAZ+DUM*DAI(3,III,ITS)
                  ENDIF
               END IF
            ENDDO
         ENDDO
         END IF
C
C
C        COMBINE THE FORCES
c        GX = FQQX   +  FAX
c        GY = FQQY   +  FAY
c        GZ = FQQZ   +  FAZ
c        write(6,'("CHK1=",3F12.8)') GX,GY,GZ
C
C
         DSOLQM(1,IAT) = DSOLQM(1,IAT)+ FQQX   +  FAX
         DSOLQM(2,IAT) = DSOLQM(2,IAT)+ FQQY   +  FAY
         DSOLQM(3,IAT) = DSOLQM(3,IAT)+ FQQZ   +  FAZ
 100  CONTINUE
C
c     DO IAT=1,NAT
c       valx=0.0D+00
c       valy=0.0D+00
c       valz=0.0D+00
c       DO ITS=1,NTS
c         valx=valx+0.5D+00*QT(ITS)*FEQ(1,IAT,ITS)
c         valy=valy+0.5D+00*QT(ITS)*FEQ(2,IAT,ITS)
c         valz=valz+0.5D+00*QT(ITS)*FEQ(3,IAT,ITS)
c       END DO
c       write(6,'("CHK2=",3F12.8)') valx,valy,valz
c     END DO
C
      RETURN
      END
C
C
C*MODULE PCMHSS  *DECK CCMBEMFMOH
C>    @brief Second order derivatiev for PCM for FMO
C>
C>    @details  Derivative of  G
C>
C>    @author Hiroya Nakata
C>
      SUBROUTINE CCMBEMFMOH(DEN,DSOLQM,FEQ,QT,VECTMP,
     *           L2,NATM,XCTS,YCTS,ZCTS,AS,QSN,QSE,QSND,QSED,Q_FS,
     *           Q_IND,Q2,XE,YE,ZE,RE,ISPHE,
     *           PCMCDR,IAGLOB,DAI,IDDAI,
     *           QSEXY,QSEMN,FEQWRK,QSEMP,
     *           EH,FD,EHQ,FDQ,nat2,MAP,imode,iwrk,natb,lbglb)
      use mx_limits, only: mxatm
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
c     MXFRG=1050, MXDFG=5, MXDPPT=MXFRG*MXDFG*12)
c     PARAMETER (ONE=1.0D+00,HALF=0.5D+00)
C
C
      COMMON /FMCOM / XX(1)
      COMMON /FMOINF/ NFG,NLAYER,NATFMO,NBDFG,NAOTYP,NBODY,NSEGM
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
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON /MLTFMO/ Q_MUL,IMLTFMO,ISWPFD,ISWNEW1,lfmoicm
      COMMON /MLTPNT/ LDELEPOT,LDELQSE,LQSETMP,LFZVEC,LMAPCMFD,LQBL2FD,
     *                LFRZQ
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /PCMDIM/ MXSP,MXTS,MEMPCM1,MEMPCM2,NTS
      COMMON /PCMPAR/ IPCM,NFT26,NFT27,IKREP,IEF,IP_F,NFMOPCM,IHET
      COMMON /SCYVEC/ IDFYVEC,LIPTYVEC,LIPTYA,LFCMWRK
      COMMON /SCZLAG/ LZVLAG,LZVWRK,LYALAG,LYAWRK,LFEQ1
C
      LOGICAL GOPARR, DSKWRK, MASWRK
      LOGICAL LFMOPCM1
C
      DIMENSION DEN(L2),DSOLQM(3,NATM),
     *          FEQ(3,NATM,*),QT(NTS),VECTMP(*),
     *          XCTS(*),YCTS(*),ZCTS(*),AS(*),QSN(*),QSE(*),
     *          QSND(*),QSED(*),Q_FS(*),
     *          Q_IND(*),Q2(*),XE(*),YE(*),ZE(*),RE(*),ISPHE(*),
     *          PCMCDR(3,*),IAGLOB(*),
     *          QSEXY(*),QSEMN(*),FEQWRK(9,*),QSEMP(*),
     *          EH(9,*),FD(L2,3,natm),
     *          EHQ(9,nat2),FDQ(L2,3,natm,*),MAP(*)
      DIMENSION DAI(3,20,MXTS),IDDAI(21,MXTS)
      DIMENSION IWRK(*)
      DIMENSION LBGLB(*)
C
      DATA ZERO/0.0D+00/
      DATA ONE/1.0D+00/
c     DATA RHF,UHF,ROHF,GVB/8HRHF     ,8HUHF     ,8HROHF    ,8HGVB     /
c     DATA NONE,RNONE/4HNONE,8HNONE    /
c     data excite,spnflp/8HEXCITE  ,8HSPNFLP  /
C---------------------------------------------------------------------
C     COMPUTE THE C-PCM AND IEF-PCM Hessian 
C---------------------------------------------------------------------
C
C     IN AN QM/EFP/PCM RUN, AN EFP(EFPA) AND ITS
C     ASC(ASCA) EXPERIENCE THE FOLLOWING FORCES/TORQUES:
C        EFPA <-- OTHER ASC
C        ASCA <-- OTHER ASC, OTHER EFP, QM ATOMS
C     FOR AN QM ATOM(ATOMB) AND ITS ASC(ASCB):
C        ATOMB <-- OTHER ASC
C        ASCB  <-- OTHER ASC, OTHER QM ATOMS, EFP
C     ***********************************************
C
C
C     ISWPFD=1: FOR FMO/PCM<1> COMPLETELY ANALYTIC GRADIENT
C               NOTE THAT ONLY FMO/CPCM[1] COMPLETELY ANALYTIC GRADIENT
C                         IS CURRENTLY AVAILABLE.
      if(nat.eq.0) write(6,*) "PCMCDR=",PCMCDR(1,1)
C
C
      LFMOPCM1 = ISWPFD.EQ.1
      LFMOPCM1 = .false.
      L1   = NUM
      L3   = L1*L1
      NEGH = 3*NAT + 9*(NAT*NAT+NAT)/2
      CALL VCLR(DSOLQM,1,NEGH)
      CALL VCLR(FD,1,3*natb*L2)
C
C
C
C     -- READ DENSITY MATRIX
C
      IF(NAT.GT.0) THEN
         CALL DAREAD(IDAF,IODA,DEN,L2,16,0)
c        IF(SCFTYP.EQ.ROHF .OR. SCFTYP.EQ.UHF
c    *                     .OR. SCFTYP.EQ.GVB) THEN
c           CALL VALFM(LOADFM)
c           LDENB = LOADFM + 1
c           LAST  = LDENB  + L2
c           NEED  = LAST - LOADFM - 1
c           CALL GETFM(NEED)
c           CALL DENDD1(DEN,XX(LDENB),L2)
c           CALL RETFM(NEED)
c        END IF
      END IF
C
c     -- COMBINE THE SURFACE CHARGES
c     write(6,'("Check charge")')
c     write(6,'("VAL="4F12.7)') 
c    *     (QSE(III),QSN(III),Q_FS(III),Q_IND(III),iii=1,nts)
C
      if(imode.ne.1) then
        DO ITS=1,NTS
           QT(ITS)=QSE(ITS)+QSN(ITS)+Q_FS(ITS)+Q_IND(ITS)
        ENDDO
      end if
      if(imode.eq.1) then
        CALL DCOPY(nts,xx(LDELQSE),1,QT,1)
        DO ITS=1,NTS
           QT(ITS+nts)=QSE(ITS)+QSN(ITS)+Q_FS(ITS)+Q_IND(ITS)
        ENDDO
c       write(6,'("wwwchk=",3F10.5)') (xx(ldelepot+iii-1),iii=1,nts)
c       write(6,'("wwwchk=",3F10.5)') (QT(iii),iii=1,nts)
      end if
C
C     -- ZERO OUT Q2 FOR IEF=10
C
      IF(IEF.EQ.10) THEN
         DO ITS=1,NTS
            Q2(ITS)=ZERO
         END DO
      END IF
C
c     DO ITS=1,nts
c       write(6,*) "ISP=",ITS,ISPHE(its)
c     END DO
      nc1 = natb * 3
      CALL VCLR(XX(LFCMWRK),1,nc1*nc1*2)
C
c     return
      if(imode.eq.2) go to 300
C
c     write(6,*) "WWW Enter PCM Hessian 2nd",nts,natfmo,natb
c     do iat=1,natfmo
c      write(6,'("www=",4I8)') iat,MAP(IAT),MAP(IAT+1)
c     end do
C
      DO 100 IAT=1,NATFMO
       NINITS=MAP(IAT)
       NENDTS=MAP(IAT+1)
       NLEN  =NENDTS-NINITS
       if(nlen.eq.0) go to 100
C      check for translation invariance
       ICENT =0
       DO JAT=1,NAT
         jzat = iaglob(jat)
         if(iat.eq.jzat) icent = jat
       END DO
 
c      write(6,'("CHK DIM=",4I4)') NINITS,NENDTS,NLEN,NTS
C
       CALL VCLR(FEQ,1,3*NAT)
       CALL VCLR(EHQ,1,9*nat2)
       CALL VCLR(FEQWRK,1,9*nat2)
       CALL VCLR(FDQ,1,3*L2*NAT)
C      TVD! One-Electron interaction
       CALL PCMTVDSPD(FEQ,EHQ,DEN,FDQ,XCTS(NINITS),YCTS(NINITS),
     *                ZCTS(NINITS),NLEN,QT(NINITS))
C
C
C      Translation invraiance is used
       CALL TRANINVHSSFMO(icent,NAT,FEQ,FDQ,EHQ,DSOLQM,FD,EH,
     *                    L2,XX(LFCMWRK+nc1*nc1),IAT,LBGLB,iwrk,nc1)
       CALL DAXPY(nc1*nc1,1.0D+00,XX(LFCMWRK+nc1*nc1),1,XX(LFCMWRK),1)
       CALL VCLR(XX(LFCMWRK+nc1*nc1),1,nc1*nc1)
C
 100  continue
C
C
 300  continue
CXC
      CALL VALFM(LOADFM)
      NDASC = 0
      LFPNT = LOADFM
      NTSPAR=(NTS-1)/nproc + 1
      CALL INDQPROD(0,NTSPAR,LFPNT,NDASC,DUM,DUM,DUM,0,0)
      LQTMP   = LFPNT   + NDASC  + 1
      LWRK    = LQTMP   + nts    
      LWRK2   = LWRK    + max(L2,natfmo*9)
      LFCMTMP = LWRK2   + max(L2,natfmo*9)
      LFEQ2   = LFCMTMP + nat      * nat * 9
      LAST    = LFEQ2   + nc1      * nts
      NEEDW   = LAST    - LOADFM   - 1
      CALL GETFM(NEEDW)
C     NUCLEAR REPULSION
      if(imode.ne.2) then
        CALL PCMNUCHSSFMO(ISPHE,XCTS,YCTS,ZCTS,NTS,QT,XX(LWRK),
     *       XX(LFMOZAN),XX(LFMOC),iaglob,XX(LFCMWRK+nc1*nc1),iwrk,nc1)
        CALL DAXPY(nc1*nc1,1.0D+00,XX(LFCMWRK+nc1*nc1),1,XX(LFCMWRK),1)
      end if
C
      CALL VCLR(feq,1,nat*3*nts)
      call dacopy(nts,one,VECTMP,1)
      CALL FELECQ(DEN,VECTMP,FEQ,L2,XCTS,YCTS,ZCTS)
C
c     if(maswrk) then
c     do iat = 1,nat
c       write(6,'("www chk 0=",3F12.8)')  (FEQ(iii,iat,1),iii=1,3)
c     end do
c     end if
C
C     Add translation Invariance
      CALL VCLR(XX(lfeq2),1,natb*3*nts)
      CALL FNUCNT2(FEQ,XCTS,YCTS,ZCTS,NTS,ISPHE,
     *     IAGLOB,XX(LFMOZAN),XX(LFMOC),1.00D+00,XX(LFEQ2),
     *     natb,iwrk,lbglb)
C
      CALL VCLR(feq,1,nat*3*nts)
C
      if(maswrk) then
        do iat = 1,nat
          do ixyz = 1,3
            izat  = lbglb(iat)
            jpnt  = LFEQ2 + nts * (3 * (izat -1) + ixyz - 1)
            CALL DCOPY(nts,XX(JPNT),1,vectmp,1)
            do ii = 1,nts
              feq(ixyz,iat,ii) = vectmp(ii)
            end do
          end do
        end do
      end if

c     if(maswrk) then
c     do iat = 1,nat
c       write(6,'("www chk 1=",3F12.8)')  (FEQ(iii,iat,1),iii=1,3)
c     end do
c     end if


c     IF(GOPARR) THEN
c       CALL DSCAL(NAT*3,ONE/NPROC,DSOLQM,1)
c       CALL DSCAL(NAT*3*NTS,ONE/NPROC,FEQ,1)
c     END IF
C
C
c     if(imode.ne.0) then 
      CALL QCHESS(NATM,0,FEQWRK,QT,VECTMP,XCTS,YCTS,ZCTS,
     *        AS,QSN,QSE,QSND,QSED,Q_FS,Q_IND,Q2,XE,YE,ZE,
     *        RE,ISPHE,DAI,IDDAI,XX(LQSETMP),QSEMP,QSEXY,QSEMN,
     *        DSOLQM,LFMOPCM1,IAGLOB,XX(LFCMTMP),FEQ,XX(LFMOC),
     *        imode)
C
C
      IF(GOPARR) then
        CALL DDI_GSUMF(2461,DSOLQM,3*NAT)
        CALL DDI_GSUMF(2462,XX(LFCMTMP),9*NAT*NAT)
        CALL DDI_GSUMF(2463,FEQ,NAT*3*NTS)
      END IF
C
c     NOCC = NA
c     NVIR = NQMT - NA
c     NROT = NOCC * NVIR
      JINI = LYAWRK      + l1*l1* NAT * 3
C
C     check
c     CALL VCLR(FD,1,3*NAT*L2)
c     do i = 1,nat
c      do j=1,3
c       val = TRACEP(DEN,FD(1,j,i),num) 
c       write(*,'("gradient0 ="2I4,F10.7)') i,j,val
c      end do
c     end do
C
      do iat=1,nat
       IZAT = lbglb(iat)
C
       do ixyz=1,3
C        FEQ --> VECTMP
         DO ITS=1,NTS
          VECTMP(ITS)=FEQ(ixyz,iat,its)
         END DO
C
         JPNT = JINI + nts * ((iat -1)*3 + ixyz -1)
         CALL DCOPY(nts,vectmp,1,XX(JPNT),1)
C
         ITER0 = 1
         CALL INDQPROD(ITER0,NTSPAR,LFPNT,NDASC,VECTMP,
     *                 XX(LQTMP),XX(LQTMP),1,0)
C
C        Fock derivative terms
         CALL PCMPOT(XX(lwrk),XX(lwrk),NTS,XX(LQTMP),XCTS,L2,1,.false.)
C
         if(izat.ne.0) 
     *     CALL DAXPY(L2,1.0D+00,xx(lwrk),1,FD(1,ixyz,izat),1)
c        Complete HSS cont is in fmopcmhssesp
         CALL PCMCMPHSS(XX(LFCMTMP),NTS,XX(LQTMP),FEQ,NAT,iat,ixyz)
       end do
      end do
C
C     check
c     CALL VCLR(FD,1,3*NAT*L2)
c     do i = 1,nat
c      do j=1,3
c       val = TRACEP(DEN,FD(1,j,i),num) 
c       write(*,'("gradient1 ="2I4,F10.7)') i,j,val
c      end do
c     end do
C
      call pcmaddhss(EH,XX(LFCMTMP),NAT)
c     end if
C
      CALL RETFM(NEEDW)
C
C
c     if(ifmostp.eq.4) then
c       CALL DCOPY(NEGH,EH,1,FEQ(1,1,2),1)
c       call pcmaddhss(FEQ(1,1,2),XX(LFCMWRK),NAT)
c       WRITE(IW,*) "wwwcheck hessian"
c       CALL HSSPRT(NAT,FEQ,FEQ(1,1,2))
c       CALL VCLR(XX(LFCMWRK),1,nc1*nc1)
c     end if
C
C     Save gradient and Hessian
      if(imode.eq.0) then
        CALL DAREAD(IDAF,IODA,FEQ,NEGH,67,0)
        CALL DAXPY(NEGH,1.0D+00,FEQ,1,DSOLQM,1)
      end if
      CALL DAWRIT(IDAF,IODA,DSOLQM,NEGH,67,0)
C     Save Fock derivative
      CALL VALFM(LOADFM)
      LVEC  = LOADFM + 1
      LWRK1 = LVEC   + L3
      LWRK2 = LWRK1  + L3
      LWRK3 = LWRK2  + L3
      LAST  = LWRK3  + L3
      NEED  = LAST   - LOADFM - 1
      CALL GETFM(NEED)
C
      CALL DAREAD(IDAF,IODA,XX(LVEC),L3,15,0)
C
      NFT18 = 18
c     NXYZ  = NAT * 3
      if(imode.eq.0) then
       CALL SEQREW(NFT18)
       DO IAT=1,NAT
c       IZAT = iaglob(iat)
        IZAT = lbglb(iat)
        if(izat.ne.0) then
C
C
C
C
         DO IXYZ=1,3
C
c         write(6,*) "www fd check 0 =",izat,ixyz
c         CALL PRTRI(FD(1,IXYZ,IZAT),L1)
C
          CALL TFTRI(XX(LWRK1),FD(1,IXYZ,IZAT),XX(LVEC),XX(LWRK2),
     *               L1,L1,L1)
          CALL SQREAD(NFT18,XX(LWRK3),L2)
          CALL VADD(XX(LWRK3),1,XX(LWRK1),1,FD(1,IXYZ,IZAT),1,L2)
         END DO
        end if
       END DO
C
       CALL SEQREW(NFT18)
C
       DO IAT=1,NAT
        IZAT = LBGLB(iat)
        if(izat.ne.0) then
         DO IXYZ=1,3
C
c         write(6,*) "www fd check =",iat,ixyz
c         CALL PRTRI(FD(1,IXYZ,IZAT),L1)
C
          CALL SQWRIT(NFT18,FD(1,IXYZ,IZAT),L2)
         END DO
        end if
       END DO
C
       DO IATM = 1, natb 
         INFRG = 0
         DO JAT = 1, NAT
            JZAT = LBGLB(jat)
            if(iatm.eq.jzat) infrg = 1
         END DO
C
         if(infrg.eq.0) then
          DO IXYZ = 1,3
           CALL TFTRI(XX(LWRK1),FD(1,IXYZ,IATM),XX(LVEC),XX(LWRK2),
     *                L1,L1,L1)
           CALL SQWRIT(NFT18,XX(LWRK1),L2)
          END DO
         end if
       END DO
      end if
C
      if(imode.eq.1) then
C
       CALL SEQREW(NFT18)
C
       DO IATM = 1, NAT
C
         IZAT = LBGLB(iatm)
C
         if(izat.ne.0) then
          DO IXYZ = 1,3
           CALL TFTRI(XX(LWRK1),FD(1,IXYZ,izat),XX(LVEC),XX(LWRK2),
     *                L1,L1,L1)
           CALL SQWRIT(NFT18,XX(LWRK1),L2)
          END DO
         end if
C
       END DO
      end if
C
C
      CALL RETFM(NEED)
C
      RETURN
      END
C
CC
c     SAVGOP=GOPARR
c     GOPARR=.false.
C
CX      IF(GOPARR) THEN
C  X        CALL DSCAL(NAT*3,ONE/NPROC,DSOLQM,1)
C  X        CALL DSCAL(NAT*3*NTS,ONE/NPROC,FEQ,1)
C  X      END IF
CXC     C^ab
CXc
CXC     - SUM UP -
CX      IF(GOPARR) THEN 
CX        CALL DDI_GSUMF(2461,DSOLQM,3*NAT)
CX        CALL DDI_GSUMF(2462,XX(LFCMTMP),9*NAT*NAT)
CX        CALL DDI_GSUMF(2463,FEQ,NAT*3*NTS)
CX      END IF
C
c     GOPARR=SAVGOP
C     - C^{-1} * V^a
C
C
C*MODULE PCMHSS  *DECK FNUCNT2
C>    @brief Second order derivatiev for PCM for FMO
C>
C>    @details  Derivative of V^a for (q = C^-1 V^a)
C>
C>    @author Hiroya Nakata
C>
C
      SUBROUTINE FNUCNT2(FEQ,XCTS,YCTS,ZCTS,NTS,ISPHE,iaglob,
     *           FMOZAN,FMOC,fda,FEQ2,natb,iwrk,lbglb)
      use mx_limits, only: mxatm
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      DIMENSION FEQ(3,nat,*)
      DIMENSION XCTS(*),YCTS(*),ZCTS(*)
      DIMENSION ISPHE(*),iaglob(*)
c     DIMENSION WRK(3,*)
      DIMENSION FMOZAN(*),FMOC(3,*)
      DIMENSION FEQ2(nts,3,natb)
      DIMENSION IWRK(*),LBGLB(*)
C
      LOGICAL GOPARR, DSKWRK, MASWRK
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK

C
c     DATA ZERO/0.0D+00/
C
c     CALL VCLR(FEQ2,1,3*NATFMO*nts)
C
      DO ITS = 1, NTS
       IAT  = ISPHE(ITS)
       IZAT = iwrk(iat)
C     
c      ICNT = 0
       do jatm = 1,NAT
        jzat = lbglb(jatm)
        if(iat.ne.iaglob(jatm)) then
C
         if(maswrk) then
          if(jzat.ne.0) then
           FEQ2(its,1,jzat) = FEQ2(its,1,jzat) + FEQ(1,jatm,its)
           FEQ2(its,2,jzat) = FEQ2(its,2,jzat) + FEQ(2,jatm,its)
           FEQ2(its,3,jzat) = FEQ2(its,3,jzat) + FEQ(3,jatm,its)
          end if
C
          if(izat.ne.0) then
           FEQ2(its,1,izat) = FEQ2(its,1,izat) - FEQ(1,jatm,its)
           FEQ2(its,2,izat) = FEQ2(its,2,izat) - FEQ(2,jatm,its)
           FEQ2(its,3,izat) = FEQ2(its,3,izat) - FEQ(3,jatm,its)
          end if
         end if
C
        end if
       end do
C
       DO JATFRG = 1, NAT
        JAT = IAGLOB(JATFRG)
        JZAT= lbglb(jatfrg)
        IF(JAT.NE.IAT) THEN
C         - OTHER NUCLEI -
          XI = XCTS(ITS)
          YI = YCTS(ITS)
          ZI = ZCTS(ITS)
          XN = FMOC(1,JAT)
          YN = FMOC(2,JAT)
          ZN = FMOC(3,JAT)
          R2 = (XN-XI)**2+(YN-YI)**2+(ZN-ZI)**2
          R  = SQRT(R2)
          R3 = R*R2
          Z  = FMOZAN(JAT) * fda
c
          DUM = Z/R3
C
          if(maswrk) then
           if(jzat.ne.0) then
            FEQ2(ITS,1,JZAT)=FEQ2(ITS,1,JZAT) - DUM*(XN-XI) 
            FEQ2(ITS,2,JZAT)=FEQ2(ITS,2,JZAT) - DUM*(YN-YI)
            FEQ2(ITS,3,JZAT)=FEQ2(ITS,3,JZAT) - DUM*(ZN-ZI)
           end if
C
           if(izat.ne.0) then
            FEQ2(ITS,1,IZAT)=FEQ2(ITS,1,IZAT) + DUM*(XN-XI) 
            FEQ2(ITS,2,IZAT)=FEQ2(ITS,2,IZAT) + DUM*(YN-YI)
            FEQ2(ITS,3,IZAT)=FEQ2(ITS,3,IZAT) + DUM*(ZN-ZI)
           end if
          END IF
C
        END IF
       END DO
C
C
C
      END DO


      RETURN
      END
C
C
C*MODULE PCMHSS  *DECK PCMNUCHSSFMO
C>    @brief Second order derivatiev for PCM for FMO
C>
C>    @details  Second order Derivative of Nuclear Contribution
C>
C>    @author Hiroya Nakata
C
      SUBROUTINE PCMNUCHSSFMO(ISPHE,XCTS,YCTS,ZCTS,NTS,QT,WRK,
     *           FMOZAN,FMOC,iaglob,FCMWRK,iwrk,nc1)
      use mx_limits, only: mxatm

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      DIMENSION QT(*)
      DIMENSION ISPHE(*),iaglob(*)
      DIMENSION XCTS(*),YCTS(*),ZCTS(*)
      DIMENSION WRK(9,*),FMOZAN(*),FMOC(3,*)
      DIMENSION FCMWRK(nc1,*)
      DIMENSION iwrk(*)
C
C
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /FMOINF/ NFG,NLAYER,NATFMO,NBDFG,NAOTYP,NBODY,NSEGM
C
C
c     NAT2   = (NAT*NAT+NAT)/2
C
C
      DO ITS = 1, NTS
       IAT  = ISPHE(ITS)
       IZAT = iwrk(iat)
       CALL VCLR(wrk,1,9*natfmo)
       DO JATFRG = 1, NAT
        JAT   = iaglob(jatfrg)
        JZAT  = iwrk(jat)
        IF(JAT.NE.IAT) THEN
C         - OTHER NUCLEI -
          XI  = XCTS(ITS)
          YI  = YCTS(ITS)
          ZI  = ZCTS(ITS)
          XN  = FMOC(1,JAT)
          YN  = FMOC(2,JAT)
          ZN  = FMOC(3,JAT)
          R2  = (XN-XI)**2+(YN-YI)**2+(ZN-ZI)**2
          R   = SQRT(R2)
          R3  = R*R2
          Z   = FMOZAN(JAT)
C
          DUM = Z*QT(ITS)/R3
c         hessian
          SCAL1=3.0D+00
          SCAL2=1.0D+00
C
          DUMXX=SCAL1*DUM*(XN-XI)*(XN-XI)/R2 - SCAL2*DUM 
          DUMXY=SCAL1*DUM*(XN-XI)*(YN-YI)/R2 
          DUMXZ=SCAL1*DUM*(XN-XI)*(ZN-ZI)/R2
          DUMYX=SCAL1*DUM*(YN-YI)*(XN-XI)/R2 
          DUMYY=SCAL1*DUM*(YN-YI)*(YN-YI)/R2 - SCAL2*DUM 
          DUMYZ=SCAL1*DUM*(YN-YI)*(ZN-ZI)/R2
          DUMZX=SCAL1*DUM*(ZN-ZI)*(XN-XI)/R2 
          DUMZY=SCAL1*DUM*(ZN-ZI)*(YN-YI)/R2
          DUMZZ=SCAL1*DUM*(ZN-ZI)*(ZN-ZI)/R2 - SCAL2*DUM 
C
          WRK(1,JAT) =  DUMXX
          WRK(2,JAT) =  DUMXY
          WRK(3,JAT) =  DUMXZ
          WRK(4,JAT) =  DUMYX
          WRK(5,JAT) =  DUMYY
          WRK(6,JAT) =  DUMYZ
          WRK(7,JAT) =  DUMZX
          WRK(8,JAT) =  DUMZY
          WRK(9,JAT) =  DUMZZ
C
          jx = (jzat-1)*3 + 1
          jy = (jzat-1)*3 + 2
          jz = (jzat-1)*3 + 3
C
          if(jzat.ne.0) then
           fcmwrk(jx,jx) = fcmwrk(jx,jx) +   DUMXX
           fcmwrk(jx,jy) = fcmwrk(jx,jy) +   DUMXY
           fcmwrk(jx,jz) = fcmwrk(jx,jz) +   DUMXZ
           fcmwrk(jy,jx) = fcmwrk(jy,jx) +   DUMYX
           fcmwrk(jy,jy) = fcmwrk(jy,jy) +   DUMYY
           fcmwrk(jy,jz) = fcmwrk(jy,jz) +   DUMYZ
           fcmwrk(jz,jx) = fcmwrk(jz,jx) +   DUMZX
           fcmwrk(jz,jy) = fcmwrk(jz,jy) +   DUMZY
           fcmwrk(jz,jz) = fcmwrk(jz,jz) +   DUMZZ
          end if
C         
C          
          WRK(1,IAT) =  WRK(1,IAT) + DUMXX
          WRK(2,IAT) =  WRK(2,IAT) + DUMXY
          WRK(3,IAT) =  WRK(3,IAT) + DUMXZ
          WRK(4,IAT) =  WRK(4,IAT) + DUMYX
          WRK(5,IAT) =  WRK(5,IAT) + DUMYY
          WRK(6,IAT) =  WRK(6,IAT) + DUMYZ
          WRK(7,IAT) =  WRK(7,IAT) + DUMZX
          WRK(8,IAT) =  WRK(8,IAT) + DUMZY
          WRK(9,IAT) =  WRK(9,IAT) + DUMZZ
C          
        END IF
       END DO
C
       ix = (izat-1)*3 + 1
       iy = (izat-1)*3 + 2
       iz = (izat-1)*3 + 3
C      Get translation invariance
c      if(.false.) then
       if(izat.ne.0) then
        DO 720 JATM=1,NATFMO
          IF (JATM.EQ.IAT) GOTO 720
          jzat = iwrk(jatm)
C
          jx = (jzat-1)*3 + 1
          jy = (jzat-1)*3 + 2
          jz = (jzat-1)*3 + 3
C
          if(jzat.ne.0) then
           fcmwrk(ix,jx) = fcmwrk(ix,jx) -  WRK(1,JATM)
           fcmwrk(ix,jy) = fcmwrk(ix,jy) -  WRK(2,JATM)
           fcmwrk(ix,jz) = fcmwrk(ix,jz) -  WRK(3,JATM)
           fcmwrk(iy,jx) = fcmwrk(iy,jx) -  WRK(4,JATM)
           fcmwrk(iy,jy) = fcmwrk(iy,jy) -  WRK(5,JATM)
           fcmwrk(iy,jz) = fcmwrk(iy,jz) -  WRK(6,JATM)
           fcmwrk(iz,jx) = fcmwrk(iz,jx) -  WRK(7,JATM)
           fcmwrk(iz,jy) = fcmwrk(iz,jy) -  WRK(8,JATM)
           fcmwrk(iz,jz) = fcmwrk(iz,jz) -  WRK(9,JATM)
C
           fcmwrk(jx,ix) = fcmwrk(jx,ix) -  WRK(1,JATM)
           fcmwrk(jy,ix) = fcmwrk(jy,ix) -  WRK(2,JATM)
           fcmwrk(jz,ix) = fcmwrk(jz,ix) -  WRK(3,JATM)
           fcmwrk(jx,iy) = fcmwrk(jx,iy) -  WRK(4,JATM)
           fcmwrk(jy,iy) = fcmwrk(jy,iy) -  WRK(5,JATM)
           fcmwrk(jz,iy) = fcmwrk(jz,iy) -  WRK(6,JATM)
           fcmwrk(jx,iz) = fcmwrk(jx,iz) -  WRK(7,JATM)
           fcmwrk(jy,iz) = fcmwrk(jy,iz) -  WRK(8,JATM)
           fcmwrk(jz,iz) = fcmwrk(jz,iz) -  WRK(9,JATM)
          end if
C
 720    CONTINUE
C
        fcmwrk(ix,ix) = fcmwrk(ix,ix) +  WRK(1,IAT)
        fcmwrk(ix,iy) = fcmwrk(ix,iy) +  WRK(2,IAT)
        fcmwrk(ix,iz) = fcmwrk(ix,iz) +  WRK(3,IAT)
        fcmwrk(iy,ix) = fcmwrk(iy,ix) +  WRK(4,IAT)
        fcmwrk(iy,iy) = fcmwrk(iy,iy) +  WRK(5,IAT)
        fcmwrk(iy,iz) = fcmwrk(iy,iz) +  WRK(6,IAT)
        fcmwrk(iz,ix) = fcmwrk(iz,ix) +  WRK(7,IAT)
        fcmwrk(iz,iy) = fcmwrk(iz,iy) +  WRK(8,IAT)
        fcmwrk(iz,iz) = fcmwrk(iz,iz) +  WRK(9,IAT)
       end if
C
      END DO
      RETURN
      END
C
C
C*MODULE PCMHSS  *DECK FMOPCMHSSESP
C>    @brief Second order derivatiev for PCM for FMO
C>
C>    @details  Second order Derivative of Coupling V^a*q^b
C>
C>    @author Hiroya Nakata
C
      SUBROUTINE FMOPCMHSSESP(L1,L2,DA,YAWRK,
     *              FMOFCM,NC1,IDAM,IDAD,imode,lfeq1,nextpcm)
      use mx_limits, only: mxatm,mxgtot,mxsh
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (half=0.5D+00)
      PARAMETER (one=1.0D+00)
      PARAMETER (ZERO=0.0D+00)
      PARAMETER (FOUR=4.0D+00)
C
      LOGICAL ISGDDI,PAROUT,INITGDDI,wasgddi,MLGDDI
      LOGICAL GOPARR, DSKWRK, MASWRK
C
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON /FMCOM / X(1)
      common /fmodim/ maxbnd,maxknd,maxcbs,maxcao,maxvec,maxl1,maxnat,
     *                maxabd,maxbas,maxbbd,maxlmo,maxslo,maxabd2,maxrij
      COMMON /FMOINF/ NFG,NLAYER,NATFMO,NBDFG,NAOTYP,NBODY,NSEGM
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
      COMMON /GDDI/   ISCOPE,NGROUPS,MYGROUP,MEGLOB,NPGLOB,NNGLOB,JBTYP,
     *                ISGDDI,PAROUT,INITGDDI,wasgddi,MLGDDI,NSUBGR,
     *                MeUniv,NPUniv,numdlb,myworld,nworlds,mogddi
      COMMON /PCMDIM/ MXSP,MXTS,MEMPCM1,MEMPCM2,NTS
      COMMON /PCMPAR/ IPCM,NFT26,NFT27,IKREP,IEF,IP_F,NFMOPCM,IHET
      COMMON /PCMPNT/ LXYZRE,LSSFE,LLIST,LALPCM,LRINPCM,LINAPCM,LINFPCM,
     *                LINIPCM,LMEPCM,LAXYZCT,LXYZCT2,LNVERT,LQSN,LQSE,
     *                LQFS,LQIND,LISPHE,LVAD,LQOR,LVECMUL,LAIPRJ,LFIPRJ,
     *                LPEL,LPCMCDR,LDAI,LIDDAI,LQSND,LQSED,LEPSPCM
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
C
      DIMENSION DA(*),YAWRK(*),FMOFCM(nc1,*)
C
      ncursh=0
C
      FDA =IDAM*half+IDAD
      if(nfmopcm.eq.2) FDA =IDAM+IDAD
      if(fda.eq.zero.and.imode.eq.0.and.ifmostp.eq.4) return
C
c     nc3= nc1 * nc1
      L3 = L1  * L1
      if(imode.eq.0) then
        CALL DAREAD(IDAF,IODA,DA(l2+1),L3,15,0)
        call dmtx2(da,DA(l2+1),na,l1,l1,na) 
      end if
C
c     DODDCOR = IAND(IXESP,1024).NE.0
C
c     ILAY    = ICURLAY
c     NATA    = NAT
C
      CALL VALFM(LOADFM)
      NDASC = 0
      if(imode.eq.0) then
        LQTMP1  = LOADFM  + 1
        LFEQ    = LQTMP1  + nts
        LWRK1   = LFEQ    + nts * 3 * max(maxnat,nat)
        LIWRK   = LWRK1   + 3       * max(maxnat,nat)
        LBGLB   = LIWRK   + natfmo
        LAST    = LBGLB   + max(maxnat,nat)
        NEEDW   = LAST    - LOADFM    - 1
C
        CALL GETFM(NEEDW)
C
      end if
      if(imode.eq.1) then 
C
c       write(6,*) "www check job ?",nts,nat,maxnat,nc1,l2,natfmo
C
        LFEQ    = LOADFM  + 1 
        LFEQ2   = LFEQ    + nts * 3 * max(maxnat,nat)
c       LQTMP1  = LFEQ2   + max(nts*nc1,l2*3*nat)
        LQTMP1  = LFEQ2   + nts*nc1
c       LQTMP1  = LWRK1   + 3       * max(maxnat,nat)
        LIWRK   = LQTMP1  + nts
        LBGLB   = LIWRK   + natfmo
c       LAST    = LBGLB   + max(nat,maxnat)
        LAST    = LBGLB   + maxnat
        NEEDW   = LAST    - LOADFM    - 1
C
        CALL GETFM(NEEDW)
C
        CALL VALFM(LOADFM)
        LFPNT   = LOADFM
        NTSPAR  = (NTS-1)/npglob+1
        CALL INDQPROD(0,NTSPAR,LFPNT,NDASC,DUM,DUM,DUM,0,0)
        LQTMP   = LFPNT   + NDASC  + 1
        LAST    = LQTMP   + nts    * 4  
        NEEDI   = LAST    - LOADFM    - 1
C
        CALL GETFM(NEEDI)
C
      end if
C
      LXPNT = LAXYZCT 
      LYPNT = LAXYZCT  + MXTS
      LZPNT = LAXYZCT  + MXTS *2
      LAST  = LAXYZCT  + MXTS *3
C
      call mapfmohess(natfmo,x(LIAGLOB),x(LINDAT),x(LLAYFRG),
     *   x(liwrk),natb,4)
C
c     write(6,*) "www check=",nat,natb,natc
C
      if(imode.eq.1) go to 300
C
      call mapfmohess(nat,x(LIAGLOB),x(LINDAT),x(LLAYFRG),
     *   x(lbglb),natc,0)
C
C     V^a * C^{-1} * V^b ||  V(ij) * C^{-1} V^a (Derivative of W)
      call dacopy(nts,fda,X(LQTMP1),1)
      CALL FELECQ(DA,X(LQTMP1),X(LFEQ),L2,X(LXPNT),X(LYPNT),X(LZPNT))
      CALL FNUCNT2(X(LFEQ),X(LXPNT),X(LYPNT),X(LZPNT),NTS,X(LISPHE),
     *     X(LIAGLOB),X(LFMOZAN),X(LFMOC),fda,X(LFEQ1),
     *     natb,x(liwrk),x(lbglb))
C
c     if(ifmstp.eq.2) then
c       lpnt = LFEQ1 +  nc1*nts + nts
c       call dacopy(nts,one,X(LQTMP1),1)
c       CALL FELECQ(DA,X(LQTMP1),X(lpnt),L2,X(LXPNT),X(LYPNT),X(LZPNT))
c       CALL FNUCNT2(X(lpnt),X(LXPNT),X(LYPNT),X(LZPNT),NTS,X(LISPHE),
c    *       X(LIAGLOB),X(LFMOZAN),X(LFMOC),one,X(LFEQ1),
c    *       natb,x(liwrk),x(lbglb))
c     end if
C
      nocc = na
      nvir = nqmt - na
      nrot = nocc * nvir
C
      if(nextpcm.le.0) then
C       it shoud be called only once
        jini = 1 + l1*l1* NAT * 3
        do iat = 1,nat
          izat = IXFTCH(X(lbglb),iat) 
          do ixyz = 1,3
C
            JPNT = JINI  + nts * ((iat  -1)*3 + ixyz -1)
            LPNT = LFEQ1 + nts * ((izat -1)*3 + ixyz -1)
C
            if(maswrk) CALL DAXPY(nts,fda,yawrk(jpnt),1,x(LPNT),1)
c           if(maswrk.and.ifmostp.eq.2) then
c             if(nfmopcm.eq.1) then
c               LPNT = LFEQ1+nts*((izat-1)*3+ixyz-1) + nts*nc1+nts
c               CALL DAXPY(nts,one,yawrk(jpnt),1,x(LPNT),1)
c             end if
c           end if
          end do
        end do
      end if
C
C
      CALL RETFM(NEEDW)

      LQTMP1  = LOADFM  + 1
      LDS     = LQTMP1  + nts
      LWRK    = LDS     + L2 * nat * 3
      LWRK2   = LWRK    + L2 + l3  * 2
      LAST    = LWRK2   + l3
      NEEDW   = LAST    - LOADFM    - 1
C
C
      CALL GETFM(NEEDW)
C
C     (4 Umi V - 2 Smi V)
c     call vclr(x(lds),1,l2*nat*3)
c     if(nfmopcm.ne.2) then
      if(nextpcm.le.0) then
       CALL EXTSDER(X(lwrk),x(lds),L1,L2,.TRUE.,.TRUE.,0)
       IF(GOPARR) CALL DDI_GSUMF(1501,x(lds),L2*3*NAT)
       do iat = 1,nat
        do iixyz = 1,3
         ixyz = (iat-1)*3 + iixyz
C        D*S^a*D
         lpnt =  lds + l2*(ixyz-1)
         CALL CPYTSQ(DA,X(lwrk2),L1,1)
         CALL TFTRI(x(lwrk),x(lpnt),X(lwrk2),x(lwrk+l2),l1,l1,l1)
         call dscal(l2,-half,x(lwrk),1)
C        C*U^a*C
         ii = 1  + NROT * (ixyz - 1)
         lcc= l2 + 1
         jj = 1  + L2 + nocc * l1 
C        C(1,nocc)*Y(nocc,nvir)*C(nvir,1)
         CALL DGEMM('N','N',l1,nvir,nocc,one,da(lcc),l1,
     *              YAWRK(ii),nocc,ZERO,x(lwrk+l2),L1)
         CALL DGEMM('N','T',l1,l1,nvir,one,x(lwrk+l2),l1,
     *              DA(jj),L1,ZERO,x(lwrk+l2+l3),l1)
         CALL SQ2TRI(L1,L1,x(lwrk+l3+l2),x(lwrk+l2),one)
         CALL DAXPY(L2,four,x(lwrk+l2),1,x(lwrk),1)
C        GET POTENTIAL
         CALL VCLR(x(LQTMP1),1,nts)
         call EPOTEN(-fda,X(LXPNT),X(LYPNT),X(LZPNT),x(LQTMP1),
     *               X(LWRK),NTS,L2)
C
         IF(GOPARR) CALL DDI_GSUMF(2461,x(LQTMP1),nts)
C
         izat = IXFTCH(X(LIAGLOB),iat) 
         lpnt = LFEQ1 + (iixyz - 1 + (izat-1) * 3)*nts 
         if(maswrk) CALL DAXPY(nts,one,x(lqtmp1),1,x(lpnt),1)
        end do
       end do
C
      end if
C
C     Additional term for  ( 1/2 V + DVIJ ) * C^{-1} * (d C/ da) * C^{-1} * V^{b,K}
c     CALL VCLR(x(LQTMP1),1,nts)
c     call pcmnup(1,x(lqtmp1))
c     call dscal(nts,fda,x(lqtmp1),1)
c     call EPOTEN(-fda,x(LXPNT),X(LYPNT),X(LZPNT),x(LQTMP1),da,NTS,L2)
c     IF(GOPARR) CALL DDI_GSUMF(2461,x(LQTMP1),nts)
c     lpnt = LFEQ1 +  nc1*nts
c     if(maswrk) CALL DAXPY(nts,one,x(lqtmp1),1,x(lpnt),1)
C
      CALL RETFM(NEEDW)
C
      if(imode.eq.0) then
       return
      end if
C
 300  CONTINUE
C
      IF(GOPARR) then 
        CALL DDI_GSUMF(2461,x(LFEQ1),nts*natb*3)
      end if
C
      if(nfmopcm.eq.2) CALL DCOPY(nts*nc1,X(LFEQ1),1,X(LFEQ2),1)
C
      do iat=1,natb
       do ixyz=1,3
c       C^-1 FEQ --> LQTMP
        IFEQ    = LFEQ1 + (ixyz - 1 + (iat-1) * 3)*nts 
        ITER0   = 1
c       write(6,*) "www get q^a"
        CALL INDQPROD(ITER0,NTSPAR,LFPNT,NDASC,X(IFEQ),
     *                X(LQTMP),X(LQTMP),1,0)
C
c       LQTMP -->  FEQ
        CALL DCOPY(nts,X(LQTMP),1,X(IFEQ),1)
C
       end do
      end do
C
c     IFEQ  = LFEQ1 +  natb*3*nts
c     ITER0 = 1
c     CALL INDQPROD(ITER0,NTSPAR,LFPNT,NDASC,X(IFEQ),
c    *              X(LQTMP),X(LQTMP),1,0)
c     CALL DCOPY(nts,X(LQTMP),1,X(IFEQ),1)
c
C
cyyy  CALL QCDFCKFMO(NAT,NAT,NTS,X(LFEQ1),X(LQTMP),X(LXPNT),X(LYPNT),
c    *           X(LZPNT),X(LAS),X(LQSN),X(LQSE),X(LQSND),X(LQSED),
c    *           X(LQFS),X(LQIND),X(LISPHE),
c    *           X(LDAI),X(LIDDAI),X(LIAGLOB),X(LEPSPCM),X(LFMOZAN),
c    *           X(LIANFRG),fmofcm,X(IFEQ),ntspar,lfpnt,ndasc,natb,
c    *           x(liwrk),x(lbglb))
C
      do iat=1,natb
       do ixyz=1,3
c       C^-1 FEQ --> LQTMP
        ii      = ixyz     + (iat-1) * 3
        IFEQ1   = LFEQ1    + (ixyz - 1 + (iat-1) * 3)*nts 
C
        do jat=1,natb
         do jxyz=1,3
           jj    = jxyz     + (jat-1) * 3
c       C^-1 FEQ --> LQTMP
           if(nfmopcm.eq.2) IFEQ2 = LFEQ2 +(jxyz-1 + (jat-1)*3)*nts 
           if(nfmopcm.eq.1) 
     *       IFEQ2 = LFEQ1 +(jxyz-1 + (jat-1)*3)*nts + nc1*nts+nts
C
           val0  = ddot(nts,x(IFEQ1),1,x(IFEQ2),1)
C
c          write(6,*) "www pcmcmphss =",ii,jj,val0
           if(maswrk) then
             fmofcm(ii,jj) = fmofcm(ii,jj) + val0 
           end if
C
         end do
        end do
       end do
      end do
C
C     
      CALL RETFM(NEEDI)
C
      if(nfmopcm.eq.2) go to 310
C
CCC      IF (ISGDDI) CALL GDDI_SCOPE(DDI_GROUP)
CCCC
CCC      KFG0    = 1
CCC      KFG1    = NFG
CCCC
CCC      IF (ISGDDI) CALL GDDICOUNT(-1,LGROUP,MYJOB)
CCCC
CCC      DO 100 KFG=KFG0,KFG1
CCC        KFGX = KFG
CCCC
CCC        KLAY  = IXFTCH(X(LLAYFRG),KFGX)
CCCC
CCC        IF(KLAY.LT.NLAYER) GOTO 100
CCCC
CCC        IF (ISGDDI) THEN
CCC          CALL GDDICOUNT(0,LGROUP,MYJOB)
CCC          IF (.NOT.MYJOB) GOTO 100
CCC        END IF
CCCC
CCC        call CLOSDA('DELETE')
CCC        CALL OPENDA(0)
CCC        CALL MAKEMOL(KFGX,0,0,KLAY,0,0,0,0,0,0,0,.true.)
CCCC
CCC        call mapfmohess(natfmo,x(LIAGLOB),x(LINDAT),x(LLAYFRG),
CCC     *     x(liwrk),natb,4)
CCC        call mapfmohess(nat,x(LIAGLOB),x(LINDAT),x(LLAYFRG),
CCC     *     x(lbglb),natc,0)
CCCC
CCCC
CCC        L1K  = NUM
CCC        L2K  = (L1K*L1K+L1K)/2
CCC        L3K  = L1K*L1K
CCCc       natk = nat
CCCC
CCC        IREC0=IDMREC0(KLAY)
CCC        IF(DODDCOR) THEN
CCC          IREC0RHF=1
CCC          IF(IREC0.EQ.1) IREC0RHF=NFG*2+1
CCC          IREC0=IREC0RHF
CCC        ENDIF
CCCC
CCC        IDMREC0K= KFG + IREC0
CCC        NQK     = IAND(NQMTFG(KFG),65535)
CCC        NAK     = ISHFT(IXFTCH(X(LNUMFRG),KFG),-16)
CCCc       MULK    = IXFTCH(X(LMULFG),KFG)
CCCc       NBK     = NAK-MULK+1
CCC        NOCK    = NAK
CCC        NVRK    = NQK - NOCK
CCC        NROTK   = NOCK* NVRK
CCCC
CCC        lenrec = L3k+l1K
CCCC       
CCC        CALL RAREADS(IDAFMO,X(LIODFMO),wrk1(1+L2K),lenrec,
CCC     *               IDMREC0K,0)
CCCC       Density Matrix
CCC        CALL DMTX2(WRK1,wrk1(1+L2K),NAK,L1K,L1K,NAK)
CCCC
CCC        CALL dawrit(IDAF,IODA,wrk1(l2k+1),l3k,15,0)
CCC        CALL dawrit(IDAF,IODA,wrk1,l2k,16,0)
CCCCC
CCC        if(nfmopcm.ne.2) then
CCC         call dacopy(nts,one,X(LQTMP1),1)
CCC         CALL FELECQ(WRK1,X(LQTMP1),X(LFEQ),L2,
CCC     *               X(LXPNT),X(LYPNT),X(LZPNT))
CCC
CCCC
CCCc        write(6,*) "Enter Fnucnt2? KFG",KFG
CCCC
CCC         CALL VCLR(X(LFEQ2),1,3*natb*nts)
CCC         CALL FNUCNT2(X(LFEQ),X(LXPNT),X(LYPNT),X(LZPNT),NTS,X(LISPHE),
CCC     *        X(LIAGLOB),X(LFMOZAN),X(LFMOC),one,X(LFEQ2),
CCC     *        natb,x(liwrk),x(lbglb))
CCCCCC
CCC         CALL VALFM(LOADFM)
CCC         LWRK2   = LOADFM + 1
CCC         LWRK3   = LWRK2  + L2K + l3k * 2
CCC         LSDER   = LWRK3  + l3k
CCC         LAST    = LSDER  + L2K * nat * 3
CCC         NEEDK   = LAST   - LOADFM -1
CCC         CALL GETFM(NEEDK)
CCCC
CCCC        - D*S^a*D + 4 Umi Vmi
CCCc        write(6,*) "NCURSH=",NCURSH
CCCc        call vclr(x(lsder),1,l2k*nat*3)
CCC         CALL EXTSDER(X(lwrk2),x(LSDER),L1K,L2K,.TRUE.,.TRUE.,0)
CCC         IF(GOPARR) CALL DDI_GSUMF(1501,x(LSDER),L2K*3*NAT)
CCC         do iat = 1,nat
CCC          do iixyz = 1,3
CCC           ixyz = (iat-1)*3 + iixyz
CCCC          D*S^a*D
CCC           lpnt =  lsder + l2k*(ixyz-1)
CCC           CALL CPYTSQ(wrk1,X(lwrk2),L1K,1)
CCC           CALL TFTRI(x(lwrk3),x(lpnt),X(lwrk2),
CCC     *                x(lwrk2+l3k),l1k,l1k,l1k)
CCC           call dscal(l2k,-half,x(lwrk3),1)
CCCC          C*U^a*C
CCC           ii = NROTK * (ixyz - 1) + iptya(kfg)
CCC           lcc= l2k + 1
CCC           jj = 1   + L2k + nock * l1k
CCCC          C(1,nocc)*Y(nocc,nvir)*C(nvir,1)
CCC           CALL DGEMM('N','N',l1k,nvrk,nock,one,wrk1(lcc),l1k,
CCC     *                YALAG(ii),nock,ZERO,x(lwrk2+l2k),l1k)
CCC           CALL DGEMM('N','T',l1k,l1k,nvrk,one,x(lwrk2+l2k),l1k,
CCC     *                wrk1(jj),l1k,ZERO,x(lwrk2+l2k+l3k),l1k)
CCC           CALL SQ2TRI(l1k,l1k,x(lwrk2+l3k+l2k),x(lwrk2+l2k),one)
CCCC
CCC           call vclr(x(lwrk3),1,l2k)
CCC           CALL DAXPY(L2K,four,x(lwrk2+l2k),1,x(lwrk3),1)
CCCC          GET POTENTIAL
CCC           CALL VCLR(x(LQTMP1),1,nts)
CCC           call EPOTEN(-one,X(LXPNT),X(LYPNT),X(LZPNT),x(LQTMP1),
CCC     *                 X(LWRK3),NTS,L2K)
CCCC
CCC           IF(GOPARR) CALL DDI_GSUMF(2461,x(LQTMP1),nts)
CCCC
CCC           izat = IXFTCH(X(LIAGLOB),iat) 
CCC           lpnt = LFEQ2 + (iixyz - 1 + (izat-1) * 3)*nts 
CCC           CALL DAXPY(nts,one,x(lqtmp1),1,x(lpnt),1)
CCC          end do
CCC         end do
CCCC
CCC         CALL RETFM(NEEDK)
CCCC
CCC         DO KATM=1,natb
CCC          DO KXYZ=1,3
CCCC          Copy tesserae
CCC           kzxyz   = (katm-1)*3 + kxyz
CCC           LKFEQ   = lfeq2  + (kzxyz-1)*nts
CCCc
CCCC          LOOP for IJ Fragments
CCC           do iat = 1,natb
CCC             do ixyz=1,3
CCC               izxyz = (iat-1)*3 + ixyz
CCC               IJFEQ = LFEQ1 + nts * (izxyz -1)
CCCC
CCC               val0  = ddot(nts,x(LKFEQ),1,x(IJFEQ),1)
CCC               if(maswrk) then
CCC                 fmofcm(izxyz,kzxyz) = fmofcm(izxyz,kzxyz) + val0
CCC                 fmofcm(kzxyz,izxyz) = fmofcm(kzxyz,izxyz) + val0
CCC               end if
CCC             end do
CCC           end do
CCC          END DO
CCC         END DO
CCC        end if
CCCC
CCC        CALL HSSCCM(nfmopcm)
CCCC
CCCc       write(6,*) "www Get memory for KFG",nat
CCCc       write(6,*) "www Get memory for KFG",l3k,l2k
CCCC
CCC        CALL VALFM(LOADFM)
CCC        LEG     = LOADFM + 1
CCC        LEH     = LEG    + nat*3
CCC        LWRK1   = LEH    + 9 * (NAT*NAT+NAT)/2
CCC        LSDER   = LWRK1  + L2K + l3k * 2
CCC        LAST    = LSDER  + L2K * nat * 3
CCCC        
CCC        NEEDK   = LAST   - LOADFM -1
CCC        CALL GETFM(NEEDK)
CCCCC
CCC        if(nfmopcm.ne.2) then
CCC         CALL EXTSDER(X(LWRK1),X(LSDER),L1K,L2K,.TRUE.,.TRUE.,0)
CCC         IF(GOPARR) CALL DDI_GSUMF(1501,x(LSDER),L2K*3*NAT)
CCC         do iat=1,nat
CCC          do ixyz=1,3
CCC           lpnt = LSDER + l2k * (iat*3+ixyz-1)
CCC           CALL TFTRI(X(LWRK1),X(lpnt),wrk1(l2k+1),X(LWRK1+l2k),
CCC     *          L1K,L1K,L1K)
CCC           CALL DCOPY(l2k,x(lwrk1),1,x(lpnt),1)
CCC          end do
CCC         end do
CCC        end if
CCCcCxxx
CCCc       write(6,*) "wwwcheck1=",LEG,lbglb,maxnat
CCCc       do iat=1,nat
CCCc        write(6,*) "wwwcheck1=",
CCCc    *     IXFTCH(X(LIAGLOB),IAT),IXFTCH(X(lbglb),IAT)
CCCc       end do
CCCCC
CCC        NEGH = 3*NAT + 9*(NAT*NAT+NAT)/2
CCC        CALL DAREAD(IDAF,IODA,x(LEG),NEGH,67,0)
CCC        call addhss(nat,natb,x(LEH),x(lfcmwrk),1,x(lbglb))
CCC        if(maswrk) CALL DAXPY(nc3,ONE,x(lfcmwrk),1,fmofcm,1)
CCCC
CCCC
CCCCC
CCC        NFT18 = 18
CCCc       NXYZ  = NAT * 3
CCCcCxxx
CCCc       do iat=1,nat
CCCc        write(6,*) "wwwcheck2=",
CCCc    *     IXFTCH(X(LIAGLOB),IAT),IXFTCH(X(lbglb),IAT)
CCCc       end do
CCC        if(nfmopcm.ne.2) then
CCCC
CCC         CALL SEQREW(NFT18)
CCC         DO IAT = 1,nat
CCCc          IZAT  = IAT
CCC           IZAT  = IXFTCH(X(lbglb),IAT)
CCCc          write(6,*) "www PCM fddchek",izat1,izat
CCC           DO IXYZ=1,3 
CCC            izxyz = (izat-1)*3 + ixyz
CCC            CALL SQREAD(NFT18,X(LWRK1+l3k),l2k)
CCCC           wwwchk
CCCc           write(6,*) "www PCM fddchek"
CCCc           CALL PRTRI(X(lwrk1+l3k),l1k)
CCCC
CCC            CALL CPYTSQ(X(lwrk1+l3k),x(lwrk1),l1k,1)
CCC            do jat=1,nat
CCC             JZAT  = IXFTCH(X(lbglb),JAT)
CCC             do jxyz=1,3
CCCC              atom
CCC               jzxyz = (jzat-1)*3 + jxyz
CCCCC             For U contribusion
CCC               ipnt= iptya(kfg) + NROTK*((jat-1)*3+jxyz-1)
CCC               val = zero
CCC               do ivir=1,nvrk
CCC                do iocc=1,nock
CCC                 jpnt = ipnt  + nock*(ivir-1) + iocc -1
CCC                 lpnt = lwrk1 + l1k*(iocc-1) + nock + ivir - 1
CCC                 val  = val   + X(lpnt)*yalag(jpnt)*four
CCC                end do
CCC               end do
CCCC              For S contribusion
CCC               ipnt = LSDER + l2k * ((jat-1)*3+jxyz-1)
CCC               CALL CPYTSQ(X(ipnt),x(lwrk1+l3k),l1k,1)
CCC               do iocc=1,nock
CCC                do jocc=1,nock
CCC                  jpnt = lwrk1 + l3k + (iocc-1)*l1k + jocc - 1
CCC                  lpnt = lwrk1 +       (iocc-1)*l1k + jocc - 1
CCC                  val  = val   - two * x(jpnt) * x(lpnt)
CCC                end do
CCC               end do
CCC               if(maswrk) then
CCC                 fmofcm(izxyz,jzxyz) = fmofcm(izxyz,jzxyz) + val
CCC                 fmofcm(jzxyz,izxyz) = fmofcm(jzxyz,izxyz) + val
CCC               end if
CCC             end do
CCC            end do
CCC          END DO
CCC         END DO
CCCC
CCCC
CCC        end if
CCCC
CCCC
CCCC
CCC        CALL RETFM(NEEDK)
CCCC
CCCC
CCCC
CCC  100 CONTINUE
CCCC
CCC      IF (ISGDDI) CALL GDDICOUNT( 1,LGROUP,MYJOB)
C
 310  continue
C
      CALL RETFM(NEEDW)
C
C

      RETURN
      END
C
C
C*MODULE PCMHSS  *DECK TRANINVHSSFMO
C>    @brief Second order derivatiev for PCM for FMO
C>
C>    @details  Transration Invaraiance of Hessian
C>
C>    @author Hiroya Nakata
C
C
      SUBROUTINE TRANINVHSSFMO(NCNTR,NAT,FEQ,FDQ,EHQ,DSOLQM,FD,EH,
     *           L2,fcmwrk,IIZAT,IAGLOB,iwrk,nc1)

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      DIMENSION FEQ(3,NAT),FDQ(L2,3,nat),EHQ(9,*)
      DIMENSION DSOLQM(3,NAT),FD(L2,3,nat),EH(9,*)
c     DIMENSION FEQWRK(9,*)
      DIMENSION FCMWRK(nc1,*)
      DIMENSION IAGLOB(*)
      DIMENSION IWRK(*)
      PARAMETER (ONE=1.0D+00)
C
      IZAT = iwrk(IIZAT)
C
C     for gradient
      DO JAT=1,NAT
       DO IXYZ=1,3
        IF(NCNTR.NE.JAT) THEN
         DSOLQM(IXYZ,JAT)  =DSOLQM(IXYZ,JAT)   + FEQ(IXYZ,JAT)
         if(NCNTR.ne.0)
     *      DSOLQM(IXYZ,NCNTR)=DSOLQM(IXYZ,NCNTR) - FEQ(IXYZ,JAT)
        END IF
       END DO
      END DO
C
C     for hessian
      DO JAT=1,NAT
       JZAT = IAGLOB(jat)
       DO IXYZ=1,3
        IF(IZAT.NE.JZAT) THEN
         if(jzat.ne.0)
     *     CALL DAXPY(L2, ONE,FDQ(1,IXYZ,JAT),1,FD(1,IXYZ,JZAT),1)
         if(izat.ne.0)
     *     CALL DAXPY(L2,-ONE,FDQ(1,IXYZ,JAT),1,FD(1,IXYZ,IZAT),1)
        END IF
       END DO
      END DO

C     W^{ab} 
C       1. Normal terms
      DO 700 JATM=1,NAT
        IF(JATM.EQ.NCNTR) GOTO 700
        DO 710 KATM=1,JATM
         IF(KATM.EQ.NCNTR) GOTO 710
         JKAT=(JATM*JATM-JATM)/2 +  KATM
         CALL DAXPY(9,ONE,EHQ(1,JKAT),1,EH(1,JKAT),1)
 710    CONTINUE
 700  CONTINUE

c     if(NCNTR.eq.0) return
C     2. Translation invariance for term 5
      DO 720 JATM=1,NAT
        IF (JATM.EQ.NCNTR) GOTO 720
        JZATM = IAGLOB(jatm)
        if(jzatm.eq.0.or.izat.eq.0) goto 720
        ix    = (izat -1)*3 + 1
        iy    = (izat -1)*3 + 2
        iz    = (izat -1)*3 + 3
        jx    = (jzatm-1)*3 + 1
        jy    = (jzatm-1)*3 + 2
        jz    = (jzatm-1)*3 + 3
        DO 730  KATM=1,NAT
         IF (KATM.EQ.NCNTR) GOTO 730
         IF(JATM.GT.KATM) THEN
           JKAT=(JATM*JATM-JATM)/2 +  KATM
C
           fcmwrk(ix,jx)  = fcmwrk(ix,jx) -EHQ(1,jkat)
           fcmwrk(ix,jy)  = fcmwrk(ix,jy) -EHQ(4,jkat)
           fcmwrk(ix,jz)  = fcmwrk(ix,jz) -EHQ(7,jkat)
           fcmwrk(iy,jx)  = fcmwrk(iy,jx) -EHQ(2,jkat)
           fcmwrk(iy,jy)  = fcmwrk(iy,jy) -EHQ(5,jkat)
           fcmwrk(iy,jz)  = fcmwrk(iy,jz) -EHQ(8,jkat)
           fcmwrk(iz,jx)  = fcmwrk(iz,jx) -EHQ(3,jkat)
           fcmwrk(iz,jy)  = fcmwrk(iz,jy) -EHQ(6,jkat)
           fcmwrk(iz,jz)  = fcmwrk(iz,jz) -EHQ(9,jkat)
 
           fcmwrk(jx,ix)  = fcmwrk(jx,ix) -EHQ(1,jkat)
           fcmwrk(jy,ix)  = fcmwrk(jy,ix) -EHQ(4,jkat)
           fcmwrk(jz,ix)  = fcmwrk(jz,ix) -EHQ(7,jkat)
           fcmwrk(jx,iy)  = fcmwrk(jx,iy) -EHQ(2,jkat)
           fcmwrk(jy,iy)  = fcmwrk(jy,iy) -EHQ(5,jkat)
           fcmwrk(jz,iy)  = fcmwrk(jz,iy) -EHQ(8,jkat)
           fcmwrk(jx,iz)  = fcmwrk(jx,iz) -EHQ(3,jkat)
           fcmwrk(jy,iz)  = fcmwrk(jy,iz) -EHQ(6,jkat)
           fcmwrk(jz,iz)  = fcmwrk(jz,iz) -EHQ(9,jkat)
C
         END IF
         IF(KATM.GE.JATM) THEN
           JKAT=(KATM*KATM-KATM)/2 + JATM
           fcmwrk(ix,jx)  = fcmwrk(ix,jx) -EHQ(1,jkat)
           fcmwrk(ix,jy)  = fcmwrk(ix,jy) -EHQ(2,jkat)
           fcmwrk(ix,jz)  = fcmwrk(ix,jz) -EHQ(3,jkat)
           fcmwrk(iy,jx)  = fcmwrk(iy,jx) -EHQ(4,jkat)
           fcmwrk(iy,jy)  = fcmwrk(iy,jy) -EHQ(5,jkat)
           fcmwrk(iy,jz)  = fcmwrk(iy,jz) -EHQ(6,jkat)
           fcmwrk(iz,jx)  = fcmwrk(iz,jx) -EHQ(7,jkat)
           fcmwrk(iz,jy)  = fcmwrk(iz,jy) -EHQ(8,jkat)
           fcmwrk(iz,jz)  = fcmwrk(iz,jz) -EHQ(9,jkat)
C
           fcmwrk(jx,ix)  = fcmwrk(jx,ix) -EHQ(1,jkat)
           fcmwrk(jy,ix)  = fcmwrk(jy,ix) -EHQ(2,jkat)
           fcmwrk(jz,ix)  = fcmwrk(jz,ix) -EHQ(3,jkat)
           fcmwrk(jx,iy)  = fcmwrk(jx,iy) -EHQ(4,jkat)
           fcmwrk(jy,iy)  = fcmwrk(jy,iy) -EHQ(5,jkat)
           fcmwrk(jz,iy)  = fcmwrk(jz,iy) -EHQ(6,jkat)
           fcmwrk(jx,iz)  = fcmwrk(jx,iz) -EHQ(7,jkat)
           fcmwrk(jy,iz)  = fcmwrk(jy,iz) -EHQ(8,jkat)
           fcmwrk(jz,iz)  = fcmwrk(jz,iz) -EHQ(9,jkat)
         END IF
 730    CONTINUE
 720  CONTINUE
      DO JATM=1,NAT
       jzat = iaglob(jatm)
       if(jzat.ne.0.and.izat.ne.0) then 
        IF (jzat.NE.izat) THEN
         ix  = (izat -1)*3 + 1
         iy  = (izat -1)*3 + 2
         iz  = (izat -1)*3 + 3
         jx  = (jzat -1)*3 + 1
         jy  = (jzat -1)*3 + 2
         jz  = (jzat -1)*3 + 3
C
         fcmwrk(ix,ix)   = fcmwrk(ix,ix) - fcmwrk(ix,jx)
         fcmwrk(ix,iy)   = fcmwrk(ix,iy) - fcmwrk(ix,jy)
         fcmwrk(ix,iz)   = fcmwrk(ix,iz) - fcmwrk(ix,jz)
         fcmwrk(iy,ix)   = fcmwrk(iy,ix) - fcmwrk(iy,jx)
         fcmwrk(iy,iy)   = fcmwrk(iy,iy) - fcmwrk(iy,jy)
         fcmwrk(iy,iz)   = fcmwrk(iy,iz) - fcmwrk(iy,jz)
         fcmwrk(iz,ix)   = fcmwrk(iz,ix) - fcmwrk(iz,jx)
         fcmwrk(iz,iy)   = fcmwrk(iz,iy) - fcmwrk(iz,jy)
         fcmwrk(iz,iz)   = fcmwrk(iz,iz) - fcmwrk(iz,jz)
        END IF 
       end if
      END DO
C


      RETURN
      END
C
C
C*MODULE PCMHSS  *DECK CPWABPCM
C>    @brief Second order derivatiev for PCM for FMO
C>
C>    @details  Contribution of U^a * V^b
C>
C>    @author Hiroya Nakata
C
C
      SUBROUTINE CPWABPCM(NOCC,NVIR,L1,YA,FCM,DA,L2,nc1)
      use mx_limits, only: mxatm
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      LOGICAL GOPARR, DSKWRK, MASWRK
C
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
C
      COMMON /FMOINF/ NFG,NLAYER,NATFMO,NBDFG,NAOTYP,NBODY,NSEGM
      COMMON /FMCOM / X(1)
C
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
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
C
      DIMENSION DA(*)
      DIMENSION YA(NOCC,NVIR,*)
      DIMENSION FCM(nc1,*)
      PARAMETER (TWO=2.0D+00)
      PARAMETER (zero=0.0D+00)
      PARAMETER (FOUR=4.0D+00)

C
c     fda    = ida
      L3     = l1*l1
      CALL DAREAD(IDAF,IODA,DA(L2+1),L3,15,0)
C
      CALL VALFM(LOADFM)
      LWRK1  = LOADFM + 1
      LDSAO  = LWRK1  + L3 * 2 + L2
      LIWRK  = LDSAO  + l2 * nat * 3
      LAST   = LIWRK  + natfmo
      NEED   = LAST   - LOADFM - 1
      CALL GETFM(NEED)
C
      call mapfmohess(nat,x(LIAGLOB),x(LINDAT),x(LLAYFRG),
     *   x(liwrk),natc,0)
C
c     write(6,*) "www check=",nat,natc
C
      CALL EXTSDER(X(lwrk1),x(ldsao),L1,L2,.TRUE.,.TRUE.,0)
      IF(GOPARR) CALL DDI_GSUMF(1501,x(ldsao),L2*3*NAT)
C
      NFT18 = 18
c     NXYZ  = NAT * 3
      CALL SEQREW(NFT18)
C     Discard
      DO IAT = 1,nat
        DO IXYZ=1,3 
         CALL SQREAD(NFT18,X(LWRK1),l2)
        END DO
      END DO
C
      DO IATM = 1,natc
       infrag = 0
       DO JAT = 1,NAT
        jzat = IXFTCH(X(liwrk),jat)
        if(jzat.eq.iatm) infrag = 1
       END DO
       if(infrag.eq.0) then
        do ixyz=1,3
         izxyz = (iatm-1)*3  + ixyz
         CALL SQREAD(NFT18,X(LWRK1+l3),l2)
         CALL CPYTSQ(X(lwrk1+l3),x(lwrk1),l1,1)
         do jat=1,nat
          JZAT  = IXFTCH(X(liwrk),JAT)
          if(jzat.ne.0) then
           do jxyz=1,3
C           atom
             jzxyz = (jzat-1)*3  + jxyz
             jx    = (jat -1)*3  + jxyz
C           For U contribusion
             val = zero
             do ivir=1,nvir
              do iocc=1,nocc
               lpnt = lwrk1 + l1 *(iocc-1)  + nocc + ivir - 1
               val  = val   + X(lpnt)*ya(iocc,ivir,jx)*four
              end do
             end do
C            For S contribusion
             lpnt = ldsao + L2 * (jx-1)
             CALL TFTRI(x(lwrk1+l3*2),X(lpnt),DA(l2+1),X(lwrk1+l3),
     *                  nqmt,l1,l1)
             CALL CPYTSQ(x(lwrk1+l3*2),x(lwrk1+l3),l1,1)
             do iocc=1,nocc
              do jocc=1,nocc
                jpnt = lwrk1 + l3  + (iocc-1)*l1  + jocc - 1
                lpnt = lwrk1 +       (iocc-1)*l1  + jocc - 1
                val  = val   - two * x(jpnt) * x(lpnt)
              end do
             end do
CC            add FCM
             if(maswrk) then
              fcm(izxyz,jzxyz) = fcm(izxyz,jzxyz) + val
              fcm(jzxyz,izxyz) = fcm(jzxyz,izxyz) + val
             end if
           end do
          end if
         end do
        end do
       end if
CC
      END DO

C
C
      CALL RETFM(NEED)
C
C
      RETURN
      END
C
C
C*MODULE PCMHSS  *DECK QCDFCKFMO
C>    @brief    Second order derivative for PCM for FMO
C>
C>    @details  Derivative of V_uv * Q^a 
C>
C>    @author   Hiroya Nakata
C
C
C
      SUBROUTINE QCDFCKFMO(NATA,NATM,NTS,FEQ,QT,XCTS,YCTS,ZCTS,
     *                  AS,QSN,QSE,QSND,QSED,Q_FS,Q_IND,
     *                  ISPHE,DAI,IDDAI,IAGLOB,EPSHET,FMOZAN,IANFRG,
     *                  fcm,QDEL,ntspar,lfpnt,ndasc,natb,iwrk,lbglb)
      use mx_limits, only: mxatm,mxfrg,mxfgpt
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C
      LOGICAL GOPARR, DSKWRK, MASWRK
C
      COMMON /FMOINF/ NFG,NLAYER,NATFMO,NBDFG,NAOTYP,NBODY,NSEGM
C
      DIMENSION FEQ(nts,3,*)
      DIMENSION XCTS(1),YCTS(1),ZCTS(1),AS(1)
      DIMENSION QT(nts,4),QSN(1),QSE(1),QSND(1),QSED(1),
     *          Q_FS(1),Q_IND(1),FMOZAN(1)
      DIMENSION ISPHE(1),DAI(3,20,*),IDDAI(21,*)
      DIMENSION IAGLOB(1)
      DIMENSION EPSHET(*)
      DIMENSION IANFRG(*),fcm(natb*3,*)
      DIMENSION QDEL(*)
      DIMENSION IWRK(*),lbglb(*)

c     COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
c    *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
c     COMMON /PCMEPS/ EPSHET(MXATM),IHET
      COMMON /PCMOPT/ RABI,RASC,REFPOL,THRSLS,DENSLS,WB,WA,ETA2,GD,EVAC,
     *                RHOW,PM,AREATL,AREAKP,BONDRY,OMEGA,RET,FRO,EPSINF,
     *                EPS,DR,RSOLV,VMOL,TCE,STEN,DSTEN,CMF,TABS,IDIRCT,
     *                IPCDER,IDP,ICOMP,IFIELD,ICAV,IDISP,IPRINT,IRETCAV,
     *                ICENT,IFAST,NEVAL,IEFPOL,KEEPSM,IMGABI,IMGASC,NADD
c     COMMON /PCMPAR/ IPCM,NFT26,NFT27,IKREP,IEF,IP_F,NFMOPCM,IHET
      COMMON /PCMPLY/ STOT,VOL,CCX,CCY,CCZ,RDIF,NESF,NESFP,NESFF,I_NESF,
     *                L_AST
c     COMMON /PCMPNT/ LXYZRE,LSSFE,LLIST,LALPCM,LRINPCM,LINAPCM,LINFPCM,
c    *                LINIPCM,LMEPCM,LAXYZCT,LXYZCT2,LNVERT,LQSN,LQSE,
c    *                LQFS,LQIND,LISPHE,LVAD,LQOR,LVECMUL,LAIPRJ,LFIPRJ,
c    *                LPEL,LPCMCDR,LDAI,LIDDAI,LQSND,LQSED,LEPSPCM
c     COMMON /PCMTMP/ IDSP(MXATM+MXFRG*5),DPCMFRG(6,MXFRG),
c    *                EGPCM(3,MXATM)

      DATA ONE /1.0D+00/
c     DATA ZERO, ONE /0.0D+00, 1.0D+00/
c     DATA PT5/0.5D+00/
c     DATA FPI/12.56637062D+00/
C
      if(nata.eq.0) write(6,*) iaglob(1),natm,qsed(1),qsnd(1),lbglb(1)
C
C     CALCULATE -1/2*Q*DC/DA*Q
C     -- COMPUTE THE VECTOR A^(-1)*S*A^(-1)*Q2 STORE IN VECTMP
c     write(6,*) "IEF=",IEF,IPCDER
C
C
C
      DO 200 IAT = 1, natfmo
C 
       IATGLB = IAT
       izat   = iwrk(iat)
       if(izat.eq.0) go to 200
       IF (IANFRG(IAT).NE.1.AND.ABS(FMOZAN(IAT)-ONE).LT.1.0D-08)
     *    GOTO 200
C
C      TERM 2.A FORCE IMPOSED ON ASC(IAT) BY OTHER ASC
       CALL VCLR(qt,1,3*nts)
C 
       IF(NESFP.EQ.0)THEN
        DO ITS = 1, NTS
         IF(ISPHE(ITS).EQ.IATGLB) THEN
          XI = XCTS(ITS)
          YI = YCTS(ITS)
          ZI = ZCTS(ITS)
          DO JTS = 1, NTS
           IF(ISPHE(JTS).NE.IATGLB) THEN
             XJ = XCTS(JTS)
             YJ = YCTS(JTS)
             ZJ = ZCTS(JTS)
             R2 = (XJ-XI)**2+(YJ-YI)**2+(ZJ-ZI)**2
             R  = SQRT(R2)
             R3 = R*R2
             EI = EPSHET(ISPHE(ITS))/(EPSHET(ISPHE(ITS))-ONE)
             EJ = EPSHET(ISPHE(JTS))/(EPSHET(ISPHE(JTS))-ONE)
c            IF(IEF.EQ.10) THEN
             QTIJTS = QDEL(JTS)
             QTIITS = QDEL(ITS)
             QITS   = EI*QTIJTS/R3
             QJTS   = EJ*QTIITS/R3
C
             QT(ITS,1)= QT(ITS,1) - (XI-XJ)*QITS
             QT(ITS,2)= QT(ITS,2) - (YI-YJ)*QITS
             QT(ITS,3)= QT(ITS,3) - (ZI-ZJ)*QITS
             QT(JTS,1)= QT(JTS,1) - (XI-XJ)*QJTS
             QT(JTS,2)= QT(JTS,2) - (YI-YJ)*QJTS
             QT(JTS,3)= QT(JTS,3) - (ZI-ZJ)*QJTS
c
           END IF
          END DO
         END IF
        END DO
       END IF
C
C
C      TERM 5. FORCE DUE TO TESSARA AREA CHANGE (FIXPVA)
       IF(IPCDER.EQ.3 .AND. NESFP.EQ.0) THEN
        SQRTPI = SQRT(3.1415926535897932D+00)
        DO ITS = 1, NTS
          DO III = 1, IDDAI(21,ITS)
            IF(IATGLB.EQ.IDDAI(III,ITS)) THEN
c            IF(IEF.EQ.10) THEN
             QTDITS=QDEL(ITS)
C
             DUM2 =-1.07D+00*SQRTPI/SQRT(AS(ITS)**3)*QTDITS
             SCALE=EPSHET(ISPHE(ITS))/(EPSHET(ISPHE(ITS))-ONE)
C
             FX   = DUM2*DAI(1,III,ITS)*SCALE 
             FY   = DUM2*DAI(2,III,ITS)*SCALE
             FZ   = DUM2*DAI(3,III,ITS)*SCALE
C
             QT(ITS,1)= QT(ITS,1) + FX
             QT(ITS,2)= QT(ITS,2) + FY
             QT(ITS,3)= QT(ITS,3) + FZ
C
            END IF
          ENDDO
        ENDDO
       END IF
C
       do ixyz=1,3
        ITER0   = 1
        CALL INDQPROD(ITER0,NTSPAR,LFPNT,NDASC,QT(1,ixyz),
     *               QT(1,4),QT(1,4),1,0)
        CALL DAXPY(NTS,ONE,QT(1,4),1,FEQ(1,ixyz,izat),1)
       end do
C
  200 CONTINUE
C
C
c     if(nfmopcm.eq.2) return
C
C
C
C     - INITIALIZE PARALLEL -
      DO ITS=1,NTS
         QT(ITS,1)=QSE(ITS)+QSN(ITS)+Q_FS(ITS)+Q_IND(ITS)
      ENDDO


      IPCOUNT = ME - 1
C
      DO 100 IAT = 1, natfmo
C
         izat = iwrk(iat)
         if(izat.eq.0) goto 100
C        GO PARALLEL!
         IF(GOPARR) THEN
           IPCOUNT = IPCOUNT + 1
           IF(MOD(IPCOUNT,NPROC).NE.0) GO TO 100
         END IF
C
         IATGLB = IAT
         ix     = (IZAT-1)*3 + 1
         iy     = (IZAT-1)*3 + 2
         iz     = (IZAT-1)*3 + 3
         IF (IANFRG(IAT).NE.1.AND.ABS(FMOZAN(IAT)-ONE).LT.1.0D-08)
     *     GOTO 100
C
C        TERM 2.A FORCE IMPOSED ON ASC(IAT) BY OTHER ASC
C
         IF(NESFP.EQ.0)THEN
          DO ITS = 1, NTS
           IF(ISPHE(ITS).EQ.IATGLB) THEN
            XI = XCTS(ITS)
            YI = YCTS(ITS)
            ZI = ZCTS(ITS)
            DO JTS = 1, NTS
             IF(ISPHE(JTS).NE.IATGLB) THEN
               XJ = XCTS(JTS)
               YJ = YCTS(JTS)
               ZJ = ZCTS(JTS)
               R2 = (XJ-XI)**2+(YJ-YI)**2+(ZJ-ZI)**2
               R  = SQRT(R2)
               R3 = R*R2
c              IF(IEF.EQ.10) THEN
               EI=EPSHET(ISPHE(ITS))/(EPSHET(ISPHE(ITS))-ONE)
               EJ=EPSHET(ISPHE(JTS))/(EPSHET(ISPHE(JTS))-ONE)
               QTIJTS= QT(JTS,1)
               QTIITS= QT(ITS,1)
c              write(6,'("QI=",2F12.8)') QTDITS,QTIITS
c              write(6,'("EPS=",2F12.8)') EI,EJ
               QITS  = EI*QTIJTS/R3
               QJTS  = EJ*QTIITS/R3
C
C
               FIX   = - (XI-XJ)*QITS
               FIY   = - (YI-YJ)*QITS
               FIZ   = - (ZI-ZJ)*QITS
C                      
               FJX   = - (XI-XJ)*QJTS
               FJY   = - (YI-YJ)*QJTS
               FJZ   = - (ZI-ZJ)*QJTS
C
               do jat =1,natb
                do jxyz = 1,3
                 jj = (jat-1)*3+jxyz
c
                 fcm(ix,jj)=fcm(ix,jj)+FIX*FEQ(its,jxyz,jat) 
                 fcm(iy,jj)=fcm(iy,jj)+FIY*FEQ(its,jxyz,jat)
                 fcm(iz,jj)=fcm(iz,jj)+FIZ*FEQ(its,jxyz,jat)
                 fcm(ix,jj)=fcm(ix,jj)+FJX*FEQ(jts,jxyz,jat) 
                 fcm(iy,jj)=fcm(iy,jj)+FJY*FEQ(jts,jxyz,jat)
                 fcm(iz,jj)=fcm(iz,jj)+FJZ*FEQ(jts,jxyz,jat)
 
                 fcm(jj,ix)=fcm(jj,ix)+FIX*FEQ(its,jxyz,jat) 
                 fcm(jj,iy)=fcm(jj,iy)+FIY*FEQ(its,jxyz,jat)
                 fcm(jj,iz)=fcm(jj,iz)+FIZ*FEQ(its,jxyz,jat)
                 fcm(jj,ix)=fcm(jj,ix)+FJX*FEQ(jts,jxyz,jat) 
                 fcm(jj,iy)=fcm(jj,iy)+FJY*FEQ(jts,jxyz,jat)
                 fcm(jj,iz)=fcm(jj,iz)+FJZ*FEQ(jts,jxyz,jat)
C
                end do
               end do
C
             END IF
            END DO
           END IF
          END DO
         END IF
C
C
C
 100  CONTINUE
C
      RETURN
      END
C
C
C
C*MODULE PCMHSS  *DECK CAVITHESS
C>    @brief    Second order derivative for PCM for FMO
C>
C>    @details  Secand order Derivative of G^cav
C>
C>    @author   Hiroya Nakata
C
      subroutine CAVITHESS(natm,fmoc,fcm,numdif,feqwrk,method)
      use mx_limits, only: mxatm
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C
      LOGICAL GOPARR, DSKWRK, MASWRK,nxt,SAVGOP
C
      COMMON /FMOINF/ NFG,NLAYER,NATFMO,NBDFG,NAOTYP,NBODY,NSEGM
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
C
      COMMON /FMCOM / X(1)
      COMMON /PCMDIM/ MXSP,MXTS,MEMPCM1,MEMPCM2,NTS
      COMMON /PCMPNT/ LXYZRE,LSSFE,LLIST,LALPCM,LRINPCM,LINAPCM,LINFPCM,
     *                LINIPCM,LMEPCM,LAXYZCT,LXYZCT2,LNVERT,LQSN,LQSE,
     *                LQFS,LQIND,LISPHE,LVAD,LQOR,LVECMUL,LAIPRJ,LFIPRJ,
     *                LPEL,LPCMCDR,LDAI,LIDDAI,LQSND,LQSED,LEPSPCM
      COMMON /OUTPUT/ NPRINT,ITOL,ICUT,NORMF,NORMP,NOPK
C
C
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
C
      dimension fmoc(3,*)
      dimension fcm(natm*3,*)
      DIMENSION D(2)
      DIMENSION FEQWRK(3,NATM)
      DIMENSION METHOD(*)
C
      NPRSAV = NPRINT
      NPRINT = 817
C
C
      NXT   = IBTYP.EQ.1
      NEXT  = -1
      KOUNT = -1
C
c     write(6,*) "www nxt =",IBTYP.EQ.1
C
      CALL VCLR(fcm,1,natm*natm*9)
C
C
      NVIB   =  2 
      VIBSIZ =  0.000001D+00
c     VIBSIZ =  0.00001D+00
      D(1)   =  VIBSIZ
      D(2)   = -VIBSIZ
c     LFIXPV = 0
c     IF(IPCDER.EQ.3) LFIXPV = 3*20*MXTS
c     LFIXPV2 = 0
c     IF(IPCDER.EQ.3) LFIXPV2= 21*MXTS/NWDVAR + 1
C     Numerical derivatives of dC/da  :  (d^2 C/da db) 
c     write(6,*) "wwww begin ",natm,mxsp
c     write(6,*) "wwww begin numerical"
      if(numdif.eq.-1234) write(6,*) numdif
      SAVGOP = GOPARR
C
      DO 390 IVIB = 1,NVIB
        DO 380 IAT = 1,NATM
          DO 370 IXYZ = 1,3
C
C
           IF(GOPARR) THEN
             KOUNT = KOUNT + 1
             IF(NXT) THEN
               IF(KOUNT.GT.NEXT) CALL DDI_DLBNEXT(NEXT)
               IF(NEXT.NE.KOUNT) GOTO 370
             ELSE
               IF(MOD(KOUNT,NPROC).NE.ME) GOTO 370
             ENDIF
           ENDIF
C
C
           JXYZ       = IAT * 3 - 3 + IXYZ
           if(nfg.ne.0) then
            fmoc(ixyz,iat) = fmoc(ixyz,iat) + D(IVIB)
           else
            C(IXYZ,IAT)= C(IXYZ,IAT) + D(IVIB)
           end if
C
c          write(6,*) "enter makcvm"
           CALL MAKCVM(1)
c          CALL TESIN(NATM,X(LXYZRE),X(LXYZRE+MXSP),X(LXYZRE+MXSP*2),
c    *                X(LALPCM),X(LRINPCM),X(LINIPCM),
c    *                X(LMEPCM))
c          write(6,*) "enter pedram"
           GOPARR = .false.
           CALL PEDRAM
C          
c          IF(IPCDER.EQ.3.AND.NFMOPCM.EQ.0) THEN
c           CALL DAREAD(IDAF,IODA,X(LDAI),LFIXPV,334,0)
c           CALL DAREAD(IDAF,IODA,X(LIDDAI),LFIXPV2,335,1)
c          ENDIF
C
           CALL VCLR(FEQWRK,1,3*NATM)
           DUM    = 0
C
c          GOPARR = .false.
C
c          CALL CAVITX(NUMDIF,X(LAXYZCT),X(LAXYZCT+MXTS),
           CALL CAVITX(1,X(LAXYZCT),X(LAXYZCT+MXTS),
     *            X(LAXYZCT+MXTS*2),
     *            X(LAXYZCT+MXTS*3),X(LNVERT),X(LXYZRE),X(LXYZRE+MXSP),
     *            X(LXYZRE+MXSP*2),X(LXYZRE+MXSP*3),X(LSSFE),X(LISPHE),
     *            X(LRINPCM),X(LDAI),X(LIDDAI),METHOD,FEQWRK,
     *            X(LALPCM),1)
c          
           goparr = savgop
c

           IF(IVIB.EQ.1) CALL DCOPY(natm*3,FEQWRK,1,FCM(1,JXYZ),1)
           IF(IVIB.EQ.2) 
     *       CALL DAXPY(3*NATM,-1.0D+00,FEQWRK,1,FCM(1,JXYZ),1)

           if(nfg.ne.0) then
            fmoc(ixyz,iat) = fmoc(ixyz,iat) - D(IVIB)
           else
            C(IXYZ,IAT)= C(IXYZ,IAT) - D(IVIB)
           end if
C
 370      CONTINUE
 380    CONTINUE
 390  CONTINUE
C
      IF(GOPARR) then
        IF(NXT) CALL DDI_DLBRESET
        CALL DDI_GSUMF(2461,FCM,9*NATM*NATM)
      end if
C
      DUM = 1.0D+00/(VIBSIZ*NVIB)
      DO 140 IAT = 1,NATM*3
         DO 130 JAT = 1,NATM*3
            FCM(IAT,JAT) = DUM*FCM(IAT,JAT)
  130    CONTINUE
  140 CONTINUE
      DO 590 IAT = 2,NATM*3
        DO 580 JAT = 1,IAT-1
          AVE = (FCM(IAT,JAT)+FCM(JAT,IAT))*0.5D+00
          FCM(JAT,IAT) = AVE
          FCM(IAT,JAT) = AVE
  580   CONTINUE
  590 CONTINUE
C
c     write(6,*) "wwww end numerical"
C
      CALL MAKCVM(0)
      CALL TESIN(NATM,X(LXYZRE),X(LXYZRE+MXSP),X(LXYZRE+MXSP*2),
     *           X(LALPCM),X(LRINPCM),X(LINIPCM),
     *           X(LMEPCM))
      CALL PEDRAM
c     IF(IPCDER.EQ.3.AND.NFMOPCM.EQ.0) THEN
c      CALL DAREAD(IDAF,IODA,X(LDAI),LFIXPV,334,0)
c      CALL DAREAD(IDAF,IODA,X(LIDDAI),LFIXPV2,335,1)
c     ENDIF
C
      NPRINT = NPRSAV

      RETURN
      END
C
C
C*MODULE PCMHSS  *DECK DISRPHSS
C>    @brief    Second order derivative for PCM for FMO
C>
C>    @details  Secand order Derivative of G^dis
C>
C>    @author   Hiroya Nakata
C
C
      SUBROUTINE DISRPHSS(NAT,C,DKT,RWT,NT,RDIFF,EPSI,SIGMA,DKA,RWA,
     *                  NTA,NTB,XCTS,YCTS,ZCTS,AS,NVERT,XE,YE,ZE,RE,
     *                  SSFE,ISPHE,RIN,ICLAV,ILJ,N,RHO,NUMDIF,DAI,
     *                  IDDAI,FEQWRK,ALPHA,DSUMTD,FCM,DSUMTR,METHOD)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C
      LOGICAL GOPARR, DSKWRK, MASWRK,nxt,SAVGOP
C
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
c     PARAMETER (MXPT=2000, MXFGPT=12000)
C
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON /FMCOM / X(1)
      COMMON /MACHIN/ NWDVAR,MAXFM,MAXSM,LIMFM,LIMSM
      COMMON /PCMDIM/ MXSP,MXTS,MEMPCM1,MEMPCM2,NTS
      COMMON /PCMOPT/ RABI,RASC,REFPOL,THRSLS,DENSLS,WB,WA,ETA2,GD,EVAC,
     *                RHOW,PM,AREATL,AREAKP,BONDRY,OMEGA,RET,FRO,EPSINF,
     *                EPS,DR,RSOLV,VMOL,TCE,STEN,DSTEN,CMF,TABS,IDIRCT,
     *                IPCDER,IDP,ICOMP,IFIELD,ICAV,IDISP,IPRINT,IRETCAV,
     *                ICENT,IFAST,NEVAL,IEFPOL,KEEPSM,IMGABI,IMGASC,NADD
      COMMON /PCMPNT/ LXYZRE,LSSFE,LLIST,LALPCM,LRINPCM,LINAPCM,LINFPCM,
     *                LINIPCM,LMEPCM,LAXYZCT,LXYZCT2,LNVERT,LQSN,LQSE,
     *                LQFS,LQIND,LISPHE,LVAD,LQOR,LVECMUL,LAIPRJ,LFIPRJ,
     *                LPEL,LPCMCDR,LDAI,LIDDAI,LQSND,LQSED,LEPSPCM
      COMMON /OUTPUT/ NPRINT,ITOL,ICUT,NORMF,NORMP,NOPK
      COMMON /PCMPAR/ IPCM,NFT26,NFT27,IKREP,IEF,IP_F,NFMOPCM,IHET
C
      DIMENSION C(3,*),EPSI(NTA),SIGMA(NTA),DKA(NTA),RWA(NTA),
     *          DKT(NTB),RWT(NTB),NT(NTB),RDIFF(NTB),XCTS(*),YCTS(*),
     *          ZCTS(*),AS(*),NVERT(*),XE(*),YE(*),ZE(*),RE(*),SSFE(*),
     *          ISPHE(*),ALPHA(*),METHOD(*)
      DIMENSION DSUMTD(3,NTA),DSUMTR(3,NTA)
      DIMENSION DAI(3,20,MXTS),IDDAI(21,MXTS)
      DIMENSION FEQWRK(3,*)
      DIMENSION FCM(nat*3,*)
      DIMENSION D(2)
      DIMENSION RIN(*)
C
C
      NPRSAV = NPRINT
      NPRINT = 817
C
      NXT    = IBTYP.EQ.1
      NEXT   = -1
      KOUNT  = -1
C
C
      CALL VCLR(fcm,1,nat*nat*9)
C
C
      NVIB   =  2 
      VIBSIZ =  0.000001D+00
c     VIBSIZ =  0.00001D+00
      D(1)   =  VIBSIZ
      D(2)   = -VIBSIZ
      LFIXPV = 0
      IF(IPCDER.EQ.3) LFIXPV = 3*20*MXTS
      LFIXPV2 = 0
      IF(IPCDER.EQ.3) LFIXPV2= 21*MXTS/NWDVAR + 1
C     Numerical derivatives of dC/da  :  (d^2 C/da db) 
c     write(6,*) "wwww begin ",natm,mxsp
      if(numdif.eq.-1234) write(6,*) numdif
      SAVGOP = GOPARR
C
      DO 390 IVIB = 1,NVIB
        DO 380 IAT = 1,NAT
          DO 370 IXYZ = 1,3
C
           IF(GOPARR) THEN
             KOUNT = KOUNT + 1
             IF(NXT) THEN
               IF(KOUNT.GT.NEXT) CALL DDI_DLBNEXT(NEXT)
               IF(NEXT.NE.KOUNT) GOTO 370
             ELSE
               IF(MOD(KOUNT,NPROC).NE.ME) GOTO 370
             ENDIF
           ENDIF
C
C
           JXYZ       = IAT * 3 - 3 + IXYZ
           C(IXYZ,IAT)= C(IXYZ,IAT) + D(IVIB)
C
           CALL MAKCVM(1)
c          CALL TESIN(NAT,XE,YE,ZE,alpha,RIN,X(LINIPCM),
c    *                X(LMEPCM))
c          write(6,*) "enter pedram"
           GOPARR = .false.
           CALL PEDRAM
C          
c          IF(IPCDER.EQ.3.AND.NFMOPCM.EQ.0) THEN
c           CALL DAREAD(IDAF,IODA,DAI,LFIXPV,334,0)
c           CALL DAREAD(IDAF,IODA,IDDAI,LFIXPV2,335,1)
c          ENDIF
C
           CALL VCLR(FEQWRK,1,3*NAT)
           DUM    = 0
C
c          GOPARR = .false.
c          
           CALL DISRPX(NAT,C,DKT,RWT,NT,RDIFF,EPSI,SIGMA,DKA,RWA,
     *            NTA,NTB,XCTS,YCTS,ZCTS,AS,NVERT,XE,YE,ZE,RE,
     *            SSFE,ISPHE,RIN,ICLAV,ILJ,N,RHO,1,DAI,
     *            IDDAI,FEQWRK,ALPHA,DSUMTD,DSUMTR,METHOD,1)
c    *            SSFE,ISPHE,RIN,ICLAV,ILJ,N,RHO,NUMDIF,DAI,
c
           goparr = savgop
C
           IF(IVIB.EQ.1) CALL DCOPY(nat*3,FEQWRK,1,FCM(1,JXYZ),1)
           IF(IVIB.EQ.2) 
     *       CALL DAXPY(3*NAT,-1.0D+00,FEQWRK,1,FCM(1,JXYZ),1)
C
           C(IXYZ,IAT)= C(IXYZ,IAT) - D(IVIB)
C
 370      CONTINUE
 380    CONTINUE
 390  CONTINUE
C
      IF(GOPARR) then
        IF(NXT) CALL DDI_DLBRESET
        CALL DDI_GSUMF(2462,FCM,9*NAT*NAT)
      end if
C
      DUM = 1.0D+00/(VIBSIZ*NVIB)
      DO 140 IAT = 1,NAT*3
         DO 130 JAT = 1,NAT*3
            FCM(IAT,JAT) = DUM*FCM(IAT,JAT)
  130    CONTINUE
  140 CONTINUE
      DO 590 IAT = 2,NAT*3
        DO 580 JAT = 1,IAT-1
          AVE = (FCM(IAT,JAT)+FCM(JAT,IAT))*0.5D+00
          FCM(JAT,IAT) = AVE
          FCM(IAT,JAT) = AVE
  580   CONTINUE
  590 CONTINUE
C
c     write(6,*) "wwww end numerical"
C
      CALL MAKCVM(0)
      CALL TESIN(NAT,XE,YE,ZE,alpha,RIN,X(LINIPCM),
     *           X(LMEPCM))
      CALL PEDRAM
      IF(IPCDER.EQ.3.AND.NFMOPCM.EQ.0) THEN
       CALL DAREAD(IDAF,IODA,DAI,LFIXPV,334,0)
       CALL DAREAD(IDAF,IODA,IDDAI,LFIXPV2,335,1)
      ENDIF
C
      NPRINT = NPRSAV
C
      RETURN
      END
