C 17 Oct 19 - HN,YN,VQV,DSK,VM,DGF - changes for FMO 5.4
C  6 Jun 18 - HN,YN,DGF - changes for FMO 5.3
C  1 Apr 16 - TN,HN,DGF - changes for FMO 5.2
c 22 Oct 14 - HN,YN,DGF,KRB,NM,SRP - changes for FMO 5.1
C 21 May 13 - DGF,HN,TN - changes for FMO 5.0
C 19 Oct 12 - MWS - synchronize FRGINF common
C 24 Jul 12 - HN,DGF - code update to finish FMO 4.3
C 23 MAR 12 - DGF,CHC - code update to finish FMO 4.2
C 28 DEC 11 - DGF,TN,CS  - FMO 4.2 and EFMO CHANGES
C 11 Aug 11 - HL  - fix for LJ term
C 15 Apr 11 - DGF,TN - frozen domain and analytic gradient related changes
C 12 Dec 10 - DGF - esdderz: fix multilayer logic for detached bonds
C  1 Oct 10 - DGF - added active atom input IACTAT to OPTFMO
C 11 Aug 10 - DGF,TN - changes for FMO 4.0 
C 25 Mar 10 - DGF - gradient projection args changed
C 14 Oct 09 - DGF - changes for FMO 3.3
C 22 May 09 - DGF - fix for MCP/FMO gradients
C  1 May 09 - NMT - new dimension for EFPBAS
C 15 Dec 08 - TN,DGF - various changes for FMO 3.2 release
C 20 Aug 07 - DGF - various changes for FMO 3.1 release
C 24 Mar 07 - MWS - pad FRGINF common block
C 22 Dec 06 - DGF - various changes for FMO 3.0 release
C  6 Nov 06 - MWS - adjust GDDI common block
C 14 Nov 05 - DGF - pad common block ENRGYS
C 19 Sep 05 - IA  - synchronize FRGINF common
C  5 Jul 05 - MWS - SELECT NEW ATOM,BASIS,EFP,PCM,DAF DIMENSIONS
C  1 Jun 05 - DGF - fixes for the 2nd release
c 15 Mar 05 - DGF - major changes for the second release
C  7 Mar 05 - IA  - fix common block FRGINF
C 13 Feb 05 - MWS - pad common block NSHEL, FRGINF
C  5 Feb 05 - MWS - pad common FMORUN
C 23 Jul 04 - MWS - GMSXYZ: use of character strings to obey f77 rules
C 19 May 04 - DGF,KK - implement Fragment Molecular Orbital (FMO) method
c
C*MODULE fmogrd  *DECK optfmox
      SUBROUTINE optfmox
      use mx_limits, only: mxfrg,mxfgpt,mxdfg,mxpt,mxshef,mxgefp,
     *     mxrt,mxatm,mxao
      USE DFTBPB_MOD,ONLY: PERIOD,LATOPT,NDIMPBC,PBCBOX,STRFMO,
     *    DFTB_LATREC,VOLPBC,LATOPA
      USE comm_PAULMO
      USE comm_EFPCT
      USE comm_FRGLJP
      USE comm_FRGTYP
      USE comm_EFPBAS
      USE comm_EFPOTD
      USE comm_FGRAD
      USE comm_FRGINF
      USE comm_POLTMP
      USE comm_EFPFMO
      USE HARMONIC_RESTRAINTS
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (nnam=19,tiny=1.0D-06,zero=0.0D+00,one=1.0D+00,
     *           three=3.0D+00,UNITS=0.52917724924D+00,diagval=THREE)
      LOGICAL DOEFPOPT,PACK2E,GOPARR,DSKWRK,MASWRK,CVGED,DSKsav,newp,
     *        qmmm,mmonly,imomm,simomm,ISGDDI,PAROUT,INITGDDI,DOSCZV,
     *        wasgddi,MLGDDI,CHGEMB,OPTLAT,DFTBFL,SCC,SRSCC,DFTB3,DAMPXH
     *       ,LCDFTB,saveasc,usefrac
      character*72 pdbcom
      character*8 cdum
      character*4 atnam
      Integer ddi_world,ddi_group
      Parameter(ddi_world=0,ddi_group=1)
      COMMON /DFTB  / DFTBFL,SCC,SRSCC,DFTB3,DAMPXH,LCDFTB
      COMMON /ENRGYS/ ENUCR,EELCT,ETOT,SZ,SZZ,ECORE,ESCF,EERD,E1,E2,
     *                VEN,VEE,EPOT,EKIN,ESTATE(MXRT),STATN,EDFT(3),EDISP
      COMMON /FMCOM / X(1)
      COMMON /GDDI/   ISCOPE,NGROUPS,MYGROUP,MEGLOB,npglob,nnglob,JBTYP,
     *                ISGDDI,PAROUT,INITGDDI,wasgddi,MLGDDI,NSUBGR,
     *                MeUniv,NPUniv,numdlb,myworld,nworlds,mogddi
      COMMON /HAMCON/ FHMCON(MXATM),SHMCON(MXATM),NHAMCON,
     *                IHMCON(3*MXATM)
      common /INT2IC/ NINTIC,iNINTIC,nxxic,lbufpic,lixic,labsix,nintix
      COMMON /INTFIL/ NINTMX,NHEX,NTUPL,PACK2E,INTTYP,IGRDTYP
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON /MACHIN/ NWDVAR,MAXFM,MAXSM,LIMFM,LIMSM
      COMMON /MACHSW/ KDIAG,ICORFL,IXDR,modio,mem10,lpnt10,mem10m
      COMMON /OPTGRD/ XYZ(3*MXATM),ES,FE(20),
     *                CONVF,FMAXT,DXMAXT,RMAX,RMIN,RLIM,
     *                EIGMAX,EIGMIN,GRDERR,FRMS,FMAX,TRMAX,TRMIN,
     *                IC(20),NSTEP,NSERCH,NPMAX,NP,IFOLOW,
     *                NNEG,IUPHSS,IEXIT,ITRUPD,IPAD,kdiagf,NPRICO
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /PCKLAB/ LABSIZ
c     COMMON /PCMDIM/ mxsp,mxts,mempcm1,mempcm2,NTS
      COMMON /PCMITR/ RCUT(2),THRES,IPCMIT,IMUL,MXDIIS,NREG,MXITR1,
     *                MXITR2,MODPAP
      COMMON /PCMOPT/ RABI,RASC,REFPOL,THRSLS,DENSLS,WB,WA,ETA2,GD,EVAC,
     *                RHOW,PM,AREATL,AREAKP,BONDRY,OMEGA,RET,FRO,EPSINF,
     *                EPS,DR,RSOLV,VMOL,TCE,STEN,DSTEN,CMF,TABS,IDIRCT,
     *                IPCDER,IDP,ICOMP,IFIELD,ICAV,IDISP,IPRINT,IRETCAV,
     *                ICENT,IFAST,NEVAL,IEFPOL,KEEPSM,IMGABI,IMGASC,NADD
      COMMON /PCMSAV/ mxtssav,ntssav(3),lisphsav,lqpcmsav
      COMMON /PCMPAR/ IPCM,NFT26,NFT27,IRPPCM,IEFPCM,IP_F,nfmopcm,IHET
      COMMON /RUNLAB/ TITLE(10),ANAM(MXATM),BNAM(MXATM),BFLAB(MXAO)
      COMMON /TINOPT/ mparti,MMONLY, QMMM
      common /QMMM1/ IMOMM,SIMOMM,NPAIR,NSEQ,CHGEMB
      COMMON /ZMTALT/ NZMAT2,NZVAR2,NVAR2,NZMTRD,ICOORD
      COMMON /FMODTB/ LGAMMA,LZREF,LISPE,LSHIFTG,LSHIFTCT,
     *                LGAMMAD,MODESD,MODGAMMA,LDCTMP,LGAMMA3,LGAMMA3D,
     *                LNBOND,LCTMUL,LESPDFTB,LESPASC,NDFTBRST,LDFTBRST
      Common /fmoinf/ nfg,nlayer,natfmo,nbdfg,naotyp,nbody,nsegm
      common /fmoopt/ espsca(9),RESPAP(2),rESPPC(2),rESDIM,restri(4),
     *                rcorsd,respct,convfg,cnvdmp,coroff,rflmo(4),
     *                orshft,orshft2,cnvafo,ascreen(4),IXESP,mxitfg,
     *                nguess,NBSSE,modorb,modpar,irststp,irstlay,nprfmo,
     *                nfmopal,modprp,maxl1c,ipieda,modgrd,modesp,ivmul,
     *                modlmo,nopden,mofock,modfd,modfmm,ncentm,ndualb,
     *                ngab,modpan
      common /fmopnt/ lichfg,lmulfg,lidmrec,lfrgnam,llayfrg,lindat,
     *                lnCBS,lfmozan,lfmoc,lfmomas,lizbas,liaglob,libdgh,
     *                liabdfg,ljabdfg,lnCAO,lidxCAO,liaprjo,ljaprjo,
     *                lCoreAO,lOccCor,lshiftb,liodfmo,lfmoda,lfmodb,
     *                lfmoespa,lfmoespb,llocfmo,lscffrg,lfmoscf,lrij,
     *                lpopmul,lpopmat,lialoc,lindbd,liatfrg,lindfrg,
     *                lindgfrg,lnatfrg,lnat0frg,lianfrg,lzanfrg,lcfrg,
     *                llibish,llibnsh,llibng,lindatg,lfmobuf(3),lfmode,
     *                lnumfrg,lloctat,liaoglob,lloadm,lfmoge,ldgrid,
     *                liodcfmo,ljob2grp,lfmopg,lemocdr,luntxyz,luntrot,
     *                lstonep,lmapsu,lfrgmul,lclmo,lialmo,lindlmo,
     *                latclmo,llmobdf,lfgflmo,lnfglmo,llfglmo,lpfglmo,
     *                LPOPDMAT,lidmpnt,liddpnt,livmpnt,liactfg,lcrfrg,
     *                lzlmfrgv,lYlmfrgv,lndtfrg,lf_mm,lg_mm,lmaxl30,
     *                libuffg,lindatp,lsmon,lsdim,ledimfed,lexcit3d,
     *                lifgfret
      COMMON /SCZVEC/ IDFZVEC,LIPTZVEC,LIPTLG,icursczv
      dimension KQNAM(nnam),QNAM(nnam),dum3(3),nfgd(2),cellpar(6)
      DATA KQNAM /3,3,1,5,1, 5,5,1,3,3, 3,21,1,1,1, 1,1,3,3/
      DATA QNAM/8HOPTTOL  ,8HSTEP    ,8HNSTEP   ,8HMETHOD  ,8HLINSER  ,
     *          8HUPDATE  ,8HHESS    ,8HIFREEZ  ,8HSTPMIN  ,8HSTPMAX  ,
     *          8HSTPFAC  ,8HNFGD    ,8HIACTAT  ,8HNPRICO  ,8HIREST   ,
     *          8HMAXNAT  ,8HIHMCON  ,8HFHMCON  ,8HSHMCON  /
     *     optfmo/8HOPTFMO  /,steep/8HSTEEP   /,cg/8HCG      /,
     *     hssupd/8HHSSUPD  /,bfgsl/8HLBFGS   /,bfgs/8HBFGS    /,
     *     dfp/8HDFP     /,rnone/8HNONE    /,GUESS/8HGUESS   /,
     *     atnam/'    '/,cdum/' '/
c
c     Main driver for genuine FMO optimisations.
c     parstat: GlobalNone
c
C        ----- READ NAMELIST $OPTFMO -----
C
c     Scan $FMOXYZ to find natfmo.
c
      qmmm=.false.
      mmonly=.false.
      optlat=.false.
      ldim=0
      if(dftbfl.and.period.and.latopt.ne.0) then
        optlat=.true.
        LDIM=NDIMPBC
      endif
      ldim3=LDIM*3
      units3=units**3
c     write(6,*) 'wwwpbc',optlat,LDIM,latopt
c
      CALL SEQREW(IR)
      CALL FNDGRP(IR,' $TINKEY',JEOF1)
      CALL SEQREW(IR)
      CALL FNDGRP(IR,' $LINK  ',JEOF2)
      if(JEOF1.eq.0 .and. JEOF2.eq.0) qmmm=.true.
      if(JEOF1.eq.0 .and. JEOF2.ne.0) mmonly=.true.
      saveasc=iand(modpap,256).ne.0.and.nfmopcm.eq.1.and.ICOMP.eq.0
c     saveasc could be adapted for PCM[2] with 2 sets of charges [1] and [2]. 
c     ICOMP.eq.2 would require storing two sets of charges (e+n).
c     However, there is no gradient for ICOMP.eq.2 now.
c     if(maswrk) write(6,*) 'wwwsave',saveasc,modpap,nfmopcm
c
      if(qmmm. or. mmonly) then
         CALL TINKIN(IR,IW)
c        force tinker coordinates to pass to LINKIN
         if(qmmm) ICOORD=5
         if(qmmm) CALL LINKIN(1)
         CALL TOPTIN
         if(qmmm) natfmo=nseq+npair
         if(mmonly) natfmo=nqmmmatoms()
      else
         call fmoxyz(' $FMOXYZ',1,natfmo,dum,dum3,idum)
      endif
c
      nc1=natfmo*3
C 
      DOEFPOPT = IEFPFMO.NE.0.AND.MOVE.NE.0
      NCVAL = NC1
      IF (DOEFPOPT) NC1=NC1+NFRG*6
      indlat=nc1
      if(optlat) nc1=nc1+ldim3
C
      nc2=(nc1*nc1+nc1)/2
      CALL VALFM(LOADFM)
      lifreez=LOADFM+1
      liactat=lifreez+(NCVAL-1)/nwdvar+1
      last=liactat+(natfmo-1)/nwdvar+1
      NEED2 = LAST- LOADFM -1
      CALL GETFM(NEED2)
c     liactat above does not provide space for freezing EFP atoms!
c
      opttol=1.0D-04
      nstep=200
      step=0.75D+00
      stepmin=zero
      stepmax=one
      stepfac=0.9D+00
c     Currently supported methods:
c     steep:  steepest descent, delta-x = - alpha * grad
c     cg:     Fletcher-Reeves conjugate gradient
c     lbfgs:  limited memory BGFS
c     hssupd: hessian update methods: BGFS or DFP updates supported.
      rmethod=hssupd
      update=bfgs
c
c     All of these can be combined with linear search.
c     Linear search options (LINSER):
c     0: do not perform 
c     1: perform 
c
      linser=0
c
c     At present hess can be read or guess.
      HESS=GUESS
c
c     IFREEZ is the array of coordinates to be brozen. Only Cartesians are
c     supported.
c     NOTE: for cell optimisations in DFTB/PBC, only fractional atomic
c     coordinates in the cell are frozen; atoms move along with 
c     expanding/shrinking cell.
c     It may not be the best thing to do; but now only that is implemented.
c
      KQNAM(8)=nc1*10+1
      KQNAM(13)=natfmo*10+1
      KQNAM(17)=mxatm*3*10+1
      KQNAM(18)=mxatm*10+3
      KQNAM(19)=mxatm*10+3
c     harmonic constraints are dimensioned MXATM here, and there is
c     a smaller limit MXRESTR, so the real limit is is MXRESTR per type,
c     i.e., in any case not more than MXRESTR*3 contraints.
      call viclr(x(lifreez),1,NCVAL)
      call viclr(x(liactat),1,natfmo)
      NHAMCON=0
      CALL VICLR(IHMCON,1,3*mxatm)
      CALL VCLR(FHMCON,1,mxatm)
      CALL VCLR(SHMCON,1,mxatm)
      nfgd(1)=0
      nfgd(2)=0
c     nfgd(1) must be set to nfg for using modpar=512 or 1024
c     nfgd(2) must be set to nfg for nguess=8
c     If either option is not used, set that element to 0.
      ipcord=0
      irestopt=0
c     0 disable restarts
c     1 enable restarts
      usefrac=.false.
c     This is a built-in option, not currently settable in the input.
c     usefrac is thought to be effective when atomic coordinates are frozen
c     and only lattice parameters are optimized. 
c
      maxnatd=0
      JRET = 0
      CALL NAMEIO(IR,JRET,optfmo,NNAM,QNAM,KQNAM,
     *            opttol,step,nstep,rmethod,
     *            linser,update,hess,x(lifreez),stepmin,stepmax,stepfac,
     *            nfgd,x(liactat),ipcord,
     *            irestopt,maxnatd,IHMCON,FHMCON,SHMCON,  0,0,0,0,0,
     *            0,0,0,0,0,  0,0,0,0,0,   0,0,0,0,0,  0,0,0,0,0,
     *            0,0,0,0,0,  0,0,0,0,0,   0,0,0,0,0,  0,0,0,0,0,
     *            0,0,0,0,0,  0)
      IF(JRET.EQ.2) THEN
         IF (MASWRK) WRITE(IW,*) 'Error reading $OPTFMO' 
         CALL ABRT
      END IF
      if(rmethod.ne.steep.and.rmethod.ne.cg.and.rmethod.ne.bfgsl.and.
     *   rmethod.ne.hssupd) then
        if(maswrk) write(iw,9800) rmethod
        call abrt
      endif
      if(rmethod.ne.hssupd) then
        update=rnone 
      else if(update.ne.dfp.and.update.ne.bfgs) then
        if(maswrk) write(iw,9810) update
        call abrt
      endif
      CALL CNTNZ(IHMCON,NHAMCON,MXATM)
      if(nhamcon.gt.mxatm) call abrtx("Too many constraints")
      ipcord0=ipcord
c     tollat=opttol
c
      if(rmethod.eq.hssupd) then
        nhblock=(nc2-1)/nproc+1
        lmin=me*nhblock+1
        lmax=min(lmin+nhblock-1,nc2)
c       triangular Hessian matrix is divided into nodes as a linear array
c       Node ME handles elements from lmin to lmax.
        call tribrk(lmin,-1,imin,jmin)
        imin=imin-1
        call tribrk(lmax,-1,imax,jmax)
        imax=imax-1
c       write(6,*) 'wwwints',nc1,nc2,nhblock,lmin,lmax,imin,jmin,imax,jmax
      else
        nhblock=0
      endif
      ndpnt=(nfgd(2)*nfgd(2)-nfgd(2))/2
      NDAR30=1+nfgd(1)*2*2
c     write(6,*) 'wwwbb',nfg,nfgd,NDAR30
c
      CALL VALFM(LOADFM)
      lfmoc=LOADFM+1
      lfmozan=lfmoc+nc1+3
      lfmomas=lfmozan+NCVAL/3
      lfmoge=lfmomas+NCVAL/3
      lizbas=lfmoge+nc1
      lp=lizbas+(NCVAL/3-1)/nwdvar+1
      lidmpnt=lp+nc1
      liddpnt=lidmpnt + (NDAR30-1)/nwdvar+1
c     write(6,*) 'wwwallocating in OPTFMOX',NDAR30
      last=liddpnt + (ndpnt-1)/nwdvar+1
c     DOSCZV = iand(modgrd,32).ne.0.and.iand(modpar,512+1024).ne.0
c     DOSCZV = .true.
      DOSCZV = nfgd(1).gt.0
c     write(6,*) 'wwwnpri0',nfgd,DOSCZV
      if (DOSCZV) then
        liptzvec = last
C       last     = liptzvec + 2*nfgd(1)/nwdvar + 1
        liptlg   = liptzvec + 2*nfgd(1)/nwdvar + 1
        last     = liptlg   +   nfgd(1)/nwdvar + 2
      else
        liptzvec = 0
        liptlg = 0
      end if
      if(saveasc) then
        mxtssav=natfmo*60
        lisphsav=last
        lqpcmsav=lisphsav+(mxtssav-1)/nwdvar + 1
        last=lqpcmsav+mxtssav*3
c       At most, three sets (ITSPHE is shared):
c       (1) monomer, SCC
c       (2) dimer, total
c       (3) monomer, SCZV
c       write(6,*) 'wwwalloc',mxtssav
      else
        lisphsav=last
        lqpcmsav=last
        mxtssav=0
      endif
      lfrc=last
      if(optlat) last=lfrc+indlat
      call viclr(ntssav,1,3)
      ndftbrst=maxnatd*nfgd(1)
      ldftbrst=last
      if(dftbfl.and.ndftbrst.ne.0) last=ldftbrst+ndftbrst*3
c     Allocate two sets (to use that set MAXNAT in OPTFMO):
c     1. Atomic charges for SCC, size ndftbrst*2
c     2. Atomic (other) charges for SCZV, size ndftbrst
c
      if(qmmm .or. mmonly) then
         lmmge=last
         last=lmmge+nc1
      endif
c
      if(rmethod.eq.bfgsl.or.rmethod.eq.hssupd.or.linser.ne.0) then
        lfmocp=last
        lfmogep=lfmocp+nc1
        last=lfmogep+nc1
        if(rmethod.eq.hssupd) then
          lhdg=last
          lu=lhdg+nc1
          last=lu+nc1
        endif
      endif
      NEED = LAST- LOADFM -1
      CALL GETFM(NEED)
c
c     See if we can surreptitiously use incore integral buffer to store the
c     Hessian. The code below needs to be synchronised with INT2EIC.
c     Nota bene: it takes about 113 megawords just for a 5000 atom molecule!
c
      iuse2e=0
      if(NINTIC.gt.0) then
        LABSIZ2=2/LABSIZ
        n2ebuf=NINTMX+NINTIC+(NINTMX+NINTIC-1)/LABSIZ2+1
c       write(6,*) 'Buffer',n2ebuf,nhblock
        if(n2ebuf.ge.nhblock) then
          lhess=LBUFPIC 
          iuse2e=1
        endif
      endif
      if(ixftch(x(liactat),1).eq.-1) then 
c       IACTAT(1)=-1 is a special input to freeze all atoms.
        do i=1,NCVAL
          call ixstor(x(lifreez),i,i)
        enddo
      else if(ixftch(x(liactat),1).ne.0) then
        call setactat(natfmo,natfmo,x(liactat),x(lifreez))
      endif
      usefrac=optlat
c
      nfrz=0
      do i=1,NCVAL
        if(ixftch(x(lifreez),i).ne.0) then
          nfrz=nfrz+1
        else
          goto 90 
        endif 
      enddo
   90 continue
      need1=0
      if(rmethod.eq.hssupd.and.iuse2e.eq.0.and.irestopt.eq.0) then
c       for irestopt=0 allocate here as the data have to be passed.
        CALL VALFM(LOADFM)
        lhess=LOADFM+1
        last=lhess+nhblock
        NEED1 = LAST- LOADFM -1
        CALL GETFM(NEED1)
      endif
      if(saveasc) then
        call viclr(x(lisphsav),1,mxtssav)
        call vclr(x(lqpcmsav),1,mxtssav*3)
      endif
c
      nfthes=38
      DSKsav=dskwrk
      if(rmethod.eq.hssupd.and.irestopt.gt.0) then
        dskwrk=.true.
        CALL SEQOPN(NFThes,'HESSIAN','UNKNOWN',.FALSE.,'UNFORMATTED')
        dskwrk=dsksav
      endif
c
c     Attempt ro read in restart data other than Hessian 
c
      frmsp=one
      irest=0
      if(irestopt.gt.0) then
        if(rmethod.eq.hssupd.or.rmethod.eq.bfgsl)
     *    call rstopt1(natfmo,x(lfmocp),x(lfmogep),irest)
        if(rmethod.eq.cg) call rstopt2(natfmo,frmsp,x(lp),irest)
      endif
c
      NSERCH=0
      CVGED=.false.
      modef=0
      e0=zero
      depre=zero
      newp=irest.eq.0
      if(linser.eq.0.and.rmethod.eq.hssupd) newp=.true.
      slopenew=zero
      slope0=zero
      stepold=zero
      stepnew=zero
      step0=zero
      pnorm=one
      lsiter=0
      rad2ang=9.0D+01/acos(zero)
      if(maswrk) then
        write(iw,9000) rmethod,update,hess,linser,nfrz,step,stepmin,
     *                 stepmax,stepfac,irestopt,ipcord,usefrac,
     *                 need+need1+need2
        IF(NFRZ.GT.0.and.iand(nprfmo,3).lt.3) 
     *    WRITE(IW,9010) (ixftch(x(lifreez),i),I=1,NFRZ)
        if(nhblock.ne.0) then
          if(iuse2e.eq.0) then
            WRITE(IW,9020) nhblock
          else
            WRITE(IW,9030) nhblock
          endif
        endif
      endif
c
      if(qmmm .or. mmonly) then
         call tin2fmo(natfmo,x(lfmozan),x(lfmoc),x(lizbas))
      endif
c     Here, natfmo is used because it is for a check of atom IDs.
      CALL hr_init(natfmo,ihmcon,fhmcon,shmcon)
      IF(maswrk.and.NHAMCON.GT.0) CALL hr_echo(iw)
c
c     Initialise the step size.
c
      stepi=step
c
  100 continue
c
c       Compute single point energy+gradient.
c       optimization of mm region will be completed first.
c
        if(qmmm. or. mmonly) then
           if(isgddi) call gddi_scope(ddi_world)
           call toysfmo(natfmo,x(lfmoc),x(lmmge),emm)
           if(isgddi) call gddi_scope(ddi_group)
        endif
        if(.not. mmonly) then
c         if(isgddi) call gddi_scope(ddi_world)
          call fmox(modef) 
          if(iand(modpar,512+1024).ne.0.and.nfgd(1).lt.nfg) then
            if(maswrk) write(iw,*) 'Set nfgd(1) in $OPTFMO to',nfg
            call abrt
          endif
          if(iand(nguess,8).ne.0.and.nfgd(2).lt.nfg) then
            if(maswrk) write(iw,*) 'Set nfgd(2) in $OPTFMO to',nfg
            call abrt
          endif
        endif
        IF(NHAMCON.GT.0) CALL HARMCON(natfmo,ETOT,x(lfmoge),x(lfmoc),
     *                                maswrk.and.iand(nprfmo,3).lt.3,iw)
c
        if(qmmm) then
c       sum up mm and qm gradients
           do i=1,nseq*3
              x(lfmoge+i-1)=x(lfmoge+i-1)+x(lmmge+i-1)
           enddo
        endif
        if(mmonly) then
           do i=1,natfmo*3
              x(lfmoge+i-1)=x(lmmge+i-1)
           enddo
        endif
        if(dftbfl.and.period.and.latopt.ne.0) then
          call dcopy(ldim3,strfmo,1,x(lfmoge+indlat),1)
          call latgconstr(latopt,latopa,ldim,x(lfmoge+indlat))
        endif
c
c       Perhaps, OPTLAT can also use IFREEZ for freezing some
c       lattice vectors (but not angles).
c
        do i=1,nfrz
          x(lfmoge-1+ixftch(x(lifreez),i))=zero
        enddo
        IF (DOEFPOPT) THEN
          II=NCVAL
          DO INF=1,NFRG
            X(LFMOGE+II  )=DEFT(1,INF)
            X(LFMOGE+II+1)=DEFT(2,INF)
            X(LFMOGE+II+2)=DEFT(3,INF)
            X(LFMOGE+II+3)=TORQ(1,INF)
            X(LFMOGE+II+4)=TORQ(2,INF)
            X(LFMOGE+II+5)=TORQ(3,INF)
            II=II+6
          END DO
        END IF
c
        if(rmethod.eq.hssupd.and.iuse2e.eq.0.and.irestopt.gt.0) then
          CALL VALFM(LOADFM)
          lhess=LOADFM+1
          last=lhess+nhblock
          NEED1 = LAST- LOADFM -1
          CALL GETFM(NEED1)
        endif
c
        NSERCH=NSERCH+1 
        modef=2
        NCVALa=NCVAL
        if(optlat) NCVALa=nc1
        call EGMAX(x(lfmoge),NCVALa,DEFT,TORQ,NFRG,FMAX,FRMS)
        if(nfrz.ne.0) then
c         rescale FRMS because the number of actual degrees of freedom
c         is not NCVALa, but NCVALa-NFRZ.
          FCVAL=NCVALa
c         write(6,*) 'rms before',FRMS,NFRZ,NCVALa
          FRMS=FRMS*SQRT(FCVAL/(NCVALa-NFRZ))
c         write(6,*) 'rms after',FRMS
        endif
c       if(optlat) 
c    *    call EGMAX(x(lfmoge+indlat),ldim3,DEFT,TORQ,0,FLMAX,FLRMS)
        CVGED= (FRMS .LT. opttol/THREE) .AND. (FMAX .LT. opttol)
c       if(optlat)
c    *    CVGED= CVGED.and.(FLRMS.LT.tollat/THREE).AND.(FLMAX.LT.tollat)
c       if(maswrk) write(iw,9200) NSERCH,etot,FMAX,FRMS
        if(maswrk) then
           if(qmmm) then
              etabmm=etot+emm
              write(iw,9201) NSERCH,etot,emm,etabmm,FMAX,FRMS
           elseif (mmonly) then
              write(iw,9202) NSERCH,emm,FMAX,FRMS
           else
c             if(optlat) then
c               write(iw,9205) NSERCH,etot,FMAX,FRMS,FLMAX,FLRMS
c             else
                write(iw,9200) NSERCH,etot,FMAX,FRMS
c             endif
           endif
        endif
c
        if(CVGED) goto 150
c
        if(rmethod.eq.hssupd) then
          if(NSERCH.eq.1.and.HESS.eq.GUESS) then 
            call hssini(x(lhess),imin,jmin,imax,jmax)
          else if(irestopt.gt.0) then
            dskwrk=.true.
            CALL SEQREW(NFThes)
            read(nfthes) (x(lhess+i-1),i=1,nhblock)
            dskwrk=dsksav
            if(NSERCH.eq.1.and.maswrk) write(iw,9250)
          endif
        endif
c
c       The search vector P is updated during linear search only if any 
c       is satisfied:
c       1. this is the first iteration and there was no restart 
c       2. the energy is lower than at the previous iteration
c       3. too many linear search iterations 
c       If no linear search is used, P is always updated except for the
c       first iteration with a restart.
c
        if(linser.ne.0) then
          if(lsiter.ge.linser) then
            newp=.true.
          else
            newp=etot.le.e0.or.irest.eq.0
          endif
        endif
c       write(6,*) 'wwwnewp chosen of',newp
c
        if(irest.eq.0.or.rmethod.eq.steep) then
c         Unless restart data have been read, the first step in all three 
c         methods is the steepest descent.
c
c         The search vector (P) is updated before
c         linear search iterations start.
c
          if(newp) then
            call dcopy(nc1,x(lfmoge),1,x(lp),1)
            call dscal(nc1,-one,x(lp),1)
          endif
        else if(rmethod.eq.cg) then
          if(newp) then
            call dscal(nc1,frms/frmsp,x(lp),1)
            call daxpy(nc1,-one,x(lfmoge),1,x(lp),1)
            frmsp=frms
          endif
        else 
c         Inverse Hessian update (bfgsl or hssupd).
c         Update for all points, including all of linear search. 

c         lfmocp and lfmogep hold previous coordinates and gradient. 
          call daxpy(NCVAL,-ONE,x(lfmoc),1,x(lfmocp),1)
          if(dftbfl.and.period.and.latopt.ne.0)
     *      call daxpy(ldim3,-ONE,pbcbox,1,x(lfmocp+indlat),1)
          call daxpy(nc1,-ONE,x(lfmoge),1,x(lfmogep),1)
          cgnorm=ddot(nc1,x(lfmocp),1,x(lfmogep),1)
          if(rmethod.eq.bfgsl) then
            if(newp) then
              cfnorm=ddot(nc1,x(lfmocp),1,x(lfmoge),1)
              gfnorm=ddot(nc1,x(lfmogep),1,x(lfmoge),1)
              ggnorm=ddot(nc1,x(lfmogep),1,x(lfmogep),1)
              dampc=gfnorm/cgnorm-(one+ggnorm/cgnorm)*cfnorm/cgnorm
              dampg=cfnorm/cgnorm
c
              call dcopy(nc1,x(lfmoge),1,x(lp),1)
              call dscal(nc1,-one,x(lp),1)
              call daxpy(nc1,dampc,x(lfmocp),1,x(lp),1)
              call daxpy(nc1,dampg,x(lfmogep),1,x(lp),1)
            endif
          else
c           hdg=[1/H]*dg
            call pMTAV(x(lhess),x(lfmogep),x(lhdg),nc1,imin,jmin,
     *                                                 imax,jmax)
            ghgnorm=ddot(nc1,x(lfmogep),1,x(lhdg),1)
            if(update.eq.bfgs) then
              call dcopy(nc1,x(lfmocp),1,x(lu),1)
              call dscal(nc1,one/cgnorm,x(lu),1)
              call daxpy(nc1,-one/ghgnorm,x(lhdg),1,x(lu),1)
            endif
            call updhss(update,x(lhess),imin,jmin,imax,jmax,
     *                  x(lfmocp),x(lhdg),x(lu),cgnorm,ghgnorm)
          endif
        endif
        if(rmethod.eq.hssupd.and.newp) then
c         hdg=[1/H]*g
          call pMTAV(x(lhess),x(lfmoge),x(lhdg),nc1,imin,jmin,imax,jmax)
c         p=-hdg=-[1/H]*g
          call vclr(x(lp),1,nc1)
          call daxpy(nc1,-one,x(lhdg),1,x(lp),1)
        endif
c
c       Determine step size using linear search 
c
        if(linser.eq.0) then 
c
          IF(NSERCH.gt.1)
     *      call stepofmo(etot,e0,depre,stepi,stepmin,stepmax,stepfac)
c
          stepnew=stepi
          e0=etot
c
c         Save previous point energy+gradient and
c         update the coordinates with the search vector.
c
          if(rmethod.eq.bfgsl.or.rmethod.eq.hssupd) then
            call dcopy(NCVAL,x(lfmoc),1,x(lfmocp),1)
            if(dftbfl.and.period.and.latopt.ne.0)
     *        call dcopy(ldim3,pbcbox,1,x(lfmocp+indlat),1)
            call dcopy(nc1,x(lfmoge),1,x(lfmogep),1)
          endif
        else
          if(etot.gt.e0.and.lsiter.lt.linser) then
            lsiter=lsiter+1
            ediff=etot-e0
            slopenew=ddot(nc1,x(lp),1,x(lfmoge),1)
            call stepsize(ediff,slope0,slopenew,step0,stepold,stepnew)
            stepnew=MAX(stepnew,stepmin)
c           Do not let step size go below stepmin.
c
c           check the angle between the search vector and the negative gradient
c           vector (=force)
c
            cosang = -slopenew/(FRMS*pnorm)
            cosang = min(one, max(-one,cosang))
            ang = rad2ang*acos(cosang)
            if(maswrk) write(iw,9300) lsiter,ang,stepnew,ediff
c
c           Revert to the linser.eq.0's coordinate+gradient.
c
            if(rmethod.eq.hssupd.or.rmethod.eq.bfgsl) then
              call daxpy(NCVAL,ONE,x(lfmoc),1,x(lfmocp),1)
              if(dftbfl.and.period.and.latopt.ne.0)
     *          call daxpy(ldim3,ONE,pbcbox,1,x(lfmocp+indlat),1)
              call daxpy(nc1,ONE,x(lfmoge),1,x(lfmogep),1)
            endif
            call dcopy(NCVAL,x(lfmocp),1,x(lfmoc),1)
            if(dftbfl.and.period.and.latopt.ne.0)
     *        call dcopy(ldim3,x(lfmocp+indlat),1,pbcbox,1)
            call dcopy(nc1,x(lfmogep),1,x(lfmoge),1)
c
c           Unless we reset the pre-LS step size (stepi), it will never go 
c           down and will only increase (by stepfac).
c           Therefore, if energy went up during linear search, reset stepi
c           to the default value. Note that stepi is not stepnew. The latter
c           is the current step size during LS. 
c
            stepi=step
          else
c           Use linear search to adjust the step size for energy raising.
c           and use stepofmo for energy lowering. 
            IF(NSERCH.gt.1.and.etot.lt.e0)
     *        call stepofmo(etot,e0,depre,stepi,stepmin,stepmax,stepfac)
            lsiter=0
            e0=etot
            stepnew=stepi
            pnorm=sqrt(ddot(nc1,x(lp),1,x(lp),1))
            slope0=ddot(nc1,x(lp),1,x(lfmoge),1)
            call dcopy(NCVAL,x(lfmoc),1,x(lfmocp),1)
            if(dftbfl.and.period.and.latopt.ne.0)
     *        call dcopy(ldim3,pbcbox,1,x(lfmocp+indlat),1)
            call dcopy(nc1,x(lfmoge),1,x(lfmogep),1)
          endif
          stepold=stepnew
        endif
c
c       Find the next set of coordinates.
c
        do i=1,nfrz
          x(lp-1+ixftch(x(lifreez),i))=zero
        enddo
c
c       write(6,*) 'the P vector is'
c       call prsq(x(lp),nc1,1,1)
c
c       This is an empirical fix. Guess Hessian produces too wide a step. 
c
        if(rmethod.eq.hssupd.and.NSERCH.eq.1.and.hess.eq.guess) 
     *    stepnew=stepnew/diagval
c
C       calculate the new coordinates of EFPs
C
C       ROTATE AND TRANSLATE EFP FRAGMENTS
C
        IF (DOEFPOPT) THEN
          NMST=0
          NPST=0
          NLJST=0
          NDPST=0
          NRST=0
          NDST=0
          NATST=0
          NORST=0
          NCTMOST=0
          II=NCVAL
C
          MXBF=0
          DO I=1,NFRG
             MXBF=MAX(MXBF,NPBF(I))
          END DO
          CALL VALFM(LOADFM)
          LPROVEC = LOADFM  + 1
          LCTVEC  = LPROVEC + MXBF*NTMO
          LAST    = LCTVEC  + MXBF*NTCTMO
          NEEDefp = LAST-LOADFM-1
          CALL GETFM(NEEDefp)
          LENPV=MXBF*NTMO
          IF (NFRG.GT.0.AND.NTMO.GT.0)
     *       CALL DAREAD(IDAF,IODA,X(LPROVEC),LENPV,267,0)
          LENCV=MXBF*NTCTMO
          IF (NFRG.GT.0.AND.NTCTMO.GT.0)
     *       CALL DAREAD(IDAF,IODA,X(LCTVEC),LENCV,269,0)
C
          DO I1=1,NFRG
c           XDEL=D(II)
c           YDEL=D(II+1)
c           ZDEL=D(II+2)
c           THETAX=D(II+3)
c           THETAY=D(II+4)
c           THETAZ=D(II+5)
            XDEL=  stepnew*X(LP+II)
            YDEL=  stepnew*X(LP+II+1)
            ZDEL=  stepnew*X(LP+II+2)
            THETAX=stepnew*X(LP+II+3)
            THETAY=stepnew*X(LP+II+4)
            THETAZ=stepnew*X(LP+II+5)
            X(LFMOCP+II  ) = -XDEL
            X(LFMOCP+II+1) = -YDEL
            X(LFMOCP+II+2) = -ZDEL
            X(LFMOCP+II+3) = -THETAX
            X(LFMOCP+II+4) = -THETAY
            X(LFMOCP+II+5) = -THETAZ
            II=II+6
C
C
            CALL EFDSPL(I1,XDEL,YDEL,ZDEL,THETAX,THETAY,THETAZ,
     *                  NMST,NPST,NLJST,NDPST,NRST,NDST,NATST,NORST,
     *                  X(LPROVEC),NCTMOST,X(LCTVEC),MXBF)
C
            NMST  = NMST  + NMPTS(I1)
            NPST  = NPST  + NPPTS(I1)
            NLJST = NLJST + NLJPTS(I1)
            NDPST = NDPST + NDPPTS(I1)
            NRST  = NRST  + NRPTS(I1)
            NDST  = NDST  + NDPTS(I1)
            NATST = NATST + NATEF(ISET(I1))
            NORST = NORST + NORBEF(I1)
            NCTMOST = NCTMOST + NCTMO(ISET(I1))
          END DO
C
          IF (NFRG.GT.0.AND.NTMO.GT.0)
     *       CALL DAWRIT(IDAF,IODA,X(LPROVEC),LENPV,267,0)
          IF (NFRG.GT.0.AND.NTCTMO.GT.0)
     *       CALL DAWRIT(IDAF,IODA,X(LCTVEC),LENCV,269,0)
C
          CALL RETFM(NEEDefp)
        END IF
        PI=3.14159265358979312D+00
        todeg=180/pi
c       write(6,*) 'update coord',NCVAL,nc1
        call daxpy(NCVAL,stepnew,x(lp),1,x(lfmoc),1)
c       nuplat=10
        if(dftbfl.and.period.and.latopt.ne.0) then
c    *     mod(NSERCH,nuplat).eq.0) then
c
c         calculate fractional coordinates with the new C and old PBCBOX.
c
          if(usefrac)
     *    CALL DGEMM('N','N',3,NATFMO,3,
     *               1.0D+00/(2.0D+00*PI),PBCBOX(1,4),3,x(lfmoc),3,
     *               0.0D+00,X(LFRC),3)
c
c         update PBCBOX
c
c         pbcmax=0.05D+00
c         steppbc=min(stepnew,pbcmax)
          steppbc=stepnew
          call daxpy(ldim3,steppbc,x(lp+indlat),1,PBCBOX,1)
c 
c         Get new cell parameters.
          call getcellpar(ldim,PBCBOX,cellpar)
c
          call DFTB_LATREC(.false.)
c         write(6,*) 'Updated the lattice (new)',pbcbox
c
c         With the new PBCBOX, scale new atomic coordinates.
c
          if(usefrac)
     *    CALL DGEMM('T','N',3,NATFMO,3,
     *       1.0D+00,PBCBOX(1,1),3,X(LFRC),3,
     *       0.0D+00,x(lfmoc),3)
c
c         Update the P vector with the new coordinates.
c         Should do?
c           call dcopy(NATFMO*3,x(lfmoc),1,x(lp),1)
        endif
c
        depre=stepnew*ddot(nc1,x(lp),1,x(lfmoge),1)
c       depre does not include Hessian corrections to the predicted energy
c       change.
  150   continue
c
c       Save restart data (Hessian and others).
c
        if(rmethod.eq.hssupd) then
          if(irestopt.gt.0) then
            dskwrk=.true.
            CALL SEQREW(NFThes)
            write(nfthes) (x(lhess+i-1),i=1,nhblock)
            dskwrk=dsksav
            if(iuse2e.eq.0) CALL RETFM(NEED1)
          endif
        endif
        irest=1
        if(irestopt.gt.0) then
          if(rmethod.eq.hssupd.or.rmethod.eq.bfgsl)
     *      call rstopt1(natfmo,x(lfmocp),x(lfmogep),irest)
          if(rmethod.eq.cg) call rstopt2(natfmo,frmsp,x(lp),irest)
c       else
c         irest=0
c         Keep irest=1 as otherwise the code is confused;
c         a fake "restart" here is fine because everything is stored in memory. 
        endif
c
        if(linser.eq.0) newp=.true.
c
c 150   continue
        if(maswrk) then
          if(ipcord0.eq.0) ipcord=nserch
          write(iw,9228) (TITLE(I),I=1,10)
          if(CVGED) then
            write(iw,9230)
          else
            if (mod(nserch,ipcord).eq.0) write(iw,9210) NSERCH+1
          endif
          if(iand(modio,512).ne.0) then
            nftxyz=104
            if(meglob.eq.0) then
              CALL SEQOPN(nftxyz,'XYZ','UNKNOWN',.FALSE.,'FORMATTED')
              write(nftxyz,9220) (i,int(x(lfmozan-1+i)+tiny),
     *                   (x(lfmoc-1+j+(i-1)*3)*UNITS,j=1,3),i=1,NCVAL/3)
              CALL SEQCLO(nftxyz,'KEEP')
            endif
          else
           if (mod(nserch,ipcord).eq.0.or.cvged.or.NSERCH.ge.nstep) then
c    *      write(iw,9220) (i,int(x(lfmozan-1+i)+tiny),
c    *                   (x(lfmoc-1+j+(i-1)*3)*UNITS,j=1,3),i=1,NCVAL/3)
            do i=1,NCVAL/3
              iz=int(x(lfmozan-1+i)+tiny)
              call zsymnum(cdum,atnam,iz)
              write(iw,9225) atnam,iz,(x(lfmoc-1+j+(i-1)*3)*UNITS,j=1,3)
            enddo
            endif
c           if(optlat) then
            if(dftbfl.and.period.and.latopt.ne.0) then
           call EGMAX(x(lfmoge+indlat),ldim3,DEFT,TORQ,NFRG,FMAXl,FRMSl)
              NSERCH1=NSERCH
              if(.not.cvged) NSERCH1=NSERCH+1
c             NSERCH1 is the next step after NSERCH except when converged.
c             (if converged, the cell is unchanged, so it is NSERCH).
              WRITE (IW,9260) NSERCH1,VOLPBC*units3,FMAXl,FRMSl,
     *                 (cellpar(i)*units,i=1,3),(cellpar(i)*todeg,i=4,6)
              WRITE (IW,'(14X,"X",19X,"Y",19X,"Z")')
              WRITE (IW,'(1X,3F20.10)') ((PBCBOX(J,I)*UNITS,I=1,3),
     *                                   J=1,NDIMPBC)
            endif
          endif
          IF (DOEFPOPT) THEN
            NATT = NAT
            NUMM = NUM
            NAT  = 0
            NUM  = 0
            CALL PRATM(IW,1)
            CALL PRATM(IP,1)
            NAT  = NATT
            NUM  = NUMM
          END IF
        endif
        if(CVGED.or.NSERCH.ge.nstep) goto 200
      goto 100
  200 continue
c
      if(rmethod.eq.hssupd.and.irestopt.gt.0) then
        dskwrk=.true.
        CALL SEQCLO(NFThes,'KEEP')
        dskwrk=dsksav
      endif
      if(rmethod.eq.hssupd.and.iuse2e.eq.0.and.irestopt.eq.0)
     *  CALL RETFM(NEED1)
c
      if(cvged) then
        if(maswrk) write(iw,9500) NSERCH
      else 
        if(maswrk) write(iw,9510) NSERCH
      endif
c     For puny output, write out constraint data.
      if(NHAMCON.GT.0.and.maswrk.and.iand(nprfmo,3).eq.3)
     *  CALL hr_status(iw)
      write(pdbcom,9240) NSERCH,etot,FMAX,FRMS
      call writepdb('$PDB    ',natfmo,x(lfmoc),pdbcom,title)
c
      DOSCZV = iand(modgrd,32).ne.0.and.iand(modpar,512+1024).ne.0
      IF (DOSCZV) CALL RETZVEC
      call retdd
c     if(isgddi) call gddi_scope(ddi_group)
      CALL RETFM(NEED)
      CALL RETFM(NEED2)
c
      RETURN
 9000 format(/1x,75(1H"),/1x,'FMO geometry optimisation driver:',
     *       /3x,'METHOD=',A8,' UPDATE=',A8,'   HESS=',A8,' LINSER=',I4,
     *           ' NFRZ=',I8,
     *       /3x,'STEP  =',F8.4,' STPMIN=',F8.4,' STPMAX=',F8.4,
     *           ' STPFAC=',F8.4,' IREST=',I2,
     *       /3x,'IPCORD=',I4,' USEFRC=',L1,
     *      /1x,'Using',I12,' words of memory.',/1x,75(1H"),/)
 9010 FORMAT(1X,'FROZEN COORD(S)=',10I6/(17X,(10I6)))
 9020 FORMAT(1x,I9,' words will be used to store inverse Hessian.',/)
 9030 FORMAT(1x,'Using',I9,' words of incore integral buffer to store',
     *          ' inverse Hessian.',/)
 9200 format(/1x,'NSERCH:',I5,' E=',F19.9,' GMAX=',F9.7,' GRMS=',F9.7/)
 9201 format(/1x,'NSERCH:',I5,' E(QM)=',F19.9,' E(MM)=',f19.9,
     *       ' E(QM+MM)=',F19.9,' GMAX(QM)=',F9.7,' GRMS(QM)=',F9.7/)
 9202 format(/1x,'NSERCH:',I5,' E(MM)=',F19.9,' GMAX=',F9.7,
     *        ' GRMS=',F9.7/)
c9205 format(/1x,'NSERCH:',I5,' E=',F19.9,' GMAX=',F9.7,' GRMS=',F9.7,
c    *           ' GLMAX=',F9.7,' GLRMS=',F9.7/)
 9210 format(/1x,'Cartesian coordinates for NSERCH =',I4,
     *           ' are (Angst.):')
 9220 format(I6,I4,3F19.12)
 9225 format(1x,A2,I4,3F19.12)
 9228 format(/1x,'Opt level: ',10A8)
 9230 format(/1x,'Cartesian coordinates for the CONVERGED geometry ',
     *           '(Angst.):')
 9240 format('N=',I4,' E=',F19.9,' GMAX=',F9.7,' GRMS=',F9.7)
 9250 format(/1x,'Hessian was read successfully.')
c9260 format(/1X,'NSTEP=',I8,' VOL=',F12.5,' LENGTHS=',3F8.3,
 9260 format(/1X,'NSTEP=',I5,' VOL=',F12.5,' GLMAX=',F9.7,' GLRMS=',F9.7
     *          ,' LENGTHS=',3F8.3,' ANGLES=',3F8.3,
     *      /1X,' REAL-SPACE CELL VECTORS (IN ANGSTROM)')
 9300 format(/8x,'linear search iteration ',i14,
     *       /8x,' angle (step/-gradient) ',f14.4,
     *       /8x,'       step size        ',f14.8
     *       /8x,'   energy difference    ',f14.10,
     *       /8x,38(1H.),/)
 9500 format(/1x,'Geometry optimisation converged in',I5,' steps.',/)
 9510 format(/1x,'Geometry optimisation DID NOT converge in',I5,
     *           ' steps.',/)
 9800 format(/1x,'Unknown FMO optimisation driver: ',A8)
 9810 format(/1x,'Unknown inverse Hessian update: ',A8)
      END
C*MODULE fmogrd  *DECK readmond
      SUBROUTINE readmond(d,orbxch,odexch,na,nb,l1,iodfmo,idmrec0
     *          ,uhfcal)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (one=1.0D+00)
      dimension d(*),iodfmo(*)
      logical orbxch,odexch,uhfcal
      common /fmorun/ espscf,e0scf(2),emp2s,IDAFMO,icurfg,jcurfg,kcurfg,
     *                icurlay,icurunt,nat1e,ncursh,ngau,icurpop,ifmostp,
     *                moncor,needr,modrst,norbproj,nunesp,iskipesp,
     *                IESDPPC,idoprop,mp2run,icurit,idmfmo,iddfmo,
     *                iddcur,nddleft,ivmfmo,nzmtfmo,ifmobas,itmfmo(2)
c
c     read monomer density, based on 3 ways of storing:
c     1. density
c     2. orbitals
c     3. density+orbitals
c     d thus must be allocated appropriate space (l2 or l2+l3).
c
      l2=(l1*l1+l1)/2
      l3=l1*l1
      if(orbxch) then
        if(odexch) then
          call rareads(IDAFMO,iodfmo,d,l2+l3,idmrec0,0)
        else
          lenrec=l3
          if(uhfcal) lenrec=l3+l3
          call rareads(IDAFMO,iodfmo,d(l2+1),lenrec,idmrec0,0)
          if(.not.uhfcal) then
            call DMTX2(d,d(l2+1),na,l1,l1,nb)
          else
            call DMTX2(d,d(l2+1),na,l1,l1,0)
            call DMTX2(d(l2+1),d(l2+1+l3),nb,l1,l1,0)
            call daxpy(l2,one,d(l2+1),1,d,1)
            call rareads(IDAFMO,iodfmo,d(l2+1),lenrec,idmrec0,0)
          endif
        endif
      else
        call rareads(IDAFMO,iodfmo,d,l2,idmrec0,0)
      endif
      RETURN
      END
C*MODULE fmogrd  *DECK readumond
      SUBROUTINE readumond(d,wrk,orbxch,odexch,na,nb,l1,l1i,
     *                     iodfmo,idmrec0,uhfcal)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (one=1.0D+00)
      dimension d(*),wrk(*),iodfmo(*) 
      logical orbxch,odexch,uhfcal
      common /fmorun/ espscf,e0scf(2),emp2s,IDAFMO,icurfg,jcurfg,kcurfg,
     *                icurlay,icurunt,nat1e,ncursh,ngau,icurpop,ifmostp,
     *                moncor,needr,modrst,norbproj,nunesp,iskipesp,
     *                IESDPPC,idoprop,mp2run,icurit,idmfmo,iddfmo,
     *                iddcur,nddleft,ivmfmo,nzmtfmo,ifmobas,itmfmo(2)
c
c     read monomer density, based on 3 ways of storing:
c     1. density (RHF)
c     2. orbitals (DFT or ROHF/ROKS, no UHF)
c     3. density+orbitals (MCSCF)
c     d thus must be allocated appropriate space (l2 or l2+l3).
c     wrk(l3)
c
      l2=(l1*l1+l1)/2 
      l3=l1*l1
      if(orbxch) then
        if(odexch) then
          call rareads(IDAFMO,iodfmo,d,l2+l3,idmrec0,0)
        else
          lenrec=l3
          if(uhfcal) lenrec=l3+l3
          call rareads(IDAFMO,iodfmo,d(l2+1),lenrec,idmrec0,0)
          if(icurunt.gt.0) then
            call rotmo(d(l2+1),l1,na,l1i,ncursh,wrk)
            call dcopy(l1*na,wrk,1,d(l2+1),1)
            if(uhfcal) then
              call rotmo(d(l2+1+l3),l1,nb,l1i,ncursh,wrk)
              call dcopy(l1*na,wrk,1,d(l2+1+l3),1)
            end if
c           Only occupied are rotated.
          endif
          if(.not.uhfcal) then
            call DMTX2(d,d(l2+1),na,l1,l1,nb)
          else
            call DMTX2(d,d(l2+1),na,l1,l1,0)
            call DMTX2(d(l2+1),d(l2+1+l3),nb,l1,l1,0)
            call daxpy(l2,one,d(l2+1),1,d,1)
            call rareads(IDAFMO,iodfmo,d(l2+1),lenrec,idmrec0,0)
          endif
        endif
      else
        call rareads(IDAFMO,iodfmo,d,l2,idmrec0,0)
      endif
      RETURN
      END
C*MODULE fmogrd  *DECK stuff123
      SUBROUTINE stuff123(l1ij,l1,ld3,mapi,mapj,map3i,map3j,maptmp,d2,d3
     *                   )
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      dimension mapi(*),mapj(*),map3i(*),map3j(*),maptmp(ld3,2),d2(*),
     *          d3(*)
c
c     Stuff dimer(IJ) density d2 into the trimer IJK density d3
c     d3 is accumulated (added to whatever there was).
c
c     invert mapi and mapj. mapi,mapj map dimer index to monomer, maptmp
c     maps monomer index to dimer.
c
c     l2=(l1*l1+l1)/2
c     call vclr(d3,1,l2)
      call viclr(maptmp,1,ld3*2)
      do i=1,l1ij
        mi=mapi(i)
        if(mi.ne.0) maptmp(mi,1)=i
        mj=mapj(i)
        if(mj.ne.0) maptmp(mj,2)=i
      enddo
      loop=0
      do i=1,l1
        do j=1,i
          loop=loop+1
          ii=map3i(i)
          ij=map3i(j)
          ji=map3j(i)
          jj=map3j(j)
          id=0
          if(ii.ne.0) then
            id=maptmp(ii,1)
          else if(ji.ne.0) then
            id=maptmp(ji,2)
          endif
          if(id.ne.0) then
            jd=0
            if(ij.ne.0) then
              jd=maptmp(ij,1)
            else if(jj.ne.0) then
              jd=maptmp(jj,2)
            endif
            if(jd.ne.0) then
              if(id.ge.jd) then
                ij=(id*id-id)/2+jd
              else
                ij=(jd*jd-jd)/2+id
              endif
              d3(loop)=d3(loop)+d2(ij)
            endif 
          endif 
        enddo
      enddo
      RETURN
      END
C*MODULE fmogrd  *DECK grdprj
      SUBROUTINE grdprj(de,c,rmass,nat)
      use mx_limits, only: mxfrg
      USE comm_FRGINF
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON /FMCOM / X(1)
      dimension DE(3,nat),c(3,*),rmass(*)
c
c     project out rot+tran for DFT gradient
cnb   use this sub in grd2a.src 
c     This subroutine is not parallelised(?) and is called for the 
c     whole molecule!
c
c     Note on DFT gradients: normally if one does a 2-body run with 2 fragments
c     the results are indentical to regular RHF (or whatever) to ALL printed
c     digits. This is not so for DFT gradient. The reason is the projection
c     of rotational degrees of freedom, that is done for fragments 1 and 2 in 
c     the order 2-1, whereas for conventional it is done in the order 1-2.
c     However, if one switches fragments 1 and 2 in the input then the DFT
c     gradient matches identically (with mental switching 1 and 2 back). 
c     Otherwise a difference of 1e-7 was observed.
c
      NCCF = 3*NAT + 6*NFRG
      CALL VALFM(LOADFM)
      LX   = LOADFM+ 1
      LGTOT= LX    + NCCF
      LRM  = LGTOT + NCCF
      LP   = LRM   + 3*NAT + 21*NFRG
      last = LP    + NCCF 
      NEEDP= LAST-LOADFM-1
      CALL GETFM(NEEDP)
c     call fmorams(nat,zan,TABM,X(LRM),1)
      CALL PRJGRD(DE,C,X(LX),X(LGTOT),X(LRM),X(LP),rmass,nat,3*NAT,NCCF,
     *            .false.)
      CALL RETFM(NEEDP)
      RETURN
      END
C*MODULE fmogrd  *DECK fmorams
      SUBROUTINE fmorams(nat,zan,tabm,rmass,iopt)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      parameter (zero=0.0D+00,one=1.0D+00)
      dimension zan(*),tabm(*),rmass(*)
c
c     a clone of RAMS working with an arbitrary number of atoms.
C     IOPT.EQ.0 RETURNS MASS OF EACH ATOM
C     IOPT.EQ.1 RETURNS TRIPLETS OF INVERSE SQUARE ROOT MASSES
C     IOPT.EQ.2 RETURNS MASS OF EACH ATOM IN EFMO ORDERING
C
C         RETURN THE NORMAL MASSES
C     
      if(iopt.eq.0) then
        DO 310 I=1,NAT
          RMASS(I) = tabm(int(zan(I)+1.0D-03))
  310   CONTINUE 
      elseif(iopt.eq.2) then
        DO 320 I=1,NAT
          RMASS(I) = tabm(I)
  320   CONTINUE
      else
C
C         TRIPLETS OF THE INVERSE SQUARE ROOTS OF THE
C         MASSES ARE USEFUL FOR MASS WEIGHTED CARTESIANS.
C
        K = 0
        DO 480 I = 1,NAT
          amass=tabm(int(zan(I)+1.0D-03))
          IF(AMASS.GT.ZERO) THEN
            DUM = ONE/SQRT(AMASS)
          ELSE
            DUM = ONE
          END IF
          DO 470 J = 1,3
            K = K+1
            RMASS(K) = DUM
  470     CONTINUE
  480   CONTINUE
      endif
      RETURN
      END
C*MODULE fmogrd  *DECK esdgrd
C>
C>     @brief Gradient of ES dimers
C>
C>     @details Calculate the gradient of ES dimers. 
C>
C>     @author Dmitri Fedorov
C>
      SUBROUTINE esdgrd(ilay,loadhf,jobgrp,loadm,layfrg,numfrg,iodfmo,
     *                  maxl30,iaglob,scffrg,ngrfmo,mannod,mastid,da,fa,
     *                  fb,wrk,fmode,fmopg,isdrec,orbxch,irec0,irec1,
     *                 dodistr,NQMTFG,IPTLG,ZVLAG,ZVWRK,dopbcmd,irmdfmo)
      use mx_limits, only: mxatm
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      Parameter (maxpst=10,half=0.5D+00,one=1.0D+00)
      integer ddi_world,ddi_group,ddi_masters
      Parameter(ddi_world=0,ddi_group=1,ddi_masters=2)
      logical isgddi,parout,INITGDDI,GOPARR,DSKWRK,MASWRK,orbxch,
     *        iodexch,myjob,dodistr,urohf,wasgddi,uhfcal,dopbcmd,MLGDDI
      COMMON /GDDI/   ISCOPE,NGROUPS,MYGROUP,MEGLOB,npglob,nnglob,JBTYP,
     *                ISGDDI,PAROUT,INITGDDI,wasgddi,MLGDDI,NSUBGR,
     *                MeUniv,NPUniv,numdlb,myworld,nworlds,mogddi
      COMMON /GRAD  / DE(3,MXATM)
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /RUNOPT/ RUNTYP,EXETYP,NEVALS,NGLEVL,NHLEVL
      COMMON /WFNOPT/ SCFTYP,VBTYP,DFTYPE,TDDFTYP,CITYP,CCTYP,
     *                MPLEVL,MPCTYP
      Common /fmoinf/ nfg,nlayer,natfmo,nbdfg,naotyp,nbody,nsegm
      common /fmoopt/ espsca(9),RESPAP(2),rESPPC(2),rESDIM,restri(4),
     *                rcorsd,respct,convfg,cnvdmp,coroff,rflmo(4),
     *                orshft,orshft2,cnvafo,ascreen(4),IXESP,mxitfg,
     *                nguess,NBSSE,modorb,modpar,irststp,irstlay,nprfmo,
     *                nfmopal,modprp,maxl1c,ipieda,modgrd,modesp,ivmul,
     *                modlmo,nopden,mofock,modfd,modfmm,ncentm,ndualb,
     *                ngab,modpan
      common /fmorun/ espscf,e0scf(2),emp2s,IDAFMO,icurfg,jcurfg,kcurfg,
     *                icurlay,icurunt,nat1e,ncursh,ngau,icurpop,ifmostp,
     *                moncor,needr,modrst,norbproj,nunesp,iskipesp,
     *                IESDPPC,idoprop,mp2run,icurit,idmfmo,iddfmo,
     *                iddcur,nddleft,ivmfmo,nzmtfmo,ifmobas,itmfmo(2)
      DIMENSION IPTLG(1),NQMTFG(1),ZVLAG(1),ZVWRK(1)
      dimension jobgrp(*),loadm(*),layfrg(*),numfrg(*),iodfmo(*),
     *          maxl30(*),iaglob(*),scffrg(*),ngrfmo(maxpst,*),mannod(*)
     *         ,mastid(*),da(*),fa(*),fb(*),wrk(*),fmode(3,natfmo,*),
     *          fmopg(3,*),isdrec(*)
      data RMC/8HMCSCF   /,UHF,ROHF/8HUHF     ,8HROHF    /,
     *     rnone/8HNONE    /
c
c     finish computing derivatives of separated dimers
c     (Lagrangian or density derivative terms). 
c     parstat: GlobalFull+GlobalFull
c
c     Here is how the separated dimer gradient is computed. All contributions
c     are scattered over all nodes!
c     <mu|V|nu> derivative contributions are computed in TVDER
c     (both derivatives of mu,nu and V) (ifmostp=2).
c     Enr (nuclear repulsion) derivative contributions are computed in
c     ESDNRDER (ifmostp=6).
c     (mu,nu|ro,sigma) derivative contributions are computed in ESD2DER 
c     (ifmostp=6).
c     Dmu,nu (density) derivative contributions are computed here (ESDGRD)
c     (ifmostp=7).
c     jobgrp is reused here with no connexion to previous usage.
c
c     write(6,*) 'Enter esdgrd'
      ifmostp=7
      icurlay=ilay
      jcurfg=0
      kcurfg=0
c     WRITE(IW,9086)
c     CALL EGOUT(fmode(1,1,2),natfmo)
c
c     Read, exchange and save monomer potentials
c     (that is, V(I)=sum(J)V(IJ), where V(IJ) is potential for monomer I
c     due to external field of separated monomer J). 
c     Those are scattered over nodes after edimer. 
c     Note that we still are and should be in the (separated) dimer group 
c     division.
c
c     First find out if some fragments have no separated external field. 
c     (in which case jobgrp(ifg) is 0).
c
      if(maswrk.or.dodistr) then
        do ifg=1,nfg
          nrec=1+ifg+nfg
          jobgrp(ifg)=isdrec(nrec)
        enddo 
      end if
c     write(6,*) 'wwwaaa',(jobgrp(i),i=1,nfg)
c
      if(isgddi.and..not.dodistr) then
         call gddi_scope(ddi_masters)
         if(maswrk) then 
           call ddi_gsumi(2419,jobgrp,nfg)
           do 100 ifg=1,nfg
             if(layfrg(ifg).lt.ilay.and.modfd.eq.0 .or.jobgrp(ifg).eq.0)
     *         goto 100
             l1=iand(numfrg(ifg),65535)
             l2=(l1*l1+l1)/2
c            The potentials are duplicated on slaves (within a group). 
c            The amount of communications here is not negligible and some 
c            more thoughts later may be helpful. (mainly it is the dilemma:
c            whether to sum on masters only and then broadcast or global
c            sum everywhere substituting zeroes on slaves. At present the
c            Gordian knot is cut the first way). Note that masters need to 
c            broadcast not everything, but only what their groups will do:
c            thus this solution seems agreeable -> unfortunately this glorious
c            approach is flawed: masters can become slaves after regrouping
c            which would lead to a total mess, thus the solution is to plain
c            broadcast by masters after global summing here.
             nrec=1+ifg+nfg
             if(isdrec(nrec).eq.0) then
               call vclr(fa,1,l2)
             else
               CALL rareads(IDAFMO,iodfmo,fa,l2,nrec,0)
             endif
c            Global sum potentials over masters.
             call ddi_gsumf(2419,fa,l2)
             imxl30=maxl30(ifg)
             CALL rawrites(IDAFMO,iodfmo,fa,imxl30,l2,nrec,0)
  100      continue 
         endif
         call gddi_scope(ddi_group)
      endif
c
c     Broadcast now, before it is too late! (before fmopre)
c  
      if(dodistr) then
        if(isgddi) call gddi_scope(ddi_world)
        call ddi_gsumi(2421,jobgrp,nfg)
        if(isgddi) call gddi_scope(ddi_group)
      else
c
        CALL DDI_BCAST(2421,'I',jobgrp,nfg,MASTER)
c
c     write(6,*) 'jobgrp=',(jobgrp(i),i=1,nfg)
      uhfcal=.false.
      if(uhfcal) write(6,*) "uhfcal"
      do 120 ifg=1,nfg
        if(layfrg(ifg).lt.ilay.and.modfd.eq.0 .or.jobgrp(ifg).eq.0) 
     *    goto 120
        l1=iand(numfrg(ifg),65535)
        l2=(l1*l1+l1)/2
        nrec=1+ifg+nfg
        if(maswrk) CALL rareads(IDAFMO,iodfmo,fa,l2,nrec,0)
        CALL DDI_BCAST(2421,'F',fa,l2,MASTER)
        imxl30=maxl30(ifg)
        if(.not.maswrk) CALL rawrites(IDAFMO,iodfmo,fa,imxl30,l2,nrec,0)
        if(scffrg(ifg).eq.uhf) uhfcal=.true.
  120 continue 
      endif
c
c     switch GDDI groups to monomers (using monomer group division).
c
      iilay=ilay
      if(modfd.ne.0) iilay=1
c     IF(MASWRK)WRITE(6,*)'SRP: calling FMOPRE FMOGRD'
c     CALL FLSHBF(6)

      if(isgddi.and.ngrfmo(1,iilay).ne.0)
     *  call fmopre(iilay,1,ngrfmo,mannod,mastid,irmdfmo)
c     save fmoc
c---  if(dopbcmd) call genfrgmic(0,-1)
      if(dopbcmd) call genimagefg(0,0)
c
c     Now add contributions to the gradient.
c     Obtain Lagrangians by transforming potentials with densities:
c     L(I)=D(I) V(I) D(I)-dagger 
c     Finally contract with derivative overlaps:
c     grad E=Tr(L*gradS).
c
      if(isgddi) call GDDICOUNT(-1,lgroup,myjob)
c
      do 200 iifg=1,nfg
        ifg=iifg
        if(loadhf.eq.1) ifg=loadm(iifg)
        if(layfrg(ifg).lt.ilay.and.modfd.eq.0 .or.jobgrp(ifg).eq.0) 
     *    goto 200
        icurfg=ifg
        if(isgddi) then
          call GDDICOUNT(0,lgroup,myjob)
          if(.not.myjob) goto 200
        endif
cpbc
c       generate fragments around the fragment ifg accoording to
c       the minimum image convention. fmoc, cfrg, and untxyz are
c       destroyed (fmoc is stored in fmoctmp), and so they must
c       be restored later using fmoctmp.
c---    if(dopbcmd) call genfrgmic(ifg,0)
        if(dopbcmd) call genimagefg(2,ifg)
c       if(dopbcmd) then
c         call makeindexiu(ilay,loadhf,loadm,x(lindat),x(lfmozan),
c    *                   x(lfmoc),x(liabdfg),x(ljabdfg),x(lscffrg),
c    *                   1,x(lrij),rwrk,mfg2,vdwrad,
c    *                   nprfrg,.false.,nimgcell,natfmo+nbdfg,
c    *                   x(luntxyz),units,x(lindxiu),.false.)
c       end if
c
        if(maswrk.and.iand(nprfmo,3).le.1) write(iw,9000) ifg
c       if(maswrk) CLOSE(UNIT=IDAF,STATUS='DELETE')
        call CLOSDA('DELETE')
        CALL OPENDA(0)
        NEVALS=0
        iilay=ilay
        if(modfd.ne.0) iilay=layfrg(ifg)
        call makemol(ifg,0,0,iilay,0,0,0,0,0,0,0,.true.)
        l1=num
        l2=(l1*l1+l1)/2
        l3=l1*l1
        iodexch=scffrg(ifg).eq.rmc
c
c       Read densities.
c
        irec=irec0
        if(modfd.ne.0.and.layfrg(ifg).eq.1) irec=irec1
        call readmond(da,orbxch,iodexch,na,nb,l1,iodfmo,ifg+irec,
     *  scffrg(ifg).eq.uhf)
c
c       For FMO/FD, compute derivatives of <mu|dV/da|nu> for fragments
c       in layer 1 (mu,nu) and potentials V in layer 2 (atom "a"). 
c
        if(modfd.ne.0.and.iilay.eq.1.and.irststp.ge.2) then
c         DUM are to be "undummied" for relativistic wave functions,
c         see STVDER. 
c         Fake ifmostp/icurlay so that TVDER computes the required derivatives.
          ifmostps=ifmostp
          ifmostp=2
          icurlays=icurlay
          icurlay=1
          mplevls=mplevl
          mplevl=0
          tdsav=tddftyp
          tddftyp=rnone
c         write(6,*) 'wwwcalling  TVDER',ifg
          CALL TVDER(da,dum,dum,L2,0)
          tddftyp=tdsav
          mplevl=mplevls
          icurlay=icurlays
          ifmostp=ifmostps
        endif
c       Read external potentials, compute the Lagrangian.
        nrec=1+ifg+nfg
        CALL rareads(IDAFMO,iodfmo,fa,l2,nrec,0)
c       write(6,*) 'Zeroing out the Lagrangian',ifg
c       call vclr(fa,1,l2)
c       write(6,*) 'Lagrangian',ifg
c       call prtril(fa,l1)
        urohf=scffrg(ifg).eq.uhf.or.scffrg(ifg).eq.rohf
        if(urohf) then
c         write(6,*) 'wwwExecuting new code'
c         -2 * 1/2 * (Da * V * Da + Db * V * Db), Da,Db is alpha/beta density 
c         Get alpha density 
          call DMTX2(da,da(l2+1),na,l1,l1,0)
          call CPYTSQ(da,wrk(1+l2),l1,1)
          call TFTRI(wrk,fa,wrk(1+l2),WRK(1+l2+l3),l1,l1,l1)
c         Get beta density 
          if(scffrg(ifg).eq.rohf) then
            call DMTX2(da,da(l2+1),nb,l1,l1,0)
          else if(scffrg(ifg).eq.uhf) then
            call DMTX2(da,da(l2+1+l3),nb,l1,l1,0)
          end if 
          call CPYTSQ(da,wrk(1+l2),l1,1)
          call TFTRI(fb,fa,wrk(1+l2),WRK(1+l2+l3),l1,l1,l1)
          call daxpy(l2,one,wrk,1,fb,1)
          call dscal(l2,-one,fb,1)
        else
c         -2 * 1/4 * D * V * D, D is total density 
          call CPYTSQ(da,da(1+l2),l1,1)
          call TFTRI(fb,fa,da(1+l2),WRK,l1,l1,l1)
          call dscal(l2,-half,fb,1)
        endif
cnb     It seems probable that FMO-MCSCF gradient is incorrect with RESDIM>0.
        IF (IAND(MODGRD,42).EQ.42.and.layfrg(ifg).ge.ilay) THEN
          if(scffrg(ifg).ne.uhf) lnum=l3
          if(scffrg(ifg).eq.uhf) lnum=l3+l3
          CALL RAREADS(IDAFMO,IODFMO,DA(L2+1),lnum,IREC+IFG,0)
          call getcnt2ubar(1,L1,L3,IFG,IPTLG,NQMTFG,DA(L2+1),
     *                     FA,ZVWRK,ZVLAG,scffrg(ifg).eq.uhf,
     *                     scffrg(ifg).eq.rohf)
        END IF
C
c       IF(OUT) THEN
c         WRITE(IW,*) 'THE LAGRANGIAN MATRIX IS'
c         CALL PRTRI(fb,L1)
c       END IF
c         WRITE(IW,*) 'THE density for SDER is',nfg
c         CALL PRTRI(da,L1)
        call vclr(de,1,NAT*3)
        call SDER(fb,fa,L1,L2,.false.)
c       Not only masters accumulate the gradient here: all nodes have
c       contributions (no global sum in SDER). 
c       WRITE(IW,*) 'FMO ES gradient: Lagrangian terms',ifg
c       CALL EGOUT(DE,NAT)
c       call fmodeg(1,fmode(1,1,2),fmopg,iaglob)
        call fmodeg(1,fmode(1,1,nbody),fmopg,iaglob)
c       write(6,*) 'wwwstopped'
c       Lump ES dimer contributions from trimers to FMO3.
        if(iand(nprfmo,3).le.1) call timit(1)
  200 continue 
cpbc
c     restore the original atomic coordinates
c---  if(dopbcmd) call genfrgmic(0,1)
      if(dopbcmd) call genimagefg(1,0)
c     if(dopbcmd) then
c       call makeindexiu(ilay,loadhf,loadm,x(lindat),x(lfmozan),
c    *                   x(lfmoc),x(liabdfg),x(ljabdfg),x(lscffrg),
c    *                   1,x(lrij),rwrk,mfg2,vdwrad,
c    *                   nprfrg,.false.,nimgcell,natfmo+nbdfg,
c    *                   x(luntxyz),units,x(lindxiu),.false.)
c     end if
c
      if(isgddi) call GDDICOUNT( 1,lgroup,myjob)
c
c     WRITE(IW,9088)
c     CALL EGOUT(fmode(1,1,2),natfmo)
c9086 FORMAT(/10X,36(1H-)/10X,'gradient of NR+2e+1e (VDER+HELFEY)',
c    *       /10X,36(1H-))
c9088 FORMAT(/10X,36(1H-)/10X,'GRADIENT OF THE ENERGY',/10X,36(1H-))
c     write(6,*) 'Exit esdgrd'
      return
 9000 format(/1x,'Finishing separated dimer gradients, ifg=',I5)
      end
C*MODULE fmogrd  *DECK ESDnrder
      SUBROUTINE ESDnrder(nat0)
      use mx_limits, only: mxatm
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      logical GOPARR,DSKWRK,MASWRK
      COMMON /GRAD  / DE(3,MXATM)
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
C
c     Add contributions to the gradient due to nuclear energy for
c     separated dimers IJ. 
c     This subroutine produces wrong results for overlapping I and J.
c     I and J can only overlap if RESDIM<0, which should never happen!
c     This subroutine is statically parallelised.
c
c     call vclr(de,1,NAT*3)
      DO 100 K=1,NAT0
        if(mod(k,nproc).ne.me) goto 100
        ZAK=ZAN(K)
        DO L=nat0+1,nat
          ZAL=ZAN(L)
          PKL=(C(1,K)-C(1,L))**2+(C(2,K)-C(2,L))**2+(C(3,K)-C(3,L))**2
          SKL=ZAK*ZAL/(PKL*sqrt(PKL))
          DO KK = 1,3
            DERKL=(C(KK,K)-C(KK,L))*SKL
            DE(KK,K)=DE(KK,K)-DERKL
            DE(KK,L)=DE(KK,L)+DERKL
         enddo
       enddo
  100 continue
c     No global sum here: it is done in ESD2DER that is called after this. 
c     WRITE(IW,9088)
c     WRITE(6,*) 'FMO ES gradient: N-N terms'
c     CALL EGOUT(DE,NAT)
      RETURN
c9088 FORMAT(/10X,36(1H-)/10X,'GRADIENT OF NUCLEAR REPULSION FORCE',
c    *       /10X,36(1H-))
      END
C*MODULE fmogrd  *DECK esdderz
      SUBROUTINE esdderz(nfg2,iat0,layfrg,indat,iabdfg,jabdfg,indbd,
     *                   fmozan,fmoc,fzcor,iactfg,ndtfrg,znuc,cx,cy,cz,
     *                   untxyz,natfmob)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      logical layman,mixdim,dofdb
      COMMON /MCPFMO/ IMCPFMO,LFZCOR,LIFMPTYP,LIFMPTYP2,LMCPSW,
     *                LIZCOR2,imp0,jmp0,icorsh0,igtf0,IECPFMO
      common /fmodim/ maxbnd,maxknd,maxcbs,maxcao,maxvec,maxl1,maxnat,
     *                maxabd,maxbas,maxbbd,maxlmo,maxslo,maxabd2,maxrij
      Common /fmoinf/ nfg,nlayer,natfmo,nbdfg,naotyp,nbody,nsegm
      common /fmoopt/ espsca(9),RESPAP(2),rESPPC(2),rESDIM,restri(4),
     *                rcorsd,respct,convfg,cnvdmp,coroff,rflmo(4),
     *                orshft,orshft2,cnvafo,ascreen(4),IXESP,mxitfg,
     *                nguess,NBSSE,modorb,modpar,irststp,irstlay,nprfmo,
     *                nfmopal,modprp,maxl1c,ipieda,modgrd,modesp,ivmul,
     *                modlmo,nopden,mofock,modfd,modfmm,ncentm,ndualb,
     *                ngab,modpan
      COMMON /FMOPMD/ FMOBOX(3),MDWPBC,NIMGCELL,IMGLVL,
     *                LTRVEC,LFMOGCTR,LFMOCTMP,LINDATMD,LWRKDSAV,LINDXIU
     *               ,IPBCFST
      common /fmorun/ espscf,e0scf(2),emp2s,IDAFMO,icurfg,jcurfg,kcurfg,
     *                icurlay,icurunt,nat1e,ncursh,ngau,icurpop,ifmostp,
     *                moncor,needr,modrst,norbproj,nunesp,iskipesp,
     *                IESDPPC,idoprop,mp2run,icurit,idmfmo,iddfmo,
     *                iddcur,nddleft,ivmfmo,nzmtfmo,ifmobas,itmfmo(2)
C
      dimension layfrg(*),indat(*),iabdfg(*),jabdfg(*),indbd(maxabd,*),
     *          fmozan(*),fmoc(3,*),fzcor(*),iactfg(*),ndtfrg(nfg2,2),
     *          t(3),untxyz(3,natfmob,0:*)
c
c     Determine partial charges for the external potentials due to 
c     the terms in the sum(J) of the separated dimer 
c     IJ energy derivative, that is, derivatives of the form <mu|Vext|nu>. 
c     This subroutine is similar to fmoatfrg. 
c     Vext are just Coulomb point charges due to atoms in separated dimers
c     (relative to the atoms in the current fragment I, mu,nu are in I as well).
c     It produces wrong results for frac2 in case of overlapping I and J.
c     I and J can only overlap if RESDIM<0, which should never happen!
c
c     An atom IAT is divided into M+1 pieces. The first (major) piece is
c     the fragment to where it belongs (indat(iat)). Then IAT can be shared
c     with M other fragments (M.ge.0) (loop over iabd). Each shared piece
c     (along with the main) can either be included in or excluded from ZNUC,
c     as determined by interfragment distances (included for separated dimers). 
c
c     The cA comments change FMO/FD to treat only those dimers in layer 2
c     where at least one monomer is active (edimer should be changed too).
c
      dofdb=iand(modfd,4).ne.0
      if(MDWPBC.ne.0.and.iat0.gt.natfmo) then
        iat=mod(iat0-1,natfmo)+1
        iunit=(iat0-1)/natfmo
        cx=untxyz(1,iat,iunit)
        cy=untxyz(2,iat,iunit)
        cz=untxyz(3,iat,iunit)
      else
        iat=iat0
        iunit=0
        cx=fmoc(1,iat)
        cy=fmoc(2,iat)
        cz=fmoc(3,iat)
      endif
      IF (IMCPFMO.EQ.1) THEN
        znuc=fmozan(iat)-fzcor(iat)
      else if (iecpfmo.eq.1) then
        znuc=fmozan(iat)-fzcor(iat)
      else
        znuc=fmozan(iat)
      endif
      znuc2=0
      ifg=icurfg
      jfg=indat(iat)
      ilay=layfrg(ifg)
      jlay=layfrg(jfg)
      layman=(ilay.eq.icurlay .or.jlay.eq.icurlay).and.
     *        ilay.ge.icurlay.and.jlay.ge.icurlay
c     accept charges only from appropriate layers.
      if(modfd.ne.0.and.icurlay.eq.1) layman=.false.
      if(iand(modfd,2).ne.0) then
        if(icurlay.eq.2) layman=iactfg(ifg).eq.1.or.iactfg(jfg).eq.1
      endif
c     stop all layer 1 ES dimers except the mixed ones set by mixdim. 
      mixdim=modfd.ne.0.and.
     *       (icurlay.gt.1.and.iactfg(ifg).eq.1.and.layfrg(jfg).eq.1.or.
     *        icurlay.eq.1.and.layfrg(ifg).eq.1.and.iactfg(jfg).eq.1)
c     if(jfg.ne.ifg) then
      if(jfg.ne.ifg.or.(MDWPBC.ne.0.and.iunit.gt.0)) then
c       note that distance R(ifg,ifg) is not defined when iunit=0!
        icurunts=icurunt
        icurunt=iunit
c       note that distance R(ifg,ifg) is not defined!
        mmdim=0
      if(modfmm.ne.0) call mmdist(ifg,0,0,jfg,t,radius,ty2z,ratio,mmdim)
        ida=1
        if(nbody.gt.2) then
          iifg=max(ifg,jfg)
          jjfg=min(ifg,jfg)
          ijfg=(iifg*iifg-3*iifg)/2+jjfg+1
          ida=ida-ndtfrg(ijfg,1)
        endif
        znuc2a=znuc2
        if(.not.dofdb.or..not.mixdim) then
        if((fmodist(ifg,0,0,jfg).gt.resdim.and.layman .or.mixdim)
     *     .and.mmdim.eq.0) then
          znuc2=znuc*ida
c         full Z 
          znuc2a=znuc
        endif
        endif
c       write(6,*) 'wwwaaac1ea',jfg,ida,iat,znuc2
c       if(nbdfg.ne.0) then
        icurunt=icurunts
c       assume no bond detachment between fragments that belong to
c       different cells.
        if(nbdfg.ne.0.and.(MDWPBC.eq.0.or.iunit.eq.0)) then
          do iabd=1,maxabd
            ibdfg=indbd(iabd,iat)
            if(ibdfg.eq.0) goto 100
            ia=abs(iabdfg(ibdfg))
            ja=abs(jabdfg(ibdfg))
            jfg=indat(ja)
            jlay=layfrg(jfg)
            layman=(ilay.eq.icurlay .or.jlay.eq.icurlay).and.
     *              ilay.ge.icurlay.and.jlay.ge.icurlay
            if(modfd.ne.0.and.icurlay.eq.1) layman=.false.
c           layman is set to true or false so that it is ignored below
c           (.or. .not..true.); instead, mixdim is used solely.
            if(iand(modfd,2).ne.0) then
            if(icurlay.eq.2) layman=iactfg(ifg).eq.1.or.iactfg(jfg).eq.1
            endif
            mixdim=modfd.ne.0.and.
     *       (icurlay.gt.1.and.iactfg(ifg).eq.1.and.layfrg(jfg).eq.1.or.
     *        icurlay.eq.1.and.layfrg(ifg).eq.1.and.iactfg(jfg).eq.1)
            if(ia.eq.iat.and.ifg.ne.jfg) then
              mmdim=0
              if(modfmm.ne.0) call mmdist(ifg,0,0,jfg,t,radius,ty2z,
     *                                    ratio,mmdim)
              ida1=1
              if(nbody.gt.2) then
                iifg=max(ifg,jfg)
                jjfg=min(ifg,jfg)
                ijfg=(iifg*iifg-3*iifg)/2+jjfg+1
                ida1=ida1-ndtfrg(ijfg,1)
                if(znuc2a.ne.0) znuc2=znuc2-ida
c               write(6,*) 'wwwaaac1eb',ifg,jfg,ida,znuc2
c               subtract BAA proton's contribution (ida from above, not ida1) 
              endif
              if(((fmodist(ifg,0,0,jfg).le.resdim.or..not.layman).and.
     *           .not.mixdim).or.mmdim.ne.0) then
                if(znuc2a.ne.0.and.nbody.eq.2) znuc2=znuc2-ida
c               subtract BAA proton's contribution. 
              else
                if(.not.dofdb.or..not.mixdim) then
                if(znuc2a.eq.0.or.nbody.ne.2) znuc2=znuc2+ida1
c               add BAA proton's contribution. 
                endif
              endif
c             The special FMO3 code should also work for FMO2.
c             write(6,*) 'wwwaaac1ec',ifg,jfg,ida1,znuc2
            endif
          enddo
  100     continue
        endif
      endif
      znuc=-znuc2
      RETURN
      END
C*MODULE fmogrd  *DECK TIN2fmo
      SUBROUTINE tin2fmo(natfmo,fmozan,fmoc,ibas)
      use mx_limits, only: MAXATM => tinker_maxatm
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      dimension fmozan(*),fmoc(3,*),ibas(*)
C
C    GAMESS
      PARAMETER (MAXLNK=100, MAXR1=2000)
      PARAMETER (ONE=1.0D+00,units=0.52917724924D+00)
C
      CHARACTER*10 NAME,xyzname
      INTEGER TYPE,TAG,CLASS,ATOMIC,VALENCE
      REAL*8 DY,DZ,X,Y,Z,MASS
      CHARACTER*20 STORY
c     CHARACTER*8 prgnam
      LOGICAL IMOMM,SIMOMM,CHGEMB
      logical GOPARR,DSKWRK,MASWRK
      logical mmonly, qmmm
      DIMENSION ITMPLNK(3*MAXLNK)
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      common /fmodim/ maxbnd,maxknd,maxcbs,maxcao,maxvec,maxl1,maxnat,
     *                maxabd,maxbas,maxbbd,maxlmo,maxslo,maxabd2,maxrij
C
c     COMMON /COORDN/ UZNUC(MXATM),CUNIQ(MXATM,3),NATREL,UATOM(MXATM)
      COMMON /ATMTYP/ MASS(MAXATM),TAG(MAXATM),CLASS(MAXATM),
     *                ATOMIC(MAXATM),VALENCE(MAXATM),NAME(MAXATM),
     *                STORY(MAXATM), xyzname(maxatm)
      COMMON /ATOMST/ X(MAXATM),Y(MAXATM),Z(MAXATM),N,TYPE(MAXATM)
C
      COMMON /QMMM1/ IMOMM,SIMOMM,NPAIR,NSEQ,CHGEMB
      COMMON /QMMM2/ IQMATM(MAXR1),ibasfmo(MAXR1)
      COMMON /QMMM3/ LINKGE(3*MAXLNK),BLQM(MAXLNK),BLMM(MAXLNK)
      COMMON /TINOPT/ mparti,MMONLY, QMMM
C
c     FMO interface, cloned from TIN2GMS by DGF.
c     MMONLY
      if(mmonly) then
        DO I=1,natfmo
          fmozan(I)=ATOMIC(I)
          fmoc(1,I)=X(I)/units
          fmoc(2,I)=Y(I)/units
          fmoc(3,I)=Z(I)/units
        ENDDO
      endif
c     QMMM
      if(qmmm) then
        DO I=1,NSEQ
          fmozan(I)=ATOMIC(IQMATM(I))
          fmoc(1,I)=X(IQMATM(I))/units
          fmoc(2,I)=Y(IQMATM(I))/units
          fmoc(3,I)=Z(IQMATM(I))/units
c         UATOM(I)=NAME(IQMATM(I))
        ENDDO
C
C       TOTAL NUMBER OF QUANTUM ATOMS
        natfmo=nseq+npair
c
C       DANGLING BONDS ARE REPLACED BY HYDROGENS
        if(npair.gt.0) then
          DO I=1,NPAIR
            LR1=LINKGE(3*I-2)
            LR2=LINKGE(3*I  )
            DX=(X(LR2)-X(LR1))
            DY=(Y(LR2)-Y(LR1))
            DZ=(Z(LR2)-Z(LR1))
            TDIST=SQRT(DX**2+DY**2+DZ**2)
            DO J=1,NSEQ
               IF (LR1.EQ.IQMATM(J)) ITMPLNK(I)=J
            ENDDO
            ZAN1=ATOMIC(LR1)
            CALL GTDIST(ZAN1,ONE,RDIST)
            fmozan(I+NSEQ)=ONE
c           UATOM(I+NSEQ)='H         '
            fmoc(1,I+NSEQ)=(X(LR1)+DX*RDIST/TDIST)/units
            fmoc(2,I+NSEQ)=(Y(LR1)+DY*RDIST/TDIST)/units
            fmoc(3,I+NSEQ)=(Z(LR1)+DZ*RDIST/TDIST)/units
          ENDDO
c         UATOM(NSEQ+NPAIR+1)='$END      '
          CALL SEQREW(IR)
          call fndgrp(IR,' $LNKXYZ',IEOF)
          if(maswrk) then
           if(ieof.eq.0) then
              call lnkxyz(' $LNKXYZ',nseq,npair,fmoc)
           endif 
          endif
          if(GOPARR) then
               CALL DDI_BCAST(2401,'F',fmoC,3*natfmo,MASTER)
          end if
          if(maswrk) then
            write(iw,990)
            do i=1,npair
              write(iw,980) i+nseq,fmoc(1,i+nseq)*units,
     *                 fmoc(2,i+nseq)*units, fmoc(3,i+nseq)*units
            enddo
          endif
C
C         ASSIGN LINKGE VALUES
C
          DO I=1,NPAIR
            LINKGE(I*3-2)=ITMPLNK(I)
          ENDDO
        endif
        mxibas=1
        do i=1,nseq+npair
          ibas(i)=ibasfmo(i)
          if(ibas(i).gt.mxibas) mxibas=ibas(i)
        enddo
      endif
      maxbas=mxibas
C
      RETURN
  990 format(/,10x,'Link atom coordinates (angstrom)',/)
  980 format(10x,i5,3f20.8)
      END
C*MODULE fmogrd  *DECK gmsxyz
C>
C>     @brief GAMESS <-> FMO coordinate exchange 
C>
C>     @details Interface for GAMESS <-> FMO coordinate exchange. 
C>
C>     @author Dmitri Fedorov
C>
      subroutine gmsxyz(modef,nfthes,nc1,natfmo,layfrg,indat,fmozan,
     *                  fmoc,fmomas,izbas,iactfg,nfcrec,nevsav,grandmas,
     *                  nactfg,fmosym,nfmosym,quit,mdoutmin)
      use mx_limits, only: mxatm,mxao
      USE comm_EFPFMO
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      logical quit,GOPARR,DSKWRK,MASWRK,LINEAR,gopsav,grandmas,MCMIN,
     *        ISGDDI,PAROUT,INITGDDI,wasgddi,doremd,MLGDDI,mdoutmin,
     *        THRDCAL
      integer ddi_world,ddi_group
      Parameter(ddi_world=0,ddi_group=1)
      character*4 atm
      character*8 cdum,nill
      character*80 fmosym
      COMMON /ECP2  / CLP(400),ZLP(400),NLP(400),KFIRST(MXATM,6),
     *                KLAST(MXATM,6),MMAX(MXATM),LPSKIP(MXATM),
     *                IZCORE(MXATM)
      COMMON /FMCOM / X(1)
      COMMON /FROZNC/ IFREEZ(3*MXATM),NFRZ
      COMMON /GDDI/   ISCOPE,NGROUPS,MYGROUP,MEGLOB,npglob,nnglob,JBTYP,
     *                ISGDDI,PAROUT,INITGDDI,wasgddi,MLGDDI,NSUBGR,
     *                MeUniv,NPUniv,numdlb,myworld,nworlds,mogddi
      COMMON /GLOFMO/ MCMIN,IFMOID(MXATM)
      COMMON /INFO  / CDUP(MXATM,3),IZAN(MXATM),NATOMS,IUNTRD,ATM(MXATM)
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON /MACHIN/ NWDVAR,MAXFM,MAXSM,LIMFM,LIMSM
      COMMON /MACHSW/ KDIAG,ICORFL,IXDR,modio,mem10,lpnt10,mem10m
      COMMON /MCPFMO/ IMCPFMO,LFZCOR,LIFMPTYP,LIFMPTYP2,LMCPSW,
     *                LIZCOR2,imp0,jmp0,icorsh0,igtf0,IECPFMO
      COMMON /OPTGRD/ XYZ(3*MXATM),ES,FE(20),
     *                CONVF,FMAXT,DXMAXT,RMAX,RMIN,RLIM,
     *                EIGMAX,EIGMIN,GRDERR,FRMS,FMAX,TRMAX,TRMIN,
     *                IC(20),NSTEP,NSERCH,NPMAX,NP,IFOLOW,
     *                NNEG,IUPHSS,IEXIT,ITRUPD,IPAD,kdiagf,NPRICO
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /RUNLAB/ TITLE(10),ANAM(MXATM),BNAM(MXATM),BFLAB(MXAO)
      COMMON /RUNOPT/ RUNTYP,EXETYP,NEVALS,NGLEVL,NHLEVL
c     COMMON /THRDER/ THRDCAL
      COMMON /ZMAT  / NZMAT,NZVAR,NVAR,NSYMC,LINEAR
      COMMON /MASSES/ ZMASS(MXATM)
      common /fmoopt/ espsca(9),RESPAP(2),rESPPC(2),rESDIM,restri(4),
     *                rcorsd,respct,convfg,cnvdmp,coroff,rflmo(4),
     *                orshft,orshft2,cnvafo,ascreen(4),IXESP,mxitfg,
     *                nguess,NBSSE,modorb,modpar,irststp,irstlay,nprfmo,
     *                nfmopal,modprp,maxl1c,ipieda,modgrd,modesp,ivmul,
     *                modlmo,nopden,mofock,modfd,modfmm,ncentm,ndualb,
     *                ngab,modpan
      common /fmorun/ espscf,e0scf(2),emp2s,IDAFMO,icurfg,jcurfg,kcurfg,
     *                icurlay,icurunt,nat1e,ncursh,ngau,icurpop,ifmostp,
     *                moncor,needr,modrst,norbproj,nunesp,iskipesp,
     *                IESDPPC,idoprop,mp2run,icurit,idmfmo,iddfmo,
     *                iddcur,nddleft,ivmfmo,nzmtfmo,ifmobas,itmfmo(2)
      dimension layfrg(*),indat(*),fmozan(*),fmoc(3,*),fmomas(*),
     *          izbas(*),iactfg(*)
      equivalence (dum,cdum),(dnill,nill)
      data OPTMIZ /8HOPTIMIZE/, SADPT  /8HSADPOINT/,GLOBOP /8HGLOBOP  /,
     *     dnill/8H        /,cdum/' '/,HSSIAN/8HHESSIAN /
      DATA FMOMD  /8HMD      /,RAMAN  /8HRAMAN   /,FMOHESS/8HFMOHESS /
      DATA AIRC   /8HIRC     /,AMEX/8HMEX     /,ANEBP/8HNEBPATH /
c
c     Interface GAMESS and FMO coordinates if needed
c     quit is set to true for SIGX runs, that at this point only read 
c     coordinates.
c
c     write(iw,1234) runefp
c1234 format(1x,'wwwhav ',2A8)
c     if(runefp.eq.AMEX.and.modef.ne.0.and.
c    *   nsubgr.gt.0.and.myworld.ne.0) then
c       if(maswrk) write(iw,*) 'Skipping dictionary file bookkeeping.'
c       Other worlds do not need F10 for MEX runs, and they have no data
c       stored too.
c       return
c     endif
      if(mdoutmin) call timit(1)
      quit=.false. 
c     if(qmmm) then
c       if(modef.ne.0) call tin2fmo(natfmo,fmozan,fmoc)
c       return
c     endif   
c
      natfmof=(natfmo-1)/nwdvar+1
c
      modeAB=0
      if(modfd.ne.0.and.runtyp.eq.HSSIAN) modeAB=1
c
      THRDCAL=.false.
      if(modef.eq.0) then
c       natbl=nat
        if(runtyp.eq.OPTMIZ.or.runtyp.eq.SADPT.or.runtyp.eq.HSSIAN.or.
     *    runefp.eq.globop.OR.RUNEFP.EQ.FMOMD.or.runtyp.eq.RAMAN.or.
     *    (runtyp.eq.FMOHESS.and.modfd.ne.0).or.runtyp.eq.airc.or.
     *    runtyp.eq.AMEX.or.(runtyp.eq.FMOHESS.and.THRDCAL).or.
     *    runefp.eq.ANEBP) then
c
c         copy FMO -> GAMESS
          call checkremd(RUNEFP,doremd)
          if(isgddi.and..not.doremd) call gddi_scope(ddi_world)
          if(modfd.eq.0) then
            nat=natfmo
            if(runefp.ne.fmomd) then
              if(nat.gt.MXATM) then
                if(maswrk) write(iw,9000) MXATM
                call abrt
              endif 
              call dcopy(3*nat,fmoc,1,C,1)
              call dcopy(nat,fmozan,1,ZAN,1)
            endif
c         FMO/MD does not need to copy data, which are taken from FMOPNT arrays.
          else
            call packat(modeAB,layfrg,indat,fmozan,fmoc,fmomas,iactfg)
c           packat sets NAT, C and ZAN.
          endif
          CALL DAWRIT(IDAF,IODA,izbas,natfmof,495,1)
c         we write to F10 to all masters here, just to enable FMO-REMD...
          if(isgddi.and.doremd) call gddi_scope(ddi_world)
c         if(isgddi) call gddi_scope(ddi_world)
c         if(isgddi) call gddi_scope(ddi_group)
c
          if(runefp.ne.fmomd.or.modfd.ne.0) then
          IF (IMCPFMO.EQ.1) call viclr(iZCORe,1,NATFMO)
c         For MCP, the nuclear charge is not split between ZAN and IZCORE.
c         Instead, the sum is stored in ZAN and IZCORE is zero.
c         It is easy to split ZAN and IZCORE, but is that necessary?
          do i=1,nat
            ian(i)=int(ZAN(i)+0.5D-02)
c           IF (IMCPFMO.EQ.1) IZCORE(I)=INT(FZCOR(I)+0.01D+00)
            call zsymnum(cdum,ATM(i),ian(i))
c           write(UNIT=cdum,FMT='(A4,I4)') ATM(i),i
            write(UNIT=cdum,FMT='(A4,A4)') ATM(i),'    '
            ANAM(i)=dum
            BNAM(i)=dnill
          enddo 
c         Generate or read Z-matrix.
          IF(NZVAR.GT.0.and.nzmtfmo.eq.0) THEN
            CALL ZMATIN(NZVAR.EQ.0)
            if(natfmo.gt.MXATM) call abrtx("Too many atoms in gmsxyz")
            call dcopy(natfmo,fmomas,1,zmass,1)
c           BANDBI must have zmass set up...
            CALL BANDBI
            if(nevsav.eq.0) CALL PZANDG(egdum,0)
          END IF
          endif
          if(runefp.eq.ANEBP) call dcopy(natfmo,fmomas,1,zmass,1)
c         if(isgddi) call gddi_scope(ddi_group)
c
c         Now comes the mean part. We rewrite BNAM for our evil purposes.
c         (namely, storing basis set ID for each atom).
c
c         do i=1,nat
c           dum=BNAM(i)
c           write(UNIT=cdum,FMT='(I8)') izbas(i)
c           BNAM(i)=dum
c         enddo
C
C         Tricky modifiction to prevent writing record 4
          IF(RUNTYP.EQ.AIRC.or.THRDCAL) THEN
            NSERCH=1
          END IF
          if(maswrk.and.mdoutmin) write(6,*) 'Copied FMO -> GAMESS'
          quit=.true. 
        endif
      else
c       running with modef=1 means being called from SIGX or MDX.
c       save the Hessian
c       copy GAMESS -> FMO
        nactat=nat-NFRZ/3
        if(isgddi) call gddi_scope(ddi_world)
c       CALL DDI_BCAST(2422,'F',c,3*nat,master)
c       write(6,*) 'wwwbcastc'
c       make sure the coordinates are exactly the same
        if(modfd.eq.0) then 
          natfmo=nat
          if(runefp.ne.fmomd) then
            call dcopy(nat,zan,1,fmoZAN,1)
            call dcopy(3*nat,c,1,fmoC,1)
          endif
        else
          call unpackat(modeAB,layfrg,indat,fmozan,fmoc,iactfg)
          nc1=nat*3
c         overwrite nc1...
c     if(nevsav.eq.1) fmoC(2,37)=fmoC(2,37)-0.005D+00/0.52917724924D+00
c     if(nevsav.eq.2) fmoC(2,37)=fmoC(2,37)+0.005D+00/0.52917724924D+00
c     nat=natfmo
c     call dcopy(3*nat,fmoC,1,C,1)
c     CALL PRATM(IW,1)
        endif
c       IF (RUNEFP.EQ.hssian.AND.modfd.eq.0) then
c         izbas is now reset to 1!!!
c         if(maswrk) write(6,*) 'izbas is reset!',natfmof
c         call icopy(natfmof,1,0,izbas,1)
c         if(maswrk) write(6,*) 'izbas reset!'
c       if(isgddi) call gddi_scope(ddi_group)
c         GOTO 555
c       endif
c         if(isgddi) call gddi_scope(ddi_world)
        CALL DAread(IDAF,IODA,izbas,natfmof,495,1)
c         if(isgddi) call gddi_scope(ddi_world)
       if(modfd.ne.0.and.grandmas.and.(RUNEFP.NE.RAMAN.or.nevsav.eq.0))
     *   call outxyz(nevsav,title,natfmo,fmozan,fmoc,izbas,nactat,
     *               nactfg,fmosym,nfmosym)
        if(isgddi) call gddi_scope(ddi_group)
c       Hessian was saved in the world scope which means only the grand
c       master has it. Moreover, the whole group must call daread.
c
c       For SCHLEGEL optimisations Hessian storage will be used to save
c       Coords+grads of NPMAX previous points.
c
C 
C            MCMIN: ENERGY+GRAD CALCULATION
C       .NOT.MCMIN: ONLY ENERGY CALCULATION
C 
        IF ((RUNEFP.EQ.GLOBOP.AND..NOT.MCMIN).or.RUNEFP.EQ.FMOMD)
     *    GOTO 555
        IF (RUNTYP.EQ.hssian.AND.modfd.eq.0) GOTO 555
        IF(runefp.eq.ANEBP) GOTO 555
c       IF(runefp.eq.AMEX) GOTO 555
c       MEX runs do not use Hessian on DAF.
c       IF(runefp.eq.AMEX.and.nsubgr.gt.0.and.myworld.ne.0) GOTO 555
c       Worlds other than 0 do not need the Hessian for AMEX.
        gopsav=goparr
        if(mem10.eq.0) goparr=.false.
        if((grandmas.or.mem10.ne.0).and.
     *     (modfd.eq.0.or.runtyp.ne.fmohess)) then
c    *    (modfd.eq.0.or.runtyp.ne.fmohess).and.runtyp.ne.hssian) then
          nc3=nc1*nc1
          nc13=3*nc1
          if(runtyp.ne.hssian) nc13=0
          nc1sch=nc1*NPMAX
          CALL VALFM(LOADFM)
          lfcm=LOADFM+1
          lddm=lfcm+max(nc3,nc1sch)
          last=lddm+nc13 
          IF(NZVAR.GT.0) THEN
            NTEMP     = NZMAT
            NZMAT     = nzmtfmo
            NCOORD    = natfmo * 3
            LFMOZMTIN = LAST
C  FOR  RECORD 39
            NDIM      = NZVAR
C  FOR  RECORD 40 
            NDIM      = NDIM  + NZMAT
C  FOR  RECORD 41 
            NDIM      = NDIM  + 3 * NAT * NZVAR
C  FOR  RECORD 46 
            NDIM      = NDIM  + NZVAR * NZVAR
C  FOR  RECORD 47 
            NDIM      = NDIM  + NVAR
C  FOR  RECORD 48 
            NDIM      = NDIM  + NZMAT
C  FOR  RECORD 49 
            NDIM      = NDIM  + NCOORD * NVAR
C  FOR  RECORD 50
            NDIM      = NDIM  + NCOORD * NVAR
            LAST      = LFMOZMTIN + NDIM
            NDIMSAV   = NDIM
          END IF
          NEEDH = LAST- LOADFM -1
          CALL GETFM(NEEDH)
c         Hessian was written in world scope divided in parallel
c         and must be read in the same scope here.
          if(isgddi.and.mem10.ne.0) call gddi_scope(ddi_world)
          CALL daread(IDAF,IODA,x(lFCM),NC3,nfcrec,0)
          if(isgddi.and.mem10.ne.0) call gddi_scope(ddi_group)
         if(runtyp.eq.hssian) CALL daread(IDAF,IODA,x(lddm),NC13,34,0)
c         Only the grand master writes to NFTHES.
          if(grandmas) then
          call seqrew(nfthes)
c         Hessian restart: if the Hessian file was there, use it rather
c         than the initial Hessian computed during this run. 
          igothes=0 
          if(nserch.eq.0.and.runtyp.ne.hssian) then 
            READ(nfthes,END=604) nserch0,nc10
c           aha! read was successful, should we use continuous nserch?
            if(nc1.ne.nc10) then
              write(6,*) 'Hessian found but the size is wrong',
     *                   nc1,nc10
              call abrt
            endif
            write(6,*) 'Hessian restart successful, ',nserch0
            igothes=1 
c           DO NOT overwrite nfthes
  604       continue
          endif 
          endif
c         Only for mem10 all cores get here.
          if(mem10.ne.0) CALL DDI_BCAST(2422,'I',igothes,1,master)
          if(igothes.eq.0) then
            IF(NZVAR.GT.0) THEN
              LPNT = LFMOZMTIN
              NDIM = NZVAR
              CALL DAREAD(IDAF,IODA,x(LPNT),NDIM,39,0)
              LPNT = LPNT + NDIM
              NDIM = NZMAT
              CALL DAREAD(IDAF,IODA,x(LPNT),NDIM,40,1)
              LPNT = LPNT + NDIM
              NDIM = 3 * NATFMO * NZVAR
              CALL DAREAD(IDAF,IODA,x(LPNT),NDIM,41,0)
              LPNT = LPNT + NDIM
              NDIM = NZVAR * NZVAR
              CALL DAREAD(IDAF,IODA,x(LPNT),NDIM,46,0)
              LPNT = LPNT + NDIM
              NDIM = NVAR
              CALL DAREAD(IDAF,IODA,x(LPNT),NDIM,47,0)
              LPNT = LPNT + NDIM
              NDIM = NZMAT
              CALL DAREAD(IDAF,IODA,x(LPNT),NDIM,48,1)
              LPNT = LPNT + NDIM
              NDIM = NCOORD * NVAR
              CALL DAREAD(IDAF,IODA,x(LPNT),NDIM,49,0)
              LPNT = LPNT + NDIM
              NDIM = NCOORD * NVAR
              CALL DAREAD(IDAF,IODA,x(LPNT),NDIM,50,0)
              LPNT = LPNT + NDIM
c             write(*,*) "finish rading ZMAT file"
              NZMAT = NTEMP     
            END IF
c           write(6,*) 'No Hessian found, using the initial one'
            if(grandmas) then
            call seqrew(nfthes)
            write(nfthes) nserch,nc1 
            write(nfthes) (x(lFCM+i),i=0,nc3-1)
            if(runtyp.eq.hssian) write(nfthes) (x(lddm+i),i=0,nc13-1)
            IF(NZVAR.GT.0) write(nfthes) (x(LFMOZMTIN+i),i=0,NDIMSAV-1)
            endif
          endif
          if(NPMAX.ne.0.and.nserch.gt.0) then
            write(6,*) 'Saving Schlegel...'
            CALL daread(IDAF,IODA,x(lFCM),nc1sch,98,0) 
            if(grandmas) write(nfthes) (x(lFCM+i),i=0,nc1sch-1)
            CALL daread(IDAF,IODA,x(lFCM),nc1sch,99,0) 
            if(grandmas) write(nfthes) (x(lFCM+i),i=0,nc1sch-1)
          endif
          CALL RETFM(NEEDH)
        else
c         call vclr(FCM,1,NC1*NC1)
        endif
        goparr=gopsav
        if(maswrk.and.mdoutmin) write(6,*) 'Copied GAMESS -> FMO',NC1
c       if(meglob.eq.0) call prsq(FCM,nc1,nc1,nc1)
  555   CONTINUE
c       do i=1,natfmo 
c         dum=BNAM(i)
c         read(UNIT=cdum,FMT='(I8)') izbas(i)
c       enddo
      endif
      if(RUNEFP.EQ.GLOBOP) then
c       This may be an overkill check, but just to be sure.
        if(nat.gt.mxatm) call abrtx("Too many atoms in gmsxyz 2")
        DO IIK=1,NAT
           IFMOID(IIK)=INDAT(IIK)
        END DO
      endif
      if(mdoutmin) call timit(1)
      return
 9000 format(/1x,'Too many atoms: ',I7, 
     *           '. Recompile with a larger MXATM.',/)
      END
C*MODULE fmogrd  *DECK gmsprop
C>
C>     @brief GAMESS <-> FMO property exchange 
C>
C>     @details Interface for GAMESS <-> FMO property exchange. 
C>
C>     @author Dmitri Fedorov
C>
      subroutine gmsprop(modef,nfthes,nc1,layfrg,indat,fmode,fmomas,
     *                  iactfg,nfcrec,grandmas0,LGLOE,dofddpcm,esolv)
      use mx_limits, only: mxatm,mxrt
      USE comm_EFPFMO
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      integer ddi_world,ddi_group
      logical grandmas,LGLOE,LINEAR,ISGDDI,PAROUT,INITGDDI,wasgddi,
     *        MLGDDI,GOPARR,DSKWRK,MASWRK,grandmas0
      logical dofddpcm
      dimension esolv(*)
      dimension layfrg(*),indat(*),fmode(3,*),fmomas(*),iactfg(*)
      Parameter(ddi_world=0,ddi_group=1)
      COMMON /ENRGYS/ ENUCR,EELCT,ETOT,SZ,SZZ,ECORE,ESCF,EERD,E1,E2,
     *                VEN,VEE,EPOT,EKIN,ESTATE(MXRT),STATN,EDFT(3),EDISP
      COMMON /FMCOM / X(1)
      COMMON /FUNCT / EHF,EG(3*MXATM)
      COMMON /GDDI/   ISCOPE,NGROUPS,MYGROUP,MEGLOB,npglob,nnglob,JBTYP,
     *                ISGDDI,PAROUT,INITGDDI,wasgddi,MLGDDI,NSUBGR,
     *                MeUniv,NPUniv,numdlb,myworld,nworlds,mogddi
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON /MACHSW/ KDIAG,ICORFL,IXDR,modio,mem10,lpnt10,mem10m
      COMMON /MASSES/ ZMASS(MXATM)
      COMMON /OPTGRD/ XDUM(3*MXATM),ES,FE(20),
     *                CONVF,FMAXT,DXMAXT,RMAX,RMIN,RLIM,
     *                EIGMAX,EIGMIN,GRDERR,FRMS,FMAX,TRMAX,TRMIN,
     *                IC(20),NSTEP,NSERCH,NPMAX,NP,IFOLOW,
     *                NNEG,IUPHSS,IEXIT,ITRUPD,IPAD,kdiagf,NPRICO
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /RUNOPT/ RUNTYP,EXETYP,NEVALS,NGLEVL,NHLEVL
      COMMON /ZMAT  / NZMAT,NZVAR,NVAR,NSYMC,LINEAR
      Common /fmoinf/ nfg,nlayer,natfmo,nbdfg,naotyp,nbody,nsegm
      common /fmoopt/ espsca(9),RESPAP(2),rESPPC(2),rESDIM,restri(4),
     *                rcorsd,respct,convfg,cnvdmp,coroff,rflmo(4),
     *                orshft,orshft2,cnvafo,ascreen(4),IXESP,mxitfg,
     *                nguess,NBSSE,modorb,modpar,irststp,irstlay,nprfmo,
     *                nfmopal,modprp,maxl1c,ipieda,modgrd,modesp,ivmul,
     *                modlmo,nopden,mofock,modfd,modfmm,ncentm,ndualb,
     *                ngab,modpan
      common /fmorun/ espscf,e0scf(2),emp2s,IDAFMO,icurfg,jcurfg,kcurfg,
     *                icurlay,icurunt,nat1e,ncursh,ngau,icurpop,ifmostp,
     *                moncor,needr,modrst,norbproj,nunesp,iskipesp,
     *                IESDPPC,idoprop,mp2run,icurit,idmfmo,iddfmo,
     *                iddcur,nddleft,ivmfmo,nzmtfmo,ifmobas,itmfmo(2)
      COMMON /SCYVEC/ IDFYVEC,LIPTYVEC,LIPTYA,LFCMWRK
      COMMON /SCZLAG/ LZVLAG,LZVWRK,LYALAG,LYAWRK,LFEQ1
      data HSSIAN/8HHESSIAN /,  FMOMD   /8HMD      /
      data FMOHESS/8HFMOHESS /,AMEX/8HMEX     /,ANEBP/8HNEBPATH /
c     data OPTMIZ /8HOPTIMIZE/, SADPT  /8HSADPOINT/
c 
c     Copy FMO properties to conventional GAMESS storage.
c  
      grandmas=grandmas0
      IF(runefp.eq.AMEX.and.nsubgr.gt.0.and.myworld.ne.0) 
     *   grandmas=.false.
      IF (LGLOE) THEN
        EHF = ETOT
        CALL CLOSDA('DELETE')
        CALL OPENDA(0)
        CALL DAWRIT(IDAF,IODA,ENUCR,MXRT+15,2,0)
        RETURN
      END IF
c     if(maswrk) CLOSE(UNIT=IDAF,STATUS='DELETE')
      call CLOSDA('DELETE')
c
      modeAB=0
      if(modfd.ne.0.and.runtyp.eq.HSSIAN) modeAB=1
c
c     Write this in the world scope for OPTIMIZE to read
c
      if(isgddi) call gddi_scope(ddi_world)
      CALL OPENDA(0)
      if(modef.eq.0) then
        if(dofddpcm) then
          lenrec=ixftch(x(liptya),1+nfg)
          CALL DAWRIT(IDAF,IODA,X(LYALAG),lenrec,570,0)
          CALL DAWRIT(IDAF,IODA,esolv,nfg,571,0)
        end if
      end if
      if(modef.ne.0) then
        EHF=etot
        if(modfd.eq.0) then
          if(runefp.ne.fmomd) then
            if(natfmo.gt.mxatm) call abrtx("Too many atoms in gmsprop")
            call dcopy(3*natfmo,fmode,1,eg,1)
            call dcopy(natfmo,fmomas,1,zmass,1)
c           write(6,*) 'gradient Copied',natfmo
          endif
        else
          if(dofddpcm) then
            lenrec=ixftch(x(liptya),1+nfg)
            CALL DAWRIT(IDAF,IODA,X(LYALAG),lenrec,570,0)
            CALL DAWRIT(IDAF,IODA,esolv,nfg,571,0)
          end if
c           write(6,*) 'gradient gmsprop',natfmo
c     call fmogout(1,.true.,indat,fmozan,fmode,DEFTF,TORQF)
          call packprp(modeAB,layfrg,indat,fmode,fmomas,iactfg,nat0)
          nc1=nat0*3
c         overwrite nc1...
        endif
        IF (RUNEFP.EQ.FMOMD) RETURN
        IF (RUNTYP.EQ.FMOHESS) RETURN
        IF (RUNTYP.EQ.hssian.AND.modfd.eq.0) RETURN
        IF (RUNEFP.EQ.ANEBP) RETURN
c
c       only the grand master has the true hessian, other nodes write zeroes.
c       It is important, because otherwise indexing of F10 will be broken.
        nc3=nc1*nc1
        nc13=3*nc1
        if(runtyp.ne.hssian) nc13=0
        nc1sch=nc1*NPMAX
        CALL VALFM(LOADFM)
        lfcm=LOADFM+1
        lddm=lfcm+max(nc3,nc1sch)
        last=lddm+nc13 
        IF(NZVAR.GT.0.and.nzmtfmo.ne.0) then
          NTEMP     = NZMAT
          NZMAT     = nzmtfmo
          NCOORD    = natfmo * 3
          LFMOZMTIN = LAST
C  FOR  RECORD 39
          NDIM      = NZVAR
C  FOR  RECORD 40 
          NDIM      = NDIM  + NZMAT
C  FOR  RECORD 41 
          NDIM      = NDIM  + 3 * NATFMO * NZVAR
C  FOR  RECORD 46 
          NDIM      = NDIM  + NZVAR * NZVAR
C  FOR  RECORD 47 
          NDIM      = NDIM  + NVAR
C  FOR  RECORD 48 
          NDIM      = NDIM  + NZMAT
C  FOR  RECORD 49 
          NDIM      = NDIM  + NCOORD * NVAR
C  FOR  RECORD 50
          NDIM      = NDIM  + NCOORD * NVAR
          LAST      = LFMOZMTIN + NDIM
          NDIMSAV   = NDIM
        END IF
        NEEDH = LAST- LOADFM -1
        CALL GETFM(NEEDH)
c       The way of reading below can possibly raise portability issues:
c       (for parallel runs, reading and writing by a single node).
c       if(meglob.eq.0) then
        if(grandmas) then
           call seqrew(nfthes)
           read(nfthes) nserch0,nc10
c          dummy
           if(nserch0.lt.0) nserch0=nserch+nc10
c          read(nfthes) (x(lFCM+i),i=0,nc3-1)
c          if(runtyp.eq.hssian) read(nfthes) (x(lddm+i),i=0,nc13-1)
           call readfcm(nfthes,nc3,nc13,x(lFCM),x(lddm))
           if(NZVAR.GT.0.and.nzmtfmo.ne.0) 
     *      read(nfthes) (x(LFMOZMTIN+i),i=0,NDIMSAV-1)
        else
           call vclr(x(lFCM),1,nc3) 
           if(nc13.ne.0) call vclr(x(lddm),1,nc13) 
        endif  
c       Tricky part: for in-core F10 need to broadcast the Hessian
c       because slaves save on F10 too; otherwise only the master does.
c       It is possible to save a bit by doing a send rather than broadcst
c       (but that would require low level processing).
        if(mem10.ne.0)
     *   CALL DDI_BCAST(2422,'F',x(lFCM),NC1*NC1,master)
        CALL DAWRIT(IDAF,IODA,x(lFCM),NC1*NC1,nfcrec,0)
        if(runtyp.eq.hssian)
     *    CALL DAWRIT(IDAF,IODA,x(lddm),NC13,34,0)
        if(NPMAX.ne.0.and.nserch.gt.0) then
          write(6,*) 'Getting Schlegel...'
          if(grandmas) read(nfthes) (x(lFCM+i),i=0,nc1sch-1)
          CALL dawrit(IDAF,IODA,x(lFCM),nc1sch,98,0)
          if(grandmas) read(nfthes) (x(lFCM+i),i=0,nc1sch-1)
          CALL dawrit(IDAF,IODA,x(lFCM),nc1sch,99,0)
        endif
c       
        if(NZVAR.GT.0.and.nzmtfmo.ne.0) THEN
            LPNT = LFMOZMTIN
            NDIM = NZVAR
            CALL DAWRIT(IDAF,IODA,x(LPNT),NDIM,39,0)
            LPNT = LPNT + NDIM
            NDIM = NZMAT
            CALL DAWRIT(IDAF,IODA,x(LPNT),NDIM,40,1)
            LPNT = LPNT + NDIM
            NDIM = 3 * NATFMO * NZVAR
            CALL DAWRIT(IDAF,IODA,x(LPNT),NDIM,41,0)
            LPNT = LPNT + NDIM
            NDIM = NZVAR * NZVAR
            CALL DAWRIT(IDAF,IODA,x(LPNT),NDIM,46,0)
            LPNT = LPNT + NDIM
            NDIM = NVAR
            CALL DAWRIT(IDAF,IODA,x(LPNT),NDIM,47,0)
            LPNT = LPNT + NDIM
            NDIM = NZMAT
            CALL DAWRIT(IDAF,IODA,x(LPNT),NDIM,48,1)
            LPNT = LPNT + NDIM
            NDIM = NCOORD * NVAR
            CALL DAWRIT(IDAF,IODA,x(LPNT),NDIM,49,0)
            LPNT = LPNT + NDIM
            NDIM = NCOORD * NVAR
            CALL DAWRIT(IDAF,IODA,x(LPNT),NDIM,50,0)
            NZMAT = NTEMP
        END IF
C
        CALL RETFM(NEEDH)
      endif
      CALL DAWRIT(IDAF,IODA,ENUCR,MXRT+15,2,0)
c     punch Cartesian hessian for unconverged runs (not ready for internals) 
c     if(nerr.gt.1.and.nfcrec.eq.4) then
c       write(6,*) 'Emergent punch of hessian enforced',NC1
c       CALL FCMPUN(x(lFCM),NC1)
c     endif
      if(isgddi) call gddi_scope(ddi_group)
      return
      END
c
C*MODULE fmogrd  *DECK fmorst
      subroutine fmorst(nfg2,nfg3,nen,nent,l1fmo,NDAR30a,NDAR30fa,iodfmo
     *                 ,layfrg,emon,edim,etrim,idmrec0,ifgfmo0,skipfini,
     *                  irestl,nts0,nts2,outpune,icanuse1,some)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      dimension iodfmo(*),layfrg(*),emon(*),edim(*),etrim(*),idmrec0(*)
      logical skipfini,outpune,some,GOPARR,DSKWRK,MASWRK,icanuse1
      COMMON /RAIOLN/ JRECLN(10),JRECST(10)
      Common /fmoinf/ nfg,nlayer,natfmo,nbdfg,naotyp,nbody,nsegm
      common /fmoopt/ espsca(9),RESPAP(2),rESPPC(2),rESDIM,restri(4),
     *                rcorsd,respct,convfg,cnvdmp,coroff,rflmo(4),
     *                orshft,orshft2,cnvafo,ascreen(4),IXESP,mxitfg,
     *                nguess,NBSSE,modorb,modpar,irststp,irstlay,nprfmo,
     *                nfmopal,modprp,maxl1c,ipieda,modgrd,modesp,ivmul,
     *                modlmo,nopden,mofock,modfd,modfmm,ncentm,ndualb,
     *                ngab,modpan
      common /fmorun/ espscf,e0scf(2),emp2s,IDAFMO,icurfg,jcurfg,kcurfg,
     *                icurlay,icurunt,nat1e,ncursh,ngau,icurpop,ifmostp,
     *                moncor,needr,modrst,norbproj,nunesp,iskipesp,
     *                IESDPPC,idoprop,mp2run,icurit,idmfmo,iddfmo,
     *                iddcur,nddleft,ivmfmo,nzmtfmo,ifmobas,itmfmo(2)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /RUNOPT/ RUNTYP,EXETYP,NEVALS,NGLEVL,NHLEVL
      DATA CHECK/8HCHECK    /
c
c     FMO restart initialisation
c
c     reserve the first record for record indexing by writing something there
c
c     skipfini is true for the F40 broadcast option when the grand master (GM)
c     broadcasts to everybody, so only the GM needs to read F40 restart info.
c     Layer restarts will probably be messed up.
c
      nts0=0
      if(irststp.eq.0.or.skipfini) then
        irstlay=1
        modrst=0
c       Setting modrst to 1 means the density will be read from the 1st rec. 
c       This is now changed to the unified modrst (in inidfmo).
        if(skipfini) modrst=1
        iodfmo(NDAR30a+1)=modorb
        iodfmo(NDAR30a+2)=l1fmo
        iodfmo(NDAR30a+3)=irstlay
        iodfmo(NDAR30a+4)=modrst
        iodfmo(NDAR30a+5)=JRECST(IDAFMO/10)
        iodfmo(NDAR30a+6)=nts2
        call icopy(nlayer,idmrec0,1,iodfmo(NDAR30a+7),1)
c       if(icanuse1) 
        call rawrites(IDAFMO,iodfmo,iodfmo,NDAR30fa,NDAR30fa,1,0)
        if(iodfmo(1).ne.1.and.EXETYP.ne.CHECK.and.icanuse1) then
c         write(6,*) 'iolib.src has been changed! No restarts.'
          if(maswrk) write(6,*) 'The F40 structure prohibits restarts.'
c         call abrt
        endif
      else
c       here is a potentially dangerous part (if iolib is later altered)
c       emulate record indices, by giving 1 as the first physical record
c       where the first logical record starts. It is assumed that no other
c       index matters (for the first record).
        if(icanuse1) iodfmo(1)=1 
        call rareads(IDAFMO,iodfmo,iodfmo,NDAR30fa,1,0)
c
c       check if the file provided is acceptable.
c
        modorb0=iodfmo(NDAR30a+1)
        l1fmo0= iodfmo(NDAR30a+2)
        irstlay=iodfmo(NDAR30a+3)
        modrst= iodfmo(NDAR30a+4)
        JRECST(IDAFMO/10)=iodfmo(NDAR30a+5)
        nts0= iodfmo(NDAR30a+6)
        call icopy(nlayer,iodfmo(NDAR30a+7),1,idmrec0,1)
c       l1fmo0 is 0 for FMO0 when l1fmo is not computed, so by writing zero
c       we skip the check here as a temporary measure. 
        if((l1fmo.ne.l1fmo0.and.l1fmo0.ne.0).or.modorb.ne.modorb0.or.
     *     irstlay.gt.nlayer) then
c         refuse to use the restart file if
c         1) basis set size for the whole molecule mismatches
c         2) orbital/density/energy flag mismatches
c         Note: it is wrong to think that for example one can read in 
c         F40 computed with energies and reuse without (in general).
c         The reason is that writing to this file allocates certain size
c         and mixing this size is not going to work. Note that if one runs
c         dimers in some cases it may work. 
c         Sometimes the error below comes from forgetting to adjust $FMOLMO
c         groups so that l1fmo is not computed correctly.
          write(6,*) 'Incompatible restart file: l1,l1, modorb,modorb=',
     *               l1fmo,l1fmo0,modorb,modorb0
c         call abrt
        endif
c       modify modrst if needed (mostly to use RHF, not MP2 densities for
c       restart jobs.
        if(iand(nguess,8192).ne.0) modrst=modrst-1
c       if(iand(nguess,65536).ne.0) irstlay=1
        if(irestl.ne.0) irstlay=irestl
c       write(6,*) 'modrst is',modrst
      endif
c     Care should be taken lest FMOENM and FMOEND be read during geometry 
c     optimisations (they should not be, except possibly the first step). 
c     During optimisations irststp is 2, so the present code seems correct.
      ilay=0
c     layer structure may be broken? irstlay.gt.0? 
      if(irststp.gt.2.or.irstlay.gt.1.or.(nbsse.eq.3.and.irststp.gt.0))
     *  call eminout(0,ilay,ifgfmo0,layfrg,emon,outpune,some)
      if((irststp.ge.4.or.irstlay.gt.1).and.nbody.gt.1) 
     *  call edin(nfg2,nen,layfrg,edim,some)
      if((irststp.ge.9.or.irstlay.gt.1).and.nbody.gt.2) 
     *  call etin(nfg3,nent,layfrg,etrim,some)
      return
      end
c
C*MODULE fmogrd  *DECK grdrst
      subroutine grdrst(imode,nfthes,nde,de,nem,em,nie,ie,lie,IDAFMO,
     *                  ndar30,niod,iodfmo,icanuse1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      dimension de(nde),em(nem),ie(nie),iodfmo(*)
      logical lie,GOPARR,DSKWRK,MASWRK,icanuse1
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /RAIOLN/ JRECLN(10),JRECST(10)
c     
c     read (imode=0) or write (imode=1) gradient restart file.
c     de gradient
c     em dipole moments 
c     ie RESDIM record indices
c     lie says if IODFMO was updated and should be saved.
c     
      CALL SEQREW(NFThes)
      if(imode.eq.0) then
        READ(nfthes,END=604) nde0,nem0,nie0
        if(nde.ne.nde0.or.nem.ne.nem0.or.nie.ne.nie0) then
          if(maswrk) write(iw,*) 'Wrong restart:',nde,nde0,nem,nem0,
     *                                            nie,nie0
          call abrt
        endif
        READ(nfthes) de
        if(nem.ne.0) READ(nfthes) em
        if(nie.ne.0) READ(nfthes) ie 
c       There is no need to read iodfmo here: it is done elsewhere.
        if(maswrk) write(iw,*) 'Gradient restart is successful!'
c       write(iw,*) de
        write(iw,*) 'wwwnie=',nie,(ie(i),i=1,nie)
  604   continue
c       write(iw,*) 'no grad restart'
      else 
        write(NFThes) nde,nem,nie
        write(NFThes) de
        if(nem.ne.0) write(nfthes) em
        if(nie.ne.0) then
          write(nfthes) ie 
          if(lie) then
            iodfmo(NDAR30+5)=JRECST(IDAFMO/10)
            if(.not.icanuse1) call abrtx("Gradient restart error")
            call rawrites(IDAFMO,iodfmo,iodfmo,niod,niod,1,0)
          endif 
c         Indexing information should be saved if IDAFMO was updated.
        endif
        if(maswrk) write(iw,*) 'Updated gradient restart file.'
c       write(iw,*) de
c       write(iw,*) 'wwwnie=',nie,(ie(i),i=1,nie)
      endif
      return
      END
C*MODULE fmogrd  *DECK hssini
      subroutine hssini(hssinv,imin,jmin,imax,jmax)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      dimension hssinv(*)
      PARAMETER (three=3.0D+00,zero=0.0D+00)
      logical GOPARR,DSKWRK,MASWRK
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
c
c     Generate initial inverse Hessian for the FMO methods.
c     The Hessian is stored as a lower triangular matrix divided over nodes.
c 
      diagval=three 
c     hssinv(iloop)=1.0D+00
      iloop=0
      jstart=jmin
      do i=imin,imax
        jend=i
        if(i.eq.imax) jend=jmax
        if(jend.eq.i) jend=jend-1
c       Off-diagonal 
        do j=jstart,jend
          iloop=iloop+1
          hssinv(iloop)=zero
        enddo 
c       diagonal 
        if(i.ne.imax.or.jmax.eq.imax) then
          iloop=iloop+1
          hssinv(iloop)=diagval
c         3 comes from 1/(1/3) where 1/3 is used in other parts of GAMESS.
        endif
        jstart=1
      enddo
      if(maswrk) write(iw,9000) diagval
c     write(6,*) 'inv hess',imin,jmin,imax,jmax
c     call prtri(hssinv,imax)
      return
 9000 format(1x,'Using diagonal inverse hessian guess of ',F4.1)
      END
c
C*MODULE fmogrd  *DECK pMTAV
      subroutine pMTAV(A,V,AV,n,imin,jmin,imax,jmax)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (zero=0.0D+00)
      dimension a(*),v(n),av(n) 
c
c     Multiply a symmetric (lower tringular) matrix A(n x n) times
c     a vector V(n) and store the results into a vector AV. 
c     A is divided over nodes.
c
      call vclr(av,1,n) 
      iloop=0
      jstart=jmin
      do i=imin,imax
        jend=i
        if(i.eq.imax) jend=jmax
        vi=v(i)
        sum=zero
        do j=jstart,jend
          iloop=iloop+1
          sum=sum+a(iloop)*v(j)
          av(j)=av(j)+a(iloop)*vi
c         write(6,*) 'wwwij',i,j,sum,a(iloop),iloop,vi,v(j)
        enddo
c       subtract double counting of the diagonal
        if(jend.eq.i) sum=sum-a(iloop)*vi
        av(i)=av(i)+sum
        jstart=1
      enddo
      call ddi_gsumf(2422,av,n)
c     write(6,*) 'inv hess',n,imin,jmin,imax,jmax
c     call prtri(a,n)
c     write(6,*) 'wwwv',(v(i),i=1,n)
c     write(6,*) 'wwwav',(av(i),i=1,n)
      return
      END
C*MODULE fmogrd  *DECK updhss
      subroutine updhss(update,hssinv,imin,jmin,imax,jmax,
     *                  cfdif,hdg,u,dxdg,dghdg)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      dimension hssinv(*),cfdif(*),hdg(*),u(*)
      data bfgs/8HBFGS    /,dfp/8HDFP     /
c
c     Numeric updates of the inverse Hessian matrix.
c     Supported types:
c     update=bfgs, ref: JCC, 18, 1011 (1997).
c     update=dfp,  ref: JCC, 18, 1011 (1997).
c
      iloop=0
      jstart=jmin
      do i=imin,imax
        jend=i
        if(i.eq.imax) jend=jmax
c
c       BFGS updating scheme
c
        if(update.eq.bfgs) then
          do j=jstart,jend
            iloop=iloop+1
            hssinv(iloop)=hssinv(iloop)+ ( cfdif(i)*cfdif(j) )/dxdg
     *                                 - ( hdg(i)*hdg(j) )/dghdg
     *                                 + dghdg* (u(i)*u(j))
          enddo
c
c       DFP updating scheme
c
        else if(update.eq.dfp) then
          do j=jstart,jend
            iloop=iloop+1
            hssinv(iloop)=hssinv(iloop)+ ( cfdif(i)*cfdif(j) )/dxdg
     *                                 - ( hdg(i)*hdg(j) )/dghdg
          enddo
        else
          call abrtx("OPTFMO method failure")
        endif
        jstart=1
      enddo
c     psb update (optional)
c       do i = 1,3*natfmo
c        do j = 1,3*natfmo
c         hesinv(i,j) = hesinv(i,j)
c    1                + ( gdif(i) - hdg(i) )*cdif(j)/dxdx
c    1                + ( gdif(j) - hdg(j) )*cdif(i)/dxdx
c    1                - cdif(i)*cdif(j)*( dxdg - dxhdx )/(dxdx*dxdx)
c         hesinv(j,i) = hesinv(i,j)
c        end do
c       end do
c     write(6,*) 'inv hess',imin,jmin,imax,jmax
c     call prtri(hssinv,imax)
      return
      END
c
C*MODULE fmogrd  *DECK stepsize
      subroutine stepsize(ediff,slope0,slopeold,step0,stepold,stepnew)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (two=2.0D+00,three=3.0D+00)
c
c     Determine the next step size during geometry optimisation,
c     using cubic interpolation methods.
c
c     The result is in stepnew, the rest are input arguments.
c
c     written by T. Ishida
c
      ramda = 0.1d+00
      rmiu = 0.5d+00
c
      aa = -two*ediff + (slope0+slopeold)*(stepold-step0)
      aa = aa/(stepold-step0)**3
      bb = three*ediff - (two*slope0+slopeold)*(stepold-step0)
      bb = bb/(stepold-step0)**2
      cc = slope0
c
      step=step0 + ( -bb + sqrt( bb*bb - three*aa*cc ) )/(three*aa)
      stepnew=max( ramda*stepold , min(step,rmiu*stepold) )
c
c         optional  .........................................................
c         step size is determind by quadratic interpolation methods
c         note : stepc is the calculated new step size
c                cc = (fb - fa - (slopea)*stepb)/(stepb**2)
c                stepc = - ( (slopea)/(2.0D+00*cc) )
c                stepc = max( ramda*stepb , min(stepc,rmiu*stepb) )
c        ....................................................................
c
      return 
      end 
C*MODULE fmogrd  *DECK rstopt1
      subroutine rstopt1(natfmo,fmocp,fmogep,irest)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      logical GOPARR,DSKWRK,MASWRK
      character*8 grpnam
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      dimension fmocp(3,natfmo),fmogep(3,natfmo)
      data grpnam/' $OPTRST'/
c
c     Read(irest=0)/write(irest=1) FMO optimisation restart. 
c     Set irest to 1 if read successfully.
c
      if(irest.eq.0) then
        CALL SEQREW(IR)
        CALL FNDGRP(IR,grpnam,IEOF)
        IF(IEOF.EQ.0) then
          if(maswrk) then
            read(ir,9020) fmocp 
            read(ir,9030) fmogep
          endif 
          if(goparr) then
            CALL DDI_BCAST(2422,'F',fmocp,3*natfmo,master)
            CALL DDI_BCAST(2422,'F',fmogep,3*natfmo,master)
          endif
          if(maswrk) write(iw,9200)
          irest=1
        endif 
      else
        write(ip,9005) grpnam
        write(ip,9020) fmocp 
        write(ip,9030) fmogep
        write(ip,9100)
      endif
      return 
 9005 format(1x,'Restart data for RUNTYP=OPTFMO, METHOD=HSSUPD or BFGSL'
     *     ,/A8)
 9020 format(6f13.8)
 9030 format(8f10.7)
 9100 format(1x,'$END')
 9200 format(/1x,'Successfully read optimisation restart data.')
      end 
C*MODULE fmogrd  *DECK rstopt2
      subroutine rstopt2(natfmo,frmsp,p,irest)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      logical GOPARR,DSKWRK,MASWRK
      character*8 grpnam
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      dimension p(3,natfmo)
      data grpnam/' $OPTRST'/
c
c     Read(irest=0)/write(irest=1) FMO optimisation restart.
c     Set irest to 1 if read successfully.
      if(irest.eq.0) then
        CALL SEQREW(IR)
        CALL FNDGRP(IR,grpnam,IEOF)
        IF(IEOF.EQ.0) then
          if(maswrk) then
            read(ir,9020) frmsp
            read(ir,9030) p
          endif
          if(goparr) then
            CALL DDI_BCAST(2422,'F',frmsp,1,master)
            CALL DDI_BCAST(2422,'F',p,3*natfmo,master)
          endif
          if(maswrk) write(iw,9200)
          irest=1
        endif
      else
        write(ip,9005) grpnam
        write(ip,9020) frmsp 
        write(ip,9030) p
        write(ip,9100)
      endif
      return
 9005 format(1x,'Restart data for RUNTYP=OPTFMO, METHOD=CG',/A8)
 9020 format(1x,f14.10)
 9030 format(8f10.7)
 9100 format(1x,'$END')
 9200 format(/1x,'Successfully read optimisation restart data.')
      end
C*MODULE fmogrd  *DECK stepofmo
      SUBROUTINE stepofmo(ENERGY,olde,depre,step,stepmin,stepmax,stepfac
     *                   )
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      logical GOPARR,DSKWRK,MASWRK
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
c     PARAMETER (ONE=1.0D+00, TWO=2.0D+00, THREE=3.0D+00, FOUR=4.0D+00,
c    *           TMONE=1.0D-01, PT75=0.75D+00)
c
      DEACT = ENERGY-OLDE
      if(DEPRE.ne.0.0D+00) then
        RATIO = DEACT/DEPRE
C
C       THE STARTING POINT FOR DXMAXT IS THE PREVIOUS VALUE
C       IF THIS WAS A GOOD STEP, INCREASE THE ALLOWED STEP LENGTH
C       REDUCE IT IF IT WAS A POOR STEP
C
c       IF(RATIO.LE.TMONE .OR. RATIO.GE.THREE)        step=step/TWO
c       IF(RATIO.GE.PT75 .AND. RATIO.LE.(FOUR/THREE)) step=step*SQRT(TWO)
c       IF(ABS(RATIO-ONE).LT.TMONE)                   step=step*SQRT(TWO)
        if(ratio.lt.-0.25D-00) step=step*stepfac
        if(ratio.gt. 1.25D+00) step=step/stepfac
      else
        RATIO=1.0D+00
      endif
      step=MAX(step,stepmin)
      step=MIN(step,stepmax)
      IF(MASWRK) WRITE(IW,9060) DEACT,DEPRE,RATIO,step
c
      return
 9060 FORMAT(10X,'   ACTUAL ENERGY CHANGE WAS',F15.10/
     *       10X,'PREDICTED ENERGY CHANGE WAS',F15.10,' RATIO=',F7.3,/,
     *       10X,'Current step size is       ',F9.4)
      end
C
C*MODULE FMOGRD  *DECK HOPDER 
C>
C>     @brief HOP derivative
C>
C>     @details Calculate HOP derivative. 
C>
C>     @author Takeshi Nagata
C>
      SUBROUTINE HOPDER(MODE,IDAM,IDAD,IDAT,L1,L2,DENS,WRK1,WRK2,FCM)
      use mx_limits, only: mxgtot,mxsh,mxatm
      USE DFTBPB_MOD,ONLY: PERIOD,STRESS
C
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (ZERO=0.0D+00,ONE=1.0D+00,thrush=1.0D+04)
      LOGICAL GOPARR,DSKWRK,MASWRK
      LOGICAL LGRAD,LHESS
      LOGICAL DFTBFL,SCC,SRSCC,DFTB3,DAMPXH,LCDFTB
      LOGICAL SAVGOP
C
      DIMENSION DENS(1),FCM(*)
      DIMENSION WRK1(1),WRK2(1),idamdt(3)
      COMMON /DFTB  / DFTBFL,SCC,SRSCC,DFTB3,DAMPXH,LCDFTB
      COMMON /FMCOM / X(1)
      COMMON /FMODIM/ MAXBND,MAXKND,MAXCBS,MAXCAO,MAXVEC,MAXL1,MAXNAT,
     *                MAXABD,MAXBAS,MAXBBD,MAXLMO,MAXSLO,MAXABD2,maxrij
      common /fmoopt/ espsca(9),RESPAP(2),rESPPC(2),rESDIM,restri(4),
     *                rcorsd,respct,convfg,cnvdmp,coroff,rflmo(4),
     *                orshft,orshft2,cnvafo,ascreen(4),IXESP,mxitfg,
     *                nguess,NBSSE,modorb,modpar,irststp,irstlay,nprfmo,
     *                nfmopal,modprp,maxl1c,ipieda,modgrd,modesp,ivmul,
     *                modlmo,nopden,mofock,modfd,modfmm,ncentm,ndualb,
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
     *                LPOPDMAT,lidmpnt,liddpnt,livmpnt,liactfg,lcrfrg,
     *                lzlmfrgv,lYlmfrgv,lndtfrg,lf_mm,lg_mm,lmaxl30,
     *                libuffg,lindatp,lsmon,lsdim,ledimfed,lexcit3d,
     *                lifgfret
      COMMON /GRAD  / DE(3,MXATM)
      Common /fmoinf/ nfg,nlayer,natfmo,nbdfg,naotyp,nbody,nsegm
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON /NSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),
     *                CF(MXGTOT),CG(MXGTOT),CH(MXGTOT),CI(MXGTOT),
     *                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),
     *                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /WFNOPT/ SCFTYP,VBTYP,DFTYPE,TDDFTYP,CITYP,CCTYP,
     *                MPLEVL,MPCTYP
      data rnone/8HNONE    /,UHF/8HUHF     /
      data RHF/8HRHF     /
C
C     Hereafter calculate HOPDER terms
C
C     MODE = 0: GRADIENT
C     MODE = 1: CPHF 
C     MODE = 2: CPHF and nder=2
C     MODE = 3: Hessian
C
      CALL DERCHK(NDER)
      valtmp = IDAD
C 
C     Call to this routine is made in grd1.src and FMOCHF calculation
C     note 1: when this routine is called in fmo.src, dscal must be 
C             applied to DE. 
C     note 2: when this routine is called in grd1.src, nothing should be 
C             done for DE (DE is already divided)
C 
      IF (MODE.EQ.0) THEN
        IF (GOPARR) CALL DSCAL(3*NAT,ONE/NPROC,DE,1)
        IF (GOPARR.AND.PERIOD) CALL DSCAL(9,ONE/NPROC,STRESS,1)
c       IF (MPLEVL.EQ.2.AND.IAND(MODGRD,32).NE.0) THEN
        IF ((MPLEVL.EQ.2.OR.TDDFTYP.NE.RNONE).AND.IAND(MODGRD,32).NE.0)
     *  THEN
          CALL DAREAD(IDAF,IODA,DENS,L2,16,0)
        ELSE
          CALL READHFD(DENS,wrk1,l2)
        END IF
c       write(6,*) 'read records in HOPDER'
c       call prtril(dens,l1) 
      END IF
C 
      LHESS = NDER.eq.2.and.mode.ne.2.and.mode.ne.4
C     Memory allocation
C 
      L3 = L1*L1
      L2NAT3 = L2*NAT*3
      CALL VALFM(LOADFM)
      LPS  = LOADFM + 1
C     LS   = LPS    + L2NAT3 
      LS   = LPS    + L2NAT3 + 3
      IF(DFTBFL.AND.mode.eq.2) LS   = LPS    + L2NAT3*2 + 3
      if(LHESS) LS = LS + natfmo * natfmo * 9
      LSS  = LS     + L2
      LDD  = LSS    + L3
      LROT = LDD    + MAXCBS*MAXCBS+MAXCBS
      LQ   = LROT   + MAXCBS*MAXCAO
      LSCR = LQ     + L3
      if(nder.eq.2) LSCR = LQ     + L3 * 4
      if(mode.eq.4.and.nder.eq.2) LSCR = LQ     + L3 * 5
      LAST = LSCR   + L1*8
C 
      LFCMWRK = LPS
      IF(LHESS) LFCMWRK = LPS + L2NAT3 + 3
C     ltest= last
C     last = ltest +  nat*3
C 
      NEED = LAST- LOADFM -1
      CALL GETFM(NEED)
C
      if(nder.eq.2) CALL DAREAD(IDAF,IODA,X(LQ+l3*3),L3,15,0)
      if(nder.eq.2.AND.mode.eq.4) 
     *    CALL DAREAD(IDAF,IODA,X(LQ+l3*4),L3,19,0)
      if(LHESS) CALL VCLR(X(LFCMWRK),1,natfmo*natfmo*9)

C     call dcopy(nat*3,de,1,x(ltest),1)
C 
C     CALCULATE OVERLAP MATRIX, X(LS) 
C 
      IF (MODE.EQ.0.or.MODE.eq.2.or.MODE.eq.4) THEN
        CALL DAREAD(IDAF,IODA,X(LS),L2,12,0)
      ELSEIF (MODE.EQ.1) THEN
        DUMMY  = ZERO 
        IDUMMY = 0
        IF (.NOT.DFTBFL) THEN
          CALL COOVLP(0,X(LS),DUMMY,L1,L1,L2,NAT,
     *                MXGTOT,NSHELL,EX,CS,CP,CD,CF,CG,CH,CI,
     *                KSTART,KATOM,KTYPE,KNG,KLOC,KMIN,KMAX,
     *                1,1,DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,
     *                DUMMY,IDUMMY,IDUMMY,IDUMMY,IDUMMY,IDUMMY,IDUMMY,
     *                IDUMMY,C,C)
        ELSE
          CALL VCLR(X(LS),1,L2)
          CALL DFTB_HAM_OVER(X(LS),X(LS),.FALSE.,.TRUE.,.FALSE.)
        END IF
      ELSE
        CALL ABRTx("Failure in HOPDER")
      END IF
      CALL CPYTSQ(X(LS),X(LSS),L1,1)
C     DSSS=TRACEP(DENS,X(LS),L1)
C     DSSS=DSSS/NA
C     if (maswrk) write(iw,*) 'TOT DS', DSSS 
      ! first arg is dummy 
      ! supposed that ncursh is 0
      if(mode.ne.2.and.mode.ne.4) then
        IF (.NOT.DFTBFL) THEN
         CALL EXTSDER(WRK1,X(LPS),L1,L2,.TRUE.,.TRUE.,0)
         IF(GOPARR)  CALL DDI_GSUMF(1666,X(LPS),L2*3*NAT)
        ELSE
         CALL DFTB_EXTSDER(X(LPS),L2,1,NAT)
        END IF
      else if(mode.eq.2.or.mode.eq.4) then
        call dcopy(L2*3*NAT,FCM,1,x(lps),1) 
        if(DFTBFL) call vclr(x(lps+L2NAT3),1,L2NAT3)
      end if
C     IF (MASWRK) THEN
C       WRITE(6,*) 'DSINT in HOPDER', ICURFG
C       WRITE(6,'(8F12.8)') (X(LPS+I-1),I=L2+1,L2*2)
C     END IF
C 
C 
C     CALCULATE HOP DERIVATIVE TERM TWO*S^a*P*S*DENS
C               (P is the hop matrix)
C      LPS: overlap derivative, LSS: overlap matrix
C     DENS: density matrix or symmetrized Z-VECTOR
C 
c     X(LPS+L2NAT3  ) = IDAM
c     X(LPS+L2NAT3+1) = IDAD
c     X(LPS+L2NAT3+2) = IDAT
      idamdt(1)=IDAM
      idamdt(2)=IDAD
      idamdt(3)=IDAT
C 
      LGRAD = .TRUE.
      CPHF  = RNONE
      if(mode.eq.2.or.mode.eq.4) then
         if(nder.ne.2) call abrtx("Internal error in HOPDER")
         CPHF  = RHF 
         if(mode.eq.4) CPHF = UHF
         CALL DAREAD(IDAF,IODA,x(LQ),L1,17,0)
         ncount = 0
         do iii = 1, NQMT
           if(x(LQ+iii-1) .gt. thrush) ncount = ncount + 1
         end do
c        X(LPS+L2NAT3  ) = ncount
         idamdt(1)=ncount
      end if
C
      SAVGOP=GOPARR
      IF(dftbfl.and.goparr.and.mode.ne.0) GOPARR=.false.
C
      CALL FMOHOP(L1,L2,X(LPS),DENS,X(LSS),X(LQ),X(LDD),
     *            X(LSCR),WRK1,WRK2,X(LIABDFG),X(LJABDFG),X(LIDXCAO)
     *           ,X(LIAGLOB),X(LNCBS),X(LNCAO),X(LIAPRJO),X(LJAPRJO),
     *            X(LSHIFTB),X(LCOREAO),X(LFMOC),X(LROT),X(LLOCFMO),
     *            NSHELL,KATOM,KTYPE,KLOC,KMIN,KMAX,.TRUE.,LGRAD,CPHF,
     *            idamdt)
C 
C     CALL DAREAD(IDAF,IODA,X(LQ),L3,15,0)
C     CALL TFSQU(WRK2,X(LSS),X(LQ),WRK1,L1,L1)
C     CALL PRSQ(WRK2,L1,L1,L1)
C 

      if(nder.eq.2.and.mode.eq.0) then
        neh = natfmo * natfmo * 9
        if(MASWRK) call daxpy(neh,valtmp,x(LFCMWRK),1,FCM,1)
      end if
C     CALL DAXPY(3*nat,-one,de,1,x(ltest),1)
C     IF(GOPARR) CALL DDI_GSUMF(1667,x(ltest),3*NAT)
      if(mode.eq.2.and.dftbfl.and.nder.eq.2) then
       do ixyz=1,nat*3
        LPNT = LPS + L2NAT3 + L2 * (ixyz-1)
        call daxpy(l2,one,x(lpnt),1,wrk1(l2*(ixyz-1)+1),1)
       end do
      end if
C     if (maswrk) then
C       if (mode.eq.0) then 
C         write(iw,*) 'CONTRIBUTION OF HOP TO GRAD'
C       else
C         write(iw,*) 'CONTRIBUTION OF HOP TO CPHF'
C       end if
C       ic = 0
C       do II = 1, NAT
C         write(iw,'(3F12.8)') x(ltest+ic),x(ltest+ic+1),x(ltest+ic+2)
C         ic = ic + 3
C       end do
C     end if
C     IF(MODE.EQ.1.AND.GOPARR) CALL DSCAL(3*NAT,ONE/NPROC,DE,1)
      IF(MODE.EQ.0.AND.GOPARR) CALL DDI_GSUMF(1666,DE,3*NAT)
      IF(MODE.EQ.0.AND.GOPARR.AND.PERIOD) CALL DDI_GSUMF(1666,STRESS,9)
C
      GOPARR=SAVGOP

      CALL RETFM(NEED)
      RETURN
      END
C*MODULE FMOGRD  *DECK HOPSDER
      SUBROUTINE HOPSDER(L1,L2,ILOC,NAO,PARRGO,SDER,DD,DENS,SS,
     *                   WRK1,WRK2,WRK3)
      use mx_limits, only: mxatm
      USE DFTBPB_MOD,ONLY: PERIOD
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (ZERO=0.0D+00,ONE=1.0D+00,TWO=2.0D+00)
      LOGICAL GOPARR,DSKWRK,MASWRK
      LOGICAL PARRGO
      DIMENSION SDER(L2,3,NAT), DD(1), DENS(1), SS(L1,*)
      DIMENSION WRK1(L1,*), WRK2(L1,*),WRK3(*)
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /GRAD  / DE(3,MXATM)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
C 
C     SDER is Overlap derivative, DD is projection matrix
C     DENS is Density matrix (density or symmetrized ZVEC)
C     SS   is Overlap matrix in the square form  
C 
      CALL CPYTSQ(DENS,WRK3,L1,1)
      CALL DGEMM('N','N',L1,NAO,L1,ONE,WRK3,L1,SS(1,ILOC),L1,ZERO,
     *           WRK1,L1)

      CALL CPYTSQ(DD,WRK3,NAO,1)
      CALL DGEMM('N','N',L1,NAO,NAO,ONE,WRK1,L1,WRK3,NAO,ZERO,
     *           WRK2(1,ILOC),L1)
C
c     old implementation
C     NLOC = ILOC + NAO - 1
C     DO I = 1, NAT
C       DO IC = 1, 3
C         CALL CPYTSQ(SDER(1,IC,I),WRK1,L1,1)
C         CONT = ZERO
C         DO IAO = ILOC, NLOC
C           DO IB = 1, L1
C             CONT = CONT + WRK1(IB,IAO)*WRK2(IB,IAO)
C           END DO
C         END DO

C         IF (PARRGO) THEN
C                      DE(IC,I) = DE(IC,I) + TWO * CONT
C         ELSE
C           IF(MASWRK) DE(IC,I) = DE(IC,I) + TWO * CONT
C         END IF
C       END DO
C     END DO
C
c     modified implementation
c     Assume that TRACEP be not parallelised (else a deadlock!)!
      NLOC = ILOC + NAO - 1
      CALL VCLR(WRK2,1,L1*(ILOC-1))
      CALL VCLR(WRK2(1,NLOC+1),1,(L1-NLOC)*L1)
      IF (PERIOD) THEN
        CALL SYMTRZE(WRK2,L1,L1)
        CALL DSCAL(L1*L1,2.0D+00,WRK2,1)
        IF (PARRGO.OR.MASWRK)
     *  CALL DFTB_ADDLAG(1,L1,WRK2,DUM,WRK1,DUM,DUM,DUM,DUM,
     *                   NAT,C,DE,IDUM)
      ELSE
        CALL SYMSQT(L1,WRK2,L1,WRK1,1)
        IF(PARRGO.or.MASWRK) then
          DO I = 1, NAT
            DO IC = 1, 3
              DE(IC,I)=DE(IC,I)+TWO*TRACEP(SDER(1,IC,I),WRK1,L1)
            END DO
c             write(6,7777) i,(TRACEP(SDER(1,IC1,I),WRK1,L1),ic1=1,3)
c7777         format(1x,'wwwgrd=',I4,3F15.8)
          END DO
        ENDIF
      END IF
C
      RETURN
      END
C*MODULE FMOGRD  *DECK HOPCODER 
      SUBROUTINE HOPCODER(L1,IRGHT,LEFT,ILOC,PARRGO,JJ,JAT,
     *                    NAO,NMO,CLCAO,LDC,ROTLCAO,LDR,NSHELL,
     *                    KATOM,KTYPE,KLOC,KMIN,SHIFTB,IBDTYP,
     *                    IAPRJO,JAPRJO,ISIDE,iDADE,DENS,SS,
     *                    ZAXIS,BOND,WRK1,WRK2,WRK3)
      use mx_limits, only: mxsh,mxatm
      USE DFTBPB_MOD,ONLY: PERIOD,STRFMO
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER (ZERO=0.0D+00,ONE=1.0D+00)
C
      LOGICAL GOPARR,DSKWRK,MASWRK
      LOGICAL PARRGO,GOPARS
      DIMENSION CLCAO(LDC,NMO),ROTLCAO(LDR,NMO),
     *          KATOM(*),KTYPE(*),KLOC(*),KMIN(*)
      DIMENSION SHIFTB(MAXCAO,*),IAPRJO(MAXCAO,*),JAPRJO(MAXCAO,*)
      DIMENSION iDADE(3),DENS(1),SS(L1,*),ZAXIS(3),BOND(3)
      DIMENSION WRK1(1), WRK2(1), WRK3(1),TMP(3)
C
      COMMON /FMCOM / X(1)
      common /fmodim/ maxbnd,maxknd,maxcbs,maxcao,maxvec,maxl1,maxnat,
     *                maxabd,maxbas,maxbbd,maxlmo,maxslo,maxabd2,maxrij
      common /fmopnt/ lichfg,lmulfg,lidmrec,lfrgnam,llayfrg,lindat,
     *                lnCBS,lfmozan,lfmoc,lfmomas,lizbas,liaglob,libdgh,
     *                liabdfg,ljabdfg,lnCAO,lidxCAO,liaprjo,ljaprjo,
     *                lCoreAO,lOccCor,lshiftb,liodfmo,lfmoda,lfmodb,
     *                lfmoespa,lfmoespb,llocfmo,lscffrg,lfmoscf,lrij,
     *                lpopmul,lpopmat,lialoc,lindbd,liatfrg,lindfrg,
     *                lindgfrg,lnatfrg,lnat0frg,lianfrg,lzanfrg,lcfrg,
     *                llibish,llibnsh,llibng,lindatg,lfmobuf(3),lfmode,
     *                lnumfrg,lloctat,liaoglob,lloadm,lfmoge,ldgrid,
     *                liodcfmo,ljob2grp,lfmopg,lemocdr,luntxyz,luntrot,
     *                lstonep,lmapsu,lfrgmul,lclmo,lialmo,lindlmo,
     *                latclmo,llmobdf,lfgflmo,lnfglmo,llfglmo,lpfglmo,
     *                LPOPDMAT,lidmpnt,liddpnt,livmpnt,liactfg,lcrfrg,
     *                lzlmfrgv,lYlmfrgv,lndtfrg,lf_mm,lg_mm,lmaxl30,
     *                libuffg,lindatp,lsmon,lsdim,ledimfed,lexcit3d,
     *                lifgfret
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /SYMTRY/ MAPSHL(MXSH,48),MAPCTR(MXATM,48),
     *                TT(432),INVT(48),NT

C 
C     COMPUTE DERIVATIVES OF HOP COEFFICIENTS
C 
      CALL VALFM(LOADFM)
      LPTRD  = LOADFM + 1
      LDTRD  = LPTRD  + 3*3
      LFTRD  = LDTRD  + 6*6
      LGTRD  = LFTRD  + 10*10
      LDROT  = LGTRD  + 15*15
      LAST   = LDROT  + LDR*NAO

      NEED  = LAST - LOADFM - 1
      CALL GETFM(NEED)

      DO IC = 1, 3
        CALL VECROTDER(ISKIP,IC,ZAXIS,BOND,TT)
        IF (ISKIP.EQ.0) THEN
          CALL TRMATDER(X(LPTRD),X(LDTRD),X(LFTRD),X(LGTRD))
          CALL VCLR(X(LDROT),1,LDR*NAO)
          ! note that nao is reset here
          CALL ROTCAO(0,JJ,0,JAT,0,X(LPTRD),X(LDTRD),X(LFTRD),X(LGTRD),
     *                NAO,NMO,CLCAO,LDC,X(LDROT),LDR,NSHELL,KATOM,KTYPE,
     *                KLOC,KMIN,.TRUE.)

C         write(6,*) 'ROT DERIVATive MATRIX', IC
C         call prsq(X(LDROT),NMO,NAO,LDR)
C 
          CALL VCLR(WRK3,1,NAO*NAO)
          DO IMO = 1, NMO
            if(iside.eq.0.and.iaprjo(imo,ibdtyp).ne.0.or.
     *         iside.ne.0.and.japrjo(imo,ibdtyp).ne.0) then

              BSHIFT = SHIFTB(imo,ibdtyp)
              CALL DGEMM('N','T',NAO,NAO,1,BSHIFT,X(LDROT+LDR*(IMO-1)),
     *                   LDR,ROTLCAO(1,IMO),LDR,ZERO,WRK2,NAO)
              CALL DAXPY(NAO*NAO,ONE,WRK2,1,WRK3,1)
            end if
          END DO
          CALL SQ2TRI(NAO,NAO,WRK3,WRK2,ONE)
          GOPARS = GOPARR
          GOPARR = .NOT.PARRGO.and.NPROC.GT.1
          CALL TFTRI0(WRK1,WRK2,SS(ILOC,1),WRK3,L1,NAO,L1)
          GOPARR = GOPARS
          IF (PARRGO.OR.MASWRK)
     *   CALL HOPDEG(L1,IC,iDADE,IRGHT,LEFT,DENS,WRK1,X(LFMODE),TMP(IC))
        END IF
      END DO
C
      IF (PERIOD.AND.(PARRGO.OR.MASWRK)) THEN
c       call vsub(X(LFMOC+3*(LEFT-1)),1,X(LFMOC+3*(IRGHT-1)),1,bond,1,3)
c       if(hoppbc) call pbcpair(bond(1),bond(2),bond(3),shpbc)
        DO K1 = 1, 3
          DO K2 = 1, 3
            STRFMO(K1,K2) = STRFMO(K1,K2) + 2.0D+00*TMP(K1)
     *       *bond(k2)
c    *       *(X(LFMOC+3*(IRGHT-1)+K2-1)-X(LFMOC+3*(LEFT-1)+K2-1))
          END DO
        END DO
      END IF
C
      CALL RETFM(NEED)
      RETURN
      END
C*MODULE FMOGRD  *DECK VECROTDER
      SUBROUTINE VECROTDER(ISKIP,IC,V1,V2,A)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      DIMENSION V1(3),V2(3),A(3,3),V12(3)
      PARAMETER(ONE=1.0D+00,SMALL=1.0D-12)
C
C     GENERATE A 3X3 UNITARY MATRIX THAT ROTATES VECTOR V1 TO V2.
C     THE VECTORS NEED NOT BE NORMALISED.
C
C 
C     V1 should be Z = (0,0,1)
C 
      ISKIP = 0
      VN1= V1(1)*V1(1)+V1(2)*V1(2)+V1(3)*V1(3)
      VN2= V2(1)*V2(1)+V2(2)*V2(2)+V2(3)*V2(3)
C     IF EITHER IS A NULL VECTOR SET THE MATRIX TO BE UNIT
c     IF(VN1.LT.SMALL.OR.VN2.LT.SMALL) THEN
C        CALL RUNITV(3,3,A)
C        if A is constant matrix, there is no contribution of HOPD
c        iskip = 1
c        RETURN
c     ENDIF
C
C     THE ORDER V1 AND V2 BELOW IS SUCH THAT THE DEFINITION OF THE
C     VECTOR PRODUCT AND THE ROTATION DIRECTION (ANTICLOCKWISE) AGREE.
C
      CALL VECPRD(V12,V2,V1)
      V12S=V1(1)*V2(1)+V1(2)*V2(2)+V1(3)*V2(3)
      COSW=V12S/SQRT(VN1*VN2)
C     IF V1 IS COLLINEAR TO V2: UNIT MATRIX
c     IF(ABS(COSW-ONE).LT.SMALL) THEN
C        CALL RUNITV(3,3,A)
c        iskip = 1
c        RETURN
c     ENDIF
      IF(ABS(COSW+ONE).LT.SMALL) THEN
C
C     IF V1 IS ANTICOLLINEAR TO V2, V12 IS A (NEARLY) ZERO VECTOR.
C       IN THIS CASE THERE IS AN INIFINITE NUMBER OF VV12 (LYING IN
C       THE SAME PLANE ORTHOGONAL TO BOTH V1 AND V2). A HARD-WIRED
C       CHOICE IS MADE THEN BY ORTHOGONALISING Z- OR Y-AXIS ORTH,
C       WHICHEVER IS LINEARLY INDEPENDENT.
C
c       iskip = 1
c       RETURN
        IF(ABS(V1(3)).GT.SMALL) THEN
          V12(1)=VN1*V1(2)*V1(1)
          V12(2)=VN1*V1(2)*V1(2)-ONE
          V12(3)=VN1*V1(2)*V1(3)
        ELSE
          V12(1)=VN1*V1(3)*V1(1)
          V12(2)=VN1*V1(3)*V1(2)
          V12(3)=VN1*V1(3)*V1(3)-ONE
        ENDIF
      ENDIF
      VV12=SQRT(V12(1)*V12(1)+V12(2)*V12(2)+V12(3)*V12(3))
      V12(1)=V12(1)/VV12
      V12(2)=V12(2)/VV12
      V12(3)=V12(3)/VV12
      OMEGA=ACOS(COSW)
      ICC = IC
      CALL MATNOMhop(ICC,OMEGA,vn2,v2,vv12,V12,A)
      RETURN
      END
C*MODULE FMOGRD  *DECK MATNOMhop
      SUBROUTINE MATNOMhop(IC,OMEGA,vn2,v2,vv12,v12,A)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER(ONE=1.0D+00,TWO=2.0D+00)
      DIMENSION v12(3),A(3,3)
      dimension v2(3),bn(3)
C 
C     dimension B(3,3), C(3,3), D(3,3)
C
C     GENERATE A 3X3 ROTATION MATRIX FROM THE N,OMEGA REPRESENTATION
C     WHERE N IS AN ORTH DEFINING A ROTATION AXIS AND OMEGA IS THE
C     ROTATION ANGLE. N IS AN ORTH AND AS SUCH MUST BE NORMALISED.
C     THIS ROUTINE GENERATES PURE ROTATION (EVEN PARITY).
C
C 
C     derivative of rotation matrix: v2 is supposed to be bond vector
C 
      CO=COS(OMEGA)
      CO1=ONE-CO
      SI=SIN(OMEGA)
C 
      VN2R = SQRT(VN2)
C 
      vv12_inv = ONE/vv12
      vn2_inv  = ONE/VN2R
      bn(1) = v2(1)*vn2_inv
      bn(2) = v2(2)*vn2_inv
      bn(3) = v2(3)*vn2_inv

      CALL VCLR(A,1,9)
C     A(1,1)=CO+CO1*N(1)*N(1)
C     A(1,2)=CO1*N(1)*N(2)-SI*N(3)
C     A(1,3)=CO1*N(1)*N(3)+SI*N(2)
C     A(2,1)=CO1*N(2)*N(1)+SI*N(3)
C     A(2,2)=CO+CO1*N(2)*N(2)
C     A(2,3)=CO1*N(2)*N(3)-SI*N(1)
C     A(3,1)=CO1*N(3)*N(1)-SI*N(2)
C     A(3,2)=CO1*N(3)*N(2)+SI*N(1)
C     A(3,3)=CO+CO1*N(3)*N(3)

      if (ic.eq.1) then
        ! x-element
        v12_2 = v12(1)*v12(1)
        A(1,1)= bn(1)*bn(1)*v12(2)*CO+two*v12_2*v12(2)*CO1
        A(1,1)= A(1,1)*vv12_inv
        A(1,2)=-bn(1)*bn(1)*v12(1)*CO-(v12_2-v12(2)*v12(2))*v12(1)*CO1
        A(1,2)= A(1,2)*vv12_inv
        A(1,3)= bn(1)*v12(2)*CO*CO-v12_2*SI
        A(1,3)= A(1,3)*vv12_inv
        A(2,1)= A(1,2)
        A(2,2)= v12(2)*bn(2)*bn(2)*CO-two*v12(2)*v12_2*CO1
        A(2,2)= A(2,2)*vv12_inv
        A(2,3)=-bn(1)*v12(1)*co*co-v12(1)*v12(2)*SI
        A(2,3)= A(2,3)*vv12_inv
        A(3,1)=-A(1,3)
        A(3,2)=-A(2,3)
        A(3,3)=-bn(1)*CO*vn2_inv

      elseif (ic.eq.2) then
        ! y-element
        v12_2 = v12(2)*v12(2)
        A(1,1)=-bn(1)*bn(1)*v12(1)*CO+two*v12(1)*v12_2*CO1
        A(1,1)= A(1,1)*vv12_inv
        A(1,2)= bn(2)*bn(2)*v12(2)*CO+(v12_2-v12(1)*v12(1))*v12(2)*CO1
        A(1,2)= A(1,2)*vv12_inv
        A(1,3)= bn(2)*v12(2)*CO*CO-v12(1)*v12(2)*SI
        A(1,3)= A(1,3)*vv12_inv
        A(2,1)= A(1,2)
        A(2,2)=-v12(1)*bn(2)*bn(2)*CO-two*CO1*v12(1)*v12(2)*v12(2)
        A(2,2)= A(2,2)*vv12_inv
        A(2,3)=-bn(2)*v12(1)*CO*CO-v12(2)*v12(2)*SI
        A(2,3)= A(2,3)*vv12_inv 
        A(3,1)=-A(1,3)
        A(3,2)=-A(2,3)
        A(3,3)=-bn(2)*CO*vn2_inv

      elseif (ic.eq.3) then
        ! z-element
        A(1,1)= bn(1)*bn(1)*vn2_inv
        A(1,2)= bn(1)*bn(2)*vn2_inv
        A(1,3)= bn(1)*CO*vn2_inv
        A(2,1)= A(1,2)
        A(2,2)= bn(2)*bn(2)*vn2_inv
        A(2,3)= bn(2)*CO*vn2_inv
        A(3,1)=-A(1,3)
        A(3,2)=-A(2,3)
        A(3,3)= SI*SI*vn2_inv

      else
        call abrtx("Confusion in MATNOMhop")
      end if
C 
C     B(1,1)=CO+CO1*v12(1)*v12(1)
C     B(1,2)=CO1*v12(1)*v12(2)-SI*v12(3)
C     B(1,3)=CO1*v12(1)*v12(3)+SI*v12(2)
C     B(2,1)=CO1*v12(2)*v12(1)+SI*v12(3)
C     B(2,2)=CO+CO1*v12(2)*v12(2)
C     B(2,3)=CO1*v12(2)*v12(3)-SI*v12(1)
C     B(3,1)=CO1*v12(3)*v12(1)-SI*v12(2)
C     B(3,2)=CO1*v12(3)*v12(2)+SI*v12(1)
C     B(3,3)=CO+CO1*v12(3)*v12(3)
C     CALL DGEMM('N','T',3,3,3,1.0D+00,A,3,B,3,0.0D+00,C,3)
C     CALL DGEMM('N','T',3,3,3,1.0D+00,B,3,A,3,0.0D+00,D,3)
C     write(6,*) 'TEST C AND D MATRIX'
C     CALL PRSQ(C,3,3,3)
C     CALL PRSQ(D,3,3,3)
C 
      RETURN
      END
C
C*MODULE FMOGRD  *DECK TRMATDDER
      SUBROUTINE TRMATDER(PTRD,DTRD,FTRD,GTRD)
      use mx_limits, only: mxsh,mxatm
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      LOGICAL GOPARR,DSKWRK,MASWRK
C
      DIMENSION PNAME(3),DNAME(6),FNAME(10),GNAME(15)
      dimension ptrd(3,*),dtrd(6,*),ftrd(10,*), gtrd(15,*)
C
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON /OUTPUT/ NPRINT,ITOL,ICUT,NORMF,NORMP,NOPK
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /SYMSPD/ PTR(3,144),DTR(6,288),FTR(10,480),GTR(15,720)
      COMMON /SYMTRY/ MAPSHL(MXSH,48),MAPCTR(MXATM,48),
     *                TT(432),INVT(48),NT
      COMMON /TRANSF/ XSMAL,YSMAL,ZSMAL,XNEW,YNEW,ZNEW,XP,YP,ZP
      COMMON /FRAME / U1,U2,U3,V1,V2,V3,W1,W2,W3,X0,Y0,Z0
C
      PARAMETER (ONE=1.0D+00, SQRT3=1.73205080756888D+00,
     *           SQRT5=2.23606797749979D+00,
     *           SQRT7=2.64575131106459D+00)
C
      DATA PNAME /4HX   ,4HY   ,4HZ   /
      DATA DNAME /4HXX  ,4HYY  ,4HZZ  ,4HXY  ,4HXZ  ,4HYZ  /
      DATA FNAME /4HXXX ,4HYYY ,4HZZZ ,4HXXY ,4HXXZ ,4HYYX ,
     *            4HYYZ ,4HZZX ,4HZZY ,4HXYZ /
      DATA GNAME /4HXXXX,4HYYYY,4HZZZZ,4HXXXY,4HXXXZ,4HYYYX,
     *            4HYYYZ,4HZZZX,4HZZZY,4HXXYY,4HXXZZ,4HYYZZ,
     *            4HXXYZ,4HYYXZ,4HZZXY/
C
C     ----- CALCULATE TRANSFORMS OF P, D, F, AND G FUNCTIONS
C           FOR ALL SYMETRY OPERATIONS.
C
C 
C     This is the clone of TRMAT and computes the derivatives of
C     HOP coefficients in FMO
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
      PTRD(1,N+1) = XP-X0
      PTRD(2,N+1) = YP-Y0
      PTRD(3,N+1) = ZP-Z0
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
      PTRD(1,N+2) = XP-X0
      PTRD(2,N+2) = YP-Y0
      PTRD(3,N+2) = ZP-Z0
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
      PTRD(1,N+3) = XP-X0
      PTRD(2,N+3) = YP-Y0
      PTRD(3,N+3) = ZP-Z0
  140 CONTINUE
      DO 640 IT = 1,NT
      NP = 3*(IT-1)
      ND = 6*(IT-1)
      NF = 10*(IT-1)
      NG = 15*(IT-1)
      DO 300 I = 1,6
      GO TO (160,180,200,220,240,260) ,I
  160 J = 1
      K = 1
      GO TO 280
C
  180 J = 2
      K = 2
      GO TO 280
C
  200 J = 3
      K = 3
      GO TO 280
C
  220 J = 1
      K = 2
      GO TO 280
C
  240 J = 1
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
  280 DTRD(1,ND+I) = PTRD(1,NP+J)*PTR(1,NP+K) + PTR(1,NP+J)*PTRD(1,NP+K)
      DTRD(2,ND+I) = PTRD(2,NP+J)*PTR(2,NP+K) + PTR(2,NP+J)*PTRD(2,NP+K)
      DTRD(3,ND+I) = PTRD(3,NP+J)*PTR(3,NP+K) + PTR(3,NP+J)*PTRD(3,NP+K)
      DTRD(4,ND+I) = PTRD(1,NP+J)*PTR(2,NP+K) + PTRD(2,NP+J)*PTR(1,NP+K)
     *             + PTR(1,NP+J)*PTRD(2,NP+K) + PTR(2,NP+J)*PTRD(1,NP+K)
      DTRD(5,ND+I) = PTRD(1,NP+J)*PTR(3,NP+K) + PTRD(3,NP+J)*PTR(1,NP+K)
     *             + PTR(1,NP+J)*PTRD(3,NP+K) + PTR(3,NP+J)*PTRD(1,NP+K)
      DTRD(6,ND+I) = PTRD(2,NP+J)*PTR(3,NP+K) + PTRD(3,NP+J)*PTR(2,NP+K)
     *             + PTR(2,NP+J)*PTRD(3,NP+K) + PTR(3,NP+J)*PTRD(2,NP+K)
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
C 420 FTR(1,NF+I)=DTR(1,ND+J)*PTR(1,NP+K)
C     FTR(2,NF+I)=DTR(2,ND+J)*PTR(2,NP+K)
C     FTR(3,NF+I)=DTR(3,ND+J)*PTR(3,NP+K)
C     FTR(4,NF+I)=DTR(1,ND+J)*PTR(2,NP+K)
C    1           +DTR(4,ND+J)*PTR(1,NP+K)
C     FTR(5,NF+I)=DTR(1,ND+J)*PTR(3,NP+K)
C    1           +DTR(5,ND+J)*PTR(1,NP+K)
C     FTR(6,NF+I)=DTR(2,ND+J)*PTR(1,NP+K)
C    1           +DTR(4,ND+J)*PTR(2,NP+K)
C     FTR(7,NF+I)=DTR(2,ND+J)*PTR(3,NP+K)
C    1           +DTR(6,ND+J)*PTR(2,NP+K)
C     FTR(8,NF+I)=DTR(3,ND+J)*PTR(1,NP+K)
C    1           +DTR(5,ND+J)*PTR(3,NP+K)
C     FTR(9,NF+I)=DTR(3,ND+J)*PTR(2,NP+K)
C    1           +DTR(6,ND+J)*PTR(3,NP+K)
C     FTR(10,NF+I)=DTR(4,ND+J)*PTR(3,NP+K)
C    1           + DTR(5,ND+J)*PTR(2,NP+K)
C    2           + DTR(6,ND+J)*PTR(1,NP+K)
  420 FTRD(1,NF+I)=DTRD(1,ND+J)*PTR(1,NP+K) + DTR(1,ND+J)*PTRD(1,NP+K)
      FTRD(2,NF+I)=DTRD(2,ND+J)*PTR(2,NP+K) + DTR(2,ND+J)*PTRD(2,NP+K)
      FTRD(3,NF+I)=DTRD(3,ND+J)*PTR(3,NP+K) + DTR(3,ND+J)*PTRD(3,NP+K)
      FTRD(4,NF+I)=DTRD(1,ND+J)*PTR(2,NP+K) + DTR(1,ND+J)*PTRD(2,NP+K)
     1            +DTRD(4,ND+J)*PTR(1,NP+K) + DTR(4,ND+J)*PTRD(1,NP+K)
      FTRD(5,NF+I)=DTRD(1,ND+J)*PTR(3,NP+K) + DTR(1,ND+J)*PTRD(3,NP+K)
     1            +DTRD(5,ND+J)*PTR(1,NP+K) + DTR(5,ND+J)*PTRD(1,NP+K)
      FTRD(6,NF+I)=DTRD(2,ND+J)*PTR(1,NP+K) + DTR(2,ND+J)*PTRD(1,NP+K)
     1            +DTRD(4,ND+J)*PTR(2,NP+K) + DTR(4,ND+J)*PTRD(2,NP+K)
      FTRD(7,NF+I)=DTRD(2,ND+J)*PTR(3,NP+K) + DTR(2,ND+J)*PTRD(3,NP+K)
     1            +DTRD(6,ND+J)*PTR(2,NP+K) + DTR(6,ND+J)*PTRD(2,NP+K)
      FTRD(8,NF+I)=DTRD(3,ND+J)*PTR(1,NP+K) + DTR(3,ND+J)*PTRD(1,NP+K)
     1            +DTRD(5,ND+J)*PTR(3,NP+K) + DTR(5,ND+J)*PTRD(3,NP+K)
      FTRD(9,NF+I)=DTRD(3,ND+J)*PTR(2,NP+K) + DTR(3,ND+J)*PTRD(2,NP+K)
     1            +DTRD(6,ND+J)*PTR(3,NP+K) + DTR(6,ND+J)*PTRD(3,NP+K)
      FTRD(10,NF+I)=DTRD(4,ND+J)*PTR(3,NP+K) + DTR(4,ND+J)*PTRD(3,NP+K)
     1           +  DTRD(5,ND+J)*PTR(2,NP+K) + DTR(5,ND+J)*PTRD(2,NP+K)
     2           +  DTRD(6,ND+J)*PTR(1,NP+K) + DTR(6,ND+J)*PTRD(1,NP+K)
  440 CONTINUE
      DO 620 I=1,15
      GO TO (460,470,480,490,500,510,520,530,
     1       540,550,560,570,580,590,600),I
  460 J=1
      K=1
      GO TO 610
  470 J=2
      K=2
      GO TO 610
  480 J=3
      K=3
      GO TO 610
  490 J=1
      K=2
      GO TO 610
  500 J=1
      K=3
      GO TO 610
  510 J=2
      K=1
      GO TO 610
  520 J=2
      K=3
      GO TO 610
  530 J=3
      K=1
      GO TO 610
  540 J=3
      K=2
      GO TO 610
  550 J=4
      K=2
      GO TO 610
  560 J=5
      K=3
      GO TO 610
  570 J=7
      K=3
      GO TO 610
  580 J=4
      K=3
      GO TO 610
  590 J=6
      K=3
      GO TO 610
  600 J=8
      K=2
C 610 GTR(1,NG+I)=FTR(1,NF+J)*PTR(1,NP+K)
C     GTR(2,NG+I)=FTR(2,NF+J)*PTR(2,NP+K)
C     GTR(3,NG+I)=FTR(3,NF+J)*PTR(3,NP+K)
C     GTR(4,NG+I)=FTR(1,NF+J)*PTR(2,NP+K)
C    1           +FTR(4,NF+J)*PTR(1,NP+K)
C     GTR(5,NG+I)=FTR(1,NF+J)*PTR(3,NP+K)
C    1           +FTR(5,NF+J)*PTR(1,NP+K)
C     GTR(6,NG+I)=FTR(2,NF+J)*PTR(1,NP+K)
C    1           +FTR(6,NF+J)*PTR(2,NP+K)
C     GTR(7,NG+I)=FTR(2,NF+J)*PTR(3,NP+K)
C    1           +FTR(7,NF+J)*PTR(2,NP+K)
C     GTR(8,NG+I)=FTR(3,NF+J)*PTR(1,NP+K)
C    1           +FTR(8,NF+J)*PTR(3,NP+K)
C     GTR(9,NG+I)=FTR(3,NF+J)*PTR(2,NP+K)
C    1           +FTR(9,NF+J)*PTR(3,NP+K)
C     GTR(10,NG+I)=FTR(4,NF+J)*PTR(2,NP+K)
C    1           + FTR(6,NF+J)*PTR(1,NP+K)
C     GTR(11,NG+I)=FTR(5,NF+J)*PTR(3,NP+K)
C    1           + FTR(8,NF+J)*PTR(1,NP+K)
C     GTR(12,NG+I)=FTR(7,NF+J)*PTR(3,NP+K)
C    1           + FTR(9,NF+J)*PTR(2,NP+K)
C     GTR(13,NG+I)=FTR(4,NF+J)*PTR(3,NP+K)
C    1           + FTR(5,NF+J)*PTR(2,NP+K)
C    2           +FTR(10,NF+J)*PTR(1,NP+K)
C     GTR(14,NG+I)=FTR(6,NF+J)*PTR(3,NP+K)
C    1           + FTR(7,NF+J)*PTR(1,NP+K)
C    2           +FTR(10,NF+J)*PTR(2,NP+K)
C     GTR(15,NG+I)=FTR(8,NF+J)*PTR(2,NP+K)
C    1           + FTR(9,NF+J)*PTR(1,NP+K)
C    2           +FTR(10,NF+J)*PTR(3,NP+K)
  610 GTRD(1,NG+I)=FTRD(1,NF+J)*PTR(1,NP+K)+FTR(1,NF+J)*PTRD(1,NP+K)
      GTRD(2,NG+I)=FTRD(2,NF+J)*PTR(2,NP+K)+FTR(2,NF+J)*PTRD(2,NP+K)
      GTRD(3,NG+I)=FTRD(3,NF+J)*PTR(3,NP+K)+FTR(3,NF+J)*PTRD(3,NP+K)
      GTRD(4,NG+I)=FTRD(1,NF+J)*PTR(2,NP+K)+FTR(1,NF+J)*PTRD(2,NP+K)
     1            +FTRD(4,NF+J)*PTR(1,NP+K)+FTR(4,NF+J)*PTRD(1,NP+K)
      GTRD(5,NG+I)=FTRD(1,NF+J)*PTR(3,NP+K)+FTR(1,NF+J)*PTRD(3,NP+K)
     1            +FTRD(5,NF+J)*PTR(1,NP+K)+FTR(5,NF+J)*PTRD(1,NP+K)
      GTRD(6,NG+I)=FTRD(2,NF+J)*PTR(1,NP+K)+FTR(2,NF+J)*PTRD(1,NP+K)
     1            +FTRD(6,NF+J)*PTR(2,NP+K)+FTR(6,NF+J)*PTRD(2,NP+K)
      GTRD(7,NG+I)=FTRD(2,NF+J)*PTR(3,NP+K)+FTR(2,NF+J)*PTRD(3,NP+K)
     1            +FTRD(7,NF+J)*PTR(2,NP+K)+FTR(7,NF+J)*PTRD(2,NP+K)
      GTRD(8,NG+I)=FTRD(3,NF+J)*PTR(1,NP+K)+FTR(3,NF+J)*PTRD(1,NP+K)
     1            +FTRD(8,NF+J)*PTR(3,NP+K)+FTR(8,NF+J)*PTRD(3,NP+K)
      GTRD(9,NG+I)=FTRD(3,NF+J)*PTR(2,NP+K)+FTR(3,NF+J)*PTRD(2,NP+K)
     1            +FTRD(9,NF+J)*PTR(3,NP+K)+FTR(9,NF+J)*PTRD(3,NP+K)
      GTRD(10,NG+I)=FTRD(4,NF+J)*PTR(2,NP+K)+FTR(4,NF+J)*PTRD(2,NP+K)
     1             +FTRD(6,NF+J)*PTR(1,NP+K)+FTR(6,NF+J)*PTRD(1,NP+K)
      GTRD(11,NG+I)=FTRD(5,NF+J)*PTR(3,NP+K)+FTR(5,NF+J)*PTRD(3,NP+K)
     1             +FTRD(8,NF+J)*PTR(1,NP+K)+FTR(8,NF+J)*PTRD(1,NP+K)
      GTRD(12,NG+I)=FTRD(7,NF+J)*PTR(3,NP+K)+FTR(7,NF+J)*PTRD(3,NP+K)
     1             +FTRD(9,NF+J)*PTR(2,NP+K)+FTR(9,NF+J)*PTRD(2,NP+K)
      GTRD(13,NG+I)=FTRD(4,NF+J)*PTR(3,NP+K)+FTR(4,NF+J)*PTRD(3,NP+K)
     1             +FTRD(5,NF+J)*PTR(2,NP+K)+FTR(5,NF+J)*PTRD(2,NP+K)
     2           +FTRD(10,NF+J)*PTR(1,NP+K)+FTR(10,NF+J)*PTRD(1,NP+K)
      GTRD(14,NG+I)=FTRD(6,NF+J)*PTR(3,NP+K)+FTR(6,NF+J)*PTRD(3,NP+K)
     1             +FTRD(7,NF+J)*PTR(1,NP+K)+FTR(7,NF+J)*PTRD(1,NP+K)
     2           +FTRD(10,NF+J)*PTR(2,NP+K)+FTR(10,NF+J)*PTRD(2,NP+K)
      GTRD(15,NG+I)=FTRD(8,NF+J)*PTR(2,NP+K)+FTR(8,NF+J)*PTRD(2,NP+K)
     1             +FTRD(9,NF+J)*PTR(1,NP+K)+FTR(9,NF+J)*PTRD(1,NP+K)
     2           +FTRD(10,NF+J)*PTR(3,NP+K)+FTR(10,NF+J)*PTRD(3,NP+K)
  620 CONTINUE
  640 CONTINUE
      IF (NORMF .EQ. 1 .AND. NORMP .EQ. 1) GO TO 750
      DO 740 IT = 1,NT
      ND = 6*(IT-1)
      NF = 10*(IT-1)
      NG = 15*(IT-1)
      DO 660 I = 1,6
      IF (I .GT. 3) GO TO 650
C     DTR(4,ND+I) = DTR(4,ND+I)/SQRT3
C     DTR(5,ND+I) = DTR(5,ND+I)/SQRT3
C     DTR(6,ND+I) = DTR(6,ND+I)/SQRT3
      DTRD(4,ND+I) = DTRD(4,ND+I)/SQRT3
      DTRD(5,ND+I) = DTRD(5,ND+I)/SQRT3
      DTRD(6,ND+I) = DTRD(6,ND+I)/SQRT3
      GO TO 660
C
C 650 DTR(1,ND+I) = DTR(1,ND+I)*SQRT3
C     DTR(2,ND+I) = DTR(2,ND+I)*SQRT3
C     DTR(3,ND+I) = DTR(3,ND+I)*SQRT3
  650 DTRD(1,ND+I) = DTRD(1,ND+I)*SQRT3
      DTRD(2,ND+I) = DTRD(2,ND+I)*SQRT3
      DTRD(3,ND+I) = DTRD(3,ND+I)*SQRT3
  660 CONTINUE
      DO 690 I=1,10
      IF(I.GT.3) GO TO 670
C     FTR(4,NF+I)=FTR(4,NF+I)/SQRT5
C     FTR(5,NF+I)=FTR(5,NF+I)/SQRT5
C     FTR(6,NF+I)=FTR(6,NF+I)/SQRT5
C     FTR(7,NF+I)=FTR(7,NF+I)/SQRT5
C     FTR(8,NF+I)=FTR(8,NF+I)/SQRT5
C     FTR(9,NF+I)=FTR(9,NF+I)/SQRT5
C     FTR(10,NF+I)=FTR(10,NF+I)/(SQRT5*SQRT3)
      FTRD(4,NF+I)=FTRD(4,NF+I)/SQRT5
      FTRD(5,NF+I)=FTRD(5,NF+I)/SQRT5
      FTRD(6,NF+I)=FTRD(6,NF+I)/SQRT5
      FTRD(7,NF+I)=FTRD(7,NF+I)/SQRT5
      FTRD(8,NF+I)=FTRD(8,NF+I)/SQRT5
      FTRD(9,NF+I)=FTRD(9,NF+I)/SQRT5
      FTRD(10,NF+I)=FTRD(10,NF+I)/(SQRT5*SQRT3)
      GO TO 690
  670 IF(I.GT.9) GO TO 680
C     FTR(1,NF+I)=FTR(1,NF+I)*SQRT5
C     FTR(2,NF+I)=FTR(2,NF+I)*SQRT5
C     FTR(3,NF+I)=FTR(3,NF+I)*SQRT5
C     FTR(10,NF+I)=FTR(10,NF+I)/SQRT3
      FTRD(1,NF+I)=FTRD(1,NF+I)*SQRT5
      FTRD(2,NF+I)=FTRD(2,NF+I)*SQRT5
      FTRD(3,NF+I)=FTRD(3,NF+I)*SQRT5
      FTRD(10,NF+I)=FTRD(10,NF+I)/SQRT3
      GO TO 690
C 680 FTR(1,NF+I)=FTR(1,NF+I)*SQRT5*SQRT3
C     FTR(2,NF+I)=FTR(2,NF+I)*SQRT5*SQRT3
C     FTR(3,NF+I)=FTR(3,NF+I)*SQRT5*SQRT3
C     FTR(4,NF+I)=FTR(4,NF+I)*SQRT3
C     FTR(5,NF+I)=FTR(5,NF+I)*SQRT3
C     FTR(6,NF+I)=FTR(6,NF+I)*SQRT3
C     FTR(7,NF+I)=FTR(7,NF+I)*SQRT3
C     FTR(8,NF+I)=FTR(8,NF+I)*SQRT3
C     FTR(9,NF+I)=FTR(9,NF+I)*SQRT3
  680 FTRD(1,NF+I)=FTRD(1,NF+I)*SQRT5*SQRT3
      FTRD(2,NF+I)=FTRD(2,NF+I)*SQRT5*SQRT3
      FTRD(3,NF+I)=FTRD(3,NF+I)*SQRT5*SQRT3
      FTRD(4,NF+I)=FTRD(4,NF+I)*SQRT3
      FTRD(5,NF+I)=FTRD(5,NF+I)*SQRT3
      FTRD(6,NF+I)=FTRD(6,NF+I)*SQRT3
      FTRD(7,NF+I)=FTRD(7,NF+I)*SQRT3
      FTRD(8,NF+I)=FTRD(8,NF+I)*SQRT3
      FTRD(9,NF+I)=FTRD(9,NF+I)*SQRT3
  690 CONTINUE
      DO 730 I=1,15
      IF(I.GT.3) GO TO 700
C     GTR(4,NG+I)=GTR(4,NG+I)/SQRT7
C     GTR(5,NG+I)=GTR(5,NG+I)/SQRT7
C     GTR(6,NG+I)=GTR(6,NG+I)/SQRT7
C     GTR(7,NG+I)=GTR(7,NG+I)/SQRT7
C     GTR(8,NG+I)=GTR(8,NG+I)/SQRT7
C     GTR(9,NG+I)=GTR(9,NG+I)/SQRT7
C     GTR(10,NG+I)=GTR(10,NG+I)*SQRT3/(SQRT5*SQRT7)
C     GTR(11,NG+I)=GTR(11,NG+I)*SQRT3/(SQRT5*SQRT7)
C     GTR(12,NG+I)=GTR(12,NG+I)*SQRT3/(SQRT5*SQRT7)
C     GTR(13,NG+I)=GTR(13,NG+I)/(SQRT5*SQRT7)
C     GTR(14,NG+I)=GTR(14,NG+I)/(SQRT5*SQRT7)
C     GTR(15,NG+I)=GTR(15,NG+I)/(SQRT5*SQRT7)
      GTRD(4,NG+I)=GTRD(4,NG+I)/SQRT7
      GTRD(5,NG+I)=GTRD(5,NG+I)/SQRT7
      GTRD(6,NG+I)=GTRD(6,NG+I)/SQRT7
      GTRD(7,NG+I)=GTRD(7,NG+I)/SQRT7
      GTRD(8,NG+I)=GTRD(8,NG+I)/SQRT7
      GTRD(9,NG+I)=GTRD(9,NG+I)/SQRT7
      GTRD(10,NG+I)=GTRD(10,NG+I)*SQRT3/(SQRT5*SQRT7)
      GTRD(11,NG+I)=GTRD(11,NG+I)*SQRT3/(SQRT5*SQRT7)
      GTRD(12,NG+I)=GTRD(12,NG+I)*SQRT3/(SQRT5*SQRT7)
      GTRD(13,NG+I)=GTRD(13,NG+I)/(SQRT5*SQRT7)
      GTRD(14,NG+I)=GTRD(14,NG+I)/(SQRT5*SQRT7)
      GTRD(15,NG+I)=GTRD(15,NG+I)/(SQRT5*SQRT7)
      GO TO 730
  700 IF(I.GT.9) GO TO 710
C     GTR(1,NG+I)=GTR(1,NG+I)*SQRT7
C     GTR(2,NG+I)=GTR(2,NG+I)*SQRT7
C     GTR(3,NG+I)=GTR(3,NG+I)*SQRT7
C     GTR(10,NG+I)=GTR(10,NG+I)*SQRT3/SQRT5
C     GTR(11,NG+I)=GTR(11,NG+I)*SQRT3/SQRT5
C     GTR(12,NG+I)=GTR(12,NG+I)*SQRT3/SQRT5
C     GTR(13,NG+I)=GTR(13,NG+I)/SQRT5
C     GTR(14,NG+I)=GTR(14,NG+I)/SQRT5
C     GTR(15,NG+I)=GTR(15,NG+I)/SQRT5
      GTRD(1,NG+I)=GTRD(1,NG+I)*SQRT7
      GTRD(2,NG+I)=GTRD(2,NG+I)*SQRT7
      GTRD(3,NG+I)=GTRD(3,NG+I)*SQRT7
      GTRD(10,NG+I)=GTRD(10,NG+I)*SQRT3/SQRT5
      GTRD(11,NG+I)=GTRD(11,NG+I)*SQRT3/SQRT5
      GTRD(12,NG+I)=GTRD(12,NG+I)*SQRT3/SQRT5
      GTRD(13,NG+I)=GTRD(13,NG+I)/SQRT5
      GTRD(14,NG+I)=GTRD(14,NG+I)/SQRT5
      GTRD(15,NG+I)=GTRD(15,NG+I)/SQRT5
      GO TO 730
  710 IF(I.GT.12) GO TO 720
C     GTR(1,NG+I)=GTR(1,NG+I)*SQRT7*SQRT5/SQRT3
C     GTR(2,NG+I)=GTR(2,NG+I)*SQRT7*SQRT5/SQRT3
C     GTR(3,NG+I)=GTR(3,NG+I)*SQRT7*SQRT5/SQRT3
C     GTR(4,NG+I)=GTR(4,NG+I)*SQRT5/SQRT3
C     GTR(5,NG+I)=GTR(5,NG+I)*SQRT5/SQRT3
C     GTR(6,NG+I)=GTR(6,NG+I)*SQRT5/SQRT3
C     GTR(7,NG+I)=GTR(7,NG+I)*SQRT5/SQRT3
C     GTR(8,NG+I)=GTR(8,NG+I)*SQRT5/SQRT3
C     GTR(9,NG+I)=GTR(9,NG+I)*SQRT5/SQRT3
C     GTR(13,NG+I)=GTR(13,NG+I)/SQRT3
C     GTR(14,NG+I)=GTR(14,NG+I)/SQRT3
C     GTR(15,NG+I)=GTR(15,NG+I)/SQRT3
      GTRD(1,NG+I)=GTRD(1,NG+I)*SQRT7*SQRT5/SQRT3
      GTRD(2,NG+I)=GTRD(2,NG+I)*SQRT7*SQRT5/SQRT3
      GTRD(3,NG+I)=GTRD(3,NG+I)*SQRT7*SQRT5/SQRT3
      GTRD(4,NG+I)=GTRD(4,NG+I)*SQRT5/SQRT3
      GTRD(5,NG+I)=GTRD(5,NG+I)*SQRT5/SQRT3
      GTRD(6,NG+I)=GTRD(6,NG+I)*SQRT5/SQRT3
      GTRD(7,NG+I)=GTRD(7,NG+I)*SQRT5/SQRT3
      GTRD(8,NG+I)=GTRD(8,NG+I)*SQRT5/SQRT3
      GTRD(9,NG+I)=GTRD(9,NG+I)*SQRT5/SQRT3
      GTRD(13,NG+I)=GTRD(13,NG+I)/SQRT3
      GTRD(14,NG+I)=GTRD(14,NG+I)/SQRT3
      GTRD(15,NG+I)=GTRD(15,NG+I)/SQRT3
      GO TO 730
C 720 GTR(1,NG+I)=GTR(1,NG+I)*SQRT7*SQRT5
C     GTR(2,NG+I)=GTR(2,NG+I)*SQRT7*SQRT5
C     GTR(3,NG+I)=GTR(3,NG+I)*SQRT7*SQRT5
C     GTR(4,NG+I)=GTR(4,NG+I)*SQRT5
C     GTR(5,NG+I)=GTR(5,NG+I)*SQRT5
C     GTR(6,NG+I)=GTR(6,NG+I)*SQRT5
C     GTR(7,NG+I)=GTR(7,NG+I)*SQRT5
C     GTR(8,NG+I)=GTR(8,NG+I)*SQRT5
C     GTR(9,NG+I)=GTR(9,NG+I)*SQRT5
C     GTR(10,NG+I)=GTR(10,NG+I)*SQRT3
C     GTR(11,NG+I)=GTR(11,NG+I)*SQRT3
C     GTR(12,NG+I)=GTR(12,NG+I)*SQRT3
  720 GTRD(1,NG+I)=GTRD(1,NG+I)*SQRT7*SQRT5
      GTRD(2,NG+I)=GTRD(2,NG+I)*SQRT7*SQRT5
      GTRD(3,NG+I)=GTRD(3,NG+I)*SQRT7*SQRT5
      GTRD(4,NG+I)=GTRD(4,NG+I)*SQRT5
      GTRD(5,NG+I)=GTRD(5,NG+I)*SQRT5
      GTRD(6,NG+I)=GTRD(6,NG+I)*SQRT5
      GTRD(7,NG+I)=GTRD(7,NG+I)*SQRT5
      GTRD(8,NG+I)=GTRD(8,NG+I)*SQRT5
      GTRD(9,NG+I)=GTRD(9,NG+I)*SQRT5
      GTRD(10,NG+I)=GTRD(10,NG+I)*SQRT3
      GTRD(11,NG+I)=GTRD(11,NG+I)*SQRT3
      GTRD(12,NG+I)=GTRD(12,NG+I)*SQRT3
  730 CONTINUE
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
C 760 WRITE (IW,9068) PNAME(I),(PTR(I,NP+J),J = 1,3)
  760 WRITE (IW,9068) PNAME(I),(PTRD(I,NP+J),J = 1,3)
      WRITE (IW,9008)
      ND = 6*(IT-1)
      WRITE (IW,9048) (DNAME(J),J = 1,6)
      WRITE (IW,9128)
      DO 780 I = 1,6
C 780 WRITE (IW,9068) DNAME(I),(DTR(I,ND+J),J = 1,6)
  780 WRITE (IW,9068) DNAME(I),(DTRD(I,ND+J),J = 1,6)
      WRITE (IW,9008)
      NF = 10*(IT-1)
      WRITE (IW,9048) (FNAME(J),J = 1,10)
      WRITE (IW,9128)
      DO 800 I = 1,10
C 800 WRITE (IW,9068) FNAME(I),(FTR(I,ND+J),J = 1,10)
  800 WRITE (IW,9068) FNAME(I),(FTRD(I,ND+J),J = 1,10)
      WRITE (IW,9008)
      NG = 15*(IT-1)
      WRITE (IW,9048) (GNAME(J),J = 1,15)
      WRITE (IW,9128)
      DO 820 I = 1,15
C 820 WRITE (IW,9068) GNAME(I),(GTR(I,ND+J),J = 1,15)
  820 WRITE (IW,9068) GNAME(I),(GTRD(I,ND+J),J = 1,15)
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
C*MODULE FMOGRD  *DECK HOPDEG
      SUBROUTINE HOPDEG(L1,IC,iDECO,IR,IL,DA,DSPS,FMODE,TMP)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      PARAMETER(ZERO=0.0D+00)
      DIMENSION iDECO(3),DA(1),DSPS(1),FMODE(3,NATFMO,*)
      COMMON /FMOINF/ NFG,NLAYER,NATFMO,NBDFG,NAOTYP,NBODY,NSEGM
C
      CONT = ZERO
      CONT = TRACEP(DA,DSPS,L1)
c     IF (ABS(DECO(1)).GE.SMALL) THEN
      IF (iDECO(1).ne.0) THEN
        TMP1 = iDECO(1)*CONT
        FMODE(IC,IR,1) = FMODE(IC,IR,1) + TMP1
        FMODE(IC,IL,1) = FMODE(IC,IL,1) - TMP1
      END IF
c     IF (ABS(DECO(2)).GE.SMALL.AND.NBODY.GT.1) THEN
      IF (iDECO(2).ne.0.AND.NBODY.GT.1) THEN
        TMP1 = iDECO(2)*CONT
        FMODE(IC,IR,2) = FMODE(IC,IR,2) + TMP1
        FMODE(IC,IL,2) = FMODE(IC,IL,2) - TMP1
      END IF
c     IF (ABS(iDECO(3)).GE.SMALL.AND.NBODY.GT.2) THEN
      IF (iDECO(3).ne.0.AND.NBODY.GT.2) THEN
        TMP1 = iDECO(3)*CONT
        FMODE(IC,IR,3) = FMODE(IC,IR,3) + TMP1
        FMODE(IC,IL,3) = FMODE(IC,IL,3) - TMP1
      END IF
      TMP = iDECO(NBODY)*CONT
      RETURN
      END
C*MODULE FMOGRD  *DECK readhfd 
      SUBROUTINE readhfd(da,db,l2)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      dimension da(l2),db(l2)
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON /WFNOPT/ SCFTYP,VBTYP,DFTYPE,TDDFTYP,CITYP,CCTYP,
     *                MPLEVL,MPCTYP
      DATA RNONE/8HNONE    /,rohf/8HROHF    /,uhf/8HUHF     /
C
c     Read RHF or ROHF total density into DA. DB is not used for RHF.
c
      if(scftyp.eq.uhf.or.scftyp.eq.rohf) then
c       read alpha and beta densities and add them up.
        IF (MPLEVL.EQ.0.AND.TDDFTYP.EQ.RNONE) THEN
          CALL DAREAD(IDAF,IODA,da,L2,16,0)
          CALL DAREAD(IDAF,IODA,db,L2,20,0)
        ELSE
          CALL DAREAD(IDAF,IODA,da,L2,418,0)
          CALL DAREAD(IDAF,IODA,db,L2,428,0)
        END IF
        CALL DAXPY(l2,1.0D+00,db,1,da,1)
      else
        IF (MPLEVL.EQ.0.AND.TDDFTYP.EQ.RNONE) THEN
          CALL DAREAD(IDAF,IODA,da,L2,16,0)
        ELSE
          CALL DAREAD(IDAF,IODA,da,L2,308,0)
        END IF
      endif
      RETURN
      END
C*MODULE FMOGRD  *DECK readfcm
      SUBROUTINE readfcm(nfthes,nc3,nc13,FCM,ddm)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      dimension fcm(nc3),ddm(nc13)
c
      read(nfthes) FCM
      if(nc13.ne.0) read(nfthes) ddm
      RETURN
      END
C*MODULE FMOGRD  *DECK clrinact
      SUBROUTINE clrinact(iactat,fmode) 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON /FMOINF/ NFG,NLAYER,NATFMO,NBDFG,NAOTYP,NBODY,NSEGM
      dimension iactat(*),fmode(3,natfmo,*)
c
      do iat=1,natfmo
        if(iactat(iat).eq.0) then
          do i=1,nbody
            do j=1,3
              fmode(j,iat,i)=0
            enddo 
          enddo 
        endif
      enddo
c     call vclr(fmode,1,3*natfmo*nbody)
c
      RETURN
      END
*MODULE FMOGRD  *DECK packat
C>
C>     @brief store atom for FMO/FD 
C>
C>     @details Store atom for FMO/FD. 
C>
C>     @author Dmitri Fedorov 
C>
      SUBROUTINE packat(mode,layfrg,indat,fmozan,fmoc,fmomas,iactfg)
      use mx_limits, only: mxatm
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      LOGICAL GOPARR,DSKWRK,MASWRK
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON /MASSES/ ZMASS(MXATM)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /FMOINF/ NFG,NLAYER,NATFMO,NBDFG,NAOTYP,NBODY,NSEGM
      dimension layfrg(*),indat(*),fmozan(*),fmoc(3,*),fmomas(*),
     *          iactfg(*)
c
c     mode=0 use B
c     mode=1 use A
c
      CALL DAWRIT(IDAF,IODA,fmozan,natfmo,493,0)
      CALL DAWRIT(IDAF,IODA,fmoc,3*natfmo,494,0)
      nat=0
      do iat=1,natfmo
c       if(iactat(iat).ne.0) then
        if(mode.eq.0.and.layfrg(indat(iat)).gt.1 .or.
     *     mode.eq.1.and.iactfg(indat(iat)).gt.0) then
          nat=nat+1
          if(nat.gt.MXATM) then
            if(maswrk) write(iw,9000) MXATM
            call abrt
          endif
          zan(nat)=fmozan(iat)
          do j=1,3
            c(j,nat)=fmoc(j,iat)
          enddo
          zmass(nat)=fmomas(iat)
        endif
      enddo
c
      RETURN
 9000 format(/1x,'Too many atoms: ',I7, 
     *           '. Recompile with a larger MXATM.',/)
      END
*MODULE FMOGRD  *DECK unpackat
      SUBROUTINE unpackat(mode,layfrg,indat,fmozan,fmoc,iactfg)
      use mx_limits, only: mxatm
C>
C>     @brief restore atom for FMO/FD 
C>
C>     @details Restore atom for FMO/FD. 
C>
C>     @author Dmitri Fedorov 
C>
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON /FMOINF/ NFG,NLAYER,NATFMO,NBDFG,NAOTYP,NBODY,NSEGM
      dimension layfrg(*),indat(*),fmozan(*),fmoc(3,*),iactfg(*)
c
c     mode=0 use B
c     mode=1 use A
c
      CALL daread(IDAF,IODA,fmozan,natfmo,493,0)
      CALL daread(IDAF,IODA,fmoc,3*natfmo,494,0)
      nat0=nat
      nat=0
      do iat=1,natfmo
c       if(iactat(iat).ne.0) then
        if(mode.eq.0.and.layfrg(indat(iat)).gt.1 .or.
     *     mode.eq.1.and.iactfg(indat(iat)).gt.0) then
          nat=nat+1
          do j=1,3
            fmoc(j,iat)=c(j,nat)
          enddo
        endif
      enddo
      if(nat.ne.nat0) call abrtx("Atom mismatch in unpackat")
c
      RETURN
      END
*MODULE FMOGRD  *DECK packprp
C>
C>     @brief Store properties for FMO/FD 
C>
C>     @details Store properties for FMO/FD. 
C>
C>     @author Dmitri Fedorov 
C>
      SUBROUTINE packprp(mode,layfrg,indat,fmode,fmomas,iactfg,nat)
      use mx_limits, only: mxatm
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      LOGICAL GOPARR,DSKWRK,MASWRK
      COMMON /FUNCT / EHF,EG(3,MXATM)
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON /MASSES/ ZMASS(MXATM)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /FMOINF/ NFG,NLAYER,NATFMO,NBDFG,NAOTYP,NBODY,NSEGM
      dimension layfrg(*),indat(*),fmode(3,*),fmomas(*),iactfg(*)
c
c     mode=0 use B
c     mode=1 use A
c
      nat=0
      do iat=1,natfmo
        if(mode.eq.0.and.layfrg(indat(iat)).gt.1 .or.
     *     mode.eq.1.and.iactfg(indat(iat)).gt.0) then
          nat=nat+1
          if(nat.gt.MXATM) then
            if(maswrk) write(iw,9000) MXATM
            call abrt
          endif
          zmass(nat)=fmomas(iat)
          do j=1,3
            eg(j,nat)=fmode(j,iat)
          enddo
c         write(6,*) 'wwwEG',nat,(eg(j,nat),j=1,3)
        endif
      enddo
c
      RETURN
 9000 format(/1x,'Too many atoms: ',I7,
     *           '. Recompile with a larger MXATM.',/)
      END
*MODULE FMOGRD  *DECK packfat
      SUBROUTINE packfat(ifact,mapaf)
      use mx_limits, only: mxatm
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      LOGICAL GOPARR,DSKWRK,MASWRK
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      dimension ifact(*),mapaf(*)
c
c     create map(full->active)
c
c     write(6,*) 'wwworig',(ifact(i),i=1,nat)
      call setmapaf(mapaf)
c
      do iat=1,nat
        if(ifact(iat).ne.0) then
          jat=mapaf(abs(ifact(iat))) 
          if(jat.eq.0) then
            if(maswrk) write(iw,*) 'Active atom',ifact(iat),
     *                             ' is in layer 1.'
            call abrt
          endif
          if(ifact(iat).lt.0) jat=-jat
          ifact(iat)=jat
        endif
      enddo
c     write(6,*) 'wwwconv',(ifact(i),i=1,nat)
c
      RETURN
      END
*MODULE FMOGRD  *DECK setmapaf
      SUBROUTINE setmapaf(mapaf)
      use mx_limits, only: mxatm
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON /FRZCRT/ IFCART(3*MXATM),NFCART,MAPFA(MXATM)
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /FMOINF/ NFG,NLAYER,NATFMO,NBDFG,NAOTYP,NBODY,NSEGM
      dimension mapaf(natfmo)
c
      call viclr(mapaf,1,natfmo)
      do iat=1,nat
        mapaf(mapfa(iat))=iat
      enddo
c     write(6,*) 'list of mapaf',(mapaf(i),i=1,natfmo)
      RETURN
      END
*MODULE FMOGRD  *DECK outxyz
      SUBROUTINE outxyz(nserch,title,natfmo,fmozan,fmoc,izbas,nactat,
     *                  nactfg,fmosym,nfmosym)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (units=0.52917724924D+00)
      character*4 atm
      character*8 cdum
      character*80 title,fmosym
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      dimension fmozan(*),fmoc(3,*),izbas(*)
c
c     comb(1:nfmosym)=fmosym(1:nfmosym)
c     comb(nfmosym+1:nfmosym+1)=' '
c     comb(nfmosym+2:80)=title(1:80-nfmosym-1)
c     Titel is now a copy of fmosym
      write(iw,9000) nserch 
c     write(iw,9010) natfmo,nserch,nactat,nactfg,comb(1:55)
      write(iw,9010) natfmo,nserch,nactat,nactfg,fmosym(1:nfmosym)
      cdum=' '
      do i=1,natfmo
        ian=int(fmoZAN(i)+0.5D-02)
        call zsymnum(cdum,ATM,ian)
        izi=izbas(i)
        if(izi.eq.1) then
          write(iw,9020) i,    atm,(fmoc(j,i)*UNITS,j=1,3) 
        else
          write(iw,9030) i,izi,atm,(fmoc(j,i)*UNITS,j=1,3) 
        endif
      enddo
      write(iw,9100)  
      RETURN
 9000 format(/1x,'Coordinates of all atoms for NSEARCH=',I5)
 9010 format(1x,'$FMOXYZ',/1x,I6,/,'NS=',I4,' NA=',i4,' NF=',i4,1x,a)
 9020 format(1x,i6,' ',1x,1x,a4,3F20.12)
 9030 format(1x,i6,'.',i1,1x,a4,3F20.12)
 9100 format(1x,'$end',/)
      END
*MODULE FMOGRD  *DECK checkremd
C>
C>     @brief Check if doing REMD. 
C>
C>     @details Read the input to see if doing REMD.
C>
C>     @author Dmitri Fedorov 
C>
      SUBROUTINE checkremd(runtyp,doremd) 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      parameter (NSTR=80)
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      character*99 STR
      logical doremd,GOPARR,DSKWRK,MASWRK,inmd
      DATA AMD/8HMD      /
c
c     God forbid from writing subroutines like this.
c
      doremd=.false.
      inmd=.false.
c
      if(RUNTYP.eq.AMD) then
        CALL SEQREW(IR)
        if(maswrk) then
  100     CONTINUE
          READ(IR,9000,END=600,ERR=600) STR
          CALL UPRCAS(STR,NSTR)
c         write(6,*) 'Processing',str
          if(str(1:4).EQ.' $MD') inmd=.true.
          if(inmd) then
            do i=1,NSTR-4
              if(str(i:i+5).eq.'MREMD=') then
                read(str(i+6:NSTR),*,END=600,ERR=600) mremd
                doremd=mremd.ne.0
                goto 600
c               Those who like to add spaces as in "MREMD = 1"
c               will not be happy. 
              else if(str(i:i+4).EQ.' $END') then
                inmd=.false.
                goto 600
              endif
            enddo
          endif
          GOTO 100
        endif
  600   continue
        if(goparr) then
          iremd=0
          if(doremd) iremd=1
          CALL DDI_BCAST(2410,'I',iremd,1,MASTER)
          doremd=iremd.ne.0 
        endif
      endif
c     write(6,*) 'wwwremd',doremd
      return
 9000 format(A80)
      END
*MODULE FMOGRD  *DECK lumpesd 
      SUBROUTINE lumpesd(nder,loadhf,nfg2,nes,nscf,ndimscf,nesdim3,
     *                   nfg2d,iifg,edisnrg,gcorrel,dodc,reducee,
     *                   ESDIMTMP,edim,loadd,savemem2,savememr,
     *                   edimsum,irij,nrij,wrk,wrk1,ichfg,mdout)
      use mx_limits, only: mxrt
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (zero=0.0D+00,one=1.0D+00)
      logical ISGDDI,PAROUT,INITGDDI,wasgddi,GOPARR,DSKWRK,MASWRK,myjob,
     *        gcorrel,dodc,reducee,DFTBFL,SCC,SRSCC,DFTB3,DAMPXH,LCDFTB,
     *        MLGDDI,savemem2,savememr,mdout,dotrunces,dotruncesd,doall
      dimension ESDIMTMP(*),edim(nfg2d,*),loadd(*),
     *         edimsum(*),irij(maxrij,*),nrij(*),wrk(*),wrk1(*),ichfg(*)
      common /cenrad/ rcutes,lcradfg,lichfgg,dotrunces,dotruncesd
      COMMON /DFTB  / DFTBFL,SCC,SRSCC,DFTB3,DAMPXH,LCDFTB
      COMMON /ENRGYS/ ENUCR,EELCT,ETOT,SZ,SZZ,ECORE,ESCF,EERD,E1,E2,
     *                VEN,VEE,EPOT,EKIN,ESTATE(MXRT),STATN,EDFT(3),EDISP
      COMMON /FMCOM / X(1)
      COMMON /GDDI/   ISCOPE,NGROUPS,MYGROUP,MEGLOB,npglob,nnglob,JBTYP,
     *                ISGDDI,PAROUT,INITGDDI,wasgddi,MLGDDI,NSUBGR,
     *                MeUniv,NPUniv,numdlb,myworld,nworlds,mogddi
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      common /fmodim/ maxbnd,maxknd,maxcbs,maxcao,maxvec,maxl1,maxnat,
     *                maxabd,maxbas,maxbbd,maxlmo,maxslo,maxabd2,maxrij
      COMMON /FMODTB/ LGAMMA,LZREF,LISPE,LSHIFTG,LSHIFTCT,
     *                LGAMMAD,MODESD,MODGAMMA,LDCTMP,LGAMMA3,LGAMMA3D,
     *                LNBOND,LCTMUL,LESPDFTB,LESPASC,NDFTBRST,LDFTBRST
      COMMON /FMOINF/ NFG,NLAYER,NATFMO,NBDFG,NAOTYP,NBODY,NSEGM
      common /fmoopt/ espsca(9),RESPAP(2),rESPPC(2),rESDIM,restri(4),
     *                rcorsd,respct,convfg,cnvdmp,coroff,rflmo(4),
     *                orshft,orshft2,cnvafo,ascreen(4),IXESP,mxitfg,
     *                nguess,NBSSE,modorb,modpar,irststp,irstlay,nprfmo,
     *                nfmopal,modprp,maxl1c,ipieda,modgrd,modesp,ivmul,
     *                modlmo,nopden,mofock,modfd,modfmm,ncentm,ndualb,
     *                ngab,modpan
c
c     Lump ESD for DFTB (add them together).
c
CDFTB CALL SYSTEM_CLOCK(ICLOCK_COUNTS_3,CLCK_RATE)
CDFTB CALL CPU_TIME(T3)
C     IIFG > NFG MEANS THAT THERE ARE NO ESD2 CALCULATIONS
CDFTB CALL SYSTEM_CLOCK(ICLOCK_COUNTS_4,CLCK_RATE)
CDFTB CALL CPU_TIME(T4)
      IF (MASWRK.AND.MDOUT) THEN
        WRITE(IW,FMT='('' ...... END OF SCF2 STEP ......'')')
        CALL TIMIT(1)
        WRITE (IW,*) ' '
      END IF
      doall=IAND(MODESD,16).ne.0
c     doall selects a full span of IFG/JFG with SAVEMEMR to avoid
c     calling loaddes at the little extra cost of looping over SCF dimers.
c     nprfmos=nprfmo
      IF (IAND(MODESD,8).EQ.0) THEN
        if (iand(nprfmo,2).ne.2) then
          if (iand(nprfmo,1).ne.1) then
            nprfmo = nprfmo + 2
          else
            nprfmo = nprfmo + 1
          end if
        end if
      END IF
      IF (MASWRK.AND.MDOUT) WRITE (IW,*) 'LUMP ES-DIM STEP FOR DFTB'
c     IF (IAND(MODESD,6).NE.0.AND.IIFG.LE.NFG) THEN
      IF (DODC) EDISP=ZERO
      ndone=0
      nskept=0
      nstep1  = 1
      IF (IAND(MODESD,6).NE.0) then
        if(IIFG.GT.NFG) goto 1001 
        if (loadhf.ne.1) then
          write (iw,*) "loadhf.ne.1 --> abrt"
          call abrt
        end if
C       CALL CPU_TIME(T0)
        call tribrk(nscf+1,-1,iifg,jjfg)
        lijfg = nscf
        NES = NFG2 - NSCF
C       global synchronization before ESD calculation.
c       IF (DODC) EDISP = ZERO
 1001   continue
        if(isgddi) then
c         call GDDICOUNT( 1,lgroup,myjob)
          call GDDICOUNT(-1,lgroup,myjob)
        end if
        if(IIFG.GT.NFG) goto 1002 
        IF (isgddi) THEN
          GI = 1.0D-01/NGROUPS
        ELSE
          GI = 1.0D-01
        END IF
        NESDDI = INT(NES*GI)+1
        ESDDII = ONE / NESDDI
        IF (IAND(MODESD,2).NE.0) THEN
          NSTART1 = IIFG
          NLAST1  = NFG
          NSTART2 = JJFG
          NLAST2  = IIFG-1
          if(doall) then
            if(.not.savememr) call abrtx("Failure in lumpesd")
c           savememr is at present the only way coded; it is possible
c           to implement MODESD=16 without it but probably not useful?!
c           The order of IFG is reversed (nstep1=-1) for parallel
c           efficiency (to do the largest IFG=NFG first).
c           There is no point at all to reverse JFG.
            NSTART1 = nfg 
            NLAST1  = 1
            nstep1  = -1
            NSTART2 = 1
            NLAST2  = nfg
c           write(6,*) 'wwwesd 16',nstep1
          endif
        ELSE IF (IAND(MODESD,4).NE.0) THEN
          NSTART1 = 1
          NLAST1  = INT(NES*ESDDII) + 1
          NSTART2 = 1
          NLAST2  = NESDDI
        END IF
C       WRITE (IW,*) "NUMBER OF PARALLEL LOOP =" , NLAST1-NSTART1+1
        inextscf=0
        do iiifg= NSTART1, NLAST1, nstep1
c         dynamic load balancing for the first index
          IF (IAND(MODESD,2).NE.0) THEN
            NLAST2 = iiifg-1
c           The same NLAST2 is good for MODESD=16. 
          ELSE IF (IAND(MODESD,4).NE.0) THEN
            NLAST2 = NESDDI
            IF (IIIFG.EQ.NLAST1) NLAST2 = NES-NESDDI*(NLAST1-1)
          END IF
          if(isgddi) then
            call GDDICOUNT(0,lgroup,myjob)
            if(.not.myjob) then
              lijfg = lijfg + (NLAST2 - NSTART2 + 1)
              goto 700
            end if
          endif
          if(doall) then
            indrij=1
            numrij=nrij(iiifg)
            inextscf=0
            ijfg=0
            if(indrij.le.numrij) inextscf=abs(irij(indrij,iiifg))
c           write(6,*) 'wwwij',iiifg,numrij
c           write(6,*) 'wwwirij',(irij(iin,iiifg),iin=1,numrij)
          endif
          do jjjfg= NSTART2, NLAST2
            lijfg=lijfg+1
c           static load balancing for the second index
            if(mod(jjjfg,nproc).eq.me) then
              if(savememr) then
                if(doall) then
                  if(jjjfg.eq.inextscf) then
c                   if iiifg,jjjfg is an SCF dimer
c                   Assume that JFG in irij are stored in ascending order!
                    indrij=indrij+1
                   if(indrij.le.numrij) inextscf=abs(irij(indrij,iiifg))
c                   Skip SCF dimers.
                    cycle
                  endif
c                 If we got here, that means iiifg,jjjfg is an ES dimer
                  ifg=iiifg
                  jfg=jjjfg
c                 ijfg is not set!
c                 write(6,*) 'wwwt',ifg,jfg,ijfg
                else
                  call loaddes(nrij,irij,lijfg-nscf,ifg,jfg,ijfg)
c                 write(6,*) 'wwwt',ifg,jfg,ijfg
                endif
              else
                ijfg=loadd(lijfg)
                call tribrk(ijfg,-1,ifg,jfg)
              endif
c             esd2=emon(ifg,memon)+emon(jfg,memon)
c
              NESDIM3 = NESDIM3 + 1
              if(dotruncesd) then
                if(ichfg(ifg).eq.0.and.ichfg(jfg).eq.0) then
                  lcradfgi=lcradfg+(ifg-1)*4
                  lcradfgj=lcradfg+(jfg-1)*4
                  if(fmodest(x(lcradfgi),x(lcradfgi+3),x(lcradfgj),
     *                       x(lcradfgj+3)).gt.rcutes) then
c                   write(6,*) 'Skipping',ifg,lfg
                    nskept=nskept+1
                    cycle
                  endif
                endif
              endif
              ndone=ndone+1
c
              CALL DFTB_ESDIM3(IFG,JFG,NDER,ESDIMTMP,WRK1,
     *          WRK1(MAXNAT*2+1),WRK,MDOUT)
              IF (DODC)
     *          CALL DFTB_ESDIM_DC_HELP(IFG,JFG,NDER,WRK)
c             if(savemem2) then
c               edimsum(1)=edimsum(1)+ESD2 
c               if(dodc) edimsum(3)=edimsum(3)+ESD2 
c             else
c               EDIM(IJFG,1) = ESD2
c               if(dodc.and..NOT.REDUCEE) EDIM(IJFG,3) = EDIM(IJFG,1)
c             endif
            endif
          end do
  700     jjfg = 1
          NSTART2 = 1
        continue
        end do
 1002   continue
        if(isgddi) call GDDICOUNT( 1,lgroup,myjob)
        if(goparr) call DDI_gsumi(2423,NESDIM3,1)
        if(maswrk) write(iw,9000) ndone,nskept,nfg2-nscf
C       nprfmo=nprfmos
c     else
c       IF (DODC) EDISP=ZERO
      end if
CDFTB CALL SYSTEM_CLOCK(ICLOCK_COUNTS_5,CLCK_RATE)
CDFTB CALL CPU_TIME(T5)
      CALL DFTB_ESDIM3(0,NFG,NDER,ESDIMTMP,WRK,WRK,WRK,MDOUT)
      if(savemem2) then
        IF (NESDIM3.EQ.0.and.DODC) EDISP = ZERO
      else
        IF (NESDIM3.EQ.0) THEN
          if(savememr) then
c           Here, find the last index to add the remainder.
c           It does not really matter which ES dimer it is,
c           do doall just calls loaddes once too.
            call loaddes(nrij,irij,nfg2-nscf,ifg,jfg,ijfg)
c           write(6,*) 'wwwaaa',ifg,jfg,ijfg
          else
            IJFG = LOADD(NFG2)
          endif
          IF (DODC) EDISP = ZERO
        END IF
      endif
      if(doall.and.ijfg.eq.0) then
c       Need to set up ijfg of the last ES dimer
        call loaddes(nrij,irij,nfg2-nscf,ifg,jfg,ijfg)
c       write(6,*) 'wwwaaa',ifg,jfg,ijfg
      endif
      IF (SCC) then
        if(savemem2) then
          edimsum(1)=edimsum(1)+ESDIMTMP(1)
          if(dodc) edimsum(3)=edimsum(3)+ESDIMTMP(1)
        else
          EDIM(IJFG,1) = EDIM(IJFG,1) + ESDIMTMP(1)
          if((gcorrel.or.dodc).and..NOT.REDUCEE)
     *      EDIM(IJFG,3) = EDIM(IJFG,3)+ESDIMTMP(1)
        endif
      endif
c     maybe separate loop is more efficient
      if(gcorrel.or.dodc) THEN
        edisnrg = edisp
        if(savemem2) then
          edimsum(1)=edimsum(1)+EDISNRG
        else
          EDIM(IJFG,1) = EDIM(IJFG,1) + EDISNRG
        endif
c       IF (.NOT.REDUCEE) THEN
c         if(savemem2) then
c           Is this necessary? These are not proper values?
c           DO LIJFG = NSCF+1, NFG2
c             ijfg=loadd(lijfg)
c             EDIM(IJFG,3) = EDIM(IJFG,1)
c           END DO
c           edimsum(3)=edimsum(3)-edisnrg
c         else
c           if(savememr) then
c           DO LIJFG = NSCF+1, NFG2
c             call loaddes(nrij,irij,LIJFG-nscf,ifg,jfg,ijfg) 
c             EDIM(IJFG,3) = EDIM(IJFG,1)
c           END DO
c           else
c           DO LIJFG = NSCF+1, NFG2
c             ijfg=loadd(lijfg)
c             EDIM(IJFG,3) = EDIM(IJFG,1)
c           END DO
c           endif
c           edim(ijfg,3)=edim(ijfg,3)-edisnrg
c         endif
c       END IF
      END IF
CDFTB CALL CPU_TIME(T6)
CDFTB CALL SYSTEM_CLOCK(ICLOCK_COUNTS_6,CLCK_RATE)
      if(maswrk.AND.MDOUT) then
        WRITE (IW,*) "NUMBERS OF SCF2 AND ESD FOR GROUP ", MYGROUP
        WRITE (IW,*) "NDIMSCF = ", NDIMSCF
        WRITE (IW,*) "NESDIM  = ", NESDIM3
        WRITE (IW,*) " "
CDFTB   WRITE (IW,*) "WALL AND CPU TIMES FOR GROUP ", MYGROUP
CDFTB   WRITE (IW,*) "PREPARATION OF EDIMER SUBROUTINE"
CDFTB   WRITE (IW,'(F10.5,3X,F10.5)') 
CDFTB* (ICLOCK_COUNTS_2-ICLOCK_COUNTS_1) / REAL(CLCK_RATE), T2-T1

CDFTB   WRITE (IW,*) "SCF2 STEP"
CDFTB   WRITE (IW,'(F10.5,3X,F10.5)') 
CDFTB* (ICLOCK_COUNTS_3-ICLOCK_COUNTS_2) / REAL(CLCK_RATE), T3-T2

CDFTB   WRITE (IW,*) "ESD MAIN LOOP"
CDFTB   WRITE (IW,'(F10.5,3X,F10.5)') 
CDFTB* (ICLOCK_COUNTS_5-ICLOCK_COUNTS_4) / REAL(CLCK_RATE), T5-T4

CDFTB   WRITE (IW,*) "ESD SUMMARY STEP"
CDFTB   WRITE (IW,'(F10.5,3X,F10.5)') 
CDFTB* (ICLOCK_COUNTS_6-ICLOCK_COUNTS_5) / REAL(CLCK_RATE), T6-T5

CDFTB   WRITE (IW,*) "CPU TIME FOR PRE-LOOP", T1-T0, MYGROUP
CDFTB   WRITE (IW,*) "CPU TIME FOR SUMMARY ", T2-T1, MYGROUP
CDFTB   WRITE (IW,'(2X,"TOTAL ESD2 ENERGY= ",F15.10)') ESDIMTMP(1)
      END IF
      return
 9000 format(1x,'LUMPESD: did',I16,' skept',I16,' total',I20)
      END
*MODULE FMOGRD  *DECK latgconstr
C>
C>     @brief Apply lattice constraints.
C>
C>     @details Use constrains.
C>
C>     @author Dmitri Fedorov
C>
      SUBROUTINE latgconstr(latopt,latopa,ldim,gelat)
      IMPLICIT NONE
C     Declarations of arguments
      INTEGER latopt,latopa(6),ldim
      DOUBLE PRECISION gelat(3,3)
C     Declarations of common blocks
      INTEGER IR,IW,IP,IS,IPK,IDAF,NAV,IODA
      INTEGER ME,MASTER,NPROC,IBTYP,IPTIM
      LOGICAL GOPARR,DSKWRK,MASWRK
C     Other declarations
      INTEGER i,j
C
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
c
c     Transform lattice gradient so that the constraints are applied.
c     Right now the only option is for a rectangular cell to remain
c     rectangular.
c
      if(latopt.ne.2) return 
      if(maswrk) write(iw,*) 'Applying rectangular lattice constraints.'
c     write(6,*) 'wwwpreg',gelat
      do i=1,ldim
        do j=1,3
          if(i.ne.j.or.latopa(i).ne.0) gelat(i,j)=0
        enddo
      enddo
c     write(6,*) 'wwwpostg',gelat
      return
      END
C
C-----------------------------------------------------------------------
*MODULE FMOGRD  *DECK getcellpar 
C>
C>     @brief   Get cell parameters.
C>
C>     @details Compute cell parameters.
C>
C>     @author  Dmitri Fedorov
C>
      SUBROUTINE getcellpar(ldim,box,cellpar)
      IMPLICIT NONE
C     Declarations of arguments
      INTEGER ldim
      DOUBLE PRECISION box(3,3),cellpar(6)
C     Other declarations
      DOUBLE PRECISION a,cosa,cosb,cosg
      INTEGER i
c
      call vclr(cellpar,1,6)
c     a, b, and c
      do i=1,ldim
        a=box(i,1)*box(i,1)+box(i,2)*box(i,2)+box(i,3)*box(i,3)
        cellpar(i)=sqrt(a)
      enddo
c     alpha, beta, gamma
      if(ldim.gt.2) then
        cosa=box(2,1)*box(3,1)+box(2,2)*box(3,2)+box(2,3)*box(3,3)
        cellpar(4)=acos(cosa/(cellpar(2)*cellpar(3)))
        cosb=box(1,1)*box(3,1)+box(1,2)*box(3,2)+box(1,3)*box(3,3)
        cellpar(5)=acos(cosb/(cellpar(1)*cellpar(3)))
      endif
      if(ldim.gt.1) then
        cosg=box(1,1)*box(2,1)+box(1,2)*box(2,2)+box(1,3)*box(2,3)
        cellpar(6)=acos(cosg/(cellpar(1)*cellpar(2)))
      endif
      return
      END
*MODULE FMOGRD  *DECK restpcm 
C>
C>     @brief Process restart data for PCM.
C>
C>     @details Do restarts.
C>
C>     @author Dmitri Fedorov
C>
      SUBROUTINE restpcm(buf,mode,instance,ires)
      IMPLICIT NONE
C     Declarations of arguments
      DOUBLE PRECISION buf(*)
      INTEGER mode,instance,ires
C     Declarations of common blocks
      DOUBLE PRECISION X
      INTEGER mxsp,mxts,mempcm1,mempcm2,NTS
      INTEGER LXYZRE,LSSFE,LLIST,LALPCM,LRINPCM,LINAPCM,LINFPCM
      INTEGER LINIPCM,LMEPCM,LAXYZCT,LXYZCT2,LNVERT,LQSN,LQSE
      INTEGER LQFS,LQIND,LISPHE,LVAD,LQOR,LVECMUL,LAIPRJ,LFIPRJ
      INTEGER mPEL,lpcmcdr,LDAI,LIDDAI,LQSND,LQSED,LEPSPCM
      INTEGER mxtssav,ntssav,lisphsav,lqpcmsav
C     Other declarations
      INTEGER lqpcmsavi
C
      COMMON /FMCOM / X(1)
      COMMON /PCMDIM/ mxsp,mxts,mempcm1,mempcm2,NTS
      COMMON /PCMPNT/ LXYZRE,LSSFE,LLIST,LALPCM,LRINPCM,LINAPCM,LINFPCM,
     *                LINIPCM,LMEPCM,LAXYZCT,LXYZCT2,LNVERT,LQSN,LQSE,
     *                LQFS,LQIND,LISPHE,LVAD,LQOR,LVECMUL,LAIPRJ,LFIPRJ,
     *                mPEL,lpcmcdr,LDAI,LIDDAI,LQSND,LQSED,LEPSPCM
      COMMON /PCMSAV/ mxtssav,ntssav(3),lisphsav,lqpcmsav
c
c     Process tessera data for restarts in PCM.
c     mode=0: save
c     mode=1: restore
c     instance can be 1, 2 or 3 (SCC, final, SCZV)
c     ISPHE is the same, but still copy.
c     Only ICOMP=0 is supported.
c     ires is set to 1 if restart was successful and 0 otherwise
c     BUF(NTS) is an input argument for mode=0 (data to be saved)
c     BUF(NTS) is an output argument for mode=1 (restored data)
c
      ires=0
      if(mxtssav.le.0) return
c
      if(instance.lt.1.or.instance.gt.3) call abrtx("PCM restart error")
      if(NTS.gt.mxtssav) call abrtx("Tessera mismatch for PCM restart")
c
      lqpcmsavi=lqpcmsav+(instance-1)*mxtssav
      if(mode.eq.0) then
c       write(6,*) 'Saving PCM data',instance,nts
        ntssav(instance)=nts
        call icopy(nts,x(LISPHE),1,x(lisphsav),1)
        call dcopy(nts,BUF,1,x(lqpcmsavi),1)
      else if(ntssav(instance).ne.0) then
c       write(6,*) 'Restoring PCM data',instance,mxtssav
        call restpcmdat(ntssav(instance),x(lisphsav),x(lqpcmsavi),nts,
     *                  x(LISPHE),BUF)
c       if(instance.eq.1) then
c         Set up PCM charges to be used in SCF
c         call dcopy(nts,BUF,1,x(LQSE),1)
c       endif
c       Nothing is done for the first call when no restart data are available.
c       The calling subroutine should treat this case. 
        ires=1
      endif
      return
      END
*MODULE FMOGRD  *DECK restpcmdat
C>
C>     @brief Map charges for PCM.
C>
C>     @details Do restarts.
C>
C>     @author Dmitri Fedorov
C>
      SUBROUTINE restpcmdat(ntssav,isphsav,qpcmsav,nts,isphe,qpcm)
      IMPLICIT NONE
C     Declarations of arguments
      INTEGER ntssav,isphsav(*)
      DOUBLE PRECISION qpcmsav(*)
      INTEGER nts,isphe(*)
      DOUBLE PRECISION qpcm(*)
C     Other declarations
      DOUBLE PRECISION, PARAMETER :: ERR=1.0D-10
      INTEGER i,iat,j,jat,ncopy
      DOUBLE PRECISION qtot,qtots
      DOUBLE PRECISION ddot
c
c     Map charges from saved to current arrays in PCM.
c     This is tricky because tessera maps ISPHSAV/ISPHE often are different.
c     This subroutine assumes that tesserae are stored in the
c     increasing order of atoms (first atom 1, then 2 etc).
c     Some atoms may not have any tesserae.
c     Spheres are loosely called "atoms" here.
c     NTS=0 and NTSSAV=0 are not possible, but still handle them.
c     Sanity check: charge conservation!
c
      if(nts.le.0) return
      call vclr(qpcm,1,nts)
      if(ntssav.le.0) return
c
      jat=-1
      j=0
c     jatp=-1
      ncopy=0
      do i=1,ntssav
        iat=isphsav(i)
   90   continue
        if(j.lt.nts) then
          jat=isphe(j+1)
c         checking for atoms is the effort to match tesserae of each
c         atom in the best possible way.
          if(iat.eq.jat.or.j.eq.0) then
c           If the same atom, just copy.
            j=j+1
            qpcm(j)=qpcmsav(i)
            ncopy=ncopy+1
c           Copying is the most direct correspondence; at least the
c           atom is preserved; for a direct tesserae-tesserae match one has
c           to match coordinates.
c           jatp=jat
          else
c           If a different atom, there are two cases.
            if(jat.gt.iat) then
c             The new list is shorter: add extra old charges to the last
c             new charge. Do not move up the new list.
c             The idea is to pile the charge of a lost tessera to the last
c             matching tessera. 
c             The piling of charges is done for the charge conservation;
c             one can argue that they should be smeared rather than piled...
c             Pile even to other atoms.
c             old: 1,2*,
c             new: 1@,3
c             The charge of * is piled to @.
c             if(iat.eq.jatp) qpcm(j)=qpcm(j)+qpcmsav(i)
              qpcm(j)=qpcm(j)+qpcmsav(i)
            else
c             The new list is longer: pad zeroes.
c             old: 1,1,2@,2
c             new: 1,1,1*,1*,2
c             Pad zeros to *.
100           continue
              j=j+1
              qpcm(j)=0
              if(j.lt.nts) then
                if(isphe(j+1).lt.iat) goto 100
              endif
c             After padding (*), now process @.
              goto 90
            endif
          endif
        else
c         The new list is done, but the old list still lasts
c         Pile new charges to the last tesserae.
c         old: 1,1,2,2*,3*,3*
c         new: 1,1,2@
c         The charge of * is piled to @.
c         if(iat.eq.jat) qpcm(j)=qpcm(j)+qpcmsav(i)
          if(j.eq.0) j=1
          qpcm(j)=qpcm(j)+qpcmsav(i)
        endif
      enddo
c     Padding zeros to * is now done in vclr above.
c     old: 1,2,3,4,4
c     new: 1,2,3,4,4,4*,5*,5*
c
c     Now, the new list should have all values set!
      qtots=ddot(ntssav,qpcmsav,1,1.0D+00,0)
      qtot=ddot(nts,qpcm,1,1.0D+00,0)
c     write(6,*) 'wwwQs:',qtots,qtot,ntssav,nts,ncopy
      if(abs(qtots-qtot).gt.ERR) call abrtx("PCM restart mismatch")
c     Not matching means the algorithm here is incorrect.
      return
      END
