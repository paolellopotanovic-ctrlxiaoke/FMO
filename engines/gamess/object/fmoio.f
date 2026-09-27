C  8 Aug 21 - DD Interface to RI-CC added in fmoccp
C 17 Oct 19 - HN,YN,VQV,DSK,VM,DGF - changes for FMO 5.4
C  6 Jun 18 - HN,YN,DK,DGF - changes for FMO 5.3
C  1 Apr 16 - TN,HN,DGF - changes for FMO 5.2
c 22 Oct 14 - HN,DGF,KRB,NM,SRP - changes for FMO 5.1
C 21 May 13 - DGF,HN,TN - changes for FMO 5.0
C 19 Oct 12 - MWS - synchronize FRGINF common
C 31 Jul 12 - DGF,CHC - LAST CHANGES FOR FMO 4.3
C 24 Jul 12 - HN,DGF - code update to finish FMO 4.3
C 21 JUN 12 - DGF - changes for FMO 4.3
C 23 MAR 12 - DGF,CHC - code update to finish FMO 4.2
C 28 DEC 11 - DGF,CS  - FMO 4.2 and EFMO changes
C  2 Nov 11 - DGF - keyword VDWRAD should be read as F.P. type
C 15 Apr 11 - DGF,TN - frozen domain, EFP and PCM related changes
C  1 Oct 10 - CS  - introduced EFMO changes
C 11 Aug 10 - DGF,TN - changes for FMO 4.0 
C 25 Mar 10 - DGF - tickles to complete FMO 3.3 release
C 14 Oct 09 - DGF - changes for FMO 3.3
C 22 May 09 - DGF - print correct FMO/MCP Mulliken charges
C 23 Jan 09 - DGF - output changes for FMO 3.2
C 15 Dec 08 - DGF,MC - various changes for FMO 3.2 release
C 18 Jul 08 - DGF - FMOPROP: all nodes should compute 3-body terms
C 11 Apr 08 - MWS - synchronize DFGRID common block
C  4 Mar 08 - DGF - FMOMINP: fix the RITRIM check
C 28 Aug 07 - DGF - small printing changes
C 20 Aug 07 - DGF - various changes for FMO 3.1 release
c 12 Jul 07 - MC  - add FMO-TDDFT arguments
C 24 Mar 07 - MWS - pad FRGINF common block
C 22 Dec 06 - TN  - synchronise EFPFMO common
C  8 Nov 06 - DGF - various changes for FMO 3.0 release
C  6 Nov 06 - MWS - adjust wavefunction and GDDI common block
C 22 Feb 06 - TN  - FMOPROP: include EFP/FMO model
C 21 Nov 05 - DGF - various changes for FMO 2.1 release
C 19 Sep 05 - IA  - synchronize FRGINF common
C  6 Jul 05 - DGF - FMOCCP: update print routine for CCSD(TQ)
C  1 Jun 05 - DGF - fixes for the 2nd release
c 31 Jan 05 - DGF,KK - add I/O module for the FMO method
c
C*MODULE fmoio   *DECK fmominp
C>
C>    @brief Reads in input parameters from $fmo
C>
C>    @author Dmitri Fedorov
C>
C>    @date October, 2012 - Colleen Bertoni
C>    - Added two additional parameters to the EFMO modefmo array so
C>      the user can turn dispersion, charge transfer, and/or
C>      exchange repulsion on or off.
C>    - Modified EFMO check code to read the two additional spots
C>      in the EFMO modefmo array and set the appropriate flags.
C>    @date Jan, 2017 - C. Bertoni
C>    - Added a variable iefmo_agrad to the input to flag if
C>      it's an analytical gradient or not
C>    @date June, 2020 - Anastasia Gunina
C>    - Added a keyword to bypass dispersion error with AFO
C>      (+2 in the third position of modefm; imodefd here)
C>    @date Dec, 2023 - Peng Xu
C>    - add more options for IMODEFD for disp7 and iso disp8
C>
C>    @param itask :
C>    @param nder :
C>    @param ifgfmo0 :
C>    @param ichfg :
C>    @param mulfg :
C>    @param frgnam :
C>    @param layfrg :
C>    @param indat :
C>    @param scffrg :
C>    @param fmoscf :
C>    @param fmoci :
C>    @param fmodft :
C>    @param fmocc :
C>    @param mpnfmo :
C>    @param fmotd :
C>    @param nacut :
C>    @param modmol :
C>    @param molfrg :
C>    @param nprfrg :
C>    @param lbody :
C>    @param gcorrel :
C>    @param fmoq :
C>    @param iexcit :
C>    @param exfid :
C>    @param modcha :
C>    @param atchrg :
C>    @param iactfg :
C>    @param nactfg :
C>
      SUBROUTINE fmominp(itask,nder,ifgfmo0,ichfg,mulfg,frgnam,layfrg,
     *                  indat,scffrg,fmoscf,fmoci,fmodft,fmocc,mpnfmo,
     *                  fmotd,nacut,modmol,molfrg,nprfrg,lbody,gcorrel,
     *                  fmoq,iexcit,exfid,modcha,atchrg,iactfg,nactfg,
     *                  rijskip,nextraind,nindatp,indatp,nmodseg,modseg,
     *                  nscfd)
      use EFP_logical
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      character*8 tstring
      logical GOPARR,DSKWRK,MASWRK,gcorrel,fmoq,
     *        DFTBFL,SCC,SRSCC,DFTB3,DAMPXH,LCDFTB
      PARAMETER (nnam=57,MaxLay=5,ten=1.0D+01)
      DIMENSION QNAM(NNAM),KQNAM(NNAM),ichfg(*),mulfg(*),frgnam(*),
     *          layfrg(*),indat(*),scffrg(*),fmoscf(*),fmoci(*),
     *          fmodft(*),fmocc(*),mpnfmo(*),fmotd(*),molfrg(*),
     *          nprfrg(*),lbody(*),iexcit(6),atchrg(*),iactfg(*),
     *          ritrim(5),indatp(*),modseg(*)
      COMMON /DFTB  / DFTBFL,SCC,SRSCC,DFTB3,DAMPXH,LCDFTB
      COMMON /EFMO  / EFMOETOT,EFMOESERG,EFMOPOLERG,EFMODISERG,
     *                EFMOREPNRG,EFMOCHTNRG,EFMOEPEN,EFMOPCMG,
     *                IEFMORUN,IMODEFP,IEFMORT,IEFMOCFRG,IEFMONFRG,
     *                IEFMODIM(2),IMODEFE,NATEFMO,IMODEFD,IMODEFER,
     *                IMODEFCT,IDIMTYP,iefmo_agrad
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON /MEXOPT/ MEXSKPGES,MEXSTATE
      COMMON /MEXPAR/ SCF1,SCF2,TGMAX,TDE,TDXMAX,TDXRMS,TGRMS,STPSZ,
     *                MULT1,MULT2,NSTEP,NRDMOS,NMOS1,NMOS2,NPRT,IMEXFG
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
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
      common /fmodim/ maxbnd,maxknd,maxcbs,maxcao,maxvec,maxl1,maxnat,
     *                maxabd,maxbas,maxbbd,maxlmo,maxslo,maxabd2,maxrij
      EQUIVALENCE (tstring,dstring)
      DATA KQNAM /1,1,1,1,1, 5,1,23,23,3, 53,1,1,1,1, 1,5,5,5,5,
     *            1,1,3,1,1, 1,1,5,3,1,    1,1,1,1,61, 5,3,1,3,43,
     *            1,1,1,1,3, 43,1,1,1,51,  1,1,21,1,1, 1,1/
      DATA QNAM/8HNFRAG   ,8HNLAYER  ,8HICHARG  ,8HMULT    ,8HINDAT   ,
     *          8HFRGNAM  ,8HLAYER   ,8HRESPAP  ,8HRESPPC  ,8HRESDIM  ,
     *          8HRITRIM  ,8HMAXKND  ,8HMAXCAO  ,8HNBSSE   ,8HMPLEVL  ,
     *          8HIXESP   ,8HDFTTYP  ,8HSCFFRG  ,8HCITYP   ,8HCCTYP   ,
     *          8HNACUT   ,8HNFRND   ,8HRCORSD  ,8HMOLFRG  ,8HNBODY   ,
     *          8HNOPFRG  ,8HMODMOL  ,8HSCFTYP  ,8HORSHFT  ,8HMODGRD  ,
     *          8HRESPCT  ,8HMAXBND  ,8HMODESP  ,8HIVMUL   ,8HIEXCIT  ,
     *          8HTDDFT   ,8HEXFID   ,8HMODCHA  ,8HATCHRG  ,8HRAFO    ,
     *          8HMODAFO  ,8HNATCHA  ,8HMAXBBD  ,8HNOPDEN  ,8HORSHF2  ,
     *          8HSCREEN  ,8HIEFMO   ,8HMODFD   ,8HIACTFG  ,8HMODEFM  ,
     *          8HMODMUL  ,8HNDUALB  ,8HMAXRIJ  ,8HIEFMOG  ,8HNGAB    ,
     *          8HINDATP  ,8HNOPSEG  /
      data fmogrp/8HFMO     /,rnone/8HNONE    /,UHF/8HUHF     /,
     *     ROHF/8HROHF    /
      dimension modefmo(5),maxrij2(2)
c
c     this subroutine reads FMO input.
c     itask=1 light version: only read dimensions and options
c     itask=2 heavy version: read everything
c     Keywords in QNAM must be arranged appropriately.
c
c     atomic population approximation 
      defval=-99.0D+00
      RESPAP(1)=0
      RESPAP(2)=0
c     point charge approximation 
      RESPPC(1)=defval
      RESPPC(2)=0
c     separated dimer approximation (energy without SCF) 
      RESDIM=defval
c     restri(1,2,3) thresholds to ignore SCF trimer contribution
c     restri(4) a threshold to ignore correlated trimer contribution
c     The condition to NOT ignore trimer SCF contribution is:
c     rmin.le.restri(1).or.rmax.le.restri(2)).and.(rmax.le.
c    *            restri(3).or.restri(3).eq.0)).or.restri(1).eq.0
c     where rmin and rmax are returned by fmodist3(ifg,jfg,kfg,rmin,rmax)
c     The condition to NOT ignore trimer correlated contribution is:
c     rmax.le.restri(4).or.restri(4).eq.0
      RESPCT=0.0D+00
c
      call dacopy(4,defval,ritrim,1)
      ritrim(5)=0
c     separated dimer approximation (energy without correlation)
      rcorsd=defval
      nbsse=0
      IXESP=0
c     bit additive
c     2 remove delta-D terms(2nd/3rd order energy correction).
c       (That is, terms of the form Tr(delta-Dij*Vij) or Tr(delta-Dijk*Vijk).
c       Useful to make energy gradient analytic by changing the 
c       energy, not the gradient!
c     4 add derivative delta-D terms (not implemented)
c     32 do not add separated dimer derivative contributions. 
c 
      NACUT=0
c     slice the molecule into fragments each having NACUT atoms,
c     if NACUT is set by the user, INDAT is ignored. All other options must be
c     set by the user.
      nfmopal=0
c     "friend option" - provide compatibility with other implementations of FMO
c     nfmopal=1 GAUSSIAN-94, round off interfragment distances to 0.1 waals. 
c     nfmopal=2 MacMolPlt support 
c
c     the B parameter in the projection operators
      orshft=1.0d+06
      orshft2=0.0d+00
c
c     n-body expansion of FMO 
      nbody=2
      if(ifgfmo0.ne.0) nbody=1
      call icopy(MaxLay,-1,0,lbody,1)
c
      modmol=0
c     Selected monomer option (molfrg), bit additive.
c     1 skip SCF mode. If set, do not perform SCF calculations of non-selected
c       n-mers otherwise, perform correlated calculations of selected and SCF
c       of all n-mers (n>1).  
c     2 intra/intermolecular (set/unset) choice. Intramolecular means
c       do n-mers within the selected group, otherwise between selected and 
c       the rest.
c     4 save memory by not allocating space for n-mer arrays (active only
c       with intramolecular mode (2) and, at present, with n=3 only).
c     "Skip SCF" is not useful with nbody=3 and the current code. Perhaps in 
c     this case one should compute all dimers and then selected trimers with at
c     most one monomer from molfrg? 
c     8 subsystem analysis (SA)
c
      modgrd=-1
c     modgrd is bit additive.
c      1 do not subtract ESP from Lagrangian. This is needed to restore the
c        old published MP2 gradient.
c      2 add ESP derivatives (T. Nagata).
c      4 read Lagrangian from data save in RHF for FMO/AFO
c      8 add Mulliken charge derivative
c     16 do not add HOP derivatives 
c     32 reserved 
c     64 revert the application of the projection of the rotational and
c        translational degrees of freedom. 
c     Options for ESP
c     modesp=0 the old uniform distances in dimers and trimers 
c     modesp=1 the n-body consistent distances except for connected n-mers 
c     modesp=2 the n-body consistent distances throughout 
c     modesp=4 additive dimer(trimer) ESP approximation (unconnected only)
c     modesp=8 use the modesp=4 approximation for connected cases.
c     modesp=16 FDIFF for ESPs
c     modesp=32 add exchange terms to ESP
c     modesp=64 screening for 1e ESP
c     modesp=128 Use FDIFF-like accelerator to SCC 
      modesp=-1
c
c     The length of the ESP multipolar expansion.
      ivmul=1
c
c     only two dimension parameters are user accessible. 
c
      call vclr(rflmo,1,4)
c
c     at most 10 sets of LMOs
      MAXKND=10
c     at most 5 MOs in LMOs. 5 is for C (1s + four sp3). 
      maxcao=5
c     maxbas=1
      maxbnd=-1
c     19 is for C of 6-311G*
c     maxbbd is maximum basis functions for the atoms in the model system
c     whose AO expansion is saved.
c     maxcbs is the upper bound of this number (as maxcbs is set for all atoms
c     not just those near BDAs which contribute in model systems.
      maxbbd=maxcbs
c
      nopden=0
c
      icurfg=0 
      jcurfg=0 
      kcurfg=0 
      ncursh=0
c     internal flags to alter computations of integrals 
c
c
c     initialize fmo-tddft options
      call viclr(iexcit,1,6)
      if(TDDFTYP.ne.rnone) iexcit(2)=1
      iexcit(4)=-1
      exfid=0
c     NOPFRG:
c     print-out level separately for each fragment (bit additive). 
c     0 no additional print-out 
c     1 set the value of NPRINT=7 (print orbitals).
c     2 set MVOQ to +6 to obtain better virtual orbitals.
c     4 generate cube file, the grid is chosen automatically. 
c     8 remove one electron, donor fragment (meaningful only for nbody=1).
c     16 add one electron, acceptor fragment (meaningful only for nbody=1).
c     32 use 10 times tighter icut for DIRSCF jobs involving this fragment. 
c     64 use fixed ATCHRG for ESP from this fragment.
c        ATCHRG may not have enough storage if > 1 bond per fragment is cut...
c        NFG words are allocated for redundant BDA atoms.
c     128 apply options 1,8 only for the last property iteration
c         (meaningful for correlation or GRADIENT only). 
c     modcha is the order of the many-body expansion of atomic chargers
c     (some charges are expensive and time can be saved by reducing the order).
      modcha=-1
      modlmo=0
      iatcha=0
      call vclr(ascreen,1,4)
c     EFMO
      IEFMORUN=0
c if this <= 0, it's not an analytic gradient
c if it's > 0, it's an analytic gradient
      iefmo_agrad=0
      call viclr(modefmo,1,5)
      IMODEFE=0
c     IMODEFE:
c     change the behavior of electrostatics from EFPs
c     in the EFMO method.
c     0'no screening of electrostatics
c     1'exponential screening of electrostatics by fixed value
c       to fit the classical potential to the QM-potential,
c       set screen(1)=-1 (experimental)
c     2'Add octupole energy to the electrostatic energy
c     4'use Hui Li's density based multipole expansion
c     8'ignore torque contributions to the gradient
c    16'generate electrostatics on bond midpoints too.
c
      IMODEFP=0
c     IMODEFP:
c     change the behavior of polarization from EFPs
c     in the EFMO method
c     0'tang-toenis type screening
c     1'do not include any polarization. at. all.
c     2
c     4'add percentage discrimination based on distance to atoms
c     8'ignore torque contributions to the gradient
c    16'use full polarization tensors
c    32'move polarizability tensors to nearest atom before induction
c    64'do not evaluate electrostatic field, induced dipoles or
c       gradient contributions from neighbouring fragments. NB,
c       this assumes fragments are made in a sequential fashion.
c   128'use ruednberg localization for localization of orbitals
c
      IMODEFD=0
c     IMODEFD
c     change the behavior of dispersion from EFPs
c     in the EFMO method
c     0 disable this feature
c     1 enable this feature
c     2 forcibly turn on dispersion with AFO (otherwise crashes
c       for non-diagonal Fock matrix)
c     4 (as for polarization) unannounced keyword, on by default
c       as the only implementation: add percentage discrimination
c       based on distance to atoms into gradient
c     8 (as for polarization) unannounced keyword, on by default
c       as the only implementation: ignore torque contributions
c       to gradient
c     16 (placeholder) enable dispersion with AFO by re-diagonalizing
c        Fock matrix of the fragments - to be implemented
c
      IMODEFCT=0
c     IMODEFCT
c     change the behavior of charge transfer from EFPs
c     in the EFMO method
c     0'disable this feature
c     1 enable this feature
c
      IMODEFER=0
c     IMODEFER
c     change the behavior of exchange repulsion from EFPs
c     in the EFMO method
c     0 disable this feature
c     1 enable this feature
c
c     1 frozen domain (FD): AF, AB, bB terms
c     2 in frozen domain with dimers (FDD): ignore bB
c     4 in frozen domain in buffer (FDB): ignore AF
c     (B=A+b)
c     Note that FD, FDD, FDB are chosen with 1, 3 and 7.
c      
      modfd=0
c
      modfmm=0
c     bit-additive
c     1 compute ES dimers (individual contributions)
c     2 compute ES dimers (sum only)
c     8 compute one-electron ESP gradients (RESPPC<=0) 
c
c     ncentm=1
c  
      ndualb=0
      maxrij2(1)=0
      maxrij2(2)=(nfg*nfg-nfg)/2
      ngab=0
c
      if(itask.ne.1.and.itask.ne.2) then
        write(6,*) 'Internal Task error',itask
        call abrt
      endif
      KQNAM(25)=MaxLay*10+1
      if(itask.gt.1) then
        KQNAM(3)=nfg*10+1
        KQNAM(4)=nfg*10+1
        KQNAM(5)=(natfmo+nfg*2+1+nextraind)*10+1
        KQNAM(6)=nfg*10+5
        KQNAM(7)=nfg*10+1
        KQNAM(15)=nlayer*10+1
        KQNAM(17)=nlayer*10+5
        KQNAM(18)=nfg*10+5
        KQNAM(19)=nlayer*10+5
        KQNAM(20)=nlayer*10+5
        KQNAM(24)=nfg*10+1
        KQNAM(26)=nfg*10+1
        KQNAM(28)=nlayer*10+5
        KQNAM(36)=nlayer*10+5
        KQNAM(39)=(natfmo+nfg)*10+3
        KQNAM(49)=10*nfg+1
        KQNAM(56)=nindatp*10+1
        KQNAM(57)=nmodseg*10+1
        do i=1,nfg
          ichfg(i)=0
          mulfg(i)=1
          layfrg(i)=1
          WRITE(UNIT=tstring,FMT='(A3,I5.5)') 'frg',i
          frgnam(i)=dstring 
          scffrg(i)=scftyp
          molfrg(i)=0
          nprfrg(i)=0
        enddo
        do i=1,nlayer
          mpnfmo(i)=mplevl
          fmoscf(i)=scftyp
          fmodft(i)=dftype
c         Although no check is done here, only grid-based DFT is supported. 
          fmocc(i)=cctyp
          fmoci(i)=cityp
          fmotd(i)=TDDFTYP
        enddo
        do i=1,natfmo+nfg*2+1
          indat(i)=1 
        enddo
        call vclr(atchrg,1,natfmo+nfg)
        call viclr(iactfg,1,nfg)
        call viclr(indatp,1,nindatp)
        call viclr(modseg,1,nmodseg)
      else
        nfg=1
        nlayer=1
        KQNAM(3)=10*1+9
        KQNAM(4)=10*1+9
        KQNAM(5)=10*1+9
        KQNAM(6)=10*5+9
        KQNAM(7)=10*1+9
        KQNAM(15)=10*1+9
        KQNAM(17)=10*5+9
        KQNAM(18)=10*5+9
        KQNAM(19)=10*5+9
        KQNAM(20)=10*5+9
        KQNAM(24)=10*1+9
        KQNAM(26)=10*1+9
        KQNAM(28)=10*5+9
        KQNAM(36)=10*5+9
        KQNAM(39)=10*3+9
        KQNAM(49)=10*1+9
        KQNAM(56)=10*1+9
        KQNAM(57)=10*1+9
      endif
      indatp(1)=0
      CALL NAMEIO(IR,JRET,FMOGRP,NNAM,QNAM,KQNAM,
     *            nfg,nlayer,ichfg,mulfg,indat,frgnam,layfrg,rESPAP,
     *            RESPPC,RESDIM,ritrim,MAXKND,maxcao,nbsse,mpnfmo,IXESP,
     *            fmodft,scffrg,fmoci,fmocc,NACUT,nfmopal,rcorsd,molfrg,
     *            lbody,nprfrg,modmol,fmoscf,orshft,modgrd,respct,
     *            maxbnd,modesp,ivmul,iexcit,fmotd,exfid,modcha,atchrg,
     *            rflmo, modlmo,iatcha,maxbbd,nopden,orshft2,ascreen,
     *            IEFMORUN,modfd,iactfg,modefmo,modfmm,ndualb,maxrij2,
     *            iefmo_agrad,ngab,indatp,modseg,
     *            0,0,0,0,0,0,0, 0,0,0,0,0,0)
c     if(maswrk) write(6,7777) fmoscf(1),(scffrg(i),i=1,4)
c7777 format(1x,A8,' wwwscf ',9A8)
c
      if(itask.eq.2) call setrmdind(nfg,natfmo,indat,ichfg)
      call dcopy(4,ritrim,1,restri,1)
      rijskip=ritrim(5)
      maxrij=maxrij2(1)
      nscfd=maxrij2(2)
c
      IF(JRET.EQ.2) THEN
         IF (MASWRK) WRITE(IW,*) 'Error reading $FMO'
         CALL ABRT
      END IF
      if(lbody(1).eq.-2) then
        lbody(1)=abs(lbody(1))
        fmoq=.true. 
      else
        fmoq=.false. 
      endif
      IF (MEXSTATE.EQ.2.and.itask.eq.2) THEN
c       Set up the multiplicity of the MEX fragment;
c       SCFTYP for it and all layers. 
        if(maswrk) write(iw,9100) mult2,scf2
        mulfg(IMEXFG)=mult2
        scffrg(IMEXFG)=scf2
        do i=1,nlayer
          fmoscf(i)=scf2
        enddo
      end if
c
      if(iand(modmol,8).ne.0.and.itask.eq.2) then
        do i=1,nfg
          if(molfrg(i).eq.0) then
            write(6,*) 'molfrg(',i,') should not be 0'
            call abrt
          endif
        enddo
      endif
C
      if(itask.eq.2) then 
       if(SCFTYP.eq.UHF.OR.SCFTYP.eq.ROHF) then
         call icopy(nfg,mulfg,1,mulfg(nfg+1),1)
         do ifg=1,nfg
           mulfg(ifg)=abs(mulfg(ifg))
         end do
       end if
      end if
c
C
      if(lbody(2).lt.0) call icopy(nlayer-1,lbody(1),0,lbody(2),1)
      if(lbody(1).lt.0) call icopy(nlayer,nbody,0,lbody,1)
      nbody=lbody(inamax(MaxLay,lbody,1))
      if(nbody.lt.0) nbody=0 
c     bizarre artefact: should use max value, not max abs. value in inamax. 
      if(modcha.eq.-1) modcha=nbody
c
      if(nbody.gt.3) then
        if(maswrk) write(iw,9050) nbody 
        call abrt
      endif
      if(nbody.gt.nfg) then
        if(maswrk) write(iw,9070)
        call abrt
      endif
      modgrd0=0 
      if(rflmo(1).ne.0) modgrd0=16
      if(modesp.eq.-1.and.dftbfl) modesp=0
      if(nbody.le.2) then
        if(RESPPC(1).eq.defval) RESPPC(1)=2
        if(modesp.eq.-1) modesp=0 
        if(modgrd.eq.-1) modgrd=10+modgrd0 
        call vclr(restri,1,4)
        if(nbody.le.1) then
          resdim=0
          rcorsd=0
        endif
        if(resdim.eq.defval) resdim=2
        if(rcorsd.eq.defval) rcorsd=2
      else
        if(RESPPC(1).eq.defval) RESPPC(1)=2.5D+00
        if(modesp.eq.-1) modesp=1
        if(modgrd.eq.-1) modgrd=10+modgrd0
        if(restri(1).eq.defval) restri(1)=1.25D+00
        if(restri(2).eq.defval) restri(2)=-1
        if(restri(3).eq.defval) restri(3)=2
        if(nder.eq.0) then
          if(restri(4).eq.defval) restri(4)=2
          if(resdim.eq.defval) resdim=restri(1)+restri(3)
          if(rcorsd.eq.defval) rcorsd=restri(1)+restri(4)
        else
          if(restri(4).eq.defval) restri(4)=restri(3)
c         restri(4)=restri(3) for correlated FMO3 gradient.
          if(resdim.eq.defval) resdim=0
          if(rcorsd.eq.defval) rcorsd=resdim
        endif
c       restri(2) is usually not used.
c       resdim and rcorsd are set so that FMO3 has all needed dimers for
c       better accuracy.
      endif
      if(iexcit(2).gt.nbody) then
        if(maswrk) write(iw,9060)
        call abrt
      endif
      if(TDDFTYP.ne.rnone.and.iexcit(1).eq.0) then 
        if(maswrk) write(iw,*) 'TDDFT fragment undefined in iexcit.'
        call abrt
      endif 
      if(iexcit(4).lt.0) then
        iexcit(4)=0
        if(iexcit(2).ge.2) iexcit(4)=2
      endif
c     FMO-TDDFT gradient requires iexcit(3)=1. 
      if(TDDFTYP.ne.rnone.and.nder.gt.0) iexcit(3)=1 
c
      if(maxbnd.le.0) maxbnd=nfg*2+1
c
      if(nbody.le.2) call vclr(restri,1,4)
c
c     note that the 2nd elements are not used in Gaussian(?)
      if(iand(nfmopal,1).ne.0) then
        itmp=int(respap(1)*ten) 
        respap(1)=itmp/ten
        itmp=int(resppc(1)*ten) 
        resppc(1)=itmp/ten
        itmp=int(resdim*ten) 
        resdim=itmp/ten
      endif
c     do i=1,3
c       if(rflmo(i).ne.0) rflmo(i)=rflmo(i)/units
c     enddo
      gcorrel=.false.  
      if(itask.gt.1) then
        do ilay=1,nlayer
c         gcorrel=gcorrel.or.mpnfmo(ilay).ne.0.or.fmoci(ilay).ne.rnone
c    *                   .or.fmocc(ilay).ne.rnone
c         Only one of these may be set. Check now if that is so.
          ncorme=0
          if(mpnfmo(ilay).ne.0) ncorme=ncorme+1
          if(fmoci(ilay).ne.rnone) ncorme=ncorme+1
          if(fmocc(ilay).ne.rnone) ncorme=ncorme+1
          if(fmotd(ilay).ne.rnone) ncorme=ncorme+1
          if(ncorme.gt.1) then
            write(6,*) 'Check MP,CI,CC and TD in layer ',ilay
            call abrt
          endif
          gcorrel=gcorrel.or.ncorme.ne.0
        enddo
        if(modfd.ne.0) then
          call icopy(nfg,iactfg,1,iactfg(nfg+1),1)
          call explist(nfg,iactfg(nfg+1),iactfg,nactfg)
c         write(6,*) 'list of act frg',(iactfg(i),i=1,nfg)
        endif
      endif 
      if(.not.gcorrel) then
        rcorsd=0.0D+00
        restri(4)=0.0D+00
      endif
c     if(gcorrel.and.rcorsd.eq.0.0D+00) then
c       write(iw,*) 'Please set rcorsd to a nonzero value.'
c       call abrt
c     endif
      if(restri(4).gt.restri(3).and.gcorrel) then
        write(iw,9040) restri(3),restri(4)
        call abrt
      endif
      if(iand(modesp,32).ne.0.and.nder.gt.0) then
        write(iw,*) 'No gradient for ESP with exchange.' 
        call abrt
      endif
      if(iand(modgrd,2).ne.0.and.nder.gt.0) then
        if(respap(1).ne.0) then
          write(iw,*) 'ESP derivative (modgrd=2) requires respap=0.'
          call abrt
        endif
        if(nbody.lt.2) then
          write(iw,*) 'ESP derivative (modgrd=2) requires nbody>1.'
          call abrt
        endif
        if(iand(modesp,3).ne.0) then
          write(iw,*) 'ESP derivative (modgrd=2) requires modesp=0.'
          call abrt
        endif
      endif
      IF( IEFMORUN.EQ.0 .AND.
     *     ((IMODEFE + IMODEFP + IMODEFD + IMODEFCT + IMODEFER)
     *     .GT.0)) THEN
         write(iw,*) 'MODEFE, MODEFP, MODEFD, MODEFCT, and MODEFER'
         write(iw,*) ' need IEFMO=1)'
        call abrt
      ENDIF
c     EFMO requires special settings
      if( IEFMORUN.GT.0 ) then
        if(nbody.gt.2) then
          write(iw,*) 'EFMO does not support nbody>2'
          call abrt
        endif
c       assign efmo options
        imodefe  = modefmo(1)
        imodefp  = modefmo(2)
        imodefd  = modefmo(3)
        imodefct = modefmo(4)
        imodefer = modefmo(5)
        if(imodefe.eq.0) imodefe=8
        if(imodefp.eq.0) imodefp=12
c       default to no torque and percentage discrimination to atoms
!        if(imodefd.eq.0 .or. imodefd.eq.1) imodefd=imodefd+12
!----PX fine control what terms on included for dispersion
        IF(IAND(IMODEFD,1).ne.0) THEN ! (4/3)*disp6
          write(6,*) 'IMODEFD fmoio=',IMODEFD
          DISP7=.false.; DISP8=.false.
          IDDDYN=.true.
          IDQDYN=.false.; IQQDYN=.false.; IDODYN=.false.
          E7DISP=.false.; E8DISP=.false.; do_aniso_disp=.false.
        ELSEIF(IAND(IMODEFD,16).ne.0) THEN ! (4/3)*disp6 + disp7
          write(6,*) 'IMODEFD in fmoio=',IMODEFD
          DISP7=.true.; DISP8=.false.
          IDDDYN=.true.
          IDQDYN=.true.; IQQDYN=.false.; IDODYN=.false.
          E7DISP=.true.; E8DISP=.false.; do_aniso_disp=.false.
        write(6,*) 'IDQDYN for efmo',IDQDYN
        ELSEIF(IAND(IMODEFD,32).ne.0) THEN ! disp6 + disp7 + iso disp8
          write(6,*) 'IMODEFD in fmoio=',IMODEFD
          DISP7=.true.; DISP8=.true.
          IDDDYN=.true.
          IDQDYN=.true.; IQQDYN=.true.; IDODYN=.false.
          E7DISP=.true.; E8DISP=.true.; do_aniso_disp=.false.
        ENDIF
!----END PX
c     change so that EFMO will run the Z-vector code.
c if iefmo_agrad is > 0, then modgrd should be set to 33.
c otherwise, do nothing to it.
        if( iefmo_agrad .gt. 0 .and. nder .gt. 0 )then
           modgrd = 33
        else
c
c     add modgrd
           modgrd = 1
           if(rflmo(1).ne.0) modgrd = modgrd + 16
        endif
        ixesp = 16416
      endif
      if(rflmo(1).ne.0.and.IAND(MODGRD,16).eq.0) then
        if(maswrk) write(iw,*) 'modgrd=16 should be set for AFO'
        call abrt 
      endif
c     if(IAND(MODGRD,32).ne.0.and.nder.gt.0.and.nbody.ne.2) then
c       if(maswrk) write(iw,*) 'modgrd=32 only works for FMO2' 
c       It seems that Takeshi did not finish this code, 
c       although he did start it.
c       call abrt 
c     endif
      if(ascreen(1).ne.0.and.resppc(1).ge.0.and.iefmorun.eq.0) then
        if(maswrk) write(iw,*) 'Set resppc=-1'
        call abrt
      endif
c     if(ascreen(1).ne.0.and.iand(modesp,64).eq.0) then
c       if(maswrk) write(iw,*) 'Add 64 to modesp (e.g., use modesp=64)'
c       call abrt
c     endif
c
c     FMO/BD always needs full allocation, since LMOs are regenerated.
c
      if(itask.gt.1.and.iand(ixesp,128).ne.0.and.maxcao.ne.maxcbs) then
        write(iw,*) 'Please set maxcao in $FMO to',maxcbs
        call abrt
      endif
c
c     if(resdim.ne.0.0D+00.and.nbody.gt.2.and.nder.ne.0) then
c       write(iw,*) '3-body gradient is not available with RESDIM.ne.0.'
c       call abrt
c       This is now implemented via the terms in TVDER, FMOESP
c       and ESP Lagrangians used in ESDGRD. Note that only the
c       FMO3 gradient is computed (i.e., the 2-body gradient has no ES
c       dimer contributions, because they are added directly to the 
c       3-body gradient).
c     endif
c     if(gcorrel.and.rcorsd.ne.0.0D+00.and.nbody.gt.2.and.nder.ne.0.and.
      if(gcorrel.and.rcorsd.ne.resdim.and.nbody.gt.2.and.nder.ne.0.and.
     *   iand(ixesp,64).eq.0) then
        write(iw,*) 'FMO3 gradient requires RCORSD=RESDIM.'
c       Implementing RCORSD!=RESDIM requires adding a loop 1,2 (=RHF,corr),
c       similar to monoscf.
        call abrt
      endif
c     if(ixesp.eq.-1) then
c       ixesp=0
c       if(nder.gt.0.and.resdim.gt.0) ixesp=32
c     endif
      if(iatcha.ne.0.and.itask.gt.1) then
        nblock=int(iatcha/(natfmo/nfg))
        if(nblock*natfmo.ne.iatcha*nfg) then
c       This is a condition to ensure that iatcha is a multiple of the first 
c       fragment in terms of the number of atoms. 
          if(maswrk) write(iw,*) 'Running uneven copy of charges'
          nblock=int(natfmo/iatcha)
          if(nblock*iatcha.ne.natfmo) call abrtx("Confusion in fmominp")
          do i=2,nblock
            call dcopy(iatcha,atchrg(1),1,atchrg((i-1)*iatcha+1),1)
          enddo
        else
          do i=2,nfg/nblock
            call dcopy(iatcha,atchrg(1),1,atchrg((i-1)*iatcha+1),1)
          enddo
        endif
        if(maswrk) write(iw,*) 'Copied charges of block',nblock,
     *                         ' as requested.'
c       write(6,*) 'aaa',(atchrg(i),i=1,natfmo)
      endif
c
c     nfg is some estimate of the cross-cuts (extra bond cuts due to X-shaped
c     crossovers in proteins; S-S links, mostly).
c     mxatmf=mxatm*nfg
      maxvec=100
      maxabd=4
c     if(rflmo(1).ne.0) maxabd=8
c     maxcbs will be set elsewhere
      return
 9040 format(1x,'Invalid trimer option ritrim(3)=',F10.5/
     *       1x,'      is smaller than ritrim(4)=',F10.5)
 9050 format(1x,'The only meaningful values for nbody are: 0,1,2,3',
     *          'and you have nbody=',I3/)
 9060 format(1x,'nbody should not be larger than iexcit(2).',/)
 9070 format(1x,'nbody should not be larger than nfg.',/)
 9100 format(/1x,'MEX sets multiplicity to',I2,' and SCFTYP to ',A8,/)
      end
c
C*MODULE fmoio   *DECK fmopinp
C>
C>     @brief $FMOPRP
C>
C>     @details Read $FMOPRP.
C>
C>     @author Dmitri Fedorov
C>
      SUBROUTINE fmopinp(itask,mannod,mcnvfmo,mconvex,mconfg,ijvec,
     *                   ngrfmo,prtdst,ncvscf,natdata,numdata,ibfconv,
     *                   mapconv,ngrmax,savgrd,loadbf,loadgr,l1dir,ngm,
     *                   vdwrad,grdpad,mpcmit,convpcm,pcmoff,imect,n0bda
     *                  ,ne0bda,r0bda,e0bda,efmo0,nefmo0,nepl0,epl0ds,
     *                   eint0,rappri,irestl,naoafo,modcnv,nlcmo,offnum,
     *                   nfmobuf,needorb,reapc,maxmand,mandist)
      use mx_limits, only: mxgrid
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      logical GOPARR,DSKWRK,MASWRK,ISGDDI,PAROUT,INITGDDI,savgrd,SG1,
     *        DIRSCF,FDIFF,wasgddi,MLGDDI,DFTBFL,SCC,SRSCC,DFTB3,DAMPXH,
     *        LCDFTB,needorb
      PARAMETER (nnam=49,MAXNZ=137,nstfmo=10,maxpst=10,ONE=1.0D+00)
      DIMENSION QNAM(NNAM),KQNAM(NNAM),mannod(*),mcnvfmo(*),mconfg(*),
     *          ijvec(5,*),ngrfmo(maxpst,*),prtdst(4),ncvscf(2),
     *          ibfconv(*),mapconv(*),loadbf(maxpst,*),loadgr(maxpst,*),
     *          ngm(3),vdwrad(MAXNZ),grdpad(3),vdwr0(MAXNZ),r0bda(*),
     *          e0bda(5,ne0bda,*),efmo0(*),epl0ds(*),eint0(4),rappri(3),
     *          nlcmo(3),nfmobuf(3),reapc(3,2,MAXNZ),rh(MAXNZ),
     *          mandist(3,maxmand)
      PARAMETER (MAXEXP=6,MXSPE=10)
      COMMON /DFTB  / DFTBFL,SCC,SRSCC,DFTB3,DAMPXH,LCDFTB
      COMMON /DFTBAO/ EXPON(MAXEXP,3*MXSPE),COEFF(3,MAXEXP,3*MXSPE),
     *                NEXP(3*MXSPE),INDSH(MXSPE+1)
      COMMON /DFGRID/ DFTTHR,DFTGTHR,SWOFF,SW0,BSLRD(137),NDFTFG,
     *                NRAD,NTHE,NPHI,NRAD0,NTHE0,NPHI0,
     *                NANGPT(MXGRID),NANGPT0(MXGRID),SG1,JANS
      COMMON /GDDI/   ISCOPE,NGROUPS,MYGROUP,MEGLOB,npglob,nnglob,JBTYP,
     *                ISGDDI,PAROUT,INITGDDI,wasgddi,MLGDDI,NSUBGR,
     *                MeUniv,NPUniv,numdlb,myworld,nworlds,mogddi
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON /MACHSW/ KDIAG,ICORFL,IXDR,modio,mem10,lpnt10,mem10m
      COMMON /OPTSCF/ DIRSCF,FDIFF
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
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
      common /fmodim/ maxbnd,maxknd,maxcbs,maxcao,maxvec,maxl1,maxnat,
     *                maxabd,maxbas,maxbbd,maxlmo,maxslo,maxabd2,maxrij
      DATA KQNAM /1,1,3,1,1, 1,1,21,1,1, 93,1,43,1,3, 1,3,1,1,1, 
     *            1,1,1,1,31, 3,33,1,3,3, 1,1,1,3,3, 3,3,43,33,1,
     *            3,1,1,31,3, 31,1,3,1/
      DATA QNAM/8HNGUESS  ,8HMAXIT   ,8HCONV    ,8HMODORB  ,8HMODPAR  ,
     *          8HMCONV   ,8HIREST   ,8HNCVSCF  ,8HIJVEC   ,8HNPRINT  ,
     *          8HESPSCA  ,8HNGRFMO  ,8HPRTDST  ,8HMANNOD  ,8HCNVDMP  ,
     *          8HMCONFG  ,8HCOROFF  ,8HIBFCON  ,8HMAPCON  ,8HMAXAOC  ,
     *          8HMODPRP  ,8HLOADBF  ,8HLOADGR  ,8HNAODIR  ,8HNGRID   ,
     *          8HVDWRAD  ,8HGRDPAD  ,8HNPCMIT  ,8HCNVPCM  ,8HPCMOFF  ,
     *          8HIMECT   ,8HIPIEDA  ,8HN0BDA   ,8HR0BDA   ,8HE0BDA   ,
     *          8HEFMO0   ,8HEPL0DS  ,8HEINT0   ,8HRAPPRI  ,8HIRESTL  ,
     *          8HCNVAFO  ,8HNAOAFO  ,8HMOFOCK  ,8HNLCMO   ,8HOFFNUM  ,
     *          8HNBUFF   ,8HMODPAN  ,8HREAPC   ,8HMANDST  /
      data fmogrp/8HFMOPRP  /,rnone/8HNONE    /
c     note: unknown parameters are set to 2.5A.
      data vdwr0/1.20d+00,1.20d+00,1.37d+00,1.45d+00,1.45d+00,1.50d+00,
     *          1.50d+00,1.40d+00,1.35d+00,1.30d+00,1.57d+00,1.36d+00,
     *          1.24d+00,1.17d+00,1.80d+00,1.75d+00,1.70d+00,12*2.5d+00,
     *          2.5d+00,4*2.5d+00,2.3d+00,102*2.5d+00/
      data rh/0.74,0,   0,   0,   1.19,1.09,1.01,0.96,0.92,0,
     *        0,   0,   0,   1.48,1.44,1.34,1.27,0,   0,   0,
     *        0,0,0,0,0, 0,0,0,0,0,
     *        0,   1.53,1.52,1.46,1.41, 0,0,0,0,0,
     *        0,0,0,0,0, 0,0,0,0,1.70,
     *        0,   1.70,1.61,0,0, 0,0,0,0,0,
     *        77*0/
c     http://www.wiredchemist.com/chemistry/data/bond_energies_lengths.html
c
c     this subroutine reads FMO property input.
c     care should be exercised about e0bda: it is used both as
c     e0bda(5,n0bda) and e0bda(5,nbdfg,nlayer)
c
      nguess=2
c     bit additive
c        1 use free molecule guess for monomers
c        2 use monomer density for dimer guess (if 0, use Huckel guess)
c        4 insert HMO projection operator in the Huckel guess 
c        8 use SCF dimer density from previous geometry if available
c       16 do RHF for each dimer before MCSCF and match orbitals
c       32 do not construct density and orbitals for n-mers in DFTB/HOP
c       64 m-mer Aufbau: for dimers (trimers), fill in the initial guess
c          according to the Aufbau principle applied to orbital energies
c          of all (2 or 3) fragments in the n-mer as 1 set;
c          if 64 is not set, then 1-mer Aufbau is used, always use occupied 
c          monomer orbitals as occupied n-mer orbitals.
c          bit 64 can lead to populating virtual monomer orbitals and
c          enforce the charge transfer states.
c      128 do not do restart during geometry optimisations (that is, always
c          use the same (Huckel) guess for monomers).
c      256 the "desperate" option: if SCF does not converge, try the
c          alternative converger (flip between SOSCF and DIIS).
c      512 manual orbital reordering, i.e. reorder initial orbitals using
c          $GUESS options (IORDER), useful for MCSCF.
c     1024 orthogonalise dimer/trimer initial orbitals (relevant for
c          DFT, ROHF, MCSCF, PIEDA(required))
c     2048 a modifier of 256: when going "desperate" do not read orbitals
c          of the previous (diverged) SCF.
c     262144 apply dimer cut orbital projection to the dimer initial guess
c     Options 4,16,32,64,262144 are effective only if some bonds are cut.
c     Options 64,128,512,1024 are obscure and almost never used.
c     Options 256 and 2048 may be useful for poor convergence.
      mxitfg=-1
c     monomer SCF convergence (energy)
      convfg=1.0D-07
c     AFO convergence
      cnvafo=1.0D-05
c     a threshold to force no static correlation (MCSCF/DFT) until energy 
c     converges up to coroff. If 0, this option is ignored.
c     This option is supposed to be useful mainly for DFT. Sensible values
c     are about 1e-2 ... 1e-4. Beside helping convergence, this option can 
c     speed the run up. 
c     coroff is not set if DFTTYP is not defined in $CONTRL.
c     This can be very confusing.  
      coroff=0.0d+00
      if(ndftfg.ne.0.and.iand(ixesp,16384).eq.0.and.iand(modesp,8).eq.0) 
     *  coroff=1.0d-03
c
c     Dual basis runs use both COROFF and SWOFF in $DFT.
c
      prtdst(1)=0.0D+00
c     all interfragment distances less then prtdst(1) will be printed.
      prtdst(2)=0.5D+00
c     nonzero interfragment distances less then prtdst(2) will be marked.
c     fragments closer than prtdst(3) will be obliged to have a fractioned
c     bond defined in $FMOBND. -1.0D+00 effectively turns off such check.
c     0.6 is a good value to check C-C bonds.
c     The cost of checking will grow rapidly with prtdst(3) (not parallelised).
      prtdst(3)=0.6D+00
      prtdst(4)=0.0D+00
c 
c     bit additive
c     MODORB=0 exchange density between fragments
c     MODORB=1 exchange MOs between fragments
c     MODORB=2 orbital energies between fragments
c     For DFT MODORB=1 must be set (i.e., MODORB=1 or MODORB=3).
c     DFT will run with either of the two, but MODORB=1 does not properly
c     reorder the initial dimer orbitals, so MODORB=3 is the only proper choice.
c     MCSCF requires both orbitals and energies (MODORB=3 only).
      MODORB=0
      if(NDFTFG.ne.0.or.TDDFTYP.ne.rnone.or.needorb) MODORB=3
c     if(NDFTFG.ne.0.or.SOSCF) MODORB=modorb+1
c     if(SOSCF)                modorb=modorb+2
c
c     MODPAR controls parallel options (bit-additive)
c     MODPAR=  1 turns on/off heavy job first strategy (aimed to reduce
c                waiting on remaining jobs at barrier points)
c     MODPAR=  2 changes ESP parallisation strategy:
c                if not set, parallelise loops over shells in each fragment,
c                otherwise parallelise loop over fragments.
c                Note that the latter way is seldom useful (i.e., many short-
c                contracted shells on slow network and few nodes). The other 
c                option is preferred, although it generates more parallel 
c                traffic to perform dynamic load balancing.
c     modpar=  4 broadcast all fragments done by a group at once rather
c                than fragment by fragment. 
c     modpar=  8 alters the behaviour of the fragment
c                initialisation: if set, fragments are always done in the
c                reverse order (nfg, nfg-1, ...1) because distance calculation 
c                costs decrease in the same order and they usually prevail over
c                making Huckel orbitals or running free monomer SCF. Note that
c                during SCC (monomer SCF) iterations the order in which monomers
c                are done is determined by MODPAR=1. 
c     modpar= 16 if set, LMO projectors (FMOHOP) will not be parallelised 
c     modpar= 32 save memory by not storing some n-mer data (n>1).
c                Note that this option requires certain other options to work.
c     modpar= 64 Distribute file F40 for restarts irest>1. The file should be
c                put only on the grand master node and it will broadcast it to
c                all other nodes. Moreover, if the file F40 is copied to slaves
c                and modpar=64, then file errors will result due to 
c                improper indexing of the existing file.
c     modpar=128 flip fragment indices in separated dimers. It improves
c                parallelisation for massively parallel runs.
c     modpar=256
c     modpar=512 Use DDI memory to store fragment densities during the monomer
c                step, using supervector (smallest memory).
c     modpar=1024 Use DDI memory to store fragment densities during the monomer
c                step, using matrix (smallest communications).
c     modpar=2048 Slaves and master do intergroup I/O for modpar=512 or 1024  
c                 (0 here means only master does it). 
c                
      MODPAR=-13
c
      call viclr(ngrfmo,1,maxpst*nlayer)
c     Set the number of GDDI groups:
c     ngrfmo(1) SCF monomers
c     ngrfmo(2) SCF dimers
c     ngrfmo(3) all trimers 
c     ngrfmo(4) correlated monomers 
c     ngrfmo(5) separated dimers
c     ngrfmo(6) SCF monomers after MCSCF monomers
c     ngrfmo(7) SCF dimers after MCSCF dimers 
c            8,9,10 reserved
c     These are defined inidivually for each layer.
c     The default is all 0s, using whatever is in $GDDI.
c  
c     mannod in $FMO similarly to $GDDI defines manual node division into 
c     groups separately for each ngrfmo(i) above.  
c     -1 defines automatic division. Note that to use manual division NGRFMO
c     must also be set.
      nprfmo=1
c     bits 1-2
c      0 normal output
c      1 reduced output
c      2 minimum output (not implemented yet)
c      4 print interfragment distances
c      8 print Mulliken charges
c     16 special test run to check for missing bonds in FMOBND. 
c        Use only with nbody=0!
c     256 dump efp information during during generation (EFMO)
c
c     the ESP from the Fock matrix and rediagonalise.
      call dacopy(9,one,espsca,1)
      espscf=one
c     scale ESP during the first three iterations by multiplying by ESPSCA.
      cnvdmp=0
c
c     Set a converger for each step in FMO
      do mfmostp=1,nstfmo
        mcnvfmo(mfmostp)=-1
      enddo
      call viclr(loadbf,1,maxpst*nlayer)
      call viclr(loadgr,1,maxpst*nlayer)
c     loadbf and loadgr define semidynamic load balancing.
c     n-mers with more AOs than loadbf are done statically on the first 
c     loadgr groups, e.g. loadbf(1) and loadgr(1) are used for monomer SCF
c     (see ngrfmo(1:maxpst) description for explanation of steps). 
c     Other groups will use dynamic load balancing. When all static tasks
c     are done, dynamic load balancing will be used by the selected loadgr
c     groups as well. The following conditions must be satisfied in order to use
c     semidynamic load balancing: 
c     a) GDDI in use 
c     b) nonzero loadbf(ifmostp)+loadgr(ifmostp), 
c     c) $GDDI BALTYP=NXTVAL 
c     d) iand(modpar,1).ne.0 . 
c     It appears that ngrfmo(6) should be 0 as well (no MCSCF selected monomer
c     group division). 
c     reapc(i,j,k):
c     i=1 k-H distance (A)
c     i=2 BDA energy (hartree)
c     i=3 BAA energy (hartree)
c     j=1,2 basis set
c     k is atomic number
c
c     rappri is in the order: rmin,rmax,rstep
      call vclr(rappri,1,3)
      if(itask.ne.1.and.itask.ne.2) then
        write(6,*) 'Internal Task error',itask
        call abrt
      endif
      KQNAM(6)=nstfmo*10+1
      KQNAM(12)=maxpst*nlayer*10+1
      KQNAM(22)=maxpst*nlayer*10+1
      KQNAM(23)=maxpst*nlayer*10+1
      KQNAM(26)=MAXNZ*10+3
      KQNAM(49)=3*maxmand*10+1
      call viclr(mandist,1,3*maxmand)
      if(itask.gt.1) then
        KQNAM(9)=5*maxvec*10+1
        KQNAM(14)=ngrmax*maxpst*10+1
        KQNAM(16)=(nfg+mconvex*4)*10+1
        KQNAM(18)=2*natdata*10+1
        KQNAM(19)=natdata*maxl1c*10+1
        KQNAM(34)=10*ne0bda+3
        KQNAM(35)=10*5*ne0bda*nlayer+3
        KQNAM(36)=10*nefmo0*nlayer+3
        KQNAM(37)=10*nepl0+3
        KQNAM(48)=3*2*MAXNZ*10+3
        do i=1,nfg+mconvex*4
          mconfg(i)=-1
        enddo
        do i=1,ngrmax*maxpst*nlayer
          mannod(i)=-1
        enddo
        call viclr(ijvec,1,5*maxvec)
        do i=1,2*natdata
          ibfconv(i)=0
        enddo
        do i=1,natdata*numdata
          mapconv(i)=-1
        enddo
c       write(6,*) 'wwwe0z?',nbdfg,ipieda
c       if(nbdfg.ne.0.and.ipieda.ne.0) then
        if(ne0bda.ne.0) then
          call vclr(r0bda,1,ne0bda)
          call vclr(e0bda,1,5*ne0bda*nlayer)
c         write(6,*) 'wwwe0z',5*ne0bda*nlayer
        endif
        if(nefmo0.ne.0) call vclr(efmo0,1,nefmo0*nlayer)
        if(nepl0.ne.0) call vclr(epl0ds,1,nepl0)
        call vclr(reapc,1,3*2*MAXNZ)
c       Set the default values of atom-H distances.
c       APC will will bomb out if a zero value is attempted to be used.
        do i=1,MAXNZ
          reapc(1,1,i)=rh(i)
          reapc(1,2,i)=rh(i)
        enddo
      else
        KQNAM(9)=10*1+9
        KQNAM(14)=10*1+9
        KQNAM(16)=10*1+9
        KQNAM(18)=10*1+9
        KQNAM(19)=10*1+9
        KQNAM(34)=10*3+9
        KQNAM(35)=10*3+9
        KQNAM(36)=10*3+9
        KQNAM(37)=10*3+9
        KQNAM(48)=10*3+9
      endif
c
c     restart options: step and where is step
c     at present only layer 1 steps can be restarted.
c     Supported options: 
c       irststp=2 resume monomer SCF 
c       irststp=4 resume dimer SCF 
c       irststp=5 add BSSE corrections
c     n-multiple of 100 can be added to skip n layers.
c     e.g. 702 means skip 7 layers (i.e. jump to layer 8), step 2.
c     all restart jobs require F40 with monomer densities.
c
      irststp=0
      irestl=0
      if(modfd.ne.0) irestl=1
c     Force restarts from layer 1. Actually, restarts should ALWAYS do that? 
c
c     switch the SCF converger after ncvscf(1) monomer SCF iterations (that is,
c     switch between SOSCF <-> FULLNR). Enforce the converger in MCONV(2) after
c     ncvscf(2), i.e. overwrite convertors in MCONFG. 
      ncvscf(1)=9999
      ncvscf(2)=9999
c  
      modprp=0
c     FMO properties:
c      1 electron density (AO-based matrix, in-core).
c      2 electron density (AO-based matrix, disk-based) - reserved.
c      4 electron density (on a grid, in-core); produces a Gaussian cube file. 
c      8 electron density (on a grid, in-core); produces a sparse cube file. 
c     16 automatically generate grid for modprp = 4 or 8.
c     64 fast PCM without separation of pairwise screening from monomer Ees.
c
      l1dir=0
c     parameter enforcing direct SCF if the current basis size is small.
c     l1dir is the maximum basis set size per 100 mln integrals.
c     The actual condition to enforce DIRSCF is:
c     (L1/l1dir)**4 * 10d+8 * k .gt. nintic*nproc
c     where L1 is the basis set size, k is 1.5 for LABSIZ=1,NWDVAR=2 or else 2.
c     nintic is the memory for the incore integrals buffer. 
c     The recommended value is 180 (based on experience).
      call viclr(ngm,1,3)
c     NGM defines the number of grid points for monomer cube files.
c     The choice of monomers for which cube files are to be computed 
c     is done in NPRFRG.
      call dcopy(MAXNZ,vdwr0,1,vdwrad,1)
c
c     Grid padding parameter, unitless, as it is used as a factor
c     on which atomic radii are multipiled, i.e., the grid padding
c     is done by adding the extra space of grdpad*Rvdw(a) around atom a. 
c     The first value determines the box size (either total grid or monomer
c     option NOPFRG=4), the second is used to reduce the total box size into 
c     smaller windows for n-mers (total or sparse dense grid). 
c     The third value is for NGRID option to reduce the space of dimers for
c     a given monomer.
c
      grdpad(1)=2.0D+00
      grdpad(2)=-2.0D+00
      grdpad(3)=-2.0D+00
c
      mpcmit=30
      convpcm=1.0D-07
      pcmoff=0
c
c     charge transfer calculation method. 0,1,2,3,4,5 are supported.
c     spin transfer for open-shell methods is always computed with method 5.
c
      imect=4
c
      ipieda=0
c
      n0bda=0
      call vclr(eint0,1,4)
      naoafo=0
c
      mofock=0
c     Fock matrix construction
c     bit additive
c     1 construct the total Fock F and overlaps S 
c     2 add exchange to the total Fock matrix (FMO/FX)
c     4 punch T (kinetic) and H (1e) integrals 
c     8 do not use resppc(2) for the recomputed ESP
c     (prior to 5.2, do not increase RESPPC by 2 for the recomputed ESP)
c     16 compute overlaps for ES dimers
      nlcmo(1)=0
c     LCMO options
c     1: do energies
c     2: add X (LCMOX)
c     4: do MOs 
c     8: do not use resppc(2) for the recomputed ESP
c     16: print Fock
c     32: print overlap
      nlcmo(2)=0
      nlcmo(3)=0
      offnum=1.0D-03
c     offset for numerical gradient in bohr
      call viclr(nfmobuf,1,3)
c     nfmobuf(1) is used to reduce the memory needed for FMO3
      modpan=0
c     modpan=1 Do partition analysis
c     modpan=2 split "rep" and "disp" terms from nones for dimers
c              (if not set, nones=rep+diso is computed)
c     modpan=4 use PL state from FMO1
c     modpan=8 split residues into backbonds and side chains (requires PDB).
c     modpan=16 define monomer and dimer REP
c               (if not set, then define partial monomer REP)
c
      CALL NAMEIO(IR,JRET,FMOGRP,NNAM,QNAM,KQNAM,
     *            NGUESS,mxitfg,convfg,modorb,modpar,mcnvfmo,irststp,
     *            ncvscf,ijvec,nprfmo,espsca,ngrfmo,prtdst,mannod,
     *            cnvdmp,mconfg,coroff,ibfconv,mapconv,maxl1c,modprp,
     *            loadbf,loadgr,l1dir,ngm,vdwrad,grdpad,mpcmit,convpcm,
     *            pcmoff,imect,ipieda,n0bda,r0bda,e0bda,efmo0,epl0ds,
     *            eint0,rappri,irestl, cnvafo,naoafo,mofock,nlcmo,
     *            offnum,nfmobuf,modpan,reapc,mandist,0, 0,0,0,0,0,
     *            0,0,0,0,0, 0,0,0,0,0,0,0,0,0, 0)
      IF(JRET.EQ.2) THEN
         IF (MASWRK) WRITE(IW,*) 'Error reading $FMOPRP'
         CALL ABRT
      END IF
      modcnv=0
      if(convfg.lt.0) then
        modcnv=1
        convfg=-convfg
      endif
      if(mxitfg.eq.-1.and.coroff.ne.0) mxitfg=40 
      if(mxitfg.eq.-1) mxitfg=30 
      if(iand(modprp,4).ne.0.and.dftbfl.and.natfmo.gt.1) then
        if(INDSH(2).eq.0) then
          if(maswrk) write(iw,*) 'Please add $DFTBAO $end to the input.'
          call abrt
        endif
      endif
      if(iand(modprp,32).eq.0) then
         if(grdpad(2).eq.-2.0D+00) grdpad(2)=grdpad(1)
         if(grdpad(3).eq.-2.0D+00) grdpad(3)=grdpad(1)
      else
c        Use large values for MEP to turn off improper screening.
         if(grdpad(2).eq.-2.0D+00) grdpad(2)=1d+10
         if(grdpad(3).eq.-2.0D+00) grdpad(3)=1d+10
      endif
      UNITS=1.0D+00/0.52917724924D+00
c     offnum=offnum*UNITS
      if(iand(modio,1024).ne.0) then
        if(nprfmo.eq.1) nprfmo=3+128
        if(modpar.eq.-13) modpar=9+1024
c       if(nguess.eq.2) nguess=nguess+32
      else
        if(modpar.eq.-13) modpar=13
      endif
      if(iand(modio,2048).ne.0) then
        if(iand(modprp,64).eq.0) modprp=modprp+64
      endif
c
c     replace level default values by the general default 
      do mfmostp=1,nstfmo
c       -1 in fmoconv resets to defaults
       if(mcnvfmo(mfmostp).eq.-1) then
         call fmoconv(mcnvfmo(mfmostp),-1,0,dum,.true.,.true.,.false.,
     *                .false.,.false.)
       else
         if(iand(mcnvfmo(mfmostp),512).eq.0.and.dirscf) then
           if(maswrk) write(iw,9030) 
           call abrt 
         endif
       endif
      enddo
      savgrd=iand(irststp,1024).ne.0
      if(irststp.ge.1024) irststp=irststp-1024
c     An option to save gradient after each dimer calculation for the purpose
c     of the following RUNTYP=GRADIENT restarts.
      if(irststp.gt.nstfmo.or.irestl.gt.nlayer) then
        write(iw,*) 'Invalid restart options',irststp,irestl
        call abrt
      endif
      if(iand(modpar,2).ne.0.and.goparr) then
        if(maswrk) write(iw,*) 'MODPAR=2 option has an unfixed bug.'
c       The bug has something to do with GOPARR set to false in parallel
c       execution; it shows up if NPROC>NFG-1??
        call abrt
      endif 
c     if(.not.parout.and.n1.ne.0.and.(n1.ne.n2.or.n1.ne.ngroups)) then
c       if(maswrk) write(iw,9020) 
c     endif 
      if(itask.gt.1) then
        loop=0
        do k=1,nlayer
          do i=1,maxpst
            nnodi=0 
            ngrik=ngrfmo(i,k)
            if(ngrik.ne.0.and.nsubgr.lt.0) then
              if(maswrk) 
     *          write(iw,*) 'Remove NGRFMO from $FMOPRP for NSUNBGR=-1'
              call abrt
            endif
            do j=1,ngrik
              loop=loop+1
              nnodi=nnodi+mannod(loop)
            enddo
c           nnodi is negative if nodes are not divided manually
            if(nnodi.gt.0.and.nnodi.ne.nnglob.or.ngrfmo(i,k)
     *         .gt.nnglob) then
              if(maswrk) write(iw,9020) i,k,ngrfmo(i,k),nnodi,nnglob
              call abrt
            endif
            if(ngrfmo(i,k).gt.0.and..not.isgddi) then
              if(maswrk) write(iw,9025) i,k,ngrfmo(i,k) 
              call abrt
            endif
          enddo 
        enddo
        do i=1,MAXNZ
          reapc(1,1,i)=reapc(1,1,i)*UNITS
          reapc(1,2,i)=reapc(1,2,i)*UNITS
          if(ipieda.ne.0) then
c           In PIEDA/APC, cap corrections have to be split,
c           which is not supported.
c           Use no cap corrections if doing PIEDA for now.
            reapc(2,1,i)=0
            reapc(3,1,i)=0
            reapc(2,2,i)=0
            reapc(3,2,i)=0
          endif
        enddo
        if(iand(ndualb,8).ne.0.and.ipieda.ne.0.and.maswrk)
     *   write(iw,*) 'All cap corrections are cleared because of PIEDA.'
      endif
C     FOR PIEDA or CPHF CALCULATION
      if(ipieda.gt.0.and.iand(modorb,1).eq.0) then
        if(maswrk) write(iw,*) 'Setting modorb=1 flag...'
        modorb=ior(modorb,1)
      endif
c     triminid cannot handle MODORB=1 for PIEDA, set to 3.
      if((IAND(MODGRD,32).NE.0.or.ipieda.gt.0.and.nbody.eq.3.or.
     *   nlcmo(1).ne.0).and.iand(modorb,3).ne.3) then
        if(maswrk.and.iand(nprfmo,3).lt.3) 
     *    write(iw,*) 'Setting modorb=3 flag...'
        modorb=ior(modorb,3)
      endif
c     if(iand(modorb,1).ne.1.and.iand(nguess,1024).eq.0) then
      if(iand(modorb,3).ne.3.and.iand(nguess,1024).eq.0) then
        nguess=nguess+1024
c       The option to orgonalise initial guess requires MOs,
c       turn it off.
      endif
c     if(itask.gt.1) then
c       do k=1,nlayer
c         if(fmodft(k).ne.rnone) then
c           if(maswrk) write(iw,*) 'Setting modorb=3 flag...'
c           modorb=ior(modorb,3)
c           goto 100
c         endif
c       enddo
c 100   continue
c     endif
c     if(ipieda.ne.0) nguess=ior(nguess,1024)
c     Force dimer initial orbital orthonormalisation for PIEDA
C
c     if(ncvscf.eq.0) ncvscf=9999
c     if(ixesp.eq.-1) then
c       ixesp=0
c       if(nder.gt.0.and.resdim.gt.0) ixesp=32
c     endif
      return
 9020 format(/1x,'Wrong group division (ngrfmo/mannod) for step',I2,
     *           ', layer',I1,':',3I6,/)
 9025 format(/1x,'ngrfmo cannot be used without GDDI (set NGROUP in ',
     *           '$GDDI).', 3I6,/)
 9030 format(/1x,'For $scf dirscf=.t., MCONV should have 512 added to', 
     *           ' it (DIRSCF),',
     *       /1x,'or 768 (recommended, DIRSCF+FDIFF).',/)
      end
c
C*MODULE fmoio   *DECK fmocinp
      SUBROUTINE fmocinp(fmozan,fmoc,fmomas,dolat)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      logical GOPARR,DSKWRK,MASWRK,dolat
      PARAMETER (nnam=14)
      DIMENSION QNAM(NNAM),KQNAM(NNAM),u(3,3),ctr(3),
     *          fmozan(*),fmoc(3,*),fmomas(*)
      COMMON /FMCOM / X(1)
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      Common /fmoinf/ nfg,nlayer,natfmo,nbdfg,naotyp,nbody,nsegm
      common /fmolat/ untang(3),untorg(3),respbc(4),abclat(3),anglat(3),
     *                symtra(3,24),symope(3,3,24),iatorg,nsymop,maxklms,
     *                ioporg(3),iopdir(3),iopabc(3),iopang(3)
      DATA KQNAM /1,33,33,43,33, 33,1,3,3,31, 31,31,31,1/
      DATA QNAM/8HIATORG  ,8HUNTANG  ,8HUNTORG  ,8HRESPBC  ,8HABCLAT  ,
     *          8HANGLAT  ,8HNSYMOP  ,8HSYMTRA  ,8HSYMOPE  ,8HIOPORG  ,
     *          8HIOPDIR  ,8HIOPABC  ,8HIOPANG  ,8HMAXKLM  /
      data fmogrp/8HFMOLAT  /
c
c     this subroutine reads lattice FMO input.
c
      tobohr=1.0D+00/0.52917724924D+00
      torad=acos(0.0D+00)/90.0D+00
c     iatorg: where to put the origin.
c              <0 at (0,0,0) 
c              =0 center of mass (default),
c              >0 at iatorg-th atom.
      KQNAM(8)=3*48*10+3
      KQNAM(9)=3*3*48*10+3
      iatorg=0
      nsymop=0
      do i=1,3
        untang(i)=0
        untorg(i)=0
        abclat(i)=0
        anglat(i)=0
        ioporg(i)=0
        iopdir(i)=0
        iopabc(i)=0
        iopang(i)=0
      enddo
      call vclr(respbc,1,4)
      call vclr(symtra,1,3*24)
      call vclr(symope,1,3*3*24)
      maxklms=10000
c
      CALL NAMEIO(IR,JRET,FMOGRP,NNAM,QNAM,KQNAM,
     *            iatorg,untang,untorg,respbc,abclat,anglat,nsymop,
     *            symtra,symope,ioporg,iopdir,iopabc,iopang,maxklms,0,
     *            0,0,0,0,0,  0,0,0,0,0, 0,0,0,0,0,
     *            0,0,0,0,0,  0,0,0,0,0, 0,0,0,0,0,
     *           0,0,0,0,0,  0,0,0,0,0, 0,0,0,0,0, 0,0,0,0, 0,0,0,0,0,0)
      IF(JRET.EQ.2) THEN
         IF (MASWRK) WRITE(IW,*) 'Error reading $FMO'
         CALL ABRT
      END IF
      dolat=JRET.EQ.0
      if(.not.dolat) then
        maxklms=0
        return
      endif
c
      call dscal(3,torad,untang,1)
      call dscal(3,torad,anglat,1)
c     if(iuntrd.eq.1) call dscal(3,tobohr,abclat,1)
      call dscal(3,tobohr,abclat,1)
      call dscal(4,tobohr,respbc,1)
c
      vol=cellvol(abclat,anglat)
      if(maswrk) write(iw,9976) vol/tobohr**3
c     orientation of molecule
      call euleru(untang(1),untang(2),untang(3),u)
c     shift center-of-mass to origin
      iwhere=0
      if(iatorg.lt.0) then
        call vclr(ctr,1,3)
      else if(iatorg.eq.0) then
        call fndcntr(natfmo,fmozan,fmoc,fmomas,iwhere,ctr(1),ctr(2),
     *               ctr(3))
      else
        call dcopy(3,fmoc(1,iatorg),1,ctr,1)
      endif
      do iat=1,natfmo
         xx=fmoc(1,i)-ctr(1)
         yy=fmoc(2,i)-ctr(2)
         zz=fmoc(3,i)-ctr(3)
         fmoc(1,i)=xx*u(1,1)+yy*u(1,2)+zz*u(1,3)
         fmoc(2,i)=xx*u(2,1)+yy*u(2,2)+zz*u(2,3)
         fmoc(3,i)=xx*u(3,1)+yy*u(3,2)+zz*u(3,3)
      enddo
c     check symmetry operations 
      CALL VALFM(LOADFM)
      lct=LOADFM+1
      last=lct+3*natfmo*2
      NEED=LAST-LOADFM-1
      CALL GETFM(NEED)
      call chksymop(natfmo,fmoc,x(lct),x(lct+3*natfmo))
      CALL RETFM(NEED)
      return
 9976 format(' cell volume (A**3) = ',f20.4)
      end
c
C*MODULE fmoio   *DECK fmoxyz
      SUBROUTINE fmoxyz(gprnam,iscan,natfmo,fmozan,fmoc,izbas)
      use mx_limits, only: mxatm
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      character*99 STR,STRF77
      character*4 atm,c4dum
      dimension iargs(5)
      character*8 atname,ATOMNM,gprnam
      dimension fmoZAN(*),fmoC(3,*),izbas(*)
      logical GOPARR,DSKWRK,MASWRK,dolend,usexyz
      COMMON /INFO  / Cdum(MXATM,3),IAN(MXATM),NATOMS,IUNTRD,ATM(MXATM)
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /MACHSW/ KDIAG,ICORFL,IXDR,modio,mem10,lpnt10,mem10m
      common /fmodim/ maxbnd,maxknd,maxcbs,maxcao,maxvec,maxl1,maxnat,
     *                maxabd,maxbas,maxbbd,maxlmo,maxslo,maxabd2,maxrij
      common /fmoopt/ espsca(9),RESPAP(2),rESPPC(2),rESDIM,restri(4),
     *                rcorsd,respct,convfg,cnvdmp,coroff,rflmo(4),
     *                orshft,orshft2,cnvafo,ascreen(4),IXESP,mxitfg,
     *                nguess,NBSSE,modorb,modpar,irststp,irstlay,nprfmo,
     *                nfmopal,modprp,maxl1c,ipieda,modgrd,modesp,ivmul,
     *                modlmo,nopden,mofock,modfd,modfmm,ncentm,ndualb,
     *                ngab,modpan
      data iargs/1,1,0,0,0/,c4dum/' '/
c
c     this subroutine reads FMO atomic coordinates.
c     parstat: GroupMasterBcast
c     iscan=0 normal reading and storing
c     iscan=1 scan to find NATFMO/maxbas, do not store fmozan, fmoc and izbas. 
c             In this case fmozan, fmoc and izbas should be declared in
c             the calling arguments as dum, dum3(3) and idum. 
c
      UNITS=1.0D+00/0.52917724924D+00
c     if(iand(nfmopal,1).ne.0) UNITS=1.889725989D+00
      if(iand(nfmopal,1).ne.0) UNITS=1.0D+00/0.529177249D+00
      if(iand(nfmopal,16).ne.0) UNITS=1.0D+00/0.52917720859D+00
      usexyz=iand(modio,512).ne.0
c
c     Try to open a file with coordinates. If it exists,
c     use it, else read from IR. 
c
      nftxyz=104
      if(usexyz) then
        CALL SEQOPN(nftxyz,'XYZ','UNKNOWN',.FALSE.,'FORMATTED')
      else
        nftxyz=ir
        CALL SEQREW(NFTXYZ)
        CALL FNDGRP(NFTXYZ,gprnam,IEOF)
        if(ieof.ne.0) then
          write(iw,9010) gprnam
          call abrt
        endif
      endif
      natfmo=0
      if(maswrk) then
  100   continue
          natfmo=natfmo+1
          natfmo0=natfmo
          if(iscan.ne.0) natfmo0=1
c         if(natfmo.gt.MXATMF) then
c           write(iw,*) 'Too many atoms',natfmo,MXATMF 
c           call abrt
c         endif
         atname='        '
         if(usexyz) then
           READ(NFTXYZ,9000,END=120,ERR=600) STR
c          No $end to look for. End if no more data found.
         else
           READ(NFTXYZ,9000,END=600,ERR=600) STR
           if(dolend(str)) goto 120
         endif
         call convsf77(80,str,strf77,5,iargs)
         READ(strf77,*,END=600,ERR=600) atname,ATOMNM,
     *                                  (fmoC(i,natfmo0),i=1,3)
c        READ(NFTXYZ,*,END=120,ERR=600) atname,ATOMNM,(fmoC(i,natfmo),i=1,3)
          call UPRCAS(atname,8)
          call UPRCAS(ATOMNM,8)
c         if(atname.eq.'$END') goto 120
          if(iuntrd.eq.1) then
            do i=1,3
              fmoC(i,natfmo0)=fmoC(i,natfmo0)*units
            enddo
          endif
          imode=0
c         write(6,*) 'Reading',natfmo
c         write(6,9999) atname,ATOMNM 
c9999     format(2A8)
          call zsymnum(ATOMNM,c4dum,imode)
          fmoZAN(natfmo0)=imode
c         write(6,*) 'wwwatm',imode,(fmoC(i,natfmo0),i=1,3)
          ipos=ifndchr(atname,8,'.')
          if(ipos.ne.0) then
            read(UNIT=atname(ipos+1:8),FMT='(I3)') ibas
            if(ibas.le.0) then
              write(6,*) 'Invalid basis number in $FMOXYZ',ibas
              call abrt
            endif
          else
            ibas=1
          endif
          if(ndualb.gt.0) ibas=2
          izbas(natfmo0)=ibas
          if(iscan.ne.0) maxbas=max(maxbas,ibas)
c         write(6,*) 'readbas',ibas,maxbas
        goto 100
      endif
  120 continue
      IF (GOPARR) CALL DDI_BCAST(2400,'I',natfmo,1,MASTER)
      natfmo=natfmo-1
      if(natfmo.le.0) then
        write(iw,*) 'No atomic coordinates found' 
        call abrt
      endif
      if(iscan.eq.0.and.GOPARR) then
        CALL DDI_BCAST(2401,'F',fmoC,3*natfmo,MASTER)
        CALL DDI_BCAST(2402,'F',fmoZAN,natfmo,MASTER)
        CALL DDI_BCAST(2403,'I',izbas,natfmo,MASTER)
      endif
      if(iscan.eq.1.and.GOPARR) then
        CALL DDI_BCAST(2401,'I',maxbas,1,MASTER)
      endif
      if(usexyz) CALL SEQCLO(nftxyz,'KEEP')
      return
  600 continue
c     call UPRCAS(atname,8)
c     if(atname.eq.'$END') goto 120
      write(6,*) 'error reading ',gprnam,', check atom',natfmo
c     now kill slave nodes by broadcasting reset natfmo
      natfmo=1
      CALL DDI_BCAST(2400,'I',natfmo,1,MASTER)
      call abrt
      return
 9000 format(A80)
 9010 format(1x,a,' input group not found')
      end
C
C*MODULE fmoio   *DECK lnkxyz
      SUBROUTINE lnkxyz(gprnam,nseq,npair,fmoc)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      Parameter (MAXLNK=100)
      character*99 STR,STRF77
      dimension iargs(3)
      character*8 gprnam
      dimension fmoC(3,*)
      logical dolend
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
c     COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      data iargs/0,0,0/
c
c     this subroutine reads LINK ATOM coordinates.
c
      UNITS=1.0D+00/0.52917724924D+00
      nattmp=0
  100   continue
          nattmp=nattmp+1
          if(nattmp.gt.MAXLNK) then
            write(iw,*) 'Too many link atoms',nattmp,MAXLNK 
            call abrt
          endif
         READ(IR,9000,END=600,ERR=600) STR
         if(dolend(str)) goto 120
         call convsf77(80,str,strf77,3,iargs)
         READ(strf77,*,END=600,ERR=600) idam1, idam2,
     *                      (fmoC(i,nattmp+nseq),i=1,3)
         if(1.lt.0) write(6,*) idam1, idam2
         do i=1,3
            fmoC(i,nattmp+nseq)=fmoC(i,nattmp+nseq)*units
         enddo
        goto 100
  120 continue
      nattmp=nattmp-1
      if(nattmp.ne.npair) then
        write(iw,9010) gprnam,nattmp,npair 
        call abrt
      endif
      return
  600 continue
c     call UPRCAS(atname,8)
c     if(atname.eq.'$END') goto 120
      write(6,*) 'error reading ',gprnam,', check atom',npair
c     now kill slave nodes by broadcasting reset natfmo
      call abrt
      return
 9000 format(A80)
 9010 format(1x,'The number of atoms in ',A8,' differs from Npair',2I8) 
      end
C
C*MODULE fmoio   *DECK fmolmo
      SUBROUTINE fmolmo(gprnam1,gprnam2,taotyp,iaprjo,japrjo,CoreAO,
     *                  OccCor,nCBS,nCAO,shiftb,libish,modQbas,nsp2hop,
     *                  nsp3shop)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      character*99 STR,STRF77
      dimension iargs(3),mCAO0(2)
      character*8 gprnam1,gprnam2,taotyp(*)
      logical GOPARR,DSKWRK,MASWRK,modQbas,dolend
      dimension iaprjo(MaxCAO,*),japrjo(MaxCAO,*),OccCor(MaxCAO,*),
     *          nCBS(*),nCAO(*),CoreAO(MaxCBS,MaxCAO,*),shiftb(MaxCAO,*)
     *         ,libish(*)
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /RUNOPT/ RUNTYP,EXETYP,NEVALS,NGLEVL,NHLEVL
      Common /fmoinf/ nfg,nlayer,natfmo,nbdfg,naotyp,nbody,nsegm
      common /fmoopt/ espsca(9),RESPAP(2),rESPPC(2),rESDIM,restri(4),
     *                rcorsd,respct,convfg,cnvdmp,coroff,rflmo(4),
     *                orshft,orshft2,cnvafo,ascreen(4),IXESP,mxitfg,
     *                nguess,NBSSE,modorb,modpar,irststp,irstlay,nprfmo,
     *                nfmopal,modprp,maxl1c,ipieda,modgrd,modesp,ivmul,
     *                modlmo,nopden,mofock,modfd,modfmm,ncentm,ndualb,
     *                ngab,modpan
      common /fmodim/ maxbnd,maxknd,maxcbs,maxcao,maxvec,maxl1,maxnat,
     *                maxabd,maxbas,maxbbd,maxlmo,maxslo,maxabd2,maxrij
      data iargs/1,0,0/,fmohmo/8HFMOHMO  /
c
c     this subroutine reads localised MOs for boundary orbitals.  
c     parstat: GroupMasterBcast
c     Two spellings of this group name are supported. 
C
c     orbdepB=iand(nguess,524288).eq.0
c     deltaB=10.0D+00
      CALL SEQREW(IR)
      CALL FNDGRP(IR,gprnam1,IEOF)
      if(ieof.ne.0) then 
        CALL SEQREW(IR)
        CALL FNDGRP(IR,gprnam2,IEOF)
        if(ieof.ne.0) return 
      endif 
      if(modQbas) call viclr(nCAO,1,maxknd*3)
      call viclr(nCBS,1,MaxKnd*2)
      ierr=0
      naotyp=0
      nsp2hop=0
      nsp3shop=0
      if(maswrk) then
  100   continue
          naotyp=naotyp+1
          if(naotyp.gt.MaxKnd) then
          write(iw,*) 'Too many LMO types (check $FMOHYB)',naotyp,MaxKnd
            call abrt
          endif
          ierr=0
          taotyp(naotyp)='        '
          READ(IR,9000,END=600,ERR=600) STR
          if(dolend(str)) goto 120
          call convsf77(80,str,strf77,3,iargs)
          READ(strf77,*,END=600,ERR=600) taotyp(naotyp),mCBS,mCAO
c         READ(IR,*,END=120,ERR=600) taotyp(naotyp),mCBS,mCAO
          call UPRCAS(taotyp(naotyp),8)
c         if(taotyp(naotyp).eq.'$END') goto 120
          if(mCBS.gt.MaxCBS) then
             write(iw,*) 'Too many AOs in $FMOHYB',mCBS,MaxCBS
             write(iw,*) 
     *           'The basis sets in $DATA and $FMOHYB are inconsistent.'
             call abrt
          endif
          if(mCAO.gt.MaxCAO) then
             write(iw,*) 'Increase maxcao in $fmo to',mCAO
             call abrt
          endif
c         The else clause covers some dormant option to be activated later.
          if(mod(nguess,2).ne.333) then
          READ(IR,*,END=600,ERR=600) (iaprjo(i,naotyp),japrjo(i,naotyp),
     *                          (CoreAO(j,i,naotyp),j=1,mCBS),i=1,mCAO)
          else
          READ(IR,*,END=600,ERR=600) (iaprjo(i,naotyp),japrjo(i,naotyp),
     *                               (CoreAO(j,i,naotyp),j=1,mCBS),
     *                                OccCor(i,naotyp),i=1,mCAO)
          endif
c         READ(IR,*,END=600,ERR=600) (iaotyp(j,naotyp),j=1,mCBS) 
          nCBS(naotyp)=mCBS
          nCAO(naotyp)=mCAO
          nCBS(naotyp+MaxKnd)=0
          if(taotyp(naotyp)(1:1).eq.'@') then
            nCBS(naotyp+MaxKnd)=1
            nsp2hop=nsp2hop+1
          endif
          if(taotyp(naotyp)(1:1).eq.'+') then
            nCBS(naotyp+MaxKnd)=2
            nsp3shop=nsp3shop+1
          endif
c         call dacopy(mCAO,orshft,shiftb(1,naotyp),1)
          do i=1,mCAO
            bi=orshft
c          if(iaprjo(i,naotyp).eq.2.or.japrjo(i,naotyp).eq.2) bi=1.0D+02
c           if(orbdepB) bi=bi+naotyp*deltaB
            shiftb(i,naotyp)=bi
c           if(iand(nguess,524288).ne.0.and.i.eq.mCAO) then
c             shiftb(i,naotyp)=orshft2
c             write(6,*) 'Setting B',i,orshft2
c           endif
          enddo
c
c         shiftb(1,naotyp)=1e+02
c
c        if(orbdepB) write(6,*) naotyp,'Bs=',(shiftb(i,naotyp),i=1,mCAO)
c         if(iand(ixesp,512).ne.0) call normhmo(iaprjo(1,naotyp),
          if(exetyp.eq.fmohmo) call normhmo(iaprjo(1,naotyp),
     *                    japrjo(1,naotyp),CoreAO(1,1,naotyp),mCBS,mCAO)
          if(modQbas) then
            ishiftl=naotyp+maxknd
            ishiftr=naotyp+maxknd*2
            mCAO0(1)=ishiftl
            mCAO0(2)=ishiftr
            call convlmo(mCBS,mCAO,mCAO0,naotyp,iaprjo(1,naotyp),
     *                   japrjo(1,naotyp),CoreAO,libish)
            nCAO(ishiftl)=mCAO0(1)
            nCAO(ishiftr)=mCAO0(2)
          endif
          ierr=1
          if(mCBS.lt.mCAO) goto 600
        goto 100
      endif
  120 continue
      IF (GOPARR) CALL DDI_BCAST(2403,'I',naotyp,1,MASTER)
      naotyp=naotyp-1 
      if(naotyp.le.0) then
        if(maswrk) write(iw,*) 'No LMOs found'
        if(naotyp.lt.0) then
          write(6,*) 'Invalid naotyp',naotyp
          call abrt
        endif
      endif
      IF (GOPARR) then
        CALL DDI_BCAST(2404,'F',taotyp,naotyp,MASTER)
        CALL DDI_BCAST(2405,'I',iaprjo,MaxCAO*naotyp,MASTER)
        CALL DDI_BCAST(2406,'I',japrjo,MaxCAO*naotyp,MASTER)
        naotyp3=naotyp
        if(modQbas) naotyp3=maxknd*3 
        CALL DDI_BCAST(2407,'F',CoreAO,MaxCBS*MaxCAO*naotyp3,MASTER)
        if(mod(nguess,2).eq.1) 
     *  CALL DDI_BCAST(2408,'F',OccCor,MaxCAO*naotyp,MASTER)
        CALL DDI_BCAST(2409,'I',nCBS,naotyp+MaxKnd,MASTER)
        CALL DDI_BCAST(2410,'I',nCAO,naotyp3,MASTER)
        CALL DDI_BCAST(2411,'F',shiftb,MaxCAO*naotyp,MASTER)
      endif
c     write(6,*) 'read ',naotyp,' LMO types'
      return
  600 continue
c     call UPRCAS(taotyp(naotyp),8)
c     if(taotyp(naotyp).eq.'$END') goto 120
      write(6,*) 'error reading LMO set ',naotyp,': check input'
      if(ierr.eq.0) write(6,*) 'format error (or insufficient data).' 
      if(ierr.eq.1) write(6,*) 'more MOs than AOs.'
c     if(ierr.eq.2) write(6,*) 'left-right assignment error.'
c     now kill slave nodes by broadcasting reset naotyp
      naotyp=0
      CALL DDI_BCAST(2403,'I',naotyp,1,MASTER)
      call abrt
      return
 9000 format(A80)
      end
C
C*MODULE fmoio   *DECK fmobon
      SUBROUTINE fmobon(gprnam,nhybnam,indat,taotyp,taotypi,iabdfg,
     *                  jabdfg,idxcao,doafo,nrothop)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      parameter (MAXL=5)
      character*99 STR,STRF77
      dimension iargs(2+MAXL*2+2)
      character*8 gprnam,taotyp(*),taotypi(*)
      logical GOPARR,DSKWRK,MASWRK,dolend,doafo
      logical mmonly,qmmm
      dimension indat(*),iabdfg(*),jabdfg(*),idxCAO(MaxBnd,*)
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      Common /fmoinf/ nfg,nlayer,natfmo,nbdfg,naotyp,nbody,nsegm
      common /fmodim/ maxbnd,maxknd,maxcbs,maxcao,maxvec,maxl1,maxnat,
     *                maxabd,maxbas,maxbbd,maxlmo,maxslo,maxabd2,maxrij
      COMMON /TINOPT/ mparti,MMONLY,QMMM
c
c     this subroutine reads boundary orbital indices.  
c     parstat: GroupMasterBcast
C
      nrothop=0
      CALL SEQREW(IR)
      CALL FNDGRP(IR,gprnam,IEOF)
      if(ieof.ne.0) return 
      iargs(1)=0
      iargs(2)=0
      do i=1,nhybnam
        iargs(2+i)=1
      enddo
      iargs(nhybnam+3)=0
      iargs(nhybnam+4)=0
      nbdfg=0
      if(maswrk) then
  100   continue
          nbdfg=nbdfg+1
          if(nbdfg.gt.MaxBnd) then 
             write(iw,*) 'Too many cut bonds',nbdfg,MaxBnd
             call abrt
          endif
          READ(IR,9000,END=600,ERR=600) STR
          if(dolend(str)) goto 120
          call convsf77(80,str,strf77,2+nhybnam+2,iargs)
          iabdfg(nbdfg)=0
          jabdfg(nbdfg)=0
          call viclr(idxcao(nbdfg,1),MaxBnd,nhybnam+2)
c         idxcao(nbdfg,nhybnam+1)=0 
c         idxcao(nbdfg,nhybnam+2)=0 
c         if(nhybnam.gt.0) then
c         READ(strf77,*,END=600,ERR=600) iabdfg(nbdfg),jabdfg(nbdfg),
c    *                                    (taotypi(i),i=1,nhybnam),
c         else
          do i=1,nhybnam
            taotypi(i)(1:1)=char(0)
          enddo
          READ(strf77,*,END=110,ERR=110) iabdfg(nbdfg),jabdfg(nbdfg),
     *        (taotypi(i),i=1,nhybnam),(idxcao(nbdfg,nhybnam+i),i=1,2)
  110 continue
c         write(6,*) 'wwwa',nhybnam,(idxcao(nbdfg,nlayer+i),i=1,2)
c         write(6,*) 'wwwb',strf77
c           read charge and multiplicity for FLMO boundaries. 
c           There is no error check now in the READ, so do it manually. 
            if(iabdfg(nbdfg).eq.0.or.jabdfg(nbdfg).eq.0) goto 600
c         endif
c         READ(IR,*,END=120,ERR=120) iabdfg(nbdfg),jabdfg(nbdfg),
c    *                               (taotypi(i),i=1,nhybnam)
          if(.not. qmmm) then
c         the check will be done later in QMMM calculation
          if(indat(abs(iabdfg(nbdfg))).eq.indat(jabdfg(nbdfg))) then
            if(maswrk) write(iw,9010) gprnam,nbdfg,iabdfg(nbdfg),
     *                                jabdfg(nbdfg),indat(jabdfg(nbdfg))
            call abrt
          endif 
          endif
          do i=1,nhybnam
            call UPRCAS(taotypi(i),8)
          enddo
c         reorder indices in the canonical way: first negative
          if(jabdfg(nbdfg).lt.0) then
             itmp=jabdfg(nbdfg)
             jabdfg(nbdfg)=iabdfg(nbdfg)
             iabdfg(nbdfg)=itmp
             if(maswrk) write(6,*) 'Warning: indices reordered',nbdfg
c            This is nothing important (is it?) 
          endif
          if(naotyp.ne.0) then
c            if(naotyp.lt.nhybnam) then
             do i=1,nhybnam
               if(taotypi(i)(1:1).eq.char(0)) then
                 if(maswrk) write(iw,9030) i,nbdfg
                 call abrt
               endif
               idxcao(nbdfg,i)=ifndtxt(naotyp,taotyp,taotypi(i))
             enddo
          else
             if(.not.doafo) then
               if(maswrk) write(iw,9020) nbdfg
               call abrt
             endif
          endif
          if(idxCAO(nbdfg,nhybnam+1).ne.0.and..not.doafo)
     *       nrothop=nrothop+1
        goto 100
      endif
  120 continue
      IF (GOPARR) then
         CALL DDI_BCAST(2410,'I',nbdfg,1,MASTER)
         CALL DDI_BCAST(2411,'I',iabdfg,nbdfg,MASTER)
         CALL DDI_BCAST(2412,'I',jabdfg,nbdfg,MASTER)
         nhybnam1=nhybnam+2
c        if(nhybnam.eq.0) nhybnam1=2
         CALL DDI_BCAST(2413,'I',idxcao,MaxBnd*nhybnam1,MASTER)
      endif
      nbdfg=nbdfg-1 
      if(nbdfg.le.0) then
        if(maswrk) write(iw,*) 'No bond cuts found.'
c       call abrt
      endif
      return
  600 continue
      write(6,*) 'error reading ',gprnam,' check entry ',nbdfg
      call abrt
      return
 9000 format(A80)
 9010 format(/1x,'Check',A8,', entree',I5,
     *          ': intrafragment fractioned bonds are not allowed.',
     *       /1x,2I5,' are in fragment',I5,/)
 9020 format(/1x,'No HMOs were read for bond',I5,
     *           ' - perhaps $FMOHYB is missing?') 
 9030 format(/1x,'HMO set',I2,' is not defined in line',I5,
     *           ' of $FMOBND.',
     *       /1x,'For multilayer runs, you need 1 set per layer.',/)
      end
c
C*MODULE fmoio   *DECK convsf77
      subroutine convsf77(lens,sinp,sout,nargs,iargs)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      character*99 sinp,sout
      character*1 quote,space
      dimension iargs(nargs) 
      data quote/''''/,space/' '/
c
c     Convert human string format into Fortran77 format by quoting all
c     string arguments listed in iargs. "99" above is obtained as 80+
c     at most 9 string arguments (since each adds 2 characters).
c     Only spaces are allowed to separate arguments (shame on tab users).
c     sinp: input string of length lens (<=80)
c     sout: output string
c     nargs: the expected number of space separated words in sinp
c     iargs(i) tells if word i is a string (iargs(i).ne.0).
c     nargs and iargs are arguments to convsf77.
c     Multiple space separation is allowed.
c
      do i=lens+1,99
        sout(i:i)=space
      enddo
      if(nargs.gt.9) call abrtx("Too many arguments read in convsf77")
      io=0
      iarg=0
      inside=0
      do i=1,lens
c
        if(sinp(i:i).eq.space.and.inside.ne.0) then
          inside=0
          if(iarg.le.nargs.and.iarg.gt.0) then
            if(iargs(iarg).ne.0) then
              io=io+1
              sout(io:io)=quote
c             end quote
            endif
          endif
        endif
c
        if(sinp(i:i).ne.space.and.inside.eq.0) then
          iarg=iarg+1
          inside=1
          if(iarg.le.nargs) then
            if(iargs(iarg).ne.0) then
              io=io+1
              sout(io:io)=quote
c             begin quote
            endif
          endif
        endif
c
        io=io+1
        sout(io:io)=sinp(i:i)
c       copy the rest, including possible NULL and garbage characters
      enddo
c     append end quote if needed at the very end
      if(inside.ne.0.and.iarg.le.nargs.and.iarg.gt.0) then
        if(iargs(iarg).ne.0) then
          io=io+1
          sout(io:io)=quote
c         end quote
        endif
      endif
      RETURN
      END
c
C*MODULE fmoio   *DECK ifndtxt
      integer function ifndtxt(n,text,texti) 
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      character*8 text(n),texti
c
      ifndtxt=0
      do i=1,n
         if(text(i).eq.texti) then
            ifndtxt=i
            return
         endif
      enddo
      write(6,9000) texti,n
 9000 format(1x,'TEXTI=',a,' not found in the stored area, n=',i4)
      call abrt
      return
      end
C
C*MODULE fmoio   *DECK zsymnum
      subroutine zsymnum(atninp,atnout,imode)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      parameter (MXEL=137)
      character*8 atninp,substring
      character*4 atnout,atmlab(MXEL)
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
     *            'LR  ','RF  ','HA  ','SG  ','NS  ','HS  ',
     *            'MT  ','110 ','111 ','112 ','113 ','114 ',
     *            '115 ','116 ','117 ','118 ','119 ','120 ',
     *            '121 ','122 ','123 ','124 ','125 ','126 ',
     *            '127 ','128 ','129 ','130 ','131 ','132 ',
     *            '133 ','134 ','135 ','136 ','137 '/
c
c     This subroutine either finds the atomic number given a text name (imode=0)
c     or returns the text name given an atomic number imode. 
C
      if(imode.eq.0) then
c       first try to see if atninp already contains an integer, possibly
c       followed by a decimal dot and some other numbers to be ignored.
        ipos=ifndchr(atninp,8,'.')
        if(ipos.eq.0) ipos=9
        substring=atninp(1:ipos-1)
        read(UNIT=substring,FMT='(I4)',ERR=100) iat
        imode=iat
        return
100     continue
        do i=1,MXEL
          if(atmlab(i).eq.atninp(1:4)) then
            imode=i
            return
          endif
        enddo
        write(6,*) 'Unknown element: ',atninp
        call abrt
      else
        atnout=atmlab(imode)
      endif
      return
      end
c
C*MODULE fmoio   *DECK fmoout
C>
C>     @brief FMO option summary
C>
C>     @details Print FMO option summary.
C>
C>     @author Dmitri Fedorov
C>
      subroutine fmoout(nder,need,lbody,fmozan,fmoscf,fmoci,fmodft,fmocc
     *                 ,mpnfmo,fmotd,modmol,molfrg,nstfmo,ncvscf,mcnvfmo
     *                 ,ngrfmo,loadbf,loadgr,nstjob,ichtot,multot,nacut,
     *                  savgrd,nxg,nyg,nzg,ngm,vdwrad,grdpad,spargrid,
     *                  nxyzg,ifgfmo0,iexcit,nafo,naoafo,ndmsiz,natbuf,
     *                 natprp,nactfg,gcorrel,fmosym,nfmosym,modcnv,docns
     *                 ,dodcesd,maxld,maxlt,fraggrid,nfraggr,nsubsys,
     *                  imect,dodos,nlcmo,IMEXFG,dofret,nfgfret,dofed,
     *                  dosczvprop,nbesp,mpcmit,doapc,mconfg,mconvex,
     *                  nsp2hop,nsp3shop,nrothop,fragnorm,IEFMORUN,dousc
     *                 ,doapcvac)
      USE DFTBPB_MOD,ONLY: PERIOD
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      parameter (one=1.0D+00,maxpst=10)  
      logical savgrd,spargrid,gcorrel,QFMM,QOPS,docns,dodcesd,fraggrid,
     *        dodos,dofret,dofed,DFTBFL,SCC,SRSCC,DFTB3,DAMPXH,LCDFTB,
     *        dosczvprop,doapc,dousc,doapcvac
      character*1 nauf
      character*2 sfock,symafo
      character*3 symfd
      character*80 fmosym
      character*8 monostr,dens,symscr,espsym,solvdip,scrnorm,apcasc
      COMMON /DFTB  / DFTBFL,SCC,SRSCC,DFTB3,DAMPXH,LCDFTB
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON /PCMITR/ RCUT(2),THRES,IPCMIT,IMUL,MXDIIS,NREG,MXITR1,
     *                MXITR2,MODPAP
      COMMON /potopt/ pcmGaussEx,pcmGaussTol,pcmGaussExSqrt,pcmGaussRCut
     *               ,modpot
      COMMON /PCMPAR/ IPCM,NFT26,NFT27,IKREP,IEF,IP_F,nfmopcm,IHET
      COMMON /QMFM  / SIZE,EPS,DPGD,QFMM,NP,NS,IWS,NPGP,MPMTHD,NUMRD,
     *                ITERMS,QOPS,ISCUT
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
      common /fmodim/ maxbnd,maxknd,maxcbs,maxcao,maxvec,maxl1,maxnat,
     *                maxabd,maxbas,maxbbd,maxlmo,maxslo,maxabd2,maxrij
      common /fmopmd/ fmobox(3),mdwpbc,nimgcell,imglvl,
     *                ltrvec,lfmogctr,lfmoctmp,lindatmd,lwrkdsav,lindxiu
     *               ,IPBCFST
      Common /fmoinf/ nfg,nlayer,natfmo,nbdfg,naotyp,nbody,nsegm
      dimension lbody(*),fmozan(*),fmoscf(*),fmoci(*),fmodft(*),fmocc(*)
     *         ,mpnfmo(*),fmotd(*),molfrg(*),ncvscf(2),mcnvfmo(*),
     *          ngrfmo(maxpst,*),loadbf(maxpst,*),loadgr(maxpst,*),
     *          nstjob(*),ngm(3),vdwrad(*),grdpad(3),iexcit(6),maxln(3),
     *          nlcmo(3),mconfg(*)
c
      write(iw,8500)
      if(ifgfmo0.le.1) write(iw,8510)
      write(iw,8520)
      monostr='internal'
      if(iand(modcnv,1).ne.0) monostr='embedded'
      maxln(1)=maxl1
      maxln(2)=maxld
      maxln(3)=maxlt
      write(iw,9000) need,nfg,nlayer,nbody,natfmo,ichtot,multot,maxbas,
     *               naotyp,nbdfg,convfg,mxitfg,monostr,
     *               (maxln(i),i=1,nbody)
      nespap=0
      if(iand(ndualb,2).ne.0) nespap=1
      if(doapc) then
        write(iw,8982)
      else
        if(ndualb.ne.0.and.iand(ndualb,4).eq.0) write(iw,8980) nespap
        if(iand(ndualb,4).ne.0) write(iw,8981)
      endif
      if(IMEXFG.ne.0) write(iw,8985) IMEXFG
      if(nbdfg.ne.0) write(iw,8990) orshft
      if(ifgfmo0.gt.0) write(iw,9005) ifgfmo0
      nmols=0
      do ifg=1,nfg
c       if(mulfg(ifg).ne.1) write(iw,8070) ifg,mulfg(ifg)
c       if(scffrg(ifg).ne.scftyp) write(iw,8060) ifg,scffrg(ifg)
c       if(molfrg(ifg).ne.0) nmols=ifg
        if(molfrg(ifg).gt.0) nmols=ifg
      enddo
      nsca=0 
      do i=1,9
        if(espsca(i).ne.one) nsca=i
      enddo
      write(iw,8000)
      write(iw,8002) (lbody(i),i=1,nlayer)
      write(iw,8005) (fmoscf(i),i=1,nlayer)
      write(iw,8010) (mpnfmo(i),i=1,nlayer)
      write(iw,8020) (fmodft(i),i=1,nlayer)
      write(iw,8030) (fmoci(i),i=1,nlayer)
      write(iw,8040) (fmocc(i),i=1,nlayer)
      write(iw,8050) (fmotd(i),i=1,nlayer)
c     do ilay=1,nlayer
c     if(nfmopcm.ne.0.and.mpcmit.ne.2) write(iw,8080) nfmopcm
c     if(nfmopcm.ne.0.and.mpcmit.eq.2) write(iw,8082) 1,nfmopcm
      write(iw,*) 'FMO method: ',fmosym(1:min(nfmosym,72))
      espsym='[1]'
      if(nbesp.eq.2.and.mpcmit.eq.2) espsym='[1(2)]'
      if(nbesp.eq.2.and.mpcmit.gt.2) espsym='[2]'
      if(nbesp.eq.2) write(iw,*) 'Embedding: ESP',espsym
c     write(iw,8080) fmosym(1:min(nfmosym,72))
      if(iexcit(1).gt.0)  write(iw,8090) (iexcit(i),i=1,4) 
      if(dofret) write(iw,8092) nfgfret
      if(dofed) write(iw,8094) nfgfret
c
      if(nbdfg.ne.0) then
        if(rflmo(1).ne.0) then
          symafo(1:1)=' '
          if((IEFMORUN.GT.0.OR.gcorrel).AND.IAND(MODLMO,128).EQ.0)
     *      symafo(1:1)='G'
          symafo(2:2)='0'
          if(nafo.eq.1) symafo(2:2)='H' 
          if(nafo.eq.2) symafo(2:2)='1'
          write(iw,8110) symafo(1:1),symafo(2:2)
          if(nafo.ge.2) write(iw,8112) cnvafo
          if(iand(modlmo,4).eq.0) write(iw,8120)
          if(iand(modlmo,4).ne.0) write(iw,8121)
          if(iand(modlmo,8).eq.0) write(iw,8122)
          if(iand(modlmo,8).ne.0) write(iw,8123)
          if(iand(modlmo,32).eq.0) write(iw,8125) 'partial'
          if(iand(modlmo,32).ne.0) write(iw,8125) 'full'
          if(naoafo.gt.0) write(iw,8130) naoafo 
          if(rflmo(4).ne.0) write(iw,8135) rflmo(4) 
        else 
          write(iw,8115)
        endif
      endif
      write(iw,8140)
      if(respap(1).gt.0) write(iw,9010) respap(1)
      if(respap(1).lt.0) write(iw,9011) 
      if(respap(2).gt.0) write(iw,9015) respap(2) 
      if(respap(2).lt.0) write(iw,9016) 
      if(resppc(1).gt.0) write(iw,9020) resppc(1)
      if(resppc(1).lt.0) write(iw,9021) 
      if(resppc(2).gt.0) write(iw,9025) resppc(2)
      if(resppc(2).lt.0) write(iw,9026) 
      if(nbody.gt.1) then  
        if(resdim.gt.0) write(iw,9030) resdim
        if(resdim.lt.0) write(iw,9031) 
        if(rcorsd.gt.0) write(iw,9035) rcorsd
        if(rcorsd.lt.0) write(iw,9036) 
        if(nbody.gt.2) then  
         if(restri(1).gt.0) write(iw,9032) restri(1),restri(2),restri(3)
         if(restri(4).gt.0) write(iw,9033) restri(4)
        endif
      endif
      if(iand(modESP,4096).ne.0) write(iw,9037) 
      if(iand(IXESP,2).ne.0) write(iw,9091) 
      if(iand(IXESP,32).ne.0) write(iw,9095)
      if(iand(modESP,7).eq.0) write(iw,9096)
      if(iand(modESP,7).eq.1) write(iw,9097) 
      if(iand(modESP,7).eq.2) write(iw,9098) 
      if(iand(modESP,8).ne.0) write(iw,9094) 
c     if(iand(modESP,64).ne.0) write(iw,9098) 
c     if(ifgdon.ne.0.and.ifgacc.ne.0) write(iw,9097) ifgdon,ifgacc
      call vdwrout(fmozan,vdwrad) 
      if(nder.gt.0) then
        if(iand(modgrd,1).ne.0) write(iw,9040) 
        if(iand(modgrd,2).ne.0) write(iw,9042) 
c       if(iand(modgrd,4).ne.0) write(iw,9043) 
        if(iand(modgrd,8).ne.0) write(iw,9044) 
        if(iand(modgrd,16).eq.0.and.nbdfg.ne.0) write(iw,9045) 
        if(iand(modgrd,32).ne.0) write(iw,9046) 
        if(iand(modgrd,64).ne.0) write(iw,9047)
      endif
      if(modfd.ne.0) then
        symfd='FD '
        if(iand(modfd,2).ne.0) symfd='FDD'
        if(iand(modfd,4).ne.0) symfd='FDB'
        write(iw,9050) symfd,nactfg,natbuf,natprp
      endif
      if(iand(modfmm,1).ne.0) write(iw,9180)
      if(iand(modfmm,2).ne.0) write(iw,9181)
c     if(iand(modfmm,4).ne.0) write(iw,9182)
      if(iand(modfmm,8).ne.0) write(iw,9183)
      if(modfmm.ne.0) write(iw,9195) NPGP/2
      if(dodcesd) write(iw,9197)
c
      write(iw,8150)
      if(mod(nguess/32,2).eq.1.and.dftbfl) then
         write(iw,9078)
      else
        if(mod(nguess,  2).eq.1) write(iw,9070) 
        if(iand(nguess,2).ne.0) write(iw,9071) 
        if(iand(nguess,4).ne.0) write(iw,9072) 
        if(iand(nguess,8).ne.0) write(iw,9073) 
        if(iand(nguess,32).ne.0) write(iw,9075)
        nauf='1'
        if(iand(nguess,64).ne.0) nauf='n'
        write(iw,9076) nauf
      endif
      if(mod(nguess/16,2).eq.1) write(iw,9074)
      if(mod(nguess/128,2).eq.1) write(iw,9077)
      if(iand(nguess,512).ne.0) write(iw,9079)
      if(iand(nguess,1024).eq.0) write(iw,9080)
      if(iand(nguess,32768).ne.0) write(iw,9081)
c     if(iand(nguess,262144).ne.0) write(iw,9076) 
      if(irststp.ne.0) write(iw,9125) irstlay,irststp,modrst
      if(savgrd) write(iw,9126)
c
      write(iw,8160)
      if(nbsse.eq.1) write(iw,9085) 
      if(nbsse.eq.2) write(iw,9086) 
      if(nbsse.eq.3) write(iw,9087) 
      if(iand(modESP,32).ne.0) write(iw,9090) 
      if(nrothop.ne.0) then
        if(nsp2hop.eq.0.and.nsp3shop.eq.0) write(iw,9092)
        if(nsp3shop.ne.0) write(iw,9088) nsp3shop 
        if(nsp2hop.ne.0) write(iw,9089) nsp2hop 
      endif
      if(iand(nfmopal,1).ne.0) write(iw,9150)
      if(iand(nfmopal,8).ne.0) write(iw,9153)
      if(ascreen(1).ne.0) write(iw,9160) (ascreen(i),i=1,2)
      if(nsubsys.ne.0) then
          write(iw,9175) nsubsys
      else if(nmols.ne.0) then
        if(iand(modmol,1).ne.0) then
          write(iw,9171) 
        else
          if(gcorrel) write(iw,9170) 
        endif
        if(iand(modmol,2).ne.0) then
          write(iw,9172) (molfrg(i),i=1,nmols)
        else
          write(iw,9173) (molfrg(i),i=1,nmols)
        endif
        if(iand(modmol,4).ne.0.and.nbody.gt.2) write(iw,9174)
      endif
      if(docns) write(iw,9093)
c
      write(iw,8170)
      if(nacut.ne.0) write(iw,9300) nacut 
c
      write(iw,8180)
      nmanscf=0
      do i=1,mconvex
        if(mconfg(nfg+(i-1)*4+1).gt.0) nmanscf=nmanscf+1
      enddo
      write(iw,9001) (mcnvfmo(i),i=1,nstfmo) 
      if(nmanscf.gt.0) write(iw,9007) nmanscf
      if(ncvscf(1).ne.9999) write(iw,9002) ncvscf(1) 
      if(ncvscf(2).ne.9999) write(iw,9003) ncvscf(2) 
      if(mod(modorb,  2).eq.0) write(iw,9100) 
      if(mod(modorb,  2).ne.0) write(iw,9101) 
      if(mod(modorb/2,2).ne.0) write(iw,9102) 
      if(nsca.gt.0) write(iw,9004) (espsca(i),i=1,nsca)
      if(coroff.ne.0) write(iw,9006) coroff
      if(iand(modESP,128).ne.0) write(iw,9099) 
c
      write(iw,8190)
      if(mod(modpar,  2).ne.0) write(iw,9110)
      if(mod(modpar/2,2).eq.0) write(iw,9111)
      if(mod(modpar/2,2).ne.0) write(iw,9112)
      if(iand(modpar,4).ne.0) write(iw,9113)
      if(iand(modpar,8).ne.0) write(iw,9114)
      if(iand(modpar,16).ne.0) write(iw,9115)
      if(iand(modpar,32).ne.0) write(iw,9109)
      if(iand(modpar,64).ne.0) write(iw,9116)
      if(iand(modpar,128).ne.0) write(iw,9117)
      if(iand(modpar,512).ne.0)  write(iw,9118) 'supervector',ndmsiz
      if(iand(modpar,1024).ne.0) write(iw,9118) '     matrix',ndmsiz
      if(iand(modpar,512+1024).ne.0.and.iand(modpar,2048).eq.0) 
     *  write(iw,9119)
      if(iand(modpar,8192).ne.0) write(iw,9120)
      dens='electron'
      if(iand(modprp,128).ne.0) dens='  spin  '
      do i=1,nlayer  
        if(ngrfmo(1,i).ne.0) write(iw,9130) i,ngrfmo(1,i),ngrfmo(4,i),
     *                                        ngrfmo(6,i)
        if(ngrfmo(2,i).ne.0) write(iw,9132) i,ngrfmo(2,i),ngrfmo(5,i),
     *                                      ngrfmo(7,i)
        if(ngrfmo(3,i).ne.0) write(iw,9134) i,ngrfmo(3,i)
        if(nstjob(1).ne.0.and.i.eq.1) write(iw,9140) i,nstjob(1),
     *                                ' monomer',loadbf(1,i),loadgr(1,i)
c       For I>1 the number is not computed correctly, but the load
c       balancing works.
      enddo
c
      write(iw,8200) imect
      if(iand(nprfmo,3).eq.0) write(iw,9200) 'Voluminous'
      if(iand(nprfmo,3).eq.1) write(iw,9200) '   Compact'
      if(iand(nprfmo,3).eq.2) write(iw,9200) '    Frugal'
      if(iand(nprfmo,3).eq.3) write(iw,9200) '      Puny'
      if(iand(nprfmo,128).ne.0) write(iw,9202)
      if(iand(nprfmo,4).ne.0) write(iw,9206) 
      if(iand(nprfmo,8).ne.0) write(iw,9208) 
      if(iand(modprp,1).ne.0) write(iw,9210)
      if(iand(modprp,4).ne.0) write(iw,9216) dens,nxg,nyg,nzg 
      if(iand(modprp,32).ne.0) write(iw,9217)
      if(iand(modprp,32).ne.0.and.nfmopcm.ne.0) 
     *  write(iw,9218) pcmGaussEx,pcmGaussTol
      if(iand(modprp,512).ne.0) write(iw,9219)
      if(iand(modprp,4096).ne.0) write(iw,9241)
      if(iand(modprp,8192).ne.0) write(iw,9242)
      if(iand(modprp,65536).ne.0) write(iw,9244)
      if(dftbfl.and.period.and.iand(modprp,16384).eq.0) write(iw,9243)
      if(spargrid) write(iw,9220) (nxg*nyg*nzg*1.0D+02)/nxyzg,nxyzg
      ngmm=ngm(1)*ngm(2)*ngm(3)
      if(ngmm.ne.0) write(iw,9222) ngm(1),ngm(2),ngm(3)
      if(ngmm.ne.0.or.fraggrid) write(iw,9223) nfraggr
      if(ngmm.ne.0.or.iand(modprp,4).ne.0.or.spargrid) 
     *  write(iw,9230) grdpad(1),grdpad(2) 
      if(fraggrid.and.grdpad(3).ne.0) write(iw,9231) grdpad(3)
      if(ipieda.eq.2) write(iw,9240) 
      if(ipieda.eq.1) write(iw,9250)
      symscr='local'
      if(iand(modpap,8).ne.0) symscr='partial'
      scrnorm='uniform'
      if(fragnorm.ne.0) scrnorm='fragment' 
      if(nfmopcm.gt.0) write(iw,9255) symscr
      if(nfmopcm.gt.0.and.iand(modpap,8).ne.0) write(iw,9257) scrnorm
      apcasc='specific'
      if(dousc) apcasc='uniform'
      if(doapcvac) apcasc='semidry'
      if(ndualb.ne.0.and.nfmopcm.gt.0) write(iw,9258) apcasc
      solvdip='local'
      if(iand(modpap,8).ne.0.and.iand(modprp,262144).eq.0) 
     *   solvdip='partial'
      if(nfmopcm.gt.0) write(iw,9256) solvdip
      if(iand(mofock,1).ne.0) then
        sfock='F '
        if(iand(modesp,32).ne.0) sfock='XF'
        if(iand(modesp,32).eq.0.and.iand(mofock,2).ne.0) sfock='FX'
        write(iw,9260) sfock
        if(iand(mofock,4).ne.0) write(iw,9261)
        if(iand(mofock,8).eq.0) write(iw,9262)
        if(iand(mofock,16).ne.0) write(iw,9263)
        if(iand(mofock,32).ne.0) write(iw,9264)
      endif
      if(iand(nlcmo(1),1).ne.0) then
        sfock='  '
        if(iand(nlcmo(1),2).ne.0) sfock='X '
        write(iw,9290) sfock,(nlcmo(i),i=1,3)
      endif
      if(mdwpbc.ne.0) then
        write(iw,8300)
         toangs = 0.52917724924d+00
         write(iw,9271) (fmobox(ixyz)*toangs,ixyz=1,3)
         write(iw,9272) nimgcell,2*imglvl+1
      end if
      if(dodos) write(iw,9280)
      if(dosczvprop) write(iw,9310)
c
      write(iw,9500) maxbnd,maxknd,maxcbs,maxcao,maxbbd,maxnat
      write(iw,9900) 
c
 8000 format(1x,'Layer electron correlation information')
 8002 format(4x,' nbody=',8(1x,I2,5x))
 8005 format(4x,'SCFTYP=',8A8)
 8010 format(4x,'MPLEVL=',8(1x,I2,5x))
 8020 format(4x,'DFTTYP=',8A8)
 8030 format(4x,' CITYP=',8A8)
 8040 format(4x,' CCTYP=',8A8/)
 8050 format(4x,' TDDFT=',8A8/)
c8060 format(4x,'Fragment ',I5,' has SCFTYP=',A8)
c8070 format(4x,'Fragment ',I5,' has MUL=',I4)
c8080 format(1x,'FMO method: ',A72)
 8090 format(1x,'Excited state options',
     *       /5x,'Excited fragment:         ',I5,
     *       /5x,'Many-body excitation level:  ',I2,
     *       /5x,'Energy option:               ',I2, 
     *       /5x,'Dimer matching scheme:       ',I2) 
 8092 format(1x,'States in Foerster resonance energy transfer ',
     *       '(FRET):',I8)
 8094 format(1x,'Fragment excitation difference (FED) set:',I8)
 8110 format(/1x,'Bond fragmentation scheme:',A1,'AFO',A1)
 8112 format(4x,'AFO1 conergence threshold:',E12.4)
 8115 format(/1x,'Bond fragmentation scheme: HOP') 
 8120 format(4x,'One AFO per bond frozen.')
 8121 format(4x,'All AFO in a bond frozen.')
 8122 format(4x,'Weak localisation.')
 8123 format(4x,'Strong localisation.')
 8125 format(4x,'Using ',A7,' Mulliken localization criterion.')
 8130 format(4x,'Max AO size for caching AFO density:',I6)
 8135 format(4x,'In the localization criterion, cut primitives with ',
     *          'exponents larger than',F10.2)
 8140 format(/1x,'FMO approximations.')
 8150 format(/1x,'FMO initial guess options.')
 8160 format(/1x,'FMO extensions.')
 8170 format(/1x,'FMO automation options.')
 8180 format(/1x,'FMO SCF converger options.')
 8190 format(/1x,'FMO parallelisation options.')
 8200 format(/1x,'FMO properties and print-out level.',
     *       /4x,'Interfragment CT method:',I5)
 8300 format(/1x,'FMO MD with Periodic Boundary Condition (PBC)')
 8500 format(/1x,71(1H-),
     *      /18x,'The Fragment Molecular Orbital (FMO) method.')
 8510 format(25x,'Version 5.5',
     *       /1x,'Contributors: N. Asada, K. R. Brorsen, ',
     *           'M. Chiba, C. Choi, D. G. Fedorov,',
     *       /1x,'D. S. Kaliakin, N. Minezawa, V. Mironov, T. Nagata, ',
     *           'H. Nakata, Y. Nishimoto,',
     *       /15x,'S. R. Pruitt, C. Steinmann, and F. Zahariev.',
     *       /7x,'All publications using the FMO method in GAMESS ',
     *           'should cite:'
     *       /10x,'D. G. Fedorov, K. Kitaura, J. Chem. Phys. ',
     *           '120, 6832 (2004).')
 8520 format(1x,71(1H-))
 9000 format(/1x,'Used memory (words): ',I12,
     *       /1x,'Number of fragments:',I7,/1x,'Number of layers:',I10,
     *       /1x,'N-body FMO method:  ',I7,/1x,'Number of atoms: ',I10,
     *       /1x,'Total charge: ',I13,
     *       /1x,'Total spin multiplicity:',I3,
     *       /1x,'Number of basis sets:',I6,
     *       /1x,'Number of LMO types:',I7,/1x,'Number of boundaries:',
     *    I6,/1x,'MonomerSCF convergence:',E13.4,
     *       /1x,'Max monomerSCF iter:',I7,
     *       /1x,'MonomerSCF criterion:     ',A8,
     *       /1x,'Max AOs per n-mer:',3I9)
 8980 format( 1x,'Dual basis:',13x,'AP',I1)
 8981 format( 1x,'Dual basis:',13x,'AE')
 8982 format( 1x,'Dual basis:',13x,'APC')
c8980 format( 1x,'Dual basis:',A23)
 8985 format( 1x,'MEX fragment:',I14)
 8990 format( 1x,'Projection operator constant:',E10.2)
 9001 format(4x,'SCF convergers:',/16I7/)
 9002 format(4x,'SCF convergers will be switched after',I3,
     *          ' monomer SCF iterations.')
 9003 format(4x,'SCF convergers will be reset to mconv(2) after',I3,
     *          ' monomer SCF iterations.')
 9004 format(4x,'SCF ESP scaling at initial iterations:',9F5.2/)
 9005 format(1x,'Free monomer:    ',I9)
 9006 format(4x,'DFT will be turned off until monomer SCF converges to',
     *          E11.4)
 9007 format(4x,'Manually set FMO convergers for n-mers:',I3)
 9010 format(4x,'Using AO pop approximation for the 2e ESP with the VdW'
     *         ,' factor',F7.2)
 9011 format(4x,'Using AO pop approximation for the 2e ESP of all ',
     *           'fragments.')
 9015 format(4x,'Using AO pop approximation for separated ',
     *          'dimers with the VdW factor',F7.2)
 9016 format(4x,'Using AO pop approximation for the energy of all ',
     *          'dimers .')
 9020 format(4x,'Using point charge approximation for the 2e ESP ',
     *           'with the VdW factor',F7.2)
 9021 format(4x,'Using point charge approximation for the 2e ESP ',
     *           'of all fragments.')
 9025 format(4x,'Using point charge approximation for the energy of ',
     *           'of separated dimers with the VdW factor',F7.2)
 9026 format(4x,'Using point charge approximation for the energy of ',
     *           'all dimers.')
 9030 format(4x,'Using ES-dimer approximation with the VdW factor',F7.2)
 9031 format(4x,'Using ES-dimer approximation for all fragments.')
 9032 format(4x,'Ignore trimers using thresholds T1, T2 and T3 if ',
     *      /4x,'RD>T1 and RM>T2 or RM>T3',3F7.2,/4x,'where RM is the ',
     *       'monomer-dimer distance and the monomer is chosen so that',
     *      /4x,'the remaining dimer has the shortest interfragment',
     *      /4x,'distance RD among all 3 dimers in the trimer.') 
 9033 format(4x,'Ignore correlation for trimers separated by the VdW ',
     *          'factor',F7.2)
 9035 format(4x,'No dynamic correlation for dimers with the VdW',
     *          ' factor ',F7.2)  
 9036 format(4x,'No dynamic correlation for all dimers.')
 9037 format(4x,'No ES dimer contributions.')
 9040 format(4x,'FMO gradient: retain ESP in the Lagrangian.')
 9042 format(4x,'FMO gradient: add ESP derivative terms.')
 9044 format(4x,'FMO gradient: add Mulliken charge derivative.')
 9045 format(4x,'FMO gradient: add HOP derivatives.')
 9046 format(4x,'FMO gradient: add response terms with SCZV.')
 9047 format(4x,'FMO gradient: reverse gradient projection.')
 9050 format(4x,'Optimisation: ',A3,' with',I4,' fragments and',I5,
     *          ' atoms in buffer (prop:',I5,' )')
 9070 format(4x,'Using free molecule guess for monomers.')
 9071 format(4x,'Using monomer density to form n-mer guess.') 
 9072 format(4x,'Adding HOP projector into Huckel guess.')
 9073 format(4x,'Using dimer electronic state from previous geometry.')
 9074 format(4x,'Run RHF to get DFT or MCSCF initial dimer orbitals.')
 9075 format(4x,'Compute MCSCF canonical orbitals based on the additive' 
     *      /4x,'monomer densities and use such orbitals for dimers.') 
c9076 format(4x,'Applying HOP projector to the dimer initial orbitals.')
c9076 format(4x,'Project out detached orbitals from initial MO space.')
 9076 format(4x,'Use ',A1,'-mer Aufbau for initial n-mer orbitals.')
 9077 format(4x,'Orbitals will not be reused during geometry ',
     *          'optimisation.')
 9078 format(4x,'Skip initial n-mer density and MOs.')
 9079 format(4x,'MCSCF orbitals will be manually reordered.')
 9080 format(4x,'Orthonormalize initial guess orbitals.')
 9081 format(4x,'Use orbitals from lower levels as initial guess.')
 9085 format(4x,'BSSE correction: many body counter poise, ',
     *          'not self-consistent.')
 9086 format(4x,'BSSE correction: many body ES counter poise, ',
     *          'not self-consistent.')
 9087 format(4x,'BSSE correction: counter poise in vaccuum, ',
     *          'not self-consistent.')
 9088 format(4x,'HOP extensions: DR(sp3o), sets',I3)
 9089 format(4x,'HOP extensions: DR(sp2), sp2 sets:',I3)
 9090 format(4x,'Adding exchange terms to ESP (ESPX).')
 9091 format(4x,'No 2nd order ESP energy terms.')
 9092 format(4x,'HOP extensions: DR(sp3)')
 9093 format(4x,'Adding CNS terms.')
 9094 format(4x,'No embedding potential (ESP) will be used.')
 9095 format(4x,'No gradient contribution from separated dimers.') 
 9096 format(4x,'The uniform distances will be used in ESPs.')
 9097 format(4x,'The consistent distances will be used in ESPs for ',
     *          'unconnected n-mers.')
 9098 format(4x,'The consistent distances will be used in ESPs for ',
     *          'all n-mers.')
 9099 format(4x,'Use FDIFF-like SCC accelerator.')
 9100 format(4x,'Fragment densities will be exchanged.')
 9101 format(4x,'Fragment orbitals will be exchanged.')
 9102 format(4x,'Fragment energies will be exchanged.')
 9109 format(4x,'Memory usage will be reduced.')
 9110 format(4x,'Smart parallel load balancing: large jobs first.')
 9111 format(4x,'ESPs will be shell parallelised.')
 9112 format(4x,'ESPs will be fragment parallelised.')
 9113 format(4x,'Broadcast all fragments at once.')
 9114 format(4x,'Do the fragment initialisation in the reverse order.')
 9115 format(4x,'Do not parallelise projector operator contributions.')
 9116 format(4x,'Initial density file will be broadcast from grand ',
     *          'master.')
 9117 format(4x,'Rearrange separated dimers.')
 9118 format(4x,'Store fragment data in DDI memory as ',A11,
     *          ', total size=',I12)
 9119 format(8x,'(limit integroup I/O to masters).') 
 9120 format(4x,'Use DDI memory for screening.')
 9125 format(4x,'Restart layer: ',I2,' step: ',I2,' record ',I9)
 9126 format(4x,'Gradient restart information will be dumped to F38.') 
 9130 format(4x,'Layer',I2,', monomer(SCF,MP2,MCSCF) groups are:',3I5)
 9132 format(4x,'Layer',I2,',   dimer(SCF,ESD,MCSCF) groups are:',3I5)
 9134 format(4x,'Layer',I2,',  trimer(SCF) groups are:          ',I5)
 9140 format(/1x,'Layer',I2,' static LB for',I6,A8,'(s) larger than ',
     *           I5,' on',I5,' group(s).',/)
 9150 format(4x,'G94-like interfragment distances: rounded off to .1A.')
 9153 format(4x,'QMCPack data will be generated.')
 9160 format(4x,'Point charge screening parameters:',2F10.4)
 9170 format(4x,'Correlation will be applied only to selected n-mers.')
 9171 format(4x,'Only selected n-mers will be computed (n>1).') 
 9172 format(4x,'Intramolecular model (within selected).',
     *      /4x,'The following monomers are selected: ',99I4)
 9173 format(4x,'Intermolecular model (between selected and the rest).',
     *   /4x,'Total interaction energies will be printed for monomers:',
     *          99I4) 
 9174 format(4x,'Economy mode: no trimer arrays.') 
 9175 format(4x,'Analysis will be printed for subsystems:',I5) 
 9180 format(4x,'Multipoles will be used to compute ES dimers.')
 9181 format(4x,'Multipoles will be used to sum over ES dimers.')
c9182 format(4x,'Multipoles will be used to compute ES dimers.')
 9183 format(4x,'Multipoles will be used for 1e ESP gradients.') 
 9195 format(8x,'The highest angular momentum is:',I3)
 9197 format(4x,'Adding dispersion to ES dimers.')
 9200 format(4x,A10,' print-out level.')
 9202 format(4x,'Skip ES dimer printout.')
 9206 format(4x,'Fragment connexion information will be printed.')
 9208 format(4x,'Atomic charges will be printed.')
 9210 format(4x,'The total electron density in AO basis will be dumped',
     *          ' to F10.')
 9216 format(4x,'The total ',A8,' density will be punched on the grid:',
     *          3I6)
 9217 format(4x,'The molecular electrostatic potential will also be ',
     *           'punched.')
 9218 FORMAT(/1X,'Using Gaussians for PCM contribution to MEP, params:',
     *           2E12.5)
 9219 format(4x,'Global matrices will be stored in DDI memory.')
 9241 format(4x,'Dispersion will not be divided into fragments.')
 9242 format(4x,'Dispersion will be computed for ES fragments.')
 9243 format(4x,'Atoms will be internally wrapped in one cell.')
 9244 format(4x,'Exchange will be computed for ES fragments.')
 9220 format(4x,'Sparse grid takes',F6.1,'% of the full grid: (',I10,
     *          ' words).')
 9222 format(4x,'Monomer grid will have dimensions of',3I6)
 9223 format(4x,'Fragments to compute data on grid',I7)
 9230 format(4x,'Total box padding parameter in vdW radii is',F9.3,
     *      /4x,'n-mer box padding parameter in vdW radii is',E9.2)
 9231 format(4x,'Interfragment threshold for grid',F7.2)
 9240 format(4x,'Pair interaction energy decomposition analysis (PIEDA)'
     *         ,' will be done (full).')
 9250 format(4x,'Pair interaction energy decomposition analysis (PIEDA)'
     *         ,' will be done (PL).')
 9255 format(4x,'Solvent screening model: ',A8)
 9256 format(4x,'Solvent dipole model: ',A8)
 9257 format(4x,'Solvent charge normalization: ',A8)
 9258 format(4x,'Solvent charges in APC: ',A8)
 9260 format(4x,'Data for the Fock matrix construction',
     *          ' will be punched: FMO/',A2)
 9261 format(7x,'Kinetic and 1e integrals will be punched.')
 9262 format(7x,'RESPPC will be raised by 2 for the total Fock matrix.')
 9263 format(7x,'Compute overlaps for ES dimers.')
 9264 format(7x,'Use the internal formulation.')
 9271 format(4x,'Box length in Angstrom:',3F12.5)
 9272 format(4x,'Number of image cells: ',i5,' (=',i2,'**3-1)')
 9280 format(4x,'DOS will be computed.')
 9290 format(4x,'The Fock matrix construction: FMO/LCMO',A2,
     *      /7x,'Options:',I12,' ; use occ/virt orbitals:',2I8)
 9300 format(4x,'Number of atoms per fragment',I5)
 9310 format(4x,'Corrections from relaxed density will be added.')
 9500 format(/1x,'Array dimensions maxbnd,maxknd,maxcbs,maxcao,maxbbd,',
     *        'maxnat are:',/17x,6I7)
 9900 format(/1x,50(1H-)//)
      return
      END
c
C*MODULE fmoio   *DECK fmoprop
C>
C>    @brief evaluates total FMO properties
C>
C>    @author Dmitri Fedorov
C>
C>    @date 25/03/13 - Casper Steinmann
C>     - Added EFMO/FD and initial support for EFMO/PCM
C>
      SUBROUTINE fmoprop(nder,ichfmo,nefmo,mulfmo,l0fmo,l1fmo,nfg2d,nfg2
     *                  ,nfg3,nenm,nend,nent,needdm,ichfg,frgnam,layfrg,
     *                   indat,iabdfg,jabdfg,itrlay,fmozan,scffrg,fmoscf
     *                  ,lcorrel,modmol,exclmol,molfrg,emolfrg,emon,edim
     *                  ,etrim,edimq,esolv,emocdr,fmode,fmoq,primul,
     *                   enucfmo,e1efmo,ekinfmo,atmulq,nfmoelm,fmoelm,
     *                   nelm,modef0,gcorrel,nedimes,nedimex,ifgfmo0,
     *                   doddcor,isgddi0,extracc,nextracc,didcc,didmp,
     *                   etotdft,nedft,ctdim,n0bda,i0bda,e0bda,efmo0,
     *                   epl0ds,eint0,dopleda,rappri,nappri,erapp,nrapp,
     *                   ext2lay,DEFTF,TORQF,dolat,nunint,nsymeq,iexcit,
     *                   eexcit,texcit,nstmono,osmd,isumd,isumt,eexfg,
     *                   modcha,e0centr,dopdc,dotd,doci,doeom,ipeam,
     *                   fzcor,totfock,numfrg,m1efmo,eigfmo,nspins,ibfmo
     *                  ,skipesd,urohf,douhf,dodc,fmohard,eaip,iactfg,
     *                   iactat,skipscc1,
     *                   nevsav,prtdst,ascat,esdi,hasgrad,fullmfmo3,
     *                   edimlow,needmd,docns,cnsdat,ALCNT2F,runtyps,
     *                   dodcesd,dopbcmd,savemem,etrimsum,reducee,ctspin
     *                  ,pmulspin,savemem2,edimsum,urospn,subsys,nsubsys
     *                  ,subprp,dofret,domipea,excit2,excit3,texcit2,
     *                   mdoutmin,dodos,ignoresd,extri,exfid,dodft,
     *                   ecorrdft,do3c,fast3loop,loop3ij,ibloop,edimmixd
     *                 ,fast3prop,loadt,net,dadb,hfretap,nfgfret,ifgfret
     *                  ,popmat,dofed,autosa,isasign,nfgasub,esatot,
     *                   samon,sadim,dosczvprop,sczvprp,iapsign,epcmap,
     *                   dofret2,szfmo,doapc)
      USE mx_limits, only: mxfrg,mxdfg,mxatm,mxrt,mxgrid
      USE comm_FRGINF
      USE comm_EFPFMO
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      logical primul,masout,needdm,exclmol,usedij,gcorrel,doddcor,GOPARR
     *       ,DSKWRK,MASWRK,isgddi,parout,INITGDDI,isgddi0,didcc,didmp,
     *        skipc,dopleda,ext2lay,fmoq,dolat,dopdc,primul2,SG1,dotd,
     *        doci,doeom,doexc,totfock,skipesd,dovir,urohf,douhf,dodc,
     *        mixdim,skipscc1,pcmprp,convSCC,wasgddi,hasgrad,fullmfmo3,
     *        MLGDDI,lowtrimer,docns,dodcesd,dopbcmd,savemem,didijk,
     *        reducee,urospn,savemem2,subsys,scfnone,dofret,dofed,
     *        mdoutmin,dodos,ignoresd,dodft,domdan,outgrd,outgrdn,
     *        fast3loop,fast3prop,autosa,autosa1,dosczvprop,domipea,
     *        dofret2,DFTBFL,SCC,SRSCC,DFTB3,DAMPXH,LCDFTB,do3c,dispdft,
     *        didmp2,doapc
      character*1 corri,corrs(4),symbe
      character*3 corstr2,apl0
      character*4 Scavcds
      character*6 corstr(3),esstr,solstr,dispstr,aapl0
      character*9 charger
      character*11 bpl0
      character*25 distr
      parameter(MXPTPT=100,sawdust=1.0D-25,
     *          TOKCAL=627.51D+00,toeV=27.21138386D+00,one=1.0D+00,
     *          zero=0.0D+00,units=0.52917724924D+00)
c    *          zero=0.0D+00,tokcals=627.509541D+00)
      common /cnsdat/ EXREF,lcnsdat,natcns,ioover
      COMMON /CDSPRT/ GCDS,AREACDS
      COMMON /DFGRID/ DFTTHR,DFTGTHR,SWOFF,SW0,BSLRD(137),NDFTFG,
     *                NRAD,NTHE,NPHI,NRAD0,NTHE0,NPHI0,
     *                NANGPT(MXGRID),NANGPT0(MXGRID),SG1,JANS
      COMMON /DFTB  / DFTBFL,SCC,SRSCC,DFTB3,DAMPXH,LCDFTB
      COMMON /EFMO  / EFMOETOT,EFMOESERG,EFMOPOLERG,EFMODISERG,
     *                EFMOREPNRG,EFMOCHTNRG,EFMOEPEN,EFMOPCMG,
     *                IEFMORUN,IMODEFP,IEFMORT,IEFMOCFRG,IEFMONFRG,
     *                IEFMODIM(2),IMODEFE,NATEFMO,IMODEFD,IMODEFER,
     *                IMODEFCT,IDIMTYP,iefmo_agrad
      COMMON /ELPROP/ ELDLOC,ELMLOC,ELPLOC,ELFLOC,
     *                IEDEN,IEMOM,IEPOT,IEFLD,MODENS,
     *                IEDOUT,IEMOUT,IEPOUT,IEFOUT,
     *                IEDINT,IEMINT,IEPINT,IEFINT
      COMMON /ENRGYS/ ENUCR,EELCT,ETOT,SZ,SZZ,ECORE,ESCF,EERD,E1,E2,
     *                VEN,VEE,EPOT,EKIN,ESTATE(MXRT),STATN,EDFT(3),EDISP
      COMMON /FMCOM / X(1)
      COMMON /FUNCT / E,EG(3,MXATM)
      COMMON /GDDI/   ISCOPE,NGROUPS,MYGROUP,MEGLOB,npglob,nnglob,JBTYP,
     *                ISGDDI,PAROUT,INITGDDI,wasgddi,MLGDDI,NSUBGR,
     *                MeUniv,NPUniv,numdlb,myworld,nworlds,mogddi
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON /LMOEDA/ GDTOLA,GJTOLA,GKTOLA,TKTOLA,VTOLA,
     *                GDTOLB,GJTOLB,GKTOLB,TKTOLB,VTOLB,
     *                ECORL,EXCOR
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /PCMOPT/ RABI,RASC,REFPOL,THRSLS,DENSLS,WB,WA,ETA2,GD,EVAC,
     *                RHOW,PM,AREATL,AREAKP,BONDRY,OMEGA,RET,FRO,EPSINF,
     *                EPS,DR,RSOLV,VMOL,TCE,STEN,DSTEN,CMF,TABS,IDIRCT,
     *                IPCDER,IDP,ICOMP,IFIELD,ICAV,IDISP,IPRINT,IRETCAV,
     *                ICENT,IFAST,NEVAL,IEFPOL,KEEPSM,IMGABI,IMGASC,NADD
      COMMON /PCMPAR/ IPCM,NFT26,NFT27,IKREP,IEF,IP_F,nfmopcm,IHET
      COMMON /PCMPLY/ STOT,VOL,CCX,CCY,CCZ,RDIF,NESF,NESFP,NESFF,I_NESF,
     *                L_AST
      COMMON /PCMPRT/ GCAVP,GCAVS,GDISP,GRP,EHFGAS
      COMMON /POINTS/ NPOINT,IPUNIT,XPOINT(MXPTPT),YPOINT(MXPTPT),
     *                ZPOINT(MXPTPT)
      COMMON /PRPOPT/ ETOLLZ,ILOCAL,IAHARD
      COMMON /SOLSMX/ SOLA,SOLB,SOLC,SOLG,SOLH,SOLN,ISMX
      COMMON /STNBUF/ STNPNT(4,2*MXATM),BIGEXP,NPTSTN,NBUFFM
      COMMON /XYZPRP/ XP,YP,ZP,
     *                DMX,DMY,DMZ,
     *                QXX,QYY,QZZ,QXY,QXZ,QYZ,
     *                QMXX,QMYY,QMZZ,QMXY,QMXZ,QMYZ,
     *                OXXX,OXXY,OXXZ,OXYY,OYYY,OYYZ,
     *                OXZZ,OYZZ,OZZZ,OXYZ,
     *                OMXXX,OMXXY,OMXXZ,OMXYY,OMYYY,
     *                OMYYZ,OMXZZ,OMYZZ,OMZZZ,OMXYZ
      COMMON /RUNOPT/ RUNTYP,EXETYP,NEVALS,NGLEVL,NHLEVL
      common /fmodim/ maxbnd,maxknd,maxcbs,maxcao,maxvec,maxl1,maxnat,
     *                maxabd,maxbas,maxbbd,maxlmo,maxslo,maxabd2,maxrij
      Common /fmoinf/ nfg,nlayer,natfmo,nbdfg,naotyp,nbody,nsegm
      COMMON /FMOMD/  efmomd(10),MODFMO,lemonmd,ledimmd,lrminmd,lindmd,
     *                lemonmds,ledimmds,letrimmd
      COMMON /MCPFMO/ IMCPFMO,LFZCOR,LIFMPTYP,LIFMPTYP2,LMCPSW,
     *                LIZCOR2,imp0,jmp0,icorsh0,igtf0,IECPFMO
      common /fmoopt/ espsca(9),RESPAP(2),rESPPC(2),rESDIM,restri(4),
     *                rcorsd,respct,convfg,cnvdmp,coroff,rflmo(4),
     *                orshft,orshft2,cnvafo,ascreen(4),IXESP,mxitfg,
     *                nguess,NBSSE,modorb,modpar,irststp,irstlay,nprfmo,
     *                nfmopal,modprp,maxl1c,ipieda,modgrd,modesp,ivmul,
     *                modlmo,nopden,mofock,modfd,modfmm,ncentm,ndualb,
     *                ngab,modpan
      common /fmopmd/ fmobox(3),mdwpbc,nimgcell,imglvl,
     *                ltrvec,lfmogctr,lfmoctmp,lindatmd,lwrkdsav,lindxiu
     *               ,IPBCFST
      common /fmorun/ espscf,e0scf(2),emp2s,IDAFMO,icurfg,jcurfg,kcurfg,
     *                icurlay,icurunt,nat1e,ncursh,ngau,icurpop,ifmostp,
     *                moncor,needr,modrst,norbproj,nunesp,iskipesp,
     *                IESDPPC,idoprop,mp2run,icurit,idmfmo,iddfmo,
     *                iddcur,nddleft,ivmfmo,nzmtfmo,ifmobas,itmfmo(2)
      COMMON /MLTFMO/ Q_MUL,IMLTFMO,ISWPFD,ISWNEW1,lfmoicm
      dimension ichfg(*),frgnam(*),layfrg(*),indat(*),iabdfg(*),
     *          jabdfg(*),itrlay(*),fmozan(*),scffrg(*),fmoscf(*),
     *          lcorrel(*),molfrg(*),emolfrg(nfg,*),emon(nfg,4,*),
     *          edim(nfg2d,*),etrim(nfg3,*),edimq(*),esolv(*),
     *          emocdr(*),fmode(3,natfmo,*),atmulq(natfmo,*),hfretap(*),
     *         fmoelm(nfmoelm,*),nelm(0:*),dumbuf(3),extracc(nextracc,*)
     *         ,etotdft(*),ctdim(*),ctmax(4),ictmax(4),e0bda(5,nbdfg,*),
     *          efmo0(nfg,2,*),epl0ds(nfg,3),eint0(4),rappri(3),
     *          erapp(nappri,2),nrapp(nappri,2),DEFTF(3,NFRG,*),
     *          TORQF(3,NFRG,*),nsymeq(*),iexcit(6),eexcit(mxrt,3),
     *          texcit(3,mxrt,2),osmd(mxrt,2),isumd(*),isumt(*),
     *          eexfg(2,nfg,mxrt),fzcor(*),numfrg(*),eigfmo(*),ibfmo(*),
     *          fmohard(natfmo,*),eaip(nfg,2),iactfg(*),iactat(*)
     *         ,prtdst(4),ascat(3,*),t(3),esdi(*),edimlow(nfg2,3),
     *          cnsdat(3,*),ALCNT2F(*),etrimsum(7),ctspin(*),
     *          pmulspin(natfmo,*),edimsum(16),subprp(*),excit2(*),
     *          excit3(*),texcit2(*),extri(*),ecorrdft(*),loop3ij(*),
     *         edimmixd(nfg2,7),loadt(*),ifgfret(*),popmat(maxnat,nfg,*)
     *         ,samon(2,*),sadim(*),sczvprp(*),epcmap(3),szfmo(3),
     *          dsoltot(3)
      data corrs/'C','N','S','M'/
      DATA RMC/8HMCSCF   /,check/8HCHECK   /,optfmo/8HOPTFMO  /,
     *     rohf/8HROHF    /,uhf/8HUHF     /,rnone/8HNONE    /,
     *     FMOMD /8HMD      /,rungrad/8HGRADIENT/
      masout=.false.
      convSCC=.false.
c
c     Compute properties for FMO.
c     parstat: GlobalNone 
c     Group masters do redundant (each group does the same) calculation
c     of properties and broadcast the results onto slaves.
c     urospn = urohf .and. nder.eq.0
c
c     ifmostp=6
      domdan=runefp.eq.fmomd.and.iand(modfmo,1).ne.0
      modef=modef0
c     a quick fix to avoid printing the gradient in MD...
      if(runefp.eq.fmomd) modef=1
      doexc=dotd.or.doci.or.doeom.or.ipeam.ne.0
      dispdft=dodc.or.gcorrel
      didmp2=didmp.and..not.dodft
      symbe='e'
      if(maswrk) then 
      if(mdoutmin) then
        call timit(1)
        write(iw,9100)
      endif
c     force output of properties for PIEDA runs
      if(isgddi0) then 
        masout=(iand(nprfmo,3).le.1.or.ipieda.gt.0).and.meglob.eq.0
        outgrd=(runefp.eq.rungrad.or.iand(nprfmo,3).le.1).and.
     *          meglob.eq.0
      else
        masout=iand(nprfmo,3).le.1.or.ipieda.gt.0.and.maswrk
        outgrd=(runefp.eq.rungrad.or.iand(nprfmo,3).le.1).and.maswrk
      endif
      primul2=primul
      if(iand(modcha,3).eq.1) primul2=.false.
      usedij=iand(IXESP,2).eq.0
      corstr(1)='uncorr' 
      apl0='DI'
      aapl0='disp0'
      bpl0='dispersion '
      if(dodft) then
        apl0='RC'
        aapl0='rc0  '
        bpl0='correlation'
      endif
      if(gcorrel) then
        if(doexc) then
          corstr(1)='exc   '
          corstr(2)='ground'
          corstr(3)='excit '
        else
          corstr(1)='corr  '
          corstr(2)='uncorr'
          corstr(3)='delta'
        endif
      else if(dodc) then
        corstr(1)='unco+D'
        corstr(2)='uncorr'
        corstr(3)='disp'
        if(do3c) corstr(3)='di+bs '
      endif
      if(nfmopcm.ne.0) then
        esstr(1:3)=corstr(1)(1:3)
        esstr(4:6)='_es'
        solstr(1:3)=corstr(1)(1:3)
        solstr(4:6)='+so'
        if(corstr(1).eq.'unco+D') then
          esstr(1:3)='unD'
          solstr(1:3)='unD'
        endif
      endif
      sz0=(mulfmo-1)/2.0D+00
      s2_0=sz0*(sz0+1)
c     TDDFT/PCM does not compute some n-mers for iexcit(3)=0, thus PCM
c     properties are not complete, and are normally suppressed.
c     At present multilayer properties are not available (one layer only).
c     pcmprp=nfmopcm.ne.0.and.
c    *       (iexcit(1).eq.0.or.iexcit(3).ne.0.or.iexcit(5).ne.0)
      pcmprp=nfmopcm.ne.0.and.
     *       .not.(nbody.lt.1.or.doexc.and.iexcit(3).eq.0)
      if(ignoresd.or.dofret.or.dofed) enucfmo=0
c     This option does not compute this property correctly.
c
c     for RHF-D:MP2 the label will be "corr", not "unco+D".
      ndiv2=0 
      ndiv3=0
c
      convSCC=.true.
      do ilay=1,nlayer
        iter=itrlay(ilay)
        if(iter.le.mxitfg) then
          if(mdoutmin) then
          if(iter.eq.0) then
            write(iw,9106) ilay
          else
            write(iw,9109) ilay,iter,mxitfg,convfg
          endif
          endif
        else
          convSCC=.false.
          write(iw,9108) ilay,mxitfg,convfg
        endif
      enddo
      corstr2=' un'
      if(doddcor) corstr2='   '
      if(mdoutmin) then
      write(iw,9220) enucfmo,natfmo,nefmo,ichfmo,mulfmo,l1fmo,l0fmo
      write(iw,9110) 
      endif
      if(nfmoelm.gt.0.and.mdoutmin) write(iw,9111) XPOINT(npoint)*units,
     *                        YPOINT(npoint)*units,ZPOINT(npoint)*units
      if(gcorrel.and.mdoutmin) write(iw,9112) corstr2
      if(didcc.and.mdoutmin) write(iw,9113)
c     if(pcmprp) then
      if(nfmopcm.ne.0) then
        if(masout.and.pcmprp) write(iw,9350)
        if(masout.and.pcmprp.and.iand(modprp,64).ne.0) write(iw,9352)
        ges=zero
        IF (IMLTFMO.EQ.1) GES = GES + Q_MUL
      endif
      if(masout.and.doapc) write(iw,9353)
      ENUCR=enucfmo
      emixdim=0
      rmixdim=1.0D+30
      imixdim=0
      jmixdim=0
c
      if(nbody.lt.1.or.doexc.and.iexcit(3).eq.0) goto 1000 
      if(dofret) then
        etotu=0
        goto 1000
      endif
c
      if(masout) then
      write(iw,9114)
c       if(ipieda.ne.0) then
        if(ipieda.gt.1) then
          if(dopleda) then
            write(iw,8117) apl0
          else
            if(nenm.gt.1.or.dodft) then
              write(iw,8115) apl0,apl0
            else
              write(iw,8116)
            endif
          endif
        else
          if(nenm.gt.1) write(iw,9115)
          if(nenm.le.1.and.iemom.ge.1) write(iw,9116)
        endif
      endif
      edftb0=0
      nmon0=1
      if(dodc) nmon0=2
      emontot=zero
      emontotu=zero
      e0tot=zero
      e0totu=zero
      if(nfmoelm.ne.0) call vclr(fmoelm(1,nfg+3),1,nfmoelm) 
      epl0s=zero
      epl0d=zero
      epld=zero
      epl0di=zero
      epld_di=zero
      epld_pld=zero
      ilay0=1
      if(modfd.ne.0) ilay0=2
      do ifg=1,nfg
        if(masout) then
        do ilay=ilay0,layfrg(ifg)
          nelm1=nelm(1)
          if(iemom.lt.1.or.ilay.lt.layfrg(ifg)) nelm1=0
          scfnone=scffrg(ifg).eq.rnone
c         if(ipieda.ne.0) then
          if(ipieda.gt.1) then
            eipl0d=epl0ds(ifg,1)
            eipl0s=epl0ds(ifg,2)
            eipl0=eipl0d+eipl0s
            if(ilay.eq.layfrg(ifg)) e0tot=e0tot+efmo0(ifg,1,ilay)
            if(gcorrel.or.dodc) then
              eipld=emon(ifg,2,ilay)-efmo0(ifg,2,ilay)
              eipl0di=epl0ds(ifg,3)
              eipld_di=emon(ifg,1,ilay)-efmo0(ifg,1,ilay)-eipld-eipl0di
              if(ilay.eq.layfrg(ifg)) e0totu=e0totu+efmo0(ifg,2,ilay)
            else if(dodft) then
              eipl0di=epl0ds(ifg,3)
              eipld_di=ecorrdft(ifg)-efmo0(ifg,2,ilay)
              eipld=emon(ifg,1,ilay)-efmo0(ifg,1,ilay)-eipld_di
              eipld_di=eipld_di-eipl0di
              if(ilay.eq.layfrg(ifg)) e0totu=e0totu+efmo0(ifg,2,ilay)
            else
              eipld=emon(ifg,1,ilay)-efmo0(ifg,1,ilay)
              eipl0di=zero
              eipld_di=zero
            endif
            eipld_pld=eipld-eipl0d
            epl0s=epl0s+eipl0s
            epl0d=epl0d+eipl0d
            epl0di=epl0di+eipl0di
            epld_di=epld_di+eipld_di
            epld_pld=epld_pld+eipld_pld 
            epld=epld+eipld
            if(.not.scfnone) then
            if(dopleda) then
              write(iw,8118) ifg,frgnam(ifg),ilay,eipl0d*tokcal,
     *                       eipl0s*tokcal,eipl0*tokcal,eipl0di*tokcal,
     *                       (fmoelm(ielm,ifg),ielm=nelm(0)+1,nelm1)
            else
            if(nenm.gt.1.or.dodft) then
              write(iw,9117) ifg,frgnam(ifg),ilay,eipl0d*tokcal,
     *                       eipl0s*tokcal,eipl0*tokcal,eipl0di*tokcal,
     *                       eipld_pld*tokcal,eipld_di*tokcal,
     *                       (fmoelm(ielm,ifg),ielm=nelm(0)+1,nelm1)
            else
              write(iw,9119) ifg,frgnam(ifg),ilay,eipl0d*tokcal,
     *                      eipl0s*tokcal,eipl0*tokcal,eipld_pld*tokcal,
     *                       (fmoelm(ielm,ifg),ielm=nelm(0)+1,nelm1)
            endif
            endif
            endif
          else 
            if(.not.scfnone) then
            if(nenm.gt.1) then
              write(iw,9118) ifg,frgnam(ifg),ilay,(emon(ifg,k,ilay),k=1,
     *                     nenm),(fmoelm(ielm,ifg),ielm=nelm(0)+1,nelm1)
            else
              write(iw,9120) ifg,frgnam(ifg),ilay,(emon(ifg,k,ilay),k=1,
     *                     nenm),(fmoelm(ielm,ifg),ielm=nelm(0)+1,nelm1)
            endif
            endif
          endif 
        enddo
          if(iemom.ge.2) 
     *      write(iw,8122) (fmoelm(ielm,ifg),ielm=nelm(1)+1,nelm(2))
          if(iemom.ge.3) 
     *      write(iw,8123) (fmoelm(ielm,ifg),ielm=nelm(2)+1,nelm(3))
        endif
        if(domdan) x(lemonmd+ifg-1)=emon(ifg,1,layfrg(ifg))
        ilay=layfrg(ifg)
        if(modfd.eq.0.or.ilay.ge.2) then
          emontot=emontot+emon(ifg,1,ilay)
          if(nenm.ge.2) emontotu=emontotu+emon(ifg,2,ilay)
          if(nfmoelm.ne.0) 
     *      call daxpy(nfmoelm,one,fmoelm(1,ifg),1,fmoelm(1,nfg+3),1)
        endif
        if(dftbfl) edftb0=edftb0+emon(ifg,nmon0,ilay)-ecorrdft(ifg)
c       write(6,*) 'wwwaaa',ifg,emon(ifg,nmon0,ilay)
      enddo
      if(dftbfl.and.masout) write(iw,9215) edftb0
c     if(nfmopcm.ne.0.and.iemom.ge.1)
c    *  call pcmmul(qse,XPOINT(npoint),YPOINT(npoint),ZPOINT(npoint),
c    *              ichfg,fmoelm(1,1))
      if(masout) then
        if(dodft.and.ipieda.ne.0) then 
          write(6,*) 'Monomer correlation energies in DFT (au)'
          do ifg=1,nfg
            write(iw,9120) ifg,frgnam(ifg),ilay,ecorrdft(ifg)
          enddo
        endif
        write(iw,8150) 
        do ifg=1,nfg
          if((modfd.eq.0.or.layfrg(ifg).ge.2).and.scffrg(ifg).ne.rnone)
     *    then
            vea=eaip(ifg,1)*toeV
            vip=eaip(ifg,2)*toeV
            vmu=-(vip+vea)/2
            veta=vip-vea
            if(veta.ne.0) then
              vs=1.0d+00/veta
              vomega=vmu*vmu/veta/2.0d+00
            else
              vs=0
              vomega=0
            endif
            write(iw,8160) ifg,vip,vea,vmu,veta,vs,vomega
          endif
        enddo
        call stabanal(layfrg,scffrg,eaip)
      endif 
c
      IF (IEFPFMO.NE.0) THEN
        emontot  = emontot  + REPNUCEFP
        emontotu = emontotu + REPNUCEFP
      END IF
c     IF (IMLTFMO.EQ.1) THEN
c       emontot  = emontot  + Q_MUL
c       emontotu = emontotu + Q_MUL
c     END IF
      etot=emontot
      etotu=emontotu
      IF (IEFPFMO.EQ.2) THEN
        WRITE(IW,8170)
        DO IFG = 1, NFG
          WRITE(IW,'(I6,F16.8)') IFG, ALCNT2F(IFG)
        END DO
      END IF
c
      if(ifgfmo0.ne.0) then
c       Save some data for FMO0, based on just fragment 1.
        ETOT=emon(1,1,layfrg(1))
        ESCF=emon(1,2,layfrg(1))
        if(dodft.or.dftbfl) ECORL=ecorrdft(1)
c       This is a despicable way to store ES moments. 
        if(nfmoelm.ne.0) call dcopy(nfmoelm,fmoelm(1,1),1,DMX,1)
      endif
      if(ipieda.eq.2) then
        write(iw,9200) symbe,corstr(1),0,e0tot
        if(gcorrel.or.dodc) then
          write(iw,9200) symbe,'uncorr',0,e0totu
          write(iw,9200) symbe,'delta ',0,e0tot-e0totu
        endif
      endif
      if(mdoutmin) write(iw,9200) symbe,corstr(1),1,emontot
c     for .not.gcorrel the data are not accumulated for uncorr. 
      if(nenm.ge.2.and.(gcorrel.or.dodc)) then
        if(mdoutmin) write(iw,9200) symbe,corstr(2),1,emontotu
        if(mdoutmin) write(iw,9200) symbe,corstr(3),1,emontot-emontotu
        if(didcc) call fmoccp(1,emontotu,extracc(1,1))
        if(didmp2) call fmompp(1,emontotu,extracc(1,1))
        if(do3c) call fmoh3c(1,extracc(1,1))
c       It is possible to implement this print-out for PCM (easy),
c       but the present code does not subtract the PCM potential properly.
      endif
      szfmo(2)=szfmo(2)+szfmo(1)
      szfmo(3)=szfmo(3)+szfmo(2)
      if(douhf.and.mdoutmin) 
     *  write(iw,9780) 1,sz0,1,szfmo(1)*(szfmo(1)+1),s2_0
c     if(ifgdon.ne.0.and.ifgacc.ne.0) then
c        ilayd=layfrg(ifgdon)
c        ilayc=layfrg(ifgacc)
c        edonacc=emon(ifgdon,1,ilayd)+emon(ifgacc,1,ilayc)
c        write(iw,9420) corstr2,edonacc
c        if(nenm.gt.1) then
c          edonacc=emon(ifgdon,2,ilayd)+emon(ifgacc,2,ilayc)
c          write(iw,9420) 'uncorr',edonacc
c        endif
c     endif
      if(iand(modcha,4).ne.0) then
        iscf=1
        if(gcorrel.or.dodc) iscf=2
        epl0dc=(emon(1,iscf,nlayer)-e0centr)*tokcal
        write(iw,9460) epl0dc,-epl0dc*2,-epl0dc
c       Assuming linear response
      endif
      if(iemom.ge.1)
     *  write(iw,9230) 1,(fmoelm(ielm,nfg+3),ielm=nelm(0)+1,nelm(1)),
     *                  DNRM2(3,fmoelm(1,nfg+3),1)
      if(iemom.ge.2) 
     *  write(iw,9231) (fmoelm(ielm,nfg+3),ielm=nelm(1)+1,nelm(2))
      if(iemom.ge.3) 
     *  write(iw,9232) (fmoelm(ielm,nfg+3),ielm=nelm(2)+1,nelm(3))
      if(ifgfmo0.eq.0.and.nfmoelm.ne.0) 
     *  call dcopy(nfmoelm,fmoelm(1,nfg+3),1,DMX,1)
c
c     if(nder.gt.0.and.modef.ne.1 .or.modfd.ne.0) then
      if(nder.gt.0.and.modef.ne.1.and.(runtyps.ne.optfmo.or.nbody.eq.1))
     *then
        if(hasgrad) then
          outgrdn=outgrd.and.(nbody.eq.1.or.iand(nprfmo,3).le.1)
          call fmogout(1,outgrdn,indat,fmozan,fmode,DEFTF,TORQF)
        else
          write(iw,8113)
        endif
      endif
c     call fmogout(1,.true.,indat,fmozan,fmode,DEFTF,TORQF)
      if(m1efmo.ne.0) call printeig(nefmo,ichfmo,mulfmo,numfrg,eigfmo,
     *                              ibfmo,dodos,m1efmo,nspins)
c
      gsolv=0
      gesfact=1
c     IF(pcmprp) then
      IF(nfmopcm.ne.0) then
        if(iswnew1.eq.1) gesfact=2
        Scavcds='Gcav'
        IF(ISMX.NE.0) Scavcds='Gcds'
        if(nbody.ge.1) then
          if(masout) write(iw,9360) Scavcds
          stotes=ddot(nfg,emocdr(6),10,one,0)
          gcavs=0
          do 610  ifg=1,nfg
            if(modfd.ne.0) then
              if(nlayer.ne.layfrg(ifg)) then
                go to 610
              end if
            end if
            gesi = ESOLV(IFG)
            gesi2= gesi-ESOLV(IFG+nfg)
c           if(masout) write(iw,9365) ifg,frgnam(ifg),esolv(ifg)*tokcal,
            gsolvi=gesi2*tokcal
            Qqi=emocdr((ifg-1)*10+7)+ichfg(ifg)
            if(abs(Qqi).gt.1.0D-06) then
              eps_eff=ichfg(ifg)/Qqi
            else
              eps_eff=0
c             to avoid printing NaNs or huge numbers. 
c             0 is chosen here as it may be easier to relay to humans
c             that something is fishy. Also, 0/0 should be resolved as 0.
            endif
            do 600 i=1,3
  600         gsolvi=gsolvi+emocdr((ifg-1)*10+i)
            if(domdan) x(lemonmd+ifg-1)=x(lemonmd+ifg-1)+gsolvi/tokcal
            if(masout) write(iw,9365) ifg,frgnam(ifg),
     *                     (emocdr((ifg-1)*10+i),i=4,6),
     *                      emocdr((ifg-1)*10+6)/stotes*1.0d+02,
     *                      emocdr((ifg-1)*10+7),eps_eff,
     *                     (emocdr((ifg-1)*10+i),i=8,10),
     *                  gesi2*tokcal,(emocdr((ifg-1)*10+i),i=1,3),gsolvi
C           ges=ges+esolv(ifg)
            ges=ges+gesi
            gsolv=gsolv+gsolvi
            gcavs=gcavs+emocdr((ifg-1)*10+1)
            ifgsa=ifg+nfgasub 
            if(autosa) 
     *       samon(2,ifgsa)=samon(2,ifgsa)+isasign*iapsign*gsolvi/tokcal
c           The signs are: isasign is from SA; iapsign is from AP.
            if(ibloop.lt.3) epcmap(1)=epcmap(1)+iapsign*gsolvi/tokcal
c           Delta-solvation energy (BS1,V)-(BS1,0).
  610     continue
          if(I_NESF.ne.natfmo.and.masout.and.icav.eq.1) then
            gcav0=gcavp-gcavs
            gsolv=gsolv+gcav0
            write(iw,9366) 0,'dummy',zero,gcav0,zero,zero,gcav0
c           Only Ecav is non-zero.
          endif
          if(masout.and.mdoutmin) write(iw,9355) 1,gsolv,1,ges*tokcal
          if(masout.and.iand(nprfmo,3).eq.0.and.nfmoelm.ne.0) then
            write(iw,9357)
            do ifg=1,nfg
c             write(6,*) 'wwwsolu',(fmoelm(j,ifg),j=1,3)
c             write(6,*) 'wwwsolv',(emocdr((ifg-1)*10+j),j=8,10)
              dsolute=sqrt(ddot(3,fmoelm(1,ifg),1,fmoelm(1,ifg),1))
              dsolvent=sqrt(ddot(3,emocdr((ifg-1)*10+8),1,
     *                      emocdr((ifg-1)*10+8),1))
           call vadd(fmoelm(1,ifg),1,emocdr((ifg-1)*10+8),1,dsoltot,1,3)
              dsolabs=sqrt(ddot(3,dsoltot,1,dsoltot,1))
              epsdeff=1
              if(dsolabs.gt.1.0D-06) epsdeff=dsolute/dsolabs
              omegaeff=0
              if(dsolute*dsolvent.gt.1.0D-06) omegaeff=
     *           acos(ddot(3,fmoelm(1,ifg),1,emocdr((ifg-1)*10+8),1)/
     *           (dsolute*dsolvent))/acos(-1.0D+00)*180.0D+00
              write(iw,9367) ifg,esolv(ifg),esolv(ifg)*tokcal,epsdeff,
     *                       omegaeff
            enddo
          endif
c         write(6,*) 'wwwpcm-esi',ges 
          if(mdoutmin) then
            write(iw,9200) symbe,esstr,1,emontot+ges
            write(iw,9200) symbe,solstr,1,emontot+gsolv/tokcal
          endif
        endif
        if(nbody.ge.2) then
          if(masout) write(iw,9370)
          ijfg=0
          rij=-one
          do ifg=1,nfg
            do 230 jfg=1,ifg-1
              ijfg=ijfg+1
              if(molfrg(ifg).eq.0.and.molfrg(jfg).eq.0.and.nbody.eq.2) 
     *        then
c               For FMO3, process all dimers fully no matter what molfrg
                go to 230
              end if
              if(modfd.ne.0) then
                if(iactfg(ifg).eq.0.and.iactfg(jfg).eq.0) then
                  go to 230
                end if
                if(layfrg(ifg).lt.nlayer.or.layfrg(jfg).lt.nlayer) then
                  go to 230
                end if
              end if
              if(needr.ne.0) then
c               separated dimers have no esolv set (since monomer 
c               values were not global summed before edimer).
                rij=fmodist(ifg,0,0,jfg)
                if(resdim.ne.0.and.rij.gt.resdim) then
                  esolv(nfg*2+ijfg)=esolv(ifg)+esolv(jfg)
c                 note that these values will be used in trimers below (FMO3)
                  if(ignoresd) goto 230
                endif
              endif
              eesiij=(esolv(nfg*2+ijfg)-esolv(ifg)-esolv(jfg))*gesfact
c             write(6,*) 'wwwesd',ifg,jfg,esolv(ifg)*tokcal,
c    *                        esolv(jfg)*tokcal,esolv(nfg*2+ijfg)*tokcal
              ges=ges+eesiij
              eesiij2=esolv(nfg*2+nfg2+ijfg)*tokcal
              eesiij3=eesiij*tokcal
              edispij=emocdr(10*nfg+(ijfg-1)*2+1)
              erepij=emocdr(10*nfg+(ijfg-1)*2+2)
              esolvij=eesiij2+eesiij3+edispij+erepij
              surfesij=emocdr((ifg-1)*10+6)+emocdr((jfg-1)*10+6)
              coverij=surfesij/stotes*1.0d+02
c             coverij=emocdr((ifg-1)*10+6)*emocdr((jfg-1)*10+6)/stotes/
c    *                stotes*1.0d+02
c             if(masout) write(iw,9375) ifg,jfg,esolv(nfg*2+ijfg)*tokcal,
              if(masout) write(iw,9375) ifg,jfg,surfesij,coverij,eesiij2
     *                                 ,eesiij3,edispij,erepij,esolvij
c             if(masout) write(iw,*) 'wwwes2+3',ifg,jfg,eesiij2+eesiij3
              gsolv=gsolv+esolvij
              esolv(nfg*2+nfg2+ijfg)=esolvij
c             We overwrite the ES2 set of values and the units are kcal/mol! 
 230        continue
          enddo
          if(iefmorun.ne.0) then
            gsolv=efmopcmg*tokcal
            ges=efmopcmg
          endif
          if(pcmprp.and.mdoutmin) write(iw,9355) 2,gsolv,2,ges*tokcal
        endif
      endif
C
c     if(nbody.lt.2.and.autosa) goto 1110
      if(nbody.lt.2.or..not.convSCC) goto 1000 
c
      looppieda3=0
c     Option to avoid computing effective PIEs in FMO3.
      if(iand(modprp,524288).ne.0.or.savemem.or.needdm.or.
     *  iand(nprfmo,3).gt.1) looppieda3=-1
c
 1100 continue
      autosa1=autosa.and.(nbody.lt.3.or.looppieda3.eq.1)
      if(masout.and..not.savemem2) then
        if(nbody.eq.3.and.looppieda3.eq.1) then
          write(iw,8131)
        else
          write(iw,9121)
          if(ipieda.ne.0) write(iw,8130)
        endif
        if(nbsse.eq.0) then
          if(ipieda.ne.0) then
            if(dopleda) then
              write(iw,8125) aapl0
            else
              if(dftbfl) then
                write(iw,9124)
              else
                distr='Erc+di    Gsol'
                if(dodc) distr=' Edisp   Gsol'
                if(dodft) distr=' Erc      Gsol'
                if(dispdft.and.dodft) distr=' Erc      Edisp   Gsol'
                if(do3c) distr=' Ebs      Edisp     Gsol'
                write(iw,9125) distr
              endif
            endif
          else 
            if(dolat.or.dopbcmd) then
              if(.not.(gcorrel.or.dodc)) then
                write(iw,8120)
              else
                write(iw,8121)
              endif
            else
              if(.not.(gcorrel.or.dodc)) then
                write(iw,9122)
              else
                write(iw,9123)
              endif
            endif
          endif 
        else if(nbsse.eq.1.or.nbsse.eq.2) then
          write(iw,9126)
        else if(nbsse.eq.3) then
          if(nend.gt.6) write(iw,9128)
          if(nend.le.6) write(iw,9129)
        endif
      endif
      ijfg=0
      e0mon=0
      e0monu=0
      edimtot=zero
      edimtotu=zero
      edimcon=zero
      edimconu=zero
      edimijq=zero
      edimijqu=zero
      eddim=zero
      eddimu=zero
      debsse=zero
      debsseu=zero
      debsseij=zero
      debsseuij=zero
      edimbs=zero
      ees=zero
      eex=zero
      ect=zero
      edi=zero
      erc=zero
      eso=zero
      eesbda=zero
      eexbda=zero
      ectbda=zero
      edibda=zero
      esobda=zero
      ebdaes=zero
      ebdaex=zero
      ebdact=zero
      ebdadi=zero
      ebdaso=zero
      rij=-one
      sepmax=zero 
      eesiij=zero
      esolvij=zero
c     esolvesij=zero
c     epsij=zero
c     ires=0
      if(subsys) then
        call vclr(emolfrg,1,nfg*nsubsys*2)
        call vclr(subprp,1,nsubsys*3+nsubsys*nsubsys)
      else if(modmol.ne.0) then
        call vclr(emolfrg,1,nfg*3)
      endif
      do ifg=1,nfg
        ilay=layfrg(ifg)
        njfg=ifg-1
        do 500 iu=0,nunint
        icurunt=iu
        if(iu.gt.0) then
          if(IPBCFST.eq.1)    goto 500
          if(nsymeq(iu).le.0) goto 500
        endif
        do jfg=1,njfg
c       do jfg=1,ifg-1
          ijfg=ijfg+1
          jlay=layfrg(jfg)
          ijlay=min(ilay,jlay)
c         if(molfrg(1).ge.0) ires=molfrg(ifg)+molfrg(jfg)
          ires=molfrg(ifg)+molfrg(jfg)
          if(exclmol) then
            skipc=ires.ge.0.and.ires.ne.2
          else
            skipc=ires.eq.0.or.ires.eq.2
          endif
          if(iand(modmol,8).ne.0.or.nbody.eq.3) skipc=.false.
          if(iand(modmol,1).ne.0.and.skipc) goto 100 
          if(modfd.ne.0) then
            if(iand(modfd,2).ne.0) then
              if(iactfg(ifg)+iactfg(jfg).eq.0) goto 100
            else
              if(layfrg(ifg)+layfrg(jfg).le.3.and.
     *           iactfg(ifg)+iactfg(jfg).eq.0) goto 100
            endif
          endif
          if(dofed.and.(iactfg(ifg).eq.0.or.iactfg(jfg).eq.0)) goto 100
          if(needr.ne.0) rij=fmodist(ifg,0,0,jfg) 
c         write(6,*) 'www',needr,rij,ifg,jfg  
c         Very radical: exit loop, because there is no monomer E to subtract
c         (nor anything to print).
          if(resdim.gt.0.and.rij.gt.resdim.and.savemem2) goto 100
          if(resdim.gt.0.and.rij.gt.resdim.and.ignoresd) goto 100
          sepmax=max(sepmax,rij)
          irij=int(rij*1.0D+02)
          ijcharge=ichfg(ifg)*ichfg(jfg)
          memon=1
c         memon1=0
          memon1=2
c         if(gcorrel) memon1=1
          if(iand(lcorrel(ijlay),1).ne.0) memon1=1
          if(modfd.ne.0.and.layfrg(ifg)+layfrg(jfg).eq.3) then
           if(gcorrel.or.dodc) memon=2
           memon1=3
          else
          if(resdim.ne.0) then
            if(rij.gt.resdim) memon1=3
            if(gcorrel) then
c             ires=molfrg(ifg)+molfrg(jfg)
              if((rij.gt.rcorsd.and.rcorsd.ne.0).or.skipc) then
                memon=2
                memon1=max(memon1,2)
              endif
            endif 
            if(rij.gt.resdim.and.dodc.and..not.dodcesd) memon=2
          endif 
          endif
          if((scffrg(ifg).eq.rmc.or.scffrg(jfg).eq.rmc).and.
     *       fmoscf(ijlay).eq.rmc.and.(resdim.eq.0.or.rij.le.resdim)) 
     *      memon1=4 
          if(doexc.and.memon1.eq.1.and.
     *      (ifg.ne.iexcit(1).and.jfg.ne.iexcit(1))) memon1=2
          memonb=memon+2
          if(dopleda) memon=memonb
c         corri=corrs(memon)
          corri=corrs(memon1) 
          if(modfmm.ne.0.and.memon1.eq.3) then
            call mmdist(ifg,0,0,jfg,t,radius,ty2z,ratio,mmdim)
            if(mmdim.ne.0) then
              corri='s'
              if(iand(modfmm,2).ne.0) goto 100
            endif 
          endif 
          mixdim=modfd.ne.0.and.layfrg(ifg)+layfrg(jfg).eq.3
          emonij=emon(ifg,memon,ijlay)+emon(jfg,memon,ijlay)
          if(mixdim) emonij=emon(ifg,memon,ilay)+emon(jfg,memon,jlay)
          emonijsav=emonij
          if(nbody.eq.2.and.(resdim.gt.0.and.rij.gt.resdim .or.mixdim))
     *    then
            emonij=0
          else
            emonijsav=0
          endif
          if(savemem2) then
c           edimij1=-emonij
            edimij1=0
          else
            edimij1=edim(ijfg,1)-emonij
          endif
          if(reducee.or.savemem2) then
            edimij2u=zero
          else
            edimij2u=edim(ijfg,2)
          endif
c         edimij2=edim(ijfg,4)
c         Exclude ES dimers.
          scfnone=scffrg(ifg).eq.rnone.or.scffrg(jfg).eq.rnone
          if(scfnone) then
            edimij1=zero
          else
            if(.not.savemem2.and.memon1.ne.3) then
              if(edim(ijfg,1).eq.zero) ndiv2=ndiv2+1
            endif
          endif
          if(mixdim) then
            emixdim=emixdim+edimij1
            if(rij.lt.rmixdim) then
              rmixdim=rij
              imixdim=ifg
              jmixdim=jfg
            endif
          endif
c
c         This is the latest fix to replace delta-D (MP2) by RHF.
c         to put back, uncontract the if below and remove the line next to it.
c
c         if(.not.doddcor.or.memon.eq.2) edimij2=edimij2u
          if(fmoq) edimijq=edimq(ijfg)
          edimijqu=edimijq
          edimij2=edimij2u
          eddim=eddim+edimij2
          eddimu=eddimu+edimij2u
c         if(nend.lt.3) then
          emonijusav=emonijsav
          if(.not.(gcorrel.or.dodc)) then
            edimij1u=edimij1
          else
            memon2=2 
            if(dopleda) memon2=memon2+2
            emoniju=emon(ifg,memon2,ijlay)+emon(jfg,memon2,ijlay)
          if(mixdim) emoniju=emon(ifg,memon2,ilay)+emon(jfg,memon2,jlay)
            emonijusav=emoniju
            if(nbody.eq.2.and.(resdim.gt.0.and.rij.gt.resdim.or.mixdim))
     *      then
              emoniju=0
            else
              emonijusav=0
            endif
            if(savemem2) then
c             edimij1u=-emoniju
              edimij1u=0
            else
              if (reducee) then
                edimij1u=edim(ijfg,1)-emoniju
              else
                edimij1u=edim(ijfg,3)-emoniju
              end if
            endif
          endif
          if(scfnone) then
            edimij1u=zero
          endif
          if(doexc.and.iexcit(2).eq.1) edimij1=edimij1u
c         No excited dimers: use ground state E'.
          if(abs(rij).lt.1.0D-08) then
            edimcon=edimcon+edimij1+edimij2+edimijq
            edimconu=edimconu+edimij1u+edimij2u+edimijqu
          endif
c         if(pcmprp) then
          if(nfmopcm.ne.0) then
c           eesiij=esolv(nfg*2+ijfg)-esolv(ifg)-esolv(jfg)
            esolvij=esolv(nfg*2+nfg2+ijfg)
c           esolvesij=esolvij
c    *            -emocdr(10*nfg+(ijfg-1)*2+1)-emocdr(10*nfg+(ijfg-1)*2+2)
          endif
c         if(nend.gt.4) then 
          if(nbsse.ne.0) then
            edimij3=edim(ijfg,5)
            edimij4=edim(ijfg,6)
            if(edimij3.eq.sawdust) edimij3=zero
            if(edimij4.eq.sawdust) edimij4=zero
            if(nend.le.6) then
              edimijcb=edim(ijfg,6)
            else
              edimijcb=edim(ijfg,7)
            endif
          endif
          if(.not.savemem2) then
            edimes=edim(ijfg,nedimes)
            edimex=edim(ijfg,nedimex)
          endif
c         for domdan, quit before printing as soon as possible...
          if(masout.or.domdan) then
          if(nbsse.eq.0) then 
            detot=edimij1+edimij2+edimijq
            detotu=edimij1u+edimij2u+edimijqu
            detot1=detot+esolvij/tokcal
            detot2=detot1
            detot1u=detotu+esolvij/tokcal
c           Remove artificial tens of hartrees for connected monomers.
c           If distances are not available, use plain energy to judge.
            if(ipieda.eq.0.and.(needr.ne.0.and.rij.eq.zero.or.
     *         needr.eq.0.and.edimij1.lt.-1.0D+01)) then
              detot=edimij2+edimijq
              detotu=edimij2u+edimijqu
            endif
            if(needr.ne.0.and.rij.eq.zero.or.
     *         needr.eq.0.and.edimij1.lt.-1.0D+01) then
              detot2=0
c             Set the unconnected contribution to 0 for connected dimers.
c             Here, -10 is chosen to cover C-C bonds (<-14).
c             detot1 any PIE
c             detot2 connected PIE
            endif
c           for EFMO, remove the dimer polarization from dimer printout
            if(iefmorun.ne.0) then
              detot=detot-edimij2
              detotu=detotu-edimij2u
            endif
            ctij=zero
c           ctsp=zero
            if(primul2) ctij=ctdim(ijfg)
c           if(urospn.and.imect.eq.5) ctsp=ctspin(ijfg)
c           rij means two fragments are connected by a bond and thus
c           edimij1 contains some meaningless (well,?) value.
c           if(nfmopcm.ne.0) epsij=-detot/esolvesij
            if(domdan) then
              x(ledimmd+ijfg-1)=detot1
              if(.not.masout.and.ipieda.eq.0) goto 400
c             skip the printing section after saving data.
            endif
            if(autosa.and.isasign.eq.1) 
     *        sadim(ijfg)=sadim(ijfg)+iapsign*detot1
c
            if(ipieda.ne.0) then
c
              edimdi=edimij1-edimij1u
              if(.not.(gcorrel.or.dodc)) edimdi=zero
              edimct=detot-edimes-edimex-edimdi
              if(dodft) then
                ecdft=ecorrdft(ijfg+nfg)-ecorrdft(ifg)-ecorrdft(jfg)
                if(resdim.ne.0.and.rij.gt.resdim) ecdft=0
                edimct=edimct-ecdft
                edimdi=edimdi+ecdft
              else
                ecdft=0
              endif
              if(do3c) then
                edimbs=ecorrdft(ijfg+nfg)-ecorrdft(ifg)-ecorrdft(jfg)
                if(resdim.ne.0.and.rij.gt.resdim) edimbs=0
              endif
              if(resdim.ne.0.and.rij.gt.resdim.and.nbody.gt.2.and.
     *          (dodft.or.do3c)) 
     *          ecorrdft(ijfg+nfg)=ecorrdft(ifg)+ecorrdft(jfg)
c               For PIEDA/3, save ES dimer value to make life easier.
              if(dftbfl) then
c               Calculate Delta-E'ES
                e1es=ecorrdft(ijfg+nfg*2)-ecorrdft(ifg)-ecorrdft(jfg)
                if(dopleda) e1es=ecorrdft(ijfg+nfg*2)
     *                          -efmo0(ifg,2,ilay)-efmo0(jfg,2,ilay)
c               Calculate DeltaE'(CT*ES) = DeltaE'ES - E_ES 
                e1ctes=e1es-edimes
                if(resdim.ne.0.and.rij.gt.resdim) e1ctes=0
c               Calculate total DeltaE(CT*ES) = DeltaE'(CT*ES) + E_V
c               and store it as edimct
                edimct=e1ctes+edimij2
c               write(6,*) 'wwwct',e1ctes*tokcal,edimij2*tokcal
c               Calculate DeltaE0 as the remainder 
c               (deltaE' - DeltaE'(CT*ES) - E_ES + DeltaErep+DeltaEP) 
c               and store it as edimex.
                edimex=detot-edimes-edimct-edimdi 
c
                if(resdim.ne.0.and.rij.gt.resdim.and.nbody.gt.2)
     *           ecorrdft(ijfg+nfg*2)=ecorrdft(ifg)+ecorrdft(jfg)+edimes
              endif
              if(dopleda.and.memon1.ne.3) then 
                epl0ij=epl0ds(ifg,1)+epl0ds(ifg,2)+
     *                 epl0ds(jfg,1)+epl0ds(jfg,2)
                if(dodft) epl0ij=epl0ij+epl0ds(ifg,3)+epl0ds(jfg,3)
                edimct=edimct-epl0ij
c               In PL0 runs polarisation is double counted.
c               For DFT, the PL*C term is double counted because PIE
c               is computed as PL0-0, which includes PL*C.
              else
                epl0ij=zero
              endif
              edimij1b=edimij1
              eesolvij=esolvij/tokcal
              if(n0bda.ne.0) then
                isign=1
                if(ndualb.ne.0.and.ibloop.eq.1) isign=-1
                e0bda(1,i0bda,1)=e0bda(1,i0bda,1)+isign*edimes
                e0bda(2,i0bda,1)=e0bda(2,i0bda,1)+isign*edimex
                e0bda(3,i0bda,1)=e0bda(3,i0bda,1)+isign*edimct
                e0bda(4,i0bda,1)=e0bda(4,i0bda,1)+isign*edimdi
                e0bda(5,i0bda,1)=e0bda(5,i0bda,1)+isign*eesolvij
              else if(nbdfg.ne.0) then
                ebdaes=ebdaes+edimes
                ebdaex=ebdaex+edimex
                ebdact=ebdact+edimct
                ebdadi=ebdadi+edimdi
                ebdaso=ebdaso+eesolvij
c               BDA corrections will not work with D in DFT/D.
                call bdasub(ifg,jfg,gcorrel.or.dodc.or.dodft,edimij1b,
     *                      detot,detotu,edimes,edimex,edimct,edimdi,
     *                      eesolvij,e0mon,e0monu,indat,iabdfg,jabdfg,
     *                      e0bda(1,1,ijlay),ires0)
                ebdaes=ebdaes-edimes
                ebdaex=ebdaex-edimex
                ebdact=ebdact-edimct
                ebdadi=ebdadi-edimdi
                ebdaso=ebdaso-eesolvij
                if(ires0.eq.0) then
                  eesbda=eesbda+edimes
                  eexbda=eexbda+edimex
                  ectbda=ectbda+edimct
                  edibda=edibda+edimdi
                  esobda=esobda+eesolvij
                endif
                esolvij=eesolvij*tokcal
              endif
              ees=ees+edimes
              eex=eex+edimex
              ect=ect+edimct
              edi=edi+edimdi
              if(domdan) then
c               The most common 5 components.
                x(ledimmd+nfg2*6+ijfg-1)=edimes
                x(ledimmd+nfg2*7+ijfg-1)=edimex
                x(ledimmd+nfg2*8+ijfg-1)=edimct
                x(ledimmd+nfg2*9+ijfg-1)=edimdi
                x(ledimmd+nfg2*10+ijfg-1)=esolvij/tokcal
c               Some PIEDAs have 6 components, and 7 can be only in AP.
                x(ledimmd+nfg2*11+ijfg-1)=ecdft
                x(ledimmd+nfg2*12+ijfg-1)=edimbs
              endif
              if((prtdst(4).eq.0.or.abs(detot)*tokcal.gt.prtdst(4)).and.
     *            .not.savemem2.and..not.scfnone) then
              if(dopleda) then
                write(iw,9132) ifg,jfg,corri,ijlay,ijcharge,rij,ctij,
     *                         edimij1b*tokcal,edimij2*tokcal,
     *                         detot*tokcal,epl0ij*tokcal,
     *                         edimes*tokcal,edimex*tokcal,edimct*tokcal
     *                        ,edimdi*tokcal
              else
                eso=eso+esolvij
c               if(dodft.and.dodc .or. do3c) then
                if(dodft.and.dispdft .or. do3c) then
c               PIEDA has 6 components:
c               DFT/D: ES, EX, CT+MIX, RC, DI, SOLV
c               HF-3c: ES, EX, CT+MIX, BS, DI, SOLV
                if(do3c) then
                  edimextra=edimbs
                else
                  edimextra=ecdft
                endif
                edi=edi-edimextra
                erc=erc+edimextra
                write(iw,9134) ifg,jfg,corri,ijlay,ijcharge,rij,ctij,
     *                         edimij1b*tokcal,edimij2*tokcal,
     *                         detot*tokcal+esolvij,edimes*tokcal,
     *                         edimex*tokcal,edimct*tokcal,
     *                       edimextra*tokcal,(edimdi-edimextra)*tokcal,
     *                         esolvij
                else
c               PIEDA has 5 components ES, EX, CT+MIX, DI+RC, SOLV
                write(iw,9134) ifg,jfg,corri,ijlay,ijcharge,rij,ctij,
     *                         edimij1b*tokcal,edimij2*tokcal,
     *                         detot*tokcal+esolvij,edimes*tokcal,
     *                         edimex*tokcal,edimct*tokcal,
     *                         edimdi*tokcal,esolvij
                endif
c               write(6,*) 'epsIJ',ifg,jfg,detot*tokcal/
c    *                              (detot*tokcal+esolvij),edimes*tokcal
              endif 
              endif 
            else  
              if((prtdst(4).eq.0.or.abs(detot)*tokcal.gt.prtdst(4)).and.
     *           .not.savemem2.and..not.scfnone) then
              edimij1p=edim(ijfg,1)+emonijsav
              edimij3p=edim(ijfg,3)+emonijusav
c             write(6,*) 'wwweee',emonijsav,emonijusav,edimij1p
                if(dolat.or.dopbcmd) then
c                 nsymeq(iu=0) is well defined?
                  if(iu.eq.0) then
                    nsymeqiu=0
                  else
                    nsymeqiu=nsymeq(iu)
                  endif
                  if(.not.(gcorrel.or.dodc)) then
                    write(iw,9161) ifg,jfg,iu,nsymeqiu,corri,ijlay,
     *                          ijcharge,rij,ctij,edimij1p,edimij1,
     *                             edimij2,detot*tokcal
                  else
                    write(iw,9163) ifg,jfg,iu,nsymeqiu,corri,ijlay,
     *                             ijcharge,rij,ctij,edimij1p,
     *                             edimij3p,edimij1,edimij1u,
     *                             edimij2,detot*tokcal
                  endif
                else
                  if(.not.(gcorrel.or.dodc)) then
                   write(iw,9131) ifg,jfg,corri,ijlay,ijcharge,rij,ctij,
     *                            edimij1p,edimij1,edimij2,
     *                            esolvij,detot*tokcal+esolvij
c    *                            edim(ijfg,1),edimij1,edimij2,edimijq,
                  else
                   edimij3p0=edimij3p
                   if (reducee) edimij3p0=zero
                   write(iw,9133) ifg,jfg,corri,ijlay,ijcharge,rij,ctij,
     *                            edimij1p,edimij3p0,edimij1,
     *                            edimij1u,edimij2,esolvij,
     *                            detot*tokcal+esolvij
c    *                            edimij1u,edimij2,edimijq,detot*tokcal
                  endif 
                endif 
              endif 
            endif 
            debsseij=zero
            if(ndualb.ne.0) edimmixd(ijfg,ibloop+4)=detot1
          else if(nbsse.eq.1.or.nbsse.eq.2) then 
            debsseij=edimij3+edimij4-edimij1-edimij2
            if(edimij3+edimij4.eq.zero) debsseij=zero
            detot=edimij1+edimij2+debsseij
            write(iw,9140) ifg,jfg,corri,ijlay,edim(ijfg,1),edimij1,
     *            edimij1u,edimij2,edimij3,edimij4,debsseij,detot*tokcal
          else if(nbsse.eq.3) then 
            debsseij=emon(ifg,memonb,ijlay)+emon(jfg,memonb,ijlay)-
     *               edimij3
            debsseuij=emon(ifg,4,ijlay)+emon(jfg,4,ijlay)-edimijcb
c           take care of separated dimers that ignore BSSE
            if(edimij3.eq.zero) debsseij=zero
            if(edimij3.eq.zero.or.nend.lt.7) debsseuij=zero
            detot=edimij1+edimij2+debsseij
            if(nend.ge.7) then 
              write(iw,9150) ifg,jfg,corri,ijlay,irij,edim(ijfg,1),
     *          edimij1,edimij1u,edimij2,debsseij,debsseuij,detot*tokcal
            else
              write(iw,9152) ifg,jfg,corri,ijlay,irij,edim(ijfg,1),
     *           edimij1,edimij2,debsseij,detot*tokcal
            endif
          endif
          endif
  400     continue
c         This continue/endif is for if(masout.and.domdan).
c
          if(dopleda.and.memon1.ne.3) then
c           Replace energies by PL0 values to get the PIEDA energy.
c           memon1=3 means separated dimers, in which case
c           edimij1 is measured against FMO0, and there is no need to
c           adjust it. 
            edimij1=edim(ijfg,1)-emon(ifg,memon-2,ijlay)
     *                          -emon(jfg,memon-2,ijlay)
            edimij2=edim(ijfg,4)
            edimij1u=edimij1
            edimij2u=edimij2
            if(gcorrel.or.dodc)
     *        edimij1u=edim(ijfg,3)-emon(ifg,2,ijlay)-emon(jfg,2,ijlay)
            if(resdim.gt.0.and.rij.gt.resdim) call abrtx("Error in prp")
c           separated dimers should not get here: the emon code was not adopted.
          endif
          edimtot=edimtot+edimij1+edimij2+edimijq
          edimtotu=edimtotu+edimij1u+edimij2u+edimijqu
c         write(6,*) 'wwwEE',ifg,jfg,edimij1,edimij1u,edimij2
          if(nappri.ne.0) then 
            if(rij.le.rappri(1)) then
              iappri=1
            else
              iappri=int((rij-rappri(1))/rappri(3))+2
            endif
            if(iappri.le.nappri) then
              erapp(iappri,1)=erapp(iappri,1)+edimij1+edimij2+edimijq
              erapp(iappri,2)=erapp(iappri,2)+edimij1u+edimij2u+edimijqu
              nrapp(iappri,1)=nrapp(iappri,1)+1 
              if(memon1.eq.1) nrapp(iappri,2)=nrapp(iappri,2)+1
            endif
          endif
          debsse=debsse+debsseij
          debsseu=debsseu+debsseuij
c         if(ires.eq.0.or.ires.eq.2) then
          if(.not.exclmol.and..not.subsys.and.ires.eq.1.and.
     *       modmol.ne.0) then
c           moldmol.ne.0 check is a double safety.
            if(molfrg(ifg).gt.0) then
              emolfrg(ifg,1)=emolfrg(ifg,1)+detot1
              emolfrg(ifg,2)=emolfrg(ifg,2)+detot1u
              if(ipieda.ne.0) emolfrg(ifg,3)=emolfrg(ifg,3)+edimes
            endif
            if(molfrg(jfg).gt.0) then
              emolfrg(jfg,1)=emolfrg(jfg,1)+detot1
              emolfrg(jfg,2)=emolfrg(jfg,2)+detot1u
              if(ipieda.ne.0) emolfrg(jfg,3)=emolfrg(jfg,3)+edimes
            endif
          endif
          if(subsys) then
            isub=molfrg(ifg)
            jsub=molfrg(jfg)
            isub1=isub+nsubsys
            jsub1=jsub+nsubsys
            emolfrg(ifg,jsub)=emolfrg(ifg,jsub)+detot1
            emolfrg(ifg,jsub1)=emolfrg(ifg,jsub1)+detot2
            emolfrg(jfg,isub)=emolfrg(jfg,isub)+detot1
            emolfrg(jfg,isub1)=emolfrg(jfg,isub1)+detot2
c           write(6,*) 'wwwij',ifg,jfg,rij,detot1,detot2
c           write(6,*) 'wwwii',ifg,emolfrg(ifg,jsub),emolfrg(ifg,jsub1)
c           write(6,*) 'wwwjj',jfg,emolfrg(jfg,isub),emolfrg(jfg,isub1)
c           ctij is from I to J, meaning J gets ctij added.
            subprp(isub)=subprp(isub)-ctij
            subprp(jsub)=subprp(jsub)+ctij
          endif
 100    continue
        enddo
        njfg=nfg
  500 continue
      enddo
c     if(nbody.eq.3.and.ipieda.gt.0.and.looppieda3.eq.0) goto 1000
c
c     Add the lump contributions from all dimers.
c
      if(savemem2) then
        edimtot=edimtot+edimsum(1)+edimsum(2)
        if(gcorrel.or.dodc) then
          edimtotu=edimtotu+edimsum(3)+edimsum(2)
        else
          edimtotu=edimtot
        endif
      endif
c
      if(molfrg(1).ge.0.and..not.exclmol.and..not.subsys) then
c       presumably, only modmol=1 can get here.
        write(iw,9160) 
        do ifg=1,nfg
          if(molfrg(ifg).gt.0) then
            if(ipieda.ne.0) then
              write(iw,9175) ifg,frgnam(ifg),(emolfrg(ifg,i)*tokcal,
     *                       i=3,1,-1)
            else
              write(iw,9170) ifg,frgnam(ifg),(emolfrg(ifg,i)*tokcal,
     *                       i=2,1,-1)
            endif
          endif
        enddo
c       No total 2-body properties are available.
c       goto 1000 
      endif
c
c1110 continue
      if(subsys) then
        write(iw,9165)
        gsolvi=0
        do ifg=1,nfg
          if(nfmopcm.gt.0) 
     *      gsolvi=ESOLV(IFG)-ESOLV(IFG+nfg)+(emocdr((ifg-1)*10+1)
     *            +emocdr((ifg-1)*10+2)+emocdr((ifg-1)*10+3))/tokcal
          ilay=layfrg(ifg)
          isub=molfrg(ifg)
          isub1=isub+nsubsys
          eparti=emon(ifg,1,ilay)+gsolvi+emolfrg(ifg,isub)/2
          ebbi=eparti-emolfrg(ifg,isub1)/2
c         write(6,*) 'wwwi',ifg,emolfrg(ifg,isub),emolfrg(ifg,isub1)
          ifgsa=ifg+nfgasub 
          if(autosa1) 
     *      samon(1,ifgsa)=samon(1,ifgsa)+isasign*iapsign*eparti
c         subtract unconnected.
c         interactions within subsystem are double counted and must be halved.
c         halve self-interactions to avoid double counting
          write(iw,9166) ifg,frgnam(ifg),molfrg(ifg),eparti,ebbi,
     *                   (emolfrg(ifg,i+nsubsys)*tokcal,i=1,nsubsys)
c    *                   gsolvi*tokcal,
          subprp(isub+nsubsys)= subprp(isub+nsubsys)+eparti
c         subprp(isub+nsubsys*2)=subprp(isub+nsubsys*2)+ebbi
          loop=nsubsys*3+(isub-1)*nsubsys
c         loop=nsubsys*3+((isub-1)*(isub-2))/2
c         do j=1,isub-1
          do j=1,nsubsys
            if(j.ne.isub) then
c             compute EU,int(i,j)
              subprp(loop+j)=subprp(loop+j)+emolfrg(ifg,j+nsubsys)
c             Compute sum(j!=i) EC,int(i,j)
              subprp(isub+nsubsys*2)=subprp(isub+nsubsys*2)
     *                            +emolfrg(ifg,j)-emolfrg(ifg,j+nsubsys)
            endif
c           to which subprp(isub+nsubsys) must be added to get BB.
          enddo
c         At this point subprp(isub) has charge due to CT only.
c         Now add original charge.
          subprp(isub)=subprp(isub)+ichfg(ifg)
        enddo
        write(iw,9167)
        do i=1,nsubsys
c         loop=nsubsys*3+((i-1)*(i-2))/2
          loop=nsubsys*3+(i-1)*nsubsys
c         The matrix is symmetric but for data processing convenience
c         do square matrix.
          einti=subprp(i+nsubsys)
          ebbi=einti+subprp(i+nsubsys*2)/2
c         ebbi=einti+subprp(i+nsubsys*2)/2-subprp(loop+i)/2
c         ebbi=einti+subprp(i+nsubsys*2)/2
c         Divide intersubsystem int,C evenly.
c         zero out self-interactions (these are included in Epart).
c         subprp(loop+i)=0 
          eparti=ebbi
          do j=1,nsubsys
            eparti=eparti+subprp(loop+j)/2
          enddo
          write(iw,9168) i,subprp(i),einti,eparti,ebbi,
     *                   (subprp(loop+j)*tokcal,j=1,nsubsys)
c    *               subprp(i+nsubsys*2)*tokcal,(subprp(loop+j),j=1,i-1)
        enddo
      endif
c
c     if(nbsse.eq.0) then
c       etot=edimtot-(nfg-2)*emontot
c     else
      if(iand(modfmm,2).ne.0) then
        esdtot=ddot(nfg,esdi,1,one,0) 
        if(maswrk) then
          if(iand(nprfmo,3).lt.3) then
            write(iw,9720)
            do ifg=1,nfg
              write(iw,9730) ifg,esdi(ifg)*tokcal 
            enddo
          endif
          write(iw,9740) esdtot*tokcal,esdtot
        endif
        edimtot=edimtot+esdtot
        edimtotu=edimtotu+esdtot
c       edimtotu=edimtot+esdtot
      endif
      etot=emontot+edimtot
      etotu=emontotu+edimtotu
c     At this point fmoelm(1,nfg+3) contains the 1st order moments
c                   fmoelm(1,nfg+2) contains the 3rd order correction 
c                   fmoelm(1,nfg+1) contains the 2nd order correction 
      if(nfmoelm.ne.0) 
     *  call daxpy(nfmoelm,one,fmoelm(1,nfg+3),1,fmoelm(1,nfg+1),1)

      if(iefmorun.ne.0) then
c       for efmo we need to add the total polarization
c       energy on both the energy and the gradient
        etot = etot + efmoetot
        etotu = etotu + efmoetot
        IEFMORT = 4
      endif

c     endif
c     enucr=ENUC(natfmo,fmoZan,fmoC)
c
c     for BSSE runs, the total energy (as will be used during optimisation)
c     includes BSSE correction, but the gradient doesn't!
c 
      if(ndiv2.ne.0.and.runtyp.ne.check) then
        write(iw,9500) ndiv2,2
        etot=0
        etotu=0
      endif
      if(skipesd) write(iw,9510)
c
      if(ipieda.ne.0.and.nbdfg.ne.0.and.n0bda.eq.0.and.e0mon.ne.0) then
        write(iw,9201) corstr(1),1,emontot+e0mon
        if(nenm.ge.2.and.(gcorrel.or.dodc)) then
          write(iw,9201) 'uncorr',1,emontotu+e0monu
          write(iw,9201) 'delta ',1,emontot-emontotu+(e0mon-e0monu)
        endif
      endif
      if(nbsse.eq.0) then
        write(iw,9200) symbe,corstr(1),2,etot
        if(ndualb.ne.0.and.looppieda3.eq.0) 
     *    epcmap(3)=epcmap(3)+iapsign*(etot+gsolv/tokcal)
        if(nfmopcm.ne.0.and.mdoutmin) then
          write(iw,9200) symbe,esstr,2,etot+ges
          write(iw,9200) symbe,solstr,2,etot+gsolv/tokcal
        endif
        if( iefmorun .gt. 0 ) write(iw,9200) symbe,'totpol',2,efmoetot
c       if(doddcor) write(iw,9200) 'corr-D',2,etot-eddim+eddimu
c       corr-D has no value because to get it one needs E'I with RHF delta-D. 
        if(nenm.gt.1) then
          if (reducee) etotu=zero
          if(.not.reducee) write(iw,9200) symbe,corstr(2),2,etotu
          if (reducee) etotu=etot
          if(.not.reducee) write(iw,9200) symbe,corstr(3),2,etot-etotu
          if(didcc) call fmoccp(2,etotu,extracc(1,2))
          if(didmp2) call fmompp(2,etotu,extracc(1,2))
          if(do3c) call fmoh3c(2,extracc(1,2))
        endif
      else
        write(iw,9210) corstr(1),etot,etot+debsse
        etot=etot+debsse
        if(nend.ge.7) then
          write(iw,9210) 'uncorr',etotu,etotu+debsseu 
          escf=escf+debsseu
        endif
      endif
      if(douhf.and.mdoutmin) 
     *  write(iw,9780) 2,sz0,2,szfmo(2)*(szfmo(2)+1),s2_0
c     if(nbody.eq.3.and.ipieda.gt.0.and.looppieda3.eq.0) goto 1000
c     EBB may be defined without distances known
c     but then the code should get the connecting
c     information from elsewhere.
      if (reducee.and.dodc) then
        emontotu = zero
        edimconu = zero
      end if
      if(nbdfg.ne.0.and.needr.ne.0.and.mdoutmin) then
cnb     for savemem2, the values of EBB are wrong - add ES?
        write(iw,9202) symbe,corstr(1),2,emontot+edimcon
        if(gcorrel.or.dodc) write(iw,9202) symbe,'uncorr',2,
     *                      emontotu+edimconu
      endif
      if(nappri.ne.0) then
        do i=2,nappri
          erapp(i,1)=erapp(i,1)+erapp(i-1,1)
          erapp(i,2)=erapp(i,2)+erapp(i-1,2)
          nrapp(i,1)=nrapp(i,1)+nrapp(i-1,1)
          nrapp(i,2)=nrapp(i,2)+nrapp(i-1,2)
        enddo
        write(iw,*) ' ' 
        do i=1,nappri
          riapp=rappri(1)+rappri(3)*(i-1)
c         sepmax=min(sepmax,resdim)
          if(riapp.gt.sepmax.and.riapp-rappri(3).le.sepmax) riapp=sepmax
          if(riapp.le.sepmax.and.sepmax.ne.0) then 
          if(gcorrel.or.dodc) then
           write(iw,9205) riapp,nrapp(i,2),nrapp(i,1),2,
     *       erapp(i,1)+emontot-erapp(i,2)-emontotu,erapp(i,2)+emontotu
          else
           write(iw,9206) riapp,nrapp(i,1),2,erapp(i,1)+emontot
          endif
          endif
        enddo
        write(iw,*) ' ' 
      endif
C
      if(reducee.and.dodc.and.masout) write (iw,9800)
c
      if(ipieda.ne.0) then
c       if(dodft.and..not.dodc) then
        if(dodft.and..not.dispdft) then
          erc=edi
          edi=0
        endif
        if(dopleda) then
c         if(dodft.and..not.dodc) eint0(4)=erc
          if(dodft.and..not.dodc) edi=erc
          eint0(1)=ees
          eint0(2)=eex
          eint0(3)=ect
          eint0(4)=edi
          epl0=epl0d+epl0s
          eib0=epl0+epl0di+ees+eex+ect+edi
          e1int=eesbda+eexbda+ectbda+edibda+esobda
          write(iw,9203) corstr(1),2,e0tot+e0mon+eib0
          if(gcorrel.or.dodc) then
            write(iw,9203) 'uncorr',2,e0totu+e0monu+eib0-epl0di-edi
          write(iw,9203) 'delta ',2,e0tot+e0mon-e0totu-e0monu+epl0di+edi
          endif
          if(ipieda.eq.2) then
            write(iw,9600) epl0d*tokcal,
     *                     epl0s*tokcal,
     *                         epl0*tokcal,
     *               bpl0,apl0,epl0di*tokcal,
     *                         ees*tokcal,
     *                         eex*tokcal,
     *                         ect*tokcal,
     *               bpl0,apl0,edi*tokcal,
     *                             eib0*tokcal,(eib0-epl0di-edi)*tokcal
            if(nbdfg.ne.0) then
             write(iw,9620) eesbda*tokcal,eexbda*tokcal,ectbda*tokcal,
     *                      bpl0,apl0,edibda*tokcal,esobda*tokcal,
     *                      e1int*tokcal,ebdaes*tokcal,ebdaex*tokcal,
     *                      ebdact*tokcal,bpl0,apl0,ebdadi*tokcal,
     *                      ebdaso*tokcal,e0mon*tokcal
            endif
          endif
        else
          ees0=eint0(1)
          eex0=eint0(2)
          ect0=eint0(3)
          edi0=eint0(4)
          epls_pls=ees-ees0-epl0s
          epl_pl=epld_pld+epls_pls
          epls_ex=eex-eex0
c         epl_plt=epl_pl
          epls_ct=ect+epl_pl-ect0
          epls_di=edi-edi0
          if(dodft) epls_di=erc-edi0
          epl0=epl0d+epl0s
          epl_di=epld_di+epls_di
c         if(dodft) epl_plt=epl_plt+epl_di
          eib0=epl0+epl0di+ees0+eex0+ect0+edi0
          eib=eib0+epls_ex+epls_ct+epl_di
c         eint=ees+eex+ect+edi
          eint=ees+eex+ect+edi+erc
          e1int=eesbda+eexbda+ectbda+edibda+esobda
          epl=epl0+epl_pl
          epldcorr=epl0d+epld_pld+epl0di+epld_di
          if(ipieda.eq.2) then
          write(iw,9203) corstr(1),2,e0tot+e0mon+eib0
          if(gcorrel.or.dodc) then
            write(iw,9203) 'uncorr',2,e0totu+e0monu+eib0-epl0di-edi0
         write(iw,9203) 'delta ',2,e0tot+e0mon-e0totu-e0monu+epl0di+edi0
          endif
            write(iw,9600) epl0d*tokcal,
     *                     epl0s*tokcal,
     *                         epl0*tokcal,
     *               bpl0,apl0,epl0di*tokcal,
     *                         ees0*tokcal,
     *                         eex0*tokcal,
     *                         ect0*tokcal,
     *               bpl0,apl0,edi0*tokcal,
     *                             eib0*tokcal,(eib0-epl0di-edi0)*tokcal
     *                            ,epls_ex*tokcal,
     *                             epls_ct*tokcal,
     *               bpl0,apl0,epld_di*tokcal,
     *               bpl0,apl0,epls_di*tokcal,
     *                   bpl0,apl0,epl_di*tokcal,
     *                                 eib*tokcal,
     *                                 (eib-epl0di-edi0-epl_di)*tokcal
c
            write(iw,9605) epld_pld*tokcal,epls_pls*tokcal,
     *                     epl_pl*tokcal,epl*tokcal,(ect0-epl_pl)*tokcal
     *                    ,epldcorr*tokcal
          endif
c
          if(dftbfl) then
            write(iw,9611) ees*tokcal,eex*tokcal,ect*tokcal,edi*tokcal,
     *                   eso,eint*tokcal+eso
            write(iw,9215) edftb0+eex
          else
              distr='Remainder correlation'
              dispstr='RC'
            if(do3c) then
              distr='Basis set correction '
              dispstr='BS'
            endif
            write(iw,9610) dispstr,ees*tokcal,eex*tokcal,ect*tokcal,
     *                     edi*tokcal,distr,dispstr,erc*tokcal,eso,
     *                     eint*tokcal+eso
          endif
          if(nbdfg.ne.0) then
             write(iw,9620) eesbda*tokcal,eexbda*tokcal,ectbda*tokcal,
     *                      bpl0,apl0,edibda*tokcal,esobda*tokcal,
     *                      e1int*tokcal,ebdaes*tokcal,ebdaex*tokcal,
     *                      ebdact*tokcal,bpl0,apl0,ebdadi*tokcal,
     *                      ebdaso*tokcal,e0mon*tokcal
             endif
        endif
      endif
c     if(nbody.eq.3.and.ipieda.gt.0.and.looppieda3.eq.0) goto 1000
      if(looppieda3.eq.1) goto 1000
c
c     write(iw,9229+i) (fmoelm(ielm,nfg+1),ielm=nelm(i-1)+1,nelm(i))
      if(iemom.ge.1) 
     *  write(iw,9230) 2,(fmoelm(ielm,nfg+1),ielm=nelm(0)+1,nelm(1)),
     *                  DNRM2(3,fmoelm(1,nfg+1),1)
      if(iemom.ge.2) 
     *  write(iw,9231) (fmoelm(ielm,nfg+1),ielm=nelm(1)+1,nelm(2))
      if(iemom.ge.3) 
     *  write(iw,9232) (fmoelm(ielm,nfg+1),ielm=nelm(2)+1,nelm(3))
      if(nfmoelm.ne.0) call dcopy(nfmoelm,fmoelm(1,nfg+1),1,DMX,1)
c
      if(primul2.and.masout) then
        write(iw,9420)
        nctsiz=min(nfg-1,4)
        ctmall=zero
        ctdall=zero
        do ifg=1,nfg
          if(modfd.eq.0.or.(layfrg(ifg).gt.1.and.
     *       (iand(modfd,2).eq.0.or.iactfg(ifg).ne.0))) then
          cttot=zero
          call vclr(ctmax,1,nctsiz)
          call viclr(ictmax,1,nctsiz)
          do jfg=1,nfg
            if(ifg.ne.jfg) then
              if(ifg.gt.jfg) then
                ctij=-ctdim(((ifg-1)*(ifg-2))/2+jfg)
                ctdall=ctdall+abs(ctij)
              else 
                ctij= ctdim(((jfg-1)*(jfg-2))/2+ifg)
              endif 
              cttot=cttot+ctij
              imin=1
              do ii=2,nctsiz
                if(abs(ctmax(ii)).lt.abs(ctmax(imin))) imin=ii
              enddo
c             Replace the smallest element in ctmax by ctij
              if(abs(ctmax(imin)).lt.abs(ctij)) then
                ctmax(imin)=ctij
                ictmax(imin)=jfg
              endif
            endif
          enddo
          call dasort(nctsiz,ctmax,ictmax)
          write(iw,9425) ifg,ichfg(ifg),cttot,
     *                   (ictmax(ii),ctmax(ii),ii=1,nctsiz)
          ctmall=ctmall+abs(cttot)
          endif
        enddo
        write(iw,9430) ctmall,ctdall
      endif
c     if(urospn.and.imect.eq.5.and.masout) then
      if(urospn.and.masout) then
        write(iw,9440)
        ijfg=0
        do ifg=1,nfg
          if(modfd.eq.0.or.(layfrg(ifg).gt.1.and.
     *       (iand(modfd,2).eq.0.or.iactfg(ifg).ne.0))) then
            do jfg=1,ifg-1
              ijfg=ijfg+1
              SCFi=scffrg(ifg)
              SCFj=scffrg(jfg)
              if(scfi.eq.uhf.or.scfi.eq.rohf.or.
     *           scfj.eq.uhf.or.scfj.eq.rohf)
     *        write(iw,9445) ifg,jfg,ctspin(ijfg)
            enddo
          endif
        enddo
      endif
c
c     if(nder.gt.0.and.modef.ne.1 .or.modfd.ne.0) then
      if(nder.gt.0.and.modef.ne.1.and.(runtyps.ne.optfmo.or.nbody.eq.2)
     *  .and.(.not.dftbfl.or.nbody.eq.2)) then
        if(hasgrad) then
          outgrdn=outgrd.and.(nbody.eq.2.or.iand(nprfmo,3).le.1)
          call fmogout(2,outgrdn,indat,fmozan,fmode,DEFTF,TORQF)
          if(resdim.ne.0.and.nbody.eq.3.and.outgrdn)
     *      write(iw,*) 'WARNING: FMO2 gradient has no ES dimer terms.'
c           They are lumped to FMO3 gradient.
C           FMO3-DFTB cannot calculate FMO2 gradient correctly...
c       else
c         write(iw,8113)
        endif
      endif
c     call fmogout(2,.true.,indat,fmozan,fmode,DEFTF,TORQF)
c
      if(nfg3.eq.0) goto 1000 
c     At present FMO3 does not work with BSSE
      if(masout) then
        write(iw,9250)
        if(.not.savemem) then
          if(needdm) then
            write(iw,9260)
          else
            if(ipieda.ne.0) then
              if(dopleda) then
                write(iw,9262)
              else
                dispstr='  Edi '
                if(do3c) dispstr='Edi+bs'
                if(dodft) dispstr='  Erc '
                if(gcorrel) dispstr='Erc+di'
                distr='Eex     Ect+mix'
                if(dftbfl) distr=' E0      Ect*es'
                write(iw,9270) distr,dispstr
              endif
            else
              write(iw,9265)
            endif
          endif
        endif
      endif 
      etrimtot1=zero
      etrimtot2=zero
      etrimtot2a=zero
      etrimtot1u=zero
      lijkfg=0
      rmin=-one
      rmax=-one
      rrij=zero
      rrik=zero
      rrjk=zero
      if(nappri.ne.0) call vclr(erapp,1,nappri*2)
      if(nappri.ne.0) call viclr(nrapp,1,nappri*2)
      sepmax=zero
      do iifg=1,nfg
        do jjfg=1,iifg-1
          if(fast3loop) then
            ijfg=((iifg-1)*(iifg-2))/2+jjfg
            if(loop3ij(ijfg).eq.0) then
              lijkfg=lijkfg+jjfg-1
              goto 900
            endif
          endif
          do kkfg=1,jjfg-1
            lijkfg=lijkfg+1
            if(fast3prop) then
               if(net.gt.0.and.lijkfg.gt.net) then
c                write(6,*) 'wwwdone all',net
                 goto 590
c                Break the trimer loop if all SCF trimers are done.
               endif
               ijkfg=loadt(lijkfg)
               call cubbrk(ijkfg,-1,ifg,jfg,kfg)
            else
               ijkfg=lijkfg
               ifg=iifg
               jfg=jjfg
               kfg=kkfg
            endif
            if(needr.ne.0) call fmodist3(ifg,jfg,kfg,rmin,rmax)
            sepmax=max(sepmax,rmax)
            ilay=layfrg(ifg)
            ijlay=min(ilay,layfrg(jfg))
            ijklay=min(ijlay,layfrg(kfg))
            memon=1
            metrim=1
            if(gcorrel) then
              ires=molfrg(ifg)+molfrg(jfg)+molfrg(kfg)
              if(exclmol) then
                skipc=ires.ge.0.and.ires.ne.3
              else
                skipc=ires.ge.0.and.ires.ne.1
              endif
              if(iand(modmol,8).ne.0) skipc=.false.
              if(rmax.gt.restri(4).and.restri(4).ne.0 .or. skipc) then
                 memon=2
                 metrim=nent
              endif
            endif
            if(ext2lay) then
              memon=2
              metrim=1
            endif
            memon1=memon
            if(iand(lcorrel(ijklay),1).eq.0) memon1=2
c           if(.not.gcorrel) memon1=2
c           Should also check memon1 vs ritrim(4)
            if(doexc.and.memon1.eq.1.and.(ifg.ne.iexcit(1).and.
     *         jfg.ne.iexcit(1).and.kfg.ne.iexcit(1))) memon1=2
            if(dodc) memon=2
            memon0=memon
            if(dopleda) memon=memon+2
            corri=corrs(memon1)
            emonijk=emon(ifg,memon,ijklay)+emon(jfg,memon,ijklay)+
     *              emon(kfg,memon,ijklay)
            emonijku=emonijk
            memon2=2
            if(dopleda) memon2=memon2+2
            if(gcorrel.or.dodc)
     *      emonijku=emon(ifg,memon2,ijklay)+emon(jfg,memon2,ijklay)+
     *              emon(kfg,memon2,ijklay)
            indij=(ifg*ifg-3*ifg)/2+jfg+1
            indik=(ifg*ifg-3*ifg)/2+kfg+1
            indjk=(jfg*jfg-3*jfg)/2+kfg+1
            esolvijk=0
            if(.not.savemem) then
              eijk=etrim(ijkfg,metrim)
c             
c             For trimers it is more difficult to identify divergence. 
              if(etrim(ijkfg,1).eq.one) then
                eijk=one
                ndiv3=ndiv3+1
              endif
              if(nfmopcm.ne.0) then
                esolvijk=esolv(nfg*2+nfg2*2+ijkfg)
                if(esolvijk.ne.zero) then
                  esolvijk=(esolvijk+esolv(ifg)+esolv(jfg)+esolv(kfg)
     *                     -esolv(nfg*2+indij)-esolv(nfg*2+indik)
     *                     -esolv(nfg*2+indjk))*gesfact
                  ges=ges+esolvijk
                  gsolv=gsolv+esolvijk*tokcal
                endif
c             write(6,*) 'wwwest',ifg,jfg,kfg,esolv(ifg)*tokcal,
c    *                   esolv(jfg)*tokcal,esolv(kfg)*tokcal,
c    *              esolv(nfg*2+indij)*tokcal,esolv(nfg*2+indik)*tokcal,
c    *                   esolv(nfg*2+indjk)*tokcal,
c    *                   esolv(nfg*2+nfg2*2+ijkfg)*tokcal
              endif
c             
              ddijk=etrim(ijkfg,2)
              didijk=eijk.ne.0.or..not.usedij
c                 if(ddijk.ne.0.or..not.usedij) then
            else
              eijk=0
              ddijk=0
              if(needr.ne.0.and.restri(1).ne.0.and.usedij) then
c               This is a non-trivial piece.
c               It should agree with ETRIMER.
                didijk=.not.(rmin.gt.restri(1).and.rmax.gt.restri(2).or.
     *                       rmax.gt.restri(3).and.restri(3).ne.0)
              else
                didijk=.true.
              endif
              if(nfmopcm.ne.0.and.didijk) then
                esolvijk=(esolv(ifg)+esolv(jfg)+esolv(kfg)
     *                   -esolv(nfg*2+indij)-esolv(nfg*2+indik)
     *                   -esolv(nfg*2+indjk))*gesfact
                ges=ges+esolvijk
                gsolv=gsolv+esolvijk*tokcal
                esolvijk=0
              endif
            endif
            if(didijk) then
              if(eijk.eq.one.and.ndualb.ne.0) then
c                if(maswrk) write(iw,*) 'Trimer',ifg,jfg,kfg,' diverged'
c                call abrt 
              endif
c             if(maswrk) write(iw,*) 'wwwTrimer',ifg,jfg,kfg
              if(eijk.eq.one) eijk=zero
c             deijk=eijk+emonijk-edim(indij,medimij)-edim(indik,medimik)-
c    *                           edim(indjk,medimjk)
c             Actually one only has to set memij,memik,memjk properly.
c             It makes no difference with dimers (medij,medik,medjk) because 
c             for them uncorrelated values are copied to the correlated location
c             so if one accesses the wrong place by mistake there is no 
c             problem. But monomers always (well, almost) have different 
c             corr/uncorr values so proper indices memij,memik,memjk must be 
c             computed.
              inisetd=1
              inisetm=1
              if(dodc) inisetd=3
              if(dodc) inisetm=2
c             For DC, use energies without Edisp, i.e., enforce no three-body
c             effects for Edisp, which is proper. 
              medij=inisetd
              medik=inisetd
              medjk=inisetd
              memij=inisetm
              memik=inisetm
              memjk=inisetm
c             The code below has not yet been changed to accomodate exclmol.
              if(gcorrel) then
                if(resdim.ne.0) then
                  rrij=fmodist(ifg,0,0,jfg)
                  rrik=fmodist(ifg,0,0,kfg)
                  rrjk=fmodist(jfg,0,0,kfg)
                endif
                ires=molfrg(ifg)+molfrg(jfg)
                if(nbody.eq.3) ires=1
                if(memon0.eq.2.or.(rrij.gt.rcorsd.and.rcorsd.ne.0).or.
     *             ires.eq.0.or.ires.eq.2) then
                  medij=3
                  memij=2
                endif
                ires=molfrg(ifg)+molfrg(kfg)
                if(nbody.eq.3) ires=1
                if(memon0.eq.2.or.(rrik.gt.rcorsd.and.rcorsd.ne.0).or.
     *             ires.eq.0.or.ires.eq.2) then
                  medik=3
                  memik=2
                endif
                ires=molfrg(jfg)+molfrg(kfg)
                if(nbody.eq.3) ires=1
                if(memon0.eq.2.or.(rrjk.gt.rcorsd.and.rcorsd.ne.0).or.
     *             ires.eq.0.or.ires.eq.2) then
                  medjk=3
                  memjk=2
                endif
              endif
              if(dopleda) then
                memij=memij+2
                memik=memik+2
                memjk=memjk+2
              endif
              lowtrimer=fullmfmo3.and.ijklay.eq.1
              deij=edim(indij,medij)-
     *               emon(ifg,memij,ijklay)-emon(jfg,memij,ijklay)
              deik=edim(indik,medik)-
     *               emon(ifg,memik,ijklay)-emon(kfg,memik,ijklay)
              dejk=edim(indjk,medjk)-
     *               emon(jfg,memjk,ijklay)-emon(kfg,memjk,ijklay)
c             write(6,*) 'wwww',jfg,kfg,edim(indjk,medjk),
c    *         emon(jfg,memjk,ijklay),emon(kfg,memjk,ijklay),medjk,memjk
              if(lowtrimer) then
                deij=deij-edim(indij,medij)+edimlow(indij,medij)
                deik=deik-edim(indik,medik)+edimlow(indik,medik)
                dejk=dejk-edim(indjk,medjk)+edimlow(indjk,medjk)
              endif
              deijk=eijk-emonijk-deij-deik-dejk
c
              if(lowtrimer) then
          ddijk=ddijk-edimlow(indij,2)-edimlow(indik,2)-edimlow(indjk,2)
              else
                ddijk=ddijk-edim(indij,2)-edim(indik,2)-edim(indjk,2)
              endif
              etrimtot1=etrimtot1+deijk
              if(gcorrel) then
                if(savemem) then
                  eijku=0
                else
                  eijku=etrim(ijkfg,nent)
                endif
                if(lowtrimer) then
                  deijku=eijku+emon(ifg,memon2,ijklay)+
     *                   emon(jfg,memon2,ijklay)+emon(kfg,memon2,ijklay)
     *               -edimlow(indij,3)-edimlow(indik,3)-edimlow(indjk,3)
                else
                  deijku=eijku+emon(ifg,memon2,ijklay)+
     *                   emon(jfg,memon2,ijklay)+emon(kfg,memon2,ijklay)
     *                  -edim(indij,3)-edim(indik,3)-edim(indjk,3)
                endif
                etrimtot1u=etrimtot1u+deijku
              else
                deijku=deijk
                if(dodc) etrimtot1u=etrimtot1u+deijk
              endif
              if(dopleda) then
c               1 comes from E'I-E'J-E'K
                npli=1
                nplj=1
                nplk=1
c               Add -1 for SCF dimers from E'IJ-E'I-E'J
                if(resdim.ne.0) then
                  if(rrij.le.resdim) npli=npli-1
                  if(rrij.le.resdim) nplj=nplj-1
                  if(rrik.le.resdim) npli=npli-1
                  if(rrik.le.resdim) nplk=nplk-1
                  if(rrjk.le.resdim) nplj=nplj-1
                  if(rrjk.le.resdim) nplk=nplk-1
                else
                  npli=npli-2
                  nplj=nplj-2
                  nplk=nplk-2
                endif
                epl0ijk=npli*(epl0ds(ifg,1)+epl0ds(ifg,2))
     *                 +nplj*(epl0ds(jfg,1)+epl0ds(jfg,2))
     *                 +nplk*(epl0ds(kfg,1)+epl0ds(kfg,2))
c               In PL0 runs polarisation is double counted.
              else
                epl0ijk=zero
              endif
              if(domdan) x(letrimmd+ijkfg-1)=deijk+ddijk+esolvijk
              if(autosa.and.isasign.eq.1) then
                etrimval=(deijk+ddijk+esolvijk)/3
                ijfg3=(ifg*ifg-3*ifg)/2+jfg+1
                ikfg3=(ifg*ifg-3*ifg)/2+kfg+1
                jkfg3=(jfg*jfg-3*jfg)/2+kfg+1
                sadim(ijfg3)=sadim(ijfg3)+iapsign*etrimval
                sadim(ikfg3)=sadim(ikfg3)+iapsign*etrimval
                sadim(jkfg3)=sadim(jkfg3)+iapsign*etrimval
              endif
c             if(masout) then
              if(needdm) then
c               Use V-consistent ddDVijk Tr(delta-delta-Density-ijk * V-ijk).
c               However, the 3-body effect that is printed as etot is computed
c               from the pure 3-body quantity ddVijk. This might be confusing,
c               but the other possibility dddijk might be equally so.
c             In the ddDVijk formula one should ignore dDij and use only dddijk;
c             but then the notation of purely three-body density effect is lost.
                dddijk=etrim(ijkfg,3) 
                detot=deijk+ddijk
                etrimtot2=etrimtot2+dddijk
                etrimtot2a=etrimtot2a+ddijk
                if(masout)
     *          write(iw,9275) ifg,jfg,kfg,corri,ijklay,rmin,rmax,eijk
     *                        ,deijk,deijku,ddijk,dddijk,esolvijk*tokcal
     *                        ,(detot+esolvijk)*tokcal
              else
                detot=deijk+ddijk
                etrimtot2=etrimtot2+ddijk
                if(ipieda.ne.0) then 
                  etrimdi=deijk-deijku
                  if(dodft) then
c                   Subtract correlation energy from E'I for each monomer.
                    etrimdi=ecorrdft(ijkfg+nfg+nfg2)-ecorrdft(indij+nfg)
     *                     -ecorrdft(indik+nfg)-ecorrdft(indjk+nfg)
     *                     +ecorrdft(ifg)+ecorrdft(jfg)+ecorrdft(kfg)
                  endif
                  if(dftbfl) then
c                   Calculate DeltaE'ES
                   e1es=ecorrdft(ijkfg+nfg*2+nfg2)-ecorrdft(indij+nfg*2)
     *                 -ecorrdft(indik+nfg*2)-ecorrdft(indjk+nfg*2)
     *                 +ecorrdft(ifg)+ecorrdft(jfg)+ecorrdft(kfg)
c                   Calculate DeltaE(CT*ES) = DeltaE'ES+Delta E_V
                    etrimct=e1es+ddijk
c                   Calculate DeltaE0 as the remainder 
                    etrimex=detot-etrimct
                  else
                    eesexij=edim(indij,nedimes)+edim(indij,nedimex)
                    eesexik=edim(indik,nedimes)+edim(indik,nedimex)
                    eesexjk=edim(indjk,nedimes)+edim(indjk,nedimex)
                    etrimex=etrim(ijkfg,3)-emonijku
     *                     -eesexij-eesexik-eesexjk
                    if(nfmopcm.ne.0) etrimex=etrimex
     *                                 -esolv(ifg)-esolv(jfg)-esolv(kfg)
                    if(dodft) etrimex=etrimex+ecorrdft(ifg)
     *                               +ecorrdft(jfg)+ecorrdft(kfg)
c                   For connected trimers in APC, EX is not defined.
                    if(doapc.and.etrim(ijkfg,3).eq.0) etrimex=0
                    etrimct=detot-etrimex-etrimdi-epl0ijk
                    if(dopleda) then
                      if(masout.and..not.savemem)
     *                write(iw,9285) ifg,jfg,kfg,corri,ijklay,rmin,rmax,
     *                               deijk*tokcal,ddijk*tokcal,
     *                               (detot+esolvijk)*tokcal,
     *                               epl0ijk*tokcal,etrimex*tokcal,
     *                               etrimct*tokcal,etrimdi*tokcal
c                     PL0 is not interfaced with PCM, so esolvijk is not done.
                    endif
                  endif
                  if(masout.and..not.savemem.and..not.dopleda)
     *              write(iw,9285) ifg,jfg,kfg,corri,ijklay,rmin,rmax,
     *                             deijk*tokcal,ddijk*tokcal,
     *                             (detot+esolvijk)*tokcal,
     *                             etrimex*tokcal,etrimct*tokcal,
     *                             etrimdi*tokcal,esolvijk*tokcal
c                  After printing, save 3-body components.
                   if(dftbfl) then
                     ecorrdft(ijkfg+nfg*2+nfg2)=etrimct
                     etrim(ijkfg,2)=etrimex
                   else
                     etrim(ijkfg,3)=etrimex
                     if(dodft) ecorrdft(ijkfg+nfg+nfg2)=etrimdi
                     if(gcorrel) etrim(ijkfg,nent)=etrimdi
                     etrim(ijkfg,2)=etrimct
                   endif
                   if(nfmopcm.ne.0) esolv(nfg*2+nfg2*2+ijkfg)=esolvijk
                   etrim(ijkfg,1)=1
c                  Setting to +1 as a stupid but simple way to check if IJK was 
c                  processed here.
                else
                  if(masout.and..not.savemem)
     *            write(iw,9280) ifg,jfg,kfg,corri,ijklay,rmin,rmax,eijk
     *                          ,deijk,deijku,ddijk,esolvijk*tokcal,
     *                           (detot+esolvijk)*tokcal
                  if(.not.savemem.and.looppieda3.eq.0) then
c                   etrim is overwritten for effective PIE
                    etrim(ijkfg,1)=deijk
                    etrim(ijkfg,2)=ddijk
                    if(gcorrel) etrim(ijkfg,3)=deijku
                    if(nfmopcm.ne.0) esolv(nfg*2+nfg2*2+ijkfg)=esolvijk
                  endif
                endif
              endif
c             endif
              if(nappri.ne.0) then
                if(rmax.le.rappri(1)) then
                  iappri=1
                else
                  iappri=int((rmax-rappri(1))/rappri(3))+2
                endif
                if(iappri.le.nappri) then
                  erapp(iappri,1)=erapp(iappri,1)+deijk+ddijk
                  erapp(iappri,2)=erapp(iappri,2)+deijku+ddijk
                  nrapp(iappri,1)=nrapp(iappri,1)+1
                  if(memon1.eq.1) nrapp(iappri,2)=nrapp(iappri,2)+1
                endif
              endif
            endif
          enddo
  900     continue
        enddo
      enddo
  590 continue
c
c     Add the lump contributions from all trimers.
c
      if(savemem) then
        etrimtot1=etrimtot1+etrimsum(1)
        etrimtot2=etrimtot2+etrimsum(2)
        if(gcorrel) then
          etrimtot1u=etrimtot1u+etrimsum(3)
        else
          etrimtot1u=etrimtot1u+etrimsum(1)
c         For D, 3-body corrections are zero, so we use no D value here.
        endif
        if(nfmopcm.ne.0) then
          ges=ges+etrimsum(7)*gesfact
          gsolv=gsolv+etrimsum(7)*gesfact*tokcal
        endif
      endif
c
      etota=etot+etrimtot1+etrimtot2a
      etotau=etotu+etrimtot1u+etrimtot2a
c     if(needr2.ne.0) then
c     As mentioned above, the ddDVijk equation should have two-body dDij
c     subtracted as it was added in the two-body sum.
      if(needdm) then
        etot=etot-eddim
        etotu=etotu-eddimu
      endif
      etotsav=etot
      etotusav=etotu
      etrimtot=etrimtot1+etrimtot2
      etrimtotu=etrimtot1u+etrimtot2
      etot1=etot+etrimtot1
      etot1u=etotu+etrimtot1u
      etot=etot+etrimtot
      etotu=etotu+etrimtotu
      if(ndiv3.ne.0.and.runtyp.ne.check) then
        write(iw,9500) ndiv3,3
        etot=0
        etotu=0
      endif
      write(iw,9200) symbe,corstr(1),3,etot
      if(nfmopcm.ne.0.and.mdoutmin) then
        write(iw,9200) symbe,esstr,3,etot+ges
        write(iw,9200) symbe,solstr,3,etot+gsolv/tokcal
      endif
      if(needdm) then
c       needdm prints three energies:
c       Euncorr is the EI'+dE'IJ+dE'IJK+ddDijkVijk
c       E no D  is the EI'+dE'IJ+dE'IJK
c       E D23   is the EI'+dE'IJ+dE'IJK+dDijVij
c                         +(dDijkVijk-dDijVij-dDikVik-dDjkVjk)
c               where dDijk=Dijk-(Di+Dj+Dk).
c       "E no D" is useless, Euncorr is proper only for no ritrim approximation.
c       "Euncorr" is equal to "E no D" without approximations.
c       "E no D" is the most sound one at present. It is used as Euncorr for 
c       needdm=.false. runs (production).
c       Its drawback is a fairly high sensitivity to ES approximations,
c       as it is based on Vij and Vijk.
c
        write(iw,9200) symbe,' no D ',3,etot1
        write(iw,9200) symbe,' D23  ',3,etota
      endif
c     if(douhf.and.mdoutmin) 
c    *  write(iw,9780) 3,sz0,3,szfmo(3)*(szfmo(3)+1),s2_0
c
c     Note that 3-body corrections to /D are zero, so none are accumulated
c     or printed (0 because they are pair additive). However, we print for
c     consistency some FMO3 properties.
c
      if(gcorrel.or.dodc) then
        write(iw,9200) symbe,corstr(2),3,etotu
        write(iw,9200) symbe,corstr(3),3,etot-etotu
        if(needdm) then
          write(iw,9200) symbe,'uncnoD',3,etot1u
          write(iw,9200) symbe,'uncD23',3,etotau
        endif
        if(didcc) call fmoccp(3,etotu,extracc(1,3))
        if(didmp2) call fmompp(3,etotu,extracc(1,3))
        if(do3c) call fmoh3c(3,extracc(1,3))
      endif
c
      if(nappri.ne.0) then
        do i=2,nappri 
          erapp(i,1)=erapp(i,1)+erapp(i-1,1)
          erapp(i,2)=erapp(i,2)+erapp(i-1,2)
          nrapp(i,1)=nrapp(i,1)+nrapp(i-1,1)
          nrapp(i,2)=nrapp(i,2)+nrapp(i-1,2)
        enddo
        write(iw,*) ' '
        do i=1,nappri 
          riapp=rappri(1)+rappri(3)*(i-1)
          if(restri(3).ne.0) sepmax=min(sepmax,restri(3))
          if(riapp.gt.sepmax.and.riapp-rappri(3).le.sepmax) riapp=sepmax
          if(riapp.le.sepmax.and.sepmax.ne.0) then 
          if(gcorrel) then
           write(iw,9205) riapp,nrapp(i,2),nrapp(i,1),3,
     *       erapp(i,1)+etotsav-erapp(i,2)-etotusav,erapp(i,2)+etotusav
          else
           write(iw,9206) riapp,nrapp(i,1),3,erapp(i,1)+etotusav
          endif
          endif
        enddo
        write(iw,*) ' '
      endif
c
c     write(iw,9200) ' + D2 ',3,etot1+eddim
      if(nfmoelm.ne.0) 
     *  call daxpy(nfmoelm,one,fmoelm(1,nfg+1),1,fmoelm(1,nfg+2),1)
      if(iemom.ge.1) 
     *  write(iw,9230) 3,(fmoelm(ielm,nfg+2),ielm=nelm(0)+1,nelm(1)),
     *                  DNRM2(3,fmoelm(1,nfg+2),1)
      if(iemom.ge.2) 
     *  write(iw,9231) (fmoelm(ielm,nfg+2),ielm=nelm(1)+1,nelm(2))
      if(iemom.ge.3) 
     *  write(iw,9232) (fmoelm(ielm,nfg+2),ielm=nelm(2)+1,nelm(3))
      if(nfmoelm.ne.0) call dcopy(nfmoelm,fmoelm(1,nfg+2),1,DMX,1)
c     Do not output the gradient for optimisation (it will be done elsewhere)
c
      if(nder.gt.0.and.modef.ne.1.and.(runtyps.ne.optfmo.or.nbody.eq.3))
     *then
        if(hasgrad) then
          call fmogout(3,outgrd,indat,fmozan,fmode,DEFTF,TORQF)
c       else
c         write(iw,8113)
        endif
      endif
      IF(pcmprp.and.mdoutmin) write(iw,9355) 3,gsolv,3,ges*tokcal
c
c     if(nbody.eq.3.and.ipieda.gt.0.and.looppieda3.eq.1) then
      if(nbody.eq.3.and.looppieda3.eq.0) then
c       a really mean jump to reprint pair values
        looppieda3=1
       call pieda3(ipieda,dftbfl,dodft,gcorrel,dodc.or.do3c,nfmopcm.ne.0
     *           ,nfg,nfg2d,nfg2,nfg3,nedimex,edim,etrim,ecorrdft,esolv)
        symbe='E'
        goto 1100
      endif
c
 1000 continue
c     
c     if(totfock) then
c       write(iw,9300) e1efmo,ekinfmo,etot-e1efmo-enucr
c     ekinfmo can be zero if FMO1 has not converged. 
      dovir=ekinfmo.ne.0.and.(.not.doexc.or.iexcit(3).ne.0).and.
     *      .not.urohf.and..not.fullmfmo3.and.needmd.gt.0.and.
     *      .not.dopbcmd.and.mdoutmin.and..not.ignoresd.and.
     *      .not.(dofret.or.dofed)
      if(dovir) then
        etotnod=etot
        if(dodc) etotnod=etotu
c       Do not include empiric dispersion into this analysis,
c       which is to be done for QM only.
        E1=e1efmo
        TKIN=ekinfmo
        E2 = ETOTnod - E1 - ENUCR
        VNE = E1 - TKIN
        VNN = ENUCR
        VEE = E2
        VTOT = VNE + VNN + VEE
        VIRIAL = -VTOT/TKIN
        WRITE(IW,9300) E1,E2,ENUCR,ETOTnod
        WRITE(IW,9305) VEE,VNE,VNN,VTOT,TKIN,VIRIAL
      endif
      if(totfock.and.meglob.eq.0) then
        write(ip,9310) enucr
        write(ip,9315) (int(fmozan(i)+0.01D+00),i=1,natfmo)
      endif
c
      if(nfmopcm.gt.0.and.iemom.ge.1.and.mdoutmin) then
        pcmq=ddot(nfg,emocdr(7),10,one,0) 
        pcmdx=ddot(nfg,emocdr(8),10,one,0) 
        pcmdy=ddot(nfg,emocdr(9),10,one,0) 
        pcmdz=ddot(nfg,emocdr(10),10,one,0) 
        pcmd=sqrt(pcmdx*pcmdx+pcmdy*pcmdy+pcmdz*pcmdz)
        idipn=nbody-1
        if(nbody.eq.1) idipn=3
        totdx=pcmdx+fmoelm(1,nfg+idipn)
        totdy=pcmdy+fmoelm(2,nfg+idipn)
        totdz=pcmdz+fmoelm(3,nfg+idipn)
        totd=sqrt(totdx*totdx+totdy*totdy+totdz*totdz)
        WRITE(IW,932) pcmq,pcmdx,pcmdy,pcmdz,pcmd,totdx,totdy,totdz,totd
      endif
      IF(pcmprp) then
        GCAV=zero
        IF(ICAV.EQ.1) gcav=gcavp 
        IF(IDISP.EQ.0) THEN
          GDISP=ZERO
          GREP=ZERO
        ELSE
          GREP=GRP
        END IF
        IF(ISMX.NE.0) then
          gint=ges+GCDS/tokcal
          etot=etot+gint
          if(gcorrel.or.dodc) etotu=etotu+gint
          if(mdoutmin) write(iw,9380) GES*tokcal,GCDS,GINT*tokcal
        else
          gint=ges+(gcavp+gdisp+grep)/tokcal
          etot=etot+gint
          if(gcorrel.or.dodc) etotu=etotu+gint
          if(mdoutmin) write(iw,9390) GES*tokcal,GCAV,GDISP,GREP,
     *                                GINT*tokcal
        endif
        write(iw,9395) corstr(1),etot,corstr(1),etot-gint
        if(mdoutmin.and.(gcorrel.or.dodc.and..not.reducee)) 
     *    write(iw,9395) corstr(2),etotu,corstr(2),etotu-gint
      END IF
c
      if(ifgfmo0.eq.0) then
        escf=etot
        if(nend.ge.7) escf=etotu
      endif 
c
      if(ndftfg.eq.1.and.mdoutmin) then
        write(iw,9450) (etotdft(i),i=1,nedft)
      endif
c     if((primul.or.urospn).and.nbody.gt.0
c     For modmol=1 or 2, there is a bookkeeping problem.
c     Perhaps, only molfrg charges are correct?!
      nbodyc=nbody
      if(primul.and.nbody.gt.0.and.(.not.doexc.or.iexcit(3).gt.0).and.
     *  (masout.or.iand(modpan,1).ne.0).and.iand(modmol,3).eq.0) then
c       charges are not summed properly in PL0
        if(iand(modcha,3).ne.0.and..not.autosa) nbodyc=iand(modcha,3)
        if(dofret.or.dofed) nbodyc=1 
c       if(dopleda) nbodyc=1 
c       Choose the charge meister.
        charger='Mulliken'
        if(NPTSTN.ne.0) charger='Stone'
        if(dopdc) charger='Potential'
        if(nfmopcm.eq.0) then
          write(iw,9400) charger
        else
          write(iw,9402) charger
        endif
        do i=1,natfmo
          if(.not.skipscc1.or.iactat(i).ne.0) then
            if(nbody.gt.1) atmulq(i,2)=atmulq(i,2)+atmulq(i,1)
            if(nbody.gt.2) atmulq(i,3)=atmulq(i,3)+atmulq(i,2)
            zi=fmozan(i)
            if(dftbfl) call dftb_nucz(zi,2)
            if(IMCPFMO.EQ.1) zi=zi-fzcor(I)
            if(nfmopcm.eq.0) then 
              write(iw,9410) i,indat(i),fmozan(i),
     *                       (zi-atmulq(i,j),j=1,nbodyc)
            else
              write(iw,9412) i,indat(i),fmozan(i),ascat(1,i),
     *                       (ascat(1,i)/ascat(2,i))*1.0D+02,ascat(3,i),
     *                       (zi-atmulq(i,j),j=1,nbodyc)
            endif
            if(dosczvprop) sczvprp(i)=zi-atmulq(i,nbodyc)
          endif
        enddo
        if(nbody.gt.2.and.ipieda.ne.0) call prifragq(IW,dftbfl,nbodyc,
     *                              IMCPFMO,indat,atmulq,fmozan,fzcor)
        if(iand(nprfmo,3).eq.0.and.maswrk.and.meglob.eq.0.and.
     *     icurpop.ge.1.and.icurpop.le.2) then
          write(ip,9480) charger,maxnat,nfg
          write(ip,9481) ((popmat(i,j,icurpop),i=1,maxnat),j=1,nfg)
        endif
      endif
      if(urospn.and.masout) then
c       write(iw,9403) charger
        write(iw,9403) 'Mulliken'
        do i=1,natfmo
          write(iw,9410) i,indat(i),fmozan(i),
     *                     (pmulspin(i,j),j=1,min(nbody,2))
        enddo
      end if
      if(iahard.gt.0) then
        write(iw,9470)
        do i=1,natfmo
          write(iw,9475) i,indat(i),fmozan(i),
     *                   (fmohard(i,j)/nefmo,j=1,nbody)
c         The factor of 2 is not included.
        enddo
      endif
c     The end of if(maswrk).
      endif
c
c     Relaxed density corrections 
c
      if(dosczvprop) then
        if(masout) write(iw,9520) nbodyc
        do i=1,natfmo
          qi=sczvprp(i)
          deltaqi=sczvprp(natfmo+i)
          if(masout) write(iw,9530) i,indat(i),fmozan(i),deltaqi,
     *                              qi+deltaqi
        enddo
        if(masout) write(iw,9540) 
        do i=1,nfg
          ind=natfmo*2+(i-1)*3
          if(masout) write(iw,9550) i,(sczvprp(ind+j),j=1,3),
     *                                (fmoelm(j,i)+sczvprp(ind+j),j=1,3)
          dmx=dmx+sczvprp(ind+1)
          dmy=dmy+sczvprp(ind+2)
          dmz=dmz+sczvprp(ind+3)
        enddo
        if(masout) write(iw,9560) dmx,dmy,dmz,
     *                            sqrt(dmx*dmx+dmy*dmy+dmz*dmz)
      endif
c
c     Slaves have no multipoles which may lead to very bizarre problems.
      if(goparr.and.nfmoelm.ne.0)
     *  CALL DDI_BCAST(2411,'F',DMX,nfmoelm,MASTER)
c
c     Print FMO-TDDFT RESULTS
c
      if(doexc) then
        if(dofret.or.dofed.or.domipea) then
          call excout2(excit2,excit3,texcit2,osmd,iexcit,dotd,doci,edim,
     *                 iactfg,nfgfret,ifgfret,dofed,modprp,doeom,ipeam,
     *                 domipea,frgnam,dofret2,hfretap,ndualb,dadb)
        else
          call excout(frgnam,eexcit,texcit,osmd,iexcit,ipeam,nstmono,
     *         isumd,isumt,eexfg,needr.ne.0,dotd,doci,doeom,extri,exfid)
        endif
      endif
c
c     Now add CNS terms.
c
      if(docns) then
        etot0=etot
        etot=etot+exref
        if(maswrk) write(iw,9750) etot0,exref,etot
        if(nbody.gt.0) then
c         print QM gradient 
          call fmogout(nbody,outgrd,indat,fmozan,fmode,DEFTF,TORQF)
c         print Wa*Exref gradient 
c         For FMO/EFP this will print some unwanted EFP terms.
          if(masout) write(iw,9755)
          call fmogout(1,outgrd,indat,fmozan,cnsdat,DEFTF,TORQF)
          call daxpy(3*natfmo,one,cnsdat,1,fmode(1,1,nbody),1)
c         Add them. The total gradient will normally be printed elsewhere.
          call fmogout(nbody,outgrd,indat,fmozan,fmode,DEFTF,TORQF)
        endif
      endif
c
c     Zero out the energy for diverged runs
c
      if(.not.convSCC.or.ndiv2+ndiv3.ne.0) then
        etot=zero 
        escf=zero
      else
        if(modfd.ne.0) then
          if(rmixdim.eq.1.0D+30) rmixdim=0.0D+00
c         That means that F domain is empty...
          write(iw,9700) nevsav,emixdim,rmixdim,imixdim,jmixdim
        endif
      endif
c
c     Now we broadcast the results to slaves: ESCF, ETOT, ENUCR.
c
      IF(GOPARR) then
        dumbuf(1)=etot
        dumbuf(2)=escf
        dumbuf(3)=enucr
        CALL DDI_BCAST(2421,'F',dumbuf,3,MASTER)
        etot=dumbuf(1)
        escf=dumbuf(2)
        enucr=dumbuf(3)
      endif
      if(domdan) call icopy(natfmo,indat,1,x(lindmd),1)
      e=etot
      if(ndualb.ne.0) escf=etotu
      if(maswrk.and.ndualb.eq.0) write(iw,9810) e
      if(autosa) then
        esatot=esatot+isasign*iapsign*e
c       for 1 fragment, SA is skept, so assign the fragment energy
        ifgsa=1+nfgasub 
        if(nbody.lt.2) samon(1,ifgsa)=samon(1,ifgsa)+isasign*iapsign*e
      endif
      if(ndualb.ne.0) epcmap(2)=epcmap(2)+iapsign*e
      if(mdoutmin) then
      if(maswrk) write(iw,9900)
      call timit(1)
      endif
      RETURN
  932 FORMAT(/1X,'PCM CHARGE=',F12.7,/1X,'PCM DIP=',4F13.7,' (DEBYE)',
     *       /1X,'TOT DIP=',4F13.7,' (DEBYE)')
 9100 format(//9x,'FMO properties',/9x,15(1H#),//)
 9106 format(1x,'No monomer SCF done in layer',I2)
 9108 format(1x,'Layer',I2,' DIVERGED (max iterations',I3,', conv ',
     *          E7.2,').')
 9109 format(1x,'Layer',I2,' converged in',I3,' iterations (max',I3,
     *          ', conv ',E7.2,').',/)
 9110 format(/1x,'Energy values (E) are given with ',
     *           'external ESP subtracted.')
 9111 format(/1x,'D(dipole), Q(quadrupole) and O(octopole) moments are',
     *           ' printed in the order:',
     *       /1x,'X,Y,Z; XX,YY,ZZ,XY,XZ,YZ; XXX,XXY,XXZ,XYY,YYY,YYZ,',
     *           'XZZ,YZZ,ZZZ,XYZ,',
     *       /1x,'computed at the point (A): ',3F13.8)
 9112 format(/1x,
     *         'Both uncorrelated and correlated energies are printed.',
     *      /1x,'Gradient is printed only as correlated values. ',
     *      /1x,'Other properties correspond to',A3,'correlated ',
     *          'wavefunction.')
 9113 format(1x,'For FMO-CC all values labelled "correlated" are for ',
     *          'the "best gradient" CC level.')
 8113 format(/1x,'Gradient is not calculated for RUNTYP=FMOHESS. ',
     *           'Use RUNTYP=GRADIENT to get it.')
 9114 format(/1x,'One-body FMO properties.',/1x,23(1H=)/)
 9115 format(26x,'E"corr',10x,'E"uncorr',7x,'DX',8x,'DY',8x,'DZ')
 9116 format(28x,'E"',10x,'DX',8x,'DY',8x,'DZ')
 8115 format(22x,'EPL0d',4x,'EPL0s',4x,'EPL0',4x,'EPL0',A2,2x,'EPLd*PLd'
     *      ,2x,'EPLd*',a2,4x,'DX',8x,'DY',8x,'DZ')
 8116 format(22x,'EPL0d',4x,'EPL0s',4x,'EPL0',3x,'EPLd*PLd',
     *        4x,'DX',8x,'DY',8x,'DZ')
 8117 format(22x,'EPL0d',4x,'EPL0s',4x,'EPL0',4x,'EPL0',a2,5x,
     *           'DX',8x,'DY',8x,'DZ')
 8118 format(1x,I4,'(',A8,',L',I1,')',4F9.3,3F10.5)
 9117 format(1x,I4,'(',A8,',L',I1,')',6F9.3,3F10.5)
 9118 format(1x,I4,'(',A8,',L',I1,')',2F17.9,3F10.5)
 9119 format(1x,I4,'(',A8,',L',I1,')',4F9.3,3F10.5)
 9120 format(1x,I4,'(',A8,',L',I1,')', F17.9,3F10.5)
 8122 format(8x,'Q=',6F10.5) 
 8123 format(8x,'O=',5F10.5,/8x,2x,5F10.5)
 8130 format(1x,'total = (EIJ-EI-EJ)+dDIJ*VIJ+Gsol = ',
     *          'Ees+Eex+Ect+mix+Edisp+Gsol (kcal/mol).')
 8131 format(/1x,'Contracted two-body FMO properties (with FMO3 ',
     *           'corrections).',/1x,59(1H=))
 9121 format(/1x,'Two-body FMO properties.',/1x,24(1H=)//,
     *        1x,'DL: D=C dynamically correlated (MP2,CI), D=N not ',
     *          'dynamically correlated',/1x,'(RHF,DFT). D=S separated',
     *          ' dimer: semiclassical ','interaction (ES), D=M MCSCF.',
     *       /1x,'L stands for layer, Z is the monomer charge product,',
     *           ' R is the interfragment',/1x,'distance relative to ',
     *         'van-der-Waals radii (-1.00 is printed if distances are',
     *       /1x,'not computed). dDIJ*VIJ is the explicit embedded ',
     *           'charge transfer energy.',/1x,'Q(I->J) is the charge ',
     *           'transfer amount, printed as zero if not available.',
     *       /1x,'Positive values correspond to I in IJ having extra ',
     *           'negative charge.')
 9122 format(/4x,'I    J DL  Z    R   Q(I->J)',8x,'E"',8x,
     *       'E"IJ-E"I-E"J   dDIJ*VIJ',6x,'Gsol',5x,'tot',/1x,90(1H-))
 9123 format(/4x,'I    J DL  Z    R   Q(I->J)',8x,'E"corr',10x,
     *       'E"uncorr',4x,'E"IJ-E"I-E"J,corr/uncorr dDIJ*VIJ,unc',3x,
     *       'Gsol',2x,'tot,corr',/1x,120(1H-))
 8120 format(/3x,'I   J  IU  SF DL  Z    R   Q(I->J)',8x,'E',11x,
     *       'EIJ-EI-EJ   dDIJ*VIJ',6x,'tot',/1x,87(1H-))
 8121 format(/3x,'I   J  IU  SF DL  Z    R   Q(I->J)',8x,'Ecorr',12x,
     *      'Euncorr',6x,'EIJ-EI-EJ,corr/uncorr  dDIJ*VIJ,unc tot,corr',
     *       /1x,118(1H-))
 9125 format(/4x,'I    J DL  Z    R   Q(I->J)  EIJ-EI-EJ dDIJ*VIJ',4x,
     *        'total',5x,'Ees ',5x,'Eex ',3x,'Ect+mix ',2x,A23,
     *       /1x,105(1H-))
 9124 format(/4x,'I    J DL  Z    R   Q(I->J)  EIJ-EI-EJ dDIJ*VIJ',4x,
     *        'total',5x,'Ees ',5x,' E0 ',3x,' Ect*es ',3x,'Edisp ',
     *       3x,'Gsol',/1x,105(1H-))
 8125 format(/4x,'I    J DL  Z    R   Q(I->J)  EIJ-EI-EJ dDIJ*VIJ',4x,
     *        'total',3x,'EPL0I+J',4x,'Ees0',5x,'Eex0',4x,'Ectmix0',3x,
     *        'E',A5,/1x,106(1H-))
 8150 format(/1x,'Frontier molecular orbital (FMO!) properties based ',
     *           'on Koopmans'' theorem.',
     *       /1x,'Electronegativity = - chemical potential. ',
     *           'All values are in eV.',
     *      //1x,' Frag  Ionization    Electron    Chemical   ',
     *           ' Global       Global  Electrophili-',
     *       /1x,'        potential    affinity    potential  ',
     *           'hardness     softness  city index',
     *       /1x,'EA and related properties are often very inaccurate.',
     *       /1x,79(1H-))
 8160 format(1x,I5,5F12.6,F14.6)
 8170 FORMAT(/1X,'FMO-EFP energies for each FMO fragment I: ',
     *           'DI*VI+0.5*SUM_J{dDIJ*VIJ}')
 9126 format(/1x,'Notation: BSSE corr. is a sum of int columns 3+4-1-2',
     *        8x,'/--- mon basis ---\\ /+++ dim basis +++\\ ',4x,
     *        'I   J DL',10x,'E',9x,'EIJ-EI-EJ   dDIJ*VIJ    EIJ-EI-EJ',
     *           ' dDIJ*VIJ BSSE corr.',4x,'tot',/1x,93(1H-))
 9128 format(/4x,'I   J DL    R',11x,'E',9x,'EIJ-EI-EJ corr/uncorr ',
     *       'dDIJ*VIJ',4x,'BSSE corr/uncorr',4x,'tot',/1x,88(1H-))
 9129 format(/4x,'I   J DL    R',11x,'E',9x,'EIJ-EI-EJ   dDIJ*VIJ',4x,
     *       'BSSE corr.',4x,'tot',/1x,79(1H-))
 9131 format(2I5,1x,A1,I1,I3,F7.2,F8.4,F17.9,F13.8,F12.8,2F9.3)
 9132 format(2I5,1x,A1,I1,I3,F7.2,F8.4,F10.3,F9.3,F10.3,F9.3,F10.3,
     *       4F9.3)
 9133 format(2I5,1x,A1,I1,I3,F7.2,F8.4,2F17.9,2F13.8,F12.8,2F9.3)
 9134 format(2I5,1x,A1,I1,I3,F7.2,F8.4,F10.3,F9.3,2F10.3,6F9.3)
 9140 format(2I5,1x,A1,I1,1x,F17.9,5F12.8,F9.3)
 9150 format(2I5,1x,A1,I1,I7,F17.9,5F12.8,F9.3)
 9152 format(2I5,1x,A1,I1,I7,F17.9,3F12.8,F9.3)
 9160 format(/1x,'Total interaction (TIE) to selected fragments in ',
     *           'kcal/mol',/1x,'(final is correlated if available).',/)
 9161 format(4I4,1x,A1,I1,I3,F7.2,F8.4,F17.9,F13.8,F12.8,F9.3)
 9163 format(4I4,1x,A1,I1,I3,F7.2,F8.4,2F17.9,2F13.8,F12.8,F9.3)
 9165 format(/1x,'Fragment-wise subsystem partition summary.',
     *       /1x,'Epart(I) = EBB(I) + EU,int(i)/2 where I is in i',
     *       /1x,'EBB(I) = E''''(I) + Esolv(I)+EC,int(i)/2',
     *       /4x,'I  Name     Sys        Epart            EBB         ',
     *           ' EU,int(i), kcal/mol',
     *       /1x,79(1H-))
c    *           ' Esolv    EU,int(i)',
 9166 format(I5,1x,A8,1x,i4,2F17.9,(8F9.3,/))
 9167 format(/1x,'Subsystem partition summary.',
     *       /1x,'E''(i) = sum(I in i){Epart(I)}', 
     *       /1x,'EBB(i)=E''(i)+sum(j){EC,int(j)}/2',
     *       /1x,'Epart(i)=EBB(i)+sum(j){EU,int(j)}/2 = ',
     *           'E''(i)+sum(j){Eint(j)}/2',
     *       /1x,'E=sum{Epart(i)} or ',
     *           'E=sum{EBB(i) + sum(j>i){EU,int(j)}}',
     *       /1x,'Also, E=sum{E''(i) + sum(j>i){Eint(j)}}, ',
     *       /1x,'where Eint(i)=EBB(i)-E''(i)+EU,int(i)',
     *       /3x,'Sys Charge        E''               Epart      ',
     *         '      EBB          EU,int(i), kcal/mol',
     *       /1x,79(1H-))
 9168 format(I5,F8.4,3F17.9,(8F9.3,/))
 9170 format(1x,I5,'(',A8,'): SCF',F11.4,', final TIE',F11.4)
 9175 format(1x,I5,'(',A8,'): ES',F11.4,', SCF',F11.4,', final TIE',
     *          F11.4)
 9200 format(/1x,'Total ',A1,'nergy of the molecule: E',A6,'(',I1,')=',
     *           F19.9)
 9201 format(/1x,'Total energy, BDA corrected E',A6,'(',I1,')=',F19.9)
 9202 format(/1x,'The backbone ',A1,'nergy EBB',A6,'(',I1,')=',F19.9)
 9203 format(/1x,'The total PIEDA energy E',A6,'(',I1,')=',F19.9)
 9205 format(1x,'R=',F6.2,' N=',2I6,' (C,N) ',I1,'-body E=',F14.9,F19.9)
 9206 format(1x,'R=',F6.2,' N=',I6,' (N) ',I1,'-body E=',F19.9)
 9210 format(//1x,'Total energy E',A6,'(2) without BSSE correction: ',
     *        F19.9/1x,24x,'with    BSSE correction: ',F19.9)
 9215 format(/1x,'The 0-th order, P and REP energy =',F19.9)
 9220 format(/1x,'Nuclear repulsion energy:   ',F22.9/,
     *        1x,'Total number of atoms:             ',I9/,
     *        1x,'Total number of electrons:         ',I9/,
     *        1x,'Total charge:                      ',I9/,
     *        1x,'Total spin multiplicity:           ',I9/,
     *        1x,'Total number of basis functions:   ',I9/,
     *        1x,'Total number of molecular orbitals:',I9)
 9230 format(1x,'Dipole moment D(xyz),DA(',I1,')=',4F13.7)
 9231 format(1x,'Quadrupole moment Q(XX,YY,ZZ)=',3F12.7/,
     *       1x,'                   (XY,XZ,YZ)=',3F12.7)
 9232 format(1x,'Octopole moment O(XXX,XXY,XXZ,XYY)=',4F11.6/,
     *       1x,'                 (YYY,YYZ,XZZ,YZZ)=',4F11.6/,
     *       1x,'                         (ZZZ,XYZ)=',2F11.6)
 9250 format(/1x,'Three-body FMO properties.',/1x,26(1H=)/)
 9260 format(/5x,'I    J    K DL   RMIN   RMAX       E"corr',6x,
     *      'deltaE"','IJK,corr/uncorr',2x,'dDIJK*VIJK',2x,'ddDIJK*VIJK'
     *       ,4x,'Gsol',5x,'tot',/1x,115(1H-))
 9262 format(/5x,'I    J    K DL   RMIN   RMAX deltaE''IJK dDIJK*VIJK',
     *        3x,'total     EPL0      Eex0     Ectmix0   Edisp0',
     *       /1x,102(1H-))
 9265 format(/5x,'I    J    K DL   RMIN   RMAX       E"corr',6x,
     *      'deltaE"','IJK,corr/uncorr',2x,'dDIJK*VIJK',5x,'Gsol',5x,
     *       'tot',/1x,103(1H-))
 9270 format(/5x,'I    J    K DL   RMIN   RMAX deltaE''IJK dDIJK*VIJK',
     *        3x,'total',6x,a15,3x,a6,6x,'Gsol',
     *       /1x,102(1H-))
 9275 format(1x,3I5,1x,A1,I1,2F7.2,F17.9,4F12.8,2F9.3)
 9280 format(1x,3I5,1x,A1,I1,2F7.2,F17.9,3F12.8,2F9.3)
 9285 format(1x,3I5,1x,A1,I1,2F7.2,7F10.3)
 9300 FORMAT(/1X,'               ONE ELECTRON ENERGY =',F24.10/
     *        1X,'               TWO ELECTRON ENERGY =',F24.10/
     *        1X,'          NUCLEAR REPULSION ENERGY =',F24.10/
     *       38X,23(1H-)/
     *        1X,'                      TOTAL ENERGY =',F24.10)
 9305 FORMAT(/1X,'ELECTRON-ELECTRON POTENTIAL ENERGY =',F24.10/
     *        1X,' NUCLEUS-ELECTRON POTENTIAL ENERGY =',F24.10/
     *        1X,'  NUCLEUS-NUCLEUS POTENTIAL ENERGY =',F24.10/
     *       38X,23(1H-)/
     *        1X,'            TOTAL POTENTIAL ENERGY =',F24.10/
     *        1X,'              TOTAL KINETIC ENERGY =',F24.10/
     *        1X,'                VIRIAL RATIO (V/T) =',F24.10)
 9310 format(1x,'TOTAL FOCK NUCLEI=',F25.13)
 9315 format(26I3)
 9350 format(/8x,'FMO-PCM properties.',/1x,38(1H=),
     *       /1x,'Ges is the energy of electrostatic interaction,',/1x,
     *       'included in the individual n-mer and total FMO energies.',
     *       /1x,'Internal energies are for solute only excluding all ',
     *           'interactions with solvent.',
     *       /1x,'The total energies in PCM are computed as:',
     *       /1x,'Gsol+Ecorr or Gsol+Euncorr.')
 9352 format(/1x,'WARNING: in this run, pairwise ES solvent screening',
     *       /1x,'is not separated from monomer ES solvent energies.',
     *       /1x,'The total properties are correct.')
 9353 format(/1x,'WARNING: EX for connected dimers/trimers unavailable')
 9355 format(/1x,'Total Gsol(',I1,')=',F10.3,' kcal/mol.',
     *       /1x,'Shift to convert internal to QM energy, Des(',I1,')=',
     *           F10.3,' kcal/mol.')
 9357 format(/1x,'Monomer shifts to convert to QM energies in',
     *           ' a.u. and kcal/mol, eps(eff) and omega')
 9360 format(/1x,'Monomer surface areas (in A**2), charges (a.u.) and',
     *           ' solute-solvent energies (kcal/mol).',
     *       //4x,'I',14x,'surf_cav',1x,'disp/rep',2x,'surf_es',1x,
     *            'cover,%',1x,'q_cav',4x,'eps_eff',4x,
     *            'DXsol',5x,'DYsol',5x,'DZsol',5x,'Ges',5x,a4,
     *            4x,'Gdisp',5x,'Grep',5x,'Gsol',/1x,143(1H-))
 9365 format(I5,' (',A8,') ',3F9.1,F6.1,F9.4,F10.3,3F10.4,5F9.3)
 9366 format(I5,' (',A8,') ',27x,6x,9x,8x,5F9.3,
     *       /1x,'Ges, Gdisp, Grep for monomers contain contributions ',
     *           'from dummy spheres.')
 9367 format(I5,F19.9,F9.3,F10.4,F8.3)
 9370 format(/1x,'Pair ES surface (A**2) and solute-solvent pair ',
     *           'interactions (kcal/mol):',
     *       /3x,'dGes is the direct pairwise solvent screening,',
     *       /3x,'dGCTes is the coupling of solvent es and solute CT',
     *       /3x,'dGdisp and dGrep are solvent corrections to pair',
     *            ' interactions,',
     *       /3x,'dGsol=dGes+dGCTes+dGdisp+dGrep.',
     *      //4x,'I    J',4x,'surf_es cover,%  dGes    dGCTes   dGdisp',
     *           '    dGrep    dGsol',/1x,71(1H-))
 9375 format(2I5,F10.1,F7.1,5F9.3)
 9380 FORMAT(/1X,'ELECTROSTATIC INTERACTION    =',F15.3,' KCAL/MOL',/,
     *        1X,'CDS INTERACTION              =',F15.3,' KCAL/MOL',/,
     *        1X,'TOTAL INTERACTION            =',F15.3,' KCAL/MOL')
 9390 FORMAT(/1X,'ELECTROSTATIC INTERACTION  =',F12.3,' KCAL/MOL'/
     *        1X,'PIEROTTI CAVITATION ENERGY =',F12.3,' KCAL/MOL'/
     *        1X,'DISPERSION FREE ENERGY     =',F12.3,' KCAL/MOL'/
     *        1X,'REPULSION FREE ENERGY      =',F12.3,' KCAL/MOL'/
     *        1X,'TOTAL INTERACTION          =',F12.3,' KCAL/MOL'/
     *     /1X,'The first energy printed below is the best in FMO/PCM.')
 9395 format(/1x,'    Free ',A6,' energy in solvent=',F19.9,
     *       /1x,'Internal ',A6,' energy in solvent=',F19.9)
 9400 format(/1x,'n-body ',A9,' atomic charges Q(n)',
     *       //4x,'IAT  IFG   Z',7x,'Q(1)',8x,'Q(2)',8x,'Q(3)')
 9403 format(/1x,'n-body ',A9,' atomic spin populations S(n)',
     *       //4x,'IAT  IFG   Z',7x,' S(1)',7x,' S(2)')
 9402 format(/1x,'Solvent q and n-body ',A9,' solute atomic charges ',
     *       'Q(n)',
     *        //4x,'IAT  IFG   Z',2x,'surface',1x,'cover,%',3x,'q(ASC)',
     *            7x,'Q(1)',8x,'Q(2)',8x,'Q(3)')
 9410 format(1x,2I5,F6.1,3F12.6)
 9412 format(1x,2I5,F6.1,F7.2,F7.1,4F12.6)
 9420 format(/1x,'Charge transfer for each fragment:',
     *      //2x,'IFG QFG  DeltaQ     and its contributions from JFG, ',
     *           'Q(JFG->IFG).',/1x,78(1H-))
 9425 format(I5,I3,F9.4,' =',4(I5,'->',F8.4))
 9430 format(/1x,'Total absolute monomer transf. charge   =',F12.6,
     *       /1x,'Total amount of absolute transf. charge =',F12.6)
c9420 format(1x,'Donor+acceptor ',A6,' energy sum is ',F19.9)
 9440 format(/1x,'Spin transfer for open-shell fragments:',
     *       /1x,'   IFG   JFG  S(I->J)')
 9445 format(1x,2I6,F9.4)
 9450 format(/1x,'DFT exchange+correlation energy=',F19.9,/,
     *       1x,'Total electron number          =',F19.9)
 9460 format(/1x,'Destabilisation polarisation energy:  ',F12.3,
     *       /1x,'Stabilisation polarisation energy:    ',F12.3,
     *       /1x,'Total polarisation energy (kcal/mol): ',F12.3,/)
 9470 format(/1x,'n-body local atomic electrostatic hardness Eta(n)',
     *       //4x,'IAT  IFG   Z',5x,'Eta(1)',6x,'Eta(2)',6x,'Eta(3)')
 9475 format(1x,2I5,F6.1,3F12.6)
 9480 format(1x,A9,' populations for',i4,'*',I7,' atoms')
 9481 format(1x,7F11.7)
 9500 format(/1x,I6,' of ',I1,'-body energies diverged!!!',
     *       ' The total properties are not meaningful.')
 9510 format(/1x,'Because separated dimers are omitted for speed,' 
     *       /1x,'the total properties are inaccurate.',
     *       /1x,'However, we only need PCM charges at this point,',
     *       /1x,'which are unaffected by omitted ES dimers.')
 9520 format(/1x,'Atomic charge corrections from relaxed density.',
     *       /4x,'IAT  IFG   Z',6x,'Q(rel)',5x,'Q(',I1,'+rel)')
 9530 format(1x,2I5,F6.1,2F12.6)
 9540 format(/1x,'Multipole corrections from relaxed density.',
     *       /5x,'IFG',5x,'dDX',7x,'dDY',7x,'dDZ',7x,'DXr',7x,'DYr',7x,
     *           'DZr')
 9550 format(1x,I7,6F10.5)
 9560 format(/1x,'Relaxed dipole moment D=',4F13.7)
 9600 format(/1x,'PIEDA summary',
     *       /5x,'Main relations:'
     *       /5x,'EIB0 = EPL0 + EPL0DI + EES0 + EEX0 + ECTmix0 + EDI0',
     *       /5x,'EIB  = EIB0 + EPL*EX + EPL*CTmix + EPL*DI',
     *       /5x,'EIB  = Eint + EPL0d + EPL0DI + EPLd*PLd + EPLd*DI',
     *       /5x,'EIB  = E(FMO2) - E(FMO0) - E(BDA)', 
     *       /5x,'EIB0 - interaction relative to free state using free',
     *           ' state densities,',
     *       /5x,'EIB  - interaction relative to free state using PL ',
     *           'state densities,',/,
     *       /1x,'Polarisation (destabilisation)  EPL0d',F14.3,
     *       /1x,'Polarisation (stabilisation)    EPL0s',F14.3,
     *       /5x,'Polarisation (total, free state)    EPL0   ',F14.3,
     *       /5x,'Pl.-induced ',a11,'(free state) EPL0',a2,' ',F14.3,
     *       /5x,'Electrostatic (free state)          EES0   ',F14.3,
     *       /5x,'Exchange (free state)               EEX0   ',F14.3,
     *       /5x,'Charge transfer+all mixing(free s.) ECTmix0',F14.3,
     *       /5x,a11,' (free state)            E',a2,'0   ',F14.3,
     *       /9x,'Internal binding (free state)        EIB0     ',F14.3,
     *       /9x,'(Uncorr. int. bind.,free s.) EIB0-EPL0DI-EDI0 ',F14.3,
     *           ')',/,
     *       /9x,'Polarisation-exchange                EPL*EX   ',F14.3,
     *       /9x,'Polarisation-charge transfer         EPL*CTmix',F14.3,
     *       /5x,'Polarisation-',a11,' (destab)   EPLd*',a2,F14.3,
     *       /5x,'Polarisation-',a11,' (stab)     EPLs*',a2,F14.3,
     *       /9x,'Polarisation-',a11,' (total)     EPL*',a2,'   ',F14.3,
     *      /13x,'Internal binding                     EIB      ',F14.3,
     *      /13x,'(Uncorr. int. bind.) EIB-EPL0DI-EDI0-EPL*DI   ',F14.3,
     *           ')',/)
 9605 format(/1x,'Polarisation-polarisation (destab)  EPLd*PLd',F14.3,
     *       /1x,'Polarisation-polarisation (stab)    EPLs*PLs',F14.3,
     *       /1x,'Polarisation-polarisation (total)    EPL*PL ',F14.3,
     *       /5x,'Polarisation (total)        EPL=EPL0+EPL*PL ',F14.3,
     *       /1x,'(EPL*PL up to two-body terms is included in ECTmix0,',
     *       /1x,'the rest of EPL*PL is in EPL*CTmix. E(CT+mix) does',
     *           ' not include any of EPL*PL.)',
     *       /1x,'Estimate of free state CT without PL*PL, E(CT0+mix0)',
     *           F14.3,
     *       /1x,'Total monomer destabilisation,   E(PLd+corr)',F14.3)
 9610 format(/1x,'Interaction energy relative to PL state:'
     *       /1x,'Eint = EES + EEX + E(CT+mix) + EDI + E',A2,' + ESOLV',
     *       /1x,'Eint = E(FMO2) - E(FMO1) - E(BDA)',
     *       /1x,'Eint - interaction relative to PL state using PL ',
     *           'state densities.',/,
     *       /5x,'Electrostatic (PL state, incl. EPLs)  EES  ',F14.3,
     *       /5x,'Exchange (PL state)                   EEX  ',F14.3,
     *       /5x,'Charge transfer (PL state)        E(CT+mix)',F14.3,
     *       /5x,'Dispersion (PL state)                 EDI  ',F14.3,
     *       /5x,A21,' (PL state)      E',A2,'  ',F14.3,
     *       /5x,'Solvent screening (PL state)          ESOLV',F14.3,
     *       /9x,'Total interaction (PL state)         Eint  ',F14.3,/)
 9611 format(/1x,'Interaction energy relative to PL state:'
     *       /1x,'Eint = EES + E0 + E(CT*ES) + EDI + ESOLV',
     *       /1x,'Eint = E(FMO2) - E(FMO1) - E(BDA)',
     *       /1x,'Eint - interaction relative to PL state using PL ',
     *           'state densities.',/,
     *       /5x,'Electrostatic (PL state, incl. EPLs)  EES  ',F14.3,
     *       /5x,'0-th order DFTB energy (PL state)      E0  ',F14.3,
     *       /5x,'Charge transfer coupled to ES (PL) E(CT*ES)',F14.3,
     *       /5x,'Dispersion (PL state)                 EDI  ',F14.3,
     *       /5x,'Solvent screening (PL state)          ESOLV',F14.3,
     *       /9x,'Total interaction (PL state)         Eint  ',F14.3,/)
 9620 format(/1x,'Unconnected interaction energy relative to PL state:',
     *       /1x,'E''int = E''ES + E''EX + E''(CT+mix) + E''DI',
     *       /1x,'E''int = E''(FMO2) - E''(FMO1)',
     *       /1x,'E''int - interaction relative to PL state using PL ',
     *           'state densities, ',
     *       /1x,'excluding contributions from dimers between which ',
     *       /1x,'a covalent bond is fractioned.',/,
     *       /5x,'Electrostatic (PL state, incl. EPLs)   E''ES',F14.3,
     *       /5x,'Exchange (PL state)                    E''EX',F14.3,
     *       /5x,'Charge transfer (PL state)       E''(CT+mix)',F14.3,
     *       /5x,a11,' (PL state)                 E''',a2,F14.3,
     *       /5x,'Solvent screening (PL state)         E''SOLV',F14.3,
     *       /9x,'Total interaction (PL state)         E''int ',F14.3,/,
     *       /1x,'BDA energies',
     *       /5x,'Electrostatic (BDA)              EESBDA',F16.3,
     *       /5x,'Exchange (BDA)                   EEXBDA',F16.3,
     *       /5x,'Charge transfer (BDA)      E(CT+mix)BDA',F16.3,
     *       /5x,a11,' (BDA)                E',a2,'BDA',F16.3,
     *       /5x,'Solvent screening (BDA)        ESOLVBDA',F16.3,
     *       /9x,'Total interaction (BDA)         EintBDA',F16.3,/)
 9700 format(/1x,'FD properties: N=',I4,' E_AF=',F19.9,' R_AF=',F7.2,
     *           ' IJ_AF=',2I6)
 9720 format(/1x,'Other ES dimer IJ energies (kcal/mol) summed over J',
     *           ' for each I.')
 9730 format(1x,'I=',I5,' Esd=',F12.3)
 9740 format(1x,'The sum of these ES dimer energies is',F9.3,
     *          ' kcal/mol or',F20.10,' (a.u.)')
 9750 format(/1x,'FMO total energy = ',F19.9, 
     *       /1x,'CNS Exref energy = ',F19.9, 
     *       /1x,'New total energy = ',F19.9)
 9755 format(/1x,'CNS Exref energy gradient is printed next.') 
 9780 FORMAT(1x,'Spin Sz(',I1,') = ',F8.3,' actual S^2(',I1,') = ',F8.3,
     *          ' should be',F8.3)
 9800 format(/1x,'In this reduced memory run',
     *       /1x,'Euncorr, Edisp, and EBBuncorr are not calculated.',
     *    /1x,'Subtract 128 from NPRINT in $FMOPRP to obtain them.',/1x)
 9810 format(/1x,'The best FMO energy is ',F25.9)
 9900 format(/1x,'Done with FMO properties.')
      END
c
C*MODULE fmoio   *DECK vdwrout
      SUBROUTINE vdwrout(fmozan,vdwrad)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      character*8 c8dum
      character*4 ATOMNM
      PARAMETER (MAXNZ=137,units=0.52917724924D+00)
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      Common /fmoinf/ nfg,nlayer,natfmo,nbdfg,naotyp,nbody,nsegm
      dimension fmozan(*),vdwrad(MAXNZ),nout(MAXNZ)
      data ATOMNM/'    '/,c8dum/' '/
c
c     Output van der Waals radii.
c
      call viclr(nout,1,MAXNZ)
c
      write(iw,9000)
      do iat=1,natfmo
        iz=int(fmozan(iat)+0.1D+00)
        if(nout(iz).eq.0) then
          call zsymnum(c8dum,ATOMNM,iz)
          write(iw,9100) ATOMNM,vdwrad(iz)*units
          nout(iz)=1
        endif
      enddo
      return
 9000 format(4x,'Van der Waals radii in Angstrom are:')
 9100 format(7x,'R(',A3,')=',F8.4) 
      END
c
C*MODULE fmoio   *DECK eminout
      subroutine eminout(imode,ilay,ifgfmo0,layfrg,emon,outpune,some)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      character*8 gprnam,end
      LOGICAL GOPARR,DSKWRK,MASWRK,outpune,some
c     parameter (sawdust=1.0D-25)
      dimension layfrg(*),emon(nfg,4,nlayer)
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      Common /fmoinf/ nfg,nlayer,natfmo,nbdfg,naotyp,nbody,nsegm
      common /fmoopt/ espsca(9),RESPAP(2),rESPPC(2),rESDIM,restri(4),
     *                rcorsd,respct,convfg,cnvdmp,coroff,rflmo(4),
     *                orshft,orshft2,cnvafo,ascreen(4),IXESP,mxitfg,
     *                nguess,NBSSE,modorb,modpar,irststp,irstlay,nprfmo,
     *                nfmopal,modprp,maxl1c,ipieda,modgrd,modesp,ivmul,
     *                modlmo,nopden,mofock,modfd,modfmm,ncentm,ndualb,
     *                ngab,modpan
c
      dimension ewrk(4)
c
      data gprnam/' $FMOENM'/,end/' $END'/
c
c     read (imode=0) or write (imode=1) monomer energies
c     ilay is ignored for imode=0.
c
c     nen=1
c     if(nbsse.eq.3) nen=nen+1
c     if(gcorrel) nen=nen+nen
      nen=4
      if(imode.eq.0) then
        CALL SEQREW(IR)
        CALL FNDGRP(IR,gprnam,IEOF)
        if(ieof.ne.0) then
          if(maswrk) write(iw,9020) gprnam
c         abort if not doing monomer SCF
          if(irststp.gt.2.or.irststp.eq.2.and.nbsse.eq.3) call abrt
c         for BSSE one needs free monomer energies
          return
        endif
        if(maswrk) then
  100   continue
          READ(IR,9010,END=120,ERR=120) ifg,ilay,(ewrk(ien),ien=1,nen)
          if(ifg.gt.nfg.or.ilay.gt.nlayer) then
            if(maswrk) write(6,*) 'bad indices read in',gprnam,ifg,ilay
            call abrt
          endif
c         this is to avoid double counting; assign all energies to the
c         grand master. Avoid parallel confusion by assigning some sawdust
c         to other masters.
c         At present no need for sawdusting monomers.
          if(some) write(iw,9000) ifg,ilay,(ewrk(ien),ien=1,nen)
c         if(mygroup.ne.0) call dacopy(nen,sawdust,ewrk,1)
          call dcopy(nen,ewrk,1,emon(ifg,1,ilay),nfg)
c         write(iw,9000) ifg,ilay,(emon(ifg,ien,ilay),ien=1,nen)
        goto 100
  120   continue 
        endif
        IF (GOPARR) CALL DDI_BCAST(2421,'F',emon,nfg*4*nlayer,MASTER)
      else
        if(outpune) then
          write(ip,9100) gprnam 
          if(ifgfmo0.ne.0) then
            if(layfrg(1).ge.ilay) write(ip,9010) ifgfmo0,ilay,
     *                                   (emon(1,ien,ilay),ien=1,nen)
          else
            do ifg=1,nfg
              if(layfrg(ifg).ge.ilay) write(ip,9010) ifg,ilay,
     *                                   (emon(ifg,ien,ilay),ien=1,nen)
            enddo
          endif
          write(ip,9100) END
        endif
      endif
      return
 9000 format(1x,'E ifg=',I5,'ilay=',I2,' E=',4F14.10)
 9010 format(1x,I5,I2,4F18.10)
 9020 format(1x,'group name=',a,' was not found')
 9100 format(a8)
      end
c
C*MODULE fmoio   *DECK edin
      subroutine edin(nfg2,nen,layfrg,edim,some)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      character*8 gprnam
      parameter (zero=0.0D+00,sawdust=1.0D-25)
      dimension layfrg(*),edim(nfg2,*)
      LOGICAL GOPARR,DSKWRK,MASWRK,isgddi,parout,INITGDDI,some,wasgddi,
     *        MLGDDI
      COMMON /GDDI/   ISCOPE,NGROUPS,MYGROUP,MEGLOB,npglob,nnglob,JBTYP,
     *                ISGDDI,PAROUT,INITGDDI,wasgddi,MLGDDI,NSUBGR,
     *                MeUniv,NPUniv,numdlb,myworld,nworlds,mogddi
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      Common /fmoinf/ nfg,nlayer,natfmo,nbdfg,naotyp,nbody,nsegm
      common /fmoopt/ espsca(9),RESPAP(2),rESPPC(2),rESDIM,restri(4),
     *                rcorsd,respct,convfg,cnvdmp,coroff,rflmo(4),
     *                orshft,orshft2,cnvafo,ascreen(4),IXESP,mxitfg,
     *                nguess,NBSSE,modorb,modpar,irststp,irstlay,nprfmo,
     *                nfmopal,modprp,maxl1c,ipieda,modgrd,modesp,ivmul,
     *                modlmo,nopden,mofock,modfd,modfmm,ncentm,ndualb,
     *                ngab,modpan
      dimension ewrk(9)
      data gprnam/' $FMOEND'/
c
c     read dimer energies
c     this does not work with dolat (change the leading dim in edim).
c
c     nen=nen0
c     if(nbsse.ne.0.and.IRSTSTP.eq.5) nen=nen-2
c     if(imode.eq.0) then
      CALL SEQREW(IR)
      CALL FNDGRP(IR,gprnam,IEOF)
      if(ieof.ne.0) then
        write(iw,9010) gprnam
 9010 format(1x,'group name',a,' not found: all dimers will be run.')
c       abort if not doing dimer SCF
        if(irststp.gt.4) then
          write(6,*) 'bad irststp',irststp
          call abrt
        endif
        return
      endif
      if(maswrk) then
      nread=0
  100 continue
c       if(nen.le.4) then
        READ(IR,9210,END=120,ERR=120) ifg,jfg,ilay,(ewrk(ien),ien=1,nen)
c       else
c       READ(IR,9215,END=120,ERR=120) ifg,jfg,ilay,(ewrk(ien),ien=1,nen)
c       endif
        if(ifg.le.1.or.ifg.gt.nfg.or.jfg.gt.nfg.or.jfg.gt.ifg.or.
     *     ilay.gt.nlayer) then
          write(6,*) 'bad indices read in',gprnam,ifg,jfg,ilay
          call abrt
        endif
        nread=nread+1
c       ignore energies for lower levels
        if(ilay.lt.layfrg(ifg).and.ilay.lt.layfrg(jfg)) goto 100
        if(some) then
          if(nen.le.3) then 
            write(iw,9000) ifg,jfg,ilay,ewrk(1)
          else
            write(iw,9000) ifg,jfg,ilay,ewrk(1),(ewrk(ii),ii=3,4)
          endif
        endif
c       this is to avoid double counting and prevent recompututing
c       if(meglob.ne.0) call dacopy(nen,sawdust,ewrk,1)
        if(mygroup.ne.0) call dacopy(nen,sawdust,ewrk,1)
c       ilay is ignored in saving data. This is because we expect the 
c       layers be given in ascending order, so that higher layers if exist
c       should replace the data for the lower ones. No check is done.
c       In fact ilay is punched for the user to check the data.
        if(ifg.lt.jfg) then
          write(6,*) 'bad index pair',ifg,jfg
          call abrt
        endif 
        ijfg=((ifg-1)*(ifg-2))/2+jfg
        call dcopy(nen,ewrk,1,edim(ijfg,1),nfg2)
        if(nbsse.ne.0.and.IRSTSTP.eq.5) then
          edim(ijfg,4)=zero
          edim(ijfg,5)=zero
cnb       correlation (6)? 
        endif
      goto 100
  120 continue 
      write(iw,9100) nread,nfg2
      endif
      IF (GOPARR) then
        CALL DDI_BCAST(2422,'F',edim,nfg2*nen,MASTER)
c       Slaves in group zero should have energies replaced by sawdust to avoid
c       double counting.
c       This code is not correct for BSSE restarts.
        if(nbsse.ne.0.and.IRSTSTP.eq.5) call abrt 
c
        if(mygroup.eq.0.and..not.maswrk) then
          do ijfg=1,nfg2
            if(edim(ijfg,1).ne.zero) 
     *        call dcopy(nen,sawdust,1,edim(ijfg,1),nfg2)
          enddo
        endif
      endif
      return
 9000 format(1x,'E ifg=',I5,' jfg=',I5,' ilay=',I2,' E=',3F14.7)
 9100 format(1x,i10,' out of',I10,' dimer energies read in.')
 9210 format(1x,2I5,I3,3F22.10,2(/14x,3F22.10))
c9215 format(1x,2I5,I3,2F22.10,/14x,2F22.10)
      end
c
C*MODULE fmoio   *DECK etin
      subroutine etin(nfg3,nent,layfrg,etrim,some)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      character*8 gprnam
      parameter (zero=0.0D+00,sawdust=1.0D-25)
      dimension layfrg(*),etrim(nfg3,*)
      LOGICAL GOPARR,DSKWRK,MASWRK,isgddi,parout,INITGDDI,some,wasgddi,
     *        MLGDDI
      COMMON /GDDI/   ISCOPE,NGROUPS,MYGROUP,MEGLOB,npglob,nnglob,JBTYP,
     *                ISGDDI,PAROUT,INITGDDI,wasgddi,MLGDDI,NSUBGR,
     *                MeUniv,NPUniv,numdlb,myworld,nworlds,mogddi
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      Common /fmoinf/ nfg,nlayer,natfmo,nbdfg,naotyp,nbody,nsegm
      common /fmoopt/ espsca(9),RESPAP(2),rESPPC(2),rESDIM,restri(4),
     *                rcorsd,respct,convfg,cnvdmp,coroff,rflmo(4),
     *                orshft,orshft2,cnvafo,ascreen(4),IXESP,mxitfg,
     *                nguess,NBSSE,modorb,modpar,irststp,irstlay,nprfmo,
     *                nfmopal,modprp,maxl1c,ipieda,modgrd,modesp,ivmul,
     *                modlmo,nopden,mofock,modfd,modfmm,ncentm,ndualb,
     *                ngab,modpan
      dimension ewrk(9)
      data gprnam/' $FMOENT'/
c
c     read trimer energies
c
      CALL SEQREW(IR)
      CALL FNDGRP(IR,gprnam,IEOF)
      if(ieof.ne.0) then
        write(iw,9010) gprnam
c       abort if not doing dimer SCF
        if(irststp.gt.9) then
          write(6,*) 'bad irststp',irststp
          call abrt
        endif
        return
      endif
      if(maswrk) then
      nread=0
  100 continue
        READ(IR,9210,END=120,ERR=120) ifg,jfg,kfg,ilay,
     *                                (ewrk(ien),ien=1,nent)
        if(ifg.le.1.or.ifg.gt.nfg.or.jfg.gt.ifg.or.kfg.gt.jfg.or.
     *     ilay.gt.nlayer) then
          write(6,*) 'bad indices read in',gprnam,ifg,jfg,kfg,ilay
          call abrt
        endif
        nread=nread+1
c       ignore energies for lower levels
        if(ilay.lt.layfrg(ifg).and.ilay.lt.layfrg(jfg)) goto 100
        if(some) then
          if(nent.le.3) then 
            write(iw,9000) ifg,jfg,kfg,ilay,ewrk(1)
          else
            write(iw,9000) ifg,jfg,kfg,ilay,ewrk(1),(ewrk(ii),ii=3,4)
          endif
        endif
c       this is to avoid double counting and prevent recompututing
        if(mygroup.ne.0) call dacopy(nent,sawdust,ewrk,1)
c       ilay is ignored in saving data. This is because we expect the 
c       layers be given in ascending order, so that higher layers if exist
c       should replace the data for the lower ones. No check is done.
c       In fact ilay is punched for the user to check the data.
        ijkfg=((ifg-1)*(ifg-2)*(ifg-3))/6+((jfg-1)*(jfg-2))/2+kfg
        call dcopy(nent,ewrk,1,etrim(ijkfg,1),nfg3)
      goto 100
  120 continue 
      write(iw,9100) nread,nfg3
      endif
      IF (GOPARR) then
        CALL DDI_BCAST(2422,'F',etrim,nfg3*nent,MASTER)
c       Slaves in group zero should have energies replaced by sawdust to avoid
c       double counting.
c
        if(mygroup.eq.0.and..not.maswrk) then
          do ijkfg=1,nfg3
            if(etrim(ijkfg,1).ne.zero) 
     *        call dcopy(nent,sawdust,1,etrim(ijkfg,1),nfg3)
          enddo
        endif
      endif
      return
 9000 format(1x,'E ifg=',I5,' jfg=',I5,' kfg=',I5,' ilay=',I2,' E=',
     *       3F14.7)
 9100 format(1x,i10,' out of',I10,' trimer energies read in.')
 9010 format(1x,'group name',a,' not found: all trimers will be run.')
 9210 format(1x,3I5,I3,2F22.10,2F16.10)
      end
c
C*MODULE fmoio   *DECK PRTRILe
      SUBROUTINE PRTRILe(D,N)
      use mx_limits, only: mxatm,mxao
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      LOGICAL GOPARR,DSKWRK,MASWRK
      DIMENSION D(*)
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON /OUTPUT/ NPRINT,ITOL,ICUT,NORMF,NORMP,NOPK
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /RUNLAB/ TITLE(10),ANAM(MXATM),BNAM(MXATM),BFLAB(MXAO)
C
C     ----- PRINT SYMMETRIC MATRIX -D- OF DIMENSION -N- -----
C     THE ROWS WILL BE LABELED WITH BASIS FUNCTION TAGS.
C
      IF (MASWRK) THEN
      MAX = 5
      IF (NPRINT .EQ. 6) MAX = 10
      MM1 = MAX - 1
      DO 120 I0=1,N,MAX
         IL = MIN(N,I0+MM1)
         WRITE(IW,9008)
         WRITE(IW,9028) (I,I=I0,IL)
         WRITE(IW,9008)
         IL = -1
         DO 100 I=I0,N
            IL=IL+1
            J0=I0+(I*I-I)/2
            JL=J0+MIN(IL,MM1)
            WRITE(IW,9048) I,BFLAB(I),(D(J),J=J0,JL)
  100    CONTINUE
  120 CONTINUE
      END IF
      RETURN
 9008 FORMAT(1X)
 9028 FORMAT(15X,10(4X,I4,3X))
c9028 FORMAT(15X,10(3X,I8,3X))
 9048 FORMAT(I5,2X,A8,10E15.6)
c9048 FORMAT(I5,2X,A8,10F12.8)
      END
C*MODULE fmoio   *DECK fmoccp
      SUBROUTINE fmoccp(ibody,escf,extracc)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*13 CMET(9)
      LOGICAL CNVR,QDR,QDCR
      DOUBLE PRECISION CCERI
      COMMON /CCENGY/ ENRG,EREF,EMP2,eccn,ETOT(6),ECORR(6),
     *                DIAGS(3),AMPMX(5,2),IAMPMX(5,4,2),XO1,XO2,
     *                DIFMAX,DIFFENG,ITER,CNVR
      COMMON /CCINFO/ TSH,NH,NP,MET,MEM,ICONV,MAXIT,IREST,IDISC
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON/T4TOT/ETT2,ETS2,EQ1A,EQ1B,EQ2A,EQ2B,EQ3A,EQ3B,EQ4A,
     *EQ4B,O4A,O4B,ECORQ(8),ETOTQ(8)
      COMMON/QUADRUPLE/QDR,QDCR
      COMMON /WFNOPT/ SCFTYP,VBTYP,DFTYPE,TDDFTYP,CITYP,CCTYP,
     *                MPLEVL,MPCTYP
      COMMON /ERIFLG/ CCERI
      DATA CMET/'LCCD,','CCD,','CCSD,','CCSD[T],','CCSD(T),',
     *          'R-CCSD[T],','R-CCSD(T),','CR-CCSD[T],','CR-CCSD(T),'/
      data CRCCL/8HCR-CCL  /
      DATA STANDARD,RI,CD/8HSTANDARD,8HRI      ,8HCD      /
      dimension extracc(*)
C
c     This is a mutilated clone of DPRINT.
c     Note that it is assumed that there is no 2nd order correlated density
c     contribution of the form Tr(deltaDcorr*V) to n-mers, which is normally
c     so unless some new method such as FMO1-RCC (relaxed CC) is developed.
c
      IMET=MET+1
      IF(MET.GT.2) IMET=3
c
      EREF=escf
      EMP2=extracc(1)
      eccn=extracc(2)
      call dcopy(6,extracc(3),1,ECORR,1)
      if(qdr) call dcopy(7,extracc(9),1,ECORQ(2),1)
c
      write(iw,9000)
      IF(CCERI .EQ. STANDARD)THEN
         WRITE(iw,1000) ibody,'MP2,         ',EMP2,eref+EMP2
         WRITE(iw,1000) ibody,CMET(IMET),eccn,eref+ECCN
      ELSE IF(CCERI .EQ. RI)THEN
         WRITE(iw,1000) ibody,'RI-MP2       ',EMP2,eref+EMP2
         WRITE(iw,1000) ibody,'RI-CCSD      ',eccn,eref+ECCN
      END IF
      IF(MET.GT.2) THEN
      IF(CCERI .EQ. STANDARD)THEN
        WRITE(iw,1000) ibody,'CCSD[T],     ',ECORR(1),eref+ECORR(1)
        WRITE(iw,1000) ibody,'CCSD(T),     ',ECORR(2),eref+ECORR(2)
      ELSE IF(CCERI .EQ. RI)THEN
        WRITE(iw,1000) ibody,'RI-CCSD(T)   ',ECORR(2),eref+ECORR(2)
      END IF
        IF(MET.GT.3) THEN
          WRITE(iw,1000) ibody,'R-CCSD[T],   ',ECORR(3),eref+ECORR(3)
          WRITE(iw,1000) ibody,'R-CCSD(T),   ',ECORR(4),eref+ECORR(4)
          IF(MET.GT.4.AND..NOT.(QDR.AND..NOT.QDCR)) THEN
            WRITE(iw,1000) ibody,'CR-CCSD[T],  ',ECORR(5),eref+ECORR(5)
            WRITE(iw,1000) ibody,'CR-CCSD(T),  ',ECORR(6),eref+ECORR(6)
          END IF
        END IF
      END IF
c
      IMET=MET+1
      IF(MET.GT.2) IMET=MET+2
      IF(QDR) THEN
        WRITE(iw,1000) ibody,'CCSD(TQ),B   ',ECORQ(2),eref+ECORQ(2)
        WRITE(iw,1000) ibody,'R1-CCSD(TQ),A',ECORQ(3),eref+ECORQ(3)
        WRITE(iw,1000) ibody,'R1-CCSD(TQ),B',ECORQ(4),eref+ECORQ(4)
        WRITE(iw,1000) ibody,'R2-CCSD(TQ),A',ECORQ(5),eref+ECORQ(5)
        WRITE(iw,1000) ibody,'R2-CCSD(TQ),B',ECORQ(6),eref+ECORQ(6)
        IF(QDCR) THEN
          WRITE(iw,1000) ibody,'CR-CCSD(TQ),A',ECORQ(7),eref+ECORQ(7)
          WRITE(iw,1000) ibody,'CR-CCSD(TQ),B',ECORQ(8),eref+ECORQ(8)
        END IF
      END IF
      if(CCTYP.EQ.CRCCL) then
        e1=eccn+extracc(16)
        e2=eccn+extracc(17)
        WRITE(iw,1000) ibody,'CCSD(2)_T,   ',e1,e1+eref
        WRITE(iw,1000) ibody,'CR-CCSD(T)_L,',e2,e2+eref
        WRITE(iw,2000)
      END IF
c
      write(iw,9000)
      return
C
 1000 FORMAT(1X,'FMO',I1,'-',A13,' corr. energy=',F15.9,
     *       ' total E=',F19.9)
 2000 FORMAT(1X,'E corr printed above is computed for CR-CCSD(T)_L')
 9000 FORMAT(1X)
      END
c
C*MODULE fmoio   *DECK fmompp
      SUBROUTINE fmompp(ibody,escf,ecorr)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON /MP2PAR/ OSPT,CODEMP,SCSPT,TOL,METHOD,NWDMP2,MEMPRI,MPPROP,
     *                NACORE,NBCORE,NOA,NOB,NO,NBF,NOMIT,MOCPHF,MAXITC
      DATA ANONE/8HNONE    /
C
c     Alternative MP2 properties are only set in CODE=IMS/DDI, RUNTYP=ENERGY
c     and CODE=SERIAL; both for RHF reference only.
c     Some efforts are invested in making ecorr exactly 0 if undefined.
c
      if(ecorr.ne.0) then
        write(iw,9000)
        IF(SCSPT.NE.ANONE) then
          WRITE(iw,1000) ibody,'MP2     ',ecorr,escf+ecorr
        else
          WRITE(iw,1000) ibody,'MP2(SCS)',ecorr,escf+ecorr
        endif
        write(iw,9000)
      endif
      return
C
 1000 FORMAT(1X,'FMO',I1,'-',A8,' corr. energy=',F15.9,
     *       ' total E=',F19.9)
 9000 FORMAT(1X)
      END
c
C*MODULE fmoio   *DECK fmoh3c
C
C     @brief Print HF-3c results.
C
C     @details Write to output the results.
C
C     @author Dmitri Fedorov
C
      SUBROUTINE fmoh3c(ibody,extracc)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      dimension extracc(*)
      parameter (AUTOKCAL=627.509541D+00)
C
c     Additional HF3c properties.
c
      WRITE(iw,1000) ibody,'DISP',extracc(1)*AUTOKCAL,extracc(1)
      WRITE(iw,1000) ibody,'GCP ',extracc(2)*AUTOKCAL,extracc(2)
      WRITE(iw,1000) ibody,'SRB ',extracc(3)*AUTOKCAL,extracc(3)
      return
C
 1000 FORMAT(1X,'FMO',I1,A8,' energy=',F20.6,' kcal/mol',F16.9,' a.u.')
      END
c
C*MODULE fmoio   *DECK fmogout
      SUBROUTINE fmogout(ibody,masout,indat,fmozan,fmode,DEFTF,TORQF)
      USE mx_limits, only: mxfrg
      USE comm_FGRAD
      USE comm_FRGINF
      USE comm_EFPFMO
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      logical masout
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      parameter (one=1.0D+00)
      Common /fmoinf/ nfg,nlayer,natfmo,nbdfg,naotyp,nbody,nsegm
      dimension indat(*),fmozan(*),fmode(3,natfmo,*),DEFTF(3,nfrg,*),
     *          TORQF(3,nfrg,*),DEFTA(3),dnam(3)
      DATA DNAM /4HE'X ,4HE'Y ,4HE'Z /
c
      write(iw,9300) ibody
c     CALL EGOUT(fmode,NATfmo)
c     EGOUT cannot be used due to ANAM and BNAM being undefined.
      if(masout) then
        WRITE(IW,9000) (DNAM(J),J=1,3)
        DO I=1,NATfmo
          WRITE(IW,9010) I,indat(i),fmozan(i),(fmode(J,I,ibody),J=1,3)
        enddo
        if(IEFPFMO.NE.0) then
          do i=1,3
            DEFTA(i)=ddot(nfrg,DEFTF(i,1,ibody),3,one,0)
          enddo
          WRITE (IW,8000) ibody
          CALL EFOUT(DEFTF(1,1,ibody),DEFTA,TORQF(1,1,ibody))
        endif
      endif
      if(IEFPFMO.NE.0) then
        call EGMAX(fmode(1,1,ibody),natfmo*3,DEFTF(1,1,ibody),
     *             TORQF(1,1,ibody),NFRG,FMAX,FRMS)
      else
        call EGMAX(fmode(1,1,ibody),natfmo*3,DEFT,TORQ,NFRG,FMAX,FRMS)
      endif
      WRITE (IW,9030) ibody,FMAX,FRMS
c
      RETURN
 8000 FORMAT(/21X,29('-')/
     *        21X,'FRAGMENT GRADIENT INFORMATION, FMO',I1,/
     *        21X,29('-')/)
 9000 FORMAT(/3X,' ATOM#   FRG#   Z',12X,A4,13X,A4,13X,A4)
 9010 FORMAT(1X,I8,I6,F6.1,3F18.9)
 9030 FORMAT(/6X,'(',I1,') ',19HMAXIMUM GRADIENT = ,F10.7,4X,
     *           'RMS GRADIENT =',F10.7)
 9300 format(/1x,'Energy gradient (hartree/bohr), no BSSE: G(',I1,')')
      END
C*MODULE fmoio   *DECK closefrg
      SUBROUTINE closefrg 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      LOGICAL GOPARR,DSKWRK,MASWRK
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON /OUTPUT/ NPRINT,ITOL,ICUT,NORMF,NORMP,NOPK
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
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
c
      if(iand(nfmopal,2).ne.0.and.maswrk.and.nprint.ne.-5)
     *  write(iw,8000) icurlay,icurfg,jcurfg,kcurfg,icurit
 8000 format(/1x,16(1H~),'  END  FMO OUTPUT for',I2,3I6,I4,1x,16(1H~),/)
      RETURN
      END
C*MODULE fmoio   *DECK printeig
      SUBROUTINE printeig(nefmo,ichfmo,mulfmo,numfrg,eigfmo,ibfmo,dodos,
     *                    m1efmo,nspins)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      parameter (toeV=27.21138386D+00)
      CHARACTER*8 sspin(2)
      logical dodos
      dimension numfrg(*),eigfmo(m1efmo,*),ibfmo(*),NABOCC(2)
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      Common /fmoinf/ nfg,nlayer,natfmo,nbdfg,naotyp,nbody,nsegm
      data sspin/'alpha','beta'/
c
      NABOCC(1)=(NEfmo-ICHfmo+MULfmo-1)/2
      NABOCC(2)=(NEfmo-ICHfmo-MULfmo+1)/2
      do ispin=1,nspins
c     compress (take away padded zeros)
      ind=0
      ind0=0
      do i=1,nfg
        l1=iand(numfrg(i),65535)
        nzero=0
        do j=l1,1,-1
          if(eigfmo(ind+j,ispin).ne.0) goto 100
          nzero=nzero+1
        enddo
  100   continue
c       write(iw,*) 'Fragment=',i,', extras=',nzero
        do j=1,l1-nzero
          eigfmo(ind0+j,ispin)=eigfmo(ind+j,ispin)
        enddo
        ind=ind+l1
        ind0=ind0+l1-nzero
      enddo
c     ind0 can be zero if FMO1 has not converged.
      if(ind0.gt.0) then
c       nafmo=nefmo/2
        nafmo=NABOCC(ispin)
        nhomo=nafmo
        nlumo=nafmo+1
c       Assume RHF
        write(iw,9000) ind0
        call DSORT(ind0,eigfmo(1,ispin),ibfmo)
        write(ip,9200) ind0,sspin(ispin),nhomo,nlumo,
     *                 (eigfmo(nlumo,ispin)-eigfmo(nhomo,ispin))*toeV
        write(ip,9210) (eigfmo(i,ispin),i=1,ind0)
        if(dodos) then
          ehomo=eigfmo(nhomo,ispin)*toeV
          elumo=eigfmo(nlumo,ispin)*toeV
          efermi=(ehomo+elumo)/2
          write(iw,9220) sspin(ispin),ehomo,elumo,efermi,elumo-ehomo
        endif
      endif
      enddo
      RETURN
 9000 Format(/1x,'Sorting',I10,' monomer orbital energies...')
 9200 Format(1x,I8,' monomer ',A5,' orbital energies in FMO1 are:',
     *       /1x,'HOMO=',I8,' , LUMO=',I8,' , gap=',F10.3,' eV')
 9210 FORMAT(6F13.6)
 9220 format(1x,A5,' FMO1 energy summary (eV):',
     *      /5x,'HOMO=  ',F12.6,
     *      /5x,'LUMO=  ',F12.6,
     *      /5x,'Fermi= ',F12.6,
     *      /5x,'gap=   ',F12.6)
      END
C*MODULE fmoio   *DECK readovd
      SUBROUTINE readovd(d,orbxch,odexch,na,l0p,l1,
     *                   iodfmo,idmrec0,modov)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      dimension d(*),iodfmo(*)
      logical orbxch,odexch
      common /fmorun/ espscf,e0scf(2),emp2s,IDAFMO,icurfg,jcurfg,kcurfg,
     *                icurlay,icurunt,nat1e,ncursh,ngau,icurpop,ifmostp,
     *                moncor,needr,modrst,norbproj,nunesp,iskipesp,
     *                IESDPPC,idoprop,mp2run,icurit,idmfmo,iddfmo,
     *                iddcur,nddleft,ivmfmo,nzmtfmo,ifmobas,itmfmo(2)
c
c     read occ+virt density, based on 3 ways of storing:
c     1. density
c     2. orbitals
c     3. density+orbitals
c     d thus must be allocated appropriate space (l2 or l2+l3).
c     For ROHF/UHF, this becomes
c     1. density(alpha)+density(beta)
c     2. the same for orbitals is not implemented (no DFT!)
c     wrk(l3)
c
      l2=(l1*l1+l1)/2
      l3=l1*l1
      if(orbxch) then
        if(odexch) then
          call rareads(IDAFMO,iodfmo,d,l2+l3,idmrec0,0)
        else
          call rareads(IDAFMO,iodfmo,d(l2+1),l3,idmrec0,0)
        endif
        if(modov.eq.1) then
          inid=1
          nda=na
        else if(modov.eq.2) then
          inid=na*l1+1
          nda=l0p-na
        else
          inid=1
          nda=l0p
        endif
        ndb=nda
c       not nb! na to have the same projectors for alpha and beta
c       This also works for MCSCF, by taking the active space as "occupied".
        call DMTX2(d,d(l2+inid),nda,l1,l1,ndb)
c       write(6,*) 'wwwid',na,l0p,l1,nda,inid
        call prsq(d(l2+inid),nda,l1,l1)
        call prtri(d,l1)
      else
        if(modov.eq.1) then
c         This is not right for ROUHF but we cannot get here then.
          call rareads(IDAFMO,iodfmo,d,l2,idmrec0,0)
        else
          call abrtx("Error in READOVD.") 
        endif
      endif
      RETURN
      END
C*MODULE fmoio   *DECK setdmpnt
      SUBROUTINE setdmpnt(ioptdm,esdder,maxl30,idmpnt,nexrst,ndmsiz)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      logical esdder
      dimension idmpnt(*),maxl30(*)
      Common /fmoinf/ nfg,nlayer,natfmo,nbdfg,naotyp,nbody,nsegm
c
c     Record 1 is used for restarts, length nexrst.
c
c     ioptdm=1 supervector
c            2 matrix 
c     on the output, ndmsiz has a different meaning
c     ioptdm=1 size of supervector (for one set) 
c            2 number of matrix columns (for one set)
c              (there are 2 or 3 sets)
c     Similarly, idmpnt refers either to 
c     ioptdm=1 supervector elements
c            2 matrix columns
c
c     define set 1
      idmpnt(1)=1
      if(ioptdm.eq.1) then
        ndmsiz=0
        do ifg=1,nfg
          idmpnt(ifg+1)=ndmsiz+1+nexrst
          ndmsiz=ndmsiz+maxl30(ifg)
c         write(6,*) 'wwwrec',ifg,maxl30(ifg)
        enddo
      else
        do ifg=1,nfg
          idmpnt(ifg+1)=ifg+nexrst
        enddo
        ndmsiz=nfg
      endif
c
c     Do mapping of records to save space.
c     The order of records is
c     set   F40        DDI memory
c     1.   dens1       dens1 
c     2.   grdV        dens2
c     3.   dens2       grdV (undefined unless esdder is true) 
c     4.   unused      undefined
c
c
c     do iset=1,4-1
c       do ifg=1,nfg
c         idmpnt(ifg+1+nfg*iset)=idmpnt(ifg+1)+ndmsiz*iset
c       enddo
c     enddo 
c
c     define set 3
      do ifg=1,nfg
        idmpnt(ifg+1+nfg*2)=idmpnt(ifg+1)+ndmsiz
      enddo
      if(esdder) then
c     define set 2, and store 0 to set 4 
      do ifg=1,nfg
        idmpnt(ifg+1+nfg)=idmpnt(ifg+1)+ndmsiz*2
        idmpnt(ifg+1+nfg*3)=0
      enddo
      else
c     store 0 to sets 3,4 
        do ifg=1,nfg
          idmpnt(ifg+1+nfg)=0
          idmpnt(ifg+1+nfg*3)=0
        enddo
      endif
c
      RETURN
      END
C*MODULE fmoio   *DECK clrvesd
      SUBROUTINE clrvesd(ioptdm,idmfmo,wrk,maxl30,idmpnt) 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      LOGICAL GOPARR,DSKWRK,MASWRK,ISGDDI,PAROUT,INITGDDI,wasgddi,
     *        MLGDDI
      integer ddi_world,ddi_group
      Parameter(ddi_world=0,ddi_group=1)
      dimension wrk(*),idmpnt(*),maxl30(*)
      COMMON /GDDI/   ISCOPE,NGROUPS,MYGROUP,MEGLOB,npglob,nnglob,JBTYP,
     *                ISGDDI,PAROUT,INITGDDI,wasgddi,MLGDDI,NSUBGR,
     *                MeUniv,NPUniv,numdlb,myworld,nworlds,mogddi
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      Common /fmoinf/ nfg,nlayer,natfmo,nbdfg,naotyp,nbody,nsegm
c
c     clear the RESDIM accumulated potentials, stored as set 2
c
      if(isgddi) call gddi_scope(ddi_world)
      do ifg=1,nfg
        if(mod(ifg,nproc).eq.ME) then
          nrec=1+ifg+nfg
          len=maxl30(ifg)
          call vclr(wrk,1,len)
          ist=idmpnt(nrec)
          if(ioptdm.eq.1) then
            iend=ist+len-1
            call ddi_put(idmfmo,1,1,ist,iend,wrk)
          else
            call ddi_put(idmfmo,1,len,ist,ist,wrk)
c           we do not zero out maxl1 here. Should be fine?
          endif
        endif 
      enddo 
      if(isgddi) call gddi_scope(ddi_group)
      RETURN
      END
C*MODULE fmoio   *DECK setddpnt
      SUBROUTINE setddpnt(orbxch,numfrg,layfrg,scffrg,iactfg,iddpnt,
     *                    nddsiz)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      logical orbxch,highdim
      dimension numfrg(*),layfrg(*),scffrg(*),iactfg(*),iddpnt(*)
      Common /fmoinf/ nfg,nlayer,natfmo,nbdfg,naotyp,nbody,nsegm
      common /fmoopt/ espsca(9),RESPAP(2),rESPPC(2),rESDIM,restri(4),
     *                rcorsd,respct,convfg,cnvdmp,coroff,rflmo(4),
     *                orshft,orshft2,cnvafo,ascreen(4),IXESP,mxitfg,
     *                nguess,NBSSE,modorb,modpar,irststp,irstlay,nprfmo,
     *                nfmopal,modprp,maxl1c,ipieda,modgrd,modesp,ivmul,
     *                modlmo,nopden,mofock,modfd,modfmm,ncentm,ndualb,
     *                ngab,modpan
      DATA RMC/8HMCSCF   /
c
      iddpnt(1)=1
      nddsiz=0
      loop=0
      do ifg=1,nfg
        do jfg=1,ifg-1
          rrij=fmodist(ifg,0,0,jfg)
          l1=iand(numfrg(ifg),65535)+iand(numfrg(jfg),65535)
c         This is an overestimate in general.
          loop=loop+1
          highdim=modfd.eq.0.or.layfrg(ifg)+layfrg(jfg).ge.4
          if(iand(modfd,2).ne.0) highdim=highdim.and.
     *                                   iactfg(ifg)+iactfg(jfg).gt.0
          if((rrij.le.resdim.or.resdim.eq.0).and.highdim) then
            if(orbxch) then
              m2=l1*l1
            else 
              m2=(l1*l1+l1)/2
            endif
            if(scffrg(ifg).eq.rmc.or.scffrg(jfg).eq.rmc) 
     *        m2=(l1*l1+l1)/2 + l1*l1
            if(scffrg(ifg).eq.rmc.or.scffrg(jfg).eq.rmc) 
     *        write(6,*) 'wwwall',ifg,jfg
            iddpnt(loop)=nddsiz+1
            nddsiz=nddsiz+m2
c           write(6,*) 'wwwidd',loop,iddpnt(loop)
          else
            iddpnt(loop)=-1
          endif
c         write(6,*) 'wwwidd',loop,iddpnt(loop)
        enddo
      enddo
c
      RETURN
      END
c
C*MODULE fmoio   *DECK setvmpnt
      SUBROUTINE setvmpnt(numfrg,ivmpnt,nvmsiz)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      dimension numfrg(*),ivmpnt(*)
      Common /fmoinf/ nfg,nlayer,natfmo,nbdfg,naotyp,nbody,nsegm
c
      nvmsiz=0
      do ifg=1,nfg
        l1=iand(numfrg(ifg),65535)
        l2=(l1*l1+l1)/2
        ivmpnt(ifg)=nvmsiz+1
        nvmsiz=nvmsiz+l2
c         write(6,*) 'wwwidd',ifg,ivmpnt(ifg)
      enddo
      do ifg=1,nfg
        ivmpnt(ifg+nfg)=ivmpnt(ifg)+nvmsiz
        ivmpnt(ifg+nfg*2)=ivmpnt(ifg)+nvmsiz*2
      enddo
c
      RETURN
      END
c
C*MODULE fmoio   *DECK getdd
      SUBROUTINE getdd(ifg,jfg,mcdim,orbxch,dd,l1,iddpnt,ires)
      use mx_limits, only: mxatm
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      dimension dd(*),iddpnt(*)
      logical orbxch,ISGDDI,PAROUT,INITGDDI,GOPARR,DSKWRK,MASWRK,wasgddi
     *       ,MLGDDI,mcdim
      Integer ddi_world 
      Parameter(ddi_world=0)
      COMMON /GDDI/   ISCOPE,NGROUPS,MYGROUP,MEGLOB,NPGLOB,NNGLOB,JBTYP,
     *                ISGDDI,PAROUT,INITGDDI,wasgddi,MLGDDI,NSUBGR,
     *                MeUniv,NPUniv,numdlb,myworld,nworlds,mogddi
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      common /fmorun/ espscf,e0scf(2),emp2s,IDAFMO,icurfg,jcurfg,kcurfg,
     *                icurlay,icurunt,nat1e,ncursh,ngau,icurpop,ifmostp,
     *                moncor,needr,modrst,norbproj,nunesp,iskipesp,
     *                IESDPPC,idoprop,mp2run,icurit,idmfmo,iddfmo,
     *                iddcur,nddleft,ivmfmo,nzmtfmo,ifmobas,itmfmo(2)
c
c     ires=0 success (got from DDI and saved on DA)
c     ires=1 failure (done nothing)
c     Failures are not fatal: backup copy will be used (sum of monomers).
c
      ii=max(ifg,jfg)-1
      jj=min(ifg,jfg)
      loop=(ii*ii-ii)/2+jj
      itmp=iscope
      ist=iddpnt(loop)
      if(ist.le.0) then
        if(maswrk) write(iw,9000) ifg,jfg
        ires=1
        return
      endif
      ires=0 
      if(orbxch) then
        len=l1*l1
        irec=15
        noff=(l1*l1+l1)/2
      else
        len=(l1*l1+l1)/2
        irec=16
c       Always save to record 16 where SCF expects it.
        noff=0
      endif
      if(mcdim) then
        len=(l1*l1+l1)/2+l1*l1
        noff=0
      endif
      iend=ist+LEN-1
      if(isgddi) call gddi_ascope(ddi_world)
      call ddi_get(iddfmo,1,1,ist,iend,dd(1+noff))
      if (isgddi.and.itmp.ne.ddi_world) call gddi_ascope(itmp)
c     write(6,*) 'wwwgetting/dawriting',ist,iend,irec
      if(mcdim) then
        len=(l1*l1+l1)/2
        CALL dawrit(IDAF,IODA,dd,len,16,0)
        CALL dawrit(IDAF,IODA,dd(len+1),l1*l1,15,0)
        len=len+l1*l1
      else
      CALL dawrit(IDAF,IODA,dd(1+noff),len,irec,0)
c
c     Next, also create and save the density.
c 
      if(orbxch) then
c       This will not work for MCSCF, but will for RHF or ROHF (+DFT or MP2).
c       (probably, not for ROHF?)
        call DMTX2(dd,dd(noff+1),na,l1,l1,nb)
        CALL dawrit(IDAF,IODA,dd,noff,16,0)
c       write(6,*) 'wwwsaved D'
      endif
      endif
c
      RETURN
 9000 FORMAT(/1X,'Dimer',2I6,' : monomer densities will be used in',
     *           ' dimer guess.')
      END
c
C*MODULE fmoio   *DECK putdd
      SUBROUTINE putdd(ifg,jfg,mcdim,irecd,orbxch,dd,l1,iddpnt,ires)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      dimension dd(*),iddpnt(*)
      logical orbxch,ISGDDI,PAROUT,INITGDDI,GOPARR,DSKWRK,MASWRK,wasgddi
     *       ,MLGDDI,mcdim
      Integer ddi_world
      Parameter(ddi_world=0)
      COMMON /GDDI/   ISCOPE,NGROUPS,MYGROUP,MEGLOB,NPGLOB,NNGLOB,JBTYP,
     *                ISGDDI,PAROUT,INITGDDI,wasgddi,MLGDDI,NSUBGR,
     *                MeUniv,NPUniv,numdlb,myworld,nworlds,mogddi
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      common /fmorun/ espscf,e0scf(2),emp2s,IDAFMO,icurfg,jcurfg,kcurfg,
     *                icurlay,icurunt,nat1e,ncursh,ngau,icurpop,ifmostp,
     *                moncor,needr,modrst,norbproj,nunesp,iskipesp,
     *                IESDPPC,idoprop,mp2run,icurit,idmfmo,iddfmo,
     *                iddcur,nddleft,ivmfmo,nzmtfmo,ifmobas,itmfmo(2)
c
c     ires=0 success (read from DA and put to DDI)
c     ires=1 failure (done nothing)
c     Failures are not fatal: backup copy will be used (sum of monomers).
c
      ii=max(ifg,jfg)-1
      jj=min(ifg,jfg)
      loop=(ii*ii-ii)/2+jj
      itmp=iscope
      ist=iddpnt(loop)
      if(orbxch) then
        len=l1*l1
        irec=15
      else
        len=(l1*l1+l1)/2
c       irec=16
        irec=irecd
c       Read SCF density.
      endif
      if(mcdim) len=(l1*l1+l1)/2+l1*l1
      if(ist.le.0) then
c       Try to allocate memory from emergency funds
        if(len.le.nddleft) then
          ist=iddcur
          iddpnt(loop)=ist
          iddcur=iddcur+len
          nddleft=nddleft-len
          if(maswrk) write(iw,9000) ifg,jfg,len,nddleft
c         write(6,*) 'wwwi',ist
        else
c         if(maswrk) write(iw,9010) ifg,jfg,len,nddleft
c         call abrt
          ires=1
          return
        endif
      endif
      ires=0
      if(mcdim) then
        len=(l1*l1+l1)/2
        CALL daread(IDAF,IODA,dd,len,16,0)
        CALL daread(IDAF,IODA,dd(len+1),l1*l1,15,0)
        len=len+l1*l1
      else
        CALL daread(IDAF,IODA,dd,len,irec,0)
      endif
      iend=ist+LEN-1
      if(isgddi) call gddi_ascope(ddi_world)
      call ddi_put(iddfmo,1,1,ist,iend,dd)
c     write(6,*) 'wwwputting/dareading',ist,iend,irec
      if (isgddi.and.itmp.ne.ddi_world) call gddi_ascope(itmp)
      RETURN
 9000 FORMAT(/1X,2I6,': allocating',i8,' from an offshore fund, ',
     *           'remainder=',i8)
c9010 FORMAT(/1X,2I6,': unable to allocate',i6,' words from an ',
c    *           'offshore fund, left',i6,
c    *       /1x,'Monomer densities will be used.')
      END
C*MODULE fmoio   *DECK retdd
      SUBROUTINE retdd 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      logical ISGDDI,PAROUT,INITGDDI,GOPARR,DSKWRK,MASWRK,wasgddi,MLGDDI
      COMMON /GDDI/   ISCOPE,NGROUPS,MYGROUP,MEGLOB,NPGLOB,NNGLOB,JBTYP,
     *                ISGDDI,PAROUT,INITGDDI,wasgddi,MLGDDI,NSUBGR,
     *                MeUniv,NPUniv,numdlb,myworld,nworlds,mogddi
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      common /fmorun/ espscf,e0scf(2),emp2s,IDAFMO,icurfg,jcurfg,kcurfg,
     *                icurlay,icurunt,nat1e,ncursh,ngau,icurpop,ifmostp,
     *                moncor,needr,modrst,norbproj,nunesp,iskipesp,
     *                IESDPPC,idoprop,mp2run,icurit,idmfmo,iddfmo,
     *                iddcur,nddleft,ivmfmo,nzmtfmo,ifmobas,itmfmo(2)
      Integer ddi_world
      Parameter(ddi_world=0)
c
      if(idmfmo.ge.0.or.iddfmo.ge.0.or.ivmfmo.ge.0) then
        if(isgddi) call gddi_scope(ddi_world)
c       The order of destruction should be exactly opposite to creation!
        if(maswrk) write(iw,*) 'Returned DDI',idmfmo,ivmfmo,iddfmo
        if(iddfmo.ge.0) then
          call ddi_destroy(iddfmo)
          iddfmo=-1
        endif
        if(ivmfmo.ge.0) then
          call ddi_destroy(ivmfmo)
          ivmfmo=-1
        endif
        if(idmfmo.ge.0) then
          call ddi_destroy(idmfmo)
          idmfmo=-1
        endif
c       if(isgddi) call gddi_scope(ddi_group)
c       This subroutine is supposed to be called from the world scope
c       outside of FMO, although this is not required, but no switch is made
c       to whatever scope it was.
      endif
      RETURN
      END
C*MODULE fmoio   *DECK takevm
      SUBROUTINE takevm(mode,ifg,irec0,vv,l2,ivmpnt)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      dimension vv(*),ivmpnt(*)
      logical ISGDDI,PAROUT,INITGDDI,GOPARR,DSKWRK,MASWRK,wasgddi,MLGDDI
      Integer ddi_world
      Parameter(ddi_world=0)
      COMMON /GDDI/   ISCOPE,NGROUPS,MYGROUP,MEGLOB,NPGLOB,NNGLOB,JBTYP,
     *                ISGDDI,PAROUT,INITGDDI,wasgddi,MLGDDI,NSUBGR,
     *                MeUniv,NPUniv,numdlb,myworld,nworlds,mogddi
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      Common /fmoinf/ nfg,nlayer,natfmo,nbdfg,naotyp,nbody,nsegm
      common /fmorun/ espscf,e0scf(2),emp2s,IDAFMO,icurfg,jcurfg,kcurfg,
     *                icurlay,icurunt,nat1e,ncursh,ngau,icurpop,ifmostp,
     *                moncor,needr,modrst,norbproj,nunesp,iskipesp,
     *                IESDPPC,idoprop,mp2run,icurit,idmfmo,iddfmo,
     *                iddcur,nddleft,ivmfmo,nzmtfmo,ifmobas,itmfmo(2)
c
c     mode=0 put
c     mode=1 get
c
      if(irec0.lt.0) then
c       one set of potentials
        nrec=nfg*2
      else
c       two sets of densities
        nrec=0
        if(irec0.gt.1) nrec=nfg
c       revert reading and writing records 
        if(mode.eq.1) nrec=nfg-nrec
      endif
      len=l2
      ist=ivmpnt(ifg+nrec)
      iend=ist+LEN-1
c     It is usual to call this in group scope, and it is mandatory that
c     vv on all nodes in the current set (group or world) is identical (mode=0).
      if(maswrk) then
        itmp=iscope
        if(isgddi.and.itmp.ne.ddi_world) call gddi_ascope(ddi_world)
c       write(6,*) 'wwwtaking',mode,ist,iend,ifg,irec0,nrec
        if(mode.eq.0) call ddi_put(ivmfmo,1,1,ist,iend,vv)
        if(mode.eq.1) call ddi_get(ivmfmo,1,1,ist,iend,vv)
        if (isgddi.and.itmp.ne.ddi_world) call gddi_ascope(itmp)
      endif 
      if(mode.eq.1) call DDI_BCAST(2422,'F',vv,len,MASTER)
      RETURN
      END
C*MODULE fmoio   *DECK qmcfout
C>
C>    @brief Output data for QMC.
C>
C>    @details Print summary for QMC. 
C>
C>    @author Dmitri Fedorov
C>
      SUBROUTINE qmcfout(nftqmc,ifg,jfg,kfg,q) 
      use mx_limits, only: mxatm,mxao
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      dimension q(*)
      Parameter(UNITS=0.52917724924D+00)
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /RUNLAB/ TITLE(10),ANAM(MXATM),BNAM(MXATM),BFLAB(MXAO)
      write(nftqmc,9000) ifg,jfg,kfg,nat,num,     mul,ne,ich,na,nb
      if(jfg.eq.0.and.kfg.eq.0) then
      write(nftqmc,9100) (ANAM(J),ZAN(J),(C(I,J)*UNITS,I=1,3),
     *                   ZAN(J)-q(j),J=1,nat)
      else
      write(nftqmc,9150) (ANAM(J),ZAN(J),(C(I,J)*UNITS,I=1,3),J=1,nat)
      endif
      write(nftqmc,9200)
      RETURN
 9000 format(/1x,60(1H-),/1x,'QMC CURRENT N-MER COORDINATES, I=',I7,
     *           ' J=',I7,' K=',I7,/1x,60(1H-),/,
     *       /1x,'QMC NUMBER OF ATOMS IN FRAGMENT                =',I10,
     *       /1x,'QMC NUMBER OF CARTESIAN ATOMIC ORBITALS        =',I10,
     *       /1x,'QMC SPIN MULTIPLICITY                          =',I10,
     *       /1x,'QMC NUMBER OF ELECTRONS                        =',I10,
     *       /1x,'QMC CHARGE OF MOLECULE                         =',I10,
     *       /1x,'QMC NUMBER OF OCCUPIED ORBITALS (ALPHA)        =',I10,
     *       /1x,'QMC NUMBER OF OCCUPIED ORBITALS (BETA)         =',I10,
     *       /1x,'QMC ATOM CHARGE',8x,'X',17x,'Y',17x,'Z (ANGST)   ',
     *           'ESP CHARGE')
c    *       /1x,'QMC TOTAL NUMBER OF MOS IN VARIATION SPACE     =',I10,
 9100 FORMAT(1X,'QMC',A4,F5.1,3F18.10,F16.12)
 9150 FORMAT(1X,'QMC',A4,F5.1,3F18.10)
 9200 FORMAT(1X,'QMC',60(1H-))
      END
C*MODULE fmoio   *DECK abortdump
C>
C>     @brief Dump info when aborting.
C>
C>     @details Write out data when aborting.
C>
C>     @author  Dmitri Fedorov
C>
      SUBROUTINE abortdump
      use mx_limits, only: mxatm
      USE comm_EFPFMO, only: runefp
      IMPLICIT NONE
C     Declarations of common blocks
      DOUBLE PRECISION X
      INTEGER nfg,nlayer,natfmo,nbdfg,naotyp,nbody,nsegm
      INTEGER lichfg,lmulfg,lidmrec,lfrgnam,llayfrg,lindat
      INTEGER lnCBS,lfmozan,lfmoc,lfmomas,lizbas,liaglob,libdgh
      INTEGER liabdfg,ljabdfg,lnCAO,lidxCAO,liaprjo,ljaprjo
      INTEGER lCoreAO,lOccCor,lshiftb,liodfmo,lfmoda,lfmodb
      INTEGER lfmoespa,lfmoespb,llocfmo,lscffrg,lfmoscf,lrij
      INTEGER lpopmul,lpopmat,lialoc,lindbd,liatfrg,lindfrg
      INTEGER lindgfrg,lnatfrg,lnat0frg,lianfrg,lzanfrg,lcfrg
      INTEGER llibish,llibnsh,llibng,lindatg,lfmobuf,lfmode
      INTEGER lnumfrg,lloctat,liaoglob,lloadm,lfmoge,ldgrid
      INTEGER liodcfmo,ljob2grp,lfmopg,lemocdr,luntxyz,luntrot
      INTEGER lstonep,lmapsu,lfrgmul,lclmo,lialmo,lindlmo
      INTEGER latclmo,llmobdf,lfgflmo,lnfglmo,llfglmo,lpfglmo
      INTEGER LPOPDMAT,lidmpnt,liddpnt,livmpnt,liactfg,lcrfrg
      INTEGER lzlmfrgv,lYlmfrgv,lndtfrg,lf_mm,lg_mm,lmaxl30
      INTEGER libuffg,lindatp,lsmon,lsdim,ledimfed,lexcit3d,lifgfret
      DOUBLE PRECISION espscf,e0scf,emp2s
      INTEGER IDAFMO,icurfg,jcurfg,kcurfg
      INTEGER icurlay,icurunt,nat1e,ncursh,ngau,icurpop,ifmostp
      INTEGER moncor,needr,modrst,norbproj,nunesp,iskipesp
      INTEGER IESDPPC,idoprop,mp2run,icurit,idmfmo,iddfmo
      INTEGER iddcur,nddleft,ivmfmo,nzmtfmo,ifmobas,itmfmo
      INTEGER ISCOPE,NGROUPS,MYGROUP,MEGLOB,npglob,nnglob,JBTYP
      INTEGER NSUBGR
      INTEGER MeUniv,NPUniv,numdlb,myworld,nworlds,mogddi
      LOGICAL ISGDDI,PAROUT,INITGDDI,wasgddi,MLGDDI
      DOUBLE PRECISION ZAN,C
      INTEGER NAT,ICH,MUL,NUM,NQMT,NE,NA,NB
      INTEGER IAN
      INTEGER IR,IW,IP,IJK,IJKT,IDAF,NAV,IODA
      INTEGER ME,MASTER,NPROC,IBTYP,IPTIM
      LOGICAL GOPARR,DSKWRK,MASWRK
C     Other declarations
      DOUBLE PRECISION, PARAMETER :: UNITS=0.52917724924D+00
      LOGICAL out
      DOUBLE PRECISION FMOMD
      INTEGER i,j
C
      COMMON /FMCOM / X(1)
      common /fmoinf/ nfg,nlayer,natfmo,nbdfg,naotyp,nbody,nsegm
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
      common /fmorun/ espscf,e0scf(2),emp2s,IDAFMO,icurfg,jcurfg,kcurfg,
     *                icurlay,icurunt,nat1e,ncursh,ngau,icurpop,ifmostp,
     *                moncor,needr,modrst,norbproj,nunesp,iskipesp,
     *                IESDPPC,idoprop,mp2run,icurit,idmfmo,iddfmo,
     *                iddcur,nddleft,ivmfmo,nzmtfmo,ifmobas,itmfmo(2)
      COMMON /GDDI/   ISCOPE,NGROUPS,MYGROUP,MEGLOB,npglob,nnglob,JBTYP,
     *                ISGDDI,PAROUT,INITGDDI,wasgddi,MLGDDI,NSUBGR,
     *                MeUniv,NPUniv,numdlb,myworld,nworlds,mogddi
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /IOFILE/ IR,IW,IP,IJK,IJKT,IDAF,NAV,IODA(950)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      DATA FMOMD /8HMD      /
c
C     Dumps data to show where the job aborts because some
c     reduced output options may obscure the cause of abort.
c     GDDI/PAROUT runs open output files on each core; write on slaves
c     too in that case.
c
      out=maswrk.or.isgddi.and.parout
      if(nfg.ne.0.and.out) then
c       For now, output only for the smallest FMO output.
        write(iw,9000) ifmostp,icurlay,icurfg,jcurfg,kcurfg,icurit,
     *                 nat,ICH,MUL,NUM,NQMT,NE,NA
c       if doing a monomer, dimer or trimer
        if(ifmostp.eq.2.or.ifmostp.eq.4.or.ifmostp.eq.9) then
         if(maswrk) write(iw,*) 'Current n-mer'
         if(maswrk) write(iw,9100) (iAN(i),(C(j,i)*units,j=1,3),i=1,nat)
        endif 
        if(runefp.eq.fmomd) then
         write(iw,9110) natfmo
         if(maswrk) write(iw,9100) (int(x(lfmozan+i-1)),
     *                    (x(lfmoc+(i-1)*3+j-1)*units,j=1,3),i=1,natfmo)
        endif
      endif
      RETURN
 9000 format(/1x,'FMO dump: step',I2,' layer',I2,' n-mer',3I8,' iter',
     *           I3,
     *      /1x,'NAT,ICH,MUL,NUM,NQMT,NE,NA=',2I3,5I6)
 9100 format(1x,I3,3F18.10)
 9110 format(1x,'All FMO atoms',I9)
      END
c
C*MODULE fmoio   *DECK pieda3
C>
C>     @brief PIEDA3 calculations.
C>
C>     @details Add corrections.
C>
C>     @author Dmitri Fedorov
C>
      SUBROUTINE pieda3(ipieda,dftbfl,dodft,gcorrel,dodisp,dopcm,nfg,
     *                nfg2d,nfg2,nfg3,nedimex,edim,etrim,ecorrdft,esolv)
      IMPLICIT NONE
      LOGICAL dftbfl,dodft,gcorrel,dodisp,dopcm
      INTEGER ipieda,nfg,nfg2d,nfg2,nfg3,nedimex,ifg,ijfg,ijkfg,ikfg,
     *        jfg,jkfg,kfg,nijk
      DOUBLE PRECISION, PARAMETER :: tokcal=627.51D+00
      DOUBLE PRECISION edim(nfg2d,*),etrim(nfg3,*),ecorrdft(*),esolv(*),
     *                 etrimct,etrimcts,etrimdi,etrimex,etrimrc,etrimso,
     *                 etrimsok,pie3
c
c     Add three-body corrections to PIEs and/or their components.
c
      etrimct=0
      etrimex=0
      etrimrc=0
      etrimdi=0
      etrimso=0
c     tokcal3=3*tokcal
      ijkfg=0
      nijk=0
      do ifg=1,nfg
        do jfg=1,ifg-1
          ijfg=(ifg*ifg-3*ifg)/2+jfg+1
          do kfg=1,jfg-1
            ijkfg=ijkfg+1
            if(ipieda.gt.0.and.etrim(ijkfg,1).eq.1 .or.
     *         ipieda.eq.0.and.etrim(ijkfg,1).ne.0) then
              nijk=nijk+1
              ikfg=(ifg*ifg-3*ifg)/2+kfg+1
              jkfg=(jfg*jfg-3*jfg)/2+kfg+1
              pie3=0
              if(ipieda.gt.0) then
                if(dftbfl) then
c               
c                 CT*ES for DFTB
c               
                  etrimct=ecorrdft(nfg*2+nfg2+ijkfg)/3
                  ecorrdft(nfg*2+ijfg)=ecorrdft(nfg*2+ijfg)+etrimct
                  ecorrdft(nfg*2+ikfg)=ecorrdft(nfg*2+ikfg)+etrimct
                  ecorrdft(nfg*2+jkfg)=ecorrdft(nfg*2+jkfg)+etrimct
                  pie3=pie3+etrimct
                  etrimcts=etrimct
                else
c               
c                 EX
c               
                  etrimex=etrim(ijkfg,3)/3
                  edim(ijfg,nedimex)=edim(ijfg,nedimex)+etrimex
                  edim(ikfg,nedimex)=edim(ikfg,nedimex)+etrimex
                  edim(jkfg,nedimex)=edim(jkfg,nedimex)+etrimex
                  pie3=pie3+etrimex
c               
                  if(dodft) then
c               
c                   RC for DFT
c               
                    etrimrc=ecorrdft(nfg+nfg2+ijkfg)/3
                    ecorrdft(ijfg+nfg)=ecorrdft(ijfg+nfg)+etrimrc
                    ecorrdft(ikfg+nfg)=ecorrdft(ikfg+nfg)+etrimrc
                    ecorrdft(jkfg+nfg)=ecorrdft(jkfg+nfg)+etrimrc
                    pie3=pie3+etrimrc
                  else if(gcorrel) then
c               
c                   DI for MP2/CC
c                   Subtract, because edim(ijfg,3) is subtracted to get DI.
c               
                    etrimdi=etrim(ijkfg,4)/3
                    edim(ijfg,3)=edim(ijfg,3)-etrimdi
                    edim(ikfg,3)=edim(ikfg,3)-etrimdi
                    edim(jkfg,3)=edim(jkfg,3)-etrimdi
                    pie3=pie3+etrimdi
                  endif
                endif
c
c               0 for DFTB, CT+MIX otherwise
c               By modifying edim(,1) the total solute PIE is corrected,
c               which also corrects for 0/CT+mix.
c               
                etrimct=etrim(ijkfg,2)/3
                pie3=pie3+etrimct
                edim(ijfg,1)=edim(ijfg,1)+pie3
                edim(ikfg,1)=edim(ikfg,1)+pie3
                edim(jkfg,1)=edim(jkfg,1)+pie3
c               
c               For MP2/CC/D also change uncorrelated edim(,3)
c               
                if(gcorrel.or.dodisp) then
                  edim(ijfg,3)=edim(ijfg,3)+pie3
                  edim(ikfg,3)=edim(ikfg,3)+pie3
                  edim(jkfg,3)=edim(jkfg,3)+pie3
                endif
                if(dftbfl) then
                  etrimex=etrimct
                  etrimct=etrimcts
                endif
              else
                pie3=etrim(ijkfg,1)/3
                edim(ijfg,1)=edim(ijfg,1)+pie3
                edim(ikfg,1)=edim(ikfg,1)+pie3
                edim(jkfg,1)=edim(jkfg,1)+pie3
                pie3=etrim(ijkfg,2)/3
                edim(ijfg,2)=edim(ijfg,2)+pie3
                edim(ikfg,2)=edim(ikfg,2)+pie3
                edim(jkfg,2)=edim(jkfg,2)+pie3
                if(gcorrel) then
                  pie3=etrim(ijkfg,3)/3
                  edim(ijfg,3)=edim(ijfg,3)+pie3
                  edim(ikfg,3)=edim(ikfg,3)+pie3
                  edim(jkfg,3)=edim(jkfg,3)+pie3
                endif
              endif
c
c             SOLV
c
              if(dopcm) then
                etrimso=esolv(nfg*2+nfg2*2+ijkfg)/3
                etrimsok=etrimso*tokcal
c               The dreaded units: esolv(nfg*2+nfg2+) are in kcal/mol
c                                  esolv(nfg*2+nfg2*2+ijkfg) is in hartree
                esolv(nfg*2+nfg2+ijfg)=esolv(nfg*2+nfg2+ijfg)+etrimsok
                esolv(nfg*2+nfg2+ikfg)=esolv(nfg*2+nfg2+ikfg)+etrimsok
                esolv(nfg*2+nfg2+jkfg)=esolv(nfg*2+nfg2+jkfg)+etrimsok
                pie3=pie3+etrimso
              endif
c             write(6,9000) ifg,jfg,kfg,pie3*tokcal3,etrimex*tokcal3,
c    *                     etrimct*tokcal3,(etrimrc+etrimdi)*tokcal3,
c    *                      etrimso*tokcal3
c9000         format(1x,'DeltaE=',3I6,': tot EX CT RC+DI SO',5F10.5)
c             for DFTB: tot E0 CT*ES DI SO
            endif
          enddo
        enddo
      enddo
      write(6,*) 'Trimer contributions:',nijk
c
      RETURN
      END
C*MODULE fmoio   *DECK prifragq
C>
C>     @brief Print fragment charges.
C>
C>     @details Print charges of fragments.
C>
C>     @author  Dmitri Fedorov
C>
      SUBROUTINE prifragq(IW,dftbfl,nbodyc,IMCPFMO,indat,atmulq,fmozan,
     *                    fzcor)
      IMPLICIT NONE
C     Declarations of arguments
      INTEGER IW
      LOGICAL dftbfl
      INTEGER nbodyc,IMCPFMO,indat(*)
      DOUBLE PRECISION atmulq(natfmo,*),fmozan(*)
      DOUBLE PRECISION fzcor(*)
C     Declarations of common blocks
      DOUBLE PRECISION X
      INTEGER nfg,nlayer,natfmo,nbdfg,naotyp,nbody,nsegm
C     Other declarations
      INTEGER i,ifg,j,last,lct,lcti,loadfm,need
      DOUBLE PRECISION qi,zi
C
      COMMON /FMCOM / X(1)
      common /fmoinf/ nfg,nlayer,natfmo,nbdfg,naotyp,nbody,nsegm
c
      write(iw,9000)
      CALL VALFM(LOADFM)
      lct=LOADFM+1
      last=lct+nfg*nbodyc
      NEED=LAST-LOADFM-1
      CALL GETFM(NEED)
      call vclr(x(lct),1,nfg*nbodyc)
      do i=1,natfmo
        zi=fmozan(i)
        if(dftbfl) call dftb_nucz(zi,2)
        if(IMCPFMO.EQ.1) zi=zi-fzcor(I)
        ifg=indat(i)
        do j=1,nbodyc
          qi=zi-atmulq(i,j)
          lcti=lct+(j-1)*nfg+ifg-1
          x(lcti)=x(lcti)+qi
        enddo
      enddo
      do i=1,nfg
        write(iw,9100) i,(x(lct+(j-1)*nfg+i-1),j=1,nbodyc)
      enddo
      CALL RETFM(NEED)
 9000 format(/1x,'Fragment charges based on INDAT:',
     *       /1x,'     IFG      FMO1        FMO2        FMO3')
 9100 format(1x,I8,3F12.6)
      RETURN
      END
C*MODULE fmoio   *DECK pares
C>
C>     @brief Print results of partition analysis.
C>
C>     @details Print PA results.
C>
C>     @author  Dmitri Fedorov
C>
      subroutine pares(monrep,totnones,mempa,mempap,mbody,nsegm,nsegm2,
     *                 natsegm,smon,sdim,rdip,tt,idipor,totq,rappa,
     *                 doascan,esolvtot,dsolv,segnam,nsegnam,indatp1,
     *                 bbsplit,nat,zan,pop,indatp,patnam,modseg,modpan,
     *                 nsubsys,eints,eintpart,eint,domdan,emonmds,
     *                 edimmds,efmomd,autosa,isasign,nsegmasub,esatot,
     *                 samons,sadims,dftbfl)
      IMPLICIT NONE
      LOGICAL monrep,totnones,doascan,bbsplit,domdan,autosa,dftbfl
      INTEGER mempa,mempap,mbody,nsegm,natsegm(*),idipor,nsegnam,indatp1
     *       ,nat,indatp(*),modseg(*),modpan,nsubsys,IR,IW,IP,IJK,IJKT,
     *        IDAF,NAV,IODA,ISMX,i,is,j,js,k,loop,modsegi,nsegm2,isasign
     *       ,nsegmasub,isa,morigin
      DOUBLE PRECISION smon(16,*),sdim(8,*),rdip(3),tt,totq,rappa(2),
     *       esolvtot,dsolv(3),zan(*),pop(*),patnam(*),eints(*),totqs,
     *       eintpart(2,*),eint(nsubsys,*),SOLA,SOLB,SOLC,SOLG,SOLH,SOLN
     *      ,eesij,enones,eparti,esolv,etots,etotss,rij,rijf,totes,qq,
     *       tote(9),totdip(4,3),emonmds(nsegm,*),edimmds(nsegm2,*),
     *       efmomd(*),esatot,samons(2,*),sadims(*),r0(3),r0s(3),
     *       rdips(3),epseff
      character*8 segnam(*),segnami
      character*4 Scavcds
      character*30 dipor,scrtyp,reptyp,nones,inptyp
      DOUBLE PRECISION, PARAMETER :: TOKCAL=627.509469D+00
      DOUBLE PRECISION, PARAMETER :: UNITS=0.52917724924D+00
      DOUBLE PRECISION, PARAMETER :: one=1.0D+00
      COMMON /IOFILE/ IR,IW,IP,IJK,IJKT,IDAF,NAV,IODA(950)
      COMMON /SOLSMX/ SOLA,SOLB,SOLC,SOLG,SOLH,SOLN,ISMX
c
      call vclr(tote,1,9)
      if(mbody.eq.0) then
        write(iw,9000)
      else
        if(dftbfl) then
          write(iw,9002) mbody
        else
          write(iw,9003) mbody
        endif
      endif
      morigin=1
c     dipole origin:
c     0 solvent origin = solute origin
c     1 solvent and solute origins are separately set as their
c       respective centers of charge.
      dipor='Center of charge'
      if(idipor.eq.1) dipor='inherited from above'
      scrtyp='total'
      if(doascan) scrtyp='partial'
      reptyp='monomer+pairwise'
      if(monrep) reptyp='monomer'
      Scavcds='cav'
      IF(ISMX.NE.0) Scavcds='cds'
      if(.not.monrep) totnones=.false.
c     For convenience, split non-es if the user wanted REP,
c     so that the analysis is as detailed as possible.
      IF(ISMX.NE.0) totnones=.true.
c     SMD cannot split them in any case, only PCM can.
      nones='split'
      if(totnones) nones='merged'
      inptyp='manual'
      if(iand(modpan,32).ne.0) inptyp='FMO'
      if(indatp1.lt.0) then
        inptyp='PDB residues'
        if(bbsplit) inptyp='PDB backbones and side chains'
      endif
      write(iw,9005) nat,nsegm,nsubsys,totq,scrtyp,nones,reptyp,mempa,
     *               mempap,dipor,morigin,'partial',inptyp,rappa(2),tt
c     local model for solvent multipoles is not coded.
c
c     Atomic properties
c
      do i=1,nat
        is=indatp(i)
        segnami='        '
        if(is.le.nsegnam) segnami=segnam(is)
        write(iw,9011) i,patnam(i),indatp(i),segnami,zan(i),pop(i)
      enddo
c
c     Monomer energies
c
      if(dftbfl) write(iw,9006) 
      if(.not.dftbfl) write(iw,9009) 
      do i=1,nsegm
        etots=0
        do j=1,8
          tote(j)=tote(j)+smon(j+8,i)
          etots=etots+smon(j+8,i)
        enddo
        esolv=smon(13,i)+smon(14,i)+smon(15,i)+smon(16,i)
        segnami='        '
        if(i.le.nsegnam) segnami=segnam(i)
        if(dftbfl) then
          write(iw,9010) i,segnami,natsegm(i),(smon(j,i),j=1,2),
     *                   (smon(j,i),j=9,12),esolv,etots
        else
          epseff=0
          if(abs(smon(2,i)).lt.1.0D-08) epseff=1
          if(abs(smon(1,i)+smon(2,i)).gt.1.0D-06)
     *      epseff=smon(1,i)/(smon(1,i)+smon(2,i))
          write(iw,9010) i,segnami,natsegm(i),(smon(j,i),j=1,2),epseff
        endif
        if(domdan) emonmds(i,1)=etots
        isa=i+nsegmasub
        if(autosa) samons(2,isa)=samons(2,isa)+isasign*esolv
      enddo
c
c     Monomer dipoles and detailed solvent components
c
      call vclr(totdip,1,4*3)
      write(iw,9007) Scavcds,Scavcds
      totqs=0
      do i=1,nsegm
        write(iw,9012) i,modseg(i+nsegm),modseg(i),(smon(j,i),j=3,8),
     *                                    (smon(j,i)*tokcal,j=13,16)
c       j runs over x,y,z
        do j=1,3
          totdip(j,1)=totdip(j,1)+smon(2+j,i)
          totdip(j,2)=totdip(j,2)+smon(5+j,i)
          totdip(j,3)=totdip(j,3)+smon(2+j,i)+smon(5+j,i)
        enddo
        totqs=totqs+smon(2,i)
      enddo
c     Total dipole
      if(.not.doascan) then
        call daxpy(3,one,dsolv,1,totdip(1,2),1)
        call daxpy(3,one,dsolv,1,totdip(1,3),1)
      endif
c     j runs over solute,solvent,total
      do j=1,3
        totdip(4,j)=sqrt(totdip(1,j)*totdip(1,j)+totdip(2,j)*totdip(2,j)
     *                  +totdip(3,j)*totdip(3,j))
      enddo
c
      write(iw,9200) (rdip(i)*units,i=1,3),((totdip(i,j),i=1,4),j=1,3)
c
c     Renormalize dipoles with a shifted origin so that the total diole is 0.
c     
      Qq=totq+totqs
c     In principle, abs(Qq) may be checked but it is prone to containing
c     numerical errors comparable to the value itself, so check solute charge.
      if(abs(totq).gt.1.0D-03) then
        write(iw,9008)
        do j=1,3
c         r0(j)=totdip(j,3)/Qq
          r0(j)=totdip(j,1)/totq
          r0s(j)=totdip(j,2)/totqs
          if(morigin.eq.0) r0s(j)=r0(j)
          rdip(j)=rdip(j)+r0(j)
          rdips(j)=rdip(j)+r0s(j)
        enddo
        call vclr(totdip,1,4*3)
        do i=1,nsegm
          do j=1,3
            smon(j+2,i)=smon(j+2,i)-smon(1,i)*r0(j)
            smon(j+5,i)=smon(j+5,i)-smon(2,i)*r0s(j)
            totdip(j,1)=totdip(j,1)+smon(2+j,i)
            totdip(j,2)=totdip(j,2)+smon(5+j,i)
            totdip(j,3)=totdip(j,3)+smon(2+j,i)+smon(5+j,i)
          enddo
          write(iw,9013) i,(smon(j,i),j=3,8)
        enddo
        do j=1,3
        totdip(4,j)=sqrt(totdip(1,j)*totdip(1,j)+totdip(2,j)*totdip(2,j)
     *                  +totdip(3,j)*totdip(3,j))
        enddo
        write(iw,9210) (rdip(i)*units,i=1,3),(rdips(i)*units,i=1,3),
     *                 ((totdip(i,j),i=1,4),j=1,3)
      endif
c
      if(.not.dftbfl) then
        write(iw,9500)
        return
      endif
c
c     write(6,*) 'wwweee',tote(5),esolvtot
c     Total es value is in esolvtot.
      if(.not.doascan) tote(5)=tote(5)+esolvtot
      do i=1,8
        tote(9)=tote(9)+tote(i)
      enddo
      totes=tote(5)+tote(6)+tote(7)+tote(8)
      write(iw,9100) (1,tote(i),i=1,4),1,totes,1,tote(9)
c
c     if(mbody.gt.1) then
      if(nsegm.gt.1) then
      if(monrep) then
        if(totnones) then
          write(iw,9016)
        else
          write(iw,9014)
        endif
      else
        write(iw,9015)
      endif
      loop=0
      do i=1,nsegm
        do j=1,i-1
          loop=loop+1
          etots=0
          do k=3,6
            tote(k-1)=tote(k-1)+sdim(k,loop)
            etots=etots+sdim(k,loop)
          enddo
          do k=7,8
            tote(k)=tote(k)+sdim(k,loop)
            etots=etots+sdim(k,loop)
          enddo
          etots=etots*tokcal
          rij=sdim(1,loop)*units
          rijf=sdim(2,loop)
          if(monrep) then
            eesij=sdim(3,loop)*tokcal
            if(totnones) then
              enones=(sdim(7,loop)+sdim(8,loop))*tokcal
           write(iw,9020) i,j,rij,rijf,eesij,(sdim(k,loop)*tokcal,k=5,6)
     *                   ,enones,etots
            else
           write(iw,9020) i,j,rij,rijf,eesij,(sdim(k,loop)*tokcal,k=5,8)
     *                   ,etots
            endif
          else
           write(iw,9020) i,j,rij,rijf,(sdim(k,loop)*tokcal,k=3,8),etots
          endif
          if(domdan) then
            edimmds(loop,1)=rij
            call dcopy(6,sdim(3,loop),1,edimmds(loop,2),nsegm2)
          endif
          if(autosa.and.isasign.eq.1) sadims(loop)=etots/tokcal
        enddo
      enddo
      tote(9)=0
      do i=1,8
        tote(9)=tote(9)+tote(i)
      enddo
      totes=tote(5)+tote(6)+tote(7)+tote(8)
      write(iw,9100) (2,tote(i),i=1,4),2,totes,2,tote(9)
      endif
      if(domdan) efmomd(7)=tote(9)
      if(autosa) esatot=esatot+isasign*tote(9)
c
c     Subsystem analysis for PA
c
      if(nsubsys.ne.0) then
        write(iw,9300)
c       First, i,J (fragment-subsystem)
        etotss=0
        call vclr(eintpart,1,2*nsubsys)
        do i=1,nsegm
          is=modseg(i+nsegm)
          call vclr(eints,1,nsubsys)
          do j=1,nsegm
            if(i.ne.j) then
              if(i.gt.j) then
                loop=(i*i-3*i+2)/2+j
              else
                loop=(j*j-3*j+2)/2+i
              endif
              js=modseg(j+nsegm)
              etots=0
              do k=3,8
                etots=etots+sdim(k,loop)
              enddo
              eints(js)=eints(js)+etots 
c             write(6,*) 'www',i,j,js,loop,etots
            endif
          enddo
          etots=0
          do k=9,16
            etots=etots+smon(k,i)
          enddo
c         Halve self-interactions for is to avoid double counting and add it
c         to the partial energy.
c         Note that other eints(js) are still doubled counted.
          eparti=etots+eints(is)/2
          eints(is)=0
          write(iw,9310) i,is,eparti,(eints(js)*tokcal,js=1,nsubsys)
          etotss=etotss+eparti
          do js=1,nsubsys
            etotss=etotss+eints(js)/2
c           halve to correct double counting
            eint(is,js)=eint(is,js)+eints(js)
            eintpart(2,is)=eintpart(2,is)+eints(js)/2
          enddo
c         (1) is internal, (2) is partial
          eintpart(1,is)=eintpart(1,is)+eparti
          eintpart(2,is)=eintpart(2,is)+eparti
          isa=i+nsegmasub
          if(autosa) samons(1,isa)=samons(1,isa)+isasign*eparti
        enddo
        write(iw,9330) etotss
c       Now, I,J (subsystem-subsystem)
        write(iw,9400)
        etotss=0
        do is=1,nsubsys
          write(iw,9410) is,(eintpart(js,is),js=1,2),
     *                      (eint(is,js)*tokcal,js=1,nsubsys)
          etotss=etotss+eintpart(2,is)
        enddo
        write(iw,9330) etotss
      endif
      RETURN
 9000 format(/1x,73(1H-),/23x,'Partition analysis for full DFTB')
 9002 format(/1x,73(1H-),/23x,'Partition analysis for FMO',I1,'-DFTB')
 9003 format(/1x,73(1H-),/19x,'Partition analysis using FMO',I1,' data')
 9005 format(1x,'Number of atoms:',I13,
     *      /1x,'Number of segments:',I10,
     *      /1x,'Number of subsystems:',I8,
     *      /1x,'Total charge:',F25.8,
     *      /1x,'Screening type:',13x,A25,
     *      /1x,'Non-es solvent pair terms:',2x,A25,
     *      /1x,'Repulsion type:',13x,A25,
     *      /1x,'Replicated memory:',I14,
     *      /1x,'Distributed memory:',I13,
     *      /1x,'Multipole origin:',11x,A25,' ',I1,
     *      /1x,'Solvent dipole model: ',6x,A7,
     *      /1x,'Segment definition:',9x,A30,
     *      /1x,'Dispersion threshold:',F12.3,
     *      /1x,'Wall-clock time, s:',F12.1,/1x,73(1H-),
     *      //5x,'Atom',1x,'Name',1x,'Segment',1x,'Seg name',5x,'Z',
     *        3x,'charge')
 9006 format(/1x,'Segment properties (a.u.); ',
     *           'E=E0+EES+EREP+EDI+Esolv',
     *      /1x,'Small and capital letters denote solvent and solute ',
     *           'values, respectively.',
     *      //8x,'i',4x,'Name',2x,'NAT',1x,'Qsolute',1x,'qsolvent',5x,
     *        'E0_i',8x,'EES_i',7x,'EREP_i',6x,'EDI_i',7x,'Esolv_i',8x,
     *       'E_i',/1x,114(1H-))
 9009 format(/1x,'Segment properties (a.u.); ',
     *      /1x,'Small and capital letters denote solvent and solute ',
     *           'values, respectively.',
     *      //8x,'i',4x,'Name',2x,'NAT',1x,'Qsolute',1x,'qsolvent ',3x,
     *      'eps(eff)',/1x,51(1H-))
 9007 format(/1x,'Dipoles (Debye) and solvation energy components ',
     *          '(kcal/mol); ',
     *          'Esolv=Ees+E',A3,'+Edisp+Erep',
     *        //8x,'i subsystem',2x,'option',2x,'solute dipole (x,y,z)',
     *          6x,'solvent dipole (x,y,z)',5x,'Ees_i',3x,'E',A3,'_i',
     *          3x,'Edisp_i',2x,'Erep_i',
     *      /1x,114(1H-))
 9008 format(/1x,'Adjusted dipoles (Debye) with shifted origin,',
     *       /8x,'i ',3x,'solute dipole (x,y,z)',
     *        6x,'solvent dipole (x,y,z)',
     *       /1x,62(1H-))
 9010 format(1x,I8,1x,A8,I4,2F8.4,F14.8,4F12.8,F15.8)
 9011 format(1x,I8,1x,A4,I8,1x,A8,1x,F6.1,F10.6)
 9012 format(1x,I8,2I8,6F9.4,4F9.3)
 9013 format(1x,I8,6F9.4)
 9014 format(1x,'Segment pair properties (kcal/mol); ',
     *          'E=EES+EDI+(Esolv=Ees+Edisp+Erep)',
     *     //8x,'i',7x,'j',4x,'R,A',3x,'R,rel',3x,'EES_ij',
     *     3x,'EDI_ij',3x,'Ees_ij',2x,'Edisp_ij',2x,'Erep_ij',3x,'E_ij',
     *      /1x,85(1H-))
 9015 format(1x,'Segment pair properties (kcal/mol); ',
     *          'E=EES+EREP+EDI+(Esolv=Ees+Edisp+Erep)',
     *     //8x,'i',7x,'j',4x,'R,A',3x,'R,rel',3x,'EES_ij',3x,'EREP_ij',
     *     2x,'EDI_ij',3x,'Ees_ij',2x,'Edisp_ij',2x,'Erep_ij',3x,'E_ij',
     *      /1x,94(1H-))
 9016 format(1x,'Segment pair properties (kcal/mol); ',
     *          'E=EES+EDI+(Esolv=Ees+Enones)',
     *     //8x,'i',7x,'j',4x,'R,A',3x,'R,rel',3x,'EES_ij',
     *     3x,'EDI_ij',3x,'Ees_ij',2x,'Enones_ij',2x,'E_ij',/1x,76(1H-))
 9020 format(1x,2I8,F8.2,F7.2,7F9.3)
 9100 FORMAT(/1X,'          NCC ELECTRON ENERGY (',I1,')=',F19.10/
     *        1X,'SCC CHARGE INTERACTION ENERGY (',I1,')=',F19.10/
     *        1X,'             REPULSION ENERGY (',I1,')=',F19.10/
     *        1X,'            DISPERSION ENERGY (',I1,')=',F19.10/
     *        1X,'               SOLVENT ENERGY (',I1,')=',F19.10/
     *        1X,'                 TOTAL ENERGY (',I1,')=',F19.10/)
 9200 FORMAT(/1X,'ORIGIN FOR DIPOLES=',3F15.8,
     *       /1X,'SOLUTE DIPOLE= ',4F12.6,/1X,'SOLVENT DIPOLE=',4F12.6,
     *       /1X,'TOTAL DIPOLE=  ',4F12.6,' (DEBYE)')
 9210 FORMAT(/1X,'ORIGIN FOR SOLUTE DIPOLES=',3F15.8,
     *       /1X,'ORIGIN FOR SOLVENT DIPOLES=',3F15.8,
     *       /1X,'SOLUTE DIPOLE= ',4F12.6,/1X,'SOLVENT DIPOLE=',4F12.6,
     *       /1X,'TOTAL DIPOLE=  ',4F12.6,' (DEBYE)')
 9300 FORMAT(/1X,29(1H-),/1X,'Subsystem partition analysis',/1X,29(1H-),
     *       /8X,'i',4x,'sys',5x,'Epart(i),h',2x,'Eint(iJ),kcal/mol ',
     *           '(i is segment, J is subsystem)')
 9310 FORMAT(1X,I8,I7,F15.8,9999(8F11.3,/,31X))
 9330 FORMAT(/1X,'Total energy (check):',F20.10,/)
 9400 FORMAT(/1X,'subsystem',4x,'E''(I),h',6x,'Epart(I),h',3x,
     *           'Eint(IJ),kcal/mol (I, J are subsystems)')
 9410 FORMAT(1X,I7,2F15.8,9999(8F11.3,/,31X))
 9500 format(/1x,'This truncated PA ends here.',/)
      END
C*MODULE fmoio   *DECK pdbind
C>
C>     @brief Read PDB to get the indices.
C>
C>     @details Collect indices in PDB.
C>
C>     @author  Dmitri Fedorov
C>
      subroutine pdbind(gprnam,indatp,segnam,natfmo,nsegm,nsegnam,
     *                  modpan,patnam,modseg,nprfmo)
c     IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      IMPLICIT NONE
      character*99 STR
      character*8 segnam(*),at,gprnam
      character*1 n1,chain,chainp
      LOGICAL ressplit,basesplit,GOPARR,DSKWRK,MASWRK,dolend,dosplit,
     *        didsplit,dosplitn,nucleo,amino,autosub,copysub,verbatim,
     *        keepA
      DOUBLE PRECISION patnam(*)
      INTEGER indatp(*),modseg(*),natfmo,nsegm,nsegnam,IR,IW,IP,IJK,
     *        IJKT,IDAF,NAV,IODA,ME,MASTER,NPROC,IBTYP,IPTIM,i,ieof,is,
     *        js,nat,nline,nsegm0,modpan,nsegme,j,nprfmo,nchain
C
      COMMON /IOFILE/ IR,IW,IP,IJK,IJKT,IDAF,NAV,IODA(950)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
c
c     Read INDATP from PDB.
c
      ressplit=iand(modpan,8).ne.0
      basesplit=iand(modpan,128).ne.0
      autosub=iand(modpan,256).ne.0
      copysub=iand(modpan,512).ne.0
      verbatim=iand(modpan,1024).ne.0
      keepA=iand(modpan,2048).ne.0
      if(autosub.and.copysub) call abrtx("modpan=256+512 is not proper")
      if(verbatim.and.(ressplit.or.basesplit)) call abrtx("modpan err.")
c     either one or the other, but not both.
      CALL SEQREW(IR)
      CALL FNDGRP(IR,gprnam,IEOF)
      if(ieof.ne.0) then
        write(iw,*) gprnam,' not found'
        call abrt
      endif
      if(maswrk) then
      if(iand(nprfmo,3).le.1)
     *  write(iw,9020) ressplit,basesplit,autosub,copysub,verbatim,keepA
      nat=0
      nline=0
      nsegm=0
      nsegm0=0
      is=-1
      didsplit=.false.
      chainp='-'
c     Presumably, an impossible chain name.
      nchain=0
c     Multiple chians with side splittings will loose chain om the name.
      do while(.true.)
        nline=nline+1
        READ(IR,9000,END=120,ERR=600) STR
        call UPRCAS(str(1:6),6)
        if(dolend(str)) exit
        if(str(1:6).eq.'ATOM  ' .or.  str(1:6).eq.'HETATM') then
          nat=nat+1
          chain=str(22:22)
          READ(str(23:26),*,END=600,ERR=600) js
          if(chain.ne.chainp) then
            nchain=nchain+1
            chainp=chain
          endif
          if(verbatim) then
            nsegm0=js
            nsegm=max(nsegm,js)
            if(nsegm0.le.nsegnam) then
              segnam(nsegm0)(1:4)=str(17:20)
              segnam(nsegm0)(5:8)=str(23:26)
            endif
          else
          if(is.ne.js.and.nsegm.gt.0) then
            do i=1,nsegm
              if(segnam(i)(4:8).eq.str(22:26)) then
                nsegm0=i
                is=js
                exit
              endif
            enddo
          endif
          if(is.ne.js) then
            nsegm=nsegm+1
            is=js
            if(nsegm.le.nsegnam) then
              if(.not.keepA.and.str(22:22).eq.'A') str(22:22)=' '
c             Remove the chain label for chain A.
              segnam(nsegm)(1:3)=str(18:20)
              segnam(nsegm)(4:8)=str(22:26)
            endif
            nsegm0=nsegm
          endif
          endif
          if(nat.gt.natfmo) then
            write(iw,*) 'Numbers of atoms mismatch in $FMOXYZ and $PDB'
            call abrt
          endif
c         write(6,7777) ,nat,is_save,js_save,nat,nsegm0
c7777     format(1x,'wwwpdb',5I4)
          indatp(nat)=nsegm0
          at=str(13:16)
c         split side chains
          dosplit=ressplit
          dosplitn=basesplit
c         modpan takes precedence over modseg.
          if(nsegm.le.nsegnam.and..not.verbatim) then
            dosplit=iand(modseg(nsegm),1).ne.0.or.ressplit
            dosplitn=iand(modseg(nsegm),2).ne.0.or.basesplit
          endif
          n1=' '
          if(str(18:19).eq.' D') n1=str(20:20)
          if(str(18:18).eq.'D' ) n1=str(19:19)
          if(str(18:19).eq.'  ') n1=str(20:20)
          nucleo=str(1:6).eq.'ATOM  '.and.(n1.eq.'A'.or.n1.eq.'C'.or.
     *           n1.eq.'G' .or.n1.eq.'T'.or.n1.eq.'U')
c         Naivite in full bloom.
          amino=str(1:6).eq.'ATOM  '.and..not.nucleo
c         split side chains in residues
          if(dosplit.and.amino.and.
     *       at.ne.' C  '.and.at.ne.' N  '.and.at.ne.' CA '.and.
     *       at.ne.' O  '.and.at.ne.' H  '.and.at(1:3).ne.' HA') then
            indatp(nat)=-nsegm0
c           write(6,*) 'atom',nat,' is side chain of ',nsegm,dosplit
            didsplit=.true.
          endif
c         split bases in nucleotides
          if(dosplitn.and.nucleo.and.
     *      (at(4:4).eq.' '.and.at(2:2).ne.'P'.or.
     *       at(2:2).eq.'H'.and.at(4:4).ne.'''')) then
            indatp(nat)=-nsegm0
c           write(6,*) 'atom',nat,' is base of ',nsegm,dosplitn
            didsplit=.true.
          endif
          patnam(nat)=transfer(at,patnam(1))
        endif
      enddo
  120 continue
      if(nat.ne.natfmo) then
        write(iw,*) 'Numbers of atoms $FMOXYZ and $PDB',natfmo,nat
        call abrt
      endif
c     write(6,*) 'indatp1=',(indatp(i),i=1,natfmo)
      nsegm0=nsegm
      nsegme=nsegm
      if(didsplit) then
        if(autosub.or.copysub) then
c         sesquipass: count segments.
          js=1
          do i=1,natfmo
            is=indatp(i)
            if(is.lt.0.and.is.ne.js) then
              nsegme=nsegme+1
              js=is
            endif
          enddo
          if(autosub) then
            do i=1,nsegm
              modseg(i+nsegme)=i
            enddo
          endif
          if(copysub) then
c           Copy old data to the new location backwards.
      write(iw,*) 'siz',nsegme,nsegm0
      write(iw,*) 'orig',(modseg(i),i=1,20)
            do i=nsegm,1,-1
              modseg(i+nsegme)=modseg(i+nsegm0)
            enddo
c           clear copied indices for tidiness
      write(iw,*) 'copied',(modseg(i),i=1,40)
            do i=nsegm0+1,nsegme
              modseg(i)=0
            enddo
      write(iw,*) 'zeroed',(modseg(i),i=1,40)
          endif
        endif
c       second pass: split side chains as new segments
        js=1
c       note that "is" is negative, so +1 for js is fine.
        write(iw,*) 'Number of segments before splitting:',nsegm
        do i=1,natfmo
          is=indatp(i)
          if(is.lt.0) then
            if(is.ne.js) then
              nsegm=nsegm+1
              js=is
              if(nsegm.le.nsegnam) then
                segnam(nsegm)=segnam(abs(is))
c               segnam(nsegm)(1:1)='s'
c               The side chain is now reflected in position 4,
c               the same one used for storing chain labels.
c               It is expected that normal chain labels are capitals,
c               and smalls are used for side chains.
                if(segnam(nsegm)(4:4).eq.' ') segnam(nsegm)(4:4)='A'
c               That is, ALA is backbone chain A; ALAa is its side chain.
c               That is, ALAB is backbone chain B; ALAb is its side chain.
                call LWRCAS(segnam(nsegm)(4:4),1)
              endif
              if(autosub) modseg(nsegm+nsegme)=abs(is)
              if(copysub) modseg(nsegm+nsegme)=modseg(abs(is)+nsegme)
            endif
            indatp(i)=nsegm
c           write(6,*) 'Atom',i,-is,nsegm
          endif
        enddo
c     write(iw,*) 'final',(modseg(i),i=1,40)
        write(iw,*) 'Split',nsegm-nsegm0,'side chains/bases.'
      endif
      write(iw,9030) nsegm,nchain,gprnam
c     nsegm.ne.nsegnam is not an error, just a warning.
      if(nsegm.ne.nsegnam) write(iw,9010) nsegnam,nsegm
c     write(iw,*) 'IND',(indatp(i),i=1,natfmo)
c     write(iw,*) 'sub',(modseg(i+nsegme),i=1,nsegm)
      if(autosub) write(iw,*) 'Automatically assigned subsystems',nsegme
      if(autosub.and.nsegme.ne.nsegm) call abrtx("modpan=256 failure")
      endif
c     Squeeze segment names
      if(nsegnam.gt.0) then
c       squeeze spaces
        do nsegm0=1,min(nsegm,nsegnam)
          j=1
          do i=1,8
            if(segnam(nsegm0)(i:i).ne.' ') then
              segnam(nsegm0)(j:j)=segnam(nsegm0)(i:i) 
              j=j+1
            endif
          enddo
          if(j.le.8) then
            do i=j,8
              segnam(nsegm0)(i:i)=' '
            enddo
          endif
        enddo
      endif
      if(goparr) then
        CALL DDI_BCAST(2400,'I',indatp,natfmo,MASTER)
        IF(autosub) CALL DDI_BCAST(2401,'I',modseg,nsegm*2,MASTER)
c       SEGNAM is not broadcast - slaves do not have it.
      endif
      return
  600 continue
      write(6,*) 'error reading ',gprnam,' check entry ',nline
      call abrt
      RETURN
 9000 format(A80)
 9010 format(1x,'Numbers of segments in INDATP(1) and $PDB',2I10)
 9020 format( 1x,'PDB parser options:',
     *       /5x,'Split residues:      ',L1,
     *       /5x,'Split bases:         ',L1,
     *       /5x,'Subsystematization:  ',L1,
     *       /5x,'Inherit subsystems:  ',L1,
     *       /5x,'Verbatim segment ID: ',L1,
     *       /5x,'Keep A chain labels: ',L1)
 9030 format(1x,'Read',I8,' segments in',I3,' chain(s) from ',A8)
      END
C*MODULE fmoio   *DECK writepdb
C>
C>     @brief Write PDB.
C>
C>     @details Write out PDB.
C>
C>     @author  Dmitri Fedorov
C>
      subroutine writepdb(grpnam,nato,c,text,fmosym)
      IMPLICIT NONE
C     Declarations of arguments
      character*8 grpnam
      INTEGER nato
      DOUBLE PRECISION c(3,nato)
      character*72 text
      character*80 fmosym
C     Declarations of common blocks
      INTEGER IR,IW,IP,IJK,IJKT,IDAF,NAV,IODA
      INTEGER ME,MASTER,NPROC,IBTYP,IPTIM
      LOGICAL GOPARR,DSKWRK,MASWRK
C     Other declarations
      DOUBLE PRECISION, PARAMETER :: UNITS=0.52917724924D+00
      INTEGER i,ieof,nat,nline
      logical dolend
      character*99 STR
C
      COMMON /IOFILE/ IR,IW,IP,IJK,IJKT,IDAF,NAV,IODA(950)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
c
c     Read PDB provided by user in the input and write
c     the same on the punch file, replacing the coordinates only.
c     "text" is some comment to add.
c
      CALL SEQREW(IR)
      CALL FNDGRP(IR,' $PDB   ',IEOF)
      if(ieof.ne.0.or..not.maswrk) return
      write(ip,*) grpnam
      write(ip,9200) fmosym(1:72),text
      nat=0
      nline=0
      do while(.true.)
        nline=nline+1
        READ(IR,9000,END=120,ERR=600) STR
        call UPRCAS(str(1:6),6)
        if(dolend(str)) exit
        if(str(1:6).eq.'ATOM  ' .or.  str(1:6).eq.'HETATM') then
          nat=nat+1
          if(nat.gt.nato) then
            write(iw,*) 'Numbers of atoms mismatch',nat,nato
            call abrt
          endif
          write(str(31:54),9100) (c(i,nat)*units,i=1,3)
        endif
        write(ip,9000) str
      enddo
  120 continue
      if(nat.ne.nato) then
        write(iw,*) 'Numbers of atoms differ in WRITEPDB',nato,nat
        call abrt
      endif
      write(ip,*) '$END'
      write(iw,*) 'Successfully wrote PDB to punch for',nat,' atoms.'
      return
  600 continue
      write(6,*) 'error reading ',grpnam,' check entry ',nline
      call abrt
      RETURN
 9000 format(A80)
 9100 format(3F8.3)
 9200 format('REMARK ',A72,/,'REMARK ',A72)
      END
C*MODULE fmoio   *DECK saout
C>
C>     @brief Write SA results. 
C>
C>     @details Print results for SA.
C>
      subroutine saout(mode,nbodym,nbodyp,nfga,nfgb,igeoma,igeomb,
     *                 frgnam,molfrg,samon,sadim,esatot)
C     IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      IMPLICIT NONE
      integer molfrg(*),mode,nfga,nfgb,igeoma,igeomb,IR,IW,IP,IJK,IJKT,
     *        IDAF,NAV,IODA,nfg,i,j,loop,jj,nbodym,nbodyp
      DOUBLE PRECISION frgnam(*),samon(2,*),sadim(*),esatot,AUTOKCAL,
     *                 epld,edesolv,eint,esa
      character*1 sysnam(2),is,js
      character*3 yesno(0:1),monstr
      character*12 unitstr,intstr
      parameter (AUTOKCAL=627.509541D+00)
      COMMON /IOFILE/ IR,IW,IP,IJK,IJKT,IDAF,NAV,IODA(950)
      data sysnam/'A','B'/,yesno/'no','yes'/
C
c     Results for auto subsystem analysis (SA) for fragments.
c     mode=0 fragments
c     mode=1 segments
c
      if(mode.eq.1) then
        unitstr='segments'
        is='i'
        js='j'
        monstr='CTd'
      else
        unitstr='fragments'
        is='I'
        js='J'
        monstr='PLd'
      endif
      intstr='none'
      if(nbodyp.eq.2) intstr='pair'
      if(nbodyp.eq.3) intstr='pair+triple'
      nfg=nfga+nfgb
      write(iw,9000) unitstr,nfg,nbodym,intstr,yesno(igeoma),
     *               yesno(igeomb)
      write(iw,9010) monstr,is,is,is,is,is,is,is
      write(iw,9020) is,monstr
      epld=0
      edesolv=0
      eint=0
      do i=1,nfg
c       samon(1,*) contains Epart. Compute Epld/Ectd here.
        samon(1,i)=samon(1,i)-samon(2,i)
        write(iw,9100) i,sysnam(molfrg(i)),frgnam(i),
     *                 (samon(j,i)*AUTOKCAL,j=1,2)
        epld=epld+samon(1,i)
        edesolv=edesolv+samon(2,i)
        esa=epld+edesolv
      enddo
      if(nbodyp.ge.2) then
        write(iw,9150) is,js
        do i=1,nfga
          do j=1,nfgb
            jj=j+nfga
            loop=(jj*jj-3*jj+2)/2+i
            write(iw,9200) i,sysnam(molfrg(i)),frgnam(i),
     *                     jj,sysnam(molfrg(jj)),frgnam(jj),
     *                     sadim(loop)*AUTOKCAL
            eint=eint+sadim(loop)
          enddo
        enddo
        esa=esa+eint
      endif
      intstr='frag'
      if(mode.eq.1) intstr='seg '
      write(iw,9300) unitstr,
     *                monstr,  intstr,epld*AUTOKCAL,
     *               'DEsolv',intstr,edesolv*AUTOKCAL,
     *               'DEint ',intstr,eint*AUTOKCAL,
     *               '  DE  ',intstr,esa*AUTOKCAL,
     *               '  DE  ',intstr,esatot*AUTOKCAL
      RETURN
 9000 format(/1x,65(1H-),
     *       /1x,'Auto subsystem analysis for A+B->AB (energies are in',
     *           ' kcal/mol)',
     *       /1x,'basic units:      ',A9, 
     *       /1x,'number of units:  ',I7, 
     *       /1x,'unit state:       CT',I1,
     *       /1x,'interaction type: ',A11, 
     *       /1x,'deformation of A: ',A3,
     *       /1x,'deformation of B: ',A3,
     *       /1x,65(1H-)) 
 9010 format(/1x,'Note that DE',A3,'(',A1,')+DEsolv(',A1,')=Epart(',A1,
     *           ',AB)-Epart(',A1,',A/B)',
     *       /1x,'whereas DEsolv(',A1,')=Esolv(',A1,',AB)-Esolv(',
     *           A1,',A/B)')
 9020 format(/7x,A1,' S',3x,'Name',9x,'DE',A2,'d',6x,'DEsolv')
 9100 format(1x,I7,1x,A1,1x,A8,2F12.3)
 9150 format(/7x,A1,' S',3x,'Name',8x,A1,' S',3x,'Name',9x,'DEint')
 9200 format(1x,2(I7,1x,A1,1x,A8),F12.3)
 9300 format(/1x,'Subsystem analysis (SA) results for ',A9,
     *       /5x,'Monomer destab. energy (DE',A3,' ,',A4,') ',F14.3, 
     *       /5x,'   Desolvation  energy (',A6,',',A4,') ',F14.3, 
     *       /5x,'   Interaction  energy (',A6,',',A4,') ',F14.3, 
     *       /5x,'        Binding energy (',A6,',',A4,') ',F14.3, 
     *       /5x,' (Check binding energy (',A6,',',A4,') ',F17.6,' )',/)
      END
