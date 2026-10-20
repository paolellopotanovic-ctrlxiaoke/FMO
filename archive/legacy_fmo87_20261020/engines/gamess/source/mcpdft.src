C*MODULE MPDFT  *DECK PDFTINP
C> @brief CHECK PDFT FUNCTIONAL TYPES
C> @detail It checks functional types among: tPBE, ftPBE, TBLYP, ftBLYP, trevPBE,ftrevPBE
      SUBROUTINE PDFTINP
      USE modmcpdft
      IMPLICIT NONE
      
      CHARACTER*10 METHOD
      INTEGER LSTRNG, LENMTH

      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      INTEGER IR,IW,IP,IS,IPK,IDAF,NAV,IODA


      F_FLAG = .FALSE.
      T_FLAG = .FALSE.
      
      IF(PDFTYP.eq.TPBE) then
        MCPDFTFUNCTIONAL=eTPBE
        T_FLAG = .TRUE.
      else if(PDFTYP.eq.FTPBE) then
        MCPDFTFUNCTIONAL=eFTPBE
        T_FLAG = .TRUE.
        F_FLAG = .TRUE.
      else if(PDFTYP.eq.TBLYP) then
        MCPDFTFUNCTIONAL=eTBLYP
        T_FLAG = .TRUE.
      else if(PDFTYP.eq.FTBLYP) then
        MCPDFTFUNCTIONAL=eFTBLYP
        T_FLAG = .TRUE.
        F_FLAG = .TRUE.
      else if(PDFTYP.eq.TREVPBE) then
        MCPDFTFUNCTIONAL=eTREVPBE
        T_FLAG = .TRUE.
      else if(PDFTYP.eq.FTREVPBE) then
        MCPDFTFUNCTIONAL=eFTREVPBE
        T_FLAG = .TRUE.
        F_FLAG = .TRUE.
      else
        write(IW,9111)
        call abrt
      end if

      WRITE(UNIT=METHOD,FMT='(A8)') PDFTYP
      LENMTH = LSTRNG(METHOD,10)
      write(IW,9000) METHOD(1:LENMTH) 

      RETURN
 
 9000 FORMAT(/10X,26(1H-)/10X,A10,' MC-PDFT CALCULATION'/10X,26(1H-))
 9111 FORMAT(///1X,'* * * * ERROR * * * *'/
     *       /1X,'INVALID MPCDFT FUNCTIONAL SPECIFIED. PLEASE SELECT:',
     *       /1X,'TPBE, FTPBE, TBLYP, FTBLYP, TREVPBE or FTREVPBE',
     *       //1X,'* * * * * * * * * * *'/)
      END



C*MODULE MCPDFT   *DECK MCPDFT
C> @brief      CALCULATE MC-PDFT ENERGY.
C> THIS IS THE MAIN MC-PDFT SUBROUTINE. 
C>    
C> @authors ANDREW SAND, PRACHI SHARMA
C>    SEE PAPERS SUCH AS
C>  J. CHEM. THEORY COMPUT. 10(9), 3669-3680(2014)
C>  J. CHEM. THEORY COMPUT  11(9), 4077-4085(2015)
C>  J. CHEM. PHYS.  146, 034101 (2017)
C>

      SUBROUTINE MCPDFT
      USE modmcpdft
      use mx_limits, only: mxatm,mxrt,mxao,mxnoro,mxgrid,mxfrz
      use comm_MAXROOTS, only: MXSTATE
      use comm_RDM_MCPDFT, only: ld1a,lxt4,nactive,nnz
      use mod_sformas, only: mulscf,mulsfci,kroot,sfoci
C
      IMPLICIT NONE

      LOGICAL GOPARR,DSKWRK,MASWRK,SVDSKW, PACK2E
      LOGICAL DDITRF,DOOOOO,DOVOOO,DOVVOO,DOVOVO,DOVVVO,DOVVVV,
     *        DOCORE,doexch
      LOGICAL NOS,CLABEL,IVANIC,SOME,MCWTS,MCCAN,MUSTCANON
      INTEGER IR,IW,IP,IS,IJKT,IDAF,NAV,IODA
      INTEGER NINTMX,NHEX,NTUPL,INTTYP,IGRDTYP
      INTEGER NFT11,NFT12,NFT13,NFT14,NFT15,NFT16,IDAF20,NEMEMX
      Integer NSTATE,NX
      Double precision  WSTAT_B(MXRT)
      Integer IWTS_B(MXRT)
      DOUBLE PRECISION METHOD,CISTEP,FINALCI,ACURCY,ENGTOL,DAMP
      INTEGER MICIT,NWORD,NORBMC,NOROT,MOFRZ
      INTEGER NPFLG,NOFO,MCFMO,IDIABAT
      LOGICAL CANONC,FCORE,FORS,EKT,LINSER, SG1
      INTEGER NPRINT,ITOL,ICUT,NORMF,NORMP,NOPK
      INTEGER ME,MASTER,NPROC,IBTYP,IPTIM
      DOUBLE PRECISION RUNTYP,EXETYP
      INTEGER NEVALS,NGLEVL,NHLEVL
      INTEGER NWDVAR,MAXFM,MAXSM,LIMFM,LIMSM
      DOUBLE PRECISION DFTTHR,DFTGTHR,SWOFF,SW0,BSLRD
      INTEGER NDFTFG,NRAD,NTHE,NPHI,NRAD0,NTHE0,NPHI0
      INTEGER NANGPT,NANGPT0
      DOUBLE PRECISION SCFTYP,VBTYP,DFTYPE,TDDFTYP,CITYP,CCTYP
      INTEGER MPLEVL,MPCTYP
      INTEGER NLEB, NLEB0
      INTEGER JANS
C
      COMMON /INTFIL/ NINTMX,NHEX,NTUPL,PACK2E,INTTYP,IGRDTYP
      COMMON /CIFILS/ NFT11,NFT12,NFT13,NFT14,NFT15,NFT16,IDAF20,NEMEMX
      DOUBLE PRECISION WSTATE,SPINS,CRIT,PRTTOL,S,SZ
      DOUBLE PRECISION GRPDET,STSYM,GLIST,DWPARM
      INTEGER NFLGDM,IWTS
      INTEGER NCORSV,NCOR,NACT,NORB
      INTEGER NA,NB,K,KST,IROOT,IPURES,MAXW1,NITER
      INTEGER MAXP,NCI,IGPDET,KSTSYM,NFTGCI,IDWEIGH
      DOUBLE PRECISION FSTATE
      INTEGER IFTS

CJI - Add intrnl to detmrp common block
      logical intrnl
      COMMON /DETMRP/ EBAR0(MXRT),EPT0(MXRT),EPT1(MXRT),WPTST(MXRT),
     *                IPTST(MXRT),NPTST,NVAL,NEXT,IFORB,NOS,IVANIC,
     *                EDSHFT,ESOLPT,intrnl
      DOUBLE PRECISION EBAR0,EPT0,EPT1,WPTST,EDSHFT,ESOLPT
      INTEGER IPTST,NPTST,NVAL,NEXT,IFORB

      COMMON /DETWFN/ WSTATE(MXRT),SPINS(MXRT),CRIT,PRTTOL,S,SZ,
     *                GRPDET,STSYM,GLIST,DWPARM,
     *                NFLGDM(MXRT),IWTS(MXRT),NCORSV,NCOR,NACT,NORB,
     *                NA,NB,K,KST,IROOT,IPURES,MAXW1,NITER,MAXP,NCI,
     *                IGPDET,KSTSYM,NFTGCI,IDWEIGH,
     *                fstate(mxrt),ifts(mxrt)
      DOUBLE PRECISION ENUCR,EELCT,ETOT,STOT,SZZ,ECORE,ESCF,EERD,E1,E2
      DOUBLE PRECISION VEN,VEE,EPOT,EKIN,ESTATE,STATN,EDFT,EDISP
      COMMON /ENRGYS/ ENUCR,EELCT,ETOT,STOT,SZZ,ECORE,ESCF,EERD,E1,E2,
     *                VEN,VEE,EPOT,EKIN,ESTATE(MXRT),STATN,EDFT(3),EDISP

      DOUBLE PRECISION X
      COMMON /FMCOM / X(1)
      INTEGER IA
      COMMON /IJPAIR/ IA(MXAO)
      INTEGER NAT,ICH,MUL,NUM,NQMT,NE,MA,MB
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,MA,MB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      DOUBLE PRECISION ZAN, C
      INTEGER IAN
      COMMON /IOFILE/ IR,IW,IP,IS,IJKT,IDAF,NAV,IODA(950)
      COMMON /MCINP / METHOD,CISTEP,FINALCI,ACURCY,ENGTOL,DAMP,
     *                MICIT,NWORD,NORBMC,NOROT(2,MXNORO),MOFRZ(MXFRZ),
     *                NPFLG(10),NOFO,MCFMO,IDIABAT,
     *                CANONC,FCORE,FORS,EKT,LINSER
      COMMON /OUTPUT/ NPRINT,ITOL,ICUT,NORMF,NORMP,NOPK
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /RUNOPT/ RUNTYP,EXETYP,NEVALS,NGLEVL,NHLEVL
      COMMON /MACHIN/ NWDVAR,MAXFM,MAXSM,LIMFM,LIMSM

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
      COMMON /DFGRID/ DFTTHR,DFTGTHR,SWOFF,SW0,BSLRD(137),NDFTFG,
     *                NRAD,NTHE,NPHI,NRAD0,NTHE0,NPHI0,
     *                NANGPT(MXGRID),NANGPT0(MXGRID),SG1,JANS
      COMMON /WFNOPT/ SCFTYP,VBTYP,DFTYPE,TDDFTYP,CITYP,CCTYP,
     *                MPLEVL,MPCTYP
      COMMON /DFLEB/  NLEB(MXGRID),NLEB0(MXGRID)
      INTEGER I,J,L,LMAX,IJ,IJKL,KL,ioff1,ioff2,IPURES_B
      INTEGER L0,L1,L2,L4,LAST,LD1,LIA,LIP,LIX,LIXX,LL1,LL2,LL3,LL4,LWRK
      INTEGER  LX,LXX,NINACT,nprttf
      INTEGER  LMOS,LNOS,LOADFM,LSINT1,LSINT2,M1,M2,M3,M4,LC
      INTEGER LSINT1c,LSINT2c,LIPc,LIXXc,LXXc
      INTEGER  NC,NO,NT,NORBS,NEED,NEED2
      DOUBLE PRECISION :: DET,DETPT,ORMAS
      INTEGER  iTri
C
C  ----- Read input defining the MRMP2 dimensions -----
C  this should not really be necessary, but the number of core orbitals
C  may be incorrect unless we repeat the -detinp- call
C
      DATA DET/8HDET     /
      DATA DETPT/8HDETPT   /
      DATA ORMAS/8HORMAS   /
      iTri(i,j) = Max(i,j)*(Max(i,j)-1)/2 + Min(i,j)
C
C  OPTION FOR ORMAS PERTURBATION THEORY
C
      IF(CISTEP.EQ.ORMAS) THEN
         CALL ORMASPT
         RETURN
      ENDIF
C
C
C      CALL DETINP(NPRINT,DET)

*      CALL PDFINP(NPRINT,DETPT)
      SOME = MASWRK  .AND.  NPRINT.NE.-5
C
      L0 = NQMT
      LL1 = NUM
      LL2 = (LL1*LL1+LL1)/2
      LL3 =  LL1*LL1
      LL4 =  (LL2*LL2+LL2)/2
C
C  If using natural orbitals  then we must exchange these with
C  canonical mos
C
      IF (NOS) THEN
         LL1 = NUM
         LL3 = LL1*LL1
         CALL VALFM(LOADFM)
         LNOS = LOADFM + 1
         LMOS = LNOS   + LL3
         LAST = LMOS   + LL3
         NEED = LAST - LOADFM - 1
         CALL GETFM(NEED)
         CALL DAREAD(IDAF,IODA,X(LMOS),LL3,15,0)
         CALL DAREAD(IDAF,IODA,X(LNOS),LL3,19,0)
         CALL DAWRIT(IDAF,IODA,X(LMOS),LL3,19,0)
         CALL DAWRIT(IDAF,IODA,X(LNOS),LL3,15,0)
         CALL RETFM(NEED)
      ENDIF
C      mxstate = K
      CALL VALFM(LOADFM)

C  Probably a good first step - transform orbitals.

C   note: we will need additional orbs transformed for the gradient.
      NCOR = NCORSV
      NINACT=NCOR
      nc= 0 
      no= NCOR + NACT
      nt= NORB
 
      NORBS=NACT
      NNZ=NCORSV+NACT

      M1= NACT
      M2=(NACT*NACT+NACT)/2
      M3 = M2*M2
      M4 = (M2*M2+M2)/2
C
      NACTIVE=M1

      L1=NNZ
      L2=(NNZ*NNZ+NNZ)/2
      L4=(L2*L2+L2)/2


      nprttf=0
      dditrf=.false.
      DOOOOO=.false.
      DOVOOO=.false.
      DOVVOO=.false.
      DOVOVO=.false.
      DOVVVO=.false.
      DOVVVV=.false.
      DOCORE=.true.
      DOEXCH=.false.
C      call TRFMCX(nprttf,nc,no,nt,.false.,doexch,

      CALL VALFM(LOADFM)
      LSINT1  = LOADFM  + 1
      LSINT1c = LSINT1  + L2
      LSINT2  = LSINT1c + M2
      LSINT2c = LSINT2  + L4
      LIP     = LSINT2c + M4
      LIPc    = LIP  + L2/NWDVAR + 1
      LXX     = LIPc + M2/NWDVAR + 1
      LXXc    = LIP     + NINTMX
      LIXX    = LXXc    + NINTMX
      LIXXc   = LIXX    + NINTMX
      LC      = LIXXc   + NINTMX  
      LWRK    = LC      + LL2      
      LMOS    = LWRK    + LL3    
      LAST    = LMOS    + LL3    
      NEED2   = LAST - LOADFM - 1
      CALL GETFM(NEED2)

C Save MCSCF orbitals and densities to VIB 0
      CALL DAREAD(IDAF,IODA,X(LMOS),LL3,15,0) 
      CALL DAWRIT(IDAF,IODA,X(LMOS),LL3,30,0) 
      CALL DAREAD(IDAF,IODA,X(LC),LL2,16,0)   
      CALL DAWRIT(IDAF,IODA,X(LC),LL2,32,0)   
C
C  Transform integrals for MCPDFT:
      call TRFMCX(nprttf,0,L1,L1,.false.,doexch,
     *            DDITRF,DOOOOO,DOVOOO,DOVVOO,DOVOVO,
     *            DOVVVO,DOVVVV,DOCORE)
      CALL RDCI12(DDITRF,IJKT,X(LSINT1),X(LSINT2),0,L1,L2,L4,
     *            X(LIP),X(LXX),X(LIXX),NINTMX)
      CALL SEQREW(IJKT)
C
!      if(CITYP.NE.RNONE.AND.MULSFCI.GT.0) then
      IF (SFOCI.EQV..TRUE.) THEN
        MXSTATE = K
      ELSE
C
C  Transform integrals for CI step
        call TRFMCX(nprttf,NCOR,L1,L1,.false.,doexch,
     *            DDITRF,DOOOOO,DOVOOO,DOVVOO,DOVOVO,
     *            DOVVVO,DOVVVV,DOCORE)
        CALL RDCI12(DDITRF,IJKT,X(LSINT1c),X(LSINT2c),0,M1,M2,M4,
     *            X(LIPc),X(LXXc),X(LIXXc),NINTMX)
        CALL SEQREW(IJKT)
      END IF
C
C  Back up WSTATE
      WSTAT_B(:) = 0.0D0
      DO 101 NSTATE=1,MXSTATE
        WSTAT_B(NSTATE)=WSTATE(NSTATE)
  101 CONTINUE
      DO 104 NSTATE=1,MXRT
        IWTS_B(NSTATE)=IWTS(NSTATE)
  104 CONTINUE
      IPURES_B = IPURES
 
        
C  Allocate memory needed for 1- and 2-RDM storage
      CALL VALFM(LOADFM)
 
      LIA  = LOADFM + 1
      LXT4 = LIA  + M1
      LX   = LXT4 + M4
      LIX  = LX   + NINTMX
      LD1  = LIX  + NINTMX 
      LD1A = LD1  + L2
      LAST = LD1A + M2
      NEED = LAST - LOADFM - 1
      CALL GETFM(NEED)

         CALL DAREAD(IDAF,IODA,X(LMOS),LL3,15,0)
C  Loop over states
      DO 100 NSTATE=1, MXSTATE


C For some reason, the orbitals on file change through the calculation,
C but they shouldn't.  So, as a dirty fix, they are rewritten to file here.

       CALL DAWRIT(IDAF,IODA,X(LMOS),LL3,15,0)

C Adjust the weights so that density matrices will be computed for the
C specific state of interest.

         WSTATE(:) = 0.0D0
         WSTATE(NSTATE)=1.0D0
         IWTS(:) = 0
         IWTS(1) = nstate
         IPURES=0

      NCOR=NCORSV
C
C     Do not need to call DETFCI+DETDM2 for SF-ORMAS-PDFT
!      if(CITYP.NE.RNONE.AND.MULSFCI.GT.0) then
      IF (SFOCI.EQV..TRUE.) THEN
        CONTINUE
      ELSE
        CALL DETFCI(NPFLG(3),.false.,DDITRF,0)
        CALL DETDM2(NPRINT)
      END IF
C
C This is implemented now only for the ALDET CI option.  We may need to check
C If CISTEP is still what it should be here.

C     IF(CISTEP.EQ.ALDECI) THEN
C         SVDSKW=DSKWRK
C         DSKWRK=.FALSE.
C         CALL DETFCI(NPFLG(3),.TRUE.,DDITRF,0)
C Included the next two calls  in the prev uncommented if loop              
C         CALL DETFCI(NPFLG(3),.false.,DDITRF,0)
C         CALL DETDM2(NPRINT)
C         DSKWRK=SVDSKW
C      ENDIF
C
C
      T_ALPHA = 0.0d0
      T_BETA = 0.0d0


C Calculate the 1-RDM and 2-RDM specific to each state
C Read and store 2-RDM and 1-RDM in memory
C The following reads the 2-RDM:
      Call DCOPY(M4,0.d0,0,X(LXT4),1)
      CALL RDTR3(NFT15,M1,M2,M4,X(LIA),X(LXT4),X(LX),X(LIX))
C  Read 1-RDM (includes core and active orbs):
      Call SQREAD(NFT15,X(LD1),L2)
C  Multiply off-diagonal elements by two:
      Call DSCAL(L2,2.0d0,X(LD1),1)
      do i=1,L1
          ij = iTri(i,i)
          X(LD1-1+ij) = X(LD1-1+ij)*0.5d0
      end do
C Calculate MCSCF conribution to the MC_PDFT energy
      CALL MCPDFT_CALCENER(E1,E2,L1,L2,LD1,LSINT1,LSINT2)
C Reshape to just active orbs:
      ioff1 = 0
      ioff2 = 0
      do i=1,NNZ
        do j=1,i
          if((i>NCORSV).and.(j>NCORSV)) then
            X(LD1A+ioff1)=X(LD1+ioff2)
            ioff1 = ioff1 + 1
          end if
          ioff2 = ioff2 + 1
        end do
      end do
    

C      IF(CISTEP.EQ.ORMAS) THEN
C some other routine for handling the density matrices
C      ENDIF



C Calculate the density functional component of the energy

      TOTELEA = 0.0d0
      TOTELEB = 0.0d0
      TOTELEAt = 0.0d0
      TOTELEBt = 0.0d0
      debcou = 0
 
      CALL RHFCL
C ASSEMBLE FINAL MC-PDFT ENERGY:
C
C   kroot is a label that stores the state picked by $CIDET 
C   in a SF-ORMAS calc. The IF statement breaks 100 Do loop
C   as SF-ORMAS-PDFT calculates one state at a time. 
!      IF(CITYP.NE.RNONE.AND.MULSFCI.GT.0) THEN
      IF (SFOCI.EQV..TRUE.) THEN
        CALL MCPDFT_ENERGY(E1,E2,EPDFT,ENUCR,KROOT,IW)
        EXIT
      END IF
c
      CALL MCPDFT_ENERGY(E1,E2,EPDFT,ENUCR,NSTATE,IW)
c      CALL ABRT
  100 CONTINUE
C  RESTORE WSTATE
      WSTATE(:) = 0.0D0
      DO 102 NSTATE=1,MXRT
        WSTATE(NSTATE)=WSTAT_B(NSTATE)
  102 CONTINUE
      CALL RETFM(NEED2) 
      CALL RETFM(NEED) 
c      CALL TIMIT(1)

C Restore MCSCF orbitals and densities from VIB0 
      CALL DAREAD(IDAF,IODA,X(LMOS),LL3,30,0)
      CALL DAWRIT(IDAF,IODA,X(LMOS),LL3,15,0)
      CALL DAREAD(IDAF,IODA,X(LC),LL2,32,0)  
      CALL DAWRIT(IDAF,IODA,X(LC),LL2,16,0)  

      return
      
c  150 continue
c      Call ABRT
      END

C*MODULE MCPDFT   *DECK MCPDFT_CALCENER
C> @brief      CALCULATE MC-PDFT ENERGY.
C> @detail IT COMPUTES MCSCF CONTRIBUTION TO THE MC-PDFT ENERGY
C>    
      SUBROUTINE MCPDFT_CALCENER(E1,E2,L1,L2,LD1,LSINT1,LSINT2)
      IMPLICIT NONE
C      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON /FMCOM / X(1)
      DOUBLE PRECISION X
      INTEGER iTri
      INTEGER KINDEX
      INTEGER I,J,L
      INTEGER, INTENT(IN)::L1,L2,LD1,LSINT1,LSINT2
      INTEGER:: LMAX,IJ,IJKL,KL,ioff1,ioff2
      DOUBLE PRECISION, INTENT(OUT):: E1,E2
      DOUBLE PRECISION val,fact,DDot
      iTri(i,j) = Max(i,j)*(Max(i,j)-1)/2 + Min(i,j)
C  Calculate the wave function components of the energy
C  First, E1=\sum_{pq} h_{pq}*D_{pq}

      E1 = 0.0d0
      E1=DDot(L2,X(LD1),1,X(LSINT1),1)

      Call DSCAL(L2,0.5d0,X(LD1),1)
      do i=1,L1
          ij = iTri(i,i)
          X(LD1-1+ij) = X(LD1-1+ij)*2.0d0
      end do

C  Second, E2= 0.5 * \sum_{pqrs} h_{pqrs}*D_{pq}D_{rs}
C  (i.e. just the classical coulomb part)
      
      E2 = 0.0d0
      do i=1,L1
        do j=1,i
          ij=iTri(i,j)
          do kindex=1,i
            if(i.eq.kindex) then
              lmax = j
            else 
              lmax = kindex
            endif
            do l=1,lmax
              kl= iTri(kindex,l)
              ijkl=iTri(ij,kl)
              fact=0.5d0
              if(ij.ne.kl) Fact=fact*2
              if(i.ne.j) Fact = fact*2
              if(kindex.ne.l) Fact = fact*2
              val= X(LSINT2-1+ijkl)*
     &        X(LD1-1+ij)*X(LD1-1+kl)*fact
              E2 = E2 + val
            end do
          end do
        end do
      end do
      
      RETURN
      END

C*MODULE MCPDFT   *DECK MCPDFT_ENERGY
C> @brief      CALCULATE MC-PDFT ENERGY.
C> @detail IT CALCULATES AND PRINTS ENERGIES AND OTEHR INGREDIENTS TO THE OUTPUT FILE
C>    
      SUBROUTINE MCPDFT_ENERGY(E1,E2,EPDFT,ENUCR,NSTATE,IW)
      USE MODMCPDFT, ONLY:TOTELEA,TOTELEB,TOTELEAT,TOTELEBT
      use mod_sformas, only: mulsfci
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: E1,E2,EPDFT,ENUCR
      INTEGER,          INTENT(IN) :: IW,NSTATE 
c      write(*,*), "k, nstate", k, nstate
c      DO 1001  I=1,k
      write(IW,9000) NSTATE
      if(MULSFCI.GT.0) then 
        write(IW,9022) E1,E2,ENUCR,E1+E2+ENUCR
        write(IW,9024) EPDFT
        write(IW,9030) NSTATE,ENUCR+E1+E2+EPDFT
      else
        write(IW,9020) E1,E2,ENUCR,E1+E2+ENUCR
        write(IW,9024) EPDFT
        write(IW,9029) NSTATE,ENUCR+E1+E2+EPDFT
      end if

      write(IW,9002)
      write(IW,9021) TOTELEA,TOTELEB,TOTELEAT,TOTELEBT
c 1001 CONTINUE
 9000 FORMAT(/18X,27(1H-)/18X,"ENERGY COMPONENTS, STATE",
     *       1X,I2,/18X,27(1H-)/)
 9001 FORMAT(/1X,65(1H-),1X/)
 9020 FORMAT(/1X,'               ONE ELECTRON ENERGY =',F19.10/
     *        1X,'               TWO ELECTRON ENERGY =',F19.10/
     *        1X,'          NUCLEAR REPULSION ENERGY =',F19.10/
     *        1X,'          TOTAL CASSCF-LIKE ENERGY =',F19.10)
 9022 FORMAT(/1X,'                 ONE ELECTRON ENERGY =',F19.10/
     *        1X,'                 TWO ELECTRON ENERGY =',F19.10/
     *        1X,'            NUCLEAR REPULSION ENERGY =',F19.10/
     *        1X,'          TOTAL SF-ORMAS-LIKE ENERGY =',F19.10)
 9024 FORMAT( 1X,'      PAIR-DENSITY FUNCTIONAL ENERGY =',F19.10) 
 9029 FORMAT(38X,18(1H-)/
     *        1X,'     STATE',1X,I2,1X,'TOTAL MC-PDFT ENERGY ='
     *       ,F19.10)
 9030 FORMAT(40X,18(1H-)/
     *        1X,' STATE',1X,I2,1X,'TOTAL SF-ORMAS-PDFT ENERGY ='
     *       ,F19.10)
 9021 FORMAT(/1X,'   ALPHA DENSITY (BEFORE TRANSLATION) =',F11.5/
     *        1X,'    BETA DENSITY (BEFORE TRANSLATION) =',F11.5/
     *        1X,'    ALPHA DENSITY (AFTER TRANSLATION) =',F11.5/
     *        1X,'     BETA DENSITY (AFTER TRANSLATION) =',F11.5)
      
 9002 FORMAT(/18X,20(1H-)/18X,20HTRANSLATED DENSITIES/18X,20(1H-)/)
      RETURN
      END


