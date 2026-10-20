c 18 Apr 16 - DGF - tweaks for FMO 5.2
C 30 SEP 13 - NDS - EFP1 EMPIRICAL DISPERSION CORRECTION (EFP1-D)
C 18 Feb 13 - MWS - EFGRAD: project gradient if PBC is being used
C 19 Oct 12 - MWS - synchronize FRGINF common
C 21 MAR 12 - TLW,ACW,MWS - REMOVED CALL VNNDER
C                           ADDED CALL EFPDERCOM
C 28 DEC 11 - CS  - EFMO CHANGES FOR COVALENT BONDS
C  2 NOV 11 - MWS - ENABLE FRAGMENT ONLY NUMERICAL GRADIENT PATHWAY
C 11 AUG 11 - AAD - QM/EFP ENERGY DECOMPOSITION ANALYSIS
C 12 DEC 10 - HL  - ADD ARTIFICIAL FORCES ON ATOMS/FRAGMENTS
C  1 OCT 10 - CS  - ADDED EFMO CHANGES
C 11 AUG 10 - DGF - ALLOCATE NEEDED STORAGE FOR PRJGRD
C 25 MAR 10 - HL  - CHANGES FOR CORRECT EFP(NO POL)+PCM JOBS
C 25 MAR 10 - DJS,DGF - PAD COMMON BLOCK PCMPNT, ARGS FOR GRADIENT PROJ.
C  1 MAY 09 - LVS - CHANGES TO ALLOW FOR FORCE POINTS
C  1 MAY 09 - NMT - NEW DIMENSION FOR EFPBAS AND FRAGMT
C 12 JAN 09 - DGF - SYNCHRONISE PCMPNT, AGAIN
C 20 NOV 08 - HL  - POLPCM: REMOVE CHARGE PENETRATION FOR EFP/PCM,
C                   USE AVERAGED INDUCED DIPOLES FOR EFP/PCM
C 23 OCT 08 - LVS - DISPERSION DAMPING, CLEAN UP JIE'S CODE
C 20 AUG 07 - LVS - PBC-EFP CLEANED
C 24 MAR 07 - LVS - PAD FRGINF,EFPPAR; PRINT FUNCTION PREFP ADDED
C 17 JAN 06 - HL  - MD SWITCHING FUNCTIONS, CHANGES TO IND.DIP.,EFP/PCM
C 14 NOV 05 - HL  - POLPCM ADDED TO ALLOW FOR EFP+PCM COMPUTATIONS
C 19 SEP 05 - HMN,IA - CHANGES FOR MOLECULAR DYNAMICS
C 15 JUL 05 - MWS - EFPPOL,EFPPOL2: RESTRICT VCLR RANGES
C 13 JUL 05 - MWS - CHANGE DYN.STORE FOR CHGTRN,EFPPOL
C  5 JUL 05 - MWS - ALLOCATE DYNAMIC MEMORY FOR DIPIT
C 27 JUN 05 - HL  - PERIODIC BOUNDARY COND, PARALLEL EFP2 DISPERSION
C  1 JUN 05 - MWS - EFSP: SUM ALL ES TERMS IN PARALLEL RUN
C 30 APR 05 - DGF - SYNCHRONISE OPTGRD COMMON BLOCK
C 14 MAR 05 - HL  - EFP CHARGE TRANSFER ENERGY, PAULI CHANGES
C  7 MAR 05 - IA  - EFP2 DISPERSION AND DYNAMIC POLARIZABILITY
C 22 FEB 05 - MWS - EFSP: CHANGE ENERGY PRINTING (AND DEBUG OUTPUTS)
C 13 FEB 05 - HMN - ADD DUMMY EFP MD HOOKS INTO EFSP,
C                   ADD EFP PARALLELIZATION INTO EFSP AND EFPPOL
C 13 FEB 05 - JS  - ADD EFP MP2 CALLS TO EFSP AND EFPPOL, ADD EFPPOL2
C 23 JUL 04 - HL  - EFSP: DO FORMAL MEMORY ALLOCATION FOR EFPPOL CALL
C  9 DEC 03 - MWS - EFGRAD: KILL NUMERICAL HESSIAN REQUEST
C  3 JUL 03 - JMM - SUPPRESS PRINTING FOR MONTE CARLO JOBS
C  7 AUG 02 - HT  - USE DYNAMIC MEMORY FOR EFP PAULI REPULSION
C 17 APR 02 - MWS - SYNCH UP FRGINF COMMON
C 13 JUN 01 - MAF - EFSP: COMPUTE CHARGE PENETRATION ENERGY
C 20 FEB 01 - PND - EFSP,EFPPOL: DIP.CONV. AND PRINTING CHANGES
C 25 MAR 00 - MWS - CHANGE ARGUMENTS TO VNNDER
C 20 NOV 97 - MWS - EFSP: USE JAN'S PAULI EXCHANGE ROUTINE
C  2 SEP 97 - MWS - EFSP: CHANGE WORDING OF ENERGY OUTPUT
C 16 MAY 97 - GNM - EFFECTIVE FRAGMENT DRIVER MODULE CREATED
C
C*MODULE EFDRVR  *DECK EFSP
C>
C>     @brief Main EFP energy driver
C>
C>     @details This routine calculates fragment-fragment interaction
C>     energies in the absence of an ab initio molecule
C>
C>     @date 10/5/12 - Spencer Pruitt
C>     - Added EFMOCHTNRG, EFMOREPNRG andD EFMOEPEN to transfer relevant
C>       energies to EFMO code.
C>     @date March, 2013 - Colleen Bertoni
C>     - Cleared LSABSAVE and LSABGRAD so they won't cause
C>       memory errors when different EFMO terms are turned on
C>       and off.
C>     @date January 2016-E.Guidez
C>     -Add printout of E6 and E7 contributions to dispersion energy
C>
C>     @date    September, 2013 - Nuwan De Silva
C>          EFP1 dispersion correction:
C>          The dispersion could be introduced into the EFP1 method
C>          in a manner that is analogous to the Grimme's dispersion
C>          correction to DFT (DFT-D) or HF (HF-D). Such a method
C>          can be called EFP1-D.
C>     @date January 2017-C.Bertoni
C>     -Add changes for EFMO gradient
C>     @date October 2018 - Peng Xu and Tosaporn Sattasathuchana
C>     -Add the call of EFPDYN_OLD to use the old overlap damping function
C>          (IDISP=2)
C>
      SUBROUTINE EFSP
C
      USE MX_LIMITS,ONLY:MXDFG,MXFRG,MXPT,MXFGPT,MXIFRQ,MXDPPT,
     *     mxao,mxatm
      use DYNPOL_DIST
      USE comm_EFPCT
      USE comm_PAULMO
      USE comm_EFPIO
      USE comm_FRGLJP
      USE comm_FRGTYP
      USE comm_EFPFRC
      USE comm_EFPMP2
      USE comm_EFPSCR
      USE comm_FGRAD
      USE comm_FRGINF
      USE comm_POLTMP
      USE comm_PCMEFP
      USE comm_FRGNRG
      USE EFP_LOGICAL,only:E8DISP
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      LOGICAL DC,DCCHG,DCABC
      LOGICAL GOPARR,DSKWRK,MASWRK
      LOGICAL IFGRAD
      LOGICAL IFEWLD,TNFOIL,EWLDPL,DEBUGOPT
C      LOGICAL IFRCPNT, efmo_exrep_grad
      LOGICAL efmo_exrep_grad
C
      PARAMETER (ZERO = 0.0D+00)
      PARAMETER (HALF=0.50D+00)
      PARAMETER (FACT=1.593601268D-03)
C
      DOUBLE PRECISION MD
C
      COMMON /CENTRL/ RQM,XQM,YQM,ZQM,RMM,XMM,YMM,ZMM,ICENFD
      COMMON /DFTDC / DC,IDCVER,DCCHG,DCABC,DCSR,DCS6,DCS8,DCALP,
     *                DCA1,DCA2
      COMMON /EFMO  / EFMOETOT,EFMOESERG,EFMOPOLERG,EFMODISERG,
     *                EFMOREPNRG,EFMOCHTNRG,EFMOEPEN,EFMOPCMG,
     *                IEFMORUN,IMODEFP,IEFMORT,IEFMOCFRG,IEFMONFRG,
     *                IEFMODIM(2),IMODEFE,NATEFMO,IMODEFD,IMODEFER,
     *                IMODEFCT,IDIMTYP,iefmo_agrad
      COMMON /EFMOPO/ LNEFMOPTS,LEFMOPTS,LEFMOPPTS,LEFMODPTS,LNEFMOBAS,
     *                LEFMOBAS,LEFMOLMO,LEFMOFM,LEFMOCV,LEFMOCF,
     *                LEFMODIMG,LEFMOTOTG,LEFMOIPT,LEFMOIGLOB,
     *                LEFMOESDER,MXZVWK,lefmo_fock_der,lefmodipder,
     *     lefmo_tran, lcpl_coeff, liexrep_offset,
     *     lefmo_tot_field, lefmo_f_resp, lefmo_scphf, lefmo_dyndisp,
     *     lefmo_scptdhf
      COMMON /FMCOM / XX(1)
      COMMON /FUNCT / E,EG(3*MXATM)
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON /OPTGRD/ Q(3*MXATM),ES,QQ(20),
     *                CONVF,FMAXT,DXMAXT,RMAX,RMIN,RLIM,
     *                EIGMAX,EIGMIN,GRDERR,FRMS,FMAX,TRMAX,TRMIN,
     *                IC(20),NSTEP,NSERCH,NPMAX,NP,IFOLOW,
     *                NNEG,IUPHSS,IEXIT,ITRUPD,IPAD,KDIAGH,NPRICO
      COMMON /OUTPUT/ NPRINT,ITOL,ICUT,NORMF,NORMP,NOPK
      COMMON /RUNOPT/ RUNTYP,EXETYP,NEVALS,NGLEVL,NHLEVL
      COMMON /SCZVEC/ IDFZVEC,LIPTZVEC,LIPTLG,icursczv
      COMMON /SCZLAG/ LZVLAG,LZVWRK,LYALAG,LYAWRK,LFEQ1
      COMMON /SIMDAT/ NACC,NREJ,IGOMIN,NRPA,IBWM,NACCT,NREJT,NRPAT,
     *                NPRTGO,IDPUNC,IGOFLG
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /PCMDIM/ MXSP,MXTS,MEMPCM1,MEMPCM2,NTS
      COMMON /PCMOPT/ RABI,RASC,REFPOL,THRSLS,DENSLS,WB,WA,ETA2,GD,EVAC,
     *                RHOW,PM,AREATL,AREAKP,BONDRY,OMEGA,RET,FRO,EPSINF,
     *                EPS,DR,RSOLV,VMOL,TCE,STEN,DSTEN,CMF,TABS,IDIRCT,
     *                IPCDER,IDP,ICOMP,IFIELD,ICAV,IDISP,IPRINT,IRETCAV,
     *                ICENT,IFAST,NEVAL,KEFPOL,KEEPSM,IMGABI,IMGASC,NADD
      COMMON /PCMPAR/ IPCM,NFT26,NFT27,IRPPCM,IEF,IP_F,NFMOPCM,IHET
      COMMON /PCMPNT/ LXYZRE,LSSFE,LLIST,LALPCM,LRINPCM,LINAPCM,LINFPCM,
     *                LINIPCM,LMEPCM,LAXYZCT,LXYZCT2,LNVERT,LQSN,LQSE,
     *                LQFS,LQIND,LISPHE,LVAD,LQOR,LVECMUL,LAIPRJ,LFIPRJ,
     *                LPEL,LPCMCDR,LDAI,LIDDAI,LQSND,LQSED,LEPSPCM
      COMMON /PCMPRT/ GCAVP,GCAVS,GDISP,GREP,EHFGAS
      COMMON /EWALD / IFEWLD,TNFOIL,BETA,ELCNST,EPSLN,CUTOFF,CUTLST,
     *                KMAX,NKVEC,LEVEL,EWLDPL
      COMMON /SELFE / ESELFCC, ESELFDD, EEXCLCC, EEXCLCD, EEXCLDD,
     *                EEXCLCQ,
     *                MPLEVL,MPCTYP
      COMMON /MACHSW/ KDIAG,ICORFL,IXDR,modio,mem10,lpnt10,mem10m
C
C
C
      DATA HSSIAN/8HHESSIAN /, SADPT/8HSADPOINT/
      DATA AIRC/8HIRC     /
      DATA MD /8HMD      /
      logical masout
      IPRT=0
      IF(IAND(MODIO,16).NE.0) IPRT=-5
      MASOUT=MASWRK .AND.IPRT.NE.-5
C
C     -------------------------------------------------------------
C     --- THIS ROUTINE CALCULATES FRAGMENT-FRAGMENT INTERACTION ---
C     --- ENERGIES IN THE ABSENCE OF AN AB INITIO MOLECULE      ---
C     -------------------------------------------------------------
C
C     FIRST, INITIALIZE SOME PARALLELIZATION INFORMATION.
C
      CALL EFPPARL
C
C     CCNRG  = CHARGE-CHARGE INTERACTION ENERGY
C     CDNRG  = CHARGE-DIPOLE INTERACTION ENERGY
C     CQNRG  = CHARGE-QUADRUPOLE INTERACTION ENERGY
C     CONRG  = CHARGE-OCTUPOLE INTERACTION ENERGY
C     DDNRG  = DIPOLE-DIPOLE INTERACTION ENERGY
C     DQNRG  = DIPOLE-QUADRUPOLE INTERACTION ENERGY
C     QQNRG  = QUADRUPOLE-QUADRUPOLE INTERACTION ENERGY
C     ESNRG  = TOTAL ELECTROSTATIC INTERACTION ENERGY
C     REPNRG = EXCHANGE-REPULSION ENERGY
C     POLNRG = POLARIZATION ENERGY
C     DISNRG = MP2/EFP1 C6 DISPERSION ENERGY
C     DISNRG8= MP2/EFP1 C8 DISPERSION ENERGY
C     EDISP  = EFP2 DAMPED DISPERSION ENERGY
C     CHTNRG = EFP2 CHARGE TRANSFER ENERGY
C     ENPCM  = PCM ELECTROSTATIC SOLVATION ENERGY
C     GCAVP  = PCM CAVITATION ENERGY
C     GDISP  = PCM DISPERSION ENERGY
C     GREP   = PCM REPULSION ENERGY
C
      E      = ZERO
      CCNRG  = ZERO
      CDNRG  = ZERO
      CQNRG  = ZERO
      CONRG  = ZERO
      DDNRG  = ZERO
      DQNRG  = ZERO
      QQNRG  = ZERO
      ESNRG  = ZERO
      REPNRG = ZERO
      POLNRG = ZERO
      DISNRG = ZERO
      DISNRG8= ZERO
      EDISP  = ZERO
      EDISP6 = ZERO
      EDISP60=ZERO
      EDISP7 = ZERO
      EDISP70 = ZERO
      EDISP0  = ZERO
      CHTNRG = ZERO
      ENPCM  = ZERO
      GCAVP  = ZERO
      GDISP  = ZERO
      GREP   = ZERO
      EPEN   = ZERO
C
C     -- EFP CENTER OF MASS --
C
      CALL EFCM
C
C     -- ELECTROSTATIC ENERGIES --
C     ELECTROSTATIC ENERGIES (AND GRADIENTS) ARE DONE IN FFELEC NOW,
C     REPLACING THE SMALL ARMY OF ROUTINES CHGCHG, ETC. ETC.
C
      MAXDER = 0
      CALL DERCHK(MAXDER)
                                          IFGRAD = .FALSE.
      IF(MAXDER.GT.0)                     IFGRAD = .TRUE.
      IF(MAXDER.EQ.1  .AND.  NGLEVL.EQ.1) IFGRAD=.FALSE.
      IF(MAXDER.EQ.2  .AND.  NHLEVL.EQ.2) IFGRAD=.FALSE.
C
      IF(IEFMORUN.GT.0) THEN
        IF(IEFMORT.NE.3) GOTO 10
        IF(IEFMORT.EQ.3.AND.IDIMTYP.EQ.-1) GOTO 10
      ENDIF
C
      CALL FFELEC(ESNRG,CCNRG,CDNRG,CQNRG,CONRG,DDNRG,DQNRG,QQNRG,
     *                  IFGRAD,.FALSE.)
      if( iefmorun .gt. 0 .and. iefmo_agrad .gt. 0) then
c     if it's a gradient run, we need to make sure the 
c     proper dimensioning is used. it's complicated, since some dimensions
c     for other terms are stored in the CT storage.
         ntctmo_temp = 0
         MXBFEFP = 0
         DO I = 1, NFRG
            MXBFEFP=MAX(MXBFEFP,NPBF(I))
            ntctmo_temp=ntctmo_temp+NCTMO(ISET(I))
         END DO
         LENCV=MXBFEFP*ntctmo_temp
      endif

      if( iefmorun .gt. 0 .and. iefmo_agrad .gt. 0 ) then
         if(maswrk) write(*,*) 'Starting multiple efelec calls'
         call efmo_es_grad( XX(LEFMOESDER), iefmodim, lstmpts, nmpts,
     *        iefmonfrg )

c     use mxzvwk as a flag to tell us whether or not to calculate the response
         if( mxzvwk .gt. 0 ) then
      if(maswrk)  write(*,*) 'Starting multiple response efelec calls'
         call efmo_es_gradz( XX(LZVWRK), iefmodim, lstmpts, nmpts,
     *        iefmonfrg, xx(LZVLAG), xx(LIPTLG),
     *        mxzvwk)

         end if
      endif

 10   CONTINUE
C
      IF (IFEWLD) THEN
C  ---- EWALD SUMS: DIRECT TERM IS CALCULATED WITHIN FFELEC
C  ---- HERE WE CALCULATE OTHER TERMS:
C  ---- SELF-ENERGY & EXCLUSION ENERGY, RECIPROCAL, AND TIN-FOIL PART
C
         ESNRG = ESNRG + ESELFCC + ESELFDD + EEXCLCC + EEXCLCD
     *         + EEXCLDD + EEXCLCQ
C
         CALL RECIPR(ERCPRCC,ERCPRCD,ERCPRDD,ERCPRCQ,ERCPRDQ,
     *                    ERCPRQQ,ERCPRCO,IFGRAD)
C         WRITE(6,988) ERCPRCC,ERCPRCD,ERCPRDD,ERCPRCQ
C  988  FORMAT(/1X,'ERCPRCC = ',F18.14/
C     *        1X,'ERCPRCD = ',F18.14/
C     *        1X,'ERCPRDD = ',F18.14/
C     *        1X,'ERCPRCQ = ',F18.14)
         ESNRG = ESNRG +
     *           ERCPRCC + ERCPRCD + ERCPRDD  + ERCPRCQ
C
         CALL TINFOIL(ETFCC,ETFCD,ETFDD, IFGRAD)
C        WRITE(6,989) ETFCC,ETFCD,ETFDD
C  989  FORMAT(/1X,'ETFCC = ',F18.14/
C     *        1X,'ETFCD = ',F18.14/
C     *        1X,'ETFDD = ',F18.14)
         ESNRG = ESNRG + ETFCC + ETFCD + ETFDD
C         WRITE(6,*) 'TOTAL CHARGE-CHARGE: ',CCNRG + ESELFCC +
C     *               EEXCLCC + ERCPRCC + ETINFOIL
C
C  ---- THIS IS DEBUGGING OPTION WHICH CALCULATES ELECTROSTATIC SUM
C  ---- DIRECTLY. FOR UNCLEAR REASON, WORKS PROPERLY ONLY FOR ATOMS NOT
C  ---- FOR MOLECULES
         DEBUGOPT=.FALSE.
         IF (DEBUGOPT) CALL DIRECTSUM
      END IF
C
C     -- EXCHANGE-REPULSION ENERGY ---
C
C      MEMORY FOR SAB MATRIX: NEED TO KEEP IT FOR DISPERSION DAMPING
C
      LSABSAVE = 0
      LSABGRAD = 0

      efmo_exrep_grad = .False.
      IF (NTMO.GT.0) THEN
         MXMO=0
         MXBF=0
         DO I = 1, NFRG
            MXMO=MAX(MXMO,NORBEF(I))
            MXBF=MAX(MXBF,NPBF(I))
         END DO
C
C     check whether or not to get the efmo exchange repulsion
c     gradient
         efmo_exrep_grad = (IEFMORUN.GT.0.AND.IEFMORT.EQ.3.AND.
     *        IDIMTYP.NE.-1 .and. iefmo_agrad .gt. 0)
         CALL VALFM(LOADFM1)
         LSABSAVE = LOADFM1 + 1
         LSABGRAD = LSABSAVE + NFRG*NFRG*MXMO*MXMO
         IF (IFGRAD) THEN
         LAST1 = LSABGRAD+ NFRG*NFRG*MXMO*MXMO*6
         ELSE
         LAST1 = LSABSAVE + NFRG*NFRG*MXMO*MXMO
         END IF
         if(iefmorun .gt. 0 .and. iefmo_agrad .gt. 0) then
            Lefmo_weight = last1
            last1 = lefmo_weight + NFRG*NFRG*MXMO*MXMO
         endif
         NEEDSAB = LAST1 - LOADFM1 -1
         CALL GETFM(NEEDSAB)
         if(iefmorun .gt.0 .and. iefmo_agrad .gt. 0) then
            call vclr( xx(lefmo_weight), 1,
     *        NFRG*NFRG*MXMO*MXMO )
         endif
      END IF
C
      IF(IREP.EQ.1 .AND. ILJP.EQ.1) THEN
         CALL ENRGLJ(REPNRG)
      ELSE IF(IREP.EQ.1 .AND. ILJP.EQ.0) THEN
         MXRPTS = 0
         DO IFRG = 1,NFRG
            MXRPTS = MAX(NRPTS(IFRG),MXRPTS)
         ENDDO
         MXRPTS = MXRPTS+1
         NDFRG2 = (NDFRG*NDFRG+NDFRG)/2
         LEN = MXRPTS*MXRPTS*NDFRG2
         CALL VALFM(LOADFM)
         LAREL = LOADFM + 1
         LCREL = LAREL + LEN
         LEFLD = LCREL + LEN
         LEFAD = LEFLD + 3*NPTTPT
         LAST  = LEFAD + 3*NPTTPT
         NEED = LAST - LOADFM - 1
         CALL GETFM(NEED)
         CALL DAREAD(IDAF,IODA,XX(LAREL),2*LEN,75,0)
C
         CALL EREPUL(REPNRG,XX(LAREL),XX(LCREL),MXRPTS,NDFRG2)
         CALL RETFM(NEED)
      ELSE IF (NTMO.GT.0) THEN
         MXMO2=(MXMO*MXMO+MXMO)/2
         LENPV=MXBF*NTMO
         LENFM=MXMO2*NFRG
C
C        IF GRADIENT IS NOT NEEDED
         IF(IFGRAD .EQV. .FALSE.) THEN
           CALL VALFM(LOADFM)
           LPROVEC = LOADFM  + 1
           LFOCKMA = LPROVEC + MXBF*NTMO
           LSMAT   = LFOCKMA + MXMO2*NFRG
           LTMAT   = LSMAT   + MXBF*MXBF
           LWRK    = LTMAT   + MXBF*MXBF
           LSIJ    = LWRK    + MXBF
           LTIJ    = LSIJ    + MXMO*MXMO
           LFASQ   = LTIJ    + MXMO*MXMO
           LFBSQ   = LFASQ   + MXMO*MXMO
           LAST    = LFBSQ   + MXMO*MXMO
           NEED=LAST-LOADFM-1
           CALL GETFM(NEED)
           CALL DAREAD(IDAF,IODA,XX(LPROVEC),LENPV,267,0)
           CALL DAREAD(IDAF,IODA,XX(LFOCKMA),LENFM,268,0)
           CALL EXREP(REPNRG,XX(LPROVEC),XX(LFOCKMA),XX(LSMAT),
     *                XX(LTMAT),XX(LWRK),XX(LSIJ),XX(LTIJ),
     *                XX(LFASQ),XX(LFBSQ),XX(LSABSAVE),MXBF,MXMO,EPEN)
           CALL RETFM(NEED)
         ELSE
C          - CALCULATE ENERGY & GRADIENTS TOGETHER
           CALL VALFM(LOADFM)
           LPROVEC = LOADFM  + 1
           LDPVEC  = LPROVEC + MXBF*NTMO
           LFOCKMA = LDPVEC  + MXBF*NTMO*3
           LSMAT   = LFOCKMA + MXMO2*NFRG
           LTMAT   = LSMAT   + MXBF*MXBF
           LWRK    = LTMAT   + MXBF*MXBF
           LSAB    = LWRK    + MXBF
           LTAB    = LSAB    + MXMO*MXMO
           LFASQ   = LTAB    + MXMO*MXMO
           LFBSQ   = LFASQ   + MXMO*MXMO
           LDSINT  = LFBSQ   + MXMO*MXMO
           LDTINT  = LDSINT  + MXBF*MXBF*6
           LDSAB   = LDTINT  + MXBF*MXBF*6
           LDTAB   = LDSAB   + MXMO*MXMO*6
           LWRK2   = LDTAB   + MXMO*MXMO*6
           LGXRFRG = LWRK2   + MXMO*MXMO
           LXRFRG  = LGXRFRG + NFRG*NFRG*6
           lweighted_da = LXRFRG  + NFRG*NFRG*6*4
c    multiply by efmo_exrep_grad so that we only allocate memory
c    for the efmo terms if we're doing an efmo run.
           if( efmo_exrep_grad ) then
              ldnew_density = lweighted_da + (MXBF*MXBF+MXBF)/2
              lctvec =  ldnew_density + MXBF*MXBF*8
              ltran     = lctvec + lencv
              last = ltran + mxmo*mxmo
           else
              last = lweighted_da
           endif
           NEED=LAST-LOADFM-1
           CALL GETFM(NEED)
           CALL DAREAD(IDAF,IODA,XX(LPROVEC),LENPV,267,0)
           CALL DAREAD(IDAF,IODA,XX(LFOCKMA),LENFM,268,0)
           CALL DERVEC(1,XX(LDPVEC),XX(LPROVEC),MXBF,NTMO)
C  DIFFERENT SUBROUTINE IF NEED TO GATHER EFP FORCES AT FORCE POINTS
           IF (IFRCPNT) THEN
           CALL GEXREP_FP(REPNRG,XX(LPROVEC),XX(LFOCKMA),XX(LSMAT),
     *                 XX(LTMAT),XX(LWRK),XX(LSAB),XX(LTAB),
     *                 XX(LFASQ),XX(LFBSQ),XX(LDPVEC),XX(LDSINT),
     *                 XX(LDTINT),XX(LDSAB),XX(LDTAB),XX(LWRK2),
     *                 XX(LGXRFRG),XX(LSABSAVE),XX(LSABGRAD),
     *                 MXBF,MXMO,XX(LXRFRG),EPEN)
           ELSE
           if( .not. efmo_exrep_grad ) then
           CALL GEXREP(REPNRG,XX(LPROVEC),XX(LFOCKMA),XX(LSMAT),
     *                 XX(LTMAT),XX(LWRK),XX(LSAB),XX(LTAB),
     *                 XX(LFASQ),XX(LFBSQ),XX(LDPVEC),XX(LDSINT),
     *                 XX(LDTINT),XX(LDSAB),XX(LDTAB),XX(LWRK2),
     *                 XX(LGXRFRG),XX(LSABSAVE),XX(LSABGRAD),
     *                 MXBF,MXMO,XX(LXRFRG),EPEN)
           endif
           if( efmo_exrep_grad ) then
           CALL DAREAD(IDAF,IODA,XX(LCTVEC),lencv,269,0)
           CALL EFMO_GEXREP(REPNRG,XX(LPROVEC),XX(LFOCKMA),XX(LSMAT),
     *                 XX(LTMAT),XX(LWRK),XX(LSAB),XX(LTAB),
     *                 XX(LFASQ),XX(LFBSQ),XX(LDPVEC),XX(LDSINT),
     *                 XX(LDTINT),XX(LDSAB),XX(LDTAB),XX(LWRK2),
     *                 XX(LGXRFRG),XX(LSABSAVE),XX(LSABGRAD),
     *                 MXBF,MXMO,XX(LXRFRG),EPEN,
     *                 xx(lweighted_da),
     *                 xx(ldnew_density), xx(lctvec),
     *                 xx(LZVLAG), xx(LIPTLG), xx(ltran) )
           end if
           END IF
           CALL RETFM(NEED)
         END IF
C
      END IF
C
C     -- ARTIFICIAL FORCE ON FRAGMENTS --
C        WATER RADIUS = 1.928 A, INPUT RMM=13.072 FOR R = 15.0
C        ONLY ENERGY IS DONE HERE, GRADIENT IS IN DECK VNNDER
C
      IF(ICENFD.EQ.1) THEN
         ECENT = ZERO
         R0    = RMM
         RX    = XMM
         RY    = YMM
         RZ    = ZMM
         POWER = 2.0D+00
         DO IFRG = 1, NFRG
            COMX  = EFCENT(1,IFRG)
            COMY  = EFCENT(2,IFRG)
            COMZ  = EFCENT(3,IFRG)
            XABS  = ABS(COMX)
            YABS  = ABS(COMY)
            ZABS  = ABS(COMZ)
            COMR2 = COMX*COMX + COMY*COMY + COMZ*COMZ
            COMR  = SQRT(COMR2)
            IF(COMR.GT.R0) THEN
               ECENT= ECENT + (COMR-R0)**POWER
            END IF
            IF(XABS.GT.RX) THEN
               ECENT= ECENT + (XABS-RX)**POWER
            END IF
            IF(YABS.GT.RY) THEN
               ECENT= ECENT + (YABS-RY)**POWER
            END IF
            IF(ZABS.GT.RZ) THEN
               ECENT= ECENT + (ZABS-RZ)**POWER
            END IF
         ENDDO
C        - THE BEST IS TO ADD ECENT TO EFP-EFP REPULSION -
         REPNRG = REPNRG + ECENT
      END IF
C
C
C     -- MP2/EFP1 DISPERSION ENERGY --
C
      IF(MP2TYP.EQ.1) THEN
         MXDPTS = 0
         DO IFRG = 1,NFRG
            MXDPTS = MAX(NDPTS(IFRG),MXDPTS)
         ENDDO
         MXDPTS = MXDPTS+1
         NDFRG2 = (NDFRG*NDFRG+NDFRG)/2
         LEN = MXDPTS*MXDPTS*NDFRG2
         DISNRG = ZERO
         DISNRG8= ZERO
C
         CALL VALFM(LOADFM)
         LAREL = LOADFM + 1
         LDREL = LAREL + LEN
         LCREL = LDREL + LEN
         LAST  = LCREL + LEN
         LBREL = LAST  + LEN
         LEREL = LBREL + LEN
         LAST1 = LEREL + LEN
         NEED  = LAST1 - LOADFM - 1
         CALL GETFM(NEED)
         CALL DAREAD(IDAF,IODA,XX(LAREL),6*LEN,81,0)
         CN=6.0D+00
         CALL EDISPR(DISNRG,XX(LAREL),XX(LDREL),XX(LCREL),
     *               MXDPTS,NDFRG2,CN)
         CN=8.0D+00
         CALL EDISPR(DISNRG8,XX(LAST),XX(LBREL),XX(LEREL),
     *               MXDPTS,NDFRG2,CN)
         CALL RETFM(NEED)
      END IF
C
C     -- POLARIZATION ENERGY ---
C     -- PCM SOLVATION ENERGY --
C
      MXPPT = 0
      DO IFRG = 1,NFRG
         MXPPT = MAX(NPPTS(IFRG),MXPPT)
      ENDDO
C
      IF(IPCM.EQ.1) THEN
          NPRTBK=NPRINT
          NPRINT=817
          IF(ICAV.EQ.1) CALL CAVITM
          IF(IDISP.EQ.1) CALL DISRPM
          NPRINT=NPRTBK
C
          CALL VALFM(LOADFM)
          LEFLD = LOADFM + 1
          LEFAD = LEFLD  + 3*NPTTPT
          LAST  = LEFAD  + 3*NPTTPT
          NEED=LAST - LOADFM -1
          CALL GETFM(NEED)
          CALL POLPCM(XX(LEFLD),XX(LEFAD),XX(LQSE),XX(LQSED),
     *                XX(LQSN),XX(LAXYZCT),XX(LAXYZCT+MXTS),
     *                XX(LAXYZCT+MXTS*2),XX(LAXYZCT+MXTS*3),
     *                XX(LXYZRE),XX(LXYZRE+MXSP),XX(LXYZRE+MXSP*2),
     *                XX(LXYZRE+MXSP*3),XX(LISPHE),XX(LLIST),POLNRG)
          CALL RETFM(NEED)
      ELSE IF (IEFP.EQ.1) THEN
          CALL VALFM(LOADFM)
          LEFLD = LOADFM + 1
          LEFAD = LEFLD + 3*NPTTPT
          LAST  = LEFAD + 3*NPTTPT
          if(iefmorun .gt. 0 .and. iefmo_agrad .gt. 0) then
             LAST  = LEFAD + 6*NPTTPT
          endif
          NEED=LAST-LOADFM-1
          CALL GETFM(NEED)

          CALL EFPPOL(XX(LEFLD),XX(LEFAD),POLNRG)

      if( iefmorun .gt. 0 .and. iefmo_agrad .gt. 0 ) then

         if(maswrk) then
            write(*,*) 'Starting multiple frgfld calls', iefmonfrg, nfrg
         endif
         call efmo_pol_grad( XX(LEFMOESDER), lstmpts, nmpts,
     *        iefmonfrg, nppts,iefmodim,iefmort,nfrg )

      endif



          CALL RETFM(NEED)
      END IF
C
C     -- EFP2 CHARGE TRANSFER ENERGY --
C
C EFP2 CHARGE TRANSFER NOT YET DERIVED FOR OPEN SHELL (ROHF) METHOD
C
C CHECK TO SEE IF ANY MULTIPLICITIES ARE GREATER THAN 1
C IF SO, SKIP CT CALCULATIONS ALTOGETHER
C
C THIS SKIPS ANY CT CALCULATION BETWEEN
C CLOSED-SHELL MOLECULES AS WELL
C
      DO I=1,NFRG
       IF (MULMAT(I) .GT. 1)  GOTO 100
      END DO
C
      IF(NTCTMO.GT.0)THEN
         MXAOEFP=0
         MXBFEFP=0
         MXCTMOEFP=0
         DO I = 1, NFRG
            MXAOEFP=MAX(MXAOEFP,NAO(ISET(I)))
            MXBFEFP=MAX(MXBFEFP,NPBF(I))
            MXCTMOEFP=MAX(MXCTMOEFP,NCTMO(ISET(I)))
         END DO
         LENCV=MXBFEFP*NTCTMO
C
C        IF NO GRADIENTS
         IF(IFGRAD .EQV. .FALSE.) THEN
           CALL VALFM(LOADFM)
           LCTVEC  = LOADFM  + 1
           LSMAT   = LCTVEC  + MXBFEFP*NTCTMO
           LTFRG   = LSMAT   + MXBFEFP*MXBFEFP
           LTMAT   = LTFRG   + MXCTMOEFP*MXCTMOEFP*MXDFG
           LWRK    = LTMAT   + MXBFEFP*MXBFEFP
           LSAB    = LWRK    + MXBFEFP
           LTAB    = LSAB    + MXCTMOEFP*MXCTMOEFP
           LVEFP   = LTAB    + MXCTMOEFP*MXCTMOEFP
           LVABEFB = LVEFP   + MXBFEFP*MXBFEFP
           LVABEFA = LVABEFB + MXCTMOEFP*MXCTMOEFP
           LVAAEFB = LVABEFA + MXCTMOEFP*MXCTMOEFP
           LVBBEFA = LVAAEFB + MXCTMOEFP*MXCTMOEFP
           LECTFRG = LVBBEFA + MXCTMOEFP*MXCTMOEFP
           LIDOIT  = LECTFRG + NFRG*NFRG
           LQUAD   = LIDOIT  + NMTTPT
           LAST    = LQUAD   + 6*NMTTPT
           NEED=LAST-LOADFM-1
           CALL GETFM(NEED)
           CALL DAREAD(IDAF,IODA,XX(LCTVEC),LENCV,269,0)
           CALL CHGTRN(CHTNRG,XX(LCTVEC),XX(LSMAT),XX(LTFRG),
     *                 XX(LTMAT),XX(LWRK),XX(LSAB),XX(LTAB),
     *                 XX(LVEFP),XX(LVABEFB),XX(LVABEFA),
     *                 XX(LVAAEFB),XX(LVBBEFA),XX(LECTFRG),
     *                 XX(LIDOIT),XX(LQUAD),MXBFEFP,MXCTMOEFP)
           CALL RETFM(NEED)
         ELSE
C          RUN ENERGY & GRADIENTS TOGETHER
           CALL VALFM(LOADFM)
           LCTVEC  = LOADFM  + 1
           LDCVEC  = LCTVEC  + MXBFEFP*NTCTMO
           LSMAT   = LDCVEC  + MXBFEFP*NTCTMO*3
           LTFRG   = LSMAT   + MXBFEFP*MXBFEFP
           LTMAT   = LTFRG   + MXCTMOEFP*MXCTMOEFP*MXDFG
           LWRK    = LTMAT   + MXBFEFP*MXBFEFP
           LSAB    = LWRK    + MXBFEFP
           LTAB    = LSAB    + MXCTMOEFP*MXCTMOEFP
           LVEFP   = LTAB    + MXCTMOEFP*MXCTMOEFP
           LVABEFB = LVEFP   + MXBFEFP*MXBFEFP
           LVABEFA = LVABEFB + MXCTMOEFP*MXCTMOEFP
           LVAAEFB = LVABEFA + MXCTMOEFP*MXCTMOEFP
           LVBBEFA = LVAAEFB + MXCTMOEFP*MXCTMOEFP
           LDSINT  = LVBBEFA + MXCTMOEFP*MXCTMOEFP
           LDTINT  = LDSINT  + MXBFEFP*MXBFEFP*6
           LDSAB   = LDTINT  + MXBFEFP*MXBFEFP*6
           LDTAB   = LDSAB   + MXCTMOEFP*MXCTMOEFP*6
           LWRK2   = LDTAB   + MXCTMOEFP*MXCTMOEFP*6
           LDABEFB = LWRK2   + MXCTMOEFP*MXCTMOEFP
           LDABEFA = LDABEFB + MXCTMOEFP*MXCTMOEFP*6
           LDAAEFB = LDABEFA + MXCTMOEFP*MXCTMOEFP*6
           LDBBEFA = LDAAEFB + MXCTMOEFP*MXCTMOEFP*6
           LDVEFP  = LDBBEFA + MXCTMOEFP*MXCTMOEFP*6
           LGCTFRG = LDVEFP  + MXBFEFP*MXBFEFP*6
           LIDOIT  = LGCTFRG + NFRG*NFRG*6
           LQUAD   = LIDOIT  + NMTTPT
           LAST    = LQUAD   + 6*NMTTPT
           NEED=LAST-LOADFM-1
           CALL GETFM(NEED)
           CALL DAREAD(IDAF,IODA,XX(LCTVEC),LENCV,269,0)
           CALL DERVEC(2,XX(LDCVEC),XX(LCTVEC),MXBFEFP,NTCTMO)
           CALL GCTANA(CHTNRG,XX(LCTVEC),XX(LSMAT),XX(LTFRG),
     *                 XX(LTMAT),XX(LWRK),XX(LSAB),XX(LTAB),
     *                 XX(LVEFP),XX(LVABEFB),XX(LVABEFA),
     *                 XX(LVAAEFB),XX(LVBBEFA),XX(LDCVEC),
     *                 XX(LDSINT),XX(LDTINT),XX(LDSAB),XX(LDTAB),
     *                 XX(LWRK2),XX(LDABEFB),XX(LDABEFA),
     *                 XX(LDAAEFB),XX(LDBBEFA),XX(LDVEFP),
     *                 XX(LGCTFRG),XX(LIDOIT),XX(LQUAD),MXBFEFP,
     *                 MXCTMOEFP)
           CALL RETFM(NEED)
         END IF
      END IF
C
C --- EFP1 GRIMME'S DISPERSION CORRECTION ---
      IF(DC) THEN
         CALL DFTD3(1,EDISP,DUMMY)
      END IF
C     -- EFP2 DYNAMIC POLARIZATION ENERGY/GRADIENTS ---
C
C DDK
C EFP2 DISPERSION NOT YET DERIVED FOR  OPEN SHELL (ROHF) METHOD
C
C CHECK TO SEE IF ANY MULTIPLICITIES ARE GREATER THAN 1
C IF SO, SKIP DISPERSION CALCULATIONS ALTOGETHER
C
C THIS SKIPS ANY DISPERSION CALCULATION BETWEEN
C CLOSED-SHELL MOLECULES AS WELL
C
      DO I=1,NFRG
       IF (MULMAT(I) .GT. 1)  GOTO 100
      END DO
      IF (NDPTTPT.GT.0) THEN
         CALL VALFM(LOADFM)
         LDYNTR = LOADFM + 1
         LAST   = LDYNTR + NFRG*MXIFRQ*MXPPT
         NEED   = LAST -LOADFM -1
         CALL GETFM(NEED)
c-----PX: if IDISP=2, call the old EFPDYN code ---
c         which uses old overlap damping formula and 
c         does not have the correct 3/pi factor
c         so that one can reproduce old results
c
        IF(IDISPS.EQ.2) THEN
         CALL EFPDYN_OLD(EDISP,EDISP0,XX(LDYNTR),XX(LSABSAVE),
     *        XX(LSABGRAD),MXPPT)
        ELSE
         CALL EFPDYN(EDISP,EDISP0,XX(LDYNTR),XX(LSABSAVE),
     *        XX(LSABGRAD),MXPPT,EDISP6,EDISP7,IFGRAD,
     *        xx(lefmo_weight),ntmo)
          
        ENDIF

         CALL RETFM(NEED)

         if(ntmo.gt.0 .and. iefmorun .gt. 0
     *        .and. iefmo_agrad .gt. 0 ) then
            call valfm(loadfm)
            lprovec = loadfm  + 1
            lsmat  = lprovec + mxbf*ntmo
            ltmat   = lsmat   + mxbf*mxbf
            lwrk   = ltmat  + mxbf*mxbf
            lsab    = lwrk    + mxbf
            ldnew_density = lsab    + mxmo*mxmo
            lctvec =  ldnew_density + mxbf*mxbf*4
            ltran     = lctvec + lencv
            last   = ltran + mxmo*mxmo*3
            need=last-loadfm-1
            call getfm(need)
            call daread(idaf,ioda,xx(lctvec),lencv,269,0)
            call daread(idaf,ioda,xx(lprovec),lenpv,267,0)
            call efmo_disp_damping(xx(lefmo_weight), xx(lprovec),
     *           xx(lsmat),xx(ltmat),xx(lwrk), xx(lsab),
     *           mxbf, mxmo,xx(ldnew_density),xx(lctvec), xx(lzvlag),
     *           xx(liptlg),xx(ltran) )
            call retfm(need)
         endif

      END IF
C
 100  CONTINUE
C
C       RETURN SAB MEMORY HERE
      IF(IREP.NE.1 .AND. NTMO.GT.0) CALL RETFM(NEEDSAB)
C
C     ******************
C     -- TOTAL ENERGY --
C     ******************
C         ESNRG = CCNRG + CDNRG + CQNRG + CONRG + DDNRG + DQNRG + QQNRG
C
      E = ESNRG+REPNRG+POLNRG+EDISP+CHTNRG+DISNRG+DISNRG8
     *  + HALF*ENPCM + FACT*(GCAVP+GDISP+GREP) + EPEN
C
C
      NEVALS = NEVALS + 1
C
C     -- PRINT OUT ENERGIES --
C
C     IF EFMO, TRANSFER ENERGIES TO COMMON AND SKIP PRINTOUT
      IF(IEFMORUN.GT.0) THEN
        EFMOESERG = ESNRG
        EFMOPOLERG = POLNRG
        EFMODISERG = EDISP
        EFMOREPNRG = REPNRG
        EFMOCHTNRG = CHTNRG
        EFMOEPEN = EPEN
C        GOTO 800
        IF(IPCM.EQ.1.AND.IEFMORT.EQ.4) EFMOPCMG=ENPCM*HALF
      ENDIF
      IF(RUNTYP.EQ.MD) GOTO 800
C
C         IF YOU ASSIGN IGOFLG=2 HERE, YOU WILL GET FULL ENERGY OUTPUT
C         FOR ALL GEOMETRIES COMPUTED, GOOD FOR DEBUGGING TERM BY TERM.
C
C      IGOFLG=2
C
      IF(MASout  .AND.  IGOFLG.EQ.2) WRITE(IW,9010)
C
C         NOTE THAT THE TOTAL ENERGY IS INSIDE THE IGOFLG=2 PORTION
C         OF THE ELSE CLAUSE, TO SUPPRESS MONTE CARLO ENERGIES,
C         BUT THE MORE BASIC RUNS ALWAYS PRODUCE TOTAL ENERGY OUTPUT.
C
C     IF((RUNTYP.EQ.OPTIMIZ .AND. NSERCH.NE.0) .OR.
      IF(
     *   (RUNTYP.EQ.SADPT   .AND. NSERCH.NE.0) .OR.
     *   (RUNTYP.EQ.AIRC)                      .OR.
     *   (RUNTYP.EQ.HSSIAN)) THEN
         IF(MASOUT  .AND.  IGOFLG.EQ.2) THEN
            WRITE(IW,9090) ESNRG
            WRITE(IW,9100) REPNRG
            WRITE(IW,9110) POLNRG
C            IF(IEFDP.EQ.1) WRITE(IW,9111) EDISP
            IF(IEFDP.EQ.1) THEN
             IF(IDISPS.EQ.2) THEN
              WRITE(IW,9150) EDISP
             ELSE
              WRITE(IW,9111) EDISP
             ENDIF
            ENDIF
C
C --- EFP1 GRIMME'S DISPERSION CORRECTION ---
            IF(IEFDP.NE.1 .AND. DC) WRITE(IW,9310) EDISP
            IF(NTCTMO.GT.0) WRITE(IW,9112) CHTNRG
            IF(MP2TYP.EQ.1) THEN
              WRITE(IW,9130) DISNRG
              WRITE(IW,9140) DISNRG8
            END IF
            WRITE(IW,9210)
         END IF
         IF(MASOUT) THEN
           IF(IPCM.EQ.1) THEN
              WRITE(IW,9300) E
           ELSE
              WRITE(IW,9200) E
           END IF
         END IF
         IGOFLG=1
      ELSE
         IF(MASOUT.AND.  IGOFLG.EQ.2) THEN
            WRITE(IW,9020) CCNRG
            WRITE(IW,9030) CDNRG
            WRITE(IW,9040) CQNRG
            WRITE(IW,9050) CONRG
            WRITE(IW,9060) DDNRG
            WRITE(IW,9070) DQNRG
            WRITE(IW,9080) QQNRG
            WRITE(IW,9089) EPEN
           WRITE(IW,9220)
C            WRITE(IW,9090) ESNRG
            WRITE(IW,9090) ESNRG+EPEN
            WRITE(IW,9100) REPNRG
            WRITE(IW,9110) POLNRG
            IF(IEFDP.EQ.1) THEN 

             IF(IDISPS.EQ.2) THEN
                WRITE(IW,9150) EDISP
             ELSE

                WRITE(IW,9111) EDISP
C            IF (E7DISP) THEN 
                WRITE(IW,9113) EDISP7
C            ENDIF
                WRITE(IW,9114) EDISP6
!---PX disp8
             IF (E8DISP) THEN
                WRITE(IW,9115) EDISP8
             ENDIF
!----END PX
             END IF !idisps=2 if block
            ENDIF
C --- EFP1 GRIMME'S DISPERSION CORRECTION ---
            IF(IEFDP.NE.1 .AND. DC) WRITE(IW,9310) EDISP
C            IF(IEFDP.EQ.1) WRITE(IW,9113) EDISP0
            IF(NTCTMO.GT.0) WRITE(IW,9112) CHTNRG
            IF(MP2TYP.EQ.1) THEN
              WRITE(IW,9130) DISNRG
              WRITE(IW,9140) DISNRG8
            END IF
            WRITE(IW,9210)
            IF(IPCM.EQ.1) THEN
               WRITE(IW,9300) E
            ELSE
               WRITE(IW,9200) E
            END IF
         END IF
C        IGOFLG=1
         IF(MASOUT) WRITE(IW,*)
         CALL EFDIPOLE
         IF(MASOUT) WRITE(IW,*)
      END IF
C
 800  CONTINUE
!      t_elec = t_elec + (t1_elec - t0_elec)
!      t_pol  = t_pol  + (t1_pol - t0_pol)
!      t_exr  = t_exr  + (t1_exr - t0_exr)
!      t_disp = t_disp + (t1_disp - t0_disp)
!      t_ct   = t_ct   + (t1_ct  - t0_ct)
!
!      if(maswrk) then
!        write(6,900) 'time for Coulomb',t_elec
!        write(6,900) 'time for polarization',t_pol
!        write(6,900) 'time for exrep',t_exr
!        write(6,900) 'time for dispersion',t_disp
!        write(6,900) 'time for charge transfer',t_ct
!      endif
! 900  FORMAT(A,F12.2)
C
      RETURN
C
 9010 FORMAT(/17X,38('-')/
     *       17X,'FRAGMENT-FRAGMENT INTERACTION ENERGIES'/
     *       17X,38('-')/)
 9020 FORMAT(17X,'CHARGE-CHARGE         = ',F14.10)
 9030 FORMAT(17X,'CHARGE-DIPOLE         = ',F14.10)
 9040 FORMAT(17X,'CHARGE-QUADRUPOLE     = ',F14.10)
 9050 FORMAT(17X,'CHARGE-OCTUPOLE       = ',F14.10)
 9060 FORMAT(17X,'DIPOLE-DIPOLE         = ',F14.10)
 9070 FORMAT(17X,'DIPOLE-QUADRUPOLE     = ',F14.10)
 9080 FORMAT(17X,'QUADRUPOLE-QUADRUPOLE = ',F14.10)
 9089 FORMAT(17X,'OVERLAP PEN. ENERGY   = ',F14.10)
 9090 FORMAT(17X,'ELECTROSTATIC ENERGY  = ',F14.10)
C 9091 FORMAT(17X,'ELECTR + PEN ENERGY   = ',F14.10)
 9100 FORMAT(17X,'REPULSION ENERGY      = ',F14.10)
 9110 FORMAT(17X,'POLARIZATION ENERGY   = ',F14.10)
 9111 FORMAT(17X,'TOTAL DISPERSION ENERGY(E6+E7+E8) = ',F14.10)
C 9113 FORMAT(17X,'DISPERSION ZERO ENERGY= ',F14.10)
 9112 FORMAT(17X,'CHARGE TRANSFER ENRGY = ',F14.10)
 9113 FORMAT(17X,'E7 DISPERSION ENERGY     = ',F14.10)
 9114 FORMAT(17X,'E6 DISPERSION ENERGY     = ',F14.10)
 9115 FORMAT(17X,'E8 DISPERSION ENERGY     = ',F14.10)
 9130 FORMAT(17X,'C6 DISPERSION ENERGY  = ',F14.10)
 9140 FORMAT(17X,'C8 DISPERSION ENERGY  = ',F14.10)
 9200 FORMAT(17X,'FINAL EFP ENERGY      = ',F14.10)
 9210 FORMAT(17X,38('-'))
 9220 FORMAT(41X,14('-'))
 9300 FORMAT(17X,'FINAL EFPPCM ENERGY   = ',F14.10)
 9310 FORMAT(17X,'GRIMME''S DISP. ENERGY = ',F14.10)
 9150 FORMAT(17X,'DISPERSION ENERGY     = ',F14.10)
C
      END
C*MODULE EFDRVR  *DECK EFPPOL
C>
C>     @brief EFP fragment-fragment polarization energy
C>
C>     @details This routine calculates fragment-fragment polarization
C>     energies in the absence of an ab initio molecule
C>
C>     @date January 2017-C.Bertoni
C>     -Add changes for EFMO gradient
C>
C>

      SUBROUTINE EFPPOL(EFLD,EFADD,POLNRG)
C
C --------------------------------------------------------------
C --- THIS ROUTINE CALCULATES FRAGMENT-FRAGMENT POLARIZATION ---
C ---   ENERGIES IN THE ABSENCE OF AN AB INITO MOLECULE      ---
C --------------------------------------------------------------
C
      USE MX_LIMITS,ONLY:MXFRG,MXFGPT,MXIFRQ,MXDFG,MXDPPT
      USE comm_EFPPAR
      USE comm_EFPAR
      USE comm_FRGINF
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      PARAMETER (ZERO=0.0D+00)
      PARAMETER (DEBYE=2.541766D+00)
C
C
      LOGICAL MINMEM
      LOGICAL GOPARR,DSKWRK,MASWRK
      LOGICAL IFEWLD,TNFOIL,EWLDPL,DEBUGOPT
C
      COMMON /EFMO  / EFMOETOT,EFMOESERG,EFMOPOLERG,EFMODISERG,
     *                EFMOREPNRG,EFMOCHTNRG,EFMOEPEN,EFMOPCMG,
     *                IEFMORUN,IMODEFP,IEFMORT,IEFMOCFRG,IEFMONFRG,
     *                IEFMODIM(2),IMODEFE,NATEFMO,IMODEFD,IMODEFER,
     *                IMODEFCT,IDIMTYP,iefmo_agrad
      COMMON /FMCOM/  XX(1)
      COMMON /FRFD  / PPEFX(MXFGPT), PPEFY(MXFGPT), PPEFZ(MXFGPT)
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON /INDDIP/ TMUX,TMUY,TMUZ,TMUXD,TMUYD,TMUZD,MINMEM
      COMMON /RUNOPT/ RUNTYP,EXETYP,NEVALS,NGLEVL,NHLEVL
      COMMON /SIMDAT/ NACC,NREJ,IGOMIN,NRPA,IBWM,NACCT,NREJT,NRPAT,
     *                NPRTGO,IDPUNC,IGOFLG
      COMMON /EWALD / IFEWLD,TNFOIL,BETA,ELCNST,EPSLN,CUTOFF,CUTLST,
     *                KMAX,NKVEC,LEVEL,EWLDPL
      COMMON /EWFLD / RCPFLD(3,MXFGPT),RCPIIFLD1(3,MXFGPT),
     *                RCPIIFLD2(3,MXFGPT), TNFLD(3),
     *                TNIIFLD1(3),TNIIFLD2(3),
     *                SELFFLD1(3,MXFGPT),SELFFLD2(3,MXFGPT),
     *                EXCLFLD(3,MXFGPT),
     *                EXCLIIFLD1(3,MXFGPT),EXCLIIFLD2(3,MXFGPT)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /MACHSW/ KDIAG,ICORFL,IXDR,modio,mem10,lpnt10,mem10m
C
c$$$      DIMENSION EFLD(3,NPTTPT),EFADD(3,NPTTPT)
      DIMENSION EFLD(3,NPTTPT),EFADD(6,NPTTPT)
C
      DATA ENERGY  /8HENERGY  /,GRADIENT/8HGRADIENT/
C
C     BE VERY CAREFUL ABOUT THE DATA UNDER COMMON FRFD, INDDIP,
C     AND EFPPAR (DIND AND DINDD), SINCE THEY ARE WRITTEN WHEN EPOL
C     IS CALCULATED!!!
C
      logical masout
      IPRT=0
      IF(IAND(MODIO,16).NE.0) IPRT=-5
      MASOUT=MASWRK .AND.IPRT.NE.-5
      CALL FRGFLD( 0, .FALSE. )
      IF (IFEWLD .AND. EWLDPL) THEN
C   ----  DEBUG: CALCULATES STATIC FIELD EXPLICITELY
         DEBUGOPT=.FALSE.
         IF (DEBUGOPT) CALL DIRECTSUM_POL
C   ---- TINFOIL CONTRIBUTION TO THE FIELD IS CALCULATED IN TNFOIL,
C   ---- TOGETHER WITH ELECTROSTATIC TINFOIL ENERGIES
C   ---- EXCLUSION CONTRIBUTION IS CALCULATED IN EXCL_FIELD AND CALLED
C   ---- ONCE IN EWALDX
C   ---- SELF CONTRIBUTION TO THE FIELD IS 0
C   ---- DIRECT CONTRIBUTION WAS CALCULATED IN FRGFLD
C
C      WRITE(6,*)'FRGFLD1 ',PPEFX(1),PPEFY(1),PPEFZ(1)
C      WRITE(6,*)'FRGFLD2 ',PPEFX(2),PPEFY(2),PPEFZ(2)
C      WRITE(6,*)'FRGFLD3 ',PPEFX(3),PPEFY(3),PPEFZ(3)
C
         CALL RECIPR_FIELD
      END IF
      DO LEFP = 1, NPTTPT
         EFLD(1,LEFP) = PPEFX(LEFP)
         EFLD(2,LEFP) = PPEFY(LEFP)
         EFLD(3,LEFP) = PPEFZ(LEFP)
         IF (IFEWLD .AND. EWLDPL) THEN
            DO I=1,3
               EFLD(I,LEFP) = EFLD(I,LEFP)+RCPFLD(I,LEFP)
     *                      +TNFLD(I)
     *                      +EXCLFLD(I,LEFP)
            END DO
C
C   ---- DIRTY THINGS
            PPEFX(LEFP) = EFLD(1,LEFP)
            PPEFY(LEFP) = EFLD(2,LEFP)
            PPEFZ(LEFP) = EFLD(3,LEFP)
         END IF
      ENDDO
C
C      WRITE(6,*)'TOTAL STATIC FIELD ',(EFLD(I,1),I=1,3),
C     *           (EFLD(I,2),I=1,3),(EFLD(I,3),I=1,3)
C
      CALL VALFM(LOADFM)
      LDIPNEW  = LOADFM   + 1
      LDIPNWD  = LDIPNEW  + NFRG
      LDINDC   = LDIPNWD  + NFRG
      LDINDDC  = LDINDC   + 3*NPTTPT
      LDIND_W  = LDINDDC  + 3*NPTTPT
      LDIND_C  = LDIND_W  + 6*NPTTPT
      LEWLD1   = LDIND_C  + 6*NPTTPT
      LEWLD2   = LEWLD1   + 3*NPTTPT
      LAST     = LEWLD2   + 3*NPTTPT
      NEED = LAST - LOADFM - 1
      CALL GETFM(NEED)
      CALL DIPIT(EFLD,EFADD,XX(LDINDC),XX(LDINDDC),
     *           XX(LDIND_W),XX(LDIND_C),
     *           NPTTPT,XX(LDIPNEW),XX(LDIPNWD),
     *           XX(LEWLD1),XX(LEWLD2),NFRG)
      CALL RETFM(NEED)
      IF(IDPUNC.EQ.1) RETURN
C
C     --- CALCULATE POLARIZATION ENERGY ---
C
      IUP=0
      DO I=1,NUP
         IUP = IUP + NPPTS(I)
      ENDDO
      ILOW=1
      DO I=1,NLOW-1
         ILOW = ILOW + NPPTS(I)
      ENDDO
      POLNRG = ZERO
      DO II=ILOW,IUP
         POLNRG = POLNRG + DIND(1,II)*PPEFX(II)
     *                   + DIND(2,II)*PPEFY(II)
     *                   + DIND(3,II)*PPEFZ(II)
      ENDDO
      POLNRG = -0.5D+00*POLNRG
      IF(GOPARR) CALL DDI_GSUMF(677,POLNRG,1)
C
C     --- GRAND INDUCED DIPOLE ---
C
      IF( IEFMORUN.GT.0 .AND. IEFMORT.NE.4 ) GOTO 800
      IF (RUNTYP.EQ.ENERGY .OR. RUNTYP.EQ.GRADIENT) THEN
         TMUX = ZERO
         TMUY = ZERO
         TMUZ = ZERO
         DO LEFP=1, NPTTPT
           TMUX = TMUX + DIND(1,LEFP)
           TMUY = TMUY + DIND(2,LEFP)
           TMUZ = TMUZ + DIND(3,LEFP)
         ENDDO
         TMUXD = TMUX*DEBYE
         TMUYD = TMUY*DEBYE
         TMUZD = TMUZ*DEBYE
         IF(MASOUT)
     *   WRITE(IW,9000) TMUX, TMUXD, TMUY, TMUYD, TMUZ, TMUZD
      ENDIF
C
  800 CONTINUE
C
      RETURN
C
 9000 FORMAT(/,17X,'INDUCED DIPOLE',4X,'ATOMIC UNITS',3X,'DEBYE',/
     *       17X,14('-'),4X,12('-'),3X,5('-')/
     *       21X,'X(IND)',9X,F8.5,3X,F8.5/
     *       21X,'Y(IND)',9X,F8.5,3X,F8.5/
     *       21X,'Z(IND)',9X,F8.5,3X,F8.5)
C
      END
C
C*MODULE EFDRVR  *DECK EFGRAD
C>
C> @brief   compute gradient for EFP-only runs
C>
C> @details compute gradient when only EFP particles are present,
C>          possibly enclosed in a PCM layer.
C>
C> @author  perhaps Paul Day
C>
C
      SUBROUTINE EFGRAD
C
C ---------------------------------------------------------------
C --- THIS ROUTINE CALCULATES THE GRADIENT ASSOCIATED WITH    ---
C --- A FRAGMENT ONLY RUN, I.E., FRAGMENTS IN THE ABSENCE OF  ---
C ---                  AN AB INITIO MOLECULE                  ---
C ---------------------------------------------------------------
C
      USE MX_LIMITS,ONLY:MXPT,MXFRG,MXFGPT,MXDFG,MXDPPT,MXSHEF,MXGEFP,
     *     mxatm
      USE comm_PAULMO
      USE comm_EFPPBC
      USE comm_EFPBAS
      USE comm_EFPMP2
      USE comm_FGRAD
      USE comm_FRGINF
      USE comm_POLTMP
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
c
      PARAMETER (ONE=1.0D+00)
C
      LOGICAL PROJGRAD,PBC
      LOGICAL GOPARR,DSKWRK,MASWRK
C
      DOUBLE PRECISION MD
C
      COMMON /EFMO  / EFMOETOT,EFMOESERG,EFMOPOLERG,EFMODISERG,
     *                EFMOREPNRG,EFMOCHTNRG,EFMOEPEN,EFMOPCMG,
     *                IEFMORUN,IMODEFP,IEFMORT,IEFMOCFRG,IEFMONFRG,
     *                IEFMODIM(2),IMODEFE,NATEFMO,IMODEFD,IMODEFER,
     *                IMODEFCT,IDIMTYP,iefmo_agrad
      COMMON /FMCOM / X(1)
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /PCMPAR/ IPCM,NFT26,NFT27,IKREP,IEF,IP_F,NFMOPCM,IHET
      COMMON /RUNOPT/ RUNTYP,EXETYP,NEVALS,NGLEVL,NHLEVL
C
      DATA MD/8HMD      /
C
      IF( IEFMORUN.GT.0 ) RETURN
C
C          NUMERICAL GRADIENT IS DONE ELSEWHERE...
C
      IF(NGLEVL.NE.0) THEN
         CALL NUMGRDX
         RETURN
      END IF
C
C          THE REST OF THIS ROUTINE DRIVES ANALYTICAL GRADIENTS.
C
      IF(MP2TYP.EQ.1) THEN
        NFRPTS=NMTTPT+NPTTPT+NRTTPT+NTPATM+NTMO+NDTTPT
      ELSE
        NFRPTS=NMTTPT+NPTTPT+NRTTPT+NTPATM+NTMO
      END IF
C
      CALL VALFM(LOADFM)
      LEF3    = LOADFM  + 1
      LWORK   = LEF3    + MAX(1,3*NFRPTS)
      LATORQN = LWORK   + MAX(1,10*NMTTPT)
      LAST    = LATORQN + MAX(1,3*NFRG)
      NEED = LAST - LOADFM - 1
      CALL GETFM(NEED)
C
      CALL VCLR(DEF,1,3*MXFGPT)
      CALL VCLR(ATORQ,1,3*NFRG)
      CALL EFPDERCOM(DEF,X(LEF3),NFRPTS,X(LWORK),X(LATORQN))
C     (NEXT CALL WAS PULLED OUT OF OLD COLLECTIVE VNNDER+EFP ROUTINE)
      IF(GOPARR) CALL DSCAL(3*NAT,ONE/NPROC,DEF,1)
      CALL EFFT
C
C       AS A BAND-AID, PROJECT RESIDUAL TORQUE/DRIFT FROM THE GRADIENT.
C       THIS EXISTS WHEN USING A FINITE PCM GRID FOR EFP+PCM,
C       AND SEEMINGLY WHEN PERIODIC BOUNDARY CONDITIONS EXIST.
C       THE THINKING IS THAT MD RUNS ARE SUPPOSED TO CONSERVE ENERGY,
C       SO THEY SHOULD AVOID ANY MANIPULATION OF THE GRADIENT HERE.
C
C       NOTE THAT NAT=0, SO THERE IS A DUMMY QM GRADIENT VECTOR,
C       WITH NO QM COORDINATES OR QM NUCLEAR MASSES.
C
      PBC = (XBOX+YBOX+ZBOX) .GT. 0.0D+00
      PROJGRAD = (IPCM.EQ.1 .OR.  PBC)  .AND. RUNTYP.NE.MD
      IF(PROJGRAD) THEN
         NCCF = 3*NAT + 6*NFRG
         CALL VALFM(LOADFM)
         LGTOT  = LOADFM + 1
         LRM    = LGTOT  + 3*NAT +  6*NFRG
         LX     = LRM    + 3*NAT + 21*NFRG
         LP     = LX     + 3*NAT +  6*NFRG
         LAST   = LP     + 3*NAT +  6*NFRG
         NEEDP  = LAST-LOADFM-1
         CALL GETFM(NEEDP)
         DUMMY =0.0D+00
         DUMMYC=0.0D+00
         DUMMYZ=0.0D+00
C---     CALL RAMS(X(LRM),1)
         CALL PRJGRD(DUMMY,DUMMYC,X(LX),X(LGTOT),X(LRM),X(LP),DUMMYZ,
     *               NAT,NAT*3,NCCF,.FALSE.)
         CALL RETFM(NEEDP)
      END IF
C
      CALL RETFM(NEED)
      RETURN
      END

C>
C*MODULE EFDRVR  *DECK EFPDYN
C>    
C>  @brief    EFP-EFP dispersion interaction and gradient
C>
C>  @author   Ivana Adamovic (energy), Lyudmila V. Slipchenko (gradient) ? 
C> 
C>  @date     1995, 2006?
C>            July 2014 -Peng Xu-
C>           - add EFP-EFP R^-7 dispersion energy
C>           - isotropic distributed E6 is corrected 
C>             with the factor of 3/pi
C>            January 2016 -E.Guidez-
C>           - Add E7 gradient terms and changed the overlap formula
C>             for F6 and F7
C> @date January 2017 -C.Bertoni
C>           -Added changes for EFMO gradient
C>
C> @param efmo_weight : Stores weight needed for dispersion damping gradient
C> @date Mar 2020 - P. Xu
C>           -Add EFP2-EFP2 R^-8 dispersion energy (disp8)
C>           -rename EDISD,EDISD0,EDISD6 to EDISP,EDISP0,EDISP6
C>           -modularize the subroutine
C> 
C> @param EDISP: damped total dispersion energy 
C> @param EDISP0: undamped total dispersion energy
C> @param EDISP6: damped R^-6 dispersion energy
C> @param EDISP7: damped R^-7 dispersion energy
C> @param IFGRAD: logical to determine if compute gradients
C> @param DYNTR: trace of dynamic dipole polarizability tensor
C>               (xx+yy+zz)
C> @param SABSAVE: overlap between LMOs of fragments A and B
C> @param SABGRAD: overlap derivative
C> @param MXPPTS: maximum number of polarizable points (?) 
C> @param efmo_weight : Stores weight needed for dispersion damping
C>                      gradient
C> @param ntmo: Total number of LMOs of all fragments in the system.
C>              This is 0 if exrep is off.
C>              If exrep is off, the overlap matrix was not computed,
C>              so overlap-based damping can not be carried out
C>
      SUBROUTINE EFPDYN(EDISP,EDISP0,DYNTR,SABSAVE,SABGRAD
     *     ,MXPPTS,EDISP6,EDISP7,IFGRAD, efmo_weight, ntmo)

      use mx_limits, only: mxatm,MXFRG,MXFGPT,MXPT,MXIFRQ,MXDFG,
     *     MXDPPT,MXPAIRS,MXGEFP,MXSHEF,
     *     MXCPUEFP

      USE comm_EFPPBC
      USE comm_EFPSCR
      USE comm_EFPPAR
      USE comm_EFPGES
      USE comm_EFPAR

      use EFP_LOGICAL
      use DYNPOL_MOL
      use DYNPOL_DIST

      USE comm_EFPFRC
      USE comm_FGRAD
      USE comm_FRGINF
      USE comm_EFPGES
      USE comm_FRGNRG, only:EDISP8,EDISP8_aniso
      USE constants

      USE EFAIDIS, only: IQMDISP

      IMPLICIT NONE
C
C      PARAMETER (B=1.5D+00)
      double precision,PARAMETER :: FACT1= 0.72086099022968040154D-02
      double precision,PARAMETER :: FACT2= 0.17697067815034886394D-01
      double precision,PARAMETER :: FACT3= 0.30660908596251749739D-01
      double precision,PARAMETER :: FACT4= 0.48381293256249884995D-01
      double precision,PARAMETER :: FACT5= 0.74878830420650517080D-01
      double precision,PARAMETER :: FACT6= 0.11806515901361630228D+00
      double precision,PARAMETER :: FACT7= 0.19535413832209084204D+00
      double precision,PARAMETER :: FACT8= 0.35055692324483221824D+00
      double precision,PARAMETER :: FACT9= 0.71577113554429568336D+00
      double precision,PARAMETER :: FACT10=1.81409759976323969729D+00
      double precision,PARAMETER :: FACT11=6.97923445114870823247D+00
      double precision,PARAMETER :: FACT12=83.2480938829658453917D+00
C
      LOGICAL :: GOPARR,DSKWRK,MASWRK
      LOGICAL :: IFGRAD
      logical :: SHIFT
      integer :: KRONBC,KRONIJ,KRONIK,KRONJK
      integer :: KRONICD,KRONICE,KRONIDE
      integer :: KRON1A,KRON1B,KRON1C,KRON2D,KRON2E,KRON2F,KRON3G
      integer :: KRON3H,KRON3I,NATOM_IM,NATOM_JM
      integer :: KRONA,KRONB,KRONC,KROND,KRONE,KRONF
      Double precision :: delta

      integer :: NEVALS,NGLEVL,NHLEVL
      double precision :: RUNTYP,EXETYP
      COMMON /RUNOPT/ RUNTYP,EXETYP,NEVALS,NGLEVL,NHLEVL
C
      double precision :: EFMOETOT,EFMOESERG,EFMOPOLERG,EFMODISERG,
     *                    EFMOREPNRG,EFMOCHTNRG,EFMOEPEN,EFMOPCMG
      integer :: IEFMORUN,IMODEFP,IEFMORT,IEFMOCFRG,IEFMONFRG,
     *           IEFMODIM,IMODEFE,NATEFMO,IMODEFD,IMODEFER,
     *           IMODEFCT,IDIMTYP,iefmo_agrad
      COMMON /EFMO  / EFMOETOT,EFMOESERG,EFMOPOLERG,EFMODISERG,
     *                EFMOREPNRG,EFMOCHTNRG,EFMOEPEN,EFMOPCMG,
     *                IEFMORUN,IMODEFP,IEFMORT,IEFMOCFRG,IEFMONFRG,
     *                IEFMODIM(2),IMODEFE,NATEFMO,IMODEFD,IMODEFER,
     *                IMODEFCT,IDIMTYP,iefmo_agrad
      integer :: IR,IW,IP,IS,IPK,IDAF,NAV,IODA
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      integer :: ME,MASTER,NPROC,IBTYP,IPTIM
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
C
      double precision :: DYNTR(MXIFRQ*MXPPTS,NFRG),
     *           SABSAVE(NFRG,NFRG,MXPPTS,MXPPTS),
     *           SABGRAD(NFRG,NFRG,MXPPTS,MXPPTS,6),DSAB(6),
     *           IFP1(10),IFP2(10),
     *     efmo_weight(nfrg,mxppts,nfrg,mxppts)
      double precision :: FACT(12),COORD(3),CRD(3),
     *           Rprime(3,MXFGPT),
     *           DYNDD_ave(12,NFRG),DUM5(3),DUM6(3),DUM7(3),DUM8(3),
     *           DER_EPAIR70(3),DER_EPAIR7(3),
     *           DER_EPAIR7_IJ(3),
     *           DER_EPAIR(6),FORDIS(3,MXFRG),TORDIS(3,MXFRG),
     *           DUMA(3),DUM3A(3),DER_T2X(9),DER_T3X(27),DER_T2Y(9),
     *           DER_T3Y(27),DER_T2Z(9),DER_T3Z(27),
     *           DER_T2TX(9),DER_T2TY(9),DER_T2TZ(9),
     *           DER_T3TX(27),DER_T3TY(27),DER_T3TZ(27),TORDIS6(3,MXFRG)
     *           ,TORDIS7(3,MXFRG),FORDIS6(3,MXFRG),FORDIS7(3,MXFRG),
     *           DER_EPAIR70_IJ(3),
     *           DER_EPAIR_FP_6(3),
     *           DER_T2TX_IJ(9), DER_T2TY_IJ(9), DER_T2TZ_IJ(9),
     *           DER_T3TX_IJ(27), DER_T3TY_IJ(27), DER_T3TZ_IJ(27),
     *           DER_EPAIR_FP_7(3)

      double precision :: TT2(3,3),TT3(3,3,3),TT4(3,3,3,3)
      double precision :: DYNQQ_AVE(MXIFRQ,NFRG)
      double precision :: DYNQQ_LMO_AVE(MXIFRQ*MXPPTS,NFRG)
      double precision :: T2(3,3),T3(3,3,3),
     *          T4(3,3,3,3),T5(3,3,3,3,3)
      double precision :: EDISP,EDISP0,EDISP60,EDISP6,
     *                    EDISP60_aniso,EDISP6_aniso,EDISP70,EDISP7,
     *                    EDISP80,EDISP80_aniso,C6,C8,
     *                    E8DQ,E8DO,E8QQ
      double precision :: EPAIR0,EPAIR,EPAIR60,EPAIR6,EPAIR60_aniso,
     *                    EPAIR6_aniso,EPAIR70,EPAIR7,EPAIR80,EPAIR8,
     *                    EPAIR80_aniso,EPAIR8_aniso
      double precision :: ALP1,ALP2,ALP3,ALP4,ALP5,ALP6,ALP7,ALP8,ALP9,
     *                    ALP10,ALP11,ALP12,ALPN1,ALPN2,ALPN3,ALPN4,
     *                    ALPN5,ALPN6,ALPN7,ALPN8,ALPN9,ALPN10,ALPN11,
     *                    ALPN12
      double precision :: COMX,COMY,COMZ,COMX1,COMX2,COMY1,COMY2,
     *                    COMZ1,COMZ2,X,Y,Z,X1,X2,Y1,Y2,Z1,Z2,
     *                    X12,Y12,Z12,R12,PIL,R12SQ,
     *                    SWF,DSWF,SWFDX,SWFDY,SWFDZ,
     *                    VUMY1,VUMY2,VUMY3,VUMY4,VUMY5,
     *                    PBCX,PBCY,PBCZ,
     *                    DUMY,DUMY1,DUMY2,DUMY3,DUMY4,DUMY5,
     *                    R,R2,R3,R4,R5,R6,R7,R8,R9,R10,
     *                    conx,cony,conz,CUT,F6,DF6,F7,DF7A,DF7B,
     *                    F8,DF8A,DF8B
      double precision :: SAB
      integer :: I,J,K,L,Ifrag,Ifreq,IM,IN,JM,JP,ID,IP1,IP2,IL,
     *           IFP1_COM,IFP2_COM,NEQUAL1,NEQUAL2
      integer :: NTMO,MXPPTS,INPPTS,JNPPTS

      IF (IFRCPNT .AND. E7DISP) THEN
      IF (MASWRK) WRITE (IW,*)
     *  'Force points cannot yet be calculated with E7'
       CALL ABRT
      ENDIF
C     ------------------------------------------------------------------
C        THIS ROUTINE CALCULATES FRAGMENT-FRAGMENT DISPERSION (EFP2)
C        ENERGY AND GRADIENTS IN THE ABSENCE OF AN AB INITO MOLECULE
C     ------------------------------------------------------------------
C
!  initialize vector FACT that stores the 12 point quadrature
      CALL VCLR(FACT,1,12)
C   E6 -- anisotropic molecular E6
C   E7 -- anisotropic molecular E7
C   E6_iso -- isotropic molecular E6
C   E60_aniso -- anisotropic distributed E6(LMO) NO damping
C   E6_aniso -- anisotropic distributed E6(LMO) damped
C   EDISP70  -- anisotropic distributed E7(LMO) NO damping
C   EDISP7   -- anisotropic distributed E7(LMO) damped
C   EDISP  -- total dispersion with damping
C   EDISP0 -- total dispersion No damping  
C   EDISP60  -- isotropic distributed E6(LMO) NO damping
C   EDISP6 -- isotropic distributed E6(LMO) damped

      IP1=0
      IF (GOPARR) THEN
         DO I=1,NLOW-1
            IP1 = IP1 + NPPTS(I)
         ENDDO
      ENDIF

      EDISP7 = ZERO 
      EDISP70 = ZERO
      EDISP = ZERO 
      EDISP0 = ZERO
      EPAIR70=ZERO
      EPAIR7=ZERO
      DF6=ZERO
              CALL VCLR(FACT,1,12)
              CALL VCLR(DSAB,1,6)
              CALL VCLR(T2,1,9)
              CALL VCLR(T3,1,27)
              CALL VCLR(DER_T2X,1,9)
              CALL VCLR(DER_T2Y,1,9)
              CALL VCLR(DER_T2Z,1,9)
              CALL VCLR(DER_T3X,1,27)
              CALL VCLR(DER_T3Y,1,27)
              CALL VCLR(DER_T3Z,1,27)
              CALL VCLR(DER_EPAIR70,1,3)
              CALL VCLR(DER_EPAIR7,1,3)
              CALL VCLR(DER_EPAIR7_IJ,1,3)
              CALL VCLR(DER_EPAIR70_IJ,1,3)
              CALL VCLR(DER_T2TX_IJ,1,9)
              CALL VCLR(DER_T2TY_IJ,1,9)
              CALL VCLR(DER_T2TZ_IJ,1,9)
              CALL VCLR(DER_T3TX_IJ,1,27)
              CALL VCLR(DER_T3TY_IJ,1,27)
              CALL VCLR(DER_T3TZ_IJ,1,27)


      IF(E7DISP.or.E8DISP.or.MOLPOL) THEN
       FACT(1) = FACT1
       FACT(2) = FACT2
       FACT(3) = FACT3
       FACT(4) = FACT4
       FACT(5) = FACT5
       FACT(6) = FACT6
       FACT(7) = FACT7
       FACT(8) = FACT8
       FACT(9) = FACT9
       FACT(10)= FACT10
       FACT(11)= FACT11
       FACT(12)= FACT12
      ELSE
       FACT(1) = ZERO
       FACT(2) = ZERO
       FACT(3) = ZERO
       FACT(4) = ZERO
       FACT(5) = ZERO
       FACT(6) = ZERO
       FACT(7) = ZERO
       FACT(8) = ZERO
       FACT(9) = ZERO
       FACT(10)= ZERO
       FACT(11)= ZERO
       FACT(12)= ZERO
      END IF


C===========================================================================
C  -- P. XU:  Dispersion using molecular polarizability 
C             E6 -- anisotropic molecular E6
C             E6_iso -- isotropic molecular E6
C             E7 -- anisotropic molecular E7
C             E8 -- anisotropic molecular E8
C             E8_iso -- isotropic molecular E8
      IF(MOLPOL) THEN

        CALL EFPDYN_MolPol(FACT)

      ENDIF  !MOLPOL
C==============================================================================

C ---  dispersion using LMO polarizability -DISTRIBUTED MODEL

      EDISP         = ZERO ! damped total dispersion E6+E7+1/3E6 
      EDISP0        = ZERO ! undamped total dispersion E6+E7+1/3E6 
      EDISP6        = ZERO ! damped isotropic E6  
      EDISP60       = ZERO ! undamped isotropic E6
      EDISP6_aniso  = ZERO ! damped anisotropic E6
      EDISP60_aniso = ZERO ! undamped anisotropic E6  
      EDISP7        = ZERO ! damped anisotropic E7
      EDISP70       = ZERO ! undamped anisotropic E7
      EDISP8_aniso  = ZERO ! damped anisotropic E8
      EDISP80_aniso = ZERO ! undamped anisotropic E8
      EDISP8        = ZERO ! damped isotropic E8
      EDISP80       = ZERO ! undamped isotropic E8

      E8DO          = ZERO ! dip-octo contribution to E8
      E8DQ          = ZERO ! dip-quad contribution to E8
      E8QQ          = ZERO ! quad-quad contribution to E8
c-----PX
!      t_Ttensor = zero
!      t_damp = zero
!      t_E6 = zero
!      t_E7 = zero
!      t_E8 = zero
!      t_aniso8 = zero
!      t_iso8 = zero
!      t_grad = zero
c-----END PX
C
C     -- TRACE OF DYNAMIC ALPHA IN THE FORM ALPHA(NDPPTS,NFRG) ---
C        used for calculating isotropic distributed E6
C
      ID=0
!$omp parallel do private(I,J) shared(DYNDD_LMO_ROT,DYNTR)
!$omp& firstprivate(ID) schedule(dynamic)
      DO I=1,NFRG
        DO J=1,NDPPTS(I)
           ID=ID+1
            DYNTR(J,I) = (ONE/THREE)*(DYNDD_LMO_ROT(1,1,J,I) +      
     *                                DYNDD_LMO_ROT(2,2,J,I) +
     *                                DYNDD_LMO_ROT(3,3,J,I))


        ENDDO
      ENDDO
!$omp end parallel do
C
c      -- Spherical LMO quad-quad polarizability --
c         for calculating isotropic distributed E8 (local)
      IF(E8DISP) THEN
      CALL VCLR(DYNQQ_LMO_AVE,1,MXIFRQ*MXPPTS*NFRG)
!$omp parallel do private(Ifrag,Ifreq,IL,I,J,K,L,INPPTS)
!$omp& shared(DYNQQ_LMO_ROT,DYNQQ_LMO_AVE,NPPTS)
!$omp& schedule(dynamic)

      Do Ifrag = 1, NFRG
         INPPTS = NPPTS(Ifrag)
       Do Ifreq = 1, MXIFRQ
        Do IL   = 1, INPPTS
           Do L = 1,3
           Do K = 1,3
           Do J = 1,3
           Do I = 1,3
           DYNQQ_LMO_AVE(IL+(Ifreq-1)*INPPTS,Ifrag)=
     *     DYNQQ_LMO_AVE(IL+(Ifreq-1)*INPPTS,Ifrag) + 
     *     (One/Five)*DYNQQ_LMO_ROT(I,J,K,L,IL+(Ifreq-1)*INPPTS,Ifrag)
     *     *((One/Two)*(delta(I,K)*delta(J,L)+delta(I,L)*delta(J,K))
     *     -(One/Three)*delta(I,J)*delta(K,L))
           Enddo
           Enddo
           Enddo
           Enddo
        Enddo
       Enddo
      Enddo
!$omp end parallel do
      ENDIF ! E8DISP

C --- P. XU: dispersion using LMO polarizability -DISTRIBUTED MODEL
      IP1=0
      IF (GOPARR) THEN
         DO I=1,NLOW-1
            IP1 = IP1 + NPPTS(I)
         ENDDO
      ENDIF

      DO 680 IM = NLOW,NUP
         INPPTS=NPPTS(IM)
        DO 670 IN=1,NPPTS(IM)
            IP1=IP1+1
            X1 = EFP(1,IP1)
            Y1 = EFP(2,IP1)
            Z1 = EFP(3,IP1)
            COMX1 = EFCENT(1,IM)
            COMY1 = EFCENT(2,IM)
            COMZ1 = EFCENT(3,IM)
            ALP1=DYNTR(IN,IM)
            ALP2=DYNTR(IN+INPPTS,IM)
            ALP3=DYNTR(IN+2*INPPTS,IM)
            ALP4=DYNTR(IN+3*INPPTS,IM)
            ALP5=DYNTR(IN+4*INPPTS,IM)
            ALP6=DYNTR(IN+5*INPPTS,IM)
            ALP7=DYNTR(IN+6*INPPTS,IM)
            ALP8=DYNTR(IN+7*INPPTS,IM)
            ALP9=DYNTR(IN+8*INPPTS,IM)
            ALP10=DYNTR(IN+9*INPPTS,IM)
            ALP11=DYNTR(IN+10*INPPTS,IM)
            ALP12=DYNTR(IN+11*INPPTS,IM)

            IF(IFRCPNT) THEN
               CALL NEARFP(IM,X1,Y1,Z1,IFP1,IFP1_COM,NEQUAL1)
            ENDIF

           IP2=0
          DO 660 JM = 1,NFRG
             IF (IM.GE.JM)THEN
                IP2=IP2+NPPTS(JM)
             !  GO TO 660
                CYCLE
             END IF
               JNPPTS=NPPTS(JM)
               COMX2 = EFCENT(1,JM)
               COMY2 = EFCENT(2,JM)
               COMZ2 = EFCENT(3,JM)
               COMX = COMX1 - COMX2
               COMY = COMY1 - COMY2
               COMZ = COMZ1 - COMZ2
               PBCX = ZERO
               PBCY = ZERO
               PBCZ = ZERO
               IF (XBOX .NE. ZERO) PBCX = XBOX * ANINT(COMX/XBOX)
               IF (YBOX .NE. ZERO) PBCY = YBOX * ANINT(COMY/YBOX)
               IF (ZBOX .NE. ZERO) PBCZ = ZBOX * ANINT(COMZ/ZBOX)
C
C   ----  DISTANCE BETWEEN FRAGMENTS
               X12  = ABS(COMX - PBCX)
               Y12  = ABS(COMY - PBCY)
               Z12  = ABS(COMZ - PBCZ)
               R12  = X12*X12+Y12*Y12+Z12*Z12
               R12SQ=SQRT(R12)

               SWF   = 1.0D+00
               SWFDX = 0.0D+00
               SWFDY = 0.0D+00
               SWFDZ = 0.0D+00
C
               IF(X12.GT.XCUT.OR.Y12.GT.YCUT.OR.Z12.GT.ZCUT)THEN
                 IP2=IP2+NPPTS(JM)
                 GOTO 660
               END IF
               R12  = X12*X12+Y12*Y12+Z12*Z12
               IF(R12.GT.RCUT)THEN
                 IP2=IP2+NPPTS(JM)
                 GOTO 660
               END IF
C
               IF(R12.GT.SWR2*SWR2) THEN
                 IP2=IP2+NPPTS(JM)
                 GOTO 660
               ELSE IF(R12.GT.SWR1*SWR1) THEN
                 IF(ISWF.EQ.1) THEN
                    DUMY1 = 1.0D+00/(SWR2*SWR2 - SWR1*SWR1)
                    DUMY2 = DUMY1*DUMY1
                    DUMY3 = DUMY2*DUMY1
                    DUMY4 = DUMY2*DUMY2
                    DUMY5 = DUMY2*DUMY3
                    VUMY1 = R12 - SWR1*SWR1
                    VUMY2 = VUMY1*VUMY1
                    VUMY3 = VUMY2*VUMY1
                    VUMY4 = VUMY2*VUMY2
                    VUMY5 = VUMY2*VUMY3
                    DSWF  = -30.0D+00*DUMY3*VUMY2
     *                      +60.0D+00*DUMY4*VUMY3
     *                      -30.0D+00*DUMY5*VUMY4
                    DSWF  = 2.0D+00*DSWF
                    SWF   = 1.0D+00 - 10.0D+00*DUMY3*VUMY3
     *                              + 15.0D+00*DUMY4*VUMY4
     *                              -  6.0D+00*DUMY5*VUMY5
                    SWFDX = DSWF*(EFCENT(1,IM)-EFCENT(1,JM)-PBCX)
                    SWFDY = DSWF*(EFCENT(2,IM)-EFCENT(2,JM)-PBCY)
                    SWFDZ = DSWF*(EFCENT(3,IM)-EFCENT(3,JM)-PBCZ)
                 ELSE IF(ISWF.EQ.2) THEN
                    R12   = SQRT(R12)
                    PIL   = 3.1415926535897932D+00/(SWR2-SWR1)
                    DUMY  = (R12 - SWR1) * PIL
                    SWF   = 0.50D+00 * (COS(DUMY) + 1.0D+00)
                    DUMY2 =-0.50D+00 * SIN(DUMY) * PIL / R12
                    SWFDX = DUMY2*(EFCENT(1,IM)-EFCENT(1,JM)-PBCX)
                    SWFDY = DUMY2*(EFCENT(2,IM)-EFCENT(2,JM)-PBCY)
                    SWFDZ = DUMY2*(EFCENT(3,IM)-EFCENT(3,JM)-PBCZ)
                 END IF
               END IF

             DO 650 JP = 1,NPPTS(JM)
                  IP2=IP2+1

                  ALPN1=DYNTR(JP,JM)
                  ALPN2=DYNTR(JP+JNPPTS,JM)
                  ALPN3=DYNTR(JP+2*JNPPTS,JM)
                  ALPN4=DYNTR(JP+3*JNPPTS,JM)
                  ALPN5=DYNTR(JP+4*JNPPTS,JM)
                  ALPN6=DYNTR(JP+5*JNPPTS,JM)
                  ALPN7=DYNTR(JP+6*JNPPTS,JM)
                  ALPN8=DYNTR(JP+7*JNPPTS,JM)
                  ALPN9=DYNTR(JP+8*JNPPTS,JM)
                  ALPN10=DYNTR(JP+9*JNPPTS,JM)
                  ALPN11=DYNTR(JP+10*JNPPTS,JM)
                  ALPN12=DYNTR(JP+11*JNPPTS,JM)

                  X2 = EFP(1,IP2)
                  Y2 = EFP(2,IP2)
                  Z2 = EFP(3,IP2)
C                   IF FORCE POINTS ARE ACTIVE
                  IF(IFRCPNT) THEN
                    CALL NEARFP(JM,X2,Y2,Z2,IFP2,IFP2_COM,NEQUAL2)
                  ENDIF

!PX   Prepare R's, where R is the distance between centroids 

                  X = X1 - X2 - PBCX
                  Y = Y1 - Y2 - PBCY
                  Z = Z1 - Z2 - PBCZ

                  COORD(1) = X
                  COORD(2) = Y
                  COORD(3) = Z

                  R2 = X*X + Y*Y + Z*Z
                  R = SQRT(R2)
                  R3 = R2*R
                  R4 = R2*R2
                  R6 = R4*R2
                  R5 = R3*R2
                  R7 = R5*R2
                  R8 = R6*R2
                  R9 = R7*R2
                  R10= R9*R

C   calculate T tensors

!--PX isotropic disp6 and disp8 evaluation 
!     does not need T tensors
                  IF(E7DISP) THEN
                     call T_tensor_2(COORD,T2)
                     call T_tensor_3(COORD,T3)
                  ENDIF
                  IF(do_aniso_disp) THEN
                    call T_tensor_4(COORD,T4)
                  ENDIF


C   Damping functions
                  CALL Damping_for_Dispersion(F6,F7,F8,DF6,DF7A,DF7B,
     *                                  DF8A,DF8B,IM,JM,IN,JP,
     *                                  IDISPS,E7DISP,E8DISP,IFGRAD,R,
     *                                  SAB,DSAB,SABSAVE,SABGRAD,MXPPTS)
!
                  EPAIR60_aniso = ZERO
                  EPAIR6_aniso  = ZERO
                  CALL Disp6_LMOpol(IN,JP,INPPTS,JNPPTS,IM,JM,FACT,
     *                          FACT1,FACT2,FACT3,FACT4,FACT5,FACT6,
     *                          FACT7,FACT8,FACT9,FACT10,FACT11,FACT12,
     *                          ALP1,ALP2,ALP3,ALP4,ALP5,ALP6,ALP7,
     *                          ALP8,ALP9,ALP10,ALP11,ALP12,ALPN1,
     *                          ALPN2,ALPN3,ALPN4,ALPN5,ALPN6,ALPN7,
     *                          ALPN8,ALPN9,ALPN10,ALPN11,ALPN12,
     *                          F6,R6,T2,C6,EPAIR60,EPAIR6,
     *                          EPAIR60_aniso,EPAIR6_aniso,
     *                          do_aniso_disp)
                  EDISP60 = EDISP60 + SWF*EPAIR60
                  EDISP6  = EDISP6  + SWF*EPAIR6

!     approximate Disp8 as 1/3 of Disp6
                  IF(.not.E8DISP) THEN
                  EPAIR=EPAIR6*FOUR/THREE
                  EPAIR0=EPAIR60*FOUR/THREE
                  ELSE ! not use 1/3 disp6 approximation
                  EPAIR=EPAIR6
                  EPAIR0=EPAIR60
                  ENDIF

                  EDISP60_aniso = EDISP60_aniso + EPAIR60_aniso
                  EDISP6_aniso  = EDISP6_aniso  + EPAIR6_aniso
!
C ---- anisotropic E7 from LMO polarizability ----

                  IF (E7DISP) THEN
                      CALL Disp7_LMOpol(IN,JP,INPPTS,JNPPTS,IM,JM,IP1,
     *                                  FACT,F7,T2,T3,COORD,R2,R4,R7,R9,
     *                                  EPAIR70,EPAIR7,IFGRAD,
     *                                  DER_EPAIR70,DER_EPAIR70_IJ)
                      EDISP70 = EDISP70 +SWF*EPAIR70
                      EDISP7  = EDISP7 +SWF*EPAIR7
                  ENDIF ! E7DISP

c---------distributed E8 from LMO polarizability
       IF(E8DISP) THEN
          CALL Disp8_LMOpol(IN,JP,INPPTS,JNPPTS,IM,JM,FACT,
     *                      F8,R8,T2,T3,T4,DYNQQ_LMO_AVE,
     *                      DYNTR,MXPPTS,EPAIR80,EPAIR8,
     *                      EPAIR80_aniso,EPAIR8_aniso,C8,
     *                      do_aniso_disp,E8DO,E8DQ,E8QQ)
!     *                      t0_aniso8,t1_aniso8,t0_iso8,t1_iso8)
       EDISP80_aniso = EDISP80_aniso + SWF*EPAIR80_aniso
       EDISP8_aniso  = EDISP8_aniso  + SWF*EPAIR8_aniso
       EDISP80   = EDISP80 + SWF*EPAIR80
       EDISP8    = EDISP8  + SWF*EPAIR8
       END IF ! E8disp

C Derivative of E6+E8 with overlap damping and TT formula
        IF (IFGRAD) THEN
           CALL GRAD_EFP_DISP(C6,F6,DF6,F7,DF7A,DF7B,C8,F8,
     *                        DF8A,DF8B,X,Y,Z,R,R6,R7,R8,R9,R10,
     *                        DSAB,conx,cony,conz,SWF,SWFDX,SWFDY,SWFDZ,
     *                        COMX,COMY,COMZ,PBCX,PBCY,PBCZ,IP1,IM,JM,
     *                        E7DISP,E8DISP,IDISPS,EPAIR,DER_EPAIR,
     *                        EPAIR70,EPAIR7,DER_EPAIR70,DER_EPAIR70_IJ,
     *                        EPAIR8)

           if( iefmorun .gt. 0 .and. iefmo_agrad .gt. 0) then
                    call efmo_store_disp( idisps, ntmo, sab, im, in,
     *                   jm, jp,df6, c6, r6, f6,swf, efmo_weight,
     *                   iefmodim,conx, cony,conz,fact1,fact2,fact3,
     *                   fact4,fact5,fact6,fact7,fact8,fact9,fact10,
     *                   fact11,fact12,alp1,alp2,alp3,alp4,alp5,
     *                   alp6,alp7,alp8,alp9,alp10,alp11,alp12,
     *                   alpn1,alpn2,alpn3,alpn4,alpn5,alpn6,alpn7,
     *                   alpn8,alpn9,alpn10,alpn11,alpn12,nfrg,mxppts)
           endif

        ENDIF !IFGRAD

C  ---- GATHER GRADS AND TORQUES IN FORCE POINTS
        IF(IFRCPNT) THEN
      CALL DISP_EFP_FP(C6,F6,DF6,X,Y,Z,R,R6,R7,R8,
     *                 DSAB,SWF,SWFDX,SWFDY,SWFDZ,
     *                 IFP1,IFP1_COM,IFP2,IFP2_COM,NEQUAL1,NEQUAL2,
     *                 COMX,COMY,COMZ,PBCX,PBCY,PBCZ,IP1,
     *                 IDISPS,EPAIR,DER_EPAIR)

        END IF !IFRCPNT

c----PX
!      t_Ttensor = t_Ttensor + (t1_Ttensor - t0_Ttensor)
!      t_damp = t_damp + (t1_damp - t0_damp)
!      t_E6 = t_E6 + (t1_disp6 - t0_disp6)
!      t_E7 = t_E7 + (t1_disp7 - t0_disp7)
!      t_E8 = t_E8 + (t1_disp8 - t0_disp8)
!      t_aniso8 = t_aniso8 + (t1_aniso8 - t0_aniso8)
!      t_iso8 = t_iso8 + (t1_iso8 - t0_iso8)
!      IF(IFGRAD) t_grad = t_grad + (t1_grad - t0_grad)
c----END PX

 650           CONTINUE
 660        CONTINUE
 670    CONTINUE
 680  CONTINUE

      IF(GOPARR) THEN
         CALL DDI_GSUMF(677,EDISP60,1)
         CALL DDI_GSUMF(677,EDISP6,1)
         CALL DDI_GSUMF(677,EDISP70,1)
         CALL DDI_GSUMF(677,EDISP7,1)
         CALL DDI_GSUMF(677,EDISP60_aniso,1)
         CALL DDI_GSUMF(677,EDISP6_aniso,1)
         IF(E8DISP) THEN
           CALL DDI_GSUMF(677,EDISP80_aniso,1)
           CALL DDI_GSUMF(677,EDISP8_aniso,1)
           CALL DDI_GSUMF(677,EDISP80,1)
           CALL DDI_GSUMF(677,EDISP8,1)
           call ddi_gsumf(677,E8DO,1)
           call ddi_gsumf(677,E8DQ,1)
           call ddi_gsumf(677,E8QQ,1)

         ENDIF
         if(IDISPS .EQ. 1 .and. ntmo.gt. 0 .and.iefmorun .gt. 0
     *        .and. iefmo_agrad .gt. 0 ) then
            CALL DDI_GSUMF(678,efmo_weight,nfrg*mxppts*nfrg*mxppts)
         endif
C
      END IF

c----PX debug printing (keep)
!      IF(MASWRK) THEN
!c    PX debug: different contributions in Disp8
!        write(6,*) 'E8DO=',E8DO
!        write(6,*) 'E8DQ=',E8DQ
!        write(6,*) 'E8QQ=',E8QQ
!c     -----------
!        write(6,*) 'distributed aniso E6=',EDISP6_aniso
!        write(6,*) 'distributed aniso E8=',EDISP8_aniso
!        write(6,*) 'distributed spher E8=',EDISP8
!      ENDIF
c-----END PX

c----PX
!      if(maswrk) then
!        write(6,101) 'time for mol. disp',t_mol
!        write(6,101) 'time for T tensors',t_Ttensor
!        write(6,101) 'time for damping',t_damp
!        write(6,101) 'time for disp6',t_E6
!        write(6,101) 'time for disp7',t_E7
!        write(6,101) 'time for disp8',t_E8
!        write(6,101) 'time for aniso8',t_aniso8
!        write(6,101) 'time for iso8',t_iso8
!        write(6,101) 'time for gradient',t_grad
!      endif
 101  FORMAT(A,F12.4)
!------PX take away 1/3 approximated disp8
!--PX: currently has to use approximated disp8 for qm-efp disp
  
      IF((.not.E8DISP).or.IQMDISP.eq.1) THEN
!      IF(.not.E8DISP) THEN
      EDISP8 = EDISP6/THREE
      EDISP80= EDISP60/THREE 
      EDISP = (EDISP6*(FOUR/THREE))+EDISP7
      EDISP0= (EDISP60*(FOUR/THREE))+EDISP70
      ELSE 
      !  write(6,*) 'EDISP8 in EFPDYN=',EDISP8
      EDISP = EDISP6+EDISP7+EDISP8
      EDISP0= EDISP60+EDISP70+EDISP80
      ENDIF

      RETURN
      END
C*MODULE EFDRVR  *DECK POLPCM
      SUBROUTINE POLPCM(EFLD,EFADD,QSE,QSED,QSN,
     *                  XCTS,YCTS,ZCTS,AS,
     *                  XE,YE,ZE,RE,ISPHE,LIST,POLNRG)
C
      USE MX_LIMITS,ONLY:MXFRG,MXFGPT,MXIFRQ,MXDFG,MXDPPT,mxatm
      USE comm_EFPPAR
      USE comm_EFMULT
      USE comm_FRGINF
      USE comm_PCMEFP
 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C
C
      LOGICAL GOPARR, DSKWRK, MASWRK
C
      COMMON /FMCOM / X(1)
      COMMON /FRFD  / PPEFX(MXFGPT), PPEFY(MXFGPT), PPEFZ(MXFGPT)
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON /OUTPUT/ NPRINT,ITOL,ICUT,NORMF,NORMP,NOPK
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /PCMDIM/ MXSP,MXTS,MEMPCM1,MEMPCM2,NTS
      COMMON /PCMITR/ RCUT(2),THRES,IPCMIT,IMUL,MXDIIS,NREG,MXITR1,
     *                MXITR2,MODPAP
      COMMON /PCMOPT/ RABI,RASC,REFPOL,THRSLS,DENSLS,WB,WA,ETA2,GD,EVAC,
     *                RHOW,PM,AREATL,AREAKP,BONDRY,OMEGA,RET,FRO,EPSINF,
     *                EPS,DR,RSOLV,VMOL,TCE,STEN,DSTEN,CMF,TABS,IDIRCT,
     *                IPCDER,IDP,ICOMP,IFIELD,ICAV,IDISP,IPRINT,IRETCAV,
     *                ICENT,IFAST,NEVAL,IEFPOL,KEEPSM,IMGABI,IMGASC,NADD
      COMMON /PCMPAR/ IPCM,NFT26,NFT27,IKREP,IEF,IP_F,NFMOPCM,IHET
      COMMON /PCMPRT/ GCAVP,GCAVS,GDISP,GREP,EHFGAS
      COMMON /SIMDAT/ NACC,NREJ,IGOMIN,NRPA,IBWM,NACCT,NREJT,NRPAT,
     *                NPRTGO,IDPUNC,IGOFLG
C
      DIMENSION  EFLD(3,NPTTPT),EFADD(3,NPTTPT),
     *           QSE(*),QSED(*),QSN(*),XCTS(*),YCTS(*),ZCTS(*),
     *           AS(*),XE(*),YE(*),
     *           ZE(*),RE(*),ISPHE(*),LIST(*)
C
      PARAMETER (ZERO = 0.00D+00)
      PARAMETER (ONE = 1.00D+00)
      PARAMETER (TOANGS=0.52917724924D+00)
      PARAMETER (ANTOAU=1.0D+00/TOANGS)
      PARAMETER (GLIM=1.0D-03)
C
C     -------------------------------------------
C       -- EFP/PCM (NO QM ATOM) CALCULATION --
C       EFP INDUCED DIPOLE AND PCM ASC ARE
C       ITERATED TO SELF-CONSISTENCY
C       HUI LI
C       IOWA STATE UNIVERSITY
C       JULY 25TH, 2005
C     -------------------------------------------
C
      NPRTBK = NPRINT
      NPRINT = 817
      CALL MAKCVM(0)
      CALL PEDRAM
      NPRINT = NPRTBK
C
      MXDII1 = MXDIIS + 1
C
      CALL VALFM(LOADFM)
      LQPOT1 = LOADFM + 1
      LQPOT2 = LQPOT1 + NTS
      LQ0    = LQPOT2 + NTS
      LQ1    = LQ0    + NTS
      LQ2    = LQ1    + NTS
      LQ3    = LQ2    + NTS
      LD0    = LQ3    + NTS
      LQA    = LD0    + NTS
      LDIMAT = LQA    + NTS
      LQREP  = LDIMAT + MXDII1*MXDII1
      LTMP0  = LQREP  + 2*NTS*MXDII1
      LTMP1  = LTMP0  + NTS*3
      LTMP2  = LTMP1  + MXDII1
      LIPVT  = LTMP2  + MXDII1*MXDII1
      LPOTTMP= LIPVT  + MXDII1
      LRMUL  = LPOTTMP+ NTS
      LAST   = LRMUL  + MXSP*10
      NEED   = LAST   - LOADFM -1
      CALL GETFM(NEED)
C
C
C     -- EFP STATIC MULTIPOLE POTENTIAL AT EACH TESSERA
C
      IF(IEFC+IEFD+IEFQ.GT.0) THEN
         DO ITS=1,NTS
            XI=XCTS(ITS)
            YI=YCTS(ITS)
            ZI=ZCTS(ITS)
            QADD=ZERO
            DO I=1,NMTTPT
               XX=EFQAD(1,I)
               YY=EFQAD(2,I)
               ZZ=EFQAD(3,I)
               XY=EFQAD(4,I)
               XZ=EFQAD(5,I)
               YZ=EFQAD(6,I)
               DUM=XX+YY+ZZ
               QXX=(3.0D+00*XX-DUM)/2.0D+00
               QYY=(3.0D+00*YY-DUM)/2.0D+00
               QZZ=(3.0D+00*ZZ-DUM)/2.0D+00
               QXY=1.50D+00*XY
               QXZ=1.50D+00*XZ
               QYZ=1.50D+00*YZ
               DX = XCTS(ITS)-EFC(1,I)
               DY = YCTS(ITS)-EFC(2,I)
               DZ = ZCTS(ITS)-EFC(3,I)
               R2 =  (EFC(1,I)-XI)**2
     *              +(EFC(2,I)-YI)**2
     *              +(EFC(3,I)-ZI)**2
               R  = SQRT(R2)
               QADD = QADD + EFCHG(1,I)/R + EFCHG(2,I)/R
     *                       +(EFDIP(1,I)*DX
     *                       +EFDIP(2,I)*DY
     *                       +EFDIP(3,I)*DZ)/(R*R2)
     *                      +(QXX*3.0D+00*DX*DX
     *                       +QYY*3.0D+00*DY*DY
     *                       +QZZ*3.0D+00*DZ*DZ
     *                       +QXY*6.0D+00*DX*DY
     *                       +QXZ*6.0D+00*DX*DZ
     *                       +QYZ*6.0D+00*DY*DZ)
     *                       /(3.0D+00*R2*R2*R)
            ENDDO
            X(LQPOT1+ITS-1)=-QADD
         ENDDO
      END IF
C
C
      CALL REGINI(ISPHE,LIST)
C
C
C     -- FIELD AT POLARIZABILITY POINTS DUE TO
C        MULTIPOLE POINTS IN OTHER FRAGMENTS
C        NOTE: THIS IS DONE ONCE AND REPEATEDLY
C              USED TO INDUCE EFP DIPOLES
C
      IF(IEFP.GT.0) CALL FRGFLD( 0, .FALSE. )
C
C
C     IF THERE IS NO DIPOLE POLARIZABILITY TENSOR,
C     NO ITERATION IS NEEDED.
C
      MXCALL = 1
      IF(IEFP.GT.0) MXCALL = 40
C
      DO 600 ICALL = 1, MXCALL
C        USUALLY 10~15 ITERATIONS ARE REQUIRED.
C
C        -- 1. ADD INDUCED DIPOLE POTENTIAL TO EACH TESSERA
         DO ITS=1,NTS
            XI=XCTS(ITS)
            YI=YCTS(ITS)
            ZI=ZCTS(ITS)
            QADD=ZERO
            DO LEFP=1,NPTTPT
               XL=EFP(1,LEFP)
               YL=EFP(2,LEFP)
               ZL=EFP(3,LEFP)
               X1 = XI - XL
               Y1 = YI - YL
               Z1 = ZI - ZL
               R2=X1*X1 + Y1*Y1 + Z1*Z1
               R = SQRT(R2)
               R3 = R*R2
               QADD= QADD + (DIND(1,LEFP)*X1
     *                    +  DIND(2,LEFP)*Y1
     *                    +  DIND(3,LEFP)*Z1)/R3/2.0D+00
     *                    + (DINDD(1,LEFP)*X1
     *                    +  DINDD(2,LEFP)*Y1
     *                    +  DINDD(3,LEFP)*Z1)/R3/2.0D+00
            ENDDO
            X(LQPOT2+ITS-1)=X(LQPOT1+ITS-1)-QADD
         ENDDO
C
C        -- 2. COMPUTE PCM SURFACE CHARGE ON EACH TESSERA
C
C        ITERATIVE CYCLE (STEP 1)
         NTSPAR=(NTS-1)/NPROC+1
         CALL ASCCYC(1,ICALL,X(LQPOT2),X(LQ0),X(LQ1),X(LD0),
     *              X(LQA),X(LDIMAT),X(LQREP),X(LTMP0),X(LTMP1),
     *              X(LTMP2),X(LIPVT),X(LRMUL),
     *              MXDII1,NTSPAR,XCTS,YCTS,ZCTS,AS,
     *              XE,YE,ZE,RE,ISPHE,LIST)
C
C        C-PCM ONLY NEEDS STEP 1.
         IF(IEF.EQ.10) THEN
           DO ITS=1,NTS
              QSE(ITS)    =X(LQ1+ITS-1)
              X(LQ0+ITS-1)=X(LQ1+ITS-1)
           ENDDO
           GOTO 500
         END IF
C
C        ITERATIVE CYCLE (STEP 2)
C
C        WRITE THE INTERMEDIATE POTENTIAL, THE STEP 2 INITIAL GUESS
C        AND COPY THE INTERMEDIATE CHARGES
C
         FACIP=-1.0D+00/(EPS-1.0D+00)
         FACIG=-1.0D+00/EPS
         DO ITS=1,NTS
            X(LQ0+ITS-1)=      X(LQ2+ITS-1)
            X(LQ2+ITS-1)=FACIP*X(LQ1+ITS-1)
            X(LQ3+ITS-1)=      X(LQ1+ITS-1)
            IF(ICALL.EQ.1) X(LQ0+ITS-1)=FACIG*X(LQ1+ITS-1)
         ENDDO
C
         CALL ASCCYC(2,ICALL,X(LQ2),X(LQ0),X(LQ1),X(LD0),
     *              X(LQA),X(LDIMAT),X(LQREP),X(LTMP0),X(LTMP1),
     *              X(LTMP2),X(LIPVT),X(LRMUL),
     *              MXDII1,NTSPAR,XCTS,YCTS,ZCTS,AS,
     *              XE,YE,ZE,RE,ISPHE,LIST)
C
C        COPY THE CHARGES
C
         DO ITS=1,NTS
            QSE(ITS)=X(LQ3+ITS-1)+X(LQ1+ITS-1)
         ENDDO
C
C        USE Q1 AND Q3 FOR NEXT GUESSES
C
         DO ITS=1,NTS
            X(LQ0+ITS-1)=X(LQ3+ITS-1)
            X(LQ2+ITS-1)=X(LQ1+ITS-1)
         ENDDO
C
 500     CONTINUE
C
C        -- SCALE C-PCM CHARGE BY (EPS-1)/EPS
         IF(IEF.EQ.10) THEN
           DO ITS = 1, NTS
              QSE(ITS) = QSE(ITS)*(EPS-ONE)/EPS
           END DO
         END IF
C
C        TEST PCM ASC FOR CONVERGENCY
C        NOTE: WE USE QSN TO SAVE THE OLD QSE
C              MAKE SURE TO CLEAR QSN AFTERMATH
C
         SDEV=ZERO
         DO ITS=1,NTS
            SDEV=SDEV+(QSE(ITS)-QSN(ITS))**2
            QSN(ITS) = QSE(ITS)
         ENDDO
         SDEV=SQRT(SDEV/NTS)
C           IF(MASWRK) WRITE(IW,*) ' THRES=',THRES,
C    *      ' SDEV=', SDEV
         IF (SDEV.EQ.ZERO) THEN
C           IF(MASWRK) WRITE(IW,*)
C    *      'POLPCM CONVERGED AFTER ',ICALL,'ITERATIONS'
            GO TO 700
         END IF
C
C
C        **********************
C        UPDATE INDUCED DIPOLES
C        **********************
C
         IF(IEFP.GT.0) THEN
C
C           -- TOTAL FIELDS AT POL-TENSORS
C
            IPCOUNT = ME-1
            DO 200 LEFP = 1, NPTTPT
C              -- GO PARALLEL!
               IF(GOPARR) THEN
                  EFLD(1,LEFP) = ZERO
                  EFLD(2,LEFP) = ZERO
                  EFLD(3,LEFP) = ZERO
                  IPCOUNT = IPCOUNT + 1
                  IF (MOD(IPCOUNT,NPROC).NE.0) GO TO 200
               END IF
C              -- SOLVENT APPARENT CHARGE FIELD AT THE POINT --
               PQ_X=ZERO
               PQ_Y=ZERO
               PQ_Z=ZERO
               IF(IEFPOL.EQ.1) THEN
                  DO 330 ITS = 1,NTS
                     RR = (EFP(1,LEFP)-XCTS(ITS))**2 +
     *                    (EFP(2,LEFP)-YCTS(ITS))**2 +
     *                    (EFP(3,LEFP)-ZCTS(ITS))**2
                     IF(RR.LT.GLIM .OR. RR.LT.(REFPOL*ANTOAU))GOTO 330
                     RR3 = RR * SQRT(RR)
                     PQ_X = PQ_X + QSE(ITS)*(EFP(1,LEFP)-XCTS(ITS))/RR3
                     PQ_Y = PQ_Y + QSE(ITS)*(EFP(2,LEFP)-YCTS(ITS))/RR3
                     PQ_Z = PQ_Z + QSE(ITS)*(EFP(3,LEFP)-ZCTS(ITS))/RR3
 330              CONTINUE
               END IF
C              -- TOTAL FIELD: EFP+PCM --
               EFLD(1,LEFP) = PQ_X+PPEFX(LEFP)
               EFLD(2,LEFP) = PQ_Y+PPEFY(LEFP)
               EFLD(3,LEFP) = PQ_Z+PPEFZ(LEFP)
 200        CONTINUE
            IF(GOPARR) THEN
               CALL DDI_GSUMF(667,EFLD, 3*NPTTPT)
            END IF
C
C           -- COMPUTE DIND/DINDD
C
            CALL VALFM(LOADFM)
            LDIPNEW  = LOADFM   + 1
            LDIPNWD  = LDIPNEW  + NFRG
            LDINDC   = LDIPNWD  + NFRG
            LDINDDC  = LDINDC   + 3*NPTTPT
            LDIND_W  = LDINDDC  + 3*NPTTPT
            LDIND_C  = LDIND_W  + 6*NPTTPT
            LEWLD1   = LDIND_C  + 6*NPTTPT
            LEWLD2   = LEWLD1   + 3*NPTTPT
            LAST     = LEWLD2   + 3*NPTTPT
            NEEED    = LAST - LOADFM - 1
            CALL GETFM(NEEED)
            CALL DIPIT(EFLD,EFADD,X(LDINDC),X(LDINDDC),
     *                 X(LDIND_W),X(LDIND_C),
     *                 NPTTPT,X(LDIPNEW),X(LDIPNWD),
     *                 X(LEWLD1),X(LEWLD2),NFRG)
            CALL RETFM(NEEED)
            IF(IDPUNC.EQ.1) RETURN
         END IF
C
 600  CONTINUE
 700  CONTINUE
      CALL VCLR(QSN,1,NTS)
C
C
C     -- COMPUTE PCM INTERACTION
C
      ENPCM  = ZERO
      QET    = ZERO
      DO ITS = 1,NTS
         QET    = QET+QSE(ITS)
         ENPCM  = ENPCM - X(LQPOT1+ITS-1) * QSE(ITS)
      ENDDO
C
C     -- COMPUTE EFP POLARIZATION ENERGY
C
      POLNRG = ZERO
      DO II=1,NPTTPT
         POLNRG = POLNRG - DIND(1,II)*PPEFX(II)
     *                   - DIND(2,II)*PPEFY(II)
     *                   - DIND(3,II)*PPEFZ(II)
     *                   - DINDD(1,II)*PPEFX(II)
     *                   - DINDD(2,II)*PPEFY(II)
     *                   - DINDD(3,II)*PPEFZ(II)
      END DO
      POLNRG = POLNRG/4.0D+00
C
C     SOLVATION OUTPUT FOR PCM MODEL
C
      IF(MASWRK) THEN
        WRITE(IW,*) ' '
        WRITE(IW,*) ' NTS    =', NTS
        WRITE(IW,'(A10,F18.10)') 'PCM ASC=',QET
        WRITE(IW,'(A10,F18.10)') 'GELEC  =',ENPCM/2.0D+00
        IF(ICAV.EQ.1)
     *  WRITE(IW,'(A10,F18.10)') 'GCAVP  =',GCAVP/627.509541D+00
        IF(IDISP.EQ.1)THEN
        WRITE(IW,'(A10,F18.10)') 'GDISP  =',GDISP/627.509541D+00
        WRITE(IW,'(A10,F18.10)') 'GREP   =', GREP/627.509541D+00
        END IF
      END IF
C
C     - QSED IS REQUIRED EVERYWHERE EVEN IF IHET=0 -
C     - QSED = QSE WHEN IHET=0
      IF(IHET.EQ.0) CALL DCOPY(NTS,QSE,1,QSED,1)
      IF(IHET.NE.0) THEN
         IF(MASWRK) WRITE(IW,*)
     *         ' ONLY HOMOGENEOUS PCM CAN BE USED FOR EFP-PCM'
         CALL ABRT
      END IF
C
      CALL RETFM(NEED)
      RETURN
      END
C
C
C*MODULE EFDRVR  *DECK PREFP
      SUBROUTINE PREFP(LU)
C
      USE MX_LIMITS,ONLY:MXFRG,MXFGPT,MXDFG,MXDPPT
      USE comm_FGRAD
      USE comm_FRGINF
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      LOGICAL GOPARR,DSKWRK,MASWRK
C
      CHARACTER*6 FRGNAM
      CHARACTER*8 PTNAM
C
      COMMON /FRAGMT/ XCRD(50,MXFRG),YCRD(50,MXFRG),ZCRD(50,MXFRG),
     *                PTNAM(50,MXFRG),FRGNAM(MXFRG)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
C
C   ---- LVS 2007: ROUTINE FOR PRINTING IN MD AND MC
C
      DO 30 I=1,NFRG
       IF (MASWRK) THEN
          WRITE(LU,35) I, FRGNAM(I)
          WRITE(LU,34) EFCENT(1,I),EFCENT(2,I),EFCENT(3,I)
          WRITE(LU,33) DEFT(1,I),DEFT(2,I),DEFT(3,I)
          WRITE(LU,31) TORQ(1,I),TORQ(2,I),TORQ(3,I)
          WRITE(LU,36)
       END IF
 30   CONTINUE
C
 35   FORMAT('FRAGMENT #',I6,'  ',A6)
 31   FORMAT('      TORQ',3F16.10)
 33   FORMAT('     FORCE',3F16.10)
 34   FORMAT('    EFCENT',3F16.10)
 36   FORMAT('    ')
C
      RETURN
      END
C
C*MODULE EFDRVR  *DECK EFPPARL
      SUBROUTINE EFPPARL
C
      USE MX_LIMITS,ONLY:MXFRG,MXDFG,MXDPPT
      USE comm_EFPAR
      USE comm_EFPAR2
      USE comm_FRGINF
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      LOGICAL GOPARR,DSKWRK,MASWRK
C
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      INTEGER DDI_NP,DDI_ME,FLO,CEIL,DIFF
C
C     PARTITION EFP'S OVER PROCESSORS
C
C     THIS SETS UP THE PARAMETERS FOR EFP PARALLELIZATION SCHEME.
C     WHEN ONE USES GDDI COMMUNICATER, THE PARALLELIZATION PARAMETER
C     FOR EFP MUST BE RESET, IN THE CORRECT COMMUNICATOR.
C
C      write(IW,*)'NFRG=',NFRG
      IF(NFRG.GT.0) THEN
        DDI_ME = ME
        DDI_NP = NPROC
        IF (MOD(NFRG,DDI_NP) .EQ. 0) THEN
          LOCNFG = NFRG / DDI_NP
          NLOW = DDI_ME * LOCNFG + 1
        ELSE
          CEIL = NFRG/DDI_NP + 1
          FLO = NFRG/DDI_NP
          DIFF = NFRG - DDI_NP*FLO
C
          IF (DDI_ME .GT. (DDI_NP-1-DIFF)) THEN
            LOCNFG = CEIL
            NLOW = DDI_ME * LOCNFG - (DDI_NP-1-DIFF)
          ELSE
            LOCNFG = FLO
            NLOW = DDI_ME * LOCNFG + 1
          END IF
        END IF
C
        NUP =  NLOW + LOCNFG - 1
        IF(GOPARR) THEN
C    THERE IS AN UNFORTUNATE COMMON BLOCK DIMENSION USED HERE.
          IF(DDI_NP.GT.MXCPUEFP) THEN
             IF(MASWRK) WRITE(IW,*)
     *          'ERROR IN PARALLEL EFP, MXCPUEFP,NPROC=',MXCPUEFP,DDI_NP
             CALL ABRT
          END IF
          CALL VICLR(NFRGAL,1,MXCPUEFP)
          CALL VICLR(NLOWAL,1,MXCPUEFP)
          CALL VICLR( NUPAL,1,MXCPUEFP)
          NFRGAL(DDI_ME+1) = LOCNFG
          NLOWAL(DDI_ME+1) = NLOW
           NUPAL(DDI_ME+1) = NUP
          CALL DDI_GSUMI(660,NFRGAL,DDI_NP)
          CALL DDI_GSUMI(661,NLOWAL,DDI_NP)
          CALL DDI_GSUMI(662, NUPAL,DDI_NP)
        ELSE
          NFRGAL(1) = NFRG
          NLOWAL(1) = 1
           NUPAL(1) = NFRG
        END IF
      END IF

C      write(IW,*)'NFRGAL=',NFRGAL
C      write(IW,*)'NLOWAL=',NLOWAL
C      write(IW,*)'NFRG=',NFRG
      RETURN
      END
C *MODULE EFDRVR  *DECK EFAIDISPDR
!> @brief   preparing memory for QMEFP2 (EFP-AI) dispersion energy
!>
!> @author  Tosaporn Sattasathuchana, Peng Xu
!> - Oct, 2019
!>
!> @param EORBS: orbital energy
!> @param L1: basis set size
!>
      SUBROUTINE EFAIDISPDR(EORBS,L1)
      USE MX_LIMITS,ONLY:MXDPPT,mxatm, MXIFRQ
      USE EFAIDIS
      USE comm_FRGINF
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /FMCOM / X(1)
      DIMENSION EORBS(NUM)
      L2 = (L1*L1+L1)/2
      NLOC = NA-MMCORE
      NVIR  = NQMT-NA
C
C    SETUP FREE MEMORY
C
      CALL VALFM(LOADFM)
      LWRKTR   = 1        + LOADFM
      LDYNTR   = LWRKTR   + MXIFRQ
      LFUM1    = LDYNTR   + MXDPPT*NFRG
      LFUM2    = LFUM1    + NLOC*NVIR
      LFUM3    = LFUM2    + NLOC*NVIR
      LSQ1     = LFUM3    + NLOC*NVIR
      LSQ2     = LSQ1     + NLOC*NLOC
      LSQ3     = LSQ2     + NLOC*NLOC
      LAVG     = LSQ3     + NLOC*NLOC
      LWORK    = LAVG     + NLOC*NLOC
      LH       = LWORK    + NA
      LUM1     = LH       + NLOC*NLOC
      LUM2     = LUM1     + NVIR*NLOC
      LUM3     = LUM2     + NVIR*NLOC
      LUM      = LUM3     + NVIR*NLOC
      LHF      = LUM      + NLOC*NVIR*3
      LTRAN    = LHF      + L2
      LE       = LTRAN    + NLOC*NLOC
      LHFP     = LE       + L1
      LCXM     = LHFP     + L1*L1
      LTM      = LCXM     + L1*L1
      LD       = LTM      + L1
      LVEC     = LD       + L2
      LCCHG    = LVEC     + L1*NA
      LCEE     = LCCHG    + 3*NLOC
      LZNU     = LCEE     + NAT*3
      LAROT    = LZNU     + NAT
      LATOTAL  = LAROT    + 3*NAT
      LFTOTAL  = LATOTAL  + 3*NAT
      LATRAN   = LFTOTAL  + 3*NAT
      LDEE     = LATRAN   + 3
      LDEETORQ = LDEE     + 3
      LXYZP    = LDEETORQ + 6
      LDIST    = LXYZP    +   3
      LAST     = LDIST    +   3
      NEEDV    = LAST - LOADFM - 1
      CALL GETFM(NEEDV)
      CALL EFAIDISP(EORBS,L1,L2,NVIR,NLOC,X(LWRKTR), X(LDYNTR),X(LFUM1),
     $      X(LFUM2),X(LFUM3), X(LSQ1), X(LSQ2), X(LSQ3), X(LAVG),
     $      X(LWORK),X(LH), X(LUM1), X(LUM2), X(LUM3), X(LUM), X(LHF),
     $      X(LTRAN), X(LHFP),X(LCXM), X(LTM), X(LD), X(LVEC), X(LCCHG),
     $      X(LCEE), X(LZNU), X(LAROT), X(LATOTAL), X(LFTOTAL),
     $      X(LATRAN), X(LDEE), X(LDEETORQ), X(LXYZP), X(LDIST))
      CALL RETFM(NEEDV)
      RETURN
      END SUBROUTINE EFAIDISPDR
C
C *MODULE EFDRVR  *DECK EFAIDISP
!> @brief   Calculates EFP-ab initio (EFP-AI/EFP-QM) dispersion energy
!>          and gradients
!>
!> @author  Tony Smith
!> - Jan 10, 2013
!> @date Oct 2019 - Peng Xu and Tosaporn Sattasathuchana
!> - adding argument in the calling lists to allow for the memory
!>   allocation. This subroutine is called by EFAIDISPDR.
!> 
!> @param   EORBS: Array of orbital energies
!> @details see Q. Smith et al, J. Chem. Phys. 136, 244107 (2012);
!>          http://dx.doi.org/10.1063/1.4729535
!> @note    Gradient code for overlap-based damping (IDISDMP=1) is
!>          incomplete!
!> @todo    Get derivatives of overlap integrals between EFP and AI
!>          orbitals to complete the overlap-based damping gradients
!> @date October 2018, Peng Xu and Tosaporn Sattasathuchana
!>  - Implement new overlap damping function and make it as default
!>       (IDISDMP=1)
!>  - bug fixing for loop index IEFP confliting with the common block
!>    variable (changed to IIEFP in do loop)
!>
!> @param EORBS: orbital energy
!> @param L1: basis set size
!> @param L2: triangular length of L1
!> @param nvir: number of virtual orbitals
!> @param nloc: number of active occupied orbitals
!> @param WRKTR: dynamic alpha polarizability tensors
!> @param DYNTR: average dynamic polarizability
!> @param FUM1: product of dipole integral (X) with the imaginary 
!>              frequency integral FR
!> @param FUM2: product of dipole integral (Y) with the imaginary 
!>              frequency integral FR
!> @param FUM3: product of dipole integral (Z) with the imaginary 
!>              frequency integral FR
!> @param SQ1: product of FUM1*UM1
!> @param SQ2: product of FUM2*UM2
!> @param SQ3: product of FUM3*UM3
!> @param AVG: average of SQ1, SQ2, and SQ3 (see comment below)
!> @param WORK: scratch space
!> @param H: transformation of AVG to localized MO basis
!> @param UM1: TRANSPOSE OF DIPOLE INTEGRAL MATRICES (X)
!> @param UM2: TRANSPOSE OF DIPOLE INTEGRAL MATRICES (Y)
!> @param UM3: TRANSPOSE OF DIPOLE INTEGRAL MATRICES (Z)
!> @param UM: DIPOLE INTEGRAL MATRICES
!> @param HF: a transformed dipole integral using CXM
!> @param TRAN: localization transformation matrix
!> @param HFP: square of the HF
!> @param CXM: alpha orbitals, read from DICTNRY 15
!> @param TM: working array
!> @param D: dipole integral
!> @param VEC: alpha orbitals during localuzation, from DICTNRY 71
!> @param CCHG: a transformed dipole integral using VEC
!> @param CEE: transpose of C
!> @param ZNU: nuclear charge
!> @param AROT: rotational force
!> @param ATOTAL: sum of ATRAN and AROT
!> @param FTOTAL: (mass-weighted) total force
!> @param ATRAN: translational force
!> @param DEE: translational gradient
!> @param DEETORQ: rotational gradient
!> @param XYZP: storage to hold center of mass data
!> @param DIST: storage for distances
!>
      SUBROUTINE EFAIDISP(EORBS,L1,L2,NVIR,NLOC,WRKTR, DYNTR, FUM1,
     $      FUM2,FUM3, SQ1, SQ2, SQ3,  AVG, WORK, H, UM1, UM2, UM3,
     $      UM, HF, TRAN, HFP, CXM, TM, D, VEC, CCHG, CEE, ZNU,
     $      AROT, ATOTAL,FTOTAL, ATRAN, DEE, DEETORQ, XYZP, DIST)

      USE MX_LIMITS,ONLY:MXDPPT,mxatm,MXIFRQ
      USE comm_EFPSCR
      USE comm_EFPPAR
      USE comm_EFPGES
      USE comm_FGRAD
      USE comm_FRGINF
      USE comm_EFPGES

      USE EFAIDIS
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON /GRAD  / DE(3,MXATM)
      COMMON /MASSES/ ZMASS(MXATM)
C
      DIMENSION EORBS(L1)
      DIMENSION WRKTR(MXIFRQ) 
      DIMENSION DYNTR(MXDPPT,NFRG)
      DIMENSION FUM1(NLOC,NVIR),FUM2(NLOC,NVIR),
     *          FUM3(NLOC,NVIR)
      DIMENSION SQ1(NLOC,NLOC),
     *          SQ2(NLOC,NLOC),
     *          SQ3(NLOC,NLOC),
     *          AVG(NLOC,NLOC),WORK(NA),
     *          H(NLOC,NLOC)
      DIMENSION UM1(NVIR,NLOC),UM2(NVIR,NLOC),
     *          UM3(NVIR,NLOC)
      DIMENSION UM(NLOC,NVIR,3),
     *          HF(L2),
     *          TRAN(NLOC,NLOC),
     *          HFP(L1,L1),CXM(L1,L1),TM(L1),
     *          D(L2)
      DIMENSION VEC(L1,NA),XYZP(3),CCHG(3,NLOC),DIST(3)
      DIMENSION CEE(NAT,3),ZNU(NAT)
      DIMENSION DEE(3),DEETORQ(6),AROT(3,NAT),
     *          ATRAN(3),ATOTAL(3,NAT),FTOTAL(3,NAT)
C
      DATA ZERO,ONE,TWO,THREE /0.0D+00,1.0D+00,2.0D+00,3.0D+00/
      DATA FOUR,FIVE,SIX /4.0D+00,5.0D+00,6.0D+00/
C
      THIRD=ONE/THREE
      FOURTHIRD=FOUR*THIRD
      MAXDER = 0
      CALL DERCHK(MAXDER)
      IF (MAXDER.NE.1) IDSGRD = 0
      NXYZ = 3
      NORB = NQMT
      NOCC = NA
      NLOC = NA - MMCORE
      NVIR = NORB - NA
      L3 = L1*L1
      DO I=1,3
         DEE(I)=ZERO
         DEETORQ(I)=ZERO
         DEETORQ(I+3)=ZERO
      ENDDO
      CALL DAREAD(IDAF,IODA,VEC,L1*NA,71,0)
      IF(NVIR .EQ. 0) RETURN
      CALL DAREAD(IDAF,IODA,TRAN,NLOC*NLOC,73,0)
      ID = 0
      DO I=1,NFRG
        DO J=1,NDPPTS(I)
           ID=ID+1
           DYNTR(J,I) = THIRD*(EFDPOL(1,ID) + EFDPOL(2,ID) +
     *                   EFDPOL(3,ID))
        ENDDO
      ENDDO
C
C  TRANSFORM DIPOLE INTEGRALS TO MO BASIS
C
      CALL DAREAD(IDAF,IODA,CXM,L3,15,0)
      IF (MMCORE.GT.0) THEN
         DO I = 1,NLOC
            CALL DCOPY(L1,VEC(1,MMCORE+I),1,VEC(1,I),1)
         END DO
      END IF
      CALL CALCOM(XCOM,YCOM,ZCOM)
      XYZP(1) = XCOM
      XYZP(2) = YCOM
      XYZP(3) = ZCOM
      CALL DIPINT(XCOM,YCOM,ZCOM,.FALSE.)
      DO 300 IXYZ=1,3
         CALL DAREAD(IDAF,IODA,D,L2,94+IXYZ,0)
         CALL TFTRI(HF,D,CXM,TM,NORB,L1,L1)
C GET (OCC-CORE) x VIR BLOCK
         CALL CPYTSQ(HF,HFP,L1,1)
         IL=0
         DO 320 IVIR=1,NVIR
            IM=0
            DO 340 IOCC=MMCORE+1,NA
               IM=IM+1
               UM(IM,IVIR,IXYZ)=HFP(IOCC,NA+IVIR)
  340       CONTINUE
  320    CONTINUE
C GET LMO CENTROIDS
         DO 360 J=1,NLOC
            CALL TFTRI(CCHG(IXYZ,J),D,VEC(1,J),TM,1,L1,L1)
            CCHG(IXYZ,J) = CCHG(IXYZ,J) + XYZP(IXYZ)
  360    CONTINUE
  300 CONTINUE
      IE=0
      SUM=0.0D+00
      DISP=0.0D+00
      JAY=0
C BEGIN LOOP OVER EFP POINTS
      DO 801 IIEFP=1,NFRG
      INPPTS=NPPTS(IIEFP)
      DO 800 JJEFP=1,INPPTS
      IE=IE+1
      DO IT=1,MXIFRQ
C CONSTRUCT ALPHA TENSOR (DYNAMIC POLARIZABILITY VALUES FOR EFP)
         WRKTR(IT)=DYNTR(JJEFP+(IT-1)*INPPTS,IIEFP)
      ENDDO
      DO 411 IOCC=1,NLOC
         DO 421 IVIR=1,NVIR
C ENERGY DIFF BETWEEN A VIRTUAL AND AN OCCUPIED AI ORBITAL
            ENRGDF=EORBS(NA+IVIR)-EORBS(MMCORE+IOCC)
C GET VALUE OF INTEGRAL OVER IMAG FREQ RANGE FOR THIS PAIR OF AI ORBITALS
C AND SET OF ALPHA VALUES
            CALL EFAIMF(FR,WRKTR,ENRGDF)
C MULTIPLY IMAG FREQ INTEGRAL WITH CORRESPONDING DIPOLE INTEGRAL
            FUM1(IOCC,IVIR)=UM(IOCC,IVIR,1)*FR
            FUM2(IOCC,IVIR)=UM(IOCC,IVIR,2)*FR
            FUM3(IOCC,IVIR)=UM(IOCC,IVIR,3)*FR
C TRANSPOSE OF DIPOLE INTEGRAL MATRICES
            UM1(IVIR,IOCC)=UM(IOCC,IVIR,1)
            UM2(IVIR,IOCC)=UM(IOCC,IVIR,2)
            UM3(IVIR,IOCC)=UM(IOCC,IVIR,3)
  421    CONTINUE
  411 CONTINUE
C MULTIPLY (IMAG FREQ*DIPOLE) BY (TRANSPOSE OF DIPOLE), OBTAINING SQUARE
C MATRICES INDEXED BY OCCUPIED (VALENCE) ORBITALS ONLY
      CALL MRARBR(FUM1,NLOC,NLOC,NVIR,UM1,NVIR,NLOC,SQ1,NLOC)
      CALL MRARBR(FUM2,NLOC,NLOC,NVIR,UM2,NVIR,NLOC,SQ2,NLOC)
      CALL MRARBR(FUM3,NLOC,NLOC,NVIR,UM3,NVIR,NLOC,SQ3,NLOC)
C  AVERAGE THE MATRICES OVER DIRECTIONAL COMPONENTS: (1/3)(Dx+Dy+Dz)
C  (The FOUR is somehow part of the integration scheme; see SUROUTINE LAPOL
C  in source file LOCPOL, from which this subroutine was adapted.)
      DO I=1,NLOC
         DO J=1,NLOC
            AVG(I,J)=SQ1(I,J)+SQ2(I,J)+SQ3(I,J)
            AVG(I,J)=AVG(I,J)*FOUR*THIRD
         ENDDO
      ENDDO
C TRANSFORM TO LOCALIZED MO BASIS
      CALL TFSQU(H,AVG,TRAN,WORK,NLOC,NLOC)
C DIAGONAL ELEMENTS OF H MATRIX NOW CORRESPOND TO C6 VALUES.
C USE THEM WITH C6*F6/R6 FORMULA:
      DO 790 I=1,NLOC
         C6=H(I,I)
         F6=ONE
         DO IXYZ=1,3
            DIST(IXYZ)=EFP(IXYZ,IE)-CCHG(IXYZ,I)
         ENDDO
         X = DIST(1)
         Y = DIST(2)
         Z = DIST(3)
C         RR=DIST(1)*DIST(1)+DIST(2)*DIST(2)+DIST(3)*DIST(3)
         R2 = X*X + Y*Y + Z*Z
         R=SQRT(R2)
         R6=R2*R2*R2
         IF (IDISDMP.EQ.2) THEN
C           OVERLAP-BASED DAMPING
            SAB=SABEFPAI(IIEFP,I,JJEFP)
            IF (ABS(SAB) .GT. 1.0D-05) THEN
               RB = -2.0D+00*LOG(ABS(SAB))
               RB2 = RB*RB
               EX1 = SAB*SAB
               F6 = ONE-(EX1*(ONE + RB + RB2/TWO ))
            ELSE
               F6 = ONE
               RB = ZERO
               RB2 = ZERO
               SAB = ZERO
            END IF
C           THE OVERLAP BASED DAMPING FOR DISPERSION
C REF. J. Phys. Chem. A, 120, 639-647 (2016)
c-----PX: new overlap damping option is 1 (default)
         ELSE IF (IDISDMP.EQ.1) THEN
            SAB=SABEFPAI(IIEFP,I,JJEFP)
            IF (ABS(SAB) .GT. 1.0D-05) THEN
               RB = -2.0D+00*LOG(ABS(SAB))
               RB2 = RB*RB
               RB3 = RB2*RB
               EX1 = SAB*SAB
               F6 = ONE-(EX1*(ONE + SQRT(RB) + (RB/TWO) +
     *       (((RB)**(THREE/TWO))/SIX) + (RB2/24.0D+00) +
     *       ((RB**(FIVE/TWO))/120.0D+00) +(RB3/720.0D+00)))
            ELSE
               F6 = ONE
               RB = ZERO
               RB2 = ZERO
               SAB = ZERO
            END IF
         ELSE IF (IDISDMP.EQ.-1) THEN
C           NO DAMPING
            F6=ONE
            RB=ZERO
            RB2 = ZERO
            SAB=ZERO
         ELSE
C           TANG-TOENNIES DAMPING
            RB =  B*R
            EX1 = EXP(-RB)
            RB2 = RB*RB
            RB3 = RB*RB2
            RB4 = RB2*RB2
            RB5 = RB*RB4
            RB6 = RB*RB5
C
            F6 = ONE-(EX1*(ONE + RB + RB2/TWO +
     *               RB3/SIX + RB4/24.0D+00 +
     *               RB5/120.0D+00 +RB6/720.0D+00))
         END IF
         SUM=SUM+C6*F6
         DISP=DISP - C6*F6/R6
C   GRADIENT
C   DE/DR=FOUR/THREE*(-6*C6*F6/R**7 + DF6/DR * C6/R**6)
         IF (IDSGRD.EQ.0) GOTO 790
         IF (IDISDMP .EQ. 1) THEN
C   NOTE: GRADIENT NOT FULLY IMPLEMENTED FOR OVERLAP-BASED
C         DAMPING. NEED "DF6" TERMS, I.E. DSAB (DERIVATIVE
C         OF OVERLAP INTEGRALS TAKEN BTWN EFP&QM ORBITALS)
            CON1=-FOURTHIRD*(C6/R6)*(SIX*F6/R2)
            CON2=-FOURTHIRD*(C6/R6)*(RB2*SAB)
         ELSE
            CON1=-FOURTHIRD*(C6/R6)*
     *           ((SIX*F6/R2)-(B*RB6*EX1)/(R*720.0D+00))
            CON2 = ZERO
         END IF
C
         CONX=CON1*X
         CONY=CON1*Y
         CONZ=CON1*Z
         DEE(1) = DEE(1) + CONX
C     *            + CON2*DSAB(1)
         DEE(2) = DEE(2) + CONY
C     *            + CON2*DSAB(2)
         DEE(3) = DEE(3) + CONZ
C     *            + CON2*DSAB(3)
         DEETORQ(1) = DEETORQ(1)+CONZ*(EFP(2,JJEFP)-EFCENT(2,IIEFP))
     *                      -CONY*(EFP(3,JJEFP)-EFCENT(3,IIEFP))
         DEETORQ(2) = DEETORQ(2)+CONX*(EFP(3,JJEFP)-EFCENT(3,IIEFP))
     *                      -CONZ*(EFP(1,JJEFP)-EFCENT(1,IIEFP))
         DEETORQ(3) = DEETORQ(3)+CONY*(EFP(1,JJEFP)-EFCENT(1,IIEFP))
     *                      -CONX*(EFP(2,JJEFP)-EFCENT(2,IIEFP))
         DEETORQ(4) = DEETORQ(4)+CONZ*(CCHG(2,I)-YCOM)
     *                      -CONY*(CCHG(3,I)-ZCOM)
         DEETORQ(5) = DEETORQ(5)+CONX*(CCHG(3,I)-ZCOM)
     *                      -CONZ*(CCHG(1,I)-XCOM)
         DEETORQ(6) = DEETORQ(6)+CONY*(CCHG(1,I)-XCOM)
     *                      -CONX*(CCHG(2,I)-YCOM)
C
C END LOOP OVER VALENCE LMOs
  790 CONTINUE
C  END LOOP OVER EFP POINTS
  800 CONTINUE
      IF (IDSGRD.EQ.0) GOTO 801
      DO J=1,3
         DESFRG(J,IIEFP) = DESFRG(J,IIEFP) + DEE(J)
         DESFRG(J+3,IIEFP) = DESFRG(J+3,IIEFP) + DEETORQ(J)
      ENDDO
C END LOOP OVER NUMBER OF FRAGMENTS
  801 CONTINUE
      IF (IDSGRD.EQ.0) GOTO 999
C     GET MOMENTS OF INERTIA AND TOTAL MASS FOR AI
      TMASS=ZERO
      DO I=1,NAT
         CEE(I,1)=C(1,I)
         CEE(I,2)=C(2,I)
         CEE(I,3)=C(3,I)
         ZNU(I)=ZAN(I)
         TMASS = TMASS + ZMASS(I)
      ENDDO
      CALL PRAXIS(CEE,ZNU,NAT,NAT,AIX,AIY,AIZ)
C     CALCULATE ROTATIONAL ACCELERATION AT ATOMS
      IF (NAT.EQ.1) THEN
         AROT(1,1)=ZERO
         AROT(2,1)=ZERO
         AROT(3,1)=ZERO
         GOTO 390
      END IF
      IF (AIX.GT.1.0D-05) THEN
         ALX=DEETORQ(4)/AIX
      ELSE
         ALX=ZERO
      END IF
      IF (AIY.GT.1.0D-05) THEN
         ALY=DEETORQ(5)/AIY
      ELSE
         ALY=ZERO
      END IF
      IF (AIZ.GT.1.0D-05) THEN
         ALZ=DEETORQ(6)/AIZ
      ELSE
         ALZ=ZERO
      END IF
      DO I=1,NAT
        AROT(1,I)=ALZ*(C(2,I)-YCOM)-ALY*(C(3,I)-ZCOM)
        AROT(2,I)=ALX*(C(3,I)-ZCOM)-ALZ*(C(1,I)-XCOM)
        AROT(3,I)=ALY*(C(1,I)-XCOM)-ALX*(C(2,I)-YCOM)
      ENDDO
C     CALCULATE TRANSLATIONAL ACCELERATION
 390  DO I=1,3
         ATRAN(I)=DEE(I)/TMASS
      ENDDO
C     CALCULATE TOTAL ACCELERATION AND TOTAL FORCE,
C     ADD TO GRADIENT "DE"
      DO I=1,NAT
         DO J=1,3
            ATOTAL(J,I)=AROT(J,I)+ATRAN(J)
            FTOTAL(J,I)=ATOTAL(J,I)*ZMASS(I)
            DE(J,I) = DE(J,I) + FTOTAL(J,I)
         ENDDO
      ENDDO
C
C  APPROXIMATE THE HIGHER ORDER TERMS IN DISPERSION ENERGY
C  EXPANSION AS 1/3 OF THE C6 CONTRIBUTION
  999 DISPENRG = DISP*FOURTHIRD
C     WRITE(IW,1000) DISPENRG
      IF (ALLOCATED(SABEFPAI)) DEALLOCATE(SABEFPAI)
      RETURN
C1000 FORMAT(1X,'EFP-AI DISPERSION ENERGY =',F20.10)
      END SUBROUTINE EFAIDISP
C
C *MODULE EFDRVR  *DECK EFAIMF
!> @brief Calculates imaginary frequency integral between
!>        pairs of AI orbitals; for use with EFAIDISP code
!>
!> @author Tony Smith
!> - Jan 10, 2013
!> @param FR : the final calculated value of the integral
!> @param ALP : alpha tensor; array of 12 values
!> @param ENRGDF : energy difference between two AI orbitals
      SUBROUTINE EFAIMF(FR,ALP,ENRGDF)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C THE 12 GAUSS-LEGENDRE WEIGHTS
      PARAMETER (FACT1= 0.72086099022968040154D-02)
      PARAMETER (FACT2= 0.17697067815034886394D-01)
      PARAMETER (FACT3= 0.30660908596251749739D-01)
      PARAMETER (FACT4= 0.48381293256249884995D-01)
      PARAMETER (FACT5= 0.74878830420650517080D-01)
      PARAMETER (FACT6= 0.11806515901361630228D+00)
      PARAMETER (FACT7= 0.19535413832209084204D+00)
      PARAMETER (FACT8= 0.35055692324483221824D+00)
      PARAMETER (FACT9= 0.71577113554429568336D+00)
      PARAMETER (FACT10=1.81409759976323969729D+00)
      PARAMETER (FACT11=6.97923445114870823247D+00)
      PARAMETER (FACT12=83.2480938829658453917D+00)
C THE 12 GAUSS-LEGENDRE ABSCISSAS A.K.A. IMAGINARY FREQ VALUES
C (THE FOLLOWING ARE ACTUALLY THE FREQUENCIES *SQUARED*)
      PARAMETER (FREQ1=0.779327022332536710825D-05)
      PARAMETER (FREQ2=0.228210717737242978745D-03)
      PARAMETER (FREQ3=0.152113192477780750688D-02)
      PARAMETER (FREQ4=0.608339199058554610325D-02)
      PARAMETER (FREQ5=0.192239670393041989464D-01)
      PARAMETER (FREQ6=0.543928293635942075334D-01)
      PARAMETER (FREQ7=0.148916688004125595980D+00)
      PARAMETER (FREQ8=0.421349037034822917794D+00)
      PARAMETER (FREQ9=1.331494010666300803437D+00)
      PARAMETER (FREQ10=5.32498192172462030801D+00)
      PARAMETER (FREQ11=3.54935126637048206534D+01)
      PARAMETER (FREQ12=1.03935828835455831714D+03)
C
      DIMENSION  ALP(12)
C
      ENRGSQ=ENRGDF*ENRGDF
C
      FR= ALP(1)*FACT1/(ENRGSQ+FREQ1)
     *  + ALP(2)*FACT2/(ENRGSQ+FREQ2)
     *  + ALP(3)*FACT3/(ENRGSQ+FREQ3)
     *  + ALP(4)*FACT4/(ENRGSQ+FREQ4)
     *  + ALP(5)*FACT5/(ENRGSQ+FREQ5)
     *  + ALP(6)*FACT6/(ENRGSQ+FREQ6)
     *  + ALP(7)*FACT7/(ENRGSQ+FREQ7)
     *  + ALP(8)*FACT8/(ENRGSQ+FREQ8)
     *  + ALP(9)*FACT9/(ENRGSQ+FREQ9)
     *  + ALP(10)*FACT10/(ENRGSQ+FREQ10)
     *  + ALP(11)*FACT11/(ENRGSQ+FREQ11)
     *  + ALP(12)*FACT12/(ENRGSQ+FREQ12)
      FR=FR*ENRGDF
      RETURN
      END SUBROUTINE EFAIMF
C*MODULE EFDRVR  *DECK EFPDYN_OLD
C>  @brief   EFP-EFP dispersion interaction and gradient
C>
C>  @author  Ivana Adamovic (energy), Lyudmila V. Slipchenko (gradient)?
C>  @details This version uses the old overlap damping for energy and
C>           gradient of EFP-EFP: The damping function is
C>           f_6 = 1 - S^2(1-2ln|S| + 2ln^2|S|)
C>
      SUBROUTINE EFPDYN_OLD(EDISP,EDISP0,DYNTR,SABSAVE,SABGRAD,MXPPTS)

      USE MX_LIMITS,ONLY:MXFRG,MXIFRQ,MXFGPT,MXDFG,MXDPPT
      use comm_EFPAR
      use comm_EFPPBC
      use comm_EFPGES
      use comm_EFPPAR
      use comm_EFPSCR
      USE comm_EFPFRC
      USE comm_FGRAD
      USE comm_FRGINF
      USE comm_EFPGES
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      PARAMETER (ZERO=0.0D+00)
      PARAMETER (THREE=3.0D+00)
      PARAMETER (FOUR=4.0D+00)
      PARAMETER (ONE=1.0D+00)
      PARAMETER (TWO=2.0D+00)
      PARAMETER (SIX=6.0D+00)
C      PARAMETER (B=1.5D+00)
      PARAMETER (FACT1= 0.72086099022968040154D-02)
      PARAMETER (FACT2= 0.17697067815034886394D-01)
      PARAMETER (FACT3= 0.30660908596251749739D-01)
      PARAMETER (FACT4= 0.48381293256249884995D-01)
      PARAMETER (FACT5= 0.74878830420650517080D-01)
      PARAMETER (FACT6= 0.11806515901361630228D+00)
      PARAMETER (FACT7= 0.19535413832209084204D+00)
      PARAMETER (FACT8= 0.35055692324483221824D+00)
      PARAMETER (FACT9= 0.71577113554429568336D+00)
      PARAMETER (FACT10=1.81409759976323969729D+00)
      PARAMETER (FACT11=6.97923445114870823247D+00)
      PARAMETER (FACT12=83.2480938829658453917D+00)
C
      LOGICAL GOPARR,DSKWRK,MASWRK
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
C
      DIMENSION  DYNTR(MXIFRQ*MXPPTS,NFRG),
     *           SABSAVE(NFRG,NFRG,MXPPTS,MXPPTS),
     *           SABGRAD(NFRG,NFRG,MXPPTS,MXPPTS,6),DSAB(6),
     *           IFP1(10),IFP2(10)
C
C
C     ------------------------------------------------------------------
C        THIS ROUTINE CALCULATES FRAGMENT-FRAGMENT DISPERSION (EFP2)
C        ENERGY AND GRADIENTS IN THE ABSENCE OF AN AB INITO MOLECULE
C     ------------------------------------------------------------------
C
      EDISP  = ZERO
      EDISP6  = ZERO
      EDISP0  = ZERO
      EDISP60  = ZERO
      CALL VCLR(DSAB,1,6)
C
C     -- TRACE OF DYNAMIC ALPHA IN THE FORM ALPHA(NDPPTS,NFRG) ---
C
      ID=0
      DO I=1,NFRG
        DO J=1,NDPPTS(I)
           ID=ID+1
           DYNTR(J,I) = (ONE/THREE)*(EFDPOL(1,ID) + EFDPOL(2,ID) +
     *                   EFDPOL(3,ID))
        ENDDO
      ENDDO
C
      ID1=0
      IF (GOPARR) THEN
         DO I=1,NLOW-1
            ID1 = ID1 + NPPTS(I)
         ENDDO
      ENDIF
C
      DO 280 IM = NLOW,NUP
         INPPTS=NPPTS(IM)
         COMX1 = EFCENT(1,IM)
         COMY1 = EFCENT(2,IM)
         COMZ1 = EFCENT(3,IM)
         DO 270 IN=1,NPPTS(IM)
            ID1=ID1+1
            ALP1=DYNTR(IN,IM)
            ALP2=DYNTR(IN+INPPTS,IM)
            ALP3=DYNTR(IN+2*INPPTS,IM)
            ALP4=DYNTR(IN+3*INPPTS,IM)
            ALP5=DYNTR(IN+4*INPPTS,IM)
            ALP6=DYNTR(IN+5*INPPTS,IM)
            ALP7=DYNTR(IN+6*INPPTS,IM)
            ALP8=DYNTR(IN+7*INPPTS,IM)
            ALP9=DYNTR(IN+8*INPPTS,IM)
            ALP10=DYNTR(IN+9*INPPTS,IM)
            ALP11=DYNTR(IN+10*INPPTS,IM)
            ALP12=DYNTR(IN+11*INPPTS,IM)
            X1 = EFP(1,ID1)
            Y1 = EFP(2,ID1)
            Z1 = EFP(3,ID1)
C                IF FORCE POINTS ARE ACTIVE
            IF(IFRCPNT) THEN
               CALL NEARFP(IM,X1,Y1,Z1,IFP1,IFP1_COM,NEQUAL1)
            ENDIF
C
            ID2=0
            DO 260 JM = 1,NFRG
               IF (IM.GE.JM)THEN
                 ID2=ID2+NPPTS(JM)
                 GO TO 260
               END IF
C
               COMX2 = EFCENT(1,JM)
               COMY2 = EFCENT(2,JM)
               COMZ2 = EFCENT(3,JM)
               COMX = COMX1 - COMX2
               COMY = COMY1 - COMY2
               COMZ = COMZ1 - COMZ2
               PBCX = ZERO
               PBCY = ZERO
               PBCZ = ZERO
               IF (XBOX .NE. ZERO) PBCX = XBOX * ANINT(COMX/XBOX)
               IF (YBOX .NE. ZERO) PBCY = YBOX * ANINT(COMY/YBOX)
               IF (ZBOX .NE. ZERO) PBCZ = ZBOX * ANINT(COMZ/ZBOX)
C
C   ----  DISTANCE BETWEEN FRAGMENTS
               X12  = ABS(COMX - PBCX)
               Y12  = ABS(COMY - PBCY)
               Z12  = ABS(COMZ - PBCZ)
               R12  = X12*X12+Y12*Y12+Z12*Z12
C
               SWF   = 1.0D+00
               SWFDX = 0.0D+00
               SWFDY = 0.0D+00
               SWFDZ = 0.0D+00
C
               IF(X12.GT.XCUT.OR.Y12.GT.YCUT.OR.Z12.GT.ZCUT)THEN
                 ID2=ID2+NPPTS(JM)
                 GOTO 260
               END IF
               R12  = X12*X12+Y12*Y12+Z12*Z12
               IF(R12.GT.RCUT)THEN
                 ID2=ID2+NPPTS(JM)
                 GOTO 260
               END IF
C
               IF(R12.GT.SWR2*SWR2) THEN
                 ID2=ID2+NPPTS(JM)
                 GOTO 260
               ELSE IF(R12.GT.SWR1*SWR1) THEN
                 IF(ISWF.EQ.1) THEN
                    DUMY1 = 1.0D+00/(SWR2*SWR2 - SWR1*SWR1)
                    DUMY2 = DUMY1*DUMY1
                    DUMY3 = DUMY2*DUMY1
                    DUMY4 = DUMY2*DUMY2
                    DUMY5 = DUMY2*DUMY3
                    VUMY1 = R12 - SWR1*SWR1
                    VUMY2 = VUMY1*VUMY1
                    VUMY3 = VUMY2*VUMY1
                    VUMY4 = VUMY2*VUMY2
                    VUMY5 = VUMY2*VUMY3
                    DSWF  = -30.0D+00*DUMY3*VUMY2
     *                      +60.0D+00*DUMY4*VUMY3
     *                      -30.0D+00*DUMY5*VUMY4
                    DSWF  = 2.0D+00*DSWF
                    SWF   = 1.0D+00 - 10.0D+00*DUMY3*VUMY3
     *                              + 15.0D+00*DUMY4*VUMY4
     *                              -  6.0D+00*DUMY5*VUMY5
                    SWFDX = DSWF*(EFCENT(1,IM)-EFCENT(1,JM)-PBCX)
                    SWFDY = DSWF*(EFCENT(2,IM)-EFCENT(2,JM)-PBCY)
                    SWFDZ = DSWF*(EFCENT(3,IM)-EFCENT(3,JM)-PBCZ)
                 ELSE IF(ISWF.EQ.2) THEN
                    R12   = SQRT(R12)
                    PIL   = 3.1415926535897932D+00/(SWR2-SWR1)
                    DUMY  = (R12 - SWR1) * PIL
                    SWF   = 0.50D+00 * (COS(DUMY) + 1.0D+00)
                    DUMY2 =-0.50D+00 * SIN(DUMY) * PIL / R12
                    SWFDX = DUMY2*(EFCENT(1,IM)-EFCENT(1,JM)-PBCX)
                    SWFDY = DUMY2*(EFCENT(2,IM)-EFCENT(2,JM)-PBCY)
                    SWFDZ = DUMY2*(EFCENT(3,IM)-EFCENT(3,JM)-PBCZ)
                 END IF
               END IF
C
               JNPPTS=NPPTS(JM)
               DO 250 JP = 1,NPPTS(JM)
                  ID2=ID2+1
                  ALPN1=DYNTR(JP,JM)
                  ALPN2=DYNTR(JP+JNPPTS,JM)
                  ALPN3=DYNTR(JP+2*JNPPTS,JM)
                  ALPN4=DYNTR(JP+3*JNPPTS,JM)
                  ALPN5=DYNTR(JP+4*JNPPTS,JM)
                  ALPN6=DYNTR(JP+5*JNPPTS,JM)
                  ALPN7=DYNTR(JP+6*JNPPTS,JM)
                  ALPN8=DYNTR(JP+7*JNPPTS,JM)
                  ALPN9=DYNTR(JP+8*JNPPTS,JM)
                  ALPN10=DYNTR(JP+9*JNPPTS,JM)
                  ALPN11=DYNTR(JP+10*JNPPTS,JM)
                  ALPN12=DYNTR(JP+11*JNPPTS,JM)
C
                  X2 = EFP(1,ID2)
                  Y2 = EFP(2,ID2)
                  Z2 = EFP(3,ID2)
C                   IF FORCE POINTS ARE ACTIVE
                  IF(IFRCPNT) THEN
                    CALL NEARFP(JM,X2,Y2,Z2,IFP2,IFP2_COM,NEQUAL2)
                  ENDIF
C
                  X = X1 - X2
                  Y = Y1 - Y2
                  Z = Z1 - Z2
                  X = X - PBCX
                  Y = Y - PBCY
                  Z = Z - PBCZ
C
                  R2 = X*X + Y*Y + Z*Z
                  R = SQRT(R2)
                  R4 = R2*R2
                  R6 = R4*R2
C
                  F6 = ONE
                  IF (IDISPS .EQ. 1) THEN
C                    THE OVERLAP BASED DAMPING FOR DISPERSION
                     SAB = SABSAVE(IM,JM,IN,JP)
                     IF (ABS(SAB) .GT. 1.0D-05) THEN
                        DO ICRD = 1,6
                            DSAB(ICRD) = SABGRAD(IM,JM,IN,JP,ICRD)
                        END DO
                        RB = -2.0D+00*LOG(ABS(SAB))
                        EX1 = SAB*SAB
                        RB2 = RB*RB
                        F6 = ONE-(EX1*(ONE + RB + RB2/TWO ))
                     ELSE
                        F6 = ONE
                        RB = ZERO
                        SAB = ZERO
                        DO ICRD = 1,6
                            DSAB(ICRD) = ZERO
                        END DO
                     END IF
                  ELSE
C                 - F6 IS DAMPING FUNCTION FOR DISPERSION BASED ON
C                   TANG AND TOENNIES FORMULA
                     RB =  B*R
                     EX1 = EXP(-RB)
                     RB2 = RB*RB
                     RB3 = RB*RB2
                     RB4 = RB2*RB2
                     RB5 = RB*RB4
                     RB6 = RB*RB5
C
                     F6 = ONE-(EX1*(ONE + RB + RB2/TWO +
     *                           RB3/SIX + RB4/24.0D+00 +
     *                           RB5/120.0D+00 +RB6/720.0D+00))
                  END IF
C
                  C6=
     *             + FACT1*ALP1*ALPN1
     *             + FACT2*ALP2*ALPN2
     *             + FACT3*ALP3*ALPN3
     *             + FACT4*ALP4*ALPN4
     *             + FACT5*ALP5*ALPN5
     *             + FACT6*ALP6*ALPN6
     *             + FACT7*ALP7*ALPN7
     *             + FACT8*ALP8*ALPN8
     *             + FACT9*ALP9*ALPN9
     *             + FACT10*ALP10*ALPN10
     *             + FACT11*ALP11*ALPN11
     *             + FACT12*ALP12*ALPN12
C
C                 - EDISP6 IS DAMPED DISPERSION ENERGY
C
                  EPAIR = - F6*C6/R6
                  EPAIR0 = - C6/R6
                  EDISP6 = EDISP6 + SWF*EPAIR
                  EDISP60 = EDISP60 + SWF*EPAIR0
                  EPAIR=EPAIR*FOUR/THREE
                  EPAIR0=EPAIR0*FOUR/THREE
C
C                 - GRADIENT OF THE DISPERSION ENERGY
C                 DE/DR=FOUR/THREE*(-6*C6*F6/R**7 + DF6/DR * C6/R**6)
C
                  IF (IDISPS .EQ. 1) THEN
                     CON1=-(FOUR/THREE)*(C6/R6)*(SIX*F6/R2)
                     CON2=-(FOUR/THREE)*(C6/R6)*(RB2*SAB)
                  ELSE
                     CON1=-(FOUR/THREE)*(C6/R6)*
     *                     ((SIX*F6/R2)-(B*RB6*EX1)/(R*720.0D+00))
                     CON2 = ZERO
                  END IF
C
                  CONX=CON1*X
                  CONY=CON1*Y
                  CONZ=CON1*Z
C
C
                  DESFRG(1,JM) = DESFRG(1,JM) - EPAIR*SWFDX
                  DESFRG(2,JM) = DESFRG(2,JM) - EPAIR*SWFDY
                  DESFRG(3,JM) = DESFRG(3,JM) - EPAIR*SWFDZ
                  DESFRG(1,IM) = DESFRG(1,IM) + EPAIR*SWFDX
                  DESFRG(2,IM) = DESFRG(2,IM) + EPAIR*SWFDY
                  DESFRG(3,IM) = DESFRG(3,IM) + EPAIR*SWFDZ
C
                  DESFRG(1,JM) = DESFRG(1,JM) + SWF*(CONX+CON2*DSAB(1))
                  DESFRG(2,JM) = DESFRG(2,JM) + SWF*(CONY+CON2*DSAB(2))
                  DESFRG(3,JM) = DESFRG(3,JM) + SWF*(CONZ+CON2*DSAB(3))
                  DESFRG(4,JM) = DESFRG(4,JM) + SWF*
     *                      (CONZ*(EFP(2,ID2)-EFCENT(2,JM))
     *                      -CONY*(EFP(3,ID2)-EFCENT(3,JM))
     *                      -CON2*DSAB(4)
     *                      +CON2*DSAB(3)*(COMY - PBCY)
     *                      -CON2*DSAB(2)*(COMZ - PBCZ))
                  DESFRG(5,JM) = DESFRG(5,JM) + SWF*
     *                      (CONX*(EFP(3,ID2)-EFCENT(3,JM))
     *                      -CONZ*(EFP(1,ID2)-EFCENT(1,JM))
     *                      -CON2*DSAB(5)
     *                      +CON2*DSAB(1)*(COMZ - PBCZ)
     *                      -CON2*DSAB(3)*(COMX - PBCX))
                  DESFRG(6,JM) = DESFRG(6,JM) + SWF*
     *                      (CONY*(EFP(1,ID2)-EFCENT(1,JM))
     *                      -CONX*(EFP(2,ID2)-EFCENT(2,JM))
     *                      -CON2*DSAB(6)
     *                      +CON2*DSAB(2)*(COMX - PBCX)
     *                      -CON2*DSAB(1)*(COMY - PBCY))
                  DESFRG(1,IM) = DESFRG(1,IM) - SWF*(CONX+CON2*DSAB(1))
                  DESFRG(2,IM) = DESFRG(2,IM) - SWF*(CONY+CON2*DSAB(2))
                  DESFRG(3,IM) = DESFRG(3,IM) - SWF*(CONZ+CON2*DSAB(3))
                  DESFRG(4,IM) = DESFRG(4,IM) - SWF*
     *                      (CONZ*(EFP(2,ID1)-EFCENT(2,IM))
     *                      -CONY*(EFP(3,ID1)-EFCENT(3,IM))
     *                      -CON2*DSAB(4))
                  DESFRG(5,IM) = DESFRG(5,IM) - SWF*
     *                      (CONX*(EFP(3,ID1)-EFCENT(3,IM))
     *                      -CONZ*(EFP(1,ID1)-EFCENT(1,IM))
     *                      -CON2*DSAB(5))
                  DESFRG(6,IM) = DESFRG(6,IM) - SWF*
     *                      (CONY*(EFP(1,ID1)-EFCENT(1,IM))
     *                      -CONX*(EFP(2,ID1)-EFCENT(2,IM))
     *                      -CON2*DSAB(6))
C
C  ---- GATHER GRADS AND TORQUES IN FORCE POINTS
             IF(IFRCPNT) THEN
                  FRCTRQ(1,IFP2_COM) = FRCTRQ(1,IFP2_COM)-EPAIR*SWFDX
                  FRCTRQ(2,IFP2_COM) = FRCTRQ(2,IFP2_COM)-EPAIR*SWFDY
                  FRCTRQ(3,IFP2_COM) = FRCTRQ(3,IFP2_COM)-EPAIR*SWFDZ
                  FRCTRQ(1,IFP1_COM) = FRCTRQ(1,IFP1_COM)+EPAIR*SWFDX
                  FRCTRQ(2,IFP1_COM) = FRCTRQ(2,IFP1_COM)+EPAIR*SWFDY
                  FRCTRQ(3,IFP1_COM) = FRCTRQ(3,IFP1_COM)+EPAIR*SWFDZ
C
                  DO I=1,NEQUAL2
                  IFP = IFP2(I)
                  COEF=ONE/NEQUAL2
                  FRCTRQ(1,IFP) = FRCTRQ(1,IFP)+
     *                            SWF*(CONX+CON2*DSAB(1))*COEF
                  FRCTRQ(2,IFP) = FRCTRQ(2,IFP)+
     *                            SWF*(CONY+CON2*DSAB(2))*COEF
                  FRCTRQ(3,IFP) = FRCTRQ(3,IFP)+
     *                            SWF*(CONZ+CON2*DSAB(3))*COEF
                  FRCTRQ(4,IFP) = FRCTRQ(4,IFP)+SWF*COEF*
     *                      (CONZ*(EFP(2,ID2)-FRCCRD(2,IFP))
     *                      -CONY*(EFP(3,ID2)-FRCCRD(3,IFP))
     *                      -CON2*DSAB(4)
     *                      +CON2*DSAB(3)*(COMY - PBCY)
     *                      -CON2*DSAB(2)*(COMZ - PBCZ))
                  FRCTRQ(5,IFP) = FRCTRQ(5,IFP)+SWF*COEF*
     *                      (CONX*(EFP(3,ID2)-FRCCRD(3,IFP))-
     *                       CONZ*(EFP(1,ID2)-FRCCRD(1,IFP))
     *                      -CON2*DSAB(5)
     *                      +CON2*DSAB(1)*(COMZ - PBCZ)
     *                      -CON2*DSAB(3)*(COMX - PBCX))
                  FRCTRQ(6,IFP) = FRCTRQ(6,IFP)+SWF*COEF*
     *                      (CONY*(EFP(1,ID2)-FRCCRD(1,IFP))
     *                      -CONX*(EFP(2,ID2)-FRCCRD(2,IFP))
     *                      -CON2*DSAB(6)
     *                      +CON2*DSAB(2)*(COMX - PBCX)
     *                      -CON2*DSAB(1)*(COMY - PBCY))
                END DO
                  DO I=1,NEQUAL1
                  IFP = IFP1(I)
                  COEF=ONE/NEQUAL1
                  FRCTRQ(1,IFP) = FRCTRQ(1,IFP)-
     *                            SWF*(CONX+CON2*DSAB(1))*COEF
                  FRCTRQ(2,IFP) = FRCTRQ(2,IFP)-
     *                            SWF*(CONY+CON2*DSAB(2))*COEF
                  FRCTRQ(3,IFP) = FRCTRQ(3,IFP)-
     *                            SWF*(CONZ+CON2*DSAB(3))*COEF
                  FRCTRQ(4,IFP) = FRCTRQ(4,IFP) - SWF*COEF*
     *                      (CONZ*(EFP(2,ID1)-FRCCRD(2,IFP))
     *                      -CONY*(EFP(3,ID1)-FRCCRD(3,IFP))
     *                      -CON2*DSAB(4))
                  FRCTRQ(5,IFP) = FRCTRQ(5,IFP) - SWF*COEF*
     *                      (CONX*(EFP(3,ID1)-FRCCRD(3,IFP))-
     *                       CONZ*(EFP(1,ID1)-FRCCRD(1,IFP))
     *                      -CON2*DSAB(5))
                  FRCTRQ(6,IFP) = FRCTRQ(6,IFP) - SWF*COEF*
     *                      (CONY*(EFP(1,ID1)-FRCCRD(1,IFP))
     *                      -CONX*(EFP(2,ID1)-FRCCRD(2,IFP))
     *                      -CON2*DSAB(6))
                END DO
              END IF
C
 250           CONTINUE
 260        CONTINUE
 270     CONTINUE
 280  CONTINUE
C
      IF(GOPARR) THEN
         CALL DDI_GSUMF(677,EDISP6,1)
         CALL DDI_GSUMF(677,EDISP60,1)
      END IF
      EDISP=EDISP6+EDISP6/THREE
      EDISP0=EDISP60+EDISP60/THREE
C
C
      RETURN
      END

C  *MODULE EFDRVR *DECK Delta
!> @brief Kronecker Delta function
!>
!> @author Peng Xu
!> - Nov 16, 2016
!
      FUNCTION Delta(I,J)
      USE constants, only: zero,one
      Implicit None 
      Integer,intent(in) :: I,J
      Double Precision :: Delta

      If (I .eq. J) then
        Delta = one
      else
        Delta = zero 
      Endif

      END FUNCTION Delta

! The following T_tensor_# subroutines compute the T tensors
! of various ranks (up to 5th rank) used for EFP-EFP dispersion
! The formula are Eq. (3.1.3)-(3.1.7) in "The Theory of Intermolecular
! Forces" by A.J. Stone
!  
C  *MODULE EFDRVR *DECK T_tensor_0
!> @brief Electrostatic rank 0 tensors, ~ O(R^-1)
!>
!> @author Peng Xu
!>  - Dec, 2016
!>
!> @param C: distance vector R between two points
!> @param T: rank 0 T tensor
!> 
      SUBROUTINE T_tensor_0(C,T)
      USE constants, only: one
      Implicit none 

      double precision :: C(3)
      Double precision :: T
      Double precision :: R,R2 

        R2 = C(1)*C(1)+C(2)*C(2)+C(3)*C(3)
        R = sqrt(R2)

        T = one/R

      RETURN
      END subroutine T_tensor_0

C  *MODULE EFDRVR *DECK T_tensor_1
!> @brief Electrostatic rank 1 tensors, ~ O(R^-2)
!>
!> @author Peng Xu
!>  - Dec, 2016
!>
!> @param C: distance vector R between two points
!> @param T: rank 1 T tensor
!
      SUBROUTINE T_tensor_1(C,T)
         
      Implicit none 
       
      integer :: I 
      double precision :: C(3)
      Double precision :: T(3)
      Double precision :: R,R2,R3

        R2 = C(1)*C(1)+C(2)*C(2)+C(3)*C(3)
        R = sqrt(R2)
        R3 = R*R2

        Do I = 1,3
           T(I)= - C(I)/abs(R3)
        Enddo

      RETURN
      END subroutine T_tensor_1

C  *MODULE EFDRVR *DECK T_tensor_2
!> @brief Electrostatic rank 2 tensors, ~ O(R^-3)
!>
!> @author Peng Xu
!>  - Dec, 2016
!>
!> @param C: distance vector R between two points
!> @param T: rank 2 T tensor
!
      SUBROUTINE T_tensor_2(C,T)
      USE constants, only:three
      Implicit none 

      integer :: I,J
      double precision :: C(3)
      Double precision :: T(3,3)
      Double precision :: R,R2,R3,R5
      Double precision :: delta

        R2 = C(1)*C(1)+C(2)*C(2)+C(3)*C(3)
        R = sqrt(R2)
        R3 = R*R2
        R5 = R2*R3

        Do I = 1,3
         Do J = 1,3
            T(I,J)= (three*C(I)*C(J)-R2*delta(I,J))/abs(R5)
         Enddo
        Enddo

      RETURN
      END subroutine T_tensor_2

C  *MODULE EFDRVR *DECK T_tensor_3
!> @brief Electrostatic rank 3 tensors, ~ O(R^-4)
!>
!> @author Peng Xu
!>  - Dec, 2016
!>
!> @param C: distance vector R between two points
!> @param T: rank 3 T tensor
!
      SUBROUTINE T_tensor_3(C,T)
      USE constants, only:three,fifteen
      Implicit none 
     
      integer :: I,J,K
      double precision :: C(3)
      Double precision :: T(3,3,3)
      Double precision :: R,R2,R3,R5,R7
      Double precision :: delta

        R2 = C(1)*C(1)+C(2)*C(2)+C(3)*C(3)
        R = sqrt(R2)
        R3 = R*R2
        R5 = R2*R3
        R7 = R2*R5
!!$omp parallel do private(I,J,K) schedule(dynamic)
!!$omp& collapse(3)
        Do I = 1,3
         Do J = 1,3
          Do K = 1,3
             T(I,J,K) = -(fifteen*C(I)*C(J)*C(K)
     *                    -three*R2*(C(I)*delta(J,K)+
     *                               C(J)*delta(I,K)+
     *                               C(K)*delta(I,J)))
     *                    /abs(R7)
c----PX: in EFPDYN, T3 had an extra negative sign
             T(I,J,K) = -T(I,J,K)
c----END PX
          Enddo
         Enddo
        Enddo
!!$omp end parallel do 
      RETURN
      END subroutine T_tensor_3

C  *MODULE EFDRVR *DECK T_tensor_4
!> @brief Electrostatic rank 4 tensors, ~ O(R^-5)
!>
!> @author Peng Xu
!>  - Dec, 2016
!>
!> @param C: distance vector R between two points
!> @param T: rank 4 T tensor
!
      SUBROUTINE T_tensor_4(C,T)
      USE constants,only: three,fifteen
      Implicit none 

      integer :: I,J,K,L
      double precision :: C(3)
      Double precision :: T(3,3,3,3)
      Double precision :: R,R2,R3,R4,R5,R7,R9
      Double precision :: delta

        R2 = C(1)*C(1)+C(2)*C(2)+C(3)*C(3)
        R = sqrt(R2)
        R3 = R*R2
        R4 = R2*R2
        R5 = R2*R3
        R7 = R2*R5
        R9 = R2*R7
!!$omp parallel do default(private) shared(C,R2,R4,R9,T)
!!$omp& schedule(dynamic) 
!!$omp& collapse(4)
        Do I = 1,3
         Do J = 1,3
          Do K = 1,3
           Do L = 1,3
           T(I,J,K,L) = ((105.0D+00)*C(I)*C(J)*C(K)*C(L)
     *      -fifteen*R2*(C(I)*C(J)*delta(K,L)+C(I)*C(K)*delta(J,L)+
     *                   C(I)*C(L)*delta(J,K)+C(J)*C(K)*delta(I,L)+
     *                   C(J)*C(L)*delta(I,K)+C(K)*C(L)*delta(I,J))
     *        +three*R4*(delta(I,J)*delta(K,L)+delta(I,K)*delta(J,L)
     *                  +delta(I,L)*delta(J,K)))/abs(R9)
           Enddo
          Enddo
         Enddo
        Enddo
!!$omp end parallel do 
      RETURN
      END subroutine T_tensor_4

C  *MODULE EFDRVR *DECK T_tensor_5
!> @brief Electrostatic rank 5 tensors
!>
!> @author Peng Xu
!>  - Dec, 2016
!>
!> @param C: distance vector R between two points
!> @param T: rank 5 T tensor
!
      SUBROUTINE T_tensor_5(C,T)
      USE constants,only: fifteen
      Implicit none 

      integer :: I,J,K,L,M
      double precision :: C(3)
      Double precision :: T(3,3,3,3,3)
      Double precision :: R,R2,R3,R4,R5,R7,R9,R11
      Double precision :: delta

        R2 = C(1)*C(1)+C(2)*C(2)+C(3)*C(3)
        R = sqrt(R2)
        R3 = R*R2
        R4 = R2*R2
        R5 = R2*R3
        R7 = R2*R5
        R9 = R2*R7
        R11= R2*R9
!!$omp parallel do schedule(dynamic) collapse(5)
!!$omp parallel do schedule(dynamic)
!!$omp&  default(private) shared(C,R2,R4,R11,T)
        Do I = 1,3
         Do J = 1,3
          Do K = 1,3
           Do L = 1,3
            Do M = 1,3
             T(I,J,K,L,M) = (-945.0D+00*C(I)*C(J)*C(K)*C(L)*C(M)
     *                     +105.0D+00*R2*(C(I)*C(J)*C(K)*delta(L,M)+
     *                                    C(I)*C(J)*C(L)*delta(K,M)+
     *                                    C(I)*C(J)*C(M)*delta(K,L)+
     *                                    C(I)*C(K)*C(L)*delta(J,M)+
     *                                    C(I)*C(L)*C(M)*delta(J,K)+
     *                                    C(I)*C(K)*C(M)*delta(J,L)+
     *                                    C(J)*C(K)*C(L)*delta(I,M)+
     *                                    C(J)*C(K)*C(M)*delta(I,L)+
     *                                    C(J)*C(L)*C(M)*delta(I,K)+
     *                                    C(K)*C(L)*C(M)*delta(I,J))
     *                     -fifteen*R4*(C(I)*delta(K,M)*delta(J,L)+
     *                                  C(I)*delta(J,M)*delta(K,L)+
     *                                  C(I)*delta(L,M)*delta(J,K)+
     *                                  C(J)*delta(K,M)*delta(I,L)+
     *                                  C(J)*delta(I,M)*delta(K,L)+
     *                                  C(J)*delta(L,M)*delta(K,I)+
     *                                  C(K)*delta(J,M)*delta(I,L)+
     *                                  C(K)*delta(I,M)*delta(J,L)+
     *                                  C(K)*delta(L,M)*delta(I,J)+
     *                                  C(L)*delta(I,M)*delta(J,K)+
     *                                  C(L)*delta(J,M)*delta(K,I)+
     *                                  C(L)*delta(K,M)*delta(I,J)+
     *                                  C(M)*delta(I,L)*delta(J,L)+
     *                                  C(M)*delta(J,L)*delta(K,I)+
     *                                  C(M)*delta(K,L)*delta(I,J))
     *                        )/abs(R11)
            Enddo
           Enddo
          Enddo
         Enddo
        Enddo
!!$omp end parallel do 
      RETURN
      END subroutine T_tensor_5

C*MODULE EFDRVR  *DECK EFPDYN_MolPol 
      SUBROUTINE EFPDYN_MolPol(FACT)
      use constants
      use mx_limits, only: MXIFRQ
      use DYNPOL_MOL
      use comm_FRGINF, only: NFRG
      use comm_EFPAR, only: NLOW,NUP
      use comm_FGRAD, only: EFCENT

      IMPLICIT NONE 

      LOGICAL GOPARR,DSKWRK,MASWRK
      integer :: ME,MASTER,NPROC,IBTYP,IPTIM
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK

      integer :: Ifrag,Ifreq,I,J,K,L,IM,JM,IA,IB,IC,ID,IE,Ii,Ij,Ik
      double precision :: delta
      double precision :: E6,E6_iso,E7,E8,E8_iso
      double precision :: CC6,CC8,COMX,COMY,COMZ,
     &                    COMX1,COMX2,COMY1,COMY2,COMZ1,COMZ2
      double precision :: DUM,DUM1,DUM2,DUM3,TMP1,TMP2,TMP3,TMP4,
     &                    TMP5,TMP6,TMP7,TMP8
      double precision :: RR,RR2,RR4,RR5,RR6,RR7,RR8
      double precision :: FACT(12),CRD(3)
      double precision :: DYNDD_AVE(MXIFRQ,NFRG),DYNQQ_AVE(MXIFRQ,NFRG)
      double precision :: TT2(3,3),TT3(3,3,3),TT4(3,3,3,3)

      E6     = ZERO
      E6_iso = ZERO
      E7     = ZERO
      E8     = ZERO
      E8_iso = ZERO

!!$omp parallel private(I,J,K,L,Ifrag,Ifreq) 
!!$omp& shared(DYNDD_ROT,DYNDD_AVE,DYNQQ_ROT,DYNQQ_AVE)

C     -- trace of dynamic molecular dip-dip polarizability --
C        used for calculating isotropic molecular E6

      call vclr(DYNDD_AVE,1,MXIFRQ*NFRG)
!$omp parallel private(I,K) 
!$omp& shared(DYNDD_ROT,DYNDD_AVE)
!$omp do schedule(dynamic) 
!      call vclr(DYNDD_AVE,1,MXIFRQ*NFRG)
      DO I = 1, NFRG
        DO K = 1, MXIFRQ
         DYNDD_AVE(K,I) = (ONE/THREE)*
     *   (DYNDD_ROT(1,1,K,I)+DYNDD_ROT(2,2,K,I)+DYNDD_ROT(3,3,K,I))
        END DO
      END DO
!$omp end do
!$omp end parallel

C    -- spherical molecular quad-quad polarizability --
C       used for calculating spherical molecular E8
      call vclr(DYNQQ_AVE,1,MXIFRQ*NFRG)
!$omp parallel private(Ifrag,Ifreq,L,K,J,I)
!$omp& shared(DYNQQ_ROT,DYNQQ_AVE)
!$omp do schedule(dynamic)
      Do Ifrag = 1, NFRG
       Do Ifreq = 1, MXIFRQ
          Do L = 1,3
          Do K = 1,3
          Do J = 1,3
          Do I = 1,3
             DYNQQ_AVE(Ifreq,Ifrag)=DYNQQ_AVE(Ifreq,Ifrag) + 
     &     (One/Five)*DYNQQ_ROT(I,J,K,L,Ifreq,Ifrag)
     &    *((One/Two)*(delta(I,K)*delta(J,L)+delta(I,L)*delta(J,K))
     &    -(One/Three)*delta(I,J)*delta(K,L))
          Enddo
          Enddo
          Enddo
          Enddo
       Enddo
      Enddo
!$omp end do
!$omp end parallel
C

      DO 95 IM = NLOW,NUP
         COMX1 = EFCENT(1,IM)
         COMY1 = EFCENT(2,IM)
         COMZ1 = EFCENT(3,IM)

         DO 90 JM = 1,NFRG
               IF (IM.GE.JM)THEN
                  CYCLE
               END IF
C
               COMX2 = EFCENT(1,JM)
               COMY2 = EFCENT(2,JM)
               COMZ2 = EFCENT(3,JM)

               COMX = COMX1 - COMX2
               COMY = COMY1 - COMY2
               COMZ = COMZ1 - COMZ2

C   ---- calculate T tensors ---- 
C   Note: these are calculated using centers of mass of the fragment
C
                CRD(1) = COMX
                CRD(2) = COMY
                CRD(3) = COMZ
                RR2=COMX*COMX+COMY*COMY+COMZ*COMZ
                RR=SQRT(RR2)
                RR4=RR2*RR2
                RR5 = RR4*RR
                RR6 = RR4*RR2
                RR7 = RR5*RR2
                RR8 = RR6*RR2
                call T_tensor_2(CRD,TT2)
                call T_tensor_3(CRD,TT3)
                call T_tensor_4(CRD,TT4)
    
C    -- the anisotropic E6 from molecular polarizability --
!-----PX
!$omp parallel do default(private)
!$omp& firstprivate(IM,JM)
!$omp& shared(TT2,FACT,DYNDD_ROT) reduction(+: E6)  
!$omp& schedule(dynamic) 

               DO 80 IFREQ = 1,MXIFRQ
                DO IA = 1,3
                 DO IB = 1,3
                  DO IC = 1,3
                   DO ID = 1,3
                      DUM = TT2(IA,IB)*TT2(IC,ID)*FACT(IFREQ)
     *                     *DYNDD_ROT(IA,IC,IFREQ,IM)
     *                     *DYNDD_ROT(IB,ID,IFREQ,JM)
                      E6 = E6 - DUM/(TWO*PI)
                   ENDDO
                  ENDDO
                 ENDDO
                ENDDO
 80            CONTINUE
!$omp end parallel do
!-----END PX
!               E6 = E6/(TWO*PI)

C    P. XU : the isotropic E6 from the molecular polarizability

               CC6 = zero
!$omp parallel do schedule(dynamic) default(private)
!$omp& shared(FACT,DYNDD_AVE)
!$omp& firstprivate(IM,JM) reduction(+: CC6)
               DO 70 IFREQ = 1, MXIFRQ
                 CC6 = CC6+(THREE/PI)*FACT(IFREQ)*
     *                DYNDD_AVE(IFREQ,IM)*DYNDD_AVE(IFREQ,JM)
 70            CONTINUE
!$omp end parallel do
               E6_iso = E6_iso - CC6/RR6


C    P. XU:    E7 from molecular polarizability
!-----PX
!$omp parallel do default(private) 
!$omp& firstprivate(IM,JM)
!$omp& shared(TT2,TT3,FACT,DYNDD_ROT,DYNDQ_ROT) 
!$omp& reduction(+:E7) 
!$omp& schedule(dynamic) 
              DO 60 IFREQ = 1,MXIFRQ
                DO IA = 1,3
                 DO IB = 1,3
                  DO IC = 1,3
                   DO ID = 1,3
                    DO IE = 1,3
                       DUM1 = DYNDD_ROT(IA,IC,IFREQ,IM)*
     *                        DYNDQ_ROT(IB,ID,IE,IFREQ,JM)
                       DUM2 = DYNDD_ROT(IB,IE,IFREQ,JM)*
     *                        DYNDQ_ROT(IA,IC,ID,IFREQ,IM)
                       DUM3 = TT2(IA,IB)*TT3(IC,ID,IE)*
     *                        FACT(IFREQ)*(DUM1-DUM2)
                       E7= E7-DUM3/(THREE*PI)
                    ENDDO
                   ENDDO
                  ENDDO
                 ENDDO
                ENDDO
 60           CONTINUE
!$omp end parallel do 
!-----END PX

C      -- anisotropic E8 from molecular polarizability --
!-----PX
!$omp parallel do default(private)
!$omp& firstprivate(IM,JM)
!$omp& shared(TT2,TT3,TT4,FACT,
!$omp$   DYNDD_ROT,DYNDQ_ROT,DYNDO_ROT,DYNQQ_ROT) 
!$omp& reduction(+:E8) 
!$omp& schedule(dynamic) 
              DO 50 IFREQ = 1, MXIFRQ
               DO IA = 1,3
                DO IB = 1,3
                 DO IC = 1,3
                  DO Ii = 1,3
                   DO Ij = 1,3
                    DO Ik = 1,3
                       TMP1 = DYNDD_ROT(Ia,Ic,Ifreq,IM)
     &                       *DYNDO_ROT(Ib,Ii,Ij,Ik,Ifreq,JM)
                       TMP2 = DYNDQ_ROT(Ia,Ic,Ii,Ifreq,IM)
     &                       *DYNDQ_ROT(Ib,Ij,Ik,Ifreq,JM)
                       TMP3 = DYNDO_ROT(Ia,Ic,Ii,Ij,Ifreq,IM)
     &                       *DYNDD_ROT(Ib,Ik,Ifreq,JM)
                       TMP4 = TT2(Ia,Ib)*TT4(Ic,Ii,Ij,Ik)*FACT(Ifreq)
     &                       *(-TMP1/five + TMP2/three -TMP3/five)
                       TMP5 = DYNDD_ROT(Ia,Ii,Ifreq,IM)
     &                       *DYNQQ_ROT(Ib,Ic,Ij,Ik,Ifreq,JM)
                       TMP6 = DYNDQ_ROT(Ii,Ia,Ib,Ifreq,IM)
     &                       *DYNDQ_ROT(Ic,Ij,Ik,Ifreq,JM)
                       TMP7 = DYNQQ_ROT(Ia,Ib,Ii,Ij,Ifreq,IM)
     &                       *DYNDD_ROT(Ic,Ik,Ifreq,JM)
                       TMP8 = TT3(Ia,Ib,Ic)*TT3(Ii,Ij,Ik)*FACT(Ifreq)
     &                       *(-TMP5/two + TMP6/three - TMP7/two)
                       E8 = E8 + (TMP4+TMP8)/(three*pi)
                    Enddo ! Ik
                   Enddo ! Ij
                  Enddo ! Ii
                 Enddo ! Ic
                Enddo ! Ib
               Enddo ! Ia   
 50           CONTINUE
!$omp end parallel do 
!-----END PX

c    -- Spherical E8 using molecular polarizability tensors --
        CC8 = zero
!$omp parallel do schedule(dynamic) default(private)
!$omp& shared(FACT,DYNDD_AVE,DYNQQ_AVE)
!$omp& firstprivate(IM,JM) reduction(+: CC8)
        Do IFREQ = 1, MXIFRQ
           CC8 = CC8 + (fifteen/pi)*FACT(Ifreq)
     &    *(DYNDD_AVE(Ifreq,IM)*DYNQQ_AVE(Ifreq,JM)+
     &      DYNDD_AVE(Ifreq,JM)*DYNQQ_AVE(Ifreq,IM))
        Enddo
!$omp end parallel do 
        E8_iso = E8_iso - CC8/RR8

  90    CONTINUE
  95   CONTINUE

      IF(GOPARR) THEN
         CALL DDI_GSUMF(670,E6,1)
         CALL DDI_GSUMF(670,E6_iso,1)
         CALL DDI_GSUMF(670,E7,1)
         CALL DDI_GSUMF(670,E8,1)
         CALL DDI_GSUMF(670,E8_iso,1)
      ENDIF

      if(maswrk) then
         write(6,*) 'anisotropic molecular E6=',E6
         write(6,*) 'isotropic molecular E6=',E6_iso
         write(6,*) 'anisotropic molecular E7=',E7
         write(6,*) 'anisotropic molecular E8=',E8
         write(6,*) 'isotropic molecular E8=',E8_iso
      endif

      RETURN
      END 


C*MODULE EFDRVR  *DECK Disp6_LMOpol
      SUBROUTINE Disp6_LMOpol(IN,JP,INPPTS,JNPPTS,IM,JM,FACT,
     *                        FACT1,FACT2,FACT3,FACT4,FACT5,FACT6,
     *                        FACT7,FACT8,FACT9,FACT10,FACT11,FACT12,
     *                        ALP1,ALP2,ALP3,ALP4,ALP5,ALP6,ALP7,
     *                        ALP8,ALP9,ALP10,ALP11,ALP12,ALPN1,
     *                        ALPN2,ALPN3,ALPN4,ALPN5,ALPN6,ALPN7,
     *                        ALPN8,ALPN9,ALPN10,ALPN11,ALPN12,
     *                        F6,R6,T2,C6,EPAIR60,EPAIR6,
     *                        EPAIR60_aniso,EPAIR6_aniso,do_aniso_disp)

      use constants, only:zero,two,three,pi
      use DYNPOL_DIST, only: DYNDD_LMO_ROT

      IMPLICIT NONE 

      LOGICAL do_aniso_disp

      integer :: IN,JP,INPPTS,JNPPTS,IM,JM,IA,IB,IC,ID,IFREQ
      double precision :: FACT1,FACT2,FACT3,FACT4,FACT5,FACT6,
     &                    FACT7,FACT8,FACT9,FACT10,FACT11,FACT12
      double precision :: ALP1,ALP2,ALP3,ALP4,ALP5,ALP6,ALP7,ALP8,
     &                    ALP9,ALP10,ALP11,ALP12,ALPN1,ALPN2,ALPN3,
     &                    ALPN4,ALPN5,ALPN6,ALPN7,ALPN8,ALPN9,
     &                    ALPN10,ALPN11,ALPN12
      double precision :: F6,R6,C6,EPAIR60,EPAIR6,EPAIR60_aniso,
     &                    EPAIR6_aniso,DUM
      double precision :: T2(3,3),FACT(12)
C C6 isotropic coefficient
C - P. XU: correction: add the factor of 3/pi
C
       EPAIR60 = ZERO
       EPAIR6  = ZERO
       C6=((THREE/PI))*(
     *   + FACT1*ALP1*ALPN1
     *   + FACT2*ALP2*ALPN2
     *   + FACT3*ALP3*ALPN3
     *   + FACT4*ALP4*ALPN4
     *   + FACT5*ALP5*ALPN5
     *   + FACT6*ALP6*ALPN6
     *   + FACT7*ALP7*ALPN7
     *   + FACT8*ALP8*ALPN8
     *   + FACT9*ALP9*ALPN9
     *   + FACT10*ALP10*ALPN10
     *   + FACT11*ALP11*ALPN11
     *   + FACT12*ALP12*ALPN12)

C
       EPAIR60 = -C6/R6
       EPAIR6  = -F6*C6/R6
C --P. XU:  anisotropic E6 from distributed LMO polarizability
C
      IF(do_aniso_disp) THEN
       EPAIR60_aniso = ZERO
       EPAIR6_aniso  = ZERO
C
!$omp parallel do default(private)
!$omp&   firstprivate(IM,JM,IN,JP,INPPTS,JNPPTS) 
!$omp$   shared(T2,FACT,DYNDD_LMO_ROT) 
!$omp&   reduction(+:EPAIR60_aniso) 
!$omp&   schedule(dynamic) 
       DO 10 IFREQ = 1,12
        DO IA = 1,3
         DO IB = 1,3
          DO IC = 1,3
           DO ID = 1,3
              DUM = T2(IA,IB)*T2(IC,ID)*FACT(IFREQ)
     *        *DYNDD_LMO_ROT(IA,IC,(IN+(IFREQ-1)*INPPTS),IM)
     *        *DYNDD_LMO_ROT(IB,ID,(JP+(IFREQ-1)*JNPPTS),JM)
              EPAIR60_aniso = EPAIR60_aniso - DUM/(TWO*PI)
           ENDDO
          ENDDO
         ENDDO
        ENDDO
 10    CONTINUE
!$omp end parallel do
       EPAIR6_aniso = F6*EPAIR60_aniso
      END IF

      RETURN
      END


C*MODULE EFDRVR  *DECK Disp7_LMOpol
      SUBROUTINE Disp7_LMOpol(IN,JP,INPPTS,JNPPTS,IM,JM,IP1,FACT,
     *                        F7,T2,T3,COORD,R2,R4,R7,R9,
     *                        EPAIR70,EPAIR7,IFGRAD,
     *                        DER_EPAIR70,DER_EPAIR70_IJ)

      use constants, only: zero,one,three,fifteen,pi
      use DYNPOL_DIST, only: DER_DYNDD,DER_DYNDQ,
     *                       DYNDD_LMO_ROT,DYNDQ_LMO_ROT
      use comm_EFPPAR, only: EFP
      use comm_FGRAD, only: EFCENT

      IMPLICIT NONE 

      LOGICAL IFGRAD 

      integer :: IN,JP,INPPTS,JNPPTS,IM,JM,IP1
      integer :: IA,IB,IC,ID,IE,IAB,IABC,IAC,IACD,IBC,IBDE,IBE,ICD,
     &           ICDE,IDE,IFREQ,KRON1A,KRON1B,KRON1C,KRON2D,KRON2E,
     &           KRON2F,KRON3G,KRON3H,KRON3I,KRONA,KRONB,KRONC,KROND,
     &           KRONE,KRONF

      double precision :: F7,DUM1,DUM2,delta
      double precision :: R2,R4,R7,R9
      double precision :: EPAIR70,EPAIR7
      double precision :: T2(3,3),T3(3,3,3),COORD(3),FACT(12)
      double precision :: DUM5(3),DUM6(3),DUM7(3),DUM8(3),
     *          DER_EPAIR70(3),DER_EPAIR7(3),
     *          DER_EPAIR7_IJ(3),
     *          DUM3A(3),DER_T2X(9),DER_T3X(27),DER_T2Y(9),
     *          DER_T3Y(27),DER_T2Z(9),DER_T3Z(27),
     *          DER_T2TX(9),DER_T2TY(9),DER_T2TZ(9),
     *          DER_T3TX(27),DER_T3TY(27),DER_T3TZ(27),
     *          DER_EPAIR70_IJ(3),
     *          DER_T2TX_IJ(9), DER_T2TY_IJ(9), DER_T2TZ_IJ(9),
     *          DER_T3TX_IJ(27), DER_T3TY_IJ(27), DER_T3TZ_IJ(27)

C  ---- EPAIR70 is undamped E7 between pairs of LMOs

      EPAIR70 = ZERO
      EPAIR7  = ZERO
      CALL VCLR(DER_EPAIR70,1,3)
      CALL VCLR(DER_EPAIR70_IJ,1,3)
      CALL VCLR(DER_T2TX_IJ,1,9)
      CALL VCLR(DER_T2TY_IJ,1,9)
      CALL VCLR(DER_T2TZ_IJ,1,9)
      CALL VCLR(DER_T3TX_IJ,1,27)
      CALL VCLR(DER_T3TY_IJ,1,27)
      CALL VCLR(DER_T3TZ_IJ,1,27)
!----PX
!$omp parallel do default(private) 
!$omp$   shared(T2,T3,FACT,DYNDD_LMO_ROT,DYNDQ_LMO_ROT) 
!$omp&   firstprivate(IM,JM,IN,JP,INPPTS,JNPPTS)
!$omp&   reduction(+:EPAIR70) 
!$omp&   schedule(dynamic) 
!!collapse(6) 

       DO 202 IFREQ=1,12
         DO 182 IA = 1,3
         DO 162 IB = 1,3
         DO 142 IC = 1,3
         DO 122 ID = 1,3
         DO 102 IE = 1,3

         DUM1= DYNDD_LMO_ROT(IA,IC,(IN+(IFREQ-1)*INPPTS),IM)*
     *         DYNDQ_LMO_ROT(IB,ID,IE,(JP+(IFREQ-1)*JNPPTS),JM)

         DUM2= DYNDD_LMO_ROT(IB,IE,(JP+(IFREQ-1)*JNPPTS),JM)*
     *         DYNDQ_LMO_ROT(IA,IC,ID,(IN+(IFREQ-1)*INPPTS),IM)

         EPAIR70=EPAIR70+((-ONE/(THREE*PI))*T2(IA,IB)*T3(IC,ID,IE)
     *                     *FACT(IFREQ)*(DUM1-DUM2))
 102     CONTINUE
 122     CONTINUE
 142     CONTINUE
 162     CONTINUE
 182     CONTINUE
 202   CONTINUE
!$omp end parallel do

       IF(IFGRAD) THEN
!$omp parallel do schedule(dynamic) default(private)
!$omp& shared(T2,T3,FACT,DYNDD_LMO_ROT,DYNDQ_LMO_ROT,EFP,EFCENT,
!$omp&        DER_DYNDD,DER_DYNDQ,COORD,R2,R4,R7,R9,
!$omp&        DER_T2X,DER_T2Y,DER_T2Z,DER_T3X,DER_T3Y,DER_T3Z)
!$omp& firstprivate(IM,JM,IN,JP,INPPTS,JNPPTS,IP1) 
!$omp& reduction(+:DER_EPAIR70,DER_EPAIR70_IJ)
       DO 200 IFREQ=1,12
         DO 180 IA = 1,3
         DO 160 IB = 1,3
                IAB = (IA-1)*3+IB
         DO 140 IC = 1,3
                IAC = (IA-1)*3+IC
                IBC=(IB-1)*3+IC
                IABC=(IA-1)*9+IBC
         DO 120 ID = 1,3
                ICD = (IC-1)*3+ID
                IACD = (IA-1)*9+ICD
         DO 100 IE = 1,3
                IBE = (IB-1)*3+IE
                IDE = (ID-1)*3+IE
                IBDE = (IB-1)*9+IDE
                ICDE = (IC-1)*9+IDE

CE.Guidez T2 and T3 tensors and their translational and rotationl
Cderivatives
               IF (COORD(1).EQ.COORD(IA)) THEN
                   KRONA=1
               ELSE
                   KRONA=0
               ENDIF
               IF (COORD(1).EQ.COORD(IB)) THEN
                   KRONB=1
               ELSE
                   KRONB=0
               ENDIF
               IF (COORD(2).EQ.COORD(IA)) THEN
                   KRONC=1
               ELSE
                   KRONC=0
               ENDIF
               IF (COORD(2).EQ.COORD(IB)) THEN
                   KROND=1
               ELSE
                   KROND=0
               ENDIF
               IF (COORD(3).EQ.COORD(IA)) THEN
                   KRONE=1
               ELSE
                   KRONE=0
               ENDIF
               IF (COORD(3).EQ.COORD(IB)) THEN
                   KRONF=1
               ELSE
                   KRONF=0
               ENDIF

C Translational derivatives
               DER_T2X(IAB)=(-(FIFTEEN*COORD(1)*COORD(IA)*COORD(IB)
     *                       -THREE*R2*(COORD(1)*delta(IA,IB)
     *                       +COORD(IA)*KRONB+COORD(IB)*KRONA))
     *                       /ABS(R7))
               DER_T2Y(IAB)=(-(FIFTEEN*COORD(2)*COORD(IA)*COORD(IB)
     *                       -THREE*R2*(COORD(2)*delta(IA,IB)
     *                       +COORD(IA)*KROND+COORD(IB)*KRONC))
     *                       /ABS(R7))
               DER_T2Z(IAB)=(-(FIFTEEN*COORD(3)*COORD(IA)*COORD(IB)
     *                       -THREE*R2*(COORD(3)*delta(IA,IB)
     *                       +COORD(IA)*KRONF+COORD(IB)*KRONE))
     *                       /ABS(R7))


               IF (COORD(1).EQ. COORD(IE)) THEN
                   KRON1A=1
               ELSE
                   KRON1A=0
               ENDIF
               IF (COORD(1).EQ. COORD(ID)) THEN
                   KRON1B=1
               ELSE
                   KRON1B=0
               ENDIF
               IF (COORD(1).EQ. COORD(IC)) THEN
                   KRON1C=1
               ELSE
                   KRON1C=0
               ENDIF
               IF (COORD(2).EQ. COORD(IE)) THEN
                   KRON2D=1
               ELSE
                   KRON2D=0
               ENDIF
               IF (COORD(2).EQ. COORD(ID)) THEN
                   KRON2E=1
               ELSE
                   KRON2E=0
               ENDIF
               IF (COORD(2).EQ.COORD(IC)) THEN
                   KRON2F=1
               ELSE
                   KRON2F=0
               ENDIF
               IF (COORD(3).EQ.COORD(IE)) THEN
                   KRON3G=1
               ELSE
                   KRON3G=0
               ENDIF
               IF (COORD(3).EQ.COORD(ID)) THEN
                   KRON3H=1
               ELSE
                   KRON3H=0
               ENDIF
               IF (COORD(3).EQ.COORD(IC)) THEN
                   KRON3I=1
               ELSE
                   KRON3I=0
               ENDIF
        DER_T3X(ICDE)=-(105.0D+00*COORD(1)*COORD(IC)*COORD(ID)*COORD(IE)
     *               -FIFTEEN*R2*(COORD(1)*COORD(IC)*delta(ID,IE)+
     *               COORD(1)*COORD(ID)*delta(IC,IE)+COORD(1)*COORD(IE)*
     *               delta(IC,ID)+COORD(IC)*COORD(ID)*KRON1A+COORD(IC)*
     *               COORD(IE)*KRON1B+COORD(ID)*COORD(IE)*KRON1C)+
     *               THREE*R4*(KRON1C*delta(ID,IE)+KRON1B*delta(IC,IE)+
     *               KRON1A*delta(IC,ID)))/ABS(R9)

       DER_T3Y(ICDE)=-(105.0D+00*COORD(2)*COORD(IC)*COORD(ID)*COORD(IE)
     *               -FIFTEEN*R2*(COORD(2)*COORD(IC)*delta(ID,IE)+
     *               COORD(2)*COORD(ID)*delta(IC,IE)+COORD(2)*COORD(IE)*
     *               delta(IC,ID)+COORD(IC)*COORD(ID)*KRON2D+COORD(IC)*
     *               COORD(IE)*KRON2E+COORD(ID)*COORD(IE)*KRON2F)+
     *               THREE*R4*(KRON2F*delta(ID,IE)+KRON2E*delta(IC,IE)+
     *               KRON2D*delta(IC,ID)))/ABS(R9)

       DER_T3Z(ICDE)=-(105.0D+00*COORD(3)*COORD(IC)*COORD(ID)*COORD(IE)
     *               -FIFTEEN*R2*(COORD(3)*COORD(IC)*delta(ID,IE)+
     *               COORD(3)*COORD(ID)*delta(IC,IE)+COORD(3)*COORD(IE)*
     *               delta(IC,ID)+COORD(IC)*COORD(ID)*KRON3G+COORD(IC)*
     *               COORD(IE)*KRON3H+COORD(ID)*COORD(IE)*KRON3I)+
     *               THREE*R4*(KRON3I*delta(ID,IE)+KRON3H*delta(IC,IE)+
     *               KRON3G*delta(IC,ID)))/ABS(R9)

                DUM1= DYNDD_LMO_ROT(IA,IC,(IN+(IFREQ-1)*INPPTS),IM)*
     *                DYNDQ_LMO_ROT(IB,ID,IE,(JP+(IFREQ-1)*JNPPTS),JM)
                DUM2= DYNDD_LMO_ROT(IB,IE,(JP+(IFREQ-1)*JNPPTS),JM)*
     *                DYNDQ_LMO_ROT(IA,IC,ID,(IN+(IFREQ-1)*INPPTS),IM)

C    Translation
               DUM3A(1)=-(ONE/(THREE*PI))*
     *             ((DER_T2X(IAB)*T3(IC,ID,IE))
     *              +(T2(IA,IB)*DER_T3X(ICDE)))
     *                  *FACT(IFREQ)*(DUM1-DUM2)
               DUM3A(2)=-(ONE/(THREE*PI))*
     *             ((DER_T2Y(IAB)*T3(IC,ID,IE))
     *             +(T2(IA,IB)*DER_T3Y(ICDE)))
     *                   *FACT(IFREQ)*(DUM1-DUM2)
               DUM3A(3)=-(ONE/(THREE*PI))*
     *             ((DER_T2Z(IAB)*T3(IC,ID,IE))
     *             +(T2(IA,IB)*DER_T3Z(ICDE)))
     *                    *FACT(IFREQ)*(DUM1-DUM2)
C    Rotation
               DUM5(1)=DER_DYNDD(1,IAC,(IN+(IFREQ-1)*INPPTS),IM)*
     *             DYNDQ_LMO_ROT(IB,ID,IE,(JP+(IFREQ-1)*JNPPTS),JM)
               DUM5(2)=DER_DYNDD(2,IAC,(IN+(IFREQ-1)*INPPTS),IM)*
     *             DYNDQ_LMO_ROT(IB,ID,IE,(JP+(IFREQ-1)*JNPPTS),JM)
               DUM5(3)=DER_DYNDD(3,IAC,(IN+(IFREQ-1)*INPPTS),IM)*
     *             DYNDQ_LMO_ROT(IB,ID,IE,(JP+(IFREQ-1)*JNPPTS),JM)
               DUM6(1)=DER_DYNDQ(1,IBDE,(JP+(IFREQ-1)*JNPPTS),JM)*
     *              DYNDD_LMO_ROT(IA,IC,(IN+(IFREQ-1)*INPPTS),IM)
               DUM6(2)=DER_DYNDQ(2,IBDE,(JP+(IFREQ-1)*JNPPTS),JM)*
     *                 DYNDD_LMO_ROT(IA,IC,(IN+(IFREQ-1)*INPPTS),IM)
               DUM6(3)=DER_DYNDQ(3,IBDE,(JP+(IFREQ-1)*JNPPTS),JM)*
     *                 DYNDD_LMO_ROT(IA,IC,(IN+(IFREQ-1)*INPPTS),IM)
               DUM7(1)=DER_DYNDD(1,IBE,(JP+(IFREQ-1)*JNPPTS),JM)*
     *                DYNDQ_LMO_ROT(IA,IC,ID,(IN+(IFREQ-1)*INPPTS),IM)
               DUM7(2)=DER_DYNDD(2,IBE,(JP+(IFREQ-1)*JNPPTS),JM)*
     *                DYNDQ_LMO_ROT(IA,IC,ID,(IN+(IFREQ-1)*INPPTS),IM)
               DUM7(3)=DER_DYNDD(3,IBE,(JP+(IFREQ-1)*JNPPTS),JM)*
     *                DYNDQ_LMO_ROT(IA,IC,ID,(IN+(IFREQ-1)*INPPTS),IM)
               DUM8(1)=DER_DYNDQ(1,IACD,(IN+(IFREQ-1)*INPPTS),IM)*
     *                DYNDD_LMO_ROT(IB,IE,(JP+(IFREQ-1)*JNPPTS),JM)
               DUM8(2)=DER_DYNDQ(2,IACD,(IN+(IFREQ-1)*INPPTS),IM)*
     *                DYNDD_LMO_ROT(IB,IE,(JP+(IFREQ-1)*JNPPTS),JM)
               DUM8(3)=DER_DYNDQ(3,IACD,(IN+(IFREQ-1)*INPPTS),IM)*
     *                DYNDD_LMO_ROT(IB,IE,(JP+(IFREQ-1)*JNPPTS),JM)
c
               DER_EPAIR70(1) =DER_EPAIR70(1)+DUM3A(1)
               DER_EPAIR70(2) =DER_EPAIR70(2)+DUM3A(2)
               DER_EPAIR70(3) =DER_EPAIR70(3)+DUM3A(3)

C T2 and T3 torques
               DER_T2TX_IJ(IAB)=-(DER_T2Z(IAB)*(
     *         EFP(2,IP1)-EFCENT(2,IM)))+
     *        (DER_T2Y(IAB)*(EFP(3,IP1)-EFCENT(3,IM)))

               DER_T2TY_IJ(IAB)=-(DER_T2X(IAB)*(
     *         EFP(3,IP1)-EFCENT(3,IM)))+
     *        (DER_T2Z(IAB)*(EFP(1,IP1)-EFCENT(1,IM)))

               DER_T2TZ_IJ(IAB)=-(DER_T2Y(IAB)*(
     *         EFP(1,IP1)-EFCENT(1,IM)))+
     *        (DER_T2X(IAB)*(EFP(2,IP1)-EFCENT(2,IM)))

               DER_T3TX_IJ(ICDE)=-(DER_T3Z(ICDE)*(
     *         EFP(2,IP1)-EFCENT(2,IM)))+
     *        (DER_T3Y(ICDE)*(EFP(3,IP1)-EFCENT(3,IM)))

               DER_T3TY_IJ(ICDE)=-(DER_T3X(ICDE)*(
     *         EFP(3,IP1)-EFCENT(3,IM)))+
     *        (DER_T3Z(ICDE)*(EFP(1,IP1)-EFCENT(1,IM)))

               DER_T3TZ_IJ(ICDE)=-(DER_T3Y(ICDE)*(
     *         EFP(1,IP1)-EFCENT(1,IM)))+
     *        (DER_T3X(ICDE)*(EFP(2,IP1)-EFCENT(2,IM)))

C E70 derivative

               DER_EPAIR70_IJ(1)=DER_EPAIR70_IJ(1)
     *          +(((-ONE/(THREE*PI))*FACT(IFREQ)*((
     *           (((DER_T2TX_IJ(IAB))*T3(IC,ID,IE))
     *          +(T2(IA,IB)*(DER_T3TX_IJ(ICDE))))*
     *           (DUM1-DUM2))+(T2(IA,IB)*T3(IC,ID,IE)*(
     *           (DUM5(1)-DUM8(1)))))))

               DER_EPAIR70_IJ(2)=DER_EPAIR70_IJ(2)
     *          +(((-ONE/(THREE*PI))*FACT(IFREQ)*((
     *           ((DER_T2TY_IJ(IAB)*T3(IC,ID,IE))
     *          +(T2(IA,IB)*(DER_T3TY_IJ(ICDE))))*
     *           (DUM1-DUM2))+(T2(IA,IB)*T3(IC,ID,IE)*(
     *           (DUM5(2)-DUM8(2)))))))

               DER_EPAIR70_IJ(3)=DER_EPAIR70_IJ(3)
     *          +(((-ONE/(THREE*PI))*FACT(IFREQ)*((
     *           ((DER_T2TZ_IJ(IAB)*T3(IC,ID,IE))
     *          +(T2(IA,IB)*(DER_T3TZ_IJ(ICDE))))*
     *           (DUM1-DUM2))+(T2(IA,IB)*T3(IC,ID,IE)*(
     *           (DUM5(3)-DUM8(3)))))))

 100     CONTINUE
 120     CONTINUE
 140     CONTINUE
 160     CONTINUE
 180     CONTINUE
 200   CONTINUE
!$omp end parallel do
       ENDIF ! IFGRAD
!----END PX
       EPAIR7  = F7*EPAIR70

       RETURN
       END


C*MODULE EFDRVR  *DECK Disp8_LMOpol
      SUBROUTINE Disp8_LMOpol(IN,JP,INPPTS,JNPPTS,IM,JM,FACT,
     *                        F8,R8,T2,T3,T4,DYNQQ_LMO_AVE,
     *                        DYNTR,MXPPTS,EPAIR80,EPAIR8,
     *                        EPAIR80_aniso,EPAIR8_aniso,C8,
     *                        do_aniso_disp,E8DO,E8DQ,E8QQ)
!     *                        t0_aniso8,t1_aniso8,t0_iso8,t1_iso8)

      use constants, only: zero,six,nine,fifteen,pi
      use mx_limits, only: MXIFRQ
      use comm_FRGINF, only: NFRG
      use DYNPOL_DIST, only: DYNDD_LMO_ROT,DYNDQ_LMO_ROT,
     *                       DYNDO_LMO_ROT,DYNQQ_LMO_ROT

      IMPLICIT NONE 

      LOGICAL do_aniso_disp

      integer :: IN,JP,INPPTS,JNPPTS,IM,JM,MXPPTS,Ia,Ib,Ic,
     &           Ii,Ij,Ik,Ifreq
      double precision :: T2(3,3),T3(3,3,3),T4(3,3,3,3),FACT(12)
      double precision :: DYNTR(MXIFRQ*MXPPTS,NFRG)
      double precision :: DYNQQ_LMO_AVE(MXIFRQ*MXPPTS,NFRG) 
      double precision :: TMP,TMP1,TMP2,TMP3,TMP4,tmpDO,tmpDQ,tmpQQ
      double precision :: F8,R8,EPAIR80,EPAIR8,EPAIR80_aniso,
     &                    EPAIR8_aniso,E8DO,E8DQ,E8QQ,C8

      EPAIR80       = ZERO
      EPAIR8        = ZERO
      EPAIR80_aniso = ZERO
      EPAIR8_aniso  = ZERO
      
      IF(do_aniso_disp) THEN

!$omp parallel do schedule(dynamic) 
!$omp& default(private) shared(DYNDD_LMO_ROT,DYNDQ_LMO_ROT,
!$omp&   DYNDO_LMO_ROT,DYNQQ_LMO_ROT,T2,T3,T4,FACT)
!$omp& firstprivate(IM,JM,IN,JP,INPPTS,JNPPTS)
!$omp& reduction(+:EPAIR80_aniso,E8DO,E8DQ,E8QQ)

       Do 203 Ifreq = 1, MXIFRQ
          Do Ia = 1,3
          Do Ib = 1,3
          Do Ic = 1,3
          Do Ii = 1,3
          Do Ij = 1,3
          Do Ik = 1,3
c TMP1 are terms containing dipole-octopole pol 
             TMP1 =DYNDD_LMO_ROT(Ia,Ic,(IN+(Ifreq-1)*INPPTS),IM)*
     &             DYNDO_LMO_ROT(Ib,Ii,Ij,Ik,(JP+(Ifreq-1)*JNPPTS),JM)
     &            +DYNDO_LMO_ROT(Ia,Ic,Ii,Ij,(IN+(Ifreq-1)*INPPTS),IM)*
     &             DYNDD_LMO_ROT(Ib,Ik,(JP+(Ifreq-1)*JNPPTS),JM)
c TMP2 and TMP3 are terms containing dipole-quadrupole pol
             TMP2 =DYNDQ_LMO_ROT(Ia,Ic,Ii,(IN+(Ifreq-1)*INPPTS),IM)*
     &             DYNDQ_LMO_ROT(Ib,Ij,Ik,(JP+(Ifreq-1)*JNPPTS),JM)
c
             TMP3 =DYNDQ_LMO_ROT(Ii,Ia,Ib,(IN+(Ifreq-1)*INPPTS),IM)*
     &             DYNDQ_LMO_ROT(Ic,Ij,Ik,(JP+(Ifreq-1)*JNPPTS),JM)
c TMP4 are terms containing quadrupole-quadrupole pol
             TMP4 =DYNDD_LMO_ROT(Ia,Ii,(IN+(Ifreq-1)*INPPTS),IM)*
     &             DYNQQ_LMO_ROT(Ib,Ic,Ij,Ik,(JP+(Ifreq-1)*JNPPTS),JM)
     &            +DYNQQ_LMO_ROT(Ia,Ib,Ii,Ij,(IN+(Ifreq-1)*INPPTS),IM)*
     &             DYNDD_LMO_ROT(Ic,Ik,(JP+(Ifreq-1)*JNPPTS),JM)

c-----PX debug
             EPAIR80_aniso = EPAIR80_aniso + (FACT(Ifreq)/pi)*
     &            (T2(Ia,Ib)*T4(Ic,Ii,Ij,Ik)
     &            *(-TMP1/fifteen+TMP2/nine)
     &            +T3(Ia,Ib,Ic)*T3(Ii,Ij,Ik)
     &            *(TMP3/nine-TMP4/six))
c------PX debug: different contributions in Disp8
       E8DO=E8DO-FACT(Ifreq)*T2(Ia,Ib)*T4(Ic,Ii,Ij,Ik)*
     &             TMP1/(fifteen*pi)
       E8DQ=E8DQ+FACT(Ifreq)*(T2(Ia,Ib)*T4(Ic,Ii,Ij,Ik)*TMP2/(nine*pi)
     &                   +T3(Ia,Ib,Ic)*T3(Ii,Ij,Ik)*TMP3/(nine*pi))
       E8QQ=E8QQ-FACT(Ifreq)*T3(Ia,Ib,Ic)*T3(Ii,Ij,Ik)*TMP4/(six*pi)
c------END PX
c             EPAIR80_aniso = EPAIR80_aniso + TMP/pi
          Enddo
          Enddo
          Enddo
          Enddo
          Enddo
          Enddo
  203  CONTINUE
!$omp end parallel do

       EPAIR8_aniso  = F8 * EPAIR80_aniso
       END IF


c         -- Distributed spherical E8 (Local)
       C8 = zero
!!!$omp parallel do schedule(dynamic) default(private)
!!!$omp& shared(FACT,DYNTR,DYNQQ_LMO_AVE)
!!!$omp& firstprivate(IM,JM,IN,JP,INPPTS,JNPPTS)
!!!$omp& reduction(+: C8) 
        Do 204 Ifreq = 1,MXIFRQ
        C8 = C8 + (fifteen/pi)*FACT(Ifreq)*
     *           (DYNTR(IN+(Ifreq-1)*INPPTS,IM)
     *           *DYNQQ_LMO_AVE(JP+(Ifreq-1)*JNPPTS,JM)+
     *            DYNTR(JP+(Ifreq-1)*JNPPTS,JM)
     *           *DYNQQ_LMO_AVE(IN+(Ifreq-1)*INPPTS,IM))

  204   Continue
!!!$omp end parallel do

       EPAIR80   = - C8/R8
       EPAIR8    = F8*EPAIR80
       
    
       RETURN
       END


C*MODULE EFDRVR  *DECK Damping_for_Dispersion
      SUBROUTINE Damping_for_Dispersion(F6,F7,F8,DF6,DF7A,DF7B,
     *                                  DF8A,DF8B,IM,JM,IN,JP,
     *                                  IDISPS,E7DISP,E8DISP,IFGRAD,R,
     *                                  SAB,DSAB,SABSAVE,SABGRAD,MXPPTS)

      use constants
      use comm_EFPSCR, only: B
      use comm_FRGINF, only: NFRG
      IMPLICIT NONE 

      LOGICAL E7DISP,E8DISP,IFGRAD
      integer :: IM,JM,IN,JP,IDISPS,MXPPTS,ICRD
      double precision :: DSAB(6),SABSAVE(NFRG,NFRG,MXPPTS,MXPPTS),
     *                    SABGRAD(NFRG,NFRG,MXPPTS,MXPPTS,6)
      double precision :: F6,F7,F8,DF6,DF7A,DF7B,DF8A,DF8B
      double precision :: SAB,EX1,R,RB,RB2,RB3,RB4,RB5,RB6,RB7,RB8

!  IDISP = 1   overlap damping
!  IDISP = 0   Tang-Toennies damping 

!  F6 : damping function for Disp6
!  F7 : damping function for Disp7
!  F8 : damping function for Disp8

!  Initialization
      F6 = ONE
      F7 = ONE
      F8 = ONE
      DF6= ZERO
      DF7A=ZERO
      DF7B=ZERO
      DF8A=ZERO
      DF8B=ZERO
      SAB = ZERO
      call VCLR(DSAB,1,6)
! overlap damping
      IF (IDISPS .EQ. 1) THEN
         SAB = SABSAVE(IM,JM,IN,JP)
         IF (ABS(SAB) .GT. 1.0D-05) THEN
             RB = -2.0D+00*LOG(ABS(SAB))
             EX1 = SAB*SAB
             RB2 = RB*RB
             RB3 = RB2*RB
             F6 = ONE-(EX1*(ONE + SQRT(RB) + (RB/TWO) +
     *     (((RB)**(THREE/TWO))/SIX) + (RB2/24.0D+00) +
     *     ((RB**(FIVE/TWO))/120.0D+00) +(RB3/720.0D+00)))
c------PX
c        write(6,*) 'overlap F6=',F6
c------END PX
          IF(E7DISP) THEN !if disp7
!             F7 = F6 - (((RB)**(SEVEN/TWO))/5040.0D+00)
             F7 = ONE-(EX1*(ONE + SQRT(RB) + (RB/TWO)
     *           + (((RB)**(THREE/TWO))/SIX)
     *           + (RB2/24.0D+00)
     *           + ((RB**(FIVE/TWO))/120.0D+00)
     *           + (RB3/720.0D+00)
     *           + (((RB)**(SEVEN/TWO))/5040.0D+00)))
c-----PX
c        write(6,*) 'overlap F7=',F7
c-----END PX
          ELSE
             F7=ZERO
          ENDIF !E7DISP

          IF(E8DISP) THEN !if disp8 
             F8 = F6 - EX1*((RB**(SEVEN/TWO))/5040.0D+00)
     *               - EX1*((RB**FOUR)/40320.0D+00)               
c-----PX        
c        write(6,*) 'overlap F8=',F8
c-----END PX
          ELSE                                                
             F8 = ZERO                                          
          END IF !e8disp                                      

         ELSE ! no significant overlap
            F6 = ONE
            RB = ZERO
            SAB = ZERO
            DO ICRD = 1,6
               DSAB(ICRD) = ZERO
            ENDDO
c
            IF (E7DISP) THEN
             F7=ONE
            ELSE
             F7=ZERO
            ENDIF
            IF(E8DISP) THEN
             F8 = ONE
            ELSE
             F8 = ZERO
            ENDIF
         END IF !SAB

      ELSE !Tang-Toennies damping
         RB =  B*R
         EX1 = EXP(-RB)
         RB2 = RB*RB
         RB3 = RB*RB2
         RB4 = RB2*RB2
         RB5 = RB*RB4
         RB6 = RB*RB5
         RB7 = RB5*RB2
         RB8 = RB7*RB

         F6 = ONE-(EX1*(ONE + RB + RB2/TWO +
     *        RB3/SIX + RB4/24.0D+00 +
     *        RB5/120.0D+00 +RB6/720.0D+00))
c----PX
c        write(6,*) 'TT F6=',F6
c-----END PX
         IF (E7DISP) THEN
            F7 = F6 - EX1*RB7/5040.0D+00
c----PX
c        write(6,*) 'TT F7=',F7
c----END PX
         ELSE
            F7=ZERO
         ENDIF !E7DISP

         IF(E8DISP) THEN
           F8 = F6 - EX1*RB7/5040.0D+00 - EX1*RB8/40320.0D+00
c-----PX
c        write(6,*) 'TT F8=',F8
c-----END PX
         ELSE
           F8 = ZERO
         END IF !E8DISP 

      END IF ! IDISP damping

!----PX debug: turn off damping 
!        F6=ONE
!        F7=ONE
!        F8=ONE
!-----END PX

! derivative of damping functions 
      IF(IFGRAD) THEN
        IF (IDISPS .EQ. 1) THEN
           IF (ABS(SAB) .GT. 1.0D-05) THEN
               DO ICRD = 1,6
                   DSAB(ICRD) = SABGRAD(IM,JM,IN,JP,ICRD)
               END DO
               RB = -2.0D+00*LOG(ABS(SAB))
               EX1 = SAB*SAB
               RB2 = RB*RB
               RB3 = RB2*RB
               DF6 =-SAB*(ONE+((THREE/TWO)*SQRT(RB))
     *             -(ONE/SQRT(RB))+((FIVE/SIX)*RB)
     *             +((7.0D+00/24.0D+00)*(RB**(THREE/TWO)))
     *             +((9.0D+00/120.0D+00)*RB2)
     *             +((ONE/60.0D+00)*(RB**(FIVE/TWO)))
     *             +((ONE/360.0D+00)*RB3))
               IF(E7DISP) THEN
                  DF7A=-SAB*(ONE+((THREE/TWO)*SQRT(RB))
     *                 -(ONE/SQRT(RB))+((FIVE/SIX)*RB)
     *                 +((7.0D+00/24.0D+00)*(RB**(THREE/TWO)))
     *                 +((9.0D+00/120.0D+00)*RB2)
     *                 +((11.0D+00/720.0D+00)*(RB**(FIVE/TWO)))
     *                 +((ONE/360.0D+00)*RB3)
     *                 +((ONE/2520.0D+00)*(RB**(7.0D+00/TWO))))
                  DF7B=ZERO
               ELSE
                  DF7A=ZERO
                  DF7B=ZERO
               ENDIF !E7DISP

               IF(E8DISP) THEN
                  DF8A =-SAB*(ONE+((THREE/TWO)*SQRT(RB))
     *                 -(ONE/SQRT(RB))+((FIVE/SIX)*RB)
     *                 +((7.0D+00/24.0D+00)*(RB**(THREE/TWO)))
     *                 +((9.0D+00/120.0D+00)*RB2)
     *                 +((11.0D+00/720.0D+00)*(RB**(FIVE/TWO)))
     *                 +((ONE/360.0D+00)*RB3)
     *                 +((ONE/2520.0D+00)*(RB**(7.0D+00/TWO)))) 
     *                 +SAB*((ONE/5040.0D+00)*(RB**THREE)      
     *                     -(ONE/20160.0D+00)*(RB**FOUR))
                  DF8B=ZERO                                          
               ELSE                                                
                  DF8A=ZERO                                          
                  DF8B=ZERO                                          
               END IF !e8disp                                      
c----PX debug
c        write(6,*) 'DF8A=',DF8A
C-----END PX
           ELSE ! no significant overlap
               RB = ZERO
               SAB = ZERO
               DO ICRD = 1,6
                  DSAB(ICRD) = ZERO
               ENDDO
C------PX
               DF6=ZERO
C--------END PX
c
               IF (E7DISP) THEN
                DF7A=ZERO
                DF7B=ZERO
               ENDIF
               IF(E8DISP) THEN
                DF8A = ZERO
                DF8B = ZERO
               ENDIF
             END IF !SAB
         
        ELSE ! Tang-Toennies damping
             RB =  B*R
             EX1 = EXP(-RB)
             RB2 = RB*RB
             RB3 = RB*RB2
             RB4 = RB2*RB2
             RB5 = RB*RB4
             RB6 = RB*RB5
             RB7 = RB5*RB2

!------PX replace with an equivalent but much simpler formula
!             DF6=-B*EX1*(-ONE+(RB*((ONE/RB)-ONE))+
!     *           ((RB2/TWO)*((TWO/RB)-ONE))+
!     *           ((RB3/SIX)*((THREE/RB)-ONE))+
!     *           ((RB4/24.0D+00)*((FOUR/RB)-ONE))+
!     *           ((RB5/120.D+00)*((FIVE/RB)-ONE))+
!     *           ((RB6/720.0D+00)*((SIX/RB)-ONE)))

             DF6=B*EX1*RB6/720.0D+00
!------END PX

             IF (E7DISP) THEN
!                DF7B=-B*EX1*(-ONE+(RB*((ONE/RB)-ONE))+
!     *                ((RB2/TWO)*((TWO/RB)-ONE))+
!     *                ((RB3/SIX)*((THREE/RB)-ONE))+
!     *                ((RB4/24.0D+00)*((FOUR/RB)-ONE))+
!     *                ((RB5/120.D+00)*((FIVE/RB)-ONE))+
!     *                ((RB6/720.0D+00)*((SIX/RB)-ONE))+
!     *                ((RB7/5040.0D+00)*((SEVEN/RB)-ONE)))
                DF7B=B*EX1*RB7/5040.0D+00
                DF7A=ZERO
              ELSE
                DF7A=ZERO
                DF7B=ZERO
             ENDIF !E7DISP

             IF(E8DISP) THEN
               DF8A=ZERO
!-------PX
!               DF8B=DF7B
!     *              -B*EX1*((RB8/40320.0D+00)*((EIGHT/RB)-ONE))
               DF8B=B*EX1*RB8/40320.0D+00
!-------END PX
             ELSE
               DF8A=ZERO
               DF8B=ZERO
             END IF !E8DISP 
        ENDIF ! IDISP
      ENDIF ! IFGRAD
!----PX debug
!        DF6 =ZERO
!        DF7A=ZERO
!        DF7B=ZERO
!        DF8A=ZERO
!        DF8B=ZERO
!----END PX
      RETURN
      END


C*MODULE EFDRVR  *DECK GRAD_EFP_DISP
      SUBROUTINE GRAD_EFP_DISP(C6,F6,DF6,F7,DF7A,DF7B,C8,F8,DF8A,DF8B,
     *                        X,Y,Z,R,R6,R7,R8,R9,R10,
     *                        DSAB,conx,cony,conz,SWF,SWFDX,SWFDY,SWFDZ,
     *                        COMX,COMY,COMZ,PBCX,PBCY,PBCZ,IP1,IM,JM,
     *                        E7DISP,E8DISP,IDISPS,EPAIR,DER_EPAIR,
     *                        EPAIR70,EPAIR7,DER_EPAIR70,DER_EPAIR70_IJ,
     *                        EPAIR8)

      use constants, only: zero,three,four,six,eight
      use mx_limits, only: MXFRG
      use comm_EFPPAR, only: EFP
      use comm_FGRAD, only: EFCENT
      use comm_EFPGES, only: DESFRG

      IMPLICIT NONE 

      LOGICAL E7DISP,E8DISP
      integer :: IP1,IM,JM,IDISPS
      double precision :: DSAB(6),DER_EPAIR(6),DER_EPAIR70(3),
     *                    DER_EPAIR70_IJ(3),DER_EPAIR7(3),
     *                    DER_EPAIR7_IJ(3),DER_EPAIR8(6),
     *                    FORDIS(3,MXFRG),FORDIS6(3,MXFRG),
     *                    FORDIS7(3,MXFRG),FORDIS8(3,MXFRG),
     *                    TORDIS(3,MXFRG),TORDIS6(3,MXFRG),
     *                    TORDIS7(3,MXFRG),TORDIS8(3,MXFRG)
      double precision :: C6,F6,DF6,F7,DF7A,DF7B,C8,F8,DF8A,DF8B,
     *                    X,Y,Z,R,R6,R7,R8,R9,R10
      double precision :: conx,cony,conz,SWF,SWFDX,SWFDY,SWFDZ,
     *                   COMX,COMY,COMZ,PBCX,PBCY,PBCZ
      double precision :: EPAIR,EPAIR70,EPAIR7,EPAIR8 
!-----PX
!      double precision :: C6_orig,DF6_orig,R6_orig,SIX_orig,R8_orig,
!     *                    F6_orig,EPAIR_orig

c---END PX
!---- initialize force and torque for dispersion
      call vclr(FORDIS,1,3*MXFRG)
      call vclr(FORDIS6,1,3*MXFRG)
      call vclr(FORDIS7,1,3*MXFRG)
      call vclr(FORDIS8,1,3*MXFRG)
      call vclr(TORDIS,1,3*MXFRG)
      call vclr(TORDIS6,1,3*MXFRG)
      call vclr(TORDIS7,1,3*MXFRG)
      call vclr(TORDIS8,1,3*MXFRG)
      call vclr(DER_EPAIR,1,6)
      IF(E7DISP) THEN
         CALL VCLR(DER_EPAIR7,1,3)
         CALL VCLR(DER_EPAIR7_IJ,1,3)
      ENDIF
      IF(E8DISP) CALL VCLR(DER_EPAIR8,1,6)


!$omp parallel  default(none)
!$omp& firstprivate(IP1,IM,JM)
!$omp& shared(IDISPS,E8DISP,E7DISP,
!$omp&        C6,F6,DF6,F7,DF7A,DF7B,C8,F8,DF8A,DF8B,
!$omp&        DSAB,R,R6,R7,R8,R9,R10,X,Y,Z,
!$omp&        EFP,EFCENT,EPAIR70,DER_EPAIR70,DER_EPAIR70_IJ,
!$omp&        SWF,SWFDX,SWFDY,SWFDZ,conx,cony,conz,
!$omp&        COMX,COMY,COMZ,PBCX,PBCY,PBCZ,EPAIR,EPAIR7,EPAIR8,
!$omp&        DER_EPAIR,DER_EPAIR7,DER_EPAIR8,DER_EPAIR7_IJ)

C Derivative of E6+E8 with overlap damping and TT formula
      IF (IDISPS.EQ.1) THEN ! overlap damping
C Translation
          IF(.not.E8DISP) THEN 
          DER_EPAIR(1)=(((FOUR/THREE)*C6*DF6*DSAB(1))/R6)-
     *                  ((EIGHT*C6*F6*X)/R8)
          DER_EPAIR(2)=(((FOUR/THREE)*C6*DF6*DSAB(2))/R6)-
     *                  ((EIGHT*C6*F6*Y)/R8)
          DER_EPAIR(3)=(((FOUR/THREE)*C6*DF6*DSAB(3))/R6)-
     *                  ((EIGHT*C6*F6*Z)/R8)
          ELSE
          DER_EPAIR(1)=((C6*DF6*DSAB(1))/R6)-
     *                  ((SIX*C6*F6*X)/R8)

          DER_EPAIR(2)=((C6*DF6*DSAB(2))/R6)-
     *                  ((SIX*C6*F6*Y)/R8)

          DER_EPAIR(3)=((C6*DF6*DSAB(3))/R6)-
     *                  ((SIX*C6*F6*Z)/R8)
          ENDIF
C Rotation
          IF(.not.E8DISP) THEN
          DER_EPAIR(4)=(((FOUR/THREE)*C6*DF6*DSAB(4))/R6)
     *                 -(EIGHT*C6*F6*((Y*(EFP(3,IP1)-
     *             EFCENT(3,IM)))-(Z*(EFP(2,IP1)-EFCENT(2,IM))))/R8)
          DER_EPAIR(5)=(((FOUR/THREE)*C6*DF6*DSAB(5))/R6)
     *                 -(EIGHT*C6*F6*((Z*(EFP(1,IP1)-
     *             EFCENT(1,IM)))-(X*(EFP(3,IP1)-EFCENT(3,IM))))/R8)
          DER_EPAIR(6)=(((FOUR/THREE)*C6*DF6*DSAB(6))/R6)
     *                 -(EIGHT*C6*F6*((X*(EFP(2,IP1)-
     *             EFCENT(2,IM)))-(Y*(EFP(1,IP1)-EFCENT(1,IM))))/R8)
          ELSE
          DER_EPAIR(4)=((C6*DF6*DSAB(4))/R6)
     *                 -(SIX*C6*F6*((Y*(EFP(3,IP1)-
     *             EFCENT(3,IM)))-(Z*(EFP(2,IP1)-EFCENT(2,IM))))/R8)

          DER_EPAIR(5)=((C6*DF6*DSAB(5))/R6)
     *                 -(SIX*C6*F6*((Z*(EFP(1,IP1)-
     *             EFCENT(1,IM)))-(X*(EFP(3,IP1)-EFCENT(3,IM))))/R8)

          DER_EPAIR(6)=((C6*DF6*DSAB(6))/R6)
     *                 -(SIX*C6*F6*((X*(EFP(2,IP1)-
     *             EFCENT(2,IM)))-(Y*(EFP(1,IP1)-EFCENT(1,IM))))/R8)
          ENDIF

c$$$ for EFMO
          conx = -((EIGHT*C6*F6*X)/R8)
          cony = -((EIGHT*C6*F6*Y)/R8)
          conz = -((EIGHT*C6*F6*Z)/R8)

      ELSE ! TT damping
C Translation
          IF(.not.E8DISP) THEN
          DER_EPAIR(1)=(X/R)*(FOUR/THREE)*C6*((DF6/R6)-
     *                 ((SIX*F6)/R7))
          DER_EPAIR(2)=(Y/R)*(FOUR/THREE)*C6*((DF6/R6)-
     *                 ((SIX*F6)/R7))
          DER_EPAIR(3)=(Z/R)*(FOUR/THREE)*C6*((DF6/R6)-
     *                 ((SIX*F6)/R7))
          ELSE
          DER_EPAIR(1)=(X/R)*C6*((DF6/R6)-((SIX*F6)/R7))
          DER_EPAIR(2)=(Y/R)*C6*((DF6/R6)-((SIX*F6)/R7))
          DER_EPAIR(3)=(Z/R)*C6*((DF6/R6)-((SIX*F6)/R7))
          ENDIF
C Rotation
          IF(.not.E8DISP) THEN
          DER_EPAIR(4)=(((Y*(EFP(3,IP1)-EFCENT(3,IM)))
     *                 -(Z*(EFP(2,IP1)-EFCENT(2,IM))))/R)*
     *                 (FOUR/THREE)*C6*((DF6/R6)-(SIX*F6/R7))
          DER_EPAIR(5)=(((Z*(EFP(1,IP1)-EFCENT(1,IM)))
     *                 -(X*(EFP(3,IP1)-EFCENT(3,IM))))/R)*
     *               (FOUR/THREE)*C6*((DF6/R6)-(SIX*F6/R7))
          DER_EPAIR(6)=(((X*(EFP(2,IP1)-EFCENT(2,IM)))
     *                 -(Y*(EFP(1,IP1)-EFCENT(1,IM))))/R)*
     *               (FOUR/THREE)*C6*((DF6/R6)-(SIX*F6/R7))
          ELSE
          DER_EPAIR(4)=(((Y*(EFP(3,IP1)-EFCENT(3,IM)))
     *                 -(Z*(EFP(2,IP1)-EFCENT(2,IM))))/R)*
     *                 C6*((DF6/R6)-(SIX*F6/R7))
          DER_EPAIR(5)=(((Z*(EFP(1,IP1)-EFCENT(1,IM)))
     *                 -(X*(EFP(3,IP1)-EFCENT(3,IM))))/R)*
     *                 C6*((DF6/R6)-(SIX*F6/R7))
          DER_EPAIR(6)=(((X*(EFP(2,IP1)-EFCENT(2,IM)))
     *                 -(Y*(EFP(1,IP1)-EFCENT(1,IM))))/R)*
     *                 C6*((DF6/R6)-(SIX*F6/R7))
          ENDIF
c$$$ for EFMO--if not overlap, turn off
          conx = 0
          cony = 0
          conz = 0
      ENDIF !IDISPS
 
c Disp7
      IF(E7DISP) THEN
!         CALL VCLR(DER_EPAIR7,1,3)
!         CALL VCLR(DER_EPAIR7_IJ,1,3)

C Translation
         DER_EPAIR7(1)=-(((DF7A*DSAB(1)+(DF7B*(X/R)))*EPAIR70)
     *                  +(F7*DER_EPAIR70(1)))
         DER_EPAIR7(2)=-(((DF7A*DSAB(2)+(DF7B*(Y/R)))*EPAIR70)
     *                  +(F7*DER_EPAIR70(2)))
         DER_EPAIR7(3)=-(((DF7A*DSAB(3)+(DF7B*(Z/R)))*EPAIR70)
     *                  +(F7*DER_EPAIR70(3)))
C Rotation 
         DER_EPAIR7_IJ(1)=(((DF7A*DSAB(4)+(DF7B*
     * ((Y*(EFP(3,IP1)-EFCENT(3,IM))-Z*((EFP(2,IP1)-EFCENT(2,IM))))/R)))
     *                 *(EPAIR70))+(F7*DER_EPAIR70_IJ(1)))
         DER_EPAIR7_IJ(2)=(((DF7A*DSAB(5)+(DF7B*
     * ((Z*(EFP(1,IP1)-EFCENT(1,IM))-X*((EFP(3,IP1)-EFCENT(3,IM))))/R)))
     *                 *(EPAIR70))+(F7*DER_EPAIR70_IJ(2)))
         DER_EPAIR7_IJ(3)=(((DF7A*DSAB(6)+(DF7B*
     * ((X*(EFP(2,IP1)-EFCENT(2,IM))-Y*((EFP(1,IP1)-EFCENT(1,IM))))/R)))
     *                 *(EPAIR70))+(F7*DER_EPAIR70_IJ(3)))
      ENDIF !E7DISP

c Disp8 
      IF(E8DISP) THEN
!         CALL VCLR(DER_EPAIR8,1,6)
! Translation
         DER_EPAIR8(1)=((C8*DF8A*DSAB(1))/R8+((X/R)*DF8B*C8)/R8)-
     *                  ((EIGHT*C8*F8*X)/R10)

         DER_EPAIR8(2)=((C8*DF8A*DSAB(2))/R8+((Y/R)*DF8B*C8)/R8)-
     *                  ((EIGHT*C8*F8*Y)/R10)

         DER_EPAIR8(3)=((C8*DF8A*DSAB(3))/R8+((Z/R)*DF8B*C8)/R8)-
     *                  ((EIGHT*C8*F8*Z)/R10)
! Rotation
         IF(IDISPS.EQ.1) THEN ! overlap damping
         DER_EPAIR8(4)=(C8*DF8A*DSAB(4)/R8)
     *                -(EIGHT*C8*F8/R10)*((Y*(EFP(3,IP1)-
     *            EFCENT(3,IM)))-(Z*(EFP(2,IP1)-EFCENT(2,IM))))

         DER_EPAIR8(5)=(C8*DF8A*DSAB(5)/R8)
     *                -(EIGHT*C8*F8/R10)*((Z*(EFP(1,IP1)-
     *            EFCENT(1,IM)))-(X*(EFP(3,IP1)-EFCENT(3,IM))))

         DER_EPAIR8(6)=(C8*DF8A*DSAB(6)/R8)
     *                -EIGHT*C8*F8/R10*((X*(EFP(2,IP1)-
     *            EFCENT(2,IM)))-(Y*(EFP(1,IP1)-EFCENT(1,IM))))

         ELSE ! TT damping
         DER_EPAIR8(4)=(((Y*(EFP(3,IP1)-EFCENT(3,IM)))
     *                -(Z*(EFP(2,IP1)-EFCENT(2,IM))))/R)*
     *                C8*((DF8B/R8)-(EIGHT*F8/R9))
         DER_EPAIR8(5)=(((Z*(EFP(1,IP1)-EFCENT(1,IM)))
     *                -(X*(EFP(3,IP1)-EFCENT(3,IM))))/R)*
     *                C8*((DF8B/R8)-(EIGHT*F8/R9))
         DER_EPAIR8(6)=(((X*(EFP(2,IP1)-EFCENT(2,IM)))
     *                -(Y*(EFP(1,IP1)-EFCENT(1,IM))))/R)*
     *                C8*((DF8B/R8)-(EIGHT*F8/R9))
         ENDIF
      ENDIF
!$omp end parallel

C E. Guidez: FOrces: Add the E7*(dSWD/dx),E7*(dSWD/dy),E7*(dSWD/dz)
C terms to the
C gradients for 1 fragment and -E7*(dSWD/dx),-E7*(dSWD/dy),-E7*(dSWD/dz)
C to the other fragment since FA(B)=-FB(A)

C FORCES
      FORDIS6(1,JM)=-((EPAIR)*SWFDX)+((DER_EPAIR(1))*SWF)
      FORDIS6(2,JM)=-((EPAIR)*SWFDY)+((DER_EPAIR(2))*SWF)
      FORDIS6(3,JM)=-((EPAIR)*SWFDZ)+((DER_EPAIR(3))*SWF)
      FORDIS6(1,IM)=-FORDIS6(1,JM)
      FORDIS6(2,IM)=-FORDIS6(2,JM)
      FORDIS6(3,IM)=-FORDIS6(3,JM)

      IF(E7DISP) THEN
         FORDIS7(1,JM)=-((EPAIR7)*SWFDX)+((DER_EPAIR7(1))*SWF)
         FORDIS7(2,JM)=-((EPAIR7)*SWFDY)+((DER_EPAIR7(2))*SWF)
         FORDIS7(3,JM)=-((EPAIR7)*SWFDZ)+((DER_EPAIR7(3))*SWF)
         FORDIS7(1,IM)=-FORDIS7(1,JM)
         FORDIS7(2,IM)=-FORDIS7(2,JM)
         FORDIS7(3,IM)=-FORDIS7(3,JM)
      ELSE
         FORDIS7(1,JM)=ZERO
         FORDIS7(2,JM)=ZERO
         FORDIS7(3,JM)=ZERO
         FORDIS7(1,IM)=ZERO
         FORDIS7(2,IM)=ZERO
         FORDIS7(3,IM)=ZERO
      ENDIF
      IF(E8DISP) THEN
         FORDIS8(1,JM)=-((EPAIR8)*SWFDX)+(DER_EPAIR8(1)*SWF)
         FORDIS8(2,JM)=-((EPAIR8)*SWFDY)+(DER_EPAIR8(2)*SWF)
         FORDIS8(3,JM)=-((EPAIR8)*SWFDZ)+(DER_EPAIR8(3)*SWF)
         FORDIS8(1,IM)=-FORDIS8(1,JM)
         FORDIS8(2,IM)=-FORDIS8(2,JM)
         FORDIS8(3,IM)=-FORDIS8(3,JM)
      ELSE
         FORDIS8(1,JM)=ZERO
         FORDIS8(2,JM)=ZERO
         FORDIS8(3,JM)=ZERO
         FORDIS8(1,IM)=ZERO
         FORDIS8(2,IM)=ZERO
         FORDIS8(3,IM)=ZERO
      ENDIF

      FORDIS(1,IM)=FORDIS6(1,IM)+FORDIS7(1,IM)+FORDIS8(1,IM)
      FORDIS(2,IM)=FORDIS6(2,IM)+FORDIS7(2,IM)+FORDIS8(2,IM)
      FORDIS(3,IM)=FORDIS6(3,IM)+FORDIS7(3,IM)+FORDIS8(3,IM)
      FORDIS(1,JM)=FORDIS6(1,JM)+FORDIS7(1,JM)+FORDIS8(1,JM)
      FORDIS(2,JM)=FORDIS6(2,JM)+FORDIS7(2,JM)+FORDIS8(2,JM)
      FORDIS(3,JM)=FORDIS6(3,JM)+FORDIS7(3,JM)+FORDIS8(3,JM)

      DESFRG(1,JM) = DESFRG(1,JM) +FORDIS(1,JM)
      DESFRG(2,JM) = DESFRG(2,JM) +FORDIS(2,JM)
      DESFRG(3,JM) = DESFRG(3,JM) +FORDIS(3,JM)
      DESFRG(1,IM) = DESFRG(1,IM) +FORDIS(1,IM)
      DESFRG(2,IM) = DESFRG(2,IM) +FORDIS(2,IM)
      DESFRG(3,IM) = DESFRG(3,IM) +FORDIS(3,IM)

C TORQUES
      TORDIS6(1,IM) =SWF*(DER_EPAIR(4))
      TORDIS6(2,IM) =SWF*(DER_EPAIR(5))
      TORDIS6(3,IM) =SWF*(DER_EPAIR(6))
      TORDIS6(1,JM) = -TORDIS6(1,IM)-((FORDIS6(2,IM)*(-COMZ+PBCZ))-
     *                (FORDIS6(3,IM)*(-COMY+PBCY)))
      TORDIS6(2,JM) = -TORDIS6(2,IM)-((FORDIS6(3,IM)*(-COMX+PBCX))-
     *                (FORDIS6(1,IM)*(-COMZ+PBCZ)))
      TORDIS6(3,JM) = -TORDIS6(3,IM)-((FORDIS6(1,IM)*(-COMY+PBCY))-
     *                (FORDIS6(2,IM)*(-COMX+PBCX)))

      IF (E7disp) then
          TORDIS7(1,IM)=-SWF*(DER_EPAIR7_IJ(1))
          TORDIS7(2,IM)=-SWF*(DER_EPAIR7_IJ(2))
          TORDIS7(3,IM)=-SWF*(DER_EPAIR7_IJ(3))
          TORDIS7(1,JM) = -TORDIS7(1,IM)
     *                    -((FORDIS7(2,IM)*(-COMZ+PBCZ))-
     *                     (FORDIS7(3,IM)*(-COMY+PBCY)))
          TORDIS7(2,JM) = -TORDIS7(2,IM)
     *                    -((FORDIS7(3,IM)*(-COMX+PBCX))-
     *                     (FORDIS7(1,IM)*(-COMZ+PBCZ)))
          TORDIS7(3,JM) = -TORDIS7(3,IM)
     *                    -((FORDIS7(1,IM)*(-COMY+PBCY))-
     *                     (FORDIS7(2,IM)*(-COMX+PBCX)))

      ELSE
          TORDIS7(1,IM)=ZERO
          TORDIS7(2,IM)=ZERO
          TORDIS7(3,IM)=ZERO
          TORDIS7(1,JM)=ZERO
          TORDIS7(2,JM)=ZERO
          TORDIS7(3,JM)=ZERO
      ENDIF
      IF(E8DISP) THEN
         TORDIS8(1,IM)=SWF*DER_EPAIR8(4)
         TORDIS8(2,IM)=SWF*DER_EPAIR8(5)
         TORDIS8(3,IM)=SWF*DER_EPAIR8(6)
         TORDIS8(1,JM)=-TORDIS8(1,IM)-((FORDIS8(2,IM)*(-COMZ+PBCZ))-
     *                  (FORDIS8(3,IM)*(-COMY+PBCY)))
         TORDIS8(2,JM)=-TORDIS8(2,IM)-((FORDIS8(3,IM)*(-COMX+PBCX))-
     *                  (FORDIS8(1,IM)*(-COMZ+PBCZ)))
         TORDIS8(3,JM)=-TORDIS8(3,IM)-((FORDIS8(1,IM)*(-COMY+PBCY))-
     *                  (FORDIS8(2,IM)*(-COMX+PBCX)))
      ELSE
         TORDIS8(1,IM)=ZERO
         TORDIS8(2,IM)=ZERO
         TORDIS8(3,IM)=ZERO
         TORDIS8(1,JM)=ZERO
         TORDIS8(2,JM)=ZERO
         TORDIS8(3,JM)=ZERO
      END IF

      TORDIS(1,IM)=TORDIS6(1,IM)+TORDIS7(1,IM)+TORDIS8(1,IM)
      TORDIS(2,IM)=TORDIS6(2,IM)+TORDIS7(2,IM)+TORDIS8(2,IM)
      TORDIS(3,IM)=TORDIS6(3,IM)+TORDIS7(3,IM)+TORDIS8(3,IM)
      TORDIS(1,JM)=TORDIS6(1,JM)+TORDIS7(1,JM)+TORDIS8(1,JM)
      TORDIS(2,JM)=TORDIS6(2,JM)+TORDIS7(2,JM)+TORDIS8(2,JM)
      TORDIS(3,JM)=TORDIS6(3,JM)+TORDIS7(3,JM)+TORDIS8(3,JM)

      DESFRG(4,JM)=DESFRG(4,JM)+TORDIS(1,JM)
      DESFRG(5,JM)=DESFRG(5,JM)+TORDIS(2,JM)
      DESFRG(6,JM)=DESFRG(6,JM)+TORDIS(3,JM)
      DESFRG(4,IM)=DESFRG(4,IM)+TORDIS(1,IM)
      DESFRG(5,IM)=DESFRG(5,IM)+TORDIS(2,IM)
      DESFRG(6,IM)=DESFRG(6,IM)+TORDIS(3,IM)

 100  FORMAT(A15,F12.8)

      RETURN
      END

C*MODULE EFDRVR  *DECK DISP_EFP_FP 
      SUBROUTINE DISP_EFP_FP(C6,F6,DF6,X,Y,Z,R,R6,R7,R8,
     *                      DSAB,SWF,SWFDX,SWFDY,SWFDZ,
     *                      IFP1,IFP1_COM,IFP2,IFP2_COM,NEQUAL1,NEQUAL2,
     *                      COMX,COMY,COMZ,PBCX,PBCY,PBCZ,IP1,
     *                      IDISPS,EPAIR,DER_EPAIR)
      
      use constants, only: one,three,four,six,eight
      use comm_EFPFRC
      use comm_EFPPAR, only: EFP
      IMPLICIT NONE 

      integer :: IFP1_COM,IFP2_COM,NEQUAL1,NEQUAL2,IP1,IDISPS,I,IFP
      integer :: IFP1(10),IFP2(10)
      double precision :: DSAB(6),DER_EPAIR(6),DER_EPAIR_FP_6(3)
      double precision :: C6,F6,DF6,X,Y,Z,R,R6,R7,R8,
     *                    SWF,SWFDX,SWFDY,SWFDZ,COMX,COMY,COMZ,
     *                    PBCX,PBCY,PBCZ,EPAIR,COEF

           FRCTRQ(1,IFP2_COM) = FRCTRQ(1,IFP2_COM)-((EPAIR)*SWFDX)
           FRCTRQ(2,IFP2_COM) = FRCTRQ(2,IFP2_COM)-((EPAIR)*SWFDY)
           FRCTRQ(3,IFP2_COM) = FRCTRQ(3,IFP2_COM)-((EPAIR)*SWFDZ)
           FRCTRQ(1,IFP1_COM) = FRCTRQ(1,IFP1_COM)+((EPAIR)*SWFDX)
           FRCTRQ(2,IFP1_COM) = FRCTRQ(2,IFP1_COM)+((EPAIR)*SWFDY)
           FRCTRQ(3,IFP1_COM) = FRCTRQ(3,IFP1_COM)+((EPAIR)*SWFDZ)

           DO I=1,NEQUAL1
              IFP = IFP1(I)
              COEF=ONE/NEQUAL1
              FRCTRQ(1,IFP) = FRCTRQ(1,IFP)+
     *                        SWF*DER_EPAIR(1)*COEF
              FRCTRQ(2,IFP) = FRCTRQ(2,IFP)+
     *                        SWF*DER_EPAIR(2)*COEF
              FRCTRQ(3,IFP) = FRCTRQ(3,IFP)+
     *                        SWF*DER_EPAIR(3)*COEF

              IF (IDISPS.EQ.1) THEN
                  DER_EPAIR_FP_6(1)=(((FOUR/THREE)*C6*DF6*DSAB(4))/R6)
     *                              -(EIGHT*C6*F6*((Y*(EFP(3,IP1)-
     *               FRCCRD(3,IFP)))-(Z*(EFP(2,IP1)-FRCCRD(2,IFP))))/R8)

                  DER_EPAIR_FP_6(2)=(((FOUR/THREE)*C6*DF6*DSAB(5))/R6)
     *                              -(EIGHT*C6*F6*((Z*(EFP(1,IP1)-
     *               FRCCRD(1,IFP)))-(X*(EFP(3,IP1)-FRCCRD(3,IFP))))/R8)

                  DER_EPAIR_FP_6(3)=(((FOUR/THREE)*C6*DF6*DSAB(6))/R6)
     *                              -(EIGHT*C6*F6*((X*(EFP(2,IP1)-
     *               FRCCRD(2,IFP)))-(Y*(EFP(1,IP1)-FRCCRD(1,IFP))))/R8)

               ELSE
                  DER_EPAIR_FP_6(1)=(((Y*(EFP(3,IP1)-
     *              FRCCRD(3,IFP)))-(Z*(EFP(2,IP1)-FRCCRD(2,IFP))))/R)*
     *               (FOUR/THREE)*C6*((DF6/R6)-(SIX*F6/R7))

                  DER_EPAIR_FP_6(2)=(((Z*(EFP(1,IP1)-
     *              FRCCRD(1,IFP)))-(X*(EFP(3,IP1)-FRCCRD(3,IFP))))/R)*
     *               (FOUR/THREE)*C6*((DF6/R6)-(SIX*F6/R7))

                  DER_EPAIR_FP_6(3)=(((X*(EFP(2,IP1)-
     *              FRCCRD(2,IFP)))-(Y*(EFP(1,IP1)-FRCCRD(1,IFP))))/R)*
     *               (FOUR/THREE)*C6*((DF6/R6)-(SIX*F6/R7))

              ENDIF

              FRCTRQ(4,IFP) = FRCTRQ(4,IFP)+SWF*COEF*(DER_EPAIR_FP_6(1))
              FRCTRQ(5,IFP) = FRCTRQ(5,IFP)+SWF*COEF*(DER_EPAIR_FP_6(2))
              FRCTRQ(6,IFP) = FRCTRQ(6,IFP)+SWF*COEF*(DER_EPAIR_FP_6(3))

          END DO

          DO I=1,NEQUAL2
             IFP = IFP2(I)
             COEF=ONE/NEQUAL2
             FRCTRQ(1,IFP) = FRCTRQ(1,IFP)+SWF*(DER_EPAIR(1)*COEF)
             FRCTRQ(2,IFP) = FRCTRQ(2,IFP)+SWF*(DER_EPAIR(2)*COEF)
             FRCTRQ(3,IFP) = FRCTRQ(3,IFP)+SWF*(DER_EPAIR(3)*COEF)

             FRCTRQ(4,IFP)= -FRCTRQ(4,IFP1(I))-(((DER_EPAIR_FP_6(2))
     *                 *(-COMZ+PBCZ))-(DER_EPAIR_FP_6(3))*(-COMY+PBCY))
             FRCTRQ(5,IFP)= -FRCTRQ(5,IFP1(I))-(((DER_EPAIR_FP_6(3))
     *                 *(-COMX+PBCX))-(DER_EPAIR_FP_6(1))*(-COMZ+PBCZ))
             FRCTRQ(6,IFP)= -FRCTRQ(6,IFP1(I))-(((DER_EPAIR_FP_6(1))
     *                 *(-COMY+PBCY))-(DER_EPAIR_FP_6(2))*(-COMX+PBCX))

          END DO

      RETURN
      END
