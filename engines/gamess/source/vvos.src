C  1 Aug 14 - MWS - patch over the rare case of negative Dj pre-norm.
C 27 Jan 14 - MWS - demonstrate direct VVO generation by SVD rotation.
C 20 Sep 13 - ACW,MWS - make VVOs for high spin ROHF and all GVB cases.
C  3 Jul 13 - MWS - VVOS: orthog atomic side, TFAO5D: fix d transf.
C  6 Feb 13 - MWS - remove small AXP compiler warning
C 05 MAY 12 - DPT - SAFLG logical keyword added to DETPAR common block
C 17 APR 12 - MWS - DEFEND AGAINST CASE OF VVO'S WITH NO EXTERNAL ORBS
C  8 Mar 12 - EAH - ADD 1ST & 2ND ROW STATE AVERAGED TRANSITION METALS,
C                   FIX THE 4th & 5th ROW MAIN GROUP WTBS: H-XE now OK.
C  7 MAR 12 - MWS - ALIGN DETWFN COMMON
C 18 Nov 11 - FZ  - adjust overlap selection arguments
C 11 AUG 10 - DGF - SYNCH COMMON BLOCK PRPOPT
C 14 AUG 09 - MWS - VVOS: HANDLE RARE CASE OF ZERO OVERLAP WITH OCC.ORBS
C  1 MAY 09 - MWS - SWITCH TO FSOCI FOR THE CI-SD STEP, WHEN >2 E-
C 15 DEC 08 - DGF - SYNCHRONISE PRPOPT
C 20 NOV 08 - MWS - CORRECT PRINTING OF EXTERNAL RANGE
C 23 OCT 08 - MWS - INITIALIZE NRNFG CORRECTLY
C 18 JUL 08 - MWS - BUILD IN MORE ATOMIC BASIS SETS, COPE WITH 5D MBS
C  7 DEC 07 - MWS - DEVELOP AUTOMATED BONDING ANALYSIS PATHWAY
C 10 JUL 06 - MWS - FIX MBS COUNTING ESTIMATES
C  8 MAY 06 - MWS - NEW MODULE FOR QUAMBO'S/VIRTUAL VALENCE ORBITALS
C
C*MODULE VVOS    *DECK VVOS
C>
C> @brief      Routine drives formation of valence virtual orbitals.
C>
C> @author     Mike Schmidt, 2007
C>
C> @details   Routine drives formation of valence virtual orbitals.
C>            W.C.LU, C.Z.WANG, T.L.CHAN, K.RUEDENBERG, K.M.HO
C>            PHYS.REV.B 70, 041101/1-4(2004)
C>
C> @date November 4, 2012-Aaron West
C> -Modified nvvos_numcor integer calls.
C>
C> @date January 16, 2013-Aaron West
C> -Implemented MCSCF VVOS for use in SVD.
C>  Copy back over occupieds since then SVD.
C>
C> @date January 21, 2013-Aaron West
C> -Modified OVLSEL arguments.
C>  We need this modification to save the correct VVOS
C>  rather than some garbage from the virtual space.
C>
C> @date February 23, 2013-Aaron West
C> -Save the MBS AO SV matrix for possible use in SVD.
C>  We do VVOS now for ISVDOP=10.
C>  So, exit out early when IVVOS.ne.0.and.ISVDOP.eq.10.
C>      We need to leave the ISVDOP option open for both
C>      true VVOS and non-true VVOS runs.
C>      i.e. true VVOS create the extra ORMAS group
C>           and replace the canonical MCSCF orbitals.
C>      So, it is fair game to assume the above setting is
C>          unique to SVD-type runs with ISVDOP.eq.10.
C>
C> @date November 1, 2020 - George Schoendorff
C> - Added PADSVD to ORNTMO common block
C>
C> @date November 2, 2020 - George Schoendorff
C> - Added support for the relativistic AAMBS and for MCPs
C>
C> @date June 7, 2024 - George Schoendorff
C> - Move IVVTYP overrides to LMOINP to correctly set IVVTYP
C>
C> @param VEC  Canonical orbitals on entry and modified orbitals
C>             on exit
C> @param FAO  FOCK operator in the AO basis on entry
C> @param EIG  Eigenvalue matrix
C> @param WRK1 L3 sized scratch space
C> @param WRK2 L3 sized scratch space
C> @param WRK3 L3 sized scratch space
C> @param SCR  Scratch space for call to GLDIAG
C> @param IWRK Scratch space for call to GLDIAG
C> @param L0   Number of basis functions
C> @param L1   Number of basis functions after removing spherical
C>             contaminants and linear dependencies, i.e., total 
C>             number of orbitals (occupied + virtual)
C> @param L2   Upper/lower triangular portion of a L1 square matrix,
C>             i.e., L2 = (L1*L1+L1)/2
C> @param L3   L1 square matrix, ie., L3 = L1*L1
C> @param NOCC Total number of occupied orbitals (NCORE + VVAL) 
C>
      SUBROUTINE VVOS(VEC,FAO,EIG,WRK1,WRK2,WRK3,SCR,IWRK,
     *                L0,L1,L2,L3,NOCC)
      use mx_limits, only: mxgtot,mxsh,mxao,mxatm,mxrt,mxnoro,mxfrz
      use mod_sformas, only: nacvvo
      use constants, only: zero,one
      use comm_orntmo
      use comm_vvopar
C
      use omp_lib
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION METHOD
C
      LOGICAL MCP
C
      DIMENSION VEC(L1,L0),FAO(L2),EIG(L1),
     *          WRK1(L3),WRK2(L3),WRK3(L3),SCR(L1,8),IWRK(L1)
      DIMENSION NBFS(4),MINF(4),MAXF(4),NANGM(4)
C
      LOGICAL GOPARR,DSKWRK,MASWRK,OK,DZ,CANONC,FCORE,FORS,EKT,LINSER,
     *        QUAMBO,EXIT
C
C
      COMMON /DETWFN/ WSTATE(MXRT),SPINS(MXRT),CRIT,PRTTOL,SDET,SZDET,
     *                GRPDET,STSYM,GLIST,DWPARM,
     *                NFLGDM(MXRT),IWTS(MXRT),NCORSV,NCOR,NACT,NORBDT,
     *                NADET,NBDET,KDET,KSTDET,IROOT,IPURES,MAXW1,NITDET,
     *                MAXP,NCIDET,IGPDET,KSTSYM,NFTGCI,IDWEIGH,
     *                fstate(mxrt),ifts(mxrt)
      COMMON /ECP2  / CLP(400),ZLP(400),NLP(400),KFIRST(MXATM,6),
     *                KLAST(MXATM,6),LMAX(MXATM),LPSKIP(MXATM),
     *                IZCORE(MXATM)
      COMMON /FMCOM / XX(1)
      COMMON /GUGWFN/ NFZC,NMCC,NDOC,NAOS,NBOS,NALP,NVAL,NEXT,NFZV,
     *                IFORS,IEXCIT,ICICI,NOIRR
      COMMON /GVBWFN/ CICOEF(2,12),FGVB(25),ALPHA(325),BETA(325),NO(10),
     *                NCO,NSETO,NOPEN,NPAIR,NORBGVB,NCONF(MXAO),NHAM
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NEX,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /INTOPT/ ISCHWZ,IECP,NECP,IEFLD
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON /MCINP / METHOD,CISTEP,FINALCI,ACURCY,ENGTOL,DAMP,
     *                MICIT,NWORD,NORB,NOROT(2,MXNORO),MOFRZ(MXFRZ),
     *                NPFLG(10),NOFO,MCFMO,IDIABAT,
     *                CANONC,FCORE,FORS,EKT,LINSER
      COMMON /NSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),
     *                CF(MXGTOT),CG(MXGTOT),CH(MXGTOT),CI(MXGTOT),
     *                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),
     *                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /PRPOPT/ ETOLLZ,ILOCAL,IAHARD
      COMMON /RELWFN/ RMETHOD,QRQMT,CLIG,CLIG2,QRTOL,TAU,
     *                IQRORD,MODQR,NESOC,NRATOM,
     *                NUMU,NQMTR,NQRDAF,MORDA,NDARELB
      COMMON /RUNOPT/ RUNTYP,EXETYP,NEVALS,NGLEVL,NHLEVL
      COMMON /SYMTRY/ MAPSHL(MXSH,48),MAPCTR(MXATM,48),
     *                TT(432),INVT(48),NT
      COMMON /WFNOPT/ SCFTYP,VBTYP,DFTYPE,TDDFTYP,CITYP,CCTYP,
     *                MPLEVL,MPCTYP
C
C         THE PERIODIC TABLE ONLY GOES UP TO SPDF:
      DATA NBFS/ 1, 3, 6,10/
      DATA MINF/ 1, 2, 5,11/
      DATA MAXF/ 1, 4,10,20/
      DATA NANGM/1, 2, 3, 4/
      DATA CHECK/8HCHECK   /
      DATA RMC/8HMCSCF   /
      DATA RHF/8HRHF     /
      DATA ROHF/8HROHF    /
      DATA GVB/8HGVB     /
      DATA RNONE/8HNONE    /
      DATA GUGA,ALDET,GENCI,ORMAS,GMCCI
     *    /8HGUGA    ,8HALDET   ,8HGENCI   ,8HORMAS   ,8HGMCCI   /
C
C         ----- GENERATE VIRTUAL VALENCE ORBITAL SPACE -----
C       SEE PARTICULARLY SECTION II.F AND V OF THE FIRST PAPER:
C
C     MOLECULE INTRINSIC MINIMAL BASIS SETS. I. EXACT RESOLUTION OF
C     AB INITIO OPTIMIZED MOLECULAR ORBITALS IN TERMS OF DEFORMED
C     MINIMAL-BASIS ORBITALS
C     W.C.LU, C.Z.WANG, M.W.SCHMIDT, L.BYTAUTAS, K.M.HO, K.RUEDENBERG
C     J.CHEM.PHYS. 120, 2629-2637(2004)
C
C     MOLECULE INTRINSIC MINIMAL BASIS SETS. II. BONDING ANALYSIS
C     FOR SI4H6 AND SI2 TO SI10.
C     W.C.LU, C.Z.WANG, M.W.SCHMIDT, L.BYTAUTAS, K.M.HO, K.RUEDENBERG
C     J.CHEM.PHYS. 120, 2638-2651(2004)
C
C     REPRESENTATION OF ELECTRONIC STRUCTURES IN CRYSTALS IN TERMS OF
C     HIGHLY LOCALIZED QUASIATOMIC MINIMAL BASIS ORBITALS.
C     W.C.LU, C.Z.WANG, T.L.CHAN, K.RUEDENBERG, K.M.HO
C     PHYS.REV.B 70, 041101/1-4(2004)
C
C     A COMPREHENSIVE ANALYSIS OF MOLECULE-INTRINSIC QUASI-ATOMIC,
C     BONDING, AND CORRELATING ORBITALS.
C     I.  HARTREE-FOCK WAVE FUNCTIONS
C     A.C.WEST, M.W.SCHMIDT, M.S.GORDON, K.RUEDENBERG
C     J.CHEM.PHYS. 139, 234107/1-21(2013)
C
C     ON ENTRY, -FAO- CONTAINS THE RHF FOCK OPERATOR IN THE AO BASIS,
C     WHILE OTHER ARRAYS ARE TREATED AS WORK SPACES.
C
C------------------------------------------------
C     CHECK THE SCF ORBITALS FOR ORTHONORMALITY.
C
      IF(ILOCAL.EQ.4) THEN
         CALL DAREAD(IDAF,IODA,WRK2,L2,12,0)
         CALL LOCAL_CHECK_ORTHOG(
     *        L0,L1,L2,
     *        WRK1,WRK2,VEC,WRK3,
     *        MASWRK,IW)
      END IF
C------------------------------------------------
C     EXTREME MEASURE TO PRESERVE ISVDOP
C     ERROR OUTS BASED ON SOME NON-SENSE ARE BETTER
      ISVDOP2=ISVDOP
C------------------------------------------------
C
C     BE SURE THE RUN IS AN ALL E- BASIS SET, WITH ONLY H-Xe 
C     (NON-RELATIVISTIC) OR H-RN (RELATIVISTIC) BECAUSE THESE ARE
C     THE ONLY KIND OF ATOMIC MINIMAL BASIS SETS THAT ARE
C     INTERNALLY STORED FOR USE IN THE PROJECTION OF AOS.
C
CGS   Model core potentials also work for VVOS
C
      OK = .TRUE.
      DO I=1,NAT
         NUCZ = INT(ZAN(I) + 0.001D+00)
         IF(NUCZ.NE.IAN(I))                    OK=.FALSE.
         IF(NUCZ.GE.87)                        OK=.FALSE.
      ENDDO
C
C     Select relativistic AAMBS for use with MCPs
C
      MCP=.FALSE.
      IF(IECP.EQ.5 .OR. IECP.EQ.6) THEN
        OK=.TRUE.
        MCP=.TRUE.
      ENDIF
C
      IF(.NOT.MCP.AND.(IVVTYP.EQ.1.AND.RMETHOD.EQ.RNONE)) THEN
         IF(MASWRK) WRITE(IW,9006)
      ENDIF
C
      IF(.NOT.OK) THEN
         IF(MASWRK) WRITE(IW,9000)
         RETURN
      END IF
C
C     Get the number of valence orbitals
C
      NVALAT=0
      ISWMBS=1
      DO I=1,NAT
         NUMVAL1=LOCAL_NUMVAL(I,ISWMBS)
         NVALAT=NVALAT+NUMVAL1
      ENDDO
C
C     NCORE   = Total number of core orbitals (all atoms)
C     NOCC    = Number of occupied orbitals (core + occupied valence)
C     NOCCVAL = Number of occupied valence orbitals
C     NVIRVAL = Number of virtual valence orbitals (VVOs)
C     NVIR    = Number of remaining virtual orbitals
C     NALFA   = Number of singly occupied orbitals (ROHF)
C     NFILL   = Number of doubly occupied valence orbitals (ROHF)
C               that are not included in the active space (MCSCF)
C     MACT    = Number of active orbitals (MCSCF)
C
      NCORE   = NVVOS_NUMCOR(0,0)
      NOCCVAL = NOCC - NCORE
      NVIRVAL = NVALAT - NOCCVAL
      NVIR    = L0 - NOCC
      IF(SCFTYP.EQ.RHF) THEN
         IF(MASWRK) WRITE(IW,9010) NCORE,NVALAT,NOCCVAL,NVIRVAL,NVIR
      END IF
      IF(SCFTYP.EQ.ROHF) THEN
         NALFA=NA-NB
         NFILL=NOCCVAL-NALFA
         IF(MASWRK) WRITE(IW,9012) NCORE,NVALAT,NOCCVAL,NFILL,NALFA,
     *                             NVIRVAL,NVIR
      END IF
      IF(SCFTYP.EQ.GVB) THEN
         IF(MASWRK) WRITE(IW,9013) NCORE,NVALAT,NCO-NCORE,
     *                             NOCC-NCO,NVIRVAL,NVIR
      END IF
      IF(SCFTYP.EQ.RMC) THEN
         IF(CISTEP.EQ.ALDET .OR.
     *      CISTEP.EQ.GENCI .OR.
     *      CISTEP.EQ.ORMAS .OR.
     *      CISTEP.EQ.GMCCI) THEN
            MACT = NACT
         END IF
         IF(CISTEP.EQ.GUGA) THEN
            MACT = NDOC+NALP+NAOS+NBOS+NVAL
         END IF
         NFILL = NOCCVAL - MACT
         IF(MASWRK) WRITE(IW,9014) NCORE,NVALAT,NOCCVAL,NFILL,MACT,
     *                             NVIRVAL,NVIR
      END IF
C
      NVVOS = NVIRVAL
C
C         ESTIMATE GAUSSIANS AND SHELLS FOR BIGGEST ATOM IN EACH ROW.
      NERR   = 0
      MXGMBS = 0
      MXSHMBS= 0
      DO IAT=1,NAT
         IF(                       IAN(IAT).LE. 2) THEN
            MXGMBS  = MXGMBS  + 8
            MXSHMBS = MXSHMBS + 1
         END IF
         IF(IAN(IAT).GE.3  .AND.  IAN(IAT).LE.18) THEN
            IF(IVVTYP.EQ.0) THEN
              IF(IAN(IAT).GE. 3  .AND.  IAN(IAT).LE.10) THEN
                MXGMBS  = MXGMBS  + 14 + 14 + 7
                MXSHMBS = MXSHMBS + 3
              END IF
              IF(IAN(IAT).GE.11  .AND.  IAN(IAT).LE.18) THEN
                MXGMBS  = MXGMBS  + 18 + 18 + 12 + 18 + 12
                MXSHMBS = MXSHMBS + 5
              END IF
            ELSEIF(IVVTYP.EQ.1) THEN
              IF(IAN(IAT).GE. 3  .AND.  IAN(IAT).LE.10) THEN
                MXGMBS  = MXGMBS  + 20 + 20 + 13
                MXSHMBS = MXSHMBS + 3
              END IF
              IF(IAN(IAT).GE.11  .AND.  IAN(IAT).LE.18) THEN
                MXGMBS  = MXGMBS  + 24 + 24 + 17 + 24 + 17
                MXSHMBS = MXSHMBS + 5
              END IF
            ELSE
               NERR=NERR+1
            ENDIF
         END IF
C           IN NEXT TWO, NO THOUGHT WAS GIVEN TO T.M. CASES
C           NOR HAVE SLIGHTLY SMALLER BASES IN GROUPS 1/2 BEEN SET UP
         IF(IAN(IAT).GE.19  .AND.  IAN(IAT).LE.36) THEN
            IF(IVVTYP.EQ.0) THEN
               MXGMBS  = MXGMBS  + 26 + 26+20 + 26+20+14 + 26+20
               MXSHMBS = MXSHMBS + 8
            ELSEIF(IVVTYP.EQ.1) THEN
               MXGMBS  = MXGMBS  + 26 + 26+21 + 26+21+17 + 26+21
               MXSHMBS = MXSHMBS + 8
            ELSE
               NERR=NERR+1
            ENDIF
         END IF
         IF(IAN(IAT).GE.37  .AND.  IAN(IAT).LE.54) THEN
            IF(IVVTYP.EQ.0) THEN
              MXGMBS  = MXGMBS  + 28 + 28+23 + 28+23+17 + 28+23+17
     *                          + 28+23
              MXSHMBS = MXSHMBS + 11
            ELSEIF(IVVTYP.EQ.1) THEN
              MXGMBS  = MXGMBS  + 30 + 30+26 + 30+26+24 + 30+26+24
     *                          + 30+26
              MXSHMBS = MXSHMBS + 11
            ELSE
               NERR=NERR+1
            ENDIF
         END IF
         IF(IAN(IAT).GE.55  .AND.  IAN(IAT).LE.56) THEN
           MXGMBS = MXGMBS + 37 + 37+32 + 37+32+29 + 37+32+29
     *                     + 37+32+29 + 37+32
           MXSHMBS = MXSHMBS + 14
         ENDIF
         IF(IAN(IAT).GE.57  .AND.  IAN(IAT).LE.86) THEN
           MXGMBS = MXGMBS + 41 + 41+38 + 41+38+30 + 41+38+30+20
     *                     + 41+38+30 + 41+38
           MXSHMBS = MXSHMBS + 15
         ENDIF
      ENDDO
C
      IF(NERR.NE.0) THEN
         IF(MASWRK) WRITE(IW,9016)
         CALL FLSHBF(IW)
         CALL ABRT
         STOP
      ENDIF
C
      CALL VALFM(LOADFM)
      LEXMBS   = LOADFM   + 1
      LCSMBS   = LEXMBS   + MXGMBS
      LCPMBS   = LCSMBS   + MXGMBS
      LCDMBS   = LCPMBS   + MXGMBS
      LCFMBS   = LCDMBS   + MXGMBS
      LCSINP   = LCFMBS   + MXGMBS
      LCPINP   = LCSINP   + MXGMBS
      LCDINP   = LCPINP   + MXGMBS
      LCFINP   = LCDINP   + MXGMBS
      LINTYP   = LCFINP   + MXGMBS
      LKSTRMBS = LINTYP   + MXSHMBS
      LKATMMBS = LKSTRMBS + MXSHMBS
      LKTYPMBS = LKATMMBS + MXSHMBS
      LKNGMBS  = LKTYPMBS + MXSHMBS
      LKLOCMBS = LKNGMBS  + MXSHMBS
      LKMINMBS = LKLOCMBS + MXSHMBS
      LKMAXMBS = LKMINMBS + MXSHMBS
      LNSHLS   = LKMAXMBS + MXSHMBS
      LAST     = LNSHLS   + NAT
      NEED1 = LAST -LOADFM -1
      CALL GETFM(NEED1)
C
C     ---- SET UP NEAR-HF QUALITY MINIMAL ATOMIC BASIS, ON THE FLY ----
C     BESIDES SETTING UP DETAILS IN EXMBS, CSMBS, CPMBS, KXXXMBS, ...,
C     THE DIMENSIONALITY IS SPECIFIED BY L1MBS,NGAUMBS,NSHLMBS
C     NOTE THAT CSINP,CPINP,INTYP,NSHLS ARE STRICTLY WORKING STORAGE.
C     CAUTION, NAT NEEDS TO COUNT 1,2,3,... AS THE BASIS IS CONSTRUCTED.
C
      IERR1   = 0
      IERR2   = 0
      L1MBS   = 0
      NSHLMBS = 0
      NGAUMBS = 0
      DZ = .FALSE.
      NATSV = NAT
C
C     Load the AAMBS
C
      DO IAT=1,NATSV
         NAT  = IAT
         NUCZ = IAN(IAT)
         CALL ETGTO(NUCZ,DZ,XX(LCSINP),XX(LCPINP),XX(LCDINP),
     *              XX(LCFINP),IERR1,IERR2,L1MBS,
     *              XX(LINTYP),NANGM,NBFS,MINF,MAXF,XX(LNSHLS),
     *              XX(LEXMBS),XX(LCSMBS),XX(LCPMBS),XX(LCDMBS),
     *              XX(LCFMBS),
     *              XX(LKSTRMBS),XX(LKATMMBS),XX(LKTYPMBS),XX(LKNGMBS),
     *              XX(LKLOCMBS),XX(LKMINMBS),XX(LKMAXMBS),
     *              NGAUMBS,NSHLMBS,MXGMBS,MXSHMBS,IAT)
         IF(MASWRK  .AND.  IERR1.NE.0) WRITE(IW,9020) MXSHMBS
         IF(MASWRK  .AND.  IERR2.NE.0) WRITE(IW,9030) MXGMBS
         IF(IERR1.NE.0  .OR.  IERR2.NE.0) CALL ABRT
      ENDDO
C
      NAT = NATSV
      L2MBS = (L1MBS*L1MBS+L1MBS)/2
      L3MBS = L1MBS*L1MBS
C
      CALL VALFM(LOADFM)
      LOVLP = LOADFM + 1
      LSMBS = LOVLP  + L1*L1MBS
      LEVAL = LSMBS  + MAX(L2MBS,L3MBS)
      LEVEC = LEVAL  + L1MBS
      LSMHF = LEVEC  + L3MBS
      LAST  = LSMHF  + L3MBS
      NEED2 = LAST - LOADFM - 1
      CALL GETFM(NEED2)
C
      IF(EXETYP.EQ.CHECK) L0MBS=L1MBS
      IF(EXETYP.EQ.CHECK) GO TO 110
C
C        STEP (0), PRE-ORTHOGONALIZE THE LARGE GTO BASIS SET SPACE
C           a) get the metric for the large GTO basis set
C           b) transform it to spherical harmonics, after
C              which it is unit matrices on the diagonal,
C              but overlaps in between the atoms.
C           c) generate S**(-1/2) orthogonalizing transformation
C        We actually apply the Lowdin orthogonalization at step (I)!
C        This step was not included in the initial implementation
C        of 2007, from 2004 equations, but was added in 6/2013.
C
      MODE=0
      DUMMY=ZERO
      IDUMMY=1
      CALL COOVLP(MODE,XX(LSMBS),DUMMY,L1MBS,IDUMMY,L2MBS,NAT,
     *            NGAUMBS,NSHLMBS,XX(LEXMBS),XX(LCSMBS),XX(LCPMBS),
     *            XX(LCDMBS),XX(LCFMBS),DUMMY,DUMMY,DUMMY,
     *            XX(LKSTRMBS),XX(LKATMMBS),XX(LKTYPMBS),XX(LKNGMBS),
     *            XX(LKLOCMBS),XX(LKMINMBS),XX(LKMAXMBS),
     *            IDUMMY,IDUMMY,DUMMY,
     *            DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,
     *            IDUMMY,IDUMMY,IDUMMY,IDUMMY,IDUMMY,IDUMMY,IDUMMY,
     *            C,DUMMY)
C
      MODE=0
      DUMMY=ZERO
      CALL TFAO5D(MODE,XX(LSMBS),DUMMY,XX(LKTYPMBS),WRK1,
     *            WRK2,WRK3,L1,NSHLMBS,L0MBS,L1MBS,L2MBS)
C
C          get S**(-1/2) for later Lowdin orthogonalization
C
      CALL GLDIAG(L0MBS,L0MBS,L0MBS,XX(LSMBS),SCR,XX(LEVAL),
     *            XX(LEVEC),IERR,IWRK)
      CALL VCLR(XX(LSMBS),1,L3MBS)
      DO I=1,L0MBS
         XX(LSMBS-1+(I-1)*L0MBS+I) = ONE/SQRT(XX(LEVAL-1+I))
      ENDDO
      CALL TRPOSQ(XX(LEVEC),L0MBS)
      CALL TFSQU(XX(LSMHF),XX(LSMBS),XX(LEVEC),SCR,L0MBS,L0MBS)
C
C        STEP (I), GENERATION OF THE MATRIX "LITTLE-A-STAR", IN -WRK1-
C
C        THE IMPLEMENTATION HERE USES A FIXED, PREVIOUSLY STORED ATOMIC
C        BASIS SET AS THE CAPITAL-A-STAR ORBITALS, RATHER THAN TAKING
C        THESE FROM AN ATOMIC COMPUTATION USING THE CURRENT MOLECULE'S
C        BASIS SET.
C
C        GET OVERLAP BETWEEN THE AO MBS AND THE BASIS SET OF THE RUN,
C        COOVLP COMPUTES <BASIS2|BASIS1>, WITH SHAPE = OVLP(L1,L1MBS),
C
      MODE=1
      DUMMY=ZERO
      IDUMMY=1
      CALL COOVLP(MODE,DUMMY,XX(LOVLP),L1MBS,L1,IDUMMY,NAT,
     *            NGAUMBS,NSHLMBS,XX(LEXMBS),XX(LCSMBS),XX(LCPMBS),
     *            XX(LCDMBS),XX(LCFMBS),DUMMY,DUMMY,DUMMY,
     *            XX(LKSTRMBS),XX(LKATMMBS),XX(LKTYPMBS),XX(LKNGMBS),
     *            XX(LKLOCMBS),XX(LKMINMBS),XX(LKMAXMBS),
     *            MXGTOT,NSHELL,EX,CS,CP,CD,CF,CG,CH,CI,
     *            KSTART,KATOM,KTYPE,KNG,KLOC,KMIN,KMAX,
     *            C,C)
C
C        POSSIBLY, TRANSFORM AUXILIARY BASIS FROM CARTESIAN TO
C        SPHERICAL HARMONICS, REDUCING -L1MBS- TO -L0MBS-.
C
      MODE=1
      DUMMY=ZERO
      CALL TFAO5D(MODE,DUMMY,XX(LOVLP),XX(LKTYPMBS),
     *            WRK1,WRK2,DUMMY,L1,NSHLMBS,L0MBS,L1MBS,L2MBS)
C
C        SAVE THE (L1,L0MBS) "SV" MATRIX FOR POSSIBLE USE IN
C        SUBSEQUENT INTERNAL SPACE SVD STAGES.
C        NOTE:  WORKS WHETHER L1MBS=L0MBS OR NOT.
C
      CALL DAWRIT(IDAF,IODA,XX(LOVLP),L1*L0MBS,533,0)
C
C          rotate the MBS side to the orthogonal AO space A* here.
C          The behavior of this code between 2007 and 2013 can be
C          recovered by avoiding the next Lowdin orthogonalization.
C
      CALL DGEMM('N','N',L1,L0MBS,L0MBS,
     *           ONE,XX(LOVLP),L1,XX(LSMHF),L0MBS,ZERO,WRK1,L1)
      CALL DCOPY(L1*L0MBS,WRK1,1,XX(LOVLP),1)
C
C        TRANSFORM FROM TRUE AO SPACE TO THE SCF MO BASIS (OCC+VIRT).
C        -WRK1- IS OVERLAP BETWEEN SCF ORBITALS AND THE AUXILIARY
C        SPHERICAL HARMONIC MINIMAL ATOMIC BASIS SET.
C
      CALL DGEMM('T','N',L0,L0MBS,L1,
     *           ONE,VEC,L1,XX(LOVLP),L1,ZERO,WRK1,L0)
C
  110 CONTINUE
      CALL RETFM(NEED2)
      CALL RETFM(NEED1)
C
C        EARLY EXITS:
C        NOTE:  YOU CAN REALLY ONLY FIGURE THIS OUT BY
C               RUNNING ...ALL... ACCUMULATED EXAMPLES.
C        --IVVOS=0 RUNS.
C          ALL RUNS REQUIRE DAF533 IN THE ORIENTED ORBITAL SECTION.
C          BEFORE, ISVDOP=0 FULL-VAL. MCSCF RUNS DID NOT USE VVOS.
C          ALSO,   ALDECI RUNS HAVE LATER ERROR OUTS FOR VVOS.
C          ALSO,   THE VVOS LOGICAL RESTRICTS CALLS TO VVOS
C                  BASED ON RUNTYPS,ETC, WHICH PREVENTS CONFLICTS.
C        --SVDCOPT
C        --CCTYP   is presently incompatible with LOCAL=SVD.
C
      IF(ISVDOP.EQ.0.AND.NVIRVAL.EQ.0
     *              .AND.SCFTYP.EQ.RMC.AND.EXTLOC.EQ.RNONE) IVVOS=0
      EXIT=.FALSE.
      IF(IVVOS.EQ.0) EXIT=.TRUE.
      IF(SVDCOPT)    EXIT=.TRUE.
      IF(CCTYP.NE.RNONE  .AND.  ILOCAL.EQ.4) EXIT=.TRUE.
      IF(EXIT) THEN
        IF(MASWRK) WRITE(IW,9035)
        CALL FLSHBF(IW)
        RETURN
      ENDIF
      IF(EXETYP.EQ.CHECK) GO TO 210
C
C        STEP (II), FORM CAPITAL-B (IN -WRK2-) AND DIAGONALIZE IT
C                   TO OBTAIN ITS -P- (NVIRVAL) LARGEST VALUES
C
C        Note that we are squaring an SVD equation, to convert it
C        into an ordinary eigenvalue problem (B is symmetric).
C        The order of the squaring means we solve for the rotation
C        only on the molecular side of the SVD, since we don't
C        intend to do anything with the atom's side of the SVD.
C
CCC no point to parallelize AND blas call inside
      IJ=0
      DO I=1,NVIR
         DO J=1,I
            IJ=IJ+1
            WRK2(IJ) = DDOT(L0MBS,WRK1(NOCC+I),L0,WRK1(NOCC+J),L0)
         ENDDO
      ENDDO
C         DIAGONALIZATION SHOULD GET EIGENVALUES IN DESCENDING ORDER
      NVIR2 = (NVIR*NVIR+NVIR)/2
CCC possible to parallelize, but very tricky, and can lose vectorization
CCC does not compile with single or master with old gfortran,
CCC ok with critical, but it is an unnecessary work
CCC leave as is for now
      DO I=1,NVIR2
         WRK2(I) = -WRK2(I)
      ENDDO
CCC need modern compiler for modern OpenMP version to support vectorization!
      CALL GLDIAG(NVIR,NVIR,NVIR,WRK2,SCR,EIG,WRK3,IERR,IWRK)
      DO I=1,NVIR
         EIG(I) = -EIG(I)
      ENDDO
C
C        SHOW THE USER THE BREAK IN EIGENVALUES FOR -B-.
C        print SQRT to correspond to the SVD problem, not its square.
C
      MINPRT=MAX(   1,NVIRVAL-4)
      MAXPRT=MIN(NVIR,NVIRVAL+5)
      IF(MASWRK) THEN
         WRITE(IW,9040)
         WRITE(IW,9050)
         WRITE(IW,9070) (SQRT(EIG(III)),III=MINPRT,NVIRVAL)
         WRITE(IW,9060)
         WRITE(IW,9070) (SQRT(EIG(III)),III=NVIRVAL+1,MAXPRT)
      END IF
C
C        Create QUAMBOs?
C
  210 CONTINUE
      QUAMBO=.TRUE.
C
C        The mathematics in the 1st 2004 paper is now known to be
C        that of the Singular Value Decomposition.  The overlap
C        between the pre-orthogonalized AAMBS orbitals A* indexed 'j'
C        and the SCF virtuals 'v' is the little-a-star generated at
C        step (I).  Its SVD would be
C           T-dagger(w,v) * a*(v,j) * U(j,k) = lambda-w * delta(w,k)
C        The 2004 paper chooses to square this equation, generating
C        the ordinary eigenvalue problem T-dagger * B * T = lambda**2,
C        for B= sum-on-j a*(v,j)a*(w,j), which avoids the need to call
C        an SVD routine here.
C
C        If our only concern is creating the pseudo-canonical VVO space,
C        we can rotate the SCF virtuals by T(v,w) (at -WRK3-) and thus
C        entirely skip QUAMBO projection.  The first few columns of T
C        correspond to the largest singular values, and the remaining
C        columns span the external space.  Since there are many zero
C        eigenvalues in the external space (high degeneracy for the
C        diagonalization), those eigenvectors might be a little dodgy?
C        We therefore choose not to branch around the procedure to
C        generate an external space by orthogonality to the internals.
C        Of course, we must still canonicalize, as a way to clean up
C        both the internal and external virtual spaces.
C
      IF(.NOT.QUAMBO) THEN
         IF(MASWRK) WRITE(IW,9200)
         CALL DGEMM('N','N',L1,NVIR,NVIR,1.0D+00,VEC(1,NOCC+1),L1,
     *              WRK3,NVIR,0.0D+00,WRK2,L1)
         CALL DCOPY(L1*NVIR,WRK2,1,VEC(1,NOCC+1),1)
         CALL DCOPY(L1*L0,VEC,1,WRK1,1)
         GO TO 600
      END IF
C
C        STEP (III), FORM CAPITAL-R, CAPITAL-D, BOTH LITTLE-A-PRIME'S
C
      CALL VALFM(LOADFM)
      LR    = LOADFM + 1
      LD    = LR     + NVIR*NVIR
      LWRK4 = LD     + L0MBS
      LWRK5 = LWRK4  + NVIR*L0MBS
      LAST  = LWRK5  + L0MBS*L0MBS
      NEED = LAST-LOADFM-1
      CALL GETFM(NEED)
      IF(EXETYP.EQ.CHECK) GO TO 800
C
C          NOTE THAT THE PROJECTOR R FORMED HERE SATISFIES R*R=R
C
      CALL DGEMM('N','T',NVIR,NVIR,NVIRVAL,
     *           ONE,WRK3,NVIR,WRK3,NVIR,ZERO,XX(LR),NVIR)
C
C          COMPUTE 2ND TERM OF D BEFORE THE 1ST TERM,
C          NOTE REUSE OF THE INTERMEDIATE -WRK4- IN THE SECOND LOOP.
C
      CALL DGEMM('N','N',NVIR,L0MBS,NVIR,
     *           ONE,XX(LR),NVIR,WRK1(NOCC+1),L0,ZERO,XX(LWRK4),NVIR)
      CALL DGEMM('T','N',L0MBS,1,NVIR,
     *           ONE,WRK1(NOCC+1),L0,XX(LWRK4),NVIR,ZERO,XX(LD),L0MBS)
      DO J=1,L0MBS
         LW1 = L0*(J-1) + 1
         XX(LD+J-1) = XX(LD+J-1) + DDOT(NOCC,WRK1(LW1),1,WRK1(LW1),1)
      ENDDO
C
C        IF THERE IS NO OVERLAP BETWEEN AN ATOMIC ORBITAL AND
C        THE OCCUPIED SPACE DUE TO SYMMETRY REASONS, AND FURTHER
C        THERE ARE ENOUGH ZERO ELEMENTS IN -R-, IT MAY HAPPEN
C        THAT AN ELEMENT OF -D- IS COMPUTED AS ZERO.  THIS
C        HAPPENS FOR CH2.  IT CAN BE PATCHED OVER BY USING THE
C        DIRECT PROJECTION OF THE AO ONTO THE VIRTUAL SPACE
C        AS THE 'PRE-NORMALIZATION'.
C        In 2014, a similar problem was found in HSiPH2, which
C        has no symmetry.  Coordinates with the Si empty p
C        function oriented almost exactly in the y direction
C        were found to have a negative value for Dj below for
C        the valence py QUAMBO.  Because skipping this pre-
C        normalization can be compensated for during the
C        following symmetric orthonormalization, a "fix" is
C        taking such normalizations as "unity".
C
      DO J=1,L0MBS
         LW1 = L0*(J+NOCC-1) + 1
         IF(ABS(XX(LD-1+J)).LT.1.0D-09)
     *          XX(LD-1+J) = DDOT(NVIR,WRK1(LW1),1,WRK1(LW1),1)
         IF(XX(LD-1+J).LT.ZERO) XX(LD-1+J)=ONE
      ENDDO
C
C          BUILD LITTLE-A-PRIME ON TOP OF LITTLE-A-STAR, IN -WRK1-
C
      DO J=1,L0MBS
         LW1 =   L0*(J-1) + 1
         LW4 = NVIR*(J-1) + LWRK4
         DFACT = ONE/SQRT(XX(LD+J-1))
         CALL DSCAL(NOCC,DFACT,WRK1(LW1),1)
         CALL DSCAL(NVIR,DFACT,  XX(LW4),1)
         CALL DCOPY(NVIR,XX(LW4),1,WRK1(LW1+NOCC),1)
      ENDDO
C
C        STEP (V), SYMMETRIC ORTHOGONALIZATION ACROSS ALL ATOMS,
C                  THUS PRODUCING THE QUAMBO'S, A''.
C
C          FORM OVERLAP MATRIX OVER THE A-PRIME ORBITALS
C
      CALL DGEMM('T','N',L0MBS,L0MBS,NOCC,
     *           ONE,WRK1,L0,WRK1,L0,ZERO,XX(LWRK5),L0MBS)
      CALL DGEMM('T','N',L0MBS,L0MBS,NVIR,
     *           ONE,WRK1(NOCC+1),L0,WRK1(NOCC+1),L0,
     *           ONE,XX(LWRK5),L0MBS)
C
C          FORM S**(-1/2) MATRIX, IN -WRK3-
C
      CALL CPYSQT(XX(LWRK5),WRK2,L0MBS,1)
      CALL GLDIAG(L0MBS,L0MBS,L0MBS,WRK2,SCR,EIG,WRK3,IERR,IWRK)
      L2MBS = (L0MBS*L0MBS+L0MBS)
      CALL VCLR(WRK2,1,L2MBS)
      II = 0
      DO I=1,L0MBS
         II=II+I
         WRK2(II) = ONE/SQRT(EIG(I))
      ENDDO
      CALL TRPOSQ(WRK3,L0MBS)
      CALL TFTRI(XX(LWRK5),WRK2,WRK3,SCR,L0MBS,L0MBS,L0MBS)
      CALL CPYTSQ(XX(LWRK5),WRK3,L0MBS,1)
C
C          COMPUTE A'' COEFFICIENTS, A'' = LITTLE-A-PRIME * S**(-1/2)
C          AND THEN COPY THEM TO -WRK1- (OVERWRITING LITTLE-A-PRIME).
C
      CALL DGEMM('N','N',L0,L0MBS,L0MBS,
     *           ONE,WRK1,L0,WRK3,L0MBS,ZERO,WRK2,L0)
      CALL DCOPY(L0*L0MBS,WRK2,1,WRK1,1)
C
C          TRANSFORM A'' FROM CANONICAL MO BASIS TO AO SPACE
C          AND THEN COPY THEM TO -WRK1-
C
      CALL DAREAD(IDAF,IODA,VEC,L3,15,0)
      CALL DGEMM('N','N',L1,L0MBS,L0,
     *           ONE,VEC,L1,WRK1,L0,ZERO,WRK2,L1)
      CALL DCOPY(L1*L0MBS,WRK2,1,WRK1,1)
C
C          DUMP THE QUAMBO'S TO THE PUNCH FILE
C
      IF(MASWRK) THEN
         WRITE(IW,9080) L0MBS,NOCC,NVIRVAL
         WRITE(IP,8000)
         WRITE(IP,8010) ' $VEC'
         CALL PUSQL(WRK1,L0MBS,L1,L1)
         WRITE(IP,8010) ' $END'
      END IF
C
C        STEP (VI), FORM CANONICAL VALENCE ORBITALS
C          THE QUAMBO'S ARE IN WRK1, AND WE ARE MISSING THE
C          EXTERNAL SET.  FORTUNATELY, THE -ORTHO- ROUTINE
C          CAN GENERATE THESE BY ORTHOGONALIZING, THUS
C          GENERATING AN EXTERNAL VIRTUAL SPACE.  WE SAVE
C          THE QUAMBO'S AND THIS RATHER ARBITRARY EXTERNAL SET.
C
  600 CONTINUE
C
C          -S- MATRIX AT WRK2, -Q- MATRIX AT WRK3
      CALL DAREAD(IDAF,IODA,WRK2,L2,12,0)
      CALL DAREAD(IDAF,IODA,WRK3,L3,45,0)
C           WE SEEK ORBITALS ORTHOGONAL TO THE -L0MBS- QUAMBO'S
      CALL ORTHO(WRK3,WRK2,WRK1,SCR,L0MBS,L0,L1,L2,L1)
      CALL TFSQB(WRK1,WRK3,SCR,L0,L1,L1)
C------
C           PERFORM 2ND ORTHOGONALIZATION FOR NUMERICAL REASONS.
      CALL DAREAD(IDAF,IODA,WRK2,L2,12,0)
      CALL DAREAD(IDAF,IODA,WRK3,L3,45,0)
      CALL ORTHO(WRK3,WRK2,WRK1,SCR,L0,L0,L1,L2,L1)
      CALL TFSQB(WRK1,WRK3,SCR,L0,L1,L1)
C------
      CALL DAWRIT(IDAF,IODA,WRK1,L3,470,0)
C
C
C       ----- TWO STEP GENERATION OF "SEMI-CANONICAL" ORBITALS -----
C          WE WANT TO DIAGONALIZE THE DIAGONAL BLOCKS OF F, ONCE
C          IN THE QUAMBO SPACE (OCCUPIED + VALENCE VIRTUAL) AND
C          ONCE IN THE EXTERNAL VIRTUAL SPACE.  WE SEEK TO HAVE
C          THE EIGENVALUES ORDERED IN EACH SPACE, BUT NOT TO MIX
C          VALENCE VIRTUALS AND EXTERNAL VIRTUALS.  WE ALSO WANT
C          THE DIAGONALIZATION TO EXPLOIT FULL POINT GROUP SYMMETRY.
C
C      --- CANONICALIZE BY FREE DIAGONALIZATION WITHIN EACH SPACE ---
C          TRANSFORM FOCK MATRIX TO QUAMBO+EXTERNAL ORBITAL SPACE,
C          THEN DIAGONALIZE THE SEPARATE BLOCKS, AND FINALLY
C          TRANSFORM THE VECTORS BACK TO THE ATOMIC BASIS SET.
C          VEC TO RECEIVE SEMI-CANONICAL QUAMBO/EXTERNAL SPACE
C
C          in case there is no external space, skip separating it.
      IF(L0.EQ.L0MBS) GO TO 700
C
      CALL TFTRI(WRK3,FAO,WRK1,SCR,L0,L1,L1)
      CALL GLDIAG(L1,L0MBS,L0MBS,WRK3,SCR,EIG,VEC,IERR,IWRK)
      CALL TFSQB(VEC,WRK1,SCR,L0MBS,L1,L1)
CCC don't parallelize, unless very good with OpenMP
      IJ=0
      DO I=L0MBS+1,L0
         IROW = (I*I-I)/2
         DO J=L0MBS+1,I
            IJ=IJ+1
            KL=IROW+J
            WRK3(IJ) = WRK3(KL)
         ENDDO
      ENDDO
      CALL GLDIAG(L1,L0-L0MBS,L0-L0MBS,WRK3,SCR,EIG(L0MBS+1),
     *            VEC(1,L0MBS+1),IERR,IWRK)
      CALL TFSQB(VEC(1,L0MBS+1),WRK1(L0MBS*L1+1),SCR,L0-L0MBS,L1,L1)
C
C         For GVB, the internal space diagonalization above will
C         not have regenerated the converged occupied orbitals,
C         but will have made an external space orthogonal to the
C         QUAMBO space.
C         Put the occupied orbitals back in, then generate the
C         unoccupied internal orbital space by orthogonality.
C         For technical reasons, the orthogonality call carries over
C         to the external space as well (it seemingly can't be skipped).
C         We are assuming that the QUAMBO space's higher eigenvectors
C         will always span the unoccupied internal orbital space.
C         It is a fairly safe assumption since all GVB cases now
C         supply the "combined Fock matrix" to the VVO procedure.
C
      IF(SCFTYP.EQ.GVB) THEN
        CALL DAREAD(IDAF,IODA,WRK3,L3,15,0)
        CALL DCOPY(L1*NOCC,WRK3,1,VEC,1)
        CALL DAREAD(IDAF,IODA,WRK3,L3,45,0)
        CALL DAREAD(IDAF,IODA,WRK2,L2,12,0)
        CALL ORTHO(WRK3,WRK2,VEC,SCR,L0MBS,L0,L1,L2,L1)
        CALL TFSQB(VEC,WRK3,SCR,L0,L1,L1)
      END IF
C
C      --- FULL SYMMETRY SEMI-CANONICALIZATION ---
C          THE ABOVE HAS PRODUCED ORBITALS ORDERED BY THEIR EIGENVALUE
C          WITHIN EACH SPACE, BUT LACKS CONTROL OF THE ORBITAL SYMMETRY.
C     A) BRING F TO PREVIOUSLY CANONICALIZED QUAMBO/EXTERNAL SPACE
C     B) ZAP ALL ELEMENTS CONNECTING THESE SPACES
C     C) TRANSFORM THIS MODIFIED FOCK OPERATOR BACK TO AO SPACE
C     D) DIAGONALIZE MODIFIED F WITH FULL Q-MATRIX SYMMETRY CONTROL
C     E) MAKE SURE THE SPACES DON'T GET MIXED UP BY OVERLAP MAXIMIZATION
C
C     THE RESULT OF THIS SHOULD BE THE FOLLOWING ORDER OF ORBITALS:
C            ORIGINAL OCCUPIED CANONICAL ORBITALS
C            VALENCE VIRTUAL ORBITALS
C            EXTERNAL VIRTUAL ORBITALS
C
  700 CONTINUE
      IF(NT.EQ.1) GO TO 800
C
C         -F IN SEMI-CANONICAL QUAMBO/EXTERNAL SPACE- AT WRK3
C         ONLY THE VIRTUAL VALENCE/EXTERNAL VIRTUAL BLOCK SHOULD
C         BE NONZERO, BUT ZAP OCCUPIED/EXTERNAL BLOCK TOO.
C
      IF(SCFTYP.EQ.RHF) THEN
        IF(MASWRK) WRITE(IW,9084)
        CALL DAREAD(IDAF,IODA,FAO,L2,14,0)
      ELSE IF(SCFTYP.EQ.ROHF) THEN
        IF(MASWRK) WRITE(IW,9086)
        CALL DAREAD(IDAF,IODA,FAO,L2,51,0)
      ELSE IF(SCFTYP.EQ.GVB) THEN
        IF(MASWRK) WRITE(IW,9087)
        CALL DAREAD(IDAF,IODA,FAO,L2,51,0)
      ELSE IF(SCFTYP.EQ.RMC) THEN
        IF(MASWRK) WRITE(IW,9082)
        CALL DAREAD(IDAF,IODA,FAO,L2,532,0)
      ELSE
        IF(MASWRK) WRITE(IW,9088)
        CALL FLSHBF(IW)
        CALL ABRT
        STOP
      ENDIF
C
      CALL TFTRI(WRK3,FAO,VEC,SCR,L0,L1,L1)
!$omp parallel do private(i,j) schedule(dynamic)
      DO I=L0MBS+1,L0
         I0 = (I*I-I)/2
         DO J=1,L0MBS
            WRK3(I0+J) = ZERO
         ENDDO
      ENDDO
!$omp end parallel do
C
C        ENSURE NO MIXING INSIDE THE GVB OCCUPIED ORBITAL SPACE,
C        OR BETWEEN OCCUPIED SPACE AND EMPTY VVO SPACE.
C        THE 'COMBINED FOCK OPERATOR' HAS SMALL ELEMENTS DUE TO
C        CONVERGENCE NOT BEING EXACT.  ZEROING THEM OUT HAS THE
C        EFFECT OF NOT LETTING THE OCCUPIED ORBITALS CHANGE HERE,
C        AT ALL.
C
      IF(SCFTYP.EQ.GVB) THEN
!$omp parallel
!$omp do schedule(dynamic) private(i,j)
         DO I=1,NOCC
            I0 = (I*I-I)/2
            DO J=1,I-1
               WRK3(I0+J) = ZERO
            ENDDO
         ENDDO
!$omp end do nowait
!$omp do schedule(dynamic) private(i,j)
         DO I=NOCC+1,L0MBS
            I0 = (I*I-I)/2
            DO J=1,NOCC
               WRK3(I0+J) = ZERO
            ENDDO
         ENDDO
!$omp end do
!$omp end parallel
      END IF
C
C         -S- AT WRK2, -BLOCKED F IN MO SPACE- AT WRK3, SCRATCH AT WRK1
      CALL TFTRIB(FAO,WRK3,WRK2,VEC,WRK1,SCR,L0,L1,L2,L3)
C         STORE COPY OF ORBITALS AT -WRK1- FOR OVERLAP SELECTION BELOW
      CALL DCOPY(L0*L1,VEC,1,WRK1,1)
C         -Q- AT WRK3 USED TO SYMMETRY BLOCK DIAG STEP OF ZAPPED -FAO-
C         SCRATCH AT WRK2
      CALL DAREAD(IDAF,IODA,WRK3,L3,45,0)
      CALL TFTRI(WRK2,FAO,WRK3,SCR,L0,L1,L1)
      CALL SYMDIA(WRK2,VEC,EIG,SCR,IWRK,L0,L2,L1)
      CALL TFSQB(VEC,WRK3,SCR,L0,L1,L1)
C         ALIGN BY MAXIMUM OVERLAP TO FIRST STEP'S ORBITAL SPACE
C         -WRK1- IS ORIGINAL ORDER, -S- AT WRK2, SCRATCH AT WRK3
      CALL DAREAD(IDAF,IODA,WRK2,L2,12,0)
      IDUM=0
      CALL OVLSEL(VEC,EIG,WRK1,WRK2,WRK3,SCR,IWRK,SCR,
     *            IDUM,IDUM,L0,L1,L2,.FALSE.,L0)
C
C         ALL DONE, SAVE ORBITALS AND SEMI-CANONICAL EIGENVALUES.
C         NOTES:
C         1. WE DO NOT NEED TO ASSIGN ORBITAL SYMMETRY LABELS,
C            WHICH IS TAKEN CARE OF IN THE HARTREE-FOCK DRIVER.
C         2. FOR IVVOS=0 AND ISVDOP=10 RUNS, WE EXIT OUT EARLY
C            TO PRESERVE THE TRUE CANONICAL MCSCF ORBITALS.
C
  800 CONTINUE
      IF(QUAMBO) CALL RETFM(NEED)
      IF(L0.GT.L0MBS) THEN
         IEXTBEG = NOCC+NVIRVAL+1
         IEXTEND = L0
      ELSE
         IEXTBEG = 0
         IEXTEND = 0
      ENDIF
      IF(MASWRK) WRITE(IW,9090) NOCC,NVIR,
     *           1,NOCC, NOCC+1,NOCC+NVIRVAL, IEXTBEG,IEXTEND
      nacvvo=nocc+nvirval
C
C     NOTES FOR ISVDOP=10 OPTION:
C     EXIT EARLY IF NOT INSERTING VVOS ORBITALS.
C     CONTINUE ON OTHERWISE...
C
      IF(IVVOS.EQ.0.AND.ISVDOP.EQ.10) THEN
        IF(MASWRK) WRITE(IW,9095)
        RETURN
      ENDIF
C
C         GVB should just recover occupied orbitals and values, so
C         only VVOs and external are generated/canonicalized above.
C         Enough care has been taken that the occupied orbitals from
C         above are OK, so their reinsertion is "extra cautious",
C         but we don't have the original GVB program's eigenvalues.
C
      IF(SCFTYP.EQ.GVB) THEN
        CALL DAREAD(IDAF,IODA,WRK1,L3,15,0)
        CALL DAREAD(IDAF,IODA,WRK2,L1,17,0)
        CALL DCOPY(L1*NOCC,WRK1,1,VEC,1)
        CALL DCOPY(   NOCC,WRK2,1,EIG,1)
      END IF
C
C         The ordering between canonicals and vec orbital set
C         is not simple.  Using L0 in the ovlsel can result in vvos
C         that leave the vvos indices.  Using L0LIM often works.
C         However, sometimes the canonical orbitals in L0LIM have
C         zero overlap with the vec orbital set in l0lim.
C         As a result, the reordering instructions that come out of
C         ovlsel will still have 0s in some cases.
C         The following keeps 1-->nocc as definitive while allowing
C         some flexibility in the vvos reordering.
C
      IF(SCFTYP.EQ.RMC) THEN
        CALL DAREAD(IDAF,IODA,WRK1,L3,15,0)
        CALL DAREAD(IDAF,IODA,WRK2,L2,12,0)
C       CREATE PROPER REORDERING INSTRUCTIONS.
C       NOTE:  OVLSEL CANNOT BE USED HERE.
        L0LIM=NCORE+NVALAT
CKEEP        IDUM=0
CKEEP        CALL OVLSEL(VEC,EIG,WRK1,WRK2,WRK3,SCR,IWRK,SCR,
CKEEP     *              IDUM,IDUM,L0,L1,L2,.FALSE.,L0LIM)
        CALL VCLR(WRK3,1,L3)
        CALL MTARBR(WRK2,L1,VEC,L0,WRK3,L1,1)
        CALL VCLR(WRK2,1,L3)
        CALL MRTRBR(WRK1,L1,L1,L0,WRK3,L1,L0,WRK2,L1)
        CALL LOCAL_MAXCOIN(1,L0LIM,IWRK,WRK2,L1)
        DO 1000 I=NCORE+1,L0LIM
           DO J=NCORE+1,L0LIM
              IF(IWRK(J).EQ.I) GO TO 1000
           ENDDO
           DO J=NCORE+1,L0LIM
              IF(IWRK(J).EQ.0) THEN
                 IWRK(J)=I
                 GO TO 1000
              ENDIF
           ENDDO
 1000   CONTINUE
        CALL ICOPY(L0LIM,IWRK,1,SCR,1)
        CALL REORDR(EIG,SCR,L0LIM,1)
        CALL ICOPY(L0LIM,IWRK,1,SCR,1)
        CALL REORDR(VEC,SCR,L0LIM,L1)
C       COPY BACK OVER OCCUPIED INFORMATION.
        CALL DAREAD(IDAF,IODA,WRK1,L3,15,0)
        CALL DCOPY(L1*NOCC,WRK1,1,VEC,1)
        CALL VCLR(WRK1,1,L3)
        CALL DAREAD(IDAF,IODA,WRK1,L1,17,0)
        CALL DCOPY(NOCC,WRK1,1,EIG,1)
        IF(NOCC.EQ.L0LIM.AND.MASWRK) WRITE(IW,9096)
      ENDIF
C
      CALL DAWRIT(IDAF,IODA,VEC,L3,15,0)
      CALL DAWRIT(IDAF,IODA,EIG,L1,17,0)
C
C     FOREVER SAVE VVOS RESULTS.
      CALL DAWRIT(IDAF,IODA,VEC,L3,538,0)
C
C     THIS NEXT WARNING IS ADDED IN CASE SOMEONE THINKS
C     THEY ARE SMART BY PUTTING IVVOS=1 WHEN ISVDOP=10.
      IF(IVVOS.EQ.0) THEN
        IF(MASWRK) WRITE(IW,9100)
        CALL FLSHBF(IW)
        IF(ISVDOP.EQ.10.OR.ISVDOP2.NE.ISVDOP) THEN
          CALL ABRT
          STOP
        ENDIF
      ENDIF
C
      RETURN
C
 8000 FORMAT(/'SYMMETRICALLY ORTHOGONALIZED QUAMBO''S,',
     *        ' EXPANDED IN THE AO BASIS SET'/
     *        'NOTE THAT QUAMBO''S ARE NOT VVO''S,',
     *        ' THOSE ARE GIVEN IN THE NEXT ORBITAL GROUP')
 8010 FORMAT(A5)
C
 9000 FORMAT(1X,'ABANDONING -VVO- GENERATION...'/
     *       1X,'THE NECESSARY ATOMIC AO SETS ONLY EXIST FOR H-RN,'/
     *       1X,'AND/OR ECPS MAY NOT BE USED.')
 9006 FORMAT(/1X,'***WARNING: A RELATIVISTIC METHOD SHOULD BE USED ',
     *           'WITH IVVTYP=1.')
C 9007 FORMAT(/1X,'THE RELATIVISTIC AAMBS MUST BE USED FOR ELEMENT ',I3,
C     *       /1X,'RESETTING IVVTYP=1')
 9010 FORMAT(/1X,'FORMING VIRTUAL VALENCE ORBITALS FOR RHF...'/
     *       1X,'THE ATOMS POSSESS',I5,' CHEMICAL CORES, AND',I5,
     *          ' VALENCE ORBITALS.'/
     *       1X,'THIS MOLECULE HAS',I5,' OCCUPIED VALENCE ORBITALS,'/
     *       1X,'SO WE ARE SEEKING',I5,'  VIRTUAL VALENCE ORBITALS',
     *          ' FROM'/
     *       1X,'WITHIN THE MOLECULE''S VIRTUAL SPACE OF DIMENSION',I5)
 9012 FORMAT(/1X,'FORMING VIRTUAL VALENCE ORBITALS FOR ROHF...'/
     *       1X,'THE ATOMS POSSESS',I5,' CHEMICAL CORES, AND',I5,
     *          ' VALENCE ORBITALS.'/
     *       1X,'THIS MOLECULE HAS',I5,' OCCUPIED VALENCE ORBITALS',
     *          ' (',I5,' FILLED,',I3,' ALPHA),'/
     *       1X,'SO WE ARE SEEKING',I5,'  VIRTUAL VALENCE ORBITALS',
     *          ' FROM'/
     *       1X,'WITHIN THE MOLECULE''S VIRTUAL SPACE OF DIMENSION',I5)
 9013 FORMAT(/1X,'FORMING VIRTUAL VALENCE ORBITALS FOR GVB...'/
     *       1X,'THE ATOMS POSSESS',I5,' CHEMICAL CORES, AND',I5,
     *          ' VALENCE ORBITALS.'/
     *       1X,'THIS MOLECULE HAS',I5,'      FILLED VALENCE ORBITALS,'/
     *       1X,'              AND',I5,' PARTLY OCC. VALENCE ORBITALS,'/
     *       1X,'SO WE ARE SEEKING',I5,'     VIRTUAL VALENCE ORBITALS',
     *          ' FROM'/
     *       1X,'WITHIN THE MOLECULE''S VIRTUAL SPACE OF DIMENSION',I5)
 9014 FORMAT(/1X,'FORMING VIRTUAL VALENCE ORBITALS FOR MCSCF...'/
     *       1X,'THE ATOMS POSSESS',I5,' CHEMICAL CORES, AND',I5,
     *          ' VALENCE ORBITALS.'/
     *       1X,'THIS MOLECULE HAS',I5,' OCCUPIED VALENCE ORBITALS',
     *          ' (',I5,' FILLED,',I3,' ACTIVE),'/
     *       1X,'SO WE ARE SEEKING',I5,'  VIRTUAL VALENCE ORBITALS',
     *          ' FROM'/
     *       1X,'WITHIN THE MOLECULE''S VIRTUAL SPACE OF DIMENSION',I5)
 9016 FORMAT(/1X,'VVOS ERROR:'
     *       /1X,'PROBLEM IN GAUSSIANS AND SHELLS MAX DEFINITIONS.')
 9020 FORMAT(1X,'TOO MANY SHELLS IN VVO ATOMIC MINIMAL BASIS, >',I8)
 9030 FORMAT(1X,'TOO MANY GAUSSIANS IN VVO ATOMIC MINIMAL BASIS, >',I8)
 9035 FORMAT(1X,'NO VVOS CAN BE GENERATED FOR SVDCOPT,',
     *       1X,'OR FOR CCTYP WITH LOCALIZATION.')
 9040 FORMAT(/1X, 'THERE SHOULD BE A CLEAR DROPOFF IN THE EIGENVALUES',
     *           ' OF -B- BETWEEN THE'/
     *           ' VALENCE VIRTUAL AND EXTERNAL VIRTUAL ORBITALS,',
     *           ' PLEASE CHECK THIS.')
 9050 FORMAT(1X,'LAST FEW SVD EIGENVALUES KEPT AS VALENCE:')
 9060 FORMAT(1X,'FIRST FEW SVD EIGENVALUES REJECTED AS EXTERNAL:')
 9070 FORMAT(1X,5F15.10)
 9080 FORMAT(/1X,'THE QUASIATOMIC MINIMAL-BASIS-SET ORBITALS',
     *          ' (QUAMBO''S) HAVE BEEN PUNCHED.'/
     *       1X,'THE',I5,' QUAMBO''S ARE MOLECULAR ORBITALS WHICH'/
     *       1X,'   A) ARE MUTUALLY ORTHONORMAL,'/
     *       1X,'   B) SPAN THE',I5,' DIMENSIONAL OCCUPIED ORBITAL',
     *          ' SPACE OF THIS MOLECULE'/
     *       1X,'      EXACTLY, AND',I5,
     *          ' ADDITIONAL VALENCE VIRTUAL ORBITALS,'/
     *       1X,'   C) ARE AS CLOSE AS POSSIBLE TO THE FREE ATOM AO''S',
     *          ' (IN THE SENSE OF '/
     *       1X,'      MAXIMUM OVERLAP), WHILE SATISFYING "A" AND "B".')
 9082 FORMAT(/1X,'USING ZAPPED MCSCF AO FOCK MATRIX TO CANONICALIZE',
     *           ' THE VVOS.')
 9084 FORMAT(/1X,'USING THE RHF FOCK MATRIX TO CANONICALIZE',
     *           ' THE VVOS.')
 9086 FORMAT(/1X,'USING ROHF COMBINED FOCK MATRIX TO CANONICALIZE',
     *           ' THE VVOS.')
 9087 FORMAT(/1X,'USING GVB COMBINED FOCK MATRIX TO CANONICALIZE',
     *           ' THE VVOS.')
 9088 FORMAT(/1X,'NEED TO SELECT A DAF FILE FOR THE AO FOCK MATRIX.',
     *       /1X,'EITHER CREATE A NEW ELSEIF-THEN OR ADD TO ONE.')
 9090 FORMAT(/1X,'VIRTUAL VALENCE ORBITAL GENERATION IS FINISHED.'/
     *       1X,'THE VIRTUAL VALENCE ORBITALS ARE A QUANTIFICATION OF',
     *          ' THE FRONTIER ORBITAL'/
     *       1X,'SPACE (I.E. LUMO, LUMO+1, ...) WITHIN THE HARTREE-',
     *          'FOCK APPROXIMATION.'/
     *       1X,'THE',I5,' OCCUPIED AND',I5,' VIRTUAL ORBITALS',
     *          ' PRINTED BELOW ARE'/
     *       1X,'     FIRST, THE OCCUPIED CANONICAL ORBITALS (',
     *           I5,'-',I5,')'/
     *       1X,'    SECOND, THE VIRTUAL VALENCE ORBITALS    (',
     *           I5,'-',I5,')'/
     *       1X,'     THIRD, THE EXTERNAL VIRTUAL ORBITALS   (',
     *           I5,'-',I5,')'/
     *       1X,'NOTE THAT THE VIRTUAL SPACE IS ONLY SEMI-',
     *          'CANONICALIZED, MEANING THAT'/
     *       1X,'THE FOCK OPERATOR IS DIAGONAL ONLY WITHIN ITS',
     *          ' VALENCE AND EXTERNAL'/
     *       1X,'VIRTUAL SUBBLOCKS.  THE VIRTUAL "EIGENVALUES" BELOW',
     *          ' ARE THESE DIAGONALS,'/
     *       1X,'WHICH ARE EXPECTATION VALUES OF THE FOCK OPERATOR,',
     *          ' NOT ITS EIGENVALUES.')
 9095 FORMAT(///1X,'FOR LOCAL_LMOSVD RUNS,',
     *        /1X,'EXIT EARLY AND NEVER TRY',
     *         1X,'TO PUT ANY VVOS ORBITALS FOR ORBITAL SET.'///)
 9096 FORMAT(///1X,'FOR LOCAL_LMOSVD RUNS,',
     *         /1X,'NOCC=NUMBER OF MBS ORBITALS.',
     *        /1X,'HENCE, WE ACTUALLY NEVER PICK UP ANY NEW',
     *        /1X,'       VVOS ORBITALS INTO THE ORBITAL SET.'///)
 9100 FORMAT(/1X,'WARNING:  IVVOS=0 FOR VVOS RUN AFTER CHANGING',
     *        1X,'DAF FILE 15,ETC.',
     *       /1X,'THERE MIGHT BE SOMETHING WRONG WITH YOUR RUN.')
 9200 FORMAT(/1X,'CARRYING OUT ONLY THE SVD VIRTUAL SIDE ROTATION',
     *           ' (NO QUAMBOS)')
      END
C
C*MODULE VVOS    *DECK NVVOS_NUMCOR
C> @brief      Routine returns number of chemical core orbitals.
C>
C> @author     Mike Schmidt
C>             -1990s
C>
C> @details    Routine accounts for semi-cores so that only the truly
C>             valence orbitals are not included in the core.
C>
C> @date November 2, 2020 - George Schoendorff
C> - Added support for the relativistic AAMBS and for MCPs
C>
C> @param IATM1 indicates which cores to return when ITYPE equals 1.
C> @param ITYPE indicates whether to return cores for all atoms or one.
C> @date November 4, 2012-Aaron West
C> -Modified routine to return cores for all or individual atoms.
C>
C> @see NUMCOR
      INTEGER FUNCTION NVVOS_NUMCOR(IATM1,ITYPE)
      use mx_limits, only: mxatm
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C
      DIMENSION NCORES(103)
C
      LOGICAL GOPARR,DSKWRK,MASWRK,MCP
C
      COMMON /ECP2  / CLP(400),ZLP(400),NLP(400),KFIRST(MXATM,6),
     *                KLAST(MXATM,6),LMAX(MXATM),LPSKIP(MXATM),
     *                IZCORE(MXATM)
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NEX,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
C
      DATA NCORES /2*0,
     *             2*1,                  6*1,
     *             2*5,                  6*5,
     *             2*9,          10*9,   6*14,
     *             2*18,         10*18,  6*23,
     *             2*27,  14*27, 10*34,  6*39,
     *             2*43,  14*43,    50/
C
C             RETURN THE NUMBER OF CHEMICAL CORE ORBITALS
C        THIS DIFFERS FROM -NUMCOR- IN DEFINING SEMI-CORE ORBS AS CORE,
C        SO THAT ONLY TRULY VALENCE ORBITALS ARE NOT INCLUDED IN CORE.
C
      NCORE = 0
      IERR = 0
      IF(ITYPE.EQ.0) THEN
        DO 100 I=1,NAT
           IZ = INT(ZAN(I)) + IZCORE(I)
           IF(IZ.LE.0) GO TO 100
           IF((IZ.GT.103) .AND. MASWRK) WRITE(IW,9000)
           IF(IZ.GT.103) IERR=1
           NCORE = NCORE + NCORES(IZ) - IZCORE(I)/2
  100   CONTINUE
      ELSEIF(ITYPE.EQ.1) THEN
           IZ = INT(ZAN(IATM1)) + IZCORE(IATM1)
           IF(IZ.LE.0) GO TO 9998
           IF((IZ.GT.103) .AND. MASWRK) WRITE(IW,9000)
           IF(IZ.GT.103) IERR=1
           NCORE = NCORE + NCORES(IZ) - IZCORE(IATM1)/2
      ELSE
        IF(MASWRK) WRITE(IW,9001)
        IERR=1
      ENDIF
C
 9998 CONTINUE
      IF(IERR.GT.0) THEN
        IWX=IW
        IF(MASWRK) WRITE(IWX,9002)
        CALL FLSHBF(IWX)
        CALL ABRT
        STOP
      ENDIF
C
      NVVOS_NUMCOR = NCORE
      RETURN
 9000 FORMAT(/1X,'NVVOS_NUMCOR: Z.GT.103 IN SETTING CORE COUNT')
 9001 FORMAT(/1X,'NVVOS_NUMCOR: BAD CHOICE FOR ITYPE')
 9002 FORMAT(/1X,'NVVOS_NUMCOR: IZ IS LESS THAN 1.  WRONG.')
      END
C*MODULE VVOS    *DECK TFAO5D
C>
C>    @brief   transform large GTO MBS to spherical harmonics.
C>
C>    @author  Mike Schmidt
C>
C>    @details transform the minimal basis set side(s) of either
C>             the overlap within the MBS, or between the virtual
C>             orbitals and the MBS.  Nothing happens at all if
C>             there are no atoms with occupied d's.
C>
C> @date October 18, 2017 - George Schoendorff
C>  -Included support for f orbitals
C>
      SUBROUTINE TFAO5D(MODE,SMBS,OVLP,KTYPE,TF,WORK,VWRK,
     *                  L1,NSHLMBS,L0MBS,L1MBS,L2MBS)
      use constants, only: zero,half=>pt5,one
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C          only one of SMBS/VWRK or OVLP should be given,
C          depending on the execution MODE.
C          WORK is L2MBS for MODE=0, and L1*L1 for MODE=1
      DIMENSION SMBS(L2MBS),OVLP(L1,L1MBS),KTYPE(NSHLMBS),
     *          TF(L1MBS,L1MBS),WORK(*),VWRK(L1MBS)
      DIMENSION IDEG(4)
      DATA IDEG/1,3,6,10/
C
C       TRANSFORM THE MINIMAL ATOMIC BASIS TO SPHERICAL HARMONICS,
C       AND -L0MBS- WILL BE SET TO THE NUMBER OF SPHERICAL HARMONICS
C
C       IF MODE=0, THE METRIC IS THAT OF THE MINIMAL BASIS SPACE,
C                  AND NEEDS TO BE TRANSFORMED TO THE SPHERICAL
C                  HARMONIC PART OF THIS SPACE.
C                  ON ENTRY, -SMBS- IS L1MBS TRI, ON EXIT L0MBS TRI
C       IF MODE=1, THE METRIC IS BETWEEN THE MBS AND THE VIRTUAL
C                  ORBITALS OF SOME MOLECULAR BASIS, SO ONLY THE
C                  MBS SIDE SHOULD BE CHANGED TO SPHERICAL HARMONICS.
C                  ON ENTRY, -OVLP- IS L1 X L1MBS, ON EXIT L1 X L0MBS
C
      V34 = SQRT(3.0D+00/4.0D+00)
      V65 = SQRT(6.0D+00/5.0D+00)
      V340 = SQRT(3.0D+00/4.0D+01)
      V38 = SQRT(3.0D+00/8.0D+00)
      V920 = SQRT(9.0D+00/2.0D+01)
      V98 = SQRT(9.0D+00/8.0D+00)
      V58 = SQRT(5.0D+00/8.0D+00)
C
      CALL VCLR(TF,1,L1MBS*L1MBS)
      L1ROW=0
      L0COL=0
C
C        NOTE THAT OUR MBS IS A LARGE GTO EXPANSION W/O ANY L SHELLS
C
      DO ISH=1,NSHLMBS
         KTYP = KTYPE(ISH)
         IF(KTYP.EQ.1) THEN
            TF(L1ROW+1,L0COL+1) = ONE
            L0COL = L0COL + 1
         END IF
         IF(KTYP.EQ.2) THEN
            TF(L1ROW+1,L0COL+1) = ONE
            TF(L1ROW+2,L0COL+2) = ONE
            TF(L1ROW+3,L0COL+3) = ONE
            L0COL = L0COL + 3
         END IF
         IF(KTYP.EQ.3) THEN
            TF(L1ROW+1,L0COL+1) =  HALF
            TF(L1ROW+2,L0COL+1) =  HALF
            TF(L1ROW+3,L0COL+1) = -ONE
            TF(L1ROW+1,L0COL+2) =  V34
            TF(L1ROW+2,L0COL+2) = -V34
            TF(L1ROW+4,L0COL+3) = ONE
            TF(L1ROW+5,L0COL+4) = ONE
            TF(L1ROW+6,L0COL+5) = ONE
            L0COL = L0COL + 5
         END IF
CGS   The transformation for the f orbitals is performed in
CGS   DH8 symmetry as opposed to Oh. This is done to preserve
CGS   chemically relevant and intuitive values of Lz for
CGS   linear molecules.
         IF(KTYP.EQ.4) THEN
            TF(L1ROW+10,L0COL+1) = ONE
            TF(L1ROW+5,L0COL+2) =  V34
            TF(L1ROW+7,L0COL+2) = -V34
            TF(L1ROW+1,L0COL+3) = -V58
            TF(L1ROW+6,L0COL+3) =  V98
            TF(L1ROW+2,L0COL+4) = -V58
            TF(L1ROW+4,L0COL+4) =  V98
            TF(L1ROW+3,L0COL+5) =  ONE
            TF(L1ROW+5,L0COL+5) = -V920
            TF(L1ROW+7,L0COL+5) = -V920
            TF(L1ROW+1,L0COL+6) = -V38
            TF(L1ROW+6,L0COL+6) = -V340
            TF(L1ROW+8,L0COL+6) =  V65
            TF(L1ROW+2,L0COL+7) = -V38
            TF(L1ROW+4,L0COL+7) = -V340
            TF(L1ROW+9,L0COL+7) =  V65
            L0COL = L0COL + 7
         END IF
         IF(KTYP.GE.5) THEN
            WRITE(6,*) 'WE DO NOT HAVE G-BLOCK ELEMENTS'
            CALL ABRT
         END IF
         L1ROW = L1ROW + IDEG(KTYP)
      ENDDO
C
      L0MBS = L0COL
C
C        IF NEED BE, TRANSFORM MINIMAL BASIS SET TO SPHERICAL HARMONICS
C
      IF(L0MBS.EQ.L1MBS) RETURN
C
      IF(MODE.EQ.0) THEN
         CALL TFTRI(WORK,SMBS,TF,VWRK,L0MBS,L1MBS,L1MBS)
         CALL DCOPY(L2MBS,WORK,1,SMBS,1)
      END IF
C
      IF(MODE.EQ.1) THEN
         CALL DGEMM('N','N',L1,L0MBS,L1MBS,ONE,OVLP,L1,TF,L1MBS,
     *              ZERO,WORK,L1)
         CALL DCOPY(L1*L0MBS,WORK,1,OVLP,1)
      END IF
      RETURN
      END
C
C*MODULE VVOS    *DECK ETGTO
C>
C> @brief         Returns AAMBS exponents and contraction coefficients
C>
C> @author        Mike Schmidt
C>
C> @date November 2, 2020 - George Schoendorff
C>  -Included support the relativistic AAMBS up to radon and
C>   support for model core potentials
C>
C> @param NUCZ    is the true nuclear charge
C> @param DZ      is a flag to determine if the outermost Gaussian
C>                should be uncontracted
C> @param CSINP   are the input s coefficients of the AAMBS
C> @param CPINP   are the input p coefficients of the AAMBS
C> @param CDINP   are the input d coefficients of the AAMBS
C> @param CFINP   are the input f coefficients of the AAMBS
C> @param IERR1   is non-zero when there are errors
C> @param IERR2   is non-zero when there are errors
C> @param LOC     is an integer that serves as a pointer to a
C>                position in KLOC
C> @param INTYP   is an array of ITYP values that specify the
C>                angular moment
C> @param NANGM   defines the angular momentum of EACH KTYPE
C> @param NBFS    is the number of basis functions for each orbital
C> @param MINF    provides the index of KMIN
C> @param MAXF    provides the index of KMAX
C> @param NSHLS   is the number of shells
C> @param EX      is an array of Gaussian exponents
C> @param CS      is an array of s coefficients
C> @param CP      is an array of p coefficients
C> @param CD      is an array of d coefficients
C> @param CF      is an array of f coefficients
C> @param KSTART  is the location of the first exponent and the first
C>                contraction coefficient contained in a particular
C>                shell
C> @param KATOM   indicates which atom the shell is centered on
C> @param KTYPE   indicates the angular momentum of the shell
C> @param KNG     is the number of Gaussian functions in the shell
C> @param KLOC    is the location of the shell in the total AO basis
C> @param KMIN    is the starting index of the shell
C> @param KMAX    is the ending index of the shell
C> @param NGAUSS  is the number of Gaussian functions in a shell
C> @param NSHELL  is the number of shells
C> @param MXGTOT  is the maximum number of Gaussian functions permitted
C> @param MXSH    is the maximum number of shells permitted
C> @param IAT     is an index specifying the current atome
C>
      SUBROUTINE ETGTO(NUCZ,DZ,CSINP,CPINP,CDINP,CFINP,IERR1,IERR2,LOC,
     *                 INTYP,NANGM,NBFS,MINF,MAXF,NSHLS,
     *                 EX,CS,CP,CD,CF,KSTART,KATOM,KTYPE,
     *                 KNG,KLOC,KMIN,KMAX,NGAUSS,NSHELL,MXGTOT,MXSH,IAT)
      use mx_limits, only: mxatm, lenilist
      use constants, only: zero,one,two,pt5,pt75,tm6,tm8,tm10,pi32
      use comm_vvopar, only: ivvtyp
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      LOGICAL DONE,DZ,GOPARR,DSKWRK,MASWRK
      LOGICAL MCP
C
      PARAMETER (LENEEX=128)
      PARAMETER (LENCOEF=590)
C
      DIMENSION CSINP(MXGTOT),CPINP(MXGTOT),CDINP(MXGTOT),
     *          CFINP(MXGTOT),
     *          INTYP(*),NANGM(4),NBFS(4),MINF(4),MAXF(4),
     *          NSHLS(*),
     *          EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),
     *          CF(MXGTOT),
     *          KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),
     *          KLOC(MXSH),KMIN(MXSH),KMAX(MXSH)
C         NEXT DIMENSIONS ARE FOR THE FIRST FIVE ROWS OF THE
C         PERIODIC TABLE, AND MAYBE FLOATING OUTER S,P IN
C         THE MAIN GROUP ELEMENTS.  IT ALSO HAS TO DO WITH THE
C         USE OF SCHMIDT/RUEDENBERG AND HUZINAGA/KLOBUKOWSKI
C         LARGE ALL ELECTRON BASIS SETS, AND THEIR PRIMITIVE COUNTS.
C         FOR EXAMPLE, IODINE WTBS HAS THESE PRIMITIVE SIZES:
C         1S=28, 2S=28, 2P=23, 3S=28, 3P=23, 4S=28, 3D=17, 4P=23,
C         5S=28, 4D=17, 5P=23 =266, FLOAT ONE S,P WOULD BE 268 GTOS.
      DIMENSION EEX(LENEEX),COEF(LENCOEF)
      DIMENSION IELIST(LENILIST),ICLIST(LENILIST),
     *          ITLIST(LENILIST),IGLIST(LENILIST)
C         4 HERE AND ABOVE IS BECAUSE PERIODIC TABLE HAS S,P,D,F OCCUP.
      CHARACTER*1 TYPE(4)
C
C
      COMMON /ECP2  / CLP(400),ZLP(400),NLP(400),KFIRST(MXATM,6),
     *                KLAST(MXATM,6),LMAX(MXATM),LPSKIP(MXATM),
     *                IZCORE(MXATM)
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /INTOPT/ ISCHWZ,IECP,NECP,IEFLD
      COMMON /OUTPUT/ NPRINT,ITOL,ICUT,NORMF,NORMP,NOPK
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
C
      DATA TYPE/'S','P','D','F'/
      DATA PT1875/1.875D+00/
C
      MCP=.FALSE.
      IF(IZCORE(IAT).GT.0.AND.(IECP.EQ.5.OR.IECP.EQ.6)) MCP=.TRUE.
C
      CALL VCLR(EEX,1,LENEEX)
      CALL VCLR(COEF,1,LENCOEF)
C
      IF(NUCZ.GT.86) THEN
         CALL BERROR(3)
         STOP
      ENDIF
C
C     ----- Get the AAMBS for the current atom -----
C
      CALL AAMBS_ATOM(EEX,COEF,NUCZ,MCP,IVVTYP)
C
C     ----- LOOP OVER EACH SHELL -----
C
      IF(MCP) THEN
         CALL ETSHMCP(NUCZ,DZ,NPASS,
     *             LENILIST,IELIST,ICLIST,ITLIST,IGLIST)
      ELSE
         CALL ETSHELL(NUCZ,DZ,NPASS,
     *             LENILIST,IELIST,ICLIST,ITLIST,IGLIST)
      ENDIF
C
      DONE  = .FALSE.
      IERR3 = 0
      IPASS = 0
  210 CONTINUE
      IPASS = IPASS+1
      IF(IPASS.GT.NPASS) DONE=.TRUE.
      IF(DONE. AND. IERR3.NE.0) CALL ABRT
      IF(DONE) RETURN
C
      NGE    = IELIST(IPASS)
      NGC    = ICLIST(IPASS)
      ITYP   = ITLIST(IPASS)
      IGAUSS = IGLIST(IPASS)
C
C     ----- DEFINE THE CURRENT SHELL -----
C
      NSHELL = NSHELL+1
      IF(NSHELL.GT.MXSH) THEN
         IERR1=1
         RETURN
      END IF
      NSHLS(NAT)     = NSHLS(NAT)+1
      KMIN(NSHELL)   = MINF(ITYP)
      KMAX(NSHELL)   = MAXF(ITYP)
      KSTART(NSHELL) = NGAUSS+1
      KATOM(NSHELL)  = NAT
      KTYPE(NSHELL)  = NANGM(ITYP)
      INTYP(NSHELL)  = ITYP
      KNG(NSHELL)    = IGAUSS
      KLOC(NSHELL)   = LOC+1
      NGAUSS         = NGAUSS+IGAUSS
      IF(NGAUSS.GT.MXGTOT) THEN
         IERR2=1
         RETURN
      END IF
      LOC   = LOC  +NBFS(ITYP)
      K1 = KSTART(NSHELL)
      K2 = K1+KNG(NSHELL)-1
      DO 440 I = 1,IGAUSS
         K = K1+I-1
         EX(K) = EEX(NGE+I)
         IF(ITYP.EQ.1) CSINP(K) = COEF(NGC+I)
         IF(ITYP.EQ.2) CPINP(K) = COEF(NGC+I)
         IF(ITYP.EQ.3) CDINP(K) = COEF(NGC+I)
         IF(ITYP.EQ.4) CFINP(K) = COEF(NGC+I)
         CS(K) = ZERO
         CP(K) = ZERO
         CD(K) = ZERO
         CF(K) = ZERO
         IF(ITYP.EQ.1) CS(K) = CSINP(K)
         IF(ITYP.EQ.2) CP(K) = CPINP(K)
         IF(ITYP.EQ.3) CD(K) = CDINP(K)
         IF(ITYP.EQ.4) CF(K) = CFINP(K)
  440 CONTINUE
C
C     ----- ALWAYS UNNORMALIZE PRIMITIVES -----
C
      DO 460 K = K1,K2
         EE = EX(K)+EX(K)
         FACS = PI32/(EE*SQRT(EE))
         FACP = PT5*FACS/EE
         FACD = PT75*FACS/(EE*EE)
         FACF = PT1875*FACS/(EE*EE*EE)
         IF(ITYP.EQ.1) CS(K) = CS(K)/SQRT(FACS)
         IF(ITYP.EQ.2) CP(K) = CP(K)/SQRT(FACP)
         IF(ITYP.EQ.3) CD(K) = CD(K)/SQRT(FACD)
         IF(ITYP.EQ.4) CF(K) = CF(K)/SQRT(FACF)
  460 CONTINUE
C
C     ----- IF(NORMF.EQ.0) NORMALIZE BASIS FUNCTIONS -----
C
      IF (NORMF .EQ. 1) GO TO 210
      FACS = ZERO
      FACP = ZERO
      FACD = ZERO
      FACF = ZERO
      DO 510 IG = K1,K2
         DO 500 JG = K1,IG
            EE = EX(IG)+EX(JG)
            FAC = EE*SQRT(EE)
            DUMS =      CS(IG)*CS(JG)/FAC
            DUMP =  PT5*CP(IG)*CP(JG)/(EE*FAC)
            DUMD = PT75*CD(IG)*CD(JG)/(EE*EE*FAC)
            DUMF = PT1875*CF(IG)*CF(JG)/(EE*EE*EE*FAC)
            IF (IG .EQ. JG) GO TO 480
               DUMS = DUMS+DUMS
               DUMP = DUMP+DUMP
               DUMD = DUMD+DUMD
               DUMF = DUMF+DUMF
  480       CONTINUE
            FACS = FACS+DUMS
            FACP = FACP+DUMP
            FACD = FACD+DUMD
            FACF = FACF+DUMF
  500    CONTINUE
  510 CONTINUE
C
      FAC=ZERO
      IF(ITYP.EQ.1 .AND. FACS.GT.TM10) FAC=ONE/SQRT(FACS*PI32)
      IF(ITYP.EQ.2 .AND. FACP.GT.TM10) FAC=ONE/SQRT(FACP*PI32)
      IF(ITYP.EQ.3 .AND. FACD.GT.TM10) FAC=ONE/SQRT(FACD*PI32)
      IF(ITYP.EQ.4 .AND. FACF.GT.TM10) FAC=ONE/SQRT(FACF*PI32)
C
C                                 VERIFY NORMALIZATION
CGS
C     Some orbitals have normalization constants that deviate
C     from 1.0 by a bit more than TM6. Most of these issues
C     were solved by reoptimization of the coefficients with 
C     the fixed precision exponents. However, the 3p orbital for
C     ruthenium remains just outside the tolerance, i.e.,
C     0.10000013E+01, so we adjust the tolerance for ruthenium
C     only.
CGS
      TOLERANCE = TM6
      IF(NUCZ.EQ.44) TOLERANCE = TWO*TM6
      TNORM = ABS(ONE-FAC)
      IF(TNORM.LT.TOLERANCE) GO TO 520
C
      IF (MASWRK) THEN
         WRITE(IW,9000) FAC,NUCZ,IPASS,TYPE(ITYP)
         WRITE(IW,9010)  (EEX(NGE+K),K=1,IGAUSS)
         WRITE(IW,9020) (COEF(NGC+K),K=1,IGAUSS)
      ENDIF
      IERR3=IERR3+1
C
  520 CONTINUE
      DO 550 IG = K1,K2
         IF(ITYP.EQ.1) CS(IG) = FAC*CS(IG)
         IF(ITYP.EQ.2) CP(IG) = FAC*CP(IG)
         IF(ITYP.EQ.3) CD(IG) = FAC*CD(IG)
         IF(ITYP.EQ.4) CF(IG) = FAC*CF(IG)
         IF(ITYP.EQ.1) CSINP(IG) = FAC*CSINP(IG)
         IF(ITYP.EQ.2) CPINP(IG) = FAC*CPINP(IG)
         IF(ITYP.EQ.3) CDINP(IG) = FAC*CDINP(IG)
         IF(ITYP.EQ.4) CFINP(IG) = FAC*CFINP(IG)
  550 CONTINUE
      GO TO 210
C
 9000 FORMAT(/1X,'ERROR!!! NORMALIZATION FACTOR=',E16.8/
     *     1X,'FOR ATOM Z=',I3,' SHELL NO.',I3,' TYPE=',A1/
     *     1X,'CHECK BUILT IN ET-GTO EXPONENTS AND CONT. COEFS')
 9010 FORMAT(5F15.5)
 9020 FORMAT(5F15.8)
      END
C
C*MODULE VVOS    *DECK ETSHELL
C>
C>    @brief          Defines the AAMBS shell parameters
C>
C>    @author         Mike Schmidt
C>
C>    @date           11/21/22 - George Schoendorff
C>                     - IVVTYP=1 now gets the contraction lengths from
C>                       AAMBS_CONTRACTION_LENGTH so that the 
C>                       contraction lengths are uniquely defined only
C>                       in AAMBS_LIMITS in mod_vvos.src
C>
C>    @date           5/27/24 - George Schoendorff
C>                     - Condensed definitions of ITLIST
C>
C>    @param NUCZ     is the true nuclear charge of the atom
C>    @param DZ       is a flag to determine if we float the
C>                    outer exponents
C>    @param NPASS    is the number of AAMBS orbitals for the atom
C>    @param LENILIST is the maximum number of shells
C>    @param IELIST   is an array specifying the starting index of
C>                    the exponents for each orbital type
C>    @param ICLIST   is an array specifying the starting index of
C>                    the coefficients for each orbital
C>    @param ITLIST   is an array specifying the orbital types
C>    @param IGLIST   is an array specifying the number of Gaussian
C>                    exponents per orbital
C>
      SUBROUTINE ETSHELL(NUCZ,DZ,NPASS,
     *                   LENILIST,IELIST,ICLIST,ITLIST,IGLIST)
      use comm_vvopar, only: ivvtyp
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION IELIST(LENILIST),ICLIST(LENILIST),
     *          ITLIST(LENILIST),IGLIST(LENILIST)
      LOGICAL DZ
      INTEGER AAMBS_CONTRACTION_LENGTH
C
      CHARACTER*4 FLOAT
      LOGICAL GOPARR,DSKWRK,MASWRK
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON /INTOPT/ ISCHWZ,IECP,NECP,IEFLD
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
C
C         NEED TO THINK ABOUT WHICH EXPONENT TO FLOAT,
C         AND TO ADD CODE IN WTBSFOUR, WTBSFIVE TO DO IT.
C
C        FILL LISTS OF EXPONENT POINTER, COEFFICIENT POINTER,
C                      TYPE OF AO, AND NUMBER OF GAUSSIANS
C
      IF(LENILIST.NE.19) THEN
         IF(MASWRK) WRITE(IW,*) 'CHANGE LENILIST SOMEWHERE'
         WRITE(IW,*) 'LENILIST=',LENILIST
         CALL FLSHBF(IW)
         CALL ABRT
         STOP
      ENDIF
C------------
      NPASS=0
      CALL VICLR(ITLIST,1,LENILIST)
      NS=0
      NP=0
      ND=0
      NF=0
      NS =  AAMBS_CONTRACTION_LENGTH(NUCZ,1,0,IVVTYP)
      IF((NUCZ.GE.3.AND.IVVTYP.EQ.1).OR.(NUCZ.GE.5.AND.IVVTYP.EQ.0))
     *   NP =  AAMBS_CONTRACTION_LENGTH(NUCZ,2,1,IVVTYP)
      IF(IVVTYP.EQ.0.AND.NUCZ.GE.21)
     *   ND = AAMBS_CONTRACTION_LENGTH(NUCZ,3,2,IVVTYP)
      IF(IVVTYP.EQ.1.AND.NUCZ.GE.19) 
     *   ND = AAMBS_CONTRACTION_LENGTH(NUCZ,3,2,IVVTYP)
      IF(NUCZ.GE.57) NF = AAMBS_CONTRACTION_LENGTH(NUCZ,4,3,IVVTYP)
      IF(NUCZ.GT. 2) GO TO 020
         NPASS=1
         ITLIST(1)=1
         IF(DZ) FLOAT='S   '
         GO TO 800
  020 IF(NUCZ.GT. 4) GO TO 030
         ITLIST(1)=1
         ITLIST(2)=1
         IF(IVVTYP.EQ.0) THEN
            NPASS=2
            IF(DZ) FLOAT='S   '
            GO TO 800
         ELSEIF(IVVTYP.EQ.1) THEN
            NPASS=3
            ITLIST(3)=2
            IF(DZ) FLOAT='SP  '
            GO TO 800
         ELSE
            GO TO 720
         ENDIF
  030 IF(NUCZ.GT.10) GO TO 040
         ITLIST(1)=1
         ITLIST(2)=1
         ITLIST(3)=2
         IF(IVVTYP.EQ.0) THEN
            NPASS=3
            IF(DZ) FLOAT='SP  '
            GO TO 800
         ELSEIF(IVVTYP.EQ.1) THEN
            NPASS=3
            IF(DZ) FLOAT='SP  '
            GO TO 800
         ELSE
            GO TO 720
         ENDIF
  040 IF(NUCZ.GT.12) GO TO 050
         ITLIST(1)=1
         ITLIST(2)=1
         ITLIST(3)=2
         ITLIST(4)=1
         IF(IVVTYP.EQ.0) THEN
            NPASS=4
            IF(DZ) FLOAT='S   '
            GO TO 800
         ELSEIF(IVVTYP.EQ.1) THEN
            NPASS=5
            ITLIST(5)=2
            IF(DZ) FLOAT='SP  '
            GO TO 800
         ELSE
            GO TO 720
         ENDIF
  050 IF(NUCZ.GT.18) GO TO 060
         NPASS=5
         ITLIST(1)=1
         ITLIST(2)=1
         ITLIST(3)=2
         ITLIST(4)=1
         ITLIST(5)=2
         IF(IVVTYP.EQ.0.OR.IVVTYP.EQ.1) THEN
            IF(DZ) FLOAT='SP  '
            GO TO 800
         ELSE
            GO TO 720
         ENDIF
C                                    K-CA
  060 IF(NUCZ.GT.20) GO TO 070
         ITLIST(1)=1
         ITLIST(2)=1
         ITLIST(3)=2
         ITLIST(4)=1
         ITLIST(5)=2
         ITLIST(6)=1
         IF(IVVTYP.EQ.0) THEN
            NPASS=6
            GO TO 800
         ELSEIF(IVVTYP.EQ.1) THEN
            NPASS=8
            ITLIST(7)=3
            ITLIST(8)=2
            GO TO 800
         ELSE
            GO TO 720
         ENDIF
C                                    Sc-Ni
  070 IF(NUCZ.GT.30) GO TO 080
         ITLIST(1)=1
         ITLIST(2)=1
         ITLIST(3)=2
         ITLIST(4)=1
         ITLIST(5)=2
         ITLIST(6)=1
         ITLIST(7)=3
         IF(IVVTYP.EQ.0) THEN
            IF(NUCZ.EQ.29 .OR. NUCZ.EQ.30) ND=14
            NPASS=7
            GO TO 800
         ELSEIF(IVVTYP.EQ.1) THEN
            NPASS=8
            ITLIST(8)=2
            GO TO 800
         ELSE
            GO TO 720
         ENDIF
C not sure what to float for the metals
C                                    GA-KR
  080 IF(NUCZ.GT.36) GO TO 090
         IF(IVVTYP.EQ.0.OR.IVVTYP.EQ.1) THEN
            NPASS=8
            ITLIST(1)=1
            ITLIST(2)=1
            ITLIST(3)=2
            ITLIST(4)=1
            ITLIST(5)=2
            ITLIST(6)=3
            ITLIST(7)=1
            ITLIST(8)=2
            IF(DZ) FLOAT='SP  '
            GO TO 800
         ELSE
            GO TO 720
         ENDIF
C                                    RB-SR
  090 IF(NUCZ.GT.38) GO TO 100
         ITLIST(1)=1
         ITLIST(2)=1
         ITLIST(3)=2
         ITLIST(4)=1
         ITLIST(5)=2
         ITLIST(6)=1
         ITLIST(7)=3
         ITLIST(8)=2
         ITLIST(9)=1
         IF(IVVTYP.EQ.0) THEN
            NPASS=9
            GO TO 800
         ELSEIF(IVVTYP.EQ.1) THEN
            NPASS=11
            ITLIST(10)=3
            ITLIST(11)=2
            GO TO 800
         ELSE
            GO TO 720
         ENDIF
C                                    Y-Cd
  100 IF(NUCZ.GT.48) GO TO 110
         ITLIST( 1)=1
         ITLIST( 2)=1
         ITLIST( 3)=2
         ITLIST( 4)=1
         ITLIST( 5)=2
         ITLIST( 6)=1
         ITLIST( 7)=3
         ITLIST( 8)=2
         ITLIST( 9)=1
         ITLIST(10)=3
         IF(IVVTYP.EQ.0) THEN
            NPASS=10
            GO TO 800
         ELSEIF(IVVTYP.EQ.1) THEN
            NPASS=11
            ITLIST(11)=2
            GO TO 800
         ELSE
            GO TO 720
         ENDIF

C not sure what to float for metals
C                                    IN-XE
  110 IF(NUCZ.GT.54) GO TO 120
         IF(IVVTYP.EQ.0.OR.IVVTYP.EQ.1) THEN
            NPASS=11
            ITLIST( 1)=1
            ITLIST( 2)=1
            ITLIST( 3)=2
            ITLIST( 4)=1
            ITLIST( 5)=2
            ITLIST( 6)=3
            ITLIST( 7)=1
            ITLIST( 8)=2
            ITLIST( 9)=3
            ITLIST(10)=1
            ITLIST(11)=2
            IF(DZ) FLOAT='SP  '
            GO TO 800
         ELSE
            GO TO 720
         ENDIF
C
C                                 CS-BA
  120 IF(NUCZ.GT.56) GO TO 130
        NPASS=14
        ITLIST( 1)=1
        ITLIST( 2)=1
        ITLIST( 3)=2
        ITLIST( 4)=1
        ITLIST( 5)=2
        ITLIST( 6)=3
        ITLIST( 7)=1
        ITLIST( 8)=2
        ITLIST( 9)=3
        ITLIST(10)=1
        ITLIST(11)=2
        ITLIST(12)=1
        ITLIST(13)=3
        ITLIST(14)=2
        IF(DZ) FLOAT='SP  '
        GO TO 800
C
C                                 LA-YB
  130 IF(NUCZ.GT.70) GO TO 140
        NPASS = 15
        ITLIST( 1)=1
        ITLIST( 2)=1
        ITLIST( 3)=2
        ITLIST( 4)=1
        ITLIST( 5)=2
        ITLIST( 6)=3
        ITLIST( 7)=1
        ITLIST( 8)=2
        ITLIST( 9)=3
        ITLIST(10)=1
        ITLIST(11)=2
        ITLIST(12)=1
        ITLIST(13)=4
        ITLIST(14)=3
        ITLIST(15)=2
        GO TO 800
C
C                                 LU-HG
  140 IF(NUCZ.GT.80) GO TO 150
        NPASS = 15
        ITLIST( 1)=1
        ITLIST( 2)=1
        ITLIST( 3)=2
        ITLIST( 4)=1
        ITLIST( 5)=2
        ITLIST( 6)=3
        ITLIST( 7)=1
        ITLIST( 8)=2
        ITLIST( 9)=3
        ITLIST(10)=4
        ITLIST(11)=1
        ITLIST(12)=2
        ITLIST(13)=1
        ITLIST(14)=3
        ITLIST(15)=2
        GO TO 800
C
C                                 TL-RN
  150 IF(NUCZ.GT.86) GO TO 720
        NPASS=15
        ITLIST( 1)=1
        ITLIST( 2)=1
        ITLIST( 3)=2
        ITLIST( 4)=1
        ITLIST( 5)=2
        ITLIST( 6)=3
        ITLIST( 7)=1
        ITLIST( 8)=2
        ITLIST( 9)=3
        ITLIST(10)=1
        ITLIST(11)=2
        ITLIST(12)=4
        ITLIST(13)=3
        ITLIST(14)=1
        ITLIST(15)=2
        IF(DZ) FLOAT='SP  '
        GO TO 800
C
  720 CONTINUE
      IF(MASWRK) WRITE(IW,*) 'IMPOSSIBLE MISTAKE IN ETSHELL'
      CALL ABRT
      STOP
C
C----------------
  800 CONTINUE
      IF(NPASS.EQ.0) THEN
         IF(MASWRK) WRITE(IW,*) 'ETSHELL ERROR:  NPASS CANNOT BE 0.'
         CALL FLSHBF(IW)
         CALL ABRT
         STOP
      ENDIF
      IF(NPASS+2.GT.LENILIST) THEN
         IF(MASWRK) THEN
            WRITE(IW,*) 'ETSHELL ERROR:'
            WRITE(IW,*) 'NPASS+2=',NPASS+2
            WRITE(IW,*) 'LENILIST=',LENILIST
            WRITE(IW,*) 'INCREASE LENILIST FOR MEMORY.  TWO PLACES.'
         ENDIF
         CALL FLSHBF(IW)
         CALL ABRT
         STOP
      ENDIF
CGS MCP potential error out
      NERR=0
      DO I=1,NPASS
         IF(ITLIST(I).EQ.0) NERR=NERR+1
      ENDDO
      IF(NERR.NE.0) THEN
         IF(MASWRK) THEN
            WRITE(IW,*) 'ETSHELL ERROR:'
            WRITE(IW,*) 'PROBLEM WITH ITLIST BEFORE NON-FLOAT.'
         ENDIF
         CALL FLSHBF(IW)
         CALL ABRT
         STOP
      ENDIF
C----------------
      IELIST(1)=0
      ICLIST(1)=0
      IGLIST(1)=NS
      DO I=2,NPASS
         IF(ITLIST(I).EQ.1) THEN
            IELIST(I) = 0
            IGLIST(I) = NS
         END IF
         IF(ITLIST(I).EQ.2) THEN
            IELIST(I) = NS
            IGLIST(I) = NP
         END IF
         IF(ITLIST(I).EQ.3) THEN
            IELIST(I) = NS+NP
            IGLIST(I) = ND
         END IF
         IF(ITLIST(I).EQ.4) THEN
            IELIST(I) = NS+NP+ND
            IGLIST(I) = NF
         END IF
         ICLIST(I) = ICLIST(I-1) + IGLIST(I-1)
      ENDDO
      IF(FLOAT.EQ.'S   ') THEN
         ITLIST(NPASS+1) = 1
         IGLIST(NPASS+1) = 1
         IELIST(NPASS+1) = 0
         ICLIST(NPASS+1) = ICLIST(NPASS) + IGLIST(NPASS)
         NPASS=NPASS+1
      END IF
      IF(FLOAT.EQ.'SP  ') THEN
         ITLIST(NPASS+1) = 1
         ITLIST(NPASS+2) = 2
         IGLIST(NPASS+1) = 1
         IGLIST(NPASS+2) = 1
         IELIST(NPASS+1) = 0
         IELIST(NPASS+2) = NS
         ICLIST(NPASS+1) = ICLIST(NPASS)   + IGLIST(NPASS)
         ICLIST(NPASS+2) = ICLIST(NPASS+1) + IGLIST(NPASS+1)
         NPASS=NPASS+2
      END IF
C
C-    WRITE(6,*) 'EXIT FROM ETSHELL ROUTINE, NPASS=',NPASS
C-    WRITE(6,899) 'IELIST',(IELIST(I),I=1,NPASS)
C-    WRITE(6,899) 'ICLIST',(ICLIST(I),I=1,NPASS)
C-    WRITE(6,899) 'ITLIST',(ITLIST(I),I=1,NPASS)
C-    WRITE(6,899) 'IGLIST',(IGLIST(I),I=1,NPASS)
C-899 FORMAT(A,'=',5I5/7X,10I5)
C
      RETURN
      END
C*MODULE VVOS    *DECK BNDANX
C> @brief      Routine drives formation and orientiation of vvos.
C>
C> @author     Mike Schmidt
C>             -2004
C>
C> @details   Routine drives calculations from energy to orientation.
C>
C> @date November 4, 2012-Aaron West
C> -Modified nvvos_numcor integer calls.
C>
      SUBROUTINE BNDANX
      use mx_limits, only: mxao,mxatm,mxrt
      use constants, only: one,two
      use comm_vvopar, only: nvvos,ivvtyp,bndden
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION LOCMO
      DIMENSION FANT(8),LFANT(8),GANT(27),LGANT(8),LOCMO(3),
     *          NRNFG(10),NPFLG(10)
      LOGICAL GOPARR,DSKWRK,MASWRK,OK,ABELSYM,ABEL,ABELPT,
     *        ANALYS,FDIRCT,QCORR,DIRSCF,FDIFF,DIRTRF,SOCI,
     *        SAFLG
      CHARACTER*5 CICODE
C
      COMMON /DETPAR/ ICLBBR,ANALYS,SAFLG,IPRTSA
      COMMON /DETWFN/ WSTATE(MXRT),SPINS(MXRT),CRIT,PRTTOL,S,SZ,
     *                GRPDET,STSYM,GLIST,DWPARM,
     *                NFLGDM(MXRT),IWTS(MXRT),NCORSV,NCOR,NACT,NORB,
     *                NADET,NBDET,K,KST,IROOT,IPURES,MAXW1,NITER,
     *                MAXP,NCI,IGPDET,KSTSYM,NFTGCI,IDWEIGH,
     *                fstate(mxrt),ifts(mxrt)
      COMMON /DESOCI/ OSPIN(MXRT),NEXT,MOSET,MAXPSO,KSO,NSOCI
      COMMON /FCCWFN/ NSPACE,MSTA(51),MNUM(51),MINI(51),MAXI(51),
     *                IAMI(51),IAMA(51),IBMI(51),IBMA(51),IDIM(51),
     *                LBST(51),NREF0,FDIRCT,QCORR,C0SQ
      COMMON /FMCOM / XX(1)
      COMMON /IJPAIR/ IA(MXAO)
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NEX,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON /OPTSCF/ DIRSCF,FDIFF
      COMMON /PRPOPT/ ETOLLZ,ILOCAL,IAHARD
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /RUNOPT/ RUNTYP,EXETYP,NEVALS,NGLEVL,NHLEVL
      COMMON /SYMMOL/ GROUP,COMPLEX,IGROUP,NAXIS,ILABMO,ABELSYM
      COMMON /WFNOPT/ SCFTYP,VBTYP,DFTYPE,TDDFTYP,CITYP,CCTYP,
     *                MPLEVL,MPCTYP
C
C      ----- INITIALIZE $BNDANA INPUT GROUP -----
C
      PARAMETER (NNAM=1)
      DIMENSION QNAM(NNAM),KQNAM(NNAM)
      DATA BNDGRP/8HBNDANA  /
      DATA   QNAM/8HDENS    /
      DATA  KQNAM/5/
C
      DATA SCF,CISD,CISDT,CISDTQ,ORMAS
     *   /8HSCF     ,8HCISD    ,8HCISDT   ,8HCISDTQ  ,8HORMAS   /
      DATA LOCMO/8HBOYS    ,8HRUEDNBRG,8HPOP     /
      DATA CHECK/8HCHECK   /
C
      DATA FANT/8HC1      ,8HCI      ,8HCS      ,8HC2      ,
     *          8HD2      ,8HC2V     ,8HC2H     ,8HD2H     /
      DATA LFANT/1,1,1,1,2,2,2,3/
      DATA GANT/8HA       ,8HAG      ,8HAU      ,8HA'      ,
     *          8HA"      ,8HA       ,8HB       ,8HA       ,
     *          8HB1      ,8HB2      ,8HB3      ,8HA1      ,
     *          8HA2      ,8HB1      ,8HB2      ,8HAG      ,
     *          8HBG      ,8HBU      ,8HAU      ,8HAG      ,
     *          8HB1G     ,8HB2G     ,8HB3G     ,8HAU      ,
     *          8HB1U     ,8HB2U     ,8HB3U     /
      DATA LGANT/0,1,3,5,7,11,15,19/
C
      DATA NRNFG /1,1,1,0,1,0,0,0,0,0/
      DATA NPFLG /0,0,0,0,0,0,0,0,0,0/
C
      DMATRX=CISD
      JRET = 0
      CALL NAMEIO(IR,JRET,BNDGRP,NNAM,QNAM,KQNAM,
     *            DMATRX,  0,0,0,
     *   0,0,0,0,0,    0,0,0,0,0,   0,0,0,0,0,   0,0,0,0,0,
     *   0,0,0,0,0,    0,0,0,0,0,   0,0,0,0,0,   0,0,0,0,0,
     *   0,0,0,0,0,    0,0,0,0,0,   0,0,0,0,0,   0,0,0,0,0,
     *   0,0,0,0,0, 0)
      IF(JRET.EQ.2) THEN
         IF(MASWRK) WRITE(IW,9000)
         CALL ABRT
      END IF
C
      OK = .FALSE.
      IF(DMATRX.EQ.SCF)    OK=.TRUE.
      IF(DMATRX.EQ.CISD)   OK=.TRUE.
      IF(DMATRX.EQ.CISDT)  OK=.TRUE.
      IF(DMATRX.EQ.CISDTQ) OK=.TRUE.
      IF(DMATRX.EQ.ORMAS)  THEN
         IF(CITYP.EQ.ORMAS) THEN
            OK=.TRUE.
         ELSE
            IF(MASWRK) WRITE(IW,9010)
         END IF
      END IF
C
      IF(OK) THEN
         BNDDEN = DMATRX
      ELSE
         IF(MASWRK) WRITE(IW,9020) DMATRX,SCF,CISD,CISDT,CISDTQ,ORMAS
         CALL ABRT
      END IF
C
C        --- PREPARE A CLOSED SHELL WAVEFUNCTION, AND ITS VVOS ---
C        IN CASE CITYP=ORMAS IS USED, THE NEXT CALL ALSO GENERATES
C        THE DENSITY MATRIX TO BE USED.  THIS CUSTOM ORMAS SHOULD
C        PROBABLY NOT BE MENTIONED IN THE DOCUMENTATION.
C
      CALL ENERGX
C
      MCORE = NVVOS_NUMCOR(0,0)
      IF(MASWRK) WRITE(IW,9030) LOCMO(ILOCAL),SCFTYP,
     *                          NA-MCORE,NVVOS,DMATRX
C
C       --- GET A DENSITY FOR THESE SEMICANONICAL VALENCE ORBITALS ---
C
      IF(DMATRX.EQ.SCF) THEN
         IF(MASWRK) WRITE(IW,9040)
      END IF
C
      IF(DMATRX.EQ.ORMAS) THEN
         IF(MASWRK) WRITE(IW,9050)
      END IF
C
C           USUAL OPTION IS A HARDWIRED GROUND STATE CI-SD WITH FSOCI,
C           BUT ORMAS CAN GENERATE HIGHER CI-SDT OR CI-SDTQ
C           THE FSOCI CODE HAS GLITCHES FOR CASES WITH JUST TWO E-,
C           BUT ORMAS EATS THOSE TRIVIAL RUNS UP.
C
      IF(DMATRX.EQ.CISD    .OR.
     *   DMATRX.EQ.CISDT   .OR.
     *   DMATRX.EQ.CISDTQ) THEN
C
         IF(DMATRX.EQ.CISD)   IEXCIT=2
         IF(DMATRX.EQ.CISDT)  IEXCIT=3
         IF(DMATRX.EQ.CISDTQ) IEXCIT=4
         SOCI = IEXCIT.EQ.2  .AND.  (NA+NB-2*MCORE).GT.2
     *                       .AND.  NPROC.EQ.1
C
C           parallel runs appear to carry out the ORMAS calculation
C           for CI-SD correctly, but the localization picks up the
C           incorrect integrals, and so doesn't localize right.
C           In 2016, stop parallel runs here, until that is addressed.
C
         IF(NPROC.GT.1) THEN
            IF(MASWRK) WRITE(IW,*) 'BONDING ANALYSIS MUST RUN SERIALLY'
            CALL ABRT
         END IF
C
C           SET UP TYPICAL DETERMINANT SPECIFICATIONS
C
         SZ     = (MUL-1)/TWO
         IPURES = 1
C
         NCOR   = MCORE
         NCORSV = MCORE
         IF(SOCI) THEN
            CICODE='FSOCI'
            NACT   = NA - MCORE
            NORB   = MCORE + NACT
         ELSE
            CICODE='ORMAS'
            NACT   = NA - MCORE + NVVOS
            NORB   = MCORE + NACT
         END IF
C
         NHIGH = INT(SZ+SZ+0.0001D+00)
         NELS  = NA + NB - 2*MCORE
         NBDET = (NELS-NHIGH)/2
         NADET = NBDET+NHIGH
C
         KSTSYM=1
         GRPDET=FANT(1)
         IF(IGROUP.EQ.1)                GRPDET = FANT(1)
         IF(IGROUP.EQ.3)                GRPDET = FANT(2)
         IF(IGROUP.EQ.2)                GRPDET = FANT(3)
         IF(IGROUP.EQ.4.AND.NAXIS.EQ.2) GRPDET = FANT(4)
         IF(IGROUP.EQ.8.AND.NAXIS.EQ.2) GRPDET = FANT(5)
         IF(IGROUP.EQ.7.AND.NAXIS.EQ.2) GRPDET = FANT(6)
         IF(IGROUP.EQ.6.AND.NAXIS.EQ.2) GRPDET = FANT(7)
         IF(IGROUP.EQ.9.AND.NAXIS.EQ.2) GRPDET = FANT(8)
         DO I=1,8
            IF (GRPDET.EQ.FANT(I)) THEN
               IGPDET=LFANT(I)
               STSYM = GANT(LGANT(I)+KSTSYM)
            ENDIF
         ENDDO
C
         IROOT=1
         CALL VCLR(WSTATE,1,MXRT)
         WSTATE(1) = ONE
         CALL VICLR(NFLGDM,1,MXRT)
         NFLGDM(1) = 1
         CALL VICLR(IWTS,1,MXRT)
         IWTS(1)=1
C
         MAXW1  = 300    ! NHGSS
         CRIT   = 1D-06  ! CVGTOL
         NITER  = 100    ! ITERMX
         K      = 1      ! NSTATE
         KST    = 1      ! NSTGSS
         MAXP   = 10     ! MXXPAN
         PRTTOL = 0.05D+00
C
         ICLBBR=0
         ANALYS=.FALSE.
C
         IF(MASWRK) WRITE(IW,9060) IEXCIT,CICODE,NCORSV,
     *                             NA-NCORSV+NVVOS,
     *                             NA-NCORSV,NVVOS,NELS
C
         IF(SOCI) THEN
C              DEFINE SINGLE REFERENCE "MR-CISD" USING FSOCI CODE
            NEXT   = NVVOS
            MOSET  = 0
            MAXPSO = MAXP
            KSO    = K
         ELSE
C              DEFINE TWO ORMAS SPACES: NA (OF INFOA) INCLUDES CORE E-
            NSPACE  = 2
            MSTA(1) = NCORSV+1
            MSTA(2) = NA+1
            MSTA(3) = NA+NVVOS+1
            MNUM(1) = NA-NCORSV
            MNUM(2) = NVVOS
            MINI(1) = NELS-IEXCIT
            MINI(2) = 0
            MAXI(1) = NELS
            MAXI(2) = IEXCIT
            FDIRCT  = .TRUE.
            QCORR   = .FALSE.
            CALL FCCHECK(IW,.FALSE.,NSPACE,MNUM,MINI,MAXI,
     *                   IAMI,IAMA,IBMI,IBMA,NADET,NBDET)
         END IF
C
                    DIRTRF=DIRSCF
         IF(GOPARR) DIRTRF=.TRUE.
C
         ABEL = ABELPT()
         IF(.NOT.ABEL  .AND.  .NOT.DIRTRF) THEN
            IF(MASWRK) WRITE(IW,9070)
            CALL SYMOFF
            CALL JANDK
            CALL SYMON
         END IF
C
C           PREPARE ORBITAL SYMMETRY LABELS
C
         L0 = NQMT
         L1 = NUM
         L2 = (L1*L1+L1)/2
         L3 =  L1*L1
C
         CALL VALFM(LOADFM)
         LMOLAB = LOADFM + 1
         LMOIRP = LMOLAB + L1
         LVEC   = LMOIRP + L1
         LS     = LVEC   + L3
         LQ     = LS     + L2
         LWRK   = LQ     + L3
         LMODEG = LWRK   + L1
         LAST   = LMODEG + L1
         NEEDL  = LAST - LOADFM - 1
         CALL GETFM(NEEDL)
C
         IF(GRPDET.EQ.FANT(1)  .OR.  EXETYP.EQ.CHECK) THEN
            CALL C1DET(XX(LMOIRP),XX(LMOLAB),L0)
         ELSE
            CALL DAREAD(IDAF,IODA,XX(LVEC),L3,15,0)
            CALL DAREAD(IDAF,IODA,XX(LS),L2,12,0)
            CALL DAREAD(IDAF,IODA,XX(LQ),L3,45,0)
            CALL TRFSYM(XX(LMOLAB),XX(LMOIRP),XX(LMODEG),XX(LQ),XX(LS),
     *                  XX(LVEC),XX(LWRK),IA,L0,L1,L0,L1)
         END IF
         CALL GAJASW(XX(LMOIRP),NUM,GRPDET)
         CALL DAWRIT(IDAF,IODA,XX(LMOIRP),L1,262,1)
         CALL RETFM(NEEDL)
C
C           CARRY OUT THE CI COMPUTATION
C
         IF(SOCI) THEN
            CALL FSODCI(NRNFG,NPFLG)
         ELSE
            CALL ORDET(NRNFG,NPFLG)
         END IF
      END IF
C
C        --- OBTAIN THE LOCALIZATION TRANSFORMATION ONTO ATOMS ---
C
      CALL LMOX
C
C        --- ROTATE DENSITY MATRIX TO THE LOCALIZED ORBITAL BASIS ---
C        -M1- INCLUDES OCCUPIED AND VIRTUAL VALENCE, BUT NOT CORE
C
      MDOC = NA - MCORE
      M1 = MDOC + NVVOS
      M2 = (M1*M1+M1)/2
      M3 = M1*M1
C
      CALL VALFM(LOADFM)
      LDCMO = LOADFM + 1
      LDLMO = LDCMO  + M2
      LTRAN = LDLMO  + M2
      LWRK  = LTRAN  + M3
      LAST  = LWRK   + M1
      NEED = LAST - LOADFM - 1
      CALL GETFM(NEED)
      CALL BNDANADEN(DMATRX,XX(LDCMO),XX(LDLMO),XX(LTRAN),XX(LWRK),
     *               M1,M2,MDOC)
      CALL RETFM(NEED)
C
C        --- ORIENT THE QUASIATOMIC ORBITALS AND PRINT ANALYSIS ---
C
      CALL DIRLMO
C
      RETURN
C
 9000 FORMAT(1X,'ERROR READING $BNDANA GROUP, TRY AGAIN')
 9010 FORMAT(1X,'*** ERROR ***'/
     *       1X,'THE USE OF DENS=ORMAS REQUIRES CITYP=ORMAS, AS WELL'/
     *       1X,'AS $CIDET AND $ORMAS INPUT GROUPS TO DEFINE THE CI.')
 9020 FORMAT(1X,'INVALID CHOICE FOR DENSITY MATRIX=',A8/
     *       1X,'PLEASE CHOOSE ONLY FROM: ',A8,4(',',A8),'.')
 9030 FORMAT(/5X,24(1H=)/5X,'BONDING ANALYSIS DETAILS'/5X,24(1H=)/
     *     1X,'THE ANALYSIS WILL BE BASED ON ',A8,
     *        ' LOCALIZATION OF ',A8,' ORBITALS,'/
     *     1X,'WHICH COMPRISE',I5,' OCCUPIED VALENCE AND',I5,
     *        ' VIRTUAL VALENCE ORBITALS,'/
     *     1X,'USING A DENSITY MATRIX OBTAINED AT THE ',A6,' LEVEL.')
 9040 FORMAT(/1X,'USING THE SIMPLISTIC SCF LEVEL DENSITY IS NOT AS',
     *           ' REVEALING AS A CI DENSITY.')
 9050 FORMAT(/1X,'THE CUSTOMIZED CI COMPUTATION IS ALREADY FINISHED,',
     *           ' AFTER THE SCF CONVERGED.')
 9060 FORMAT(/1X,'PERFORMING THE CI CALCULATION AT EXCITATION',
     *           ' LEVEL',I2,', USING ',A,','/
     *        1X,'WITH',I5,' CHEMICAL CORE AND',I5,' VALENCE ORBITALS',
     *           ' (OCC=',I5,', VIRT=',I5,'),'/
     *        1X,'CORRELATING ALL',I5,' VALENCE ELECTRONS.')
 9070 FORMAT(/1X,'REGENERATING AO INTEGRALS WITHOUT EXPLOITING POINT',
     *           ' GROUP SYMMETRY'/
     *        1X,'FOR USE BY THE NON-ABELIAN GROUP INTEGRAL',
     *           ' TRANSFORMATION')
      END
C
C*MODULE VVOS    *DECK BNDANADEN
      SUBROUTINE BNDANADEN(DMATRX,DCMO,DLMO,TRAN,WRK,M1,M2,MDOC)
      use constants, only: two
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION DCMO(M2),DLMO(M2),TRAN(M1,M1),WRK(M1)
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON /RUNOPT/ RUNTYP,EXETYP,NEVALS,NGLEVL,NHLEVL
C
      DATA SCF/8HSCF     /
      DATA CHECK/8HCHECK   /
C
C         PREPARE THE USE OF JOE'S ORIENTATION PROGRAM BY GETTING
C         THE CORRECT DENSITY MATRIX TO THE CORRECT DISK RECORD.
C
      IF(EXETYP.EQ.CHECK) GO TO 800
C
C         READ THE TRANSFORMATION TO LOCALIZED ORBITALS
C
      CALL DAREAD(IDAF,IODA,TRAN,M1*M1,73,0)
C
C         READ (OR GENERATE) THE DENSITY MATRIX
C
      IF(DMATRX.EQ.SCF) THEN
         CALL VCLR(DCMO,1,M2)
         DO I=1,MDOC
            II = (I*I+I)/2
            DCMO(II) = TWO
         ENDDO
      ELSE
         CALL DAREAD(IDAF,IODA,DCMO,M2,320,0)
      END IF
C
C          ROTATE THE DENSITY TO THE LOCALIZED ORBITAL BASIS
C
      CALL TFTRI(DLMO,DCMO,TRAN,WRK,M1,M1,M1)
C
C          WRITE OUT THE LOCALIZED DENSITY TO DAF RECORD 285
C
  800 CONTINUE
      IF(EXETYP.EQ.CHECK) CALL VCLR(DLMO,1,M2)
      CALL DAWRIT(IDAF,IODA,DLMO,M2,285,0)
      RETURN
      END
C*MODULE VVOS    *DECK ETSHMCP
!>
!>    @brief   This is a clone of ETSHELL adapted for use with
!>             model core potentials
!>
!>    @details This subroutine is called instead of ETSHELL when a
!>             model core potential is used. This sets up information
!>             about each AAMBS orbital. Since this is a clone of
!>             much older code from vvos.src, some degree of older
!>             style of programming is retained such as GOTO
!>             statments. I see no sense in breaking something that
!>             works.
!>
!>    @author  George Schoendorff
!>     - October 18, 2017
!>
!>    @date 11/21/22 - George Schoendorff
!>     - Contractions lengths obtained via AAMBS_CONTRACTION_LENGTH
!>
!>    @param NUCZ     is the true nuclear charge of the atom
!>    @param DZ       is a flag to determine if we float the
!>                    outer exponents
!>    @param NPASS    is the number of AAMBS orbitals for the atom
!>    @param LENILIST is the maximum number of shells
!>    @param IELIST   is an array specifying the starting index of
!>                    the exponents for each orbital type
!>    @param ICLIST   is an array specifying the starting index of
!>                    the coefficients for each orbital
!>    @param ITLIST   is an array specifying the orbital types
!>    @param IGLIST   is an array specifying the number of Gaussian
!>                    exponents per orbital
!>
      SUBROUTINE ETSHMCP(NUCZ,DZ,NPASS,
     *                   LENILIST,IELIST,ICLIST,ITLIST,IGLIST)
C
      USE MX_LIMITS, ONLY: MXGTOT,MXSH,MXAO,MXATM,MXRT,MXNORO
      USE comm_VVOPAR, only: ivvtyp
C
      IMPLICIT NONE
C
      INTEGER NUCZ,NPASS,LENILIST,IELIST,ICLIST,ITLIST,IGLIST
      INTEGER I,NS,NP,ND,NF,NERR
      INTEGER AAMBS_CONTRACTION_LENGTH
      DIMENSION IELIST(LENILIST),ICLIST(LENILIST),
     *          ITLIST(LENILIST),IGLIST(LENILIST)
C
      LOGICAL DZ,PBLOCK
      LOGICAL SUBVAL
C
      CHARACTER(LEN=4) FLOAT
C      CHARACTER*4 FLOAT
C
      INTEGER IW,IDAF,IODA,IP,IPK,IR,IS,NAV
C
      INTEGER ME,MASTER,NPROC,IBTYP,IPTIM
      LOGICAL GOPARR,DSKWRK,MASWRK
C
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
C
C------------
      IF(LENILIST.NE.19) THEN
         IF(MASWRK) WRITE(IW,*) 'CHANGE LENILIST SOMEWHERE'
         CALL FLSHBF(IW)
         CALL ABRT
         STOP
      ENDIF
C------------
      NPASS=0
      CALL VICLR(ITLIST,1,LENILIST)
      NS=0
      NP=0
      ND=0
      NF=0
      NS =  AAMBS_CONTRACTION_LENGTH(NUCZ,1,0,1)
      IF((NUCZ.GE.3.AND.IVVTYP.EQ.1).OR.(NUCZ.GE.5.AND.IVVTYP.EQ.0))
     *   NP =  AAMBS_CONTRACTION_LENGTH(NUCZ,2,1,IVVTYP)
      IF(NUCZ.GE.21) ND = AAMBS_CONTRACTION_LENGTH(NUCZ,3,2,1)
      IF(NUCZ.GE.57) NF = AAMBS_CONTRACTION_LENGTH(NUCZ,4,3,1)
C
C     H - He (using the nonrelativistic AAMBS)
C
      IF(NUCZ.LE. 2) THEN
         NS =  8
         NP =  0
         NPASS=1
         ITLIST(1)=1
         IF(DZ) FLOAT='S   '
C
C     Li - Be
C
      ELSEIF(NUCZ.GE.3.AND.NUCZ.LE. 4) THEN
         NPASS=2
         ITLIST(1)=1
         ITLIST(2)=2
         IF(DZ) FLOAT='SP  '
C
C     B - Ne
C
      ELSEIF(NUCZ.GE.5.AND.NUCZ.LE.10) THEN
         NPASS=2
         ITLIST(1)=1
         ITLIST(2)=2
         IF(DZ) FLOAT='SP  '
C
C     Na - Mg
C
      ELSEIF(NUCZ.GE.11.AND.NUCZ.LE.12) THEN
         NPASS=2
         ITLIST(1)=1
         ITLIST(2)=2
         IF(DZ) FLOAT='SP  '
C
C     Al - Ar
C
      ELSEIF(NUCZ.GE.13.AND.NUCZ.LE.18) THEN
         NPASS=2
         ITLIST(1)=1
         ITLIST(2)=2
         IF(DZ) FLOAT='SP  '
C
C     K - Ca
C
C     Omit the 3d AAMBS set here since the MCP-DZP
C     valence basis set does not include d functions
C
      ELSEIF(NUCZ.GE.19.AND.NUCZ.LE.20) THEN
         NPASS=2
         ITLIST(1)=1
         ITLIST(2)=2
C                                    Sc-Ni
      ELSEIF(NUCZ.GE.21.AND.NUCZ.LE.30) THEN
         NPASS=5
         ITLIST(1)=1
         ITLIST(2)=2
         ITLIST(3)=1
         ITLIST(4)=3
         ITLIST(5)=2
C                                    GA-KR
      ELSEIF(NUCZ.GE.31.AND.NUCZ.LE.36) THEN
         NPASS=3
         ITLIST(1)=3
         ITLIST(2)=1
         ITLIST(3)=2
         IF(DZ) FLOAT='SP  '
C
C                                    RB-SR
C
C     Omit the 4d AAMBS set here since the MCP-DZP
C     valence basis set does not include d functions
C
      ELSEIF(NUCZ.GE.37.AND.NUCZ.LE.38) THEN
         NPASS=4
         ITLIST(1)=1
         ITLIST(2)=2
         ITLIST(3)=1
         ITLIST(4)=2
C                                    Y-Cd
      ELSEIF(NUCZ.GE.39.AND.NUCZ.LE.48) THEN
         NPASS=5
         ITLIST(1)=1
         ITLIST(2)=2
         ITLIST(3)=1
         ITLIST(4)=3
         ITLIST(5)=2
C                                    IN-XE
      ELSEIF(NUCZ.GE.49.AND.NUCZ.LE.54) THEN
         NPASS=3
         ITLIST(1)=3
         ITLIST(2)=1
         ITLIST(3)=2
         IF(DZ) FLOAT='SP  '
C
C                                 CS-BA
      ELSEIF(NUCZ.GE.55.AND.NUCZ.LE.56) THEN
        NPASS=4
        ITLIST(1)=1
        ITLIST(2)=2
        ITLIST(3)=1
        ITLIST(4)=2
        IF(DZ) FLOAT='SP  '
C
C                                 LA-YB
      ELSEIF(NUCZ.GE.57.AND.NUCZ.LE.70) THEN
        NPASS = 6
        ITLIST(1)=1
        ITLIST(2)=2
        ITLIST(3)=1
        ITLIST(4)=4
        ITLIST(5)=3
        ITLIST(6)=2
C
C                                 LU-HG
      ELSEIF(NUCZ.GE.71.AND.NUCZ.LE.80) THEN
        NPASS = 5
        ITLIST(1)=1
        ITLIST(2)=2
        ITLIST(3)=1
        ITLIST(4)=3
        ITLIST(5)=2
C
C                                 TL-RN
      ELSEIF(NUCZ.GE.81.AND.NUCZ.LE.86) THEN
        NPASS=3
        ITLIST(1)=3
        ITLIST(2)=1
        ITLIST(3)=2
        IF(DZ) FLOAT='SP  '
      ELSE
        IF(MASWRK) WRITE(IW,*) 'IMPOSSIBLE MISTAKE IN ETSHMCP'
        CALL ABRT
        STOP
      ENDIF
C----------------
      IF(NPASS.EQ.0) THEN
         IF(MASWRK) WRITE(IW,*) 'ETSHMCP ERROR: NPASS CANNOT BE 0.'
         CALL FLSHBF(IW)
         CALL ABRT
         STOP
      ENDIF
      IF(NPASS+2.GT.LENILIST) THEN
         IF(MASWRK) THEN
            WRITE(IW,*) 'ETSHELL ERROR:'
            WRITE(IW,*) 'NPASS+2=',NPASS+2
            WRITE(IW,*) 'LENILIST=',LENILIST
            WRITE(IW,*) 'INCREASE LENILIST FOR MEMORY.  TWO PLACES.'
         ENDIF
         CALL FLSHBF(IW)
         CALL ABRT
         STOP
      ENDIF
CGS MCP potential error out
      NERR=0
      DO I=1,NPASS
         IF(ITLIST(I).EQ.0) NERR=NERR+1
      ENDDO
      IF(NERR.NE.0) THEN
         IF(MASWRK) THEN
            WRITE(IW,*) 'ETSHMCP ERROR:'
            WRITE(IW,*) 'PROBLEM WITH ITLIST BEFORE NON-FLOAT.'
         ENDIF
         CALL FLSHBF(IW)
         CALL ABRT
         STOP
      ENDIF
C----------------
      IELIST(1)=0
      ICLIST(1)=0
CGS
      PBLOCK=.FALSE.
      IF(NUCZ.GE.31.AND.NUCZ.LE.36) PBLOCK=.TRUE.
      IF(NUCZ.GE.49.AND.NUCZ.LE.54) PBLOCK=.TRUE.
      IF(NUCZ.GE.81.AND.NUCZ.LE.86) PBLOCK=.TRUE.
      IF(PBLOCK) THEN
        IGLIST(1)=ND
        DO I=2,NPASS
          IF(ITLIST(I).EQ.3) THEN
             IELIST(I) = 0
             IGLIST(I) = ND
          END IF
          IF(ITLIST(I).EQ.1) THEN
             IELIST(I) = ND
             IGLIST(I) = NS
          END IF
          IF(ITLIST(I).EQ.2) THEN
             IELIST(I) = ND+NS
             IGLIST(I) = NP
          END IF
          IF(ITLIST(I).EQ.4) THEN
            IF(MASWRK) THEN
              WRITE(IW,*) 'ETSHMCP ERROR:'
              WRITE(IW,*) 'NO F ORBITALS IN P-BLOCK MCP'
            ENDIF
            CALL FLSHBF(IW)
            CALL ABRT
            STOP
          ENDIF
          ICLIST(I) = ICLIST(I-1) + IGLIST(I-1)
        ENDDO
      ELSE
        IGLIST(1)=NS
        DO I=2,NPASS
          IF(ITLIST(I).EQ.1) THEN
             IELIST(I) = 0
             IGLIST(I) = NS
          END IF
          IF(ITLIST(I).EQ.2) THEN
             IELIST(I) = NS
             IGLIST(I) = NP
          END IF
          IF(ITLIST(I).EQ.3) THEN
             IELIST(I) = NS+NP
             IGLIST(I) = ND
          END IF
          IF(ITLIST(I).EQ.4) THEN
             IELIST(I) = NS+NP+ND
             IGLIST(I) = NF
          END IF
          ICLIST(I) = ICLIST(I-1) + IGLIST(I-1)
        ENDDO
      ENDIF
CGS
      IF(FLOAT.EQ.'S   ') THEN
         ITLIST(NPASS+1) = 1
         IGLIST(NPASS+1) = 1
         IELIST(NPASS+1) = 0
         ICLIST(NPASS+1) = ICLIST(NPASS) + IGLIST(NPASS)
         NPASS=NPASS+1
      END IF
      IF(FLOAT.EQ.'SP  ') THEN
         ITLIST(NPASS+1) = 1
         ITLIST(NPASS+2) = 2
         IGLIST(NPASS+1) = 1
         IGLIST(NPASS+2) = 1
         IELIST(NPASS+1) = 0
         IELIST(NPASS+2) = NS
         ICLIST(NPASS+1) = ICLIST(NPASS)   + IGLIST(NPASS)
         ICLIST(NPASS+2) = ICLIST(NPASS+1) + IGLIST(NPASS+1)
         NPASS=NPASS+2
      END IF
C
C-    WRITE(6,*) 'EXIT FROM ETSHMCP ROUTINE, NPASS=',NPASS
C-    WRITE(6,899) 'IELIST',(IELIST(I),I=1,NPASS)
C-    WRITE(6,899) 'ICLIST',(ICLIST(I),I=1,NPASS)
C-    WRITE(6,899) 'ITLIST',(ITLIST(I),I=1,NPASS)
C-    WRITE(6,899) 'IGLIST',(IGLIST(I),I=1,NPASS)
C-899 FORMAT(A,'=',5I5/7X,10I5)
C
      RETURN
      END
C

      subroutine print_array(length,array)
      implicit none
      integer, intent(in) :: length
      double precision, intent(in) :: array(*)
      integer i
      INTEGER IW,IDAF,IODA,IP,IPK,IR,IS,NAV
      INTEGER ME,MASTER,NPROC,IBTYP,IPTIM
      LOGICAL GOPARR,DSKWRK,MASWRK
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      if(maswrk) write(iw,9999) (array(i),',',i=1,length-1),array(i)
 9999 format((5000(ES15.8,A1)),ES15.8,'/')
      return
      end
