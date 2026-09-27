C*MODULE REKS    *DECK REKSSCF
C>
C> @brief      REKS SCF driver.
C>
C> @author     Michael Filatov and Seunghoon Lee
C>
C> @date       2021, Oct
C>
C> @details    create ensemble Fock matrices, add DFT contributions,
C>             solve HF equations.
C>
      SUBROUTINE REKSSCF
C
      USE comm_SSR, ONLY: S2SARE,WS2SA, S3SARE, WS3SA, CLXGRD
      USE comm_REKSCM, ONLY: NMICRO, MTTYP, WPPS, WOSS, G1, DNR,
     * DNS, DELTA, FR, FS
      USE comm_REXOPT, ONLY: REXTYPE, REXTARGET, REXSHIFT, REXDIIS,
     * REXLDL, RLXDEN, REXEKT, EKTEA, REXCG, RXCGIT, RXCGTH
      USE metaGGA, ONLY: PRTTAU => printtau
      USE camdft, ONLY: ALPHAC => cam_alpha, BETAC => cam_beta,
     *  CAMMU => cam_mu, CAMFLAG
      USE lrcdft, ONLY: LCFLAG, EMU, EMU2
      USE MX_LIMITS, ONLY: mxsh, mxgtot, mxg2, mxatm, mxao, mxrt, mxgrid
      USE comm_FRGINF
      USE constants, ONLY: zero, half, one, two, ten
      USE prec, ONLY: dp
      IMPLICIT NONE
C
      INTEGER, PARAMETER :: NUMRA = 23
      REAL(KIND=dp), PARAMETER :: PT2 = 0.2D+00
      REAL(KIND=dp), PARAMETER :: TWOPT2 = 2.2D+00
C
C
      REAL(KIND=dp), DIMENSION(MXATM) :: ANAM, BNAM, ZAN
      REAL(KIND=dp), DIMENSION(4) :: ASCREEN, RESTRI, RFLMO
      REAL(KIND=dp), DIMENSION(MXAO) :: BFLAB
      REAL(KIND=dp) :: BFSZGB, CCTYP, CITYP, CNVAFO, CNVDMP, CONVFG,    &
     &                 CONVHF, COROFF, DCA1, DCA2, DCALP, DCS6, DCS8,   &
     &                 DCSR, DENTOL, DFTGTHR, DFTTHR, DFTYPE,           &
     &                 DIELEC, DIFF, DIRTHR, DMPCUT, DMPTOL,            &
     &                 DPGD, E0NBO, E1, E2, EBORN, ECORE, ECORL, EDISP, &
     &                 EELCT, EERD, EHF, EHF0, EHFGAS, EKIN, ELRD10,    &
     &                 ELRD6, ELRD8, EMP2S, EMULT, EN, ENERGY, ENUCR,   &
     &                 EPOT, EPSLN, ESCF, ESPSCF, ETHRSH, ETOLLZ, ETOT, &
     &                 ETOT1, EXCOR, EXENA, EXENB, EXENC, EXETYP,       &
     &                 EXTTOL, FSHIFT, GCAVP, GCAVS, GDISP,             &
     &                 GDTOLA, GDTOLB, GJTOLA, GJTOLB, GKTOLA, GKTOLB,  &
     &                 GNUCF, GREP, GVIR, GZRF, ORSHFT
      REAL(KIND=dp), DIMENSION(137) :: BSLRD
      REAL(KIND=dp), DIMENSION(3,MXATM) :: C, EGRAD
      REAL(KIND=dp), DIMENSION(MXGTOT) :: CD, CF, CG, CH, CI, CP, CS, EX
      LOGICAL :: DC, DCABC, DCCHG, DENFIT, DIRSCF, DOLRD, DSKWRK,       &
     &           FDIFF, GOPARR, GOSMP, INITGDDI, ISGDDI, LRDFLG, LRINT, &
     &           LVCLN, MASWRK, MINMEM, MLGDDI, MLTINT, OTHAUX, PACK2E, &
     &           PAROUT, QFMM, QOPS, SG1, USEDM, VIROK, VTSCAL, WASGDDI
      REAL(KIND=dp), DIMENSION(20) :: DFTTYP
      REAL(KIND=dp), DIMENSION(2) :: E0SCF, RESPAP, RESPPC
      REAL(KIND=dp), DIMENSION(3) :: EDFT, FIND
      REAL(KIND=dp), DIMENSION(9) :: ESPSCA
      REAL(KIND=dp), DIMENSION(MXRT) :: ESTATE
      REAL(KIND=dp), DIMENSION(4,48) :: EULANG
      REAL(KIND=dp), DIMENSION(48,48) :: GAM
      INTEGER, DIMENSION(MXAO) :: IA, IRMON, LBUFF
      INTEGER :: IAHARD, IAUXBF, IBTYP, IBWM, ICALC, ICALP, ICBET, ICH, &
     &           ICURFG, ICURIT, ICURLAY, ICURPOP, ICURUNT, ICUT, IDAF, &
     &           IDAFMO, IDCVER, IDDCUR, IDDFMO, IDFT34, IDMFMO, IDOCHG,&
     &           IDOCMAP, IDOLJ, IDOPOL, IDOPROP, IDPUNC, IECP, IEFLD,  &
     &           IESDPPC, IEXTIN, IFMOBAS, IFMOSTP, IGOFLG, IGOMIN,     &
     &           IGRDTYP, ILENG, ILOCAL, ININTIC, INTLOC, INTTYP,       &
     &           INVCUT, IP, IPIEDA, IPK, IPTIM, IQP, IQRAF, IR, IRAF,  &
     &           IREST, IRSTLAY, IRSTSTP, IRTCUT, IS, ISCHWZ, ISCOPE,   &
     &           ISCUT, ISKIPESP, IST, ISVP, ITER, ITERMS, ITOL, IVMFMO,&
     &           IVMTD, IVMUL, IW, IWS, IXESP, IZRF, JANS, JBTYP,       &
     &           JCURFG, JST, KCURFG, KST, LABSIX, LBF, LBUFPIC, NSEGM
      INTEGER, DIMENSION(MXATM) :: IAN
      INTEGER, DIMENSION(10) :: IAPOL
      INTEGER, DIMENSION(950) :: IODA
      INTEGER, DIMENSION(NUMRA) :: IORA
      INTEGER, DIMENSION(2) :: ITMFMO
      INTEGER, DIMENSION(MXSH) :: KATOM, KLOC, KMAX, KMIN, KNG, KSTART, &
     &                            KTYPE
      INTEGER :: LIXIC, LST, MASTER, MAXDII,                            &
     &           MAXGRD, MAXIT, MAXL1C, MAXNYP, MAXVT, MAXWS,           &
     &           MCONV, ME, MEGLOB, MEMSH, MEUNIV, MEXSKPGES, MEXSTATE, &
     &           MIJKL, MODESP, MODFD, MODFMM, MODGRD, MODLMO, MODORB,  &
     &           MODPAN, MODPAR, MODPRP, MODRST, MOFOCK, MONCOR, MP2RUN,&
     &           MPCTYP, MPLEVL, MPMTHD, MTHSVP, MUL, MXITFG,           &
     &           MYGROUP, MYWORLD, N1213J, N14J, NA, NACC, NACCT, NANGL,&
     &           NAOTYP, NAPOL, NAT, NAT1E, NATFMO, NAUXFUN, NAUXSHL,   &
     &           NAV, NB, NBDFG, NBODY, NBOND, NBSSE, NBUFMO, NCAV,     &
     &           NCENTM, NCMAP, NCURSH, NCXYZ, NDDLEFT, NDFTFG, NDIHB,  &
     &           NDIHR, NDUALB, NE, NECP, NEEDR, NELERM, NEORUN, NEVALS
      INTEGER, DIMENSION(MXGRID) :: NANGPT, NANGPT0, NLEB, NLEB0
      INTEGER :: NFFAT, NFG, NFMOPAL, NGAB, NGAU, NGLEVL, NGROUPS,      &
     &           NGUESS, NHEX, NHLEVL, NINTIC, NINTIX, NINTMX, NLAYER,  &
     &           NLKQMM, NNGLOB, NOPDEN, NOPK, NORBPROJ, NORMF, NORMP,  &
     &           NP, NPGLOB, NPGP, NPHI, NPHI0, NPRFMO, NPRINT, NPROC,  &
     &           NPRTGO, NPUNCH, NPUNIV, NQMT, NRAD, NRAD0, NREC, NREJ, &
     &           NREJT, NRPA, NRPAT, NS, NSHELL, NSUBGR, NSVP, NTBOX,   &
     &           NTHE, NTHE0, NTMPL, NTUPL, NUM, NUMDLB, NUMRD, NUNESP, &
     &           NVLPL, NWAGG, NWORLDS, NXXIC, NZMTFMO, MOGDDI
      INTEGER, DIMENSION(4) :: NPREO
      REAL(KIND=dp) :: ORSHFT2, RCORSD, RESDIM, RESPCT,                 &
     &                 RHOMIN, RRSHFT, RUNTYP, SCALTE, SCALTT,          &
     &                 SCFTYP, SHIFTO, SHIFTV, SIZE, SOGTOL, STATN,     &
     &                 STOL, SW0, SWDIIS, SWOFF, SZ1, SZZ, TDDFTYP,     &
     &                 TIMLIM, TKTOLA, TKTOLB, TMUX, TMUXD, TMUY, TMUYD,&
     &                 TMUZ, TMUZD, VBTYP, VEE, VEN, VSHTOL, VTCONV,    &
     &                 VTOL, VTOLA, VTOLB
      REAL(KIND=dp), DIMENSION(10) :: POLCHG, TITLE
      REAL(KIND=dp), DIMENSION(1) :: X
      COMMON /ACONV / RRSHFT, EXTTOL, DMPTOL, VSHTOL, IEXTIN
      COMMON /CONV  / DENTOL, EN, ETOT, EHF, EHF0, DIFF, ITER, ICALP,   &
     &                ICBET
      COMMON /DFGRID/ DFTTHR, DFTGTHR, SWOFF, SW0, BSLRD, NDFTFG, NRAD, &
     &                NTHE, NPHI, NRAD0, NTHE0, NPHI0, NANGPT, NANGPT0, &
     &                SG1, JANS
      COMMON /DFLEB / NLEB, NLEB0
      COMMON /DFTDC / DC, IDCVER, DCCHG, DCABC, DCSR, DCS6, DCS8, DCALP,&
     &                DCA1, DCA2
      COMMON /DFTPAR/ DFTTYP, EXENA, EXENB, EXENC, IDFT34, NAUXFUN,     &
     &                NAUXSHL
      COMMON /DIISSO/ SOGTOL, ETHRSH, MAXDII, IRAF
      COMMON /DMPING/ SHIFTO, SHIFTV, DMPCUT, SWDIIS, DIRTHR
      COMMON /EFPBUF/ POLCHG, NBUFMO, LBUFF, LBF, NAPOL, IAPOL
      COMMON /ENEDA / E0NBO
      COMMON /ENRGYS/ ENUCR, EELCT, ETOT1, SZ1, SZZ, ECORE, ESCF, EERD, &
     &                E1, E2, VEN, VEE, EPOT, EKIN, ESTATE, STATN, EDFT,&
     &                EDISP
      COMMON /FFPARM/ NFFAT, NBOND, NANGL, NDIHR, NDIHB, NCMAP, NWAGG,  &
     &                N1213J, N14J, NLKQMM, IDOCHG, IDOPOL, IDOLJ,      &
     &                IDOCMAP
      COMMON /FMCOM / X
      COMMON /FMMDER/ MAXWS, NTBOX, NCXYZ, NTMPL, MAXNYP
      COMMON /FMOINF/ NFG, NLAYER, NATFMO, NBDFG, NAOTYP, NBODY, NSEGM
      COMMON /FMOOPT/ ESPSCA, RESPAP, RESPPC, RESDIM, RESTRI, RCORSD,   &
     &                RESPCT, CONVFG, CNVDMP, COROFF, RFLMO, ORSHFT,    &
     &                ORSHFT2, CNVAFO, ASCREEN, IXESP, MXITFG, NGUESS,  &
     &                NBSSE, MODORB, MODPAR, IRSTSTP, IRSTLAY, NPRFMO,  &
     &                NFMOPAL, MODPRP, MAXL1C, IPIEDA, MODGRD, MODESP,  &
     &                IVMUL, MODLMO, NOPDEN, MOFOCK, MODFD, MODFMM,     &
     &                NCENTM, NDUALB, NGAB, MODPAN
      COMMON /FMORUN/ ESPSCF, E0SCF, EMP2S, IDAFMO, ICURFG, JCURFG,     &
     &                KCURFG, ICURLAY, ICURUNT, NAT1E, NCURSH, NGAU,    &
     &                ICURPOP, IFMOSTP, MONCOR, NEEDR, MODRST, NORBPROJ,&
     &                NUNESP, ISKIPESP, IESDPPC, IDOPROP, MP2RUN,       &
     &                ICURIT, IDMFMO, IDDFMO, IDDCUR, NDDLEFT, IVMFMO,  &
     &                NZMTFMO, IFMOBAS, ITMFMO
      COMMON /FUNCT / ENERGY, EGRAD
      COMMON /GDDI  / ISCOPE, NGROUPS, MYGROUP, MEGLOB, NPGLOB, NNGLOB, &
     &                JBTYP, ISGDDI, PAROUT, INITGDDI, WASGDDI, MLGDDI, &
     &                NSUBGR, MEUNIV, NPUNIV, NUMDLB, MYWORLD, NWORLDS, &
     &                MOGDDI
      COMMON /IJPAIR/ IA
      COMMON /INDDIP/ TMUX, TMUY, TMUZ, TMUXD, TMUYD, TMUZD, MINMEM
      COMMON /INFGRD/ RHOMIN, ILENG, MAXGRD
      COMMON /INFOA / NAT, ICH, MUL, NUM, NQMT, NE, NA, NB, ZAN, C, IAN
      COMMON /INT2IC/ NINTIC, ININTIC, NXXIC, LBUFPIC, LIXIC, LABSIX,   &
     &                NINTIX
      COMMON /INTFIL/ NINTMX, NHEX, NTUPL, PACK2E, INTTYP, IGRDTYP
      COMMON /INTOPT/ ISCHWZ, IECP, NECP, IEFLD
      COMMON /IOFILE/ IR, IW, IP, IS, IPK, IDAF, NAV, IODA

      COMMON /LMOEDA/ GDTOLA, GJTOLA, GKTOLA, TKTOLA, VTOLA, GDTOLB,    &
     &                GJTOLB, GKTOLB, TKTOLB, VTOLB, ECORL, EXCOR
      COMMON /LRDISP/ ELRD6, ELRD8, ELRD10, EMULT, LRDFLG, MLTINT, DOLRD
      COMMON /MEXOPT/ MEXSKPGES, MEXSTATE
      COMMON /NEOJOB/ NEORUN, NELERM
      COMMON /NLRCF / LRINT
      COMMON /NSHEL / EX, CS, CP, CD, CF, CG, CH, CI, KSTART, KATOM,    &
     &                KTYPE, KNG, KLOC, KMIN, KMAX, NSHELL
      COMMON /OPTSCF/ DIRSCF, FDIFF
      COMMON /OUTPUT/ NPRINT, ITOL, ICUT, NORMF, NORMP, NOPK
      COMMON /PAR   / ME, MASTER, NPROC, IBTYP, IPTIM, GOPARR, DSKWRK,  &
     &                MASWRK
      COMMON /PCMPRT/ GCAVP, GCAVS, GDISP, GREP, EHFGAS
      COMMON /PRPOPT/ ETOLLZ, ILOCAL, IAHARD
      COMMON /QFMMRA/ IORA, IQRAF, MIJKL
      COMMON /QMFM  / SIZE, EPSLN, DPGD, QFMM, NP, NS, IWS, NPGP,       &
     &                MPMTHD, NUMRD, ITERMS, QOPS, ISCUT
      COMMON /RESTAR/ TIMLIM, IREST, NREC, INTLOC, IST, JST, KST, LST
      COMMON /RUNLAB/ TITLE, ANAM, BNAM, BFLAB
      COMMON /RUNOPT/ RUNTYP, EXETYP, NEVALS, NGLEVL, NHLEVL
      COMMON /SCFOPT/ CONVHF, MAXIT, MCONV, NPUNCH, NPREO, FSHIFT
      COMMON /SIMDAT/ NACC, NREJ, IGOMIN, NRPA, IBWM, NACCT, NREJT,     &
     &                NRPAT, NPRTGO, IDPUNC, IGOFLG
      COMMON /SOSYM / EULANG, GAM, IRMON
      COMMON /SVPOPT/ ISVP, NSVP, MTHSVP, NCAV, NVLPL, IQP
      COMMON /VIRIAL/ SCALTE, SCALTT, GVIR, VTCONV, MAXVT, VTSCAL,      &
     &                VIROK, LVCLN
      COMMON /WFNOPT/ SCFTYP, VBTYP, DFTYPE, TDDFTYP, CITYP, CCTYP,     &
     &                MPLEVL, MPCTYP
      COMMON /ZRFPAR/ GZRF, FIND, GNUCF, EBORN, DIELEC, IZRF, ICALC
C
      LOGICAL :: ABINIT, CVDENS, CVDIIS, CVFOCK, CVGED, CVGING,         &
     &           CVGSVP, CVGSVPDN, CVGSVPDN1, CVGSVPNR, DAMPH, DBG,     &
     &           DOAFO, DODIIS, EXTRAH, FMOEX, FMOPL, GMS_CCHEM_RI,     &
     &           IPCFP, IZRFLS, LCFLAGS, LOCOPT, LRINTS, MASPRT, NOTOPN,&
     &           PRDENS, RESET, RSTRCT, SCHWRZ, SOME, SOSCF, SOSCF0,    &
     &           SWGROK, TDSKWRK, VSHIFT
      REAL(KIND=dp) :: AFON, AMINB2, APLUSB, AU2DEBYE, BFON,            &
     &                 CFACT, CFTDES, CFTOSS, CFTPPS, CFTPPSOSS, CFTRAN,&
     &                 CFTRAN2, CFTTRN, CFTTRN1, CFTTRN2, CMTMP, CNVLIM,&
     &                 CONST, DAMP, DAMP0, DEAVG, DELE, DELE0, DET,     &
     &                 DIFFA, DIFFATMP, DIFFB, DIFFP, DIISDMP, DIITOL,  &
     &                 DIPDESX, DIPDESY, DIPDESZ, DIPOSSX, DIPOSSY,     &
     &                 DIPOSSZ, DIPPPSX, DIPPPSY, DIPPPSZ, DTOL2, DUMMY,&
     &                 EDES, ENGTOL, EOSS, EPPS, ERDIIS, EREKS, ERRORC, &
     &                 ERRVAL, ETHNEW, ETHSAV, ETOL2, ETOT0, ETOTD,     &
     &                 ETOTTMP, ETRP, FONA, FONB, GRIDTYP, GRMS, HFSSAV,&
     &                 OFFFOCK, OFFTOL, ORBGRD, REXENERGY, S2,          &
     &                 SOGNEW, SOGSAV, SWDFT, SZ, TEMP, TFITER, TFOCK,  &
     &                 TFOCK1, TFOCKN, THRESH, TIM0, TIM1
      REAL(KIND=dp) :: DDOT, REXEM2EE, REXOFFFOC, TRACEP
      LOGICAL, SAVE :: FT15OP
      INTEGER :: I, ICAB, ICCFLG, IDFTMEM, IDFTSWITCH, IDFTSWITCH_SAVE, &
     &           IFPK, IGRDDFT, II, IJL, IJSPP, IJTBL, IJTBLP,          &
     &           IJTPP, INDX, INFO, INFOSI, IOUT, ISWOFF, ITDIIS,       &
     &           ITERLV, ITERV, ITMP, ITNOND, ITSO, J, JJ, JL1, K, KK,  &
     &           KLSPP, KLTPP, L, L0, L1, L2, L3, LADIIS, LAFA, LAFB,   &
     &           LAOC, LAST, LAST1, LAST2, LBDIIS, LBOC, LCEIG,         &
     &           LCLM, LCM, LCXYZ, LDA, LDB, LDDIJ, LDISPLO, LDLIST,    &
     &           LDLN, LDSH, LDT, LDT1, LDT2, LDT3, LDT4, LDTTMP, LEA,  &
     &           LEB, LEBOX, LEDFA, LEDFB, LEEXC, LEEXCTMP, LEM,        &
     &           LENLINE, LENPRNT, LENSVP, LF, LFA, LFAM, LFAMMO
      INTEGER :: LFAMMOTMP, LFAMTMP, LFATMP, LFA_ACCUM, LFB, LFBM,      &
     &           LFBMMO, LFBMMOTMP, LFBMTMP, LFBTMP, LFB_ACCUM, LFLM,   &
     &           LFXCA, LFXCB, LG, LGHOND, LGRADA, LGRADB, LGRADO, LH1, &
     &           LHESSO, LIBS, LIDXBOX, LIDXIJK, LIDXSHL, LIDXWS, LINDX,&
     &           LINDX2, LINDX3, LIODII, LIPP, LISP, LITPP, LITPP2,     &
     &           LITSP, LITSP2, LIWRK, LIXPK, LIYP, LJT, LJT1, LJT2,    &
     &           LJT3, LJT4, LJTT, LKT, LKT1, LKT2, LKT3, LKT4, LKTT,   &
     &           LMAX, LMAXIJ, LMBOX, LMLIST, LMLPNT, LNBR, LNUMWS,     &
     &           LOADFM, LPDIIS, LPGRADO, LPNTR, LPTBL, LPTOT2, LPUN,   &
     &           LQ, LS, LSCR, LSLIST, LSLN, LSVP, LTMPGPL, LTMPGPS,    &
     &           LTS, LTWOEI, LVA, LVB, LWRK1, LWRK2, LWRK3
      INTEGER :: LWRK4, LWRK5, LWRK6, LWRK7, LXDIIS, LXDINT, LXINTS,    &
     &           LXK, LXP, LYDINT, LYP, LYZPNT, LYZTBL, LZDINT, LZLL,   &
     &           LZP, MAXG, MAXIT2, MAXITREX, MEMDFT, MINTMX, MNEARJ,   &
     &           MPUNCH, MQBOX, MTIJ, NANGM, NBOX, NDAF, NDAFA, NDAFB,  &
     &           NDER, NDFTSV, NEED, NEED1, NEED2, NEEDI, NFT15, NFT16, &
     &           NGOTMX, NGRDMEM, NINT, NMIC, NOCP, NONDMX, NPRA, NPRB, &
     &           NPRO, NPRU, NSAVE, NSCHWZ, NSH2, NSHL2
      INTEGER, DIMENSION(20) :: NSBOX
      CHARACTER(4) :: REXPREFIX
      REAL(KIND=dp), DIMENSION(3) :: TIMSTR
      REAL(KIND=dp) :: TITER0, TITER1, TITER2, TLEFT, TLEFTS, TMP,      &
     &                 TOTELE, TOTKIN, TRDIPT, TRDIPX, TRDIPY, TRDIPZ,  &
     &                 TSITER, TSOLVE, VSHFTVAL, WRS, XCOM, YCOM, ZCOM
C arrays for SI-SA-REKS
      CHARACTER(256) :: DASHES, DDASHES, CONVERT                         !define strings for printing
      CHARACTER(LEN=:), ALLOCATABLE :: LINE
      CHARACTER(LEN=:), ALLOCATABLE :: LINE1
      REAL(KIND=dp), DIMENSION(5) :: TMS2SA                              !temp array for 2SI-2SA-REKS
      REAL(KIND=dp), DIMENSION(9) :: TMS3SA                              !temp array for 3SI-2SA-REKS
      LOGICAL :: cllCPR                                                  !whether to do CP-REKS (for gradient)
C
C
      DATA FT15OP/.FALSE./
      REAL(KIND=dp), PARAMETER :: CHECK = transfer('CHECK   ',1.0d0)
      REAL(KIND=dp), PARAMETER :: DEBUG = transfer('DEBUG   ',1.0d0)
      REAL(KIND=dp), PARAMETER :: DBUGME = transfer('UHFOP   ',1.0d0)
      REAL(KIND=dp), PARAMETER :: ROHF = transfer('ROHF    ',1.0d0)
      REAL(KIND=dp), PARAMETER :: UHF = transfer('UHF     ',1.0d0)
      REAL(KIND=dp), PARAMETER :: RMP = transfer('RMP     ',1.0d0)
      REAL(KIND=dp), PARAMETER :: ZAPT = transfer('ZAPT    ',1.0d0)
c
c       keywords used to call cp-reks routine; more can be added
      REAL(KIND=dp), PARAMETER :: OPTMZE = transfer('OPTIMIZE',1.0d0)
      REAL(KIND=dp), PARAMETER :: HESS = transfer('HESSIAN ',1.0d0)
      REAL(KIND=dp), PARAMETER :: SADPT = transfer('SADPOINT',1.0d0)
      REAL(KIND=dp), PARAMETER :: GRAD = transfer('GRADIENT',1.0d0)
      REAL(KIND=dp), PARAMETER :: AMEX = transfer('MEX     ',1.0d0)
      REAL(KIND=dp), PARAMETER :: MAKEFP = transfer('MAKEFP  ',1.0d0)
      REAL(KIND=dp), PARAMETER :: MD = transfer('MD      ',1.0d0)
      REAL(KIND=dp), PARAMETER :: G3MP2 = transfer('G3MP2   ',1.0d0)
      REAL(KIND=dp), PARAMETER :: ANONE = transfer('NONE    ',1.0d0)
      REAL(KIND=dp), PARAMETER :: COARSE = transfer('COARSE  ',1.0d0)
      REAL(KIND=dp), PARAMETER :: FINE = transfer('FINE    ',1.0d0)
      REAL(KIND=dp), PARAMETER :: COMP = transfer('COMP    ',1.0d0)
      REAL(KIND=dp), PARAMETER :: SFCIS = transfer('SFCIS   ',1.0d0)
      REAL(KIND=dp), PARAMETER :: SPNFLP = transfer('SPNFLP  ',1.0d0)
      REAL(KIND=dp), PARAMETER :: REKS = transfer('REKS    ',1.0d0)
      REAL(KIND=dp), PARAMETER :: NONE = transfer('NONE    ',1.0d0)

C
      NMICRO = reshape( (/2,1,2,1, 1,2,1,2, 2,1,1,2, 2,2,1,1 /),        &
     &                  (/ 2,2,4 /) )
C
      cllCPR = RUNTYP.EQ.OPTMZE.OR.RUNTYP.EQ.HESS.OR.RUNTYP.EQ.SADPT
     *         .OR.RUNTYP.EQ.GRAD
C
      vshftVal = rexShift                                                !set the level shift value for REKS
      DELTA = 0.4D+00                                                    !the value of the δ parameter in the REKS coupling term
      WOSS  = 1.d0 - WPPS                                                !ensemble (state averaging) weighting factor for the OSS configuration
      ILENG = 200
      DODIIS=rexDIIS
C
      dashes='-'
      ddashes='='
      do i = 1,255                                                       !populate strings
       dashes=trim(dashes)//'-'
       ddashes=trim(ddashes)//'='
      enddo
C
      if(maswrk)then                                                     !do only if master process
         line = "Performing a "
         if (rexType.eq.2) line = line // "3SI-2"
         if (rexType.eq.1) line = line // "2SI-2"
         !do nothing
         !if (rexType.eq.0) line = line // ""
         if (abs(wpps-1.d0).gt.1.d-6) line = line // "SA-"
         line = line // "REKS(2,2) calculation"
         write (iw,'(/1x,a)') ddashes(:80)
         if (abs(wpps-1.d0).gt.1.d-6) then
            line = line // " for the "
            if (rexTarget.eq.0) line = line // "averaged state"
            if (rexType.eq.0) then
               if (rexTarget.eq.1) line = line // "PPS state"
               if (rexTarget.eq.2) line = line // "OSS state"
            else
               if (rexTarget.eq.1) line = line // "S0 state"
               if (rexTarget.eq.2) line = line // "S1 state"
            endif                                                        !if(rexType.eq.0)then
         endif                                                           !if(abs(wpps-1.d0).gt.1.d-6)then
         write (iw,'(/10x,a)') trim(line)

         if (abs(wpps-1.d0).gt.1.d-6)
     *      write(iw,'(/10x,"Weights: WPPS =",F12.6,", WOSS =",F12.6)')
     *         WPPS,1.d0-WPPS
C
         write (convert(1:),'(a,F6.3,a)') "Applying",rexShift,
     *      " orbital eigenvalue shift, "
         line = trim(convert)
         if (rexDIIS) then
            line = line // "using"
         else
            line = line // "not using"
         endif
         line = line // " DIIS"
         write (iw,'(/10x,a)') trim(line)
C
         if (rexLdL.eq.1) then
            write (iw,'(/10x,"Localizing the REKS active orbitals")')
         elseif (rexLdL.eq.2) then
            write (iw,'(/10x,"De-localizing the REKS active orbitals")')
         endif
C
         write (iw,'(/1x,a)') ddashes(:80)
      endif                                                              !if(maswrk)then
C
clsh  initial CM
      DNR   = 1.0D+00                                                    !set up initial values of the REKS orbital FONs
      DNS   = 1.0D+00-DNR
cclsh test closed
C
      line = " "
      if (denfit) then                                                   !case: density fitting
         line = line // "DF-"
      endif
      if(abs(WPPS-1.d0).GT.1.d-6)then                                    !case: SA-REKS
         line = line // "SA-"
      endif
      line = line // 'REKS(2,2)'                                         !REKS(2,2)
      if (dftype.ne.anone) then
         write(convert(1:),'(a)') dftype
         line = line // "-" // trim(convert)
      endif
      line = trim(adjustl(line))
      lenline = len(line)
C
      IZRFLS = IZRF.NE.0
C
      DBG    = NPRINT.EQ. 5   .AND.   MASWRK  .AND.  NPRTGO.NE.2
      SOME   = NPRINT.NE.-5   .AND.   MASWRK  .AND.  NPRTGO.NE.2
      PRDENS = NPRINT.GT. 8   .AND.   MASWRK  .AND.  NPRTGO.NE.2
      MASPRT = MASWRK  .AND.  NPRTGO.NE.2 .AND.
     *   (NFG.EQ.0.OR.IAND(NPRFMO,3).LE.1).AND. (RUNTYP.NE.MD)
      IF(EXETYP.EQ.DEBUG  .OR.  EXETYP.EQ.DBUGME) THEN
         DBG    = .TRUE. .AND. MASWRK
         SOME   = .TRUE. .AND. MASWRK
         PRDENS = .TRUE. .AND. MASWRK
      END IF
C
      IF(MASPRT  .AND.  (NFG.EQ.0.OR.NPRINT.NE.-5)) THEN
       write(iw,'(/10x,a)')dashes(:lenline)
       write(iw,'(10x,a)')line(:lenline)
       write(iw,'(10x,a)')dashes(:lenline)
      END IF
C
      TIM0 = ZERO                                                        !START SCF CLOCK
      TIM1 = ZERO
      CALL TSECND(TIM0)
      TLEFTS = TIMLIM - TIM0
C
c     EXTRAH = MOD(MCONV,2)  .EQ. 1
c     DAMPH  = MOD(MCONV,4)  .GE. 2
c     RSTRCT = MOD(MCONV,16) .GE. 8
c     SOSCF0 = MOD(MCONV,128).GE.64
c     LOCOPT = MOD(MCONV,256).GE.128
c     reset  = MOD(MCONV,512).GE.256
c     SOSCF  = SOSCF0 .OR.  SWDIIS.GT.0.0D+00
      EXTRAH = .false.
      DAMPH  = .false.
      RSTRCT = .false.
      VSHIFT = .false.
      SOSCF0 = .false.
      LOCOPT = .false.
      reset  = .false.
      SOSCF  = .false.
C
      L0 = NQMT
      L1 = NUM
      L2 = (L1*L1+L1)/2
      L3 = L1*L1
C
      MAXIT2=(MAXIT*MAXIT+MAXIT)/2
C
      NSH2 = (NSHELL*NSHELL+NSHELL)/2
c
c     fmo
      FMOEX=.false.
      FMOPL=.false.
      DOAFO=.false.
c     reks
      NOCP=4
      NMIC=4
c
      L0 = NQMT
      L1 = NUM
      L2 = (L1*L1+L1)/2
      L3 = L1*L1
C
      MAXIT2=(MAXIT*MAXIT+MAXIT)/2
C
      IFPK  = 1                                                          !SET MEMORY
      IF(NOPK.EQ.1) IFPK=0
C
      CALL VALFM(LOADFM)
      CALL GOTFM(NGOTMX)
c     scratch memory
      LSCR  = LOADFM + 1
      LIWRK = LSCR   + 1+6*L1+3*L3                                       !set new length for the -LSCR- array
      LWRK1 = LIWRK  + L1
      LWRK2 = LWRK1  + MAX(L2,L3)
      LWRK3 = LWRK2  + MAX(L2,L3)
      LWRK4 = LWRK3  + MAX(L2,L3)
      LAST  = LWRK4  + MAX(L2,L3)
c     scf memory
      LDT   = LAST                                                       ! c,l,r,d density mat in AO
      LJT   = LDT    + L2*NOCP                                           ! c,l,r,d coulomb mat in AO
      LKT   = LJT    + L2*NOCP                                           ! c,l,r,d exchange mat in AO
      LFAM  = LKT    + L2*NOCP                                           ! microstates alpha fock mat in AO
      LFBM  = LFAM   + L2*NMIC                                           ! microstates beta fock mat in AO
      LFAMMO= LFBM   + L2*NMIC                                           ! microstates alpha fock mat in MO
      LFBMMO= LFAMMO + L3*NMIC                                           ! microstates beta fock mat in MO
      LEM   = LFBMMO + L3*NMIC                                           ! microstates energy
      LCM   = LEM    + NMIC                                              ! microstates C_L's (weighting factors)
      LVA   = LCM    + NMIC                                              ! MO eigenvectors
      LVB   = LVA                                                        ! set bogus vaue for β MOs (no β's here)
      LDA   = LVB    + L3                                                ! ensemble alpha density mat
      LDB   = LDA    + L2                                                ! ensemble beta density mat
      LFA   = LDB    + L2                                                ! ensemble Fock mat
      LFB   = LFA                                                        ! set bogus value for β
      LEA   = LFB    + L2                                                ! MO eigenvalues
      LEB   = LEA                                                        ! set bogus for β
      LAOC  = LEB    + L1                                                ! α occupation numbers
      LBOC  = LAOC   + L1                                                ! β occupation numbers
      LEEXC = LBOC   + L1                                                ! exchange correlation energy
      LEDFA = LEEXC  + NMIC                                              ! exchange correlation energy
      LEDFB = LEDFA  + NMIC                                              ! exchange correlation energy
      LAST  = LEDFB  + NMIC                                              ! end of it

c     integral memory
      NSH2 = (NSHELL*NSHELL+NSHELL)/2
      CALL BASCHK(LMAX)
                    NANGM =  4
      IF(LMAX.EQ.2) NANGM =  6
      IF(LMAX.EQ.3) NANGM = 10
      IF(LMAX.EQ.4) NANGM = 15
      IF(LMAX.EQ.5) NANGM = 21
      IF(LMAX.EQ.6) NANGM = 28
      MAXG = NANGM**4
C
      MINTMX=NINTMX
      IF(NINTIC.NE.0) MINTMX=0
C
C     MEMORY IS ALLOCATED ELSEWHERE
C     FANCY OPTIONS (LOOK FOR NINTMX BELOW) WILL NOT WORK WITH IN
C     CORE INTEGRAL STORAGE
C
      IF(DIRSCF) THEN
         LXP   = LAST
         LXK   = LAST
         LIXPK = LAST
         LXINTS= LAST
         LGHOND= LXINTS + NSH2
         LDSH  = LGHOND + MAXG
         LDDIJ = LDSH   + NSH2
         LAST  = LDDIJ  + 49*MXG2
      ELSE
         LXP   = LAST
         LXK   = LXP    + MINTMX
         LIXPK = LXK    + MINTMX * IFPK
         LAST  = LIXPK  + MINTMX
         LXINTS= LAST
         LGHOND= LAST
         LDSH  = LAST
         LDDIJ = LAST
      END IF
      IF(NINTIC.NE.0) THEN
         LXP  = LBUFPIC
         LIXPK= LIXIC
      END IF
      LTWOEI = LAST
      LPTOT2 = LAST
C
      IF(DODIIS) THEN
         LADIIS = LAST
         LBDIIS = LADIIS + MAXDII*MAXDII
         LPDIIS = LBDIIS + MAXIT2
         LXDIIS = LPDIIS + MAXIT
         LIODII = LXDIIS + MAXIT
         LAST   = LIODII + 4*MAXDII
      ELSE
         LADIIS = LAST
         LBDIIS = LAST
         LPDIIS = LAST
         LXDIIS = LAST
         LIODII = LAST
      END IF
C
       IPCFP=.FALSE.
       LSVP=LAST                                                         !block svp solvent model
       LENSVP=0
C
C         FOR SOSCF CONVERGER
C
      NFT15=15
      NFT16=16
      NPRA=(L0-NA)*NA
      NPRB=(L0-NB)*NB
      NPRU=NPRA+NPRB
      NPRO=NPRA+(NA-NB)*NB
      ITSO=0
      ORBGRD=ZERO
      LGRADA = LAST
      LGRADB = LAST
      LGRADO = LAST
      LHESSO = LAST
      LPGRADO= LAST
      LDISPLO= LAST

      LCXYZ= LAST
      LIYP = LAST
      LIBS = LAST
      LISP = LAST
      LIPP = LAST
      LIDXWS = LAST
      LIDXIJK= LAST
      LPTBL  = LAST
      LYZTBL = LAST
      LCLM   = LAST
      LFLM   = LAST
      LPNTR  = LAST
      LG     = LAST
      LF     = LAST
      LZLL   = LAST
      LINDX  = LAST
      LINDX2 = LAST
      LINDX3 = LAST
      LDLIST = LAST
      LSLIST = LAST
      LDLN   = LAST
      LSLN   = LAST
      LMLIST = LAST
      LMLPNT = LAST
      LTS    = LAST
      LMAXIJ = LAST
      LIDXSHL= LAST
C
      ISWOFF=0                                                           !ADD IN GRID-BASED DFT MEMORY
      IF(NDFTFG.EQ.1) THEN
         LFXCA    = LAST
         LFXCB    = LFXCA + L2
         IGRDDFT  = LFXCB + L2
C
         IF (CAMFLAG) THEN
            LAFA    = IGRDDFT
            LAFB    = LAFA + L2
            IGRDDFT = LAFB + L2
         ELSE
            LAFA    = IGRDDFT
            LAFB    = IGRDDFT
         END IF
C
         CALL GRDDFT(L2,NGRDMEM)
         LAST     = IGRDDFT + NGRDMEM
         IF((NEVALS.EQ.0 .OR. NFG.NE.0)  .AND.  SWOFF.NE.0) ISWOFF=1
      END IF
      HFSSAV=DFTTYP(3)
      NDFTSV=NDFTFG
      IDFTSWITCH=0
      LCFLAGS=LCFLAG
      LRINTS=LRINT
C
      MEMDFT=LAST                                                        !ADD IN GRID-FREE DFT MEMORY
      IF(DFTTYP(1) .NE. 0.0D+00) THEN
         CALL DFTMEM(IDFTMEM)
         LAST=LAST+IDFTMEM
         NEED  = LAST -LOADFM -1
      END IF
C
      LWRK5=LAST                                                         !allocate memory for cuhf
      LWRK6=LAST
      LWRK7=LAST
      LCEIG=LAST
C
C        WE CAN REDUCE MEMORY NEED BY 2N**2 BY DOING MORE DAREADS,
C        TO OVERLAP -H1- AND -WRK1-, -Q- AND -WRK2-, -S- AND -WRK3-.
C
      NEED = LAST - LOADFM - 1
      NSAVE = 2*L2 + L3
      NEEDI = NEED + NSAVE

      IF((NEEDI.LT.NGOTMX)  .OR.
     *   (VSHIFT .OR. IZRFLS .OR. NFFAT.GT.0)) THEN
         MINMEM = .FALSE.
         LH1   = LAST
         LQ    = LH1    + L2
         LS    = LQ     + L3
         LAST  = LS     + L2
         NEED  = LAST -LOADFM -1
      ELSE
         MINMEM = .TRUE.
         IF(SOME) WRITE(IW,9005) NSAVE
         LH1 = LWRK1
         LQ  = LWRK2
         LS  = LWRK3
      END IF
C
      CALL GETFM(NEED)                                                   !GET MEMORY
C
      II = 0
      DO I = 1,L1
         IA(I) = II                                                      !this array actually saves nothing; better: "ij=0//do i=...//do j==...//ij=ij+1//..."
         II = II + I
      END DO
C
C     ----- INITIALIZE VARIABLES -----
C
      ERDIIS = ZERO
      ITDIIS = 1
      ETHSAV = ETHRSH
      ETHNEW = ETHRSH
      SOGSAV = SOGTOL
      SOGNEW = SOGTOL
      DIISDMP= 0.9D+00
      NONDMX = 5
      ITNOND = NONDMX
      NOTOPN = .TRUE.
C
      IF(MAXIT.LE.0) MAXIT = 30
C
C        LOOSEN CONVERGENCE FOR DIRECT SCF BY JUST A BIT
C
      IF(DIRSCF) THEN
         DENTOL = TWO * CONVHF
         ENGTOL = 1.0D-08
         DIITOL = 1.0D-06
         OFFTOL = DENTOL
      ELSE
         DENTOL = CONVHF
         ENGTOL = 1.0D-09
         DIITOL = 1.0D-07
         OFFTOL = DENTOL
      END IF
      IF(NDFTFG.NE.1) THEN
c        SW0=ZERO
c        For multilayer FMO comment out.
      ELSE
         IF(SW0.LT.CONVHF) THEN
            IF(MASWRK) WRITE(IW,9006) SW0,CONVHF
            CALL ABRT
         END IF
      END IF
      SWDFT = SW0
C
C        AT EARLY POINTS IN GEOMETRY OPTS, WE CAN BE A BIT COARSER,
C        WHENEVER THE LAST GEOMETRY'S GRADIENT WAS STILL BIG.
C
      IF((RUNTYP.EQ.OPTMZE  .OR. RUNTYP.EQ.COMP)
     *   .AND. (NEVALS.GT.0.AND.NFG.EQ.0)) THEN
         GRMS = DDOT(3*NAT,EGRAD,1,EGRAD,1)
         GRMS = SQRT(GRMS/(3*NAT))
                               CFACT =  1.0D+00
         IF(GRMS.GT.0.005D+00) CFACT =  5.0D+00
         IF(GRMS.GT.0.020D+00) CFACT = 20.0D+00
         IF(GRMS.GT.0.100D+00) CFACT = 50.0D+00
         DENTOL = CFACT*DENTOL
         ENGTOL = CFACT*ENGTOL
         DIITOL = CFACT*ENGTOL
         SWDFT  = CFACT*SWDFT
         OFFTOL = DENTOL
      END IF
C
      DTOL2 = TWO * DENTOL
      ETOL2 = TEN * ENGTOL
      CVDIIS = .FALSE.
      CVGING = .FALSE.
      CVGED  = .FALSE.
C
      CVGSVP    = .FALSE.                                                !CVGSVP KEEPS TRACK OF WHETHER SVP CHARGES/DIPOLES APPEAR CONVERGED
      CVGSVPDN  = .FALSE.
      CVGSVPDN1 = .FALSE.
      CVGSVPNR  = .FALSE.
C
      ICCFLG=0
      MPUNCH = NPUNCH
      EHF = ZERO
      ETOT= ZERO
      ECORL=ZERO
      ITERV= 0
      ICALP = 0
      ICBET = 0
      ICAB  = 0
      DAMP  = ZERO
      DAMP0 = ZERO
      IF(DMPCUT .LE. ZERO) DMPCUT = ZERO
      RRSHFT = ZERO
      ITERLV = 0
      DIFF = ZERO
      DIFFP= ZERO
      DIFFA= ZERO
      DIFFB= ZERO
      DELE  = ZERO
      DEAVG = ZERO
      NDAFA = 23
      NDAFB = 26
c     IF(DAMPH .OR. VSHIFT) DAMP = ONE
      SWGROK=.TRUE.
C
      CALL REXCM(X(LCM),DNR,DNS,DELTA,WPPS,WOSS,NMIC)                    !initiate C_L's
C
      IF(MASPRT.AND.(NFG.EQ.0.OR.NPRINT.NE.-5)) THEN                     !PRINT NUCLEAR ENERGY, CONTROL PARAMETERS
         IF(SOME) WRITE(IW,9020) EN,MAXIT,MPUNCH,MUL,EXTRAH,DAMPH,
     *                  VSHIFT,RSTRCT,DODIIS,SOSCF
         IF(NDFTFG.NE.1) THEN
            IF(MASPRT) WRITE(IW,9012) DENTOL,OFFTOL
         ELSE
            IF(MASPRT) WRITE(IW,9013) DENTOL,OFFTOL,SWDFT,SWOFF
         END IF
         IF(SOME  .AND.  VTSCAL) WRITE(IW,9015) VTSCAL,MAXVT,VTCONV
         IF(SOME) WRITE(IW,9040) NEED
      END IF
C
      NINT=0
      NSCHWZ=0
      IF(DIRSCF.AND.ABINIT) THEN                                         !EXCHANGE INTEGRALS FOR DIRECT SCF THRESHOLD TESTS
         SCHWRZ = ISCHWZ.GT.0
         IF(MASPRT) WRITE(IW,9045) SCHWRZ,FDIFF
         IF(SCHWRZ) CALL EXCHNG(X(LXINTS),X(LGHOND),X(LDDIJ),
     *                          NSH2,MAXG,INTTYP)
      ELSE
         SCHWRZ=.FALSE.
      END IF
C
C     L0 = NUMBER OF CANONICAL ORTHONORMAL VECTORS KEPT.
C
      IF(.NOT.MINMEM) THEN
         CALL DAREAD(IDAF,IODA,X(LH1),L2,11,0)                           !READ CORE HAMILTONIAN
         CALL DAREAD(IDAF,IODA,X(LS) ,L2,12,0)                           !READ OVERLAP MATRIX
         CALL DAREAD(IDAF,IODA,X(LQ) ,L3,45,0)                           !READ SYMMETRY ADDAPTED Q MATRIX
      END IF
C
      CALL DAREAD(IDAF,IODA,X(LDA),L2,16,0)                              !READ INITIAL DENSITY
C
      IF(DIRSCF  .AND.  FDIFF) THEN                                      !PREVIOUS DENSITY IS ZERO, PREVIOUS FOCK MAT. IS HCORE
         IF(MINMEM) CALL DAREAD(IDAF,IODA,X(LH1),L2,11,0)
         CALL DAWRIT(IDAF,IODA,X(LH1),L2,14,0)
         CALL DAWRIT(IDAF,IODA,X(LH1),L2,18,0)
         CALL VCLR(X(LWRK1),1,L2)
         CALL VCLR(X(LWRK2),1,L2)
      END IF
C
      IF(IZRFLS) THEN                                                    !IF SCRF: READ IN DIPOLE INTEGRALS
         CALL DAREAD(IDAF,IODA,X(LXDINT),L2,95,0)
         CALL DAREAD(IDAF,IODA,X(LYDINT),L2,96,0)
         CALL DAREAD(IDAF,IODA,X(LZDINT),L2,97,0)
      END IF
C
      IF(IZRFLS  .OR.  IEFP.EQ.1  .OR.  IPCFP  .OR.  ISVP.EQ.1
     *   .OR. NFFAT.GT.0) THEN                                           !IF SCRF: SAVE PRISTINE ONE ELEC. HAM.
         CALL DAWRIT(IDAF,IODA,X(LH1),L2,87,0)
      END IF
C
      IF (NDFTFG.EQ.1) THEN                                              !SETTING FOR GRID DFT
         IF(SG1  .OR.  JANS.GT.0) THEN
            GRIDTYP=FINE
         ELSE
            IF (NLEB(1).NE.0) THEN
               IF(NLEB(1)*NRAD.GT.NRAD0*NLEB0(1)) THEN
                  IF(MASPRT) WRITE(IW,9335) NRAD,NLEB(1),NRAD0,NLEB0(1)
                  IF(EXETYP.NE.CHECK) THEN
                     CALL SWGRID
                     GRIDTYP=COARSE
                  ELSE
                     GRIDTYP=FINE
                  END IF
               ELSE
                  GRIDTYP=FINE
               END IF
            ELSE
               IF(NRAD*NPHI*NTHE.GT.NRAD0*NPHI0*NTHE0) THEN
                  IF(MASPRT) WRITE(IW,9340) NRAD,NTHE,NPHI,NRAD0,NTHE0,
     *                 NPHI0
                  IF(EXETYP.NE.CHECK) THEN                               !POSSIBLE TDDFT CHECK RUN WILL NEED TO SEE TRUE GRID SIZE
                     CALL SWGRID
                     GRIDTYP=COARSE
                  ELSE
                     GRIDTYP=FINE
                  END IF
               ELSE
                  GRIDTYP=FINE
               END IF
            ENDIF
         ENDIF                                                           !IF (NLEB(1).NE.0) THEN
         CALL DFTSET(X(IGRDDFT),1,.FALSE.)
         IF(ISWOFF.GT.0) THEN                                            !SWITCH OFF DFT AT THIS POINT IF WE INITIATE BY HF ITERS
            IF(MASWRK  .AND.  NDFTFG.NE.0) WRITE(IW,9350)
            DFTTYP(3)=1.0D+00
            NDFTFG=0
            LCFLAG=.FALSE.
            LRINT=.FALSE.
         END IF
      END IF                                                             !IF (NDFTFG.EQ.1) THEN
C
      IF (QFMM) THEN                                                     !SETTING FOR QFMM
         CALL QFMMUP(L2,LS,NSHL2,LSLIST,NSHELL,LSLN,LTS,
     *     NCXYZ,LCXYZ,LIYP,LIBS,LISP,LIPP,LIDXWS,LINDX,MAXWS,
     *     MAXNYP,NSH2,LIDXSHL,LPNTR,NBOX,NTBOX,LCLM,LFLM,
     *     LAST,LIDXBOX,LMBOX,LEBOX,LNBR,LNUMWS,LYZPNT,LTMPGPS,
     *     LTMPGPL,NEED1,NEED,LAST1,LINDX2,LINDX3,
     *     LIDXIJK,LPTBL,LYZTBL,NSBOX,NTMPL,
     *     MTIJ,MIJKL,LYP,LZP,LITSP,LITPP,LITSP2,
     *     LITPP2,IJTPP,KLTPP,IJSPP,KLSPP,IJTBLP,IJTBL,LAST2,NEED2)
C
C           THERE ARE ADDITIONAL MEMORY ALLOCATIONS INSIDE ROUTINES
C           -QFMMBOX- AND -NEARJ-, WHICH THE CHECK RUN WILL BRANCH
C           AROUND W/O CALLING WHEN BRANCHING PAST THE ITERS BELOW.
C           THESE MEMORY NEEDS ARE NOT SIMULTANEOUSLY NECESSARY.
         MQBOX = 3*NCXYZ + 2*NCXYZ+ 3*2**NS+ 3
         MNEARJ= 6*MIJKL + 2
         IF(EXETYP.EQ.CHECK) THEN
            NEED2 = NEED2 + MAX(MQBOX,MNEARJ)
         ELSE
            IF(MASPRT) WRITE(IW,9035) NEED2 + MAX(MQBOX,MNEARJ)
         END IF
         CALL GETFM(NEED2)
         NEED=NEED+NEED2
      END IF                                                             !IF (QFMM) THEN
C
      IF(EXETYP.NE.CHECK) THEN                                           !go to line 1998
C
      SOSCF = SOSCF0                                                     !save SOSCF flag; SOSCF is not used here
C
      IF(MASPRT) THEN                                                    !PRINT ITERATION HEADER
         IF(NBUFMO.GT.0) THEN
           WRITE(IW,9039)
         END IF
       if(dirscf) then
        write(iw,'(1x,a)')ddashes(:132)
        write(iw,9055)
        write(iw,'(/1x,a)')ddashes(:132)
       else
        write(iw,'(1x,a)')ddashes(:116)
        write(iw,9050)
        write(iw,'(/1x,a)')ddashes(:116)
       endif
      END IF                                                             !IF(MASPRT) THEN
C
C           *************************
C       ***** START REKS ITERATIONS *****
C           *************************
C
      TFOCK = ZERO
      TSOLVE = ZERO

      DO ITER=1,MAXIT
C
C     ----- CONSTRUCT SKELETON 2 ELECTRON FOCK MATRIX -----
C     FOR DIRECT SCF, EITHER -LWRK1- AND -LWRK2- SHOULD BE THE CHANGE IN THE DENSITY,
C     WITH THE CHANGE IN THE FOCK MATRICES TO BE COMPUTED AT -LFA- AND -LFB-,
C     OR THEY SHOULD CONTAIN THE TOTAL DENSITY WITH THE ENTIRE TWO ELECTRON FOCK
C     OPERATORS TO BE FORMED.
C
  310    CALL TSECND(TITER0)

         IF(.NOT.(DIRSCF))THEN
C
C        ORDINARY INTEGRAL CONTRIBUTIONS TO FOCK MATRIX,
C        WHICH ARE ALSO NEEDED FOR RANGE-SEPARATED DFT.
C
            IF(ITER.EQ.1) CALL DAREAD(IDAF,IODA,X(LVA),L3,15,0)          !read in the eigenvectors
            CALL REXDMAT(X(LDT),X(LVA),IA,L1,L2,L3,NA)                   !build the REKS density matrices
C
            IF(DBG) THEN
               WRITE(IW,*) 'GUESS MO'
               CALL PREVS(X(LVA),X(LEA),X(LIWRK),L0,L1,L1,1)
               DO L=1,NOCP
                  LDTTMP=LDT + (L-1)*L2
                  WRITE(IW,'(1X,I1,A)') L,'TH OCCUP'
                  WRITE(IW,*) 'DENSITY MAT'
                  CALL PRTRIL(X(LDTTMP),L1)
               ENDDO
            END IF
C
            TDSKWRK = DSKWRK                                             !save the current DSKWRK and set it to .true.
            DSKWRK  = .TRUE.
            CALL SEQREW(IS)
C
            LDT1 = LDT
            LDT2 = LDT1 + L2
            LDT3 = LDT2 + L2
            LDT4 = LDT3 + L2
            CALL VCLR(X(LJT),1,L2*NOCP)
            CALL VCLR(X(LKT),1,L2*NOCP)
            LJT1 = LJT
            LJT2 = LJT1 + L2
            LJT3 = LJT2 + L2
            LJT4 = LJT3 + L2
            LKT1 = LKT
            LKT2 = LKT1 + L2
            LKT3 = LKT2 + L2
            LKT4 = LKT3 + L2
c           J_L, c_HF*K_L
            if(LCFLAG) then
               LRINT=.TRUE.
               CALL HSTARREX(X(LDT1),X(LJT1),X(LKT1),X(LDT2),X(LJT2),
     *                       X(LKT2),X(LDT3),X(LJT3),X(LKT3),
     *                       X(LDT4),X(LJT4),X(LKT4),X(LXP),X(LIXPK),
     *                       NINTMX,IA,NOPK,.FALSE.)
               LRINT=.FALSE.
            endif
            IF(CAMFLAG) THEN
               LRINT     = .TRUE.
               EMU       = CAMMU
               EMU2      = CAMMU*CAMMU
               DFTTYP(3) = BETAC
               CALL HSTARREX(X(LDT1),X(LJT1),X(LKT1),X(LDT2),X(LJT2),
     *                       X(LKT2),X(LDT3),X(LJT3),X(LKT3),
     *                       X(LDT4),X(LJT4),X(LKT4),X(LXP),X(LIXPK),
     *                       NINTMX,IA,NOPK,.FALSE.)
               DFTTYP(3) = ALPHAC
               LRINT     = .FALSE.
            endif
            CALL HSTARREX(X(LDT1),X(LJT1),X(LKT1),X(LDT2),X(LJT2),
     *             X(LKT2),X(LDT3),X(LJT3),X(LKT3),
     *             X(LDT4),X(LJT4),X(LKT4),X(LXP),X(LIXPK),
     *             NINTMX,IA,NOPK,.FALSE.)
            DSKWRK  = TDSKWRK
         ENDIF                                                           !IF(DIRSCF) THEN
c        F_L = F_L + J_L - c_HF*K_L
         CALL REX2E(X(LFAM),X(LFBM),X(LJT),X(LKT),IA,L1,L2,NMIC)

         IF(GOPARR) THEN
            CALL DDI_GSUMF(1000,X(LFAM),L2*NMIC)
            CALL DDI_GSUMF(1001,X(LFBM),L2*NMIC)
c           CALL DDI_GSUMI(1002,NINT  ,1)                                !we don't have direct scf; these are not needed for disk-based scf
c           CALL DDI_GSUMI(1003,NSCHWZ,1)
         END IF
C
         IF(DBG) THEN
            DO L=1,NOCP
               LJTT=LFAM + (L-1)*L2
               LKTT=LFBM + (L-1)*L2
               WRITE(IW,'(1X,I1,A)') L,'TH MICROSTATE'
               WRITE(IW,*) 'ALPHA SKELETON FOCK MATRIX'
               CALL PRTRIL(X(LJTT),L1)
               WRITE(IW,*) 'BETA SKELETON FOCK MATRIX'
               CALL PRTRIL(X(LKTT),L1)
            ENDDO
         END IF
C
         IF (NDFTFG.EQ.1) THEN                                           !ADD GRID DFT CONTRIBUTION TO 2E- FOCK OPERATOR
           call stopwa(4,0)
           SWGROK=.TRUE.
           IF(ITER.EQ.1) CALL DAREAD(IDAF,IODA,X(LVA),L3,15,0)
C
           DO L=1,NMIC                                                   ! loop over microstates
              MTTYP   = L
              LFATMP  = LFAM  + (L-1)*L2
              LFBTMP  = LFBM  + (L-1)*L2
              LEEXCTMP= LEEXC +  L-1
              CALL TRPOSE(X(LVA),X(LWRK1),L1,L1,1)                       ! TRANSPOSING ORBITALS ALLOWS UNIT STRIDE INSIDE THE DFT CODES
              DOLRD=.FALSE.
              IF(LRDFLG .AND. CVGING) DOLRD=.TRUE.
C
              if(dbg) WRITE(IW,'(1X,I1,A)') L,'TH MICRO DFT'
              CALL DFTEXCOR(X(IGRDDFT),X(LFXCA),X(LFXCB),
     *                      X(LVA),X(LVB),L1,L2,X(LEEXCTMP),
     *                      TOTELE,TOTKIN)
C
              CALL TRPOSE(X(LVA),X(LWRK1),L1,L1,1)                       !TRANSPOSE THE ORBITALS BACK
              CALL SYMH(X(LFXCA),X(LWRK1),IA)
              if(na.eq.1.and.l.eq.4)then                                 !if there are no core electrons, then
                 call vclr(x(lfxcb),1,l2)                                !zero β-Fock matrix XC contribution for the last microstate (triplet)
              else
                 CALL SYMH(X(LFXCB),X(LWRK1),IA)                         !go on as usual otherwise
              endif
              IF(GOPARR) THEN
                 CALL DDI_GSUMF(2319,ECORL,1)
                 CALL DDI_GSUMF(2310,X(LEEXCTMP),1)
                 CALL DDI_GSUMF(2311,TOTELE,1)
                 CALL DDI_GSUMF(1004,X(LFXCA),L2)
                 CALL DDI_GSUMF(1005,X(LFXCB),L2)
              END IF

              IF(DBG) THEN
                WRITE(IW,'(1X,I1,A)') L,'TH MICROSTATE'
                WRITE(IW,*) 'ALPHA DFT FOCK MATRIX'
                CALL PRTRIL(X(LFXCA),L1)
                WRITE(IW,*) 'BETA  DFT FOCK MATRIX'
                CALL PRTRIL(X(LFXCB),L1)
              END IF

C
              IF (PRTTAU) THEN
                IF(GOPARR) CALL DDI_GSUMF(2317,TOTKIN,1)
                WRITE(IW,9311) TOTKIN
              END IF
              IF(.NOT.(ABINIT.AND.DIRSCF.AND.FDIFF)) THEN
                IF(L.EQ.1 .OR. L.EQ.3) ITMP=2
                IF(L.EQ.2)             ITMP=3
                IF(L.EQ.4)             ITMP=4
                LDTTMP=LDT+(ITMP-1)*L2
                X(LEDFA+L-1)=TRACEP(X(LDTTMP),X(LFXCA),L1)/TWO
                IF(L.EQ.1)             ITMP=2
                IF(L.EQ.2 .OR. L.EQ.3) ITMP=3
                IF(L.EQ.4)             ITMP=1
                LDTTMP=LDT+(ITMP-1)*L2
                X(LEDFB+L-1)=TRACEP(X(LDTTMP),X(LFXCB),L1)/TWO
                CALL VADD(X(LFATMP),1,X(LFXCA),1,X(LFATMP),1,L2)
                CALL VADD(X(LFBTMP),1,X(LFXCB),1,X(LFBTMP),1,L2)
              END IF
           ENDDO                                                         !end of loop over microstates
           call stopwa(4,1)
         END IF                                                          !IF (NDFTFG.EQ.1) THEN
C
C           AT THIS POINT, THE SKELETON 2E- FOCK OPERATOR IS FINISHED
C
      IF(DBG) THEN
         WRITE(IW,*) '2E- FOCK'
         DO L=1,NMIC                                                     ! loop over microstates
            LFAMTMP=LFAM + (L-1)*L2
            LFBMTMP=LFBM + (L-1)*L2
            WRITE(IW,'(1X,I1,A)') L,'TH MICROSTATE'
            WRITE(IW,*) 'ALPHA SKELETON FOCK MATRIX'
            CALL PRTRIL(X(LFAMTMP),L1)
            WRITE(IW,*) 'BETA SKELETON FOCK MATRIX'
            CALL PRTRIL(X(LFBMTMP),L1)
         ENDDO
      END IF
C
C     ----- COMPLETE THE FOCK MATRIX -----
C     WE EITHER ADD THE FOCK OPERATOR OF THE PREVIOUS ITERATION TO
C     THE CHANGE IN THE FOCK MATRIX, OR WE ADD THE ONE ELECTRON
C     HAMILTONIAN TO THE FULLY COMPUTED TWO ELECTRON OPERATOR.
C
C
C     FOR DIRECT SCF-COSMO WHERE THE SURFACE CHARGES ARE UPDATED EACH
C     SCF ITERATION, WE NEED THE CHANGE IN H BETWEEN THE CURRENT AND
C     THE PREVIOUS SCF ITERATION, SO THE CURRENT H IS ADDED FROM
C     SECTION 11 AND THE PREVIOUS H SUBTRACTED FROM SECTION 87
C
      IF(.NOT.(ABINIT  .AND.  DIRSCF  .AND.  FDIFF)) THEN
         IF(MINMEM) CALL DAREAD(IDAF,IODA,X(LH1),L2,11,0)
         DO L=1,NMIC                                                     ! loop over microstates
            LFAMTMP=LFAM + (L-1)*L2
            LFBMTMP=LFBM + (L-1)*L2
            CALL VADD(X(LFAMTMP),1,X(LH1),1,X(LFAMTMP),1,L2)
            CALL VADD(X(LFBMTMP),1,X(LH1),1,X(LFBMTMP),1,L2)
         ENDDO
      END IF
      IF(DBG) THEN
         WRITE(IW,*) 'FOCK'
         DO L=1,NMIC                                                     ! loop over microstates
            LFAMTMP=LFAM + (L-1)*L2
            LFBMTMP=LFBM + (L-1)*L2
            WRITE(IW,'(1X,I1,A)') L,'TH MICROSTATE'
c            WRITE(IW,*) 'ALPHA FOCK MATRIX'
c            CALL PRTRIL(X(LFAMTMP),L1)
c            WRITE(IW,*) 'BETA FOCK MATRIX'
c            CALL PRTRIL(X(LFBMTMP),L1)

            CALL CANTOGEN(X(LFAMTMP),X(LWRK1),L1)
            WRITE(IW,*) 'ALPHA FOCK MATRIX'
            CALL PRSQL(X(LWRK1),L1,L1,L1)
            CALL CANTOGEN(X(LFBMTMP),X(LWRK1),L1)
            WRITE(IW,*) 'BETA FOCK MATRIX'
            CALL PRSQL(X(LWRK1),L1,L1,L1)
         ENDDO
      END IF
C
      CALL REXEM(X(LEM),X(LDT),X(LFAM),X(LFBM),X(LH1),X(LEEXC),          !calculate E_L's
     *           EXENA,EXENB,EXENC,EN,X(LEDFA),X(LEDFB),
     *           L1,L2,NMIC,DBG)
      IF(DBG) THEN
         WRITE(IW,*) 'EL: ENERGY OF MICROSTATES'
         DO L=1,NMIC                                                     !microstates
            WRITE(IW,'(1X,I1,A,f15.5)') L,'TH MICROSTATE',X(LEM+L-1)
         ENDDO
         WRITE(IW,*) 'PREIVOUS INFO'
         WRITE(IW,*) 'nr/2, ns/2', DNR,DNS
         WRITE(IW,*) 'CL: WEIGHT FACTOR OF MICROSTATES'
         DO L=1,NMIC                                                     !microstates
            WRITE(IW,'(1X,I1,A,f15.5)') L,'TH MICROSTATE',X(LCM+L-1)
         ENDDO
      END IF
C
      THRESH=1.0D-8
      CNVLIM=1.0D-12
      MAXITREX=20
      IOUT=IW
      CALL REXSOLVER(X(LEM),DELTA,THRESH,CNVLIM,DNR,EREKS,MAXITREX,      !update FONs
     *               IOUT,CONST,G1,DBG)
C
      if(abs(dnr).lt.0.5d0.and..not.DODIIS)then                          !swap the active orbitals, if dnr<dns; don't do that, when DIIS is on; need to figure out why
         if(MASPRT.AND.(NFG.EQ.0.OR.NPRINT.NE.-5))                       !write, if master proc
     *   write(iw,'(1x,"DNR =",F14.9,"<0.5; Swap the active orbitals")')
     *   dnr
         DNR = one - DNR
         do kk = 1,l1
            temp = x(lva+l1*(na-1)+kk-1)
            x(lva+l1*(na-1)+kk-1) = x(lva+l1*na+kk-1)
            x(lva+l1*na+kk-1) = temp
         enddo
         temp = x(lea+l1*(na-1)+kk-1)
         x(lea+l1*(na-1)+kk-1) = x(lea+l1*na+kk-1)
         x(lea+l1*na+kk-1) = temp
C
         temp=x(lem)                                                     !swap energies of the 1st and 2nd microstates
         x(lem)=x(lem+1)
         x(lem+1)=temp
C
         do kk=1,L2
            temp=x(lfam + kk - 1)                                        !swap the 1st and the 2nd α-Fock matrices
            x(lfam + kk - 1) = x(lfam + L2 + kk -1)
            x(lfam + L2 + kk -1) = temp
            temp=x(lfbm + kk - 1)                                        !swap the 1st and the 2nd β-Fock matrices
            x(lfbm + kk - 1) = x(lfbm + L2 + kk -1)
            x(lfbm + L2 + kk -1) = temp
         enddo
C
c        IF(DODIIS) THEN                                                 !reset DIIS, if orbitals swapped
c          ITDIIS=1
c          IF(.NOT.NOTOPN) CALL RACLOS(20,'DELETE')
c          NOTOPN=.TRUE.
c        END IF
      endif                                                              !if(abs(dnr).lt.0.5d0.and..not.DODIIS)then
C
      DNS=ONE-DNR
      CALL REXCM(X(LCM),DNR,DNS,DELTA,WPPS,WOSS,NMIC)                    !update C_L's
C
      IF(DBG) THEN
         WRITE(IW,*) 'UPDATED INFO'
         WRITE(IW,*) 'nr/2, ns/2', DNR,DNS
         WRITE(IW,*) 'CL: WEIGHT FACTOR OF MICROSTATES'
         DO L=1,NMIC                                                     !microstates
            WRITE(IW,'(1X,I1,A,f15.5)') L,'TH MICROSTATE',X(LCM+L-1)
         ENDDO
      END IF
      CALL VCLR(X(LFA),1,L2)                                             !clear the array for REKS Fock matrix
      ETOTtmp = ZERO
      FR=DNR*WPPS+HALF*WOSS                                              !compute the active orbitals FONs; for SA-REKS too
      FS=DNS*WPPS+HALF*WOSS
C
      DO L=1,NMIC                                                        !transform Fock matrices from AO to MO basis
         LFAMTMP  =LFAM  + (L-1)*L2
         LFAMMOTMP=LFAMMO+ (L-1)*L3
         CALL CANTOGEN(X(LFAMTMP),X(LWRK1),L1)
         CALL DGEMM('T','N',L1,L1,L1,ONE,X(LVA),L1,X(LWRK1),L1,ZERO,
     *              X(LWRK2),L1)
         CALL DGEMM('N','N',L1,L1,L1,ONE,X(LWRK2),L1,X(LVA),L1,ZERO,
     *              X(LFAMMOTMP),L1)

         LFBMTMP  =LFBM  + (L-1)*L2
         LFBMMOTMP=LFBMMO+ (L-1)*L3
         CALL CANTOGEN(X(LFBMTMP),X(LWRK1),L1)
         CALL DGEMM('T','N',L1,L1,L1,ONE,X(LVB),L1,X(LWRK1),L1,ZERO,
     *              X(LWRK2),L1)
         CALL DGEMM('N','N',L1,L1,L1,ONE,X(LWRK2),L1,X(LVB),L1,ZERO,
     *              X(LFBMMOTMP),L1)
      ENDDO

      IF(DBG) THEN
         WRITE(IW,*) 'FOCK in MO'
         DO L=1,NMIC                                                     !microstates
            LFAMTMP=LFAMMO + (L-1)*L3
            LFBMTMP=LFBMMO + (L-1)*L3
            WRITE(IW,'(1X,I1,A)') L,'TH MICROSTATE'
            WRITE(IW,*) 'ALPHA FOCK MATRIX'
            CALL PRSQL(X(LFAMTMP),L1,L1,L1)
            WRITE(IW,*) 'BETA FOCK MATRIX'
            CALL PRSQL(X(LFBMTMP),L1,L1,L1)
         ENDDO
      END IF
C
      WRS = 0.d0                                                         !set r-s Lagrangian to zero
      CALL VCLR(X(LWRK3),1,L3)                                           !clear the array for SA-REKS Lagrangian
      DO L=1,NMIC                                                        !build REKS Fock matrix
         LFAMTMP=LFAMMO + (L-1)*L3
         LFBMTMP=LFBMMO + (L-1)*L3
         CALL GENTOCAN(X(LFAMTMP),X(LWRK1),L1)
         CALL GENTOCAN(X(LFBMTMP),X(LWRK2),L1)
         CMTMP =X(LCM+L-1)
         IF(L.GE.3) CMTMP=TWO*CMTMP
         CALL REXFM2FE(L,CMTMP,FR,FS,L1,L2,L3,NA-1,NA+2,
     *                 X(LWRK1),X(LWRK2),X(LFA),X(LWRK3),WRS,IA)         !SA-REKS Lagrangian (lower triangular) goes to -LWRK3-
      ENDDO
      CALL DAWRIT(IDAF,IODA,X(LWRK3),L2,36,0)                            !save SA-REKS Lagrangian to disk (in MO rep.)
C
      OFFFOCK = REXOFFFOC(X(LFA),L1,NA-1)                                !FIND MAX OFF-DIAGONAL ELEMENT OF THE REKS FOCK MATRIX
C
      ETOT0= ETOT
      ETOT = REXEM2EE(X(LCM),X(LEM),NMIC)                                !compute REKS total energy
C
      IF(ITER.EQ.1) E0NBO = ETOT
      DELE0= DELE
      DELE = ETOT-ETOT0
      IF(ITER.EQ.1) DEAVG = ZERO
      IF(ITER.EQ.2) DEAVG = ABS(DELE)
      IF(ITER.GE.3) DEAVG = (ABS(DELE)+ABS(DELE0)+PT2*DEAVG)/TWOPT2
C
      CALL DAWRIT(IDAF,IODA,X(LFA),L2,14,0)                              !save REKS Fock matrix to disk
      CALL DAWRIT(IDAF,IODA,X(LFA),L2,18,0)                              !need to make bogus write down to avoid further error messages
c
      CALL TSECND(TITER1)
C
      IF(DBG) THEN
         CALL CANTOGEN(X(LFA),X(LWRK1),L1)
         WRITE(IW,*) 'ALPHA REKS FOCK MATRIX'
         CALL PRSQL(X(LWRK1),L1,L1,L1)
      END IF
C
      IF(.NOT.(CVGED  .AND.  ICCFLG.GT.0)) THEN                          !go to line 1518
C
C     -LWRK1,LWRK2,LWRK3,LWRK4- ARE USED AS SQUARE STORAGE.
      IF(DODIIS) THEN                                                    !PERFORM DIIS INTERPOLATION
        CALL EXPND(X(LFA),X(LWRK3),L1,1)                                 !expand lower triangular M on the full square M
        CALL REXDIISER(X(LWRK3),FR,FS,L1,NA-1,errval)                    !do the REKS DIIS error vector in LWRK3

C convert DIIS error to AO basis: S*C*F(MO)*C^{\dagger}*S
      CALL DGEMM('n','N',L1,L1,L1,ONE,X(LVA),L1,X(LWRK3),L1,ZERO,        !take the product C*F
     *           X(LWRK1),L1)
      CALL TRPOSQ(X(LWRK1),L1)                                           !and transpose it to F*C^{\dagger}
      CALL DGEMM('n','N',L1,L1,L1,ONE,X(LVA),L1,X(LWRK1),L1,ZERO,        !multiply by C: C*(F*C^{\dagger})
     *           X(LWRK3),L1)
      IF(MINMEM) CALL DAREAD(IDAF,IODA,X(LS),L2,12,0)
      CALL EXPND(X(LS),X(LWRK2),L1,0)
      CALL DGEMM('n','N',L1,L1,L1,ONE,X(LWRK2),L1,X(LWRK3),L1,ZERO,      !multiply by S; first on the left
     *           X(LWRK1),L1)
      CALL DGEMM('n','N',L1,L1,L1,ONE,X(LWRK1),L1,X(LWRK2),L1,ZERO,      !now, on the right
     *           X(LWRK3),L1)
C convert Fock matrix to AO basis; this should be done, because the MO phase is arbitrary; in AO everything is constant
      CALL EXPND(X(LFA),X(LWRK1),L1,0)                                   !unpack F.lower onto F.square
      CALL DGEMM('n','N',L1,L1,L1,ONE,X(LVA),L1,X(LWRK1),L1,ZERO,
     *           X(LWRK2),L1)
      CALL TRPOSQ(X(LWRK2),L1)
      CALL DGEMM('n','N',L1,L1,L1,ONE,X(LVA),L1,X(LWRK2),L1,ZERO,
     *           X(LWRK1),L1)
      IF(MINMEM) CALL DAREAD(IDAF,IODA,X(LS),L2,12,0)
      CALL EXPND(X(LS),X(LWRK2),L1,0)
      CALL DGEMM('n','N',L1,L1,L1,ONE,X(LWRK2),L1,X(LWRK1),L1,ZERO,
     *           X(LWRK4),L1)
      CALL DGEMM('n','N',L1,L1,L1,ONE,X(LWRK4),L1,X(LWRK2),L1,ZERO,
     *           X(LWRK1),L1)
C
      CALL CPYSQT(X(LWRK1),X(LFA),L1,1)                                  !pack F.square back to F.lower
cclsh test
c        write(iw,*) 'Fock in MO'
c        CALL PRSQ(X(LWRK3),L1,L1,L1)

C
      IF(MINMEM) CALL DAREAD(IDAF,IODA,X(LQ),L3,45,0)
      CALL DIIS(SCFTYP,IW,ITDIIS,X(LQ),X(LFA),X(LFB),X(LWRK3),
     *          X(LWRK4),X(LWRK1),X(LADIIS),X(LXDIIS),X(LPDIIS),
     *          X(LBDIIS),X(LIODII),X(LSCR),L1,L2,L3,MAXIT,MAXIT2,
     *          4*MAXDII,ERDIIS,NOTOPN,MASWRK)
C convert Fock back to MO: C^{\dagger}*F(AO)*C
      CALL EXPND(X(LFA),X(LWRK1),L1,0)
      CALL DGEMM('t','N',L1,L1,L1,ONE,X(LVA),L1,X(LWRK1),L1,ZERO,
     *           X(LWRK2),L1)
      CALL DGEMM('n','N',L1,L1,L1,ONE,X(LWRK2),L1,X(LVA),L1,ZERO,
     *           X(LWRK1),L1)
      CALL CPYSQT(X(LWRK1),X(LFA),L1,1)
C
      END IF                                                             !IF(DODIIS) THEN
C
      do k = na,l1                                                       !apply level shift
       x(lfa-1+k*(k+1)/2) = x(lfa-1+k*(k+1)/2) + vshftVal
      enddo
      do k = na+1,l1
       x(lfa-1+k*(k+1)/2) = x(lfa-1+k*(k+1)/2) + vshftVal
      enddo
      do k = na+2,l1
       x(lfa-1+k*(k+1)/2) = x(lfa-1+k*(k+1)/2) + vshftVal
      enddo
C
      IF(DBG) THEN
         write(iw,*) 'vector test'
         WRITE(IW,*) 'previous MO'
         CALL PREVS(X(LVA),X(LEA),X(LIWRK),L0,L1,L1,1)
      END IF
C
      CALL EXPND(X(LFA),X(LWRK1),L1,0)                                   !copy Fock matrix to full square matrix
      info = 0
      call dsyev('V','L',L1,X(LWRK1),L1,X(LEA),X(LSCR),1+6*L1+3*L3,info) !diagonalize
      if (Info.ne.0) then
        write(iw,'(1x,"Error in DSYEV in REKS; INFO =",i10)')info
        call ABRT
      endif
      do j = 1,L1                                                        !fixing orbital phase before multiplying by the old vectors
       jl1 = (j - 1) * L1
        if(X(LWRK1 + jl1 + j - 1).lt.0.d0)                               !only diagonal elements are checked; the eigenvectors converge
     *  call dscal(L1,-1.d0,X(LWRK1 + jl1),1)                            !to unit matrix anyway
      enddo
      CALL DCOPY(L3,X(LVA),1,X(LWRK2),1)                                 !copy old vectors in LWRK2
      CALL DGEMM('N','N',L1,L1,L1,ONE,X(LWRK2),L1,X(LWRK1),L1,ZERO,      !multiply the new and old vectors
     *           X(LVA),L1)

      do k = na,l1                                                       !remove the level shift; from the orbital energies
       x(lea-1+k) = x(lea-1+k) - vshftVal
      enddo
      do k = na+1,l1
       x(lea-1+k) = x(lea-1+k) - vshftVal
      enddo
      do k = na+2,l1
       x(lea-1+k) = x(lea-1+k) - vshftVal
      enddo
C
      if(rexLdL.eq.1)then
      call localz(x(lwrk1),x(lva+l1*(na-1)),x(lva+l1*na),l1,iw)          !localize a pair of fractionally occupied orbitals
      elseif(rexLdL.eq.2)then
      call delclz(x(lwrk1),x(lva+l1*(na-1)),x(lva+l1*na),l1,iw)          !de-localize a pair of fractionally occupied orbitals
      endif

C
      IF(DBG) THEN
         WRITE(IW,*) 'current MO'
         CALL PREVS(X(LVA),X(LEA),X(LIWRK),L0,L1,L1,1)

         write(iw,*) 'vector test'
         diffatmp=0.0
         DO I=1,L1
         DO J=1,L1
            ijl=(i-1)*l3+j
            tmp=abs(x(lva+ijl-1)-x(lwrk2+ijl-1))
            diffatmp=diffatmp+tmp
            write(iw,'(2i3,3f10.5)')
     *         i,j,x(lva+ijl-1),x(lwrk2+ijl-1),tmp
         ENDDO
         ENDDO
         write(iw,*) 'diffav',diffatmp
         WRITE(IW,*) 'NEW ALPHA ORBITALS'
         CALL PREV(X(LVA),X(LEA),L0,L1,L1)
      END IF
C
      CALL DCOPY(L2*NOCP,X(LDT),1,X(LWRK1),1)                            !copy previous density matrix to -lwrk1-
      CALL REXDMAT(X(LDT),X(LVA),IA,L1,L2,L3,NA)                         !build new density matrix
      CALL DDIFF(X(LWRK1),X(LDT),L2*NOCP,DIFFA)                          !compute density matrices difference

      IF(DBG) THEN
      DO L=1,4
      diffatmp=0.0
      DO I=1,L2
         ijl=(l-1)*l2+i
         tmp=abs(x(ldt+ijl-1)-x(lwrk1+ijl-1))
         diffatmp=diffatmp+tmp
         write(iw,'(2i3,3f10.5)')
     *   l,i,x(ldt+ijl-1),x(lwrk1+ijl-1),tmp
      ENDDO
      write(iw,*) 'diffa1',diffatmp
      diffatmp=0.0

      DO I=1,L2
         ijl=(l-1)*l2+i
         tmp=abs(x(ldt+ijl-1))-abs(x(lwrk1+ijl-1))
         diffatmp=diffatmp+tmp
         write(iw,'(2i3,3f10.5)')
     *   l,i,x(ldt+ijl-1),x(lwrk1+ijl-1),tmp
      ENDDO
      write(iw,*) 'diffa2',diffatmp
      ENDDO

         WRITE(IW,*) 'NEW DENSITY MATRIX BASED ON OCCUPATION NUMBERS'
         DO L=1,NOCP
            LDTTMP=LDT + (L-1)*L2
            WRITE(IW,'(1X,I1,A)') L,'TH OCCUP'
            WRITE(IW,*) 'DENSITY MAT'
            CALL PRTRIL(X(LDTTMP),L1)
         ENDDO
         WRITE(IW,*) '*** END OF DEBUG OUTPUT FOR ITERATION',ITER
      END IF

      DIFFB=ZERO
      DIFF  = DIFFA+DIFFB
C
      CALL TSECND(TITER2)
      TFOCK  = TFOCK  + (TITER1-TITER0)
      TSOLVE = TSOLVE + (TITER2-TITER1)
      IF(ITER.EQ.1) TFOCK1 = TITER1 - TITER0
      TFOCKN = TITER1 - TITER0
C
      IF(MASPRT) THEN                                                    !PRINT CURRENT ITERATION'S RESULTS
            ERRORC = ERDIIS
         IF(DIRSCF) THEN
           WRITE(IW,9070) ITER,ICAB,DNR,ETOT,OFFFOCK,DELE,DIFF,
     *                    ERRORC,vshftVal,NINT,NSCHWZ
         ELSE
           WRITE(IW,9070) ITER,ICAB,DNR,ETOT,OFFFOCK,DELE,DIFF,
     *                    ERRORC,vshftVal
         END IF
      CALL FLSHBF(IW)
      END IF
C
      ICALP = ICALP+1
      ICBET = ICBET+1
      ICAB  = ICAB+1
C
      CVDENS = (DIFF.LT.DENTOL)                                          !convergence on density matrix
      CVFOCK = (OFFFOCK.LT.OFFTOL)                                       !convergence on fock matrix off-diag element
      IF(DODIIS) CVDIIS = ERDIIS.LT.DIITOL  .AND.  DIFF.LT.DTOL2         !convergence on DIIS error
      CVGED  = CVGING  .AND.  (CVDENS.OR.CVFOCK.OR.CVDIIS)
      CVGED  = CVGED  .AND.  ISWOFF.EQ.0
      IF(NDFTFG.EQ.1) CVGED = CVGED  .AND.  GRIDTYP.EQ.FINE
      CVGING = CVDENS.OR.CVFOCK.OR.CVDIIS
C
      IF(ISVP.EQ.1) THEN
         CVGED  = CVGED   .AND.  CVGSVP
         CVGING = CVGING  .AND.  CVGSVP
         CVGSVPNR = DIFF  .LT.  100.0D+00 * DENTOL
         CVGSVPDN = CVGING
         CVGSVPDN1 = CVGING
      END IF
C
      IDFTSWITCH_SAVE=IDFTSWITCH
      IF(ISWOFF.GT.0.AND.DIFF.LT.SWOFF) THEN                             !TURN ON THE DFT IF PRELIMINARY SCF IS BECOMING CONVERGED
        IF(MASWRK) WRITE(IW,9355)
        DFTTYP(3)=HFSSAV
        NDFTFG=NDFTSV
        LCFLAG=LCFLAGS
        LRINT=LRINTS
        ISWOFF=0
        IDFTSWITCH=1
      END IF
C
      IF(NDFTFG.EQ.1.AND.DIFF.LT.SWDFT  .AND.  SWGROK) THEN              !DFT MUST SWITCH TO A TIGHTER GRID AS IT NEARS CONVERGENCE
         IF(NRAD*NANGPT(1).LT.NRAD0*NANGPT0(1)  .AND.
     *       .NOT.(SG1.OR.JANS.GT.0)) THEN
            IF(MASPRT) WRITE(IW,9330)
            CALL SWGRID
            GRIDTYP = FINE
            CALL DFTSET(X(IGRDDFT),0,.FALSE.)
            IDFTSWITCH=2
         END IF
      END IF
C
      IF(IDFTSWITCH.NE.0) THEN                                           !ANY SWITCH IN GRIDS OR CONVERGERS MUST RESET CONVERGERS, AND PERHAPS RESET THE DIFFERENTIAL FOCK FORMATION.
         if(.not.locopt) then
            EXTRAH = .FALSE.
            DAMPH  = .FALSE.
            VSHIFT = .FALSE.
            endif
         IF(DODIIS) THEN
           ITDIIS=1
           IF(.NOT.NOTOPN) CALL RACLOS(20,'DELETE')
           NOTOPN=.TRUE.
         END IF
         IF(ABINIT .AND. DIRSCF .AND. FDIFF .AND. IDFTSWITCH.NE.4) THEN
           IF(MINMEM) CALL DAREAD(IDAF,IODA,X(LH1),L2,11,0)
           CALL DAWRIT(IDAF,IODA,X(LH1),L2,14,0)
           CALL DAWRIT(IDAF,IODA,X(LH1),L2,18,0)
           CALL VCLR(X(LWRK1),1,L2)
           CALL VCLR(X(LWRK2),1,L2)
         END IF
         IDFTSWITCH=0
      END IF
C
      IF(CVGED  .AND.  CCTYP.NE.ANONE) THEN
         ICCFLG=ICCFLG+1
         IF(MASPRT) WRITE(IW,9105)
         IF(ABINIT  .AND.  DIRSCF  .AND.  FDIFF) THEN
            IF(MINMEM) CALL DAREAD(IDAF,IODA,X(LH1),L2,11,0)
            CALL DAWRIT(IDAF,IODA,X(LH1),L2,14,0)
            CALL DAWRIT(IDAF,IODA,X(LH1),L2,18,0)
            CALL VCLR(X(LWRK1),1,L2)
            CALL VCLR(X(LWRK2),1,L2)
         END IF
         GO TO 310
      END IF
      END IF                                                             !IF(.NOT.(CVGED  .AND.  ICCFLG.GT.0))
C
      IF(CVGED .AND. CVDENS) THEN
         IF(MASPRT) WRITE(IW,9080)
         GO TO 600
      END IF
      IF(CVGED .AND. CVFOCK) THEN
         IF(MASPRT) WRITE(IW,9090)
         GO TO 600
      END IF
      IF(CVGED .AND. CVDIIS) THEN
         IF(MASPRT) WRITE(IW,9100)
         GO TO 600
      END IF
C
      CALL TSECND(TIM1)                                                  !EXIT IN CASE OF TIME LIMIT
      TLEFT = TIMLIM - TIM1
      TIM0 = TIM1
      IF(TLEFT .GT. 2.0D+00*(TLEFTS-TLEFT)/ITER) CYCLE
      IF(CVGING) THEN
         IF(MASPRT) WRITE(IW,9110)
         GO TO 600
      ELSE
         IF(MASPRT) WRITE(IW,9120)
         GO TO 550
      END IF
      END DO
C           **********************
C     ***** END OF REKS ITERATIONS *****
C           **********************
      IF(MASPRT) WRITE(IW,9130)                                          !FALLING OUT OF LOOP 500 MEANS CONVERGENCE FAILURE
      ITER = MAXIT
C
  550 CONTINUE
      ETOT = ZERO                                                        !SET ENERGY TO ZERO IF FAILED TO CONVERGE
      EHF = -EN
C
  600 CONTINUE                                                           !BRANCH TO HERE ON SUCCESFUL CONVERGENCE
C
      call dawrit(idaf,ioda,x(lfam),l2*nmic,700,0)
      call dawrit(idaf,ioda,x(lfbm),l2*nmic,701,0)
      call dawrit(idaf,ioda,x(lem),nmic,702,0)
      call dawrit(idaf,ioda,x(lcm),nmic,703,0)
C
      TFITER = TFOCK/ITER
      TSITER = TSOLVE/ITER
      IF(MASPRT) THEN
          IF(DIRSCF) THEN
             WRITE(IW,9400) TFOCK,TFITER,TFOCK1,TFOCKN,TSOLVE,TSITER
          ELSE
             WRITE(IW,9410) TFOCK,TFITER,TSOLVE,TSITER
          END IF
      END IF
C
      IF (NDFTFG.EQ.1.AND.ABINIT.AND.DIRSCF.AND.FDIFF) THEN              !SAVE GRID DFT FOCK MATRICES
            CALL DAWRIT(IDAF,IODA,X(LFA),L2,14,0)
            CALL DAWRIT(IDAF,IODA,X(LFB),L2,18,0)
      END IF
C
      IF(NFG.NE.0.AND.NDFTFG.EQ.1.AND.                                   !IF SCF IS UNCONVERGED, GRID MAY HAVE BEEN LEFT IN REVERSE ORDER.
     *     (NRAD*NANGPT(1).LT.NRAD0*NANGPT0(1)  .AND.                    !WE MAY LET SCF BE UNCONVERGED HOPING IT WILL CONVERGE IN THE NEXT MONOMER SCF ITERATION.
     *       .NOT.(SG1.OR.JANS.GT.0))) THEN                              !NOW THE PROBLEM WITH REVERSE ORDER IS THAT THE NEXT DFT RUN WILL NOT RUN! (MEMORY INADEQUATE)!
         CALL SWGRID
         GRIDTYP = FINE
      END IF
C
      IF(DC) THEN                                                        !DISPERSION CORRECTIONS
         CALL DFTD3(1,EDISP,DUMMY)
         ETOTD = ETOT
         IF(ETOT.NE.ZERO) ETOT = ETOT + EDISP
         E0NBO = E0NBO + EDISP
      ELSE
         EDISP=0.0D+00
      END IF
C
      IF(LRDFLG) THEN                                                    !LOCAL RESPONSE DISPERSION CORRECTIONS
         EDISP = ELRD6 + ELRD8 + ELRD10 + EMULT
         ETOTD = ETOT
         IF(ETOT.NE.ZERO) ETOT  = ETOT + EDISP
         CALL PRTLRD
         E0NBO = E0NBO + EDISP
      END IF
C
      IF(MASPRT) THEN                                                    !PRINT FINAL RESULTS
      write(iw,'(/1x,a6,a,a10,F20.10,a6,I4,a11)')'FINAL ',
     *     line(:lenline),' ENERGY IS',ETOT,' AFTER',ITER,' ITERATIONS'
            IF(DFTTYP(1) .NE. 0.0D+00) THEN
               WRITE(IW,9210) EXENA
               WRITE(IW,9220) EXENB
            END IF
            IF(DC  .AND.  MASPRT) THEN                                   !PRINT FOR DISPERSION CORRECTIONS
               WRITE(IW,9510) EDISP
               WRITE(IW,9500) ETOTD
            END IF
C
            IF(LRDFLG) THEN
               WRITE(IW,9550) ETOTD
            END IF
C
      SZ = (NA-NB)/2.0                                                   !SPIN EXPECTATION VALUE
      S2 = SZ*(SZ+1.0)
c     WRITE(IW,9148) SZ,S2
      write(iw,'(1x,a,":  FON(",i3,") =",F9.6," FON(",i3,") =",F9.6)')
     *line(:lenline),NA,DNR+DNR,NA+1,DNS+DNS
C
      if(dirscf) then
       write(iw,'(/1x,a)')ddashes(:132)
      else
       write(iw,'(/1x,a)')ddashes(:116)
      endif
      END IF
C
      if(abs(WPPS-1.d0).LE.1.d-6)then                                    !skip the following section, if single state REKS
C compute the total density matrix for single state REKS
      call vclr(x(lwrk1),1,l3)                                           !clear -LWRK1- array; this will be for the total density matrix
      call dgemm('n','t',l1,l1,na-1,2.d0,x(lva),l1,x(lva),l1,0.d0,       !compute density matrix
     *          x(lwrk1),l1)                                             !take product C_{core}*C_{core}^{\dagger}; this is the core density matrix
      call dger(l1,l1,2.d0*fr,x(lva+l1*(na-1)),1,x(lva+l1*(na-1)),1,
     *          x(lwrk1),l1)                                             !take product afon*C_r*C_r^{\dagger} and add to -LDA-
      call dger(l1,l1,2.d0*fs,x(lva+l1*na),1,x(lva+l1*na),1,
     *          x(lwrk1),l1)                                             !take product bfon*C_s*C_s^{\dagger} and add to -LDA-
      call cpysqt(x(lwrk1),x(lda),L1,1)                                  !save the density matrix to lower triangular matrix -LDA-
        goto 650
      endif
C SA/SSR final printout; do transition dipole here
      rexprefix = ' SA:'
      if(rexType.ne.0) rexprefix = 'SSR:'                                !set prefix for SA/SSR printing

      if(masprt)
     *write(iw,'(/1x,a,1x,a)')rexprefix,'final printout'
      CALL REXCM(X(LCM),DNR,DNS,DELTA,1.d0,0.d0,NMIC)                    !compute new C_L's for pure PPS state
      ETRP = X(LEM+3)
      EPPS = X(LCM)*X(LEM) + X(LCM+1)*X(LEM+1) + 2.d0*X(LCM+2)*X(LEM+2)
     *     + 2.d0*X(LCM+3)*X(LEM+3)
      EDES = X(LCM+1)*X(LEM) + X(LCM)*X(LEM+1) - 2.d0*X(LCM+2)*X(LEM+2)
     *     - 2.d0*X(LCM+3)*X(LEM+3)
      EOSS = 2.d0*X(LEM+2) - X(LEM+3)
      lenprnt=lenline+43
      if(masprt)
     *write(iw,'(/1x,a)')ddashes(:lenprnt)
      line1='.'
      do kk=1,lenprnt-55
       line1=line1//'.'
      enddo
      if(masprt)
     *write(iw,'(4(/1x,a,a,F22.12)//1x,a,a,F22.12)')
     *'Triplet:.......................',line1,ETRP,
     *'Doubly excited singlet (DES):..',line1,EDES,
     *'Open Shell Singlet (OSS):......',line1,EOSS,
     *'Perfectly Paired Singlet (PPS):',line1,EPPS,
     *'Lagrangian Wrs:................',line1,WRS
C do 2SI-2SA-REKS energies
      if(rexType.eq.1) then
      if(masprt)
     *write(iw,'(1x,a)')ddashes(:66)
      s2sare(1,1)=EPPS
      s2sare(2,2)=EOSS
      s2sare(1,2)=WRS*(sqrt(dnr)-sqrt(dns))*sqrt(2.d0)
      s2sare(2,1)=s2sare(1,2)
      infosi = 0
      call dsyev('V','L',2,s2sare,2,ws2sa,tms2sa,5,infosi)
      if (Infosi.ne.0) then
        write(iw,'(1x,"Error in DSYEV in REKS; INFO =",i10)')infosi
        call ABRT
      endif
      if(masprt)
     *write(iw,'(1x,a,a,a,/23x,"E_k",15x,"C_{PPS}",9x,"C_{OSS}")')
     *'2SI-2',line(:lenline),' states:'
      if(masprt)
     *write(iw,'(1x,a)')dashes(:66)
      do kk=1,2
      if(masprt)
     *write(iw,'(1x,"SSR state",i2,F22.12,2F16.8)')
     *kk-1,ws2sa(kk),(s2sare(jj,kk),jj=1,2)
      enddo
      if(masprt)
     *write(iw,'(1x,a)')dashes(:66)
      elseif(rexType.eq.2) then
C do 3SI-2SA-REKS energies
      if(masprt)
     *write(iw,'(1x,a)')ddashes(:82)
      s3sare(1,1)=EPPS
      s3sare(2,2)=EOSS
      s3sare(3,3)=EDES
      s3sare(1,2)=WRS*(sqrt(dnr)-sqrt(dns))*sqrt(2.d0)
      s3sare(2,1)=s3sare(1,2)
      s3sare(2,3)=WRS*(sqrt(dnr)+sqrt(dns))*sqrt(2.d0)
      s3sare(3,2)=s3sare(2,3)
      s3sare(1,3)=0.d0
      s3sare(3,1)=0.d0
      infosi = 0
      call dsyev('V','L',3,s3sare,3,ws3sa,tms3sa,9,infosi)
      if (Infosi.ne.0) then
        write(iw,'(1x,"Error in DSYEV in REKS; INFO =",i10)')infosi
        call ABRT
      endif
      if(masprt)
     *write(iw,'(1x,a,a,a,/23x,"E_k",15x,"C_{PPS}",9x,"C_{OSS}",
     *      9x,"C_{DES}")')
     *'3SI-2',line(:lenline),' states:'
      if(masprt)
     *write(iw,'(1x,a)')dashes(:82)
      do kk=1,3
      if(masprt)
     *write(iw,'(1x,"SSR state",i2,F22.12,3F16.8)')
     *kk-1,ws3sa(kk),(s3sare(jj,kk),jj=1,3)
      enddo
      if(masprt)
     *write(iw,'(1x,a)')dashes(:82)
      endif                                                              !if..elseif(rexType.eq.2)
C compute FONs and transition densities
      if(rexTarget.ne.0) then
      if(rexTarget.eq.1) then
      if(rexType.eq.1) then
       rexenergy = ws2sa(1)                                              !set the first eigenvalue as the reks energy
       afon = dnr*s2sare(1,1)**2 + 0.5d0*s2sare(2,1)**2                  !CL(1)*a_11^2 + 0.5*a_21^2
       bfon = dns*s2sare(1,1)**2 + 0.5d0*s2sare(2,1)**2                  !CL(2)*a_11^2 + 0.5*a_21^2
       cftran = s2sare(1,1)*s2sare(2,1)*
     *         (sqrt(0.5d0*dnr)-sqrt(0.5d0*dns))                         !0.5*(sqrt(nr)-sqrt(ns))*a_11*a_21
      elseif(rexType.eq.2) then
       rexenergy = ws3sa(1)                                              !set the first eigenvalue as the reks energy
       afon = dnr*s3sare(1,1)**2 + 0.5d0*s3sare(2,1)**2                  !CL(1)*a_11^2 + CL(2)*a_31^2 + 0.5*a_21^2 + a_11*a_31*sqrt(CL(1)*CL(2))
     *      + dns*s3sare(3,1)**2 + s3sare(1,1)*s3sare(3,1)*
     *        sqrt(dnr*dns)
       bfon = dns*s3sare(1,1)**2 + 0.5d0*s3sare(2,1)**2                  !CL(2)*a_11^2 + CL(1)*a_31^2 + 0.5*a_21^2 - a_11*a_31*sqrt(CL(1)*CL(2))
     *      + dnr*s3sare(3,1)**2 - s3sare(1,1)*s3sare(3,1)*
     *        sqrt(dnr*dns)
       cftran = s3sare(1,1)*s3sare(2,1)*
     *         (sqrt(0.5d0*dnr)-sqrt(0.5d0*dns))                         !0.5*(sqrt(nr)-sqrt(ns))*a_21*a_11 + 0.5*(sqrt(nr)+sqrt(ns))*a_21*a_31
     *        + s3sare(2,1)*s3sare(3,1)*
     *         (sqrt(0.5d0*dnr)+sqrt(0.5d0*dns))
      elseif(rexType.eq.0) then
       rexenergy = EPPS
       afon = dnr
       bfon = dns
       cftran = 0.d0
      endif                                                              !if..elseif(rexType.eq.2)
      if(masprt)
     * write(iw,'(/1x,a,1x,a)')rexprefix,'the ground state is reported'
      elseif(rexTarget.eq.2) then
      if(rexType.eq.1) then
       rexenergy = ws2sa(2)
       afon = dnr*s2sare(1,2)**2 + 0.5d0*s2sare(2,2)**2                  !CL(1)*a_12^2 + 0.5*a_22^2
       bfon = dns*s2sare(1,2)**2 + 0.5d0*s2sare(2,2)**2                  !CL(2)*a_12^2 + 0.5*a_22^2
       cftran = s2sare(1,2)*s2sare(2,2)*
     *          (sqrt(0.5d0*dnr)-sqrt(0.5d0*dns))                        !0.5*(sqrt(nr)-sqrt(ns))*a_12*a_22
      elseif(rexType.eq.2) then
       rexenergy = ws3sa(2)                                              !set the second eigenvalue as the reks energy
       afon = dnr*s3sare(1,2)**2 + 0.5d0*s3sare(2,2)**2                  !CL(1)*a_12^2 + CL(2)*a_32^2 + 0.5*a_22^2 + a_12*a_32*sqrt(CL(1)*CL(2))
     *      + dns*s3sare(3,2)**2 + s3sare(1,2)*s3sare(3,2)*
     *        sqrt(dnr*dns)
       bfon = dns*s3sare(1,2)**2 + 0.5d0*s3sare(2,2)**2                  !CL(2)*a_12^2 + CL(1)*a_32^2 + 0.5*a_22^2 - a_12*a_32*sqrt(CL(1)*CL(2))
     *      + dnr*s3sare(3,2)**2 - s3sare(1,2)*s3sare(3,2)*
     *        sqrt(dnr*dns)
       cftran = s3sare(1,2)*s3sare(2,2)*
     *         (sqrt(0.5d0*dnr)-sqrt(0.5d0*dns))                         !0.5*(sqrt(nr)-sqrt(ns))*a_22*a_12 + 0.5*(sqrt(nr)+sqrt(ns))*a_22*a_32
     *        + s3sare(2,2)*s3sare(3,2)*
     *         (sqrt(0.5d0*dnr)+sqrt(0.5d0*dns))
      elseif(rexType.eq.0) then
       rexenergy = EOSS
       afon = 0.5d0
       bfon = 0.5d0
       cftran = 0.d0
      endif                                                              !if..elseif(rexType.eq.2)
      if(masprt)
     * write(iw,'(/1x,a,1x,a)')rexprefix,
     *                       'the first excited state is reported'
      endif                                                              !if...elseif(rexTarget.eq.2)
       AplusB = afon + bfon                                              !find actual FONs by diagonalizing
       AminB2 = (afon - bfon)**2                                         ! |  afon  | cftran |
       det = sqrt(4.d0*cftran**2 + AminB2)                               ! | cftran |  bfon  |
       fona = AplusB + det
       fonb = AplusB - det
      if(masprt)
     * write(iw,'(1x,a,1x,a,
     *       ":  FON(",i3,") =",F9.6," FON(",i3,") =",F9.6)')
     *       rexprefix,'unrelaxed occupation numbers',NA,fona,NA+1,fonb
      elseif(rexTarget.eq.0) then                                        !if(rexTarget.ne.0); case: SA-REKS averaged state
      afon=fr
      bfon=fs
      cftran=0.d0
      rexenergy=etot
      if(masprt)
     *write(iw,'(/1x,a,1x,a)')rexprefix,'the averaged state is reported'
      endif                                                              !if(rexTarget.ne.0)
      if(masprt)
     *write(iw,'(/1x,a," final energy =",F22.12)')rexprefix,rexenergy    !print the REKS/SA/SSR energy
C compute density matrix of the target SA/SSR state
      call vclr(x(lwrk1),1,l3)                                           !clear -LWRK1- array; this will be for the total density matrix
      if(na.ne.1)                                                        !precaution: don't do this if there are no core electrons
     *call dgemm('n','t',l1,l1,na-1,2.d0,x(lva),l1,x(lva),l1,0.d0,
     *          x(lwrk1),l1)                                             !take product C_{core}*C_{core}^{\dagger}; this is the core density matrix
      call dcopy(l3,x(lwrk1),1,x(lwrk2),1)                               !copy -LWRK1- to -LWRK2- for further use
      call dger(l1,l1,2.d0*afon,x(lva+l1*(na-1)),1,x(lva+l1*(na-1)),1,
     *          x(lwrk1),l1)                                             !take product afon*C_r*C_r^{\dagger} and add to -LWRK1-
      call dger(l1,l1,2.d0*bfon,x(lva+l1*na),1,x(lva+l1*na),1,
     *          x(lwrk1),l1)                                             !take product bfon*C_s*C_s^{\dagger} and add to -LWRK1-
      call dger(l1,l1,+2.d0*cftran,x(lva+l1*(na-1)),1,x(lva+l1*na),1,    !check out the sign
     *          x(lwrk1),l1)                                             !take product cft*C_r*C_s^{\dagger} and add to -LWRK1-
      call dger(l1,l1,+2.d0*cftran,x(lva+l1*na),1,x(lva+l1*(na-1)),1,
     *          x(lwrk1),l1)                                             !take product cft*C_s*C_r^{\dagger} and add to -LWRK1-
      CALL CPYSQT(X(LWRK1),X(LDA),L1,1)                                  !save the density matrix to lower triangular matrix -LDA-
c     CALL CPYSQT(X(LWRK1),X(LDB),L1,1)                                  !save the density matrix to lower triangular matrix -LDB-
      call vclr(x(ldb),1,l2)                                             !clear β density matrix
c     CALL PRSQL(X(LWRK1),L1,L1,L1)

C compute transition dipole moment between the S0 and S1 states, SA or SSR
      call vclr(x(lwrk1),1,l3)                                           !clear -LWRK1- array; this will be for the transition density matrix
c     cftran = (sqrt(2.d0*dnr)-sqrt(2.d0*dns))                           !tran. dens. coefficient (\sqrt(n_r/2) - \sqrt(n_s/2))
      cftran = (sqrt(1.d0*dnr)-sqrt(1.d0*dns))*sqrt(0.5d0)               !tran. dens. coefficient (\sqrt(n_r) - \sqrt(n_s))/2
      call dger(l1,l1,1.d0,x(lva+l1*(na-1)),1,x(lva+l1*na),1,            !compute transition density matrix; to be multiplied by cftran later
     *          x(lwrk1),l1)                                             !take product cft*C_r*C_s^{\dagger} and put to -LWRK1-
      call dger(l1,l1,1.d0,x(lva+l1*na),1,x(lva+l1*(na-1)),1,
     *          x(lwrk1),l1)                                             !take product cft*C_s*C_r^{\dagger} and add to -LWRK1-
C compute density matrices of PPS, OSS, and (possibly) DES states
      if(rexType.ne.0)then                                               !do density matrices of PPS, OSS, and DES in -LWRK2-, -LWRK3, and -LWRK4-
      call dcopy(l3,x(lwrk2),1,x(lwrk3),1)                               !core density in -LWRK2-; copy to -LWRK3- for OSS
      if(rexType.eq.2) call dcopy(l3,x(lwrk2),1,x(lwrk4),1)              !if need DES, copy to -LWRK4-
      call dger(l1,l1,2.d0*dnr,x(lva+l1*(na-1)),1,x(lva+l1*(na-1)),1,    !build PPS density matrix in -LWRK2-
     *          x(lwrk2),l1)                                             !take product afon*C_r*C_r^{\dagger} and add to -LWRK2-
      call dger(l1,l1,2.d0*dns,x(lva+l1*na),1,x(lva+l1*na),1,
     *          x(lwrk2),l1)                                             !take product bfon*C_s*C_s^{\dagger} and add to -LWRK2-
      call dger(l1,l1,1.d0,x(lva+l1*(na-1)),1,x(lva+l1*(na-1)),1,        !build OSS density matrix in -LWRK3-
     *          x(lwrk3),l1)                                             !take product afon*C_r*C_r^{\dagger} and add to -LWRK3-
      call dger(l1,l1,1.d0,x(lva+l1*na),1,x(lva+l1*na),1,
     *          x(lwrk3),l1)                                             !take product bfon*C_s*C_s^{\dagger} and add to -LWRK3-
      if(rexType.eq.2) then
      call dger(l1,l1,2.d0*dns,x(lva+l1*(na-1)),1,x(lva+l1*(na-1)),1,    !build DES density matrix in -LWRK4-
     *          x(lwrk4),l1)                                             !take product afon*C_r*C_r^{\dagger} and add to -LWRK4-
      call dger(l1,l1,2.d0*dnr,x(lva+l1*na),1,x(lva+l1*na),1,
     *          x(lwrk4),l1)                                             !take product bfon*C_s*C_s^{\dagger} and add to -LWRK4-
      endif                                                              !if(rexType.eq.2)
      endif                                                              !if(rexType.ne.0)then
C
      call calcom(xcom,ycom,zcom)                                        !calculate center of mass coordinates
      call dipint(xcom,ycom,zcom,.false.)                                !calculate the dipole integrals w.r.t. C.o.M.
      CALL DAREAD(IDAF,IODA,X(LFBM),L2,95,0)                             !read in the X dipole integrals; to -LFBM-, which is free
      CALL EXPND(X(LFBM),X(LFBMMO),L1,0)                                 !and expand it to -LFBMMO- to a full square form
      trdipX = ddot(l3,x(lwrk1),1,x(lfbmmo),1)                           !take trace of the product D^{tran}*Dip_x
      if(rexType.ne.0)then                                               !do electronic dipole moments of PPS, OSS, abd (possibly) DES states
       dipPPSX = ddot(l3,x(lwrk2),1,x(lfbmmo),1)*0.5d0                   !divide by 2; the core den mat was done for 2e occupations
       dipOSSX = ddot(l3,x(lwrk3),1,x(lfbmmo),1)*0.5d0
      if(rexType.eq.2)
     * dipDESX = ddot(l3,x(lwrk4),1,x(lfbmmo),1)*0.5d0
      endif                                                              !if(rexType.ne.0)then
      CALL DAREAD(IDAF,IODA,X(LFBM),L2,96,0)                             !read in the Y dipole integrals to -LFBM-
      CALL EXPND(X(LFBM),X(LFBMMO),L1,0)                                 !and expand it to -LFBMMO- to a full square form
      trdipY = ddot(l3,x(lwrk1),1,x(lfbmmo),1)                           !take trace of the product D^{tran}*Dip_y
      if(rexType.ne.0)then                                               !do electronic dipole moments of PPS, OSS, abd (possibly) DES states
       dipPPSY = ddot(l3,x(lwrk2),1,x(lfbmmo),1)*0.5d0
       dipOSSY = ddot(l3,x(lwrk3),1,x(lfbmmo),1)*0.5d0
      if(rexType.eq.2)
     * dipDESY = ddot(l3,x(lwrk4),1,x(lfbmmo),1)*0.5d0
      endif                                                              !if(rexType.ne.0)then
      CALL DAREAD(IDAF,IODA,X(LFBM),L2,97,0)                             !read in the Z dipole integrals to -LFBM-
      CALL EXPND(X(LFBM),X(LFBMMO),L1,0)                                 !and expand it to -LFBMMO- to a full square form
      trdipZ = ddot(l3,x(lwrk1),1,x(lfbmmo),1)                           !take trace of the product D^{tran}*Dip_z
      if(rexType.ne.0)then                                               !do electronic dipole moments of PPS, OSS, abd (possibly) DES states
       dipPPSZ = ddot(l3,x(lwrk2),1,x(lfbmmo),1)*0.5d0
       dipOSSZ = ddot(l3,x(lwrk3),1,x(lfbmmo),1)*0.5d0
      if(rexType.eq.2)
     * dipDESZ = ddot(l3,x(lwrk4),1,x(lfbmmo),1)*0.5d0
      endif                                                              !if(rexType.ne.0)then
      au2debye = 2.541766d0                                              !conversion from au to debye; common in gamess-us
      if(rexType.eq.0)then                                               !coeff before trans.dipole; cftran
       trdipX = cftran*trdipX
       trdipY = cftran*trdipY
       trdipZ = cftran*trdipZ
      elseif(rexType.eq.1)then                                           !if(rexType.eq.0)then
       cftTrn = cftran*(s2sare(1,1)**2 - s2sare(1,2)**2)                 !coeff before trans. dipole; (a_11^2 - a_12^2)*cftran
       cftPPSOSS = s2sare(1,1)*s2sare(1,2)                               !coeff before PPS/OSS dipoles
       trdipX = cftTrn*trdipX + cftPPSOSS*(dipPPSX - dipOSSX)
       trdipY = cftTrn*trdipY + cftPPSOSS*(dipPPSY - dipOSSY)
       trdipZ = cftTrn*trdipZ + cftPPSOSS*(dipPPSZ - dipOSSZ)
      elseif(rexType.eq.2)then                                           !if(rexType.eq.0)then
       cftran2 = (sqrt(1.d0*dnr)+sqrt(1.d0*dns))*sqrt(0.5d0)             !coefficient before transition density btw OSS and DES states
       cftTrn1 = s3sare(1,1)*s3sare(2,2)+s3sare(1,2)*s3sare(2,1)         !for PPS to OSS tran dipole
       cftTrn2 = s3sare(2,1)*s3sare(3,2)+s3sare(2,2)*s3sare(3,1)         !for OSS to DES tran dipole
       cftTrn = cftTrn1*cftran + cftran2*cftTrn2
       cftPPS = s3sare(1,1)*s3sare(1,2)
       cftOSS = s3sare(2,1)*s3sare(2,2)
       cftDES = s3sare(3,1)*s3sare(3,2)
      trdipX=cftTrn*trdipX+cftPPS*dipPPSX+cftOSS*dipOSSX+cftDES*dipDESX
      trdipY=cftTrn*trdipY+cftPPS*dipPPSY+cftOSS*dipOSSY+cftDES*dipDESY
      trdipZ=cftTrn*trdipZ+cftPPS*dipPPSZ+cftOSS*dipOSSZ+cftDES*dipDESZ
      endif                                                              !if(rexType.eq.0)then
      trdipX = trdipX*au2debye
      trdipY = trdipY*au2debye
      trdipZ = trdipZ*au2debye
      trdipT = sqrt(trdipX**2 + trdipY**2 + trdipZ**2)
      if(rexType.ne.0)then
      if(masprt)
     *write(iw,
     *'(/1x,a," transition dipole between the S1 and S0 states")')
     * rexprefix
      else
      if(masprt)
     *write(iw,
     *'(/1x,a," transition dipole between the OSS and PPS states")')
     * rexprefix
      endif
      if(masprt)
     *write(iw,'(8x,"Dipole = (",F12.8,",",F12.8,",",F12.8,"); |D| ="
     *,F12.8," Debye")') trdipX,trdipY,trdipZ,trdipT
C
      if(dirscf) then
      if(masprt)
     * write(iw,'(/1x,a)')ddashes(:132)
      else
      if(masprt)
     * write(iw,'(/1x,a)')ddashes(:116)
      endif
C
      ETOT = rexenergy                                                   !put rexEnergy to ETOT
C
  650 NDAF = 15                                                          !SAVE ORBITALS AND DENSITIES AND ENERGIES; leap here, if single state REKS
      CALL SCFSAV(X(LVA),X(LDA),X(LEA),NDAF,L1,L2,L3)
      NDAF = 19
      CALL SCFSAV(X(LVA),X(LDA),X(LEA),NDAF,L1,L2,L3)                    !save α orbitals as β
C   do SA-REKS Lagrangian                                                !the Lagrangian was checked by 1-e gradient; compared to RHF
      CALL DAREAD(IDAF,IODA,X(LWRK2),L2,36,0)                            !read the SA-REKS Lagrangian in -LWRK2- and convert to AO rep.
      CALL EXPND(X(LWRK2),X(LWRK1),L1,0)                                 !unpack lower to square in -LWRK1-
c       if(MASWRK)then
c         write(iw,*) 'test L.MO'
c         CALL PRSQ(X(LWRK1),L1,L1,L1)                                   !test print
c       endif
      CALL DAWRIT(IDAF,IODA,X(LWRK1),L3,704,0)                           !save the SA Lagrangian in MO to 704 for further use
      CALL DGEMM('N','N',L1,L1,L1,ONE,X(LVA),L1,X(LWRK1),L1,ZERO,        !here, a half transform is done; it's C*L.MO*C^†, not S*C*L.MO*C^†*S
     *           X(LWRK2),L1)
      CALL TRPOSQ(X(LWRK2),L1)
      CALL DGEMM('N','N',L1,L1,L1,-ONE,X(LVA),L1,X(LWRK2),L1,ZERO,       !the Lagrangian in gamess should be negative (-)
     *           X(LWRK1),L1)
      CALL CPYSQT(X(LWRK1),X(LWRK2),L1,1)                                !pack square back to lower
      CALL DAWRIT(IDAF,IODA,X(LWRK2),L2,36,0)                            !save it back to where it belongs
C
      IF(DBG) THEN
         write(IW,'(/9x,a,/10x,a,/9x,a)')
     *   dashes(:17),'REKS Lagrangian',dashes(:17)
         CALL PRTRI(X(LWRK2),L1)
      END IF
C
      CALL DAREAD(IDAF,IODA,X(LVA),L3,15,0)                              !PRINT ALPHA ORBITALS
      CALL DAREAD(IDAF,IODA,X(LDA),L2,16,0)
      CALL DAREAD(IDAF,IODA,X(LEA),L1,17,0)
      IF(MINMEM) THEN
         CALL DAREAD(IDAF,IODA,X(LQ),L3,45,0)
         CALL DAREAD(IDAF,IODA,X(LS),L2,12,0)
      END IF
      CALL SYMMOS(X(LIWRK),X(LQ),X(LS),X(LVA),X(LSCR),L0,L1,L0,L1)       !CALLING -SYMMOS- DESTROYS THE -Q- MATRIX
      CALL DAWRIT(IDAF,IODA,X(LIWRK),L1,255,1)
      CALL DAWRIT(IDAF,IODA,IRMON,L1,356,1)
      IF(ETOLLZ.NE.0) CALL LZMOS(X(LEA),X(LVA),L0,L1,L1,L2,L3)
      IF(SOME) THEN
         WRITE(IW,9170)
         CALL PREVS(X(LVA),X(LEA),X(LIWRK),L0,L1,L1,1)
      END IF
      IF(PRDENS) THEN
         WRITE(IW,9180)
         CALL PRTRIL(X(LDA),L1)
      END IF
C
c 700 CONTINUE
      IF(MPUNCH.NE.0  .AND.  MASPRT) THEN                                !PUNCH THE ORBITALS
         CALL TMDATE(TIMSTR)
         IF(NFG.NE.0) THEN
            WRITE(IP,8005) TIMSTR,ICURFG,JCURFG,KCURFG,line(:lenline),
     *                     ETOT,EN,ITER
         ELSE
         WRITE(IP,8000) TIMSTR,TITLE,line(:lenline),ETOT,EN,ITER
         ENDIF
         LPUN = NA
         IF(MPUNCH.EQ.2) LPUN = L0
         IF(RUNTYP.NE.AMEX) THEN
                              WRITE(IP,8010) ' $VEC   '
         ELSE
            IF(MEXSTATE.EQ.1) WRITE(IP,8010) ' $VEC1  '
            IF(MEXSTATE.EQ.2) WRITE(IP,8010) ' $VEC2  '
         END IF
         CALL PUSQL(X(LVA),LPUN,L1,L1)
                           WRITE(IP,8010) ' $END   '
      END IF
C
      IF(SOME .AND. NEORUN.EQ.1) THEN                                    !IF NEO RUN OUPUT NUCLEAR MOS
         CALL NMOOUT(MPUNCH)
      END IF
C
      END IF                                                             !IF(EXETYP.NE.CHECK)
      IF(CVGED) IREST = 0
      ETHRSH = ETHSAV
      SOGTOL = SOGSAV
C
      IF(FT15OP) THEN
         TDSKWRK = DSKWRK
         DSKWRK  = ISGDDI
         CALL SEQCLO(NFT15,'DELETE')
         DSKWRK = TDSKWRK
         FT15OP=.FALSE.
      END IF
      IF(.NOT.NOTOPN) CALL RACLOS(20,'DELETE')
      NOTOPN=.TRUE.
C
      IF (QFMM) THEN                                                     !SAVE FMM DATA FOR DERIVATIVE CALCULATIONS
         CALL DERCHK(NDER)
         IF (NDER.GT.0) THEN
            CALL DENWZP(SCFTYP,NCXYZ,X(LIYP),X(LINDX2),
     *         X(LIDXIJK),X(LIDXWS),X(LCXYZ),X(LIBS),X(LYZTBL),NTMPL,
     *         X(LYP),X(LZP),L2,X(LDA),X(LDB),NTBOX,MAXWS,X(LYZPNT),
     *         X(LF),X(LG),X(LZLL),X(LCLM),X(LFLM),X(LTMPGPS),NSBOX,
     *         X(LISP),X(LIPP))
            CALL SHLDEN(SCFTYP,X(LDA),X(LDB),DUMMY,X(LDSH),IA,L1,L2,
     *         NSH2,1)
            CALL WRTFMM(X(LYP),X(LZP),NP,NTMPL,NCXYZ,X(LIYP),
     *         X(LINDX2),X(LIDXIJK),X(LIDXWS),X(LIBS),
     *         X(LYZTBL),NTBOX,MAXWS,X(LYZPNT),X(LCLM),
     *         X(LIDXBOX),NSBOX,X(LISP),X(LIPP),X(LEBOX),
     *         X(LMBOX),X(LINDX),X(LNUMWS),X(LDSH),X(LSLIST),
     *         X(LIDXSHL),X(LSLN),X(LXINTS),NSHELL,NSH2,NSHL2,NS)
         END IF
      END IF
C
      CALL RETFM(NEED)                                                   !RETURN ALL FAST MEMORY
      if(masprt)
     *WRITE(IW,FMT='('' ...... END OF REKS SCF CALCULATION ......'')')   !BANNER THAT ENDS THE CALCULATION
      CALL TIMIT(1)
C
c       if(cllCPR.and.rexEKT)then
c     write(iw,*)
c    *' ==> SA/SSR gradient is currently incompatible with EKT!'
c     call ABRT
c       endif
      if(cllCPR.OR.rlxDen.OR.rexEKT)                                     !call CP-REKS if gradient, or relaxed density matrix, or EKT requested
     *  call cpreks(cllCPR)                                              !test call
C
      RETURN
C
 8000 FORMAT('--- OPEN SHELL ORBITALS --- GENERATED AT ',3A8/10A8/
     *    'E(',A,')=',F20.10,', E(NUC)=',F16.10,',',I5,' ITERS')
 8005 FORMAT('--- OPEN SHELL ORBITALS --- GENERATED AT ',3A8/,
     *       'FMO ORBITALS',3I6/,
     *       'E(',A,')=',F20.10,', E(NUC)=',F16.10,',',I5,' ITERS')
 8010 FORMAT(A8)
C
 9000 FORMAT(/17X,20(1H-)/10X,A11,' SCF CALCULATION'/17X,20(1H-))
 9005 FORMAT(/1X,'OVERLAPPING STORAGE FOR H,S,Q MATRICES,',
     *        1X,'TO SAVE',I10,' WORDS.'/1X,'ADDING THIS AMOUNT',
     *        1X,'MAY IMPROVE EFFICIENCY SOMEWHAT.')
 9006 FORMAT(1X,'* * * ERROR * * *'/
     *       1X,'IT IS ILLOGICAL FOR THE GRID SWITCHING VALUE',
     *          ' $DFT SWITCH=',1P,E7.1/
     *       1X,'    TO BE SMALLER THAN THE DENSITY CONVERGEN',
     *          'CE $SCF CONV=',1P,E7.1)
 9012 FORMAT(5X,'DENSITY MATRIX CONV=',1P,E10.2/
     *       5X,'   FOCK MATRIX CONV=',1P,E10.2)
 9013 FORMAT(5X,'    DENSITY MATRIX CONVERGENCE THRESHOLD=',1P,E10.2/
     *       5X,'       FOCK MATRIX CONVERGENCE THRESHOLD=',1P,E10.2/
     *       5X,'COARSE -> FINE DFT GRID SWITCH THRESHOLD=',1P,E10.2,
     *          ' (SWITCH IN $DFT)'/
     *       5X,'              HF -> DFT SWITCH THRESHOLD=',1P,E10.2,
     *          ' (SWOFF IN $DFT)')
 9015 FORMAT(5X,'VTSCAL=',L1,'  MAXVT=',I3,'  VTCONV=',F12.8)
 9017 FORMAT(5X,'SOSCF WILL OPTIMIZE',I8,' ORBITAL ROTATION ANGLES.',
     *          ' SOGTOL=',1P,E10.3)
 9018 FORMAT(5X,'SOSCF WILL OPTIMIZE',I8,' ALPHA AND',I8,
     *          ' BETA ROTATION ANGLES.'/5X,'SOGTOL=',1P,E10.3)
 9020 FORMAT(/5X,'NUCLEAR ENERGY = ',F20.10/
     *       5X,'MAXIT =',I3,5X,'NPUNCH=',I3,5X,'MULT=',I3/
     *       5X,'EXTRAP=',L1,'  DAMP=',L1,'  SHIFT=',L1,
     *          '  RSTRCT=',L1,'  DIIS=',L1,'  SOSCF=',L1)
 9030 FORMAT(5X,'ROHF CANONICALIZATION PARAMETERS'/
     *       5X,8X,'C-C',5X,'O-O',5X,'V-V'/
     *       5X,'ALPHA',3F8.4/5X,'BETA ',3F8.4)
 9035 FORMAT(1X,'QFMM WILL USE ADDITIONAL MEMORY OF',I10,' WORDS')
 9039 FORMAT(/1X,'THE BUFFER ZONE DENSITY WILL BE EXCLUDED FROM THE ',
     *        'INTERACTION WITH THE EFP')
 9040 FORMAT(5X,'MEMORY REQUIRED FOR REKS ITERS=',I10,' WORDS.')
 9045 FORMAT(/1X,'DIRECT SCF CALCULATION, SCHWRZ=',L1,'   FDIFF=',L1)
 9048 FORMAT(/' ITER EX      TOTAL ENERGY        E CHANGE  DENSITY ',
     *        'CHANGE    ORB. GRAD       VIR. SHIFT       DAMPING')
 9049 FORMAT(/112X,'NONZERO     BLOCKS'/
     *        ' ITER EX      TOTAL ENERGY        E CHANGE  DENSITY ',
     *        'CHANGE    ORB. GRAD       VIR. SHIFT       DAMPING',
     *        '        INTEGRALS    SKIPPED')
 9050 FORMAT(/' ITER EX     NR/2           TOTAL ENERGY    OFFDIAG ',
     *        'FOCK      E CHANGE    DENSITY ',
     *        'CHANGE  DIIS ERROR     VIR. SHIFT')
 9055 FORMAT(/112X,'NONZERO     BLOCKS'/
     *        ' ITER EX     NR/2           TOTAL ENERGY        E ',
     *        'CHANGE  DENSITY CHANGE    DIIS ERROR      VIR. SHIFT',
     *        '        INTEGRALS    SKIPPED')
 9058 FORMAT(/' ITER EX      TOTAL ENERGY        E CHANGE  DENSITY ',
     *        'CHANGE    ORB. GRAD ')
 9059 FORMAT(/80X,'NONZERO     BLOCKS'/
     *        ' ITER EX      TOTAL ENERGY        E CHANGE  DENSITY ',
     *        'CHANGE    ORB. GRAD       INTEGRALS    SKIPPED')
 9060 FORMAT(/' ITER EX      TOTAL ENERGY   OFFDIAG FOCK',
     *        '      E CHANGE     DENSITY CHANGE  DIIS ERROR')
 9065 FORMAT(/80X,'NONZERO     BLOCKS'/
     *        ' ITER EX      TOTAL ENERGY        E CHANGE  DENSITY ',
     *        'CHANGE    DIIS ERROR      INTEGRALS    SKIPPED')
 9070 FORMAT(1X,2I3,F14.9,F20.10,F14.9,F17.10,2F14.9,F16.9,I15,I11)
 9075 FORMAT(1X,2I3,F20.10,F14.9,F17.10,2F14.9,F16.9,I15,I11)
 9080 FORMAT(/10X,17(1H-)/10X,17HDENSITY CONVERGED/10X,17(1H-))
 9090 FORMAT(/10X,21(1H-)/10X,21HFOCK MATRIX CONVERGED/10X,21(1H-))
 9100 FORMAT(/10X,14(1H-)/10X,14HDIIS CONVERGED/10X,14(1H-))
 9105 FORMAT(/1X,'REKS HAS CONVERGED, NOW COMPUTING EXACT ALPHA,BETA',
     *           ' FOCK MATRICES'/
     *        1X,'FOR USE DURING THE COUPLED CLUSTER CALCULATION',
     *           ' THAT FOLLOWS.')
 9110 FORMAT(1X,'... SCF HAS ALMOST CONVERGED BUT ROUTINE STOPPED',
     *       1X,' FOR TIMLIM BEFORE ULTIMATE CYCLE ...')
 9120 FORMAT(/1X,'SCF IS UNCONVERGED, TOO LITTLE TIME')
 9130 FORMAT(/1X,'SCF IS UNCONVERGED, TOO MANY ITERATIONS')
 9135 FORMAT(/10X,'INDUCED DIPOLE',5X,'ATOMIC UNITS',10X,'DEBYE',
     *       /10X,'      X(IND)  ',5X,F12.5,5X,F12.5,
     *       /10X,'      Y(IND)  ',5X,F12.5,5X,F12.5,
     *       /10X,'      Z(IND)  ',5X,F12.5,5X,F12.5)
 9140 FORMAT(/1X,'FINAL ',A,' ENERGY IS',F20.10,' AFTER',I4,
     *           ' ITERATIONS')
 9141 FORMAT(/1X,'PURE DFT ',A,' ENERGY IS',F20.10,' AFTER',I4,
     *           ' ITERATIONS')
 9142 FORMAT(1X,'THIS IS ONLY THE DFT PART OF THE TOTAL ENERGY,')
 9143 FORMAT(1X,'LOOK AFTER FOR THE TOTAL ENERGY WITH MP2 ADDITION')
 9145 FORMAT(/1X,'HEAT OF FORMATION IS',F15.5,' KCAL/MOL')
 9148 FORMAT(/10X,20(1H-)/10X,12HSPIN SZ   = ,F8.3/
     *        10X,12HS-SQUARED = ,F8.3/10X,20(1H-))
 9510 FORMAT(1X,'GRIMME''S DISPERSION ENERGY                     =',
     *           F20.10)
 9500 FORMAT(1X,'ENERGY WITHOUT GRIMME''S DISPERSION CORRECTION IS',
     *           F20.10)
 9550 FORMAT(1X,'ENERGY WITHOUT LRD CORRECTION IS',F20.10)
 9150 FORMAT(/1X,' ----- ALPHA SET ----- ')
 9160 FORMAT(/1X,' ----- BETA SET ----- ')
 9170 FORMAT(/10X,12(1H-)/10X,12HEIGENVECTORS/10X,12(1H-))
 9180 FORMAT(/10X,14(1H-)/10X,14HDENSITY MATRIX/10X,14(1H-))
 9190 FORMAT(/10X,43(1H-)/10X,'UHF NATURAL ORBITALS AND OCCUPATION',
     *        ' NUMBERS'/10X,43(1H-))
 9200 FORMAT(10X,15(1H-),'START SECOND ORDER SCF',15(1H-))
 9210 FORMAT(1X,'ALPHA DFT EXCHANGE + CORRELATION ENERGY = ',F20.10)
 9220 FORMAT(1X,' BETA DFT EXCHANGE + CORRELATION ENERGY = ',F20.10)
 9235 FORMAT(1X,'NEXT ENERGY RISES: SWITCHING DIIS OFF FOR',I3,
     *          ' ITERS, NEW ETHRSH=',1P,E9.2)
 9236 FORMAT(1X,'NEXT ENERGY RISES: SWITCHING SOSCF OFF FOR',I3,
     *          ' ITERS, NEW SOGTOL=',1P,E9.2)
 9230 FORMAT(/1X,'THE CONVERGED ORBITALS WILL UNDERGO GUEST/SAUNDERS'/
     *        1X,'CANONICALIZATION FOR ZAPT PERTURBATION THEORY.')
 9310 FORMAT(1X,'TOTAL ELECTRON NUMBER             = ',F20.10)
 9311 FORMAT(/10X,'TOTAL KINETIC ENERGY DENSITY      = ',F20.10,1X/)
 9312 FORMAT(1X,'TOTAL KINETIC ENERGY DENSITY      = ',F20.10,1X)
 9320 FORMAT(1X,'DFT EXCHANGE + CORRELATION ENERGY = ',F20.10)
 9330 FORMAT(1X,'DFT CODE IS SWITCHING BACK TO THE FINE GRID')
 9335 FORMAT(/1X,'DFT CODE IS SWITCHING FROM THE FINE GRID NRAD=',I3,
     *           ',  NLEB=',I5/
     *       22X,'TO THE COARSE GRID NRAD0=',I3,', NLEB0=',I5)
 9340 FORMAT(/1X,'DFT CODE IS SWITCHING FROM THE FINE GRID NRAD=',I3,
     *          ',  NTHE, NPHI=',2I3/
     *       22X,'TO THE COARSE GRID NRAD0=',I3,
     *          ', NTHE0,NPHI0=',2I3)
 9350 FORMAT(1X,'DFT IS SWITCHED OFF, PERFORMING PURE SCF UNTIL SWOFF',
     *          ' THRESHOLD IS REACHED.')
 9355 FORMAT(1X,'CONVERGED TO SWOFF, SO DFT CALCULATION IS NOW',
     *          ' SWITCHED ON.')
 9380 FORMAT(9X,'SWDIIS THRESHOLD HAS BEEN REACHED:',
     *          ' SWITCHING DIIS OFF.')
 9390 FORMAT(1X,'ULTIMATE SCF CYCLE OF MAKEFP SWITCHES TO THE',
     *          ' ROOTHAAN CANONICALIZATION...')
 9400 FORMAT(5X,'TIME TO FORM FOCK OPERATORS=',F10.1,' SECONDS (',
     *          F10.1,' SEC/ITER)'/
     *       5X,'FOCK TIME ON FIRST ITERATION=',F10.1,
     *          ', LAST ITERATION=',F10.1/
     *       5X,'TIME TO SOLVE SCF EQUATIONS=',F10.1,' SECONDS (',
     *          F10.1,' SEC/ITER)')
 9410 FORMAT(5X,'TIME TO FORM FOCK OPERATORS=',F10.1,' SECONDS (',
     *          F10.1,' SEC/ITER)'/
     *       5X,'TIME TO SOLVE SCF EQUATIONS=',F10.1,' SECONDS (',
     *          F10.1,' SEC/ITER)')
C
 1005 FORMAT(1X,A10,F18.10)
 1000 FORMAT(1X,'CALCULATED ',F10.5,' THEORETICAL ',
     * F10.5,' NOT RENORMALIZED',F10.5)
 1100 FORMAT(1X,'CALCULATED ',F10.5,' THEORETICAL ',
     * F10.5,' RENORMALIZED',F10.5)
 1110 FORMAT(1X,'CALCULATED',F10.5,' (THEOR= ',
     * F10.5,')  ESCAPED',F10.5,' FINAL', F10.5)
C
      END SUBROUTINE REKSSCF
C*MODULE REKS    *DECK REXDMAT
C>
C> @author     Seunghoon Lee
C>
C> @date       2021, Oct
C>
      SUBROUTINE REXDMAT(D,V,IA,L1,L2,L3,N)
      USE prec, ONLY: dp
      IMPLICIT NONE
C
      INTEGER :: L1, L2, L3, N
      REAL(KIND=dp), DIMENSION(L2,2,2) :: D
      INTEGER, DIMENSION(L1) :: IA
      REAL(KIND=dp), DIMENSION(L1,L1) :: V
C
      INTEGER :: I, IJ, J, K
C
      CALL VCLR(D,1,L2*4)
C
      DO I=1,L1                                                          !build core density matrix D00
      DO J=1,I
         IJ=IA(I)+J
         DO K=1,N-1
            D(IJ,1,1)=D(IJ,1,1)+V(I,K)*V(J,K)
         ENDDO
      ENDDO
      ENDDO
C
      CALL DCOPY(L2,D(1,1,1),1,D(1,2,1),1)
      CALL DCOPY(L2,D(1,1,1),1,D(1,1,2),1)
C
      DO I=1,L1
      DO J=1,I
         IJ=IA(I)+J
         D(IJ,2,1)=D(IJ,2,1)+V(I,N)*V(J,N)                               !build D10
         D(IJ,1,2)=D(IJ,1,2)+V(I,N+1)*V(J,N+1)                           !build D01
      ENDDO
      ENDDO
C
      CALL DCOPY(L2,D(1,2,1),1,D(1,2,2),1)
C
      DO I=1,L1
      DO J=1,I
         IJ=IA(I)+J
         D(IJ,2,2)=D(IJ,2,2)+V(I,N+1)*V(J,N+1)                           !build D11
      ENDDO
      ENDDO
      RETURN
      END SUBROUTINE REXDMAT
C*MODULE REKS    *DECK REX2E
C>
C> @author     Seunghoon Lee
C>
C> @date       2021, Oct
C>
      SUBROUTINE REX2E(FAM,FBM,DJT,DKT,IA,L1,L2,NMIC)
      USE comm_REKSCM, ONLY: NMICRO, MTTYP, WPPS, WOSS, G1, DNR,
     * DNS, DELTA, FR, FS
      USE prec, ONLY: dp
      IMPLICIT NONE
C
C
      INTEGER :: L1, L2, NMIC
      REAL(KIND=dp), DIMENSION(L2,2,2) :: DJT, DKT
      REAL(KIND=dp), DIMENSION(L2,NMIC) :: FAM, FBM
      INTEGER, DIMENSION(*) :: IA
C
      INTEGER :: I, IJ, IRA, IRB, ISA, ISB, J, L
C
      DO L=1,NMIC  ! L     - microstates
         IRA=NMICRO(1,1,L)
         ISA=NMICRO(2,1,L)
         IRB=NMICRO(1,2,L)
         ISB=NMICRO(2,2,L)
C
         DO I=1,L1
         DO J=1,I
            IJ=IA(I)+J
            FAM(IJ,L)=DJT(IJ,IRA,ISA)+DJT(IJ,IRB,ISB)-DKT(IJ,IRA,ISA)
            FBM(IJ,L)=DJT(IJ,IRA,ISA)+DJT(IJ,IRB,ISB)-DKT(IJ,IRB,ISB)
         ENDDO
         ENDDO
      ENDDO
      RETURN
      END SUBROUTINE REX2E
C*MODULE REKS    *DECK HSTARREX
C>
C> @brief      Clone of HSTAR
C>
C> @author     Seunghoon Lee
C>
C> @date       2021, Oct
C
      SUBROUTINE HSTARREX(D1,DJ1,DK1,D2,DJ2,DK2,D3,DJ3,DK3,D4,DJ4,DK4,  &
     &                    XX,IX,NINTMX,IA,NOPK,DUPAO)
C
      USE camdft, ONLY: CAMFLAG
      USE lrcdft, ONLY: LCFLAG, LRFILE
      USE mx_limits, ONLY: MXATM
      USE constants, ONLY: half
      USE prec, ONLY: dp
      IMPLICIT NONE
C
      REAL(KIND=dp), DIMENSION(3,MXATM) :: C
      REAL(KIND=dp), DIMENSION(20) :: DFTTYP
      LOGICAL :: DSKWRK, GOPARR, MASWRK
      REAL(KIND=dp) :: EXENA, EXENB, EXENC
      INTEGER, DIMENSION(MXATM) :: IAN
      INTEGER :: IBTYP, ICH, IDAF, IDFT34, ININTIC, IP, IPK, IPRI,      &
     &           IPTIM, IR, IS, IW, LABSIX, LABSIZ, LBUFPIC, LDOS,      &
     &           LIXIC, MASTER, ME, MUL, NA, NASPIN, NAT, NATMM,        &
     &           NAUXFUN, NAUXSHL, NAV, NB, NBNDAB, NCT, NE, NELONG,    &
     &           NINTIC, NINTIX, NPROC, NQMT, NTMLB, NUM, NXXIC
      INTEGER, DIMENSION(950) :: IODA
      INTEGER, DIMENSION(7) :: NORDER
      REAL(KIND=dp), DIMENSION(MXATM) :: ZAN
      COMMON /DFTPAR/ DFTTYP, EXENA, EXENB, EXENC, IDFT34, NAUXFUN,     &
     &                NAUXSHL
      COMMON /ELGPMT/ NELONG, NATMM, NASPIN, NCT, NBNDAB, NTMLB, IPRI,  &
     &                LDOS
      COMMON /INFOA / NAT, ICH, MUL, NUM, NQMT, NE, NA, NB, ZAN, C, IAN
      COMMON /INT2IC/ NINTIC, ININTIC, NXXIC, LBUFPIC, LIXIC, LABSIX,   &
     &                NINTIX
      COMMON /IOFILE/ IR, IW, IP, IS, IPK, IDAF, NAV, IODA
        LOGICAL :: LRINT
      COMMON /NLRCF / LRINT
      COMMON /ORDOPT/ NORDER
      COMMON /PAR   / ME, MASTER, NPROC, IBTYP, IPTIM, GOPARR, DSKWRK,  &
     &                MASWRK
      COMMON /PCKLAB/ LABSIZ
C
      LOGICAL :: DUPAO
      INTEGER :: NINTMX, NOPK
      REAL(KIND=dp), DIMENSION(*) :: D1, D2, D3, D4, DJ1, DJ2, DJ3, DJ4, &
     &                              DK1, DK2, DK3, DK4
      INTEGER, DIMENSION(*) :: IA, IX
      REAL(KIND=dp), DIMENSION(NINTMX) :: XX
C
      REAL(KIND=dp) :: CSCALT, HFSCAL, VAL, VAL2, VAL4
      INTEGER :: I, II, IPACK, IPCOUNT, J, JPACK, K, KPACK, L,          &
     &           LPACK, M, NIJ, NIK, NIL, NINT, NJK, NJL, NKL,          &
     &           NPACK, NPACKIC, NRECRD, NUMTRI, NXX, NFILE,            &
     &           LABEL, NPAC
C
C         ----- FORM THE SKELETON CLOSED SHELL FOCK MATRIX -----
C     THIS MEANS ONLY SYMMETRY UNIQUE 2E- INTEGRALS ARE PROCESSED, SO
C     THE RESULTING MATRIX MUST BE SYMMETRIZED BY -SYMH- AND THEN WILL
C     NEED THE 1 ELECTRON TERMS ADDED IN ORDER TO FINISH A FOCK MATRIX.
C
C     INDICES IN AO INTEGRAL LABELS ARE IN STANDARD ORDER,
C               I.GE.J , K.GE.L , (IJ).GE.(KL)
C     THIS ROUTINE ALSO HAS AN OPTION TO PROCESS A -P- INTEGRAL LIST,
C     AND THERE IS SPECIAL CODE FOR FORMING A MATRIX WITH A DESIRED
C     ADMIXTURE OF THE EXCHANGE TERM (FOR PURPOSES OF DFT MATRICES).
C
C     IN CASE A DUPLICATED AO INTEGRAL LIST, AS OPPOSED TO A DISTRIBUTED
C     LIST, IS TO BE PROCESSED, WE MUST SKIP RECORDS IN ORDER TO MAKE THE
C     CPU WORK OF A FOCK BUILD BE PARALLELIZED.  DUPAO=.TRUE. IN CALL.
C
      if(.not.LRINT) CALL SEQREW(IS)
      if(     LRINT) CALL SEQREW(LRFILE)
C
      IF(NELONG.GT.1) THEN
         NRECRD = 1
      END IF
C
      NUMTRI = (NUM*NUM+NUM)/2
      I = 0
      J = 0
      K = 0
      L = 0
      NXX = 0
C
      IF(DUPAO) IPCOUNT = ME
C
C     ----- INTEGRALS ARE IN STANDARD -J- LIST (NOPK=.TRUE.) -----
C     ----- DOING SPECIAL DFT CASE WHERE EXCHANGE CONTRIBUTION  -----
C     ----- IS MULTIPLIED BY SOME NUMBER BETWEEN ZERO AND ONE   -----
C
      HFSCAL=DFTTYP(3)
      CSCALT=1.0D+00
      IF(LCFLAG) THEN
         IF(LRINT) THEN
            HFSCAL=1.0D+00
            CSCALT=0.0D+00
         ELSE
            HFSCAL=0.0D+00
            CSCALT=1.0D+00
         ENDIF
      ENDIF
      IF(CAMFLAG.AND.LRINT) CSCALT=0.0D+00
C     FIRST, PROCESS IN CORE INTEGRALS
      NPACKIC=0
C
  410 IF(DUPAO) THEN
         IPCOUNT = IPCOUNT + 1
         IF (MOD(IPCOUNT,NPROC).NE.0) THEN
            if(.not.LRINT) READ (IS)     NXX
            if(     LRINT) READ (LRFILE) NXX
            GO TO 480
         END IF
      END IF
      IF(NPACKIC.EQ.0.AND.NINTIC.NE.0) THEN
        NXX=NXXIC
        NINT=NXX
      ELSE
C       IF ALL IN CORE BAIL OUT
        IF(NINTIC.NE.0.AND.NXXIC.LT.NINTIC) GO TO 490
        if(.not.LRINT)
     *    CALL PREAD(IS,XX(NINTIC+1),IX(ININTIC+1),NXX,NINTMX)
        if(LRINT)
     *    CALL PREAD(LRFILE,XX(NINTIC+1),IX(ININTIC+1),NXX,NINTMX)
        IF(NXX .EQ. 0) GO TO 490
        NINT = IABS(NXX)
        IF(NINT .GT. NINTMX) CALL ABRT
      END IF
C
      DO M = 1,NINT
         NPACK = NPACKIC+M
         IF (LABSIZ .EQ. 2) THEN
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
         I = IPACK
         J = JPACK
         K = KPACK
         L = LPACK
C
         VAL = XX(NPACK)
         NIJ = IA(I)+J
         NKL = IA(K)+L
C
C           USING SQUARE CANONICAL INTEGRAL FILE
C
         IF(NORDER(7).EQ.1) THEN
            IF(NKL .GT. NIJ) CYCLE
            IF(I  .EQ.  J) VAL=VAL*HALF
            IF(K  .EQ.  L) VAL=VAL*HALF
            IF(NIJ.EQ.NKL) VAL=VAL*HALF
         END IF
C
         NIK = IA(I)+K
         NIL = IA(I)+L
         IF(J .LT. K) THEN
            NJK = IA(K)+J
            IF(J .LT. L) THEN
               NJL = IA(L)+J
            ELSE
               NJL = IA(J)+L
            END IF
         ELSE
            NJK = IA(J)+K
            NJL = IA(J)+L
         END IF
C
         VAL2     = VAL  +  VAL
         VAL4     = VAL2 + VAL2
         DJ1(NIJ) = DJ1(NIJ)+VAL4*D1(NKL)*CSCALT
         DJ1(NKL) = DJ1(NKL)+VAL4*D1(NIJ)*CSCALT
         DK1(NIK) = DK1(NIK)+VAL2*D1(NJL)*HFSCAL
         DK1(NIL) = DK1(NIL)+VAL2*D1(NJK)*HFSCAL
         DK1(NJK) = DK1(NJK)+VAL2*D1(NIL)*HFSCAL
         DK1(NJL) = DK1(NJL)+VAL2*D1(NIK)*HFSCAL
C
         DJ2(NIJ) = DJ2(NIJ)+VAL4*D2(NKL)*CSCALT
         DJ2(NKL) = DJ2(NKL)+VAL4*D2(NIJ)*CSCALT
         DK2(NIK) = DK2(NIK)+VAL2*D2(NJL)*HFSCAL
         DK2(NIL) = DK2(NIL)+VAL2*D2(NJK)*HFSCAL
         DK2(NJK) = DK2(NJK)+VAL2*D2(NIL)*HFSCAL
         DK2(NJL) = DK2(NJL)+VAL2*D2(NIK)*HFSCAL
C
         DJ3(NIJ) = DJ3(NIJ)+VAL4*D3(NKL)*CSCALT
         DJ3(NKL) = DJ3(NKL)+VAL4*D3(NIJ)*CSCALT
         DK3(NIK) = DK3(NIK)+VAL2*D3(NJL)*HFSCAL
         DK3(NIL) = DK3(NIL)+VAL2*D3(NJK)*HFSCAL
         DK3(NJK) = DK3(NJK)+VAL2*D3(NIL)*HFSCAL
         DK3(NJL) = DK3(NJL)+VAL2*D3(NIK)*HFSCAL
C
         DJ4(NIJ) = DJ4(NIJ)+VAL4*D4(NKL)*CSCALT
         DJ4(NKL) = DJ4(NKL)+VAL4*D4(NIJ)*CSCALT
         DK4(NIK) = DK4(NIK)+VAL2*D4(NJL)*HFSCAL
         DK4(NIL) = DK4(NIL)+VAL2*D4(NJK)*HFSCAL
         DK4(NJK) = DK4(NJK)+VAL2*D4(NIL)*HFSCAL
         DK4(NJL) = DK4(NJL)+VAL2*D4(NIK)*HFSCAL
C
      END DO
C
  480 IF(NPACKIC.EQ.0) NPACKIC=NINTIC
C     SWITCH TO INTEGRALS ON DISK
      IF(NXX .GT. 0) GO TO 410
C
C     DONE CONSTRUCTING MATRIX, OFF DIAGONAL ELEMENTS NEED HALVING
C
  490 if(.not.LRINT) then
        CALL DSCAL(NUMTRI,HALF,DJ1,1)
        CALL DSCAL(NUMTRI,HALF,DJ2,1)
        CALL DSCAL(NUMTRI,HALF,DJ3,1)
        CALL DSCAL(NUMTRI,HALF,DJ4,1)
        CALL DSCAL(NUMTRI,HALF,DK1,1)
        CALL DSCAL(NUMTRI,HALF,DK2,1)
        CALL DSCAL(NUMTRI,HALF,DK3,1)
        CALL DSCAL(NUMTRI,HALF,DK4,1)
        II = 0
        DO I=1,NUM
           II = II + I
           DJ1(II) = DJ1(II)+DJ1(II)
           DJ2(II) = DJ2(II)+DJ2(II)
           DJ3(II) = DJ3(II)+DJ3(II)
           DJ4(II) = DJ4(II)+DJ4(II)
           DK1(II) = DK1(II)+DK1(II)
           DK2(II) = DK2(II)+DK2(II)
           DK3(II) = DK3(II)+DK3(II)
           DK4(II) = DK4(II)+DK4(II)
        ENDDO
      endif
      if(.not.LRINT) CALL SEQREW(IS)
      if(     LRINT) CALL SEQREW(LRFILE)
      RETURN
      END SUBROUTINE HSTARREX
C*MODULE REKS    *DECK REXEM
C>
C> @author     Seunghoon Lee
C>
C> @date       2021, Oct
C>
      SUBROUTINE REXEM(EM,DT,FAM,FBM,H1,EXC,EXENA,EXENB,EXENC,EN,EDFA,  &
     &                 EDFB,L1,L2,NMIC,DBG)
C
      USE comm_REKSCM, ONLY: NMICRO, MTTYP, WPPS, WOSS, G1, DNR,
     * DNS, DELTA, FR, FS
      USE constants, ONLY: half
      USE prec, ONLY: dp
      IMPLICIT NONE
C
C
      LOGICAL :: DBG
      REAL(KIND=dp) :: EN, EXENA, EXENB, EXENC
      INTEGER :: L1, L2, NMIC
      REAL(KIND=dp), DIMENSION(L2,2,2) :: DT
      REAL(KIND=dp), DIMENSION(NMIC) :: EDFA, EDFB, EM, EXC
      REAL(KIND=dp), DIMENSION(L2,NMIC) :: FAM, FBM
      REAL(KIND=dp), DIMENSION(L2) :: H1
C
      REAL(KIND=dp) :: EFAM, EFBM, EH1A, EH1B, EMA, EMB
      INTEGER :: IW, L, LRA, LRB, LSA, LSB
      REAL(KIND=dp) :: TRACEP
C
      DO L=1,NMIC
         LRA=NMICRO(1,1,L)
         LSA=NMICRO(2,1,L)
         LRB=NMICRO(1,2,L)
         LSB=NMICRO(2,2,L)
C
         EH1A=TRACEP(DT(1,LRA,LSA),H1,L1)
         EH1B=TRACEP(DT(1,LRB,LSB),H1,L1)
         EFAM=TRACEP(DT(1,LRA,LSA),FAM(1,L),L1)
         EFBM=TRACEP(DT(1,LRB,LSB),FBM(1,L),L1)

         EMA=EH1A+EFAM
         EMB=EH1B+EFBM
         EM(L)=HALF*(EMA+EMB)+EXC(L)+EXENA+EXENB+EXENC
     *         -EDFA(L)-EDFB(L)+EN
         if(dbg) write(iw,'(a,i3,4f10.5)')
     *    'l,h1a,fam,h1b,fbm   ',l,eh1a,efam,eh1b,efbm
         if(dbg) write(iw,'(a,5f10.5)')
     *    'exc,exena,exenb,exenc,en',EXC(L),EXENA,EXENB,EXENC,EN
         if(dbg) write(iw,'(a,4f10.5)')
     *    'edfa,edfb',EDFA(L),EDFB(L)
      ENDDO

      RETURN
      END SUBROUTINE REXEM
C*MODULE REKS    *DECK RexSlv4x4
C>
C> @brief      REKS solver 4x4. 
C>             The following functions 
C> sgvb01,dsgvb01dx1,dsgvb01dx2,dsgvb01dy1,dsgvb01dy2,d2sgvb01dx12,
C> d2sgvb01dx22,d2sgvb01dy12,d2sgvb01dy22,d2sgvb01dx1x2,d2sgvb01dy1y2,
C> d2sgvb01dx1y1,d2sgvb01dx2y2,d2sgvb01dx1y2,d2sgvb01dx2y1,degvbdx1,
C> degvbdx2,d2egvbdx12,d2egvbdx22,d2egvbdx1x2,hgvb01,dhgvb01dx1,dhgvb01dx2,
C> dhgvb01dy1,dhgvb01dy2,d2hgvb01dx1y1,d2hgvb01dx1y2,d2hgvb01dx2y1,
C> d2hgvb01dx2y2,d2hgvb01dx1x2,d2hgvb01dy1y2,d2hgvb01dx12,d2hgvb01dx22,
C> d2hgvb01dy12,d2hgvb01dy22
C>             are used to optimize the fill numbers 
C>             in general and particular cases.
C>
C> @author     Michael Filatov
C>
C> @date       2021, Oct
C>
      SUBROUTINE RexSlv4x4(Eab,Eac,Ebd,Ecd,Dad,Dbc,Delta,CnvLim,x1,x2,  &
     &                     EREKS,IOut,MaxIt,dbg)
      USE constants, only: zero, one, two, five
      USE prec, ONLY: dp
      IMPLICIT NONE
C
      INTEGER :: IR, IW, IP, IJK, IJKT, IDAF, NAV 
      INTEGER, DIMENSION(950) :: IODA
      COMMON /IOFILE/ IR, IW, IP, IJK, IJKT, IDAF, NAV, IODA             !common block with I/O files
      LOGICAL :: DSKWRK, GOPARR, MASWRK
      INTEGER :: IBTYP, IPTIM, MASTER, ME, NPROC
      COMMON /PAR   / ME, MASTER, NPROC, IBTYP, IPTIM, GOPARR, DSKWRK,  &
     &                MASWRK
C
      REAL(KIND=dp) :: CNVLIM, DAD, DBC, DELTA, EAB, EAC, EBD, ECD,     &
     &                 EREKS, X1, X2
      LOGICAL :: DBG
      INTEGER :: IOUT, MAXIT
C
      REAL(KIND=dp) :: ARG1, ARG2, CC, CFT1, CFT2, DARG1, DARG2, DER11, &
     &                 DER12, DER21, DER22, DETH, DX1, DX2, DXNORM,     &
     &                 EGVB0, EGVB1, FACTOR, GOLDEN,  X10, X20
      REAL(KIND=dp) :: DERCFT0Z, DERCFT1Z, DERCFT2Z, EGVB
      REAL(KIND=dp), DIMENSION(4,4) :: H
      INTEGER :: I, INFO, ITER, J, MITER
      REAL(KIND=dp), DIMENSION(2) :: VEC
      REAL(KIND=dp), DIMENSION(4) :: W
      REAL(KIND=dp), DIMENSION(16) :: WORK
C
C Minimizes the REKS GVB-PP(2/4) energy w.r.t. the x1 and x2 parameters
C
C E_REKS[x1,x2] = x1*x2*Eab + x1*(1-x2)*Eac+ (1-x1)*x2*Ebd + (1-x1)*(1-x2)*Ecd
C               - Cft[x1]*Δad - Cft[x2]*Δbc
C
C by Newton-Raphson method
C The following set of equations is solved
C   Δx = H^-1 * V, where
C Δx is a vector of x1 and x2 increments
C H is the Hessian matrix
C
C  |  Cft"[x2]*Δbc  | E_REKS"[x1,x2] |
C  | E_REKS'[x1,x2] |  Cft'[x1]*Δad  |
C
C and V is a vector
C
C  | -Cft'[x1]*Δad + E_REKS"[x1,x2]*x2 + Eac -Ecd |
C  | -Cft'[x2]*Δbc + E_REKS"[x1,x2]*x1 + Ebd -Ecd |
C
C==============================================================================
C Input: Eab,Eac,Ebd,Ecd,Dad,Dbc - energies of the microstates (E's) and
C                                  coupling elements (D's)
C        Delta - parameter in the REKS function before the coupling
C                elements, Cft[x]; x1 - before Dad, x2 - before Dbc
C        CnvLim - convergence threshold in the NR solver
C        IOut   - output file
C        MaxIt  - maximum number of iterations in the NR solver
C Output:
C        x1, x2 - optimal values of x1 and x2
C        EREKS - REKS energy
C==============================================================================
C The function Cft[x] is given by
C        Cft[x] = Y**(1-(1/2)*(Y+Delta)/(1+Delta),
C        where Y = 4*x*(1-x)
C Its first derivative is
C        dCft[x]/dx = (2*(1 - 2*x)/(1 + d))*(Y**(-(d + Y)/(2*(1 + d)))*
C                   * (2 + d - Y - Y*Log[Y])
C Its second derivative is
C        d2Cft[x]/dx2 = (4/(1 + d)**2) Y**(-1-(d + Y)/(2*(1 + d)))*
C                       ((-d**2 - Y*(8 + (Y - 8)*Y) +
C                       d*(-2 + 5*(Y - 1)*Y)) -
C       Y*Log[Y]*(4 + d*(2 - 3*Y) + Y*(-7 + 2*Y) + (-1 + Y)*Y*Log[Y]))
C
C==============================================================================
      CC = Eab + Ecd - Eac - Ebd                                         !second derivative of the GVB-PP-REKS energy w.r.t. x1 and x2
C==============================================================================
C Setting the initial values of x1 and x2 by diagonalization of the
C Hamiltonian matrix:
C
C  | Eab |-Dbc |-Dad |  0  |
C  |-Dbc | Eac |  0  |-Dad |
C  |-Dad |  0  | Ebd |-Dbc |
C  |  0  |-Dad |-Dbc | Ecd |
C
C==============================================================================
      do i=1,4
       do j=1,4
        H(i,j)=0.d0                                                      !zero it out first
       enddo
      enddo
      H(1,1)=Eab                                                         !set up the Hamiltonian
      H(2,2)=Eac
      H(3,3)=Ebd
      H(4,4)=Ecd
      H(2,1)=Dbc
      H(3,1)=Dad
      H(4,2)=Dad
      H(4,3)=Dbc
      info = 0
      call dsyev('V','L',4,H,4,W,Work,16,Info)                           !diagonalize
      if (Info.ne.0) then
        write(iw,'(1x,"Error in DSYEV in REKS; INFO =",i10)')info
        call ABRT
      else
c     write(IOut,'(2x,F12.6,2x,4(2x,F12.6))')W(1),(H(j,1),j=1,4)
         x1 = H(1,1)**2 + H(2,1)**2                                      !compute occupations: x1 = C11^2 + C21^2
         x2 = H(1,1)**2 + H(3,1)**2                                      !                     x2 = C11^2 + C31^2
C
         if(abs(x1-one).le.0.1d0) x1=0.9                                 !if close to one, set it a bit apart
         if(abs(x2-one).le.0.1d0) x2=0.9
C
         if(abs(Dad).le.1.d-8) x1=1.d0                                   !set x1 to one if Δad is zero
         if(abs(Dbc).le.1.d-8) x1=1.d0                                   !set x2 to one if Δbc is zero
         egvb0=egvb(Eab,Eac,Ebd,Ecd,Dad,Dbc,delta,x1,x2)
         if(dbg) write(IOut,'(2x,a,F12.6,a,F12.6,a,F14.8)')
     $           'x1 =',x1,', x2 =',x2,', eGVB =',egvb0
         iter=0                                                          !iteration count
C
         if(dbg)write(IOut,'(a)')
     $   " iter     x1         x2            |Δ|                  eGVB"
         if(dbg)write(IOut,'(i4,4x,F10.8,1x,F10.8,4x,F14.12,1x,F20.12)')
     $          iter,x1,x2,one,egvb0
      end if
C
   10 iter=iter+1                                                        !increment iteration count
      x10=x1                                                             !save the current x values to x0
      x20=x2
C Now calculate the derivatives of Cft[x]
      arg1 = x1*(one-x1)
      arg2 = x2*(one-x2)
      darg1= one-two*x1
      darg2= one-two*x2
C
      Cft1 = derCft0z(arg1,delta)
      Cft2 = derCft0z(arg2,delta)
      Der11 = derCft1z(arg1,delta)*darg1
      Der12 = derCft1z(arg2,delta)*darg2
      Der21 = -two*derCft1z(arg1,delta)
     $      + derCft2z(arg1,delta)*darg1**2
      Der22 = -two*derCft1z(arg2,delta)
     $      + derCft2z(arg2,delta)*darg2**2
C
      detH = Der21*Der22*Dad*Dbc - CC*CC                                 !determinant of the Hessian matrix
c     write(IOut,'(2x,a,F12.6)')'detH = ',detH
      vec(1) = -Der11*Dad + CC*x2 + Eac - Ecd                            !-Cft'[x1]*Δad + E_REKS"[x1,x2]*x2 + Eac -Ecd
      vec(2) = -Der12*Dbc + CC*x1 + Ebd - Ecd                            !-Cft'[x2]*Δbc + E_REKS"[x1,x2]*x1 + Ebd -Ecd
      H(1,1) = Der22*Dbc                                                 !Inverse Hessian matrix
      H(2,2) = Der21*Dad                                                 ! |  Cft"[x2]*Δbc  | E_REKS"[x1,x2] |
      H(1,2) = CC                                                        ! | E_REKS"[x1,x2] |  Cft"[x1]*Δad  |/detH
      H(2,1) = H(1,2)
      if(abs(Dad).le.1.d-8.or.abs(Dbc).le.1.d-8)then                     !check if one of the occupations is fixed
         dx1 = zero
         dx2 = zero
         if(abs(H(2,2)).ge.CnvLim)dx1 = vec(1)/H(2,2)
         if(Dad.lt.zero) dx1 = -dx1                                      !inverse the sign if Δad is negative; should be not, however...
         if(abs(H(1,1)).ge.CnvLim)dx2 = vec(2)/H(1,1)
         if(Dbc.lt.zero) dx2 = -dx2
      else                                                               !do both occupations; none is fixed
         dx1 = H(1,1)*vec(1) + H(1,2)*vec(2)                             !the increment vector
         dx2 = H(2,1)*vec(1) + H(2,2)*vec(2)
         if(detH.gt.CnvLim) then                                         !divide by detH
            dx1 = dx1/detH
            dx2 = dx2/detH
         endif
      endif
c     write(IOut,'(2x,a,F12.6,a,F12.6)')'δx1 =',dx1,', δx2 =',dx2
      factor=one                                                         !line search factor
      golden=(one+sqrt(five))/two                                        !golden ratio
      dxnorm = sqrt(dx1*dx1 + dx2*dx2)                                   !norm of the increment vector
      miter=0
      do
      miter=miter+1                                                      !top of line search micro iterations
      x1 = x10 + dx1*factor                                              !new X vector
      x2 = x20 + dx2*factor
      if(x1.gt.one) x1=one                                               !cap x1 from above
      if(x1.lt.zero)x1=zero                                              !cap x1 from below
      if(x2.gt.one) x2=one                                               !cap x2 from above
      if(x2.lt.zero)x2=zero                                              !cap x2 from below
C
      if(abs(Dad).le.1.d-8) x1=1.d0                                      !set x1 to one if Δad is zero
      if(abs(Dbc).le.1.d-8) x1=1.d0                                      !set x2 to one if Δbc is zero
C
      egvb1=egvb(Eab,Eac,Ebd,Ecd,Dad,Dbc,delta,x1,x2)                    !current GVB energy
c     write(IOut,'(2x,a,i4,a,F12.8,a,F20.12)')
c    $     'micro_iter =',miter,', λ =',factor,', eGVB =',egvb1
      if(egvb1.lt.egvb0)goto 30                                          !if energy goes down, continue iterations
      if(abs(factor*dxnorm).lt.CnvLim)goto 30                            !if no lower energy is found, bail out
      if(miter.gt.MaxIt)goto 30
      factor=factor/(golden*golden)                                      !if energy rises, do line search by golden ratio
      if(miter.gt.MaxIt/2)factor=factor/(golden*golden)                  !speed up the division when many micro iters were done
      end do                                                             !loop the line search
   30 egvb0=egvb1
      dx1=x1-x10                                                         !now do the actual displacement
      dx2=x2-x20
      dxnorm = sqrt(dx1*dx1 + dx2*dx2)                                   !norm of the increment vector
C
      if(dbg) write(IOut,'(i4,4x,F10.8,1x,F10.8,4x,F14.12,1x,F20.12)')
     $ iter,x1,x2,dxnorm,egvb0
C
      if(dxnorm.gt.CnvLim) then                                          !convergence test
         if(iter.le.MaxIt) goto 10
         if(iter.gt.MaxIt) then
            if(dbg)  write(IOut,'(2x,a,i4,2x,a)')
     $                         'Not converged within',MaxIt,'iterations'
         endif
      else
         if(dbg)  write(IOut,'(2x,a,i4,2x,a)')
     *          'Converged within',iter,'iterations'
      endif
      EREKS = egvb(Eab,Eac,Ebd,Ecd,Dad,Dbc,delta,x1,x2)                  !do the E_REKS energy
      return
      end subroutine rexslv4x4
C
C © Michael Filatov 2021
      FUNCTION derCft0z(z,delta)
      USE constants, only: one, two, four
      USE prec, ONLY: dp
      IMPLICIT NONE
C
      REAL(KIND=dp) :: DELTA, Z
      REAL(KIND=dp) :: DERCFT0Z
C
      REAL(KIND=dp) :: POW, RAT, VAR
      if(z.le.1.d-10) z = 1.d-10                                         !derivatives diverge at z = 0; cap the z value
      var = four*z
      pow = one-(var+delta)/(two+two*delta)
      derCft0z = var**pow
      return
      end function
      function derCft1z(z,delta)
      USE constants, only: one, two, four
      USE prec, ONLY: dp
      IMPLICIT NONE
C
      REAL(KIND=dp) :: DELTA, Z
C
      REAL(KIND=dp) :: DERCFT1Z
      REAL(KIND=dp) :: POW, RAT, VAR
      if(z.le.1.d-10) z = 1.d-10                                         !derCft1z diverges at z = 0; cap the z values
      var = four*z
      pow = (var+delta)/(two+two*delta)
      rat = (two+delta-var*(one+log(var)))/(one+delta)
      derCft1z = two*rat/var**pow
      return
      end function
      function derCft2z(z,delta)
      USE constants, only: one, two, four
      USE prec, ONLY: dp
      IMPLICIT NONE
C
      REAL(KIND=dp) :: DELTA, Z
C
      REAL(KIND=dp) :: DERCFT2Z
      REAL(KIND=dp) :: POW, RAT, VAR
      if(z.le.1.d-10) z = 1.d-10                                         !derCft2z diverges at z = 0; cap the z values
      var = four*z
      pow = (var+delta)/(two+two*delta)
      rat = ((two+delta-var*(one+log(var)))**2
     $    -two*(one+delta)*(two+delta+var))/(one+delta)**2
      derCft2z = rat/(z*var**pow)
      return
      end function derCft2z
C
C © Michael Filatov 2021
      FUNCTION sgvb01(delta,x1,x2,y1,y2)
      USE constants, only:one, four
      USE prec, ONLY: dp
      IMPLICIT NONE
C
      REAL(KIND=dp) :: delta, x1, x2, y1, y2
      REAL(KIND=dp) :: sgvb01
C
      REAL(KIND=dp) :: dercft0z, dercft1z, dercft2z
      sgvb01 = derCft0z(x1*x2,delta)*derCft0z(y1*y2,delta)/four
     $       + derCft0z((one-x1)*(one-x2),delta)
     $       * derCft0z((one-y1)*(one-y2),delta)/four
      return
      end function
      function dsgvb01dx1(delta,x1,x2,y1,y2)
      USE constants, only:one, four
      USE prec, ONLY: dp
      IMPLICIT NONE
C
      REAL(KIND=dp) :: delta, x1, x2, y1, y2
      REAL(KIND=dp) :: dsgvb01dx1
C
      REAL(KIND=dp) :: dercft0z, dercft1z, dercft2z
      dsgvb01dx1 = derCft1z(x1*x2,delta)*derCft0z(y1*y2,delta)*x2/four
     $ - (one-x2)*derCft0z((one-y1)*(one-y2),delta)
     $ * derCft1z((one-x1)*(one-x2),delta)/four
      return
      end function
      function dsgvb01dx2(delta,x1,x2,y1,y2)
      USE constants, only:one, four
      USE prec, ONLY: dp
      IMPLICIT NONE
C
      REAL(KIND=dp) :: delta, x1, x2, y1, y2
      REAL(KIND=dp) :: dsgvb01dx2
C
      REAL(KIND=dp) :: dercft0z, dercft1z, dercft2z
      dsgvb01dx2 = derCft1z(x1*x2,delta)*derCft0z(y1*y2,delta)*x1/four
     $ - (one-x1)*derCft0z((one-y1)*(one-y2),delta)
     $ * derCft1z((one-x1)*(one-x2),delta)/four
      return
      end function
      function dsgvb01dy1(delta,x1,x2,y1,y2)
      USE constants, only:one, four
      USE prec, ONLY: dp
      IMPLICIT NONE
C
      REAL(KIND=dp) :: delta, x1, x2, y1, y2
      REAL(KIND=dp) :: dsgvb01dy1
C
      REAL(KIND=dp) :: dercft0z, dercft1z, dercft2z
      dsgvb01dy1 = derCft0z(x1*x2,delta)*derCft1z(y1*y2,delta)*y2/four
     $ - (one-y2)*derCft1z((one-y1)*(one-y2),delta)
     $ * derCft0z((one-x1)*(one-x2),delta)/four
      return
      end function
      function dsgvb01dy2(delta,x1,x2,y1,y2)
      USE constants, only:one, four
      USE prec, ONLY: dp
      IMPLICIT NONE
C
      REAL(KIND=dp) :: delta, x1, x2, y1, y2
      REAL(KIND=dp) :: dsgvb01dy2
C
      REAL(KIND=dp) :: dercft0z, dercft1z, dercft2z
      dsgvb01dy2 = derCft0z(x1*x2,delta)*derCft1z(y1*y2,delta)*y1/four
     $ - (one-y1)*derCft1z((one-y1)*(one-y2),delta)
     $ * derCft0z((one-x1)*(one-x2),delta)/four
      return
      end function
      function d2sgvb01dx12(delta,x1,x2,y1,y2)
      USE constants, only:one, four
      USE prec, ONLY: dp
      IMPLICIT NONE
C
      REAL(KIND=dp) :: delta, x1, x2, y1, y2
      REAL(KIND=dp) :: d2sgvb01dx12
C
      REAL(KIND=dp) :: dercft0z, dercft1z, dercft2z
      d2sgvb01dx12 = derCft0z(y1*y2,delta)*derCft2z(x1*x2,delta)
     $ * x2*x2/four
     $ + (one-x2)**2*derCft2z((one-x1)*(one-x2),delta)
     $ * derCft0z((one-y1)*(one-y2),delta)/four
      return
      end function
      function d2sgvb01dx22(delta,x1,x2,y1,y2)
      USE constants, only:one, four
      USE prec, ONLY: dp
      IMPLICIT NONE
C
      REAL(KIND=dp) :: delta, x1, x2, y1, y2
      REAL(KIND=dp) :: d2sgvb01dx22
C
      REAL(KIND=dp) :: dercft0z, dercft1z, dercft2z
      d2sgvb01dx22 = derCft0z(y1*y2,delta)*derCft2z(x1*x2,delta)
     $ * x1*x1/four
     $ + (one-x1)**2*derCft2z((one-x1)*(one-x2),delta)
     $ * derCft0z((one-y1)*(one-y2),delta)/four
      return
      end function
      function d2sgvb01dy12(delta,x1,x2,y1,y2)
      USE constants, only:one, four
      USE prec, ONLY: dp
      IMPLICIT NONE
C
      REAL(KIND=dp) :: delta, x1, x2, y1, y2
      REAL(KIND=dp) :: d2sgvb01dy12
C
      REAL(KIND=dp) :: dercft0z, dercft1z, dercft2z
      d2sgvb01dy12 = derCft0z(x1*x2,delta)*derCft2z(y1*y2,delta)
     $ * y2*y2/four
     $ + (one-y2)**2*derCft2z((one-y1)*(one-y2),delta)
     $ * derCft0z((one-x1)*(one-x2),delta)/four
      return
      end function
      function d2sgvb01dy22(delta,x1,x2,y1,y2)
      USE constants, only:one, four
      USE prec, ONLY: dp
      IMPLICIT NONE
C
      REAL(KIND=dp) :: delta, x1, x2, y1, y2
      REAL(KIND=dp) :: d2sgvb01dy22
C
      REAL(KIND=dp) :: dercft0z, dercft1z, dercft2z
      d2sgvb01dy22 = derCft0z(x1*x2,delta)*derCft2z(y1*y2,delta)
     $ * y1*y1/four
     $ + (one-y1)**2*derCft2z((one-y1)*(one-y2),delta)
     $ * derCft0z((one-x1)*(one-x2),delta)/four
      return
      end function
      function d2sgvb01dx1x2(delta,x1,x2,y1,y2)
      USE constants, only:one, four
      USE prec, ONLY: dp
      IMPLICIT NONE
C
      REAL(KIND=dp) :: delta, x1, x2, y1, y2
      REAL(KIND=dp) :: d2sgvb01dx1x2
C
      REAL(KIND=dp) :: dercft0z, dercft1z, dercft2z
      d2sgvb01dx1x2 = derCft0z((one-y1)*(one-y2),delta)
     $ * (derCft1z((one-x1)*(one-x2),delta)+(one-x1)*(one-x2)
     $ * derCft2z((one-x1)*(one-x2),delta))/four
     $ + derCft0z(y1*y2,delta)*(derCft1z(x1*x2,delta)
     $ + x1*x2*derCft2z(x1*x2,delta))/four
      return
      end function
      function d2sgvb01dy1y2(delta,x1,x2,y1,y2)
      USE constants, only:one, four
      USE prec, ONLY: dp
      IMPLICIT NONE
C
      REAL(KIND=dp) :: delta, x1, x2, y1, y2
      REAL(KIND=dp) :: d2sgvb01dy1y2
C
      REAL(KIND=dp) :: dercft0z, dercft1z, dercft2z
      d2sgvb01dy1y2 = derCft0z((one-x1)*(one-x2),delta)
     $ * (derCft1z((one-y1)*(one-y2),delta)+(one-y1)*(one-y2)
     $ * derCft2z((one-y1)*(one-y2),delta))/four
     $ + derCft0z(x1*x2,delta)*(derCft1z(y1*y2,delta)
     $ + y1*y2*derCft2z(y1*y2,delta))/four
      return
      end function
      function d2sgvb01dx1y1(delta,x1,x2,y1,y2)
      USE constants, only:one, four
      USE prec, ONLY: dp
      IMPLICIT NONE
C
      REAL(KIND=dp) :: delta, x1, x2, y1, y2
      REAL(KIND=dp) :: d2sgvb01dx1y1
C
      REAL(KIND=dp) :: dercft0z, dercft1z, dercft2z
      d2sgvb01dx1y1 =
     $   (one-x2)*(one-y2)*derCft1z((one-x1)*(one-x2),delta)
     $ * derCft1z((one-y1)*(one-y2),delta)/four
     $ + x2*y2*derCft1z(x1*x2,delta)*derCft1z(y1*y2,delta)/four
      return
      end function
      function d2sgvb01dx2y2(delta,x1,x2,y1,y2)
      USE constants, only:one, four
      USE prec, ONLY: dp
      IMPLICIT NONE
C
      REAL(KIND=dp) :: delta, x1, x2, y1, y2
      REAL(KIND=dp) :: d2sgvb01dx2y2
C
      REAL(KIND=dp) :: dercft0z, dercft1z, dercft2z
      d2sgvb01dx2y2 =
     $   (one-x1)*(one-y1)*derCft1z((one-x1)*(one-x2),delta)
     $ * derCft1z((one-y1)*(one-y2),delta)/four
     $ + x1*y1*derCft1z(x1*x2,delta)*derCft1z(y1*y2,delta)/four
      return
      end function
      function d2sgvb01dx1y2(delta,x1,x2,y1,y2)
      USE constants, only:one, four
      USE prec, ONLY: dp
      IMPLICIT NONE
C
      REAL(KIND=dp) :: delta, x1, x2, y1, y2
      REAL(KIND=dp) :: d2sgvb01dx1y2
C
      REAL(KIND=dp) :: dercft0z, dercft1z, dercft2z
      d2sgvb01dx1y2 =
     $   (one-x2)*(one-y1)*derCft1z((one-x1)*(one-x2),delta)
     $ * derCft1z((one-y1)*(one-y2),delta)/four
     $ + x2*y1*derCft1z(x1*x2,delta)*derCft1z(y1*y2,delta)/four
      return
      end function
      function d2sgvb01dx2y1(delta,x1,x2,y1,y2)
      USE constants, only:one, four
      USE prec, ONLY: dp
      IMPLICIT NONE
C
      REAL(KIND=dp) :: delta, x1, x2, y1, y2
      REAL(KIND=dp) :: d2sgvb01dx2y1
C
      REAL(KIND=dp) :: dercft0z, dercft1z, dercft2z
      d2sgvb01dx2y1 =
     $   (one-x1)*(one-y2)*derCft1z((one-x1)*(one-x2),delta)
     $ * derCft1z((one-y1)*(one-y2),delta)/four
     $ + x1*y2*derCft1z(x1*x2,delta)*derCft1z(y1*y2,delta)/four
      return
      end function d2sgvb01dx2y1
C
C GVB(2/4)-REKS energy and derivatives
C © Michael Filatov 2021
      FUNCTION egvb(Eab,Eac,Ebd,Ecd,Dad,Dbc,delta,x1,x2)
      USE constants, only:one, two
      USE prec, ONLY: dp
      IMPLICIT NONE
C
      REAL(KIND=dp) :: Eab, Eac, Ebd, Ecd, Dad, Dbc, delta, x1, x2
      REAL(KIND=dp) :: egvb
      REAL(KIND=dp) :: derCft0z
      egvb =
     $ x1*x2*Eab+x1*(one-x2)*Eac+(one-x1)*x2*Ebd+(one-x1)*(one-x2)*Ecd
     $-derCft0z(x1*(one-x1),delta)*Dad-derCft0z(x2*(one-x2),delta)*Dbc
      return
      end function egvb
      function degvbdx1(Eab,Eac,Ebd,Ecd,Dad,Dbc
     $                                ,delta,x1,x2)
      USE constants, only:one, two
      USE prec, ONLY: dp
      IMPLICIT NONE
C
      REAL(KIND=dp) :: Eab, Eac, Ebd, Ecd, Dad, Dbc, delta, x1, x2
      REAL(KIND=dp) :: degvbdx1
      REAL(KIND=dp) :: derCft1z
      degvbdx1 =
     $ (Eac-Ecd)+(Eab-Eac-Ebd+Ecd)*x2
     $-Dad*(one-two*x1)*derCft1z((one-x1)*x1,delta)
      return
      end function
      function degvbdx2(Eab,Eac,Ebd,Ecd,Dad,Dbc
     $                                ,delta,x1,x2)
      USE constants, only:one, two
      USE prec, ONLY: dp
      IMPLICIT NONE
C
      REAL(KIND=dp) :: Eab, Eac, Ebd, Ecd, Dad, Dbc, delta, x1, x2
      REAL(KIND=dp) :: degvbdx2
      REAL(KIND=dp) :: derCft1z
      degvbdx2 =
     $ (Ebd-Ecd)+(Eab-Eac-Ebd+Ecd)*x1
     $-Dbc*(one-two*x2)*derCft1z((one-x2)*x2,delta)
      return
      end function
      function d2egvbdx12(Eab,Eac,Ebd,Ecd,Dad,Dbc
     $                                ,delta,x1,x2)
      USE constants, only:one, two
      USE prec, ONLY: dp
      IMPLICIT NONE
C
      REAL(KIND=dp) :: Eab, Eac, Ebd, Ecd, Dad, Dbc, delta, x1, x2
      REAL(KIND=dp) :: d2egvbdx12
      REAL(KIND=dp) :: derCft1z, derCft2z
      d2egvbdx12 =
     $ two*Dad*derCft1z((one-x1)*x1,delta)
     $-Dad*(one-two*x1)**2*derCft2z((one-x1)*x1,delta)
      return
      end function
      function d2egvbdx22(Eab,Eac,Ebd,Ecd,Dad,Dbc
     $                                ,delta,x1,x2)
      USE constants, only:one, two
      USE prec, ONLY: dp
      IMPLICIT NONE
C
      REAL(KIND=dp) :: Eab, Eac, Ebd, Ecd, Dad, Dbc, delta, x1, x2
      REAL(KIND=dp) :: d2egvbdx22
      REAL(KIND=dp) :: derCft1z, derCft2z
      d2egvbdx22 =
     $ two*Dbc*derCft1z((one-x2)*x2,delta)
     $-Dbc*(one-two*x2)**2*derCft2z((one-x2)*x2,delta)
      return
      end function
      function d2egvbdx1x2(Eab,Eac,Ebd,Ecd,Dad,Dbc
     $                                ,delta,x1,x2)
      USE constants, only:one, two
      USE prec, ONLY: dp
      IMPLICIT NONE
C
      REAL(KIND=dp) :: Eab, Eac, Ebd, Ecd, Dad, Dbc, delta, x1, x2
      REAL(KIND=dp) :: d2egvbdx1x2
      d2egvbdx1x2 = Eab-Eac-Ebd+Ecd
      return
      end function d2egvbdx1x2
C
C GRVB(2/4)-REKS coupling element and derivatives
C © Michael Filatov 2021
      FUNCTION hgvb01(Eab,Ecd,Dab,Dac,Dad,Dbc,Dbd,Dcd,delta,x1,x2,y1,y2)
      USE constants, only: one, four
      USE prec, ONLY: dp
      IMPLICIT NONE
      REAL(KIND=dp) :: Eab, Ecd, Dab, Dac, Dad, Dbc, Dbd, Dcd, delta,
     &                 x1, x2, y1, y2
      REAL(KIND=dp) :: dercft0z, dercft1z, dercft2z
      REAL(KIND=dp) :: hgvb01
      hgvb01 = (Eab*derCft0z(x1*x2,delta)*derCft0z(y1*y2,delta)
     $+ Ecd*derCft0z((one-x1)*(one-x2),delta)
     $* derCft0z((one-y1)*(one-y2),delta))/four
     $+ (Dab*(derCft0z((one-x1)*x2,delta)*derCft0z(y1*(one-y2),delta)
     $+ derCft0z(x1*(one-x2),delta)*derCft0z((one-y1)*y2,delta)))/four
     $- (Dac*((derCft0z((one-x1)*(one-x2),delta)
     $* derCft0z(y1*(one-y2),delta))
     $+ derCft0z(x1*x2,delta)*derCft0z((one-y1)*y2,delta)))/four
     $- (Dad*((derCft0z(x1*(one-x2),delta)
     $* derCft0z((one-y1)*(one-y2),delta))
     $+ derCft0z((one-x1)*x2,delta)*derCft0z(y1*y2,delta)))/four
     $- (Dbc*((derCft0z((one-x1)*x2,delta)
     $* derCft0z((one-y1)*(one-y2),delta))
     $+ derCft0z(x1*(one-x2),delta)*derCft0z(y1*y2,delta)))/four
     $- (Dbd*((derCft0z(x1*x2,delta)*derCft0z(y1*(one-y2),delta))
     $+ derCft0z((one-x1)*(one-x2),delta)
     $* derCft0z((one-y1)*y2,delta)))/four
     $+ (Dcd*(derCft0z(x1*(one-x2),delta)*derCft0z(y1*(one-y2),delta)
     $+ derCft0z((one-x1)*x2,delta)*derCft0z((one-y1)*y2,delta)))/four
      return
      end function
      function dhgvb01dx1(Eab,Ecd,Dab,Dac,Dad,Dbc,Dbd,Dcd
     $                ,delta,x1,x2,y1,y2)
      USE constants, only: one, four
      USE prec, ONLY: dp
      IMPLICIT NONE
      REAL(KIND=dp) :: Eab, Ecd, Dab, Dac, Dad, Dbc, Dbd, Dcd, delta,
     &                 x1, x2, y1, y2
      REAL(KIND=dp) :: dercft0z, dercft1z, dercft2z
      REAL(KIND=dp) :: dhgvb01dx1
      dhgvb01dx1 = (Eab*x2*derCft1z(x1*x2,delta)*derCft0z(y1*y2,delta)
     $- Ecd*(one-x2)*derCft1z((one-x1)*(one-x2),delta)
     $* derCft0z((one-y1)*(one-y2),delta))/four
     $+ Dab*((one-x2)*derCft1z(x1*(one-x2),delta)
     $* derCft0z((one-y1)*y2,delta)
     $- x2*derCft1z((one-x1)*x2,delta)*derCft0z(y1*(one-y2),delta))/four
     $- Dac*(-(one-x2)*derCft1z((one-x1)*(one-x2),delta)
     $* derCft0z(y1*(one-y2),delta)
     $+ x2*derCft1z(x1*x2,delta)*derCft0z((one-y1)*y2,delta))/four
     $- Dad*((one-x2)*derCft1z(x1*(one-x2),delta)
     $* derCft0z((one-y1)*(one-y2),delta)
     $- x2*derCft1z((one-x1)*x2,delta)*derCft0z(y1*y2,delta))/four
     $- Dbc*(-x2*derCft1z((one-x1)*x2,delta)
     $* derCft0z((one-y1)*(one-y2),delta)
     $+ (one-x2)*derCft1z(x1*(one-x2),delta)*derCft0z(y1*y2,delta))/four
     $- Dbd*(x2*derCft1z(x1*x2,delta)*derCft0z(y1*(one-y2),delta)
     $- (one-x2)*derCft1z((one-x1)*(one-x2),delta)
     $* derCft0z((one-y1)*y2,delta))/four
     $+ Dcd*((one-x2)*derCft1z(x1*(one-x2),delta)
     $* derCft0z(y1*(one-y2),delta)
     $- x2*derCft1z((one-x1)*x2,delta)*derCft0z((one-y1)*y2,delta))/four
      return
      end function
      function dhgvb01dx2(Eab,Ecd,Dab,Dac,Dad,Dbc,Dbd,Dcd
     $                ,delta,x1,x2,y1,y2)
      USE constants, only: one, four
      USE prec, ONLY: dp
      IMPLICIT NONE
      REAL(KIND=dp) :: Eab, Ecd, Dab, Dac, Dad, Dbc, Dbd, Dcd, delta,
     &                 x1, x2, y1, y2
      REAL(KIND=dp) :: dercft0z, dercft1z, dercft2z
      REAL(KIND=dp) :: dhgvb01dx2
      dhgvb01dx2 = (Eab*x1*derCft1z(x1*x2,delta)*derCft0z(y1*y2,delta)
     $- Ecd*(one-x1)*derCft1z((one-x1)*(one-x2),delta)
     $* derCft0z((one-y1)*(one-y2),delta))/four
     $+ Dab*(-x1*derCft1z(x1*(one-x2),delta)*derCft0z((one-y1)*y2,delta)
     $+ (one-x1)*derCft1z((one-x1)*x2,delta)*derCft0z(y1*(one-y2),delta)
     $  )/four
     $- Dac*(-(one-x1)*derCft1z((one-x1)*(one-x2),delta)
     $* derCft0z(y1*(one-y2),delta)
     $+ x1*derCft1z(x1*x2,delta)*derCft0z((one-y1)*y2,delta))/four
     $- Dad*(-x1*derCft1z(x1*(one-x2),delta)
     $* derCft0z((one-y1)*(one-y2),delta)
     $+ (one-x1)*derCft1z((one-x1)*x2,delta)*derCft0z(y1*y2,delta))/four
     $- Dbc*((one-x1)*derCft1z((one-x1)*x2,delta)
     $* derCft0z((one-y1)*(one-y2),delta)
     $- x1*derCft1z(x1*(one-x2),delta)*derCft0z(y1*y2,delta))/four
     $- Dbd*(x1*derCft1z(x1*x2,delta)*derCft0z(y1*(one-y2),delta)
     $- (one-x1)*derCft1z((one-x1)*(one-x2),delta)
     $* derCft0z((one-y1)*y2,delta))/four
     $+ Dcd*(-x1*derCft1z(x1*(one-x2),delta)*derCft0z(y1*(one-y2),delta)
     $+ (one-x1)*derCft1z((one-x1)*x2,delta)*derCft0z((one-y1)*y2,delta)
     $  )/four
      return
      end function
      function dhgvb01dy1(Eab,Ecd,Dab,Dac,Dad,Dbc,Dbd,Dcd
     $                ,delta,x1,x2,y1,y2)
      USE constants, only: one, four
      USE prec, ONLY: dp
      IMPLICIT NONE
      REAL(KIND=dp) :: Eab, Ecd, Dab, Dac, Dad, Dbc, Dbd, Dcd, delta,
     &                 x1, x2, y1, y2
      REAL(KIND=dp) :: dercft0z, dercft1z, dercft2z
      REAL(KIND=dp) :: dhgvb01dy1
      dhgvb01dy1 = (Eab*y2*derCft0z(x1*x2,delta)*derCft1z(y1*y2,delta)
     $- Ecd*(one-y2)*derCft0z((one-x1)*(one-x2),delta)
     $* derCft1z((one-y1)*(one-y2),delta))/four
     $+ Dab*(-y2*derCft0z(x1*(one-x2),delta)*derCft1z((one-y1)*y2,delta)
     $+ (one-y2)*derCft0z((one-x1)*x2,delta)*derCft1z(y1*(one-y2),delta)
     $  )/four
     $- Dac*((one-y2)*derCft0z((one-x1)*(one-x2),delta)
     $* derCft1z(y1*(one-y2),delta)
     $- y2*derCft0z(x1*x2,delta)*derCft1z((one-y1)*y2,delta))/four
     $- Dad*(-(one-y2)*derCft0z(x1*(one-x2),delta)
     $* derCft1z((one-y1)*(one-y2),delta)
     $+ y2*derCft0z((one-x1)*x2,delta)*derCft1z(y1*y2,delta))/four
     $- Dbc*(-(one-y2)*derCft0z((one-x1)*x2,delta)
     $* derCft1z((one-y1)*(one-y2),delta)
     $+ y2*derCft0z(x1*(one-x2),delta)*derCft1z(y1*y2,delta))/four
     $- Dbd*((one-y2)*derCft0z(x1*x2,delta)*derCft1z(y1*(one-y2),delta)
     $- y2*derCft0z((one-x1)*(one-x2),delta)*derCft1z((one-y1)*y2,delta)
     $  )/four
     $+ Dcd*((one-y2)*derCft0z(x1*(one-x2),delta)
     $* derCft1z(y1*(one-y2),delta)
     $- y2*derCft0z((one-x1)*x2,delta)*derCft1z((one-y1)*y2,delta))/four
      return
      end function
      function dhgvb01dy2(Eab,Ecd,Dab,Dac,Dad,Dbc,Dbd,Dcd
     $                ,delta,x1,x2,y1,y2)
      USE constants, only: one, four
      USE prec, ONLY: dp
      IMPLICIT NONE
      REAL(KIND=dp) :: Eab, Ecd, Dab, Dac, Dad, Dbc, Dbd, Dcd, delta,
     &                 x1, x2, y1, y2
      REAL(KIND=dp) :: dercft0z, dercft1z, dercft2z
      REAL(KIND=dp) :: dhgvb01dy2
      dhgvb01dy2 = (Eab*y1*derCft0z(x1*x2,delta)*derCft1z(y1*y2,delta)
     $- Ecd*(one-y1)*derCft0z((one-x1)*(one-x2),delta)
     $* derCft1z((one-y1)*(one-y2),delta))/four
     $+ Dab*((one-y1)*derCft0z(x1*(one-x2),delta)
     $* derCft1z((one-y1)*y2,delta)
     $- y1*derCft0z((one-x1)*x2,delta)*derCft1z(y1*(one-y2),delta)
     $  )/four
     $- Dac*(-y1*derCft0z((one-x1)*(one-x2),delta)
     $* derCft1z(y1*(one-y2),delta)
     $+ (one-y1)*derCft0z(x1*x2,delta)*derCft1z((one-y1)*y2,delta))/four
     $- Dad*(-(one-y1)*derCft0z(x1*(one-x2),delta)
     $* derCft1z((one-y1)*(one-y2),delta)
     $+ y1*derCft0z((one-x1)*x2,delta)*derCft1z(y1*y2,delta))/four
     $- Dbc*(-(one-y1)*derCft0z((one-x1)*x2,delta)
     $* derCft1z((one-y1)*(one-y2),delta)
     $+ y1*derCft0z(x1*(one-x2),delta)*derCft1z(y1*y2,delta))/four
     $- Dbd*(-y1*derCft0z(x1*x2,delta)*derCft1z(y1*(one-y2),delta)
     $+ (one-y1)*derCft0z((one-x1)*(one-x2),delta)
     $* derCft1z((one-y1)*y2,delta))/four
     $+ Dcd*(-y1*derCft0z(x1*(one-x2),delta)*derCft1z(y1*(one-y2),delta)
     $+ (one-y1)*derCft0z((one-x1)*x2,delta)*derCft1z((one-y1)*y2,delta)
     $  )/four
      return
      end function
      function d2hgvb01dx1y1(Eab,Ecd,Dab,Dac,Dad,Dbc,Dbd,Dcd
     $                ,delta,x1,x2,y1,y2)
      USE constants, only: one, four
      USE prec, ONLY: dp
      IMPLICIT NONE
      REAL(KIND=dp) :: Eab, Ecd, Dab, Dac, Dad, Dbc, Dbd, Dcd, delta,
     &                 x1, x2, y1, y2
      REAL(KIND=dp) :: dercft0z, dercft1z, dercft2z
      REAL(KIND=dp) :: d2hgvb01dx1y1
      d2hgvb01dx1y1 = (Eab*x2*y2*derCft1z(x1*x2,delta)
     $* derCft1z(y1*y2,delta)
     $+ Ecd*(one-x2)*(one-y2)*derCft1z((one-x1)*(one-x2),delta)
     $* derCft1z((one-y1)*(one-y2),delta))/four
     $+ Dab*(-(one-x2)*y2*derCft1z(x1*(one-x2),delta)
     $* derCft1z((one-y1)*y2,delta)
     $- x2*(one-y2)*derCft1z((one-x1)*x2,delta)
     $* derCft1z(y1*(one-y2),delta))/four
     $- Dac*(-(one-x2)*(one-y2)*derCft1z((one-x1)*(one-x2),delta)
     $* derCft1z(y1*(one-y2),delta)
     $- x2*y2*derCft1z(x1*x2,delta)*derCft1z((one-y1)*y2,delta))/four
     $- Dad*(-(one-x2)*(one-y2)*derCft1z(x1*(one-x2),delta)
     $* derCft1z((one-y1)*(one-y2),delta)
     $- x2*y2*derCft1z((one-x1)*x2,delta)*derCft1z(y1*y2,delta))/four
     $- Dbc*(x2*(one-y2)*derCft1z((one-x1)*x2,delta)
     $* derCft1z((one-y1)*(one-y2),delta)
     $+ (one-x2)*y2*derCft1z(x1*(one-x2),delta)
     $* derCft1z(y1*y2,delta))/four
     $- Dbd*(x2*(one-y2)*derCft1z(x1*x2,delta)
     $* derCft1z(y1*(one-y2),delta)
     $+ (one-x2)*y2*derCft1z((one-x1)*(one-x2),delta)
     $* derCft1z((one-y1)*y2,delta))/four
     $+ Dcd*((one-x2)*(one-y2)*derCft1z(x1*(one-x2),delta)
     $* derCft1z(y1*(one-y2),delta)
     $+ x2*y2*derCft1z((one-x1)*x2,delta)*derCft1z((one-y1)*y2,delta)
     $  )/four
      return
      end function
      function d2hgvb01dx1y2(Eab,Ecd,Dab,Dac,Dad,Dbc,Dbd,Dcd
     $                ,delta,x1,x2,y1,y2)
      USE constants, only: one, four
      USE prec, ONLY: dp
      IMPLICIT NONE
      REAL(KIND=dp) :: Eab, Ecd, Dab, Dac, Dad, Dbc, Dbd, Dcd, delta,
     &                 x1, x2, y1, y2
      REAL(KIND=dp) :: dercft0z, dercft1z, dercft2z
      REAL(KIND=dp) :: d2hgvb01dx1y2
      d2hgvb01dx1y2 = (Eab*x2*y1*derCft1z(x1*x2,delta)
     $* derCft1z(y1*y2,delta)
     $+ Ecd*(one-x2)*(one-y1)*derCft1z((one-x1)*(one-x2),delta)
     $* derCft1z((one-y1)*(one-y2),delta))/four
     $+ Dab*((one-x2)*(one-y1)*derCft1z(x1*(one-x2),delta)
     $* derCft1z((one-y1)*y2,delta)
     $+ x2*y1*derCft1z((one-x1)*x2,delta)
     $* derCft1z(y1*(one-y2),delta))/four
     $- Dac*((one-x2)*y1*derCft1z((one-x1)*(one-x2),delta)
     $* derCft1z(y1*(one-y2),delta)
     $+ x2*(one-y1)*derCft1z(x1*x2,delta)*derCft1z((one-y1)*y2,delta)
     $  )/four
     $- Dad*(-(one-x2)*(one-y1)*derCft1z(x1*(one-x2),delta)
     $* derCft1z((one-y1)*(one-y2),delta)
     $- x2*y1*derCft1z((one-x1)*x2,delta)*derCft1z(y1*y2,delta))/four
     $- Dbc*(x2*(one-y1)*derCft1z((one-x1)*x2,delta)
     $* derCft1z((one-y1)*(one-y2),delta)
     $+ (one-x2)*y1*derCft1z(x1*(one-x2),delta)
     $* derCft1z(y1*y2,delta))/four
     $- Dbd*(-x2*y1*derCft1z(x1*x2,delta)
     $* derCft1z(y1*(one-y2),delta)
     $- (one-x2)*(one-y1)*derCft1z((one-x1)*(one-x2),delta)
     $* derCft1z((one-y1)*y2,delta))/four
     $+ Dcd*(-(one-x2)*y1*derCft1z(x1*(one-x2),delta)
     $* derCft1z(y1*(one-y2),delta)
     $- x2*(one-y1)*derCft1z((one-x1)*x2,delta)
     $* derCft1z((one-y1)*y2,delta))/four
      return
      end function
      function d2hgvb01dx2y1(Eab,Ecd,Dab,Dac,Dad,Dbc,Dbd,Dcd
     $                ,delta,x1,x2,y1,y2)
      USE constants, only: one, four
      USE prec, ONLY: dp
      IMPLICIT NONE
      REAL(KIND=dp) :: Eab, Ecd, Dab, Dac, Dad, Dbc, Dbd, Dcd, delta,
     &                 x1, x2, y1, y2
      REAL(KIND=dp) :: dercft0z, dercft1z, dercft2z
      REAL(KIND=dp) :: d2hgvb01dx2y1
      d2hgvb01dx2y1 = (Eab*x1*y2*derCft1z(x1*x2,delta)
     $* derCft1z(y1*y2,delta)
     $+ Ecd*(one-x1)*(one-y2)*derCft1z((one-x1)*(one-x2),delta)
     $* derCft1z((one-y1)*(one-y2),delta))/four
     $+ Dab*(x1*y2*derCft1z(x1*(one-x2),delta)
     $* derCft1z((one-y1)*y2,delta)
     $+ (one-x1)*(one-y2)*derCft1z((one-x1)*x2,delta)
     $* derCft1z(y1*(one-y2),delta))/four
     $- Dac*(-(one-x1)*(one-y2)*derCft1z((one-x1)*(one-x2),delta)
     $* derCft1z(y1*(one-y2),delta)
     $- x1*y2*derCft1z(x1*x2,delta)*derCft1z((one-y1)*y2,delta))/four
     $- Dad*(x1*(one-y2)*derCft1z(x1*(one-x2),delta)
     $* derCft1z((one-y1)*(one-y2),delta)
     $+ (one-x1)*y2*derCft1z((one-x1)*x2,delta)*derCft1z(y1*y2,delta)
     $  )/four
     $- Dbc*(-(one-x1)*(one-y2)*derCft1z((one-x1)*x2,delta)
     $* derCft1z((one-y1)*(one-y2),delta)
     $- x1*y2*derCft1z(x1*(one-x2),delta)*derCft1z(y1*y2,delta))/four
     $- Dbd*(x1*(one-y2)*derCft1z(x1*x2,delta)
     $* derCft1z(y1*(one-y2),delta)
     $+ (one-x1)*y2*derCft1z((one-x1)*(one-x2),delta)
     $* derCft1z((one-y1)*y2,delta))/four
     $+ Dcd*(-x1*(one-y2)*derCft1z(x1*(one-x2),delta)
     $* derCft1z(y1*(one-y2),delta)
     $- (one-x1)*y2*derCft1z((one-x1)*x2,delta)
     $* derCft1z((one-y1)*y2,delta))/four
      return
      end function
      function d2hgvb01dx2y2(Eab,Ecd,Dab,Dac,Dad,Dbc,Dbd,Dcd
     $                ,delta,x1,x2,y1,y2)
      USE constants, only: one, four
      USE prec, ONLY: dp
      IMPLICIT NONE
      REAL(KIND=dp) :: Eab, Ecd, Dab, Dac, Dad, Dbc, Dbd, Dcd, delta,
     &                 x1, x2, y1, y2
      REAL(KIND=dp) :: dercft0z, dercft1z, dercft2z
      REAL(KIND=dp) :: d2hgvb01dx2y2
      d2hgvb01dx2y2 = (Eab*x1*y1*derCft1z(x1*x2,delta)
     $* derCft1z(y1*y2,delta)
     $+ Ecd*(one-x1)*(one-y1)*derCft1z((one-x1)*(one-x2),delta)
     $* derCft1z((one-y1)*(one-y2),delta))/four
     $+ Dab*(-x1*(one-y1)*derCft1z(x1*(one-x2),delta)
     $* derCft1z((one-y1)*y2,delta)
     $- (one-x1)*y1*derCft1z((one-x1)*x2,delta)
     $* derCft1z(y1*(one-y2),delta))/four
     $- Dac*((one-x1)*y1*derCft1z((one-x1)*(one-x2),delta)
     $* derCft1z(y1*(one-y2),delta)
     $+ x1*(one-y1)*derCft1z(x1*x2,delta)*derCft1z((one-y1)*y2,delta)
     $  )/four
     $- Dad*(x1*(one-y1)*derCft1z(x1*(one-x2),delta)
     $* derCft1z((one-y1)*(one-y2),delta)
     $+ (one-x1)*y1*derCft1z((one-x1)*x2,delta)*derCft1z(y1*y2,delta)
     $  )/four
     $- Dbc*(-(one-x1)*(one-y1)*derCft1z((one-x1)*x2,delta)
     $* derCft1z((one-y1)*(one-y2),delta)
     $- x1*y1*derCft1z(x1*(one-x2),delta)*derCft1z(y1*y2,delta))/four
     $- Dbd*(-x1*y1*derCft1z(x1*x2,delta)*derCft1z(y1*(one-y2),delta)
     $- (one-x1)*(one-y1)*derCft1z((one-x1)*(one-x2),delta)
     $* derCft1z((one-y1)*y2,delta))/four
     $+ Dcd*(x1*y1*derCft1z(x1*(one-x2),delta)
     $* derCft1z(y1*(one-y2),delta)
     $+ (one-x1)*(one-y1)*derCft1z((one-x1)*x2,delta)
     $* derCft1z((one-y1)*y2,delta))/four
      return
      end function
      function d2hgvb01dx1x2(Eab,Ecd,Dab,Dac,Dad,Dbc,Dbd,Dcd
     $                ,delta,x1,x2,y1,y2)
      USE constants, only: one, four
      USE prec, ONLY: dp
      IMPLICIT NONE
      REAL(KIND=dp) :: Eab, Ecd, Dab, Dac, Dad, Dbc, Dbd, Dcd, delta,
     &                 x1, x2, y1, y2
      REAL(KIND=dp) :: dercft0z, dercft1z, dercft2z
      REAL(KIND=dp) :: d2hgvb01dx1x2
      d2hgvb01dx1x2 =
     $  (Ecd*(derCft0z((one-y1)*(one-y2),delta)
     $*derCft1z((one-x1)*(one-x2),delta)
     $+ (one-x1)*(one-x2)*derCft0z((one-y1)*(one-y2),delta)
     $*derCft2z((one-x1)*(one-x2),delta)))/four +
     $  (Dbc*(derCft0z(y1*y2,delta)*(derCft1z(x1*(one-x2),delta)
     $+ x1*(one-x2)*derCft2z(x1*(one-x2),delta)) +
     $  derCft0z((one-y1)*(one-y2),delta)*(derCft1z((one-x1)*x2,delta)
     $+ (one-x1)*x2*derCft2z((one-x1)*x2,delta))))/four +
     $  (Dab*(-derCft0z((one-y1)*y2,delta)*(derCft1z(x1*(one-x2),delta)
     $+ x1*(one-x2)*derCft2z(x1*(one-x2),delta)) -
     $  derCft0z(y1*(one-y2),delta)*(derCft1z((one-x1)*x2,delta)
     $+ (one-x1)*x2*derCft2z((one-x1)*x2,delta))))/four +
     $  (Dcd*(-derCft0z(y1*(one-y2),delta)*(derCft1z(x1*(one-x2),delta)
     $+ x1*(one-x2)*derCft2z(x1*(one-x2),delta)) -
     $  derCft0z((one-y1)*y2,delta)*(derCft1z((one-x1)*x2,delta)
     $+ (one-x1)*x2*derCft2z((one-x1)*x2,delta))))/four +
     $  (Dad*(derCft0z((one-y1)*(one-y2),delta)
     $*(derCft1z(x1*(one-x2),delta)
     $+ x1*(one-x2)*derCft2z(x1*(one-x2),delta)) +
     $  derCft0z(y1*y2,delta)*(derCft1z((one-x1)*x2,delta)
     $+ (one-x1)*x2*derCft2z((one-x1)*x2,delta))))/four +
     $  (Eab*(derCft0z(y1*y2,delta)*derCft1z(x1*x2,delta)
     $+ x1*x2*derCft0z(y1*y2,delta)*derCft2z(x1*x2,delta)))/four +
     $  (Dbd*(-(derCft0z((one-y1)*y2,delta)
     $*(derCft1z((one-x1)*(one-x2),delta)
     $+ (one-x1)*(one-x2)*derCft2z((one-x1)*(one-x2),delta))) -
     $  derCft0z(y1*(one-y2),delta)*(derCft1z(x1*x2,delta)
     $+ x1*x2*derCft2z(x1*x2,delta))))/four +
     $  (Dac*(-derCft0z(y1*(one-y2),delta)
     $*(derCft1z((one-x1)*(one-x2),delta)
     $+ (one-x1)*(one-x2)*derCft2z((one-x1)*(one-x2),delta)) -
     $  derCft0z((one-y1)*y2,delta)*(derCft1z(x1*x2,delta)
     $+ x1*x2*derCft2z(x1*x2,delta))))/four
      return
      end function
      function d2hgvb01dy1y2(Eab,Ecd,Dab,Dac,Dad,Dbc,Dbd,Dcd
     $                ,delta,x1,x2,y1,y2)
      USE constants, only: one, four
      USE prec, ONLY: dp
      IMPLICIT NONE
      REAL(KIND=dp) :: Eab, Ecd, Dab, Dac, Dad, Dbc, Dbd, Dcd, delta,
     &                 x1, x2, y1, y2
      REAL(KIND=dp) :: dercft0z, dercft1z, dercft2z
      REAL(KIND=dp) :: d2hgvb01dy1y2
      d2hgvb01dy1y2 =
     $ (Ecd*(derCft0z((one-x1)*(one-x2),delta)
     $*derCft1z((one-y1)*(one-y2),delta)
     $+ (one-y1)*(one-y2)*derCft0z((one-x1)*(one-x2),delta)
     $*derCft2z((one-y1)*(one-y2),delta)))/four +
     $  (Dbd*(derCft0z(x1*x2,delta)*(derCft1z(y1*(one-y2),delta)
     $+ y1*(one-y2)*derCft2z(y1*(one-y2),delta)) +
     $  derCft0z((one-x1)*(one-x2),delta)*(derCft1z((one-y1)*y2,delta)
     $+ (one-y1)*y2*derCft2z((one-y1)*y2,delta))))/four +
     $  (Dab*(-(derCft0z((one-x1)*x2,delta)*(derCft1z(y1*(one-y2),delta)
     $+ y1*(one-y2)*derCft2z(y1*(one-y2),delta))) -
     $  derCft0z(x1*(one-x2),delta)*(derCft1z((one-y1)*y2,delta)
     $+ (one-y1)*y2*derCft2z((one-y1)*y2,delta))))/four +
     $  (Dcd*(-(derCft0z(x1*(one-x2),delta)*(derCft1z(y1*(one-y2),delta)
     $+ y1*(one-y2)*derCft2z(y1*(one-y2),delta))) -
     $  derCft0z((one-x1)*x2,delta)*(derCft1z((one-y1)*y2,delta)
     $+ (one-y1)*y2*derCft2z((one-y1)*y2,delta))))/four +
     $ (Dac*(derCft0z((one-x1)*(one-x2),delta)
     $*(derCft1z(y1*(one-y2),delta)
     $+ y1*(one-y2)*derCft2z(y1*(one-y2),delta)) +
     $  derCft0z(x1*x2,delta)*(derCft1z((one-y1)*y2,delta)
     $+ (one-y1)*y2*derCft2z((one-y1)*y2,delta))))/four +
     $  (Eab*(derCft0z(x1*x2,delta)*derCft1z(y1*y2,delta)
     $+ y1*y2*derCft0z(x1*x2,delta)*derCft2z(y1*y2,delta)))/four +
     $ (Dbc*(-(derCft0z((one-x1)*x2,delta)
     $*(derCft1z((one-y1)*(one-y2),delta)
     $+ (one-y1)*(one-y2)*derCft2z((one-y1)*(one-y2),delta))) -
     $  derCft0z(x1*(one-x2),delta)*(derCft1z(y1*y2,delta)
     $+ y1*y2*derCft2z(y1*y2,delta))))/four +
     $ (Dad*(-(derCft0z(x1*(one-x2),delta)
     $*(derCft1z((one-y1)*(one-y2),delta)
     $+ (one-y1)*(one-y2)*derCft2z((one-y1)*(one-y2),delta))) -
     $  derCft0z((one-x1)*x2,delta)*(derCft1z(y1*y2,delta)
     $+ y1*y2*derCft2z(y1*y2,delta))))/four
      return
      end function
      function d2hgvb01dx12(Eab,Ecd,Dab,Dac,Dad,Dbc,Dbd,Dcd
     $                ,delta,x1,x2,y1,y2)
      USE constants, only: one, four
      USE prec, ONLY: dp
      IMPLICIT NONE
      REAL(KIND=dp) :: Eab, Ecd, Dab, Dac, Dad, Dbc, Dbd, Dcd, delta,
     &                 x1, x2, y1, y2
      REAL(KIND=dp) :: dercft0z, dercft1z, dercft2z
      REAL(KIND=dp) :: d2hgvb01dx12
      d2hgvb01dx12 =
     $  (Eab*x2**2*derCft0z(y1*y2,delta)*derCft2z(x1*x2,delta))/four +
     $ (Ecd*(one-x2)**2*derCft0z((one-y1)*(one-y2),delta)
     $*derCft2z((one-x1)*(one-x2),delta))/four +
     $ (Dbc*(-((one-x2)**2*derCft0z(y1*y2,delta)
     $*derCft2z(x1*(one-x2),delta)) -
     $x2**2*derCft0z((one-y1)*(one-y2),delta)
     $*derCft2z((one-x1)*x2,delta)))/four +
     $ (Dab*((one-x2)**2*derCft0z((one-y1)*y2,delta)
     $*derCft2z(x1*(one-x2),delta) + x2**2*derCft0z(y1*(one-y2),delta)
     $*derCft2z((one-x1)*x2,delta)))/four +
     $ (Dcd*((one-x2)**2*derCft0z(y1*(one-y2),delta)
     $*derCft2z(x1*(one-x2),delta) + x2**2*derCft0z((one-y1)*y2,delta)
     $*derCft2z((one-x1)*x2,delta)))/four +
     $ (Dad*(-((one-x2)**2*derCft0z((one-y1)*(one-y2),delta)
     $*derCft2z(x1*(one-x2),delta)) - x2**2*derCft0z(y1*y2,delta)
     $*derCft2z((one-x1)*x2,delta)))/four +
     $ (Dbd*(-((one-x2)**2*derCft0z((one-y1)*y2,delta)
     $*derCft2z((one-x1)*(one-x2),delta)) -
     $x2**2*derCft0z(y1*(one-y2),delta)*derCft2z(x1*x2,delta)))/four +
     $ (Dac*(-((one-x2)**2*derCft0z(y1*(one-y2),delta)
     $*derCft2z((one-x1)*(one-x2),delta)) -
     $x2**2*derCft0z((one-y1)*y2,delta)*derCft2z(x1*x2,delta)))/four
      return
      end function
      function d2hgvb01dx22(Eab,Ecd,Dab,Dac,Dad,Dbc,Dbd,Dcd
     $                ,delta,x1,x2,y1,y2)
      USE constants, only: one, four
      USE prec, ONLY: dp
      IMPLICIT NONE
      REAL(KIND=dp) :: Eab, Ecd, Dab, Dac, Dad, Dbc, Dbd, Dcd, delta,
     &                 x1, x2, y1, y2
      REAL(KIND=dp) :: dercft0z, dercft1z, dercft2z
      REAL(KIND=dp) :: d2hgvb01dx22
      d2hgvb01dx22 =
     $  (Eab*x1**2*derCft0z(y1*y2,delta)*derCft2z(x1*x2,delta))/four +
     $ (Ecd*(one-x1)**2*derCft0z((one-y1)*(one-y2),delta)
     $*derCft2z((one-x1)*(one-x2),delta))/four +
     $  (Dbc*(-(x1**2*derCft0z(y1*y2,delta)
     $*derCft2z(x1*(one-x2),delta)) -
     $(one-x1)**2*derCft0z((one-y1)*(one-y2),delta)
     $*derCft2z((one-x1)*x2,delta)))/four +
     $ (Dab*(x1**2*derCft0z((one-y1)*y2,delta)
     $*derCft2z(x1*(one-x2),delta) +
     $(one-x1)**2*derCft0z(y1*(one-y2),delta)
     $*derCft2z((one-x1)*x2,delta)))/four +
     $ (Dcd*(x1**2*derCft0z(y1*(one-y2),delta)
     $*derCft2z(x1*(one-x2),delta) +
     $(one-x1)**2*derCft0z((one-y1)*y2,delta)
     $*derCft2z((one-x1)*x2,delta)))/four +
     $ (Dad*(-(x1**2*derCft0z((one-y1)*(one-y2),delta)
     $*derCft2z(x1*(one-x2),delta)) -
     $(one-x1)**2*derCft0z(y1*y2,delta)
     $*derCft2z((one-x1)*x2,delta)))/four +
     $ (Dbd*(-((one-x1)**2*derCft0z((one-y1)*y2,delta)
     $*derCft2z((one-x1)*(one-x2),delta)) -
     $x1**2*derCft0z(y1*(one-y2),delta)*derCft2z(x1*x2,delta)))/four +
     $ (Dac*(-((one-x1)**2*derCft0z(y1*(one-y2),delta)
     $*derCft2z((one-x1)*(one-x2),delta)) -
     $x1**2*derCft0z((one-y1)*y2,delta)*derCft2z(x1*x2,delta)))/four
      return
      end function
      function d2hgvb01dy12(Eab,Ecd,Dab,Dac,Dad,Dbc,Dbd,Dcd
     $                ,delta,x1,x2,y1,y2)
      USE constants, only: one, four
      USE prec, ONLY: dp
      IMPLICIT NONE
      REAL(KIND=dp) :: Eab, Ecd, Dab, Dac, Dad, Dbc, Dbd, Dcd, delta,
     &                 x1, x2, y1, y2
      REAL(KIND=dp) :: dercft0z, dercft1z, dercft2z
      REAL(KIND=dp) :: d2hgvb01dy12
      d2hgvb01dy12 =
     $  (Eab*y2**2*derCft0z(x1*x2,delta)*derCft2z(y1*y2,delta))/four +
     $ (Ecd*(one-y2)**2*derCft0z((one-x1)*(one-x2),delta)
     $*derCft2z((one-y1)*(one-y2),delta))/four +
     $ (Dbd*(-((one-y2)**2*derCft0z(x1*x2,delta)
     $*derCft2z(y1*(one-y2),delta)) -
     $y2**2*derCft0z((one-x1)*(one-x2),delta)
     $*derCft2z((one-y1)*y2,delta)))/four +
     $ (Dab*((one-y2)**2*derCft0z((one-x1)*x2,delta)
     $*derCft2z(y1*(one-y2),delta) +
     $y2**2*derCft0z(x1*(one-x2),delta)
     $*derCft2z((one-y1)*y2,delta)))/four +
     $ (Dcd*((one-y2)**2*derCft0z(x1*(one-x2),delta)
     $*derCft2z(y1*(one-y2),delta) +
     $y2**2*derCft0z((one-x1)*x2,delta)
     $*derCft2z((one-y1)*y2,delta)))/four +
     $ (Dac*(-((one-y2)**2*derCft0z((one-x1)*(one-x2),delta)
     $*derCft2z(y1*(one-y2),delta)) -
     $y2**2*derCft0z(x1*x2,delta)*derCft2z((one-y1)*y2,delta)))/four +
     $ (Dbc*(-((one-y2)**2*derCft0z((one-x1)*x2,delta)
     $*derCft2z((one-y1)*(one-y2),delta)) -
     $y2**2*derCft0z(x1*(one-x2),delta)*derCft2z(y1*y2,delta)))/four +
     $ (Dad*(-((one-y2)**2*derCft0z(x1*(one-x2),delta)
     $*derCft2z((one-y1)*(one-y2),delta)) -
     $y2**2*derCft0z((one-x1)*x2,delta)*derCft2z(y1*y2,delta)))/four
      return
      end function
      function d2hgvb01dy22(Eab,Ecd,Dab,Dac,Dad,Dbc,Dbd,Dcd
     $                ,delta,x1,x2,y1,y2)
      USE constants, only: one, four
      USE prec, ONLY: dp
      IMPLICIT NONE
      REAL(KIND=dp) :: Eab, Ecd, Dab, Dac, Dad, Dbc, Dbd, Dcd, delta,
     $                 x1, x2, y1, y2
      REAL(KIND=dp) :: dercft0z, dercft1z, dercft2z
      REAL(KIND=dp) :: d2hgvb01dy22
      d2hgvb01dy22 =
     $  (Eab*y1**2*derCft0z(x1*x2,delta)*derCft2z(y1*y2,delta))/four +
     $ (Ecd*(one-y1)**2*derCft0z((one-x1)*(one-x2),delta)
     $*derCft2z((one-y1)*(one-y2),delta))/four +
     $  (Dbd*(-(y1**2*derCft0z(x1*x2,delta)*derCft2z(y1*(one-y2),delta))
     $ - (one-y1)**2*derCft0z((one-x1)*(one-x2),delta)
     $*derCft2z((one-y1)*y2,delta)))/four +
     $ (Dab*(y1**2*derCft0z((one-x1)*x2,delta)
     $*derCft2z(y1*(one-y2),delta) +
     $(one-y1)**2*derCft0z(x1*(one-x2),delta)
     $*derCft2z((one-y1)*y2,delta)))/four +
     $ (Dcd*(y1**2*derCft0z(x1*(one-x2),delta)
     $*derCft2z(y1*(one-y2),delta) +
     $(one-y1)**2*derCft0z((one-x1)*x2,delta)
     $*derCft2z((one-y1)*y2,delta)))/four +
     $ (Dac*(-(y1**2*derCft0z((one-x1)*(one-x2),delta)
     $*derCft2z(y1*(one-y2),delta)) -
     $(one-y1)**2*derCft0z(x1*x2,delta)
     $*derCft2z((one-y1)*y2,delta)))/four +
     $ (Dbc*(-((one-y1)**2*derCft0z((one-x1)*x2,delta)
     $*derCft2z((one-y1)*(one-y2),delta)) -
     $y1**2*derCft0z(x1*(one-x2),delta)*derCft2z(y1*y2,delta)))/four +
     $ (Dad*(-((one-y1)**2*derCft0z(x1*(one-x2),delta)
     $*derCft2z((one-y1)*(one-y2),delta)) -
     $y1**2*derCft0z((one-x1)*x2,delta)*derCft2z(y1*y2,delta)))/four
      return
      end function d2hgvb01dy22
C
C*MODULE REKS    *DECK REXSOLVER
C>
C> @author     Seunghoon Lee
C>
C> @date       2021, Oct
C>
C> @details    Assemble Microstate Fock matrices into the REKS fock
C
      SUBROUTINE REXSOLVER(EL,DELTA,THRESHOLD,CNVLIM,X,ENREKS,MAXITER,  &
     &                     IOUT,CONST,G1,DBG)
      USE constants, only: zero, one, two
      USE prec, ONLY: dp
      IMPLICIT NONE
C
      LOGICAL :: DSKWRK, GOPARR, MASWRK
      INTEGER :: IBTYP, IPTIM, MASTER, ME, NPROC
      COMMON /PAR   / ME, MASTER, NPROC, IBTYP, IPTIM, GOPARR, DSKWRK,  &
     &                MASWRK
C
      REAL(KIND=dp) :: CNVLIM, CONST, DELTA, ENREKS, G1, THRESHOLD, X
      LOGICAL :: DBG
      INTEGER :: IOUT, MAXITER
      REAL(KIND=dp), DIMENSION(5) :: EL
C
      REAL(KIND=dp) :: ARG, CONSTDEN, DAB, DAC, DAD, DARG, DBC, DBD,    &
     &                 DCD, DER2X, EAB, EAC, EAD, EBC, EBD, ECD, X1
      REAL(KIND=dp) :: dercft1z, dercft2z
C
C uses rexslv4x4 to find the reks(2,2) orbital occupations
C only the (b,c) pair is optimized; the (a,d) pair is kept fixed
C the electronic configurations involved are (aa)(bb) and (aa)(cc)
C the orbital a is always doubly occupied
C the orbital d is empty
      Dab=zero        !set these coupling parameters to zero
      Dac=zero
      Dad=zero
      Dbd=zero
      Dcd=zero
      Ead=zero        !these energies are set to zero
      Ebc=zero
      Ebd=zero
      Ecd=zero
      Eab=El(1)       !energy of the (aa)(bb) configuration
      Eac=El(2)       !energy of the (aa)(cc) configuration
      Dbc=El(3)-El(4) !coupling parameter between the (bb) and (cc) configurations
c     call the reks(4,4) occupations solver with the parameters adjusted for reks(2,2)
      call RexSlv4x4(Eab,Eac,Ebd,Ecd,Dad,Dbc,Delta,CnvLim,
     *               x1,x,EnREKS,IOut,MaxIter,dbg)
C
      ConstDen = -Dbc
      if(abs(ConstDen).gt.1.d-6) then     !do Const and G1 needed for further use in CP-REKS
         Const = (El(2) - El(1))/ConstDen
         if(Const.gt.zero) Const = -Const   !Const is always negative
         arg = x*(one - x)
         darg = one - two*x
         der2x = -two*derCft1z(arg,delta) +
     $            derCft2z(arg,delta)*darg*darg
         G1 = One/(ConstDen*der2x)
      else                                 !special case when the denominator too small
         if(El(2).lt.El(1)) x = One - x
         Const = Zero
         G1 = Zero
      endif
      if(dbg) write(IOut,'(a,F12.8,1x,a,F12.8)')"Const:",Const,"G1:",G1
      return
      end subroutine rexsolver
C*MODULE REKS    *DECK REXCM
C>
C> @author     Seunghoon Lee
C>
C> @date       2021, Oct
C>
      SUBROUTINE REXCM(CM,DNR,DNS,DELTA,WPPS,WOSS,NMIC)
C
      USE constants, only: half, one, four
      USE prec, ONLY: dp
      IMPLICIT NONE
C
      REAL(KIND=dp) :: DELTA, DNR, DNS, WOSS, WPPS
      INTEGER :: NMIC
      REAL(KIND=dp), DIMENSION(NMIC) :: CM
C
      REAL(KIND=dp) :: F, TMP
      REAL(KIND=dp) :: REXCONVF
      LOGICAL :: STATAVG
C
      STATAVG = abs(WPPS - ONE).GT.1.D-8

      TMP = FOUR*DNR*DNS
      F   = WPPS*REXCONVF(TMP,DELTA)
C
      IF(STATAVG) THEN
         CM(1)= WPPS*DNR
         CM(2)= WPPS*DNS
         CM(3)= WOSS-HALF*F
         CM(4)= HALF*F - HALF*WOSS
      ELSE
         CM(1)= DNR
         CM(2)= DNS
         CM(3)= -HALF*F
         CM(4)= -CM(3)
      ENDIF
      RETURN
      END SUBROUTINE REXCM
C*MODULE REKS    *DECK REXCONVF
      FUNCTION REXCONVF(DNRNS,DELTA)
      USE constants, only:half, one
      USE prec, ONLY: dp
      IMPLICIT NONE
C
      REAL(KIND=dp) :: DELTA, DNRNS
      REAL(KIND=dp) :: REXCONVF
C
      REAL(KIND=dp) :: TMP
C © Michael Filatov 2021
C
      TMP      = HALF*(DNRNS+DELTA)/(ONE+DELTA)
      REXCONVF = DNRNS**(ONE-TMP)
      RETURN
      END FUNCTION REXCONVF
C*MODULE REKS    *DECK REXFM2FE
C>
C> @author     Seunghoon Lee
C>
C> @date       2021, Oct
C>
      SUBROUTINE REXFM2FE(L,CL,FR,FS,L1,L2,L3,NCORE,NVIR,FA,FB,FREX,    &
     &                    WREX,WRS,IA)
C
      USE comm_REKSCM, ONLY: NMICRO
      USE constants, only: zero, half, one
      USE prec, ONLY: dp
      IMPLICIT NONE
C
      REAL(KIND=dp), PARAMETER :: LOWLIM = 1.0D-08
C
C
      REAL(KIND=dp) :: CL, FR, FS, WRS
      INTEGER :: L, L1, L2, L3, NCORE, NVIR
      REAL(KIND=dp), DIMENSION(L2) :: FA, FB, FREX, WREX
      INTEGER, DIMENSION(*) :: IA
C
      REAL(KIND=dp) :: NRA, NRB, NSA, NSB, SIGNRS, WC, WRA, WRB, WRCA,  &
     &                 WRCB, WSA, WSB, WSCA, WSCB, WSRA, WSRB,          &
     &                 WVRA, WVRB, WVSA, WVSB
      INTEGER :: I, IJ, J
C
C     ASSIGN OCCUPATION
c     IRA=NMICRO(1,1,L)
c     IF(IRA.EQ.1) NRA=ZERO
c     IF(IRA.EQ.2) NRA=ONE
c     ISA=NMICRO(2,1,L)
c     IF(ISA.EQ.1) NSA=ZERO
c     IF(ISA.EQ.2) NSA=ONE
c     IRB=NMICRO(1,2,L)
c     IF(IRB.EQ.1) NRB=ZERO
c     IF(IRB.EQ.2) NRB=ONE
c     ISB=NMICRO(2,2,L)
c     IF(ISB.EQ.1) NSB=ZERO
c     IF(ISB.EQ.2) NSB=ONE
C
      NRA = dble(NMICRO(1,1,L) - 1)
      NSA = dble(NMICRO(2,1,L) - 1)
      NRB = dble(NMICRO(1,2,L) - 1)
      NsB = dble(NMICRO(2,2,L) - 1)
C     ASSIGN WEIGHTING FACTORS
      WC  =HALF*CL
      WRA =WC*NRA
      IF(FR.GT.LOWLIM) WRA =WRA/FR
      WRB =WC*NRB
      IF(FR.GT.LOWLIM) WRB =WRB/FR
      WSA =WC*NSA
      IF(FS.GT.LOWLIM) WSA =WSA/FS
      WSB =WC*NSB
      IF(FS.GT.LOWLIM) WSB =WSB/FS
      WRCA=WC*(ONE-NRA)
      IF(ONE-FR.GT.LOWLIM) WRCA=WRCA/(ONE-FR)
      WRCB=WC*(ONE-NRB)
      IF(ONE-FR.GT.LOWLIM) WRCB=WRCB/(ONE-FR)
      WSCA=WC*(ONE-NSA)
      IF(ONE-FS.GT.LOWLIM) WSCA=WSCA/(ONE-FS)
      WSCB=WC*(ONE-NSB)
      IF(ONE-FS.GT.LOWLIM) WSCB=WSCB/(ONE-FS)
      WVRA=WRA
      WVRB=WRB
      WVSA=WSA
      WVSB=WSB
      SIGNRS=ONE
      IF(FR-FS.LT.ZERO) SIGNRS=-SIGNRS
      WSRA=CL*(NRA-NSA)*SIGNRS
      WSRB=CL*(NRB-NSB)*SIGNRS
C
C     CALCULATE REKS FOCK FOR THIS MICROSTATE
C
C     FC
      DO I=1,NCORE
      DO J=1,I
         IJ=IA(I)+J
         FREX(IJ)=FREX(IJ)+WC*(FA(IJ)+FB(IJ))
         WREX(IJ)=WREX(IJ)+CL*(FA(IJ)+FB(IJ))
      ENDDO
      ENDDO
C     FV
      DO I=NVIR,L1
      DO J=NVIR,I
         IJ=IA(I)+J
         FREX(IJ)=FREX(IJ)+WC*(FA(IJ)+FB(IJ))
      ENDDO
      ENDDO
C     FVC
      DO I=NVIR,L1
      DO J=1,NCORE
         IJ=IA(I)+J
         FREX(IJ)=FREX(IJ)+WC*(FA(IJ)+FB(IJ))
      ENDDO
      ENDDO
C
C     FR
      I=NCORE+1
      J=NCORE+1
      IJ=IA(I)+J
      FREX(IJ)=FREX(IJ)+WRA*FA(IJ)+WRB*FB(IJ)
      WREX(IJ)=WREX(IJ)+CL*(NRA*FA(IJ)+NRB*FB(IJ))
C     FS
      I=NCORE+2
      J=NCORE+2
      IJ=IA(I)+J
      FREX(IJ)=FREX(IJ)+WSA*FA(IJ)+WSB*FB(IJ)
      WREX(IJ)=WREX(IJ)+CL*(NSA*FA(IJ)+NSB*FB(IJ))
C
C     FRC
      I=NCORE+1
      DO J=1,NCORE
         IJ=IA(I)+J
         FREX(IJ)=FREX(IJ)+WRCA*FA(IJ)+WRCB*FB(IJ)
         WREX(IJ)=WREX(IJ)+CL*(NRA*FA(IJ)+NRB*FB(IJ))
      ENDDO
C     FSC
      I=NCORE+2
      DO J=1,NCORE
         IJ=IA(I)+J
         FREX(IJ)=FREX(IJ)+WSCA*FA(IJ)+WSCB*FB(IJ)
         WREX(IJ)=WREX(IJ)+CL*(NSA*FA(IJ)+NSB*FB(IJ))
      ENDDO
C
C     FVR
      J=NCORE+1
      DO I=NVIR,L1
         IJ=IA(I)+J
         FREX(IJ)=FREX(IJ)+WVRA*FA(IJ)+WVRB*FB(IJ)
      ENDDO
C     FVS
      J=NCORE+2
      DO I=NVIR,L1
         IJ=IA(I)+J
         FREX(IJ)=FREX(IJ)+WVSA*FA(IJ)+WVSB*FB(IJ)
      ENDDO
C
C     FSR
      I=NCORE+2
      J=NCORE+1
      IJ=IA(I)+J
      FREX(IJ)=FREX(IJ)+WSRA*FA(IJ)+WSRB*FB(IJ)
      WRS = WRS + CL*(NRA*FA(IJ) + NRB*FB(IJ))
      WREX(IJ)=WRS
c     WREX(IJ)=WREX(IJ)+CL*(NSA*FA(IJ)+NSB*FB(IJ))
C
      RETURN
      END SUBROUTINE REXFM2FE
C*MODULE REKS    *DECK REXEM2EE
      FUNCTION REXEM2EE(CM,EM,NMIC)
      USE constants, only: zero, one, two
      USE prec, ONLY: dp
      IMPLICIT NONE
C
      INTEGER :: NMIC
      REAL(KIND=dp), DIMENSION(NMIC) :: CM, EM
      REAL(KIND=dp) :: REXEM2EE
C
      REAL(KIND=dp) :: FACT
      INTEGER :: I
C © Michael Filatov 2021
C
      REXEM2EE = ZERO
      FACT=ONE
      DO I=1,NMIC
         IF(I.GE.3) FACT=TWO
         REXEM2EE = REXEM2EE + FACT*CM(I)*EM(I)
      ENDDO
      RETURN
      END FUNCTION REXEM2EE

C*MODULE REKS    *DECK REXOFFFOC
      FUNCTION REXOFFFOC(F,L1,NCORE)
      USE prec, ONLY: dp
      IMPLICIT NONE
C
      INTEGER :: L1, NCORE
      REAL(KIND=dp), DIMENSION(*) :: F
      REAL(KIND=dp) :: REXOFFFOC
C
      INTEGER :: I, IJ, IMAX, J, JMAX
C © Michael Filatov 2021
C
      REXOFFFOC = 0.d0
C
c           CALL PRTRIL(F,L1)
C
      ij = ncore*(ncore+1)/2                                             !skip all core orbitals
      do j = 1,ncore                                                     !loop over core-r block
        ij = ij + 1
        if(abs(f(ij)).gt.rexofffoc) then
           rexofffoc = abs(f(ij))
           imax = ncore+1
           jmax = j
        endif
      enddo
C
      ij = ij +1                                                         !increment ij for r-r element
      do j = 1,ncore+1                                                   !loop over core-s and r-s blocks
        ij = ij + 1
        if(abs(f(ij)).gt.rexofffoc) then
           rexofffoc = abs(f(ij))
           imax = ncore+2
           jmax = j
        endif
      enddo
C
      ij = ij + 1                                                        !increment ij for s-s block
      do i = 1,l1 - ncore - 2                                            !loop over core-virt, r-virt, and s-virt blocks; note, i \in [1,L1-ncore-2], but i \in [ncore+3,L1]
       do j = 1,ncore+2
         ij = ij + 1
         if(abs(f(ij)).gt.rexofffoc) then
            rexofffoc = abs(f(ij))
            imax = i + ncore + 2
            jmax = j
         endif
       enddo
       ij = ij + i + ncore + 3 - j                                       !increment ij for missing cells, i+1-j
      enddo
C
c     write(6,*)' imax, jmax =',imax,jmax
      RETURN
      END FUNCTION REXOFFFOC
C*MODULE REKS    *DECK REXDIISER
C>
C> @author     Michael Filatov
C>
C> @date       2021, Oct
C>
      SUBROUTINE REXDIISER(ERR,DNR,DNS,L1,NCORE,errval)
      USE prec, ONLY: dp
      IMPLICIT NONE
C
      REAL(KIND=dp) :: DNR, DNS, ERRVAL
      INTEGER :: L1, NCORE
      REAL(KIND=dp), DIMENSION(L1,L1) :: ERR
C
      REAL(KIND=dp) :: DCMR, DCMS, DRMS, DRMV, DSMV
      INTEGER :: I, J

      DCMR = 1.d0 - DNR
      DCMS = 1.d0 - DNS
      DRMV = DNR
      DSMV = DNS
      DRMS = DNR - DNS

c     DCMR = 1.d0
c     DCMS = 1.d0
c     DRMV = 1.d0
c     DSMV = 1.d0
c     DRMS = 1.d0
C v-v
      DO I=NCORE+3,L1
      DO J=NCORE+3,L1
         ERR(I,J) = 0.d0
      ENDDO
      ENDDO
C c-c
      DO I=1,NCORE
      DO J=1,NCORE
         ERR(I,J) = 0.d0
      ENDDO
      ENDDO
C c-r
      do i=1,ncore
        Err(ncore+1,i) = err(ncore+1,i)*dcmr
        err(i,ncore+1) = err(i,ncore+1)*dcmr
      enddo
C c_s
      do i=1,ncore
        Err(ncore+2,i) = err(ncore+2,i)*dcms
        err(i,ncore+2) = err(i,ncore+2)*dcms
      enddo
C r-v
      do i=ncore+3,l1
        err(ncore+1,i) = err(ncore+1,i)*drmv
        err(i,ncore+1) = err(i,ncore+1)*drmv
      enddo
C s-v
      do i=ncore+3,l1
        err(ncore+2,i) = err(ncore+2,i)*dsmv
        err(i,ncore+2) =err(i,ncore+2)*dsmv
      enddo
C r-s
      err(ncore+1,ncore+2) = err(ncore+1,ncore+2)*drms
      err(ncore+2,ncore+1) = err(ncore+2,ncore+1)*drms
C r-r, s-s
      err(ncore+1,ncore+1) = 0.d0
      err(ncore+2,ncore+2) = 0.d0
C find max value of the error
      errval =  0.d0
      do i = 1,L1
       do j = 1,L1
        if(abs(err(i,j)).gt.errval) errval = abs(err(i,j))
       enddo
      enddo
C
      RETURN
      END SUBROUTINE REXDIISER
C*MODULE REKS    *DECK LOCALZ
C>
C> @author     Michael Filatov
C>
C> @date       2021, Oct
C>
C> @details    Localize a pair of fractionally occupied orbitals
C>             minimizes ff = \sum_{i,j} [v1(i)*S(i,j)*v2(j)]^2
C>
      SUBROUTINE localz(ovlp,v1,v2,num,iw)
      USE prec, ONLY: dp
      IMPLICIT NONE
C
      LOGICAL :: DSKWRK, GOPARR, MASWRK
      INTEGER :: IBTYP, IPTIM, MASTER, ME, NPROC
      COMMON /PAR   / ME, MASTER, NPROC, IBTYP, IPTIM, GOPARR, DSKWRK,  &
     &                MASWRK
C
      INTEGER :: IW, NUM
      REAL(KIND=dp), DIMENSION(*) :: OVLP, V1, V2                        !S(i,j) is in lower triangular form
C
      REAL(KIND=dp) :: A, B, CS, DEN, F1, F2, FF, OVLP2, PHI, PHINEW,    &
     &                PIOVER4, SS, TEMP, TEMP0, TEMP1, TEMP2, TEST0,    &
     &                TEST1, TMP21, TMP22, ULIM, VD, VMAX
      INTEGER :: I, IJ, IMAX, ITR, J, JJ, MAXIT
      REAL(KIND=dp), DIMENSION(num) :: V1OLD, V1V1, V1V2, V2OLD, V2V2

      ulim = 1.d-8                                                       !convergence criterion
      piover4 = asin(1.d0)/2.d0
      phi = 0.d0
      maxit = 20

      call dcopy(num,v1,1,v1old,1)                                       !save the old vectors
      call dcopy(num,v2,1,v2old,1)

      test0 = 0.d0
      do j = 1,num                                                       !compute -ff- btw the old vectors
       temp = 0.d0
       jj = (j*j-j)/2
       do i = 1,num
        ij = (i*i-i)/2 + j                                               !lower triangular index
        if(i.lt.j)ij = jj + i                                            !use S(i,j) symmetry if in the upper triangle
        temp = temp + (v1old(i) * ovlp(ij))**2
       enddo
       test0 = test0 + temp*v2old(j)**2
      enddo

      vmax = 0.d0
      imax = 0
      do i = 1,num                                                       !find the largest element of -v1- and its index -imax-
       vd = abs(v1(i))
        if(vd.gt.vmax)then
         vmax = vd
         imax = i
        endif
      enddo

      test1 = 0.d0
      if(test0.gt.ulim)then

      a = vmax                                                           !make v2(imax) to zero by orthogonal rotation
      b = v2(imax)
      den = 1.d0/sqrt(a*a + b*b)
      cs = a*den
      ss = b*den*sign(1.d0,v1(imax))
      phinew = asin(ss)                                                  !this is the angle

      itr = 0
  100 phi = phinew                                                       !top of the loop; Newton-Raphson iterations

      do i  = 1,num                                                      !apply the orthogonal rotation
       v1(i) = v1old(i)*cs + v2old(i)*ss
       v2(i) =-v1old(i)*ss + v2old(i)*cs
       v1v2(i) = v1(i)*v2(i)                                             !do intermediate variables
       v1v1(i) = v1(i)**2
       v2v2(i) = v2(i)**2
      enddo

      ff = 0.d0                                                          !-ff- is the target function
      f1 = 0.d0                                                          !-f1- is its first derivative, w.r.t. rotation angle -phi-
      f2 = 0.d0                                                          !-f2- is its second derivative
      do j = 1,num
       temp0 = 0.d0
       temp1 = 0.d0
       temp2 = 0.d0
       tmp21 = 0.d0
       tmp22 = 0.d0
       jj = (j*j-j)/2
       do i = 1,num
        ij = (i*i-i)/2 + j
        if(i.lt.j)ij = jj + i
        ovlp2 = ovlp(ij)**2                                              !square the overlap
        temp0 = temp0 + v1v1(i)*ovlp2
        temp1 = temp1 + (v2v2(i) - v1v1(i))*ovlp2
        temp2 = temp2 + v1v2(i)*ovlp2
        tmp21 = temp0
        tmp22 = tmp22 + v2v2(i)*ovlp2
       enddo
        ff = ff + temp0*v2v2(j)                                          !ff = \sum_{i,j} [v1(i)*S(i,j)*v2(j)]^2
        f1 = f1 + temp1*v1v2(j)                                          !f1 = 2*\sum_{i,j} [v2(i)^2 - v1(i)^2]*S(i,j)**2*v1(j)*v2(j)
        f2 = f2 + 2.d0*tmp21*v1v1(j) + 2.d0*tmp22*v2v2(j)                !f2 = 2*\sum{i,j} [v1(i)^2*S(i,j)^2*v1(j)^2 + v2(i)^2*S(i,j)^2*v2(j)^2]
     *          - 8.d0*temp2*v1v2(j)                                     !   - 8*\sum_{i,j} v1(i)*v2(i)*S(i,j)^2*v1(j)*v2(j) - 4*ff
      enddo
      f2 = f2 - 4.d0*ff
      if(f2.lt.0.d0) f2 = -f2                                            !inverse sign of -f2- if it's negative
      test1 = ff

      phinew = phi - f1/f2                                               !Newton-Raphson increment
      if(abs(phinew).gt.piover4)phinew=-sign(piover4,phinew)
      ss = sin(phinew)
      cs = cos(phinew)

      itr = itr + 1

      if(abs(f1).gt.ulim.and.itr.le.maxit)goto 100                       !if not converged, go to the top of the loop

      if(abs(phi).gt.1.d-6.and.test1.lt.test0)then                       !if the absolute overlap is smaller after rotation
      if(maswrk)
     *write(iw,'(1x,"REKS localization: before ",F12.8,", after ",F12.8,
     *", angle ",F10.6," rad.")')test0,test1,phi
      else                                                               !if failed to minimize -ff-
      call dcopy(num,v1old,1,v1,1)
      call dcopy(num,v2old,1,v2,1)
      if(maswrk)
     *write(iw,'(1x,"REKS localization failed: before ",F12.8,", after "
     *,F12.8,", angle ",F10.6," rad. Using old vectors")')test0,test1
     *,phi
      endif

      endif                                                              !if(test0.gt.ulim)then

      return
      end subroutine localz
C*MODULE REKS    *DECK DELCLZ
C>
C> @author     Michael Filatov
C>
C> @date       2021, Oct
C>
C> @details    De-localize a pair of fractionally occupied orbitals
C>             maximizes ff = \sum_{i,j} [v1(i)*S(i,j)*v2(j)]^2
C>
      SUBROUTINE delclz(ovlp,v1,v2,num,iw)
      USE prec, ONLY: dp
      IMPLICIT NONE
C
      LOGICAL :: DSKWRK, GOPARR, MASWRK
      INTEGER :: IBTYP, IPTIM, MASTER, ME, NPROC
      COMMON /PAR   / ME, MASTER, NPROC, IBTYP, IPTIM, GOPARR, DSKWRK,  &
     &                MASWRK
C
      INTEGER :: IW, NUM
      REAL(KIND=dp), DIMENSION(*) :: OVLP, V1, V2                        !S(i,j) is in lower triangular form
C
      REAL(KIND=dp) :: CS, F1, F2, FF, OVLP2, PHI, PHINEW, PIOVER4, SIG, &
     &                 SS, STEP, TEMP0, TEMP1, TEMP2, TEST0, TEST1,      &
     &                 TMP21, TMP22, ULIM
      INTEGER :: I, IJ, ITR, J, JJ, MAXIT
      REAL(KIND=dp), DIMENSION(num) :: V1OLD, V1V1, V1V2, V2OLD, V2V2

      ulim = 1.d-8                                                       !convergence criterion
      piover4 = asin(1.d0)/2.d0
      phi = 0.d0
      maxit = 20

      call dcopy(num,v1,1,v1old,1)                                       !save the old vectors
      call dcopy(num,v2,1,v2old,1)

      test0 = 0.d0
      f1 = 0.d0
      f2 = 0.d0
      do j = 1,num
       temp0 = 0.d0
       temp1 = 0.d0
       temp2 = 0.d0
       tmp21 = 0.d0
       tmp22 = 0.d0
       jj = (j*j-j)/2
       do i = 1,num
        ij = (i*i-i)/2 + j                                               !lower triangular index
        if(i.lt.j)ij = jj + i                                            !use S(i,j) symmetry if in the upper triangle
        ovlp2 = ovlp(ij)**2                                              !square the overlap
        temp0 = temp0 + v1old(i)**2 * ovlp2
        temp1 = temp1 + (v2old(i)*v2old(i) - v1old(i)*v1old(i))*ovlp2
        temp2 = temp2 + v1old(i)*v2old(i)*ovlp2
        tmp21 = temp0
        tmp22 = tmp22 + v2old(i)*v2old(i)*ovlp2
       enddo
       test0 = test0 + temp0*v2old(j)**2
       f1 = f1 + temp1*v1old(j)*v2old(j)
       f2 = f2 + 2.d0*tmp21*v1old(j)*v1old(j)
     *         + 2.d0*tmp22*v2old(j)*v2old(j)
     *         - 8.d0*temp2*v1old(j)*v2old(j)
      enddo
      f2 = f2 - 4.d0*test0

      if(abs(f1).le.ulim) return                                         !exit if the first derivative is already small

      phi = 0.d0
      step = piover4
      sig = 0.5d0
      if(f2.lt.-0.1d0) goto 200                                          !go to Newton-Raphson iterations; skip the next step

      itr = 0
  100 phi = phi + sig*step                                               !top of the loop; locate negative second derivative
      ss = sin(phi)
      cs = cos(phi)

      do i  = 1,num                                                      !apply the orthogonal rotation
       v1(i) = v1old(i)*cs + v2old(i)*ss
       v2(i) =-v1old(i)*ss + v2old(i)*cs
       v1v2(i) = v1(i)*v2(i)                                             !do intermediate variables
       v1v1(i) = v1(i)**2
       v2v2(i) = v2(i)**2
      enddo

      ff = 0.d0                                                          !-ff- is the target function
      f1 = 0.d0                                                          !-f1- is its first derivative, w.r.t. rotation angle -phi-
      f2 = 0.d0                                                          !-f2- is its second derivative
      do j = 1,num
       temp0 = 0.d0
       temp1 = 0.d0
       temp2 = 0.d0
       tmp21 = 0.d0
       tmp22 = 0.d0
       jj = (j*j-j)/2
       do i = 1,num
        ij = (i*i-i)/2 + j
        if(i.lt.j)ij = jj + i
        ovlp2 = ovlp(ij)**2                                              !square the overlap
        temp0 = temp0 + v1v1(i)*ovlp2
        temp1 = temp1 + (v2v2(i) - v1v1(i))*ovlp2
        temp2 = temp2 + v1v2(i)*ovlp2
        tmp21 = temp0
        tmp22 = tmp22 + v2v2(i)*ovlp2
       enddo
        ff = ff + temp0*v2v2(j)                                          !ff = \sum_{i,j} [v1(i)*S(i,j)*v2(j)]^2
        f1 = f1 + temp1*v1v2(j)                                          !f1 = 2*\sum_{i,j} [v2(i)^2 - v1(i)^2]*S(i,j)**2*v1(j)*v2(j)
        f2 = f2 + 2.d0*tmp21*v1v1(j) + 2.d0*tmp22*v2v2(j)                !f2 = 2*\sum{i,j} [v1(i)^2*S(i,j)^2*v1(j)^2 + v2(i)^2*S(i,j)^2*v2(j)^2]
     *          - 8.d0*temp2*v1v2(j)                                     !   - 8*\sum_{i,j} v1(i)*v2(i)*S(i,j)^2*v1(j)*v2(j) - 4*ff
      enddo
      f2 = f2 - 4.d0*ff
      test1 = ff
      if(f1.lt.0.d0.and.sig.gt.0.d0) sig = -sig                          !inverse the increment, if first derivative becomes negative

      step = step*0.5d0                                                  !divide step by 2

      itr = itr + 1

      if(f2.ge.-0.1d0.and.itr.le.maxit)goto 100                          !if not converged, go to the top of the loop

      if(f2.ge.0.d0.and.itr.gt.maxit)then                                !failed to locate negative -f2-
      if(maswrk)
     *write(iw,'(1x,"Unable to find maximum after",i4," iterations")')
     *      itr
      call dcopy(num,v1old,1,v1,1)                                       !revert to the old vectors
      call dcopy(num,v2old,1,v2,1)
      return                                                             !and exit
      endif

  200 itr = 0
      phinew = phi

  300 phi = phinew                                                       !top of the loop; Newton-Raphson iterations
      ss = sin(phi)
      cs = cos(phi)

      do i  = 1,num                                                      !apply the orthogonal rotation
       v1(i) = v1old(i)*cs + v2old(i)*ss
       v2(i) =-v1old(i)*ss + v2old(i)*cs
       v1v2(i) = v1(i)*v2(i)                                             !do intermediate variables
       v1v1(i) = v1(i)**2
       v2v2(i) = v2(i)**2
      enddo

      ff = 0.d0                                                          !-ff- is the target function
      f1 = 0.d0                                                          !-f1- is its first derivative, w.r.t. rotation angle -phi-
      f2 = 0.d0                                                          !-f2- is its second derivative
      do j = 1,num
       temp0 = 0.d0
       temp1 = 0.d0
       temp2 = 0.d0
       tmp21 = 0.d0
       tmp22 = 0.d0
       jj = (j*j-j)/2
       do i = 1,num
        ij = (i*i-i)/2 + j
        if(i.lt.j)ij = jj + i
        ovlp2 = ovlp(ij)**2                                              !square the overlap
        temp0 = temp0 + v1v1(i)*ovlp2
        temp1 = temp1 + (v2v2(i) - v1v1(i))*ovlp2
        temp2 = temp2 + v1v2(i)*ovlp2
        tmp21 = temp0
        tmp22 = tmp22 + v2v2(i)*ovlp2
       enddo
        ff = ff + temp0*v2v2(j)                                          !ff = \sum_{i,j} [v1(i)*S(i,j)*v2(j)]^2
        f1 = f1 + temp1*v1v2(j)                                          !f1 = 2*\sum_{i,j} [v2(i)^2 - v1(i)^2]*S(i,j)**2*v1(j)*v2(j)
        f2 = f2 + 2.d0*tmp21*v1v1(j) + 2.d0*tmp22*v2v2(j)                !f2 = 2*\sum{i,j} [v1(i)^2*S(i,j)^2*v1(j)^2 + v2(i)^2*S(i,j)^2*v2(j)^2]
     *          - 8.d0*temp2*v1v2(j)                                     !   - 8*\sum_{i,j} v1(i)*v2(i)*S(i,j)^2*v1(j)*v2(j) - 4*ff
      enddo
      f2 = f2 - 4.d0*ff
      test1 = ff

      phinew = phi - f1/f2                                               !Newton-Raphson increment

      itr = itr + 1

      if(abs(f1).gt.ulim.and.itr.le.maxit)goto 300                       !if not converged, go to the top of the loop
      if(maswrk)
     *write(iw,'(1x,"REKS de-localization: before ",F12.8,", after ",
     *F12.8,", angle ",F10.6," rad.")')test0,test1,phi

      return
      end subroutine delclz
C*MODULE REKS    *DECK REXINPUT
C>
C> @author     Michael Filatov
C>
C> @date       2021, Oct
C>
      SUBROUTINE rexinput
      USE comm_REKSCM, ONLY: NMICRO, MTTYP, WPPS, WOSS, G1, DNR,
     * DNS, DELTA, FR, FS
      USE comm_REXOPT, ONLY: REXTYPE, REXTARGET, REXSHIFT, REXDIIS,
     * REXLDL, RLXDEN, REXEKT, EKTEA, REXCG, RXCGIT, RXCGTH
      USE prec, ONLY: dp
      IMPLICIT NONE
C
      LOGICAL :: DSKWRK, GOPARR, MASWRK
      INTEGER :: IBTYP, IDAF, IJK, IJKT, IP, IPTIM, IR, IW, MASTER, ME  &
     &           , NAV, NPROC
      INTEGER, DIMENSION(950) :: IODA
      COMMON /IOFILE/ IR, IW, IP, IJK, IJKT, IDAF, NAV, IODA             !common block with I/O files
      COMMON /PAR   / ME, MASTER, NPROC, IBTYP, IPTIM, GOPARR, DSKWRK,  &!common block for parallel execution parameters
     &                MASWRK
C
      CHARACTER(2) :: BLNK
      CHARACTER(16) :: DUMMY
      LOGICAL :: FOUND, INVEC, IROPND
      INTEGER :: IINDEX, LENDUM, LENSTR
      CHARACTER(256) :: LINE, STRING                                     !declare strings for reading
      blnk = " |"                                                        !terminating character for string; "|" does no harm
C
C Default values for REKS/SA-REKS/SSR calculation
C rexType   = 0.....SA-REKS(2,2); SA: PPS+OSS
C             1.....2SI-2SA-REKS(2,2); SA: PPS+OSS; SI: PPS+OSS
C             2.....3SI-2SA-REKS(2,2); SA: PPS+OSS; SI: PPS+OSS+DES
C
C rexTarget = 1.....either PPS (of SA) or S0 (of SSR) is reported
C             2.....either OSS of S1 is reported
C
C WPPS      = float.....weighting factor of PPS state in SA-REKS(2,2) ensemble (0.<=WPPS<-0.5)
C
C rexShift  = float.....orbital level shift used in SA-REKS calculation
C
C rexDIIS   = yes/no....to do or not to do DIIS with REKS; default: to do
C
C rexLdL    = 0.....do nothing
C             1.....localize the REKS active orbitals by minimizing the square of their absolute overlap
C             2.....delocalize the REKS active orbitals
C
C rlxDen    = yes/no....to do/or not to do the relaxed density matrix for the target state
C
C rexEKT    = yes/no....to compute/not to compute the Extended Koopmans' Theorem energies and orbitals
C
C rexCG     = yes/no....to use conjugate gradient for solving CP-REKS equations; default yes; if no, then use direct solver
C
C rxCGit    = N....max number of CG iterations; default 100
C
C rxCGth    = float....convergence criterion for CG
C
      rexType = 0                                                        !set default to SA-REKS(2,2)
      rexTarget = 0                                                      !set default to averaged state
      wpps = 0.5d0                                                       !set default to equiensemble of PPS and OSS states
      rexShift = 0.3d0                                                   !set default to 0.3 a.u.
      rexDIIS = .true.
      rlxDen = .false.
      rexEKT = .false.
      EKTEA = .false.
      rexLdL = 0                                                         !by default: do nothing
      rexCG = .true.                                                     !use CG solver by default
      rxCGit = 100                                                       !max CG iterations
      rxCGth = 1.d-6                                                     !CG convergence criterion
C
      if(maswrk)then                                                     !do on master proc and broadcast to slaves
C
      inquire(unit=ir,opened=iropnd)                                     !inquire if the inut file is opened
      if(.not.iropnd)open(unit=ir,status='old')                          !open the input file for reading, if it's not opened already
      rewind(unit=ir)                                                    !rewind it
      invec=.false.                                                      !inside $VEC block; this block is skipped
      found=.false.                                                      !if the $REKS block is found
      do                                                                 !endless loop
  101 read(ir,'(a)',end=999)string                                       !read in a string from the input file
      lenstr = len(trim(string))                                         !get length of the string
      call uprcas(string,lenstr)                                         !to UPPER case; internal gamess feature
      if(index(string,"$VEC").ne.0)then
      do
      read(ir,'(a)',end=999)dummy                                        !shorter reads to "dummy" (c16)
      lendum = len(trim(dummy))                                          !get length of "dummy"
      call uprcas(dummy,lendum)                                          !convert to UPPER case
      if(index(dummy,"$END").ne.0)goto 101                               !break out the loop
      end do                                                             !make an endless loop
      endif
      if(index(string,"$REKS").ne.0)then                                 !search for opening $REKS keyword
       found=.true.                                                      !make it .true., if found
       if(index(string,"$END").ne.0)then                                 !if $END is in the same string
       goto 999                                                          !break out of the loop if $REKS and $END in the same string
       else                                                              !if $END is not in the same string
C
  111  read(ir,'(a)',end=999)line                                        !make another endless loop; read the next string in
       lenstr = len(trim(line))
       call uprcas(line,lenstr)
       string = trim(string)//blnk                                       !paranoia; append blank at the end of the string
       string = trim(string)//trim(line)                                 !concatenate with the previous string
       if(index(string,"$END").ne.0)then                                 !if $END found, bail out
        goto 999
       else                                                              !if not
        goto 111                                                         !continue looping
       endif                                                             !if(index(string,"$END").ne.0)
C
       endif                                                             !if(index(string,"$END").ne.0)
      endif                                                              !if(index(string,"$REKS").ne.0)
      enddo                                                              !end of the endless loop
  999 continue                                                           !break from the endless loop
      if(.not.found)then                                                 !what to do, if $REKS list not found
       write(iw,'(3(/1x,a))')
     * '====> REKS calcuation requested, but $REKS list not found!',
     * '====> Running default SA-REKS(2,2) calculation!',
     * '====> Not happy, call tall free 1-800-pissoff'
      goto 1000
      endif                                                              !if(.not.found)then
      write(iw,'(/1x,a,/1x,a,/1x)')
     *      '=====> REKS control string found <=====',
     *      trim(string)
      lenstr = len(trim(string))
      call rmcommas(string,lenstr)                                       !replace "," with " "
C
      if(index(string,"REXTYPE").ne.0)then                               !found REXTYPE keyword
       iindex = index(string,"REXTYPE")
       iindex = iindex + 8
       call skpblnks(string,iindex)                                      !move -iindex- to first non-blank character ("=") after the keyword
       read(string(iindex:),*)rexType
       write(iw,'(1x,a,i4)')'======> set rexType to',rexType
      endif
      if(index(string,"REXTARGET").ne.0)then                             !found REXTARGET keyword
       iindex = index(string,"REXTARGET")
       iindex = iindex + 10
       call skpblnks(string,iindex)
       read(string(iindex:),*)rexTarget
       write(iw,'(1x,a,i4)')'======> set rexTarget to',rexTarget
      endif
      if(index(string,"WPPS").ne.0)then                                  !found REXTARGET keyword
       iindex = index(string,"WPPS")
       iindex = iindex + 5
       call skpblnks(string,iindex)
       read(string(iindex:),*)wpps
       write(iw,'(1x,a,f14.8)')'======> set WPPS to',wpps
      endif
      if(index(string,"REXSHIFT").ne.0)then                              !found REXSHIFT keyword
       iindex = index(string,"REXSHIFT")
       iindex = iindex + 9
       call skpblnks(string,iindex)
       read(string(iindex:),*)rexShift
       write(iw,'(1x,a,f14.8)')'======> set REXSHIFT to',rexShift
      endif
      if(index(string,"REXDIIS").ne.0)then                               !found REXDIIS keyword
       iindex = index(string,"REXDIIS")
       iindex = iindex + 8
       call skpblnks(string,iindex)
       read(string(iindex:),*)line
       lenstr = len(trim(line))
       call uprcas(line,lenstr)
       if(index(line,"NO").ne.0)rexDIIS=.false.
       write(iw,'(1x,a,l8)')'======> set REXDIIS to',rexDIIS
      endif
      if(index(string,"RLXDEN").ne.0)then                                !found RLXDEN keyword
       iindex = index(string,"RLXDEN")
       iindex = iindex + 7
       call skpblnks(string,iindex)
       read(string(iindex:),*)line
       lenstr = len(trim(line))
       call uprcas(line,lenstr)
       if(index(line,"YES").ne.0)rlxDen=.true.
       write(iw,'(1x,a,l8)')'======> set RLXDEN to',rlxDen
      endif
      if(index(string,"REXEKT").ne.0)then                                !found REXEKT keyword
       iindex = index(string,"REXEKT")
       iindex = iindex + 7
       call skpblnks(string,iindex)
       read(string(iindex:),*)line
       lenstr = len(trim(line))
       call uprcas(line,lenstr)
       if(index(line,"YES").ne.0)rexEKT=.true.
       write(iw,'(1x,a,l8)')'======> set REXEKT to',rexEKT
      endif
      if(index(string,"EKTEA").ne.0)then                                 !found EKTEA keyword
       iindex = index(string,"EKTEA")
       iindex = iindex + 6
       call skpblnks(string,iindex)
       read(string(iindex:),*)line
       lenstr = len(trim(line))
       call uprcas(line,lenstr)
       if(index(line,"YES").ne.0)EKTEA=.true.
       write(iw,'(1x,a,l8)')'======> set EKTEA to',EKTEA
      endif
      if(index(string,"REXCG").ne.0)then                                 !found REXCG keyword
       iindex = index(string,"REXCG")
       iindex = iindex + 6
       call skpblnks(string,iindex)
       read(string(iindex:),*)line
       lenstr = len(trim(line))
       call uprcas(line,lenstr)
       if(index(line,"NO").ne.0)rexCG=.false.
       write(iw,'(1x,a,l8)')'======> set REXCG to',rexCG
      endif
      if(index(string,"RXCGIT").ne.0)then                                !found RXCGIT keyword
       iindex = index(string,"RXCGIT")
       iindex = iindex + 7
       call skpblnks(string,iindex)
       read(string(iindex:),*)rxCGit
       write(iw,'(1x,a,i4)')'======> set rxCGit to',rxCGit
      endif
      if(index(string,"RXCGTH").ne.0)then                                !found RXCGTH keyword
       iindex = index(string,"RXCGTH")
       iindex = iindex + 7
       call skpblnks(string,iindex)
       read(string(iindex:),*)rxCGth
       write(iw,'(1x,a,f14.8)')'======> set RXCGTH to',rxCGth
      endif
      if(index(string,"REXLDL").ne.0)then                                !found REXLDL keyword
       iindex = index(string,"REXLDL")
       iindex = iindex + 7
       call skpblnks(string,iindex)
       read(string(iindex:),*)rexLdL
       write(iw,'(1x,a,i4)')'======> set rexLdL to',rexLdL
      endif
C
      write(iw,'(/1x,a)')'=====> End of REKS control string <====='
      if(.not.iropnd)close(unit=ir)                                      !zip it, if it was not opened
      if(abs(wpps-1.d0).le.1.d-6)then                                    !case of a single state REKS calculation; reset the type and target
      rexType = 0
      rexTarget = 0
      endif
C
 1000 continue
      endif                                                              !if(maswrk)then
C
      if(goparr)then                                                     !broadcast to slaves, if parallel
      call ddi_bcast(1234,'F',WPPS,1,MASTER)
      call ddi_bcast(2134,'F',rexShift,1,MASTER)
      call ddi_bcast(2314,'I',rexType,1,MASTER)
      call ddi_bcast(2341,'I',rexTarget,1,MASTER)
      call ddi_bcast(3241,'I',rexLdL,1,MASTER)
      call ddi_bcast(3421,'I',rexDIIS,1,MASTER)
      call ddi_bcast(1342,'I',rlxDen,1,MASTER)
      call ddi_bcast(2342,'I',rexEKT,1,MASTER)
      call ddi_bcast(3452,'I',EKTEA,1,MASTER)
      call ddi_bcast(2345,'I',rexCG,1,MASTER)
      call ddi_bcast(5234,'I',rxCGit,1,MASTER)
      call ddi_bcast(4523,'F',rxCGth,1,MASTER)
      endif                                                              !if(goparr)then
C
      return
      end subroutine rexinput
C*MODULE REKS    *DECK RMCOMMAS
C>
C> @author     Michael Filatov
C>
C> @date       2021, Oct
C>
      SUBROUTINE rmcommas(string,lenstr)                                 !removes commas from the string; replaced by blanks
      USE prec, ONLY: dp
      IMPLICIT NONE
C
      INTEGER :: LENSTR, I
C
      CHARACTER(256) :: LINE, STRING

      line = ' '
      do i = 1,lenstr
         if(string(i:i).ne.',') line(i:i) = string(i:i)
      enddo
      string = line
      return
      end subroutine rmcommas
C*MODULE REKS    *DECK SKPBLNKS
C>
C> @author     Michael Filatov
C>
C> @date       2021, Oct
C>
      SUBROUTINE skpblnks(string,iindex)                                 !skips blanks before "=" in the -string-; returns new -iindex-
      USE prec, ONLY: dp
      IMPLICIT NONE
C
      INTEGER :: IINDEX
      CHARACTER(256) :: STRING
C
      DO
      IF (STRING(IINDEX-1:IINDEX-1).EQ.' ') THEN
        IINDEX=IINDEX+1
        CYCLE
      END IF
      EXIT
      END DO
      RETURN
      END SUBROUTINE SKPBLNKS
C*MODULE REKS    *DECK CPREKS
C>
C> @author     Michael Filatov
C>
C> @date       2021, Oct
C>
      SUBROUTINE CPREKS(grdflg)                                          !grdflg - flag for gradient calculation; some quantities are needed only for gradient
C
      USE comm_SSR, ONLY: S2SARE,WS2SA, S3SARE, WS3SA, CLXGRD
      USE comm_REKSCM, ONLY: NMICRO, MTTYP, WPPS, WOSS, G1, DNR,
     * DNS, DELTA, FR, FS
      USE comm_REXOPT, ONLY: REXTYPE, REXTARGET, REXSHIFT, REXDIIS,
     * REXLDL, RLXDEN, REXEKT, EKTEA, REXCG, RXCGIT, RXCGTH
      USE camdft, ONLY: CAMFLAG
      USE lrcdft, ONLY: LCFLAG
      USE mx_limits, only: mxgrid, mxsh, mxatm
      USE constants, only: zero, one
      USE prec, ONLY: dp
      IMPLICIT NONE
C
      LOGICAL :: GRDFLG
      LOGICAL :: ALPHKWD, BETAKWD, DAMPXH, DFTB3, DFTBFL, DIRSCF,       &
     &           DSKWRK, GOPARR, LCDFTB, LRINT, MASWRK, MRDEA, MREKT,   &
     &           PACK2E, SCC, SG1, SG1T, SRSCC, TAMMD, TPA, TRIPLET
      REAL(KIND=dp), DIMENSION(137) :: BSLRD
      REAL(KIND=dp), DIMENSION(3,MXATM) :: C
      REAL(KIND=dp) :: CCTYP, CITYP, CNVTOL, DFTGTHR, DFTTHR,           &
     &                 DFTYPE, EXETYP, FDIFF , RHOMIN, RUNTYP, SCFTYP,  &
     &                 SW0, SWOFF, TDDFTYP, VBTYP
      INTEGER, DIMENSION(MXATM) :: IAN
      INTEGER :: IBTYP, ICH, IDAF, IGRDTYP, ILENG, INTTYP, IP, IPK,     &
     &           IPTIM, IR, IRECTD, IS, ITDFG, ITDPRP, IW, JANS, JANST, &
     &           LABSIZ, MASTER, MAXGRD, MAXVEC, ME, MODTD, MPCTYP,     &
     &           MPLEVL, MUL, MTHST, MULTD, NA, NAT, NAV, NB, NDFTFG,   &
     &           NE, NEVALS, NGLEVL, NHEX, NHLEVL, NINTMX, NLEBT,       &
     &           NONEQR, NPHI, NPHI0, NPHIT, NPROC, NQMT, NRAD, NRAD0,  &
     &           NRADT, NSTAT, NT, NTHE, NTHE0, NTHET, NTHST, NTRIAL,   &
     &           NTUPL, NUM
      INTEGER, DIMENSION(4) :: IFEDAT
      INTEGER, DIMENSION(48) :: INVT
      INTEGER, DIMENSION(950) :: IODA
      INTEGER, DIMENSION(MXATM,48) :: MAPCTR
      INTEGER, DIMENSION(MXSH,48) :: MAPSHL
      INTEGER, DIMENSION(MXGRID) :: NANGPT, NANGPT0, NLEB, NLEB0
      REAL(KIND=dp), DIMENSION(2) :: PFREQ
      REAL(KIND=dp), DIMENSION(3) :: SPCP
      REAL(KIND=dp), DIMENSION(432) :: T
      REAL(KIND=dp), DIMENSION(1) :: X
      REAL(KIND=dp), DIMENSION(MXATM) :: ZAN
      COMMON /DFGRID/ DFTTHR, DFTGTHR, SWOFF, SW0, BSLRD, NDFTFG, NRAD, &
     &                NTHE, NPHI, NRAD0, NTHE0, NPHI0, NANGPT, NANGPT0, &
     &                SG1, JANS
      COMMON /DFLEB / NLEB, NLEB0
      COMMON /DFTB  / DFTBFL, SCC, SRSCC, DFTB3, DAMPXH, LCDFTB
      COMMON /FMCOM / X
      COMMON /INFGRD/ RHOMIN, ILENG, MAXGRD
      COMMON /INFOA / NAT, ICH, MUL, NUM, NQMT, NE, NA, NB, ZAN, C, IAN
      COMMON /INFOTD/ CNVTOL, PFREQ, MODTD, JANST, NRADT, NTHET, NPHIT, &
     &                NLEBT, NSTAT, NTRIAL, MAXVEC, NTHST, IRECTD,      &
     &                ITDFG, ITDPRP, TRIPLET, SG1T, NONEQR, TAMMD, TPA, &
     &                ALPHKWD, BETAKWD, SPCP, MULTD, MREKT, MRDEA,      &
     &                MTHST, IFEDAT
      COMMON /INTFIL/ NINTMX, NHEX, NTUPL, PACK2E, INTTYP, IGRDTYP
      COMMON /IOFILE/ IR, IW, IP, IS, IPK, IDAF, NAV, IODA
      COMMON /NLRCF / LRINT
      COMMON /OPTSCF/ DIRSCF, FDIFF
      COMMON /PAR   / ME, MASTER, NPROC, IBTYP, IPTIM, GOPARR, DSKWRK,  &
     &                MASWRK
      COMMON /PCKLAB/ LABSIZ
      COMMON /RUNOPT/ RUNTYP, EXETYP, NEVALS, NGLEVL, NHLEVL
      COMMON /SYMTRY/ MAPSHL, MAPCTR, T, INVT, NT
      COMMON /WFNOPT/ SCFTYP, VBTYP, DFTYPE, TDDFTYP, CITYP, CCTYP,     &
     &                MPLEVL, MPCTYP
C
      REAL(KIND=dp) :: ATHRSH, CFT1, CFT2, CFT3, CFT4, CGTHRSH,         &
     &                 CLSA, CLX, CLXD, CORRFCTR, DFOCCA, DFOCCB, ELNEG,&
     &                 ELTOT, ERS, FCTR, FI, FJ, FOCCA, FOCCB, FONA,    &
     &                 FONB, OCCIA, OCCIB, OCCJA, OCCJB, PKAPK,         &
     &                 RKSQRD, RKZK, RKZK1, RSR, RSR12, SAIL, SQRTNR,   &
     &                 SQRTNS, THRESH, XMU, XTAU, XTHRSH, ZOMEGA,       &
     &                 ZPQ
      LOGICAL :: DBGAMAT, SOME
      LOGICAL, SAVE :: DBUG, RMDEEX
      REAL(KIND=dp) :: DDOT
      INTEGER :: I, IEND, II, IJ, IJ1, IJ2, IJIJ, IJRR, IJSR, IJSS, IK, &
     &           INDX, INFO, INFOF, INFOI, IORB, IRS, IRX, ISRART,      &
     &           ISTART, ISUM, ISX, ITERCG, J, JJ, JK, JRX, JSX, KK, KL,&
     &           KL1, KT, L, L0, L1, L2, L3, L7, LAMO, LAMOTEST, LAOMAX,&
     &           LAST, LBSMA, LBSMB, LBUF, LBVAO, LBVEC,                &
     &           LCMPPS, LCMSA, LCOEF, LDA, LDB, LDCH, LDM, LEA, LEB,   &
     &           LEC, LEC0, LEM, LEX, LEX0, LFAM, LFAMMO, LFAMTMP, LFBM,&
     &           LFBMMO, LFBMTMP, LFXC, LGMO, LGRD, LIAO, LIBUF, LINXOV,&
     &           LIPRGA, LIPRGB, LIRMOA, LIRMOB, LIWM
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IJCOPY                       !temporary array for orbital pairs indexing (dynamic allocation; works with ifort)
      REAL(8), ALLOCATABLE, DIMENSION(:) :: SCR1,SCR2,SCR3,SCR4,CMOSS
      INTEGER :: LJ2EA, LJ2EB, LK2EA, LK2EB, LOADFM, LOMATA, LOMATB,    &
     &           LOMG, LRHOI, LSCR, LSMAT,                              &
     &           LTA, LTAUI, LTB, LTRAI, LVA, LVALGA, LVB, LVHA, LVHB,  &
     &           LWGT, LWRK1, LWRK2, LWRK3, LX, LXM, MAXITCPRX,         &
     &           MXVEC, NCORE, NEED, NGOTMX, NLEBS, NMAX,               &
     &           NMIC, NNEG, NOCR, NOCS, NPHIS, NRADS, NRAL, NRBL, NRSV,&
     &           NSAL, NSBL, NTHES, NV, NVEC, NVIR, PQ
      CHARACTER(4) :: REXPREFIX
      LOGICAL :: ssr321,ssr320                                           !flags for SSR(3,2): SI3-SA2-REKS (PPS,OSS,+DES)
C
c       dimension ijcopy(MXAO*4)                                         !temporary array for orbital pairs indexing
C
      DATA DBUG/.false./
      DATA RMDEEX/.true./
c     DATA DBGAMAT/.true./                                               !set temporarily to .true. to enforce solving cp-reks by matrix inversion
c     DATA DBGAMAT/.false./                                              !set temporarily to .false. to enforce solving cp-reks by CG method
      REAL(KIND=dp), PARAMETER :: RNONE = transfer('NONE    ',1.0d0)
      REAL(KIND=dp), PARAMETER :: CHECK = transfer('CHECK   ',1.0d0)
C
      DBGAMAT = .not.rexCG                                               !use direct CP-REKS solver, if not CG solver
      maxitcprx = rxCGit                                                 !set max iterations for CP-REKS
      cgthrsh = rxCGth                                                   !set convergence criterion for CP-REKS
C
      if(rexTarget.eq.0)then
         if(MASWRK)write(iw,'(1x,"SA-REKS: no CP-REKS required")')
         RETURN
      endif
C
      ssr320 = .false.
      ssr321 = .false.
c     if(rexType.gt.1)then
      if(rexType.gt.2)then
        if(MASWRK)write(iw,'(1x,"CP-REKS for SSR > 3,2 NYI!")')
        call ABRT
      else
        if(MASWRK)then
          write(iw,'(1x,"==> Solving CP-REKS equations")')
          if(rexTarget.eq.2.and.rexType.eq.0)then
            write(iw,'(1x,"==> for OSS state")')
          elseif(rexTarget.eq.1.and.rexType.eq.0)then
            write(iw,'(1x,"==> for PPS state")')
          elseif(rexTarget.eq.2.and.rexType.eq.1)then
            write(iw,'(1x,"==> for SSR22 state S1")')
          elseif(rexTarget.eq.1.and.rexType.eq.1)then
            write(iw,'(1x,"==> for SSR22 state S0")')
          elseif(rexTarget.eq.2.and.rexType.eq.2)then
                ssr321 = .true.
            write(iw,'(1x,"==> for SSR32 state S1")')
          elseif(rexTarget.eq.1.and.rexType.eq.2)then
                ssr320 = .true.
            write(iw,'(1x,"==> for SSR32 state S0")')
          else
            write(iw,'(1x,"==> unknown target state,",i4)')rexTarget
            call ABRT
          endif
        endif
      endif
C
      SOME = MASWRK
C
C     ----- TURN OFF SYMMETRY -----
C
C
C     ----- PREPARE DFT GRID INFORMATION -----
C
      ITDFG=1
      IF(DFTYPE.NE.RNONE) NDFTFG = 1
C
      IF(NT.GT.1.AND..NOT.DFTBFL) THEN
         CALL SYMOFF
         IF(.NOT.DIRSCF) THEN
            IF(LCFLAG.OR.CAMFLAG) THEN
               LRINT=.TRUE.
               CALL JANDK
               LRINT=.FALSE.
            ENDIF
            CALL JANDK
         END IF
         CALL SYMON
      END IF
C
C     ----- CHANGE GRID SIZE
C
      IF(NDFTFG.EQ.1) THEN
         NRADS=NRAD
         NTHES=NTHE
         NPHIS=NPHI
         NLEBS=NLEB(1)
C
         NRADT  = 48
         NTHET  =  0
         NPHIT  =  0
         NLEBT  = 110
C
         NRAD=NRADT
         NTHE=NTHET
         NPHI=NPHIT
         NLEB(1)=NLEBT
C
         IF(.NOT.DFTBFL) THEN
            CALL SYMOFF
            CALL MEMGRD
            CALL SYMON
         END IF
      END IF
C here the parallel execution resumes; it appears that ssr320 and ssr321 need to be bradcasted to slaves
      if(goparr)then                                                     !broadcast to slaves, if parallel
         call ddi_bcast(3422,'I',ssr320,1,MASTER)
         call ddi_bcast(3423,'I',ssr321,1,MASTER)
      endif
C
C     ----- DEFINE # OF OCCUPIED AND VIRTUAL ORBITALS -----
C
      LX    = NQMT                                                       ! # OF MO's
      NCORE = NA-1                                                       ! # OF OCCUPIED CORE ORBITALS
      NOCR  = NA                                                         ! # OF ACTIVE R ORBITAL
      NOCS  = NA+1                                                       ! # OF ACTIVE S ORBITAL
      NVIR  = LX - NOCS                                                  ! # OF VIRTUAL  ORBITALS
      NRSV  = LX - NCORE                                                 ! # OF ACTIVE AND VIRTUAL ORBITALS
C
      L0 = NINTMX
      L1 = NUM
      L2 = (L1*(L1+1))/2
      L3 = L1*L1
C
      L7 = NOCS*NRSV                                                     !length of the Z vector
C
C     ----- DEFINE MXVEC & NVEC:NUMBER OF INITIAL VECTORS -----
C
      NMIC=4                                                             !number of REKS micristates

      IF(DBGAMAT) THEN
         NVEC=L7                                                         !number of "trial" vectors for direct method
         MXVEC = L7                                                      !max number of "trial" vectors
         IEND=L7
         ISTART=1
      ELSE
         NVEC = 1                                                        !only one vector is needed for CG solver
         MXVEC = 1                                                       !set MXVEC to 1, as only one vector is needed
c        MXVEC = L7                                                      !NEED TO CHECK THIS SETTING; IS IT REALLY NECESSARY?
         ISRART = 1
         IEND = 1
      ENDIF
C
      NMAX   = MAX(2,NVEC)                                               !for scratch arrays
C
      CALL GOTFM(NGOTMX)
      CALL VALFM(LOADFM)
C
C     ----- MEMORY ALLOCATION -----
C
      LVA    = LOADFM + 1
      LVB    = LVA
      LEA    = LVB    + L1*LX
      LEB    = LEA
      LINXOV = LEB    + LX
      LXM    = LINXOV + L7*2                                             ! X vector
      LDM    = LXM    + L7                                               ! XΔ vector
      if(ssr320.or.ssr321) LDM = LDM + L7                                ! book one more array for X.DES vector in LXM [X,X.DES]
      LBVEC  = LDM    + L7                                               ! trial vector in MO
      if(ssr320.or.ssr321) LBVEC = LBVEC + L7                            ! book one more array for XΔ12 vector in LDM [XΔ,XΔ12]
      LWRK1  = LBVEC  + L7*MXVEC                                         ! scratch 1
      LWRK2  = LWRK1  + L3*NMIC*NMAX                                     ! scratch 2
      LWRK3  = LWRK2  + L3*NMIC*NMAX                                     ! scratch 3
      LJ2EA  = LWRK3  + L3*4                                             ! scratch J2e
      LJ2EB  = LJ2EA  + L3*NMIC*NMAX                                     ! scratch J2e
      LK2EA  = LJ2EB  + L3*NMIC*NMAX                                     ! scratch K2e
      LK2EB  = LK2EA  + L3*NMIC*NMAX                                     ! scratch K2e
      LBVAO  = LK2EB  + L3*NMIC*NMAX                                     ! trial vector in AO
      LBSMA  = LBVAO  + L3*NMAX
      LBSMB  = LBSMA  + L3*NMIC*NMAX
      LFAM   = LBSMB  + L3*NMIC*NMAX                                     ! microstates alpha fock mat in AO
      LFBM   = LFAM   + L2*NMIC                                          ! microstates beta fock mat in AO
      LFAMMO = LFBM   + L2*NMIC                                          ! microstates alpha fock mat in MO
      LFBMMO = LFAMMO + LX*LX*NMIC                                       ! microstates beta fock mat in MO
      LEM    = LFBMMO + LX*LX*NMIC                                       ! microstates energy
      LCMSA  = LEM    + NMIC                                             ! microstates SA weight
      LIWM   = LCMSA  + NMIC                                             ! microstates Omega weight
      LCMPPS = LIWM   + NMIC                                             ! microstates PPS weight
      allocate(CMOSS(NMIC))                                              ! microstates OSS weight
      LOMG   = LCMPPS + NMIC                                             ! Omega
      LBUF   = LOMG   + L7
      LIBUF  = LBUF   + L0
      LAMO   = LIBUF  + L0*LABSIZ                                        ! A * Bvec
      LAMOtest= LAMO  + L7*MXVEC                                         ! tmp for A * Bvec
      LSCR  = LAMOtest+ L7*MXVEC
      LIRMOA = LSCR   + L7
      LIRMOB = LIRMOA + LX
      LDA    = LIRMOB + LX
      LDB    = LDA    + L2
      LTA    = LDB    + L2
      LTB    = LTA    + L2
      LSMAT  = LTB    + L2
      LOMATA = LSMAT  + L2                                               ! (n^{alpha}_p - n^{alpha}_q)/2.0D+00
      LOMATB = LOMATA + LX*LX*NMIC                                       ! (n^{beta }_p - n^{beta }_q)/2.0D+00
      LAST   = LOMATB + LX*LX*NMIC

C MEMORY FOR DFT GRID
      IF(NDFTFG.EQ.1) THEN
         LGRD   = LAST                                                   ! XYZGRD(MAXGRD*3)
         LWGT   = LGRD   + MAXGRD*3                                      ! XYZWGT(MAXGRD)
         LDCH   = LWGT   + MAXGRD                                        ! XYZDCH
         LRHOI  = LDCH   + 4*L1*ILENG                                    ! DRHOI_(MAXGRD,4,2)
         LTAUI  = LRHOI  + 8*MAXGRD                                      ! DTAUI_(MAXGRD,2)
         LTRAI  = LTAUI  + 2*MAXGRD                                      ! TRAI (5,2,ILENG)
         LCOEF  = LTRAI  + 10*ILENG                                      ! COEF (8,4,ILENG)
         LEX    = LCOEF  + 32*ILENG                                      ! EX (ILENG,18), IF 2ND DERIVS
         LEX0   = LEX    + 18*ILENG                                      ! EX0(ILENG)
         LEC    = LEX0   + ILENG                                         ! EC (ILENG,35), IF 2ND DERIVS
         LEC0   = LEC    + 35*ILENG                                      ! EC0(ILENG)
         LAOMAX = LEC0   + ILENG                                         ! AOMAX(L1)
         LVHA   = LAOMAX + L1                                            ! VHA(L1,L1)
         LVHB   = LVHA   + L3                                            ! VHB(L1,L1)
         LGMO   = LVHB   + L3                                            ! GMO (L1,4,2) FOR ALPHA AND BETA
         LVALGA = LGMO   + L1*8                                          ! VALGA (L1,2) FOR ALPHA AND BETA
         LFXC   = LVALGA + L1*2                                          ! FXC(L2,2)
         LIAO   = LFXC   + L2*2                                          ! IAO(L1)
         LIPRGA = LIAO   + L1                                            ! IPRGRDA(2,L3)
         LIPRGB = LIPRGA + L3*2                                          ! IPRGRDB(2,L3)
         LAST   = LIPRGB + L3*2
      END IF
C
      NEED = LAST - LOADFM - 1

      IF(DBUG.AND.MASWRK) THEN
         write (iw,*) 'L0', L0
         write (iw,*) 'L1', L1
         write (iw,*) 'L2', L2
         write (iw,*) 'L3', L3
         write (iw,*) 'LX', LX
         write (iw,*) 'L7', L7
         write (iw,*) 'MXVEC', MXVEC
         write (iw,*) 'NMIC', NMIC
         write (iw,*) 'NMAX', NMAX
         write (iw,*) 'LABSIZ', LABSIZ
         write (iw,*) 'MAXGRD', MAXGRD
         write (iw,*) 'ILENG', ILENG
      END IF

      IF(MASWRK)
     * WRITE(IW,'(1x,"==> Memory required for CP-REKS is",i10
     *           ," words")')NEED
      CALL GETFM(NEED)
C
      IF(EXETYP.NE.CHECK) THEN                                           !go to line 6259
C
C     ----- READ MO AND ORBITAL ENERGIES -----
C
      CALL DAREAD(IDAF,IODA,X(LVA),L1*LX,15,0)                           !α-orbitals
      CALL DAREAD(IDAF,IODA,X(LEA),LX   ,17,0)                           !α-energies
      call dcopy(L1*LX,X(LVA),1,X(LVB),1)                                !copy to β-orbitals
      call dcopy(LX,X(LEA),1,X(LEB),1)                                   !and energies
C
      IF(DBUG) THEN
         WRITE(IW,*) 'ORBITAL ENERGIES, ALPHA/BETA'
         DO IORB=1,LX
            WRITE(IW,*) IORB,X(LEA+IORB-1),X(LEB+IORB-1)
         END DO
      END IF
C
C     ----- READ REKS INFORMATION -----
C
      CALL DAREAD(IDAF,IODA,X(LFAM) ,L2*NMIC,700,0)                      !α-Fock matrices; all 4
      CALL DAREAD(IDAF,IODA,X(LFBM) ,L2*NMIC,701,0)                      !β-Fock matrices
      CALL DAREAD(IDAF,IODA,X(LEM)  ,NMIC   ,702,0)                      !energies of microstates
      CALL DAREAD(IDAF,IODA,X(LCMSA),NMIC   ,703,0)                      !SA-REKS weighting factors
C SCRATCH MEMORY
      allocate(scr1(l3), scr2(l3), scr3(l3), scr4(l3))
C======================================================================================================================
C       build the right hand side of REKS response equations; the X-vector
C======================================================================================================================
C
      CALL TSRFAO2MO(X(LFAMMO),X(LFBMMO),X(LFAM),X(LFBM),X(LVA),         !transform Fock matrices to MO representation
     *               X(LWRK1),X(LWRK2),L1,L2,L3,LX,NMIC)
C
      CALL REXCM(X(LCMPPS),DNR,DNS,DELTA,1.d0,0.d0,NMIC)                 !do PPS weighting factors
      CMOSS(1:4)   = (/ 0.d0, 0.d0, 1.d0, -0.5d0 /)                      !do OSS weighting factors
      CALL VCLR(X(LXM),1,L7)                                             !clear array for X-vector
        if(ssr320.or.ssr321)CALL VCLR(X(LXM+L7),1,L7)                    !clear array for X.DES vector
          DO L=1,NMIC                                                    !loop over microstates
             LFAMTMP=LFAMMO + (L-1)*LX*LX                                !pointer to L-th microstate Fα
             LFBMTMP=LFBMMO + (L-1)*LX*LX                                !pointer to L-th microstate Fβ
c            WRITE(IW,'(1X,I1,A)') L,'TH MICROSTATE'
c            WRITE(IW,*) 'ALPHA SKELETON FOCK MATRIX'
c            CALL PRSQ(X(LFAMTMP),L1,L1,L1)
c            WRITE(IW,*) 'BETA SKELETON FOCK MATRIX'
c            CALL PRSQ(X(LFBMTMP),L1,L1,L1)
        nraL = nmicro(1,1,L) - 1                                         !occupation of -r- orbital in L-th microstate; α-spin
        nrbL = nmicro(1,2,L) - 1                                         !occupation of -r- orbital in L-th microstate; β-spin
        nsaL = nmicro(2,1,L) - 1                                         !occupation of -s- orbital in L-th microstate; α-spin
        nsbL = nmicro(2,2,L) - 1                                         !occupation of -s- orbital in L-th microstate; β-spin
        if(rexTarget.eq.2)then
        CLX = CMOSS(L)                                                   !weighting factor of L-th microstate in OSS state
        elseif(rexTarget.eq.1)then
        CLX = X(LCMPPS+L-1)                                              !weighting factor of L-th microstate in PPS state
        else
        write(iw,'(1x,"==> Unknown target state,",i4)')rexTarget
        call ABRT
        endif
        if(L.gt.2)CLX = CLX + CLX                                        !double the weights for BS and T microstates
        if(ssr320.or.ssr321)then
        CLXD = -2.d0*X(LCMPPS+L-1)                                       !invert the sign for DES; valid only for L=3,4
        if(L.eq.1)CLXD = X(LCMPPS+L)
        if(L.eq.2)CLXD = X(LCMPPS+L-2)
        endif
c               if(MASWRK)write(iout,*)CLX,CLXD
        ij = 0
        do j = NOCR,LX                                                   !loop over orbital pairs, i<j; r,s,virt first
         occja = 0.d0
         if(j.eq.NOCR)occja = dble(nraL)
         if(j.eq.NOCS)occja = dble(nsaL)
         occjb = 0.d0
         if(j.eq.NOCR)occjb = dble(nrbL)
         if(j.eq.NOCS)occjb = dble(nsbL)
          do i = 1,NOCS                                                  !now core,r,s orbitals
            occia = 1.d0
            if(i.eq.NOCR)occia = dble(nraL)
            if(i.eq.NOCS)occia = dble(nsaL)
            occib = 1.d0
            if(i.eq.NOCR)occib = dble(nrbL)
            if(i.eq.NOCS)occib = dble(nsbL)
            ij = ij + 1                                                  !indexing X vector
            pq = (j - 1)*LX + i - 1                                      !indexing Fock matrices
            X(LXM+ij-1) = X(LXM+ij-1)                                    !X_{ij} = iε_{ij} - jε_{ij}
     *    + CLX * (occia - occja) * X(LFAMTMP + pq)
     *    + CLX * (occib - occjb) * X(LFBMTMP + pq)
            if(j.eq.NOCR.and.i.eq.NOCS)X(LXM+ij-1) = 0.d0                !zero -sr- element
        if(ssr320.or.ssr321)then                                         !for ssr32, build X.DES
            X(LXM+ij+L7-1) = X(LXM+ij+L7-1)                              !X_{ij} = iε_{ij} - jε_{ij}
     *    + CLXD * (occia - occja) * X(LFAMTMP + pq)
     *    + CLXD * (occib - occjb) * X(LFBMTMP + pq)
            if(j.eq.NOCR.and.i.eq.NOCS)X(LXM+ij+L7-1) = 0.d0             !zero -sr- element
        endif
          enddo
        enddo
          ENDDO                                                          !end of loop over microstates
c       if(MASWRK)then
c       write(iw,*)' test X-vector'
c       write(iw,'(/10(1x,F10.6))')(x(LXM+i-1),i=1,L7)
c       write(iw,*)' test X.DES-vector'
c       write(iw,'(/10(1x,F10.6))')(x(LXM+i+L7-1),i=1,L7)
c       endif
c       call ABRT
C=============== begin XDelta =========================================================================================
        if(rexType.gt.0)then                                             !begin building XΔ-vector for SSR; the vector is already multiplied by 2
        call VCLR(X(LOMATA),1,LX*LX)
        call VCLR(X(LOMATB),1,LX*LX)
        do L = 1,NMIC                                                    !build rε_{ij} and sε_{ij} matrices for building XΔ vector; use -LOMATA- and -LOMATB-
        CLSA = X(LCMSA+L-1)                                              !weights for the SA state
        if(L.gt.2)CLSA = CLSA + CLSA
             LFAMTMP=LFAMMO + (L-1)*LX*LX                                !pointer to L-th microstate Fα
             LFBMTMP=LFBMMO + (L-1)*LX*LX                                !pointer to L-th microstate Fβ
        nraL = nmicro(1,1,L) - 1                                         !occupation of -r- orbital in L-th microstate; α-spin
        nrbL = nmicro(1,2,L) - 1                                         !occupation of -r- orbital in L-th microstate; β-spin
        nsaL = nmicro(2,1,L) - 1                                         !occupation of -s- orbital in L-th microstate; α-spin
        nsbL = nmicro(2,2,L) - 1                                         !occupation of -s- orbital in L-th microstate; β-spin
      call daxpy(LX*LX,0.5d0*dble(nraL)*CLSA,X(LFAMTMP),1,X(LOMATA),1) !rε matrix, α-spins
      call daxpy(LX*LX,0.5d0*dble(nrbL)*CLSA,X(LFBMTMP),1,X(LOMATA),1) !rε matrix, β-spins
      call daxpy(LX*LX,0.5d0*dble(nsaL)*CLSA,X(LFAMTMP),1,X(LOMATB),1) !sε matrix, α-spins
      call daxpy(LX*LX,0.5d0*dble(nsbL)*CLSA,X(LFBMTMP),1,X(LOMATB),1) !sε matrix, β-spins
        enddo
c       write(iw,*)' ==> REpsilon'
c       call PRSQ(X(LOMATA),LX,LX,LX)
c       write(iw,*)' ==> SEpsilon'
c       call PRSQ(X(LOMATB),LX,LX,LX)
                                                                         !compute Lagrangian part of XΔ-vector
      CALL VCLR(X(LDM),1,L7)                                             !clear array for XΔ-vector
        sqrtNr = sqrt(2.d0*DNR)
        sqrtNs = sqrt(2.d0*DNS)
        ij = 0
        do j = NOCR,LX                                                   !loop over orbital pairs, i<j; r,s,virt first
          do i = 1,NOCS                                                  !now core,r,s orbitals
            ij = ij + 1                                                  !indexing XΔ-vector
            if(i.eq.NOCR.and.i.lt.j)then                                 !case δ_pr in eq. (55)
            jsx = (j - 1)*LX + NOCS - 1                                  !pointer to rε and sε elements (q,s); 0.5 is already included in ε
            X(LDM+ij-1) = X(LDM+ij-1) - sqrtNr*X(LOMATA+jsx)
     *                                + sqrtNs*X(LOMATB+jsx)
            endif
            if(i.eq.NOCS.and.i.lt.j)then                                 !case δ_ps in eq. (55)
            jrx = (j - 1)*LX + NOCR - 1                                  !elements (q,r)
            X(LDM+ij-1) = X(LDM+ij-1) - sqrtNr*X(LOMATA+jrx)
     *                                + sqrtNs*X(LOMATB+jrx)
            endif
            if(j.eq.NOCR.and.i.lt.j)then                                 !case δ_qr in eq. (55)
            isx = (NOCS - 1)*LX + i - 1                                  !elements (p,s)
            X(LDM+ij-1) = X(LDM+ij-1) + sqrtNr*X(LOMATA+isx)
     *                                - sqrtNs*X(LOMATB+isx)
            endif
            if(j.eq.NOCS.and.i.lt.j)then                                 !case δ_qs in eq. (55)
            irx = (NOCR - 1)*LX + i - 1                                  !elements (p,r)
            X(LDM+ij-1) = X(LDM+ij-1) + sqrtNr*X(LOMATA+irx)
     *                                - sqrtNs*X(LOMATB+irx)
            endif
          enddo
        enddo
C
        if(ssr320.or.ssr321)then                                         !for ssr32, build XΔ12
      CALL VCLR(X(LDM+L7),1,L7)                                          !clear array for XΔ12-vector
        ij = 0
        do j = NOCR,LX                                                   !loop over orbital pairs, i<j; r,s,virt first
          do i = 1,NOCS                                                  !now core,r,s orbitals
            ij = ij + 1                                                  !indexing XΔ-vector
            if(i.eq.NOCR.and.i.lt.j)then                                 !case δ_pr in eq. (55)
            jsx = (j - 1)*LX + NOCS - 1                                  !pointer to rε and sε elements (q,s)
            X(LDM+ij+L7-1) = X(LDM+ij+L7-1) - sqrtNr*X(LOMATA+jsx)
     *                                      - sqrtNs*X(LOMATB+jsx)
            endif
            if(i.eq.NOCS.and.i.lt.j)then                                 !case δ_ps in eq. (55)
            jrx = (j - 1)*LX + NOCR - 1                                  !elements (q,r)
            X(LDM+ij+L7-1) = X(LDM+ij+L7-1) - sqrtNr*X(LOMATA+jrx)
     *                                      - sqrtNs*X(LOMATB+jrx)
            endif
            if(j.eq.NOCR.and.i.lt.j)then                                 !case δ_qr in eq. (55)
            isx = (NOCS - 1)*LX + i - 1                                  !elements (p,s)
            X(LDM+ij+L7-1) = X(LDM+ij+L7-1) + sqrtNr*X(LOMATA+isx)
     *                                      + sqrtNs*X(LOMATB+isx)
            endif
            if(j.eq.NOCS.and.i.lt.j)then                                 !case δ_qs in eq. (55)
            irx = (NOCR - 1)*LX + i - 1                                  !elements (p,r)
            X(LDM+ij+L7-1) = X(LDM+ij+L7-1) + sqrtNr*X(LOMATA+irx)
     *                                      + sqrtNs*X(LOMATB+irx)
            endif
          enddo
        enddo
c       if(MASWRK)then
c       write(iw,*)' test XΔ-vector'
c       write(iw,'(/10(1x,F10.6))')(x(LDM+i-1),i=1,L7)
c       write(iw,*)' test XΔ12-vector'
c       write(iw,'(/10(1x,F10.6))')(x(LDM+i+L7-1),i=1,L7)
c       endif
c       call ABRT
        endif                                                            !if(ssr320.or.ssr321)then
C========== build ΔQ1 ===================================================
      CALL VCLR(SCR1,1,L3)                                               !clear array for ΔQ1 matrix for coupling term in SSR Lagrangian; it's miltiplied by -1
        ij = 0                                                           !it's needed anyway in EKT and in gradient
        do j = 1,LX                                                      !loop over orbital pairs, i<j; r,s,virt first
          do i = 1,LX                                                    !now core,r,s orbitals
            ij = ij + 1                                                  !indexing ΔQ1 matrix; the matrix is symmetrized
            if(i.eq.NOCR)then                                            !case δ_pr in eq. (55)
            jsx = (j - 1)*LX + NOCS - 1                                  !pointer to rε and sε elements (q,s); 0.5 is already in ε
            SCR1(ij) = SCR1(ij) - sqrtNr*X(LOMATA+jsx)
     *                          + sqrtNs*X(LOMATB+jsx)
            endif
            if(i.eq.NOCS)then                                            !case δ_ps in eq. (55)
            jrx = (j - 1)*LX + NOCR - 1                                  !elements (q,r)
            SCR1(ij) = SCR1(ij) - sqrtNr*X(LOMATA+jrx)
     *                          + sqrtNs*X(LOMATB+jrx)
            endif
            if(j.eq.NOCR)then                                            !case δ_qr in eq. (55)
            isx = (NOCS - 1)*LX + i - 1                                  !elements (p,s)
            SCR1(ij) = SCR1(ij) - sqrtNr*X(LOMATA+isx)
     *                          + sqrtNs*X(LOMATB+isx)
            endif
            if(j.eq.NOCS)then                                            !case δ_qs in eq. (55)
            irx = (NOCR - 1)*LX + i - 1                                  !elements (p,r)
            SCR1(ij) = SCR1(ij) - sqrtNr*X(LOMATA+irx)
     *                          + sqrtNs*X(LOMATB+irx)
            endif
          enddo
        enddo
        call dscal(LX*LX,0.5d0,SCR1,1)                                   !scale by 0.5; this is from symmetrization
c       write(iw,*)' ==> ΔQ1 matrix'
c       call PRSQ(SCR1,LX,LX,LX)
        call DAWRIT(IDAF,IODA,SCR1,LX*LX,705,0)                          !save ΔQ1 matrix to 705 for further use
C
C========== build Δ12Q1 =================================================
        if(ssr320.or.ssr321)then                                         !for ssr32, build Δ12Q1 matrix
      CALL VCLR(SCR1,1,L3)                                               !clear array for Δ12Q1 matrix; it's miltiplied by -1
        ij = 0
        do j = 1,LX                                                      !loop over orbital pairs, i<j; r,s,virt first
          do i = 1,LX                                                    !now core,r,s orbitals
            ij = ij + 1                                                  !indexing Δ12Q1 matrix; the matrix is symmetrized
            if(i.eq.NOCR)then                                            !case δ_pr in eq. (55)
            jsx = (j - 1)*LX + NOCS - 1                                  !pointer to rε and sε elements (q,s); 0.5 is already in ε
            SCR1(ij) = SCR1(ij) - sqrtNr*X(LOMATA+jsx)
     *                          - sqrtNs*X(LOMATB+jsx)
            endif
            if(i.eq.NOCS)then                                            !case δ_ps in eq. (55)
            jrx = (j - 1)*LX + NOCR - 1                                  !elements (q,r)
            SCR1(ij) = SCR1(ij) - sqrtNr*X(LOMATA+jrx)
     *                          - sqrtNs*X(LOMATB+jrx)
            endif
            if(j.eq.NOCR)then                                            !case δ_qr in eq. (55)
            isx = (NOCS - 1)*LX + i - 1                                  !elements (p,s)
            SCR1(ij) = SCR1(ij) - sqrtNr*X(LOMATA+isx)
     *                          - sqrtNs*X(LOMATB+isx)
            endif
            if(j.eq.NOCS)then                                            !case δ_qs in eq. (55)
            irx = (NOCR - 1)*LX + i - 1                                  !elements (p,r)
            SCR1(ij) = SCR1(ij) - sqrtNr*X(LOMATA+irx)
     *                          - sqrtNs*X(LOMATB+irx)
            endif
          enddo
        enddo
        call dscal(LX*LX,0.5d0,SCR1,1)                                   !scale by 0.5; this is from symmetrization
c       write(iw,*)' ==> Δ12Q1 matrix'
c       call PRSQ(SCR1,LX,LX,LX)
        call DAWRIT(IDAF,IODA,SCR1,LX*LX,708,0)                          !save Δ12Q1 matrix to 708 for further use
        endif                                                            !if(ssr320.or.ssr321)then
C
        Ers = X(LOMATA + (NOCR - 1)*LX + NOCS - 1)                       !element (r,s) of rε (or sε)
        Rsr = 2.d0/sqrtNr*Ers                                            !Rsr first term; already multiplied by 2
        if(sqrtNs.gt.1.d-8) Rsr = Rsr + 2.d0/sqrtNs*Ers                  !cap the 1/n_s term; if n_s -> 0, limiting value of this term times G1 is zero
        if(ssr320.or.ssr321)then                                         !for ssr32, build Rsr12
        Rsr12 = 2.d0/sqrtNr*Ers                                          !Rsr12 first term; already multiplied by 2
        if(sqrtNs.gt.1.d-8) Rsr12 = Rsr12 - 2.d0/sqrtNs*Ers              !sign is different compared to Rsr
        endif
        endif                                                            !if(rexType.gt.0)then
C=============== end XDelta ===========================================================================================
C
        call VCLR(X(LOMATA),1,LX*LX)
        call VCLR(X(LOMATB),1,LX*LX)
C OCCUPUTAION NUMBER DIFFERENCE
      CALL TSROMAT(X(LOMATA),X(LOMATB),NCORE,NOCR,NOCS,LX,NMIC,RMDEEX)
C I_L WEIGHT
      CALL TDSRIM(X(LIWM),X(LEM),NMIC)
C
C=============== begin XDelta =========================================================================================
        if(rexType.gt.0)then                                             !continue building XΔ
        irs = (NOCR - 1)*LX + NOCS - 1                                   !pointer to (r,s) element
        do L = 1,NMIC
        nraL = nmicro(1,1,L) - 1                                         !occupation of -r- orbital in L-th microstate; α-spin
        nrbL = nmicro(1,2,L) - 1                                         !occupation of -r- orbital in L-th microstate; β-spin
        nsaL = nmicro(2,1,L) - 1                                         !occupation of -s- orbital in L-th microstate; α-spin
        nsbL = nmicro(2,2,L) - 1                                         !occupation of -s- orbital in L-th microstate; β-spin
        LFAMTMP=LFAMMO + (L-1)*LX*LX                                     !pointer to L-th microstate Fα
        LFBMTMP=LFBMMO + (L-1)*LX*LX                                     !pointer to L-th microstate Fβ
        SAIL = X(LIWM+L-1)                                               !I_L coefficient for L-th microstate
        if(L.gt.2) SAIL = SAIL + SAIL                                    !double it for BS and T microstates
        Rsr = Rsr + WPPS * SAIL
     *      * (sqrtNr*dble(nraL)*X(LFAMTMP+irs)
     *      -  sqrtNs*dble(nsaL)*X(LFAMTMP+irs)
     *      +  sqrtNr*dble(nrbL)*X(LFBMTMP+irs)
     *      -  sqrtNs*dble(nsbL)*X(LFBMTMP+irs))
        if(ssr320.or.ssr321)then                                         !for ssr32, build Rsr12
        Rsr12 = Rsr12 + WPPS * SAIL
     *      * (sqrtNr*dble(nraL)*X(LFAMTMP+irs)
     *      +  sqrtNs*dble(nsaL)*X(LFAMTMP+irs)
     *      +  sqrtNr*dble(nrbL)*X(LFBMTMP+irs)
     *      +  sqrtNs*dble(nsbL)*X(LFBMTMP+irs))
        endif
        enddo
        endif                                                            !if(rexType.gt.0)then
C=============== end XDelta ===========================================================================================
C OMEGA MATRIX
      CALL VCLR(X(LOMG),1,L7)
      CALL TSRWSUM(X(LOMG),X(LFAMMO),X(LFBMMO),LX,X(LOMATA),X(LOMATB),
     *             X(LIWM),NCORE,NOCS,LX,NMIC,L7)
C
C=============== begin XDelta =========================================================================================
      if(rexType.gt.0)then                                               !continue with XΔ
      call daxpy(L7,G1*Rsr,X(LOMG),1,X(LDM),1)                           !add G1*Rsr*Ω_pq term to XΔ-vector
        if(ssr320.or.ssr321)                                             !for ssr32
     *call daxpy(L7,G1*Rsr12,X(LOMG),1,X(LDM+L7),1)                      !add G1*Rsr12*Ω_pq term to XΔ12-vector
c       if(MASWRK)then
c       write(iw,*)' ==> test X'
c       write(iw,'(/10(1x,F10.6))')(x(LXM+i-1),i=1,L7)
c       write(iw,*)' ==> test XΔ.1e'
c       write(iw,'(/10(1x,F10.6))')(x(LDM+i-1),i=1,L7)
c       write(iw,*)' ==> test XΔ12.1e'
c       write(iw,'(/10(1x,F10.6))')(x(LDM+i+L7-1),i=1,L7)
c       endif
c       call ABRT
C LXM - X vector, LXM+L7 - X.DES vector, LDM - XΔ vector, LDM+L7 - XΔ12 vector
      call SSRX2E(X(LDM),X(LCMSA),X(LWRK1),X(LWRK2),X(LOMATA),X(LOMATB), !compute 2e part of XΔ-vector
     *            X(LBSMA),X(LBSMB),X(LVA),X(LBUF),X(LIBUF),SCR1,
     *     X(LGRD),X(LWGT),X(LDCH),X(LRHOI),X(LTAUI),X(LAOMAX),X(LGMO),  !arrays for DFT XC kernels
     *     X(LFXC),X(LTRAI),X(LCOEF),X(LEX),X(LEC),X(LEX0),X(LEC0),      !arrays for DFT XC kernels
     *     X(LIAO),X(LIPRGA),X(LIPRGB),                                  !arrays for DFT XC kernels
     *              L1,L2,L3,LX,L7,NOCR,NOCS,NCORE,
     *              NMIC,.false.)
c       if(MASWRK)then
c       write(iw,*)' ==> test XΔ.2e'
c       write(iw,'(/10(1x,F10.6))')(x(LDM+i-1),i=1,L7)
c       endif
        if(ssr320.or.ssr321)then                                         !for ssr32 !!! this logical switch interferes with parallel execution!?
      call SSRX2E(X(LDM+L7),X(LCMSA),X(LWRK1),X(LWRK2),X(LOMATA),        !compute 2e part of XΔ12-vector
     *     X(LOMATB),X(LBSMA),X(LBSMB),X(LVA),X(LBUF),X(LIBUF),SCR1,
     *     X(LGRD),X(LWGT),X(LDCH),X(LRHOI),X(LTAUI),X(LAOMAX),X(LGMO),  !arrays for DFT XC kernels
     *     X(LFXC),X(LTRAI),X(LCOEF),X(LEX),X(LEC),X(LEX0),X(LEC0),      !arrays for DFT XC kernels
     *     X(LIAO),X(LIPRGA),X(LIPRGB),                                  !arrays for DFT XC kernels
     *              L1,L2,L3,LX,L7,NOCR,NOCS,NCORE,
     *              NMIC,.true.)

c       if(MASWRK)then
c       write(iw,*)' ==> test XΔ12.2e'
c       write(iw,'(/10(1x,F10.6))')(x(LDM+L7+i-1),i=1,L7)
c       endif
        endif                                                            !if(ssr320.or.ssr321)then
        if(.not.(ssr320.or.ssr321))then                                  !SSR(2,2)
                                                                         !exploit the fact that X^PPS = -X^OSS; a_{1m}^2*X^PPS+a_{2m}^2*X^OSS-2*a_{1m}*a_{2m}*XΔ
        cft1 = s2sare(1,1)**2 - s2sare(2,1)**2                           !if X^PPS is computed; SSR state 0
        if(rexTarget.eq.2)cft1 = s2sare(2,2)**2 - s2sare(1,2)**2         !if X^OSS is computed; SSR state 1
        cft2 = -2.d0*s2sare(1,1)*s2sare(2,1)                             !coefficient before XΔ; if state 0 is requested; check out "minus" sign
        if(rexTarget.eq.2)cft2 = -2.d0*s2sare(1,2)*s2sare(2,2)           !if state 1 is requested
        call dscal(L7,cft1,X(LXM),1)                                     !scale X^PPS/OSS by cft1
        call daxpy(L7,cft2,X(LDM),1,X(LXM),1)                            !cft1 * X^PPS/OSS + cft2*XΔ
c       write(iw,*)' ==> test X.SSR'
c       write(iw,'(/10(1x,F10.6))')(x(LXM+i-1),i=1,L7)
        else
                                                                         !SSR(3,2); note X^PPS = -X^OSS
        cft1 = s3sare(1,1)**2 - s3sare(2,1)**2                           !if X^PPS is computed; SSR state 0
        if(rexTarget.eq.2)cft1 = s3sare(2,2)**2 - s3sare(1,2)**2         !if X^OSS is computed; SSR state 1
        cft2 = s3sare(3,1)**2                                            !coefficient before X.DES; state 0
        if(rexTarget.eq.2)cft2 = s3sare(3,2)**2                          !state 1
        cft3 = -2.d0*s3sare(1,1)*s3sare(2,1)                             !coefficient before XΔ; if state 0 is requested; check out "minus" sign
        if(rexTarget.eq.2)cft3 = -2.d0*s3sare(1,2)*s3sare(2,2)           !if state 1 is requested
        cft4 = -2.d0*s3sare(2,1)*s3sare(3,1)                             !coefficient before XΔ12; state 0
        if(rexTarget.eq.2)cft4 = -2.d0*s3sare(2,2)*s3sare(3,2)           !state 1
        call dscal(L7,cft1,X(LXM),1)                                     !scale X^PPS/OSS by cft1
        call daxpy(L7,cft2,X(LXM+L7),1,X(LXM),1)                         !cft1 * X^PPS/OSS + cft2*X.DES
        call daxpy(L7,cft3,X(LDM),1,X(LXM),1)                            !cft1 * X^PPS/OSS + cft2*X.DES + cft3*XΔ
        call daxpy(L7,cft4,X(LDM+L7),1,X(LXM),1)                         !cft1 * X^PPS/OSS + cft2*X.DES + cft3*XΔ + cft4*XΔ12
c       write(iw,*)' ==> test X.SSR32'
c       write(iw,'(/10(1x,F10.6))')(x(LXM+i-1),i=1,L7)
        endif                                                            !if(.not.(ssr320.or.ssr321))
      endif                                                              !if(rexType.gt.0)then
C=============== end XDelta ===========================================================================================
C
      CALL VCLR(X(LAMO),1,L7*MXVEC)
C
      ISTART  = 1
      IEND    = NVEC
C
C======================================================================================================================
C       solve CP-REKS response equations:
C======================================================================================================================
C
        allocate(ijcopy(1:(LX-NOCR+1)*NOCS))                             !allocate ijcopy array
        ij = 0
        isum = 0
        do j = NOCR,LX                                                   !prepare to remove zero rows and columns from Hessian A and Z and other vectors
          do i = 1,NOCS                                                  !the original orbital pair vector has pairs with i>=j;these should be removed
           ij = ij + 1
            ijcopy(ij) = 1
            if(i.ge.j)ijcopy(ij) = 0                                     !these orbital pairs have i>=j
           isum = isum + ijcopy(ij)                                      !length of new orbital pairs array (i<j); keep the -isum-
          enddo
        enddo
C
      IF(DBGAMAT) THEN                                                   !compute the whole A matrix from TSRABVEC
C======================================================================================================================
C       by matrix inversion: Z = A^-1 * X
C======================================================================================================================
      if(MASWRK)
     *write(iw,'(1x,"Solving CP-REKS equations by direct method")')

c        CALL VCLR(X(LBVEC),1,L7*L7)                                     !trial vector is a L7*L7 unit matrix
         CALL VCLR(X(LBVEC),1,L7*MXVEC)                                  !trial vector is a L7*L7 unit matrix
         IJRR = (NOCR-NCORE-1)*NOCS+NOCR
         IJSR = (NOCR-NCORE-1)*NOCS+NOCS
         IJSS = (NOCS-NCORE-1)*NOCS+NOCS

         DO IJ1=ISTART,IEND
         DO IJ2=ISTART,IEND
            IJIJ=(IJ1-1)*L7+IJ2
            X(LBVEC+IJIJ-1)=ZERO
            IF(RMDEEX) THEN
               IF(IJ1.EQ.IJ2 .AND. IJ1.NE.IJRR .AND.
     *            IJ1.NE.IJSR.AND. IJ1.NE.IJSS) X(LBVEC+IJIJ-1)=ONE
            ELSE
               IF(IJ1.EQ.IJ2 .AND. IJ1.NE.IJRR .AND.
     *            IJ1.NE.IJSS) X(LBVEC+IJIJ-1)=ONE
            ENDIF
         ENDDO
         ENDDO
C
c         write(iw,*) 'test bvec', iter
c         CALL PRSQ(X(LBVEC),L7,L7,L7)
C
      NV = IEND - ISTART + 1
C
C     ----- MATRIX-VECTOR PRODUCTS: A_{SA}*BVEC -----
C
      CALL TSRABVEC(0,X(LAMO),X(LBVEC),X(LCMSA),X(LWRK1),X(LWRK2),
     *              X(LOMATA),X(LOMATB),X(LBSMA),X(LBSMB),X(LVA),
     *              X(LOMG),X(LBVAO),X(LFAM),X(LFBM),X(LBUF),X(LIBUF),
     *              X(LAMOTEST),SCR1,SCR2,SCR3,SCR4,
     *              X(LJ2EA),X(LJ2EB),X(LK2EA),X(LK2EB),
     *              X(LGRD),X(LWGT),X(LDCH),X(LRHOI),X(LTAUI),
     *              X(LAOMAX),X(LGMO),X(LFXC),X(LTRAI),X(LCOEF),
     *              X(LEX),X(LEC),X(LEX0),X(LEC0),X(LVALGA),
     *              X(LIAO),X(LIPRGA),X(LIPRGB),
     *              L1,L2,L3,LX,L7,NOCR,NOCS,NCORE,
     *              ISTART,IEND,NMAX,NV,MXVEC,NMIC)
C
        call dcopy(L7*L7,X(LAMO),1,X(LAMOtest),1)
c         write(iw,*) 'test A'
c         CALL PRSQ(X(LAMOtest),L7,L7,L7)
C
        call dcopy(L7,X(LXM),1,X(LWRK1),1)
        CALL VCLR(X(LXM),1,L7)
        CALL VCLR(X(LAMO),1,L7*L7)
        ij1 = 0                                                          !pointer to rows of "clean" A and X matrices
        do ij = 1,L7                                                     !now, remove zero rows and columns from A and X;loop over rows
         if(ijcopy(ij).eq.1)then
          ij1 = ij1 + 1
          X(LXM+ij1-1) = X(LWRK1+ij-1)                                   !move proper (i<j) elements of X from -LWRK1- to -LXM-
         endif
          kl1 = 0                                                        !pointer to columns of "clean" A matrix
          do kl = 1,L7
           if(ijcopy(kl).eq.1)kl1 = kl1 + 1
            if(ijcopy(ij).eq.1.and.ijcopy(kl).eq.1)
     *      X(LAMO+(kl1-1)*isum+ij1-1) = X(LAMOtest+(kl-1)*L7+ij-1)      !move proper (i<j;k<l) elements of A from -LAMOtest- to -LAMO-
          enddo
        enddo
c       write(iw,'(/10(1x,F10.6))')(x(LXM+i-1),i=1,isum)
c         write(iw,*) 'test A'
c         CALL PRSQ(X(LAMO),isum,isum,isum)
C
        call DGETRF(isum,isum,X(LAMO),isum,ijcopy,INFOF)                 !make LU factorization of A
        call DGETRI(isum,X(LAMO),isum,ijcopy,X(LWRK1),isum,INFOI)        !inverse A
        call dgemm('n','n',isum,1,isum,1.d0,X(LAMO),isum,X(LXM),isum,    !multiply A^-1 * X in -LWRK1-; this is the Z-vector
     *             0.d0,X(LWRK1),isum)
c       write(iw,'(/10(1x,F10.6))')(x(LWRK1+i-1),i=1,isum)
C
        ij = 0
        ij1 = 0
        CALL VCLR(X(LBVEC),1,L7)                                         !clear array for expanded Z-vector
        do j = NOCR,LX                                                   !expand the Z-vector back to the initial ordering; with i>=j pairs
          do i = 1,NOCS
           if(i.ge.j)then                                                !these orbital pairs have i>=j
            ij1 = ij1 + 1                                                !increment pointer for the i>=j pairs
           else
            ij = ij + 1                                                  !increment pointer for i<j pairs
            ij1 = ij1 + 1                                                !increment pointer for the i>=j pairs
            X(LBVEC+ij1-1) = X(LWRK1+ij-1)                               !copy from compressed to expanded Z-vector
           endif
          enddo
        enddo
C
      ELSE                                                               !IF(DBGAMAT) THEN; solve by CG method
C======================================================================================================================
C       by conjugate gradient method
C======================================================================================================================
C
      if(MASWRK)
     *write(iw,'(1x,"Solving CP-REKS equations by CG method",
     *          /1x,"convergence threshold",F18.12," within maximum of",
     *           i6," iterations")')cgthrsh,maxitcprx
      xthrsh = 1.d-16                                                    !threshold value for small denominators in GC method
      athrsh = 1.d-12                                                    !threshold value for small A1e.diag element
      itercg = 0                                                         !set iteration counter
      NV = 1                                                             !set the number of trial vectors (1)
c     CALL VCLR(X(LBVEC),1,L7*L7)                                        !clear the whole array -LBVEC-
      CALL VCLR(X(LBVEC),1,L7*MXVEC)                                     !clear the whole array -LBVEC-
C
      CALL VCLR(X(LSCR),1,L7)                                            !clear array for diagonal part of 1e A matrix; A^{1e}_{ij,ij}
      DO L=1,NMIC                                                        !loop over microstates
             LFAMTMP=LFAMMO + (L-1)*LX*LX                                !pointer to L-th microstate Fα
             LFBMTMP=LFBMMO + (L-1)*LX*LX                                !pointer to L-th microstate Fβ
        nraL = nmicro(1,1,L) - 1                                         !occupation of -r- orbital in L-th microstate; α-spin
        nrbL = nmicro(1,2,L) - 1                                         !occupation of -r- orbital in L-th microstate; β-spin
        nsaL = nmicro(2,1,L) - 1                                         !occupation of -s- orbital in L-th microstate; α-spin
        nsbL = nmicro(2,2,L) - 1                                         !occupation of -s- orbital in L-th microstate; β-spin
        CLX = X(LCMSA+L-1)*0.5d0                                         !weighting factor of L-th microstate in SA state; factor 1/2 included
        if(L.gt.2)CLX = CLX + CLX                                        !double the weights for BS and T microstates
        ij = 0
        do j = NOCR,LX                                                   !loop over orbital pairs, i<j; r,s,virt first
         occja = 0.d0
         if(j.eq.NOCR)occja = dble(nraL)
         if(j.eq.NOCS)occja = dble(nsaL)
         occjb = 0.d0
         if(j.eq.NOCR)occjb = dble(nrbL)
         if(j.eq.NOCS)occjb = dble(nsbL)
          do i = 1,NOCS                                                  !now core,r,s orbitals
            occia = 1.d0
            if(i.eq.NOCR)occia = dble(nraL)
            if(i.eq.NOCS)occia = dble(nsaL)
            occib = 1.d0
            if(i.eq.NOCR)occib = dble(nrbL)
            if(i.eq.NOCS)occib = dble(nsbL)
            ij = ij + 1                                                  !indexing X vector
            jj = (j - 1)*LX + j - 1                                      !<j|F|j>
            ii = (i - 1)*LX + i - 1                                      !<i|F|i>
            X(LSCR+ij-1) = X(LSCR+ij-1)                                  !A_{ij,ij}^1e; Ω contribution is neglected here; it's a preconditioner
     *    + CLX * (occia - occja) * (X(LFAMTMP + jj) - X(LFAMTMP + ii))
     *    + CLX * (occib - occjb) * (X(LFBMTMP + jj) - X(LFBMTMP + ii))
            if(j.eq.NOCR.and.i.eq.NOCS)X(LSCR+ij-1) = 0.d0               !zero -sr- element
          enddo
        enddo
      ENDDO                                                              !end of loop over microstates
        call dcopy(L7,X(LXM),1,X(LBVEC),1)                               !copy X to Z; this is the starting guess for Z
        call dcopy(L7,X(LBVEC),1,X(LDM),1)                               !save Z in -LDM-; -LDM- is used to store Z_k
C
      CALL TSRABVEC(0,X(LAMO),X(LBVEC),X(LCMSA),X(LWRK1),X(LWRK2),
     *              X(LOMATA),X(LOMATB),X(LBSMA),X(LBSMB),X(LVA),
     *              X(LOMG),X(LBVAO),X(LFAM),X(LFBM),X(LBUF),X(LIBUF),
     *              X(LAMOTEST),SCR1,SCR2,SCR3,SCR4,
     *              X(LJ2EA),X(LJ2EB),X(LK2EA),X(LK2EB),
     *              X(LGRD),X(LWGT),X(LDCH),X(LRHOI),X(LTAUI),
     *              X(LAOMAX),X(LGMO),X(LFXC),X(LTRAI),X(LCOEF),
     *              X(LEX),X(LEC),X(LEX0),X(LEC0),X(LVALGA),
     *              X(LIAO),X(LIPRGA),X(LIPRGB),
     *              L1,L2,L3,LX,L7,NOCR,NOCS,NCORE,
     *              ISTART,IEND,NMAX,NV,MXVEC,NMIC)
        call daxpy(L7,-1.d0,X(LAMO),1,X(LBVEC),1)                        !r = X - A*X; X in -LBVEC-
        call dcopy(L7,X(LBVEC),1,X(LXM),1)                               !copy r in -LXM-
c         rksqrd = ddot(L7,X(LXM),1,x(LXM),1)                            !calculate ||rk||^2
        rksqrd = dot_product(x(lxm:lxm+l7-1),x(lbvec:lbvec+l7-1))
        do ij = 1,L7                                                     !multiply the current Z0 vector by A1e.diag^-1
         if(ijcopy(ij).eq.1)then                                         !use only non-zero elements
          fctr = sign(1.d0/athrsh,X(LSCR+ij-1))                          !if A1e.diag is below the threshold, cap it
          if(abs(X(LSCR+ij-1)).gt.athrsh) fctr = 1.d0/X(LSCR+ij-1)       !inverse of A1e.diag
          X(LBVEC+ij-1) = X(LBVEC+ij-1)*fctr                             !z0 elements
         endif
        enddo
c        rkzk = ddot(L7,X(LBVEC),1,x(LXM),1)                             !calculate r_k.z_k
        rkzk = dot_product(x(lxm:lxm+l7-1),x(lbvec:lbvec+l7-1))
      if(MASWRK)
     *write(iw,'(1x,"CG: iteration",i6,", residue",F18.12)')
     *      itercg,sqrt(rksqrd)
      if(sqrt(rksqrd).le.cgthrsh)then
      if(MASWRK)
     * write(iw,'(1x,"CG: Residue (",F18.12,
     *               ") is less than the threshold (",F18.12,")")')
     *       sqrt(rksqrd),cgthrsh
       goto 200
      endif
C
  100 itercg = itercg + 1
      CALL VCLR(X(LAMO),1,L7)
      CALL TSRABVEC(0,X(LAMO),X(LBVEC),X(LCMSA),X(LWRK1),X(LWRK2),
     *              X(LOMATA),X(LOMATB),X(LBSMA),X(LBSMB),X(LVA),
     *              X(LOMG),X(LBVAO),X(LFAM),X(LFBM),X(LBUF),X(LIBUF),
     *              X(LAMOTEST),SCR1,SCR2,SCR3,SCR4,
     *              X(LJ2EA),X(LJ2EB),X(LK2EA),X(LK2EB),
     *              X(LGRD),X(LWGT),X(LDCH),X(LRHOI),X(LTAUI),
     *              X(LAOMAX),X(LGMO),X(LFXC),X(LTRAI),X(LCOEF),
     *              X(LEX),X(LEC),X(LEX0),X(LEC0),X(LVALGA),
     *              X(LIAO),X(LIPRGA),X(LIPRGB),
     *              L1,L2,L3,LX,L7,NOCR,NOCS,NCORE,
     *              ISTART,IEND,NMAX,NV,MXVEC,NMIC)
        pkApk = ddot(L7,X(LBVEC),1,x(LAMO),1)                            !calculate p_k.(Ap_k)
        if(abs(pkApk).gt.xthrsh)then
        xmu = rkzk/pkApk
        else
      if(MASWRK)
     *  write(iw,'(1x,"CG: Warning! Small value p_k*(A*p_k) = ",
     *             G18.12)')pkApk
        endif
        call daxpy(L7,xmu,X(LBVEC),1,X(LDM),1)                           !Z_k = Z_{k-1} + μ * p_k
        call daxpy(L7,-xmu,X(LAMO),1,X(LXM),1)                           !r_k = r_{k-1} - μ * (A*p_k)
        rksqrd = ddot(L7,X(LXM),1,x(LXM),1)                              !calculate ||rk||^2
      if(MASWRK)
     *write(iw,'(1x,"CG: iteration",i6,", residue",F18.12)')
     *      itercg,sqrt(rksqrd)
      if(sqrt(rksqrd).le.cgthrsh)then                                    !convergence test
      if(MASWRK)
     * write(iw,'(1x,"CG: Residue (",F18.12,
     *               ") is less than the threshold (",F18.12,")")')
     *       sqrt(rksqrd),cgthrsh
       goto 200
      endif
        call VCLR(X(LWRK1),1,L7)                                         !clear -LWRK1- for z_k
        do ij = 1,L7                                                     !multiply the current r_k vector by A1e.diag^-1
         if(ijcopy(ij).eq.1)then                                         !use only non-zero elements
          fctr = sign(1.d0/athrsh,X(LSCR+ij-1))                          !if A1e.diag is below the threshold, cap it
          if(abs(X(LSCR+ij-1)).gt.athrsh) fctr = 1.d0/X(LSCR+ij-1)       !inverse of A1e.diag
          X(LWRK1+ij-1) = X(LXM+ij-1)*fctr                               !z_k elements
         endif
        enddo
        rkzk1 = ddot(L7,X(LWRK1),1,x(LXM),1)                             !calculate  current r_k.z_k
        if(abs(rkzk).gt.xthrsh)then
        xtau = rkzk1/rkzk                                                !τ = (r_k.z_k)/(r_{k-1}.z_{k-1})
        else
      if(MASWRK)
     *  write(iw,'(1x,"CG: Warning! Small value r_k * z_k = ",G18.12)')
     *  rkzk
        endif
        rkzk = rkzk1                                                     !copy current r_k*z_k to old r_{k-1}*z_{k-1}
        call daxpy(L7,xtau,X(LBVEC),1,X(LWRK1),1)                        !p_k = z_{k-1} + τ * p_{k-1}
        call dcopy(L7,X(LWRK1),1,X(LBVEC),1)
C
        if(itercg.lt.maxitcprx)goto 100
      if(MASWRK)
     *write(iw,'(1x,"CG: Convergence not reached after",i6," iterations"
     *     )')itercg
  200 continue                                                           !bail out the loop
        call dcopy(L7,X(LDM),1,X(LBVEC),1)                               !copy Z-vector in -LBVEC-
C
      ENDIF                                                              !IF(DBGAMAT) THEN
C
        call DAWRIT(IDAF,IODA,X(LBVEC),L7,706,0)                         !save Z-vector to 706
c       write(iw,*)' test Z-vector'
c       write(iw,'(/10(1x,F10.6))')(x(LBVEC+i-1),i=1,L7)
C======================================================================================================================
C       build full Lagrangian (with relaxation terms), for the gradient calculation
C======================================================================================================================
      if(grdflg.or.rexEKT)then                                           !build XQ matrix for the gradient (and EKT)
        call VCLR(X(LWRK2),1,LX*LX)                                      !clear -LWRK2- for Z-vector expanded in LX*LX matrix
        ij = 0
        do j = NOCR,LX                                                   !expand Z vector to square array; LX x LX
          do i = 1,NOCS
           ij = ij + 1
           indx = (j-1)*LX + i                                           !indexing the upper triangle (i<j)
           X(LWRK2+indx-1) = X(LBVEC+ij-1)
          enddo
        enddo
c         write(iw,*) 'test Z'
c         CALL PRSQ(X(LWRK2),L1,L1,L1)
        call VCLR(X(LWRK1),1,LX*LX)                                      !clear -LWRK1- for XQ1 matrix in MO representation LX*LX
      DO L=1,NMIC                                                        !loop over microstates
             LFAMTMP=LFAMMO + (L-1)*L3                                   !pointer to L-th microstate Fα
             LFBMTMP=LFBMMO + (L-1)*L3                                   !pointer to L-th microstate Fβ
        nraL = nmicro(1,1,L) - 1                                         !occupation of -r- orbital in L-th microstate; α-spin
        nrbL = nmicro(1,2,L) - 1                                         !occupation of -r- orbital in L-th microstate; β-spin
        nsaL = nmicro(2,1,L) - 1                                         !occupation of -s- orbital in L-th microstate; α-spin
        nsbL = nmicro(2,2,L) - 1                                         !occupation of -s- orbital in L-th microstate; β-spin
        CLX = X(LCMSA+L-1)                                               !weighting factor of L-th microstate in SA state
        if(L.gt.2)CLX = CLX + CLX                                        !double the weights for BS and T microstates
        ij = 0
        do j = 1,LX                                                      !loop over orbitals q
         occja = 0.d0
         if(j.lt.NOCR)occja = 1.d0
         if(j.eq.NOCR)occja = dble(nraL)
         if(j.eq.NOCS)occja = dble(nsaL)
         occjb = 0.d0
         if(j.lt.NOCR)occjb = 1.d0
         if(j.eq.NOCR)occjb = dble(nrbL)
         if(j.eq.NOCS)occjb = dble(nsbL)
          do i = 1,LX                                                    !loop over orbitals p
            occia = 0.d0
            if(i.lt.NOCR)occia = 1.d0
            if(i.eq.NOCR)occia = dble(nraL)
            if(i.eq.NOCS)occia = dble(nsaL)
            occib = 0.d0
            if(i.lt.NOCR)occib = 1.d0
            if(i.eq.NOCR)occib = dble(nrbL)
            if(i.eq.NOCS)occib = dble(nsbL)
            pq = (j - 1)*LX + i - 1                                      !indexing Z vector (expanded); p < q
            if(i.gt.j) pq = (i - 1)*LX + j - 1                           !Z: p > q
            Zpq = 0.5d0 * X(LWRK2+pq) * CLX
            if(i.gt.j)Zpq = -Zpq                                         !change sign for p > q
            focca = (occia - occja)*Zpq
            foccb = (occib - occjb)*Zpq
            dfocca = abs(focca)
            dfoccb = abs(foccb)
           do kt = 1,LX                                                  !loop over orbitals t
            kk = (kt - 1)*LX - 1
            ik = kk + i                                                  !pt pointer
            jk = kk + j                                                  !qt pointer
            if(dfocca.gt.1.d-12)
     *      X(LWRK1+ik) = X(LWRK1+ik) + focca * X(LFAMTMP + jk)
            if(dfoccb.gt.1.d-12)
     *      X(LWRK1+ik) = X(LWRK1+ik) + foccb * X(LFBMTMP + jk)
           enddo
          enddo
        enddo
      ENDDO                                                              !DO L=1,NMIC
c         write(iw,*) 'test Q1(unsymm)'
c         CALL PRSQ(X(LWRK1),LX,LX,LX)
        call SYMTRZE(X(LWRK1),LX,LX)                                     !symmetrize Q1
c         write(iw,*) 'test Q1(symm)'
c         CALL PRSQ(X(LWRK1),LX,LX,LX)
C
c     call VCLR(X(LWRK1),1,LX*LX)                                        !test Q2
      call SSRQ2E(X(LWRK1),X(LBVEC),X(LCMSA),X(LOMATA),X(LOMATB),        !add Q2 matrix
     *            X(LBSMA),X(LBSMB),X(LVA),X(LBUF),X(LIBUF),SCR1,
     *     X(LGRD),X(LWGT),X(LDCH),X(LRHOI),X(LTAUI),X(LAOMAX),X(LGMO),  !arrays for DFT XC kernels
     *     X(LFXC),X(LTRAI),X(LCOEF),X(LEX),X(LEC),X(LEX0),X(LEC0),      !arrays for DFT XC kernels
     *     X(LIAO),X(LIPRGA),X(LIPRGB),                                  !arrays for DFT XC kernels
     *            L1,L2,L3,LX,L7,NOCR,NOCS,NCORE,NMIC)
c         write(iw,*) 'test Q1+Q2'
c         CALL PRSQ(X(LWRK1),LX,LX,LX)
C
      endif                                                              !if(grdflg.and.(.not.rexEKT))then
C======================================================================================================================
C       compute relaxed density matrix
C======================================================================================================================
                                                                         !this part is valid for both SSR(2,2) and SSR(3,2)
        zomega = ddot(L7,X(LBVEC),1,x(LOMG),1)                           !take trace of the product Z * Ω
        CALL VCLR(X(LWRK3),1,L7)                                         !clear array for occ. number weighted Z-vector
        ij = 0
        do j = NOCR,LX                                                   !make occ. number weighted Z; occ. numbers for SA-REKS
         fj = 0.d0                                                       !assume vacant orbitals
         if(j.eq.NOCR) fj = fr                                           !if not, then R or S
         if(j.eq.NOCS) fj = fs
          do i = 1,NOCS
           fi = 1.d0                                                     !assume core orbitals
           if(i.eq.NOCR) fi = fr                                         !if not, R or S
           if(i.eq.NOCS) fi = fs
           ij = ij + 1
           X(LWRK3+ij-1) = X(LBVEC+ij-1)*(fi - fj)                       ! Z_{pq} * (f_p - f_q); p<q
          enddo
        enddo
c       write(iw,*)' fr, fs =',fr,fs
c       write(iw,'(/10(1x,F10.6))')(x(LWRK3+i-1),i=1,L7)
        call VCLR(X(LFAM),1,LX*LX)
        ij = 0
        do j = NOCR,LX                                                   !expand occ.number weighted Z to square array; LX x LX
          do i = 1,NOCS
           ij = ij + 1
           indx = (i-1)*LX + j                                           !indexing the lower triangle
           X(LFAM+indx-1) = X(LWRK3+ij-1)
          enddo
        enddo
c         write(iw,*) 'test Z'
c         CALL PRSQ(X(LFAM),LX,LX,LX)
        call VCLR(X(LFBM),1,L1*L1)                                       !clear array for scratch
        call dgemm('n','t',LX,L1,LX,1.d0,X(LFAM),LX,X(LVA),L1,0.d0,      !take the right product: Z * C^{\dagger}; -LVA- is a L1 x LX array (AO x MO)
     *             X(LFBM),LX)
        call VCLR(X(LFAM),1,L1*L1)                                       !clear array for o.n.w. Z transformed to AO rep.
        call dgemm('n','n',L1,L1,LX,1.d0,X(LVA),L1,X(LFBM),LX,0.d0,      !take the left product: C * [Z*C^{\dagger}]
     *             X(LFAM),L1)
        call SYMTRZE(X(LFAM),L1,L1)                                      !symmetrize the resulting contribution to relaxed den. mat.
        call dscal(L3,-2.d0,X(LFAM),1)                                   !scale by -2
c         write(iw,*) 'test Pr1'
c         CALL PRSQ(X(LFAM),L1,L1,L1)
C
        fctr = 2.d0*WPPS*G1*zomega                                       !add the second contribution to the relaxed den. mat.
        call dger(l1,l1,fctr,x(lva+l1*(na-1)),1,x(lva+l1*(na-1)),1,
     *          x(lfam),l1)                                              !take product afon*C_r*C_r^{\dagger} and add to -LFAM-
        call dger(l1,l1,-fctr,x(lva+l1*na),1,x(lva+l1*na),1,
     *          x(lfam),l1)                                              !take product bfon*C_s*C_s^{\dagger} and add to -LFAM-
c         write(iw,*) 'test Pr2'
c         CALL PRSQ(X(LFAM),L1,L1,L1)
C
                                                                         !for SSR(2,2) and SSR(3,2) the density matrix was tested by the dipole moment calculation
        if(rexType.gt.0)then                                             !add SSR contribution ΔP^n; this part is different for SSR(2,2) and SSR(3,2)
c        call VCLR(X(LFAM),1,L1*L1)                                      !test
         cft1 = 0.5d0*WPPS*((2.d0*DNR-1.d0)*sqrt(2.d0*DNR)
     *               -(2.d0*DNS-1.d0)*sqrt(2.d0*DNS))                    !first factor in ΔP^n: this is common for (2,2) and (3,2)
         cft2 = -2.d0*G1*Rsr                                             !second factor in ΔP^n
         call VCLR(X(LFBM),1,L1*L1)                                      !build in -LFBM-
         call dger(L1,L1,cft1,X(LVA+L1*(NA-1)),1,X(LVA+L1*NA),1,         !cft1 * C_r * C_s^{\dagger}
     *             X(LFBM),L1)
         call dger(L1,L1,cft1,X(LVA+L1*NA),1,X(LVA+L1*(NA-1)),1,         !cft1 * C_s * C_r^{\dagger}
     *             X(LFBM),L1)
         call dger(L1,L1,cft2,X(LVA+L1*(NA-1)),1,X(LVA+L1*(NA-1)),1,     !cft2 * C_r * C_r^{\dagger}
     *             X(LFBM),L1)
         call dger(L1,L1,-cft2,X(LVA+L1*NA),1,X(LVA+L1*NA),1,            !-cft2 * C_s * C_s^{\dagger}
     *             X(LFBM),L1)
      if(ssr320.or.ssr321)then                                           !this is for SSR(3,2)
        fctr = +2.d0*s3sare(1,1)*s3sare(2,1)                             !coefficient before ΔP^n; if state 0 is requested; check out the sign
        if(rexTarget.eq.2)fctr = +2.d0*s3sare(1,2)*s3sare(2,2)           !if state 1 is requested
         call dscal(L3,fctr,X(LFBM),1)                                   !scale ΔP^n
         call VADD(X(LFBM),1,X(LFAM),1,X(LFAM),1,L3)                     !add to relaxation contribution
         cft3 = 0.5d0*WPPS*((2.d0*DNR-1.d0)*sqrt(2.d0*DNR)               !second contribution to SSR(3,2) ΔP^n
     *               +(2.d0*DNS-1.d0)*sqrt(2.d0*DNS))                    !first factor in ΔP^n
         cft4 = -2.d0*G1*Rsr12                                           !second factor in ΔP^n
         call VCLR(X(LFBM),1,L1*L1)                                      !build in -LFBM-
         call dger(L1,L1,cft3,X(LVA+L1*(NA-1)),1,X(LVA+L1*NA),1,         !cft3 * C_r * C_s^{\dagger}
     *             X(LFBM),L1)
         call dger(L1,L1,cft3,X(LVA+L1*NA),1,X(LVA+L1*(NA-1)),1,         !cft3 * C_s * C_r^{\dagger}
     *             X(LFBM),L1)
         call dger(L1,L1,cft4,X(LVA+L1*(NA-1)),1,X(LVA+L1*(NA-1)),1,     !cft4 * C_r * C_r^{\dagger}
     *             X(LFBM),L1)
         call dger(L1,L1,-cft4,X(LVA+L1*NA),1,X(LVA+L1*NA),1,            !-cft4 * C_s * C_s^{\dagger}
     *             X(LFBM),L1)
        fctr = +2.d0*s3sare(2,1)*s3sare(3,1)                             !coefficient before ΔP^n; if state 0 is requested
        if(rexTarget.eq.2)fctr = +2.d0*s3sare(2,2)*s3sare(3,2)           !if state 1 is requested
         call dscal(L3,fctr,X(LFBM),1)                                   !scale ΔP^n
         call VADD(X(LFBM),1,X(LFAM),1,X(LFAM),1,L3)                     !add to relaxation contribution
      else                                                               !this is for SSR(2,2) only
        fctr = +2.d0*s2sare(1,1)*s2sare(2,1)                             !coefficient before ΔP^n; if state 0 is requested; check out the sign
        if(rexTarget.eq.2)fctr = +2.d0*s2sare(1,2)*s2sare(2,2)           !if state 1 is requested
         call dscal(L3,fctr,X(LFBM),1)                                   !scale ΔP^n
         call VADD(X(LFBM),1,X(LFAM),1,X(LFAM),1,L3)                     !add to relaxation contribution
      endif                                                              !if(ssr320.or.ssr321)then
c         write(iw,*) 'test ΔP^n'
c         CALL PRSQ(X(LFAM),L1,L1,L1)
c        call VCLR(X(LFAM),1,L1*L1)                                      !test
        endif                                                            !if(rexType.gt.0)then
C
      CALL DAREAD(IDAF,IODA,X(LDA),L2,16,0)                              !read in unrelaxed den. mat.
        call VCLR(X(LFBM),1,L1*L1)
        call EXPND(X(LDA),X(LFBM),L1,0)                                  !expand to square
        call VADD(X(LFAM),1,X(LFBM),1,X(LFAM),1,L3)                      !add relaxation contribution
c         write(iw,*) 'test P(relaxed)'
c         CALL PRSQ(X(LFAM),L1,L1,L1)
        call CPYSQT(X(LFAM),X(LDA),L1,1)                                 !compress to lower triangular
      CALL DAWRIT(IDAF,IODA,X(LDA),L2,16,0)                              !save REKS relaxed density matrix to record 16
      CALL DAWRIT(IDAF,IODA,X(LDA),L2,20,0)                              !and to 20
c       call ABRT
C======================================================================================================================
C       build new weignting factors including the relaxation terms and finish building the Lagrangian
C======================================================================================================================
        if(rexEKT.or.grdflg)then                                         !bypass this part, if only relaxed density was needed
        do L = 1,nmic                                                    !do new SA/SSR weighting factors for the gradient
        if(rexTarget.eq.2)then
        CLX = CMOSS(L)                                                   !weighting factor of L-th microstate in OSS state (SA)
        if(rexType.gt.0)                                                 !weighting factor of L-th microstate in S1 state (SSR)
     *  CLX = s2sare(1,2)**2 * X(LCMPPS+L-1)                             !PPS state
     *      + s2sare(2,2)**2 * CMOSS(L)                                  !OSS state
     *      - 2.d0*s2sare(2,2)*s2sare(1,2) * G1 * Rsr * X(LIWM+L-1)      !I_L coefficient for L-th microstate
        if(ssr321)then
        CLXD = -X(LCMPPS+L-1)                                            !invert the sign for DES; valid only for L=3,4
        if(L.eq.1)CLXD = X(LCMPPS+L)
        if(L.eq.2)CLXD = X(LCMPPS+L-2)
        CLX = s3sare(1,2)**2 * X(LCMPPS+L-1)                             !PPS state
     *      + s3sare(2,2)**2 * CMOSS(L)                                  !OSS state
     *      + s3sare(3,2)**2 * CLXD                                      !DES state
     *      - 2.d0*s3sare(2,2)*s3sare(1,2) * G1 * Rsr * X(LIWM+L-1)      !from Δ
     *      - 2.d0*s3sare(2,2)*s3sare(3,2) * G1 * Rsr12 * X(LIWM+L-1)    !from Δ12
        endif                                                            !if(ssr321)then
        elseif(rexTarget.eq.1)then
        CLX = X(LCMPPS+L-1)                                              !weighting factor of L-th microstate in PPS state (SA)
        if(rexType.gt.0)                                                 !weighting factor of L-th microstate in S0 state (SSR)
     *  CLX = s2sare(1,1)**2 * X(LCMPPS+L-1)                             !PPS state
     *      + s2sare(2,1)**2 * CMOSS(L)                                  !OSS state
     *      - 2.d0*s2sare(2,1)*s2sare(1,1) * G1 * Rsr * X(LIWM+L-1)      !I_L coefficient for L-th microstate
        if(ssr320)then
        CLXD = -X(LCMPPS+L-1)                                            !invert the sign for DES; valid only for L=3,4
        if(L.eq.1)CLXD = X(LCMPPS+L)
        if(L.eq.2)CLXD = X(LCMPPS+L-2)
        CLX = s3sare(1,1)**2 * X(LCMPPS+L-1)                             !PPS state
     *      + s3sare(2,1)**2 * CMOSS(L)                                  !OSS state
     *      + s3sare(3,1)**2 * CLXD                                      !DES state
     *      - 2.d0*s3sare(2,1)*s3sare(1,1) * G1 * Rsr * X(LIWM+L-1)      !from Δ
     *      - 2.d0*s3sare(2,1)*s3sare(3,1) * G1 * Rsr12 * X(LIWM+L-1)    !from Δ12
        endif                                                            !if(ssr320)then
        else
        write(iw,'(1x,"==> Unknown target state,",i4)')rexTarget
        call ABRT
        endif
          CLXgrd(L) = CLX
     *              + WPPS*G1*X(LIWM+L-1)*zomega                         !finish the new weighting factors
        enddo                                                            !do L = 1,nmic
c       write(iw,*)' ==> new C(L) coefficients'
c       write(iw,'(4(1x,F14.8))')(CLXgrd(i),i=1,nmic)
C
        call dscal(LX*LX,-1.d0,X(LWRK1),1)                               !scale Q1+Q2 by -1; Q1 & Q2 are the same for SSR22 and SSR32
c       call VCLR(X(LWRK1),1,LX*LX)                                      !test
          DO L=1,NMIC                                                    !loop over microstates
             LFAMTMP=LFAMMO + (L-1)*L3                                   !pointer to L-th microstate Fα
             LFBMTMP=LFBMMO + (L-1)*L3                                   !pointer to L-th microstate Fβ
        nraL = nmicro(1,1,L) - 1                                         !occupation of -r- orbital in L-th microstate; α-spin
        nrbL = nmicro(1,2,L) - 1                                         !occupation of -r- orbital in L-th microstate; β-spin
        nsaL = nmicro(2,1,L) - 1                                         !occupation of -s- orbital in L-th microstate; α-spin
        nsbL = nmicro(2,2,L) - 1                                         !occupation of -s- orbital in L-th microstate; β-spin
        CLX = CLXgrd(L)
        if(L.gt.2)CLX = CLX + CLX                                        !double the weights for BS and T microstates
        ij = 0
        do j = 1,LX                                                      !loop over orbitals |ket>
         occja = 0.d0
         if(j.lt.NOCR)occja = 1.d0
         if(j.eq.NOCR)occja = dble(nraL)
         if(j.eq.NOCS)occja = dble(nsaL)
         occjb = 0.d0
         if(j.lt.NOCR)occjb = 1.d0
         if(j.eq.NOCR)occjb = dble(nrbL)
         if(j.eq.NOCS)occjb = dble(nsbL)
          do i = 1,LX                                                    !loop over orbitals <bra|
            occia = 0.d0
            if(i.lt.NOCR)occia = 1.d0
            if(i.eq.NOCR)occia = dble(nraL)
            if(i.eq.NOCS)occia = dble(nsaL)
            occib = 0.d0
            if(i.lt.NOCR)occib = 1.d0
            if(i.eq.NOCR)occib = dble(nrbL)
            if(i.eq.NOCS)occib = dble(nsbL)
            pq = (j - 1)*LX + i - 1                                      !indexing Lagrangian matrices; Largangians of PPS and OSS are symmetrized
            X(LWRK1+pq) = X(LWRK1+pq)                                    !F_{pq} = 1/2*(qε_{pq} + pε_{pq})
     *    + 0.5d0*CLX * (occia + occja) * X(LFAMTMP + pq)
     *    + 0.5d0*CLX * (occib + occjb) * X(LFBMTMP + pq)
          enddo
        enddo
          ENDDO                                                          !DO L=1,NMIC
C
        if(rexType.gt.0)then                                             !for SSR
        call VCLR(X(LWRK2),1,LX*LX)
        call DAREAD(IDAF,IODA,X(LWRK2),LX*LX,705,0)                      !read ΔQ1 matrix from 705
      call SSRDQ2E(X(LWRK2),X(LCMSA),X(LOMATA),X(LOMATB),
     *             X(LBSMA),X(LBSMB),X(LVA),X(LBUF),X(LIBUF),SCR1,
     *     X(LGRD),X(LWGT),X(LDCH),X(LRHOI),X(LTAUI),X(LAOMAX),X(LGMO),  !arrays for DFT XC kernels
     *     X(LFXC),X(LTRAI),X(LCOEF),X(LEX),X(LEC),X(LEX0),X(LEC0),      !arrays for DFT XC kernels
     *     X(LIAO),X(LIPRGA),X(LIPRGB),                                  !arrays for DFT XC kernels
     *             L1,L2,L3,LX,L7,NOCR,NOCS,NCORE,NMIC,.false.)
c         write(iw,*) 'test ΔQ2'
c         CALL PRSQ(X(LWRK2),LX,LX,LX)
C
        if(.not.(ssr320.or.ssr321))then                                  !switch between ssr22 and ssr32
        fctr = -2.d0*s2sare(1,1)*s2sare(2,1)                             !ssr22: coefficient before ΔQ contribution to Lagrangian
        if(rexTarget.eq.2)fctr = -2.d0*s2sare(1,2)*s2sare(2,2)           !if state 1 is requested
        else                                                             !for ssr32, use s3sare matrix
        fctr = -2.d0*s3sare(1,1)*s3sare(2,1)
        if(rexTarget.eq.2)fctr = -2.d0*s3sare(1,2)*s3sare(2,2)
        endif                                                            !if(ssr320.or.ssr321)then
         call dscal(LX*LX,fctr,X(LWRK2),1)                               !scale ΔQ
         call daxpy(LX*LX,1.d0,X(LWRK2),1,X(LWRK1),1)                    !and add L.ΔQ to the Lagrangian
C
        if(ssr320.or.ssr321)then                                         !for ssr32, read Δ12Q1 matrix
        call VCLR(X(LWRK2),1,LX*LX)
        call DAREAD(IDAF,IODA,X(LWRK2),LX*LX,708,0)                      !read Δ12Q1 matrix from 708
      call SSRDQ2E(X(LWRK2),X(LCMSA),X(LOMATA),X(LOMATB),
     *             X(LBSMA),X(LBSMB),X(LVA),X(LBUF),X(LIBUF),SCR1,
     *     X(LGRD),X(LWGT),X(LDCH),X(LRHOI),X(LTAUI),X(LAOMAX),X(LGMO),  !arrays for DFT XC kernels
     *     X(LFXC),X(LTRAI),X(LCOEF),X(LEX),X(LEC),X(LEX0),X(LEC0),      !arrays for DFT XC kernels
     *     X(LIAO),X(LIPRGA),X(LIPRGB),                                  !arrays for DFT XC kernels
     *             L1,L2,L3,LX,L7,NOCR,NOCS,NCORE,NMIC,.true.)
c         write(iw,*) 'test Δ12Q2'
c         CALL PRSQ(X(LWRK2),LX,LX,LX)
        fctr = -2.d0*s3sare(2,1)*s3sare(3,1)
        if(rexTarget.eq.2)fctr = -2.d0*s3sare(2,2)*s3sare(3,2)
         call dscal(LX*LX,fctr,X(LWRK2),1)                               !scale ΔQ
         call daxpy(LX*LX,1.d0,X(LWRK2),1,X(LWRK1),1)                    !and add L.ΔQ to the Lagrangian
        endif                                                            !if(ssr320.or.ssr321)then
        endif                                                            !if(rexType.gt.0)then
C
c         write(iw,*) 'test L.MO'
c         CALL PRSQ(X(LWRK1),LX,LX,LX)
      CALL DAWRIT(IDAF,IODA,X(LWRK1) ,L3,704,0)                          !save full SA/SSR Lagrangian in MO to 704; this is for EKT
      call VCLR(X(LWRK2),1,LX*LX)
      CALL DGEMM('N','N',L1,L1,L1,1.d0,X(LVA),L1,X(LWRK1),L1,0.d0,       !here, a half transform is done; it's C*L.MO*C^†, not S*C*L.MO*C^†*S
     *           X(LWRK2),L1)
      CALL TRPOSQ(X(LWRK2),L1)
      CALL DGEMM('N','N',L1,L1,L1,-1.d0,X(LVA),L1,X(LWRK2),L1,0.d0,      !the Lagrangian in gamess should be negative (-)
     *           X(LWRK1),L1)
      CALL CPYSQT(X(LWRK1),X(LWRK2),L1,1)                                !pack square back to lower
      CALL DAWRIT(IDAF,IODA,X(LWRK2),L2,36,0)                            !save it to where it belongs
C======================================================================================================================
C       test part used for debugging; build relaxed density matrix alternative way, using CLXgrd, etc.
C======================================================================================================================
c       call VCLR(X(LWRK3),1,LX*LX)                                      !clear aray for diagonal part of dm
c       DO L=1,NMIC
c       nraL = nmicro(1,1,L) - 1                                         !occupation of -r- orbital in L-th microstate; α-spin
c       nrbL = nmicro(1,2,L) - 1                                         !occupation of -r- orbital in L-th microstate; β-spin
c       nsaL = nmicro(2,1,L) - 1                                         !occupation of -s- orbital in L-th microstate; α-spin
c       nsbL = nmicro(2,2,L) - 1                                         !occupation of -s- orbital in L-th microstate; β-spin
c       CLX = CLXgrd(L)
c       if(L.gt.2)CLX = CLX + CLX                                        !double the weights for BS and T microstates
c       do i = 1,LX
c        occia = 0.d0
c        if(i.lt.NOCR)occia = 1.d0                                       !one-electron DM; scaled by 2 later
c        if(i.eq.NOCR)occia = dble(nraL)
c        if(i.eq.NOCS)occia = dble(nsaL)
c        occib = 0.d0
c        if(i.lt.NOCR)occib = 1.d0
c        if(i.eq.NOCR)occib = dble(nrbL)
c        if(i.eq.NOCS)occib = dble(nsbL)
c        occi = 0.5d0*(occia + occib)
c        pq = (i - 1)*LX + i - 1                                         !indexing diagonal part
c       X(LWRK3+pq) = X(LWRK3+pq) + CLX * occi
c       enddo
c       ENDDO                                                            !DO L=1,NMIC
c       call VCLR(X(LFAM),1,LX*LX)
c       ij = 0
c       do j = NOCR,LX                                                   !expand occ.number weighted Z to square array; LX x LX
c        fj = 0.d0                                                       !assume vacant orbitals
c        if(j.eq.NOCR) fj = fr                                           !if not, then R or S
c        if(j.eq.NOCS) fj = fs
c         do i = 1,NOCS
c          fi = 1.d0                                                     !assume core orbitals
c          if(i.eq.NOCR) fi = fr                                         !if not, R or S
c          if(i.eq.NOCS) fi = fs
c          ij = ij + 1
c          indx = (i-1)*LX + j                                           !indexing the lower triangle
c          X(LFAM+indx-1) = X(LBVEC+ij-1)*(fi - fj)                      !occ. number weighted Z; occ. numbers for SA-REKS
c         enddo
c       enddo
c       if(rexType.gt.0)then                                             !for SSR
c       DRsr = sqrt(2.d0*DNR)*fr - sqrt(2.d0*DNS)*fs
c       fctr = -2.d0*s2sare(1,1)*s2sare(2,1)                             !coefficient before -sr- contribution from ΔR_L
c       if(rexTarget.eq.2)fctr = -2.d0*s2sare(1,2)*s2sare(2,2)           !if state 1 is requested
c       indx = (NOCS-1)*LX + NOCR - 1
c       X(LFAM+indx) = X(LFAM+indx) + fctr * DRsr
c       endif                                                            !if(rexType.gt.0)then
c       call SYMTRZE(X(LFAM),LX,LX)                                      !symmetrize the resulting contribution to relaxed den. mat.
c       call daxpy(LX*LX,-1.d0,X(LFAM),1,X(LWRK3),1)                     !scale by -1 and add to the diagonal DM
c         write(iw,*) 'test DM.1'
c         CALL PRSQ(X(LWRK3),LX,LX,LX)
C
c     call VCLR(X(LFBM),1,L1*L1)                                         !clear array for scratch
c     call dgemm('n','t',LX,L1,LX,2.d0,X(LWRK3),LX,X(LVA),L1,0.d0,       !take the right product: DM * C^{\dagger}; -LVA- is a L1 x LX array (AO x MO)
c    *           X(LFBM),LX)                                             !scaled by 2 here
c     call VCLR(X(LFAM),1,L1*L1)                                         !clear array for DM in AO rep.
c     call dgemm('n','n',L1,L1,LX,1.d0,X(LVA),L1,X(LFBM),LX,0.d0,        !take the left product: C * [DM*C^{\dagger}]
c    *           X(LFAM),L1)
c     call CPYSQT(X(LFAM),X(LDA),L1,1)                                   !compress to lower triangular
c     CALL DAWRIT(IDAF,IODA,X(LDA),L2,16,0)                              !save REKS relaxed density matrix to record 16
c     CALL DAWRIT(IDAF,IODA,X(LDA),L2,20,0)                              !and to 20
C
        endif                                                            !if(rexEKT.or.grdflg)then
        if(rexEKT) then                                                  !bail out, if no EKT requested
C======================================================================================================================
C       EKT: F * x = P * x * λ; F - Fock matrix; P - relaxed density matrix; λ - EKT ionization energies (negatives)
C======================================================================================================================
c       EKTEA=.true.                                                     !test setting
        if(EKTEA)then                                                    !do electron affinities (EA) via the EKT
c       if(MASWRK)then
c         write(iw,*) 'test P.AO'
c         CALL PRTRIL(X(LDA),L1)
c       endif

      CALL BUILDFOCK(X(LDA),X(LWRK1),X(LWRK2),X(LWRK3),X(LGRD),X(LBUF),
     *                                                 X(LIBUF),.false.) !build closed-shell Fock matrix with the current (relaxed) DM (in -LDA-)
c       if(MASWRK)then
c         write(iw,*) 'test F.AO'
c         CALL PRTRIL(X(LWRK1),L1)
c       endif
        call VCLR(X(LWRK2),1,L1*L1)                                      !clear -LWRK2-
        call EXPND(X(LWRK1),X(LWRK2),L1,0)                               !expand to full square array
        call VCLR(X(LWRK1),1,L1*L1)                                      !clear -LWRK1-
        CALL DGEMM('n','n',L1,LX,L1,1.d0,X(LWRK2),L1,X(LVA),L1,          !transform to MO basis; F.AO * C
     *             0.d0,X(LWRK1),L1)
        CALL DGEMM('t','n',LX,LX,L1,1.d0,X(LVA),L1,X(LWRK1),L1,          !F.MO = C^{\dagger} * (F.AO * C)
     *             0.d0,X(LWRK2),L1)
      CALL DAWRIT(IDAF,IODA,X(LWRK2),LX*LX,707,0)                        !save the c-s FM to record -707-, for further use
c       if(MASWRK)then
c         write(iw,*) 'test F.MO'
c         CALL PRSQ(X(LWRK2),LX,LX,LX)                                   !test print
c       endif
        endif                                                            !if(EKTEA)then
C
                                                                         !transform the relaxed DM to MO basis: P.MO = (S*C)^† * P.AO * S*C
      CALL DAREAD(IDAF,IODA,X(LSMAT),L2,12,0)                            !read in the overlap matrix
      CALL EXPND(X(LSMAT),X(LWRK1),L1,0)                                 !expand in -LWRK1-
        call VCLR(X(LWRK2),1,L1*L1)                                      !clear -LWRK2-
        call dgemm('n','n',L1,LX,L1,1.d0,X(LWRK1),L1,X(LVA),L1,0.d0,     !take S*C in -LWRK2-; L1.LX
     *             X(LWRK2),L1)
        call dcopy(L1*LX,X(LWRK2),1,X(LWRK1),1)                          !copy -LWRK2- to -LWRK1-
        call VCLR(X(LWRK2),1,L1*L1)                                      !clear -LWRK2-
        call dgemm('n','n',L1,LX,L1,1.0d0,X(LFAM),L1,X(LWRK1),L1,0.d0,   !take P(AO) * (SC) in -LWRK2-
     *             X(LWRK2),L1)
        call VCLR(X(LFAM),1,LX*LX)
        call dgemm('t','n',LX,LX,L1,0.5d0,X(LWRK1),L1,X(LWRK2),L1,0.d0,  !take (SC)^\dagger * -LWRK2- in -LFAM-; scale by 0.5 to make one-electron matrix
     *             X(LFAM),LX)
C
        if(EKTEA)then
                call VCLR(X(LWRK2),1,LX*LX)
                do kk=1,LX
                X(LWRK2+(LX+1)*(kk-1))=1.d0
                enddo
        call dscal(LX*LX,-1.d0,X(LFAM),1)                                !scale DM by -1.0
        call VADD(X(LWRK2),1,X(LFAM),1,X(LFAM),1,LX*LX)
        endif                                                            !if(EKTEA)then
C
c       if(MASWRK)then
c         write(iw,*) 'test P.MO'
c         CALL PRSQ(X(LFAM),LX,LX,LX)                                    !test print
c       endif
C
        call dscal(LX*LX,-1.d0,X(LFAM),1)                                !make it negative; then the occupied orbitals will be the first
        info = 0
        call dsyev('V','L',LX,X(LFAM),LX,X(LEA),X(LWRK2),1+6*LX+3*L3,
     *             info)                                                 !diagonalize
        if (Info.ne.0) then
          write(iw,'(1x,"Error in DSYEV in REKS; INFO =",i10)')info
          call ABRT
        endif
C
c       if(maswrk)then                                                   !test print
c               write(iw,*)'---> P.diag:'
c       write(iw,'(5(1x,F14.8))')(X(LEA-1+KK),KK=1,LX)
c       endif
C
        nneg = 0                                                         !find out the negative eigenvalues; these are the occupied orbitals
        thresh = 1.d-6                                                   !threshold for absolute value of zero eigenvalues
        eltot = 0.d0
        elneg = 0.d0
        do i = 1,LX
         eltot = eltot + X(LEA+i-1)
          if(X(LEA+i-1).le.-thresh)then
           elneg = elneg + X(LEA+i-1)
            nneg = nneg + 1
          endif
        enddo
        call VCLR(X(LOMG),1,LX)
        call dcopy(nneg,X(LEA),1,X(LOMG),1)                              !save NO occupations to -LOMG-
        corrfctr = 1.d0
        if(nneg.gt.0)corrfctr = sqrt(eltot/elneg)                        !to correct the orbital occupation numbers by the total number of electrons
        rexprefix = ' SA:'
        if(rexType.ne.0) rexprefix = 'SSR:'                              !set prefix for SA/SSR printing
        fona=-2.d0*X(LOMG+NOCR-1)
        fonb=-2.d0*X(LOMG+NOCS-1)
                if(EKTEA)then                                            !print NO occupations for the "hole" states; EKT for EA's
        fona=2.d0*(1.d0+X(LOMG+LX-NOCR))
        fonb=2.d0*(1.d0+X(LOMG+LX-NOCS))
                endif
        if(MASWRK)then
        write(iw,'(/1x,a,1x,a,
     *       ":  FON(",i3,") =",F9.6," FON(",i3,") =",F9.6)')
     *       rexprefix,'relaxed occupation numbers',NOCR,fona,NOCS,fonb
        if(EKTEA)then
        write(iw,'(/1x,
     *  "Extended Koopmans'' Theorem for Electron Affinities")')
        else
        write(iw,'(/1x,
     *  "Extended Koopmans'' Theorem for Ionization Energies")')
        endif
       write(iw,'(/1x,"==> EKT: there are",i4," Dyson''s orbitals with",
     *             " non-zero strengths")')nneg
        endif
        call VCLR(X(LWRK3),1,nneg*nneg)                                  !prepare to calculate inverse square root of density matrix
        do i = 1,nneg
          pq = (i - 1)*nneg + i - 1
        X(LWRK3+pq) = sqrt(-1.d0/X(LEA+i-1))                             !put 1/Sqrt[eig_k] on the main diagonal of a matrix
        enddo
        call dgemm('n','n',LX,nneg,nneg,1.0d0,X(LFAM),LX,X(LWRK3),nneg,  !C*(eig)^-1/2
     *             0.d0,X(LWRK2),LX)
        call VCLR(X(LWRK1),1,LX*LX)
        call dgemm('n','t',LX,LX,nneg,1.0d0,X(LWRK2),LX,X(LFAM),LX,      !C*(eig)^-1/2 * C^†  in -LWRK1-
     *             0.d0,X(LWRK1),LX)
c       call DAWRIT(IDAF,IODA,X(LWRK1),LX*LX,706,0)                      !save P^-1/2 matrix to 706 for further use
C
        call VCLR(X(LWRK3),1,nneg*nneg)                                  !prepare to restore the density matrix using non-zero eigenvalues only
        do i = 1,nneg
          pq = (i - 1)*nneg + i - 1
        X(LWRK3+pq) = -X(LEA+i-1)                                        !copy eig_k on the main diagonal of a matrix
        enddo
        call dgemm('n','n',LX,nneg,nneg,1.0d0,X(LFAM),LX,X(LWRK3),nneg,  !C*(eig)
     *             0.d0,X(LWRK2),LX)
        call VCLR(X(LWRK3),1,LX*LX)
        call dgemm('n','t',LX,LX,nneg,1.0d0,X(LWRK2),LX,X(LFAM),LX,      !C*(eig) * C^†  in -LWRK3-
     *             0.d0,X(LWRK3),LX)
        call VCLR(X(LOMATA),1,LX*LX)
        call dcopy(LX*LX,X(LWRK3),1,X(LWRK2),1)                          !copy -LWRK3- to -LWRK2-
        call dgemm('n','n',LX,LX,LX,1.0d0,X(LWRK2),LX,X(LWRK3),LX,       !take square of the density matrix in -LOMATA-
     *             0.d0,X(LOMATA),LX)                                    !the square is needed to compute pole strength
C
c       if(MASWRK)then
c         write(iw,*) 'test P.sqrd.new'
c         CALL PRSQ(X(LOMATA),LX,LX,LX)
c       endif
C
        call dcopy(LX*LX,X(LWRK1),1,X(LFAM),1)                           !copy P^-1/2 to -LFAM-
C
c       if(MASWRK)then
c         write(iw,*) 'test P.MO.new'
c         CALL PRSQ(X(LWRK1),LX,LX,LX)
c       endif
C
        goto 702                                                         !bypass the following block; here the Lagrangian with relaxation terms is used
C====================================== here, SA/SSR Lagrangian without relaxation terms is computed ==================
c       call VCLR(X(LWRK1),1,LX*LX)                                      !clear -LWRK1- for Lagrangian matrix in MO representation LX*LX
c         DO L=1,NMIC                                                    !loop over microstates
c            LFAMTMP=LFAMMO + (L-1)*L3                                   !pointer to L-th microstate Fα
c            LFBMTMP=LFBMMO + (L-1)*L3                                   !pointer to L-th microstate Fβ
c       nraL = nmicro(1,1,L) - 1                                         !occupation of -r- orbital in L-th microstate; α-spin
c       nrbL = nmicro(1,2,L) - 1                                         !occupation of -r- orbital in L-th microstate; β-spin
c       nsaL = nmicro(2,1,L) - 1                                         !occupation of -s- orbital in L-th microstate; α-spin
c       nsbL = nmicro(2,2,L) - 1                                         !occupation of -s- orbital in L-th microstate; β-spin
c       if(rexTarget.eq.2)then
c       CLX = X(LCMOSS+L-1)                                              !weighting factor of L-th microstate in OSS state
c       elseif(rexTarget.eq.1)then
c       CLX = X(LCMPPS+L-1)                                              !weighting factor of L-th microstate in PPS state
c       else
c       write(iw,'(1x,"==> Unknown target state,",i4)')rexTarget
c       call ABRT
c       endif
c       if(L.gt.2)CLX = CLX + CLX                                        !double the weights for BS and T microstates
c       ij = 0
c       do j = 1,LX                                                      !loop over orbitals |ket>
c        occja = 0.d0
c        if(j.lt.NOCR)occja = 1.d0
c        if(j.eq.NOCR)occja = dble(nraL)
c        if(j.eq.NOCS)occja = dble(nsaL)
c        occjb = 0.d0
c        if(j.lt.NOCR)occjb = 1.d0
c        if(j.eq.NOCR)occjb = dble(nrbL)
c        if(j.eq.NOCS)occjb = dble(nsbL)
c         do i = 1,LX                                                    !loop over orbitals <bra|
c           occia = 0.d0
c           if(i.lt.NOCR)occia = 1.d0
c           if(i.eq.NOCR)occia = dble(nraL)
c           if(i.eq.NOCS)occia = dble(nsaL)
c           occib = 0.d0
c           if(i.lt.NOCR)occib = 1.d0
c           if(i.eq.NOCR)occib = dble(nrbL)
c           if(i.eq.NOCS)occib = dble(nsbL)
c           pq = (j - 1)*LX + i - 1                                      !indexing Lagrangian matrices; Largangians of PPS and OSS are symmetrized
c           X(LWRK1+pq) = X(LWRK1+pq)                                    !F_{pq} = 1/2*(qε_{pq} + pε_{pq})
c    *    + 0.25d0*CLX * (occia + occja) * X(LFAMTMP + pq)
c    *    + 0.25d0*CLX * (occib + occjb) * X(LFBMTMP + pq)
c         enddo
c       enddo
c         ENDDO                                                          !DO L=1,NMIC
C
c       if(rexType.gt.0)then
c        CALL DAREAD(IDAF,IODA,X(LWRK2) ,L3,704,0)                       !read SA Lagrangian in -LWRK2-; it's multiplied by 1/2
c        cft1 = s2sare(1,1)**2 - s2sare(2,1)**2                          !L.SA = L.PPS + L.OSS; L.SSR0 = (a_11^2 - a_21^2)*L.PPS + a_21^2*L.SA
c        if(rexTarget.eq.2)cft1 = s2sare(2,2)**2 - s2sare(1,2)**2        !L.SSR1 = (a_22^2 - a_12^2)*L.OSS + a_12^2*L.SA
c        cft2 = s2sare(2,1)**2
c        if(rexType.eq.2)cft2 = s2sare(1,2)**2
c        call dscal(LX*LX,cft1,X(LWRK1),1)                               !scale the calculated L by cft1
c        call daxpy(LX*LX,cft2,X(LWRK2),1,X(LWRK1),1)                    !and add L.SA scaled by cft2
c        call dscal(LX*LX,2.d0,X(LWRK1),1)                               !scale the calculated L by 2
C
c       call VCLR(X(LWRK2),1,LX*LX)
c       call DAREAD(IDAF,IODA,X(LWRK2),LX*LX,705,0)                      !read ΔQ1 matrix from 705
c     call SSRDQ2E(X(LWRK2),X(LCMSA),X(LFAMMO),X(LFBMMO),
c    *            X(LBSMA),X(LBSMB),X(LVA),X(LBUF),X(LIBUF),
c    *            X(LSCR1),L1,L2,L3,LX,L7,NOCR,NOCS,NCORE,NMIC)
c       fctr = -2.d0*s2sare(1,1)*s2sare(2,1)                             !coefficient before ΔQ contribution to Lagrangian
c       if(rexTarget.eq.2)fctr = -2.d0*s2sare(1,2)*s2sare(2,2)           !if state 1 is requested
c        call dscal(LX*LX,fctr,X(LWRK2),1)                               !scale ΔQ
c        call daxpy(LX*LX,1.d0,X(LWRK2),1,X(LWRK1),1)                    !and add L.ΔQ to the Lagrangian
c         write(iw,*) 'test ΔQ2'
c         CALL PRSQ(X(LWRK2),LX,LX,LX)
c        call dscal(LX*LX,0.5d0,X(LWRK1),1)                              !scale L back by 1/2
c       endif                                                            !if(rexTyle.gt.0)then
C
c         write(iw,*) 'test L.MO'
c         CALL PRSQ(X(LWRK1),LX,LX,LX)
C======================================================================================================================
C
 702     CALL DAREAD(IDAF,IODA,X(LWRK1) ,L3,704,0)                       !read SA/SSR Lagrangian in MO in -LWRK1-
         call dscal(LX*LX,0.5d0,X(LWRK1),1)                              !scale L by 1/2
C
        if(EKTEA)then
         CALL DAREAD(IDAF,IODA,X(LWRK2),LX*LX,707,0)                     !read the c-s FM from record -707-
         call dscal(LX*LX,-1.d0,X(LWRK1),1)                              !scale L by -1.0
         call VADD(X(LWRK2),1,X(LWRK1),1,X(LWRK1),1,LX*LX)
        endif                                                            !if(EKTEA)then
C
c       if(MASWRK)then
c         write(iw,*) 'test L.MO'
c         CALL PRSQ(X(LWRK1),LX,LX,LX)                                   !test print
c       endif
C
        call VCLR(X(LWRK2),1,LX*LX)                                      !multiply by P^-1/2 on both sides
        call dgemm('n','n',LX,LX,LX,1.0d0,X(LWRK1),LX,X(LFAM),LX,        !on the right
     *             0.d0,X(LWRK2),LX)
        call VCLR(X(LWRK1),1,LX*LX)
        call dgemm('t','n',LX,LX,LX,1.0d0,X(LFAM),LX,X(LWRK2),LX,        !and on the left
     *             0.d0,X(LWRK1),LX)
        info = 0
        call dsyev('V','L',LX,X(LWRK1),LX,X(LEA),X(LWRK2),1+6*L1+3*L3
     *             ,info)
        if (Info.ne.0) then
          write(iw,'(1x,"Error in DSYEV in REKS; INFO =",i10)')info
          call ABRT
        endif
        call dgemm('n','n',LX,LX,LX,1.0d0,X(LFAM),LX,X(LWRK1),LX,        !multiply eigenvectors by P^-1/2
     *             0.d0,X(LWRK2),LX)                                     !now they are normalized on P; X^† * P * X = 1
C
        call dgemm('n','n',LX,LX,LX,1.0d0,X(LOMATA),LX,X(LWRK2),LX,      !compute pole strength X^† * P**2 * X; right product
     *             0.d0,X(LWRK3),LX)
        call dgemm('t','n',LX,LX,LX,1.0d0,X(LWRK2),LX,X(LWRK3),LX,       !left product; put to -LOMATA-
     *             0.d0,X(LOMATA),LX)
C
        call dcopy(LX*LX,X(LWRK1),1,X(LWRK2),1)                          !copy -LWRK1- to -LWRK2-; X.MO orthonormalized
        call dgemm('n','n',L1,LX,LX,1.0d0,X(LVA),L1,X(LWRK2),LX,         !multiply AO eigenvectors and MO eigenvectors C.AO * X.MO
     *             0.d0,X(LWRK1),L1)                                     !now they are in the AO representation; not normalized on S!!!
C
        call VCLR(X(LWRK2),1,L1*L1)
        call VCLR(X(LWRK3),1,LX)
        call VCLR(ijcopy,1,LX)
        indx = 0
        do j = 1,LX                                                      !loop over MOs
         if(abs(X(LEA+j-1)).gt.thresh)then                               !copy only eigenvectors for which eigenvalue is non-zero
          indx = indx + 1
          ijcopy(indx) = j                                               !prepare indexing of pole strengths; to be copied below
          X(LWRK3+indx-1) = X(LEA+j-1)
          do i = 1,LX                                                    !loop over AOs
            pq = (j - 1)*LX + i - 1                                      !pointer to elements of the original orbital set
            ij = (indx - 1)*LX + i - 1                                   !pointer to elements of new orbital set
            X(LWRK2+ij) = X(LWRK1+pq)                                    !copy from old to new
          enddo
         endif
        enddo
        call VCLR(X(LWRK1),1,L3)                                         !determine symmetry labels of EKT orbitals
        call DAREAD(IDAF,IODA,X(LWRK1),L3,45,0)                          !read Q matrix to -LWRK1-
        call DAREAD(IDAF,IODA,X(LFBM),L2,12,0)                           !read overlap to -LFBM-
        call SYMMOS(X(LEA),X(LWRK1),X(LFBM),X(LWRK2),X(LSCR),LX,L1,indx, !determine symmetries; labels go to -LEA-
     *              L1)
C
c       do i = 1,indx
c         pq = (i - 1)*LX + i - 1
c         X(LFBM+i-1) = X(LOMATA + pq)                                   !copy pole strengths to -LFBM-; diag of 2D -> 1D array
c       enddo
        kk=0
        do i = 1,LX                                                      !copy pole strengths to -LFBM-; diag of 2D -> 1D array
        if(ijcopy(i).ne.0)then                                           !only copy the strengths for non-sero eigenvalues of the EKT equation
        kk=kk+1
        pq = (ijcopy(kk) - 1)*LX + ijcopy(kk) - 1
        X(LFBM+kk-1) = X(LOMATA + pq)
        endif
        enddo
C
      if(MASWRK)then                                                     !print teh Dyson's orbitals
         write(iw,'(/10x,46(1H-)/12x,
     *       "EKT orbitals, energies, and pole strengths"
     *       /10x,46(1H-))')
         call PREKT(X(LWRK2),X(LWRK3),X(LFBM),X(LEA),nneg,L1,L1,1)
      endif
C
      END IF                                                             !if(rexEKT)
      END IF                                                             !IF(EXETYP.NE.CHECK)
C
      IF(MASWRK) WRITE(IW,'(/1x,"==> CP-REKS done")')
      if(allocated(ijcopy)) deallocate(ijcopy)                           !deallocate ijcopy, if it was allocated
      CALL RETFM(NEED)
      CALL TIMIT(1)
C
      RETURN
C
      END SUBROUTINE CPREKS
C*MODULE REKS    *DECK PREKT
C>
C> @author     Michael Filatov
C>
C> @date       2021, Oct
C>
      SUBROUTINE PREKT(V,E,P,LABMO,NMO,NAO,LDV,ISTMO)                    !hacked PREVS; to print EKT orbitals with symmetry labels and norms
C
      USE mx_limits, only: mxatm, mxao
      USE prec, ONLY: dp
      IMPLICIT NONE
C
      REAL(KIND=dp), DIMENSION(MXATM) :: ANAM, BNAM
      REAL(KIND=dp), DIMENSION(MXAO) :: BFLAB
      LOGICAL :: DSKWRK, GOPARR, MASWRK
      INTEGER :: IBTYP, ICUT, IDAF, IJK, IJKT, IP, IPTIM, IR, ITOL, IW, &
     &           MASTER, ME, NAV, NOPK, NORMF, NORMP, NPRINT, NPROC
      INTEGER, DIMENSION(950) :: IODA
      REAL(KIND=dp), DIMENSION(10) :: TITLE
      COMMON /IOFILE/ IR, IW, IP, IJK, IJKT, IDAF, NAV, IODA
      COMMON /OUTPUT/ NPRINT, ITOL, ICUT, NORMF, NORMP, NOPK
      COMMON /PAR   / ME, MASTER, NPROC, IBTYP, IPTIM, GOPARR, DSKWRK,  &
     &                MASWRK
      COMMON /RUNLAB/ TITLE, ANAM, BNAM, BFLAB
C
      INTEGER :: ISTMO, LDV, NAO, NMO
      REAL(KIND=dp), DIMENSION(NMO) :: E, P
      INTEGER, DIMENSION(NMO) :: LABMO
      REAL(KIND=dp), DIMENSION(LDV,NMO) :: V
C
      INTEGER :: I, IMAX, IMIN, J, MAX
C
C     ----- PRINT OUT EIGENDATA, WITH MO SYMMETRY LABELS -----
C     THE ROWS ARE LABELED WITH THE BASIS FUNCTION NAMES.
C
      MAX = 5
      IF (NPRINT .EQ. 6) MAX = 10
      IMAX = ISTMO-1
C
      IF (MASWRK) THEN
C
  100 IMIN = IMAX+1
      IMAX = IMAX+MAX
      IF (IMAX .GT. NMO) IMAX = NMO
      WRITE (IW,'(1X)')
      WRITE (IW,'(15X,10(4X,I4,3X))')       (I,       I=IMIN,IMAX)
      WRITE (IW,'(5X,"ENERGY",4X,10F11.6)') (E(I),    I=IMIN,IMAX)
      WRITE (IW,'(5X,"STRENGTH",2X,10F11.6)') (P(I),    I=IMIN,IMAX)
      WRITE (IW,'(16X,10(5X,A4,2X))')       (LABMO(I),I=IMIN,IMAX)
      DO J = 1,NAO
      WRITE (IW,'(I5,2X,A8,10F11.6)') J,BFLAB(J),(V(J,I),I = IMIN,IMAX)
      ENDDO
      IF (IMAX .LT. NMO) GO TO 100
      END IF
      RETURN
      END SUBROUTINE PREKT
C*MODULE REKS    *DECK SSRX2E
C>
C> @author     Michael Filatov
C>
C> @date       2021, Oct
C>
      SUBROUTINE SSRX2E(X2E,CM,WRK1,WRK2,OMATA,OMATB,
     &                  BNMA,BNMB,VA,BUF,IBUF,SCR,
     &                  GRD,WGT,DCH,RHOI,TAUI,AOMAX,GMO,                 !arrays for DFT XC kernels
     &                  FXC,TRAI,COEF,EX,EC,EX0,EC0,IAO,VPRGA,VPRGB,     !arrays for DFT XC kernels
     &                  L1,L2,L3,LX,L7,NOCR,NOCS,NCORE,
     &                  NMIC,SWTC12)
c X2E                   output: - sqrt(n_r) (rs | ft^{r(p-q)}_{Hxc} | pq)
c                               + sqrt(n_s) (rs | ft^{s(p-q)}_{Hxc} | pq)
c if SWTC12=.TRUE.,then
c X2E                   output: - sqrt(n_r) (rs | ft^{r(p-q)}_{Hxc} | pq)
c                               - sqrt(n_s) (rs | ft^{s(p-q)}_{Hxc} | pq)
c CM                    C_L
c WRK1,WRK2             scratch
c OMATA,OMATB           were defined before
c BNMA,BNMB             scratch
c VA                    MOs
c BUF, IBUF             scratch for HF
c L1,L2,L3,LX,L7        the same
c NOCR,NOCS,NCORE       the same
c NMIC                  the number of microstates
c SWTC12                switch between XΔ(.false.) and XΔ12 (.true.)
      USE comm_REKSCM, ONLY: NMICRO, MTTYP, WPPS, WOSS, G1, DNR,
     * DNS, DELTA, FR, FS
      USE mx_limits, only: mxgrid
      USE constants, only: zero, one, two
      USE prec, ONLY: dp
      IMPLICIT NONE
C
      REAL(KIND=dp), DIMENSION(137) :: BSLRD
      REAL(KIND=dp) :: DFTGTHR, DFTTHR, RHOMIN, SW0, SWOFF
      INTEGER :: IGRDTYP, ILENG, INTTYP, JANS, MAXGRD, NDFTFG,          &
     &           NHEX, NINTMX, NPHI, NPHI0, NRAD, NRAD0, NTHE, NTHE0,   &
     &           NTUPL
      INTEGER, DIMENSION(MXGRID) :: NANGPT, NANGPT0
      LOGICAL :: PACK2E, SG1
      COMMON /DFGRID/ DFTTHR, DFTGTHR, SWOFF, SW0, BSLRD, NDFTFG, NRAD, &
     &                NTHE, NPHI, NRAD0, NTHE0, NPHI0, NANGPT, NANGPT0, &
     &                SG1, JANS
      COMMON /INFGRD/ RHOMIN, ILENG, MAXGRD
      COMMON /INTFIL/ NINTMX, NHEX, NTUPL, PACK2E, INTTYP, IGRDTYP
C
      LOGICAL :: SWTC12                                                  !switch to build 2e part of XΔ12, instead of XΔ
      INTEGER :: L1, L2, L3, L7, LX, NCORE, NMIC, NOCR, NOCS
      REAL(KIND=dp), DIMENSION(L1) :: AOMAX
      REAL(KIND=dp), DIMENSION(L1,L1,NMIC) :: BNMA, BNMB
      REAL(KIND=dp), DIMENSION(NINTMX) :: BUF
      REAL(KIND=dp), DIMENSION(NMIC) :: CM
      REAL(KIND=dp), DIMENSION(32*ILENG) :: COEF
      REAL(KIND=dp), DIMENSION(4*L1*ILENG) :: DCH
      REAL(KIND=dp), DIMENSION(35*ILENG) :: EC
      REAL(KIND=dp), DIMENSION(ILENG) :: EC0, EX0
      REAL(KIND=dp), DIMENSION(18*ILENG) :: EX
      REAL(KIND=dp), DIMENSION(L2,2) :: FXC
      REAL(KIND=dp), DIMENSION(L1*8) :: GMO
      REAL(KIND=dp), DIMENSION(MAXGRD*3) :: GRD
      INTEGER, DIMENSION(*) :: IAO
      INTEGER, DIMENSION(NINTMX) :: IBUF
      REAL(KIND=dp), DIMENSION(LX,LX,NMIC) :: OMATA, OMATB
      REAL(KIND=dp), DIMENSION(8*MAXGRD) :: RHOI
      REAL(KIND=dp), DIMENSION(L3) :: SCR
      REAL(KIND=dp), DIMENSION(2*MAXGRD) :: TAUI
      REAL(KIND=dp), DIMENSION(10*ILENG) :: TRAI
      REAL(KIND=dp), DIMENSION(L1,LX) :: VA
      REAL(KIND=dp), DIMENSION(*) :: VPRGA, VPRGB
      REAL(KIND=dp), DIMENSION(MAXGRD) :: WGT
      REAL(KIND=dp), DIMENSION(L3,NMIC) :: WRK1, WRK2
      REAL(KIND=dp), DIMENSION(L7) :: X2E
C
      LOGICAL :: DBGAMAT
      INTEGER :: I, J, L, NPTGRD, NRA, NRB, NSA, NSB
      REAL(KIND=dp) :: RHO, SQNR, SQNS
C
      CALL VCLR(BNMA,1,L3*NMIC)
      CALL VCLR(BNMB,1,L3*NMIC)
      SQNR=sqrt(TWO*DNR)
      SQNS=sqrt(TWO*DNS)

      DO L=1,NMIC
         nra = nmicro(1,1,L) - 1                                         !occupation of -r- orbital in L-th microstate; α-spin
         nrb = nmicro(1,2,L) - 1                                         !occupation of -r- orbital in L-th microstate; β-spin
         nsa = nmicro(2,1,L) - 1                                         !occupation of -s- orbital in L-th microstate; α-spin
         nsb = nmicro(2,2,L) - 1                                         !occupation of -s- orbital in L-th microstate; β-spin
      if(.not.SWTC12)then                                                !do XΔ
         DO I=1,L1
         DO J=1,L1
            BNMA(I,J,L)=                                                 !do (\sqrt(n_s)*n_{s,L}^α - \sqrt(n_r)*n_{r,L}^α)*C_μr*C_νs
     *                  VA(I,NOCR)*VA(J,NOCS)*
     *                 (SQNS*dble(NSA)-SQNR*dble(NRA))
            BNMB(I,J,L)=                                                 !the same for β-spins
     *                  VA(I,NOCR)*VA(J,NOCS)*
     *                 (SQNS*dble(NSB)-SQNR*dble(NRB))
         ENDDO
         ENDDO
      else                                                               !do XΔ12
         DO I=1,L1
         DO J=1,L1
            BNMA(I,J,L)=                                                 !do (-\sqrt(n_s)*n_{s,L}^α - \sqrt(n_r)*n_{r,L}^α)*C_μr*C_νs
     *                  VA(I,NOCR)*VA(J,NOCS)*
     *                 (-SQNS*dble(NSA)-SQNR*dble(NRA))
            BNMB(I,J,L)=                                                 !the same for β-spins
     *                  VA(I,NOCR)*VA(J,NOCS)*
     *                 (-SQNS*dble(NSB)-SQNR*dble(NRB))
         ENDDO
         ENDDO
      endif                                                              !if(.not.SWTC12)then
      ENDDO
C
      CALL VCLR(WRK1,1,L3*NMIC)
      CALL VCLR(WRK2,1,L3*NMIC)
      CALL SSR2E(BNMA,BNMB,WRK1,WRK2,BUF,IBUF,L1,NMIC)                   !do HF part in AO; WRK1/2 - α/β-spin
C
      if(NDFTFG.eq.1)then                                                !do the DFT XC part in AO
      do L=1,4                                                           !loop over microstates
         MTTYP = L                                                       !L is communicated via MTTYP to UDENCNST
         NPTGRD = MAXGRD                                                 !grid setting for the current microstate
         CALL UTDDFTSET(GRD,WGT,DCH,VA,VA,RHOI,TAUI,AOMAX,GMO,           !prepares the density of the L-th microstate; L is communicated to UDENCNST via the MTTYP constant
     *                  ILENG,NPTGRD,L1)
         CALL VCLR(FXC,1,L2*2)
         CALL UTDFXCP2(FXC,RHO,GRD,WGT,DCH,BNMA(1,1,L),BNMB(1,1,L),      !computes integrals of the XC kernel and contracts with the BNMA/BNMB matrices; 1 - α-spin, 2 -β-spin
     *                 RHOI,TAUI,TRAI,COEF,EX,EC,EX0,
     *                 EC0,AOMAX,VPRGA,VPRGB,IAO,ILENG,NPTGRD,
     *                 L1,L2,2,.FALSE.,0)
         call EXPND(FXC(1,1),SCR,L1,0)
         call VADD(WRK1(1,L),1,SCR,1,WRK1(1,L),1,L3)
         call EXPND(FXC(1,2),SCR,L1,0)
         call VADD(WRK2(1,L),1,SCR,1,WRK2(1,L),1,L3)
      enddo                                                              !do L=1,4
      endif                                                              !if(NDFTFG.eq.1)then
C
      DO L=1,NMIC                                                        !convert to MO basis
         CALL DGEMM('N','N',L1,LX,L1,ONE,WRK1(1,L),L1,VA,L1,
     *              ZERO,SCR,L1)
         CALL DGEMM('T','N',LX,LX,L1,ONE,VA,L1,SCR,L1,
     *              ZERO,WRK1(1,L),L1)
         CALL DGEMM('N','N',L1,LX,L1,ONE,WRK2(1,L),L1,VA,L1,
     *              ZERO,SCR,L1)
         CALL DGEMM('T','N',LX,LX,L1,ONE,VA,L1,SCR,L1,
     *              ZERO,WRK2(1,L),L1)
      ENDDO
C
      CALL TSRWSUM(X2E,WRK1,WRK2,L1,OMATA,                               !multiply by occupation number differences and take summation over L
     *             OMATB,CM,NCORE,NOCS,LX,NMIC,L7)
C
      RETURN
      END SUBROUTINE SSRX2E
C*MODULE REKS    *DECK SSR2E
C>
C> @author     Michael Filatov
C>
C> @date       2021, Oct
C>
      SUBROUTINE SSR2E(BSMA,BSMB,ABA,ABB,BUF,IBUF,NBF,NMIC)
      USE mx_limits, only: mxsh, mxgtot, mxao
      USE prec, ONLY: dp
      USE camdft, ONLY: ALPHAC => cam_alpha, BETAC => cam_beta,
     *  CAMMU => cam_mu, CAMFLAG
      USE lrcdft, ONLY: LCFLAG, EMU, EMU2
      IMPLICIT NONE
C
        REAL(KIND=dp), DIMENSION(20)  :: DFTTYP
        REAL(KIND=dp) :: EXENA, EXENB, EXENC
        INTEGER       :: IDFT34, NAUXFUN, NAUXSHL
      COMMON /DFTPAR/ DFTTYP, EXENA, EXENB, EXENC, IDFT34, NAUXFUN,     &
     &                NAUXSHL
        INTEGER, DIMENSION(MXAO) :: IA
      COMMON /IJPAIR/ IA
        INTEGER :: IGRDTYP, INTTYP, NHEX, NINTMX, NTUPL
        LOGICAL :: PACK2E
      COMMON /INTFIL/ NINTMX, NHEX, NTUPL, PACK2E, INTTYP, IGRDTYP
        INTEGER :: IECP, IEFLD, ISCHWZ, NECP
      COMMON /INTOPT/ ISCHWZ, IECP, NECP, IEFLD
        INTEGER :: IDAF, IJK, IPK, IP, IR, IW, NAV
        INTEGER, DIMENSION(950) :: IODA
      COMMON /IOFILE/ IR, IW, IP, IJK, IPK, IDAF, NAV, IODA
        LOGICAL :: LRINT
      COMMON /NLRCF / LRINT
        REAL(KIND=dp), DIMENSION(MXGTOT) :: CD, CF, CG, CH, CI, CP, CS, &
     &                                      EX
        INTEGER, DIMENSION(MXSH) :: KATOM, KLOC, KMAX, KMIN, KNG,       &
     &                              KSTART, KTYPE
        INTEGER :: NSHELL
      COMMON /NSHEL / EX, CS, CP, CD, CF, CG, CH, CI, KSTART, KATOM,    &
     &                KTYPE, KNG, KLOC, KMIN, KMAX, NSHELL
        LOGICAL :: DIRSCF, FDIFF
      COMMON /OPTSCF/ DIRSCF, FDIFF
        INTEGER :: ICUT, ITOL, NOPK, NORMF, NORMP, NPRINT
      COMMON /OUTPUT/ NPRINT, ITOL, ICUT, NORMF, NORMP, NOPK
        LOGICAL :: DSKWRK, GOPARR, MASWRK
        INTEGER :: IBTYP, IPTIM, MASTER, ME, NPROC
      COMMON /PAR   / ME, MASTER, NPROC, IBTYP, IPTIM, GOPARR, DSKWRK,  &
     &                MASWRK
C
      INTEGER :: NBF, NMIC
      REAL(KIND=dp), DIMENSION(NBF,NBF,NMIC) :: ABA, ABB
      REAL(KIND=dp), DIMENSION(NBF*NBF*NMIC) :: BSMA, BSMB
      REAL(KIND=dp), DIMENSION(NINTMX) :: BUF
      INTEGER, DIMENSION(NINTMX) :: IBUF
C
      INTEGER :: NBF3, NINT, NSCHWZ
      REAL(KIND=dp), SAVE :: RHF
      LOGICAL :: SCHWRZ, TDSKWRK
C
C     --- FORM SQUARE NON-SYMMETRIC FOCK-LIKE MATRIX ---
C         DIRECT METHOD = RECOMPUTE 2E- AO INTEGRALS
C         STANDARD METHOD = PROCESS INTEGRALS FROM DISK
C
      NINT   = 0
      NSCHWZ = 0
      NBF3   = NBF*NBF

      IF(DIRSCF) THEN
clsh I have to include direct scf
c         SCHWRZ=ISCHWZ.GT.0
c         NSH2 = (NSHELL*NSHELL+NSHELL)/2
cC
c         CALL BASCHK(LMAX)
c                       NANGM =  4
c         IF(LMAX.EQ.2) NANGM =  6
c         IF(LMAX.EQ.3) NANGM = 10
c         IF(LMAX.EQ.4) NANGM = 15
c         IF(LMAX.EQ.5) NANGM = 21
c         IF(LMAX.EQ.6) NANGM = 28
c         MAXG = NANGM**4
cC
c         CALL VALFM(LOADFM)
c         LGHOND = LOADFM + 1
c         LXINTS = LGHOND + MAXG
c         LDSH   = LXINTS + NSH2
c         LAST   = LDSH     + NSH2
c         NEED   = LAST - LOADFM - 1
c         CALL GETFM(NEED)
cC
c         IF(SCHWRZ) THEN
c            DUMMY = 0.0D+00
c            CALL SHLTD(RHF,BALL,DUMMY,XX(LDSH),NBF,NSH2,1)
c            CALL DAREAD(IDAF,IODA,XX(LXINTS),NSH2,54,0)
cC
c         END IF
c         CALL MRSFTWOEI(SCHWRZ,NINT,NSCHWZ,NBF,XX(LXINTS),NSH2,
c     *       XX(LGHOND),MAXG,IA,BALL,AGDLR,XX(LDSH),1,.TRUE.,1,
c     *                  BO2V,BO1V,BCO1,BCO2,BALL,
c     *                  CO2V,CO1V,CCO1,CCO2,AGDLR,ACO2V,ACO1V,
c     *                  ACCO1,ACCO2,ADO2V,ADO1V,ADCO1,ADCO2)
c         CALL RETFM(NEED)
      ELSE
         TDSKWRK = DSKWRK
         DSKWRK  = .TRUE.

         if(LCFLAG) then
          LRINT=.TRUE.
          CALL SSRADISK(BSMA,BSMB,ABA,ABB,NBF,BUF,IBUF,
     *                 NINTMX,NOPK,NMIC)
          LRINT=.FALSE.
         endif
         IF(CAMFLAG) THEN
           LRINT     = .TRUE.
           EMU       = CAMMU
           EMU2      = CAMMU*CAMMU
           DFTTYP(3) = BETAC
           CALL SSRADISK(BSMA,BSMB,ABA,ABB,NBF,BUF,IBUF,
     *                   NINTMX,NOPK,NMIC)
           DFTTYP(3) = ALPHAC
           LRINT     = .FALSE.
         endif
         CALL SSRADISK(BSMA,BSMB,ABA,ABB,NBF,BUF,IBUF,
     *                 NINTMX,NOPK,NMIC)
         DSKWRK  = TDSKWRK
      END IF
C
C     --- SUM UP PARTIAL FOCK-LIKE MATRICES ---
C
      IF(GOPARR) THEN
         CALL DDI_GSUMI(2311,NINT,1)
         CALL DDI_GSUMI(2312,NSCHWZ,1)
         CALL DDI_GSUMF(2313,ABA,NBF3*NMIC)
         CALL DDI_GSUMF(2314,ABB,NBF3*NMIC)
      END IF
C
      RETURN
      END SUBROUTINE SSR2E
C*MODULE REKS    *DECK SSRADISK
C>
C> @author     Michael Filatov
C>
C> @date       2021, Oct
C>
      SUBROUTINE SSRADISK(BSMA,BSMB,ABA,ABB,L1,XX,IX,NINTMX,NOPK,NMIC)
C
      USE camdft, ONLY: CAMFLAG
      USE lrcdft, ONLY: LCFLAG, LRFILE
      USE mx_limits, only: mxao
      USE constants, only: half
      USE prec, ONLY: dp
      IMPLICIT NONE
C
      REAL(KIND=dp), DIMENSION(20) :: DFTTYP
      REAL(KIND=dp) :: EXENA, EXENB, EXENC
      INTEGER, DIMENSION(MXAO) :: IA
      INTEGER :: IDAF, IDFT34, IP, IPK, IR, IS, IW, LABSIZ, NAUXFUN,    &
     &           NAUXSHL, NAV
      INTEGER, DIMENSION(950) :: IODA
      INTEGER, DIMENSION(7) :: NORDER
      COMMON /DFTPAR/ DFTTYP, EXENA, EXENB, EXENC, IDFT34, NAUXFUN,     &
     &                NAUXSHL
      COMMON /IJPAIR/ IA
      COMMON /IOFILE/ IR, IW, IP, IS, IPK, IDAF, NAV, IODA
      LOGICAL ::  LRINT
      COMMON /NLRCF / LRINT
      COMMON /ORDOPT/ NORDER
      COMMON /PCKLAB/ LABSIZ
C
      INTEGER :: L1, NINTMX, NMIC, NOPK
      REAL(KIND=dp), DIMENSION(L1,L1,NMIC) :: ABA, ABB, BSMA, BSMB
      INTEGER, DIMENSION(*) :: IX
      REAL(KIND=dp), DIMENSION(NINTMX) :: XX
C
      REAL(KIND=dp) :: CSCALT, DIJ, DIKA, DIKB, DILA, DILB, DJI, DJKA,  &
     &                 DJKB, DJLA, DJLB, DKIA, DKIB, DKJA, DKJB, DKL,   &
     &                 DLIA, DLIB, DLJA, DLJB, DLK, HFSCAL, VAL, VALEX, &
     &                 VALES
      INTEGER :: I, IMIC, IPACK, J, JPACK, K, KPACK, L, LPACK, M,       &
     &           NIJ, NINT, NKL, NPACK, NXX, LABEL
C
      if(.not.LRINT) CALL SEQREW(IS)
      if(     LRINT) CALL SEQREW(LRFILE)
      HFSCAL=DFTTYP(3)
      CSCALT=1.0D+00
      IF(LCFLAG) THEN
         IF(LRINT) THEN
            HFSCAL=1.0D+00
            CSCALT=0.0D+00
         ELSE
            HFSCAL=0.0D+00
            CSCALT=1.0D+00
         ENDIF
      ENDIF
      IF(CAMFLAG.AND.LRINT) CSCALT=0.0D+00
c      write(iw,*) 'hfscal',hfscal
      I = 0
      J = 0
      K = 0
      L = 0
      NXX = 0
      IF(NOPK.NE.1) THEN
         WRITE(IW,*) 'NOPK.NE.1'
         WRITE(IW,*) 'REKS CALCULATION DOES NOT SUPPORT ',
     *               'THIS INTEGRAL TYPE'
         CALL ABRT
      END IF
C
C     ----- INTEGRALS ARE NOT IN SUPERMATRIX FORM (NOPK=.TRUE.) -----
C
  620 if(.not.LRINT) CALL PREAD(IS,    XX,IX,NXX,NINTMX)
      if(     LRINT) CALL PREAD(LRFILE,XX,IX,NXX,NINTMX)
      IF(NXX .NE. 0) THEN
         NINT = IABS(NXX)
         IF(NINT .GT. NINTMX) CALL ABRT
         DO IMIC=1,NMIC                                                  ! microstate
            DO M = 1,NINT
C
               NPACK = M
               IF (LABSIZ .EQ. 2) THEN
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
               I = IPACK
               J = JPACK
               K = KPACK
               L = LPACK
C            
               VAL = XX(M)
               NIJ = IA(I)+J
               NKL = IA(K)+L
C            
C                 USING SQUARE CANONICAL INTEGRAL FILE
C            
               IF(NORDER(7) .EQ. 1) THEN
                  IF(NKL .GT. NIJ) CYCLE
                  IF(I .EQ. J) VAL=VAL*HALF
                  IF(K .EQ. L) VAL=VAL*HALF
                  IF(NIJ .EQ. NKL) VAL=VAL*HALF
               END IF
             
               VALES  =  VAL * CSCALT
               VALEX  =  VAL * HFSCAL * HALF
C            
               DKJA=BSMA(K,J,IMIC)
               DLJA=BSMA(L,J,IMIC)
               DKIA=BSMA(K,I,IMIC)
               DLIA=BSMA(L,I,IMIC)
               DJLA=BSMA(J,L,IMIC)
               DILA=BSMA(I,L,IMIC)
               DJKA=BSMA(J,K,IMIC)
               DIKA=BSMA(I,K,IMIC)
C            
               ABA(I,L,IMIC)=ABA(I,L,IMIC)-VALEX*DKJA! (IJ|LK)
               ABA(I,K,IMIC)=ABA(I,K,IMIC)-VALEX*DLJA! (IJ|LK)
               ABA(J,L,IMIC)=ABA(J,L,IMIC)-VALEX*DKIA! (IJ|LK)
               ABA(J,K,IMIC)=ABA(J,K,IMIC)-VALEX*DLIA! (IJ|LK)
               ABA(K,I,IMIC)=ABA(K,I,IMIC)-VALEX*DJLA! (IJ|LK)
               ABA(K,J,IMIC)=ABA(K,J,IMIC)-VALEX*DILA! (IJ|LK)
               ABA(L,I,IMIC)=ABA(L,I,IMIC)-VALEX*DJKA! (IJ|LK)
               ABA(L,J,IMIC)=ABA(L,J,IMIC)-VALEX*DIKA! (IJ|LK)
C            
               ABA(I,L,IMIC)=ABA(I,L,IMIC)-VALEX*DJKA! (IJ|LK)
               ABA(I,K,IMIC)=ABA(I,K,IMIC)-VALEX*DJLA! (IJ|LK)
               ABA(J,L,IMIC)=ABA(J,L,IMIC)-VALEX*DIKA! (IJ|LK)
               ABA(J,K,IMIC)=ABA(J,K,IMIC)-VALEX*DILA! (IJ|LK)
               ABA(K,I,IMIC)=ABA(K,I,IMIC)-VALEX*DLJA! (IJ|LK)
               ABA(K,J,IMIC)=ABA(K,J,IMIC)-VALEX*DLIA! (IJ|LK)
               ABA(L,I,IMIC)=ABA(L,I,IMIC)-VALEX*DKJA! (IJ|LK)
               ABA(L,J,IMIC)=ABA(L,J,IMIC)-VALEX*DKIA! (IJ|LK)
C            
               DKJB=BSMB(K,J,IMIC)
               DLJB=BSMB(L,J,IMIC)
               DKIB=BSMB(K,I,IMIC)
               DLIB=BSMB(L,I,IMIC)
               DJLB=BSMB(J,L,IMIC)
               DILB=BSMB(I,L,IMIC)
               DJKB=BSMB(J,K,IMIC)
               DIKB=BSMB(I,K,IMIC)
C            
               ABB(I,L,IMIC)=ABB(I,L,IMIC)-VALEX*DKJB! (IJ|LK)
               ABB(I,K,IMIC)=ABB(I,K,IMIC)-VALEX*DLJB! (IJ|LK)
               ABB(J,L,IMIC)=ABB(J,L,IMIC)-VALEX*DKIB! (IJ|LK)
               ABB(J,K,IMIC)=ABB(J,K,IMIC)-VALEX*DLIB! (IJ|LK)
               ABB(K,I,IMIC)=ABB(K,I,IMIC)-VALEX*DJLB! (IJ|LK)
               ABB(K,J,IMIC)=ABB(K,J,IMIC)-VALEX*DILB! (IJ|LK)
               ABB(L,I,IMIC)=ABB(L,I,IMIC)-VALEX*DJKB! (IJ|LK)
               ABB(L,J,IMIC)=ABB(L,J,IMIC)-VALEX*DIKB! (IJ|LK)
C            
               ABB(I,L,IMIC)=ABB(I,L,IMIC)-VALEX*DJKB! (IJ|LK)
               ABB(I,K,IMIC)=ABB(I,K,IMIC)-VALEX*DJLB! (IJ|LK)
               ABB(J,L,IMIC)=ABB(J,L,IMIC)-VALEX*DIKB! (IJ|LK)
               ABB(J,K,IMIC)=ABB(J,K,IMIC)-VALEX*DILB! (IJ|LK)
               ABB(K,I,IMIC)=ABB(K,I,IMIC)-VALEX*DLJB! (IJ|LK)
               ABB(K,J,IMIC)=ABB(K,J,IMIC)-VALEX*DLIB! (IJ|LK)
               ABB(L,I,IMIC)=ABB(L,I,IMIC)-VALEX*DKJB! (IJ|LK)
               ABB(L,J,IMIC)=ABB(L,J,IMIC)-VALEX*DKIB! (IJ|LK)
C            
               DKL=BSMA(K,L,IMIC)+BSMB(K,L,IMIC)
               DLK=BSMA(L,K,IMIC)+BSMB(L,K,IMIC)
               DIJ=BSMA(I,J,IMIC)+BSMB(I,J,IMIC)
               DJI=BSMA(J,I,IMIC)+BSMB(J,I,IMIC)
C            
               ABA(I,J,IMIC)=ABA(I,J,IMIC)+VALES*(DKL+DLK)! (IJ|LK)
               ABA(J,I,IMIC)=ABA(J,I,IMIC)+VALES*(DKL+DLK)! (IJ|LK)
               ABA(K,L,IMIC)=ABA(K,L,IMIC)+VALES*(DIJ+DJI)! (IJ|LK)
               ABA(L,K,IMIC)=ABA(L,K,IMIC)+VALES*(DIJ+DJI)! (IJ|LK)
C            
               ABB(I,J,IMIC)=ABB(I,J,IMIC)+VALES*(DKL+DLK)! (IJ|LK)
               ABB(J,I,IMIC)=ABB(J,I,IMIC)+VALES*(DKL+DLK)! (IJ|LK)
               ABB(K,L,IMIC)=ABB(K,L,IMIC)+VALES*(DIJ+DJI)! (IJ|LK)
               ABB(L,K,IMIC)=ABB(L,K,IMIC)+VALES*(DIJ+DJI)! (IJ|LK)
C
            END DO
         END DO
C
         IF(NXX .GT. 0) GO TO 620
      END IF
C
C
      if(.not.LRINT) CALL SEQREW(IS)
      if(     LRINT) CALL SEQREW(LRFILE)
C
      RETURN
      END SUBROUTINE SSRADISK
C*MODULE REKS    *DECK SSRDQ2E
C>
C> @author     Michael Filatov
C>
C> @date       2021, Oct
C>
      SUBROUTINE SSRDQ2E(X2E,CM,WRK1,WRK2,
     *                   BNMA,BNMB,VA,BUF,IBUF,SCR,
     *                   GRD,WGT,DCH,RHOI,TAUI,AOMAX,GMO,                !arrays for DFT XC kernels
     *                   FXC,TRAI,COEF,EX,EC,EX0,EC0,IAO,VPRGA,VPRGB,    !arrays for DFT XC kernels
     *                   L1,L2,L3,LX,L7,NOCR,NOCS,NCORE,
     *                   NMIC,SWTC12)
C> X2E                   output: - sqrt(n_r) (rs | ft^{r(p+q)}_{Hxc} | pq)
C>                               + sqrt(n_s) (rs | ft^{s(p+q)}_{Hxc} | pq)
C> if SWTC12=.TRUE.,then
C> X2E                   output: - sqrt(n_r) (rs | ft^{r(p+q)}_{Hxc} | pq)
C>                               - sqrt(n_s) (rs | ft^{s(p+q)}_{Hxc} | pq)
C> CM                    C_L
C> WRK1,WRK2             scratch
C> OMATA,OMATB           were defined before
C> BNMA,BNMB             scratch
C> VA                    MOs
C> BUF, IBUF             scratch for HF
C> L1,L2,L3,LX,L7        the same
C> NOCR,NOCS,NCORE       the same
C> NMIC                  the number of microstates
C> SWTC12                switch between ΔQ2(.false.) and Δ12Q2 (.true.)
C>
      USE comm_REKSCM, ONLY: NMICRO, MTTYP, WPPS, WOSS, G1, DNR,
     * DNS, DELTA, FR, FS
      USE mx_limits, only: mxgrid
      USE constants, only: zero, half, one, two
      USE prec, ONLY: dp
      IMPLICIT NONE
C
      REAL(KIND=dp), DIMENSION(137) :: BSLRD
      REAL(KIND=dp) :: DFTGTHR, DFTTHR, RHOMIN, SW0, SWOFF
      INTEGER :: IGRDTYP, ILENG, INTTYP, JANS, MAXGRD, NDFTFG,          &
     &           NHEX, NINTMX, NPHI, NPHI0, NRAD, NRAD0, NTHE, NTHE0,   &
     &           NTUPL
      INTEGER, DIMENSION(MXGRID) :: NANGPT, NANGPT0
      LOGICAL :: PACK2E, SG1
      COMMON /DFGRID/ DFTTHR, DFTGTHR, SWOFF, SW0, BSLRD, NDFTFG, NRAD, &
     &                NTHE, NPHI, NRAD0, NTHE0, NPHI0, NANGPT, NANGPT0, &
     &                SG1, JANS
      COMMON /INFGRD/ RHOMIN, ILENG, MAXGRD
      COMMON /INTFIL/ NINTMX, NHEX, NTUPL, PACK2E, INTTYP, IGRDTYP
C
      INTEGER :: L1, L2, L3, L7, LX, NCORE, NMIC, NOCR, NOCS
      LOGICAL :: SWTC12                                                  !switch to build 2e part of Δ12Q2, instead of ΔQ2
      REAL(KIND=dp), DIMENSION(L1) :: AOMAX
      REAL(KIND=dp), DIMENSION(L1,L1,NMIC) :: BNMA, BNMB
      REAL(KIND=dp), DIMENSION(NINTMX) :: BUF
      REAL(KIND=dp), DIMENSION(NMIC) :: CM
      REAL(KIND=dp), DIMENSION(32*ILENG) :: COEF
      REAL(KIND=dp), DIMENSION(4*L1*ILENG) :: DCH
      REAL(KIND=dp), DIMENSION(35*ILENG) :: EC
      REAL(KIND=dp), DIMENSION(ILENG) :: EC0, EX0
      REAL(KIND=dp), DIMENSION(18*ILENG) :: EX
      REAL(KIND=dp), DIMENSION(L2,2) :: FXC
      REAL(KIND=dp), DIMENSION(L1*8) :: GMO
      REAL(KIND=dp), DIMENSION(MAXGRD*3) :: GRD
      INTEGER, DIMENSION(*) :: IAO
      INTEGER, DIMENSION(NINTMX) :: IBUF
      REAL(KIND=dp), DIMENSION(8*MAXGRD) :: RHOI
      REAL(KIND=dp), DIMENSION(L3) :: SCR
      REAL(KIND=dp), DIMENSION(2*MAXGRD) :: TAUI
      REAL(KIND=dp), DIMENSION(10*ILENG) :: TRAI
      REAL(KIND=dp), DIMENSION(L1,LX) :: VA
      REAL(KIND=dp), DIMENSION(*) :: VPRGA, VPRGB
      REAL(KIND=dp), DIMENSION(MAXGRD) :: WGT
      REAL(KIND=dp), DIMENSION(L3,NMIC) :: WRK1, WRK2
      REAL(KIND=dp), DIMENSION(LX*LX) :: X2E
C
      LOGICAL :: DBGAMAT
      REAL(KIND=dp) :: OCCIA, OCCIB, OCCJA, OCCJB, RHO, SQNR,           &
     &                SQNS, TMPA, TMPB, WMTMP
      INTEGER :: I, IJ, J, JLX, L, NPTGRD, NRA, NRB, NSA, NSB
C
      CALL VCLR(BNMA,1,L3*NMIC)
      CALL VCLR(BNMB,1,L3*NMIC)
C
      SQNR=sqrt(TWO*DNR)
      SQNS=sqrt(TWO*DNS)

      DO L=1,NMIC
         nra = nmicro(1,1,L) - 1                                         !occupation of -r- orbital in L-th microstate; α-spin
         nrb = nmicro(1,2,L) - 1                                         !occupation of -r- orbital in L-th microstate; β-spin
         nsa = nmicro(2,1,L) - 1                                         !occupation of -s- orbital in L-th microstate; α-spin
         nsb = nmicro(2,2,L) - 1                                         !occupation of -s- orbital in L-th microstate; β-spin

      if(.not.SWTC12)then                                                !do ΔQ2
         DO I=1,L1
         DO J=1,L1
            BNMA(I,J,L)=
     *                  VA(I,NOCR)*VA(J,NOCS)*
     *                 (SQNS*dble(NSA)-SQNR*dble(NRA))
            BNMB(I,J,L)=
     *                  VA(I,NOCR)*VA(J,NOCS)*
     *                 (SQNS*dble(NSB)-SQNR*dble(NRB))
         ENDDO
         ENDDO
      else                                                               !do Δ12Q2
         DO I=1,L1
         DO J=1,L1
            BNMA(I,J,L)=
     *                  VA(I,NOCR)*VA(J,NOCS)*
     *                 (-SQNS*dble(NSA)-SQNR*dble(NRA))
            BNMB(I,J,L)=
     *                  VA(I,NOCR)*VA(J,NOCS)*
     *                 (-SQNS*dble(NSB)-SQNR*dble(NRB))
         ENDDO
         ENDDO
      endif                                                              !if(.not.SWTC12)then
      ENDDO
C
      CALL VCLR(WRK1,1,L3*NMIC)
      CALL VCLR(WRK2,1,L3*NMIC)
      CALL SSR2E(BNMA,BNMB,WRK1,WRK2,BUF,IBUF,L1,NMIC)
C
      if(NDFTFG.eq.1)then                                                !do the DFT XC part in AO
      do L=1,4                                                           !loop over microstates
         MTTYP = L                                                       !L is communicated via MTTYP to UDENCNST
         NPTGRD = MAXGRD                                                 !grid setting for the current microstate
         CALL UTDDFTSET(GRD,WGT,DCH,VA,VA,RHOI,TAUI,AOMAX,GMO,           !prepares the density of the L-th microstate; L is communicated to UDENCNST via the MTTYP constant
     *                  ILENG,NPTGRD,L1)
         CALL VCLR(FXC,1,L2*2)
         CALL UTDFXCP2(FXC,RHO,GRD,WGT,DCH,BNMA(1,1,L),BNMB(1,1,L),      !computes integrals of the XC kernel and contracts with the BNMA/BNMB matrices; 1 - α-spin, 2 -β-spin
     *                 RHOI,TAUI,TRAI,COEF,EX,EC,EX0,
     *                 EC0,AOMAX,VPRGA,VPRGB,IAO,ILENG,NPTGRD,
     *                 L1,L2,2,.FALSE.,0)
         call EXPND(FXC(1,1),SCR,L1,0)
         call VADD(WRK1(1,L),1,SCR,1,WRK1(1,L),1,L3)
         call EXPND(FXC(1,2),SCR,L1,0)
         call VADD(WRK2(1,L),1,SCR,1,WRK2(1,L),1,L3)
      enddo                                                              !do L=1,4
      endif                                                              !if(NDFTFG.eq.1)then
C
      DO L=1,NMIC                                                        !convert to MO
         CALL DGEMM('N','N',L1,LX,L1,ONE,WRK1(1,L),L1,VA,L1,
     *              ZERO,SCR,L1)
         CALL DGEMM('T','N',LX,LX,L1,ONE,VA,L1,SCR,L1,
     *              ZERO,WRK1(1,L),L1)
         CALL DGEMM('N','N',L1,LX,L1,ONE,WRK2(1,L),L1,VA,L1,
     *              ZERO,SCR,L1)
         CALL DGEMM('T','N',LX,LX,L1,ONE,VA,L1,SCR,L1,
     *              ZERO,WRK2(1,L),L1)
      ENDDO
C
      DO L=1,NMIC
         WMTMP = CM(L)
         IF(L.LE.2) WMTMP = WMTMP*HALF

         nra = nmicro(1,1,L) - 1                                         !occupation of -r- orbital in L-th microstate; α-spin
         nrb = nmicro(1,2,L) - 1                                         !occupation of -r- orbital in L-th microstate; β-spin
         nsa = nmicro(2,1,L) - 1                                         !occupation of -s- orbital in L-th microstate; α-spin
         nsb = nmicro(2,2,L) - 1                                         !occupation of -s- orbital in L-th microstate; β-spin

         IJ=0
         DO J=1,LX
         occja = 1.d0
         if(j.eq.NOCR)occja = dble(nra)
         if(j.eq.NOCS)occja = dble(nsa)
         if(j.gt.NOCS)occja = 0.d0
         occjb = 1.d0
         if(j.eq.NOCR)occjb = dble(nrb)
         if(j.eq.NOCS)occjb = dble(nsb)
         if(j.gt.NOCS)occjb = 0.d0
         JLX = (J - 1) * LX
            DO I=1,LX
            occia = 1.d0
            if(i.eq.NOCR)occia = dble(nra)
            if(i.eq.NOCS)occia = dble(nsa)
            if(i.gt.NOCS)occia = 0.d0
            occib = 1.d0
            if(i.eq.NOCR)occib = dble(nrb)
            if(i.eq.NOCS)occib = dble(nsb)
            if(i.gt.NOCS)occib = 0.d0
               IJ = IJ+1
               TMPA = WRK1(JLX+I,L)*(occia+occja)
               TMPB = WRK2(JLX+I,L)*(occib+occjb)
C
               X2E(IJ)=X2E(IJ)+WMTMP*(TMPA+TMPB)
            END DO
         END DO
      END DO
C
      RETURN
      END SUBROUTINE SSRDQ2E
C*MODULE REKS    *DECK SSRQ2E
C>
C> @author     Michael Filatov
C>
C> @date       2021, Oct
C>
      SUBROUTINE SSRQ2E(X2E,ZX,CM,WRK1,WRK2,
     *                  BNMA,BNMB,VA,BUF,IBUF,SCR,
     *                  GRD,WGT,DCH,RHOI,TAUI,AOMAX,GMO,                 !arrays for DFT XC kernels
     *                  FXC,TRAI,COEF,EX,EC,EX0,EC0,IAO,VPRGA,VPRGB,     !arrays for DFT XC kernels
     *                  L1,L2,L3,LX,L7,NOCR,NOCS,NCORE,
     *                  NMIC)
      USE comm_REKSCM, ONLY: NMICRO, MTTYP, WPPS, WOSS, G1, DNR,
     * DNS, DELTA, FR, FS
C>
C> X2E                  output: (1/2) \sum_{p<q} Z_{pq} (pq | ft^{(p-q)(i+j)}_{Hxc} | ij)
C> ZX                   input: Z-vector for state X
C> CM                   C_L
C> WRK1,WRK2            scratch
C> BNMA,BNMB            scratch
C> VA                   MOs
C> BUF, IBUF            scratch for DFT
C> L1,L2,L3,LX,L7       the same
C> NOCR,NOCS,NCORE      the same
C> NMIC                 the number of microstates
C>
      USE mx_limits, only: mxgrid
      USE constants, only: zero, half, one
      USE prec, ONLY: dp
      IMPLICIT NONE
C
      REAL(KIND=dp), DIMENSION(137) :: BSLRD
      REAL(KIND=dp) :: DFTGTHR, DFTTHR, RHOMIN, SW0, SWOFF
      INTEGER :: IGRDTYP, ILENG, INTTYP, JANS, MAXGRD, NDFTFG,          &
     &           NHEX, NINTMX, NPHI, NPHI0, NRAD, NRAD0, NTHE, NTHE0,   &
     &           NTUPL
      INTEGER, DIMENSION(MXGRID) :: NANGPT, NANGPT0
      LOGICAL :: PACK2E, SG1
      COMMON /DFGRID/ DFTTHR, DFTGTHR, SWOFF, SW0, BSLRD, NDFTFG, NRAD, &
     &                NTHE, NPHI, NRAD0, NTHE0, NPHI0, NANGPT, NANGPT0, &
     &                SG1, JANS
      COMMON /INFGRD/ RHOMIN, ILENG, MAXGRD
      COMMON /INTFIL/ NINTMX, NHEX, NTUPL, PACK2E, INTTYP, IGRDTYP
C
      INTEGER :: L1, L2, L3, L7, LX, NCORE, NMIC, NOCR, NOCS
      REAL(KIND=dp), DIMENSION(L1) :: AOMAX
      REAL(KIND=dp), DIMENSION(L1,L1,NMIC) :: BNMA, BNMB
      REAL(KIND=dp), DIMENSION(NINTMX) :: BUF
      REAL(KIND=dp), DIMENSION(NMIC) :: CM
      REAL(KIND=dp), DIMENSION(32*ILENG) :: COEF
      REAL(KIND=dp), DIMENSION(4*L1*ILENG) :: DCH
      REAL(KIND=dp), DIMENSION(35*ILENG) :: EC
      REAL(KIND=dp), DIMENSION(ILENG) :: EC0, EX0
      REAL(KIND=dp), DIMENSION(18*ILENG) :: EX
      REAL(KIND=dp), DIMENSION(L2,2) :: FXC
      REAL(KIND=dp), DIMENSION(L1*8) :: GMO
      REAL(KIND=dp), DIMENSION(MAXGRD*3) :: GRD
      INTEGER, DIMENSION(*) :: IAO
      INTEGER, DIMENSION(NINTMX) :: IBUF
      REAL(KIND=dp), DIMENSION(8*MAXGRD) :: RHOI
      REAL(KIND=dp), DIMENSION(L3) :: SCR
      REAL(KIND=dp), DIMENSION(2*MAXGRD) :: TAUI
      REAL(KIND=dp), DIMENSION(10*ILENG) :: TRAI
      REAL(KIND=dp), DIMENSION(L1,LX) :: VA
      REAL(KIND=dp), DIMENSION(*) :: VPRGA, VPRGB
      REAL(KIND=dp), DIMENSION(MAXGRD) :: WGT
      REAL(KIND=dp), DIMENSION(L3,NMIC) :: WRK1, WRK2
      REAL(KIND=dp), DIMENSION(LX*LX) :: X2E
      REAL(KIND=dp), DIMENSION(L7) :: ZX
C
      LOGICAL :: DBGAMAT
      REAL(KIND=dp) :: OCCIA, OCCIB, OCCJA, OCCJB, RHO, TMPA,           &
     &                 TMPB, WMTMP
      INTEGER :: I, IJ, J, JLX, L, NPTGRD, NRA, NRB, NSA, NSB
C
      CALL VCLR(BNMA,1,L3*NMIC)
      CALL VCLR(BNMB,1,L3*NMIC)
C
      DO L=1,NMIC
         nra = nmicro(1,1,L) - 1                                         !occupation of -r- orbital in L-th microstate; α-spin
         nrb = nmicro(1,2,L) - 1                                         !occupation of -r- orbital in L-th microstate; β-spin
         nsa = nmicro(2,1,L) - 1                                         !occupation of -s- orbital in L-th microstate; α-spin
         nsb = nmicro(2,2,L) - 1                                         !occupation of -s- orbital in L-th microstate; β-spin

        call VCLR(SCR,1,L3)                                              !clear SCR for occ. number weighted Z = Z * (n_{p,L}^σ - n_{q,L}^σ)
        ij = 0
        do j = NOCR,LX                                                   !make occ. number weighted Z
         occja = 0.d0
         if(j.eq.NOCR)occja = dble(nra)
         if(j.eq.NOCS)occja = dble(nsa)
          JLX = (J - 1) * LX
          do i = 1,NOCS
            occia = 1.d0
            if(i.eq.NOCR)occia = dble(nra)
            if(i.eq.NOCS)occia = dble(nsa)
            ij = ij + 1
            SCR(JLX+I) = ZX(ij)*(occia - occja)                          ! Z_{pq} * (n_{p,L}^α - n_{q,L}^α); p<q
          enddo
        enddo
        CALL VCLR(WRK1,1,L3*NMIC)                                        !transform to AO basis; use WRK1 for scratch
        call dgemm('n','t',LX,L1,LX,1.d0,SCR,LX,VA,L1,
     *             0.d0,WRK1,LX)
        call dgemm('n','n',L1,L1,LX,1.d0,VA,L1,WRK1,LX,
     *             0.d0,BNMA(1,1,L),L1)
C
        call VCLR(SCR,1,L3)                                              !clear SCR for occ. number weighted Z = Z * (n_{p,L}^σ - n_{q,L}^σ)
        ij = 0
        do j = NOCR,LX                                                   !make occ. number weighted Z
         occjb = 0.d0
         if(j.eq.NOCR)occjb = dble(nrb)
         if(j.eq.NOCS)occjb = dble(nsb)
          JLX = (J - 1) * LX
          do i = 1,NOCS
            occib = 1.d0
            if(i.eq.NOCR)occib = dble(nrb)
            if(i.eq.NOCS)occib = dble(nsb)
            ij = ij + 1
            SCR(JLX+I) = ZX(ij)*(occib - occjb)                          ! Z_{pq} * (n_{p,L}^β - n_{q,L}^β); p<q
          enddo
        enddo
        CALL VCLR(WRK1,1,L3*NMIC)                                        !transform to AO basis; use WRK1 for scratch
        call dgemm('n','t',LX,L1,LX,1.d0,SCR,LX,VA,L1,
     *             0.d0,WRK1,LX)
        call dgemm('n','n',L1,L1,LX,1.d0,VA,L1,WRK1,LX,
     *             0.d0,BNMB(1,1,L),L1)
      ENDDO
C
      CALL VCLR(WRK1,1,L3*NMIC)
      CALL VCLR(WRK2,1,L3*NMIC)
      CALL SSR2E(BNMA,BNMB,WRK1,WRK2,BUF,IBUF,L1,NMIC)
C
      if(NDFTFG.eq.1)then                                                !do the DFT XC part in AO
      do L=1,4                                                           !loop over microstates
         MTTYP = L                                                       !L is communicated via MTTYP to UDENCNST
         NPTGRD = MAXGRD                                                 !grid setting for the current microstate
         CALL UTDDFTSET(GRD,WGT,DCH,VA,VA,RHOI,TAUI,AOMAX,GMO,           !prepares the density of the L-th microstate; L is communicated to UDENCNST via the MTTYP constant
     *                  ILENG,NPTGRD,L1)
         CALL VCLR(FXC,1,L2*2)
         CALL UTDFXCP2(FXC,RHO,GRD,WGT,DCH,BNMA(1,1,L),BNMB(1,1,L),      !computes integrals of the XC kernel and contracts with the BNMA/BNMB matrices; 1 - α-spin, 2 -β-spin
     *                 RHOI,TAUI,TRAI,COEF,EX,EC,EX0,
     *                 EC0,AOMAX,VPRGA,VPRGB,IAO,ILENG,NPTGRD,
     *                 L1,L2,2,.FALSE.,0)
         call EXPND(FXC(1,1),SCR,L1,0)
         call VADD(WRK1(1,L),1,SCR,1,WRK1(1,L),1,L3)
         call EXPND(FXC(1,2),SCR,L1,0)
         call VADD(WRK2(1,L),1,SCR,1,WRK2(1,L),1,L3)
      enddo                                                              !do L=1,4
      endif                                                              !if(NDFTFG.eq.1)then
C
      DO L=1,NMIC                                                        !convert to MO
         CALL DGEMM('N','N',L1,LX,L1,ONE,WRK1(1,L),L1,VA,L1,
     *              ZERO,SCR,L1)
         CALL DGEMM('T','N',LX,LX,L1,ONE,VA,L1,SCR,L1,
     *              ZERO,WRK1(1,L),L1)
         CALL DGEMM('N','N',L1,LX,L1,ONE,WRK2(1,L),L1,VA,L1,
     *              ZERO,SCR,L1)
         CALL DGEMM('T','N',LX,LX,L1,ONE,VA,L1,SCR,L1,
     *              ZERO,WRK2(1,L),L1)
      ENDDO
C
      DO L=1,NMIC
         WMTMP = CM(L)
         IF(L.LE.2) WMTMP = WMTMP*HALF

         nra = nmicro(1,1,L) - 1                                         !occupation of -r- orbital in L-th microstate; α-spin
         nrb = nmicro(1,2,L) - 1                                         !occupation of -r- orbital in L-th microstate; β-spin
         nsa = nmicro(2,1,L) - 1                                         !occupation of -s- orbital in L-th microstate; α-spin
         nsb = nmicro(2,2,L) - 1                                         !occupation of -s- orbital in L-th microstate; β-spin

         IJ=0
         DO J=1,LX                                                       !assemble Q2 matrix in MO basis
         occja = 1.d0
         if(j.eq.NOCR)occja = dble(nra)
         if(j.eq.NOCS)occja = dble(nsa)
         if(j.gt.NOCS)occja = 0.d0
         occjb = 1.d0
         if(j.eq.NOCR)occjb = dble(nrb)
         if(j.eq.NOCS)occjb = dble(nsb)
         if(j.gt.NOCS)occjb = 0.d0
          JLX = (J - 1) * LX
            DO I=1,LX
            occia = 1.d0
            if(i.eq.NOCR)occia = dble(nra)
            if(i.eq.NOCS)occia = dble(nsa)
            if(i.gt.NOCS)occia = 0.d0
            occib = 1.d0
            if(i.eq.NOCR)occib = dble(nrb)
            if(i.eq.NOCS)occib = dble(nsb)
            if(i.gt.NOCS)occib = 0.d0
               IJ = IJ+1
               TMPA = WRK1(JLX+I,L)*(occia+occja)
               TMPB = WRK2(JLX+I,L)*(occib+occjb)
C
               X2E(IJ)=X2E(IJ)+WMTMP*(TMPA+TMPB)
            END DO
         END DO

      END DO
C
      RETURN
      END SUBROUTINE SSRQ2E
C
C*MODULE REKS    *DECK BUILDFOCK
C>
C> @author     Michael Filatov
C>
C> @date       2021, Oct
C>
      SUBROUTINE BUILDFOCK(DM,FM,SCR1,SCR2,XDFT,BUFF,IBUF,HLFXC)
C>
C> DM   - density matrix (input); L2
C> FM   - Fock matrix (output); L2
C> SCR1 - scratch; L2
C> SCR2 - scratch; L2
C> BUFF - scratch: MINTMX
C> IBUF - scratch; MINTMX
C>
      USE camdft, ONLY:ALPHAC => cam_alpha, CAMMU => cam_mu, CAMFLAG
      USE lrcdft, ONLY:LCFLAG, EMU, EMU2, LRFILE
      USE mx_limits, only:mxatm, mxao, mxgrid
      USE constants, only:half
      USE prec, ONLY: dp
      IMPLICIT NONE
C
      REAL(KIND=dp), DIMENSION(137) :: BSLRD
      REAL(KIND=dp) :: BUFRAD, BUFTYP, DFTGTHR, DFTTHR, ELRD10, ELRD6,   &
     &                ELRD8, EMULT, EXENA, EXENB, EXENC, PACK2E, SUBLNG,&
     &                SUBTYP, SW0, SWOFF
      REAL(KIND=dp), DIMENSION(3,MXATM) :: C
      LOGICAL :: DCFLG, DOLRD, DSKWRK, GOPARR, LRDFLG, MASWRK, MLTINT,  &
     &           SG1
      REAL(KIND=dp), DIMENSION(20) :: DFTTYP
      INTEGER, DIMENSION(MXAO) :: IA
      INTEGER, DIMENSION(MXATM) :: IAN
      INTEGER :: IBTYP, ICH, ICUT, IDAF, IDENAO, IDFT34, IGRDTYP,       &
     &           ININTIC, INTTYP, IP, IPK, IPTIM, IR, IS, ITOL, IW,     &
     &           JANS, LABSIX, LBUFPIC, LIXIC, MASTER, ME, MUL, NA, NAT,&
     &           NAUXFUN, NAUXSHL, NAV, NB, NDCPRT, NDFTFG, NE, NHEX,   &
     &           NINTIC, NINTIX, NINTMX, NOPK, NORMF, NORMP, NPHI,      &
     &           NPHI0, NPRINT, NPROC, NQMT, NRAD, NRAD0, NSUBS, NTHE,  &
     &           NTHE0, NTUPL, NUM, NXXIC
      INTEGER, DIMENSION(950) :: IODA
      INTEGER, DIMENSION(MXGRID) :: NANGPT, NANGPT0
      REAL(KIND=dp), DIMENSION(MXATM) :: ZAN
      COMMON /DCOPT / SUBTYP, BUFTYP, SUBLNG, BUFRAD, NDCPRT, NSUBS,    &
     &                DCFLG
      COMMON /DFGRID/ DFTTHR, DFTGTHR, SWOFF, SW0, BSLRD, NDFTFG, NRAD, &
     &                NTHE, NPHI, NRAD0, NTHE0, NPHI0, NANGPT, NANGPT0, &
     &                SG1, JANS
      COMMON /DFTPAR/ DFTTYP, EXENA, EXENB, EXENC, IDFT34, NAUXFUN,     &
     &                NAUXSHL
      COMMON /DNSAO / IDENAO
      COMMON /IJPAIR/ IA
      COMMON /INFOA / NAT, ICH, MUL, NUM, NQMT, NE, NA, NB, ZAN, C, IAN
      COMMON /INT2IC/ NINTIC, ININTIC, NXXIC, LBUFPIC, LIXIC, LABSIX,   &
     &                NINTIX
      COMMON /INTFIL/ NINTMX, NHEX, NTUPL, PACK2E, INTTYP, IGRDTYP
      COMMON /IOFILE/ IR, IW, IP, IS, IPK, IDAF, NAV, IODA
      COMMON /LRDISP/ ELRD6, ELRD8, ELRD10, EMULT, LRDFLG, MLTINT, DOLRD
      COMMON /OUTPUT/ NPRINT, ITOL, ICUT, NORMF, NORMP, NOPK
      COMMON /PAR   / ME, MASTER, NPROC, IBTYP, IPTIM, GOPARR, DSKWRK,  &
     &                MASWRK
C
      LOGICAL :: HLFXC
      REAL(KIND=dp) :: XDFT
      REAL(KIND=dp), DIMENSION(*) :: BUFF, DM, FM, SCR1, SCR2
      INTEGER, DIMENSION(*) :: IBUF
C
      REAL(KIND=dp) :: EEXC, TOTELE, TOTKIN
      LOGICAL :: DUPAO, TDSKWRK
      INTEGER :: I, IDENAOSAV, II, L0, L1, L2, L3, MINTMX
C
C     Build the closed-shell Fock matrix for the given density matrix
C
C     L0 = NUMBER OF CANONICAL ORTHONORMAL VECTORS KEPT.
C
      L0 = NQMT
      L1 = NUM
      L2 = (L1*L1+L1)/2
      L3 = L1*L1
C
      MINTMX=NINTMX
      IF(NINTIC.NE.0) MINTMX=0
C
C do the 2-e part of the HF Fock matrix with the integrals from hard drive
C
      TDSKWRK = DSKWRK
      DSKWRK  = .TRUE.
C
      IF(CAMFLAG  .AND.  NDFTFG.EQ.1) THEN
         CALL SEQREW(LRFILE)
         EMU = CAMMU
         EMU2 = CAMMU*CAMMU
         CALL VCLR(SCR1,1,L2)
         CALL HSTARLC(DM,SCR1,BUFF,IBUF,NINTMX,IA)
         CALL DSCAL(L2,HALF,SCR1,1)                                      !off-diagonal elements need to be halved
         II=0
         DO I=1,L1
            II = II + I
            SCR1(II) = SCR1(II) + SCR1(II)                               !restore the diagonal elements
         ENDDO
      ENDIF
C
      IF(LCFLAG) THEN
         CALL SEQREW(LRFILE)
         CALL HSTARLC(DM,FM,BUFF,IBUF,NINTMX,IA)
      END IF
C
      IF(CAMFLAG  .AND.  NDFTFG.EQ.1) DFTTYP(3) = ALPHAC
C
C               ORDINARY INTEGRAL CONTRIBUTIONS TO FOCK MATRIX,
C               WHICH ARE ALSO NEEDED FOR RANGE-SEPARATED DFT.
C
      CALL SEQREW(IS)
      CALL HSTAR(DM,FM,BUFF,IBUF,NINTMX,IA,NOPK,.FALSE.)
C
      IF(CAMFLAG  .AND.  NDFTFG.EQ.1) CALL VADD(FM,1,SCR1,1,FM,1,L2)
C
      DSKWRK = TDSKWRK
C
      IF (NDFTFG.EQ.1) THEN
C
         DOLRD=.FALSE.
         IDENAOSAV=IDENAO
         IDENAO=1                                                        !set as if it's divide&conquer; however, DCFLG is false
C
         CALL DSCAL(L2,HALF,DM,1)                                        !scale the DM by 1/2; should be DM for α-spin only
C
         CALL DFTSET(XDFT,0,.FALSE.)
         CALL DFTEXCOR(XDFT,SCR1,SCR1,DM,DM,L1,L2,EEXC,TOTELE,TOTKIN)    !with DCFLG=.F., it does the XC part with the density matrix only
         CALL SYMH(SCR1,SCR2,IA)
         IDENAO=IDENAOSAV
C
C          BY ADDING FXC NOW WE SAVE ONE CALL TO DDI_GSUM.
C
         CALL DSCAL(L2,HALF,SCR1,1)                                      ! VAP scale the XC_SCR1 by 1/2
         CALL VADD(FM,1,SCR1,1,FM,1,L2)
C         CALL PRTRIL(SCR1,L1)
C
      END IF
C
      IF(GOPARR) THEN
         CALL DDI_GSUMF(1000,FM,L2)
      END IF
C
      CALL DAREAD(IDAF,IODA,SCR1,L2,11,0)                                !read the 1-e core hamiltonian
      CALL VADD(FM,1,SCR1,1,FM,1,L2)
C      CALL PRTRIL(SCR1,L1)
C
      RETURN
      END SUBROUTINE BUILDFOCK
