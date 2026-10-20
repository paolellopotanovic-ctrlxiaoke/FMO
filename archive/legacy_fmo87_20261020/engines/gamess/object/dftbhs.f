C  6 Jun 18 - YN -  updates for FMO 5.3
C 18 Apr 16 - YN  - updates for FMO 5.2
C 29 Sep 14 - YN  - Implement Density-Functional Tight-Binding Method
C
C*MODULE DFTBHS    *DECK DFTB_HESSIAN
C>
C>    @brief   Hessian of DFTB
C>
C>    @details Caculate Higher-Order Derivatives for DFTB.
C>             Details of implementations:
C>             - Hessian with fractional occupation number
C>               Nishimoto, Y.; Irle, S.
C>               Chem. Phys. Lett. 2017, 667, 317-321.
C>             - Third-order derivatives
C>               Nishimoto, Y.
C>               J. Chem. Phys. 2017, 146, 084101.
C>
C>    @author  Yoshio Nishimoto
C>             - Sep, 2014 - Subroutine written
C>    @date    Feb, 2016 -  Hiroya Nakata
C>             - Added FMO-DFTB Hessian
C>    @date    Feb, 2016 - Yoshio Nishimoto
C>             - Performance improvement + DFTB/PCM Hessian
C>    @date    Jun, 2016 - Yoshio Nishimoto
C>             - Bug fix for zero electron
C>    @date    Feb, 2018 - Yoshio Nishimoto
C>             - Performance improvement
C>             - Hyperpolarizability + polarizability derivative
C>             - (Higher-order coupled-perturbed)
C>             - (Experimental code of AO-based CP)
C>    @date    Oct, 2019 - Yoshio Nishimoto
C>             - Added LC-DFTB Hessian
C>
C>           --- OUTPUT ---
C>    @param ALPHA   Polarizability
C>    @param ALPDER  Polarizability derivative
C>
C
      SUBROUTINE DFTB_HESSIAN(ALPHA,ALPDER)
      USE lrcdft, ONLY: LCFLAG, EMU, EMU2, LRFILE
      use mx_limits, only: mxatm
      USE DFTBPB_MOD,ONLY: PERIOD
C
C     CONVERT DEBYE**2/AMU-ANGSTROM**2 UNIT TO KM/MOLE WITH 42.255
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION ALPHA(6),ALPDER(6,3,NAT)
C
      PARAMETER (MXSPE=10,ZERO=0.0D+00,ONE=1.0D+00)
      INTEGER, PARAMETER :: DDI_WORLD=0
      LOGICAL DFTBFL,SCC,SRSCC,DFTB3,DAMPXH,LCDFTB,ARAMAN
      LOGICAL ISGDDI,PAROUT,INITGDDI,WASGDDI,SAVGDDI,MLGDDI
C
      COMMON /DFTB  / DFTBFL,SCC,SRSCC,DFTB3,DAMPXH,LCDFTB
      COMMON /DFTBHS/ DEGTHR(2),PFREQ(2),CPCONV,OCCTHR,MODHSS,
     *                MODPRJ,MXCPIT,ARAMAN
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
      COMMON /FMCOM / X(1)
      common /fmoinf/ nfg,nlayer,natfmo,nbdfg,naotyp,nbody,nsegm
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
      COMMON /FUNCT / ENERGY,EGRAD(3,MXATM)
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON /MACHIN/ NWDVAR,MAXFM,MAXSM,LIMFM,LIMSM
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      LOGICAL GOPARR,DSKWRK,MASWRK,SAVGOP
      COMMON /GDDI/   ISCOPE,NGROUPS,MYGROUP,MEGLOB,npglob,nnglob,JBTYP,
     *                ISGDDI,PAROUT,INITGDDI,wasgddi,MLGDDI,NSUBGR,
     *                MeUniv,NPUniv,numdlb,myworld,nworlds,mogddi
      COMMON /MACHSW/ KDIAG,ICORFL,IXDR,modio,mem10,lpnt10,mem10m
      COMMON /OUTPUT/ NPRINT,ITOL,ICUT,NORMF,NORMP,NOPK
      COMMON /PCMDIM/ MXSP,MXTS,MEMPCM1,MEMPCM2,NTS
      COMMON /PCMPAR/ IPCM,NFT26,NFT27,IKREP,IEF,IP_F,NFMOPCM,IHET
      COMMON /PCMOPT/ RABI,RASC,REFPOL,THRSLS,DENSLS,WB,WA,ETA2,GD,EVAC,
     *                RHOW,PM,AREATL,AREAKP,BONDRY,OMEGA,RET,FRO,EPSINF,
     *                EPS,DR,RSOLV,VMOL,TCE,STEN,DSTEN,CMF,TABS,IDIRCT,
     *                IPCDER,IDP,ICOMP,IFIELD,ICAV,IDISP,IPRINT,IRETCAV,
     *                ICENT,IFAST,NEVAL,IEFPOL,KEEPSM,IMGABI,IMGASC,NADD
      COMMON /RUNOPT/ RUNTYP,EXETYP,NEVALS,NGLEVL,NHLEVL
      COMMON /WFNOPT/ SCFTYP,VBTYP,DFTYPE,TDDFTYP,CITYP,CCTYP,
     *                MPLEVL,MPCTYP
C
      DATA CHECK/8HCHECK   /, DEBUG/8HDEBUG   /, DBUGME/8HRHFCL   /
      DATA ROHFC/8HROHF    /
      DATA UHFC/8HUHF     /
      DATA RAMAN/8HRAMAN   /,HESSIAN/8HHESSIAN /,OPTIMIZE/8HOPTIMIZE/,
     *     SADPOINT/8HSADPOINT/
C
      CHARACTER(11) :: CHARA
      LOGICAL :: DC,FON(3),ROHF,UHF,ISGDDI0,SOME,MASSAV,FULL3RD,CPF,CPA,
     *  CPF2,CPA2,CPF3
      INTEGER :: INDSH(NAT),INDSAO(NUM+1),NDOCC(2),
     *  NPART(2),NVIRT(2),NOCC(2)
C     DOUBLE PRECISION :: HESTMP(3*NAT*3*NAT)
C
      IF(MASWRK) THEN
        SOME = NPRINT.NE.-5 .OR. EXETYP.EQ.DEBUG .OR. EXETYP.EQ.DBUGME
        SOME = .TRUE.
        IF (NFG.NE.0) SOME = NPRINT.NE.-5 .AND. NPRFMO.EQ.0
      END IF
      MASSAV = MASWRK
C
      IF (ARAMAN.AND.RUNTYP.NE.RAMAN) ARAMAN=.FALSE.
      IF (ARAMAN) THEN
        CALL VCLR(ALPHA,1,6)
        CALL VCLR(ALPDER,1,18*NAT)
      END IF
C
C     CHECK WHAT TO DO
C
      FULL3RD= .FALSE.
      IF (ARAMAN.AND.MODHSS.EQ.999) FULL3RD = .TRUE.
      IF (FULL3RD) MODHSS = 1+2+4+8+16
      IF (MODHSS.EQ.-1.AND.ARAMAN) MODHSS=7
      IF (MODHSS.EQ.-1.AND.(RUNTYP.EQ.HESSIAN.OR.RUNTYP.EQ.OPTIMIZE
     *    .OR.RUNTYP.EQ.SADPOINT)) MODHSS=0
C
      CPF  = .FALSE.
      CPA  = .FALSE.
      CPF2 = .FALSE.
      CPA2 = .FALSE.
      CPF3 = .FALSE.
      ncpf = 0
C     cpf=.true.
C     cpf2=.true.
      IF (MODHSS.NE.0.AND.MODHSS.NE.999.AND.MODHSS.NE.9999) THEN
      IF (.NOT.ARAMAN) NCPF = 3
      IF (NCPF.GT.8) CALL ABRT
      IF (MASWRK.AND.SOME) THEN
        WRITE (IW,*)
        WRITE (IW,'(" FOLLOWING PROPERTIES ARE COMPUTED (MODHSS=",I4,
     *    ")")') MODHSS
      END IF
      IF (IAND(MODHSS,  1).NE.0) THEN
C       IF (MASWRK.AND.SOME) WRITE (IW,*) "   - POLARIZABILITY"
C       CPF  = .TRUE.
        IF (MASWRK.AND.SOME) WRITE (IW,*) "   - HESSIAN"
        CPA  = .TRUE.
      END IF
      IF (IAND(MODHSS,  2).NE.0) THEN
        IF (MASWRK.AND.SOME) WRITE (IW,*) "   - HYPER-POLARIZABILITY"
        CPF  = .TRUE.
        IF (NCPF.LT.1) NCPF = 1
        IF (.NOT.ARAMAN) THEN
          CPF2 = .TRUE.
          IF (NCPF.LT.2) NCPF = 2
        END IF
      END IF
      IF (IAND(MODHSS,  4).NE.0) THEN
        IF (MASWRK.AND.SOME)
     *    WRITE (IW,*) "   - POLARIZABILITY DERIVATIVE"
        CPF  = .TRUE.
        CPA  = .TRUE.
        IF (NCPF.LT.1) NCPF = 1
      END IF
      IF (IAND(MODHSS,  8).NE.0) THEN
        IF (MASWRK.AND.SOME)
     *    WRITE (IW,*) "   - ELECTRICAL ANHARMONICITY?"
        CPF  = .TRUE.
        CPA  = .TRUE.
        IF (NCPF.LT.1) NCPF = 1
      END IF
      IF (IAND(MODHSS, 16).NE.0) THEN
        IF (MASWRK.AND.SOME) WRITE (IW,*) "   - CUBIC FORCE CONSTANT"
        CPA  = .TRUE.
      END IF
      IF (IAND(MODHSS, 32).NE.0) THEN
        IF (MASWRK.AND.SOME)
     *    WRITE (IW,*) "   - SECOND HYPER-POLARIZABILITY"
        CPF  = .TRUE.
        CPF2 = .TRUE.
        IF (NCPF.LT.3) NCPF = 3
      END IF
      IF (IAND(MODHSS, 64).NE.0) THEN
        IF (MASWRK.AND.SOME) WRITE (IW,*) "   - HYPER-RAMAN"
        CPF  = .TRUE.
        CPF2 = .TRUE.
        CPF3 = .TRUE.
        IF (NCPF.LT.3) NCPF = 3
      END IF
      IF (IAND(MODHSS,128).NE.0) THEN
        IF (MASWRK.AND.SOME) WRITE (IW,*) "   - QUARTIC FORCE CONSTANT"
        CPA  = .TRUE.
        CPA2 = .TRUE.
      END IF
C
      IF (MODHSS.GE.1.AND.MASWRK.AND.SOME) THEN
        WRITE (IW,*)
        WRITE (IW,*) "COUPLED-PERTURBED (CP) EQUATIONS TO BE SOLVED"
        NVEC = 0
        IF (CPF ) THEN
          WRITE (IW,'("   1ST-ORDER CP WRT ELECTRIC FIELD:",I5)') 3
          NVEC = NVEC + 3
        END IF
        IF (CPA ) THEN
          WRITE (IW,'("   1ST-ORDER CP WRT NUCLEAR COORD.:",I5)') 3*NAT
          NVEC = NVEC + 3*NAT
        END IF
        IF (CPF2) THEN
          WRITE (IW,'("   2ND-ORDER CP WRT ELECTRIC FIELD:",I5)') 6
          NVEC = NVEC + 6
        END IF
        IF (CPA2) THEN
          WRITE (IW,'("   2ND-ORDER CP WRT NUCLEAR COORD.:",I5)')
     *      3*NAT*(3*NAT+1)/2
          NVEC = NVEC + 3*NAT*(3*NAT+1)/2
        END IF
        IF (CPF3) THEN
          WRITE (IW,'("   3RD-ORDER CP WRT ELECTRIC FIELD:",I5)') 10
          NVEC = NVEC + 10
        END IF
        IVEC = 10
        DO ICP = 4, NCPF
          IVEC = IVEC + (ICP+1)
          WRITE (IW,'("  ",I2,"TH-ORDER CP WRT ELECTRIC FIELD:",I5)') 
     *      ICP,IVEC
          NVEC = NVEC + IVEC
        END DO
        WRITE (IW,'("   -------------------------------------")')
        WRITE (IW,'("   NUMBER OF TOTAL RESPONSE EQS.  :",I5)') NVEC
        IF (IAND(MODHSS,24).NE.0) FULL3RD=.TRUE.
      END IF
C
      ELSE IF (ARAMAN.AND.MODHSS.EQ.0) THEN
        MODHSS=7 !! POLAR + HYPER + POLDER
        CPF = .TRUE.
        CPA = .TRUE.
      ELSE
        CPA = .TRUE.
      END IF
      if (cpf3 .or. ncpf.ge.3) then
        if (maswrk) write (iw,*) "higher-order CP is out of service"
        call abrt
      end if
C
      IF (IAND(MODHSS,7).NE.7) ARAMAN=.FALSE.
C
      IF (MASWRK.AND.SOME) THEN
        WRITE (IW,*)
        WRITE (IW,'(5X,"---------------------------",5X)')
        IF (IAND(MODHSS,224).NE.0) THEN
        WRITE (IW,'(5X,"2,3,4TH DERIVATIVE FOR DFTB",5X)')
        ELSE IF (ARAMAN) THEN
        WRITE (IW,'(5X,"2ND+3RD DERIVATIVE FOR DFTB",5X)')
        ELSE
        WRITE (IW,'(5X,"SECOND  DERIVATIVE FOR DFTB",5X)')
        END IF
        WRITE (IW,'(5X,"---------------------------",5X)')
        WRITE (IW,*)
        WRITE (IW,'(5X," NOTE: PLEASE CHECK WHETHER ATOMIC WEIGHTS ARE",
     *                 " CORRECT.")')
        WRITE (IW,'(5X,"       DFTB MAY USE DIFFERENT VALUES")')
        IF (CPF) THEN
        WRITE (IW,'(5X," NOTE: ELECTRIC FIELD DERIVATIVE RESULTS (E.G.",
     *                 " POLARIZABILITIES) ARE")')
        WRITE (IW,'(5X,"       SHOWN IN A.U.")')
        END IF
      END IF
C
      FON(1) = .FALSE.
      FON(2) = .FALSE.
      FON(3) = .FALSE.
      ROHF   = .FALSE.
      UHF    = .FALSE.
C
      SAVGDDI= ISGDDI
      if(nfg.ne.0) ISGDDI=.false.
C
      ISGDDI0= ISGDDI
      IF (EXETYP.EQ.CHECK) ISGDDI0 = .FALSE.
C
      IF (SCFTYP.EQ.UHFC) THEN
        UHF = .TRUE.
        IF (.NOT.SCC) THEN
          IF (MASWRK)
     *      WRITE (IW,*) "UNRESTRICTED NCC-DFTB IS NOT SUPPORTED"
          CALL ABRT
        END IF
      END IF
C
      IF (PERIOD) THEN
        IF (MASWRK) WRITE (IW,*) "DFTB/PBC HESSIAN IS NOT YET"
        CALL ABRT
      END IF
C
      IF (MEM10.NE.0.AND.ISGDDI.AND.CPF) THEN
        IF (MASWRK) WRITE (IW,*) "MEM10 OPTION HAS TO BE REMOVED FOR ",
     *                           "CP-DFTB WRT ELECTRIC FIELDS"
        CALL ABRT
      END IF
C
      NSPIN = 1
      IF (UHF) NSPIN = 2
      NSHELL = 0
      DO I = 1, NAT
        NSHELL = NSHELL + MAXANG(ISPE(I))
      END DO
      IF (UHF.OR.SRSCC) THEN
        INDSH(1) = 0
        DO I = 1, NAT-1
          INDSH(I+1) = INDSH(I) + MAXANG(ISPE(I))
        END DO
        ISH = 0
        INDSAO(1) = 0
        DO I = 1, NAT
          DO J = 1, MAXANG(ISPE(I))
            ISH = ISH + 1
            IF (J.EQ.1) INDSAO(ISH+1) = INDSAO(ISH) + 1
            IF (J.EQ.2) INDSAO(ISH+1) = INDSAO(ISH) + 3
            IF (J.EQ.3) INDSAO(ISH+1) = INDSAO(ISH) + 5
          END DO
        END DO
      END IF
C
      L0 = NQMT
      L1 = NUM
      L2 = NUM*(NUM+1)/2
      L3 = NUM*NUM
C
      IF (L0.NE.L1) THEN
        IF (MASWRK) WRITE (IW,*) "CANNOT CALCULATE HESSIAN IF LINEARLY",
     *    "DEPENDENT ORBITALS HAVE BEEN REMOVED..."
        CALL ABRT
      END IF
C
C     CHECK THE NUMBERS OF DOUBLY AND PARTLY (PARTIALLY?) OCCUPIED AND
C     UNOCCUPIED ORBITALS
C
      CALL GOTFM(NGOTMX)
      CALL VALFM(LOADFM)
      LOCC = LOADFM + 1
      LEIG = LOCC   + L1*NSPIN
      LAST = LEIG   + L1*NSPIN
      IF (SCFTYP.EQ.ROHFC) THEN
        LWRK = LAST
        LAST = LWRK + L1
      END IF
      NEED = LAST - LOADFM - 1
      CALL GETFM(NEED)
C
      DC = .FALSE.
      IF (EXETYP.NE.CHECK) THEN
        CALL DAREAD(IDAF,IODA,X(LOCC),L1,562,0)
        CALL DAREAD(IDAF,IODA,X(LEIG),L1,17,0)
        IF (UHF) THEN
          CALL DAREAD(IDAF,IODA,X(LOCC+L1),L1,563,0)
          CALL DAREAD(IDAF,IODA,X(LEIG+L1),L1,21,0)
        END IF
        NOCC(1) = 0
        NOCC(2) = 0
        CALL DFTB_HESS_CHKOCC(X(LOCC),X(LEIG),NDOCC,NPART,
     *    NVIRT,NSPIN,L1,IW,FON,ROHF,MASWRK,SOME,SCFTYP,OCCTHR)
        IF (FON(1).OR.FON(2)) FON(3) = .TRUE.
        CALL INIFAC(10)
        NOCC(1) = NDOCC(1) + NPART(1)
        IF (NSPIN.EQ.2) NOCC(2) = NDOCC(2) + NPART(2)
C
        IF (MASWRK.AND.SOME) THEN
          WRITE (IW,*)
          DO ISPIN = 1, NSPIN
            IF (NSPIN.EQ.2) THEN
              IF (ISPIN.EQ.1) WRITE (IW,'( " *** ALPHA ***")')
              IF (ISPIN.EQ.2) WRITE (IW,'(/" *** BETA  ***")')
            END IF
            WRITE (IW,'(" NUMBER OF        OCCUPIED ORBITALS = ", I5)')
     *        NDOCC(ISPIN)
            IF (FON(ISPIN))
     *       WRITE (IW,'(" NUMBER OF PARTLY OCCUPIED ORBITALS = ", I5)')
     *          NPART(ISPIN)
            WRITE (IW,'(" NUMBER OF        VIRTUAL  ORBITALS = ", I5)')
     *        NVIRT(ISPIN)
          END DO
          WRITE (IW,'(" ------------------------------------------")')
         WRITE (IW,'(" NUMBER OF TOTAL BASIS FUNCTIONS    = ", I5)') NUM
          IF (FON(3).AND.ROHF) THEN
            WRITE
     *      (IW,'(" ... FRACTIONAL OCCUPATION NUMBERS (FON) ARE USED")')
          WRITE (IW,'(" BUT ZERO ELECTRONIC TEMPERATURE. ASSUME ROHF")')
          END IF
          WRITE (IW,*)
        END IF
      ELSE
        !! THESE VARIABLES ARE ONLY GUESS, BECAUSE WE DON'T KNOW THESE
        !! VALUES UNLESS SCF IS DONE
        NDOCC(1) = NA
        NPART(1) = 0
        NVIRT(1) = NUM - NA
        IF (UHF) THEN
          NDOCC(2) = NB
          NPART(2) = 0
          NVIRT(2) = NUM - NB
        END IF
      END IF
      CALL RETFM(NEED)
C
C     ----- ESTIMATE REQUIREMENT FOR AO DERIVATIVE STORAGE -----
C
      CALL VALFM(LOADFM)
      LWRK = LOADFM + 1
      LAST = LWRK   + NAT
      NEED = LAST - LOADFM - 1
      CALL GETFM(NEED)
      MAXBND = 0
      CALL DFTB_DERMAT(1,NAT,C,INDSH,INDSAO,MAXANG,ISPE,
     *                 X(LWRK),SKCUT2,LNMAT,MAXBND,X(LAST),X(LAST),
     *                 SRSCC,UHF)
      CALL RETFM(NEED)
C
C     ESTIMATE MEMORY REQUIREMENT
C
      NDIM = NAT
      IF (SRSCC) NDIM = NSHELL
C
      MEMNCC = NAT*NAT*9 + NAT*9 + L1*2 + L2*2 + L3*3 +
     *  NAT*(NAT+1)/2 + NDIM*(NDIM+1)/2*2
      MEMNCC = MEMNCC + NAT !! LMULA ... is it needed for NCC?
C
      IF (SCC) THEN
        MEMSCC = NDIM*(NDIM+1)/2*2 + NAT
        IF (NSPIN.EQ.2) MEMSCC = MEMSCC + L1*2 + L2*1 + L3*3
        IF (UHF) MEMSCC = MEMSCC + 5*NSHELL
        IF (SRSCC.AND..NOT.UHF) MEMSCC = MEMSCC + NSHELL*4
        !! ESTIMATED MEMORY USED FOR BROYDEN MIXING
        NDIMBR = NAT
        IF (SRSCC) NDIMBR = NSHELL
        IF (FON(3)) NDIMBR = NDIMBR + 1
        LL = NDIMBR*MXCPIT !! MXITER = 50
        LMX1 = MXCPIT
        LMX2 = (LMX1*LMX1+LMX1)/2
        LMX3 = LMX1*LMX1
        MEMBR = 7*NDIMBR+2*LL+2*LMX1+2*LMX2+2*LMX3
      ELSE
        MEMSCC = 0
        MEMBR  = 0
      END IF
C
      IF (DFTB3) THEN
        MEMDFTB3 = NDIM*NDIM*2
      ELSE
        MEMDFTB3 = 0
      END IF
C
      !! MEMORY FOR CP-DFTB
      IF (SRSCC.OR.UHF) THEN
        NDIMSH = NSHELL
      ELSE
        NDIMSH = NAT
      END IF
      MAXOCC = MAX(NOCC(1),NOCC(2))
      MAXVIR = MAX(NPART(1)+NVIRT(1),NPART(2)+NVIRT(2))
      MEMCP = (L1+L2*7+NDIMSH*2)*NSPIN + L2*1 + L3*4
     *  + MAXOCC*MAXVIR*NSPIN
      IF (SCC) MEMCP = MEMCP + NDIMSH*(2+NSPIN) + L2
      IF (FON(3)) MEMCP = MEMCP + L1*2
      IF (UHF) MEMCP = MEMCP + NDIMSH
      IF (DFTB3.AND.NFG.NE.0) MEMCP = MEMCP + NAT*2
C
      !! MEMORY FOR PCM
      IF (IPCM.NE.0) THEN
        NDASC = 0
        if (isgddi) then
          NTSPAR=(NTS-1)/npglob+1
        else
          NTSPAR=(NTS-1)/nproc+1
        end if
        LFPNT = LAST
        CALL INDQPROD(0,NTSPAR,LFPNT,NDASC,DUM,DUM,DUM,0,0)
        MEMPCM  = NDASC + NTS*6
        IF (IPCDER.EQ.3) MEMPCM = MEMPCM + 60*MXTS + 21*MXTS/NWDVAR+1
        MEMPCM = MAX(MEMPCM,3*NAT+9*(NAT*NAT+NAT)/2)
      ELSE
        MEMPCM = 0
      END IF
C
      IF (ARAMAN) THEN
        MEMAR = (NDIMSH+L1*2+L2*2+L3)*NSPIN*3 + L3*NSPIN
      ELSE
        MEMAR = 0
      END IF
C
      IF (FULL3RD) THEN
        MEMF3 = L2*NSPIN*NAT*3
      ELSE
        MEMF3 = 0
      END IF
C
      !! MEMORY FOR AO DERIVATIVE
      MEMDIR  = NAT*(MAXBND+1)*4 + LNMAT*2
      IF (NFG.NE.0.AND.NBDFG.NE.0) MEMDIR = MEMDIR + L2*NAT*NSPIN*6
C
      MEMCORE = MEMNCC + MEMBR + MEMCP + MEMSCC + MEMDFTB3
     *        + MEMPCM + MEMDIR+ MEMF3
C
      IF (MASWRK.AND.SOME) THEN
        WRITE(IW,'(" ESTIMATED MEMORY REQUIREMENT (UNIT IN WORDS)")')
        WRITE(IW,'("   FOR DFTB1                  = ",I12)') MEMNCC
        IF (SCC) THEN
        WRITE(IW,'("   FOR DFTB2                  = ",I12)') MEMSCC
        IF (DFTB3)
     *  WRITE(IW,'("   FOR DFTB3                  = ",I12)') MEMDFTB3
        END IF
        WRITE(IW,'("   FOR AO DERIVATIVE MATRIX   = ",I12)') MEMDIR
        WRITE(IW,'("   FOR COUPLED PERTURBED DFTB = ",I12)') MEMCP
        IF (IPCM.NE.0)
     *  WRITE(IW,'("   FOR PCM                    = ",I12)') MEMPCM
        IF (ARAMAN)
     *  WRITE(IW,'("   FOR 3RD-ORDER DERIVATIVES  = ",I12)') MEMAR
        IF (FULL3RD)
     *  WRITE(IW,'("   FOR FULL 3RD  DERIVATIVES  = ",I12)') MEMF3
        WRITE(IW,'("   FOR BROYDEN MIXING         = ",I12)') MEMBR
        WRITE(IW,'(" -------------------------------------------")')
        WRITE(IW,'("   MEMORY REQUIRED            = ",I12)') MEMCORE
        WRITE(IW,'("   CURRENT AVAILABLE MEMORY   = ",I12,
     *             " WORDS")') NGOTMX
      END IF
C
      IF (MASWRK.AND.SOME)  WRITE(IW,*)
C     DISK = .FALSE.
      IF (FULL3RD.AND.MEMCORE.GT.NGOTMX) THEN
        MEMCORE = MEMCORE - MEMF3
C       DISK = .TRUE.
        write (*,*) "not yet"
        call abrt
      END IF
      IF (MEMCORE.GT.NGOTMX) THEN
        IF (MASWRK) WRITE(IW,*) "NOT ENOUGH MEMORY"
        CALL ABRT
      END IF
C
      CALL VALFM(LOADFM)
      NDIM = NAT
      IF (SRSCC) NDIM = NSHELL
      LHESSIAN     = LOADFM   + 1
      LDIPDER      = LHESSIAN + NAT*NAT*9
      LOCC         = LDIPDER  + NAT*9
      LEIG         = LOCC     + L1*NSPIN
      LD           = LEIG     + L1*NSPIN
      LD2          = LD       + L2*NSPIN
      LS           = LD2      + L2*NSPIN
      LVEC         = LS       + L2
      LWRK         = LVEC     + L3*NSPIN
      LWRK1        = LWRK     + L3*NSPIN !! U MATRIX
      LDIST        = LWRK1    + L3*NSPIN !! UINIT
      LGRAD        = LDIST    + NAT*(NAT+1)/2
      LHESS        = LGRAD    + NDIM*(NDIM+1)/2
      LAST         = LHESS    + NDIM*(NDIM+1)/2
      IF (UHF) THEN
        LSHIFTS    = LAST
        LSHIFTSH   = LSHIFTS  + NSHELL
        LSHIFTSA   = LSHIFTSH + NSHELL
        LSHIFTSB   = LSHIFTSA + NSHELL
        LMULS      = LSHIFTSB + NSHELL
        LAST       = LMULS    + NSHELL
      END IF
      LMULA        = LAST
      LAST         = LMULA    + NAT
C
      IF (SCC) THEN
        LGAMMA     = LAST
        LGAMDER    = LGAMMA   + NDIM*(NDIM+1)/2
        LSHIFT     = LGAMDER  + NDIM*(NDIM+1)/2
        LAST       = LSHIFT   + NAT
        IF (SRSCC.AND..NOT.UHF) THEN
          LSHIFTSH = LAST
          LSHIFTSA = LSHIFTSH + NSHELL
          LSHIFTSB = LSHIFTSA + NSHELL
          LMULS    = LSHIFTSB + NSHELL
          LAST     = LMULS    + NSHELL
          LSHIFTS  = LAST
        END IF
        IF (DFTB3) THEN
          LGAMMA3  = LAST
          LGAMDER3 = LGAMMA3  + NDIM*NDIM
          LAST     = LGAMDER3 + NDIM*NDIM
        ELSE
          LGAMMA3  = LAST
          LGAMDER3 = LAST
        END IF
      ELSE
        LGAMMA     = LAST
        LGAMDER    = LAST
        LSHIFT     = LAST
        LGAMMA3    = LAST
        LGAMDER3   = LAST
      END IF
C
      IF (LCDFTB) THEN
        LGAMLC     = LAST
        LGAMDERLC  = LGAMLC    + NDIM*NDIM
        LAST       = LGAMDERLC + NDIM*(NDIM+1)/2
      ELSE
        LGAMLC     = LAST
        LGAMDERLC  = LAST
      END IF
C
      LINDBND    = LAST
      LDERMAT    = LINDBND  + NAT*(MAXBND+1)*4
      LAST       = LDERMAT  + LNMAT*2
      IF (NFG.NE.0.AND.NBDFG.NE.0) THEN
        LHOPAO   = LAST
        LAST     = LHOPAO   + L2*NAT*3*NSPIN
      END IF
C
      NEED         = LAST - LOADFM - 1
C
C     ----- GET MEMORY -----
C
      CALL GETFM(NEED)
      IF (MASWRK.AND.SOME) THEN
        WRITE (IW,'(" ...... ALLOCATED MEMORY = ",I10," WORDS")') NEED
      END IF
C
      IF (EXETYP.EQ.CHECK) THEN
        IF (MASWRK) THEN
          WRITE (IW,'(" SUGGESTED MWORDS = ",I5)')
     *      INT(MEMCORE/DBLE(1.0D+06))+1
        END IF
        GO TO 20
      END IF
C
C     ----- RESTORE VECTORS AND ARRAYS OBTAINED IN SCF -----
C
      !! EIGENVECTORS
      CALL DAREAD(IDAF,IODA,X(LVEC),L3,15,0)
      !! DENSITY MATRIX
      CALL DAREAD(IDAF,IODA,X(LD),L2,16,0)
      IF (ROHF.OR.SCFTYP.EQ.ROHFC) THEN
        CALL DAREAD(IDAF,IODA,X(LWRK1),L2,20,0)
        CALL VADD(X(LD),1,X(LWRK1),1,X(LD),1,L2)
      END IF
      !! EIGENVALUES
      CALL DAREAD(IDAF,IODA,X(LEIG),L1,17,0)
      !! OCCUPATION NUMBER
      CALL DAREAD(IDAF,IODA,X(LOCC),L1,562,0)
      !! OVERLAP
      CALL DAREAD(IDAF,IODA,X(LS),L2,12,0)
      !! ATOMIC MULLIKEN POPULATION
      CALL DAREAD(IDAF,IODA,X(LMULA),NAT,556,0)
      !! SHIFT MATRIX
      IF (SCC) THEN
        CALL DAREAD(IDAF,IODA,X(LSHIFT),NAT,559,0)
        IF (SRSCC) CALL DAREAD(IDAF,IODA,X(LSHIFTSH),NSHELL,560,0)
        IF (NFG.NE.0) THEN
          !! Due to the change on October, 2019
          !! Subtract the ESP shift
          CALL DAREAD(IDAF,IODA,X(LWRK),NAT,566,0) !! ESP for DFTB2
          IF (DFTB3) THEN
            CALL DAREAD(IDAF,IODA,X(LWRK1),NAT*2,567,0) !! ESP for DFTB3
            DO I = 1, NAT
              DQ = X(LMULA+I-1) - ZREF(I)
              X(LWRK+I-1) = X(LWRK+I-1)
     *                    + 2.0D+00*X(LWRK1+I-1)*DQ + X(LWRK1+I-1+NAT)
            END DO
          END IF
          CALL DAXPY(NAT,-ONE,X(LWRK),1,X(LSHIFT),1)
        END IF
      END IF
      IF (UHF.OR.SRSCC) THEN
        CALL DAREAD(IDAF,IODA,X(LMULS),NSHELL,557,0)
        ISH = 0
        DO I = 1, NAT
          DO J = 1, MAXANG(ISPE(I))
            ISH = ISH + 1
            X(LMULS+ISH-1) = X(LMULS+ISH-1) - QREFL(J,ISPE(I))
          END DO
        END DO
      END IF
      IF (UHF) THEN
        CALL DAREAD(IDAF,IODA,X(LVEC+L3),L3,19,0)
        CALL DAREAD(IDAF,IODA,X(LD+L2),L2,20,0)
        CALL DAREAD(IDAF,IODA,X(LEIG+L1),L1,21,0)
        CALL DAREAD(IDAF,IODA,X(LOCC+L1),L1,563,0)
        CALL DAREAD(IDAF,IODA,X(LSHIFTS),NSHELL,561,0)
      END IF
C
      DC = .FALSE.
      IF (IDFTBD.EQ.1.OR.IDFTBD.EQ.2.OR.IDFTBD.EQ.4) DC = .TRUE.
C
      CALL VCLR(X(LDIST),1,NAT*(NAT+1)/2)
      CALL VCLR(X(LGRAD),1,NAT*(NAT+1)/2)
      CALL VCLR(X(LHESS),1,NAT*(NAT+1)/2)
      CALL VCLR(X(LHESSIAN),1,3*NAT*3*NAT)
      CALL VCLR(X(LDIPDER),1,9*NAT)
C
C     ----- RESTORE PCM GRADIENT AND HESSIAN -----
C
      IF (IPCM.NE.0) THEN
        CALL GOTFM(NGOTMX)
        CALL VALFM(LOADFM)
        NEGH = 3*NAT + 9*(NAT*NAT+NAT)/2
        LEGH = LOADFM + 1
        LEG  = LEGH
        LEH  = LEG   + 3*NAT
        LAST = LEH   + 9*(NAT*NAT+NAT)/2
        NEED2= LAST - LOADFM - 1
        CALL GETFM(NEED2)
        !! gradient
        CALL DAREAD(IDAF,IODA,X(LEGH),NEGH,67,0)
        CALL DAXPY(3*NAT,1.0D+00,X(LEGH),1,EGRAD,1)
C       CALL DFTB_PRJGRD(EGRAD,C,RUNTYP,NAT)
        CALL DAWRIT(IDAF,IODA,EGRAD,3*NAT,3,0)
        !! hessian
        CALL CPWAB0(X(LHESSIAN),X(LEH),NAT)
C       CALL CPFCM(X(LHESSIAN),NAT)
        IF (GOPARR) CALL DSCAL(9*NAT*NAT,ONE/DBLE(NPROC),X(LHESSIAN),1)
        CALL RETFM(NEED2)
      END IF
C
C     ----- CALCULATE INTER-ATOMIC DISTANCES -----
C
      NSEQ = 0
      DO I = 1, NAT
        DO J = 1, I
          NSEQ = NSEQ + 1
          DIST = SQRT( (C(1,I)-C(1,J))**2
     &               + (C(2,I)-C(2,J))**2
     &               + (C(3,I)-C(3,J))**2 )
          X(LDIST+NSEQ-1) = DIST
        END DO
      END DO
C
      SAVGOP=GOPARR
      IF(GOPARR.AND.NFG.NE.0) GOPARR=.false.
C
C     ----- CALCULATE GAMMA AND ITS DERIVATIVE OF DFTB2 AND DFTB3 -----
C
      IF (SCC) THEN
        CALL DFTB_GAMMA(DAMPXHE,HUBBL,X(LGAMMA),X(LDIST),NAT,NSPE,
     *    ISPE,MAXANG,SPE,DAMPXH,SRSCC)
        !! CALCULATE GAMMA DERIVATIVES
        IF (SRSCC) THEN
          NN = 0
          DO I = 1, NSPE
            NN = NN + MAXANG(I)
          END DO
          NNN = NSHELL
        ELSE
          NN = NSPE
          NNN = NAT
        END IF
        CALL DFTB_GAMDER(NAT,NSPE,HUBBL,ISPE,X(LGAMDER),DAMPXH,
     *    DAMPXHE,SPE,X(LDIST),DFTB3,X(LGAMDER3),HUBDER,SRSCC,MAXANG,
     *    NNN,NN)
        IF (DFTB3) THEN
          !! DFTB_GAMDER GIVES GAMDER3 VALUES SCALED BY 1/3
          CALL DSCAL(NNN*NNN,3.0D+00,X(LGAMDER3),1)
          !! NOW, GAMMA3 VALUES ARE NOT SCALED>
          CALL DFTB_GAMMA3(NAT,HUBBL,X(LGAMMA3),X(LDIST))
        END IF
      END IF
C
      IF (LCDFTB) THEN
        CALL DFTB_LCGAMMADER(EMU,HUBB,NAT,ISPE,X(LDIST),
     *                       X(LWRK),X(LGAMDERLC))
        CALL CPYTSQ(X(LWRK),X(LGAMLC),NAT,1)
      END IF
C
C     ----- NOTHING SHOULD BE DONE IF THERE ARE NO ELECTRONS -----
C
      IF (NA.EQ.0) THEN
        DO I = 1, NAT
          X(LMULA+I-1) = X(LMULA+I-1) - ZREF(I)
        END DO
C
        !! DIPOLE DERIVATIVE
        DO NI = 1, NAT
          DO KI = 1, 3
            X(LDIPDER+KI-1+(KI-1)*3+(NI-1)*NAT*3) = - X(LMULA+NI-1)
          END DO
        END DO
        CALL DSCAL(9*NAT,4.803204749D+00,X(LDIPDER),1)
        GO TO 30
      END IF
C
      DO I = 1, NAT
        X(LMULA+I-1) = X(LMULA+I-1) - ZREF(I)
      END DO
C
C     ----- CALCULATE AO DERIVATIVES -----
C
      CALL DFTB_DERMAT(2,NAT,C,INDSH,INDSAO,MAXANG,ISPE,
     *                 X(LWRK),SKCUT2,LNMAT,MAXBND,X(LINDBND),
     *                 X(LDERMAT),SRSCC,UHF)
      IF (EXETYP.EQ.CHECK) GO TO 20
C
C     ----- CONVERT ATOM-RESOLVED SHIFT TO SHELL-RESOLVED SHIFT -----
C
      IF (.NOT.SRSCC.AND.UHF) THEN
        CALL DFTB_SHIFT_ATOSH(X(LSHIFT),X(LSHIFTSH),NAT,NSHELL,NSPE,
     *    ISPE,MAXANG)
      END IF
C
C     ----- MAKE SHIFT MATRICES INCLUDING SPIN CONTRIBUTION -----
C     SHIFTSA = (SHELL-RESOLVED SHIFT) + (SPIN SHIFT)
C     SHIFTSB = (SHELL-RESOLVED SHIFT) - (SPIN SHIFT)
C
      IF (UHF) THEN
        CALL DCOPY(NSHELL,X(LSHIFTSH),1,X(LSHIFTSA),1)
        CALL DCOPY(NSHELL,X(LSHIFTSH),1,X(LSHIFTSB),1)
        CALL DAXPY(NSHELL, ONE,X(LSHIFTS),1,X(LSHIFTSA),1)
        CALL DAXPY(NSHELL,-ONE,X(LSHIFTS),1,X(LSHIFTSB),1)
      ELSE IF (SRSCC) THEN
        CALL DCOPY(NSHELL,X(LSHIFTSH),1,X(LSHIFTSA),1)
        CALL VCLR(X(LSHIFTSB),1,NSHELL)
      END IF
C
C     ----- HOP CONTRIBUTIONS FOR FMO-DFTB -----
C
      IF(NFG.ne.0.and.NBDFG.ne.0) then
c       IF(MASWRK)  write(*,*) "Entering HOP DER in CPHF"
        CALL VCLR(X(LHOPAO),1,L2*NAT*3*NSPIN)
        CALL VALFM(LOADFM)
        LPS   = LOADFM + 1
        LAST  = LPS    + L2*NAT*3
        NEED2 = LAST   - LOADFM - 1
        CALL GETFM(NEED2)
        CALL DFTB_EXTSDER(X(LPS),L2,1,NAT)
        GOPARR = SAVGOP
        CALL HOPDER(2,0,0,0,L1,L2,
     *              x(lfmoda),X(LHOPAO),DUM,X(LPS))
        CALL DAREAD(IDAF,IODA,X(LPS),L1,17,0)
        IF(GOPARR) GOPARR=.false.
C
        NCOUNT = 0
        THRUSH = 1.0D+04
        DO III = 1, NQMT
          IF(X(LPS+III-1) .GT. THRUSH) NCOUNT = NCOUNT + 1
        END DO
C
        !! NEED TO PROJECT HOP MO CONTRIBUTIONS OUT
        !! CONSTRUCT (SCC)_ma = S*C*C^T
        !! = \sum_{n} \sum_{i}^{occ+virt} Smn*Cni*Cai
        !! where m,n,a = AO, i = MO
        CALL DFTB_CPYTSQ(X(LS),X(LWRK),L1,L1)
        CALL DGEMM('N','T',L1,L1,L1-NCOUNT,ONE,X(LVEC),L1,X(LVEC),L1,
     *             ZERO,X(LWRK1),L1)
        CALL DGEMM('N','N',L1,L1,L1,ONE,X(LWRK),L1,X(LWRK1),L1,
     *             ZERO,X(LPS),L1)
        !! \tilde{P}mn = \sum_{ab} (SCC)ma*Pab*(SCC)nb^T
        DO I = 1, 3*NAT
          CALL CPYTSQ(X(LHOPAO+L2*(I-1)),X(LWRK),L1,1)
          CALL DGEMM('N','N',L1,L1,L1,ONE,
     *               X(LPS),L1,X(LWRK),L1,
     *               ZERO,X(LWRK1),L1)
          CALL DGEMM('N','T',L1,L1,L1,ONE,
     *               X(LWRK1),L1,X(LPS),L1,ZERO,
     *               X(LWRK),L1)
          CALL CPYSQT(X(LWRK),X(LHOPAO+L2*(I-1)),L1,1)
        END DO
        CALL RETFM(NEED2)
      end if
C
C     ----- SOLVE COUPLED PERTURBED DFTB EQUATIONS -----
C
   20 CONTINUE
      NGAM = NAT
      IF (SRSCC) NGAM = NSHELL
      NDIM = NAT
      IF (SRSCC.OR.UHF) NDIM = NSHELL
C     CALL DCOPY(3*NAT*3*NAT,X(LHESSIAN),1,HESTMP,1)
      CALL DFTB_CP(L0,L1,L2,NDOCC,NVIRT,NPART,NOCC,
     *  X(LWRK),X(LVEC),X(LEIG),X(LOCC),X(LSHIFT),
     *  X(LGAMMA),X(LGAMDER),X(LS),X(LMULA),X(LDIST),X(LWRK1),
     *  X(LHESSIAN),X(LD),X(LHOPAO),X(LGAMMA3),
     *  X(LGAMDER3),X(LDIPDER),FON,ROHF,UHF,NGAM,
     *  X(LSHIFTSA),X(LSHIFTSB),INDSH,
     *  INDSAO,NSPIN,X(LMULS),NDIM,HUBB,HUBBL,
     *  SAVGOP,SOME,ALPHA,ALPDER,
     *  MAXBND,X(LINDBND),X(LDERMAT),FULL3RD,CPF,CPA,CPF2,CPA2,CPF3,
     *  NCPF,NDIMSH,X(LGAMLC),X(LGAMDERLC),NSHELL)
      IF (EXETYP.EQ.CHECK) GO TO 100
C
      IF (ISGDDI0) CALL GDDI_SCOPE(DDI_WORLD)
      IF (.NOT.CPA) GO TO 100
C
C     WRITE (*,*) "AFTER DFTB_CP"
C     DO I = 1, 3*3*NAT*NAT
C       WRITE (*,'(i3,x,f20.10)') i,x(lhessian+i-1)
C     END DO
C     WRITE (*,*) "DFTB_CP CONTRIBUTION"
C     DO I = 1, 3*3*NAT*NAT
C       WRITE (*,'(i3,x,f20.10)') i,x(lhessian+i-1)-hestmp(i)
C     END DO
C
C     CALL DCOPY(3*NAT*3*NAT,X(LHESSIAN),1,HESTMP,1)
C     IF (MASWRK.AND.SOME) THEN
C       WRITE (IW,'(" CALCULATE DERIVATIVES OF H AND S ...")',
C    *  ADVANCE='NO')
C       CALL CPU_TIME(T0)
C       CALL SYSTEM_CLOCK(IC0,ICR,ICM)
C     END IF
C
      !! CALCULATE ENERGY-WEIGHTED DENSITY MATRIX
      CALL VCLR(X(LD2),1,L2)
      CALL DFTB_DMTEIG(X(LD2),X(LVEC),X(LOCC),X(LEIG),NA,L1,L1)
      IF (UHF) THEN
        CALL VCLR(X(LD2+L2),1,L2)
        CALL DFTB_DMTEIG(X(LD2+L2),X(LVEC+L3),X(LOCC+L1),X(LEIG+L1),
     *    NB,L1,L1)
      END IF
      IF (LCDFTB) THEN
        !! Destroy X(LGAMMA)
        CALL DFTB_LCGAMMA(EMU,HUBB,NAT,ISPE,X(LDIST),X(LGAMMA))
C
        CALL DAREAD(IDAF,IODA,X(LWRK1),L2,16,0)
        CALL VCLR(X(LWRK),1,L2)
        CALL DFTB_MAKEREFD(NAT,X(LWRK))
        CALL DAXPY(L2,-ONE/DBLE(NSPIN),X(LWRK),1,X(LWRK1),1)
        CALL CPYTSQ(X(LWRK1),X(LWRK),L1,1)
C
        CALL DFTB_LCGRAD1(0,NAT,L1,ISPE,MAXANG,IND,X(LWRK),X(LWRK),
     *                    X(LS),X(LGAMMA),X(LWRK1))
        CALL DAXPY(L2,-DBLE(NSPIN),X(LWRK1),1,X(LD2),1)
C
        IF (NSPIN.EQ.2) THEN
          CALL DAREAD(IDAF,IODA,X(LWRK1),L2,20,0)
          CALL VCLR(X(LWRK),1,L2)
          CALL DFTB_MAKEREFD(NAT,X(LWRK))
          CALL DAXPY(L2,-ONE/DBLE(NSPIN),X(LWRK),1,X(LWRK1),1)
          CALL CPYTSQ(X(LWRK1),X(LWRK),L1,1)
C
          CALL DFTB_LCGRAD1(0,NAT,L1,ISPE,MAXANG,IND,X(LWRK),X(LWRK),
     *                      X(LS),X(LGAMMA),X(LWRK1))
          CALL DAXPY(L2,-DBLE(NSPIN),X(LWRK1),1,X(LD2+L2),1)
        END IF
      END IF
C
      CALL DFTB_D2HS(X(LVEC),X(LVEC+L2),X(LD),X(LD2),
     *  X(LSHIFT),X(LSHIFTSA),X(LSHIFTSB),
     *  C,X(LHESSIAN),NAT,L2,NSPIN,IND,NSHELL,NSPE,ISPE,MAXANG,
     *  INDSAO,SCC,SRSCC,UHF)
C
C     IF (MASWRK.AND.SOME) THEN
C       CALL CPU_TIME(T1)
C       CALL SYSTEM_CLOCK(IC1,ICR,ICM)
C       WRITE (IW,'(" FINISHED IN ",F8.2," (",F8.2,") SECONDS")')
C    *  T1-T0,DBLE(IC1-IC0)/DBLE(ICR)
C     END IF
C     WRITE (*,*) "AFTER D2HS"
C     DO I = 1, 3*3*NAT*NAT
C       WRITE (*,'(i3,x,f20.10)') i,x(lhessian+i-1)
C     END DO
C     WRITE (*,*) "D2HS CONTRIBUTION"
C     DO I = 1, 3*3*NAT*NAT
C     WRITE (*,'(i3,x,f20.10)') i,x(lhessian+i-1)-hestmp(i)
C     END DO
C     CALL CPU_TIME(T6)
C
C     CALL DCOPY(3*NAT*3*NAT,X(LHESSIAN),1,HESTMP,1)
   30 CONTINUE
      IF (SCC) THEN
C       IF (MASWRK.AND.SOME) THEN
C         WRITE (IW,'(" CALCULATE GAMMA DERIVATIVE TERM  ...")',
C    *      ADVANCE='NO')
C         CALL CPU_TIME(T0)
C         CALL SYSTEM_CLOCK(IC0,ICR,ICM)
C       END IF
        IF (SRSCC) THEN
          NN = 0
          DO I = 1, NSPE
            NN = NN + MAXANG(I)
          END DO
          CALL DFTB_GAM_DER2(NAT,NSPE,HUBBL,C,X(LDIST),X(LGRAD),
     *      X(LHESS),X(LWRK),X(LWRK1),X(LMULS),X(LHESSIAN),DFTB3,
     *      HUBDER,DAMPXHE,ISPE,MAXANG,NSHELL,NSHELL,NN,SPE,SRSCC,
     *      DAMPXH,X(LMULA))
        ELSE
          NN = NSPE
          CALL DFTB_GAM_DER2(NAT,NSPE,HUBBL,C,X(LDIST),X(LGRAD),
     *      X(LHESS),X(LWRK),X(LWRK1),X(LMULA),X(LHESSIAN),DFTB3,
     *      HUBDER,DAMPXHE,ISPE,MAXANG,NSHELL,NAT,NN,SPE,SRSCC,
     *      DAMPXH,X(LMULA))
        END IF
C       IF (MASWRK.AND.SOME) THEN
C         CALL CPU_TIME(T1)
C         CALL SYSTEM_CLOCK(IC1,ICR,ICM)
C         WRITE (IW,'(" FINISHED IN ",F8.2," (",F8.2,") SECONDS")')
C    *    T1-T0,DBLE(IC1-IC0)/DBLE(ICR)
C       END IF
      END IF
C     END IF
C     WRITE (*,*) "AFTER GAM_DER2"
C     DO I = 1, 3*3*NAT*NAT
C       WRITE (*,'(i3,x,f20.10)') i,x(lhessian+i-1)
C     END DO
C     WRITE (*,*) "GAM_DER2 CONTRIBUTION"
C     DO I = 1, 3*3*NAT*NAT
C       WRITE (*,'(i3,x,f20.10)') i,x(lhessian+i-1)-hestmp(i)
C     END DO
C
C     ----- SECOND-ORDER DERIVATIVE OF LC-GAMMA -----
C
      IF (LCDFTB) THEN
        CALL DAREAD(IDAF,IODA,X(LWRK1),L2,16,0)
        CALL VCLR(X(LWRK),1,L2)
        CALL DFTB_MAKEREFD(NAT,X(LWRK))
        CALL DAXPY(L2,-ONE/DBLE(NSPIN),X(LWRK),1,X(LWRK1),1)
        CALL CPYTSQ(X(LWRK1),X(LWRK),L1,1)
C
        CALL DAREAD(IDAF,IODA,X(LS),L2,12,0)
        CALL CPYTSQ(X(LS),X(LWRK1),L1,1)
C
C       call dcopy(9*nat*nat,x(lhessian),1,hestmp,1)
        CALL DFTB_LCGRAD2(2,DBLE(NSPIN),NAT,L1,ISPE,MAXANG,IND,
     *                    X(LWRK),X(LWRK),X(LWRK1),X(LWRK1),X(LDIST),
     *                    X(LAST),C,X(LAST),X(LHESSIAN))
        IF (NSPIN.EQ.2) THEN
          CALL DAREAD(IDAF,IODA,X(LWRK1),L2,20,0)
          CALL VCLR(X(LWRK),1,L2)
          CALL DFTB_MAKEREFD(NAT,X(LWRK))
          CALL DAXPY(L2,-ONE/DBLE(NSPIN),X(LWRK),1,X(LWRK1),1)
          CALL CPYTSQ(X(LWRK1),X(LWRK),L1,1)
C
          CALL DAREAD(IDAF,IODA,X(LS),L2,12,0)
          CALL CPYTSQ(X(LS),X(LWRK1),L1,1)
C
          CALL DFTB_LCGRAD2(2,DBLE(NSPIN),NAT,L1,ISPE,MAXANG,IND,
     *                      X(LWRK),X(LWRK),X(LWRK1),X(LWRK1),X(LDIST),
     *                      X(LAST),C,X(LAST),X(LHESSIAN))
        END IF
C       WRITE (*,*) "GAM_DER2 CONTRIBUTION"
C       DO I = 1, 3*3*NAT*NAT
C         WRITE (*,'(i3,x,f20.10)') i,x(lhessian+i-1)-hestmp(i)
C       END DO
      END IF
C
C     ----- ADD EREP CONTRIBUTION -----
C
C     CALL DCOPY(3*NAT*3*NAT,X(LHESSIAN),1,HESTMP,1)
C     IF (MASWRK.AND.SOME) THEN
C       WRITE (IW,'(" CALCULATE DERIVATIVES OF EREP    ...")',
C    *  ADVANCE='NO')
C       CALL CPU_TIME(T0)
C       CALL SYSTEM_CLOCK(IC0,ICR,ICM)
C     END IF
      CALL DFTB_EREP_HESSIAN(2,X(LDIST),X(LGRAD),X(LHESS),X(LWRK),
     *                       X(LHESSIAN),X(LREPINTV(1,1)),
     *                       X(LREPSHORT(1,1)),X(LREPCOEFF(1,1)))
C     IF (MASWRK.AND.SOME) THEN
C       CALL CPU_TIME(T1)
C       CALL SYSTEM_CLOCK(IC1,ICR,ICM)
C       WRITE (IW,'(" FINISHED IN ",F8.2," (",F8.2,") SECONDS")')
C    *  T1-T0,DBLE(IC1-IC0)/DBLE(ICR)
C     END IF
C     goto 100
C     WRITE (*,*) "AFTER EREP"
C     DO I = 1, 3*3*NAT*NAT
C       WRITE (*,'(i3,x,f20.10)') i,x(lhessian+i-1)
C     END DO
C     WRITE (*,*) "EREP CONTRIBUTION"
C     DO I = 1, 3*3*NAT*NAT
C       WRITE (*,'(i3,x,f20.10)') i,x(lhessian+i-1)-hestmp(i)
C     END DO
C
C     ----- ADD DISPERSION CONTRIBUTION -----
C
      IF (DC) THEN
        IF (IDFTBD.NE.3) THEN
C         IF (MASWRK.AND.SOME) THEN
C           WRITE (IW,'(" CALCULATE DERIVATIVES OF EDISP   ...")',
C    *      ADVANCE='NO')
C           CALL CPU_TIME(T0)
C           CALL SYSTEM_CLOCK(IC0,ICR,ICM)
C         END IF
          CALL DFTB_DISP_HESS(NAT,NSPE,ISPE,DFTBDP,ZREF,QREF,C,X(LDIST),
     *    X(LGRAD),X(LHESS),X(LHESSIAN),IDFTBD,X(LWRK))
C         IF (MASWRK.AND.SOME) THEN
C           CALL CPU_TIME(T1)
C           CALL SYSTEM_CLOCK(IC1,ICR,ICM)
C           WRITE (IW,'(" FINISHED IN ",F8.2," (",F8.2,") SECONDS")')
C    *      T1-T0,DBLE(IC1-IC0)/DBLE(ICR)
C         END IF
        END IF
      END IF
C
  100 CONTINUE
C
      ISGDDI=SAVGDDI
C
      IF (EXETYP.EQ.CHECK) THEN
        CALL VCLR(X(LHESSIAN),1,3*NAT*3*NAT)
        CALL VCLR(X(LDIPDER),1,9*NAT)
      END IF
      IF (IPCM.NE.0) THEN
        IMODE=MODPRJ
        IF (MODPRJ.EQ.1.OR.MODPRJ.EQ.2)
     *  CALL DFTB_PRJGRD_HESS(IMODE,EGRAD,X(LHESSIAN),C,RUNTYP,NAT,
     *                        X(LWRK),X(LWRK1))
      END IF
      CALL DAWRIT(IDAF,IODA,X(LHESSIAN),3*NAT*3*NAT,4,0)
      CALL DAWRIT(IDAF,IODA,X(LDIPDER),9*NAT,34,0)
      IF (IDFTBD.EQ.3) CALL DFTD3(3,DUMMY,DUMMY)
      MASWRK = MASSAV
C
      CALL RETFM(NEED)
C
      IF(MASWRK.AND.SOME) THEN
        CHARA(1:11) = '           '
        CHARA(1:4)  = 'DFTB'
        NST = 5
        IF (.NOT.SCC)           CHARA(NST:NST) = '1'
        IF (SCC.AND..NOT.DFTB3) CHARA(NST:NST) = '2'
        IF (SCC.AND.DFTB3)      CHARA(NST:NST) = '3'
        NST = NST+1
        IF (DC) THEN
          IF (IDFTBD.EQ.3) THEN
            CHARA(NST:NST+1) = '-D3'
            NST = NST+2
          ELSE
            CHARA(NST:NST+1) = '-D'
            NST = NST+2
          END IF
        END IF
        IF (IPCM.NE.0) THEN
          CHARA(NST:NST+3) = '/PCM'
          NST = NST+4
        END IF
        WRITE(IW,'(" ...... END OF ",A," SECOND DERIVATIVE ",
     * "......")') CHARA(1:NST-1)
      END IF
C
      GOPARR=SAVGOP
C
      CALL TIMIT(1)
C
      END SUBROUTINE DFTB_HESSIAN
C
C-----------------------------------------------------------------------
C*MODULE DFTBHS    *DECK DFTB_QDER
C>
C>    @brief Hessian of DFTB
C>
C>    @details Caculate Hessian contribution of \Delta q terms
C>
C>    @author Yoshio Nishimoto
C>    - Sep, 2014- Subroutine written
C>    @date Feb, 2016 - Yoshio Nishimoto
C>    - Changes due to PCM and SRSCC
C>
C>           --- INPUT ---
C>    @param NJ Index of atom to be displaced
C>    @param KJ Index of axis to be displaced
C>    @param NAT Number of atoms
C>    @param C Coordinates
C>    @param QDER Delta Q
C>    @param GAMDER Derivative of DFTB2 gamma
C>    @param GAMDER3 Derivative of DFTB3 gamma
C>    @param CHAMUL Mulliken charges
C>    @param DISTMAT Matrix of inter-atomic distances
C>    @param DFTB3 Whether DFTB3 or not
C>    @parma SRSCC Whether Shell-resolved SCC or not
C>    @param NDIMSH If (SRSCC); NDIMSH=NSHELL; else NDIMSH=NAT
C>    @param NSPE Number of species
C>    @param ISPE Index of species
C>    @param MAXANG Maximum angular momentum of each species
C>           --- OUTPUT ---
C>    @param HESSIAN Hessian matrix
C>
C
      SUBROUTINE DFTB_QDER(NJ,KJ,NAT,C,QDER,GAMDER,GAMDER3,CHAMUL,
     * DISTMAT,HESSIAN,DFTB3,SRSCC,NDIMSH,NSPE,ISPE,MAXANG,LPCM,QDERAT,
     * CHAMULAT,NGAMMA)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DOUBLE PRECISION :: QDER(NAT),GAMDER(*),CHAMUL(NDIMSH),DISTMAT(*),
     *  HESSIAN(3*NAT,3*NAT),C(3,*),GAMDER3(NGAMMA,NGAMMA),QDERAT(*),
     *  CHAMULAT(*)
      INTEGER :: NSPE,ISPE(NAT),MAXANG(NSPE)
      LOGICAL :: DFTB3,SRSCC,LPCM
C
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      LOGICAL GOPARR,DSKWRK,MASWRK,PARALL2
C
      DOUBLE PRECISION,PARAMETER :: ONE=1.0D+00,ONE_THIRD=ONE/3.0D+00
      DOUBLE PRECISION :: V(3)
C
      JR = 3*(NJ-1)+KJ
      PARALL2 = GOPARR !.AND.  M.GT.MXSQN2()
      IF (PARALL2.AND..NOT.LPCM) CALL VCLR(HESSIAN(1,JR),1,3*NAT)
      IPCOUNT = ME - 1
C
      IF (SRSCC) THEN
        ISH0 = 0
        DO NI = 1, NAT
          IF(PARALL2) THEN
            IPCOUNT = IPCOUNT + 1
            IF (MOD(IPCOUNT,NPROC).NE.0) THEN
              ISH0 = ISH0 + MAXANG(ISPE(NI))
              CYCLE
            END IF
          END IF
          CHAMULI = CHAMULAT(NI)
          DO I = 1, MAXANG(ISPE(NI))
            ISH = ISH0 + I
            CHAMUL1 = CHAMUL(ISH)
            KSH0 = 0
            DO NK = 1, NAT
              IF (NI.EQ.NK) THEN
                KSH0 = KSH0 + MAXANG(ISPE(NK))
                CYCLE
              END IF
              CHAMULK = CHAMULAT(NK)
              CALL DFTB_CNVSQ(NI,NK,NSEQ)
              DIST = DISTMAT(NSEQ)
              DISTI = ONE/DIST
              DO K = 1, MAXANG(ISPE(NK))
                KSH = KSH0 + K
                CHAMUL2 = CHAMUL(KSH)
                CALL DFTB_CNVSQ(ISH,KSH,NSEQ)
                GAMMAV = GAMDER(NSEQ)
                IF (DFTB3) THEN
                  GAMMAV3A = GAMDER3(ISH,KSH)
                  GAMMAV3B = GAMDER3(KSH,ISH)
                END IF
                DO KI = 1, 3
                  V(KI) = C(KI,NI) - C(KI,NK)
                  V(KI) = V(KI)*DISTI
                  IR = 3*(NI-1)+KI
                  GAMTMP = GAMMAV*V(KI)
                  HESSIAN(IR,JR) = HESSIAN(IR,JR)
     *              + GAMTMP*(QDER(ISH)*CHAMUL2 + CHAMUL1*QDER(KSH))
                  IF (DFTB3) THEN
                    GAMTMP3A = GAMMAV3A*V(KI)
                    GAMTMP3B = GAMMAV3B*V(KI)
                    HESSIAN(IR,JR) = HESSIAN(IR,JR) + ONE_THIRD*
     *          (CHAMULI*GAMTMP3A*(QDER(ISH)*CHAMUL2+CHAMUL1*QDER(KSH))
     *         + CHAMULK*GAMTMP3B*(QDER(ISH)*CHAMUL2+CHAMUL1*QDER(KSH))
     *       +(QDERAT(NI)*GAMTMP3A+QDERAT(NK)*GAMTMP3B)*CHAMUL1*CHAMUL2)
                  END IF
                END DO
              END DO
              KSH0 = KSH0 + MAXANG(ISPE(NK))
            END DO
          END DO
          ISH0 = ISH0 + MAXANG(ISPE(NI))
        END DO
      ELSE
        DO NI = 1, NAT !! A
          IF(PARALL2) THEN
            IPCOUNT = IPCOUNT + 1
C           IF (MOD(IPCOUNT,NPROC).NE.0) CYCLE
            IF (MOD(NI,NPROC).NE.ME) CYCLE
          END IF
          CHAMUL1 = CHAMUL(NI)
          DO NK = 1, NAT !! C
            IF (NI.EQ.NK) CYCLE
            CHAMUL2 = CHAMUL(NK)
            CALL DFTB_CNVSQ(NI,NK,NSEQ)
            DIST = DISTMAT(NSEQ)
            DISTI = ONE/DIST
            GAMMAV = GAMDER(NSEQ)
            IF (DFTB3) THEN
              GAMMAV3A = GAMDER3(NI,NK)
              GAMMAV3B = GAMDER3(NK,NI)
            END IF
            DO KI = 1, 3
              V(KI) = C(KI,NI) - C(KI,NK)
              V(KI) = V(KI)*DISTI
              IR = 3*(NI-1)+KI
              GAMTMP = GAMMAV*V(KI)
              HESSIAN(IR,JR) = HESSIAN(IR,JR) + GAMTMP*(QDER(NI)*CHAMUL2
     +          + CHAMUL1*QDER(NK))
C             write (*,*) "val = ", GAMTMP*(QDER(NI)*CHAMUL2+
C    *        CHAMUL1*QDER(NK))
              IF (DFTB3) THEN
                GAMTMP3A = GAMMAV3A*V(KI)
                GAMTMP3B = GAMMAV3B*V(KI)
                HESSIAN(IR,JR) = HESSIAN(IR,JR)
     *            + ONE_THIRD
     *              *(CHAMUL2*CHAMUL2*QDER(NI)*GAMTMP3B
     *              + CHAMUL1*CHAMUL1*QDER(NK)*GAMTMP3A
     *              + 2.0D+00*CHAMUL1*CHAMUL2*(QDER(NK)*GAMTMP3B
     *                                       + QDER(NI)*GAMTMP3A))
              END IF
            END DO
          END DO
        END DO
      END IF
C
      RETURN
C
      END SUBROUTINE DFTB_QDER
C
C-----------------------------------------------------------------------
C*MODULE DFTBHS    *DECK DFTB_CALC_DERSHIFT
C>
C>    @brief Hessian of DFTB
C>
C>    @details Caculate non-response terms in d(Omega)/dR
C>
C>    @author Yoshio Nishimoto
C>    - Sep, 2014- Subroutine written
C>    @date Feb, 2016 - Yoshio Nishimoto
C>    - Changes due to FMO-DFTB Raman
C>
C>           --- INPUT ---
C>    @param NAT Number of atoms
C>    @param NI Index of atom to be displaced
C>    @param KI Index of axis to be displaced
C>    @param C Coordinates
C>    @param GAMDER Derivative of DFTB2 gamma
C>    @param GAMDER3 Derivative of DFTB3 gamma
C>    @param CHAMUL Mulliken charges
C>    @param DFTB3 Whether DFTB3 or not
C>           --- OUTPUT ---
C>    @param DERSHIFT Shift contribution
C
      SUBROUTINE DFTB_CALC_DERSHIFT(NAT,NI,KI,C,DISTMAT,GAMDER,GAMDER3,
     *  CHAMUL,DERSHIFT,DFTB3)
C
      IMPLICIT NONE
C
      COMMON /EFLDC / EVEC(3),EFLDL
      LOGICAL EFLDL
      DOUBLE PRECISION EVEC
C
      INTEGER, INTENT(IN) :: NAT,NI,KI
      DOUBLE PRECISION, INTENT(IN) :: C(3,*),DISTMAT(*),GAMDER(*),
     *  GAMDER3(NAT,NAT),CHAMUL(NAT)
      DOUBLE PRECISION, INTENT(OUT) :: DERSHIFT(NAT)
      LOGICAL, INTENT(IN) :: DFTB3
C
      INTEGER :: NC,ND,NSEQ
      DOUBLE PRECISION :: V(3),CHAMULI,CHAMULJ,VAL,GAMMAV,GAMMAV3,VAL1,
     *  VAL2,DIST
      DOUBLE PRECISION, PARAMETER :: ONE_THIRD=1.0D+00/3.0D+00,
     *  ONE_SIX=1.0D+00/6.0D+00,TWO=2.0D+00
C
      CALL VCLR(DERSHIFT,1,NAT)
C
C     CALCULATE NON-RESPONSE TERMS IN d(OMEGA)/dR.
C     -NI- CORRESPONDS TO da IN EQUATIONS.
C     -NC- CORRESPONDS TO C IN EQUATIONS.
C     -ND- CORRESPONDS TO C IN EQUATIONS WHEN NI.EQ.NC
C
      DO NC = 1, NAT
        CHAMULI = CHAMUL(NC) !! \Delta q_A
        IF (NI.EQ.NC) THEN
          DO ND = 1, NAT
C           ND = C, NI = NC = a (A)
            IF (NC.EQ.ND) CYCLE
            CHAMULJ = CHAMUL(ND) !! \Delta q_C
            CALL DFTB_CNVSQ(NC,ND,NSEQ)
            V(1) = C(1,NC) - C(1,ND) !! A-C
            V(2) = C(2,NC) - C(2,ND)
            V(3) = C(3,NC) - C(3,ND)
            DIST = DISTMAT(NSEQ)
            GAMMAV = GAMDER(NSEQ)
            VAL = GAMMAV / DIST
            DERSHIFT(NC) = DERSHIFT(NC) + V(KI)*VAL*CHAMULJ
            IF (DFTB3) THEN
              GAMMAV3 = GAMDER3(NC,ND)
              VAL1 = GAMMAV3 / DIST
              GAMMAV3 = GAMDER3(ND,NC)
              VAL2 = GAMMAV3 / DIST
              DERSHIFT(NC) = DERSHIFT(NC)
     *          + ONE_THIRD*VAL1*CHAMULI*CHAMULJ*V(KI)*TWO
     *          + ONE_SIX*VAL2*CHAMULJ*CHAMULJ*V(KI)*TWO
            END IF
          END DO
        ELSE
C         d(gamma)/dR is not zero only when ni=a
C         NC = a (A), NI = C
          CHAMULJ = CHAMUL(NI) !! \Delta q_C
          CALL DFTB_CNVSQ(NC,NI,NSEQ)
          V(1) = C(1,NI) - C(1,NC) !! C-A
          V(2) = C(2,NI) - C(2,NC)
          V(3) = C(3,NI) - C(3,NC)
          DIST = DISTMAT(NSEQ)
          GAMMAV = GAMDER(NSEQ)
          VAL = GAMMAV / DIST
          DERSHIFT(NC) = DERSHIFT(NC) + V(KI)*VAL*CHAMULJ
          IF (DFTB3) THEN
            GAMMAV3 = GAMDER3(NC,NI)
C           GAMMAV3 = GAMDER3(NI,NC)
            VAL1 = GAMMAV3 / DIST
            GAMMAV3 = GAMDER3(NI,NC)
C           GAMMAV3 = GAMDER3(NC,NI)
            VAL2 = GAMMAV3 / DIST
            DERSHIFT(NC) = DERSHIFT(NC)
     *        + ONE_THIRD*VAL1*CHAMULI*CHAMULJ*V(KI)*TWO
     *        + ONE_SIX*VAL2*CHAMULJ*CHAMULJ*V(KI)*TWO
          END IF
        END IF
      END DO
C
      IF (EFLDL) THEN
        DERSHIFT(NI) = DERSHIFT(NI) + EVEC(KI)
      END IF
C
      RETURN
C
      END SUBROUTINE DFTB_CALC_DERSHIFT
C
C-----------------------------------------------------------------------
C*MODULE DFTBHS    *DECK DFTB_CALC_DERSHIFT_SHELL
C>
C>    @brief Hessian of DFTB
C>
C>    @details Caculate non-response terms in d(Omega)/dR for SRSCC
C>
C>    @author Yoshio Nishimoto
C>    - Sep, 2014- Subroutine written
C>    @date Feb, 2016 - Yoshio Nishimoto
C>    - Changes due to SRSCC
C>
C>           --- INPUT ---
C>    @param NAT Number of atoms
C>    @param NI Index of atom to be displaced
C>    @param KI Index of axis to be displaced
C>    @param C Coordinates
C>    @param DISTMAT Matrix of inter-atomic distances
C>    @param GAMDER Derivative of DFTB2 gamma
C>    @param GAMDER3 Derivative of DFTB3 gamma
C>    @param CHAMUL Mulliken charges
C>    @param DFTB3 Whether DFTB3 or not
C>    @param NSHELL Number of shells
C>    @param INDSH Index of shell for each atom
C>    @param ISPE Index of species
C>    @param MAXANG Maximum angular momentum of each species
C>           --- OUTPUT ---
C>    @param DERSHIFT Shift contribution
C>
C>    @todo merged into DFTB_CALC_DERSHIFT
C
      SUBROUTINE DFTB_CALC_DERSHIFT_SHELL(NAT,NI,KI,C,DISTMAT,GAMDER,
     *  GAMDER3,CHAMUL,DERSHIFT,DFTB3,NSHELL,INDSH,NSPE,ISPE,MAXANG,
     *  CHAMULAT)
C
      IMPLICIT NONE
C
      INTEGER, INTENT(IN) :: NAT,NI,KI,NSHELL,INDSH(NAT),NSPE,ISPE(NAT),
     *  MAXANG(NSPE)
      DOUBLE PRECISION, INTENT(IN) :: C(3,*),DISTMAT(*),GAMDER(*),
     *  GAMDER3(NSHELL,NSHELL),CHAMUL(NSHELL),CHAMULAT(NAT)
      DOUBLE PRECISION, INTENT(OUT) :: DERSHIFT(NSHELL)
      LOGICAL, INTENT(IN) :: DFTB3
C
      INTEGER:: ISH0,ISH,JSH0,JSH,KSH0,KSH,I,J,NC,ND,NSEQ,LL,KK,LSH
      DOUBLE PRECISION :: V(3),DIST,DISTI,GAMMAV,CHAMULJ,VAL,
     *  VAL1,VAL2
      DOUBLE PRECISION, PARAMETER :: ONE_THIRD=1.0D+00/3.0D+00,
     *  ONE=1.0D+00
C
      CALL VCLR(DERSHIFT,1,NSHELL)
C
C     CALCULATE SHELL-RESOLVED NON-RESPONSE TERMS IN d(OMEGA)/dR.
C
      ISH0 = INDSH(NI)
      DO I = 1, MAXANG(ISPE(NI)) !! SHELL OF NI
        ISH = ISH0 + I
        JSH0 = 0
        DO NC = 1, NAT
          IF (NI.EQ.NC) THEN
            KSH0 = 0
            DO ND = 1, NAT
C             ND = C, NI = NC = a (A)
              IF (NC.EQ.ND) THEN
                KSH0 = KSH0 + MAXANG(ISPE(ND))
                CYCLE
              END IF
              V(1) = C(1,NC) - C(1,ND) !! A-C
              V(2) = C(2,NC) - C(2,ND)
              V(3) = C(3,NC) - C(3,ND)
              CALL DFTB_CNVSQ(NC,ND,NSEQ)
              DIST = DISTMAT(NSEQ)
              DISTI = ONE/DIST
              DO J = 1, MAXANG(ISPE(ND)) !! SHELL OF ND
                KSH = KSH0 + J
                CALL DFTB_CNVSQ(ISH,KSH,NSEQ)
                GAMMAV = GAMDER(NSEQ)
                CHAMULJ = CHAMUL(KSH) !! \Delta q_C
                VAL = GAMMAV*DISTI
                DERSHIFT(ISH) = DERSHIFT(ISH) + V(KI)*VAL*CHAMULJ
                IF (DFTB3) THEN
                  VAL1 = GAMDER3(ISH,KSH)*DISTI
                  VAL2 = GAMDER3(KSH,ISH)*DISTI
                  DERSHIFT(ISH) = DERSHIFT(ISH) + ONE_THIRD*(
     *                VAL1*CHAMULAT(NC)*CHAMULJ
     *              + VAL2*CHAMULAT(ND)*CHAMULJ)*V(KI)
                  DO LL = 1, MAXANG(ISPE(NI))
                    LSH = ISH0 + LL
                    DERSHIFT(ISH) = DERSHIFT(ISH)
     *                + ONE_THIRD*GAMDER3(LSH,KSH)*DISTI
     *                 *CHAMUL(LSH)*CHAMULJ*V(KI)
                  END DO
                END IF
              END DO
              KSH0 = KSH0 + MAXANG(ISPE(ND))
            END DO
          ELSE
            CHAMULJ = CHAMUL(ISH) !! \Delta q_A
            V(1) = C(1,NI) - C(1,NC) !! C-A
            V(2) = C(2,NI) - C(2,NC)
            V(3) = C(3,NI) - C(3,NC)
            CALL DFTB_CNVSQ(NI,NC,NSEQ)
            DIST = DISTMAT(NSEQ)
            DISTI = ONE/DIST
            DO J = 1, MAXANG(ISPE(NC)) !! SHELL OF NC
              JSH = JSH0 + J
              CALL DFTB_CNVSQ(ISH,JSH,NSEQ)
              GAMMAV = GAMDER(NSEQ)
              VAL = GAMMAV*DISTI
              DERSHIFT(JSH) = DERSHIFT(JSH) + V(KI)*VAL*CHAMULJ
              IF (DFTB3) THEN
                VAL1 = GAMDER3(JSH,ISH)*DISTI
                VAL2 = GAMDER3(ISH,JSH)*DISTI
                DERSHIFT(JSH) = DERSHIFT(JSH) + ONE_THIRD*(
     *              VAL1*CHAMULAT(NC)*CHAMULJ
     *            + VAL2*CHAMULAT(NI)*CHAMULJ)*V(KI)
                DO KK = 1, MAXANG(ISPE(NC))
                  KSH = JSH0 + KK
                  DERSHIFT(JSH) = DERSHIFT(JSH)
     *            + ONE_THIRD*GAMDER3(KSH,ISH)*DISTI
     *             *CHAMUL(KSH)*CHAMULJ*V(KI)
                END DO
              END IF
            END DO
          END IF
          JSH0 = JSH0 + MAXANG(ISPE(NC))
        END DO
      END DO
C
      RETURN
C
      END SUBROUTINE DFTB_CALC_DERSHIFT_SHELL
C
C-----------------------------------------------------------------------
C*MODULE DFTBHS    *DECK DFTB_CP_NR
C>
C>    @brief Hessian of DFTB
C>
C>    @details Caculate non-response terms in CP-DFTB equation
C>
C>    @author Yoshio Nishimoto
C>    - Sep, 2014- Subroutine written
C>    @date Feb, 2016 - Yoshio Nishimoto
C>    - Performance improvement?
C>    @date Feb, 2018 - Yoshio Nishimoto
C>    - For frequency-dependent derivative
C>
C>           --- INPUT ---
C>    @param L1     Number of basis functions
C>    @param L2     L1*(L1+1)/2
C>    @param EIG    Eigenvalues of each orbital
C>    @param OCC    Occupation number of each orbital
C>    @param HDERIJ H_{ij}^a + V_{ij}^a
C>    @param SDERIJ S_{ij}^a
C>    @param NCP    Order of coupled-perturbed
C>    @param FREQ   Flag for frequency-dependent derivative
C>           --- OUTPUT ---
C>    @param U      U matrix (geometrical derivative of MO coefficients)
C>
C
      SUBROUTINE DFTB_CP_NR(L1,L2,NDOCC,NPART,EIG,OCC,U,HDERIJ,SDERIJ,
     *                      EIGINV,MAXVIR,MAXOCC,BETA,EFERMI,DEGTHR,
     *                      NCP,FREQ)
C
      IMPLICIT NONE
C
      INTEGER, INTENT(IN) :: L1,L2,NDOCC,NPART,MAXVIR,MAXOCC,ncp
      DOUBLE PRECISION, INTENT(IN) :: EIG(L1),OCC(L1),HDERIJ(L2),
     *  SDERIJ(L2),EIGINV(MAXVIR,MAXOCC),BETA,EFERMI,DEGTHR
      DOUBLE PRECISION, INTENT(OUT) :: U(L1,L1)
C
      INTEGER :: IO,JO,NSEQ
      DOUBLE PRECISION :: EI,EJ,EIJ,OCCI,OCCJ,VAL,VALI,VALJ,SCAL
      LOGICAL :: PARTI,FREQ
      DOUBLE PRECISION, PARAMETER :: ZERO=0.0D+00,HALF=0.5D+00
C
      NSEQ = 0
      DO IO = 1, L1
        EI = EIG(IO)
        OCCI = OCC(IO)
        PARTI = IO.GT.NDOCC .AND. IO.LE.NDOCC+NPART
        DO JO = 1, IO !! L1
          NSEQ = NSEQ + 1
          IF (IO.LE.NDOCC) CYCLE
          IF (JO.GT.NDOCC+NPART) CYCLE
          IF (IO.LT.JO) CYCLE
          EJ = EIG(JO)
          OCCJ = OCC(JO)
C         IF (OCCI.EQ.OCCJ) CYCLE !! FOR FON
          IF (IO.EQ.JO) THEN
            U(IO,JO) = -HALF*SDERIJ(NSEQ)
            IF (NCP.EQ.2) U(IO,JO) = ZERO
          ELSE
            VAL = HDERIJ(NSEQ)
            VALI = VAL - EJ*SDERIJ(NSEQ)
            VALJ = VAL - EI*SDERIJ(NSEQ)
            IF (NCP.EQ.2) THEN
              VALI = VAL
              VALJ = VAL
            END IF
            SCAL = EIGINV(IO-NDOCC,JO)
            IF (PARTI) THEN
              IF (FREQ) THEN
                U(IO,JO) = VAL*SCAL
                U(JO,IO) = ZERO
              ELSE IF (ABS(EJ-EI).GE.DEGTHR) THEN
                !! AVOID NEARLY DEGENERATED P.O-P.O COUPLING ELEMENTS
                U(IO,JO) = -OCCI*SDERIJ(NSEQ)/OCCJ + VALI*SCAL
                IF (NCP.EQ.2) U(IO,JO) = VALI*SCAL
                U(JO,IO) = ZERO
              ELSE
                !! USE LIMIT OR MACLAURIN EXPANSION
                U(IO,JO) = -OCCI*SDERIJ(NSEQ)/OCCJ + VALI*SCAL
                IF (NCP.EQ.2) U(IO,JO) = VALI*SCAL
                U(JO,IO) = ZERO
              END IF
            ELSE
C             EIJ = ONE/(EJ-EI)
              EIJ = SCAL
              U(IO,JO) = VALI*EIJ
              U(JO,IO) = -VALJ*EIJ !ZERO
            END IF
          END IF
        END DO
      END DO
C
      RETURN
C
      END SUBROUTINE DFTB_CP_NR
C
C-----------------------------------------------------------------------
C*MODULE DFTBHS    *DECK SHIFT_TRIANGLE_DER_CLEAR
C>
C>    @brief Hessian of DFTB
C>
C>    @details VCLR only for values NI-th atom concerns
C>
C>    @author Yoshio Nishimoto
C>    - Sep, 2014- Subroutine written
C>
C>           --- INPUT ---
C>    @param NAT Number of aotms
C>    @param NI Atom to differentiated
C>    @param L2 L1*(L1+1)/2
C>    @param IND Index of AO
C>           --- IN/OUTPUT ---
C>    @param OUTPUT Output (all zero matrix)
C>
C
      SUBROUTINE SHIFT_TRIANGLE_DER_CLEAR(NAT,NI,L2,IND,OUTPUT)
C
      IMPLICIT NONE
C
      INTEGER, INTENT(IN) :: NAT,NI,L2,IND(*)
      DOUBLE PRECISION, INTENT(OUT) :: OUTPUT(L2)
C
      INTEGER :: NJ,M,L,MU,NU,NSEQ
      DOUBLE PRECISION, PARAMETER :: ZERO=0.0D+00
C
      IF (NI.EQ.1) GO TO 100
C
      NSEQ = IND(NI)*(IND(NI)+1)/2
      DO M = 1, IND(NI+1)-IND(NI)
        DO NJ = 1, NI-1
          DO L = 1, IND(NJ+1) - IND(NJ)
            NSEQ = NSEQ + 1
            OUTPUT(NSEQ) = ZERO
          END DO
        END DO
        NSEQ = NSEQ + M
      END DO
C
100   IF (NI.EQ.NAT) RETURN

      NSEQ = IND(NI+1)*(IND(NI+1)+1)/2 + IND(NI)
      DO NJ = NI+1, NAT
        DO L = 1, IND(NJ+1) - IND(NJ)
          NU = IND(NJ)+L
          DO M = 1, IND(NI+1)-IND(NI)
            MU = IND(NI)+M
            NSEQ = NSEQ + 1
            OUTPUT(NSEQ) = ZERO
          END DO
C         MU = IND(NI+1)+1
          NSEQ = NSEQ + IND(NI) + NU - MU
        END DO
      END DO
C
      RETURN
C
      END SUBROUTINE SHIFT_TRIANGLE_DER_CLEAR
C
C-----------------------------------------------------------------------
C*MODULE DFTBHS    *DECK SHIFT_TRIANGLE_DER_SHELL
C>
C>    @brief Hessian of DFTB
C>
C>    @details Multiply shift contribution with a derivative matrix
C>             in AO basis.  Because derivative matrix is sparse
C>             (only columns of NI is non-trivial), this is faster
C>             than simple multiplication.  SHIFT should be
C>             shell-resolved.
C>
C>    @author Yoshio Nishimoto
C>    - Sep, 2014- Subroutine written
C>
C>           --- INPUT ---
C>    @param NAT Number of atoms
C>    @param NI Index of atom to be displaced
C>    @param L2 Size of upper triangular AO matrix
C>    @param IND AO index of each atom
C>    @param INDSH Index of shell for each atom
C>    @param INDSAO AO index of each shell
C>    @param NSPE Number of species
C>    @param ISPE Index of species
C>    @param MAXANG Maximum angular momentum of each species
C>    @param SHIFT Shift contribution (NSHELL size)
C>    @param INPUT Input matrix, should be AO derivative matrix
C>           --- OUTPUT ---
C>    @param OUTPUT Output
C>
C
      SUBROUTINE SHIFT_TRIANGLE_DER_SHELL(NAT,NI,L2,IND,INDSH,INDSAO,
     *  NSPE,ISPE,MAXANG,SHIFT,INPUT,OUTPUT)
C
      IMPLICIT NONE
C
      INTEGER, INTENT(IN) :: NAT,NI,L2,IND(*),INDSH(*),INDSAO(*),NSPE,
     *  ISPE(NAT),MAXANG(NSPE)
      DOUBLE PRECISION, INTENT(IN) :: SHIFT(*),INPUT(L2)
      DOUBLE PRECISION, INTENT(OUT) :: OUTPUT(L2)
C
      INTEGER :: NSEQ,I,J,ISH,JSH,L,M,MU,NU,NJ
      DOUBLE PRECISION :: SHIFTI,SHIFTJ,SHIFTIJ
      DOUBLE PRECISION, PARAMETER :: HALF=0.5D+00
C
      IF (NI.EQ.1) GO TO 200
C
      NSEQ = IND(NI)*(IND(NI)+1)/2
      DO I = 1, MAXANG(ISPE(NI))
        ISH = INDSH(NI)+I
        SHIFTI = SHIFT(ISH)
        DO M = 1, INDSAO(ISH+1)-INDSAO(ISH)
          MU = INDSAO(ISH)+M
          DO NJ = 1, NI-1
            DO J = 1, MAXANG(ISPE(NJ))
              JSH = INDSH(NJ)+J
              SHIFTJ = SHIFT(JSH)
              SHIFTIJ = HALF*(SHIFTI+SHIFTJ)
              DO L = 1, INDSAO(JSH+1) - INDSAO(JSH)
                NU = INDSAO(JSH)+L
                NSEQ = NSEQ + 1
                OUTPUT(NSEQ) = INPUT(NSEQ)*SHIFTIJ
              END DO
            END DO
          END DO
          NSEQ = NSEQ + MU - NU
        END DO
      END DO
C
200   IF (NI.EQ.NAT) RETURN
C
      NSEQ = IND(NI+1)*(IND(NI+1)+1)/2 + IND(NI)
      DO NJ = NI+1, NAT
        DO J = 1, MAXANG(ISPE(NJ))
          JSH = INDSH(NJ)+J
          SHIFTJ = SHIFT(JSH)
          DO L = 1, INDSAO(JSH+1)-INDSAO(JSH)
            NU = INDSAO(JSH)+L
            DO I = 1, MAXANG(ISPE(NI))
              ISH = INDSH(NI)+I
              SHIFTI = SHIFT(ISH)
              SHIFTIJ = HALF*(SHIFTI+SHIFTJ)
              DO M = 1, INDSAO(ISH+1)-INDSAO(ISH)
                MU = INDSAO(ISH)+M
                NSEQ = NSEQ + 1
                OUTPUT(NSEQ) = INPUT(NSEQ)*SHIFTIJ
              END DO
            END DO
            NSEQ = NSEQ + IND(NI) + NU - MU
          END DO
        END DO
      END DO
C
      RETURN
C
      END SUBROUTINE SHIFT_TRIANGLE_DER_SHELL
C
C-----------------------------------------------------------------------
C*MODULE DFTBHS    *DECK SHIFT_TRIANGLE_DER
C>
C>    @brief Hessian of DFTB
C>
C>    @details Multiply shift contribution with a derivative matrix
C>             in AO basis.  Because derivative matrix is sparse
C>             (only columns of NI is non-trivial), this is faster
C>             than simple multiplication.  SHIFT should be
C>             atom-resolved.
C>
C>    @author Yoshio Nishimoto
C>    - Sep, 2014- Subroutine written
C>
C>           --- INPUT ---
C>    @param NAT Number of atoms
C>    @param NI Index of atom to be displaced
C>    @param L2 Size of upper triangular AO matrix
C>    @param IND AO index of each atom
C>    @param SHIFT Shift contribution (NAT size)
C>    @param INPUT Input matrix, should be AO derivative matrix
C>           --- OUTPUT ---
C>    @param OUTPUT Output
C
      SUBROUTINE SHIFT_TRIANGLE_DER(NAT,NI,L2,IND,SHIFT,INPUT,OUTPUT)
C
      IMPLICIT NONE
C
      INTEGER, INTENT(IN) :: NAT,NI,L2,IND(*)
      DOUBLE PRECISION, INTENT(IN) :: SHIFT(*),INPUT(L2)
      DOUBLE PRECISION, INTENT(OUT) :: OUTPUT(L2)
C
      INTEGER :: NSEQ,NJ,L,M,MU,NU
      DOUBLE PRECISION :: SHIFTI,SHIFTJ,SHIFTIJ
      DOUBLE PRECISION, PARAMETER :: HALF=0.5D+00
C
      SHIFTI = SHIFT(NI)
C
      IF (NI.EQ.1) GO TO 200
C
      NSEQ = IND(NI)*(IND(NI)+1)/2
      DO M = 1, IND(NI+1)-IND(NI)
        DO NJ = 1, NI-1
          SHIFTJ = SHIFT(NJ)
          SHIFTIJ = HALF*(SHIFTI+SHIFTJ)
          DO L = 1, IND(NJ+1) - IND(NJ)
            NSEQ = NSEQ + 1
            OUTPUT(NSEQ) = INPUT(NSEQ)*SHIFTIJ
          END DO
        END DO
        NSEQ = NSEQ + M
      END DO
C
200   IF (NI.EQ.NAT) RETURN

      NSEQ = IND(NI+1)*(IND(NI+1)+1)/2 + IND(NI)
      DO NJ = NI+1, NAT
        SHIFTJ = SHIFT(NJ)
        SHIFTIJ = HALF*(SHIFTI+SHIFTJ)
        DO L = 1, IND(NJ+1) - IND(NJ)
          NU = IND(NJ)+L
          DO M = 1, IND(NI+1)-IND(NI)
            MU = IND(NI)+M
            NSEQ = NSEQ + 1
            OUTPUT(NSEQ) = INPUT(NSEQ)*SHIFTIJ
          END DO
          NSEQ = NSEQ + IND(NI) + NU - MU
        END DO
      END DO
C
      RETURN
C
      END SUBROUTINE SHIFT_TRIANGLE_DER
C
C-----------------------------------------------------------------------
C*MODULE DFTBHS    *DECK SHIFT_TRIANGLE_SHELL
C>
C>    @brief Hessian of DFTB
C>
C>    @details Multiply shift contribution with a matrix in AO basis.
C>             SHIFT should be shell-resolved.
C>
C>
C>    @author Yoshio Nishimoto
C>    - Sep, 2014- Subroutine written
C>
C>           --- INPUT ---
C>    @param SHIFT Shift contribution (NSHELL size)
C>    @param INPUT Input matrix, should be AO matrix
C>    @param NAT Number of atoms
C>    @param NSPE Number of species
C>    @param NSHELL Number of shells
C>    @param L2 Size of upper triangular AO matrix
C>    @param ISPE Index of species
C>    @param MAXANG Maximum angular momentum of each species
C>           --- OUTPUT ---
C>    @param OUTPUT Output
C
      SUBROUTINE SHIFT_TRIANGLE_SHELL(SHIFT,INPUT,OUTPUT,
     *  NAT,NSPE,NSHELL,L2,ISPE,MAXANG)
C
      IMPLICIT NONE
C
      INTEGER, INTENT(IN) :: L2,NSHELL
      DOUBLE PRECISION, INTENT(IN) :: INPUT(L2),SHIFT(NSHELL)
      DOUBLE PRECISION, INTENT(OUT) :: OUTPUT(L2)
      INTEGER, INTENT(IN) :: NAT,NSPE,ISPE(NAT),MAXANG(NSPE)
C
      DOUBLE PRECISION :: SHIFTI,SHIFTJ,SHIFTIJ
      DOUBLE PRECISION, PARAMETER :: HALF=0.5D+00
      INTEGER :: NSEQ,MU,NU,I,II,III,J,JJ,JJJ,M,N,ISH,JSH
C
      NSEQ = 0
      MU = 0
      ISH = 0
      DO I = 1, NAT
        DO II = 1, MAXANG(ISPE(I))
          ISH = ISH + 1
          SHIFTI = SHIFT(ISH)
          IF (II.EQ.1) III = 1
          IF (II.EQ.2) III = 3
          IF (II.EQ.3) III = 5
          DO M = 1, III
            MU = MU + 1
            JSH = 0
            NU = 0
            DO J = 1, I
              DO JJ = 1, MAXANG(ISPE(J))
                JSH = JSH + 1
                SHIFTJ = SHIFT(JSH)
                SHIFTIJ = HALF*(SHIFTI+SHIFTJ)
                IF (JJ.EQ.1) JJJ = 1
                IF (JJ.EQ.2) JJJ = 3
                IF (JJ.EQ.3) JJJ = 5
                DO N = 1, JJJ
                  NU = NU + 1
                  IF (NU.GT.MU) EXIT
                  CALL DFTB_CNVSQ(MU,NU,NSEQ)
C                 write (*,*) mu, nu
C                 write (*,*) shifti,shiftj,shiftij
                  OUTPUT(NSEQ) = INPUT(NSEQ)*SHIFTIJ
C                 write (*,*) input(nseq),"->", output(nseq)
                END DO
              END DO
            END DO
          END DO
        END DO
      END DO
C
      RETURN
C
      END SUBROUTINE SHIFT_TRIANGLE_SHELL
C
C-----------------------------------------------------------------------
C*MODULE DFTBHS    *DECK SHIFT_TRIANGLE
C>
C>    @brief Hessian of DFTB
C>
C>    @details Multiply shift contribution with a matrix in AO basis.
C>             SHIFT should be atom-resolved.
C>
C>
C>    @author Yoshio Nishimoto
C>    - Sep, 2014- Subroutine written
C>
C>           --- INPUT ---
C>    @param NAT Number of atoms
C>    @param L2 Size of upper triangular AO matrix
C>    @param IND Index of AOs
C>    @param SHIFT Shift contribution (NSHELL size)
C>    @param INPUT Input matrix, should be AO matrix
C>           --- OUTPUT ---
C>    @param OUTPUT Output
C
      SUBROUTINE SHIFT_TRIANGLE(NAT,L2,IND,SHIFT,INPUT,OUTPUT)
C
      IMPLICIT NONE
C
      INTEGER, INTENT(IN) :: NAT,L2,IND(*)
      DOUBLE PRECISION, INTENT(IN) :: SHIFT(*),INPUT(L2)
      DOUBLE PRECISION, INTENT(OUT) :: OUTPUT(L2)
C
      INTEGER :: NSEQ,MU,NU,L,M,NM,NN
      DOUBLE PRECISION :: SHIFTI,SHIFTJ,SHIFTIJ
      DOUBLE PRECISION, PARAMETER :: HALF=0.5D+00
C
      NSEQ = 0
      NU = 0
      DO NM = 1, NAT
        SHIFTI = SHIFT(NM)
        DO L = 1, IND(NM+1) - IND(NM)
          NU = NU + 1
          MU = 0
          DO NN = 1, NM
            SHIFTJ = SHIFT(NN)
            SHIFTIJ = HALF*(SHIFTI+SHIFTJ)
            DO M = 1, IND(NN+1) - IND(NN)
              MU = MU + 1
              IF (MU.GT.NU) EXIT
              NSEQ = NSEQ + 1
              OUTPUT(NSEQ) = INPUT(NSEQ)*SHIFTIJ
            END DO
          END DO
        END DO
      END DO
C
      RETURN
C
      END SUBROUTINE SHIFT_TRIANGLE
C
C-----------------------------------------------------------------------
C*MODULE DFTBHS    *DECK DFTB_CP
C>
C>    @brief Hessian of DFTB
C>
C>    @details Solve CP-DFTB equation for each vector.
C>
C>    @author Yoshio Nishimoto
C>    - Sep, 2014- Subroutine written
C>    @date Feb, 2016 - Hiroya Nakata
C>    - Added FMO-DFTB Hessian
C>    @date Feb, 2016 - Yoshio Nishimoto
C>    - Performance improvement + DFTB/PCM
C>    @date Feb, 2018 - Yoshio Nishimoto
C>    - Almost rewritten to improve performance
C>    - Third-order derivatives
C>    @date Oct, 2019 - Yoshio Nishimoto
C>    - Added LC-DFTB Hessian
C>    - Bug fix for ROHF-type DFTB Hessian
C>
C>           --- INPUT ---
C>    @param L0      LQMT
C>    @param L1      Number of basis functions
C>    @param L2      L1*(L1+1)/2
C>    @param NDOCC   Number of doubly occupied orbitals
C>    @param NVIRT   Number of virtual orbitals
C>    @param NPART   Number of partially occupied orbitals
C>    @param U       U matrix (derivative of MO coefficients)
C>    @parma VEC     MO coefficients
C>    @param EIG     Eigenvalue of each MO
C>    @param OCC     Occupation number of each MO
C>    @param SHIFT   Shift contribution of each atom
C>    @param GAMMA   Matrix of DFTB2 gamma
C>    @param GAMDER  Derivative matrix of DFTB2 gamma
C>    @param S       Overlap matrix
C>    @param CHAMUL  Mulliken charges of each atom
C>    @param DISTMAT Matrix of inter-atomic distances
C>    @param UINIT   Initial U matrix (only with non-response terms)
C>    @param D       Density matrix
C>    @param HOPAO   HOP contributions in AO (only for FMO)
C>    @param GAMMA3  Matrix of DFTB3 gamma
C>    @param GAMDER3 Derivative matrix of DFTB3 gamma
C>    @param FON     Whether fractional occupation is used or not
C>    @param ROHF    Whether (same-electron) RO-DFTB or not
C>                   RO-DFTB with canonicalization has not been implemen
C>    @param UHF     Whether U-DFTB or not
C>    @param NGAMMA  Dimension of gamma matrices (NAT or NSHELL)
C>    @param SHIFTSA Shift contribution for alpha orbitals (shell-res.)
C>    @param SHIFTSB Shift contribution for beta  orbitals (shell-res.)
C>    @param INDSH   Shell index of each atom
C>    @param INDSAO  AO index of each shell
C>    @param NSPIN   Number of spins (1=R-DFTB,2=U-DFTB)
C>    @param CHAMULS Shell-resolved Mulliken charges
C>    @param NDIMSH  If (SRSCC); NDIMSH=NSHELL; else NDIMSH=NAT
C>    @param HUBBL   Hubbard values of each angular momentum
C>    @param SAVGOP  Saved GOPARR
C>    @param SOME    Print option
C>    @param SPNCST  Spin constants
C>    @param MAXBND  Number of maximum interacting atoms
C>    @param INDBND  Index of interacting atoms
C>    @param DERMAT  AO derivatives of F and S
C>    @param FULL3RD For future
C>    @param CPF     Flag for first-order CP wrt electric field
C>    @param CPA     Flag for first-order CP wrt coordinate
C>    @param CPF2    Flag for second-order CP wrt electric field
C>    @param CPA2    Flag for second-order CP wrt coordinate
C>    @param CPF3    Flag for third-order CP wrt electric field
C>    @param NCPF    Maximum order of CP wrt electric field
C>    @param NSHELL  Number of shells
C>           --- OUTPUT ---
C>    @param HESSIAN Hessian
C>    @param DIPDER  Dipole derivative
C>    @param POLAR   Polarizability
C>    @param ALPDER  Polarizability derivative
C>
C
      SUBROUTINE DFTB_CP(L0,L1,L2,NDOCC,NVIRT,NPART,NOCC,U,
     *  VEC,EIG,OCC,SHIFT,GAMMA,GAMDER,S,CHAMUL,DISTMAT,UINIT,
     *  HESSIAN,D,HOPAO,GAMMA3,GAMDER3,DIPDER,FON,
     *  ROHF,UHF,NGAMMA,SHIFTSA,SHIFTSB,INDSH,INDSAO,
     *  NSPIN,CHAMULS,NDIMSH,HUBB,HUBBL,
     *  SAVGOP,SOME,POLAR,ALPDER,MAXBND,INDBND,DERMAT,FULL3RD,
     *  CPF,CPA,CPF2,CPA2,CPF3,NCPF,NSHELL,GAMLC,GAMDERLC,NSHELL2)
      use mx_limits, only: mxatm
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER (MXSPE=10)
      LOGICAL DFTBFL,SCC,SRSCC,DFTB3,DAMPXH,LCDFTB
      LOGICAL SAVGOP,EFLDL,ARAMAN
C
      DOUBLE PRECISION :: U(L1,L1,NSPIN),VEC(L1,L1,NSPIN),
     *  EIG(L1,NSPIN),OCC(L1,NSPIN),SHIFT(NAT),
     *  UINIT(L1,L1,NSPIN),GAMMA(*),GAMDER(*),S(L2),CHAMUL(NAT),
     *  DISTMAT(*),HESSIAN(3*NAT,3*NAT),HOPAO(L2,NSPIN,3,NAT),
     *  DIPDER(3,3,NAT),D(L2,NSPIN),
     *  GAMMA3(NGAMMA,NGAMMA),GAMDER3(NGAMMA,NGAMMA),
     *  SHIFTSA(NDIMSH),SHIFTSB(NDIMSH),CHAMULS(NDIMSH),
     *  HUBB(*),HUBBL(*),GAMLC(*),GAMDERLC(*)
      DIMENSION INDBND(MAXBND+1,4,NAT),DERMAT(2,*)
      DIMENSION POLAR(6),ALPDER(6,3,NAT)
C
      COMMON /DFTB  / DFTBFL,SCC,SRSCC,DFTB3,DAMPXH,LCDFTB
      COMMON /DFTBHS/ DEGTHR(2),PFREQ(2),CPCONV,OCCTHR,MODHSS,
     *                MODPRJ,MXCPIT,ARAMAN
      COMMON /DFTBMX/ SIMPLE,ITYPMX,TBROYDEN
      COMMON /DFTBPR/ ETEMP,DFTBDP(MXSPE*14),DAMPXHE,HUBDER(MXSPE),
     *                ZREF(MXATM),SPNCST(6,MXSPE),SPE(MXATM),NSPE,
     *                MAXANG(MXATM),ISPE(MXATM),IND(MXATM+1),IDFTBD,
     *                PARAMDIR
      COMMON /EFLDC / EVEC(3),EFLDL
      COMMON /FMCOM / X(1)
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON /MACHIN/ NWDVAR,MAXFM,MAXSM,LIMFM,LIMSM
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      LOGICAL GOPARR,DSKWRK,MASWRK,PARALL3
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
      COMMON /GDDI/   ISCOPE,NGROUPS,MYGROUP,MEGLOB,npglob,nnglob,JBTYP,
     *                ISGDDI,PAROUT,INITGDDI,wasgddi,MLGDDI,NSUBGR,
     *                MeUniv,NPUniv,numdlb,myworld,nworlds,mogddi
C     for FMO-DFTB Hessian
      common /fmoinf/ nfg,nlayer,natfmo,nbdfg,naotyp,nbody,nsegm
      common /fmorun/ espscf,e0scf(2),emp2s,IDAFMO,icurfg,jcurfg,kcurfg,
     *                icurlay,icurunt,nat1e,ncursh,ngau,icurpop,ifmostp,
     *                moncor,needr,modrst,norbproj,nunesp,iskipesp,
     *                IESDPPC,idoprop,mp2run,icurit,idmfmo,iddfmo,
     *                iddcur,nddleft,ivmfmo,nzmtfmo,ifmobas,itmfmo(2)
      COMMON /SCZLAG/ LZVLAG,LZVWRK,LYALAG,LYAWRK,LFEQ1
C
      LOGICAL MYJOB,ISGDDI,PAROUT,INITGDDI,WASGDDI,MLGDDI
      COMMON /RUNOPT/ RUNTYP,EXETYP,NEVALS,NGLEVL,NHLEVL
      DATA  CHECK/8HCHECK   /
C
      INTEGER :: O,INDSH(*),INDSAO(*),NDOCC(2),NPART(2),NVIRT(2),
     *  NOCC(2),NDEG2(3),NDEG3(3),INDPAR(NCPF)
      INTEGER, PARAMETER :: DDI_WORLD=0,DDI_GROUP=1
      DOUBLE PRECISION,PARAMETER ::
     *  ZERO=0.0D+00,ONE=1.0D+00,TWO=2.0D+00,
     *  DEBYE_ANGSTROM=2.541746675D+00*1.889725989D+00,
     *  BOLTZ=3.1668114D-06,ONE_THIRD=1.0D+00/3.0D+00
      DOUBLE PRECISION :: V(3),AUP(9,9),BUP(9,9),FEDERTMP(4)
C
      LOGICAL :: FON(3),ROHF,UHF,SKP,SOME,SAVDCQ,FULL3RD,CPF,CPA,CPF2,
     *           CPA2,CPF3,FREQ !! first
      CHARACTER(8) :: CHARA
C
      PARALL3 = .FALSE.
      IF (MASWRK.AND.SOME) THEN
        WRITE (IW,*)
        WRITE (IW,*) "       *** START COUPLED-PERTURBED DFTB ***"
        WRITE (IW,*)
      END IF
C     IF (ISGDDI) CALL GDDI_SCOPE(DDI_GROUP)
C     IF (ISGDDI) CALL GDDI_SCOPE(DDI_WORLD)
      FREQ = .FALSE.
      DO I = 1, 2 !! NCPF
        IF (ABS(PFREQ(I)).NE.ZERO) THEN
          FREQ = .TRUE.
          EXIT
        END IF
      END DO
      NFREQ = I
      IF (I.EQ.3 .OR. .NOT.CPF) NFREQ = 0
      if (freq .or. cpf3) then
        if (maswrk) write (iw,*) "PFREQ is not yet implemented"
        call abrt
      end if
      L3 = L1*L1
      CALL VALFM(LOADFM)
      CALL GOTFM(NGOTMX)
C
      LDEDB        = LOADFM     + 1
      LSDER        = LDEDB      + L1*NSPIN
      LHDERIJ      = LSDER      + L2
      LSDERIJ      = LHDERIJ    + L2*NSPIN
      LSOMGDIJ     = LSDERIJ    + L2*NSPIN
      LQDER_INIT   = LSOMGDIJ   + L2*NSPIN
      LDOMEGA      = LQDER_INIT + NDIMSH*NSPIN
      LWRK         = LDOMEGA    + L2*NSPIN
      LWRK1        = LWRK       + MAX(L3,3*NAT)
      LWRK2        = LWRK1      + L3
      LWRK3        = LWRK2      + MAX(L3,L2*NSPIN)
      LQDER        = LWRK3      + L3
      IF (FON(3)) THEN
      LDDER        = LQDER      + NDIMSH*NSPIN+1
      ELSE
      LDDER        = LQDER      + NDIMSH*NSPIN
      END IF
      IF (FREQ) THEN
      LWDER        = LDDER      + L3*NSPIN
      ELSE
      LWDER        = LDDER      + L2*NSPIN
      END IF
      LDTMP        = LWDER      + L2*NSPIN
      LEIGINV      = LDTMP      + L2*NSPIN
        MAXOCC     = MAX(NOCC(1),NOCC(2))
        MAXVIR     = MAX(NPART(1)+NVIRT(1),NPART(2)+NVIRT(2))
      LAST         = LEIGINV    + MAX(MAXOCC*MAXVIR*NSPIN,L1+1) !! L1+1 IS FOR SPARSE
      LFDERIJ = LWRK !! may be not needed actually
C
      IF (SCC) THEN
        LDERSHIFT  = LAST
        LDSH       = LDERSHIFT  + NDIMSH
        LQDER_PREV = LDSH       + NDIMSH
        IF (FON(3)) THEN
        LSDEROMG   = LQDER_PREV + NDIMSH*NSPIN+1
        ELSE
        LSDEROMG   = LQDER_PREV + NDIMSH*NSPIN
        END IF
        LAST       = LSDEROMG   + L2
      END IF
C
      IF (FON(3)) THEN
        !! UHF has not been implemented...
        LDEDB_INIT = LAST
        LDNDR      = LDEDB_INIT + L1
        LAST       = LDNDR      + L1
      ELSE
        LDEDB_INIT = LAST
        LDNDR      = LAST
      END IF
C
      IF (UHF) THEN
        LSHIFTSP   = LAST
        LAST       = LSHIFTSP   + NDIMSH
      ELSE
        LSHIFTSP   = LAST
      END IF
C
      IF (DFTB3.AND.NFG.NE.0) THEN
        LSHIFT3    = LAST
        LAST       = LSHIFT3    + NAT*2
      END IF
C
      MEMPCM = 0
      IF (IPCM.NE.0.OR.NFMOPCM.GT.0) THEN
        NDASC = 0
        LFPNT = LAST
        if (isgddi) then
          !! Silly strategy for ASCDII
          !! If Broyden is used, no need to do such
          NTSPAR0=(NTS-1)/npglob+1
          NTSPAR =(NTS-1)/nproc+1
        else
          NTSPAR0=(NTS-1)/nproc+1
          NTSPAR =(NTS-1)/nproc+1
        end if
        CALL INDQPROD(0,NTSPAR,LFPNT,NDASC,DUM,DUM,DUM,0,0)
        LAST    = LFPNT   + NDASC  + 1
C
        LEPOT   = LAST
        LPCMTMP = LEPOT      + NTS
        LQSEDER = LPCMTMP    + NTS
        LQTMP   = LQSEDER    + NTS
        LAST    = LQTMP      + NTS
C
C       MUST BE CHECKED FOR FMO
        LFIXPV = 0
        LFIXPV2 = 0
        IF(NFMOPCM.EQ.0) THEN
          IF(IPCDER.EQ.3) LFIXPV = 3*20*MXTS
          IF(IPCDER.EQ.3) LFIXPV2= 21*MXTS/NWDVAR + 1
          LDAI    = LAST
          LIDDAI  = LDAI        + LFIXPV
          LAST    = LIDDAI      + LFIXPV2
        ENDIF
        MEMPCM  = NDASC+NTS*6+LFIXPV+LFIXPV2
      END IF
C
      NVECF = 0
      NVECF2= 0
      IF (CPF.OR.CPF2) THEN
C       NCPF=1: NN=3
C       NCPF=2: NN=9
C       NCPF=3: NN=19
        NN = 0
        IVEC = 1
        DO ICP = 1, NCPF
          IF (FREQ) THEN
            NN = NN + 3**ICP
          ELSE
            IVEC = IVEC + (ICP+1)
            NN = NN + IVEC
          END IF
        END DO
        NVECF = NN
        !! 3**NCPF * (1+2*(NCPF-1))
C       IF (NFREQ.GT.0) THEN
        NFREQ1 = 0
        IF (NFREQ.GT.0) THEN
          IF (IAND(MODHSS,2+32).NE.0) THEN
            DO ICP = 1, NCPF
              NN = NN + 3**ICP * (1+2*(ICP-1))
              NFREQ1 = NFREQ1 + 2+2*(ICP-1)
            END DO
          END IF
          IF (IAND(MODHSS,64).NE.0) THEN
            !! FO: (-2w)
            !! SO: (-2w,w)
            !! TO: (-2w,w,w)
            NN = NN + 3 + 9 + 27
            NFREQ1 = NFREQ1 + 3
          END IF
        END IF
C                     nfreq1 = nfreq1 + 1
        NVECF2 = NN
        LEDERF  = LAST
        LQDERF  = LEDERF  +     L1*NSPIN*3
        LFONDER = LQDERF  + NDIMSH*NSPIN*NN !! Q_A^F
        LFDERF  = LFONDER +     L1*NSPIN*3  !! derivative of FON
        LDDERF  = LFDERF  +     L2*NSPIN*NN !! F_{ij}^F
        IF (FREQ) THEN
        LEIGINV2= LDDERF  +     L3*NSPIN*NN !! D_{mn}^F
        ELSE
        LEIGINV2= LDDERF  +     L2*NSPIN*NN !! D_{mn}^F
        END IF
        LFDERFSC= LEIGINV2+     L3*NSPIN
        LFEDER  = LFDERFSC+     L3*NSPIN*3 
        LAST    = LFEDER  +        NSPIN*NN
        IF (FREQ) THEN
        LINDFRQ = LAST
        LAST    = LINDFRQ + (NCPF+1)*NVECF2/NWDVAR+1
        END IF
C
        LFDERSC = LSOMGDIJ !! NEED L3*NSPIN (< 2*L2*NSPIN)
        LSDERSC = LDDER
      ELSE
        LEDERF  = LAST
        LQDERF  = LAST
        LFONDER = LAST
        LDDERF  = LAST
        LEIGINV2= LAST
        IF (CPA2) LAST = LEIGINV2 + L3*NSPIN
        LFEDER  = LAST
C
        LFDERFSC= LWRK
        LFDERF  = LWRK2
      END IF
C
      IF (IAND(MODHSS,8+16+128).NE.0) THEN
        NVEC = 3*NAT
        IF (CPA2) NVEC = NVEC + 3*NAT*(3*NAT+1)/2
        LDSAVE  = LAST
        LFSAVE  = LDSAVE  + L2*NSPIN*NVEC
        LSSAVE  = LFSAVE  + L2*NSPIN*NVEC
        LFESAVE = LSSAVE  + L2*NSPIN*NVEC
        LAST    = LFESAVE + NSPIN*NVEC
      END IF
C
      IF (CPA2) THEN
        LGAMHES = LAST
        LAST    = LGAMHES + NGAMMA*(NGAMMA+1)/2
        IF (DFTB3) THEN
          LGAMHES3 = LAST
          LAST     = LGAMHES3 + NGAMMA*NGAMMA
        ELSE
          LGAMHES3 = LAST
        END IF
      END IF
C
      IF (CPF) THEN
        NFDER2 = 0
        NFDER3 = 0
        NFDER4 = 0
        IF (MASWRK.AND.SOME .AND. IAND(MODHSS,2+32).NE.0) THEN
          IF (PFREQ(1).EQ.ZERO) THEN
            IF (NCPF.GE.1) NFDER2 = 9
            IF (NCPF.GE.2) NFDER3 = 27
            IF (NCPF.GE.3) NFDER4 = 81
          ELSE
            IF (NCPF.GE.1) NFDER2 = 9*2
            IF (NCPF.GE.2) NFDER3 = 27*4
            IF (NCPF.GE.3) NFDER4 = 81*6
          END IF
        END IF
        LFDER2 = LAST
        LFDER3 = LFDER2 + NFDER2
        LFDER4 = LFDER3 + NFDER3
        LAST   = LFDER4 + NFDER4
      END IF
C
      IF (LCDFTB) THEN
        LDELD  = LAST
        LFAOLC = LDELD  + L2*NSPIN
        LFAOLC2= LFAOLC + L2*NSPIN
        LAST   = LFAOLC2+ L2*NSPIN
      END IF
C
      IF (ITYPMX.LT.0 .OR. LCDFTB) THEN
        LDDER_PREV = LAST
        LAST       = LDDER_PREV + L2*NSPIN
      END IF
      NEED         = LAST - LOADFM - 1
C
C     ----- GET MEMORY -----
C
      CALL GETFM(NEED)
      IF (EXETYP.EQ.CHECK) GO TO 1000
      IF (UHF.OR.SRSCC) NSHELL = NDIMSH
      IF (MASWRK.AND.SOME) THEN
        WRITE (IW,'(" ...... ALLOCATED MEMORY IN DFTB_CP",
     * " = ",I10," WORDS")') NEED - MEMPCM
      END IF
C
      CALL GOTFM(NGOTMX)
C
      IF (EXETYP.EQ.CHECK) GO TO 1000
C
      CALL VCLR(X(LDEDB),1,L1*NSPIN)
      IF (FON(3)) THEN
        CALL VCLR(X(LDEDB_INIT),1,L1)
        CALL VCLR(X(LDNDR),1,L1)
        IF (SCC) THEN
          CALL DFTB_ETEMP(L0,NE,EIG,OCC,ETEMP,EMERMIN,EFERMI,OCCTHR)
          IF (ETEMP.EQ.0.0D+00) THEN
            BETA=ONE
          ELSE
            BETA=ONE/(ETEMP*BOLTZ)
          END IF
          CALL VCLR(FEDERTMP,1,4)
        END IF
      END IF
C
      IF (FON(3)) THEN
        !! CALCULATE PARTIAL DENSITY MATRIX
        !! D'_{\mu \nu} = \sum_i^{NDOCC} n_i c_{\mu i} c_{\nu i}
        IF (NDOCC(1).GE.1) THEN
          DO ISPIN = 1, NSPIN
            CALL DFTB_DPART(L1,L2,NDOCC(ISPIN),OCC(1,ISPIN),
     *        X(LDTMP+L2*(ISPIN-1)),VEC(1,1,ISPIN))
          END DO
        ELSE
          CALL VCLR(X(LDTMP),1,L2*NSPIN)
        END IF
      ELSE
        CALL DCOPY(L2*NSPIN,D,1,X(LDTMP),1)
      END IF
C
C     ONE/(EI-EJ)
C
      CALL DFTB_EIGINV(1,L1,NSPIN,NDOCC,NPART,MAXVIR,MAXOCC,X(LEIGINV),
     *                 EIG,OCC,DEGTHR(1),BETA,FREQ,PFREQ(1))
C
      IF (CPA2) THEN
        IF (SRSCC) THEN
          NN = 0
          DO I = 1, NSPE
            NN = NN + MAXANG(I)
          END DO
        ELSE
          NN = NSPE
        END IF
        CALL DFTB_CALC_GAM_DER2(2,X(LWRK),X(LGAMHES),X(LAST),X(LWRK),
     *                          X(LGAMHES3),X(LAST),DISTMAT,HUBBL,
     *                          HUBDER,DAMPXHE,ISPE,MAXANG,NAT,NSPE,
     *                          NSHELL,NGAMMA,NN,SPE,DFTB3,SRSCC,DAMPXH)
      END IF
C
C     ----- RESTORE ESP CONTRIBUTION, ONLY WHEN DFTB3 -----
C
      IF (DFTB3.AND.NFG.NE.0) THEN
        GOPARR = SAVGOP
        CALL DAREAD(IDAF,IODA,X(LSHIFT3),NAT*2,567,0)
        IF (GOPARR) GOPARR = .FALSE.
      END IF
C
      MXITER = MXCPIT
      ALPHA = 1.0D+00
C     IF (SRSCC.OR.UHF) THRES = 1.0D-08
      THRES = CPCONV
      NITER = 0
      NITER2= 0
      NEEDAR = 0
      IF ((FON(3).AND.SOME.AND.MASWRK).OR.ARAMAN.OR.CPF.OR.CPA2) THEN
        CALL DFTB_HESS_CHKDEG(NSPIN,L1,EIG,OCC,DEGTHR(2),NDOCC,
     *                        NPART,IW,NDEG2,NDEG3,SOME.AND.MASWRK,
     *                        X(LFDERFSC),X(LFDERF))
C         write (*,*) ndeg2(1),ndeg2(2),ndeg2(3)
C         write (*,*) ndeg3(1),ndeg3(2),ndeg3(3)
        IF (ARAMAN.OR.CPF.OR.CPA2) THEN
          IF (NDEG2(3)*2.GT.L3*NSPIN*3.OR.NDEG3(3)*3.GT.L2*NSPIN*6) THEN
            if (maswrk) write (iw,*) "check chkdeg"
            if (maswrk) write (iw,*) "ndeg2(3)*2 = ", ndeg2(3)*2
            if (maswrk) write (iw,*) "ndeg3(3)*3 = ", ndeg3(3)*3
            call abrt
          END IF
          CALL VALFM(LOADFM)
          LLSTDEG2 = LOADFM + 1
          LLSTDEG3 = LLSTDEG2 + NDEG2(3)*2
          LAST     = LLSTDEG3 + NDEG3(3)*3
          NEEDAR   = LAST - LOADFM - 1
          CALL GETFM(NEEDAR)
C
          CALL DCOPY(NDEG2(3)*2,X(LFDERFSC),1,X(LLSTDEG2),1)
          CALL DCOPY(NDEG3(3)*3,X(LFDERF),1,X(LLSTDEG3),1)
        END IF
      END IF
C
      IF (NEEDAR.EQ.0) THEN
        NDEG2(1) = 0
        NDEG2(2) = 0
        NDEG2(3) = 0
        NDEG3(1) = 0
        NDEG3(2) = 0
        NDEG3(3) = 0
        LLSTDEG2 = LAST
        LLSTDEG3 = LAST
      END IF
C     PREPARE PCM INFORMATION (IN WORLD SCOPE)
C
      SAVDCQ = .FALSE.
      IF (IPCM.NE.0.OR.NFMOPCM.GT.0) THEN
        CALL DFTB_QPOTEN(NAT,NTS,C,X(LAXYZCT),X(LAXYZCT+MXTS),
     *                   X(LAXYZCT+MXTS*2),CHAMUL,X(LEPOT))
C       CALL INDQPROD(1,NTSPAR,LFPNT,NDASC,X(LEPOT),X(LQSE),X(LQSE),1,0)
       CALL INDQPROD(1,NTSPAR0,LFPNT,NDASC,X(LEPOT),X(LQSE),X(LQSE),1,0)
C       call dscal(nts,-one,x(lepot),1)
C
        IF(IPCDER.EQ.3.AND.NFMOPCM.EQ.0) THEN
          CALL DAREAD(IDAF,IODA,X(LDAI),LFIXPV,334,0)
          CALL DAREAD(IDAF,IODA,X(LIDDAI),LFIXPV2,335,1)
        ENDIF
C
        !! Store dC/da*q matrix, if memory is available
        CALL GOTFM(NGOTMX)
        CALL VALFM(LOADFM)
        LDCDAQ = LOADFM + 1
        LAST   = LDCDAQ + NTS*NAT*3
        NEED3  = LAST - LOADFM - 1
        IF (NEED3.LE.NGOTMX) THEN
          CALL GETFM(NEED3)
          SAVDCQ = .TRUE.
          CALL VCLR(X(LDCDAQ),1,NTS*NAT*3)
          IPCOUNT = ME - 1
          DO NI = 1, NAT
            IF(GOPARR) THEN
              IPCOUNT = IPCOUNT + 1
              IF (MOD(IPCOUNT,NPROC).NE.0) CYCLE
            END IF
            DO KI = 1, 3
              CALL DFTB_CALC_DCDA(NTS,NI,KI,X(LAXYZCT),
     *             X(LAXYZCT+MXTS),X(LAXYZCT+MXTS*2),
     *             X(LAXYZCT+MXTS*3),X(LISPHE),
     *             X(LEPSPCM),X(LIDDAI),X(LDAI),X(LQSE),
     *             X(LDCDAQ+3*NTS*(NI-1)+NTS*(KI-1)))
            END DO
          END DO
          IF(GOPARR) CALL DDI_GSUMF(2418,X(LDCDAQ),NTS*NAT*3)
        ELSE
          NEED3 = 0
        END IF
        IF (MASWRK.AND.SOME) THEN
          IF (SAVDCQ)
     *    WRITE (IW,'(" ENOUGH MEMORY TO STORE DC/DA*Q FOR PCM",
     *     " ( ",I10," WORDS)")') NEED3
          WRITE (IW,'(" ...... ALLOCATED MEMORY FOR PCM = ",I10,
     *      " WORDS")') MEMPCM + NEED3
        END IF
      END IF
C
      IF (LCDFTB) THEN
        CALL VCLR(X(LWRK),1,L2)
        CALL DFTB_MAKEREFD(NAT,X(LWRK))
        CALL DCOPY(L2,D,1,X(LDELD),1)
        CALL DAXPY(L2,-ONE/DBLE(NSPIN),X(LWRK),1,X(LDELD),1)
        IF (NSPIN.EQ.2) THEN
          CALL DCOPY(L2,D(1,2),1,X(LDELD+L2),1)
          CALL DAXPY(L2,-ONE/DBLE(NSPIN),X(LWRK),1,X(LDELD+L2),1)
        END IF
      END IF
C
C
C     ----- CONSTRUCT X(LINDFRQ); CONTAINS THE INFORMATION OF W -----
C
C     IF (FREQ) THEN
C       CALL DFTB_CNSTINDFRQ(NCPF,MODHSS,NFREQ1,X(LINDFRQ))
C     END IF
      IF (ISGDDI) CALL GDDI_SCOPE(DDI_GROUP)
      PARALL3 = GOPARR
C     IPCOUNT = ME - 1
      IF (PARALL3) THEN !! .AND..NOT.ISGDDI) THEN
C       GOPARR = .FALSE. !! AVOID PARALLEL CALCULATION IN TFTRI0
        IF (MASWRK.AND.SOME) THEN
          WRITE (IW,'(" PARALLEL CALCULATION WITH ",I3," CORES")')NPROC
        ELSE
C         DSKWRK = .FALSE.
        END IF
      END IF
C
      IF (ISGDDI) THEN
        IF (MASWRK.AND.SOME) THEN
          WRITE (IW,'(" PARALLEL (GDDI) CALCULATION WITH ",
     *      I3," GROUPS")') NGROUPS
        END IF
        CALL GDDI_SCOPE(DDI_GROUP)
c       CALL DDI_COMMID(DDI_WORLD,IDWORLD)
        CALL GDDICOUNT(-1,LGROUP,MYJOB)
      END IF
C
      NI0 = 1
      IF (ARAMAN.OR.CPF.OR.CPF2.OR.CPA2) THEN
        IF (CPF.OR.CPF2) THEN
          NI0 = 0
          CALL DFTB_DIPINT(IND,ZERO,ZERO,ZERO,S,.FALSE.)
        END IF
        CALL DFTB_EIGINV(2,L1,NSPIN,NDOCC,NPART,L1,L1,X(LEIGINV2),
     *                   EIG,OCC,DEGTHR(1),BETA,.FALSE.,ZERO)
C       TAR  = ZERO
C       TARW = ZERO
C       CALL CPU_TIME(TAR0)
C       CALL SYSTEM_CLOCK(ICAR0,ICR,ICM)
      END IF
                NAT0 = 0
      IF (CPA)  NAT0 = NAT
      IF (CPA2) NAT0 = NAT + NAT*(NAT+1)/2
      IFREQ2 = 0
C     FIRST = .TRUE.
      DO NI = NI0, NAT0
        NCP = 1
        NKI = 3 ! x,y,z
C       IF (NI.EQ.0.AND.CPF2) NKI = 9 ! x,y,z,xx,xy,...
C       IF (NI.EQ.0.AND.CPF3) NKI = 19
        IF (NI.EQ.0.AND.NCPF.GT.1) NKI = NVECF
        IF (NI.GT.NAT) THEN
          NCP = 2
C
          NI1 = NI
          NI2 = 1
          DO
            NI1 = NI1 - (NAT-NI2+1)
            IF (NI1.LE.NAT) EXIT
            NI2 = NI2 + 1
          END DO
          IF (NI1.EQ.NI2) THEN
            NKI = 6
          ELSE
            NKI = 9
          END IF
        END IF
        DO KI = 1, NKI
        NFREQ2 = 1
        NFREQ3 = 0
        IF (NI.EQ.0 .AND. NFREQ.EQ.1) THEN
          NTMP = 0
          DO ICPF = 1, NCPF
            IF (KI.LE.NTMP+3**ICPF) EXIT
            NTMP = NTMP + 3**ICPF
            NFREQ3 = NFREQ3 + 2*ICPF
          END DO
          NCP = ICPF
          NFREQ2 = 2*NCP
          IF (IAND(MODHSS,64).NE.0) THEN
            IF (NCP.LE.3) THEN
              NFREQ2 = NFREQ2 + 1
              NFREQ3 = NFREQ3 + NCP-1
            ELSE
              NFREQ3 = NFREQ3 + 3
            END IF
          END IF
        END IF
        DO IFREQ = 1, NFREQ2
          IFREQ2 = IFREQ2 + 1
          IF (ISGDDI.AND.NI.GE.1) THEN
            call GDDICOUNT(0,LGROUP,MYJOB)
            IF (.NOT.MYJOB) CYCLE
          END IF
          IF (NI.EQ.0) THEN
            IF (MASWRK.AND.SOME) THEN
              IF (KI.EQ.1 .AND. IFREQ.EQ.1) WRITE (IW,*)
              IF (IFREQ.EQ.1)
     *        WRITE (IW,'(" CP-DFTB WITH ELECTRIC FIELD")',ADVANCE='NO')
            END IF
            NTMP = 0
            IVEC = 1
            DO ICPF = 1, NCPF
C             IF (FREQ) THEN
C               IF (KI.LE.NTMP+3**ICPF) EXIT
C               NTMP = NTMP + 3**ICPF
C             ELSE
                IVEC = IVEC + (ICPF+1)
                IF (KI.LE.NTMP+IVEC) EXIT
                NTMP = NTMP + IVEC
C             END IF
            END DO
            IF (KI.GT.NTMP) THEN
C             KITMP = KI-NTMP
C             ICPFTMP = ICPF
C             CALL VICLR(INDPAR,1,NCPF)
C             CALL DFTB_SETKI(1,KITMP,ICPFTMP,1,3,INDPAR,FREQ)
                indpar(1) = ki
              CHARA(1:8) = '        '
              DO I = 1, NCPF
                J = INDPAR(I)
                IF (J.EQ.1) CHARA(I:I) = 'X'
                IF (J.EQ.2) CHARA(I:I) = 'Y'
                IF (J.EQ.3) CHARA(I:I) = 'Z'
              END DO
              IF (IFREQ.EQ.1 .AND. MASWRK .AND. SOME)
     *        WRITE (IW,'(" ALONG ",A8," AXIS:")',ADVANCE='NO') CHARA
            END IF
            IF (FREQ .AND. KI.GE.4 .AND. KI.LE.12) THEN
              KI1 = INDPAR(1)
              KI2 = INDPAR(2)
              NCP = 2
            ELSE IF (KI.GE.4.AND.KI.LE.9) THEN
              CALL DFTB_KI12(KI-3,KI1,KI2)
              NCP = 2
            END IF
            NCP = ICPF
            IF (NCP.EQ.1) THEN
              CALL DAREAD(IDAF,IODA,X(LWRK),L2,94+KI,0)
              CALL DFTB_CPYTSQ(X(LWRK),X(LWRK1),L1,L1)
              DO ISPIN = 1, NSPIN
                CALL DGEMM('T','N',L1,L1,L1,ONE,
     *                     VEC(1,1,ISPIN),L1,X(LWRK1),L1,
     *                     ZERO,X(LWRK2),L1)
                CALL DGEMM('N','N',L1,L1,L1,ONE,
     *                     X(LWRK2),L1,VEC(1,1,ISPIN),L1,ZERO,
     *                     X(LWRK),L1)
                CALL CPYSQT(X(LWRK),X(LHDERIJ+L2*(ISPIN-1)),L1,1)
              END DO
            ELSE
              NI1 = 0
              NI2 = 0
              CALL VCLR(X(LHDERIJ),1,L2*NSPIN)
            END IF
            IF (FREQ) THEN
C             CALL DFTB_CNVOMG(NCP,PFREQ,IFREQ,INDOMG,FRQOMG)
C             CALL ICOPY(NCPF,X(LINDFRQ+(NCPF+1)*(IFREQ+NFREQ3-1)+1),1,
C    *                   INDOMG,1)
C             FRQOMG = PFREQ(1)*DBLE(SUM(INDOMG,NCPF))
C             CALL DFTB_EIGINV(1,L1,NSPIN,NDOCC,NPART,MAXVIR,MAXOCC,
C    *          X(LEIGINV),EIG,OCC,DEGTHR(1),BETA,FREQ,FRQOMG)
            END IF
            CALL VCLR(X(LSDERIJ),1,L2*NSPIN)
          ELSE
            IF (NCP.EQ.1) THEN
              NI1 = NI
              KI1 = KI
              NI2 = 0
              KI2 = 0
            ELSE IF (NCP.EQ.2) THEN
              IF (NI1.EQ.NI2) THEN
                CALL DFTB_KI12(KI,KI1,KI2)
              ELSE
                KI1 = (KI-1)/3+1
                KI2 = MOD(KI,3)
                IF (KI2.EQ.0) KI2 = 3
              END IF
            END IF
            CALL DFTB_HESS_DIRECT(X(LHDERIJ),X(LSDERIJ),X(LSDER),
     *        HOPAO(1,1,KI,NI),C,VEC,SHIFT,SHIFTSA,SHIFTSB,X(LWRK),
     *        X(LWRK1),X(LWRK2),X(LWRK3),NI1,NI2,KI1,KI2,NAT,L0,L1,L2,
     *        NDOCC,NPART,IND,NSHELL,NSPE,NSPIN,INDSH,INDSAO,ISPE,
     *        MAXANG,SCC,SRSCC,UHF,ARAMAN.OR.CPA2.OR.FULL3RD,
     *        NFG.NE.0.AND.NBDFG.NE.0,INDBND,MAXBND)
          END IF
C
C         CALCULATE INITIAL DENSITY MATRIX DERIVATIVE FOR SO-CP
C
          FEDER = ZERO
C         IF (NCP.GE.2) THEN
C           CALL DFTB_DDERINIT(NI1,NI2,KI1,KI2,L1,L2,NAT,NSPIN,NDIMSH,
C    *           NDOCC,NOCC,IND,INDSAO,ISPE,MAXANG,
C    *           NDEG2,NDEG3,X(LLSTDEG2),X(LLSTDEG3),
C    *           OCC,EIG,VEC,BETA,DEGTHR(2),
C    *           X(LSDERIJ),X(LFDERF),X(LFSAVE),X(LSSAVE),FEDERTMP,
C    *           X(LFESAVE),X(LEIGINV2),
C    *           U,X(LWRK),X(LWRK1),X(LWRK2),X(LWRK3),
C    *           FON,SRSCC,FREQ,
C    *           X(LDDER),X(LQDER_INIT),S,NCP,FEDER0)
C                FEDER = FEDER0
C         END IF
C
C         CALCULATE GAMDER_SHIFT FIRST
C
          IF (SCC.AND.NI.GE.1) THEN
C           IF (NCP.EQ.1) THEN
              IF (SRSCC) THEN
               CALL DFTB_CALC_DERSHIFT_SHELL(NAT,NI,KI,C,DISTMAT,GAMDER,
     *            GAMDER3,CHAMULS,X(LDSH),DFTB3,NSHELL,INDSH,NSPE,ISPE,
     *            MAXANG,CHAMUL)
                !! LWRK = LDERSH (OLD)
              ELSE
                CALL DFTB_CALC_DERSHIFT(NAT,NI,KI,C,DISTMAT,GAMDER,
     *            GAMDER3,CHAMUL,X(LDSH),DFTB3)
              END IF
C           ELSE
C             CALL DFTB_FCKDER(NCP,NAT,L1,L2,NI1,NI2,KI1,KI2,
C    *             NDIMSH,NGAMMA,NSPIN,IND,MAXANG,ISPE,INDSAO,
C    *             C,S,D,DISTMAT,CHAMUL,CHAMULS,
C    *             GAMMA,GAMMA3,GAMDER,GAMDER3,X(LGAMHES),X(LGAMHES3),
C    *             X(LWRK),X(LWRK1),X(LDSAVE),X(LDSH),X(LSDEROMG),
C    *             X(LQDER_INIT),X(LSDER),
C    *             DFTB3,SRSCC)
C           END IF
C
            IF (IPCM.NE.0.OR.NFMOPCM.GT.0) THEN
              IF (NCP.EQ.2) THEN
                WRITE (IW,*) "NCP=2: NOT YET FOR PCM"
                CALL ABRT
              END IF
              CALL VCLR(X(LPCMTMP),1,NTS)
              CALL VCLR(X(LWRK),1,NAT)
              !! explicit term 2
              !! also prepare implicit term 2-a
              CALL DFTB_HESS_PCM1(NAT,NI,KI,NTS,C,X(LAXYZCT),
     *             X(LAXYZCT+MXTS),X(LAXYZCT+MXTS*2),X(LQSE),X(LWRK),
     *             CHAMUL,X(LPCMTMP),X(LISPHE))
              !! implicit term 1 
              IF (SAVDCQ) THEN
                CALL DAXPY(NTS,ONE,X(LDCDAQ+3*NTS*(NI-1)+NTS*(KI-1)),
     *                     1,X(LPCMTMP),1)
              ELSE
                CALL DFTB_CALC_DCDA(NTS,NI,KI,X(LAXYZCT),
     *               X(LAXYZCT+MXTS),X(LAXYZCT+MXTS*2),
     *               X(LAXYZCT+MXTS*3),X(LISPHE),
     *               X(LEPSPCM),X(LIDDAI),X(LDAI),X(LQSE),X(LPCMTMP))
              END IF
              CALL INDQPROD(1,NTSPAR,LFPNT,NDASC,X(LPCMTMP),
     *                      X(LQTMP),X(LQTMP),1,0)
              CALL DCOPY(NTS,X(LQTMP),1,X(LQSEDER),1)
              IF (SRSCC) THEN
                CALL DFTB_QSHIFT(NAT,NTS,C,X(LAXYZCT),X(LAXYZCT+MXTS),
     *               X(LAXYZCT+MXTS*2),X(LQTMP),X(LWRK))
                CALL DFTB_SHIFT_ATOSH(X(LWRK),X(LWRK1),NAT,
     *               NSHELL,NSPE,ISPE,MAXANG)
                CALL DAXPY(NSHELL,ONE,X(LWRK1),1,X(LDSH),1)
              ELSE
                CALL DFTB_QSHIFT(NAT,NTS,C,X(LAXYZCT),X(LAXYZCT+MXTS),
     *               X(LAXYZCT+MXTS*2),X(LQTMP),X(LWRK))
                CALL DAXPY(NAT,ONE,X(LWRK),1,X(LDSH),1)
              END IF
            END IF
C
            IF (SRSCC) THEN
              CALL SHIFT_TRIANGLE_SHELL(X(LDSH),S,X(LWRK),NAT,NSPE,
     *          NSHELL,L2,ISPE,MAXANG)
            ELSE
              CALL SHIFT_TRIANGLE(NAT,L2,IND,X(LDSH),S,X(LWRK))
            END IF
            IF (NCP.EQ.2) CALL DAXPY(L2,ONE,X(LSDEROMG),1,X(LWRK),1)
C
            IF (LCDFTB) THEN
              CALL DFTB_LCHESS1(NI,KI,X(LFAOLC),X(LSDER),S,X(LDELD),
     *                          GAMLC,GAMDERLC,DISTMAT)
              IF (NSPIN.EQ.2) THEN
                CALL DFTB_LCHESS1(NI,KI,X(LFAOLC+L2),X(LSDER),S,
     *                            X(LDELD+L2),GAMLC,GAMDERLC,DISTMAT)
                CALL DSCAL(2*L2,TWO,X(LFAOLC),1)
              END IF
              CALL DCOPY(L2*NSPIN,X(LFAOLC),1,X(LFAOLC2),1)
            END IF
            DO ISPIN = 1, NSPIN
              IF (LCDFTB) CALL DAXPY(L2,ONE,X(LFAOLC+L2*(ISPIN-1)),1,
     *                               X(LWRK),1)
              CALL DFTB_CPYTSQ(X(LWRK),X(LWRK1),L1,L1)
              CALL DGEMM('T','N',L1-NDOCC(ISPIN),L1,L1,ONE,
     *                   VEC(1,NDOCC(ISPIN)+1,ISPIN),L1,X(LWRK1),L1,
     *                   ZERO,X(LWRK2),L1)
              CALL DGEMM('N','N',L1-NDOCC(ISPIN),NOCC(ISPIN),L1,ONE,
     *                   X(LWRK2),L1,VEC(1,1,ISPIN),L1,ZERO,
     *                   X(LWRK1+NDOCC(ISPIN)),L1)
              CALL CPYSQT(X(LWRK1),X(LSOMGDIJ+L2*(ISPIN-1)),L1,1)
              IF (LCDFTB) CALL DAXPY(L2,-ONE,X(LFAOLC+L2*(ISPIN-1)),1,
     *                               X(LWRK),1)
            END DO
            !! X(LSOMGDIJ) = H_{IJ}^a + non-response d(Omega)/da*S_{IJ}
            CALL VADD(X(LSOMGDIJ),1,X(LHDERIJ),1,X(LSOMGDIJ),1,L2*NSPIN)
            IF (NCP.EQ.2) THEN
              CALL DCOPY(L2*NSPIN,X(LSOMGDIJ),1,X(LHDERIJ),1)
              CALL VCLR(X(LDSH),1,NDIMSH)
            END IF
C
C         CALCULATE NON-RESEPONT TERMS FOR QDER
C
            IF (NCP.EQ.1) THEN
              IF (SRSCC.OR.UHF) THEN
                DO ISPIN = 1, NSPIN
                  CALL DFTB_PREP_QDER(NDIMSH,L1,L2,
     *              INDSAO,X(LQDER_INIT+NDIMSH*(ISPIN-1)),
     *              D(1,ISPIN),X(LDTMP+L2*(ISPIN-1)),S,X(LSDER),X(LWRK),
     *             X(LWRK1),X(LWRK2),X(LWRK3),UHF,X(LDDER+L2*(ISPIN-1)))
                END DO
              ELSE
                CALL DFTB_PREP_QDER(NAT,L1,L2,IND,
     *            X(LQDER_INIT),D,X(LDTMP),S,X(LSDER),X(LWRK),X(LWRK1),
     *            X(LWRK2),X(LWRK3),UHF,X(LDDER))
              END IF
C             write (iw,*) "QDER_INIT"
C             do i = 1, ndimsh
C               tmp = x(lqder_init+i-1)
C               if (uhf) tmp = tmp + x(lqder_init+i-1+ndimsh)
C               write (iw,'(i3,x,f20.10)') i,tmp
C             end do
            END IF
          ELSE
            CALL DCOPY(L2*NSPIN,X(LHDERIJ),1,X(LSOMGDIJ),1)
            CALL VCLR(UINIT,1,L3*NSPIN)
            IF (NCP.EQ.1) THEN
              IF (FREQ) THEN
                CALL VCLR(X(LDDER),1,L3*NSPIN)
              ELSE
                CALL VCLR(X(LDDER),1,L2*NSPIN)
              END IF
              CALL VCLR(X(LQDER_INIT),1,NDIMSH*NSPIN)
            END IF
          END IF
C
C         CALCULATE NON-RESEPONCE TERMS
C
          DO ISPIN = 1, NSPIN
            CALL DFTB_CP_NR(L1,L2,NDOCC(ISPIN),NPART(ISPIN),
     *        EIG(1,ISPIN),OCC(1,ISPIN),UINIT(1,1,ISPIN),
     *        X(LSOMGDIJ+L2*(ISPIN-1)),X(LSDERIJ+L2*(ISPIN-1)),
     *        X(LEIGINV+MAXVIR*MAXOCC*(ISPIN-1)),MAXVIR,MAXOCC,
     *        BETA,EFERMI,DEGTHR(1),NCP,FREQ)
            IF(NFG.NE.0) THEN
              IFG = icurfg
              JFG = jcurfg
              KFG = kcurfg
c             write(6,'("ICHK=",I3)') NDOCC(1)
              CALL FMOESP_CP_DFTB(NI,KI,UINIT,L1,L2,IFG,JFG,KFG,
     *             X(LSDER), S,VEC,EIG,CHAMUL,SAVGOP)
            END IF
C           write (*,*) "hderij"
C           call prtril(x(lhderij),l1)
C           write (*,*) "sderij"
C           call prtril(x(lsderij),l1)
C           write (*,*) "lsomgdij"
C           call prtril(x(lsomgdij),l1)
C      write (iw,*) "Ki, ni = ", ki, ni
C      call prsq(uinit(1,1,ispin),l1,l1,l1)
          END DO
          CALL DCOPY(L1*L1*NSPIN,UINIT,1,U,1)
C
C         PREPARE FOR DE/DR CALCULATION
C
          IF (FON(3).AND..NOT.ROHF) THEN
            !! X(LHDERIJ) IS DUMMY
            IF (NCP.EQ.1) THEN
            CALL DFTB_CALC_DNDR(0,L1,L2,NDOCC(1),NPART(1),ETEMP,
     *        X(LDEDB_INIT),X(LDEDB),X(LDNDR),OCC,EIG,X(LSOMGDIJ),
     *        X(LSDERIJ),X(LHDERIJ),FEDER,SCC)
            END IF
          END IF
          IF (.NOT.SCC)  then
            IF (FON(3).AND..NOT.ROHF) THEN
                CALL DFTB_CALC_DNDR(1,L1,L2,NDOCC(1),NPART(1),ETEMP,
     *            X(LDEDB_INIT),X(LDEDB),X(LDNDR),OCC,EIG,X(LHDERIJ),
     *            X(LSDERIJ),X(LWRK),FEDER,SCC)
            END IF
            GO TO 100
          END IF
C
C         START CP-DFTB ITERATION
C
          MEMMIX = 0
          CALL DCOPY(NDIMSH*NSPIN,X(LQDER_INIT),1,X(LQDER_PREV),1)
          IF (FON(3)) X(LQDER_PREV+NDIMSH) = ZERO
          IF (ITYPMX.LT.0.OR.LCDFTB)
     *      CALL DCOPY(L2*NSPIN,X(LDDER),1,X(LDDER_PREV),1)
          DO ITER = 1, MXITER
            DO ISPIN = 1, NSPIN
              CALL DFTB_CALC_QDER(NCP,NAT,L1,IND,NDOCC(ISPIN),
     *        NPART(ISPIN),X(LQDER+NDIMSH*(ISPIN-1)),
     *        U(1,1,ISPIN),OCC(1,ISPIN),EIG(1,ISPIN),
     *        X(LQDER_INIT+NDIMSH*(ISPIN-1)),
     *        X(LDNDR+L1*(ISPIN-1)),FEDER,BETA,
     *        VEC(1,1,ISPIN),FON(ISPIN),SRSCC.OR.UHF,NSPE,
     *        ISPE,MAXANG,NSHELL,X(LWRK),X(LWRK1),S,INDSAO)
C    *        X(LSOMGDIJ),X(LSDERIJ+L2*(ISPIN-1)),FRQOMG)
              IF (LCDFTB) THEN
                CALL DCOPY(L2,X(LWRK1),1,X(LWRK2+L2*(ISPIN-1)),1)
                !! X(LWRK2) has the derivative of D from independent pairs
                CALL DAXPY(L2,ONE,X(LDDER+L2*(ISPIN-1)),1,
     *                     X(LWRK2+L2*(ISPIN-1)),1)
                !! X(LWRK2) now has the full density matrix derivative
              END IF
            END DO
C
C           CALCULATE CHANGE OF QDER FROM PREVIOUS ITERATION
C
            N = 0
            VMAX = ZERO
            DO ISPIN = 1, NSPIN
              IF (LCDFTB) THEN
                DO I = 1, L2
                  N = N + 1
                  TMP = ABS(X(LWRK2+N-1)-X(LDDER_PREV+N-1))
                  VMAX = MAX(VMAX,TMP)
                END DO
              ELSE
                DO I = 1, NDIMSH
                  N = N + 1
                  TMP = ABS(X(LQDER+N-1)-X(LQDER_PREV+N-1))
                  VMAX = MAX(VMAX,TMP)
                END DO
              END IF
            END DO
C
C           BROYDEN MIXING FOR QDER
C
            IF (ITYPMX.GE.0 .AND. .NOT.LCDFTB) THEN
              IF (FON(3)) THEN
                X(LQDER+NDIMSH) = FEDER
                CALL DFTB_MIX_BROYDEN(ITER,NDIMSH*NSPIN+1,MXITER,MEMMIX,
     *            ALPHA,X(LQDER_PREV),X(LQDER),ERDIIS)
                CALL DCOPY(NDIMSH*NSPIN+1,X(LQDER),1,X(LQDER_PREV),1)
                FEDER = X(LQDER_PREV+NDIMSH)
              ELSE
                CALL DFTB_MIX_BROYDEN(ITER,NDIMSH*NSPIN,MXITER,MEMMIX,
     *            ALPHA,X(LQDER_PREV),X(LQDER),ERDIIS)
                CALL DCOPY(NDIMSH*NSPIN,X(LQDER),1,X(LQDER_PREV),1)
              END IF
            ELSE
              CALL DCOPY(NDIMSH*NSPIN,X(LQDER),1,X(LQDER_PREV),1)
              NNN = L2*NSPIN
              IF (FON(3)) THEN
                X(LWRK2+L2*NSPIN) = FEDER
                NNN = NNN +1
              END IF
C     call cpu_time(t0)
              CALL DFTB_MIX_BROYDEN(ITER,NNN,MXITER,MEMMIX,
     *          ALPHA,X(LDDER_PREV),X(LWRK2),ERDIIS)
C     call cpu_time(t1)
C     t_broyden = t_broyden + t1-t0
              CALL DCOPY(NNN,X(LWRK2),1,X(LDDER_PREV),1)
              IF (FON(3)) THEN
                FEDER = X(LWRK2+L2*NSPIN)
              END IF
C
              CALL DAXPY(L2*NSPIN,-ONE,X(LDDER),1,X(LWRK2),1)
              CALL DFTB_MULLIKEN(L1,L2,X(LWRK2),S,X(LWRK))
              CALL DFTB_MULSA(X(LWRK),X(LWRK1),X(LWRK3),L1,NSHELL2,NAT)
              CALL DCOPY(NDIMSH,X(LQDER_INIT),1,X(LQDER),1)
              IF (SRSCC.OR.UHF) THEN
                CALL DAXPY(NDIMSH,ONE,X(LWRK1),1,X(LQDER),1)
              ELSE
                CALL DAXPY(NDIMSH,ONE,X(LWRK3),1,X(LQDER),1)
              END IF
              IF (NSPIN.EQ.2) THEN
                CALL DFTB_MULLIKEN(L1,L2,X(LWRK2+L2),S,X(LWRK))
               CALL DFTB_MULSA(X(LWRK),X(LWRK1),X(LWRK3),L1,NSHELL2,NAT)
                CALL DCOPY(NDIMSH,X(LQDER_INIT+NDIMSH),1,
     *                     X(LQDER+NDIMSH),1)
                IF (SRSCC.OR.UHF) THEN
                  CALL DAXPY(NDIMSH,ONE,X(LWRK1),1,X(LQDER+NDIMSH),1)
                ELSE
                  CALL DAXPY(NDIMSH,ONE,X(LWRK3),1,X(LQDER+NDIMSH),1)
                END IF
              END IF
            END IF
C
            CALL DCOPY(NDIMSH,X(LQDER),1,X(LWRK),1)
            IF (UHF) THEN
              CALL VADD(X(LWRK),1,X(LQDER+NDIMSH),1,X(LWRK),1,NDIMSH)
              IF (.NOT.SRSCC) THEN
                CALL DFTB_SHIFT_SHTOA(X(LWRK),X(LWRK1),NAT,NSHELL,NSPE,
     *            ISPE,MAXANG)
                CALL DCOPY(NAT,X(LWRK1),1,X(LWRK),1)
              END IF
            END IF
C
C           CONSTRUCT A NEW SHIFT MATRIX WITH THE UPDATED QDER
C
            CALL VCLR(X(LDERSHIFT),1,NGAMMA)
            DO NN = 1, NGAMMA
              DO NC = 1, NGAMMA
                CALL DFTB_CNVSQ(NN,NC,NSEQ)
                X(LDERSHIFT+NN-1)
     *            = X(LDERSHIFT+NN-1) + GAMMA(NSEQ)*X(LWRK+NC-1)
              END DO
            END DO
C
            IF (DFTB3) THEN
              !! GENERATE ATOM-RESOLVED MULLIKEN DER.
              IF (SRSCC) THEN
                CALL DCOPY(NDIMSH,X(LQDER),1,X(LWRK1),1)
                IF (UHF) THEN
                  CALL DAXPY(NDIMSH,ONE,X(LQDER+NDIMSH),1,X(LWRK1),1)
                END IF
                CALL DFTB_SHIFT_SHTOA(X(LWRK1),X(LWRK+NSHELL+1),NAT,
     *            NSHELL,NSPE,ISPE,MAXANG)
C
                ISH0 = 0
                DO I = 1, NAT
                  JSH0 = 0
                  DERMULI = X(LWRK+NSHELL+ I)
                  DO J = 1, NAT
                    DERMULJ = X(LWRK+NSHELL+ J)
                    DO II = 1, MAXANG(ISPE(I))
                      ISH = ISH0 + II
                      DO JJ = 1, MAXANG(ISPE(J))
                        JSH = JSH0 + JJ
                        DERSHLJ = X(LWRK+JSH-1)
                        X(LDERSHIFT+ISH-1) = X(LDERSHIFT+ISH-1)
     *                    + ONE_THIRD*
     *                     (DERMULI*CHAMULS(JSH)*GAMMA3(ISH,JSH)
     *                    + DERSHLJ*CHAMUL(I)   *GAMMA3(ISH,JSH)
     *                    + DERMULJ*CHAMULS(JSH)*GAMMA3(JSH,ISH)
     *                    + DERSHLJ*CHAMUL(J)   *GAMMA3(JSH,ISH))
                        DO KK = 1, MAXANG(ISPE(I))
                          KSH = ISH0 + KK
                          DERSHLK = X(LWRK+KSH-1)
                          X(LDERSHIFT+ISH-1) = X(LDERSHIFT+ISH-1)
     *                      + ONE_THIRD*
     *                       (DERSHLK*CHAMULS(JSH)*GAMMA3(KSH,JSH)
     *                      + DERSHLJ*CHAMULS(KSH)*GAMMA3(KSH,JSH))
                        END DO
                      END DO
                    END DO
                    JSH0 = JSH0 + MAXANG(ISPE(J))
                  END DO
                  ISH0 = ISH0 + MAXANG(ISPE(I))
                END DO
              ELSE
                DO NN = 1, NGAMMA
                  DO NC = 1, NGAMMA
                    X(LDERSHIFT+NN-1) = X(LDERSHIFT+NN-1) + TWO*(
     *               + ONE_THIRD*GAMMA3(NN,NC)*CHAMUL(NN)*X(LWRK+NC-1)
     *               + ONE_THIRD*GAMMA3(NC,NN)*CHAMUL(NC)*X(LWRK+NC-1)
     *               + ONE_THIRD*GAMMA3(NN,NC)*CHAMUL(NC)*X(LWRK+NN-1))
                  END DO
                END DO
              END IF
              IF (NFG.NE.0) THEN
                DO NN = 1, NGAMMA !! NAT
                  X(LDERSHIFT+NN-1) = X(LDERSHIFT+NN-1)
     *             + TWO*X(LSHIFT3+NN-1)*X(LWRK+NN-1)
                END DO
              END IF
            END IF
C
            IF (LCDFTB) THEN
              CALL VCLR(X(LFAOLC),1,L2*NSPIN)
              CALL DFTB_LCSHIFT2(1,NAT,NSPIN,L1,ISPE,MAXANG,IND,
     *                          X(LDDER_PREV),DUM,S,GAMLC,X(LFAOLC),ELC)
            END IF
C           ----- CALCULATE CHARGES INDUCED BY MULLIKEN DER. -----
C           START INDQPROD WITH ITER+1, BECAUSE WE ALREADY HAVE
C           AN INITIAL GUESS WITHOUT ASCS INDUCED BY MULLIKEN
C           DERIVATIVES
C
            IF (IPCM.NE.0.OR.NFMOPCM.GT.0) THEN
              !! implicit term 2-b
              !! use dq/da (X(LWRK))
              IF (SRSCC.OR.UHF) THEN
                !! Convert shell-resolved X(LWRK) to atom-resolved
                CALL DCOPY(NDIMSH,X(LQDER),1,X(LWRK1),1)
                IF (UHF) THEN
                  CALL DAXPY(NDIMSH,ONE,X(LQDER+NDIMSH),1,X(LWRK1),1)
                END IF
                CALL DFTB_SHIFT_SHTOA(X(LWRK1),X(LWRK),NAT,NSHELL,NSPE,
     *            ISPE,MAXANG)
              END IF
              CALL DFTB_QPOTEN(NAT,NTS,C,X(LAXYZCT),X(LAXYZCT+MXTS),
     *                         X(LAXYZCT+MXTS*2),X(LWRK),X(LEPOT))
              CALL INDQPROD(ITER+1,NTSPAR,LFPNT,NDASC,X(LEPOT),
     *                      X(LQSEDER),X(LQSEDER),1,0)
              IF (SRSCC) THEN
                CALL VCLR(X(LWRK1),1,NAT)
                CALL DFTB_QSHIFT(NAT,NTS,C,X(LAXYZCT),X(LAXYZCT+MXTS),
     *               X(LAXYZCT+MXTS*2),X(LQSEDER),X(LWRK1))
                CALL DFTB_SHIFT_ATOSH(X(LWRK1),X(LWRK2),NAT,NSHELL,
     *               NSPE,ISPE,MAXANG)
                CALL DAXPY(NGAMMA,ONE,X(LWRK2),1,X(LDERSHIFT),1)
              ELSE
                CALL DFTB_QSHIFT(NAT,NTS,C,X(LAXYZCT),X(LAXYZCT+MXTS),
     *               X(LAXYZCT+MXTS*2),X(LQSEDER),X(LDERSHIFT))
              END IF
            END IF
C
C           CONVERT ATOM-RESOLVED SHIFT TO SHELL-RESOLVED SHIFT,
C           AND CALCULATE SPIN CONTRIBUTION
C
            IF (UHF) THEN
              IF (.NOT.SRSCC) THEN
                CALL DCOPY(NAT,X(LDERSHIFT),1,X(LWRK),1)
                CALL DFTB_SHIFT_ATOSH(X(LWRK),X(LDERSHIFT),NAT,NSHELL,
     *            NSPE,ISPE,MAXANG)
              END IF
C             FORM d(q^a-q^b)/dR (DERIVATIVE OF SPIN DENSITY) AS X(LWRK)
              CALL DCOPY(NSHELL,X(LQDER),1,X(LWRK),1)
              CALL VSUB(X(LQDER+NDIMSH),1,X(LWRK),1,
     *          X(LWRK),1,NSHELL)
C
              CALL VCLR(X(LSHIFTSP),1,NDIMSH)
              CALL DFTB_SPIN_SHIFT(X(LSHIFTSP),X(LWRK),SPNCST,NAT,
     *          NSHELL,NSPE,ISPE,MAXANG)
              CALL DCOPY(NDIMSH,X(LDERSHIFT),1,X(LWRK1),1)
            END IF
C
C           CALCULATE RESPONSE PART OF S*d(OMEGA)/dR AS SDSHIJ
C
            IF (KI.EQ.1) CALL VCLR(X(LSDEROMG),1,L2)
            DO ISPIN = 1, NSPIN
              IF (SRSCC.OR.UHF) THEN
                IF (UHF) THEN
                  IF (ISPIN.EQ.1) THEN
                   CALL VADD(X(LSHIFTSP),1,X(LDERSHIFT),1,X(LDERSHIFT),
     *               1,NDIMSH)
                  ELSE
                    CALL DCOPY(NDIMSH,X(LWRK1),1,X(LDERSHIFT),1)
                    CALL VSUB(X(LSHIFTSP),1,X(LDERSHIFT),1,X(LDERSHIFT),
     *                1,NDIMSH)
                  END IF
                END IF
                CALL SHIFT_TRIANGLE_SHELL(X(LDERSHIFT),S,X(LSDEROMG),
     *            NAT,NSPE,NSHELL,L2,ISPE,MAXANG)
              ELSE
               CALL SHIFT_TRIANGLE(NAT,L2,IND,X(LDERSHIFT),S,
     *           X(LSDEROMG))
              END IF
              IF (LCDFTB) THEN
               CALL DAXPY(L2,ONE,X(LFAOLC+L2*(ISPIN-1)),1,X(LSDEROMG),1)
              END IF
              CALL DFTB_CPYTSQ(X(LSDEROMG),X(LWRK),L1,L1)
              CALL DGEMM('T','N',L1-NDOCC(ISPIN),L1,L1,ONE,
     *                   VEC(1,NDOCC(ISPIN)+1,ISPIN),L1,X(LWRK),L1,
     *                   ZERO,X(LWRK2),L1)
              CALL DGEMM('N','N',L1-NDOCC(ISPIN),NOCC(ISPIN),L1,ONE,
     *                   X(LWRK2),L1,VEC(1,1,ISPIN),L1,ZERO,
     *                   X(LWRK+NDOCC(ISPIN)),L1)
C             NSEQ=1
C             DO IO = 1, L1
C               CALL DCOPY(IO,X(LSDEROMG+NSEQ-1),1,X(LWRK+L1*(IO-1)),1)
C               NSEQ = NSEQ + IO
C             END DO
C             CALL DSYMM('L','U',L1,NOCC(ISPIN),ONE,
C    *                   X(LWRK),L1,VEC(1,1,ISPIN),L1,
C    *                   ZERO,X(LWRK2),L1)
C             CALL DGEMM('T','N',L1-NDOCC(ISPIN),NOCC(ISPIN),L1,ONE,
C    *                   VEC(1,NDOCC(ISPIN)+1,ISPIN),L1,X(LWRK2),L1,
C    *                   ZERO,X(LWRK+NDOCC(ISPIN)),L1)
C
C             UPDATE de/dR, de_F/dR, and dn/dR VECTORS
C
              IF (FON(3).AND..NOT.ROHF) THEN
C               IF (NCP.EQ.1) THEN
                  CALL DFTB_CALC_DNDR(1,L1,L2,NDOCC(1),NPART(1),ETEMP,
     *              X(LDEDB_INIT),X(LDEDB),X(LDNDR),OCC,EIG,
     *              X(LHDERIJ),X(LSDERIJ),X(LWRK),FEDER,SCC)
                  CALL CPYSQT(X(LWRK),X(LSOMGDIJ),L1,1)
                  CALL DAXPY(L2,ONE,X(LHDERIJ),1,X(LSOMGDIJ),1)
C               ELSE
C                 !! X(LWRK) contains only d2q/dadb contribution
C                 CALL CPYSQT(X(LWRK),X(LSOMGDIJ),L1,1)
C                 !! Now, X(LOMGDIJ) corresponds to F_{ij}^{ab}
C                 CALL DFTB_FEDER(L1,L2,NSPIN,NDOCC,NOCC,FEDER0,FEDER,
C    *                            OCC,BETA,X(LSOMGDIJ))
C               END IF
              END IF
C
C             CONSTRUCT NEW U MATRIX
C
              CALL DFTB_HESS_UUPDATE(U(1,1,ISPIN),UINIT(1,1,ISPIN),
     *          X(LWRK),L1,NDOCC(ISPIN),NPART(ISPIN),MAXOCC,MAXVIR,
     *          X(LEIGINV+MAXVIR*MAXOCC*(ISPIN-1)))
            END DO !! END OF SPIN LOOP
            IF (ITER.GE.2.AND.VMAX.LE.THRES) THEN
              IF (NI.GE.1.AND.NI.LE.NAT) NITER = NITER + ITER
              IF (NI.GT.NAT) NITER2 = NITER2 + ITER
C             WRITE (*,'(" ROOT",I3," - ",I1," CONVERGED IN", I3,
C    *          " CYCLES: DELTA Q_DER = ",E10.3)') NI,KI,ITER,VMAX
              EXIT
            END IF
          END DO !! END OF CP ITERATION
C         WRITE (*,*) "U MATRIX"
C         CALL PRSQ(U,L1,L1,L1)
C         WRITE (*,*)
C         WRITE (*,*) "SDER MATRIX"
C         CALL PRTRIL(X(LSDERIJ),L1)
C         WRITE (*,*)
C         WRITE (*,*) "DEGENERATED ORBITALS"
C         DO IO = 1, L1
C           DO JO = IO+1, L1
C             IF (ABS(EIG(IO)-EIG(JO)).LE.1.0D-06) THEN
C               CALL DFTB_CNVSQ(IO,JO,NSEQ)
C               WRITE (*,'(I3,X,I3,4(X,F20.10))')
C    *            IO,JO,U(IO,JO),U(JO,IO),U(IO,JO)+U(JO,IO),
C    *            X(LSDERIJ+NSEQ-1)
C             END IF
C           END DO
C         END DO
          CALL RETFM(MEMMIX)
C
C         FOR UHF CALCULATION:
C         QDER(1,1) :: TOTAL CHARGE DERIVATIVE
C         QDER(1,2) :: SPIN DERIVATIVE
C
          IF (UHF) THEN
            CALL DCOPY(NDIMSH,X(LQDER),1,X(LWRK),1)
            CALL VADD(X(LQDER+NDIMSH),1,X(LQDER),1,X(LQDER),1,NDIMSH)
            CALL DSCAL(NDIMSH,-ONE,X(LQDER+NDIMSH),1)
            CALL VADD(X(LWRK),1,X(LQDER+NDIMSH),1,
     *        X(LQDER+NDIMSH),1,NDIMSH)
          END IF
C
C         CALCULATE d(OMEGA)/dR
C         d(OMEGA)/dR IS SYMMETRIC
C
          IF (UHF.AND..NOT.SRSCC) THEN
            CALL DCOPY(NDIMSH,X(LQDER),1,X(LWRK),1)
            CALL DFTB_SHIFT_SHTOA(X(LWRK),
     *        X(LQDER),NAT,NSHELL,NSPE,ISPE,MAXANG)
          END IF
          IF (NI.EQ.0) CALL VCLR(X(LDSH),1,NGAMMA)
          DO NM = 1, NGAMMA
            TMP1 = X(LDSH+NM-1)
            TMP2 = ZERO
            DO NC = 1, NGAMMA
              CALL DFTB_CNVSQ(NM,NC,NSEQ)
              TMP2 = TMP2 + GAMMA(NSEQ)*X(LQDER+NC-1)
            END DO
            X(LDOMEGA+NM-1) = TMP1 + TMP2
          END DO
C
          IF (DFTB3) THEN
            IF (SRSCC) THEN
              CALL DCOPY(NDIMSH,X(LQDER),1,X(LWRK1),1)
              IF (ARAMAN) CALL DCOPY(NDIMSH,X(LQDER),1,X(LQDER_INIT),1)
              CALL DFTB_SHIFT_SHTOA(X(LWRK1),X(LWRK+NSHELL+1),NAT,
     *          NSHELL,NSPE,ISPE,MAXANG)
C
              ISH0 = 0
              DO I = 1, NAT
                JSH0 = 0
                DERMULI = X(LWRK+NSHELL+ I)
                DO J = 1, NAT
                  DERMULJ = X(LWRK+NSHELL+ J)
                  DO II = 1, MAXANG(ISPE(I))
                    ISH = ISH0 + II
                    DO JJ = 1, MAXANG(ISPE(J))
                      JSH = JSH0 + JJ
                      DERSHLJ = X(LQDER+JSH-1)
                      X(LDOMEGA+ISH-1) = X(LDOMEGA+ISH-1)
     *                  + ONE_THIRD*
     *                   (DERMULI*CHAMULS(JSH)*GAMMA3(ISH,JSH)
     *                  + DERSHLJ*CHAMUL(I)   *GAMMA3(ISH,JSH)
     *                  + DERMULJ*CHAMULS(JSH)*GAMMA3(JSH,ISH)
     *                  + DERSHLJ*CHAMUL(J)   *GAMMA3(JSH,ISH))
                      DO KK = 1, MAXANG(ISPE(I))
                        KSH = ISH0 + KK
                        DERSHLK = X(LQDER+KSH-1)
                        X(LDOMEGA+ISH-1) = X(LDOMEGA+ISH-1)
     *                    + ONE_THIRD*
     *                     (DERSHLK*CHAMULS(JSH)*GAMMA3(KSH,JSH)
     *                    + DERSHLJ*CHAMULS(KSH)*GAMMA3(KSH,JSH))
                      END DO
                    END DO
                  END DO
                  JSH0 = JSH0 + MAXANG(ISPE(J))
                END DO
                ISH0 = ISH0 + MAXANG(ISPE(I))
              END DO
            ELSE
              DO NN = 1, NAT
                DO NC = 1, NAT
                  X(LDOMEGA+NN-1) = X(LDOMEGA+NN-1) + TWO*(
     *             + ONE_THIRD*GAMMA3(NN,NC)*CHAMUL(NN)*X(LQDER+NC-1)
     *             + ONE_THIRD*GAMMA3(NC,NN)*CHAMUL(NC)*X(LQDER+NC-1)
     *             + ONE_THIRD*GAMMA3(NN,NC)*CHAMUL(NC)*X(LQDER+NN-1))
                END DO
              END DO
            END IF
          END IF
C
C         ----- PCM CONTRIBUTIONS -----
C         QSEDER HAS ASC INDUCED BY MULLIKEN DERIVATIVE.
C         SHIFT CONTRIBUTIONS OTHER THAN QSEDER HAVE BEEN ADDED AS
C         X(LDSH) ABOVE, SO SHIFT HERE IS DONE BEFORE FORMING TRUE
C         QSEDER (\bar{q}^a)
C
          IF (IPCM.NE.0.OR.NFMOPCM.GT.0) THEN
            IF (SRSCC) THEN
              CALL VCLR(X(LWRK),1,NAT)
              CALL DFTB_QSHIFT(NAT,NTS,C,X(LAXYZCT),X(LAXYZCT+MXTS),
     *             X(LAXYZCT+MXTS*2),X(LQSEDER),X(LWRK))
              CALL DFTB_SHIFT_ATOSH(X(LWRK),X(LWRK1),NAT,NSHELL,NSPE,
     *             ISPE,MAXANG)
              CALL DAXPY(NSHELL,ONE,X(LWRK1),1,X(LDOMEGA),1)
            ELSE
              CALL DFTB_QSHIFT(NAT,NTS,C,X(LAXYZCT),X(LAXYZCT+MXTS),
     *             X(LAXYZCT+MXTS*2),X(LQSEDER),X(LDOMEGA))
            END IF
            CALL DAXPY(NTS,ONE,X(LQTMP),1,X(LQSEDER),1)
            !! THIS X(LQSEDER) IS TRUE \bar{q}^a
          END IF
C
          IF (UHF) THEN
            IF (.NOT.SRSCC) THEN
              CALL DCOPY(NAT,X(LDOMEGA),1,X(LWRK),1)
              CALL DFTB_SHIFT_ATOSH(X(LWRK),X(LDOMEGA),NAT,NSHELL,NSPE,
     *          ISPE,MAXANG)
            END IF
            CALL VCLR(X(LSHIFTSP),1,NDIMSH)
            CALL DFTB_SPIN_SHIFT(X(LSHIFTSP),X(LQDER+NDIMSH),SPNCST,NAT,
     *        NSHELL,NSPE,ISPE,MAXANG)
            CALL DCOPY(NSHELL,X(LDOMEGA),1,X(LDOMEGA+NSHELL),1)
            CALL VADD(X(LSHIFTSP),1,X(LDOMEGA),1,X(LDOMEGA),1,NSHELL)
            CALL VSUB(X(LSHIFTSP),1,X(LDOMEGA+NSHELL),1,
     *        X(LDOMEGA+NSHELL),1,NSHELL)
          END IF
C
C         WE HAVE OBTAINED U MATRIX AND QDER NOW.
C
C         ADD QDER HESSIAN TERMS
C
          IF (NI.GE.1.AND.NCP.EQ.1) THEN
            IF (SRSCC) THEN
              !! CONVERT FROM SHELL-RESOLVED TO ATOM-RESOLVED QDER
              CALL DFTB_SHIFT_SHTOA(X(LQDER),X(LWRK),NAT,NSHELL,NSPE,
     *          ISPE,MAXANG)
              CALL DFTB_QDER(NI,KI,NAT,C,X(LQDER),GAMDER,GAMDER3,
     *          CHAMULS,DISTMAT,HESSIAN,DFTB3,SRSCC,NDIMSH,NSPE,ISPE,
     *          MAXANG,IPCM.NE.0.OR.NFMOPCM.GT.0,X(LWRK),CHAMUL,NGAMMA)
              CALL DCOPY(NAT,X(LWRK),1,X(LQDER),1)
            ELSE
              CALL DFTB_QDER(NI,KI,NAT,C,X(LQDER),GAMDER,GAMDER3,
     *          CHAMUL,DISTMAT,HESSIAN,DFTB3,SRSCC,NDIMSH,NSPE,ISPE,
     *          MAXANG,IPCM.NE.0.OR.NFMOPCM.GT.0,X(LQDER),CHAMUL,NGAMMA)
            END IF
          END IF
C
C         d(OMEGA)/dR + d(EIG(IO))/dR TERMS
C
          !! X(LDOMEGA) d(OMEGA)/dR ON ATOM BASIS
          DO ISPIN = 1, NSPIN
            IF (SRSCC.OR.UHF) THEN
              CALL SHIFT_TRIANGLE_SHELL(X(LDOMEGA+NDIMSH*(ISPIN-1)),
     *          S,X(LWRK2),NAT,NSPE,NSHELL,L2,ISPE,MAXANG)
            ELSE
              CALL SHIFT_TRIANGLE(NAT,L2,IND,X(LDOMEGA),S,X(LWRK2))
            END IF
            IF (LCDFTB) THEN
             CALL DAXPY(L2,ONE,X(LFAOLC +L2*(ISPIN-1)),1,X(LWRK2),1)
             CALL DAXPY(L2,ONE,X(LFAOLC2+L2*(ISPIN-1)),1,X(LWRK2),1)
            END IF
            CALL DFTB_CPYTSQ(X(LWRK2),X(LWRK1),L1,L1)
            !! OCC-OCC SPACE IS ALSO NEEDED
            !! possibly full space?
            NOCCTMP = NOCC(ISPIN)
            IF (ARAMAN.OR.CPF.OR.CPA2.OR.FULL3RD) NOCCTMP = L1
            CALL DGEMM('T','N',L1,L1,L1,ONE,
     *                 VEC(1,1,ISPIN),L1,X(LWRK1),L1,
     *                 ZERO,X(LWRK2),L1)
            CALL DGEMM('N','N',L1,NOCCTMP,L1,ONE,
     *                 X(LWRK2),L1,VEC(1,1,ISPIN),L1,ZERO,
     *                 X(LWRK),L1)
            CALL CPYSQT(X(LWRK),X(LSOMGDIJ+L2*(ISPIN-1)),L1,1)
          END DO
          CALL VADD(X(LSOMGDIJ),1,X(LHDERIJ),1,X(LHDERIJ),1,L2*NSPIN)
C
C         ----- ADD PCM CONTRIBUTIONS TO HESSIAN -----
C         LQSE   : ASC INDUCED BY D
C         LQSEDER: dq/da
C
          IF (IPCM.NE.0.OR.NFMOPCM.GT.0.AND.NI.GE.1) THEN
            CALL DFTB_HESS_PCM2(NAT,NTS,C,X(LAXYZCT),X(LAXYZCT+MXTS),
     *                          X(LAXYZCT+MXTS*2),X(LQSE),X(LQSEDER),
     *                          X(LQDER),CHAMUL,X(LAXYZCT+MXTS*3),
     *                          X(LISPHE),X(LEPSPCM),X(LIDDAI),X(LDAI),
     *                          X(LPCMTMP),HESSIAN(1,3*(NI-1)+KI),
     *                          X(LDCDAQ),SAVDCQ)
          END IF
C
C         ADD U-DEPENDENT TERM AND DE/DB DEPENDENT TERMS
C
  100     CONTINUE
          CALL DFTB_CALC_FDER_DEDB(L1,L2,NSPIN,X(LFDERIJ),X(LDEDB),
     *      X(LHDERIJ),X(LSDERIJ),EIG)
          CALL DFTB_CALC_DWDER(NCP,L1,L2,NSPIN,NDOCC,NPART,EIG,OCC,VEC,
     *                         X(LDEDB),X(LDNDR),X(LHDERIJ),U,
     *                         X(LDDER),X(LWDER),X(LWRK),X(LWRK1),
     *                         X(LWRK2),X(LSDERIJ),FON,SCC,BETA,FEDER)
C       if (ncp.eq.2) call prtril(x(ldder),l1)
          IF (NI.EQ.0) THEN
            CALL DCOPY(L2*NSPIN,X(LHDERIJ),1,
     *                          X(LFDERF+L2*NSPIN*(IFREQ2-1)),1)
            IF (FREQ) THEN
              CALL DCOPY(L3*NSPIN,X(LDDER),1,
     *                            X(LDDERF+L3*NSPIN*(IFREQ2-1)),1)
            ELSE
              CALL DCOPY(L2*NSPIN,X(LDDER),1,
     *                            X(LDDERF+L2*NSPIN*(IFREQ2-1)),1)
            END IF
            IF (FON(3)) THEN
              IF (KI.LE.3 .AND. IFREQ.EQ.1) FEDERTMP(KI) = FEDER
              FEDERTMP(4)  = EFERMI
              X(LFEDER+IFREQ2-1) = FEDER
            END IF
            CALL DCOPY(NDIMSH*NSPIN,X(LQDER),1,
     *                 X(LQDERF+NDIMSH*NSPIN*(KI-1)),1)
C
            IF (NCP.EQ.1 .AND. IFREQ.EQ.1) THEN
              DO ISPIN = 1, NSPIN
                LOC = L3*NSPIN*(KI-1) + L3*(ISPIN-1)
                CALL DFTB_CPYTSQ(X(LHDERIJ+L2*(ISPIN-1)),
     *            X(LFDERFSC+LOC),L1,L1)
                DO I = 1, L3
                  X(LFDERFSC+I-1+LOC)
     *              = X(LFDERFSC+I-1+LOC)*X(LEIGINV2+I-1+L3*(ISPIN-1))
                END DO
                CALL DFTB_DEGCLR(L1,NDEG2,ISPIN,X(LLSTDEG2),
     *                           X(LFDERFSC+LOC))
              END DO
              CALL DCOPY(L1*NSPIN,X(LDEDB),1,
     *                            X(LEDERF+L1*NSPIN*(KI-1)),1)
              IF (FON(3)) THEN
                CALL DCOPY(L1*NSPIN,X(LDNDR),1,
     *                              X(LFONDER+L1*NSPIN*(KI-1)),1)
              END IF
            END IF
C
            IF (MASWRK.AND.SOME) THEN
              WRITE (IW,'(" CONVERGED AFTER ",I2," ITERATIONS")') ITER
              IF (.NOT.SRSCC) THEN
                DO KJ = 1, 3
                  V(KJ) = -DDOT(NAT,X(LQDER),1,C(KJ,1),3)
                END DO
              ELSE
                DO KJ = 1, 3
                  TMP = ZERO
                  ISH0 = 0
                  DO NN = 1, NAT
                    DO II = 1, MAXANG(ISPE(NI))
                      ISH = ISH0 + II
                      TMP = TMP + X(LQDER+ISH-1)*C(KJ,NN)
                    END DO
                    ISH0 = ISH0 + MAXANG(ISPE(NI))
                  END DO
                  V(KJ) = -TMP
                END DO
              END IF
C             DO KJ = 1, 3
C               WRITE (IW,'(" HYPERPOLARIZABILITY ALONG ",A1,"=",
C    *            F20.10)') ".",V(KJ)
C             END DO
              IF (KI.EQ.1) CHARA(1:1) = 'X'
              IF (KI.EQ.2) CHARA(1:1) = 'Y'
              IF (KI.EQ.3) CHARA(1:1) = 'Z'
              DO KJ = 1, 3
                IF (KJ.EQ.1) CHARA(NCP+1:NCP+1) = 'X'
                IF (KJ.EQ.2) CHARA(NCP+1:NCP+1) = 'Y'
                IF (KJ.EQ.3) CHARA(NCP+1:NCP+1) = 'Z'
C               IF (NCP.EQ.1) THEN
                  WRITE (IW,'(" POLARIZABILITY      ",A2,"  =",
     *              F20.10)') CHARA(1:2),V(KJ)
C               ELSE IF (NCP.EQ.2) THEN
C                 WRITE (IW,'(" HYPERPOLARIZABILITY ",A3, " =",
C    *              F20.10)') CHARA(1:3),V(KJ)
C               END IF
              END DO
C             IF (KI.EQ.3) THEN
C               CALL SYSTEM_CLOCK(ICAR1,ICR,ICM)
C               CALL CPU_TIME(TAR1)
C               TAR  = TAR1-TAR0
C               TARW = DBLE(ICAR1-ICAR0)/DBLE(ICR)
C               WRITE (IW,'(X,"CP-DFTB FOR E.FIELD TOOK ",F10.2," (",
C    *            F10.2,") CPU (WALL) SECONDS")') TAR,TARW
C               TAR  = ZERO
C               TARW = ZERO
C             END IF
            END IF
            CYCLE
          END IF
C
C         ----- COMPUTE D^a*F^b - W^a*S^b -----
C
          IF (NCP.EQ.1) THEN
            CALL DFTB_DEDB(HESSIAN(1,3*(NI-1)+KI),HOPAO,
     *                     X(LDOMEGA),SHIFT,SHIFTSA,
     *                     SHIFTSB,C,D,
     *                     NAT,L1,L2,NSPIN,NDIMSH,
     *                     NSPE,ISPE,MAXANG,IND,INDSH,INDSAO,
     *                     SCC,SRSCC,UHF,MASWRK,
     *                     X(LDDER),X(LWDER),MAXBND,INDBND,DERMAT)
            !! Compute S^a*S^b, S^a*D^b, S^a*g^b, g^a*D^b
            IF (LCDFTB) THEN
              CALL DFTB_LCHESS2(NI,KI,NSPIN,X(LFAOLC),X(LSDER),S,
     *                          X(LDELD),X(LDDER),GAMLC,GAMDERLC,
     *                          DISTMAT,MAXBND,INDBND,DERMAT,
     *                          MASWRK,HESSIAN,NPROC)
              IF (NSPIN.EQ.2)
     *        CALL DFTB_LCHESS2(NI,KI,NSPIN,X(LFAOLC+L2),X(LSDER),S,
     *                          X(LDELD+L2),X(LDDER+L2),GAMLC,GAMDERLC,
     *                          DISTMAT,MAXBND,INDBND,DERMAT,
     *                          MASWRK,HESSIAN,NPROC)
            END IF
          END IF
          IF (IAND(MODHSS,8+16+128).NE.0) THEN
            LOC = L2*NSPIN*(3*(NI-1)+KI-1)
            IF (NCP.EQ.2) THEN
              LOC = L2*NSPIN*3*NAT
     *            + L2*NSPIN*(6*(NI2-1)
     *            + 9*(NAT*(NI2-1)-NI2*(NI2-1)/2) + KI-1)
              IF (NI1.NE.NI2) THEN
                LOC = LOC + L2*NSPIN*6 + L2*NSPIN*(NI1-NI2-1)*9
              END IF
            END IF
C                write (*,*) "loc = ", loc/l2
            CALL DCOPY(L2*NSPIN,X(LDDER),1,X(LDSAVE+LOC),1)
            CALL DCOPY(L2*NSPIN,X(LHDERIJ),1,X(LFSAVE+LOC),1)
            CALL DCOPY(L2*NSPIN,X(LSDERIJ),1,X(LSSAVE+LOC),1)
C       if (ncp.eq.1) then
C       write (*,*) "dsave"
C       call prtril(x(ldsave+loc),l1)
C       write (*,*) "fsave"
C       call prtril(x(lfsave+loc),l1)
C       write (*,*) "ssave"
C       call prtril(x(lssave+loc),l1)
C       end if
C           IF (FON(3)) X(LFESAVE+NSPIN*(3*(NI-1)+KI-1)) = FEDER
            IF (FON(3)) X(LFESAVE+LOC/L2) = FEDER
          END IF
C
C         GET U^a (Vir-Occ)
          IF(NFG.NE.0.and.nbody.gt.1) THEN
            NVOC = NDOCC(1)*(L1-NDOCC(1))
c           WRITE (6,'("CHECK=",I3)') NVOC
            IXYZ = (NI-1)*3 + KI
            LPNT = LYAWRK   + NVOC*(ixyz-1)
            DO ivir=1,L1-NDOCC(1)
              DO jocc=1,NDOCC(1)
                X(LPNT) = U(ivir+NDOCC(1),jocc,1)
                LPNT    = LPNT  + 1
              END DO
            END DO
            LPNT = LYAWRK   + NVOC*NAT*3 + NAT*(ixyz-1)
            CALL DCOPY(NAT,X(LQDER),1,X(LPNT),1)
          END IF
C         WRITE (IW,*) "QDER (ALPHA)"
C         DO IO = 1, NAT
C           WRITE (IW,'(I3,1(X,F20.10))') IO,X(LQDER+IO-1)
C         END DO
C         WRITE (*,*)
C         call abrt
C         WRITE (*,*) "QDER,DEDB,DNDR"
C         DO IO = 1, NAT
C           WRITE (*,'(I3,3(X,F20.10))')
C    *      IO,X(LDEDB+IO-1),X(LDNDR+IO-1),X(LQDER+IO-1)
C         END DO
C         DO IO = NAT+1, L1
C           WRITE (*,'(I3,3(X,F20.10))')
C    *      IO,X(LDEDB+IO-1),X(LDNDR+IO-1),ZERO
C         END DO
C         WRITE (*,'(" CONVERGED DE_F/DR = ",F20.10)') FEDER
CD        call cpu_time(t5)
C
C         CALCULATE QDER IF NCC TO CALCULATE DIPOLE DERIVATIVE
C
          IF (.NOT.SCC) THEN
            IF (FON(3)) THEN
              IF (NDOCC(1).GE.1)
     *          CALL DFTB_DPART(L1,L2,NDOCC,OCC,X(LDTMP),VEC)
            ELSE
              CALL DCOPY(L2,D,1,X(LDTMP),1)
            END IF
            CALL VCLR(X(LSDER),1,L2)
            !! CALCULATE OVERLAP DERIVATIVE IN AO
            CALL VCLR(BUP,1,81)
            DO NJ = 1, NAT
              IF (NI.EQ.NJ) CYCLE
              V(1) = C(1,NJ)-C(1,NI)
              V(2) = C(2,NJ)-C(2,NI)
              V(3) = C(3,NJ)-C(3,NI)
              CALL DFTB_SKODE_HELPA(1,NI,NJ,KI,0,0,0,V,AUP,BUP,.FALSE.,
     *          .TRUE.,SKP)
              IF (SKP) CYCLE
              DO L = 1, IND(NJ+1) - IND(NJ)
                O = IND(NJ)+L
                DO M = 1, IND(NI+1) - IND(NI)
                  N = IND(NI)+M
                  CALL DFTB_CNVSQ(N,O,NSEQ)
                  X(LSDER+NSEQ-1) = BUP(M,L)
                END DO
              END DO
            END DO
            CALL DFTB_PREP_QDER(NAT,L1,L2,IND,
     *        X(LQDER),D,X(LDTMP),S,X(LSDER),X(LWRK),X(LWRK1),
     *        X(LWRK2),X(LWRK3),UHF,X(LDDER+L2*(ISPIN-1)))
            CALL DFTB_CALC_QDER(1,NAT,L1,IND,NDOCC,NPART,
     *        X(LQDER),U,OCC,EIG,X(LQDER),X(LDNDR),DUM,DUM,VEC,FON,
     *        .FALSE.,NSPE,ISPE,MAXANG,NSHELL,X(LWRK),X(LWRK1),S,INDSAO)
C    *        X(LAST),X(LAST),FRQOMG)
          END IF
C
C         EXPLICIT ELECTRIC FIELD CONTRIBUTIONS
C
          IF (EFLDL.AND.NCP.EQ.1) THEN
            J = 3*(NI-1)+KI
            DO I = 1, NAT
              DO K = 1, 3
                HESSIAN(3*(I-1)+K,J)
     *            = HESSIAN(3*(I-1)+K,J) + EVEC(K)*X(LQDER+I-1)
              END DO
            END DO
          END IF
C
C         CALCULATE DIPOLE DERIVATIVE (TENSOR)
C
          IF (NCP.EQ.1) THEN
            DO I = 1, NAT
              DO J = 1, 3
                DIPDER(J,KI,NI) = DIPDER(J,KI,NI) - X(LQDER+I-1)*C(J,I)
              END DO
            END DO
            DIPDER(KI,KI,NI) = DIPDER(KI,KI,NI) - CHAMUL(NI)
            DO I = 1, 3
              DIPDER(I,KI,NI) = DIPDER(I,KI,NI) * DEBYE_ANGSTROM
            END DO
          END IF
C
C         CALCULATE POLARIZABILITY DERIVATIVE (d3E/dadFdG)
C
          IF (ARAMAN.AND.IAND(MODHSS,4).NE.0) THEN
C           CALL CPU_TIME(TAR0)
C           CALL SYSTEM_CLOCK(ICAR0,ICR,ICM)
            IF (DFTB3.AND.SRSCC)
     *        CALL DCOPY(NDIMSH,X(LQDER_INIT),1,X(LQDER),1)
            DO ISPIN = 1, NSPIN
              LOC = L3*(ISPIN-1)
              CALL DFTB_CPYTSQ(X(LHDERIJ+L2*(ISPIN-1)),
     *                         X(LFDERSC+LOC),L1,L1)
              CALL DFTB_CPYTSQ(X(LSDERIJ+L2*(ISPIN-1)),
     *                         X(LSDERSC+LOC),L1,L1)
              DO I = 1, L3
                X(LFDERSC+I-1+LOC)
     *            = X(LFDERSC+I-1+LOC)*X(LEIGINV2+I-1+LOC)
                X(LSDERSC+I-1+LOC)
     *            = X(LSDERSC+I-1+LOC)*X(LEIGINV2+I-1+LOC)
              END DO
             CALL DFTB_DEGCLR(L1,NDEG2,ISPIN,X(LLSTDEG2),X(LFDERSC+LOC))
             CALL DFTB_DEGCLR(L1,NDEG2,ISPIN,X(LLSTDEG2),X(LSDERSC+LOC))
            END DO
            CALL DFTB_CALC_POLDER(NI,KI,OCC,EIG,NAT,IND,L1,L2,NSPIN,
     *        NDOCC,NPART,C,
     *        FON,BETA,X(LDEDB),FEDER,X(LEDERF),X(LFDERF),
     *        FEDERTMP,ALPDER(1,KI,NI),X(LQDERF),X(LQDER),DISTMAT,
     *        GAMMA,
     *        GAMDER,X(LSDER),X(LDDERF),X(LWRK),X(LHDERIJ),X(LSDERIJ),
     * X(LFONDER),X(LFDERFSC),X(LFDERSC),X(LSDERSC),
     * NDEG2,NDEG3,X(LLSTDEG2),X(LLSTDEG3),DEGTHR(2),
     * X(LWRK1),X(LWRK2),X(LWRK3),
     * DFTB3,GAMMA3,GAMDER3,CHAMUL,MAXANG,ISPE,NDIMSH,INDSH,SRSCC,UHF,
     * SPNCST,NSHELL,NSPE,NGAMMA,CHAMULS,LCDFTB,S,GAMLC,GAMDERLC)
C           CALL SYSTEM_CLOCK(ICAR1,ICR,ICM)
C           CALL CPU_TIME(TAR1)
C           TAR  = TAR  + TAR1-TAR0
C           TARW = TARW + DBLE(ICAR1-ICAR0)/DBLE(ICR)
          END IF
C
CD        IF (FIRST) THEN
CD          CALL CPU_TIME(TIT1)
CD          CALL SYSTEM_CLOCK(ICLOCK_COUNTS_2,ICLOCK_RATE)
CD          FIRST = .FALSE.
CD          TT = TIT1-TIT0
CD          TTW = (ICLOCK_COUNTS_2-ICLOCK_COUNTS_1)/DBLE(ICLOCK_RATE)
CD          IF (MASWRK.AND.SOME) THEN
CD            IF (ISGDDI) THEN
CD              NP = NGROUPS
CD            ELSE
CD              NP = NPROC
CD            END IF
CD            WRITE (IW,'(/"   CALCULATION OF THE FIRST VECTOR TOOK
CD   * ",F8.2," CPU TIME (SECONDS)")') TT
CD            WRITE (IW,'("   ESTIMATED TIME FOR ",I4," VECTORS IS ",
CD   *          F10.2," SECONDS (",F7.2," HOURS)")')
CD   *          3*NAT,3*NAT*TT/DBLE(NP),3*NAT*TT/(3.6D+03*DBLE(NP))
CD            WRITE (IW,'( "   CALCULATION OF THE FIRST VECTOR TOOK
CD   * ",F7.2," WALL TIME (SECONDS)")') TTW
CD            WRITE (IW,'("   ESTIMATED TIME FOR ",I4," VECTORS IS ",
CD   *          F10.2," SECONDS (",F7.2," HOURS)")')
CD   *          3*NAT,3*NAT*TTW/DBLE(NP),3*NAT*TTW/(3.6D+03*DBLE(NP))
CD          END IF
CD        END IF
        END DO
        END DO
      END DO
C
      IF (ISGDDI) CALL GDDICOUNT( 1,LGROUP,MYJOB)
C
      IF (.NOT.MASWRK) THEN
        NITER = 0
        NITER2= 0
        IF (NFG.EQ.0) CALL VCLR(DIPDER,1,9*NAT)
      END IF
C
      IF (ISGDDI) CALL GDDI_SCOPE(DDI_WORLD)
C
      IF (PARALL3.OR.ISGDDI) THEN
        CALL DDI_GSUMI(1,NITER,1)
        CALL DDI_GSUMI(1,NITER2,1)
      END IF
C
      IF (MASWRK.AND.SOME) THEN
C       CALL CPU_TIME(T1)
C       CALL SYSTEM_CLOCK(IC1,ICR,ICM)
        IF (CPA) THEN
          WRITE (IW,*)
          IF (.NOT.CPA2) THEN
            WRITE (IW,'(X,I4," ROOTS CONVERGED AFTER ",F5.1,
     *        " CYCLES ON AVERAGE")') 3*NAT,DBLE(NITER)/DBLE(3*NAT)
          ELSE
            WRITE (IW,'(" FOR FIRST-ORDER COUPLED-PERTURBED:")')
            WRITE (IW,'(X,I4," ROOTS CONVERGED AFTER ",F5.1,
     *        " CYCLES ON AVERAGE")') 3*NAT,DBLE(NITER)/DBLE(3*NAT)
            WRITE (IW,'(" FOR SECOND-ORDER COUPLED-PERTURBED:")')
            WRITE (IW,'(X,I4," ROOTS CONVERGED AFTER ",F5.1,
     *       " CYCLES ON AVERAGE")')
     *       NVEC-3*NAT,DBLE(NITER2)/DBLE(NVEC-3*NAT)
          END IF
C         WRITE (IW,'(X,"SOLUTION OF CP-DFTB TOOK ",F10.2," (",F10.2,
C    *      ") CPU (WALL) SECONDS")')
C    *      T1-T0-TAR,DBLE(IC1-IC0)/DBLE(ICR)-TARW
        END IF
C       IF (ARAMAN.AND.IAND(MODHSS,4).NE.0) THEN
C         WRITE (IW,'(X,"POLARIZABILITY DER. TOOK ",F10.2," (",F10.2,
C    *      ") CPU (WALL) SECONDS")') TAR,TARW
C       END IF
        WRITE (IW,*)
      END IF
C
      IF (ARAMAN.OR.IAND(MODHSS,30).NE.0) THEN
        IF (IAND(MODHSS,4).NE.0) CALL DSCAL(18*NAT,-ONE,ALPDER,1)
C
C       ----- CALCULATE ALPHA POLARIZABILITY IF ANALYTIC RAMAN -----
C       d2E/dFdG = dq/dF * Cg
C
C       CALL CPU_TIME(TAR0)
C       CALL SYSTEM_CLOCK(ICAR0,ICR,ICM)
        IF (MASWRK) THEN
          IF (ARAMAN) THEN
            CALL DFTB_CALC_POLAR(POLAR,X(LQDERF),C,NAT,NDIMSH,NSPIN,
     *                           MAXANG,ISPE,SRSCC,PFREQ(1).NE.ZERO,
     *                           IAND(MODHSS,64).NE.0)
            TMP = POLAR(3)
            POLAR(3) = POLAR(4)
            POLAR(4) = TMP
          ELSE IF (CPF2) THEN
            CALL DFTB_CALC_POLAR(X(LWRK),X(LQDERF),C,NAT,NDIMSH,NSPIN,
     *                           MAXANG,ISPE,SRSCC,PFREQ(1).NE.ZERO,
     *                           IAND(MODHSS,64).NE.0)
            WRITE (IW,'(" POLARIZABILITY FOR X-X =",F25.10)') X(LWRK+0)
            WRITE (IW,'(" POLARIZABILITY FOR X-Y =",F25.10)') X(LWRK+1)
            WRITE (IW,'(" POLARIZABILITY FOR X-Z =",F25.10)') X(LWRK+2)
            WRITE (IW,'(" POLARIZABILITY FOR Y-Y =",F25.10)') X(LWRK+3)
            WRITE (IW,'(" POLARIZABILITY FOR Y-Z =",F25.10)') X(LWRK+4)
            WRITE (IW,'(" POLARIZABILITY FOR Z-Z =",F25.10)') X(LWRK+5)
          END IF
        END IF
C
C       ----- ALSO, CALCULATE HYPER-POLARIZABILITY -----
C       d3E/dFdGdH
        IF (IAND(MODHSS,2).NE.0)
     *  CALL DFTB_CALC_HYPER(OCC,EIG,IW,NAT,L1,L2,NSPIN,NDOCC,NPART,FON,
     *                       BETA,X(LEDERF),X(LFDERF),FEDERTMP,
     *                       X(LFDERFSC),NDEG2,NDEG3,
     *                       X(LLSTDEG2),X(LLSTDEG3),X(LWRK),X(LWRK1),
     *                       X(LWRK2),DEGTHR(2),DFTB3,
     *                       GAMMA3,X(LQDERF),NDIMSH,SRSCC,INDSH,
     *                       ISPE,MAXANG,PFREQ(1).NE.ZERO,
     *                       IAND(MODHSS,64).NE.0)
C       CALL SYSTEM_CLOCK(ICAR1,ICR,ICM)
C       CALL CPU_TIME(TAR1)
C       TAR  = TAR1-TAR0
C       TARW = DBLE(ICAR1-ICAR0)/DBLE(ICR)
C       IF (MASWRK.AND.SOME.AND.(CPF.OR.CPF2)) THEN
C         WRITE (IW,'(X,"POLARIZABILITY      TOOK ",F10.2," (",F10.2,
C    *      ") CPU (WALL) SECONDS")') TAR,TARW
C       END IF
C
C       ----- CALCULATE d3E/dadbdF AND d3E/dadbdc -----
C
        IF (IAND(MODHSS,16).NE.0) THEN
C         CALL CPU_TIME(TAR0)
C         CALL SYSTEM_CLOCK(ICAR0,ICR,ICM)
C         CALL DFTB_FULL3RD
C         CALL DFTB_FULL3RD(IW,L1,L2,L3,NAT,NA,NB,NSPE,NDIMSH,NGAMMA,
C    *                      NSHELL,NSPIN,NDOCC,NPART,NOCC,IND,INDSH,
C    *                      ISPE,MAXANG,NDEG2,NDEG3,X(LLSTDEG2),
C    *                      X(LLSTDEG3),MAXBND,INDBND,
C    *                      C,OCC,EIG,VEC,S,D,DISTMAT,CHAMUL,CHAMULS,
C    *                      GAMMA,GAMMA3,GAMDER,GAMDER3,
C    *                      HUBBL,HUBDER,DAMPXHE,SPE,BETA,SHIFT,
C    *                      DEGTHR(2),
C    *                      X(LEDERF),X(LFDERF),FEDERTMP,X(LFDERFSC),
C    *                      X(LDSAVE),X(LFSAVE),X(LSSAVE),X(LFESAVE),
C    *                      X(LWRK),X(LWRK1),X(LWRK2),X(LWRK3),efermi,
C    *                      FON,DFTB3,SRSCC,DAMPXH)
C         CALL SYSTEM_CLOCK(ICAR1,ICR,ICM)
C         CALL CPU_TIME(TAR1)
C         TAR  = TAR1-TAR0
C         TARW = DBLE(ICAR1-ICAR0)/DBLE(ICR)
C         IF (MASWRK.AND.SOME) THEN
C           WRITE (IW,'(X,"FULL 3RD DERIVATIVE TOOK ",F10.2," (",F10.2,
C    *        ") CPU (WALL) SECONDS")') TAR,TARW
C         END IF
        END IF
      END IF
C
C     IF (PFREQ(1).NE.ZERO.AND.CPF.AND.MASWRK.AND.SOME
C    *  .AND.IAND(MODHSS,2+32).NE.0) THEN
C       CALL DFTB_FDERPRINT(NCPF,X(LFDER2),X(LFDER3),X(LFDER4),
C    *                      PFREQ(1))
C     END IF
C
      IF (NEEDAR.NE.0) CALL RETFM(NEEDAR)
C
C     ----- RELEASE ALLOCATED MEMORY -----
C
1000  CONTINUE
      IF (SAVDCQ) CALL RETFM(NEED3)
      CALL RETFM(NEED)
C
      IF (PARALL3.OR.ISGDDI) THEN
C       GOPARR = .TRUE.
C       CALL DDI_GSUMF(1,HESSIAN,9*NAT*NAT)
        CALL DDI_GSUMF(1,DIPDER,9*NAT)
        IF (ARAMAN) THEN
          IF (PARALL3) THEN
            CALL GDDI_SCOPE(DDI_GROUP)
            CALL DSCAL(18*NAT,1.0D+00/DBLE(NPROC),ALPDER,1)
            CALL GDDI_SCOPE(DDI_WORLD)
          END IF
          CALL DDI_GSUMF(1,ALPDER,18*NAT)
        END IF
      END IF
C
      RETURN
C
      END SUBROUTINE DFTB_CP
C
C-----------------------------------------------------------------------
C*MODULE DFTBHS    *DECK DFTB_CALC_DWDER
C>
C>    @brief   Calculate quasi-degenerated coefficient
C>
C>    @details Compute binomial coefficient (nCm)
C>
C>    @author  Yoshio Nishimoto
C>    - Feb, 2018- Subroutine written
C>
C>           --- INPUT ---
C>    @param NCP     Order of coupled-perturbed
C>    @param L1      Number of basis functions
C>    @param L2      L1*(L1+1)/2
C>    @param NSPIN   Number of spins (1=R-DFTB,2=U-DFTB)
C>    @param NDOCC   Number of doubly occupied orbitals
C>    @param NPART   Number of partially occupied orbitals
C>    @param EIG     Eigenvalue of each MO
C>    @param OCC     Occupation number of each MO
C>    @parma VEC     MO coefficients
C>    @parma EDER    Eigenvalue derivative
C>    @parma FONDER  FON derivative
C>    @parma FDER    Fock matrix derivative in MO
C>    @parma U       Rotation of MOs
C>    @param WRK     Work space
C>    @param WRK1    Work space
C>    @param WRK2    Work space
C>    @param SDERIJ  Overlap derivative in MO
C>    @param FON     Whether fractional occupation is used or not
C>    @param SCC     SCC or not
C>    @param BETA    1/kTe
C>    @param FEDER   Fermi energy derivative
C>           --- OUTPUT ---
C>    @parma DDER    Density matrix derivative
C>    @parma WDER    Energy-weighted density matrix derivative
C>
C
C
      SUBROUTINE DFTB_CALC_DWDER(NCP,L1,L2,NSPIN,NDOCC,NPART,EIG,OCC,
     *                           VEC,EDER,FONDER,FDER,U,DDER,WDER,WRK,
     *                           WRK1,WRK2,SDERIJ,FON,SCC,BETA,FEDER)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      COMMON /DFTBHS/ DEGTHR(2),PFREQ(2),CPCONV,OCCTHR,MODHSS,
     *                MODPRJ,MXCPIT,ARAMAN
C
      PARAMETER (ZERO=0.0D+00,ONE=1.0D+00,TWO=2.0D+00,HALF=0.5D+00)
C
      DIMENSION NDOCC(2),NPART(2),EIG(L1,*),OCC(L1,*),VEC(L1,L1,*),
     *          EDER(L1,*),FONDER(L1,*),FDER(L2,*),U(L1,L1,*),
     *          DDER(L2,*),WDER(L2,*),WRK(L1,*),WRK1(L1,*),WRK2(L1,*),
     *          SDERIJ(L2,*)
      LOGICAL   ARAMAN,FON(3),SCC
C
C     L0=L1
      L3=L1*L1
      DO ISPIN = 1, NSPIN
        NOCC = NDOCC(ISPIN) + NPART(ISPIN)
        !! independent pairs
        IF (NDOCC(ISPIN).GE.1) THEN
          CALL DFTB_CPYTSQ(SDERIJ(1,ISPIN),WRK1,L1,L1)
C
          !! independent pairs for W derivative
          CALL VCLR(WRK2,1,L3)
          DO IO = 1, NDOCC(ISPIN)
            DO JO = 1, NDOCC(ISPIN)
              CALL DFTB_CNVSQ(IO,JO,NSEQ)
              WRK2(IO,JO) = FDER(NSEQ,ISPIN)
     *          -(EIG(IO,ISPIN)+EIG(JO,ISPIN))*WRK1(IO,JO)
              IF (IO.EQ.JO) WRK2(IO,JO) = -EIG(IO,ISPIN)*WRK1(IO,JO)
            END DO
          END DO
          CALL DGEMM('N','N',L1,NDOCC(ISPIN),NDOCC(ISPIN),ONE,
     *               VEC(1,1,ISPIN),L1,WRK2,L1,ZERO,WRK,L1)
          CALL DGEMM('N','T',L1,L1,NDOCC(ISPIN),ONE,WRK,L1,
     *               VEC(1,1,ISPIN),L1,ZERO,WRK2,L1)
          CALL SYMSQT(L1,WRK2,L1,WDER(1,ISPIN),1)
          IF (NSPIN.EQ.1) CALL DSCAL(L2,TWO,WDER(1,ISPIN),1)
          IF (NSPIN.EQ.2) CALL DSCAL(L2,ONE,WDER(1,ISPIN),1)
C
          !! independent pairs for D derivative
          IF (.NOT.SCC) THEN
            DO IO = 1, NDOCC(ISPIN)
              CALL DSCAL(NDOCC(ISPIN),OCC(IO,ISPIN),WRK1(1,IO),1)
            END DO
            CALL DGEMM('N','N',L1,NDOCC(ISPIN),NDOCC(ISPIN),ONE,
     *                 VEC(1,1,ISPIN),L1,WRK1,L1,ZERO,WRK2,L1)
            CALL DGEMM('N','T',L1,L1,NDOCC(ISPIN),ONE,WRK2,L1,
     *                 VEC(1,1,ISPIN),L1,ZERO,WRK,L1)
            CALL SYMSQT(L1,WRK,L1,DDER(1,ISPIN),1)
            CALL DSCAL(L2,-ONE,DDER(1,ISPIN),1)
          END IF
        ELSE
          IF (.NOT.SCC) CALL VCLR(DDER(1,ISPIN),1,L2)
          CALL VCLR(WDER(1,ISPIN),1,L2)
        END IF
C
        !! C*U
        IF (FON(ISPIN)) THEN
          CALL DCOPY(L3,U(1,1,ISPIN),1,WRK,1)
          CALL DGEMM('N','N',L1,NOCC,L1-NDOCC(ISPIN),TWO,
     *               VEC(1,NDOCC(ISPIN)+1,ISPIN),L1,
     *               WRK(NDOCC(ISPIN)+1,1),L1,ZERO,WRK1,L1)
          CALL DSCAL(L1*NDOCC(ISPIN),OCC(1,ISPIN),WRK1,1)
          DO IO = NDOCC(ISPIN)+1, NOCC
            CALL DSCAL(L1,OCC(IO,ISPIN),WRK1(1,IO),1)
            IF (NCP.EQ.1) THEN
              CALL DAXPY(L1,FONDER(IO,ISPIN),VEC(1,IO,ISPIN),1,
     *                   WRK1(1,IO),1)
            ELSE IF (NCP.GE.2) THEN
              OCCI = OCC(IO,ISPIN)*HALF
              SCAL = -OCCI*(ONE-OCCI)*BETA
              NSEQ = IO*(IO+1)/2
              EI = EIG(IO,ISPIN)
C             CALL DAXPY(L1,SCAL*(FDER(NSEQ,ISPIN)
C    *                      -EI*SDERIJ(NSEQ,ISPIN)-FEDER)*TWO,
C    *                   VEC(1,IO,ISPIN),1,WRK1(1,IO),1)
              CALL DAXPY(L1,SCAL*(FDER(NSEQ,ISPIN)-FEDER)*TWO,
     *                   VEC(1,IO,ISPIN),1,WRK1(1,IO),1)
            END IF
          END DO
        ELSE
          CALL DGEMM('N','N',L1,NOCC,L1-NDOCC(ISPIN),TWO,
     *               VEC(1,NDOCC(ISPIN)+1,ISPIN),L1,
     *               U(NDOCC(ISPIN)+1,1,ISPIN),L1,ZERO,WRK1,L1)
          CALL DSCAL(L1*NDOCC(ISPIN),OCC(1,ISPIN),WRK1,1)
        END IF
        !! non-independent pairs + df/dR for D derivative
        CALL DGEMM('N','T',L1,L1,NOCC,ONE,WRK1,L1,
     *             VEC(1,1,ISPIN),L1,ZERO,WRK2,L1)
        CALL SYMSQT(L1,WRK2,L1,WRK,1)
        CALL DAXPY(L2,ONE,WRK,1,DDER(1,ISPIN),1)
C
        !! non-independent pairs + f*de/dR + df/dR*e for W derivative
        IF (FON(ISPIN)) THEN
          CALL DCOPY(L3,U(1,1,ISPIN),1,WRK,1)
          DO IO = 1, NDOCC(ISPIN)
            CALL VCLR(WRK(1,IO),1,NDOCC(ISPIN))
          END DO
          DO IO = NOCC+1, L1
            CALL VCLR(WRK(NOCC+1,IO),1,L1-NOCC)
          END DO
C
          DO IO = 1, NOCC
            EI = EIG(IO,ISPIN)
            OCCI = OCC(IO,ISPIN)
            DO JO = NDOCC(ISPIN)+1, NOCC
              IF (IO.GE.JO) CYCLE
              EJ = EIG(JO,ISPIN)
              OCCJ = OCC(JO,ISPIN)
              !! simply use U, if io~=jo
              IF (ABS(EI-EJ).GE.DEGTHR(1)) THEN
                CALL DFTB_CNVSQ(IO,JO,NSEQ)
                SCAL = (OCCJ*EJ-OCCI*EI)/((EJ-EI)*OCCI*EI)
                WRK(JO,IO) = -SDERIJ(NSEQ,ISPIN)
     *            + SCAL*(FDER(NSEQ,ISPIN)-EJ*SDERIJ(NSEQ,ISPIN))
              ELSE
                CALL DFTB_CNVSQ(IO,JO,NSEQ)
                !! use clcdeg
C               SCAL = (OCCJ*EJ-OCCI*EI)/((EJ-EI)*OCCI*EI)
                SCAL = ONE/EI - (ONE-OCCJ*HALF)*BETA
                WRK(JO,IO) = -SDERIJ(NSEQ,ISPIN)
     *            + SCAL*(FDER(NSEQ,ISPIN)-EJ*SDERIJ(NSEQ,ISPIN))
              END IF
            END DO
          END DO
          CALL DGEMM('N','N',L1,NOCC,L1-NDOCC(ISPIN),TWO,
     *               VEC(1,NDOCC(ISPIN)+1,ISPIN),L1,
     *               WRK(NDOCC(ISPIN)+1,1),L1,ZERO,WRK1,L1)
          DO IO = 1, NOCC
            CALL DSCAL(L1,OCC(IO,ISPIN)*EIG(IO,ISPIN),WRK1(1,IO),1)
            SCAL = EDER(IO,ISPIN)*OCC(IO,ISPIN)
     *           + FONDER(IO,ISPIN)*EIG(IO,ISPIN)
            CALL DAXPY(L1,SCAL,VEC(1,IO,ISPIN),1,WRK1(1,IO),1)
          END DO
        ELSE
          DO IO = 1, NDOCC(ISPIN)
            CALL DSCAL(L1,EIG(IO,ISPIN),WRK1(1,IO),1)
            SCAL = EDER(IO,ISPIN)*OCC(IO,ISPIN)
            CALL DAXPY(L1,SCAL,VEC(1,IO,ISPIN),1,WRK1(1,IO),1)
          END DO
        END IF
        CALL DGEMM('N','T',L1,L1,NOCC,ONE,WRK1,L1,
     *             VEC(1,1,ISPIN),L1,ZERO,WRK2,L1)
        CALL SYMSQT(L1,WRK2,L1,WRK,1)
        CALL DAXPY(L2,ONE,WRK,1,WDER(1,ISPIN),1)
      END DO
C     call cpu_time(t1)
C     write (*,*) "cpu time = ", t1-t0
C     write (*,*) "dder"
C     call prtril(dder(1,1),l1)
C     write (*,*) "wder"
C     call prtril(wder(1,1),l1)
C     if (nspin.eq.2) then
C       write (*,*) "dder(beta)"
C       call prtril(dder(1,2),l1)
C       write (*,*) "wder(beta)"
C       call prtril(wder(1,2),l1)
C     end if
C
      RETURN
C
      END SUBROUTINE DFTB_CALC_DWDER
C
C-----------------------------------------------------------------------
C*MODULE DFTBHS    *DECK DFTB_CALC_FDER_DEDB
C>
C>    @brief   DFTB Hessian
C>
C>    @details Calculate F_{ij}^a and de_i/da
C>             Hessian
C>
C>    @author Yoshio Nishimoto
C>    - Feb, 2016- Subroutine written
C>
C>           --- INPUT ---
C>    @param L1        Number of AOs
C>    @param L2        (L1*L1+L1)/2
C>    @param NSPIN     Number of spins (1: RHF, 2:UHF)
C>    @param NDIMSH    If (SRSCC); NDIMSH=NSHELL; else NDIMSH=NAT
C>    @param NDIM      Number of occupied orbitals?
C>    @param HDERMO    H_{ij}^{0,a}
C>    @param SDERMO    S_{ij}^a
C>    @param EIG       e_i
C>           --- OUTPUT ---
C>    @param FDERMO    F_{ij}^a
C>    @param DEDB      de_i/da
C>
C
      SUBROUTINE DFTB_CALC_FDER_DEDB(L1,L2,NSPIN,FDERMO,DEDB,HDERMO,
     *  SDERMO,EIG)
C
      IMPLICIT NONE
C
      DOUBLE PRECISION, INTENT(IN) :: HDERMO(L2,NSPIN),SDERMO(L2,NSPIN),
     *  EIG(L1,NSPIN)
      DOUBLE PRECISION, INTENT(OUT) :: FDERMO(L2,NSPIN),DEDB(L1,NSPIN)
      INTEGER, INTENT(IN) :: L1,L2,NSPIN
C
      INTEGER :: ISPIN,NSEQ,IO,JO
      DOUBLE PRECISION :: EI,EJ
C
      common /fmoinf/ nfg,nlayer,natfmo,nbdfg,naotyp,nbody,nsegm
      INTEGER :: nfg,nlayer,natfmo,nbdfg,naotyp,nbody,nsegm
C
      DOUBLE PRECISION, PARAMETER :: TWO=2.0D+00
C      
C     CALCULATE DEDB (DERIVATIVE OF EIGENVALUES) AND FDERIJ TERMS, WHERE
C     FDERIJ_{IJ} = \SUM_{\MU \NU} C_{\MU I} C_{\NU J} dF_{\MU \NU}/dB
C     HERE WRITES FDERMO.
C     "F" IS A KIND OF FOCK MATRIX (SIMPLY HAMILTONIAN IN DFTB)
C
      DO ISPIN = 1, NSPIN
        NSEQ = 0
        IF (NFG.EQ.0) THEN
          DO IO = 1, L1
            EI = EIG(IO,ISPIN)
            DO JO = 1, IO-1
              EJ = EIG(JO,ISPIN)
              NSEQ = NSEQ + 1
CORIG         FDERMO(NSEQ,ISPIN)
CORIG*          = HDERMO(NSEQ,ISPIN) - (EI+EJ)*SDERMO(NSEQ,ISPIN)
              FDERMO(NSEQ,ISPIN)
     *          = HDERMO(NSEQ,ISPIN)-(EI+EJ)*SDERMO(NSEQ,ISPIN)*0.5d+00
            END DO
            NSEQ = NSEQ + 1
CORIG       FDERMO(NSEQ,ISPIN) = -EI*SDERMO(NSEQ,ISPIN)
           FDERMO(NSEQ,ISPIN) = HDERMO(NSEQ,ISPIN)-EI*SDERMO(NSEQ,ISPIN)
C           FDERMO(NSEQ,ISPIN) = -SDERMO(NSEQ,ISPIN)*0.5d+00
C           FDERMO(NSEQ,ISPIN) = 0.0d+00
C           FDERMO(NSEQ,ISPIN) = HDERMO(NSEQ,ISPIN)
            DEDB(IO,ISPIN)
     *        = HDERMO(NSEQ,ISPIN)-EI*SDERMO(NSEQ,ISPIN)
          END DO
        ELSE
          DO IO = 1, L1
            EI = EIG(IO,ISPIN)
            DO JO = 1, IO-1
              EJ = EIG(JO,ISPIN)
              NSEQ = NSEQ + 1
              FDERMO(NSEQ,ISPIN)
     *          = TWO*HDERMO(NSEQ,ISPIN)-(EI+EJ)*SDERMO(NSEQ,ISPIN)
            END DO
            NSEQ = NSEQ + 1
            FDERMO(NSEQ,ISPIN)
     *        = HDERMO(NSEQ,ISPIN)-EI*SDERMO(NSEQ,ISPIN)
C
            DEDB(IO,ISPIN)
     *        = HDERMO(NSEQ,ISPIN)-EI*SDERMO(NSEQ,ISPIN)
          END DO
        END IF
      END DO
C
      RETURN
C 
      END SUBROUTINE DFTB_CALC_FDER_DEDB
C
C-----------------------------------------------------------------------
C*MODULE DFTBHS    *DECK DFTB_HESS_UUPDATE
C>
C>    @brief Hessian of DFTB
C>
C>    @details Update U matrix, by calculating response terms
C>
C>    @author Yoshio Nishimoto
C>    - Sep, 2014- Subroutine written
C>    @date Feb, 2016 - Yoshio Nishimoto
C>    - Avoid degenerate problem
C>
C>           --- INPUT ---
C>    @param UINIT   Initial U matrix (only with non-response terms)
C>    @param L1      Number of basis functions
C>    @param NDOCC   Number of doubly occupied orbitals
C>    @param NPART   Number of partially occupied orbitals
C>    @param MAXOCC  Number of occupied + partially occupied orbitals
C>    @param MAXVIR  Number of partially occupied + virtual orbitals
C>    @param EIGINV  1/(ei-ej)
C>           --- OUTPUT ---
C>    @param U       U matrix (derivatives of MO coefficient)
C>
C
      SUBROUTINE DFTB_HESS_UUPDATE(U,UINIT,SDSHIJ,L1,
     *  NDOCC,NPART,MAXOCC,MAXVIR,EIGINV)
C
      IMPLICIT NONE
C
      INTEGER, INTENT(IN) :: L1,NDOCC,NPART,MAXOCC,MAXVIR
      DOUBLE PRECISION, INTENT(INOUT) :: U(L1,L1)
      DOUBLE PRECISION, INTENT(IN) :: UINIT(L1,L1),SDSHIJ(L1,*),
     *  EIGINV(MAXVIR,MAXOCC)
C
      INTEGER :: IO,JO
C
      DO IO = 1, NDOCC+NPART
        DO JO = NDOCC+1, L1
          IF (IO.GE.JO) CYCLE
          U(JO,IO) = UINIT(JO,IO) + EIGINV(JO-NDOCC,IO)*SDSHIJ(JO,IO)
        END DO
      END DO
C
      RETURN
C
      END SUBROUTINE DFTB_HESS_UUPDATE
C
C-----------------------------------------------------------------------
C*MODULE DFTBHS    *DECK DFTB_DPART
C>
C>    @brief Hessian of DFTB
C>
C>    @details Calculate partial density matrix.  The index of "i" runs
C>             over only doubly occupied orbitals.
C>
C>    @author Yoshio Nishimoto
C>    - Sep, 2014- Subroutine written
C>
C>           --- INPUT ---
C>    @param L1 Number of basis functions
C>    @param L2 L1*(L1+1)/2
C>    @param NDOCC Number of doubly occupied orbitals
C>    @param OCC Occupation number of each MO
C>    @param VEC MO coefficients
C>           --- OUTPUT ---
C>    @param D Partial density matrix
C>
C
      SUBROUTINE DFTB_DPART(L1,L2,NDOCC,OCC,D,VEC)
C
      IMPLICIT NONE
C
      INTEGER, INTENT(IN) :: L1,L2,NDOCC
      DOUBLE PRECISION ,INTENT(IN) :: OCC(L1),VEC(L1,L1)
      DOUBLE PRECISION, INTENT(OUT) :: D(L2)
C
      INTEGER :: I,IJ,K
C
      CALL VCLR(D,1,L2)
C
      DO K = 1, NDOCC
        IJ = 1
        DO I = 1, L1
          CALL DAXPY(I,OCC(K)*VEC(I,K),VEC(1,K),1,D(IJ),1)
          IJ = IJ + I
        END DO
      END DO
C
      RETURN
C
      END SUBROUTINE DFTB_DPART
C
C-----------------------------------------------------------------------
C*MODULE DFTBHS    *DECK DFTB_CALC_DNDR
C>
C>    @brief Hessian of DFTB
C>
C>    @details Calculate derivatives of occupation numbers.
C>
C>    @author Yoshio Nishimoto
C>    - Sep, 2014- Subroutine written
C>
C>           --- INPUT ---
C>    @param IMODE See below
C>    @param L1 Number of basis functions
C>    @param L2 L1*(L1+1)/2
C>    @param NDOCC Number of doubly occupied orbitals
C>    @param NPART Number of partially occupied orbitals
C>    @param ETEMP Electronic temperature
C>    @param OCC Occupatio number of each MO
C>    @param EIG Eigenvalues of each MO
C>    @param HDER H_{ij}^b + V_{ij}^b
C>    @param SDER S_{ij}^b
C>    @param SOMGD Either non-response or response contribution of
C>           the derivative of shift contribution
C>    @param VEC MO coefficients
C>           --- OUTPUT ---
C>    @parma DEDB_INIT Non-response contributions of de/db (derivatives
C>           of eigenvalue)
C>    @param DEDB Derivatives of eigenvalue
C>    @param DNDR Derivatiees of occupation number
C>    @param FEDER Derivative of fermi-level
C>
C
      SUBROUTINE DFTB_CALC_DNDR(IMODE,L1,L2,NDOCC,NPART,ETEMP,
     *  DEDB_INIT,DEDB,DNDR,OCC,EIG,HDER,SDER,SOMGD,FEDER,SCC)
C
      IMPLICIT NONE
C
      INTEGER, INTENT(IN) :: IMODE,L1,L2,NDOCC,NPART
      DOUBLE PRECISION, INTENT(IN) :: ETEMP,OCC(L1),EIG(L1),HDER(L2),
     *  SDER(L2),SOMGD(L1,L1)
      DOUBLE PRECISION, INTENT(INOUT) :: DEDB_INIT(L1),DEDB(L1),
     *  DNDR(L1),FEDER
      LOGICAL, INTENT(IN) :: SCC
C
      INTEGER :: NSEQ,IO
      DOUBLE PRECISION :: EI,OCCI,TMP,BETA
      DOUBLE PRECISION,PARAMETER :: ZERO=0.0D+00,HALF=0.5D+00,
     *  ONE=1.0D+00,TWO=2.0D+00,BOLTZ=3.1668114D-06
C
C     IMODE = 0
C     USED FOR FON CALCULATION.  JUST ADD NON-RESPONSE CONTRIBUTIONS TO
C     DEDB_INIT MATRIX.
C
C     IMODE = 1
C     USED FOR FON CALCULATION, AND ADD RESPONSE CONTRIBUTIONS
C
      NSEQ = 0
      IF (IMODE.EQ.0) THEN
C       WRITE (*,*) "DEDB_INIT"
        DO IO = 1, L1
          NSEQ = NSEQ + IO
          EI = EIG(IO)
C         WRITE (*,'(I3,2(X,F20.10))') NSEQ,HDER(NSEQ),SDER(NSEQ)
C         WRITE (*,'(3X,2(X,F20.10))') EI,SOMGD(NSEQ)
          TMP = HDER(NSEQ) - EI*SDER(NSEQ)
          DEDB_INIT(IO) = TMP
C         WRITE (*,'(I3,X,F20.10)') IO,DEDB_INIT(IO)
        END DO
        CALL VCLR(DNDR,1,L1)
        RETURN
      ELSE
        DO IO = 1, L1
          IF (SCC) THEN
            TMP = SOMGD(IO,IO)
          ELSE
            NSEQ = NSEQ + IO
            EI = EIG(IO)
            TMP = HDER(NSEQ) - EI*SDER(NSEQ)
          END IF
          DEDB(IO) = DEDB_INIT(IO) + TMP
        END DO
      END IF
C
C       ----- CALCULATE DERVIATIVE OF FERMI ENERGY FIRST -----
C
      CALL VCLR(DNDR,1,L1)
      FEDER = ZERO
      TMP = ZERO
      DO IO = NDOCC+1, NDOCC+NPART
        OCCI = OCC(IO)*HALF
        FEDER = FEDER + OCCI*(ONE-OCCI)*DEDB(IO)
        TMP = TMP + OCCI*(ONE-OCCI)
      END DO
      FEDER = FEDER/TMP
      BETA = ONE/(BOLTZ*ETEMP)
C
C       ----- THEN CALCULATE DERIVATIVE OF OCCUPATION NUMBERS -----
C
      DO IO = NDOCC+1, NDOCC+NPART
        OCCI = OCC(IO)*HALF
        DNDR(IO) = -OCCI*(ONE-OCCI)*(DEDB(IO)-FEDER) !! *BETA
      END DO
      CALL DSCAL(NPART,TWO*BETA,DNDR(NDOCC+1),1)
C
      END SUBROUTINE DFTB_CALC_DNDR
C
C-----------------------------------------------------------------------
C*MODULE DFTBHS    *DECK DFTB_PREP_QDER
C>
C>    @brief Hessian of DFTB
C>
C>    @details Calculate non-response contributions to the derivatives
C>             of Mulliken population
C>
C>    @author Yoshio Nishimoto
C>    - Sep, 2014- Subroutine written
C>    @date Feb, 2018 - Yoshio Nishimoto
C>    - Modified to use DGEMM and give initial D derivative
C>
C>           --- INPUT ---
C>    @param NAT       Number of atoms
C>    @param L1        Number of basis functions
C>    @param L2        L1*(L1+1)/2
C>    @param IND       Index of AO
C>    @param D1        Normal density matrix
C>    @param D2        Partial density matrix (for FON)
C>    @param S         Overlap matrix
C>    @param SDER      Matrix of overlap matrix in AO basis
C>    @param UHF Whether U-DFTB or not
C>           --- OUTPUT ---
C>    @param QDER_INIT Non-response contribution to Mulliken derivative
C>    @param WRK,WRK1,W,TMP Working array
C>    @param DDER      Initial density matrix derivative
C>
C
      SUBROUTINE DFTB_PREP_QDER(NAT,L1,L2,IND,QDER_INIT,
     *  D1,D2,S,SDER,WRK,WRK1,W,TMP,UHF,DDER)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      INTEGER :: IND(*)
      DOUBLE PRECISION :: QDER_INIT(NAT),D1(L2),D2(L2),S(L2),SDER(L2),
     *  WRK(*),WRK1(*),W(*),TMP(L1),DDER(L2)
      LOGICAL :: UHF
      DOUBLE PRECISION,PARAMETER :: ZERO=0.0D+00,HALF=0.5D+00,
     *  ONE=1.0D+00
C
C     D1 :: NORMAL DENSITY MATRIX
C     D2 :: PARTIAL DENSITY MATRIX (IF FON)
C
      !! CALCULATE LAGRANGIAN-LIKE TERM
      CALL DFTB_CPYTSQ(D2,WRK,L1,L1)
      CALL DFTB_CPYTSQ(SDER,WRK1,L1,L1)
      CALL DGEMM('T','N',L1,L1,L1,ONE,WRK,L1,WRK1,L1,ZERO,W   ,L1)
      CALL DGEMM('N','N',L1,L1,L1,ONE,W  ,L1,WRK ,L1,ZERO,WRK1,L1)
      CALL SYMSQT(L1,WRK1,L1,W,1)
      IF (.NOT.UHF) CALL DSCAL(L2,HALF,W,1)
      CALL DCOPY(L2,W,1,DDER,1)
      CALL DSCAL(L2,-ONE,DDER,1)
C
C     CALCULATE S*W (BOTH MATRICES ARE SYMMETRIC)
C
      CALL DFTB_MULLIKEN(L1,L2,S,W,TMP)
C
C     CALCULATE D*(dS/dR), BY DESTROYING W, WHICH IS NO LONGER NEEDED.
C     ONLY DIAGONAL ELEMENT IS NEEDED, AND USE THE NATURE OF
C     SYMMETRY OF D AND dS/dR.
C     
      CALL DFTB_MULLIKEN(L1,L2,D1,SDER,W)
      CALL VSUB(TMP,1,W,1,TMP,1,L1)
C
      DO NC = 1, NAT
        VAL = ZERO
        DO L = 1, IND(NC+1) - IND(NC)
          MU = IND(NC)+L
C         VAL = VAL + TMP1(MU) - TMP2(MU)
          VAL = VAL + TMP(MU)
        END DO
        QDER_INIT(NC) = VAL
      END DO
C
      RETURN
C
      END SUBROUTINE DFTB_PREP_QDER
C
C-----------------------------------------------------------------------
C*MODULE DFTBHS    *DECK DFTB_CALC_QDER
C>
C>    @brief Hessian of DFTB
C>
C>    @details Calculate non-response contributions to the derivatives
C>             of Mulliken population.
C>
C>    @author Yoshio Nishimoto
C>    - Sep, 2014- Subroutine written
C>    @data Oct, 2014 - Yoshio Nishimoto
C>    - Performance improvement
C>    @data Feb, 2018 - Yoshio Nishimoto
C>    - Performance improvement
C>
C>           --- INPUT ---
C>    @param NCP       Order of CP
C>    @param NAT       Number of atoms
C>    @param L1        Number of basis functions
C>    @param IND       Index of AO
C>    @param NDOCC     Number of doubly occupied orbitals
C>    @param NPART     Number of partially occupied orbitals
C>    @param U         U matrix (derivative of MO coefficients)
C>    @param OCC       Occupation number of each MO
C>    @param EIG       Eigenvalue of each MO
C>    @param QDER_INIT Non-response Mulliken derivative
C>    @param DNDR      Derivative of occupation numbers
C>    @param FEDER     Derivative of Fermi energy
C>    @param BETA      1/kTe
C>    @param VEC       MO coefficients
C>    @param FON       Whether fractional occupation or not
C>    @param SRSCC     Whether Shell-resolved SCC or not
C>    @param NSPE      Number of species
C>    @param ISPE      Index of species
C>    @param MAXANG    Maximum angular momentum of each species
C>    @param NSHELL    Number of shells
C>    @param WRK,WRK1  Working arrays
C>    @param S         Overlap matrix in AO
C>    @param INDSAO    AO index of shell
C>           --- OUTPUT ---
C>    @param QDER      Derivative of Mulliken charges
C>
C
      SUBROUTINE DFTB_CALC_QDER(NCP,NAT,L1,IND,NDOCC,NPART,
     *  QDER,U,OCC,EIG,QDER_INIT,DNDR,FEDER,BETA,VEC,FON,SRSCC,
     *  NSPE,ISPE,MAXANG,NSHELL,WRK,WRK1,S,INDSAO)
C    *  FDERIJ,SDERIJ,FRQOMG)
C
      IMPLICIT NONE
C
      INTEGER, INTENT(IN) :: NCP,NAT,L1,IND(*),NDOCC,NPART,
     *  NSPE,ISPE(NAT),MAXANG(NSPE),NSHELL,INDSAO(*)
      DOUBLE PRECISION, INTENT(IN)  :: U(L1,L1),OCC(L1),EIG(L1),
     *  QDER_INIT(*),VEC(L1,L1),DNDR(L1),S(*),FEDER,BETA
C    *  FDERIJ(*),SDERIJ(*),FRQOMG
      DOUBLE PRECISION, INTENT(OUT) :: QDER(*)
      DOUBLE PRECISION, INTENT(INOUT) :: WRK(L1,L1),WRK1(L1,L1)
      LOGICAL, INTENT(IN) :: FON,SRSCC
C
      INTEGER :: ISH0,ISH,I,IO,MM,NN,MU,NOCC,L2
      DOUBLE PRECISION :: VAL
      DOUBLE PRECISION,PARAMETER :: ZERO=0.0D+00,ONE=1.0D+00,TWO=2.0D+00
C
      L2 = L1*(L1+1)/2
      NOCC = NDOCC + NPART
C
C     CALCULATE DENSITY MATRIX DERIVATIVE DUE TO ORBITAL RESPONSE.
C     NON-INDEPENDENT PAIRS HAVE BEEN CALCULATED ELSEWHERE.
C
      IF (FON) THEN
        CALL DGEMM('N','N',L1,NOCC,L1-NDOCC,TWO,VEC(1,NDOCC+1),L1,
     *             U(NDOCC+1,1),L1,ZERO,WRK1,L1)
C
        IF (NDOCC.GE.1) CALL DSCAL(L1*NDOCC,OCC(1),WRK1,1)
        DO IO = NDOCC+1, NOCC
          CALL DSCAL(L1,OCC(IO),WRK1(1,IO),1)
          IF (NCP.EQ.1) CALL DAXPY(L1,DNDR(IO),VEC(1,IO),1,WRK1(1,IO),1)
        END DO
C       IF (NCP.GE.2) THEN
C         DO IO = NDOCC+1, NOCC
C           VAL = OCC(IO)*0.5D+00
C           VAL = -VAL*(ONE-VAL)*BETA
C           CALL DAXPY(L1,VAL*(FDERIJ(IO*(IO+1)/2)-FEDER)*TWO,
C    *        VEC(1,IO),1,WRK1(1,IO),1)
C         END DO
C       END IF
      ELSE
        CALL DGEMM('N','N',L1,NOCC,L1-NDOCC,TWO,
     *             VEC(1,NDOCC+1),L1,
     *             U(NDOCC+1,1),L1,ZERO,WRK1,L1)
        IF (NDOCC.GE.1) CALL DSCAL(L1*NDOCC,OCC,WRK1,1)
      END IF
      CALL DGEMM('N','T',L1,L1,NOCC,ONE,WRK1,L1,VEC,L1,ZERO,WRK,L1)
      CALL SYMSQT(L1,WRK,L1,WRK1,1)
C
      !! AO-RESOLVED MULLIKEN DERIVATIVE (WRK)
      CALL DFTB_MULLIKEN(L1,L2,WRK1,S,WRK)
C
      IF (SRSCC) THEN
        CALL DCOPY(NSHELL,QDER_INIT,1,QDER,1)
        ISH0 = 0
        DO NN = 1, NAT
          DO I = 1, MAXANG(ISPE(NN))
            ISH = ISH0 + I
            VAL = ZERO
            DO MM = 1, INDSAO(ISH+1)-INDSAO(ISH)
              MU = INDSAO(ISH)+MM
              VAL = VAL + WRK(MU,1)
            END DO
            QDER(ISH) = QDER(ISH) + VAL
          END DO
          ISH0 = ISH0 + MAXANG(ISPE(NN))
        END DO
      ELSE
        CALL DCOPY(NAT,QDER_INIT,1,QDER,1)
        DO NN = 1, NAT
          VAL = ZERO
          DO MM = 1, IND(NN+1)-IND(NN)
            MU = IND(NN)+MM
            VAL = VAL + WRK(MU,1)
          END DO
          QDER(NN) = QDER(NN) + VAL
        END DO
      END IF
C
      RETURN
C
      END SUBROUTINE DFTB_CALC_QDER
C
C-----------------------------------------------------------------------
C*MODULE DFTBHS    *DECK DFTB_PREP_GAMMA_HESS
C>
C>    @brief Hessian of DFTB
C>
C>    @details Prepare the calculation of second-order derivative of
C>             gamma values, making TMP matrix
C>
C>    @author Yoshio Nishimoto
C>    - Sep, 2014- Subroutine written
C>
C>           --- INPUT ---
C>    @param HUBA,HUBB Hubbard values of atom A and B
C>    @param DFTB3 Whether DFTB3 or not
C>    @param DAMPXH Whether DAMPXH is true or not
C>           --- OUTPUT ---
C>    @param TMP Constants
C>
C
      SUBROUTINE DFTB_PREP_GAMMA_HESS(HUBA,HUBB,TMP,DFTB3,DAMPXH)
C
      IMPLICIT NONE
C
      DOUBLE PRECISION, INTENT(IN) :: HUBA,HUBB
      LOGICAL, INTENT(IN) :: DFTB3,DAMPXH
      DOUBLE PRECISION, INTENT(OUT) :: TMP(11)
C
      DOUBLE PRECISION :: TAUA,TAUA2,TAUA3,TAUA4,TAUA6,TAUB,TAUB2,
     *  TAUB3,TAUB4,TAUB6,TAUAB,TAUABI,TAUABI2,TAUABI3,TAUABI4
      DOUBLE PRECISION,PARAMETER :: HUBTOL=1.0D-04,
     *  HALF=0.5D+00,ONE=1.0D+00,TWO=2.0D+00,THREE=3.0D+00,
     *  TWELVE=1.2D+01,ONE_15=ONE/1.5D+01,ONE_48=ONE/4.8D+01
C
      TAUA  = 3.2D+00*HUBA
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
        TMP(1) = 6.875D-01*TAUA  !! 11/16 a
        TMP(2) = 1.875D-01*TAUA2 !! 3/16  a^2
        TMP(3) = ONE_48*TAUA3    !! 1/48  a^3
        TMP(4) = TWO*TMP(3)  !! 1/24  a^3
        TMP(5) = 3.125D-01*TAUA3 !! 5/16  a^3
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
      END SUBROUTINE DFTB_PREP_GAMMA_HESS
C
C-----------------------------------------------------------------------
C*MODULE DFTBHS    *DECK DFTB_GAMHESSF
C>
C>    @brief Hessian of DFTB
C>
C>    @details Actual function which calculates the second-order
C>             derivative of gamma values
C>
C>    @author Yoshio Nishimoto
C>    - Sep, 2014- Subroutine written
C>    @date Feb, 2018 - Yoshio Nishimoto
C>    - Third-order derivative
C>
C>           --- INPUT ---
C>    @param NDER      Order of derivatives
C>    @param HUBA,HUBB Hubbard values of atom A and B
C>    @param HUBDER    Hubbard derivatives
C>    @param DIST      Inter-atomic distance between atom A and B
C>    @param DAMPXHE   Exponent of X-H damping function
C>    @param TMP Constants produced in DFTB_PREP_GAMMA_HESS subroutine
C>    @param DFTB3 Whether DFTB3 or not
C>    @param HDAMP Whether the pair contains H or not
C>           --- OUTPUT ---
C>    @param GRAD2 1st-order derivative of DFTB2 gamma
C>    @param HESS2 2nd-order derivative of DFTB2 gamma
C>    @param GRAD3 1st-order derivative of DFTB3 gamma
C>    @param HESS3 2nd-order derivative of DFTB3 gamma
C>    @param DER32 3rd-order derivative of DFTB2 gamma
C>    @param DER33 3rd-order derivative of DFTB3 gamma
C>
C
      SUBROUTINE DFTB_GAMHESSF(NDER,GRAD2,HESS2,DER32,GRAD3,HESS3,DER33,
     *  HUBA,HUBB,HUBDER,DIST,DAMPXHE,TMP,DFTB3,HDAMP)
C
      IMPLICIT NONE
C
      INTEGER NDER
      DOUBLE PRECISION, INTENT(IN) :: HUBA,HUBB,HUBDER,DIST,DAMPXHE,
     *  TMP(11)
      DOUBLE PRECISION, INTENT(OUT) :: GRAD2,HESS2,DER32,GRAD3,HESS3,
     *  DER33
      LOGICAL, INTENT(IN) :: DFTB3,HDAMP
C
      DOUBLE PRECISION :: TAUA,TAUA2,TAUB,TAUB2,DIST2,DIST3,DISTI,
     *  DISTI2,DISTI3,EXPA,GAMTMP,GAMTMP2,GAMTMP3,GA,DGA_DR,D2GA_DR,S,
     *  DS_DR,D2S_DR,D2S_DR2,HUBMEAN,HUBMEANEXP,EXPH,DH_DR,D2H_DR2,
     *  DGA_DA,D2GA_DADR,D2GA_DR2,D3GA_DADR2,D2S_DUDR,D3S_DUDR2,
     *  HUBMEANEXP1,DS_DU,DH_DU,D2H_DUDR,D3H_DUDR2,EXPB,FAB,FBA,
     *  DFAB_DR,DFBA_DR,D2FAB_DR,D2FBA_DR,DFAB_DA,DFBA_DA,D2FAB_DR2,
     *  D2FAB_DADR,D2FBA_DADR,D3FAB_DADR2,D3FBA_DADR2
      DOUBLE PRECISION DISTI4,D3GA_DR,D3S_DR,D3FAB_DR,D3FBA_DR
C     DOUBLE PRECISION DISTI5,D4GA_DR,D4S_DR
      DOUBLE PRECISION, PARAMETER :: DISTTOL=1.0D-08,HUBTOL=1.0D-04,
     *  ZERO=0.0D+00,ONE=1.0D+00,TWO=2.0D+00,HALF=0.5D+00
C
      IF (DIST.LT.DISTTOL) THEN
        GRAD2 = ZERO
        HESS2 = ZERO
        GRAD3 = ZERO
        HESS3 = ZERO
        IF (NDER.GE.3) THEN
          DER32 = ZERO
          DER33 = ZERO
        END IF
      ELSE
        TAUA = 3.2D+00*HUBA
        TAUA2 = TAUA**2
        TAUB = 3.2D+00*HUBB
        TAUB2 = TAUB**2
        DISTI = ONE/DIST
        DIST2 = DIST**2
        DISTI2 = DISTI**2
        DISTI3 = DISTI**3
        IF (NDER.GE.3) DISTI4 = DISTI2*DISTI2
C       IF (NDER.EQ.4) DISTI5 = DISTI4*DISTI
        EXPA = EXP(-TAUA*DIST)
        GAMTMP2 = ZERO
        GAMTMP3 = ZERO
        IF (ABS(HUBA-HUBB).LE.HUBTOL) THEN
          GA = DISTI+TMP(1)+TMP(2)*DIST+TMP(3)*DIST2
          DGA_DR = -DISTI2+TMP(2)+TMP(4)*DIST
          D2GA_DR = TWO*DISTI3+TMP(4)
          S = EXPA*GA
          DS_DR = EXPA*(DGA_DR-TAUA*GA)
          D2S_DR = EXPA*(D2GA_DR-TWO*TAUA*DGA_DR+TAUA2*GA)
          GAMTMP = -DS_DR
          D2S_DR2 = D2S_DR
          IF (NDER.GE.3) THEN
            D3GA_DR = -6.0D+00*DISTI4
            D3S_DR = -TAUA*D2S_DR
     *             + EXPA*(D3GA_DR-TWO*TAUA*D2GA_DR+TAUA2*DGA_DR)
C           IF (NDER.EQ.4) THEN
C             D4GA_DR = 2.4D+01*DISTI5
C             D4S_DR = -TAUA*D3S_DR
C    *               -TAUA*EXPA*(D3GA_DR-TWO*TAUA*D2GA_DR+TAUA2*DGA_DR)
C    *               + EXPA*(D4GA_DR-TWO*TAUA*D3GA_DR+TAUA2*D2GA_DR)
C           END IF
          END IF
          IF (HDAMP) THEN
            HUBMEAN = HUBA
            HUBMEANEXP = HUBMEAN**DAMPXHE
            EXPH = EXP(-HUBMEANEXP*DIST2)
            DH_DR = -TWO*DIST*HUBMEANEXP*EXPH
            D2H_DR2 = -TWO*HUBMEANEXP
     `        *(ONE-TWO*DIST2*HUBMEANEXP)*EXPH
            GAMTMP = GAMTMP*EXPH !! -h*dSg/dR
            GAMTMP2 = -DH_DR*S
            !! D2S_DR = d2S/dR2 + 2*(dS/dR)*(dh/dR)
            D2S_DR = D2S_DR*EXPH+TWO*DS_DR*DH_DR+S*D2H_DR2
          END IF
          GRAD2 = -DISTI2+GAMTMP+GAMTMP2
          HESS2 = TWO*DISTI3-D2S_DR
          IF (NDER.GE.3) DER32 = -6.0D+00*DISTI4-D3S_DR
C         IF (NDER.EQ.4) DER42 =  2.4D+01*DISTI5-D4S_DR !! ???
C         write (*,*) disti2,gamtmp,gamtmp2
C         write (*,*) disti3,d2s_dr
          IF (DFTB3) THEN
            DIST3 = DIST2*DIST
            DGA_DA =
     *        6.875D-01+3.75D-01*TAUA*DIST+6.25D-02*TAUA2*DIST2
            D2GA_DADR = 3.75D-01*TAUA+1.25D-01*TAUA2*DIST
            D2GA_DR2 = TWO*DISTI3+TMP(4)
            D3GA_DADR2 = 1.25D-01*TAUA2
            GAMTMP = (TAUA*DIST-ONE)*GA-TAUA*DGA_DA+D2GA_DADR
     *        - DIST*DGA_DR
            D2S_DUDR = GAMTMP*EXPA*3.2D+00
            GAMTMP = -D2S_DUDR
            GAMTMP3 = -TAUA*(TAUA*DIST-TWO)*GA
     *        +TWO*(TAUA*DIST-ONE)*DGA_DR+TAUA2*DGA_DA
     *        -TWO*TAUA*D2GA_DADR+D3GA_DADR2-DIST*D2GA_DR2
            D3S_DUDR2 = GAMTMP3*EXPA*3.2D+00
            GAMTMP3 = -D3S_DUDR2
            IF (HDAMP) THEN !! damp correction with 3rd
              HUBMEANEXP1 = HUBMEAN**(DAMPXHE-ONE)
              DS_DU = EXPA*(-ONE+TMP(9)*DIST
     *          + TMP(10)*DIST2+TMP(11)*DIST3)
              DH_DU = -HALF*DIST2*DAMPXHE*HUBMEANEXP1*EXPH
              D2H_DUDR = DAMPXHE*DIST*HUBMEANEXP1
     *                   *(DIST2*HUBMEANEXP-ONE)*EXPH
              D3H_DUDR2 = DAMPXHE*HUBMEANEXP1*EXPH
     *          *(DIST2*HUBMEANEXP
     *             *(5.0D+00-TWO*DIST2*HUBMEANEXP)-ONE)
C             WRITE (*,*) GAMTMP*EXPH,-DS_DU*DH_DR,-DS_DR*DH_DU
C             WRITE (*,*) -S*D2H_DUDR
              GAMTMP = -D2S_DUDR*EXPH-DS_DU*DH_DR-DS_DR*DH_DU
     *          -S*D2H_DUDR
C             D2S_DR2 = D2S_DR
              GAMTMP3 = -D3S_DUDR2*EXPH-TWO*D2S_DUDR*DH_DR
     *          -TWO*DS_DR*D2H_DUDR-DS_DU*D2H_DR2-D2S_DR2*DH_DU
     *          -S*D3H_DUDR2
            END IF
C           GAMTMP = 3.2D+00*GAMTMP !! /3.0D+00
C           GAMTMP3 = 3.2D+00*GAMTMP3 !! /3.0D+00
            GRAD3 = GAMTMP*HUBDER
            HESS3 = GAMTMP3*HUBDER
C           WRITE (*,*) "i,j = ",i ,j
C           write (*,*) gamtmp*HUBDER(NSP1),GAMTMP3*HUBDER(NSP1)
          END IF
        ELSE
          EXPB = EXP(-TAUB*DIST)
          FAB = TMP(2)-TMP(1)*DISTI
          FBA = TMP(4)-TMP(3)*DISTI
          DFAB_DR = TMP(1)*DISTI2
          DFBA_DR = TMP(3)*DISTI2
          D2FAB_DR = -DFAB_DR*TWO*DISTI
          D2FBA_DR = -DFBA_DR*TWO*DISTI
          S = EXPA*FAB+EXPB*FBA
          DS_DR = EXPA*(DFAB_DR-TAUA*FAB)+EXPB*(DFBA_DR-TAUB*FBA)
          GAMTMP = -DS_DR
          D2S_DR = EXPA*(D2FAB_DR-TWO*TAUA*DFAB_DR+TAUA2*FAB)
     *           + EXPB*(D2FBA_DR-TWO*TAUB*DFBA_DR+TAUB2*FBA)
          D2S_DR2 = D2S_DR
          IF (NDER.GE.3) THEN
            D3FAB_DR =  -D2FAB_DR*TWO*DISTI + DFAB_DR*TWO*DISTI2
            D3FBA_DR =  -D2FBA_DR*TWO*DISTI + DFBA_DR*TWO*DISTI2
            D3S_DR =
     *             - TAUA*EXPA*(D2FAB_DR-TWO*TAUA*DFAB_DR+TAUA2*FAB)
     *             + EXPA*(D3FAB_DR-TWO*TAUA*D2FAB_DR+TAUA2*DFAB_DR)
     *             - TAUB*EXPB*(D2FBA_DR-TWO*TAUB*DFBA_DR+TAUB2*FBA)
     *             + EXPB*(D3FBA_DR-TWO*TAUB*D2FBA_DR+TAUB2*DFBA_DR)
          END IF
          IF (HDAMP) THEN
            HUBMEAN = HALF*(HUBA+HUBB)
            HUBMEANEXP = HUBMEAN**DAMPXHE
            EXPH = EXP(-HUBMEANEXP*DIST2)
            DH_DR = -TWO*DIST*HUBMEANEXP*EXPH
            D2H_DR2 = -TWO*HUBMEANEXP
     `        *(ONE-TWO*DIST2*HUBMEANEXP)*EXPH
            GAMTMP = GAMTMP*EXPH !! -h*dSg/dR
            GAMTMP2 = -S*DH_DR
            !! D2S_DR = d2S/dR2 + 2*(dS/dR)*(dh/dR)
            D2S_DR = D2S_DR*EXPH+TWO*DS_DR*DH_DR+S*D2H_DR2
          END IF
          GRAD2 = -DISTI2+GAMTMP+GAMTMP2
C         write (*,*) disti2,gamtmp,gamtmp2
          HESS2 = TWO*DISTI3-D2S_DR
C         write (*,*) disti3,d2s_dr
          IF (NDER.GE.3) DER32 = -6.0D+00*DISTI4-D3S_DR
C         IF (NDER.GE.4) DER42 = 2.4D+01*DISTI5-D4S_DR
          IF (DFTB3) THEN
            DFAB_DA = TMP(5) + TMP(6)*DISTI
            DFBA_DA = TMP(7) + TMP(8)*DISTI
            D2FAB_DR2 = -DFAB_DR*TWO*DISTI
            D2FAB_DADR = -TMP(6)*DISTI2
            D2FBA_DADR = -D2FAB_DADR
            D3FAB_DADR2 = -D2FAB_DADR*TWO*DISTI
            D3FBA_DADR2 = -D2FBA_DADR*TWO*DISTI
            GAMTMP = EXPA*((TAUA*DIST-ONE)*FAB-TAUA*DFAB_DA
     *        + D2FAB_DADR - DIST*DFAB_DR)
     *        + EXPB*(D2FBA_DADR - TAUB*DFBA_DA)
            D2S_DUDR = GAMTMP*3.2D+00
            GAMTMP = -D2S_DUDR
            GAMTMP3 = EXPA*(-TAUA*(TAUA*DIST-TWO)*FAB
     *        +TWO*(TAUA*DIST-ONE)*DFAB_DR+TAUA2*DFAB_DA
     *        -TWO*TAUA*D2FAB_DADR+D3FAB_DADR2
     *        -DIST*D2FAB_DR2)
     *        +EXPB*(TAUB2*DFBA_DA-TWO*TAUB*D2FBA_DADR
     *        +D3FBA_DADR2)
            D3S_DUDR2 = GAMTMP3*3.2D+00
            GAMTMP3 = -D3S_DUDR2
C           WRITE (*,*) D2S_DUDR
C           WRITE (*,*) "D3S_DUDR2 = ", D3S_DUDR2
            IF (HDAMP) THEN  !! damp correction with 3rd
              HUBMEANEXP1 = HUBMEAN**(DAMPXHE-ONE)
              DS_DU = EXPA*(DFAB_DA-DIST*FAB)+EXPB*DFBA_DA
              DS_DU = DS_DU*3.2D+00
              DH_DR = -TWO*DIST*HUBMEANEXP*EXPH
              DH_DU = -HALF*DIST2*DAMPXHE*HUBMEANEXP1*EXPH
              D2H_DUDR = DAMPXHE*DIST*HUBMEANEXP1
     *                   *(DIST2*HUBMEANEXP-ONE)*EXPH
C             WRITE (*,*) D2H_DUDR
              D3H_DUDR2 = DAMPXHE*HUBMEANEXP1*EXPH
     *          *(DIST2*HUBMEANEXP
     *             *(5.0D+00-TWO*DIST2*HUBMEANEXP)-ONE)
C             WRITE (*,*) "D3H_DUDR2 = ", D3H_DUDR2
C             WRITE (*,*) GAMTMP*EXPH,-DS_DU*DH_DR*3.2D+00,
C    *        -DS_DR*DH_DU
C             WRITE (*,*) -S*D2H_DUDR
              GAMTMP = -D2S_DUDR*EXPH-DS_DU*DH_DR
     *          -DS_DR*DH_DU-S*D2H_DUDR
C             D2S_DR2 = D2S_DR
              GAMTMP3 = -D3S_DUDR2*EXPH-TWO*D2S_DUDR*DH_DR
     *          -TWO*DS_DR*D2H_DUDR-DS_DU*D2H_DR2-D2S_DR2*DH_DU
     *          -S*D3H_DUDR2
            END IF
C           GAMTMP = 3.2D+00*GAMTMP !! /3.0D+00
C           GAMTMP3 = 3.2D+00*GAMTMP3 !! /3.0D+00
            GRAD3 = GAMTMP*HUBDER
            HESS3 = GAMTMP3*HUBDER
C           WRITE (*,*) "i,j = ",i ,j
C           write (*,*) gamtmp*HUBDER(NSP1),GAMTMP3*HUBDER(NSP1)
          END IF
        END IF
C       WRITE (*,'(I3,X,I3,X,F12.6,X,F12.6)')
C    *    1,1,GRAD2,HESS2
      END IF
C
      RETURN
C
      END SUBROUTINE DFTB_GAMHESSF
C
C-----------------------------------------------------------------------
C*MODULE DFTBHS    *DECK DFTB_CALC_GAM_DER2
C>
C>    @brief Hessian of DFTB
C>
C>    @details Calculate second-order geometrical derivative of gamma
C>
C>    @author Yoshio Nishimoto
C>    - Sep, 2014- Subroutine written
C>    @date Feb, 2016 - Yoshio Nishimoto
C>    - Bug fix
C>    @date Feb, 2018 - Yoshio Nishimoto
C>    - Third-order derivative
C>
C>           --- INPUT ---
C>    @param NDER    Order of the derivative
C>    @param DISTMAT Matrix of inter-atomic distances
C>    @param HUBBL Shell-resolved Hubbar values
C>    @param HUBDER Hubbard derivatives
C>    @param DAMPXHE Exponent of X-H damping function
C>    @param ISPE Index of species
C>    @parma MAXANG Maximum angular momentum of each species
C>    @param NAT Number of atoms
C>    @param Number of species
C>    @param Number of shells
C>    @param NDIMGAM Dimension of gamma
C>    @param NN 
C>    @param SPE Name of species
C>    @param DFTB3 Whether DFTB3 or not
C>    @param SRSCC Whether shell-resolved SCC or not
C>    @param DAMPXH Whether H-X damping is used or not
C>           --- OUTPUT ---
C>    @param GRAD2 1st-order derivative of DFTB2 gamma
C>    @param HESS2 2nd-order derivative of DFTB2 gamma
C>    @param GRAD3 1st-order derivative of DFTB3 gamma
C>    @param HESS3 2nd-order derivative of DFTB3 gamma
C>    @param DER32 3rd-order derivative of DFTB2 gamma
C>    @param DER33 3rd-order derivative of DFTB3 gamma
C>
C
      SUBROUTINE DFTB_CALC_GAM_DER2(NDER,GRAD2,HESS2,DER32,GRAD3,HESS3,
     *  DER33,DISTMAT,HUBBL,HUBDER,DAMPXHE,ISPE,MAXANG,NAT,NSPE,NSHELL,
     *  NDIMGAM,NN,SPE,DFTB3,SRSCC,DAMPXH)
C
      IMPLICIT NONE
C
      DOUBLE PRECISION, INTENT(INOUT) :: GRAD2(*),HESS2(*),DER32(*),
     *  GRAD3(NDIMGAM,NDIMGAM),HESS3(NDIMGAM,NDIMGAM),
     *  DER33(NDIMGAM,NDIMGAM)
      DOUBLE PRECISION, INTENT(IN) :: DISTMAT(*),HUBBL(3,NSPE),DAMPXHE,
     *  HUBDER(*),SPE(*)
      INTEGER, INTENT(IN) :: NDER,ISPE(NAT),MAXANG(NSPE),NAT,NSPE,
     *  NSHELL,NDIMGAM,NN
      LOGICAL, INTENT(IN) :: DFTB3,SRSCC,DAMPXH
C
      INTEGER :: I,II,ISH,J,JJ,JSH,NSP1,NSP2,NSEQ,ISHELL(NSHELL),
     *  IND(NSPE),ISH0,JSH0
      DOUBLE PRECISION :: HUBA,HUBB,HUBDERTMP,DIST,CHARA1,CHARA2,
     *  TMP(11,NN,NN)
      LOGICAL :: HDAMP
      DOUBLE PRECISION H/8HH       /
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
                CALL DFTB_PREP_GAMMA_HESS(HUBA,HUBB,TMP(1,ISH,JSH),
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
            NSP2 = ISPE(J)
            CHARA2 = SPE(ISPE(J))
            CALL DFTB_CNVSQ(I,J,NSEQ)
            DIST = DISTMAT(NSEQ)
            IF (DAMPXH.AND.(CHARA1.EQ.H.OR.CHARA2.EQ.H)) THEN
              HDAMP = .TRUE.
              CALL DFTB_DMPCHK(NSP1,NSP2,HDAMP)
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
                CALL DFTB_GAMHESSF(NDER,GRAD2(NSEQ),HESS2(NSEQ),
     *            DER32(NSEQ),GRAD3(ISH,JSH),HESS3(ISH,JSH),
     *            DER33(ISH,JSH),HUBA,HUBB,HUBDERTMP,DIST,DAMPXHE,
     *            TMP(1,NSP1,NSP2),DFTB3,HDAMP)
C      write (*,'(6i3,f20.10)') i,j,ish,jsh,nsp1,nsp2,gamma(nseq)
C       WRITE (*,'(I3,X,I3,X,I3,F20.10,X,F20.10)')
C    *    ISH,JSH,NSEQ,GRAD2(NSEQ),HESS2(NSEQ)
              END DO
            END DO
            JSH0 = JSH0 + MAXANG(ISPE(J))
          END DO
          ISH0 = ISH0 + MAXANG(ISPE(I))
        END DO
      ELSE
        DO I = 1, NSPE
          HUBA  = HUBBL(1,I)
          DO J = 1, NSPE
            HUBB = HUBBL(1,J)
C           WRITE (*,*) I,J
            CALL DFTB_PREP_GAMMA_HESS(HUBA,HUBB,TMP(1,I,J),DFTB3,DAMPXH)
C             write (*,'(x,f20.10)') tmp(1,i,j)
C             write (*,'(x,f20.10)') tmp(2,i,j)
C             write (*,'(x,f20.10)') tmp(3,i,j)
C             write (*,'(x,f20.10)') tmp(4,i,j)
C             write (*,'(x,f20.10)') tmp(5,i,j)
C             write (*,'(x,f20.10)') tmp(6,i,j)
C             write (*,'(x,f20.10)') tmp(7,i,j)
C             write (*,'(x,f20.10)') tmp(8,i,j)
C             write (*,'(x,f20.10)') tmp(9,i,j)
C             write (*,'(x,f20.10)') tmp(10,i,j)
C             write (*,'(x,f20.10)') tmp(11,i,j)
C             write (*,*)
          END DO
        END DO
C
        DO I = 1, NAT
          NSP1 = ISPE(I)
          HUBA = HUBBL(1,NSP1)
          CHARA1 = SPE(ISPE(I))
          HUBDERTMP = HUBDER(NSP1)
          DO J = 1, NAT !! I
C           WRITE (*,*) I,J
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
            CALL DFTB_GAMHESSF(NDER,GRAD2(NSEQ),HESS2(NSEQ),DER32(NSEQ),
     *        GRAD3(I,J),HESS3(I,J),DER33(I,J),HUBA,HUBB,HUBDERTMP,
     *        DIST,DAMPXHE,TMP(1,NSP1,NSP2),DFTB3,HDAMP)
C       WRITE (*,'(I3,X,I3,X,F12.6,X,F12.6)')
C    *    I,J,GRAD2(NSEQ),HESS2(NSEQ)
          END DO
        END DO
      END IF
C
      END SUBROUTINE DFTB_CALC_GAM_DER2
C
C-----------------------------------------------------------------------
C*MODULE DFTBHS    *DECK DFTB_GAM_DER2
C>
C>    @brief Hessian of DFTB
C>
C>    @details Calculate Hessian contribution of the second-order
C>             derivative of gamma values
C>
C>    @author Yoshio Nishimoto
C>    - Sep, 2014- Subroutine written
C>    @date Feb, 2016 - Yoshio Nishimoto
C>    - Added DFTB3 with SRSCC
C>
C>           --- INPUT ---
C>    @param NAT Number of atoms
C>    @param NSPE Number of species
C>    @param HUBBL Shell-resolved Hubbar values
C>    @param C Coordinates
C>    @param DISTMAT Matrix of inter-atomic distances
C>    @param CHAMUL Atom-reslve mulliken charges
C>    @parma HESSIAN Hessian
C>    @param DFTB3 Whether DFTB3 or not
C>    @param HUBDER Hubbard derivatives
C>    @param DAMPXHE Exponent of X-H damping function
C>    @param ISPE Index of species
C>    @parma MAXANG Maximum angular momentum of each species
C>    @param NSHELL Number of shells
C>    @param NDIMGAM Dimension of gamma
C>    @param NN
C>    @param SPE Name of species
C>    @param SRSCC Whether shell-resolved SCC or not
C>    @param DAMPXH Whether H-X damping is used or not
C>           --- OUTPUT ---
C>    @param GRADMAT 1st-order derivative of DFTB2 gamma
C>    @param HESSMAT 2nd-order derivative of DFTB2 gamma
C>    @param GRADMAT3 1st-order derivative of DFTB3 gamma
C>    @param HESSMAT3 2nd-order derivative of DFTB3 gamma
C>
C
      SUBROUTINE DFTB_GAM_DER2(NAT,NSPE,HUBBL,C,DISTMAT,GRADMAT,HESSMAT,
     *  GRADMAT3,HESSMAT3,CHAMUL,HESSIAN,DFTB3,HUBDER,DAMPXHE,ISPE,
     *  MAXANG,NSHELL,NDIMGAM,NN,SPE,SRSCC,DAMPXH,CHAMULAT)
C
      IMPLICIT NONE
C
      INTEGER, INTENT(IN) :: NAT,NSPE,ISPE(NAT),MAXANG(NSPE),NSHELL,
     *  NDIMGAM,NN
      DOUBLE PRECISION, INTENT(IN) :: HUBBL(3,NSPE),C(3,*),DISTMAT(*),
     *  CHAMUL(*),HUBDER(NSPE),DAMPXHE,SPE(*),CHAMULAT(*)
      DOUBLE PRECISION, INTENT(OUT) :: GRADMAT(*),HESSMAT(*),
     *  GRADMAT3(NDIMGAM,NDIMGAM),HESSMAT3(NDIMGAM,NDIMGAM)
      DOUBLE PRECISION, INTENT(INOUT) :: HESSIAN(3*NAT,3*NAT)
      LOGICAL, INTENT(IN) :: DFTB3,SRSCC,DAMPXH
C
      INTEGER :: I,J,K,NI,NJ,NK,ISH,ISH0,JSH,JSH0,KSH,KSH0,KK,KI,KJ,IR,
     *  JR,NSEQ
      DOUBLE PRECISION :: VEC(3),CHAMUL1,CHAMUL2,HESS_TMP,DIST,DISTI,
     *  GRAD,HESS,GRAD3C,GRAD3A,HESS3C,HESS3A,HESS_TMPA,HESS_TMPC
      DOUBLE PRECISION, PARAMETER :: ONE=1.0D+00,ONE_THIRD=ONE/3.0D+00
C
C     CALL DFTB_CALC_GAM_DER2(HUBBL,DISTMAT,GRADMAT,HESSMAT,GRADMAT3,
C    *  HESSMAT3)
      CALL DFTB_CALC_GAM_DER2(2,GRADMAT,HESSMAT,HESSMAT,GRADMAT3,
     *  HESSMAT3,HESSMAT3,DISTMAT,HUBBL,HUBDER,DAMPXHE,ISPE,MAXANG,
     *  NAT,NSPE,NSHELL,NDIMGAM,NN,SPE,DFTB3,SRSCC,DAMPXH)
C     do i = 1, nshell
C       do j = 1, nshell
C        call DFTB_CNVSQ(i,j,nseq)
C       write (*,'(3i3,2f20.10)') i,j,nseq,gradmat(nseq),hessmat(nseq)
C     end do
C     end do
C
      IF (SRSCC) THEN
        ISH0 = 0
        DO NI = 1, NAT !! A
          DO I = 1, MAXANG(ISPE(NI))
            ISH = ISH0 + I
            CHAMUL1 = CHAMUL(ISH)
            JSH0 = 0
            DO NJ = 1, NAT !! B
              IF (NI.EQ.NJ) THEN
                KSH0 = 0
                DO NK = 1, NAT !! C
                  IF (NI.EQ.NK) THEN
                    KSH0 = KSH0 + MAXANG(ISPE(NK))
                    CYCLE
                  END IF
                  CALL DFTB_CNVSQ(NI,NK,NSEQ)
                  DIST = DISTMAT(NSEQ)
                  DISTI = ONE/DIST
                  DO K = 1, MAXANG(ISPE(NK))
                    KSH = KSH0 + K
                    CHAMUL2 = CHAMUL(KSH)
                    CALL DFTB_CNVSQ(ISH,KSH,NSEQ)
                    GRAD = GRADMAT(NSEQ)
                    HESS = HESSMAT(NSEQ)
C                   write (*,'(3i3,2f20.10)') ish,ksh,nseq,grad,hess
                    IF (DFTB3) THEN
                      GRAD3C = GRADMAT3(KSH,ISH)
                      GRAD3A = GRADMAT3(ISH,KSH)
                      HESS3C = HESSMAT3(KSH,ISH)
                      HESS3A = HESSMAT3(ISH,KSH)
                    END IF
                    DO KK = 1, 3
                      VEC(KK) = C(KK,NI) - C(KK,NK)
                      VEC(KK) = VEC(KK)*DISTI
                    END DO
                    DO KI = 1, 3
                      IR = 3*(NI-1)+KI
                      DO KJ = 1, 3
                        JR = 3*(NJ-1)+KJ
                        IF (KI.EQ.KJ) THEN
                          HESS_TMP = VEC(KI)*VEC(KJ)*HESS
     *                      + (ONE-VEC(KI)*VEC(KJ))*DISTI*GRAD
                          IF (DFTB3) THEN
                            HESS_TMPA = VEC(KI)*VEC(KJ)*HESS3A
     *                        + (ONE-VEC(KI)*VEC(KJ))*DISTI*GRAD3A
                            HESS_TMPC = VEC(KI)*VEC(KJ)*HESS3C
     *                        + (ONE-VEC(KI)*VEC(KJ))*DISTI*GRAD3C
                          END IF
                        ELSE
                          HESS_TMP = VEC(KI)*VEC(KJ)*HESS
     *                      - VEC(KI)*VEC(KJ)*DISTI*GRAD
                          IF (DFTB3) THEN
                            HESS_TMPA = VEC(KI)*VEC(KJ)*HESS3A
     *                        - VEC(KI)*VEC(KJ)*DISTI*GRAD3A
                            HESS_TMPC = VEC(KI)*VEC(KJ)*HESS3C
     *                        - VEC(KI)*VEC(KJ)*DISTI*GRAD3C
                          END IF
                        END IF
                        HESSIAN(IR,JR) = HESSIAN(IR,JR)
     *                    + HESS_TMP*CHAMUL1*CHAMUL2
C                       write (*,'(i3,i3,3f20.10)')
C    *                  ir,jr,hess_tmp,chamul1,chamul2
C                       write (*,'(i3,i3,3f20.10)')
C    *                  ir,jr,dist,grad,hess
                        IF (DFTB3) THEN
                          HESSIAN(IR,JR)= HESSIAN(IR,JR)+CHAMUL1*CHAMUL2
     *                      *(CHAMULAT(NI)*HESS_TMPA
     *                      + CHAMULAT(NK)*HESS_TMPC)*ONE_THIRD
C                        HESSIAN(IR,JR) = HESSIAN(IR,JR)+CHAMUL1*CHAMUL2
C    *                    *(CHAMUL1*HESS_TMPA+CHAMUL2*HESS_TMPC)/3.0D+00
C                         WRITE (*,*) HESS_TMPA,HESS_TMPC
                        END IF
                      END DO
                    END DO
                  END DO
                  KSH0 = KSH0 + MAXANG(ISPE(NK))
                END DO
              ELSE
                CALL DFTB_CNVSQ(NI,NJ,NSEQ)
                DIST = DISTMAT(NSEQ)
                DISTI = ONE/DIST
                DO J = 1, MAXANG(ISPE(NJ))
                  JSH = JSH0 + J
                  CHAMUL2 = CHAMUL(JSH)
                  CALL DFTB_CNVSQ(ISH,JSH,NSEQ)
                  GRAD = GRADMAT(NSEQ)
                  HESS = HESSMAT(NSEQ)
                  IF (DFTB3) THEN
                    GRAD3C = GRADMAT3(JSH,ISH)
                    GRAD3A = GRADMAT3(ISH,JSH)
                    HESS3C = HESSMAT3(JSH,ISH)
                    HESS3A = HESSMAT3(ISH,JSH)
                  END IF
                  DO KK = 1, 3
                    VEC(KK) = C(KK,NI) - C(KK,NJ)
                    VEC(KK) = VEC(KK)*DISTI
                  END DO
                  DO KI = 1, 3
                    IR = 3*(NI-1)+KI
                    DO KJ = 1, 3
                      JR = 3*(NJ-1)+KJ
                      IF (KI.EQ.KJ) THEN
                        HESS_TMP = -VEC(KI)*VEC(KJ)*HESS
     *                           - (ONE-VEC(KI)*VEC(KJ))*DISTI*GRAD
                        IF (DFTB3) THEN
                          HESS_TMPA = -VEC(KI)*VEC(KJ)*HESS3A
     *                              - (ONE-VEC(KI)*VEC(KJ))*DISTI*GRAD3A
                          HESS_TMPC = -VEC(KI)*VEC(KJ)*HESS3C
     *                              - (ONE-VEC(KI)*VEC(KJ))*DISTI*GRAD3C
                        END IF
                      ELSE
                        HESS_TMP = -VEC(KI)*VEC(KJ)*HESS
     *                           + VEC(KI)*VEC(KJ)*DISTI*GRAD
                        IF (DFTB3) THEN
                          HESS_TMPA = -VEC(KI)*VEC(KJ)*HESS3A
     *                              + VEC(KI)*VEC(KJ)*DISTI*GRAD3A
                          HESS_TMPC = -VEC(KI)*VEC(KJ)*HESS3C
     *                              + VEC(KI)*VEC(KJ)*DISTI*GRAD3C
                        END IF
                      END IF
                      HESSIAN(IR,JR) = HESSIAN(IR,JR)
     *                  + HESS_TMP*CHAMUL1*CHAMUL2
C                       write (*,'(i3,i3,3f20.10)')
C    *                  ir,jr,hess_tmp,chamul1,chamul2
                      IF (DFTB3) THEN
                        HESSIAN(IR,JR)= HESSIAN(IR,JR)+CHAMUL1*CHAMUL2
     *                    *(CHAMULAT(NI)*HESS_TMPA
     *                    + CHAMULAT(NJ)*HESS_TMPC)*ONE_THIRD
C                       HESSIAN(IR,JR) = HESSIAN(IR,JR)+CHAMUL1*CHAMUL2
C    *                  *(CHAMUL1*HESS_TMPA+CHAMUL2*HESS_TMPC)*ONE_THIRD
C                         WRITE (*,*) HESS_TMPA,HESS_TMPC
                      END IF
                    END DO
                  END DO
                END DO
              END IF
              JSH0 = JSH0 + MAXANG(ISPE(NJ))
            END DO
          END DO
          ISH0 = ISH0 + MAXANG(ISPE(NI))
        END DO
      ELSE
        DO NI = 1, NAT !! A
          CHAMUL1 = CHAMUL(NI)
          DO NJ = 1, NAT !! B
            IF (NI.EQ.NJ) THEN
              DO NK = 1, NAT !! C
                IF (NI.EQ.NK) CYCLE
                CHAMUL2 = CHAMUL(NK)
                CALL DFTB_CNVSQ(NI,NK,NSEQ)
                DIST = DISTMAT(NSEQ)
                GRAD = GRADMAT(NSEQ)
                HESS = HESSMAT(NSEQ)
                IF (DFTB3) THEN
                  GRAD3C = GRADMAT3(NK,NI)
                  GRAD3A = GRADMAT3(NI,NK)
                  HESS3C = HESSMAT3(NK,NI)
                  HESS3A = HESSMAT3(NI,NK)
                END IF
                DISTI = ONE/DIST
                DO KK = 1, 3
                  VEC(KK) = C(KK,NI) - C(KK,NK)
                  VEC(KK) = VEC(KK)*DISTI
                END DO
                DO KI = 1, 3
                  IR = 3*(NI-1)+KI
                  DO KJ = 1, 3
                    JR = 3*(NJ-1)+KJ
                    IF (KI.EQ.KJ) THEN
                      HESS_TMP = VEC(KI)*VEC(KJ)*HESS
     *                         + (ONE-VEC(KI)*VEC(KJ))*DISTI*GRAD
                      IF (DFTB3) THEN
                        HESS_TMPA = VEC(KI)*VEC(KJ)*HESS3A
     *                            + (ONE-VEC(KI)*VEC(KJ))*DISTI*GRAD3A
                        HESS_TMPC = VEC(KI)*VEC(KJ)*HESS3C
     *                            + (ONE-VEC(KI)*VEC(KJ))*DISTI*GRAD3C
                      END IF
                    ELSE
                      HESS_TMP = VEC(KI)*VEC(KJ)*HESS
     *                         - VEC(KI)*VEC(KJ)*DISTI*GRAD
                      IF (DFTB3) THEN
                        HESS_TMPA = VEC(KI)*VEC(KJ)*HESS3A
     *                            - VEC(KI)*VEC(KJ)*DISTI*GRAD3A
                        HESS_TMPC = VEC(KI)*VEC(KJ)*HESS3C
     *                            - VEC(KI)*VEC(KJ)*DISTI*GRAD3C
                      END IF
                    END IF
                    HESSIAN(IR,JR) = HESSIAN(IR,JR)
     *                + HESS_TMP*CHAMUL1*CHAMUL2
C                       write (*,'(i3,i3,3f20.10)')
C    *                  ir,jr,hess_tmp,chamul1,chamul2
                    IF (DFTB3) THEN
                      HESSIAN(IR,JR) = HESSIAN(IR,JR) + CHAMUL1*CHAMUL2
     *                  *(CHAMUL1*HESS_TMPA+CHAMUL2*HESS_TMPC)*ONE_THIRD
C                     WRITE (*,*) HESS_TMPA,HESS_TMPC
                    END IF
                  END DO
                END DO
              END DO
            ELSE
              CHAMUL2 = CHAMUL(NJ)
              CALL DFTB_CNVSQ(NI,NJ,NSEQ)
              DIST = DISTMAT(NSEQ)
              GRAD = GRADMAT(NSEQ)
              HESS = HESSMAT(NSEQ)
              DISTI = ONE/DIST
              IF (DFTB3) THEN
                GRAD3C = GRADMAT3(NJ,NI)
                GRAD3A = GRADMAT3(NI,NJ)
                HESS3C = HESSMAT3(NJ,NI)
                HESS3A = HESSMAT3(NI,NJ)
              END IF
              DO KK = 1, 3
                VEC(KK) = C(KK,NI) - C(KK,NJ)
                VEC(KK) = VEC(KK)*DISTI
              END DO
              DO KI = 1, 3
                IR = 3*(NI-1)+KI
                DO KJ = 1, 3
                  JR = 3*(NJ-1)+KJ
                  IF (KI.EQ.KJ) THEN
                    HESS_TMP = -VEC(KI)*VEC(KJ)*HESS
     *                       - (ONE-VEC(KI)*VEC(KJ))*DISTI*GRAD
                    IF (DFTB3) THEN
                      HESS_TMPA = -VEC(KI)*VEC(KJ)*HESS3A
     *                          - (ONE-VEC(KI)*VEC(KJ))*DISTI*GRAD3A
                      HESS_TMPC = -VEC(KI)*VEC(KJ)*HESS3C
     *                          - (ONE-VEC(KI)*VEC(KJ))*DISTI*GRAD3C
                    END IF
                  ELSE
                    HESS_TMP = -VEC(KI)*VEC(KJ)*HESS
     *                       + VEC(KI)*VEC(KJ)*DISTI*GRAD
                    IF (DFTB3) THEN
                      HESS_TMPA = -VEC(KI)*VEC(KJ)*HESS3A
     *                          + VEC(KI)*VEC(KJ)*DISTI*GRAD3A
                      HESS_TMPC = -VEC(KI)*VEC(KJ)*HESS3C
     *                          + VEC(KI)*VEC(KJ)*DISTI*GRAD3C
                    END IF
                  END IF
                  HESSIAN(IR,JR) = HESSIAN(IR,JR)
     *              + HESS_TMP*CHAMUL1*CHAMUL2
C                       write (*,'(i3,i3,3f20.10)')
C    *                  ir,jr,hess_tmp,chamul1,chamul2
                  IF (DFTB3) THEN
                    HESSIAN(IR,JR) = HESSIAN(IR,JR) + CHAMUL1*CHAMUL2
     *                *(CHAMUL1*HESS_TMPA+CHAMUL2*HESS_TMPC)*ONE_THIRD
C                     WRITE (*,*) HESS_TMPA,HESS_TMPC
                  END IF
                END DO
              END DO
            END IF
          END DO
        END DO
      END IF
C
      END SUBROUTINE DFTB_GAM_DER2
C
C-----------------------------------------------------------------------
C*MODULE DFTBHS    *DECK DFTB_DEDB
C>
C>    @brief   Hessian of DFTB
C>
C>    @details Calculate Hessian contribution of U-dependent terms and
C>             de/db dependent terms
C>
C>    @author  Yoshio Nishimoto
C>             Sep, 2014 - Subroutine written
C>    @date    Feb, 2016 - Yoshio Nishimoto and Hiroya Nakata
C>             - Performance improvement
C>             - Implementation of FMO-DFTB Hessian
C>    @date    Feb, 2018 - Yoshio Nishimoto
C>             - Almost rewritten to improve performance
C>
C>           --- INPUT ---
C>    @param SDERMO  Contains all S_{ij}^a
C>    @param DOMEGA  Derivatives of shift contribution in atom basis
C>    @parma SHIFT   Shift contribution (atom-resolved)
C>    @param SHIFTSA Shift contribution of alpha orbitals (shell-res.)
C>    @param SHIFTSB Shift contribution of beta  orbitals (shell-res.)
C>    @param C       Coordinates
C>    @param D       Density matrix
C>    @param NAT     Number of atoms
C>    @param L1      Number of basis functions
C>    @param L2      L1*(L1+1)/2
C>    @param NSPIN   Number of spins (1=R-DFTB,2=U-DFTB)
C>    @param NDIMSH  If (SRSCC); NDIMSH=NSHELL; else NDIMSH=NAT
C>    @param NSPE    Number of species
C>    @param ISPE    Index of species
C>    @parma MAXANG  Maximum angular momentum of each species
C>    @param IND     AO index of each atom
C>    @param INDSH   Index of shell for each atom
C>    @param INDSAO  AO index of each shell
C>    @param SCC     Whether SCC or not
C>    @param SRSCC   Whether shell-resolved SCC or not
C>    @param UHF     Whether U-DFTB or not
C>    @param MASWRK  MASWRK
C>    @param DDER    Derivative of density matrix
C>    @param WDER    Derivative of energy-weighted density matrix
C>    @param MAXBND  Maximum number of effective bonds
C>    @param INDBND  Index of bond and AO derivatives
C>    @param DERMAT  Derivative of AO integral
C>           --- OUTPUT ---
C>    @param HESSIAN Hessian
C>
C
      SUBROUTINE DFTB_DEDB(HESSIAN,HOPAO,DOMEGA,SHIFT,SHIFTSA,
     *  SHIFTSB,C,D,NAT,L1,L2,NSPIN,NDIMSH,
     *  NSPE,ISPE,MAXANG,IND,INDSH,INDSAO,
     *  SCC,SRSCC,UHF,MASWRK,
     *  DDER,WDER,MAXBND,INDBND,DERMAT)
C
      IMPLICIT NONE
C
      DOUBLE PRECISION, INTENT(INOUT) :: HESSIAN(3*NAT,*)
      DOUBLE PRECISION, INTENT(IN) :: HOPAO(L2,NSPIN,3,NAT),
     *  DOMEGA(NDIMSH,NSPIN),SHIFT(NAT),SHIFTSA(NDIMSH),
     *  SHIFTSB(NDIMSH),C(3,NAT),D(L2,*),DDER(L2,*),WDER(L2,*),
     *  DERMAT(2,*)
      INTEGER, INTENT(IN) :: NAT,L1,L2,NSPIN,NDIMSH,MAXBND,
     *  NSPE,ISPE(NAT),MAXANG(NSPE),IND(*),INDSH(*),INDSAO(*),
     *  INDBND(MAXBND+1,4,NAT)
      LOGICAL, INTENT(IN) :: SCC,SRSCC,UHF,MASWRK
C
      INTEGER :: NSEQ,IR,NI,KI,N,MU,NU,M,NJ,ISH,ISH0,JSH,JSH0,
     *  NSEQ2,IBND
      DOUBLE PRECISION :: VAL,TRACEP,SHIFTI,SHIFTJ,SHIFTIJ,SHIFTID,
     *  SHIFTJD,SHIFTIJD,DHDR,DSDR,SHIFTIA,SHIFTIB,SHIFTIAD,SHIFTIBD,
     *  SHIFTJAD,SHIFTJBD,SHIFTIJA,SHIFTIJB,SHIFTIJAD,SHIFTIJBD,SHIFTJA,
     *  SHIFTJB
      DOUBLE PRECISION, PARAMETER :: ZERO=0.0D+00,TWO=2.0D+00,
     *  HALF=0.5D+00
C
      common /fmoinf/ nfg,nlayer,natfmo,nbdfg,naotyp,nbody,nsegm
      INTEGER :: nfg,nlayer,natfmo,nbdfg,naotyp,nbody,nsegm
C
C     CALCULATE CONTRIBUTIONS OF U-DEPENDENT TERMS AND DE/DB DEPENDENT
C     TERMS
C
      IR = 0
      DO NI = 1, NAT
        IF (SCC) THEN
          SHIFTI  = SHIFT(NI)
          SHIFTID = DOMEGA(NI,1)
        ELSE
          SHIFTI  = ZERO
          SHIFTID = ZERO
        END IF
        DO KI = 1, 3
          IR = IR + 1
          VAL = ZERO
          DO IBND = 1, INDBND(1,1,NI)
            NJ = INDBND(IBND+1,1,NI)
            NSEQ2 = INDBND(IBND+1,KI+1,NI)
C           write (*,*) "nj = ", nj,nseq2
            IF (SRSCC.OR.UHF) THEN
              DO JSH0 = 1, MAXANG(ISPE(NJ))
                JSH = INDSH(NJ)+JSH0
                SHIFTJA  = SHIFTSA(JSH)
                SHIFTJAD = DOMEGA(JSH,1)
                IF (UHF) THEN
                  SHIFTJB  = SHIFTSB(JSH)
                  SHIFTJBD = DOMEGA(JSH,2)
                END IF
                DO N = 1, INDSAO(JSH+1)-INDSAO(JSH)
                  NU = INDSAO(JSH)+N
                  DO ISH0 = 1, MAXANG(ISPE(NI))
                    ISH = INDSH(NI)+ISH0
                    SHIFTIA   = SHIFTSA(ISH)
                    SHIFTIAD  = DOMEGA(ISH,1)
                    SHIFTIJA  = (SHIFTIA+SHIFTJA)*HALF
                    SHIFTIJAD = (SHIFTIAD+SHIFTJAD)*HALF
                    IF (UHF) THEN
                      SHIFTIB   = SHIFTSB(ISH)
                      SHIFTIBD  = DOMEGA(ISH,2)
                      SHIFTIJB  = (SHIFTIB+SHIFTJB)*HALF
                      SHIFTIJBD = (SHIFTIBD+SHIFTJBD)*HALF
                    END IF
                    DO M = 1, INDSAO(ISH+1)-INDSAO(ISH)
                      MU = INDSAO(ISH)+M
                      CALL DFTB_CNVSQ(MU,NU,NSEQ)
                      DHDR = DERMAT(1,NSEQ2)
                      DSDR = DERMAT(2,NSEQ2)
                      VAL = VAL
     *                    + DDER(NSEQ,1)*(DHDR+SHIFTIJA*DSDR)
     *                    - WDER(NSEQ,1)*DSDR
     *                    + D(NSEQ,1)*SHIFTIJAD*DSDR
                      IF (UHF) THEN
                      VAL = VAL
     *                    + DDER(NSEQ,2)*(DHDR+SHIFTIJB*DSDR)
     *                    - WDER(NSEQ,2)*DSDR
     *                    + D(NSEQ,2)*SHIFTIJBD*DSDR
                      END IF
                      NSEQ2 = NSEQ2 + 1
                    END DO
                  END DO
                END DO
              END DO
            ELSE
              IF (SCC) THEN
                SHIFTJ  = SHIFT(NJ)
                SHIFTJD = DOMEGA(NJ,1)
                SHIFTIJ = (SHIFTI +SHIFTJ )*HALF
                SHIFTIJD= (SHIFTID+SHIFTJD)*HALF
              ELSE
                SHIFTIJ = ZERO
                SHIFTIJD= ZERO
              END IF
              DO N = 1, IND(NJ+1) - IND(NJ)
                NU = IND(NJ)+N
                DO M = 1, IND(NI+1) - IND(NI)
                  MU = IND(NI)+M
                  CALL DFTB_CNVSQ(MU,NU,NSEQ)
                  DHDR = DERMAT(1,NSEQ2)
                  DSDR = DERMAT(2,NSEQ2)
                  VAL = VAL
     *                + DDER(NSEQ,1)*(DHDR+SHIFTIJ*DSDR) !! dF/da*dD/db
     *                - WDER(NSEQ,1)*DSDR !! -dS/da*dW/db
     *                + D(NSEQ,1)*SHIFTIJD*DSDR !! D*dF/dadb
                  NSEQ2 = NSEQ2 + 1
                END DO
              END DO
            END IF
          END DO
          VAL = VAL*TWO
          !! HOP CONTRIBUTION
          IF (NFG.NE.0.AND.NBDFG.NE.0)
     *      VAL = VAL + TRACEP(DDER,HOPAO(1,1,KI,NI),L1)
          IF (MASWRK.OR.NFG.NE.0) HESSIAN(IR,1) = HESSIAN(IR,1) + VAL
        END DO
      END DO
C
      RETURN
C
      END SUBROUTINE DFTB_DEDB
C
C-----------------------------------------------------------------------
C*MODULE DFTBHS    *DECK DDOT_DER
C>
C>    @brief Hessian of DFTB
C>
C>    @details Calculate vector multiplication.  Either MAT1 or MAT2
C>             should be a derivative matrix, then most matrix elements
C>             will be trivial, and this subroutine will hopefully be
C>             faster than ordinary DDOT (function).
C>
C>    @author  Yoshio Nishimoto
C>             - Sep, 2014 - Subroutine written
C>    @date    Feb, 2016 - Hiroya Nakata
C>             - Added FMO-DFTB Hessian
C>
C>           --- INPUT ---
C>    @param NAT Number of atoms
C>    @param NI Atom differentiated
C>    @param L2 L1*(L1+1)/2
C>    @param IND AO index of each atom
C>    @param MAT1,MAT2 Vecotrs to be multiplied
C>           --- OUTPUT ---
C>    @param VALUE Result of dot product
C>
C
      SUBROUTINE DDOT_DER(NAT,NI,L2,IND,VALUE,MAT1,MAT2)
C
      IMPLICIT NONE
C
      INTEGER, INTENT(IN) :: NAT,NI,L2,IND(*)
      DOUBLE PRECISION, INTENT(IN) :: MAT1(L2),MAT2(L2)
      DOUBLE PRECISION, INTENT(INOUT) :: VALUE
C
      INTEGER :: L,M,MU,NU,NSEQ,NJ
C
      IF (NI.EQ.1) GO TO 100
C
      NSEQ = IND(NI)*(IND(NI)+1)/2
      DO M = 1, IND(NI+1)-IND(NI)
        MU = IND(NI)+M
        DO NJ = 1, NI-1
          DO L = 1, IND(NJ+1) - IND(NJ)
            NU = IND(NJ)+L
            NSEQ = NSEQ + 1
            VALUE = VALUE + MAT1(NSEQ)*MAT2(NSEQ)
          END DO
        END DO
        NSEQ = NSEQ + M
      END DO
C
100   IF (NI.EQ.NAT) GO TO 200
C
      NSEQ = IND(NI+1)*(IND(NI+1)+1)/2 + IND(NI)
      DO NJ = NI+1, NAT
        DO L = 1, IND(NJ+1) - IND(NJ)
          NU = IND(NJ)+L
          DO M = 1, IND(NI+1)-IND(NI)
            MU = IND(NI)+M
            NSEQ = NSEQ + 1
            VALUE = VALUE + MAT1(NSEQ)*MAT2(NSEQ)
          END DO
          NSEQ = NSEQ + IND(NI) + NU - MU
        END DO
      END DO
C
200   RETURN
C
      END SUBROUTINE DDOT_DER
C
C-----------------------------------------------------------------------
C*MODULE DFTBHS    *DECK DFTB_DER2
C>
C>    @brief Hessian of DFTB
C>
C>    @details Calculate second-order geometrical derivatives of
C>             non-perturbed Hamiltonian and overlap matrix (AO basis)
C>
C>    @author Yoshio Nishimoto
C>    - Sep, 2014- Subroutine written
C>
C>           --- INPUT ---
C>    @param L2 L1*(L1+1)/2
C>    @param NI Atom at the begining of vector
C>    @param NJ Atom at the end of vector
C>    @param KI First axis displacement is added
C>    @param KJ Second axis displacement is added
C>    @param NAT Number of atoms
C>    @param IND AO index of each atom
C>    @param C Coordinates
C>           --- OUTPUT ---
C>    @param HDER2 2nd der. of Hamiltonian
C>    @param SDER2 2nd der. of overlap
C>
C
      SUBROUTINE DFTB_DER2(L2,NI,KI,NJ,KJ,NAT,IND,C,HDER2,SDER2,DOH2)
C
      IMPLICIT NONE
C
      INTEGER,INTENT(IN) :: L2,NI,KI,NJ,KJ,NAT,IND(*)
      DOUBLE PRECISION, INTENT(IN) :: C(3,*)
      DOUBLE PRECISION, INTENT(OUT) :: HDER2(L2),SDER2(L2)
C
      INTEGER :: NK,L,M,MU,NU,NSEQ
      DOUBLE PRECISION :: V(3),AUMM(9,9),BUMM(9,9),DHDR,DSDR
      LOGICAL :: SKP,DOH2
C
      IF (NI.EQ.NJ) THEN
        IF (DOH2) CALL SHIFT_TRIANGLE_DER_CLEAR(NAT,NI,L2,IND,HDER2)
        CALL SHIFT_TRIANGLE_DER_CLEAR(NAT,NI,L2,IND,SDER2)
      END IF
C
C     write (*,*) "ni,ki = ", ni, ki
C     write (*,*) "nj,kj = ", kj, kj
      CALL VCLR(AUMM,1,81)
      CALL VCLR(BUMM,1,81)
      IF (NI.EQ.NJ .AND. KI.EQ.KJ) THEN
        DO NK = 1, NAT
          IF (NI.EQ.NK) CYCLE
C         write (*,*) "nk = ", nk
          V(1) = C(1,NK) - C(1,NI)
          V(2) = C(2,NK) - C(2,NI)
          V(3) = C(3,NK) - C(3,NI)
          CALL DFTB_SKODE_HELPA(2,NI,NK,KI,KI,0,0,V,AUMM,BUMM,DOH2,
     *      .TRUE.,SKP)
          IF (SKP) CYCLE
          DO L = 1, IND(NK+1) - IND(NK)
            NU = IND(NK)+L
            DO M = 1, IND(NI+1) - IND(NI)
              MU = IND(NI)+M
              CALL DFTB_CNVSQ(MU,NU,NSEQ)
              DHDR = AUMM(M,L)
              DSDR = BUMM(M,L)
c             write (6,'(I3,2(X,F20.10))'),NSEQ,DHDR,DSDR
              IF(DOH2) HDER2(NSEQ) = DHDR
              SDER2(NSEQ) = DSDR
            END DO
          END DO
        END DO
      ELSE IF (NI.EQ.NJ .AND. KI.NE.KJ) THEN
        DO NK = 1, NAT
          IF (NI.EQ.NK) CYCLE
C         write (*,*) "nk = ", nk
          V(1) = C(1,NK) - C(1,NI)
          V(2) = C(2,NK) - C(2,NI)
          V(3) = C(3,NK) - C(3,NI)
          CALL DFTB_SKODE_HELPA(2,NI,NK,KI,KJ,0,0,V,AUMM,BUMM,DOH2,
     *      .TRUE.,SKP)
          IF (SKP) CYCLE
          DO L = 1, IND(NK+1) - IND(NK)
            NU = IND(NK)+L
            DO M = 1, IND(NI+1) - IND(NI)
              MU = IND(NI)+M
              CALL DFTB_CNVSQ(MU,NU,NSEQ)
              DHDR = AUMM(M,L)
              DSDR = BUMM(M,L)
C             write (*,'(I3,2(X,F20.10))'),NSEQ,DHDR,DSDR
              IF(DOH2) HDER2(NSEQ) = DHDR
              SDER2(NSEQ) = DSDR
            END DO
          END DO
        END DO
      ELSE
        V(1) = C(1,NJ) - C(1,NI)
        V(2) = C(2,NJ) - C(2,NI)
        V(3) = C(3,NJ) - C(3,NI)
        CALL DFTB_SKODE_HELPA(2,NI,NJ,KI,KJ,0,0,V,AUMM,BUMM,DOH2,
     *    .TRUE.,SKP)
        IF (SKP) THEN
          CALL VCLR(AUMM,1,81)
          CALL VCLR(BUMM,1,81)
        END IF
        DO L = 1, IND(NJ+1) - IND(NJ)
          NU = IND(NJ)+L
          DO M = 1, IND(NI+1) - IND(NI)
            MU = IND(NI)+M
            CALL DFTB_CNVSQ(MU,NU,NSEQ)
            DHDR = AUMM(M,L)
            DSDR = BUMM(M,L)
c           write (6,'("CHK=",I3,2(X,F20.10))'),NSEQ,DHDR,DSDR
            IF(DOH2) HDER2(NSEQ) = DHDR
            SDER2(NSEQ) = DSDR
          END DO
        END DO
      END IF
C
      RETURN
C
      END SUBROUTINE DFTB_DER2
C
C-----------------------------------------------------------------------
C*MODULE DFTBHS    *DECK DFTB_D2HS
C>
C>    @brief   Hessian of DFTB
C>
C>    @details Calculate Hessian contribution of second-order
C>             geometrical derivatives of non-perturbed Hamiltonian
C>             and overlap matrix
C>
C>    @author  Yoshio Nishimoto
C>             - Sep, 2014 - Subroutine written
C>    @date    Feb, 2016 - Hiroya Nakata
C>             - Added FMO-DFTB Hessian
C>    @date    Oct, 2019 - Yoshio Nishimoto
C>             - Bug fix for UDFTB parallel calculation
C>
C>           --- INPUT ---
C>    @param D Density matrix
C>    @param DEIG Energy-weighted density matrix
C>    @param SHIFT Shift contribution (atom-resolved)
C>    @param SHIFTSA Shift contribution of alpha orbitals (shell-res.)
C>    @param SHIFTSB Shift contribution of beta  orbitals (shell-res.)
C>    @param C Coordinates
C>    @param NAT Number of atoms
C>    @param L2 L1*(L1+1)/2
C>    @param NSPIN Number of spins (1=R-DFTB,2=U-DFTB)
C>    @param IND AO index of each atom
C>    @param NSPE Number of species
C>    @param ISPE Index of species
C>    @param MAXANG Maximum angular momentum of each species
C>    @param INDSAO AO index of each shell
C>    @param SCC Whether SCC or not
C>    @param SRSCC Whether shell-resolved SCC or not
C>    @param UHF Whether U-DFTB or not
C>           --- OUTPUT ---
C>    @param HDER2 2nd der. of Hamiltonian
C>    @param SDER2 2nd der. of overlap
C>    @param HESSIAN Hessian
C>
C
      SUBROUTINE DFTB_D2HS(HDER2,SDER2,D,DEIG,SHIFT,SHIFTSA,SHIFTSB,
     *  C,HESSIAN,NAT,L2,NSPIN,IND,NSHELL,NSPE,ISPE,
     *  MAXANG,INDSAO,SCC,SRSCC,UHF)
C
      IMPLICIT NONE
C
      INTEGER,INTENT(IN) :: L2,NAT,NSPIN,NSHELL
      DOUBLE PRECISION,INTENT(IN) :: D(L2,NSPIN),DEIG(L2,NSPIN),
     *  SHIFT(NAT),SHIFTSA(NSHELL),SHIFTSB(NSHELL),
     *  C(3,NAT)
      DOUBLE PRECISION,INTENT(INOUT) :: HDER2(L2),SDER2(L2),
     *  HESSIAN(3*NAT,3*NAT)
      INTEGER,INTENT(IN) :: IND(*),NSPE,ISPE(NAT),MAXANG(NSPE),INDSAO(*)
      LOGICAL,INTENT(IN) :: SCC,SRSCC,UHF
C
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      INTEGER ME,MASTER,NPROC,IBTYP,IPTIM,NSUBGR,MeUniv,NPUniv,numdlb,
     *        myworld,nworlds,mogddi
      LOGICAL GOPARR,DSKWRK,MASWRK
      COMMON /GDDI/   ISCOPE,NGROUPS,MYGROUP,MEGLOB,npglob,nnglob,JBTYP,
     *                ISGDDI,PAROUT,INITGDDI,wasgddi,MLGDDI,NSUBGR,
     *                MeUniv,NPUniv,numdlb,myworld,nworlds,mogddi
      INTEGER ISCOPE,NGROUPS,MYGROUP,MEGLOB,NPGLOB,NNGLOB,JBTYP
      LOGICAL MYJOB,ISGDDI,PAROUT,INITGDDI,WASGDDI,MLGDDI
C
      INTEGER :: NI,NJ,NK,KI,KJ,I,J,K,L,M,MU,NU,ISH,ISH0,JSH,JSH0,
     *  KSH,KSH0,NSEQ,LGROUP
      DOUBLE PRECISION :: TMP_HESS(3,3),SHIFTI,SHIFTJ,SHIFTIJ,
     *  SHIFTIA,SHIFTIB,SHIFTJA,SHIFTJB,SHIFTIJA,SHIFTIJB,VAL,
     *  TMP,D2HDR,D2SDR
C
      INTEGER, PARAMETER :: DDI_WORLD=0,DDI_GROUP=1
      DOUBLE PRECISION,PARAMETER :: ZERO=0.0D+00,TWO=2.0D+00,
     *  HALF=0.5D+00
C
      IF (ISGDDI) THEN
        CALL GDDI_SCOPE(DDI_GROUP)
        CALL GDDICOUNT(-1,LGROUP,MYJOB)
      END IF
C
      IF (UHF.OR.SRSCC) THEN
        ISH0 = 0
        DO NI = 1, NAT
          IF (ISGDDI) THEN
            CALL GDDICOUNT(0,LGROUP,MYJOB)
            IF (.NOT.MYJOB) THEN
              ISH0 = ISH0 + MAXANG(ISPE(NI))
              CYCLE
            END IF
          ELSE IF (GOPARR) THEN
            IF (MOD(NI,NPROC).NE.ME) THEN
              ISH0 = ISH0 + MAXANG(ISPE(NI))
              CYCLE
            END IF
          END IF
          JSH0 = 0
          DO NJ = 1, NI
            IF (ISGDDI.AND.GOPARR) THEN
              IF (MOD(NJ,NPROC).NE.ME) THEN
                JSH0 = JSH0 + MAXANG(ISPE(NJ))
                CYCLE
              END IF
            END IF
            DO KI = 1, 3
              DO KJ = 1, KI
                CALL DFTB_DER2(L2,NI,KI,NJ,KJ,NAT,IND,C,HDER2,SDER2,
     *               .true.)
                CALL VCLR(TMP_HESS,1,9)
                VAL = ZERO
                IF (NI.EQ.NJ) THEN
                  DO I = 1, MAXANG(ISPE(NI))
                    ISH = ISH0 + I
                    SHIFTIA = SHIFTSA(ISH)
                    SHIFTIB = SHIFTSB(ISH)
                    KSH0 = 0
                    DO NK = 1, NAT
                      DO K = 1, MAXANG(ISPE(NK))
                        KSH = KSH0 + K
                        IF (NI.EQ.NK) CYCLE
                        SHIFTJA = SHIFTSA(KSH)
                        SHIFTJB = SHIFTSB(KSH)
                        SHIFTIJA = (SHIFTIA+SHIFTJA)*HALF
                        SHIFTIJB = (SHIFTIB+SHIFTJB)*HALF
                        DO L = 1, INDSAO(KSH+1) - INDSAO(KSH)
                          NU = INDSAO(KSH)+L
                          DO M = 1, INDSAO(ISH+1) - INDSAO(ISH)
                            MU = INDSAO(ISH)+M
                            CALL DFTB_CNVSQ(MU,NU,NSEQ)
                            D2HDR = HDER2(NSEQ)
                            D2SDR = SDER2(NSEQ)
                            IF (UHF) THEN
                              TMP = D(NSEQ,1)*(D2HDR+SHIFTIJA*D2SDR)
     *                            + D(NSEQ,2)*(D2HDR+SHIFTIJB*D2SDR)
     *                            - (DEIG(NSEQ,1)+DEIG(NSEQ,2))*D2SDR
                            ELSE
                              TMP = D(NSEQ,1)*(D2HDR+SHIFTIJA*D2SDR)
     *                            - DEIG(NSEQ,1)*D2SDR
                            END IF
                            VAL = VAL + TMP
                          END DO
                        END DO
                      END DO
                      KSH0 = KSH0 + MAXANG(ISPE(NK))
                    END DO
                  END DO
                  TMP_HESS(KI,KJ) = TMP_HESS(KI,KJ) - VAL*TWO
                  IF (KI.NE.KJ)
     *              TMP_HESS(KJ,KI) = TMP_HESS(KJ,KI) - VAL*TWO
                ELSE
                  DO I = 1, MAXANG(ISPE(NI))
                    ISH = ISH0 + I
                    SHIFTIA = SHIFTSA(ISH)
                    SHIFTIB = SHIFTSB(ISH)
                    DO J = 1, MAXANG(ISPE(NJ))
                      JSH = JSH0 + J
                      SHIFTJA = SHIFTSA(JSH)
                      SHIFTJB = SHIFTSB(JSH)
                      SHIFTIJA = (SHIFTIA+SHIFTJA)*HALF
                      SHIFTIJB = (SHIFTIB+SHIFTJB)*HALF
                      DO L = 1, INDSAO(JSH+1) - INDSAO(JSH)
                        NU = INDSAO(JSH)+L
                        DO M = 1, INDSAO(ISH+1) - INDSAO(ISH)
                          MU = INDSAO(ISH)+M
                          CALL DFTB_CNVSQ(MU,NU,NSEQ)
                          D2HDR = HDER2(NSEQ)
                          D2SDR = SDER2(NSEQ)
                          IF (UHF) THEN
                            TMP = D(NSEQ,1)*(D2HDR+SHIFTIJA*D2SDR)
     *                          + D(NSEQ,2)*(D2HDR+SHIFTIJB*D2SDR)
     *                          - (DEIG(NSEQ,1)+DEIG(NSEQ,2))*D2SDR
                          ELSE
                            TMP = D(NSEQ,1)*(D2HDR+SHIFTIJA*D2SDR)
     *                          - DEIG(NSEQ,1)*D2SDR
                          END IF
                          VAL = VAL + TMP
                        END DO
                      END DO
                    END DO
                  END DO
                  TMP_HESS(KI,KJ) = TMP_HESS(KI,KJ) + VAL*TWO
                  IF (KI.NE.KJ)
     *              TMP_HESS(KJ,KI) = TMP_HESS(KJ,KI) + VAL*TWO
                END IF
                CALL ADD_HESS(NI,NJ,NAT,TMP_HESS,HESSIAN)
              END DO
            END DO
            JSH0 = JSH0 + MAXANG(ISPE(NJ))
          END DO
          ISH0 = ISH0 + MAXANG(ISPE(NI))
        END DO
      ELSE
        DO NI = 1, NAT
          IF (ISGDDI) THEN
            CALL GDDICOUNT(0,LGROUP,MYJOB)
            IF (.NOT.MYJOB) CYCLE
          ELSE IF (GOPARR) THEN
            IF (MOD(NI,NPROC).NE.ME) CYCLE
          END IF
          SHIFTI = SHIFT(NI)
          DO NJ = 1, NI
            IF (ISGDDI.AND.GOPARR) THEN
              IF (MOD(NJ,NPROC).NE.ME) CYCLE
            END IF
            DO KI = 1, 3
              DO KJ = 1, KI
                CALL DFTB_DER2(L2,NI,KI,NJ,KJ,NAT,IND,C,HDER2,SDER2,
     *               .true.)
                VAL = ZERO
                CALL VCLR(TMP_HESS,1,9)
                IF (NI.EQ.NJ) THEN
                  DO NK = 1, NAT
                    IF (NI.EQ.NK) CYCLE
                    SHIFTJ = SHIFT(NK)
                    SHIFTIJ = (SHIFTI+SHIFTJ)*HALF
                    DO L = 1, IND(NK+1) - IND(NK)
                      NU = IND(NK)+L
                      DO M = 1, IND(NI+1) - IND(NI)
                        MU = IND(NI)+M
                        CALL DFTB_CNVSQ(MU,NU,NSEQ)
                        D2HDR = HDER2(NSEQ)
                        D2SDR = SDER2(NSEQ)
                        IF (SCC) THEN
                          TMP = D(NSEQ,1)*(D2HDR+SHIFTIJ*D2SDR)
     *                        - DEIG(NSEQ,1)*D2SDR
                        ELSE
                          TMP = D(NSEQ,1)*D2HDR - DEIG(NSEQ,1)*D2SDR
                        END IF
                        VAL = VAL + TMP
                      END DO
                    END DO
                  END DO
                  TMP_HESS(KI,KJ) = TMP_HESS(KI,KJ) - VAL*TWO
                  IF (KI.NE.KJ)
     *              TMP_HESS(KJ,KI) = TMP_HESS(KJ,KI) - VAL*TWO
                ELSE
                  SHIFTJ = SHIFT(NJ)
                  SHIFTIJ = (SHIFTI+SHIFTJ)*HALF
                  DO L = 1, IND(NJ+1) - IND(NJ)
                    NU = IND(NJ)+L
                    DO M = 1, IND(NI+1) - IND(NI)
                      MU = IND(NI)+M
                      CALL DFTB_CNVSQ(MU,NU,NSEQ)
                      D2HDR = HDER2(NSEQ)
                      D2SDR = SDER2(NSEQ)
                      IF (SCC) THEN
                        TMP = D(NSEQ,1)*(D2HDR+SHIFTIJ*D2SDR)
     *                      - DEIG(NSEQ,1)*D2SDR
                      ELSE
                        TMP = D(NSEQ,1)*D2HDR - DEIG(NSEQ,1)*D2SDR
                      END IF
                      VAL = VAL + TMP
                    END DO
                  END DO
                  TMP_HESS(KI,KJ) = TMP_HESS(KI,KJ) + VAL*TWO
                  IF (KI.NE.KJ)
     *              TMP_HESS(KJ,KI) = TMP_HESS(KJ,KI) + VAL*TWO
                END IF
                CALL ADD_HESS(NI,NJ,NAT,TMP_HESS,HESSIAN)
              END DO
            END DO
          END DO
        END DO
      END IF
C
      IF (ISGDDI) THEN
        CALL GDDICOUNT( 1,LGROUP,MYJOB)
      END IF
C
C     IF (.NOT.MASWRK) THEN
C       CALL VCLR(HESSIAN,1,9*NAT*NAT)
C     END IF
C
      IF (GOPARR.OR.ISGDDI) THEN
        IF (ISGDDI) CALL GDDI_SCOPE(DDI_WORLD)
        CALL DDI_GSUMF(1,HESSIAN,9*NAT*NAT)
      END IF
C
      RETURN
C
      END SUBROUTINE DFTB_D2HS
C
C-----------------------------------------------------------------------
C*MODULE DFTBHS    *DECK DFTB_SLKODE_HELPA
C>
C>    @brief Hessian of DFTB
C>
C>    @details Helps analytical Slater-Koster transformation ...
C>
C>    @author  Yoshio Nishimoto
C>             - Sep, 2014 - Subroutine written
C>    @date    Feb, 2016 - Yoshio Nishimoto
C>             - Changes due to the DFTBSK common block
C>    @date    Feb, 2018 - Yoshio Nishimoto
C>             - Third- and fourth-order derivative
C>
C>           --- INPUT ---
C>    @param NDER Order of derivative to be calculated
C>    @param NI   Atom at the begining of vector
C>    @param NJ   Atom at the end of vector
C>    @param KI   First axis displacement is added
C>    @param KJ   Second axis displacement is added
C>    @param KK   Third axis displacement is added
C>    @param KL   Fourth axis displacement is added
C>    @param VEC  Vector
C>    @param T1   Whether Hamiltonian derivative is needed or not
C>    @param T2   Whether overlap derivative is needed or not
C>           --- OUTPUT ---
C>    @param MAT1,MAT2 NDER-th derivative of Hamiltonian and overlap
C>
C
      SUBROUTINE DFTB_SKODE_HELPA(NDER,NI,NJ,KI,KJ,KK,KL,VEC,
     *                            MAT1,MAT2,T1,T2,SKP)
      use mx_limits, only: mxatm
C
      IMPLICIT NONE

C
      INTEGER,PARAMETER :: MXSPE=10
      INTEGER, INTENT(IN) :: NDER,NI,NJ,KI,KJ,KK,KL
      DOUBLE PRECISION, INTENT(INOUT) :: VEC(3),MAT1(9,9),MAT2(9,9)
      LOGICAL, INTENT(INOUT) :: T1,T2,SKP
C
      COMMON /FMCOM / X(1)
      DOUBLE PRECISION :: X
C
      COMMON /DFTBPR/ ETEMP,DFTBDP(MXSPE*14),DAMPXHE,HUBDER(MXSPE),
     *                ZREF(MXATM),SPNCST(6,MXSPE),SPE(MXATM),NSPE,
     *                MAXANG(MXATM),ISPE(MXATM),IND(MXATM+1),IDFTBD,
     *                PARAMDIR
      DOUBLE PRECISION :: ETEMP,DFTBDP,DAMPXHE,HUBDER,ZREF,SPNCST,SPE
      INTEGER :: NSPE,MAXANG,ISPE,IND,IDFTBD
C
      COMMON /DFTBSK/ SKDIM(MXSPE,MXSPE),SKSPIN(MXSPE),QREFL(3,MXSPE),
     *                HUBBL(3,MXSPE),QREF(MXSPE),HUBB(MXSPE),SKCUT2,
     *                NEEDSK,LSKHTAB(MXSPE,MXSPE),LSKSTAB(MXSPE,MXSPE),
     *                LSKGRID(MXSPE,MXSPE),LSKSELF(MXSPE)
      DOUBLE PRECISION SKDIM,SKSPIN,QREFL,HUBBL,QREF,HUBB,SKCUT2,
     *                 PARAMDIR
      INTEGER NEEDSK,LSKHTAB,LSKSTAB,LSKGRID,LSKSELF
C
C     IF (T1) CALL VCLR(MAT1,1,81)
C     IF (T2) CALL VCLR(MAT2,1,81)
C     CALL DFTB_SKODEA(NI.EQ.NJ,ISPE(NI),ISPE(NJ),NDER,KI,KJ,KK,KL,
      CALL DFTB_SKODEA(.FALSE.,ISPE(NI),ISPE(NJ),NDER,KI,KJ,KK,KL,
     &  MAXANG(ISPE(NI)),MAXANG(ISPE(NJ)),VEC,
     &  X(LSKGRID(ISPE(NI),ISPE(NJ))),X(LSKGRID(ISPE(NJ),ISPE(NI))),
     &  SKDIM(ISPE(NI),ISPE(NJ)),SKDIM(ISPE(NJ),ISPE(NI)),
     &  X(LSKHTAB(ISPE(NI),ISPE(NJ))),X(LSKHTAB(ISPE(NJ),ISPE(NI))),
     &  X(LSKSTAB(ISPE(NI),ISPE(NJ))),X(LSKSTAB(ISPE(NJ),ISPE(NI))),
     &  X(LSKSELF(ISPE(NI))),SKCUT2,MAT1,MAT2,T1,T2,SKP)
C
      RETURN
C
      END SUBROUTINE DFTB_SKODE_HELPA
C
C-----------------------------------------------------------------------
C*MODULE DFTBHS    *DECK DFTB_EREP_HESSIAN
C>
C>    @brief Hessian of DFTB
C>
C>    @details Calculate Hessian contribution from Erep
C>
C>    @author Yoshio Nishimoto
C>    - Sep, 2014- Subroutine written
C>    @date    Feb, 2018 - Yoshio Nishimoto
C>             - Third-order derivative
C>
C>           --- INPUT ---
C>    @param NDER     Order of derivative
C>    @param DISTMAT  Matrix of inter-atomic distances
C>    @param GRADMAT  Matrix of gradient
C>    @param HESSMAT  Matrix of Hessian
C>    @param GAMMAT   Matrix of third-order derivatives
C>    @param REPINTV  Intervals of discretized repulsive energy
C>    @param REPSHORT Coefficients for short-range repulsive energy
C>    @param REPCOEFF Coefficients for repulsive energy
C>           --- OUTPUT ---
C>    @param HESSIAN  Hessian
C>
C
      SUBROUTINE DFTB_EREP_HESSIAN(NDER,DISTMAT,GRADMAT,HESSMAT,
     *  GAMMAT,HESSIAN,REPINTV,REPSHORT,REPCOEFF)
      use mx_limits, only: mxatm
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER (MXSPE=10,MXINT=150)
C
      COMMON /DFTBPR/ ETEMP,DFTBDP(MXSPE*14),DAMPXHE,HUBDER(MXSPE),
     *                ZREF(MXATM),SPNCST(6,MXSPE),SPE(MXATM),NSPE,
     *                MAXANG(MXATM),ISPE(MXATM),IND(MXATM+1),IDFTBD,
     *                PARAMDIR
      COMMON /DFTBRE/ REPCUT(MXSPE,MXSPE),NUMREP(MXSPE,MXSPE),
     *                LREPINTV(MXSPE,MXSPE),LREPSHORT(MXSPE,MXSPE),
     *                LREPCOEFF(MXSPE,MXSPE)
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
C
      DOUBLE PRECISION :: DISTMAT(*),GRADMAT(*),HESSMAT(*),GAMMAT(*),
     *  HESSIAN(3*NAT,3*NAT),REPINTV(2,MXINT,NSPE,NSPE),
     *  REPSHORT(3,NSPE,NSPE),REPCOEFF(6,MXINT,NSPE,NSPE)
      DOUBLE PRECISION, PARAMETER :: ZERO=0.0D+00,ONE=1.0D+00
C
      DO I = 1, NAT
        ISP1 = ISPE(I)
        DO J = 1, I
          ISP2 = ISPE(J)
          GRAD = ZERO
          HESS = ZERO
          IF (NDER.EQ.3) DER3 = ZERO
          CALL DFTB_CNVSQ(I,J,NSEQ)
          DIST = DISTMAT(NSEQ)
          IF (DIST.LT.1.0D-02) THEN
C         ----- VAL = 0.0 IF DISTANCE IS LESS THAN 0.01 AU
            GRAD = ZERO
            HESS = ZERO
            IF (NDER.EQ.3) DER3 = ZERO
          ELSE IF (DIST.GT.REPCUT(ISP1,ISP2)) THEN
C         ----- VAL = 0.0 IF DISTANCE IS MORE THAN CUTOFF
            GRAD = ZERO
            HESS = ZERO
            IF (NDER.EQ.3) DER3 = ZERO
          ELSE
            IF (DIST.LT.REPINTV(1,1,ISP1,ISP2)) THEN
C            --- SHORT RANGE REPULSIVE PART ---
              GRAD = -REPSHORT(1,ISP1,ISP2)
     &            * EXP(-REPSHORT(1,ISP1,ISP2)*DIST
     &                 + REPSHORT(2,ISP1,ISP2))
              HESS = (REPSHORT(1,ISP1,ISP2)**2)
     &            * EXP(-REPSHORT(1,ISP1,ISP2)*DIST
     &                 + REPSHORT(2,ISP1,ISP2))
              IF (NDER.EQ.3)
     *        DER3 = (REPSHORT(1,ISP1,ISP2)**3)
     &            * EXP(-REPSHORT(1,ISP1,ISP2)*DIST
     &                 + REPSHORT(2,ISP1,ISP2))
            ELSE
              DO K = 1, NUMREP(ISP1,ISP2)
                IF (DIST.GE.REPINTV(1,K,ISP1,ISP2)
     &              .AND. DIST.LE.REPINTV(2,K,ISP1,ISP2)) EXIT
              END DO
              X1  = DIST - REPINTV(1,K,ISP1,ISP2)
              X0  = ONE
              X00 = ONE
              X000= ONE
              IF (K.LT.NUMREP(ISP1,ISP2)) THEN
C            --- MEDIUM RANGE REPULSIVE PART ---
                DO L = 2, 4
                  GRAD = GRAD + (L-1)*REPCOEFF(L,K,ISP1,ISP2) * X0
                  HESS = HESS + (L-1)*(L-2)*REPCOEFF(L,K,ISP1,ISP2)*X00
                  IF (L.NE.2) X00 = X00 * X1
                  IF (NDER.EQ.3) THEN
                    DER3 = DER3
     *                + (L-1)*(L-2)*(L-3)*REPCOEFF(L,K,ISP1,ISP2)*X000
                  END IF
                  X0  = X0 * X1
                END DO
              ELSE
C            --- LONG RANGE REPULSIVE PART ---
                DO L = 2, 6
                  GRAD = GRAD + (L-1)*REPCOEFF(L,K,ISP1,ISP2) * X0
                  HESS = HESS + (L-1)*(L-2)*REPCOEFF(L,K,ISP1,ISP2)*X00
                  IF (L.NE.2) X00 = X00 * X1
                  IF (NDER.EQ.3) THEN
                    DER3 = DER3
     *                + (L-1)*(L-2)*(L-3)*REPCOEFF(L,K,ISP1,ISP2)*X000
                    IF (L.GE.4) X000 = X000*X1
                  END IF
                  X0  = X0 * X1
                END DO
              END IF
            END IF
          END IF
C         WRITE (*,'(I3,X,I3,X,F12.6,X,F12.6)') I,J,GRAD,HESS
          GRADMAT(NSEQ) = GRAD
          HESSMAT(NSEQ) = HESS
          IF (NDER.EQ.3) GAMMAT(NSEQ) = DER3
        END DO
      END DO
C
C     ----- ADD CONTRIBUTIONS (GRADMAT AND HESSMAT) TO HESSIAN -----
C
      IF (NDER.EQ.2)
     *CALL DFTB_MAT_HESS(NAT,C,DISTMAT,GRADMAT,HESSMAT,HESSIAN)
C     IF (NDER.EQ.3)
C    *CALL DFTB_MAT_GAMMA(NAT,C,DISTMAT,GRADMAT,HESSMAT,GAMMAT,HESSIAN)
C
      END SUBROUTINE DFTB_EREP_HESSIAN
C
C-----------------------------------------------------------------------
C*MODULE DFTBHS    *DECK ADD_HESS
C>
C>    @brief Hessian of DFTB
C>
C>    @details Add Hessian contribution to NI and NJ blocks
C>
C>    @author Yoshio Nishimoto
C>    - Sep, 2014- Subroutine written
C>
C>           --- INPUT ---
C>    @param NI,NJ Atom block to be added
C>    @param NAT Nymber of atoms
C>    @param TMP Block Hessian
C>           --- OUTPUT ---
C>    @param HESSIAN Hessian
C>
C
      SUBROUTINE ADD_HESS(NI,NJ,NAT,TMP,HESSIAN)
C
      IMPLICIT NONE
C
      INTEGER, INTENT(IN) :: NI,NJ,NAT
      DOUBLE PRECISION, INTENT(INOUT) :: TMP(3,3)
      DOUBLE PRECISION, INTENT(OUT) :: HESSIAN(3*NAT,3*NAT)
C
      INTEGER :: I,ITMP,JTMP
      DOUBLE PRECISION :: VTMP,ONE=1.0D+00
C
      IF (NI.EQ.NJ) THEN
        ITMP = 3*(NI-1)+1
        JTMP = 3*(NI-1)
        DO I = 1, 3
          JTMP = JTMP + 1
          CALL DAXPY(3,ONE,TMP(1,I),1,HESSIAN(ITMP,JTMP),1)
        END DO
      ELSE
        ITMP = 3*(NJ-1)+1
        JTMP = 3*(NI-1)
        DO I = 1, 3
          JTMP = JTMP + 1
          CALL DAXPY(3,ONE,TMP(1,I),1,HESSIAN(ITMP,JTMP),1)
        END DO
        VTMP = TMP(2,1)
        TMP(2,1) = TMP(1,2)
        TMP(1,2) = VTMP
        VTMP = TMP(3,1)
        TMP(3,1) = TMP(1,3)
        TMP(1,3) = VTMP
        VTMP = TMP(3,2)
        TMP(3,2) = TMP(2,3)
        TMP(2,3) = VTMP
        ITMP = 3*(NI-1)+1
        JTMP = 3*(NJ-1)
        DO I = 1, 3
          JTMP = JTMP + 1
          CALL DAXPY(3,ONE,TMP(1,I),1,HESSIAN(ITMP,JTMP),1)
        END DO
      END IF
C
      RETURN
C
      END SUBROUTINE ADD_HESS
C
C-----------------------------------------------------------------------
C*MODULE DFTBHS    *DECK DFTB_DISP_HESS
C>
C>    @brief Hessian of DFTB
C>
C>    @details Calculate Hessian contribution of dispersion correction.
C>
C>    @author Yoshio Nishimoto
C>    - Sep, 2014- Subroutine written
C>    @date Nov, 2014 - Yoshio Nishimoto
C>    - Added SKHP
C>    @date Apr, 2018 - Yoshio Nishimoto
C>    - Bug fix for SK and SKHP
C>
C>           --- INPUT ---
C>    @param NAT Nymber of atoms
C>    @param NSPE Number of species
C>    @param ISPE Index of species
C>    @param DFTBDP Parameters of dispersion correction
C>    @param ZREF Nuclear charges of a given atom.  In FMO-DFTB, ZREF
C>           of BDA or BAA is differet from QREF
C>    @param QREF Nuclear charges of a given species.  This is taken
C>           from Slater-Koster files
C>    @param C Coordinates
C>    @param DISTMAT Matrix of inter-atomic distances
C>    @param IDFTBD (1) UFF (2) SK
C>    @param IWRK
C>           --- OUTPUT ---
C>    @param GRADMAT Matrix of 1st-order derivatives of dispersion
C>    @param HESSMAT Matrix of 2nd-order derivatives of dispersion
C>    @param HESSIAN Hessian
C>
C
      SUBROUTINE DFTB_DISP_HESS(NAT,NSPE,ISPE,DFTBDP,ZREF,QREF,C,
     *  DISTMAT,GRADMAT,HESSMAT,HESSIAN,IDFTBD,IWRK)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DOUBLE PRECISION, PARAMETER :: ZERO=0.0D+00,ONE=1.0D+00
C
      INTEGER, INTENT(IN) :: NAT,NSPE,ISPE(*),IDFTBD
      DOUBLE PRECISION :: RDMAT(3,NSPE,NSPE)
      DOUBLE PRECISION, INTENT(IN) :: DFTBDP(*),ZREF(*),QREF(*),
     *  DISTMAT(*)
      DOUBLE PRECISION, INTENT(INOUT) :: GRADMAT(*),HESSMAT(*),
     *  HESSIAN(3*NAT,3*NAT)
      INTEGER, INTENT(INOUT) :: IWRK(*)
C
      DOUBLE PRECISION :: NEA, NEB
C
      NN = NAT*(NAT+1)/2
      CALL VCLR(GRADMAT,1,NN)
      CALL VCLR(HESSMAT,1,NN)
C
      NSP1 = 0
      NSP2 = 0
      IF (IDFTBD.EQ.1) THEN
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
          DO J = 1, I
            CALL DFTB_CNVSQ(I,J,NSEQ)
            IF (I.EQ.J) THEN
              GRADMAT(NSEQ) = ZERO
              HESSMAT(NSEQ) = ZERO
              CYCLE
            END IF
            NSP2 = ISPE(J)
            CJ = ZREF(J)/QREF(NSP2)
            GRAD = ZERO
            HESS = ZERO
            R0  = RDMAT(1,NSP1,NSP2)
            RIJ = RDMAT(2,NSP1,NSP2)
            DIJ = RDMAT(3,NSP1,NSP2)
            DIST = DISTMAT(NSEQ)
            DISTI = ONE/DIST
            DIJ = DIJ * CI * CJ
            IF (DIST.GE.R0) THEN !! ATTRACTIVE POTENTIAL
              TMP = RIJ/DIST
              R6 = TMP**6
              R12 = R6**2
              GRAD = 1.2D+01*DIJ*(R6-R12)*DISTI
              HESS = 1.2D+01*DIJ*(-7.0D+00*R6+1.3D+01*R12)*DISTI*DISTI
            ELSE !! REPULSIVE POTENTIAL
              TMP = DIST/RIJ
              R5 = TMP**5
              R10 = R5**2
              GRAD = DIJ*(-5.0D+00*C1*R5+1.0D+01*C2*R10)*DISTI
              HESS = DIJ*(-2.0D+01*C1*R5+9.0D+01*C2*R10)*DISTI*DISTI
            END IF
            GRADMAT(NSEQ) = GRAD
            HESSMAT(NSEQ) = HESS
          END DO
        END DO
      ELSE IF (IDFTBD.EQ.2) THEN
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
            C6AB = C6AB * CI * CJ
            CALL DFTB_CNVSQ(I,J,NSEQ)
            DIST = DISTMAT(NSEQ)
            DISTI = ONE/DIST
            DIST7  = DIST**7
            DIST7I = DISTI**7
            R7R07 = DIST7*R0AB7I
            VEXP = EXP(-3.0D+00*R7R07)
            DAMP1 = ONE - VEXP
            DAMP2 = DAMP1 * DAMP1
            DAMP3 = DAMP2 * DAMP1
            VAL = 8.4D+01*VEXP*R7R07-6.0D+00*DAMP1
            GRAD = -VAL*C6AB*DAMP3*DIST7I
            VAL = 2.52D+02*DAMP1*R7R07*(2.0D+00-7.0D+00*R7R07)*VEXP
     *          + 5.292D+03*R7R07*R7R07*VEXP*VEXP
     *          - 1.008D+03*DAMP1*R7R07*VEXP + 4.2D+01*DAMP2
            HESS = -VAL*C6AB*DAMP2*DIST7I*DISTI
            GRADMAT(NSEQ) = GRAD
            HESSMAT(NSEQ) = HESS
          END DO
        END DO
      ELSE IF (IDFTBD.EQ.4) THEN
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
        !! CALCULATE DISPERSION CORRECTION
        NSEQ = 0
        DO I = 1, NAT
          NSP1 = ISPE(I)
          PA   = DFTBDP(NSP1*14-12+IWRK(I))
          R0A  = DFTBDP(NSP1*14- 6+IWRK(I))
          NEA  = DFTBDP(NSP1*14)
          C6A  = 7.5D-01*SQRT(NEA*PA*PA*PA)
          CI = ZREF(I)/QREF(NSP1)
          DO J = 1, I-1 !! I+1, NAT
            NSP2 = ISPE(J)
            PB   = DFTBDP(NSP2*14-12+IWRK(J))
            R0B  = DFTBDP(NSP2*14- 6+IWRK(J))
            NEB  = DFTBDP(NSP2*14)
            C6B  = 7.5D-01*SQRT(NEB*PB*PB*PB)
            CJ = ZREF(J)/QREF(NSP2)
            C6AB = (2.0D+00*C6A*C6B*PA*PB) / (PA*PA*C6B+PB*PB*C6A)
            R0AB = (R0A**3 + R0B**3) / (R0A**2 + R0B**2)
            R0AB7I = 1.0D+00 / R0AB**7
            C6AB = C6AB * CI * CJ
            NSEQ = NSEQ + 1
            DIST = DISTMAT(NSEQ)
            DISTI = ONE/DIST
            DIST7  = DIST**7
            DIST7I = DISTI**7
            R7R07 = DIST7*R0AB7I
            VEXP = EXP(-3.0D+00*R7R07)
            DAMP1 = ONE - VEXP
            DAMP2 = DAMP1 * DAMP1
            DAMP3 = DAMP2 * DAMP1
            VAL = 8.4D+01*VEXP*R7R07-6.0D+00*DAMP1
            GRAD = -VAL*C6AB*DAMP3*DIST7I
C           VAL = 2.52D+02*R7R07*VEXP*(-2.0D+00*DAMP1
C    *                 +7.0D+00*R7R07*(3.0D+00-4.0D+00*DAMP1))
C    *          + 4.2D+01*DAMP2
C           HESS = VAL*C6AB*DAMP2*DIST7I*DISTI
            VAL = 2.52D+02*DAMP1*R7R07*(2.0D+00-7.0D+00*R7R07)*VEXP
     *          + 5.292D+03*R7R07*R7R07*VEXP*VEXP
     *          - 1.008D+03*DAMP1*R7R07*VEXP + 4.2D+01*DAMP2
            HESS = -VAL*C6AB*DAMP2*DIST7I*DISTI
            GRADMAT(NSEQ) = GRAD
            HESSMAT(NSEQ) = HESS
          END DO
          NSEQ = NSEQ + 1
        END DO
      END IF
C
C     ----- ADD CONTRIBUTIONS (GRADMAT AND HESSMAT) TO HESSIAN -----
C
      CALL DFTB_MAT_HESS(NAT,C,DISTMAT,GRADMAT,HESSMAT,HESSIAN)
C
      RETURN
C
      END SUBROUTINE DFTB_DISP_HESS
C
C-----------------------------------------------------------------------
C*MODULE DFTBHS    *DECK DFTB_MAT_HESS
C>
C>    @brief Hessian of DFTB
C>
C>    @details Calculate Hessian contribution of classical terms from
C>             1st- and 2nd-order derivatives
C>
C>    @author Yoshio Nishimoto
C>    - Sep, 2014- Subroutine written
C>
C>           --- INPUT ---
C>    @param NAT Nymber of atoms
C>    @param C Coordinates
C>    @param DISTMAT Matrix of inter-atomic distances
C>    @param GRADMAT Matrix of 1st-order derivatives of dispersion
C>    @param HESSMAT Matrix of 2nd-order derivatives of dispersion
C>           --- OUTPUT ---
C>    @param HESSIAN Hessian
C>
C
      SUBROUTINE DFTB_MAT_HESS(NAT,C,DISTMAT,GRADMAT,HESSMAT,HESSIAN)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DOUBLE PRECISION, PARAMETER :: ZERO=0.0D+00,ONE=1.0D+00
C
      INTEGER :: NAT
      DOUBLE PRECISION, INTENT(IN) :: C(3,*),DISTMAT(*),GRADMAT(*),
     *  HESSMAT(*)
      DOUBLE PRECISION, INTENT(INOUT) :: HESSIAN(3*NAT,3*NAT)
      DOUBLE PRECISION :: TMP(3,3),VEC(3)
C
C     NI :: DISPLACEMENT OF ATOM A
C     KI :: DIRECTION ALONG X, Y, OR Z (ATOM A)
C     NJ :: DISPLACEMENT OF ATOM B
C     KJ :: DIRECTION ALONG X, Y, OR Z (ATOM B)
C     NK :: SUM FOR WHEN NI = NJ
C
      DO NI = 1, NAT
        DO NJ = 1, NI
          CALL VCLR(TMP,1,9)
          IF (NI.EQ.NJ) THEN
            DO NK = 1, NAT
              IF (NI.EQ.NK) CYCLE
              CALL DFTB_CNVSQ(NI,NK,NSEQ)
              GRAD = GRADMAT(NSEQ)
              HESS = HESSMAT(NSEQ)
              IF (GRAD.EQ.ZERO.AND.HESS.EQ.ZERO) CYCLE
              DIST = DISTMAT(NSEQ)
              DISTI = ONE/DIST
              DO KK = 1, 3
                VEC(KK) = C(KK,NI) - C(KK,NK)
                VEC(KK) = VEC(KK)*DISTI
              END DO
              DO KI = 1, 3
                DO KJ = 1, KI
                  IF (KI.EQ.KJ) THEN
                    HESS_TMP = VEC(KI)*VEC(KJ)*HESS
     *                       + (ONE-VEC(KI)*VEC(KJ))*DISTI*GRAD
                    TMP(KI,KJ) = TMP(KI,KJ) + HESS_TMP
                  ELSE
                    HESS_TMP = VEC(KI)*VEC(KJ)*HESS
     *                       - VEC(KI)*VEC(KJ)*DISTI*GRAD
                    TMP(KI,KJ) = TMP(KI,KJ) + HESS_TMP
                    TMP(KJ,KI) = TMP(KJ,KI) + HESS_TMP
                  END IF
                END DO
              END DO
            END DO
          ELSE
            CALL DFTB_CNVSQ(NI,NJ,NSEQ)
            GRAD = GRADMAT(NSEQ)
            HESS = HESSMAT(NSEQ)
            IF (GRAD.EQ.ZERO.AND.HESS.EQ.ZERO) CYCLE
            DIST = DISTMAT(NSEQ)
            DISTI = ONE/DIST
            DO KK = 1, 3
              VEC(KK) = C(KK,NI) - C(KK,NJ)
              VEC(KK) = VEC(KK)*DISTI
            END DO
            DO KI = 1, 3
              DO KJ = 1, KI
                IF (KI.EQ.KJ) THEN
                  HESS_TMP = -VEC(KI)*VEC(KJ)*HESS
     *                     - (ONE-VEC(KI)*VEC(KJ))*DISTI*GRAD
                  TMP(KI,KJ) = TMP(KI,KJ) + HESS_TMP
                ELSE
                  HESS_TMP = -VEC(KI)*VEC(KJ)*HESS
     *                     + VEC(KI)*VEC(KJ)*DISTI*GRAD
                  TMP(KI,KJ) = TMP(KI,KJ) + HESS_TMP
                  TMP(KJ,KI) = TMP(KJ,KI) + HESS_TMP
                END IF
              END DO
            END DO
          END IF
          CALL ADD_HESS(NI,NJ,NAT,TMP,HESSIAN)
        END DO
      END DO
C
      RETURN
C
      END SUBROUTINE DFTB_MAT_HESS
C
C-----------------------------------------------------------------------
C*MODULE DFTBHS    *DECK DFTB_TFTRI0
C>
C>    @brief Hessian of DFTB
C>
C>    @details See below
C>
C>    @author Yoshio Nishimoto
C>    - Sep, 2014- Subroutine written
C>
C
      SUBROUTINE DFTB_TFTRI0(NDOCC,NPART,H,F,T,WRK,M,N,LDT)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      DIMENSION H(*),F(*),T(LDT,M),WRK(N)
C
      LOGICAL GOPARR,DSKWRK,MASWRK,PARALL3
C
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
C
      PARAMETER (ZERO=0.0D+00,ONE=1.0D+00)
C
C     ----- TRANSFORM THE TRIANGULAR MATRIX F USING VECTORS T -----
C                      H = T-DAGGER * F * T
C     THE ORDER OF THE TRIANGULAR MATRIX -H- IS -M-   AND   -F- IS -N-
C     THIS SUBROUTINE CALCULATES ONLY INDEPENDENT (OCCUPIED-VIRTUAL
C     COUPLING) TERMS.
C
      M2 = (M*M+M)/2
C
      PARALL3 = GOPARR  .AND.  M.GT.MXSQN3()
      IF (PARALL3) CALL VCLR(H,1,M2)
      IPCOUNT = ME - 1
C
C        THE COMPUTATION HERE IS H = T-DAGGER * (F * T),
C        WITH THE -DSPMV- FIRST PRODUCING ONE COLUMN OF F*T,
C        THEN THE -DGEMV- GENERATES AN ENTIRE ROW -J- OF -H-.
C
      DO 100 J = NDOCC+1, M
         IF(PARALL3) THEN
           IPCOUNT = IPCOUNT + 1
           IF (MOD(IPCOUNT,NPROC).NE.0) GO TO 100
         END IF
         IJ = (J*J-J)/2
         CALL DSPMV('U',N,ONE,F,T(1,J),1,ZERO,WRK,1)
         CALL DGEMV('T',M,MIN(J,NDOCC+NPART),ONE,T,LDT,WRK,1,ZERO,
     *              H(IJ+1),1)
  100 CONTINUE
C
      IF(PARALL3) CALL DDI_GSUMF(520,H,M2)
C
      RETURN
      end
C
C-----------------------------------------------------------------------
C*MODULE DFTBHS    *DECK DFTB_HESS_CHKOCC
C>
C>    @brief Hessian of DFTB
C>
C>    @details Check occupation of all orbitals in alpha and beta
C>             space.
C>
C>    @author Yoshio Nishimoto
C>    - Sep, 2014- Subroutine written
C>
C>           --- INPUT ---
C>    @param OCC Occupation number of each MO
C>    @param EIG Eigenvalues of each MO
C>    @param NSPIN Number of spins (1=R-DFTB,2=U-DFTB)
C>    @param L1 Number of basis functions
C>    @param IW 6
C>    @param MASWRK
C>    @param SOME Some output
C>    @param SCFTYP
C>    @param OCCTHR  Threshold for orbital occupation
C>           --- OUTPUT ---
C>    @param NDOCC Numer of doubly occupied orbitals in R-DFTB, but
C>           of singly occupied orbitals in U-DFTB.
C>    @param NPART Number of partially occupied orbitals
C>    @param NVIRT Number of virtual orbitals
C>    @param FON Whether fractional occupation is used or not
C>    @param ROHF Whether (same-electron) RO-DFTB or not
C>
C
      SUBROUTINE DFTB_HESS_CHKOCC(OCC,EIG,NDOCC,NPART,NVIRT,NSPIN,L1,IW,
     *  FON,ROHF,MASWRK,SOME,SCFTYP,OCCTHR)
C
      IMPLICIT NONE
C
      INTEGER, INTENT(IN) :: NSPIN,L1
      DOUBLE PRECISION, INTENT(IN) :: OCC(L1,NSPIN),EIG(L1,NSPIN),
     *                                SCFTYP,OCCTHR
      INTEGER, INTENT(OUT) :: NDOCC(2),NPART(2),NVIRT(2)
      INTEGER, INTENT(IN) :: IW
      LOGICAL, INTENT(OUT) :: FON(3),ROHF
      LOGICAL, INTENT(IN) :: MASWRK,SOME
C
      INTEGER :: IO,NHALF,ISPIN
      DOUBLE PRECISION :: DOC,OCCI
      DOUBLE PRECISION, PARAMETER :: ONE=1.0D+00,TWO=2.0D+00
      DOUBLE PRECISION RHF/8HRHF     /,ROHFC/8HROHF    /,UHF/8HUHF     /
C
      !! NHALF SHOULD BE USED, MAYBE
      !! FIRST ONLY FOR ALPHA
C
      IF (SCFTYP.EQ.RHF)   DOC = TWO
      IF (SCFTYP.EQ.UHF)   DOC = ONE
      IF (SCFTYP.EQ.ROHFC) DOC = TWO
      NDOCC(1) = 0
      NDOCC(2) = 0
      NPART(1) = 0
      NPART(2) = 0
      NVIRT(1) = 0
      NVIRT(2) = 0
C
      DO ISPIN = 1, NSPIN
        FON(ISPIN) = .FALSE.
        NHALF = 0
        DO IO = 1, L1
          OCCI = OCC(IO,ISPIN)
          IF (OCCI+OCCTHR.GT.DOC) THEN
            NDOCC(ISPIN) = NDOCC(ISPIN) + 1
          ELSE IF (OCCI.LE.OCCTHR) THEN
            NVIRT(ISPIN) = NVIRT(ISPIN) + 1
          ELSE
            NPART(ISPIN) = NPART(ISPIN) + 1
            IF (OCCI.EQ.ONE) NHALF = NHALF + 1
            IF (.NOT.FON(ISPIN).AND.MASWRK.AND.SOME) THEN
              IF (ISPIN.EQ.1) THEN
                WRITE (IW,'(/" CHECK PARTLY OCCUPIED ALPHA-ORBITALS"/)')
              ELSE
                WRITE (IW,'(/" CHECK PARTLY OCCUPIED  BETA-ORBITALS"/)')
              END IF
              WRITE (IW,*) " #ORB       EIGENVALUES(AU)      OCCUPATION
     * NUM."
            END IF
            IF (MASWRK.AND.SOME) WRITE (IW,'(2X,I4,2(X,F20.10))')
     *        IO,EIG(IO,ISPIN),OCC(IO,ISPIN)
            FON(ISPIN) = .TRUE.
          END IF
        END DO
        !! ASSUME ROHF WAVE FUNCTION, IF THE OCCUPATION OF ALL PARTLY
        !! OCCUPIED ARE ONE.
      END DO
      IF (NSPIN.EQ.1.AND.NHALF.NE.0.AND.NHALF.EQ.NPART(1)) ROHF=.TRUE.
C
      RETURN
C
      END SUBROUTINE DFTB_HESS_CHKOCC
C
C-----------------------------------------------------------------------
C*MODULE DFTBHS    *DECK DFTB_HESS_DIRECT
C>
C>    @brief Hessian of DFTB
C>
C>    @details Calculate derivative integrals in integral direct way.
C>
C>    @author Yoshio Nishimoto
C>            Sep, 2014 - Subroutine written
C>    @date   Feb, 2018 - Yoshio Nishimoto
C>            - Almost rewritten to improve performance.
C>
C>           --- INPUT ---
C>    @param HOP     HOP contributions in AO (only for FMO)
C>    @param C       Coordinates
C>    @param SHIFT   Shift contribution (atom-resolved)
C>    @param SHIFTSA Shift contribution of alpha orbitals (shell-res.)
C>    @param SHIFTSB Shift contribution of beta  orbitals (shell-res.)
C>    @param WRKH    Working array for AO derivative of Hmn
C>    @param WRKS    Working array for AO derivative of Smn
C>    @param WRK     Work space
C>    @param WRK1    Work space
C>    @param NI1     First atom to be differenciated
C>    @param NI2     Second atom to be differenciated
C>    @param KI1     First axis to be differenciated
C>    @param KI2     Second axis to be differenciated
C>    @param NAT     Number of atoms
C>    @param L0      LQMT
C>    @parma L1      Number of basis functions
C>    @param L2      L1*(L1+1)/2
C>    @param NDOCC   Number of doubly occupied orbitals
C>    @param NPART   Number of partially occupied orbitals
C>    @param IND     AO index of each atom
C>    @param NSHELL  Number of shells
C>    @param NSPE    Number of species
C>    @param NSPIN   Number of spins (1=R-DFTB,2=U-DFTB)
C>    @param INDSH   Index of shell for each atom
C>    @param INDSAO  AO index of each shell
C>    @param ISPE    Index of species
C>    @param MAXANG  Maximum angular momentum of each species
C>    @param SCC     Whether SCC or not
C>    @param SRSCC   Whether shell-resolved SCC or not
C>    @param UHF     Whether U-DFTB or not
C>    @param ARAMAN  Analytic Raman
C>    @param ADDHOP  HOP is added or not
C>    @param INDBND  Index of interacting atoms
C>    @param MAXBND  Number of maximum interacting atoms
C>           --- OUTPUT ---
C>    @param HDERMO  H_{ij}^a + V_{ij}^a
C>    @param SDERMO  S_{ij}^a
C>    @param SDERAO  S_{\mu \nu}^a
C>
C
      SUBROUTINE DFTB_HESS_DIRECT(HDERMO,SDERMO,SDERAO,HOPAO,C,VEC,
     *  SHIFT,SHIFTSA,SHIFTSB,WRKH,WRKS,WRK,WRK1,NI1,NI2,KI1,KI2,NAT,L0,
     *  L1,L2,NDOCC,NPART,IND,NSHELL,NSPE,NSPIN,INDSH,INDSAO,ISPE,
     *  MAXANG,SCC,SRSCC,UHF,ARAMAN,ADDHOP,INDBND,MAXBND)
C
      IMPLICIT NONE
C
      DOUBLE PRECISION, PARAMETER :: ZERO=0.0D+00,ONE=1.0D+00
C
      INTEGER, INTENT(IN) :: L0,L1,L2,NDOCC(2),NPART(2),NAT,NSHELL,NSPIN
      DOUBLE PRECISION, INTENT(OUT) :: HDERMO(L2,NSPIN),
     *  SDERMO(L2,NSPIN),SDERAO(L2)
      DOUBLE PRECISION, INTENT(IN) :: HOPAO(*),C(3,NAT),
     *  VEC(L1,L1,NSPIN),SHIFT(NAT),SHIFTSA(NSHELL),SHIFTSB(NSHELL)
      DOUBLE PRECISION, INTENT(INOUT) :: WRKH(L2),WRKS(L2),WRK(L1,L1),
     *  WRK1(L1,L1)
      INTEGER, INTENT(IN) :: NI1,NI2,KI1,KI2,IND(*),NSPE,
     *                       INDSH(*),INDSAO(*),ISPE(NAT),MAXANG(NSPE),
     *                       INDBND(MAXBND+1,4,*),MAXBND
      LOGICAL, INTENT(IN) :: SCC,SRSCC,UHF,ARAMAN
C
      DOUBLE PRECISION :: AU(9,9),BU(9,9),V(3)
      INTEGER :: NSEQ,NJ,M,N,MU,NU,ISPIN,NOCC,NI,KI,NAO,NNN
      LOGICAL :: SKP,ADDHOP
C
      CALL VCLR(WRKH,1,L2)
      CALL VCLR(WRKS,1,L2)
      IF (NI2.EQ.0.AND.KI2.EQ.0) THEN
      !! FO-CP
        NI = NI1
        KI = KI1
        CALL VCLR(AU,1,81)
        CALL VCLR(BU,1,81)
        DO NJ = 1, NAT
          IF (NI.EQ.NJ) CYCLE
          V(1) = C(1,NJ)-C(1,NI)
          V(2) = C(2,NJ)-C(2,NI)
          V(3) = C(3,NJ)-C(3,NI)
          CALL DFTB_SKODE_HELPA(1,NI,NJ,KI,0,0,0,V,AU,BU,.TRUE.,.TRUE.,
     *                          SKP)
          IF (SKP) CYCLE
          DO N = 1, IND(NJ+1) - IND(NJ)
            NU = IND(NJ)+N
            DO M = 1, IND(NI+1) - IND(NI)
              MU = IND(NI)+M
              CALL DFTB_CNVSQ(MU,NU,NSEQ)
              WRKH(NSEQ) = AU(M,N)
              WRKS(NSEQ) = BU(M,N)
            END DO
          END DO
        END DO
      ELSE
      !! SO-CP
        CALL DFTB_DER2(L2,NI1,KI1,NI2,KI2,NAT,IND,C,WRKH,WRKS,.TRUE.)
        IF (NI1.EQ.NI2) THEN
          CALL DSCAL(L2,-ONE,WRKH,1)
          CALL DSCAL(L2,-ONE,WRKS,1)
        END IF
        NI = NI1 !! correct?
C       write (*,*) "d2H/dadb"
C       call prtril(wrkh,l1)
C       write (*,*) "d2S/dadb"
C       call prtril(wrks,l1)
      END IF
      IF (ADDHOP) CALL DAXPY(L2,ONE,HOPAO,1,WRKH,1)
C
      !! Set parameters for sparse matrix multiplication
      NAO = MAXANG(ISPE(NI))**2
      NNN = IND(NI)+1
C
      !! SDERMO
      DO ISPIN = 1, NSPIN
        NOCC = NDOCC(ISPIN) + NPART(ISPIN)
        IF (ARAMAN) NOCC = L1
        CALL DFTB_CPYTSQ(WRKS,WRK,L1,L1)
        !! sparse matrix multiplication
        !! 1) (S^a*C) = \sum_{\nu} dS_{\mu \nu}/da * C_{\nu j}
        CALL DGEMM('N','N',L1,NOCC,NAO,ONE,WRK(1,NNN),L1,
     *             VEC(NNN,1,ISPIN),L1,ZERO,WRK1,L1)
        CALL DGEMM('N','N',NAO,NOCC,L1,ONE,WRK(NNN,1),L1,
     *             VEC(1,1,ISPIN),L1,ZERO,WRK1(NNN,1),L1)
        !! 2) C^T*(S^a*C) = \sum_{\mu} C_{\mu i} * (S^a*C)_{\mu j}
        DO N = 0, INDBND(1,1,NI)
          IF (N.EQ.0) THEN
            MU = IND(NI)+1
            NU = IND(NI+1)
            CALL DGEMM('T','N',L1,NOCC,NU-MU+1,ONE,VEC(MU,1,ISPIN),L1,
     *                 WRK1(MU,1),L1,ZERO,WRK,L1)
          ELSE
            NJ = INDBND(N+1,1,NI)
            IF (NU.NE.-1) MU = IND(NJ)+1
            NU = IND(NJ+1)
            IF (N.LT.INDBND(1,1,NI).AND.INDBND(N+2,1,NI).EQ.NJ+1) THEN
              NU = -1
              CYCLE
            END IF
            CALL DGEMM('T','N',L1,NOCC,NU-MU+1,ONE,VEC(MU,1,ISPIN),L1,
     *                 WRK1(MU,1),L1,ONE,WRK,L1)
          END IF
        END DO
        CALL CPYSQT(WRK,SDERMO(1,ISPIN),L1,1)
      END DO
C
      !! SAVE SDERAO
      IF (SCC) CALL DCOPY(L2,WRKS,1,SDERAO,1)
C
      !! HDERMO
      DO ISPIN = 1, NSPIN
        CALL DCOPY(L2,WRKH,1,WRK1,1)
        IF (SCC) THEN
          IF (UHF.OR.SRSCC) THEN
            IF (ISPIN.EQ.1) THEN
              CALL SHIFT_TRIANGLE_DER_SHELL(NAT,NI,L2,IND,INDSH,
     *          INDSAO,NSPE,ISPE,MAXANG,SHIFTSA,SDERAO,WRKS)
            ELSE
              CALL SHIFT_TRIANGLE_DER_SHELL(NAT,NI,L2,IND,INDSH,
     *          INDSAO,NSPE,ISPE,MAXANG,SHIFTSB,SDERAO,WRKS)
            END IF
          ELSE
            CALL SHIFT_TRIANGLE_DER(NAT,NI,L2,IND,SHIFT,
     *           SDERAO,WRKS)
            CALL SHIFT_TRIANGLE(NAT,L2,IND,SHIFT,SDERAO,WRKS)
          END IF
          CALL DAXPY(L2,ONE,WRKS,1,WRK1,1)
        END IF
        NOCC = NDOCC(ISPIN) + NPART(ISPIN)
        IF (ARAMAN) NOCC = L1
        CALL DFTB_CPYTSQ(WRK1,WRK,L1,L1)
C
        IF (ADDHOP) THEN
          !! If HOP contributions are added, WRK is no longer sparse
          CALL DGEMM('T','N',L1,L1,L1,ONE,VEC(1,1,ISPIN),L1,WRK,L1,
     *               ZERO,WRK1,l1)
          CALL DGEMM('N','N',L1,NOCC,L1,ONE,WRK1,L1,
     *               VEC(1,1,ISPIN),L1,ZERO,WRK,L1)
        ELSE
          !! sparse matrix multiplication
          CALL DGEMM('N','N',L1,NOCC,NAO,ONE,WRK(1,NNN),L1,
     *               VEC(NNN,1,ISPIN),L1,ZERO,WRK1,L1)
          CALL DGEMM('N','N',NAO,NOCC,L1,ONE,WRK(NNN,1),L1,
     *               VEC(1,1,ISPIN),L1,ZERO,WRK1(NNN,1),L1)
          DO N = 0, INDBND(1,1,NI)
            IF (N.EQ.0) THEN
              MU = IND(NI)+1
              NU = IND(NI+1)
              CALL DGEMM('T','N',L1,NOCC,NU-MU+1,ONE,VEC(MU,1,ISPIN),L1,
     *                   WRK1(MU,1),L1,ZERO,WRK,L1)
            ELSE
              NJ = INDBND(N+1,1,NI)
              IF (NU.NE.-1) MU = IND(NJ)+1
              NU = IND(NJ+1)
              IF (N.LT.INDBND(1,1,NI).AND.INDBND(N+2,1,NI).EQ.NJ+1) THEN
                NU = -1
                CYCLE
              END IF
              CALL DGEMM('T','N',L1,NOCC,NU-MU+1,ONE,VEC(MU,1,ISPIN),L1,
     *                   WRK1(MU,1),L1,ONE,WRK,L1)
            END IF
          END DO
        END IF
        CALL CPYSQT(WRK,HDERMO(1,ISPIN),L1,1)
      END DO
C
      END SUBROUTINE DFTB_HESS_DIRECT
C
C-----------------------------------------------------------------------
C*MODULE DFTBHS    *DECK DFTB_DISP_HESS_ESDIM
C>
C>    @brief Hessian of FMO-DFTB
C>
C>    @details Calculate dispersion Hessian of FMO-DFTB for separated
C>
C>    @author Yoshio Nishimoto
C>    - Feb, 2016- Subroutine written
C>
C>           --- INPUT ---
C>    @param IFG     I fragment index
C>    @param JFG     J fragment index
C>    @param NATI    Number of atoms in IFG
C>    @param NATJ    Number of atoms in JFG
C>    @param NC1     NATFMO*3
C>    @param NATFMO  Number of atoms in the whole system
C>    @param DISTMAT Contains distance matrix
C>    @param GRADMAT Working array for first derivative
C>    @param HESSMAT Working array for second derivative
C>    @param NBOND   Number of bonds
C>    @param INDFRG  First atom of each fragment in terms of NATFMOB?
C>    @param IATFRG  Global atom index with NBDFG?
C>    @param ISPEFMO ISPE for all atoms (NATFMO)
C>    @param DCTMP   Constants for dispersion correction, see DFTB_DCFMO
C>    @param ZREFFMO ZREF for all atoms (NATFMO)
C>    @param FMOC    Coordinates of allatoms (NATFMO)
C>    @param SCAL    Scaling factor of contributions
C>           --- IN/OUTPUT ---
C>    @param HESSIAN Hessian
C>
C
      SUBROUTINE DFTB_DISP_HESS_ESDIM(IFG,JFG,NATI,NATJ,NC1,NATFMO,
     *  DISTMAT,GRADMAT,HESSMAT,HESSIAN,NBOND,
     *  INDFRG,IATFRG,ISPEFMO,DCTMP,ZREFFMO,FMOC,SCAL)
      use mx_limits, only: mxatm

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      INTEGER,PARAMETER :: MXSPE=10
      DOUBLE PRECISION, PARAMETER :: ZERO=0.0D+00,ONE=1.0D+00
C
      COMMON /DFTBPR/ ETEMP,DFTBDP(MXSPE*14),DAMPXHE,HUBDER(MXSPE),
     *                ZREF(MXATM),SPNCST(6,MXSPE),SPE(MXATM),NSPE,
     *                MAXANG(MXATM),ISPE(MXATM),IND(MXATM+1),IDFTBD,
     *                PARAMDIR
      COMMON /DFTBSK/ SKDIM(MXSPE,MXSPE),SKSPIN(MXSPE),QREFL(3,MXSPE),
     *                HUBBL(3,MXSPE),QREF(MXSPE),HUBB(MXSPE),SKCUT2,
     *                NEEDSK,LSKHTAB(MXSPE,MXSPE),LSKSTAB(MXSPE,MXSPE),
     *                LSKGRID(MXSPE,MXSPE),LSKSELF(MXSPE)
C
      INTEGER, INTENT(IN) :: IFG,JFG,NATI,NATJ,NC1,NATFMO,
     *  NBOND(*),INDFRG(*),IATFRG(*),ISPEFMO(*)
      DOUBLE PRECISION :: RDMAT(3,NSPE,NSPE),VEC(3)
      DOUBLE PRECISION, INTENT(IN) :: DISTMAT(*),
     *  DCTMP(*),ZREFFMO(*),FMOC(3,NATFMO)
      DOUBLE PRECISION, INTENT(INOUT) :: GRADMAT(*),HESSMAT(*),
     *  HESSIAN(3*NATFMO,3*NATFMO)
C
      DOUBLE PRECISION :: NEA, NEB
C
      CALL VCLR(GRADMAT,1,NATI*NATJ)
      CALL VCLR(HESSMAT,1,NATI*NATJ)
C
      INDI = INDFRG(IFG)
      INDJ = INDFRG(JFG)
      CALL DCOPY(3*NSPE*NSPE,DCTMP,1,RDMAT,1)
      IF (IDFTBD.EQ.1) THEN
        C1 = 4.78947150872246D+01 !! 2^(5/6)*672/25
        C2 = 3.50498152274578D+01 !! 2^(2/3)*552/25
C       CR0 = 8.908987181D-01
        !! CALCULATE DISPERSION CORRECTION
        NSEQ = 0
        DO I = 1, NATI
          IAG = IATFRG(INDI+I-1)
          ISP = ISPEFMO(IAG)
          CI = ZREFFMO(INDI+I-1)/QREF(ISP)
          DO J = 1, NATJ
            NSEQ = NSEQ + 1
            JAG = IATFRG(INDJ+J-1)
            JSP = ISPEFMO(JAG)
            CJ = ZREFFMO(INDJ+J-1)/QREF(JSP)
            R0  = RDMAT(1,ISP,JSP)
            RIJ = RDMAT(2,ISP,JSP)
            DIJ = RDMAT(3,ISP,JSP)
            DIST = DISTMAT(NSEQ)
            DISTI = ONE/DIST
            DIJ = DIJ * CI * CJ
            IF (DIST.GE.R0) THEN !! ATTRACTIVE POTENTIAL
              TMPV = RIJ/DIST
              R6 = TMPV**6
              R12 = R6**2
              GRAD = 1.2D+01*DIJ*(R6-R12)*DISTI
              HESS = 1.2D+01*DIJ*(-7.0D+00*R6+1.3D+01*R12)*DISTI*DISTI
            ELSE !! REPULSIVE POTENTIAL
              TMPV = DIST/RIJ
              R5 = TMPV**5
              R10 = R5**2
              GRAD = DIJ*(-5.0D+00*C1*R5+1.0D+01*C2*R10)*DISTI
              HESS = DIJ*(-2.0D+01*C1*R5+9.0D+01*C2*R10)*DISTI*DISTI
            END IF
            GRADMAT(NSEQ) = GRAD
            HESSMAT(NSEQ) = HESS
          END DO
        END DO
      ELSE IF (IDFTBD.EQ.2) THEN
        NSEQ = 0
        DO I = 1, NATI
          IAG = IATFRG(INDI+I-1)
          ISP = ISPEFMO(IAG)
          CI = ZREFFMO(INDI+I-1)/QREF(ISP)
          DO J = 1, NATJ
            JAG = IATFRG(INDJ+J-1)
            JSP = ISPEFMO(JAG)
            CJ = ZREFFMO(INDJ+J-1)/QREF(JSP)
            C6AB   = RDMAT(1,ISP,JSP)
            R0AB7I = RDMAT(2,ISP,JSP)
            C6AB = C6AB * CI * CJ
            NSEQ = NSEQ + 1
            DIST = DISTMAT(NSEQ)
            DISTI = ONE/DIST
            DIST7  = DIST**7
            DIST7I = DISTI**7
            R7R07 = DIST7*R0AB7I
            VEXP = EXP(-3.0D+00*R7R07)
            DAMP1 = ONE - VEXP
            DAMP2 = DAMP1 * DAMP1
            DAMP3 = DAMP2 * DAMP1
            VAL = 8.4D+01*VEXP*R7R07-6.0D+00*DAMP1
            GRAD = -VAL*C6AB*DAMP3*DIST7I
            VAL = 2.52D+02*R7R07*VEXP*(-2.0D+00*DAMP1
     *                 +7.0D+00*R7R07*(3.0D+00-4.0D+00*DAMP1))
     *          + 4.2D+01*DAMP2
            HESS = VAL*C6AB*DAMP2*DIST7I*DISTI
            GRADMAT(NSEQ) = GRAD
            HESSMAT(NSEQ) = HESS
          END DO
        END DO
      ELSE IF (IDFTBD.EQ.4) THEN
        !! CALCULATE DISPERSION CORRECTION
        NSEQ = 0
        DO I = 1, NATI
          IAG = IATFRG(INDI+I-1)
          ISP = ISPEFMO(IAG)
          PA  = DFTBDP(ISP*14-12+NBOND(IAG))
          R0A = DFTBDP(ISP*14- 6+NBOND(IAG))
          NEA = DFTBDP(ISP*14)
          C6A = 7.5D-01*SQRT(NEA*PA*PA*PA)
          CI  = ZREFFMO(INDI+I-1)/QREF(ISP)
          DO J = 1, NATJ
            JAG = IATFRG(INDJ+J-1)
            JSP = ISPEFMO(JAG)
            PB  = DFTBDP(JSP*14-12+NBOND(JAG))
            R0B = DFTBDP(JSP*14- 6+NBOND(JAG))
            NEB = DFTBDP(JSP*14)
            C6B = 7.5D-01*SQRT(NEB*PB*PB*PB)
            CJ  = ZREFFMO(INDJ+J-1)/QREF(JSP)
            C6AB = (2.0D+00*C6A*C6B*PA*PB) / (PA*PA*C6B+PB*PB*C6A)
            R0AB = (R0A**3 + R0B**3) / (R0A**2 + R0B**2)
            R0AB7I = 1.0D+00 / R0AB**7
            C6AB = C6AB * CI * CJ
            NSEQ = NSEQ + 1
            DIST = DISTMAT(NSEQ)
            DISTI = ONE/DIST
            DIST7  = DIST**7
            DIST7I = DISTI**7
            R7R07 = DIST7*R0AB7I
            VEXP = EXP(-3.0D+00*R7R07)
            DAMP1 = ONE - VEXP
            DAMP2 = DAMP1 * DAMP1
            DAMP3 = DAMP2 * DAMP1
            VAL = 8.4D+01*VEXP*R7R07-6.0D+00*DAMP1
            GRAD = -VAL*C6AB*DAMP3*DIST7I
            VAL = 2.52D+02*R7R07*VEXP*(-2.0D+00*DAMP1
     *                 +7.0D+00*R7R07*(3.0D+00-4.0D+00*DAMP1))
     *          + 4.2D+01*DAMP2
            HESS = VAL*C6AB*DAMP2*DIST7I*DISTI
            GRADMAT(NSEQ) = GRAD
            HESSMAT(NSEQ) = HESS
          END DO
        END DO
      END IF
C
      IF (SCAL.NE.ONE) THEN
        CALL DSCAL(NATI*NATJ,SCAL,GRADMAT,1)
        CALL DSCAL(NATI*NATJ,SCAL,HESSMAT,1)
      END IF
C
C     ----- ADD CONTRIBUTIONS (GRADMAT AND HESSMAT) TO HESSIAN -----
C
      NSEQ = 0
      DO NI = 1, NATI
        IAG = IATFRG(INDI+NI-1)
        IGZ = (IAG-1)*3
        DO NJ = 1, NATJ
          JAG = IATFRG(INDJ+NJ-1)
          JGZ = (JAG-1)*3
          NSEQ = NSEQ + 1
          GRAD = GRADMAT(NSEQ)
          HESS = HESSMAT(NSEQ)
          IF (GRAD.EQ.ZERO.AND.HESS.EQ.ZERO) CYCLE
          DIST = DISTMAT(NSEQ)
          DISTI = ONE/DIST
          DO KK = 1, 3
            VEC(KK) = FMOC(KK,IAG) - FMOC(KK,JAG)
            VEC(KK) = VEC(KK)*DISTI
          END DO
          DO KI = 1, 3
            DO KJ = 1, 3
              VAL = VEC(KI)*VEC(KJ)*(HESS-GRAD*DISTI)
              IF (KI.EQ.KJ) VAL = VAL + GRAD*DISTI
              HESSIAN(IGZ+KI,IGZ+KJ) = HESSIAN(IGZ+KI,IGZ+KJ) + VAL
              HESSIAN(JGZ+KI,JGZ+KJ) = HESSIAN(JGZ+KI,JGZ+KJ) + VAL
              HESSIAN(IGZ+KI,JGZ+KJ) = HESSIAN(IGZ+KI,JGZ+KJ) - VAL
              HESSIAN(JGZ+KI,IGZ+KJ) = HESSIAN(JGZ+KI,IGZ+KJ) - VAL
            END DO
          END DO
        END DO
      END DO
C
      RETURN
C
      END SUBROUTINE DFTB_DISP_HESS_ESDIM
C
C-----------------------------------------------------------------------
C*MODULE DFTBHS    *DECK DFTB_DISP_HESS_ESDIM
C>
C>    @brief   DFTB Hessian
C>
C>    @details Check the degeneracy of orbitals
C>
C>    @author  Yoshio Nishimoto
C>    - Feb, 2016- Subroutine written
C>    @date   Feb, 2018 - Yoshio Nishimoto
C>            - Triply quasi-degenerated orbitals are considered
C>
C>           --- INPUT ---
C>    @param NSPIN   Number of spins
C>    @param L1      Number of atomic orbitals
C>    @param EIG     e_i
C>    @param OCC     f_i
C>    @param DEGTHR  Quasi-degenerate threshold
C>    @param NDOCC   Number of doubly occupied orbitals
C>    @param NPART   Number of partially occupied orbitals
C>    @param IW      6
C>    @param NDEG2   Number of doubly quasi-degenerate orbitals
C>    @param NDEG3   Number of triply quasi-degenerate orbitals
C>    @param MASWRK  MASWRK
C>           --- INPUT ---
C>    @param IWRK2   Index of doubly quasi-degenerate orbitals
C>    @param IWRK3   Index of triply quasi-degenerate orbitals
C>
C
      SUBROUTINE DFTB_HESS_CHKDEG(NSPIN,L1,EIG,OCC,DEGTHR,
     *                            NDOCC,NPART,IW,NDEG2,NDEG3,MASWRK,
     *                            IWRK2,IWRK3)
C
      IMPLICIT NONE
C
      INTEGER, INTENT(IN) :: NSPIN,L1,NDOCC(2),NPART(2),IW
      INTEGER, INTENT(OUT) :: NDEG2(3),NDEG3(3),IWRK2(2,*),IWRK3(3,*)
      DOUBLE PRECISION,INTENT(IN) :: DEGTHR,OCC(L1,*),EIG(L1,*)
      LOGICAL, INTENT(IN) :: MASWRK
C
      INTEGER :: IO,JO,KO,ISPIN,IDEG
      DOUBLE PRECISION :: EI,EJ,EK,OCCI,OCCJ
      LOGICAL :: FIRST
C
      FIRST = .TRUE.
      NDEG2(1) = 0
      NDEG2(2) = 0
      NDEG2(3) = 0
      NDEG3(1) = 0
      NDEG3(2) = 0
      NDEG3(3) = 0
C
      DO ISPIN = 1, NSPIN
C       DO IO = NDOCC(ISPIN)+1, NDOCC(ISPIN)+NPART(ISPIN)-1
        DO IO = 1, L1 !! temp for araman
          EI   = EIG(IO,ISPIN)
          OCCI = OCC(IO,ISPIN)
C         DO JO = IO+1, NDOCC(ISPIN)+NPART(ISPIN)
          DO JO = IO+1, L1 !! temp for araman
            EJ   = EIG(JO,ISPIN)
            OCCJ = OCC(JO,ISPIN)
C           IF (ABS(EI-EJ).LT.DEGTHR2.OR.ABS(OCCI-OCCJ).LT.DEGTHR2) THEN
            IF (ABS(EI-EJ).LT.DEGTHR) THEN
              IF (FIRST) THEN
                IF (MASWRK) THEN
                  WRITE (IW,'(" * FOLLOWING PAIRS OF ORBITALS ARE",
     *                        " DEGENERATED")')
                WRITE (IW,'("   THRESHOLD: ABS(EIG(I)-EIG(J)) < ",
     *                         E8.2)') DEGTHR
                  WRITE (IW,'(" #ORB(I)  #ORB(J)  EIG(I)-EIG(J)",
     *                        "  OCC(I)-OCC(J)")')
                  IF (NSPIN.EQ.2) WRITE (IW,'( " *** ALPHA ***")')
                END IF
                FIRST = .FALSE.
              END IF
              IF (MASWRK) WRITE (IW,'(3X,I4,5X,I4,4X,E11.4,4X,E11.4)')
     *          IO,JO,EI-EJ,OCCI-OCCJ
              NDEG2(ISPIN) = NDEG2(ISPIN) + 1
              NDEG2(3) = NDEG2(3) + 1
              IWRK2(1,NDEG2(3)) = IO
              IWRK2(2,NDEG2(3)) = JO
            END IF
          END DO
        END DO
        IF (.NOT.FIRST.AND.NSPIN.EQ.2.AND.ISPIN.EQ.1)
     *    WRITE (IW,'(/" *** BETA  ***")')
      END DO
C
      DO ISPIN = 1, NSPIN
        DO IDEG = 1, NDEG2(ISPIN)
          IO = IWRK2(1,IDEG)
          EI = EIG(IO,ISPIN)
          DO JO = IO-2, IO+2
            EJ = EIG(JO,ISPIN)
            IF (ABS(EI-EJ).GE.DEGTHR) CYCLE
            DO KO = IO-2, IO+2
              IF (IO.EQ.JO.AND.IO.EQ.KO) CYCLE
              EK = EIG(KO,ISPIN)
              IF (ABS(EI-EK).GE.DEGTHR) CYCLE
              NDEG3(ISPIN) = NDEG3(ISPIN) + 1
              NDEG3(3) = NDEG3(3) + 1
              IWRK3(1,NDEG3(3)) = IO
              IWRK3(2,NDEG3(3)) = JO
              IWRK3(3,NDEG3(3)) = KO
            END DO
          END DO
          IO = IWRK2(2,IDEG)
          EI = EIG(IO,ISPIN)
          DO JO = IO-2, IO+2
            EJ = EIG(JO,ISPIN)
            IF (ABS(EI-EJ).GE.DEGTHR) CYCLE
            DO KO = IO-2, IO+2
              IF (IO.EQ.JO.AND.IO.EQ.KO) CYCLE
              EK = EIG(KO,ISPIN)
              IF (ABS(EI-EK).GE.DEGTHR) CYCLE
              NDEG3(ISPIN) = NDEG3(ISPIN) + 1
              NDEG3(3) = NDEG3(3) + 1
              IWRK3(1,NDEG3(3)) = IO
              IWRK3(2,NDEG3(3)) = JO
              IWRK3(3,NDEG3(3)) = KO
            END DO
          END DO
        END DO
      END DO
      NDEG3(1) = 0
      NDEG3(2) = 0
      NDEG3(3) = 0
      DO ISPIN = 1, NSPIN
        DO IO = 1, L1
          EI   = EIG(IO,ISPIN)
          DO JO = IO, L1
            EJ   = EIG(JO,ISPIN)
            IF (ABS(EI-EJ).LT.DEGTHR) THEN
              DO KO = JO, L1
                IF (IO.EQ.KO) CYCLE
                EK   = EIG(KO,ISPIN)
                IF (ABS(EI-EK).LT.DEGTHR) THEN
                  NDEG3(ISPIN) = NDEG3(ISPIN) + 1
                  NDEG3(3) = NDEG3(3) + 1
                  IWRK3(1,NDEG3(3)) = IO
                  IWRK3(2,NDEG3(3)) = JO
                  IWRK3(3,NDEG3(3)) = KO
                  IF (IO.EQ.JO) THEN
                    IWRK3(1,NDEG3(3)+1) = IO
                    IWRK3(2,NDEG3(3)+1) = KO
                    IWRK3(3,NDEG3(3)+1) = IO
                    IWRK3(1,NDEG3(3)+2) = KO
                    IWRK3(2,NDEG3(3)+2) = IO
                    IWRK3(3,NDEG3(3)+2) = IO
                    NDEG3(ISPIN) = NDEG3(ISPIN) + 2
                    NDEG3(3) = NDEG3(3) + 2
                  ELSE IF (JO.EQ.KO) THEN
                    IWRK3(1,NDEG3(3)+1) = KO
                    IWRK3(2,NDEG3(3)+1) = IO
                    IWRK3(3,NDEG3(3)+1) = KO
                    IWRK3(1,NDEG3(3)+2) = KO
                    IWRK3(2,NDEG3(3)+2) = KO
                    IWRK3(3,NDEG3(3)+2) = IO
                    NDEG3(ISPIN) = NDEG3(ISPIN) + 2
                    NDEG3(3) = NDEG3(3) + 2
                  ELSE
                    IWRK3(1,NDEG3(3)+1) = IO
                    IWRK3(2,NDEG3(3)+1) = KO
                    IWRK3(3,NDEG3(3)+1) = JO
                    IWRK3(1,NDEG3(3)+2) = JO
                    IWRK3(2,NDEG3(3)+2) = IO
                    IWRK3(3,NDEG3(3)+2) = KO
                    IWRK3(1,NDEG3(3)+3) = JO
                    IWRK3(2,NDEG3(3)+3) = KO
                    IWRK3(3,NDEG3(3)+3) = IO
                    IWRK3(1,NDEG3(3)+4) = KO
                    IWRK3(2,NDEG3(3)+4) = IO
                    IWRK3(3,NDEG3(3)+4) = JO
                    IWRK3(1,NDEG3(3)+5) = KO
                    IWRK3(2,NDEG3(3)+5) = JO
                    IWRK3(3,NDEG3(3)+5) = IO
                    NDEG3(ISPIN) = NDEG3(ISPIN) + 5
                    NDEG3(3) = NDEG3(3) + 5
                  END IF
                END IF
              END DO
            END IF
          END DO
        END DO
      END DO
C
C
      IF (MASWRK) THEN
        IF (FIRST) THEN
          WRITE (IW,*) "THERE ARE NO DEGENERATE ORBITALS IN THE ",
     *                 "PARTIALLY OCCUPIED REGION"
        ELSE
          WRITE (IW,*) "IN ORDER TO PREVENT DIVERGENCE OF U MATRIX,"
          WRITE (IW,*) "MACLAURIN EXPANSION IS USED FOR THESE PAIRS"
        END IF
      END IF
C
      RETURN
C
      END SUBROUTINE DFTB_HESS_CHKDEG
C
C-----------------------------------------------------------------------
C*MODULE DFTBHS    *DECK DFTB_HESS_PCM1
C>
C>    @brief   DFTB/PCM Hessian
C>
C>    @details Calculate the derivative of shift contributions that come
C>             from PCM (SHIFT). dq/da contribution is calculated
C>             somewhere. Additionally, this calculates some of
C>             derivative contributions of ESP on tesserae (PCMTMP).
C>
C>    @author  Yoshio Nishimoto
C>    - Feb, 2016- Subroutine written
C>
C>           --- INPUT ---
C>    @param NAT    Number of atoms
C>    @param NI     Atom to be differenciated
C>    @param KI     Axis to be differenciated
C>    @param NTS    Number of tesserae
C>    @param C      Coordinate of the solute
C>    @param XCTS   X coordinate of tesserae
C>    @param YCTS   Y coordinate of tesserae
C>    @param ZCTS   Z coordinate of tesserae
C>    @param QSE    Induced charges by the pristine density matrix
C>    @param CHAMUL Mulliken charges
C>    @param ISPHE  Which atom the tessera belongs
C>           --- IN/OUTPUT ---
C>    @param SHIFT  Shift contributions
C>    @param PCMTMP Derivative contributions of ESP on tesserae
C>
C
      SUBROUTINE DFTB_HESS_PCM1(NAT,NI,KI,NTS,C,XCTS,YCTS,ZCTS,QSE,
     *                          SHIFT,CHAMUL,PCMTMP,ISPHE)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      COMMON /PCMPLY/ STOT,VOL,CCX,CCY,CCZ,RDIF,NESF,NESFP,NESFF,I_NESF,
     *                L_AST
C
      INTEGER, INTENT(IN) :: NAT,NI,KI,ISPHE(*)
      DOUBLE PRECISION, INTENT(IN) :: C(3,NAT),XCTS(NTS),YCTS(NTS),
     *  ZCTS(NTS),QSE(NTS),CHAMUL(NAT)
      DOUBLE PRECISION, INTENT(INOUT) :: SHIFT(NAT),PCMTMP(NTS)
C
      CX = C(1,NI)
      CY = C(2,NI)
      CZ = C(3,NI)
      CC = C(KI,NI)
      DO IPT = 1, NTS
        XP = XCTS(IPT)
        YP = YCTS(IPT)
        ZP = ZCTS(IPT)
        IF (KI.EQ.1) PP = XP
        IF (KI.EQ.2) PP = YP
        IF (KI.EQ.3) PP = ZP
C
        IF (ISPHE(IPT).NE.NI.OR.NESFP.GT.0) THEN
          DIST = SQRT((CX-XP)*(CX-XP)+(CY-YP)*(CY-YP)+(CZ-ZP)*(CZ-ZP))
          DISTI = 1.0D+00/DIST
          DISTI3 = DISTI*DISTI*DISTI
          !! explicit term 2
          !! negative sign, because of -S()
          SHIFT(NI) = SHIFT(NI) - QSE(IPT)*(PP-CC)*DISTI3
          !! for implicit term 2-a
          !! negative sign, because of the formulation
          !! -Dq_a*dR/da
          PCMTMP(IPT) = -CHAMUL(NI)*(PP-CC)*DISTI3
        ELSE IF (ISPHE(IPT).EQ.NI.AND.NESFP.EQ.0) THEN
C
C       DUE TO ASC DISPLACEMENT
C
          DO NJ = 1, NAT
            IF (NI.EQ.NJ) CYCLE
            CXJ = C(1,NJ)
            CYJ = C(2,NJ)
            CZJ = C(3,NJ)
            CCJ = C(KI,NJ)
            DIST = SQRT((CXJ-XP)*(CXJ-XP)+(CYJ-YP)*(CYJ-YP)
     *                 +(CZJ-ZP)*(CZJ-ZP))
            DISTI = 1.0D+00/DIST
            DISTI3 = DISTI*DISTI*DISTI
            PCMTMP(IPT) = PCMTMP(IPT) + CHAMUL(NJ)*(PP-CCJ)*DISTI3
  
            SHIFT(NJ) = SHIFT(NJ) + QSE(IPT)*(PP-CCJ)*DISTI3
          END DO
        END IF
      END DO
C
      RETURN
C
      END SUBROUTINE DFTB_HESS_PCM1
C
C-----------------------------------------------------------------------
C*MODULE DFTBHS    *DECK DFTB_QPOTEN
C>
C>    @brief   DFTB/PCM Hessian
C>
C>    @details Calculate ESP, exerted by CHAMUL, on tesserae
C>
C>    @author  Yoshio Nishimoto
C>    - Feb, 2016- Subroutine written
C>
C>           --- INPUT ---
C>    @param NAT    Number of atoms
C>    @param NTS    Number of tesserae
C>    @param C      Coordinate of the solute
C>    @param XCTS   X coordinate of tesserae
C>    @param YCTS   Y coordinate of tesserae
C>    @param ZCTS   Z coordinate of tesserae
C>           --- OUTPUT ---
C>    @param PCMTMP Contributions of ESP on tesserae
C>
C
      SUBROUTINE DFTB_QPOTEN(NAT,NTS,C,XCTS,YCTS,ZCTS,CHAMUL,PCMPOT)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      INTEGER, INTENT(IN) :: NAT,NTS
      DOUBLE PRECISION, INTENT(IN) :: C(3,NAT),XCTS(NTS),YCTS(NTS),
     *  ZCTS(NTS),CHAMUL(NTS)
      DOUBLE PRECISION, INTENT(INOUT) :: PCMPOT(NAT)
C
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      LOGICAL GOPARR,DSKWRK,MASWRK
C
      IF (GOPARR) CALL VCLR(PCMPOT,1,NTS)
      IPCOUNT = ME - 1
C
      DO ITS = 1, NTS
        IF(GOPARR) THEN
          IPCOUNT = IPCOUNT + 1
          IF(MOD(IPCOUNT,NPROC).NE.0) CYCLE
        END IF
        XP = XCTS(ITS)
        YP = YCTS(ITS)
        ZP = ZCTS(ITS)
        TMP = 0.0D+00
        DO IAT = 1, NAT
          CX = C(1,IAT)
          CY = C(2,IAT)
          CZ = C(3,IAT)
          DIST = SQRT((CX-XP)*(CX-XP)+(CY-YP)*(CY-YP)+(CZ-ZP)*(CZ-ZP))
          TMP = TMP - CHAMUL(IAT)/DIST
        END DO
        PCMPOT(ITS) = TMP
      END DO
C
      IF(GOPARR) CALL DDI_GSUMF(2418,PCMPOT,NTS)
C
      RETURN
C
      END SUBROUTINE DFTB_QPOTEN
C
C-----------------------------------------------------------------------
C*MODULE DFTBHS    *DECK DFTB_QSHIFT
C>
C>    @brief   DFTB/PCM Hessian
C>
C>    @details Calculate shift contributions due to ASCs (QSE)
C>
C>    @author  Yoshio Nishimoto
C>    - Feb, 2016- Subroutine written
C>
C>           --- INPUT ---
C>    @param NAT    Number of atoms
C>    @param NTS    Number of tesserae
C>    @param C      Coordinate of the solute
C>    @param XCTS   X coordinate of tesserae
C>    @param YCTS   Y coordinate of tesserae
C>    @param ZCTS   Z coordinate of tesserae
C>    @param QSE    Induced charges on tesserae
C>           --- IN/OUTPUT ---
C>    @param SHIFT  Shift contribution
C>
C
      SUBROUTINE DFTB_QSHIFT(NAT,NTS,C,XCTS,YCTS,ZCTS,QSE,SHIFT)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      INTEGER, INTENT(IN) :: NAT,NTS
      DOUBLE PRECISION, INTENT(IN) :: C(3,NAT),XCTS(NTS),YCTS(NTS),
     *  ZCTS(NTS),QSE(NTS)
      DOUBLE PRECISION, INTENT(INOUT) :: SHIFT(NAT)
C
      DO IAT = 1, NAT
        TMP = 0.0D+00
        CX = C(1,IAT)
        CY = C(2,IAT)
        CZ = C(3,IAT)
        DO ITS = 1, NTS
          XP = XCTS(ITS)
          YP = YCTS(ITS)
          ZP = ZCTS(ITS)
          DIST = SQRT((CX-XP)*(CX-XP)+(CY-YP)*(CY-YP)+(CZ-ZP)*(CZ-ZP))
          TMP = TMP + QSE(ITS)/DIST
        END DO
        SHIFT(IAT) = SHIFT(IAT) - TMP
      END DO
C
      RETURN
C
      END SUBROUTINE DFTB_QSHIFT
C
C-----------------------------------------------------------------------
C*MODULE DFTBHS    *DECK DFTB_HESS_PCM2
C>
C>    @brief   DFTB/PCM Hessian
C>
C>    @details Calculate remaining Hessian contributions due to PCM
C>
C>    @author  Yoshio Nishimoto
C>    - Feb, 2016- Subroutine written
C>
C>           --- INPUT ---
C>    @param NAT     Number of atoms
C>    @param NTS     Number of tesserae
C>    @param C       Coordinate of the solute
C>    @param XCTS    X coordinate of tesserae
C>    @param YCTS    Y coordinate of tesserae
C>    @param ZCTS    Z coordinate of tesserae
C>    @param QSE     Induced charges on tesserae
C>    @param QSEDER  dq/da (where q is ASC)
C>    @param QDER    dq/da (where q is Mulliken charge)
C>    @parma CHAMUL  Mulliken charges
C>    @param AS      Area of tessera
C>    @param ISPHE   Which atom the tessera belongs
C>    @param EPSHET  Epsilon of tesserae
C>    @param IDDAI   ?? used in FIXPVA
C>    @param DAI     ?? used in FIXPVA
C>    @param DCDAQ   dC/da*\bar{q} if available
C>    @param DCDAQ   dC/da*\bar{q} is stored or not
C>           --- IN/OUTPUT ---
C>    @param HESSIAN Hessian
C>
C
      SUBROUTINE DFTB_HESS_PCM2(NAT,NTS,C,XCTS,YCTS,ZCTS,QSE,QSEDER,
     *                          QDER,CHAMUL,AS,ISPHE,EPSHET,IDDAI,
     *                          DAI,PCMTMP,HESSIAN,DCDAQ,SAVDCQ)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      INTEGER, INTENT(IN) :: NAT,ISPHE(*),IDDAI(21,*)
      DOUBLE PRECISION, INTENT(IN) :: C(3,NAT),XCTS(NTS),YCTS(NTS),
     *  ZCTS(NTS),QSE(NTS),QSEDER(NTS),QDER(NAT),CHAMUL(NAT),AS(*),
     *  EPSHET(*),DAI(3,20,*),DCDAQ(NTS,3,NAT)
      DOUBLE PRECISION, INTENT(INOUT) :: PCMTMP(*),HESSIAN(3*NAT)
      LOGICAL, INTENT(IN) :: SAVDCQ
C
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      LOGICAL GOPARR,DSKWRK,MASWRK,PARALL2
      COMMON /PCMPLY/ STOT,VOL,CCX,CCY,CCZ,RDIF,NESF,NESFP,NESFF,I_NESF,
     *                L_AST
C
      DIMENSION VEC(3)
C
C     MUST BE PARALLELIZED!
C
      PARALL2 = GOPARR !.AND.  M.GT.MXSQN2()
      IPCOUNT = ME - 1
C
      DO NI = 1, NAT
        IF(PARALL2) THEN
          IPCOUNT = IPCOUNT + 1
          IF (MOD(IPCOUNT,NPROC).NE.0) CYCLE
        END IF
        IR = 3*(NI-1)+1
        DO IPT = 1, NTS
          IF (ISPHE(IPT).NE.NI.OR.NESFP.GT.0) THEN
C           Due to other ASCs
            VEC(1) = C(1,NI) - XCTS(IPT)
            VEC(2) = C(2,NI) - YCTS(IPT)
            VEC(3) = C(3,NI) - ZCTS(IPT)
            DIST = SQRT(VEC(1)*VEC(1)+VEC(2)*VEC(2)+VEC(3)*VEC(3))
            DISTI = 1.0D+00/DIST
            DISTI3 = DISTI*DISTI*DISTI
            TMP = (QSE(IPT)*QDER(NI)+QSEDER(IPT)*CHAMUL(NI))*DISTI3
            HESSIAN(IR  ) = HESSIAN(IR  ) + TMP*VEC(1)
            HESSIAN(IR+1) = HESSIAN(IR+1) + TMP*VEC(2)
            HESSIAN(IR+2) = HESSIAN(IR+2) + TMP*VEC(3)
          ELSE IF (ISPHE(IPT).EQ.NI.AND.NESFP.EQ.0) THEN
C           Due to ASC displacement
            DO NJ = 1, NAT
              IF (NI.EQ.NJ) CYCLE
              VEC(1) = C(1,NJ) - XCTS(IPT)
              VEC(2) = C(2,NJ) - YCTS(IPT)
              VEC(3) = C(3,NJ) - ZCTS(IPT)
              DIST = SQRT(VEC(1)*VEC(1)+VEC(2)*VEC(2)+VEC(3)*VEC(3))
              DISTI = 1.0D+00/DIST
              DISTI3 = DISTI*DISTI*DISTI
              TMP = (QSE(IPT)*QDER(NJ)+QSEDER(IPT)*CHAMUL(NJ))*DISTI3
              HESSIAN(IR  ) = HESSIAN(IR  ) - TMP*VEC(1)
              HESSIAN(IR+1) = HESSIAN(IR+1) - TMP*VEC(2)
              HESSIAN(IR+2) = HESSIAN(IR+2) - TMP*VEC(3)
            END DO
          END IF
        END DO
        DO KI = 1, 3
          IF (SAVDCQ) THEN
            TMP = DDOT(NTS,DCDAQ(1,KI,NI),1,QSEDER,1)
          ELSE
            CALL VCLR(PCMTMP,1,NTS)
            CALL DFTB_CALC_DCDA(NTS,NI,KI,XCTS,YCTS,ZCTS,AS,ISPHE,
     *                          EPSHET,IDDAI,DAI,QSE,PCMTMP)
            TMP = DDOT(NTS,PCMTMP,1,QSEDER,1)
          END IF
          IR = 3*(NI-1)+KI
          HESSIAN(IR) = HESSIAN(IR) + TMP
        END DO
      END DO
C
      RETURN
C
      END SUBROUTINE DFTB_HESS_PCM2
C-----------------------------------------------------------------------
C*MODULE DFTBHS    *DECK DFTB_CALC_DCDA
C>
C>    @brief   DFTB/PCM Hessian
C>
C>    @details Calculate dC/da*\bar{q}. This is taken from QCDERQ.
C>             This does not work for EFP, TD, etc.
C>
C>    @author  Yoshio Nishimoto
C>    - Feb, 2016- Subroutine written
C>
C>           --- INPUT ---
C>    @param NTS     Number of tesserae
C>    @param NI     Atom to be differenciated
C>    @param KI     Axis to be differenciated
C>    @param XCTS    X coordinate of tesserae
C>    @param YCTS    Y coordinate of tesserae
C>    @param ZCTS    Z coordinate of tesserae
C>    @param AS      Area of tessera
C>    @param ISPHE   Which atom the tessera belongs
C>    @param EPSHET  Epsilon of tesserae
C>    @param IDDAI   ?? used in FIXPVA
C>    @param DAI     ?? used in FIXPVA
C>    @param QSE     Induced charges on tesserae
C>           --- IN/OUTPUT ---
C>    @param PCMTMP  dC/da*\bar{q}
C>
C
      SUBROUTINE DFTB_CALC_DCDA(NTS,NI,KI,XCTS,YCTS,ZCTS,AS,ISPHE,
     *                          EPSHET,IDDAI,DAI,QSE,PCMTMP)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      COMMON /PCMOPT/ RABI,RASC,REFPOL,THRSLS,DENSLS,WB,WA,ETA2,GD,EVAC,
     *                RHOW,PM,AREATL,AREAKP,BONDRY,OMEGA,RET,FRO,EPSINF,
     *                EPS,DR,RSOLV,VMOL,TCE,STEN,DSTEN,CMF,TABS,IDIRCT,
     *                IPCDER,IDP,ICOMP,IFIELD,ICAV,IDISP,IPRINT,IRETCAV,
     *                ICENT,IFAST,NEVAL,IEFPOL,KEEPSM,IMGABI,IMGASC,NADD
      COMMON /PCMPLY/ STOT,VOL,CCX,CCY,CCZ,RDIF,NESF,NESFP,NESFF,I_NESF,
     *                L_AST
C
      DIMENSION XCTS(NTS),YCTS(NTS),ZCTS(NTS),AS(NTS),ISPHE(*),
     *          EPSHET(*),IDDAI(21,*),DAI(3,20,*),QSE(NTS),PCMTMP(NTS)
C
      DATA ONE/1.0D+00/
C
      !THIS IS FOR FMO?
C     IF (IAN(NI).NE.1.AND.ABS(ZAN(NI)-ONE).LT.1.0D-08) GOTO 100
C
      IF(NESFP.EQ.0)THEN
        DO ITS = 1, NTS
          IF(ISPHE(ITS).EQ.NI) THEN
            XI = XCTS(ITS)
            YI = YCTS(ITS)
            ZI = ZCTS(ITS)
            DO JTS = 1, NTS
              IF(ISPHE(JTS).NE.NI) THEN
                XJ = XCTS(JTS)
                YJ = YCTS(JTS)
                ZJ = ZCTS(JTS)
                R2 = (XJ-XI)**2+(YJ-YI)**2+(ZJ-ZI)**2
                R  = SQRT(R2)
                R3 = R*R2
                EI=EPSHET(ISPHE(ITS))/(EPSHET(ISPHE(ITS))-ONE)
                EJ=EPSHET(ISPHE(JTS))/(EPSHET(ISPHE(JTS))-ONE)
                FT = 0.5D+00*(EJ+EI)/R3
                IF (KI.EQ.1) THEN
                  PCMTMP(ITS) = PCMTMP(ITS) - (XI-XJ)*FT*QSE(JTS)
                  PCMTMP(JTS) = PCMTMP(JTS) - (XI-XJ)*FT*QSE(ITS)
                ELSE IF (KI.EQ.2) THEN 
                  PCMTMP(ITS) = PCMTMP(ITS) - (YI-YJ)*FT*QSE(JTS)
                  PCMTMP(JTS) = PCMTMP(JTS) - (YI-YJ)*FT*QSE(ITS)
                ELSE IF (KI.EQ.3) THEN
                  PCMTMP(ITS) = PCMTMP(ITS) - (ZI-ZJ)*FT*QSE(JTS)
                  PCMTMP(JTS) = PCMTMP(JTS) - (ZI-ZJ)*FT*QSE(ITS)
                END IF
              END IF
            END DO
          END IF
        END DO
      END IF
C
      IF(IPCDER.EQ.3 .AND. NESFP.EQ.0) THEN
        SQRTPI = SQRT(3.1415926535897932D+00)
        DO ITS = 1, NTS
          DO III = 1, IDDAI(21,ITS)
            IF(NI.EQ.IDDAI(III,ITS)) THEN
              DUM=-1.07D+00*SQRTPI*QSE(ITS)/SQRT(AS(ITS)**3)
              SCALE=EPSHET(ISPHE(ITS))/(EPSHET(ISPHE(ITS))-ONE)
              PCMTMP(ITS) = PCMTMP(ITS) + DUM*DAI(KI,III,ITS)*SCALE
            END IF
          ENDDO
        ENDDO
      END IF
C
      RETURN
C
      END SUBROUTINE DFTB_CALC_DCDA
C
C-----------------------------------------------------------------------
C*MODULE DFTBHS    *DECK DFTB_PRJGRD_HESS
C>
C>    @brief   DFTB/PCM Hessian
C>
C>    @details Expeimental code to project rotation mode out in PCM
C>
C>    @author Yoshio Nishimoto
C>    - Feb, 2016- Subroutine written
C>
C>           --- INPUT ---
C>    @param IMODE   Choose method?
C>    @param EGRAD   Gradient?
C>    @param KI      Axis to be differenciated
C>    @param C       Coordinate of the solute
C>    @param NAT     Number of atoms
C>    @param WRK     Working array 1
C>    @param WRK2    Working array 2
C>           --- IN/OUTPUT ---
C>    @param HESSIAN Hessian
C>
C
      SUBROUTINE DFTB_PRJGRD_HESS(IMODE,EGRAD,HESSIAN,C,RUNTYP,NAT,WRK,
     *                            WRK2)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER (ONE=1.0D+00)
C
      DIMENSION EGRAD(3,NAT),HESSIAN(3*NAT,3*NAT),C(3,NAT),WRK(*),
     *          WRK2(*)
      DIMENSION D(2)
      DATA GRADIENT/8HGRADIENT/
C
      IF (IMODE.EQ.1) THEN
        DO I = 1, 3*NAT
          CALL DFTB_PRJGRD(HESSIAN(1,I),C,RUNTYP,NAT)
        END DO
        CALL TRPOSQ(HESSIAN,3*NAT)
        DO I = 1, 3*NAT
          CALL DFTB_PRJGRD(HESSIAN(1,I),C,RUNTYP,NAT)
        END DO
      ELSE IF (IMODE.EQ.2) THEN
        NVIB   = 2 
        VIBSIZ = 1.0D-04
        D(1)   =  VIBSIZ
        D(2)   = -VIBSIZ
        SCAL   = ONE/(VIBSIZ*DBLE(NVIB))
        CALL DCOPY(3*NAT,EGRAD,1,WRK,1)
        DO NI = 1, NAT
          DO KI = 1, 3
            DO IVIB = 1, NVIB
              C(KI,NI)= C(KI,NI) + D(IVIB)
              CALL DFTB_PRJGRD(EGRAD,C,GRADIENT,NAT)
              IF (IVIB.EQ.1) CALL DCOPY(3*NAT,EGRAD,1,WRK2,1)
              IF (IVIB.EQ.2) CALL DAXPY(3*NAT,-ONE,EGRAD,1,WRK2,1)
              C(KI,NI)= C(KI,NI) - D(IVIB)
              CALL DCOPY(3*NAT,WRK,1,EGRAD,1)
            END DO
            CALL DSCAL(3*NAT,SCAL,WRK2,1)
            CALL DAXPY(3*NAT,ONE,WRK2,1,HESSIAN(1,3*(NI-1)+KI),1)
          END DO
        END DO
      END IF
C
      RETURN
C
      END SUBROUTINE DFTB_PRJGRD_HESS
C
C-----------------------------------------------------------------------
C*MODULE DFTBHS    *DECK DFTB_DERMAT
C>
C>    @brief   DFTB Hessian
C>
C>    @details Calculate derivatives of the Fock(-like) and overlap
C>             integrals in the AO and store them on memory in a sparse
C>             format.
C>
C>    @author  Yoshio Nishimoto
C>    - Feb, 2018- Subroutine written
C>
C>           --- INPUT ---
C>    @param MODE    See the comment in the code
C>    @param NAT     Number of atoms
C>    @param C       Coordinate of the solute
C>    @param INDSAO  Index of AO with respect to the index of shell
C>    @param MAXANG  Maximum angular momentum of each species
C>    @param ISPE    Index of species
C>    @param IWRK    Work space
C>    @param SKCUT2  Cutoff of Slater--Koster tables
C>    @parma SRSCC   Whether Shell-resolved SCC or not
C>    @parma UHF     Whether U-DFTB or not
C>           --- IN/OUTPUT ---
C>    @param MAXBND  Maximum number of interacting atoms
C>           --- OUTPUT ---
C>    @param LNMAT   Number of effective bonds (interacting pairs)
C>    @param INDBND  Index of interacting atoms
C>    @param DERMAT  Derivative of F and S in the AO
C>
C
      SUBROUTINE DFTB_DERMAT(MODE,NAT,C,INDSH,INDSAO,MAXANG,
     *                       ISPE,IWRK,SKCUT2,LNMAT,MAXBND,INDBND,
     *                       DERMAT,SRSCC,UHF)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION C(3,NAT),INDSH(*),INDSAO(*),MAXANG(*),ISPE(*),
     *          IWRK(NAT),INDBND(MAXBND+1,4,NAT),DERMAT(2,*)
C
      DOUBLE PRECISION  AUP(9,9),BUP(9,9),V(3)
      LOGICAL SRSCC,UHF,SKP
C
C     MODE
C     1: COUNT THE MAXIMUM NUMBER OF NEIGHBORING BONDS
C     2: MAKE NEIGHBORING BOND LIST AND CALCULATE DERIVATIVES
C
      IF (MODE.EQ.1) THEN
        !! INITIALIZE IWRK ARRAY
        DO NI = 1, NAT
          IWRK(NI) = 0
        END DO
        !! COUNT THE NUMBER OF EFFECTIVE BONDS
        LNMAT = 0
        DO NI = 1, NAT
          DO NJ = NI+1, NAT
            C1 = C(1,NI)-C(1,NJ)
            C2 = C(2,NI)-C(2,NJ)
            C3 = C(3,NI)-C(3,NJ)
            DIST2 = C1*C1 + C2*C2 + C3*C3
            IF (DIST2.LE.SKCUT2) THEN
              IWRK(NI) = IWRK(NI) + 1
              IWRK(NJ) = IWRK(NJ) + 1
              IORB = MAXANG(ISPE(NI))*MAXANG(ISPE(NI))
              JORB = MAXANG(ISPE(NJ))*MAXANG(ISPE(NJ))
              LNMAT = LNMAT + IORB*JORB*2*3
            END IF
          END DO
        END DO
C
        !! GET A MAXIMUM NUMBER
        MAXBND = 0
        DO NI = 1, NAT
          IF (IWRK(NI).GE.MAXBND) MAXBND = IWRK(NI)
        END DO
      ELSE IF (MODE.EQ.2) THEN
        !! INITIALIZE FIRST INDEX
        DO NI = 1, NAT
          INDBND(1,1,NI) = 0
        END DO
C
        !! MAKE A BOND LIST
        DO NI = 1, NAT
          DO NJ = NI+1, NAT
            C1 = C(1,NI)-C(1,NJ)
            C2 = C(2,NI)-C(2,NJ)
            C3 = C(3,NI)-C(3,NJ)
            DIST2 = C1*C1 + C2*C2 + C3*C3
            IF (DIST2.LE.SKCUT2) THEN
              INDBND(1,1,NI) = INDBND(1,1,NI) + 1
              INDBND(1,1,NJ) = INDBND(1,1,NJ) + 1
              IBND = INDBND(1,1,NI)+1
              JBND = INDBND(1,1,NJ)+1
              INDBND(IBND,1,NI) = NJ
              INDBND(JBND,1,NJ) = NI
            END IF
          END DO
        END DO
C
        NSEQ = 1
        CALL VCLR(AUP,1,81)
        CALL VCLR(BUP,1,81)
        DO NI = 1, NAT
          IORB = MAXANG(ISPE(NI))*MAXANG(ISPE(NI))
          DO KI = 1, 3
            DO IBND = 1, INDBND(1,1,NI)
              NJ = INDBND(IBND+1,1,NI)
              V(1) = C(1,NJ)-C(1,NI)
              V(2) = C(2,NJ)-C(2,NI)
              V(3) = C(3,NJ)-C(3,NI)
              CALL DFTB_SKODE_HELPA(1,NI,NJ,KI,0,0,0,V,AUP,BUP,.TRUE.,
     *             .TRUE.,SKP)
              INDBND(IBND+1,KI+1,NI) = NSEQ
              !! They should actually generate identical DERMAT
              IF (SRSCC.OR.UHF) THEN
                NN0 = INDSAO(INDSH(NI)+1)
                MM0 = INDSAO(INDSH(NJ)+1)
                DO JSH0 = 1, MAXANG(ISPE(NJ))
                  JSH = INDSH(NJ)+JSH0
                  DO MM = INDSAO(JSH)+1, INDSAO(JSH+1)
                    M = MM - MM0
                    DO ISH0 = 1, MAXANG(ISPE(NI))
                      ISH = INDSH(NI)+ISH0
                      DO NN = INDSAO(ISH)+1, INDSAO(ISH+1)
                        N = NN - NN0
                        DERMAT(1,NSEQ) = AUP(N,M)
                        DERMAT(2,NSEQ) = BUP(N,M)
                        NSEQ = NSEQ + 1
                      END DO
                    END DO
                  END DO
                END DO
              ELSE
                JORB = MAXANG(ISPE(NJ))*MAXANG(ISPE(NJ))
                DO M = 1, JORB !! IND(NJ+1) - IND(NJ)
                  DO N = 1, IORB !! IND(NI+1) - IND(NI)
                    DERMAT(1,NSEQ) = AUP(N,M)
                    DERMAT(2,NSEQ) = BUP(N,M)
                    NSEQ = NSEQ + 1
                  END DO
                END DO
              END IF
            END DO
          END DO
        END DO
        IF (NSEQ-1.GT.LNMAT) THEN
          write (*,*) "not good in dftb_dermat"
          call abrt
        END IF
      END IF
C
      RETURN
C
      END SUBROUTINE DFTB_DERMAT
C
C-----------------------------------------------------------------------
C*MODULE DFTBHS    *DECK DFTB_EIGINV
C>
C>    @brief   DFTB Hessian
C>
C>    @details Calculate 1/(ei-ej) factor for orbital rotation.
C>
C>    @author  Yoshio Nishimoto
C>    - Feb, 2018- Subroutine written
C>    @date    Mar, 2020 - Yoshio Nishimoto
C>             - Fixed the -fpe0 issue
C>
C>           --- INPUT ---
C>    @param MODE    1: For orbital rotation
C>                   2: For third-order derivatives (antisymmetric)
C>    @param L1      Number of molecular orbitals
C>    @param NDOCC   Number of doubly occupied orbitals
C>    @param NPART   Number of partially occupied orbitals
C>    @param MAXVIR  Number of partially occupied + virtual orbitals
C>    @param MAXOCC  Number of occupied + partially occupied orbitals
C>    @param EIG     Eigenvalue of each MO
C>    @param OCC     Occupation number of each MO
C>    @param DEGTHR  Degeneracy threshold
C>    @param BETA    1/kTe
C>    @param FREQ    Flag for frequency-dependent derivative
C>    @param OMEGA   Frequency of the electric field
C>           --- OUTPUT ---
C>    @param EIGINV  1/(ei-ej)
C>
C
      SUBROUTINE DFTB_EIGINV(MODE,L1,NSPIN,NDOCC,NPART,MAXVIR,MAXOCC,
     *                       EIGINV,EIG,OCC,DEGTHR,BETA,FREQ,OMEGA)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER (HALF=0.5D+00,ONE=1.0D+00)
C
      DIMENSION NDOCC(2),NPART(2),EIGINV(MAXVIR,MAXOCC,*),EIG(L1,*),
     *          OCC(L1,*)
      LOGICAL   FREQ,PARTI,PARTM
      DIMENSION EIGTMP(2),OCCTMP(2)
C
      CALL VCLR(EIGINV,1,MAXVIR*MAXOCC*NSPIN)
      DO ISPIN = 1, NSPIN
        DO IO = 1, MAXOCC
          EI = EIG(IO,ISPIN)
          OCCI = OCC(IO,ISPIN)
          IF (MODE.EQ.1) MO0 = NDOCC(ISPIN)+1
          IF (MODE.EQ.2) MO0 = 1
          IF (MODE.EQ.2.AND.NSPIN.EQ.1) OCCI = OCCI*HALF
          PARTI = .FALSE.
          IF (IO.GE.NDOCC(ISPIN)+1.AND.IO-NDOCC(ISPIN).LE.NPART(ISPIN))
     *       PARTI = .TRUE.
          DO MO = MO0, L1
            EM = EIG(MO,ISPIN)
            PARTM = .FALSE.
            IF (MO.GE.NDOCC(ISPIN)+1.AND.
     *        MO-NDOCC(ISPIN).LE.NPART(ISPIN)) PARTM = .TRUE.
            IF (MODE.EQ.1) THEN
              IF (MO.LT.IO) CYCLE
              IF (IO.EQ.MO) THEN
                SCAL = 0.0D+00
              ELSE IF (PARTM) THEN
                IF (ABS(EI-EM).GE.DEGTHR) THEN
C                 SCAL = ONE/((EI-EM)*OCCI)
                  SCAL = (OCCI-OCC(MO,ISPIN))/((EI-EM)*OCCI)
                ELSE
C                 SCAL = ONE/OCCI
                  EIGTMP(1) = EI
                  EIGTMP(2) = EM
                  OCCTMP(1) = OCCI
                  IF (NSPIN.EQ.1) OCCTMP(1)=OCCTMP(1)*HALF
                  OCCTMP(2) = OCC(MO,ISPIN)
                  IF (NSPIN.EQ.1) OCCTMP(2)=OCCTMP(2)*HALF
                  SCAL = DFTB_CLCDEG(1,0,0,BETA,EIGTMP,OCCTMP)/OCCI
                END IF
              ELSE
                IF (ABS(EI-EM).LE.1.0D-12) THEN
                  SCAL = 0.0D+00
                ELSE
                  SCAL = ONE/(EI-EM)
                END IF
              END IF
              EIGINV(MO-NDOCC(ISPIN),IO,ISPIN) = SCAL
              IF (FREQ) THEN
                !! (U(omega)+U(-omega))/2
                OCCM = OCC(MO,ISPIN)
                EIGINV(MO-NDOCC(ISPIN),IO,ISPIN) = (EI-EM)*(OCCI-OCCM)
     *            /(((EI-EM)*(EI-EM)-OMEGA*OMEGA)*OCCI)
                IF (IO.EQ.MO) EIGINV(MO-NDOCC(ISPIN),IO,ISPIN) = 0.0D+00
              END IF
            ELSE IF (MODE.EQ.2) THEN
              IF (PARTI.OR.PARTM) THEN
                IF (ABS(EI-EM).GE.DEGTHR) THEN
                  SCAL = ONE/(EI-EM)
                ELSE
C                 SCAL = -OCCI*(ONE-OCCI)*BETA
                  EIGTMP(1) = EI
                  EIGTMP(2) = EM
                  OCCTMP(1) = OCCI
                  OCCTMP(2) = OCC(MO,ISPIN)
                  IF (NSPIN.EQ.1) OCCTMP(2)=OCCTMP(2)*HALF
                  SCAL = DFTB_CLCDEG(1,0,0,BETA,EIGTMP,OCCTMP)
                END IF
              ELSE
                IF (IO.EQ.MO) THEN
                  SCAL = 0.0D+00
                ELSE
                  SCAL = ONE/(EI-EM)
                END IF
              END IF
              EIGINV(MO,IO,ISPIN) = SCAL
            END IF
          END DO
        END DO
      END DO
C
      RETURN
C
      END SUBROUTINE DFTB_EIGINV
C
C-----------------------------------------------------------------------
C*MODULE DFTBHS    *DECK DFTB_KI12R
C>
C>    @brief   DFTB third-order derivatives
C>
C>    @details Identify the index of the polarizability from KI1 and KI2
C>             (= x ,y, or z) 
C>
C>    @author  Yoshio Nishimoto
C>    - Feb, 2018- Subroutine written
C>
C>           --- INPUT ---
C>    @param KI1     The first index of the polarizability (x, y, or z)
C>    @param KI2     The second index of the polarizability
C>           --- OUTPUT ---
C>    @param KI      The index of the polarizability
C>
C
      SUBROUTINE DFTB_KI12R(KI1,KI2,KI)
      INTEGER KI1,KI2,KI,KI1T,KI2T
      KI1T = MIN(KI1,KI2)
      KI2T = MAX(KI1,KI2)
      IF (KI1T.EQ.1.AND.KI2T.EQ.1) KI = 1
      IF (KI1T.EQ.1.AND.KI2T.EQ.2) KI = 2
      IF (KI1T.EQ.2.AND.KI2T.EQ.2) KI = 4 ! 3
      IF (KI1T.EQ.1.AND.KI2T.EQ.3) KI = 3 ! 4
      IF (KI1T.EQ.2.AND.KI2T.EQ.3) KI = 5
      IF (KI1T.EQ.3.AND.KI2T.EQ.3) KI = 6
      END SUBROUTINE DFTB_KI12R
C
C-----------------------------------------------------------------------
C*MODULE DFTBHS    *DECK DFTB_CALC_POLAR
C>
C>    @brief   Polarizability in the dipole approximation
C>
C>    @details Calculate n-th order electric field derivative of the
C>             total energy using the (n-1)-th order derivatives of the
C>             Mulliken charge. Note that this is valid only when the
C>             dipole approximation holds true.
C>
C>    @author  Yoshio Nishimoto
C>    - Feb, 2018- Subroutine written
C>
C>           --- INPUT ---
C>    @param QFD     (n-1)-th order Mulliken charge derivative
C>    @param C       Coordinates
C>    @param NAT     Number of atoms
C>    @param NDIMSH  If (SRSCC); NDIMSH=NSHELL; else NDIMSH=NAT
C>    @param NSPIN   Number of spins (1=R-DFTB,2=U-DFTB)
C>    @param MAXANG  Maximum angular momentum of each species
C>    @param ISPE    Index of species
C>    @param SRSCC   Whether Shell-resolved SCC or not
C>    @param FREQ    Flag for frequency-dependent derivative
C>    @param HRAMAN  Fourh-order is computing?
C>           --- OUTPUT ---
C>    @param ALPHA   n-th order electric field derivative
C>
C
      SUBROUTINE DFTB_CALC_POLAR(ALPHA,QFD,C,NAT,NDIMSH,NSPIN,
     *                           MAXANG,ISPE,SRSCC,FREQ,HRAMAN)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION ALPHA(6),QFD(NDIMSH,NSPIN,*),C(3,*),MAXANG(*),ISPE(*)
      LOGICAL SRSCC,FREQ,HRAMAN
C
      IF (.NOT.SRSCC) THEN
        DO KI = 1, 6
          CALL DFTB_KI12(KI,KI1,KI2)
          IF (HRAMAN.AND.FREQ) THEN
            IF (KI1.EQ.2) KI1 = 4
            IF (KI1.EQ.3) KI1 = 7
          ELSE IF (FREQ) THEN
            IF (KI1.EQ.3) KI1 = 5
            IF (KI1.EQ.2) KI1 = 3
          END IF
          ALPHA(KI) = -DDOT(NAT,QFD(1,1,KI1),1,C(KI2,1),3)
        END DO
      ELSE IF (SRSCC) THEN
        DO KI = 1, 6
          CALL DFTB_KI12(KI,KI1,KI2)
          IF (HRAMAN.AND.FREQ) THEN
            IF (KI1.EQ.2) KI1 = 4
            IF (KI1.EQ.3) KI1 = 7
          ELSE IF (FREQ) THEN
            IF (KI1.EQ.3) KI1 = 5
            IF (KI1.EQ.2) KI1 = 3
          END IF
          TMP = 0.0D+00
          ISH0 = 0
          DO NI = 1, NAT
            DO II = 1, MAXANG(ISPE(NI))
              ISH = ISH0 + II
              TMP = TMP + QFD(ISH,1,KI1)*C(KI2,NI)
            END DO
            ISH0 = ISH0 + MAXANG(ISPE(NI))
          END DO
          ALPHA(KI) = -TMP
        END DO
      END IF
C
      END SUBROUTINE DFTB_CALC_POLAR
C
C-----------------------------------------------------------------------
C*MODULE DFTBHS    *DECK DFTB_KI12
C>
C>    @brief   DFTB third-order derivatives
C>
C>    @details Identify the first and second index of the polarizability
C>             from an index
C>
C>    @author  Yoshio Nishimoto
C>    - Feb, 2018- Subroutine written
C>
C>           --- INPUT ---
C>    @param KI      The index of the polarizability
C>           --- OUTPUT ---
C>    @param KI1     The first index of the polarizability (x, y, or z)
C>    @param KI2     The second index of the polarizability
C>
C
      SUBROUTINE DFTB_KI12(KI,KI1,KI2)
C
      INTEGER KI,KI1,KI2
C
      IF (KI.EQ.1) THEN
        KI1 = 1
        KI2 = 1
      ELSE IF (KI.EQ.2) THEN
        KI1 = 1
        KI2 = 2
      ELSE IF (KI.EQ.3) THEN
        KI1 = 1
        KI2 = 3
      ELSE IF (KI.EQ.4) THEN
        KI1 = 2
        KI2 = 2
      ELSE IF (KI.EQ.5) THEN
        KI1 = 2
        KI2 = 3
      ELSE IF (KI.EQ.6) THEN
        KI1 = 3
        KI2 = 3
      END IF
      END SUBROUTINE DFTB_KI12
C
C-----------------------------------------------------------------------
C*MODULE DFTBHS    *DECK DFTB_CALC_HYPER
C>
C>    @brief   Calculate hyperpolarizability
C>
C>    @details Calculate hyperpolarizabilities using Wigner's 2n+1 rule.
C>             The third-order derivatives are computed with first-order
C>             response. No output. This only prints out the results.
C>
C>    @author  Yoshio Nishimoto
C>    - Feb, 2018- Subroutine written
C>
C>           --- INPUT ---
C>    @param EIG     Eigenvalues of each orbital
C>    @param OCC     Occupation number of each orbital
C>    @param IW      6
C>    @param NAT     Number of atoms
C>    @param L1      Number of basis functions
C>    @param L2      L1*(L1+1)/2
C>    @param NSPIN   Number of spins (1=R-DFTB,2=U-DFTB)
C>    @param NDOCC   Number of doubly occupied orbitals
C>    @param NPART   Number of partially occupied orbitals
C>    @param FON     Whether fractional occupation is used or not
C>    @param BETA    1/kTe
C>    @param EDER    Derivatives of the orbital energy
C>    @param FDER    Derivatives of the Fock matrix in MO wrt electric f
C>    @param FDER    Derivatives of the Fermi energy
C>    @param FDERFSC Square and scaled FDER
C>    @param NDEG2   Number of doubly quasi-degenerated pairs
C>    @param NDEG3   Number of triply quasi-degenerated pairs
C>    @param LSTDEG2 Index of doubly quasi-degenerated pairs
C>    @param LSTDEG3 Index of triply quasi-degenerated pairs
C>    @param WRK1    Work space
C>    @param WRK2    Work space
C>    @param WRK3    Work space
C>    @param DEGTHR  Degeneracy threshold
C>    @param DFTB3   Whether DFTB3 or not
C>    @param GAMMA3  Matrix of DFTB3 gamma
C>    @param QDERF   Mulliken charge derivatives wrt electric field
C>    @param NDIMSH  If (SRSCC); NDIMSH=NSHELL; else NDIMSH=NAT
C>    @param SRSCC   Whether Shell-resolved SCC or not
C>    @param INDSH   Index of shell for each atom
C>    @param ISPE    Index of species
C>    @parma MAXANG  Maximum angular momentum of each species
C>    @param FREQ    Flag for frequency-dependent derivative
C>    @param HRAMAN  Fourh-order is computing?
C>
C
C     Loss of trailing digits occurs
      SUBROUTINE DFTB_CALC_HYPER(OCC,EIG,IW,NAT,L1,L2,NSPIN,NDOCC,NPART,
     *                           FON,BETA,EDER,FDER,FEDER,FDERFSC,
     *                           NDEG2,NDEG3,LSTDEG2,LSTDEG3,
     *                           WRK1,WRK2,WRK3,DEGTHR,
     *                           DFTB3,GAMMA3,QDERF,NDIMSH,SRSCC,INDSH,
     *                           ISPE,MAXANG,FREQ,HRAMAN)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER (ZERO=0.0D+00,ONE=1.0D+00,TWO=2.0D+00,HALF=0.5D+00)
      INTEGER, PARAMETER :: DDI_WORLD=0,DDI_GROUP=1
C
      COMMON /GDDI/   ISCOPE,NGROUPS,MYGROUP,MEGLOB,npglob,nnglob,JBTYP,
     *                ISGDDI,PAROUT,INITGDDI,wasgddi,MLGDDI,NSUBGR,
     *                MeUniv,NPUniv,numdlb,myworld,nworlds,mogddi
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
C
      DIMENSION NDOCC(2),NPART(2),OCC(L1,NSPIN),EIG(L1,NSPIN),
     *          EDER(L1,NSPIN,3),FDER(L2,NSPIN,3),FEDER(4),
     *          FDERFSC(L1,L1,NSPIN,3),LSTDEG2(2,*),LSTDEG3(3,*),
     *          NDEG2(3),NDEG3(3),
     *          WRK1(L1,*),WRK2(L1,*),WRK3(L1,*),
     *          GAMMA3(NDIMSH,NDIMSH),QDERF(NDIMSH,NSPIN,3),INDSH(*),
     *          ISPE(*),MAXANG(*)
      LOGICAL FON(3),DFTB3,SRSCC,FREQ,HRAMAN
      LOGICAL MYJOB,ISGDDI,PAROUT,INITGDDI,WASGDDI,MLGDDI
      LOGICAL GOPARR,DSKWRK,MASWRK
C
      DIMENSION VCHARA(3),HYPER(10),EIGTMP(3),OCCTMP(3)
      DATA VCHARA/8HX       ,8HY       ,8HZ       /
C
      IF (ISGDDI) CALL GDDI_SCOPE(DDI_WORLD)
C
      ISPIN = 1
      L3 = L1*L1
      THRES=DEGTHR
      NOCC = NDOCC(ISPIN) + NPART(ISPIN)
      CALL VCLR(HYPER,1,10)
      IHYPER = 0
C
      IF (ISGDDI) THEN
        CALL GDDI_SCOPE(DDI_GROUP)
        CALL GDDICOUNT(-1,LGROUP,MYJOB)
      END IF
C
      DO KI1 = 1, 3
      DO KI2 = 1, KI1
      DO KI3 = 1, KI2
      IHYPER = IHYPER + 1
      IF (ISGDDI) THEN
        call GDDICOUNT(0,LGROUP,MYJOB)
        IF (.NOT.MYJOB) CYCLE
      END IF
      CALL DFTB_KI12R(KI1,KI2,KI4)
      CALL DFTB_KI12R(KI2,KI3,KI5)
      CALL DFTB_KI12R(KI3,KI1,KI6)
C
      VAL=ZERO
C
      !! loop (a) and (b) and (d)
      DO KI0 = 1, 3
        IF (KI0.EQ.1) THEN
          KI4 = KI1
          KI5 = KI2
          KI6 = KI3
        ELSE IF (KI0.EQ.2) THEN
          KI4 = KI2
          KI5 = KI3
          KI6 = KI1
        ELSE IF (KI0.EQ.3) THEN
          KI4 = KI3
          KI5 = KI1
          KI6 = KI2
        END IF
        IF (FREQ.AND.HRAMAN) THEN
          IF (KI4.EQ.1) KI4A = 1
          IF (KI4.EQ.2) KI4A = 4
          IF (KI4.EQ.3) KI4A = 7
C         IF (KI5.EQ.1) KI5A = 1
C         IF (KI5.EQ.2) KI5A = 4
C         IF (KI5.EQ.3) KI5A = 7
          IF (KI6.EQ.1) KI6A = 1
          IF (KI6.EQ.2) KI6A = 4
          IF (KI6.EQ.3) KI6A = 7
        ELSE IF (FREQ) THEN
          IF (KI4.EQ.1) KI4A = 1
          IF (KI4.EQ.3) KI4A = 5
          IF (KI4.EQ.2) KI4A = 3
C         IF (KI5.EQ.1) KI5A = 1
C         IF (KI5.EQ.3) KI5A = 5
C         IF (KI5.EQ.2) KI5A = 3
          IF (KI6.EQ.1) KI6A = 1
          IF (KI6.EQ.3) KI6A = 5
          IF (KI6.EQ.2) KI6A = 3
        ELSE
          KI4A = KI4
C         KI5A = KI5
          KI6A = KI6
        END IF
        DO ISPIN = 1, NSPIN
          NOCC = NDOCC(ISPIN) + NPART(ISPIN)
          CALL DFTB_POLDER1(1,NSPIN,ISPIN,L1,L2,L3,NDOCC(ISPIN),NOCC,
     *                      NDEG2,LSTDEG2,VAL,OCC(1,ISPIN),
     *                      EIG(1,ISPIN),WRK1,WRK2,WRK3,
     *                      FDER(1,ISPIN,KI4A),FDERFSC(1,1,ISPIN,KI5),
     *                      FDERFSC(1,1,ISPIN,KI6),FDER(1,ISPIN,KI6A),
     *                      WRK1,FON(ISPIN))
C
          IF (FON(3)) THEN
            !! WRK1 is generated in DFTB_POLDER1
            CALL DFTB_CPYTSQ(FDER(1,ISPIN,KI6A),WRK2,L1,L1)
            CALL DGEMM('T','N',NPART(ISPIN),NPART(ISPIN),L1,ONE,
     *                 FDERFSC(1,NDOCC(ISPIN)+1,ISPIN,KI5),L1,
     *                 WRK2(1,NDOCC(ISPIN)+1),L1,ZERO,
     *                 WRK3(NDOCC(ISPIN)+1,NDOCC(ISPIN)+1),L1)
            DO IDEG0 = 1, NDEG2(ISPIN)
              IF (ISPIN.EQ.1) IDEG = IDEG0
              IF (ISPIN.EQ.2) IDEG = IDEG0 + NDEG2(1)
              IO = LSTDEG2(1,IDEG)
              JO = LSTDEG2(2,IDEG)
              OCCI = OCC(IO,ISPIN)
              IF (NSPIN.EQ.1) OCCI = OCCI*HALF
              OCCJ = OCC(JO,ISPIN)
              IF (NSPIN.EQ.1) OCCJ = OCCJ*HALF
C             SCAL= -OCCI*(ONE-OCCI)*BETA
              EIGTMP(1) = EIG(IO,ISPIN)
              EIGTMP(2) = EIG(JO,ISPIN)
              OCCTMP(1) = OCCI
              OCCTMP(2) = OCCJ
              SCAL = DFTB_CLCDEG(1,0,0,BETA,EIGTMP,OCCTMP)
C             SCAL= -OCCJ*(ONE-OCCJ)*BETA
              VAL = VAL + WRK3(IO,JO)*WRK1(IO,JO)*SCAL
              VAL = VAL + WRK3(JO,IO)*WRK1(JO,IO)*SCAL
            END DO
            DO IO = NDOCC(ISPIN)+1, NOCC
              OCCI = OCC(IO,ISPIN)
              IF (NSPIN.EQ.1) OCCI = OCCI*HALF
              SCAL= -OCCI*(ONE-OCCI)*BETA
              VAL = VAL + WRK3(IO,IO)*WRK1(IO,IO)*SCAL
            END DO
          END IF
        END DO
      END DO
      VAL1 = VAL
C
      KI4 = KI1
      KI5 = KI2
      KI6 = KI3
      IF (FREQ.AND.HRAMAN) THEN
        IF (KI4.EQ.2) KI4 = 4
        IF (KI4.EQ.3) KI4 = 7
        IF (KI5.EQ.2) KI5 = 4
        IF (KI5.EQ.3) KI5 = 7
        IF (KI6.EQ.2) KI6 = 4
        IF (KI6.EQ.3) KI6 = 7
      ELSE IF (FREQ) THEN
        IF (KI4.EQ.3) KI4 = 5
        IF (KI4.EQ.2) KI4 = 3
        IF (KI5.EQ.3) KI5 = 5
        IF (KI5.EQ.2) KI5 = 3
        IF (KI6.EQ.3) KI6 = 5
        IF (KI6.EQ.2) KI6 = 3
      END IF
      VAL2 = ZERO
      IF (FON(3)) THEN
      !! loop (c)
      !! triply quasi-degenerated
        TMP = ZERO
        DO ISPIN = 1, NSPIN
          DO IDEG0 = 1, NDEG3(ISPIN)
            IF (ISPIN.EQ.1) IDEG = IDEG0
            IF (ISPIN.EQ.2) IDEG = IDEG0 + NDEG3(1)
            IO = LSTDEG3(1,IDEG)
            JO = LSTDEG3(2,IDEG)
            KO = LSTDEG3(3,IDEG)
C           OCCI = OCC(IO,ISPIN)
C           IF (NSPIN.EQ.1) OCCI = OCCI*HALF
C           SCAL = -TWO*OCCI*(ONE-OCCI)*(TWO*OCCI-ONE)*BETA*BETA
            EIGTMP(1) = EIG(IO,ISPIN)
            EIGTMP(2) = EIG(JO,ISPIN)
            EIGTMP(3) = EIG(KO,ISPIN)
            OCCTMP(1) = OCC(IO,ISPIN)*HALF
            OCCTMP(2) = OCC(JO,ISPIN)*HALF
            OCCTMP(3) = OCC(KO,ISPIN)*HALF
            SCAL = DFTB_CLCDEG(2,0,0,BETA,EIGTMP,OCCTMP)!*TWO
            CALL DFTB_CNVSQ(IO,JO,NSEQ1)
            CALL DFTB_CNVSQ(JO,KO,NSEQ2)
            CALL DFTB_CNVSQ(KO,IO,NSEQ3)
            EDERI = FDER(NSEQ1,ISPIN,KI4)
            EDERJ = FDER(NSEQ2,ISPIN,KI5)
            EDERK = FDER(NSEQ3,ISPIN,KI6)
            AAA = EDERI*EDERJ*EDERK
            TMP = TMP + SCAL*AAA
          END DO
        END DO
        VAL2=TMP
        val2=val2*two
      END IF
C
      !! term (2)+(3)-(4) eF
      VAL4 = ZERO
      TMP = ZERO
      IF (FON(3)) THEN
        DO ISPIN = 1, NSPIN
          DO IO = 1, L1 !NDOCC(ISPIN)+1, NOCC
            OCCI = OCC(IO,ISPIN)
            IF (NSPIN.EQ.1) OCCI = OCCI*HALF
            EI   = EIG(IO,ISPIN)
            DO MO = 1, L1
              EM   = EIG(MO,ISPIN)
              CALL DFTB_CNVSQ(MO,IO,NSEQ)
              IF (ABS(EI-EM).LE.THRES) THEN
                EIGTMP(1) = EIG(IO,ISPIN)
                EIGTMP(2) = EIG(MO,ISPIN)
                OCCTMP(1) = OCC(IO,ISPIN)*HALF
                OCCTMP(2) = OCC(MO,ISPIN)*HALF
                SCAL = -DFTB_CLCDEG(1,1,0,BETA,EIGTMP,OCCTMP)
              ELSE
                SCAL = ONE/(EI-EM)
                OCCM = OCC(MO,ISPIN)
                IF (NSPIN.EQ.1) OCCM = OCCM*HALF
                SCAL = SCAL*(OCCI*(ONE-OCCI)-OCCM*(ONE-OCCM))*BETA
              END IF
              FONDER1 = FEDER(KI1)
              FONDER2 = FEDER(KI2)
              FONDER3 = FEDER(KI3)
              AAA = FONDER1*FDER(NSEQ,1,KI5)*FDER(NSEQ,1,KI6)
     *            + FONDER2*FDER(NSEQ,1,KI6)*FDER(NSEQ,1,KI4)
     *            + FONDER3*FDER(NSEQ,1,KI4)*FDER(NSEQ,1,KI5)
              TMP = TMP + AAA*SCAL
C             WRK1(IO,1) = WRK1(IO,1) + AAA*SCAL*two
            END DO
          END DO
        END DO
        IF (NSPIN.EQ.1) TMP = TMP*TWO
        VAL4 = TMP
C
        TMP = ZERO
        ISPIN=1
        DO IO = NDOCC(ISPIN)+1, NOCC
          OCCI = OCC(IO,ISPIN)
          IF (NSPIN.EQ.1) OCCI = OCCI*HALF
          SCAL = -TWO*OCCI*(ONE-OCCI)*(TWO*OCCI-ONE)*BETA*BETA
          CALL DFTB_CNVSQ(IO,IO,NSEQ)
          EDERI = EDER(IO,1,KI1)
          EDERJ = EDER(IO,1,KI2)
          EDERK = EDER(IO,1,KI3)
          AAA = EDERI*FEDER(KI2)*FEDER(KI3)
     *        + FEDER(KI1)*EDERJ*FEDER(KI3)
     *        + FEDER(KI1)*FEDER(KI2)*EDERK
     *        - FEDER(KI1)*FEDER(KI2)*FEDER(KI3)
     *        + EDERI*EDERJ*EDERK !! possible to calculate above?
          TMP = TMP + SCAL*AAA
C           WRK1(IO,2) = WRK1(IO,2) + AAA*SCAL
C      wrk1(io,6) = wrk1(io,6) + ederi*feder(ki2)*feder(ki3)*scal
C      wrk1(io,7) = wrk1(io,7) + FEDER(KI1)*EDERJ*FEDER(KI3)*scal
C      wrk1(io,8) = wrk1(io,8) + FEDER(KI1)*FEDER(KI2)*EDERK*scal
        END DO
        VAL4 = VAL4 + TMP
C       TMP = ZERO
C       DO IO = 1, L1
C         TMP = TMP + WRK1(IO,1) + WRK1(IO,2)
C       END DO
C       VAL4 = TMP

      END IF
C
      VAL5 = ZERO
      IF (DFTB3) THEN
        IF (SRSCC) THEN
          DO NI = 1, NAT
            DO ISH0 = 1, MAXANG(ISPE(NI))
              ISH = INDSH(NI)+ISH0
              DO NJ = 1, NAT
                DO JSH0 = 1, MAXANG(ISPE(NJ))
                  JSH = INDSH(NJ)+JSH0
                  GAMI = GAMMA3(ISH,JSH)
C                 GAMJ = GAMMA3(JSH,ISH)
                  DO KK = 1, MAXANG(ISPE(NI))
                    KSH = INDSH(NI)+KK
                    VAL5 = VAL5
     *          +QDERF(ISH,1,KI4)*QDERF(JSH,1,KI5)*QDERF(KSH,1,KI6)*GAMI
     *          +QDERF(ISH,1,KI5)*QDERF(JSH,1,KI6)*QDERF(KSH,1,KI4)*GAMI
     *          +QDERF(ISH,1,KI6)*QDERF(JSH,1,KI4)*QDERF(KSH,1,KI5)*GAMI
                  END DO
                END DO
              END DO
            END DO
          END DO
          VAL5=VAL5*TWO
        ELSE
          DO NI = 1, NAT
            DO NJ = 1, NAT
              GAMI = GAMMA3(NI,NJ)
              GAMJ = GAMMA3(NJ,NI)
              VAL5 = VAL5 + QDERF(NI,1,KI5)     *QDERF(NJ,1,KI6)
     *              * (GAMI*QDERF(NI,1,KI4)+GAMJ*QDERF(NJ,1,KI4))
     *                    + QDERF(NI,1,KI6)     *QDERF(NJ,1,KI4)
     *              * (GAMI*QDERF(NI,1,KI5)+GAMJ*QDERF(NJ,1,KI5))
     *                    + QDERF(NI,1,KI4)     *QDERF(NJ,1,KI5)
     *              * (GAMI*QDERF(NI,1,KI6)+GAMJ*QDERF(NJ,1,KI6))
            END DO
          END DO
        END IF
      END IF
C
      VAL = VAL2*TWO
      VAL = VAL + VAL4
      VAL = VAL + VAL1*TWO
      IF (NSPIN.EQ.1) VAL = VAL + VAL1*TWO
      IF (DFTB3) VAL = VAL + VAL5/3.0D+00
C
      HYPER(IHYPER) = -VAL
C
      END DO
      END DO
      END DO
C
      IF (ISGDDI) THEN
        IF (GOPARR) CALL DSCAL(10,1.0D+00/DBLE(NPROC),HYPER,1)
        CALL GDDICOUNT( 1,LGROUP,MYJOB)
        CALL GDDI_SCOPE(DDI_WORLD)
        CALL DDI_GSUMF(2418,HYPER,10)
      END IF
C
      IF (MASWRK) THEN
        IHYPER = 0
        DO KI1 = 1, 3
        DO KI2 = 1, KI1
        DO KI3 = 1, KI2
          IHYPER = IHYPER + 1
          WRITE (IW,'(" HYPER-POLARIZABILITY FOR ",A1,"-",A1,
     *    "-",A1," = ", F25.10)') VCHARA(KI1),VCHARA(KI2),VCHARA(KI3),
     *    HYPER(IHYPER)
        END DO
        END DO
        END DO
        WRITE (IW,*)
      END IF
C
C
      END SUBROUTINE DFTB_CALC_HYPER
C
C-----------------------------------------------------------------------
C*MODULE DFTBHS    *DECK DFTB_CALC_POLDER
C>
C>    @brief   Calculate polarizability derivative
C>
C>    @details Calculate polarizability derivatives using Wigner's 2n+1
C>             rule. They are the mixed derivative of two electric
C>             fields and one nuclear coordinate. The third-order
C>             derivatives are computed with first-order response.
C>
C>    @author  Yoshio Nishimoto
C>    - Feb, 2018- Subroutine written
C>    @date    Mar, 2020 - Yoshio Nishimoto
C>             - Added LC-DFTB polarizability derivative
C>
C>           --- INPUT ---
C>    @param NI      Index of the atomic coordinate to be differentiated
C>    @param KI0     Index of the axis to be differentiated
C>    @param EIG     Eigenvalues of each orbital
C>    @param OCC     Occupation number of each orbital
C>    @param NAT     Number of atoms
C>    @param IND     Index of AO
C>    @param L1      Number of basis functions
C>    @param L2      L1*(L1+1)/2
C>    @param NSPIN   Number of spins (1=R-DFTB,2=U-DFTB)
C>    @param NDOCC   Number of doubly occupied orbitals
C>    @param NPART   Number of partially occupied orbitals
C>    @param C       Coordinates
C>    @param FON     Whether fractional occupation is used or not
C>    @param BETA    1/kTe
C>    @param EDERA   Derivatives of the orbital energy wrt a coordinate
C>    @param FEDERA  Derivatives of the Fermi energy wrt a coordinate
C>    @param EDER    Derivatives of the orbital energy wrt electric fiel
C>    @param FDER    Derivatives of the Fock matrix in MO wrt electric f
C>    @param FEDER   Derivatives of the Fermi energy wrt electric field
C>    @param QDERF   Derivatives of the Mulliken charge wrt electric fie
C>    @param QDER    Derivatives of the Mulliken charge wrt a coordinate
C>    @param DISTMAT Matrix of inter-atomic distances
C>    @param GAMMA   Matrix of DFTB2 gamma
C>    @param GAMDER  Derivative of DFTB2 gamma
C>    @param SDERAO  Overlap derivative in AO wrt a coordinate
C>    @param DDER    Density matrix derivative in AO wrt a coordinate
C>    @param WRK     Work space
C>    @param HDERIJ  Derivatives of the Fock matrix in MO wrt a coordina
C>    @param SDERIJ  Derivatives of the overlap matrix in MO wrt a coord
C>    @param FONDER  Derivatives of the occupation number wrt a coordina
C>    @param FDERFSC Square and scaled FDER
C>    @param FDERSC  Square and scaled HDERIJ
C>    @param SDERSC  Square and scaled SDERIJ
C>    @param NDEG2   Number of doubly quasi-degenerated pairs
C>    @param NDEG3   Number of triply quasi-degenerated pairs
C>    @param LSTDEG2 Index of doubly quasi-degenerated pairs
C>    @param LSTDEG3 Index of triply quasi-degenerated pairs
C>    @param DEGTHR  Degeneracy threshold
C>    @param WRK1    Work space
C>    @param WRK2    Work space
C>    @param WRK3    Work space
C>    @param DFTB3   Whether DFTB3 or not
C>    @param GAMMA3  Matrix of DFTB3 gamma
C>    @param GAMDER3 Derivative of DFTB3 gamma
C>    @param CHAMUL  Mulliken charge
C>    @parma MAXANG  Maximum angular momentum of each species
C>    @param ISPE    Index of species
C>    @param NDIMSH  If (SRSCC); NDIMSH=NSHELL; else NDIMSH=NAT
C>    @param INDSH   Index of shell for each atom
C>    @param SRSCC   Whether Shell-resolved SCC or not
C>    @param UHF     Whether U-DFTB or not
C>    @param SPNCST  Spin constants
C>    @param NSHELL  Number of shells
C>    @param NSPE    Number of species
C>    @param NGAMMA  Dimension of gamma matrices (NAT or NSHELL)
C>    @param CHAMULS Shell-resolved Mulliken charges
C>    @param LCDFTB  Flag for LC-DFTB
C>    @param S       Overlap matrix in AO
C>    @param GAMLC   Long-range gamma
C>    @param GAMDERLC Long-range gamma derivative
C>           --- OUTPUT ---
C>    @param ALPDER  Polarizability derivative
C>
C
      SUBROUTINE DFTB_CALC_POLDER(NI,KI0,OCC,EIG,NAT,IND,L1,L2,NSPIN,
     *                            NDOCC,NPART,C,
     *                            FON,BETA,EDERA,FEDERA,EDER,FDER,
     *                            FEDER,ALPDER,QDERF,QDER,DISTMAT,
     *                            GAMMA,GAMDER,SDERAO,DDER,WRK,HDERIJ,
     *                            SDERIJ,FONDER,FDERFSC,FDERSC,SDERSC,
     *                            NDEG2,NDEG3,LSTDEG2,LSTDEG3,DEGTHR,
     *                            WRK1,WRK2,WRK3,DFTB3,GAMMA3,GAMDER3,
     *                            CHAMUL,MAXANG,ISPE,NDIMSH,INDSH,SRSCC,
     *                            UHF,SPNCST,NSHELL,NSPE,NGAMMA,CHAMULS,
     *                            LCDFTB,S,GAMLC,GAMDERLC)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER (ZERO=0.0D+00,ONE=1.0D+00,TWO=2.0D+00,HALF=0.5D+00)
      PARAMETER (ONE_THIRD=ONE/3.0D+00)
C
      DIMENSION OCC(L1,NSPIN),EIG(L1,NSPIN),C(3,*),EDERA(L1,NSPIN),
     *          EDER(L1,NSPIN,3),FDER(L2,NSPIN,3),
     *          FEDER(4),ALPDER(6),QDERF(NDIMSH,NSPIN,3),
     *          QDER(NDIMSH),DISTMAT(*),GAMMA(*),GAMDER(*),
     *          SDERAO(*),DDER(L2,NSPIN,3),WRK(*),HDERIJ(L2,NSPIN),
     *          SDERIJ(L2,NSPIN),FONDER(L1,NSPIN,3),
     * FDERFSC(L1,L1,NSPIN,3),FDERSC(L1,L1,NSPIN),
     * SDERSC(L1,L1,NSPIN),
     * LSTDEG2(2,*),LSTDEG3(3,*),
     * WRK1(L1,L1),WRK2(L1,L1),WRK3(L1,L1),
     * GAMMA3(NGAMMA,NGAMMA),GAMDER3(NGAMMA,NGAMMA),CHAMUL(NAT),
     * MAXANG(*),ISPE(*),INDSH(*),SPNCST(6,*),CHAMULS(NDIMSH),
     * S(*),GAMLC(*),GAMDERLC(*)
      DIMENSION NDOCC(2),NPART(2),NDEG2(3),NDEG3(3)
      LOGICAL FON(3),DFTB3,SRSCC,UHF,LCDFTB
C
      DIMENSION SHIFT(NDIMSH)
C
      ISPIN = 1
      THRES = DEGTHR
C     write (*,*) "ndeg2,ndeg3 = ", ndeg2,ndeg3
      NOCC = NDOCC(ISPIN) + NPART(ISPIN)
      L3=L1*L1
C
      DO KI = 1, 6
      CALL DFTB_KI12(KI,KI1,KI2)
C
      VAL = ZERO
      ALPDER(KI) = ZERO
C
      IF (.false.) THEN
      val0=zero
      val1=zero
      DO ISPIN = 1, NSPIN
      DO IO = 1, L1
      OCCI = OCC(IO,ISPIN)*HALF
      EI   = EIG(IO,ISPIN)
      DO JO = 1, L1
      OCCJ = OCC(JO,ISPIN)*HALF
      EJ   = EIG(JO,ISPIN)
      CALL DFTB_CNVSQ(IO,JO,NSEQ1)
      HDERTMP = HDERIJ(NSEQ1,ISPIN)
      SDERTMP = SDERIJ(NSEQ1,ISPIN)
      DO KO = 1, L1
      OCCK = OCC(KO,ISPIN)*HALF
      EK   = EIG(KO,ISPIN)
C
      CALL DFTB_CNVSQ(JO,KO,NSEQ2)
      CALL DFTB_CNVSQ(KO,IO,NSEQ3)
      AAA  = FDER(NSEQ2,ISPIN,KI1)*FDER(NSEQ3,ISPIN,KI2)
      BBB  = SDERTMP*AAA
      AAA  = HDERTMP*AAA
      IF      (ABS(EI-EJ).LE.THRES.AND.ABS(EI-EK).LE.THRES) THEN
      SCALA= OCCI*(ONE-OCCI)*(TWO*OCCI-ONE)*BETA*BETA*HALF
      SCALB= SCALA*EI !! may not be correct
      ELSE IF (ABS(EI-EJ).LE.THRES) THEN
      SCALA= (OCCK-OCCI+OCCI*(ONE-OCCI)*(EK-EI)*BETA)/((EJ-EK)*(EK-EI))
      SCALB= (EK*(OCCK-OCCI)+EJ*(EK-EI)*OCCI*(ONE-OCCI)*BETA)/
     *       ((EJ-EK)*(EK-EI))
      ELSE IF (ABS(EJ-EK).LE.THRES) THEN
      SCALA= (OCCI-OCCJ+OCCJ*(ONE-OCCJ)*(EI-EJ)*BETA)/((EK-EI)*(EI-EJ))
      SCALB= (EI*(OCCI-OCCJ)+EK*(EI-EJ)*OCCJ*(ONE-OCCJ)*BETA)/
     *       ((EK-EI)*(EI-EJ))
      ELSE IF (ABS(EK-EI).LE.THRES) THEN
      SCALA= (OCCJ-OCCK+OCCK*(ONE-OCCK)*(EJ-EK)*BETA)/((EI-EJ)*(EJ-EK))
      SCALB= (EJ*(OCCJ-OCCK)+EI*(EJ-EK)*OCCK*(ONE-OCCK)*BETA)/
     *       ((EI-EJ)*(EJ-EK))
      ELSE
      SCALA= (OCCI*(EJ-EK)+OCCJ*(EK-EI)+OCCK*(EI-EJ))/
     *       ((EI-EJ)*(EJ-EK)*(EK-EI))
      SCALB= (EI*OCCI*(EJ-EK)+EJ*OCCJ*(EK-EI)+EK*OCCK*(EI-EJ))/
     *       ((EI-EJ)*(EJ-EK)*(EK-EI))
      END IF
C     write (*,'(3i3,2f20.10)') io,jo,ko, aaa*scala,-bbb*scalb
      VAL = VAL + AAA*SCALA-BBB*SCALB
C     end if
C       write (*,*) io,jo,ko
      END DO
      END DO
      END DO
      END DO
      VAL = -VAL*TWO
      IF (NSPIN.EQ.1) VAL = VAL*TWO
      write (*,*) "val = ", val
      ALPDER(KI) = ALPDER(KI) + VAL
      END IF
      val = zero
      val0 = zero
      val1 = zero
      tmp = alpder(ki)
C
C     term (a)-(c)
C
      !! term (a): afg
      DO ISPIN = 1, NSPIN
      NOCC = NDOCC(ISPIN) + NPART(ISPIN)
      CALL DFTB_POLDER1(1,NSPIN,ISPIN,L1,L2,L3,NDOCC(ISPIN),NOCC,NDEG2,
     *                  LSTDEG2,VAL,OCC(1,ISPIN),EIG(1,ISPIN),WRK1,WRK2,
     *                  WRK3,HDERIJ(1,ISPIN),FDERFSC(1,1,ISPIN,KI1),
     *                  FDERFSC(1,1,ISPIN,KI2),FDER(1,ISPIN,KI2),
     *                  WRK1,FON(ISPIN))
      !! term (b)+(c)+(d): afg
      CALL DFTB_POLDER1(2,NSPIN,ISPIN,L1,L2,L3,NDOCC(ISPIN),NOCC,NDEG2,
     *                  LSTDEG2,VAL,OCC(1,ISPIN),EIG(1,ISPIN),WRK1,WRK2,
     *                  WRK3,SDERIJ(1,ISPIN),FDERFSC(1,1,ISPIN,KI1),
     *                  FDERFSC(1,1,ISPIN,KI2),FDER(1,ISPIN,KI2),
     *                  WRK1,FON(ISPIN))
C
      !! term (a): fga
      CALL DFTB_POLDER1(3,NSPIN,ISPIN,L1,L2,L3,NDOCC(ISPIN),NOCC,NDEG2,
     *                  LSTDEG2,VAL,OCC(1,ISPIN),EIG(1,ISPIN),WRK1,WRK2,
     *                  WRK3,FDER(1,ISPIN,KI1),FDERFSC(1,1,ISPIN,KI2),
     *                  FDERSC(1,1,ISPIN),HDERIJ(1,ISPIN),
     *                  SDERSC(1,1,ISPIN),FON(ISPIN))
      !! term (b)+(c)+(d): fga
      CALL DFTB_POLDER1(4,NSPIN,ISPIN,L1,L2,L3,NDOCC(ISPIN),NOCC,NDEG2,
     *                  LSTDEG2,VAL,OCC(1,ISPIN),EIG(1,ISPIN),WRK1,WRK2,
     *                  WRK3,FDER(1,ISPIN,KI1),FDERFSC(1,1,ISPIN,KI2),
     *                  SDERSC(1,1,ISPIN),SDERIJ(1,ISPIN),
     *                  SDERSC(1,1,ISPIN),FON(ISPIN))
C
      !! term (a): gaf
      CALL DFTB_POLDER1(5,NSPIN,ISPIN,L1,L2,L3,NDOCC(ISPIN),NOCC,NDEG2,
     *                  LSTDEG2,VAL,OCC(1,ISPIN),EIG(1,ISPIN),WRK1,WRK2,
     *                  WRK3,FDER(1,ISPIN,KI2),FDERSC(1,1,ISPIN),
     *                  FDERFSC(1,1,ISPIN,KI1),FDER(1,ISPIN,KI1),
     *                  SDERSC(1,1,ISPIN),FON(ISPIN))
      !! term (b)+(c)+(d): gaf
      CALL DFTB_POLDER1(6,NSPIN,ISPIN,L1,L2,L3,NDOCC(ISPIN),NOCC,NDEG2,
     *                  LSTDEG2,VAL,OCC(1,ISPIN),EIG(1,ISPIN),WRK1,WRK2,
     *                  WRK3,FDER(1,ISPIN,KI2),SDERSC(1,1,ISPIN),
     *                  FDERFSC(1,1,ISPIN,KI1),FDER(1,ISPIN,KI1),
     *                  SDERSC(1,1,ISPIN),FON(ISPIN))
      END DO
C
      ALPDER(KI) = VAL
C
      IF (FON(3)) THEN
        DO ISPIN = 1, NSPIN
        !! (fi-fi)/(ei-ei) * F_{ii}^a * F_{mi}^f * F_{mi}^g / (ei-em)
        !! afg
        CALL DFTB_CPYTSQ(FDER(1,ISPIN,KI2),WRK1,L1,L1)
        DO IO = NDOCC(ISPIN)+1, NOCC
          OCCI = OCC(IO,ISPIN)
          IF (NSPIN.EQ.1) OCCI = OCCI*HALF
          EI = EIG(IO,ISPIN)
          SCAL= -OCCI*(ONE-OCCI)*BETA
          CALL DFTB_CNVSQ(IO,IO,NSEQ)
          VAL = VAL + (HDERIJ(NSEQ,ISPIN)-EI*SDERIJ(NSEQ,ISPIN))
     *        * SCAL*DDOT(L1,FDERFSC(1,IO,ISPIN,KI1),1,WRK1(1,IO),1)
        END DO
        DO IDEG0 = 1, NDEG2(ISPIN)
          IF (ISPIN.EQ.1) IDEG = IDEG0
          IF (ISPIN.EQ.2) IDEG = IDEG0 + NDEG2(1)
          IO = LSTDEG2(1,IDEG)
          JO = LSTDEG2(2,IDEG)
          CALL DFTB_CNVSQ(IO,JO,NSEQ)
          OCCI = OCC(IO,ISPIN)
          IF (NSPIN.EQ.1) OCCI = OCCI*HALF
          EI = EIG(IO,ISPIN)
          SCAL= -OCCI*(ONE-OCCI)*BETA
          VAL = VAL + (HDERIJ(NSEQ,ISPIN)-EI*SDERIJ(NSEQ,ISPIN))
     *        * SCAL*DDOT(L1,FDERFSC(1,IO,ISPIN,KI1),1,WRK1(1,JO),1)
          OCCJ = OCC(JO,ISPIN)
          IF (NSPIN.EQ.1) OCCJ = OCCJ*HALF
          EJ = EIG(JO,ISPIN)
          SCAL= -OCCJ*(ONE-OCCJ)*BETA
          VAL = VAL + (HDERIJ(NSEQ,ISPIN)-EJ*SDERIJ(NSEQ,ISPIN))
     *        * SCAL*DDOT(L1,FDERFSC(1,JO,ISPIN,KI1),1,WRK1(1,IO),1)
        END DO
C
        !! fga
        CALL DFTB_CPYTSQ(HDERIJ(1,ISPIN),WRK1,L1,L1)
        CALL DFTB_CPYTSQ(SDERIJ(1,ISPIN),WRK2,L1,L1)
        DO IO = NDOCC(ISPIN)+1, NOCC
          OCCI = OCC(IO,ISPIN)
          IF (NSPIN.EQ.1) OCCI = OCCI*HALF
          EI = EIG(IO,ISPIN)
          SCAL= -OCCI*(ONE-OCCI)*BETA
          CALL DFTB_CNVSQ(IO,IO,NSEQ)
          VAL = VAL + FDER(NSEQ,ISPIN,KI1)*SCAL*
     *         (DDOT(L1,FDERFSC(1,IO,ISPIN,KI2),1,WRK1(1,IO),1)
     *      -EI*DDOT(L1,FDERFSC(1,IO,ISPIN,KI2),1,WRK2(1,IO),1))
        END DO
        DO IDEG0 = 1, NDEG2(ISPIN)
          IF (ISPIN.EQ.1) IDEG = IDEG0
          IF (ISPIN.EQ.2) IDEG = IDEG0 + NDEG2(1)
          IO = LSTDEG2(1,IDEG)
          JO = LSTDEG2(2,IDEG)
          CALL DFTB_CNVSQ(IO,JO,NSEQ)
          OCCI = OCC(IO,ISPIN)
          IF (NSPIN.EQ.1) OCCI = OCCI*HALF
          EI = EIG(IO,ISPIN)
          SCAL= -OCCI*(ONE-OCCI)*BETA
          VAL = VAL + FDER(NSEQ,ISPIN,KI1)*SCAL*
     *         (DDOT(L1,FDERFSC(1,IO,ISPIN,KI2),1,WRK1(1,JO),1)
     *      -EI*DDOT(L1,FDERFSC(1,IO,ISPIN,KI2),1,WRK2(1,JO),1))
          OCCJ = OCC(JO,ISPIN)
          IF (NSPIN.EQ.1) OCCJ = OCCJ*HALF
          EJ = EIG(JO,ISPIN)
          SCAL= -OCCJ*(ONE-OCCJ)*BETA
          VAL = VAL + FDER(NSEQ,ISPIN,KI1)*SCAL*
     *         (DDOT(L1,FDERFSC(1,JO,ISPIN,KI2),1,WRK1(1,IO),1)
     *      -EJ*DDOT(L1,FDERFSC(1,JO,ISPIN,KI2),1,WRK2(1,IO),1))
        END DO
C
        !! gaf
        CALL DFTB_CPYTSQ(FDER(1,ISPIN,KI1),WRK1,L1,L1)
        DO IO = NDOCC(ISPIN)+1, NOCC
          OCCI = OCC(IO,ISPIN)
          IF (NSPIN.EQ.1) OCCI = OCCI*HALF
          EI = EIG(IO,ISPIN)
          SCAL= -OCCI*(ONE-OCCI)*BETA
          CALL DFTB_CNVSQ(IO,IO,NSEQ)
          VAL = VAL + FDER(NSEQ,ISPIN,KI2)*SCAL*
     *         (DDOT(L1,FDERSC(1,IO,ISPIN),1,WRK1(1,IO),1)
     *      -EI*DDOT(L1,SDERSC(1,IO,ISPIN),1,WRK1(1,IO),1))
        END DO
        DO IDEG0 = 1, NDEG2(ISPIN)
          IF (ISPIN.EQ.1) IDEG = IDEG0
          IF (ISPIN.EQ.2) IDEG = IDEG0 + NDEG2(1)
          IO = LSTDEG2(1,IDEG)
          JO = LSTDEG2(2,IDEG)
          CALL DFTB_CNVSQ(IO,JO,NSEQ)
          OCCI = OCC(IO,ISPIN)
          IF (NSPIN.EQ.1) OCCI = OCCI*HALF
          SCAL= -OCCI*(ONE-OCCI)*BETA
          VAL = VAL + FDER(NSEQ,ISPIN,KI2)*SCAL*
     *         (DDOT(L1,FDERSC(1,IO,ISPIN),1,WRK1(1,JO),1)
     *      -EI*DDOT(L1,SDERSC(1,IO,ISPIN),1,WRK1(1,JO),1))
          OCCJ = OCC(JO,ISPIN)
          IF (NSPIN.EQ.1) OCCJ = OCCJ*HALF
          SCAL= -OCCJ*(ONE-OCCJ)*BETA
          VAL = VAL + FDER(NSEQ,ISPIN,KI2)*SCAL*
     *         (DDOT(L1,FDERSC(1,JO,ISPIN),1,WRK1(1,IO),1)
     *      -EJ*DDOT(L1,SDERSC(1,JO,ISPIN),1,WRK1(1,IO),1))
        END DO
        END DO
      END IF
C
      IF (FON(3)) THEN
      !! loop (c)
      !! triply degenerated
        TMP = ZERO
C     val0=val
        tmp0 = zero
        tmp1 = zero
        ISPIN = 1
        DO IDEG0 = 1, NDEG3(ISPIN)
          IF (ISPIN.EQ.1) IDEG = IDEG0
          IF (ISPIN.EQ.2) IDEG = IDEG0 + NDEG3(1)
          IO = LSTDEG3(1,IDEG)
          JO = LSTDEG3(2,IDEG)
          KO = LSTDEG3(3,IDEG)
          OCCI = OCC(IO,ISPIN)
          IF (NSPIN.EQ.1) OCCI = OCCI*HALF
          SCAL = -TWO*OCCI*(ONE-OCCI)*(TWO*OCCI-ONE)*BETA*BETA
          CALL DFTB_CNVSQ(IO,JO,NSEQ1)
          CALL DFTB_CNVSQ(JO,KO,NSEQ2)
          CALL DFTB_CNVSQ(KO,IO,NSEQ3)
          EDERI = HDERIJ(NSEQ1,ISPIN)
          EDERJ = FDER(NSEQ2,ISPIN,KI1)
          EDERK = FDER(NSEQ3,ISPIN,KI2)
          AAA = EDERJ*EDERK
          TMP = TMP + SCAL*AAA*EDERI
C
C         IF (IO.NE.JO.OR.IO.NE.KO) SCAL=SCAL-OCCI*(ONE-OCCI)*BETA*TWO
          EI = EIG(IO,ISPIN)
          EDERI = SDERIJ(NSEQ1,ISPIN)
          TMP = TMP - SCAL*AAA*EDERI*EI
        END DO
        DO IO = NDOCC(ISPIN)+1, NOCC
          OCCI = OCC(IO,ISPIN)*HALF
          SCAL = -TWO*OCCI*(ONE-OCCI)*(TWO*OCCI-ONE)*BETA*BETA
          CALL DFTB_CNVSQ(IO,IO,NSEQ)
          EDERI = HDERIJ(NSEQ,ISPIN)
          EDERJ = FDER(NSEQ,ISPIN,KI1)
          EDERK = FDER(NSEQ,ISPIN,KI2)
          AAA = EDERJ*EDERK
          TMP = TMP + SCAL*AAA*ederi
C
          EI = EIG(IO,ISPIN)
          EDERI = SDERIJ(NSEQ,ISPIN)
          TMP = TMP - SCAL*AAA*EDERI*EI
        END DO
        TMP = TMP*HALF
        IF (NSPIN.EQ.1) TMP = TMP*HALF !! ??? not checked
        VAL=VAL+TMP
        VAL2=TMP
      val0 = val0 + tmp0*half*half
      val1 = val1 + tmp1*half*half
C     write (*,'(a,f20.10)') "val0 = ", val0*two*two
C     write (*,'(a,f20.10)') "val1 = ", val1*two*two
C     write (*,'(a,f20.10)') "tmp0 = ", tmp0
C     write (*,'(a,f20.10)') "tmp1 = ", tmp1
      END IF
      VAL=VAL*TWO
      IF (NSPIN.EQ.1) VAL=VAL*TWO
C     write (*,'(a,f20.10)') "val  = ", val
      ALPDER(KI) = VAL
C
      IF (FON(3)) THEN
        VAL = ZERO
        TMP = ZERO
        !! Fermi energy derivative term (1)
        ISPIN = 1
        DO IO = NDOCC(ISPIN)+1, NOCC
          OCCI = OCC(IO,ISPIN)*HALF
          EI   = EIG(IO,ISPIN)
          IF (IO.LE.NDOCC(ISPIN)) NNN = NDOCC(ISPIN)+1
          IF (IO.GT.NDOCC(ISPIN)) NNN = 1
          DO MO = NNN, L1
            EM   = EIG(MO,ISPIN)
            CALL DFTB_CNVSQ(MO,IO,NSEQ)
            IF (ABS(EI-EM).LE.THRES) THEN
              SCAL = OCCI*(ONE-OCCI)*(TWO*OCCI-ONE)*BETA*BETA*HALF
            ELSE
              SCAL = OCCI*(ONE-OCCI)*BETA/(EI-EM)
            END IF
            FONDER1 = FEDERA
            FONDER2 = FEDER(KI1)
            FONDER3 = FEDER(KI2)
            AAA = FONDER1*FDER(NSEQ,1,KI1)*FDER(NSEQ,1,KI2)
     *          + FONDER2*FDER(NSEQ,1,KI2)*HDERIJ(NSEQ,1)
     *          + FONDER3*HDERIJ(NSEQ,1)  *FDER(NSEQ,1,KI1)
            BBB = FONDER2*FDER(NSEQ,1,KI2)*SDERIJ(NSEQ,1)
     *          + FONDER3*SDERIJ(NSEQ,1)  *FDER(NSEQ,1,KI1)
            TMP = TMP + SCAL*(AAA-EI*BBB)
          END DO
        END DO
        TMP = TMP*TWO
        IF (NSPIN.EQ.1) TMP = TMP*TWO
        VAL = VAL + TMP
C
        TMP = ZERO
        !! Fermi energy derivative term (2)+(3)
        DO IO = NDOCC(ISPIN)+1, NOCC
          OCCI = OCC(IO,1)*HALF
          SCAL = -OCCI*(ONE-OCCI)*(TWO*OCCI-ONE)*BETA*BETA
          CALL DFTB_CNVSQ(IO,IO,NSEQ)
          EDERI = EDERA(IO,1)
          EDERJ = EDER(IO,1,KI1)
          EDERK = EDER(IO,1,KI2)
          AAA = EDERI *FEDER(KI1)*FEDER(KI2)
     *        + FEDERA*EDERJ     *FEDER(KI2)
     *        + FEDERA*FEDER(KI1)*EDERK
     *        - FEDERA*FEDER(KI1)*FEDER(KI2)
          TMP = TMP + SCAL*AAA
        END DO
        VAL = VAL + TMP*TWO
C
        !! unique term
        TMP = ZERO
        DO IO = NDOCC(ISPIN)+1, NOCC
          CALL DFTB_CNVSQ(IO,IO,NSEQ)
          AAA = FONDER(IO,1,KI1)*FDER(NSEQ,1,KI2)
     *        + FONDER(IO,1,KI2)*FDER(NSEQ,1,KI1)
          TMP = TMP + AAA*SDERIJ(NSEQ,1)
        END DO
        VAL = VAL - TMP
C
        ALPDER(KI) = ALPDER(KI) + VAL
      END IF
C
C     REMAINING TERMS ARE ZERO FOR HYPER-POLARIZABILITY
C
      !! E^{af}*C^g
      IF (KI1.EQ.KI0) THEN
        IF (SRSCC) THEN
          DO ISH0 = 1, MAXANG(ISPE(NI))
            ISH = INDSH(NI)+ISH0
            ALPDER(KI) = ALPDER(KI) + QDERF(ISH,1,KI2)
          END DO
        ELSE
          ALPDER(KI) = ALPDER(KI) + QDERF(NI,1,KI2)
        END IF
      END IF
      !! E^{ag}*C^f
      IF (KI2.EQ.KI0) THEN
        IF (SRSCC) THEN
          DO ISH0 = 1, MAXANG(ISPE(NI))
            ISH = INDSH(NI)+ISH0
            ALPDER(KI) = ALPDER(KI) + QDERF(ISH,1,KI1)
          END DO
        ELSE
          ALPDER(KI) = ALPDER(KI) + QDERF(NI,1,KI1)
        END IF
      END IF
C
      !! E^a*C^f*C^g (gamma derivative)
      IF (SRSCC) THEN
        DO ISH0 = 1, MAXANG(ISPE(NI))
          ISH = INDSH(NI)+ISH0
          QFA = QDERF(ISH,1,KI1)
          QGA = QDERF(ISH,1,KI2)
          DO NJ = 1, NAT
            IF (NI.EQ.NJ) CYCLE
            CALL DFTB_CNVSQ(NI,NJ,NSEQ)
            DIST = DISTMAT(NSEQ)
            DISTI = ONE/DIST
            VAL = ZERO
            DO JSH0 = 1, MAXANG(ISPE(NJ))
              JSH = INDSH(NJ)+JSH0
              CALL DFTB_CNVSQ(ISH,JSH,NSEQ)
              GAMMAV = GAMDER(NSEQ)
              QFB = QDERF(JSH,1,KI1)
              QGB = QDERF(JSH,1,KI2)
              VAL = VAL + (QFA*QGB + QGA*QFB)*GAMMAV
              IF (DFTB3) THEN
                DQI = CHAMULS(ISH)
                DQJ = CHAMULS(JSH)
                GAMI = GAMDER3(ISH,JSH) !! G_ac
                GAMJ = GAMDER3(JSH,ISH) !! G_ca
                DO KSH0 = 1, MAXANG(ISPE(NI))
                  KSH = INDSH(NI)+KSH0
                  QFC = QDERF(KSH,1,KI1)
                  QGC = QDERF(KSH,1,KI2)
                  DQK = CHAMULS(KSH)
                  VAL2 = GAMI*(
     *                QFA*QGB*DQK+QGA*QFB*DQK
     *              + QFA*DQJ*QGC+QGA*DQJ*QFC
     *              + DQI*QFB*QGC+DQI*QGB*QFC)
                  VAL = VAL + VAL2*ONE_THIRD
                END DO
                DO KSH0 = 1, MAXANG(ISPE(NJ))
                  KSH = INDSH(NJ)+KSH0
                  QFC = QDERF(KSH,1,KI1)
                  QGC = QDERF(KSH,1,KI2)
                  DQK = CHAMULS(KSH)
                  VAL2 = GAMJ*(
     *                QFA*QGB*DQK+QGA*QFB*DQK
     *              + QFA*DQJ*QGC+QGA*DQJ*QFC
     *              + DQI*QFB*QGC+DQI*QGB*QFC)
                  VAL = VAL + VAL2*ONE_THIRD
                END DO
              END IF
            END DO
            VAL = VAL*DISTI*(C(KI0,NI)-C(KI0,NJ))
            ALPDER(KI) = ALPDER(KI) + VAL
          END DO
        END DO
      ELSE
        QFA = QDERF(NI,1,KI1)
        QGA = QDERF(NI,1,KI2)
        DO NJ = 1, NAT
          IF (NI.EQ.NJ) CYCLE
          CALL DFTB_CNVSQ(NI,NJ,NSEQ)
          DIST = DISTMAT(NSEQ)
          DISTI = ONE/DIST
          GAMMAV = GAMDER(NSEQ)
          QFB = QDERF(NJ,1,KI1)
          QGB = QDERF(NJ,1,KI2)
          VAL = (QFA*QGB + QGA*QFB)*GAMMAV
          IF (DFTB3) THEN
            GAMI = GAMDER3(NI,NJ) !! G_ac
            GAMJ = GAMDER3(NJ,NI) !! G_ca
            DQI = CHAMUL(NI)
            DQJ = CHAMUL(NJ)
            VAL2 = TWO*GAMI*(QFA*QGA*DQJ + QFA*DQI*QGB + QGA*DQI*QFB)
     *           + TWO*GAMJ*(QFB*QGB*DQI + QFA*DQJ*QGB + QGA*DQJ*QFB)
            VAL = VAL + VAL2*ONE_THIRD
          END IF
          VAL = VAL*DISTI*(C(KI0,NI)-C(KI0,NJ))
          ALPDER(KI) = ALPDER(KI) + VAL
        END DO
      END IF
C
      IF (UHF) THEN
        CALL VCLR(WRK2,1,NSHELL)
        CALL DFTB_SPIN_SHIFT(WRK2,QDERF(1,2,KI1),SPNCST,NAT,
     *                       NSHELL,NSPE,ISPE,MAXANG)
      END IF
      DO ISPIN = 1, NSPIN
        !! C^g*E^a*C^f (Smn derivative)
        CALL VCLR(SHIFT,1,NDIMSH)
        CALL DFTB_2ND_SHIFT(SHIFT,QDERF(1,1,KI1),GAMMA,NGAMMA)
        CALL DFTB_3RDAR(NAT,NGAMMA,NSHELL,NSPE,ISPE,MAXANG,SHIFT,
     *                  QDERF(1,1,KI1),CHAMUL,CHAMULS,GAMMA3,WRK,
     *                  DFTB3,SRSCC)
        IF (SRSCC) THEN
          DO NJ = 1, NAT
            DO ISH0 = 1, MAXANG(ISPE(NJ))
              ISH = INDSH(NJ)+ISH0
              SHIFT(ISH) = SHIFT(ISH) + C(KI1,NJ)
            END DO
          END DO
        ELSE
          CALL DAXPY(NAT,ONE,C(KI1,1),3,SHIFT,1) !! E^{af}*C^g
        END IF
        IF (UHF) THEN
          IF (.NOT.SRSCC) THEN
            !! ATOM-RESOLVED --> SHELL-RESOLVED
            CALL DCOPY(NAT,SHIFT,1,WRK,1)
            CALL DFTB_SHIFT_ATOSH(WRK,SHIFT,NAT,NSHELL,NSPE,ISPE,MAXANG)
          END IF
          IF (ISPIN.EQ.1) CALL DAXPY(NSHELL, ONE,WRK2,1,SHIFT,1)
          IF (ISPIN.EQ.2) CALL DAXPY(NSHELL,-ONE,WRK2,1,SHIFT,1)
        END IF
        IF (SRSCC.OR.UHF) THEN
          CALL SHIFT_TRIANGLE_SHELL(SHIFT,SDERAO,WRK,
     *      NAT,NSPE,NSHELL,L2,ISPE,MAXANG)
        ELSE
          CALL SHIFT_TRIANGLE(NAT,L2,IND,SHIFT,SDERAO,WRK)
        END IF
        ALPDER(KI) = ALPDER(KI) + TRACEP(DDER(1,ISPIN,KI2),WRK,L1)
      END DO
C
      IF (UHF) THEN
        CALL VCLR(WRK2,1,NSHELL)
        CALL DFTB_SPIN_SHIFT(WRK2,QDERF(1,2,KI2),SPNCST,NAT,
     *                       NSHELL,NSPE,ISPE,MAXANG)
      END IF
      DO ISPIN = 1, NSPIN
        !! C^f*E^a*C^g (Smn derivative)
        CALL VCLR(SHIFT,1,NDIMSH)
        CALL DFTB_2ND_SHIFT(SHIFT,QDERF(1,1,KI2),GAMMA,NGAMMA)
        CALL DFTB_3RDAR(NAT,NGAMMA,NSHELL,NSPE,ISPE,MAXANG,SHIFT,
     *                  QDERF(1,1,KI2),CHAMUL,CHAMULS,GAMMA3,WRK,
     *                  DFTB3,SRSCC)
        IF (SRSCC) THEN
          DO NJ = 1, NAT
            DO ISH0 = 1, MAXANG(ISPE(NJ))
              ISH = INDSH(NJ)+ISH0
              SHIFT(ISH) = SHIFT(ISH) + C(KI2,NJ)
            END DO
          END DO
        ELSE
          CALL DAXPY(NAT,ONE,C(KI2,1),3,SHIFT,1) !! E^{ag}*C^f
        END IF
        IF (UHF) THEN
          IF (.NOT.SRSCC) THEN
            CALL DCOPY(NAT,SHIFT,1,WRK,1)
            CALL DFTB_SHIFT_ATOSH(WRK,SHIFT,NAT,NSHELL,NSPE,ISPE,MAXANG)
          END IF
          IF (ISPIN.EQ.1) CALL DAXPY(NSHELL, ONE,WRK2,1,SHIFT,1)
          IF (ISPIN.EQ.2) CALL DAXPY(NSHELL,-ONE,WRK2,1,SHIFT,1)
        END IF
        IF (SRSCC.OR.UHF) THEN
          CALL SHIFT_TRIANGLE_SHELL(SHIFT,SDERAO,WRK,
     *      NAT,NSPE,NSHELL,L2,ISPE,MAXANG)
        ELSE
          CALL SHIFT_TRIANGLE(NAT,L2,IND,SHIFT,SDERAO,WRK)
        END IF
        ALPDER(KI) = ALPDER(KI) + TRACEP(DDER(1,ISPIN,KI1),WRK,L1)
      END DO
C
      IF (LCDFTB) THEN
        IF (NSPIN.EQ.1) SCAL = 2.0D+00
        IF (NSPIN.EQ.2) SCAL = 4.0D+00
        DO ISPIN = 1, NSPIN
          CALL CPYTSQ(DDER(1,ISPIN,KI1),WRK1,L1,1)
          CALL CPYTSQ(DDER(1,ISPIN,KI2),WRK2,L1,1)
          !! Smn derivative
          CALL CPYSQT(GAMLC,WRK,NAT,1)
          CALL DFTB_LCGRAD1(1,NAT,L1,ISPE,MAXANG,IND,WRK1,WRK2,
     *                      S,WRK,WRK3)
          ALPDER(KI) = ALPDER(KI) + TRACEP(SDERAO,WRK3,L1)*SCAL
          !! gamma derivative
          !! DFTB_LCGRAD2 computes for all nuclear coordinates.
          !! One can write a subroutine specific for the polarizability
          !! derivative, but I think not much benefit is obtained, because
          !! the matrix multiplication is more computationally demanding
          !! than calculating contributions for all nuclear coordinates.
          CALL CPYTSQ(S,WRK3,L1,1)
          CALL VCLR(WRK,1,3*NAT)
          CALL DFTB_LCGRAD2(1,SCAL,NAT,L1,ISPE,MAXANG,IND,WRK1,
     *                      WRK2,WRK3,WRK3,DISTMAT,
     *                      GAMDERLC,C,WRK,DUM)
          ALPDER(KI) = ALPDER(KI) + WRK(KI0+3*(NI-1))
        END DO
      END IF
C
      !! q^a*q^f*q^g for DFTB3
      IF (DFTB3) THEN
        TMP=ZERO
        IF (SRSCC) THEN
          DO NJ = 1, NAT
            DO JSH0 = 1, MAXANG(ISPE(NJ))
              JSH = INDSH(NJ)+JSH0
              DO NK = 1, NAT
                DO KSH0 = 1, MAXANG(ISPE(NK))
                  KSH = INDSH(NK)+KSH0
                  GAMJ = GAMMA3(JSH,KSH)
                  GAMK = GAMMA3(KSH,JSH)
                  DO LL = 1, MAXANG(ISPE(NJ))
                    LSH = INDSH(NJ)+LL
                    TMP = TMP
     *          +QDERF(JSH,1,KI1)*QDERF(KSH,1,KI2)*QDER (LSH)      *GAMJ
     *          +QDERF(JSH,1,KI2)*QDER (KSH)      *QDERF(LSH,1,KI1)*GAMJ
     *          +QDER (JSH)      *QDERF(KSH,1,KI1)*QDERF(LSH,1,KI2)*GAMJ
                  END DO
                  DO LL = 1, MAXANG(ISPE(NK))
                    LSH = INDSH(NK)+LL
                    TMP = TMP
     *          +QDERF(JSH,1,KI1)*QDERF(KSH,1,KI2)*QDER (LSH)      *GAMK
     *          +QDERF(JSH,1,KI2)*QDER (KSH)      *QDERF(LSH,1,KI1)*GAMK
     *          +QDER (JSH)      *QDERF(KSH,1,KI1)*QDERF(LSH,1,KI2)*GAMK
                  END DO
                END DO
              END DO
            END DO
          END DO
        ELSE
          DO NJ = 1, NAT
            DO NK = 1, NAT
              GAMJ = GAMMA3(NJ,NK)
              GAMK = GAMMA3(NK,NJ)
              TMP = TMP + QDERF(NJ,1,KI1)     *QDERF(NK,1,KI2)
     *            * (GAMJ*QDER (NJ)      +GAMK*QDER (NK))
     *                  + QDERF(NJ,1,KI2)     *QDER (NK)
     *            * (GAMJ*QDERF(NJ,1,KI1)+GAMK*QDERF(NK,1,KI1))
     *                  + QDER (NJ)           *QDERF(NK,1,KI1)
     *            * (GAMJ*QDERF(NJ,1,KI2)+GAMK*QDERF(NK,1,KI2))
            END DO
          END DO
        END IF
        ALPDER(KI) = ALPDER(KI) + TMP*ONE_THIRD
      END IF
C
      END DO
      TMP = ALPDER(3)
      ALPDER(3) = ALPDER(4)
      ALPDER(4) = TMP
C
      END SUBROUTINE DFTB_CALC_POLDER
C
C-----------------------------------------------------------------------
C*MODULE DFTBHS    *DECK DFTB_POLDER1
C>
C>    @brief   Calculate polarizability derivative
C>
C>    @details Calculate the first term in Eq. (30) and the first and
C>             second terms in Eq. (34) in the third-order derivative
C>             paper, although this subroutine does not compute all
C>             the contributions.
C>
C>    @author  Yoshio Nishimoto
C>    - Feb, 2018- Subroutine written
C>
C>           --- INPUT ---
C>    @param MODE    This specifies permutation of a, F, and G
C>    @param NSPIN   Number of spins (1=R-DFTB,2=U-DFTB)
C>    @param ISPIN   Current index of spin
C>    @param L1      Number of basis functions
C>    @param L2      L1*(L1+1)/2
C>    @param L3      L1*L1
C>    @param NDOCC   Number of doubly occupied orbitals
C>    @param NOCC    Number of doubly + partially occupied orbitals
C>    @param NDEG2   Number of doubly quasi-degenerated pairs
C>    @param LSTDEG2 Index of doubly quasi-degenerated pairs
C>    @param OCC     Occupation number of each orbital
C>    @param EIG     Eigenvalues of each orbital
C>    @param WRK1    Work space
C>    @param WRK2    Work space
C>    @param WRK3    Work space
C>    @param FA      Derivatives of Fock matrix in MO (1)
C>    @param FB      Derivatives of symmetric and scaled Fock matrix (2)
C>    @param FC      Derivatives of symmetric and scaled Fock matrix (3)
C>    @param FDER    Derivatives of Fock matrix (3)
C>    @param FBC2    FB*FC
C>    @param FON     Whether fractional occupation is used or not
C>           --- OUTPUT ---
C>    @param VAL     Computed result
C>
C
      SUBROUTINE DFTB_POLDER1(MODE,NSPIN,ISPIN,L1,L2,L3,NDOCC,NOCC,
     *                        NDEG2,LSTDEG2,VAL,OCC,EIG,WRK1,WRK2,WRK3,
     *                        FA,FB,FC,FDER,FBC2,FON)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (ZERO=0.0D+00,HALF=0.5D+00,ONE=1.0D+00)
C
      DIMENSION OCC(L1),EIG(L1),WRK1(L1,L1),WRK2(L1,L1),
     *          WRK3(L1,L1),FA(L2),FB(L1,L1),FC(L1,L1),FDER(L2),
     *          FBC2(L1,L1)
      DIMENSION NDEG2(3),LSTDEG2(2,*)
      LOGICAL   FON
C
      TMP = VAL
      VAL = ZERO
      NPART = NOCC-NDOCC
C
      !! WRK1 should not be modified between MODE=3-4 and MODE=5-6.
      IF (MODE.LE.3.OR.MODE.EQ.5) CALL DFTB_CPYTSQ(FA,WRK1,L1,L1)
C
      IF (MODE.EQ.1.OR.MODE.EQ.2) THEN
        !! F_{mi}^b * F_{ni}^c = (FF)_{mn}^bc
        CALL DCOPY(L1*NOCC,FB,1,WRK2,1)
        DO IO = 1, NOCC
          OCCI = OCC(IO)
          IF (NSPIN.EQ.1) OCCI = OCCI*HALF
          IF (MODE.EQ.2) OCCI = OCCI*EIG(IO)
          CALL DSCAL(L1,OCCI,WRK2(1,IO),1)
        END DO
        CALL DGEMM('N','T',L1,L1,NOCC,ONE,WRK2,L1,FC,L1,ZERO,WRK3,L1)
        IF (FON.AND.NDOCC.GE.1)
     *      CALL DGEMM('N','T',NDOCC,NDOCC,NPART,ONE,
     *                 WRK2(1,NDOCC+1),L1,FC(1,NDOCC+1),L1,
     *                 ZERO,WRK3,L1)
      
C
        !! F_{mn}^a * (FF)_{mn}^bc
        IF (FON) THEN
          VAL = VAL + DDOT(L3,WRK1,1,WRK3,1)
        ELSE
          DO NO = 1, NDOCC
            VAL = VAL
     *          + DDOT(L1-NDOCC,WRK1(NDOCC+1,NO),1,WRK3(NDOCC+1,NO),1)
          END DO
          VAL = VAL + DDOT(L1*(L1-NDOCC),
     *                     WRK1(1,NDOCC+1),1,WRK3(1,NDOCC+1),1)
        END IF
C       if (mode.eq.2) write (*,*) "val1 = ", val
      ELSE IF (MODE.EQ.3) THEN
        !! F_{mn}^f * F_{mi}^g = (FF)_{ni}^fg
        CALL DGEMM('T','N',L1,NOCC,L1,ONE,WRK1,L1,FB,L1,ZERO,WRK3,L1)
        !! (FF)_{ni}^fg * F_{ni}^a - e_i * (FF)_{ni}^fg * S_{ni}^a
        DO IO = 1, NOCC
          OCCI = OCC(IO)
          IF (NSPIN.EQ.1) OCCI = OCCI*HALF
          EI = EIG(IO)
          VAL = VAL +    OCCI*DDOT(L1,WRK3(1,IO),1,FC  (1,IO),1)
     *              - EI*OCCI*DDOT(L1,WRK3(1,IO),1,FBC2(1,IO),1)
        END DO
C
        !! remove occ-occ-occ part
        CALL DGEMM('T','N',NDOCC,NDOCC,NDOCC,ONE,WRK1,L1,FB,L1,
     *             ZERO,WRK3,L1)
        DO IO = 1, NDOCC
          OCCI = OCC(IO)
          IF (NSPIN.EQ.1) OCCI = OCCI*HALF
          EI = EIG(IO)
          VAL = VAL -    OCCI*DDOT(NDOCC,WRK3(1,IO),1,FC  (1,IO),1)
     *              + EI*OCCI*DDOT(NDOCC,WRK3(1,IO),1,FBC2(1,IO),1)
        END DO
      ELSE IF (MODE.EQ.4) THEN
        CONTINUE !! ALREADY DONE IN MODE=3
      ELSE IF (MODE.EQ.5) THEN
        !! F_{mn}^g (sym) * F_{ni}^f = (FF)_{mi}^gf
        CALL DGEMM('T','N',L1,NOCC,L1,ONE,WRK1,L1,FC,L1,ZERO,WRK3,L1)
        !! (FF)_{mi}^gf * F_{mi}^a - e_i * (FF)_{mi}^gf * S_{mi}^a
        DO IO = 1, NOCC
          OCCI = OCC(IO)
          IF (NSPIN.EQ.1) OCCI = OCCI*HALF
          EI = EIG(IO)
          VAL = VAL +    OCCI*DDOT(L1,WRK3(1,IO),1,FB  (1,IO),1)
     *              - EI*OCCI*DDOT(L1,WRK3(1,IO),1,FBC2(1,IO),1)
        END DO
C
        !! remove occ-occ-occ part
        CALL DGEMM('T','N',NDOCC,NDOCC,NDOCC,ONE,WRK1,L1,FC,L1,
     *             ZERO,WRK3,L1)
        DO IO = 1, NDOCC
          OCCI = OCC(IO)
          IF (NSPIN.EQ.1) OCCI = OCCI*HALF
          EI = EIG(IO)
          VAL = VAL -    OCCI*DDOT(NDOCC,WRK3(1,IO),1,FB  (1,IO),1)
     *              + EI*OCCI*DDOT(NDOCC,WRK3(1,IO),1,FBC2(1,IO),1)
        END DO
      ELSE IF (MODE.EQ.6) THEN
        CONTINUE !! ALREADY DONE IN MODE=5
      END IF
C
C
C
      !! - F_{ii}^a * (F_{mi}^b * F_{mi}^c)
      DO IO = 1, NDOCC
        OCCI = OCC(IO)
        IF (NSPIN.EQ.1) OCCI = OCCI*HALF
        IF (MODE.EQ.2.OR.MODE.EQ.4.OR.MODE.EQ.6) OCCI = OCCI*EIG(IO)
        VAL = VAL - WRK1(IO,IO)*OCCI*DDOT(L1-NDOCC,
     *              FB(NDOCC+1,IO),1,FC(NDOCC+1,IO),1)
      END DO
      DO IO = NDOCC+1, NOCC
        OCCI = OCC(IO)
        IF (NSPIN.EQ.1) OCCI = OCCI*HALF
        IF (MODE.EQ.2.OR.MODE.EQ.4.OR.MODE.EQ.6) OCCI = OCCI*EIG(IO)
        VAL = VAL - WRK1(IO,IO)*OCCI*DDOT(L1,
     *              FB(1,IO),1,FC(1,IO),1)
      END DO
C       if (mode.eq.2) write (*,*) "val2 = ", val
      IF (MODE.EQ.2.OR.MODE.EQ.4.OR.MODE.EQ.6) THEN
        NSEQ = 1
        DO IO = 1, NDOCC
          OCCI = OCC(IO)
          IF (NSPIN.EQ.1) OCCI = OCCI*HALF
          CALL DFTB_CNVSQ(IO,IO,NSEQ)
          VAL = VAL + OCCI*FDER(NSEQ)*
     *                DDOT(L1-NDOCC,WRK1(NDOCC+1,IO),1,FB(NDOCC+1,IO),1)
        END DO
        DO IO = NDOCC+1, NOCC
          OCCI = OCC(IO)
          IF (NSPIN.EQ.1) OCCI = OCCI*HALF
          CALL DFTB_CNVSQ(IO,IO,NSEQ)
          VAL = VAL + OCCI*FDER(NSEQ)*
     *                DDOT(L1,WRK1(1,IO),1,FB(1,IO),1)
        END DO
      END IF
C
      DO IDEG0 = 1, NDEG2(ISPIN)
        IF (ISPIN.EQ.1) IDEG = IDEG0
        IF (ISPIN.EQ.2) IDEG = IDEG0 + NDEG2(1)
        IO = LSTDEG2(1,IDEG)
        JO = LSTDEG2(2,IDEG)
        IF (IO.LE.NDOCC.AND.JO.LE.NDOCC) CYCLE
        IF (IO.GT.NOCC.AND.JO.GT.NOCC) CYCLE
        IF (MODE.EQ.1.OR.MODE.EQ.3.OR.MODE.EQ.5) THEN
          OCCI = OCC(IO)
          IF (NSPIN.EQ.1) OCCI = OCCI*HALF
          VAL = VAL - DDOT(L1,FB(1,IO),1,FC(1,JO),1)*WRK1(IO,JO)*OCCI
          OCCJ = OCC(JO)
          IF (NSPIN.EQ.1) OCCJ = OCCJ*HALF
          VAL = VAL - DDOT(L1,FB(1,JO),1,FC(1,IO),1)*WRK1(JO,IO)*OCCJ
        ELSE IF (MODE.EQ.2.OR.MODE.EQ.4.OR.MODE.EQ.6) THEN
          OCCI = OCC(IO)
          IF (NSPIN.EQ.1) OCCI = OCCI*HALF
          EI = EIG(IO)
          VAL = VAL - DDOT(L1,FB(1,IO),1,FC(1,JO),1)*WRK1(IO,JO)*OCCI*EI
          OCCJ = OCC(JO)
          IF (NSPIN.EQ.1) OCCJ = OCCJ*HALF
          EJ = EIG(JO)
          VAL = VAL - DDOT(L1,FB(1,JO),1,FC(1,IO),1)*WRK1(JO,IO)*OCCJ*EJ
          CALL DFTB_CNVSQ(IO,JO,NSEQ)
          VAL = VAL + OCCI*FDER(NSEQ)*DDOT(L1,WRK1(1,IO),1,FB(1,JO),1)
          VAL = VAL + OCCJ*FDER(NSEQ)*DDOT(L1,WRK1(1,JO),1,FB(1,IO),1)
        END IF
      END DO
C
      IF (MODE.EQ.1.OR.MODE.EQ.3.OR.MODE.EQ.5) THEN
        VAL = TMP + VAL
      ELSE IF (MODE.EQ.2.OR.MODE.EQ.4.OR.MODE.EQ.6) THEN
        VAL = TMP - VAL
      END IF
C
      RETURN
C
      END SUBROUTINE DFTB_POLDER1
C
C-----------------------------------------------------------------------
C*MODULE DFTBHS    *DECK DFTB_DEGCLR
C>
C>    @brief   Calculate third-order derivatives
C>
C>    @details Put zero for quasi-degenerated pairs
C>
C>    @author  Yoshio Nishimoto
C>    - Feb, 2018- Subroutine written
C>
C>           --- INPUT ---
C>    @param L1      Number of basis functions
C>    @param NDEG    Number of quasi-degenerated pairs
C>    @param ISPIN   Current index of spin
C>    @param LSTDEG  Index of quasi-degenerated pairs
C>           --- IN/OUTPUT ---
C>    @param WRK     The matrix to be put zero
C>
C
      SUBROUTINE DFTB_DEGCLR(L1,NDEG,ISPIN,LSTDEG,WRK)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (ZERO=0.0D+00)
      DIMENSION NDEG(3),LSTDEG(2,*),WRK(L1,*)
C
      IF (ISPIN.EQ.1) THEN
        DO IDEG = 1, NDEG(1)
          IO = LSTDEG(1,IDEG)
          JO = LSTDEG(2,IDEG)
          WRK(JO,IO) = ZERO
          WRK(IO,JO) = ZERO
        END DO
      ELSE
        DO IDEG = 1, NDEG(2)
          IO = LSTDEG(1,NDEG(1)+IDEG)
          JO = LSTDEG(2,NDEG(1)+IDEG)
          WRK(JO,IO) = ZERO
          WRK(IO,JO) = ZERO
        END DO
      END IF
C
      DO IO = 1, L1
        WRK(IO,IO) = ZERO
      END DO
C
      END SUBROUTINE DFTB_DEGCLR
C
C-----------------------------------------------------------------------
C*MODULE DFTBHS    *DECK DFTB_CPYTSQ
C>
C>    @brief   Similar to CPYTSQ
C>
C>    @details This is similar to CPYTSQ, but only written with DCOPY,
C>             resulting in approximately 30% speed up (?). In addition,
C>             it is possible to restrict the number MOs to be copied.
C>
C>    @author  Yoshio Nishimoto
C>    - Feb, 2018- Subroutine written
C>
C>           --- INPUT ---
C>    @param A       Triangular matrix
C>    @param NA      The size (leading dimension) of A and B
C>    @param NB      Number of MOs to be copied.
C>           --- OUTPUT ---
C>    @param B       Square matrix
C
      SUBROUTINE DFTB_CPYTSQ(A,B,NA,NB)
C     CLONE OF CPYTSQT, BUT THIS TAKES 2/3
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION A(*),B(NA,NA)
      IJ=1
      DO I = 1, NB
        CALL DCOPY(I,A(IJ),1,B(1,I),1)
        IJ = IJ + I
      END DO
      DO I = 1, NB-1
        CALL DCOPY(NB-I,B(I,I+1),NA,B(I+1,I),1)
      END DO
      RETURN
      END
C
C-----------------------------------------------------------------------
C*MODULE DFTBHS    *DECK DFTB_3RDAR
C>
C>    @brief   Calculate third-order derivatives
C>
C>    @details Compute the contributions of three first-order response.
C>             This is only for DFTB3.
C>
C>    @author  Yoshio Nishimoto
C>    - Feb, 2018- Subroutine written
C>
C>           --- INPUT ---
C>    @param NAT     Number of atoms
C>    @param NDIM    Number of atoms or shells (SRSCC=T)
C>    @param NSPE    Number of species
C>    @param ISPE    Index of species
C>    @parma MAXANG  Maximum angular momentum of each species
C>    @param DERMUL  Derivative of Mulliken charges
C>    @param CHAMUL  Atom-resolved Mulliken charges
C>    @param CHAMULS Shell-resolved Mulliken charges
C>    @param GAMMA   Matrix of DFTB3 gamma
C>    @param WRK     Work space
C>    @param DFTB3   Whether DFTB3 or not
C>    @param SRSCC   Whether Shell-resolved SCC or not
C>           --- OUTPUT ---
C>    @param SHIFT   Shift contribution of each atom or shell
C
      SUBROUTINE DFTB_3RDAR(NAT,NDIM,NSHELL,NSPE,ISPE,MAXANG,SHIFT,
     *                      DERMUL,CHAMUL,CHAMULS,GAMMA,WRK,
     *                      DFTB3,SRSCC)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER (ONE_THIRD=1.0D+00/3.0D+00,TWO_THIRD=2.0D+00/3.0D+00)
C
      DIMENSION ISPE(*),MAXANG(*),SHIFT(*),DERMUL(*),CHAMUL(*),
     *          CHAMULS(*),GAMMA(NDIM,NDIM),WRK(*)
      LOGICAL   DFTB3,SRSCC
C
      IF (.NOT.DFTB3) RETURN
C
      IF (SRSCC) THEN
        !! NEED ATOM-RESOLVED
        CALL DCOPY(NDIM,DERMUL,1,WRK(NDIM+1),1)
        CALL DFTB_SHIFT_SHTOA(WRK(NDIM+1),WRK,NAT,NSHELL,NSPE,ISPE,
     *                        MAXANG)
C
        ISH0 = 0
        DO I = 1, NAT
          JSH0 = 0
          DERMULI = WRK(I)
          DO J = 1, NAT
            DERMULJ = WRK(J)
            DO II = 1, MAXANG(ISPE(I))
              ISH = ISH0 + II
              DO JJ = 1, MAXANG(ISPE(J))
                JSH = JSH0 + JJ
                DERSHLJ = DERMUL(JSH)
                SHIFT(ISH) = SHIFT(ISH)
     *            + ONE_THIRD*
     *             (DERMULI*CHAMULS(JSH)*GAMMA(ISH,JSH)
     *            + DERSHLJ*CHAMUL(I)   *GAMMA(ISH,JSH)
     *            + DERMULJ*CHAMULS(JSH)*GAMMA(JSH,ISH)
     *            + DERSHLJ*CHAMUL(J)   *GAMMA(JSH,ISH))
                DO KK = 1, MAXANG(ISPE(I))
                  KSH = ISH0 + KK
                  DERSHLK = DERMUL(KSH)
                  SHIFT(ISH) = SHIFT(ISH)
     *              + ONE_THIRD*
     *               (DERSHLK*CHAMULS(JSH)*GAMMA(KSH,JSH)
     *              + DERSHLJ*CHAMULS(KSH)*GAMMA(KSH,JSH))
                END DO
              END DO
            END DO
            JSH0 = JSH0 + MAXANG(ISPE(J))
          END DO
          ISH0 = ISH0 + MAXANG(ISPE(I))
        END DO
      ELSE
        DO I = 1, NDIM !! a
          DO J = 1, NDIM !! c
            GAMI = GAMMA(I,J) !! G_ac
            GAMJ = GAMMA(J,I) !! G_ca
            SHIFT(I) = SHIFT(I)
     *        + TWO_THIRD*(GAMI*DERMUL(I)*CHAMUL(J)
     *          +GAMI*DERMUL(J)*CHAMUL(I)+GAMJ*DERMUL(J)*CHAMUL(J))
          END DO
        END DO
      END IF
C
      RETURN
C
      END SUBROUTINE DFTB_3RDAR
C
C-----------------------------------------------------------------------
C*MODULE DFTBHS    *DECK DFTB_CLCDEG
C>
C>    @brief   Calculate quasi-degenerated coefficient
C>
C>    @details Compute quasi-degenerated coefficient, such as
C>             (fi*ei-fj*ej)/(ei-ej), by calling DFTB_CLCDEG2 subroutine
C>             recursively.
C>
C>    @author  Yoshio Nishimoto
C>    - Feb, 2018- Subroutine written
C>
C>           --- INPUT ---
C>    @param N0      Order of degeneracy
C>    @param L0      Derivative order of fi
C>    @param M0      Exponent of ei
C>    @param BETA    1/kTe
C>    @param EIG     Eigenvalues of each orbital
C>    @param OCC     Occupation number of each orbital
C
      FUNCTION DFTB_CLCDEG(N0,L0,M0,BETA,EIG,OCC)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER (NADDMAX=3)
      DIMENSION DEGTMP((L0+1)*(N0+NADDMAX)**2),EIG(*),OCC(*)
      INTEGER   DFTB_ICOEFF
C
      DFTB_CLCDEG = 0.0D+00
C
      OCCMIN = 1.0D+00
      OCCMAX = 0.0D+00
      DO IO = 1, N0+1
        IF (OCC(IO).GT.OCCMAX) OCCMAX = OCC(IO)
        IF (OCC(IO).LT.OCCMIN) OCCMIN = OCC(IO)
      END DO
      NADD = NADDMAX
      IF (OCCMAX-OCCMIN.LE.1.0D-09) NADD = 2 - N0
      IF (OCCMAX-OCCMIN.LE.1.0D-15) NADD = 1 - N0
C
      DO LL = L0+1, 1, -1
        NORDER = N0
        CALL DFTB_CLCDEG2(TMP,NORDER,LL,M0,BETA,EIG,OCC,DEGTMP,.TRUE.,
     *                    L0+1,N0+NADD)
        DFTB_CLCDEG = DFTB_CLCDEG + DBLE(DFTB_ICOEFF(L0,LL))*TMP
      END DO
C
      IF (L0.NE.0) DFTB_CLCDEG = DFTB_CLCDEG * BETA**L0
C
      RETURN
C
      END FUNCTION DFTB_CLCDEG
C
C-----------------------------------------------------------------------
C*MODULE DFTBHS    *DECK DFTB_CLCDEG2
C>
C>    @brief   Calculate quasi-degenerated coefficient
C>
C>    @details Compute quasi-degenerated coefficient, such as
C>             (fi*ei-fj*ej)/(ei-ej).
C>
C>    @author  Yoshio Nishimoto
C>    - Feb, 2018- Subroutine written
C>
C>           --- INPUT ---
C>    @param N0      Order of degeneracy
C>    @param LL      Derivative order of fi
C>    @param MM      Exponent of ei
C>    @param BETA    1/kTe
C>    @param EIG     Eigenvalues of each orbital
C>    @param OCC     Occupation number of each orbital
C>    @param DEGTMP  Temporary contribution
C>    @param FIRST   Is this the initial call or not
C>    @param LMAX    Specifies the dimension of DEGTMP
C>    @param NMAX    Specifies the dimension of DEGTMP
C>           --- OUTPUT ---
C>    @param VAL     Result
C>
C
      RECURSIVE SUBROUTINE DFTB_CLCDEG2(VAL,N0,LL,MM,BETA,EIG,OCC,
     *                                  DEGTMP,FIRST,LMAX,NMAX)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER (ZERO=0.0D+00,ONE=1.0D+00)
C
      COMMON /FACTS / FACT(0:100),NFACT
C
      DIMENSION EIG(N0+1),OCC(N0+1),DEGTMP(NMAX,LMAX,*)
C
      LOGICAL FIRST
C
      VALTMP = ZERO
C
      NORDER = N0
      IF (FIRST) THEN
        NORDER = N0+1
        M0 = MM
      END IF
C
      IF (NORDER.EQ.1) THEN
        OCCX = ONE
        IF (LL.NE.0) OCCX = OCC(1)**LL
        EIGX = ONE
        IF (MM.NE.0) EIGX = EIG(1)**MM
        VALTMP = OCCX*EIGX
      ELSE
        IF ((FIRST.AND.MM.EQ.M0 .OR. .NOT.FIRST.AND.MM.EQ.0)
     *    .AND.LL.EQ.LMAX) THEN
          DO II = 1, NMAX
            DO KK = LMAX, 0, -1
              CALL DFTB_CLCDEG2(DEGTMP(II,KK+1,1),NORDER-1,KK,II-1,BETA,
     *                          EIG,OCC,DEGTMP(1,1,2),.FALSE.,LMAX,NMAX)
            END DO
          END DO
        END IF

        EIGX = ONE
        DO II = 0, MM-1
          TMP = DEGTMP(MM-II,LL+1,1)
          VALTMP = VALTMP + TMP*EIGX
          EIGX = EIGX*EIG(NORDER)
        END DO
  
        DO KK = 0, LL-1
          OCCX = ONE
          IF (KK.NE.0) OCCX = OCC(NORDER)**KK
          DO II = 1, NMAX
            SCAL = -OCCX*(ONE-OCC(NORDER)) * BETA**II * EIGX / FACT(II)
            EIGX2 = ONE
            DO JJ = 0, II-1
              TMP = DEGTMP(II-JJ,LL-KK+1,1)
              CALL DFTB_BINOM(NVAL,II-1,JJ)
              CC = DBLE(NVAL)
              VALTMP = VALTMP + TMP * SCAL * CC * EIGX2
              EIGX2 = -EIGX2*EIG(NORDER)
            END DO
          END DO
        END DO
      END IF
      VAL = VALTMP
C
      END SUBROUTINE DFTB_CLCDEG2
C
C-----------------------------------------------------------------------
C*MODULE DFTBHS    *DECK DFTB_ICOEFF
C>
C>    @brief   Calculate quasi-degenerated coefficient
C>
C>    @details Compute the coefficient of the quasi-degenerate coeff.
C>             for L0 > 0
C>
C>    @author  Yoshio Nishimoto
C>    - Feb, 2018- Subroutine written
C>
C>           --- INPUT ---
C>    @param N0      Order of degeneracy
C>    @param LL      Derivative order of fi
C>    @param MM      Exponent of ei
C>    @param BETA    1/kTe
C>    @param EIG     Eigenvalues of each orbital
C>    @param OCC     Occupation number of each orbital
C>    @param DEGTMP  Temporary contribution
C>    @param FIRST   Is this the initial call or not
C>    @param LMAX    Specifies the dimension of DEGTMP
C>    @param NMAX    Specifies the dimension of DEGTMP
C>           --- OUTPUT ---
C>    @param VAL     Result
C>
C
      RECURSIVE INTEGER FUNCTION DFTB_ICOEFF(NN,MM) RESULT(II)
C
      INTEGER NN,MM,NFACT
      DOUBLE PRECISION FACT
C
      COMMON /FACTS / FACT(0:100),NFACT
C
      !! NN > 0
      !! 1 <= M <= N+1
      IF (MM.EQ.NN+1) THEN
        II = INT(FACT(NN))
      ELSE IF (MM.EQ.1) THEN
        II = (-1)**NN
      ELSE
        II = (MM-1)*DFTB_ICOEFF(NN-1,MM-1) - MM*DFTB_ICOEFF(NN-1,MM)
      END IF
C
      RETURN
C
      END FUNCTION
C
C-----------------------------------------------------------------------
C*MODULE DFTBHS    *DECK DFTB_BINOM
C>
C>    @brief   Calculate quasi-degenerated coefficient
C>
C>    @details Compute binomial coefficient (nCm)
C>
C>    @author  Yoshio Nishimoto
C>    - Feb, 2018- Subroutine written
C>
C>           --- INPUT ---
C>    @param N       Left index of nCm
C>    @param M       Right index of nCm
C>           --- OUTPUT ---
C>    @param NVAL    Result
C>
C
      SUBROUTINE DFTB_BINOM(NVAL,N,M)
      INTEGER NVAL,N,M,NN
      IF (N.EQ.0.OR.M.EQ.0) THEN
        NVAL = 1
      ELSE IF (N.LT.M) THEN
        NVAL = 0
      ELSE
        NVAL = 1
        DO NN = M, 1, -1
          NVAL = NVAL*(N+1-NN)
        END DO
        DO NN = M, 1, -1
          NVAL = NVAL/(M+1-NN)
        END DO
      END IF
      END SUBROUTINE DFTB_BINOM
C
C-----------------------------------------------------------------------
C*MODULE DFTBHS    *DECK DFTB_LCHESS
C>
C>    @brief   Calculate LC terms.
C>
C>    @details Compute LC terms.
C>
C>    @author  Yoshio Nishimoto
C>
C
      SUBROUTINE DFTB_LCHESS(IAT,JAT,DIST,GRAD,HESS)
      use mx_limits, only: MXATM
      USE lrcdft, ONLY: LCFLAG, EMU, EMU2, LRFILE
C
      IMPLICIT NONE
C     Declarations of arguments
      DOUBLE PRECISION DIST,GRAD,HESS
      INTEGER IAT,JAT
C     Declarations of common blocks
      DOUBLE PRECISION ETEMP,DFTBDP,DAMPXHE,HUBDER
      DOUBLE PRECISION ZREF,SPNCST,SPE
      DOUBLE PRECISION PARAMDIR
      INTEGER NSPE
      INTEGER MAXANG,ISPE,IND,IDFTBD
      DOUBLE PRECISION SKDIM,SKSPIN,QREFL
      DOUBLE PRECISION HUBBL,QREF,HUBB,SKCUT2
      INTEGER NEEDSK,LSKHTAB,LSKSTAB
      INTEGER LSKGRID,LSKSELF
C     Other declarations
      DOUBLE PRECISION aa,bab,bb,bba,f0,f1,f2,f3,gab,gba
      INTEGER ii
      DOUBLE PRECISION ta,tau,tb,tmpg,tmph,ww
C
      INTEGER,PARAMETER :: MXSPE=10
C
      COMMON /DFTBPR/ ETEMP,DFTBDP(MXSPE*14),DAMPXHE,HUBDER(MXSPE),
     *                ZREF(MXATM),SPNCST(6,MXSPE),SPE(MXATM),NSPE,
     *                MAXANG(MXATM),ISPE(MXATM),IND(MXATM+1),IDFTBD,
     *                PARAMDIR
      COMMON /DFTBSK/ SKDIM(MXSPE,MXSPE),SKSPIN(MXSPE),QREFL(3,MXSPE),
     *                HUBBL(3,MXSPE),QREF(MXSPE),HUBB(MXSPE),SKCUT2,
     *                NEEDSK,LSKHTAB(MXSPE,MXSPE),LSKSTAB(MXSPE,MXSPE),
     *                LSKGRID(MXSPE,MXSPE),LSKSELF(MXSPE)
C
C     VVV  = 0.0D+00
      GRAD = 0.0D+00
      HESS = 0.0D+00
      IF (IAT.EQ.JAT .OR. DIST.LE.1.0D-05) RETURN
C     if (ww.eq.0.0d+00) return
C
      TA = HUBB(ISPE(IAT))*3.2D+00
      TB = HUBB(ISPE(JAT))*3.2D+00
C
      DO II = 1, 2
        IF (II.EQ.1) WW = 0.0D+00
        IF (II.EQ.2) WW = EMU
        IF (ABS(TA-TB).LT.0.3125D-05) THEN
          TAU = 0.5D+00*(TA+TB)
          AA = TAU**8/((TAU*TAU-WW*WW)**4)
          BB = TAU**3/(4.8D+01*(TAU*TAU-WW*WW)**4)
          F3 = TAU*TAU*(TAU*TAU-WW*WW)**3
          F2 = 3.0D+00*TAU*(TAU*TAU-WW*WW)**2*(3.0D+00*TAU*TAU-WW*WW)
          F1 = 3.0D+00*(TAU*TAU-WW*WW)
     *         *(1.1D+01*TAU**4-4.0D+00*TAU*TAU*WW*WW+WW**4)
          F0 = 4.8D+01*TAU**5
C
C         TMPV = AA*EXP(-WW*DIST)/DIST
C    *         - BB*EXP(-TAU*DIST)/DIST*(
C    *            F3*DIST**3 + F2*DIST**2
C    *           +F1*DIST + F0)
          TMPG = -AA*EXP(-WW*DIST)/DIST*(WW+1.0D+00/DIST)
     *         + BB*EXP(-TAU*DIST)/DIST*(F3*TAU*DIST**3
     *           +(F2*TAU-2.0D+00*F3)*DIST**2
     *           +(F1*TAU-F2)*DIST
     *           + F0*TAU + F0/DIST)
          TMPH = AA*EXP(-WW*DIST)/DIST
     *           *(WW*WW+2.0D+00*WW/DIST+2.0D+00/DIST**2)
     *         - BB*EXP(-TAU*DIST)/DIST*(
     *            F3*TAU*TAU*DIST**3
     *           +(F2*TAU-4.0D+00*F3)*TAU*DIST**2
     *           +(F1*TAU*TAU-2.0D+00*F2*TAU+2.0D+00*F3)*DIST
     *           +(F0*TAU*TAU) + 2.0D+00*F0*TAU/DIST
     *           +2.0D+00*F0/DIST**2)
        ELSE
          AA = TA**4*TB**4/((TA*TA-WW*WW)**2*(TB*TB-WW*WW)**2)
          BAB= TA*TA/(TA*TA-WW*WW)*TA*TB**4/((TB*TB-TA*TA)**2)*0.5D+00
          BBA= TB*TB/(TB*TB-WW*WW)*TB*TA**4/((TA*TA-TB*TB)**2)*0.5D+00
          GAB= TA**4/((WW*WW-TA*TA)**2)
     *         *(TB**6-3.0D+00*TA**2*TB**4+2.0D+00*WW*WW*TB**4)
     *         /((TA*TA-TB*TB)**3)
          GBA= TB**4/((WW*WW-TB*TB)**2)
     *         *(TA**6-3.0D+00*TB**2*TA**4+2.0D+00*WW*WW*TA**4)
     *         /((TB*TB-TA*TA)**3)
C
C         TMPV = EXP(-WW*DIST)/DIST*AA
C    *         - EXP(-TA*DIST)*(BAB-GAB/DIST)
C    *         - EXP(-TB*DIST)*(BBA-GBA/DIST)
          TMPG = -EXP(-WW*DIST)/DIST*(WW+1.0D+00/DIST)*AA
     *         + EXP(-TA*DIST)*(TA*BAB-TA*GAB/DIST-GAB/DIST**2)
     *         + EXP(-TB*DIST)*(TB*BBA-TB*GBA/DIST-GBA/DIST**2)
          TMPH = EXP(-WW*DIST)/DIST
     *           *(WW*WW+2.0D+00*WW/DIST+2.0D+00/DIST**2)*AA
     *         - EXP(-TA*DIST)*(TA*TA*BAB-TA*TA*GAB/DIST
     *           -2.0D+00*TA*GAB/DIST**2-2.0D+00*GAB/DIST**3)
     *         - EXP(-TB*DIST)*(TB*TB*BBA-TB*TB*GBA/DIST
     *           -2.0D+00*TB*GBA/DIST**2-2.0D+00*GBA/DIST**3)
        END IF
        IF (II.EQ.1) THEN
C         VVV  = TMPV
          GRAD = TMPG
          HESS = TMPH
        ELSE IF (II.EQ.2) THEN
C         VVV  = VVV  - TMPV
          GRAD = GRAD - TMPG
          HESS = HESS - TMPH
        END IF
      END DO
C
      RETURN
C
      END SUBROUTINE DFTB_LCHESS
C
C-----------------------------------------------------------------------
C*MODULE DFTBHS    *DECK DFTB_LCHESS1
C>
C>    @brief   Calculate LC terms.
C>
C>    @details Compute LC terms.
C>
C>    @author  Yoshio Nishimoto
C>
C
      SUBROUTINE DFTB_LCHESS1(NI,KI,FAO,SDER,S,DELD,GAMLC,GAMDERLC,DIST)
      use mx_limits, only: MXATM
      use constants, only: ZERO, ONE
C
      IMPLICIT NONE
C     Declarations of arguments
      DOUBLE PRECISION FAO(*),SDER(*),S(*),DELD(*),GAMLC(NAT,NAT),
     *                 GAMDERLC(*),DIST(*)
      INTEGER NI,KI
C     Declarations of common blocks
      DOUBLE PRECISION ETEMP,DFTBDP,DAMPXHE,HUBDER
      DOUBLE PRECISION ZREF,SPNCST,SPE
      DOUBLE PRECISION PARAMDIR
      INTEGER NSPE
      INTEGER MAXANG,ISPE,IND,IDFTBD
      DOUBLE PRECISION X
      DOUBLE PRECISION ZAN, C
      INTEGER NAT,ICH,MUL,NUM,NQMT,NE,NA,NB
      INTEGER IAN
C     Other declarations
      INTEGER, PARAMETER :: MXSPE=10
      DOUBLE PRECISION gamder
      INTEGER i,iat,l1,l2,l3,last,ldsq,lgamma,loadfm,lss,lssder
      INTEGER lwrk,lwrk1,lwrk2,mu,need,ngotmx,nj,nseq,nu
C
      COMMON /DFTBPR/ ETEMP,DFTBDP(MXSPE*14),DAMPXHE,HUBDER(MXSPE),
     *                ZREF(MXATM),SPNCST(6,MXSPE),SPE(MXATM),NSPE,
     *                MAXANG(MXATM),ISPE(MXATM),IND(MXATM+1),IDFTBD,
     *                PARAMDIR
      COMMON /FMCOM / X(1)
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
C
      L1=NUM
      L2=L1*(L1+1)/2
      L3=L1*L1
      CALL VALFM(LOADFM)
      CALL GOTFM(NGOTMX)
      LSS     = LOADFM  + 1
      LSSDER  = LSS     + L3
      LDSQ    = LSSDER  + L3
      LWRK    = LDSQ    + L3
      LWRK1   = LWRK    + L3
      LWRK2   = LWRK1   + L3
      LGAMMA  = LWRK2   + L3
      LAST    = LGAMMA  + L3
      NEED    = LAST - LOADFM - 1
      CALL GETFM(NEED)
C
      CALL CPYTSQ(S,X(LSS),L1,1)
      CALL CPYTSQ(SDER,X(LSSDER),L1,1)
      CALL CPYTSQ(DELD,X(LDSQ),L1,1)
C
C     ----- gamma derivative contributions -----
C
      !! Prepare gamma derivative matrix
      CALL VCLR(X(LGAMMA),1,L3)
      DO NJ = 1, NAT
        IF (NI.EQ.NJ) CYCLE
        CALL DFTB_CNVSQ(NI,NJ,NSEQ)
        GAMDER = GAMDERLC(NSEQ)*(C(KI,NI)-C(KI,NJ))/DIST(NSEQ)
        DO MU = IND(NI)+1, IND(NI)+MAXANG(ISPE(NI))**2
          DO NU = IND(NJ)+1, IND(NJ)+MAXANG(ISPE(NJ))**2
            X(LGAMMA+MU-1+L1*(NU-1)) = GAMDER
            X(LGAMMA+NU-1+L1*(MU-1)) = GAMDER
          END DO
        END DO
      END DO
C
      !! 1) dgamma_{AB}/da
      CALL DGEMM('N','N',L1,L1,L1,ONE,X(LSS),L1,X(LDSQ),L1,
     *           ZERO,X(LWRK1),L1)
      CALL DGEMM('N','N',L1,L1,L1,ONE,X(LWRK1),L1,X(LSS),L1,
     *           ZERO,X(LWRK),L1)
      DO I = 1, L3
        X(LWRK+I-1) = X(LWRK+I-1)*X(LGAMMA+I-1)
      END DO
      !! 2) dgamma_{AD}/da + dgamma_{BC}/da
      DO I = 1, L3
        X(LWRK1+I-1) = X(LWRK1+I-1)*X(LGAMMA+I-1)*2.0D+00
      END DO
      !! 4) dgamma_{BC}/da
      DO I = 1, L3
        X(LDSQ+I-1) = X(LDSQ+I-1)*X(LGAMMA+I-1)
      END DO
      CALL DGEMM('N','N',L1,L1,L1,ONE,X(LSS),L1,X(LDSQ),L1,
     *           ONE,X(LWRK1),L1)
      CALL DGEMM('N','N',L1,L1,L1,ONE,X(LWRK1),L1,X(LSS),L1,
     *           ONE,X(LWRK),L1)
C
C     ---- S derivative contributions -----
C
      CALL CPYTSQ(DELD,X(LDSQ),L1,1)
      !! Prepare gamma matrix
      CALL VCLR(X(LGAMMA),1,L3)
      DO IAT = 1, NAT
        DO NJ = 1, NAT
          CALL DFTB_CNVSQ(IAT,NJ,NSEQ)
          DO MU = IND(IAT)+1, IND(IAT)+MAXANG(ISPE(IAT))**2
            DO NU = IND(NJ)+1, IND(NJ)+MAXANG(ISPE(NJ))**2
              !! Scale with 2, because S^a*D*S + S*D*S^a
              X(LGAMMA+MU-1+L1*(NU-1)) = GAMLC(IAT,NJ)*2.0D+00
            END DO
          END DO
        END DO
      END DO
C
      !! 1) gamma_{AB}/da
      CALL DGEMM('N','N',L1,L1,L1,ONE,X(LSSDER),L1,X(LDSQ),L1,
     *           ZERO,X(LWRK1),L1)
      CALL DGEMM('N','N',L1,L1,L1,ONE,X(LWRK1),L1,X(LSS),L1,
     *           ZERO,X(LWRK2),L1)
      DO I = 1, L3
        X(LWRK+I-1) = X(LWRK+I-1) + X(LWRK2+I-1)*X(LGAMMA+I-1)
      END DO
      !! 2) gamma_{AD}/da
      DO I = 1, L3
        X(LWRK1+I-1) = X(LWRK1+I-1)*X(LGAMMA+I-1)
      END DO
      !! 4) gamma_{BC}/da
      DO I = 1, L3
        X(LDSQ+I-1) = X(LDSQ+I-1)*X(LGAMMA+I-1)
      END DO
      CALL DGEMM('N','N',L1,L1,L1,ONE,X(LSSDER),L1,X(LDSQ),L1,
     *           ONE,X(LWRK1),L1)
      CALL DGEMM('N','N',L1,L1,L1,ONE,X(LWRK1),L1,X(LSS),L1,
     *           ONE,X(LWRK),L1)
      !! 3) gamma_{BC}/da
      CALL CPYTSQ(DELD,X(LDSQ),L1,1)
      CALL DGEMM('N','N',L1,L1,L1,ONE,X(LDSQ),L1,X(LSS),L1,
     *           ZERO,X(LWRK1),L1)
      DO I = 1, L3
        X(LWRK1+I-1) = X(LWRK1+I-1)*X(LGAMMA+I-1)
      END DO
      CALL DGEMM('N','N',L1,L1,L1,ONE,X(LSSDER),L1,X(LWRK1),L1,
     *           ONE,X(LWRK),L1)
C
      CALL SYMSQT(L1,X(LWRK),L1,FAO,1)
      CALL DSCAL(L2,-0.125D+00,FAO,1)
C
      CALL RETFM(NEED)
C
      RETURN
C
      END SUBROUTINE DFTB_LCHESS1
C
C-----------------------------------------------------------------------
C*MODULE DFTBHS    *DECK DFTB_LCHESS2
C>
C>    @brief   Calculate LC terms.
C>
C>    @details Compute LC terms.
C>
C>    @author  Yoshio Nishimoto
C>
C
      SUBROUTINE DFTB_LCHESS2(NI0,KI0,NSPIN,FAO,SDER,S,DELD,DDER,GAMLC,
     *                        GAMDERLC,DIST,MAXBND,INDBND,
     *                        DERMAT,MASWRK,HESSIAN,NPROC)
      use mx_limits, only: MXATM
      use constants, only: ZERO, ONE
C
      IMPLICIT NONE
C     Declarations of arguments
      INTEGER NI0,KI0,NSPIN
      DOUBLE PRECISION FAO,SDER,S,DELD,DDER,GAMLC
      DOUBLE PRECISION GAMDERLC,DIST
      INTEGER MAXBND,INDBND
      DOUBLE PRECISION DERMAT
      LOGICAL MASWRK
      DOUBLE PRECISION HESSIAN
      INTEGER NPROC
C     Declarations of common blocks
      DOUBLE PRECISION ETEMP,DFTBDP,DAMPXHE,HUBDER
      DOUBLE PRECISION ZREF,SPNCST,SPE
      DOUBLE PRECISION PARAMDIR
      INTEGER NSPE
      INTEGER MAXANG,ISPE,IND,IDFTBD
      DOUBLE PRECISION X
      DOUBLE PRECISION ZAN, C
      INTEGER NAT,ICH,MUL,NUM,NQMT,NE,NA,NB
      INTEGER IAN
C     Other declarations
      DOUBLE PRECISION, PARAMETER :: TWO=2.0D+00
      INTEGER, PARAMETER :: MXSPE=10
      DOUBLE PRECISION dsdr,gamder
      INTEGER ibnd,ir,jr,ki,l1,l2,l3,last,lddsq,ldsq,lgamma
      INTEGER loadfm,lss,lssder,lwrk,lwrk1,lwrk2,m,mu,n,need
      INTEGER ngotmx,ni,nj,nseq,nseq2,nu
      DOUBLE PRECISION val,val2
C
      COMMON /DFTBPR/ ETEMP,DFTBDP(MXSPE*14),DAMPXHE,HUBDER(MXSPE),
     *                ZREF(MXATM),SPNCST(6,MXSPE),SPE(MXATM),NSPE,
     *                MAXANG(MXATM),ISPE(MXATM),IND(MXATM+1),IDFTBD,
     *                PARAMDIR
      COMMON /FMCOM / X(1)
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
C
      DIMENSION FAO(*),SDER(*),S(*),DELD(*),DDER(*),GAMLC(NAT,NAT),
     *          GAMDERLC(*),DIST(*),INDBND(MAXBND+1,4,NAT),
     *          DERMAT(2,*),HESSIAN(3*NAT,3*NAT)
C
      JR = 3*(NI0-1)+KI0
      L1=NUM
      L2=L1*(L1+1)/2
      L3=L1*L1
      CALL VALFM(LOADFM)
      CALL GOTFM(NGOTMX)
      LSS     = LOADFM  + 1
      LSSDER  = LSS     + L3
      LDSQ    = LSSDER  + L3
      LDDSQ   = LDSQ    + L3
      LWRK    = LDDSQ   + L3
      LWRK1   = LWRK    + L3
      LWRK2   = LWRK1   + L3
      LGAMMA  = LWRK2   + L3
      LAST    = LGAMMA  + NAT*(NAT+1)/2
      NEED    = LAST - LOADFM - 1
      CALL GETFM(NEED)
C
      CALL CPYTSQ(S,X(LSS),L1,1)
      CALL CPYTSQ(SDER,X(LSSDER),L1,1)
      CALL CPYTSQ(DELD,X(LDSQ),L1,1)
      CALL CPYTSQ(DDER,X(LDDSQ),L1,1)
C
      !! Prepare gamma derivative matrix
      CALL VCLR(X(LGAMMA),1,NAT*(NAT+1)/2)
      DO NJ = 1, NAT
        IF (NI0.EQ.NJ) CYCLE
        CALL DFTB_CNVSQ(NI0,NJ,NSEQ)
        GAMDER = GAMDERLC(NSEQ)*(C(KI0,NI0)-C(KI0,NJ))/DIST(NSEQ)
        X(LGAMMA+NSEQ-1) = GAMDER
      END DO
C
C     ----- S^a * gamma^b ----- 4 terms
C     a: variable
C     b: fixed
C
      CALL DFTB_LCGRAD1(1,NAT,L1,ISPE,MAXANG,IND,X(LDSQ),X(LDSQ),
     *                  S,X(LGAMMA),X(LWRK))
      CALL DSCAL(L2,ONE,X(LWRK),1)
C
C     ----- S^b * gamma^a ----- 4 terms
C
      CALL DFTB_LCGRAD2(3,ONE*DBLE(NSPIN)/NPROC,NAT,L1,ISPE,MAXANG,IND,
     *                  X(LDSQ),X(LDSQ),X(LSS),X(LSSDER),DIST,GAMDERLC,
     *                  C,HESSIAN(1,JR),X(LAST))
      CALL DFTB_LCGRAD2(3,ONE*DBLE(NSPIN)/NPROC,NAT,L1,ISPE,MAXANG,IND,
     *                  X(LDSQ),X(LDSQ),X(LSSDER),X(LSS),DIST,GAMDERLC,
     *                  C,HESSIAN(1,JR),X(LAST))
C
C     ----- D^b terms (1) ----- 2 terms
C     D^b*S*D*S*gamma^a
C
      CALL DFTB_LCGRAD2(1,TWO*DBLE(NSPIN)/NPROC,NAT,L1,ISPE,MAXANG,IND,
     *                  X(LDDSQ),X(LDSQ),X(LSS),X(LSS),DIST,GAMDERLC,
     *                  C,HESSIAN(1,JR),X(LAST))
C
C     ----- D^b terms (2) ----- 4 terms
C     D^b*S^a*D*S*gamma -> S^a*(...)
C
      CALL SYMSQT(NAT,GAMLC,NAT,X(LGAMMA),1)
      CALL DFTB_LCGRAD1(1,NAT,L1,ISPE,MAXANG,IND,X(LDSQ),X(LDDSQ),
     *                  S,X(LGAMMA),X(LWRK1))
      CALL DAXPY(L2,TWO,X(LWRK1),1,X(LWRK),1)
C
C     ----- S^a * S^b terms ----- 2 terms
C     D*S^a*D*S^b*gamma -> S^a*(...)
C
      CALL DFTB_LCGRAD1(0,NAT,L1,ISPE,MAXANG,IND,X(LDSQ),X(LDSQ),
     *                  SDER,X(LGAMMA),X(LWRK1))
      CALL DAXPY(L2,ONE,X(LWRK1),1,X(LWRK),1)
C
      IR = 0
      DO NI = 1, NAT
        DO KI = 1, 3
          IR = IR + 1
          VAL = ZERO
          VAL2= ZERO
          DO IBND = 1, INDBND(1,1,NI)
            NJ = INDBND(IBND+1,1,NI)
            NSEQ2 = INDBND(IBND+1,KI+1,NI)
            DO N = 1, IND(NJ+1) - IND(NJ)
              NU = IND(NJ)+N
              DO M = 1, IND(NI+1) - IND(NI)
                MU = IND(NI)+M
                CALL DFTB_CNVSQ(MU,NU,NSEQ)
                DSDR = DERMAT(2,NSEQ2)
                VAL = VAL + X(LWRK+NSEQ-1)*DSDR
                NSEQ2 = NSEQ2 + 1
              END DO
            END DO
          END DO
          VAL = VAL*TWO
          IF (NSPIN.EQ.2) VAL = VAL*TWO
          IF (MASWRK) HESSIAN(IR,JR) = HESSIAN(IR,JR) + VAL
        END DO
      END DO
C
      CALL RETFM(NEED)
C
      RETURN
C
      END SUBROUTINE DFTB_LCHESS2
