c 17 Oct 19 - DGF - patch FPE division by 0
C  6 Jun 18 - YN -  updates for FMO 5.3
C 18 Apr 16 - YN  - updates for FMO 5.2 
C 29 Sep 14 - YN  - Implement Density-Functional Tight-Binding Method
C
C*MODULE DFTBFO    *DECK DFTB_INPUT_FRAG
C>
C>    @brief Make fragment information for FMO-DFTB
C>
C>    @details Put information for FMO-DFTB, such as the index of
C>             species, index of atomic orbitals, number of basis
C>             functions, labels, nuclear charges (consider BDA),
C>             and number of electrons
C>
C>    @author Yoshio Nishimoto
C>    - Sep, 2014- Subroutine written
C>
C>           --- INPUT ---
C>    @param SPNCST Spin constant array
C>    @param NSPE Number of species
C>    @param MAXANG Maximum angular momentum of each species
C>    @param IW (6)
C>    @param SPE Name of each species
C
C     Possible References (specific to FMO-DFTB):
C     - FMO-DFTB2 (original) + DFTB in GAMESS (?)
C       Nishimoto, Y.; Fedorov, D. G.; Irle, S.
C       J. Chem. Theory Comput. 2014, 10, 4801-4812.
C     - FMO-DFTB3
C       Nishimoto, Y.; Fedorov, D. G.; Irle, S.
C       Chem. Phys. Lett. 2015, 636, 90-96.
C     - FMO-DFTB/SCZV and FMO-DFTB/MD
C       Nishimoto, Y.; Nakata, H.; Fedorov, D. G.; Irle, S.
C       J. Phys. Chem. Lett. 2015, 6, 5034-5039.
C     - FMO-DFTB/Hessian
C       Nakata, H.; Nishimoto, Y.; Fedorov, D. G.
C       J. Chem. Phys. 2016, 145, 044113.
C     - FMO-DFTB/PCM
C       Nishimoto, Y.; Fedorov, D. G.
C       Phys. Chem. Chem. Phys. 2016, 18, 22047-22061.
C     - FMO3-DFTB
C       Nishimoto, Y.; Fedorov, D. G.
C       J. Comput. Chem. 2017, 38, 406-418.
C     - FMO-DFTB/AFO
C       Nishimoto, Y.; Fedorov, D. G.
C       J. Chem. Phys. 2018, 148, 064115.
C     - FMO-LC-DFTB
C       Vuong, V. O. et al.
C       J. Chem. Theory Comput. 2019, 15, 3008-3020.
C
      SUBROUTINE DFTB_INPUT_FRAG(NAT0)
      use mx_limits, only: mxatm,mxao,mxgtot,mxsh
C
      IMPLICIT NONE
C
      INTEGER, INTENT(IN) :: NAT0
      INTEGER,PARAMETER :: MXSPE=10
C
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
      DOUBLE PRECISION SKDIM,SKSPIN,QREFL,HUBBL,QREF,HUBB,SKCUT2
      INTEGER NEEDSK,LSKHTAB,LSKSTAB,LSKGRID,LSKSELF
C
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     &  ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      INTEGER :: NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,IAN
      DOUBLE PRECISION :: ZAN,C,PARAMDIR
C
      COMMON /NSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),
     *                CF(MXGTOT),CG(MXGTOT),CH(MXGTOT),CI(MXGTOT),
     *                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),
     *                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
      DOUBLE PRECISION :: EX,CS,CP,CD,CF,CG,CH,CI
      INTEGER :: KSTART,KATOM,KTYPE,KNG,KLOC,KMIN,KMAX,NSHELL
C
      COMMON /RUNLAB/ TITLE(10),ANAM(MXATM),B(MXATM),BFLAB(MXAO)
      DOUBLE PRECISION :: TITLE,ANAM,B,BFLAB
C     character*8 a(MXATM) 
C     EQUIVALENCE (a,anam)
C
      INTEGER :: I, J, IZAN
      CHARACTER(4) :: DUM4
C
      DO I = NAT0+1, NAT
        DO J = 1, NSPE
          IF (ANAM(I)==SPE(J)) THEN
            ISPE(I) = J
            EXIT
          END IF
        END DO
      END DO
C
C       CALCULATION OF INDICES FOR MATRICES H AND S
C
      IND(1) = 0
      DO I = NAT0+1, NAT
        IND(I+1) = IND(I) + MAXANG(ISPE(I))**2
      END DO
C
      NUM  = IND(NAT+1) !! NUMBER OF CARTESIAN BASIS FUNCTIONS
      NQMT = IND(NAT+1)
C
      NSHELL = 0
      DO I = NAT0+1, NAT
        NSHELL = NSHELL + MAXANG(ISPE(I))
      END DO
C
C     CALL DFTB_SKTAB
      I = 1
      CALL DFTB_SETLAB(I,DUM4)
C
C         SET NUCLEAR CHARGES FOR ALL ATOMS
C
      DO I = NAT0+1, NAT
        ZREF(I) = QREF(ISPE(I))
c     adjust for the BDA
        izan=int(zan(i)+1.0D-02)
        if(ian(i).ne.izan) then
          if(izan.eq.1) then
            zref(i)=1.0D+00 
          else
            zref(i)=zref(i)-1.0D+00 
          endif
        endif
c     adjust for ghost atoms
        if(izan.eq.0) ZREF(I) = 0
      END DO
      CALL DCOPY(NAT,ZREF,1,ZAN,1)
C
C     ----- CALCULATE NUMBER OF REFERENCE ELECTRONS -----
C
      NE = 0
      DO I = NAT0+1, NAT
        NE = NE + INT(ZREF(I) + 1.0D-02)
      END DO
C
      END SUBROUTINE DFTB_INPUT_FRAG
C
C-----------------------------------------------------------------------
C*MODULE DFTBFO    *DECK DFTB_PREVMUL
C>
C>    @brief Manage Mulliken populations
C>
C>    @details Save Mulliken charges after energy calculation to .F23
C>             file by calling with IMODE=1, and  restore Mulliken
C>             charges before monomer SCF calculations from .F23 to
C>             use them as initial guess by calling with IMODE=0
C>
C>    @author Yoshio Nishimoto
C>    - Sep, 2014- Subroutine written
C>    @date   Feb, 2016 - Yoshio Nishimoto
C>    - GDDI behavior may be fixed?
C>
C>           --- INPUT ---
C>    @param NFTDFTB Unit No. for storage of Mulliken population (F23)
C>    @param MATSIZ Size of POPMAT matrix
C>    @param IMODE (1) Save (2) Restore Mulliken charges
C>    @param POPMAT Matrix of Mulliken populations
C>    @param RUNTYP
C
      SUBROUTINE DFTB_PREVMUL(NFTDFTB,MATSIZ,IMODE,POPMAT,RUNTYP,mdout)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      INTEGER, INTENT(IN) :: NFTDFTB,MATSIZ,IMODE
      DOUBLE PRECISION, INTENT(INOUT) :: POPMAT(MATSIZ)
      CHARACTER(*), INTENT(IN) :: RUNTYP
      PARAMETER (MXUNIT=299)
      COMMON /FMCOM / X(1)
      COMMON /FMODTB/ LGAMMA,LZREF,LISPE,LSHIFTG,LSHIFTCT,
     *                LGAMMAD,MODESD,MODGAMMA,LDCTMP,LGAMMA3,LGAMMA3D,
     *                LNBOND,LCTMUL,LESPDFTB,LESPASC,NDFTBRST,LDFTBRST
      COMMON /IOFILE/ IR,IW,IP,IJK,IJKT,IDAF,NAV,IODA(950)
      COMMON /OPNNFT/ NFTOPN(MXUNIT),NODEXT(MXUNIT),IOSMP(2)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      LOGICAL :: OPT,GRAD
      INTEGER, SAVE :: NSERCH
      LOGICAL GOPARR,DSKWRK,MASWRK,DSKSAV,mdout
C
C     CALL WITH IMODE=0 BEFORE MONOSCF SUBROUTINE TO RESTORE PREVIOUS
C     FMO-DFTB CHARGES.
C     CALL WITH IMODE=1 AFTER ALL CALCULATIONS TO SAVE FMO-DFTB
C     CHARGES TO THE DICTIONATRY FILE
C
C     WRITE MULLIKEN CHARGES INTO $SCR/$JOB.F23 FILE
C
      OPT  = RUNTYP(1:8).EQ."OPTFMO  ".OR.RUNTYP(1:8).EQ."OPTIMIZE".OR.
     *       RUNTYP(1:8).EQ."SADPOINT".OR.RUNTYP(1:8).EQ."SURFACE ".OR.
     *       RUNTYP(1:8).EQ."MD      "
      GRAD = OPT.OR.RUNTYP(1:8).EQ."GRADIENT"
C
C     IN-MEMORY RESTART IF MEMORY IS ALLOCATED
C     NFTDFTB.EQ. 23 and 24 NEED NDFTBRST*2 and NDFTBRST MEMORY.
      MATID=NFTDFTB-22
      IF(OPT.AND.NDFTBRST*2.GE.MATSIZ.AND.(MATID.EQ.1.OR.MATID.EQ.2)) 
     *  THEN
        LDFTBRSTI=LDFTBRST+(MATID-1)*NDFTBRST*2
        IF(IMODE.EQ.1) CALL DCOPY(MATSIZ,POPMAT,1,X(LDFTBRSTI),1)
        IF(IMODE.EQ.0.AND.IAND(NSERCH,MATID).NE.0) 
     *    CALL DCOPY(MATSIZ,X(LDFTBRSTI),1,POPMAT,1)
        IF (IMODE.EQ.1) NSERCH = IOR(NSERCH,MATID)
        IF (MASWRK.AND.(IMODE.EQ.0.OR.MDOUT))
     *    WRITE (IW,*) 'DFTB RESTARTED DATA IN MEMORY FOR',IMODE,NFTDFTB
        RETURN
      ENDIF
C
      DSKSAV = DSKWRK
      DSKWRK = .TRUE.
      IF (GRAD) THEN
        IF (IMODE.EQ.1) THEN
          IF (OPT) THEN
C           IF (NSERCH.EQ.0) THEN
            IF (NFTOPN(NFTDFTB).EQ.0) THEN
           CALL SEQOPN(NFTDFTB,'JKFILE','UNKNOWN',.FALSE.,'UNFORMATTED')
              NSERCH = 0
            END IF
            CALL SEQREW(NFTDFTB)
            CALL SQWRIT(NFTDFTB,POPMAT,MATSIZ)
            IF (MASWRK.and.mdout) THEN
             WRITE (IW,'(" SAVE MULLIKEN IN .F",I2," FOR NSERCH =",I8)')
     *         NFTDFTB,NSERCH
            END IF
            NSERCH = NSERCH + 1
          END IF
        ELSE IF (IMODE.EQ.0) THEN
          IF (NSERCH.GT.0) THEN
            CALL SEQREW(NFTDFTB)
            CALL SQREAD(NFTDFTB,POPMAT,MATSIZ)
            IF (MASWRK) THEN
           WRITE (IW,'(" READ MULLIKEN FROM .F",I2," FOR NSERCH =",I8)')
     *         NFTDFTB,NSERCH
            END IF
c           write(6,*) (POPMAT(I),I=1,MATSIZ)
          END IF
        END IF
      END IF
      DSKWRK = DSKSAV
C
      END SUBROUTINE DFTB_PREVMUL
C
C-----------------------------------------------------------------------
C*MODULE DFTBFO    *DECK DFTB_ESP
C>
C>    @brief Calculate ESP in FMO-DFTB
C>
C>    @details Calculate electrostatic potential (ESP) in FMO-DFTB.
C>
C>    @author Yoshio Nishimoto
C>    - Sep, 2014- Subroutine written
C>    @date   Feb, 2018 - Yoshio Nishimoto
C>    - Added an option to omit doing for dimers/trimers
C>
C>           --- INPUT ---
C>    @param IAGLOB  Global atom index
C>    @param NATFRG  Number of atoms in each fragment
C>    @param INDFRG  First atom of each fragment in terms of NATFMOB?
C>    @param IATFRG  Global atom index with NBDFG?
C>    @param POPMAT  Matrix of Mulliken population
C>    @param SS      Overlap matrix
C>    @param FAO     Fock (Hamilton) matrix
C>    @param ISPEFMO ISPE (Index of species) for all NATFMO
C>    @param FMOC    Coordinates for all NATFMO
C>    @param GAMMA3  Gamma of DFTB3
C>    @param INDAT   Converts global atomic location to fragment index
C>    @param IALOC   Converts global to local (fragment) atomic location
C>    @param NN      See the comment below
C>    @param NN2     See the comment below
C
      SUBROUTINE DFTB_ESP(iaglob,natfrg,indfrg,iatfrg,popmat,
     *  ss,fao,ISPEFMO,FMOC,gamma3,INDAT,IALOC,NN,NN2,MODE)
      use mx_limits, only: mxatm
      USE DFTBPB_MOD,ONLY: PERIOD,PME,DFTB_GAMFMO_PBC,DFTB_ESP_PME
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      INTEGER,PARAMETER :: MXSPE=10
      LOGICAL DFTBFL,SCC,SRSCC,DFTB3,DAMPXH0,LCDFTB,dotrunces,dotruncesd
      DOUBLE PRECISION, PARAMETER :: ZERO=0.0D+00,HALF=5.0D-01,
     * ONE=1.0D+00,ONE_THIRD=ONE/3.0D+00
C
      COMMON /DFTB  / DFTBFL,SCC,SRSCC,DFTB3,DAMPXH0,LCDFTB
      COMMON /DFTBPR/ ETEMP,DFTBDP(MXSPE*14),DAMPXHE,HUBDER(MXSPE),
     *                ZREF(MXATM),SPNCST(6,MXSPE),SPE(MXATM),NSPE,
     *                MAXANG(MXATM),ISPE(MXATM),IND(MXATM+1),IDFTBD,
     *                PARAMDIR
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /IOFILE/ IR,IW,IP,IJK,IJKT,IDAF,NAV,IODA(950)
      COMMON /FMCOM / X(1)
      COMMON /FMOINF/ NFG,NLAYER,NATFMO,NBDFG,NAOTYP,NBODY,NSEGM
      COMMON /FMODTB/ LGAMMA,LZREF,LISPE,LSHIFTG,LSHIFTCT,
     *                LGAMMAD,MODESD,MODGAMMA,LDCTMP,LGAMMA3,LGAMMA3D,
     *                LNBOND,LCTMUL,LESPDFTB,LESPASC,NDFTBRST,LDFTBRST
      common /fmodim/ maxbnd,maxknd,maxcbs,maxcao,maxvec,maxl1,maxnat,
     *                maxabd,maxbas,maxbbd,maxlmo,maxslo,maxabd2,maxrij
      common /fmorun/ espscf,e0scf(2),emp2s,IDAFMO,icurfg,jcurfg,kcurfg,
     *                icurlay,icurunt,nat1e,ncursh,ngau,icurpop,ifmostp,
     *                moncor,needr,modrst,norbproj,nunesp,iskipesp,
     *                IESDPPC,idoprop,mp2run,icurit,idmfmo,iddfmo,
     *                iddcur,nddleft,ivmfmo,nzmtfmo,ifmobas,itmfmo(2)
      common /cenrad/ rcutes,lcradfg,lichfgg,dotrunces,dotruncesd
C
      double precision popmat(maxnat,nfg),ss(*),fao(*),cencur(3)
C     integer layfrg(*)
c     MODE=0 update FAO
c     MODE=1 do not update FAO
C
      INTEGER :: IAGLOB(*),NATFRG(*),INDFRG(*),IATFRG(*),ISPEFMO(*),
     *  INDAT(*),IALOC(*)
      !! NN/NN2 = 1/NATFMO  (IAND(MODGAM,1).NE.0)
      !! NN/NN2 = 2/NATFMOB (IAND(MODGAM,2).NE.0)
      DOUBLE PRECISION :: SHIFT(NAT),FMOC(3,NATFMO),TMP(12,NSPE,NSPE),
     * SHIFT2(NAT*2),TMP3(10,NSPE,NSPE),GAMMA3(NN2,NN2,*)
      DOUBLE PRECISION :: CHATMP(MAXNAT*2)
      LOGICAL GAMDIF,DAMPXH
C
C     ------ ADD ESP for FMO-DFTB corresponding to the
C     SECOND-ORDER HAMILTONIAN term in DFTB -----
C
      ifg=icurfg
      jfg=jcurfg
      kfg=kcurfg
C
      CALL VCLR(SHIFT,1,NAT)
      IF (DFTB3) CALL VCLR(SHIFT2,1,NAT*2)
C
      GAMDIF=.FALSE.
C     write (*,*) "modgamma = ", modgamma
      IF (IAND(MODGAMMA,8).NE.0) GAMDIF=.TRUE.
      IF (IAND(MODGAMMA,2).NE.0) GAMDIF=.TRUE.
      IF (PERIOD.AND.PME)        GAMDIF=.TRUE.
C     gamdif = gamdif.and.mode.eq.0
      NATFMOB = NATFMO + NBDFG
      IF ((GAMDIF.AND.(JFG.NE.0.OR.KFG.NE.0)).OR.(PERIOD.AND.PME)) THEN
       !IF (JFG.EQ.0.AND.KFG.EQ.0) THEN
       !  INDI = INDFRG(IFG)
       !  CALL DCOPY(NAT,X(LESPDFTB+INDI-1),1,SHIFT,1)
       !  CALL DSCAL(NAT,-ONE,SHIFT,1)
       !  IF (DFTB3) THEN
       !    CALL DCOPY(NAT,X(LESPDFTB+NATFMOB+INDI-1),1,SHIFT2,1)
       !    CALL DSCAL(NAT,-ONE,SHIFT2,1)
       !  END IF
       !ELSE IF (JFG.NE.0.OR.KFG.NE.0) THEN
       !IF (JFG.NE.0.OR.KFG.NE.0) THEN
C       IF (IAND(MODGAMMA,2).NE.2.OR.(JFG.NE.0.OR.KFG.NE.0)) THEN
        !! NO NEED TO ADD BDA CONTRIBUTIONS; IT IS JUST DOUBLE-COUNTING
          DO IAT = 1, NAT
            IAG = IAGLOB(IAT)
            IATFG = INDAT(IAG)
            IATLOC=IAND(IALOC(IAG),65535)
            INDI = INDFRG(IATFG)
            ILOC = INDI + IATLOC - 1
            SHIFT(IAT) = SHIFT(IAT) + X(LESPDFTB+ILOC-1)
            IF (DFTB3) THEN
              SHIFT2(IAT) = SHIFT2(IAT) + X(LESPDFTB+NATFMOB+ILOC-1)
            END IF
          END DO
          CALL DSCAL(NAT,-ONE,SHIFT,1)
          IF (DFTB3) CALL DSCAL(NAT,-ONE,SHIFT2,1)
       !END IF
       !END IF
      ELSE IF (GAMDIF.AND.IAND(MODGAMMA,2).NE.0) THEN
      ELSE
        GAMDIF=.FALSE.
      END IF
C
C     LOOP FOR ENVIRONMENTAL FRAGMENTS
C
      IF (IAND(MODGAMMA,3).EQ.0) THEN
        CALL DCOPY(12*NSPE*NSPE,X(LGAMMA),1,TMP,1)
        IF (DFTB3) CALL DCOPY(10*NSPE*NSPE,X(LGAMMA3),1,TMP3,1)
        NFG0 = NFG
        IF (GAMDIF) THEN
          IF (IFG.NE.0) NFG0 = 1
          IF (JFG.NE.0) NFG0 = 2
          IF (KFG.NE.0) NFG0 = 3
        END IF
        if(dotrunces) call radcent(cencur,radcur)
        DO LFG0 = 1, NFG0
          LFG = LFG0
          IF (GAMDIF) LFG=-1
          IF (LFG.NE.IFG.AND.LFG.NE.JFG.AND.LFG.NE.KFG) THEN
            IF (GAMDIF) THEN
              IF (LFG0.EQ.1) THEN
                LFG = IFG
              ELSE IF (LFG0.EQ.2) THEN
                LFG = JFG
              ELSE IF (LFG0.EQ.3) THEN
                LFG = KFG
              END IF
            END IF
            if(dotrunces.and.ich.eq.0) then
              lcradfgl=lcradfg+(lfg-1)*4
              lch=ixftch(x(lichfgg),lfg)  
              if(lch.eq.0.and.fmodest(cencur,radcur,x(lcradfgl),
     *            x(lcradfgl+3)).gt.rcutes) then
c               write(6,*) 'Skipping',ifg,jfg,lfg
                cycle
              endif
            endif
            NATJ = NATFRG(LFG)
            INDJ = INDFRG(LFG)
            JJ = INDJ - 1
            DO J = 1, NATJ !! ENVIRONMENTAL FRAGMENTS
              JJ = JJ + 1
              CHAMUL = POPMAT(J,LFG) - X(LZREF+JJ-1)
              IF (CHAMUL.EQ.ZERO) CYCLE
              JAG = IATFRG(JJ)
              JSP = ISPEFMO(JAG)
              DO I = 1, NAT !! FRAGMNET TO BE CONSIDERED
                IAG = IAGLOB(I)
                IATFG = INDAT(IAG)
                ISP = ISPEFMO(IAG)
                DAMPXH = DAMPXH0
                CALL DFTB_DMPCHK(ISP,JSP,DAMPXH)
                IF (PERIOD) THEN
                  CALL DFTB_GAMFMO_PBC(IAG,JAG,FMOC(1,IAG),FMOC(1,JAG),
     *                                 GAMMAV,GAMI,GAMJ,TMP(1,ISP,JSP),
     *                                 TMP3(1,ISP,JSP),TMP3(1,JSP,ISP),
     *                                 DAMPXH,ISP,JSP)
                ELSE
                  XX = FMOC(1,IAG)-FMOC(1,JAG)
                  YY = FMOC(2,IAG)-FMOC(2,JAG)
                  ZZ = FMOC(3,IAG)-FMOC(3,JAG)
                  DIST = XX*XX + YY*YY + ZZ*ZZ
                  CALL DFTB_GAM23(0,DIST,GAMMAV,GAMI,GAMJ,
     *                 TMP(1,ISP,JSP),TMP3(1,ISP,JSP),TMP3(1,JSP,ISP),
     *                 DAMPXH,DAMPXHE,SPE(ISP),SPE(JSP),HUBDER(ISP),
     *                 HUBDER(JSP),DFTB3)
                END IF
                SHIFT(I) = SHIFT(I) + CHAMUL * GAMMAV
                IF (DFTB3) THEN
                  SHIFT2(I) = SHIFT2(I)
     *              + CHAMUL * GAMI !! 2/3 coefficiets later multiplied
                  SHIFT2(I+NAT) = SHIFT2(I+NAT)
     *              + CHAMUL * CHAMUL * GAMJ
                END IF
              END DO
            END DO
          END IF
        END DO
      END IF
C
      IF (IAND(MODGAMMA,1).NE.0) THEN
        NFG0 = NFG
        IF (GAMDIF) THEN
          IF (IFG.NE.0) NFG0 = 1
          IF (JFG.NE.0) NFG0 = 2
          IF (KFG.NE.0) NFG0 = 3
        END IF
        DO LFG0 = 1, NFG0
          LFG = LFG0
          IF (GAMDIF) LFG=-1
          IF (LFG.NE.IFG.AND.LFG.NE.JFG.AND.LFG.NE.KFG) THEN
            IF (GAMDIF) THEN
              IF (LFG0.EQ.1) THEN
                LFG = IFG
              ELSE IF (LFG0.EQ.2) THEN
                LFG = JFG
              ELSE IF (LFG0.EQ.3) THEN
                LFG = KFG
              END IF
            END IF
            NATJ = NATFRG(LFG)
            INDJ = INDFRG(LFG)
            JJ = INDJ - 1
            DO J = 1, NATJ !! ENVIRONMENTAL FRAGMENTS
              JJ = JJ + 1
              CHAMUL = POPMAT(J,LFG) - X(LZREF+JJ-1)
              IF (CHAMUL.EQ.ZERO) CYCLE
              JAG = IATFRG(JJ)
              DO I = 1, NAT !! FRAGMNET TO BE CONSIDERED
                IAG = IAGLOB(I)
c               CALL DFTB_CNVSQ(IAG,JAG,NSEQ)
                IF (IAG.GE.JAG) THEN
                  NSEQ = IAG*(IAG-1)/2 + JAG
                ELSE
                  NSEQ = JAG*(JAG-1)/2 + IAG
                END IF
                GAMMAV = X(LGAMMA+NSEQ-1)
                SHIFT(I) = SHIFT(I) + CHAMUL * GAMMAV
                IF (DFTB3) THEN
                  GAMI = GAMMA3(IAG,JAG,1)
                  GAMJ = GAMMA3(JAG,IAG,1)
                  SHIFT2(I) = SHIFT2(I)
     *              + CHAMUL * GAMI !! 2/3 coefficiets later multiplied
                  SHIFT2(I+NAT) = SHIFT2(I+NAT)
     *              + CHAMUL * CHAMUL * GAMJ
                END IF
              END DO
            END DO
          END IF
        END DO
      END IF
C
      IF (IAND(MODGAMMA,2).NE.0) THEN
        NATFMOB = NATFMO + NBDFG
        NFG0 = NFG
        IF (KFG.NE.0) THEN
          NFG0 = 3
        ELSE IF (JFG.NE.0) THEN
          NFG0 = 2
        ELSE
          NFG0 = 1
          NATI = NATFRG(IFG)
          INDI = INDFRG(IFG)
          CALL DGEMV('T',NN2,NATI,1.0D+00,X(LGAMMA+NN2*(INDI-1)),NN2,
     *               X(LSHIFTG),1,0.0D+00,SHIFT,1)
          CALL DSCAL(NATI,-ONE,SHIFT,1)
          IF (DFTB3) THEN
            CALL DGEMV('T',NN2,NATI,1.0D+00,GAMMA3(1,INDI,2),NN2,
     *                 X(LSHIFTG),1,0.0D+00,SHIFT2,1)
            CALL DGEMV('T',NN2,NATI,1.0D+00,GAMMA3(1,INDI,1),NN2,
     *                 X(LSHIFTG+NN2),1,0.0D+00,SHIFT2(NATI+1),1)
            CALL DSCAL(2*NATI,-ONE,SHIFT2,1)
          END IF
        END IF
        DO LFG0 = 1, NFG0
          IF (LFG0.EQ.1) THEN
            LFG = IFG
          ELSE IF (LFG0.EQ.2) THEN
            LFG = JFG
          ELSE IF (LFG0.EQ.3) THEN
            LFG = KFG
          END IF
          NATJ = NATFRG(LFG)
          INDJ = INDFRG(LFG)
          DO I = 1, NAT !! FRAGMNET TO BE CONSIDERED
            IAG = IAGLOB(I)
            IATFG = INDAT(IAG)
            IATLOC = IAND(IALOC(IAG),65535)
            INDI = INDFRG(IATFG)
            NSEQ = INDJ-1 + NN2*(INDI+IATLOC-2)
            SHIFT(I) = SHIFT(I)  
     *        + DDOT(NATJ,X(LSHIFTG+INDJ-1),1,X(LGAMMA+NSEQ),1)
            IF (DFTB3) THEN
              II = INDI+IATLOC-1
              SHIFT2(I) = SHIFT2(I)
     *          + DDOT(NATJ,X(LSHIFTG+INDJ-1),1,GAMMA3(INDJ,II,2),1)
              SHIFT2(I+NAT) = SHIFT2(I+NAT)
     *          + DDOT(NATJ,X(LSHIFTG+NN2+INDJ-1),1,GAMMA3(INDJ,II,1),1)
            END IF
          END DO
        END DO
      END IF
C
      IF (PERIOD.AND.PME.AND.MODE.EQ.0) THEN
        CALL DFTB_ESP_PME(IFG,JFG,KFG,NAT,POPMAT,SHIFT)
      END IF
C
      IF (GAMDIF) THEN
!       IF (JFG.EQ.0.AND.KFG.EQ.0) THEN
!         NATI = NATFRG(IFG)
!         INDI = INDFRG(IFG)
!         CALL DAXPY(NATI,ONE,SHIFT,1,X(LESPDFTB+INDI-1),1)
!         CALL DCOPY(NATI,SHIFT,1,X(LESPDFTB+INDI-1),1)
!       ELSE
          CALL DSCAL(NAT,-ONE,SHIFT,1)
          IF (DFTB3) CALL DSCAL(NAT*2,-ONE,SHIFT2,1)
!       END IF
      END IF
C
      IF (IAND(MODGAMMA,4096).NE.0) THEN
        CALL VCLR(SHIFT,1,NAT)
        IF (DFTB3) CALL VCLR(SHIFT2,1,NAT*2)
      END IF
C
C     ----- SAVE THE ARRAY FOR GRADIENT CALCULATION -----
C
C     SAVE THE SHIFT ARRAY.  THIS WILL BE USED TO CALCULATE THE NCC
C     ENERGY IN DFTBX SUBROUTINE TO SUBTRACT THE SCC CONTRIBUTION.
C     SHIFT MATRIX IS THEN USED TO CALCULATE ESP ENERGY DURING
C     SCF CALCULATIONS.
C
      IF (DFTB3) THEN
        !! COEFFICIENTS ARE MULTIPLIED LATER IN dftbx
C       DO I = 1, NAT*2
C         SHIFT2(I) = SHIFT2(I) * ONE_THIRD
C       END DO
        CALL DSCAL(NAT*2,ONE_THIRD,SHIFT2,1)
        CALL DAWRIT(IDAF,IODA,SHIFT2,NAT*2,567,0)
      END IF
      CALL DAWRIT(IDAF,IODA,SHIFT,NAT,566,0)
C
      IF(MODE.EQ.0) THEN
      DO I = 1, NAT
        SHIFTI = SHIFT(I)
        DO J = 1, I
          SHIFTJ = SHIFT(J)
          SHIFTIJ = (SHIFTI+SHIFTJ)*HALF
          DO K = 1, MAXANG(ISPE(I))**2
            M = IND(I)+K
            DO L = 1, MAXANG(ISPE(J))**2
              N = IND(J)+L
              CALL DFTB_CNVSQ(M,N,NSEQ)
              FAO(NSEQ)=FAO(NSEQ) + SS(NSEQ)*SHIFTIJ
              IF (I.EQ.J .AND. K.EQ.L) EXIT
            END DO
          END DO
        END DO
      END DO
      ENDIF
C
      RETURN
C
      END SUBROUTINE DFTB_ESP
C
C-----------------------------------------------------------------------
C*MODULE DFTBFO    *DECK DFTB_ESP3
C>
C>    @brief Calculate ESP in FMO-DFTB
C>
C>    @details Calculate electrostatic potential (ESP) in FMO-DFTB.
C>
C>    @author Yoshio Nishimoto
C>    - Xxx, 2020- Subroutine written
C>
C>           --- INPUT ---
C>    @param NATFRG  Number of atoms in each fragment
C>    @param INDFRG  First atom of each fragment in terms of NATFMOB?
C>    @param IATFRG  Global atom index with NBDFG?
C>    @param POPMAT  Matrix of Mulliken population
C>    @param ISPEFMO ISPE (Index of species) for all NATFMO
C>    @param FMOC    Coordinates for all NATFMO
C
      SUBROUTINE DFTB_ESP3(NATFMOB,NATFRG,INDFRG,IATFRG,POPMAT,
     *                     ISPEFMO,FMOC,ZREFFMO,TMP,TMP3,SHIFT)
      use mx_limits, only: mxatm
      USE DFTBPB_MOD,ONLY: PERIOD,DFTB_GAMFMO_PBC
C
      IMPLICIT NONE
C     Declarations of arguments
      DOUBLE PRECISION POPMAT(maxnat,nfg)
      DOUBLE PRECISION FMOC(3,NATFMO),ZREFFMO(*),TMP(12,NSPE,NSPE),
     *                 TMP3(10,NSPE,NSPE),SHIFT(NATFMOB,2)
      INTEGER NATFMOB,NATFRG(*),INDFRG(*),IATFRG(*)
      INTEGER ISPEFMO(*)
C     Declarations of common blocks
      LOGICAL DFTBFL,SCC,SRSCC,DFTB3,DAMPXH0,LCDFTB
      DOUBLE PRECISION ETEMP,DFTBDP,DAMPXHE,HUBDER
      DOUBLE PRECISION ZREF,SPNCST,SPE
      DOUBLE PRECISION PARAMDIR
      INTEGER NSPE
      INTEGER MAXANG,ISPE,IND,IDFTBD
      DOUBLE PRECISION X
      INTEGER NFG,NLAYER,NATFMO,NBDFG,NAOTYP,NBODY,NSEGM
      INTEGER maxbnd,maxknd,maxcbs,maxcao,maxvec,maxl1,maxnat
      INTEGER maxabd,maxbas,maxbbd,maxlmo,maxslo,maxabd2,maxrij
      INTEGER ISCOPE,NGROUPS,MYGROUP,MEGLOB,npglob,nnglob,JBTYP
      INTEGER NSUBGR
      LOGICAL ISGDDI,PAROUT,INITGDDI,wasgddi,MLGDDI
      INTEGER MeUniv,NPUniv,numdlb,myworld,nworlds,mogddi
C     Other declarations
      INTEGER, PARAMETER :: MXSPE=10
      DOUBLE PRECISION, PARAMETER :: ONE=1.0D+00,ONE_THIRD=ONE/3.0D+00
      DOUBLE PRECISION TWO_THIRD
      INTEGER IFG,LGROUP
      LOGICAL MYJOB
      INTEGER NATI,INDI,JFG,NATJ,INDJ,IAT,IAG,ISP
      DOUBLE PRECISION CHAMULI
      INTEGER JATMAX,JAT,JAG,JSP
      LOGICAL DAMPXH
      DOUBLE PRECISION CHAMULJ,GAMMAV,GAMI,GAMJ,XX,YY,ZZ,DIST,DISTI
      DOUBLE PRECISION DFTB_GAM2,VALI,VALJ
C
      COMMON /DFTB  / DFTBFL,SCC,SRSCC,DFTB3,DAMPXH0,LCDFTB
      COMMON /DFTBPR/ ETEMP,DFTBDP(MXSPE*14),DAMPXHE,HUBDER(MXSPE),
     *                ZREF(MXATM),SPNCST(6,MXSPE),SPE(MXATM),NSPE,
     *                MAXANG(MXATM),ISPE(MXATM),IND(MXATM+1),IDFTBD,
     *                PARAMDIR
      COMMON /FMCOM / X(1)
      COMMON /FMOINF/ NFG,NLAYER,NATFMO,NBDFG,NAOTYP,NBODY,NSEGM
      common /fmodim/ maxbnd,maxknd,maxcbs,maxcao,maxvec,maxl1,maxnat,
     *                maxabd,maxbas,maxbbd,maxlmo,maxslo,maxabd2,maxrij
      COMMON /GDDI/   ISCOPE,NGROUPS,MYGROUP,MEGLOB,npglob,nnglob,JBTYP,
     *                ISGDDI,PAROUT,INITGDDI,wasgddi,MLGDDI,NSUBGR,
     *                MeUniv,NPUniv,numdlb,myworld,nworlds,mogddi
C
      TWO_THIRD = 2.0D+00*ONE_THIRD
      DO IFG = 1, NFG
        IF (ISGDDI) THEN
          CALL GDDICOUNT(0,LGROUP,MYJOB)
          IF (.NOT.MYJOB) CYCLE
        END IF
        NATI = NATFRG(IFG)
        INDI = INDFRG(IFG)
        DO JFG = 1, IFG
          NATJ = NATFRG(JFG)
          INDJ = INDFRG(JFG)
          DO IAT = 1, NATI
            IAG = IATFRG(INDI+IAT-1)
            ISP = ISPEFMO(IAG)
            CHAMULI = POPMAT(IAT,IFG) - ZREFFMO(INDI+IAT-1)
            JATMAX = NATJ
            IF (IFG.EQ.JFG) JATMAX = IAT
            DO JAT = 1, JATMAX
              JAG = IATFRG(INDJ+JAT-1)
              JSP = ISPEFMO(JAG)
              DAMPXH = DAMPXH0
              CALL DFTB_DMPCHK(ISP,JSP,DAMPXH)
              CHAMULJ = POPMAT(JAT,JFG) - ZREFFMO(INDJ+JAT-1)
              IF (PERIOD) THEN
              CALL DFTB_GAMFMO_PBC(IAG,JAG,FMOC(1,IAG),FMOC(1,JAG),
     *                             GAMMAV,GAMI,GAMJ,TMP(1,ISP,JSP),
     *                             TMP3(1,ISP,JSP),TMP3(1,JSP,ISP),
     *                             DAMPXH,ISP,JSP)
              ELSE
                XX = FMOC(1,IAG)-FMOC(1,JAG)
                YY = FMOC(2,IAG)-FMOC(2,JAG)
                ZZ = FMOC(3,IAG)-FMOC(3,JAG)
                DIST = XX*XX + YY*YY + ZZ*ZZ
                IF (DFTB3) THEN
                  CALL DFTB_GAM23(0,DIST,GAMMAV,GAMI,GAMJ,
     *                 TMP(1,ISP,JSP),TMP3(1,ISP,JSP),TMP3(1,JSP,ISP),
     *                 DAMPXH,DAMPXHE,SPE(ISP),SPE(JSP),HUBDER(ISP),
     *                 HUBDER(JSP),DFTB3)
                ELSE
                  IF (DIST.GE.1.0D-15) THEN
                    DIST = SQRT(DIST)
                    DISTI = ONE/DIST
                  ELSE
                    DISTI = 0
                  END IF
                  GAMMAV = DFTB_GAM2(DIST,DISTI,TMP(1,ISP,JSP),DAMPXH,
     *                               DAMPXHE,SPE(ISP),SPE(JSP))
                END IF
              END IF
              VALI = CHAMULJ*GAMMAV
              VALJ = CHAMULI*GAMMAV
              IF (DFTB3) THEN
                VALI = VALI + ONE_THIRD*CHAMULJ*CHAMULJ*GAMJ
                VALJ = VALJ + ONE_THIRD*CHAMULI*CHAMULI*GAMI
              END IF
              IF (IFG.EQ.JFG.AND.IAT.EQ.JAT) THEN
                VALI = VALI*0.5D+00
                VALJ = VALJ*0.5D+00
              END IF
              SHIFT(INDI+IAT-1,1) = SHIFT(INDI+IAT-1,1) + VALI
              SHIFT(INDJ+JAT-1,1) = SHIFT(INDJ+JAT-1,1) + VALJ
              IF (DFTB3) THEN
                if (ifg.eq.jfg.and.iat.eq.jat) then
                  gami = gami*0.5d+00
                  gamj = gamj*0.5d+00
                end if
                SHIFT(INDI+IAT-1,2) = SHIFT(INDI+IAT-1,2) + CHAMULJ*GAMI
                SHIFT(INDJ+JAT-1,2) = SHIFT(INDJ+JAT-1,2) + CHAMULI*GAMJ
              END IF
            END DO
          END DO
        END DO
      END DO
C
C
      RETURN
C
      END SUBROUTINE DFTB_ESP3
C
C-----------------------------------------------------------------------
C*MODULE DFTBFO    *DECK DFTB_GAMFMO
C>
C>    @brief Calculate gamma values in FMO-DFTB
C>
C>    @details Calculate gamma (DFTB2 and DFTB3) values in FMO-DFTB.
C>             IF MODGAMMA.EQ.0, ONLY "TMP" ARRAY IS CALCULATED.  THE
C>             ARRAY WILL BE USED WHEN GAMMA VALUES ARE CALCULATED.
C>             THIS IMPLMENTATION IS FOR THE PROBLEM OF MEMORY.
C>             IF MODGAMMA.EQ.1, ALL GAMMA VALUES ARE CALCULATED IN THIS
C>             SUBROUTINE.  THE PROBLEM IS THAT TOO MUCH MEMORY IS USED
C>             FOR TRULY LARGE MOLECULES.
C>             ALSO, IF NDER=1 (GRADIENT), "TMPD" MATRIX IS CALCULATED.
C>             MODGAMMA=4096: Ignore ESP; only for development purpose
C>
C>    @author Yoshio Nishimoto
C>    - Sep, 2014- Subroutine written
C>    @date   Feb, 2016 - Yoshio Nishimoto
C>    - Added GDDI synchronization, which is rarely used
C>
C>           --- INPUT ---
C>    @param MODGAMMA Which values are stored (see above)
C>    @param NDER     Order of derivatives to be calculated
C>    @param LAYFRG   Layer of each fragment
C>    @param INDFRG   First atom of each fragment in terms of NATFMOB?
C>    @param IANFRG   Global atom index
C>    @param IATFRG   Global atom index with NBDFG?
C>    @param NATFRG   Number of atoms in each fragment
C>    @param INDFRG   First atom of each fragment in terms of NATFMOB?
C>    @param FMOZAN   Nuclear charges of all atoms
C>    @param ZANFRG   Nuclear charges in terms of NATFMOB?
C>    @param FMOC     Coordinates for all NATFMO
C>           --- OUTPUT ---
C>    @param ISPEFMO  ISPE (Index of species) for all NATFMO
C>    @param GAMMAFMO Array of (MODGAMMA=0) TMP (MODGAMMA=1) DFTB2 gamma
C>    @param GAMMAFMOD Array of TMPD for DFTB2
C>    @param ZREFMO   DFTB nuclear charges in terms of NATFMOB?
C>    @param GAMMA3   Array of (MODGAMMA=0,1) TMP (MODGAMMA=2) DFTB3 gamma
C>    @param GAMMA3D  Array of TMPD for DFTB3
C>    @param NBOND    Number of bonds for SKHP
C
      SUBROUTINE DFTB_GAMFMO(MODGAMMA,NDER,INDFRG,IANFRG,IATFRG,
     *  NATFRG,ISPEFMO,FMOZAN,ZANFRG,GAMMAFMO,GAMMAFMOD,ZREFFMO,FMOC,
     *  GAMMA3,GAMMA3D,NBOND,MDOUT)
      use mx_limits, only: mxatm
      USE DFTBPB_MOD,ONLY: PERIOD,DFTB_GAMFMO_PBC,PBCBOX
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      INTEGER,PARAMETER :: MXSPE=10
      LOGICAL MDOUT
      LOGICAL DFTBFL,SCC,SRSCC,DFTB3,DAMPXH0,LCDFTB
      LOGICAL GOPARR,DSKWRK,MASWRK
      LOGICAL MYJOB,ISGDDI,PAROUT,INITGDDI,WASGDDI,MLGDDI
C
      INTEGER,PARAMETER :: DDI_WORLD=0,DDI_GROUP=1
C
      COMMON /DFTB  / DFTBFL,SCC,SRSCC,DFTB3,DAMPXH0,LCDFTB
      COMMON /DFTBPR/ ETEMP,DFTBDP(MXSPE*14),DAMPXHE,HUBDER(MXSPE),
     *                ZREF(MXATM),SPNCST(6,MXSPE),SPE(MXATM),NSPE,
     *                MAXANG(MXATM),ISPE(MXATM),IND(MXATM+1),IDFTBD,
     *                PARAMDIR
      COMMON /DFTBSK/ SKDIM(MXSPE,MXSPE),SKSPIN(MXSPE),QREFL(3,MXSPE),
     *                HUBBL(3,MXSPE),QREF(MXSPE),HUBB(MXSPE),SKCUT2,
     *                NEEDSK,LSKHTAB(MXSPE,MXSPE),LSKSTAB(MXSPE,MXSPE),
     *                LSKGRID(MXSPE,MXSPE),LSKSELF(MXSPE)
      COMMON /IOFILE/ IR,IW,IP,IJK,IJKT,IDAF,NAV,IODA(950)
      COMMON /FMOINF/ NFG,NLAYER,NATFMO,NBDFG,NAOTYP,NBODY,NSEGM
      COMMON /GDDI/   ISCOPE,NGROUPS,MYGROUP,MEGLOB,npglob,nnglob,JBTYP,
     *                ISGDDI,PAROUT,INITGDDI,wasgddi,MLGDDI,NSUBGR,
     *                MeUniv,NPUniv,numdlb,myworld,nworlds,mogddi
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /TIMING/ CPU,WALL
C
      INTEGER INDFRG(*),IANFRG(*),IATFRG(*),
     * NATFRG(*),ISPEFMO(*),NBOND(*)
      DOUBLE PRECISION GAMMAFMO(*),ZREFFMO(*),FMOC(3,*),FMOZAN(*),
     * ZANFRG(*),ZREFTMP,TMP(12,NSPE,NSPE),TMPD(10,NSPE,NSPE),
     * GAMMAFMOD(*),GAMMA3(NATFMO,NATFMO),GAMMA3D(NATFMO,NATFMO),
     * TMP3(10,NSPE,NSPE),TMPD3(11,NSPE,NSPE)
      DOUBLE PRECISION VEC(3)
C     CHARACTER(4) :: ATOMTMP
      CHARACTER(8) :: C8DUM = '        '!, ATOMNM
      DOUBLE PRECISION C8/8H        /
C
      DOUBLE PRECISION, PARAMETER :: ZERO=0.0D+00,HALF=0.5D+00,
     * ONE=1.0D+00,TWO=2.0D+00,THREE=3.0D+00,ONE_48=ONE/4.8D+01
      DOUBLE PRECISION, PARAMETER :: HUBTOL = 1.0D-04
      LOGICAL :: TSQ,DDIGAM,DAMPXH
C
      ATOMNM = C8
      DO I = 1, NATFMO
        IMODE = INT(FMOZAN(I))
        CALL ZSYMNUM(C8DUM,ATOMNM,IMODE)
C       CALL ZSYMNUM(C8DUM,ATOMTMP,IMODE)
C       WRITE (UNIT=ATOMNM,FMT='(A4,A4)') ATOMTMP,'    '
        DO J = 1, NSPE
          IF (ATOMNM==SPE(J)) THEN
            ISPEFMO(I) = J
            EXIT
          END IF
        END DO
      END DO
C
      DO IFG = 1, NFG
        INDI = INDFRG(IFG)
        NATI = NATFRG(IFG)
        DO I = 1, NATI
          IAG = IATFRG(INDI+I-1)
          ZREFTMP = QREF(ISPEFMO(IAG))
          ZANTMP = ZANFRG(INDI+I-1)
          izan = int(ZANTMP+1.0D-02)
          IANTMP = IANFRG(INDI+I-1)
          if (IANTMP.ne.izan) then
            if (izan.eq.1) then
              ZREFTMP=1.0D+00 
            else
              ZREFTMP=ZREFTMP-1.0D+00 
            endif
          endif
          ZREFFMO(INDI+I-1) = ZREFTMP
        END DO
      END DO
C
      IF (.NOT.SCC) GO TO 104
C
C     CALL CPU_TIME(T0)
      DO ISP = 1, NSPE
        DO JSP = ISP, NSPE
          TAUA    = 3.2D+00 * HUBBL(1,ISP)
          TAUB    = 3.2D+00 * HUBBL(1,JSP)
          TAUMEAN = HALF * (TAUA + TAUB)
          TAUA2   = TAUA**2
          TAUB2   = TAUB**2
          TAUA4   = TAUA2**2
          TAUB4   = TAUB2**2
          TAUA6   = TAUA2*TAUA4
          TAUB6   = TAUB2*TAUB4
          IF (ISP.EQ.JSP) THEN
            TAUAB2I = ZERO
            TAUAB3I = ZERO
          ELSE
            TAUAB2I = HALF / (TAUA2-TAUB2)**2
            TAUAB3I = ONE / (TAUA2-TAUB2)**3
          END IF
          TMP( 1,ISP,JSP) = TAUA
          TMP( 2,ISP,JSP) = TAUB
          TMP( 3,ISP,JSP) = TAUMEAN
          TMP( 4,ISP,JSP) = 0.6875D+00*TAUMEAN
          TMP( 5,ISP,JSP) = 0.1875D+00*(TAUMEAN**2)
          TMP( 6,ISP,JSP) = 0.0208333333333D+00*(TAUMEAN**3)
          TMP( 7,ISP,JSP) = TAUB4*TAUA*TAUAB2I
          TMP( 8,ISP,JSP) = TAUA4*TAUB*TAUAB2I
          TMP( 9,ISP,JSP) = TAUB6-THREE*TAUB4*TAUA2
          TMP(10,ISP,JSP) = TAUA6-THREE*TAUA4*TAUB2
          TMP(11,ISP,JSP) = TAUAB3I
          TMP(12,ISP,JSP) = 0.15625D+00*(TAUA+TAUB)
          IF (ISP.NE.JSP) THEN
            CALL DCOPY(12,TMP(1,ISP,JSP),1,TMP(1,JSP,ISP),1)
            TMP( 1,JSP,ISP) =  TMP( 2,ISP,JSP)
            TMP( 2,JSP,ISP) =  TMP( 1,ISP,JSP)
            TMP( 7,JSP,ISP) =  TMP( 8,ISP,JSP)
            TMP( 8,JSP,ISP) =  TMP( 7,ISP,JSP)
            TMP( 9,JSP,ISP) =  TMP(10,ISP,JSP)
            TMP(10,JSP,ISP) =  TMP( 9,ISP,JSP)
            TMP(11,JSP,ISP) = -TMP(11,ISP,JSP)
          END IF
        END DO
      END DO
C
      IF (DFTB3) THEN
        DO I = 1, NSPE
          TAUA = 3.2D+00*HUBB(I)
          TAUA2 = TAUA**2
          TAUA3 = TAUA*TAUA2
          TAUA4 = TAUA2**2
          TAUA6 = TAUA2*TAUA4
          DO J = 1, NSPE
            TAUB = 3.2D+00*HUBB(J)
            TAUB2 = TAUB**2
            TAUB3 = TAUB*TAUB2
            TAUB4 = TAUB2**2
            TAUB6 = TAUB2*TAUB4
            IF (ABS(HUBB(I)-HUBB(J)).LE.HUBTOL) THEN
              TMP3( 1,I,J) = TAUA
              TMP3( 2,I,J) = TAUB
              TMP3( 3,I,J) = -3.125D-01
              TMP3( 4,I,J) = -3.125D-01*TAUA
              TMP3( 5,I,J) = -1.25D-01*TAUA2
              TMP3( 6,I,J) = -ONE_48*TAUA3
              TMP3( 7,I,J) = 6.875D-01*TAUA
              TMP3( 8,I,J) = 1.875D-01*TAUA2
              TMP3( 9,I,J) = ONE_48*TAUA3
              TMP3(10,I,J) = ZERO
            ELSE
              TAUAB = TAUA2 - TAUB2
              TAUABI = ONE/TAUAB
              TAUABI2 = TAUABI**2
              TAUABI3 = TAUABI**3
              TAUABI4 = TAUABI**4
              TMP3( 1,I,J) = TAUA
              TMP3( 2,I,J) = TAUB
              TMP3( 3,I,J) = -HALF*(TAUB6+THREE*TAUA2*TAUB4)*TAUABI3
              TMP3( 4,I,J) = -1.2D+01*TAUB4*TAUA3*TAUABI4
              TMP3( 5,I,J) = HALF*TAUA*TAUB4*TAUABI2
              TMP3( 6,I,J) = -(TAUB6-THREE*TAUA2*TAUB4)*TAUABI3
              TMP3( 7,I,J) = -TWO*TAUA3*TAUB3*TAUABI3
              TMP3( 8,I,J) = 1.2D+01*TAUB4*TAUA3*TAUABI4
              TMP3( 9,I,J) = HALF*TAUB*TAUA4*TAUABI2
              TMP3(10,I,J) = (TAUA6-THREE*TAUB2*TAUA4)*TAUABI3
            END IF
          END DO
        END DO
      END IF
C
C     CALCULATE "TMPD" ARRAY IF GRADIENT CALCULATION IS SPECIFIED.
C
      IF (NDER.GT.0) THEN
        DO ISP = 1, NSPE
          DO JSP = ISP, NSPE
            TAUA     = 3.2D+00 * HUBBL(1,ISP)
            TAUB     = 3.2D+00 * HUBBL(1,JSP)
            TAUMEAN  = HALF * (TAUA + TAUB)
            TAUMEAN2 = TAUMEAN*TAUMEAN
            TAUMEAN3 = TAUMEAN*TAUMEAN2
            TAUA2   = TAUA**2
            TAUB2   = TAUB**2
            TAUA4   = TAUA2**2
            TAUB4   = TAUB2**2
            TAUA6   = TAUA2*TAUA4
            TAUB6   = TAUB2*TAUB4
            IF (ISP.EQ.JSP) THEN
              TAUAB2I = ZERO
              TAUAB3I = ZERO
            ELSE
              TAUAB2I = HALF / (TAUA2-TAUB2)**2
              TAUAB3I = ONE / (TAUA2-TAUB2)**3
            END IF
            TMPD( 1,ISP,JSP) = TAUA
            TMPD( 2,ISP,JSP) = TAUB
            TMPD( 3,ISP,JSP) = TAUMEAN
            TMPD( 4,ISP,JSP) = 0.6875D+00*TAUMEAN
            TMPD( 5,ISP,JSP) = 0.1875D+00*TAUMEAN2
            TMPD( 6,ISP,JSP) = 0.0208333333333D+00*TAUMEAN3
            TMPD( 7,ISP,JSP) = (TAUB6-THREE*TAUA2*TAUB4)*TAUAB3I
            TMPD( 8,ISP,JSP) = -(TAUA6-THREE*TAUB2*TAUA4)*TAUAB3I
            TMPD( 9,ISP,JSP) = TAUB4*TAUA*TAUAB2I
            TMPD(10,ISP,JSP) = TAUA4*TAUB*TAUAB2I
            IF (ISP.NE.JSP) THEN
              CALL DCOPY(10,TMPD(1,ISP,JSP),1,TMPD(1,JSP,ISP),1)
              TMPD( 1,JSP,ISP) = TMPD( 2,ISP,JSP)
              TMPD( 2,JSP,ISP) = TMPD( 1,ISP,JSP)
              TMPD( 7,JSP,ISP) = TMPD( 8,ISP,JSP)
              TMPD( 8,JSP,ISP) = TMPD( 7,ISP,JSP)
              TMPD( 9,JSP,ISP) = TMPD(10,ISP,JSP)
              TMPD(10,JSP,ISP) = TMPD( 9,ISP,JSP)
            END IF
          END DO
        END DO
        CALL DCOPY(10*NSPE*NSPE,TMPD,1,GAMMAFMOD,1)
C       WRITE (*,*) GAMMAFMOD
C
        IF (DFTB3) THEN
          TWELVE = 1.2D+01
          ONE_15 = ONE/1.5D+01
          ONE_24 = ONE/2.4D+01
          DAMPXH = DAMPXH0
          DO ISP = 1, NSPE
C           HUBA = HUBBL(1,ISP)
C           TAUA = 3.2D+00*HUBA
            TAUA = 3.2D+00*HUBBL(1,ISP)
            TAUA2 = TAUA**2
            TAUA3 = TAUA*TAUA2
            TAUA4 = TAUA2**2
            TAUA6 = TAUA2*TAUA4
            DO JSP = 1, NSPE
C             HUBB = HUBBL(1,JSP)
C             TAUB = 3.2D+00*HUBB
              TAUB = 3.2D+00*HUBBL(1,JSP)
              TAUB2 = TAUB**2
              TAUB3 = TAUB*TAUB2
              TAUB4 = TAUB2**2
              TAUB6 = TAUB2*TAUB4
              IF (ABS(HUBBL(1,ISP)-HUBBL(1,JSP)).LE.HUBTOL) THEN
                TMPD3(1,ISP,JSP) = 1.875D-01*TAUA2       !! 3/16 a^2
                TMPD3(2,ISP,JSP) = ONE_24*TAUA3          !! 1/24 a^3
                TMPD3(3,ISP,JSP) = 6.875D-01*TAUA        !! 11/16 a
                TMPD3(4,ISP,JSP) = HALF*TMPD3(2,ISP,JSP) !! 1/48 a^3
                TMPD3(5,ISP,JSP) = -6.0D-01*TAUA
                TMPD3(6,ISP,JSP) = 2.0D-01*TAUA2
                TMPD3(7,ISP,JSP) = 2.0D-01*TAUA3
                TMPD3(8,ISP,JSP) = TAUA4*ONE_15
                IF (DAMPXH) THEN
                  TMPD3(9,ISP,JSP)  = -TAUA
                  TMPD3(10,ISP,JSP) = -4.0D-01*TAUA2
                  TMPD3(11,ISP,JSP) = -TAUA3*ONE_15
                END IF
              ELSE
                TAUAB = TAUA2 - TAUB2
                TAUABI = ONE/TAUAB
                TAUABI2 = TAUABI**2
                TAUABI3 = TAUABI**3
                TAUABI4 = TAUABI**4
                TMPD3(1,ISP,JSP) =  (TAUB6-THREE*TAUA2*TAUB4)*TAUABI3
                TMPD3(2,ISP,JSP) = HALF*TAUA*TAUB4*TAUABI2
                TMPD3(3,ISP,JSP) = -(TAUA6-THREE*TAUB2*TAUA4)*TAUABI3
                TMPD3(4,ISP,JSP) = HALF*TAUB*TAUA4*TAUABI2
                TMPD3(5,ISP,JSP)=-HALF*(TAUB6+THREE*TAUA2*TAUB4)*TAUABI3
                TMPD3(6,ISP,JSP) = -TWELVE*TAUA3*TAUB4*TAUABI4
                TMPD3(7,ISP,JSP) = -TWO*TAUB3*TAUA3*TAUABI3
                TMPD3(8,ISP,JSP) =  TWELVE*TAUB4*TAUA3*TAUABI4
                TMPD3(9,ISP,JSP) =  ONE/((TAUA+TAUB)*1.5625D-01)
              END IF
            END DO
          END DO
          CALL DCOPY(11*NSPE*NSPE,TMPD3,1,GAMMA3D,1)
        END IF
      END IF
C
      TSQ  = IAND(MODGAMMA,2).NE.0
C
      NATFMOB = NATFMO+NBDFG
      DDIGAM = .FALSE.
      IF (MASWRK.AND.MDOUT)
     *  WRITE (IW,'(" PREPARE CONSTANTS FOR GAMMA CALCULATION")')
      IF (IAND(MODGAMMA,3).EQ.0) THEN
        CALL DCOPY(12*NSPE*NSPE,TMP,1,GAMMAFMO,1)
        IF (DFTB3) CALL DCOPY(10*NSPE*NSPE,TMP3,1,GAMMA3,1)
      ELSE
        IF (MASWRK.AND.MDOUT)
     *    WRITE (IW,'(" ALL GAMMA VALUES ARE CALCULATED HERE")')
        IF (ISGDDI.AND.IAND(MODGAMMA,4).NE.0) THEN
          DDIGAM=.TRUE.
          CALL GDDI_SCOPE(DDI_GROUP)
          CALL GDDICOUNT(-1,LGROUP,MYJOB)
          IF (MASWRK.AND.MDOUT) WRITE (IW,'(" USE GDDI")')
        END IF
        IF (TSQ) THEN
          DO IFG = 1, NFG
            IF (DDIGAM) THEN
              CALL GDDICOUNT(0,LGROUP,MYJOB)
              IF (.NOT.MYJOB) CYCLE
            END IF
            NATI = NATFRG(IFG)
            INDI = INDFRG(IFG)
            DO JFG = 1, IFG
              NATJ = NATFRG(JFG)
              INDJ = INDFRG(JFG)
              DO I = 1, NATI
                IAG = IATFRG(INDI+I-1)
                ISP = ISPEFMO(IAG)
                DO J = 1, NATJ
                  JAG = IATFRG(INDJ+J-1)
                  JSP = ISPEFMO(JAG)
                  DAMPXH = DAMPXH0
                  CALL DFTB_DMPCHK(ISP,JSP,DAMPXH)
                  IF (PERIOD) THEN
                  CALL DFTB_GAMFMO_PBC(IAG,JAG,FMOC(1,IAG),FMOC(1,JAG),
     *                                 GAMMAV,GAMI,GAMJ,TMP(1,ISP,JSP),
     *                                 TMP3(1,ISP,JSP),TMP3(1,JSP,ISP),
     *                                 DAMPXH,ISP,JSP)
                  ELSE
                    XX = FMOC(1,IAG)-FMOC(1,JAG)
                    YY = FMOC(2,IAG)-FMOC(2,JAG)
                    ZZ = FMOC(3,IAG)-FMOC(3,JAG)
                    DIST = XX*XX + YY*YY + ZZ*ZZ
                    IF (DFTB3) THEN
                      CALL DFTB_GAM23(0,DIST,GAMMAV,GAMI,GAMJ,
     *                   TMP(1,ISP,JSP),TMP3(1,ISP,JSP),TMP3(1,JSP,ISP),
     *                     DAMPXH,DAMPXHE,SPE(ISP),SPE(JSP),HUBDER(ISP),
     *                     HUBDER(JSP),DFTB3)
                    ELSE
                      IF (DIST.GE.1.0D-15) THEN
                        DIST = SQRT(DIST)
                        DISTI = ONE/DIST
                      ELSE
                        DISTI = 0
                      END IF
                    GAMMAV = DFTB_GAM2(DIST,DISTI,TMP(1,ISP,JSP),DAMPXH,
     *                  DAMPXHE,SPE(ISP),SPE(JSP))
                    END IF
                  END IF
                  CALL DFTB_PUTGAM(INDI+I-1,INDJ+J-1,NATFMOB,GAMMAV,
     *                             GAMI,GAMJ,GAMMAFMO,GAMMA3,DFTB3)
                END DO
              END DO
            END DO
          END DO
        ELSE
          DO IAG = 1, NATFMO
            IF (DDIGAM) THEN
              CALL GDDICOUNT(0,LGROUP,MYJOB)
              IF (.NOT.MYJOB) CYCLE
            END IF
            ISP = ISPEFMO(IAG)
            DO JAG = 1, IAG
              JSP = ISPEFMO(JAG)
              DAMPXH = DAMPXH0
              CALL DFTB_DMPCHK(ISP,JSP,DAMPXH)
              IF (PERIOD) THEN
                CALL DFTB_GAMFMO_PBC(IAG,JAG,FMOC(1,IAG),FMOC(1,JAG),
     *                               GAMMAV,GAMI,GAMJ,TMP(1,ISP,JSP),
     *                               TMP3(1,ISP,JSP),TMP3(1,JSP,ISP),
     *                               DAMPXH,ISP,JSP)
              ELSE
                XX = FMOC(1,IAG)-FMOC(1,JAG)
                YY = FMOC(2,IAG)-FMOC(2,JAG)
                ZZ = FMOC(3,IAG)-FMOC(3,JAG)
                DIST = XX*XX + YY*YY + ZZ*ZZ
                IF (DFTB3) THEN
                  CALL DFTB_GAM23(0,DIST,GAMMAV,GAMI,GAMJ,
     *                 TMP(1,ISP,JSP),TMP3(1,ISP,JSP),TMP3(1,JSP,ISP),
     *                 DAMPXH,DAMPXHE,SPE(ISP),SPE(JSP),HUBDER(ISP),
     *                 HUBDER(JSP),DFTB3)
                ELSE
                  IF (DIST.GE.1.0D-15) THEN
                    DIST = SQRT(DIST)
                    DISTI = ONE/DIST
                  ELSE
                    DISTI = 0
                  END IF
                  GAMMAV = DFTB_GAM2(DIST,DISTI,TMP(1,ISP,JSP),DAMPXH,
     *              DAMPXHE,SPE(ISP),SPE(JSP))
                END IF
              END IF
              CALL DFTB_CNVSQ(IAG,JAG,NSEQ)
              GAMMAFMO(NSEQ) = GAMMAV
              IF (DFTB3) THEN
C             write (*,'(2i3,2f20.10)') iag,jag,gami,gamj
                GAMMA3(IAG,JAG) = GAMI
                GAMMA3(JAG,IAG) = GAMJ
              END IF
            END DO
          END DO
        END IF
C       CALL CPU_TIME(T1)
C       IF (MASWRK.AND.MDOUT)
C    *    WRITE (IW,'(" CPU TIME FOR GAMMA PREPARATION0= ",F10.5)')T2-T0
C       IF (MASWRK.AND.MDOUT)
C    *    WRITE (IW,'(" CPU TIME FOR GAMMA PREPARATION = ",F10.5)')T1-T0
      END IF
C
      IF (DDIGAM) THEN
        CALL GDDICOUNT( 1,LGROUP,MYJOB)
        CALL TSECND(TIMe0)
        IF (NPROC.GT.1) THEN
          IF (TSQ) THEN
            NN1 = NATFMOB*NATFMOB
            NN2 = 2*NATFMOB*NATFMOB
          ELSE
            NN1 = NATFMO*(NATFMO+1)/2
            NN2 = NATFMO*NATFMO
          END IF
          CALL DSCAL(NN1,1.0D+00/NPROC,GAMMAFMO,1)
          IF (DFTB3) CALL DSCAL(NN2,1.0D+00/NPROC,GAMMA3,1)
        END IF
        CALL GDDI_SCOPE(DDI_WORLD)
        wall0=wall 
        IF (TSQ) THEN
          CALL DDI_GSUMF(2424,GAMMAFMO,NATFMOB*NATFMOB)
          IF (DFTB3) CALL DDI_GSUMF(2424,GAMMA3,NATFMOB*NATFMOB*2)
          ndata = natfmob*natfmob
          if (dftb3) ndata = ndata + natfmob*natfmob*2
        ELSE
          CALL DDI_GSUMF(2424,GAMMAFMO,NATFMO*(NATFMO+1)/2)
          IF (DFTB3) CALL DDI_GSUMF(2424,GAMMA3,NATFMO*NATFMO)
          ndata = natfmo*(natfmo+1)/2
          if (dftb3) ndata = ndata + natfmob*natfmob
        END IF
C       IF (IDFTBD.EQ.4) CALL DDI_GSUMI(2424,NBOND,NATFMOB)
        CALL TSECND(TIMe1)
        if(maswrk.and.mdout) write(iw,9600) ndata,wall-wall0
c       if(maswrk) write(iw,9600) ndata,wall-wall0
 9600 format(/1x,'Data exchange of',I10,' words took',F8.1,' s.')
        CALL GDDI_SCOPE(DDI_GROUP)
      END IF
C       CALL CPU_TIME(T3)
C       IF (MASWRK.AND.MDOUT)
C    *    WRITE (IW,'(" CPU TIME FOR SYNCRONIZATION= ",F10.5)')T3-T1
C     if (maswrk) write (iw,*) "finish synchronization"
C     call timit(1)
C
C     ----- COUNT THE NUMBER OF BONDS PER ATOM -----
C     HERE ASSUMES THAT THERE ARE NO INTER-FRAGMENT COVALENT BONDS
C     OTHER THAN BDA-BAA BONDS.
C
  104 IF (IDFTBD.EQ.4) THEN
        call viclr(nbond,1,natfmob)
        DO IFG = 1, NFG
          NATI = NATFRG(IFG)
          INDI = INDFRG(IFG)
C         DO JFG = 1, NFG
C           NATJ = NATFRG(JFG)
C           INDJ = INDFRG(JFG)
            DO I = 1, NATI
              IAG = IATFRG(INDI+I-1)
              ISP = ISPEFMO(IAG)
C             DO J = 1, NATJ
              DO J = 1, I !! NATI
C               JAG = IATFRG(INDJ+J-1)
C               JSP = ISPEFMO(JAG)
                JAG = IATFRG(INDI+J-1)
                JSP = ISPEFMO(JAG)
                IF (PERIOD) THEN
                  ITRVEC = 1
                  ITRVECX = ITRVEC
                  ITRVECY = ITRVEC
                  ITRVECZ = ITRVEC
                  DO IXVEC = -ITRVECX, ITRVECX
                  DO IYVEC = -ITRVECY, ITRVECY
                  DO IZVEC = -ITRVECZ, ITRVECZ
                    VEC(1) = FMOC(1,IAG)-FMOC(1,JAG)
     *                     - PBCBOX(1,1)*IXVEC
     *                     - PBCBOX(2,1)*IYVEC
     *                     - PBCBOX(3,1)*IZVEC
                    VEC(2) = FMOC(2,IAG)-FMOC(2,JAG)
     *                     - PBCBOX(1,2)*IXVEC
     *                     - PBCBOX(2,2)*IYVEC
     *                     - PBCBOX(3,2)*IZVEC
                    VEC(3) = FMOC(3,IAG)-FMOC(3,JAG)
     *                     - PBCBOX(1,3)*IXVEC
     *                     - PBCBOX(2,3)*IYVEC
     *                     - PBCBOX(3,3)*IZVEC
                    DIST = VEC(1)*VEC(1)+VEC(2)*VEC(2)+VEC(3)*VEC(3)
                    DIST = SQRT(DIST)
                    R0A = DFTBDP(ISP*14-13) !! COVALENT RADIUS
                    R0B = DFTBDP(JSP*14-13)
                    IF (DIST.LT.(R0A+R0B).AND.DIST.GT.1.0D-08) THEN
                        NBOND(IAG) = NBOND(IAG) + 1
                        NBOND(JAG) = NBOND(JAG) + 1
                    END IF
                  END DO
                  END DO
                  END DO
                ELSE
                  IF (IAG.EQ.JAG) CYCLE
                  DIST = (FMOC(1,IAG)-FMOC(1,JAG))**2
     *                 + (FMOC(2,IAG)-FMOC(2,JAG))**2
     *                 + (FMOC(3,IAG)-FMOC(3,JAG))**2
                  DIST = SQRT(DIST)
                  R0A = DFTBDP(ISP*14-13) !! COVALENT RADIUS
                  R0B = DFTBDP(JSP*14-13)
                  IF (DIST.LT.(R0A+R0B)) THEN
                      NBOND(IAG) = NBOND(IAG) + 1
                      NBOND(JAG) = NBOND(JAG) + 1
                  END IF
                END IF
              END DO
            END DO
C         END DO
        END DO
      END IF
C
      IF (IDFTBD.EQ.4) THEN
        DO I = 1, NATFMOB
          IF (NBOND(I).GT.5) NBOND(I) = 5 !! MAX NEIGHBORS SHOULD BE 5
        END DO
      END IF
C
      RETURN
C
      END SUBROUTINE DFTB_GAMFMO
C
C-----------------------------------------------------------------------
C*MODULE DFTBFO    *DECK DFTB_GAM23
C>
C>    @brief   Calculate gamma values in FMO-DFTB
C>
C>    @details Calculate gamma (DFTB2 and DFTB3) values in FMO-DFTB
C>             simultaneously to save computational timing.
C>
C>    @author Yoshio Nishimoto
C>    - Feb, 2016- Subroutine written
C>
C>           --- INPUT ---
C>    @param IM       Mode: calculate SQRT(DIST) or not
C>    @param DIST     Distance
C>    @param VALUE2   Some constants for DFTB2 gamma
C>    @param VALUE3A  Some constants for DFTB3(AB) gamma
C>    @param VALUE3B  Some constants for DFTB3(BA) gamma
C>    @param DAMPXH   Damping is used or not
C>    @param DAMPXHE  Damping exponent
C>    @param SPE1     Name of atom A
C>    @param SPE2     Name of atom B
C>    @param HUBDERI  Hubbard derivative of atom A
C>    @param HUBDERJ  Hubbard derivative of atom B
C>    @param DFTB3    DFTB3 or not
C>           --- OUTPUT ---
C>    @param GAM2     DFTB2 gamma (\gamma_{AB})
C>    @param GAM3I    DFTB3 gamma (\Gamma_{AB})
C>    @param GAM3J    DFTB3 gamma (\Gamma_{BA})
C
C
      SUBROUTINE DFTB_GAM23(IM,DIST,GAM2,GAM3I,GAM3J,VALUE2,VALUE3A,
     *  VALUE3B,DAMPXH,DAMPXHE,SPE1,SPE2,HUBDERI,HUBDERJ,DFTB3)
C
      IMPLICIT NONE
C
      DOUBLE PRECISION, PARAMETER :: ZERO=0.0D+00,DISTTOL=1.0D-08,
     *  HUBTOL=1.0D-04,HALF=0.5D+00,ONE=1.0D+00
C
      DOUBLE PRECISION, INTENT(IN) :: VALUE2(12),VALUE3A(10),
     *  VALUE3B(10),DAMPXHE,SPE1,SPE2,HUBDERI,HUBDERJ
      DOUBLE PRECISION, INTENT(OUT) :: GAM2,GAM3I,GAM3J
      DOUBLE PRECISION, INTENT(INOUT) :: DIST
      INTEGER, INTENT(IN) :: IM
      LOGICAL, INTENT(IN) :: DAMPXH,DFTB3
C
      DOUBLE PRECISION :: TAUA,TAUB,TAUMEAN,TAUAB3I,DISTI,
     * EXPA,EXPB,DIST2,DIST3,
     * DSAB_DR,DSBA_DR,DFAB_DA,DFBA_DA,DFAB_DB,DFBA_DB,FAB,FBA,HUBMEAN,
     * HUBMEANEXP,HUBMEANEXP1,SAB,DH_DU
      DOUBLE PRECISION H/8HH       /
C
      IF (DIST.LE.DISTTOL) THEN
        GAM2 = VALUE2(12)
        IF (DFTB3) THEN
          GAM3I = HALF*HUBDERI
          GAM3J = HALF*HUBDERJ
        END IF
        RETURN
      ELSE
        IF (IM.EQ.0) DIST = SQRT(DIST)
        DISTI = ONE/DIST
        TAUA = VALUE2(1)
        TAUB = VALUE2(2)
        DH_DU = ZERO
        DIST2 = DIST*DIST
        IF (DAMPXH.AND.(SPE1.EQ.H.OR.SPE2.EQ.H)) THEN
          IF (ABS(TAUA-TAUB).LE.HUBTOL) THEN
            TAUMEAN = VALUE2(3)
            HUBMEAN = VALUE2(12)
            HUBMEANEXP1 = HUBMEAN**(DAMPXHE-ONE)
            HUBMEANEXP  = HUBMEANEXP1*HUBMEAN
            EXPA = EXP(-TAUMEAN*DIST-HUBMEANEXP*DIST2)
            SAB = EXPA*(DISTI+VALUE2(4)+VALUE2(5)*DIST+VALUE2(6)*DIST2)
            GAM2 = DISTI - SAB
            IF (DFTB3) THEN
              !! DS^g/Da
              DIST3 = DIST2*DIST
              DSAB_DR = EXPA*(VALUE3A(3)+VALUE3A(4)*DIST
     *          + VALUE3A(5)*DIST2+VALUE3A(6)*DIST3)
              DSBA_DR = DSAB_DR
              DH_DU = -DIST2*DAMPXHE*HUBMEANEXP1*HALF
              GAM3I = -HUBDERI*(DSAB_DR*3.2D+00 + SAB*DH_DU)
              GAM3J = -HUBDERJ*(DSBA_DR*3.2D+00 + SAB*DH_DU)
            END IF
          ELSE
            TAUAB3I = DISTI*VALUE2(11)
            HUBMEAN = VALUE2(12)
            HUBMEANEXP1 = HUBMEAN**(DAMPXHE-ONE)
            HUBMEANEXP  = HUBMEANEXP1*HUBMEAN
            EXPA = EXP(-TAUA*DIST-HUBMEANEXP*DIST2)
            EXPB = EXP(-TAUB*DIST-HUBMEANEXP*DIST2)
            SAB = EXPA*(VALUE2(7)-VALUE2( 9)*TAUAB3I)
     &          + EXPB*(VALUE2(8)+VALUE2(10)*TAUAB3I)
            GAM2 = DISTI - SAB
            IF (DFTB3) THEN
              !! DS^f/Da
              DFAB_DA = VALUE3A( 3) + VALUE3A( 4)*DISTI
              FAB     = VALUE3A( 5) + VALUE3A( 6)*DISTI
              DFBA_DA = VALUE3A( 7) + VALUE3A( 8)*DISTI
              FBA     = VALUE3A( 9) + VALUE3A(10)*DISTI
              DFAB_DB = VALUE3B( 7) + VALUE3B( 8)*DISTI
              DFBA_DB = VALUE3B( 3) + VALUE3B( 4)*DISTI
              DSAB_DR = EXPA*(DFAB_DA-DIST*FAB) + EXPB*DFBA_DA
              DSBA_DR = EXPB*(DFBA_DB-DIST*FBA) + EXPA*DFAB_DB
              DH_DU = -DIST2*DAMPXHE*HUBMEANEXP1*HALF
              GAM3I = -HUBDERI*(DSAB_DR*3.2D+00 + SAB*DH_DU)
              GAM3J = -HUBDERJ*(DSBA_DR*3.2D+00 + SAB*DH_DU)
            END IF
          END IF
        ELSE
          IF (ABS(TAUA-TAUB).LE.HUBTOL) THEN
            TAUMEAN = VALUE2(3)
            EXPA = EXP(-TAUMEAN*DIST)
            SAB = EXPA*(DISTI+VALUE2(4)+VALUE2(5)*DIST+VALUE2(6)*DIST2)
            GAM2 = DISTI - SAB
            IF (DFTB3) THEN
              !! DS^g/Da
              DIST3 = DIST2*DIST
              DSAB_DR = EXPA*(VALUE3A(3)+VALUE3A(4)*DIST
     *          + VALUE3A(5)*DIST2+VALUE3A(6)*DIST3)
              DSBA_DR = DSAB_DR
              GAM3I = -HUBDERI*DSAB_DR*3.2D+00
              GAM3J = -HUBDERJ*DSBA_DR*3.2D+00
            END IF
          ELSE
            TAUAB3I = DISTI*VALUE2(11)
            EXPA = EXP(-TAUA*DIST)
            EXPB = EXP(-TAUB*DIST)
            SAB = EXPA*(VALUE2(7)-VALUE2( 9)*TAUAB3I)
     &          + EXPB*(VALUE2(8)+VALUE2(10)*TAUAB3I)
            GAM2 = DISTI - SAB
            IF (DFTB3) THEN
              !! DS^f/Da
              DFAB_DA = VALUE3A( 3) + VALUE3A( 4)*DISTI
              FAB     = VALUE3A( 5) + VALUE3A( 6)*DISTI
              DFBA_DA = VALUE3A( 7) + VALUE3A( 8)*DISTI
              FBA     = VALUE3A( 9) + VALUE3A(10)*DISTI
              DFAB_DB = VALUE3B( 7) + VALUE3B( 8)*DISTI
              DFBA_DB = VALUE3B( 3) + VALUE3B( 4)*DISTI
              DSAB_DR = EXPA*(DFAB_DA-DIST*FAB) + EXPB*DFBA_DA
              DSBA_DR = EXPB*(DFBA_DB-DIST*FBA) + EXPA*DFAB_DB
              GAM3I = -HUBDERI*DSAB_DR*3.2D+00
              GAM3J = -HUBDERJ*DSBA_DR*3.2D+00
            END IF
          END IF
        END IF
      END IF
C
      RETURN
C
C     IF (DAMPXH) THEN
C       IF (SPE1.EQ.H.OR.SPE2.EQ.H) THEN
C         HUBMEAN = VALUE2(12)
C         HUBMEANEXP1 = HUBMEAN**(DAMPXHE-ONE)
C         HUBMEANEXP  = HUBMEANEXP1*HUBMEAN
C         EXPH = EXP(-HUBMEANEXP*DIST2)
C         IF (DFTB3) THEN
C           DH_DU = -DIST2*DAMPXHE*HUBMEANEXP1*EXPH*HALF
C         END IF
C       END IF
C     END IF
C
C     GAM2 = DISTI - SAB*EXPH
C     IF (DFTB3) THEN
C       GAM3I = -HUBDERI*(DSAB_DR*EXPH*3.2D+00 + SAB*DH_DU)
C       GAM3J = -HUBDERJ*(DSBA_DR*EXPH*3.2D+00 + SAB*DH_DU)
C     END IF
C
C     RETURN
C
      END SUBROUTINE DFTB_GAM23
C
C-----------------------------------------------------------------------
C
      SUBROUTINE DFTB_PUTGAM(IAG,JAG,NATFMOB,GAM2,GAMI,GAMJ,GAMMA2,
     *  GAMMA3,DFTB3)
C
      IMPLICIT NONE
C
      INTEGER, INTENT(IN) :: IAG,JAG,NATFMOB !! NOT TRUE IAG/JAG
      DOUBLE PRECISION, INTENT(IN) :: GAM2,GAMI,GAMJ
      DOUBLE PRECISION, INTENT(INOUT) :: GAMMA2(NATFMOB,*),
C    *  GAMMA3(2,NATFMOB,*)
     *  GAMMA3(NATFMOB,NATFMOB,*)
      LOGICAL, INTENT(IN) :: DFTB3
C
      GAMMA2(IAG,JAG) = GAM2
      GAMMA2(JAG,IAG) = GAM2
      IF (DFTB3) THEN
C       GAMMA3(1,IAG,JAG) = GAMI
C       GAMMA3(2,IAG,JAG) = GAMJ
C       GAMMA3(1,JAG,IAG) = GAMJ
C       GAMMA3(2,JAG,IAG) = GAMI
        GAMMA3(IAG,JAG,1) = GAMI
        GAMMA3(IAG,JAG,2) = GAMJ
        GAMMA3(JAG,IAG,1) = GAMJ
        GAMMA3(JAG,IAG,2) = GAMI
      END IF
C
      RETURN
C
      END SUBROUTINE DFTB_PUTGAM
C
C-----------------------------------------------------------------------
C*MODULE DFTBFO    *DECK DFTB_DCFMO
C>
C>    @brief Dispersion correction for FMO-DFTB
C>
C>    @details Prepare dispersion correction for FMO-DFTB
C>             If called with MODGAM>=0, only constants are calculated
C>             here and returned them as DCTMP.
C>             If called with MODGAM<0, dispersion interaction between
C>             each atom is calculated here, however this does not
C>             improve performance and efficiency.
C>
C>    @author Yoshio Nishimoto
C>    - Sep, 2014- Subroutine written
C>
C>           --- INPUT ---
C>    @param MODGAM Whether only constants are calculated or each
C>           atom pair interaction is calculated.
C>    @param NATFMO Number of atoms
C>    @param ISPEFMO ISPE (Index of species) for all NATFMO,
C>           prepared in DFTB_GAMFMO subroutine
C>    @param C Coordinate
C>           --- OUTPUT ---
C>    @param DCTMP (MODGAM>=0) Constants used for dispersion correction
C>                 (MODGAM<0) each potential
C>
C
      SUBROUTINE DFTB_DCFMO(MODGAM,NATFMO,DCTMP,ISPEFMO,C)
      use mx_limits, only: mxatm
C
      IMPLICIT NONE
C
      INTEGER,PARAMETER :: MXSPE=10
      DOUBLE PRECISION, PARAMETER :: TWO=2.0D+00
C
      COMMON /DFTBPR/ ETEMP,DFTBDP(MXSPE*14),DAMPXHE,HUBDER(MXSPE),
     *                ZREF(MXATM),SPNCST(6,MXSPE),SPE(MXATM),NSPE,
     *                MAXANG(MXATM),ISPE(MXATM),IND(MXATM+1),IDFTBD,
     *                PARAMDIR
      DOUBLE PRECISION :: ETEMP,DFTBDP,DAMPXHE,HUBDER,ZREF,SPNCST,SPE
      INTEGER :: NSPE,MAXANG,ISPE,IND,IDFTBD
C
      INTEGER :: MODGAM,NATFMO,I,J,NSEQ,NSP1,NSP2,ISPEFMO(*)
      DOUBLE PRECISION :: CR0,RI,RJ,DI,DJ,RIJ,DIJ,R0,PA,R0A,NEA,C6A,
     * PB,R0B,NEB,C6B,C6AB,R0AB,R0AB7I,DIST,R6,R12,TMP,R5,R10,DISTI,
     * DIST6I,DIST7,VEXP,DAMP1,DAMP4,C0,C1,C2
      DOUBLE PRECISION :: WRK(3,NSPE,NSPE),DCTMP(*),C(3,*),PARAMDIR
C
      IF (IDFTBD.EQ.1) THEN
        DO I = 1, NSPE
          CR0 = 8.908987181D-01
          RI = DFTBDP(I*2-1)
          DI = DFTBDP(I*2  )
          DO J = I, NSPE
            RJ = DFTBDP(J*2-1)
            DJ = DFTBDP(J*2  )
            RIJ = SQRT(RI*RJ) !! van der Waals distance
            DIJ = SQRT(DI*DJ) !! well depth
            R0 = CR0*RIJ
            WRK(1,I,J) = R0
            WRK(2,I,J) = RIJ
            WRK(3,I,J) = DIJ
            WRK(1,J,I) = R0
            WRK(2,J,I) = RIJ
            WRK(3,J,I) = DIJ
          END DO
        END DO
      ELSE IF (IDFTBD.EQ.2) THEN
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
            WRK(1,I,J) = C6AB
            WRK(2,I,J) = R0AB7I
            WRK(1,J,I) = C6AB
            WRK(2,J,I) = R0AB7I
          END DO
        END DO
      END IF
C
      IF (MODGAM.GE.0) THEN
        CALL DCOPY(3*NSPE*NSPE,WRK,1,DCTMP,1)
        RETURN
      END IF
C
      IF (IDFTBD.EQ.1) THEN
        !! CALCULATE DISPERSION CORRECTION
        C0 = 1.584D+01 !! 396/25
        C1 = 4.78947150872246D+01 !! 2^(5/6)*672/25
        C2 = 3.50498152274578D+01 !! 2^(2/3)*552/25
        DO I = 1, NATFMO
          NSP1 = ISPEFMO(I)
          DO J = I+1, NATFMO
            NSP2 = ISPEFMO(J)
            R0  = WRK(1,NSP1,NSP2)
            RIJ = WRK(2,NSP1,NSP2)
            DIJ = WRK(3,NSP1,NSP2)
            DIST = SQRT( (C(1,I)-C(1,J))**2
     &                 + (C(2,I)-C(2,J))**2
     &                 + (C(3,I)-C(3,J))**2 )
            IF (DIST.GE.R0) THEN !! ATTRACTIVE POTENTIAL
              TMP = RIJ/DIST
              R6 = TMP**6
              R12 = R6**2
              TMP = DIJ*(-TWO*R6+R12)
            ELSE !! REPULSIVE POTENTIAL
              TMP = DIST/RIJ
              R5 = TMP**5
              R10 = R5**2
              TMP = DIJ*(C0-C1*R5+C2*R10)
            END IF
            CALL DFTB_CNVSQ(I,J,NSEQ)
            DCTMP(NSEQ) = TMP
          END DO
        END DO
      ELSE IF (IDFTBD.EQ.2) THEN
        DO I = 1, NATFMO
          NSP1 = ISPEFMO(I)
          DO J = I+1, NATFMO
            NSP2 = ISPEFMO(J)
            C6AB   = WRK(1,NSP1,NSP2)
            R0AB7I = WRK(2,NSP1,NSP2)
            DIST = SQRT( (C(1,I)-C(1,J))**2
     &                 + (C(2,I)-C(2,J))**2
     &                 + (C(3,I)-C(3,J))**2 )
            DISTI = 1.0D+00 / DIST
            DIST6I = DISTI**6
            DIST7  = DIST**7
            VEXP = EXP(-3.0D+00*R0AB7I*DIST7)
            DAMP1 = 1.0D+00 - VEXP
            DAMP4 = DAMP1**4 !! THIS IS THE DAMPING FUNCTION (f(R))
            TMP = -1.0D+00 * DAMP4 * C6AB * DIST6I
            CALL DFTB_CNVSQ(I,J,NSEQ)
            DCTMP(NSEQ) = TMP
          END DO
        END DO
      END IF
C
      END SUBROUTINE DFTB_DCFMO
C
C-----------------------------------------------------------------------
C*MODULE DFTBFO    *DECK DFTB_ESDIM
C>
C>    @brief ES-DIM for FMO-DFTB
C>
C>    @details Calculate electrostatic dimer (ES-DIM) contribution to
C>             energy at FMO-DFTB.  This is used only for MODESD=0,
C>             which is slowest ES-DIM option, but useful to get
C>             complete dimer energy table.
C>
C>    @author Yoshio Nishimoto
C>    - Sep, 2014- Subroutine written
C>
C>           --- INPUT ---
C>    @param NDER       Order of dervative to be calculated
C>    @param IFG,JFG    I and J fragment index
C>    @param POPMAT     Mulliken population
C>    @param ILAY       Index of layer?
C>    @param ESDREP     Repulsive energy is calculated or not
C>    @param POPMATI    Mulliken populations of IFG
C>    @param POPMATJ    Mulliken populations of JFG
C>    @param DIST       Distance matrix between IFG-JFG
C>           --- OUTPUT ---
C>    @param DFTB_ESDIM ES-DIM energy contribution
C>
C
C     DOUBLE PRECISION FUNCTION DFTB_ESDIM(nder,ifg,jfg,LAYFRG,
C    *  IAGLOB,POPMAT)
      DOUBLE PRECISION FUNCTION DFTB_ESDIM(nder,ifg,jfg,POPMAT,
     *  ilay,ESDREP,POPMATI,POPMATJ,DIST)
      use mx_limits, only: mxatm,mxrt
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      INTEGER,PARAMETER :: MXSPE=10
C
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
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
      COMMON /ENRGYS/ ENUCR,EELCT,ETOT,SZ,SZZ,ECORE,ESCF,EERD,E1,E2,
     *                VEN,VEE,EPOT,EKIN,ESTATE(MXRT),STATN,EDFT(3),EDISP
      COMMON /FMCOM / X(1)
      common /fmodim/ maxbnd,maxknd,maxcbs,maxcao,maxvec,maxl1,maxnat,
     *                maxabd,maxbas,maxbbd,maxlmo,maxslo,maxabd2,maxrij
      COMMON /FMODTB/ LGAMMA,LZREF,LISPE,LSHIFTG,LSHIFTCT,
     *                LGAMMAD,MODESD,MODGAMMA,LDCTMP,LGAMMA3,LGAMMA3D,
     *                LNBOND,LCTMUL,LESPDFTB,LESPASC,NDFTBRST,LDFTBRST
      COMMON /FMOINF/ NFG,NLAYER,NATFMO,NBDFG,NAOTYP,NBODY,NSEGM
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
      common /GRAD  / DE(3,MXATM)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      LOGICAL GOPARR,DSKWRK,MASWRK
C
      LOGICAL DFTBFL,SCC,SRSCC,DFTB3,DAMPXH,LCDFTB
      DOUBLE PRECISION popmat(maxnat,nfg)
      LOGICAL, INTENT(IN) :: ESDREP
      DOUBLE PRECISION :: POPMATI(*),POPMATJ(*),DIST(*)
C     integer layfrg(*)
C
C     ----- FIRST CALCULATE EREP -----
C
      IF (ESDREP) THEN
        call makemol(ifg,jfg,0,ilay,0,0,0,0,0,0,0,.false.)
        CALL DFTB_EREP(NAT,C(1:3,1:MXATM),ISPE,NSPE,NUMREP,REPCUT,
     *         X(LREPINTV(1,1)),X(LREPSHORT(1,1)),X(LREPCOEFF(1,1)),
     *         ENUCR,DIST,ZREF,QREF)
      END IF
C
C     ----- CALCULATE ESD ENERGY -----
C
      DFTB_ESDIM = 0.0D+00
C
C     DE HAS THE GRADIENT IN THE PREVIOUS CYCLE, SO IT HAS TO BE
C     INITIALIZED.
C
C     ----- NOTHING TO BE DONE FOR NCC CALCULATION -----
C
      IF (.NOT.SCC) THEN
        CALL VCLR(DE,1,3*NAT)
        RETURN
      END IF
C
C     ----- CALCULATE ESP ENERGY FOR SEPARATED DIMER -----
C
      IF (IAND(MODGAMMA,2).NE.0) THEN
        NN = 2
        NN2 = NATFMO+NBDFG
      ELSE
        NN = 1
        NN2 = NATFMO
      END IF
      CALL DFTB_ESDIM2(IFG,JFG,MODGAMMA,X(LNATFRG),X(LINDFRG),
     *  X(LIATFRG),MAXNAT,X(LISPE),NSPE,X(LGAMMA),X(LZREF),POPMAT,
     *  X(LFMOC),DFTB_ESDIM,DFTB3,HUBDER,X(LGAMMA3),NATFMO,HUBB,DAMPXH,
     *  DAMPXHE,SPE,NATFMO+NBDFG,NN,NN2,POPMATI,POPMATJ,DIST)
C
      IF (NDER.EQ.0) RETURN
C
C     ----- CALCULATE GRADIENTS -----
C
      CALL VCLR(DE,1,3*NAT)
      if (nbody.eq.3) then
        iifg=max(ifg,jfg)
        jjfg=min(ifg,jfg)
        ijfg=(iifg*iifg-3*iifg)/2+jjfg+1
        ida=1
        ida=ida-ixftch(x(lndtfrg),ijfg)
        SCAL = DBLE(IDA)
      else
        SCAL = 1.0D+00
      end if
      IF (MASWRK) CALL DFTB_ESDGRAD(IFG,JFG,X(LISPE),
     *  X(LINDFRG),X(LNATFRG),X(LIATFRG),X(LGAMMAD),
     *  X(LSHIFTG+NATFMO+NBDFG),X(LFMOC),X(LFMODE),X(LGAMMA3D),POPMATI,
     *  POPMATJ,X(LGAMMA3),DIST,NN,NN2,SCAL)
C
      RETURN

      END FUNCTION DFTB_ESDIM
C
C-----------------------------------------------------------------------
C*MODULE DFTBFO    *DECK DFTB_ESDIM_DC
C>
C>    @brief ES-DIM for FMO-DFTB
C>
C>    @details Calculate dispersion correction in calculating the
C>             electrostatic dimer (ES-DIM) contribution to
C>             energy at FMO-DFTB.  This is called from
C>             DFTB_ESDIM_DC_HELP subroutine.
C>             For MODESD=0, initial EDISP is zero, however in other
C>             cases, EDISP should have certain value.
C>
C>    @author Yoshio Nishimoto
C>    - Sep, 2014- Subroutine written
C>    @date   Feb, 2016 - Yoshio Nishimoto
C>    - Added Slater-Kirkwoord (Hybrid-Polarized)
C>
C>           --- INPUT ---
C>    @param IFG,JFG I and J fragment index
C>    @param NDER    Order of dervative to be calculated
C>    @param NATFMO  Number of atoms in the system
C>    @param MODGAM  How gamma values are calculated
C>    @param NATFRG  Number of atoms in each fragment
C>    @param INDFRG  First atom of each fragment in terms of NATFMOB?
C>    @param IATFRG  Global atom index with NBDFG?
C>    @param ISPEFMO ISPE for all atoms (NATFMO)
C>    @param DCTMP   Constants for dispersion correction, see DFTB_DCFMO
C>    @param FMOC    Coordinates of allatoms (NATFMO)
C>    @param ZREFFMO ZREF for all atoms (NATFMO)
C>    @param DISTMAT Distance matrix between IFG-JFG
C>    @param NBOND   Number of bonds
C>           --- IN/OUTPUT ---
C>    @param EDISP   Dispersion energy
C>
C
      SUBROUTINE DFTB_ESDIM_DC(IFG,JFG,NDER,NATFMO,MODGAM,NATFRG,INDFRG,
     *  IATFRG,ISPEFMO,DCTMP,FMOC,FMODE,ZREFFMO,DISTMAT,NBOND,EDISP,
     *  SCAL,NBODY)
      use mx_limits, only: mxatm
C
      IMPLICIT NONE
C
      INTEGER,PARAMETER :: MXSPE=10
      DOUBLE PRECISION, PARAMETER :: TWO=2.0D+00
C
      COMMON /DFTB  / DFTBFL,SCC,SRSCC,DFTB3,DAMPXH,LCDFTB
      LOGICAL DFTBFL,SCC,SRSCC,DFTB3,DAMPXH,LCDFTB
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
      INTEGER, INTENT(IN) :: IFG,JFG,NDER,NATFRG(*),INDFRG(*),
     *  IATFRG(*),ISPEFMO(*),NATFMO,MODGAM,NBOND(*),NBODY
      DOUBLE PRECISION, INTENT(IN) :: DCTMP(*),FMOC(3,*),
     *  ZREFFMO(*),SCAL
      DOUBLE PRECISION, INTENT(INOUT) :: FMODE(3,NATFMO,*),EDISP,
     *  DISTMAT(*)
C
      INTEGER :: NATI,NATJ,INDI,INDJ,IAG,JAG,I,J,K,ISP,JSP,NSEQ
      DOUBLE PRECISION :: R0,RIJ,DIJ,DIST,TMP,R6,R12,R5,R10,VAL,
     * C0,C1,C2,CI,CJ,C6AB,R0AB7I,DISTI,DIST6I,DIST7I,DIST7,DAMP1,
     * DAMP3,DAMP4,VEXP,R0A,R0B,PA,PB,C6A,C6B,R0AB,NEA,NEB
      DOUBLE PRECISION :: RDMAT(3,NSPE,NSPE)
C
      NATI = NATFRG(IFG)
      NATJ = NATFRG(JFG)
      INDI = INDFRG(IFG)
      INDJ = INDFRG(JFG)
C
      IF (MODGAM.LT.0) THEN
        CI = 1.0D+00
        CJ = 1.0D+00
        DO I = 1, NATI
          IAG = IATFRG(INDI+I-1)
          ISP = ISPEFMO(IAG)
          IF (QREF(ISP).NE.ZREFFMO(INDI+I-1)) THEN
            CI = ZREFFMO(INDI+I-1)/QREF(ISP)
          ELSE
            CI = 1.0D+00
          END IF
          DO J = 1, NATJ
            JAG = IATFRG(INDJ+J-1)
            JSP = ISPEFMO(JAG)
            IF (QREF(JSP).NE.ZREFFMO(INDJ+J-1)) THEN
              CJ = ZREFFMO(INDJ+J-1)/QREF(JSP)
            ELSE
              CJ = 1.0D+00
            END IF
            CALL DFTB_CNVSQ(IAG,JAG,NSEQ)
            TMP = DCTMP(NSEQ)
            EDISP = EDISP + TMP*CI*CJ
          END DO
        END DO
        RETURN
      END IF
C
      IF (IAND(ABS(MODGAM),3).NE.0.OR..NOT.SCC) THEN
        NSEQ = 0
        DO I = 1, NATI
          IAG = IATFRG(INDI+I-1)
          DO J = 1, NATJ
            JAG = IATFRG(INDJ+J-1)
            DIST = SQRT( (FMOC(1,IAG)-FMOC(1,JAG))**2
     &                 + (FMOC(2,IAG)-FMOC(2,JAG))**2
     &                 + (FMOC(3,IAG)-FMOC(3,JAG))**2 )
            NSEQ = NSEQ + 1
            DISTMAT(NSEQ) = DIST
          END DO
        END DO
      END IF
C
C     EDISP = 0.0D+00 !! DO NOT INITIALIZE INSIDE
      CALL DCOPY(3*NSPE*NSPE,DCTMP,1,RDMAT,1)
      IF (IDFTBD.EQ.1) THEN
        C0 = 1.584D+01 !! 396/25
        C1 = 4.78947150872246D+01 !! 2^(5/6)*672/25
        C2 = 3.50498152274578D+01 !! 2^(2/3)*552/25
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
C           DIST = SQRT( (FMOC(1,IAG)-FMOC(1,JAG))**2
C    &                 + (FMOC(2,IAG)-FMOC(2,JAG))**2
C    &                 + (FMOC(3,IAG)-FMOC(3,JAG))**2 )
            DIST = DISTMAT(NSEQ)
            DIJ = DIJ * CI * CJ
            IF (DIST.GE.R0) THEN !! ATTRACTIVE POTENTIAL
              TMP = RIJ/DIST
              R6 = TMP**6
              R12 = R6**2
              TMP = DIJ*(-TWO*R6+R12)
              IF (NDER.NE.0) THEN
                VAL = 1.2D+01*DIJ*(R6-R12)/(DIST**2)
                DO K = 1, 3
                  FMODE(K,IAG,2) = FMODE(K,IAG,2)
     *                             + (FMOC(K,IAG)-FMOC(K,JAG))*VAL
                  FMODE(K,JAG,2) = FMODE(K,JAG,2)
     *                             - (FMOC(K,IAG)-FMOC(K,JAG))*VAL
                  IF (NBODY.EQ.3) THEN
                    FMODE(K,IAG,3) = FMODE(K,IAG,3)
     *                + (FMOC(K,IAG)-FMOC(K,JAG))*VAL*SCAL
                    FMODE(K,JAG,3) = FMODE(K,JAG,3)
     *                - (FMOC(K,IAG)-FMOC(K,JAG))*VAL*SCAL
                  END IF
                END DO
              END IF
            ELSE !! REPULSIVE POTENTIAL
              TMP = DIST/RIJ
              R5 = TMP**5
              R10 = R5**2
              TMP = DIJ*(C0-C1*R5+C2*R10)
              IF (NDER.NE.0) THEN
                VAL = DIJ*(-C1*R5*5.0D+00+C2*R10*1.0D+01)/(DIST**2)
                DO K = 1, 3
                  FMODE(K,IAG,2) = FMODE(K,IAG,2)
     *                             + (FMOC(K,IAG)-FMOC(K,JAG))*VAL
                  FMODE(K,JAG,2) = FMODE(K,JAG,2)
     *                             - (FMOC(K,IAG)-FMOC(K,JAG))*VAL
                  IF (NBODY.EQ.3) THEN
                    FMODE(K,IAG,3) = FMODE(K,IAG,3)
     *                + (FMOC(K,IAG)-FMOC(K,JAG))*VAL*SCAL
                    FMODE(K,JAG,3) = FMODE(K,JAG,3)
     *                - (FMOC(K,IAG)-FMOC(K,JAG))*VAL*SCAL
                  END IF
                END DO
              END IF
            END IF
            EDISP = EDISP + TMP
          END DO
        END DO
      ELSE IF (IDFTBD.EQ.2) THEN
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
            C6AB   = RDMAT(1,ISP,JSP)
            R0AB7I = RDMAT(2,ISP,JSP)
C           DIST = SQRT( (FMOC(1,IAG)-FMOC(1,JAG))**2
C    &                 + (FMOC(2,IAG)-FMOC(2,JAG))**2
C    &                 + (FMOC(3,IAG)-FMOC(3,JAG))**2 )
            DIST = DISTMAT(NSEQ)
            C6AB = C6AB * CI * CJ
            DISTI = 1.0D+00 / DIST
            DIST6I = DISTI**6
            DIST7  = DIST**7
            DIST7I = DIST6I*DISTI
            VEXP = EXP(-3.0D+00*R0AB7I*DIST7)
            DAMP1 = 1.0D+00 - VEXP
            DAMP4 = DAMP1**4 !! THIS IS THE DAMPING FUNCTION (f(R))
            TMP = -1.0D+00 * DAMP4 * C6AB * DIST6I
            EDISP = EDISP + TMP
            IF (NDER.NE.0) THEN
              DAMP3 = DAMP1**3
              VAL = 8.4D+01 * R0AB7I * VEXP - 6.0D+00 * DAMP1 * DIST7I
              VAL = -VAL * C6AB * DAMP3 * DISTI
              DO K = 1, 3
                FMODE(K,IAG,2) = FMODE(K,IAG,2)
     *                             + (FMOC(K,IAG)-FMOC(K,JAG))*VAL
                FMODE(K,JAG,2) = FMODE(K,JAG,2)
     *                             - (FMOC(K,IAG)-FMOC(K,JAG))*VAL
                IF (NBODY.EQ.3) THEN
                  FMODE(K,IAG,3) = FMODE(K,IAG,3)
     *              + (FMOC(K,IAG)-FMOC(K,JAG))*VAL*SCAL
                  FMODE(K,JAG,3) = FMODE(K,JAG,3)
     *              - (FMOC(K,IAG)-FMOC(K,JAG))*VAL*SCAL
                END IF
              END DO
            END IF
          END DO
        END DO
      ELSE IF (IDFTBD.EQ.4) THEN
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
            NSEQ = NSEQ + 1
            DIST = DISTMAT(NSEQ)
            C6AB = C6AB * CI * CJ
            DISTI = 1.0D+00 / DIST
            DIST6I = DISTI**6
            DIST7  = DIST**7
            VEXP = EXP(-3.0D+00*R0AB7I*DIST7)
            DAMP1 = 1.0D+00 - VEXP
            DAMP4 = DAMP1**4 !! THIS IS THE DAMPING FUNCTION (f(R))
            TMP = -1.0D+00 * DAMP4 * C6AB * DIST6I
            EDISP = EDISP + TMP
            IF (NDER.NE.0) THEN
              DIST7I = DIST6I*DISTI
              DAMP3 = DAMP1**3
              VAL = 8.4D+01*VEXP*R0AB7I*DIST7-6.0D+00*DAMP1
              VAL = -VAL*C6AB*DAMP3*DIST7I*DISTI
              DO K = 1, 3
                FMODE(K,IAG,2) = FMODE(K,IAG,2)
     *                             + (FMOC(K,IAG)-FMOC(K,JAG))*VAL
                FMODE(K,JAG,2) = FMODE(K,JAG,2)
     *                             - (FMOC(K,IAG)-FMOC(K,JAG))*VAL
                IF (NBODY.EQ.3) THEN
                  FMODE(K,IAG,3) = FMODE(K,IAG,3)
     *              + (FMOC(K,IAG)-FMOC(K,JAG))*VAL*SCAL
                  FMODE(K,JAG,3) = FMODE(K,JAG,3)
     *              - (FMOC(K,IAG)-FMOC(K,JAG))*VAL*SCAL
                END IF
              END DO
            END IF
          END DO
        END DO
      END IF
C
      RETURN
C
      END SUBROUTINE DFTB_ESDIM_DC
C
C-----------------------------------------------------------------------
C*MODULE DFTBFO    *DECK DFTB_ESDIM2
C>
C>    @brief ES-DIM for FMO-DFTB
C>
C>    @details Calculate electrostatic dimer (ES-DIM) contribution to
C>             energy at FMO-DFTB.  This is actual subroutine.
C>
C>    @author Yoshio Nishimoto
C>    - Sep, 2014- Subroutine written
C>    @date   Feb, 2016 - Yoshio Nishimoto
C>    - Performance improvement
C>
C>           --- INPUT ---
C>    @param IFG,JFG  I and J fragment index
C>    @param MODGAMMA How gamma values are calculated
C>    @param NATFRG   Number of atoms in each fragment
C>    @param IAGLOB   Global atom index
C>    @param INDFRG   First atom of each fragment in terms of NATFMOB?
C>    @param IATFRG   Global atom index with NBDFG?
C>    @param MAXNAT   Max atom
C>    @param ISPEFMO  ISPE for all atoms (NATFMO)
C>    @param NSPE     Number of species in the system
C>    @param GAMMA    DFTB2 gamma values
C>    @param ZREF     Nuclear charges of a given atom.  In FMO-DFTB, ZREF
C>                    of BDA or BAA is differet from QREF
C>    @param POPMAT   Mulliken population
C>    @param FMOC     Coordinates of allatoms (NATFMO)
C>    @parma DFTB3    Whether DFTB3 or not
C>    @param HUBDER   Hubbard derivatives
C>    @param GAMMA3   DFTB3 gamma values
C>    @param NATFMO   Number of atoms in the system
C>    @param HUBB     Atom-resolved Hubbard values
C>    @param DAMPXH0  Whether X-H damping is used or not
C>    @param DAMPXHE  Exponent of X-H damping
C>    @param SPE      Name of species
C>    @param NATFMOB  NATFMO+NBDFG
C>    @param NN
C>    @param NN2
C>           --- IN/OUTPUT ---
C>    @param ESDIM    ES-DIM contribution
C>    @param CHAMULI  Mulliken population of IFG
C>    @param CHAMULJ  Mulliken population of IFG
C>    @param DISTMAT  Distance matrix between IFG-JFG
C>
C
      SUBROUTINE DFTB_ESDIM2(IFG,JFG,MODGAMMA,NATFRG,INDFRG,
     *  IATFRG,MAXNAT,ISPEFMO,NSPE,GAMMA,ZREF,POPMAT,FMOC,ESDIM,
     *  DFTB3,HUBDER,GAMMA3,NATFMO,HUBB,DAMPXH0,DAMPXHE,SPE,NATFMOB,
     *  NN,NN2,CHAMULI,CHAMULJ,DISTMAT)
C
      USE DFTBPB_MOD,ONLY: PERIOD,DFTB_GAMFMO_PBC
C
      IMPLICIT NONE
C
      DOUBLE PRECISION, PARAMETER :: ZERO = 0.0D+00
C
      INTEGER, INTENT(IN) :: IFG,JFG,NATFRG(*),INDFRG(*),
     *  IATFRG(*),MAXNAT,MODGAMMA,ISPEFMO(*),NSPE,NATFMO,NATFMOB,NN,NN2
      DOUBLE PRECISION, INTENT(IN) :: GAMMA(*),ZREF(*),
     *  POPMAT(MAXNAT,*),FMOC(3,*),HUBDER(*),GAMMA3(NN2,NN2,*),HUBB(*),
     *  DAMPXHE,SPE(*)
      DOUBLE PRECISION, INTENT(INOUT) :: ESDIM,CHAMULI(*),CHAMULJ(*),
     *  DISTMAT(*)
      LOGICAL,INTENT(IN) :: DFTB3,DAMPXH0
C
      INTEGER :: NATI,NATJ,INDI,INDJ,IAG,JAG,I,J,NSEQ,ISP,JSP,II
      DOUBLE PRECISION :: CHAMUL1, CHAMUL2, GAMTMP, TMP(12,NSPE,NSPE),
     * DFTB_GAM2,DIST,TMP3(10,NSPE,NSPE),GAMI,GAMJ,ESDIM3,DISTI,XX,YY,
     * ZZ,ONE_THIRD
      DOUBLE PRECISION :: DDOT
      DOUBLE PRECISION :: XDOT
      LOGICAL DAMPXH
C
      NATI = NATFRG(IFG)
      NATJ = NATFRG(JFG)
      INDI = INDFRG(IFG)
      INDJ = INDFRG(JFG)
      ESDIM3=ZERO
C
      CALL DCOPY(NATI,POPMAT(1,IFG),1,CHAMULI,1)
      CALL VSUB(ZREF(INDI),1,CHAMULI,1,CHAMULI,1,NATI)
      CALL DCOPY(NATJ,POPMAT(1,JFG),1,CHAMULJ,1)
      CALL VSUB(ZREF(INDJ),1,CHAMULJ,1,CHAMULJ,1,NATJ)
      IF (DFTB3) THEN
        CALL DCOPY(NATI,CHAMULI,1,CHAMULI(NATI+1),1)
        DO I = 1, NATI
          CHAMULI(NATI+I) = CHAMULI(NATI+I)*CHAMULI(NATI+I)
        END DO
        CALL DCOPY(NATJ,CHAMULJ,1,CHAMULJ(NATJ+1),1)
        DO J = 1, NATJ
          CHAMULJ(NATJ+J) = CHAMULJ(NATJ+J)*CHAMULJ(NATJ+J)
        END DO
      END IF
C
      IF (IAND(MODGAMMA,3).EQ.0) THEN
        CALL DCOPY(12*NSPE*NSPE,GAMMA,1,TMP,1)
        IF (DFTB3) CALL DCOPY(10*NSPE*NSPE,GAMMA3,1,TMP3,1)
        NSEQ = 0
        dist = 0.0d+00
        ONE_THIRD = 1.0D+00/3.0D+00
        DO I = 1, NATI
          CHAMUL1 = CHAMULI(I)
C         IF (CHAMUL1.EQ.ZERO) CYCLE
          IAG = IATFRG(INDI+I-1)
          ISP = ISPEFMO(IAG)
          DO J = 1, NATJ
            NSEQ = NSEQ + 1
            CHAMUL2 = CHAMULJ(J)
C           IF (CHAMUL2.EQ.ZERO) CYCLE
            JAG = IATFRG(INDJ+J-1)
            JSP = ISPEFMO(JAG)
            DAMPXH = DAMPXH0
            CALL DFTB_DMPCHK(ISP,JSP,DAMPXH)
            IF (PERIOD) THEN
              CALL DFTB_GAMFMO_PBC(IAG,JAG,FMOC(1,IAG),FMOC(1,JAG),
     *                             GAMTMP,GAMI,GAMJ,TMP(1,ISP,JSP),
     *                             TMP3(1,ISP,JSP),TMP3(1,JSP,ISP),
     *                             DAMPXH,ISP,JSP)
            ELSE
              XX = FMOC(1,IAG)-FMOC(1,JAG)
              YY = FMOC(2,IAG)-FMOC(2,JAG)
              ZZ = FMOC(3,IAG)-FMOC(3,JAG)
              DIST = XX*XX + YY*YY + ZZ*ZZ
              CALL DFTB_GAM23(0,DIST,GAMTMP,GAMI,GAMJ,
     *             TMP(1,ISP,JSP),TMP3(1,ISP,JSP),TMP3(1,JSP,ISP),
     *             DAMPXH,DAMPXHE,SPE(ISP),SPE(JSP),HUBDER(ISP),
     *             HUBDER(JSP),DFTB3)
            END IF
            IF (DFTB3) 
     *        GAMTMP = GAMTMP + ONE_THIRD*(CHAMUL1*GAMI+CHAMUL2*GAMJ)
c           write(6,*) 'www3es',i,j,CHAMUL1,CHAMUL2,GAMI,GAMJ

            ESDIM = ESDIM + GAMTMP*CHAMUL1*CHAMUL2
            DISTMAT(NSEQ) = DIST
          END DO
        END DO
      ELSE IF (IAND(MODGAMMA,1).NE.0) THEN
        DO I = 1, NATI
          CHAMUL1 = CHAMULI(I)
          IF (CHAMUL1.EQ.ZERO) CYCLE
          IAG = IATFRG(INDI+I-1)
          DO J = 1, NATJ
            CHAMUL2 = CHAMULJ(J)
            IF (CHAMUL2.EQ.ZERO) CYCLE
            JAG = IATFRG(INDJ+J-1)
            CALL DFTB_CNVSQ(IAG,JAG,NSEQ)
            GAMTMP = GAMMA(NSEQ)
            ESDIM = ESDIM + GAMTMP*CHAMUL1*CHAMUL2
            IF (DFTB3) THEN
              GAMI = GAMMA3(IAG,JAG,1) !! G_ac
              GAMJ = GAMMA3(JAG,IAG,1) !! G_ca
              ESDIM3 = ESDIM3
     *          + CHAMUL1*CHAMUL1*CHAMUL2*GAMI
     *          + CHAMUL1*CHAMUL2*CHAMUL2*GAMJ
            END IF
          END DO
        END DO
      ELSE IF (IAND(MODGAMMA,2).NE.0) THEN
        DO I = 1, NATI
          CHAMUL1 = CHAMULI(I)
          IF (CHAMUL1.EQ.ZERO) CYCLE
          NSEQ = INDJ + NATFMOB*(INDI+I-2)
          ESDIM = ESDIM + CHAMUL1*DDOT(NATJ,CHAMULJ,1,GAMMA(NSEQ),1)
          IF (DFTB3) THEN
            II = INDI+I-1
            ESDIM3 = ESDIM3
     *        + CHAMUL1*CHAMUL1*DDOT(NATJ,CHAMULJ,1,GAMMA3(INDJ,II,2),1)
            ESDIM3 = ESDIM3
     *        + CHAMUL1*DDOT(NATJ,CHAMULJ(NATJ+1),1,GAMMA3(INDJ,II,1),1)
          END IF
        END DO
      END IF
C
      IF (DFTB3) ESDIM = ESDIM + ESDIM3/3.0D+00
C
      RETURN
C
      END SUBROUTINE DFTB_ESDIM2
C
C-----------------------------------------------------------------------
C*MODULE DFTBFO    *DECK DFTB_ESDIM_DC_HELP
C>
C>    @brief ES-DIM for FMO-DFTB
C>
C>    @details Helper of dispersion correction during ES-DIM calculation
C>
C>    @author Yoshio Nishimoto
C>    - Sep, 2014- Subroutine written
C>
C>           --- INPUT ---
C>    @param IFG,JFG I and J fragment index
C>    @param NDER    Order of derivative to be calculated
C>    @param DISTMAT Distance matrix between IFG-JFG
C>
C
      SUBROUTINE DFTB_ESDIM_DC_HELP(IFG,JFG,NDER,DISTMAT)
      use mx_limits, only: mxatm,mxrt
      USE DFTBPB_MOD,ONLY: PERIOD,DFTB_ESDIM_DCPBC
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      INTEGER,PARAMETER :: MXSPE=10
C
      INTEGER :: IFG,JFG,NDER
      DOUBLE PRECISION :: DISTMAT(*)
C
      COMMON /DFTBPR/ ETEMP,DFTBDP(MXSPE*14),DAMPXHE,HUBDER(MXSPE),
     *                ZREF(MXATM),SPNCST(6,MXSPE),SPE(MXATM),NSPE,
     *                MAXANG(MXATM),ISPE(MXATM),IND(MXATM+1),IDFTBD,
     *                PARAMDIR
      COMMON /ENRGYS/ ENUCR,EELCT,ETOT,SZ,SZZ,ECORE,ESCF,EERD,E1,E2,
     *                VEN,VEE,EPOT,EKIN,ESTATE(MXRT),STATN,EDFT(3),EDISP
      COMMON /FMCOM / X(1)
      COMMON /FMODTB/ LGAMMA,LZREF,LISPE,LSHIFTG,LSHIFTCT,
     *                LGAMMAD,MODESD,MODGAMMA,LDCTMP,LGAMMA3,LGAMMA3D,
     *                LNBOND,LCTMUL,LESPDFTB,LESPASC,NDFTBRST,LDFTBRST
      COMMON /FMOINF/ NFG,NLAYER,NATFMO,NBDFG,NAOTYP,NBODY,NSEGM
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
C
      IF (IDFTBD.EQ.3) RETURN
C
        if (nbody.eq.3) then
          iifg=max(ifg,jfg)
          jjfg=min(ifg,jfg)
          ijfg=(iifg*iifg-3*iifg)/2+jjfg+1
          ida=1
          ida=ida-ixftch(x(lndtfrg),ijfg)
          SCAL = DBLE(IDA)
        ELSE
          SCAL = 1.0D+00
        END IF
C
      IF (PERIOD) THEN
        CALL DFTB_ESDIM_DCPBC(IFG,JFG,NDER,MODGAMMA,X(LDCTMP),
     *                        X(LFMODE),EDISP,SCAL)
      ELSE
        CALL DFTB_ESDIM_DC(IFG,JFG,NDER,NATFMO,MODGAMMA,X(LNATFRG),
     *    X(LINDFRG),X(LIATFRG),X(LISPE),X(LDCTMP),
     *    X(LFMOC),X(LFMODE),X(LZREF),DISTMAT,X(LNBOND),EDISP,SCAL,
     *    NBODY)
      END IF
C
      RETURN
C
      END SUBROUTINE DFTB_ESDIM_DC_HELP
C
C-----------------------------------------------------------------------
C*MODULE DFTBFO    *DECK DFTB_ESDIM3
C>
C>    @brief ES-DIM for FMO-DFTB
C>
C>    @details Helper of dispersion correction during ES-DIM calculation
C>             This subroutine is used when MODESD=1 or 2.  In terms of
C>             algorithm, the calculation is similar to that of shift
C>             contribution in normal DFTB calculation.  This only
C>             calculates SHIFT(A)=\sum_I \sum_B \gamma_{AB} \Delta q_B,
C>             and in the end by calling DFTB_ESDIM5 with IFG=0, actual
C>             energy contribution, \sum_A \Delta q_A * shift(A) will be
C>             calculated.
C>
C>    @author Yoshio Nishimoto
C>    - Sep, 2014- Subroutine written
C>
C>           --- INPUT ---
C>    @param IFG,JFG I and J fragment index
C>    @param NDER    Order of derivative to be calculated
C>    @param CHAMULI Mulliken population of IFG
C>    @param CHAMULJ Mulliken population of IFG
C>    @param DIST    Distance matrix between IFG-JFG
C>           --- OUTPUT ---
C>    @param ESDIMTMP ES-DIM contribution
C>
C
      SUBROUTINE DFTB_ESDIM3(IFG,JFG,NDER,ESDIMTMP,CHAMULI,CHAMULJ,DIST,
     *                       MDOUT)
      use mx_limits, only: mxatm,mxrt
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      INTEGER,PARAMETER :: MXSPE=10
C
      INTEGER :: IFG, JFG
      DOUBLE PRECISION :: ESDIMTMP(*),CHAMULI(*),CHAMULJ(*),DIST(*)
      LOGICAL :: DC,DCCHG,DCABC,MDOUT
C
      COMMON /DFTB  / DFTBFL,SCC,SRSCC,DFTB3,DAMPXH,LCDFTB
      LOGICAL DFTBFL,SCC,SRSCC,DFTB3,DAMPXH,LCDFTB
      COMMON /DFTBPR/ ETEMP,DFTBDP(MXSPE*14),DAMPXHE,HUBDER(MXSPE),
     *                ZREF(MXATM),SPNCST(6,MXSPE),SPE(MXATM),NSPE,
     *                MAXANG(MXATM),ISPE(MXATM),IND(MXATM+1),IDFTBD,
     *                PARAMDIR
      COMMON /DFTBSK/ SKDIM(MXSPE,MXSPE),SKSPIN(MXSPE),QREFL(3,MXSPE),
     *                HUBBL(3,MXSPE),QREF(MXSPE),HUBB(MXSPE),SKCUT2,
     *                NEEDSK,LSKHTAB(MXSPE,MXSPE),LSKSTAB(MXSPE,MXSPE),
     *                LSKGRID(MXSPE,MXSPE),LSKSELF(MXSPE)
      COMMON /DFTDC / DC,IDCVER,DCCHG,DCABC,DCSR,DCS6,DCS8,DCALP,DCA1,
     *                DCA2
      COMMON /ENRGYS/ ENUCR,EELCT,ETOT,SZ,SZZ,ECORE,ESCF,EERD,E1,E2,
     *                VEN,VEE,EPOT,EKIN,ESTATE(MXRT),STATN,EDFT(3),EDISP
      COMMON /FMCOM / X(1)
      COMMON /FMODTB/ LGAMMA,LZREF,LISPE,LSHIFTG,LSHIFTCT,
     *                LGAMMAD,MODESD,MODGAMMA,LDCTMP,LGAMMA3,LGAMMA3D,
     *                LNBOND,LCTMUL,LESPDFTB,LESPASC,NDFTBRST,LDFTBRST
      common /fmodim/ maxbnd,maxknd,maxcbs,maxcao,maxvec,maxl1,maxnat,
     *                maxabd,maxbas,maxbbd,maxlmo,maxslo,maxabd2,maxrij
      COMMON /FMOINF/ NFG,NLAYER,NATFMO,NBDFG,NAOTYP,NBODY,NSEGM
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
c     common /fmorun/ espscf,e0scf(2),emp2s,IDAFMO,icurfg,jcurfg,kcurfg,
c    *                icurlay,icurunt,nat1e,ncursh,ngau,icurpop,ifmostp,
c    *                moncor,needr,modrst,norbproj,nunesp,iskipesp,
c    *                IESDPPC,idoprop,mp2run,icurit,idmfmo,iddfmo,
c    *                iddcur,nddleft,ivmfmo,nzmtfmo,ifmobas,itmfmo(2)
      common /GRAD  / DE(3,MXATM)
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      LOGICAL GOPARR,DSKWRK,MASWRK
C
C     ----- NOTHING TO BE DONE FOR NCC CALCULATION -----
C
      IF (.NOT.SCC) THEN
        IF (DC) THEN
C         CALL DFTB_ESDIM_DC(IFG,JFG,NDER,NATFMO,MODGAMMA,X(LNATFRG),
C    *    X(LINDFRG),X(LIATFRG),X(LISPE),X(LDCTMP),
C    *    X(LFMOC),X(LFMODE),X(LZREF),EDISP)
          IF (IFG.EQ.0) THEN
            WRITE (IW,'(  2X,"TOTAL EDISP      = ",F15.10,/)') EDISP
            IF(MDOUT) THEN
            IF(MASWRK)
     *      WRITE(IW,FMT='('' ...... END OF LUMP ESD ......'')')
            CALL TIMIT(1)
            END IF
          END IF
        END IF
        RETURN
      END IF
C
c     LPOPMATI=LPOPMAT+(MAXNAT*NFG)*(ICURPOP-1)
      IF (IFG.EQ.0) THEN
        CALL DFTB_ESDIM5(JFG,X(LNATFRG),X(LINDFRG),
     *         MAXNAT,X(LPOPMAT),X(LZREF),ESDIMTMP)
c    *         MAXNAT,X(LPOPMATI),X(LZREF),ESDIMTMP)
c       IF (IAND(NPRFMO,3).NE.3) THEN
        if (MDOUT) THEN
          IF (DC) THEN
            WRITE (IW,'(/,2X,"TOTAL ESD2 ENERGY= ",F15.10)') ESDIMTMP(1)
            WRITE (IW,'(2X,"TOTAL EDISP      = ",F15.10,/)') EDISP
          ELSE
          WRITE (IW,'(/,2X,"TOTAL ESD2 ENERGY= ",F15.10,/)') ESDIMTMP(1)
          END IF
c         IF(MDOUT) THEN
          IF(MASWRK)
     *      WRITE(IW,FMT='('' ...... END OF LUMP ESD ......'')')
          CALL TIMIT(1)
c         END IF
        END IF
        RETURN
      END IF
C
      IF (IAND(MODGAMMA,2).NE.0) THEN
        NN = 2
        NN2 = NATFMO+NBDFG
      ELSE
        NN = 1
        NN2 = NATFMO
      END IF
c    * X(LINDFRG),X(LIATFRG),MAXNAT,X(LISPE),NSPE,X(LPOPMATI),X(LGAMMA),
      CALL DFTB_ESDIM4(IFG,JFG,MODGAMMA,X(LNATFRG),
     * X(LINDFRG),X(LIATFRG),MAXNAT,X(LISPE),NSPE,X(LPOPMAT),X(LGAMMA),
     * X(LZREF),X(LFMOC),ESDIMTMP,DFTB3,HUBDER,X(LGAMMA3),NATFMO,HUBB,
     * DAMPXH,DAMPXHE,SPE,NATFMO+NBDFG,NN,NN2,CHAMULI,CHAMULJ,DIST)
C
C     IF (DC) THEN
C       CALL DFTB_ESDIM_DC(IFG,JFG,NDER,NATFMO,MODGAMMA,X(LNATFRG),
C    *  X(LINDFRG),X(LIATFRG),X(LISPE),X(LDCTMP),
C    *  X(LFMOC),X(LFMODE),X(LZREF),EDISP)
C     END IF
C
      IF (NDER.GT.0) THEN
        CALL VCLR(DE,1,3*NAT)
        if (nbody.eq.3) then
          iifg=max(ifg,jfg)
          jjfg=min(ifg,jfg)
          ijfg=(iifg*iifg-3*iifg)/2+jjfg+1
          ida=1
          ida=ida-ixftch(x(lndtfrg),ijfg)
          SCAL = DBLE(IDA)
        ELSE
          SCAL = 1.0D+00
        END IF
      CALL DFTB_ESDGRAD(IFG,JFG,X(LISPE),
     *  X(LINDFRG),X(LNATFRG),X(LIATFRG),X(LGAMMAD),
     *  X(LSHIFTG+NATFMO+NBDFG),X(LFMOC),X(LFMODE),X(LGAMMA3D),CHAMULI,
     *  CHAMULJ,X(LGAMMA3),DIST,NN,NN2,SCAL)
      END IF
C
      RETURN
C
      END SUBROUTINE DFTB_ESDIM3
C
C-----------------------------------------------------------------------
C*MODULE DFTBFO    *DECK DFTB_ESDIM4
C>
C>    @brief ES-DIM for FMO-DFTB
C>
C>    @details Calculate electrostatic dimer (ES-DIM) contribution to
C>             energy at FMO-DFTB.  This calculates only shift-like
C>             terms.
C>
C>    @author Yoshio Nishimoto
C>    - Sep, 2014- Subroutine written
C>    @date   Feb, 2016 - Yoshio Nishimoto
C>    - Performance improvement
C>
C>           --- INPUT ---
C>    @param IFG,JFG I and J fragment index
C>    @param MODGAMMA How gamma values are calculated
C>    @param NATFRG Number of atoms in each fragment
C>    @param INDFRG First atom of each fragment in terms of NATFMOB?
C>    @param IATFRG Global atom index with NBDFG?
C>    @param MAXNAT Max atom
C>    @param ISPEFMO ISPE for all atoms (NATFMO)
C>    @param NSPE Number of species in the system
C>    @param POPMAT Mulliken population
C>    @param GAMMA DFTB2 gamma values
C>    @param ZREF Nuclear charges of a given atom.  In FMO-DFTB, ZREF
C>           of BDA or BAA is differet from QREF
C>    @param FMOC Coordinates of allatoms (NATFMO)
C>    @parma DFTB3 Whether DFTB3 or not
C>    @param HUBDER Hubbard derivatives
C>    @param GAMMA3 DFTB3 gamma values
C>    @param NATFMO Number of atoms in the system
C>    @param HUBB Atom-resolved Hubbard values
C>    @param DAMPXH0 Whether X-H damping is used or not
C>    @param DAMPXHE Exponent of X-H damping
C>    @param SPE Name of species
C>    @param NATFMOB NATFMO+NBDFG
C>    @param NN
C>    @param NN2
C>           --- IN/OUTPUT ---
C>    @param ESDIM ES-DIM contribution
C>    @param CHAMULI Mulliken population of IFG
C>    @param CHAMULJ Mulliken population of JFG
C>
C
      SUBROUTINE DFTB_ESDIM4(IFG,JFG,MODGAMMA,NATFRG,INDFRG,
     *  IATFRG,MAXNAT,ISPEFMO,NSPE,POPMAT,GAMMA,ZREF,FMOC,ESDIMTMP,
     *  DFTB3,HUBDER,GAMMA3,NATFMO,HUBB,DAMPXH0,DAMPXHE,SPE,NATFMOB,
     *  NN,NN2,CHAMULI,CHAMULJ,DISTMAT)
C
      IMPLICIT NONE
C
      DOUBLE PRECISION, PARAMETER :: ZERO = 0.0D+00
C
      INTEGER, INTENT(IN) :: IFG,JFG,NATFRG(*),INDFRG(*),
     *  IATFRG(*),MAXNAT,MODGAMMA,ISPEFMO(*),NSPE,NATFMO,NATFMOB,NN,NN2
      DOUBLE PRECISION, INTENT(IN) :: POPMAT(MAXNAT,*),GAMMA(*),
     *  ZREF(*),FMOC(3,*),HUBDER(*),GAMMA3(NN2,NN2,*),HUBB(*),DAMPXHE,
     *  SPE(*)
      DOUBLE PRECISION, INTENT(INOUT) :: ESDIMTMP(*),CHAMULI(*),
     *  CHAMULJ(*),DISTMAT(*)
      LOGICAL,INTENT(IN) :: DFTB3,DAMPXH0
C
      INTEGER :: NATI,NATJ,INDI,INDJ,IAG,JAG,I,J,NSEQ,ISP,JSP,II
      DOUBLE PRECISION :: CHAMUL1,CHAMUL2,GAMTMP,TMPSH,DFTB_GAM2,DIST,
     *  TMP(12,NSPE,NSPE),TMP3(10,NSPE,NSPE),GAMI,GAMJ,TMPSH3,DISTI,
     *  XX,YY,ZZ
      DOUBLE PRECISION :: DDOT
      DOUBLE PRECISION :: XDOT
      LOGICAL DAMPXH
C
C     ALWAYS IFG > JFG
      NATI = NATFRG(IFG)
      NATJ = NATFRG(JFG)
      INDI = INDFRG(IFG)
      INDJ = INDFRG(JFG)
C
      CALL DCOPY(NATI,POPMAT(1,IFG),1,CHAMULI,1)
      CALL VSUB(ZREF(INDI),1,CHAMULI,1,CHAMULI,1,NATI)
      CALL DCOPY(NATJ,POPMAT(1,JFG),1,CHAMULJ,1)
      CALL VSUB(ZREF(INDJ),1,CHAMULJ,1,CHAMULJ,1,NATJ)
      IF (DFTB3) THEN
        CALL DCOPY(NATI,CHAMULI,1,CHAMULI(NATI+1),1)
        DO I = 1, NATI
          CHAMULI(NATI+I) = CHAMULI(NATI+I)*CHAMULI(NATI+I)
        END DO
        CALL DCOPY(NATJ,CHAMULJ,1,CHAMULJ(NATJ+1),1)
        DO J = 1, NATJ
          CHAMULJ(NATJ+J) = CHAMULJ(NATJ+J)*CHAMULJ(NATJ+J)
        END DO
      END IF
C
      IF (IAND(MODGAMMA,3).EQ.0) THEN
        CALL DCOPY(12*NSPE*NSPE,GAMMA,1,TMP,1)
        IF (DFTB3) CALL DCOPY(10*NSPE*NSPE,GAMMA3,1,TMP3,1)
        NSEQ = 0
        DO I = 1, NATI
          TMPSH = ZERO
          TMPSH3 = ZERO
          IAG = IATFRG(INDI+I-1)
          ISP = ISPEFMO(IAG)
          DO J = 1, NATJ
            CHAMUL2 = CHAMULJ(J)
C           IF (CHAMUL2.EQ.ZERO) CYCLE
            JAG = IATFRG(INDJ+J-1)
            JSP = ISPEFMO(JAG)
            XX = FMOC(1,IAG)-FMOC(1,JAG)
            YY = FMOC(2,IAG)-FMOC(2,JAG)
            ZZ = FMOC(3,IAG)-FMOC(3,JAG)
            DIST = XX*XX + YY*YY + ZZ*ZZ
            DAMPXH = DAMPXH0
            CALL DFTB_DMPCHK(ISP,JSP,DAMPXH)
            IF (DFTB3) THEN
              CHAMUL1 = CHAMULI(I)
              CALL DFTB_GAM23(0,DIST,GAMTMP,GAMI,GAMJ,
     *             TMP(1,ISP,JSP),TMP3(1,ISP,JSP),TMP3(1,JSP,ISP),
     *             DAMPXH,DAMPXHE,SPE(ISP),SPE(JSP),HUBDER(ISP),
     *             HUBDER(JSP),DFTB3)
              TMPSH3 = TMPSH3 + CHAMUL1*CHAMUL2*GAMI
     *                        + CHAMUL2*CHAMUL2*GAMJ
            ELSE
              DIST = SQRT(DIST)
              DISTI = 0
              IF(DIST.NE.0) DISTI = 1.0D+00/DIST
              GAMTMP = DFTB_GAM2(DIST,DISTI,TMP(1,ISP,JSP),DAMPXH,
     *          DAMPXHE,SPE(ISP),SPE(JSP))
            END IF
            TMPSH = TMPSH + GAMTMP * CHAMUL2
            NSEQ = NSEQ + 1
            DISTMAT(NSEQ) = DIST
          END DO
          TMPSH = TMPSH + TMPSH3/3.0D+00
          ESDIMTMP(INDI+I-1) = ESDIMTMP(INDI+I-1) + TMPSH
        END DO
      ELSE IF (IAND(MODGAMMA,1).NE.0) THEN
        DO I = 1, NATI
          TMPSH = ZERO
          TMPSH3 = ZERO
          IAG = IATFRG(INDI+I-1)
          DO J = 1, NATJ
            CHAMUL2 = CHAMULJ(J)
            IF (CHAMUL2.EQ.ZERO) CYCLE
            JAG = IATFRG(INDJ+J-1)
            CALL DFTB_CNVSQ(IAG,JAG,NSEQ)
            GAMTMP = GAMMA(NSEQ)
            TMPSH = TMPSH + GAMTMP * CHAMUL2
            IF (DFTB3) THEN
              CHAMUL1 = CHAMULI(I)
              GAMI = GAMMA3(IAG,JAG,1) !! G_ac
              GAMJ = GAMMA3(JAG,IAG,1) !! G_ca
              TMPSH3 = TMPSH3 + CHAMUL1*CHAMUL2*GAMI
     *                        + CHAMUL2*CHAMUL2*GAMJ
            END IF
          END DO
          TMPSH = TMPSH + TMPSH3/3.0D+00
          ESDIMTMP(INDI+I-1) = ESDIMTMP(INDI+I-1) + TMPSH
        END DO
      ELSE IF (IAND(MODGAMMA,2).NE.0) THEN
        DO I = 1, NATI
          TMPSH = ZERO
          TMPSH3 = ZERO
          NSEQ = INDJ + NATFMOB*(INDI+I-2)
          TMPSH = TMPSH + DDOT(NATJ,CHAMULJ,1,GAMMA(NSEQ),1)
          IF (DFTB3) THEN
            CHAMUL1 = CHAMULI(I)
            II = INDI+I-1
            TMPSH3 = TMPSH3
     *        + CHAMUL1*DDOT(NATJ,CHAMULJ,1,GAMMA3(INDJ,II,2),1)
     *        + DDOT(NATJ,CHAMULJ(NATJ+1),1,GAMMA3(INDJ,II,1),1)
          END IF
          TMPSH = TMPSH + TMPSH3/3.0D+00
          ESDIMTMP(INDI+I-1) = ESDIMTMP(INDI+I-1) + TMPSH
        END DO
      END IF
C
      RETURN
C
      END SUBROUTINE DFTB_ESDIM4
C
C-----------------------------------------------------------------------
C*MODULE DFTBFO    *DECK DFTB_ESDIM5
C>
C>    @brief ES-DIM for FMO-DFTB
C>
C>    @details Calculate electrostatic dimer (ES-DIM) contribution to
C>             energy at FMO-DFTB.  This should be calculated in the
C>             end of ES-DIM calculation.
C>
C>    @author Yoshio Nishimoto
C>    - Sep, 2014- Subroutine written
C>
C>           --- INPUT ---
C>    @param NFG Number of fragments
C>    @param NATFRG Number of atoms in each fragment
C>    @param INDFRG First atom of each fragment in terms of NATFMOB?
C>    @param MAXNAT Max atom
C>    @param POPMAT Mulliken population
C>    @param ZREF Nuclear charges of a given atom.  In FMO-DFTB, ZREF
C>           of BDA or BAA is differet from QREF
C>           --- IN/OUTPUT ---
C>    @param ESDIM ES-DIM contribution.  On entry, this has
C>           atom-resolve shift contribution, and on exit, the first
C>           element will contain the sum of all ES-DIM contribution.
C>
C
      SUBROUTINE DFTB_ESDIM5(NFG,NATFRG,INDFRG,MAXNAT,
     *  POPMAT,ZREF,ESDIMTMP)
C
      IMPLICIT NONE
C
      DOUBLE PRECISION, PARAMETER :: ZERO = 0.0D+00
C
      INTEGER :: NFG, NATFRG(*), INDFRG(*), MAXNAT
      DOUBLE PRECISION :: POPMAT(MAXNAT,*), ZREF(*), ESDIMTMP(*)
C
      INTEGER :: NATI, INDI, I, LFG
      DOUBLE PRECISION :: CHAMUL1, TMP
C
      TMP = ZERO
      DO LFG = 1, NFG
        NATI = NATFRG(LFG)
        INDI = INDFRG(LFG)
        DO I = 1, NATI
          IF (ESDIMTMP(INDI+I-1).EQ.ZERO) CYCLE
          CHAMUL1 = POPMAT(I,LFG) - ZREF(INDI+I-1)
          TMP = TMP + ESDIMTMP(INDI+I-1) * CHAMUL1
C         write (*,*) esdimtmp(indi+i-1),chamul1
        END DO
      END DO
      ESDIMTMP(1) = TMP
C
      RETURN
C
      END SUBROUTINE DFTB_ESDIM5
C
C-----------------------------------------------------------------------
C*MODULE DFTBFO    *DECK DFTB_FMOGRAD
C>
C>    @brief Gradient in FMO-DFTB
C>
C>    @details Helper of gradient calculation in FMO-DFTB.  Actual
C>             calculation is done in DFTB_FMOGRAD2.
C>
C>    @author Yoshio Nishimoto
C>    - Sep, 2014- Subroutine written
C>    @date   Hiroya Nakata
C>    - Feb, 2016- Added FMO-DFTB Hessian
C>
C>           --- INPUT ---
C>    @param LAYFRG Layer of each fragment
C>    @param FMODE Cartesian geometrical gradient of monomers and dimers
C>    @param IDMREC0
C>    @param NUMFRG Nuber of basis functions in each fragment
C>    @param IODFMO 
C>    @param ORBXCH Whether eigenvectors are exchanged or not
C>    @param FCM    Force constant
C>    @param MDOUT  Print option
C>
C
      SUBROUTINE DFTB_FMOGRAD(LAYFRG,FMODE,IDMREC0,NUMFRG,IODFMO,ORBXCH,
     *           FCM,MDOUT)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      
      DOUBLE PRECISION, PARAMETER :: ZERO = 0.0D+00
C
      COMMON /FMCOM / X(1)
      COMMON /FMOINF/ NFG,NLAYER,NATFMO,NBDFG,NAOTYP,NBODY,NSEGM
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
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      LOGICAL GOPARR,DSKWRK,MASWRK
      LOGICAL ORBXCH,MDOUT

      INTEGER :: LAYFRG(*),IDMREC0(*),IODFMO(*)
      DOUBLE PRECISION :: FMODE(3,NATFMO,*),FCM(*)
      COMMON /FMODTB/ LGAMMA,LZREF,LISPE,LSHIFTG,LSHIFTCT,
     *                LGAMMAD,MODESD,MODGAMMA,LDCTMP,LGAMMA3,LGAMMA3D,
     *                LNBOND,LCTMUL,LESPDFTB,LESPASC,NDFTBRST,LDFTBRST

      TIM0 = ZERO
C     TIM1 = ZERO
      CALL TSECND(TIM0)
C     TLEFTS = TIMLIM - TIM0
C
      CALL DFTB_FMOGRAD2(LAYFRG,X(LNATFRG),
     *  X(LIATFRG),X(LINDFRG),
     *  FMODE,IDMREC0,NUMFRG,IODFMO,X(LISPE),X(LFMOC),
     *  ORBXCH,FCM,X(LPOPMAT))
C
      IF(MASWRK.AND.MDOUT) THEN
        WRITE (IW,*)
        WRITE(IW,FMT='('' ...... END OF FMO-DFTB GRADIENT ......'')')
      END IF
      IF (MDOUT) CALL TIMIT(1)
C
      END SUBROUTINE DFTB_FMOGRAD
C
C-----------------------------------------------------------------------
C*MODULE DFTBFO    *DECK DFTB_FMOGRAD2
C>
C>    @brief Gradient in FMO-DFTB
C>
C>    @details Gradient calculation in FMO-DFTB.  It uses two shift
C>             contributions, X(LSHIFTG) and X(LSHIFTCT).  Some details
C>             can be found as comments below.
C>
C>    @author Yoshio Nishimoto
C>    - Sep, 2014- Subroutine written
C>    @date   Feb, 2016 - Yoshio Nishimoto
C>    - Added FMO-DFTB SCZV
C>    @date   Feb, 2016 - Hiroya Nakata
C>    - Added FMO-DFTB Hessian
C>    @date   Jun, 2016 - Hiroya Nakata
C>    - Bug fix for zero electron
C>    @date   Feb, 2018 - Yoshio Nishimoto
C>    - Added FMO-DFTB/AFO etc.
C>
C>           --- INPUT ---
C>    @param LAYFRG Layer of each fragment
C>    @param NATFRG Number of atoms in each fragment
C>    @param IATFRG Global atom index with NBDFG?
C>    @param INDFRG First atom of each fragment in terms of NATFMOB?
C>    @param IDMREC0
C>    @param NUMFRG Nuber of basis functions in each fragment
C>    @param IODFMO 
C>    @param ISPEFMO ISPE for all atoms (NATFMO)
C>    @param FMOC Coordinates for all NATFMO
C>    @param ORBXCH Whether eigenvectors are exchanged or not
C>    @param FCM    Force constant matrix
C>    @param POPMAT Mulliken population
C>           --- IN/OUTPUT ---
C>    @param FMODE Cartesian geometrical gradient of monomers and dimers
C>
C
      SUBROUTINE DFTB_FMOGRAD2(LAYFRG,NATFRG,IATFRG,INDFRG,
     * FMODE,IDMREC0,NUMFRG,IODFMO,ISPEFMO,FMOC,ORBXCH,FCM,
     * POPMAT)
      use mx_limits, only: mxatm,mxsh,mxao,mxgtot
      USE DFTBPB_MOD,ONLY: MAXINT,DFTB_INTPAIR
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      INTEGER,PARAMETER :: MXSPE=10
      DOUBLE PRECISION, PARAMETER:: ZERO=0.0D+00,ONE=1.0D+00,TWO=2.0D+00
C
      LOGICAL DFTBFL,SCC,SRSCC,DFTB3,DAMPXH,LCDFTB
      logical GOPARR,DSKWRK,MASWRK,ISGDDI,PAROUT,INITGDDI,wasgddi,MLGDDI
      integer ddi_world,ddi_group
      PARAMETER (ddi_world=0,ddi_group=1)
      logical myjob,MFRZ,SKIPOR
C
      COMMON /DFTB  / DFTBFL,SCC,SRSCC,DFTB3,DAMPXH,LCDFTB
      COMMON /FMCOM / X(1)
      common /fmodim/ maxbnd,maxknd,maxcbs,maxcao,maxvec,maxl1,maxnat,
     *                maxabd,maxbas,maxbbd,maxlmo,maxslo,maxabd2,maxrij
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
      INTEGER INDFRG(*),IATFRG(*),NATFRG(*),LAYFRG(*),IDMREC0(*),
     *  NUMFRG(*)
      DOUBLE PRECISION FMODE(3,NATFMO,*),FMOC(3,NATFMO),FCM(*)
      COMMON /DFTBPR/ ETEMP,DFTBDP(MXSPE*14),DAMPXHE,HUBDER(MXSPE),
     *                ZREF(MXATM),SPNCST(6,MXSPE),SPE(MXATM),NSPE,
     *                MAXANG(MXATM),ISPE(MXATM),IND(MXATM+1),IDFTBD,
     *                PARAMDIR
C
      COMMON /FMOAFO/ LINDFRZ,LGMK,LGMKSAV
      COMMON /FMOINF/ NFG,NLAYER,NATFMO,NBDFG,NAOTYP,NBODY,NSEGM
      COMMON /FMODTB/ LGAMMA,LZREF,LISPE,LSHIFTG,LSHIFTCT,
     *                LGAMMAD,MODESD,MODGAMMA,LDCTMP,LGAMMA3,LGAMMA3D,
     *                LNBOND,LCTMUL,LESPDFTB,LESPASC,NDFTBRST,LDFTBRST
      COMMON /FMOOPT/ ESPSCA(9),RESPAP(2),RESPPC(2),RESDIM,RESTRI(4),
     *                RCORSD,RESPCT,CONVFG,CNVDMP,COROFF,RFLMO(4),
     *                ORSHFT,ORSHFT2,CNVAFO,ASCREEN(4),IXESP,MXITFG,
     *                NGUESS,NBSSE,MODORB,MODPAR,IRSTSTP,IRSTLAY,NPRFMO,
     *                NFMOPAL,MODPRP,MAXL1C,IPIEDA,MODGRD,MODESP,IVMUL,
     *                MODLMO,NOPDEN,MOFOCK,modfd,modfmm,ncentm,ndualb,
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
      COMMON /MFRPAR/ MFRZ,NUMFRZ,NORFRZ,IFRZ(MXAO)
      COMMON /NSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),
     *                CF(MXGTOT),CG(MXGTOT),CH(MXGTOT),CI(MXGTOT),
     *                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),
     *                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
      COMMON /SCYVEC/ IDFYVEC,LIPTYVEC,LIPTYA,LFCMWRK
      COMMON /SCZVEC/ IDFZVEC,LIPTZVEC,LIPTLG,icursczv
      COMMON /SCZLAG/ LZVLAG,LZVWRK,LYALAG,LYAWRK,LFEQ1
C
      INTEGER :: ISPEFMO(*),IODFMO(*)
      LOGICAL ORBXCH
      DOUBLE PRECISION :: TMPD(10,NSPE,NSPE),
     * TMP(12,NSPE,NSPE),TMP3(10,NSPE,NSPE),TMPD3(11,NSPE,NSPE)
      DOUBLE PRECISION :: POPMAT(MAXNAT,*)
C
C     double precision detmp(3,natfmo)
C
C     THIS SUBROUTINE PERFORMS THE LOOP FOR ALL MONOMERS TO CALCULATE
C     GRADIENTS.  ALL DIMER GRADIENTS ARE CALCULATED SOMEWHERE.
C
      L1  = MAXL1
      L2  = (L1*L1+L1)/2
      L3  = L1 * L1
C
      CALL DERCHK(NDER)
C     ----- ALLOCATE MEMORY -----
C
      CALL VALFM(LOADFM)
      CALL GOTFM(NGOTMX)
      LWRK    = LOADFM  + 1
      LVEC    = LWRK    + L1
      LWRK1   = LVEC    + L3
      LWRK2   = LWRK1   + L3
      LOCC    = LWRK2   + L3
      LD      = LOCC    + L1
      LS      = LD      + L2
      LDS     = LS      + L2
      LPS     = LDS     + L3 !! D SQUARE
      LSHIFT  = LPS     + L2*MAXNAT*3
      LSHIFT2 = LSHIFT  + MAXNAT
      LAST    = LSHIFT2 + MAXNAT
C
      IF (RFLMO(1).NE.0) THEN
        LVECL = LAST
        LSCR  = LVECL + L3
        LE    = LSCR  + L1*8
        LIWRK = LE    + L1
        LAST  = LIWRK + L1
C
        LFAO  = LD
        LWRK3 = LDS
      END IF
      NEED    = LAST - LOADFM - 1
C
C     ----- GET MEMORY -----
C
      CALL GETFM(NEED)
C
C     ----- RESET GLOBAL COUNTER -----
C
      NATFMOB = NATFMO+NBDFG
      IF (ISGDDI) THEN
C       if(nproc.gt.1) then
C         call dscal((NATFMO+NBDFG)*2,1.0D+00/NPROC,X(LSHIFTG),1)
C         call dscal(NATFMO+NBDFG,1.0D+00/NPROC,X(LSHIFTCT),1)
C       endif
        CALL GDDI_SCOPE(DDI_WORLD)
        CALL DDI_GSUMF(2424,X(LSHIFTG),NATFMOB*2)
        CALL DDI_GSUMF(2424,X(LSHIFTCT),NATFMOB)
        IF (DFTB3) THEN
          CALL DDI_GSUMF(2424,X(LCTMUL),NATFMOB*2)
        ELSE
          CALL DDI_GSUMF(2424,X(LCTMUL),NATFMOB)
        END IF
        CALL GDDI_SCOPE(DDI_GROUP)
C       write(6,1234) (X(LSHIFTG+i-1),i=1,(NATFMO+NBDFG)*2)
C       write(6,1234) (X(LSHIFTCT+i-1),i=1,NATFMO+NBDFG)
C1234 format(6F12.8)
        CALL GDDICOUNT(-1,LGROUP,MYJOB)
      END IF
      if (iand(modgamma,2).ne.0) call vclr(x(lshiftg),1,natfmob*2)
C     if (maswrk) then
C     write (6,*) "shiftct"
C     do i = 1, natfmo+nbdfg
C       write (6,*) i,x(lshiftct+i-1)
C     end do
C     write (6,*) "shiftg(1)"
C     do i = 1, natfmo+nbdfg
C       write (6,*) i,x(lshiftg+i-1)
C     end do
C     write (6,*) "shiftg(2)"
C     do i = 1, natfmo+nbdfg
C       write (6,*) i,x(lshiftg+i-1+natfmob)
C     end do
C     write (6,*) "ctmul"
C     do i = 1, natfmo
C       write (6,*) i,x(lctmul+i-1)
C     end do
C     if (dftb3) then
C     do i = 1, natfmo
C       write (6,*) i,x(lctmul+i-1+natfmo)
C     end do
C     end if
C     write (6,*)
C     call valfm(loadfm)
C     lfmodewrk = loadfm    + 1
C     last      = lfmodewrk + natfmo*3
C     need1     = last - loadfm - 1
C     call getfm(need1)
C     call dcopy(natfmo*3,fmode(1,1,nbody),1,x(lfmodewrk),1)
C     if (isgddi) then
C       call gddi_scope(ddi_world)
C       call ddi_gsumf(2424,x(lfmodewrk),natfmo*3)
C       call gddi_scope(ddi_group)
C       call gddicount(-1,lgroup,myjob)
C     end if
C     write (6,*) "fmode"
C     do i = 1, natfmo
C       write (6,'(i3,3f20.10)') i,(fmode(k,i,2),k=1,3)
C     end do
C     call retfm(need1)
C     end if
C
C     ----- PREPARE FOR GAMMA DERIVATIVE -----
C
      IF (IAND(MODGAMMA,3).EQ.0) THEN
        CALL DCOPY(12*NSPE*NSPE,X(LGAMMA),1,TMP,1)
        IF (DFTB3) CALL DCOPY(10*NSPE*NSPE,X(LGAMMA3),1,TMP3,1)
      END IF
      CALL DCOPY(10*NSPE*NSPE,X(LGAMMAD),1,TMPD,1)
      IF (DFTB3) CALL DCOPY(11*NSPE*NSPE,X(LGAMMA3D),1,TMPD3,1)
C
      IF (IAND(MODGAMMA,2).NE.0) THEN
        NN = 2
        NN2 = NATFMO+NBDFG
      ELSE IF (IAND(MODGAMMA,1).NE.0) THEN
        NN = 1
        NN2 = NATFMO
      END IF
C
C     ----- START LOOP -----
C
      DO IFG = 1, NFG
        IF (ISGDDI) THEN
          call GDDICOUNT(0,LGROUP,MYJOB)
          IF (.NOT.MYJOB) CYCLE
        END IF
C
C     ----- RESTORE DENSITY MATRIX FOR MONOMER IFG -----
C
        CALL VCLR(X(LD),1,L2)
        ilay=layfrg(ifg)
        L1 = IAND(NUMFRG(IFG),65535)
        L2 = (L1*L1+L1)/2
        L3 = L1*L1
        NAT = NATFRG(IFG)
        NA = ISHFT(NUMFRG(IFG),-16)
        NB = NA
        irec0=idmrec0(ilay)
        IF (ORBXCH) THEN
          !! RESTORE VECTOR (PUT INTO X(LVEC)), THEN CONVERT TO
          !! DENSITY MATRIX
          CALL rareads(IDAFMO,iodfmo,X(LVEC),L3,irec0+ifg,0)
          CALL VCLR(X(LOCC),1,L1)
          IF (NA.EQ.NB) THEN
            DO I = 1, NA
              X(LOCC-1+I) = TWO
            END DO
          ELSE
            DO I = 1, NB
              X(LOCC-1+I) = TWO
            END DO
            DO I = NB+1, NA
              X(LOCC-1+I) = ONE
            END DO
          END IF
          CALL DFTB_DMTOCC(X(LD),X(LVEC),X(LOCC),X(LWRK2),L1,L1,1.0D-09)
        ELSE
          CALL RAREADS(IDAFMO,IODFMO,X(LD),L2,IREC0+IFG,0)
        END IF
C       write (*,*) "ifg = ",ifg
C       CALL PRTRIL(X(LD),L1L)
C
C     ----- PREPARE INFORMATION (INSTEAD OF CALLING MAKEMOL) -----
C
        INDI = INDFRG(IFG)
        IND(1) = 0
        DO I = 1, NAT
          IAG = IATFRG(INDI+I-1)
          CALL DCOPY(3,FMOC(1,IAG),1,C(1,I),1)
          ISPE(I) = ISPEFMO(IAG)
          IND(I+1) = IND(I) + MAXANG(ISPE(I))**2
          CALL IXSTOR(X(LIAGLOB),I,IAG)
        END DO
        IF (MAXINT.GT.0) CALL DFTB_INTPAIR(C,NAT,ISPE,MAXANG,NAOPBC)
C
C     ----- CALCULATE GAMMA DERIVATIVE TERMS FROM EMBEDDING CHARGES -----
C
        CALL DFTB_FMOGRAD_EMBED(NAT,NATFMO,NBDFG,ISPEFMO,IATFRG(INDI),
     *    X(LINDAT),X(LIALOC),X(LINDFRG),NN,NN2,INDI,X(LSHIFTCT+INDI-1),
     *    FMODE(1,1,NBODY),X(LCTMUL),FMOC,POPMAT(1,IFG),X(LZREF+INDI-1),
     *    TMP,TMP3,TMPD,TMPD3,X(LGAMMA),X(LGAMMA3),MASWRK,
     *    NATFRG,IATFRG,NFG,ich)
C
C     ----- MAKE SHIFT AND SHIFT2 ARRAYS -----
C
C     X(LSHIFTG ~ LSHIFTG+NATFMO+NBDFG) CONTAINS THE SUM OF CHARGE
C     CONTRIBUTIONS FROM ENVIRONMENTAL FRAGMENTS.  IN TERMS OF
C     EQUATIONS, IT IS USED FOR CALCULATING THE GRADIENTS FOR DIMER
C     IFG AND JFG WITHOUT ESD APPROXIMATION.  ALIAS IS THE ELECTRON
C     EMBEDDING, WHICH IS CALCULATED IN DFTBESP SUBROUTINE.
C     \SUM^K_{K \NE I,J} \SUM^{M_K}_{C \IN K} \GAMMA_{AC} \DELTA q_C^K
C
C     X(LSHIFTG+NATFMO+NBDFG ~ LSHIFTG+2(NATFMO+NBDFG)) CONTAINS THE
C     CONTRIBUTIONS WHICH COME FROM ESD APPROXIMATION.  WE CAN SIMPLY
C     ADD THE CONTRIBUTION TO THE X(LSHIFT) ARRAY, AND CAN BE
C     CALCULATED IN ONE TIME.
C     \SUM^{M_J}_{B \NE J} \GAMMA_{AB} \DELTA q_B^J
C
C     ABOVE TWO ARRAYS ARE USED FOR D*(DS/DR) CONTRIBUTION
C
C     WE NEED ANOTHER TERM, THE LAGRANGIAN-LIKE TERM (WS).
C     WE NEED ONE MORE SHIFT CONTRIBUTION, X(LSHIFTG), TO CALCULATE
C     THE WS TERM.  THIS COMES FROM THE DIMER CALCULATIONS, BUT IT
C     HAS NOT CALCULATED YET WHEN \DELTA WS IS CALCULATED.  THAT IS,
C       \DELTA WS(IJ) = WS(IJ) - WS(I) (+) WS(J)
C     THE SECOND AND THIRD TERMS ARE CALCULATED HERE, BUT THE FIRST
C     TERM HAS BEEN CALCULATED IN DFTB_LAG SUBROUTINE.
C     THE CONTRIBUTION HAS TO BE SUBTRACTED, BECAUSE THE SIGN IN ABOVE
C     EQUATION IS MINUS.
C
        CALL VCLR(X(LSHIFT),1,MAXNAT)
        CALL DCOPY(NAT,X(LSHIFTCT+INDI-1),1,X(LSHIFT),1)
        CALL DAXPY(NAT,ONE,X(LSHIFTG+NATFMOB+INDI-1),1,X(LSHIFT),1)
        CALL DCOPY(NAT,X(LSHIFT),1,X(LSHIFT2),1)
        CALL DAXPY(NAT,-ONE,X(LSHIFTG+INDI-1),1,X(LSHIFT),1)
C
C     ----- CALCULATE OVERLAP DERIVATIVES -----
C
        CALL DFTB_EXTSDER(X(LPS),L2,1,NAT)
C
C     ----- CALCULATE OVERLAP INTEGRAL FOR MONOMER IFG -----
C
        CALL VCLR(X(LS),1,L2)
        CALL DFTB_HAM_OVER(X(LWRK2),X(LS),.FALSE.,.TRUE.,.FALSE.)
C
C     ----- STORE L_{mi}^{IJ,K} FOR SCZV -----
C
        IF (IAND(MODGRD,32).NE.0.AND.NDER.NE.2) THEN
C         CALL DCOPY(L3,X(LVEC),1,X(LWRK1),1)
C         CALL DFTB_ZVLAGGET(X(LSHIFTCT+INDI-1),X(LWRK1),X(LS),L1,
C    *         NAT,X(LIPTLG),X(LZVWRK),X(LZVLAG),NA,L1-NA,IFG)

          CALL VCLR(X(LWRK1),1,L2)
          CALL DFTB_SHIFT_FOCK(X(LSHIFTCT+INDI-1),X(LWRK1),X(LS),NAT,L2,
     *         NSPE,ISPE,IND,MAXANG)
          CALL CPYTSQ(X(LWRK1),X(LWRK2),L1,1)
          CALL DGEMM('T','N',L1-NA,L1,L1,ONE,X(LVEC+L1*NA),L1,
     *               X(LWRK2),L1,ZERO,X(LWRK1),L1-NA)
          CALL DGEMM('N','N',L1-NA,NA,L1,ONE,X(LWRK1),L1-NA,
     *               X(LVEC),L1,ZERO,X(LZVWRK),L1-NA)
          CALL ZVLGMULT(NA*(L1-NA),1,IFG,X(LIPTLG),
     *         IAND(MODPAR,512+1024).NE.0,X(LZVLAG),X(LZVWRK),X(LWRK1))
        END IF
C
        IF(NDER.EQ.2.AND.NA.GT.0) THEN
          CALL FMO2DER_2NDMULL(NAT,X(LVEC),X(LPS),X(LD),L1,L2,
     *         X(LYALAG),X(LIPTYA),NA,L1-NA,FCM,INDFRG,IATFRG,IFG,
     *         X(LSHIFTCT+INDI-1),X(LWRK2))
        END IF
C
C     ----- CALCULATE LAGRANGIAN-LIKE TERMS -----
C
      call dftb_addlag(0,l1,x(ld),x(ls),x(lds),x(lwrk1),x(lwrk2),
     *                 x(lshift),x(lshift2),nat,c,fmode(1,1,nbody),
     *                 iatfrg(indi))
      ! CALL CPYTSQ(X(LD),X(LDS),L1,1)
      ! CALL CPYTSQ(X(LS),X(LWRK2),L1,1)
      ! !! X(LWRK1) = X(LWRK2) * X(LDS) = LS * LD
      ! CALL MRARBR(X(LWRK2),L1,L1,L1,X(LDS),L1,L1,X(LWRK1),L1) !! (SD)
      ! CALL DS_SHIFT(L1,NAT,X(LDS),X(LSHIFT))
      ! !! X(LWRK2) = X(LDS) * X(LWRK1) = (LD*SHIFT) * (LS*LD)
      ! CALL MRARBR(X(LDS),L1,L1,L1,X(LWRK1),L1,L1,X(LWRK2),L1) !! (DSD+SHIFT)
C
C     ! ----- CALCULATE \DELTA D * (DS/DR) TERMS FOR MONOMER -----
C
      ! CALL CPYTSQ(X(LD),X(LDS),L1,1)
      ! CALL DS_SHIFT(L1,NAT,X(LDS),X(LSHIFT2))
      ! CALL DAXPY(L3,-0.5D+00,X(LWRK2),1,X(LDS),1)
      ! IF (MASWRK) CALL ADD_LAG_GRAD(L1,L2,NAT,NATFMO,IATFRG(INDI),
     *!                   FMODE(1,1,NBODY),X(LPS),X(LDS),X(LWRK2))
C
        IF (RFLMO(1).NE.0) THEN
          L0=L1
          ICURFG = IFG
          IODA(12)=0
          IODA(15)=0
          IODA(318)=0
          IODA(569)=0
C
          CALL DCOPY(NAT,X(LZREF+INDI-1),1,ZAN,1)
          CALL ICOPY(NAT,X(LIANFRG+INDI-1),1,IAN,1)
          NSHELL = 0
          KLOC(1) = 1
          DO I = 1, NAT
            DO J = 1, MAXANG(ISPE(I))
              NSHELL = NSHELL + 1
              KATOM(NSHELL) = I
              KTYPE(NSHELL) = J
              IF (J.EQ.1) THEN
                KMIN(NSHELL) = 1
                KMAX(NSHELL) = 1
              ELSE IF (J.EQ.2) THEN
                KMIN(NSHELL) = 2
                KMAX(NSHELL) = 4
              ELSE IF (J.EQ.3) THEN
                KMIN(NSHELL) = 5
                KMAX(NSHELL) = 9
              END IF
              IF (NSHELL.NE.1) KLOC(NSHELL) = KLOC(NSHELL-1)
     *          + KMAX(NSHELL-1) - KMIN(NSHELL-1) + 1
            END DO
          END DO
C
          CALL DAWRIT(IDAF,IODA,X(LS),L2,12,0)
          CALL QMTSYM(X(LS),X(LWRK1),X(LWRK2),X(LE),X(LSCR),
     *                X(LIWRK),L0,L1,L2,L3,.FALSE.)
          CALL FMOORBS(X(LVECL),X(LS),X(LE),X(LWRK2),X(LSCR),
     *                 X(LIWRK),L0,L1,L2,L3,SKIPOR,1)
C
          IF (NUMFRZ.EQ.0) CYCLE
C
          NFRZ = MAX(NUMFRZ,NORFRZ)
          CALL VALFM(LOADFM)
          LUROT   = LOADFM  + 1
          LRROT   = LUROT   + NFRZ*NFRZ
          LAST    = LRROT   + NFRZ*NFRZ
          NEED2   = LAST - LOADFM - 1
          CALL GETFM(NEED2)
C
          CALL DAREAD(IDAF,IODA,X(LUROT),NFRZ*NFRZ,569,0)
          CALL DFTB_RROT(L1,NFRZ,NUMFRZ,X(LUROT),X(LRROT),X(LS),
     *                   X(LVECL),X(LVEC),X(LWRK1),X(LWRK2))
          CALL VCLR(X(LFAO),1,L2)
          CALL DFTB_SHIFT_FOCK(X(LSHIFTCT+INDI-1),X(LFAO),X(LS),NAT,L2,
     *         NSPE,ISPE,IND,MAXANG)
          CALL VCLR(X(LGMK),1,L1*NFRZ)
          CALL DFTB_CNSTGMK(L1,NA,NUMFRZ,NFRZ,ONE,
     *                     X(LGMK),X(LS),X(LFAO),X(LVEC),X(LRROT),
     *                     X(LWRK1),X(LWRK2),X(LWRK3))
C
          JCURFG = -1
          CALL DFTB_PUTGMK(L1,1,X(LGMK),X(LGMKSAV),X(LINDFRZ),X(LIALMO),
     *                     X(LIAGLOB),X(LAST),X(LAST),X(LAST),.FALSE.,
     *                     X(LAST),X(LAST),X(LAST))
          JCURFG = 0
C
          CALL RETFM(NEED2)
        END IF
      END DO
C
      IF (ISGDDI) CALL GDDICOUNT( 1,LGROUP,MYJOB)
C
      CALL RETFM(NEED)
C
C     call valfm(loadfm)
C     lfmodewrk = loadfm    + 1
C     last      = lfmodewrk + natfmo*3
C     need1     = last - loadfm - 1
C     call getfm(need1)
C     call dcopy(natfmo*3,fmode(1,1,nbody),1,x(lfmodewrk),1)
C     if (isgddi) then
C       call gddi_scope(ddi_world)
C       call ddi_gsumf(2424,x(lfmodewrk),natfmo*3)
C       call gddi_scope(ddi_group)
C       call gddicount(-1,lgroup,myjob)
C     end if
C     if (maswrk) then
C       write (6,*) "fmode after fmograd2"
C       do i = 1, natfmo
C         write (6,'(i3,3f20.10)') i,(x(lfmodewrk+k-1+3*(i-1)),k=1,3)
C       end do
C     end if
C     call retfm(need1)
C
      RETURN
C
      END SUBROUTINE DFTB_FMOGRAD2
C
C-----------------------------------------------------------------------
C*MODULE DFTBFO    *DECK DFTB_ESPGRAD
C>
C>    @brief Gradient in FMO-DFTB
C>
C>    @details Calculate gradient of ESP contribution in FMO-DFTB.
C>
C>    @author Yoshio Nishimoto
C>    - Sep, 2014- Subroutine written
C>    @date   Nov, 2014- Yoshio Nishimoto
C>    - Changes for FMO-DFTB3 gradient
C>    @date   Feb, 2016- Yoshio Nishimoto
C>    - Performance improvement
C>    @date   Feb, 2018- Yoshio Nishimoto
C>    - Added FMO-DFTB/AFO etc.
C>
C>           --- INPUT ---
C>    @param DM      Difference density matrix (\Delta D_{\mu \nu}^{IJ})
C>    @param DM2     Current density matrix (D_{\mu \nu}^{IJ})
C>    @param CTIJ    Charge transfer (CTIJ(A) = q_A^{IJ} - q_A^I - q_A^J)
C>    @param POPMAT  Mulliken population of monomer
C>    @param IAGLOB  Global atom index
C>    @param INDFRG  First atom of each fragment in terms of NATFMOB?
C>    @param ZREFFMO ZREF for all atoms (NATFMOB length)
C>    @param GAMMA   DFTB2 gamma values
C>    @param FMOC    Coordinates for all NATFMO
C>    @param ISPEFMO ISPE for all atoms (NATFMOB length)
C>    @param INDBD   Something about BDA or BAA
C>    @param IALOC   Converts global to local (fragment) atomic location
C>    @param INDAT   Converts global atomic location to fragment index
C>    @param CTIJ2   Square charge transfer used in FMO-DFTB3
C>    @param CHAMUL  Mulliken charges of dimer IJ (for FMO-DFTB3)
C>    @param GAMMA2  DFTB2 Gamma for internal fragment?
C>    @param GAMMA3  DFTB3 Gamma for internal fragment?
C>    @param FMOAFO  FMO/AFO or not
C>           --- IN/OUTPUT ---
C>    @param FMODE   Cartesian geometrical gradient of monomers and dimers
C>    @param SHIFTCT Shift contribution of charge transfer
C>    @param SHIFTG  Global shift contribution
C>
C
      SUBROUTINE DFTB_ESPGRAD(DM,DM2,FMODE,CTIJ,POPMAT,iaglob,INDFRG,
     * SHIFTCT,ZREFFMO,GAMMA,GAMMAD,SHIFTG,
     * FMOC,ISPEFMO,INDBD,IALOC,INDAT,CTIJ2,CHAMUL,GAMMA2,GAMMA3,FMOAFO)
      use mx_limits, only: mxatm
      USE DFTBPB_MOD,ONLY: PERIOD,DFTB_GAMDER_PBC,STRFMO
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      INTEGER,PARAMETER :: MXSPE=10
      DOUBLE PRECISION, PARAMETER :: ZERO=0.0D+00,ONE=1.0D+00,
     * TWO=2.0D+00,ONE_THIRD=1.0D+00/3.0D+00
C
      LOGICAL DFTBFL,SCC,SRSCC,DFTB3,DAMPXH0,LCDFTB
      LOGICAL GOPARR,DSKWRK,MASWRK
      LOGICAL FMOAFO
C
      COMMON /FMCOM / X(1)
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      common /fmorun/ espscf,e0scf(2),emp2s,IDAFMO,icurfg,jcurfg,kcurfg,
     *                icurlay,icurunt,nat1e,ncursh,ngau,icurpop,ifmostp,
     *                moncor,needr,modrst,norbproj,nunesp,iskipesp,
     *                IESDPPC,idoprop,mp2run,icurit,idmfmo,iddfmo,
     *                iddcur,nddleft,ivmfmo,nzmtfmo,ifmobas,itmfmo(2)
      common /fmodim/ maxbnd,maxknd,maxcbs,maxcao,maxvec,maxl1,maxnat,
     *                maxabd,maxbas,maxbbd,maxlmo,maxslo,maxabd2,maxrij
      COMMON /DFTB  / DFTBFL,SCC,SRSCC,DFTB3,DAMPXH0,LCDFTB
      COMMON /DFTBPR/ ETEMP,DFTBDP(MXSPE*14),DAMPXHE,HUBDER(MXSPE),
     *                ZREF(MXATM),SPNCST(6,MXSPE),SPE(MXATM),NSPE,
     *                MAXANG(MXATM),ISPE(MXATM),IND(MXATM+1),IDFTBD,
     *                PARAMDIR
      COMMON /DFTBSK/ SKDIM(MXSPE,MXSPE),SKSPIN(MXSPE),QREFL(3,MXSPE),
     *                HUBBL(3,MXSPE),QREF(MXSPE),HUBB(MXSPE),SKCUT2,
     *                NEEDSK,LSKHTAB(MXSPE,MXSPE),LSKSTAB(MXSPE,MXSPE),
     *                LSKGRID(MXSPE,MXSPE),LSKSELF(MXSPE)
      COMMON /FMOINF/ NFG,NLAYER,NATFMO,NBDFG,NAOTYP,NBODY,NSEGM
      COMMON /FMODTB/ LGAMMA,LZREF,LISPE,LSHIFTG,LSHIFTCT,
     *                LGAMMAD,MODESD,MODGAMMA,LDCTMP,LGAMMA3,LGAMMA3D,
     *                LNBOND,LCTMUL,LESPDFTB,LESPASC,NDFTBRST,LDFTBRST
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
C
      double precision popmat(maxnat,nfg),fmode(3,NATfmo,*),CTIJ(MXATM),
     * CTIJ2(*),SHIFTCT(*),GAMMA(*),GAMMAD(*),SHIFT(NAT),SHIFTG(*),
     * SHIFT3(NAT*2),FMOC(3,NATFMO),ZREFFMO(*),CHAMUL(*),
     *                 DM(*), !! DIFFERENCE DENSITY MATRIX
     *                 DM2(*),!! CURRENT DENSITY MATRIX
     * GAMMA2(*),GAMMA3(NAT,NAT)
      integer iaglob(*),indfrg(*),
     * INDBD(MAXABD,*),IALOC(*),INDAT(*),ISPEFMO(*)
C
      DOUBLE PRECISION :: VEC(3),SHIFTMP(NAT),TMPD(10,NSPE,NSPE),
     * TMP(12,NSPE,NSPE),TMP3(10,NSPE,NSPE),TMPD3(11,NSPE,NSPE)
      LOGICAL DAMPXH
      !! Gamma derivative vector and stress tensor for PBC
      DIMENSION GDV2(3),GDV3I(3),GDV3J(3),
     *          GDST2(3,3),GDST3I(3,3),GDST3J(3,3)
C
C     ------ ADD GRADIENT for FMO-DFTB corresponding to the
C     SECOND-ORDER HAMILTONIAN term in DFTB -----
C
      IF (.NOT.SCC) RETURN
      IF (IAND(MODGAMMA,4096).NE.0) RETURN
C
C     TIM0 = ZERO
C     TIM1 = ZERO
C     CALL TSECND(TIM0)
C     TLEFTS = TIMLIM - TIM0
C
      ifg=icurfg
      jfg=jcurfg
      kfg=kcurfg
      SCAL = ONE
      if (nbody.eq.3.and.kfg.eq.0) then
        iifg=max(ifg,jfg)
        jjfg=min(ifg,jfg)
        ijfg=(iifg*iifg-3*iifg)/2+jjfg+1
        ida=1
        ida=ida-ixftch(x(lndtfrg),ijfg)
        SCAL = DBLE(IDA)
      end if
C
      IF (.NOT.MASWRK) GO TO 100
C
C     ------ ACCUMULATE CHARGE TRANSFER CONTRIBUTION -----
C     HOWEVER, SUBTRACT THE INTRA-FRAGMENT GRADIENT CONTRIBUTION
C     FOR I-J TO AVOID DOUBLE COUNTING IN THE SUMMATION PART.
C
      CALL DCOPY(12*NSPE*NSPE,GAMMA,1,TMP,1)
      IF (DFTB3) CALL DCOPY(10*NSPE*NSPE,X(LGAMMA3),1,TMP3,1)
      CALL DCOPY(10*NSPE*NSPE,GAMMAD,1,TMPD,1)
      IF (DFTB3) CALL DCOPY(11*NSPE*NSPE,X(LGAMMA3D),1,TMPD3,1)
C
      VALTMP = ZERO
      TMPSH2 = ZERO
      TMPSH3 = ZERO
      DO I = 1, NAT
        IAG = IAGLOB(I)
        ISP = ISPEFMO(IAG)
        IATFG = INDAT(IAG)
        INDI = INDFRG(IATFG)
        CHAMUL1 = CTIJ(I)
        IATLOC = IAND(IALOC(IAG),65535)
        ILOC = INDI + IATLOC
        X(LCTMUL+ILOC-2) = X(LCTMUL+ILOC-2) + CTIJ(I)*SCAL
        IF (DFTB3) X(LCTMUL+NATFMO+NBDFG+ILOC-2)
     *               = X(LCTMUL+NATFMO+NBDFG+ILOC-2) + CTIJ2(I)*SCAL
        DO J = 1, NAT
          JAG = IAGLOB(J)
          JATFG = INDAT(JAG)
          INDJ = INDFRG(JATFG)
          JSP = ISPEFMO(JAG)
          JATLOC = IAND(IALOC(JAG),65535)
          JLOC = INDJ + JATLOC
          CHAMUL2 = POPMAT(JATLOC,JATFG) - ZREFFMO(JLOC-1)
          VEC(1) = FMOC(1,IAG) - FMOC(1,JAG)
          VEC(2) = FMOC(2,IAG) - FMOC(2,JAG)
          VEC(3) = FMOC(3,IAG) - FMOC(3,JAG)
          DIST = SQRT(VEC(1)*VEC(1)+VEC(2)*VEC(2)+VEC(3)*VEC(3))
          DISTI = 0
          IF(DIST.NE.0) DISTI = ONE/DIST
          CALL DFTB_CNVSQ(I,J,NSEQ)
          GAMTMP = GAMMA2(NSEQ)
          TMPSH = GAMTMP*CHAMUL1
          IF (DFTB3) THEN
            GAMI = GAMMA3(I,J)
            GAMJ = GAMMA3(J,I)
            TMPSH2 = CTIJ2(I)*GAMI
            TMPSH3 = CHAMUL2*CHAMUL1*GAMJ
            GAMTMPJ = GAMJ
          END IF
          VALTMP3 = TMPSH+ONE_THIRD*(TMPSH2+TMPSH3*TWO)
          SHIFTCT(JLOC-1) = SHIFTCT(JLOC-1) - VALTMP3*SCAL
          IF (IAG.NE.JAG.OR.PERIOD) THEN
            DAMPXH = DAMPXH0
            CALL DFTB_DMPCHK(ISP,JSP,DAMPXH)
            IF (PERIOD) THEN
              CALL DFTB_GAMDER_PBC(I,J,FMOC(1,IAG),FMOC(1,JAG),GDV2,
     *                             GDV3I,GDV3J,GDST2,GDST3I,GDST3J,
     *                             TMPD(1,ISP,JSP),TMPD3(1,ISP,JSP),
     *                             TMPD3(1,JSP,ISP),DAMPXH,ISP,JSP)
              !! Stress tensor
              VAL = CHAMUL1*CHAMUL2
              CALL DAXPY(9,-SCAL*VAL,GDST2,1,STRFMO,1)
              IF (DFTB3) THEN
                VAL = ONE_THIRD*CTIJ2(I)*CHAMUL2
                CALL DAXPY(9,-SCAL*VAL,GDST3I,1,STRFMO,1)
                VAL = ONE_THIRD*CHAMUL1*CHAMUL2*CHAMUL2
                CALL DAXPY(9,-SCAL*VAL,GDST3J,1,STRFMO,1)
              END IF
            ELSE
              CALL DFTB_GAM23D(DIST,DISTI,GAMMAV,GAMI,GAMJ,
     *                         TMPD(1,ISP,JSP),TMPD3(1,ISP,JSP),
     *                         TMPD3(1,JSP,ISP),DAMPXH,DAMPXHE,SPE(ISP),
     *                         SPE(JSP),HUBDER(ISP),HUBDER(JSP),DFTB3)
              DO K = 1, 3
                GDV2(K) = GAMMAV*VEC(K)*DISTI
                IF (DFTB3) THEN
                  GDV3I(K) = GAMI*VEC(K)*DISTI
                  GDV3J(K) = GAMJ*VEC(K)*DISTI
                END IF
              END DO
            END IF
C           VAL = CHAMUL1*CHAMUL2*GAMMAV*DISTI
C           IF (DFTB3) THEN
C             VAL = VAL + ONE_THIRD*(CTIJ2(I)*GAMI
C    *          +CHAMUL1*CHAMUL2*GAMJ)*CHAMUL2*DISTI
C           END IF
C           VALTMP2 = VAL
C           VAL = VAL*SCAL
C           FMODE(1,IAG,NBODY) = FMODE(1,IAG,NBODY) + VEC(1)*VAL
C           FMODE(2,IAG,NBODY) = FMODE(2,IAG,NBODY) + VEC(2)*VAL
C           FMODE(3,IAG,NBODY) = FMODE(3,IAG,NBODY) + VEC(3)*VAL
C           FMODE(1,JAG,NBODY) = FMODE(1,JAG,NBODY) - VEC(1)*VAL
C           FMODE(2,JAG,NBODY) = FMODE(2,JAG,NBODY) - VEC(2)*VAL
C           FMODE(3,JAG,NBODY) = FMODE(3,JAG,NBODY) - VEC(3)*VAL
            DO K = 1, 3
              VAL = CHAMUL1*CHAMUL2*GDV2(K)
              IF (DFTB3) THEN
                VAL = VAL
     *            + ONE_THIRD*(CTIJ2(I)*GDV3I(K)
     *                        +CHAMUL1*CHAMUL2*GDV3J(K))*CHAMUL2
              END IF
              IF (NBODY.EQ.3.AND.KFG.EQ.0) VAL = VAL*SCAL
              FMODE(K,IAG,NBODY) = FMODE(K,IAG,NBODY) - VAL
              FMODE(K,JAG,NBODY) = FMODE(K,JAG,NBODY) + VAL
            END DO
          END IF
          IF (NBDFG.NE.0) THEN
            DO IABD = 1, MAXABD
              IBDFG = INDBD(IABD,JAG)
              IF (IBDFG.EQ.0) CYCLE
              KATFG = ISHFT(IALOC(NATFMO+IBDFG),-16)
              KATLOC = IAND(IALOC(NATFMO+IBDFG),65535)
              INDK = INDFRG(KATFG)
              KLOC = INDK + KATLOC
              IF ((KATFG.EQ.IFG.OR.KATFG.EQ.JFG.OR.KATFG.EQ.KFG)
     *            .AND.JATFG.NE.KATFG)THEN
                CHAMUL3 = POPMAT(KATLOC,KATFG) - ZREFFMO(KLOC-1)
                IF (DFTB3) TMPSH3 = CHAMUL3*CHAMUL1*GAMTMPJ
                VAL = TMPSH+ONE_THIRD*(TMPSH2+TMPSH3*TWO)
                SHIFTCT(KLOC-1) = SHIFTCT(KLOC-1) - VAL*SCAL
                IF (IAG.NE.JAG.OR.PERIOD) THEN
                  DO K = 1, 3
                    VAL = CHAMUL1*CHAMUL3*GDV2(K)
                    IF (DFTB3) THEN
                      VAL = VAL
     *                  + ONE_THIRD*(CTIJ2(I)*GDV3I(K)
     *                              +CHAMUL1*CHAMUL3*GDV3J(K))*CHAMUL3
                    END IF
                    IF (NBODY.EQ.3.AND.KFG.EQ.0) VAL = VAL*SCAL
                    FMODE(K,IAG,NBODY) = FMODE(K,IAG,NBODY) - VAL
                    FMODE(K,JAG,NBODY) = FMODE(K,JAG,NBODY) + VAL
                  END DO
                  IF (PERIOD) THEN
                    !! Stress tensor
                    VAL = CHAMUL1*CHAMUL3
                    CALL DAXPY(9,-SCAL*VAL,GDST2,1,STRFMO,1)
                    IF (DFTB3) THEN
                      VAL = ONE_THIRD*CTIJ2(I)*CHAMUL3
                      CALL DAXPY(9,-SCAL*VAL,GDST3I,1,STRFMO,1)
                      VAL = ONE_THIRD*CHAMUL1*CHAMUL3*CHAMUL3
                      CALL DAXPY(9,-SCAL*VAL,GDST3J,1,STRFMO,1)
                    END IF
                  END IF
                END IF
                !! AVOID DOUBLE COUNT
                IF (JATFG.NE.IFG.AND.JATFG.NE.JFG.AND.JATFG.NE.KFG) THEN
                  SHIFTCT(JLOC-1) = SHIFTCT(JLOC-1) + VALTMP3*SCAL
                  IF (IAG.NE.JAG.or.period) THEN
                    DO K = 1, 3
                      VAL = CHAMUL1*CHAMUL2*GDV2(K)
                      IF (DFTB3) THEN
                        VAL = VAL
     *                    + ONE_THIRD*(CTIJ2(I)*GDV3I(K)
     *                                +CHAMUL1*CHAMUL2*GDV3J(K))*CHAMUL2
                      END IF
                      IF (NBODY.EQ.3.AND.KFG.EQ.0) VAL = VAL*SCAL
                      FMODE(K,IAG,NBODY) = FMODE(K,IAG,NBODY) + VAL
                      FMODE(K,JAG,NBODY) = FMODE(K,JAG,NBODY) - VAL
                    END DO
                    IF (PERIOD) THEN
                      !! Stress tensor
                      VAL = CHAMUL1*CHAMUL2
                      CALL DAXPY(9,SCAL*VAL,GDST2,1,STRFMO,1)
                      IF (DFTB3) THEN
                        VAL = ONE_THIRD*CTIJ2(I)*CHAMUL2
                        CALL DAXPY(9,SCAL*VAL,GDST3I,1,STRFMO,1)
                        VAL = ONE_THIRD*CHAMUL1*CHAMUL2*CHAMUL2
                        CALL DAXPY(9,SCAL*VAL,GDST3J,1,STRFMO,1)
                      END IF
                    END IF
                  END IF
                END IF
              END IF
            END DO
          END IF
        END DO
      END DO
C     
  100 CONTINUE
      RETURN
c     CALL VCLR(SHIFT,1,NAT)
c     CALL DAREAD(IDAF,IODA,SHIFT,NAT,566,0)
c     IF (NBODY.EQ.3) CALL DSCAL(NAT,SCAL,SHIFT,1)
c     IF (DFTB3) THEN
c       CALL DAREAD(IDAF,IODA,SHIFT3,NAT*2,567,0)
c       IF (NBODY.EQ.3)  CALL DSCAL(NAT*2,SCAL,SHIFT3,1)
c       CALL VCLR(SHIFTMP,1,NAT)
c       DO I = 1, NAT
c         SHIFTMP(I) = TWO*CHAMUL(I)*SHIFT3(I)
c       END DO
c       CALL DAXPY(NAT,1.0D+00,SHIFTMP,1,SHIFT,1)
c       DO I = 1, NAT
c         SHIFT(I) = SHIFT(I) + SHIFT3(NAT+I)
c       END DO
c     END IF
c     !! ON ENTRY, SHIFT(I) = SHIFT(I)+2*SHIFT3(I)*Q(I)+SHIFT3(NAT+I)
c     CALL DFTB_LAG(IAGLOB,FMODE(1,1,NBODY),SHIFT,SHIFTMP,DM,DM2,
c    *              MASWRK,DFTB3,FMOAFO)
c     !! ON EXIT,  SHIFT(I) = SHIFT(I)                 +SHIFT3(NAT+I)
c     IF (.NOT.MASWRK) RETURN
C
C     ----- ACCUMULATE SHIFT CONTRIBUTIONS -----
C     THE SIGN OF DFTB3 CONTRIBUTION IS OPPOSITE.  THE LATTER HALF OF
C     SHIFTG IS ORIGINALLY USED FOR ESD GRADIENT WHOSE SIGN IS PLUS.
C     NOW, DFTB3 CONTRIBUTION IS ONLY APPLIED FOR - (W(I) (+) W(J))
C     AND - (DS(I) (+) DS(J)) WHOSE SIGNS ARE MINUS.  THEREFORE,
C     THE SIGN OF DFTB3 CONTRIBUTION BELOW SHOULD BE MINUS.
C
c     NATFMOB = NATFMO + NBDFG
c     DO I = 1, NAT
c       IAG = IAGLOB(I)
c       IATFG = INDAT(IAG)
c       IATLOC = IAND(IALOC(IAG),65535)
c       INDI = INDFRG(IATFG)
c       ILOC = INDI + IATLOC
c       SHIFTG(ILOC-1) = SHIFTG(ILOC-1) + SHIFT(I)
c       IF (DFTB3) SHIFTG(NATFMOB+ILOC-1) = SHIFTG(NATFMOB+ILOC-1)
c    *    - TWO*SHIFT3(I)*(POPMAT(IATLOC,IATFG)-ZREFFMO(ILOC-1))
c       IF (NBDFG.NE.0) THEN
c         DO IABD = 1, MAXABD
c           IBDFG = INDBD(IABD,IAG)
c           IF (IBDFG.EQ.0) CYCLE
c           JATFG = ISHFT(IALOC(NATFMO+IBDFG),-16)
c           INDJ =  INDFRG(JATFG)
c           JATLOC = IAND(IALOC(NATFMO+IBDFG),65535)
c           IF ((JATFG.EQ.IFG.OR.JATFG.EQ.JFG.OR.JATFG.EQ.KFG)
c    *          .AND.JATFG.NE.IATFG) THEN
c             JLOC = INDJ + JATLOC
c             SHIFTG(JLOC-1) = SHIFTG(JLOC-1) + SHIFT(I)
c             IF (DFTB3) SHIFTG(NATFMOB+JLOC-1) = SHIFTG(NATFMOB+JLOC-1)
c    *          - TWO*SHIFT3(I)*(POPMAT(JATLOC,JATFG)-ZREFFMO(JLOC-1))
c             !! AVOID DOUBLE COUNT
c             IF (IATFG.NE.IFG.AND.IATFG.NE.JFG.AND.IATFG.NE.KFG) THEN
c               SHIFTG(ILOC-1) = SHIFTG(ILOC-1) - SHIFT(I)
c               IF (DFTB3) THEN
c                 SHIFTG(NATFMOB+ILOC-1) = SHIFTG(NATFMOB+ILOC-1)
c    *            + TWO*SHIFT3(I)*(POPMAT(IATLOC,IATFG)-ZREFFMO(ILOC-1))
c               END IF
c             END IF
c           END IF
c         END DO
c       END IF
c     END DO
C
C     IF(MASWRK)
C    *   WRITE(IW,FMT='('' ...... END OF ESP GRADIENT ......'')')
C     CALL TIMIT(1)
C
c     RETURN

      END SUBROUTINE DFTB_ESPGRAD
C
C-----------------------------------------------------------------------
C*MODULE DFTBFO    *DECK DFTB_LAG
C>
C>    @brief Gradient in FMO-DFTB
C>
C>    @details Calculate gradient of contributions of Lagrangian-like
C>             terms in FMO-DFTB.
C>
C>    @author Yoshio Nishimoto
C>    - Sep, 2014- Subroutine written
C>    @date   Nov, 2014- Yoshio Nishimoto
C>    - Changes for FMO-DFTB3 gradient
C>    @date   Feb, 2018- Yoshio Nishimoto
C>    - Added FMO-DFTB/AFO etc.
C>
C>           --- INPUT ---
C>    @param IAGLOB  Global atom index
C>    @param SHIFT   Shift contribution of current SCF calculation
C>    @param SHIFTMP A part of DFTB3 contribution, 2/3*\Delta q_A^IJ*
C>                   \Delta q_C^K* \Gamma_{AC}, see DFTB_ESPGRAD.
C>    @param DM      Difference density matrix (\Delta D_{\mu \nu}^{IJ})
C>    @param DM2     Current density matrix
C>    @param MASWRK  Master or not
C>    @param DFTB3   DFTB3 or not
C>    @param FMOAFO  FMO/AFO or not
C>           --- IN/OUTPUT ---
C>    @param FMODE   Cartesian geometrical gradient of monomers and dimers
C>
C
      SUBROUTINE DFTB_LAG(IAGLOB,FMODE,SHIFT,SHIFTMP,DM,DM2,
     *  MASWRK,DFTB3,FMOAFO)
      use mx_limits, only: mxatm
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      COMMON /FMCOM / X(1)
      COMMON /FMOINF/ NFG,NLAYER,NATFMO,NBDFG,NAOTYP,NBODY,NSEGM
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
C
      INTEGER :: IAGLOB(*)
      DOUBLE PRECISION :: FMODE(3,NATFMO,*),SHIFT(*),SHIFTMP(*),
     *                    DM(*), !! DIFFERENCE DENSITY MATRIX
     *                    DM2(*) !! CURRENT DENSITY MATRIX
      LOGICAL :: MASWRK,DFTB3,FMOAFO
C
C     THIS SUBROUTINE CALCULATES PARTS OF FMO-DFTB DIMER GRADIENTS.
C     IT CALCULATES ONLY TWO TERMS:
C     W_{\MU\NU}^{IJ,\ALPHA} S_{\MU\NU}^{IJ}
C        AND
C     \DELTA D_{\MU\NU}^{IJ} \FRAC{\PARTIAL S_{\MU\NU}^{IJ}}{\PARTIAL R}
C     NOTE THAT \DELTA WS = WS(IJ) - W(I) (+) W(J)
C
C     OTHER TERMS SUCH AS THE GRADIENTS OF ENVIRONMENTAL MONOMERS AND
C     MONOMER -WS TERMS, NEEDED FOR \DELTA WS, WILL BE CALCULATED IN
C     THE DFTB_FMOGRAD SUBROUTINE.
C
      L1  = NUM
      L2  = (L1*L1+L1)/2
      L3  = L1 * L1
C
C     ----- ALLOCATE MEMORY -----
C
      CALL VALFM(LOADFM)
      CALL GOTFM(NGOTMX)
      LWRK1   = LOADFM  + 1
      LWRK2   = LWRK1   + L3
      LS      = LWRK2   + L3
      LDS     = LS      + L2
      LPS     = LDS     + L3 !! D SQUARE
      LAST    = LPS     + L2*NAT*3
      IF (FMOAFO.AND.DFTB3) THEN
        LWRK3 = LAST
        LAST  = LWRK3   + L3
      END IF
      NEED    = LAST - LOADFM - 1
C
C     ----- GET MEMORY -----
C
      CALL GETFM(NEED)
C
C     ----- CALCULATE OVERLAP DERIVATIVES -----
C
      CALL DFTB_EXTSDER(X(LPS),L2,1,NAT)
C
C     ----- RESTOR OVERLAP INTEGRAL FOR FRAGMENT I OR J FROM DICT. -----
C
      CALL VCLR(X(LS),1,L2)
      CALL DAREAD(IDAF,IODA,X(LS),L2,12,0)
C
C     ----- CALCULATE LAGRANGIAN TERM GRADIENTS -----
C       \DELTA W = W(IJ) - W(I) (+) W(J)
C     ONLY THE FIRST W(IJ) TERM IS CALCULATED HERE.  THE LATTER TERM IS
C     CALCULATED IN DFTB_FMOGRAD2.  DFTB_ESPGRAD ACCUMULATES
C     CONTRIBUTION.
C
      CALL CPYTSQ(DM2,X(LDS),L1,1)
      CALL CPYTSQ(X(LS),X(LWRK2),L1,1)
      CALL MRARBR(X(LWRK2),L1,L1,L1,X(LDS),L1,L1,X(LWRK1),L1) !! (SD)
      CALL DS_SHIFT(L1,NAT,X(LDS),SHIFT)
      CALL MRARBR(X(LDS),L1,L1,L1,X(LWRK1),L1,L1,X(LWRK2),L1) !! (DSD+SHIFT)
C
C     ----- CALCULATE \DELTA D * (DS/DR) TERMS -----
C     ONLY FOR I-J DIMER.  
C     \DELTA D (DM) IS COPIED TO X(LDS), AND THEN INSIDE DS_SHIFT,
C     SHIFT MATRIX IS MULTIPLIED.
C
C     FOR DFTB3, BELOW CALCULATES THE FIRST TERM OF
C       \DELTA DS = DS(IJ) - DS(I) (+) DS(J).
C     OTHER TERMS ARE CALCULATED LATER IN DFTB_FMOGRAD2.
C     DFTB_ESPGRAD ACCUMULATES THE TERM NEEDED TO CALCULATE MONOMER
C     CONTRIBUTIONS AS THE LATTER PART OF SHIFTG MATRIX.
C
      IF (DFTB3) THEN
        CALL DAXPY(NAT,-1.0D+00,SHIFTMP,1,SHIFT,1)
      END IF
      CALL CPYTSQ(DM,X(LDS),L1,1) !! USE DIFFERENCE DM
      CALL DS_SHIFT(L1,NAT,X(LDS),SHIFT)
      IF (DFTB3) THEN
        CALL CPYTSQ(DM2,X(LWRK1),L1,1) !! CURRENT DM
        IF (FMOAFO) THEN
          CALL DAREAD(IDAF,IODA,X(LWRK3),L3,476,0)
          CALL SYMTRZE(X(LWRK3),L1,L1)
          CALL DAXPY(L3,1.0D+00,X(LWRK3),1,X(LWRK1),1)
        END IF
        CALL DS_SHIFT(L1,NAT,X(LWRK1),SHIFTMP)
        CALL VADD(X(LDS),1,X(LWRK1),1,X(LDS),1,L3)
      END IF
      CALL DAXPY(L3,-0.5D+00,X(LWRK2),1,X(LDS),1)
      IF (MASWRK) CALL ADD_LAG_GRAD(L1,L2,NAT,NATFMO,IAGLOB,
     *  FMODE,X(LPS),X(LDS),X(LWRK1))
C
      CALL RETFM(NEED)
C
      RETURN
C
      END SUBROUTINE DFTB_LAG
C
C-----------------------------------------------------------------------
C*MODULE DFTBFO    *DECK DFTB_ESDGRAD
C>
C>    @brief Gradient in FMO-DFTB
C>
C>    @details Calculate gradient of ES-DIM contribution in FMO-DFTB.
C>
C>    @author Yoshio Nishimoto
C>    - Sep, 2014- Subroutine written
C>    @date   Feb, 2016- Yoshio Nishimoto
C>    - Some improvements
C>    @date   Jun, 2016- Yoshio Nishimoto
C>    - Added FMO3-DFTB
C>
C>           --- INPUT ---
C>    @param IFG,JFG I and J fragment index
C>    @param ISPEFMO ISPE for all atoms (NATFMO)
C>    @param INDFRG  First atom of each fragment in terms of NATFMOB?
C>    @param NATFRG  Number of atoms in each fragment
C>    @param IATFRG  Global atom index with NBDFG?
C>    @param GAMMAD  Derivatives of DFTB2 gamma values
C>    @param FMOC    Coordinates for all NATFMO
C>    @param GAMMAD3 Derivatives of DFTB3 gamma values
C>    @param CHAMULI Mulliken population of IFG
C>    @param CHAMULJ Mulliken population of IFG
C>    @param DISTMAT Distance matrix between IFG-JFG
C>    @param GAMMA3  GAMMA3...?
C>    @param NN
C>    @param NN2
C>           --- OUTPUT ---
C>    @param SHIFTG Global shift contribution
C>    @param FMODE Cartesian geometrical gradient of monomers and dimers
C>
C
      SUBROUTINE DFTB_ESDGRAD(IFG,JFG,ISPEFMO,
     * INDFRG,NATFRG,IATFRG,GAMMAD,SHIFTG,FMOC,FMODE,
     * GAMMAD3,CHAMULI,CHAMULJ,GAMMA3,DISTMAT,NN,NN2,SCAL)
      use mx_limits, only: mxatm
      USE DFTBPB_MOD,ONLY: PERIOD,DFTB_GAMDER_PBC,STRFMO
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER,PARAMETER :: MXSPE=10
      DOUBLE PRECISION, PARAMETER :: ZERO=0.0D+00,ONE=1.0D+00,
     *  TWO=2.0D+00
C
      COMMON /FMCOM / X(1)
      COMMON /DFTB  / DFTBFL,SCC,SRSCC,DFTB3,DAMPXH0,LCDFTB
      LOGICAL DFTBFL,SCC,SRSCC,DFTB3,DAMPXH0,LCDFTB,dotrunces,dotruncesd
      COMMON /DFTBPR/ ETEMP,DFTBDP(MXSPE*14),DAMPXHE,HUBDER(MXSPE),
     *                ZREF(MXATM),SPNCST(6,MXSPE),SPE(MXATM),NSPE,
     *                MAXANG(MXATM),ISPE(MXATM),IND(MXATM+1),IDFTBD,
     *                PARAMDIR
      COMMON /DFTBSK/ SKDIM(MXSPE,MXSPE),SKSPIN(MXSPE),QREFL(3,MXSPE),
     *                HUBBL(3,MXSPE),QREF(MXSPE),HUBB(MXSPE),SKCUT2,
     *                NEEDSK,LSKHTAB(MXSPE,MXSPE),LSKSTAB(MXSPE,MXSPE),
     *                LSKGRID(MXSPE,MXSPE),LSKSELF(MXSPE)
C
      COMMON /FMOINF/ NFG,NLAYER,NATFMO,NBDFG,NAOTYP,NBODY,NSEGM
      COMMON /FMODTB/ LGAMMA,LZREF,LISPE,LSHIFTG,LSHIFTCT,
     *                LGAMMAD,MODESD,MODGAMMA,LDCTMP,LGAMMA3,LGAMMA3D,
     *                LNBOND,LCTMUL,LESPDFTB,LESPASC,NDFTBRST,LDFTBRST
      common /cenrad/ rcutes,lcradfg,lichfgg,dotrunces,dotruncesd
C
      DOUBLE PRECISION :: VEC(3), DIST
      DOUBLE PRECISION :: CHAMULI(*),CHAMULJ(*),GAMMA3(NN,NN2,*)
      DOUBLE PRECISION :: TMPD(10,NSPE,NSPE),
     * GAMMAD(*),SHIFTG(*),FMOC(3,*),
     * FMODE(3,NATFMO,*),TMP(12,NSPE,NSPE),TMP3(10,NSPE,NSPE),
     * TMPD3(11,NSPE,NSPE),GAMMAD3(NN,NN2,*),DISTMAT(*)
      INTEGER :: ISPEFMO(*),NATFRG(*),IATFRG(*),INDFRG(*)
      LOGICAL DAMPXH
      !! Gamma derivative vector and stress tensor for PBC
      DIMENSION GDV2(3),GDV3I(3),GDV3J(3),
     *          GDST2(3,3),GDST3I(3,3),GDST3J(3,3)
C
C     CALCULATE THE GRADIENT OF ELECTROSTATIC DIMER APPROXIMATION.
C     THIS SUBROUTINE CALCULATES PART OF GRADIENTS.  THE FIRST BLOCK
C     CALCULATES THE INTER-FRAGMENT GRADIENT, WHICH IS CALCULATED WITH
C     DERIVATIVES OF GAMMA.  ANOTHER TERMS ARE CALCULATED IN THE
C     DFTB_FMOGRAD SUBROUTINE.  LATTER BLOCK CALCULATES THE SHIFT
C     (\DELTA Q * \GAMMA) AND THIS WILL USED LATER IN DFTB_FMOGRAD.
C
C     ----- CALCULATE INTER-FRAGMENT GRADIENTS OF ENERGY -----
c
      if(dotruncesd) then
        ich=ixftch(x(lichfgg),ifg)
        jch=ixftch(x(lichfgg),jfg)
        if(ich.eq.0.and.jch.eq.0) then
          lcradfgi=lcradfg+(ifg-1)*4
          lcradfgj=lcradfg+(jfg-1)*4
          if(fmodest(x(lcradfgi),x(lcradfgi+3),
     *               x(lcradfgj),x(lcradfgj+3)).gt.rcutes) then
c           write(6,*) 'Skipping',ifg,lfg
            return
          endif
        endif
      endif
C
C     CALL VCLR(DE,1,3*NAT)
C
      NATI = NATFRG(IFG)
      NATJ = NATFRG(JFG)
      INDI = INDFRG(IFG)
      INDJ = INDFRG(JFG)
      ONE_THIRD = 1.0D+00/3.0D+00
C
      CALL DCOPY(10*NSPE*NSPE,GAMMAD,1,TMPD,1)
      IF (DFTB3) CALL DCOPY(11*NSPE*NSPE,GAMMAD3,1,TMPD3,1)
C
      NATFMOB = NATFMO + NBDFG
      NSEQ = 0
      DO I = 1, NATI
        IAG = IATFRG(INDI+I-1)
        ISP = ISPEFMO(IAG)
        CHAMUL1 = CHAMULI(I) !! CHARGE OF I FRAGMENT
        DO J = 1, NATJ
          JAG = IATFRG(INDJ+J-1)
          JSP = ISPEFMO(JAG)
          VALTMP = ZERO
          CHAMUL2 = CHAMULJ(J) !! CHARGE OF J FRAGMENT
          DAMPXH = DAMPXH0
          CALL DFTB_DMPCHK(ISP,JSP,DAMPXH)
          IF (PERIOD) THEN
            CALL DFTB_GAMDER_PBC(I,J,FMOC(1,IAG),FMOC(1,JAG),GDV2,
     *                           GDV3I,GDV3J,GDST2,GDST3I,GDST3J,
     *                           TMPD(1,ISP,JSP),TMPD3(1,ISP,JSP),
     *                           TMPD3(1,JSP,ISP),DAMPXH,ISP,JSP)
            !! Stress tensor
            VAL = CHAMUL1*CHAMUL2
            CALL DAXPY(9,SCAL*VAL,GDST2,1,STRFMO,1)
            IF (DFTB3) THEN
              VAL = ONE_THIRD*CHAMUL1*CHAMUL1*CHAMUL2
              CALL DAXPY(9,SCAL*VAL,GDST3I,1,STRFMO,1)
              VAL = ONE_THIRD*CHAMUL1*CHAMUL2*CHAMUL2
              CALL DAXPY(9,SCAL*VAL,GDST3J,1,STRFMO,1)
            END IF
            !! gradient
            IF (DFTB3) THEN
              GDV2(1) = GDV2(1)
     *          + ONE_THIRD*(CHAMUL1*GDV3I(1)+CHAMUL2*GDV3J(1))
              GDV2(2) = GDV2(2)
     *          + ONE_THIRD*(CHAMUL1*GDV3I(2)+CHAMUL2*GDV3J(2))
              GDV2(3) = GDV2(3)
     *          + ONE_THIRD*(CHAMUL1*GDV3I(3)+CHAMUL2*GDV3J(3))
            END IF
            VAL = CHAMUL1*CHAMUL2
            VEC(1) = GDV2(1)*VAL
            VEC(2) = GDV2(2)*VAL
            VEC(3) = GDV2(3)*VAL
          ELSE
            VEC(1) = FMOC(1,IAG) - FMOC(1,JAG)
            VEC(2) = FMOC(2,IAG) - FMOC(2,JAG)
            VEC(3) = FMOC(3,IAG) - FMOC(3,JAG)
            DIST = SQRT(VEC(1)*VEC(1)+VEC(2)*VEC(2)+VEC(3)*VEC(3))
            DISTI = ONE/DIST
            CALL DFTB_GAM23D(DIST,DISTI,GAMMAV,GAMI,GAMJ,
     *                       TMPD(1,ISP,JSP),TMPD3(1,ISP,JSP),
     *                       TMPD3(1,JSP,ISP),DAMPXH,DAMPXHE,SPE(ISP),
     *                       SPE(JSP),HUBDER(ISP),HUBDER(JSP),DFTB3)
            IF (DFTB3) GAMMAV = GAMMAV
     *                        + ONE_THIRD*(CHAMUL1*GAMI+CHAMUL2*GAMJ)
            VAL = CHAMUL1*CHAMUL2*GAMMAV*DISTI
            VEC(1) = VEC(1)*VAL
            VEC(2) = VEC(2)*VAL
            VEC(3) = VEC(3)*VAL
          END IF
          FMODE(1,IAG,2) = FMODE(1,IAG,2) + VEC(1)
          FMODE(2,IAG,2) = FMODE(2,IAG,2) + VEC(2)
          FMODE(3,IAG,2) = FMODE(3,IAG,2) + VEC(3)
          FMODE(1,JAG,2) = FMODE(1,JAG,2) - VEC(1)
          FMODE(2,JAG,2) = FMODE(2,JAG,2) - VEC(2)
          FMODE(3,JAG,2) = FMODE(3,JAG,2) - VEC(3)
          IF (NBODY.EQ.3) THEN
            FMODE(1,IAG,3) = FMODE(1,IAG,3) + VEC(1)*SCAL
            FMODE(2,IAG,3) = FMODE(2,IAG,3) + VEC(2)*SCAL
            FMODE(3,IAG,3) = FMODE(3,IAG,3) + VEC(3)*SCAL
            FMODE(1,JAG,3) = FMODE(1,JAG,3) - VEC(1)*SCAL
            FMODE(2,JAG,3) = FMODE(2,JAG,3) - VEC(2)*SCAL
            FMODE(3,JAG,3) = FMODE(3,JAG,3) - VEC(3)*SCAL
          END IF
        END DO
      END DO
      return
C
C     ----- PUT SHIFT -----
C     THE FIRST LOCATION OF SHIFTG IS X(LSHIFTG+NATFMO+NBDFG).
C     mysterious behavior, somehow i and j seem to be opposite?
C     However the gradient is correct, compared with ESP calculations.
C     Maybe definition of direction is differet.
C
      IF (IAND(MODGAMMA,3).EQ.0) THEN
        CALL DCOPY(12*NSPE*NSPE,X(LGAMMA),1,TMP,1)
        IF (DFTB3) CALL DCOPY(10*NSPE*NSPE,X(LGAMMA3),1,TMP3,1)
        NSEQ = 0
        DO I = 1, NATI
          TMPI = ZERO
          IAG = IATFRG(INDI+I-1)
          ISP = ISPEFMO(IAG)
          CHAMUL1 = CHAMULI(I) !! CHARGE OF I FRAGMENT
          DO J = 1, NATJ
            NSEQ = NSEQ + 1
            JAG = IATFRG(INDJ+J-1)
            JSP = ISPEFMO(JAG)
            CHAMUL2 = CHAMULJ(J) !! CHARGE OF J FRAGMENT
            DIST = DISTMAT(NSEQ)
            DAMPXH = DAMPXH0
            CALL DFTB_DMPCHK(ISP,JSP,DAMPXH)
            IF (DFTB3) THEN
              CALL DFTB_GAM23(1,DIST,GAMTMP,GAMI,GAMJ,
     *             TMP(1,ISP,JSP),TMP3(1,ISP,JSP),TMP3(1,JSP,ISP),
     *             DAMPXH,DAMPXHE,SPE(ISP),SPE(JSP),HUBDER(ISP),
     *             HUBDER(JSP),DFTB3)
              GAMTMPI = GAMTMP
     *                + ONE_THIRD*(TWO*CHAMUL1*GAMI+CHAMUL2*GAMJ)
              GAMTMPJ = GAMTMP
     *                + ONE_THIRD*(CHAMUL1*GAMI+TWO*CHAMUL2*GAMJ)
            ELSE
              DISTI = 0
              IF(DIST.NE.0) DISTI = ONE/DIST
              GAMTMPI = DFTB_GAM2(DIST,DISTI,TMP(1,ISP,JSP),DAMPXH,
     *          DAMPXHE,SPE(ISP),SPE(JSP))
              GAMTMPJ = GAMTMPI
            END IF
            TMPJ = GAMTMPJ*CHAMUL1
            TMPI = TMPI + GAMTMPI*CHAMUL2
            SHIFTG(INDJ+J-1) = SHIFTG(INDJ+J-1) + TMPJ*SCAL
          END DO
          SHIFTG(INDI+I-1) = SHIFTG(INDI+I-1) + TMPI*SCAL
        END DO
      ELSE IF (IAND(MODGAMMA,1).NE.0) THEN
        DO J = 1, NATJ
          TMPJ = ZERO
          JAG = IATFRG(INDJ+J-1)
          JSP = ISPEFMO(JAG)
          CHAMUL2 = CHAMULJ(J) !! CHARGE OF J FRAGMENT
          DO I = 1, NATI
            IAG = IATFRG(INDI+I-1)
            CHAMUL1 = CHAMULI(I) !! CHARGE OF I FRAGMENT
            CALL DFTB_CNVSQ(IAG,JAG,NSEQ)
            GAMTMPI = X(LGAMMA+NSEQ-1)
            GAMTMPJ = GAMTMPI
            IF (DFTB3) THEN
              GAMI = X(LGAMMA3+IAG-1+(JAG-1)*NATFMO) !! G_ac
              GAMJ = X(LGAMMA3+JAG-1+(IAG-1)*NATFMO) !! G_ca
              GAMTMPI = GAMTMPI
     *                + ONE_THIRD*(TWO*CHAMUL1*GAMI+CHAMUL2*GAMJ)
              GAMTMPJ = GAMTMPJ
     *                + ONE_THIRD*(CHAMUL1*GAMI+TWO*CHAMUL2*GAMJ)
            END IF
            TMPI = GAMTMPI*CHAMUL2
            TMPJ = TMPJ + GAMTMPJ*CHAMUL1
            SHIFTG(INDI+I-1) = SHIFTG(INDI+I-1) + TMPI*SCAL
          END DO
          SHIFTG(INDJ+J-1) = SHIFTG(INDJ+J-1) + TMPJ*SCAL
        END DO
      ELSE IF (IAND(MODGAMMA,2).NE.0) THEN
        NATFMOB = NATFMO + NBDFG
        !! POTENTIAL FOR I-TH FRAGMENT
        NSEQ = INDJ-1 + NATFMOB*(INDI-1)
        DO I = 1, NATI
          SHIFTG(INDI+I-1) = SHIFTG(INDI+I-1)
     *      + DDOT(NATJ,CHAMULJ,1,X(LGAMMA+NSEQ),1)*SCAL
          NSEQ = NSEQ + NATFMOB
        END DO
        IF (DFTB3) THEN
          DO I = 1, NATI
            II = INDI+I-1
            SHIFTG(INDI+I-1) = SHIFTG(INDI+I-1) + ONE_THIRD*(
     *        TWO*CHAMULI(I)*DDOT(NATJ,CHAMULJ,1,
     *                                 GAMMA3(2,INDJ,II),2)
     *        + DDOT(NATJ,CHAMULJ(NATJ+1),1,GAMMA3(1,INDJ,II),2))*SCAL
          END DO
        END IF
        !! POTENTIAL FOR J-TH FRAGMENT
        NSEQ = INDI-1 + NATFMOB*(INDJ-1)
        DO J = 1, NATJ
          SHIFTG(INDJ+J-1) = SHIFTG(INDJ+J-1)
     *      + DDOT(NATI,CHAMULI,1,X(LGAMMA+NSEQ),1)*SCAL
          NSEQ = NSEQ + NATFMOB
        END DO
        IF (DFTB3) THEN
          DO J = 1, NATJ
            JJ = INDJ+J-1
            SHIFTG(INDJ+J-1) = SHIFTG(INDJ+J-1) + ONE_THIRD*(
     *        DDOT(NATI,CHAMULI(NATI+1),1,GAMMA3(1,INDI,JJ),2)
     *        + TWO*CHAMULJ(J)*DDOT(NATI,CHAMULI,1,
     *                                   GAMMA3(2,INDI,JJ),2))*SCAL
          END DO
        END IF
      END IF
C
      RETURN
C
      END SUBROUTINE DFTB_ESDGRAD
C
C-----------------------------------------------------------------------
C*MODULE DFTBFO    *DECK ADD_LAG_GRAD
C>
C>    @brief Gradient in FMO-DFTB
C>
C>    @details Add Lagrangian-like term contribution in FMO-DFTB.
C>
C>    @author Yoshio Nishimoto
C>    - Sep, 2014- Subroutine written
C>    @date   Dec, 2014- Yoshio Nishimoto
C>    - Unnecessary arguments are removed
C>    - Performance improvement
C>
C>           --- INPUT ---
C>    @param L1     Number of basis functions
C>    @param L2     L1*(L1+1)/2
C>    @param NAT    Number of atoms
C>    @param NATFRG Number of atoms in each fragment
C>    @param IAGLOB Global atom index
C>    @param DSINT  Maybe dS_{\mu \nu}/da
C>    @param DSD    Lagrangian-like term
C>    @param WRK    Maybe working array
C>    @param ITYPE  Decide coefficients
C>           --- IN/OUTPUT ---
C>    @param FMODE  FMO2 gradient
C>
C
      SUBROUTINE ADD_LAG_GRAD(L1,L2,NAT,NATFMO,IAGLOB,FMODE,
     *                        DSINT,DSD,WRK)
      use mx_limits, only: mxatm
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER (MXSPE=10)
      COMMON /DFTBPR/ ETEMP,DFTBDP(MXSPE*14),DAMPXHE,HUBDER(MXSPE),
     *                ZREF(MXATM),SPNCST(6,MXSPE),SPE(MXATM),NSPE,
     *                MAXANG(MXATM),ISPE(MXATM),IND(MXATM+1),IDFTBD,
     *                PARAMDIR
C
      DOUBLE PRECISION :: DSINT(L2,3,*),DSD(L1,L1),FMODE(3,NATFMO),
     * WRK(*)
      INTEGER :: IAGLOB(*)
C
C     ----- CALCULATE (DSD)*(DS/DR) AND ADD GRADIENT -----
C
C     FIRST SYMMETRIZE DSD MATRIX, SCALING BY TWO
      NSEQ = 0
      DO I = 1, L1
        DO J = 1, I
          NSEQ = NSEQ + 1
          WRK(NSEQ) = DSD(I,J) + DSD(J,I)
        END DO
      END DO
C
C     THEN ADD TO GRADIENT BY CALCULATING (DSD)*(dS/dR)
      DO II = 1, NAT
        IG = IAGLOB(II)
        SUM1 = 0.0D+00
        SUM2 = 0.0D+00
        SUM3 = 0.0D+00
        CALL DDOT_DER(NAT,II,L2,IND,SUM1,WRK,DSINT(1,1,II))
        CALL DDOT_DER(NAT,II,L2,IND,SUM2,WRK,DSINT(1,2,II))
        CALL DDOT_DER(NAT,II,L2,IND,SUM3,WRK,DSINT(1,3,II))
        FMODE(1,IG) = FMODE(1,IG) + SUM1
        FMODE(2,IG) = FMODE(2,IG) + SUM2
        FMODE(3,IG) = FMODE(3,IG) + SUM3
      END DO
C
      RETURN
C
      END SUBROUTINE ADD_LAG_GRAD
C
C-----------------------------------------------------------------------
C*MODULE DFTBFO    *DECK GAMD_FMO2
C>
C>    @brief Gradient in FMO-DFTB
C>
C>    @details Calculate derivatives of DFTB2 gamma values
C>
C>    @author Yoshio Nishimoto
C>    - Sep, 2014- Subroutine written
C>
C>           --- INPUT ---
C>    @param DIST Distance between atom A and B
C>    @param VALUE Constants
C>    @param DAMPXH Whether X-H damping is used or not
C>    @param DAMPXHE Exponent of X-H damping
C>    @param SPE1,SPE2 Name of atom A and B
C>
C
      DOUBLE PRECISION FUNCTION GAMD_FMO2(DIST,DISTI,VALUE,DAMPXH,
     * DAMPXHE,SPE1,SPE2)
C
      IMPLICIT NONE
C
      DOUBLE PRECISION,INTENT(IN) :: DIST,DISTI,VALUE(10),DAMPXHE,
     *  SPE1,SPE2
      LOGICAL,INTENT(IN) :: DAMPXH
C
      DOUBLE PRECISION :: TAUA,TAUB,TAUMEAN,EXP_UR,
     * DISTI2,DIST2,HUBMEAN,GA,DGA_DR,DS_DR,GAMTMP,
     * GAMTMP2,HUBMEANEXP,EXPH,EXPA,EXPB,FAB,FBA,DFAB_DR,DFBA_DR
      DOUBLE PRECISION, PARAMETER :: DISTTOL=1.0D-08,HUBTOL=1.0D-04,
     * TWO=2.0D+00
      DOUBLE PRECISION H/8HH       /
C
C     --- Here calculates THE DERIVATIVE OF \gamma_{\alpha \beta} values ---
C
      TAUA = VALUE(1)*1.0D+00
      TAUB = VALUE(2)*1.0D+00
      GAMTMP=0.0D+00
      GAMTMP2=0.0D+00
C     write (*,*) "taua, taub = ", taua,taub
      IF (DIST.LE.DISTTOL) THEN
        GAMD_FMO2 = 0.0D+00
        RETURN
      ELSE IF (ABS(TAUA-TAUB).LE.HUBTOL) THEN
        TAUMEAN = VALUE(3)
        EXP_UR = EXP(-TAUMEAN*DIST)
        DIST2 = DIST*DIST
C       DISTI = 1.0D+00/DIST
        DISTI2 = DISTI*DISTI
        GA = DISTI+VALUE(4)+VALUE(5)*DIST+VALUE(6)*DIST2
        DGA_DR = -DISTI2+VALUE(5)+TWO*VALUE(6)*DIST
        DS_DR = EXP_UR*DGA_DR - TAUA*EXP_UR*GA
        GAMTMP = -DS_DR
        IF (DAMPXH) THEN
          IF (SPE1.EQ.H.OR.SPE2.EQ.H) THEN
            HUBMEAN = TAUA*3.125D-01
            HUBMEANEXP = HUBMEAN**DAMPXHE
            EXPH = EXP(-HUBMEANEXP*DIST2)
            GAMTMP = GAMTMP*EXPH
            GAMTMP2 = EXP_UR*GA
            GAMTMP2 = GAMTMP2*HUBMEANEXP*TWO*DIST*EXPH
          END IF
        END IF
      ELSE
C       DISTI = 1.0D+00/DIST
        DISTI2 = DISTI*DISTI
        EXPA = EXP(-TAUA*DIST)
        EXPB = EXP(-TAUB*DIST)
        FAB = VALUE(9) - VALUE(7)*DISTI
        FBA = VALUE(10) - VALUE(8)*DISTI
        DFAB_DR = VALUE(7)*DISTI2
        DFBA_DR = VALUE(8)*DISTI2
        DS_DR = EXPA*DFAB_DR - TAUA*EXPA*FAB
     *        + EXPB*DFBA_DR - TAUB*EXPB*FBA
        GAMTMP = -DS_DR
        IF (DAMPXH) THEN
          IF (SPE1.EQ.H.OR.SPE2.EQ.H) THEN
            DIST2 = DIST*DIST
            HUBMEAN = (TAUA+TAUB)*1.5625D-01
            HUBMEANEXP = HUBMEAN**DAMPXHE
            EXPH = EXP(-HUBMEANEXP*DIST2)
            GAMTMP = GAMTMP * EXPH !! -h*dSf/dR
            GAMTMP2 = EXPA*FAB + EXPB*FBA
            GAMTMP2 = GAMTMP2*HUBMEANEXP*TWO*DIST*EXPH
          END IF
        END IF
      END IF
      GAMD_FMO2 = -DISTI2 + GAMTMP + GAMTMP2
C     WRITE (*,'(2I4,2F20.10)') I,J,GAMD_FMO,-1.0/(DIST**2)
C
      RETURN

      END FUNCTION GAMD_FMO2
C
C-----------------------------------------------------------------------
C*MODULE DFTBFO    *DECK DS_SHIFT
C>
C>    @brief Gradient in FMO-DFTB
C>
C>    @details Add shift contribution toe array D*S
C>
C>    @author Yoshio Nishimoto
C>    - Sep, 2014- Subroutine written
C>    @date   Dec, 2014- Yoshio Nishimoto
C>    - Performance improvement
C>
C>           --- INPUT ---
C>    @param L1 Number of basis functions
C>    @param L2 L1*(L1+1)/2
C>    @param NSTART,NLAST Lower and upper bound to be calculated
C>    @param IAGLOB Global atom index
C>    @param SHIFT Shiftcontribution
C>           --- IN/OUTPUT ---
C>    @param DS (density matrix)*(overlap)
C>
C
      SUBROUTINE DS_SHIFT(L1,NAT,DS,SHIFT)
      use mx_limits, only: mxatm
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION DS(L1,L1),SHIFT(*)
C
      PARAMETER (MXSPE=10)
      COMMON /DFTBPR/ ETEMP,DFTBDP(MXSPE*14),DAMPXHE,HUBDER(MXSPE),
     *                ZREF(MXATM),SPNCST(6,MXSPE),SPE(MXATM),NSPE,
     *                MAXANG(MXATM),ISPE(MXATM),IND(MXATM+1),IDFTBD,
     *                PARAMDIR
C
C     THIS SUBROUTINE ADDS THE SHIFT CONTRIBUTION TO THE ARRAY DS.
C     NOTE THAT SHIFT IS THE CONTRIBUTION ON AN ATOM. 
C
      K = 1
      DO I = 1, NAT
        NVEC = IND(I+1)-IND(I)
        CALL DSCAL(L1*NVEC,SHIFT(I),DS(1,K),1)
        K = K + NVEC
      END DO
C
      RETURN
C
      END SUBROUTINE DS_SHIFT
C
C-----------------------------------------------------------------------
C*MODULE DFTBFO    *DECK DFTB_ZVLAGGET
C>
C>    @brief   Fully analytic FMO-DFTB gradient
C>
C>    @details Get Lagrangian for FMO-DFTB SCZV
C>
C>    @author  Yoshio Nishimoto
C>    - Feb, 2016- Subroutine written
C>
C>           --- INPUT ---
C>    @param SHIFTCT Shift due to charge transfer
C>    @param VEC     Orbitals
C>    @param S       Overlap in AO
C>    @param L1      Number of AOs
C>    @param NAT     Number of atoms
C>    @param IPTLG   Some index
C>    @param ZVWRK   Working arrays?
C>    @param NOCC    Number of occupied orbitals
C>    @param NVIR    Number of virtual orbitals
C>    @param IFG     Current fragment
C>           --- OUTPUT ---
C>    @param ZVLAG   Calculated lagrangian?
C>
C
      SUBROUTINE DFTB_ZVLAGGET(SHIFTCT,VEC,S,L1,NAT,IPTLG,ZVWRK,ZVLAG,
     *                         NOCC,NVIR,IFG)
      use mx_limits, only: mxatm
C      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER (MXSPE=10)
      DIMENSION SHIFTCT(*), VEC(L1,*), S(*)
      DIMENSION IPTLG(*),ZVWRK(NVIR,NOCC),ZVLAG(*)
      LOGICAL   dodistr
C
      COMMON /DFTBPR/ ETEMP,DFTBDP(MXSPE*14),DAMPXHE,HUBDER(MXSPE),
     *                ZREF(MXATM),SPNCST(6,MXSPE),SPE(MXATM),NSPE,
     *                MAXANG(MXATM),ISPE(MXATM),IND(MXATM+1),IDFTBD,
     *                PARAMDIR
      common /fmoopt/ espsca(9),RESPAP(2),rESPPC(2),rESDIM,restri(4),
     *                rcorsd,respct,convfg,cnvdmp,coroff,rflmo(4),
     *                orshft,orshft2,cnvafo,ascreen(4),IXESP,mxitfg,
     *                nguess,NBSSE,modorb,modpar,irststp,irstlay,nprfmo,
     *                nfmopal,modprp,maxl1c,ipieda,modgrd,modesp,ivmul,
     *                modlmo,nopden,mofock,modfd,modfmm,ncentm,ndualb,
     *                ngab,modpan
C
      dodistr=iand(modpar,512+1024).ne.0
C
      DO MO = 1, NVIR
        DO IO = 1, NOCC
          TMP = 0.0D+00
          DO NC = 1, NAT
            TMP2 = 0.0D+00
            DO M = 1, IND(NC+1)-IND(NC)
              MU = IND(NC)+M
              DO NU = 1, L1
                CALL DFTB_CNVSQ(MU,NU,NSEQ)
                TMP2 = TMP2 + (VEC(MU,MO+NOCC)*VEC(NU,IO)
     *                      +  VEC(MU,IO)*VEC(NU,MO+NOCC))*S(NSEQ)
              END DO
            END DO
            TMP = TMP + SHIFTCT(NC)*TMP2
          END DO
          ZVWRK(MO,IO) = TMP*0.5D+00
        END DO
      END DO
C
      CALL ZVLGMULT(NOCC*NVIR,1,IFG,IPTLG,DODISTR,ZVLAG,ZVWRK,VEC)
C
      RETURN
C
      END SUBROUTINE DFTB_ZVLAGGET
C
C-----------------------------------------------------------------------
C*MODULE DFTBFO    *DECK DFTB_DIAORBHSS
C>
C>    @brief   Fully analytic FMO-DFTB gradient
C>
C>    @details Calculate diagonal block of orbital rotation Hessian
C>
C>    @author  Yoshio Nishimoto
C>    - Feb, 2016- Subroutine written
C>    @date    Feb, 2018 - Yoshio Nishimoto
C>    - Added an option to omit doing for dimers/trimers
C>
C>           --- INPUT ---
C>    @param L2      (NUM*NUM+NUM)/2
C>    @param SHIFT   Shift due to induced vector?
C>    @param S       Overlap in AO
C>    @param Q       Induced charges
C>    @param GAMMA   DFTB2 gamma
C>    @param GAMMA3  DFTB3 gamma
C>    @param DFTB3   DFTB or not
C>    @param CHAMUL  Mulliken population
C>    @param IFG     Current fragment
C>    @param IAGLOB  Global index of atom
C>    @param NATFRG  Number of atoms in fragment
C>    @param INDFRG  First atom of each fragment in terms of NATFMOB?
C>    @param IATFRG  Global atom index with NBDFG?
C>    @param POPMAT  Matrix of Mulliken population
C>    @param FMOC    Coordinates of the whole system
C>    @param INDAT   Converts global atomic location to fragment index
C>    @param IALOC   Converts global to local (fragment) atomic location
C>           --- OUTPUT ---
C>    @param WRT1    Orbital rotation Hessian in AO
C>
C
      SUBROUTINE DFTB_DIAORBHSS(L1,L2,SHIFT,S,Q,GAMMA,GAMMA3,WRT1,DFTB3,
     *  CHAMUL,IFG,IAGLOB,NATFRG,INDFRG,IATFRG,POPMAT,FMOC,INDAT,IALOC,
     *  ZZ,LCDFTB,GAMLC,SMAX,INTLST)
      use mx_limits, only: mxatm
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER (MXSPE=10)
      COMMON /DFTBPR/ ETEMP,DFTBDP(MXSPE*14),DAMPXHE,HUBDER(MXSPE),
     *                ZREF(MXATM),SPNCST(6,MXSPE),SPE(MXATM),NSPE,
     *                MAXANG(MXATM),ISPE(MXATM),IND(MXATM+1),IDFTBD,
     *                PARAMDIR
      COMMON /FMCOM / X(1)
      common /fmodim/ maxbnd,maxknd,maxcbs,maxcao,maxvec,maxl1,maxnat,
     *                maxabd,maxbas,maxbbd,maxlmo,maxslo,maxabd2,maxrij
      COMMON /FMODTB/ LGAMMA,LZREF,LISPE,LSHIFTG,LSHIFTCT,
     *                LGAMMAD,MODESD,MODGAMMA,LDCTMP,LGAMMA3,LGAMMA3D,
     *                LNBOND,LCTMUL,LESPDFTB,LESPASC,NDFTBRST,LDFTBRST
      COMMON /FMOINF/ NFG,NLAYER,NATFMO,NBDFG,NAOTYP,NBODY,NSEGM
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /OUTPUT/ NPRINT,ITOL,ICUT,NORMF,NORMP,NOPK
C
      INTEGER, INTENT(IN) :: L1,L2,IFG,IAGLOB(*),NATFRG(*),INDFRG(*),
     *  IATFRG(*),INDAT(*),IALOC(*),INTLST(*)
      DOUBLE PRECISION, INTENT(IN) :: S(L2),Q(NAT),GAMMA(*),
     *  GAMMA3(NAT,NAT),CHAMUL(NAT),POPMAT(MAXNAT,*),FMOC(3,*),ZZ(L2),
     *  SMAX(*),GAMLC(*)
      DOUBLE PRECISION, INTENT(INOUT) :: SHIFT(NAT),WRT1(L2)
      LOGICAL, INTENT(IN) :: DFTB3
      logical LCDFTB
C
C     THIS SUBROUTINE SHOULD RETURN 1/4 ORBIAL HESSIAN
C
C     ----- CONSTRUCT SHIFT MATRIX WITH QINDZ CHARGES -----
C     NOTE THAT QIND = \sum \sum (czc + czc)*S = 1/2 Q_C
C
      IF (DFTB3) THEN
        IF (IAND(MODGAMMA,4096).EQ.0) THEN
          !! INTER-FRAGMENT INTERACTION OF 3RD ORDER
          IF (IAND(MODGAMMA,8).EQ.0) THEN
            IF (IAND(MODGAMMA,2).NE.0) THEN
              NN = 2
              NN2 = NATFMO+NBDFG
            ELSE
              NN = 1
              NN2 = NATFMO
            END IF
            CALL DFTB_OFFORBHSS(0,IAGLOB,NATFRG,INDFRG,IATFRG,
     *           POPMAT,S,WRT1,X(LISPE),FMOC,X(LGAMMA3),
     *           SHIFT,Q,POPMAT,.FALSE.,INDAT,IALOC,NN,NN2)
            !! INTRA-FRAGMENT INTERACTION OF 3RD ORDER
            DO I = 1, NAT
              DO J = 1, NAT
                GAMI = GAMMA3(I,J)
                GAMJ = GAMMA3(J,I)
                SHIFT(I) = SHIFT(I)
     *            + 2.0D+00/3.0D+00*Q(J)*CHAMUL(I)*GAMI
     *            + 2.0D+00/3.0D+00*Q(J)*CHAMUL(J)*GAMJ
              END DO
            END DO
          ELSE
            !! INTRA-FRAGMENT INTERACTION OF 3RD ORDER
            INDI = INDFRG(IFG)
            NATFMOB = NATFMO + NBDFG
            DO I = 1, NAT
              SHIFT(I) = 0.0D+00
              DO J = 1, NAT
                GAMI = GAMMA3(I,J)
                GAMJ = GAMMA3(J,I)
                SHIFT(I) = SHIFT(I)
     *            + 2.0D+00/3.0D+00*Q(J)*CHAMUL(I)*GAMI
     *            + 2.0D+00/3.0D+00*Q(J)*CHAMUL(J)*GAMJ
     *            + 2.0D+00/3.0D+00*Q(I)*CHAMUL(J)*GAMI
              END DO
              !! INTER-FRAGMENT INTERACTION OF 3RD ORDER
              SHIFT(I) = SHIFT(I)
     *          + 2.0D+00*Q(I)*X(LESPDFTB+2*NATFMOB+INDI+I-2)
            END DO
          END IF
C      write (*,*) "shift = ", shift
        ELSE
          CALL VCLR(SHIFT,1,NAT)
          CALL DFTB_TD_3RD_SHIFT(SHIFT,SHIFT,Q,Q,
     *                           CHAMUL,GAMMA3,NAT,.TRUE.)
        END IF
      ELSE
        CALL VCLR(SHIFT,1,NAT)
      END IF
      !! INTRA-FRAGMENT INTERACTION OF 2ND ORDER
      CALL DFTB_2ND_SHIFT(SHIFT,Q,GAMMA,NAT)
C
C     ----- CONSTRUCT WRT1 MATRIX (2ND TERM OF EQ.(8)) -----
C
      !! DO LC-DFTB FIRST
      CALL VCLR(WRT1,1,L2)
      IF (LCDFTB) THEN
        IF (ICUT.GE.0) THEN
          CALL VALFM(LOADFM)
          LWRK1  = LOADFM + 1
          LWRK2  = LWRK1  + L1*L1
          LSUMDDA= LWRK2  + L2
          LWRK3  = LSUMDDA+ NAT
          LSSQ   = LWRK3  + NAT*NAT
          LAST   = LSSQ   + L1*L1
          NEED   = LAST - LOADFM - 1
          CALL GETFM(NEED)
          CALL CPYTSQ(S,X(LSSQ),L1,1)
          CALL DFTB_LCSHIFT(1,NAT,L1,GAMLC,X(LSSQ),X(LAST),
     *                      X(LAST),ZZ,X(LWRK1),X(LAST),WRT1,
     *                      X(LWRK2),SMAX,X(LSUMDDA),
     *                      X(LWRK3),ELC,.FALSE.,ICUT,INTLST)
          CALL RETFM(NEED)
        ELSE
          CALL DFTB_LCSHIFT2(1,NAT,1,L1,ISPE,MAXANG,IND,ZZ,DUM,
     *                       S,GAMLC,WRT1,ELC)
        END IF
      END IF
      CALL DFTB_SHIFT_FOCK(SHIFT,WRT1,S,NAT,L2,NSPE,ISPE,IND,MAXANG)
C
      RETURN
C
      END SUBROUTINE DFTB_DIAORBHSS
C
C-----------------------------------------------------------------------
C*MODULE DFTBFO    *DECK DFTB_OFFORBHSS
C>
C>    @brief   Fully analytic FMO-DFTB gradient
C>
C>    @details Calculate off-diagonal block of orbital rotation Hessian
C>
C>    @author  Yoshio Nishimoto
C>    - Feb, 2016- Subroutine written
C>
C>           --- INPUT ---
C>    @param IFG     Current fragment
C>    @param IAGLOB  Global index of atom
C>    @param NATFRG  Number of atoms in fragment
C>    @param INDFRG  First atom of each fragment in terms of NATFMOB?
C>    @param IATFRG  Global atom index with NBDFG?
C>    @param QINDZ   Induced charges of Z
C>    @param SS      Overlap in AO
C>    @param ISPEFMO ISPE for all atoms (NATFMOB length)
C>    @param FMOC    Coordinates of the whole system
C>    @param GAMMA3  DFTB3 gamma
C>    @param SHIFTOUT Some shift contributions used in gradient
C>    @param CHAMULIFG Mulliken population of IFG
C>    @param DFTB2   DFTB2 or not
C>    @param INDAT   Converts global atomic location to fragment index
C>    @param IALOC   Converts global to local (fragment) atomic location
C>    @param NN      Species the dimension of array
C>    @param NN2     Species the dimension of array
C>           --- OUTPUT ---
C>    @param FAO     Orbital rotation Hessian in AO
C>
C
      SUBROUTINE DFTB_OFFORBHSS(IFG,iaglob,natfrg,indfrg,iatfrg,QINDZ,
     *  ss,fao,ISPEFMO,FMOC,GAMMA3,SHIFTOUT,CHAMULIFG,POPMAT,DFTB2,
     *  INDAT,IALOC,NN,NN2)
      use mx_limits, only: mxatm
      USE DFTBPB_MOD,ONLY: PERIOD,DFTB_GAMFMO_PBC
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      INTEGER,PARAMETER :: MXSPE=10
      LOGICAL DFTBFL,SCC,SRSCC,DFTB3,DAMPXH0,LCDFTB,dotrunces,dotruncesd
      DOUBLE PRECISION, PARAMETER :: ZERO=0.0D+00,TWO=2.0D+00,
     * ONE_THIRD=1.0D+00/3.0D+00
C
      COMMON /DFTB  / DFTBFL,SCC,SRSCC,DFTB3,DAMPXH0,LCDFTB
      COMMON /DFTBPR/ ETEMP,DFTBDP(MXSPE*14),DAMPXHE,HUBDER(MXSPE),
     *                ZREF(MXATM),SPNCST(6,MXSPE),SPE(MXATM),NSPE,
     *                MAXANG(MXATM),ISPE(MXATM),IND(MXATM+1),IDFTBD,
     *                PARAMDIR
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /FMCOM / X(1)
      COMMON /FMOINF/ NFG,NLAYER,NATFMO,NBDFG,NAOTYP,NBODY,NSEGM
      COMMON /FMODTB/ LGAMMA,LZREF,LISPE,LSHIFTG,LSHIFTCT,
     *                LGAMMAD,MODESD,MODGAMMA,LDCTMP,LGAMMA3,LGAMMA3D,
     *                LNBOND,LCTMUL,LESPDFTB,LESPASC,NDFTBRST,LDFTBRST
      common /fmodim/ maxbnd,maxknd,maxcbs,maxcao,maxvec,maxl1,maxnat,
     *                maxabd,maxbas,maxbbd,maxlmo,maxslo,maxabd2,maxrij
      common /cenrad/ rcutes,lcradfg,lichfgg,dotrunces,dotruncesd
C 
      double precision QINDZ(maxnat,nfg),ss(*),fao(*),CHAMULIFG(*),
     *  POPMAT(MAXNAT,NFG)
C     integer layfrg(*)
C
      INTEGER :: IAGLOB(*),NATFRG(*),INDFRG(*),IATFRG(*),ISPEFMO(*),
     *  INDAT(*),IALOC(*)
      DOUBLE PRECISION :: SHIFT(NAT),FMOC(3,NATFMO),TMP(12,NSPE,NSPE),
     * SHIFT2(NAT*2),TMP3(10,NSPE,NSPE),GAMMA3(NN2,NN2,*),
     * SHIFTOUT(NAT)
      DOUBLE PRECISION :: CHATMP(MAXNAT*2)
      LOGICAL :: DFTB2,DAMPXH
C
C     ------ CALCULATE OFF-DIAGONAL A_{ij,kl}^{I,K} MATRIX -----
C     THIS SUBROUTINE SHOULD RETURN 1/4 ORBIAL HESSIAN
C     QINDZ CONTAINS ONLY 1/2(czc + czc)*S
C
      CALL VCLR(SHIFT,1,NAT)
      IF (DFTB3) CALL VCLR(SHIFT2,1,NAT*2)
      IF (IAND(MODGAMMA,4096).NE.0) THEN
        CALL VCLR(SHIFTOUT,1,NAT)
        CALL VCLR(FAO,1,NUM*(NUM+1)/2)
        RETURN
      END IF
C
C     LOOP FOR ENVIRONMENTAL FRAGMENTS
C     THIS SUBROUTINE IS MOSTLY A COPY OF DFTB_ESP
C
      IF (IAND(MODGAMMA,3).EQ.0) THEN
        CALL DCOPY(12*NSPE*NSPE,X(LGAMMA),1,TMP,1)
        IF (DFTB3)  CALL DCOPY(10*NSPE*NSPE,X(LGAMMA3),1,TMP3,1)
        DO LFG = 1, NFG
          IF (LFG.EQ.IFG) CYCLE
          if(dotrunces.and.ich.eq.0) then
            lcradfgi=lcradfg+(ifg-1)*4
            lcradfgl=lcradfg+(lfg-1)*4
            lch=ixftch(x(lichfgg),lfg)  
            if(lch.eq.0.and.fmodest(x(lcradfgi),x(lcradfgi+3),
     *         x(lcradfgl),x(lcradfgl+3)).gt.rcutes/10) then
c             write(6,*) 'Skipping',ifg,lfg
              cycle
            endif
          endif
          NATJ = NATFRG(LFG)
          INDJ = INDFRG(LFG)
          JJ = INDJ - 1
          DO J = 1, NATJ !! ENVIRONMENTAL FRAGMENTS
            JJ = JJ + 1
            CHAMUL = QINDZ(J,LFG)
            IF (IFG.EQ.0) CHAMUL = POPMAT(J,LFG)-X(LZREF+INDJ+J-2)
            IF (CHAMUL.EQ.ZERO) CYCLE
            JAG = IATFRG(JJ)
            JSP = ISPEFMO(JAG)
            DO I = 1, NAT !! FRAGMNET TO BE CONSIDERED
              IAG = IAGLOB(I)
              ISP = ISPEFMO(IAG)
              DAMPXH = DAMPXH0
              CALL DFTB_DMPCHK(ISP,JSP,DAMPXH)
              IF (PERIOD) THEN
                CALL DFTB_GAMFMO_PBC(IAG,JAG,FMOC(1,IAG),FMOC(1,JAG),
     *                               GAMMAV,GAMI,GAMJ,TMP(1,ISP,JSP),
     *                               TMP3(1,ISP,JSP),TMP3(1,JSP,ISP),
     *                               DAMPXH,ISP,JSP)
              ELSE
               XX = FMOC(1,IAG)-FMOC(1,JAG)
               YY = FMOC(2,IAG)-FMOC(2,JAG)
               ZZ = FMOC(3,IAG)-FMOC(3,JAG)
               DIST = XX*XX + YY*YY + ZZ*ZZ
               CALL DFTB_GAM23(0,DIST,GAMMAV,GAMI,GAMJ,
     *              TMP(1,ISP,JSP),TMP3(1,ISP,JSP),TMP3(1,JSP,ISP),
     *              DAMPXH,DAMPXHE,SPE(ISP),SPE(JSP),HUBDER(ISP),
     *              HUBDER(JSP),DFTB3)
              END IF
              SHIFT(I) = SHIFT(I) + CHAMUL * GAMMAV
              IF (DFTB3) THEN
                CHAMULI = CHAMULIFG(I)
                !! MULTIPLY 1/3 LATER WHEN ADDED TO SHIFT
                IF (.NOT.DFTB2) THEN
                  SHIFT2(I) = SHIFT2(I)
     *              + CHAMUL*CHAMULI*GAMI*TWO
                ELSE IF (IFG.EQ.0) THEN
                  SHIFT2(I) = SHIFT2(I)
     *              + CHAMUL*(CHAMULI*GAMI*TWO+CHAMUL*GAMJ)
                ELSE
                  CHAMULJ = POPMAT(J,LFG)-X(LZREF+INDJ+J-2)
                  SHIFT2(I) = SHIFT2(I)
     *              + TWO*CHAMUL*(CHAMULI*GAMI+CHAMULJ*GAMJ)
                END IF
              END IF
            END DO
          END DO
        END DO
      END IF
C
      IF (IAND(MODGAMMA,1).NE.0) THEN
        DO LFG = 1, NFG
          IF (LFG.EQ.IFG) CYCLE
          NATJ = NATFRG(LFG)
          INDJ = INDFRG(LFG)
          JJ = INDJ - 1
          DO J = 1, NATJ !! ENVIRONMENTAL FRAGMENTS
            JJ = JJ + 1
            CHAMUL = QINDZ(J,LFG)
            IF (IFG.EQ.0) CHAMUL = QINDZ(J,LFG)-X(LZREF+INDJ+J-2)
            IF (CHAMUL.EQ.ZERO) CYCLE
            JAG = IATFRG(JJ)
            DO I = 1, NAT !! FRAGMNET TO BE CONSIDERED
              IAG = IAGLOB(I)
              IF (IAG.GE.JAG) THEN
                NSEQ = IAG*(IAG-1)/2 + JAG
              ELSE
                NSEQ = JAG*(JAG-1)/2 + IAG
              END IF
              GAMMAV = X(LGAMMA+NSEQ-1)
              SHIFT(I) = SHIFT(I) + CHAMUL * GAMMAV
              IF (DFTB3) THEN
                CHAMULI = CHAMULIFG(I)
                GAMI = GAMMA3(IAG,JAG,1)
                GAMJ = GAMMA3(JAG,IAG,1)
                IF (.NOT.DFTB2) THEN
                  SHIFT2(I) = SHIFT2(I)
     *              + CHAMUL*CHAMULI*GAMI*TWO
                ELSE IF (IFG.EQ.0) THEN
                  SHIFT2(I) = SHIFT2(I)
     *              + CHAMUL*(CHAMULI*GAMI*TWO+CHAMUL*GAMJ)
                ELSE
                  CHAMULJ = POPMAT(J,LFG)-X(LZREF+INDJ+J-2)
                  SHIFT2(I) = SHIFT2(I)
     *              + TWO*CHAMUL*(CHAMULI*GAMI+CHAMULJ*GAMJ)
                END IF
              END IF
            END DO
          END DO
        END DO
      END IF
C
      IF (IAND(MODGAMMA,2).NE.0) THEN
        natfmob = natfmo + nbdfg
        nati = natfrg(ifg)
        indi = indfrg(ifg)
        call dgemv('t',nn2,nati,1.0d+00,x(lgamma+nn2*(indi-1)),nn2,
     *             x(lshiftg),1,0.0d+00,shift,1)
        if (dftb3) then
          call dgemv('t',nn2,nati,2.0d+00,gamma3(1,indi,2),nn2,
     *               x(lshiftg),1,0.0d+00,shift2(nat+1),1)
          call dgemv('t',nn2,nati,2.0d+00,gamma3(1,indi,1),nn2,
     *               x(lshiftg+nn2),1,0.0d+00,shift2,1)
        end if
        call dgemv('t',nati,nati,-1.0d+00,
     *             x(lgamma+indi-1+nn2*(indi-1)),nn2,
     *             x(lshiftg+indi-1),1,1.0d+00,shift,1)
        IF (DFTB3) THEN
          call dgemv('t',nati,nati,-2.0d+00,gamma3(indi,indi,2),nn2,
     *               x(lshiftg+indi-1),1,1.0d+00,shift2(nat+1),1)
          call dgemv('t',nati,nati,-2.0d+00,gamma3(indi,indi,1),nn2,
     *               x(lshiftg+nn2+indi-1),1,1.0d+00,shift2,1)
          do j = 1, nati
           shift2(j) = shift2(j) + shift2(j+nat)*chamulifg(j)
          end do
        END IF
C       NATFMOB = NATFMO + NBDFG
C       DO LFG = 1, NFG
C         IF (LFG.EQ.IFG) CYCLE
C         NATJ = NATFRG(LFG)
C         INDJ = INDFRG(LFG)
C         CALL DCOPY(NATJ,QINDZ(1,LFG),1,CHATMP,1)
C         IF (DFTB3) THEN
C           IF (IFG.EQ.0) THEN
C             CALL VSUB(X(LZREF+INDJ-1),1,CHATMP,1,CHATMP,1,NATJ)
C             CALL DCOPY(NATJ,CHATMP,1,CHATMP(NATJ+1),1)
C             DO J = 1, NATJ
C               CHATMP(NATJ+J) = CHATMP(NATJ+J)*CHATMP(NATJ+J)
C             END DO
C           ELSE
C             CALL DCOPY(NATJ,POPMAT(1,LFG),1,CHATMP(NATJ+1),1)
C             CALL VSUB(X(LZREF+INDJ-1),1,CHATMP(NATJ+1),1,
C    *                                    CHATMP(NATJ+1),1,NATJ)
C             DO J = 1, NATJ
C               CHATMP(NATJ+J) = CHATMP(J)*CHATMP(NATJ+J)
C             END DO
C           END IF
C         ELSE
C           IF (IFG.EQ.0)
C    *        CALL VSUB(X(LZREF+INDJ-1),1,CHATMP,1,CHATMP,1,NATJ)
C         END IF
C         DO I = 1, NAT !! FRAGMNET TO BE CONSIDERED
C           IAG = IAGLOB(I)
C           IATFG = INDAT(IAG)
C           IATLOC = IAND(IALOC(IAG),65535)
C           INDI = INDFRG(IATFG)
C           NSEQ = INDJ-1 + NATFMOB*(INDI+IATLOC-2)
C           SHIFT(I) = SHIFT(I)
C    *        + DDOT(NATJ,CHATMP,1,X(LGAMMA+NSEQ),1)
C           IF (DFTB3) THEN
C             CHAMULI = CHAMULIFG(I)
C             II = INDI+IATLOC-1
C             IF (.NOT.DFTB2) THEN
C               SHIFT2(I) = SHIFT2(I)
C    *            + TWO*CHAMULI*DDOT(NATJ,CHATMP,1,GAMMA3(2,INDJ,II),2)
C             ELSE IF (IFG.EQ.0) THEN
C               SHIFT2(I) = SHIFT2(I)
C    *            + TWO*CHAMULI*DDOT(NATJ,CHATMP,1,GAMMA3(2,INDJ,II),2)
C    *            + DDOT(NATJ,CHATMP(NATJ+1),1,GAMMA3(1,INDJ,II),2)
C             ELSE
C               CHAMULJ = POPMAT(J,LFG)-X(LZREF+INDJ+J-2)
C               SHIFT2(I) = SHIFT2(I)
C    *            + TWO*CHAMULI*DDOT(NATJ,CHATMP,1,GAMMA3(2,INDJ,II),2)
C    *            + TWO*DDOT(NATJ,CHATMP(NATJ+1),1,GAMMA3(1,INDJ,II),2)
C             END IF
C           END IF
C         END DO
C       END DO
      END IF
C
      if (iand(modgamma,2).eq.2) then
        IF (DFTB3) CALL DAXPY(NAT,ONE_THIRD,SHIFT2,1,SHIFT,1)
        CALL DCOPY(NAT,SHIFT,1,SHIFTOUT,1)
      else
      IF (DFTB3) THEN
        IF (DFTB2) THEN
          CALL DAXPY(NAT,ONE_THIRD,SHIFT2,1,SHIFT,1)
        ELSE
          CALL DCOPY(NAT,SHIFT2,1,SHIFT,1)
          CALL DSCAL(NAT,ONE_THIRD,SHIFT,1)
        END IF
      END IF
      end if
C
      IF (IFG.EQ.0.OR..NOT.DFTB2) THEN
        CALL DCOPY(NAT,SHIFT,1,SHIFTOUT,1)
        RETURN
      END IF
C
C     ----- CALCULATE SHIFT CONTRIBUTIONS AS FAO (WRT1) ----
C     FAO HAS TO BE IN AO BASIS
C
      CALL VCLR(FAO,1,NUM*(NUM+1)/2)
      DO I = 1, NAT
        SHIFTI = SHIFT(I)
        DO K = 1, MAXANG(ISPE(I))**2
          M = IND(I)+K
          DO J = 1, I
            SHIFTJ = SHIFT(J)
            SHIFTIJ = SHIFTI + SHIFTJ
            SHIFTIJ = 0.5D+00*SHIFTIJ
            DO L = 1, MAXANG(ISPE(J))**2
              N = IND(J)+L
              CALL DFTB_CNVSQ(M,N,NSEQ)
              FAO(NSEQ) = SS(NSEQ)*SHIFTIJ
            END DO
          END DO
        END DO
      END DO
c
      RETURN

      END SUBROUTINE DFTB_OFFORBHSS
C
C-----------------------------------------------------------------------
C*MODULE DFTBFO    *DECK DFTB_ESDCPHF
C>
C>    @brief   Fully analytic FMO-DFTB gradient
C>
C>    @details Calculate actual gradient contributions due to Z-vector
C>
C>    @author  Yoshio Nishimoto
C>    - Feb, 2016- Subroutine written
C>    @date    Feb, 2018 - Yoshio Nishimoto
C>    - Added FMO-DFTB/AFO
C>
C>           --- INPUT ---
C>    @param IXFG    Current fragment
C>    @param IXLAYG  Current layer
C>    @param TRI1    Z_{\mu \nu}
C>    @param TRI2    \tilde{Z}_{\mu \nu} (energy-weighted)
C>    @param L1      Number of AOs
C>    @param L2      (L1*L1+L1)/2
C>    @param SCFFRG  SCF type of fragment
C>    @param WRT1    Maybe working array
C>    @param WRK1    Working array
C>    @param WRK2    Working array
C>    @param HOPGRD  HOP gradient should be added or not
C>    @param VEC     MO coefficient
C>    @param EIG     Eigenvalue of each MO
C>    @param Z1      Z-vector in MO
C>    @param VECL    LMO coefficient
C>    @param SL      S*L
C>    @param FMOAFO  FMO/AFO or not
C>           --- IN/OUTPUT ---
C>    @param FMODE   Gradient
C>    @param SHIFTSAV Shift contributions by induced charges
C>
C
      SUBROUTINE DFTB_ESDCPHF(IXFG,IXLAY,TRI1,TRI2,L1,L2,FMODE,
     *  SCFFRG,WRT1,WRK1,WRK2,HOPGRD,VEC,EIG,Z1,VECL,SL,SHIFTSAV,FMOAFO,
     *  GAMLC)
      USE lrcdft, ONLY: LCFLAG, EMU, EMU2, LRFILE
      use mx_limits, only: mxatm,mxsh,mxgtot,mxao
      USE DFTBPB_MOD,ONLY: PERIOD,DFTB_GAMMA_PBC,DFTB_GAMDER_PBC,
     *                     MAXINT,INDBND,PBCBOX,DFTB_FMOSTR,STRESS
C     SUBROUTINE DFTB_ESDCPHF(IXFG,IXLAY,DX,TRI1,TRI2,L1,L2,FMODE,
C    *  LAYFRG,SCFFRG,IDMREC0,IDAM,IDAD,IDAT,WRT1,WRK1,WRK2,WRK3,HOPGRD,
C    *  NAI,NBI,ZVEC1,NQMTFG)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER (ONE=1.0D+00)
      INTEGER,PARAMETER :: MXSPE=10
C
      LOGICAL GOPARR,DSKWRK,MASWRK
      LOGICAL HOPGRD,urohf,EFLDL
      LOGICAL DFTBFL,SCC,SRSCC,DFTB3,DAMPXH,LCDFTB
      LOGICAL MFRZ,FMOAFO
C
      DIMENSION SCFFRG(*)
      DIMENSION TRI1(1), TRI2(1)
      DIMENSION WRT1(*),WRK1(*),WRK2(*)
      DIMENSION FMODE(3,NATFMO)
      DIMENSION VEC(*),EIG(*),Z1(*),VECL(L1,*),SL(*),SHIFTSAV(*),
     *          GAMLC(MAXNAT**2,*)
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
      COMMON /EFLDC / EVEC(3),EFLDL
C     COMMON /EFPFMO/ RUNEFP,REFFMO,REPNUCEFP,IEFPFMO,ISWENUC,EFPEFP
      COMMON /FMCOM / X(1)
      COMMON /FMOAFO/ LINDFRZ,LGMK,LGMKSAV
      COMMON /FMODIM/ MAXBND,MAXKND,MAXCBS,MAXCAO,MAXVEC,MAXL1,MAXNAT,
     *                MAXABD,MAXBAS,MAXBBD,MAXLMO,MAXSLO,MAXABD2,maxrij
      COMMON /FMODTB/ LGAMMAGLOB,LZREF,LISPE,LSHIFTG,LSHIFTCT,
     *                LGAMMAD,MODESD,MODGAMMA,LDCTMP,LGAMMA3,LGAMMA3D,
     *                LNBOND,LCTMUL,LESPDFTB,LESPASC,NDFTBRST,LDFTBRST
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
     *                LPOPDMAT,LIDMPNT,LIDDPNT,LIVMPNT,liactfg,lcrfrg,
     *                lzlmfrgv,lYlmfrgv,lndtfrg,lf_mm,lg_mm,lmaxl30,
     *                libuffg,lindatp,lsmon,lsdim,ledimfed,lexcit3d,
     *                lifgfret
C     COMMON /FMOPMD/ FMOBOX(3),MDWPBC,NIMGCELL,IMGLVL,
C    *                LTRVEC,LFMOGCTR,LFMOCTMP,LINDATMD,LWRKDSAV,LINDXIU
C    *               ,IPBCFST
      COMMON /GRAD  / DE(3,MXATM)
C     COMMON /SYMTRY/ MAPSHL(MXSH,48),MAPCTR(MXATM,48),
C    *                TT(432),INVT(48),NT
      COMMON /FMORUN/ ESPSCF,E0SCF(2),EMP2S,IDAFMO,ICURFG,JCURFG,KCURFG,
     *                ICURLAY,ICURUNT,NAT1E,NCURSH,NGAU,ICURPOP,IFMOSTP,
     *                MONCOR,NEEDR,MODRST,NORBPROJ,NUNESP,ISKIPESP,
     *                IESDPPC,IDOPROP,MP2RUN,ICURIT,IDMFMO,IDDFMO,
     *                IDDCUR,NDDLEFT,IVMFMO,nzmtfmo,ifmobas,itmfmo(2)
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
C     COMMON /INTOPT/ ISCHWZ,IECP,NECP,IEFLD
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON /MFRPAR/ MFRZ,NUMFRZ,NORFRZ,IFRZ(MXAO)
      COMMON /MLTPNT/ LDELEPOT,LDELQSE,LQSETMP,LFZVEC,LMAPCMFD,LQBL2FD
      COMMON /NSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),
     *                CF(MXGTOT),CG(MXGTOT),CH(MXGTOT),CI(MXGTOT),
     *                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),
     *                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
      COMMON /OUTPUT/ NPRINT,ITOL,ICUT,NORMF,NORMP,NOPK
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /PCMDIM/ MXSP,MXTS,MEMPCM1,MEMPCM2,NTS
      COMMON /PCMPAR/ IPCM,NFT26,NFT27,IRPPCM,IEF,IP_F,NFMOPCM,IHET
      COMMON /PCMPNT/ LXYZRE,LSSFE,LLIST,LALPCM,LRINPCM,LINAPCM,LINFPCM,
     *                LINIPCM,LMEPCM,LAXYZCT,LXYZCT2,LNVERT,LQSN,LQSE,
     *                LQFS,LQIND,LISPHE,LVAD,LQOR,LVECMUL,LAIPRJ,LFIPRJ,
     *                LPEL,LPCMCDR,LDAI,LIDDAI,LQSND,LQSED,LEPSPCM
C
      DATA UHF/8HUHF     /,ROHF/8HROHF    /
C
      DOUBLE PRECISION :: AUP(9,9),BUP(9,9),V(3)
      LOGICAL :: SKP
C
C     THIS SUBROUTINE IS CALLED INSIDE IFG LOOP.
C
      IF (IFMOSTP.EQ.6) THEN
        WRITE (IW,*) 'IFMOSTP=6 IS NOT ALLOWED IN FMOESDER'
        CALL ABRT
      END IF
      IF (IFMOSTP.NE.2.AND.IFMOSTP.NE.4.AND.IFMOSTP.NE.9) RETURN
      IESDPPC = 1   ! ALREADY SET IN FMOESDER?
C
C     ORBXCH=MOD(MODORB,2).NE.0
C     .TRUE.  DIVIDE FRAGMENTS
C     .FALSE. DIVIDE SHELLS
C     LWRKESP  = LFMOBUF(2)  ! LWRK2 -> DIJ
C
C     FOR MONOMER SCF THERE IS NO NEED TO ADJUST IREC0 BELOW (AND MPPROP MAY
C     NOT HAVE BEEN READ ANYWAY).
C
      L0       = NQMT
      IFG      = IXFG
      ILAY     = IXLAY
      IXLAY = ILAY !! ftnchk...
      L3       = L1*L1
      INDI     = IXFTCH(X(LINDFRG),IFG)
      NFRZ     = MAX(NUMFRZ,NORFRZ)
      NATFMOB  = NATFMO+NBDFG
C
      IF(IFG.EQ.0) RETURN
      urohf = scffrg(ifg).eq.uhf.or.scffrg(ifg).eq.rohf
C
C     ----- CLEAR DE FIRST -----
C 
      CALL VCLR(DE,1,NAT*3)
      IF (PERIOD) CALL VCLR(STRESS,1,9)
C
C     ----- ALLOCATE MEMORY -----
C
      CALL VALFM(LOADFM)
      LWRK         = LOADFM   + 1
      LDIST        = LWRK     + MAX(L3,NATFMO+NBDFG,L2*2)
      LGAMMA       = LDIST    + NAT*(NAT+1)/2
      LGAMDER      = LGAMMA   + NAT*(NAT+1)/2
      LSHIFT       = LGAMDER  + NAT*(NAT+1)/2
      LSHIFTESP    = LSHIFT   + MAX(NAT,NSHELL)
      LDERSHIFT    = LSHIFTESP+ NAT
      LS           = LDERSHIFT+ NAT
      LFDERAO      = LS       + L2
      LCHAMUL      = LFDERAO  + L2
      LQZ          = LCHAMUL  + NAT
      LAST         = LQZ      + NAT
      IF (DFTB3) THEN
        LGAM3      = LAST
        !! DO NOT USE LGAMMA3, WHICH CONFLICTS WITH COMMON BLOCK
        !! USE LGAM3 AS A WORKING ARRAY IN KFG
        LGAMDER3   = LGAM3    + MAXNAT*MAXNAT*2
        LAST       = LGAMDER3 + NAT*NAT
      ELSE
        LGAM3      = LAST
        LGAMDER3   = LAST
      END IF
      IF (LCDFTB) THEN
        LGAMLC     = LAST
        LGAMDERLC  = LGAMLC   + NAT*(NAT+1)/2
        LAST       = LGAMDERLC+ NAT*(NAT+1)/2
      ELSE
        LGAMLC     = LAST
        LGAMDERLC  = LAST
      END IF
C
      IF (FMOAFO) THEN
        LFAO       = LAST
        LEIGL      = LFAO     + L2
        LUROT      = LEIGL    + L2
        LRROT      = LUROT    + NFRZ*NFRZ
        LZTMP      = LRROT    + NFRZ*NFRZ
        LWRK1      = LZTMP    + L3
        LAST       = LWRK1    + L3
      END IF
C
      IF (FMOAFO.OR.LCDFTB) THEN
        LD         = LAST
        LAST       = LD       + L2
      END IF
C
      NEED         = LAST - LOADFM - 1
      CALL GETFM(NEED)
C
C     ----- CALCULATE OVERLAP MATRIX -----
C
      IF (FMOAFO) THEN
        CALL VCLR(X(LFAO),1,L2)
        CALL VCLR(X(LS),1,L2)
        CALL DFTB_HAM_OVER(X(LFAO),X(LS),.TRUE.,.TRUE.,.FALSE.)
        CALL DCOPY(NAT,SHIFTSAV(INDI),1,X(LSHIFT),1)
        CALL DFTB_SHIFT_FOCK(X(LSHIFT),X(LFAO),X(LS),NAT,L2,NSPE,
     *       ISPE,IND,MAXANG)
        CALL DAWRIT(IDAF,IODA,X(LFAO),L2,317,0)
      ELSE
        CALL VCLR(X(LS),1,L2)
        CALL DFTB_HAM_OVER(X(LS),X(LS),.FALSE.,.TRUE.,.FALSE.)
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
C     ----- CALCULATE GAMMA AND ITS DERIVATIVE OF DFTB2 AND DFTB3 -----
C
      IF (PERIOD) THEN
        IF (SRSCC) THEN
          NNN = 0
          DO I = 1, NSPE
            NNN = NNN + MAXANG(I)
          END DO
        ELSE
          NNN = NSPE
        END IF
        CALL VALFM(LOADFM)
        LTMPGAM  = LOADFM + 1
        LAST     = LTMPGAM + 11*NNN*NNN
        NEEDGAM  = LAST  - LOADFM - 1
        CALL GETFM(NEEDGAM)
        CALL DFTB_GAMMA_PBC(0,0,0,NAT,NSPE,NSHELL,NNN,ISPE,MAXANG,
     *                      C,SPE,HUBBL,DAMPXHE,HUBDER,
     *                      DAMPXH,SRSCC,DFTB3,X(LGAMMA),
     *                      X(LGAM3),X(LTMPGAM),X(LAST),X(LAST),
     *                      X(LAST),X(LAST),X(LAST),.FALSE.,
     *                      DUM,X(LAST))
        CALL RETFM(NEEDGAM)
       !DO I = 1, NAT
       !  ISP = ISPE(I)
       !  DO J = 1, I
       !    JSP = ISPE(J)
       !    NSEQIJ = ISP-1+NNN*(JSP-1)
       !    NSEQJI = JSP-1+NNN*(ISP-1)
       !    IF (I.EQ.J) THEN
       !      DIST = 0.0D+00
       !    ELSE
       !      DIST = SQRT((C(1,I)-C(1,J))**2
     * !                 +(C(2,I)-C(2,J))**2
     * !                 +(C(3,I)-C(3,J))**2)
       !    END IF
       !    CALL DFTB_GAMDER_PBC(I,J,C(1,I),C(1,J),GAMMAV,
     * !                         GAMI,GAMJ,X(LGAMMAD+10*NSEQIJ),
     * !                         X(LGAMMA3D+11*NSEQIJ),
     * !                         X(LGAMMA3D+11*NSEQJI),
     * !                         DAMPXH,ISP,JSP,CHAMUL1,CHAMUL2,
     * !                         CHAMULSQ)
       !    GAMMAV = GAMMAV*DIST
       !    CALL DFTB_CNVSQ(I,J,NSEQ)
       !    X(LGAMMA+NSEQ-1) = GAMMAV
       !    IF (DFTB3) THEN
       !      GAMI = GAMI*DIST
       !      GAMJ = GAMJ*DIST
       !      X(LGAMDER3+I-1+NAT*(J-1)) = GAMI
       !      X(LGAMDER3+J-1+NAT*(I-1)) = GAMJ
       !    END IF
       !  END DO
       !END DO
      ELSE
        CALL DFTB_GAMMA(DAMPXHE,HUBBL,X(LGAMMA),X(LDIST),NAT,NSPE,
     *    ISPE,MAXANG,SPE,DAMPXH,SRSCC)
        NN = NSPE
        CALL DFTB_GAMDER(NAT,NSPE,HUBBL,ISPE,X(LGAMDER),DAMPXH,
     *    DAMPXHE,SPE,X(LDIST),DFTB3,X(LGAMDER3),HUBDER,SRSCC,MAXANG,
     *    NAT,NN)
        IF (DFTB3) THEN
          !! DFTB_GAMDER GIVES GAMDER3 VALUES SCALED BY 1/3
          CALL DSCAL(NAT*NAT,3.0D+00,X(LGAMDER3),1)
          !! NOW, GAMMA3 VALUES ARE NOT SCALED
          CALL DFTB_GAMMA3(NAT,HUBBL,X(LGAM3),X(LDIST))
        END IF
      END IF
      !! GET MULLIKEN CHARGES OF I-TH FRAGMENT AND DEFINE LPOPMATI
      CALL DFTB_ZVEC_MUL(NAT,LPOPMAT,MAXNAT,IFG,X(LCHAMUL),LPOPMATI)
C
      IF (FMOAFO.OR.LCDFTB) THEN
        CALL VCLR(WRK1,1,L1)
        DO I = 1, NA
          WRK1(I)=2.0D+00
        END DO
        CALL DFTB_DMTOCC(X(LD),VEC,WRK1,WRK2,L1,L1,1.0D-09)
      END IF
C
      IF (IAND(MODGAMMA,2).NE.0) THEN
        NN = 2
        NN2 = NATFMO+NBDFG
      ELSE
        NN = 1
        NN2 = NATFMO
      END IF
      IF (FMOAFO) THEN
        NNUM = (L1-NA)*NA
        NORBPRJ = NFRZ - NUMFRZ
        CALL DAREAD(IDAF,IODA,X(LUROT),NFRZ*NFRZ,569,0)
        CALL DFTB_RROT(L1,NFRZ,NUMFRZ,X(LUROT),X(LRROT),X(LS),
     *                 VECL,VEC,WRK1,WRK2)
        IF (NORBPRJ.GT.0) THEN
C         CALL DCOPY(L1*NORBPRJ,VECL(1,NUMFRZ+1),1,X(LWRK1),1)
C         CALL DCOPY(L1*(L1-NUMFRZ),VECL(1,NORFRZ+1),1,
C    *               VECL(1,NUMFRZ+1),1)
C         CALL DCOPY(L1*NORBPRJ,X(LWRK1),1,VECL(1,L1-NORBPRJ+1),1)
          CALL DCOPY(L1*NORBPRJ,VECL(1,NUMFRZ+1),1,
     *               VECL(1,L1-NORBPRJ+1),1)
        END IF
        CALL DFTB_ZVEC_ZB(1,L1,NNUM,NFRZ,Z1,TRI1,TRI2,X(LCHAMUL),X(LQZ),
     *                    X(LSHIFT),X(LAST),EIG,X(LS),X(LD),X(LFAO),
     *                    X(LEIGL),VEC,VECL,X(LUROT),X(LRROT),X(LGMK),
     *                    X(LZTMP),SL,X(LGAMMA),X(LGAM3),X(LGAMMA3),
     *                    WRT1,WRK1,WRK2,X(LWRK),X(LWRK1),.TRUE.,DFTB3,
     *                    MODGAMMA,NATFMO,NBDFG,X(LIAGLOB),X(LNATFRG),
     *                    X(LINDFRG),X(LIATFRG),X(LPOPMATI),X(LISPE),
     *                    X(LFMOC),X(LINDAT),X(LIALOC),
     *                    X(LESPDFTB+2*NATFMOB+INDI-1))
        IF (NORBPRJ.GT.0) THEN
          CALL DCOPY(L1*NORBPRJ,VECL(1,L1-NORBPRJ+1),1,
     *               VECL(1,NUMFRZ+1),1)
        END IF
        CALL DSCAL(L2,-ONE,TRI2,1)
        JCURFG = -1
        CALL DFTB_PUTGMK(L1,1,X(LGMK),X(LGMKSAV),X(LINDFRZ),X(LIALMO),
     *                   X(LIAGLOB),X(LAST),X(LAST),X(LAST),.FALSE.,
     *                   X(LAST),X(LAST),X(LAST))
        JCURFG = 0
      ELSE
C
C       ----- CALCULATE Q_A^Z CHARGES AS X(LQZ) -----
C
        CALL DFTB_MULLIKEN(L1,L2,TRI1,X(LS),X(LWRK))
        !! X(LSHIFT) IS USED AS A WORKING ARRAY
        CALL DFTB_MULSA(X(LWRK),X(LSHIFT),X(LQZ),L1,NSHELL,NAT)
C
C       ----- CONSTRUCT INTRA-FRAGMENT SHIFTZ ARRAY -----
C
        CALL VCLR(X(LSHIFT),1,NAT)
        IF (DFTB3) THEN
          IF (IAND(MODGAMMA,4096).EQ.0) THEN
            IF (IAND(MODGAMMA,8).EQ.0) THEN
              CALL DFTB_OFFORBHSS(0,X(LIAGLOB),X(LNATFRG),X(LINDFRG),
     *          X(LIATFRG),X(LPOPMATI),X(LS),WRT1,X(LISPE),X(LFMOC),
     *          X(LGAMMA3),X(LSHIFT),X(LQZ),X(LPOPMATI),.FALSE.,
     *          X(LINDAT),X(LIALOC),NN,NN2)
              DO I = 1, NAT 
                DO J = 1, NAT
                  GAMI = X(LGAM3+I-1+(J-1)*NAT)
                  GAMJ = X(LGAM3+J-1+(I-1)*NAT)
                  X(LSHIFT+I-1) = X(LSHIFT+I-1)
     *              + 2.0D+00/3.0D+00*X(LQZ+J-1)*X(LCHAMUL+I-1)*GAMI
     *              + 2.0D+00/3.0D+00*X(LQZ+J-1)*X(LCHAMUL+J-1)*GAMJ
                END DO
              END DO
            ELSE
              DO I = 1, NAT 
                DO J = 1, NAT
                  GAMI = X(LGAM3+I-1+(J-1)*NAT)
                  GAMJ = X(LGAM3+J-1+(I-1)*NAT)
                  X(LSHIFT+I-1) = X(LSHIFT+I-1)
     *              + 2.0D+00/3.0D+00*X(LQZ+J-1)*X(LCHAMUL+I-1)*GAMI
     *              + 2.0D+00/3.0D+00*X(LQZ+J-1)*X(LCHAMUL+J-1)*GAMJ
     *              + 2.0D+00/3.0D+00*X(LQZ+I-1)*X(LCHAMUL+J-1)*GAMI
                END DO
                X(LSHIFT+I-1) = X(LSHIFT+I-1)
     *            + 2.0D+00*X(LQZ+I-1)*X(LESPDFTB+2*NATFMOB+INDI+I-2)
              END DO
            END IF
          ELSE
            CALL DFTB_TD_3RD_SHIFT(X(LSHIFT),X(LSHIFT),X(LQZ),X(LQZ),
     *                             X(LCHAMUL),X(LGAM3),NAT,.TRUE.)
          END IF
        END IF
        CALL DFTB_2ND_SHIFT(X(LSHIFT),X(LQZ),X(LGAMMA),NAT)
      END IF
      !! This contribution comes from intra-fragment
      IF (MASWRK) CALL DAXPY(NAT,ONE,X(LSHIFT),1,X(LSHIFTCT+INDI-1),1)
C
C     ----- CONSTRUCT SHIFT ARRAY -----
C     CALLING WITH 0 FOR THE FIRST ARGUMENT, SHIFT CONTRIBUTIONS
C     FROM ALL ATOMS INCLUDING IFG WILL BE CALCULATED.  SHIFT IS
C     RETURNED AS X(LSHIFTESP).
C
      !! X(LSHIFTESP) SHOULD BE IDENTICAL TO SHIFT+SHIFT2 IN DFTBX.
      IF (IAND(MODGAMMA,10).NE.0) THEN
        CALL DCOPY(NAT,X(LESPDFTB+INDI-1),1,X(LSHIFTESP),1)
        IF (DFTB3) THEN
          DO I = 1, NAT
            X(LSHIFTESP+I-1) = X(LSHIFTESP+I-1)
     *      +X(LESPDFTB+NATFMOB+INDI+I-2)*X(LCHAMUL+I-1)*2.0D+00/3.0D+00
          END DO
        END IF
      ELSE
        CALL DFTB_OFFORBHSS(0,X(LIAGLOB),X(LNATFRG),X(LINDFRG),
     *       X(LIATFRG),X(LPOPMATI),X(LS),WRT1,X(LISPE),X(LFMOC),
     *       X(LGAMMA3),X(LSHIFTESP),X(LCHAMUL),X(LPOPMATI),.TRUE.,
     *       X(LINDAT),X(LIALOC),NN,NN2)
      END IF
      IF (IAND(MODGAMMA,4096).NE.0) THEN
        CALL DFTB_2ND_SHIFT(X(LSHIFTESP),X(LCHAMUL),X(LGAMMA),NAT)
        IF (DFTB3) THEN
          CALL DFTB_3RD_SHIFT(X(LSHIFTESP),X(LCHAMUL),X(LGAM3),NAT)
        END IF
      END IF
      IF (EFLDL) CALL DFTB_EFIELD_SHIFT(X(LSHIFTESP),EVEC,C,NAT)
C     IF (MASWRK) CALL DCOPY(NAT,X(LQZ),1,X(LESPDFTB+INDI-1),1)
C
C     ---- PCM<1> CONTRIBUTIONS (V(Z)^a TERM) -----
C     ADD ONLY Z*dS/da CONTRIBUTIONS. dr_Al/da, D*dS/da, and q*dC/da*q
C     TERMS ARE COMPUTED IN THE CCMBEM AND DFTB_SDER SUBROUTINE.
C
      IF (NFMOPCM.NE.0) THEN
        CALL VCLR(X(LWRK),1,NAT)
        CALL DFTB_QSHIFT(NAT,NTS,C,X(LAXYZCT),X(LAXYZCT+MXTS),
     *                   X(LAXYZCT+MXTS*2),X(LQSE),X(LWRK))
        CALL DAXPY(NAT,ONE,X(LWRK),1,X(LSHIFTESP),1)
        !! for geometry-dependence
        IF (MASWRK.AND.FMOAFO) THEN
          CALL VCLR(X(LWRK),1,NAT)
          CALL DFTB_QSHIFT(NAT,NTS,C,X(LAXYZCT),X(LAXYZCT+MXTS),
     *                    X(LAXYZCT+MXTS*2),X(LQSETMP),X(LWRK))
          CALL DAXPY(NAT,ONE,X(LWRK),1,X(LSHIFTG+INDI-1),1)
        END IF
      END IF
C
C     ----- ADD HOP TERMS -----
C     THIS (SPECIFICALLY, FMOHOP) HAS TO BE CALLED BEFORE DFTB_MAKEREFD
C     FOR LC-DFTB, BECAUSE ROTCHOP IS NOT DEFINED OTHERWISE.
C
      IF (HOPGRD) THEN
        I2D = 0
        I3D = 0
        IF (NBODY.EQ.2) I2D = 1 
        IF (NBODY.EQ.3) I3D = 1 
CZ      WRITE(6,*) 'HOPDER IN SCZV'
        CALL HOPDER(1,0,I2D,I3D,L1,L2,TRI1,WRK1,WRK2,DUM)
        if(urohf)  then
          CALL HOPDER(1,0,I2D,I3D,L1,L2,TRI1(L2+1),WRK1,WRK2,DUM)
        end if
      END IF
C
C     ----- EXCHANGE CONTRIBUTIONS FOR LC-DFTB -----
C
      IF (LCDFTB) THEN
        CALL VCLR(WRK1,1,L2)
        IF (ICUT.GE.0) THEN
          !! GAMLC(1,2): SMAX
          !! GAMLC(1,3): INTLST
          CALL VALFM(LOADFM)
          LWRK1  = LOADFM + 1
          LWRK2  = LWRK1  + L1*L1
          LSUMDDA= LWRK2  + L2
          LWRK3  = LSUMDDA+ NAT
          LAST   = LWRK3  + NAT*NAT
          NEED2  = LAST - LOADFM - 1
          CALL GETFM(NEED2)
          CALL CPYTSQ(X(LS),WRK2,L1,1)
          CALL DFTB_LCSHIFT(1,NAT,L1,GAMLC,WRK2,X(LAST),
     *                      X(LAST),TRI1,X(LWRK1),X(LAST),WRK1,
     *                      X(LWRK2),GAMLC(1,2),X(LSUMDDA),
     *                      X(LWRK3),ELC,.FALSE.,ICUT,GAMLC(1,3))
          CALL RETFM(NEED2)
        ELSE
          CALL DFTB_LCSHIFT2(1,NAT,1,L1,ISPE,MAXANG,IND,TRI1,DUM,
     *                       X(LS),GAMLC,WRK1,ELC)
C         IF (UROHF) THEN
C           CALL DCOPY(L2,WRK1(1+L2),1,X(LWRK),1)
C           CALL DCOPY(L2,X(LWRK),1,WRK1(1+L3),1)
C         END IF
        END IF
      ! !! GAMLC(1,2): SMAX
      ! !! GAMLC(1,3): INTLST
      ! CALL VALFM(LOADFM)
      ! LWRK1  = LOADFM + 1
      ! LWRK2  = LWRK1  + L1*L1
      ! LSUMDDA= LWRK2  + L2
      ! LWRK3  = LSUMDDA+ NAT
      ! LAST   = LWRK3  + NAT*NAT
      ! NEED2  = LAST - LOADFM - 1
      ! CALL GETFM(NEED2)
C
      ! !! Exchange Fock due to Z
      ! CALL CPYTSQ(X(LS),WRK2,L1,1)
      ! CALL VCLR(WRK1,1,L2)
      ! write (*,*) "6.1"
      ! CALL DFTB_LCSHIFT2(1,NAT,1,L1,ISPE,MAXANG,IND,TRI1,DUM,
     *!                    X(LS),GAMLC,WRK1,ELC)
C     ! CALL DFTB_LCSHIFT(1,NAT,L1,GAMLC,WRK2,X(LAST),
C    *!                   X(LAST),TRI1,X(LWRK1),X(LAST),WRK1,
C    *!                   X(LWRK2),GAMLC(1,2),X(LSUMDDA),
C    *!                   X(LWRK3),ELC,.FALSE.,ICUT,GAMLC(1,3))
      ! CALL RETFM(NEED2)
      ! write (*,*) "6.2"
C
        CALL CPYTSQ(X(LD),WRK2,L1,1)
        !! D*F^Z*D
        CALL TFTRI(X(LWRK),WRK1,WRK2,X(LWRK+L2),L0,L1,L1)
        !! Add to W
        CALL DAXPY(L2,-0.5D+00,X(LWRK),1,TRI2,1)
C
        !! Construct \Delta D
        CALL VCLR(X(LWRK),1,L2)
        CALL DFTB_MAKEREFD(NAT,X(LWRK))
        CALL DCOPY(L2,X(LD),1,WRK1,1)
        CALL DAXPY(L2,-1.0D+00,X(LWRK),1,WRK1,1)
        CALL CPYTSQ(WRK1,WRK2,L1,1) !! \Delta D
        CALL CPYTSQ(TRI1,WRK1,L1,1) !! Z
C
        CALL DFTB_LCGAMMADER(EMU,HUBB,NAT,ISPE,X(LDIST),
     *                       X(LGAMLC),X(LGAMDERLC))
C
        !! Prepare for S derivative contributions in (ab|cd)^a
        !! Z and \Delta D
        CALL DFTB_LCGRAD1(1,NAT,L1,ISPE,MAXANG,IND,WRK1,WRK2,
     *                    X(LS),X(LGAMLC),X(LWRK))
        CALL DAXPY(L2,+2.0D+00,X(LWRK),1,TRI2,1)
C
        !! gamma derivative contributions in (ab|cd)^a
        IF (MASWRK) THEN
          CALL CPYTSQ(X(LS),X(LWRK),L1,1)
          CALL DFTB_LCGRAD2(1,2.0D+00,NAT,L1,ISPE,MAXANG,
     *                      IND,WRK1,WRK2,X(LWRK),X(LWRK),X(LDIST),
     *                      X(LGAMDERLC),C,DE,X(LAST))
C         CALL DFTB_LCGRAD2(1,NAT,L1,ISPE,MAXANG,IND,WRK1,WRK2,
C    *                      X(LWRK),X(LDIST),X(LGAMDERLC),C,DE)
C         CALL DFTB_LCGRAD2(1,NAT,L1,ISPE,MAXANG,IND,WRK2,WRK1,
C    *                      X(LWRK),X(LDIST),X(LGAMDERLC),C,DE)
        END IF
      END IF
C
      CALL VCLR(WRK1,1,L2)
      CALL VCLR(WRK2,1,L2)
      CALL VCLR(AUP,1,81)
      CALL VCLR(BUP,1,81)
c     DO NI = 1, 0 ! NAT
c       DO KI = 1, 3
C
C         ----- CALCULATE HDERAO AND SDERAO TERM -----
C         HDERAO: H_{\mu \nu}^{0,a}
C         SDERAO: S_{\mu \nu}^{a}
C
c          IF (PERIOD) THEN
c            CALL VCLR(WRK1,1,L2)
C            CALL VCLR(WRK2,1,L2)
C            NNN = 4210
C            CALL VCLR(AUP,1,81)
C            CALL VCLR(BUP,1,81)
C            NINTPAIR = INDBND(1,1,NI)
c            DO JJ = 2, NINTPAIR !! skip JJ=1, because the same atom
c              NJ = INDBND(1,JJ+1,NI)
C              IF (NI.EQ.NJ) CYCLE
C              NNN = INDBND(2,JJ+1,NI)
C              IXTRV =  NNN/400-10
C              IYTRV = (NNN-(IXTRV+10)*400)/20-10
C              IZTRV = (NNN-(IXTRV+10)*400-(IYTRV+10)*20)-10
c              V(1) = C(1,NJ) - C(1,NI) + PBCBOX(1,1)*IXTRV
c     *             + PBCBOX(2,1)*IYTRV + PBCBOX(3,1)*IZTRV
C              V(2) = C(2,NJ) - C(2,NI) + PBCBOX(1,2)*IXTRV
C     *             + PBCBOX(2,2)*IYTRV + PBCBOX(3,2)*IZTRV
C              V(3) = C(3,NJ) - C(3,NI) + PBCBOX(1,3)*IXTRV
C     *             + PBCBOX(2,3)*IYTRV + PBCBOX(3,3)*IZTRV
C              CALL DFTB_SKODE_HELPA(1,NI,NJ,KI,0,0,0,V,AUP,BUP,
c     *                              .TRUE.,.TRUE.,SKP)
c              IF (SKP) CYCLE
C              DO M = 1, IND(NJ+1)-IND(NJ)
C                MU = IND(NJ)+M
C                DO N = 1, IND(NI+1)-IND(NI)
C                  NU = IND(NI)+N
C                  CALL DFTB_CNVSQ(MU,NU,NSEQ)
c                  WRK1(NSEQ) = WRK1(NSEQ) + AUP(N,M) !! HDERAO
c                  WRK2(NSEQ) = WRK2(NSEQ) + BUP(N,M) !! SDERAO
C                END DO
C              END DO
C            END DO
C          ELSE
C            DO NJ = 1, NAT
c              IF (NI.EQ.NJ) CYCLE
c              V(1) = C(1,NJ)-C(1,NI)
C              V(2) = C(2,NJ)-C(2,NI)
C              V(3) = C(3,NJ)-C(3,NI)
C              CALL DFTB_SKODE_HELPA(1,NI,NJ,KI,0,0,0,V,AUP,BUP,.TRUE.,
C     *             .TRUE.,SKP)
C              IF (SKP) CYCLE
c              DO M = 1, IND(NJ+1) - IND(NJ)
c                MU = IND(NJ)+M
C                DO N = 1, IND(NI+1) - IND(NI)
C                  NU = IND(NI)+N
C                  CALL DFTB_CNVSQ(MU,NU,NSEQ)
C                  WRK1(NSEQ) = AUP(N,M) !! HDERAO
C                  WRK2(NSEQ) = BUP(N,M) !! SDERAO
c                END DO
c              END DO
c            END DO
c          END IF
C
C         ----- CALCULATE F_{\mu \nu}^a TERMS FOR IFG -----
C         ALL DS+WS TERMS ARE CALCULATED IN DFTB_ZVEC_KGRAD SUBROUTINE
C
C          CALL VCLR(X(LFDERAO),1,L2)
c          !! S_{\mu \nu}/da * \Omega_{AB}
c          CALL SHIFT_TRIANGLE_DER(NAT,NI,L2,IND,X(LSHIFTESP),WRK2,
c     *         X(LFDERAO))
c          !! H_{\mu \nu}^0/da
C          CALL DAXPY(L2,ONE,WRK1,1,X(LFDERAO),1)
C          !! S_{\mu \nu} * \Omega_{AB}/da
C      call vclr(x(ldershift),1,nat)
C         CALL DFTB_CALC_DERSHIFT(NAT,NI,KI,C,X(LDIST),X(LGAMDER),
C    *      X(LGAMDER3),X(LCHAMUL),X(LDERSHIFT),DFTB3)
C         CALL DFTB_SHIFT_FOCK(X(LDERSHIFT),X(LFDERAO),X(LS),NAT,L2,
C    *      NSPE,ISPE,IND,MAXANG)
C
C         ----- CALCULATE Z_{\mu \nu} F_{\mu \nu}^a TERM -----
C
        ! IF (MASWRK) DE(KI,NI) = DE(KI,NI) + TRACEP(TRI1,X(LFDERAO),L1)
C
C         ----- CALCULATE \tilde{Z}_{\mu \nu} S_{\mu \nu}^a TERM -----
C
c        ! IF (MASWRK) DE(KI,NI) = DE(KI,NI) + TRACEP(TRI2,WRK2,L1)
c        END DO
c        CALL SHIFT_TRIANGLE_DER_CLEAR(NAT,NI,L2,IND,WRK1)
c        CALL SHIFT_TRIANGLE_DER_CLEAR(NAT,NI,L2,IND,WRK2)
c      END DO
C
      NNN = 4210
      CALL VCLR(AUP,1,81)
      CALL VCLR(BUP,1,81)
      IPCOUNT = ME - 1
      DO NI = 1, NAT
        IF(GOPARR) THEN
          IPCOUNT = IPCOUNT + 1
          IF(MOD(IPCOUNT,NPROC).NE.0) CYCLE
        END IF
        SHIFTI = X(LSHIFTESP+NI-1)
        IF (MAXINT.EQ.0) THEN
          NINTPAIR = NAT
        ELSE
          NINTPAIR = INDBND(1,1,NI)
        END IF
        DO JJ = 1, NINTPAIR
          IF (MAXINT.EQ.0) THEN
            NJ = JJ
          ELSE
            NJ = INDBND(1,JJ+1,NI)
          END IF
          IF (PERIOD) THEN
            NNN = INDBND(2,JJ+1,NI)
            IXTRV =  NNN/400-10
            IYTRV = (NNN-(IXTRV+10)*400)/20-10
            IZTRV = (NNN-(IXTRV+10)*400-(IYTRV+10)*20)-10
          END IF
          IF (NI.EQ.NJ .AND. NNN.EQ.4210) CYCLE
          SHIFTJ = X(LSHIFTESP+NJ-1)
          SHIFTIJ = (SHIFTI+SHIFTJ)*0.5D+00
          V(1) = C(1,NJ) - C(1,NI)
          V(2) = C(2,NJ) - C(2,NI)
          V(3) = C(3,NJ) - C(3,NI)
          IF (PERIOD) THEN
            V(1) = V(1) + PBCBOX(1,1)*IXTRV
     *           + PBCBOX(2,1)*IYTRV + PBCBOX(3,1)*IZTRV
            V(2) = V(2) + PBCBOX(1,2)*IXTRV
     *           + PBCBOX(2,2)*IYTRV + PBCBOX(3,2)*IZTRV
            V(3) = V(3) + PBCBOX(1,3)*IXTRV
     *           + PBCBOX(2,3)*IYTRV + PBCBOX(3,3)*IZTRV
          END IF
          DO KI = 1, 3
            CALL DFTB_SKODE_HELPA(1,NI,NJ,KI,0,0,0,V,AUP,BUP,
     *                            .TRUE.,.TRUE.,SKP)
            IF (SKP) CYCLE
            STRTMP = 0.0D+00
            DO M = 1, IND(NJ+1)-IND(NJ)
              MU = IND(NJ)+M
              DO N = 1, IND(NI+1)-IND(NI)
                NU = IND(NI)+N
                CALL DFTB_CNVSQ(MU,NU,NSEQ)
                DHDR = 2.0D+00*AUP(N,M)
                DSDR = 2.0D+00*BUP(N,M)
                TMP = DHDR*TRI1(NSEQ)
     *              + DSDR*(SHIFTIJ*TRI1(NSEQ)+TRI2(NSEQ))
C               DE(KI,NI) = DE(KI,NI) + TMP
                STRTMP = STRTMP + TMP
              END DO
            END DO
            DE(KI,NI) = DE(KI,NI) + STRTMP
            IF (PERIOD) THEN
              DO K2 = 1, 3
                STRESS(KI,K2) = STRESS(KI,K2) - STRTMP*V(K2)
              END DO
            END IF
          END DO
        END DO
      END DO
C
      IF (GOPARR) CALL DDI_GSUMF(2418,DE,NAT*3)
      IF (GOPARR.AND.PERIOD) CALL DDI_GSUMF(2418,STRESS,9)
C
C     ----- LOOP FOR K-TH FRAGMENT -----
C
      IF (IAND(MODGAMMA,4096).EQ.0) THEN
        CALL DFTB_ZVEC_KFG(IFG,X(LNATFRG),X(LIATFRG),
     *    X(LINDFRG),FMODE,X(LISPE),
     *    X(LFMOC),X(LS),TRI1,X(LPOPMATI),NSHELL,
     *    X(LSHIFTCT))
C       IF (MASWRK) CALL DAXPY(NAT,ONE,X(LSHIFT),1,X(LSHIFTCT+INDI-1),1)
      END IF
C
C     THE COLLECTED DE IS NOW COPIED TO FMODE
C 
      IF (MASWRK) THEN
        IF     (IFMOSTP.EQ.4) THEN
          CALL FMODEG(1,FMODE,X(LFMOPG),X(LIAGLOB))
        ELSEIF (IFMOSTP.EQ.9) THEN
          CALL FMODEG(1,FMODE,X(LFMOPG),X(LIAGLOB))
        ENDIF
        if(period) call dftb_fmostr(1,0,0,1)
      ENDIF
C
      CALL RETFM(NEED)
      IESDPPC = 0
C
      RETURN
C
      END SUBROUTINE DFTB_ESDCPHF
C
C-----------------------------------------------------------------------
C*MODULE DFTBFO    *DECK DFTB_ZVEC_KFG
C>
C>    @brief   Fully analytic FMO-DFTB gradient
C>
C>    @details Calculate gradient contributions in KFG loop
C>
C>    @author  Yoshio Nishimoto
C>    - Feb, 2016- Subroutine written
C>    @date    Feb, 2018 - Yoshio Nishimoto
C>    - Added FMO-DFTB/AFO
C>
C>           --- INPUT ---
C>    @param IIFG    Current fragment
C>    @param NATFRG  Number of atoms in each fragment
C>    @param IATFRG  Global atom index with NBDFG?
C>    @param INDFRG  First atom of each fragment in terms of NATFMOB?
C>    @param ISPEFMO ISPE for all atoms (NATFMO)
C>    @param FMOC    Coordinates for all NATFMO
C>    @param S       Overlap matrix in AO
C>    @param Z       Z_{\mu \nu}
C>    @param POPMAT  Mulliken population
C>    @param NSHELL  Number of shells
C>    @param GAMMA3  DFTB3 gamma
C>           --- IN/OUTPUT ---
C>    @param FMODE   Gradient
C>    @param SHIFTZ  Shift due to induced charges by Z
C>
C
      SUBROUTINE DFTB_ZVEC_KFG(IIFG,NATFRG,IATFRG,
     *  INDFRG,FMODE,ISPEFMO,FMOC,
     *  S,Z,POPMAT,NSHELL,SHIFTZ)
      use mx_limits, only: mxatm
      USE DFTBPB_MOD,ONLY: PERIOD,DFTB_GAMFMO_PBC,DFTB_GAMDER_PBC,STRFMO
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      INTEGER,PARAMETER :: MXSPE=10
      DOUBLE PRECISION, PARAMETER :: ONE=1.0D+00,TWO=2.0D+00,
     *  ZERO=0.0D+00
      LOGICAL DFTBFL,SCC,SRSCC,DFTB3,DAMPXH0,LCDFTB,dotrunces,dotruncesd
      LOGICAL GOPARR,DSKWRK,MASWRK
C
      COMMON /DFTB  / DFTBFL,SCC,SRSCC,DFTB3,DAMPXH0,LCDFTB
      COMMON /DFTBPR/ ETEMP,DFTBDP(MXSPE*14),DAMPXHE,HUBDER(MXSPE),
     *                ZREF(MXATM),SPNCST(6,MXSPE),SPE(MXATM),NSPE,
     *                MAXANG(MXATM),ISPE(MXATM),IND(MXATM+1),IDFTBD,
     *                PARAMDIR
      COMMON /DFTBSK/ SKDIM(MXSPE,MXSPE),SKSPIN(MXSPE),QREFL(3,MXSPE),
     *                HUBBL(3,MXSPE),QREF(MXSPE),HUBB(MXSPE),SKCUT2,
     *                NEEDSK,LSKHTAB(MXSPE,MXSPE),LSKSTAB(MXSPE,MXSPE),
     *                LSKGRID(MXSPE,MXSPE),LSKSELF(MXSPE)
      COMMON /FMCOM / X(1)
      COMMON /FMOINF/ NFG,NLAYER,NATFMO,NBDFG,NAOTYP,NBODY,NSEGM
      COMMON /FMODTB/ LGAMMAGLOB,LZREF,LISPE,LSHIFTG,LSHIFTCT,
     *                LGAMMAD,MODESD,MODGAMMA,LDCTMP,LGAMMA3,LGAMMA3D,
     *                LNBOND,LCTMUL,LESPDFTB,LESPASC,NDFTBRST,LDFTBRST
      common /fmodim/ maxbnd,maxknd,maxcbs,maxcao,maxvec,maxl1,maxnat,
     *                maxabd,maxbas,maxbbd,maxlmo,maxslo,maxabd2,maxrij
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      common /cenrad/ rcutes,lcradfg,lichfgg,dotrunces,dotruncesd
C
      INTEGER INDFRG(*),IATFRG(*),NATFRG(*),
     *  ISPEFMO(*)
      DOUBLE PRECISION FMODE(3,NATFMO),FMOC(3,NATFMO),
     *  TMP(12,NSPE,NSPE),TMP3(10,NSPE,NSPE),
     *  TMPD3(11,NSPE,NSPE),VEC(3),TMPD(10,NSPE,NSPE)
      DOUBLE PRECISION S(*),Z(*),POPMAT(MAXNAT,*),SHIFTZ(*)
      LOGICAL DAMPXH
      !! Gamma derivative vector and stress tensor for PBC
      DIMENSION GDV2(3),GDV3I(3),GDV3J(3),
     *          GDST2(3,3),GDST3I(3,3),GDST3J(3,3)
C
      IF (.NOT.MASWRK) RETURN
C
      IF (IAND(MODGAMMA,3).EQ.0) THEN
        CALL DCOPY(12*NSPE*NSPE,X(LGAMMAGLOB),1,TMP,1)
        IF (DFTB3) CALL DCOPY(10*NSPE*NSPE,X(LGAMMA3),1,TMP3,1)
      END IF
C     IF (IAND(MODGAMMA,12).EQ.0) THEN
        CALL DCOPY(10*NSPE*NSPE,X(LGAMMAD),1,TMPD,1)
        IF (DFTB3) CALL DCOPY(11*NSPE*NSPE,X(LGAMMA3D),1,TMPD3,1)
C     END IF
C
C
      IFG = IIFG
      L1 = MAXL1
      L2 = (L1*L1+L1)/2
      L3 = L1*L1
C
C     ----- ALLOCATE MEMORY -----
C     USE THE SIZE OF THE BIGGEST FRAGMENT
C
      CALL VALFM(LOADFM)
      CALL GOTFM(NGOTMX)
      LWRK    = LOADFM  + 1
      LWRK1   = LWRK    + L1
      LWRK2   = LWRK1   + L3
      LWRK3   = LWRK2   + L3
      LOCC    = LWRK3   + L3
      LD      = LOCC    + L1
      LS      = LD      + L2
      LDS     = LS      + L2
      LQZ     = LDS     + L3 !! D SQUARE
      LQTMP   = LQZ     + MAXNAT
      LGAMMA  = LQTMP   + MAXNAT
      LSDERAO = LGAMMA  + MAXNAT*MAXNAT
      LAST    = LSDERAO + L2
      NEED    = LAST - LOADFM - 1
C
C     ----- GET MEMORY -----
C
      CALL GETFM(NEED)
C
C     ----- CALCULATE Q^Z POPULATIONS -----
C     Q_A^Z = 2 \sum_{\mu \in A} \sum_\nu Z_{\mu \nu}' S_{\mu \nu}
C
C     USE "NUM" WHOSE DIMENSION CORRESPONDS TO IFG'S NUM
      CALL DFTB_MULLIKEN(NUM,NUM*(NUM+1)/2,Z,S,X(LWRK))
      CALL DFTB_MULSA(X(LWRK),X(LWRK1),X(LQZ),NUM,NSHELL,NAT)
C
      NATI = NATFRG(IFG)
      INDI = INDFRG(IFG)
      ONE_THIRD = 1.0D+00/3.0D+00
C
      NATFMOB = NATFMO + NBDFG
      DO JFG = 1, NFG !! CORRESPONDS TO KFG IN OUR EQUATIONS
C       IF (JFG.EQ.IFG) CYCLE
        if(dotrunces) then
          lcradfgi=lcradfg+(ifg-1)*4
          lcradfgj=lcradfg+(jfg-1)*4
          iich=ixftch(x(lichfgg),ifg)
          jjch=ixftch(x(lichfgg),jfg)
          if(iich.eq.0.and.jjch.eq.0.and.fmodest(x(lcradfgi),
     *       x(lcradfgi+3),x(lcradfgj),x(lcradfgj+3)).gt.rcutes/10) then
c           write(6,*) 'Skipping',ifg,jfg
            cycle
          endif
        endif
        INDJ = INDFRG(JFG)
        NATJ = NATFRG(JFG)
C
C       ----- CALCULATE GAMMA DERIVATIVE TERMS -----
C
        NSEQ = 0
        NSEQ2 = 0
        DO I = 1, NATI
          IAG = IATFRG(INDI+I-1)
          ISP = ISPEFMO(IAG)
          CHAMUL1 = X(LQZ+I-1)
          IF (DFTB3) CHNI = POPMAT(I,IFG) - X(LZREF+INDI+I-2)
          DO J = 1, NATJ
            NSEQ = NSEQ + 1
            JAG = IATFRG(INDJ+J-1)
            JSP = ISPEFMO(JAG)
            VALTMP = ZERO
            CHAMUL2 = POPMAT(J,JFG) - X(LZREF+INDJ+J-2) !! CHARGE OF J FRAGMENT
            VEC(1) = FMOC(1,IAG) - FMOC(1,JAG)
            VEC(2) = FMOC(2,IAG) - FMOC(2,JAG)
            VEC(3) = FMOC(3,IAG) - FMOC(3,JAG)
            DIST = SQRT(VEC(1)**2+VEC(2)**2+VEC(3)**2)
            IF (IAG.NE.JAG) DISTI = ONE/DIST
            DAMPXH = DAMPXH0
            CALL DFTB_DMPCHK(ISP,JSP,DAMPXH)
            IF (IAND(MODGAMMA,3).EQ.0) THEN
              IF (PERIOD) THEN
                CALL DFTB_GAMFMO_PBC(IAG,JAG,FMOC(1,IAG),FMOC(1,JAG),
     *                               GAMTMP,GAMI,GAMJ,TMP(1,ISP,JSP),
     *                               TMP3(1,ISP,JSP),TMP3(1,JSP,ISP),
     *                               DAMPXH,ISP,JSP)
              ELSE
                GAMTMP = DFTB_GAM2(DIST,DISTI,TMP(1,ISP,JSP),
     *                             DAMPXH,DAMPXHE,SPE(ISP),SPE(JSP))
              END IF
            ELSE IF (IAND(MODGAMMA,1).NE.0) THEN
              CALL DFTB_CNVSQ(IAG,JAG,NSEQ2)
              GAMTMP = X(LGAMMAGLOB+NSEQ2-1)
            ELSE IF (IAND(MODGAMMA,2).NE.0) THEN
              GAMTMP = X(LGAMMAGLOB+INDJ+J-2+NATFMOB*(INDI+I-2))
            END IF
            IF (DFTB3) THEN
              IF (IAND(MODGAMMA,3).EQ.0) THEN
                IF (PERIOD) THEN
                ELSE
                  CALL DFTB_GAM3(DIST,DISTI,GAMI,GAMJ,HUBB(ISP),
     *              HUBB(JSP),TMP3(1,ISP,JSP),TMP3(1,JSP,ISP),DAMPXH,
     *              DAMPXHE,SPE(ISP),SPE(JSP))
                  GAMI = GAMI*HUBDER(ISP)
                  GAMJ = GAMJ*HUBDER(JSP)
                END IF
              ELSE IF (IAND(MODGAMMA,1).NE.0) THEN
                GAMI = X(LGAMMA3+IAG-1+(JAG-1)*NATFMO) !! G_ac
                GAMJ = X(LGAMMA3+JAG-1+(IAG-1)*NATFMO) !! G_ca
              ELSE IF (IAND(MODGAMMA,2).NE.0) THEN
                II = INDI+I-2
                JJ = INDJ+J-2
                IIJJ = II+NATFMOB*JJ
                gami = x(lgamma3+iijj)
                gamj = x(lgamma3+iijj+natfmob*natfmob)
              END IF
              GAMTMP = GAMTMP + TWO*ONE_THIRD*(GAMI*CHNI+GAMJ*CHAMUL2)
            END IF
            IF (IFG.NE.JFG)
     *        SHIFTZ(INDJ+J-1) = SHIFTZ(INDJ+J-1) + GAMTMP*CHAMUL1
            IF (IAG.EQ.JAG.AND..NOT.PERIOD) CYCLE
C           DIST = SQRT(VEC(1)**2+VEC(2)**2+VEC(3)**2)
            IF (PERIOD) THEN
              CALL DFTB_GAMDER_PBC(IAG,JAG,FMOC(1,IAG),FMOC(1,JAG),
     *                             GDV2,GDV3I,GDV3J,GDST2,GDST3I,GDST3J,
     *                             TMPD(1,ISP,JSP),TMPD3(1,ISP,JSP),
     *                             TMPD3(1,JSP,ISP),DAMPXH,ISP,JSP)
              !! Stress tensor
              VAL = CHAMUL1*CHAMUL2
              CALL DAXPY(9,VAL,GDST2,1,STRFMO,1)
              IF (DFTB3) THEN
                VAL = TWO*ONE_THIRD*CHNI*CHAMUL1*CHAMUL2
                CALL DAXPY(9,VAL,GDST3I,1,STRFMO,1)
                VAL = ONE_THIRD*CHAMUL1*CHAMUL2*CHAMUL2
                CALL DAXPY(9,VAL,GDST3J,1,STRFMO,1)
              END IF
            ELSE
C             GAMMAV = GAMD_FMO2(DIST,DISTI,TMPD(1,ISP,JSP),DAMPXH,
C    *                           DAMPXHE,SPE(ISP),SPE(JSP))
C             IF (DFTB3) THEN
C              CALL DFTB_GAMD3(DIST,DISTI,GAMI,GAMJ,HUBB(ISP),HUBB(JSP),
C    *            HUBDER(ISP),HUBDER(JSP),TMPD3(1,ISP,JSP),
C    *            TMPD3(1,JSP,ISP),DAMPXH,DAMPXHE,SPE(ISP),SPE(JSP))
C               VALTMP = ONE_THIRD*(CHNI*GAMI*TWO+CHAMUL2*GAMJ)
C               GAMMAV = GAMMAV + VALTMP
C             END IF
              CALL DFTB_GAM23D(DIST,DISTI,GAMMAV,GAMI,GAMJ,
     *                         TMPD(1,ISP,JSP),TMPD3(1,ISP,JSP),
     *                         TMPD3(1,JSP,ISP),DAMPXH,DAMPXHE,
     *                         SPE(ISP),SPE(JSP),HUBDER(ISP),
     *                         HUBDER(JSP),DFTB3)
              DO K = 1, 3
                GDV2(K) = GAMMAV*VEC(K)*DISTI
                IF (DFTB3) THEN
                  GDV3I(K) = GAMI*VEC(K)*DISTI
                  GDV3J(K) = GAMJ*VEC(K)*DISTI
                END IF
              END DO
            END IF
            DO K = 1, 3
C             VAL = CHAMUL1 * CHAMUL2 * GAMMAV*DISTI
              VAL = CHAMUL1*CHAMUL2*GDV2(K)
              IF (DFTB3) THEN
                VAL = VAL + ONE_THIRD*(CHNI*GDV3I(K)*TWO
     *                      +CHAMUL2*GDV3J(K))*CHAMUL1*CHAMUL2
              END IF
              FMODE(K,IAG) = FMODE(K,IAG) + VAL
              FMODE(K,JAG) = FMODE(K,JAG) - VAL
            END DO
C           VAL = CHAMUL1 * CHAMUL2 * GAMMAV*DISTI
C           DO K = 1, 3
C             FMODE(K,IAG) = FMODE(K,IAG) - VEC(K)*VAL
C             FMODE(K,JAG) = FMODE(K,JAG) + VEC(K)*VAL
C           END DO
          END DO
        END DO
      END DO !! END OF KFG LOOP
C
      CALL RETFM(NEED)
C
      RETURN
C
      END SUBROUTINE DFTB_ZVEC_KFG
C
C-----------------------------------------------------------------------
C*MODULE DFTBFO    *DECK DFTB_ZVEC_MUL
C>
C>    @brief   Fully analytic FMO-DFTB gradient
C>
C>    @details Calculate Mulliken charges in FMO-DFTB SCZV
C>
C>    @author  Yoshio Nishimoto
C>    - Feb, 2016- Subroutine written
C>
C>           --- INPUT ---
C>    @param NAT     Number of atoms
C>    @param LPOPMAT Where POPMAT is
C>    @param MAXNAT  Maximum number of atoms in fragment
C>    @param IFG     Current fragment
C>           --- IN/OUTPUT ---
C>    @param CHAMUL  Mulliken population
C>    @param LCHAMULI ...
C>
C
      SUBROUTINE DFTB_ZVEC_MUL(NAT,LPOPMAT,MAXNAT,IFG,CHAMUL,LPOPMATI)
      use mx_limits, only: mxatm
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      INTEGER,PARAMETER :: MXSPE=10
C
      COMMON /DFTBPR/ ETEMP,DFTBDP(MXSPE*14),DAMPXHE,HUBDER(MXSPE),
     *                ZREF(MXATM),SPNCST(6,MXSPE),SPE(MXATM),NSPE,
     *                MAXANG(MXATM),ISPE(MXATM),IND(MXATM+1),IDFTBD,
     *                PARAMDIR
      COMMON /FMCOM / X(1)
      COMMON /FMOINF/ NFG,NLAYER,NATFMO,NBDFG,NAOTYP,NBODY,NSEGM
      common /fmorun/ espscf,e0scf(2),emp2s,IDAFMO,icurfg,jcurfg,kcurfg,
     *                icurlay,icurunt,nat1e,ncursh,ngau,icurpop,ifmostp,
     *                moncor,needr,modrst,norbproj,nunesp,iskipesp,
     *                IESDPPC,idoprop,mp2run,icurit,idmfmo,iddfmo,
     *                iddcur,nddleft,ivmfmo,nzmtfmo,ifmobas,itmfmo(2)
C
      DIMENSION CHAMUL(NAT)
C
      lpopmati=lpopmat+(maxnat*nfg)*(icurpop-1)
      !! MULLIKEN CHARGES IN LCHAMUL
      CALL DCOPY(NAT,X(LPOPMATI+(IFG-1)*MAXNAT),1,CHAMUL,1)
      DO I = 1, NAT
        CHAMUL(I) = CHAMUL(I) - ZREF(I)
      END DO
C
      RETURN
C
      END SUBROUTINE DFTB_ZVEC_MUL
C
C-----------------------------------------------------------------------
C*MODULE DFTBFO    *DECK DFTB_ZVEC_KGRAD
C>
C>    @brief   Fully analytic FMO-DFTB gradient
C>
C>    @details Calculate gradient contributions in KFG loop
C>
C>    @author  Yoshio Nishimoto
C>    - Feb, 2016- Subroutine written
C>    @date    Feb, 2018 - Yoshio Nishimoto
C>    - Added FMO-DFTB/AFO
C>
C>           --- INPUT ---
C>    @param NATFRG  Number of atoms in each fragment
C>    @param IATFRG  Global atom index with NBDFG?
C>    @param INDFRG  First atom of each fragment in terms of NATFMOB?
C>    @param ISPEFMO ISPE for all atoms (NATFMO)
C>    @parma IDMREC0
C>    @param NUMFRG  Nuber of AOs in fragment
C>    @param IODFMO 
C>    @param FMOC    Coordinates for all NATFMO
C>    @param MODORB  Do we have orbitals?
C>    @param MDOUT   Controls output length
C>    @param IAGLOB  Global atom index
C>    @param IATFRG  Global atom index with NBDFG?
C>    @param IALMO   Index of adaptive LMO?
C>    @param FMOAFO  FMO/AFO or not
C>    @param ISWNEW1 PCM<1> option
C>           --- IN/OUTPUT ---
C>    @param FMODE   Gradient
C>
C
      SUBROUTINE DFTB_ZVEC_KGRAD(LAYFRG,NATFRG,IATFRG,
     *  INDFRG,FMODE,IDMREC0,NUMFRG,IODFMO,ISPEFMO,FMOC,MODORB,MDOUT,
     *  IAGLOB,IANFRG,IALMO,FMOAFO,ISWNEW1)
      use mx_limits, only: mxatm,mxsh,mxao,mxgtot
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      INTEGER,PARAMETER :: MXSPE=10
      DOUBLE PRECISION, PARAMETER :: ONE=1.0D+00,TWO=2.0D+00,
     *  ZERO=0.0D+00
      integer ddi_world,ddi_group
      PARAMETER (ddi_world=0)
      PARAMETER (ddi_group=1)
C
      logical GOPARR,DSKWRK,MASWRK
      LOGICAL MYJOB,ISGDDI,PAROUT,INITGDDI,WASGDDI,MLGDDI
      LOGICAL MFRZ,SKIPOR
C
      COMMON /DFTBPR/ ETEMP,DFTBDP(MXSPE*14),DAMPXHE,HUBDER(MXSPE),
     *                ZREF(MXATM),SPNCST(6,MXSPE),SPE(MXATM),NSPE,
     *                MAXANG(MXATM),ISPE(MXATM),IND(MXATM+1),IDFTBD,
     *                PARAMDIR
      COMMON /FMCOM / X(1)
      COMMON /FMOAFO/ LINDFRZ,LGMK,LGMKSAV
      COMMON /FMODTB/ LGAMMAGLOB,LZREF,LISPE,LSHIFTG,LSHIFTCT,
     *                LGAMMAD,MODESD,MODGAMMA,LDCTMP,LGAMMA3,LGAMMA3D,
     *                LNBOND,LCTMUL,LESPDFTB,LESPASC,NDFTBRST,LDFTBRST
      common /fmodim/ maxbnd,maxknd,maxcbs,maxcao,maxvec,maxl1,maxnat,
     *                maxabd,maxbas,maxbbd,maxlmo,maxslo,maxabd2,maxrij
      COMMON /FMOINF/ NFG,NLAYER,NATFMO,NBDFG,NAOTYP,NBODY,NSEGM
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
      COMMON /MFRPAR/ MFRZ,NUMFRZ,NORFRZ,IFRZ(MXAO)
      COMMON /NSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),
     *                CF(MXGTOT),CG(MXGTOT),CH(MXGTOT),CI(MXGTOT),
     *                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),
     *                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
C
      INTEGER INDFRG(*),IATFRG(*),NATFRG(*),LAYFRG(*),IDMREC0(*),
     *  NUMFRG(*),ISPEFMO(*),IODFMO(*),MODORB,IAGLOB(*),IANFRG(*),
     *  IALMO(*)
      DOUBLE PRECISION FMODE(3,NATFMO),FMOC(3,NATFMO)
      LOGICAL ORBXCH,MDOUT,FMOAFO
C
      TIM0 = ZERO
      CALL TSECND(TIM0)
C
      L1 = MAXL1
      L2 = (L1*L1+L1)/2
      L3 = L1*L1
C
C     ----- ALLOCATE MEMORY -----
C     USE THE SIZE OF THE BIGGEST FRAGMENT
C
      CALL VALFM(LOADFM)
      CALL GOTFM(NGOTMX)
      LWRK    = LOADFM  + 1
      LVEC    = LWRK    + L1
      LWRK1   = LVEC    + L3
      LWRK2   = LWRK1   + L3
      LOCC    = LWRK2   + L3
      LD      = LOCC    + L1
      LS      = LD      + L2
      LDS     = LS      + L2
      LPS     = LDS     + L3 !! D SQUARE
      LSHIFT  = LPS     + L2*MAXNAT*3
      LAST    = LSHIFT  + MAXNAT
C
      IF (FMOAFO) THEN
        LVECL = LAST
        LSCR  = LVECL + L3
        LE    = LSCR  + L1*8
        LIWRK = LE    + L1
        LAST  = LIWRK + L1
C
        LWRK3 = LDS
        LFAO  = LD
      END IF
C
      NEED    = LAST - LOADFM - 1
C
C     ----- GET MEMORY -----
C
      CALL GETFM(NEED)
C
C     ----- RESET GLOBAL COUNTER -----
C
      IF (ISGDDI) THEN
        CALL GDDI_SCOPE(DDI_WORLD)
        CALL DDI_GSUMF(2424,X(LSHIFTCT),NATFMO+NBDFG)
        IF (ISWNEW1.EQ.1) CALL DDI_GSUMF(2424,X(LSHIFTG),NATFMO+NBDFG)
        CALL GDDI_SCOPE(DDI_GROUP)
        CALL GDDICOUNT(-1,LGROUP,MYJOB)
      END IF
      ORBXCH=MOD(MODORB,2).NE.0
C
      DO IFG = 1, NFG !! CORRESPONDS TO KFG IN OUR EQUATIONS
        IF (ISGDDI) THEN
          call GDDICOUNT(0,LGROUP,MYJOB)
          IF (.NOT.MYJOB) CYCLE
        END IF
        ILAY=LAYFRG(IFG)
        NUM = IAND(NUMFRG(IFG),65535)
        L1 = NUM
        L2 = (L1*L1+L1)/2
        L3 = L1*L1
C
        NAT = NATFRG(IFG)
        NAL = ISHFT(NUMFRG(IFG),-16)
        NBL = NAL
        NA  = NAL
        IREC0=IDMREC0(ILAY)
        CALL VCLR(X(LD),1,L2)
        IF (ORBXCH) THEN
          CALL RAREADS(IDAFMO,IODFMO,X(LVEC),L3,IREC0+IFG,0)
          CALL VCLR(X(LOCC),1,L1)
          IF (NAL.EQ.NBL) THEN
            DO I = 1, NAL
              X(LOCC-1+I) = TWO
            END DO
          ELSE
            DO I = 1, NBL
              X(LOCC-1+I) = TWO
            END DO
            DO I = NBL+1, NAL
              X(LOCC-1+I) = ONE
            END DO
          END IF
          CALL DFTB_DMTOCC(X(LD),X(LVEC),X(LOCC),X(LWRK2),L1,L1,1.0D-09)
        ELSE
          CALL RAREADS(IDAFMO,IODFMO,X(LD),L2,IREC0+IFG,0)
        END IF
C
C       ----- PREPARE INFORMATION (INSTEAD OF CALLING MAKEMOL) -----
C
        INDI = INDFRG(IFG)
        NATI = NATFRG(IFG)
        NAT  = NATI
        IND(1) = 0
        DO I = 1, NATI
          IAG = IATFRG(INDI+I-1)
          CALL DCOPY(3,FMOC(1,IAG),1,C(1,I),1)
          ISPE(I) = ISPEFMO(IAG)
          IND(I+1) = IND(I) + MAXANG(ISPE(I))**2
          IAGLOB(I) = IAG
        END DO
C
C       ----- MAKE SHIFT ARRAY -----
C
        CALL VCLR(X(LSHIFT),1,MAXNAT)
        CALL DCOPY(NAT,X(LSHIFTCT+INDI-1),1,X(LSHIFT),1)
C
C       ----- CALCULATE OVERLAP DERIVATIVES -----
C
        CALL DFTB_EXTSDER(X(LPS),L2,1,NAT)
C
C       ----- CALCULATE OVERLAP INTEGRAL FOR MONOMER IFG -----
C
        CALL VCLR(X(LS),1,L2)
        CALL DFTB_HAM_OVER(X(LWRK2),X(LS),.FALSE.,.TRUE.,.FALSE.)
C
C       ----- PREPARE DS+WS TERM -----
C
        call dftb_addlag(0,l1,x(ld),x(ls),x(lds),x(lwrk1),x(lwrk2),
     *                   x(lshift),x(lshift),nat,c,fmode,iatfrg(indi))
      ! CALL CPYTSQ(X(LD),X(LDS),L1,1)
      ! CALL CPYTSQ(X(LS),X(LWRK2),L1,1)
      ! !! X(LWRK1) = X(LWRK2) * X(LDS) = LS * LD
      ! CALL MRARBR(X(LWRK2),L1,L1,L1,X(LDS),L1,L1,X(LWRK1),L1) !! (SD)
      ! CALL DS_SHIFT(L1,NAT,X(LDS),X(LSHIFT))
      ! !! X(LWRK2) = X(LDS) * X(LWRK1) = (LD*SHIFT) * (LS*LD)
      ! CALL MRARBR(X(LDS),L1,L1,L1,X(LWRK1),L1,L1,X(LWRK2),L1) !! (DSD+SHIFT)
      ! CALL DAXPY(L3,-0.5D+00,X(LWRK2),1,X(LDS),1)
C
C     ! ----- CALCULATE GRADIENT OF KFG -----
C
      ! IF (MASWRK) CALL ADD_LAG_GRAD(L1,L2,NAT,NATFMO,IATFRG(INDI),
     *!                   FMODE,X(LPS),X(LDS),X(LWRK1))
C
        IF (FMOAFO) THEN
          L0=L1
          ICURFG = IFG
          IODA(12)=0
          IODA(15)=0
          IODA(318)=0
          IODA(569)=0
C
          CALL DCOPY(NAT,X(LZREF+INDI-1),1,ZAN,1)
          CALL ICOPY(NAT,IANFRG(INDI),1,IAN,1)
          NSHELL = 0
          KLOC(1) = 1
          DO I = 1, NAT
            DO J = 1, MAXANG(ISPE(I))
              NSHELL = NSHELL + 1
              KATOM(NSHELL) = I
              KTYPE(NSHELL) = J
              IF (J.EQ.1) THEN
                KMIN(NSHELL) = 1
                KMAX(NSHELL) = 1
              ELSE IF (J.EQ.2) THEN
                KMIN(NSHELL) = 2
                KMAX(NSHELL) = 4
              ELSE IF (J.EQ.3) THEN
                KMIN(NSHELL) = 5
                KMAX(NSHELL) = 9
              END IF
              IF (NSHELL.NE.1) KLOC(NSHELL) = KLOC(NSHELL-1)
     *          + KMAX(NSHELL-1) - KMIN(NSHELL-1) + 1
            END DO
          END DO
          CALL DAWRIT(IDAF,IODA,X(LS),L2,12,0)
          CALL QMTSYM(X(LS),X(LWRK1),X(LWRK2),X(LE),X(LSCR),
     *                X(LIWRK),L0,L1,L2,L3,.FALSE.)
          CALL FMOORBS(X(LVECL),X(LS),X(LE),X(LWRK2),X(LSCR),
     *                 X(LIWRK),L0,L1,L2,L3,SKIPOR,1)
C
          IF (NUMFRZ.EQ.0) CYCLE
C
          NFRZ = MAX(NUMFRZ,NORFRZ)
          CALL VALFM(LOADFM)
          LUROT   = LOADFM  + 1
          LRROT   = LUROT   + NFRZ*NFRZ
          LAST    = LRROT   + NFRZ*NFRZ
          NEED2   = LAST - LOADFM - 1
          CALL GETFM(NEED2)
C
          CALL DAREAD(IDAF,IODA,X(LUROT),NFRZ*NFRZ,569,0)
          CALL DFTB_RROT(L1,NFRZ,NUMFRZ,X(LUROT),X(LRROT),X(LS),
     *                   X(LVECL),X(LVEC),X(LWRK1),X(LWRK2))
C
          CALL VCLR(X(LFAO),1,L2)
          IF (ISWNEW1.EQ.1)
     *      CALL DAXPY(NAT,ONE,X(LSHIFTG+INDI-1),1,X(LSHIFT),1)
          CALL DFTB_SHIFT_FOCK(X(LSHIFT),X(LFAO),X(LS),NAT,L2,
     *         NSPE,ISPE,IND,MAXANG)
          CALL VCLR(X(LGMK),1,L1*NFRZ)
          CALL DFTB_CNSTGMK(L1,NA,NUMFRZ,NFRZ,ONE,
     *                     X(LGMK),X(LS),X(LFAO),X(LVEC),X(LRROT),
     *                     X(LWRK1),X(LWRK2),X(LWRK3))
C
          JCURFG = -1
          CALL DFTB_PUTGMK(L1,1,X(LGMK),X(LGMKSAV),X(LINDFRZ),IALMO,
     *                     IAGLOB,X(LAST),X(LAST),X(LAST),.FALSE.,
     *                     X(LAST),X(LAST),X(LAST))
          JCURFG = 0
C
          CALL RETFM(NEED2)
        END IF
      END DO !! END OF IFG LOOP
C
      IF (ISGDDI) CALL GDDICOUNT( 1,LGROUP,MYJOB)
      CALL RETFM(NEED)
C
      IF (MDOUT) THEN
        IF(MASWRK) THEN
          IW = 6
          WRITE (IW,*)
          WRITE(IW,FMT='('' ...... END OF FMO-DFTB SCZV ......'')')
        END IF
        CALL TIMIT(1)
      END IF
C
      RETURN
C
      END SUBROUTINE DFTB_ZVEC_KGRAD
C
C-----------------------------------------------------------------------
C*MODULE DFTBFO    *DECK DFTB_FMOGRAD_EMBED
C>
C>    @brief   FMO-DFTB gradient
C>
C>    @details Calculate gradient contributions of embedding charges
C>
C>    @author  Yoshio Nishimoto
C>    - Feb, 2016- Subroutine written
C>
C>           --- INPUT ---
C>    @param NAT     Number of atoms
C>    @param NATFMO  Number of atoms in the whole system
C>    @param NBDFG   Number of BDA-BAA pairs
C>    @param ISPEFMO ISPE for all atoms (NATFMO)
C>    @param IAGLOB  Global atom index
C>    @param INDAT   Converts global atomic location to fragment index
C>    @param IALOC   Converts global to local (fragment) atomic location
C>    @param INDFRG  First atom of each fragment in terms of NATFMOB?
C>    @param NN      Length of an array
C>    @param NN2     Length of an array
C>    @param INDJ    INDFRG(JFG)?
C>    @param SHIFTCT Shift contributions due to charge transfer
C>    @param CTMUL   Sum of charge transfer
C>    @param FMOC    Coordinates for all NATFMO
C>    @param POPMAT  Mulliken population
C>    @param ZREFFMO ZREF for all atoms (NATFMOB length)
C>    @param TMP     Some constants for DFTB2 gamma
C>    @param TMP3    Some constants for DFTB3 gamma
C>    @param TMPD    Some constants for derivative of DFTB2 gamma
C>    @param TMPD3   Some constants for derivative of DFTB3 gamma
C>    @param GAMMA   DFTB2 gamma
C>    @param GAMMA3  DFTB3 gamma
C>    @param MASWRK  Master!
C>    @param NATFRG  Number of atoms in each fragment
C>    @param IATFRG  Global atom index with NBDFG?
C>    @param NFG     Number of fragments
C>           --- IN/OUTPUT ---
C>    @param FMODE   Gradient
C>
C
      SUBROUTINE DFTB_FMOGRAD_EMBED(NAT,NATFMO,NBDFG,ISPEFMO,IAGLOB,
     *  INDAT,IALOC,INDFRG,NN,NN2,INDJ,SHIFTCT,FMODE,CTMUL,FMOC,POPMAT,
     *  ZREFFMO,TMP,TMP3,TMPD,TMPD3,GAMMA,GAMMA3,MASWRK,
     *  NATFRG,IATFRG,NFG,ich)
      use mx_limits, only: mxatm
      USE DFTBPB_MOD,ONLY: PERIOD,DFTB_GAMFMO_PBC,DFTB_GAMDER_PBC,STRFMO
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      INTEGER,PARAMETER :: MXSPE=10
      DOUBLE PRECISION, PARAMETER :: ONE=1.0D+00,TWO=2.0D+00,
     *  ZERO=0.0D+00,ONE_THIRD=ONE/3.0D+00
C
      LOGICAL DFTBFL,SCC,SRSCC,DFTB3,DAMPXH0,LCDFTB,dotrunces,dotruncesd
C
      INTEGER, INTENT(IN) :: ISPEFMO(NATFMO),IAGLOB(NAT),INDAT(*),
     *  IALOC(*),INDFRG(*),NATFRG(*),IATFRG(*)
      DOUBLE PRECISION, INTENT(IN) :: CTMUL(*),FMOC(3,*),
     *  POPMAT(NAT),ZREFFMO(NAT),GAMMA(*),GAMMA3(NN2,NN2,*)
      DOUBLE PRECISION, INTENT(INOUT) :: SHIFTCT(NAT),FMODE(3,NATFMO)
C
      COMMON /DFTB  / DFTBFL,SCC,SRSCC,DFTB3,DAMPXH0,LCDFTB
      COMMON /DFTBPR/ ETEMP,DFTBDP(MXSPE*14),DAMPXHE,HUBDER(MXSPE),
     *                ZREF(MXATM),SPNCST(6,MXSPE),SPE(MXATM),NSPE,
     *                MAXANG(MXATM),ISPE(MXATM),IND(MXATM+1),IDFTBD,
     *                PARAMDIR
      COMMON /DFTBSK/ SKDIM(MXSPE,MXSPE),SKSPIN(MXSPE),QREFL(3,MXSPE),
     *                HUBBL(3,MXSPE),QREF(MXSPE),HUBB(MXSPE),SKCUT2,
     *                NEEDSK,LSKHTAB(MXSPE,MXSPE),LSKSTAB(MXSPE,MXSPE),
     *                LSKGRID(MXSPE,MXSPE),LSKSELF(MXSPE)
      COMMON /FMCOM / X(1)
      COMMON /FMODTB/ LGAMMA,LZREF,LISPE,LSHIFTG,LSHIFTCT,
     *                LGAMMAD,MODESD,MODGAMMA,LDCTMP,LGAMMA3,LGAMMA3D,
     *                LNBOND,LCTMUL,LESPDFTB,LESPASC,NDFTBRST,LDFTBRST
      common /cenrad/ rcutes,lcradfg,lichfgg,dotrunces,dotruncesd
C
      DOUBLE PRECISION :: VEC(3),TMPD(10,NSPE,NSPE),TMP(12,NSPE,NSPE),
     *  TMP3(10,NSPE,NSPE),TMPD3(11,NSPE,NSPE),cencur(3)
      LOGICAL :: MASWRK,DAMPXH
      !! Gamma derivative vector and stress tensor for PBC
      DIMENSION GDV2(3),GDV3I(3),GDV3J(3),
     *          GDST2(3,3),GDST3I(3,3),GDST3J(3,3)
C
C     IAGLOB IS ONLY FOR THE CURRENT FRAGMENT
C
C     SCAL=ONE
C     IF (.NOT.MASWRK) SCAL=ZERO
      NATFMOB = NATFMO + NBDFG
C
      TMPSH2 = ZERO
      TMPSH3 = ZERO
      VALTMP = ZERO
      IF (IAND(MODGAMMA,3).EQ.0) THEN
        if(dotrunces) call radcent(cencur,radcur)
        DO LFG = 1, NFG
          if(dotrunces.and.ich.eq.0) then
            lcradfgl=lcradfg+(lfg-1)*4
            lch=ixftch(x(lichfgg),lfg)
            if(lch.eq.0.and.fmodest(cencur,radcur,x(lcradfgl),
     *          x(lcradfgl+3)).gt.rcutes) then
c             write(6,*) 'Skipping',ifg,jfg,lfg
              cycle
            endif
          endif
          NATI = NATFRG(LFG)
          INDI = INDFRG(LFG)
          DO IAT = 1, NATI
            IAG = IATFRG(INDI+IAT-1)
            CHAMUL1 = CTMUL(INDI+IAT-1)
            IF (DFTB3) CHAMULSQ = CTMUL(NATFMOB+INDI+IAT-1)
            ISP = ISPEFMO(IAG)
            DO J = 1, NAT
              JAG = IAGLOB(J)
              JSP = ISPEFMO(JAG)
              CHAMUL2 = POPMAT(J) - ZREFFMO(J)
              VEC(1) = FMOC(1,IAG) - FMOC(1,JAG)
              VEC(2) = FMOC(2,IAG) - FMOC(2,JAG)
              VEC(3) = FMOC(3,IAG) - FMOC(3,JAG)
              DAMPXH = DAMPXH0
              CALL DFTB_DMPCHK(ISP,JSP,DAMPXH)
              IF (PERIOD) THEN
                CALL DFTB_GAMFMO_PBC(IAG,JAG,FMOC(1,IAG),FMOC(1,JAG),
     *                               GAMTMP,GAMI,GAMJ,TMP(1,ISP,JSP),
     *                               TMP3(1,ISP,JSP),TMP3(1,JSP,ISP),
     *                               DAMPXH,ISP,JSP)
              ELSE
                DIST = SQRT(VEC(1)*VEC(1)+VEC(2)*VEC(2)+VEC(3)*VEC(3))
                IF (DIST.GE.1.0D-15) DISTI = ONE/DIST
                CALL DFTB_GAM23(1,DIST,GAMTMP,GAMI,GAMJ,
     *               TMP(1,ISP,JSP),TMP3(1,ISP,JSP),TMP3(1,JSP,ISP),
     *               DAMPXH,DAMPXHE,SPE(ISP),SPE(JSP),HUBDER(ISP),
     *               HUBDER(JSP),DFTB3)
              END IF
              SHIFTCT(J) = SHIFTCT(J) + GAMTMP*CHAMUL1
              IF (DFTB3) THEN
                SHIFTCT(J) = SHIFTCT(J)
     *            + ONE_THIRD*(CHAMULSQ*GAMI+TWO*CHAMUL1*CHAMUL2*GAMJ)
              END IF
C
              IF (.NOT.MASWRK .OR. IAG.EQ.JAG.AND..NOT.PERIOD) CYCLE
C
              IF (PERIOD) THEN
                CALL DFTB_GAMDER_PBC(IAG,JAG,FMOC(1,IAG),FMOC(1,JAG),
     *                               GDV2,GDV3I,GDV3J,GDST2,GDST3I,
     *                               GDST3J,TMPD(1,ISP,JSP),
     *                               TMPD3(1,ISP,JSP),TMPD3(1,JSP,ISP),
     *                               DAMPXH,ISP,JSP)
                !! Stress tensor
                VAL = CHAMUL1*CHAMUL2
                CALL DAXPY(9,VAL,GDST2,1,STRFMO,1)
                IF (DFTB3) THEN
                  VAL = ONE_THIRD*CHAMULSQ*CHAMUL2
                  CALL DAXPY(9,VAL,GDST3I,1,STRFMO,1)
                  VAL = ONE_THIRD*CHAMUL1*CHAMUL2*CHAMUL2
                  CALL DAXPY(9,VAL,GDST3J,1,STRFMO,1)
                END IF
                IF (IAG.EQ.JAG) CYCLE
              ELSE
                CALL DFTB_GAM23D(DIST,DISTI,GAMMAV,GAMI,GAMJ,
     *                           TMPD(1,ISP,JSP),TMPD3(1,ISP,JSP),
     *                           TMPD3(1,JSP,ISP),DAMPXH,DAMPXHE,
     *                           SPE(ISP),SPE(JSP),HUBDER(ISP),
     *                           HUBDER(JSP),DFTB3)
                DO K = 1, 3
                  GDV2(K) = GAMMAV*VEC(K)*DISTI
                  IF (DFTB3) THEN
                    GDV3I(K) = GAMI*VEC(K)*DISTI
                    GDV3J(K) = GAMJ*VEC(K)*DISTI
                  END IF
                END DO
              END IF
              DO K = 1, 3
                VAL = CHAMUL1*CHAMUL2*GDV2(K)
                IF (DFTB3) THEN
                  VAL = VAL
     *              + ONE_THIRD*(CHAMULSQ*GDV3I(K)
     *                          +CHAMUL1*CHAMUL2*GDV3J(K))*CHAMUL2
                END IF
                FMODE(K,IAG) = FMODE(K,IAG) + VAL
                FMODE(K,JAG) = FMODE(K,JAG) - VAL
              END DO
            END DO
          END DO
        END DO
      ELSE IF (IAND(MODGAMMA,1).NE.0) THEN
        DO LFG = 1, NFG
          NATI = NATFRG(LFG)
          INDI = INDFRG(LFG)
          DO IAT = 1, NATI
            IAG = IATFRG(INDI+IAT-1)
            CHAMUL1 = CTMUL(INDI+IAT-1)
            IF (DFTB3) CHAMULSQ = CTMUL(NATFMOB+INDI+IAT-1)
            ISP = ISPEFMO(IAG)
            DO J = 1, NAT
              JAG = IAGLOB(J)
              JSP = ISPEFMO(JAG)
              CHAMUL2 = POPMAT(J) - ZREFFMO(J)
              VEC(1) = FMOC(1,IAG) - FMOC(1,JAG)
              VEC(2) = FMOC(2,IAG) - FMOC(2,JAG)
              VEC(3) = FMOC(3,IAG) - FMOC(3,JAG)
              CALL DFTB_CNVSQ(IAG,JAG,NSEQ)
              GAMTMP = GAMMA(NSEQ)
              TMPSH  = GAMTMP*CHAMUL1
              DAMPXH = DAMPXH0
              CALL DFTB_DMPCHK(ISP,JSP,DAMPXH)
              IF (DFTB3) THEN
                GAMI = GAMMA3(IAG,JAG,1)
                GAMJ = GAMMA3(JAG,IAG,1)
                TMPSH2 = CHAMULSQ*GAMI
                TMPSH3 = CHAMUL2*CHAMUL1*GAMJ
                TMPSH = TMPSH
     *            + ONE_THIRD*(CHAMULSQ*GAMI+TWO*CHAMUL1*CHAMUL2*GAMJ)
              END IF
              SHIFTCT(J) = SHIFTCT(J) + TMPSH
C
              IF (.NOT.MASWRK .OR. IAG.EQ.JAG.AND..NOT.PERIOD) CYCLE
C
              IF (PERIOD) THEN
                CALL DFTB_GAMDER_PBC(IAG,JAG,FMOC(1,IAG),FMOC(1,JAG),
     *                               GDV2,GDV3I,GDV3J,GDST2,GDST3I,
     *                               GDST3J,TMPD(1,ISP,JSP),
     *                               TMPD3(1,ISP,JSP),TMPD3(1,JSP,ISP),
     *                               DAMPXH,ISP,JSP)
                !! Stress tensor
                VAL = CHAMUL1*CHAMUL2
                CALL DAXPY(9,VAL,GDST2,1,STRFMO,1)
                IF (DFTB3) THEN
                  VAL = ONE_THIRD*CHAMULSQ*CHAMUL2
                  CALL DAXPY(9,VAL,GDST3I,1,STRFMO,1)
                  VAL = ONE_THIRD*CHAMUL1*CHAMUL2*CHAMUL2
                  CALL DAXPY(9,VAL,GDST3J,1,STRFMO,1)
                END IF
                IF (IAG.EQ.JAG) CYCLE
              ELSE
                DIST = SQRT(VEC(1)*VEC(1)+VEC(2)*VEC(2)+VEC(3)*VEC(3))
                IF (DIST.GE.1.0D-15) DISTI = ONE/DIST
                CALL DFTB_GAM23D(DIST,DISTI,GAMMAV,GAMI,GAMJ,
     *                           TMPD(1,ISP,JSP),TMPD3(1,ISP,JSP),
     *                           TMPD3(1,JSP,ISP),DAMPXH,DAMPXHE,
     *                           SPE(ISP),SPE(JSP),HUBDER(ISP),
     *                           HUBDER(JSP),DFTB3)
                DO K = 1, 3
                  GDV2(K) = GAMMAV*VEC(K)*DISTI
                  IF (DFTB3) THEN
                    GDV3I(K) = GAMI*VEC(K)*DISTI
                    GDV3J(K) = GAMJ*VEC(K)*DISTI
                  END IF
                END DO
              END IF
C             GAMMAV = GAMD_FMO2(DIST,DISTI,TMPD(1,ISP,JSP),DAMPXH,
C    *                           DAMPXHE,SPE(ISP),SPE(JSP))
C             VAL = SCAL*(CHAMUL1*CHAMUL2*GAMMAV+VALTMP)*DISTI
C             FMODE(1,IAG) = FMODE(1,IAG) - VEC(1)*VAL
C             FMODE(2,IAG) = FMODE(2,IAG) - VEC(2)*VAL
C             FMODE(3,IAG) = FMODE(3,IAG) - VEC(3)*VAL
C             FMODE(1,JAG) = FMODE(1,JAG) + VEC(1)*VAL
C             FMODE(2,JAG) = FMODE(2,JAG) + VEC(2)*VAL
C             FMODE(3,JAG) = FMODE(3,JAG) + VEC(3)*VAL
              DO K = 1, 3
                VAL = CHAMUL1*CHAMUL2*GDV2(K)
                IF (DFTB3) THEN
                  VAL = VAL
     *              + ONE_THIRD*(CHAMULSQ*GDV3I(K)
     *                          +CHAMUL1*CHAMUL2*GDV3J(K))*CHAMUL2
                END IF
                FMODE(K,IAG) = FMODE(K,IAG) + VAL
                FMODE(K,JAG) = FMODE(K,JAG) - VAL
              END DO
            END DO
          END DO
        END DO
      ELSE IF (IAND(MODGAMMA,2).NE.0) THEN
        !! THIS IS NOT VERY EFFICIENT
        DO LFG = 1, NFG
          NATI = NATFRG(LFG)
          INDI = INDFRG(LFG)
          DO J = 1, NAT
            JAG = IAGLOB(J)
            JSP = ISPEFMO(JAG)
            NSEQ = INDI + NATFMOB*(INDJ+J-2)
            CHAMUL2 = POPMAT(J) - ZREFFMO(J)
            SHIFTCT(J) = SHIFTCT(J)
     *          + DDOT(NATI,CTMUL(INDI),1,GAMMA(NSEQ),1)
            IF (DFTB3) THEN
              JJ = INDJ+J-1
              SHIFTCT(J) = SHIFTCT(J) + ONE_THIRD*(
C    *          DDOT(NATI,CTMUL(INDI),1,GAMMA3(2,INDI,JJ),2)*TWO*CHAMUL2
C    *          + DDOT(NATI,CTMUL(NATFMOB+INDI),1,GAMMA3(1,INDI,JJ),2))
     *          DDOT(NATI,CTMUL(INDI),1,GAMMA3(INDI,JJ,2),1)*TWO*CHAMUL2
     *          + DDOT(NATI,CTMUL(NATFMOB+INDI),1,GAMMA3(INDI,JJ,1),1))
            END IF
            DO IAT = 1, NATI
              IAG = IATFRG(INDI+IAT-1)
              IF (.NOT.MASWRK .OR. IAG.EQ.JAG.AND..NOT.PERIOD) CYCLE
              CHAMUL1 = CTMUL(INDI+IAT-1)
              IF (DFTB3) CHAMULSQ = CTMUL(NATFMOB+INDI+IAT-1)
              ISP = ISPEFMO(IAG)
              VEC(1) = FMOC(1,IAG) - FMOC(1,JAG)
              VEC(2) = FMOC(2,IAG) - FMOC(2,JAG)
              VEC(3) = FMOC(3,IAG) - FMOC(3,JAG)
              DAMPXH = DAMPXH0
              CALL DFTB_DMPCHK(ISP,JSP,DAMPXH)
!             IF (DFTB3) THEN
!              CALL DFTB_GAMD3(DIST,DISTI,GAMI,GAMJ,HUBB(ISP),HUBB(JSP),
!    *            HUBDER(ISP),HUBDER(JSP),TMPD3(1,ISP,JSP),
!    *            TMPD3(1,JSP,ISP),DAMPXH,DAMPXHE,SPE(ISP),SPE(JSP))
!               VALTMP = ONE_THIRD
!    *            *CHAMUL2*(CHAMULSQ*GAMI+CHAMUL1*CHAMUL2*GAMJ)
!             END IF
!             GAMMAV = GAMD_FMO2(DIST,DISTI,TMPD(1,ISP,JSP),DAMPXH,
!    *                           DAMPXHE,SPE(ISP),SPE(JSP))
!             VAL = SCAL*(CHAMUL1*CHAMUL2*GAMMAV+VALTMP)*DISTI
!             FMODE(1,IAG) = FMODE(1,IAG) - VEC(1)*VAL
!             FMODE(2,IAG) = FMODE(2,IAG) - VEC(2)*VAL
!             FMODE(3,IAG) = FMODE(3,IAG) - VEC(3)*VAL
!             FMODE(1,JAG) = FMODE(1,JAG) + VEC(1)*VAL
!             FMODE(2,JAG) = FMODE(2,JAG) + VEC(2)*VAL
!             FMODE(3,JAG) = FMODE(3,JAG) + VEC(3)*VAL
              IF (PERIOD) THEN
                CALL DFTB_GAMDER_PBC(IAG,JAG,FMOC(1,IAG),FMOC(1,JAG),
     *                               GDV2,GDV3I,GDV3J,GDST2,GDST3I,
     *                               GDST3J,TMPD(1,ISP,JSP),
     *                               TMPD3(1,ISP,JSP),TMPD3(1,JSP,ISP),
     *                               DAMPXH,ISP,JSP)
                !! Stress tensor
                VAL = CHAMUL1*CHAMUL2
                CALL DAXPY(9,VAL,GDST2,1,STRFMO,1)
                IF (DFTB3) THEN
                  VAL = ONE_THIRD*CHAMULSQ*CHAMUL2
                  CALL DAXPY(9,VAL,GDST3I,1,STRFMO,1)
                  VAL = ONE_THIRD*CHAMUL1*CHAMUL2*CHAMUL2
                  CALL DAXPY(9,VAL,GDST3J,1,STRFMO,1)
                END IF
                IF (IAG.EQ.JAG) CYCLE
              ELSE
                DIST = SQRT(VEC(1)*VEC(1)+VEC(2)*VEC(2)+VEC(3)*VEC(3))
                IF (DIST.GE.1.0D-15) DISTI = ONE/DIST
                CALL DFTB_GAM23D(DIST,DISTI,GAMMAV,GAMI,GAMJ,
     *                           TMPD(1,ISP,JSP),TMPD3(1,ISP,JSP),
     *                           TMPD3(1,JSP,ISP),DAMPXH,DAMPXHE,
     *                           SPE(ISP),SPE(JSP),HUBDER(ISP),
     *                           HUBDER(JSP),DFTB3)
                DO K = 1, 3
                  GDV2(K) = GAMMAV*VEC(K)*DISTI
                  IF (DFTB3) THEN
                    GDV3I(K) = GAMI*VEC(K)*DISTI
                    GDV3J(K) = GAMJ*VEC(K)*DISTI
                  END IF
                END DO
              END IF
C             GAMMAV = GAMD_FMO2(DIST,DISTI,TMPD(1,ISP,JSP),DAMPXH,
C    *                           DAMPXHE,SPE(ISP),SPE(JSP))
C             VAL = SCAL*(CHAMUL1*CHAMUL2*GAMMAV+VALTMP)*DISTI
C             FMODE(1,IAG) = FMODE(1,IAG) - VEC(1)*VAL
C             FMODE(2,IAG) = FMODE(2,IAG) - VEC(2)*VAL
C             FMODE(3,IAG) = FMODE(3,IAG) - VEC(3)*VAL
C             FMODE(1,JAG) = FMODE(1,JAG) + VEC(1)*VAL
C             FMODE(2,JAG) = FMODE(2,JAG) + VEC(2)*VAL
C             FMODE(3,JAG) = FMODE(3,JAG) + VEC(3)*VAL
              DO K = 1, 3
                VAL = CHAMUL1*CHAMUL2*GDV2(K)
                IF (DFTB3) THEN
                  VAL = VAL
     *              + ONE_THIRD*(CHAMULSQ*GDV3I(K)
     *                          +CHAMUL1*CHAMUL2*GDV3J(K))*CHAMUL2
                END IF
                FMODE(K,IAG) = FMODE(K,IAG) + VAL
                FMODE(K,JAG) = FMODE(K,JAG) - VAL
              END DO
            END DO
          END DO
        END DO
      END IF
C
      RETURN
C
      END SUBROUTINE DFTB_FMOGRAD_EMBED
C
C-----------------------------------------------------------------------
C*MODULE DFTBFO    *DECK DFTB_PCMPOT
C>
C>    @brief   FMO-DFTB/PCM energy and gradient
C>
C>    @details This does the same to PCMPOT
C>
C>    @author Yoshio Nishimoto
C>    - Feb, 2016- Subroutine written
C>
C>           --- INPUT ---
C>    @param QSE     ASCs
C>    @param XYZ     Coordinates of tesserae
C>    @param L2      (NUM*NUM+NUM)/2
C>    @param MODE    PCM shift is saved or not
C>           --- OUTPUT ---
C>    @param VPCM    Contributions of PCM into Fock in AO basis
C>
C
      SUBROUTINE DFTB_PCMPOT(VPCM,QSE,XYZ,L2,MODE)
      use mx_limits, only: mxatm
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      DIMENSION VPCM(L2),QSE(*),XYZ(MXTS,3)
      LOGICAL GOPARR,DSKWRK,MASWRK
C
      PARAMETER (MXSPE=10)
      PARAMETER (ZERO=0.0D+00)
C
      COMMON /DFTBPR/ ETEMP,DFTBDP(MXSPE*14),DAMPXHE,HUBDER(MXSPE),
     *                ZREF(MXATM),SPNCST(6,MXSPE),SPE(MXATM),NSPE,
     *                MAXANG(MXATM),ISPE(MXATM),IND(MXATM+1),IDFTBD,
     *                PARAMDIR
      COMMON /FMCOM / X(1)
      COMMON /FMODTB/ LGAMMA,LZREF,LISPE,LSHIFTG,LSHIFTCT,
     *                LGAMMAD,MODESD,MODGAMMA,LDCTMP,LGAMMA3,LGAMMA3D,
     *                LNBOND,LCTMUL,LESPDFTB,LESPASC,NDFTBRST,LDFTBRST
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
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /PCMDIM/ mxsp,mxts,mempcm1,mempcm2,NTS
C
      CALL VALFM(LOADFM)
      LS     = LOADFM + 1
      LSHIFT = LS     + L2
      LAST   = LSHIFT + NAT
      NEED   = LAST - LOADFM - 1
      CALL GETFM(NEED)
C     IF (IODA(12).EQ.L2) THEN
C       CALL DAREAD(IDAF,IODA,X(LS),L2,12,0)
C     ELSE
        CALL VCLR(X(LS),1,L2)
        CALL DFTB_HAM_OVER(X(LS),X(LS),.FALSE.,.TRUE.,.FALSE.)
C     END IF
C
      IF (MODE.NE.0.OR.IAND(MODGAMMA,8).EQ.0) THEN
        IF (GOPARR) CALL VCLR(X(LSHIFT),1,NAT)
C
        IPCOUNT = ME - 1
        DO NM = 1, NAT
          IF(GOPARR) THEN
            IPCOUNT = IPCOUNT + 1
            IF(MOD(IPCOUNT,NPROC).NE.0) CYCLE
          END IF
          CX = C(1,NM)
          CY = C(2,NM)
          CZ = C(3,NM)
          TMP = ZERO
          DO ITS = 1, NTS
            DIST = SQRT((XYZ(ITS,1)-CX)*(XYZ(ITS,1)-CX)
     *                + (XYZ(ITS,2)-CY)*(XYZ(ITS,2)-CY)
     *                + (XYZ(ITS,3)-CZ)*(XYZ(ITS,3)-CZ))
            TMP = TMP + QSE(ITS)/DIST
          END DO
          X(LSHIFT+NM-1) = -TMP
        END DO
      ELSE
        IF (MASWRK) THEN
          DO NM = 1, NAT
            iag = ixftch(x(liaglob),nm)
            x(LSHIFT+NM-1) = X(LESPASC+IAG-1)
          END DO
        ELSE
          CALL VCLR(X(LSHIFT),1,NAT)
        END IF
      END IF
C     DO ITS = 1, NTS
C       XPP = XYZ(ITS,1)
C       YPP = XYZ(ITS,2)
C       ZPP = XYZ(ITS,3)
C       DO NM = 1, NAT
C         DIST = SQRT((C(1,NM)-XPP)*(C(1,NM)-XPP)
C    *              + (C(2,NM)-YPP)*(C(2,NM)-YPP)
C    *              + (C(3,NM)-ZPP)*(C(3,NM)-ZPP))
C         X(LSHIFT+NM-1) = X(LSHIFT+NM-1) - QSE(ITS)/DIST
C       END DO
C     END DO
c     write (6,*) "SHIFT IN PCMPOT"
c     do i = 1, nat
c       write (6,'(i3,f20.10)') i,x(lshift+i-1)
c     end do
c     write (6,*) "QSE IN PCMPOT"
c     write (6,'(10f8.2)') (qse(its),its=1,nts)
C
      CALL DFTB_SHIFT_FOCK(X(LSHIFT),VPCM,X(LS),NAT,L2,NSPE,ISPE,IND,
     *                     MAXANG)
c     write (6,*) "over ",(X(LS+i-1),i=1,l2)
c     write (6,*) "fock ",(VPCM(i),i=1,l2)
c     CALL PRTRIL(VPCM,NUM)
C
      CALL RETFM(NEED)
C
      RETURN
C
      END
C
C-----------------------------------------------------------------------
C*MODULE DFTBFO    *DECK DFTB_SDER
C>
C>    @brief FMO-DFTB/PCM gradient
C>
C>    @details Calculate dSmn/da contributions during FMO-DFTB/PCM
C>             gradient calculation. Although FELECQ formally calculates
C>             both d(r_Al)^-1/da and dSmn/da terms (V^a terms) in HF/DFT
C>             calculations, only d(r_Al)^-1/da terms are calculated
C>             in DFTB, so this subroutine has to calculate dSmn/da,
C>             which comes from V^a, too. In full (w/o fragments) DFTB,
C>             dSmn/da is calculated as shift in DFTB_GRAD.
C>
C>    @author Yoshio Nishimoto
C>    - Feb, 2016- Subroutine written
C>
C>           --- INPUT ---
C>    @param EPS     Xmn?
C>    @param DSINT   Temporary array for dSmn/da
C>    @param DA      Density matrix
C>    @param POT     Potential... see below
C>    @param L1      NUM
C>    @param L2      (NUM*NUM+NUM)/2
C>
C
      SUBROUTINE DFTB_SDER(EPS,DSINT,DA,POT,L1,L2)
      use mx_limits, only: mxatm
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION EPS(L2),DSINT(L2),DA(L2),POT(L2)
      LOGICAL GOPARR,DSKWRK,MASWRK
C
      PARAMETER (MXSPE=10)
C
      COMMON /DFTBPR/ ETEMP,DFTBDP(MXSPE*14),DAMPXHE,HUBDER(MXSPE),
     *                ZREF(MXATM),SPNCST(6,MXSPE),SPE(MXATM),NSPE,
     *                MAXANG(MXATM),ISPE(MXATM),IND(MXATM+1),IDFTBD,
     *                PARAMDIR
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /GRAD  / DE(3,MXATM)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
C
      DIMENSION AUP(9,9),BUP(9,9),VEC(3),SHIFT(NAT)
      LOGICAL SKP
C
      NSEQ = 0
C     POT(L2) = 1/2 Smn \sum_l (1/r_Al + 1/r_Bl)*\bar{q}_l
C     where q = (\Delta \bar{q}_l + \bar{q}_l^Z)
      DO IAT = 1, NAT
        DO K = 1, MAXANG(ISPE(IAT))*MAXANG(ISPE(IAT))
          MU = IND(IAT) + K
          NSEQ = NSEQ + MU
          IF (K.EQ.1) SHIFT(IAT) = POT(NSEQ)
        END DO
      END DO
C
      CALL VCLR(POT,1,L2)
      IPCOUNT = ME - 1
      CALL VCLR(BUP,1,81)
      DO NI = 1, NAT
        IF(GOPARR) THEN
          IPCOUNT = IPCOUNT + 1
          IF(MOD(IPCOUNT,NPROC).NE.0) CYCLE
        END IF
        CALL VCLR(DSINT,1,L2)
        DO KI = 1, 3
          DO NJ = 1, NAT
            IF (NI.EQ.NJ) CYCLE
            VEC(1) = C(1,NJ)-C(1,NI)
            VEC(2) = C(2,NJ)-C(2,NI)
            VEC(3) = C(3,NJ)-C(3,NI)
            CALL DFTB_SKODE_HELPA(1,NI,NJ,KI,0,0,0,VEC,AUP,BUP,
     *                            .FALSE.,.TRUE.,SKP)
            IF (SKP) CYCLE
            DO M = 1, IND(NJ+1)-IND(NJ)
              MU = IND(NJ)+M
              DO N = 1, IND(NI+1)-IND(NI)
                NU = IND(NI)+N
                CALL DFTB_CNVSQ(MU,NU,NSEQ)
                DSINT(NSEQ) = BUP(N,M) !! SDERAO
              END DO
            END DO
          END DO
C         HERE CALCULATES S_ij*W_ij^\Delta TERMS
C         \sum_ij Sij*W_ij^\Delta = \sum_mn dSmn/da*Xmn
C         Xmn = -1/2 \sum_op Dmo Wop^\Delta Dpn
          DE(KI,NI) = DE(KI,NI) + TRACEP(EPS,DSINT,L1)
C         HERE CALCULATES Dmn*dSmn/da TERMS
          CALL SHIFT_TRIANGLE_DER(NAT,NI,L2,IND,SHIFT,DSINT,POT)
          DE(KI,NI) = DE(KI,NI) + TRACEP(DA,POT,L1)
        END DO
        CALL SHIFT_TRIANGLE_DER_CLEAR(NAT,NI,L2,IND,POT)
      END DO
C
      END SUBROUTINE DFTB_SDER
C
C-----------------------------------------------------------------------
C*MODULE DFTBFO    *DECK DFTB_FMOPRC
C>
C>    @brief   FMO-DFTB cube
C>
C>    @details This does FMOPRC for DFTB.
C>
C>    @author Yoshio Nishimoto
C>    - Feb, 2016- Subroutine written
C>
C>           --- INPUT ---
C>    @param FDA     Scaling factor
C>    @param IBLOCK
C>    @param L1      NUM
C>    @param L2      (NUM*NUM+NUM)/2
C>    @param DA      Density matrix
C>    @param MINX    Min index for x coordinate
C>    @param NXG     Max index for x coordinate
C>    @param MINY    Min index for y coordinate
C>    @param NYG     Max index for y coordinate
C>    @param MINZ    Min index for z coordinate
C>    @param NZG     Max index for z coordinate
C>    @param IXMIN   Min index for x to be calculated
C>    @param IXMAX   Max index for x to be calculated
C>    @param IYMIN   Min index for y to be calculated
C>    @param IYMAX   Max index for y to be calculated
C>    @param IZMIN   Min index for z to be calculated
C>    @param IZMAX   Max index for z to be calculated
C>    @param JXMIN   Min range for x
C>    @param JXMAX   Max range for x
C>    @param JYMIN   Min range for y
C>    @param JYMAX   Max range for y
C>    @param JZMIN   Min range for z
C>    @param JZMAX   Max range for z
C>    @param ORIGIN  Origin of grid
C>    @param UX      Grid length along x
C>    @param UY      Grid length along y
C>    @param UZ      Grid length along z
C>    @param NSETGRID 1 or 2
C>    @param GRIDDISTR Grid is distributed or not
C>    @param LDGRIDL Memory location of grid
C>    @param NSHELL  Number of shells
C>           --- IN/OUTPUT ---
C>    @param DGRID   Density (probability?) on grids
C>
      SUBROUTINE DFTB_FMOPRC(FDA,IBLOCK,L1,L2,DA,DGRID,MINX,NXG,MINY,
     *                       NYG,MINZ,NZG,IXMIN,IXMAX,IYMIN,IYMAX,IZMIN,
     *                       IZMAX,JXMIN,JXMAX,JYMIN,JYMAX,JZMIN,JZMAX,
     *                       ORIGIN,UX,UY,UZ,NSETGRID,GRIDDISTR,
     *                       LDGRIDL,NSHELL)
      use mx_limits, only: mxatm
      USE DFTBPB_MOD,ONLY: PERIOD
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL GOPARR,DSKWRK,MASWRK,GOPSAV
C
      COMMON /ELPROP/ ELDLOC,ELMLOC,ELPLOC,ELFLOC,
     *                IEDEN,IEMOM,IEPOT,IEFLD,MODENS,
     *                IEDOUT,IEMOUT,IEPOUT,IEFOUT,
     *                IEDINT,IEMINT,IEPINT,IEFINT
      COMMON /FMCOM / X(1)
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /XYZPRP/ XP,YP,ZP
     *               ,DMX,DMY,DMZ
     *               ,QXX,QYY,QZZ,QXY,QXZ,QYZ
     *               ,QMXX,QMYY,QMZZ,QMXY,QMXZ,QMYZ
     *               ,OXXX,OXXY,OXXZ,OXYY,OYYY,OYYZ
     *               ,OXZZ,OYZZ,OZZZ,OXYZ
     *               ,OMXXX,OMXXY,OMXXZ,OMXYY,OMYYY
     *               ,OMYYZ,OMXZZ,OMYZZ,OMZZZ,OMXYZ
C
      DIMENSION DA(L2),DGRID(MINZ:NZG,MINY:NYG,MINX:NXG,NSETGRID),
     *          ORIGIN(3),UX(3),UY(3),UZ(3)
      LOGICAL GRIDDISTR
C
      DATA ELDEN  /8HELDEN   /
C
      MODENS=0
      CALL VALFM (LOADFM)
      LWRK = LOADFM + 1
      LAST = LWRK   + L2
      IF (NSETGRID.GT.1) THEN
        LMULA = LAST
        LMUL  = LMULA + NAT
        LMULS = LMUL  + L1
        LS    = LMULS + NSHELL
        LAST  = LS    + L2
      END IF
      NEED = LAST - LOADFM - 1
      CALL GETFM(NEED)
      IF (NSETGRID.GT.1) THEN
        !! for noprfg=4
C       IF (JCURFG.EQ.0) THEN
          CALL VCLR(X(LS),1,L2)
          CALL DFTB_HAM_OVER(X(LS),X(LS),.FALSE.,.TRUE.,.FALSE.)
          !! use original density matrix
          LOOP=0
          DO I=1,L1
            DO J=1,I-1
              LOOP=LOOP+1 
              DA(LOOP)=DA(LOOP)*0.5D+00
            ENDDO
C           SKIP THE DIAGONAL
            LOOP=LOOP+1
          ENDDO
          CALL DFTB_MULLIKEN(L1,L2,DA,X(LS),X(LMUL))
          LOOP=0
          DO I=1,L1
            DO J=1,I-1
              LOOP=LOOP+1
              DA(LOOP)=DA(LOOP)*2.0D+00
            ENDDO
C           SKIP THE DIAGONAL
            LOOP=LOOP+1
          ENDDO
          CALL DFTB_MULSA(X(LMUL),X(LMULS),X(LMULA),NUM,NSHELL,NAT)
C       ELSE
C         CALL DAREAD(IDAF,IODA,X(LMULA),NAT,556,0)
          CALL DAXPY(NAT,-1.0D+00,ZAN,1,X(LMULA),1)
C       END IF
        CALL DSCAL(NAT,FDA,X(LMULA),1)
      END IF
C
      NGLOCX=IXMAX-IXMIN+1
      NGLOCY=IYMAX-IYMIN+1
      NGLOCZ=IZMAX-IZMIN+1
      IF (GOPARR) FDA = FDA/NPROC
      GOPSAV = GOPARR
      GOPARR = .FALSE.
C
C     For PBC, wrap grid to home in the elementary cell
c     (grid spilling occurs because of having to account for adjacent cells).
C     For now, griddistr does not work with PBC.
C
      DO MX = IXMIN, IXMAX
        MX0=MX-1
        MX1=MX-IXMIN+1
C       IF(TESTYZ) THEN
C         XP = ORIGIN(1) + MX0*UX(1)
C         DUM=DUM0 + AA*((XP-AX)**2)
C         IF(DUM .GT. TOL) GO TO 490
C       ENDIF
        MXG=MX
c       IF(PERIOD) THEN
c         MXG=MOD(MX,NXG)
c         IF(MXG.LT.0) MXG=MXG+NXG
c         IF(MXG.EQ.0) MXG=NXG
c       ENDIF
        DO MY = IYMIN, IYMAX
          MY0=MY-1 
          MY1=MY-IYMIN+1 
C         IF(TESTZ2) THEN
C           XP = ORIGIN(1) + MX0*UX(1) + MY0*UY(1)
C           YP = ORIGIN(2) + MX0*UX(2) + MY0*UY(2)
C           DUM=DUM0 + AA*((XP-AX)**2 + (YP-AY)**2)
C           IF(DUM .GT. TOL) GO TO 480
C         ENDIF
          MYG=MY
c         IF(PERIOD) THEN
c           MYG=MOD(MY,NYG)
c           IF(MYG.LT.0) MYG=MYG+NYG
c           IF(MYG.EQ.0) MYG=NYG
c         ENDIF
          DO MZ = IZMIN, IZMAX
            IF(IBLOCK.GT.1 .AND. 
     *             MX.GE.JXMIN.AND.MX.LE.JXMAX. AND.
     *             MY.GE.JYMIN.AND.MY.LE.JYMAX. AND.
     *             MZ.GE.JZMIN.AND.MZ.LE.JZMAX) CYCLE
            MZ0=MZ-1
            MZ1=MZ-IZMIN+1
            MZG=MZ
c           IF(PERIOD) THEN
c             MZG=MOD(MZ,NZG)
c             IF(MZG.LT.0) MZG=MZG+NZG
c             IF(MZG.EQ.0) MZG=NZG
c           ENDIF
            XP = ORIGIN(1) + MX0*UX(1) + MY0*UY(1) + MZ0*UZ(1)
            YP = ORIGIN(2) + MX0*UX(2) + MY0*UY(2) + MZ0*UZ(2)
            ZP = ORIGIN(3) + MX0*UX(3) + MY0*UY(3) + MZ0*UZ(3)
C
C           DUM : ELECTRON DENSITY
C           DUM2: MOLECULAR ELECTROSTATIC POTENTIAL
C
            CALL DFTB_PRCALC(ELDEN,X(LWRK),L2)
            DUM = FDA*DDOT(L2,X(LWRK),1,DA,1)
            IF (NSETGRID.GT.1) THEN
              DUM2 = 0.0D+00
              DO IAT = 1, NAT
                DIST = SQRT((C(1,IAT)-XP)**2+(C(2,IAT)-YP)**2
     *                     +(C(3,IAT)-ZP)**2)
                IF (DIST.LE.1.0D-08) CYCLE
                !! fda is multiplied above
                DUM2 = DUM2 + X(LMULA+IAT-1)/DIST
              END DO
            END IF
C
            IF(GRIDDISTR) THEN
              LDGRIDLI=LDGRIDL+MZ1-1+(MY1-1+(MX1-1)*NGLOCY)*NGLOCZ
              X(LDGRIDLI)=X(LDGRIDLI)+DUM
              IF(NSETGRID.GT.1) THEN
                LDGRIDLI=LDGRIDLI+NGLOCZ*NGLOCY*NGLOCX
                X(LDGRIDLI)=X(LDGRIDLI)-DUM2
              ENDIF
            ELSE
c             write(6,7777) mx,my,mz,mxg,myg,mzg
c7777         format(1x,'wwwinds',3I3,' to',3I3)
              DGRID(MZG,MYG,MXG,1)=DGRID(MZG,MYG,MXG,1)+DUM
              IF(NSETGRID.GT.1)
     *         DGRID(MZG,MYG,MXG,2)=DGRID(MZG,MYG,MXG,2)-DUM2
            ENDIF
          END DO
        END DO
      END DO
C
      CALL RETFM(NEED)
      GOPARR = GOPSAV
C
      RETURN
C
      END SUBROUTINE DFTB_FMOPRC
C
C-----------------------------------------------------------------------
C*MODULE DFTBFO    *DECK DFTB_ESP2
C>
C>    @brief   FMO-DFTB electrostatic potential
C>
C>    @details Compute electrostatic potential in FMO-DFTB using GDDI,
C>             store it in memory in monomer iteration, and reuse it in
C>             monomer, dimer, and trimer SCF calculations.
C>
C>    @author Yoshio Nishimoto
C>    - Feb, 2018- Subroutine written
C>
C>           --- INPUT ---
C>    @param IDAF    10
C>    @param IODA    Length of unit
C>           --- OUTPUT ---
C>    @param SHIFT   Computed ESP in the current monomer SCC cycle.
C>
C
      SUBROUTINE DFTB_ESP2(IDAF,IODA,SHIFT,CONV)
      use mx_limits, only: mxatm
      USE DFTBPB_MOD,ONLY: PERIOD,PME,DFTB_GAMMA_PBC,DFTB_ESP_PME
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DDI_WORLD,DDI_GROUP
C
      PARAMETER (MXSPE=10,ONE=1.0D+00,THREE=3.0D+00)
      PARAMETER (DDI_WORLD=0,DDI_GROUP=1)
C
      LOGICAL DFTBFL,SCC,SRSCC,DFTB3,DAMPXH,LCDFTB
      LOGICAL GOPARR,DSKWRK,MASWRK
      LOGICAL MYJOB,ISGDDI,PAROUT,INITGDDI,WASGDDI,MLGDDI
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
      COMMON /FMODTB/ LGAMMA,LZREF,LISPE,LSHIFTG,LSHIFTCT,
     *                LGAMMAD,MODESD,MODGAMMA,LDCTMP,LGAMMA3,LGAMMA3D,
     *                LNBOND,LCTMUL,LESPDFTB,LESPASC,NDFTBRST,LDFTBRST
      common /fmodim/ maxbnd,maxknd,maxcbs,maxcao,maxvec,maxl1,maxnat,
     *                maxabd,maxbas,maxbbd,maxlmo,maxslo,maxabd2,maxrij
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
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
C
      DIMENSION IODA(*),SHIFT(*)
      DIMENSION SHIFT3(MAXNAT*2),CHAMUL(MAXNAT)
      LOGICAL   CONV
C
      NATFMOB = NATFMO + NBDFG
      IF (DFTB3) THEN
        NDIM = NATFMOB*3
      ELSE
        NDIM = NATFMOB
      END IF
      CALL DCOPY(NDIM,0.0D+00,0,X(LESPDFTB),1)
      IF (ISGDDI) THEN
        CALL GDDI_SCOPE(DDI_GROUP)
        CALL GDDICOUNT(-1,LGROUP,MYJOB)
      END IF
C
      LPOPMATI=LPOPMAT+(MAXNAT*NFG)*(ICURPOP-1)
C
      IF (.NOT.CONV.AND.IAND(MODGAMMA,3).EQ.0) THEN
        CALL DFTB_ESP3(NATFMOB,X(LNATFRG),X(LINDFRG),X(LIATFRG),
     *                 X(LPOPMATI),X(LISPE),X(LFMOC),X(LZREF),
     *                 X(LGAMMA),X(LGAMMA3),X(LESPDFTB))
      ELSE IF (.NOT.CONV.AND.IAND(MODGAMMA,2).EQ.2) THEN
        IF (ISGDDI) CALL GDDI_SCOPE(DDI_WORLD)
        NVCORE = NATFMOB/NPROC
        IF (MOD(NATFMOB,NPROC).NE.0) NVCORE = NVCORE + 1
        INDL = 1 + ME*NVCORE
        INDH = (ME+1)*NVCORE
        IF (INDH.GT.NATFMOB) INDH = NATFMOB
        NVEC = INDH-INDL+1
c       NVEC can be negative!
C       
        if(nvec.gt.0)
     *  CALL DGEMV('T',NATFMOB,NVEC,
     *             1.0D+00,X(LGAMMA+NATFMOB*(INDL-1)),NATFMOB,
     *                     X(LSHIFTG),1,
     *             0.0D+00,X(LESPDFTB+INDL-1),1)
        IF (DFTB3) THEN
        if(nvec.gt.0)
     *    CALL DGEMV('T',NATFMOB,NVEC,
     *               1.0D+00,X(LGAMMA3+NATFMOB**2+NATFMOB*(INDL-1)),
     *                       NATFMOB,X(LSHIFTG),1,
     *               0.0D+00,X(LESPDFTB+NATFMOB+INDL-1),1)
C         CALL DGEMV('T',NATFMOB,NVEC,
C    *               1.0D+00/3.0D+00,X(LGAMMA3+NATFMOB*(INDL-1)),
C    *                               NATFMOB,X(LSHIFTG+NATFMOB),1,
C    *               1.0D+00,X(LESPDFTB+INDL-1),1)
        if(nvec.gt.0)
     *    CALL DGEMV('T',NATFMOB,NVEC,
     *               1.0D+00,X(LGAMMA3+NATFMOB*(INDL-1)),
     *                       NATFMOB,X(LSHIFTG+NATFMOB),1,
     *               1.0D+00,X(LESPDFTB+NATFMOB*2+INDL-1),1)
          CALL DAXPY(NATFMOB,1.0D+00/3.0D+00,X(LESPDFTB+NATFMOB*2),1,
     *               X(LESPDFTB),1)
        END IF
        IF (GOPARR.AND..NOT.ISGDDI)
     *    CALL DDI_GSUMF(2423,X(LESPDFTB),NDIM)
      ELSE
        NN = 1
        NN2 = NATFMO
        IF (IAND(MODGAMMA,2).NE.0) THEN
          NN = NN + 1
          NN2 = NN2 + NBDFG
        END IF
C
        CALL VALFM(LOADFM)
        LDIST = LOADFM + 1
        LGAM2 = LDIST + MAXNAT*(MAXNAT+1)/2
        LGAM3 = LGAM2 + MAXNAT*(MAXNAT+1)/2
        LAST  = LGAM3 + MAXNAT*MAXNAT
        NEED  = LAST - LOADFM - 1
        CALL GETFM(NEED)
C
        DO LFG = 1, NFG
          IF (ISGDDI) THEN
            CALL GDDICOUNT(0,LGROUP,MYJOB)
            IF (.NOT.MYJOB) CYCLE
          END IF
          ICURFG = LFG
C    
          NAT  = IXFTCH(X(LNATFRG),LFG)
          INDI = IXFTCH(X(LINDFRG),LFG)
          DO I = 1, NAT
            IAG = IXFTCH(X(LIATFRG),INDI+I-1)
            CALL IXSTOR(X(LIAGLOB),I,IAG)
            CALL DCOPY(3,X(LFMOC+3*(IAG-1)),1,C(1,I),1)
            ISPE(I) = IXFTCH(X(LISPE),IAG)
          END DO
          IODA(566) = 0
          IODA(567) = 0
C    
          modgamsav = modgamma
          if (iand(modgamma,8).eq.8.and.iand(modgamma,2).eq.0) then
            modgamma = modgamma - 8
          end if
          CALL DFTB_ESP(X(LIAGLOB),X(LNATFRG),X(LINDFRG),X(LIATFRG),
     *                  X(LPOPMATI),SHIFT,SHIFT,X(LISPE),X(LFMOC),
     *                  X(LGAMMA3),X(LINDAT),X(LIALOC),NN,NN2,1)
          modgamma = modgamsav
C    
          IF (PERIOD) THEN
            IF (SRSCC) THEN
              NNN = 0
              DO I = 1, NSPE
                NNN = NNN + MAXANG(I)
              END DO
              NSHELL=1 !! needs the correct value
            ELSE
              NNN = NSPE
              NSHELL=nat !! NSPE
            END IF
            CALL VALFM(LOADFM)
            LTMPGAM  = LOADFM + 1
            LAST     = LTMPGAM + 11*NNN*NNN
            NEEDGAM  = LAST  - LOADFM - 1
            CALL GETFM(NEEDGAM)
            CALL DFTB_GAMMA_PBC(0,0,0,NAT,NSPE,NSHELL,NNN,ISPE,MAXANG,
     *                          C,SPE,HUBBL,DAMPXHE,HUBDER,
     *                          DAMPXH,SRSCC,DFTB3,X(LGAM2),
     *                          X(LGAM3),X(LTMPGAM),X(LAST),X(LAST),
     *                          X(LAST),X(LAST),X(LAST),.FALSE.,
     *                          DUM,X(LAST))
            CALL RETFM(NEEDGAM)
          ELSE
            NSEQ = 0
            DO I = 1, NAT
              DO J = 1, I
                DIST = SQRT((C(1,I)-C(1,J))**2
     *                     +(C(2,I)-C(2,J))**2
     *                     +(C(3,I)-C(3,J))**2)
                NSEQ = NSEQ + 1
                X(LDIST+NSEQ-1) = DIST
              END DO
            END DO
            CALL DFTB_GAMMA(DAMPXHE,HUBBL,X(LGAM2),X(LDIST),NAT,NSPE,
     *                      ISPE,MAXANG,SPE,DAMPXH,SRSCC)
          END IF
C    
          CALL DCOPY(NAT,X(LPOPMATI+MAXNAT*(LFG-1)),1,CHAMUL,1)
          DO I = 1, NAT
            CHAMUL(I) = CHAMUL(I) - X(LZREF+INDI+I-2)
          END DO
          CALL DAREAD(IDAF,IODA,SHIFT,NAT,566,0)
          CALL DFTB_2ND_SHIFT(SHIFT,CHAMUL,X(LGAM2),NAT)
C    
          IF (DFTB3) THEN
            CALL DAREAD(IDAF,IODA,SHIFT3,NAT*2,567,0)
            !! For SCZV
            CALL DCOPY(NAT,SHIFT3,1,X(LESPDFTB+2*NATFMOB+INDI-1),1)
            IF (.NOT.PERIOD) THEN
              CALL VCLR(X(LGAM3),1,NAT*NAT)
              CALL DFTB_GAMMA3(NAT,HUBBL,X(LGAM3),X(LDIST))
            END IF
            DO I = 1, NAT
              SHIFT(I) = SHIFT(I) + SHIFT3(NAT+I)
              SHIFT3(I) = SHIFT3(I)*THREE
              DO J = 1, NAT
                GAMI = X(LGAM3+I-1+NAT*(J-1))
                GAMJ = X(LGAM3+J-1+NAT*(I-1))
                SHIFT(I) = SHIFT(I) + CHAMUL(J)*CHAMUL(J)*GAMJ/THREE
                SHIFT3(I) = SHIFT3(I) + CHAMUL(J)*GAMI
              END DO
            END DO
C    
            CALL DCOPY(NAT,SHIFT3,1,X(LESPDFTB+NATFMOB+INDI-1),1)
          END IF
          CALL DCOPY(NAT,SHIFT,1,X(LESPDFTB+INDI-1),1)
        END DO
C
        CALL RETFM(NEED)
      END IF
C
      IF (PERIOD.AND.PME) THEN
        CALL DFTB_ESP_PME(0,0,0,0,X(LPOPMATI),X(LESPDFTB))
      END IF
C     write (*,*) "aaa"
C     do i = 1, natfmob
C       write (*,'(i3,f20.10)') i,x(lespdftb+i-1)
C     end do
C
      IF (ISGDDI) THEN
        CALL GDDICOUNT(1,LGROUP,MYJOB)
C       IF (GOPARR.AND.IAND(MODGAMMA,2).NE.2)
C    *    CALL DSCAL(NDIM,ONE/NPROC,X(LESPDFTB),1)
        IF (GOPARR)  CALL DSCAL(NDIM,ONE/NPROC,X(LESPDFTB),1)
        CALL GDDI_SCOPE(DDI_WORLD)
        CALL DDI_GSUMF(2423,X(LESPDFTB),NDIM)
C       CALL DDI_BCAST(2423,'F',X(LESPDFTB),NDIM,MASTER)
        CALL GDDI_SCOPE(DDI_GROUP)
      END IF
C     if (maswrk) then
C     write (*,*) "shiftesp"
C       do i = 1, ndim
C         write (*,'(i3,f20.10)') i,x(lespdftb+i-1)
C       end do
C     end if
C 
      RETURN
C
      END SUBROUTINE DFTB_ESP2
C
C-----------------------------------------------------------------------
C*MODULE DFTBFO    *DECK DFTB_PCMPOT2
C>
C>    @brief   FMO-DFTB/PCM ASC contribution
C>
C>    @details Compute ASC contributions in FMO-DFTB/PCM using DDI,
C>             store it in memory in monomer iteration, and reuse it in
C>             monomer, dimer, and trimer SCF calculations.
C>
C>    @author Yoshio Nishimoto
C>    - Feb, 2018- Subroutine written
C>
C>           --- INPUT ---
C>    @param QSE     Induced ASCs
C>    @param XYZ     Coordinates of tesserae
C>
C
      SUBROUTINE DFTB_PCMPOT2(QSE,XYZ)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION QSE(NTS),XYZ(MXTS,3)
C
      LOGICAL GOPARR,DSKWRK,MASWRK
C
      PARAMETER (ZERO=0.0D+00)
C
      COMMON /FMCOM / X(1)
      COMMON /FMODTB/ LGAMMA,LZREF,LISPE,LSHIFTG,LSHIFTCT,
     *                LGAMMAD,MODESD,MODGAMMA,LDCTMP,LGAMMA3,LGAMMA3D,
     *                LNBOND,LCTMUL,LESPDFTB,LESPASC,NDFTBRST,LDFTBRST
      COMMON /FMOINF/ NFG,NLAYER,NATFMO,NBDFG,NAOTYP,NBODY,NSEGM
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
      COMMON /PCMDIM/ mxsp,mxts,mempcm1,mempcm2,NTS
C
      IF (IAND(MODGAMMA,8).EQ.0) RETURN
C
      IF (GOPARR) CALL VCLR(X(LESPASC),1,NATFMO)
C
      IPCOUNT = ME - 1
      DO NM = 1, NATFMO
        IF(GOPARR) THEN
          IPCOUNT = IPCOUNT + 1
          IF(MOD(IPCOUNT,NPROC).NE.0) CYCLE
        END IF
        CX = X(LFMOC+3*(NM-1))
        CY = X(LFMOC+3*(NM-1)+1)
        CZ = X(LFMOC+3*(NM-1)+2)
        TMP = ZERO
        DO ITS = 1, NTS
          DIST = SQRT((XYZ(ITS,1)-CX)*(XYZ(ITS,1)-CX)
     *              + (XYZ(ITS,2)-CY)*(XYZ(ITS,2)-CY)
     *              + (XYZ(ITS,3)-CZ)*(XYZ(ITS,3)-CZ))
          TMP = TMP + QSE(ITS)/DIST
        END DO
        X(LESPASC+NM-1) = -TMP
      END DO
C
      IF (GOPARR) CALL DDI_GSUMF(914,X(LESPASC),NATFMO)
c     write(6,*) 'wwwset',x(lespasc)
C
      RETURN
C
      END SUBROUTINE DFTB_PCMPOT2
C
C-----------------------------------------------------------------------
C*MODULE DFTBFO    *DECK DFTB_GAM23D
C>
C>    @brief   Gamma function for DFTB
C>
C>    @details Compute Gamma function for DFTB.
C>
C>    @author  Yoshio Nishimoto
C>
C
      SUBROUTINE DFTB_GAM23D(DIST,DISTI,GAM2,GAM3I,GAM3J,VALUE2,VALUE3A,
     *  VALUE3B,DAMPXH,DAMPXHE,SPE1,SPE2,HUBDERI,HUBDERJ,DFTB3)
C
      IMPLICIT NONE
C     Declarations of arguments
      DOUBLE PRECISION DIST,DISTI,GAM2,GAM3I,GAM3J,VALUE2(10),
     *                 VALUE3A(11)
      DOUBLE PRECISION VALUE3B(11),DAMPXHE,SPE1,SPE2,HUBDERI,HUBDERJ
      LOGICAL DAMPXH,DFTB3
C     Other declarations
      DOUBLE PRECISION, PARAMETER :: ZERO=0.0D+00,DISTTOL=1.0D-08,
     *  HUBTOL=1.0D-04,HALF=0.5D+00,ONE=1.0D+00,TWO=2.0D+00
C
      DOUBLE PRECISION H/8HH       /
      DOUBLE PRECISION d2fab_dadr,d2fab_dbdr,d2fba_dadr,d2fba_dbdr
      DOUBLE PRECISION d2h_dudr,dfab_da,dfab_db,dfab_dr,dfba_da,dfba_db
      DOUBLE PRECISION dfba_dr,dga_dr,dh_dr,dh_du,dist2,dist3,disti2
      DOUBLE PRECISION ds_dr,ds_drba,ds_duab,ds_duba,exp_ur,expa,expb
      DOUBLE PRECISION fab,fba,ga,gamtmp3,gamtmp3ab,gamtmp3ba,hubmean
      DOUBLE PRECISION hubmeanexp,hubmeanexp1,s,taua,taub,taumean,tmp
C
      IF (DIST.LE.DISTTOL) THEN
        GAM2  = 0.0D+00
        GAM3I = 0.0D+00
        GAM3J = 0.0D+00
        RETURN
      ELSE
        TAUA = VALUE2(1)
        TAUB = VALUE2(2)
        IF (DAMPXH.AND.(SPE1.EQ.H.OR.SPE2.EQ.H)) THEN
          IF (ABS(TAUA-TAUB).LE.HUBTOL) THEN
            TAUMEAN = VALUE2(3)
            HUBMEAN = TAUA*3.125D-01
            HUBMEANEXP = HUBMEAN**DAMPXHE
            DIST2 = DIST*DIST
            EXP_UR = EXP(-TAUMEAN*DIST-HUBMEANEXP*DIST2)
            DISTI2 = DISTI*DISTI
            GA = DISTI+VALUE2(4)+VALUE2(5)*DIST+VALUE2(6)*DIST2
            DGA_DR = -DISTI2+VALUE2(5)+TWO*VALUE2(6)*DIST
            DS_DR = EXP_UR*DGA_DR - TAUA*EXP_UR*GA
C
            GAM2 = -DISTI2 - DS_DR + EXP_UR*GA*HUBMEANEXP*TWO*DIST
            IF (DFTB3) THEN
              GA = DISTI+VALUE3A(3)+VALUE3A(1)*DIST+VALUE3A(4)*DIST2
              GAM3I =-EXP_UR*(VALUE3A(6) + VALUE3A(7)*DIST
     *              + VALUE3A(8)*DIST2)*DIST
C
              HUBMEANEXP = HUBMEANEXP*DIST2
              GAM3I = GAM3I + HUBMEANEXP
     *          *(TWO*EXP_UR*(-DISTI+VALUE3A(9)
     *              + VALUE3A(10)*DIST+VALUE3A(11)*DIST2)
     *          - DISTI*DAMPXHE/HUBMEAN
     *            *(-HALF*DIST*DS_DR + EXP_UR*GA*(HUBMEANEXP-ONE)))
              GAM3I = GAM3I*HUBDERI
              GAM3J = GAM3I
            END IF
          ELSE
            DIST2 = DIST*DIST
            HUBMEAN = (TAUA+TAUB)*1.5625D-01
            HUBMEANEXP = HUBMEAN**DAMPXHE
            EXPA = EXP(-TAUA*DIST-HUBMEANEXP*DIST2)
            EXPB = EXP(-TAUB*DIST-HUBMEANEXP*DIST2)
            FAB = VALUE2(9) - VALUE2(7)*DISTI
            FBA = VALUE2(10) - VALUE2(8)*DISTI
            DISTI2 = DISTI*DISTI
            DFAB_DR = VALUE2(7)*DISTI2
            DFBA_DR = VALUE2(8)*DISTI2
            DS_DR = EXPA*(DFAB_DR-TAUA*FAB)
     *            + EXPB*(DFBA_DR-TAUB*FBA)
            GAM2 = -DISTI2 - DS_DR
     *           + (EXPA*FAB + EXPB*FBA)*HUBMEANEXP*TWO*DIST
            IF (DFTB3) THEN
              FAB = VALUE3A(2)-VALUE3A(1)*DISTI
              DFAB_DA = VALUE3A(5) + VALUE3A(6)*DISTI
              DFBA_DA = VALUE3A(7) + VALUE3A(8)*DISTI
              D2FAB_DADR = -VALUE3A(6)*DISTI2
              GAMTMP3AB = EXPA*((TAUA*DIST-ONE)*FAB-TAUA*DFAB_DA
     *                  + D2FAB_DADR - DISTI*VALUE3A(1))
     *                  + EXPB*(-D2FAB_DADR - TAUB*DFBA_DA)
              DS_DUAB = EXPA*(DFAB_DA-DIST*FAB)+EXPB*DFBA_DA
              HUBMEANEXP1 = HUBMEANEXP*VALUE3A(9)
              DH_DR = -TWO*DIST*HUBMEANEXP
              GAM3I = -3.2D+00*(GAMTMP3AB+DS_DUAB*DH_DR)
C
              FBA = VALUE3A(4)-VALUE3A(3)*DISTI
              DFAB_DB = VALUE3B(7) + VALUE3B(8)*DISTI
              DFBA_DB = VALUE3B(5) + VALUE3B(6)*DISTI
              D2FBA_DBDR = -VALUE3B(6)*DISTI2
              GAMTMP3BA = EXPB*((TAUB*DIST-ONE)*FBA-TAUB*DFBA_DB
     *                  + D2FBA_DBDR - DISTI*VALUE3A(3))
     *                  + EXPA*(-D2FBA_DBDR - TAUA*DFAB_DB)
              DS_DUBA = EXPB*(DFBA_DB-DIST*FBA)+EXPA*DFAB_DB
              GAM3J = -3.2D+00*(GAMTMP3BA+DS_DUBA*DH_DR)
C
              DH_DU = -HALF*DIST2*DAMPXHE*HUBMEANEXP1
              D2H_DUDR = DAMPXHE*DIST*HUBMEANEXP1
     *                   *(DIST2*HUBMEANEXP-ONE)
              S = EXPA*FAB + EXPB*FBA
              TMP = DS_DR*DH_DU + S*D2H_DUDR
              GAM3I = (GAM3I-TMP)*HUBDERI
              GAM3J = (GAM3J-TMP)*HUBDERJ
            END IF
          END IF
        ELSE
          IF (ABS(TAUA-TAUB).LE.HUBTOL) THEN
            TAUMEAN = VALUE2(3)
            EXP_UR = EXP(-TAUMEAN*DIST)
            DIST2 = DIST*DIST
            DISTI2 = DISTI*DISTI
            GA = DISTI+VALUE2(4)+VALUE2(5)*DIST+VALUE2(6)*DIST2
            DGA_DR = -DISTI2+VALUE2(5)+TWO*VALUE2(6)*DIST
            DS_DR = EXP_UR*DGA_DR - TAUA*EXP_UR*GA
            GAM2 = -DISTI2 - DS_DR
            IF (DFTB3) THEN
              GA = DISTI+VALUE3A(3)+VALUE3A(1)*DIST+VALUE3A(4)*DIST2
              DGA_DR = -DISTI2+VALUE3A(1)+VALUE3A(2)*DIST
              DS_DR = EXP_UR*DGA_DR-TAUA*EXP_UR*GA
              DIST3 = DIST2*DIST
              GAMTMP3 = EXP_UR*(VALUE3A(6)*DIST + VALUE3A(7)*DIST2
     *                + VALUE3A(8)*DIST3)
              GAMTMP3 = -GAMTMP3
              GAM3I = GAMTMP3 * HUBDERI
              GAM3J = GAM3I
            END IF
          ELSE
            DISTI2 = DISTI*DISTI
            EXPA = EXP(-TAUA*DIST)
            EXPB = EXP(-TAUB*DIST)
            FAB = VALUE2(9) - VALUE2(7)*DISTI
            FBA = VALUE2(10) - VALUE2(8)*DISTI
            DFAB_DR = VALUE2(7)*DISTI2
            DFBA_DR = VALUE2(8)*DISTI2
            DS_DR = EXPA*DFAB_DR - TAUA*EXPA*FAB
     *            + EXPB*DFBA_DR - TAUB*EXPB*FBA
            GAM2 = -DISTI2 - DS_DR
            IF (DFTB3) THEN
              !! CALCULATE HUBA NE HUBB
              FAB = VALUE3A(2)-VALUE3A(1)*DISTI
              FBA = VALUE3A(4)-VALUE3A(3)*DISTI
              DFAB_DR = VALUE3A(1)*DISTI2
              DFBA_DR = VALUE3A(3)*DISTI2
              !! CALCULATE -dSf/dR
              DS_DR   = EXPA*DFAB_DR - TAUA*EXPA*FAB
     *                + EXPB*DFBA_DR - TAUB*EXPB*FBA
              DS_DRBA = EXPB*DFBA_DR - TAUB*EXPB*FBA
     *                + EXPA*DFAB_DR - TAUA*EXPA*FAB
              DFAB_DA = VALUE3A(5) + VALUE3A(6)*DISTI
              DFBA_DA = VALUE3A(7) + VALUE3A(8)*DISTI
              DFAB_DB = VALUE3B(7) + VALUE3B(8)*DISTI
              DFBA_DB = VALUE3B(5) + VALUE3B(6)*DISTI
              D2FAB_DADR = -VALUE3A(6)*DISTI2
              D2FBA_DADR = -D2FAB_DADR
              D2FBA_DBDR = -VALUE3B(6)*DISTI2
              D2FAB_DBDR = -D2FBA_DBDR
              GAMTMP3 = EXPA*((TAUA*DIST-ONE)*FAB-TAUA*DFAB_DA
     *          + D2FAB_DADR - DIST*DFAB_DR)
     *          + EXPB*(D2FBA_DADR - TAUB*DFBA_DA)
              GAMTMP3 = -3.2D+00*GAMTMP3
              GAMTMP3BA = EXPB*((TAUB*DIST-ONE)*FBA-TAUB*DFBA_DB
     *          + D2FBA_DBDR - DIST*DFBA_DR)
     *          + EXPA*(D2FAB_DBDR - TAUA*DFAB_DB)
              GAMTMP3BA = -3.2D+00*GAMTMP3BA
              GAM3I = GAMTMP3 * HUBDERI
              GAM3J = GAMTMP3BA * HUBDERJ
            END IF
          END IF
        END IF
      END IF
C
      RETURN
C
      END SUBROUTINE DFTB_GAM23D
C
C-----------------------------------------------------------------------
C*MODULE DFTBFO    *DECK dftb_addlag
C>
C>    @brief   Add Lagrangian for DFTB.
C>
C>    @details Compute Lagrangian for DFTB.
C>
C>    @author  Yoshio Nishimoto
C>
C
      subroutine dftb_addlag(mode,l1,dd,ss,wrk1,wrk2,wrk3,shift,shift2,
     *                       nat,c,fmode,iaglob)
c
      USE MX_LIMITS, ONLY: MXATM
      USE DFTBPB_MOD,ONLY: PERIOD,INDBND,MAXINT,PBCBOX,STRFMO,STRESS
C
      IMPLICIT NONE
C     Declarations of arguments
      DOUBLE PRECISION dd,ss,wrk1,wrk2,wrk3,shift,shift2
      INTEGER  mode,l1
      DOUBLE PRECISION c,fmode
      INTEGER nat,iaglob
C     Declarations of common blocks
      DOUBLE PRECISION ETEMP,DFTBDP,DAMPXHE,HUBDER
      DOUBLE PRECISION ZREF,SPNCST,SPE
      DOUBLE PRECISION PARAMDIR
      INTEGER NSPE
      INTEGER MAXANG,ISPE,IND,IDFTBD
      DOUBLE PRECISION REPCUT
      INTEGER NUMREP
      INTEGER LREPINTV,LREPSHORT,LREPCOEFF
      DOUBLE PRECISION SKDIM,SKSPIN,QREFL
      DOUBLE PRECISION HUBBL,QREF,HUBB,SKCUT2
      INTEGER NEEDSK,LSKHTAB,LSKSTAB
      INTEGER LSKGRID,LSKSELF
      INTEGER ME,MASTER,NPROC,IBTYP,IPTIM
      LOGICAL GOPARR,DSKWRK,MASWRK
C     Other declarations
      INTEGER, PARAMETER :: MXSPE=10
      DOUBLE PRECISION aum,bum
      INTEGER i,ig,ixtrv,iytrv,iztrv,j,jj,k,k2,l,l2,l3,m,nintpair,nnn
      DOUBLE PRECISION tmp,vec
C
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
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
c
      dimension dd(*),ss(*),wrk1(l1,l1),wrk2(*),wrk3(*),shift(*),
     *          shift2(*),c(3,*),fmode(3,*)
      dimension iaglob(*)
      dimension vec(3),aum(9,9),bum(9,9)
c
      logical   skp
      integer   o
c
      l2 = l1*(l1+1)/2
      l3 = l1*l1
c
      if (mode.eq.0) then
        call vclr(wrk2,1,l2)
        call dftb_shift_fock(shift,wrk2,ss,nat,l2,nspe,ispe,ind,maxang)
        call cpytsq(wrk2,wrk1,l1,1)
        call cpytsq(dd,wrk2,l1,1)
C
        call dgemm('T','N',l1,l1,l1,
     *             1.0d+00,wrk2,l1,wrk1,l1,
     *             0.0d+00,wrk3,l1)
        call dgemm('N','N',l1,l1,l1,
     *             1.0d+00,wrk3,l1,wrk2,l1,
     *             0.0d+00,wrk1,l1)
        call dscal(l3,-0.5d+00,wrk1,1)
C
        call vclr(wrk2,1,l2)
        call dftb_shift_fock(shift2,wrk2,dd,nat,l2,nspe,ispe,ind,maxang)
        call cpytsq(wrk2,wrk3,l1,1)
        call daxpy(l3,1.0d+00,wrk3,1,wrk1,1)
      else if (mode.eq.1) then
        call dcopy(l3,dd,1,wrk1,1)
      end if
C
      NNN = 4210
      CALL VCLR(BUM,1,81)
      DO I = 1, NAT
        IF (MODE.EQ.0) THEN
          IG = IAGLOB(I)
        ELSE IF (MODE.EQ.1) THEN
          IG = I
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
            CALL DFTB_SKODE_HELPA(1,I,J,K,0,0,0,VEC,AUM,BUM,
     *                            .FALSE.,.TRUE.,SKP)
            IF (SKP) CYCLE
            TMP = 0.0D+00
            DO L = 1, IND(J+1) - IND(J)
              DO M = 1, IND(I+1) - IND(I)
                TMP = TMP + 2.0D+00*BUM(M,L)*WRK1(IND(I)+M,IND(J)+L)
              END DO
            END DO
            IF ((MODE.EQ.0.AND.MASWRK) .OR. MODE.EQ.1)
     *        FMODE(K,IG) = FMODE(K,IG) + TMP
            IF (PERIOD) THEN
              IF (MODE.EQ.0) THEN
                if (maswrk) then
                DO K2 = 1, 3
                  STRFMO(K,K2) = STRFMO(K,K2) - TMP*VEC(K2)
                END DO
                end if
              ELSE IF (MODE.EQ.1) THEN
                DO K2 = 1, 3
                  STRESS(K,K2) = STRESS(K,K2) - TMP*VEC(K2)
                END DO
              END IF
            END IF
          END DO
        END DO
      END DO
C       if (maswrk.and.mode.eq.1) then
C         write (*,*) (stress(i,1),i=1,9)
C       end if
c
      return
c
      end subroutine dftb_addlag
C
C-----------------------------------------------------------------------
C*MODULE DFTBFO    *DECK DFTB_PBCPTR_FMO
C>
C>    @brief   Save pointers for DFTB.
C>
C>    @details Store array pointers.
C>
C>    @author  Yoshio Nishimoto
C>
C
      SUBROUTINE DFTB_PBCPTR_FMO(NFG,MAXNAT,NATFMO,NBDFG,NATFRG,INDFRG,
     *                           IATFRG,IAGLOB,INDAT,IALOC,ISPEFMO,
     *                           NBOND,
     *                           FMOC,ZREF)
      USE DFTBPB_MOD,ONLY:DFTB_PBCPTR_FMO2
      IMPLICIT NONE
C     Declarations of arguments
      INTEGER NFG,MAXNAT,NATFMO,NBDFG,NATFRG,INDFRG
      INTEGER IATFRG,IAGLOB,INDAT,IALOC,ISPEFMO
      INTEGER NBOND
      DOUBLE PRECISION FMOC,ZREF
      DIMENSION NATFRG(*),INDFRG(*),IATFRG(*),IAGLOB(*),INDAT(*),
     *          IALOC(*),ISPEFMO(*),NBOND(*),FMOC(*),ZREF(*)
      CALL DFTB_PBCPTR_FMO2(NFG,MAXNAT,NATFMO,NBDFG,NATFRG,INDFRG,
     *                      IATFRG,IAGLOB,INDAT,IALOC,ISPEFMO,NBOND,
     *                      FMOC,ZREF)
      END SUBROUTINE DFTB_PBCPTR_FMO
C
C-----------------------------------------------------------------------
C*MODULE DFTBFO    *DECK dftb_cnstchamul
C>
C>    @brief   Process charges in DFTB.
C>
C>    @details Store charges in DFTB.
C>
C>    @author  Yoshio Nishimoto
C>
C
      subroutine dftb_cnstchamul(mode,mode2,popmat,popzat,
     *                           shiftg,zref,nfg,natfrg,iatfrg,
     *                           ispefmo,maxnat,natfmob,dftb3)
      use mx_limits, only: mxatm
c
      IMPLICIT NONE
C     Declarations of arguments
      DOUBLE PRECISION popmat(maxnat,*),popzat(maxnat,*)
      INTEGER mode,mode2
      DOUBLE PRECISION shiftg(natfmob,2),zref(*)
      INTEGER nfg,natfrg(*),iatfrg(*)
      INTEGER ispefmo(*),maxnat,natfmob
      LOGICAL dftb3
C     Declarations of common blocks
      DOUBLE PRECISION ETEMP,DFTBDP,DAMPXHE,HUBDER
      DOUBLE PRECISION zz,SPNCST,SPE
      DOUBLE PRECISION PARAMDIR
      INTEGER NSPE
      INTEGER MAXANG,ISPE,IND,IDFTBD
      DOUBLE PRECISION X
C     Other declarations
C
      INTEGER, PARAMETER :: MXSPE=10
      INTEGER i,iat,ifg,indi,nati

c
      common /dftbpr/ etemp,dftbdp(mxspe*14),dampxhe,hubder(mxspe),
     *                zz(mxatm),spncst(6,mxspe),spe(mxatm),nspe,
     *                maxang(mxatm),ispe(mxatm),ind(mxatm+1),idftbd,
     *                PARAMDIR
      common /fmcom / x(1)
c
      if (mode2.eq.0) then
        indi = 1
        do ifg = 1, nfg
          nati = natfrg(ifg)
          !! check the third index of popmat
          call dcopy(nati,popzat(1,ifg),1,shiftg(indi,1),1)
          indi = indi + nati
        end do
      else
        call dcopy(natfmob,popzat,1,shiftg,1)
      end if
      if (mode.eq.0) call daxpy(natfmob,-1.0d+00,zref,1,shiftg,1)
      if (dftb3) then
        if (mode.eq.0) then
          do i = 1, natfmob
            shiftg(i,2) = shiftg(i,1)**2
          end do
        else
          indi = 1
          do ifg = 1, nfg
            nati = natfrg(ifg)
            do iat = 1, nati
              shiftg(indi,2)
     *          = shiftg(indi,1)*(popmat(iat,ifg)-zref(indi))
              indi = indi + 1
            end do
          end do
        end if
      end if
C     write (*,*) "cnstchamul"
C     do i = 1, natfmob
C       if (dftb3) then
C         write (*,'(i3,2f20.10)') i,shiftg(i,1),shiftg(i,2)
C       else
C         write (*,'(i3, f20.10)') i,shiftg(i,1)
C       end if
C     end do
c
      return
c
      end subroutine dftb_cnstchamul
