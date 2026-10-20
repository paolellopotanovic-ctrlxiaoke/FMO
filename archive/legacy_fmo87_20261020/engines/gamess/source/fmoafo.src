C  6 Jun 18 - YN  - new module for FMO/AFO gradient
C
C*MODULE FMOAFO    *DECK DFTB_UROT
C>
C>    @brief   FMO(-DFTB)/AFO analytic derivative
C>
C>    @details Compute U matrix which orthonormalizes LMOs adopted to
C>             the fragment. The U matrix is stored in dictrionary file
C>
C>    @author Yoshio Nishimoto
C>    - Feb, 2018- Subroutine written
C>
C>           --- INPUT ---
C>    @param L1      Number of AOs
C>    @param NMOS    Number of frozen occupied + projecte virtual orbs
C>    @param VEC0    L_{\mu i} before orthonormalization
C>    @param VEC     L_{\mu i} after orthonormalization
C>
C
C     See the FMO-DFTB/AFO reference for more details:
C       Nishimoto, Y.; Fedorov, D. G.
C       J. Chem. Phys. 2018, 148, 064115.
C
C     Some subroutines may be used in HF/DFT. DFTB_LOCZVEC, DFTB_PUTGMK,
C     and DFTB_ZVEC_ZB have to be modified at least.
C
      SUBROUTINE DFTB_UROT(L1,NMOS,VEC0,VEC)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER (ZERO=0.0D+00,ONE=1.0D+00)
C
      COMMON /FMCOM / X(1)
      COMMON /IOFILE/ IR,IW,IP,IJK,IJKT,IDAF,NAV,IODA(950)
C
      DIMENSION VEC0(L1,*),VEC(L1,*)
      DIMENSION WRK(1)
C
      IF (NMOS.EQ.0) THEN
        RETURN
      ELSE IF (NMOS.EQ.1) THEN
        !! Here, only normalization has been performed, so it's easy to
        !! obtain \widetilde{U}_{ki}
        DNORM  = ONE/DNRM2(L1,VEC0,1)
        CALL DAWRIT(IDAF,IODA,DNORM,1,569,0)
      ELSE
        !! QUERY THE WORKSPACE
        CALL VALFM(LOADFM)
        LAST = LOADFM
        CALL DGESVD('A','A',L1,NMOS,VEC0,L1,X(LAST),X(LAST),L1,X(LAST),
     *              NMOS,WRK,-1,INFO)
        LWORK = INT(WRK(1))
C
        LS    = LOADFM  + 1
        LU    = LS    + NMOS
        LV    = LU    + L1*L1
        LWRK  = LV    + NMOS*NMOS
        LWRK1 = LWRK  + LWORK
        LAST  = LWRK1 + L1*L1
        NEED  = LAST - LOADFM - 1
        CALL GETFM(NEED)
C
        !! L = L^0*X
        !! X = (L^0)-1*L
        !! L_{\mu i} = \sum_i L^0{\mu m} X_{mi}
C
        !! L^0 = U*S*V^T (SINGULAR VALUE DECOMPOSITION)
        !! X = (L^0)^-1 * L
        !!   = V * S^-1 * U^T * L
        !! L^0_{\mu i} = \sum_m U_{\mu m} \sigma_m V_{im}
        !! U: L1*L1
        !! S: NMOS
        !! V: NMOS*NMOS
        CALL DCOPY(L1*NMOS,VEC0,1,X(LWRK1),1)
        CALL DGESVD('A','A',L1,NMOS,VEC0,L1,X(LS),X(LU),L1,X(LV),NMOS,
     *              X(LWRK),LWORK,INFO)
        CALL DCOPY(L1*NMOS,X(LWRK1),1,VEC0,1)
C
        !! V
        CALL TRPOSQ(X(LV),NMOS)
        !! V * S^-1
        DO I = 1, NMOS
          CALL DSCAL(NMOS,ONE/X(LS+I-1),X(LV+NMOS*(I-1)),1)
        END DO
        !! WRK1 = (L^0)^-1 = V * S^-1 * U^T
        CALL DGEMM('N','T',NMOS,L1,NMOS,ONE,X(LV),NMOS,X(LU),L1,
     *             ZERO,X(LWRK1),NMOS)
        !! X = (L^0)^-1 * L
        CALL DGEMM('N','N',NMOS,NMOS,L1,ONE,X(LWRK1),NMOS,VEC,L1,
     *             ZERO,X(LV),NMOS)
        CALL DAWRIT(IDAF,IODA,X(LV),NMOS*NMOS,569,0)
C       !! L^0 * X
C       CALL DGEMM('N','N',L1,NMOS,NMOS,ONE,VEC0,L1,
C    *             X(LV),NMOS,ZERO,X(LWRK1),L1)
        CALL RETFM(NEED)
      END IF
C
      RETURN
C
      END SUBROUTINE DFTB_UROT
C
C-----------------------------------------------------------------------
C*MODULE FMOAFO    *DECK DFTB_RROT
C>
C>    @brief   FMO(-DFTB)/AFO analytic derivative
C>
C>    @details Compute R matrix from U and LSC matrices. See Eq. (27) in
C>             the FMO-DFTB/AFO paper.
C>
C>    @author Yoshio Nishimoto
C>    - Feb, 2018- Subroutine written
C>
C>           --- INPUT ---
C>    @param L1      Number of AOs
C>    @param NFRZ    Number of frozen occupied + projecte virtual orbs
C>    @param NUMFRZ  Number of frozen occupied orbitals
C>    @param UROT    U matrix
C>    @param S       S_{\mu \nu} (Overlap matrix in AO)
C>    @param VECL    L_{\mu i} (LMOs)
C>    @param VEC     C_{\mu i} (MO coefficient)
C>    @param WRK1    Work space
C>    @param WRK2    Work space
C>           --- OUTPUT ---
C>    @param RROT    R matrix
C>
C
      SUBROUTINE DFTB_RROT(L1,NFRZ,NUMFRZ,UROT,RROT,S,VECL,VEC,
     *                     WRK1,WRK2)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER (ZERO=0.0D+00,ONE=1.0D+00)
C
      DIMENSION UROT(NFRZ,NFRZ),RROT(NFRZ,NFRZ),S(*),VECL(L1,*),
     *          VEC(L1,*),WRK1(*),WRK2(L1,*)
C
      IF (NFRZ.EQ.0) RETURN
      NORBPRJ = NFRZ - NUMFRZ
C
      IF (NFRZ.GT.1) THEN
        CALL CPYTSQ(S,WRK2,L1,1)
        !! L^T*S
        CALL DGEMM('T','N',NFRZ,L1,L1,ONE,VECL,L1,
     *             WRK2,L1,ZERO,WRK1,NFRZ)
        CALL DCOPY(L1*NUMFRZ,VEC,1,WRK2,1)
        IF (NORBPRJ.GT.0) THEN
          CALL DCOPY(L1*NORBPRJ,VEC(1,L1-NORBPRJ+1),1, !! use L0?
     *               WRK2(1,NUMFRZ+1),1)
        END IF
        !! L^T*S*C
        CALL DGEMM('N','N',NFRZ,NFRZ,L1,ONE,WRK1,NFRZ,
     *             WRK2,L1,ZERO,RROT,NFRZ)
        !! WRK1 = RROT = UROT*(L^T*S*C)
        CALL DGEMM('N','N',NFRZ,NFRZ,NFRZ,ONE,UROT,NFRZ,
     *             RROT,NFRZ,ZERO,WRK1,NFRZ)
        CALL DCOPY(NFRZ*NFRZ,WRK1,1,RROT,1)
      ELSE
        RROT(1,1) = UROT(1,1)
      END IF
C
      RETURN
C
      END SUBROUTINE DFTB_RROT
C
C-----------------------------------------------------------------------
C*MODULE FMOAFO    *DECK DFTB_LOCZVEC
C>
C>    @brief   FMO(-DFTB)/AFO analytic derivative
C>
C>    @details Solve two types of Z-vector equations for model system
C>             At present, only Pipek--Mezey localization is implemented
C>
C>    @author Yoshio Nishimoto
C>    - Feb, 2018- Subroutine written
C>
C>           --- INPUT ---
C>    @param MODE    1: Perform SCF, 2: SCF has already been performed
C>    @param GMKZ    G_{\mu k}^\Xi
C>    @param MDOUT   Controls the amount of output.
C>    @param UROT    U matrix
C>    @param S       S_{\mu \nu} (Overlap matrix in AO)
C>    @param VECL    L_{\mu i} (LMOs)
C>    @param VEC     C_{\mu i} (MO coefficient)
C>    @param WRK1    Work space
C>    @param WRK2    Work space
C>           --- OUTPUT ---
C>    @param RROT    R matrix
C>
C
      SUBROUTINE DFTB_LOCZVEC(MODE,GMK,MDOUT)
      use mx_limits, only: mxatm,mxao
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      LOGICAL DFTBFL,SCC,SRSCC,DFTB3,DAMPXH,LCDFTB
      LOGICAL MDOUT,MFRZ,GOPARR,DSKWRK,MASWRK
C
      PARAMETER (ZERO=0.0D+00,HALF=0.5D+00,ONE=1.0D+00)
      PARAMETER (MXSPE=10)
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
      COMMON /FMODTB/ LGAMMAGLOB,LZREF,LISPE,LSHIFTG,LSHIFTCT,
     *                LGAMMAD,MODESD,MODGAMMA,LDCTMP,LGAMMA3G,LGAMMA3D,
     *                LNBOND,LCTMUL,LESPDFTB,LESPASC,NDFTBRST,LDFTBRST
      COMMON /FUNCT / ENERGY,EGRAD(3,MXATM)
      COMMON /GRAD  / DE(3,MXATM)
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON /MFRPAR/ MFRZ,NUMFRZ,NORFRZ,IFRZ(MXAO)
      COMMON /OPTLOC/ CVGLOC,MAXLOC,IPRTLO,ISYMLO,IFCORE,NOUTA,NOUTB,
     *                MOOUTA(MXAO),MOOUTB(MXAO),IBOYAL
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /PRPOPT/ ETOLLZ,ILOCAL,IAHARD
C
      DIMENSION GMK(NUM,*)
      DIMENSION AU(9,9),BU(9,9),V(3)
      LOGICAL SKP
C
      if (.not.dftbfl) return !! only for DFTB is implemented
C
C     MOVE TO A MODEL SYSTEM
C
C     IF (ILOCAL.EQ.0) GO TO 100
      IF (MODE.EQ.0) THEN
        ITMP1  = NUMFRZ
        ITMP2  = NORFRZ
        ITMP3  = NOUTA
        NUMFRZ = 0
        NORFRZ = 0
        NOUTA  = 0
        MFRZ   = .FALSE.
        CALL DFTBX(DUM,DUM) !ENERGX
        CALL LMOLOC(ILOCAL)
        CALL LMOPSI(1,MCORE,MDOC,MACT,NUMLOC)
        NLOC   = NUMLOC
        NUMFRZ = ITMP1
        NORFRZ = ITMP2
        NOUTA  = ITMP3
        MFRZ   = .TRUE.
      ELSE IF (MODE.EQ.1) THEN
        CALL VCLR(EGRAD,1,3*NAT)
      END IF
C
      L0 = NQMT
      L1 = NUM
      L2 = L1*(L1+1)/2
      L3 = L1*L1
C
      NLOC  = NA
      NVIRT = L1-NA
      NOCC  = NA
      NNUM  = MAX(NLOC*NLOC,NVIRT*NOCC)
C
      CALL VALFM(LOADFM)
      LWRK      = LOADFM    + 1
      LEIG      = LWRK      + L1
      LS        = LEIG      + L1
      LD        = LS        + L2
      LW        = LD        + L2
      LFAO      = LW        + L2
      LWRK1     = LFAO      + L2
      LWRK2     = LWRK1     + L3
      LWRK3     = LWRK2     + L3
      LLAG1     = LWRK3     + L3
      LLAG2     = LLAG1     + MAX(L2,NLOC*NLOC)
      LZTMP     = LLAG2     + MAX(L2,NOCC*NVIRT)
      LVEC      = LZTMP     + L3
      LVECL     = LVEC      + L3
      LAST      = LVECL     + L3
C
      LCHAMUL   = LAST
      LSHIFT    = LCHAMUL   + NAT
      LSHIFTZ   = LSHIFT    + NAT
      LDERSHIFT = LSHIFTZ   + NAT
      LDIST     = LDERSHIFT + NAT
      LGAMMA    = LDIST     + NAT*(NAT+1)/2
      LAST      = LGAMMA    + NAT*(NAT+1)/2
      IF (DFTB3) THEN
        LGAMMA3 = LAST
        LAST    = LGAMMA3   + NAT*NAT
      END IF
C
      LR0  = LAST
      LZ0  = LR0    + NNUM
      LP0  = LZ0    + NNUM
      LAP  = LP0    + NNUM
      LZ1  = LAP    + NNUM
      LPREC= LZ1    + NNUM
      LAST = LPREC  + NNUM
C
      LQQ       = LAST
      LQBAR     = LQQ       + NAT*NLOC*NLOC
      LABC      = LQBAR     + NAT*NLOC*L1
      LMM       = LABC      + NLOC*NLOC*L1*L1
      LAST      = LMM       + NAT*L1*L1
C
      NEED  = LAST - LOADFM - 1
      CALL GETFM(NEED)
C
      CALL DAREAD(IDAF,IODA,X(LS)   ,L2,12,0)
C     CALL DFTB_HAM_OVER(X(LS),X(LS),.FALSE.,.TRUE.,.FALSE.)
      CALL DAREAD(IDAF,IODA,X(LVEC) ,L3,15,0)
      CALL DAREAD(IDAF,IODA,X(LD)   ,L2,16,0)
      CALL DAREAD(IDAF,IODA,X(LEIG) ,L1,17,0)
      CALL DAREAD(IDAF,IODA,X(LVECL),L3,71,0)
      !! WRK1 = L_{\mu i}
C
      NSEQ = 0
      DO I = 1, NAT
        DO J = 1, I
          NSEQ = NSEQ + 1
          DIST = SQRT( (C(1,I)-C(1,J))**2
     *               + (C(2,I)-C(2,J))**2
     *               + (C(3,I)-C(3,J))**2 )
          X(LDIST+NSEQ-1) = DIST
        END DO
      END DO
      CALL DFTB_GAMMA(DAMPXHE,HUBBL,X(LGAMMA),X(LDIST),NAT,NSPE,
     *                ISPE,MAXANG,SPE,DAMPXH,SRSCC)
      IF (DFTB3) THEN
        CALL VCLR(X(LGAMMA3),1,NAT*NAT)
        CALL DFTB_GAMMA3(NAT,HUBBL,X(LGAMMA3),X(LDIST))
      END IF
      CALL DAREAD(IDAF,IODA,X(LCHAMUL),NAT,556,0)
      DO I = 1, NAT
        X(LCHAMUL+I-1) = X(LCHAMUL+I-1) - ZREF(I)
      END DO
C
C     ----- COMPUTE A NON-RESPONSE TERM -----
C
      CALL CPYTSQ(X(LD),X(LWRK2),L1,1)
      !! D_{rho mu}*G_{mu k}
      CALL DGEMM('T','N',L1,L1,L1,ONE,X(LWRK2),L1,
     *           GMK,L1,ZERO,X(LWRK1),L1)
      !! D_{rho mu}*G_{mu k}*L_{sigma k}
      CALL DGEMM('N','T',L1,L1,NLOC,-ONE,X(LWRK1),L1,
     *           X(LVECL),L1,ZERO,X(LWRK2),L1)
      CALL SYMSQT(L1,X(LWRK2),L1,X(LW),1)
      !! Now, X(LW) is Eq. (75)
C
C     ----- CONSTRUCT TWO L TERMS -----
C
      !! L for U^a
      CALL DAREAD(IDAF,IODA,X(LWRK2),NLOC*NLOC,73,0)
      !! C_{mu m}*G_{mu k}
      CALL DGEMM('T','N',L1-NA,NLOC,L1,ONE,X(LVEC+L1*NA),L1,
     *           GMK,L1,ZERO,X(LWRK1),L1-NA)
      !! C_{mu m}*G_{mu k}*W_{n k}
      CALL DGEMM('N','T',L1-NA,NA,NLOC,ONE,X(LWRK1),L1-NA,
     *           X(LWRK2),NLOC,ZERO,X(LLAG2),L1-NA)
      !! 4*C_{mu m}*G_{mu k}*W_{n k}
      CALL DSCAL(NA*(L1-NA),4.0D+00,X(LLAG2),1)
C
      !! L for V^a
      CALL DGEMM('T','N',NLOC,NLOC,L1,ONE,X(LVECL),L1,
     *           GMK,L1,ZERO,X(LLAG1),NLOC)
      CALL DSCAL(NLOC*NLOC,8.0D+00,X(LLAG1),1)
      CALL ASYMTRZE(X(LLAG1),NLOC,NLOC,ASYM)
      !! X(LLAG1) is asymmetrized Eq. (77)
C
C     SOLVE Z-VECTOR USING COUPLED-PERTURBED LOCALIZATION
C
C     (A) CONSTRUCT INITIAL RESIDUE R0
      NNUM = NLOC*NLOC
      CALL DCOPY(NNUM,X(LLAG1),1,X(LR0),1)
C
      CALL DFTB_CNSTQ(1,NAT,NLOC,L1,IND,X(LS),X(LVEC),X(LVECL),X(LQQ),
     *                X(LQBAR))
      CALL DFTB_CNSTC(NAT,NLOC,X(LABC),X(LQQ))
C
C     (B) CONSTRUCT Z0 FROM R0 AND MAKE P0=Z0
C     (C) PREPARE PRECONDITIONER
      !! Preconditioning is not described in the paper...
      CALL DFTB_PRECPL(NAT,NLOC,X(LQQ),X(LPREC))
      DO I = 1, NNUM
        X(LZ0+I-1) = X(LR0+I-1)*X(LPREC+I-1)
      END DO
C
C     ----- START ITERATION -----
C
      !! Eq. (87)
      IDUM=0
      CALL DFTB_ZVEC(3,L1,L2,NNUM,NA,NAT,NLOC,X(LR0),
     *               X(LZ0),X(LP0),X(LAP),X(LZ1),X(LPREC),X(LZTMP),
     *               X(LSHIFT),X(LWRK1),X(LWRK2),X(LABC),
     *               X(LVEC),X(LFAO),X(LS),X(LEIG),X(LGAMMA),
     *               X(LGAMMA3),X(LCHAMUL),
     *               NSPE,ISPE,IND,MAXANG,IW,MASWRK,X(LVEC),IDUM,
     *               MDOUT,DFTB3,X(LAST),MODGAMMA)
C
C     ----- COMPUTE Z*A AND ADD TO ENERGY-WEIGHTED DENSITY MATRIX  -----
C
      !! Eq. (81)
      CALL DFTB_CNSTA(NAT,NA,L1,NLOC,X(LABC),X(LMM),X(LQQ),X(LVECL),
     *   IND,X(LZ1),X(LWRK1),X(LW))
      !! X(LW) is now the first term in Eq. (93) without overlap derivat
C
C     ----- COMPUTE B*Z AND ADD TO L^U -----
C
      CALL DAREAD(IDAF,IODA,X(LWRK1),NLOC*NLOC,73,0)
      !! Eq. (82)
      CALL DFTB_CNSTB(NAT,NA,L1,NLOC,X(LABC),X(LQQ),X(LQBAR),X(LWRK1),
     *     X(LLAG2),X(LZ1))
      !! X(LLAG2) is now the LHS of Eq. (92)
C
C     SOLVE Z-VECTOR USING STANDARD COUPLED-PERTURBED EQUATION
C
      NNUM = (L1-NA)*NA
C     (A) CONSTRUCT INITIAL RESIDUE R0
      CALL DCOPY(NNUM,X(LLAG2),1,X(LR0),1)
C     (B) CONSTRUCT Z0 FROM R0 AND MAKE P0=Z0
C     (C) PREPARE PRECONDITIONER
      NSEQ = 0
      DO JO = 1, NA
        DO IO = NA+1, L1
          NSEQ = NSEQ + 1
          X(LPREC+NSEQ-1) = ONE/(X(LEIG+JO-1)-X(LEIG+IO-1))
          X(LZ0+NSEQ-1) = X(LR0+NSEQ-1)*X(LPREC+NSEQ-1)
        END DO
      END DO
C
C     ----- START ITERATION -----
C
      !! Eq. (92)
      CALL DFTB_ZVEC(2,L1,L2,NNUM,NA,NAT,0,X(LR0),
     *               X(LZ0),X(LP0),X(LAP),X(LZ1),X(LPREC),X(LZTMP),
     *               X(LSHIFT),X(LWRK1),X(LWRK2),X(LWRK3),
     *               X(LVEC),X(LFAO),X(LS),X(LEIG),X(LGAMMA),
     *               X(LGAMMA3),X(LCHAMUL),
     *               NSPE,ISPE,IND,MAXANG,IW,MASWRK,X(LVEC),IDUM,
     *               MDOUT,DFTB3,X(LAST),MODGAMMA)
C
C     ----- COMPUTE B*Z -----
C
      CALL SYMSQT(L1,X(LZTMP),L1,X(LWRK1),1)
      CALL DCOPY(L2,X(LWRK1),1,X(LZTMP),1)
      DO NI = 1, NAT
        TMP = ZERO
        DO MM = 1, IND(NI+1)-IND(NI)
          MU = IND(NI)+MM
          DO NU = 1, L1
            CALL DFTB_CNVSQ(MU,NU,NSEQ)
            TMP = TMP + X(LZTMP+NSEQ-1)*X(LS+NSEQ-1)
          END DO
        END DO
        X(LWRK2+NI-1) = TMP
      END DO
      CALL VCLR(X(LSHIFTZ),1,NAT)
      CALL DFTB_2ND_SHIFT(X(LSHIFTZ),X(LWRK2),X(LGAMMA),NAT)
      IF (DFTB3) CALL DFTB_TD_3RD_SHIFT(X(LSHIFTZ),X(LSHIFTZ),X(LWRK2),
     *                X(LWRK2),X(LCHAMUL),X(LGAMMA3),NAT,.TRUE.)
      CALL VCLR(X(LFAO),1,L2)
      CALL DFTB_SHIFT_FOCK(X(LSHIFTZ),X(LFAO),X(LS),NAT,L2,NSPE,ISPE,
     *  IND,MAXANG)
      CALL CPYTSQ(X(LD),X(LWRK2),L1,1)
      CALL TFTRI(X(LWRK1),X(LFAO),X(LWRK2),X(LWRK),L0,L1,L1)
      CALL DAXPY(L2,-HALF,X(LWRK1),1,X(LW),1)
C
      !! construct W^Z
      CALL DCOPY(NNUM,X(LZ1),1,X(LWRK1),1)
      DO JO = 1, NA
        CALL DSCAL(L1-NA,X(LEIG+JO-1),X(LWRK1+(L1-NA)*(JO-1)),1)
      END DO
      CALL DGEMM('N','N',L1,NA,L1-NA,ONE,X(LVEC+L1*NA),L1,
     *           X(LWRK1),L1-NA,ZERO,X(LWRK2),L1)
      CALL DGEMM('N','T',L1,L1,NA,ONE,X(LWRK2),L1,
     *           X(LVEC),L1,ZERO,X(LWRK1),L1)
      CALL SYMSQT(L1,X(LWRK1),L1,X(LWRK2),1)
      CALL DAXPY(L2,-ONE,X(LWRK2),1,X(LW),1)
C
      ! construct SHIFTES with complete DQES
      CALL DAREAD(IDAF,IODA,X(LSHIFT),NAT,559,0)
C
      NN=NSPE
      CALL DFTB_GAMDER(NAT,NSPE,HUBBL,ISPE,X(LGAMMA),DAMPXH,
     *  DAMPXHE,SPE,X(LDIST),DFTB3,X(LGAMMA3),HUBDER,SRSCC,MAXANG,
     *  NAT,NN)
      IF (DFTB3) CALL DSCAL(NAT*NAT,3.0D+00,X(LGAMMA3),1)
C
      CALL VCLR(AU,1,81)
      CALL VCLR(BU,1,81)
      !! Finally, compute Eq. (93)
      DO NI = 1, NAT 
        DO KI = 1, 3
C
C         ----- CALCULATE HDERAO AND SDERAO TERM -----
C         HDERAO: H_{\mu \nu}^{0,a}
C         SDERAO: S_{\mu \nu}^{a}
C
          CALL VCLR(X(LLAG1),1,L2)
          CALL VCLR(X(LLAG2),1,L2)
          DO NJ = 1, NAT
            IF (NI.EQ.NJ) CYCLE
            V(1) = C(1,NJ)-C(1,NI)
            V(2) = C(2,NJ)-C(2,NI)
            V(3) = C(3,NJ)-C(3,NI)
            CALL DFTB_SKODE_HELPA(1,NI,NJ,KI,0,0,0,V,AU,BU,.TRUE.,
     *           .TRUE.,SKP)
            IF (SKP) CYCLE
            DO M = 1, IND(NJ+1) - IND(NJ)
              MU = IND(NJ)+M
              DO N = 1, IND(NI+1) - IND(NI)
                NU = IND(NI)+N
                CALL DFTB_CNVSQ(MU,NU,NSEQ)
                X(LLAG1+NSEQ-1) = AU(N,M) !! HDERAO
                X(LLAG2+NSEQ-1) = BU(N,M) !! SDERAO
              END DO
            END DO
          END DO
C
C         ----- CALCULATE F_{\mu \nu}^a TERMS FOR IFG -----
C
          !! S_{\mu \nu}/da * \Omega_{AB}
          CALL DFTB_SHIFT_FOCK(X(LSHIFT),X(LLAG1),X(LLAG2),NAT,L2,NSPE,
     *      ISPE,IND,MAXANG)
          !! S_{\mu \nu} * \Omega_{AB}/da
          CALL DFTB_CALC_DERSHIFT(NAT,NI,KI,C,X(LDIST),X(LGAMMA),
     *      X(LGAMMA3),X(LCHAMUL),X(LDERSHIFT),DFTB3)
          CALL DFTB_SHIFT_FOCK(X(LDERSHIFT),X(LLAG1),X(LS),NAT,L2,
     *      NSPE,ISPE,IND,MAXANG)
C
C         ----- CALCULATE Z_{\mu \nu} F_{\mu \nu}^a TERM -----
C
          EGRAD(KI,NI) = EGRAD(KI,NI) + TRACEP(X(LZTMP),X(LLAG1),L1)
C
          CALL VCLR(X(LFAO),1,L2)
          CALL DFTB_SHIFT_FOCK(X(LSHIFTZ),X(LFAO),X(LLAG2),NAT,L2,NSPE,
     *      ISPE,IND,MAXANG)
          EGRAD(KI,NI) = EGRAD(KI,NI) + TRACEP(X(LD),X(LFAO),L1)
C
C         ----- CALCULATE \tilde{Z}_{\mu \nu} S_{\mu \nu}^a TERM -----
C
          EGRAD(KI,NI) = EGRAD(KI,NI) + TRACEP(X(LW),X(LLAG2),L1)
        END DO
      END DO
C
      IF (MODE.EQ.1) CALL DCOPY(3*NAT,EGRAD,1,DE,1)
C     write (*,*) "DEPLMO contribution"
C     do i = 1, nat
C       write (*,'(i3,3f20.10)') i,(egrad(j,i),j=1,3)
C     end do
      CALL RETFM(NEED)
C
      RETURN
C
      END SUBROUTINE DFTB_LOCZVEC
C
C-----------------------------------------------------------------------
C*MODULE FMOAFO    *DECK DFTB_CNSTGMK
C>
C>    @brief   FMO(-DFTB)/AFO analytic derivative
C>
C>    @details Construct G_{\mu k}.
C>             See Eq. (61) in the FMO-DFTB/AFO paper.
C>
C>    @author Yoshio Nishimoto
C>    - Feb, 2018- Subroutine written
C>
C>           --- INPUT ---
C>    @param L1       Number of AOs
C>    @param NA       Number of alpha electrons
C>    @param NUMFRZ   Number of frozen occupied orbitals
C>    @param NFRZ     Number of frozen occupied + projecte virtual orbs
C>    @param SCAL     Scaling factor
C>    @param S        Overlap matrix in AO
C>    @param FAO      Fock matrix in AO
C>    @param VEC       MO coefficient
C>    @param RROT     R matrix (see DFTB_RROT)
C>    @param WRK1     Work space
C>    @param WRK2     Work space
C>    @param WRK3     Work space
C>           --- IN/OUTPUT ---
C>    @param GMK      G_{\mu k} in this fragment
C>
C
      SUBROUTINE DFTB_CNSTGMK(L1,NA,NUMFRZ,NFRZ,SCAL,
     *                        GMK,S,FAO,VEC,RROT,WRK1,WRK2,WRK3)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER (ZERO=0.0D+00,ONE=1.0D+00)
C
      DIMENSION GMK(L1,NFRZ),S(*),FAO(*),VEC(L1,L1),RROT(NFRZ,NFRZ),
     *          WRK1(L1,L1),WRK2(L1,L1),WRK3(L1,L1)
C
      CALL CPYTSQ(FAO,WRK1,L1,1)
      !! F_{mu nu}*C_{nu j}
      CALL DGEMM('N','N',L1,NA,L1,ONE,WRK1,L1,VEC,L1,ZERO,WRK2,L1)
      DO KO = 1, NFRZ
        DO MO = 1, NUMFRZ
          CALL DAXPY(L1,SCAL*RROT(KO,MO),WRK2(1,MO),1,GMK(1,KO),1)
        END DO
      END DO
C
      !! construct full F_{mn}
      CALL DGEMM('T','N',NUMFRZ,NA,L1,ONE,VEC,L1,
     *           WRK2,L1,ZERO,WRK1,NUMFRZ)
      !! R_{km}*F_{mn}
      CALL DGEMM('N','N',NFRZ,NA,NUMFRZ,ONE,RROT,NFRZ,
     *           WRK1,NUMFRZ,ZERO,WRK2,NFRZ)
      !! S_{mu nu}*C_{nu n}
      CALL CPYTSQ(S,WRK1,L1,1)
      CALL DGEMM('N','N',L1,NA,L1,ONE,WRK1,L1,VEC,L1,ZERO,WRK3,L1)
      !! (S_{mu nu}*C_{nu n})*(R_{km}*F_{mn})^T
      CALL DGEMM('N','T',L1,NFRZ,NA,-SCAL,WRK3,L1,WRK2,NFRZ,ONE,GMK,L1)
C
      RETURN
C
      END SUBROUTINE DFTB_CNSTGMK
C
C-----------------------------------------------------------------------
C*MODULE FMOAFO    *DECK DFTB_PUTGMK
C>
C>    @brief   FMO(-DFTB)/AFO analytic derivative
C>
C>    @details Compute G_{\mu k} and store it in appropriate place.
C>             Perhaps this corresponds to applying mapping (Eq. (73))
C>             In addition, comput F_{kl}^{\prime \prime K} and store.
C>             This is the second term in Eq. (70)
C>
C>    @author Yoshio Nishimoto
C>    - Feb, 2018- Subroutine written
C>    @date   Yoshio Nishimoto
C>    - Oct, 2019- Small modification for the shift vector
C>
C>           --- INPUT ---
C>    @param L1       Number of AOs
C>    @param ISC      Scaling factor (sometimes 1-N_{md}^I)
C>    @param GMK      G_{\mu k} (work space)
C>    @param INDFRZ   Index of frozen orbitals
C>    @param IALMO    Index of adaptive LMO?
C>    @param IAGLOB   Global atom index
C>    @param WRK      Work space
C>    @param FAO      Work space
C>    @param IPTLG    Index of Lagrangian for SCZV
C>    @param DODISTR  Data are sotred in GDDI memory pool or not
C>    @param ZVWRK    Work space
C>           --- IN/OUTPUT ---
C>    @param GMKSAV   G_{\mu k} for model systems (Eq. (70))
C>    @param ZVLAG    Calculated Lagrangian
C>    @param SHIFTSAV ESP contributions for each atom.
C>
C
      SUBROUTINE DFTB_PUTGMK(L1,ISC,GMK,GMKSAV,INDFRZ,IALMO,IAGLOB,
     *                       WRK,FAO,IPTLG,DODISTR,ZVLAG,ZVWRK,
     *                       SHIFTSAV)
      use mx_limits, only: mxsh,mxgtot,mxatm,mxao
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER (ZERO=0.0D+00,ONE=1.0D+00)
c     PARAMETER (MXSPE=10)
C
      LOGICAL DFTBFL,SCC,SRSCC,DFTB3,DAMPXH0,LCDFTB
      LOGICAL GOPARR,DSKWRK,MASWRK
C
      COMMON /DFTB  / DFTBFL,SCC,SRSCC,DFTB3,DAMPXH0,LCDFTB
c     COMMON /DFTBPR/ ETEMP,DFTBDP(MXSPE*14),DAMPXHE,HUBDER(MXSPE),
c    *                ZREF(MXATM),SPNCST(6,MXSPE),SPE(MXATM),NSPE,
c    *                MAXANG(MXATM),ISPE(MXATM),IND(MXATM+1),IDFTBD
      COMMON /FMCOM / X(1)
      COMMON /FMODIM/ MAXBND,MAXKND,MAXCBS,MAXCAO,MAXVEC,MAXL1,MAXNAT,
     *                MAXABD,MAXBAS,MAXBBD,MAXLMO,MAXSLO,MAXABD2,MAXRIJ
      COMMON /FMOOPT/ ESPSCA(9),RESPAP(2),RESPPC(2),RESDIM,RESTRI(4),
     *                RCORSD,RESPCT,CONVFG,CNVDMP,COROFF,RFLMO(4),
     *                ORSHFT,ORSHFT2,CNVAFO,ASCREEN(4),IXESP,MXITFG,
     *                NGUESS,NBSSE,MODORB,MODPAR,IRSTSTP,IRSTLAY,NPRFMO,
     *                NFMOPAL,MODPRP,MAXL1C,IPIEDA,MODGRD,MODESP,IVMUL,
     *                MODLMO,NOPDEN,MOFOCK,MODFD,MODFMM,NCENTM,ndualb,
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
      common /fmorun/ espscf,e0scf(2),emp2s,IDAFMO,icurfg,jcurfg,kcurfg,
     *                icurlay,icurunt,nat1e,ncursh,ngau,icurpop,ifmostp,
     *                moncor,needr,modrst,norbproj,nunesp,iskipesp,
     *                IESDPPC,idoprop,mp2run,icurit,idmfmo,iddfmo,
     *                iddcur,nddleft,ivmfmo,nzmtfmo,ifmobas,itmfmo(2)
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
      DIMENSION GMK(L1,*),GMKSAV(MAXBBD,MAXABD2,MAXLMO,*),
     *          INDFRZ(2,*),IALMO(2,MAXABD2,MAXLMO,*),IAGLOB(*),
     *          WRK(*),FAO(*),IPTLG(*),ZVLAG(*),ZVWRK(L1-NA,*),
     *          SHIFTSAV(*)
      LOGICAL   MFRZ,DODISTR
C
      NFRZ = MAX(NUMFRZ,NORFRZ)
      IF (NFRZ.EQ.0) RETURN
      SCAL = DBLE(ISC)
C
      if (.not.dftbfl) return !! only for DFTB is implemented
C
      IF (JCURFG.EQ.0.AND.IAND(MODGRD,32).NE.0) THEN
        L2=L1*(L1+1)/2
        L3=L1*L1
C
        INDI = IXFTCH(X(LINDFRG),ICURFG)
        CALL DAREAD(IDAF,IODA,WRK,NAT,559,0) !! internal shift (+PCM)
        IF (MASWRK) CALL DAXPY(NAT,ONE,WRK,1,SHIFTSAV(INDI),1)
        !! October 2019
        !! Change due to a different definition of the internal gradient
        !! The shift in 559 already contains the ESP shift
C       CALL DAREAD(IDAF,IODA,WRK,NAT,566,0) !! ESP for DFTB2
C       IF (MASWRK) CALL DAXPY(NAT,ONE,WRK,1,SHIFTSAV(INDI),1)
C       IF (DFTB3) THEN
C         CALL DAREAD(IDAF,IODA,FAO,NAT,556,0)
C         CALL DAREAD(IDAF,IODA,WRK,NAT*2,567,0) !! ESP for DFTB3
C         DO I = 1, NAT
C           DQ = FAO(I) - ZREF(I)
C           IF (MASWRK) SHIFTSAV(INDI+I-1) = SHIFTSAV(INDI+I-1)
C    *        + 2.0D+00*WRK(I)*DQ
C         END DO
C         IF (MASWRK) CALL DAXPY(NAT,ONE,WRK(NAT+1),1,SHIFTSAV(INDI),1)
C       END IF
C
        CALL VALFM(LOADFM)
        LS    = LOADFM  + 1
        LWRK1 = LS      + L2
        LWRK2 = LWRK1   + L3
        LWRK3 = LWRK2   + L3
        LVEC  = LWRK3   + L3
        LVECL = LVEC    + L3
        LUROT = LVECL   + L3
        LRROT = LUROT   + NFRZ*NFRZ
        LAST  = LRROT   + NFRZ*NFRZ
        NEED  = LAST - LOADFM - 1
        CALL GETFM(NEED)
C
        CALL DAREAD(IDAF,IODA,X(LVEC),L3,15,0)
        CALL DAREAD(IDAF,IODA,FAO,L2,317,0)
        CALL CPYTSQ(FAO,X(LWRK1),L1,1)
        CALL DGEMM('T','N',L1-NA,L1,L1,ONE,X(LVEC+L1*NA),L1,X(LWRK1),L1,
     *             ZERO,X(LWRK2),L1-NA)
        CALL DGEMM('N','N',L1-NA,NA,L1,ONE,X(LWRK2),L1-NA,
     *             X(LVEC),L1,ZERO,ZVWRK,L1-NA)
        !! Second term in Eq. (70)
        CALL ZVLGMULT(NA*(L1-NA),ISC,ICURFG,IPTLG,
     *       DODISTR,ZVLAG,ZVWRK,X(LWRK1))
C
C       ----- CONSTRUCT GMK FOR MONOMER (Eq. (61)) -----
C
        CALL VCLR(GMK,1,NFRZ*L1)
        IF (NUMFRZ.GT.0) THEN
          CALL DAREAD(IDAF,IODA,X(LS),L2,12,0)
          CALL DAREAD(IDAF,IODA,X(LVECL),L3,318,0)
          CALL DAREAD(IDAF,IODA,X(LUROT),NFRZ*NFRZ,569,0)
          CALL DFTB_RROT(L1,NFRZ,NUMFRZ,X(LUROT),X(LRROT),X(LS),
     *                   X(LVECL),X(LVEC),X(LWRK1),X(LWRK2))
          CALL DFTB_CNSTGMK(L1,NA,NUMFRZ,NFRZ,ONE,
     *                      GMK,X(LS),FAO,X(LVEC),X(LRROT),X(LWRK1),
     *                      X(LWRK2),X(LWRK3))
        END IF
C
        CALL RETFM(NEED)
      END IF
C
      IF (.NOT.MASWRK) RETURN
C
C     WRITE (*,*) "G_{\MU K}"
C     CALL PRSQ(GMK,NFRZ,L1,L1)
      MAXIA = MAXABD2
      DO IO = 1, NFRZ
        IBDFG = INDFRZ(1,IO)
        IMO   = INDFRZ(2,IO)
        DO IA = 1, MAXIA
          IBDA=IALMO(1,IA,IMO,IBDFG)
          IF(IBDA.NE.0) THEN
            NBBDA=IALMO(2,IA,IMO,IBDFG)
            DO JA = 1, NAT
              IF (IAGLOB(JA).EQ.IBDA) THEN
               DO I=1,NSHELL
                 IF (KATOM(I).EQ.JA) THEN
                   CALL DAXPY(NBBDA,SCAL,GMK(KLOC(I),IO),1,
     *                              GMKSAV(1,IA,IMO,IBDFG),1)
                   GO TO 100
                 END IF
               END DO
              END IF
            END DO
  100       CONTINUE
          END IF
        END DO
      END DO
C
      RETURN
C
      END SUBROUTINE DFTB_PUTGMK
C
C-----------------------------------------------------------------------
C*MODULE FMOAFO    *DECK DFTB_PUTGMK2
C>
C>    @brief   FMO(-DFTB)/AFO analytic derivative
C>
C>    @details Gather G_{\mu k} for the current model system.
C>
C>    @author Yoshio Nishimoto
C>    - Feb, 2018- Subroutine written
C>
C>           --- INPUT ---
C>    @param L1       Number of AOs
C>    @param GMKSAV   G_{\mu k} for all the model systems (Eq. (70))
C>    @param ILMO     Index of LMOs
C>    @param NLMO     Number of LMOs
C>    @param IBDFG    Index of splitted covalent bonds
C>           --- OUTPUT ---
C>    @param GMK      G_{\mu k} for the present model system
C>
C
      SUBROUTINE DFTB_PUTGMK2(L1,GMK,GMKSAV,ILMO,NLMO,IBDFG)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      COMMON /FMODIM/ MAXBND,MAXKND,MAXCBS,MAXCAO,MAXVEC,MAXL1,MAXNAT,
     *                MAXABD,MAXBAS,MAXBBD,MAXLMO,MAXSLO,MAXABD2,MAXRIJ
C
      DIMENSION GMK(L1,*),GMKSAV(MAXBBD,MAXABD2,MAXLMO,*),ILMO(*)
C
      DO IO = 1, NLMO
        JLMO = ILMO(IO)
        DO IA = 1, MAXABD2
          CALL MOCOCP(1,IA,GMK(1,JLMO),IDUM,GMKSAV(1,IA,IO,IBDFG))
        END DO
      END DO
C
      RETURN
C
      END SUBROUTINE DFTB_PUTGMK2
C
C-----------------------------------------------------------------------
C*MODULE FMOAFO    *DECK DFTB_ZVEC_ZB
C>
C>    @brief   FMO(-DFTB)/AFO analytic derivative
C>
C>    @details Compute Z*B and Z*G. B and G are defined in Eqs. (39) and
C>             (40), respectively.
C>
C>    @author Yoshio Nishimoto
C>    - Feb, 2018- Subroutine written
C>
C>           --- INPUTT ---
C>    @param MODE      MODE=0: called from dimers/trimers calculation
C>                     MODE=1: called after self-consistent Z-vector
C>    @param L1        Number of AOs
C>    @param NNUM      (No. of virtual)*(No. of non-frozen occupied)
C>    @param NFRZ      Number of frozen occupied + projecte virtual orbs
C>    @param Z1        Z after solving Z-vector
C>    @param WW        W_{\mu \nu} contributions
C>    @param DQ        Mulliken charge (not population)
C>    @param DQES      Charges induced by Z
C>    @param SHIFTES   Shift contributions from the charges induced by Z
C>    @param SHIFT3    DFTB3 shift
C>    @param EIG       Orbital energy
C>    @param S         Overlap matrix in AO
C>    @param D         Density matrix in AO
C>    @param FAO       Fock matrix in AO
C>    @param EIGL      Fock matrix in LMO
C>    @param VEC        MO coefficient
C>    @param VECL      LMO coefficient
C>    @param UROT      U matrix (see DFTB_UROT)
C>    @param RROT      R matrix (see DFTB_RROT)
C>    @param ZTMP      Work space for various Z
C>    @param SL        Work space (specifically S*L is stored)
C>    @param GAMMA     gamma for DFTB2
C>    @param GAMMA3    Gamma for DFTB3
C>    @param GAMMA3FMO Gamma for FMO-DFTB3
C>    @param WRK       Work space
C>    @param WRK1      Work space
C>    @param WRK2      Work space
C>    @param WRK3      Work space
C>    @param WRK4      Work space
C>    @param DEPLMO    Whether LMOs dependent on geometry or not
C>    @param DFTB3     Whether DFTB3 or not
C>    @param MODGAMMA  How gamma/Gamma are stored
C>    @param NATFMO    Number of atoms in the whole system
C>    @param NBDFG     Number of BDA-BAA pairs
C>    @param IAGLOB    Global atom index
C>    @param NATFRG    Number of atoms in each fragment
C>    @param INDFRG    First atom of each fragment in terms of NATFMOB?
C>    @param IATFRG    Global atom index with NBDFG?
C>    @param POPMATI   Mulliken populations of the fragment
C>    @param ISPEFMO   ISPE for all atoms (NATFMO)
C>    @param FMOC      Coordinates for all NATFMO
C>    @param INDAT     Converts global atomic location to fragment index
C>    @param IALOC     Converts global to local (fragment) atomic locati
C>           --- OUTPUT ---
C>    @param ZZ       D on entry and D+Z on exit
C>    @param WW       W_{\mu \nu} contributions
C>    @param GMK      G_{\mu k} in this fragment
C>
C
      SUBROUTINE DFTB_ZVEC_ZB(MODE,L1,NNUM,NFRZ,Z1,ZZ,WW,DQ,DQES,
     *                        SHIFTES,SHIFT3,EIG,S,D,FAO,EIGL,VEC,VECL,
     *                        UROT,RROT,GMK,ZTMP,SL,GAMMA,GAMMA3,
     *                        GAMMA3FMO,WRK,WRK1,WRK2,WRK3,WRK4,DEPLMO,
     *                        DFTB3,MODGAMMA,NATFMO,NBDFG,IAGLOB,NATFRG,
     *                        INDFRG,IATFRG,POPMATI,ISPEFMO,FMOC,INDAT,
     *                        IALOC,ESPDFTB)
      use mx_limits, only: mxatm,mxao
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER (ZERO=0.0D+00,ONE=1.0D+00,TWO=2.0D+00,HALF=0.5D+00)
      PARAMETER (MXSPE=10)
      LOGICAL MFRZ
C
      COMMON /DFTBPR/ ETEMP,DFTBDP(MXSPE*14),DAMPXHE,HUBDER(MXSPE),
     *                ZREF(MXATM),SPNCST(6,MXSPE),SPE(MXATM),NSPE,
     *                MAXANG(MXATM),ISPE(MXATM),IND(MXATM+1),IDFTBD,
     *                PARAMDIR
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
C
      DIMENSION Z1(*),ZZ(*),WW(*),DQ(*),DQES(*),SHIFTES(*),SHIFT3(*),
     *          EIG(L1),S(*),D(*),FAO(*),EIGL(*),VEC(L1,*),VECL(L1,*),
     *          UROT(NFRZ,NFRZ),RROT(NFRZ,NFRZ),GMK(L1,*),ZTMP(L1,L1),
     *          SL(L1,*),GAMMA(*),GAMMA3(NAT,NAT),GAMMA3FMO(*),WRK(*),
     *          WRK1(L1,*),WRK2(L1,*),WRK3(L1,*),WRK4(L1,*),IAGLOB(*),
     *          NATFRG(*),INDFRG(*),IATFRG(*),POPMATI(*),ISPEFMO(*),
     *          FMOC(*),INDAT(*),IALOC(*),ESPDFTB(*)
      LOGICAL DEPLMO,DFTB3
C
C     Z*B: Eqs. (462)--(476) in frzorb.pdf (without ESP)
C     Z*G: Eqs. (588) in frzorb.pdf
C
      L2 = L1*(L1+1)/2
      L3 = L1*L1
      NORBPRJ = NFRZ-NUMFRZ
      IF (MODE.EQ.1) CALL VCLR(GMK,1,L1*NFRZ)
C
      !! F_{\mu \nu}^{(a)} terms. Actual computation is elsewhere.
      !! WRK2 should contain S*L
      IF (MODE.EQ.1) CALL CPYTSQ(ZZ,ZTMP,L1,1)
      CALL DCOPY(L3,ZTMP,1,WRK4,1)
      CALL DGEMM('T','N',L1,L1,L1,ONE,SL,L1,
     *           ZTMP,L1,ZERO,WRK1,L1)
      CALL DGEMM('N','N',L1,L1,L1,ONE,WRK1,L1,
     *           SL,L1,ZERO,WRK3,L1)
C
      !! Now, 6th term
      CALL DCOPY(L1*NORBPRJ,VECL(1,NUMFRZ+1),1,WRK1,1)
      DO KO0 = 1, NORBPRJ
        KO = KO0 + NUMFRZ
        CALL DSCAL(L1,WRK3(KO,KO),WRK1(1,KO0),1)
      END DO
      CALL DGEMM('N','T',L1,L1,NORBPRJ,-ONE,WRK1,L1,
     *           VECL(1,NUMFRZ+1),L1,ONE,WRK4,L1)
      IF (MODE.EQ.0) CALL DAWRIT(IDAF,IODA,WRK4,L3,476,0) !! P
C
      !! WRK1 CONTAINS Z_{\mu \nu}^\prime
      CALL SYMSQT(L1,WRK4,L1,WRK1,1)
C
      CALL VCLR(DQES,1,NAT)
      DO NI = 1, NAT
        TMP = ZERO
        DO MM = 1, IND(NI+1)-IND(NI)
          MU = IND(NI)+MM
          DO NU = 1, L1
            CALL DFTB_CNVSQ(MU,NU,NSEQ)
            TMP = TMP + WRK1(NSEQ,1)*S(NSEQ)
          END DO
        END DO
        DQES(NI) = TMP
      END DO
      !! W term from Z*F (= W^{(2),Z})
      CALL VCLR(SHIFTES,1,NAT)
      IF (MODE.EQ.1 .AND. DFTB3) THEN
        IF (IAND(MODGAMMA,4096).EQ.0) THEN
          IF (IAND(MODGAMMA,8).EQ.0) THEN
            IF (IAND(MODGAMMA,2).NE.0) THEN
              NN = 2
              NN2 = NATFMO+NBDFG
            ELSE
              NN = 1
              NN2 = NATFMO
            END IF
            CALL DFTB_OFFORBHSS(0,IAGLOB,NATFRG,INDFRG,
     *        IATFRG,POPMATI,S,WRK1,ISPEFMO,FMOC,
     *        GAMMA3FMO,SHIFTES,DQES,POPMATI,.FALSE.,
     *        INDAT,IALOC,NN,NN2)
            DO I = 1, NAT 
              DO J = 1, NAT
                GAMI = GAMMA3(I,J)
                GAMJ = GAMMA3(J,I)
                SHIFTES(I) = SHIFTES(I)
     *            + 2.0D+00/3.0D+00*DQES(J)*DQ(I)*GAMI
     *            + 2.0D+00/3.0D+00*DQES(J)*DQ(J)*GAMJ
              END DO
            END DO
          ELSE
            DO I = 1, NAT
              DO J = 1, NAT
                GAMI = GAMMA3(I,J)
                GAMJ = GAMMA3(J,I)
                SHIFTES(I) = SHIFTES(I)
     *            + 2.0D+00/3.0D+00*DQES(J)*DQ(I)*GAMI
     *            + 2.0D+00/3.0D+00*DQES(J)*DQ(J)*GAMJ
     *            + 2.0D+00/3.0D+00*DQES(I)*DQ(J)*GAMI
              END DO
              SHIFTES(I) = SHIFTES(I)
     *          + 2.0D+00*DQES(I)*ESPDFTB(I)
            END DO
          END IF
        ELSE
          CALL DFTB_TD_3RD_SHIFT(SHIFTES,SHIFTES,DQES,DQES,
     *                           DQ,GAMMA3,NAT,.TRUE.)
        END IF
      END IF
      CALL DFTB_2ND_SHIFT(SHIFTES,DQES,GAMMA,NAT)
      IF (MODE.EQ.0 .AND. DFTB3) THEN
        IF (IAND(MODGAMMA,4096).EQ.0) THEN
          !! ADD INTER-FRAGMENT
          DO NI = 1, NAT
            SHIFTES(NI) = SHIFTES(NI) + TWO*SHIFT3(NI)*DQES(NI)
          END DO
        END IF
        CALL DFTB_TD_3RD_SHIFT(SHIFTES,SHIFTES,DQES,DQES,DQ,GAMMA3,
     *       NAT,.TRUE.)
      END IF
      CALL VCLR(FAO,1,L2)
      CALL DFTB_SHIFT_FOCK(SHIFTES,FAO,S,NAT,L2,NSPE,ISPE,IND,MAXANG)
      !! WW is computed from scratch below for MODE=1.
      IF (MODE.EQ.0) THEN
        CALL CPYTSQ(D,WRK2,L1,1)
        CALL TFTRI(WRK4,FAO,WRK2,WRK,L1,L1,L1)
        CALL DAXPY(L2,HALF,WRK4,1,WW,1)
      ELSE
        CALL VCLR(ZZ,1,L2)
        CALL VCLR(WW,1,L2)
      END IF
C     !! COMPLETE D = D+Z^\prime
      CALL DAXPY(L2,ONE,WRK1,1,ZZ,1)
C
      IF (DEPLMO.AND.NUMFRZ.GT.0 .AND. MODE.EQ.0) THEN
        !! ZG term (1): 1st and 2nd terms in Eq. (40)
        CALL DFTB_CNSTGMK(L1,NA,NUMFRZ,NFRZ,-ONE,
     *                    GMK,S,FAO,VEC,RROT,WRK1,
     *                    WRK2,WRK4)
C       WRITE (*,*) "ZG1"
C       CALL PRSQ(GMK,NFRZ,L1,L1)
      END IF
C
C
C
      !! W^Z: the second term in Eq. (39)
      !! Z-vector in dimers/trimers are solved for (L1-NA)*(NA-NUMFRZ)
      !! Z-vector in SCZV are for NOCC*NVIR
      CALL DCOPY(NNUM,Z1,1,WRK1,1)
      IF (MODE.EQ.0) THEN
        DO JO = NUMFRZ+1, NA
          CALL DSCAL(L1-NA,EIG(JO),WRK1(1+(L1-NA)*(JO-1-NUMFRZ),1),1)
        END DO
        CALL DGEMM('N','N',L1,NA-NUMFRZ,L1-NA,ONE,VEC(1,NA+1),L1,
     *             WRK1,L1-NA,ZERO,WRK2,L1)
        CALL DGEMM('N','T',L1,L1,NA-NUMFRZ,ONE,WRK2,L1,
     *             VEC(1,NUMFRZ+1),L1,ZERO,WRK1,L1)
      ELSE
        CALL VCLR(WRK1,1,NUMFRZ*(L1-NA))
        DO JO = 1, NA
          CALL DSCAL(L1-NA,EIG(JO),WRK1(1+(L1-NA)*(JO-1),1),1)
        END DO
        CALL DGEMM('N','N',L1,NA,L1-NA,ONE,VEC(1,NA+1),L1,
     *             WRK1,L1-NA,ZERO,WRK2,L1)
        CALL DGEMM('N','T',L1,L1,NA,ONE,WRK2,L1,
     *             VEC(1,1),L1,ZERO,WRK1,L1)
      END IF
      CALL SYMSQT(L1,WRK1,L1,WRK4,1)
      CALL DAXPY(L2,ONE,WRK4,1,WW,1)
C
C
C
      CALL SYMTRZE(ZTMP,L1,L1)
      CALL DSCAL(L3,TWO,ZTMP,1) !! Z0 = \widetilde{Z}_{\mu \nu}
      !! Now, X(LZTMP) = \widetilde{Z}_{\mu \nu}
      !!               = (1+\tau_{ij}) Zij*(CI*CJ + CJ*CI)
      !! Z(1): 3rd term.
      CALL DAREAD(IDAF,IODA,FAO,L2,317,0)
      IF (NUMFRZ.GT.0) THEN
        ! D^fo
        CALL DGEMM('N','T',L1,L1,NUMFRZ,ONE,VEC,L1,
     *             VEC,L1,ZERO,WRK1,L1)
        CALL CPYTSQ(FAO,WRK2,L1,1)
        ! D^fo*F
        CALL DGEMM('T','N',L1,L1,L1,ONE,WRK1,L1,
     *             WRK2,L1,ZERO,WRK4,L1)
      ELSE
        CALL VCLR(WRK4,1,L3)
      END IF
      !! Z(2): 4th term?
      CALL DAREAD(IDAF,IODA,VECL,L3,318,0)
      IF (NORBPRJ.GT.0) THEN
        CALL DCOPY(L1*NORBPRJ,VECL(1,NUMFRZ+1),1,WRK1,1)
        CALL DCOPY(L1*(L1-NUMFRZ),VECL(1,NORFRZ+1),1,
     *             VECL(1,NUMFRZ+1),1)
        CALL DCOPY(L1*NORBPRJ,WRK1,1,VECL(1,L1-NORBPRJ+1),1)
      END IF
      ! first, construct M_{\mu \nu}
      CALL TFTRI(EIGL,FAO,VECL,WRK,L1,L1,L1)
      CALL DCOPY(L3,VECL,1,WRK2,1)
      DO KO = L1-NORBPRJ+1, L1
        NSEQ = KO*(KO-1)/2+KO
        CALL DSCAL(L1,EIGL(NSEQ)-ORSHFT,WRK2(1,KO),1)
      END DO
      CALL DGEMM('N','T',L1,L1,NORBPRJ,ONE,VECL(1,L1-NORBPRJ+1),L1,
     *           WRK2(1,L1-NORBPRJ+1),L1,ZERO,WRK1,L1)
      CALL CPYTSQ(S,WRK2,L1,1)
      ! M*S
      CALL DGEMM('T','N',L1,L1,L1,ONE,WRK1,L1,
     *           WRK2,L1,ONE,WRK4,L1)
      ! complete the contribution: (D^fo*F + M*S)*\widetilde{Z}
      CALL DGEMM('N','N',L1,L1,L1,ONE,WRK4,L1,
     *           ZTMP,L1,ZERO,WRK2,L1)
      CALL SYMSQT(L1,WRK2,L1,WRK4,1)
      CALL DAXPY(L2,ONE,WRK4,1,WW,1)
C
      IF (DEPLMO) THEN
        IF (NUMFRZ.GT.0) THEN
          !! ZG term (2): 3rd term in G
          IF (MODE.EQ.1) CALL DAREAD(IDAF,IODA,FAO,L2,317,0)
          CALL CPYTSQ(FAO,WRK1,L1,1)
          !! F_{mu nu}*C_{nu m}
          CALL DGEMM('T','N',L1,NUMFRZ,L1,ONE,WRK1,L1,
     *               VEC,L1,ZERO,WRK2,L1)
          !! F_{mu nu}*C_{nu m}*R_{km}
          CALL DGEMM('N','T',L1,NFRZ,NUMFRZ,ONE,WRK2,L1,
     *               RROT,NFRZ,ZERO,WRK4,L1)
        ELSE
          CALL VCLR(WRK4,1,L1*NFRZ)
        END IF
C
        !! ZG term (3): 4th term in G
        !! L_{nu m}^\prime = L_{nu m}*(epsilon_m-B)
        CALL DCOPY(L1*NORBPRJ,VECL(1,L1-NORBPRJ+1),1,WRK1,1)
        DO MO0 = 1, NORBPRJ
          MO = MO0 + L1-NORBPRJ
          NSEQ = MO*(MO-1)/2+MO
          CALL DSCAL(L1,EIGL(NSEQ)-ORSHFT,WRK1(1,MO0),1)
        END DO
        CALL CPYTSQ(S,WRK2,L1,1)
        !! S_{mu nu}*L_{nu m}^\prime
        CALL DGEMM('T','N',L1,NORBPRJ,L1,ONE,WRK2,L1,
     *             WRK1,L1,ZERO,WRK3,L1)
        CALL DGEMM('N','T',L1,NFRZ,NORBPRJ,ONE,WRK3,L1,
     *             UROT(1,NUMFRZ+1),NFRZ,ONE,WRK4,L1)
C
        !! S_{mu nu}*Z_{nu sigma}
        CALL DGEMM('T','N',L1,L1,L1,ONE,WRK2,L1,
     *             ZTMP,L1,ZERO,WRK1,L1)
        !! complete ZG term (2) + (3)
        CALL DGEMM('N','N',L1,NFRZ,L1,0.25D+00,WRK1,L1,
     *             WRK4,L1,ONE,GMK,L1)
C       WRITE (*,*) "ZG2"
C       CALL PRSQ(GMK,NFRZ,L1,L1)
      END IF
C
C
C
      !! Z(3)+Z(5): maybe 5th and 7th terms
      !! Zmk*(ek^L-B) + Zkk*emk^L
      CALL CPYTSQ(S,WRK1,L1,1)
      CALL DGEMM('N','N',L1,L1,L1,ONE,WRK1,L1,
     *           VECL,L1,ZERO,WRK2,L1) !! WRK2 = SL
      CALL DGEMM('T','N',L1,L1,L1,ONE,WRK2,L1,
     *           ZTMP,L1,ZERO,WRK1,L1)
      CALL DGEMM('N','N',L1,L1,L1,ONE,WRK1,L1,
     *           WRK2,L1,ZERO,ZTMP,L1)
      !! ZTMP = \widetilde{\widetilde{Z}}_{ij}
      IF (NUMFRZ.GT.0) THEN
        ! first, fv--fo region
        DO KO0 = 1, NORBPRJ
          KO = KO0 + L1-NORBPRJ
          DO MO0 = 1, NUMFRZ
            MO = MO0
            CALL DFTB_CNVSQ(KO,MO,NSEQ)
            WRK1(MO0+NUMFRZ*(KO0-1),1)
     *      = ZTMP(MO,KO)*(EIGL(KO*(KO-1)/2+KO)-ORSHFT)
     *      + ZTMP(KO,KO)*EIGL(NSEQ)
          END DO
        END DO
        IF (DEPLMO) THEN
          !! ZG (5)+(6): fo--fv region (5th term)
          !! the summation for k is likely redundant
          ! latter
          CALL DGEMM('N','N',L1,NORBPRJ,NUMFRZ,ONE,WRK2,L1,
     *               WRK1,NUMFRZ,ZERO,WRK3(1,NUMFRZ+1),L1)
          ! former
          CALL DGEMM('N','T',L1,NUMFRZ,NORBPRJ,ONE,
     *               WRK2(1,L1-NORBPRJ+1),L1,
     *               WRK1,NUMFRZ,ZERO,WRK3,L1)
          CALL DGEMM('N','T',L1,NFRZ,NFRZ,-0.25D+00,WRK3,L1,
     *               UROT,NFRZ,ONE,GMK,L1)
C         write (*,*) "ZG3"
C         CALL PRSQ(GMK,NFRZ,L1,L1)
        END IF
        CALL DGEMM('N','N',L1,L1,NUMFRZ,ONE,VECL,L1,
     *             WRK1,NUMFRZ,ZERO,WRK2,L1)
        CALL DGEMM('N','T',L1,L1,NORBPRJ,ONE,WRK2,L1,
     *             VECL(1,L1-NORBPRJ+1),L1,ZERO,WRK4,L1)
      ELSE
        CALL VCLR(WRK4,1,L3)
      END IF
      ! next, fv--fv region
      CALL VCLR(WRK1,1,NORBPRJ*NORBPRJ)
      DO KO0 = 1, NORBPRJ
        KO = KO0 + L1-NORBPRJ
        DO MO0 = 1, KO0
          MO = MO0 + L1-NORBPRJ
          CALL DFTB_CNVSQ(KO,MO,NSEQ)
          WRK1(MO0+NORBPRJ*(KO0-1),1)
     *      = ZTMP(MO,KO)*(EIGL(KO*(KO-1)/2+KO)-ORSHFT)
     *      + ZTMP(KO,KO)*EIGL(NSEQ)
          IF (MO.NE.KO) WRK1(MO0+NORBPRJ*(KO0-1),1) !! why!?
     *        =  WRK1(MO0+NORBPRJ*(KO0-1),1)*TWO
        END DO
      END DO
      IF (DEPLMO) THEN
        !! ZG (5)+(6): fv--fv region (5th term)
        CALL CPYTSQ(S,WRK3,L1,1)
        CALL DGEMM('N','N',L1,L1,L1,ONE,WRK3,L1,
     *             VECL,L1,ZERO,WRK2,L1) !! WRK2 = SL
        !! latter
        CALL DGEMM('N','N',L1,NORBPRJ,NORBPRJ,ONE,
     *             WRK2(1,L1-NORBPRJ+1),L1,
     *             WRK1,NORBPRJ,ZERO,WRK3,L1)
        ! former
        CALL DGEMM('N','T',L1,NORBPRJ,NORBPRJ,ONE,
     *             WRK2(1,L1-NORBPRJ+1),L1,
     *             WRK1,NORBPRJ,ONE,WRK3,L1)
        CALL DGEMM('N','T',L1,NFRZ,NORBPRJ,-0.125D+00,WRK3,L1,
     *             UROT(1,NUMFRZ+1),NFRZ,ONE,GMK,L1)
      END IF
      CALL DGEMM('N','N',L1,L1,NORBPRJ,ONE,VECL(1,L1-NORBPRJ+1),L1,
     *           WRK1,NORBPRJ,ZERO,WRK2,L1)
      CALL DGEMM('N','T',L1,L1,NORBPRJ,HALF,WRK2,L1,
     *           VECL(1,L1-NORBPRJ+1),L1,ONE,WRK4,L1)
      CALL SYMSQT(L1,WRK4,L1,WRK1,1)
      CALL DAXPY(L2,-ONE,WRK1,1,WW,1)
C
      IF (DEPLMO) THEN
        !! ZG term (4): the latter term of the fourth term?
        CALL CPYTSQ(FAO,WRK1,L1,1)
        CALL DCOPY(L1*NORBPRJ,VECL(1,L1-NORBPRJ+1),1,WRK2,1)
        DO MO0 = 1, NORBPRJ
          MO = MO0 + L1-NORBPRJ
          CALL DSCAL(L1,ZTMP(MO,MO),WRK2(1,MO0),1)
        END DO
        CALL DGEMM('T','N',L1,NORBPRJ,L1,ONE,WRK1,L1,
     *             WRK2,L1,ZERO,WRK4,L1)
        CALL DGEMM('N','T',L1,NFRZ,NORBPRJ,0.25D+00,WRK4,L1,
     *             UROT(1,NUMFRZ+1),NFRZ,ONE,GMK,L1)
C       WRITE (*,*) "ZG4"
C       CALL PRSQ(GMK,NFRZ,L1,L1)
      END IF
C
      IF (DEPLMO) CALL DSCAL(L1*NFRZ,-ONE,GMK,1)
C
      RETURN
C
      END SUBROUTINE DFTB_ZVEC_ZB
C
C-----------------------------------------------------------------------
C*MODULE FMOAFO    *DECK FRFOCK4
C>
C>    @brief   FMO(-DFTB)/AFO in SCF
C>
C>    @details Manipulate the Fock matrix in AFO, see Eq. (15) etc.
C>             Some calibrations are required, see dftbx.src.
C>             DSYMM may be used if LAPACK is linked.
C>
C>    @author Yoshio Nishimoto
C>    - Feb, 2018- Subroutine written
C>
C>           --- INPUT ---
C>    @param L1       Number of AOs
C>    @param NUMFRZ   Number of frozen occupied orbitals
C>    @param NORFRZ   Number of frozen occupied + projecte virtual orbs
C>    @param VEC      MO coefficient
C>    @param SL       S*L
C>    @param SLL1     S*L for i \in fo
C>    @param SLL2     S*L for i \in no+nv+fv
C>    @param BTMP     SL*SL*B
C>    @param WRK1     Work space
C>    @param WRK2     Work space
C>    @param WRK3     Work space
C>           --- IN/OUTPUT ---
C>    @param FAO      F_{\mu \nu}^{\prime \prime} on entry,
C>                    and F_{\mu \nu} on exit
C>
C
      SUBROUTINE FRFOCK4(L1,NUMFRZ,NORFRZ,FAO,VEC,SL,SLL1,SLL2,BTMP,
     *                   WRK1,WRK2,WRK3)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER (ZERO=0.0D+00,ONE=1.0D+00,TWO=2.0D+00)
C
      DIMENSION FAO(*),VEC(L1,*),SL(L1,*),SLL1(*),SLL2(*),BTMP(*),
     *          WRK1(L1,*),WRK2(L1,*),WRK3(*)
C
      L2 = L1*(L1+1)/2
      NORBPRJ = NORFRZ-NUMFRZ
C
      CALL CPYTSQ(FAO,WRK1,L1,1)
      CALL DGEMM('N','N',L1,L1,L1,ONE,SLL1,L1,
     *           WRK1,L1,ZERO,WRK2,L1)
CLAP
C     DO IO = 1, L1
C       NSEQ = IO*(IO-1)/2+1
C       CALL DCOPY(IO,FAO(NSEQ),1,WRK1(1,IO),1)
C     END DO
C     CALL DSYMM('R','U',L1,L1,ONE,WRK1,L1,
C    *           SLL1,L1,ZERO,WRK2,L1)
CLAP
      CALL DGEMM('N','T',L1,L1,L1,TWO,WRK2,L1,
     *           SLL2,L1,ZERO,WRK3,L1)
 
      IF (NORBPRJ.GT.0) THEN
        CALL DAXPY(L2,+ONE,BTMP,1,FAO,1)
C
        CALL DGEMM('T','N',L1,NORBPRJ,L1,ONE,WRK1,L1,
     *             VEC(1,NUMFRZ+1),L1,ZERO,WRK2,L1)
CLAP
C       CALL DSYMM('L','U',L1,NORBPRJ,ONE,WRK1,L1,
C    *             VEC(1,NUMFRZ+1),L1,ZERO,WRK2,L1)
CLAP
        CALL DCOPY(L1*NORBPRJ,SL(1,NUMFRZ+1),1,WRK1,1)
        DO KO = 1, NORBPRJ
          TMP = DDOT(L1,WRK2(1,KO),1,VEC(1,NUMFRZ+KO),1)
          CALL DSCAL(L1,TMP,WRK1(1,KO),1)
        END DO
        CALL DGEMM('N','T',L1,L1,NORBPRJ,ONE,WRK1,L1,
     *             SL(1,NUMFRZ+1),L1,ONE,WRK3,L1)
      END IF
      CALL SYMSQT(L1,WRK3,L1,WRK1,1)
      CALL DAXPY(L2,-ONE,WRK1,1,FAO,1)
C
      RETURN
C
      END SUBROUTINE FRFOCK4
C
C-----------------------------------------------------------------------
C*MODULE FMOAFO    *DECK DFTB_CNSTQ
C>
C>    @brief   FMO(-DFTB)/AFO for analytic derivative
C>
C>    @details Compute Q and QBAR (Eqs. (85) and (86), respectively).
C>             Brute-force implementation. There must be a better approach.
C>
C>    @author Yoshio Nishimoto
C>    - Feb, 2018- Subroutine written
C>
C>           --- INPUT ---
C>    @param IMODE    With IMODE=0, only Q is computed.
C>    @param NAT      Number of atoms
C>    @param NLOC     Number of localized orbitals
C>    @param L1       Number of AOs
C>    @param IND      AO index of each atom
C>    @param S        Overlap matrix
C>    @param VEC      MO coefficient
C>    @param VECL     LMO coefficient
C>           --- OUTPUT ---
C>    @param Q        Q_{ij}^{A,\Xi}, Eq. (85)
C>    @param QBAR     \bar{Q}_{ij}^{A,\Xi}, Eq. (86)
C>
C
      SUBROUTINE DFTB_CNSTQ(IMODE,NAT,NLOC,L1,IND,S,VEC,VECL,Q,QBAR)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER (ZERO=0.0D+00)
C
      DIMENSION IND(*),S(*),VEC(L1,*),VECL(L1,*),Q(NAT,NLOC,NLOC),
     *          QBAR(NAT,L1,NLOC)
C
      DO IAT = 1, NAT
        DO IO = 1, NLOC
          DO JO = 1, NLOC
            TMP1 = ZERO
            DO MU0 = 1, IND(IAT+1)-IND(IAT)
              MU = IND(IAT)+MU0
              DO NU = 1, L1
                CALL DFTB_CNVSQ(MU,NU,NSEQ)
                SS = S(NSEQ)
                TMP1 = TMP1
     *            + (VECL(MU,IO)*VECL(NU,JO)+VECL(MU,JO)*VECL(NU,IO))*SS
              END DO
            END DO
            Q(IAT,IO,JO) = TMP1
          END DO
        END DO
      END DO
C
      IF (IMODE.EQ.0) RETURN
C
      DO IAT = 1, NAT
        DO IO = 1, L1
          DO JO = 1, NLOC
            TMP2 = ZERO
            DO MU0 = 1, IND(IAT+1)-IND(IAT)
              MU = IND(IAT)+MU0
              DO NU = 1, L1
                CALL DFTB_CNVSQ(MU,NU,NSEQ)
                SS = S(NSEQ)
                TMP2 = TMP2
     *            + (VEC(MU,IO)*VECL(NU,JO)+VECL(MU,JO)*VEC(NU,IO))*SS
              END DO
            END DO
            QBAR(IAT,IO,JO) = TMP2
          END DO
        END DO
      END DO
C
      RETURN
C
      END SUBROUTINE DFTB_CNSTQ
C
C-----------------------------------------------------------------------
C*MODULE FMOAFO    *DECK DFTB_CNSTA
C>
C>    @brief   FMO(-DFTB)/AFO for analytic derivative
C>
C>    @details Compute C matrix (Eq. (83)). Brute-force implementation.
C>             There must be a better approach.
C>
C>    @author Yoshio Nishimoto
C>    - Feb, 2018- Subroutine written
C>
C>           --- INPUT ---
C>    @param NAT      Number of atoms
C>    @param NLOC     Number of localized orbitals
C>    @param Q        Q_{ij}^{A,\Xi}, Eq. (85)
C>           --- OUTPUT ---
C>    @param C        C matrix (\widetilde{C}_{ij,kl}^\Xi; Eq. (83))
C>
C
      SUBROUTINE DFTB_CNSTC(NAT,NLOC,C,Q)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER (ZERO=0.0D+00,HALF=0.5D+00,ONE=1.0D+00)
C
      DIMENSION C(NLOC,NLOC,*),Q(NAT,NLOC,NLOC)
C
      DO IO0 = 1, NLOC
        DO JO0 = 1, NLOC !IO0
          CALL DFTB_CNVSQ(IO0,JO0,NSEQ1)
          DO KO0 = 1, NLOC
            DO LO0 = 1, KO0
              CALL DFTB_CNVSQ(KO0,LO0,NSEQ2)
              C(IO0,JO0,NSEQ2) = ZERO
              DO IP = 1, 4
                IF (IP.EQ.1) THEN
                  IO = IO0
                  JO = JO0
                  KO = KO0
                  LO = LO0
                  SCAL = ONE
                ELSE IF (IP.EQ.2) THEN
                  IO = JO0
                  JO = IO0
                  KO = KO0
                  LO = LO0
                  SCAL = -ONE
                ELSE IF (IP.EQ.3) THEN
                  IO = IO0
                  JO = JO0
                  KO = LO0
                  LO = KO0
                  SCAL = -ONE
                ELSE IF (IP.EQ.4) THEN
                  IO = JO0
                  JO = IO0
                  KO = LO0
                  LO = KO0
                  SCAL = ONE
                END IF
                IF (IO.NE.LO) CYCLE
C
                TMP = ZERO
                DO IAT = 1, NAT
                  TMP = TMP + HALF*Q(IAT,KO,JO)
     *                *(Q(IAT,IO,IO)+Q(IAT,LO,LO)-2.0D+00*Q(IAT,JO,JO))
     *                + 2.0D+00*Q(IAT,IO,JO)*Q(IAT,KO,LO)
                END DO
                if (ko.eq.lo) scal=zero
                TMP = TMP*SCAL
                C(IO0,JO0,NSEQ2) = C(IO0,JO0,NSEQ2) + TMP
              END DO
            END DO
          END DO
        END DO
      END DO
C
      RETURN
C
      END SUBROUTINE DFTB_CNSTC
C
C-----------------------------------------------------------------------
C*MODULE FMOAFO    *DECK DFTB_PRECPL
C>
C>    @brief   FMO(-DFTB)/AFO for analytic derivative
C>
C>    @details Compute preconditioner for coupled-perturbed localization
C>             I'm not sure the correct formalism actually. I chose the
C>             one that gives fastest convergence.
C>
C>    @author Yoshio Nishimoto
C>    - Feb, 2018- Subroutine written
C>
C>           --- INPUT ---
C>    @param NAT      Number of atoms
C>    @param NLOC     Number of localized orbitals
C>    @param QQ       Q_{ij}^{A,\Xi}
C>           --- OUTPUT ---
C>    @param PREC     Precondition
C>
C
      SUBROUTINE DFTB_PRECPL(NAT,NLOC,QQ,PREC)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER (ZERO=0.0D+00)
C
      DIMENSION QQ(NAT,NLOC,NLOC),PREC(NLOC,NLOC)
C
      DO IO = 1, NLOC
        DO JO = 1, IO-1
          TMP = ZERO
          DO IAT = 1, NAT
            ! This leads to fastest convergence, somehow.
            TMP = TMP - 4.0D+00*QQ(IAT,IO,JO)*QQ(IAT,IO,JO)
     *          + (QQ(IAT,IO,IO)-QQ(IAT,JO,JO))**2.0D+00
            ! this reproduces 1/C
C           TMP = TMP + 2.0D+00*(QQ(IAT,IO,IO)-QQ(IAT,JO,JO)
C    *          -2.0d+00*QQ(IAT,IO,JO))*QQ(IAT,IO,JO)
C    *          + (QQ(IAT,IO,IO)-QQ(IAT,JO,JO))**2.0D+00
            ! this expression should be correct? (AFO No.46)
            ! diagonal precondition
C           TMP = TMP + 2.0D+00*(QQ(IAT,IO,IO)+QQ(IAT,JO,JO)
C    *          -2.0d+00*QQ(IAT,IO,JO))*QQ(IAT,IO,JO)
C    *          + (QQ(IAT,IO,IO)-QQ(IAT,JO,JO))**2.0D+00
          END DO
          PREC(IO,JO) = 1.0D+00/TMP
          PREC(JO,IO) = PREC(IO,JO)
        END DO
      END DO
C
      RETURN
C
      END SUBROUTINE DFTB_PRECPL
C
C-----------------------------------------------------------------------
C*MODULE FMOAFO    *DECK DFTB_CNSTA
C>
C>    @brief   FMO(-DFTB)/AFO for analytic derivative
C>
C>    @details Compute temporary A matrix (Eq. (81)) and the Lagrangian
C>             of the first term in Eq. (88). Brute-force
C>             implementation. There must be a better approach.
C>
C>    @author Yoshio Nishimoto
C>    - Feb, 2018- Subroutine written
C>
C>           --- INPUT ---
C>    @param NAT      Number of atoms
C>    @param NA       Number of alpha electrons
C>    @param L1       Number of AOs
C>    @param NLOC     Number of localized orbitals
C>    @param VALM     M_{\mu i}^{A,\Xi} (Eq. (84), as a work space)
C>    @param QQ       Q_{ij}^{A,\Xi}, Eq. (85)
C>    @param VECL     LMO coefficient
C>    @param IND      AO index of each atom
C>    @param ZZ       \mathbf{Z}^{\mathrm{V},\Xi} (solution of Eq. (87))
C>    @param WRK1     Work space
C>           --- OUTPUT ---
C>    @param AA       A matrix (\widetilde{A}_{ij,mn}^\Xi; Eq. (81))
C>    @param WW       Energy-weighted density matrix-like contribution
C>                    (the first term in Eq. (88))
C>
C
      SUBROUTINE DFTB_CNSTA(NAT,NA,L1,NLOC,AA,VALM,QQ,VECL,IND,
     *  ZZ,WRK1,WW)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER (ZERO=0.0D+00,ONE=1.0D+00)
C
      DIMENSION AA(L1,L1),VALM(NAT,L1,NLOC),QQ(NAT,NLOC,NLOC),
     *          VECL(L1,*),IND(*)
      DIMENSION ZZ(NLOC,*),WRK1(*),WW(*)
C
      DO IAT = 1, NAT
        DO MU = 1, L1
          DO IO = 1, NLOC
            TMP = ZERO
            DO MM0 = 1, IND(IAT+1)-IND(IAT)
              MM=IND(IAT)+MM0
              IF (MU.EQ.MM) THEN
                TMP = TMP + VECL(MM,IO)
              END IF
            END DO
            DO KO = 1, NA
              TMP = TMP - 0.5D+00*VECL(MU,KO)*QQ(IAT,KO,IO)
            END DO
            VALM(IAT,MU,IO) = TMP
          END DO
        END DO
      END DO
C
      DO MU = 1, L1
      DO NU = 1, L1
        TMP0 = ZERO
        DO IO = 1, NLOC
        DO JO = 1, IO-1 !NLOC
          TMP = ZERO
          DO IAT = 1, NAT
            TMP = TMP
     *        + VALM(IAT,MU,IO)*VECL(NU,JO)*QQ(IAT,IO,IO)
     *        + VALM(IAT,MU,JO)*VECL(NU,IO)*QQ(IAT,IO,IO)
     *        + 2.0D+00*VALM(IAT,MU,IO)*VECL(NU,IO)*QQ(IAT,IO,JO)
     *        - VALM(IAT,MU,JO)*VECL(NU,IO)*QQ(IAT,JO,JO)
     *        - VALM(IAT,MU,IO)*VECL(NU,JO)*QQ(IAT,JO,JO)
     *        - 2.0D+00*VALM(IAT,MU,JO)*VECL(NU,JO)*QQ(IAT,JO,IO)
          END DO
          TMP0 = TMP0 + TMP*ZZ(IO,JO)
        END DO
        END DO
        AA(MU,NU) = TMP0
      END DO
      END DO
C
      CALL SYMSQT(L1,AA,L1,WRK1,1)
      L2=L1*(L1+1)/2
      CALL DAXPY(L2,-ONE,WRK1,1,WW,1)
C
      RETURN
C
      END SUBROUTINE DFTB_CNSTA
C
C-----------------------------------------------------------------------
C*MODULE FMOAFO    *DECK DFTB_CNSTB
C>
C>    @brief   FMO(-DFTB)/AFO for analytic derivative
C>
C>    @details Compute B matrix (Eq. (82)) and the second term in
C>             Eq. (92). Brute-force implementation. There must be a
C>             better approach.
C>
C>    @author Yoshio Nishimoto
C>    - Feb, 2018- Subroutine written
C>
C>           --- INPUT ---
C>    @param NAT      Number of atoms
C>    @param NA       Number of alpha electrons
C>    @param L1       Number of AOs
C>    @param NLOC     Number of localized orbitals
C>    @param QQ       Q_{ij}^{A,\Xi}, Eq. (85)
C>    @param QBAR     \bar{Q}_{ij}^{A,\Xi}, Eq. (86)
C>    @param WW       W_{mi}^\Xi, see Eq. (79)
C>    @param ZZ       \mathbf{Z}^{\mathrm{V},\Xi} (solution of Eq. (87))
C>           --- OUTPUT ---
C>    @param BB       B matrix (\widetilde{B}_{ij,kl}^\Xi; Eq. (82))
C>    @param VALLAG   Lagrangian contribution (R.H.S. of Eq. (92))
C>
C
      SUBROUTINE DFTB_CNSTB(NAT,NA,L1,NLOC,BB,QQ,QBAR,WW,VALLAG,ZZ)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION BB(NLOC,NLOC,L1,L1),QQ(NAT,NLOC,NLOC),QBAR(NAT,L1,NLOC),
     *          WW(NLOC,NLOC),VALLAG(L1-NA,NA),ZZ(NLOC,NLOC)
C
      DO IO = 1, NLOC
      DO JO = 1, NLOC
      DO KO = NA+1, L1
      DO LO = 1, NA
        TMP = 0.0D+00 !ZERO
        DO IAT = 1, NAT
          TMP = TMP
     *      + WW(LO,IO)*QBAR(IAT,KO,JO)*(QQ(IAT,IO,IO)-QQ(IAT,JO,JO))
     *      + 2.0D+00*QQ(IAT,IO,JO)*WW(LO,IO)*QBAR(IAT,KO,IO)
        END DO
        BB(IO,JO,KO,LO) = TMP
      END DO
      END DO
      END DO
      END DO
C
      DO MO0 = 1, L1-NA
        MO = MO0 + NA
        DO NO = 1, NA
          VALLAG(MO0,NO) = VALLAG(MO0,NO)
     *      - DDOT(NLOC*NLOC,BB(1,1,MO,NO),1,ZZ,1)
        END DO
      END DO
C
      END SUBROUTINE DFTB_CNSTB
C
C-----------------------------------------------------------------------
C*MODULE FMOAFO    *DECK DFTB_CCC
C>
C>    @brief   FMO(-DFTB)/AFO for analytic derivative
C>
C>    @details Calculate Z*C in solving Z-vector equations for localized
C>             orbitals, Eq. (87). This subroutine may actually computes
C>             C*Z, but C is symmetric.
C>
C>    @author Yoshio Nishimoto
C>    - Feb, 2018- Subroutine written
C>
C>           --- INPUT ---
C>    @param NLOC     Number of localized orbitals
C>    @param C        C matrix constructed in DFTB_CNSTC
C>    @param P        Residual vector
C>           --- OUTPUT ---
C>    @param AP       C*P
C>
C
      SUBROUTINE DFTB_CCC(NLOC,C,P,AP)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION C(NLOC,NLOC,*),P(NLOC,NLOC),AP(NLOC,NLOC)
C
      DO IO = 1, NLOC
        DO JO = 1, NLOC
          CALL DFTB_CNVSQ(IO,JO,NSEQ1)
          TMP=0.0D+00
          DO KO = 1, NLOC
            DO LO = 1, KO-1
              CALL DFTB_CNVSQ(KO,LO,NSEQ2)
              TMP = TMP + C(IO,JO,NSEQ2)*P(KO,LO)
            END DO
          END DO
          AP(IO,JO) = TMP
        END DO
      END DO
C
      RETURN
C
      END SUBROUTINE DFTB_CCC
