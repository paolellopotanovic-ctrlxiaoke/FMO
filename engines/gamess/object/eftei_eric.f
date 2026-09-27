C*MODULE EFTEI_ERIC   *DECK ERIC_QMEFP
C>
C>    @brief   ERIC: Electron Repulsion Integral Calculator for QMEFP
C>
C>    @author  Peng Xu and Tosaporn Sattasathuchana 
C>             Jan 2020
C>    @details This subroutine prepares the necessary shell information
C>             for QM/EFP ERI (ERIC) calculation and 
C>             is based on the subroutine ERIC in int2c.src
C>
C>    @param   JM: index for the type of EFP potentials 
C>    @param   JNAT: counter for number of EFP atoms 
C>    @param   MODE: different type of TEI involve in QM/EFP EXREP
C>                 MODE=2 two-center on QM, two-center on EFP
C>                 MODE=3 three-center on QM, One-center on EFP
C>
      SUBROUTINE ERIC_QMEFP(ISH,JSH,KSH,LSH,ERI,JM,JNAT,MODE)
      use mx_limits, only: mxatm,mxgtot,mxsh,mxgsh,mxg2,MXGEFP
      USE comm_EFPBAS
C ----------------------------------------------------------------------
C
C              ELECTRON REPULSION INTEGRAL CALCULATOR
C
C ----------------------------------------------------------------------
C
C  METHOD:
C     A) FORM SCALED, CONTRACTED 1-CENTER PRECURSOR INTEGRALS.
C        CONVERT THESE TO 4-CENTER INTEGRALS OVER CARTESIAN
C        GAUSSIANS USING,
C     B) PRECURSOR-HERMITE TRANSFER EQUATION (PTE)
C     C) CONTRACTED TRANSFER  EQUATION (CTE)
C     D) HORIZONTAL RECURSION RELATION (HRR)
C
C
C  REFERENCES
C
C     PTE:
C     "Recursion Formula for Electron Repulsion Integrals Over
C     Hermite Polynomials"
C         G.D.Fletcher  Int.J.Quantum Chem. 106, 355-360(2006)
C
C     HRR:
C     M. HEAD-GORDON & J. A. POPLE,
C     J. CHEM. PHYS., 89, 5777-5786 (1988).
C
C     CTE:
C     P. M. W. GILL, M. HEAD-GORDON, & J. A. POPLE,
C     INT. J. Q. CHEM., SYMP. 23, 269-280 (1989).
C
C     INTERPOLATION METHOD:
C     P. M. W. GILL, B. G. JOHNSON, & J. A. POPLE,
C     INT. J. Q. CHEM., 40, 745-752 (1991).
C
C
C  THIS IS THE DRIVER ROUTINE THAT INTERFACES TO GAMESS COMMON NSHEL
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      DOUBLE PRECISION ERI(*)
C
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /NSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),
     *                CF(MXGTOT),CG(MXGTOT),CH(MXGTOT),CI(MXGTOT),
     *                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),
     *                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
C
C  /ERIPRM/ SYMBOLS:
C     EXI,J,K,L   = PRIMITIVE EXPONENTS OF I,J,K,L SHELLS
C     CCI,J,K,L   = CONTRACTION COEFFICIENTS OF ISH,JSH,KSH,LSH
C     XAB,YAB,ZAB = EXPONENT-WEIGHTED DISTANCE FOR BRA
C     XCD,YCD,ZCD = EXPONENT-WEIGHTED DISTANCE FOR KET
C     CCBRA,CCKET = PRIMITIVE CHARGE-CLOUD FACTOR FOR BRA,KET
C     SLBRA,SLKET = FACTOR CONVERTING SMALL-T EXPRESSION TO LARGE-T
C
      COMMON /ERIPRM/ EXI(MXGSH),EXJ(MXGSH),EXK(MXGSH),EXL(MXGSH),
     *                CCI(MXGSH),CCJ(MXGSH),CCK(MXGSH),CCL(MXGSH),
     *                XAB(MXG2),YAB(MXG2),ZAB(MXG2),
     *                XCD(MXG2),YCD(MXG2),ZCD(MXG2),
     *                CCBRA(MXG2),CCKET(MXG2),RXB(MXG2),
     *                SLBRA(MXG2),SLKET(MXG2),RXK(MXG2)
C
      COMMON /ERIDAT/ LEN1,LEN2,LEN3,LEN4
      COMMON /ERIOUT/ INW,JNW,KNW,LNW,LSTRI,LSTRJ,LSTRK,LSTRL
C
      COMMON /SHLEXC/ NORGSH(3),NORGSP(3),IEXCH,NANGM,NGTH(4)
      COMMON /SHLNOS/ QQ4,LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,
     2                MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL,
     3                NIJ,IJ,KL,IJKL

      COMMON /FMCOM /XX(1)
C
      LOGICAL IEQJ,KEQL
      LOGICAL TEST1,TEST2,TEST3
C
      PARAMETER (PT5=0.5D+00)
      PARAMETER (ONE=1.0D+00)
      PARAMETER (TWO=2.0D+00)
C
C     SR3 <==> SQRT( 3)   <==> 1.7320508075688773D+00
C     SR5 <==> SQRT( 5)   <==> 2.2360679774997897D+00
C     SR7 <==> SQRT( 7)   <==> 2.6457513110645906D+00
C     S15 <==> SQRT(15)   <==> 3.8729833462074169D+00
C     S35 <==> SQRT(35)   <==> 5.9160797830996160D+00
C     S53 <==> SQRT(35/3) <==> 3.4156502553198661D+00
C
      PARAMETER (SR3=1.7320508075688773D+00)
      PARAMETER (SR5=2.2360679774997897D+00)
      PARAMETER (SR7=2.6457513110645906D+00)
      PARAMETER (S15=SR3*SR5)
      PARAMETER (S35=SR5*SR7)
      PARAMETER (S53=S35/SR3)
C
      PARAMETER (PI254 =5.9149671727956129D+00)
      PARAMETER (PI214 =0.94139626377671481D+00)
C
      DIMENSION IORD(35),ANGL(35)
      DATA IORD/
     1       1,
     2       2,  3,  4,
     3       5,  7, 10,  6,  8,  9,
     4      11, 14, 20, 12, 15, 13, 17, 18, 19, 16,
     5      21, 25, 35, 22, 26, 24, 29, 33, 34, 23, 30, 32, 27, 28, 31/
      DATA ANGL/
     1     ONE,
     2     ONE,ONE,ONE,
     3     ONE,SR3,ONE,SR3,SR3,ONE,
     4     ONE,SR5,SR5,ONE,SR5,S15,SR5,SR5,SR5,ONE,
     5     ONE,SR7,S53,SR7,ONE,SR7,S35,S35,SR7,S53,S35,S53,SR7,SR7,ONE/
C
C  FAST CODES
C
C  RE-ORDER CHARGE-CLOUD SHELLS
C  (HRR DATA ONLY STORED FOR LI>=LJ)
C
      INW = ISH
      JNW = JSH
      KNW = KSH
      LNW = LSH
C
      LSTRI = NGTH(1)
      LSTRJ = NGTH(2)
      LSTRK = NGTH(3)
      LSTRL = NGTH(4)
C
      LANGI = KTYPEF(INW,JM) - 1
      LANGJ = KTYPE(JNW) - 1
      IF(MODE.EQ.2) THEN
      LANGK = KTYPEF(KNW,JM) - 1
      ELSE
      LANGK = KTYPE(KNW) - 1
      ENDIF
      LANGL = KTYPE(LNW) - 1

      LBRA  = LANGI + LANGJ
      LKET  = LANGK + LANGL

      TEST1 = LBRA.GT.LKET
      TEST2 = LBRA.EQ.LKET

      IF((LANGI.LT.LANGJ).AND.(LANGK.LT.LANGL)) THEN

        INW = JSH
        JNW = ISH
        KNW = LSH
        LNW = KSH

        LANGI = KTYPE(INW)
        LANGJ = KTYPEF(JNW,JM)
        LANGK = KTYPE(KNW)
        IF(MODE.EQ.2) THEN
        LANGL = KTYPEF(LNW,JM) 
        ELSE
        LANGL = KTYPE(LNW)
        ENDIF

        LNGIJ = (LANGI*LANGI-LANGI)/2 + LANGJ
        LNGKL = (LANGK*LANGK-LANGK)/2 + LANGL

        TEST3 = LNGKL.GT.LNGIJ

        IF(TEST1) THEN
          CALL ERIC_QMEFP_CASE_8(MODE,JM,JNAT,ISH,JSH,KSH,LSH,
     &         LANGI,LANGJ,LANGK,LANGL,IPRIM,JPRIM,KPRIM,LPRIM,
     &         IEQJ,KEQL,XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD,
     &         MINI,MAXI,MINJ,MAXJ,MINK,MAXK,MINL,MAXL)
        ELSEIF(TEST2.and.TEST3) THEN
          CALL ERIC_QMEFP_CASE_8(MODE,JM,JNAT,ISH,JSH,KSH,LSH,
     &         LANGI,LANGJ,LANGK,LANGL,IPRIM,JPRIM,KPRIM,LPRIM,
     &         IEQJ,KEQL,XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD,
     &         MINI,MAXI,MINJ,MAXJ,MINK,MAXK,MINL,MAXL)
        ELSE
          CALL ERIC_QMEFP_CASE_4(MODE,JM,JNAT,ISH,JSH,KSH,LSH,
     &         LANGI,LANGJ,LANGK,LANGL,IPRIM,JPRIM,KPRIM,LPRIM,
     &         IEQJ,KEQL,XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD,
     &         MINI,MAXI,MINJ,MAXJ,MINK,MAXK,MINL,MAXL)
        ENDIF          

      ELSEIF((LANGI.LT.LANGJ).AND.(LANGK.GE.LANGL)) THEN

        INW = JSH
        JNW = ISH
        KNW = KSH
        LNW = LSH

        LANGI = KTYPE(INW)
        LANGJ = KTYPEF(JNW,JM)
        IF(MODE.EQ.2) THEN
        LANGK = KTYPEF(KNW,JM)
        ELSE
        LANGK = KTYPE(KNW)
        ENDIF
        LANGL = KTYPE(LNW)

        LNGIJ = (LANGI*LANGI-LANGI)/2 + LANGJ
        LNGKL = (LANGK*LANGK-LANGK)/2 + LANGL

        TEST3 = LNGKL.GT.LNGIJ

        IF(TEST1) THEN
          CALL ERIC_QMEFP_CASE_6(MODE,JM,JNAT,ISH,JSH,KSH,LSH,
     &         LANGI,LANGJ,LANGK,LANGL,IPRIM,JPRIM,KPRIM,LPRIM,
     &         IEQJ,KEQL,XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD,
     &         MINI,MAXI,MINJ,MAXJ,MINK,MAXK,MINL,MAXL)
        ELSEIF(TEST2.and.TEST3) THEN
          CALL ERIC_QMEFP_CASE_6(MODE,JM,JNAT,ISH,JSH,KSH,LSH,
     &         LANGI,LANGJ,LANGK,LANGL,IPRIM,JPRIM,KPRIM,LPRIM,
     &         IEQJ,KEQL,XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD,
     &         MINI,MAXI,MINJ,MAXJ,MINK,MAXK,MINL,MAXL)
        ELSE
          CALL ERIC_QMEFP_CASE_2(MODE,JM,JNAT,ISH,JSH,KSH,LSH,
     &         LANGI,LANGJ,LANGK,LANGL,IPRIM,JPRIM,KPRIM,LPRIM,
     &         IEQJ,KEQL,XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD,
     &         MINI,MAXI,MINJ,MAXJ,MINK,MAXK,MINL,MAXL)
        ENDIF

      ELSEIF((LANGI.GE.LANGJ).AND.(LANGK.LT.LANGL)) THEN

        INW = ISH
        JNW = JSH
        KNW = LSH
        LNW = KSH
        
        LANGI = KTYPEF(INW,JM)
        LANGJ = KTYPE(JNW)
        LANGK = KTYPE(KNW)
        IF(MODE.EQ.2) THEN
        LANGL = KTYPEF(LNW,JM)
        ELSE
        LANGL = KTYPE(LNW)
        ENDIF

        LNGIJ = (LANGI*LANGI-LANGI)/2 + LANGJ
        LNGKL = (LANGK*LANGK-LANGK)/2 + LANGL

        TEST3 = LNGKL.GT.LNGIJ

        IF(TEST1) THEN
          CALL ERIC_QMEFP_CASE_7(MODE,JM,JNAT,ISH,JSH,KSH,LSH,
     &         LANGI,LANGJ,LANGK,LANGL,IPRIM,JPRIM,KPRIM,LPRIM,
     &         IEQJ,KEQL,XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD,
     &         MINI,MAXI,MINJ,MAXJ,MINK,MAXK,MINL,MAXL)
        ELSEIF(TEST2.and.TEST3) THEN
          CALL ERIC_QMEFP_CASE_7(MODE,JM,JNAT,ISH,JSH,KSH,LSH,
     &         LANGI,LANGJ,LANGK,LANGL,IPRIM,JPRIM,KPRIM,LPRIM,
     &         IEQJ,KEQL,XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD,
     &         MINI,MAXI,MINJ,MAXJ,MINK,MAXK,MINL,MAXL)
        ELSE
          CALL ERIC_QMEFP_CASE_3(MODE,JM,JNAT,ISH,JSH,KSH,LSH,
     &         LANGI,LANGJ,LANGK,LANGL,IPRIM,JPRIM,KPRIM,LPRIM,
     &         IEQJ,KEQL,XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD,
     &         MINI,MAXI,MINJ,MAXJ,MINK,MAXK,MINL,MAXL)
        ENDIF

      ELSEIF((LANGI.GE.LANGJ).AND.(LANGK.GE.LANGL)) THEN
       
        LANGI = KTYPEF(INW,JM)
        LANGJ = KTYPE(JNW)
        IF(MODE.EQ.2) THEN
        LANGK = KTYPEF(KNW,JM)
        ELSE
        LANGK = KTYPE(KNW)
        ENDIF
        LANGL = KTYPE(LNW)

        LNGIJ = (LANGI*LANGI-LANGI)/2 + LANGJ
        LNGKL = (LANGK*LANGK-LANGK)/2 + LANGL

        TEST3 = LNGKL.GT.LNGIJ

        IF(TEST1) THEN
          CALL ERIC_QMEFP_CASE_5(MODE,JM,JNAT,ISH,JSH,KSH,LSH,
     &         LANGI,LANGJ,LANGK,LANGL,IPRIM,JPRIM,KPRIM,LPRIM,
     &         IEQJ,KEQL,XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD,
     &         MINI,MAXI,MINJ,MAXJ,MINK,MAXK,MINL,MAXL)
        ELSEIF(TEST2.and.TEST3) THEN
          CALL ERIC_QMEFP_CASE_5(MODE,JM,JNAT,ISH,JSH,KSH,LSH,
     &         LANGI,LANGJ,LANGK,LANGL,IPRIM,JPRIM,KPRIM,LPRIM,
     &         IEQJ,KEQL,XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD,
     &         MINI,MAXI,MINJ,MAXJ,MINK,MAXK,MINL,MAXL)
        ELSE
          CALL ERIC_QMEFP_CASE_1(MODE,JM,JNAT,ISH,JSH,KSH,LSH,
     &         LANGI,LANGJ,LANGK,LANGL,IPRIM,JPRIM,KPRIM,LPRIM,
     &         IEQJ,KEQL,XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD,
     &         MINI,MAXI,MINJ,MAXJ,MINK,MAXK,MINL,MAXL)
        ENDIF

      ELSE
        write(6,*) 'SHOULD NOT GET HERE!'
      ENDIF

      RAB = (XA-XB)**2 + (YA-YB)**2 + (ZA-ZB)**2
      IJ = 0
      DO JJ = 1, JPRIM
         EJ = EXJ(JJ)
         CJ = CCJ(JJ)
         ITOP = IPRIM
         IF(IEQJ) ITOP = JJ
         DO II = 1, ITOP
            EI = EXI(II)
            CIX = CCI(II)
            EIJ = ONE/(EI+EJ)
            IJ = IJ + 1
            CCFAC = PI254*CIX*CJ*EIJ*EXP( -EI*EJ*RAB*EIJ )
            IF(IEQJ .AND. II.NE.JJ) CCFAC = CCFAC*TWO
            CCBRA(IJ) = CCFAC
            SLBRA(IJ) = PI214*SQRT(EIJ)
            XAB(IJ) = (EI*XA + EJ*XB)*EIJ
            YAB(IJ) = (EI*YA + EJ*YB)*EIJ
            ZAB(IJ) = (EI*ZA + EJ*ZB)*EIJ
            RXB(IJ) = EIJ*PT5
         END DO
      END DO

      RCD = (XC-XD)**2 + (YC-YD)**2 + (ZC-ZD)**2
      KL = 0
      DO LL = 1, LPRIM
         EL = EXL(LL)
         CL = CCL(LL)
         KTOP = KPRIM
         IF(KEQL) KTOP = LL
         DO KK = 1, KTOP
            EK = EXK(KK)
            CK = CCK(KK)
            EKL = ONE/(EK+EL)
            KL = KL + 1
            CCFAC = PI254*CK*CL*EKL*EXP( -EK*EL*RCD*EKL )
            IF(KEQL .AND. KK.NE.LL) CCFAC = CCFAC*TWO
            CCKET(KL) = CCFAC
            SLKET(KL) = PI214*SQRT(EKL)
            XCD(KL) = (EK*XC + EL*XD)*EKL
            YCD(KL) = (EK*YC + EL*YD)*EKL
            ZCD(KL) = (EK*ZC + EL*ZD)*EKL
            RXK(KL) = EKL*PT5
         END DO
      END DO
C
      CALL VALFM(LOADFM)
      IPHI = LOADFM + 1
C
C  ANGULAR MOMENTUM 4-INDEX
C
      LBGT  = MAX(LANGI,LANGJ)
      LBLT  = MIN(LANGI,LANGJ)
      LNGIJ = (LBGT*LBGT+LBGT)/2 + LBLT
      LKGT  = MAX(LANGK,LANGL)
      LKLT  = MIN(LANGK,LANGL)
      LNGKL = (LKGT*LKGT+LKGT)/2 + LKLT
      LQGT  = MAX(LNGKL,LNGIJ)
      LQLT  = MIN(LNGKL,LNGIJ)
      LIJKL = (LQGT*LQGT+LQGT)/2 + LQLT
C
      IF(LIJKL.LE.5) THEN
C
C  SP CASES
C
         IF(LIJKL.EQ.0) THEN
C  SSSS
            IWK1 = IPHI +    1     ! LPHI
            IWK2 = IWK1 +    1     ! LWK1
            LAST = IWK2 +    1     ! LWK1
            NEED = LAST - LOADFM
            CALL GETFM(NEED)
            IDIM =           1
            CALL SSSS (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL,
     *                 XX(IPHI),XX(IWK2),IDIM)
            IOFF = IWK2 +    0
         ELSE IF(LIJKL.EQ.1) THEN
C  PSSS
            IWK1 = IPHI +   21
            IWK2 = IWK1 +    1
            LAST = IWK2 +    4
            NEED = LAST - LOADFM
            CALL GETFM(NEED)
            IDIM =           4
            CALL PSSS (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL,
     *                 XC,YC,ZC,XD,YD,ZD,
     *                 XX(IPHI),XX(IWK2),IDIM)
            IOFF = IWK2 +    1
         ELSE IF(LIJKL.EQ.2) THEN
C  PSPS
            IWK1 = IPHI +   65
            IWK2 = IWK1 +    4
            LAST = IWK2 +   12
            NEED = LAST - LOADFM
            CALL GETFM(NEED)
            IDIM =           4
            CALL PSPS (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL
     *,                XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD
     *,                XX(IPHI),XX(IWK1),XX(IWK2),IDIM)
            IOFF = IWK2 +    1
         ELSE IF(LIJKL.EQ.3) THEN
C  PPSS
            IWK1 = IPHI +   65
            IWK2 = IWK1 +    1
            LAST = IWK2 +   24
            NEED = LAST - LOADFM
            CALL GETFM(NEED)
            IDIM =          24
            CALL PPSS (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL,
     *                 XC,YC,ZC,XD,YD,ZD,
     *                 XX(IPHI),XX(IWK2),IDIM)
            IOFF = IWK2 +   15
         ELSE IF(LIJKL.EQ.4) THEN
C  PPPS
            IWK1 = IPHI +  183
            IWK2 = IWK1 +    4
            LAST = IWK2 +   72
            NEED = LAST - LOADFM
            CALL GETFM(NEED)
            IDIM =          24
            CALL PPPS (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL
     *,                XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD
     *,                XX(IPHI),XX(IWK1),XX(IWK2),IDIM)
            IOFF = IWK2 +   15
         ELSE IF(LIJKL.EQ.5) THEN
C  PPPP
            IWK1 = IPHI +  574
            IWK2 = IWK1 +   24
            LAST = IWK2 +  216
            NEED = LAST - LOADFM
            CALL GETFM(NEED)
            IDIM =          24
            CALL PPPP (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL
     *,                XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD
     *,                XX(IPHI),XX(IWK1),XX(IWK2),IDIM)
            IOFF = IWK2 +   15
         END IF
C
C  END OF SP CASES
C
      ELSE IF(LIJKL.GE.6 .AND. LIJKL.LE.20) THEN
C
C  D CASES
C
         IF(LIJKL.EQ.6) THEN
C  DSSS
            IWK1 = IPHI +   58
            IWK2 = IWK1 +    1
            LAST = IWK2 +   11
            NEED = LAST - LOADFM
            CALL GETFM(NEED)
            IDIM =          11
            CALL DSSS (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL,
     *                 XC,YC,ZC,XD,YD,ZD,
     *                 XX(IPHI),XX(IWK2),IDIM)
            IOFF = IWK2 +    5
         ELSE IF(LIJKL.EQ.7) THEN
C  DSPS
            IWK1 = IPHI +  159
            IWK2 = IWK1 +    4
            LAST = IWK2 +   33
            NEED = LAST - LOADFM
            CALL GETFM(NEED)
            IDIM =          11
            CALL DSPS (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL
     *,                XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD
     *,                XX(IPHI),XX(IWK1),XX(IWK2),IDIM)
            IOFF = IWK2 +    5
         ELSE IF(LIJKL.EQ.8) THEN
C  DSPP
            IWK1 = IPHI +  453
            IWK2 = IWK1 +   11
            LAST = IWK2 +  144
            NEED = LAST - LOADFM
            CALL GETFM(NEED)
            IDIM =          24
            CALL DSPP (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL
     *,                XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD
     *,                XX(IPHI),XX(IWK1),XX(IWK2),IDIM)
            IOFF = IWK2 +   15
         ELSE IF(LIJKL.EQ.9) THEN
C  DSDS
            IWK1 = IPHI +  389
            IWK2 = IWK1 +   11
            LAST = IWK2 +   66
            NEED = LAST - LOADFM
            CALL GETFM(NEED)
            IDIM =          11
            CALL DSDS (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL
     *,                XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD
     *,                XX(IPHI),XX(IWK1),XX(IWK2),IDIM)
            IOFF = IWK2 +    5
         ELSE IF(LIJKL.EQ.10) THEN
C  DPSS
            IWK1 = IPHI +  149
            IWK2 = IWK1 +    1
            LAST = IWK2 +   53
            NEED = LAST - LOADFM
            CALL GETFM(NEED)
            IDIM =          53
            CALL DPSS (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL,
     *                 XC,YC,ZC,XD,YD,ZD,
     *                 XX(IPHI),XX(IWK2),IDIM)
            IOFF = IWK2 +   35
         ELSE IF(LIJKL.EQ.11) THEN
C  DPPS
            IWK1 = IPHI +  390
            IWK2 = IWK1 +    4
            LAST = IWK2 +  159
            NEED = LAST - LOADFM
            CALL GETFM(NEED)
            IDIM =          53
            CALL DPPS (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL
     *,                XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD
     *,                XX(IPHI),XX(IWK1),XX(IWK2),IDIM)
            IOFF = IWK2 +   35
         ELSE IF(LIJKL.EQ.12) THEN
C  DPPP
            IWK1 = IPHI + 1181
            IWK2 = IWK1 +   24
            LAST = IWK2 +  477
            NEED = LAST - LOADFM
            CALL GETFM(NEED)
            IDIM =          53
            CALL DPPP (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL
     *,                XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD
     *,                XX(IPHI),XX(IWK1),XX(IWK2),IDIM)
            IOFF = IWK2 +   35
         ELSE IF(LIJKL.EQ.13) THEN
C  DPDS
            IWK1 = IPHI +  928
            IWK2 = IWK1 +   11
            LAST = IWK2 +  318
            NEED = LAST - LOADFM
            CALL GETFM(NEED)
            IDIM =          53
            CALL DPDS (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL
     *,                XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD
     *,                XX(IPHI),XX(IWK1),XX(IWK2),IDIM)
            IOFF = IWK2 +   35
         ELSE IF(LIJKL.EQ.15) THEN
C  DDSS
            IWK1 = IPHI +  319
            IWK2 = IWK1 +    1
            LAST = IWK2 +  165
            NEED = LAST - LOADFM
            CALL GETFM(NEED)
            IDIM =         165
            CALL DDSS (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL,
     *                 XC,YC,ZC,XD,YD,ZD,
     *                 XX(IPHI),XX(IWK2),IDIM)
            IOFF = IWK2 +  129
         ELSE IF(LIJKL.EQ.16) THEN
C  DDPS
            IWK1 = IPHI +  802
            IWK2 = IWK1 +    4
            LAST = IWK2 +  495
            NEED = LAST - LOADFM
            CALL GETFM(NEED)
            IDIM =         165
            CALL DDPS (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL
     *,                XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD
     *,                XX(IPHI),XX(IWK1),XX(IWK2),IDIM)
            IOFF = IWK2 +  129
         END IF
C
C  END OF D CASES
C
      ELSE IF(LIJKL.GE.21 .AND. LIJKL.LE.54) THEN
C
C  F CASES
C
         IF(LIJKL.EQ.21) THEN
C  FSSS
            IWK1 = IPHI +  129
            IWK2 = IWK1 +    1
            LAST = IWK2 +   24
            NEED = LAST - LOADFM
            CALL GETFM(NEED)
            IDIM =          24
            CALL FSSS (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL,
     *                 XC,YC,ZC,XD,YD,ZD,
     *                 XX(IPHI),XX(IWK2),IDIM)
            IOFF = IWK2 +   14
         ELSE IF(LIJKL.EQ.22) THEN
C  FSPS
            IWK1 = IPHI +  328
            IWK2 = IWK1 +    4
            LAST = IWK2 +   72
            NEED = LAST - LOADFM
            CALL GETFM(NEED)
            IDIM =          24
            CALL FSPS (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL
     *,                XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD
     *,                XX(IPHI),XX(IWK1),XX(IWK2),IDIM)
            IOFF = IWK2 +   14
         ELSE IF(LIJKL.EQ.23) THEN
C  FSPP
            IWK1 = IPHI +  976
            IWK2 = IWK1 +   24
            LAST = IWK2 +  216
            NEED = LAST - LOADFM
            CALL GETFM(NEED)
            IDIM =          24
            CALL FSPP (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL
     *,                XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD
     *,                XX(IPHI),XX(IWK1),XX(IWK2),IDIM)
            IOFF = IWK2 +   14
         ELSE IF(LIJKL.EQ.24) THEN
C  FSDS
            IWK1 = IPHI +  768
            IWK2 = IWK1 +   11
            LAST = IWK2 +  144
            NEED = LAST - LOADFM
            CALL GETFM(NEED)
            IDIM =          24
            CALL FSDS (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL
     *,                XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD
     *,                XX(IPHI),XX(IWK1),XX(IWK2),IDIM)
            IOFF = IWK2 +   14
         ELSE IF(LIJKL.EQ.28) THEN
C  FPSS
            IWK1 = IPHI +  299
            IWK2 = IWK1 +    1
            LAST = IWK2 +  100
            NEED = LAST - LOADFM
            CALL GETFM(NEED)
            IDIM =         100
            CALL FPSS (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL,
     *                 XC,YC,ZC,XD,YD,ZD,
     *                 XX(IPHI),XX(IWK2),IDIM)
            IOFF = IWK2 +   70
         ELSE IF(LIJKL.EQ.29) THEN
C  FPPS
            IWK1 = IPHI +  740
            IWK2 = IWK1 +    4
            LAST = IWK2 +  300
            NEED = LAST - LOADFM
            CALL GETFM(NEED)
            IDIM =         100
            CALL FPPS (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL
     *,                XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD
     *,                XX(IPHI),XX(IWK1),XX(IWK2),IDIM)
            IOFF = IWK2 +   70
         ELSE IF(LIJKL.EQ.36) THEN
C  FDSS
            IWK1 = IPHI +  589
            IWK2 = IWK1 +    1
            LAST = IWK2 +  285
            NEED = LAST - LOADFM
            CALL GETFM(NEED)
            IDIM =         285
            CALL FDSS (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL,
     *                 XC,YC,ZC,XD,YD,ZD,
     *                 XX(IPHI),XX(IWK2),IDIM)
            IOFF = IWK2 +  225
         END IF
C
C  END OF F CASES
C
      ELSE IF(LIJKL.GE.55 .AND. LIJKL.LE.119) THEN
C
C  G CASES
C
         IF(LIJKL.EQ.55) THEN
C  GSSS
            IWK1 = IPHI +  255
            IWK2 = IWK1 +    1
            LAST = IWK2 +   46
            NEED = LAST - LOADFM
            CALL GETFM(NEED)
            IDIM =          46
            CALL GSSS (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL,
     *                 XC,YC,ZC,XD,YD,ZD,
     *                 XX(IPHI),XX(IWK2),IDIM)
            IOFF = IWK2 +   31
         ELSE IF(LIJKL.EQ.56) THEN
C  GSPS
            IWK1 = IPHI +  612
            IWK2 = IWK1 +    4
            LAST = IWK2 +  138
            NEED = LAST - LOADFM
            CALL GETFM(NEED)
            IDIM =          46
            CALL GSPS (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL
     *,                XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD
     *,                XX(IPHI),XX(IWK1),XX(IWK2),IDIM)
            IOFF = IWK2 +   31
         ELSE IF(LIJKL.EQ.66) THEN
C  GPSS
            IWK1 = IPHI +  545
            IWK2 = IWK1 +    1
            LAST = IWK2 +  171
            NEED = LAST - LOADFM
            CALL GETFM(NEED)
            IDIM =         171
            CALL GPSS (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL,
     *                 XC,YC,ZC,XD,YD,ZD,
     *                 XX(IPHI),XX(IWK2),IDIM)
            IOFF = IWK2 +  126
         END IF
C
C  END OF G CASES
C
      END IF    ! LIJKL
C
C  ANGULAR NORMALIZATION AND SAVE TO OUTPUT ARRAY
C  WITH HONDO INDEXING AND REORDERING FOR GAMESS
C
      LENI = MAXI - MINI + 1
      LENK = MAXK - MINK + 1

      II = 1
      DO I = MINI, MAXI
         IO = IORD(I)
         AI = ANGL(IO)
         IJ = II
         DO J = MINJ, MAXJ
            JO = IORD(J)
            AIJ = ANGL(JO)*AI
            JC = ((JO-MINJ)*LENI + IO-MINI)*IDIM + IOFF 
            IJK = IJ
            DO K = MINK, MAXK
               KO = IORD(K)
               AIJK = ANGL(KO)*AIJ
               IJKL = IJK
               DO L = MINL, MAXL
                  LO = IORD(L)
                  AIJKL = ANGL(LO)*AIJK
                  IR = (LO-MINL)*LENK + KO-MINK
                  ERI(IJKL) = XX(JC+IR)*AIJKL
                  IJKL = IJKL + LSTRL
               END DO
               IJK = IJK + LSTRK
            END DO
            IJ = IJ + LSTRJ
         END DO
         II = II + LSTRI
      END DO
      CALL RETFM(NEED)
C
      INW = ISH
      JNW = JSH
      KNW = KSH
      LNW = LSH
      LSTRI = NGTH(1)
      LSTRJ = NGTH(2)
      LSTRK = NGTH(3)
      LSTRL = NGTH(4)

      RETURN
      END

C*MODULE EFTEI_ERIC  *DECK SETUP_SHELLEFP
C>    @brief  prepare EFP shells
C>    @author  Peng Xu and Tosaporn Sattasathuchana 
C>             Jan 2020
C>    @detail prepare EFP shells for ERIC_QMEFP subroutine
C>
C>    @param  INW: index of the shell
C>    @param  LANGI: type of basis function 
C>    @param  IPRIM: number of gaussian function in shell INW
C>    @param  EXI: exponenent of basis function
C>    @param  CCI: coefficient of basis function
C>    @param  XA,YA,ZA: coordinates of atoms associated with the shell 
C>    @param  JM: index for the type of EFP potentials
C>    @param  JNAT: counter for number of EFP atoms 
C>    @param  MINI: starting counter of shell I
C>    @param  MAXI: ending counter of shell I

      SUBROUTINE SETUP_SHELLEFP(INW,LANGI,IPRIM,
     &   EXI,CCI,XA,YA,ZA,MINI,MAXI,JM,JNAT)
C
      use mx_limits, only: mxatm,mxgtot,mxsh,mxgsh,mxg2
      USE comm_EFPBAS

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      DIMENSION EXI(MXGSH), CCI(MXGSH)

      LANGI = KTYPEF(INW,JM) - 1
      IPRIM = KNGEF(INW,JM)
      MINI = KMINEF(INW,JM)
      MAXI = KMAXEF(INW,JM)
C
      K1 = KSTREF(INW,JM)
      K2 = K1 + IPRIM - 1
      II = 0
      DO I = K1, K2
         II = II + 1
         EXI(II) = EXEF(I,JM)
      END DO
      ITYP = KTYPEF(INW,JM)
      IF(ITYP.EQ.1) THEN
         II = 0
         DO I = K1, K2
            II = II + 1
            CCI(II) = CSEF(I,JM)
         END DO
      ELSE IF(ITYP.EQ.2) THEN
         II = 0
         DO I = K1, K2
            II = II + 1
            CCI(II) = CPEF(I,JM)
         END DO
      ELSE IF(ITYP.EQ.3) THEN
         II = 0
         DO I = K1, K2
            II = II + 1
            CCI(II) = CDEF(I,JM)
         END DO
      ELSE IF(ITYP.EQ.4) THEN
         II = 0
         DO I = K1, K2
            II = II + 1
            CCI(II) = CFEF(I,JM)
         END DO
      ELSE IF(ITYP.EQ.5) THEN
         II = 0
         DO I = K1, K2
            II = II + 1
            CCI(II) = CGEF(I,JM)
         END DO
      END IF
C
C  2-INDEX PARAMETERS
C
      IATM = KATMEF(INW,JM)
      XA   = PRCORD(1,IATM+JNAT)
      YA   = PRCORD(2,IATM+JNAT)
      ZA   = PRCORD(3,IATM+JNAT)
      RETURN
      END

C*MODULE EFTEI_ERIC  *DECK SETUP_SHELLAI
C>    @brief  prepare QM shells
C>    @author  Peng Xu and Tosaporn Sattasathuchana 
C>             Jan 2020
C>    @detail prepare QM shells for ERIC_QMEFP subroutine
C>
C>    @param  JNW: index of the shell
C>    @param  LANGJ: type of basis function 
C>    @param  JPRIM: number of gaussian function in shell INW
C>    @param  EXJ: exponenent of basis function
C>    @param  CCJ: coefficient of basis function
C>    @param  XB,YB,ZB: coordinates of atoms associated with the shell 
C>    @param  MINJ: starting counter of shell J 
C>    @param  MAXJ: ending counter of shell J

      SUBROUTINE SETUP_SHELLAI(JNW,LANGJ,JPRIM,
     &   EXJ,CCJ,XB,YB,ZB,MINJ,MAXJ)
      use mx_limits, only: mxatm,mxgtot,mxsh,mxgsh,mxg2

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /NSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),
     *                CF(MXGTOT),CG(MXGTOT),CH(MXGTOT),CI(MXGTOT),
     *                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),
     *                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL

      DIMENSION EXJ(MXGSH), CCJ(MXGSH)

      LANGJ = KTYPE(JNW) - 1
      JPRIM = KNG(JNW)
      MINJ = KMIN(JNW)
      MAXJ = KMAX(JNW)
C
      K1 = KSTART(JNW)
      K2 = K1 + JPRIM - 1
      JJ = 0
      DO J = K1, K2
         JJ = JJ + 1
         EXJ(JJ) = EX(J)
      END DO
      JTYP = KTYPE(JNW)
      IF(JTYP.EQ.1) THEN
         JJ = 0
         DO J = K1, K2
            JJ = JJ + 1
            CCJ(JJ) = CS(J)
         END DO
      ELSE IF(JTYP.EQ.2) THEN
         JJ = 0
         DO J = K1, K2
            JJ = JJ + 1
            CCJ(JJ) = CP(J)
         END DO
      ELSE IF(JTYP.EQ.3) THEN
         JJ = 0
         DO J = K1, K2
            JJ = JJ + 1
            CCJ(JJ) = CD(J)
         END DO
      ELSE IF(JTYP.EQ.4) THEN
         JJ = 0
         DO J = K1, K2
            JJ = JJ + 1
            CCJ(JJ) = CF(J)
         END DO
      ELSE IF(JTYP.EQ.5) THEN
         JJ = 0
         DO J = K1, K2
            JJ = JJ + 1
            CCJ(JJ) = CG(J)
         END DO
      END IF
C
      JATM = KATOM(JNW)
      XB   = C(1,JATM)
      YB   = C(2,JATM)
      ZB   = C(3,JATM)
      RETURN
      END


C*MODULE EFTEI_ERIC  *DECK ERIC_QMEFP_CASE_1 
C>    @brief  ERIC scheme for case 1 (IJ|KL) 
C>    @author  Peng Xu and Tosaporn Sattasathuchana 
C>             Jan 2020
C>    @detail Case 1: no switching of I,J,K and L
C>    @param  MODE: different type of TEI involve in QM/EFP EXREP
C>                 MODE=2 two-center on QM, two-center on EFP
C>                 MODE=3 three-center on QM, One-center on EFP
C>    @param   JM: index for the type of EFP potentials
C>    @param   JNAT: counter for number of EFP atoms
C>    @param   ISH,JSH,KSH,LSH: counter of shells I,J,K and L
C>    @param   LANGI,LANGJ,LANGK,LANGL: type of basis function of shells 
C>    @param   IPRIM,JPRIM,KPRIM,LPRIM: number of gaussian function in
C>                                      shells
C>    @param   IEQJ,KEQL: logicals for pair IJ and KL, respectively 
C>    @param   XA,YA,ZA: xyz coordinate of atoms associated with shell I
C>    @param   XB,YB,ZB: xyz coordinate of atoms associated with shell J 
C>    @param   XC,YC,ZC: xyz coordinate of atoms associated with shell K
C>    @param   XD,YD,ZD: xyz coordinate of atoms associated with shell L
C>    @param   MINI,MINJ,MINK,MINL: starting counter of shells I,J,K & L 
C>    @param   MAXI,MAXJ,MAXK,MAXL: ending counter of shells I,J,K & L 

      SUBROUTINE ERIC_QMEFP_CASE_1(MODE,JM,JNAT,ISH,JSH,KSH,LSH,
     &                             LANGI,LANGJ,LANGK,LANGL,
     &                             IPRIM,JPRIM,KPRIM,LPRIM,
     &                             IEQJ,KEQL,XA,YA,ZA,XB,YB,ZB,
     &                             XC,YC,ZC,XD,YD,ZD,
     &                             MINI,MAXI,MINJ,MAXJ,MINK,MAXK,
     &                             MINL,MAXL)

      use mx_limits, only: mxatm,mxgtot,mxsh,mxgsh,mxg2,MXGEFP
      USE comm_EFPBAS

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON /ERIPRM/ EXI(MXGSH),EXJ(MXGSH),EXK(MXGSH),EXL(MXGSH),
     *                CCI(MXGSH),CCJ(MXGSH),CCK(MXGSH),CCL(MXGSH),
     *                XAB(MXG2),YAB(MXG2),ZAB(MXG2),
     *                XCD(MXG2),YCD(MXG2),ZCD(MXG2),
     *                CCBRA(MXG2),CCKET(MXG2),RXB(MXG2),
     *                SLBRA(MXG2),SLKET(MXG2),RXK(MXG2)
C
      COMMON /ERIDAT/ LEN1,LEN2,LEN3,LEN4
      COMMON /ERIOUT/ INW,JNW,KNW,LNW,LSTRI,LSTRJ,LSTRK,LSTRL
C
      COMMON /SHLEXC/ NORGSH(3),NORGSP(3),IEXCH,NANGM,NGTH(4)
      LOGICAL IEQJ,KEQL

      INW = ISH
      JNW = JSH
      KNW = KSH
      LNW = LSH
      LSTRI = NGTH(1)
      LSTRJ = NGTH(2)
      LSTRK = NGTH(3)
      LSTRL = NGTH(4)
       IEQJ = .false.
       KEQL = .false.

      CALL SETUP_SHELLEFP(INW,LANGI,IPRIM,EXI,CCI,XA,YA,ZA,MINI,MAXI,
     & JM,JNAT)

      CALL SETUP_SHELLAI(JNW,LANGJ,JPRIM,EXJ,CCJ,XB,YB,ZB,MINJ,MAXJ)

      IF(MODE.EQ.2) THEN
      CALL SETUP_SHELLEFP(KNW,LANGK,KPRIM,EXK,CCK,XC,YC,ZC,MINK,MAXK,
     & JM,JNAT)
      ELSE
      CALL SETUP_SHELLAI(KNW,LANGK,KPRIM,EXK,CCK,XC,YC,ZC,MINK,MAXK)
      KEQL = KNW.EQ.LNW
      ENDIF

      CALL SETUP_SHELLAI(LNW,LANGL,LPRIM,EXL,CCL,XD,YD,ZD,MINL,MAXL)

      RETURN
      END


C*MODULE EFTEI_ERIC  *DECK ERIC_QMEFP_CASE_2 
C>    @brief  ERIC scheme for case 2 (JI|KL) 
C>    @author  Peng Xu and Tosaporn Sattasathuchana 
C>             Jan 2020
C>    @detail Case 2: IJ switch
C>    @see    parameters same as ERIC_QMEFP_CASE_1 


      SUBROUTINE ERIC_QMEFP_CASE_2(MODE,JM,JNAT,ISH,JSH,KSH,LSH,
     &                             LANGI,LANGJ,LANGK,LANGL,
     &                             IPRIM,JPRIM,KPRIM,LPRIM,
     &                             IEQJ,KEQL,XA,YA,ZA,XB,YB,ZB,
     &                             XC,YC,ZC,XD,YD,ZD,
     &                             MINI,MAXI,MINJ,MAXJ,MINK,MAXK,
     &                             MINL,MAXL)

      use mx_limits, only: mxatm,mxgtot,mxsh,mxgsh,mxg2,MXGEFP
      USE comm_EFPBAS

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      COMMON /ERIPRM/ EXI(MXGSH),EXJ(MXGSH),EXK(MXGSH),EXL(MXGSH),
     *                CCI(MXGSH),CCJ(MXGSH),CCK(MXGSH),CCL(MXGSH),
     *                XAB(MXG2),YAB(MXG2),ZAB(MXG2),
     *                XCD(MXG2),YCD(MXG2),ZCD(MXG2),
     *                CCBRA(MXG2),CCKET(MXG2),RXB(MXG2),
     *                SLBRA(MXG2),SLKET(MXG2),RXK(MXG2)
C
      COMMON /ERIDAT/ LEN1,LEN2,LEN3,LEN4
      COMMON /ERIOUT/ INW,JNW,KNW,LNW,LSTRI,LSTRJ,LSTRK,LSTRL
C
      COMMON /SHLEXC/ NORGSH(3),NORGSP(3),IEXCH,NANGM,NGTH(4)
      LOGICAL IEQJ,KEQL

      INW = JSH
      JNW = ISH
      KNW = KSH
      LNW = LSH

      LSTRI = NGTH(2)
      LSTRJ = NGTH(1)
      LSTRK = NGTH(3)
      LSTRL = NGTH(4)

      IEQJ = .false.
      KEQL = .false.

      CALL SETUP_SHELLAI(INW,LANGI,IPRIM,EXI,CCI,XA,YA,ZA,MINI,MAXI)
      CALL SETUP_SHELLEFP(JNW,LANGJ,JPRIM,EXJ,CCJ,XB,YB,ZB,MINJ,MAXJ,
     & JM,JNAT)

      IF(MODE.EQ.2) THEN
      CALL SETUP_SHELLEFP(KNW,LANGK,KPRIM,EXK,CCK,XC,YC,ZC,MINK,MAXK,
     & JM,JNAT)
      ELSE
      CALL SETUP_SHELLAI(KNW,LANGK,KPRIM,EXK,CCK,XC,YC,ZC,MINK,MAXK)
      KEQL = KNW.EQ.LNW
      ENDIF

      CALL SETUP_SHELLAI(LNW,LANGL,LPRIM,EXL,CCL,XD,YD,ZD,MINL,MAXL)

      RETURN
      END


C*MODULE EFTEI_ERIC  *DECK ERIC_QMEFP_CASE_3 
C>    @brief  ERIC scheme for case 3 (IJ|LK) 
C>    @author  Peng Xu and Tosaporn Sattasathuchana 
C>             Jan 2020
C>    @detail Case 2: KL switch
C>    @see    parameters same as ERIC_QMEFP_CASE_1 

      SUBROUTINE ERIC_QMEFP_CASE_3(MODE,JM,JNAT,ISH,JSH,KSH,LSH,
     &                             LANGI,LANGJ,LANGK,LANGL,
     &                             IPRIM,JPRIM,KPRIM,LPRIM,
     &                             IEQJ,KEQL,XA,YA,ZA,XB,YB,ZB,
     &                             XC,YC,ZC,XD,YD,ZD,
     &                             MINI,MAXI,MINJ,MAXJ,MINK,MAXK,
     &                             MINL,MAXL)

      use mx_limits, only: mxatm,mxgtot,mxsh,mxgsh,mxg2,MXGEFP
      USE comm_EFPBAS

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      COMMON /ERIPRM/ EXI(MXGSH),EXJ(MXGSH),EXK(MXGSH),EXL(MXGSH),
     *                CCI(MXGSH),CCJ(MXGSH),CCK(MXGSH),CCL(MXGSH),
     *                XAB(MXG2),YAB(MXG2),ZAB(MXG2),
     *                XCD(MXG2),YCD(MXG2),ZCD(MXG2),
     *                CCBRA(MXG2),CCKET(MXG2),RXB(MXG2),
     *                SLBRA(MXG2),SLKET(MXG2),RXK(MXG2)
C
      COMMON /ERIDAT/ LEN1,LEN2,LEN3,LEN4
      COMMON /ERIOUT/ INW,JNW,KNW,LNW,LSTRI,LSTRJ,LSTRK,LSTRL
C
      COMMON /SHLEXC/ NORGSH(3),NORGSP(3),IEXCH,NANGM,NGTH(4)
      LOGICAL IEQJ,KEQL

      INW = ISH
      JNW = JSH
      KNW = LSH
      LNW = KSH

      LSTRI = NGTH(1)
      LSTRJ = NGTH(2)
      LSTRK = NGTH(4)
      LSTRL = NGTH(3)

      IEQJ = .false.
      KEQL = .false.

      CALL SETUP_SHELLEFP(INW,LANGI,IPRIM,EXI,CCI,XA,YA,ZA,MINI,MAXI,
     & JM,JNAT)

      CALL SETUP_SHELLAI(JNW,LANGJ,JPRIM,EXJ,CCJ,XB,YB,ZB,MINJ,MAXJ)

      CALL SETUP_SHELLAI(KNW,LANGK,KPRIM,EXK,CCK,XC,YC,ZC,MINK,MAXK)

      IF(MODE.EQ.2) THEN
      CALL SETUP_SHELLEFP(LNW,LANGL,LPRIM,EXL,CCL,XD,YD,ZD,MINL,MAXL,
     & JM,JNAT)
      ELSE
      CALL SETUP_SHELLAI(LNW,LANGL,LPRIM,EXL,CCL,XD,YD,ZD,MINL,MAXL)
      KEQL = KNW.EQ.LNW
      ENDIF

      RETURN
      END

C*MODULE EFTEI_ERIC  *DECK ERIC_QMEFP_CASE_4 
C>    @brief  ERIC scheme for case 4 (JI|LK) 
C>    @author  Peng Xu and Tosaporn Sattasathuchana 
C>             Jan 2020
C>    @detail Case 4: IJ switch and KL switch
C>    @see    parameters same as ERIC_QMEFP_CASE_1 

      SUBROUTINE ERIC_QMEFP_CASE_4(MODE,JM,JNAT,ISH,JSH,KSH,LSH,
     &                             LANGI,LANGJ,LANGK,LANGL,
     &                             IPRIM,JPRIM,KPRIM,LPRIM,
     &                             IEQJ,KEQL,XA,YA,ZA,XB,YB,ZB,
     &                             XC,YC,ZC,XD,YD,ZD,
     &                             MINI,MAXI,MINJ,MAXJ,MINK,MAXK,
     &                             MINL,MAXL)

      use mx_limits, only: mxatm,mxgtot,mxsh,mxgsh,mxg2,MXGEFP
      USE comm_EFPBAS

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      COMMON /ERIPRM/ EXI(MXGSH),EXJ(MXGSH),EXK(MXGSH),EXL(MXGSH),
     *                CCI(MXGSH),CCJ(MXGSH),CCK(MXGSH),CCL(MXGSH),
     *                XAB(MXG2),YAB(MXG2),ZAB(MXG2),
     *                XCD(MXG2),YCD(MXG2),ZCD(MXG2),
     *                CCBRA(MXG2),CCKET(MXG2),RXB(MXG2),
     *                SLBRA(MXG2),SLKET(MXG2),RXK(MXG2)
C
      COMMON /ERIDAT/ LEN1,LEN2,LEN3,LEN4
      COMMON /ERIOUT/ INW,JNW,KNW,LNW,LSTRI,LSTRJ,LSTRK,LSTRL
C
      COMMON /SHLEXC/ NORGSH(3),NORGSP(3),IEXCH,NANGM,NGTH(4)
      LOGICAL IEQJ,KEQL

      INW = JSH
      JNW = ISH
      KNW = LSH
      LNW = KSH

      LSTRI = NGTH(2)
      LSTRJ = NGTH(1)
      LSTRK = NGTH(4)
      LSTRL = NGTH(3)

      IEQJ = .false.
      KEQL = .false.

      CALL SETUP_SHELLAI(INW,LANGI,IPRIM,EXI,CCI,XA,YA,ZA,MINI,MAXI)

      CALL SETUP_SHELLEFP(JNW,LANGJ,JPRIM,EXJ,CCJ,XB,YB,ZB,MINJ,MAXJ,
     & JM,JNAT)

      CALL SETUP_SHELLAI(KNW,LANGK,KPRIM,EXK,CCK,XC,YC,ZC,MINK,MAXK)

      IF(MODE.EQ.2) THEN
      CALL SETUP_SHELLEFP(LNW,LANGL,LPRIM,EXL,CCL,XD,YD,ZD,MINL,MAXL,
     & JM,JNAT)
      ELSE
      CALL SETUP_SHELLAI(LNW,LANGL,LPRIM,EXL,CCL,XD,YD,ZD,MINL,MAXL)
      KEQL = KNW.EQ.LNW
      ENDIF

      RETURN
      END

C*MODULE EFTEI_ERIC  *DECK ERIC_QMEFP_CASE_5 
C>    @brief  ERIC scheme for case 5 (KL|IJ) 
C>    @author  Peng Xu and Tosaporn Sattasathuchana 
C>             Jan 2020
C>    @detail Case 5: pairs IJ and KL switch 
C>    @see    parameters same as ERIC_QMEFP_CASE_1 

      SUBROUTINE ERIC_QMEFP_CASE_5(MODE,JM,JNAT,ISH,JSH,KSH,LSH,
     &                             LANGI,LANGJ,LANGK,LANGL,
     &                             IPRIM,JPRIM,KPRIM,LPRIM,
     &                             IEQJ,KEQL,XA,YA,ZA,XB,YB,ZB,
     &                             XC,YC,ZC,XD,YD,ZD,
     &                             MINI,MAXI,MINJ,MAXJ,MINK,MAXK,
     &                             MINL,MAXL)

      use mx_limits, only: mxatm,mxgtot,mxsh,mxgsh,mxg2,MXGEFP
      USE comm_EFPBAS

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      COMMON /ERIPRM/ EXI(MXGSH),EXJ(MXGSH),EXK(MXGSH),EXL(MXGSH),
     *                CCI(MXGSH),CCJ(MXGSH),CCK(MXGSH),CCL(MXGSH),
     *                XAB(MXG2),YAB(MXG2),ZAB(MXG2),
     *                XCD(MXG2),YCD(MXG2),ZCD(MXG2),
     *                CCBRA(MXG2),CCKET(MXG2),RXB(MXG2),
     *                SLBRA(MXG2),SLKET(MXG2),RXK(MXG2)
C
      COMMON /ERIDAT/ LEN1,LEN2,LEN3,LEN4
      COMMON /ERIOUT/ INW,JNW,KNW,LNW,LSTRI,LSTRJ,LSTRK,LSTRL
C
      COMMON /SHLEXC/ NORGSH(3),NORGSP(3),IEXCH,NANGM,NGTH(4)
      LOGICAL IEQJ,KEQL

      INW = KSH
      JNW = LSH
      KNW = ISH
      LNW = JSH

      LSTRI = NGTH(3)
      LSTRJ = NGTH(4)
      LSTRK = NGTH(1)
      LSTRL = NGTH(2)

      IEQJ = .false.
      KEQL = .false.

      IF(MODE.EQ.2) THEN
      CALL SETUP_SHELLEFP(INW,LANGI,IPRIM,EXI,CCI,XA,YA,ZA,MINI,MAXI,
     & JM,JNAT)
      ELSE
      CALL SETUP_SHELLAI(INW,LANGI,IPRIM,EXI,CCI,XA,YA,ZA,MINI,MAXI)
      IEQJ = INW.EQ.JNW
      ENDIF

      CALL SETUP_SHELLAI(JNW,LANGJ,JPRIM,EXJ,CCJ,XB,YB,ZB,MINJ,MAXJ)

      CALL SETUP_SHELLEFP(KNW,LANGK,KPRIM,EXK,CCK,XC,YC,ZC,MINK,MAXK,
     & JM,JNAT)

      CALL SETUP_SHELLAI(LNW,LANGL,LPRIM,EXL,CCL,XD,YD,ZD,MINL,MAXL)

      RETURN
      END


C*MODULE EFTEI_ERIC  *DECK ERIC_QMEFP_CASE_6 
C>    @brief  ERIC scheme for case 6 (KL|JI) 
C>    @author  Peng Xu and Tosaporn Sattasathuchana 
C>             Jan 2020
C>    @detail Case 6: pairs IJ and KL switch followed by KL switch
C>    @see    parameters same as ERIC_QMEFP_CASE_1 

      SUBROUTINE ERIC_QMEFP_CASE_6(MODE,JM,JNAT,ISH,JSH,KSH,LSH,
     &                             LANGI,LANGJ,LANGK,LANGL,
     &                             IPRIM,JPRIM,KPRIM,LPRIM,
     &                             IEQJ,KEQL,XA,YA,ZA,XB,YB,ZB,
     &                             XC,YC,ZC,XD,YD,ZD,
     &                             MINI,MAXI,MINJ,MAXJ,MINK,MAXK,
     &                             MINL,MAXL)

      use mx_limits, only: mxatm,mxgtot,mxsh,mxgsh,mxg2,MXGEFP
      USE comm_EFPBAS

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      COMMON /ERIPRM/ EXI(MXGSH),EXJ(MXGSH),EXK(MXGSH),EXL(MXGSH),
     *                CCI(MXGSH),CCJ(MXGSH),CCK(MXGSH),CCL(MXGSH),
     *                XAB(MXG2),YAB(MXG2),ZAB(MXG2),
     *                XCD(MXG2),YCD(MXG2),ZCD(MXG2),
     *                CCBRA(MXG2),CCKET(MXG2),RXB(MXG2),
     *                SLBRA(MXG2),SLKET(MXG2),RXK(MXG2)
C
      COMMON /ERIDAT/ LEN1,LEN2,LEN3,LEN4
      COMMON /ERIOUT/ INW,JNW,KNW,LNW,LSTRI,LSTRJ,LSTRK,LSTRL
C
      COMMON /SHLEXC/ NORGSH(3),NORGSP(3),IEXCH,NANGM,NGTH(4)
      LOGICAL IEQJ,KEQL

      INW = KSH
      JNW = LSH
      KNW = JSH
      LNW = ISH

      LSTRI = NGTH(3)
      LSTRJ = NGTH(4)
      LSTRK = NGTH(2)
      LSTRL = NGTH(1)

      IEQJ = .false.
      KEQL = .false.

      IF(MODE.EQ.2) THEN
      CALL SETUP_SHELLEFP(INW,LANGI,IPRIM,EXI,CCI,XA,YA,ZA,MINI,MAXI,
     & JM,JNAT)
      ELSE
      CALL SETUP_SHELLAI(INW,LANGI,IPRIM,EXI,CCI,XA,YA,ZA,MINI,MAXI)
      IEQJ = INW.EQ.JNW
      ENDIF

      CALL SETUP_SHELLAI(JNW,LANGJ,JPRIM,EXJ,CCJ,XB,YB,ZB,MINJ,MAXJ)

      CALL SETUP_SHELLAI(KNW,LANGK,KPRIM,EXK,CCK,XC,YC,ZC,MINK,MAXK)

      CALL SETUP_SHELLEFP(LNW,LANGL,LPRIM,EXL,CCL,XD,YD,ZD,MINL,MAXL,
     & JM,JNAT)

      RETURN
      END


C*MODULE EFTEI_ERIC  *DECK ERIC_QMEFP_CASE_7 
C>    @brief  ERIC scheme for case 7 (LK|IJ) 
C>    @author  Peng Xu and Tosaporn Sattasathuchana 
C>             Jan 2020
C>    @detail Case 7: pairs IJ and KL switch, followed by IJ switch
C>    @see    parameters same as ERIC_QMEFP_CASE_1 

      SUBROUTINE ERIC_QMEFP_CASE_7(MODE,JM,JNAT,ISH,JSH,KSH,LSH,
     &                             LANGI,LANGJ,LANGK,LANGL,
     &                             IPRIM,JPRIM,KPRIM,LPRIM,
     &                             IEQJ,KEQL,XA,YA,ZA,XB,YB,ZB,
     &                             XC,YC,ZC,XD,YD,ZD,
     &                             MINI,MAXI,MINJ,MAXJ,MINK,MAXK,
     &                             MINL,MAXL)

      use mx_limits, only: mxatm,mxgtot,mxsh,mxgsh,mxg2,MXGEFP
      USE comm_EFPBAS

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      COMMON /ERIPRM/ EXI(MXGSH),EXJ(MXGSH),EXK(MXGSH),EXL(MXGSH),
     *                CCI(MXGSH),CCJ(MXGSH),CCK(MXGSH),CCL(MXGSH),
     *                XAB(MXG2),YAB(MXG2),ZAB(MXG2),
     *                XCD(MXG2),YCD(MXG2),ZCD(MXG2),
     *                CCBRA(MXG2),CCKET(MXG2),RXB(MXG2),
     *                SLBRA(MXG2),SLKET(MXG2),RXK(MXG2)
C
      COMMON /ERIDAT/ LEN1,LEN2,LEN3,LEN4
      COMMON /ERIOUT/ INW,JNW,KNW,LNW,LSTRI,LSTRJ,LSTRK,LSTRL
C
      COMMON /SHLEXC/ NORGSH(3),NORGSP(3),IEXCH,NANGM,NGTH(4)
      LOGICAL IEQJ,KEQL

      INW = LSH
      JNW = KSH
      KNW = ISH
      LNW = JSH

      LSTRI = NGTH(4)
      LSTRJ = NGTH(3)
      LSTRK = NGTH(1)
      LSTRL = NGTH(2)

      IEQJ = .false.
      KEQL = .false.

      CALL SETUP_SHELLAI(INW,LANGI,IPRIM,EXI,CCI,XA,YA,ZA,MINI,MAXI)

      IF(MODE.EQ.2) THEN
      CALL SETUP_SHELLEFP(JNW,LANGJ,JPRIM,EXJ,CCJ,XB,YB,ZB,MINJ,MAXJ,
     &                    JM,JNAT)
      ELSE
      CALL SETUP_SHELLAI(JNW,LANGJ,JPRIM,EXJ,CCJ,XB,YB,ZB,MINJ,MAXJ)
      ENDIF

      CALL SETUP_SHELLEFP(KNW,LANGK,KPRIM,EXK,CCK,XC,YC,ZC,MINK,MAXK,
     &                    JM,JNAT)

      CALL SETUP_SHELLAI(LNW,LANGL,LPRIM,EXL,CCL,XD,YD,ZD,MINL,MAXL)

      RETURN
      END

C*MODULE EFTEI_ERIC  *DECK ERIC_QMEFP_CASE_8 
C>    @brief  ERIC scheme for case 8 (LK|JI) 
C>    @author  Peng Xu and Tosaporn Sattasathuchana 
C>             Jan 2020
C>    @detail Case 8: pairs IJ and KL switch followed by IJ switch
C>                    and KL switch
C>    @see    parameters same as ERIC_QMEFP_CASE_1 

      SUBROUTINE ERIC_QMEFP_CASE_8(MODE,JM,JNAT,ISH,JSH,KSH,LSH,
     &                             LANGI,LANGJ,LANGK,LANGL,
     &                             IPRIM,JPRIM,KPRIM,LPRIM,
     &                             IEQJ,KEQL,XA,YA,ZA,XB,YB,ZB,
     &                             XC,YC,ZC,XD,YD,ZD,
     &                             MINI,MAXI,MINJ,MAXJ,MINK,MAXK,
     &                             MINL,MAXL)

      use mx_limits, only: mxatm,mxgtot,mxsh,mxgsh,mxg2,MXGEFP
      USE comm_EFPBAS

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      COMMON /ERIPRM/ EXI(MXGSH),EXJ(MXGSH),EXK(MXGSH),EXL(MXGSH),
     *                CCI(MXGSH),CCJ(MXGSH),CCK(MXGSH),CCL(MXGSH),
     *                XAB(MXG2),YAB(MXG2),ZAB(MXG2),
     *                XCD(MXG2),YCD(MXG2),ZCD(MXG2),
     *                CCBRA(MXG2),CCKET(MXG2),RXB(MXG2),
     *                SLBRA(MXG2),SLKET(MXG2),RXK(MXG2)
C
      COMMON /ERIDAT/ LEN1,LEN2,LEN3,LEN4
      COMMON /ERIOUT/ INW,JNW,KNW,LNW,LSTRI,LSTRJ,LSTRK,LSTRL
C
      COMMON /SHLEXC/ NORGSH(3),NORGSP(3),IEXCH,NANGM,NGTH(4)
      LOGICAL IEQJ,KEQL
 
      INW = LSH
      JNW = KSH
      KNW = JSH
      LNW = ISH

      LSTRI = NGTH(4)
      LSTRJ = NGTH(3)
      LSTRK = NGTH(2)
      LSTRL = NGTH(1)

      IEQJ = .false.
      KEQL = .false.

      CALL SETUP_SHELLAI(INW,LANGI,IPRIM,EXI,CCI,XA,YA,ZA,MINI,MAXI)

      IF(MODE.EQ.2) THEN
      CALL SETUP_SHELLEFP(JNW,LANGJ,JPRIM,EXJ,CCJ,XB,YB,ZB,MINJ,MAXJ,
     &                    JM,JNAT)
      ELSE
      CALL SETUP_SHELLAI(JNW,LANGJ,JPRIM,EXJ,CCJ,XB,YB,ZB,MINJ,MAXJ)
      IEQJ = INW.EQ.JNW
      ENDIF

      CALL SETUP_SHELLAI(KNW,LANGK,KPRIM,EXK,CCK,XC,YC,ZC,MINK,MAXK)

      CALL SETUP_SHELLEFP(LNW,LANGL,LPRIM,EXL,CCL,XD,YD,ZD,MINL,MAXL,
     &                    JM,JNAT)

      RETURN
      END

 
C*MODULE EFTEI_ERIC   *DECK ERIC_QMEFP_TS
C>
C>    @brief   ERIC: Electron Repulsion Integral Calculator
C>
C>    @details a thread-safe version of ERIC_QMEFP subroutine
C>
C>    @author  Peng Xu and Tosaporn Sattasathuchana 
C>             Jan 2020
C>
C>    @param   JM: index for the type of EFP potentials 
C>    @param   JNAT: counter for number of EFP atoms 
C>    @param   MODE: different type of TEI involve in QM/EFP EXREP
C>                 MODE=2 two-center on QM, two-center on EFP
C>                 MODE=3 three-center on QM, One-center on EFP
C>
      SUBROUTINE ERIC_QMEFP_TS(ISH,JSH,KSH,LSH,ERI,JM,JNAT,MODE)
      use mx_limits, only: mxatm,mxgtot,mxsh,mxgsh,mxg2,MXGEFP
      USE comm_EFPBAS
C ----------------------------------------------------------------------
C
C              ELECTRON REPULSION INTEGRAL CALCULATOR
C
C ----------------------------------------------------------------------
C
C  METHOD:
C     A) FORM SCALED, CONTRACTED 1-CENTER PRECURSOR INTEGRALS.
C        CONVERT THESE TO 4-CENTER INTEGRALS OVER CARTESIAN
C        GAUSSIANS USING,
C     B) PRECURSOR-HERMITE TRANSFER EQUATION (PTE)
C     C) CONTRACTED TRANSFER  EQUATION (CTE)
C     D) HORIZONTAL RECURSION RELATION (HRR)
C
C
C  REFERENCES
C
C     PTE:
C     "Recursion Formula for Electron Repulsion Integrals Over
C     Hermite Polynomials"
C         G.D.Fletcher  Int.J.Quantum Chem. 106, 355-360(2006)
C
C     HRR:
C     M. HEAD-GORDON & J. A. POPLE,
C     J. CHEM. PHYS., 89, 5777-5786 (1988).
C
C     CTE:
C     P. M. W. GILL, M. HEAD-GORDON, & J. A. POPLE,
C     INT. J. Q. CHEM., SYMP. 23, 269-280 (1989).
C
C     INTERPOLATION METHOD:
C     P. M. W. GILL, B. G. JOHNSON, & J. A. POPLE,
C     INT. J. Q. CHEM., 40, 745-752 (1991).
C
C
C  THIS IS THE DRIVER ROUTINE THAT INTERFACES TO GAMESS COMMON NSHEL
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      DOUBLE PRECISION ERI(*)
C
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /NSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),
     *                CF(MXGTOT),CG(MXGTOT),CH(MXGTOT),CI(MXGTOT),
     *                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),
     *                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
C
C  /ERIPRM/ SYMBOLS:
C     EXI,J,K,L   = PRIMITIVE EXPONENTS OF I,J,K,L SHELLS
C     CCI,J,K,L   = CONTRACTION COEFFICIENTS OF ISH,JSH,KSH,LSH
C     XAB,YAB,ZAB = EXPONENT-WEIGHTED DISTANCE FOR BRA
C     XCD,YCD,ZCD = EXPONENT-WEIGHTED DISTANCE FOR KET
C     CCBRA,CCKET = PRIMITIVE CHARGE-CLOUD FACTOR FOR BRA,KET
C     SLBRA,SLKET = FACTOR CONVERTING SMALL-T EXPRESSION TO LARGE-T
C
      COMMON /ERIPRM/ EXI(MXGSH),EXJ(MXGSH),EXK(MXGSH),EXL(MXGSH),
     *                CCI(MXGSH),CCJ(MXGSH),CCK(MXGSH),CCL(MXGSH),
     *                XAB(MXG2),YAB(MXG2),ZAB(MXG2),
     *                XCD(MXG2),YCD(MXG2),ZCD(MXG2),
     *                CCBRA(MXG2),CCKET(MXG2),RXB(MXG2),
     *                SLBRA(MXG2),SLKET(MXG2),RXK(MXG2)
C
      COMMON /ERIDAT/ LEN1,LEN2,LEN3,LEN4
      COMMON /ERIOUT/ INW,JNW,KNW,LNW,LSTRI,LSTRJ,LSTRK,LSTRL
C
      COMMON /SHLEXC/ NORGSH(3),NORGSP(3),IEXCH,NANGM,NGTH(4)
      COMMON /SHLNOS/ QQ4,LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,
     2                MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL,
     3                NIJ,IJ,KL,IJKL

      COMMON /FMCOM /XX(1)
C
      LOGICAL IEQJ,KEQL
      LOGICAL TEST1,TEST2,TEST3
C
      PARAMETER (PT5=0.5D+00)
      PARAMETER (ONE=1.0D+00)
      PARAMETER (TWO=2.0D+00)
C
C     SR3 <==> SQRT( 3)   <==> 1.7320508075688773D+00
C     SR5 <==> SQRT( 5)   <==> 2.2360679774997897D+00
C     SR7 <==> SQRT( 7)   <==> 2.6457513110645906D+00
C     S15 <==> SQRT(15)   <==> 3.8729833462074169D+00
C     S35 <==> SQRT(35)   <==> 5.9160797830996160D+00
C     S53 <==> SQRT(35/3) <==> 3.4156502553198661D+00
C
      PARAMETER (SR3=1.7320508075688773D+00)
      PARAMETER (SR5=2.2360679774997897D+00)
      PARAMETER (SR7=2.6457513110645906D+00)
      PARAMETER (S15=SR3*SR5)
      PARAMETER (S35=SR5*SR7)
      PARAMETER (S53=S35/SR3)
C
      PARAMETER (PI254 =5.9149671727956129D+00)
      PARAMETER (PI214 =0.94139626377671481D+00)
 
      DIMENSION xxiphi(1182), xxiwk1(25), xxiwk2(496)
      SAVE xxiphi, xxiwk1, xxiwk2
!$omp threadprivate(xxiphi, xxiwk1, xxiwk2)     
C
      DIMENSION IORD(35),ANGL(35)
      DATA IORD/
     1       1,
     2       2,  3,  4,
     3       5,  7, 10,  6,  8,  9,
     4      11, 14, 20, 12, 15, 13, 17, 18, 19, 16,
     5      21, 25, 35, 22, 26, 24, 29, 33, 34, 23, 30, 32, 27, 28, 31/
      DATA ANGL/
     1     ONE,
     2     ONE,ONE,ONE,
     3     ONE,SR3,ONE,SR3,SR3,ONE,
     4     ONE,SR5,SR5,ONE,SR5,S15,SR5,SR5,SR5,ONE,
     5     ONE,SR7,S53,SR7,ONE,SR7,S35,S35,SR7,S53,S35,S53,SR7,SR7,ONE/
C
C  FAST CODES
C
C  RE-ORDER CHARGE-CLOUD SHELLS
C  (HRR DATA ONLY STORED FOR LI>=LJ)
C
      INW = ISH
      JNW = JSH
      KNW = KSH
      LNW = LSH
C
      LSTRI = NGTH(1)
      LSTRJ = NGTH(2)
      LSTRK = NGTH(3)
      LSTRL = NGTH(4)
C
      LANGI = KTYPEF(INW,JM) - 1
      LANGJ = KTYPE(JNW) - 1
      IF(MODE.EQ.2) THEN
      LANGK = KTYPEF(KNW,JM) - 1
      ELSE
      LANGK = KTYPE(KNW) - 1
      ENDIF
      LANGL = KTYPE(LNW) - 1

      LBRA  = LANGI + LANGJ
      LKET  = LANGK + LANGL

      TEST1 = LBRA.GT.LKET
      TEST2 = LBRA.EQ.LKET

      IF((LANGI.LT.LANGJ).AND.(LANGK.LT.LANGL)) THEN

        INW = JSH
        JNW = ISH
        KNW = LSH
        LNW = KSH

        LANGI = KTYPE(INW)
        LANGJ = KTYPEF(JNW,JM)
        LANGK = KTYPE(KNW)
        IF(MODE.EQ.2) THEN
        LANGL = KTYPEF(LNW,JM) 
        ELSE
        LANGL = KTYPE(LNW)
        ENDIF

        LNGIJ = (LANGI*LANGI-LANGI)/2 + LANGJ
        LNGKL = (LANGK*LANGK-LANGK)/2 + LANGL

        TEST3 = LNGKL.GT.LNGIJ

        IF(TEST1) THEN
          CALL ERIC_QMEFP_CASE_8(MODE,JM,JNAT,ISH,JSH,KSH,LSH,
     &         LANGI,LANGJ,LANGK,LANGL,IPRIM,JPRIM,KPRIM,LPRIM,
     &         IEQJ,KEQL,XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD,
     &         MINI,MAXI,MINJ,MAXJ,MINK,MAXK,MINL,MAXL)
        ELSEIF(TEST2.and.TEST3) THEN
          CALL ERIC_QMEFP_CASE_8(MODE,JM,JNAT,ISH,JSH,KSH,LSH,
     &         LANGI,LANGJ,LANGK,LANGL,IPRIM,JPRIM,KPRIM,LPRIM,
     &         IEQJ,KEQL,XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD,
     &         MINI,MAXI,MINJ,MAXJ,MINK,MAXK,MINL,MAXL)
        ELSE
          CALL ERIC_QMEFP_CASE_4(MODE,JM,JNAT,ISH,JSH,KSH,LSH,
     &         LANGI,LANGJ,LANGK,LANGL,IPRIM,JPRIM,KPRIM,LPRIM,
     &         IEQJ,KEQL,XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD,
     &         MINI,MAXI,MINJ,MAXJ,MINK,MAXK,MINL,MAXL)
        ENDIF          

      ELSEIF((LANGI.LT.LANGJ).AND.(LANGK.GE.LANGL)) THEN

        INW = JSH
        JNW = ISH
        KNW = KSH
        LNW = LSH

        LANGI = KTYPE(INW)
        LANGJ = KTYPEF(JNW,JM)
        IF(MODE.EQ.2) THEN
        LANGK = KTYPEF(KNW,JM)
        ELSE
        LANGK = KTYPE(KNW)
        ENDIF
        LANGL = KTYPE(LNW)

        LNGIJ = (LANGI*LANGI-LANGI)/2 + LANGJ
        LNGKL = (LANGK*LANGK-LANGK)/2 + LANGL

        TEST3 = LNGKL.GT.LNGIJ

        IF(TEST1) THEN
          CALL ERIC_QMEFP_CASE_6(MODE,JM,JNAT,ISH,JSH,KSH,LSH,
     &         LANGI,LANGJ,LANGK,LANGL,IPRIM,JPRIM,KPRIM,LPRIM,
     &         IEQJ,KEQL,XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD,
     &         MINI,MAXI,MINJ,MAXJ,MINK,MAXK,MINL,MAXL)
        ELSEIF(TEST2.and.TEST3) THEN
          CALL ERIC_QMEFP_CASE_6(MODE,JM,JNAT,ISH,JSH,KSH,LSH,
     &         LANGI,LANGJ,LANGK,LANGL,IPRIM,JPRIM,KPRIM,LPRIM,
     &         IEQJ,KEQL,XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD,
     &         MINI,MAXI,MINJ,MAXJ,MINK,MAXK,MINL,MAXL)
        ELSE
          CALL ERIC_QMEFP_CASE_2(MODE,JM,JNAT,ISH,JSH,KSH,LSH,
     &         LANGI,LANGJ,LANGK,LANGL,IPRIM,JPRIM,KPRIM,LPRIM,
     &         IEQJ,KEQL,XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD,
     &         MINI,MAXI,MINJ,MAXJ,MINK,MAXK,MINL,MAXL)
        ENDIF

      ELSEIF((LANGI.GE.LANGJ).AND.(LANGK.LT.LANGL)) THEN

        INW = ISH
        JNW = JSH
        KNW = LSH
        LNW = KSH
        
        LANGI = KTYPEF(INW,JM)
        LANGJ = KTYPE(JNW)
        LANGK = KTYPE(KNW)
        IF(MODE.EQ.2) THEN
        LANGL = KTYPEF(LNW,JM)
        ELSE
        LANGL = KTYPE(LNW)
        ENDIF

        LNGIJ = (LANGI*LANGI-LANGI)/2 + LANGJ
        LNGKL = (LANGK*LANGK-LANGK)/2 + LANGL

        TEST3 = LNGKL.GT.LNGIJ

        IF(TEST1) THEN
          CALL ERIC_QMEFP_CASE_7(MODE,JM,JNAT,ISH,JSH,KSH,LSH,
     &         LANGI,LANGJ,LANGK,LANGL,IPRIM,JPRIM,KPRIM,LPRIM,
     &         IEQJ,KEQL,XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD,
     &         MINI,MAXI,MINJ,MAXJ,MINK,MAXK,MINL,MAXL)
        ELSEIF(TEST2.and.TEST3) THEN
          CALL ERIC_QMEFP_CASE_7(MODE,JM,JNAT,ISH,JSH,KSH,LSH,
     &         LANGI,LANGJ,LANGK,LANGL,IPRIM,JPRIM,KPRIM,LPRIM,
     &         IEQJ,KEQL,XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD,
     &         MINI,MAXI,MINJ,MAXJ,MINK,MAXK,MINL,MAXL)
        ELSE
          CALL ERIC_QMEFP_CASE_3(MODE,JM,JNAT,ISH,JSH,KSH,LSH,
     &         LANGI,LANGJ,LANGK,LANGL,IPRIM,JPRIM,KPRIM,LPRIM,
     &         IEQJ,KEQL,XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD,
     &         MINI,MAXI,MINJ,MAXJ,MINK,MAXK,MINL,MAXL)
        ENDIF

      ELSEIF((LANGI.GE.LANGJ).AND.(LANGK.GE.LANGL)) THEN
       
        LANGI = KTYPEF(INW,JM)
        LANGJ = KTYPE(JNW)
        IF(MODE.EQ.2) THEN
        LANGK = KTYPEF(KNW,JM)
        ELSE
        LANGK = KTYPE(KNW)
        ENDIF
        LANGL = KTYPE(LNW)

        LNGIJ = (LANGI*LANGI-LANGI)/2 + LANGJ
        LNGKL = (LANGK*LANGK-LANGK)/2 + LANGL

        TEST3 = LNGKL.GT.LNGIJ

        IF(TEST1) THEN
          CALL ERIC_QMEFP_CASE_5(MODE,JM,JNAT,ISH,JSH,KSH,LSH,
     &         LANGI,LANGJ,LANGK,LANGL,IPRIM,JPRIM,KPRIM,LPRIM,
     &         IEQJ,KEQL,XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD,
     &         MINI,MAXI,MINJ,MAXJ,MINK,MAXK,MINL,MAXL)
        ELSEIF(TEST2.and.TEST3) THEN
          CALL ERIC_QMEFP_CASE_5(MODE,JM,JNAT,ISH,JSH,KSH,LSH,
     &         LANGI,LANGJ,LANGK,LANGL,IPRIM,JPRIM,KPRIM,LPRIM,
     &         IEQJ,KEQL,XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD,
     &         MINI,MAXI,MINJ,MAXJ,MINK,MAXK,MINL,MAXL)
        ELSE
          CALL ERIC_QMEFP_CASE_1(MODE,JM,JNAT,ISH,JSH,KSH,LSH,
     &         LANGI,LANGJ,LANGK,LANGL,IPRIM,JPRIM,KPRIM,LPRIM,
     &         IEQJ,KEQL,XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD,
     &         MINI,MAXI,MINJ,MAXJ,MINK,MAXK,MINL,MAXL)
        ENDIF

      ELSE
        write(6,*) 'SHOULD NOT GET HERE!'
      ENDIF

      RAB = (XA-XB)**2 + (YA-YB)**2 + (ZA-ZB)**2
      IJ = 0
      DO JJ = 1, JPRIM
         EJ = EXJ(JJ)
         CJ = CCJ(JJ)
         ITOP = IPRIM
         IF(IEQJ) ITOP = JJ
         DO II = 1, ITOP
            EI = EXI(II)
            CIX = CCI(II)
            EIJ = ONE/(EI+EJ)
            IJ = IJ + 1
            CCFAC = PI254*CIX*CJ*EIJ*EXP( -EI*EJ*RAB*EIJ )
            IF(IEQJ .AND. II.NE.JJ) CCFAC = CCFAC*TWO
            CCBRA(IJ) = CCFAC
            SLBRA(IJ) = PI214*SQRT(EIJ)
            XAB(IJ) = (EI*XA + EJ*XB)*EIJ
            YAB(IJ) = (EI*YA + EJ*YB)*EIJ
            ZAB(IJ) = (EI*ZA + EJ*ZB)*EIJ
            RXB(IJ) = EIJ*PT5
         END DO
      END DO

      RCD = (XC-XD)**2 + (YC-YD)**2 + (ZC-ZD)**2
      KL = 0
      DO LL = 1, LPRIM
         EL = EXL(LL)
         CL = CCL(LL)
         KTOP = KPRIM
         IF(KEQL) KTOP = LL
         DO KK = 1, KTOP
            EK = EXK(KK)
            CK = CCK(KK)
            EKL = ONE/(EK+EL)
            KL = KL + 1
            CCFAC = PI254*CK*CL*EKL*EXP( -EK*EL*RCD*EKL )
            IF(KEQL .AND. KK.NE.LL) CCFAC = CCFAC*TWO
            CCKET(KL) = CCFAC
            SLKET(KL) = PI214*SQRT(EKL)
            XCD(KL) = (EK*XC + EL*XD)*EKL
            YCD(KL) = (EK*YC + EL*YD)*EKL
            ZCD(KL) = (EK*ZC + EL*ZD)*EKL
            RXK(KL) = EKL*PT5
         END DO
      END DO
C
      CALL VALFM(LOADFM)
      IPHI = LOADFM + 1
C
C  ANGULAR MOMENTUM 4-INDEX
C
      LBGT  = MAX(LANGI,LANGJ)
      LBLT  = MIN(LANGI,LANGJ)
      LNGIJ = (LBGT*LBGT+LBGT)/2 + LBLT
      LKGT  = MAX(LANGK,LANGL)
      LKLT  = MIN(LANGK,LANGL)
      LNGKL = (LKGT*LKGT+LKGT)/2 + LKLT
      LQGT  = MAX(LNGKL,LNGIJ)
      LQLT  = MIN(LNGKL,LNGIJ)
      LIJKL = (LQGT*LQGT+LQGT)/2 + LQLT
C
      IF(LIJKL.LE.5) THEN
C
C  SP CASES
C
         IF(LIJKL.EQ.0) THEN
C  SSSS
            IDIM =           1
            IOFF = 0
            CALL SSSS (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL,
     *                 XXIPHI,XXIWK2,IDIM)
         ELSE IF(LIJKL.EQ.1) THEN
C  PSSS
            IDIM =           4
            IOFF = 1
            CALL PSSS (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL,
     *                 XC,YC,ZC,XD,YD,ZD,
     *                 XXIPHI,XXIWK2,IDIM)
         ELSE IF(LIJKL.EQ.2) THEN
C  PSPS
            IDIM =           4
            IOFF = 1
            CALL PSPS (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL
     *,                XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD
     *,                XXIPHI,XXIWK1,XXIWK2,IDIM)
         ELSE IF(LIJKL.EQ.3) THEN
C  PPSS
            IDIM =          24
            IOFF = 15
            CALL PPSS (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL,
     *                 XC,YC,ZC,XD,YD,ZD,
     *                 XXIPHI,XXIWK2,IDIM)
         ELSE IF(LIJKL.EQ.4) THEN
C  PPPS
            IDIM =          24
            IOFF = 15
            CALL PPPS (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL
     *,                XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD
     *,                XXIPHI,XXIWK1,XXIWK2,IDIM)
         ELSE IF(LIJKL.EQ.5) THEN
C  PPPP
            IDIM =          24
            IOFF = 15
            CALL PPPP (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL
     *,                XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD
     *,                XXIPHI,XXIWK1,XXIWK2,IDIM)
         END IF
C
C  END OF SP CASES
C
      ELSE IF(LIJKL.GE.6 .AND. LIJKL.LE.20) THEN
C
C  D CASES
C
         IF(LIJKL.EQ.6) THEN
C  DSSS
            IDIM =          11
            IOFF = 5
            CALL DSSS (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL,
     *                 XC,YC,ZC,XD,YD,ZD,
     *                 XXIPHI,XXIWK2,IDIM)
         ELSE IF(LIJKL.EQ.7) THEN
C  DSPS
            IDIM =          11
            IOFF = 5
            CALL DSPS (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL
     *,                XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD
     *,                XXIPHI,XXIWK1,XXIWK2,IDIM)
         ELSE IF(LIJKL.EQ.8) THEN
C  DSPP
            IDIM =          24
            IOFF = 15
            CALL DSPP (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL
     *,                XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD
     *,                XXIPHI,XXIWK1,XXIWK2,IDIM)
         ELSE IF(LIJKL.EQ.9) THEN
C  DSDS
            IDIM =          11
            IOFF = 5
            CALL DSDS (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL
     *,                XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD
     *,                XXIPHI,XXIWK1,XXIWK2,IDIM)
         ELSE IF(LIJKL.EQ.10) THEN
C  DPSS
            IDIM =          53
            IOFF = 35
            CALL DPSS (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL,
     *                 XC,YC,ZC,XD,YD,ZD,
     *                 XXIPHI,XXIWK2,IDIM)
         ELSE IF(LIJKL.EQ.11) THEN
C  DPPS
            IDIM =          53
            IOFF = 35
            CALL DPPS (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL
     *,                XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD
     *,                XXIPHI,XXIWK1,XXIWK2,IDIM)
         ELSE IF(LIJKL.EQ.12) THEN
C  DPPP
            IDIM =          53
            IOFF = 35
            CALL DPPP (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL
     *,                XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD
     *,                XXIPHI,XXIWK1,XXIWK2,IDIM)
         ELSE IF(LIJKL.EQ.13) THEN
C  DPDS
            IDIM =          53
            IOFF = 35
            CALL DPDS (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL
     *,                XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD
     *,                XXIPHI,XXIWK1,XXIWK2,IDIM)
         ELSE IF(LIJKL.EQ.15) THEN
C  DDSS
            IDIM =         165
            IOFF = 129
            CALL DDSS (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL,
     *                 XC,YC,ZC,XD,YD,ZD,
     *                 XXIPHI,XXIWK2,IDIM)
         ELSE IF(LIJKL.EQ.16) THEN
C  DDPS
            IDIM =         165
            IOFF = 129
            CALL DDPS (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL
     *,                XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD
     *,                XXIPHI,XXIWK1,XXIWK2,IDIM)
         END IF
C
C  END OF D CASES
C
      ELSE IF(LIJKL.GE.21 .AND. LIJKL.LE.54) THEN
C
C  F CASES
C
         IF(LIJKL.EQ.21) THEN
C  FSSS
            IDIM =          24
            IOFF = 14
            CALL FSSS (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL,
     *                 XC,YC,ZC,XD,YD,ZD,
     *                 XXIPHI,XXIWK2,IDIM)
         ELSE IF(LIJKL.EQ.22) THEN
C  FSPS
            IDIM =          24
            IOFF = 14
            CALL FSPS (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL
     *,                XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD
     *,                XXIPHI,XXIWK1,XXIWK2,IDIM)
         ELSE IF(LIJKL.EQ.23) THEN
C  FSPP
            IDIM =          24
            IOFF = 14
            CALL FSPP (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL
     *,                XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD
     *,                XXIPHI,XXIWK1,XXIWK2,IDIM)
         ELSE IF(LIJKL.EQ.24) THEN
C  FSDS
            IDIM =          24
            IOFF = 14
            CALL FSDS (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL
     *,                XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD
     *,                XXIPHI,XXIWK1,XXIWK2,IDIM)
         ELSE IF(LIJKL.EQ.28) THEN
C  FPSS
            IDIM =         100
            IOFF = 70
            CALL FPSS (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL,
     *                 XC,YC,ZC,XD,YD,ZD,
     *                 XXIPHI,XXIWK2,IDIM)
         ELSE IF(LIJKL.EQ.29) THEN
C  FPPS
            IDIM =         100
            IOFF = 70
            CALL FPPS (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL
     *,                XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD
     *,                XXIPHI,XXIWK1,XXIWK2,IDIM)
         ELSE IF(LIJKL.EQ.36) THEN
C  FDSS
            IDIM =         285
            IOFF = 225
            CALL FDSS (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL,
     *                 XC,YC,ZC,XD,YD,ZD,
     *                 XXIPHI,XXIWK2,IDIM)
         END IF
C
C  END OF F CASES
C
      ELSE IF(LIJKL.GE.55 .AND. LIJKL.LE.119) THEN
C
C  G CASES
C
         IF(LIJKL.EQ.55) THEN
C  GSSS
            IDIM =          46
            IOFF = 31
            CALL GSSS (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL,
     *                 XC,YC,ZC,XD,YD,ZD,
     *                 XXIPHI,XXIWK2,IDIM)
         ELSE IF(LIJKL.EQ.56) THEN
C  GSPS
            IDIM =          46
            IOFF = 31
            CALL GSPS (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL
     *,                XA,YA,ZA,XB,YB,ZB,XC,YC,ZC,XD,YD,ZD
     *,                XXIPHI,XXIWK1,XXIWK2,IDIM)
         ELSE IF(LIJKL.EQ.66) THEN
C  GPSS
            IDIM =         171
            IOFF = 126
            CALL GPSS (IPRIM,JPRIM,KPRIM,LPRIM,IEQJ,KEQL,
     *                 XC,YC,ZC,XD,YD,ZD,
     *                 XXIPHI,XXIWK2,IDIM)
         END IF
C
C  END OF G CASES
C
      END IF    ! LIJKL
C
C  ANGULAR NORMALIZATION AND SAVE TO OUTPUT ARRAY
C  WITH HONDO INDEXING AND REORDERING FOR GAMESS
C
      LENI = MAXI - MINI + 1
      LENK = MAXK - MINK + 1

      II = 1
      DO I = MINI, MAXI
         IO = IORD(I)
         AI = ANGL(IO)
         IJ = II
         DO J = MINJ, MAXJ
            JO = IORD(J)
            AIJ = ANGL(JO)*AI
            JC = ((JO-MINJ)*LENI + IO-MINI)*IDIM + IOFF + 1 
            IJK = IJ
            DO K = MINK, MAXK
               KO = IORD(K)
               AIJK = ANGL(KO)*AIJ
               IJKL = IJK
               DO L = MINL, MAXL
                  LO = IORD(L)
                  AIJKL = ANGL(LO)*AIJK
                  IR = (LO-MINL)*LENK + KO-MINK
                  ERI(IJKL) = XXIWK2(JC+IR)*AIJKL
                  IJKL = IJKL + LSTRL
               END DO
               IJK = IJK + LSTRK
            END DO
            IJ = IJ + LSTRJ
         END DO
         II = II + LSTRI
      END DO
!      CALL RETFM(NEED)

        INW = ISH
        JNW = JSH
        KNW = KSH
        LNW = LSH
        LSTRI = NGTH(1)
        LSTRJ = NGTH(2)
        LSTRK = NGTH(3)
        LSTRL = NGTH(4)

      RETURN
      END
