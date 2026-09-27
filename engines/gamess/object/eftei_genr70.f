C*MODULE EFTEI_GENR70   *DECK GENR70_QMEFP
C>
C>    @brief   rotated axis integration involving s,p,L shells for QMEFP
C>
C>    @author  Peng Xu and Tosaporn Sattasathuchana 
C>             Jan 2020
C>    @details This subroutine prepares the necessary shell information
C>             for QM/EFP ERI (rotated axis) calculation and 
C>             is based on the subroutine GENR70 in int2b.src
C>
C>    @param   JM: index for the type of EFP potentials 
C>    @param   JNAT: counter for number of EFP atoms 
C>    @param   MODE: different type of TEI involve in QM/EFP EXREP
C>                 MODE=2 two-center on QM, two-center on EFP
C>                 MODE=3 three-center on QM, One-center on EFP

      SUBROUTINE GENR70_QMEFP(JM,JNAT,MODE)
      use mx_limits, only: mxgtot,mxsh,mxgsh,mxg2
      USE comm_EFPBAS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      COMMON /GEOMPQ/ R12,RAB,X34,X43,AQZ,QPR,QPS,
     2                TX12(MXG2),TX21(MXG2),TY01(MXG2),TY02(MXG2),
     3                D00P(MXG2),D01P(MXG2),D10P(MXG2),D11P(MXG2),
     4                NGANGB
      COMMON /GOUT  / GPOPLE(256*3),NORG
      COMMON /INTAC2/ EI1,EI2,CUX
      COMMON /MAXC  / CMAX(MXGTOT),CMAXA(MXGSH),CMAXB(MXGSH),
     2                CMAXC(MXGSH),CMAXD(MXGSH),ISMLP(MXG2),ISMLQ
      COMMON /NSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),
     2                CF(MXGTOT),CG(MXGTOT),CH(MXGTOT),CI(MXGTOT),
     3                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),
     4                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
      COMMON /POPOUT/ LPOPI,LPOPJ,LPOPK,LPOPL
      COMMON /SHLG70/ ISHELL,JSHELL,KSHELL,LSHELL,INEW,JNEW,KNEW,LNEW
      COMMON /SHLLFO/ NGA,LA,EXA(MXGSH),CSA(MXGSH),CPA(MXGSH),
     2                NGB,LB,EXB(MXGSH),CSB(MXGSH),CPB(MXGSH),
     3                NGC,LC,EXC(MXGSH),CSC(MXGSH),CPC(MXGSH),
     4                NGD,LD,EXD(MXGSH),CSD(MXGSH),CPD(MXGSH)
C
C     ==================================================================
C
CJMS  LABELLED COMMON JMSGYH DEFINED FOR COMPUTATIONAL EFFICIENCY.
CJMS  IT IS ONLY USED IN THIS MODULE INT2B AND IN MODULE INT2R.
C
      COMMON /JMSGYH/ SQ(0:1,0:1)
      COMMON /KI2 / ACY,ACY2,AQX,AQX2,AQXY,Y03,Y04
      COMMON /KI3 / R00(25),R01(120),R02(156),R03(80),R04(15)
C
      DIMENSION P12(3,3),P34(3,3),P(3,3),T(3)
C
      PARAMETER (ZER=0.0D+00)
      PARAMETER (PT5=0.5D+00)
      PARAMETER (ONE=1.0D+00)
      PARAMETER (PT7=0.7D+00)
      PARAMETER (PT9=0.9D+00)
C        ACYCUT CHANGED FROM A VALUE OF 1E-4 IN 2004,
C        FOR ACCURACY WHEN USING DIFFUSE EXPONENTS.
      PARAMETER (ACYCUT=1.0D-10)
      PARAMETER (TENM12=1.0D-12)
C
C                PITO52=(PI +PI )*PI * SQRT(PI )
      PARAMETER (PITO52=34.986836655249726D+00)
C
      PIF= PITO52
C
      LAT= KTYPEF(ISHELL,JM)-1
      LBT= KTYPE(JSHELL)-1
      IF(MODE.EQ.2) THEN 
      LCT= KTYPEF(KSHELL,JM)-1
      ELSE !MODE=3
      LCT= KTYPE(KSHELL)-1
      ENDIF
      LDT= KTYPE(LSHELL)-1
      ITYPE= 1+LDT+2*(LCT+2*(LBT+2*LAT))

      R12= ZER
      R34= ZER
C
      SELECT CASE(ITYPE)
C
C TYPES 0000,0001,0011,0101,0111,1111 ARE UNALTERED
C
         CASE(1,2,4,6,8,16)
         CALL GENR70_QMEFP_INTYPE_1(MODE,JM,JNAT,
     *                         LAT,LBT,LCT,LDT,P12,P34,R34)
C
C TYPES 1001,1011 HAVE IJ SWITCHED
C
         CASE(10,12)
         CALL GENR70_QMEFP_INTYPE_2(MODE,JM,JNAT,
     *                         LAT,LBT,LCT,LDT,P12,P34,R34)
C
C TYPES 0010,0110 HAVE KL SWITCHED
C
         CASE(3,7)
         CALL GENR70_QMEFP_INTYPE_3(MODE,JM,JNAT,
     *                         LAT,LBT,LCT,LDT,P12,P34,R34)
C
C TYPES 0100,1100,1101 HAVE PAIRS IJ AND KL SWITCHED
C
         CASE(5,13,14)
         CALL GENR70_QMEFP_INTYPE_4(MODE,JM,JNAT,
     *                         LAT,LBT,LCT,LDT,P12,P34,R34)
C
C TYPE 1010 HAS IJ SWITCHED AND KL SWITCHED
C
         CASE(11)
         CALL GENR70_QMEFP_INTYPE_5(MODE,JM,JNAT,
     *                         LAT,LBT,LCT,LDT,P12,P34,R34)
C
C TYPE 1000  HAS PAIRS IJ AND KL SWITCHED FOLLOWED BY KL SWITCH
C
         CASE(9)
         CALL GENR70_QMEFP_INTYPE_6(MODE,JM,JNAT,
     *                         LAT,LBT,LCT,LDT,P12,P34,R34)
C
C TYPE 1110 HAS PAIRS IJ AND KL SWITCHED FOLLOWED BY IJ SWITCH
C
         CASE(15)
         CALL GENR70_QMEFP_INTYPE_7(MODE,JM,JNAT,
     *                         LAT,LBT,LCT,LDT,P12,P34,R34)
       END SELECT
C
C ONLY 6 STANDARD TYPES REMAIN. 0000,0001,0011,0101,0111,1111
C SPECIFY THESE BY JTYPE
C
      IF(ITYPE.EQ. 1) THEN
         JTYPE= 1
      ELSEIF(ITYPE.EQ. 2 .OR. ITYPE.EQ. 3 .OR.
     2       ITYPE.EQ. 5 .OR. ITYPE.EQ. 9) THEN
         JTYPE= 2
      ELSEIF(ITYPE.EQ. 4 .OR. ITYPE.EQ.13) THEN
         JTYPE= 3
      ELSEIF(ITYPE.EQ. 6 .OR. ITYPE.EQ. 7 .OR.
     2       ITYPE.EQ.10 .OR. ITYPE.EQ.11) THEN
         JTYPE= 4
      ELSEIF(ITYPE.EQ. 8 .OR. ITYPE.EQ.12 .OR.
     2       ITYPE.EQ.14 .OR. ITYPE.EQ.15) THEN
         JTYPE= 5
      ELSEIF(ITYPE.EQ.16) THEN
         JTYPE= 6
      ENDIF
C
C EMPTY INTEGRAL SUMMATION STORAGE
C
      IKL= 0
      IF(JTYPE.EQ. 1) THEN
         R00(1)= ZER
      ELSEIF(JTYPE.EQ. 2) THEN
         CALL INTK2(IKL)
      ELSEIF(JTYPE.EQ. 3) THEN
         CALL INTK3(IKL)
      ELSEIF(JTYPE.EQ. 4) THEN
         CALL INTK4(IKL)
      ELSEIF(JTYPE.EQ. 5) THEN
         CALL INTK5(IKL)
      ELSEIF(JTYPE.EQ. 6) THEN
         CALL INTK6(IKL)
      ENDIF
C
C EMPTY COMMON GOUT
C
C---  NGOUT=256
C---  CALL VCLR(GPOPLE(1+NORG),1,NGOUT)
C
C
C FIND DIRECTION COSINES OF PENULTIMATE AXES FROM COORDINATES OF AB
C P(1,1),P(1,2),... ARE DIRECTION COSINES OF AXES AT P.  Z-AXIS ALONG AB
C T(1),T(2),T(3)... ARE DIRECTION COSINES OF AXES AT Q.  Z-AXIS ALONG CD
C
C FIND DIRECTION COSINES OF AB AND CD. THESE ARE LOCAL Z-AXES.
C IF INDETERMINATE TAKE ALONG SPACE Z-AXIS
C
      P(1,3)= ZER
      P(2,3)= ZER
      P(3,3)= ONE
      RAB= ZER
      IF(R12.NE.ZER) THEN
         RAB= SQRT(R12)
         TMP= ONE/RAB
         P(1,3)= P12(1,3)*TMP
         P(2,3)= P12(2,3)*TMP
         P(3,3)= P12(3,3)*TMP
      ENDIF
C
      T(1)= ZER
      T(2)= ZER
      T(3)= ONE
      RCD= ZER
      IF(R34.NE.ZER) THEN
         RCD= SQRT(R34)
         TMP= ONE/RCD
         T(1)= P34(1,3)*TMP
         T(2)= P34(2,3)*TMP
         T(3)= P34(3,3)*TMP
      ENDIF
C
C FIND LOCAL Y-AXIS AS COMMON PERPENDICULAR TO AB AND CD
C IF INDETERMINATE TAKE PERPENDICULAR TO AB AND SPACE Z-AXIS
C IF STILL INDETERMINATE TAKE PERPENDICULAR TO AB AND SPACE X-AXIS
C
      COSG= T(1)*P(1,3)+T(2)*P(2,3)+T(3)*P(3,3)
      COSG= MIN( ONE,COSG)
      COSG= MAX(-ONE,COSG)
C     SING= SQRT(ONE-COSG*COSG)
C
C MODIFIED ROTATION TESTING.
C THIS FIX CURES THE SMALL ANGLE PROBLEM.
C
      P(1,2)= T(3)*P(2,3)-T(2)*P(3,3)
      P(2,2)= T(1)*P(3,3)-T(3)*P(1,3)
      P(3,2)= T(2)*P(1,3)-T(1)*P(2,3)
      IF( ABS(COSG).GT.PT9) THEN
         SING= SQRT(P(1,2)*P(1,2)+P(2,2)*P(2,2)+P(3,2)*P(3,2))
      ELSE
         SING= SQRT(ONE-COSG*COSG)
      ENDIF
      IF( ABS(COSG).LE.PT9 .OR. SING.GE.TENM12) THEN
         TMP= ONE/SING
         P(1,2)= P(1,2)*TMP
         P(2,2)= P(2,2)*TMP
         P(3,2)= P(3,2)*TMP
      ELSE
         I=3
         IF( ABS(P(1,3)).LE.PT7) I=1
         TMP = P(I,3)*P(I,3)
         TMP = MIN( ONE,TMP)
         TMP = SQRT(ONE-TMP)
         IF(TMP.NE.ZER) TMP= ONE/TMP
         IF( ABS(P(1,3)).LE.PT7) THEN
            P(1,2)= ZER
            P(2,2)= P(3,3)*TMP
            P(3,2)=-P(2,3)*TMP
         ELSE
            P(1,2)= P(2,3)*TMP
            P(2,2)=-P(1,3)*TMP
            P(3,2)= ZER
         ENDIF
      ENDIF
C
C FIND DIRECTION COSINES OF LOCAL X-AXES
C
      P(1,1)= P(2,2)*P(3,3)-P(3,2)*P(2,3)
      P(2,1)= P(3,2)*P(1,3)-P(1,2)*P(3,3)
      P(3,1)= P(1,2)*P(2,3)-P(2,2)*P(1,3)
C
C FIND COORDINATES OF C RELATIVE TO LOCAL AXES AT A
C
      T(1)= P34(1,1)-P12(1,1)
      T(2)= P34(2,1)-P12(2,1)
      T(3)= P34(3,1)-P12(3,1)
      ACX = T(1)*P(1,1)+T(2)*P(2,1)+T(3)*P(3,1)
      ACY = T(1)*P(1,2)+T(2)*P(2,2)+T(3)*P(3,2)
      ACZ = T(1)*P(1,3)+T(2)*P(2,3)+T(3)*P(3,3)
C
C SET ACY= 0  IF CLOSE
C
      IF( ABS(ACY).LE.ACYCUT) THEN
         ACY = ZER
         ACY2= ZER
      ELSE
         ACY2= ACY*ACY
      ENDIF
C
C DIRECTION COSINES OF CD LOCAL AXES WITH RESPECT TO AB LOCAL AXES
C (COSG,0,-SING)  (0,1,0)  (SING,0,COSG)
C
C PRELIMINARY P LOOP
C
C FILL GEOMPQ WITH INFORMATION ABOUT P IN PRELIMINARY P-LOOP
C
      JI= 1
      DO 170 I=1,NGA
         X01= EXA(I)
         DO 170 J=1,NGB
            X02= EXB(J)
            X12= X01+X02
            X21= ONE/X12
            Y01= X01*X21
            Y02= ONE-Y01
            Y12= Y01*X02
            TX12(JI)= X12
            TX21(JI)= X21*PT5
            TY02(JI)= Y02*RAB
CJMS        TY01(JI)=-Y01*RAB
            TY01(JI)= TY02(JI)-RAB
C           IF(JTYPE.EQ.4 .OR. JTYPE.EQ.5) TY01(JI)= TY01(JI)*X12
            R12Y12= R12*Y12
            IF(R12Y12.GT.CUX) THEN
               ISMLP(JI)=2
               GO TO 160
            ENDIF
            E12= X21* EXP(-R12Y12)
            TST= E12*CMAXA(I)*CMAXB(J)
            ISMLP(JI)=0
            IF(TST.LE.EI1) ISMLP(JI)=1
            IF(TST.LE.EI2) ISMLP(JI)=2
            E12= PIF*E12
C
C FOR TYPES 0000,0001,0011 ONLY D00P NEEDED
C
            D00P(JI)= E12*CSA(I)*CSB(J)
            IF(JTYPE.GT.3) THEN
               D01P(JI)= E12*CSA(I)*CPB(J)
               IF(JTYPE.LE.5) THEN
                  IF(D01P(JI).NE.ZER) THEN
                     D00P(JI)= D00P(JI)/D01P(JI)
                  ENDIF
               ELSE
                  D10P(JI)= E12*CPA(I)*CSB(J)
                  D11P(JI)= E12*CPA(I)*CPB(J)
                  IF(D11P(JI).NE.ZER) THEN
                     TMP = ONE/D11P(JI)
                     D00P(JI)= D00P(JI)*TMP
                     D01P(JI)= D01P(JI)*TMP
                     D10P(JI)= D10P(JI)*TMP
                  ENDIF
               ENDIF
            ENDIF
  160       CONTINUE
  170 JI=JI+1
C
C BEGIN Q LOOP
C
      IKL= 0
      DO 190 K=1,NGC
         X03= EXC(K)
         DO 190 L=1,NGD
            X04= EXD(L)
            X34= X03+X04
            X43= ONE/X34
            Y03= X03*X43
            Y04= ONE-Y03
            Y34= Y03*X04
            R34Y34= R34*Y34
            IF(R34Y34.GT.CUX) GO TO 180
            E34= X43* EXP(-R34Y34)
            TST= E34*CMAXC(K)*CMAXD(L)
            IF(TST.LE.EI2) GO TO 180
            ISMLQ= 0
            IF(TST.LE.EI1) ISMLQ= 1
C
C CQX = COMPONENT OF CQ ALONG PENULTIMATE X-AXIS
C CQZ = COMPONENT OF CQ ALONG PENULTIMATE Z-AXIS
C
            CQ = RCD*Y04
            CQX= CQ*SING
            CQZ= CQ*COSG
C
C FIND COORDINATES OF Q RELATIVE TO AXES AT A
C QPR IS PERPENDICULAR FROM Q TO AB
C
            AQX= ACX+CQX
            AQX2=AQX*AQX
            AQXY=AQX*ACY
            AQZ= ACZ+CQZ
            QPS= AQX2+ACY2
C
            IF(JTYPE.NE.1) THEN
               SQ(0,0)= E34*CSC(K)*CSD(L)
               SQ(1,0)= E34*CSC(K)*CPD(L)
               SQ(0,1)= E34*CPC(K)*CSD(L)
               SQ(1,1)= E34*CPC(K)*CPD(L)
            ELSE
               SQ(0,0)= E34*CSC(K)*CSD(L)
            ENDIF
C
C USE SPECIAL FAST ROUTINE FOR INNER LOOPS FOR 0000 ... 1111
C
CJMS  ZEROING OF THE FQx (x=0..4) ARRAYS OF LABELLED COMMON /FQ04/
CJMS  TAKES PLACE IN THE INTJx (x=1..6) SUBROUTINES
C
CJMS  ZEROING OF THE R0x (x=0..4) ARRAYS OF LABELLED COMMON /KI3 /
CJMS  TOOK PLACE IN THE CALL (FOR IKL=0) OF THE INTKx (x=2..6)
CJMS  SUBROUTINES
C
            IKL= IKL+1
            CALL SP0S1S(JTYPE,IKL)
C
  180       CONTINUE
  190 CONTINUE
C
      IG= 1+NORG
      IF(JTYPE.EQ.1) THEN
         GPOPLE(IG)= R00(1)
         GO TO 999
      ENDIF
C
      QX= RCD*SING
      QZ= RCD*COSG
      IF(JTYPE.EQ. 2)THEN
         CALL MCDV2(GPOPLE(IG),QX,QZ)
      ELSEIF(JTYPE.EQ. 3) THEN
         CALL MCDV3(GPOPLE(IG),QX,QZ)
      ELSEIF(JTYPE.EQ. 4) THEN
         CALL MCDV4(GPOPLE(IG),QX,QZ)
      ELSEIF(JTYPE.EQ. 5) THEN
         CALL MCDV5(GPOPLE(IG),QX,QZ)
      ELSEIF(JTYPE.EQ. 6) THEN
         CALL MCDV6(GPOPLE(IG),QX,QZ)
      ENDIF
C
CJMS  NOW, THE TRANSPOSE OF P TO BE USED FOR COMPUTATIONAL EFFICIENCY
C
      DO 195 J=1,2
         DO 195 I=J+1,3
            TMP= P(I,J)
            P(I,J)= P(J,I)
            P(J,I)= TMP
  195 CONTINUE
C
      CALL R30S1S(JTYPE,GPOPLE(IG),P)
C
  999 CONTINUE
      RETURN
      END


C*MODULE EFTEI_GENR70   *DECK SETSP_SHELLEFP
C>    @brief  prepare EFP shells
C>    @author  Peng Xu and Tosaporn Sattasathuchana 
C>             Jan 2020
C>    @detail prepare EFP shells for GENR70_QMEFP subroutine
C>
C>    @param  ISHELL: index of the shell
C>    @param  LAT: type of basis function 
C>    @param  INEW: updated shell index
C>    @param  LA: updated basis type
C>    @param  NGA: number of gaussian functions in a shell
C>    @param  CMAXA: maximum coefficent 
C>    @param  EXA: exponenent of basis function
C>    @param  CSA: coefficient of s basis function
C>    @param  CPA: coefficient of p basis function
C>    @param  CRDI: coordinates of atoms associated with the shell
C>    @param  JM: index for the type of EFP potentials
C>    @param  JNAT: counter for number of EFP atoms 
      SUBROUTINE SETSP_SHELLEFP(ISHELL,LAT,INEW,LA,NGA,
     *                          CMAXA,EXA,CSA,CPA,CRDI,JM,JNAT)
      use mx_limits, only: mxgtot,mxsh,mxgsh,mxg2,MXSH
      USE comm_EFPBAS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION CRDI(MXSH,3)
      DIMENSION EXA(MXGSH),CSA(MXGSH),CPA(MXGSH),CMAXA(MXGSH)

      INEW = ISHELL
      LA   = LAT
C NUMBERS OF GAUSSIAN FUNCTIONS IN SHELLS INEW
      NGA  = KNGEF(INEW,JM)
C STARTING LOCATIONS OF SHELLS INEW
      I    = KSTREF(INEW,JM)-1
C ATOM COUNTER ASSOCIATED WITH SHELLS INEW
      ICC  = KATMEF(INEW,JM)

C LOOP OVER GAUSSIANS IN SHELL INEW
      DO NI=1,NGA
         N=I+NI
         CMAXA(NI)= CMAXEF(N,JM)
         EXA(NI)= EXEF(N,JM)
         CSA(NI)= CSEF(N,JM)
         CPA(NI)= CPEF(N,JM)
      ENDDO

C COORDINATES OF ATOMS ASSOCIATED WITH SHELLS INEW
      DO N=1,3
         CRDI(INEW,N) = PRCORD(N,ICC+JNAT)
      ENDDO

      RETURN
      END

C*MODULE EFTEI_GENR70   *DECK SETSP_SHELLAI
C>    @brief  prepare QM shells
C>    @author  Peng Xu and Tosaporn Sattasathuchana 
C>             Jan 2020
C>    @detail prepare QM shells for GENR70_QMEFP subroutine
C>
C>    @param  JSHELL: index of the shell
C>    @param  LBT: type of basis function 
C>    @param  JNEW: updated shell index
C>    @param  LB: updated basis type
C>    @param  NGB: number of gaussian functions in a shell
C>    @param  CMAX: maximum coefficent 
C>    @param  CMAXB: updated maximum coefficient of each shell
C>    @param  EXB: exponenent of basis function
C>    @param  CSB: coefficient of s basis function
C>    @param  CPB: coefficient of p basis function
C>    @param  CRDJ: coordinates of atoms associated with the shell

      SUBROUTINE SETSP_SHELLAI(JSHELL,LBT,JNEW,LB,CMAX,
     *                         NGB,CMAXB,EXB,CSB,CPB,CRDJ)
      use mx_limits, only: mxgtot,mxsh,mxgsh,mxg2,mxatm
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION CMAX(MXGTOT)
      DIMENSION CRDJ(MXSH,3)
      DIMENSION EXB(MXGSH),CSB(MXGSH),CPB(MXGSH),CMAXB(MXGSH)

      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /NSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),
     2                CF(MXGTOT),CG(MXGTOT),CH(MXGTOT),CI(MXGTOT),
     3                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),
     4                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
C
      JNEW = JSHELL
      LB   = LBT
C NUMBERS OF GAUSSIAN FUNCTIONS IN SHELLS JNEW
      NGB  = KNG(JNEW)
C STARTING LOCATIONS OF SHELLS JNEW
      J    = KSTART(JNEW)-1

C LOOP OVER GAUSSIANS IN SHELL JNEW
      DO NJ=1,NGB
         N=J+NJ
         CMAXB(NJ)= CMAX(N)
         EXB(NJ)= EX(N)
         CSB(NJ)= CS(N)
         CPB(NJ)= CP(N)
      ENDDO

C ATOM COUNTER ASSOCIATED WITH SHELLS JNEW
      JATM=KATOM(JNEW)
C COORDINATES OF ATOMS ASSOCIATED WITH SHELLS JNEW
      DO N=1,3
         CRDJ(JNEW,N)= C(N,JATM)
      ENDDO

      RETURN
      END


C*MODULE EFTEI_GENR70   *DECK SETUP_COORD 
C>    @brief  prepare atomic coordinates and distances associated
C>            with shells 
C>    @author  Peng Xu and Tosaporn Sattasathuchana 
C>             Jan 2020
C>    @param  CRDI,CRDJ,CRDK,CRDL: coordinates of atoms associated with
C>                                 shells I,J,K and L
C>    @param  INEW,JNEW,KNEW,LNEW: shell indices
C>    @param  P12: coordinates of atoms associated with I&J shells and
C>                 the distance between them
C>    @param  P34: coordinates of atoms associated with K&L shells and
C>                 the distance between them
C>    @param  R12: R^2, where R is the distance between atoms associated
C>                 with I and J shells 
C>    @param  R34: R^2, where R is the distance between atoms associated
C>                 with K and L shells 

      SUBROUTINE SETUP_COORD(CRDI,CRDJ,CRDK,CRDL,INEW,JNEW,KNEW,LNEW,
     &                       P12,P34,R12,R34)
      use mx_limits, only:mxsh
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DIMENSION P12(3,3),P34(3,3)
      DIMENSION CRDI(MXSH,3),CRDJ(MXSH,3),CRDK(MXSH,3),CRDL(MXSH,3)

      DO N=1,3
         P12(N,1)= CRDI(INEW,N)
         P12(N,2)= CRDJ(JNEW,N)
         P12(N,3)= P12(N,2)-P12(N,1)
      R12= R12+P12(N,3)*P12(N,3)
         P34(N,1)= CRDK(KNEW,N)
         P34(N,2)= CRDL(LNEW,N)
         P34(N,3)= P34(N,2)-P34(N,1)
      R34= R34+P34(N,3)*P34(N,3)
      ENDDO
      
      RETURN
      END

C*MODULE EFTEI_GENR70   *DECK GENR70_QMEFP_INTYPE_1 
C>    @brief  rotated axis (sp) for case 1 (IJ|KL) 
C>    @author  Peng Xu and Tosaporn Sattasathuchana 
C>             Jan 2020
C>    @detail Case 1 integral includes types:
C>            0000,0001,0011,0101,0111,1111
C>            no switching of I, J, K, L 
C>
C>    @param  MODE: different type of TEI involve in QM/EFP EXREP
C>                 MODE=2 two-center on QM, two-center on EFP
C>                 MODE=3 three-center on QM, One-center on EFP
C>    @param   JM: index for the type of EFP potentials
C>    @param   JNAT: counter for number of EFP atoms
C>
C>    @param  LAT: type of basis frunction for shell I
C>    @param  LBT: type of basis function  for shell J
C>    @param  LCT: type of basis function  for shell K
C>    @param  LDT: type of basis function  for shell L
C>    @param  P12: coordinates of atoms associated with I&J shells and
C>                 the distance between them
C>    @param  P34: coordinates of atoms associated with K&L shells and
C>                 the distance between them
C>    @param  R34: R^2, where R is the distance between atoms associated
C>                 with K and L shells 

      SUBROUTINE GENR70_QMEFP_INTYPE_1(MODE,JM,JNAT,
     *                         LAT,LBT,LCT,LDT,P12,P34,R34)
      use mx_limits, only: mxgtot,mxsh,mxgsh,mxg2
      USE comm_EFPBAS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      COMMON /GEOMPQ/ R12,RAB,X34,X43,AQZ,QPR,QPS,
     2                TX12(MXG2),TX21(MXG2),TY01(MXG2),TY02(MXG2),
     3                D00P(MXG2),D01P(MXG2),D10P(MXG2),D11P(MXG2),
     4                NGANGB
      COMMON /INTAC2/ EI1,EI2,CUX
      COMMON /MAXC  / CMAX(MXGTOT),CMAXA(MXGSH),CMAXB(MXGSH),
     2                CMAXC(MXGSH),CMAXD(MXGSH),ISMLP(MXG2),ISMLQ
      COMMON /NSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),
     2                CF(MXGTOT),CG(MXGTOT),CH(MXGTOT),CI(MXGTOT),
     3                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),
     4                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
      COMMON /POPOUT/ LPOPI,LPOPJ,LPOPK,LPOPL
      COMMON /SHLG70/ ISHELL,JSHELL,KSHELL,LSHELL,INEW,JNEW,KNEW,LNEW
      COMMON /SHLLFO/ NGA,LA,EXA(MXGSH),CSA(MXGSH),CPA(MXGSH),
     2                NGB,LB,EXB(MXGSH),CSB(MXGSH),CPB(MXGSH),
     3                NGC,LC,EXC(MXGSH),CSC(MXGSH),CPC(MXGSH),
     4                NGD,LD,EXD(MXGSH),CSD(MXGSH),CPD(MXGSH)
C
      DIMENSION P12(3,3),P34(3,3)
      DIMENSION CRDI(MXSH,3),CRDJ(MXSH,3),CRDK(MXSH,3),CRDL(MXSH,3)
C
C
C TYPES 0000,0001,0011,0101,0111,1111 ARE UNALTERED
C
      LPOPI = 64
      LPOPJ = 16
      LPOPK =  4
      LPOPL =  1

      CALL SETSP_SHELLEFP(ISHELL,LAT,INEW,LA,NGA,
     *                    CMAXA,EXA,CSA,CPA,CRDI,JM,JNAT)
      CALL SETSP_SHELLAI(JSHELL,LBT,JNEW,LB,CMAX,
     *                   NGB,CMAXB,EXB,CSB,CPB,CRDJ)
      IF(MODE.EQ.2) THEN
        CALL SETSP_SHELLEFP(KSHELL,LCT,KNEW,LC,NGC,
     *                      CMAXC,EXC,CSC,CPC,CRDK,JM,JNAT)
      ELSE
        CALL SETSP_SHELLAI(KSHELL,LCT,KNEW,LC,CMAX,
     *                     NGC,CMAXC,EXC,CSC,CPC,CRDK)
      ENDIF
      CALL SETSP_SHELLAI(LSHELL,LDT,LNEW,LD,CMAX,
     *                   NGD,CMAXD,EXD,CSD,CPD,CRDL)

      NGANGB = NGA * NGB

      CALL SETUP_COORD(CRDI,CRDJ,CRDK,CRDL,INEW,JNEW,KNEW,LNEW,
     &                       P12,P34,R12,R34)    
      RETURN
      END


C*MODULE EFTEI_GENR70   *DECK GENR70_QMEFP_INTYPE_2 
C>    @brief  rotated axis (sp) for case 2 (JI|KL) 
C>    @author  Peng Xu and Tosaporn Sattasathuchana 
C>             Jan 2020
C>    @detail Case 2 integral includes types:
C>            1001,1011 (IJ SWITCHED)
C>    @see    parameters same as GENR70_QMEFP_INTYPE_1

      SUBROUTINE GENR70_QMEFP_INTYPE_2(MODE,JM,JNAT,
     *                         LAT,LBT,LCT,LDT,P12,P34,R34)
      use mx_limits, only: mxgtot,mxsh,mxgsh,mxg2
      USE comm_EFPBAS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      COMMON /GEOMPQ/ R12,RAB,X34,X43,AQZ,QPR,QPS,
     2                TX12(MXG2),TX21(MXG2),TY01(MXG2),TY02(MXG2),
     3                D00P(MXG2),D01P(MXG2),D10P(MXG2),D11P(MXG2),
     4                NGANGB
      COMMON /INTAC2/ EI1,EI2,CUX
      COMMON /MAXC  / CMAX(MXGTOT),CMAXA(MXGSH),CMAXB(MXGSH),
     2                CMAXC(MXGSH),CMAXD(MXGSH),ISMLP(MXG2),ISMLQ
      COMMON /NSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),
     2                CF(MXGTOT),CG(MXGTOT),CH(MXGTOT),CI(MXGTOT),
     3                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),
     4                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
      COMMON /POPOUT/ LPOPI,LPOPJ,LPOPK,LPOPL
      COMMON /SHLG70/ ISHELL,JSHELL,KSHELL,LSHELL,INEW,JNEW,KNEW,LNEW
      COMMON /SHLLFO/ NGA,LA,EXA(MXGSH),CSA(MXGSH),CPA(MXGSH),
     2                NGB,LB,EXB(MXGSH),CSB(MXGSH),CPB(MXGSH),
     3                NGC,LC,EXC(MXGSH),CSC(MXGSH),CPC(MXGSH),
     4                NGD,LD,EXD(MXGSH),CSD(MXGSH),CPD(MXGSH)
C
      DIMENSION P12(3,3),P34(3,3)
      DIMENSION CRDI(MXSH,3),CRDJ(MXSH,3),CRDK(MXSH,3),CRDL(MXSH,3)
C
C
C TYPES 1001,1011 HAVE IJ SWITCHED
C
      LPOPI = 16
      LPOPJ = 64
      LPOPK =  4
      LPOPL =  1

      CALL SETSP_SHELLAI(JSHELL,LBT,INEW,LA,CMAX,
     *                   NGA,CMAXA,EXA,CSA,CPA,CRDI)


      CALL SETSP_SHELLEFP(ISHELL,LAT,JNEW,LB,NGB,
     *                    CMAXB,EXB,CSB,CPB,CRDJ,JM,JNAT)

      IF(MODE.EQ.2) THEN
        CALL SETSP_SHELLEFP(KSHELL,LCT,KNEW,LC,NGC,
     *                      CMAXC,EXC,CSC,CPC,CRDK,JM,JNAT)

      ELSE
        CALL SETSP_SHELLAI(KSHELL,LCT,KNEW,LC,CMAX,
     *                     NGC,CMAXC,EXC,CSC,CPC,CRDK)

      ENDIF

      CALL SETSP_SHELLAI(LSHELL,LDT,LNEW,LD,CMAX,
     *                   NGD,CMAXD,EXD,CSD,CPD,CRDL)


      NGANGB = NGA * NGB

      CALL SETUP_COORD(CRDI,CRDJ,CRDK,CRDL,INEW,JNEW,KNEW,LNEW,
     &                       P12,P34,R12,R34)

      RETURN
      END

C*MODULE EFTEI_GENR70   *DECK GENR70_QMEFP_INTYPE_3 
C>    @brief  rotated axis (sp) for case 3 (IJ|LK) 
C>    @author  Peng Xu and Tosaporn Sattasathuchana 
C>             Jan 2020
C>    @detail Case 3 integral includes types:
C>            0010,0110 (KL SWITCHED)
C>    @see    parameters same as GENR70_QMEFP_INTYPE_1

      SUBROUTINE GENR70_QMEFP_INTYPE_3(MODE,JM,JNAT,
     *                         LAT,LBT,LCT,LDT,P12,P34,R34)
      use mx_limits, only: mxgtot,mxsh,mxgsh,mxg2
      USE comm_EFPBAS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      COMMON /GEOMPQ/ R12,RAB,X34,X43,AQZ,QPR,QPS,
     2                TX12(MXG2),TX21(MXG2),TY01(MXG2),TY02(MXG2),
     3                D00P(MXG2),D01P(MXG2),D10P(MXG2),D11P(MXG2),
     4                NGANGB
      COMMON /INTAC2/ EI1,EI2,CUX
      COMMON /MAXC  / CMAX(MXGTOT),CMAXA(MXGSH),CMAXB(MXGSH),
     2                CMAXC(MXGSH),CMAXD(MXGSH),ISMLP(MXG2),ISMLQ
      COMMON /NSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),
     2                CF(MXGTOT),CG(MXGTOT),CH(MXGTOT),CI(MXGTOT),
     3                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),
     4                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
      COMMON /POPOUT/ LPOPI,LPOPJ,LPOPK,LPOPL
      COMMON /SHLG70/ ISHELL,JSHELL,KSHELL,LSHELL,INEW,JNEW,KNEW,LNEW
      COMMON /SHLLFO/ NGA,LA,EXA(MXGSH),CSA(MXGSH),CPA(MXGSH),
     2                NGB,LB,EXB(MXGSH),CSB(MXGSH),CPB(MXGSH),
     3                NGC,LC,EXC(MXGSH),CSC(MXGSH),CPC(MXGSH),
     4                NGD,LD,EXD(MXGSH),CSD(MXGSH),CPD(MXGSH)
C
      DIMENSION P12(3,3),P34(3,3)
      DIMENSION CRDI(MXSH,3),CRDJ(MXSH,3),CRDK(MXSH,3),CRDL(MXSH,3)
C
C
C TYPES 0010,0110 HAVE KL SWITCHED
C
      LPOPI = 64
      LPOPJ = 16
      LPOPK =  1
      LPOPL =  4

      CALL SETSP_SHELLEFP(ISHELL,LAT,INEW,LA,NGA,
     *                    CMAXA,EXA,CSA,CPA,CRDI,JM,JNAT)

      CALL SETSP_SHELLAI(JSHELL,LBT,JNEW,LB,CMAX,
     *                   NGB,CMAXB,EXB,CSB,CPB,CRDJ)

      CALL SETSP_SHELLAI(LSHELL,LDT,KNEW,LC,CMAX,
     *                   NGC,CMAXC,EXC,CSC,CPC,CRDK)

      IF(MODE.EQ.2) THEN
        CALL SETSP_SHELLEFP(KSHELL,LCT,LNEW,LD,NGD,
     *                      CMAXD,EXD,CSD,CPD,CRDL,JM,JNAT)
      ELSE
        CALL SETSP_SHELLAI(KSHELL,LCT,LNEW,LD,CMAX,
     *                     NGD,CMAXD,EXD,CSD,CPD,CRDL)
      ENDIF

      NGANGB = NGA * NGB

      CALL SETUP_COORD(CRDI,CRDJ,CRDK,CRDL,INEW,JNEW,KNEW,LNEW,
     &                       P12,P34,R12,R34)

      RETURN
      END

C*MODULE EFTEI_GENR70   *DECK GENR70_QMEFP_INTYPE_4 
C>    @brief  rotated axis (sp) for case 4 (KL|IJ)
C>    @author  Peng Xu and Tosaporn Sattasathuchana 
C>             Jan 2020
C>    @detail Case 4 integral includes types:
C>            0100,1100,1101 (pairs IJ and KL SWITCHED)
C>    @see    parameters same as GENR70_QMEFP_INTYPE_1

      SUBROUTINE GENR70_QMEFP_INTYPE_4(MODE,JM,JNAT,
     *                         LAT,LBT,LCT,LDT,P12,P34,R34)
      use mx_limits, only: mxgtot,mxsh,mxgsh,mxg2
      USE comm_EFPBAS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      COMMON /GEOMPQ/ R12,RAB,X34,X43,AQZ,QPR,QPS,
     2                TX12(MXG2),TX21(MXG2),TY01(MXG2),TY02(MXG2),
     3                D00P(MXG2),D01P(MXG2),D10P(MXG2),D11P(MXG2),
     4                NGANGB
      COMMON /INTAC2/ EI1,EI2,CUX
      COMMON /MAXC  / CMAX(MXGTOT),CMAXA(MXGSH),CMAXB(MXGSH),
     2                CMAXC(MXGSH),CMAXD(MXGSH),ISMLP(MXG2),ISMLQ
      COMMON /NSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),
     2                CF(MXGTOT),CG(MXGTOT),CH(MXGTOT),CI(MXGTOT),
     3                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),
     4                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
      COMMON /POPOUT/ LPOPI,LPOPJ,LPOPK,LPOPL
      COMMON /SHLG70/ ISHELL,JSHELL,KSHELL,LSHELL,INEW,JNEW,KNEW,LNEW
      COMMON /SHLLFO/ NGA,LA,EXA(MXGSH),CSA(MXGSH),CPA(MXGSH),
     2                NGB,LB,EXB(MXGSH),CSB(MXGSH),CPB(MXGSH),
     3                NGC,LC,EXC(MXGSH),CSC(MXGSH),CPC(MXGSH),
     4                NGD,LD,EXD(MXGSH),CSD(MXGSH),CPD(MXGSH)
C
      DIMENSION P12(3,3),P34(3,3)
      DIMENSION CRDI(MXSH,3),CRDJ(MXSH,3),CRDK(MXSH,3),CRDL(MXSH,3)
C
C
C TYPES 0100,1100,1101 HAVE PAIRS IJ AND KL SWITCHED
C
      LPOPI =  4
      LPOPJ =  1
      LPOPK = 64
      LPOPL = 16

      IF(MODE.EQ.2) THEN
        CALL SETSP_SHELLEFP(KSHELL,LCT,INEW,LA,NGA,
     *                      CMAXA,EXA,CSA,CPA,CRDI,JM,JNAT)
      ELSE
        CALL SETSP_SHELLAI(KSHELL,LCT,INEW,LA,CMAX,
     *                     NGA,CMAXA,EXA,CSA,CPA,CRDI)
      ENDIF

      CALL SETSP_SHELLAI(LSHELL,LDT,JNEW,LB,CMAX,
     *                   NGB,CMAXB,EXB,CSB,CPB,CRDJ)

      CALL SETSP_SHELLEFP(ISHELL,LAT,KNEW,LC,NGC,
     *                    CMAXC,EXC,CSC,CPC,CRDK,JM,JNAT)

      CALL SETSP_SHELLAI(JSHELL,LBT,LNEW,LD,CMAX,
     *                   NGD,CMAXD,EXD,CSD,CPD,CRDL)


      NGANGB = NGA * NGB

      CALL SETUP_COORD(CRDI,CRDJ,CRDK,CRDL,INEW,JNEW,KNEW,LNEW,
     &                       P12,P34,R12,R34)

      RETURN
      END


C*MODULE EFTEI_GENR70   *DECK GENR70_QMEFP_INTYPE_5 
C>    @brief  rotated axis (sp) for case 5 (JI|LK) 
C>    @author  Peng Xu and Tosaporn Sattasathuchana 
C>             Jan 2020
C>    @detail Case 5 integral includes types:
C>            1010 (IJ SWITCHED and KL SWITCHED)
C>    @see    parameters same as GENR70_QMEFP_INTYPE_1

      SUBROUTINE GENR70_QMEFP_INTYPE_5(MODE,JM,JNAT,
     *                         LAT,LBT,LCT,LDT,P12,P34,R34)
      use mx_limits, only: mxgtot,mxsh,mxgsh,mxg2
      USE comm_EFPBAS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      COMMON /GEOMPQ/ R12,RAB,X34,X43,AQZ,QPR,QPS,
     2                TX12(MXG2),TX21(MXG2),TY01(MXG2),TY02(MXG2),
     3                D00P(MXG2),D01P(MXG2),D10P(MXG2),D11P(MXG2),
     4                NGANGB
      COMMON /INTAC2/ EI1,EI2,CUX
      COMMON /MAXC  / CMAX(MXGTOT),CMAXA(MXGSH),CMAXB(MXGSH),
     2                CMAXC(MXGSH),CMAXD(MXGSH),ISMLP(MXG2),ISMLQ
      COMMON /NSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),
     2                CF(MXGTOT),CG(MXGTOT),CH(MXGTOT),CI(MXGTOT),
     3                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),
     4                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
      COMMON /POPOUT/ LPOPI,LPOPJ,LPOPK,LPOPL
      COMMON /SHLG70/ ISHELL,JSHELL,KSHELL,LSHELL,INEW,JNEW,KNEW,LNEW
      COMMON /SHLLFO/ NGA,LA,EXA(MXGSH),CSA(MXGSH),CPA(MXGSH),
     2                NGB,LB,EXB(MXGSH),CSB(MXGSH),CPB(MXGSH),
     3                NGC,LC,EXC(MXGSH),CSC(MXGSH),CPC(MXGSH),
     4                NGD,LD,EXD(MXGSH),CSD(MXGSH),CPD(MXGSH)
C
      DIMENSION P12(3,3),P34(3,3)
      DIMENSION CRDI(MXSH,3),CRDJ(MXSH,3),CRDK(MXSH,3),CRDL(MXSH,3)
C
C
C TYPE 1010 HAS IJ SWITCHED AND KL SWITCHED
C
      LPOPI = 16
      LPOPJ = 64
      LPOPK =  1
      LPOPL =  4

      CALL SETSP_SHELLAI(JSHELL,LBT,INEW,LA,CMAX,
     *                   NGA,CMAXA,EXA,CSA,CPA,CRDI)

      CALL SETSP_SHELLEFP(ISHELL,LAT,JNEW,LB,NGB,
     *                    CMAXB,EXB,CSB,CPB,CRDJ,JM,JNAT)

      CALL SETSP_SHELLAI(LSHELL,LDT,KNEW,LC,CMAX,
     *                   NGC,CMAXC,EXC,CSC,CPC,CRDK)

      IF(MODE.EQ.2) THEN
        CALL SETSP_SHELLEFP(KSHELL,LCT,LNEW,LD,NGD,
     *                      CMAXD,EXD,CSD,CPD,CRDL,JM,JNAT)
      ELSE
        CALL SETSP_SHELLAI(KSHELL,LCT,LNEW,LD,CMAX,
     *                     NGD,CMAXD,EXD,CSD,CPD,CRDL)
      ENDIF

      NGANGB = NGA * NGB

      CALL SETUP_COORD(CRDI,CRDJ,CRDK,CRDL,INEW,JNEW,KNEW,LNEW,
     &                       P12,P34,R12,R34)

      RETURN
      END

C*MODULE EFTEI_GENR70   *DECK GENR70_QMEFP_INTYPE_6 
C>    @brief  rotated axis (sp) for case 6 (KL|JI) 
C>    @author  Peng Xu and Tosaporn Sattasathuchana 
C>             Jan 2020
C>    @detail Case 6 integral includes types:
C>            1000 (PAIRS IJ and KL SWITCHED,followed by KL SWITCH)
C>    @see    parameters same as GENR70_QMEFP_INTYPE_1

      SUBROUTINE GENR70_QMEFP_INTYPE_6(MODE,JM,JNAT,
     *                         LAT,LBT,LCT,LDT,P12,P34,R34)
      use mx_limits, only: mxgtot,mxsh,mxgsh,mxg2
      USE comm_EFPBAS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      COMMON /GEOMPQ/ R12,RAB,X34,X43,AQZ,QPR,QPS,
     2                TX12(MXG2),TX21(MXG2),TY01(MXG2),TY02(MXG2),
     3                D00P(MXG2),D01P(MXG2),D10P(MXG2),D11P(MXG2),
     4                NGANGB
      COMMON /INTAC2/ EI1,EI2,CUX
      COMMON /MAXC  / CMAX(MXGTOT),CMAXA(MXGSH),CMAXB(MXGSH),
     2                CMAXC(MXGSH),CMAXD(MXGSH),ISMLP(MXG2),ISMLQ
      COMMON /NSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),
     2                CF(MXGTOT),CG(MXGTOT),CH(MXGTOT),CI(MXGTOT),
     3                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),
     4                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
      COMMON /POPOUT/ LPOPI,LPOPJ,LPOPK,LPOPL
      COMMON /SHLG70/ ISHELL,JSHELL,KSHELL,LSHELL,INEW,JNEW,KNEW,LNEW
      COMMON /SHLLFO/ NGA,LA,EXA(MXGSH),CSA(MXGSH),CPA(MXGSH),
     2                NGB,LB,EXB(MXGSH),CSB(MXGSH),CPB(MXGSH),
     3                NGC,LC,EXC(MXGSH),CSC(MXGSH),CPC(MXGSH),
     4                NGD,LD,EXD(MXGSH),CSD(MXGSH),CPD(MXGSH)
C
      DIMENSION P12(3,3),P34(3,3)
      DIMENSION CRDI(MXSH,3),CRDJ(MXSH,3),CRDK(MXSH,3),CRDL(MXSH,3)
C
C
C TYPE 1000  HAS PAIRS IJ AND KL SWITCHED FOLLOWED BY KL SWITCH
C
      LPOPI =  1
      LPOPJ =  4
      LPOPK = 64
      LPOPL = 16

      IF(MODE.EQ.2) THEN
        CALL SETSP_SHELLEFP(KSHELL,LCT,INEW,LA,NGA,
     *                      CMAXA,EXA,CSA,CPA,CRDI,JM,JNAT)
      ELSE
        CALL SETSP_SHELLAI(KSHELL,LCT,INEW,LA,CMAX,
     *                     NGA,CMAXA,EXA,CSA,CPA,CRDI)
      ENDIF

      CALL SETSP_SHELLAI(LSHELL,LDT,JNEW,LB,CMAX,
     *                   NGB,CMAXB,EXB,CSB,CPB,CRDJ)

      CALL SETSP_SHELLAI(JSHELL,LBT,KNEW,LC,CMAX,
     *                   NGC,CMAXC,EXC,CSC,CPC,CRDK)
 
      CALL SETSP_SHELLEFP(ISHELL,LAT,LNEW,LD,NGD,
     *                    CMAXD,EXD,CSD,CPD,CRDL,JM,JNAT)

      NGANGB = NGA * NGB

      CALL SETUP_COORD(CRDI,CRDJ,CRDK,CRDL,INEW,JNEW,KNEW,LNEW,
     &                       P12,P34,R12,R34)

      RETURN
      END


C*MODULE EFTEI_GENR70   *DECK GENR70_QMEFP_INTYPE_7 
C>    @brief  rotated axis (sp) for case 7 (LK|IJ) 
C>    @author  Peng Xu and Tosaporn Sattasathuchana 
C>             Jan 2020
C>    @detail Case 7 integral includes types:
C>            1110 (PAIRS IJ and KL SWITCHED followed by IJ SWITCH)
C>    @see    parameters same as GENR70_QMEFP_INTYPE_1

      SUBROUTINE GENR70_QMEFP_INTYPE_7(MODE,JM,JNAT,
     *                         LAT,LBT,LCT,LDT,P12,P34,R34)
      use mx_limits, only: mxgtot,mxsh,mxgsh,mxg2
      USE comm_EFPBAS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      COMMON /GEOMPQ/ R12,RAB,X34,X43,AQZ,QPR,QPS,
     2                TX12(MXG2),TX21(MXG2),TY01(MXG2),TY02(MXG2),
     3                D00P(MXG2),D01P(MXG2),D10P(MXG2),D11P(MXG2),
     4                NGANGB
      COMMON /INTAC2/ EI1,EI2,CUX
      COMMON /MAXC  / CMAX(MXGTOT),CMAXA(MXGSH),CMAXB(MXGSH),
     2                CMAXC(MXGSH),CMAXD(MXGSH),ISMLP(MXG2),ISMLQ
      COMMON /NSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),
     2                CF(MXGTOT),CG(MXGTOT),CH(MXGTOT),CI(MXGTOT),
     3                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),
     4                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
      COMMON /POPOUT/ LPOPI,LPOPJ,LPOPK,LPOPL
      COMMON /SHLG70/ ISHELL,JSHELL,KSHELL,LSHELL,INEW,JNEW,KNEW,LNEW
      COMMON /SHLLFO/ NGA,LA,EXA(MXGSH),CSA(MXGSH),CPA(MXGSH),
     2                NGB,LB,EXB(MXGSH),CSB(MXGSH),CPB(MXGSH),
     3                NGC,LC,EXC(MXGSH),CSC(MXGSH),CPC(MXGSH),
     4                NGD,LD,EXD(MXGSH),CSD(MXGSH),CPD(MXGSH)
C
      DIMENSION P12(3,3),P34(3,3)
      DIMENSION CRDI(MXSH,3),CRDJ(MXSH,3),CRDK(MXSH,3),CRDL(MXSH,3)
C
C
C TYPE 1110 HAS PAIRS IJ AND KL SWITCHED FOLLOWED BY IJ SWITCH
C
      LPOPI =  4
      LPOPJ =  1
      LPOPK = 16
      LPOPL = 64

      CALL SETSP_SHELLAI(LSHELL,LDT,INEW,LA,CMAX,
     *                   NGA,CMAXA,EXA,CSA,CPA,CRDI)

      IF(MODE.EQ.2) THEN
        CALL SETSP_SHELLEFP(KSHELL,LCT,JNEW,LB,NGB,
     *                      CMAXB,EXB,CSB,CPB,CRDJ,JM,JNAT)
      ELSE
        CALL SETSP_SHELLAI(KSHELL,LCT,JNEW,LB,CMAX,
     *                     NGB,CMAXB,EXB,CSB,CPB,CRDJ)
      ENDIF

      CALL SETSP_SHELLEFP(ISHELL,LAT,KNEW,LC,NGC,
     *                    CMAXC,EXC,CSC,CPC,CRDK,JM,JNAT)

      CALL SETSP_SHELLAI(JSHELL,LBT,LNEW,LD,CMAX,
     *                   NGD,CMAXD,EXD,CSD,CPD,CRDL)


      NGANGB = NGA * NGB

      CALL SETUP_COORD(CRDI,CRDJ,CRDK,CRDL,INEW,JNEW,KNEW,LNEW,
     &                       P12,P34,R12,R34)

      RETURN
      END

