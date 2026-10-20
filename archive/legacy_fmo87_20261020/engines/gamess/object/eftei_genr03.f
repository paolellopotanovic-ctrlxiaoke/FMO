C*MODULE EFTEI_GENR03  *DECK GENR03_QMEFP 
C>
C>    @brief   rotated axis integration involving s,p,L,d shells for QMEFP
C>
C>    @author  Peng Xu and Tosaporn Sattasathuchana 
C>             Jan 2020
C>    @details This subroutine prepares the necessary shell information
C>             for QM/EFP ERI (rotated axis) calculation and 
C>             is based on the subroutine GENR03 in int2r.src
C>
C>    @param   JM: index for the type of EFP potentials 
C>    @param   JNAT: counter for number of EFP atoms 
C>    @param   MODE: different type of TEI involve in QM/EFP EXREP
C>                 MODE=2 two-center on QM, two-center on EFP
C>                 MODE=3 three-center on QM, One-center on EFP

      SUBROUTINE GENR03_QMEFP(GHONDO,JM,JNAT,MODE)
      use mx_limits, only: mxgtot,mxsh,mxgsh,mxg2
      USE comm_EFPBAS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION GHONDO(*)
C
C
      COMMON /FLIPS / IB(4,3)
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
      COMMON /SHLG70/ ISHELL,JSHELL,KSHELL,LSHELL,INEW,JNEW,KNEW,LNEW
      COMMON /SHLLFO/ NGA,LA,EXA(MXGSH),CSA(MXGSH),CPA(MXGSH),
     2                NGB,LB,EXB(MXGSH),CSB(MXGSH),CPB(MXGSH),
     3                NGC,LC,EXC(MXGSH),CSC(MXGSH),CPC(MXGSH),
     4                NGD,LD,EXD(MXGSH),CSD(MXGSH),CPD(MXGSH)
      COMMON /SHLNOS/ QQ4,LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,
     2                MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL,
     3                NIJ,IJ,KL,IJKL
      COMMON /SHLSPD/ CDA(MXGSH),CDB(MXGSH),CDC(MXGSH),CDD(MXGSH),
     2                D02D(MXG2),D12D(MXG2),D22D(MXG2)
C
C     LABELLED COMMON JMSGYH DEFINED FOR COMPUTATIONAL EFFICIENCY.
C     IT IS ONLY USED IN THIS MODULE INT2R AND IN MODULE INT2B.
C
      COMMON /JMSGYH/ SQ(4)
      COMMON /KI2 / ACY,ACY2,AQX,AQX2,AQXY,Y03,Y04
C
      DIMENSION P12(3,3),P34(3,3),P(3,3),T(3)
C
      PARAMETER (ZER=0.0D+00)
      PARAMETER (PT5=0.5D+00)
      PARAMETER (ONE=1.0D+00)
      PARAMETER (PT7=0.7D+00)
      PARAMETER (PT9=0.9D+00)
      PARAMETER (ACYCUT=1.0D-10)
      PARAMETER (TENM12=1.0D-12)
C
C                PITO52=(PI +PI )*PI * SQRT(PI )
      PARAMETER (PITO52=34.98683665524973D+00)
C
      PIF= PITO52
C
      IEXCH= 1
      LAT= KTYPEF(ISHELL,JM)-1
      LBT= KTYPE(JSHELL)-1

      IF(MODE.EQ.2) THEN
      LCT= KTYPEF(KSHELL,JM)-1
      ELSE ! MODE=3 
      LCT= KTYPE(KSHELL)-1
      ENDIF

      LDT= KTYPE(LSHELL)-1
      ITYPE= 1+LDT+3*(LCT+3*(LBT+3*LAT))
      LTYPE= LAT+LBT+LCT+LDT

      R12= ZER
      R34= ZER
C
      IF(LTYPE.EQ.2) THEN
         JTYPE= 7
         IF(ITYPE.EQ.3) THEN
            GO TO 1010
         ELSEIF(ITYPE.EQ. 7) THEN
            GO TO 1030
         ELSEIF(ITYPE.EQ.19) THEN
            GO TO 1040
         ELSEIF(ITYPE.EQ.55) THEN
            GO TO 1060
         ENDIF
C
      ELSEIF(LTYPE.EQ.3) THEN
         IF(ITYPE.EQ. 6) THEN
            JTYPE= 8
            GO TO 1010
         ELSEIF(ITYPE.EQ. 8) THEN
            JTYPE= 8
            GO TO 1030
         ELSEIF(ITYPE.EQ.46) THEN
            JTYPE= 8
            GO TO 1040
         ELSEIF(ITYPE.EQ.64) THEN
            JTYPE= 8
            GO TO 1060
         ELSEIF(ITYPE.EQ.12) THEN
            JTYPE= 9
            GO TO 1010
         ELSEIF(ITYPE.EQ.30) THEN
            JTYPE= 9
            GO TO 1020
         ELSEIF(ITYPE.EQ.16) THEN
            JTYPE= 9
            GO TO 1030
         ELSEIF(ITYPE.EQ.20) THEN
            JTYPE= 9
            GO TO 1040
         ELSEIF(ITYPE.EQ.34) THEN
            JTYPE= 9
            GO TO 1050
         ELSEIF(ITYPE.EQ.56) THEN
            JTYPE= 9
            GO TO 1060
         ELSEIF(ITYPE.EQ.22) THEN
            JTYPE= 9
            GO TO 1070
         ELSEIF(ITYPE.EQ.58) THEN
            JTYPE= 9
            GO TO 1080
         ENDIF
C
      ELSEIF(LTYPE.EQ.4) THEN
         IF(ITYPE.EQ. 9) THEN
            JTYPE= 10
            GO TO 1010
         ELSEIF(ITYPE.EQ.73) THEN
            JTYPE= 10
            GO TO 1040
         ELSEIF(ITYPE.EQ.15) THEN
            JTYPE= 11
            GO TO 1010
         ELSEIF(ITYPE.EQ.33) THEN
            JTYPE= 11
            GO TO 1020
         ELSEIF(ITYPE.EQ.17) THEN
            JTYPE= 11
            GO TO 1030
         ELSEIF(ITYPE.EQ.47) THEN
            JTYPE= 11
            GO TO 1040
         ELSEIF(ITYPE.EQ.35) THEN
            JTYPE= 11
            GO TO 1050
         ELSEIF(ITYPE.EQ.65) THEN
            JTYPE= 11
            GO TO 1060
         ELSEIF(ITYPE.EQ.49) THEN
            JTYPE= 11
            GO TO 1070
         ELSEIF(ITYPE.EQ.67) THEN
            JTYPE= 11
            GO TO 1080
         ELSEIF(ITYPE.EQ.21) THEN
            JTYPE= 12
            GO TO 1010
         ELSEIF(ITYPE.EQ.57) THEN
            JTYPE= 12
            GO TO 1020
         ELSEIF(ITYPE.EQ.25) THEN
            JTYPE= 12
            GO TO 1030
         ELSEIF(ITYPE.EQ.61) THEN
            JTYPE= 12
            GO TO 1050
         ELSEIF(ITYPE.EQ.39) THEN
            JTYPE= 13
            GO TO 1010
         ELSEIF(ITYPE.EQ.43) THEN
            JTYPE= 13
            GO TO 1030
         ELSEIF(ITYPE.EQ.23) THEN
            JTYPE= 13
            GO TO 1040
         ELSEIF(ITYPE.EQ.59) THEN
            JTYPE= 13
            GO TO 1060
         ENDIF
C
      ELSEIF(LTYPE.EQ.5) THEN
         IF(ITYPE.EQ.18) THEN
            JTYPE= 14
            GO TO 1010
         ELSEIF(ITYPE.EQ.36) THEN
            JTYPE= 14
            GO TO 1020
         ELSEIF(ITYPE.EQ.74) THEN
            JTYPE= 14
            GO TO 1040
         ELSEIF(ITYPE.EQ.76) THEN
            JTYPE= 14
            GO TO 1070
         ELSEIF(ITYPE.EQ.24) THEN
            JTYPE= 15
            GO TO 1010
         ELSEIF(ITYPE.EQ.60) THEN
            JTYPE= 15
            GO TO 1020
         ELSEIF(ITYPE.EQ.26) THEN
            JTYPE= 15
            GO TO 1030
         ELSEIF(ITYPE.EQ.48) THEN
            JTYPE= 15
            GO TO 1040
         ELSEIF(ITYPE.EQ.62) THEN
            JTYPE= 15
            GO TO 1050
         ELSEIF(ITYPE.EQ.66) THEN
            JTYPE= 15
            GO TO 1060
         ELSEIF(ITYPE.EQ.52) THEN
            JTYPE= 15
            GO TO 1070
         ELSEIF(ITYPE.EQ.70) THEN
            JTYPE= 15
            GO TO 1080
         ELSEIF(ITYPE.EQ.42) THEN
            JTYPE= 16
            GO TO 1010
         ELSEIF(ITYPE.EQ.44) THEN
            JTYPE= 16
            GO TO 1030
         ELSEIF(ITYPE.EQ.50) THEN
            JTYPE= 16
            GO TO 1040
         ELSEIF(ITYPE.EQ.68) THEN
            JTYPE= 16
            GO TO 1060
         ENDIF
C
      ELSEIF(LTYPE.EQ.6) THEN
         IF(ITYPE.EQ.27) THEN
            JTYPE= 17
            GO TO 1010
         ELSEIF(ITYPE.EQ.63) THEN
            JTYPE= 17
            GO TO 1020
         ELSEIF(ITYPE.EQ.75) THEN
            JTYPE= 17
            GO TO 1040
         ELSEIF(ITYPE.EQ.79) THEN
            JTYPE= 17
            GO TO 1070
         ELSEIF(ITYPE.EQ.45) THEN
            JTYPE= 18
            GO TO 1010
         ELSEIF(ITYPE.EQ.77) THEN
            JTYPE= 18
            GO TO 1040
         ELSEIF(ITYPE.EQ.51) THEN
            JTYPE= 19
            GO TO 1010
         ELSEIF(ITYPE.EQ.69) THEN
            JTYPE= 19
            GO TO 1020
         ELSEIF(ITYPE.EQ.53) THEN
            JTYPE= 19
            GO TO 1030
         ELSEIF(ITYPE.EQ.71) THEN
            JTYPE= 19
            GO TO 1050
         ENDIF
C
      ELSEIF(LTYPE.EQ.7) THEN
         JTYPE= 20
         IF(ITYPE.EQ.54) THEN
            GO TO 1010
         ELSEIF(ITYPE.EQ.72) THEN
            GO TO 1020
         ELSEIF(ITYPE.EQ.78) THEN
            GO TO 1040
         ELSEIF(ITYPE.EQ.80) THEN
            GO TO 1070
         ENDIF
C
      ELSEIF(LTYPE.EQ.8) THEN
         JTYPE= 21
         GO TO 1010
      ENDIF

1010  CONTINUE 
         CALL GENR03_QMEFP_INTYPE_1(MODE,JM,JNAT,
     *        LAT,LBT,LCT,LDT,P12,P34,R34)
      GO TO 1090

1020  CONTINUE 
         CALL GENR03_QMEFP_INTYPE_2(MODE,JM,JNAT,
     *        LAT,LBT,LCT,LDT,P12,P34,R34)
      GO TO 1090

1030  CONTINUE 
         CALL GENR03_QMEFP_INTYPE_3(MODE,JM,JNAT,
     *        LAT,LBT,LCT,LDT,P12,P34,R34)
      GO TO 1090

1040  CONTINUE 
         CALL GENR03_QMEFP_INTYPE_4(MODE,JM,JNAT,
     *        LAT,LBT,LCT,LDT,P12,P34,R34)
      GO TO 1090

1050  CONTINUE 
         CALL GENR03_QMEFP_INTYPE_5(MODE,JM,JNAT,
     *        LAT,LBT,LCT,LDT,P12,P34,R34)
      GO TO 1090

1060  CONTINUE 
         CALL GENR03_QMEFP_INTYPE_6(MODE,JM,JNAT,
     *        LAT,LBT,LCT,LDT,P12,P34,R34)
      GO TO 1090

1070  CONTINUE 
         CALL GENR03_QMEFP_INTYPE_7(MODE,JM,JNAT,
     *        LAT,LBT,LCT,LDT,P12,P34,R34)
      GO TO 1090

1080  CONTINUE 
         CALL GENR03_QMEFP_INTYPE_8(MODE,JM,JNAT,
     *        LAT,LBT,LCT,LDT,P12,P34,R34)
      GO TO 1090
1090  CONTINUE
C
C
C EMPTY INTEGRAL SUMMATION STORAGE
C
      IKL= 0
      IF(JTYPE.EQ. 7) THEN
         CALL INTK07(IKL)
      ELSEIF(JTYPE.EQ. 8) THEN
         CALL INTK08(IKL)
      ELSEIF(JTYPE.EQ. 9) THEN
         CALL INTK09(IKL)
      ELSEIF(JTYPE.EQ.10) THEN
         CALL INTK10(IKL)
      ELSEIF(JTYPE.EQ.11) THEN
         CALL INTK11(IKL)
      ELSEIF(JTYPE.EQ.12) THEN
         CALL INTK12(IKL)
      ELSEIF(JTYPE.EQ.13) THEN
         CALL INTK13(IKL)
      ELSEIF(JTYPE.EQ.14) THEN
         CALL INTK14(IKL)
      ELSEIF(JTYPE.EQ.15) THEN
         CALL INTK15(IKL)
      ELSEIF(JTYPE.EQ.16) THEN
         CALL INTK16(IKL)
      ELSEIF(JTYPE.EQ.17) THEN
         CALL INTK17(IKL)
      ELSEIF(JTYPE.EQ.18) THEN
         CALL INTK18(IKL)
      ELSEIF(JTYPE.EQ.19) THEN
         CALL INTK19(IKL)
      ELSEIF(JTYPE.EQ.20) THEN
         CALL INTK20(IKL)
      ELSEIF(JTYPE.EQ.21) THEN
         CALL INTK21(IKL)
      ENDIF
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
            TY01(JI)= TY02(JI)-RAB
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
            IF(LB.EQ.2) THEN
               D02D(JI)= E12*CSA(I)*CDB(J)
               IF(LA.EQ.1) THEN
                  D12D(JI)= E12*CPA(I)*CDB(J)
                  IF(D12D(JI).NE.ZER) D02D(JI)=D02D(JI)/D12D(JI)
               ELSEIF(LA.EQ.2) THEN
                  D22D(JI)= E12*CDA(I)*CDB(J)
               ENDIF
            ELSEIF(LB.EQ.1) THEN
               D00P(JI)= E12*CSA(I)*CSB(J)
               D01P(JI)= E12*CSA(I)*CPB(J)
               IF(LA.EQ.0) THEN
                  IF(D01P(JI).NE.ZER) D00P(JI)= D00P(JI)/D01P(JI)
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
            ELSE
               D00P(JI)= E12*CSA(I)*CSB(J)
            ENDIF
  160       CONTINUE
  170 JI=JI+1
C
C BEGIN Q LOOP
C
      IKL= 0
      DO 190 K=1,NGC
         X03= EXC(K)
         DO 180 L=1,NGD
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
            SQ(1)= E34*CSC(K)*CDD(L)
            SQ(2)= E34*CPC(K)*CDD(L)
            SQ(3)= E34*CDC(K)*CDD(L)
C
C USE SPECIAL FAST ROUTINE FOR INNER LOOPS FOR 0000 ... 1111
C
CJMS  ZEROING OF THE FQDx (x=0..8) ARRAYS OF LABELLED COMMON /FQ08/
CJMS  TAKES PLACE IN THE INTJxx (xx=07..21) ROUTINES
C
CJMS  ZEROING OF THE  RDx (x=0..8) ARRAYS OF LABELLED COMMON /KI4 /
CJMS  TOOK PLACE IN THE CALL (FOR IKL=0) OF THE INTKxx (xx=07..20)
CJMS  ROUTINES
C
            IKL= IKL+1
            CALL SPDGEN(JTYPE,IKL)
C
  180    CONTINUE
  190 CONTINUE
C
      QX= RCD*SING
      QZ= RCD*COSG
      IF(JTYPE.EQ. 7) THEN
         CALL MCDV07(GHONDO,QX,QZ)
      ELSEIF(JTYPE.EQ. 8) THEN
         CALL MCDV08(GHONDO,QX,QZ)
      ELSEIF(JTYPE.EQ. 9) THEN
         CALL MCDV09(GHONDO,QX,QZ)
      ELSEIF(JTYPE.EQ.10) THEN
         CALL MCDV10(GHONDO,QX,QZ)
      ELSEIF(JTYPE.EQ.11) THEN
         CALL MCDV11(GHONDO,QX,QZ)
      ELSEIF(JTYPE.EQ.12) THEN
         CALL MCDV12(GHONDO,QX,QZ)
      ELSEIF(JTYPE.EQ.13) THEN
         CALL MCDV13(GHONDO,QX,QZ)
      ELSEIF(JTYPE.EQ.14) THEN
         CALL MCDV14(GHONDO,QX,QZ)
      ELSEIF(JTYPE.EQ.15) THEN
         CALL MCDV15(GHONDO,QX,QZ)
      ELSEIF(JTYPE.EQ.16) THEN
         CALL MCDV16(GHONDO,QX,QZ)
      ELSEIF(JTYPE.EQ.17) THEN
         CALL MCDV17(GHONDO,QX,QZ)
      ELSEIF(JTYPE.EQ.18) THEN
         CALL MCDV18(GHONDO,QX,QZ)
      ELSEIF(JTYPE.EQ.19) THEN
         CALL MCDV19(GHONDO,QX,QZ)
      ELSEIF(JTYPE.EQ.20) THEN
         CALL MCDV20(GHONDO,QX,QZ)
      ELSEIF(JTYPE.EQ.21) THEN
         CALL MCDV21(GHONDO,QX,QZ)
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
      CALL R30S1D(JTYPE,GHONDO,P)
C
      RETURN
      END



C*MODULE EFTEI_GENR03   *DECK SETSPD_SHELLEFP
C>    @brief  prepare EFP shells
C>    @author  Peng Xu and Tosaporn Sattasathuchana 
C>             Jan 2020
C>    @detail prepare EFP shells for GENR03_QMEFP subroutine
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
C>    @param  MINI: starting counter of shell I
C>    @param  MAXI: ending counter of shell I
C>    @param  LOCI: location of the basis function in shell I 
      SUBROUTINE SETSPD_SHELLEFP(JM,JNAT,ISHELL,LAT,INEW,
     *  LA,NGA,CMAXA,EXA,CSA,CPA,CDA,CRDI,MINI,MAXI,LOCI)
C
      use mx_limits, only: mxgtot,mxsh,mxgsh,mxg2,MXSH
      USE comm_EFPBAS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION CRDI(MXSH,3)
      DIMENSION EXA(MXGSH),CSA(MXGSH),CPA(MXGSH),CMAXA(MXGSH)
      DIMENSION CDA(MXGSH)

      INEW= ISHELL
      LA= LAT
      NGA  = KNGEF(INEW,JM)
      I    = KSTREF(INEW,JM)-1
      ICC  = KATMEF(INEW,JM)
      MINI = KMINEF(INEW,JM)
      MAXI = KMAXEF(INEW,JM)
      LOCI = KLOCEF(INEW,JM)-MINI

      DO NI=1,NGA
         N=I+NI
         CMAXA(NI)= CMAXEF(N,JM)
         EXA(NI)= EXEF(N,JM)
         CSA(NI)= CSEF(N,JM)
         CPA(NI)= CPEF(N,JM)
         CDA(NI)= CDEF(N,JM)
      ENDDO

      DO N=1,3
         CRDI(INEW,N) = PRCORD(N,ICC+JNAT)
      ENDDO
      RETURN
      END
C

C*MODULE EFTEI_GENR03   *DECK SETSPD_SHELLAI
C>    @brief  prepare QM shells
C>    @author  Peng Xu and Tosaporn Sattasathuchana 
C>             Jan 2020
C>    @detail prepare QM shells for GENR03_QMEFP subroutine
C>
C>    @param  JSHELL: index of the shell
C>    @param  LBT: type of basis function 
C>    @param  JNEW: updated shell index
C>    @param  LB: updated basis type
C>    @param  NGB: number of gaussian functions in a shell
C>    @param  CMAX: maximum coefficient
C>    @param  CMAXB: updated maximum coefficent 
C>    @param  EXB: exponenent of basis function
C>    @param  CSB: coefficient of s basis function
C>    @param  CPB: coefficient of p basis function
C>    @param  CRDJ: coordinates of atoms associated with the shell
C>    @param  MINJ: starting counter of shell J 
C>    @param  MAXJ: ending counter of shell J
C>    @param  LOCJ: location of the basis function in shell J 
      SUBROUTINE SETSPD_SHELLAI(JSHELL,LBT,CMAX,
     *  JNEW,LB,NGB,CMAXB,EXB,CSB,CPB,CDB,CRDJ,MINJ,MAXJ,LOCJ)
      use mx_limits, only: mxgtot,mxsh,mxgsh,mxg2,mxatm
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION CMAX(MXGTOT)
      DIMENSION CRDJ(MXSH,3)
      DIMENSION EXB(MXGSH),CSB(MXGSH),CPB(MXGSH),CMAXB(MXGSH)
      DIMENSION CDB(MXGSH)

      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /NSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),
     2                CF(MXGTOT),CG(MXGTOT),CH(MXGTOT),CI(MXGTOT),
     3                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),
     4                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
C
      JNEW= JSHELL
      LB= LBT
      NGB  = KNG(JNEW)
      J    = KSTART(JNEW)-1
      MINJ = KMIN(JNEW)
      MAXJ = KMAX(JNEW)
      LOCJ = KLOC(JNEW)-MINJ

      DO NJ=1,NGB
         N=J+NJ
         CMAXB(NJ)= CMAX(N)
         EXB(NJ)= EX(N)
         CSB(NJ)= CS(N)
         CPB(NJ)= CP(N)
         CDB(NJ)= CD(N)
      ENDDO

      JATM=KATOM(JNEW)
      DO N=1,3
         CRDJ(JNEW,N)= C(N,JATM)
      ENDDO

      RETURN
      END

C*MODULE EFTEI_GENR03   *DECK GENR03_QMEFP_INTYPE_1 
C>    @brief  rotated axis (spd) for case 1 (IJ|KL) 
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

      SUBROUTINE GENR03_QMEFP_INTYPE_1(MODE,JM,JNAT,
     *        LAT,LBT,LCT,LDT,P12,P34,R34)

      use mx_limits, only: mxgtot,mxsh,mxgsh,mxg2
      USE comm_EFPBAS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      COMMON /FLIPS / IB(4,3)
      COMMON /GEOMPQ/ R12,RAB,X34,X43,AQZ,QPR,QPS,
     2                TX12(MXG2),TX21(MXG2),TY01(MXG2),TY02(MXG2),
     3                D00P(MXG2),D01P(MXG2),D10P(MXG2),D11P(MXG2),
     4                NGANGB
      COMMON /MAXC  / CMAX(MXGTOT),CMAXA(MXGSH),CMAXB(MXGSH),
     2                CMAXC(MXGSH),CMAXD(MXGSH),ISMLP(MXG2),ISMLQ
      COMMON /NSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),
     2                CF(MXGTOT),CG(MXGTOT),CH(MXGTOT),CI(MXGTOT),
     3                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),
     4                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
      COMMON /SHLG70/ ISHELL,JSHELL,KSHELL,LSHELL,INEW,JNEW,KNEW,LNEW
      COMMON /SHLLFO/ NGA,LA,EXA(MXGSH),CSA(MXGSH),CPA(MXGSH),
     2                NGB,LB,EXB(MXGSH),CSB(MXGSH),CPB(MXGSH),
     3                NGC,LC,EXC(MXGSH),CSC(MXGSH),CPC(MXGSH),
     4                NGD,LD,EXD(MXGSH),CSD(MXGSH),CPD(MXGSH)
      COMMON /SHLNOS/ QQ4,LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,
     2                MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL,
     3                NIJ,IJ,KL,IJKL
      COMMON /SHLSPD/ CDA(MXGSH),CDB(MXGSH),CDC(MXGSH),CDD(MXGSH),
     2                D02D(MXG2),D12D(MXG2),D22D(MXG2)

      DIMENSION P12(3,3),P34(3,3)
      DIMENSION CRDI(MXSH,3), CRDJ(MXSH,3), CRDK(MXSH,3), CRDL(MXSH,3)
C
      IEXCH = 1
      IB(1,IEXCH)= 1
      IB(2,IEXCH)= 2
      IB(3,IEXCH)= 3
      IB(4,IEXCH)= 4

      CALL SETSPD_SHELLEFP(JM,JNAT,ISHELL,LAT,
     *  INEW,LA,NGA,CMAXA,EXA,CSA,CPA,CDA,CRDI,MINI,MAXI,LOCI)
C
      CALL SETSPD_SHELLAI(JSHELL,LBT,CMAX,
     *  JNEW,LB,NGB,CMAXB,EXB,CSB,CPB,CDB,CRDJ,MINJ,MAXJ,LOCJ)
C
      IF(MODE.EQ.2) THEN
       CALL SETSPD_SHELLEFP(JM,JNAT,KSHELL,LCT,
     *  KNEW,LC,NGC,CMAXC,EXC,CSC,CPC,CDC,CRDK,MINK,MAXK,LOCK)
C
      ELSE
       CALL SETSPD_SHELLAI(KSHELL,LCT,CMAX,
     *  KNEW,LC,NGC,CMAXC,EXC,CSC,CPC,CDC,CRDK,MINK,MAXK,LOCK)
C
      ENDIF

      CALL SETSPD_SHELLAI(LSHELL,LDT,CMAX,
     *  LNEW,LD,NGD,CMAXD,EXD,CSD,CPD,CDD,CRDL,MINL,MAXL,LOCL)

      NGANGB = NGA * NGB

      CALL SETUP_COORD(CRDI,CRDJ,CRDK,CRDL,INEW,JNEW,KNEW,LNEW,
     &                       P12,P34,R12,R34)

      RETURN
      END

C*MODULE EFTEI_GENR03   *DECK GENR03_QMEFP_INTYPE_2 
C>    @brief  rotated axis (spd) for case 2 (JI|KL) 
C>    @author  Peng Xu and Tosaporn Sattasathuchana 
C>             Jan 2020
C>    @detail Case 2 integral includes types:
C>            1001,1011 (IJ SWITCHED)
C>    @see    parameters same as GENR03_QMEFP_INTYPE_1

      SUBROUTINE GENR03_QMEFP_INTYPE_2(MODE,JM,JNAT,
     *        LAT,LBT,LCT,LDT,P12,P34,R34)

      use mx_limits, only: mxgtot,mxsh,mxgsh,mxg2
      USE comm_EFPBAS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      COMMON /FLIPS / IB(4,3)
      COMMON /GEOMPQ/ R12,RAB,X34,X43,AQZ,QPR,QPS,
     2                TX12(MXG2),TX21(MXG2),TY01(MXG2),TY02(MXG2),
     3                D00P(MXG2),D01P(MXG2),D10P(MXG2),D11P(MXG2),
     4                NGANGB
      COMMON /MAXC  / CMAX(MXGTOT),CMAXA(MXGSH),CMAXB(MXGSH),
     2                CMAXC(MXGSH),CMAXD(MXGSH),ISMLP(MXG2),ISMLQ
      COMMON /NSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),
     2                CF(MXGTOT),CG(MXGTOT),CH(MXGTOT),CI(MXGTOT),
     3                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),
     4                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
      COMMON /SHLG70/ ISHELL,JSHELL,KSHELL,LSHELL,INEW,JNEW,KNEW,LNEW
      COMMON /SHLLFO/ NGA,LA,EXA(MXGSH),CSA(MXGSH),CPA(MXGSH),
     2                NGB,LB,EXB(MXGSH),CSB(MXGSH),CPB(MXGSH),
     3                NGC,LC,EXC(MXGSH),CSC(MXGSH),CPC(MXGSH),
     4                NGD,LD,EXD(MXGSH),CSD(MXGSH),CPD(MXGSH)
      COMMON /SHLNOS/ QQ4,LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,
     2                MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL,
     3                NIJ,IJ,KL,IJKL
      COMMON /SHLSPD/ CDA(MXGSH),CDB(MXGSH),CDC(MXGSH),CDD(MXGSH),
     2                D02D(MXG2),D12D(MXG2),D22D(MXG2)

      DIMENSION P12(3,3),P34(3,3)
      DIMENSION CRDI(MXSH,3), CRDJ(MXSH,3), CRDK(MXSH,3), CRDL(MXSH,3)

      IEXCH=1
      IB(1,IEXCH)= 2
      IB(2,IEXCH)= 1
      IB(3,IEXCH)= 3
      IB(4,IEXCH)= 4

      CALL SETSPD_SHELLAI(JSHELL,LBT,CMAX,
     *  INEW,LA,NGA,CMAXA,EXA,CSA,CPA,CDA,CRDI,MINI,MAXI,LOCI)
C
      CALL SETSPD_SHELLEFP(JM,JNAT,ISHELL,LAT,
     *  JNEW,LB,NGB,CMAXB,EXB,CSB,CPB,CDB,CRDJ,MINJ,MAXJ,LOCJ)
C
      IF(MODE.EQ.2) THEN
       CALL SETSPD_SHELLEFP(JM,JNAT,KSHELL,LCT,
     *  KNEW,LC,NGC,CMAXC,EXC,CSC,CPC,CDC,CRDK,MINK,MAXK,LOCK)
C
      ELSE
       CALL SETSPD_SHELLAI(KSHELL,LCT,CMAX,
     *  KNEW,LC,NGC,CMAXC,EXC,CSC,CPC,CDC,CRDK,MINK,MAXK,LOCK)
C
      ENDIF

      CALL SETSPD_SHELLAI(LSHELL,LDT,CMAX,
     *  LNEW,LD,NGD,CMAXD,EXD,CSD,CPD,CDD,CRDL,MINL,MAXL,LOCL)

      NGANGB = NGA * NGB

      CALL SETUP_COORD(CRDI,CRDJ,CRDK,CRDL,INEW,JNEW,KNEW,LNEW,
     &                       P12,P34,R12,R34)

      RETURN
      END

C*MODULE EFTEI_GENR03   *DECK GENR03_QMEFP_INTYPE_3 
C>    @brief  rotated axis (spd) for case 3 (IJ|LK) 
C>    @author  Peng Xu and Tosaporn Sattasathuchana 
C>             Jan 2020
C>    @detail Case 3 integral includes types:
C>            0010,0110 (KL SWITCHED)
C>    @see    parameters same as GENR03_QMEFP_INTYPE_1

      SUBROUTINE GENR03_QMEFP_INTYPE_3(MODE,JM,JNAT,
     *        LAT,LBT,LCT,LDT,P12,P34,R34)

      use mx_limits, only: mxgtot,mxsh,mxgsh,mxg2
      USE comm_EFPBAS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      COMMON /FLIPS / IB(4,3)
      COMMON /GEOMPQ/ R12,RAB,X34,X43,AQZ,QPR,QPS,
     2                TX12(MXG2),TX21(MXG2),TY01(MXG2),TY02(MXG2),
     3                D00P(MXG2),D01P(MXG2),D10P(MXG2),D11P(MXG2),
     4                NGANGB
      COMMON /MAXC  / CMAX(MXGTOT),CMAXA(MXGSH),CMAXB(MXGSH),
     2                CMAXC(MXGSH),CMAXD(MXGSH),ISMLP(MXG2),ISMLQ
      COMMON /NSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),
     2                CF(MXGTOT),CG(MXGTOT),CH(MXGTOT),CI(MXGTOT),
     3                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),
     4                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
      COMMON /SHLG70/ ISHELL,JSHELL,KSHELL,LSHELL,INEW,JNEW,KNEW,LNEW
      COMMON /SHLLFO/ NGA,LA,EXA(MXGSH),CSA(MXGSH),CPA(MXGSH),
     2                NGB,LB,EXB(MXGSH),CSB(MXGSH),CPB(MXGSH),
     3                NGC,LC,EXC(MXGSH),CSC(MXGSH),CPC(MXGSH),
     4                NGD,LD,EXD(MXGSH),CSD(MXGSH),CPD(MXGSH)
      COMMON /SHLNOS/ QQ4,LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,
     2                MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL,
     3                NIJ,IJ,KL,IJKL
      COMMON /SHLSPD/ CDA(MXGSH),CDB(MXGSH),CDC(MXGSH),CDD(MXGSH),
     2                D02D(MXG2),D12D(MXG2),D22D(MXG2)

      DIMENSION P12(3,3),P34(3,3)
      DIMENSION CRDI(MXSH,3), CRDJ(MXSH,3), CRDK(MXSH,3), CRDL(MXSH,3)

      IEXCH = 1
      IB(1,IEXCH)= 1
      IB(2,IEXCH)= 2
      IB(3,IEXCH)= 4
      IB(4,IEXCH)= 3

      CALL SETSPD_SHELLEFP(JM,JNAT,ISHELL,LAT,
     *  INEW,LA,NGA,CMAXA,EXA,CSA,CPA,CDA,CRDI,MINI,MAXI,LOCI)
C
      CALL SETSPD_SHELLAI(JSHELL,LBT,CMAX,
     *  JNEW,LB,NGB,CMAXB,EXB,CSB,CPB,CDB,CRDJ,MINJ,MAXJ,LOCJ)
C
      CALL SETSPD_SHELLAI(LSHELL,LDT,CMAX,
     *  KNEW,LC,NGC,CMAXC,EXC,CSC,CPC,CDC,CRDK,MINK,MAXK,LOCK)
C
      IF(MODE.EQ.2) THEN
       CALL SETSPD_SHELLEFP(JM,JNAT,KSHELL,LCT,
     *  LNEW,LD,NGD,CMAXD,EXD,CSD,CPD,CDD,CRDL,MINL,MAXL,LOCL)
C
      ELSE
       CALL SETSPD_SHELLAI(KSHELL,LCT,CMAX,
     *  LNEW,LD,NGD,CMAXD,EXD,CSD,CPD,CDD,CRDL,MINL,MAXL,LOCL)
C
      ENDIF

      NGANGB = NGA * NGB

      CALL SETUP_COORD(CRDI,CRDJ,CRDK,CRDL,INEW,JNEW,KNEW,LNEW,
     &                       P12,P34,R12,R34)


      RETURN
      END

C*MODULE EFTEI_GENR03   *DECK GENR03_QMEFP_INTYPE_4 
C>    @brief  rotated axis (spd) for case 4 (KL|IJ)
C>    @author  Peng Xu and Tosaporn Sattasathuchana 
C>             Jan 2020
C>    @detail Case 4 integral includes types:
C>            0100,1100,1101 (pairs IJ and KL SWITCHED)
C>    @see    parameters same as GENR03_QMEFP_INTYPE_1

      SUBROUTINE GENR03_QMEFP_INTYPE_4(MODE,JM,JNAT,
     *        LAT,LBT,LCT,LDT,P12,P34,R34)

      use mx_limits, only: mxgtot,mxsh,mxgsh,mxg2
      USE comm_EFPBAS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      COMMON /FLIPS / IB(4,3)
      COMMON /GEOMPQ/ R12,RAB,X34,X43,AQZ,QPR,QPS,
     2                TX12(MXG2),TX21(MXG2),TY01(MXG2),TY02(MXG2),
     3                D00P(MXG2),D01P(MXG2),D10P(MXG2),D11P(MXG2),
     4                NGANGB
      COMMON /MAXC  / CMAX(MXGTOT),CMAXA(MXGSH),CMAXB(MXGSH),
     2                CMAXC(MXGSH),CMAXD(MXGSH),ISMLP(MXG2),ISMLQ
      COMMON /NSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),
     2                CF(MXGTOT),CG(MXGTOT),CH(MXGTOT),CI(MXGTOT),
     3                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),
     4                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
      COMMON /SHLG70/ ISHELL,JSHELL,KSHELL,LSHELL,INEW,JNEW,KNEW,LNEW
      COMMON /SHLLFO/ NGA,LA,EXA(MXGSH),CSA(MXGSH),CPA(MXGSH),
     2                NGB,LB,EXB(MXGSH),CSB(MXGSH),CPB(MXGSH),
     3                NGC,LC,EXC(MXGSH),CSC(MXGSH),CPC(MXGSH),
     4                NGD,LD,EXD(MXGSH),CSD(MXGSH),CPD(MXGSH)
      COMMON /SHLNOS/ QQ4,LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,
     2                MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL,
     3                NIJ,IJ,KL,IJKL
      COMMON /SHLSPD/ CDA(MXGSH),CDB(MXGSH),CDC(MXGSH),CDD(MXGSH),
     2                D02D(MXG2),D12D(MXG2),D22D(MXG2)

      DIMENSION P12(3,3),P34(3,3)
      DIMENSION CRDI(MXSH,3), CRDJ(MXSH,3), CRDK(MXSH,3), CRDL(MXSH,3)

      IEXCH = 1
      IB(1,IEXCH)= 3
      IB(2,IEXCH)= 4
      IB(3,IEXCH)= 1
      IB(4,IEXCH)= 2

      IF(MODE.EQ.2) THEN
       CALL SETSPD_SHELLEFP(JM,JNAT,KSHELL,LCT,
     *  INEW,LA,NGA,CMAXA,EXA,CSA,CPA,CDA,CRDI,MINI,MAXI,LOCI)
C
      ELSE
       CALL SETSPD_SHELLAI(KSHELL,LCT,CMAX,
     *  INEW,LA,NGA,CMAXA,EXA,CSA,CPA,CDA,CRDI,MINI,MAXI,LOCI)
C
      ENDIF

      CALL SETSPD_SHELLAI(LSHELL,LDT,CMAX,
     *  JNEW,LB,NGB,CMAXB,EXB,CSB,CPB,CDB,CRDJ,MINJ,MAXJ,LOCJ)
C
      CALL SETSPD_SHELLEFP(JM,JNAT,ISHELL,LAT,
     *  KNEW,LC,NGC,CMAXC,EXC,CSC,CPC,CDC,CRDK,MINK,MAXK,LOCK)
C
      CALL SETSPD_SHELLAI(JSHELL,LBT,CMAX,
     *  LNEW,LD,NGD,CMAXD,EXD,CSD,CPD,CDD,CRDL,MINL,MAXL,LOCL)

      NGANGB = NGA * NGB

      CALL SETUP_COORD(CRDI,CRDJ,CRDK,CRDL,INEW,JNEW,KNEW,LNEW,
     &                       P12,P34,R12,R34)

      RETURN
      END

C*MODULE EFTEI_GENR03   *DECK GENR03_QMEFP_INTYPE_5 
C>    @brief  rotated axis (spd) for case 5 (JI|LK) 
C>    @author  Peng Xu and Tosaporn Sattasathuchana 
C>             Jan 2020
C>    @detail Case 5 integral includes types:
C>            1010 (IJ SWITCHED and KL SWITCHED)
C>    @see    parameters same as GENR03_QMEFP_INTYPE_1

      SUBROUTINE GENR03_QMEFP_INTYPE_5(MODE,JM,JNAT,
     *        LAT,LBT,LCT,LDT,P12,P34,R34)

      use mx_limits, only: mxgtot,mxsh,mxgsh,mxg2
      USE comm_EFPBAS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      COMMON /FLIPS / IB(4,3)
      COMMON /GEOMPQ/ R12,RAB,X34,X43,AQZ,QPR,QPS,
     2                TX12(MXG2),TX21(MXG2),TY01(MXG2),TY02(MXG2),
     3                D00P(MXG2),D01P(MXG2),D10P(MXG2),D11P(MXG2),
     4                NGANGB
      COMMON /MAXC  / CMAX(MXGTOT),CMAXA(MXGSH),CMAXB(MXGSH),
     2                CMAXC(MXGSH),CMAXD(MXGSH),ISMLP(MXG2),ISMLQ
      COMMON /NSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),
     2                CF(MXGTOT),CG(MXGTOT),CH(MXGTOT),CI(MXGTOT),
     3                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),
     4                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
      COMMON /SHLG70/ ISHELL,JSHELL,KSHELL,LSHELL,INEW,JNEW,KNEW,LNEW
      COMMON /SHLLFO/ NGA,LA,EXA(MXGSH),CSA(MXGSH),CPA(MXGSH),
     2                NGB,LB,EXB(MXGSH),CSB(MXGSH),CPB(MXGSH),
     3                NGC,LC,EXC(MXGSH),CSC(MXGSH),CPC(MXGSH),
     4                NGD,LD,EXD(MXGSH),CSD(MXGSH),CPD(MXGSH)
      COMMON /SHLNOS/ QQ4,LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,
     2                MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL,
     3                NIJ,IJ,KL,IJKL
      COMMON /SHLSPD/ CDA(MXGSH),CDB(MXGSH),CDC(MXGSH),CDD(MXGSH),
     2                D02D(MXG2),D12D(MXG2),D22D(MXG2)

      DIMENSION P12(3,3),P34(3,3)
      DIMENSION CRDI(MXSH,3), CRDJ(MXSH,3), CRDK(MXSH,3), CRDL(MXSH,3)

      IEXCH = 1
      IB(1,IEXCH)= 2
      IB(2,IEXCH)= 1
      IB(3,IEXCH)= 4
      IB(4,IEXCH)= 3

      CALL SETSPD_SHELLAI(JSHELL,LBT,CMAX,
     *  INEW,LA,NGA,CMAXA,EXA,CSA,CPA,CDA,CRDI,MINI,MAXI,LOCI)
C
      CALL SETSPD_SHELLEFP(JM,JNAT,ISHELL,LAT,
     *  JNEW,LB,NGB,CMAXB,EXB,CSB,CPB,CDB,CRDJ,MINJ,MAXJ,LOCJ)
C
      CALL SETSPD_SHELLAI(LSHELL,LDT,CMAX,
     *  KNEW,LC,NGC,CMAXC,EXC,CSC,CPC,CDC,CRDK,MINK,MAXK,LOCK)
C
      IF(MODE.EQ.2) THEN
       CALL SETSPD_SHELLEFP(JM,JNAT,KSHELL,LCT,
     *  LNEW,LD,NGD,CMAXD,EXD,CSD,CPD,CDD,CRDL,MINL,MAXL,LOCL)

      ELSE
       CALL SETSPD_SHELLAI(KSHELL,LCT,CMAX,
     *  LNEW,LD,NGD,CMAXD,EXD,CSD,CPD,CDD,CRDL,MINL,MAXL,LOCL)
C
      ENDIF

      NGANGB = NGA * NGB

      CALL SETUP_COORD(CRDI,CRDJ,CRDK,CRDL,INEW,JNEW,KNEW,LNEW,
     &                       P12,P34,R12,R34)

      RETURN
      END

C*MODULE EFTEI_GENR03   *DECK GENR03_QMEFP_INTYPE_6 
C>    @brief  rotated axis (spd) for case 6 (KL|JI) 
C>    @author  Peng Xu and Tosaporn Sattasathuchana 
C>             Jan 2020
C>    @detail Case 6 integral includes types:
C>            1000 (PAIRS IJ and KL SWITCHED,followed by KL SWITCH)
C>    @see    parameters same as GENR03_QMEFP_INTYPE_1

      SUBROUTINE GENR03_QMEFP_INTYPE_6(MODE,JM,JNAT,
     *        LAT,LBT,LCT,LDT,P12,P34,R34)

      use mx_limits, only: mxgtot,mxsh,mxgsh,mxg2
      USE comm_EFPBAS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      COMMON /FLIPS / IB(4,3)
      COMMON /GEOMPQ/ R12,RAB,X34,X43,AQZ,QPR,QPS,
     2                TX12(MXG2),TX21(MXG2),TY01(MXG2),TY02(MXG2),
     3                D00P(MXG2),D01P(MXG2),D10P(MXG2),D11P(MXG2),
     4                NGANGB
      COMMON /MAXC  / CMAX(MXGTOT),CMAXA(MXGSH),CMAXB(MXGSH),
     2                CMAXC(MXGSH),CMAXD(MXGSH),ISMLP(MXG2),ISMLQ
      COMMON /NSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),
     2                CF(MXGTOT),CG(MXGTOT),CH(MXGTOT),CI(MXGTOT),
     3                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),
     4                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
      COMMON /SHLG70/ ISHELL,JSHELL,KSHELL,LSHELL,INEW,JNEW,KNEW,LNEW
      COMMON /SHLLFO/ NGA,LA,EXA(MXGSH),CSA(MXGSH),CPA(MXGSH),
     2                NGB,LB,EXB(MXGSH),CSB(MXGSH),CPB(MXGSH),
     3                NGC,LC,EXC(MXGSH),CSC(MXGSH),CPC(MXGSH),
     4                NGD,LD,EXD(MXGSH),CSD(MXGSH),CPD(MXGSH)
      COMMON /SHLNOS/ QQ4,LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,
     2                MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL,
     3                NIJ,IJ,KL,IJKL
      COMMON /SHLSPD/ CDA(MXGSH),CDB(MXGSH),CDC(MXGSH),CDD(MXGSH),
     2                D02D(MXG2),D12D(MXG2),D22D(MXG2)

      DIMENSION P12(3,3),P34(3,3)
      DIMENSION CRDI(MXSH,3), CRDJ(MXSH,3), CRDK(MXSH,3), CRDL(MXSH,3)

      IEXCH = 1
      IB(1,IEXCH)= 4
      IB(2,IEXCH)= 3
      IB(3,IEXCH)= 1
      IB(4,IEXCH)= 2

      IF(MODE.EQ.2) THEN
       CALL SETSPD_SHELLEFP(JM,JNAT,KSHELL,LCT,
     *  INEW,LA,NGA,CMAXA,EXA,CSA,CPA,CDA,CRDI,MINI,MAXI,LOCI)
C
      ELSE
       CALL SETSPD_SHELLAI(KSHELL,LCT,CMAX,
     *  INEW,LA,NGA,CMAXA,EXA,CSA,CPA,CDA,CRDI,MINI,MAXI,LOCI)
C
      ENDIF

      CALL SETSPD_SHELLAI(LSHELL,LDT,CMAX,
     *  JNEW,LB,NGB,CMAXB,EXB,CSB,CPB,CDB,CRDJ,MINJ,MAXJ,LOCJ)
C
      CALL SETSPD_SHELLAI(JSHELL,LBT,CMAX,
     *  KNEW,LC,NGC,CMAXC,EXC,CSC,CPC,CDC,CRDK,MINK,MAXK,LOCK)
C
      CALL SETSPD_SHELLEFP(JM,JNAT,ISHELL,LAT,
     *  LNEW,LD,NGD,CMAXD,EXD,CSD,CPD,CDD,CRDL,MINL,MAXL,LOCL)

      NGANGB = NGA * NGB

      CALL SETUP_COORD(CRDI,CRDJ,CRDK,CRDL,INEW,JNEW,KNEW,LNEW,
     &                       P12,P34,R12,R34)

      RETURN
      END

C*MODULE EFTEI_GENR03   *DECK GENR03_QMEFP_INTYPE_7 
C>    @brief  rotated axis (spd) for case 7 (LK|IJ) 
C>    @author  Peng Xu and Tosaporn Sattasathuchana 
C>             Jan 2020
C>    @detail Case 7 integral includes types:
C>            1110 (PAIRS IJ and KL SWITCHED followed by IJ SWITCH)
C>    @see    parameters same as GENR03_QMEFP_INTYPE_1

      SUBROUTINE GENR03_QMEFP_INTYPE_7(MODE,JM,JNAT,
     *        LAT,LBT,LCT,LDT,P12,P34,R34)

      use mx_limits, only: mxgtot,mxsh,mxgsh,mxg2
      USE comm_EFPBAS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      COMMON /FLIPS / IB(4,3)
      COMMON /GEOMPQ/ R12,RAB,X34,X43,AQZ,QPR,QPS,
     2                TX12(MXG2),TX21(MXG2),TY01(MXG2),TY02(MXG2),
     3                D00P(MXG2),D01P(MXG2),D10P(MXG2),D11P(MXG2),
     4                NGANGB
      COMMON /MAXC  / CMAX(MXGTOT),CMAXA(MXGSH),CMAXB(MXGSH),
     2                CMAXC(MXGSH),CMAXD(MXGSH),ISMLP(MXG2),ISMLQ
      COMMON /NSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),
     2                CF(MXGTOT),CG(MXGTOT),CH(MXGTOT),CI(MXGTOT),
     3                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),
     4                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
      COMMON /SHLG70/ ISHELL,JSHELL,KSHELL,LSHELL,INEW,JNEW,KNEW,LNEW
      COMMON /SHLLFO/ NGA,LA,EXA(MXGSH),CSA(MXGSH),CPA(MXGSH),
     2                NGB,LB,EXB(MXGSH),CSB(MXGSH),CPB(MXGSH),
     3                NGC,LC,EXC(MXGSH),CSC(MXGSH),CPC(MXGSH),
     4                NGD,LD,EXD(MXGSH),CSD(MXGSH),CPD(MXGSH)
      COMMON /SHLNOS/ QQ4,LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,
     2                MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL,
     3                NIJ,IJ,KL,IJKL
      COMMON /SHLSPD/ CDA(MXGSH),CDB(MXGSH),CDC(MXGSH),CDD(MXGSH),
     2                D02D(MXG2),D12D(MXG2),D22D(MXG2)

      DIMENSION P12(3,3),P34(3,3)
      DIMENSION CRDI(MXSH,3), CRDJ(MXSH,3), CRDK(MXSH,3), CRDL(MXSH,3)

      IEXCH = 1
      IB(1,IEXCH)= 3
      IB(2,IEXCH)= 4
      IB(3,IEXCH)= 2
      IB(4,IEXCH)= 1

      CALL SETSPD_SHELLAI(LSHELL,LDT,CMAX,
     *  INEW,LA,NGA,CMAXA,EXA,CSA,CPA,CDA,CRDI,MINI,MAXI,LOCI)
C
      IF(MODE.EQ.2) THEN
       CALL SETSPD_SHELLEFP(JM,JNAT,KSHELL,LCT,
     *  JNEW,LB,NGB,CMAXB,EXB,CSB,CPB,CDB,CRDJ,MINJ,MAXJ,LOCJ)
C
      ELSE
       CALL SETSPD_SHELLAI(KSHELL,LCT,CMAX,
     *  JNEW,LB,NGB,CMAXB,EXB,CSB,CPB,CDB,CRDJ,MINJ,MAXJ,LOCJ)
C
      ENDIF

      CALL SETSPD_SHELLEFP(JM,JNAT,ISHELL,LAT,
     *  KNEW,LC,NGC,CMAXC,EXC,CSC,CPC,CDC,CRDK,MINK,MAXK,LOCK)
C

      CALL SETSPD_SHELLAI(JSHELL,LBT,CMAX,
     *  LNEW,LD,NGD,CMAXD,EXD,CSD,CPD,CDD,CRDL,MINL,MAXL,LOCL)

      NGANGB = NGA * NGB

      CALL SETUP_COORD(CRDI,CRDJ,CRDK,CRDL,INEW,JNEW,KNEW,LNEW,
     &                       P12,P34,R12,R34)

      RETURN
      END


C*MODULE EFTEI_GENR03   *DECK GENR03_QMEFP_INTYPE_8 
C>    @brief  rotated axis (spd) for case 8 (LK|JI) 
C>    @author  Peng Xu and Tosaporn Sattasathuchana 
C>             Jan 2020
C>    @detail Case 8 integral includes types:
C>            PAIRS IJ and KL SWITCHED followed by IJ SWITCH
C>            and KL SWITCH
C>    @see    parameters same as GENR03_QMEFP_INTYPE_1

      SUBROUTINE GENR03_QMEFP_INTYPE_8(MODE,JM,JNAT,
     *        LAT,LBT,LCT,LDT,P12,P34,R34)

      use mx_limits, only: mxgtot,mxsh,mxgsh,mxg2
      USE comm_EFPBAS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      COMMON /FLIPS / IB(4,3)
      COMMON /GEOMPQ/ R12,RAB,X34,X43,AQZ,QPR,QPS,
     2                TX12(MXG2),TX21(MXG2),TY01(MXG2),TY02(MXG2),
     3                D00P(MXG2),D01P(MXG2),D10P(MXG2),D11P(MXG2),
     4                NGANGB
      COMMON /MAXC  / CMAX(MXGTOT),CMAXA(MXGSH),CMAXB(MXGSH),
     2                CMAXC(MXGSH),CMAXD(MXGSH),ISMLP(MXG2),ISMLQ
      COMMON /NSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),
     2                CF(MXGTOT),CG(MXGTOT),CH(MXGTOT),CI(MXGTOT),
     3                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),
     4                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
      COMMON /SHLG70/ ISHELL,JSHELL,KSHELL,LSHELL,INEW,JNEW,KNEW,LNEW
      COMMON /SHLLFO/ NGA,LA,EXA(MXGSH),CSA(MXGSH),CPA(MXGSH),
     2                NGB,LB,EXB(MXGSH),CSB(MXGSH),CPB(MXGSH),
     3                NGC,LC,EXC(MXGSH),CSC(MXGSH),CPC(MXGSH),
     4                NGD,LD,EXD(MXGSH),CSD(MXGSH),CPD(MXGSH)
      COMMON /SHLNOS/ QQ4,LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,
     2                MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL,
     3                NIJ,IJ,KL,IJKL
      COMMON /SHLSPD/ CDA(MXGSH),CDB(MXGSH),CDC(MXGSH),CDD(MXGSH),
     2                D02D(MXG2),D12D(MXG2),D22D(MXG2)

      DIMENSION P12(3,3),P34(3,3)
      DIMENSION CRDI(MXSH,3), CRDJ(MXSH,3), CRDK(MXSH,3), CRDL(MXSH,3)

      IEXCH = 1
      IB(1,IEXCH)= 4
      IB(2,IEXCH)= 3
      IB(3,IEXCH)= 2
      IB(4,IEXCH)= 1

      CALL SETSPD_SHELLAI(LSHELL,LDT,CMAX,
     *  INEW,LA,NGA,CMAXA,EXA,CSA,CPA,CDA,CRDI,MINI,MAXI,LOCI)
C
      IF(MODE.EQ.2) THEN
       CALL SETSPD_SHELLEFP(JM,JNAT,KSHELL,LCT,
     *  JNEW,LB,NGB,CMAXB,EXB,CSB,CPB,CDB,CRDJ,MINJ,MAXJ,LOCJ)
C
      ELSE
       CALL SETSPD_SHELLAI(KSHELL,LCT,CMAX,
     *  JNEW,LB,NGB,CMAXB,EXB,CSB,CPB,CDB,CRDJ,MINJ,MAXJ,LOCJ)
C
      ENDIF
      CALL SETSPD_SHELLAI(JSHELL,LBT,CMAX,
     *  KNEW,LC,NGC,CMAXC,EXC,CSC,CPC,CDC,CRDK,MINK,MAXK,LOCK)

      CALL SETSPD_SHELLEFP(JM,JNAT,ISHELL,LAT,
     *  LNEW,LD,NGD,CMAXD,EXD,CSD,CPD,CDD,CRDL,MINL,MAXL,LOCL)

      NGANGB = NGA * NGB

      CALL SETUP_COORD(CRDI,CRDJ,CRDK,CRDL,INEW,JNEW,KNEW,LNEW,
     &                       P12,P34,R12,R34)

      RETURN
      END

