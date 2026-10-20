C 8 DEC 16 -YNS - META-GGA FUNCTIONALS FOR TDDFT:MN12/MN15 MODULE
C
C
C*MODULE TDDXCE  *DECK TD_MN12CF_SecD
      SUBROUTINE TD_MN12CF_SecD(EC0,EC,DRHOI_,DTAUI_,
     *               XYZWGT,SCALE,RHOMIN,ILENG,NPTGRD,IRU,IST,IEND,IJZY)
!>
      USE CDERIX, ONLY: IRA, IRB, IGA, IGB, IGC, ITA, ITB, IRARA,
     *  IRARB, IRAGA, IRAGB, IRAGC, IRATA, IRATB, IRBRB, IRBGA, IRBGB,
     *  IRBGC, IRBTA, IRBTB, IGAGA, IGAGB, IGAGC, IGATA, IGATB, IGBGB,
     *  IGBGC, IGBTA, IGBTB, IGCGC, IGCTA, IGCTB, ITATA, ITATB, ITBTB
!> @brief evaluates the second derivatives of the correlation
!> part of the MN12 and MN15 suites of functionals
!>
!> @detail
!> ijzy - 1 MN12-L
!> ijzy - 2 MN12-SX
!> ijzy - 3 MN15-L
!> ijzy - 4 MN15
!> Hence, suppose total grids is 1012, the routine every time
!> will evaluate 200 grids, so the first time, IST=1, IEND=200
!> next time, IST=201, IEND=400, so on so forth, until the fifth
!> time, IST=801, IEND=1000, the ILENG always equals 200
!> at the 6th time, IST=1001, IEND=1012, ILENG=200 as well
!>
!> RHOA=DRHOI_(IPT,1,1), RHOB=DRHOI_(IPT,1,IRU)
!> GRADXA=DRHOI_(IPT,2,1);GRADYA=DRHOI_(IPT,3,1);GRADZA=DRHOI_(IPT,3,1)*
!> TauA=DTAUI_(IPT,1), TauB=DTAUI_(IPT,IRU)
!>
!> OUTPUT:
!> XGRD    - EXCHANGE FUNCTIONAL VALUES
!> DF1     - FIRST DERIVATIVES WITH RESPECT TO RA, RB, GRADA, GRADB, TAUA, TAUB
!>
!> INPUT:
!> RA,B    - SPIN DENSITIES
!> GRADA,B - SPIN DENSITY GRADIENTS
!> GRDAA   = GRADXA*GRADXA+GRADYA*GRADYA+GRADZA*GRADZA
!> GRDBB   = GRADXB*GRADXB+GRADYB*GRADYB+GRADZB*GRADZB
!> GRDAB   = GRADXA*GRADXB+GRADYA*GRADYB+GRADZA*GRADZB
!> TAUA,B  - SPIN KINETIC ENERGY DENSITIES
!> NGRID   - NUMBER OF GRIDS
!> @author Yinan Shu
!> @param EX0 is the energy
!> @param EX is a vector that shows the first and second derivatives of energy
!> @param DRHOI_ is a vector of Rho A and Rho B
!> @param DTAUI_ is a vector of Tau A and Tau B
!> @param IRU is a indicator of RHF and UHF, for RHF, IRU=1, for UHF, IRU=2
!> @param XYZWGT is the weight of grid
!> @param EX0 and XYZWGT has length ILENG (I believe it equals 200)
!> @param NPTGRD is the total number of grids
!> @param IST is the first number of current grids
!> @param IEND is the last number of current grids
!>
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION DRHOI_(NPTGRD,4,IRU)
      DIMENSION DTAUI_(NPTGRD,IRU)
      DIMENSION EC0(*)
      DIMENSION EC(ILENG,*)
      DIMENSION XYZWGT(*)
      DOUBLE PRECISION GRDAA,GRDBB, GRDAB
      DOUBLE PRECISION GRADXA,GRADYA,GRADZA,GRADXB,GRADYB,GRADZB
      DOUBLE PRECISION LSDA
      DOUBLE PRECISION ECF


      DATA ZERO/0.0D0/, F1/1.0D0/,  F2/2.0D0/,  F3/3.0D0/,
     >       F4/4.0D0/, F5/5.0D0/,  F6/6.0D0/,  F7/7.0D0/,
     >       F8/8.0D0/, F9/9.0D0/,F10/10.0D0/, F11/11.0D0/
      Data One/1.0d0/,Two/2.0d0/,Three/3.0d0/,Four/4.0d0/,Five/5.0d0/,
     $  Six/6.0d0/,Seven/7.0d0/,Eight/8.0d0/

      PARAMETER( PI = 3.1415926535897932384626433832795D0 )

      IF (IJZY.EQ.1) THEN
C      PARAMETERS FOR MN12-L
        AT0 =    8.844610D-01
        AT1 =   -2.202279D-01
        AT2 =    5.701372D+00
        AT3 =   -2.562378D+00
        AT4 =   -9.646827D-01
        AT5 =    1.982183D-01
        AT6 =    1.019976D+01
        AT7 =    9.789352D-01
        AT8 =   -1.512722D+00
        AT9 =    0.000000D+00
        AT10=    0.000000D+00
        AT11=    0.000000D+00
C
        BT0 =    5.323948D-01
        BT1 =   -5.831909D+00
        BT2 =    3.882386D+00
        BT3 =    5.878488D+00
        BT4 =    1.493228D+01
        BT5 =   -1.374636D+01
        BT6 =   -8.492327D+00
        BT7 =   -2.486548D+00
        BT8 =   -1.822346D+01
        BT9 =    0.000000D+00
        BT10=    0.000000D+00
        BT11=    0.000000D+00
C
       ELSEIF (IJZY.EQ.2) THEN
C     PARAMETERS FOR MN12-SX
        AT0 =    7.171161D-01
        AT1 =   -2.380914D+00
        AT2 =    5.793565D+00
        AT3 =   -1.243624D+00
        AT4 =    1.364920D+01
        AT5 =   -2.110812D+01
        AT6 =   -1.598767D+01
        AT7 =    1.429208D+01
        AT8 =    6.149191D+00
        AT9 =    0.000000D+00
        AT10=    0.000000D+00
        AT11=    0.000000D+00
C
        BT0 =    4.663699D-01
        BT1 =   -9.110685D+00
        BT2 =    8.705051D+00
        BT3 =   -1.813949D+00
        BT4 =   -4.147211D-01
        BT5 =   -1.021527D+01
        BT6 =    8.240270D-01
        BT7 =    4.993815D+00
        BT8 =   -2.563930D+01
        BT9 =    0.000000D+00
        BT10=    0.000000D+00
        BT11=    0.000000D+00
C
       ELSEIF (IJZY.EQ.3) THEN
C     PARAMETERS FOR MN15-L
        AT0 =    0.952058087D+00
        AT1 =   -0.756954364D+00
        AT2 =    5.677396094D+00
        AT3 =   -5.017104782D+00
        AT4 =   -5.10654071D+00
        AT5 =   -4.812053335D+00
        AT6 =    3.397640087D+00
        AT7 =    1.980041517D+00
        AT8 =    10.1231046D+00
        AT9 =    0.000000D+00
        AT10=    0.000000D+00
        AT11=    0.000000D+00
C
        BT0 =    0.819504932D+00
        BT1 =   -7.689358913D+00
        BT2 =   -0.70532663D+00
        BT3 =   -0.600096421D+00
        BT4 =    11.03332527D+00
        BT5 =    5.861969337D+00
        BT6 =    8.913865465D+00
        BT7 =    5.74529876D+00
        BT8 =    4.254880837D+00
        BT9 =    0.000000D+00
        BT10=    0.000000D+00
        BT11=    0.000000D+00
C
       ELSEIF (IJZY.EQ.4) THEN
C     PARAMETERS FOR MN15
        AT0 =    1.093250748D+00
        AT1 =   -0.269735037D+00
        AT2 =    6.368997613D+00
        AT3 =   -0.245337101D+00
        AT4 =   -1.587103441D+00
        AT5 =    0.124698862D+00
        AT6 =    1.605819855D+00
        AT7 =    0.466206031D+00
        AT8 =    3.484978654D+00
        AT9 =    0.000000D+00
        AT10=    0.000000D+00
        AT11=    0.000000D+00
C
        BT0 =    1.427424993D+00
        BT1 =   -3.57883682D+00
        BT2 =    7.398727547D+00
        BT3 =    3.927810559D+00
        BT4 =    2.789804639D+00
        BT5 =    4.988320462D+00
        BT6 =    3.079464318D+00
        BT7 =    3.521636859D+00
        BT8 =    4.769671992D+00
        BT9 =    0.000000D+00
        BT10=    0.000000D+00
        BT11=    0.000000D+00
C
      ENDIF

C     ECF    = INITIALIZE THE ENERGY OF THE CORRELATION FUNCTION TO ZERO
      ECF    = ZERO

C Useful variables
      Tiny   = 1.0D-8
      F1o3 = F1/F3
      F5o3 = F5/F3
      F2o3 = F2/F3
      PI34 = F3/(F4*PI)


      DO IPT=IST,IEND

         IIPT = IPT-IST+1
         RHOA = DRHOI_(IPT,1,1)
         RHOB = DRHOI_(IPT,1,IRU)
         TauA = DTAUI_(IPT,1)
         TauB = DTAUI_(IPT,IRU)
         RHO  = RHOA + RHOB
         Tau  = TauA + TauB


      IF((RHO.GT.TINY).AND.(TAU.GT.TINY)) THEN

        GRADXA=DRHOI_(IPT,2,1)
        GRADYA=DRHOI_(IPT,3,1)
        GRADZA=DRHOI_(IPT,4,1)
        GRADXB=DRHOI_(IPT,2,IRU)
        GRADYB=DRHOI_(IPT,3,IRU)
        GRADZB=DRHOI_(IPT,4,IRU)
        GRDAA  = GRADXA*GRADXA+GRADYA*GRADYA+GRADZA*GRADZA
        GRDBB  = GRADXB*GRADXB+GRADYB*GRADYB+GRADZB*GRADZB
        GRDAB  = GRADXA*GRADXB+GRADYA*GRADYB+GRADZA*GRADZB
        Y    = GRDAA + GRDBB + 2.00D+00*GRDAB

C     WRITE(IW,*)'RHOA,RHOB,TAUA,TAUB',RHOA,RHOB,TAUA,TAUB

C     WRITE(IW,*)'GRADXA,GRADYA,GRADZA,GRADXB,GRADYB,GRADZB'
C     $   ,GRADXA,GRADYA,GRADZA,GRADXB,GRADYB,GRADZB

        RS     =(PI34/RHO)**F1O3
        ZETA   =(RHOA-RHOB)/RHO
        TAUUEG = F3*(F3*PI*PI)**F2O3*RHO**F5O3/F10
        TSIG   = TAUUEG/TAU
        WSIG   =(TSIG - F1)/(TSIG + F1)
        W1=Wsig
        W2=Wsig*W1
        W3=Wsig*W2
        W4=Wsig*W3
        W5=Wsig*W4
        W6=Wsig*W5
        W7=Wsig*W6
        W8=Wsig*W7
        W9=Wsig*W8
        W10=Wsig*W9
        W11=Wsig*W10

        FSIG1  =(AT0 + WSIG*(AT1 + WSIG*(AT2 + WSIG*(AT3 + WSIG*(
     >           AT4 + WSIG*(AT5 + WSIG*(AT6 + WSIG*(AT7 + WSIG*(
     >           AT8 + WSIG*(AT9 + WSIG*(AT10+ WSIG*AT11)))))))))))
C
        FSIG2  =(BT0 + WSIG*(BT1 + WSIG*(BT2 + WSIG*(BT3 + WSIG*(
     >           BT4 + WSIG*(BT5 + WSIG*(BT6 + WSIG*(BT7 + WSIG*(
     >           BT8 + WSIG*(BT9 + WSIG*(BT10+WSIG*BT11)))))))))))
C
        GRHO   = SQRT(Y)
        CALL LSDAC(RS,ZETA,POTLC,DLDS,DLDZ,D2LDSS,D2LDSZ,D2LDZZ)
        LSDA = RHO*POTLC
        CALL H0TERM(2,Rho,GRho,Rs,Zeta,PotLC,dLdS,dLdZ,
     $      d2LdSS,d2LdSZ,d2LdZZ,d3LdSSS,d3LdSSZ,d3LdSZZ,d3LdZZZ,
     $      H,dHdR,dHdG,dHdZ,d2HdRR,d2HdRG,d2HdGG,d2HdRZ,d2HdGZ,d2HdZZ,
     $      d3HdRRR,d3HdRRG,d3HdRGG,d3HdGGG,d3HdRRZ,d3HdRGZ,d3HdGGZ,
     $      d3HdRZZ,d3HdGZZ,d3HdZZZ)
        GGA  = RHO*H
        E1   = LSDA*FSIG1
        E2   = GGA*FSIG2
      ENERGY = E1 + E2

C energies are computed

C     WRITE(IW,*) ECF,ENERGY,XYZWGT(IPT)

        ECF             = ECF             + SCALE*ENERGY*XYZWGT(IPT)
        EC0(IIPT)       = EC0(IIPT)       + SCALE*ENERGY

C        WRITE(IW,*) 'ENERGY, XYZWGT(IPT), EC0(IIPT), ECF',ENERGY,
C     $  XYZWGT(IPT), EC0(IIPT), ECF

C
C First Derivative
C
         RSP = -RS/(F3*Rho)
         dZdA = (F1-Zeta)/Rho
         dZdB = (-F1-Zeta)/Rho
         dLdRA = dLdS*RSP + dLdZ*dZdA
         dLdRB = dLdS*RSP + dLdZ*dZdB
         dF1dW = (at1 + F2*at2*W1 + F3*at3*W2
     $             + F4*at4*W3 + F5*at5*W4
     $              + F6*at6*W5 + F7*at7*W6
     $              + F8*at8*W7 + F9*at9*W8
     $              + F10*at10*W9 + F11*at11*W10)
         dF2dW = (bt1 + F2*bt2*W1 + F3*bt3*W2
     $             + F4*bt4*W3 + F5*bt5*W4
     $              + F6*bt6*W5 + F7*bt7*W6
     $              + F8*bt8*W7 + F9*bt9*W8
     $              + F10*bt10*W9 + F11*bt11*W10)

         dWdT = F2/((F1 + Tsig)**F2)
         dTdR = Tsig*F5/(F3*Rho)
         dTdTau = -Tsig/Tau
         dF1dR = dF1dW*dWdT*dTdR
         dF1dTau=dF1dW*dWdT*dTdTau
         dF2dR = dF2dW*dWdT*dTdR
         dF2dTau=dF2dW*dWdT*dTdTau
         dLDdRA = PotLC + Rho*dLdRA
         dLDdRB = PotLC + Rho*dLdRB
         dHdRA = dHdR + dHdZ*dZdA
         dHdRB = dHdR + dHdZ*dZdB
         dGRhodY = F1/(F2*GRho)
         dHdY = dHdG * dGRhodY
         dHdGA = dHdY
         dHdGB = dHdY
         dHdGC = dHdY*F2
         dGGAdRA = H + Rho*dHdRA
         dGGAdRB = H + Rho*dHdRB
         dGGAdGA = Rho*dHdGA
         dGGAdGB = dGGAdGA
         dGGAdGC = Rho*dHdGC
C
         dE1dRA = dLDdRA*Fsig1 + LSDA*dF1dR
         dE1dRB = dLDdRB*Fsig1 + LSDA*dF1dR
         dE1dKA = LSDA*dF1dTau
         dE1dKB = dE1dKA
C
         dE2dRA = dGGAdRA*Fsig2 + GGA*dF2dR
         dE2dRB = dGGAdRB*Fsig2 + GGA*dF2dR
         dE2dKA = GGA*dF2dTau
         dE2dKB = dE2dKA
         dE2dGA = dGGAdGA*Fsig2
         dE2dGB = dGGAdGB*Fsig2
         dE2dGC = dGGAdGC*Fsig2
C
C     STORE FIRST FUNCTIONAL DERIVATIVE
C
            EC (IIPT,IRA) = EC (IIPT,IRA) + SCALE*(DE1DRA + DE2DRA)
            EC (IIPT,IRB) = EC (IIPT,IRB) + SCALE*(DE1DRB + DE2DRB)
            EC (IIPT,IGA) = EC (IIPT,IGA) + SCALE*(DE2DGA)
            EC (IIPT,IGB) = EC (IIPT,IGB) + SCALE*(DE2DGB)
            EC (IIPT,IGC) = EC (IIPT,IGC) + SCALE*(DE2DGC)
            EC (IIPT,ITA) = EC (IIPT,ITA) + SCALE*((DE1DKA + DE2DKA))
            EC (IIPT,ITB) = EC (IIPT,ITB) + SCALE*((DE1DKB + DE2DKA))

C      WRITE(IW,*)'EC(IIPT,IRA),EC(IIPT,IGA),EC(IIPT,IGC),EC(IIPT,ITA)'
C     $  ,EC(IIPT,IRA),EC(IIPT,IGA),EC(IIPT,IGC),EC(IIPT,ITA)

C
C Second Derivatives
C
          d2RSdR = -RSP*F4/(F3*Rho)
          RhoCub = Rho*Rho*Rho
          d2ZdAA = -F4*RhoB / RhoCub
          d2ZdAB = F2*Zeta/(Rho*Rho)
          d2ZdBB = F4*RhoA / RhoCub
C d2PotL
          d2LdAA = d2LdSS*RSP**F2 + F2*d2LdSZ*RSP*dZdA +
     $          dLdS*d2RSdR + d2LdZZ*dZdA**F2 + dLdZ*d2ZdAA
          d2LdAB = d2LdSS*RSP**F2 + d2LdSZ*RSP*(dZdA+dZdB) +
     $          dLdS*d2RSdR + d2LdZZ*dZdA*dZdB + dLdZ*d2ZdAB
          d2LdBB = d2LdSS*RSP**F2 + F2*d2LdSZ*RSP*dZdB +
     $          dLdS*d2RSdR + d2LdZZ*dZdB*dZdB + dLdZ*d2ZdBB
C d2H
          d2HdRAA = d2HdRR + F2*d2HdRZ*dZdA +
     $          d2HdZZ*dZdA*dZdA + dHdZ*d2ZdAA
          d2HdRAB = d2HdRR + d2HdRZ*(dZdA+dZdB) +
     $          d2HdZZ*dZdA*dZdB + dHdZ*d2ZdAB
          d2HdRBB = d2HdRR + F2*d2HdRZ*dZdB +
     $          d2HdZZ*dZdB*dZdB + dHdZ*d2ZdBB
          d2HdRAGA = (d2HdRG + d2HdGZ*dZdA)*dGRhodY
          d2HdRBGA = (d2HdRG + d2HdGZ*dZdB)*dGRhodY
          d2HdRAGB = d2HdRAGA
          d2HdRBGB = d2HdRBGA
          d2HdRAGC = F2*d2HdRAGA
          d2HdRBGC = F2*d2HdRBGA
          d2GRhodY = -F1*dGRhodY/(F2*Y)
          d2HdGAA = d2HdGG*dGRhodY*dGRhodY + dHdG*d2GRhodY
          d2HdGAB = d2HdGAA
          d2HdGBB = d2HdGAA
          d2HdGAC = F2*d2HdGAA
          d2HdGBC = d2HdGAC
          d2HdGCC = F4*d2HdGAA


C d2Fsig
          d2TdTau=-F2*dTdTau/Tau
          d2TdR=dTdR*F2/(F3*Rho)
          d2TdRdTau=-dTdR/Tau
          d2WdT=-F2*dWdT/(F1+Tsig)

          d2F1dW=(F2*at2+F6*at3*W1 + F3*F4*at4*W2
     $         + F4*F5*at5*W3 + F5*F6*at6*W4
     $         + F6*F7*at7*W5 + F7*F8*at8*W6
     $         + F8*F9*at9*W7 + F9*F10*at10*W8
     $         + F10*F11*at11*W9)
          d2F2dW=(F2*bt2+F6*bt3*W1 + F3*F4*bt4*W2
     $         + F4*F5*bt5*W3 + F5*F6*bt6*W4
     $         + F6*F7*bt7*W5 + F7*F8*bt8*W6
     $         + F8*F9*bt9*W7 + F9*F10*bt10*W8
     $         + F10*F11*bt11*W9)

          d2F1dTau = dF1dW*dWdT*d2TdTau
     $              +dF1dW*d2WdT*((dTdTau)*(dTdTau))
     $              +d2F1dW*((dWdT*dTdTau)*(dWdT*dTdTau))
          d2F1dR = dF1dW*dWdT*d2TdR
     $             +dF1dW*d2WdT*((dTdR)**F2)
     $             +d2F1dW*((dWdT*dTdR)**F2)
          d2F1dRdTau = dF1dW*dWdT*d2TdRdTau
     $                +dF1dW*d2WdT*dTdR*dTdTau
     $                +d2F1dW*((dWdT)*(dWdT))*dTdR*dTdTau

          d2F2dTau = dF2dW*dWdT*d2TdTau
     $               +dF2dW*d2WdT*((dTdTau)*(dTdTau))
     $               +d2F2dW*((dWdT*dTdTau)*(dWdT*dTdTau))
          d2F2dR = dF2dW*dWdT*d2TdR
     $             +dF2dW*d2WdT*((dTdR)**F2)
     $             +d2F2dW*((dWdT*dTdR)**F2)
          d2F2dRdTau = dF2dW*dWdT*d2TdRdTau
     $             +dF2dW*d2WdT*dTdR*dTdTau
     $             +d2F2dW*((dWdT)*(dWdT))*dTdR*dTdTau
C d2LSDA
          d2LDdRAA = F2*dLdRA + Rho*d2LdAA
          d2LDdRBB = F2*dLdRB + Rho*d2LdBB
          d2LDdRAB = dLdRA + dLdRB + Rho*d2LdAB
C d2GGA
          d2GGAdRAA = F2*dHdRA + Rho*d2HdRAA
          d2GGAdRBB = F2*dHdRB + Rho*d2HdRBB
          d2GGAdRAB = dHdRA + dHdRB + Rho*d2HdRAB
          d2GGAdRAGA = dHdGA + Rho*d2HdRAGA
          d2GGAdRBGA = dHdGA + Rho*d2HdRBGA
          d2GGAdRAGB = dHdGB + Rho*d2HdRAGB
          d2GGAdRBGB = dHdGB + Rho*d2HdRBGB
          d2GGAdRAGC = dHdGC + Rho*d2HdRAGC
          d2GGAdRBGC = dHdGC + Rho*d2HdRBGC
          d2GGAdGAA = Rho*d2HdGAA
          d2GGAdGAB = Rho*d2HdGAB
          d2GGAdGBB = Rho*d2HdGBB
          d2GGAdGAC = Rho*d2HdGAC
          d2GGAdGBC = Rho*d2HdGBC
          d2GGAdGCC = Rho*d2HdGCC

C d2E1
          d2E1dRAA = d2LDdRAA*Fsig1 + F2*dLDdRA*dF1dR + LSDA*d2F1dR
          d2E1dRBB = d2LDdRBB*Fsig1 + F2*dLDdRB*dF1dR + LSDA*d2F1dR
          d2E1dRAB = d2LDdRAB*Fsig1 + dLDdRA*dF1dR
     $             + dLDdRB*dF1dR + LSDA*d2F1dR
          d2E1dRAKA = dLDdRA*dF1dTau + LSDA*d2F1dRdTau
          d2E1dRBKA = dLDdRB*dF1dTau + LSDA*d2F1dRdTau
          d2E1dRAKB = d2E1dRAKA
          d2E1dRBKB = d2E1dRBKA
          d2E1dKAA = LSDA*d2F1dTau
          d2E1dKBB = d2E1dKAA
          d2E1dKAB = d2E1dKAA
C d2E2
          d2E2dRAA=d2GGAdRAA*Fsig2+F2*dGGAdRA*dF2dR+GGA*d2F2dR
          d2E2dRBB=d2GGAdRBB*Fsig2+F2*dGGAdRB*dF2dR+GGA*d2F2dR
          d2E2dRAB=d2GGAdRAB*Fsig2+dGGAdRA*dF2dR
     $             +dGGAdRB*dF2dR+GGA*d2F2dR
          d2E2dRAKA = dGGAdRA*dF2dTau + GGA*d2F2dRdTau
          d2E2dRBKA = dGGAdRB*dF2dTau + GGA*d2F2dRdTau
          d2E2dRAKB = d2E2dRAKA
          d2E2dRBKB = d2E2dRBKA
          d2E2dKAA = GGA*d2F2dTau
          d2E2dKBB = d2E2dKAA
          d2E2dKAB = d2E2dKAA

C
C     STORE SECOND FUNCTIONAL DERIVATIVE
C
            EC (IIPT,IRARA) = EC (IIPT,IRARA)
     $                        + SCALE*(d2E1dRAA+d2E2dRAA)
            EC (IIPT,IRARB) = EC (IIPT,IRARB)
     $                        + SCALE*(d2E1dRAB+d2E2dRAB)
            EC (IIPT,IRAGA) = EC (IIPT,IRAGA)
     $                        + SCALE*(Fsig2*d2GGAdRAGA+dGGAdGA*dF2dR)
            EC (IIPT,IRAGB) = EC (IIPT,IRAGB)
     $                        + SCALE*(Fsig2*d2GGAdRAGB+dGGAdGB*dF2dR)
            EC (IIPT,IRAGC) = EC (IIPT,IRAGC)
     $                        + SCALE*(Fsig2*d2GGAdRAGC+dGGAdGC*dF2dR)
            EC (IIPT,IRATA) = EC (IIPT,IRATA)
     $                        + SCALE*(d2E1dRAKA+d2E2dRAKA)
            EC (IIPT,IRATB) = EC (IIPT,IRATB)
     $                        + SCALE*(d2E1dRAKB+d2E2dRAKB)
            EC (IIPT,IRBRB) = EC (IIPT,IRBRB)
     $                        + SCALE*(d2E1dRBB+d2E2dRBB)
            EC (IIPT,IRBGA) = EC (IIPT,IRBGA)
     $                        + SCALE*(Fsig2*d2GGAdRBGA+dGGAdGA*dF2dR)
            EC (IIPT,IRBGB) = EC (IIPT,IRBGB)
     $                        + SCALE*(Fsig2*d2GGAdRBGB+dGGAdGB*dF2dR)
            EC (IIPT,IRBGC) = EC (IIPT,IRBGC)
     $                        + SCALE*(Fsig2*d2GGAdRBGC+dGGAdGC*dF2dR)
            EC (IIPT,IRBTA) = EC (IIPT,IRBTA)
     $                        + SCALE*(d2E1dRBKA+d2E2dRBKA)
            EC (IIPT,IRBTB) = EC (IIPT,IRBTB)
     $                        + SCALE*(d2E1dRBKB+d2E2dRBKB)
            EC (IIPT,IGAGA) = EC (IIPT,IGAGA)
     $                        + SCALE*(Fsig2*d2GGAdGAA)
            EC (IIPT,IGAGB) = EC (IIPT,IGAGB)
     $                        + SCALE*(Fsig2*d2GGAdGAB)
            EC (IIPT,IGAGC) = EC (IIPT,IGAGC)
     $                        + SCALE*(Fsig2*d2GGAdGAC)
            EC (IIPT,IGATA) = EC (IIPT,IGATA)
     $                        + SCALE*(dGGAdGA*dF2dTau)
            EC (IIPT,IGATB) = EC (IIPT,IGATB)
     $                        + SCALE*(dGGAdGA*dF2dTau)
            EC (IIPT,IGBGB) = EC (IIPT,IGBGB)
     $                        + SCALE*(Fsig2*d2GGAdGBB)
            EC (IIPT,IGBGC) = EC (IIPT,IGBGC)
     $                        + SCALE*(Fsig2*d2GGAdGBC)
            EC (IIPT,IGBTA) = EC (IIPT,IGBTA)
     $                        + SCALE*(dGGAdGB*dF2dTau)
            EC (IIPT,IGBTB) = EC (IIPT,IGBTB)
     $                        + SCALE*(dGGAdGB*dF2dTau)
            EC (IIPT,IGCGC) = EC (IIPT,IGCGC)
     $                        + SCALE*(Fsig2*d2GGAdGCC)
            EC (IIPT,IGCTA) = EC (IIPT,IGCTA)
     $                        + SCALE*(dGGAdGC*dF2dTau)
            EC (IIPT,IGCTB) = EC (IIPT,IGCTB)
     $                        + SCALE*(dGGAdGC*dF2dTau)
            EC (IIPT,ITATA) = EC (IIPT,ITATA)
     $                        + SCALE*(d2E1dKAA+d2E2dKAA)
            EC (IIPT,ITATB) = EC (IIPT,ITATB)
     $                        + SCALE*(d2E1dKAB+d2E2dKAB)
            EC (IIPT,ITBTB) = EC (IIPT,ITBTB)
     $                        + SCALE*(d2E1dKBB+d2E2dKBB)


        END IF
      END DO
      RETURN
      END

C*MODULE TDDXCE  *DECK TD_MN12XF_SecD
      SUBROUTINE TD_MN12XF_SecD(EX0,EX,DRHOI_,DTAUI_,
     *               XYZWGT,SCALE,RHOMIN,ILENG,NPTGRD,IRU,IST,IEND,IJZY)

      USE XDERIX, ONLY: KRA, KRB, KGA, KGB, KTA, KTB, KRARA, KRAGA,
     *  KRATA, KGAGA, KGATA, KTATA, KRBRB, KRBGB, KRBTB, KGBGB, KGBTB,
     *  KTBTB
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION DRHOI_(NPTGRD,4,IRU)
      DIMENSION DTAUI_(NPTGRD,IRU)
      DIMENSION EX0(*)
      DIMENSION EX(ILENG,*)
      DIMENSION XYZWGT(*)
      DOUBLE PRECISION GRDAA,GRDBB
      DOUBLE PRECISION GRADXA,GRADYA,GRADZA,GRADXB,GRADYB,GRADZB
      DOUBLE PRECISION EXF



      DATA ZERO/0.0D0/, F1/1.0D0/,  F2/2.0D0/,  F3/3.0D0/,
     >       F4/4.0D0/, F5/5.0D0/,  F6/6.0D0/,  F7/7.0D0/,
     >       F8/8.0D0/, F9/9.0D0/,F10/10.0D0/, F11/11.0D0/
      Data One/1.0d0/,Two/2.0d0/,Three/3.0d0/,Four/4.0d0/,Five/5.0d0/,
     $  Six/6.0d0/,Seven/7.0d0/,Eight/8.0d0/

      PARAMETER( PI = 3.1415926535897932384626433832795D0 )

      If (ijzy.eq.1) then
c
c       MN12-L
c
        CC000 =  6.735981D-01
        CC001 = -2.270598D+00
        CC002 = -2.613712D+00
        CC003 =  3.993609D+00
        CC004 =  4.635575D+00
        CC005 =  1.250676D+00
        CC010 =  8.444920D-01
        CC011 = -1.301173D+01
        CC012 = -1.777730D+01
        CC013 = -4.627211D+00
        CC014 =  5.976605D+00
        CC020 =  1.142897D+00
        CC021 = -2.040226D+01
        CC022 = -2.382843D+01
        CC023 =  7.119109D+00
        CC030 = -2.335726D+01
        CC031 = -1.622633D+01
        CC032 =  1.482732D+01
        CC100 =  1.449285D+00
        CC101 =  1.020598D+01
        CC102 =  4.407450D+00
        CC103 = -2.008193D+01
        CC104 = -1.253561D+01
        CC110 = -5.435031D+00
        CC111 =  1.656736D+01
        CC112 =  2.000229D+01
        CC113 = -2.513105D+00
        CC120 =  9.658436D+00
        CC121 = -3.825281D+00
        CC122 = -2.500000D+01
        CC200 = -2.070080D+00
        CC201 = -9.951913D+00
        CC202 =  8.731211D-01
        CC203 =  2.210891D+01
        CC210 =  8.822633D+00
        CC211 =  2.499949D+01
        CC212 =  2.500000D+01
        CC300 =  6.851693D-01
        CC301 = -7.406948D-02
        CC302 = -6.788000D-01

      Else If (ijzy.eq.2) then
c
c       MN12-SX
c
        CC000 =  5.226556D-01
        CC001 = -2.681208D-01
        CC002 = -4.670705D+00
        CC003 =  3.067320D+00
        CC004 =  4.095370D+00
        CC005 =  2.653023D+00
        CC010 =  5.165969D-01
        CC011 = -2.035442D+01
        CC012 = -9.946472D+00
        CC013 =  2.938637D+00
        CC014 =  1.131100D+01
        CC020 =  4.752452D+00
        CC021 = -3.061331D+00
        CC022 = -2.523173D+01
        CC023 =  1.710903D+01
        CC030 = -2.357480D+01
        CC031 = -2.727754D+01
        CC032 =  1.603291D+01
        CC100 =  1.842503D+00
        CC101 =  1.927120D+00
        CC102 =  1.107987D+01
        CC103 = -1.182087D+01
        CC104 = -1.117768D+01
        CC110 = -5.821000D+00
        CC111 =  2.266545D+01
        CC112 =  8.246708D+00
        CC113 = -4.778364D+00
        CC120 =  5.329122D-01
        CC121 = -6.666755D+00
        CC122 =  1.671429D+00
        CC200 = -3.311409D+00
        CC201 =  3.415913D-01
        CC202 = -6.413076D+00
        CC203 =  1.038584D+01
        CC210 =  9.026277D+00
        CC211 =  1.929689D+01
        CC212 =  2.669232D+01
        CC300 =  1.517278D+00
        CC301 = -3.442503D+00
        CC302 =  1.100161D+00

      Else If (ijzy.eq.3) then
c
c       MN15-L
c
        CC000 =  0.670864162d0
        CC001 = -0.822003903d0
        CC002 = -1.022407046d0
        CC003 =  1.689460986d0
        CC004 = -0.00562032d0
        CC005 = -0.110293849d0
        CC010 =  0.972245178d0
        CC011 = -6.697641991d0
        CC012 = -4.322814495d0
        CC013 = -6.786641376d0
        CC014 = -5.687461462d0
        CC020 =  9.419643818d0
        CC021 =  11.83939406d0
        CC022 =  5.086951311d0
        CC023 =  4.302369948d0
        CC030 = -8.07344065d0
        CC031 =  2.429988978d0
        CC032 =  11.09485698d0
        CC100 =  1.247333909d0
        CC101 =  3.700485291d0
        CC102 =  0.867791614d0
        CC103 = -0.591190518d0
        CC104 = -0.295305435d0
        CC110 = -5.825759145d0
        CC111 =  2.537532196d0
        CC112 =  3.143390933d0
        CC113 =  2.939126332d0
        CC120 =  0.599342114d0
        CC121 =  2.241702738d0
        CC122 =  2.035713838d0
        CC200 = -1.525344043d0
        CC201 = -2.325875691d0
        CC202 =  1.141940663d0
        CC203 = -1.563165026d0
        CC210 =  7.882032871d0
        CC211 =  11.93400684d0
        CC212 =  9.852928303d0
        CC300 =  0.584030245d0
        CC301 = -0.720941131d0
        CC302 = -2.836037078d0

      Else If (ijzy.eq.4) then
c
c        MN15
        CC000  =  0.073852235d0
        CC001  = -0.839976156d0
        CC002  = -3.082660125d0
        CC003  = -1.02881285d0
        CC004  = -0.811697255d0
        CC005  = -0.063404387d0
        CC010  =  2.54805518d0
        CC011  = -5.031578906d0
        CC012  =  0.31702159d0
        CC013  =  2.981868205d0
        CC014  = -0.749503735d0
        CC020  =  0.231825661d0
        CC021  =  1.261961411d0
        CC022  =  1.665920815d0
        CC023  =  7.483304941d0
        CC030  = -2.544245723d0
        CC031  =  1.384720031d0
        CC032  =  6.902569885d0
        CC100  =  1.657399451d0
        CC101  =  2.98526709d0
        CC102  =  6.89391326d0
        CC103  =  2.489813993d0
        CC104  =  1.454724691d0
        CC110  = -5.054324071d0
        CC111  =  2.35273334d0
        CC112  =  1.299104132d0
        CC113  =  1.203168217d0
        CC120  =  0.121595877d0
        CC121  =  8.048348238d0
        CC122  =  21.91203659d0
        CC200  = -1.852335832d0
        CC201  = -3.4722735d0
        CC202  = -1.564591493d0
        CC203  = -2.29578769d0
        CC210  =  3.666482991d0
        CC211  =  10.87074639d0
        CC212  =  9.696691388d0
        CC300  =  0.630701064d0
        CC301  = -0.505825216d0
        CC302  = -3.562354535d0

       End If
C
C     EXF    = INITIALIZE THE ENERGY OF THE EXCHANGE FUNCTION TO ZERO
      EXF    = ZERO

C Working variables
C
        Tiny = 1.0d-8
        G   = 0.004d+0
        ome = 2.5d+0
        F28   = Four * Seven
        F2o3  = Two / Three
        F3o2  = Three / Two
        F3o5  = Three/Five
        F1o3  = One / Three
        F4o3  = Four / Three
        F5o3  = Five/Three
        F8o3  = Eight / Three
        F28o9 = F28 / Three / Three
        Ax    = -F3o2*(F4o3*PI)**(-F1o3)

      DO IPT=IST,IEND

         IIPT = IPT-IST+1
         RHOA = DRHOI_(IPT,1,1)
         RHOB = DRHOI_(IPT,1,IRU)
         TauA = DTAUI_(IPT,1)*F2
         TauB = DTAUI_(IPT,IRU)*F2

        GRADXA=DRHOI_(IPT,2,1)
        GRADYA=DRHOI_(IPT,3,1)
        GRADZA=DRHOI_(IPT,4,1)
        GRADXB=DRHOI_(IPT,2,IRU)
        GRADYB=DRHOI_(IPT,3,IRU)
        GRADZB=DRHOI_(IPT,4,IRU)
        GRDAA  = GRADXA*GRADXA+GRADYA*GRADYA+GRADZA*GRADZA
        GRDBB  = GRADXB*GRADXB+GRADYB*GRADYB+GRADZB*GRADZB


C       WRITE(IW,*)'RHOA,RHOB,TAUA,TAUB',RHOA,RHOB,TAUA,TAUB

C       WRITE(IW,*)'GRDAA,GRDBB',GRDAA,GRDBB

C       WRITE(IW,*)'RHOMIN',RHOMIN


      IF(IRU.EQ.1) THEN

      IF(RHOA.GT.RHOMIN) THEN
          pX    = RHOA
          GamX2 = GRDAA
          S2    = GamX2*pX**(-F8o3)
          E = Ax*pX**F4o3
          Tau    = TAUA
          TauUEG =F3o5*((F6*PI*PI)**F2o3)*(pX**F5o3)
          Tsig =TauUEG/Tau
          FV = ome*pX**F1o3/(One+ome*pX**F1o3)
          FU = G*S2/(One+G*S2)
          FW = (Tsig - One)/(Tsig + One)

          FMN12 = CC000 + FW*CC001 + FW**2*CC002 + FW**3*CC003 +
     $          FW**4*CC004 + FW**5*CC005 +
     $          FU*CC010 + FW*FU*CC011 + FW**2*FU*CC012 +
     $          FW**3*FU*CC013 +
     $          FW**4*FU*CC014 + FU**2*CC020 +
     $          FW*FU**2*CC021 + FW**2*FU**2*CC022 +
     $          FW**3*FU**2*CC023 +
     $          FU**3*CC030 + FW*FU**3*CC031 +
     $          FW**2*FU**3*CC032 + FV*CC100 + FV*FW*CC101 +
     $          FV*FW**2*CC102 + FV*FW**3*CC103 +
     $          FV*FW**4*CC104 + FV*FU*CC110 + FV*FW*FU*CC111 +
     $          FV*FW**2*FU*CC112 + FV*FW**3*FU*CC113 +
     $          FV*FU**2*CC120 + FV*FW*FU**2*CC121 +
     $          FV*FW**2*FU**2*CC122 +
     $          FV**2*CC200 + FV**2*FW*CC201 +
     $          FV**2*FW**2*CC202 + FV**2*FW**3*CC203 +
     $          FV**2*FU*CC210 +
     $          FV**2*FW*FU*CC211 +
     $          FV**2*FW**2*FU*CC212 + FV**3*CC300 +
     $          FV**3*FW*CC301 +
     $          FV**3*FW**2*CC302

        ENERGY = E * FMN12

C energies are computed

        EXF          = EXF          + SCALE*ENERGY*F2*XYZWGT(IPT)
        EX0(IIPT)    = EX0(IIPT)    + SCALE*ENERGY*F2

C TESTVALUE=ENERGY*F2

C        WRITE(IW,*) 'TESTVALUE, XYZWGT(IPT), EX0(IIPT),EXF',TESTVALUE,
C     $  XYZWGT(IPT),EX0(IIPT), EXF
C        WRITE(IW,*) 'EXF,ENERGY,XYZWGT(IPT)',EXF,ENERGY,XYZWGT(IPT)

C First Derivatives
          ER = F4o3*E/pX
          S    = Sqrt(S2)
          GamX = Sqrt(GamX2)
          SR   = -F4o3*S/pX
          SG   = S/GamX
          TR  = ((F6*PI*PI)**F2o3)*(pX**F2o3)/Tau
          TK  = -TauUEG/Tau**2
          FUS = Two*G*S/((One+G*S*S)**2)
          FVR =  (ome/(Three*pX**F2o3))*
     $           (One+ome*pX**F1o3)**(-Two)
          FWT = Two/((One + Tsig)**2)
          FUR = FUS*SR
          FUG = FUS*SG
          FWR = FWT*TR
          FWK = FWT*TK

C        WRITE(IW,*)'FUS,SG,FWT,TK',FUS,SG,FWT,TK
C        WRITE(IW,*)'S,GAMX,TSIG',S,GAMX,TSIG
C        WRITE(IW,*)'Im Calling MN12XF now computing these parameters:',
C     $  ER,S,GAMX,SR,SG,TR,TK,FUS,FVR,FWT,FUR,FUG,FWR,FWK




          FMN12FV = CC100+FW*CC101+FW**2*CC102+FW**3*CC103+
     $              FW**4*CC104+FU*CC110+
     $              FW*FU*CC111+FW**2*FU*CC112+FW**3*FU*CC113+
     $              FU**2*CC120+FW*FU**2*CC121+
     $              FW**2*FU**2*CC122+Two*FV*CC200+Two*FV*FW*CC201+
     $              Two*FV*FW**2*CC202+Two*FV*FW**3*CC203+
     $              Two*FV*FU*CC210+Two*FV*FW*FU*CC211+
     $              Two*FV*FW**2*FU*CC212+
     $              Three*FV**2*CC300+Three*FV**2*FW*CC301+
     $              Three*FV**2*FW**2*CC302
          FMN12FU = CC010+FW*CC011+FW**2*CC012+FW**3*CC013+
     $              FW**4*CC014+Two*FU*CC020+
     $              Two*FW*FU*CC021+Two*FW**2*FU*CC022+
     $              Two*FW**3*FU*CC023+
     $              Three*FU**2*CC030+Three*FW*FU**2*CC031+
     $              Three*FW**2*FU**2*CC032+FV*CC110+FV*FW*CC111+
     $              FV*FW**2*CC112+FV*FW**3*CC113+
     $              Two*FV*FU*CC120+Two*FV*FW*FU*CC121+
     $              Two*FV*FW**2*FU*CC122+
     $              FV**2*CC210+FV**2*FW*CC211+
     $              FV**2*FW**2*CC212
          FMN12FW = CC001+Two*FW*CC002+Three*FW**2*CC003+
     $              Four*FW**3*CC004+Five*FW**4*CC005+FU*CC011+
     $              Two*FW*FU*CC012+Three*FW**2*FU*CC013+
     $              Four*FW**3*FU*CC014+
     $              FU**2*CC021+Two*FW*FU**2*CC022+
     $              Three*FW**2*FU**2*CC023+FU**3*CC031+
     $              Two*FW*FU**3*CC032+
     $              FV*CC101+Two*FV*FW*CC102+
     $              Three*FV*FW**2*CC103+Four*FV*FW**3*CC104+
     $              FV*FU*CC111+
     $              Two*FV*FW*FU*CC112+Three*FV*FW**2*FU*CC113+
     $              FV*FU**2*CC121+Two*FV*FW*FU**2*CC122+FV**2*CC201+
     $              Two*FV**2*FW*CC202+
     $              Three*FV**2*FW**2*CC203+FV**2*FU*CC211+
     $              Two*FV**2*FW*FU*CC212+FV**3*CC301+
     $              Two*FV**3*FW*CC302
c
          FMN12R  = FMN12FV*FVR+FMN12FU*FUR+FMN12FW*FWR
          FMN12G  = FMN12FU*FUG
          FMN12K  = FMN12FW*FWK
C
C     WRITE(IW,*)'FMN12R,FMN12G,FMN12K',FMN12R,FMN12G,FMN12K
C     WRITE(IW,*)'FMN12FU,FUG,FMN12FW,FWK',FMN12FU,FUG,FMN12FW,FWK

C
C     STORE FIRST FUNCTIONAL DERIVATIVE
C

            EX (IIPT,KRA) = EX (IIPT,KRA) + SCALE*(ER*FMN12+E*FMN12R)
            EX (IIPT,KGA) = EX (IIPT,KGA) + SCALE*E*FMN12G/(Two*GamX)
            EX (IIPT,KTA) = EX (IIPT,KTA) + SCALE*E*FMN12K*F2
            EX (IIPT,KRB) = EX (IIPT,KRB) + SCALE*(ER*FMN12+E*FMN12R)
            EX (IIPT,KGB) = EX (IIPT,KGB) + SCALE*E*FMN12G/(Two*GamX)
            EX (IIPT,KTB) = EX (IIPT,KTB) + SCALE*E*FMN12K*F2

C      WRITE(IW,*)'EX (IIPT,KRA),EX (IIPT,KGA),EX (IIPT,KTA)'
C     $  ,EX (IIPT,KRA),EX (IIPT,KGA),EX (IIPT,KTA)

C Second Derivatives
          ERR = ER/(Three*pX)
          SRR = F28o9*S/(pX*pX)
          SRG = SR/GamX
          TRR = F2o3*TR/pX
          TRK = -TR/tau
          TKK = -Two*TK/tau
          FWTT = -Four/(Tsig+One)**3
          FUSS = Two*G*(One-Three*G*S*S)/((One+G*S*S)**3)
          FVRR = (-Two*ome*(One + Two*ome*pX**F1o3))/
     $           (F9*(One + ome*pX**F1o3)**3*pX**F5o3)
          FURR = FUSS*SR*SR +
     $           FUS*SRR
          FUGG = FUSS*SG*SG
          FURG = FUSS*SR*SG +
     $           FUS*SRG
          FWRR = FWTT*TR*TR +
     $           FWT*TRR
          FWRK = FWTT*TK*TR +
     $           FWT*TRK
          FWKK = FWTT*TK*TK +
     $           FWT*TKK
c
          FMN12FVFV = Two*CC200+Two*FW*CC201+Two*FW**2*CC202+
     $                Two*FW**3*CC203+Two*FU*CC210+Two*FU*FW*CC211+
     $                Two*FU*FW**2*CC212+Six*FV*CC300+
     $                Six*FV*FW*CC301+Six*FV*FW**2*CC302
          FMN12FUFU = Two*CC020+Two*FW*CC021+Two*FW**2*CC022+
     $                Two*FW**3*CC023+Six*FU*CC030+Six*FU*FW*CC031+
     $                Six*FU*FW**2*CC032+Two*FV*CC120+Two*FV*FW*CC121+
     $                Two*FV*FW**2*CC122
          FMN12FWFW = Two*CC002+Six*FW*CC003+12*FW**2*CC004+
     $                20.0D+00*FW**3*CC005+Two*FU*CC012+
     $                Six*FU*FW*CC013+
     $                12.0D+00*FU*FW**2*CC014+Two*FU**2*CC022+
     $                Six*FU**2*FW*CC023+Two*FU**3*CC032+Two*FV*CC102+
     $                Six*FV*FW*CC103+12.0D+00*FV*FW**2*CC104+
     $                Two*FU*FV*CC112+Six*FU*FV*FW*CC113+
     $                Two*FU**2*FV*CC122+Two*FV**2*CC202+
     $                Six*FV**2*FW*CC203+Two*FU*FV**2*CC212+
     $                Two*FV**3*CC302
          FMN12FVFU = CC110+FW*CC111+FW**2*CC112+FW**3*CC113+
     $                Two*FU*CC120+Two*FU*FW*CC121+
     $                Two*FU*FW**2*CC122+Two*FV*CC210+
     $                Two*FV*FW*CC211+Two*FV*FW**2*CC212
          FMN12FVFW = CC101+Two*FW*CC102+Three*FW**2*CC103+
     $                Four*FW**3*CC104+FU*CC111+Two*FU*FW*CC112+
     $                Three*FU*FW**2*CC113+FU**2*CC121+
     $                Two*FU**2*FW*CC122+Two*FV*CC201+Four*FV*FW*CC202+
     $                Six*FV*FW**2*CC203+Two*FU*FV*CC211+
     $                Four*FU*FV*FW*CC212+Three*FV**2*CC301+
     $                Six*FV**2*FW*CC302
          FMN12FUFW = CC011+Two*FW*CC012+Three*FW**2*CC013+
     $                Four*FW**3*CC014+Two*FU*CC021+Four*FU*FW*CC022+
     $                Six*FU*FW**2*CC023+Three*FU**2*CC031+
     $                Six*FU**2*FW*CC032+FV*CC111+Two*FV*FW*CC112+
     $                Three*FV*FW**2*CC113+Two*FU*FV*CC121+
     $                Four*FU*FV*FW*CC122+FV**2*CC211+
     $                Two*FV**2*FW*CC212
c
          FMN12RR = FMN12FVFV*FVR*FVR +
     $              FMN12FVFU*FVR*FUR +
     $              FMN12FVFW*FVR*FWR +
     $              FMN12FV*FVRR +
     $              FMN12FVFU*FUR*FVR +
     $              FMN12FUFU*FUR*FUR +
     $              FMN12FUFW*FUR*FWR +
     $              FMN12FU*FURR +
     $              FMN12FVFW*FWR*FVR +
     $              FMN12FUFW*FWR*FUR +
     $              FMN12FWFW*FWR*FWR +
     $              FMN12FW*FWRR
          FMN12RG = FMN12FVFU*FUG*FVR +
     $              FMN12FUFU*FUG*FUR +
     $              FMN12FUFW*FUG*FWR +
     $              FMN12FU*FURG
          FMN12RK = FMN12FVFW*FWK*FVR +
     $              FMN12FUFW*FWK*FUR +
     $              FMN12FWFW*FWK*FWR +
     $              FMN12FW*FWRK
          FMN12GG = FMN12FUFU*FUG*FUG+
     $              FMN12FU*FUGG
          FMN12KK = FMN12FWFW*FWK*FWK+
     $              FMN12FW*FWKK
          FMN12GK = FMN12FUFW*FWK*FUG
c
C
C     STORE SECOND FUNCTIONAL DERIVATIVE
C
            EX (IIPT,KRARA) = EX (IIPT,KRARA)
     $          + SCALE*(ERR*FMN12+Two*ER*FMN12R+E*FMN12RR)
            EX (IIPT,KRAGA) = EX (IIPT,KRAGA)
     $          + SCALE*(ER*FMN12G+E*FMN12RG)/(Two*GamX)
            EX (IIPT,KRATA) = EX (IIPT,KRATA)
     $          + SCALE*(ER*FMN12K+E*FMN12RK)*F2
            EX (IIPT,KGAGA) = EX (IIPT,KGAGA)
     $      + SCALE*(E/(Four*GamX2*GamX))*(GamX*FMN12GG-FMN12G)
            EX (IIPT,KGATA) = EX (IIPT,KGATA)
     $          + SCALE*E*FMN12GK/(Two*GamX)*F2
            EX (IIPT,KTATA) = EX (IIPT,KTATA)
     $          + SCALE*(E*FMN12KK)*F4
            EX (IIPT,KRBRB) = EX (IIPT,KRBRB)
     $          + SCALE*(ERR*FMN12+Two*ER*FMN12R+E*FMN12RR)
            EX (IIPT,KRBGB) = EX (IIPT,KRBGB)
     $          + SCALE*(ER*FMN12G+E*FMN12RG)/(Two*GamX)
            EX (IIPT,KRBTB) = EX (IIPT,KRBTB)
     $          + SCALE*(ER*FMN12K+E*FMN12RK)*F2
            EX (IIPT,KGBGB) = EX (IIPT,KGBGB)
     $      + SCALE*(E/(Four*GamX2*GamX))*(GamX*FMN12GG - FMN12G)
            EX (IIPT,KGBTB) = EX (IIPT,KGBTB)
     $          + SCALE*(E*FMN12GK)/(Two*GamX)*F2
            EX (IIPT,KTBTB) = EX (IIPT,KTBTB)
     $          + SCALE*(E*FMN12KK)*F4



C      WRITE(IW,*)'RARA,RAHA,RATA,GAGA,GATA,TATA',EX (IIPT,KRARA),
C     $  EX (IIPT,KRAGA),EX (IIPT,KRATA),EX (IIPT,KGAGA),
C     $  EX (IIPT,KGATA),EX (IIPT,KTATA)


      ENDIF
      ENDIF



      IF(IRU.EQ.2) THEN

C ---Alpha component---

      IF((RHOA.GT.TINY).AND.(TAUA.GT.TINY)) THEN
          pX    = RHOA
          GamX2 = GRDAA
          S2    = GamX2*pX**(-F8o3)
          E = Ax*pX**F4o3
          Tau    = TAUA
          TauUEG =F3o5*((F6*PI*PI)**F2o3)*(pX**F5o3)
          Tsig =TauUEG/Tau
          FV = ome*pX**F1o3/(One+ome*pX**F1o3)
          FU = G*S2/(One+G*S2)
          FW = (Tsig - One)/(Tsig + One)

          FMN12 = CC000 + FW*CC001 + FW**2*CC002 + FW**3*CC003 +
     $          FW**4*CC004 + FW**5*CC005 +
     $          FU*CC010 + FW*FU*CC011 + FW**2*FU*CC012 +
     $          FW**3*FU*CC013 +
     $          FW**4*FU*CC014 + FU**2*CC020 +
     $          FW*FU**2*CC021 + FW**2*FU**2*CC022 +
     $          FW**3*FU**2*CC023 +
     $          FU**3*CC030 + FW*FU**3*CC031 +
     $          FW**2*FU**3*CC032 + FV*CC100 + FV*FW*CC101 +
     $          FV*FW**2*CC102 + FV*FW**3*CC103 +
     $          FV*FW**4*CC104 + FV*FU*CC110 + FV*FW*FU*CC111 +
     $          FV*FW**2*FU*CC112 + FV*FW**3*FU*CC113 +
     $          FV*FU**2*CC120 + FV*FW*FU**2*CC121 +
     $          FV*FW**2*FU**2*CC122 +
     $          FV**2*CC200 + FV**2*FW*CC201 +
     $          FV**2*FW**2*CC202 + FV**2*FW**3*CC203 +
     $          FV**2*FU*CC210 +
     $          FV**2*FW*FU*CC211 +
     $          FV**2*FW**2*FU*CC212 + FV**3*CC300 +
     $          FV**3*FW*CC301 +
     $          FV**3*FW**2*CC302

      ENERGY = E * FMN12
C energies are computed
        EXF          = EXF          + SCALE*ENERGY*XYZWGT(IPT)
        EX0(IIPT)    = EX0(IIPT)    + SCALE*ENERGY

C First Derivatives
          ER = F4o3*E/pX
          S    = Sqrt(S2)
          GamX = Sqrt(GamX2)
          SR   = -F4o3*S/pX
          SG   = S/GamX
          TR  = ((F6*PI*PI)**F2o3)*(pX**F2o3)/Tau
          TK  = -TauUEG/Tau**2
          FUS = Two*G*S/((One+G*S*S)**2)
          FVR =  (ome/(Three*pX**F2o3))*
     $           (One+ome*pX**F1o3)**(-Two)
          FWT = Two/((One + Tsig)**2)
          FUR = FUS*SR
          FUG = FUS*SG
          FWR = FWT*TR
          FWK = FWT*TK
          FMN12FV = CC100+FW*CC101+FW**2*CC102+FW**3*CC103+
     $              FW**4*CC104+FU*CC110+
     $              FW*FU*CC111+FW**2*FU*CC112+FW**3*FU*CC113+
     $              FU**2*CC120+FW*FU**2*CC121+
     $              FW**2*FU**2*CC122+Two*FV*CC200+Two*FV*FW*CC201+
     $              Two*FV*FW**2*CC202+Two*FV*FW**3*CC203+
     $              Two*FV*FU*CC210+Two*FV*FW*FU*CC211+
     $              Two*FV*FW**2*FU*CC212+
     $              Three*FV**2*CC300+Three*FV**2*FW*CC301+
     $              Three*FV**2*FW**2*CC302
          FMN12FU = CC010+FW*CC011+FW**2*CC012+FW**3*CC013+
     $              FW**4*CC014+Two*FU*CC020+
     $              Two*FW*FU*CC021+Two*FW**2*FU*CC022+
     $              Two*FW**3*FU*CC023+
     $              Three*FU**2*CC030+Three*FW*FU**2*CC031+
     $              Three*FW**2*FU**2*CC032+FV*CC110+FV*FW*CC111+
     $              FV*FW**2*CC112+FV*FW**3*CC113+
     $              Two*FV*FU*CC120+Two*FV*FW*FU*CC121+
     $              Two*FV*FW**2*FU*CC122+
     $              FV**2*CC210+FV**2*FW*CC211+
     $              FV**2*FW**2*CC212
          FMN12FW = CC001+Two*FW*CC002+Three*FW**2*CC003+
     $              Four*FW**3*CC004+Five*FW**4*CC005+FU*CC011+
     $              Two*FW*FU*CC012+Three*FW**2*FU*CC013+
     $              Four*FW**3*FU*CC014+
     $              FU**2*CC021+Two*FW*FU**2*CC022+
     $              Three*FW**2*FU**2*CC023+FU**3*CC031+
     $              Two*FW*FU**3*CC032+
     $              FV*CC101+Two*FV*FW*CC102+
     $              Three*FV*FW**2*CC103+Four*FV*FW**3*CC104+
     $              FV*FU*CC111+
     $              Two*FV*FW*FU*CC112+Three*FV*FW**2*FU*CC113+
     $              FV*FU**2*CC121+Two*FV*FW*FU**2*CC122+FV**2*CC201+
     $              Two*FV**2*FW*CC202+
     $              Three*FV**2*FW**2*CC203+FV**2*FU*CC211+
     $              Two*FV**2*FW*FU*CC212+FV**3*CC301+
     $              Two*FV**3*FW*CC302
c
          FMN12R  = FMN12FV*FVR+FMN12FU*FUR+FMN12FW*FWR
          FMN12G  = FMN12FU*FUG
          FMN12K  = FMN12FW*FWK
C
C     STORE FIRST FUNCTIONAL DERIVATIVE
C

            EX (IIPT,KRA) = EX (IIPT,KRA) + SCALE*(ER*FMN12+E*FMN12R)
            EX (IIPT,KGA) = EX (IIPT,KGA) + SCALE*E*FMN12G/(Two*GamX)
            EX (IIPT,KTA) = EX (IIPT,KTA) + SCALE*E*FMN12K*F2
C Second Derivatives
          ERR = ER/(Three*pX)
          SRR = F28o9*S/(pX*pX)
          SRG = SR/GamX
          TRR = F2o3*TR/pX
          TRK = -TR/tau
          TKK = -Two*TK/tau
          FWTT = -Four/(Tsig+One)**3
          FUSS = Two*G*(One-Three*G*S*S)/((One+G*S*S)**3)
          FVRR = (-Two*ome*(One + Two*ome*pX**F1o3))/
     $           (F9*(One + ome*pX**F1o3)**3*pX**F5o3)
          FURR = FUSS*SR*SR +
     $           FUS*SRR
          FUGG = FUSS*SG*SG
          FURG = FUSS*SR*SG +
     $           FUS*SRG
          FWRR = FWTT*TR*TR +
     $           FWT*TRR
          FWRK = FWTT*TK*TR +
     $           FWT*TRK
          FWKK = FWTT*TK*TK +
     $           FWT*TKK
c
          FMN12FVFV = Two*CC200+Two*FW*CC201+Two*FW**2*CC202+
     $                Two*FW**3*CC203+Two*FU*CC210+Two*FU*FW*CC211+
     $                Two*FU*FW**2*CC212+Six*FV*CC300+
     $                Six*FV*FW*CC301+Six*FV*FW**2*CC302
          FMN12FUFU = Two*CC020+Two*FW*CC021+Two*FW**2*CC022+
     $                Two*FW**3*CC023+Six*FU*CC030+Six*FU*FW*CC031+
     $                Six*FU*FW**2*CC032+Two*FV*CC120+Two*FV*FW*CC121+
     $                Two*FV*FW**2*CC122
          FMN12FWFW = Two*CC002+Six*FW*CC003+12*FW**2*CC004+
     $                20.0D+00*FW**3*CC005+Two*FU*CC012+
     $                Six*FU*FW*CC013+
     $                12.0D+00*FU*FW**2*CC014+Two*FU**2*CC022+
     $                Six*FU**2*FW*CC023+Two*FU**3*CC032+Two*FV*CC102+
     $                Six*FV*FW*CC103+12.0D+00*FV*FW**2*CC104+
     $                Two*FU*FV*CC112+Six*FU*FV*FW*CC113+
     $                Two*FU**2*FV*CC122+Two*FV**2*CC202+
     $                Six*FV**2*FW*CC203+Two*FU*FV**2*CC212+
     $                Two*FV**3*CC302
          FMN12FVFU = CC110+FW*CC111+FW**2*CC112+FW**3*CC113+
     $                Two*FU*CC120+Two*FU*FW*CC121+
     $                Two*FU*FW**2*CC122+Two*FV*CC210+
     $                Two*FV*FW*CC211+Two*FV*FW**2*CC212
          FMN12FVFW = CC101+Two*FW*CC102+Three*FW**2*CC103+
     $                Four*FW**3*CC104+FU*CC111+Two*FU*FW*CC112+
     $                Three*FU*FW**2*CC113+FU**2*CC121+
     $                Two*FU**2*FW*CC122+Two*FV*CC201+Four*FV*FW*CC202+
     $                Six*FV*FW**2*CC203+Two*FU*FV*CC211+
     $                Four*FU*FV*FW*CC212+Three*FV**2*CC301+
     $                Six*FV**2*FW*CC302
          FMN12FUFW = CC011+Two*FW*CC012+Three*FW**2*CC013+
     $                Four*FW**3*CC014+Two*FU*CC021+Four*FU*FW*CC022+
     $                Six*FU*FW**2*CC023+Three*FU**2*CC031+
     $                Six*FU**2*FW*CC032+FV*CC111+Two*FV*FW*CC112+
     $                Three*FV*FW**2*CC113+Two*FU*FV*CC121+
     $                Four*FU*FV*FW*CC122+FV**2*CC211+
     $                Two*FV**2*FW*CC212
c
          FMN12RR = FMN12FVFV*FVR*FVR +
     $              FMN12FVFU*FVR*FUR +
     $              FMN12FVFW*FVR*FWR +
     $              FMN12FV*FVRR +
     $              FMN12FVFU*FUR*FVR +
     $              FMN12FUFU*FUR*FUR +
     $              FMN12FUFW*FUR*FWR +
     $              FMN12FU*FURR +
     $              FMN12FVFW*FWR*FVR +
     $              FMN12FUFW*FWR*FUR +
     $              FMN12FWFW*FWR*FWR +
     $              FMN12FW*FWRR
          FMN12RG = FMN12FVFU*FUG*FVR +
     $              FMN12FUFU*FUG*FUR +
     $              FMN12FUFW*FUG*FWR +
     $              FMN12FU*FURG
          FMN12RK = FMN12FVFW*FWK*FVR +
     $              FMN12FUFW*FWK*FUR +
     $              FMN12FWFW*FWK*FWR +
     $              FMN12FW*FWRK
          FMN12GG = FMN12FUFU*FUG*FUG+
     $              FMN12FU*FUGG
          FMN12KK = FMN12FWFW*FWK*FWK+
     $              FMN12FW*FWKK
          FMN12GK = FMN12FUFW*FWK*FUG
c
C
C     STORE SECOND FUNCTIONAL DERIVATIVE
C
            EX (IIPT,KRARA) = EX (IIPT,KRARA)
     $          + SCALE*(ERR*FMN12+Two*ER*FMN12R+E*FMN12RR)
            EX (IIPT,KRAGA) = EX (IIPT,KRAGA)
     $          + SCALE*(ER*FMN12G+E*FMN12RG)/(Two*GamX)
            EX (IIPT,KRATA) = EX (IIPT,KRATA)
     $          + SCALE*(ER*FMN12K+E*FMN12RK)*F2
            EX (IIPT,KGAGA) = EX (IIPT,KGAGA)
     $          + SCALE*(E/(Four*GamX2*GamX))*(GamX*FMN12GG - FMN12G)
            EX (IIPT,KGATA) = EX (IIPT,KGATA)
     $          + SCALE*E*FMN12GK/(Two*GamX)*F2
            EX (IIPT,KTATA) = EX (IIPT,KTATA)
     $          + SCALE*(E*FMN12KK)*F4


      ENDIF

C ---Beta component---

      IF((RHOB.GT.TINY).AND.(TAUB.GT.TINY)) THEN
          pX    = RHOB
          GamX2 = GRDBB
          S2    = GamX2*pX**(-F8o3)
          E = Ax*pX**F4o3
          Tau    = TAUB
          TauUEG =F3o5*((F6*PI*PI)**F2o3)*(pX**F5o3)
          Tsig =TauUEG/Tau
          FV = ome*pX**F1o3/(One+ome*pX**F1o3)
          FU = G*S2/(One+G*S2)
          FW = (Tsig - One)/(Tsig + One)

          FMN12 = CC000 + FW*CC001 + FW**2*CC002 + FW**3*CC003 +
     $          FW**4*CC004 + FW**5*CC005 +
     $          FU*CC010 + FW*FU*CC011 + FW**2*FU*CC012 +
     $          FW**3*FU*CC013 +
     $          FW**4*FU*CC014 + FU**2*CC020 +
     $          FW*FU**2*CC021 + FW**2*FU**2*CC022 +
     $          FW**3*FU**2*CC023 +
     $          FU**3*CC030 + FW*FU**3*CC031 +
     $          FW**2*FU**3*CC032 + FV*CC100 + FV*FW*CC101 +
     $          FV*FW**2*CC102 + FV*FW**3*CC103 +
     $          FV*FW**4*CC104 + FV*FU*CC110 + FV*FW*FU*CC111 +
     $          FV*FW**2*FU*CC112 + FV*FW**3*FU*CC113 +
     $          FV*FU**2*CC120 + FV*FW*FU**2*CC121 +
     $          FV*FW**2*FU**2*CC122 +
     $          FV**2*CC200 + FV**2*FW*CC201 +
     $          FV**2*FW**2*CC202 + FV**2*FW**3*CC203 +
     $          FV**2*FU*CC210 +
     $          FV**2*FW*FU*CC211 +
     $          FV**2*FW**2*FU*CC212 + FV**3*CC300 +
     $          FV**3*FW*CC301 +
     $          FV**3*FW**2*CC302

      ENERGY = E * FMN12
C energies are computed
        EXF          = EXF          + SCALE*ENERGY*XYZWGT(IPT)
        EX0(IIPT)    = EX0(IIPT)    + SCALE*ENERGY

C First Derivatives
          ER = F4o3*E/pX
          S    = Sqrt(S2)
          GamX = Sqrt(GamX2)
          SR   = -F4o3*S/pX
          SG   = S/GamX
          TR  = ((F6*PI*PI)**F2o3)*(pX**F2o3)/Tau
          TK  = -TauUEG/Tau**2
          FUS = Two*G*S/((One+G*S*S)**2)
          FVR =  (ome/(Three*pX**F2o3))*
     $           (One+ome*pX**F1o3)**(-Two)
          FWT = Two/((One + Tsig)**2)
          FUR = FUS*SR
          FUG = FUS*SG
          FWR = FWT*TR
          FWK = FWT*TK


          FMN12FV = CC100+FW*CC101+FW**2*CC102+FW**3*CC103+
     $              FW**4*CC104+FU*CC110+
     $              FW*FU*CC111+FW**2*FU*CC112+FW**3*FU*CC113+
     $              FU**2*CC120+FW*FU**2*CC121+
     $              FW**2*FU**2*CC122+Two*FV*CC200+Two*FV*FW*CC201+
     $              Two*FV*FW**2*CC202+Two*FV*FW**3*CC203+
     $              Two*FV*FU*CC210+Two*FV*FW*FU*CC211+
     $              Two*FV*FW**2*FU*CC212+
     $              Three*FV**2*CC300+Three*FV**2*FW*CC301+
     $              Three*FV**2*FW**2*CC302
          FMN12FU = CC010+FW*CC011+FW**2*CC012+FW**3*CC013+
     $              FW**4*CC014+Two*FU*CC020+
     $              Two*FW*FU*CC021+Two*FW**2*FU*CC022+
     $              Two*FW**3*FU*CC023+
     $              Three*FU**2*CC030+Three*FW*FU**2*CC031+
     $              Three*FW**2*FU**2*CC032+FV*CC110+FV*FW*CC111+
     $              FV*FW**2*CC112+FV*FW**3*CC113+
     $              Two*FV*FU*CC120+Two*FV*FW*FU*CC121+
     $              Two*FV*FW**2*FU*CC122+
     $              FV**2*CC210+FV**2*FW*CC211+
     $              FV**2*FW**2*CC212
          FMN12FW = CC001+Two*FW*CC002+Three*FW**2*CC003+
     $              Four*FW**3*CC004+Five*FW**4*CC005+FU*CC011+
     $              Two*FW*FU*CC012+Three*FW**2*FU*CC013+
     $              Four*FW**3*FU*CC014+
     $              FU**2*CC021+Two*FW*FU**2*CC022+
     $              Three*FW**2*FU**2*CC023+FU**3*CC031+
     $              Two*FW*FU**3*CC032+
     $              FV*CC101+Two*FV*FW*CC102+
     $              Three*FV*FW**2*CC103+Four*FV*FW**3*CC104+
     $              FV*FU*CC111+
     $              Two*FV*FW*FU*CC112+Three*FV*FW**2*FU*CC113+
     $              FV*FU**2*CC121+Two*FV*FW*FU**2*CC122+FV**2*CC201+
     $              Two*FV**2*FW*CC202+
     $              Three*FV**2*FW**2*CC203+FV**2*FU*CC211+
     $              Two*FV**2*FW*FU*CC212+FV**3*CC301+
     $              Two*FV**3*FW*CC302
c
          FMN12R  = FMN12FV*FVR+FMN12FU*FUR+FMN12FW*FWR
          FMN12G  = FMN12FU*FUG
          FMN12K  = FMN12FW*FWK
C
C     STORE FIRST FUNCTIONAL DERIVATIVE
C
            EX (IIPT,KRB) = EX (IIPT,KRB) + SCALE*(ER*FMN12+E*FMN12R)
            EX (IIPT,KGB) = EX (IIPT,KGB) + SCALE*E*FMN12G/(Two*GamX)
            EX (IIPT,KTB) = EX (IIPT,KTB) + SCALE*E*FMN12K*F2
C Second Derivatives
          ERR = ER/(Three*pX)
          SRR = F28o9*S/(pX*pX)
          SRG = SR/GamX
          TRR = F2o3*TR/pX
          TRK = -TR/tau
          TKK = -Two*TK/tau
          FWTT = -Four/(Tsig+One)**3
          FUSS = Two*G*(One-Three*G*S*S)/((One+G*S*S)**3)
          FVRR = (-Two*ome*(One + Two*ome*pX**F1o3))/
     $           (F9*(One + ome*pX**F1o3)**3*pX**F5o3)
          FURR = FUSS*SR*SR +
     $           FUS*SRR
          FUGG = FUSS*SG*SG
          FURG = FUSS*SR*SG +
     $           FUS*SRG
          FWRR = FWTT*TR*TR +
     $           FWT*TRR
          FWRK = FWTT*TK*TR +
     $           FWT*TRK
          FWKK = FWTT*TK*TK +
     $           FWT*TKK
c
          FMN12FVFV = Two*CC200+Two*FW*CC201+Two*FW**2*CC202+
     $                Two*FW**3*CC203+Two*FU*CC210+Two*FU*FW*CC211+
     $                Two*FU*FW**2*CC212+Six*FV*CC300+
     $                Six*FV*FW*CC301+Six*FV*FW**2*CC302
          FMN12FUFU = Two*CC020+Two*FW*CC021+Two*FW**2*CC022+
     $                Two*FW**3*CC023+Six*FU*CC030+Six*FU*FW*CC031+
     $                Six*FU*FW**2*CC032+Two*FV*CC120+Two*FV*FW*CC121+
     $                Two*FV*FW**2*CC122
          FMN12FWFW = Two*CC002+Six*FW*CC003+12*FW**2*CC004+
     $                20.0D+00*FW**3*CC005+Two*FU*CC012+
     $                Six*FU*FW*CC013+
     $                12.0D+00*FU*FW**2*CC014+Two*FU**2*CC022+
     $                Six*FU**2*FW*CC023+Two*FU**3*CC032+Two*FV*CC102+
     $                Six*FV*FW*CC103+12.0D+00*FV*FW**2*CC104+
     $                Two*FU*FV*CC112+Six*FU*FV*FW*CC113+
     $                Two*FU**2*FV*CC122+Two*FV**2*CC202+
     $                Six*FV**2*FW*CC203+Two*FU*FV**2*CC212+
     $                Two*FV**3*CC302
          FMN12FVFU = CC110+FW*CC111+FW**2*CC112+FW**3*CC113+
     $                Two*FU*CC120+Two*FU*FW*CC121+
     $                Two*FU*FW**2*CC122+Two*FV*CC210+
     $                Two*FV*FW*CC211+Two*FV*FW**2*CC212
          FMN12FVFW = CC101+Two*FW*CC102+Three*FW**2*CC103+
     $                Four*FW**3*CC104+FU*CC111+Two*FU*FW*CC112+
     $                Three*FU*FW**2*CC113+FU**2*CC121+
     $                Two*FU**2*FW*CC122+Two*FV*CC201+Four*FV*FW*CC202+
     $                Six*FV*FW**2*CC203+Two*FU*FV*CC211+
     $                Four*FU*FV*FW*CC212+Three*FV**2*CC301+
     $                Six*FV**2*FW*CC302
          FMN12FUFW = CC011+Two*FW*CC012+Three*FW**2*CC013+
     $                Four*FW**3*CC014+Two*FU*CC021+Four*FU*FW*CC022+
     $                Six*FU*FW**2*CC023+Three*FU**2*CC031+
     $                Six*FU**2*FW*CC032+FV*CC111+Two*FV*FW*CC112+
     $                Three*FV*FW**2*CC113+Two*FU*FV*CC121+
     $                Four*FU*FV*FW*CC122+FV**2*CC211+
     $                Two*FV**2*FW*CC212
C

          FMN12RR = FMN12FVFV*FVR*FVR +
     $              FMN12FVFU*FVR*FUR +
     $              FMN12FVFW*FVR*FWR +
     $              FMN12FV*FVRR +
     $              FMN12FVFU*FUR*FVR +
     $              FMN12FUFU*FUR*FUR +
     $              FMN12FUFW*FUR*FWR +
     $              FMN12FU*FURR +
     $              FMN12FVFW*FWR*FVR +
     $              FMN12FUFW*FWR*FUR +
     $              FMN12FWFW*FWR*FWR +
     $              FMN12FW*FWRR

          FMN12RG = FMN12FVFU*FUG*FVR +
     $              FMN12FUFU*FUG*FUR +
     $              FMN12FUFW*FUG*FWR +
     $              FMN12FU*FURG
          FMN12RK = FMN12FVFW*FWK*FVR +
     $              FMN12FUFW*FWK*FUR +
     $              FMN12FWFW*FWK*FWR +
     $              FMN12FW*FWRK
          FMN12GG = FMN12FUFU*FUG*FUG+
     $              FMN12FU*FUGG
          FMN12KK = FMN12FWFW*FWK*FWK+
     $              FMN12FW*FWKK
          FMN12GK = FMN12FUFW*FWK*FUG
c
C
C     STORE SECOND FUNCTIONAL DERIVATIVE
C
            EX (IIPT,KRBRB) = EX (IIPT,KRBRB)
     $          + SCALE*(ERR*FMN12+Two*ER*FMN12R+E*FMN12RR)
            EX (IIPT,KRBGB) = EX (IIPT,KRBGB)
     $          + SCALE*(ER*FMN12G+E*FMN12RG)/(Two*GamX)
            EX (IIPT,KRBTB) = EX (IIPT,KRBTB)
     $          + SCALE*(ER*FMN12K+E*FMN12RK)*F2
            EX (IIPT,KGBGB) = EX (IIPT,KGBGB)
     $      + SCALE*(E/(Four*GamX2*GamX))*(GamX*FMN12GG - FMN12G)
            EX (IIPT,KGBTB) = EX (IIPT,KGBTB)
     $          + SCALE*(E*FMN12GK)/(Two*GamX)*F2
            EX (IIPT,KTBTB) = EX (IIPT,KTBTB)
     $          + SCALE*(E*FMN12KK)*F4



      ENDIF
      ENDIF

      ENDDO
      RETURN
      END

C*MODULE TDDXCE  *DECK H0TERM
      Subroutine H0TERM(MaxDrF,Rho,GRho,rS,Zeta,PotLC,dLdS,dLdZ,dLdSS,
     $  dLdSZ,dLdZZ,dLdSSS,dLdSSZ,dLdSZZ,dLdZZZ,H0,d0dR,d0dG,d0dZ,d0dRR,
     $  d0dRG,d0dGG,d0dRZ,d0dGZ,d0dZZ,d0dRRR,d0dRRG,d0dRGG,d0dGGG,
     $  d0dRRZ,d0dRGZ,d0dGGZ,d0dRZZ,d0dGZZ,d0dZZZ)
      Implicit Real*8(A-H,O-Z)
      Save XNu,CC0,F1,F2,F3,F4,F6,F7,F8,F9,F12,F18,F27,F28,F36,F60
      Data XNu/15.75592D0/,CC0/0.004235D0/,F1/1.0d0/,F2/2.0d0/,
     $  F3/3.0d0/,F4/4.0d0/,F6/6.0d0/,F7/7.0d0/,F8/8.0d0/,F9/9.0d0/,
     $  F12/12.0d0/,F18/18.0d0/,F27/27.0d0/,F28/28.0d0/,F36/36.0d0/,
     $  F60/60.0d0/
      F13 = F1 / F3
      F16 = F13 / F2
      F23 = F13 + F13
      F43 = F23 + F23
      F76 = F7 / F6
      Beta = XNu*CC0
      Pi = F4*ATan(F1)
      GPBE = (F1-Log(F2))/Pi**2
      Alpha =  Beta**2/(F2*GPBE)
      Delta = F2*Alpha/Beta
      BOvD = Beta/Delta
      FKFac = (F3*Pi*Pi)**F13
      SKFac = Sqrt(FKFac/Pi) * F2
C
      Rho16 = Rho**F16
      Rho76 = Rho*Rho16
      Call EvGZet(Zeta,G,dGdZ,dGdZZ,dGdZZZ)
      T = GRho/(F2*SKFac*Rho76*G)
      G2 = G*G
      G3 = G*G2
      PON = -Delta*PotLC/(G3*Beta)
      ExpPON = Exp(PON)
      B = Delta/(ExpPON-F1)
      B2 = B*B
      T2 = T*T
      T4 = T2*T2
      Q4 = F1 + B*T2
      Q5 = Q4 + B2*T4
      A  = F1 + Delta*Q4*T2/Q5
      H0X = Log(A)
      H0 = G3*BOvD*H0X
      If(MaxDrF.eq.0) Return
C
C     1st derivatives.
C
      drSdR = (-(F13*rS/Rho))
      dPodR = (dLdS*drSdR*PON/PotLC)
      dPodZ = (-F3*dGdZ*PON/G + dLdZ*PON/PotLC)
      dBdR = (B*dPodR*ExpPON/(F1 - ExpPON))
      dBdZ = (B*dPodZ*ExpPON/(F1 - ExpPON))
      dT2dR = (-F2*F76*T2/Rho)
      dT2dG = (F2*T2/GRho)
      dT2dZ = (-F2*dGdZ*T2/G)
      dQ4dR = (B*dT2dR + dBdR*T2)
      dQ4dG = (B*dT2dG)
      dQ4dZ = (B*dT2dZ + dBdZ*T2)
      dQ5dR = (dQ4dR + F2*B**2*dT2dR*T2 + F2*B*dBdR*T2**2)
      dQ5dG = (dQ4dG + F2*B**2*dT2dG*T2)
      dQ5dZ = (dQ4dZ + F2*B**2*dT2dZ*T2 + F2*B*dBdZ*T2**2)
      dAdR = Delta*(dT2dR*Q4 + T2*dQ4dR - T2*Q4*dQ5dR/Q5 )/Q5
      dAdG = Delta*(dT2dG*Q4 + T2*dQ4dG - T2*Q4*dQ5dG/Q5 )/Q5
      dAdZ = Delta*(dT2dZ*Q4 + T2*dQ4dZ - T2*Q4*dQ5dZ/Q5 )/Q5
      d0dR = (BOvD*dAdR*G**3/A)
      d0dG = (BOvD*dAdG*G**3/A)
      d0dZ = (BOvD*G**2*(dAdZ*G + F3*A*dGdZ*Log(A))/A)
      If(MaxDrF.le.1) Return
C
C     2nd derivatives.
C
      drSdRR = (-F43*drSdR/Rho)
      dPodRR = ((dLdSS*drSdR**2 + dLdS*drSdRR)*PON/PotLC)
      dPodRZ = ((-F3*dGdZ*dLdS*drSdR + dLdSZ*drSdR*G)*PON/(G*PotLC))
      dPodZZ = (PON*(-F6*dGdZ*dLdZ*G + dLdZZ*G**2 + F12*dGdZ**2*PotLC -
     $  F3*dGdZZ*G*PotLC)/(G**2*PotLC))
      dBdRR = (B*ExpPON*(dPodR**2 + dPodRR + dPodR**2*ExpPON -
     $  dPodRR*ExpPON)/(-F1 + ExpPON)**2)
      dBdRZ = (B*ExpPON*(dPodRZ + dPodR*dPodZ - dPodRZ*ExpPON +
     $  dPodR*dPodZ*ExpPON)/(-F1 + ExpPON)**2)
      dBdZZ = (B*ExpPON*(dPodZ**2 + dPodZZ + dPodZ**2*ExpPON -
     $  dPodZZ*ExpPON)/(-F1 + ExpPON)**2)
      dT2dRR = (F2*F76*(F1 + F2*F76)*T2/Rho**2)
      dT2dRG = (-F4*F76*T2/(GRho*Rho))
      dT2dRZ = (F4*dGdZ*F76*T2/(G*Rho))
      dT2dGG = (F2*T2/GRho**2)
      dT2dGZ = (-F4*dGdZ*T2/(G*GRho))
      dT2dZZ = (F2*(F3*dGdZ**2 - dGdZZ*G)*T2/G**2)
      dQ4dRR = (F2*dBdR*dT2dR + B*dT2dRR + dBdRR*T2)
      dQ4dRG = (dBdR*dT2dG + B*dT2dRG)
      dQ4dRZ = (dBdZ*dT2dR + B*dT2dRZ + dBdR*dT2dZ + dBdRZ*T2)
      dQ4dGG = (B*dT2dGG)
      dQ4dGZ = (dBdZ*dT2dG + B*dT2dGZ)
      dQ4dZZ = (F2*dBdZ*dT2dZ + B*dT2dZZ + dBdZZ*T2)
      dQ5dRR = (dQ4dRR + F2*B**2*dT2dR**2 + F8*B*dBdR*dT2dR*T2 +
     $  F2*B**2*dT2dRR*T2 + F2*dBdR**2*T2**2 + F2*B*dBdRR*T2**2)
      dQ5dRG = (dQ4dRG + F2*B**2*dT2dG*dT2dR + F4*B*dBdR*dT2dG*T2 +
     $  F2*B**2*dT2dRG*T2)
      dQ5dRZ = (dQ4dRZ + F2*B**2*dT2dR*dT2dZ + F4*B*dBdZ*dT2dR*T2 +
     $  F2*B**2*dT2dRZ*T2 + F4*B*dBdR*dT2dZ*T2 + F2*B*dBdRZ*T2**2 +
     $  F2*dBdR*dBdZ*T2**2)
      dQ5dGG = (dQ4dGG + F2*B**2*dT2dG**2 + F2*B**2*dT2dGG*T2)
      dQ5dGZ = (dQ4dGZ + F2*B**2*dT2dG*dT2dZ + F4*B*dBdZ*dT2dG*T2 +
     $  F2*B**2*dT2dGZ*T2)
      dQ5dZZ = (dQ4dZZ + F2*B**2*dT2dZ**2 + F8*B*dBdZ*dT2dZ*T2 +
     $  F2*B**2*dT2dZZ*T2 + F2*dBdZ**2*T2**2 + F2*B*dBdZZ*T2**2)
      dAdRR = Delta*( dT2dRR*Q4 + F2*dT2dR*dQ4dR + T2*dQ4dRR - (
     $  F2*dT2dR*Q4*dQ5dR + F2*T2*dQ4dR*dQ5dR + T2*Q4*dQ5dRR -
     $  F2*T2*Q4*dQ5dR*dQ5dR/Q5 ) / Q5 ) / Q5
      dAdRG = Delta*( dT2dRG*Q4 + dT2dR*dQ4dG + dT2dG*dQ4dR + T2*dQ4dRG
     $  - ( dT2dR*Q4*dQ5dG + dT2dG*Q4*dQ5dR + T2*dQ4dR*dQ5dG +
     $  T2*dQ4dG*dQ5dR + T2*Q4*dQ5dRG - F2*T2*Q4*dQ5dR*dQ5dG/Q5 ) / Q5
     $  ) / Q5
      dAdRZ = Delta*( dT2dRZ*Q4 + dT2dR*dQ4dZ + dT2dZ*dQ4dR + T2*dQ4dRZ
     $  - ( dT2dR*Q4*dQ5dZ + dT2dZ*Q4*dQ5dR + T2*dQ4dR*dQ5dZ +
     $  T2*dQ4dZ*dQ5dR + T2*Q4*dQ5dRZ - F2*T2*Q4*dQ5dR*dQ5dZ/Q5 ) / Q5
     $  ) / Q5
      dAdGG = Delta*( dT2dGG*Q4 + F2*dT2dG*dQ4dG + T2*dQ4dGG - (
     $  F2*dT2dG*Q4*dQ5dG + F2*T2*dQ4dG*dQ5dG + T2*Q4*dQ5dGG -
     $  F2*T2*Q4*dQ5dG*dQ5dG/Q5 ) / Q5 ) / Q5
      dAdGZ = Delta*( dT2dGZ*Q4 + dT2dG*dQ4dZ + dT2dZ*dQ4dG + T2*dQ4dGZ
     $  - ( dT2dG*Q4*dQ5dZ + dT2dZ*Q4*dQ5dG + T2*dQ4dG*dQ5dZ +
     $  T2*dQ4dZ*dQ5dG + T2*Q4*dQ5dGZ - F2*T2*Q4*dQ5dG*dQ5dZ/Q5 ) / Q5
     $  ) / Q5
      dAdZZ = Delta*( dT2dZZ*Q4 + F2*dT2dZ*dQ4dZ + T2*dQ4dZZ - (
     $  F2*dT2dZ*Q4*dQ5dZ + F2*T2*dQ4dZ*dQ5dZ + T2*Q4*dQ5dZZ -
     $  F2*T2*Q4*dQ5dZ*dQ5dZ/Q5 ) / Q5 ) / Q5
      d0dRR = (BOvD*(-dAdR**2 + A*dAdRR)*G**3/A**2)
      d0dRG = (BOvD*(-(dAdG*dAdR) + A*dAdRG)*G**3/A**2)
      d0dRZ = (BOvD*G**2*(F3*A*dAdR*dGdZ + A*dAdRZ*G -
     $  dAdR*dAdZ*G)/A**2)
      d0dGG = (BOvD*(-dAdG**2 + A*dAdGG)*G**3/A**2)
      d0dGZ = (BOvD*G**2*(F3*A*dAdG*dGdZ + A*dAdGZ*G -
     $  dAdG*dAdZ*G)/A**2)
      d0dZZ = (BOvD*G*(F6*A*dAdZ*dGdZ*G - dAdZ**2*G**2 +
     $  A*dAdZZ*G**2 + F6*A**2*dGdZ**2*Log(A) +
     $  F3*A**2*dGdZZ*G*Log(A))/A**2)
      If(MaxDrF.le.2) Return
C
C     3rd derivatives.
C
      rSdRRR = -F28*rS/(F27*Rho**3)
      T2dRRR = -F7/F3*( dT2dRR/Rho - F2*dT2dR/Rho**2 + F2*T2/Rho**3 )
      T2dRRG = -F7/F3*( dT2dRG/Rho - dT2dG/Rho**2 )
      T2dRRZ = -F7/F3*( dT2dRZ/Rho - dT2dZ/Rho**2 )
      T2dRGG = +F2*( dT2dRG/GRho - dT2dR/GRho**2 )
      T2dRGZ = -F7/F3*dT2dGZ/Rho
      T2dRZZ = -F2*( dT2dRZ*( dGdZ/G ) + dT2dR*( dGdZZ/G -
     $  dGdZ**2/G**2 ) )
      T2dGGG = +F2*( dT2dGG/GRho - F2*dT2dG/GRho**2 + F2*T2/GRho**3 )
      T2dGGZ = +F2*( dT2dGZ/GRho - dT2dZ/GRho**2 )
      T2dGZZ = -F2*( dT2dGZ*( dGdZ/G ) + dT2dG*( dGdZZ/G -
     $  dGdZ**2/G**2 ) )
      T2dZZZ = -F2*( dT2dZZ*( dGdZ/G ) + F2*dT2dZ*( dGdZZ/G -
     $  dGdZ**2/G**2 ) + T2*( dGdZZZ/G - F3*dGdZ*dGdZZ/G**2 +
     $  F2*dGdZ**3/G**3 ) )
      PodRRR = -Delta/Beta*( dLdSSS*drSdR**3 + F3*dLdSS*drSdRR*drSdR
     $  + dLdS*rSdRRR )/G**3
      PodRRZ = -Delta/Beta*( G*dLdSSZ*drSdR**2 + G*dLdSZ*drSdRR -
     $  F3*dLdSS*drSdR**2*dGdZ - F3*dLdS*drSdRR*dGdZ )/G**4
      PodRZZ = -Delta/Beta*( G**2*dLdSZZ*drSdR -
     $  F6*G*dLdSZ*drSdR*dGdZ + F12*dLdS*drSdR*dGdZ**2 -
     $  F3*G*dLdS*drSdR*dGdZZ )/G**5
      PodZZZ = -Delta/Beta*( dLdZZZ/G**3 - F9*dLdZZ*dGdZ/G**4 -
     $  F9*dLdZ*dGdZZ/G**4 - F3*PotLC*dGdZZZ/G**4 -
     $  F60*PotLC*dGdZ**3/G**6 + F36*dLdZ*dGdZ**2/G**5 +
     $  F36*PotLC*dGdZZ*dGdZ/G**5 )
      dBdRRR = (F1-ExpPON/(ExpPON-F1))*dBdRR*dPodR -
     $  ExpPON/(ExpPON-F1)*( dPodR*dBdRR + F2*dPodRR*dBdR +
     $  B*PodRRR + (F1-ExpPON/(ExpPON-F1))*( dPodR**2*dBdR -
     $  B*ExpPON/(ExpPON-F1)*dPodR**3 + F2*B*dPodRR*dPodR ) )
      dBdRRZ = (F1-ExpPON/(ExpPON-F1))*dBdRR*dPodZ -
     $  ExpPON/(ExpPON-F1)*( dPodRZ*dBdR + dPodR*dBdRZ + dPodRR*dBdZ +
     $  B*PodRRZ + (F1-ExpPON/(ExpPON-F1))*( dPodR**2*dBdZ -
     $  B*ExpPON/(ExpPON-F1)*dPodR**2*dPodZ + F2*B*dPodRZ*dPodR ) )
      dBdRZZ = (F1-ExpPON/(ExpPON-F1))*dBdRZ*dPodZ -
     $  ExpPON/(ExpPON-F1)*( dPodR*dBdZZ + F2*dPodRZ*dBdZ +
     $  B*PodRZZ + (F1-ExpPON/(ExpPON-F1))*( dPodR*dPodZ*dBdZ -
     $  B*ExpPON/(ExpPON-F1)*dPodR*dPodZ**2 + B*( dPodRZ*dPodZ +
     $  dPodR*dPodZZ ) ) )
      dBdZZZ = (F1-ExpPON/(ExpPON-F1))*dBdZZ*dPodZ -
     $  ExpPON/(ExpPON-F1)*( dPodZ*dBdZZ + F2*dPodZZ*dBdZ +
     $  B*PodZZZ + (F1-ExpPON/(ExpPON-F1))*( dPodZ**2*dBdZ -
     $  B*ExpPON/(ExpPON-F1)*dPodZ**3 + F2*B*dPodZZ*dPodZ ) )
      Q4dRRR = dBdRRR*T2 + F3*dBdRR*dT2dR + F3*dBdR*dT2dRR +
     $  B*T2dRRR
      Q4dRRG = dBdRR*dT2dG + F2*dBdR*dT2dRG + B*T2dRRG
      Q4dRRZ = dBdRRZ*T2 + dBdRR*dT2dZ + F2*dBdRZ*dT2dR +
     $  F2*dBdR*dT2dRZ + dBdZ*dT2dRR + B*T2dRRZ
      Q4dRGG = dBdR*dT2dGG + B*T2dRGG
      Q4dRGZ = dBdRZ*dT2dG + dBdR*dT2dGZ + dBdZ*dT2dRG + B*T2dRGZ
      Q4dRZZ = dBdRZZ*T2 + dBdZZ*dT2dR + F2*dBdRZ*dT2dZ +
     $  F2*dBdZ*dT2dRZ + dBdR*dT2dZZ + B*T2dRZZ
      Q4dGGG = B*T2dGGG
      Q4dGGZ = dBdZ*dT2dGG + B*T2dGGZ
      Q4dGZZ = dBdZZ*dT2dG + F2*dBdZ*dT2dGZ + B*T2dGZZ
      Q4dZZZ = dBdZZZ*T2 + F3*dBdZZ*dT2dZ + F3*dBdZ*dT2dZZ +
     $  B*T2dZZZ
      Q5dRRR = (F1+F2*B*T2)*Q4dRRR + F6*dQ4dRR*dQ4dR
      Q5dRRG = (F1+F2*B*T2)*Q4dRRG + F2*dQ4dRR*dQ4dG +
     $  F4*dQ4dRG*dQ4dR
      Q5dRRZ = (F1+F2*B*T2)*Q4dRRZ + F2*dQ4dRR*dQ4dZ +
     $  F4*dQ4dRZ*dQ4dR
      Q5dRGG = (F1+F2*B*T2)*Q4dRGG + F4*dQ4dRG*dQ4dG +
     $  F2*dQ4dGG*dQ4dR
      Q5dRGZ = (F1+F2*B*T2)*Q4dRGZ + F2*(dQ4dRG*dQ4dZ +
     $  dQ4dRZ*dQ4dG + dQ4dGZ*dQ4dR )
      Q5dRZZ = (F1+F2*B*T2)*Q4dRZZ + F4*dQ4dRZ*dQ4dZ +
     $  F2*dQ4dZZ*dQ4dR
      Q5dGGG = (F1+F2*B*T2)*Q4dGGG + F6*dQ4dGG*dQ4dG
      Q5dGGZ = (F1+F2*B*T2)*Q4dGGZ + F2*dQ4dGG*dQ4dZ +
     $  F4*dQ4dGZ*dQ4dG
      Q5dGZZ = (F1+F2*B*T2)*Q4dGZZ + F4*dQ4dGZ*dQ4dZ +
     $  F2*dQ4dZZ*dQ4dG
      Q5dZZZ = (F1+F2*B*T2)*Q4dZZZ + F6*dQ4dZZ*dQ4dZ
      dAdRRR = Delta*( T2dRRR*Q4 + F3*dT2dRR*dQ4dR + F3*dT2dR*dQ4dRR +
     $  T2*Q4dRRR - ( F3*dT2dRR*Q4*dQ5dR + F6*dT2dR*dQ4dR*dQ5dR +
     $  F3*dT2dR*Q4*dQ5dRR + F3*T2*dQ4dRR*dQ5dR + F3*T2*dQ4dR*dQ5dRR  +
     $  T2*Q4*Q5dRRR - ( F6*dT2dR*Q4*dQ5dR**2 + F6*T2*dQ4dR*dQ5dR**2 +
     $  F6*T2*Q4*dQ5dRR*dQ5dR - F6*T2*Q4*dQ5dR**3/Q5 )/Q5 )/Q5 )/Q5
      dAdRRG = Delta*( T2dRRG*Q4 + dT2dRR*dQ4dG + F2*dT2dRG*dQ4dR +
     $  F2*dT2dR*dQ4dRG + dT2dG*dQ4dRR + T2*Q4dRRG - ( dT2dRR*Q4*dQ5dG
     $  + dT2dRG*Q4*dQ5dR + dT2dRG*Q4*dQ5dR + F2*dT2dR*dQ4dR*dQ5dG +
     $  F2*dT2dG*dQ4dR*dQ5dR + F2*dT2dR*dQ4dG*dQ5dR +
     $  F2*dT2dR*Q4*dQ5dRG + dT2dG*Q4*dQ5dRR + T2*dQ4dRR*dQ5dG +
     $  F2*T2*dQ4dRG*dQ5dR + F2*T2*dQ4dR*dQ5dRG + T2*dQ4dG*dQ5dRR +
     $  T2*Q4*Q5dRRG - ( F4*dT2dR*Q4*dQ5dR*dQ5dG + F2*dT2dG*Q4*dQ5dR**2
     $  + F4*T2*dQ4dR*dQ5dR*dQ5dG + F2*T2*dQ4dG*dQ5dR**2 +
     $  F2*T2*Q4*dQ5dRR*dQ5dG + F4*T2*Q4*dQ5dRG*dQ5dR -
     $  F6*T2*Q4*dQ5dR**2*dQ5dG/Q5 )/Q5 )/Q5 )/Q5
      dAdRRZ = Delta*( T2dRRZ*Q4 + dT2dRR*dQ4dZ + F2*dT2dRZ*dQ4dR +
     $  F2*dT2dR*dQ4dRZ + dT2dZ*dQ4dRR + T2*Q4dRRZ - ( dT2dRR*Q4*dQ5dZ
     $  + dT2dRZ*Q4*dQ5dR + dT2dRZ*Q4*dQ5dR + F2*dT2dR*dQ4dR*dQ5dZ +
     $  F2*dT2dZ*dQ4dR*dQ5dR + F2*dT2dR*dQ4dZ*dQ5dR +
     $  F2*dT2dR*Q4*dQ5dRZ + dT2dZ*Q4*dQ5dRR + T2*dQ4dRR*dQ5dZ +
     $  F2*T2*dQ4dRZ*dQ5dR + F2*T2*dQ4dR*dQ5dRZ + T2*dQ4dZ*dQ5dRR +
     $  T2*Q4*Q5dRRZ - ( F4*dT2dR*Q4*dQ5dR*dQ5dZ + F2*dT2dZ*Q4*dQ5dR**2
     $  + F4*T2*dQ4dR*dQ5dR*dQ5dZ + F2*T2*dQ4dZ*dQ5dR**2 +
     $  F2*T2*Q4*dQ5dRR*dQ5dZ + F4*T2*Q4*dQ5dRZ*dQ5dR -
     $  F6*T2*Q4*dQ5dR**2*dQ5dZ/Q5 )/Q5 )/Q5 )/Q5
      dAdRGG = Delta*( T2dRGG*Q4 + F2*dT2dRG*dQ4dG + dT2dGG*dQ4dR +
     $  dT2dR*dQ4dGG + F2*dT2dG*dQ4dRG + T2*Q4dRGG - (
     $  F2*dT2dRG*Q4*dQ5dG + dT2dGG*Q4*dQ5dR + F2*dT2dR*dQ4dG*dQ5dG +
     $  dT2dR*Q4*dQ5dGG + F2*dT2dG*dQ4dR*dQ5dG + F2*dT2dG*dQ4dG*dQ5dR +
     $  F2*dT2dG*Q4*dQ5dRG + F2*T2*dQ4dRG*dQ5dG + T2*dQ4dGG*dQ5dR +
     $  T2*dQ4dR*dQ5dGG + F2*T2*dQ4dG*dQ5dRG + T2*Q4*Q5dRGG - (
     $  F2*dT2dR*Q4*dQ5dG*dQ5dG + F4*dT2dG*Q4*dQ5dR*dQ5dG +
     $  F2*T2*dQ4dR*dQ5dG*dQ5dG + F4*T2*dQ4dG*dQ5dR*dQ5dG +
     $  F4*T2*Q4*dQ5dRG*dQ5dG + F2*T2*Q4*dQ5dR*dQ5dGG -
     $  F6*T2*Q4*dQ5dR*dQ5dG*dQ5dG/Q5 )/Q5 )/Q5 )/Q5
      dAdRGZ = Delta*( T2dRGZ*Q4 + dT2dRG*dQ4dZ + dT2dRZ*dQ4dG +
     $  dT2dGZ*dQ4dR + dT2dR*dQ4dGZ + dT2dG*dQ4dRZ + dT2dZ*dQ4dRG +
     $  T2*Q4dRGZ - ( dT2dRG*Q4*dQ5dZ + dT2dRZ*Q4*dQ5dG +
     $  dT2dGZ*Q4*dQ5dR + dT2dR*dQ4dG*dQ5dZ + dT2dR*dQ4dZ*dQ5dG +
     $  dT2dR*Q4*dQ5dGZ + dT2dG*dQ4dR*dQ5dZ + dT2dG*dQ4dZ*dQ5dR +
     $  dT2dG*Q4*dQ5dRZ + dT2dZ*dQ4dR*dQ5dG + dT2dZ*dQ4dG*dQ5dR +
     $  dT2dZ*Q4*dQ5dRG + T2*dQ4dRG*dQ5dZ + T2*dQ4dRZ*dQ5dG +
     $  T2*dQ4dGZ*dQ5dR + T2*dQ4dR*dQ5dGZ + T2*dQ4dG*dQ5dRZ +
     $  T2*dQ4dZ*dQ5dRG + T2*Q4*Q5dRGZ - ( F2*dT2dR*Q4*dQ5dG*dQ5dZ +
     $  F2*dT2dG*Q4*dQ5dR*dQ5dZ + F2*dT2dZ*Q4*dQ5dR*dQ5dG +
     $  F2*T2*dQ4dR*dQ5dG*dQ5dZ + F2*T2*dQ4dG*dQ5dR*dQ5dZ +
     $  F2*T2*dQ4dZ*dQ5dR*dQ5dG + F2*T2*Q4*dQ5dRG*dQ5dZ +
     $  F2*T2*Q4*dQ5dRZ*dQ5dG + F2*T2*Q4*dQ5dR*dQ5dGZ -
     $  F6*T2*Q4*dQ5dR*dQ5dG*dQ5dZ/Q5 )/Q5 )/Q5 )/Q5
      dAdRZZ = Delta*( T2dRZZ*Q4 + F2*dT2dRZ*dQ4dZ + dT2dZZ*dQ4dR +
     $  dT2dR*dQ4dZZ + F2*dT2dZ*dQ4dRZ + T2*Q4dRZZ - (
     $  F2*dT2dRZ*Q4*dQ5dZ + dT2dZZ*Q4*dQ5dR + F2*dT2dR*dQ4dZ*dQ5dZ +
     $  dT2dR*Q4*dQ5dZZ + F2*dT2dZ*dQ4dR*dQ5dZ + F2*dT2dZ*dQ4dZ*dQ5dR +
     $  F2*dT2dZ*Q4*dQ5dRZ + F2*T2*dQ4dRZ*dQ5dZ + T2*dQ4dZZ*dQ5dR +
     $  T2*dQ4dR*dQ5dZZ + F2*T2*dQ4dZ*dQ5dRZ + T2*Q4*Q5dRZZ - (
     $  F2*dT2dR*Q4*dQ5dZ*dQ5dZ + F4*dT2dZ*Q4*dQ5dR*dQ5dZ +
     $  F2*T2*dQ4dR*dQ5dZ*dQ5dZ + F4*T2*dQ4dZ*dQ5dR*dQ5dZ +
     $  F4*T2*Q4*dQ5dRZ*dQ5dZ + F2*T2*Q4*dQ5dR*dQ5dZZ -
     $  F6*T2*Q4*dQ5dR*dQ5dZ*dQ5dZ/Q5 )/Q5 )/Q5 )/Q5
      dAdGGG = Delta*( T2dGGG*Q4 + F3*dT2dGG*dQ4dG + F3*dT2dG*dQ4dGG +
     $  T2*Q4dGGG - ( F3*dT2dGG*Q4*dQ5dG + F6*dT2dG*dQ4dG*dQ5dG +
     $  F3*dT2dG*Q4*dQ5dGG + F3*T2*dQ4dGG*dQ5dG + F3*T2*dQ4dG*dQ5dGG  +
     $  T2*Q4*Q5dGGG - ( F6*dT2dG*Q4*dQ5dG**2 + F6*T2*dQ4dG*dQ5dG**2 +
     $  F6*T2*Q4*dQ5dGG*dQ5dG - F6*T2*Q4*dQ5dG**3/Q5 )/Q5 )/Q5 )/Q5
      dAdGGZ = Delta*( T2dGGZ*Q4 + dT2dGG*dQ4dZ + F2*dT2dGZ*dQ4dG +
     $  F2*dT2dG*dQ4dGZ + dT2dZ*dQ4dGG + T2*Q4dGGZ - ( dT2dGG*Q4*dQ5dZ
     $  + dT2dGZ*Q4*dQ5dG + dT2dGZ*Q4*dQ5dG + F2*dT2dG*dQ4dG*dQ5dZ +
     $  F2*dT2dZ*dQ4dG*dQ5dG + F2*dT2dG*dQ4dZ*dQ5dG +
     $  F2*dT2dG*Q4*dQ5dGZ + dT2dZ*Q4*dQ5dGG + T2*dQ4dGG*dQ5dZ +
     $  F2*T2*dQ4dGZ*dQ5dG + F2*T2*dQ4dG*dQ5dGZ + T2*dQ4dZ*dQ5dGG +
     $  T2*Q4*Q5dGGZ - ( F4*dT2dG*Q4*dQ5dG*dQ5dZ + F2*dT2dZ*Q4*dQ5dG**2
     $  + F4*T2*dQ4dG*dQ5dG*dQ5dZ + F2*T2*dQ4dZ*dQ5dG**2 +
     $  F2*T2*Q4*dQ5dGG*dQ5dZ + F4*T2*Q4*dQ5dGZ*dQ5dG -
     $  F6*T2*Q4*dQ5dG**2*dQ5dZ/Q5 )/Q5 )/Q5 )/Q5
      dAdGZZ = Delta*( T2dGZZ*Q4 + F2*dT2dGZ*dQ4dZ + dT2dZZ*dQ4dG +
     $  dT2dG*dQ4dZZ + F2*dT2dZ*dQ4dGZ + T2*Q4dGZZ - (
     $  F2*dT2dGZ*Q4*dQ5dZ + dT2dZZ*Q4*dQ5dG + F2*dT2dG*dQ4dZ*dQ5dZ +
     $  dT2dG*Q4*dQ5dZZ + F2*dT2dZ*dQ4dG*dQ5dZ + F2*dT2dZ*dQ4dZ*dQ5dG +
     $  F2*dT2dZ*Q4*dQ5dGZ + F2*T2*dQ4dGZ*dQ5dZ + T2*dQ4dZZ*dQ5dG +
     $  T2*dQ4dG*dQ5dZZ + F2*T2*dQ4dZ*dQ5dGZ + T2*Q4*Q5dGZZ - (
     $  F2*dT2dG*Q4*dQ5dZ*dQ5dZ + F4*dT2dZ*Q4*dQ5dG*dQ5dZ +
     $  F2*T2*dQ4dG*dQ5dZ*dQ5dZ + F4*T2*dQ4dZ*dQ5dG*dQ5dZ +
     $  F4*T2*Q4*dQ5dGZ*dQ5dZ + F2*T2*Q4*dQ5dG*dQ5dZZ -
     $  F6*T2*Q4*dQ5dG*dQ5dZ*dQ5dZ/Q5 )/Q5 )/Q5 )/Q5
      dAdZZZ = Delta*( T2dZZZ*Q4 + F3*dT2dZZ*dQ4dZ + F3*dT2dZ*dQ4dZZ +
     $  T2*Q4dZZZ - ( F3*dT2dZZ*Q4*dQ5dZ + F6*dT2dZ*dQ4dZ*dQ5dZ +
     $  F3*dT2dZ*Q4*dQ5dZZ + F3*T2*dQ4dZZ*dQ5dZ + F3*T2*dQ4dZ*dQ5dZZ  +
     $  T2*Q4*Q5dZZZ - ( F6*dT2dZ*Q4*dQ5dZ**2 + F6*T2*dQ4dZ*dQ5dZ**2 +
     $  F6*T2*Q4*dQ5dZZ*dQ5dZ - F6*T2*Q4*dQ5dZ**3/Q5 )/Q5 )/Q5 )/Q5
      d0dRRR = (BOvD*(F2*dAdR**3 - F3*A*dAdR*dAdRR + A**2*dAdRRR
     $  )*G**3/A**3)
      d0dRRG = (BOvD*(F2*dAdR**2*dAdG - A*(F2*dAdR*dAdRG + dAdRR*dAdG )
     $  + A**2*dAdRRG )*G**3/A**3)
      d0dRRZ = (BOvD*(F3*A/G*dGdZ*( A*dAdRR - dAdR**2 ) + A**2*dAdRRZ +
     $  dAdZ*( F2*dAdR**2 - A*dAdRR ) - F2*A*dAdR*dAdRZ )*G**3/A**3)
      d0dRGG = (BOvD*(F2*dAdR*dAdG**2 - A*(F2*dAdRG*dAdG + dAdR*dAdGG )
     $  + A**2*dAdRGG )*G**3/A**3)
      d0dRGZ = (BOvD*(F3*A/G*dGdZ*( A*dAdRG - dAdR*dAdG ) +
     $  A**2*dAdRGZ + dAdZ*( F2*dAdR*dAdG - A*dAdRG ) - A*( dAdRZ*dAdG +
     $  dAdR*dAdGZ ) )*G**3/A**3)
      d0dRZZ = (BOvD*(F6*A**2*dGdZ**2*dAdR + F3*G*A**2*dGdZZ*dAdR -
     $  F6*G*A*dGdZ*dAdR*dAdZ + F6*G*A**2*dGdZ*dAdRZ +
     $  F2*G**2*dAdR*dAdZ**2 - F2*G**2*A*dAdZ*dAdRZ - G**2*A*dAdR*dAdZZ
     $  + G**2*A**2*dAdRZZ )*G/A**3)
      d0dGGG = (BOvD*(F2*dAdG**3 - F3*A*dAdG*dAdGG + A**2*dAdGGG
     $  )*G**3/A**3)
      d0dGGZ = (BOvD*(F3*A/G*dGdZ*( A*dAdGG - dAdG**2 ) + A**2*dAdGGZ +
     $  dAdZ*( F2*dAdG**2 - A*dAdGG ) - F2*A*dAdG*dAdGZ )*G**3/A**3)
      d0dGZZ = (BOvD*(F6*A**2*dGdZ**2*dAdG + F3*G*A**2*dGdZZ*dAdG -
     $  F6*G*A*dGdZ*dAdG*dAdZ + F6*G*A**2*dGdZ*dAdGZ +
     $  F2*G**2*dAdG*dAdZ**2 - F2*G**2*A*dAdZ*dAdGZ - G**2*A*dAdG*dAdZZ
     $  + G**2*A**2*dAdGZZ )*G/A**3)
      d0dZZZ = (BOvD*(F6*A**3*Log(A)*dGdZ**3 +
     $  F18*G*A**3*Log(A)*dGdZ*dGdZZ +
     $  F18*G*A**2*dGdZ**2*dAdZ - F9*G**2*A*dGdZ*dAdZ**2 +
     $  F2*G**3*dAdZ**3 + F9*G**2*A**2*( dGdZZ*dAdZ + dGdZ*dAdZZ ) -
     $  F3*G**3*A*dAdZ*dAdZZ + F3*G**2*A**3*Log(A)*dGdZZZ +
     $  G**3*A**2*dAdZZZ )/A**3)
      Return
      End

C*MODULE TDDXCE  *Deck EvGZet
      Subroutine EvGZet(Zeta,GZeta,dGZdz,d2GZdz,d3GZdz)
      Implicit Real*8(A-H,O-Z)
C
C     Evaluate G(Zeta) and its derivatives for DFT.
C
      Real*8 Nine
      Save Zero, One, Two, Three, Four, Nine, F27
      Data Zero/0.0d0/, One/1.0d0/, Two/2.0d0/, Three/3.0d0/,
     $  Four/4.0d0/, Nine/9.0d0/, F27/27.0D0/
C
      Small = 1.d-12
      GZeta = Zero
      dGZdz = Zero
      d2GZdz = Zero
      d3GZdz = Zero
      OMZ = One - Zeta
      OPZ = One + Zeta
      OMZ2 = OMZ**2
      OPZ2 = OPZ**2
      F13 = One / Three
      F19 = One / Nine
      F427 = Four / F27
      If(OMZ.gt.Small) then
        OMZ3 = OMZ ** (-F13)
        GZeta = GZeta + OMZ*OMZ3
        dGZdz = dGZdz - OMZ3
        d2GZdz = d2GZdz - OMZ3/OMZ
        d3GZdz = d3GZdz - OMZ3/OMZ2
        endIf
      If(OPZ.gt.Small) then
        OPZ3 = OPZ ** (-F13)
        GZeta = GZeta + OPZ*OPZ3
        dGZdz = dGZdz + OPZ3
        d2GZdz = d2GZdz - OPZ3/OPZ
        d3GZdz = d3GZdz + OPZ3/OPZ2
        endIf
      GZeta = GZeta / Two
      dGZdz = dGZdz * F13
      d2GZdz = d2GZdz * F19
      d3GZdz = d3GZdz * F427
      Return
      End
