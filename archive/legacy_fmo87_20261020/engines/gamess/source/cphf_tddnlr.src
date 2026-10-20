C*MODULE CPUHF   *DECK AOCPCG_DYN
      SUBROUTINE AOCPCG_DYN(WAX,YA,
     *                  VEC,EIG,DTEMP,FTEMP,WRK1,WRK2,
     *                  NOCP,GHONDO,XINTS,DSH,DDIJ,XX,IX,
     *                  NNXYZ,NUNIQ,NFOCK,NROT,NOCC,NVIR,
     *                  L1,L2,L3,NSH2,MAXG,MXG2,FTEMP2,
     *                  HF,YAL,HFL,TRAN,VCOC,DYNDD_LMO,HQ,
     *                  HQL,DYNDD,DYNDQ,DYNDQ_LMO,NA,MCORE,NLOC,
     *                  WAXH2,Q1,Q2,
     *                  RESID,ZRES,PDIR,PRCND,BNORM,BKNUM,BKDEN,WAXA)
       
      use mx_limits, only: mxao

      USE EFP_LOGICAL
      USE MAKEFP_CPHF, only:makefp_print,cpmakefp,makefp_print_default
      use constants, only: zero,one,two,four,half
      IMPLICIT NONE
C
      DIMENSION WAX(NROT,NNXYZ),YA(NROT,NNXYZ),
     *          VEC(L1,L1),EIG(L1),
     *          WRK1(L1,L1),WRK2(L3),NOCP(NNXYZ),
     *          GHONDO(MAXG),XINTS(NSH2),DSH(NSH2),DDIJ(49*MXG2),
     *          XX(*),IX(*)
      DIMENSION YAL(NVIR,3),HF(L2,3),HFL(NVIR,3),TRAN(NLOC,NLOC),
     *          VCOC(3*NA)
      double precision :: HQ(L2,9),HQL(NVIR,9)
      double precision :: DYNDD(9,NDPFREQ),DYNDQ(27,NDPFREQ)
      double precision :: DYNDD_LMO(9,NLOC,NDPFREQ)
      double precision :: DYNDQ_LMO(27,NLOC,NDPFREQ)

      double precision :: WAXH2(NROT,NNXYZ)
      double precision :: Q1(NROT,NNXYZ), Q2(NROT,NNXYZ)
      double precision :: DTEMP(L2,NUNIQ), FTEMP(L2,NUNIQ)
      double precision :: FTEMP2(L2,NUNIQ)
      DIMENSION RESID(NROT,NUNIQ),ZRES(NROT,NUNIQ),PDIR(NROT,NUNIQ)
     *         ,PRCND(NROT),BNORM(NNXYZ),BKNUM(NNXYZ),BKDEN(NNXYZ)

      double precision :: WAXA(NROT)

      LOGICAL GOPARR,DSKWRK,MASWRK,DIRSCF,FDIFF,CVGING
C
      LOGICAL MOIDON,EDCOMP,DIPDCM,QADDCM,DEPRNT,ZDO,
     *        POLDCM,POLANG,POLAPP,KMIDPT,POLDYN,POLAR
C
      LOGICAL second_order_flag
C
      DOUBLE PRECISION DDIJ,DSH,EIG,GHONDO,HF,VEC,WAX,WRK1,WRK2,XINTS
      DOUBLE PRECISION XX,YA,YAL,IX
      DOUBLE PRECISION BKDEN,BKNUM,BNORM,CHFSLV,CPTYPE,HFL,PDIR,PRCND
      DOUBLE PRECISION RESID,TOLCP,TOLSCZV,TRAN,VCOC,ZRES
      DOUBLE PRECISION DPFREQ,OCCUP,ZIJ,ZMO
      DOUBLE PRECISION AK,AKDEN,ERR,ERRSAVE,FREQ,TEMP1,TEMP2,TEMP3,TOL
      DOUBLE PRECISION TEST
C
      INTEGER L1,L2,L3,MAXG,MXG2,NNXYZ,NOCP,NROT,NSH2,NUNIQ,NVIR
      INTEGER IJKO,IP,IR,IW,NA,NCPHF,NDPFREQ,NLOC,NSCZV,NSPLIT,NWDCHF
      INTEGER IBTYP,IDAF,IJKT,IJMO,IODA,IPTIM,MASTER,ME,MOIDNO,MOIJ
      INTEGER NAOTYP,NATFMO,NAV,NBDFG,NBODY,NFG,NLAYER,NMO,NMOAT,NMOIJ
      INTEGER NPROC
      INTEGER IA,ICUT,IPROT,ITOL,NOPK,NORMF,NORMP,NPRINT,NPROT
      INTEGER IFREQ,IOCC,IROT,ITER,IUNIQ,IVIR,IXYZ,MAXCP,MAXCP2,MCORE
      INTEGER MODE,NFOCK,NINT,NOCC,NSCHWZ,NXYZF
      INTEGER NSTEP
C
      COMMON /CHFINF/ TOLCP,TOLSCZV(2),CHFSLV,CPTYPE,NWDCHF,POLAR,
     *                NSPLIT,NCPHF,NSCZV(2)
      COMMON /IOFILE/ IR,IW,IP,IJKO,IJKT,IDAF,NAV,IODA(950)
      COMMON /OPTSCF/ DIRSCF,FDIFF
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /FMOINF/ NFG,NLAYER,NATFMO,NBDFG,NAOTYP,NBODY
C
      PARAMETER (NMO=500)
C
      COMMON /EDCMP / ZIJ(NMO),ZMO(5,NMO),OCCUP(NMO),DPFREQ(50),
     *                MOIDNO(5,NMO),IJMO(2,NMO),MOIJ(NMO),NMOIJ(NMO),
     *                NMOAT(NMO),NDPFREQ,IPROT(5),NPROT,
     *                MOIDON,EDCOMP,DIPDCM,DEPRNT,QADDCM,ZDO,POLDCM,
     *                POLANG,POLAPP,KMIDPT,POLDYN
      COMMON /IJPAIR/ IA(MXAO)
      logical masout
      COMMON /OUTPUT/ NPRINT,ITOL,ICUT,NORMF,NORMP,NOPK

      TOL=5.0D-05
      IF(cpmakefp) TOL=TOLCP

      masout=maswrk.and.makefp_print_default
      !masout=maswrk
      if(masout) WRITE(IW,9000) TOL

C Construct 2*(H2)*P
       call constr_waxh2(WAX,WAXH2,VEC,EIG,NNXYZ,L1,L2,L3
     &         ,DTEMP,FTEMP2,DSH,DDIJ,NUNIQ,NROT,NOCC,NVIR,NOCP
     &         ,WRK1,WRK2,XX,IX,MAXG,MXG2,NSH2,NINT,NSCHWZ,IDAF,IODA)
C
C        SET UP PRECONDITIONER AND INITIAL GUESS OF RESPONSES
C
      CALL DSCAL(NROT*NNXYZ,TWO,YA,1)
C    ------- LOOP OVER FREQUENCIES
      DO Ifreq=1,NDPFREQ
C Freq is an imaginary freq.
         FREQ=DPFREQ(Ifreq)
      IF(MASOUT) WRITE(IW,9031) FREQ
      call iniv_cphf_dyn(nnxyz,nrot,NUNIQ,bknum,bnorm,
     $               bkden,RESID,PDIR,ZRES,NOCP)
C PRECONDITIONER
      IROT = 0
      DO IVIR = NOCC+1,NOCC+NVIR
         DO IOCC= 1,NOCC
            IROT = IROT + 1
            TEMP1 = (EIG(IVIR) - EIG(IOCC))
            TEMP2 = 4.0D+00*(TEMP1**2 + FREQ*FREQ)
            TEMP3 = 4.0D+00*TEMP1**2
            PRCND(IROT) = 1.0D+00/(TEMP2)
            DO 110 IXYZ = 1,NNXYZ
            if(freq.lt.0.1D00) then
               YA(IROT,IXYZ) = TEMP3*YA(IROT,IXYZ)*PRCND(IROT)
            else
               YA(IROT,IXYZ) = -WAXH2(IROT,IXYZ)*PRCND(IROT)
            endif
  110       CONTINUE
         ENDDO !IOCC
      ENDDO !IVIR

      NFOCK = 0
      MODE  = 1
      ITER  = 1
      NXYZF = 0
C
C         THE INITIAL ITERATION IS DONE BEFORE THE MAIN LOOP,
C         CALCULATE INITIAL ITERATION'S RESIDUAL -RESID-
C
      CALL AOCPTD_DYN(VEC,YA,RESID,DTEMP,FTEMP,FTEMP2,EIG,WRK1,WRK2,
     *              XX,IX,GHONDO,XINTS,DSH,DDIJ,NOCP,MODE,
     *              NNXYZ,NUNIQ,NXYZF,NFOCK,NINT,NSCHWZ,
     *              NROT,NOCC,NVIR,L1,L2,L3,NSH2,MAXG,MXG2,
     *              Q1,Q2,FREQ)

      IUNIQ=0
      DO IXYZ=1,NNXYZ
         IF(NOCP(IXYZ).NE.1) IUNIQ=IUNIQ+1
         IF(NOCP(IXYZ).NE.0) GO TO 160
         DO IROT=1,NROT
            RESID(IROT,IUNIQ) = -WAXH2(IROT,IXYZ) - RESID(IROT,IUNIQ)
         ENDDO ! IROT
  160 CONTINUE
      ENDDO ! IXYZ
C
C        CALCULATE INITIAL PSEUDORESIDUAL -ZRES-, CHECK CONVERGENCE
C
        call zres_cphf_dyn1(nnxyz,nrot,NUNIQ,ERR,tol,WAXH2,bnorm,
     $          RESID,ZRES,NOCP,PRCND,WAXA)
C
      IF(masout) THEN
         IF(DIRSCF) THEN
            WRITE(IW,9010)
            WRITE(IW,9020) ITER,ERR,NXYZF,NINT,NSCHWZ
         ELSE
            WRITE(IW,9030)
            WRITE(IW,9020) ITER,ERR,NXYZF
         END IF
         CALL FLSHBF(IW)
      END IF
C
      IF(ERR.LT.TOL) GO TO 800

C        THE REMAINING CG ITERATIONS (2,3,...,MAXCP) START NOW...
C
      MODE=2
      BKDEN(1)  = 1.0D+00
C        ITERATION LIMIT, IF WE AREN'T MAKING PROGRESS, IS SMALLISH.
C        PRIOR TO 6/2012, ONLY 50 TOTAL ITERATIONS ALLOWED HERE,
C        relaxed to allow more iters so long as progress is made.
      MAXCP = 50
      MAXCP2= 299
      NFOCK = 0
      ITER=1
       ERR = ZERO
  300 CONTINUE
      ITER=ITER+1
      ERRSAVE = ERR
C
C           CALCULATE COEFFICIENT -BK- AND DIRECTION VECTOR -PDIR-
C
      CALL pdir_cphf_dyn1(NNXYZ,NROT,NUNIQ,iter,tol,bknum,bkden,
     $     NOCP,RESID,PDIR,ZRES,IW,masout)

      CALL AOCPTD_DYN(VEC,PDIR,ZRES,DTEMP,FTEMP,FTEMP2,EIG,WRK1,WRK2,
     *              XX,IX,GHONDO,XINTS,DSH,DDIJ,NOCP,MODE,
     *              NNXYZ,NUNIQ,NXYZF,NFOCK,NINT,NSCHWZ,
     *              NROT,NOCC,NVIR,L1,L2,L3,NSH2,MAXG,MXG2,
     *              Q1,Q2,FREQ)

         ERR = ZERO
         IUNIQ=0
         DO IXYZ=1,NNXYZ
            IF(NOCP(IXYZ).NE.1) IUNIQ=IUNIQ+1
            IF(NOCP(IXYZ).NE.0) GO TO 340
            call akden_cphf_dyn1(NROT,NUNIQ,PDIR,ZRES,akden,
     &                tol,IUNIQ,IW,masout)
            IF(ABS(BKNUM(IXYZ)).GT.1.0D-35 .AND. AKDEN.NE.0.0D+00) THEN
               AK = BKNUM(IXYZ)/AKDEN
            ELSE
               AK = 0.0D+00
            ENDIF
C ALPHA
            CALL DAXPY(NROT, AK,PDIR(1,IUNIQ),1,   YA(1, IXYZ),1)
            CALL DAXPY(NROT,-AK,ZRES(1,IUNIQ),1,RESID(1,IUNIQ),1)
            DO IROT=1,NROT
               ZRES(IROT,IUNIQ)  = PRCND(IROT)*RESID(IROT,IUNIQ)
            ENDDO
C CHECK CONVERGENCE
          call conv_cphf_dyn1(NROT,NUNIQ,ERR,test,bnorm(ixyz),IUNIQ,
     $              RESID,PRCND,WAXA)
            IF(TEST.LT.TOL) NOCP(IXYZ)=2
  340    CONTINUE
         ENDDO ! IXYZ
         CVGING = ERR.LT.TWO*ERRSAVE
         IF(masout) THEN
            IF(DIRSCF) THEN
               WRITE(IW,9020) ITER,ERR,NXYZF,NINT,NSCHWZ
            ELSE
               WRITE(IW,9020) ITER,ERR,NXYZF
            END IF
            CALL FLSHBF(IW)
         END IF
C
C              exit if converged
C
      IF(ERR.LT.TOL) GO TO 800
C
C        Another iteration?
C
      IF(MASWRK  .AND.  ITER.EQ.MAXCP) WRITE(IW,9090)
      IF(ITER.EQ.MAXCP2) CVGING=.FALSE.
      IF(ITER.LE.MAXCP) GO TO 300
      IF(CVGING)        GO TO 300
C
      IF(MASWRK) WRITE(IW,9070) MAXCP,MAXCP2
      IF(MASWRK) WRITE(IW,9060)
      CALL ABRT
      STOP
C
C        PRINT CONVERGENCE MESSAGE, RESTORE ORIGINAL -NOCP- ARRAY
C
  800 CONTINUE
      IF(masout) WRITE(IW,9080) ITER,NFOCK,NUNIQ
      DO IXYZ=1,NNXYZ
         IF(NOCP(IXYZ).EQ.2) NOCP(IXYZ)=0
      ENDDO
      NSTEP=ifreq

      !IF(MASWRK) THEN
      !   WRITE(IW,*) 'CPHF RESPONSE VECTORS -U- DYN', NLOC
      !   CALL PRSQ(YA,NNXYZ,NROT,NROT)
      !END IF

      CALL LAPOL_DYN(YA,YAL,HF,HFL,TRAN,DYNDD_LMO(1,1,Ifreq),
      !CALL LAPOL(YA,YAL,HF,HFL,TRAN,DYNDD_LMO(1,1,Ifreq),VCOC,
     *          IA,NNXYZ,L2,NOCC,NVIR,MCORE,NLOC,NSTEP,DYNDD(1,Ifreq))
      IF(IDQDYN .or. IDQ) THEN
         CALL LDQPOL(YA,YAL,HQ,HQL,TRAN,DYNDQ_LMO(1,1,Ifreq),VCOC,
     *               IA,NNXYZ,L2,NOCC,NVIR,MCORE,NLOC,NSTEP,
     *               DYNDQ(1,Ifreq))
      END IF
C
         !CALL TIMIT(1)
      ENDDO ! end Ifreq loop
      if(maswrk) write(iw,*) ' DONE FOR DYNAMIC POLARIZABILITY '
      RETURN
C
 9000 FORMAT(1X,'PRECONDITIONED CONJUGATE GRADIENT SOLVER',5X,
     *          'CONV. TOLERANCE=',1P,E8.2)
 9010 FORMAT(1X,12X,'MAXIMUM',10X,'RESPONSES',8X,'NONZERO',5X,'BLOCKS'/
     *       1X,'ITER',4X,'RESPONSE ERROR',8X,'IMPROVED',3X,
     *          'AO INTEGRALS',4X,'SKIPPED')
 9020 FORMAT(1X,I3,5X,1P,E13.5,0P,8X,I6,3X,I15,I11)
 9030 FORMAT(1X,'ITER',4X,'RESPONSE ERROR',8X,'IMPROVED')
 9031 FORMAT(/1X,'COMPUTING DYNAMIC POLARIZABILTY FOR IMAGINARY',
     *           ' FREQUENCY',F12.7,'I (A.U.)')
 9040 FORMAT(1X,'THE PRECONDITIONER IS NOT POSITIVE DEFINITE.'/
     *       1X,'THIS MAY BE DUE TO PECULIAR ORBITAL ENERGIES.')
c9050 FORMAT(//1X,'THE ORBITAL HESSIAN IS NOT POSITIVE DEFINITE.'//)
 9060 FORMAT(1X,'MOST OFTEN THIS IS',
     *          ' CAUSED BY USE OF AN INAPPROPRIATE WAVEFUNCTION.'/
     *       1X,'CHANGE SCFTYP, OR CHECK HOMO/LUMO FILLING ORDER.'//)
 9070 FORMAT(//1X,'*** TOO MANY ITERATIONS IN AOCPCG *** MAX CPHF=',I5/
     *       1X,'AS LONG AS EQUATIONS REMAIN CONVERGENT,',
     *          ' GRACE LIMIT=',I4/
     *       1X,'THIS VALUE CANNOT BE RAISED BY THE INPUT BECAUSE'/
     *       1X,'MORE CPHF ITERATIONS ARE UNLIKELY TO HELP.'//)
 9080 FORMAT(1X,'THE CPHF HAS CONVERGED AFTER',I3,' ITERATIONS.'/
     *       1X,'IT REQUIRED',I6,' FOCK-LIKE BUILDS TO FIND THE',I4,
     *          ' SYMMETRY UNIQUE RESPONSES.')
 9090 FORMAT(28X,'*** WARNING ***'/
     *       1X,'THE RESPONSE EQUATIONS HAVE NOW REACHED THE NORMAL',
     *          ' ITERATION LIMIT,'/
     *       1X,'BUT THE SOLVER WILL CONTINUE AS LONG AS THE',
     *          ' EQUATIONS ARE CONVERGING.')
      END
C
      subroutine constr_waxh2(WAX,WAXH2,VEC,EIG,NNXYZ,L1,L2,L3
     &            ,DTEMP,FTEMP2,DSH,DDIJ,NUNIQ,NROT,NOCC,NVIR,NOCP
     &            ,WRK1,WRK2,XX,IX,MAXG,MXG2,NSH2,NINT,NSCHWZ,IDAF,IODA)
      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(INOUT) :: WAXH2(NROT,NNXYZ)
     &         ,WAX(NROT,NNXYZ)

      INTEGER, INTENT(IN) :: NNXYZ,L1,L2,L3,NUNIQ,NROT,NOCC,NVIR
     &         ,MAXG,MXG2,NSH2,NINT,NSCHWZ, IDAF
      INTEGER, INTENT(IN) :: NOCP(NNXYZ), IX(*),IODA(950)
      DOUBLE PRECISION, INTENT(IN) :: VEC(L1,L1), EIG(L1)
     &         ,DTEMP(*), FTEMP2(*), DSH(*), DDIJ(*), WRK1(*), WRK2(*)
     &         ,XX(*)

      INTEGER :: MODE, NXYZF, IXYZ, lx
      DOUBLE PRECISION, PARAMETER :: HALF=0.5D+00
      DOUBLE PRECISION, PARAMETER :: ONE=1.0D+00
      DOUBLE PRECISION, PARAMETER :: TWO=2.0D+00
      DOUBLE PRECISION, PARAMETER :: FOUR=4.0D+00


C     ----- READ IN MO DIPOLE INTEGARLS -------
      WAX(:,:)=0.0D+00
      WAXH2(:,:)=0.0D+00
      call poldb_wax(wax,l2,nnxyz,nrot,nvir,nocc,l1)
      MODE  = 1
      NXYZF = 0
      DO IXYZ = 1,NNXYZ
         IF (NOCP(IXYZ).EQ.0) NXYZF = NXYZF+1
      ENDDO
C Construct WAXH2= 2* H(2) * P
      CALL cpdyn_AMBX(WAX,WAXH2,VEC,EIG,NNXYZ,L1,L2,L3,DTEMP,FTEMP2
     &      ,DSH,DDIJ,NUNIQ,NROT,NOCC,NVIR,NOCP,WRK1,WRK2,XX,IX
     &      ,NXYZF,MODE,MAXG,MXG2,NSH2,NINT,NSCHWZ)

C Scale with 2
      CALL DSCAL(NROT*NNXYZ,TWO,WAXH2,1)
      CALL DSCAL(NROT*NNXYZ,four,WAXH2,1)
      !write(iw,*) 'waxh2'
      !call prsq(waxh2,nnxyz,nrot,nrot)
      RETURN
      END
C
      SUBROUTINE poldb_wax(wax,l2,nnxyz,nrot,nvir,nocc,l1)

      use mx_limits, only: mxao
      IMPLICIT NONE
      
      double precision :: WAX(NROT,NNXYZ)
      DOUBLE PRECISION,DIMENSION(:,:), ALLOCATABLE :: HF

      COMMON /IJPAIR/ IA(MXAO)
      COMMON /IOFILE/ IR,IW,IP,IJKO,IJKT,IDAF,NAV,IODA(950)
C
      INTEGER IA,IDAF,IJKO,IJKT,IODA,IP,IR,IW,NAV,NNXYZ,NROT
      INTEGER IJ,IOCC,IROT,IVIR,L1,L2,NOCC,NVIR
C
      ALLOCATE(HF(L2,NNXYZ))
C
C     ----- READ IN MO DIPOLE INTEGARLS -------
C
      CALL DAREAD(IDAF,IODA,HF(1,1),L2,252,0)
      CALL DAREAD(IDAF,IODA,HF(1,2),L2,253,0)
      CALL DAREAD(IDAF,IODA,HF(1,3),L2,254,0)

      IROT = 0 
      DO IVIR=1,NVIR
       DO IOCC=1,NOCC
         IROT=IROT+1
         IJ=IA(IVIR+NOCC)+IOCC
         WAX(IROT,1) =  HF(IJ,1)
         WAX(IROT,2) =  HF(IJ,2)
         WAX(IROT,3) =  HF(IJ,3)
        ENDDO
      ENDDO
      DEALLOCATE(HF)
      RETURN
      END
C
      SUBROUTINE cpdyn_AMBX(YA,RHS,C,EIG,NNXYZ,L1,L2,L3,DTEMP,FTEMP2
     &               ,DSH,DDIJ,NUNIQ,NROT,NOCC,NVIR,NOCP,WRK1,WRK2
     &               ,XX,IX,NXYZF,MODE,MAXG,MXG2,NSH2,NINT,NSCHWZ)
      IMPLICIT NONE
      
      double precision :: YA(NROT,NNXYZ), C(L1,L1), EIG(L1)
      double precision :: RHS(NROT,NNXYZ)
      double precision :: DTEMP(L2,NUNIQ)
      double precision :: FTEMP2(L2,NUNIQ)
      double precision :: WRK1(L1,L1),WRK2(L3)
     &          ,XX(*),DSH(NSH2),DDIJ(49*MXG2)
      double precision :: GHONDO(MAXG),XINTS(NSH2)
      integer :: NOCP(NNXYZ),IX(*)
      integer :: MODE,MAXG,MXG2,NSH2,NINT,NSCHWZ
      integer :: NXYZF,NUNIQ,NROT,NOCC,NVIR,NNXYZ, L1,L2,L3
      DOUBLE PRECISION, PARAMETER :: HALF=0.5D+00
      DOUBLE PRECISION, PARAMETER :: ONE=1.0D+00
      DOUBLE PRECISION, PARAMETER :: TWO=2.0D+00

      COMMON /iofile/ ir,iw,ip,is,ipk,idaf,nav,ioda(950)
        INTEGER :: ir,iw,ip,is,ipk,idaf,nav,ioda

      INTEGER :: IXYZF,IUNIQ, I, J, IJ, IKOL, IXYZ, IVIR, IOCC
      double precision :: dummy=0.D0
      double precision :: ediff, temp1

      integer :: iprint=0

      CALL VCLR(DTEMP ,1,NUNIQ*L2)
      CALL VCLR(FTEMP2 ,1,NUNIQ*L2)

      IXYZF = 0
      IUNIQ = 0
      DO IXYZ = 1,NNXYZ
         IF(NOCP(IXYZ).NE.1) IUNIQ=IUNIQ+1
         IF(NOCP(IXYZ).NE.0) GO TO 270
         IXYZF = IXYZF+1
                       IKOL = IXYZ
         IF(MODE.EQ.2) IKOL = IUNIQ
         CALL MRARTR(YA(1,IKOL),NOCC,NOCC,NVIR,C(1,NOCC+1),
     *               L1,L1,WRK2,L1)
         CALL MRARBR(C,L1,L1,NOCC,WRK2,L1,L1,WRK1,L1)

         IJ = 0
         DO I=1,L1
            DO J=1,I
              ij= ij+1
              DTEMP(IJ,ixyzf) =-WRK1(I,J) + WRK1(J,I)
              !wrk2(ij) = -(WRK1(J,I) - WRK1(I,J))
            ENDDO ! J loop
         ENDDO ! I loop


       !call CPYSQT(wrk2,dtemp(1,ixyzf),l1,1)


       if(ixyz.eq.1.and.iprint.eq.1) then
         write(iw,*) 'dtemp a-b wrk2'
         call prsq(wrk2,l1,l1,l1)
         write(iw,*) 'dtemp a-b'
         call prtri(dtemp(1,ixyzf),l1)
       endif
  270 CONTINUE
      ENDDO ! IXYZ loop

C compute H2*YA in AO :  (A-B)*X
       CALL TWOEI_CPHFDYN_AMB(NINT,NSCHWZ,L1,L2,1,XINTS,NSH2,
     *                 GHONDO,MAXG,DDIJ,DTEMP,FTEMP2,DSH,NXYZF)
      IUNIQ = 0
      IXYZF=0
      DO IXYZ = 1,NNXYZ
         IF(NOCP(IXYZ).NE.1) IUNIQ=IUNIQ+1
         IF(NOCP(IXYZ).NE.0) GO TO 390
         IXYZF = IXYZF+1
C WRK2 = (A-B)X
         CALL DCOPY(L2,FTEMP2(1,IXYZF),1,WRK1,1)
         CALL EXPND(WRK1,WRK2,L1,1)

       if(ixyz.eq.1.and.iprint.eq.1) then
         write(iw,*) 'ftemp a-b wrk2'
         call prsq(wrk2,l1,l1,l1)
         endif

         CALL MRARBR(WRK2,L1,L1,L1,C(1,NOCC+1),L1,NVIR,WRK1,L1)
         CALL MRTRBR(C,L1,L1,NOCC,WRK1,L1,NVIR,WRK2,NOCC)

       if(ixyz.eq.1.and.iprint.eq.1) then
         write(iw,*) 'ftemp a-b wrk2 mo'
         call prsq(wrk2,l1,l1,l1)
         endif

            IJ = 0
                          IKOL = IXYZ
         IF(MODE.EQ.2) IKOL = IUNIQ
            DO IVIR = NOCC+1,NOCC+NVIR
               DO IOCC= 1,NOCC
                  IJ = IJ + 1
                  EDIFF= (EIG(IVIR) - EIG(IOCC))
C RHS = WRK2 + EDIFF*YA
                  RHS(IJ,IKOL)= 2.0D+00*WRK2(IJ) + EDIFF*YA(IJ,IKOL) 
               ENDDO ! IOCC loop
            ENDDO ! IVIR loop
  390 CONTINUE
      ENDDO ! IXYZ loop
      RETURN
      END
C
      SUBROUTINE cpdyn_APBX(YA,RHS,C,EIG,NNXYZ,L1,L2,L3,DTEMP,FTEMP
     &               ,DSH,DDIJ,NUNIQ,NROT,NOCC,NVIR,NOCP,WRK1,WRK2
     &               ,XX,IX,NXYZF,MODE,MAXG,MXG2,NSH2,NINT,NSCHWZ)
      IMPLICIT NONE
      
      double precision :: YA(NROT,NNXYZ), C(L1,L1), EIG(L1)
      double precision :: RHS(NROT,NNXYZ)
      double precision :: DTEMP(L2,NUNIQ)
      double precision :: FTEMP(L2,NUNIQ)
      double precision :: WRK1(L1,L1),WRK2(L3)
     &          ,XX(*),DSH(NSH2),DDIJ(49*MXG2)
      double precision :: GHONDO(MAXG),XINTS(NSH2)
      integer :: NOCP(NNXYZ), IX(*)
      integer :: MODE,MAXG,MXG2,NSH2,NINT,NSCHWZ
      integer :: NXYZF,NUNIQ,NROT,NOCC,NVIR,NNXYZ, L1,L2,L3
      DOUBLE PRECISION, PARAMETER :: HALF=0.5D+00
      DOUBLE PRECISION, PARAMETER :: ONE=1.0D+00
      DOUBLE PRECISION, PARAMETER :: TWO=2.0D+00

      INTEGER :: IXYZF,IUNIQ, I, J, IJ, IKOL, IXYZ, IVIR, IOCC
      double precision :: dummy=0.D0
      double precision :: ediff

      integer :: iprint

      CALL VCLR(DTEMP ,1,NUNIQ*L2)
      CALL VCLR(FTEMP ,1,NUNIQ*L2)

      IXYZF = 0
      IUNIQ = 0
      DO IXYZ = 1,NNXYZ
         IF(NOCP(IXYZ).NE.1) IUNIQ=IUNIQ+1
         IF(NOCP(IXYZ).NE.0) GO TO 270
         IXYZF = IXYZF+1
                       IKOL = IXYZ
         IF(MODE.EQ.2) IKOL = IUNIQ
         CALL MRARTR(YA(1,IKOL),NOCC,NOCC,NVIR,C(1,NOCC+1),
     *               L1,L1,WRK2,L1)
         CALL MRARBR(C,L1,L1,NOCC,WRK2,L1,L1,WRK1,L1)

         IJ = 0
         DO I=1,L1
            DO J=1,I
              ij= ij+1
                DTEMP(IJ,ixyzf) = WRK1(I,J) + WRK1(J,I)
            ENDDO ! J loop
         ENDDO ! I loop

  270 CONTINUE
      ENDDO ! IXYZ loop

C compute H1*YA in AO :  (A+B)*X
       CALL TWOEI_CPHFDYN_APB(NINT,NSCHWZ,L1,L2,1,XINTS,NSH2,
     *                 GHONDO,MAXG,DDIJ,DTEMP,FTEMP,DSH,NXYZF)

      call vclr(wrk1,1,L3)
      IUNIQ = 0
      IXYZF = 0
      DO IXYZ = 1,NNXYZ
         IF(NOCP(IXYZ).NE.1) IUNIQ=IUNIQ+1
         IF(NOCP(IXYZ).NE.0) GO TO 390
         IXYZF = IXYZF+1
C WRK2 = (A+B)X
         CALL DCOPY(L2,FTEMP(1,IXYZF),1,WRK1,1)
         CALL EXPND(WRK1,WRK2,L1,0)
         CALL MRARBR(WRK2,L1,L1,L1,C(1,NOCC+1),L1,NVIR,WRK1,L1)
         CALL MRTRBR(C,L1,L1,NOCC,WRK1,L1,NVIR,WRK2,NOCC)

            IJ = 0
                          IKOL = IXYZ
            IF(MODE.EQ.2) IKOL = IUNIQ
            DO IVIR = NOCC+1,NOCC+NVIR
               DO IOCC= 1,NOCC
                  IJ = IJ + 1
                  EDIFF= 4.0D+00*(EIG(IVIR) - EIG(IOCC))
C RHS = WRK2 + EDIFF*YA
                  RHS(IJ,IKOL)= 8.0D+00*WRK2(IJ) + EDIFF*YA(IJ,IKOL) 
               ENDDO ! IOCC loop
            ENDDO ! IVIR loop
  390 CONTINUE
      ENDDO ! IXYZ loop
      RETURN
      END

      SUBROUTINE AOCPTD_DYN(VEC,YA,RHS,DTEMP,FTEMP1,FTEMP2,EIG,WRK1,WRK2
     *              ,XX,IX,GHONDO,XINTS,DSH,DDIJ,NOCP,MODE
     *              ,NNXYZ,NUNIQ,NXYZF,NFOCK,NINT,NSCHWZ
     *              ,NROT,NOCC,NVIR,L1,L2,L3,NSH2,MAXG,MXG2
     *              ,Q1,Q2,FREQ)
      use constants, only: two,four
      IMPLICIT NONE
C
      DIMENSION VEC(L1,L1),YA(NROT,NNXYZ),
     *          RHS(NROT,NUNIQ),
     *          FTEMP1(L2,NUNIQ),
     *          WRK1(L1,L1),WRK2(L3),EIG(L1),NOCP(NNXYZ),
     *          GHONDO(MAXG),XINTS(NSH2),DSH(NSH2),DDIJ(49*MXG2),
     *          XX(*),IX(*),
     *          Q1(NROT,NNXYZ),Q2(NROT,NNXYZ),
     *          DTEMP(L2,NUNIQ),FTEMP2(L2,NUNIQ)
C
      DOUBLE PRECISION DDIJ,DSH,DTEMP,EIG,FTEMP1,FTEMP2,GHONDO,RHS,VEC
      DOUBLE PRECISION WRK1,WRK2,XINTS,XX,YA
      DOUBLE PRECISION Q1,Q2
      DOUBLE PRECISION FREQ,FREQ_I
C
      INTEGER IX,L1,L2,L3,MAXG,MXG2,NNXYZ,NOCP,NROT,NSH2,NUNIQ
      INTEGER IPRINT,IROT,IUNIQ,IW,IXYZ,MODE,NFOCK,NINT,NOCC,NSCHWZ
      INTEGER NVIR,NXYZF
C
      NXYZF = 0
      DO IXYZ = 1,NNXYZ
         IF (NOCP(IXYZ).EQ.0) NXYZF = NXYZF+1
      ENDDO
      NFOCK = NFOCK + 2*NXYZF
      CALL VCLR(Q1,1,NROT*NNXYZ)
      CALL VCLR(Q2,1,NROT*NNXYZ)

      iprint=0

      CALL cpdyn_APBX(YA,Q1,VEC,EIG,NNXYZ,L1,L2,L3,DTEMP,FTEMP1,DSH,DDIJ
     &      ,NUNIQ,NROT,NOCC,NVIR,NOCP,WRK1,WRK2,XX,IX,NXYZF,MODE
     &      ,MAXG,MXG2,NSH2,NINT,NSCHWZ)
       if(iprint.eq.1) then
       write(iw,*) 'q1 test'
       call prsq(q1,nnxyz,nrot,nrot)
       endif

      CALL cpdyn_AMBX(Q1,Q2,VEC,EIG,NNXYZ,L1,L2,L3,DTEMP,FTEMP2,DSH,DDIJ
     &      ,NUNIQ,NROT,NOCC,NVIR,NOCP,WRK1,WRK2,XX,IX,NXYZF,MODE
     &      ,MAXG,MXG2,NSH2,NINT,NSCHWZ)

       if(iprint.eq.1) then
       write(iw,*) 'q2'
       call prsq(q2,nnxyz,nrot,nrot)
       endif
      IUNIQ=0
      DO IXYZ=1,NNXYZ
         IF(NOCP(IXYZ).NE.1) IUNIQ=IUNIQ+1
         IF(NOCP(IXYZ).NE.0) GO TO 160
         DO IROT=1,NROT
            FREQ_I= 4.0D+00*FREQ*FREQ
            RHS(IROT,IUNIQ) = Q2(IROT,IXYZ) + FREQ_I* YA(IROT,IUNIQ)
         ENDDO !IROT
 160  CONTINUE
      ENDDO !IXYZ
       if(iprint.eq.1) then
       write(iw,*) 'rhs'
       call prsq(rhs,nnxyz,nrot,nrot)
       endif
      RETURN
      END

      subroutine iniv_cphf_dyn(nnxyz,nrot,NUNIQ,bknum,bnorm,
     $               bkden,RESID,PDIR,ZRES,NOCP)

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(INOUT) ::
     &   RESID(NROT,NUNIQ),
     &   ZRES(NROT,NUNIQ), 
     &   PDIR(NROT,NUNIQ), 
     &   BNORM(NNXYZ), BKNUM(NNXYZ), BKDEN(NNXYZ)

      INTEGER, INTENT(INOUT) ::  NOCP(NNXYZ)

      INTEGER NNXYZ,NROT,NUNIQ,IXYZ

      call vclr(bknum ,1,nnxyz)
      call vclr(bnorm ,1,nnxyz)
      call vclr(bkden ,1,nnxyz)
      call vclr(pdir  ,1,NUNIQ*nrot)
      call vclr(RESID ,1,NUNIQ*nrot)
      call vclr(ZRES,1,NUNIQ*nrot)
      do ixyz=1,nnxyz
       nocp(ixyz)=0
      enddo
      return
      end
C
      subroutine zres_cphf_dyn1(nnxyz,nrot,NUNIQ,ERR,tol,WAXH2,bnorm,
     $               RESID,ZRES,NOCP,PRCND,WAXA)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NNXYZ, nrot,NUNIQ

      INTEGER, INTENT(INOUT) :: NOCP(NNXYZ)
      DOUBLE PRECISION, INTENT(INOUT) ::
     &   ZRES(NROT,NUNIQ), BNORM(NNXYZ)

      DOUBLE PRECISION, INTENT(IN) :: TOL
      DOUBLE PRECISION, INTENT(IN) :: 
     &   RESID(NROT,NUNIQ), PRCND(NROT), WAXH2(NROT,NNXYZ)

      DOUBLE PRECISION, INTENT(OUT) :: ERR

      INTEGER :: IUNIQ, IXYZ, IROT
      DOUBLE PRECISION :: DDOT, TEST
      DOUBLE PRECISION :: WAXA(NROT)

      DOUBLE PRECISION, PARAMETER :: ZERO=0.0D+00

      ERR = ZERO
      IUNIQ=0
      DO IXYZ=1,NNXYZ
       IF(NOCP(IXYZ).NE.1) IUNIQ=IUNIQ+1
       IF(NOCP(IXYZ).NE.0) GO TO 210
       DO IROT=1,NROT
          ZRES(IROT,IUNIQ) = PRCND(IROT)*RESID(IROT,IUNIQ)
       ENDDO

       DO IROT=1,NROT
          WAXA(IROT) =  PRCND(IROT) *WAXH2(IROT,IUNIQ)
       ENDDO

       !BNORM(IXYZ) = DDOT(NROT, WAXH2(1,IXYZ),1, WAXH2(1, IXYZ),1)
       BNORM(IXYZ) = DDOT(NROT, WAXA,1, WAXA,1)
       BNORM(IXYZ) = SQRT(BNORM(IXYZ))
C check convergence
       call conv_cphf_dyn1(NROT,NUNIQ,ERR,test,bnorm(ixyz),IUNIQ,
     $           RESID,PRCND,WAXA)

         IF(TEST.LT.TOL) NOCP(IXYZ)=2
  210 CONTINUE
      ENDDO
      return
      end
      subroutine conv_cphf_dyn1(NROT,NUNIQ,ERR,test,bnorm,IUNIQ,
     $               RESID,PRCND,TRESID)

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(INOUT) :: ERR, BNORM

      INTEGER, INTENT(IN) :: NROT,NUNIQ,IUNIQ
      DOUBLE PRECISION, INTENT(IN) :: 
     &   RESID(NROT,NUNIQ), PRCND(NROT)

      DOUBLE PRECISION, INTENT(INOUT) :: TEST
      DOUBLE PRECISION :: RNORM, DDOT, RNORMWAX
      DOUBLE PRECISION :: TRESID(NROT)
      INTEGER :: IROT

C note
C LET A= (prcnd*resid)^T * (prcnd*resid)
C and B= (prcnd*wax)^T * (prcnd*wax)
C IF RNORM  = A/B, the its PCG convergence can be faster 
C than |resid|/|wax|. However, the numerical results between these two
C convergence conditions will be slightly different.

C COMPUTE A= (prcnd*resid)^T * (prcnd*resid)
       DO IROT=1,NROT
           TRESID(IROT) =  PRCND(IROT) * RESID(IROT,IUNIQ)
       ENDDO
C B= BNORM

       !RNORM       = DDOT(NROT,RESID(1,IUNIQ),1,RESID(1,IUNIQ),1)
       RNORM       = DDOT(NROT,TRESID,1,TRESID,1)
       RNORM       = SQRT(RNORM)
       TEST = RNORM/BNORM
C       write(iw,901) 'bnorm', BNORM, RNORM, TEST
       ERR = MAX(ERR,TEST)
      RETURN
 900  format(7F15.7)
 901  format(A,5F15.7)
      END
      subroutine pdir_cphf_dyn1(NNXYZ,NROT,NUNIQ,iter,tol,bknum,bkden
     $               ,NOCP,RESID,PDIR,ZRES,IW,MASWRK)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NNXYZ, NROT,NUNIQ,IW,ITER
      DOUBLE PRECISION, INTENT(IN) :: TOL
      DOUBLE PRECISION, INTENT(IN) :: 
     &   RESID(NROT,NUNIQ)

      DOUBLE PRECISION, INTENT(INOUT) ::
     &   ZRES(NROT,NUNIQ), PDIR(NROT,NUNIQ),
     &   BKNUM(NNXYZ), BKDEN(NNXYZ)
      INTEGER, INTENT(INOUT) :: NOCP(NNXYZ)

      LOGICAL, INTENT(IN) :: MASWRK

      INTEGER :: IUNIQ, IXYZ, IROT
      DOUBLE PRECISION :: DDOT, TEST, BK
      DOUBLE PRECISION, PARAMETER :: ZERO=0.0D+00

C      write(iw,*) 'iter', iter
      IUNIQ=0
      DO IXYZ=1,NNXYZ
         IF(NOCP(IXYZ).NE.1) IUNIQ=IUNIQ+1
         IF(NOCP(IXYZ).NE.0) GO TO 320

         BK = ZERO
         BKNUM(IXYZ) =DDOT(NROT,ZRES(1,IUNIQ),1,RESID(1,IUNIQ),1)
         IF(BKNUM(IXYZ).LE.-TOL) THEN
            IF(MASWRK) WRITE(IW,9040)
            CALL ABRT
         END IF
         IF(ITER.EQ.2) THEN
            CALL DCOPY(NROT,ZRES(1,IUNIQ),1,PDIR(1,IUNIQ),1)
         ELSE
            BK = BKNUM(IXYZ)/BKDEN(IXYZ)
            DO IROT=1,NROT
               PDIR(IROT,IUNIQ) =    ZRES(IROT,IUNIQ)
     *                          +  BK*PDIR(IROT,IUNIQ)
            ENDDO
         END IF
         BKDEN(IXYZ)  = BKNUM(IXYZ)
  320 CONTINUE
      ENDDO
      RETURN
 9040 FORMAT(1X,'THE PRECONDITIONER IS NOT POSITIVE DEFINITE.'/
     *       1X,'THIS MAY BE DUE TO PECULIAR ORBITAL ENERGIES.')
      END
C
      subroutine akden_cphf_dyn1(NROT,NUNIQ,PDIR,ZRES,akden,
     &                 tol,IUNIQ,IW,MASWRK)

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(OUT) :: AKDEN

      INTEGER, INTENT(IN) :: NROT,NUNIQ,IUNIQ,IW
      DOUBLE PRECISION, INTENT(IN) :: TOL
      DOUBLE PRECISION, INTENT(IN) :: 
     &   ZRES(NROT,NUNIQ), PDIR(NROT,NUNIQ)

      LOGICAL,INTENT(IN) :: MASWRK

      DOUBLE PRECISION :: DDOT

        AKDEN = DDOT(NROT,PDIR(1,IUNIQ),1,ZRES(1,IUNIQ),1)
C
C          THIS TERMINATION IS SUPPOSED TO BE DUE TO THE
C          MATRIX -A-'S NOT BEING POSITIVE DEFINITE.  A TEST
C          CALCULATION SHOWED THAT A CASE WITH UHF INSTABILITY
C          DID NOT GET INTO THIS ERROR TERMINATION SECTION,
C          SO AN INSTABILITY IN THE ORBITAL HESSIAN MAY BE A
C          NECESSARY BUT NOT SUFFICIENT CONDITION TO GET HERE.
C
        IF(AKDEN.LE.-TOL) THEN
          IF(MASWRK) WRITE(IW,9060)
          CALL ABRT
        END IF

      RETURN
 9060 FORMAT(1X,'MOST OFTEN THIS IS',
     *          ' CAUSED BY USE OF AN INAPPROPRIATE WAVEFUNCTION.'/
     *       1X,'CHANGE SCFTYP, OR CHECK HOMO/LUMO FILLING ORDER.'//)
      END
C
C*MODULE LOCPOL    *DECK LAPOL_DYN
C> @brief  calculate molecular dipole polarizability and LMO contribution
C>
C> @author unknown (probably Simon Webb)
C>
C> @date July 13, 2014 -Peng Xu-
C>      add: calculate molecular dipole-dipole polarizability
C> @date January 2017 - C.Bertoni
C>      Changes to store coefficients for EFMO gradient
C>
C> @param NSTEP : flag for whether this calculates the distributed
C>        LMO polarizability for the static or dynamic polarizability
C>        if this is 0, it's static polarizability. if it's not 0,
C>        it's one of the frequencies for a corresponding the dynamic
C>        polarizability
C>
C> @param TPOL : molecular dipole-dipole polarizability
C>
      SUBROUTINE LAPOL_DYN(U,UL,HF,HFL,TRAN,DLPOL,IA,NXYZ,
     *                  NUM2,NOCC,NVIR,MCORE,NLOC,NSTEP,tpol)
      use mx_limits, only: mxatm,mxao
      use constants, only: zero,two,three,four
      USE EFP_LOGICAL
      USE MAKEFP_CPHF, only:makefp_print
C
      IMPLICIT NONE
C
      LOGICAL GOPARR,DSKWRK,MASWRK,MOIDON,EDCOMP,DIPDCM,QADDCM,
     *        DEPRNT,ZDO,POLDCM,POLANG,POLAPP,KMIDPT,POLDYN
C	
      INTEGER, PARAMETER :: NMO=500
C
      DIMENSION U(NOCC,NVIR,NXYZ),
     *          UL(NVIR,3),HF(NUM2,3),HFL(NVIR,3),TRAN(NLOC,NLOC),
     *          IA(NUM),CAPOL(9),APOLANG(9),DLPOL(9,NLOC)

      double precision :: tpol(9)
C
      COMMON /EDCMP / ZIJ(NMO),ZMO(5,NMO),OCCUP(NMO),DPFREQ(50),
     *                MOIDNO(5,NMO),IJMO(2,NMO),MOIJ(NMO),NMOIJ(NMO),
     *                NMOAT(NMO),NDPFREQ,IPROT(5),NPROT,
     *                MOIDON,EDCOMP,DIPDCM,DEPRNT,QADDCM,ZDO,POLDCM,
     *                POLANG,POLAPP,KMIDPT,POLDYN
      COMMON /EFMO  / EFMOETOT,EFMOESERG,EFMOPOLERG,EFMODISERG,
     *                EFMOREPNRG,EFMOCHTNRG,EFMOEPEN,EFMOPCMG,
     *                IEFMORUN,IMODEFP,IEFMORT,IEFMOCFRG,IEFMONFRG,
     *                IEFMODIM(2),IMODEFE,NATEFMO,IMODEFD,IMODEFER,
     *                IMODEFCT,IDIMTYP,iefmo_agrad
      COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,
     *                ZAN(MXATM),C(3,MXATM),IAN(MXATM)
      COMMON /IOFILE/ IR,IW,IP,IJK,IJKT,IDAF,NAV,IODA(950)
      COMMON /OPTLOC/ CVGLOC,MAXLOC,IPRTLO,ISYMLO,IFCORE,NOUTA,NOUTB,
     *                MOOUTA(MXAO),MOOUTB(MXAO),IBOYAL
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      logical masout
      COMMON /OUTPUT/ NPRINT,ITOL,ICUT,NORMF,NORMP,NOPK
C
      DOUBLE PRECISION, PARAMETER :: UNIT=0.52917724924D+00
C
      DOUBLE PRECISION APOLANG,CAPOL,DLPOL,DPFREQ,HF,HFL,OCCUP,TRAN,U
      DOUBLE PRECISION UL,ZIJ,ZMO
      DOUBLE PRECISION EFMOCHTNRG,EFMODISERG,EFMOEPEN,EFMOESERG
      DOUBLE PRECISION EFMOETOT,EFMOPCMG,EFMOPOLERG,EFMOREPNRG
      DOUBLE PRECISION C,CVGLOC,ZAN
      DOUBLE PRECISION ALPHAM,DUM
C
      INTEGER IA,IJMO,MOIDNO,MOIJ,NLOC,NMOAT,NMOIJ,NOCC,NUM,NUM2
      INTEGER NVIR,NXYZ
      INTEGER IDIMTYP,IEFMO_AGRAD,IEFMOCFRG,IEFMODIM,IEFMONFRG
      INTEGER IEFMORT,IEFMORUN,IMODEFCT,IMODEFD,IMODEFE,IMODEFER
      INTEGER IMODEFP,IPROT,NAT,NATEFMO,NDPFREQ,NPROT
      INTEGER IAN,ICH,IDAF,IFCORE,IJK,IJKT,IODA,IP,IPRTLO,IR,ISYMLO
      INTEGER IW,MAXLOC,MOOUTA,MUL,NA,NAV,NB,NE,NOUTA,NOUTB,NQMT
      INTEGER IBOYAL,IBTYP,ICUT,IPTIM,ITOL,MASTER,ME,MOOUTB,NOPK,NORMF
      INTEGER NORMP,NPRINT,NPROC
      INTEGER I,IJ,IL,ILO,IOCC,IVIR,J,K,KL,L,LOC,MCORE,NSTEP
C
      masout=maswrk.and.makefp_print
C
      CALL DAREAD(IDAF,IODA,HF(1,1),NUM2,252,0)
      CALL DAREAD(IDAF,IODA,HF(1,2),NUM2,253,0)
      CALL DAREAD(IDAF,IODA,HF(1,3),NUM2,254,0)
C
C     ---- CALCULATE ALPHA POLARIZABILITY TENSOR FOR THE CORE ----
C
      IF (IFCORE.EQ.1) THEN
      DO 400 K=1,3
         DO 410 L=1,3
         KL=(K-1)*3 + L
         CAPOL(KL)=ZERO
            DO 420 IVIR=1,NVIR
               DO 430 IOCC=1,MCORE
                  IJ=IA(IVIR+NOCC)+IOCC
                     IF(NSTEP.EQ.0) THEN
                     CAPOL(KL)=CAPOL(KL)-(FOUR*U(IOCC,IVIR,K)*HF(IJ,L))
                     ELSE
                     CAPOL(KL)=CAPOL(KL)-(TWO*U(IOCC,IVIR,K)*HF(IJ,L))
                     END IF
  430          CONTINUE
  420       CONTINUE
  410    CONTINUE
  400 CONTINUE
C
C     ----- CONVERT UNITS FROM BOHRS TO ANGTROMS ----
C
      DO 435 I=1,9
         IF(POLANG) THEN
            APOLANG(I)=CAPOL(I)*UNIT**THREE
         ELSE
            APOLANG(I)=CAPOL(I)
         END IF
  435 CONTINUE
C
C     ----- PRINT CORE ALPHA POLARIZABILITY MATRIX -----
C
      IF (masout) THEN
         IF(POLANG) THEN
            WRITE(IW,9020)
         ELSE
            WRITE(IW,9025)
         END IF
         WRITE(IW,9040)
         WRITE(IW,9060) (APOLANG(I),I=1,3)
         WRITE(IW,9080) (APOLANG(I),I=4,6)
         WRITE(IW,9100) (APOLANG(I),I=7,9)
      END IF
C
C     ---- CALCULATE THE MEAN CORE  POLARIZABILITY ----
C
      ALPHAM=(APOLANG(1)+APOLANG(5)+APOLANG(9))/THREE
      IF(masout) WRITE(IW,9110) ALPHAM
C
      ELSE
         CALL VCLR(CAPOL,1,9)
      END IF
C
C P. Xu:  -- calculate molecular dip-dip polarizability:TPOL --
C
      IF(IDD) THEN
      DO 700 K=1,3
         DO 710 L=1,3
         KL=(K-1)*3 + L
         TPOL(KL)=ZERO
            DO 720 IVIR=1,NVIR
               DO 730 IOCC=MCORE+1,NOCC
                  IJ=IA(IVIR+NOCC)+IOCC
                     IF(NSTEP.EQ.0) THEN
                     TPOL(KL)=TPOL(KL)-(FOUR*U(IOCC,IVIR,K)*HF(IJ,L))
                     ELSE
                     TPOL(KL)=TPOL(KL)-(TWO*U(IOCC,IVIR,K)*HF(IJ,L))
                     END IF
  730          CONTINUE
  720       CONTINUE
  710    CONTINUE
  700 CONTINUE
      END IF

C     -- CALCULATE ALPHA POLARIZABILITY TENSOR FOR LOCALIZED ORBITALS --
C
      CALL DAREAD(IDAF,IODA,TRAN,NLOC*NLOC,73,0)
C
c$$$ nstep=0 means that it's the static polarizability
         if( iefmorun .gt. 0 .and. masout .and.
     *     iefmo_agrad .gt. 0 ) then
            call efmo_store_ug( IEFMOCFRG, u, nocc, nvir,nstep+1 )
          endif

      DO 440 LOC=1,NLOC
         CALL VCLR(HFL,1,3*NVIR)
         CALL VCLR(UL,1,3*NVIR)
         DO 450 IVIR=1,NVIR
            IL=0
            DO 460 IOCC=MCORE+1,NOCC
               IL=IL+1
               IJ=IA(IVIR+NOCC)+IOCC
               HFL(IVIR,1)=HFL(IVIR,1)+(HF(IJ,1)*TRAN(IL,LOC))
               HFL(IVIR,2)=HFL(IVIR,2)+(HF(IJ,2)*TRAN(IL,LOC))
               HFL(IVIR,3)=HFL(IVIR,3)+(HF(IJ,3)*TRAN(IL,LOC))
               UL(IVIR,1)=UL(IVIR,1)+(U(IOCC,IVIR,1)*TRAN(IL,LOC))
               UL(IVIR,2)=UL(IVIR,2)+(U(IOCC,IVIR,2)*TRAN(IL,LOC))
               UL(IVIR,3)=UL(IVIR,3)+(U(IOCC,IVIR,3)*TRAN(IL,LOC))
  460       CONTINUE
  450    CONTINUE
C
         DO 470 K=1,3
            DO 480 L=1,3
               KL=(K-1)*3 + L
               DLPOL(KL,LOC)=ZERO
               DO 490 IVIR=1,NVIR
                  IF(NSTEP.EQ.0) THEN
                  DUM=(FOUR*UL(IVIR,K)*HFL(IVIR,L))
                  ELSE
                  DUM=(TWO*UL(IVIR,K)*HFL(IVIR,L))
                  END IF
                  DLPOL(KL,LOC)=DLPOL(KL,LOC)-DUM
                  CAPOL(KL)=CAPOL(KL)-DUM
  490          CONTINUE
  480       CONTINUE
  470    CONTINUE
  440 CONTINUE
C
      DO 510 I=1,NLOC
         DO 520 J=1,9
            IF(POLANG) THEN
               APOLANG(J)=DLPOL(J,I)*UNIT**THREE
            ELSE
               APOLANG(J)=DLPOL(J,I)
            END IF
  520    CONTINUE
C
         ILO=I+MCORE
         IF (masout) THEN
            IF(NMOAT(ILO).EQ.1) WRITE(IW,9000) MOIDNO(1,ILO)
            IF(NMOAT(ILO).EQ.2) WRITE(IW,9005)
     *                          MOIDNO(1,ILO),MOIDNO(2,ILO)
            WRITE(IW,9050)
            WRITE(IW,9060) (APOLANG(K),K=1,3)
            WRITE(IW,9080) (APOLANG(K),K=4,6)
            WRITE(IW,9100) (APOLANG(K),K=7,9)
            ALPHAM=(APOLANG(1)+APOLANG(5)+APOLANG(9))/THREE
            WRITE(IW,9110) ALPHAM
         END IF
  510 CONTINUE
C
C     ----- PRINT TOTAL ALPHA POLARIZABILITY MATRIX -----
C
      DO 550 I=1,9
         IF(POLANG) THEN
            APOLANG(I)=CAPOL(I)*UNIT**THREE
         ELSE
            APOLANG(I)=CAPOL(I)
         END IF
  550 CONTINUE
C
      IF (masout) THEN
         WRITE(IW,9055)
         WRITE(IW,9060) (APOLANG(I),I=1,3)
         WRITE(IW,9080) (APOLANG(I),I=4,6)
         WRITE(IW,9100) (APOLANG(I),I=7,9)
         ALPHAM=(APOLANG(1)+APOLANG(5)+APOLANG(9))/THREE
         WRITE(IW,9110) ALPHAM
      END IF
      RETURN
C
 9020 FORMAT(/10X,35(1H-)/
     *        10X,' LOCALIZED ALPHA POLARIZABILITIES '/
     *        10X,'       IN ANGSTROMS**3            '/
     *        10X,35(1H-)/)
 9025 FORMAT(//10X,35(1H-)/
     *        10X,' LOCALIZED ALPHA POLARIZABILITIES '/
     *        10X,'       IN ATOMIC UNITS            '/
     *        10X,35(1H-)/)
 9000 FORMAT(//1X,' LMO ALPHA POLARIZABILITY TENSOR ',
     *        'FOR CORE OR LONE PAIR ON ATOM',I3)
 9005 FORMAT(//1X,' LMO ALPHA POLARIZABILITY TENSOR ',
     *        'FOR BOND BETWEEN ATOM',I3,' AND ATOM',I3)
 9010 FORMAT(//10X,47(1H-)/
     *        10X,' LOCALIZED ALPHA POLARIZABILITY TENSORS FOR '/
     *        10X,' LONE PAIRS AND BONDS PROJECTED ONTO THE Z-AXIS '/
     *        10X,47(1H-))
 9040 FORMAT(10X,' CORE ALPHA POLARIZABILITY TENSOR '/
     *       22X,'UX',13X,'UY',13X,'UZ')
 9050 FORMAT(22X,'UX',13X,'UY',13X,'UZ')
 9055 FORMAT(//15X,35(1H-)/
     *        15X,' TOTAL ALPHA POLARIZABILITY TENSOR '/
     *        15X,35(1H-)/
     *        22X,'UX',13X,'UY',13X,'UZ')
 9060 FORMAT(10X,     ' UX ',3F15.9)
 9080 FORMAT(10X,     ' UY ',3F15.9)
 9100 FORMAT(10X,     ' UZ ',3F15.9)
 9110 FORMAT(/10X,     ' MEAN ALPHA POLARIZABILITY = ',3F15.9)
      END

      SUBROUTINE TWOEI_CPHFDYN_APB(NINT,NSCHWZ,L1,L2,NINTMX,XINTS,NSH2,
     *                 GHONDO,MAXG,DDIJ,DA,FA,DSH,NXYZ)
      use mx_limits, only: mxsh,mxgtot,mxatm,MXAO
      use constants, only: zero,half
C$    USE params, ONLY: intomp
C$    use TWOEI_CPHFDYN_omp
      use MOD_OFFLOAD, only: offload_makefp
      USE MAKEFP_CPHF, only: makefp_print_default
C
      IMPLICIT NONE
C
      LOGICAL OUT,SCHWRZ,SCHSKP,GOPARR,DSKWRK,MASWRK,DLB,SLB,C1GRP
      LOGICAL SKIPA,SKIPB,SKIPC,NPSYM
      LOGICAL PK,PANDK,NOTPK,BLOCK,GPSAVE,SCREEN
C
      DIMENSION XINTS(NSH2),
     *          GHONDO(MAXG),
     *          DSH(NSH2),DDIJ(*)
      DIMENSION FA(L2,NXYZ), DA(L2,NXYZ)
      DIMENSION MI(48),MJ(48),MK(48),M0(48)
C
      COMMON /IJPAIR/ IA(MXAO)
      COMMON /INTOPT/ ISCHWZ,IECP,NECP,IEFLD
      COMMON /INT2IC/ NINTIC,ININTIC,NXXIC,LBUFPIC,LIXIC,LABSIX,NINTIX
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON /NSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),
     *                CF(MXGTOT),CG(MXGTOT),CH(MXGTOT),CI(MXGTOT),
     *                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),
     *                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
      COMMON /OUTPUT/ NPRINT,ITOL,ICUT,NORMF,NORMP,NOPK
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /PKFIL / PK,PANDK,BLOCK
      COMMON /RESTAR/ TIMLIM,IREST,NREC,INTLOC,IST,JST,KST,LST
      COMMON /SCINP / VLAMB,SCREEN
      COMMON /SHLEXC/ NORGSH(3),NORGSP(3),IEXCH,NANGM,NGTH(4)
      COMMON /SHLG70/ ISH,JSH,KSH,LSH,IJKLXX(4)
      COMMON /SHLNOS/ QQ4,LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,
     *                MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL,
     *                NIJ,IJ,KL,IJKL
      COMMON /SHLT  / TOL,CUTOFF,ICOUNT,OUT
      COMMON /SYMTRY/ MAPSHL(MXSH,48),MAPCTR(MXATM,48),
     *                T(432),INVT(48),NT
      COMMON /ERIOUT/ INW,JNW,KNW,LNW,LSTRI,LSTRJ,LSTRK,LSTRL
      COMMON /DFTPAR/ DFTTYP(20),EXENA,EXENB,EXENC,
     *                IDFT34,NAUXFUN,NAUXSHL
      !COMMON /MACHSW/ KDIAG,ICORFL,IXDR,modio,mem10,lpnt10,mem10m
      logical not_gpu, masout
C
      DOUBLE PRECISION DA,DDIJ,DSH,FA,GHONDO,XINTS
      DOUBLE PRECISION CD,CF,CG,CH,CI,CP,CS,EX
      DOUBLE PRECISION TIMLIM,VLAMB
      DOUBLE PRECISION QQ4
      DOUBLE PRECISION CUTOFF,DFTTYP,EXENA,EXENB,EXENC,T,TOL
      DOUBLE PRECISION CSCALT,CUTSV,DENMAX,HFSCAL,TIM
      DOUBLE PRECISION TEST
C
      INTEGER IA,IECP,IEFLD,ININTIC,ISCHWZ,L2,LABSIX,LBUFPIC,LIXIC,M0
      INTEGER MAXG,MI,MJ,MK,NECP,NINTIC,NSH2,NXXIC,NXYZ
      INTEGER IDAF,IODA,IP,IPK,IR,IS,IW,KATOM,KLOC,KMAX,KMIN,KNG,KSTART
      INTEGER KTYPE,NAV,NINTIX,NSHELL
      INTEGER IBTYP,ICUT,IEXCH,INTLOC,IPTIM,IREST,IST,ITOL,JST,KST,LST
      INTEGER MASTER,ME,NANGM,NGTH,NOPK,NORGSH,NORGSP,NORMF,NORMP
      INTEGER NPRINT,NPROC,NREC
      INTEGER IJ,IJKLXX,ISH,JSH,KL,KSH,LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK
      INTEGER LOCL,LSH,MAXI,MAXJ,MAXK,MAXL,MINI,MINJ,MINK,MINL,NIJ
      INTEGER ICOUNT,IDFT34,IJKL,INVT,INW,JNW,KNW,LNW,LSTRI,LSTRJ,LSTRK
      INTEGER LSTRL,MAPCTR,MAPSHL,NAUXFUN,NAUXSHL,NT
      INTEGER I,II,IJIJ,IPCOUNT,JJ,JORK,KK,KLKL,L1,LL,LMAX,MINE,NEXT
      INTEGER NINT,NINTMX,NSCHWZ,NTTEMP
      INTEGER J
C
      MASOUT=MASWRK.and.makefp_print_default
C
C          ----- TWO-ELECTRON INTEGRALS -----
C     THIS VERSION CAN HANDLE S,P,D,F,G,H,I AND L SHELLS
C
      TIM = ZERO
      CALL TSECND(TIM)
C
      QQ4=1.0D+00
      ICOUNT=  0
      GPSAVE = GOPARR
      CUTSV  = CUTOFF
      CUTOFF = MIN(CUTOFF,1.0D-10)
      NTTEMP = NT
      NT = 1
      HFSCAL=DFTTYP(3)
      CSCALT=1.0D+00
C      IF(DIRTRF) GOPARR=.FALSE.
C
C     ----- INITIALIZATION FOR PARALLEL WORK -----
C     BOTH STATIC AND DYNAMIC LOAD BALANCING ARE IMPLEMENTED BELOW
C
      SLB = GOPARR  .AND.  IBTYP.EQ.0
      DLB = GOPARR  .AND.  IBTYP.EQ.1
      NEXT = -1
      MINE = -1
      IPCOUNT = ME - 1
C
      C1GRP = NT.EQ.1
C
      CALL BASCHK(LMAX)
                    NANGM =  4
      IF(LMAX.EQ.2) NANGM =  6
      IF(LMAX.EQ.3) NANGM = 10
      IF(LMAX.EQ.4) NANGM = 15
      IF(LMAX.EQ.5) NANGM = 21
      IF(LMAX.EQ.6) NANGM = 28
      NGTH(4) = 1
      NGTH(3) = NGTH(4) * NANGM
      NGTH(2) = NGTH(3) * NANGM
      NGTH(1) = NGTH(2) * NANGM
      IF(NOPK.EQ.0) THEN
         NORGSH(1) = 0
         NORGSH(2) = NORGSH(1) + NANGM**4
         NORGSH(3) = NORGSH(2) + NANGM**4
         NORGSP(1) = 0
         NORGSP(2) = 256
         NORGSP(3) = 512
      ELSE
         DO I=1,3
            NORGSH(I) = 0
            NORGSP(I) = 0
         ENDDO
      END IF
C
      NOTPK = .NOT.PK
      NINT  = 0
      NSCHWZ= 0
      SCHSKP=.FALSE.
      DENMAX = ZERO
         SCHWRZ = ISCHWZ.EQ.1
         IF(SCHWRZ) THEN
            CALL DAREAD(IDAF,IODA,XINTS,NSH2,54,0)
            CALL AOSHLDxyz(DA,DSH,IA,L1,L2,NSH2,NXYZ)
         END IF
C========================      OPENMP CODE       ======================
!$    if(intomp.gt.0) then
!$    CALL  TWOEI_CPHFDYN_APB_omp(SCHWRZ,NINT,NSCHWZ,L1,L2,XINTS,
!$   &        NSH2,MAXG,IA,DA,FA,DSH,NXYZ,cutoff,.false.,HFSCAL,CSCALT)
!$    goto 930
!$    endif
C========================       END OPENMP       ======================
C
C        NOW WE ARE READY TO LOOP OVER ALL NSHELL**4 SHELL QUARTETS
C
C     ----- I SHELL -----
C
      DO 920 II = 1,NSHELL
      DO 900 JJ = 1,II
C
C     ----- GO PARALLEL! -----
C
      IF (DLB) THEN
         MINE = MINE + 1
         IF (MINE.GT.NEXT) CALL DDI_DLBNEXT(NEXT)
         IF (NEXT.NE.MINE) GO TO 900
      END IF
      DO 880 KK = 1,II
      jork=kk
      if(ii.eq.kk) jork=jj
C
C     ----- GO PARALLEL! -----
C
      IF(SLB) THEN
         IPCOUNT = IPCOUNT + 1
         IF (MOD(IPCOUNT,NPROC).NE.0) GO TO 880
      END IF
      DO 860 LL = 1,jork
      ISH=II
      JSH=JJ
      KSH=KK
      LSH=LL
C
C     APPLY THE SCHWARZ INEQUALITY TO SCREEN OUT SMALL INTEGRALS,
C     SEE, FOR EXAMPLE, J.L.WHITTEN, J.CHEM.PHYS. 58,4496-4501(1973)
C
      IF(SCHWRZ) THEN
         IJIJ = (ISH*ISH-ISH)/2 + JSH
         KLKL = (KSH*KSH-KSH)/2 + LSH
         TEST = QQ4*XINTS(IJIJ)*XINTS(KLKL)
            !DENMAX = SCHWDN(DSH,ISH,JSH,KSH,LSH,IA)
            !TEST = TEST*DENMAX
         SCHSKP = TEST.LT.CUTOFF
         IF(SCHSKP) NSCHWZ = NSCHWZ + 1
      END IF
      IF(SCHSKP) GO TO 820
C
C        ----- ELECTRON REPULSION INTEGRAL CALCULATION -----
C     THIS MAY USE ROTATED AXIS, ERIC, OR RYS QUADRATURE METHODS
C
C update shell index to ISH,JSH,KSH,LSH in /SHLG70/
      CALL SHELLQUART(ISH,JSH,KSH,LSH,GHONDO)
C need to get shell index from /ERIOUT/
         CALL DIRFCK_CPHFDYN_APB(IA,DA,FA,NXYZ,INW,JNW,KNW,LNW,
     &                  GHONDO,LSTRI,LSTRJ,LSTRK,LSTRL,
     &                  CUTOFF,NINT,HFSCAL,CSCALT,L1,L2)
C
  820 CONTINUE
C
C     ----- END OF SHELL LOOPS -----
C
  860 CONTINUE
  880 CONTINUE
  900 CONTINUE
  920 CONTINUE
  930 CONTINUE
      IF(GOPARR) THEN
         CALL DDI_GSUMF(1000,FA ,NXYZ*L2)
         CALL DDI_GSUMI(1001,NINT  ,1)
         CALL DDI_GSUMI(1002,NSCHWZ,1)
      END IF
C A+B
         CALL DSCAL(NXYZ*L2,HALF,FA,1)
            II=0
            DO I=1,L1
               II = II+I
               DO J=1,NXYZ
                  FA(II,J)  = FA(II,J)  + FA(II,J)
               ENDDO ! J loop
            ENDDO ! I loop
C
      IF(DLB.and.(.not.offload_makefp)) CALL DDI_DLBRESET
C
C     ----- OUTPUT THE LAST BITS OF INTEGRALS -----
C
      GOPARR = GPSAVE
      CUTOFF = CUTSV
      NT = NTTEMP
      RETURN
C
 9010 FORMAT(1X,'II,JST,KST,LST =',4I3,' NREC =',I10,' INTLOC =',I5)
 9015 FORMAT(1X,'II,JST,KST,LST =',4I3,' IN CORE, INTLOC =',I12)
 9020 FORMAT(1X,'SCHWARZ INEQUALITY TEST SKIPPED',I12,
     *        ' INTEGRAL BLOCKS.')
 9030 FORMAT(//1X,'*** THIS JOB HAS EXHAUSTED ITS CPU TIME ***'/
     *         1X,'     (WHILE COMPUTING 2E- INTEGRALS)'///)
      END
C*MODULE SCFLIB  *DECK DIRFCK_CPHFDYN
      SUBROUTINE DIRFCK_CPHFDYN_APB(IA,DA,FA,NXYZ,ISH,JSH,KSH,LSH,
     &                  GHONDO,ISTRIDE,JSTRIDE,KSTRIDE,LSTRIDE,
     &                  CUTOFF,NINT,HFSCAL,CSCALT,LL1,LL2)
      use mx_limits, only: mxsh,mxgtot
      USE CONSTANTS, ONLY: HALF,TWO,FOUR
C
      IMPLICIT NONE
C
      DIMENSION IA(*),GHONDO(*),
     &  DA(LL2,NXYZ), FA(LL2,NXYZ)
C
      LOGICAL IANDJ,KANDL,SAME
C
      DOUBLE PRECISION CD,CF,CG,CH,CI,CP,CS,DA,EX,FA,GHONDO
      DOUBLE PRECISION CSCALT,CUTOFF,HFSCAL
      DOUBLE PRECISION VAL
      DOUBLE PRECISION VAL1,VAL4
C
      INTEGER IA,KATOM,KLOC,KMAX,KMIN,KNG,KSTART,KTYPE,LL2,NSHELL,NXYZ
      INTEGER ISH,ISTRIDE,JSH,JSTRIDE,KSH,KSTRIDE,LL1,LOCI,LOCJ,LOCK
      INTEGER LOCL,LSH,LSTRIDE,MAXI,MAXJ,MAXK,MAXL,MINI,MINJ,MINK,MINL
      INTEGER NINT,I,I1,I2,I_INDEX,II,IJ_INDEX,IJK_INDEX,IJKL_INDEX
      INTEGER ITMP,J,J1,J2,JJ,K,K1,K2,KK,L,L1,L2,MAXJ2,MAXL2,NIJ,NKL
      INTEGER II2,IJ,IK,IL,IXYZ,JJ2,JK,JL,KK2,KL,LL
C
      COMMON /NSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),
     *                CF(MXGTOT),CG(MXGTOT),CH(MXGTOT),CI(MXGTOT),
     *                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),
     *                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
C
      SAME  = ISH.EQ.KSH.AND.JSH.EQ.LSH
      IANDJ = ISH.EQ.JSH
      KANDL = KSH.EQ.LSH
C
      MINI = KMIN(ISH)
      MINJ = KMIN(JSH)
      MINK = KMIN(KSH)
      MINL = KMIN(LSH)
      MAXI = KMAX(ISH)
      MAXJ = KMAX(JSH)
      MAXK = KMAX(KSH)
      MAXL = KMAX(LSH)
      LOCI = KLOC(ISH)-MINI
      LOCJ = KLOC(JSH)-MINJ
      LOCK = KLOC(KSH)-MINK
      LOCL = KLOC(LSH)-MINL
C
      NIJ = 0
      MAXJ2 = MAXJ
      I_INDEX = 1
      DO 360 I = MINI,MAXI
         IF (IANDJ) MAXJ2 = I
C
         I1 = I+LOCI
C
         IJ_INDEX = I_INDEX
         I_INDEX = I_INDEX + ISTRIDE
C
         DO 340 J = MINJ,MAXJ2
            NIJ = NIJ+1
            MAXL2 = MAXL
C
            J1 = J+LOCJ
            I2 = I1
            J2 = J1
            IF (I1.LT.J1) THEN ! SORT <IJ|
               I2 = J1
               J2 = I1
            ENDIF
C
            IJK_INDEX = IJ_INDEX
            IJ_INDEX = IJ_INDEX + JSTRIDE
C
            NKL = NIJ
C
            DO 320 K =  MINK,MAXK
               IF (KANDL) MAXL2 = K
C
               K1 = K + LOCK
C
               IF(SAME) THEN ! ACCOUNT FOR NON-UNIQUE PERMUTATIONS
                  ITMP = MIN(MAXL2-MINL+1,NKL)
                  IF (ITMP.EQ.0) GOTO 340
                  MAXL2 = MINL + ITMP - 1
                  NKL = NKL - ITMP
               ENDIF
C
               IJKL_INDEX = IJK_INDEX
               IJK_INDEX = IJK_INDEX + KSTRIDE
C
               DO 300 L=MINL,MAXL2
C
                  VAL = GHONDO( IJKL_INDEX )
                  IJKL_INDEX = IJKL_INDEX + LSTRIDE
                  IF(ABS(VAL).LT.CUTOFF) GOTO 300
                  NINT = NINT + 1
C
                  L1 = L + LOCL
                  K2 = K1
                  L2 = L1
C
                  IF (K2.LT.L2) THEN ! SORT |KL>
                     K2 = L1
                     L2 = K1
                  ENDIF
C
                  II = I2
                  JJ = J2
                  KK = K2
                  LL = L2
C
                  IF (II.LT.KK) THEN ! SORT <IJ|KL>
                     II = K2
                     JJ = L2
                     KK = I2
                     LL = J2
                  ELSE IF (II.EQ.KK.AND.JJ.LT.LL) THEN ! SORT <IJ|IL>
                     JJ = L2
                     LL = J2
                  ENDIF
C
                  II2 = IA(II)
                  JJ2 = IA(JJ)
                  KK2 = IA(KK)
C
                  IJ = II2 + JJ
                  IK = II2 + KK
                  IL = II2 + LL
                  JK = JJ2 + KK
                  JL = JJ2 + LL
                  KL = KK2 + LL
                  IF (JJ.LT.KK) JK = KK2 + JJ
                  IF (JJ.LT.LL) JL = IA(LL) + JJ
C
C       ACCOUNT FOR IDENTICAL PERMUTATIONS.
C
                  IF(II.EQ.JJ) VAL = VAL*HALF
                  IF(KK.EQ.LL) VAL = VAL*HALF
                  IF(II.EQ.KK.AND.JJ.EQ.LL) VAL = VAL*HALF
                  VAL1 = VAL
                  VAL4 = VAL*FOUR
C A+B
                  DO IXYZ=1,NXYZ
                   FA(IJ,IXYZ) = FA(IJ,IXYZ) + VAL4*DA(KL,IXYZ)
                   FA(KL,IXYZ) = FA(KL,IXYZ) + VAL4*DA(IJ,IXYZ)
                   FA(IK,IXYZ) = FA(IK,IXYZ) - VAL1*DA(JL,IXYZ)
                   FA(JL,IXYZ) = FA(JL,IXYZ) - VAL1*DA(IK,IXYZ)
                   FA(IL,IXYZ) = FA(IL,IXYZ) - VAL1*DA(JK,IXYZ)
                   FA(JK,IXYZ) = FA(JK,IXYZ) - VAL1*DA(IL,IXYZ)
                  ENDDO

 300           ENDDO
 320        ENDDO
 340     ENDDO
 360  ENDDO
C
      RETURN
      END
      SUBROUTINE TWOEI_CPHFDYN_AMB(NINT,NSCHWZ,L1,L2,NINTMX,XINTS,NSH2,
     *                 GHONDO,MAXG,DDIJ,DA,FB,DSH,NXYZ)
      use mx_limits, only: mxsh,mxgtot,mxatm,MXAO
C$    USE params, ONLY: intomp
C$    use TWOEI_CPHFDYN_omp
      use MOD_OFFLOAD, only: offload_makefp
      USE CONSTANTS, ONLY: ZERO,HALF
      USE MAKEFP_CPHF, only: makefp_print_default
C
      IMPLICIT NONE
C
      LOGICAL OUT,SCHWRZ,SCHSKP,GOPARR,DSKWRK,MASWRK,DLB,SLB,C1GRP
      LOGICAL SKIPA,SKIPB,SKIPC,NPSYM
      LOGICAL PK,PANDK,NOTPK,BLOCK,GPSAVE,SCREEN
C
      DIMENSION XINTS(NSH2),
     *          GHONDO(MAXG),
     *          DSH(NSH2),DDIJ(*)
      DIMENSION FB(L2,NXYZ), DA(L2,NXYZ)
      DIMENSION MI(48),MJ(48),MK(48),M0(48)
C
      COMMON /IJPAIR/ IA(MXAO)
      COMMON /INTOPT/ ISCHWZ,IECP,NECP,IEFLD
      COMMON /INT2IC/ NINTIC,ININTIC,NXXIC,LBUFPIC,LIXIC,LABSIX,NINTIX
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON /NSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),
     *                CF(MXGTOT),CG(MXGTOT),CH(MXGTOT),CI(MXGTOT),
     *                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),
     *                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
      COMMON /OUTPUT/ NPRINT,ITOL,ICUT,NORMF,NORMP,NOPK
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /PKFIL / PK,PANDK,BLOCK
      COMMON /RESTAR/ TIMLIM,IREST,NREC,INTLOC,IST,JST,KST,LST
      COMMON /SCINP / VLAMB,SCREEN
      COMMON /SHLEXC/ NORGSH(3),NORGSP(3),IEXCH,NANGM,NGTH(4)
      COMMON /SHLG70/ ISH,JSH,KSH,LSH,IJKLXX(4)
      COMMON /SHLNOS/ QQ4,LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK,LOCL,
     *                MINI,MINJ,MINK,MINL,MAXI,MAXJ,MAXK,MAXL,
     *                NIJ,IJ,KL,IJKL
      COMMON /SHLT  / TOL,CUTOFF,ICOUNT,OUT
      COMMON /SYMTRY/ MAPSHL(MXSH,48),MAPCTR(MXATM,48),
     *                T(432),INVT(48),NT
      COMMON /ERIOUT/ INW,JNW,KNW,LNW,LSTRI,LSTRJ,LSTRK,LSTRL
      COMMON /DFTPAR/ DFTTYP(20),EXENA,EXENB,EXENC,
     *                IDFT34,NAUXFUN,NAUXSHL
      !COMMON /MACHSW/ KDIAG,ICORFL,IXDR,modio,mem10,lpnt10,mem10m
      logical not_gpu, masout
C
      DOUBLE PRECISION DA,DDIJ,DSH,FB,GHONDO,XINTS
      DOUBLE PRECISION CD,CF,CG,CH,CI,CP,CS,EX
      DOUBLE PRECISION TIMLIM,VLAMB
      DOUBLE PRECISION QQ4
      DOUBLE PRECISION CUTOFF,DFTTYP,EXENA,EXENB,EXENC,T,TOL
      DOUBLE PRECISION CSCALT,CUTSV,DENMAX,HFSCAL,TEST,TIM
C
      INTEGER I,IA,IECP,IEFLD,ININTIC,ISCHWZ,L2,LABSIX,LBUFPIC,LIXIC,M0
      INTEGER MAXG,MI,MJ,MK,NECP,NINTIC,NSH2,NXXIC,NXYZ
      INTEGER IDAF,IODA,IP,IPK,IR,IS,IW,KATOM,KLOC,KMAX,KMIN,KNG,KSTART
      INTEGER KTYPE,NAV,NINTIX,NSHELL
      INTEGER IBTYP,ICUT,IEXCH,INTLOC,IPTIM,IREST,IST,ITOL,JST,KST,LST
      INTEGER MASTER,ME,NANGM,NGTH,NOPK,NORGSH,NORGSP,NORMF,NORMP
      INTEGER NPRINT,NPROC,NREC
      INTEGER IJ,IJKLXX,ISH,JSH,KL,KSH,LIT,LJT,LKT,LLT,LOCI,LOCJ,LOCK
      INTEGER LOCL,LSH,MAXI,MAXJ,MAXK,MAXL,MINI,MINJ,MINK,MINL,NIJ
      INTEGER ICOUNT,IDFT34,IJKL,INVT,INW,JNW,KNW,LNW,LSTRI,LSTRJ,LSTRK
      INTEGER LSTRL,MAPCTR,MAPSHL,NAUXFUN,NAUXSHL,NT
      INTEGER II,IJIJ,IPCOUNT,JJ,JORK,KK,KLKL,L1,LL,LMAX,MINE,NEXT,NINT
      INTEGER NINTMX,NSCHWZ,NTTEMP
      INTEGER J
C
      MASOUT=MASWRK.and.makefp_print_default
C
C          ----- TWO-ELECTRON INTEGRALS -----
C     THIS VERSION CAN HANDLE S,P,D,F,G,H,I AND L SHELLS
C
      TIM = ZERO
      CALL TSECND(TIM)
C
      QQ4=1.0D+00
      ICOUNT=  0
      GPSAVE = GOPARR
      CUTSV  = CUTOFF
      CUTOFF = MIN(CUTOFF,1.0D-10)
      NTTEMP = NT
      NT = 1
      HFSCAL=DFTTYP(3)
      CSCALT=1.0D+00
C
C     ----- INITIALIZATION FOR PARALLEL WORK -----
C     BOTH STATIC AND DYNAMIC LOAD BALANCING ARE IMPLEMENTED BELOW
C
      SLB = GOPARR  .AND.  IBTYP.EQ.0
      DLB = GOPARR  .AND.  IBTYP.EQ.1
      NEXT = -1
      MINE = -1
      IPCOUNT = ME - 1
C
      C1GRP = NT.EQ.1
C
      CALL BASCHK(LMAX)
                    NANGM =  4
      IF(LMAX.EQ.2) NANGM =  6
      IF(LMAX.EQ.3) NANGM = 10
      IF(LMAX.EQ.4) NANGM = 15
      IF(LMAX.EQ.5) NANGM = 21
      IF(LMAX.EQ.6) NANGM = 28
      NGTH(4) = 1
      NGTH(3) = NGTH(4) * NANGM
      NGTH(2) = NGTH(3) * NANGM
      NGTH(1) = NGTH(2) * NANGM
      IF(NOPK.EQ.0) THEN
         NORGSH(1) = 0
         NORGSH(2) = NORGSH(1) + NANGM**4
         NORGSH(3) = NORGSH(2) + NANGM**4
         NORGSP(1) = 0
         NORGSP(2) = 256
         NORGSP(3) = 512
      ELSE
         DO I=1,3
            NORGSH(I) = 0
            NORGSP(I) = 0
         ENDDO
      END IF
C
      NOTPK = .NOT.PK
      NINT  = 0
      NSCHWZ= 0
      SCHSKP=.FALSE.
      DENMAX = ZERO
         SCHWRZ = ISCHWZ.EQ.1
         IF(SCHWRZ) THEN
            CALL DAREAD(IDAF,IODA,XINTS,NSH2,54,0)
            CALL AOSHLDxyz(DA,DSH,IA,L1,L2,NSH2,NXYZ)
         END IF
C========================      OPENMP CODE       ======================
!$    if(intomp.gt.0) then
!$    CALL  TWOEI_CPHFDYN_AMB_omp(SCHWRZ,NINT,NSCHWZ,L1,L2,XINTS,
!$   &        NSH2,MAXG,IA,DA,FB,DSH,NXYZ,cutoff,.false.,HFSCAL,CSCALT)
!$    goto 930
!$    endif
C========================       END OPENMP       ======================
C
C        NOW WE ARE READY TO LOOP OVER ALL NSHELL**4 SHELL QUARTETS
C
C     ----- I SHELL -----
C
      DO 920 II = 1,NSHELL
      DO 900 JJ = 1,II
C
C     ----- GO PARALLEL! -----
C
      IF (DLB) THEN
         MINE = MINE + 1
         IF (MINE.GT.NEXT) CALL DDI_DLBNEXT(NEXT)
         IF (NEXT.NE.MINE) GO TO 900
      END IF
      DO 880 KK = 1,II
      jork=kk
      if(ii.eq.kk) jork=jj
C
C     ----- GO PARALLEL! -----
C
      IF(SLB) THEN
         IPCOUNT = IPCOUNT + 1
         IF (MOD(IPCOUNT,NPROC).NE.0) GO TO 880
      END IF
      DO 860 LL = 1,jork
      ISH=II
      JSH=JJ
      KSH=KK
      LSH=LL
C
C     APPLY THE SCHWARZ INEQUALITY TO SCREEN OUT SMALL INTEGRALS,
C     SEE, FOR EXAMPLE, J.L.WHITTEN, J.CHEM.PHYS. 58,4496-4501(1973)
C
      IF(SCHWRZ) THEN
         IJIJ = (ISH*ISH-ISH)/2 + JSH
         KLKL = (KSH*KSH-KSH)/2 + LSH
         TEST = QQ4*XINTS(IJIJ)*XINTS(KLKL)
            !DENMAX = SCHWDN(DSH,ISH,JSH,KSH,LSH,IA)
            !TEST = TEST*DENMAX
         SCHSKP = TEST.LT.CUTOFF
         IF(SCHSKP) NSCHWZ = NSCHWZ + 1
      END IF
      IF(SCHSKP) GO TO 820
C
C        ----- ELECTRON REPULSION INTEGRAL CALCULATION -----
C     THIS MAY USE ROTATED AXIS, ERIC, OR RYS QUADRATURE METHODS
C
C update shell index to ISH,JSH,KSH,LSH in /SHLG70/
      CALL SHELLQUART(ISH,JSH,KSH,LSH,GHONDO)
C need to get shell index from /ERIOUT/
         CALL DIRFCK_CPHFDYN_AMB(IA,DA,FB,NXYZ,INW,JNW,KNW,LNW,
     &                  GHONDO,LSTRI,LSTRJ,LSTRK,LSTRL,
     &                  CUTOFF,NINT,HFSCAL,CSCALT,L1,L2)
C
  820 CONTINUE
C
C     ----- END OF SHELL LOOPS -----
C
  860 CONTINUE
  880 CONTINUE
  900 CONTINUE
  920 CONTINUE
  930 CONTINUE
      IF(GOPARR) THEN
         CALL DDI_GSUMF(1000,FB ,NXYZ*L2)
         CALL DDI_GSUMI(1001,NINT  ,1)
         CALL DDI_GSUMI(1002,NSCHWZ,1)
      END IF

C A-B
         CALL DSCAL(NXYZ*L2,HALF,FB,1)
            II=0
            DO I=1,L1
               II = II+I
               DO J=1,NXYZ
                  FB(II,J)  = FB(II,J)  + FB(II,J)
               ENDDO ! J loop
            ENDDO ! I loop
C
      IF(DLB.and.(.not.offload_makefp)) CALL DDI_DLBRESET
C
C     ----- OUTPUT THE LAST BITS OF INTEGRALS -----
C
      GOPARR = GPSAVE
      CUTOFF = CUTSV
      NT = NTTEMP
      RETURN
C
 9010 FORMAT(1X,'II,JST,KST,LST =',4I3,' NREC =',I10,' INTLOC =',I5)
 9015 FORMAT(1X,'II,JST,KST,LST =',4I3,' IN CORE, INTLOC =',I12)
 9020 FORMAT(1X,'SCHWARZ INEQUALITY TEST SKIPPED',I12,
     *        ' INTEGRAL BLOCKS.')
 9030 FORMAT(//1X,'*** THIS JOB HAS EXHAUSTED ITS CPU TIME ***'/
     *         1X,'     (WHILE COMPUTING 2E- INTEGRALS)'///)
      END
C*MODULE SCFLIB  *DECK DIRFCK_CPHFDYN
      SUBROUTINE DIRFCK_CPHFDYN_AMB(IA,DA,FB,NXYZ,ISH,JSH,KSH,LSH,
     &                  GHONDO,ISTRIDE,JSTRIDE,KSTRIDE,LSTRIDE,
     &                  CUTOFF,NINT,HFSCAL,CSCALT,LL1,LL2)
      use mx_limits, only: mxsh,mxgtot
      USE CONSTANTS, ONLY: HALF,TWO,FOUR
C
      IMPLICIT NONE
C
      DIMENSION IA(*),GHONDO(*),
     &  DA(LL2,NXYZ), FB(LL2,NXYZ)
C
      LOGICAL IANDJ,KANDL,SAME
C
C
      COMMON /NSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),
     *                CF(MXGTOT),CG(MXGTOT),CH(MXGTOT),CI(MXGTOT),
     *                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),
     *                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
C
      DOUBLE PRECISION CD,CF,CG,CH,CI,CP,CS,DA,EX,FB,GHONDO
      DOUBLE PRECISION CSCALT,CUTOFF,HFSCAL
      DOUBLE PRECISION VAL
      DOUBLE PRECISION DIK,DIL,DJK,DJL,VAL1,VAL2
C
      INTEGER IA,KATOM,KLOC,KMAX,KMIN,KNG,KSTART,KTYPE,LL2,NSHELL,NXYZ
      INTEGER ISH,ISTRIDE,JSH,JSTRIDE,KSH,KSTRIDE,LL1,LOCI,LOCJ,LOCK
      INTEGER LOCL,LSH,LSTRIDE,MAXI,MAXJ,MAXK,MAXL,MINI,MINJ,MINK,MINL
      INTEGER NINT
      INTEGER I,I1,I2,I_INDEX,II,IJ_INDEX,IJK_INDEX,IJKL_INDEX,ITMP,J
      INTEGER J1,J2,JJ,K,K1,K2,KK,L,L1,L2,MAXJ2,MAXL2,NIJ,NKL
      INTEGER II2,IJ,IK,IL,IXYZ,JJ2,JK,JL,KK2,KL,LL
C
      SAME  = ISH.EQ.KSH.AND.JSH.EQ.LSH
      IANDJ = ISH.EQ.JSH
      KANDL = KSH.EQ.LSH
C
      MINI = KMIN(ISH)
      MINJ = KMIN(JSH)
      MINK = KMIN(KSH)
      MINL = KMIN(LSH)
      MAXI = KMAX(ISH)
      MAXJ = KMAX(JSH)
      MAXK = KMAX(KSH)
      MAXL = KMAX(LSH)
      LOCI = KLOC(ISH)-MINI
      LOCJ = KLOC(JSH)-MINJ
      LOCK = KLOC(KSH)-MINK
      LOCL = KLOC(LSH)-MINL
C
      NIJ = 0
      MAXJ2 = MAXJ
      I_INDEX = 1
      DO 360 I = MINI,MAXI
         IF (IANDJ) MAXJ2 = I
C
         I1 = I+LOCI
C
         IJ_INDEX = I_INDEX
         I_INDEX = I_INDEX + ISTRIDE
C
         DO 340 J = MINJ,MAXJ2
            NIJ = NIJ+1
            MAXL2 = MAXL
C
            J1 = J+LOCJ
            I2 = I1
            J2 = J1
            IF (I1.LT.J1) THEN ! SORT <IJ|
               I2 = J1
               J2 = I1
            ENDIF
C
            IJK_INDEX = IJ_INDEX
            IJ_INDEX = IJ_INDEX + JSTRIDE
C
            NKL = NIJ
C
            DO 320 K =  MINK,MAXK
               IF (KANDL) MAXL2 = K
C
               K1 = K + LOCK
C
               IF(SAME) THEN ! ACCOUNT FOR NON-UNIQUE PERMUTATIONS
                  ITMP = MIN(MAXL2-MINL+1,NKL)
                  IF (ITMP.EQ.0) GOTO 340
                  MAXL2 = MINL + ITMP - 1
                  NKL = NKL - ITMP
               ENDIF
C
               IJKL_INDEX = IJK_INDEX
               IJK_INDEX = IJK_INDEX + KSTRIDE
C
               DO 300 L=MINL,MAXL2
C
                  VAL = GHONDO( IJKL_INDEX )
                  IJKL_INDEX = IJKL_INDEX + LSTRIDE
                  IF(ABS(VAL).LT.CUTOFF) GOTO 300
                  NINT = NINT + 1
C
                  L1 = L + LOCL
                  K2 = K1
                  L2 = L1
C
                  IF (K2.LT.L2) THEN ! SORT |KL>
                     K2 = L1
                     L2 = K1
                  ENDIF
C
                  II = I2
                  JJ = J2
                  KK = K2
                  LL = L2
C
                  IF (II.LT.KK) THEN ! SORT <IJ|KL>
                     II = K2
                     JJ = L2
                     KK = I2
                     LL = J2
                  ELSE IF (II.EQ.KK.AND.JJ.LT.LL) THEN ! SORT <IJ|IL>
                     JJ = L2
                     LL = J2
                  ENDIF
C
                  II2 = IA(II)
                  JJ2 = IA(JJ)
                  KK2 = IA(KK)
C
                  IJ = II2 + JJ
                  IK = II2 + KK
                  IL = II2 + LL
                  JK = JJ2 + KK
                  JL = JJ2 + LL
                  KL = KK2 + LL
                  IF (JJ.LT.KK) JK = KK2 + JJ
                  IF (JJ.LT.LL) JL = IA(LL) + JJ
C
C       ACCOUNT FOR IDENTICAL PERMUTATIONS.
C
                  IF(II.EQ.JJ) VAL = VAL*HALF
                  IF(KK.EQ.LL) VAL = VAL*HALF
                  IF(II.EQ.KK.AND.JJ.EQ.LL) VAL = VAL*HALF
                  VAL1 = VAL
                  VAL2 = VAL
! A-B
                  DO IXYZ=1,NXYZ
                   DJL = DA(JL,IXYZ)
                   DIK = DA(IK,IXYZ)
                   DJK = DA(JK,IXYZ)
                   DIL = DA(IL,IXYZ)
                IF (JJ.GE.LL) VAL1 = -VAL
                IF (JJ.GE.KK) VAL2 = -VAL
                   FB(IK,IXYZ) = FB(IK,IXYZ) + VAL1*DJL
                   FB(JL,IXYZ) = FB(JL,IXYZ) + VAL1*DIK
                   FB(IL,IXYZ) = FB(IL,IXYZ) + VAL2*DJK
                   FB(JK,IXYZ) = FB(JK,IXYZ) + VAL2*DIL
                  ENDDO

 300           ENDDO
 320        ENDDO
 340     ENDDO
 360  ENDDO
C
      RETURN
      END
!*MODULE CPHF    *DECK AOSHLDxyz
      SUBROUTINE AOSHLDxyz(D,DSH,IA,L1,L2,NSH2,NFO)
      use mx_limits, only: mxsh,mxgtot
      use constants, only: zero
!
      IMPLICIT NONE
!
      integer, intent(in) :: L1,L2,NSH2,NFO, IA(L1)
      double precision, intent(in) :: D(NFO,L2)
      double precision, intent(inout) :: DSH(NSH2)
!
      COMMON /NSHEL / EX(MXGTOT),CS(MXGTOT),CP(MXGTOT),CD(MXGTOT),
     *                CF(MXGTOT),CG(MXGTOT),CH(MXGTOT),CI(MXGTOT),
     *                KSTART(MXSH),KATOM(MXSH),KTYPE(MXSH),KNG(MXSH),
     *                KLOC(MXSH),KMIN(MXSH),KMAX(MXSH),NSHELL
      INTEGER :: kstart,katom,ktype,kng,kloc,kmin,kmax,nshell
      double precision :: ex,cs,cp,cd,cf,cg,ch,ci

      integer :: IJSH,ISH,JSH, MINI,MAXI, MINJ,MAXJ, IFO, I,J,IJ
      double precision :: DMAX
!
!     ----- COMPRESS DENSITY OVER AOS TO DENSITY OVER SHELLS -----
!
      IJSH=0
      DO 240 ISH=1,NSHELL
         MINI = KLOC(ISH)
         MAXI = MINI + KMAX(ISH) - KMIN(ISH)
         DO 230 JSH=1,ISH
            MINJ = KLOC(JSH)
            MAXJ = MINJ + KMAX(JSH) - KMIN(JSH)
            IJSH = IJSH+1
            DMAX = ZERO
            DO 130 I=MINI,MAXI
               IF(ISH.EQ.JSH) MAXJ=I
               DO 120 J=MINJ,MAXJ
                  IJ = IA(I) + J
                  DO 110 IFO=1,NFO
                     IF(ABS(D(IJ,IFO)).GT.DMAX) DMAX = ABS(D(IJ,IFO))
  110             CONTINUE
  120          CONTINUE
  130       CONTINUE
            DSH(IJSH) = DMAX
  230    CONTINUE
  240 CONTINUE
      RETURN
      END
