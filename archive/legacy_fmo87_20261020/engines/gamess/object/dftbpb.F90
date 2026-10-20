!  1 Dec 20 - YN  - new module for PBC
!
MODULE DFTBPB_MOD
!
  IMPLICIT NONE
!
  INTEGER, PUBLIC :: LINDBND,& ! bond index
                     NKPTS,  & ! Number of k-points
                     LATOPT, & ! Option for lattice optimization
                     LATOPA(6), & ! Which cell parameters are optimized.
                     MAXINT, & ! Maximum number of interacting atoms (default=200)
                     L2PBC, &  ! dimension of the vector for non-G-point approx
                     NSPILL    ! Number of cells for spilled atoms 
  ! PBCBOX(1:3,1:3): real-space PBC box (=BOXREAL)
  ! PBCBOX(1:3,4:6): reciprocal-space PBC box (=BOXREC)
  DOUBLE PRECISION, TARGET, PUBLIC :: PBCBOX(1:3,1:6)
  DOUBLE PRECISION, PUBLIC :: STRESS(3,3), & ! stress tensor
                              STRFMO(3,3), & ! stress tensor used in FMO
                              VOLPBC         ! Volume of the PBC box
!
  LOGICAL, PUBLIC :: PERIOD, & ! Flag for PBC
                     PME,    & ! Flag for particle mesh Ewald (not merged)
                     DIPCOR    ! Flag for dipole moments using Berry phase
  INTEGER, PUBLIC, PARAMETER :: MXKPTS=10000 ! Maximum number of k-points
  INTEGER, PUBLIC, POINTER, DIMENSION(:,:,:) :: INDBND => NULL()
  INTEGER, PUBLIC :: NDIMPBC   ! Dimension of PBC (=1,2,3)
!
! The remaining parameters are private
!
  DOUBLE PRECISION EWALPHA, VOLINV,THREW,THREWD,EWAINP
  INTEGER LVKPTS,LVKWGHT,IERFTP,NBODY,NKFFT,NCBSP
  !! For general PBC
  PRIVATE VOLINV, & ! inverse of VOLPBC
          EWALPHA,& ! Alpha parameter for Ewald summation
          LVKPTS, & ! coordinate of the k-point sampling
          LVKWGHT,& ! weight of k-points
!         NDIMPBC,& ! Dimension of PBC (=1,2,3)
          THREW,  & ! convergence threshold of Ewald summation
          THREWD, & ! convergence threshold of Ewald summation of dispersion
          IERFTP, & ! Type of erf and erfc: 0 (slow exact) and 1 (fast approx)
          EWAINP, & ! ewald alpha parmeters in the input
          NBODY     ! nbody used in FMO-DFTB
  !! For PME
  PRIVATE NKFFT,NCBSP
!
  DOUBLE PRECISION, PRIVATE :: KPTS(4,25), & ! Used only input
                               SCALEXP       ! 1/(4*alpha*alpha)
  DOUBLE PRECISION, PRIVATE, POINTER :: &
    BOXREAL(:,:) => NULL(), & ! PBC box in the real space
    BOXREC(:,:)  => NULL(), & ! PBC box in the reciprocal space
    VKPTS(:,:)   => NULL(), & ! positions of k-points (not merged)
    VKWGHT(:)    => NULL(), & ! weights of k-points (not merged)
    BSQ(:)       => NULL()    ! B^2 for PME (may not be merged)
  INTEGER, PRIVATE, POINTER :: &
    IW           => NULL(), &
    IDAF         => NULL(), &
    IODA(:)      => NULL()
  LOGICAL, PRIVATE, POINTER :: &
    MASWRK       => NULL(), &
    GOPARR       => NULL()
!
  !! For FMO
  DOUBLE PRECISION, PRIVATE, POINTER :: &
    FMOC(:,:)    => NULL(), &
    ZREFFMO(:)   => NULL()
  INTEGER, PRIVATE, POINTER :: &
    NFG          => NULL(), &
    MAXNAT       => NULL(), &
    NATFMO       => NULL(), &
    NBDFG        => NULL(), &
    NATFRG(:)    => NULL(), &
    INDFRG(:)    => NULL(), &
    IATFRG(:)    => NULL(), &
    IAGLOB(:)    => NULL(), &
    INDAT(:)     => NULL(), &
    IALOC(:)     => NULL(), &
    ISPEFMO(:)   => NULL(), &
    NBOND(:)     => NULL()
  INTEGER, PRIVATE :: NATFMOB
!
  DOUBLE PRECISION, PRIVATE, PARAMETER :: &
    PI=3.14159265358979312D+00, &
    PISQ=SQRT(PI), &
    ONE_THIRD=1.0D+00/3.0D+00, &
    ANG2AU=1.8897259877D+00
!
CONTAINS
!
!-------------------------------------------------------------------
!*MODULE DFTBPB    *DECK DFTB_PBCINIT
!>
!>     @brief Initialize PBC.
!>
!>     @details Set initial values in PBC.
!>
!>     @author Yoshio Nishimoto
!>
SUBROUTINE DFTB_PBCINIT(NAT,MAXINT_,NDIMPBC_,EWALPHA_,KPTS_,PME_, &
                        NKFFT_,NCBSP_,DIPCOR_,THREW_,THREWD_,IERFTP_)
!
  IMPLICIT NONE
!
  LOGICAL DSKWRK
  LOGICAL, TARGET :: GOPARR_,MASWRK_
  INTEGER, TARGET :: IW_,IDAF_,IODA_
  INTEGER IR,IP,IJK,IJKT,NAV, ME,MASTER,NPROC,IBTYP,IPTIM,I,NPBC,J
  INTEGER IERFTP_,NAT,NDER
  DOUBLE PRECISION THREW_,THREWD_,VLEN,VMIN,DDOT,dist,derfinv
!
  COMMON /IOFILE/ IR,IW_,IP,IJK,IJKT,IDAF_,NAV,IODA_(950)
  COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR_,DSKWRK,MASWRK_
!
  DOUBLE PRECISION, TARGET :: EWALPHA_,KPTS_(100)
  INTEGER,          TARGET :: MAXINT_,NDIMPBC_,NKFFT_,NCBSP_
  LOGICAL,          TARGET :: PME_,DIPCOR_
!
  MASWRK      => MASWRK_
  GOPARR      => GOPARR_
  IW          => IW_
  IDAF        => IDAF_
  IODA        => IODA_
!
  MAXINT  = MAXINT_
  NDIMPBC = NDIMPBC_
  EWALPHA = EWALPHA_
  CALL DCOPY(100,KPTS_,1,KPTS,1)
  PME     = PME_
  NKFFT   = NKFFT_
  NCBSP   = NCBSP_
  DIPCOR  = DIPCOR_
  BOXREAL(1:3,1:3) => PBCBOX(:,1:3)
  BOXREC (1:3,1:3) => PBCBOX(:,4:6)
  THREW   = THREW_
  THREWD  = THREWD_
  IERFTP  = IERFTP_
!
  NBODY = 0
!
  PERIOD = .FALSE.
  DO I = 1, 9
    IF (PBCBOX(I,1).NE.0.0D+00) THEN
      PERIOD = .TRUE.
      EXIT
    ENDIF
  END DO
!
  IF (.NOT.PERIOD) THEN
    NDIMPBC = 0
    RETURN
  END IF
!
  IF (NDIMPBC.LE.0 .OR. NDIMPBC.GE.4) NDIMPBC = 3
  CALL DSCAL(9,ANG2AU,PBCBOX,1)
  IF (MAXINT.EQ.0) MAXINT = 200 !! USE AS DEFAULT
  !! FIND RECIPROCAL LATTICE VECTORS AND COMPUTE VOLPBC
  CALL DFTB_RECVEC
!
  IF(THREWD.LT.0) THEN
!   THREWD = 1.0D-15
!   IF (NDIMPBC.LE.2 .OR. IDFTBD.NE.1) THREWD=1.0D-10
    THREWD = THREW
  ENDIF
!
  !! Jackson, R. A.; Catlow, C. R. A. Mol. Simul. 1998, 1, 207-224.
  !! Eq. (A.30)
  !! Derivatives of the alpha value is not considered, because
  !! the alpha value should not affect the result of the summation,
  !! if the series has converged.
  EWAINP = EWALPHA
  IF (EWALPHA.EQ.0.0D+00) THEN
    EWALPHA=(NAT*PI*PI*PI/(VOLPBC*VOLPBC))**(1.0D+00/6.0D+00)
    IF (NDIMPBC.EQ.2) EWALPHA=EWALPHA**3 !! ?
    IF (NDIMPBC.EQ.1) EWALPHA=1.0D+00/EWALPHA !! ?
    IF (PME.AND.NDIMPBC.EQ.3) THEN
      VMIN = 1.0D+06
      DO I = 1, 3
        VLEN = SQRT(DDOT(3,PBCBOX(1,I),1,PBCBOX(1,I),1))
        IF (VLEN.LE.VMIN) VMIN = VLEN
      END DO
      !! 0.32 is a very rough lower bound of 16U/5 (U: Hubbard value)
      dist = -(log(threw)*2.302585d+00-1.0d+00)*0.32d+00
      ewalpha = derfinv(1.0d+00-threw*0.1d+00*dist)/dist
    END IF
  END IF
  IF (MASWRK) THEN
    WRITE (IW,*)
    WRITE (IW,'(X,I1,"-DIMENSIONAL PBC CALCULATION")') NDIMPBC
    WRITE (IW,'(X,"REAL-SPACE CELL VECTORS (IN AU)")')
    WRITE (IW,'(14X,"X",19X,"Y",19X,"Z")')
    WRITE (IW,'(X,"A1",3F20.10)') (PBCBOX(I,1),I=1,3)
    IF (NDIMPBC.GE.2) WRITE (IW,'(X,"A2",3F20.10)') (PBCBOX(I,2),I=1,3)
    IF (NDIMPBC.GE.3) WRITE (IW,'(X,"A3",3F20.10)') (PBCBOX(I,3),I=1,3)
    WRITE (IW,*)
    WRITE (IW,'(X,"RECIPROCAL-SPACE VECTORS (IN AU^-1)")')
    WRITE (IW,'(14X,"X",19X,"Y",19X,"Z")')
    WRITE (IW,'(X,"B1",3F20.10)') (PBCBOX(I,4),I=1,3)
    IF (NDIMPBC.GE.2) WRITE (IW,'(X,"B2",3F20.10)') (PBCBOX(I,5),I=1,3)
    IF (NDIMPBC.GE.3) WRITE (IW,'(X,"B3",3F20.10)') (PBCBOX(I,6),I=1,3)
    WRITE (IW,*)
    IF (NDIMPBC.EQ.3) THEN
      WRITE (IW,'(X,"VOLUME OF THE UNIT CELL = ", F12.5," AU^3")') VOLPBC
    ELSE IF (NDIMPBC.EQ.2) THEN
      WRITE (IW,'(X,"AREA   OF THE UNIT CELL = ", F12.5," AU^2")') VOLPBC
    ELSE IF (NDIMPBC.EQ.1) THEN
      WRITE (IW,'(X,"LENGTH OF THE UNIT CELL = ", F12.5," AU")') VOLPBC
    END IF
    WRITE (IW,'(X,"OPTIMUM EWALD PARAMETER = ", F12.5  )') EWALPHA
    WRITE (IW,'(X,"EWALD TRUNC. THRESHOLDS = ",2ES9.2E2)') THREW,THREWD
    WRITE (IW,'(X,"ERF TYPE = ",I2)') IERFTP
  END IF
  SCALEXP = 1.0D+00/(4.0D+00*EWALPHA*EWALPHA)
!
  !! CHECK CELL PARAMETERS FOR 1,2-DIMENSIONAL PBC
  IF (NDIMPBC.NE.3) THEN
    NPBC = 0
    DO I = 1, 3
      DO J = 1, 3
        IF (PBCBOX(I,J).NE.0.0D+00) EXIT
      END DO
      IF (J.NE.4) NPBC = NPBC + 1
    END DO
    IF (NPBC.NE.NDIMPBC) THEN
      IF (MASWRK) THEN
        WRITE (IW,'(/,X, &
        & "THERE IS INCONSISTENCY BETWEEN NDPBC AND PBCBOX")')
        WRITE (IW,'(X,"IN ",I1,"-DIMENSIONAL PBC, ",I1, &
        & " AXES HAVE TO BE ZERO,")') NDIMPBC,3-NDIMPBC
        WRITE (IW,'(X,"WHILE THE PRESENT PBCBOX CONTAINS ",I1, &
        & " ZERO AXES.",/)') 3-NPBC
      END IF
      CALL ABRT
    END IF
  END IF
!
  IF (PME) THEN
    IF (NDIMPBC.NE.3) THEN
      IF (MASWRK) WRITE (IW, &
        '(" PME FOR 2D AND 1D HAS NOT BEEN IMPLEMENTED")')
      CALL ABRT
    END IF
    IF (MASWRK) WRITE (IW,'(" PARTICLE MESH EWALD IS EMPLOYED", &
    & " (NKFFT = ",I4,", NCBSP = ",I2,")")') NKFFT,NCBSP
    IF (NKFFT.NE.  2.AND.NKFFT.NE.  4.AND.NKFFT.NE.  8.AND. &
        NKFFT.NE. 16.AND.NKFFT.NE. 32.AND.NKFFT.NE. 64.AND. &
        NKFFT.NE.128.AND.NKFFT.NE.256) THEN
      !! Fixed memory is allocated in DFTB_SKMEM
      IF (MASWRK) WRITE (IW,'(" NKFFT CAN BE 2,4,8,16,32,64,128,256")')
      CALL ABRT
    END IF
  END IF
!
  CALL DERCHK(NDER)
  IF (NDER.GE.1) THEN
    IF (LATOPT.EQ.999 .OR. LATOPT.LT.0) LATOPT = 0
    IF (MASWRK) THEN
      WRITE (IW,'(" LATTICE OPTIMIZATION LEVEL = ",I3, &
      & " | A,B,C,ALPHA,BETA,GAMMA =",6I2)') LATOPT,(LATOPA(i),i=1,6)
    END IF
  END IF
  IF (DIPCOR) LATOPT = LATOPT + 100
  IF (VOLPBC.LE.1.0D-08) THEN
    IF (MASWRK) THEN
      WRITE (IW,*) "(VOLUME/AREA/LENGTH OF THE UNIT CELL) < 0"
      WRITE (IW,*) "CELL PARAMETERS MAY BE LINEARLY DEPENDENT"
    END IF
    CALL ABRT
  END IF
!
  CALL TRPOSQ(PBCBOX(1,1),3)
  CALL TRPOSQ(PBCBOX(1,4),3)
!
END SUBROUTINE DFTB_PBCINIT
!
!-------------------------------------------------------------------
!*MODULE DFTBPB    *DECK DFTB_PBCPTR
!>
!>     @brief Initialize PBC pointers.
!>
!>     @details Set initial pointers in PBC.
!>
!>     @author Yoshio Nishimoto
!>
SUBROUTINE DFTB_PBCPTR(INDBND_,VKPTS_,VKWGHT_,BSQ_,MXATM)
!
  IMPLICIT NONE
!
  INTEGER MXATM
  DOUBLE PRECISION, TARGET :: VKPTS_(*),VKWGHT_(*),BSQ_(*)
  INTEGER, TARGET :: INDBND_(*)
!
  IF (MAXINT.NE.0) THEN
    INDBND(1:3,1:MAXINT+1,1:MXATM) => INDBND_(1:3*(MAXINT+1)*MXATM)
  END IF
!
  IF (PERIOD) THEN
    VKPTS(1:3,1:MXKPTS) => VKPTS_(1:3*MXKPTS)
    VKWGHT(1:MXKPTS) => VKWGHT_(1:MXKPTS)
    IF (PME) THEN
      BSQ(0:NKFFT) => BSQ_(1:NKFFT+1)
    END IF
  END IF
!
  RETURN
!
END SUBROUTINE DFTB_PBCPTR
!
!-------------------------------------------------------------------
!*MODULE DFTBPB    *DECK DFTB_PBCPTR
!>
!>     @brief Initialize FMO/PBC pointers.
!>
!>     @details Set initial pointers in FMO/PBC.
!>
!>     @author Yoshio Nishimoto
!>
SUBROUTINE DFTB_PBCPTR_FMO2(NFG_,MAXNAT_,NATFMO_,NBDFG_,NATFRG_, &
                            INDFRG_,IATFRG_,IAGLOB_,INDAT_,IALOC_, &
                            ISPEFMO_,NBOND_, &
                            FMOC_,ZREF_)
!
  IMPLICIT NONE
!
  INTEGER, TARGET :: NFG_,MAXNAT_,NATFMO_,NBDFG_,NATFRG_(*), &
                     INDFRG_(*),IATFRG_(*),IAGLOB_(*),INDAT_(*), &
                     IALOC_(*),ISPEFMO_(*),NBOND_(*)
! DOUBLE PRECISION, TARGET :: FMOC_(3*NATFMO),ZREF_(*)
  DOUBLE PRECISION, TARGET :: FMOC_(*),ZREF_(*)
!
  NFG     => NFG_
  MAXNAT  => MAXNAT_
  NATFMO  => NATFMO_
  NBDFG   => NBDFG_
  NATFMOB = NATFMO + NBDFG
!
  NATFRG(1:NFG)      => NATFRG_(1:NFG)
  INDFRG(1:NFG)      => INDFRG_(1:NFG)
  IATFRG(1:NATFMOB)  => IATFRG_(1:NATFMOB)
  IAGLOB(1:MAXNAT)   => IAGLOB_(1:MAXNAT)
  INDAT(1:NATFMO)    => INDAT_(1:NATFMO)
  IALOC(1:NATFMOB)   => IALOC_(1:NATFMOB)
  ISPEFMO(1:NATFMO)  => ISPEFMO_(1:NATFMO)
  NBOND(1:NATFMOB)   => NBOND_(1:NATFMOB)
!
  FMOC(1:3,1:NATFMO) => FMOC_(1:3*NATFMO)
  ZREFFMO(1:NATFMOB) => ZREF_(1:NATFMOB)
!
  RETURN
!
END SUBROUTINE DFTB_PBCPTR_FMO2
!
!-------------------------------------------------------------------
!*MODULE DFTBPB    *DECK DFTB_RECVEC
!>
!>     @brief Compute reciprocal vectors.
!>
!>     @details Set reciprocal vectors in PBC.
!>
!>     @author Yoshio Nishimoto
!>
  SUBROUTINE DFTB_RECVEC
!
  IMPLICIT NONE
!
  DOUBLE PRECISION X,ZERO
  INTEGER LS,LU,LV,LWRK,LWRK1,LAST,NEED,I,LWORK,NNN,INFO,LOADFM
  COMMON /FMCOM / X(1)
!
  DOUBLE PRECISION WRK(1),DDOT
!
  IF (NDIMPBC.EQ.3) THEN
    !! b1 = 2pi/V * (a2*a3)
    !! b2 = 2pi/V * (a3*a1)
    !! b3 = 2pi/V * (a1*a2)
    !! First, a2*a3, a3*a1, a1*a2
    CALL VECPRD(BOXREC(1,1),BOXREAL(1,2),BOXREAL(1,3))
    CALL VECPRD(BOXREC(1,2),BOXREAL(1,3),BOXREAL(1,1))
    CALL VECPRD(BOXREC(1,3),BOXREAL(1,1),BOXREAL(1,2))
!
    !! V = a1.(a2*a3)
    VOLPBC = DDOT(3,BOXREAL,1,BOXREC,1)
  ELSE IF (NDIMPBC.GE.2) THEN
    CALL VECPRD(BOXREAL(1,3),BOXREAL(1,1),BOXREAL(1,2))
    VOLPBC = SQRT(BOXREAL(1,3)**2+BOXREAL(2,3)**2+BOXREAL(3,3)**2)
    CALL VCLR(BOXREAL(1,3),1,3)
!
    !! a*b^T = 2\pi
    !! Obtain b = 2*\pi*(a^-1)^T
    CALL VALFM(LOADFM)
    LAST = LOADFM
    CALL DGESVD('A','A',3,2,BOXREAL,3,X(LAST),X(LAST),3, &
                X(LAST),2,WRK,-1,INFO)
    LWORK = INT(WRK(1))
!
    LS    = LOADFM  + 1
    LU    = LS    + 2
    LV    = LU    + 3*3
    LWRK  = LV    + 2*2
    LWRK1 = LWRK  + LWORK
    LAST  = LWRK1 + 3*3
    NEED  = LAST - LOADFM - 1
    CALL GETFM(NEED)
!
    CALL DCOPY(3*2,BOXREAL,1,X(LWRK1),1)
    CALL DGESVD('A','A',3,2,X(LWRK1),3,X(LS),X(LU),3,X(LV),2, &
                X(LWRK),LWORK,INFO)
!
    !! V
    CALL TRPOSQ(X(LV),2)
    !! V * S^-1
    DO I = 1, 2
      CALL DSCAL(2,1.0D+00/X(LS+I-1),X(LV+2*(I-1)),1)
    END DO
    !! WRK1 = (L^0)^-1 = V * S^-1 * U^T
    ZERO=0.0D+00
    CALL DGEMM('N','T',2,3,2,VOLPBC,X(LV),2,X(LU),3,ZERO,X(LWRK1),2)
    BOXREC(1,1) = X(LWRK1+0)
    BOXREC(2,1) = X(LWRK1+2)
    BOXREC(3,1) = X(LWRK1+4)
    BOXREC(1,2) = X(LWRK1+1)
    BOXREC(2,2) = X(LWRK1+3)
    BOXREC(3,2) = X(LWRK1+5)
    CALL RETFM(NEED)
  ELSE IF (NDIMPBC.EQ.1) THEN
    CALL VCLR(BOXREAL(1,2),1,6)
    CALL DCOPY(3,BOXREAL(1,1),1,BOXREC(1,1),1)
    VOLPBC = BOXREAL(1,1)**2+BOXREAL(2,1)**2+BOXREAL(3,1)**2
  END IF
!
  !! COMPLETE RECIPROCAL VECTORS
  NNN = NDIMPBC*3
  CALL DSCAL(NNN,2.0D+00*PI/VOLPBC,BOXREC,1)
!
  VOLPBC = ABS(VOLPBC)
  IF (NDIMPBC.EQ.1) VOLPBC = SQRT(VOLPBC)
  VOLINV = 1.0D+00/VOLPBC
!
END SUBROUTINE DFTB_RECVEC
!
!-----------------------------------------------------------------------
!*MODULE DFTBPB    *DECK DFTB_SETKPTS
!>
!>     @brief Set K-points.
!>
!>     @details Set K-points in PBC.
!>
!>     @author Yoshio Nishimoto
!>
SUBROUTINE DFTB_SETKPTS
!
  IMPLICIT NONE
!
  DOUBLE PRECISION VEC(3)
!
  LOGICAL KEXP,KAUTO,BAND,DUPLIC
  INTEGER IK,IK2,I,J,NKX,NKY,NKZ,NKTMP,IKX,IKY,IKZ,IKPTS,K
  DOUBLE PRECISION VAL,VX1,VY1,VZ1,VX2,VY2,VZ2,DD2,TMP,DASUM,VALMAX
!
  COMPLEX(KIND(0D0)) TMPZ
!
  KEXP  = .FALSE.
  BAND  = .FALSE.
  KAUTO = .FALSE.
!
  !! For some reasons, B vector for PME is constructed here
  !! Eq. (38) in JCC 2016, 37, 2701-2711.
  !! B(m) = conjg(B(-m))
  !! Here computes B(m)*B(-m)
  IF (PME) THEN
    DO IK = 0, NKFFT
      TMPZ = DCMPLX(0.0D+00,0.0D+00)
      DO IK2 = 0, NCBSP-2
        VAL = 2.0D+00*PI*DBLE(IK*IK2)/DBLE(NKFFT)
        TMPZ = TMPZ &
          + DFTB_CBSP(NCBSP,DBLE(IK2+1))*DCMPLX(COS(VAL),SIN(VAL))
      END DO
      TMPZ = 1.0D+00/TMPZ
      VAL = 2.0D+00*PI*DBLE((NCBSP-1)*IK)/DBLE(NKFFT)
      TMPZ = TMPZ*DCMPLX(COS(VAL),SIN(VAL))
      BSQ(IK) = REAL(TMPZ*DCONJG(TMPZ))
    END DO
  END IF
!
  !! should be improved!
  DO I = 1, 100
    IF (KPTS(I,1).NE.0.0D+00) EXIT
  END DO
  IF (I.EQ.101) THEN
    NKPTS = 0
    RETURN
  END IF
!
  DO I = 1, 25
    DO J = 1, 3
      IF (ABS(KPTS(J,I)).GT.1.0D+00) THEN
        if (i.ne.1) write (*,*) "should be error"
        KAUTO = .TRUE.
      ELSE IF (ABS(KPTS(J,I)).LE.1.0D+00) THEN
        IF (KPTS(4,I).GT.1.0D+00 .OR. KPTS(4,I).EQ.0.0D+00) THEN
          BAND = .TRUE.
        ELSE
          KEXP = .TRUE.
        END IF
      ELSE
        write (*,*) "wrong k-point specification"
        call abrt
      END IF
    END DO
  END DO
  if (band) kexp=.false.
!
! NKPTS = 0
! !! If none of them are true, employ Gamma point approximation
! IF (.NOT.KEXP .AND. .NOT.BAND .AND. .NOT.KAUTO) RETURN
!
  !! If both are turned on, use Monkhorst--Pack grid
  IF (KAUTO) THEN
    DUPLIC = .TRUE.
    IF (KPTS(1,1).LT.0.0D+00 .OR. KPTS(2,1).LT.0.0D+00 &
        .OR.KPTS(3,1).LT.0.0D+00) THEN 
      DUPLIC = .FALSE.
      KPTS(1,1) = ABS(KPTS(1,1))
      KPTS(2,1) = ABS(KPTS(2,1))
      KPTS(3,1) = ABS(KPTS(3,1))
    END IF
    NKX = INT(KPTS(1,1))
    NKY = INT(KPTS(2,1))
    NKZ = INT(KPTS(3,1))
    DO IKX = 1, NKX
      DO IKY = 1, NKY
        DO IKZ = 1, NKZ
          VEC(1) = DBLE(2*IKX-NKX-KPTS(4,1))/DBLE(2*NKX)
          VEC(2) = DBLE(2*IKY-NKY-KPTS(1,2))/DBLE(2*NKY)
          VEC(3) = DBLE(2*IKZ-NKZ-KPTS(2,2))/DBLE(2*NKZ)
          VX1 = VEC(1)
          VY1 = VEC(2)
          VZ1 = VEC(3)
          IF (DUPLIC) THEN
            DO IKPTS = 1, NKPTS
              VX2 = VKPTS(1,IKPTS)
              VY2 = VKPTS(2,IKPTS)
              VZ2 = VKPTS(3,IKPTS)
              IF (VX2.EQ.0.5D+00) VX2 = -0.5D+00
              IF (VY2.EQ.0.5D+00) VY2 = -0.5D+00
              IF (VZ2.EQ.0.5D+00) VZ2 = -0.5D+00
              DD2 = VX2*VX2+VY2*VY2+VZ2*VZ2
              IF (VX1.EQ.-VX2.AND.VY1.EQ.-VY2.AND.VZ1.EQ.-VZ2) THEN
                VKWGHT(IKPTS) = VKWGHT(IKPTS) + 1.0D+00
                EXIT
              END IF
            END DO
            IF (IKPTS.EQ.NKPTS+1) THEN
              NKPTS = NKPTS + 1
              VKPTS(1,NKPTS) = VEC(1)
              VKPTS(2,NKPTS) = VEC(2)
              VKPTS(3,NKPTS) = VEC(3)
              VKWGHT(NKPTS) = 1.0D+00
            END IF
          ELSE
            NKPTS = NKPTS + 1
            VKPTS(1,NKPTS) = VEC(1)
            VKPTS(2,NKPTS) = VEC(2)
            VKPTS(3,NKPTS) = VEC(3)
            VKWGHT(NKPTS) = 1.0D+00
          END IF
        END DO
      END DO
    END DO
!
!
!
  ELSE IF (KEXP) THEN
    DO I = 1, 25
      IF (KPTS(1,I).NE.0.0D+00 .OR. KPTS(2,I).NE.0.0D+00 &
      .OR.KPTS(3,I).NE.0.0D+00 .OR. KPTS(4,I).NE.0.0D+00) THEN
        NKPTS = NKPTS + 1
        CALL DCOPY(3,KPTS(1,I),1,VKPTS(1,NKPTS),1)
        VKWGHT(NKPTS) = KPTS(4,I)
      END IF
    END DO
!
!
!
  ELSE IF (BAND) THEN
!   DO I = 1, 25
    DO I = 1, 24
      IF (KPTS(1,I).NE.0.0D+00 .OR. KPTS(2,I).NE.0.0D+00 &
      .OR.KPTS(3,I).NE.0.0D+00 .OR. KPTS(4,I).NE.0.0D+00) THEN
        NKPTS = NKPTS + 1
        CALL DCOPY(3,KPTS(1,I),1,VKPTS(1,NKPTS),1)
        IF (KPTS(1,I+1).NE.0.0D+00 .OR. KPTS(2,I+1).NE.0.0D+00 &
        .OR.KPTS(3,I+1).NE.0.0D+00 .OR. KPTS(4,I+1).NE.0.0D+00) THEN
!
          NKTMP = KPTS(4,I+1)
          VEC(1) = (KPTS(1,I+1)-KPTS(1,I))/DBLE(NKTMP)
          VEC(2) = (KPTS(2,I+1)-KPTS(2,I))/DBLE(NKTMP)
          VEC(3) = (KPTS(3,I+1)-KPTS(3,I))/DBLE(NKTMP)
          DO J = 1, NKTMP-1
            NKPTS = NKPTS + 1
            DO K = 1, 3
              VKPTS(K,NKPTS) = VKPTS(K,NKPTS-1) + VEC(K)!*J
            END DO
          END DO
        END IF
      END IF
    END DO
    DO I = 1, NKPTS
      VKWGHT(I) = 1.0D+00
    END DO
  END IF
!
! !! REMOVE EQUIVALENT K-POINTS
! IF (KAUTO .AND. KPTS(7,1).NE.0.0D+00) THEN
! END IF
!
  !! Normalize the weight
  TMP = DASUM(NKPTS,VKWGHT,1)
  CALL DSCAL(NKPTS,1.0D+00/TMP,VKWGHT,1)
! do i = 1, nkpts
!   vkpts(1,i) = vkpts(1,i) - 0.0001d+00
! end do
!
  IF (MASWRK) THEN
    WRITE (IW,*)
    WRITE (IW,'(" ********************")')
    WRITE (IW,'(" * K-POINT SAMPLING *")')
    WRITE (IW,'(" ********************")')
    WRITE (IW,*)
    IF (KAUTO) THEN
      WRITE (IW,'(" METHOD = MONKHORST--PACK GRID")')
    ELSE IF (KEXP) THEN
      WRITE (IW,'(" METHOD = EXPLICIT K-POINT SPECIFICATION")')
    ELSE IF (BAND) THEN
      WRITE (IW,'(" METHOD = INTERPOLATION (FOR BAND STRUCTURE)")')
    END IF
!
    WRITE (IW,'(" NUMBER OF KPOINTS = ", I4)') NKPTS
    DO I = 1, NKPTS
      WRITE (IW,'(" K(",I4,") = (",3F12.6,") : WEIGHT =",F9.6)') &
             I,(VKPTS(J,I),J=1,3),VKWGHT(I)
    END DO
  END IF
!
  IF (NKPTS.GT.MXKPTS) THEN
    IF (MASWRK) THEN
      WRITE (IW,'(" THE NUMBER OF K-POINTS IS GREATER THAN THE "&
      &"INTERNAL MAXIMUM NUMBER = ",I5)') MXKPTS
      WRITE (IW,'(" INCREASE MXKPTS AND COMPILE AGAIN, OR DECREASE "&
      &"THE NUMBER OF K-POINTS IN THE CALCULATION")') 
    END IF
    CALL ABRT
  END IF
!
  RETURN
!
END SUBROUTINE DFTB_SETKPTS
!
!-----------------------------------------------------------------------
!*MODULE DFTBPB    *DECK DFTB_INTPAIR
!>
!>     @brief Set interacting pairs.
!>
!>     @details Set interacting pairs in PBC.
!>
!>     @author Yoshio Nishimoto
!>
SUBROUTINE DFTB_INTPAIR(C,NAT,ISPE,MAXANG,NAOPBC)
!
  IMPLICIT NONE
!
  INTEGER,PARAMETER :: MXSPE=10
!
  INTEGER NUMREP,LREPINTV,LREPSHORT,LREPCOEFF,NEEDSK,LSKHTAB,LSKSTAB
  INTEGER LSKGRID,LSKSELF,I,NSEQ,IAT,NBND,NAOI,JAT,NAOJ,IXVEC,IYVEC
  INTEGER IZVEC,J,MINI,MAXXTR,MAXYTR,MAXZTR,ITMP,NAOPBC,NAT
  DOUBLE PRECISION REPCUT,SKDIM,SKSPIN,QREFL,HUBBL,QREF,HUBB,SKCUT2
  DOUBLE PRECISION SKCUT,REPMAX,R2,DISTI,TMP,VALMAX
  COMMON /DFTBRE/ REPCUT(MXSPE,MXSPE),NUMREP(MXSPE,MXSPE),      &
                  LREPINTV(MXSPE,MXSPE),LREPSHORT(MXSPE,MXSPE), &
                  LREPCOEFF(MXSPE,MXSPE)
  COMMON /DFTBSK/ SKDIM(MXSPE,MXSPE),SKSPIN(MXSPE),QREFL(3,MXSPE),  &
                  HUBBL(3,MXSPE),QREF(MXSPE),HUBB(MXSPE),SKCUT2,    &
                  NEEDSK,LSKHTAB(MXSPE,MXSPE),LSKSTAB(MXSPE,MXSPE), &
                  LSKGRID(MXSPE,MXSPE),LSKSELF(MXSPE)
!
  DOUBLE PRECISION C(3,NAT)
  INTEGER ISPE(*),MAXANG(*)
!
  DOUBLE PRECISION VEC1(3),VEC2(3),DISTTMP(MAXINT)
!
  !! Check the size of the PBC box
  SKCUT = SQRT(SKCUT2)
  VALMAX = 0.0D+00
  DO I = 1, 3
    IF (ABS(PBCBOX(1,I)).GE.VALMAX) VALMAX = ABS(PBCBOX(1,I))
  END DO
  MAXXTR = FLOOR(SKCUT/VALMAX) + NSPILL
  VALMAX = 0.0D+00
  DO I = 1, 3
    IF (ABS(PBCBOX(2,I)).GE.VALMAX) VALMAX = ABS(PBCBOX(2,I))
  END DO
  MAXYTR = FLOOR(SKCUT/VALMAX) + NSPILL
  VALMAX = 0.0D+00
  DO I = 1, 3
    IF (ABS(PBCBOX(3,I)).GE.VALMAX) VALMAX = ABS(PBCBOX(3,I))
  END DO
  MAXZTR = FLOOR(SKCUT/VALMAX) + NSPILL
  REPMAX = 0.0D+00
  DO I = 1, MXSPE*MXSPE
    IF (REPCUT(I,1).GE.REPMAX) REPMAX = REPCUT(I,1)
  END DO
  REPMAX = REPMAX*REPMAX
  IF (NDIMPBC.LE.0) MAXXTR = 0
  IF (NDIMPBC.LE.1) MAXYTR = 0
  IF (NDIMPBC.LE.2) MAXZTR = 0
  NSEQ = 0
! write (*,*) "skcut2 = ", skcut2
! write (*,*) "repmax = ", repmax
! write (*,*) "max translational vectors"
! write (*,*) maxxtr, maxytr, maxztr
  DO IAT = 1, NAT
    CALL DCOPY(3,C(1,IAT),1,VEC1,1)
    NBND = 0
    NAOI = MAXANG(ISPE(IAT))*MAXANG(ISPE(IAT))
    DO JAT = 1, NAT
      NAOJ = MAXANG(ISPE(JAT))*MAXANG(ISPE(JAT))
      DO IXVEC = -MAXXTR, MAXXTR
        DO IYVEC = -MAXYTR, MAXYTR
          DO IZVEC = -MAXZTR, MAXZTR
            VEC2(1) = C(1,JAT) + PBCBOX(1,1)*IXVEC   &
              + PBCBOX(2,1)*IYVEC + PBCBOX(3,1)*IZVEC
            IF (ABS(VEC1(1)-VEC2(1)).GE.SKCUT) CYCLE
            VEC2(2) = C(2,JAT) + PBCBOX(1,2)*IXVEC   &
              + PBCBOX(2,2)*IYVEC + PBCBOX(3,2)*IZVEC
            IF (ABS(VEC1(2)-VEC2(2)).GE.SKCUT) CYCLE
            VEC2(3) = C(3,JAT) + PBCBOX(1,3)*IXVEC   &
              + PBCBOX(2,3)*IYVEC + PBCBOX(3,3)*IZVEC
            IF (ABS(VEC1(3)-VEC2(3)).GE.SKCUT) CYCLE
            R2 = (VEC1(1)-VEC2(1))*(VEC1(1)-VEC2(1)) &
               + (VEC1(2)-VEC2(2))*(VEC1(2)-VEC2(2)) &
               + (VEC1(3)-VEC2(3))*(VEC1(3)-VEC2(3))
            IF (R2.LE.SKCUT2) THEN
              NBND = NBND + 1
              IF (NBND.GT.MAXINT) THEN
                if (maswrk) write (iw,*) "increase maxint",maxint
                call abrt
              END IF
!             write (*,*) nbnd,r2
!             write (*,*) ixvec,iyvec,izvec
!             INDBND(1,1,IAT) = NBND
              INDBND(1,NBND+1,IAT) = JAT
              INDBND(2,NBND+1,IAT) = (IXVEC+10)*400 &
                                   + (IYVEC+10)*20  &
                                   + (IZVEC+10)*1
              INDBND(3,NBND+1,IAT) = NSEQ
              NSEQ = NSEQ + NAOI*NAOJ
              DISTTMP(NBND) = R2
            END IF
          END DO
        END DO
      END DO
    END DO
    !! REORDER ACCORDING TO THE DISTANCE
    INDBND(1,1,IAT) = NBND
    DO I = 1, NBND-1
      DISTI = DISTTMP(I)
      MINI = I
      DO J = I+1, NBND
        IF (DISTTMP(J).LT.DISTI) THEN
          DISTI = DISTTMP(J)
          MINI = J
        END IF
      END DO
      TMP = DISTTMP(I)
      DISTTMP(I) = DISTTMP(MINI)
      DISTTMP(MINI) = TMP
      ITMP = INDBND(1,I+1,IAT)
      INDBND(1,I+1,IAT) = INDBND(1,MINI+1,IAT)
      INDBND(1,MINI+1,IAT) = ITMP
      ITMP = INDBND(2,I+1,IAT)
      INDBND(2,I+1,IAT) = INDBND(2,MINI+1,IAT)
      INDBND(2,MINI+1,IAT) = ITMP
      ITMP = INDBND(3,I+1,IAT)
      INDBND(3,I+1,IAT) = INDBND(3,MINI+1,IAT)
      INDBND(3,MINI+1,IAT) = ITMP
    END DO
!   write (*,*) "after reordering"
!   write (*,*) "iat = ", iat
!   write (*,*) "iat1 = ", iat1
!   write (*,*) "nbnd = ", nbnd
    DO I = 1, NBND
      IF (DISTTMP(I).GT.REPMAX) EXIT
    END DO
    INDBND(2,1,IAT) = I-1
!   write (*,*) "nbnd2 = ", indbnd(2,1,iat)
!   do i = 1, nbnd
!     write (*,*) i,indbnd(1,i+1,iat),indbnd(2,i+1,iat), &
!                   indbnd(3,i+1,iat),disttmp(i)
!     nnn = indbnd(2,i+1,iat)
!     nseq2 = indbnd(3,i+1,iat)
!     ix = nnn/400-10
!     iy = (nnn-(ix+10)*400)/20-10
!     iz = (nnn-(ix+10)*400-(iy+10)*20)-10
!     write (*,*) ix,iy,iz
!  write (*,*) "nseq2 = ", i,nseq2
!   end do
!   write (*,*) "---"
  END DO
! call abrt
! write (*,*) "nseq = ", nseq
  NAOPBC = NSEQ
! INDBND(3,1,1) = NAOPBC
!
  RETURN
!
  END SUBROUTINE DFTB_INTPAIR
!
!-----------------------------------------------------------------------
!*MODULE DFTBPB    *DECK DFTB_GAMMA_PBC
!>
!>     @brief Calculate Gamma.
!>
!>     @details Calculate Gamma in PBC.
!>
!>     @author Yoshio Nishimoto
!>
SUBROUTINE DFTB_GAMMA_PBC(NDER,NI,KI0,NAT,NSPE,NSHELL,NNN,ISPE,MAXANG, &
                          C,SPE,HUBBL,DAMPXHE,HUBDER,DAMPXH,SRSCC,     &
                          DFTB3,GAMMA2,GAMMA3,TMPGAM,EGRAD,SHIFT,DQ,   &
                          DQES,DQXY,TDDFTB,EPME,HESSIAN)
!
  IMPLICIT NONE
  INTEGER I,ISH,JSH,ISHMAX,J,II,JSHMAX,JJ,NCOLUMN,NPBC,NONPBC,IN2,IP1
  INTEGER IN1,IP2,LOADFM,LQQ,LFQ,LPLAN,LAST,LMUL1A,NEED,NKFFT3,IAT,JSH0
  INTEGER NSTART,NLAST,JAT,KJ,N,IZ,JZ,KK,KI1,KJ1,IP0,IN0,NNN,NFACT
  DOUBLE PRECISION ALPHA,X,VAL,CHARA1,CHARA2,GRAD2,HESS2,DER32,GRAD3A
  DOUBLE PRECISION HESS3A,GRAD3B,HESS3B,DEVMAX,TMPCHG,ZZ,YY,DD,VALS
  DOUBLE PRECISION OLDVAL,GMSERFC,XX,XX2,GAMINC,TMPOLD,XINV,SCAL2
  DOUBLE PRECISION FACT,VALM1,VALM2,VALN1,VALN2,SCOS,SSIN,GAMID
  DOUBLE PRECISION GSUM,GDX,GDYZ,GSUMO,GDOX,GDOYZ,SCALX,SCALYZ,VALG,XI
  DOUBLE PRECISION VSIN,VCOS,GRAD,HESS,GSUM0,GSUM1,GSUM2,GDSUM0,GDSUM1
  DOUBLE PRECISION GDSUM2,SCAL0,SCAL1,V1,V2,BS2,EI
!
!PME  include '/home/nisimoto/lib/fftw/fftw-3.3.7/api/fftw3.f'
!
  DOUBLE PRECISION, PARAMETER :: ZERO=0.0D+00,ONE=1.0D+00, &
                                 TWO=2.0D+00
!
  COMMON /FACTS / FACT(0:100),NFACT
  COMMON /FMCOM / X(1)
!
  INTEGER, INTENT(IN) :: NDER,NI,KI0,NAT,NSPE,NSHELL,ISPE(*),MAXANG(*)
  DOUBLE PRECISION, INTENT(IN) :: C(3,*),SPE(*),HUBBL(3,*),DAMPXHE, &
                                  HUBDER(*),DQ(*),DQES(*),DQXY(*)
  LOGICAL, INTENT(IN) :: DAMPXH,SRSCC,DFTB3,TDDFTB
!
  DOUBLE PRECISION, INTENT(INOUT) :: GAMMA2(*),GAMMA3(*), &
                                     TMPGAM(11,NNN,NNN), &
                                     EGRAD(3,NAT),SHIFT(*), &
                                     HESSIAN(3*NAT,3*NAT)
!
  INTEGER ISHELL(NSHELL),IND(NSPE*3),KI,ISH0,ITRVECX,ITRVECY,ITRVECZ
  INTEGER IGZ,JGZ,ITRVEC,IXVEC,IYVEC,IZVEC,NSP1,NSP2,NSEQ,K,L
  DOUBLE PRECISION VEC(3),VEC1(3),VEC2(3),VECREC(4),GAMTMP(54,2)
  DOUBLE PRECISION HUBA,HUBB,EPME,SCAL,VOL,DIST,HUBDER1,HUBDER2,VAL2
  DOUBLE PRECISION DFTB_GAMMA2F,GMSERF,VAL3A,VAL3B,DISTI,TMP
  DOUBLE PRECISION H/8HH       /
  LOGICAL HDAMP
  COMPLEX(KIND(0D0)) TMPZ
!
!     2D-Ewald
!     The stress tensor implemented here is based on the third ref.
!     - Harris, F. E. Int. J. Quantum Chem. 1998, 68, 385--404.
!     - Grzybowski, A. et al Phys. Rev. B 2000, 61, 6706--6712.
!     - Kawata, M. et al Chem. Phys. Lett. 2001, 340, 157--164.
!
!     1D-Ewald
!     - Porto, M. J. Phys. A: Math. Gen. 2000, 33, 6211--6218.
!       The sign of the last term in Eq. (13) should be opposite.
!     - Tornberg, A.-K. Adv. Comput. Math 2016, 42, 227--248.
!
  IF (NDER.EQ.0.AND.NI.EQ.0) THEN
    IF (SRSCC) THEN
      NCOLUMN = NSHELL
    ELSE
      NCOLUMN = NAT
    END IF
    CALL VCLR(GAMMA2,1,NCOLUMN*(NCOLUMN+1)/2)
    IF (DFTB3) CALL VCLR(GAMMA3,1,NCOLUMN*NCOLUMN)
  ELSE IF (NDER.EQ.1) THEN
    CALL VCLR(TMPGAM,1,11*NNN*NNN)
  ELSE IF (NDER.EQ.2) THEN
    CALL VCLR(TMPGAM,1,11*NNN*NNN)
  END IF
  DD=0
!
  !! PREPARE PARAMETERS
  IND(1) = 0
  DO I = 2, NSPE
    IF (SRSCC) THEN
      IND(I) = IND(I-1) + MAXANG(I-1)
    ELSE
      IND(I) = IND(I-1) + 1
    END IF
  END DO
  ISH = 0
  JSH = 0
  DO I = 1, NAT
    JSH = IND(ISPE(I))
    ISHMAX = 1
    IF (SRSCC) ISHMAX = MAXANG(ISPE(I))
    DO J = 1, ISHMAX
      ISH = ISH + 1
      ISHELL(ISH) = JSH + J
    END DO
  END DO
  ISH = 0
  DO I = 1, NSPE
    ISHMAX = 1
    IF (SRSCC) ISHMAX = MAXANG(I)
    DO II = 1, ISHMAX
      ISH = ISH + 1
      HUBA = HUBBL(II,I)
      JSH = 0
      DO J = 1, NSPE
        JSHMAX = 1
        IF (SRSCC) JSHMAX = MAXANG(J)
        DO JJ = 1, JSHMAX
          JSH = JSH + 1
          HUBB = HUBBL(JJ,J)
          IF (NDER.EQ.0) THEN
            CALL DFTB_PREP_GAMMA3(HUBA,HUBB,TMPGAM(1,ISH,JSH))
          ELSE IF (NDER.EQ.1) THEN
            CALL DFTB_PREP_GAMMA_GRAD(HUBA,HUBB,TMPGAM(1,ISH,JSH), &
                 DFTB3,DAMPXH)
          ELSE IF (NDER.EQ.2) THEN
            CALL DFTB_PREP_GAMMA_HESS(HUBA,HUBB,TMPGAM(1,ISH,JSH), &
                 DFTB3,DAMPXH)
          END IF
        END DO
      END DO
    END DO
  END DO
!
  !! IDENTIFY PERIODIC AND NON-PERIODIC AXES
  IF (NDIMPBC.NE.3) THEN
    NPBC = 0
    NONPBC = 0
    DO I = 1, 3
      DO J = 1, 3
        IF (BOXREC(J,I).NE.ZERO) EXIT
      END DO
      IF (J.EQ.4) THEN
        NONPBC = NONPBC + 1
        IF (NDIMPBC.EQ.2) THEN
          IN1 = I
        ELSE IF (NDIMPBC.EQ.1) THEN
          IF (NONPBC.EQ.1) IN1 = I
          IF (NONPBC.EQ.2) IN2 = I
        END IF
      ELSE
        NPBC = NPBC + 1
        IF (NDIMPBC.EQ.2) THEN
          IF (NPBC.EQ.1) IP1 = I
          IF (NPBC.EQ.2) IP2 = I
        ELSE IF (NDIMPBC.EQ.1) THEN
          IP1 = I
        END IF
      END IF
    END DO
    IF (NDIMPBC.EQ.1) CALL INIFAC(100)
!   write (*,*) "npbc,nonpbc= ", npbc,nonpbc
!   write (*,*) "ip1,ip2 = ",ip1, ip2
!   write (*,*) "in1,in2 = ",in1, in2
!   if (npbc.ne.ndimpbc) then
!     write (iw,*) "inconsistency of unit cell and ndimpbc"
!     call abrt
!   end if
  END IF
!
! PREPARE PARTICLE MESH EWALD
!
  ALPHA = EWALPHA
  IF (PME.AND.(NI.NE.0.OR.NDER.NE.0)) THEN
    CALL VALFM(LOADFM)
    LQQ   = LOADFM + 1
    LFQ   = LQQ    + NKFFT*NKFFT*NKFFT
    LPLAN = LFQ    + 2*(NKFFT/2+1)*NKFFT*NKFFT
    LAST  = LPLAN  + 8
!
    LMUL1A= LAST
    LAST  = LMUL1A + NAT
!
    NEED  = LAST   - LOADFM - 1
    CALL GETFM(NEED)
    NKFFT3 = NKFFT*NKFFT*NKFFT
!
    CALL VCLR(X(LQQ),1,NKFFT3)
    IF (NFG.NE.0.AND.IODA(556).NE.-1) THEN
    ! CALL DAREAD(IDAF,IODA,X(LMUL1A),NAT,556,0)
    ! do i = 1, nat
    ! write (*,'(i3,f20.10)') i,x(lmul1a+i-1)
    ! end do
!   ! call daxpy(nat,-1.0d+00,zref,1,x(lmul1a),1)
    ! x(lmul1a+1-1) = x(lmul1a+1-1) - 6.0d+00
    ! x(lmul1a+2-1) = x(lmul1a+2-1) - 1.0d+00
    ! x(lmul1a+3-1) = x(lmul1a+3-1) - 1.0d+00
    ! CALL DAXPY(NAT,-1.0D+00,X(LMUL1A),1,DQ,1)
    END IF
    !! Eq. (42)
    DO IAT = 1, NAT
      CALL DFTB_PMEQ(1,NKFFT,NCBSP,BOXREC,C(1,IAT),DQ(IAT), &
                     X(LQQ),VAL,VEC)
    END DO
!
    !! Q -> F(Q)
!PME  CALL DFFTW_PLAN_DFT_R2C_3D(X(LPLAN),NKFFT,NKFFT,NKFFT,X(LQQ), &
!PME                             X(LFQ),FFTW_ESTIMATE)
!PME  CALL DFFTW_EXECUTE_DFT_R2C(X(LPLAN),X(LQQ),X(LFQ))
!PME  CALL DFFTW_DESTROY_PLAN(X(LPLAN))
!
    !! Evaluate theta = B*C
    CALL DFTB_PMEBC(NDER,X(LFQ))
!
    !! F^-1(theta*Q)
!PME  CALL DFFTW_PLAN_DFT_C2R_3D(X(LPLAN),NKFFT,NKFFT,NKFFT,X(LFQ), &
!PME                             X(LQQ),FFTW_ESTIMATE)
!PME  CALL DFFTW_EXECUTE_DFT_C2R(X(LPLAN),X(LFQ),X(LQQ))
!PME  CALL DFFTW_DESTROY_PLAN(X(LPLAN))
!
    IF (NFG.NE.0.AND.IODA(556).NE.-1) THEN
    ! CALL DAXPY(NAT,+1.0D+00,X(LMUL1A),1,DQ,1)
    END IF
    IF (NDER.EQ.0) THEN
      EPME = ZERO
      DO IAT = 1, NAT
        CALL DFTB_PMEQ(2,NKFFT,NCBSP,BOXREC,C(1,IAT),ONE,X(LQQ),VAL,VEC)
        SHIFT(IAT) = SHIFT(IAT) + VAL*4.0D+00*PI*VOLINV
        EPME = EPME + VAL*4.0D+00*PI*VOLINV * DQ(IAT)
      END DO
      EPME = EPME * 0.5d+00
    ELSE IF (NDER.EQ.1) THEN
      DO IAT = 1, NAT
        SCAL = 4.0D+00*PI*VOLINV * DQ(IAT)
        CALL DFTB_PMEQ(3,NKFFT,NCBSP,BOXREC,C(1,IAT),ONE,X(LQQ),VAL,VEC)
        DO KI = 1, 3
          EGRAD(KI,IAT) = EGRAD(KI,IAT) + SCAL*VEC(KI)
        END DO
      END DO
    END IF
!
    CALL RETFM(NEED)
    IF (NDER.EQ.0) RETURN
  END IF
!
  VOL = VOLPBC
  ALPHA = EWALPHA
!
  ISH0 = 0
  ITRVECX = 0
  ITRVECY = 0
  ITRVECZ = 0
  DO IAT = 1, NAT
    CHARA1 = SPE(ISPE(IAT))
    JSH0 = 0
    VEC1(1) = C(1,IAT)
    VEC1(2) = C(2,IAT)
    VEC1(3) = C(3,IAT)
    IF (NI.EQ.0 .AND. KI0.EQ.0) THEN
      NSTART = 1
      NLAST  = IAT
    ELSE
      NSTART = NI
      NLAST  = NI
      DO JAT = 1, NSTART-1
        IF (SRSCC) THEN
          JSHMAX = MAXANG(ISPE(JAT))
        ELSE
          JSHMAX = 1
        END IF
        JSH0 = JSH0 + JSHMAX
      END DO
    END IF
    DO JAT = NSTART, NLAST
      CHARA2 = SPE(ISPE(JAT))
      IF (DAMPXH.AND.(CHARA1.EQ.H.OR.CHARA2.EQ.H)) THEN
        HDAMP = .TRUE.
        CALL DFTB_DMPCHK(ISPE(IAT),ISPE(JAT),HDAMP)
      ELSE
        HDAMP = .FALSE.
      END IF
      IF (SRSCC) THEN
        ISHMAX = MAXANG(ISPE(IAT))
        JSHMAX = MAXANG(ISPE(JAT))
      ELSE
        ISHMAX = 1
        JSHMAX = 1
      END IF
      VEC2(1) = C(1,JAT)
      VEC2(2) = C(2,JAT)
      VEC2(3) = C(3,JAT)
      IGZ = 3*(IAT-1)
      JGZ = 3*(JAT-1)
!
!     ---- SHORT-RANGE + REAL-SPACE CONTRIBUTIONS -----
!
      CALL VCLR(GAMTMP,1,54)
      DO ITRVEC = 0, 100
        IF (NDIMPBC.GE.1) ITRVECX = ITRVEC
        IF (NDIMPBC.GE.2) ITRVECY = ITRVEC
        IF (NDIMPBC.GE.3) ITRVECZ = ITRVEC
        DO IXVEC = -ITRVECX, ITRVECX
        DO IYVEC = -ITRVECY, ITRVECY
        DO IZVEC = -ITRVECZ, ITRVECZ
        IF (ABS(IXVEC).NE.ITRVEC .AND. ABS(IYVEC).NE.ITRVEC .AND. &
            ABS(IZVEC).NE.ITRVEC) CYCLE
!
        IF (PERIOD) THEN
          VEC(1) = VEC1(1) - VEC2(1)  - BOXREAL(1,1)*IXVEC &
                 - BOXREAL(2,1)*IYVEC - BOXREAL(3,1)*IZVEC
          VEC(2) = VEC1(2) - VEC2(2)  - BOXREAL(1,2)*IXVEC &
                 - BOXREAL(2,2)*IYVEC - BOXREAL(3,2)*IZVEC
          VEC(3) = VEC1(3) - VEC2(3)  - BOXREAL(1,3)*IXVEC &
                 - BOXREAL(2,3)*IYVEC - BOXREAL(3,3)*IZVEC
          DIST = SQRT(VEC(1)*VEC(1)+VEC(2)*VEC(2)+VEC(3)*VEC(3))
        ELSE
!         CALL DFTB_CNVSQ(IAT,JAT,NSEQ)
!         DIST = DISTMAT(NSEQ)
        END IF
        IF (NDER.EQ.1 .AND. DIST.LT.1.0D-16) CYCLE
        IF (NDER.EQ.2 .AND. DIST.LT.1.0D-16) CYCLE
!
        DO II = 1, ISHMAX
          ISH = ISH0 + II
          HUBA = HUBBL(II,ISPE(IAT))
          NSP1 = ISHELL(ISH)
          HUBDER1 = HUBDER(NSP1)
          DO JJ = 1, JSHMAX
            JSH = JSH0 + JJ
            HUBB = HUBBL(JJ,ISPE(JAT))
            NSP2 = ISHELL(JSH)
            HUBDER2 = HUBDER(NSP2)
            IF (NDER.EQ.0) THEN
              VAL2 = DFTB_GAMMA2F(DIST,HUBA,HUBB,HUBA*3.2D+00, &
                                  HUBB*3.2D+00,DAMPXHE,HDAMP)
              CALL DFTB_CNVSQ(ISH,JSH,NSEQ)
              GAMMA2(NSEQ) = GAMMA2(NSEQ) + VAL2
              IF (PERIOD .AND. DIST.NE.0.0D+00) THEN
                !! 1/R IS ALREADY ADDED IN DFTB_GAMMA2F
                GAMMA2(NSEQ) = GAMMA2(NSEQ) - GMSERF(ALPHA*DIST)/DIST
              END IF
              IF (DFTB3) THEN
                !! NO 1/R TERM IN DFTB3 GAMMA, SO ONLY SHORT-RANGE
                !! CONTRIBUTIONS ARE CONSIDERED.
                CALL DFTB_GAMMA3F(VAL3A,VAL3B, &
                     HUBA,HUBB,HUBDER1,HUBDER2,DIST,DAMPXHE, &
                     TMPGAM(1,NSP1,NSP2),TMPGAM(1,NSP2,NSP1),HDAMP)
                GAMMA3(ISH+NCOLUMN*(JSH-1)) &
                  = GAMMA3(ISH+NCOLUMN*(JSH-1))+VAL3A
                IF (ISH.NE.JSH) &
                GAMMA3(JSH+NCOLUMN*(ISH-1)) &
                  = GAMMA3(JSH+NCOLUMN*(ISH-1))+VAL3B
              END IF
            ELSE IF (NDER.EQ.1) THEN
              DISTI = 1.0D+00/DIST
              CALL DFTB_GAMDERF(VAL2,VAL3A,HUBA,HUBB,HUBDER1, &
                   DIST,DAMPXHE,TMPGAM(1,NSP1,NSP2),DFTB3,HDAMP)
              VAL = VAL2
!     write (*,*) "iat,jat = ", iat,jat
!  write (*,*) val2, - (2.0D+00*ALPHA/SQRT(PI) &
!                  *EXP(-ALPHA*ALPHA*DIST*DIST) &
!         - GMSERF(ALPHA*DIST)*DISTI)*DISTI
              IF (PERIOD) THEN
                VAL = VAL - (2.0D+00*ALPHA/SQRT(PI) &
                              *EXP(-ALPHA*ALPHA*DIST*DIST) &
                     - GMSERF(ALPHA*DIST)*DISTI)*DISTI
              END IF
              val2=val
              IF (DFTB3) THEN
                CALL DFTB_GAMDERF(TMP,VAL3B,HUBB,HUBA,HUBDER2,DIST, &
                     DAMPXHE,TMPGAM(1,NSP2,NSP1),DFTB3,HDAMP)
                VAL = VAL + DQ(ISH)*VAL3A + DQ(JSH)*VAL3B
              END IF
              IF (NI.EQ.0 .AND. KI0.EQ.0) THEN
                VAL = VAL*DQ(ISH)*DQ(JSH)
                IF (TDDFTB) THEN
                  VAL = VAL &
                      + VAL2*(DQES(ISH)*DQ(JSH) + DQ(ISH)*DQES(JSH) &
                      + 2.0D+00*DQXY(ISH)*DQXY(JSH))
                  IF (DFTB3) THEN
                    VAL = VAL                                        &
                        + VAL3A*DQES(ISH)*DQ(ISH)*DQ(JSH)*2.0D+00    &
                        + VAL3B*DQES(ISH)*DQ(JSH)*DQ(JSH)            &
                        + VAL3B*DQES(JSH)*DQ(JSH)*DQ(ISH)*2.0D+00    &
                        + VAL3A*DQES(JSH)*DQ(ISH)*DQ(ISH)            &
                        +(VAL3A*DQ(ISH)*DQXY(ISH)*DQXY(JSH)*2.0D+00  &
                        + VAL3B*DQ(JSH)*DQXY(ISH)*DQXY(JSH)*2.0D+00  &
                        + VAL3A*DQ(JSH)*DQXY(ISH)*DQXY(JSH)          &
                        + VAL3B*DQ(ISH)*DQXY(JSH)*DQXY(JSH))*2.0D+00
                  END IF
                END IF
                VAL = VAL*DISTI
                SCAL = 2.0D+00
                !! does this occur?
                IF (IAT.EQ.JAT .AND. ISH.EQ.JSH) SCAL = 1.0D+00
                DO K = 1, 3
                  EGRAD(K,IAT) = EGRAD(K,IAT) + VEC(K)*VAL
                  EGRAD(K,JAT) = EGRAD(K,JAT) - VEC(K)*VAL
                  DO L = 1, 3
                    STRESS(K,L) = STRESS(K,L) + VAL*VEC(K)*VEC(L)*SCAL
                  END DO
                END DO
              ELSE
                SHIFT(JSH) = SHIFT(JSH) - VEC(KI0)*VAL*DISTI*DQ(ISH)
                SHIFT(ISH) = SHIFT(ISH) - VEC(KI0)*VAL*DISTI*DQ(JSH)
                !! Some missing terms for DFTB3
                IF (DFTB3) THEN
                  SHIFT(JSH) = SHIFT(JSH) &
                    - VEC(KI0)*DISTI*VAL3B*DQ(JSH)*DQ(ISH)
                  SHIFT(ISH) = SHIFT(ISH) &
                    - VEC(KI0)*DISTI*VAL3A*DQ(ISH)*DQ(JSH)
                END IF
              END IF
            ELSE IF (NDER.EQ.2) THEN
              DISTI = ONE/DIST
              CALL DFTB_GAMHESSF(2,GRAD2,HESS2,DER32,GRAD3A,HESS3A, &
                                 TMP,HUBA,HUBB,HUBDER1,DIST,        &
                                 DAMPXHE,TMPGAM(1,NSP1,NSP2),DFTB3, &
                                 HDAMP)
              IF (PERIOD) THEN
                HESS2 = HESS2 - TWO*(GMSERF(ALPHA*DIST)/DIST          &
                   -TWO*ALPHA/SQRT(PI)*EXP(-ALPHA*ALPHA*DIST*DIST) &
                   *(ALPHA*ALPHA*DIST*DIST+1.0D+00))/(DIST*DIST)
                GRAD2 = GRAD2 &
                 - (TWO*ALPHA/SQRT(PI)*EXP(-ALPHA*ALPHA*DIST*DIST) &
                    -GMSERF(ALPHA*DIST)/DIST)/DIST
              END IF
              IF (DFTB3) THEN
                CALL DFTB_GAMHESSF(2,TMP,TMP,TMP,GRAD3B,HESS3B, &
                                   TMP,HUBB,HUBA,HUBDER2,DIST,  &
                                   DAMPXHE,TMPGAM(1,NSP2,NSP1), &
                                   DFTB3,HDAMP)
                GRAD2= GRAD2+(GRAD3A*DQ(ISH)+GRAD3B*DQ(JSH))/3.0D+00
                HESS2= HESS2+(HESS3A*DQ(ISH)+HESS3B*DQ(JSH))/3.0D+00
              END IF
              DO KI = 1, 3
                DO KJ = 1, 3
               VAL = VEC(KI)*VEC(KJ)*(HESS2-DISTI*GRAD2)*DISTI*DISTI
                  IF (KI.EQ.KJ) VAL = VAL + DISTI*GRAD2
                  VAL = VAL*DQ(ISH)*DQ(JSH)
                  HESSIAN(IGZ+KI,IGZ+KJ)= HESSIAN(IGZ+KI,IGZ+KJ)+VAL
                  HESSIAN(JGZ+KI,JGZ+KJ)= HESSIAN(JGZ+KI,JGZ+KJ)+VAL
                  HESSIAN(IGZ+KI,JGZ+KJ)= HESSIAN(IGZ+KI,JGZ+KJ)-VAL
                  HESSIAN(JGZ+KI,IGZ+KJ)= HESSIAN(JGZ+KI,IGZ+KJ)-VAL
                END DO
              END DO
            END IF
            IF (IAT.EQ.JAT.AND.II.EQ.JJ) EXIT
          END DO
        END DO
!
        END DO
        END DO
        END DO
!       write (*,*) "itrvec = ", itrvec, gamma2(nseq)
!
        !! CHECK THE CONVERGENCE
        N = 0
        IF (NDER.EQ.0) THEN
          DO II = 1, ISHMAX
            ISH = ISH0 + II
            DO JJ = 1, JSHMAX
              JSH = JSH0 + JJ
              CALL DFTB_CNVSQ(ISH,JSH,NSEQ)
              N = N + 1
              GAMTMP(N,2) = GAMMA2(NSEQ)
              IF (DFTB3) THEN
                N = N + 1
                GAMTMP(N,2) = GAMMA3(ISH+NCOLUMN*(JSH-1))
                N = N + 1
                GAMTMP(N,2) = GAMMA3(JSH+NCOLUMN*(ISH-1))
              END IF
              IF (IAT.EQ.JAT.AND.II.EQ.JJ) EXIT
            END DO
          END DO
        ELSE IF (NDER.EQ.1) THEN
          IF (NI.EQ.0 .AND. KI0.EQ.0) THEN
            N = 15
            CALL DCOPY(3,EGRAD(1,IAT),1,GAMTMP(1,2),1)
            CALL DCOPY(3,EGRAD(1,JAT),1,GAMTMP(4,2),1)
            CALL DCOPY(9,STRESS,1,GAMTMP(7,2),1)
          ELSE
            N = 2
            GAMTMP(1,2) = SHIFT(ISH)
            GAMTMP(2,2) = SHIFT(JSH)
          END IF
        ELSE IF (NDER.EQ.2) THEN
          N = 1
          DO II = 1, 3
            IF (II.EQ.1) THEN
              IZ = IGZ
              JZ = IGZ
            ELSE IF (II.EQ.2) THEN
              IZ = JGZ
              JZ = JGZ
            ELSE IF (II.EQ.3) THEN
              IZ = IGZ
              JZ = JGZ
            END IF
            DO KI = 1, 3
              DO KJ = 1, 3
                GAMTMP(N,2) = HESSIAN(IZ+KI,JZ+KJ)
                N=N+1
              END DO
            END DO
          END DO
        END IF
        DEVMAX = 0.0D+00
        DO I = 1, N
          VAL = ABS((GAMTMP(I,1)-GAMTMP(I,2)))
          IF (VAL.GE.DEVMAX) DEVMAX = VAL
        END DO
        IF (DEVMAX.LE.THREW.AND.ITRVEC.GE.NSPILL) EXIT
!       IF (DEVMAX.LE.THREW.AND.ITRVEC.GT.0) EXIT
        CALL DCOPY(N,GAMTMP(1,2),1,GAMTMP(1,1),1)
      END DO
!     write (*,*) "converged short-range at", itrvec
!
!     ---- LONG-RANGE (RECIPROCAL) AND SELF-ENERGY CONTRIBUTIONS ---
!
      IF (PERIOD) THEN
        TMPCHG = ZERO
        DO II = 1, ISHMAX
          ISH = ISH0 + II
          DO JJ = 1, JSHMAX
            JSH = JSH0 + JJ
            TMPCHG = TMPCHG + DQ(ISH)*DQ(JSH)
            IF (TDDFTB) THEN
              TMPCHG = TMPCHG &
                  + DQES(ISH)*DQ(JSH) + DQ(ISH)*DQES(JSH) &
                  + 2.0D+00*DQXY(ISH)*DQXY(JSH)
            END IF
          END DO
        END DO
        IF (PME) THEN
          VAL = 0.0D+00
          VEC(1) = 0.0D+00
          VEC(2) = 0.0D+00
          VEC(3) = 0.0D+00
        ELSE
          IF (NDIMPBC.EQ.2) THEN
            ZZ = VEC1(IN1) - VEC2(IN1)
          ELSE IF (NDIMPBC.EQ.1) THEN
            YY = VEC1(IN1) - VEC2(IN1)
            ZZ = VEC1(IN2) - VEC2(IN2)
            DD = YY*YY + ZZ*ZZ
          END IF
          VAL = 0.0D+00
          VALS= 0.0D+00
          IF (NDER.EQ.1) THEN
            CALL VCLR(VEC,1,3)
            CALL VCLR(GAMTMP(1,1),1,3)
            CALL VCLR(GAMTMP(1,2),1,3)
          END IF
          OLDVAL = 0.0D+00
          DO ITRVEC = 1, 100
            IF (NDIMPBC.GE.1) ITRVECX = ITRVEC
            IF (NDIMPBC.GE.2) ITRVECY = ITRVEC
            IF (NDIMPBC.GE.3) ITRVECZ = ITRVEC
            DO IXVEC = -ITRVECX, ITRVECX
            DO IYVEC = -ITRVECY, ITRVECY
            DO IZVEC = -ITRVECZ, ITRVECZ
            IF (ABS(IXVEC).NE.ITRVEC .AND. ABS(IYVEC).NE.ITRVEC .AND. &
                ABS(IZVEC).NE.ITRVEC) CYCLE
!       
            VECREC(1) = BOXREC(1,1)*IXVEC &
                      + BOXREC(2,1)*IYVEC &
                      + BOXREC(3,1)*IZVEC
            VECREC(2) = BOXREC(1,2)*IXVEC &
                      + BOXREC(2,2)*IYVEC &
                      + BOXREC(3,2)*IZVEC
            VECREC(3) = BOXREC(1,3)*IXVEC &
                      + BOXREC(2,3)*IYVEC &
                      + BOXREC(3,3)*IZVEC
            VECREC(4) = VECREC(1)**2 + VECREC(2)**2 + VECREC(3)**2
            !! i*sin disappears
            IF (NDER.EQ.0) THEN
              IF (NDIMPBC.EQ.3) THEN
                VAL = VAL + EXP(-SCALEXP*VECREC(4))/VECREC(4) &
                    * COS(VECREC(1)*(VEC1(1)-VEC2(1))         &
                         +VECREC(2)*(VEC1(2)-VEC2(2))         &
                         +VECREC(3)*(VEC1(3)-VEC2(3)))
              ELSE IF (NDIMPBC.EQ.2) THEN
                VECREC(4) = SQRT(VECREC(4))
                SCAL= COS(VECREC(IP1)*(VEC1(IP1)-VEC2(IP1))   &
                         +VECREC(IP2)*(VEC1(IP2)-VEC2(IP2)))
                VAL = VAL + SCAL/VECREC(4)                    &
                    * (EXP( VECREC(4)*ZZ)                     &
                      *GMSERFC(VECREC(4)/(TWO*ALPHA)+ALPHA*ZZ)   &
                    +  EXP(-VECREC(4)*ZZ)                     &
                      *GMSERFC(VECREC(4)/(TWO*ALPHA)-ALPHA*ZZ))
              ELSE IF (NDIMPBC.EQ.1) THEN
                XX = VECREC(4)/(4.0D+00*ALPHA*ALPHA)
                XX2= -0.25D+00*VECREC(4)*DD
                CALL CALCEI(XX,GAMINC,2)
                IF (DD.LE.1.0D-08) THEN
                  TMP = GAMINC
                ELSE
                  !! kappa = 0
                  TMP = GAMINC
                  SCAL = ONE
                  TMPOLD = TMP
                  !! kappa = 1
                  SCAL = XX2*SCAL
                  XINV = ONE/XX
                  SCAL2= EXP(-XX)*XINV
                  GAMINC = SCAL2 - GAMINC
                  DO KK = 1, 100
                    TMP = TMP + SCAL*GAMINC/FACT(KK)
                    IF (ABS(TMP-TMPOLD).LE.1.0D-15) EXIT
                    SCAL = XX2*SCAL
                    SCAL2= SCAL2*XINV
                    GAMINC = -(GAMINC-SCAL2)/(KK+1)
                    TMPOLD = TMP
                  END DO
                END IF
                VAL = VAL + TMP*COS(VECREC(IP1)*(VEC1(IP1)-VEC2(IP1)))
              END IF
            ELSE IF (NDER.EQ.1) THEN
              IF (NDIMPBC.EQ.3) THEN
                VAL = EXP(-SCALEXP*VECREC(4))/VECREC(4) &
                    * SIN(VECREC(1)*(VEC1(1)-VEC2(1))   &
                         +VECREC(2)*(VEC1(2)-VEC2(2))   &
                         +VECREC(3)*(VEC1(3)-VEC2(3)))
                VEC(1) = VEC(1) + VECREC(1)*VAL
                VEC(2) = VEC(2) + VECREC(2)*VAL
                VEC(3) = VEC(3) + VECREC(3)*VAL
              ELSE IF (NDIMPBC.EQ.2) THEN
                VECREC(4) = SQRT(VECREC(4))
                VALM1 = EXP( VECREC(4)*ZZ) &
                      * GMSERFC(VECREC(4)/(TWO*ALPHA)+ALPHA*ZZ)
                VALM2 = EXP(-VECREC(4)*ZZ) &
                      * GMSERFC(VECREC(4)/(TWO*ALPHA)-ALPHA*ZZ)
                VALN1 = EXP( VECREC(4)*ZZ &
                            -(VECREC(4)/(TWO*ALPHA)+ALPHA*ZZ)**2)
                VALN2 = EXP(-VECREC(4)*ZZ &
                            -(VECREC(4)/(TWO*ALPHA)-ALPHA*ZZ)**2)
                SCOS = COS(VECREC(IP1)*(VEC1(IP1)-VEC2(IP1)) &
                          +VECREC(IP2)*(VEC1(IP2)-VEC2(IP2)))
                SSIN = SIN(VECREC(IP1)*(VEC1(IP1)-VEC2(IP1)) &
                          +VECREC(IP2)*(VEC1(IP2)-VEC2(IP2)))
                VAL = SSIN/VECREC(4)*(VALM1+VALM2)
                VAL = VAL*0.25D+00 !! because multiplied by -4 later
                VEC(IP1) = VEC(IP1) + VECREC(IP1)*VAL
                VEC(IP2) = VEC(IP2) + VECREC(IP2)*VAL
!        
                VAL = SCOS &
                    *((VALM1-VALM2)-TWO*ALPHA/(SQRT(PI)*VECREC(4)) &
                    *(VALN1-VALN2))
                VEC(IN1) = VEC(IN1) - VAL*0.25D+00
              ELSE IF (NDIMPBC.EQ.1) THEN
                XX = VECREC(4)/(4.0D+00*ALPHA*ALPHA)
                XX2= -0.25D+00*VECREC(4)
                CALL CALCEI(XX,GAMINC,2)
                GAMID = TWO*EXP(-XX)
                IF (DD.LE.1.0D-08) THEN
                  GSUM = GAMINC
                  GDX = GAMID
                  GDYZ = ZERO
                ELSE
                  !! kappa = 0
                  GSUM = GAMINC
                  GDX  = GAMID
                  GDYZ = ZERO
                  GSUMO= GSUM
                  GDOX = GDX
                  GDOYZ= GDYZ
                  !! kappa = 1
                  SCALX = XX2*DD
                  SCALYZ= XX2
                  XINV  = ONE/XX
                  SCAL2 = EXP(-XX)*XINV
                  GAMINC= SCAL2 - GAMINC
                  GAMID = -(GAMID-TWO*(ONE+XX)*SCAL2)
                  XX2 = XX2*DD
                  DO KK = 1, 100
                    GSUM = GSUM + SCALX*GAMINC/FACT(KK)
                    GDX  = GDX  + SCALX*GAMID/FACT(KK) &
                         - TWO*KK*SCALX*GAMINC/FACT(KK)
                    GDYZ = GDYZ+ TWO*KK*SCALYZ*GAMINC/FACT(KK)
                    IF (ABS(GSUM-GSUMO).LE.1.0D-15 .AND. &
                        ABS(GDX-GDOX).LE.1.0D-15 .AND. &
                        ABS(GDYZ-GDOYZ).LE.1.0D-15) EXIT
                    SCALX = XX2*SCALX
                    SCALYZ= XX2*SCALYZ
                    SCAL2 = SCAL2*XINV
                    GAMINC = -(GAMINC-SCAL2)/(KK+1)
                    GAMID = -(GAMID-TWO*(KK+ONE+XX)*SCAL2)/(KK+1)
                    GSUMO = GSUM
                    GDOX  = GDX
                    GDOYZ = GDYZ
                  END DO
                END IF
                VEC(IP1) = VEC(IP1) &
                  - GSUM*SIN(VECREC(IP1)*(VEC1(IP1)-VEC2(IP1)))*VECREC(IP1)
                VEC(IN1) = VEC(IN1) &
                  + GDYZ*YY*COS(VECREC(IP1)*(VEC1(IP1)-VEC2(IP1)))
                VEC(IN2) = VEC(IN2) &
                  + GDYZ*ZZ*COS(VECREC(IP1)*(VEC1(IP1)-VEC2(IP1)))
              END IF
!        
              !! FOR STRESS TENSOR
              IF (NDIMPBC.EQ.3) THEN
                SCAL = 4.0D+00*PI*TMPCHG*EXP(-SCALEXP*VECREC(4)) &
                     *COS(VECREC(1)*(VEC1(1)-VEC2(1)) &
                         +VECREC(2)*(VEC1(2)-VEC2(2)) &
                         +VECREC(3)*(VEC1(3)-VEC2(3)))/(VOL*VECREC(4))
!               if (srscc) then
!                 write (*,*) "have to modify DFTB_GAMMA_PBC"
!                 write (*,*) "see stress computation in reciprocal"
!                 call abrt
!               end if
                IF (IAT.NE.JAT) SCAL=SCAL+SCAL
                DO K = 1, 3
                  DO L = 1, 3
                    VALS = -2.0D+00*(1.0D+00+SCALEXP*VECREC(4)) &
                           *VECREC(K)*VECREC(L)/VECREC(4)
                    IF (K.EQ.L) VALS = VALS + 1.0D+00
                    STRESS(K,L) = STRESS(K,L) - VALS*SCAL
                  END DO
                END DO
              ELSE IF (NDIMPBC.EQ.2) THEN
                TMP = TMPCHG
                IF (IAT.NE.JAT) TMP = TMP*TWO
                VALG = (VALM1+VALM2)/VECREC(4)
                XI = ZZ*(VALM1-VALM2) &
                   - (VALM1+VALM2)/VECREC(4) &
                   - (VALN1+VALN2)/(ALPHA*SQRT(PI))
                STRESS(IP1,IP1) = STRESS(IP1,IP1) &
                  - PI/VOL*TMP*SCOS &
                  *(VALG+VECREC(IP1)*VECREC(IP1)/(VECREC(4)**2)*XI)
                STRESS(IP1,IP2) = STRESS(IP1,IP2) &
                  - PI/VOL*TMP*SCOS &
                  *(    +VECREC(IP1)*VECREC(IP2)/(VECREC(4)**2)*XI)
                STRESS(IP2,IP1) = STRESS(IP2,IP1) &
                  - PI/VOL*TMP*SCOS &
                  *(    +VECREC(IP2)*VECREC(IP1)/(VECREC(4)**2)*XI)
                STRESS(IP2,IP2) = STRESS(IP2,IP2) &
                  - PI/VOL*TMP*SCOS &
                  *(VALG+VECREC(IP2)*VECREC(IP2)/(VECREC(4)**2)*XI)
                STRESS(IP1,IN1) = STRESS(IP1,IN1) &
                  - PI/VOL*TMP*VECREC(IP1)*ZZ*SSIN*VALG
                STRESS(IP2,IN1) = STRESS(IP2,IN1) &
                  - PI/VOL*TMP*VECREC(IP2)*ZZ*SSIN*VALG
                STRESS(IN1,IP1) = STRESS(IN1,IP1) &
                  - PI/VOL*TMP*VECREC(IP1)*ZZ*SSIN*VALG
                STRESS(IN1,IP2) = STRESS(IN1,IP2) &
                  - PI/VOL*TMP*VECREC(IP2)*ZZ*SSIN*VALG
                STRESS(IN1,IN1) = STRESS(IN1,IN1) &
                  - PI/VOL*TMP*ZZ*SCOS &
                  *(VALM2-VALM1-TWO*ALPHA/(SQRT(PI)*VECREC(4)) &
                  *(VALN2-VALN1))
              ELSE IF (NDIMPBC.EQ.1) THEN
                TMP = TMPCHG
                IF (IAT.NE.JAT) TMP = TMP*TWO
                STRESS(IP1,IP1) = STRESS(IP1,IP1)+ TMP*(-GSUM+GDX)/VOL &
                    *COS(VECREC(IP1)*(VEC1(IP1)-VEC2(IP1)))
                SCAL = -TMP*GSUM/VOL &
                        *SIN(VECREC(IP1)*(VEC1(IP1)-VEC2(IP1)))
                STRESS(IP1,IN1) = STRESS(IP1,IN1)+ SCAL*VECREC(IP1)*YY
                STRESS(IP1,IN2) = STRESS(IP1,IN2)+ SCAL*VECREC(IP1)*ZZ
                STRESS(IN1,IP1) = STRESS(IN1,IP1)+ SCAL*YY*VECREC(IP1)
                STRESS(IN2,IP1) = STRESS(IN2,IP1)+ SCAL*ZZ*VECREC(IP1)
                SCAL = TMP*GDYZ/VOL*COS(VECREC(IP1)*(VEC1(IP1)-VEC2(IP1)))
                STRESS(IN1,IN1) = STRESS(IN1,IN1) + SCAL*YY*YY
                STRESS(IN2,IN2) = STRESS(IN2,IN2) + SCAL*ZZ*ZZ
                IF (DD.GE.1.0D-08) THEN
                  STRESS(IN1,IN2) = STRESS(IN1,IN2) + SCAL*YY*ZZ
                  STRESS(IN2,IN1) = STRESS(IN2,IN1) + SCAL*ZZ*YY
                END IF
              END IF
            ELSE IF (NDER.EQ.2) THEN
              IF (NDIMPBC.EQ.3) THEN
                SCAL = EXP(-SCALEXP*VECREC(4))/VECREC(4)
                VSIN = SIN(VECREC(1)*(VEC1(1)-VEC2(1)) &
                          +VECREC(2)*(VEC1(2)-VEC2(2)) &
                          +VECREC(3)*(VEC1(3)-VEC2(3)))
                VCOS = COS(VECREC(1)*(VEC1(1)-VEC2(1)) &
                          +VECREC(2)*(VEC1(2)-VEC2(2)) &
                          +VECREC(3)*(VEC1(3)-VEC2(3)))
                GRAD = -4.0D+00*PI/VOL*SCAL*VSIN
                HESS = -4.0D+00*PI/VOL*SCAL*VCOS
                DO KI = 1, 3
                  DO KJ = 1, 3
                     VAL= VECREC(KI)*VECREC(KJ) &
                          *(HESS-DISTI*GRAD*DISTI*DISTI)
                    IF (KI.EQ.KJ) VAL = VAL + DISTI*GRAD
                    VAL = VAL*TMPCHG
                    HESSIAN(IGZ+KI,IGZ+KJ) = HESSIAN(IGZ+KI,IGZ+KJ)+VAL
                    HESSIAN(JGZ+KI,JGZ+KJ) = HESSIAN(JGZ+KI,JGZ+KJ)+VAL
                    HESSIAN(IGZ+KI,JGZ+KJ) = HESSIAN(IGZ+KI,JGZ+KJ)-VAL
                    HESSIAN(JGZ+KI,IGZ+KJ) = HESSIAN(JGZ+KI,IGZ+KJ)-VAL
                  END DO
                END DO
!               VEC(1) = VEC(1) + VECREC(1)*VAL
!               VEC(2) = VEC(2) + VECREC(2)*VAL
!               VEC(3) = VEC(3) + VECREC(3)*VAL
              ELSE IF (NDIMPBC.EQ.2) THEN
                VECREC(4) = SQRT(VECREC(4))
                VALM1 = EXP( VECREC(4)*ZZ) &
                      * GMSERFC(VECREC(4)/(TWO*ALPHA)+ALPHA*ZZ)
                VALM2 = EXP(-VECREC(4)*ZZ) &
                      * GMSERFC(VECREC(4)/(TWO*ALPHA)-ALPHA*ZZ)
                VALN1 = EXP( VECREC(4)*ZZ  &
                            -(VECREC(4)/(TWO*ALPHA)+ALPHA*ZZ)**2)
                VALN2 = EXP(-VECREC(4)*ZZ  &
                            -(VECREC(4)/(TWO*ALPHA)-ALPHA*ZZ)**2)
                SCOS = COS(VECREC(IP1)*(VEC1(IP1)-VEC2(IP1))  &
                          +VECREC(IP2)*(VEC1(IP2)-VEC2(IP2)))
                SSIN = SIN(VECREC(IP1)*(VEC1(IP1)-VEC2(IP1))  &
                          +VECREC(IP2)*(VEC1(IP2)-VEC2(IP2)))
                SCAL = PI/VOL*TMPCHG
                DO KI1 = 1, 3
                  IF (KI1.EQ.1) KI = IP1
                  IF (KI1.EQ.2) KI = IP2
                  IF (KI1.EQ.3) KI = IN1
                  DO KJ1 = 1, 3
                    IF (KJ1.EQ.1) KJ = IP1
                    IF (KJ1.EQ.2) KJ = IP2
                    IF (KJ1.EQ.3) KJ = IN1
                    IF (KI1.LE.2 .AND. KJ1.LE.2) THEN
                      VAL= -VECREC(KI)*VECREC(KJ)*SCOS/VECREC(4) &
                         *(VALM1+VALM2)
                    ELSE IF (KI1.LE.2 .OR. KJ1.LE.2) THEN
                      IF (KI1.LE.2) IP0 = KI
                      IF (KJ1.LE.2) IP0 = KJ
                      VAL = -VECREC(IP0)*SSIN/VECREC(4) &
                          *(VECREC(4)*(VALM1-VALM2) &
                           -TWO*ALPHA/SQRT(PI)*(VALN1-VALN2))
                    ELSE IF (KI1.EQ.3 .AND. KJ1.EQ.3) THEN
                      VAL= SCOS/VECREC(4) &
                       *(VECREC(4)*VECREC(4)*(VALM1+VALM2) &
                        -TWO*ALPHA*VECREC(4)/SQRT(PI)*(VALN1+VALN2) &
                       +TWO*TWO*ALPHA*ALPHA*ZZ/SQRT(PI)*(VALN1-VALN2))
                    END IF
                    VAL = VAL*SCAL
                    HESSIAN(IGZ+KI,IGZ+KJ) = HESSIAN(IGZ+KI,IGZ+KJ)+VAL
                    HESSIAN(JGZ+KI,JGZ+KJ) = HESSIAN(JGZ+KI,JGZ+KJ)+VAL
                    HESSIAN(IGZ+KI,JGZ+KJ) = HESSIAN(IGZ+KI,JGZ+KJ)-VAL
                    HESSIAN(JGZ+KI,IGZ+KJ) = HESSIAN(JGZ+KI,IGZ+KJ)-VAL
                  END DO
                END DO
              ELSE IF (NDIMPBC.EQ.1) THEN
                XX = VECREC(4)/(4.0D+00*ALPHA*ALPHA)
                XX2= -0.25D+00*VECREC(4)
                CALL CALCEI(XX,GAMINC,2)
                IF (DD.LE.1.0D-08) THEN
                  GSUM0 = GAMINC
                  GSUM1 = ZERO
                  GSUM2 = ZERO
                ELSE
                  !! kappa = 0
                  GSUM0 = GAMINC
                  GSUM1 = ZERO
                  GSUM2 = ZERO
                  GDSUM0= GSUM0
                  GDSUM1= GSUM1
                  GDSUM2= GSUM2
                  !! kappa = 1
                  SCAL0 = XX2*DD
                  SCAL1 = XX2
                  SCAL2 = XX2/DD
                  XINV  = ONE/XX
                  SCAL  = EXP(-XX)*XINV
                  GAMINC= SCAL  - GAMINC
                  XX2 = XX2*DD
                  DO KK = 1, 100
                    GSUM0 = GSUM0 + SCAL0*GAMINC/FACT(KK)
                    GSUM1 = GSUM1 + TWO*KK*SCAL1*GAMINC/FACT(KK)
                    GSUM2 = GSUM2 + TWO*TWO*KK*(KK-1)*SCAL2*GAMINC/FACT(KK)
                    IF (ABS(GDSUM0-GSUM0).LE.1.0D-15 .AND. &
                        ABS(GDSUM1-GSUM1).LE.1.0D-15 .AND. &
                        ABS(GDSUM2-GSUM2).LE.1.0D-15) EXIT
                    SCAL0 = XX2*SCAL0
                    SCAL1 = XX2*SCAL1
                    SCAL2 = XX2*SCAL2
                    SCAL  = SCAL *XINV
                    GAMINC = -(GAMINC-SCAL)/(KK+1)
                    GDSUM0 = GSUM0
                    GDSUM1 = GSUM1
                    GDSUM2 = GSUM2
                  END DO
                END IF
!        
                SSIN = SIN(VECREC(IP1)*(VEC1(IP1)-VEC2(IP1)))
                SCOS = COS(VECREC(IP1)*(VEC1(IP1)-VEC2(IP1)))
                SCAL = TMPCHG/VOL
                DO KI1 = 1, 3
                  IF (KI1.EQ.1) KI = IP1
                  IF (KI1.EQ.2) KI = IN1
                  IF (KI1.EQ.3) KI = IN2
                  DO KJ1 = 1, 3
                    IF (KJ1.EQ.1) KJ = IP1
                    IF (KJ1.EQ.2) KJ = IN1
                    IF (KJ1.EQ.3) KJ = IN2
                    IF (KI1.EQ.1 .AND. KJ1.EQ.1) THEN
                      VAL = -VECREC(IP1)*VECREC(IP1)*SCOS*GSUM0
                    ELSE IF (KI1.EQ.1 .OR. KJ1.EQ.1) THEN
                      IF (KI1.EQ.1) IN0 = KJ
                      IF (KJ1.EQ.1) IN0 = KI
                      VAL = -VECREC(IP1)*(VEC1(IN0)-VEC2(IN0))*SSIN*GSUM1
                    ELSE IF (KI1.GE.2 .AND. KJ1.GE.2) THEN
                      VAL = SCOS*GSUM2 &
                          *(VEC1(KI1)-VEC2(KI1))*(VEC1(KJ1)-VEC2(KJ1))
                      IF (KI1.EQ.KJ1) THEN
                        VAL = VAL + SCOS*GSUM1
                      END IF
                    END IF
                    VAL = VAL*SCAL
                    HESSIAN(IGZ+KI,IGZ+KJ) = HESSIAN(IGZ+KI,IGZ+KJ)+VAL
                    HESSIAN(JGZ+KI,JGZ+KJ) = HESSIAN(JGZ+KI,JGZ+KJ)+VAL
                    HESSIAN(IGZ+KI,JGZ+KJ) = HESSIAN(IGZ+KI,JGZ+KJ)-VAL
                    HESSIAN(JGZ+KI,IGZ+KJ) = HESSIAN(JGZ+KI,IGZ+KJ)-VAL
                  END DO
                END DO
              END IF
            END IF
            END DO
            END DO
            END DO
!           write (*,*) "itrvec = ", itrvec, val
!        
            !! CHECK THE CONVERGENCE
            IF (NDER.EQ.0) THEN
              IF (4.0D+00*PI*ABS(VAL-OLDVAL)/VOL.LE.THREW) EXIT
              OLDVAL = VAL
            ELSE IF (NDER.EQ.1) THEN
              CALL DCOPY(3,VEC,1,GAMTMP(1,2),1)
              CALL DCOPY(9,STRESS,1,GAMTMP(4,2),1)
              DEVMAX = 0.0D+00
              DO I = 1, 12
                VAL = ABS((GAMTMP(I,1)-GAMTMP(I,2)))
                IF (VAL.GE.DEVMAX) DEVMAX = VAL
              END DO
              DEVMAX = DEVMAX*4.0D+00*PI/VOL
              IF (DEVMAX.LE.THREW) EXIT
              CALL DCOPY(12,GAMTMP(1,2),1,GAMTMP(1,1),1)
            ELSE IF (NDER.EQ.2) THEN
              N = 1
              DO II = 1, 3
                IF (II.EQ.1) THEN
                  IZ = IGZ
                  JZ = IGZ
                ELSE IF (II.EQ.2) THEN
                  IZ = JGZ
                  JZ = JGZ
                ELSE IF (II.EQ.3) THEN
                  IZ = IGZ
                  JZ = JGZ
                END IF
                DO KI = 1, 3
                  DO KJ = 1, 3
                    GAMTMP(N,2) = HESSIAN(IZ+KI,JZ+KJ)
                    N=N+1
                  END DO
                END DO
              END DO
              DEVMAX = 0.0D+00
              DO I = 1, N
                VAL = ABS((GAMTMP(I,1)-GAMTMP(I,2)))
                IF (VAL.GE.DEVMAX) DEVMAX = VAL
              END DO
              IF (DEVMAX.LE.THREW) EXIT
              CALL DCOPY(N,GAMTMP(1,2),1,GAMTMP(1,1),1)
            END IF
          END DO
        END IF !! end of PME brunch
!
        !! Added k-independent terms
        IF (NDER.EQ.1) THEN
          IF (NDIMPBC.EQ.2) THEN
            VEC(IN1) = VEC(IN1) + GMSERF(ALPHA*ZZ)*0.5D+00
          ELSE IF (NDIMPBC.EQ.1) THEN
            IF (DD.GE.1.0D-08) THEN
              VEC(IN1) = VEC(IN1) + TWO*YY*(EXP(-ALPHA*ALPHA*DD)-ONE)/DD
              VEC(IN2) = VEC(IN2) + TWO*ZZ*(EXP(-ALPHA*ALPHA*DD)-ONE)/DD
            END IF
            !! because multiplied by -4pi later
            VEC(IP1) = VEC(IP1) / (-4.0D+00*PI)
            VEC(IN1) = VEC(IN1) / (-4.0D+00*PI)
            VEC(IN2) = VEC(IN2) / (-4.0D+00*PI)
          END IF
        ELSE IF (NDER.EQ.2) THEN
          IF (NDIMPBC.EQ.2) THEN
            VAL =-TWO*TWO*ALPHA*SQRT(PI)/VOL*EXP(-ALPHA*ALPHA*ZZ*ZZ)
            VAL = VAL*TMPCHG
            HESSIAN(IGZ+IN1,IGZ+IN1) = HESSIAN(IGZ+IN1,IGZ+IN1)+VAL
            HESSIAN(JGZ+IN1,JGZ+IN1) = HESSIAN(JGZ+IN1,JGZ+IN1)+VAL
            HESSIAN(IGZ+IN1,JGZ+IN1) = HESSIAN(IGZ+IN1,JGZ+IN1)-VAL
            HESSIAN(JGZ+IN1,IGZ+IN1) = HESSIAN(JGZ+IN1,IGZ+IN1)-VAL
          ELSE IF (NDIMPBC.EQ.1 .AND. DD.GE.1.0D-08) THEN
            SCAL = TWO*TMPCHG/(VOL*DD*DD)
            DO KI1 = 1, 2
              IF (KI1.EQ.1) KI = IN1
              IF (KI1.EQ.1) V1 = YY
              IF (KI1.EQ.2) KI = IN2
              IF (KI1.EQ.2) V1 = ZZ
              DO KJ1 = 1, 2
                IF (KJ1.EQ.1) KJ = IN1
                IF (KJ1.EQ.1) V2 = YY
                IF (KJ1.EQ.2) KJ = IN2
                IF (KJ1.EQ.2) V2 = ZZ
                VAL = -TWO*ALPHA*ALPHA*V1*V2*EXP(-ALPHA*ALPHA*DD)*DD &
                      -TWO*V1*V2*(EXP(-ALPHA*ALPHA*DD)-ONE)
                IF (KI1.EQ.KJ1) VAL = VAL + (EXP(-ALPHA*ALPHA*DD)-ONE)*DD
                VAL = VAL*SCAL
                HESSIAN(IGZ+KI,IGZ+KJ) = HESSIAN(IGZ+KI,IGZ+KJ)+VAL
                HESSIAN(JGZ+KI,JGZ+KJ) = HESSIAN(JGZ+KI,JGZ+KJ)+VAL
                HESSIAN(IGZ+KI,JGZ+KJ) = HESSIAN(IGZ+KI,JGZ+KJ)-VAL
                HESSIAN(JGZ+KI,IGZ+KJ) = HESSIAN(JGZ+KI,IGZ+KJ)-VAL
              END DO
            END DO
          END IF
        END IF
!
        DO II = 1, ISHMAX
          ISH = ISH0 + II
          DO JJ = 1, JSHMAX
            JSH = JSH0 + JJ
            IF (NDER.EQ.0) THEN
              CALL DFTB_CNVSQ(ISH,JSH,NSEQ)
              IF (NDIMPBC.EQ.3) THEN
                !! RECIPROCAL
                GAMMA2(NSEQ) = GAMMA2(NSEQ) + 4.0D+00*PI*VAL/VOL
                !! CORRECTION FOR CHARGED SYSTEMS
                !! (sum_A q_A)^2 = sum_AB q_A q_B
                !! should it be added for cell deformation?
                GAMMA2(NSEQ) = GAMMA2(NSEQ) - PI/(VOL*ALPHA*ALPHA)
              ELSE IF (NDIMPBC.EQ.2) THEN
                !! RECIPROCAL
                GAMMA2(NSEQ) = GAMMA2(NSEQ) + PI*VAL/VOL
                !! CORRECTION FOR CHARGED SYSTEMS (K=0)
                GAMMA2(NSEQ) = GAMMA2(NSEQ) &
                  - (SQRT(PI)*EXP(-ALPHA*ALPHA*ZZ*ZZ)/ALPHA &
                    +PI*ZZ*GMSERF(ALPHA*ZZ))/VOL*2.0D+00
              ELSE IF (NDIMPBC.EQ.1) THEN
                GAMMA2(NSEQ) = GAMMA2(NSEQ) + VAL/VOL
                IF (DD.GE.1.0D-08) THEN
                  BS2 = DD*ALPHA*ALPHA
                  CALL CALCEI(BS2,EI,2)
                  GAMMA2(NSEQ) = GAMMA2(NSEQ) &
                    - (0.5772156649D+00+EI+LOG(BS2))/VOL
                END IF
              END IF
              !! SELF-ENERGY CONTRIBUTIONS
              IF (IAT.EQ.JAT) &
                GAMMA2(NSEQ) = GAMMA2(NSEQ) - 2.0D+00*ALPHA/SQRT(PI)
            ELSE IF (NDER.EQ.1) THEN
              TMP = DQ(ISH)*DQ(JSH)
              IF (TDDFTB) THEN
                TMP = TMP &
                    + DQES(ISH)*DQ(JSH) + DQ(ISH)*DQES(JSH) &
                    + 2.0D+00*DQXY(ISH)*DQXY(JSH)
              END IF
              VAL = -4.0D+00*PI*TMP/VOL
              IF (NI.EQ.0 .AND. KI0.EQ.0) THEN
                DO K = 1, 3
                  EGRAD(K,IAT) = EGRAD(K,IAT) + VEC(K)*VAL
                  EGRAD(K,JAT) = EGRAD(K,JAT) - VEC(K)*VAL
                END DO
              ELSE
                VAL = -4.0D+00*PI/VOL
                SHIFT(JSH) = SHIFT(JSH) - VEC(KI0)*VAL*DQ(ISH)
                SHIFT(ISH) = SHIFT(ISH) - VEC(KI0)*VAL*DQ(JSH)
              END IF
              !! STRESS CONTRIBUTIONS FROM K=0 IN 2D-EWALD
              IF (NDIMPBC.EQ.2) THEN
                STRESS(IP1,IP1) = STRESS(IP1,IP1) &
                  + TWO*PI/VOL*TMP &
                    *(EXP(-ALPHA*ALPHA*ZZ*ZZ)/(ALPHA*SQRT(PI)) &
                     +ZZ*GMSERF(ALPHA*ZZ))*SCAL
                STRESS(IP2,IP2) = STRESS(IP2,IP2) &
                  + TWO*PI/VOL*TMP &
                    *(EXP(-ALPHA*ALPHA*ZZ*ZZ)/(ALPHA*SQRT(PI)) &
                     +ZZ*GMSERF(ALPHA*ZZ))*SCAL
                STRESS(IN1,IN1) = STRESS(IN1,IN1) &
                  - TWO*PI/VOL*TMP*ZZ*GMSERF(ALPHA*ZZ)*SCAL
              ELSE IF (NDIMPBC.EQ.1 .AND. DD.GE.1.0D-08) THEN
                SCAL = 2.0D+00
                IF (IAT.EQ.JAT) SCAL = 1.0D+00
                BS2 = DD*ALPHA*ALPHA
                CALL CALCEI(BS2,EI,2)
                STRESS(IP1,IP1) = STRESS(IP1,IP1) &
                  + TMP*(0.5772156649D+00+EI+LOG(BS2))/VOL*SCAL
                SCAL = SCAL*TMP*TWO*(EXP(-ALPHA*ALPHA*DD)-ONE) &
                       /(VOL*DD)
                STRESS(IN1,IN1) = STRESS(IN1,IN1) + SCAL*YY*YY
                STRESS(IN1,IN2) = STRESS(IN1,IN2) + SCAL*YY*ZZ
                STRESS(IN2,IN1) = STRESS(IN2,IN1) + SCAL*ZZ*YY
                STRESS(IN2,IN2) = STRESS(IN2,IN2) + SCAL*ZZ*ZZ
              END IF
            END IF
!           IF (IAT.EQ.JAT.AND.II.EQ.JJ.AND.NDER.le.1) EXIT
            IF (IAT.EQ.JAT.AND.II.EQ.JJ) EXIT
          END DO
        END DO
      END IF
      JSH0 = JSH0 + JSHMAX !! MAXANG(ISPE(J))
    END DO
    ISH0 = ISH0 + ISHMAX !! MAXANG(ISPE(I))
  END DO
! call prtril(gamma2,nat)
  RETURN
!
END SUBROUTINE DFTB_GAMMA_PBC
!
!-----------------------------------------------------------------------
!*MODULE DFTBPB    *DECK DFTB_PBCSTR
!>
!>     @brief Calculate stress tensor.
!>
!>     @details Calculate stress tensor in PBC.
!>
!>     @author Yoshio Nishimoto
!>
SUBROUTINE DFTB_PBCSTR(TMP,NDIMPBC_,SOME)
!
  IMPLICIT NONE
!
  DOUBLE PRECISION TMP(3,3),PRESS
  INTEGER I,J,NDIMPBC_
  LOGICAL SOME
!
  CALL DSCAL(9,-0.5D+00/VOLPBC,STRESS,1)
  PRESS =  (STRESS(1,1)+STRESS(2,2)+STRESS(3,3))/3.0D+00
  IF (MASWRK.AND.SOME) THEN
    WRITE (IW,'(X,"STRESS TENSOR IN AU")')
    DO I = 1, 3
      WRITE (IW,'(3F20.10)') (STRESS(I,J),J=1,3)
    END DO
    WRITE (IW,'(X,"PRESSURE = ", F20.10," AU")') PRESS
  END IF
!
  CALL DGEMM('N','T',3,3,3, &
             -VOLPBC/(2.0D+00*PI),STRESS,3,PBCBOX(1,4),3, &
             0.0D+00,TMP,3)
  CALL TRPOSQ(TMP,3)
  IF (MASWRK.AND.SOME) THEN
    WRITE (IW,'(X,"DERIVATIVE WITH RESPECT TO CELL DEFORMATION")')
    DO I = 1, 3
      WRITE (IW,'(3F20.10)') (TMP(I,J),J=1,3)
    END DO
  END IF
!
! CALL DFTB_PRJLAT(LATOPT,TMP,PBCBOX,PRESS)
! IF (MASWRK) THEN
!     WRITE (IW,'(X,"PROJECTED DERIVATIVES")')
!     DO I = 1, 3
!       WRITE (IW,'(3F20.10)') (TMP(I,J),J=1,3)
!     END DO
! END IF
  NDIMPBC_ = NDIMPBC
!
  RETURN
!
END SUBROUTINE DFTB_PBCSTR
!
!-----------------------------------------------------------------------
!*MODULE DFTBPB    *DECK DFTB_GAMFMO_PBC
!>
!>     @brief Calculate gamma for FMO.
!>
!>     @details Calculate gamma for FMO/PBC.
!>
!>     @author Yoshio Nishimoto
!>
SUBROUTINE DFTB_GAMFMO_PBC(IAT,JAT,C1,C2,GAM2,GAM3I,GAM3J,TMP2,TMP3I, &
                           TMP3J,DAMPXH,ISP,JSP)
USE MX_LIMITS, ONLY: MXATM
!
  IMPLICIT NONE
!
  INTEGER MXSPE
  PARAMETER (MXSPE=10)
!
  INTEGER ISH0,JSH0,ISHMAX,JSHMAX,MAXANG,ISP,JSP,ITRVEC,IXVEC,IYVEC
  INTEGER IZVEC,II,NSPE,IND,IDFTB,NEEDSK,LSKHTAB,LSKSTAB,LSKGRID,LSKSELF
  INTEGER JJ,IAT,JAT,N,I,MAXVEC,ISPE,IDFTBD
  LOGICAL DFTBFL,SCC,SRSCC,DFTB3,DAMPXH0,LCDFTB
  DOUBLE PRECISION SPEI,SPEJ,GAM2,GAM3I,GAM3J,DIST,DISTI,DITERM
  DOUBLE PRECISION GMSERF,HUBA,HUBDERI,ETEMP,DFTBDP,DAMPXHE,HUBDER,ZREF
  DOUBLE PRECISION SPNCST,SPE,SKDIM,SKSPIN,QREFL,HUBBL,QREF,HUBBS,SKCUT2
  DOUBLE PRECISION HUBB,HUBDERJ,GAM2T,GAM3IT,GAM3JT,DFTB_GAM2,DEVMAX
  DOUBLE PRECISION VAL,OLDVAL,GAMK,PARAMDIR
!
  COMMON /DFTB  / DFTBFL,SCC,SRSCC,DFTB3,DAMPXH0,LCDFTB
  COMMON /DFTBPR/ ETEMP,DFTBDP(MXSPE*14),DAMPXHE,HUBDER(MXSPE), &
                  ZREF(MXATM),SPNCST(6,MXSPE),SPE(MXATM),NSPE,  &
                  MAXANG(MXATM),ISPE(MXATM),IND(MXATM+1),IDFTBD, &
                  PARAMDIR
  COMMON /DFTBSK/ SKDIM(MXSPE,MXSPE),SKSPIN(MXSPE),QREFL(3,MXSPE),  &
                  HUBBL(3,MXSPE),QREF(MXSPE),HUBBS(MXSPE),SKCUT2,   &
                  NEEDSK,LSKHTAB(MXSPE,MXSPE),LSKSTAB(MXSPE,MXSPE), &
                  LSKGRID(MXSPE,MXSPE),LSKSELF(MXSPE)
!
  DOUBLE PRECISION C1(3),C2(3),TMP2(*),TMP3I(*),TMP3J(*)
  DOUBLE PRECISION VEC(3),VECREC(4),GAMTMP(27,2)
!
  LOGICAL DAMPXH
!
! DOUBLE PRECISION H/8HH       /
!
  ISH0 = 0
  JSH0 = 0
  IF (SRSCC) THEN
    ISHMAX = MAXANG(ISP)
    JSHMAX = MAXANG(JSP)
    !! ish0,jsh0,ishell has to be defined
    call abrt
  ELSE
    ISHMAX = 1
    JSHMAX = 1
  END IF
  SPEI = SPE(ISP)
  SPEJ = SPE(JSP)
!
! ---- SHORT-RANGE + REAL-SPACE CONTRIBUTIONS -----
!
  GAM2 = 0.0D+00
  GAM3I= 0.0D+00
  GAM3J= 0.0D+00
  CALL VCLR(GAMTMP,1,54)
  DO ITRVEC = 0, 100
    DO IXVEC = -ITRVEC, ITRVEC
    DO IYVEC = -ITRVEC, ITRVEC
    DO IZVEC = -ITRVEC, ITRVEC
    IF (ABS(IXVEC).NE.ITRVEC .AND. ABS(IYVEC).NE.ITRVEC .AND. &
        ABS(IZVEC).NE.ITRVEC) CYCLE
!
    IF (PERIOD) THEN
      VEC(1) = C1(1) - C2(1) &
             - BOXREAL(1,1)*IXVEC  &
             - BOXREAL(2,1)*IYVEC  &
             - BOXREAL(3,1)*IZVEC
      VEC(2) = C1(2) - C2(2) &
             - BOXREAL(1,2)*IXVEC  &
             - BOXREAL(2,2)*IYVEC  &
             - BOXREAL(3,2)*IZVEC
      VEC(3) = C1(3) - C2(3) &
             - BOXREAL(1,3)*IXVEC  &
             - BOXREAL(2,3)*IYVEC  &
             - BOXREAL(3,3)*IZVEC
    ELSE
      VEC(1) = C1(1) - C2(1)
      VEC(2) = C1(2) - C2(2)
      VEC(3) = C1(3) - C2(3)
    END IF
    DIST = VEC(1)*VEC(1)+VEC(2)*VEC(2)+VEC(3)*VEC(3)
    IF (DIST.GE.1.0D-15) THEN
      DIST = SQRT(DIST)
      DISTI= 1.0D+00/DIST
    END IF
!
    DITERM = 0.0D+00
    IF (PERIOD .AND. DIST.NE.0.0D+00) THEN
      !! 1/R IS ALREADY ADDED IN DFTB_GAMMA2F
      DITERM = -GMSERF(EWALPHA*DIST)*DISTI
    END IF
!
    DO II = 1, ISHMAX
!     ISH = ISH0 + II
      HUBA = HUBBL(II,ISP)
!     NSP1 = ISHELL(ISH)
      HUBDERI = HUBDER(ISP)
      DO JJ = 1, JSHMAX
!       JSH = JSH0 + JJ
        HUBB = HUBBL(JJ,JSP)
!       NSP2 = ISHELL(JSH)
        HUBDERJ = HUBDER(JSP)
        IF (DFTB3) THEN
          CALL DFTB_GAM23(1,DIST,GAM2T,GAM3IT,GAM3JT,TMP2,TMP3I,TMP3J, &
               DAMPXH,DAMPXHE,SPEI,SPEJ,HUBDERI,HUBDERJ,DFTB3)
        ELSE
          GAM2T = DFTB_GAM2(DIST,DISTI,TMP2,DAMPXH,DAMPXHE,SPEI,SPEJ)
        END IF
!       CALL DFTB_CNVSQ(ISH,JSH,NSEQ)
!       GAM2(II,JJ) = GAM2(II,JJ) + GAM2T + DITERM
!       IF (DFTB3) THEN
!         GAM3I(II,JJ) = GAM3I(II,JJ) + GAM3IT
!         GAM3J(II,JJ) = GAM3J(II,JJ) + GAM3JT
!       END IF
        GAM2 = GAM2 + GAM2T + DITERM
        IF (DFTB3) THEN
          GAM3I = GAM3I + GAM3IT
          GAM3J = GAM3J + GAM3JT
        END IF
        IF (IAT.EQ.JAT.AND.II.EQ.JJ) EXIT
      END DO
    END DO
!
    END DO
    END DO
    END DO
!   write (*,*) "itrvec = ", itrvec, gamma2(nseq)
!
    !! CHECK THE CONVERGENCE
    N = 0
!   DO II = 1, ISHMAX
!     DO JJ = 1, JSHMAX
!       N = N + 1
!       GAMTMP(N,2) = GAM2(II,JJ)
!       IF (DFTB3) THEN
!         N = N + 1
!         GAMTMP(N,2) = GAM3I(II,JJ)
!         N = N + 1
!         GAMTMP(N,2) = GAM3J(II,JJ)
!       END IF
!       IF (IAT.EQ.JAT.AND.II.EQ.JJ) EXIT
!     END DO
!   END DO
    N = N + 1
    GAMTMP(N,2) = GAM2
    IF (DFTB3) THEN
      N = N + 1
      GAMTMP(N,2) = GAM3I
      N = N + 1
      GAMTMP(N,2) = GAM3J
    END IF
    DEVMAX = 0.0D+00
    DO I = 1, N
      VAL = ABS((GAMTMP(I,1)-GAMTMP(I,2)))
      IF (VAL.GE.DEVMAX) DEVMAX = VAL
    END DO
    IF (DEVMAX.LE.THREW.AND.ITRVEC.GE.NSPILL) EXIT
!   IF (DEVMAX.LE.THREW.AND.ITRVEC.GT.0) EXIT
    CALL DCOPY(N,GAMTMP(1,2),1,GAMTMP(1,1),1)
  END DO
!     write (*,*) "converged short-range at", itrvec
!
!     ---- LONG-RANGE (RECIPROCAL) AND SELF-ENERGY CONTRIBUTIONS ---
!
  IF (PERIOD) THEN
    VAL = 0.0D+00
    OLDVAL = 0.0D+00
    VEC(1) = C1(1) - C2(1)
    VEC(2) = C1(2) - C2(2)
    VEC(3) = C1(3) - C2(3)
    MAXVEC = 100
    IF (PME) MAXVEC = 0
    DO ITRVEC = 1, MAXVEC
      DO IXVEC = -ITRVEC, ITRVEC
      DO IYVEC = -ITRVEC, ITRVEC
      DO IZVEC = -ITRVEC, ITRVEC
      IF (ABS(IXVEC).NE.ITRVEC .AND. ABS(IYVEC).NE.ITRVEC .AND. &
          ABS(IZVEC).NE.ITRVEC) CYCLE
!
      VECREC(1) = BOXREC(1,1)*IXVEC &
                + BOXREC(2,1)*IYVEC &
                + BOXREC(3,1)*IZVEC
      VECREC(2) = BOXREC(1,2)*IXVEC &
                + BOXREC(2,2)*IYVEC &
                + BOXREC(3,2)*IZVEC
      VECREC(3) = BOXREC(1,3)*IXVEC &
                + BOXREC(2,3)*IYVEC &
                + BOXREC(3,3)*IZVEC
      VECREC(4) = VECREC(1)**2 + VECREC(2)**2 + VECREC(3)**2
      !! i*sin disappears
      VAL = VAL + EXP(-SCALEXP*VECREC(4))/VECREC(4) &
          * COS(VECREC(1)*VEC(1) &
               +VECREC(2)*VEC(2) &
               +VECREC(3)*VEC(3))
      END DO
      END DO
      END DO
!
      !! CHECK THE CONVERGENCE
      IF (4.0D+00*PI*ABS(VAL-OLDVAL)/VOLPBC.LE.THREW) EXIT
      OLDVAL = VAL
    END DO
!     write (*,*) "converged  long-range at", itrvec
!
    GAMK = 4.0D+00*PI*VAL/VOLPBC   &
         - PI/(VOLPBC*EWALPHA*EWALPHA)
    IF (IAT.EQ.JAT) GAMK = GAMK - 2.0D+00*EWALPHA/SQRT(PI)
!   DO II = 1, ISHMAX
!     DO JJ = 1, JSHMAX
!       GAM2(II,JJ) = GAM2(II,JJ) + GAMK
!       IF (IAT.EQ.JAT.AND.II.EQ.JJ) EXIT
!     END DO
!   END DO
    GAM2 = GAM2 + GAMK
  END IF
!          write (*,'("refgam:",i1,i1,f20.10)') iat,jat,gam2
!
  RETURN
!
END SUBROUTINE DFTB_GAMFMO_PBC
!
!-----------------------------------------------------------------------
!*MODULE DFTBPB    *DECK DFTB_GAMDER_PBC
!>
!>     @brief Calculate gamma derivative.
!>
!>     @details Calculate gamma derivative for PBC.
!>
!>     @author Yoshio Nishimoto
!>
SUBROUTINE DFTB_GAMDER_PBC(IAT,JAT,C1,C2,GAM2,GAM3I,GAM3J,  &
                           GAM2ST,GAM3IST,GAM3JST,          &
                           TMP2,TMP3I,TMP3J,DAMPXH,ISP,JSP)
USE MX_LIMITS, ONLY: MXATM
!
  IMPLICIT NONE
!
  INTEGER MXSPE
  PARAMETER (MXSPE=10)
!
  LOGICAL DFTBFL,SCC,SRSCC,DFTB3,DAMPXH0,LCDFTB
!
  COMMON /DFTB  / DFTBFL,SCC,SRSCC,DFTB3,DAMPXH0,LCDFTB
  COMMON /DFTBPR/ ETEMP,DFTBDP(MXSPE*14),DAMPXHE,HUBDER(MXSPE), &
                  ZREF(MXATM),SPNCST(6,MXSPE),SPE(MXATM),NSPE,  &
                  MAXANG(MXATM),ISPE(MXATM),IND(MXATM+1),IDFTBD, &
                  PARAMDIR
  COMMON /DFTBSK/ SKDIM(MXSPE,MXSPE),SKSPIN(MXSPE),QREFL(3,MXSPE),  &
                  HUBBL(3,MXSPE),QREF(MXSPE),HUBBS(MXSPE),SKCUT2,   &
                  NEEDSK,LSKHTAB(MXSPE,MXSPE),LSKSTAB(MXSPE,MXSPE), &
                  LSKGRID(MXSPE,MXSPE),LSKSELF(MXSPE)
!
  INTEGER ISH0,JSH0,ISHMAX,JSHMAX,MAXANG,ISP,JSP,ITRVEC,IXVEC,IYVEC
  INTEGER IZVEC,II,NSPE,ISPE,IND,IDFTBD,NEEDSK,LSKHTAB,LSKSTAB,LSKGRID
  INTEGER LSKSELF,JJ,K,L,IAT,JAT,N,I,NBODY
  DOUBLE PRECISION C1(3),C2(3),TMP2(10),TMP3I(11),TMP3J(11)
  DOUBLE PRECISION VEC(3),VECREC(4),GAMTMP(36,2)
  DOUBLE PRECISION GAM2(3),GAM3I(3),GAM3J(3),GAM2ST(3,3),GAM3IST(3,3)
  DOUBLE PRECISION GAM3JST(3,3),SPEI,SPEJ,DIST,DISTI,DITERM,GMSERF
  DOUBLE PRECISION HUBA,HUBDERI,ETEMP,DFTBDP,DAMPXHE,HUBDER,ZREF,SPNCST,SPE
  DOUBLE PRECISION SKDIM,SKSPIN,QREFL,HUBBL,QREF,HUBBS,SKCUT2,HUBB,HUBDERJ
  DOUBLE PRECISION GAM2T,GAM3IT,GAM3JT,SCAL,DEVMAX,VAL,OLDVAL,EXPREC
  DOUBLE PRECISION RECINV,VECTMP,SINREC,VALS,PARAMDIR
!
  LOGICAL DAMPXH
!
  DOUBLE PRECISION H/8HH       /
!
  CALL DCOPY(3,0.0D+00,0,GAM2,1)
  CALL DCOPY(9,0.0D+00,0,GAM2ST,1)
  IF (DFTB3) THEN
    CALL DCOPY(3,0.0D+00,0,GAM3I,1)
    CALL DCOPY(3,0.0D+00,0,GAM3J,1)
    CALL DCOPY(9,0.0D+00,0,GAM3IST,1)
    CALL DCOPY(9,0.0D+00,0,GAM3JST,1)
  END IF
!
  ISH0 = 0
  JSH0 = 0
  IF (SRSCC) THEN
    ISHMAX = MAXANG(ISP)
    JSHMAX = MAXANG(JSP)
    !! ish0,jsh0,ishell has to be defined
    call abrt
  ELSE
    ISHMAX = 1
    JSHMAX = 1
  END IF
  SPEI = SPE(ISP)
  SPEJ = SPE(JSP)
!
! ---- SHORT-RANGE + REAL-SPACE CONTRIBUTIONS -----
!
  CALL VCLR(GAMTMP,1,72)
  DO ITRVEC = 0, 100
    DO IXVEC = -ITRVEC, ITRVEC
    DO IYVEC = -ITRVEC, ITRVEC
    DO IZVEC = -ITRVEC, ITRVEC
    IF (ABS(IXVEC).NE.ITRVEC .AND. ABS(IYVEC).NE.ITRVEC .AND. &
        ABS(IZVEC).NE.ITRVEC) CYCLE
!
    IF (PERIOD) THEN
      VEC(1) = C1(1) - C2(1) &
             - BOXREAL(1,1)*IXVEC  &
             - BOXREAL(2,1)*IYVEC  &
             - BOXREAL(3,1)*IZVEC
      VEC(2) = C1(2) - C2(2) &
             - BOXREAL(1,2)*IXVEC  &
             - BOXREAL(2,2)*IYVEC  &
             - BOXREAL(3,2)*IZVEC
      VEC(3) = C1(3) - C2(3) &
             - BOXREAL(1,3)*IXVEC  &
             - BOXREAL(2,3)*IYVEC  &
             - BOXREAL(3,3)*IZVEC
    ELSE
      VEC(1) = C1(1) - C2(1)
      VEC(2) = C1(2) - C2(2)
      VEC(3) = C1(3) - C2(3)
    END IF
    DIST = VEC(1)*VEC(1)+VEC(2)*VEC(2)+VEC(3)*VEC(3)
    IF (DIST.LT.1.0D-15) CYCLE
    DIST = SQRT(DIST)
    DISTI= 1.0D+00/DIST
!
    DITERM = 0.0D+00
    IF (PERIOD) THEN
      DITERM = DITERM - (2.0D+00*EWALPHA/PISQ &
                    *EXP(-EWALPHA*EWALPHA*DIST*DIST) &
           - GMSERF(EWALPHA*DIST)*DISTI)*DISTI
    END IF
!
    DO II = 1, ISHMAX
!     ISH = ISH0 + II
      HUBA = HUBBL(II,ISP)
!     NSP1 = ISHELL(ISH)
      HUBDERI = HUBDER(ISP)
      DO JJ = 1, JSHMAX
!       JSH = JSH0 + JJ
        HUBB = HUBBL(JJ,JSP)
!       NSP2 = ISHELL(JSH)
        HUBDERJ = HUBDER(JSP)
        CALL DFTB_GAM23D(DIST,DISTI,GAM2T,GAM3IT,GAM3JT, &
                         TMP2,TMP3I,TMP3J, &
                         DAMPXH,DAMPXHE,SPEI,SPEJ,HUBDERI,HUBDERJ,DFTB3)
!   write (*,*) gam2t,diterm
!       GAM2 = GAM2 + GAM2T + DITERM
!       VAL = CT1I*DQJ*(GAM2T+DITERM)
!       IF (DFTB3) THEN
!         GAM3I = GAM3I + GAM3IT
!         GAM3J = GAM3J + GAM3JT
!         VAL = VAL + (CT2I*GAM3IT+CT1I*DQJ*GAM3JT)
!       END IF
!
        SCAL = 2.0D+00
    !   IF (IAT.EQ.JAT) SCAL = 1.0D+00
!       scal = 1.0d+00
        DO K = 1, 3
          GAM2(K) = GAM2(K) + (GAM2T+DITERM)*VEC(K)*DISTI
          IF (DFTB3) THEN
            GAM3I(K) = GAM3I(K) + GAM3IT*VEC(K)*DISTI
            GAM3J(K) = GAM3J(K) + GAM3JT*VEC(K)*DISTI
          END IF
          DO L = 1, 3
            GAM2ST(K,L) = GAM2ST(K,L) + (GAM2T+DITERM)*VEC(K)*VEC(L)*DISTI*SCAL
            IF (DFTB3) THEN
              GAM3IST(K,L) = GAM3IST(K,L) + GAM3IT*VEC(K)*VEC(L)*DISTI*SCAL
              GAM3JST(K,L) = GAM3JST(K,L) + GAM3JT*VEC(K)*VEC(L)*DISTI*SCAL
            END IF
          END DO
        END DO
        IF (IAT.EQ.JAT.AND.II.EQ.JJ) EXIT
      END DO
    END DO
!
    END DO
    END DO
    END DO
!   write (*,*) "itrvec = ", itrvec, gamma2(nseq)
!
    !! CHECK THE CONVERGENCE
    N = 12
    CALL DCOPY(3,GAM2,1,GAMTMP(1,2),1)
    CALL DCOPY(9,GAM2ST,1,GAMTMP(4,2),1)
    IF (DFTB3) THEN
      N = 36
      CALL DCOPY(3,GAM3I,1,GAMTMP(13,2),1)
      CALL DCOPY(3,GAM3J,1,GAMTMP(16,2),1)
      CALL DCOPY(9,GAM3IST,1,GAMTMP(19,2),1)
      CALL DCOPY(9,GAM3JST,1,GAMTMP(28,2),1)
    END IF
    DEVMAX = 0.0D+00
    DO I = 1, N
      VAL = ABS((GAMTMP(I,1)-GAMTMP(I,2)))
      IF (VAL.GE.DEVMAX) DEVMAX = VAL
    END DO
    IF (DEVMAX.LE.THREW.AND.ITRVEC.GE.NSPILL) EXIT
!   IF (DEVMAX.LE.THREW.AND.ITRVEC.GT.0) EXIT
    CALL DCOPY(N,GAMTMP(1,2),1,GAMTMP(1,1),1)
  END DO
!     write (*,*) "converged short-range at", itrvec
!
!     ---- LONG-RANGE (RECIPROCAL) AND SELF-ENERGY CONTRIBUTIONS ---
!
  IF (PERIOD) THEN
    VAL = 0.0D+00
    OLDVAL = 0.0D+00
    CALL VCLR(VEC,1,3)
    CALL VCLR(GAMTMP,1,72)
!   CALL VCLR(GAMTMP,1,3)
!   CALL DCOPY(9,GAM2ST,1,GAMTMP(4,1),1)
    DO ITRVEC = 1, 100
      DO IXVEC = -ITRVEC, ITRVEC
      DO IYVEC = -ITRVEC, ITRVEC
      DO IZVEC = -ITRVEC, ITRVEC
      IF (ABS(IXVEC).NE.ITRVEC .AND. ABS(IYVEC).NE.ITRVEC .AND. &
          ABS(IZVEC).NE.ITRVEC) CYCLE
!
      VECREC(1) = BOXREC(1,1)*IXVEC &
                + BOXREC(2,1)*IYVEC &
                + BOXREC(3,1)*IZVEC
      VECREC(2) = BOXREC(1,2)*IXVEC &
                + BOXREC(2,2)*IYVEC &
                + BOXREC(3,2)*IZVEC
      VECREC(3) = BOXREC(1,3)*IXVEC &
                + BOXREC(2,3)*IYVEC &
                + BOXREC(3,3)*IZVEC
      VECREC(4) = VECREC(1)**2 + VECREC(2)**2 + VECREC(3)**2
      EXPREC    = EXP(-SCALEXP*VECREC(4))
      RECINV    = 1.0D+00/VECREC(4)
      VECTMP    = VECREC(1)*(C1(1)-C2(1))   &
                 +VECREC(2)*(C1(2)-C2(2))   &
                 +VECREC(3)*(C1(3)-C2(3))
      SINREC    = SIN(VECTMP)
!!    SINREC    = SIN(VECREC(1)*(C1(1)-C2(1))   &
!!                   +VECREC(2)*(C1(2)-C2(2))   &
!!                   +VECREC(3)*(C1(3)-C2(3)))
!     COSREC    = SQRT(1.0D+00-SINREC*SINREC)
!
!     VAL = EXP(-SCALEXP*VECREC(4))/VECREC(4) &
!         * SIN(VECREC(1)*(C1(1)-C2(1))   &
!              +VECREC(2)*(C1(2)-C2(2))   &
!              +VECREC(3)*(C1(3)-C2(3)))
      VAL = EXPREC*RECINV*SINREC
      VEC(1) = VEC(1) + VECREC(1)*VAL
      VEC(2) = VEC(2) + VECREC(2)*VAL
      VEC(3) = VEC(3) + VECREC(3)*VAL
      !! For stress tensor
!     SCAL = 4.0D+00*PI*EXP(-SCALEXP*VECREC(4)) &
!          *COS(VECREC(1)*(C1(1)-C2(1)) &
!              +VECREC(2)*(C1(2)-C2(2)) &
!              +VECREC(3)*(C1(3)-C2(3)))/(VOLPBC*VECREC(4))
!     SCAL = 4.0D+00*PI*EXPREC*COSREC*VOLINV*RECINV
      SCAL = 4.0D+00*PI*EXPREC*VOLINV*RECINV &
           *COS(VECTMP)
!!         *COS(VECREC(1)*(C1(1)-C2(1)) &
!!             +VECREC(2)*(C1(2)-C2(2)) &
!!             +VECREC(3)*(C1(3)-C2(3)))
      IF (IAT.NE.JAT) SCAL=SCAL+SCAL
      !!asdf
      if (iat.eq.jat) scal=scal+scal
!     VAL = SCAL*SCFMO
      DO K = 1, 3
        DO L = 1, 3
!         VALS = -2.0D+00*(1.0D+00+SCALEXP*VECREC(4)) &
!                *VECREC(K)*VECREC(L)/VECREC(4)
          VALS = -2.0D+00*(1.0D+00+SCALEXP*VECREC(4)) &
                 *VECREC(K)*VECREC(L)*RECINV
          IF (K.EQ.L) VALS = VALS + 1.0D+00
          GAM2ST(K,L) = GAM2ST(K,L) - VALS*SCAL
        END DO
      END DO
      END DO
      END DO
      END DO
!
      !! CHECK THE CONVERGENCE
      CALL DCOPY(3,VEC,1,GAMTMP(1,2),1)
      CALL DCOPY(9,GAM2ST,1,GAMTMP(4,2),1)
      DEVMAX = 0.0D+00
      DO I = 1, 12
        VAL = ABS((GAMTMP(I,1)-GAMTMP(I,2)))
        IF (VAL.GE.DEVMAX) DEVMAX = VAL
      END DO
      DEVMAX = DEVMAX*4.0D+00*PI*VOLINV
!     DEVMAX = DEVMAX*4.0D+00*PI/VOLPBC
      IF (DEVMAX.LE.THREW) EXIT
      CALL DCOPY(12,GAMTMP(1,2),1,GAMTMP(1,1),1)
    END DO
!
    VAL = -4.0D+00*PI*VOLINV
!   VAL = -4.0D+00*PI/VOLPBC
    DO K = 1, 3
      GAM2(K) = GAM2(K) + VEC(K)*VAL
    END DO
    !! For the charged system term
!   do k = 1, 3
!     gam2st(k,k) = gam2st(k,k) + pi/(ewalpha*ewalpha)*volinv
!   end do
  END IF
!
  if (iat.ne.jat) then
!   call dscal(9,0.5d+00,gam2st,1)
!   if (dftb3) then
!     call dscal(9,0.5d+00,gam3ist,1)
!     call dscal(9,0.5d+00,gam3jst,1)
!   end if
  end if
!
  RETURN
!
END SUBROUTINE DFTB_GAMDER_PBC
!
!-----------------------------------------------------------------------
!*MODULE DFTBPB    *DECK DFTB_FMOSTR
!>
!>     @brief Calculate stress tensor in FMO.
!>
!>     @details Calculate stress tensor in FMO/PBC.
!>
!>     @author Yoshio Nishimoto
!>
SUBROUTINE DFTB_FMOSTR(IDAM,IDAD,IDAT,NBODY)
!
  IMPLICIT NONE
  INTEGER IDAM,IDAD,IDAT,NBODY
  DOUBLE PRECISION DAM,DAD,DAT
!
!     WRITE (IW,'(X,"STRESS TENSOR IN FMOSTR")')
!     DO I = 1, 3
!       WRITE (IW,'(3F20.10)') (STRESS(I,J),J=1,3)
!     END DO
  IF (NBODY.EQ.1) THEN
    DAM = IDAM
    CALL DAXPY(9,DAM,STRESS,1,STRFMO,1)
  ELSE IF (NBODY.EQ.2) THEN
    DAM = IDAM
    DAD = IDAD
    CALL DAXPY(9,DAM+DAD,STRESS,1,STRFMO,1)
  ELSE IF (NBODY.EQ.3) THEN
    DAM = IDAM
    DAD = IDAD
    DAT = IDAT
    CALL DAXPY(9,DAM+DAD+DAT,STRESS,1,STRFMO,1)
  END IF
!
END SUBROUTINE DFTB_FMOSTR
!
!-----------------------------------------------------------------------
!*MODULE DFTBPB    *DECK DFTB_EWALPHA
!>
!>     @brief Calculate alpha for Ewald.
!>
!>     @details Calculate alpha for Ewald PBC.
!>
!>     @author Yoshio Nishimoto
!>
SUBROUTINE DFTB_EWALPHA(NAT,NBODY_)
!
  IMPLICIT NONE
  INTEGER NAT,NBODY_
!
  NBODY = NBODY_
  IF (EWAINP.NE.0.0D+00) RETURN
! IF (MASWRK) WRITE (IW,'(X,"RECOMPUTE ALPHA PARAMETER FOR FMO")')
  EWALPHA=(NAT*PI*PI*PI/(VOLPBC*VOLPBC))**(1.0D+00/6.0D+00)
  IF (MASWRK) WRITE (IW,'(X,"OPTIMUM EWALD PARAMETER = ", F12.5)') EWALPHA
!
  SCALEXP = 1.0D+00/(4.0D+00*EWALPHA*EWALPHA)
!
END SUBROUTINE DFTB_EWALPHA
!
!-----------------------------------------------------------------------
!
!*MODULE DFTBPB    *DECK DFTB_DISP_PBC
!>
!>     @brief Calculate dispersion.
!>
!>     @details Calculate dispersion for PBC.
!>
!>     @author Yoshio Nishimoto
!>
SUBROUTINE DFTB_DISP_PBC(NDER,EDISP,GRAD,HESS)
!
  USE MX_LIMITS, ONLY: MXATM
!
  IMPLICIT NONE
!
  INTEGER MXSPE,NDER
  DOUBLE PRECISION ZERO,ONE,TWO,HALF,EDISP,ETEMP,DFTBDP,DAMPXHE,HUBDER
  DOUBLE PRECISION ZREF,SPNCST,SPE,SKDIM,SKSPIN,QREFL,HUBBL,QREF,HUBB
  DOUBLE PRECISION SKCUT2,C0,C1,C2,CR0,RI,DI,RJ,DJ,RIJ,DIJ,R0,PA,R0A
  DOUBLE PRECISION C6A,PB,R0B,C6B,C6AB,R0AB,R0AB7I,X,ZAN,C,DIST
  DOUBLE PRECISION ALPHA,ALPHAI,ALPHA3,ALPHA6,VOL,CI,CJ,SCAL,RIJ6,BIJ
  DOUBLE PRECISION THRES,DISTI,TMP,R6,R12,DIST2,DISTI2,DISTI4,DISTI6
  DOUBLE PRECISION R5,R10,DIST7,DAMP4,DIST6I,VAL,DIST7I,R7R07,VEXP
  DOUBLE PRECISION DAMP1,DAMP3,VALG,VALH,DAMP2,DEVMAX,BB,BBI,ERFCBB
  DOUBLE PRECISION EXP_BBBB,GMSERFC,DDOT,SCOS,VAL1,VAL2
  INTEGER NSPE,MAXANG,ISPE,IND,IDFTBD,NEEDSK,LSKHTAB,LSKSTAB,LSKGRID
  INTEGER LSKSELF,I,J,LOADFM,LIWRK,LAST,NEED,IAT,NSP1,JAT,NSP2
  INTEGER ITRVEC,ITRVECX,ITRVECY,ITRVECZ,IXVEC,IYVEC,IZVEC,IGZ,JGZ
  INTEGER NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,IAN,IXFTCH,N,NBI,NBJ,KI,NSEQ,KJ
  PARAMETER (MXSPE=10)
  PARAMETER (ZERO=0.0D+00,ONE=1.0D+00,TWO=2.0D+00,HALF=0.5D+00)
!
  COMMON /DFTBPR/ ETEMP,DFTBDP(MXSPE*14),DAMPXHE,HUBDER(MXSPE),     &
                  ZREF(MXATM),SPNCST(6,MXSPE),SPE(MXATM),NSPE,      &
                  MAXANG(MXATM),ISPE(MXATM),IND(MXATM+1),IDFTBD,    &
                  PARAMDIR
  COMMON /DFTBSK/ SKDIM(MXSPE,MXSPE),SKSPIN(MXSPE),QREFL(3,MXSPE),  &
                  HUBBL(3,MXSPE),QREF(MXSPE),HUBB(MXSPE),SKCUT2,    &
                  NEEDSK,LSKHTAB(MXSPE,MXSPE),LSKSTAB(MXSPE,MXSPE), &
                  LSKGRID(MXSPE,MXSPE),LSKSELF(MXSPE)
  COMMON /FMCOM / X(1)
  COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB, &
                  ZAN(MXATM),C(3,MXATM),IAN(MXATM)
!
  DOUBLE PRECISION NEA,NEB
  DOUBLE PRECISION GRAD(3,NAT),HESS(3*NAT,3*NAT)
  DOUBLE PRECISION VALTMP(12,2),RDMAT(3,NSPE,NSPE),VEC(3),VEC1(3)
  DOUBLE PRECISION VEC2(3),VECREC(4),PARAMDIR
  LOGICAL FMO
!
  FMO = .FALSE.
  IF (NBODY.NE.0) FMO = .TRUE.
!
! Only 3D-UFF dispersion is computed using Ewald sum.
! Other dimensional or dispersion is evaluated directly with
! a pruned threshold.
!
  IF (NDER.EQ.0) EDISP = ZERO
  IF (IDFTBD.EQ.1) THEN
    C0 = 1.584D+01 !! 396/25
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
!
        C6AB = (2.0D+00*C6A*C6B*PA*PB) / (PA*PA*C6B+PB*PB*C6A)
        R0AB = (R0A**3 + R0B**3) / (R0A**2 + R0B**2)
        R0AB7I = 1.0D+00 / R0AB**7
        RDMAT(1,I,J) = C6AB
        RDMAT(2,I,J) = R0AB7I
        RDMAT(1,J,I) = C6AB
        RDMAT(2,J,I) = R0AB7I
      END DO
    END DO
  ELSE IF (IDFTBD.EQ.4) THEN
    CALL VALFM(LOADFM)
    LIWRK = LOADFM  + 1
    LAST  = LIWRK + NAT
    NEED  = LAST - LOADFM - 1
    CALL GETFM(NEED)
!
    !! COUNT NUMBER OF BONDS FOR ALL ATOMS
    CALL VICLR(X(LIWRK),1,NAT)
    DO IAT = 1, NAT
      NSP1 = ISPE(IAT)
      R0A = DFTBDP(NSP1*14-13) !! COVALENT RADIUS
      DO JAT = 1, IAT
        NSP2 = ISPE(JAT)
        R0B = DFTBDP(NSP2*14-13)
        ITRVEC = 1
        IF (NDIMPBC.GE.1) ITRVECX = ITRVEC
        IF (NDIMPBC.GE.2) ITRVECY = ITRVEC
        IF (NDIMPBC.GE.3) ITRVECZ = ITRVEC
        DO IXVEC = -ITRVECX, ITRVECX
        DO IYVEC = -ITRVECY, ITRVECY
        DO IZVEC = -ITRVECZ, ITRVECZ
          VEC(1) = C(1,IAT) - C(1,JAT) &
                 - PBCBOX(1,1)*IXVEC   &
                 - PBCBOX(2,1)*IYVEC   &
                 - PBCBOX(3,1)*IZVEC
          VEC(2) = C(2,IAT) - C(2,JAT) &
                 - PBCBOX(1,2)*IXVEC   &
                 - PBCBOX(2,2)*IYVEC   &
                 - PBCBOX(3,2)*IZVEC
          VEC(3) = C(3,IAT) - C(3,JAT) &
                 - PBCBOX(1,3)*IXVEC   &
                 - PBCBOX(2,3)*IYVEC   &
                 - PBCBOX(3,3)*IZVEC
          DIST = VEC(1)*VEC(1)+VEC(2)*VEC(2)+VEC(3)*VEC(3)
          IF (DIST.LT.(R0A+R0B)**2 .AND. DIST.GT.1.0D-08) THEN
            CALL IXSTOR(X(LIWRK),IAT,IXFTCH(X(LIWRK),IAT)+1)
            IF (IAT.NE.JAT) &
              CALL IXSTOR(X(LIWRK),JAT,IXFTCH(X(LIWRK),JAT)+1)
          END IF
        END DO
        END DO
        END DO
      END DO
    END DO
    DO IAT = 1, NAT
      IF (IXFTCH(X(LIWRK),IAT).GT.5) CALL IXSTOR(X(LIWRK),IAT,5)
    END DO
  ELSE
    if (maswrk) write (iw,*) "no such dispersion"
    call abrt
  END IF
!
  ALPHA  = EWALPHA
! alpha = 1.0d-03
  ALPHAI = ONE/ALPHA
  ALPHA3 = ALPHA*ALPHA*ALPHA
  ALPHA6 = ALPHA3*ALPHA3
  VOL = VOLPBC
  IF (NDER.EQ.0) N = 1
  IF (NDER.EQ.1) N = 12
  IF (NDER.EQ.2) N = 9
! THRES = 1.0D-15
! IF (NDIMPBC.LE.2 .OR. IDFTBD.NE.1) THRES=1.0D-10
  THRES = THREWD
!
  DO IAT = 1, NAT
    NSP1 = ISPE(IAT)
    IF (FMO) CI = ZREF(IAT)/QREF(NSP1)
    VEC1(1) = C(1,IAT)
    VEC1(2) = C(2,IAT)
    VEC1(3) = C(3,IAT)
    DO JAT = 1, IAT
      NSP2 = ISPE(JAT)
      IF (FMO) CJ = ZREF(JAT)/QREF(NSP2)
      VEC2(1) = C(1,JAT)
      VEC2(2) = C(2,JAT)
      VEC2(3) = C(3,JAT)
!
      SCAL = 1.0D+00
      IF (IAT.EQ.JAT) SCAL = 0.5D+00
      IGZ = 3*(IAT-1)
      JGZ = 3*(JAT-1)
      IF (IDFTBD.EQ.1) THEN
        R0  = RDMAT(1,NSP1,NSP2)
        RIJ = RDMAT(2,NSP1,NSP2)
        DIJ = RDMAT(3,NSP1,NSP2)
        RIJ6= RIJ**6
        BIJ = -TWO*DIJ*RIJ6
!       IF (FMO) DIJ = DIJ * CI * CJ
      ELSE IF (IDFTBD.EQ.2) THEN
        C6AB   = RDMAT(1,NSP1,NSP2)
        R0AB7I = RDMAT(2,NSP1,NSP2)
!       IF (FMO) C6AB = C6AB * CI * CJ
      ELSE IF (IDFTBD.EQ.4) THEN
        NBI  = IXFTCH(X(LIWRK),IAT)
        PA   = DFTBDP(NSP1*14-12+NBI)
        R0A  = DFTBDP(NSP1*14- 6+NBI)
        NEA  = DFTBDP(NSP1*14)
        C6A  = 7.5D-01*SQRT(NEA*PA*PA*PA)
        NBJ  = IXFTCH(X(LIWRK),JAT)
        PB   = DFTBDP(NSP2*14-12+NBJ)
        R0B  = DFTBDP(NSP2*14- 6+NBJ)
        NEB  = DFTBDP(NSP2*14)
        C6B  = 7.5D-01*SQRT(NEB*PB*PB*PB)
        C6AB = (2.0D+00*C6A*C6B*PA*PB) / (PA*PA*C6B+PB*PB*C6A)
        R0AB = (R0A**3 + R0B**3) / (R0A**2 + R0B**2)
        R0AB7I = 1.0D+00 / R0AB**7
!       IF (FMO) C6AB = C6AB * CI * CJ
      END IF
!
!     ----- REAL-SPACE
!
      CALL VCLR(VALTMP,1,2*12)
      DO ITRVEC = 0, 100
        IF (NDIMPBC.GE.1) ITRVECX = ITRVEC
        IF (NDIMPBC.GE.2) ITRVECY = ITRVEC
        IF (NDIMPBC.GE.3) ITRVECZ = ITRVEC
        DO IXVEC = -ITRVECX, ITRVECX
        DO IYVEC = -ITRVECY, ITRVECY
        DO IZVEC = -ITRVECZ, ITRVECZ
        IF (ABS(IXVEC).NE.ITRVEC .AND. ABS(IYVEC).NE.ITRVEC .AND. &
            ABS(IZVEC).NE.ITRVEC) CYCLE
!
        VEC(1) = VEC1(1) - VEC2(1) &
               - PBCBOX(1,1)*IXVEC &
               - PBCBOX(2,1)*IYVEC &
               - PBCBOX(3,1)*IZVEC
        VEC(2) = VEC1(2) - VEC2(2) &
               - PBCBOX(1,2)*IXVEC &
               - PBCBOX(2,2)*IYVEC &
               - PBCBOX(3,2)*IZVEC
        VEC(3) = VEC1(3) - VEC2(3) &
               - PBCBOX(1,3)*IXVEC &
               - PBCBOX(2,3)*IYVEC &
               - PBCBOX(3,3)*IZVEC
        DIST = SQRT(VEC(1)*VEC(1)+VEC(2)*VEC(2)+VEC(3)*VEC(3))
        IF (DIST.LE.1.0D-08) CYCLE
        DISTI  = ONE/DIST
!
        IF (NDER.EQ.0) THEN
          IF (IDFTBD.EQ.1) THEN
            IF (DIST.GE.R0) THEN !! ATTRACTIVE POTENTIAL
              TMP = RIJ*DISTI
              R6 = TMP**6
              R12 = R6*R6
              IF (NDIMPBC.LE.2) THEN
                VALTMP(1,1) = VALTMP(1,1) + DIJ*(-TWO*R6+R12)
                CYCLE
              END IF
              DIST = DIST*ALPHA
              DIST2 = DIST*DIST
              DISTI2 = ONE/DIST2
              DISTI4 = DISTI2*DISTI2
              DISTI6 = DISTI4*DISTI2
              VALTMP(1,1) = VALTMP(1,1) + DIJ*R12 &
                + BIJ*ALPHA6*(DISTI6+DISTI4+HALF*DISTI2)*EXP(-DIST2)
            ELSE !! REPULSIVE POTENTIAL
              TMP = DIST/RIJ
              R5 = TMP**5
              R10 = R5**2
              VALTMP(1,1) = VALTMP(1,1) + DIJ*(C0-C1*R5+C2*R10)
!
              !! Compensate the R^-6 term
              TMP = RIJ*DISTI
              R6 = TMP**6
              VALTMP(1,1) = VALTMP(1,1) + 2.0D+00*DIJ*R6
!
              DIST = DIST*ALPHA
              DIST2 = DIST*DIST
              DISTI2 = ONE/DIST2
              DISTI4 = DISTI2*DISTI2
              DISTI6 = DISTI4*DISTI2
              VALTMP(1,1) = VALTMP(1,1) &
                + BIJ*ALPHA6*(DISTI6+DISTI4+HALF*DISTI2)*EXP(-DIST2)
            END IF
          ELSE IF (IDFTBD.EQ.2 .OR. IDFTBD.EQ.4) THEN
            DIST7  = DIST**7
            DAMP4 = (ONE-EXP(-3.0D+00*R0AB7I*DIST7))**4
            DIST6I = DISTI**6
            VALTMP(1,1) = VALTMP(1,1) - DAMP4*C6AB*DIST6I
          END IF
        ELSE IF (NDER.EQ.1) THEN
          !! GRADIENT
          IF (IDFTBD.EQ.1) THEN
            IF (DIST.GE.R0) THEN !! ATTRACTIVE POTENTIAL
              TMP = RIJ*DISTI
              R6 = TMP**6
              R12 = R6*R6
              IF (NDIMPBC.LE.2) THEN
                VAL = 1.2D+01*DIJ*(R6-R12)*DISTI
              ELSE
                VAL = -1.2D+01*DIJ*R12*DISTI
!
                DIST   = DIST*ALPHA
                DIST2  = DIST*DIST
                DISTI2 = ONE/DIST2
                DISTI4 = DISTI2*DISTI2
                DISTI6 = DISTI4*DISTI2
                VAL = VAL &
!               - BIJ*ALPHA6*(6.0D+00*DISTI6+4.0D+00*DISTI4+DISTI2) &
!                *EXP(-DIST2)*DISTI &
!               - 2.0D+00*ALPHA*DIST &
!                *BIJ*ALPHA6*(DISTI6+DISTI4+HALF*DISTI2)*EXP(-DIST2)
                    - BIJ*ALPHA6*(6.0D+00*DISTI6+6.0D+00*DISTI4 &
                                 +3.0d+00*DISTI2+1.0D+00) &
                      *EXP(-DIST2)*DISTI
              END IF
            ELSE !! REPULSIVE POTENTIAL
              TMP = DIST/RIJ
              R5 = TMP**5
              R10 = R5**2
              VAL = DIJ*(-C1*R5*5.0D+00+C2*R10*1.0D+01)*DISTI
!
              !! Compensate the R^-6 term
              TMP = RIJ*DISTI
              R6 = TMP**6
              VAL = VAL - 1.2D+01*DIJ*R6*DISTI
!
              DIST   = DIST*ALPHA
              DIST2  = DIST*DIST
              DISTI2 = ONE/DIST2
              DISTI4 = DISTI2*DISTI2
              DISTI6 = DISTI4*DISTI2
              VAL = VAL &
                  - BIJ*ALPHA6*(6.0D+00*DISTI6+6.0D+00*DISTI4 &
                               +3.0d+00*DISTI2+1.0D+00) &
                    *EXP(-DIST2)*DISTI
            END IF
          ELSE IF (IDFTBD.EQ.2 .OR. IDFTBD.EQ.4) THEN
            DIST6I = DISTI**6
            DIST7  = DIST**7
            DIST7I = DIST6I*DISTI
            R7R07 = DIST7*R0AB7I
            VEXP = EXP(-3.0D+00*R7R07)
            DAMP1 = 1.0D+00 - VEXP
            DAMP3 = DAMP1**3
            VAL = 8.4D+01*VEXP*R7R07-6.0D+00*DAMP1
            VAL = -VAL*C6AB*DAMP3*DIST7I
          END IF
!
          !! STRESS TENSOR
          DO KI = 1, 3
            VALTMP(KI,1) = VALTMP(KI,1) + VEC(KI)*VAL*DISTI
          END DO
          NSEQ = 0
          DO KJ = 1, 3
            DO KI = 1, 3
              NSEQ = NSEQ + 1
                VALTMP(3+NSEQ,1) = VALTMP(3+NSEQ,1) &
                  + VEC(KI)*VEC(KJ)*VAL*DISTI
            END DO
          END DO
        ELSE IF (NDER.EQ.2) THEN
          IF (IDFTBD.EQ.1) THEN
            IF (DIST.GE.R0) THEN
              TMP = RIJ*DISTI
              R6 = TMP**6
              R12 = R6**2
              IF (NDIMPBC.LE.2) THEN
                VALG = 1.2D+01*DIJ*(R6-R12)*DISTI
                VALH = 1.2D+01*DIJ*(-7.0D+00*R6+1.3D+01*R12)*DISTI*DISTI
              ELSE
                DIST = DIST*ALPHA
                DIST2 = DIST*DIST
                DISTI2 = ONE/DIST2
                DISTI4 = DISTI2*DISTI2
                DISTI6 = DISTI4*DISTI2
                VALG = -1.2D+01*DIJ*R12*DISTI &
                     - BIJ*ALPHA6*(6.0D+00*DISTI6+6.0D+00*DISTI4 &
                                  +3.0d+00*DISTI2+1.0D+00) &
                       *EXP(-DIST2)*DISTI
                VALH = 1.56D+02*DIJ*R12*DISTI*DISTI &
                     + BIJ*ALPHA6*(4.2D+01*DISTI6+4.2D+01*DISTI4 &
                         +2.1D+01*DISTI2+7.0D+00+2.0D+00*DIST2) &
                       *DISTI*DISTI*EXP(-DIST2)
              END IF
            ELSE
              TMP = DIST/RIJ
              R5 = TMP**5
              R10 = R5**2
              VALG = DIJ*(-5.0D+00*C1*R5+1.0D+01*C2*R10)*DISTI
              VALH = DIJ*(-2.0D+01*C1*R5+9.0D+01*C2*R10)*DISTI*DISTI
              !! Compensate the R^-6 term
              TMP = RIJ/DIST
              R6 = TMP**6
!
              DIST = DIST*ALPHA
              DIST2 = DIST*DIST
              DISTI2 = ONE/DIST2
              DISTI4 = DISTI2*DISTI2
              DISTI6 = DISTI4*DISTI2
              VALG = VALG - 1.2D+01*DIJ*R6*DISTI &
                   - BIJ*ALPHA6*(6.0D+00*DISTI6+6.0D+00*DISTI4 &
                                +3.0d+00*DISTI2+1.0D+00) &
                     *EXP(-DIST2)*DISTI
              VALH = VALH + 8.4D+01*DIJ*R6*DISTI*DISTI &
                   + BIJ*ALPHA6*(4.2D+01*DISTI6+4.2D+01*DISTI4 &
                       +2.1D+01*DISTI2+7.0D+00+2.0D+00*DIST2) &
                     *DISTI*DISTI*EXP(-DIST2)
            END IF
          ELSE IF (IDFTBD.EQ.2 .OR. IDFTBD.EQ.4) THEN
            DIST7  = DIST**7
            DIST7I = DISTI**7
            R7R07 = DIST7*R0AB7I
            VEXP = EXP(-3.0D+00*R7R07)
            DAMP1 = ONE - VEXP
            DAMP2 = DAMP1 * DAMP1
            DAMP3 = DAMP2 * DAMP1
            VAL = 8.4D+01*VEXP*R7R07-6.0D+00*DAMP1
            VALG = -VAL*C6AB*DAMP3*DIST7I
            VAL = 2.52D+02*DAMP1*R7R07*(2.0D+00-7.0D+00*R7R07)*VEXP &
                + 5.292D+03*R7R07*R7R07*VEXP*VEXP &
                - 1.008D+03*DAMP1*R7R07*VEXP + 4.2D+01*DAMP2
            VALH = -VAL*C6AB*DAMP2*DIST7I*DISTI
          END IF
!
          DO KI = 1, 3
            DO KJ = 1, 3
              VAL = VEC(KI)*VEC(KJ)*(VALH-DISTI*VALG)*DISTI*DISTI
              IF (KI.EQ.KJ) VAL = VAL + DISTI*VALG
              VALTMP(KI+3*(KJ-1),1) = VALTMP(KI+3*(KJ-1),1) + VAL
            END DO
          END DO
        END IF

        END DO
        END DO
        END DO
!
        !! CHECK CONVERGENCE
        DEVMAX = ZERO
        DO I = 1, N
          VAL = ABS((VALTMP(I,1)-VALTMP(I,2)))
          IF (VAL.GE.DEVMAX) DEVMAX = VAL
        END DO
        IF (DEVMAX.LE.THRES .AND. ITRVEC.GE.1) EXIT
        CALL DCOPY(N,VALTMP(1,1),1,VALTMP(1,2),1)
      END DO
!     write (*,*) "converged short-range at", itrvec
!
      !! Finally, add short-range contributions
      CALL DSCAL(N,SCAL,VALTMP,1)
      IF (FMO) CALL DSCAL(N,CI*CJ,VALTMP,1)
      IF (NDER.EQ.0) THEN
        EDISP = EDISP + VALTMP(1,1)
      ELSE IF (NDER.EQ.1) THEN
        DO KI = 1, 3
          GRAD(KI,IAT) = GRAD(KI,IAT) + VALTMP(KI,1)
          GRAD(KI,JAT) = GRAD(KI,JAT) - VALTMP(KI,1)
        END DO
        CALL DAXPY(9,TWO,VALTMP(4,1),1,STRESS,1)
      ELSE IF (NDER.EQ.2) THEN
        DO KI = 1, 3
          DO KJ = 1, 3
            VAL = VALTMP(KI+3*(KJ-1),1)
            HESS(IGZ+KI,IGZ+KJ) = HESS(IGZ+KI,IGZ+KJ) + VAL
            HESS(JGZ+KI,JGZ+KJ) = HESS(JGZ+KI,JGZ+KJ) + VAL
            HESS(IGZ+KI,JGZ+KJ) = HESS(IGZ+KI,JGZ+KJ) - VAL
            HESS(JGZ+KI,IGZ+KJ) = HESS(JGZ+KI,IGZ+KJ) - VAL
          END DO
        END DO
      END IF
!
      IF (NDIMPBC.LE.2 .OR. IDFTBD.NE.1) CYCLE
!
!     ----- RECIPROCAL-SPACE
!
      VEC(1) = VEC1(1) - VEC2(1)
      VEC(2) = VEC1(2) - VEC2(2)
      VEC(3) = VEC1(3) - VEC2(3)
      CALL VCLR(VALTMP,1,2*12)
      DO ITRVEC = 1, 100
        IF (NDIMPBC.GE.1) ITRVECX = ITRVEC
        IF (NDIMPBC.GE.2) ITRVECY = ITRVEC
        IF (NDIMPBC.GE.3) ITRVECZ = ITRVEC
        DO IXVEC = -ITRVECX, ITRVECX
        DO IYVEC = -ITRVECY, ITRVECY
        DO IZVEC = -ITRVECZ, ITRVECZ
        IF (ABS(IXVEC).NE.ITRVEC .AND. ABS(IYVEC).NE.ITRVEC .AND. &
            ABS(IZVEC).NE.ITRVEC) CYCLE
!     
        VECREC(1) = PBCBOX(1,4)*IXVEC &
                  + PBCBOX(2,4)*IYVEC &
                  + PBCBOX(3,4)*IZVEC
        VECREC(2) = PBCBOX(1,5)*IXVEC &
                  + PBCBOX(2,5)*IYVEC &
                  + PBCBOX(3,5)*IZVEC
        VECREC(3) = PBCBOX(1,6)*IXVEC &
                  + PBCBOX(2,6)*IYVEC &
                  + PBCBOX(3,6)*IZVEC
        VECREC(4) = SQRT(VECREC(1)**2 + VECREC(2)**2 + VECREC(3)**2)
        BB = HALF*VECREC(4)*ALPHAI
        BBI= ONE/BB
!       if(IERFTP.NE.0) then
          ERFCBB=GMSERFC(BB)
!        else
!          ERFCBB=ERFC(BB)
!        endif
        EXP_BBBB=EXP(-BB*BB)
!
        IF (NDER.EQ.0) THEN
          VALTMP(1,1) = VALTMP(1,1) &
            + COS(DDOT(3,VEC,1,VECREC,1)) &
              *VECREC(4)*VECREC(4)*VECREC(4) &
              *(PISQ*ERFCBB+(HALF*BBI*BBI-ONE)*BBI*EXP_BBBB)
        ELSE IF (NDER.EQ.1) THEN
          !! GRADIENT
          VAL = -SIN(DDOT(3,VEC,1,VECREC,1)) &
                *VECREC(4)*VECREC(4)*VECREC(4) &
                *(PISQ*ERFCBB+(HALF*BBI*BBI-ONE)*BBI*EXP_BBBB)
          DO KI = 1, 3
            VALTMP(KI,1) = VALTMP(KI,1) + VAL*VECREC(KI)
          END DO
!
          !! STRESS TENSOR
          SCOS = COS(DDOT(3,VEC,1,VECREC,1))
          VAL1 = VECREC(4)*VECREC(4)*VECREC(4) &
                *(PISQ*ERFCBB+(HALF*BBI*BBI-ONE)*BBI*EXP_BBBB)
          VAL2 = 3.0D+00*VECREC(4)*(PISQ*ERFCBB-BBI*EXP_BBBB)
          NSEQ = 0
          DO KJ = 1, 3
            DO KI = 1, 3
              NSEQ = NSEQ + 1
              VAL = VAL2*VECREC(KI)*VECREC(KJ)
              IF (KI.EQ.KJ) VAL = VAL + VAL1
              VALTMP(3+NSEQ,1) = VALTMP(3+NSEQ,1) - SCOS*VAL
            END DO
          END DO
        ELSE IF (NDER.EQ.2) THEN
          VAL = -COS(DDOT(3,VEC,1,VECREC,1)) &
                *VECREC(4)*VECREC(4)*VECREC(4) &
                *(PISQ*ERFCBB+(HALF*BBI*BBI-ONE)*BBI*EXP_BBBB)
          DO KI = 1, 3
            DO KJ = 1, 3
              VALTMP(KI+3*(KJ-1),1) = VALTMP(KI+3*(KJ-1),1) &
                + VAL*VECREC(KI)*VECREC(KJ)
            END DO
          END DO
        END IF
!
        END DO
        END DO
        END DO
!
        !! CHECK CONVERGENCE
        DEVMAX = ZERO
        DO I = 1, N
          VAL = ABS((VALTMP(I,1)-VALTMP(I,2)))
          IF (VAL.GE.DEVMAX) DEVMAX = VAL
        END DO
        IF (DEVMAX.LE.THRES) EXIT
        CALL DCOPY(N,VALTMP(1,1),1,VALTMP(1,2),1)
      END DO
!     write (*,*) "converged  long-range at", itrvec
!
      !! Finally, add reciprocal-space contributions
      CALL DSCAL(N,SCAL*BIJ*PI*PISQ/(1.2D+01*VOL),VALTMP,1)
      IF (FMO) CALL DSCAL(N,CI*CJ,VALTMP,1)
      IF (NDER.EQ.0) THEN
        EDISP = EDISP + VALTMP(1,1)
      ELSE IF (NDER.EQ.1) THEN
        DO KI = 1, 3
          GRAD(KI,IAT) = GRAD(KI,IAT) + VALTMP(KI,1)
          GRAD(KI,JAT) = GRAD(KI,JAT) - VALTMP(KI,1)
        END DO
        CALL DAXPY(9,TWO,VALTMP(4,1),1,STRESS,1)
      ELSE IF (NDER.EQ.2) THEN
        DO KI = 1, 3
          DO KJ = 1, 3
            VAL = VALTMP(KI+3*(KJ-1),1)
            HESS(IGZ+KI,IGZ+KJ) = HESS(IGZ+KI,IGZ+KJ) + VAL
            HESS(JGZ+KI,JGZ+KJ) = HESS(JGZ+KI,JGZ+KJ) + VAL
            HESS(IGZ+KI,JGZ+KJ) = HESS(IGZ+KI,JGZ+KJ) - VAL
            HESS(JGZ+KI,IGZ+KJ) = HESS(JGZ+KI,IGZ+KJ) - VAL
          END DO
        END DO
      END IF
!
!     ----- K = 0 FOR 1/R6 TERM
!
      IF (FMO) BIJ = BIJ*CI*CJ
      IF (NDER.EQ.0) THEN
        EDISP = EDISP + BIJ*ALPHA3*PI*PISQ/(3.0D+00*VOL)*SCAL
        IF (IAT.EQ.JAT) THEN
          EDISP = EDISP - BIJ*ALPHA6/1.2D+01
        END IF
      ELSE IF (NDER.EQ.1) THEN
        DO KI = 1, 3
          STRESS(KI,KI) = STRESS(KI,KI) &
            - BIJ*ALPHA3*PI*PISQ/(3.0D+00*VOL)*SCAL*TWO
        END DO
      ELSE IF (NDER.EQ.2) THEN
        !! Nothing
      END IF
    END DO
  END DO
!
  IF (IDFTBD.EQ.4) CALL RETFM(NEED)
!
  RETURN
!
END SUBROUTINE DFTB_DISP_PBC
!
!-----------------------------------------------------------------------
!
!*MODULE DFTBPB    *DECK DFTB_ESDIM_DCPBC
!>
!>     @brief Calculate dispersion.
!>
!>     @details Calculate dispersion for ES dimers in FMO.
!>
!>     @author Yoshio Nishimoto
!>
SUBROUTINE DFTB_ESDIM_DCPBC(IFG,JFG,NDER,MODGAM,DCTMP,FMODE,EDISP,SCAL)
!
  USE MX_LIMITS, ONLY: MXATM
!
  IMPLICIT NONE
!
  INTEGER MXSPE
  PARAMETER (MXSPE=10)
  DOUBLE PRECISION ZERO,ONE,TWO,HALF
  PARAMETER (ZERO=0.0D+00,ONE=1.0D+00,TWO=2.0D+00,HALF=0.5D+00)
!
  DOUBLE PRECISION ETEMP,DFTBDP,DAMPXHE,HUBDER,ZREF,SPNCS,SPE,SKDIM
  DOUBLE PRECISION SKSPIN,QREFL,HUBBL,QREF,HUBB,SKCUT2,SPNCST
  INTEGER NSPE,MAXANG,ISPE,IND,IDFTBD,NEEDSK,LSKHTAB,LSKSTAB,LSKGRID
  INTEGER LSKSELF,NATI,NATJ,INDI,INDJ,IFG,JFG,NDER,MODGAM,I,J,N
  INTEGER IAT,IAG,ISP,JAT,JAG,JSP,NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,IAN
  INTEGER NBI,NBJ,ITRVEC,ITRVECX,ITRVECY,ITRVECZ,IXVEC,IYVEC,IZVEC
  INTEGER KI,NSEQ,KJ
  COMMON /DFTBPR/ ETEMP,DFTBDP(MXSPE*14),DAMPXHE,HUBDER(MXSPE),     &
                  ZREF(MXATM),SPNCST(6,MXSPE),SPE(MXATM),NSPE,      &
                  MAXANG(MXATM),ISPE(MXATM),IND(MXATM+1),IDFTBD,    &
                  PARAMDIR
  COMMON /DFTBSK/ SKDIM(MXSPE,MXSPE),SKSPIN(MXSPE),QREFL(3,MXSPE),  &
                  HUBBL(3,MXSPE),QREF(MXSPE),HUBB(MXSPE),SKCUT2,    &
                  NEEDSK,LSKHTAB(MXSPE,MXSPE),LSKSTAB(MXSPE,MXSPE), &
                  LSKGRID(MXSPE,MXSPE),LSKSELF(MXSPE)
  COMMON /FMCOM / X(1)
  COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB, &
                  ZAN(MXATM),C(3,MXATM),IAN(MXATM)
!
  DOUBLE PRECISION NEA,NEB,PARAMDIR
  DOUBLE PRECISION VALTMP(13,2),RDMAT(3,NSPE,NSPE),VEC(3),VEC1(3)
  DOUBLE PRECISION VEC2(3),VECREC(4),FMODE(3,NATFMO,*),DCTMP(*)
  DOUBLE PRECISION C0,C1,C2,CR0,RI,DI,RJ,DJ,RIJ,DIJ,R0,PA,R0A,C6A,PB
  DOUBLE PRECISION R0B,C6B,C6AB,R0AB,R0AB7I,ALPHA,ALPHAI,ALPHA3,ALPHA6
  DOUBLE PRECISION VOL,THRES,ZAN,C,CI,CJ,RIJ6,BIJ,DIST,DISTI,TMP,R6,R12
  DOUBLE PRECISION DIST2,DISTI2,DISTI4,DISTI6,R5,R10,DIST7,DAMP4,DIST6I
  DOUBLE PRECISION VAL,DIST7I,R7R07,VEXP,DAMP1,DAMP3,DEVMAX,EDISP,X,VAL2
  DOUBLE PRECISION SCAL,BB,BBI,ERFCBB,GMSERFC,EXP_BBBB,DDOT,SCOS,VAL1
!
  NATI = NATFRG(IFG)
  NATJ = NATFRG(JFG)
  INDI = INDFRG(IFG)
  INDJ = INDFRG(JFG)
  IF (NDER.EQ.1) CALL DCOPY(9,0.0D+00,0,STRESS,1)
!
! Only 3D-UFF dispersion is computed using Ewald sum.
! Other dimensional or dispersion is evaluated directly with
! a pruned threshold.
!
  CALL DCOPY(3*NSPE*NSPE,DCTMP,1,RDMAT,1)
  IF (IDFTBD.EQ.1) THEN
    C0 = 1.584D+01 !! 396/25
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
!
        C6AB = (2.0D+00*C6A*C6B*PA*PB) / (PA*PA*C6B+PB*PB*C6A)
        R0AB = (R0A**3 + R0B**3) / (R0A**2 + R0B**2)
        R0AB7I = 1.0D+00 / R0AB**7
        RDMAT(1,I,J) = C6AB
        RDMAT(2,I,J) = R0AB7I
        RDMAT(1,J,I) = C6AB
        RDMAT(2,J,I) = R0AB7I
      END DO
    END DO
  ELSE IF (IDFTBD.EQ.4) THEN
  ! LIWRK = LOADFM  + 1
  ! LAST  = LIWRK + NAT
  ! NEED  = LAST - LOADFM - 1
  ! CALL GETFM(NEED)
!
    !! COUNT NUMBER OF BONDS FOR ALL ATOMS
  ! CALL VICLR(X(LIWRK),1,NAT)
  ! DO IAT = 1, NAT
  !   NSP1 = ISPE(IAT)
  !   R0A = DFTBDP(NSP1*14-13) !! COVALENT RADIUS
  !   DO JAT = 1, IAT
  !     NSP2 = ISPE(JAT)
  !     R0B = DFTBDP(NSP2*14-13)
  !     ITRVEC = 1
  !     IF (NDIMPBC.GE.1) ITRVECX = ITRVEC
  !     IF (NDIMPBC.GE.2) ITRVECY = ITRVEC
  !     IF (NDIMPBC.GE.3) ITRVECZ = ITRVEC
  !     DO IXVEC = -ITRVECX, ITRVECX
  !     DO IYVEC = -ITRVECY, ITRVECY
  !     DO IZVEC = -ITRVECZ, ITRVECZ
  !       VEC(1) = C(1,IAT) - C(1,JAT) &
  !              - PBCBOX(1,1)*IXVEC   &
  !              - PBCBOX(2,1)*IYVEC   &
  !              - PBCBOX(3,1)*IZVEC
  !       VEC(2) = C(2,IAT) - C(2,JAT) &
  !              - PBCBOX(1,2)*IXVEC   &
  !              - PBCBOX(2,2)*IYVEC   &
  !              - PBCBOX(3,2)*IZVEC
  !       VEC(3) = C(3,IAT) - C(3,JAT) &
  !              - PBCBOX(1,3)*IXVEC   &
  !              - PBCBOX(2,3)*IYVEC   &
  !              - PBCBOX(3,3)*IZVEC
  !       DIST = VEC(1)*VEC(1)+VEC(2)*VEC(2)+VEC(3)*VEC(3)
  !       IF (DIST.LT.(R0A+R0B)**2 .AND. DIST.GT.1.0D-08) THEN
  !         CALL IXSTOR(X(LIWRK),IAT,IXFTCH(X(LIWRK),IAT)+1)
  !         IF (IAT.NE.JAT) &
  !           CALL IXSTOR(X(LIWRK),JAT,IXFTCH(X(LIWRK),JAT)+1)
  !       END IF
  !     END DO
  !     END DO
  !     END DO
  !   END DO
  ! END DO
  ! DO IAT = 1, NAT
  !   IF (IXFTCH(X(LIWRK),IAT).GT.5) CALL IXSTOR(X(LIWRK),IAT,5)
  ! END DO
  ELSE
    if (maswrk) write (iw,*) "no such dispersion"
    call abrt
  END IF
!
  ALPHA  = EWALPHA
! alpha = 1.0d-03
  ALPHAI = ONE/ALPHA
  ALPHA3 = ALPHA*ALPHA*ALPHA
  ALPHA6 = ALPHA3*ALPHA3
  VOL = VOLPBC
  N = 1
  IF (NDER.EQ.1) N = 12+1
! THRES = 1.0D-15
! IF (NDIMPBC.LE.2 .OR. IDFTBD.NE.1) THRES=1.0D-10
  THRES = THREWD
!
! write (*,*) "asdf"
  DO IAT = 1, NATI
    IAG = IATFRG(INDI+IAT-1)
    ISP = ISPEFMO(IAG)
    CI = ZREFFMO(INDI+IAT-1)/QREF(ISP)
    VEC1(1) = FMOC(1,IAG)
    VEC1(2) = FMOC(2,IAG)
    VEC1(3) = FMOC(3,IAG)
    DO JAT = 1, NATJ
      JAG = IATFRG(INDJ+JAT-1)
      JSP = ISPEFMO(JAG)
      CJ = ZREFFMO(INDJ+JAT-1)/QREF(JSP)
      VEC2(1) = FMOC(1,JAG)
      VEC2(2) = FMOC(2,JAG)
      VEC2(3) = FMOC(3,JAG)
!
      IF (IDFTBD.EQ.1) THEN
        R0  = RDMAT(1,ISP,JSP)
        RIJ = RDMAT(2,ISP,JSP)
        DIJ = RDMAT(3,ISP,JSP)
        RIJ6= RIJ**6
        BIJ = -TWO*DIJ*RIJ6
      ELSE IF (IDFTBD.EQ.2) THEN
        C6AB   = RDMAT(1,ISP,JSP)
        R0AB7I = RDMAT(2,ISP,JSP)
      ELSE IF (IDFTBD.EQ.4) THEN
        NBI  = NBOND(IAG)
        PA   = DFTBDP(ISP*14-12+NBI)
        R0A  = DFTBDP(ISP*14- 6+NBI)
        NEA  = DFTBDP(ISP*14)
        C6A  = 7.5D-01*SQRT(NEA*PA*PA*PA)
        NBJ  = NBOND(JAG)
        PB   = DFTBDP(JSP*14-12+NBJ)
        R0B  = DFTBDP(JSP*14- 6+NBJ)
        NEB  = DFTBDP(JSP*14)
        C6B  = 7.5D-01*SQRT(NEB*PB*PB*PB)
        if (pa*pa*c6b+pb*pb*c6a.eq.0.0d+00) cycle
        C6AB = (2.0D+00*C6A*C6B*PA*PB) / (PA*PA*C6B+PB*PB*C6A)
        R0AB = (R0A**3 + R0B**3) / (R0A**2 + R0B**2)
        R0AB7I = 1.0D+00 / R0AB**7
      END IF
!
!     ----- REAL-SPACE
!
      CALL VCLR(VALTMP,1,2*13)
      DO ITRVEC = 0, 100
        IF (NDIMPBC.GE.1) ITRVECX = ITRVEC
        IF (NDIMPBC.GE.2) ITRVECY = ITRVEC
        IF (NDIMPBC.GE.3) ITRVECZ = ITRVEC
        DO IXVEC = -ITRVECX, ITRVECX
        DO IYVEC = -ITRVECY, ITRVECY
        DO IZVEC = -ITRVECZ, ITRVECZ
        IF (ABS(IXVEC).NE.ITRVEC .AND. ABS(IYVEC).NE.ITRVEC .AND. &
            ABS(IZVEC).NE.ITRVEC) CYCLE
!
        VEC(1) = VEC1(1) - VEC2(1) &
               - PBCBOX(1,1)*IXVEC &
               - PBCBOX(2,1)*IYVEC &
               - PBCBOX(3,1)*IZVEC
        VEC(2) = VEC1(2) - VEC2(2) &
               - PBCBOX(1,2)*IXVEC &
               - PBCBOX(2,2)*IYVEC &
               - PBCBOX(3,2)*IZVEC
        VEC(3) = VEC1(3) - VEC2(3) &
               - PBCBOX(1,3)*IXVEC &
               - PBCBOX(2,3)*IYVEC &
               - PBCBOX(3,3)*IZVEC
        DIST = SQRT(VEC(1)*VEC(1)+VEC(2)*VEC(2)+VEC(3)*VEC(3))
        IF (DIST.LE.1.0D-08) CYCLE
        DISTI  = ONE/DIST
!
        IF (IDFTBD.EQ.1) THEN
          IF (DIST.GE.R0) THEN !! ATTRACTIVE POTENTIAL
            TMP = RIJ*DISTI
            R6 = TMP**6
            R12 = R6*R6
            IF (NDIMPBC.LE.2) THEN
              VALTMP(1,1) = VALTMP(1,1) + DIJ*(-TWO*R6+R12)
              CYCLE
            END IF
            DIST = DIST*ALPHA
            DIST2 = DIST*DIST
            DISTI2 = ONE/DIST2
            DISTI4 = DISTI2*DISTI2
            DISTI6 = DISTI4*DISTI2
            VALTMP(1,1) = VALTMP(1,1) + DIJ*R12 &
              + BIJ*ALPHA6*(DISTI6+DISTI4+HALF*DISTI2)*EXP(-DIST2)
          ELSE !! REPULSIVE POTENTIAL
            TMP = DIST/RIJ
            R5 = TMP**5
            R10 = R5**2
            VALTMP(1,1) = VALTMP(1,1) + DIJ*(C0-C1*R5+C2*R10)
!
            !! Compensate the R^-6 term
            TMP = RIJ*DISTI
            R6 = TMP**6
            VALTMP(1,1) = VALTMP(1,1) + 2.0D+00*DIJ*R6
!
            DIST = DIST*ALPHA
            DIST2 = DIST*DIST
            DISTI2 = ONE/DIST2
            DISTI4 = DISTI2*DISTI2
            DISTI6 = DISTI4*DISTI2
            VALTMP(1,1) = VALTMP(1,1) &
              + BIJ*ALPHA6*(DISTI6+DISTI4+HALF*DISTI2)*EXP(-DIST2)
          END IF
        ELSE IF (IDFTBD.EQ.2 .OR. IDFTBD.EQ.4) THEN
          DIST7  = DIST**7
          DAMP4 = (ONE-EXP(-3.0D+00*R0AB7I*DIST7))**4
          DIST6I = DISTI**6
          VALTMP(1,1) = VALTMP(1,1) - DAMP4*C6AB*DIST6I
        END IF
!
        IF (NDER.EQ.1) THEN
          IF (IDFTBD.EQ.1) THEN
            DIST = DIST*ALPHAI
            IF (DIST.GE.R0) THEN !! ATTRACTIVE POTENTIAL
              TMP = RIJ*DISTI
              R6 = TMP**6
              R12 = R6*R6
              IF (NDIMPBC.LE.2) THEN
                VAL = 1.2D+01*DIJ*(R6-R12)*DISTI
              ELSE
                VAL = -1.2D+01*DIJ*R12*DISTI
                VAL = VAL &
                    - BIJ*ALPHA6*(6.0D+00*DISTI6+6.0D+00*DISTI4 &
                                 +3.0d+00*DISTI2+1.0D+00) &
                      *EXP(-DIST2)*DISTI
              END IF
            ELSE !! REPULSIVE POTENTIAL
              TMP = DIST/RIJ
              R5 = TMP**5
              R10 = R5**2
              VAL = DIJ*(-C1*R5*5.0D+00+C2*R10*1.0D+01)*DISTI
!
              !! Compensate the R^-6 term
              TMP = RIJ*DISTI
              R6 = TMP**6
              VAL = VAL - 1.2D+01*DIJ*R6*DISTI
              VAL = VAL &
                  - BIJ*ALPHA6*(6.0D+00*DISTI6+6.0D+00*DISTI4 &
                               +3.0d+00*DISTI2+1.0D+00) &
                    *EXP(-DIST2)*DISTI
            END IF
          ELSE IF (IDFTBD.EQ.2 .OR. IDFTBD.EQ.4) THEN
            DIST6I = DISTI**6
            DIST7  = DIST**7
            DIST7I = DIST6I*DISTI
            R7R07 = DIST7*R0AB7I
            VEXP = EXP(-3.0D+00*R7R07)
            DAMP1 = 1.0D+00 - VEXP
            DAMP3 = DAMP1**3
            VAL = 8.4D+01*VEXP*R7R07-6.0D+00*DAMP1
            VAL = -VAL*C6AB*DAMP3*DIST7I
          END IF
!
          !! STRESS TENSOR
          DO KI = 1, 3
            VALTMP(KI+1,1) = VALTMP(KI+1,1) + VEC(KI)*VAL*DISTI
          END DO
          NSEQ = 0
          DO KJ = 1, 3
            DO KI = 1, 3
              NSEQ = NSEQ + 1
                VALTMP(NSEQ+4,1) = VALTMP(NSEQ+4,1) &
                  + VEC(KI)*VEC(KJ)*VAL*DISTI
            END DO
          END DO
        END IF

        END DO
        END DO
        END DO
!
        !! CHECK CONVERGENCE
        DEVMAX = ZERO
        DO I = 1, N
          VAL = ABS((VALTMP(I,1)-VALTMP(I,2)))
          IF (VAL.GE.DEVMAX) DEVMAX = VAL
        END DO
        IF (DEVMAX.LE.THRES .AND. ITRVEC.GE.1) EXIT
        CALL DCOPY(N,VALTMP(1,1),1,VALTMP(1,2),1)
      END DO
!     write (*,*) "converged short-range at", itrvec
!
      !! Finally, add short-range contributions
      CALL DSCAL(N,CI*CJ,VALTMP,1)
      EDISP = EDISP + VALTMP(1,1)
      IF (NDER.EQ.1) THEN
        DO KI = 1, 3
          FMODE(KI,IAG,2) = FMODE(KI,IAG,2) + VALTMP(KI+1,1)
          FMODE(KI,JAG,2) = FMODE(KI,JAG,2) - VALTMP(KI+1,1)
          IF (NBODY.EQ.3) THEN
            FMODE(KI,IAG,3) = FMODE(KI,IAG,3) + VALTMP(KI+1,1)*SCAL
            FMODE(KI,JAG,3) = FMODE(KI,JAG,3) - VALTMP(KI+1,1)*SCAL
          END IF
        END DO
        CALL DAXPY(9,TWO,VALTMP(5,1),1,STRESS,1)
      END IF
!
      IF (NDIMPBC.LE.2 .OR. IDFTBD.NE.1) CYCLE
!
!     ----- RECIPROCAL-SPACE
!
      VEC(1) = VEC1(1) - VEC2(1)
      VEC(2) = VEC1(2) - VEC2(2)
      VEC(3) = VEC1(3) - VEC2(3)
      CALL VCLR(VALTMP,1,2*13)
      DO ITRVEC = 1, 100
        IF (NDIMPBC.GE.1) ITRVECX = ITRVEC
        IF (NDIMPBC.GE.2) ITRVECY = ITRVEC
        IF (NDIMPBC.GE.3) ITRVECZ = ITRVEC
        DO IXVEC = -ITRVECX, ITRVECX
        DO IYVEC = -ITRVECY, ITRVECY
        DO IZVEC = -ITRVECZ, ITRVECZ
        IF (ABS(IXVEC).NE.ITRVEC .AND. ABS(IYVEC).NE.ITRVEC .AND. &
            ABS(IZVEC).NE.ITRVEC) CYCLE
!     
        VECREC(1) = PBCBOX(1,4)*IXVEC &
                  + PBCBOX(2,4)*IYVEC &
                  + PBCBOX(3,4)*IZVEC
        VECREC(2) = PBCBOX(1,5)*IXVEC &
                  + PBCBOX(2,5)*IYVEC &
                  + PBCBOX(3,5)*IZVEC
        VECREC(3) = PBCBOX(1,6)*IXVEC &
                  + PBCBOX(2,6)*IYVEC &
                  + PBCBOX(3,6)*IZVEC
        VECREC(4) = SQRT(VECREC(1)**2 + VECREC(2)**2 + VECREC(3)**2)
        BB = HALF*VECREC(4)*ALPHAI
        BBI= ONE/BB
!        if(IERFTP.NE.0) then
          ERFCBB=GMSERFC(BB)
!        else
!          ERFCBB=ERFC(BB)
!        endif
        EXP_BBBB=EXP(-BB*BB)
!
        VALTMP(1,1) = VALTMP(1,1) &
          + COS(DDOT(3,VEC,1,VECREC,1)) &
            *VECREC(4)*VECREC(4)*VECREC(4) &
            *(PISQ*ERFCBB+(HALF*BBI*BBI-ONE)*BBI*EXP_BBBB)
        IF (NDER.EQ.1) THEN
          !! GRADIENT
          VAL = -SIN(DDOT(3,VEC,1,VECREC,1)) &
                *VECREC(4)*VECREC(4)*VECREC(4) &
                *(PISQ*ERFCBB+(HALF*BBI*BBI-ONE)*BBI*EXP_BBBB)
          DO KI = 1, 3
            VALTMP(KI+1,1) = VALTMP(KI+1,1) + VAL*VECREC(KI)
          END DO
!
          !! STRESS TENSOR
          SCOS = COS(DDOT(3,VEC,1,VECREC,1))
          VAL1 = VECREC(4)*VECREC(4)*VECREC(4) &
                *(PISQ*ERFCBB+(HALF*BBI*BBI-ONE)*BBI*EXP_BBBB)
          VAL2 = 3.0D+00*VECREC(4)*(PISQ*ERFCBB-BBI*EXP_BBBB)
          NSEQ = 0
          DO KJ = 1, 3
            DO KI = 1, 3
              NSEQ = NSEQ + 1
              VAL = VAL2*VECREC(KI)*VECREC(KJ)
              IF (KI.EQ.KJ) VAL = VAL + VAL1
              VALTMP(NSEQ+4,1) = VALTMP(NSEQ+4,1) - SCOS*VAL
            END DO
          END DO
        END IF
!
        END DO
        END DO
        END DO
!
        !! CHECK CONVERGENCE
        DEVMAX = ZERO
        DO I = 1, N
          VAL = ABS((VALTMP(I,1)-VALTMP(I,2)))
          IF (VAL.GE.DEVMAX) DEVMAX = VAL
        END DO
        IF (DEVMAX.LE.THRES) EXIT
        CALL DCOPY(N,VALTMP(1,1),1,VALTMP(1,2),1)
      END DO
!     write (*,*) "converged  long-range at", itrvec
!
      !! Finally, add reciprocal-space contributions
      CALL DSCAL(N,BIJ*PI*PISQ/(1.2D+01*VOL),VALTMP,1)
      CALL DSCAL(N,CI*CJ,VALTMP,1)
      EDISP = EDISP + VALTMP(1,1)
      IF (NDER.EQ.1) THEN
        DO KI = 1, 3
          FMODE(KI,IAG,2) = FMODE(KI,IAG,2) + VALTMP(KI+1,1)
          FMODE(KI,JAG,2) = FMODE(KI,JAG,2) - VALTMP(KI+1,1)
          IF (NBODY.EQ.3) THEN
            FMODE(KI,IAG,3) = FMODE(KI,IAG,3) + VALTMP(KI+1,1)*SCAL
            FMODE(KI,JAG,3) = FMODE(KI,JAG,3) - VALTMP(KI+1,1)*SCAL
          END IF
        END DO
        CALL DAXPY(9,TWO,VALTMP(5,1),1,STRESS,1)
      END IF
!
!     ----- K = 0 FOR 1/R6 TERM
!
      BIJ = BIJ*CI*CJ
      EDISP = EDISP + BIJ*ALPHA3*PI*PISQ/(3.0D+00*VOL)
      IF (NDER.EQ.1) THEN
        DO KI = 1, 3
          STRESS(KI,KI) = STRESS(KI,KI) &
            - BIJ*ALPHA3*PI*PISQ/(3.0D+00*VOL)*TWO
        END DO
      END IF
    END DO
  END DO
!
  ! IF (IDFTBD.EQ.4) CALL RETFM(NEED)
! CALL DFTB_FMOSTR
!
  RETURN
!
END SUBROUTINE DFTB_ESDIM_DCPBC
!
!-----------------------------------------------------------------------
!*MODULE DFTBPB    *DECK DFTB_GDPBC_PREP
!>
!>     @brief Prepare dispersion data.
!>
!>     @details Prepare dispersion data for PBC.
!>
!>     @author Yoshio Nishimoto
!>
SUBROUTINE DFTB_GDPBC_PREP(TMPD,TMPD3,NSPE,HUBBL,DFTB3,DAMPXH0)
!
  IMPLICIT NONE
!
  INTEGER ISP,JSP,NSPE
  DOUBLE PRECISION HALF,ZERO,ONE,THREE,TWO,TAUA,TAUB,TAUMEAN,TAUMEAN2
  DOUBLE PRECISION TAUMEAN3,TAUA2,TAUB2,TAUA4,TAUB4,TAUA6,TAUB6,TAUAB2I
  DOUBLE PRECISION TAUAB3I,hubtol,TWELVE,ONE_15,ONE_24,TAUB3
  DOUBLE PRECISION TAUAB,TAUABI,TAUABI2,TAUABI3,TAUABI4,TAUA3
  PARAMETER (HALF=0.5D+00)
  PARAMETER (ZERO=0.0D+00)
  PARAMETER (ONE=1.0D+00)
  PARAMETER (THREE=3.0D+00)
  PARAMETER (TWO=2.0D+00)
!
  DOUBLE PRECISION TMPD(10,NSPE,NSPE),TMPD3(11,NSPE,NSPE),HUBBL(3,*)
  LOGICAL DFTB3,DAMPXH,DAMPXH0
!
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
! CALL DCOPY(10*NSPE*NSPE,TMPD,1,GAMMAFMOD,1)
!
  IF (DFTB3) THEN
    hubtol = 1.0d-04
    TWELVE = 1.2D+01
    ONE_15 = ONE/1.5D+01
    ONE_24 = ONE/2.4D+01
    DAMPXH = DAMPXH0
    DO ISP = 1, NSPE
      TAUA = 3.2D+00*HUBBL(1,ISP)
      TAUA2 = TAUA**2
      TAUA3 = TAUA*TAUA2
      TAUA4 = TAUA2**2
      TAUA6 = TAUA2*TAUA4
      DO JSP = 1, NSPE
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
!   CALL DCOPY(11*NSPE*NSPE,TMPD3,1,GAMMA3D,1)
  END IF
!
END SUBROUTINE DFTB_GDPBC_PREP
!
!-------------------------------------------------------------------
!
!*MODULE DFTBPB    *DECK DFTB_NUMGRDX
!>
!>     @brief Calculate numerical gradient.
!>
!>     @details Calculate numerical gradient for DFTB.
!>
!>     @author Yoshio Nishimoto
!>
SUBROUTINE DFTB_NUMGRDX(NAT,C,FRC,WRK)
!
  IMPLICIT NONE
  INTEGER NAT,I,IAT,J
!
  DOUBLE PRECISION C(3,*),FRC(3,*),WRK(*)
!
  CALL DSCAL(9,ANG2AU,PBCBOX,1)
  !! FIND RECIPROCAL LATTICE VECTORS AND COMPUTE VOLPBC
  CALL DFTB_RECVEC
!
  CALL DGEMM('N','N',3,NAT,3,             &
             1.0D+00,PBCBOX(1,1),3,FRC,3, &
             0.0D+00,C,3)
!
  IF (MASWRK) THEN
    WRITE (IW,*)
    WRITE (IW,'(X,"REAL-SPACE CELL VECTORS (IN AU)")')
    WRITE (IW,'(14X,"X",19X,"Y",19X,"Z")')
    WRITE (IW,'(X,"A1",3F20.10)') (PBCBOX(I,1),I=1,3)
    WRITE (IW,'(X,"A2",3F20.10)') (PBCBOX(I,2),I=1,3)
    WRITE (IW,'(X,"A3",3F20.10)') (PBCBOX(I,3),I=1,3)
    WRITE (IW,*)
    !!!
    WRITE (IW,'(X,"REAL-SPACE CELL VECTORS (IN ANGSTROM)")')
    WRITE (IW,'(14X,"X",19X,"Y",19X,"Z")')
    WRITE (IW,'(X,"A1",3F20.10)') (PBCBOX(I,1)/ANG2AU,I=1,3)
    WRITE (IW,'(X,"A2",3F20.10)') (PBCBOX(I,2)/ANG2AU,I=1,3)
    WRITE (IW,'(X,"A3",3F20.10)') (PBCBOX(I,3)/ANG2AU,I=1,3)
    WRITE (IW,*)
    !!!
    WRITE (IW,'(X,"RECIPROCAL-SPACE VECTORS (IN AU^-1)")')
    WRITE (IW,'(14X,"X",19X,"Y",19X,"Z")')
    WRITE (IW,'(X,"B1",3F20.10)') (PBCBOX(I,4),I=1,3)
    WRITE (IW,'(X,"B2",3F20.10)') (PBCBOX(I,5),I=1,3)
    WRITE (IW,'(X,"B3",3F20.10)') (PBCBOX(I,6),I=1,3)
    WRITE (IW,*)
    WRITE (IW,'(X,"VOLUME OF THE UNIT CELL = ", F12.5," AU^3")') VOLPBC
    WRITE (IW,'(X,"OPTIMUM EWALD PARAMETER = ", F12.5)') EWALPHA
    WRITE (IW,*)
    WRITE (IW,'(X,"TRANSFORMED COORDINATES")')
    DO IAT = 1, NAT
      WRITE (IW,'(X,I4,3F20.10)') IAT,(C(J,IAT),J=1,3)
    END DO
  END IF
  SCALEXP = 1.0D+00/(4.0D+00*EWALPHA*EWALPHA)
!
  CALL TRPOSQ(PBCBOX(1,1),3)
  CALL TRPOSQ(PBCBOX(1,4),3)
!
END SUBROUTINE DFTB_NUMGRDX
!
!-----------------------------------------------------------------------
!
!>
!>     @brief Calculate CBSP.
!>
!>     @details Calculate CBSP for DFTB.
!>
!>     @author Yoshio Nishimoto
!>
RECURSIVE FUNCTION DFTB_CBSP(NB,XX) RESULT(CBSP)
!
  IMPLICIT NONE
  DOUBLE PRECISION XX,CBSP
  INTEGER NB
!
  CBSP = 0.0D+00
  IF (XX.LT.0.0D+00 .OR. XX.GT.DBLE(NB)) THEN
    CBSP = 0.0D+00
  ELSE IF (NB.EQ.2) THEN
    CBSP = 1.0D+00 - ABS(XX-1.0D+00)
  ELSE 
    CBSP = XX*DFTB_CBSP(NB-1,XX)/DBLE(NB-1) &
         + DBLE(NB-XX)/DBLE(NB-1)*DFTB_CBSP(NB-1,XX-1.0D+00)
  END IF
!
END FUNCTION DFTB_CBSP
!
!-----------------------------------------------------------------------
!
!*MODULE DFTBPB    *DECK DFTB_PMEQ
!>
!>     @brief Calculate PMEQ.
!>
!>     @details Calculate PMEQ for DFTB.
!>
!>     @author Yoshio Nishimoto
!>
SUBROUTINE DFTB_PMEQ(MODE,NK,NB,BOXREC,VEC,CHAMUL,QQ,VAL,GRAD)
!
  IMPLICIT NONE
!
  INTEGER MODE,NB
  DOUBLE PRECISION CHAMUL,VAL,TPINV,U01,U02,U03,U1,VM1,VM11
  INTEGER K1S,K1L,K2S,K2L,K3S,K3L,K1,NK,K1T,K2,K2T,K3,K3T
  DOUBLE PRECISION VEC(3),BOXREC(3,3),QQ(NK,NK,NK),GRAD(3)
  DOUBLE PRECISION U2,VM2,VM21,U3,VM3,VM31,SCAL
!
  TPINV = 1.0D+00/(2.0D+00*PI)
  U01 = NK*(BOXREC(1,1)*VEC(1)+BOXREC(2,1)*VEC(2) &
           +BOXREC(3,1)*VEC(3))*TPINV
! U01 = DDOT(3,BOXREC(1,1),1,VEC,1)*TPINV*NK
  IF (U01.LT.0) U01 = U01 + DBLE(NK)
  U02 = NK*(BOXREC(1,2)*VEC(1)+BOXREC(2,2)*VEC(2) &
           +BOXREC(3,2)*VEC(3))*TPINV
! U02 = DDOT(3,BOXREC(1,2),1,VEC,1)*TPINV*NK
  IF (U02.LT.0) U02 = U02 + DBLE(NK)
  U03 = NK*(BOXREC(1,3)*VEC(1)+BOXREC(2,3)*VEC(2) &
           +BOXREC(3,3)*VEC(3))*TPINV
! U03 = DDOT(3,BOXREC(1,3),1,VEC,1)*TPINV*NK
  IF (U03.LT.0) U03 = U03 + DBLE(NK)
  K1S = U01 - NB
  K1L = K1S + NB + 1
  K2S = U02 - NB
  K2L = K2S + NB + 1
  K3S = U03 - NB
  K3L = K3S + NB + 1
!
  IF (MODE.EQ.2) THEN
    VAL = 0.0D+00
  ELSE IF (MODE.EQ.3) THEN
    GRAD(1) = 0.0D+00
    GRAD(2) = 0.0D+00
    GRAD(3) = 0.0D+00
  END IF
!
  DO K1 = K1S, K1L
    U1 = DBLE(U01-K1)
    VM1 = DFTB_CBSP(NB,U1)
    IF (MODE.EQ.3) VM11 = DFTB_CBSP(NB-1,U1)-DFTB_CBSP(NB-1,U1-1.0D+00)
    K1T = K1+1
    IF (K1T.LE.0) THEN
      K1T = K1T + NK*(ABS(K1T)/NK+1)
    ELSE IF (K1T.NE.NK) THEN
      K1T = MOD(K1T,NK)
    END IF
    DO K2 = K2S, K2L
      U2 = DBLE(U02-K2)
      VM2 = DFTB_CBSP(NB,U2)
      IF (MODE.EQ.3) VM21 = DFTB_CBSP(NB-1,U2)-DFTB_CBSP(NB-1,U2-1.0D+00)
      K2T = K2+1
      IF (K2T.LE.0) THEN
        K2T = K2T + NK*(ABS(K2T)/NK+1)
      ELSE IF (K2T.NE.NK) THEN
        K2T = MOD(K2T,NK)
      END IF
      DO K3 = K3S, K3L
        U3 = DBLE(U03-K3)
        VM3 = DFTB_CBSP(NB,U3)
        IF (MODE.EQ.3) VM31 = DFTB_CBSP(NB-1,U3)-DFTB_CBSP(NB-1,U3-1.0D+00)
        K3T = K3+1
        IF (K3T.LE.0) THEN
          K3T = K3T + NK*(ABS(K3T)/NK+1)
        ELSE IF (K3T.NE.NK) THEN
          K3T = MOD(K3T,NK)
        END IF
!
        IF (MODE.EQ.1) THEN
          QQ(K1T,K2T,K3T) = QQ(K1T,K2T,K3T) + VM1*VM2*VM3*CHAMUL
        ELSE IF (MODE.EQ.2) THEN
          VAL = VAL + VM1*VM2*VM3*QQ(K1T,K2T,K3T)
        ELSE IF (MODE.EQ.3) THEN
          SCAL = QQ(K1T,K2T,K3T)*TPINV*NK
          GRAD(1) = GRAD(1) &
           + SCAL*(BOXREC(1,1)+BOXREC(1,2)+BOXREC(1,3))*VM11*VM2*VM3
          GRAD(2) = GRAD(2) &
           + SCAL*(BOXREC(2,1)+BOXREC(2,2)+BOXREC(2,3))*VM1*VM21*VM3
          GRAD(3) = GRAD(3) &
           + SCAL*(BOXREC(3,1)+BOXREC(3,2)+BOXREC(3,3))*VM1*VM2*VM31
        END IF
!
      END DO
    END DO
  END DO
!
END SUBROUTINE DFTB_PMEQ
!
!-----------------------------------------------------------------------
!
!*MODULE DFTBPB    *DECK DFTB_PMEBC
!>
!>     @brief Calculate PMEBC.
!>
!>     @details Calculate PMEBC for DFTB.
!>
!>     @author Yoshio Nishimoto
!>
SUBROUTINE DFTB_PMEBC(NDER,FQ)
!
  IMPLICIT NONE
!
  INTEGER NDER,NSEQ,IZVEC,IZVEC3,IYVEC,IYVEC2,IYVEC3,IZVEC2,IXVEC,IXVEC2
  INTEGER IXVEC3,L,K
  DOUBLE PRECISION FQ(2,NKFFT/2+1,NKFFT,NKFFT),VECREC(4)
  DOUBLE PRECISION BBZ,BBY,BBX,VECRECI,EXPFACT,SCAL,VAL0,VAL
!
  COMPLEX(KIND(0D0)) TMPZ
!
  !! Evaluate theta = B*C
  NSEQ = 0
  DO IZVEC = 1, NKFFT
    IZVEC2 = IZVEC-1
    BBZ = BSQ(IZVEC2) ! X(LBSQ+IZVEC2)
    IZVEC3 = IZVEC2
    IF (IZVEC2.GT.NKFFT/2) IZVEC3 = IZVEC3 - NKFFT
    DO IYVEC = 1, NKFFT
      IYVEC2 = IYVEC-1
      BBY = BSQ(IYVEC2) ! X(LBSQ+IYVEC2)
      IYVEC3 = IYVEC2
      IF (IYVEC2.GT.NKFFT/2) IYVEC3 = IYVEC3 - NKFFT
      DO IXVEC = 1, NKFFT/2+1
        NSEQ = NSEQ + 1
        IF (NSEQ.EQ.1) CYCLE !! Zero for F(Q)(0,0,0)
        IXVEC2 = IXVEC-1
        BBX = BSQ(IXVEC2) ! X(LBSQ+IXVEC2)
        IXVEC3 = IXVEC2
        VECREC(1) = BOXREC(1,1)*IXVEC3 &
                  + BOXREC(2,1)*IYVEC3 &
                  + BOXREC(3,1)*IZVEC3
        VECREC(2) = BOXREC(1,2)*IXVEC3 &
                  + BOXREC(2,2)*IYVEC3 &
                  + BOXREC(3,2)*IZVEC3
        VECREC(3) = BOXREC(1,3)*IXVEC3 &
                  + BOXREC(2,3)*IYVEC3 &
                  + BOXREC(3,3)*IZVEC3
        VECREC(4) = VECREC(1)**2 + VECREC(2)**2 + VECREC(3)**2
        VECRECI   = 1.0D+00/VECREC(4)
        !! EXPFACT = EXP(-pi^2*m^2/beta^2)/m^2*Bx(m)^2*By(m)^2*Bz(m)^2
        !! Note that the reciprocal lattice vector is defined
        !! differently. In the reference, a*b = 1, but a*b=2*pi
        !! in the implementation.
        EXPFACT   = EXP(-SCALEXP*VECREC(4))*VECRECI*BBX*BBY*BBZ
        IF (NDER.EQ.1) THEN
          !! Compute the stress tensor
          !! Eq. (2.7) in JCP 1995, 103, 8577-8593.
          TMPZ = DCMPLX(FQ(1,IXVEC,IYVEC,IZVEC),FQ(2,IXVEC,IYVEC,IZVEC))
!         TMPZ = DCMPLX(X(LFQ+2*(NSEQ-1)),X(LFQ+2*(NSEQ-1)+1))
          !! S(m) = conjug(S(-m)),
          !! so S(m)*S(-m)) = F(Q)*conjg(F(Q))*Bx^2*By^2*Bz^2
          !! Additional scaling by two is unclear.
          SCAL = 4.0D+00*PI*VOLINV*TMPZ*DCONJG(TMPZ)*EXPFACT
          IF (IXVEC.NE.1) SCAL=SCAL+SCAL
          VAL0 = -2.0D+00*(1.0D+00+SCALEXP*VECREC(4))*VECRECI
          DO L = 1, 3
            DO K = 1, 3
              VAL = VAL0*VECREC(K)*VECREC(L)
              IF (K.EQ.L) VAL = VAL + 1.0D+00
              STRESS(K,L) = STRESS(K,L) - SCAL*VAL
            END DO
          END DO
        END IF
        FQ(1,IXVEC,IYVEC,IZVEC) = FQ(1,IXVEC,IYVEC,IZVEC)*EXPFACT
        FQ(2,IXVEC,IYVEC,IZVEC) = FQ(2,IXVEC,IYVEC,IZVEC)*EXPFACT
      END DO
    END DO
  END DO
  FQ(1,1,1,1) = 0.0D+00
  FQ(2,1,1,1) = 0.0D+00
!
  RETURN
END SUBROUTINE DFTB_PMEBC
!
!-----------------------------------------------------------------------
!
!*MODULE DFTBPB    *DECK DFTB_ESP_PME
!>
!>     @brief Calculate ESP with PME.
!>
!>     @details Calculate ESP with PME for DFTB.
!>
!>     @author Yoshio Nishimoto
!>
SUBROUTINE DFTB_ESP_PME(IFG,JFG,KFG,NAT,POPMAT,SHIFT)
!
  IMPLICIT NONE
!
!PME  include '/home/nisimoto/lib/fftw/fftw-3.3.7/api/fftw3.f'
!
  COMMON /FMCOM / X(1)
!
  DOUBLE PRECISION POPMAT(MAXNAT,NFG),SHIFT(*)
  DOUBLE PRECISION VEC(3),X,VAL
  INTEGER LOADFM,LAST,LQQ,LFQ,LPLAN,NEED,NKFFT3,IFG,JFG,KFG,NAT,NFG0
  INTEGER LFG0,NATI,INDI,IAT,IAG,LFG
!
  CALL VALFM(LOADFM)
  LAST  = LOADFM
  LQQ   = LOADFM + 1
  LFQ   = LQQ    + NKFFT*NKFFT*NKFFT
  LPLAN = LFQ    + 2*(NKFFT/2+1)*NKFFT*NKFFT
  LAST  = LPLAN  + 8
  NEED  = LAST   - LOADFM - 1
  CALL GETFM(NEED)
  NKFFT3 = NKFFT*NKFFT*NKFFT
!
  CALL VCLR(X(LQQ),1,NKFFT3)
  IF (IFG.EQ.0) NFG0 = NFG
  IF (IFG.NE.0) NFG0 = 1
  IF (JFG.NE.0) NFG0 = 2
  IF (KFG.NE.0) NFG0 = 3
!
  DO LFG0 = 1, NFG0
    IF (IFG.EQ.0) THEN
      LFG = LFG0
    ELSE
      IF (LFG0.EQ.1) LFG = IFG
      IF (LFG0.EQ.2) LFG = JFG
      IF (LFG0.EQ.3) LFG = KFG
    END IF
    NATI = NATFRG(LFG)
    INDI = INDFRG(LFG)
    DO IAT = 1, NATI
      IAG = IATFRG(INDI+IAT-1)
      CALL DFTB_PMEQ(1,NKFFT,NCBSP,BOXREC,FMOC(1:3,IAG), &
                     POPMAT(IAT,LFG)-ZREFFMO(INDI+IAT-1), &
                     X(LQQ),VAL,VEC)
    END DO
  END DO
!
  !! Q -> F(Q)
!PME  CALL DFFTW_PLAN_DFT_R2C_3D(X(LPLAN),NKFFT,NKFFT,NKFFT,X(LQQ), &
!PME                             X(LFQ),FFTW_ESTIMATE)
!PME  CALL DFFTW_EXECUTE_DFT_R2C(X(LPLAN),X(LQQ),X(LFQ))
!PME  CALL DFFTW_DESTROY_PLAN(X(LPLAN))
!
  !! Evaluate theta = B*C
  CALL DFTB_PMEBC(0,X(LFQ))
!
  !! F^-1(theta*Q)
!PME  CALL DFFTW_PLAN_DFT_C2R_3D(X(LPLAN),NKFFT,NKFFT,NKFFT,X(LFQ), &
!PME                             X(LQQ),FFTW_ESTIMATE)
!PME  CALL DFFTW_EXECUTE_DFT_C2R(X(LPLAN),X(LFQ),X(LQQ))
!PME  CALL DFFTW_DESTROY_PLAN(X(LPLAN))
!
  IF (IFG.EQ.0) THEN
    DO LFG = 1, NFG
      NATI = NATFRG(LFG)
      INDI = INDFRG(LFG)
      DO IAT = 1, NATI
        IAG = IATFRG(INDI+IAT-1)
        CALL DFTB_PMEQ(2,NKFFT,NCBSP,BOXREC,FMOC(1:3,IAG),1.0D+00,X(LQQ),VAL,VEC)
        SHIFT(INDI+IAT-1) = SHIFT(INDI+IAT-1) + VAL*4.0D+00*PI*VOLINV
      END DO
    END DO
  ELSE
    DO IAT = 1, NAT
      IAG = IAGLOB(IAT)
      CALL DFTB_PMEQ(2,NKFFT,NCBSP,BOXREC,FMOC(1:3,IAG),1.0D+00,X(LQQ),VAL,VEC)
      SHIFT(IAT) = SHIFT(IAT) + VAL*4.0D+00*PI*VOLINV
    END DO
  END IF
!
  CALL RETFM(NEED)
! write (*,*) "*************PME*****************"
!
  RETURN
!
END SUBROUTINE DFTB_ESP_PME
!
!-----------------------------------------------------------------------
!
!*MODULE DFTBPB    *DECK dftb_getgam
!>
!>     @brief Get gamma.
!>
!>     @details Get gamma for DFTB.
!>
!>     @author Yoshio Nishimoto
!>
subroutine dftb_getgam(nat,gam2fmo,gam3fmo,gamma2,gamma3,nn,dftb3)
!
  IMPLICIT NONE
!
  INTEGER LGAMMA,LZREF,LISPE,LSHIFTG,LSHIFTCT, &
                  LGAMMAD,MODESD,MODGAMMA,LDCTMP,LGAMMA3,LGAMMA3D, &
                  LNBOND,LCTMUL,LESPDFTB,LESPASC,NDFTBRST,LDFTBRST
  INTEGER nseq,iat,iag,jat,jag,nseq2,iatfg,iatloc,indi,jatfg,jatloc,indj
  INTEGER nat,nn
!
  COMMON /FMODTB/ LGAMMA,LZREF,LISPE,LSHIFTG,LSHIFTCT, &
                  LGAMMAD,MODESD,MODGAMMA,LDCTMP,LGAMMA3,LGAMMA3D, &
                  LNBOND,LCTMUL,LESPDFTB,LESPASC,NDFTBRST,LDFTBRST
!
  double precision gam2fmo(nn,nn),gam3fmo(nn,nn,*),gamma2(*), &
            gamma3(nat,nat)
  logical   dftb3
!
  if (iand(modgamma,3).eq.1) then
    nseq = 0
    do iat = 1, nat
      iag = iaglob(iat)
      do jat = 1, iat
        jag = iaglob(jat)
        call dftb_cnvsq(iag,jag,nseq2)
        nseq = nseq + 1
        gamma2(nseq) = gam2fmo(nseq2,1)
        if (dftb3) then
          gamma3(iat,jat) = gam3fmo(iag,jag,1)
          gamma3(jat,iat) = gam3fmo(jag,iag,1)
        end if
      end do
    end do
  else if (iand(modgamma,3).eq.2) then
    nseq = 0
    do iat = 1, nat
      iag = iaglob(iat)
      iatfg = indat(iag)
      iatloc = iand(ialoc(iag),65535)
      indi = indfrg(iatfg)
      do jat = 1, iat
        jag = iaglob(jat)
        jatfg = indat(jag)
        jatloc = iand(ialoc(jag),65535)
        indj = indfrg(jatfg)
        nseq = nseq + 1
        gamma2(nseq) = gam2fmo(indi+iatloc-1,indj+jatloc-1)
        if (dftb3) then
          gamma3(iat,jat) = gam3fmo(indi+iatloc-1,indj+jatloc-1,1)
          gamma3(jat,iat) = gam3fmo(indi+iatloc-1,indj+jatloc-1,2)
        end if
      end do
    end do
  else
    if (maswrk) write (iw,'(" this should not happen in dftb_getgam")')
    call abrt
  end if
!
  return
!
end subroutine dftb_getgam
!
!-----------------------------------------------------------------------
!
!*MODULE DFTBPB    *DECK DFTB_LATREC
!>
!>     @brief Print lattice.
!>
!>     @details Print lattice in PBC.
!>
!>     @author Yoshio Nishimoto
!>
  SUBROUTINE DFTB_LATREC(OUT)
  USE MX_LIMITS, ONLY: MXATM
!
  IMPLICIT NONE
  LOGICAL OUT
  INTEGER I,NAT,ICH,MUL,NUM,NQMT,NE,NA,NB,IAN
  DOUBLE PRECISION ZAN,C
  COMMON /INFOA / NAT,ICH,MUL,NUM,NQMT,NE,NA,NB, &
                  ZAN(MXATM),C(3,MXATM),IAN(MXATM)
!
  IF (.NOT.PERIOD.OR.LATOPT.EQ.0) RETURN
!
  CALL TRPOSQ(BOXREAL,3)
  CALL DFTB_RECVEC
  IF (OUT) THEN
    WRITE (IW,'(X,"REAL-SPACE CELL PARAMETERS (IN ANGSTROM)")')
    WRITE (IW,'(14X,"X",19X,"Y",19X,"Z")')
    WRITE (IW,'(X,"A1",3F20.10)') (BOXREAL(I,1)/ANG2AU,I=1,3)
    IF (NDIMPBC.GE.2) &
       WRITE (IW,'(X,"A2",3F20.10)') (BOXREAL(I,2)/ANG2AU,I=1,3)
    IF (NDIMPBC.GE.3) &
       WRITE (IW,'(X,"A3",3F20.10)') (BOXREAL(I,3)/ANG2AU,I=1,3)
  END IF
  CALL TRPOSQ(BOXREAL,3)
  CALL TRPOSQ(BOXREC,3)
  IF (NFG.EQ.0) THEN
    CALL DFTB_EWALPHA(NAT,NBODY)
  ELSE
    CALL DFTB_EWALPHA(NATFMO,NBODY)
  END IF
!
  RETURN
!
  end subroutine dftb_latrec
!
!

END MODULE DFTBPB_MOD
