!> @brief Harmonic restraints (HRs)
!> @details Harmonic restraint are defined as \f$ E_{HR} = k\cdot(x-x_0)^2 \f$, where
!>  \f$ E_{HR} \f$ - HR energy, \f$ k \f$ - force constant, \f$ x, x_0 \f$ - current and
!>  equilibrium value of the HR respectively
!> @note HRs are defined as in CHARMM force field without 1/2 prefactor
!> @author Vladimir Mironov
!> @date Jul 2015
! This code needs refactoring
MODULE harmonic_restraints

 USE prec, ONLY: fp
 IMPLICIT NONE

 INTEGER, PARAMETER :: MXRESTR = 100 !< Maximum number of harmonic restraints of each kind

 REAL(KIND=fp), PARAMETER :: BOHR = 1.889726125_fp       !< 1 bohr in angstroms
 REAL(KIND=fp), PARAMETER :: ANGSTROM  = 1.0_fp/BOHR     !< 1 angstrom in bohrs
 REAL(KIND=fp), PARAMETER :: HARTREE = 0.001593601446_fp !< 1 a.u. in kcal/mol
 REAL(KIND=fp), PARAMETER :: PI = 3.1415926535897932_fp  !< \f$ \pi \f$ constant
 REAL(KIND=fp), PARAMETER :: TWOPI = 2*PI                !< \f$ 2\pi \f$
 REAL(KIND=fp), PARAMETER :: RADIAN = PI/180.0_fp        !< 1 radian in degrees
 REAL(KIND=fp), PARAMETER :: HB2 = HARTREE/BOHR/BOHR     !< 1 \f$ hartree/bohr^2 \f$ in \f$ kcal/(mol\cdot A^2) \f$

 REAL(KIND=fp), PARAMETER :: GMAX = 10.0 !< gradient limit for bonds in cartesian coordinates
 REAL(KIND=fp), PARAMETER :: DEMAX = 1.0 !< gradient limit for angles and dihedrals in internal coordinates

 INTEGER, SAVE :: nBndR = 0 !< Number of bond HRs
 INTEGER, SAVE :: nAngR = 0 !< Number of angle HRs
 INTEGER, SAVE :: nDihR = 0 !< Number of dihedral HRs

 REAL(KIND=fp), DIMENSION(MXRESTR), SAVE :: bndR0 !< Equilibrium bond values of HRs
 REAL(KIND=fp), DIMENSION(MXRESTR), SAVE :: angR0 !< Equilibrium angle values of HRs
 REAL(KIND=fp), DIMENSION(MXRESTR), SAVE :: dihR0 !< Equilibrium dihedral values of HRs

 REAL(KIND=fp), DIMENSION(MXRESTR), SAVE :: bndRK = 0 !< Bond HRs force constants
 REAL(KIND=fp), DIMENSION(MXRESTR), SAVE :: angRK = 0 !< Angle HRs force constants
 REAL(KIND=fp), DIMENSION(MXRESTR), SAVE :: dihRK = 0 !< Dihedral HRs force constants


 INTEGER, DIMENSION(2,MXRESTR), SAVE :: IDBndR = 0 !< Atomic IDs for bond HRs
 INTEGER, DIMENSION(3,MXRESTR), SAVE :: IDAngR = 0 !< Atomic IDs for angle HRs
 INTEGER, DIMENSION(4,MXRESTR), SAVE :: IDDihR = 0 !< Atomic IDs for dihedral HRs

! Data for printout
 REAL(KIND=fp), DIMENSION(MXRESTR), SAVE :: bndE = 0 !< Bond HR energies
 REAL(KIND=fp), DIMENSION(MXRESTR), SAVE :: angE = 0 !< Angle HR energies
 REAL(KIND=fp), DIMENSION(MXRESTR), SAVE :: dihE = 0 !< Dihedral HR energies
! Differences from equilibria
 REAL(KIND=fp), DIMENSION(MXRESTR), SAVE :: bndDelta = 0 !< Bond HR deltas
 REAL(KIND=fp), DIMENSION(MXRESTR), SAVE :: angDelta = 0 !< Angle HR deltas
 REAL(KIND=fp), DIMENSION(MXRESTR), SAVE :: dihDelta = 0 !< Dihedral HR deltas

! Current values of restrained coordinates
 REAL(KIND=fp), DIMENSION(MXRESTR), SAVE :: bndCur = 0 !< Bond HR current values
 REAL(KIND=fp), DIMENSION(MXRESTR), SAVE :: angCur = 0 !< Angle HR current values
 REAL(KIND=fp), DIMENSION(MXRESTR), SAVE :: dihCur = 0 !< Dihedral HR current values

 PRIVATE
 PUBLIC  hr_init
 PUBLIC  hr_echo
 PUBLIC  hr_status
 PUBLIC  hr_bnd
 PUBLIC  hr_ang
 PUBLIC  hr_dih

CONTAINS

!-------------------------------------------------------------------------------
!> @brief Init HR arrays
!> @param[in]  natm         number of atoms
!> @param[in]  ihmcon       array of HR types and atomic indices, format is same to IZMAT
!> @param[in]  fhmcon       array of HR force constants (\f$ k \f$ in \f$ E = k(x-x_0)^2 \f$)
!> @param[in]  shmcon       array of HR equilibrium values (\f$ x_0 \f$ in \f$ E = k(x-x_0)^2 \f$)
 SUBROUTINE hr_init(natm, ihmcon, fhmcon, shmcon)

     INTEGER, INTENT(IN) :: natm, ihmcon(:)
     REAL(KIND=fp), INTENT(IN) :: fhmcon(:), shmcon(:)

     INTEGER :: i,n

     REAL(KIND=fp) :: phi, sine, cosine

     CHARACTER(LEN=*), PARAMETER :: &
         fErrRestType = '(1X,"ERROR! RESTRAINT TYPE COULD BE &
                             &ONLY 1, 2 AND 3,",/&
                         &1X,"BUT YOU SPECIFIED IHMCON(",I5,")=",I5)', &

         fErrAtomID   = '(1X,"ERROR! SOME OF ATOMIC IDS IN IHMCON ARRAY ARE &
                            &OUT OF ALLOWED RANGE!",/&
                         &1X,"PLEASE, REVISE YOU $STATPT INPUT SECTION")', &

         fErrForceVal = '(1X,"ERROR! INCORRECT VALUE OF FORCE &
                            &CONSTANT SPECIFIED",/&
                         &1X,"PLEASE, REVISE YOU $STATPT INPUT SECTION")'

     i = 1
     n = 1
     DO WHILE (ihmcon(i)/=0)
         SELECT CASE (ihmcon(i))
         CASE (1)
             nBndR = nBndR + 1
             IDBndR(1:2,nBndR) = ihmcon(i+1:i+2)
             bndR0(nBndR) = BOHR*shmcon(n)
             bndRK(nBndR) = HB2*fhmcon(n)
             i = i + 3
             n = n + 1
         CASE (2)
             nAngR = nAngR + 1
             IDAngR(1:3,nAngR) = ihmcon(i+1:i+3)
             phi = RADIAN*shmcon(n)
             cosine = cos(phi)
             angR0(nAngR) = acos(cosine)
             angRK(nAngR) = HARTREE*fhmcon(n)
             i = i + 4
             n = n + 1
         CASE (3)
             nDihR = nDihR + 1
             IDDihR(1:4,nDihR) = ihmcon(i+1:i+4)
             phi = RADIAN*shmcon(n)
             sine = sin(phi)
             cosine = cos(phi)
             dihR0(nDihR) = atan2(sine, cosine)
             dihRK(nDihR) = HARTREE*fhmcon(n)
             i = i + 5
             n = n + 1
         CASE DEFAULT
             WRITE(*,fErrRestType) i, ihmcon(i)
             CALL ABRT
         END SELECT
     END DO

     ! Check atomic IDs for correctness
     IF (any(IDBndR(:,1:nBndR)>natm).OR.&
         any(IDAngR(:,1:nAngR)>natm).OR.&
         any(IDDihR(:,1:nDihR)>natm).OR.&
         any(IDBndR(:,1:nBndR)<=0).OR.&
         any(IDAngR(:,1:nAngR)<=0).OR.&
         any(IDDihR(:,1:nDihR)<=0)) THEN

         WRITE(*,fErrAtomID)
         CALL ABRT
     END IF

     ! Check that force constants are positive
     IF (ANY(bndRK(1:nBndR)<0).OR.&
         ANY(angRK(1:nAngR)<0).OR.&
         ANY(dihRK(1:nDihR)<0)) THEN

         WRITE(*,fErrForceVal)
         CALL ABRT
     END IF

 END SUBROUTINE hr_init

!-------------------------------------------------------------------------------

!> @brief Print HR parameters
!> @param[in] iw    output file descriptor
 SUBROUTINE hr_echo(iw)

     INTEGER, INTENT(IN) :: iw

     INTEGER :: i

     CHARACTER(LEN=*), PARAMETER :: &
         fPrtHdrBnd = '(//10X,"HARMONIC RESTRAINTS: BONDS")', &
         fPrtDatBnd = '(1X,"I=",I5," J=",I5,16X,&
                            &" R0=",F6.4," K=",F9.3)', &

         fPrtHdrAng = '(//10X,"HARMONIC RESTRAINTS: ANGLES")', &
         fPrtDatAng = '(1X,"I=",I5," J=",I5," K=",I5,8X,&
                            &" A0=",F6.3," K=",F9.3)', &

         fPrtHdrDih = '(//10X,"HARMONIC RESTRAINTS: DIHEDRALS")', &
         fPrtDatDih = '(1X,"I=",I5," J=",I5," K=",I5," L=",I5,&
                            &" D0=",F6.3," K=",F9.3)'

     IF (nBndR>0) THEN
         WRITE(iw,fPrtHdrBnd)
         DO i = 1, nBndR
             WRITE(iw,fPrtDatBnd) IDBndR(:,i), bndR0(i), bndRK(i)
         END DO
     END IF

     IF (nAngR>0) THEN
         WRITE(iw,fPrtHdrAng)
         DO i = 1, nAngR
             WRITE(iw,fPrtDatAng) IDAngR(:,i), angR0(i), angRK(i)
         END DO
     END IF

     IF (nDihR>0) THEN
         WRITE(iw,fPrtHdrDih)
         DO i = 1, nDihR
             WRITE(iw,fPrtDatDih) IDDihR(:,i), dihR0(i), dihRK(i)
         END DO
     END IF

     WRITE(iw,*)


 END SUBROUTINE hr_echo

!-------------------------------------------------------------------------------

!> @brief Print current status of HRs (value, energy, etc.)
!> @param[in] iw    output file descriptor
 SUBROUTINE hr_status(iw)

     INTEGER, INTENT(IN) :: iw

     INTEGER :: i

     REAL(KIND=fp) :: ener

     CHARACTER(LEN=*), PARAMETER :: &
         fPrtHdr    = '(//16X,"-----------------",&
                       &/16X,"RESTRAINTS STATUS",&
                       &/16X,"-----------------",&
                       &/1X,"TYPE",1X,"| ATOMIC INDICES",10X,&
                       &1X,"| VALUE(A,DEG)",&
                       &1X,"| DELTA(A,DEG)",&
                       &1X,"| ENERGY(A.U.)")', &

         fPrtDatBnd = '(1X,"BND",4X,2(I5,1X),14X,2(F13.5,2X),F12.8)', &

         fPrtDatAng = '(1X,"ANG",4X,3(I5,1X),8X,2(F13.5,2X),F12.8)', &

         fPrtDatDih = '(1X,"DIH",4X,4(I5,1X),2X,2(F13.5,2X),F12.8)', &

         fPrtTotalE = '(32X,"TOTAL RESTRAINT ENERGY (A.U.) = ",F12.8,/)'

     WRITE(iw,fPrtHdr)

     DO i = 1, nBndR
         WRITE(iw,fPrtDatBnd) IDBndR(:,i), bndCur(i), bndDelta(i), bndE(i)
     END DO

     DO i = 1, nAngR
         WRITE(iw,fPrtDatAng) IDAngR(:,i), angCur(i), angDelta(i), angE(i)
     END DO

     DO i = 1, nDihR
         WRITE(iw,fPrtDatDih) IDDihR(:,i), dihCur(i), dihDelta(i), dihE(i)
     END DO

     ener = sum(bndE(1:nBndR)) + &
            sum(angE(1:nAngR)) + &
            sum(dihE(1:nDihR))

     WRITE(iw,fPrtTotalE) ener

 END SUBROUTINE hr_status

!-------------------------------------------------------------------------------

!> @brief Compute contribution of bond restraints to energy and gradient
!> @param[inout] ener       energy of the system (a.u.)
!> @param[inout] grad       gradient of the system (a.u.)
!> @param[in]    coor       coordinates of all atoms
!> @param[in]    limit      if `.TRUE.` max. gradient contribution from HR will be limited by `GMAX`
 SUBROUTINE hr_bnd(ener, grad, coor, limit)

     REAL(KIND=fp), INTENT(INOUT) :: ener, grad(:)
     REAL(KIND=fp), INTENT(IN) :: coor(:,:)
     LOGICAL, INTENT(IN) :: limit

     INTEGER :: i, ia, ib

     REAL(KIND=fp) :: &
         a(3), b(3), ab(3), de(3), &
         Fk, r0, delta, rab

     DO i = 1, nBndR

         ia = IDBndR(1,i)
         ib = IDBndR(2,i)
         Fk  = bndRK(i)
         r0 = bndR0(i)
         a(1:3) = coor(1:3,ia)
         b(1:3) = coor(1:3,ib)

         rab = hr_veclen(a(1:3)-b(1:3))
         delta = rab - r0
         ab(1:3) = (a(1:3) - b(1:3))/rab
         de(1:3) = 2.0_fp * delta * Fk * ab(1:3)
         IF (limit) CALL hr_vscale(de,gmax)

!        Data for detailed printout
         bndCur(i)   = rab*ANGSTROM
         bndDelta(i) = delta*ANGSTROM
         bndE(i)     = Fk * delta * delta

         ener = ener + Fk * delta * delta

         grad((ia-1)*3+1:(ia-1)*3+3) = grad((ia-1)*3+1:(ia-1)*3+3) + de(1:3)
         grad((ib-1)*3+1:(ib-1)*3+3) = grad((ib-1)*3+1:(ib-1)*3+3) - de(1:3)

     END DO

 END SUBROUTINE hr_bnd

!-------------------------------------------------------------------------------

!> @brief Compute contribution of angle restraints to energy and gradient
!> @param[inout] ener       energy of the system (a.u.)
!> @param[inout] grad       gradient of the system (a.u.)
!> @param[in]    coor       coordinates of all atoms
!> @param[in]    limit      if `.TRUE.` max. gradient contribution from HR will be limited by `DEMAX`
 SUBROUTINE hr_ang(ener, grad, coor, limit)

     REAL(KIND=fp), INTENT(INOUT) :: ener, grad(:)
     REAL(KIND=fp), INTENT(IN) :: coor(:,:)
     LOGICAL, INTENT(IN) :: limit

     INTEGER :: &
         i, ia, ib, ic

     REAL(KIND=fp), DIMENSION(3):: &
         a, b, c, &
         ab, cb, p, &
         deda, dedb, dedc

     REAL(KIND=fp) :: &
         Fk, r0, &
         rab, rcb, rp, &
         cosine, ang, &
         delta, de, af, cf

     DO i = 1, nAngR

         ia = IDAngR(1,i)
         ib = IDAngR(2,i)
         ic = IDAngR(3,i)
         Fk = angRK(i)
         r0 = angR0(i)

         a(1:3) = coor(1:3,ia)
         b(1:3) = coor(1:3,ib)
         c(1:3) = coor(1:3,ic)

         ab(:) = a(:) - b(:)
         cb(:) = c(:) - b(:)

         rab = hr_veclen(ab)
         rcb = hr_veclen(cb)

         p = hr_crossprod(cb,ab)
         rp = hr_veclen(p)

         cosine = sum(ab(:)*cb(:))/(rab*rcb)
         ang    = acos(cosine)
         delta  = ang - r0

!        Data for detailed printout
         angCur(i)   = ang/PI*180
         angDelta(i) = delta/PI*180
         angE(i)     = Fk * delta * delta


         ener = ener + Fk * delta * delta
         de = 2.0_fp * Fk * delta

         IF (limit) THEN
             IF (abs(de)>DEMAX) THEN
                 de = sign(DEMAX,de)
             END IF
         END IF

         af = -de/(rab*rab*rp)
         cf =  de/(rcb*rcb*rp)
         deda = af * hr_crossprod(ab,p)
         dedc = cf * hr_crossprod(cb,p)
         dedb = -(deda+dedc)

         grad((ia-1)*3+1:(ia-1)*3+3) = grad((ia-1)*3+1:(ia-1)*3+3) + &
                                         deda(1:3)

         grad((ib-1)*3+1:(ib-1)*3+3) = grad((ib-1)*3+1:(ib-1)*3+3) + &
                                         dedb(1:3)

         grad((ic-1)*3+1:(ic-1)*3+3) = grad((ic-1)*3+1:(ic-1)*3+3) + &
                                         dedc(1:3)

     END DO

 END SUBROUTINE hr_ang

!-------------------------------------------------------------------------------

!> @brief contribution of dihedral restraints to energy and gradient
!> @param[inout] ener       energy of the system (a.u.)
!> @param[inout] grad       gradient of the system (a.u.)
!> @param[in]    coor       coordinates of all atoms
!> @param[in]    limit      if `.TRUE.` max. gradient contribution from HR will be limited by `DEMAX`
 SUBROUTINE hr_dih(ener, grad, coor, limit)

     REAL(KIND=fp), INTENT(INOUT) :: ener, grad(:)
     REAL(KIND=fp), INTENT(IN) :: coor(:,:)
     LOGICAL, INTENT(IN) :: limit

     INTEGER :: &
         i, ia, ib, ic, id

     REAL(KIND=fp), DIMENSION(3):: &
         a, b, c, d, &
         ba, cb, dc, ca, db, &
         t, u, v, &
         dedt, dedu, &
         deda, dedb, dedc, dedd

     REAL(KIND=fp) :: &
         Fk, r0, &
         rt, ru, rv, rcb, &
         sine, cosine, ang, &
         delta, de

     DO i = 1, nDihR

         ia = IDDihR(1,i)
         ib = IDDihR(2,i)
         ic = IDDihR(3,i)
         id = IDDihR(4,i)
         Fk = dihRK(i)
         r0 = dihR0(i)

         a(1:3) = coor(1:3,ia)
         b(1:3) = coor(1:3,ib)
         c(1:3) = coor(1:3,ic)
         d(1:3) = coor(1:3,id)

         ba(:) = b(:) - a(:)
         cb(:) = c(:) - b(:)
         dc(:) = d(:) - c(:)

         ca(:) = c(:) - a(:)
         db(:) = d(:) - b(:)

         rcb=hr_veclen(cb)

         t = hr_crossprod(ba, cb)
         u = hr_crossprod(cb, dc)
         v = hr_crossprod(cb,  t)

         rt = hr_veclen(t)
         ru = hr_veclen(u)
         rv = hr_veclen(v)

         cosine = sum(t(:)*u(:))/(rt*ru)
         sine   = sum(v(:)*u(:))/(rv*ru)

         ang = atan2(sine, cosine)

         delta = ang - r0

         IF (delta<-PI) THEN
             delta = delta + TWOPI
         ELSE IF (delta>PI) THEN
             delta = delta - TWOPI
         END IF

!        Data for detailed printout
         dihCur(i)   = ang/PI*180
         dihDelta(i) = delta/PI*180
         dihE(i)     = Fk * delta * delta

         ener = ener + Fk * delta * delta
         de = 2.0_fp * Fk * delta

         IF (limit) THEN
             IF (abs(de)>DEMAX) THEN
                 de = sign(DEMAX,de)
             END IF
         END IF

         dedt(1:3) =  de * hr_crossprod(t,cb) /(rt*rt*rcb)
         dedu(1:3) = -de * hr_crossprod(u,cb) /(ru*ru*rcb)

         deda(1:3) = hr_crossprod(dedt,cb)
         dedb(1:3) = hr_crossprod(ca,dedt) + hr_crossprod(dedu,dc)
         dedc(1:3) = hr_crossprod(dedt,ba) + hr_crossprod(db,dedu)
         dedd(1:3) = hr_crossprod(dedu,cb)

         grad((ia-1)*3+1:(ia-1)*3+3) = grad((ia-1)*3+1:(ia-1)*3+3) + &
                                         deda(1:3)

         grad((ib-1)*3+1:(ib-1)*3+3) = grad((ib-1)*3+1:(ib-1)*3+3) + &
                                         dedb(1:3)

         grad((ic-1)*3+1:(ic-1)*3+3) = grad((ic-1)*3+1:(ic-1)*3+3) + &
                                         dedc(1:3)

         grad((id-1)*3+1:(id-1)*3+3) = grad((id-1)*3+1:(id-1)*3+3) + &
                                         dedd(1:3)

     END DO

 END SUBROUTINE hr_dih

!-------------------------------------------------------------------------------

!> @brief Scale vector so that its norm does not exceed vmax
!> @param[inout]    vec     vector to scale
!> @param[in]       vmax    vector norm limit
 SUBROUTINE hr_vscale(vec, vmax)
!dir$ attributes forceinline :: hr_vscale
     REAL(KIND=fp), INTENT(INOUT) :: vec(:)
     REAL(KIND=fp), INTENT(IN) :: vmax

     REAL(KIND=fp) :: fact

     fact = sqrt(sum(vec(:)**2))/vmax
     IF (fact>1.0_fp) vec(:) = vec(:)/fact

 END SUBROUTINE hr_vscale

!-------------------------------------------------------------------------------

!> @brief Calculate cross product (exterial product) of two 3D vectors
!> @param[in]    vec1    first vector
!> @param[in]    vec2    second vector
 FUNCTION hr_crossprod(vec1, vec2) RESULT(prod)
!dir$ attributes forceinline :: hr_crossprod
     REAL(KIND=fp), INTENT(IN) :: vec1(3), vec2(3)
     REAL(KIND=fp) :: prod(3)

     prod(1) = vec1(2)*vec2(3) - vec2(2)*vec1(3)
     prod(2) = vec1(3)*vec2(1) - vec2(3)*vec1(1)
     prod(3) = vec1(1)*vec2(2) - vec2(1)*vec1(2)

 END FUNCTION hr_crossprod

!-------------------------------------------------------------------------------

!> @brief Calculate vector length
!> @param[in]    vec    vector
 FUNCTION hr_veclen(vec)
!dir$ attributes forceinline :: hr_veclen
     REAL(KIND=fp), INTENT(IN) :: vec(:)
     REAL(KIND=fp) :: hr_veclen

     hr_veclen = sqrt(sum(vec(:)**2))

 END FUNCTION hr_veclen

END MODULE harmonic_restraints
