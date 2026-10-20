 MODULE mod_dft_gridint

    USE prec, ONLY: fp
    USE params, ONLY: dft_wt_der
    USE mx_limits, ONLY: mxatm, mxgtot, mxgrid, mxgridtyp
    USE mod_nosp_basis, ONLY: nosp_basis, bas_norm_matrix
    IMPLICIT NONE

    INTEGER, PARAMETER :: FLAG_SYNC = 23800

    REAL(KIND=fp), PARAMETER :: &
        ZERO = 0.0D+00, TWO = 2.0D+00, HALF = 0.5D+00, RLN10 = 2.30258D+00

    INTEGER, PARAMETER :: &
        IJX(84) = (/&
            1, &
            2, 1, 1, &
            3, 1, 1, 2, 2, 1, &
            4, 1, 1, 3, 3, 2, 1, 2, 1, 2, &
            5, 1, 1, 4, 4, 2, 1, 2, 1, 3, 3, 1, 3, 2, 2, &
            6, 1, 1, 5, 5, 2, 1, 2, 1, 4, 4, 3, 1, 3, 1, 4, 2, 2, 3, 3, 2, &
            7, 1, 1, 6, 6, 2, 1, 2, 1, 5, 5, 3, 1, 3, 1, 5, 2, 2, 4, 4, 1, 4, 4, 3, 2, 3, 2, 3/), &

        IJY(84) = (/&
            1, &
            1, 2, 1, &
            1, 3, 1, 2, 1, 2, &
            1, 4, 1, 2, 1, 3, 3, 1, 2, 2, &
            1, 5, 1, 2, 1, 4, 4, 1, 2, 3, 1, 3, 2, 3, 2, &
            1, 6, 1, 2, 1, 5, 5, 1, 2, 3, 1, 4, 4, 1, 3, 2, 4, 2, 3, 2, 3, &
            1, 7, 1, 2, 1, 6, 6, 1, 2, 3, 1, 5, 5, 1, 3, 2, 5, 2, 4, 1, 4, 3, 2, 4, 4, 2, 3, 3/), &

        IJZ(84) = (/&
            1, &
            1, 1, 2, &
            1, 1, 3, 1, 2, 2, &
            1, 1, 4, 1, 2, 1, 2, 3, 3, 2, &
            1, 1, 5, 1, 2, 1, 2, 4, 4, 1, 3, 3, 2, 2, 3, &
            1, 1, 6, 1, 2, 1, 2, 5, 5, 1, 3, 1, 3, 4, 4, 2, 2, 4, 2, 3, 3, &
            1, 1, 7, 1, 2, 1, 2, 6, 6, 1, 3, 1, 3, 5, 5, 2, 2, 5, 1, 4, 4, 2, 3, 2, 3, 4, 4, 3/)

    TYPE :: xc_engine_t
        PRIVATE
        REAL(KIND=fp), ALLOCATABLE :: compdata(:)
        !REAL(KIND=fp), ALLOCATABLE :: densdata(:)
        REAL(KIND=fp), ALLOCATABLE :: xyzw(:,:)
        REAL(KIND=fp), ALLOCATABLE :: rsqrd(:)
        REAL(KIND=fp), ALLOCATABLE :: angvl(:,:,:)

        REAL(KIND=fp), CONTIGUOUS, POINTER :: &
              aoMem(:,:,:)  => NULL() & !< AO memory
            , moMemA(:,:,:) => NULL() & !< MO memory (alpha)
            , moMemB(:,:,:) => NULL() & !< MO memory (beta)
            , densdata(:)   => NULL() & !< wts, rho, grad,..
            , tmp(:,:)      => NULL() & !< Temporary storage

            , aoV(:,:)      => NULL() & !< AO values
            , moVA(:,:)     => NULL() & !< MO values (alpha)
            , moVB(:,:)     => NULL() & !< MO values (beta)
            , aoG1(:,:,:)   => NULL() & !< AO gradient
            , aoG2(:,:,:)   => NULL() & !< AO 2nd der.
            , moG1A(:,:,:)  => NULL() & !< MO gradient (alpha)
            , moG2A(:,:,:)  => NULL() & !< MO 2nd der. (alpha)
            , moG1B(:,:,:)  => NULL() & !< MO gradient (beta)
            , moG2B(:,:,:)  => NULL() & !< MO 2nd der. (beta)
            , wts(:)        => NULL() & !< weights
!            , rhos(:)       => NULL() & !< density
!            , gradxc(:,:)   => NULL() & !< density gradient
!            , tauxc(:,:)    => NULL() & !< density Laplacian
            , dedr(:,:)     => NULL() & !< E(XC) LDA values
            , deds(:,:)     => NULL() & !< E(XC) GGA values
            , dedt(:,:)     => NULL() & !< E(XC) MGGA values

            , wfAlpha(:,:)  => NULL() & !< alpha spin wavefunction
            , wfBeta(:,:)   => NULL()   !< beta spin wavefunction

        LOGICAL :: isGGA     = .FALSE.
        LOGICAL :: needTau   = .FALSE.
        LOGICAL :: hasBeta   = .FALSE.
        LOGICAL :: isWFVecs  = .TRUE.   !< .TRUE.  - wfA and wfB are MO vectors
                                        !< .FALSE. - wfA and wfB are densities
        INTEGER :: numAOs    = 0
        INTEGER :: numPts    = 0
        INTEGER :: numAtoms  = 0
        INTEGER :: maxAngMom = 0
        INTEGER :: nDer      = 0
        INTEGER :: funTyp    = 0   !< 0 - LDA, 1 - GGA, 2 - MGGA

        INTEGER :: numAOVecs = 0

        INTEGER :: numOccAlpha = 0
        INTEGER :: numOccBeta  = 0

        REAL(KIND=fp) :: tolAO

        PROCEDURE(compute_density),      POINTER, PASS :: compRho   => NULL()
        PROCEDURE(compute_density_grad), POINTER, PASS :: compDRho  => NULL()
        PROCEDURE(compute_density_grad), POINTER, PASS :: compD2Rho => NULL()

    CONTAINS
        PROCEDURE :: init
        PROCEDURE :: resetOrbPointers
        PROCEDURE :: resetXCPointers
        PROCEDURE :: compAOs
        PROCEDURE :: compMOs
        PROCEDURE :: compDist
        PROCEDURE :: compXC
        PROCEDURE :: compXCgrad
        PROCEDURE :: updKS
        PROCEDURE :: compAtGradRho
        PROCEDURE :: compAtGradDRho
        PROCEDURE :: compAtGradTau
        PROCEDURE :: echo => echoVars
    END TYPE xc_engine_t

    ABSTRACT INTERFACE
        SUBROUTINE compute_density(self, iPt, rhoa, rhob)
            IMPORT
            CLASS(xc_engine_t) :: self
            INTEGER, INTENT(IN) :: iPt
            REAL(KIND=fp), INTENT(OUT) :: rhoa, rhob
        END SUBROUTINE
        SUBROUTINE compute_density_grad(self, iPt, gxa, gya, gza, gxb, gyb, gzb)
            IMPORT
            CLASS(xc_engine_t) :: self
            INTEGER, INTENT(IN) :: iPt
            REAL(KIND=fp), INTENT(OUT) :: gxa, gya, gza, gxb, gyb, gzb
        END SUBROUTINE
    END INTERFACE

    PRIVATE
    PUBLIC dmatd_blk
    PUBLIC derexc_blk
 CONTAINS

!> @brief Print parameters of the xc_engine_t instance
!> @author Vladimir Mironov
 SUBROUTINE echoVars(self)
    CLASS(xc_engine_t) :: self
    write(*,*) 'isGGA=', self%isGGA
    write(*,*) 'needTau =', self%needTau
    write(*,*) 'hasBeta =', self%hasBeta
    write(*,*) 'numAOs  =', self%numAOs
    write(*,*) 'numPts  =', self%numPts
    write(*,*) 'nDer    =', self%nDer
 END SUBROUTINE

!> @brief Initialize xc_engine_t instance
!> @param[in] numAOs    number of atomic orbitals in a basis
!> @param[in] nAt       number of atoms in a system
!> @param[in] maxAngMom maximum angular momentum of basis functions and their derivatives
!> @param[in] maxPts    maximum known number of non-zero points in a slice
!> @param[in] limPts    maximum possible number of points (i.e. max(nRad*nAng)) in a slice
!> @param[in] nDer      degree of energy derivative needed
!> @param[in] hasBeta   .TRUE. if open-shell calculation
!> @param[in] isGGA     .TRUE. if GGA/metaGGA functional
!> @param[in] needTau   .TRUE. if metaGGA functional
!> @param[in] tolAO     AO/derivative calculation cutoff
!> @param[in] vec_or_dens .TRUE./.FALSE. - wavefunction is MO vectors/density
!> @param[in] nOccAlpha number of occupied orbitals, alpha spin
!> @param[in] nOccBeta  number of occupied orbitals, beta spin
!> @param[in] wfAlpha   wavefunction, alpha spin
!> @param[in] wfBeta    wavefunction, beta spin
!> @author Vladimir Mironov
 SUBROUTINE init(self, numAOs, nAt, maxAngMom, maxPts, limPts, nDer, &
                 hasBeta, isGGA, needTau, tolAO, vec_or_dens, &
                 nOccAlpha, wfAlpha, &
                 nOccBeta, wfBeta)
    CLASS(xc_engine_t), TARGET :: self
    INTEGER, PARAMETER :: NDENS = 11
!   Will be possibly needed to use LibXC:
!   INTEGER, PARAMETER :: NDENS = 28
    INTEGER, INTENT(IN) :: numAOs, nAt, maxAngMom, nDer, nOccAlpha
    INTEGER, INTENT(IN) :: maxPts, limPts
    INTEGER, INTENT(IN), OPTIONAL :: nOccBeta
    LOGICAL, INTENT(IN) :: hasBeta, isGGA, needTau, vec_or_dens
    REAL(KIND=fp), INTENT(IN) :: tolAO
    REAL(KIND=fp), TARGET, CONTIGUOUS, INTENT(IN) :: wfAlpha(:,:)
    REAL(KIND=fp), TARGET, CONTIGUOUS, INTENT(IN), OPTIONAL :: wfBeta(:,:)

    INTEGER :: numVecs, numTmpVec

    self%funTyp = 0
    IF (isGGA)   self%funTyp = self%funTyp + 1
    IF (needTau) self%funTyp = self%funTyp + 1

    self%nDer = 0
    IF (isGGA) self%nDer = 1
    self%nDer = self%nDer + nDer

!    write(*,*) numAOs, nAt, maxAngMom, hasBeta, isGGA, needTau

!   Find out the amount of memory needed
    SELECT CASE (self%nDer)
    CASE (0)
        self%numAOVecs = 1
    CASE (1)
        self%numAOVecs = 4
    CASE (2)
        self%numAOVecs = 10
    CASE DEFAULT
        write(*,*) 'Invalid grad level in xc_engine_t % INIT'
        CALL abrt
    END SELECT

    numVecs = 2*self%numAOVecs
    IF (hasBeta) numVecs = 3*self%numAOVecs

    numTmpVec = 1
    IF (needTau) numTmpVec = 3

!   Allocate memory for XC calculations
    ALLOCATE(self%compdata((numAOs*(numVecs+numTmpVec)+NDENS)*maxPts))

    self%tmp(1:numAOs,1:numTmpVec*maxPts) => &
            self%compdata( numAOs*numVecs*maxPts+1:&
                           numAOs*(numVecs+numTmpVec)*maxPts)

    self%densdata(1:maxPts*NDENS) => &
            self%compdata( numAOs*(numVecs+numTmpVec)*maxPts+1:&
                           numAOs*(numVecs+numTmpVec)*maxPts+maxPts*NDENS)

!   Allocate memory for grid points storage
    ALLOCATE(self%xyzw(limPts,4))

!   Allocate memory for point-atom distances
    ALLOCATE(self%rsqrd(nAt))

!   ANGXVL=X**I, ANGYVL=Y**J, ANGZVL=Z**K, for needed values of I,J,K
!   Where (X,Y,Z) is the centre of a DFT grid point Ri
!   minus an atomic centre Ra: XYZ = (Ri - Ra)
!    write(*,*) 'allocating angvl'
    self%maxAngMom = maxAngMom + self%nDer
    ALLOCATE(self%angvl(self%maxAngMom,3,nAt), source=0.0_fp)

!   Set up common angvl values
    self%angvl(2,1:3,:) = 1.0_fp

!   Set up other runtime options
    self%numAOs    = numAOs
    self%numAtoms  = nAt

    self%isGGA   = isGGA
    self%needTau = needTau
    self%hasBeta = hasBeta

    self%tolAO = tolAO

!   Manage density/MO vectors
    self%isWFVecs = vec_or_dens
    self%wfAlpha => wfAlpha

    self%numOccAlpha = nOccAlpha

    IF (vec_or_dens) THEN
        self%compRho   => compRhoMO
        self%compDRho  => compDRhoMO
        self%compD2Rho => compD2RhoMO
    ELSE
        self%compRho   => compRhoAO
        self%compDRho  => compDRhoAO
        self%compD2Rho => compD2RhoAO
    END IF

    IF (present(nOccBeta).AND.self%hasBeta) THEN
        self%numOccBeta =  nOccBeta
        self%wfBeta     => wfBeta
    ELSE
        self%numOccBeta  =  nOccAlpha
        self%wfBeta      => wfAlpha
    END IF

 END SUBROUTINE

!> @brief Adjust internal AO/MO memory storage for a given
!>  number of grid points
!> @param[in] numPts    number of grid points
!> @author Vladimir Mironov
 SUBROUTINE resetOrbPointers(self, numPts)
    CLASS(xc_engine_t), TARGET :: self
    INTEGER, INTENT(IN) :: numPts
    INTEGER :: ind

    ind = 0
    self%aoMem(1:self%numAOs,1:numPts,1:self%numAOVecs) => &
            self%compdata(ind+1:ind+self%numAOs*numPts*self%numAOVecs)

    ind = ind+self%numAOs*numPts*self%numAOVecs

    self%moMemA(1:self%numAOs,1:numPts,1:self%numAOVecs) => &
            self%compdata(ind+1:ind+self%numAOs*numPts*self%numAOVecs)

    ind = ind+self%numAOs*numPts*self%numAOVecs

    IF (self%hasBeta) THEN
        self%moMemB(1:self%numAOs,1:numPts,1:self%numAOVecs) => &
                self%compdata(ind+1:ind+self%numAOs*numPts*self%numAOVecs)
        ind = ind+self%numAOs*numPts*self%numAOVecs
    END IF


    SELECT CASE (self%nDer)
    CASE (0)
        self%aoV    => self%aoMem(:,:,1)
        self%moVA    => self%moMemA(:,:,1)

        IF (self%hasBeta) THEN
            self%moMemB(1:self%numAOs,1:numPts,1:1) => &
                    self%compdata(self%numAOs*numPts*2+1:self%numAOs*numPts*3)
            self%moVB   => self%moMemB(:,:,1)
        END IF
    CASE (1)
        self%aoV    => self%aoMem(:,:,1)
        self%aoG1   => self%aoMem(:,:,2:4)

        self%moVA    => self%moMemA(:,:,1)
        self%moG1A   => self%moMemA(:,:,2:4)

        IF (self%hasBeta) THEN
            self%moVB   => self%moMemB(:,:,1)
            self%moG1B  => self%moMemB(:,:,2:4)
        END IF
    CASE (2)
        self%aoV    => self%aoMem(:,:,1)
        self%aoG1   => self%aoMem(:,:,2:4)
        self%aoG2   => self%aoMem(:,:,5:10)

        self%moVA    => self%moMemA(:,:,1)
        self%moG1A   => self%moMemA(:,:,2:4)
        self%moG2A   => self%moMemA(:,:,5:10)

        IF (self%hasBeta) THEN
            self%moVB   => self%moMemB(:,:,1)
            self%moG1B  => self%moMemB(:,:,2:4)
            self%moG2B  => self%moMemB(:,:,5:10)
        END IF
    END SELECT
    self%wts(1:numPts) => self%densdata(1:numPts)
 END SUBROUTINE

!> @brief Adjust internal XC memory storage for a given
!>  number of grid points
!> @param[in] numPts    number of relevant grid points
!> @author Vladimir Mironov
 SUBROUTINE resetXCPointers(self, numPts)
    CLASS(xc_engine_t), TARGET :: self
    INTEGER, INTENT(IN) :: numPts
    self%numPts = numPts
    self%wts(1:numPts)      => self%densdata(          1:   numPts)

    self%dedr(1:numPts,1:2) => self%densdata( 1*numPts+1: 3*numPts)

    !IF (self%isGGA.OR.self%needTau) THEN
    IF (self%funTyp>0) THEN
      self%deds(1:numPts,1:6)    => self%densdata( 3*numPts+1: 9*numPts)
      self%dedt(1:numPts,1:2)    => self%densdata( 9*numPts+1:11*numPts)
!     Reserve for LibXC:
      !self%rhos(1:numPts*2)      => self%densdata(11*numPts+1:13*numPts)
      !self%gradxc(1:numPts,1:9)  => self%densdata(13*numPts+1:22*numPts)
      !self%tauxc(1:numPts,1:6)   => self%densdata(22*numPts+1:28*numPts)
    END IF
 END SUBROUTINE

!> @brief Compute atomic orbital values/gradient/hessian in a grid point
!> @param[in]  atmxvec    array of atomic X coordinates
!> @param[in]  atmyvec    array of atomic Y coordinates
!> @param[in]  atmzvec    array of atomic Z coordinates
!> @param[in]  iPtIn      index of the point in self%xyzw array
!> @param[in]  iPtOut     index of the point in AO/MO arrays
!> @param[in]  iAtom      index of the atom who owns the point
!> @param[out] nnz        number of non-zero AOs in the point
!> @author Vladimir Mironov
 SUBROUTINE compAOs(self, &
                 atmxvec, atmyvec, atmzvec, &
                 iPtIn, iPtOut, iAtom, nnz)
    CLASS(xc_engine_t) :: self
    INTEGER, INTENT(IN) :: iAtom
    INTEGER, INTENT(IN) :: iPtIn, iPtOut
    REAL(KIND=fp), CONTIGUOUS, INTENT(IN) :: &
            atmxvec(:,:), atmyvec(:,:), atmzvec(:,:)
    INTEGER, INTENT(OUT) :: nnz

    CALL self%compDist(atmxvec, atmyvec, atmzvec, iPtIn, iAtom)

    SELECT CASE (self%nDer)
    CASE (0)
        CALL compAOv(nosp_basis, self%tolAO, &
                self%angvl, self%rsqrd, nnz, &
                self%aoV(:,iPtOut))
    CASE (1)
        CALL compAOvg(nosp_basis, self%tolAO, &
                 self%angvl, self%rsqrd, nnz, &

                 self%aoV(:,iPtOut), &

                 self%aoG1(:,iPtOut,1), &
                 self%aoG1(:,iPtOut,2), &
                 self%aoG1(:,iPtOut,3))

    CASE (2)
        CALL compAOvgg(nosp_basis, self%tolAO, &
                 self%angvl, self%rsqrd, nnz, &

                 self%aoV(:,iPtOut), &

                 self%aoG1(:,iPtOut,1), &
                 self%aoG1(:,iPtOut,2), &
                 self%aoG1(:,iPtOut,3), &

                 self%aoG2(:,iPtOut,1), &
                 self%aoG2(:,iPtOut,2), &
                 self%aoG2(:,iPtOut,3), &
                 self%aoG2(:,iPtOut,4), &
                 self%aoG2(:,iPtOut,5), &
                 self%aoG2(:,iPtOut,6))

    CASE DEFAULT
        WRITE(*,*) 'Invalid grad level in xc_engine_t % COMPAOS'
        CALL abrt
    END SELECT

 END SUBROUTINE

!> @brief Compute distances between the point and atoms in a system
!> @param[in]  atmxvec    array of atomic X coordinates
!> @param[in]  atmyvec    array of atomic Y coordinates
!> @param[in]  atmzvec    array of atomic Z coordinates
!> @param[in]  iPt        index of the point in self%xyzw array
!> @param[in]  iAtom      index of the atom who owns the point
!> @author Vladimir Mironov
 SUBROUTINE compDist(self, &
                 atmxvec, atmyvec, atmzvec, &
                 iPt, iAtom)

    CLASS(xc_engine_t) :: self
    INTEGER, INTENT(IN) :: ipt, iAtom
    REAL(KIND=fp), CONTIGUOUS, INTENT(IN) :: &
            atmxvec(:,:), atmyvec(:,:), atmzvec(:,:)

    REAL(KIND=fp) :: xyzcd(3)
    INTEGER :: i, j

    DO j = 1, self%numAtoms
        ASSOCIATE ( d2 => self%rsqrd(j), av => self%angvl(:,:,j), &
                    xd => self%xyzw(iPt,1), yd => self%xyzw(iPt,2), &
                    zd => self%xyzw(iPt,3))
            xyzcd(1) = xd-atmxvec(j,iAtom)
            xyzcd(2) = yd-atmyvec(j,iAtom)
            xyzcd(3) = zd-atmzvec(j,iAtom)
            d2 = sum(xyzcd**2)
            av(3,1:3) = xyzcd
            DO i = 4, self%maxAngMom
                av(i,1:3) = av(i-1,1:3) * xyzcd(1:3)
            END DO
        END ASSOCIATE
    END DO

 END SUBROUTINE

!> @brief Transform AOs to "MOs"
!> @details Multiply AO vector to the MO coefficient matrix or density matrix.
!>  True MOs are only obtained in the former case.
!> @author Vladimir Mironov
 SUBROUTINE compMOs(self)
    CLASS(xc_engine_t) :: self
    INTEGER :: nAlpha, nBeta, nVecs

    nVecs = min(self%numAOVecs, 4) ! Don't transform second derivatives

    IF (self%isWFVecs) THEN
        nAlpha = self%numOccAlpha
        nBeta  = self%numOccBeta
        CALL dgemm('N', 'N', &
            nAlpha, nVecs*ubound(self%aoMem,2), self%numAOs, &
            1.0_fp, self%wfAlpha(1,1), self%numAOs, &
                    self%aoV(1,1),     self%numAOs, &
            0.0_fp, self%moVA(1,1),    self%numAOs)
        IF (self%hasBeta) &
            CALL dgemm('N', 'N', &
                nBeta, nVecs*ubound(self%aoMem,2), self%numAOs, &
                1.0_fp, self%wfBeta(1,1), self%numAOs, &
                        self%aoV(1,1),    self%numAOs, &
                0.0_fp, self%moVB(1,1),   self%numAOs)
    ELSE
        CALL dsymm('L', 'U', &
            self%numAOs, nVecs*ubound(self%aoMem,2), &
            1.0_fp, self%wfAlpha(1,1), self%numAOs, &
                    self%aoV(1,1),     self%numAOs, &
            0.0_fp, self%moVA(1,1),    self%numAOs)
        IF (self%hasBeta) &
            CALL dsymm('L', 'U', &
                self%numAOs, nVecs*ubound(self%aoMem,2), &
                1.0_fp, self%wfBeta(1,1), self%numAOs, &
                        self%aoV(1,1),    self%numAOs, &
                0.0_fp, self%moVB(1,1),   self%numAOs)
    END IF
 END SUBROUTINE

!!   Will be possibly needed to use LibXC:
! SUBROUTINE compRhoAll(self, nOccAlpha, numPts, totele, totgradx, totgrady, totgradz, totkin)
!    CLASS(xc_engine_t) :: self
!    INTEGER, INTENT(IN) :: nOccAlpha, numPts
!    REAL(KIND=fp) :: totgradx, totgrady, totgradz, totele, totkin
!    INTEGER :: i
!
!    ! Form density
!    IF (self%isGGA.OR.self%needTau) THEN
!        DO i = 1, numPts
!            self%rhos(i)     =   dot_product(self%moVA(:nOccAlpha,i), self%moVA(:nOccAlpha,i))
!            self%gradxc(i,1) = 2*dot_product(self%moG1A(:nOccAlpha,i,1), self%moVA(:nOccAlpha,i))
!            self%gradxc(i,2) = 2*dot_product(self%moG1A(:nOccAlpha,i,2), self%moVA(:nOccAlpha,i))
!            self%gradxc(i,3) = 2*dot_product(self%moG1A(:nOccAlpha,i,3), self%moVA(:nOccAlpha,i))
!
!            self%rhos(i+numPts) = self%rhos(i) !rob
!            self%gradxc(i,4:6) = self%gradxc(i,1:3) ! grad-b
!
!            totele = totele + self%wts(i)*(self%rhos(i)+self%rhos(i+numPts))
!        END DO
!!       Form gradient invariant (GRAD dot GRAD)
!        DO i = 1, numPts
!            self%gradxc(i,7) = sum(self%gradxc(i,1:3)**2)
!            self%gradxc(i,8) = sum(self%gradxc(i,4:6)**2)
!            self%gradxc(i,9) = sum(self%gradxc(i,1:3)*self%gradxc(i,4:6))
!
!!           The total electron density gradient
!            totgradx = totgradx + self%wts(i)*self%gradxc(i,7)
!            totgrady = totgrady + self%wts(i)*self%gradxc(i,8)
!            totgradz = totgradz + self%wts(i)*self%gradxc(i,9)
!        END DO
!        IF (self%needTau) THEN
!            DO i = 1, numPts
!                self%tauxc(i,1) = dot_product(self%moG1A(:nOccAlpha,i,1), self%moG1A(:nOccAlpha,i,1))
!                self%tauxc(i,2) = dot_product(self%moG1A(:nOccAlpha,i,2), self%moG1A(:nOccAlpha,i,2))
!                self%tauxc(i,3) = dot_product(self%moG1A(:nOccAlpha,i,3), self%moG1A(:nOccAlpha,i,3))
!                self%tauxc(i,4:6) = self%tauxc(i,1:3)
!                totkin = totkin + 0.5_fp*self%wts(i)*sum(self%tauxc(i,:))
!            END DO
!        END IF
!    ELSE
!        DO i = 1, numPts
!            self%rhos(i) = dot_product(self%moVA(:nOccAlpha,i), &
!                                       self%moVA(:nOccAlpha,i))
!            self%rhos(i+numPts) = self%rhos(i) !rob
!!           The total electron density
!            totele = totele + self%wts(i)*(self%rhos(i)+self%rhos(i+numPts))
!        END DO
!    END IF
!
! END SUBROUTINE

!> @brief Compute XC contribution to the gradient
!> @param[inout]  bfGrad     array of gradient contributinos per basis function
!> @param[inout]  exec       XC energy integral
!> @param[inout]  ecorl      correlation energy integral
!> @param[inout]  totele     density integral == number of electrons
!> @param[inout]  totkin     kinetic energy integral
!> @author Vladimir Mironov
 SUBROUTINE compXCgrad(self, bfGrad, exec, ecorl, totele, totkin)
    USE funclib, ONLY: funcl, funfl
    use libxc, only: use_libxc, libxc_calc
    CLASS(xc_engine_t) :: self
    REAL(KIND=fp), CONTIGUOUS, INTENT(INOUT) :: bfGrad(:,:)
    REAL(KIND=fp), INTENT(OUT) :: exec, ecorl, totele, totkin

    INTEGER :: i
    REAL(KIND=fp) :: xalpha, xgrd, ecf
    REAL(KIND=fp) :: &
            rhoa, gradxa, gradya, gradza, tauxa, tauya, tauza, &
            rhob, gradxb, gradyb, gradzb, tauxb, tauyb, tauzb, &
            grdaa, grdbb, grdab, &
            vxca, vxcb, dsa(3), dsb(3), dmgga, dmggb
!   Compute MOs
    !write(*,*) 'comp MOs'
    CALL self%compMOs

!   Reserve for LibXC:
    !write(*,*) 'comp dens'
    !CALL self%compRhoAll(nOccAlpha, numPts, totele, totgradx, totgrady, totgradz, totkin)

    DO i = 1, self%numPts
    ASSOCIATE ( wt    => self%wts(i) )

        vxca  = ZERO
        vxcb  = ZERO
        IF (self%isGGA) THEN
            dsa  = ZERO
            dsb  = ZERO
        END IF

        IF (self%needTau) THEN
            dmgga = ZERO
            dmggb = ZERO
        END IF

        xalpha = ZERO
        xgrd   = ZERO
        ecf    = ZERO

        !LDA will not assign any values to density gradients
        grdaa = ZERO
        grdbb = ZERO
        grdab = ZERO

        ! electronic density
        CALL self%compRho(i, rhoa, rhob)

        totele = totele + wt*(rhoa+rhob)

        IF (self%isGGA) THEN
            ! electronic density 1st derivative
            CALL self%compDRho(i, gradxa, gradya, gradza, gradxb, gradyb, gradzb)

            grdaa  = gradxa*gradxa + gradya*gradya + gradza*gradza
            grdbb  = gradxb*gradxb + gradyb*gradyb + gradzb*gradzb
            grdab  = gradxa*gradxb + gradya*gradyb + gradza*gradzb

        END IF

        IF (self%needTau) THEN
            ! electronic density 2nd derivative
            CALL self%compD2Rho(i, tauxa, tauya, tauza, tauxb, tauyb, tauzb)
            totkin = totkin + 0.5_fp*wt*(tauxa+tauya+tauza+tauxb+tauyb+tauzb)
        END IF

        ! compute XC in point
        IF (use_libxc) THEN
           call libxc_calc(wt, &
                           rhoa,rhob, &
                           grdaa,grdab,grdbb, &
                           gradxa,gradya,gradza,gradxb,gradyb,gradzb, &
                           tauxa,tauya,tauza,tauxb,tauyb,tauzb, &
                           xalpha,xgrd,ecf, &
                           vxca,vxcb, &
                           dsa(1),dsa(2),dsa(3),dsb(1),dsb(2),dsb(3), &
                           dmgga,dmggb)
        ELSE IF (funcl) THEN
            CALL ccalcexc(rhoa,rhob,wt,grdaa,grdbb,grdab, &
                          gradxa,gradya,gradza,gradxb,gradyb,gradzb, &
                          xalpha,xgrd,vxca,dsa(1),dsa(2),dsa(3), &
                          vxcb,dsb(1),dsb(2),dsb(3),ecf)
        ELSE IF (funfl) THEN
            CALL fcalcexc(rhoa,rhob,wt,grdaa,grdbb,grdab, &
                          gradxa,gradya,gradza,gradxb,gradyb,gradzb, &
                          xalpha,xgrd,vxca,dsa(1),dsa(2),dsa(3), &
                          vxcb,dsb(1),dsb(2),dsb(3),ecf)
        ELSE
            CALL calcexc(rhoa,rhob,wt,grdaa,grdbb,grdab, &
                         gradxa,gradya,gradza,gradxb,gradyb,gradzb, &
                         xalpha,xgrd,vxca,dsa(1),dsa(2),dsa(3), &
                         vxcb,dsb(1),dsb(2),dsb(3),ecf, &
                         tauxa,tauya,tauza,dmgga, &
                         tauxb,tauyb,tauzb,dmggb)
        END IF

        exec = exec + xalpha + xgrd + ecf
        ecorl = ecorl + ecf

!       Scale values by quadrature weights and compute gradients

        ! LDA gradient
        vxca  = vxca * wt
        vxcb  = vxcb * wt
        CALL self%compAtGradRho(i, bfGrad, vxca, vxcb)

        IF (self%isGGA) THEN
            ! GGA gradient
            dsa  = dsa * wt
            dsb  = dsb * wt
            CALL self%compAtGradDRho(i, bfGrad, dsa, dsb)
        END IF

        IF (self%needTau) THEN
            ! meta-GGA gradient
            dmgga = dmgga * wt
            dmggb = dmggb * wt
            CALL self%compAtGradTau(i, bfGrad, dmgga, dmggb)
        END IF

    END ASSOCIATE
    END DO

 END SUBROUTINE

!> @brief Compute XC contribution to the energy
!> @param[inout]  bfGrad     array of gradient contributinos per basis function
!> @param[inout]  exec       XC energy integral
!> @param[inout]  ecorl      correlation energy integral
!> @param[inout]  totele     density integral == number of electrons
!> @param[inout]  totkin     kinetic energy integral
!> @param[inout]  togradxyz  density gradient integral
!> @author Vladimir Mironov
 SUBROUTINE compXC(self, exec, ecorl, totele, totkin, totgradxyz)
    USE funclib, ONLY: funcl, funfl
    use libxc, only: use_libxc, libxc_calc
    CLASS(xc_engine_t) :: self
    REAL(KIND=fp), INTENT(OUT) :: exec, ecorl, totele, totkin
    REAL(KIND=fp), INTENT(OUT), OPTIONAL :: totgradxyz(:)

    INTEGER :: i
    REAL(KIND=fp) :: xalpha, xgrd, ecf
    REAL(KIND=fp) :: &
            rhoa, gradxa, gradya, gradza, tauxa, tauya, tauza, &
            rhob, gradxb, gradyb, gradzb, tauxb, tauyb, tauzb, &
            grdaa, grdbb, grdab

!   Compute MOs
    CALL self%compMOs

!   Reserve for LibXC:
    !CALL self%compRhoAll(nOccAlpha, numPts, totele, totgradx, totgrady, totgradz, totkin)

    self%dedr = 0.0_fp
    IF (self%isGGA)   self%deds = 0.0_fp
    IF (self%needTau) self%dedt = 0.0_fp

    DO i = 1, self%numPts
    ASSOCIATE ( wt    => self%wts(i) &
              , vxca  => self%dedr(i,1) &
              , vxcb  => self%dedr(i,2) &
              , dsax  => self%deds(i,1) &
              , dsay  => self%deds(i,2) &
              , dsaz  => self%deds(i,3) &
              , dsbx  => self%deds(i,4) &
              , dsby  => self%deds(i,5) &
              , dsbz  => self%deds(i,6) &
              , dmgga => self%dedt(i,1) &
              , dmggb => self%dedt(i,2) &
              )

!        vxca  = ZERO
!        vxcb  = ZERO
!        IF (self%isGGA) THEN
!            dsax  = ZERO
!            dsay  = ZERO
!            dsaz  = ZERO
!            dsbx  = ZERO
!            dsby  = ZERO
!            dsbz  = ZERO
!        END IF
!
!        IF (self%needTau) THEN
!            dmgga = ZERO
!            dmggb = ZERO
!        END IF

        xalpha = ZERO
        xgrd   = ZERO
        ecf    = ZERO

        !LDA will not assign any values to density gradients
        grdaa = ZERO
        grdbb = ZERO
        grdab = ZERO

        ! electron density
        CALL self%compRho(i, rhoa, rhob)

        totele = totele + wt*(rhoa+rhob)

        IF (self%isGGA) THEN
            ! electron density 1st derivative
            CALL self%compDRho(i, gradxa, gradya, gradza, gradxb, gradyb, gradzb)

            grdaa  = gradxa*gradxa + gradya*gradya + gradza*gradza
            grdbb  = gradxb*gradxb + gradyb*gradyb + gradzb*gradzb
            grdab  = gradxa*gradxb + gradya*gradyb + gradza*gradzb

!           The total electron density gradient
            IF (present(totgradxyz)) THEN
                totgradxyz(1)= totgradxyz(1) + wt*grdaa
                totgradxyz(2)= totgradxyz(2) + wt*grdbb
                totgradxyz(3)= totgradxyz(3) + wt*grdab
            END IF

        END IF

        IF (self%needTau) THEN
            ! electron density 2nd derivative
            CALL self%compD2Rho(i, tauxa, tauya, tauza, tauxb, tauyb, tauzb)
            totkin = totkin + 0.5_fp*wt*(tauxa+tauya+tauza+tauxb+tauyb+tauzb)
        END IF

        IF (use_libxc) THEN
            CALL libxc_calc(wt, &
                            rhoa,rhob, &
                            grdaa,grdab,grdbb, &
                            gradxa,gradya,gradza,gradxb,gradyb,gradzb, &
                            tauxa,tauya,tauza,tauxb,tauyb,tauzb, &
                            xalpha,xgrd,ecf, &
                            vxca,vxcb, &
                            dsax,dsay,dsaz,dsbx,dsby,dsbz, &
                            dmgga,dmggb)
        ELSE IF (funcl) THEN
            CALL ccalcexc(rhoa,rhob,wt,grdaa,grdbb,grdab, &
                          gradxa,gradya,gradza,gradxb,gradyb,gradzb, &
                          xalpha,xgrd,vxca,dsax,dsay,dsaz, &
                          vxcb,dsbx,dsby,dsbz,ecf)
        ELSE IF (funfl) THEN
            CALL fcalcexc(rhoa,rhob,wt,grdaa,grdbb,grdab, &
                          gradxa,gradya,gradza,gradxb,gradyb,gradzb, &
                          xalpha,xgrd,vxca,dsax,dsay,dsaz, &
                          vxcb,dsbx,dsby,dsbz,ecf)
        ELSE
            CALL calcexc(rhoa,rhob,wt,grdaa,grdbb,grdab, &
                         gradxa,gradya,gradza,gradxb,gradyb,gradzb, &
                         xalpha,xgrd,vxca,dsax,dsay,dsaz, &
                         vxcb,dsbx,dsby,dsbz,ecf, &
                         tauxa,tauya,tauza,dmgga, &
                         tauxb,tauyb,tauzb,dmggb)
        END IF

!       Scale values by quadrature weights
        vxca  = vxca * wt
        vxcb  = vxcb * wt
        IF (self%isGGA) THEN
            dsax  = dsax * wt
            dsay  = dsay * wt
            dsaz  = dsaz * wt
            dsbx  = dsbx * wt
            dsby  = dsby * wt
            dsbz  = dsbz * wt
        END IF

        IF (self%needTau) THEN
            dmgga = dmgga * wt
            dmggb = dmggb * wt
        END IF

        exec = exec + xalpha + xgrd + ecf
        ecorl = ecorl + ecf

    END ASSOCIATE
    END DO

 END SUBROUTINE

! SUBROUTINE compRho(self, iPt, rhoa, rhob)
!    CLASS(xc_engine_t) :: self
!    INTEGER, INTENT(IN) :: iPt
!    REAL(KIND=fp), INTENT(OUT) :: rhoa, rhob
!
!    IF (self%isWFVecs) THEN
!        CALL compRhoMO(self, iPt, rhoa, rhob)
!    ELSE
!        CALL compRhoAO(self, iPt, rhoa, rhob)
!    END IF
!
! END SUBROUTINE

!> @brief Compute electronic density in a grid point, AO-driven calculation
!> @param[in]  iPt      index of a grid point
!> @param[out] rhoa     electronic density, alpha-spin
!> @param[out] rhob     electronic density, beta-spin
!> @author Vladimir Mironov
 SUBROUTINE compRhoAO(self, iPt, rhoa, rhob)

    CLASS(xc_engine_t) :: self
    INTEGER, INTENT(IN) :: iPt
    REAL(KIND=fp), INTENT(OUT) :: rhoa, rhob

    IF (self%hasBeta) THEN
        rhoa = dot_product(self%aoV(:,iPt), self%moVA(:,iPt))
        rhob = dot_product(self%aoV(:,iPt), self%moVB(:,iPt))
    ELSE
        rhoa = 0.5_fp*dot_product(self%aoV(:,iPt), self%moVA(:,iPt))
        rhob = rhoa
    END IF

 END SUBROUTINE

!> @brief Compute electronic density in a grid point, MO-driven calculation
!> @param[in]  iPt      index of a grid point
!> @param[out] rhoa     electronic density, alpha-spin
!> @param[out] rhob     electronic density, beta-spin
!> @author Vladimir Mironov
 SUBROUTINE compRhoMO(self, iPt, rhoa, rhob)

    CLASS(xc_engine_t) :: self
    INTEGER, INTENT(IN) :: iPt
    REAL(KIND=fp), INTENT(OUT) :: rhoa, rhob

    ASSOCIATE ( noa => self%numOccAlpha, nob => self%numOccBeta)

        rhoa = dot_product(self%moVA(:noa,iPt), self%moVA(:noa,iPt))
        IF (self%hasBeta) THEN
            rhob = dot_product(self%moVB(:nob,iPt), self%moVB(:nob,iPt))
        ELSE
            rhob = rhoa
        END IF

    END ASSOCIATE

 END SUBROUTINE

! SUBROUTINE compDRho(self, iPt, gradxa, gradya, gradza, gradxb, gradyb, gradzb)
!    CLASS(xc_engine_t) :: self
!    INTEGER, INTENT(IN) :: iPt
!    REAL(KIND=fp), INTENT(OUT) :: &
!            gradxa, gradya, gradza, &
!            gradxb, gradyb, gradzb
!
!    IF (self%isWFVecs) THEN
!        CALL compDRhoMO(self, iPt, gradxa, gradya, gradza, gradxb, gradyb, gradzb)
!    ELSE
!        CALL compDRhoAO(self, iPt, gradxa, gradya, gradza, gradxb, gradyb, gradzb)
!    END IF
!
! END SUBROUTINE

!> @brief Compute electronic density gradient in a grid point, AO-driven calculation
!> @param[in]  iPt      index of a grid point
!> @param[out] gradxa   dRho/dX, alpha-spin
!> @param[out] gradya   dRho/dY, alpha-spin
!> @param[out] gradza   dRho/dZ, alpha-spin
!> @param[out] gradxb   dRho/dX, beta-spin
!> @param[out] gradyb   dRho/dY, beta-spin
!> @param[out] gradzb   dRho/dZ, beta-spin
!> @author Vladimir Mironov
 SUBROUTINE compDRhoAO(self, iPt, &
            gradxa, gradya, gradza, &
            gradxb, gradyb, gradzb)

    CLASS(xc_engine_t) :: self
    INTEGER, INTENT(IN) :: iPt
    REAL(KIND=fp), INTENT(OUT) :: &
            gradxa, gradya, gradza, &
            gradxb, gradyb, gradzb

    IF (self%hasBeta) THEN
        gradxa = 2*dot_product(self%aoG1(:,iPt,1), self%moVA(:,iPt))
        gradya = 2*dot_product(self%aoG1(:,iPt,2), self%moVA(:,iPt))
        gradza = 2*dot_product(self%aoG1(:,iPt,3), self%moVA(:,iPt))

        gradxb = 2*dot_product(self%aoG1(:,iPt,1), self%moVB(:,iPt))
        gradyb = 2*dot_product(self%aoG1(:,iPt,2), self%moVB(:,iPt))
        gradzb = 2*dot_product(self%aoG1(:,iPt,3), self%moVB(:,iPt))
    ELSE
        gradxa = dot_product(self%aoG1(:,iPt,1), self%moVA(:,iPt))
        gradya = dot_product(self%aoG1(:,iPt,2), self%moVA(:,iPt))
        gradza = dot_product(self%aoG1(:,iPt,3), self%moVA(:,iPt))

        gradxb = gradxa
        gradyb = gradya
        gradzb = gradza
    END IF

 END SUBROUTINE

!> @brief Compute electronic density gradient in a grid point, MO-driven calculation
!> @param[in]  iPt      index of a grid point
!> @param[out] gradxa   dRho/dX, alpha-spin
!> @param[out] gradya   dRho/dY, alpha-spin
!> @param[out] gradza   dRho/dZ, alpha-spin
!> @param[out] gradxb   dRho/dX, beta-spin
!> @param[out] gradyb   dRho/dY, beta-spin
!> @param[out] gradzb   dRho/dZ, beta-spin
!> @author Vladimir Mironov
 SUBROUTINE compDRhoMO(self, iPt, &
            gradxa, gradya, gradza, &
            gradxb, gradyb, gradzb)

    CLASS(xc_engine_t) :: self
    INTEGER, INTENT(IN) :: iPt
    REAL(KIND=fp), INTENT(OUT) :: &
            gradxa, gradya, gradza, &
            gradxb, gradyb, gradzb

    ASSOCIATE ( noa => self%numOccAlpha, nob => self%numOccBeta)
        gradxa = 2*dot_product(self%moG1A(:noa,iPt,1), self%moVA(:noa,iPt))
        gradya = 2*dot_product(self%moG1A(:noa,iPt,2), self%moVA(:noa,iPt))
        gradza = 2*dot_product(self%moG1A(:noa,iPt,3), self%moVA(:noa,iPt))
        IF (self%hasBeta) THEN
            gradxb = 2*dot_product(self%moG1B(:nob,iPt,1), self%moVB(:nob,iPt))
            gradyb = 2*dot_product(self%moG1B(:nob,iPt,2), self%moVB(:nob,iPt))
            gradzb = 2*dot_product(self%moG1B(:nob,iPt,3), self%moVB(:nob,iPt))
        ELSE
            gradxb = gradxa
            gradyb = gradya
            gradzb = gradza
        END IF
    END ASSOCIATE

 END SUBROUTINE

! SUBROUTINE compD2Rho(self, iPt, tauxa, tauya, tauza, tauxb, tauyb, tauzb)
!    CLASS(xc_engine_t) :: self
!    INTEGER, INTENT(IN) :: iPt
!    REAL(KIND=fp), INTENT(OUT) :: tauxa, tauya, tauza, tauxb, tauyb, tauzb
!
!    IF (self%isWFVecs) THEN
!        CALL compD2RhoMO(self, iPt, tauxa, tauya, tauza, tauxb, tauyb, tauzb)
!    ELSE
!        CALL compD2RhoAO(self, iPt, tauxa, tauya, tauza, tauxb, tauyb, tauzb)
!    END IF
!
! END SUBROUTINE

!> @brief Compute electronic density 2nd derivatives in a grid point, AO-driven calculation
!> @param[in]  iPt     index of a grid point
!> @param[out] tauxa   d^2(Rho)/dX^2, alpha-spin
!> @param[out] tauya   d^2(Rho)/dY^2, alpha-spin
!> @param[out] tauza   d^2(Rho)/dZ^2, alpha-spin
!> @param[out] tauxb   d^2(Rho)/dX^2, beta-spin
!> @param[out] tauyb   d^2(Rho)/dY^2, beta-spin
!> @param[out] tauzb   d^2(Rho)/dZ^2, beta-spin
!> @author Vladimir Mironov
 SUBROUTINE compD2RhoAO(self, iPt, tauxa, tauya, tauza, tauxb, tauyb, tauzb)

    CLASS(xc_engine_t) :: self
    INTEGER, INTENT(IN) :: iPt
    REAL(KIND=fp), INTENT(OUT) :: tauxa, tauya, tauza, tauxb, tauyb, tauzb

    IF (self%hasBeta) THEN
        tauxa  = dot_product(self%aoG1(:,iPt,1), self%moG1A(:,iPt,1))
        tauya  = dot_product(self%aoG1(:,iPt,2), self%moG1A(:,iPt,2))
        tauza  = dot_product(self%aoG1(:,iPt,3), self%moG1A(:,iPt,3))

        tauxb  = dot_product(self%aoG1(:,iPt,1), self%moG1B(:,iPt,1))
        tauyb  = dot_product(self%aoG1(:,iPt,2), self%moG1B(:,iPt,2))
        tauzb  = dot_product(self%aoG1(:,iPt,3), self%moG1B(:,iPt,3))
    ELSE
        tauxa  = 0.5_fp*dot_product(self%aoG1(:,iPt,1), self%moG1A(:,iPt,1))
        tauya  = 0.5_fp*dot_product(self%aoG1(:,iPt,2), self%moG1A(:,iPt,2))
        tauza  = 0.5_fp*dot_product(self%aoG1(:,iPt,3), self%moG1A(:,iPt,3))

        tauxb  = tauxa
        tauyb  = tauya
        tauzb  = tauza
    END IF

 END SUBROUTINE

!> @brief Compute electronic density 2nd derivatives in a grid point, MO-driven calculation
!> @param[in]  iPt     index of a grid point
!> @param[out] tauxa   d^2(Rho)/dX^2, alpha-spin
!> @param[out] tauya   d^2(Rho)/dY^2, alpha-spin
!> @param[out] tauza   d^2(Rho)/dZ^2, alpha-spin
!> @param[out] tauxb   d^2(Rho)/dX^2, beta-spin
!> @param[out] tauyb   d^2(Rho)/dY^2, beta-spin
!> @param[out] tauzb   d^2(Rho)/dZ^2, beta-spin
!> @author Vladimir Mironov
 SUBROUTINE compD2RhoMO(self, iPt, tauxa, tauya, tauza, tauxb, tauyb, tauzb)

    CLASS(xc_engine_t) :: self
    INTEGER, INTENT(IN) :: iPt
    REAL(KIND=fp), INTENT(OUT) :: tauxa, tauya, tauza, tauxb, tauyb, tauzb

    ASSOCIATE ( noa => self%numOccAlpha, nob => self%numOccBeta)
    tauxa  = dot_product(self%moG1A(:noa,iPt,1), self%moG1A(:noa,iPt,1))
    tauya  = dot_product(self%moG1A(:noa,iPt,2), self%moG1A(:noa,iPt,2))
    tauza  = dot_product(self%moG1A(:noa,iPt,3), self%moG1A(:noa,iPt,3))
    IF (self%hasBeta) THEN
        tauxb  = dot_product(self%moG1B(:nob,iPt,1), self%moG1B(:nob,iPt,1))
        tauyb  = dot_product(self%moG1B(:nob,iPt,2), self%moG1B(:nob,iPt,2))
        tauzb  = dot_product(self%moG1B(:nob,iPt,3), self%moG1B(:nob,iPt,3))
    ELSE
        tauxb  = tauxa
        tauyb  = tauya
        tauzb  = tauza
    END IF
    END ASSOCIATE

 END SUBROUTINE

!> @brief Compute XC contributions to the gradient from a grid point, LDA part
!> @param[in]    iPt      index of a grid point
!> @param[inout] bfGrad   array of gradient contributinos per basis function
!> @param[in]    dedra    XC energy per unit volume, alpha-spin
!> @param[in]    dedrb    XC energy per unit volume, beta-spin
!> @author Vladimir Mironov
 SUBROUTINE compAtGradRho(self, iPt, bfGrad, dedra, dedrb)
    CLASS(xc_engine_t) :: self
    INTEGER, INTENT(IN) :: iPt
    REAL(KIND=fp), INTENT(IN) :: dedra, dedrb
    REAL(KIND=fp), CONTIGUOUS, INTENT(INOUT) :: bfGrad(:,:)

    REAL(KIND=fp) :: f

    ASSOCIATE ( aoG1  => self%aoG1 &
              , moVA => self%moVA &
              , moVB => self%moVB )

    f = dedra*2.0_fp

    bfGrad(:,1) = bfGrad(:,1) + (aoG1(:,iPt,1)*moVA(:,iPt))*f
    bfGrad(:,2) = bfGrad(:,2) + (aoG1(:,iPt,2)*moVA(:,iPt))*f
    bfGrad(:,3) = bfGrad(:,3) + (aoG1(:,iPt,3)*moVA(:,iPt))*f

    IF (.NOT.self%hasBeta) RETURN

    f = dedrb*2.0_fp

    bfGrad(:,1) = bfGrad(:,1) + (aoG1(:,iPt,1)*moVB(:,iPt))*f
    bfGrad(:,2) = bfGrad(:,2) + (aoG1(:,iPt,2)*moVB(:,iPt))*f
    bfGrad(:,3) = bfGrad(:,3) + (aoG1(:,iPt,3)*moVB(:,iPt))*f

    END ASSOCIATE
 END SUBROUTINE

!> @brief Compute XC contributions to the gradient from a grid point, GGA part
!> @param[in]    iPt      index of a grid point
!> @param[inout] bfGrad   array of gradient contributinos per basis function
!> @param[in]    dedsa    XC energy, GGA contribution, alpha-spin
!> @param[in]    dedsb    XC energy, GGA contribution, beta-spin
!> @author Vladimir Mironov
 SUBROUTINE compAtGradDRho(self, iPt, bfGrad, dedsa, dedsb)
    CLASS(xc_engine_t) :: self
    INTEGER, INTENT(IN) :: iPt
    REAL(KIND=fp), CONTIGUOUS, INTENT(IN) :: dedsa(:), dedsb(:)
    REAL(KIND=fp), CONTIGUOUS, INTENT(INOUT) :: bfGrad(:,:)

    REAL(KIND=fp) :: f(3)

    ASSOCIATE ( aoG1    => self%aoG1 &
              , moVA    => self%moVA &
              , moVB    => self%moVB &
              , moG1A   => self%moG1A &
              , moG1B   => self%moG1B &
              , aoG2xx  => self%aoG2(:,iPt,1) &
              , aoG2yy  => self%aoG2(:,iPt,2) &
              , aoG2zz  => self%aoG2(:,iPt,3) &
              , aoG2xy  => self%aoG2(:,iPt,4) &
              , aoG2yz  => self%aoG2(:,iPt,5) &
              , aoG2xz  => self%aoG2(:,iPt,6) &
              )

    f = dedsa(1:3)*2.0_fp

    bfGrad(:,1) = bfGrad(:,1) + &
            (aoG2xx(:)*moVA(:,iPt) + aoG1(:,iPt,1)*moG1A(:,iPt,1))*f(1) &
          + (aoG2xy(:)*moVA(:,iPt) + aoG1(:,iPt,1)*moG1A(:,iPt,2))*f(2) &
          + (aoG2xz(:)*moVA(:,iPt) + aoG1(:,iPt,1)*moG1A(:,iPt,3))*f(3)

    bfGrad(:,2) = bfGrad(:,2) + &
            (aoG2xy(:)*moVA(:,iPt) + aoG1(:,iPt,2)*moG1A(:,iPt,1))*f(1) &
          + (aoG2yy(:)*moVA(:,iPt) + aoG1(:,iPt,2)*moG1A(:,iPt,2))*f(2) &
          + (aoG2yz(:)*moVA(:,iPt) + aoG1(:,iPt,2)*moG1A(:,iPt,3))*f(3)

    bfGrad(:,3) = bfGrad(:,3) + &
            (aoG2xz(:)*moVA(:,iPt) + aoG1(:,iPt,3)*moG1A(:,iPt,1))*f(1) &
          + (aoG2yz(:)*moVA(:,iPt) + aoG1(:,iPt,3)*moG1A(:,iPt,2))*f(2) &
          + (aoG2zz(:)*moVA(:,iPt) + aoG1(:,iPt,3)*moG1A(:,iPt,3))*f(3)

    IF (.NOT.self%hasBeta) RETURN

    f = dedsb(1:3)*2.0_fp

    bfGrad(:,1) = bfGrad(:,1) + &
            (aoG2xx(:)*moVB(:,iPt) + aoG1(:,iPt,1)*moG1B(:,iPt,1))*f(1) &
          + (aoG2xy(:)*moVB(:,iPt) + aoG1(:,iPt,1)*moG1B(:,iPt,2))*f(2) &
          + (aoG2xz(:)*moVB(:,iPt) + aoG1(:,iPt,1)*moG1B(:,iPt,3))*f(3)

    bfGrad(:,2) = bfGrad(:,2) + &
            (aoG2xy(:)*moVB(:,iPt) + aoG1(:,iPt,2)*moG1B(:,iPt,1))*f(1) &
          + (aoG2yy(:)*moVB(:,iPt) + aoG1(:,iPt,2)*moG1B(:,iPt,2))*f(2) &
          + (aoG2yz(:)*moVB(:,iPt) + aoG1(:,iPt,2)*moG1B(:,iPt,3))*f(3)

    bfGrad(:,3) = bfGrad(:,3) + &
            (aoG2xz(:)*moVB(:,iPt) + aoG1(:,iPt,3)*moG1B(:,iPt,1))*f(1) &
          + (aoG2yz(:)*moVB(:,iPt) + aoG1(:,iPt,3)*moG1B(:,iPt,2))*f(2) &
          + (aoG2zz(:)*moVB(:,iPt) + aoG1(:,iPt,3)*moG1B(:,iPt,3))*f(3)


    END ASSOCIATE
 END SUBROUTINE


!> @brief Compute XC contributions to the gradient from a grid point, mGGA part
!> @param[in]    iPt      index of a grid point
!> @param[inout] bfGrad   array of gradient contributinos per basis function
!> @param[in]    dedta    XC energy, mGGA contribution, alpha-spin
!> @param[in]    dedtb    XC energy, mGGA contribution, beta-spin
!> @author Vladimir Mironov
 SUBROUTINE compAtGradTau(self, iPt, bfGrad, dedta, dedtb)
    CLASS(xc_engine_t) :: self
    INTEGER, INTENT(IN) :: iPt
    REAL(KIND=fp), INTENT(IN) :: dedta, dedtb
    REAL(KIND=fp), CONTIGUOUS, INTENT(INOUT) :: bfGrad(:,:)

    REAL(KIND=fp) :: f

    ASSOCIATE ( aoG1    => self%aoG1 &
              , moVA    => self%moVA &
              , moVB    => self%moVB &
              , moG1A   => self%moG1A &
              , moG1B   => self%moG1B &
              , aoG2xx  => self%aoG2(:,iPt,1) &
              , aoG2yy  => self%aoG2(:,iPt,2) &
              , aoG2zz  => self%aoG2(:,iPt,3) &
              , aoG2xy  => self%aoG2(:,iPt,4) &
              , aoG2yz  => self%aoG2(:,iPt,5) &
              , aoG2xz  => self%aoG2(:,iPt,6) &
              )

    f = dedta*2.0_fp

    bfGrad(:,1) = bfGrad(:,1) + &
            (  aoG2xx(:)*moG1A(:,iPt,1) &
             + aoG2xy(:)*moG1A(:,iPt,2) &
             + aoG2xz(:)*moG1A(:,iPt,3) )*f

    bfGrad(:,2) = bfGrad(:,2) + &
            (  aoG2xy(:)*moG1A(:,iPt,1) &
             + aoG2yy(:)*moG1A(:,iPt,2) &
             + aoG2yz(:)*moG1A(:,iPt,3) )*f

    bfGrad(:,3) = bfGrad(:,3) + &
            (  aoG2xz(:)*moG1A(:,iPt,1) &
             + aoG2yz(:)*moG1A(:,iPt,2) &
             + aoG2zz(:)*moG1A(:,iPt,3) )*f


    IF (.NOT.self%hasBeta) RETURN

    f = dedtb*2.0_fp

    bfGrad(:,1) = bfGrad(:,1) + &
            (  aoG2xx(:)*moG1B(:,iPt,1) &
             + aoG2xy(:)*moG1B(:,iPt,2) &
             + aoG2xz(:)*moG1B(:,iPt,3) )*f

    bfGrad(:,2) = bfGrad(:,2) + &
            (  aoG2xy(:)*moG1B(:,iPt,1) &
             + aoG2yy(:)*moG1B(:,iPt,2) &
             + aoG2yz(:)*moG1B(:,iPt,3) )*f

    bfGrad(:,3) = bfGrad(:,3) + &
            (  aoG2xz(:)*moG1B(:,iPt,1) &
             + aoG2yz(:)*moG1B(:,iPt,2) &
             + aoG2zz(:)*moG1B(:,iPt,3) )*f

    END ASSOCIATE
 END SUBROUTINE

!> @brief Add XC contribution to the Kohn-Sham matrix
!> @param[inout] fa2   Kohn-Sham matrix, full format, alpha-spin
!> @param[inout] fb2   Kohn-Sham matrix, full format, beta-spin
!> @author Vladimir Mironov
 SUBROUTINE updKS(self,fa2,fb2)
    CLASS(xc_engine_t) :: self
    REAL(KIND=fp), CONTIGUOUS, INTENT(INOUT) :: fa2(:), fb2(:)

    INTEGER :: i, j, k

    IF (.NOT.self%isGGA) THEN

        DO i = 1, self%numPts
            self%tmp(:,i) = 0.5 * self%dedr(i,1) &
                                * self%aoV(:,i)! &
                             !* self%wts(i)
        END DO
        CALL dsyr2k('U', 'N', self%numAOs, self%numPts, 1.0_fp, &
                     self%aoV(1,1),  self%numAOs, &
                     self%tmp(1,1),  self%numAOs, &
                     1.0_fp, fa2(1), self%numAOs)

        IF (self%hasBeta) THEN

            DO i = 1, self%numPts
                self%tmp(:,i) = 0.5 * self%dedr(i,2) &
                                    * self%aoV(:,i)! &
                                 !* self%wts(i)
            END DO
            CALL dsyr2k('U', 'N', self%numAOs, self%numPts, 1.0_fp, &
                         self%aoV(1,1),  self%numAOs, &
                         self%tmp(1,1),  self%numAOs, &
                         1.0_fp, fb2(1), self%numAOs)
        END IF

    ELSE !IF (self%isGGA) THEN
        DO i = 1, self%numPts
            self%tmp(:,i) = &
                      0.5_fp*self%dedr(i,1)*self%aoV(:,i) &
                           + self%deds(i,1)*self%aoG1(:,i,1) &
                           + self%deds(i,2)*self%aoG1(:,i,2) &
                           + self%deds(i,3)*self%aoG1(:,i,3)
            !self%tmp(:,i) = self%tmp(:,i) * self%wts(i)
        END DO

        CALL dsyr2k('U', 'N', self%numAOs, self%numPts, 1.0_fp, &
                     self%aoV(1,1),  self%numAOs, &
                     self%tmp(1,1),  self%numAOs, &
                     1.0_fp, fa2(1), self%numAOs)

        IF (self%hasBeta) THEN

            DO i = 1, self%numPts
                self%tmp(:,i) = &
                          0.5_fp*self%dedr(i,2)*self%aoV(:,i) &
                               + self%deds(i,4)*self%aoG1(:,i,1) &
                               + self%deds(i,5)*self%aoG1(:,i,2) &
                               + self%deds(i,6)*self%aoG1(:,i,3)
                !self%tmp(:,i) = self%tmp(:,i) * self%wts(i)
            END DO

            CALL dsyr2k('U', 'N', self%numAOs, self%numPts, 1.0_fp, &
                         self%aoV(1,1),  self%numAOs, &
                         self%tmp(1,1),  self%numAOs, &
                         1.0_fp, fb2(1), self%numAOs)
        END IF

    END IF
    IF (self%needTau) THEN
        DO i = 1, self%numPts
            self%tmp(:,i              ) = &
                    self%dedt(i,1)*self%aoG1(:,i,1)!*self%wts(i)
            self%tmp(:,i+  self%numPts) = &
                    self%dedt(i,1)*self%aoG1(:,i,2)!*self%wts(i)
            self%tmp(:,i+2*self%numPts) = &
                    self%dedt(i,1)*self%aoG1(:,i,3)!*self%wts(i)
        END DO

        IF (ubound(self%aoG1,2)==self%numPts) THEN
            CALL dsyr2k('U', 'N', self%numAOs, 3*self%numPts, 0.5_fp, &
                         self%aoG1(1,1,1), self%numAOs, &
                         self%tmp(1,1),    self%numAOs, &
                         1.0_fp, fa2(1),   self%numAOs)
         ELSE
            DO j = 1, 3
                k = (j-1)*self%numPts+1
                CALL dsyr2k('U', 'N', self%numAOs, self%numPts, 0.5_fp, &
                    self%aoG1(1,1,j), self%numAOs, &
                    self%tmp(1,k),    self%numAOs, &
                    1.0_fp, fa2(1),   self%numAOs)
            END DO
        END IF

        IF (self%hasBeta) THEN
            DO i = 1, self%numPts
                self%tmp(:,i              ) = &
                        self%dedt(i,2)*self%aoG1(:,i,1)!*self%wts(i)
                self%tmp(:,i+  self%numPts) = &
                        self%dedt(i,2)*self%aoG1(:,i,2)!*self%wts(i)
                self%tmp(:,i+2*self%numPts) = &
                        self%dedt(i,2)*self%aoG1(:,i,3)!*self%wts(i)
            END DO

            IF (ubound(self%aoG1,2)==self%numPts) THEN
                CALL dsyr2k('U', 'N', self%numAOs, 3*self%numPts, 0.5_fp, &
                             self%aoG1(1,1,1), self%numAOs, &
                             self%tmp(1,1),    self%numAOs, &
                             1.0_fp, fb2(1),   self%numAOs)
             ELSE
                DO j = 1, 3
                    k = (j-1)*self%numPts+1
                    CALL dsyr2k('U', 'N', self%numAOs, self%numPts, 0.5_fp, &
                                 self%aoG1(1,1,j), self%numAOs, &
                                 self%tmp(1,k),    self%numAOs, &
                                 1.0_fp, fb2(1),   self%numAOs)
                END DO
            END IF
        END IF
    END IF
 END SUBROUTINE

! SUBROUTINE updGrad(self,gradSh,n,da2,db2)
!    CLASS(xc_engine_t) :: self
!    REAL(KIND=fp) :: da2(:),db2(:), gradSh(:)
!    INTEGER :: n
!    INTEGER :: i
!
!!   Compute sum(j) d(i,j) f(j)
!!           sum(j) d(i,j) g[k](j)
!    DRa (1:num) = matmul(Da(1:num,1:num)* AOx(1:num))
!    DRxa(1:num) = matmul(Da(1:num,1:num)*gAOx(1:num))
!    DRya(1:num) = matmul(Da(1:num,1:num)*gAOy(1:num))
!    DRza(1:num) = matmul(Da(1:num,1:num)*gAOz(1:num))
!
!
!    gradx_lda = gaox(:)*DRa(:)*duma
!    grady_lda = gaoy(:)*DRa(:)*duma
!    gradz_lda = gaoz(:)*DRa(:)*duma
!
!
!    gradx_gga = (g2aoxx(:)*DRa(:)+gaox(:)*DRxa(:))*dumax + &
!                (g2aoxy(:)*DRa(:)+gaox(:)*DRya(:))*dumay + &
!                (g2aoxz(:)*DRa(:)+gaox(:)*DRza(:))*dumaz
!
!    grady_gga = (g2aoxy(:)*DRa(:)+gaoy(:)*DRxa(:))*dumax + &
!                (g2aoyy(:)*DRa(:)+gaoy(:)*DRya(:))*dumay + &
!                (g2aoyz(:)*DRa(:)+gaoy(:)*DRza(:))*dumaz
!
!    gradz_gga = (g2aoxz(:)*DRa(:)+gaoz(:)*DRxa(:))*dumax + &
!                (g2aoyz(:)*DRa(:)+gaoz(:)*DRya(:))*dumay + &
!                (g2aozz(:)*DRa(:)+gaoz(:)*DRza(:))*dumaz
!
!
!
!    gradx_mgga = (g2aoxx(:)*DRxa(:) + &
!                  g2aoxy(:)*DRya(:) + &
!                  g2aoxz(:)*DRza(:))*dmga
!
!    grady_mgga = (g2aoxy(:)*DRxa(:) + &
!                  g2aoyy(:)*DRya(:) + &
!                  g2aoyz(:)*DRza(:))*dmga
!
!    gradz_mgga = (g2aoxz(:)*DRxa(:) + &
!                  g2aoyz(:)*DRya(:) + &
!                  g2aozz(:)*DRza(:))*dmga
!
!
!
!    dedftsh(1,:) = dedftsh(1,:) - TWO*ftotwt*(gradx_lda+gradx_gga+gradx_mgga)
!    dedftsh(2,:) = dedftsh(2,:) - TWO*ftotwt*(grady_lda+grady_gga+grady_mgga)
!    dedftsh(3,:) = dedftsh(3,:) - TWO*ftotwt*(gradz_lda+gradz_gga+gradz_mgga)
!
!    dedft(1,n) = dedft(1,n) + TWO*ftotwt*sum(gradx_lda+gradx_gga+gradx_mgga)
!    dedft(2,n) = dedft(2,n) + TWO*ftotwt*sum(grady_lda+grady_gga+grady_mgga)
!    dedft(3,n) = dedft(3,n) + TWO*ftotwt*sum(gradz_lda+gradz_gga+gradz_mgga)
!
! END SUBROUTINE

!> @brief Compute AO values in a point
!> @param[in]    basis  atomic basis set
!> @param[in]    tol    threshold value
!> @param[in]    dr     array of (atom[xyz]-point[xyz])**(0..maxAng)
!> @param[in]    rsqrd  array of square distances between point and atoms
!> @param[out]   naos   number of significant AOs
!> @param[out]   aov    AO values
!> @author Vladimir Mironov
 SUBROUTINE compAOv(basis,tol,dr,rsqrd,naos, &
                 aov)

    USE prec, ONLY: fp
    USE mx_limits, ONLY: mxsh, mxgtot
    USE mod_nosp_basis, ONLY: basis_set
    IMPLICIT NONE

    TYPE(basis_set) :: basis

    REAL(KIND=fp), CONTIGUOUS, INTENT(OUT) :: aov(:)
    INTEGER, INTENT(OUT) :: naos

    REAL(KIND=fp), CONTIGUOUS, INTENT(IN) :: dr(:,:,:) ! dim: (nat,*,3)
    REAL(KIND=fp), CONTIGUOUS, INTENT(IN) :: rsqrd(:) !dim: (nat)

    REAL(KIND=fp), INTENT(IN) :: tol

    REAL(KIND=fp) :: &
        vexp, vexp1, vexp2, dum
    INTEGER :: &
        k2, loci, &
        ishell, imomfct, ix, iy, iz, ityp

    naos = 0
    DO ishell = 1, basis%nshell
        ASSOCIATE ( &
            iatm  =>  basis%katom(ishell), &
            kloc  =>   basis%kloc(ishell), &
            k1    => basis%kstart(ishell), &
            kng   =>    basis%kng(ishell), &
            ktype =>  basis%ktype(ishell), &
            mini  =>   basis%kmin(ishell), &
            maxi  =>   basis%kmax(ishell) )
!       Angular intermediates for the density and gradient

        vexp1 = 0.0_fp
        vexp2 = 0.0_fp

        k2 = kng - 1 + k1
        DO imomfct = k1, k2
            dum = basis%ex(imomfct)*rsqrd(iatm)
            IF (dum>tol) CYCLE
            vexp = exp(-dum)*basis%cc(imomfct)
            vexp1 = vexp1 + vexp
            vexp2 = vexp2 + 2*vexp*basis%ex(imomfct)
        END DO

        IF (vexp1==0.0_fp) THEN
            aov(kloc:kloc+maxi-mini) = 0.0_fp
            CYCLE
        END IF
        naos = naos + 1

        SELECT CASE (ktype)
        CASE (1)
            aov(kloc) = vexp1
        CASE (2)
            aov(kloc  ) = vexp1*dr(3,1,iatm)
            aov(kloc+1) = vexp1*dr(3,2,iatm)
            aov(kloc+2) = vexp1*dr(3,3,iatm)
        CASE DEFAULT
            loci  = kloc - mini
            DO ityp = mini, maxi
                ix = ijx(ityp)
                iy = ijy(ityp)
                iz = ijz(ityp)

!               Compute AO value at a grid point
                aov(loci+ityp) = vexp1 * dr(ix+1,1,iatm) &
                                       * dr(iy+1,2,iatm) &
                                       * dr(iz+1,3,iatm)
             END DO
         END SELECT
         END ASSOCIATE

    END DO

 END SUBROUTINE

!> @brief Compute AO values and gradient in a point
!> @param[in]    basis  atomic basis set
!> @param[in]    tol    threshold value (not used now)
!> @param[in]    dr     array of (atom[xyz]-point[xyz])**(0..maxAng)
!> @param[in]    rsqrd  array of square distances between point and atoms
!> @param[out]   naos   number of significant AOs
!> @param[out]   aov    AO values
!> @param[out]   aogx   AO gradient, X component
!> @param[out]   aogy   AO gradient, Y component
!> @param[out]   aogz   AO gradient, Z component
!> @author Vladimir Mironov
 SUBROUTINE compAOvg(basis,tol,dr,rsqrd,naos, &
                 aov,aogx,aogy,aogz)

    USE prec, ONLY: fp
    USE mx_limits, ONLY: mxsh, mxgtot
    USE mod_nosp_basis, ONLY: basis_set
    USE mod_dft_fuzzycell, ONLY: prim_mx_dist2, shell_mx_dist2
    IMPLICIT NONE

    TYPE(basis_set) :: basis

    INTEGER, INTENT(OUT) :: naos
    REAL(KIND=fp), CONTIGUOUS, INTENT(IN) :: rsqrd(:) !dim: (nat)
    REAL(KIND=fp), INTENT(IN) :: tol

    REAL(KIND=fp), CONTIGUOUS, INTENT(OUT) :: aov(:)
    REAL(KIND=fp), CONTIGUOUS, INTENT(OUT) :: aogx(:),aogy(:),aogz(:)
    REAL(KIND=fp), CONTIGUOUS, INTENT(IN)  :: dr(:,:,:)

    INTEGER :: &
        ishell, loci, &!, iatm, mini, maxi, loci0, & !ifct
        ityp, ix, iy, iz

    REAL(KIND=fp) :: &
        vexp1, vexp2, xy, xz, yz, &
        xm, ym, zm, &
        xp, yp, zp

    REAL(KIND=fp) :: &
        dum, vexp
    INTEGER :: &
        k2, imomfct

    naos = 0

    DO ishell = 1, basis%nshell
        ASSOCIATE ( &
            iatm  => basis%katom(ishell), &
            mini  =>  basis%kmin(ishell), &
            maxi  =>  basis%kmax(ishell), &
            kloc  =>  basis%kloc(ishell), &
            k1    => basis%kstart(ishell), &
            kng   =>    basis%kng(ishell), &
            ktype => basis%ktype(ishell) )

        vexp1 = 0.0_fp
        vexp2 = 0.0_fp

        IF (rsqrd(iatm)<=shell_mx_dist2(ishell)) THEN
            k2 = kng - 1 + k1
            DO imomfct = k1, k2
                IF (rsqrd(iatm)>prim_mx_dist2(imomfct)) CYCLE
                dum = basis%ex(imomfct)*rsqrd(iatm)
                !IF (dum>tol) CYCLE
                vexp = exp(-dum)*basis%cc(imomfct)
                vexp1 = vexp1 + vexp
                vexp2 = vexp2 + 2*vexp*basis%ex(imomfct)
            END DO
!        END IF
!        IF (vexp1==0.0_fp) THEN
        ELSE
            aov(kloc:kloc+maxi-mini)  = 0.0_fp
            aogx(kloc:kloc+maxi-mini) = 0.0_fp
            aogy(kloc:kloc+maxi-mini) = 0.0_fp
            aogz(kloc:kloc+maxi-mini) = 0.0_fp
            CYCLE
        END IF
        naos = naos + 1

        loci  = kloc-mini
        SELECT CASE (ktype)
        CASE (1)
!           Special fast code for S functions
            aov(kloc)  = vexp1
            aogx(kloc) = -vexp2*dr(3,1,iatm)
            aogy(kloc) = -vexp2*dr(3,2,iatm)
            aogz(kloc) = -vexp2*dr(3,3,iatm)
        CASE(2)
!           Special fast code for P functions
            xy = vexp2*dr(3,1,iatm)*dr(3,2,iatm)
            xz = vexp2*dr(3,1,iatm)*dr(3,3,iatm)
            yz = vexp2*dr(3,2,iatm)*dr(3,3,iatm)

            aov(kloc  ) = vexp1*dr(3,1,iatm)
            aov(kloc+1) = vexp1*dr(3,2,iatm)
            aov(kloc+2) = vexp1*dr(3,3,iatm)

            aogx(kloc  ) = vexp1 - vexp2*dr(4,1,iatm)
            aogx(kloc+1) = -xy
            aogx(kloc+2) = -xz

            aogy(kloc  ) = -xy
            aogy(kloc+1) = vexp1 - vexp2*dr(4,2,iatm)
            aogy(kloc+2) = -yz

            aogz(kloc  ) = -xz
            aogz(kloc+1) = -yz
            aogz(kloc+2) = vexp1 - vexp2*dr(4,3,iatm)
        CASE DEFAULT
            DO ityp = mini, maxi
                ix = ijx(ityp)
                iy = ijy(ityp)
                iz = ijz(ityp)

                aov(loci+ityp) = vexp1 * dr(ix+1,1,iatm) &
                                       * dr(iy+1,2,iatm) &
                                       * dr(iz+1,3,iatm)

!               Compute gradient AO value at a grid point
!               Gradient is by the electron (not nuclear) coordinates
                xy = dr(ix+1,1,iatm)*dr(iy+1,2,iatm)
                xz = dr(ix+1,1,iatm)*dr(iz+1,3,iatm)
                yz = dr(iy+1,2,iatm)*dr(iz+1,3,iatm)

!               Gradient minus one component
                xm = (ix-1)*dr(ix,1,iatm)
                ym = (iy-1)*dr(iy,2,iatm)
                zm = (iz-1)*dr(iz,3,iatm)

!               Gradient plus one component
                xp = dr(ix+2,1,iatm)*vexp2
                yp = dr(iy+2,2,iatm)*vexp2
                zp = dr(iz+2,3,iatm)*vexp2

                aogx(loci+ityp) =  ( -xp + vexp1*xm ) * yz
                aogy(loci+ityp) =  ( -yp + vexp1*ym ) * xz
                aogz(loci+ityp) =  ( -zp + vexp1*zm ) * xy
            END DO
        END SELECT

        END ASSOCIATE

    END DO

 END SUBROUTINE

!> @brief Compute AO values, and their 1st and 2nd derivatives in a point
!> @param[in]    basis    atomic basis set
!> @param[in]    tol      threshold value (not used now)
!> @param[in]    dr       array of (atom[xyz]-point[xyz])**(0..maxAng)
!> @param[in]    rsqrd    array of square distances between point and atoms
!> @param[out]   naos     number of significant AOs
!> @param[out]   aov      AO values
!> @param[out]   aogx     AO gradient, X component
!> @param[out]   aogy     AO gradient, Y component
!> @param[out]   aogz     AO gradient, Z component
!> @param[out]   aog2xx   AO 2nd derivative, XX component
!> @param[out]   aog2yy   AO 2nd derivative, YY component
!> @param[out]   aog2zz   AO 2nd derivative, ZZ component
!> @param[out]   aog2xy   AO 2nd derivative, XY component
!> @param[out]   aog2yz   AO 2nd derivative, YZ component
!> @param[out]   aog2xz   AO 2nd derivative, XZ component
!> @author Vladimir Mironov
 SUBROUTINE compAOvgg(basis,tol,dr,rsqrd,naos,&
                 aov,aogx,aogy,aogz, &
                 aog2xx,aog2yy,aog2zz,aog2xy,aog2yz,aog2xz)

    USE prec, ONLY: fp
    USE mx_limits, ONLY: mxsh, mxgtot
    USE mod_nosp_basis, ONLY: basis_set
    USE mod_dft_fuzzycell, ONLY: prim_mx_dist2, shell_mx_dist2
    IMPLICIT NONE

    TYPE(basis_set) :: basis

    INTEGER, INTENT(OUT) :: naos
    REAL(KIND=fp), CONTIGUOUS, INTENT(IN) :: rsqrd(:) !dim: (nat)
    REAL(KIND=fp), INTENT(IN) :: tol

    REAL(KIND=fp), CONTIGUOUS, INTENT(OUT) :: aov(:)
    REAL(KIND=fp), CONTIGUOUS, INTENT(OUT) :: aogx(:),aogy(:),aogz(:)
    REAL(KIND=fp), CONTIGUOUS, INTENT(OUT) :: &
            aog2xx(:), aog2yy(:), aog2zz(:), &
            aog2xy(:), aog2yz(:), aog2xz(:)
    REAL(KIND=fp), CONTIGUOUS, INTENT(IN)  :: dr(:,:,:)

    INTEGER :: &
        ishell, loci, &!, iatm, mini, maxi, loci0, & !ifct
        ityp, ix, iy, iz
    INTEGER :: &
        ixm, iym, izm

    REAL(KIND=fp) :: &
        vexp1, vexp2, vexp3, xy, xz, yz, xyz, &
        xm, ym, zm, &
        xp, yp, zp, &
        x, y, z

    REAL(KIND=fp) :: &
        xmxm, ymym, zmzm, &
        xpxp, ypyp, zpzp, &
        xmym, xmzm, ymzm, &
        xpyp, xpzp, ypzp, &
        xpym, xmyp, xpzm, xmzp, ypzm, ymzp

    REAL(KIND=fp) :: &
        dum, vexp, tmp(3)
    INTEGER :: &
        k2, imomfct

    naos = 0

    DO ishell = 1, basis%nshell
      ASSOCIATE ( &
          iatm  => basis%katom(ishell), &
          mini  =>  basis%kmin(ishell), &
          maxi  =>  basis%kmax(ishell), &
          kloc  =>  basis%kloc(ishell), &
          k1    => basis%kstart(ishell), &
          kng   =>    basis%kng(ishell), &
          ktype => basis%ktype(ishell) )

      vexp1 = 0.0_fp
      vexp2 = 0.0_fp
      vexp3 = 0.0_fp

      IF (rsqrd(iatm)<=shell_mx_dist2(ishell)) THEN
          k2 = kng - 1 + k1
          DO imomfct = k1, k2
              IF (rsqrd(iatm)>prim_mx_dist2(imomfct)) CYCLE
              dum = basis%ex(imomfct)*rsqrd(iatm)
              !IF (dum>tol) CYCLE
              vexp = exp(-dum)*basis%cc(imomfct)
              vexp1 = vexp1 + vexp
              vexp2 = vexp2 + 2*vexp*basis%ex(imomfct)
              vexp3 = vexp3 + 4*vexp*basis%ex(imomfct)*basis%ex(imomfct)
          END DO
!      END IF
!      IF (vexp1==0.0_fp) THEN
      ELSE
        aov(kloc:kloc+maxi-mini)    = 0.0_fp
        aogx(kloc:kloc+maxi-mini)   = 0.0_fp
        aogy(kloc:kloc+maxi-mini)   = 0.0_fp
        aogz(kloc:kloc+maxi-mini)   = 0.0_fp
        aog2xx(kloc:kloc+maxi-mini) = 0.0_fp
        aog2yy(kloc:kloc+maxi-mini) = 0.0_fp
        aog2zz(kloc:kloc+maxi-mini) = 0.0_fp
        aog2xy(kloc:kloc+maxi-mini) = 0.0_fp
        aog2yz(kloc:kloc+maxi-mini) = 0.0_fp
        aog2xz(kloc:kloc+maxi-mini) = 0.0_fp
        CYCLE
      END IF
      naos = naos + 1

      loci  = kloc-mini
      SELECT CASE (ktype)
      CASE (1)
!       Special fast code for S functions
        aov(kloc)  =  vexp1

        aogx(kloc) = -vexp2*dr(3,1,iatm)
        aogy(kloc) = -vexp2*dr(3,2,iatm)
        aogz(kloc) = -vexp2*dr(3,3,iatm)

        xy = dr(3,1,iatm)*dr(3,2,iatm)
        xz = dr(3,1,iatm)*dr(3,3,iatm)
        yz = dr(3,2,iatm)*dr(3,3,iatm)

        aog2xx(kloc) = vexp3*dr(4,1,iatm) - vexp2
        aog2yy(kloc) = vexp3*dr(4,2,iatm) - vexp2
        aog2zz(kloc) = vexp3*dr(4,3,iatm) - vexp2
        aog2xy(kloc) = vexp3*xy
        aog2yz(kloc) = vexp3*yz
        aog2xz(kloc) = vexp3*xz

      CASE(2)
!       Special fast code for P functions
        xy  = dr(3,1,iatm)*dr(3,2,iatm)
        xz  = dr(3,1,iatm)*dr(3,3,iatm)
        yz  = dr(3,2,iatm)*dr(3,3,iatm)
        xyz = dr(3,1,iatm)*yz

        aov(kloc  ) = vexp1*dr(3,1,iatm)
        aov(kloc+1) = vexp1*dr(3,2,iatm)
        aov(kloc+2) = vexp1*dr(3,3,iatm)

        aogx(kloc  ) = vexp1 - vexp2*dr(4,1,iatm)
        aogx(kloc+1) =       - vexp2*xy
        aogx(kloc+2) =       - vexp2*xz

        aogy(kloc  ) =       - vexp2*xy
        aogy(kloc+1) = vexp1 - vexp2*dr(4,2,iatm)
        aogy(kloc+2) =       - vexp2*yz

        aogz(kloc  ) =       - vexp2*xz
        aogz(kloc+1) =       - vexp2*yz
        aogz(kloc+2) = vexp1 - vexp2*dr(4,3,iatm)

        tmp = vexp3*dr(4,1:3,iatm) - vexp2

        aog2xx(kloc  ) = (tmp(1) - 2.0_fp*vexp2)*dr(3,1,iatm)
        aog2xx(kloc+1) =  tmp(1)                *dr(3,2,iatm)
        aog2xx(kloc+2) =  tmp(1)                *dr(3,3,iatm)

        aog2yy(kloc  ) =  tmp(2)                *dr(3,1,iatm)
        aog2yy(kloc+1) = (tmp(2) - 2.0_fp*vexp2)*dr(3,2,iatm)
        aog2yy(kloc+2) =  tmp(2)                *dr(3,3,iatm)

        aog2zz(kloc  ) =  tmp(3)                *dr(3,1,iatm)
        aog2zz(kloc+1) =  tmp(3)                *dr(3,2,iatm)
        aog2zz(kloc+2) = (tmp(3) - 2.0_fp*vexp2)*dr(3,3,iatm)

        aog2xy(kloc  ) =       aog2xx(kloc+1)
        aog2xy(kloc+1) =       aog2yy(kloc  )
        aog2xy(kloc+2) = vexp3*xyz

        aog2yz(kloc  ) = vexp3*xyz
        aog2yz(kloc+1) =       aog2yy(kloc+2)
        aog2yz(kloc+2) =       aog2zz(kloc+1)

        aog2xz(kloc  ) =       aog2xx(kloc+2)
        aog2xz(kloc+1) = vexp3*xyz
        aog2xz(kloc+2) =       aog2zz(kloc  )

      CASE DEFAULT
        DO ityp = mini, maxi
          ix = ijx(ityp)
          iy = ijy(ityp)
          iz = ijz(ityp)

!         0 components:
          x = dr(ix+1,1,iatm)
          y = dr(iy+1,2,iatm)
          z = dr(iz+1,3,iatm)

          xz = x*z
          xy = x*y
          yz = y*z
          xyz = x*y*z

!         +1 components
          xp = dr(ix+2,1,iatm)
          yp = dr(iy+2,2,iatm)
          zp = dr(iz+2,3,iatm)

!         -1 components
          xm = (ix-1)*dr(ix,1,iatm)
          ym = (iy-1)*dr(iy,2,iatm)
          zm = (iz-1)*dr(iz,3,iatm)

!         +2 components
          xpxp = dr(ix+3,1,iatm)
          ypyp = dr(iy+3,2,iatm)
          zpzp = dr(iz+3,3,iatm)

          xpyp = dr(ix+2,1,iatm)*dr(iy+2,2,iatm)
          xpzp = dr(ix+2,1,iatm)*dr(iz+2,3,iatm)
          ypzp = dr(iy+2,2,iatm)*dr(iz+2,3,iatm)

!         -2 components
          ixm = max(ix-1,1)
          iym = max(iy-1,1)
          izm = max(iz-1,1)
          xmxm = (ix-1)*(ix-2)*dr(ixm,1,iatm)
          ymym = (iy-1)*(iy-2)*dr(iym,2,iatm)
          zmzm = (iz-1)*(iz-2)*dr(izm,3,iatm)

          xmym = (ix-1)*(iy-1)*dr(ix,1,iatm)*dr(iy,2,iatm)
          xmzm = (ix-1)*(iz-1)*dr(ix,1,iatm)*dr(iz,3,iatm)
          ymzm = (iy-1)*(iz-1)*dr(iy,2,iatm)*dr(iz,3,iatm)

!         +1/-1  components
          xpym = (iy-1)*dr(ix+2,1,iatm)*dr(iy  ,2,iatm)
          xmyp = (ix-1)*dr(ix  ,1,iatm)*dr(iy+2,2,iatm)

          xpzm = (iz-1)*dr(ix+2,1,iatm)*dr(iz  ,3,iatm)
          xmzp = (ix-1)*dr(ix  ,1,iatm)*dr(iz+2,3,iatm)

          ypzm = (iz-1)*dr(iy+2,2,iatm)*dr(iz  ,3,iatm)
          ymzp = (iy-1)*dr(iy  ,2,iatm)*dr(iz+2,3,iatm)

!         AO value:
          aov(loci+ityp) = vexp1 * xyz

!         First AO derivatives:
          aogx(loci+ityp) =  (-vexp2*xp + vexp1*xm) * yz
          aogy(loci+ityp) =  (-vexp2*yp + vexp1*ym) * xz
          aogz(loci+ityp) =  (-vexp2*zp + vexp1*zm) * xy

!         Second AO derivatives:
          aog2xx(loci+ityp) = (vexp3*xpxp - vexp2*(2*ix-1)*x + vexp1*xmxm)*yz
          aog2yy(loci+ityp) = (vexp3*ypyp - vexp2*(2*iy-1)*y + vexp1*ymym)*xz
          aog2zz(loci+ityp) = (vexp3*zpzp - vexp2*(2*iz-1)*z + vexp1*zmzm)*xy

          aog2xy(loci+ityp) = (vexp3*xpyp - vexp2*(xpym+xmyp) + vexp1*xmym)*z
          aog2yz(loci+ityp) = (vexp3*ypzp - vexp2*(ypzm+ymzp) + vexp1*ymzm)*x
          aog2xz(loci+ityp) = (vexp3*xpzp - vexp2*(xpzm+xmzp) + vexp1*xmzm)*y

        END DO
      END SELECT

      END ASSOCIATE

    END DO

 END SUBROUTINE

!> @brief Original algorithm of KS matrix update
!>  Reserved for further use
!> @author Vladimir Mironov
 SUBROUTINE dftfock_sparse(needgr,needtau,urohf,&
                     ftotwt,duma,dumb,dumax,dumay, &
                     dumaz,dumbx,dumby,dumbz,valg,aox,gaox, &
                     gaoy,gaoz,fa,fb,cutoff,l1,dmgga,dmggb)
    USE prec, ONLY: fp
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: urohf, needgr, needtau
    REAL(KIND=fp) :: aox(*),gaox(*),gaoy(*),gaoz(*),fa(*),fb(*),valg(*)

    REAL(KIND=fp) :: ftotwt, cutoff, dmgga, dmggb, &
            duma, dumax, dumay, dumaz, &
            dumb, dumbx, dumby, dumbz

    INTEGER :: l1

    INTEGER :: ind, i
    REAL(KIND=fp) :: ffa, ffb, aomax, amaxval, bmaxval, c1, c2

    IF (.NOT.needgr) THEN

        DO i = 1, l1
            ind = (i-1)*l1
            ffa = ftotwt*duma*aox(i)
            IF (abs(ffa)>cutoff) &
                    fa(ind+1:ind+i) = fa(ind+1:ind+i) + ffa*aox(1:i)
        END DO

        IF (.NOT.urohf) RETURN
        DO i = 1, l1
            ind = (i-1)*l1
            ffb = ftotwt*dumb*aox(i)
            IF (abs(ffb)>cutoff) &
                    fb(ind+1:ind+i) = fb(ind+1:ind+i) + ffb*aox(1:i)
        END DO

    ELSE IF (.NOT.needtau) THEN ! old GGA code
        aomax   = 0.0
        amaxval = 0.0
        bmaxval = 0.0
        DO i = 1, l1
            ind = (i-1)*l1
            c1 = aox(i)
            c2 = ftotwt*(0.5_fp*duma * aox(i) &
                              + dumax*gaox(i) &
                              + dumay*gaoy(i) &
                              + dumaz*gaoz(i))
            valg(i) = c2
            aomax   = max(aomax,   abs(c1))
            amaxval = max(amaxval, abs(c2))
            IF (abs(c1*amaxval) + abs(c2*aomax) > cutoff) THEN
                fa(ind+1:ind+i) = fa(ind+1:ind+i) &
                        + c1*valg(1:i) &
                        + c2*aox(1:i)
            END IF
        END DO

        IF (.NOT.urohf) RETURN

        DO i = 1, l1
            ind = (i-1)*l1
            c1 = aox(i)
            c2 = ftotwt*(0.5_fp*dumb * aox(i) &
                               + dumbx*gaox(i) &
                               + dumby*gaoy(i) &
                               + dumbz*gaoz(i))
            valg(i) = c2
            aomax   = max(aomax,  abs(c1))
            bmaxval = max(bmaxval,abs(c2))
            IF (abs(c1*bmaxval) + abs(c2)*aomax > cutoff) &
                    fb(ind+1:ind+i) = fb(ind+1:ind+i) &
                                    + c1*valg(1:i) &
                                    + c2*aox(1:i)
        END DO
    ELSE ! Meta-GGA belongs here for the `ELSE` case
        !old code
        aomax   = 0.0
        amaxval = 0.0
        bmaxval = 0.0
        DO i = 1, l1
            ind = (i-1)*l1
            c1 = aox(i)
            c2 = ftotwt*(0.5_fp*duma * aox(i) &
                              + dumax*gaox(i) &
                              + dumay*gaoy(i) &
                              + dumaz*gaoz(i))
            valg(i) = c2
            aomax   = max(aomax,   abs(c1))
            amaxval = max(amaxval, abs(c2))
            IF (abs(c1*amaxval) + abs(c2*aomax) > cutoff) THEN
                fa(ind+1:ind+i) = fa(ind+1:ind+i) &
                        + c1*valg(1:i) &
                        + c2*aox(1:i) &
                        + ftotwt*dmgga*(gaox(i)*gaox(1:i) &
                                      + gaoy(i)*gaoy(1:i) &
                                      + gaoz(i)*gaoz(1:i))
            END IF
        END DO

        IF (.NOT.urohf) RETURN
        DO i = 1, l1
            ind = (i-1)*l1
            c1 = aox(i)
            c2 = ftotwt*(0.5_fp*dumb * aox(i) &
                               + dumbx*gaox(i) &
                               + dumby*gaoy(i) &
                               + dumbz*gaoz(i))
            valg(i) = c2
            aomax   = max(aomax,  abs(c1))
            bmaxval = max(bmaxval,abs(c2))
            IF (abs(c1*bmaxval) + abs(c2)*aomax > cutoff) &
                    fb(ind+1:ind+i) = fb(ind+1:ind+i) &
                                    + c1*valg(1:i) &
                                    + c2*aox(1:i) &
                                    + ftotwt*dmggb*(gaox(i)*gaox(1:i) &
                                                  + gaoy(i)*gaoy(1:i) &
                                                  + gaoz(i)*gaoz(1:i))
        END DO
        ! end old code
    END IF

 END SUBROUTINE

!> @brief Compute grid XC contribution to the Kohn-Sham matrix
!> @param[in]    coeffa     MO coefficients, alpha-spin
!> @param[in]    coeffb     MO coefficients, beta-spin
!> @param[inout] fa         KS matrix, alpha-spin
!> @param[inout] fb         KS matrix, beta-spin
!> @param[out]   eexc       XC energy
!> @param[out]   totele     electronic denisty integral
!> @param[out]   totkin     kinetic energy integral
!> @param[in]    atmxvec    array of atomic X coordinates
!> @param[in]    atmyvec    array of atomic Y coordinates
!> @param[in]    atmzvec    array of atomic Z coordinates
!> @param[in]    mxAngMom   max. needed ang. mom. value (incl. derivatives)
!> @param[in]    l1         basis set size
!> @param[in]    isGGA      .TRUE. if GGA/mGGA functional used
!> @param[in]    urohf      .TRUE. if open-shell calculation
!> @author Vladimir Mironov
 SUBROUTINE dmatd_blk(coeffa,coeffb,fa,fb, &
                   eexc,totele,totkin, &
                   atmxvec,atmyvec,atmzvec, &
                   mxAngMom,l1,isGGA,urohf)
    USE metaGGA, ONLY: needtau
    USE omp_lib, ONLY: omp_get_num_threads
    USE mod_dft_molgrid, ONLY: molGrid
    USE mod_grid_storage, ONLY: grid_3d_t
    IMPLICIT NONE

    LOGICAL, INTENT(IN) :: urohf, isGGA
    INTEGER, INTENT(IN) :: mxAngMom, l1
    REAL(KIND=fp), CONTIGUOUS, INTENT(IN) :: &
            atmxvec(:,:), atmyvec(:,:), atmzvec(:,:)
    REAL(KIND=fp), INTENT(INOUT) :: eexc, totele, totkin
    REAL(KIND=fp), INTENT(INOUT) :: coeffa(l1,l1), coeffb(l1,l1)
    REAL(KIND=fp), INTENT(INOUT) :: fa(*),fb(*)

    COMMON /DFGRID/ dftthr,dftgthr,swoff,sw0,bslrd(137),ndftfg, &
                    nrad,nthe,nphi,nrad0,nthe0,nphi0, &
                    nangpt(mxgrid),nangpt0(mxgrid),sg1,jans
        REAL(KIND=fp) :: dftthr,dftgthr,swoff,sw0,bslrd
        INTEGER :: ndftfg,nrad,nthe,nphi,nrad0,nthe0,nphi0, &
                   nangpt,nangpt0,jans
        LOGICAL :: sg1

    COMMON /DFPRUN/ prunerads(mxgrid,mxgridtyp), &
                    pruneatoms(2,mxgridtyp), &
                    iprunecuts(mxatm),ntotgridpoints(mxatm), &
                    ngrids,maxang,ngridtyps
        REAL(KIND=fp) :: prunerads, pruneatoms
        INTEGER :: iprunecuts,ntotgridpoints,ngrids,maxang,ngridtyps

    COMMON /INFOA / nat,ich,mul,num,nqmt,ne,na,nb, &
                    zan(mxatm),c(3,mxatm),ian(mxatm)
        INTEGER :: nat,ich,mul,num,nqmt,ne,na,nb,ian
        REAL(KIND=fp) :: zan,c

    COMMON /IOFILE/ ir,iw,ip,is,ipk,idaf,nav,ioda(950)
        INTEGER :: ir,iw,ip,is,ipk,idaf,nav,ioda

    COMMON /LMOEDA/ gdtola,gjtola,gktola,tktola,vtola, &
                    gdtolb,gjtolb,gktolb,tktolb,vtolb, &
                    ecorl,excor
        REAL(KIND=fp) :: gdtola,gjtola,gktola,tktola,vtola, &
            gdtolb,gjtolb,gktolb,tktolb,vtolb,ecorl,excor

    COMMON /PAR   / me,master,nproc,ibtyp,iptim,goparr,dskwrk,maswrk
        INTEGER :: me,master,nproc,ibtyp,iptim
        LOGICAL :: dskwrk, maswrk, goparr

    COMMON /SCFOPT/ convhf,maxit,mconv,npunch,npreo(4),fshift
        REAL(KIND=fp) :: convhf,fshift
        INTEGER :: maxit,mconv,npunch,npreo

    COMMON /OUTPUT/ nprint,itol,icut,normf,normp,nopk
      INTEGER :: nprint,itol,icut,normf,normp,nopk

    LOGICAL :: dlb

    REAL(KIND=fp) :: exec, ecorl1, totgradxyz(3), &
        dftthrs, wcutoff, rcutoff, ccutoff, tol

    INTEGER :: noa, nob, npt, next, ipt, i0
    INTEGER :: iSlice, iAng, iRad
    INTEGER :: nPtSlc, nAngPts, nRadPts, wtStart, iAtom

    REAL(KIND=fp), ALLOCATABLE :: fa2(:), fb2(:)

    TYPE(xc_engine_t), ALLOCATABLE :: integrator

    INTEGER :: i, j
    INTEGER :: numNzPts
    INTEGER :: numMolPts
    INTEGER :: nNzSh
    LOGICAL :: done

    INTEGER :: needWtSync

    INTEGER :: myjob
    INTEGER :: iChunk, chunkSize
    INTEGER :: numThreads
    INTEGER :: maxPts, limPts

    tol = RLN10*itol/TWO

!   ITOL sets threshold for integrals. Here it is applied to
!   wavefunction, thus it is halved. It is found that this is a
!   rather conservative setting and one gets about 11 sig. figures
!   of accuracy in the total energy with the default ITOL=20.
!   See DFTGAO and DFTG2AO.
    !write(*,*) 'isGGA, needtau=', isGGA, needtau

    noa = na
    nob = nb
    npt = nrad*maxang

!   Set cut-offs for the density RCUTOFF and weight WCUTOFF
!   RCUTOFF is set to depend upon SCF density conv. and the grid size
!   WCUTOFF is a cell volume and we set it to a fixed value.
!   Most cells have large volume (about 97% have volume >1e-14)
    dftthrs = dftthr
    IF (dftthr==ZERO) dftthr = 1.0d-04/(npt*nat)
    wcutoff = 1.0d-08/(npt*nat)
    rcutoff = convhf/(npt*nat)
    ccutoff = 1.0d-03/(npt*nat)
    IF (dftthr<1.1D-15) THEN
       wcutoff = 1.0d-15
       rcutoff = 1.0d-15
       ccutoff = 1.0d-15
    END IF

    ALLOCATE(fa2(l1*l1), source=0.0_fp)
    IF (urohf) THEN
        ALLOCATE(fb2(l1*l1), source=0.0_fp)
    ELSE
        ALLOCATE(fb2(1), source=0.0_fp)
    END IF

    DO j = 1, l1
        coeffa(:,j) = coeffa(:,j) * nosp_basis%bfnrm(j)
    END DO
    IF (urohf) THEN
        DO j = 1, l1
            coeffb(:,j) = coeffb(:,j) * nosp_basis%bfnrm(j)
        END DO
    END IF

    dlb = goparr .AND. ibtyp==1

    exec     = ZERO
    ecorl1   = ZERO
    totele   = ZERO
    totkin   = ZERO
    totgradxyz = ZERO

    next  = -1
    myjob = -1

    numMolPts = 0

    maxPts = molGrid%maxSlicePts
    limPts = molGrid%maxNRadTimesNAng

!$omp parallel &
!$omp   private(iChunk, chunkSize, numThreads, done) &
!$omp   private(iSlice, numNzPts, integrator) &
!$omp   private(iAng, iRad, iPt, nNzSh) &
!$omp   private(nPtSlc, nAngPts, nRadPts, wtStart, iAtom), &
!$omp   reduction(+:numMolPts, exec, ecorl1, totele, totgradxyz, totkin, fa2, fb2) &
!$omp   reduction(+:needWtSync)

    numThreads = 1
!$  numThreads = omp_get_num_threads()
    chunkSize = max(1,molGrid%nSlices/(nproc*4))
    IF (chunkSize/numThreads>40) THEN
        chunkSize = 40*numThreads
    END IF
    ALLOCATE(integrator)
    CALL integrator%init(l1,nat,mxAngMom,maxPts,limPts, &
            0,urohf,isGGA,needtau,tol,.TRUE.,&
            noa,coeffa,nob,coeffb)
!    CALL integrator%echo

    numMolPts = 0

    done = .FALSE.
    needWtSync = 0

    DO iChunk = 1, molGrid%nSlices, chunkSize
        IF (dlb) THEN
            !$omp barrier
            !$omp master
            myjob = myjob + 1
            IF (myjob>next) CALL ddi_dlbnext(next)
            !$omp end master
            !$omp flush(myjob,next)
            !$omp barrier
            IF (myjob/=next) CYCLE
        ELSE
            IF (mod(iChunk/chunkSize,nproc)/=me) CYCLE
        END IF

!$omp do schedule(dynamic)
slc: DO iSlice = iChunk, min(molGrid%nSlices, iChunk-1+chunkSize)

        CALL molGrid%getSliceData(iSlice, integrator%xyzw)

        !write(*,*) iSlice
        !IF (nPtSlc==0) CYCLE

        nPtSlc  = molGrid%nTotPts(iSlice)
        nAngPts = molGrid%nAngPts(iSlice)
        nRadPts = molGrid%nRadPts(iSlice)
        iAtom   = molGrid%idOrigin(iSlice)
        wtStart = molGrid%wtStart(iSlice)
        !write(*,*) 'reset orb ptr'
        CALL integrator%resetOrbPointers(nPtSlc)

        !write(*,*) 'comp AOs'
        numNzPts = 1
        DO iAng = 1, nAngPts
            DO iRad = 1, nRadPts
                iPt = (iAng-1)*nRadPts+iRad
                IF ( integrator%xyzw(iPt,4)==0.0D0 ) EXIT
                IF (abs(integrator%xyzw(iPt,4))<wcutoff) CYCLE

                integrator%wts(numNzPts) = integrator%xyzw(iPt,4)

                CALL integrator%compAOs( &
                  atmxvec, atmyvec, atmzvec, &
                  iPt, numNzPts, iAtom, nNzSh)

                IF (nNzSh==0) THEN
                    ! the ugly way to clear variable in parallel
                    molGrid%totWts(wtStart+iPt-1,iAtom) = &
                        -(nproc-1)*molGrid%totWts(wtStart+iPt-1,iAtom)
                    needWtSync = 1
                    CYCLE
                END IF

                numNzPts = numNzPts + 1

            END DO
        END DO

        IF (numNzPts==1) THEN
            ! the ugly way to clear variable in parallel
            molGrid%nTotPts(iSlice) = -(nproc-1)*nPtSlc
            CYCLE
        END IF

        ! the ugly way to set variable in parallel
        molGrid%nTotPts(iSlice) = -(nproc-1)*nPtSlc + nproc*numNzPts
        numMolPts = numMolPts + numNzPts
        numNzPts = numNzPts - 1

        ! (Re)map pointers
        CALL integrator%resetXCPointers(numNzPts)

        !write(*,*) 'comp xc'
        CALL integrator%compXC(exec, ecorl1, totele, totkin, totgradxyz)

        !write(*,*) 'upd ks'
!       Uptate Kohn-Sham matrix
        CALL integrator%updKS(fa2, fb2)

        !write(*,*) 'end'
    END DO slc
!$omp end do nowait
    END DO

    DEALLOCATE(integrator)
!$omp end parallel

    IF (dlb) CALL ddi_dlbreset

    CALL molGrid%syncNumPts(numMolPts)
    CALL ddi_gsumi(FLAG_SYNC,needWtSync,1)
    IF (needWtSync>0) CALL molGrid%syncWts

    eexc = eexc + exec
    ecorl = ecorl + ecorl1

    dftthr = dftthrs

    !write(*,*) 'new totele=', totele
    !write(*,*) 'new totkin=', totkin

    DO j = 1, l1
        coeffa(:,j) = coeffa(:,j) / nosp_basis%bfnrm(j)
    END DO

    i0 = 0
    DO i = 1, l1
        fa(i0+1:i0+i) = fa(i0+1:i0+i) + &
                  fa2((i-1)*l1+1:(i-1)*l1+i) &
                  * nosp_basis%bfnrm(i) &
                  * nosp_basis%bfnrm(1:i)
        i0 = i0 + i
    END DO
    DEALLOCATE(fa2)

    IF (urohf) THEN
        DO j = 1, l1
            coeffb(:,j) = coeffb(:,j) / nosp_basis%bfnrm(j)
        END DO

        i0 = 0
        DO i = 1, l1
            fb(i0+1:i0+i) = fb(i0+1:i0+i) + &
                    fb2((i-1)*l1+1:(i-1)*l1+i) &
                    * nosp_basis%bfnrm(i) &
                    * nosp_basis%bfnrm(1:i)
            i0 = i0 + i
        END DO
    END IF
    DEALLOCATE(fb2)

 END SUBROUTINE

!> @brief Compute grid XC contribution to the nuclear gradient
!> @note  Weight derivatives are not applied here. The gradient seems
!>  to be good enough even using fairly poor grids. However, I do
!>  not recomment to use it for numerical Hessian calculation until
!>  weight derivatives are implemented
!> @param[in]    da        density matrix, alpha-spin
!> @param[in]    db        density matrix, beta-spin
!> @param[inout] dedft     nuclear gradient
!> @param[out]   totele    electronic denisty integral
!> @param[out]   totkin    kinetic energy integral
!> @param[in]    atmxvec   array of atomic X coordinates
!> @param[in]    atmyvec   array of atomic Y coordinates
!> @param[in]    atmzvec   array of atomic Z coordinates
!> @param[in]    mxAngMom  max. needed ang. mom. value (incl. derivatives)
!> @param[in]    l1        basis set size
!> @param[in]    isGGA     .TRUE. if GGA/mGGA functional used
!> @param[in]    urohf     .TRUE. if open-shell calculation
!> @author Vladimir Mironov
 SUBROUTINE derexc_blk(da,db, dedft, &
                   totele,totkin, &
                   atmxvec,atmyvec,atmzvec, &
                   mxAngMom,l1,isGGA,urohf)
    USE metaGGA, ONLY: needtau
    USE omp_lib, ONLY: omp_get_num_threads
    USE mod_dft_molgrid
    USE mod_grid_storage, ONLY: grid_3d_t
    IMPLICIT NONE

    LOGICAL, INTENT(IN) :: urohf, isGGA
    INTEGER, INTENT(IN) :: mxAngMom, l1
    REAL(KIND=fp), CONTIGUOUS, INTENT(IN) :: &
            atmxvec(:,:), atmyvec(:,:), atmzvec(:,:)
    REAL(KIND=fp), INTENT(INOUT) :: totele, totkin
    REAL(KIND=fp), INTENT(INOUT) :: da(l1,*), db(l1,*), dedft(3,nat)

    COMMON /DFGRID/ dftthr,dftgthr,swoff,sw0,bslrd(137),ndftfg, &
                    nrad,nthe,nphi,nrad0,nthe0,nphi0, &
                    nangpt(mxgrid),nangpt0(mxgrid),sg1,jans
        REAL(KIND=fp) :: dftthr,dftgthr,swoff,sw0,bslrd
        INTEGER :: ndftfg,nrad,nthe,nphi,nrad0,nthe0,nphi0, &
                   nangpt,nangpt0,jans
        LOGICAL :: sg1

    COMMON /DFPRUN/ prunerads(mxgrid,mxgridtyp), &
                    pruneatoms(2,mxgridtyp), &
                    iprunecuts(mxatm),ntotgridpoints(mxatm), &
                    ngrids,maxang,ngridtyps
        REAL(KIND=fp) :: prunerads, pruneatoms
        INTEGER :: iprunecuts,ntotgridpoints,ngrids,maxang,ngridtyps

    COMMON /INFOA / nat,ich,mul,num,nqmt,ne,na,nb, &
                    zan(mxatm),c(3,mxatm),ian(mxatm)
        INTEGER :: nat,ich,mul,num,nqmt,ne,na,nb,ian
        REAL(KIND=fp) :: zan,c

    COMMON /IOFILE/ ir,iw,ip,is,ipk,idaf,nav,ioda(950)
        INTEGER :: ir,iw,ip,is,ipk,idaf,nav,ioda

    COMMON /LMOEDA/ gdtola,gjtola,gktola,tktola,vtola, &
                    gdtolb,gjtolb,gktolb,tktolb,vtolb, &
                    ecorl,excor
        REAL(KIND=fp) :: gdtola,gjtola,gktola,tktola,vtola, &
            gdtolb,gjtolb,gktolb,tktolb,vtolb,ecorl,excor

    COMMON /PAR   / me,master,nproc,ibtyp,iptim,goparr,dskwrk,maswrk
        INTEGER :: me,master,nproc,ibtyp,iptim
        LOGICAL :: dskwrk, maswrk, goparr

    COMMON /SCFOPT/ convhf,maxit,mconv,npunch,npreo(4),fshift
        REAL(KIND=fp) :: convhf,fshift
        INTEGER :: maxit,mconv,npunch,npreo

    COMMON /OUTPUT/ nprint,itol,icut,normf,normp,nopk
      INTEGER :: nprint,itol,icut,normf,normp,nopk

    LOGICAL :: dlb

    REAL(KIND=fp) :: exec, ecorl1, totgradxyz(3), &
        dftthrs, wcutoff, rcutoff, ccutoff, tol, eexc

    INTEGER :: noa, nob, npt, next, ipt
    INTEGER :: iSlice, iAng, iRad
    INTEGER :: nPtSlc, nAngPts, nRadPts, wtStart, iAtom

    TYPE(xc_engine_t), ALLOCATABLE, TARGET :: integrator

    INTEGER :: j
    INTEGER :: numNzPts
    INTEGER :: numMolPts
    INTEGER :: maxPts, limPts
    INTEGER :: nNzSh
    LOGICAL :: done

    INTEGER :: needWtSync

    INTEGER :: myjob
    INTEGER :: iChunk, chunkSize
    INTEGER :: numThreads

    REAL(KIND=fp), ALLOCATABLE :: bfGrad(:,:), bfGrad_priv(:,:)
    REAL(KIND=fp), ALLOCATABLE :: da2(:,:), db2(:,:)

    tol = RLN10*itol/TWO

!   ITOL sets threshold for integrals. Here it is applied to
!   wavefunction, thus it is halved. It is found that this is a
!   rather conservative setting and one gets about 11 sig. figures
!   of accuracy in the total energy with the default ITOL=20.
!   See DFTGAO and DFTG2AO.
    !write(*,*) 'isGGA, needtau=', isGGA, needtau

    noa = na
    nob = nb
    npt = nrad*maxang

!   Set cut-offs for the density RCUTOFF and weight WCUTOFF
!   RCUTOFF is set to depend upon SCF density conv. and the grid size
!   WCUTOFF is a cell volume and we set it to a fixed value.
!   Most cells have large volume (about 97% have volume >1e-14)
    dftthrs = dftthr
    IF (dftthr==ZERO) dftthr = 1.0d-04/(npt*nat)
    wcutoff = 1.0d-08/(npt*nat)
    rcutoff = convhf/(npt*nat)
    ccutoff = 1.0d-03/(npt*nat)
    IF (dftthr<1.1D-15) THEN
       wcutoff = 1.0d-15
       rcutoff = 1.0d-15
       ccutoff = 1.0d-15
    END IF

    ALLOCATE(da2(l1,l1))
    DO j = 1, l1
        da2(:,j) = da(:,j) * nosp_basis%bfnrm(j) &
                           * nosp_basis%bfnrm(1:l1)
    END DO
    IF (urohf) THEN
        ALLOCATE(db2(l1,l1))
        DO j = 1, l1
            db2(:,j) = db(:,j) * nosp_basis%bfnrm(j) &
                               * nosp_basis%bfnrm(1:l1)
        END DO
    END IF

    dlb = goparr .AND. ibtyp==1

    exec     = ZERO
    ecorl1   = ZERO
    totele   = ZERO
    totkin   = ZERO
    totgradxyz = ZERO

    next  = -1
    myjob = -1

    numMolPts = 0

    maxPts = molGrid%maxSlicePts
    limPts = molGrid%maxNRadTimesNAng

    ALLOCATE(bfGrad(l1,3), source=0.0_fp)

!$omp parallel &
!$omp   private(iChunk, chunkSize, numThreads, done) &
!$omp   private(iSlice, numNzPts, integrator) &
!$omp   private(iAng, iRad, iPt, nNzSh, j) &
!$omp   private(nPtSlc, nAngPts, nRadPts, wtStart, iAtom) &
!$omp   private(bfGrad_priv) &
!$omp   reduction(+:bfGrad, dedft) &
!$omp   reduction(+:numMolPts, exec, ecorl1, totele, totgradxyz, totkin) &
!$omp   reduction(+:needWtSync)

    numThreads = 1
!$  numThreads = omp_get_num_threads()
    chunkSize = max(1,molGrid%nSlices/(nproc*4))
    IF (chunkSize/numThreads>40) THEN
        chunkSize = 40*numThreads
    END IF
    ALLOCATE(bfGrad_priv(l1,3), source=0.0_fp)
    ALLOCATE(integrator)
    CALL integrator%init(l1,nat,mxAngMom,maxPts, limPts, &
            1,urohf,isGGA,needtau,tol,.FALSE.,&
            noa,da2,nob,db2)
!    CALL integrator%echo

    numMolPts = 0

    done = .FALSE.
    needWtSync = 0

    DO iChunk = 1, molGrid%nSlices, chunkSize
        IF (dlb) THEN
            !$omp barrier
            !$omp master
            myjob = myjob + 1
            IF (myjob>next) CALL ddi_dlbnext(next)
            !$omp end master
            !$omp flush(myjob,next)
            !$omp barrier
            IF (myjob/=next) CYCLE
        ELSE
            IF (mod(iChunk/chunkSize,nproc)/=me) CYCLE
        END IF

!$omp do schedule(dynamic)
slc: DO iSlice = iChunk, min(molGrid%nSlices, iChunk-1+chunkSize)

        CALL molGrid%getSliceData(iSlice, integrator%xyzw)

        !write(*,*) iSlice
        !IF (nPtSlc==0) CYCLE

        nPtSlc  = molGrid%nTotPts(iSlice)
        nAngPts = molGrid%nAngPts(iSlice)
        nRadPts = molGrid%nRadPts(iSlice)
        iAtom   = molGrid%idOrigin(iSlice)
        wtStart = molGrid%wtStart(iSlice)
        !write(*,*) 'reset orb ptr'
        CALL integrator%resetOrbPointers(nPtSlc)

        !write(*,*) 'comp AOs'
        numNzPts = 1
        DO iAng = 1, nAngPts
            DO iRad = 1, nRadPts
                iPt = (iAng-1)*nRadPts+iRad
                IF ( integrator%xyzw(iPt,4)==0.0 ) EXIT
                IF (abs(integrator%xyzw(iPt,4))<wcutoff) CYCLE

                integrator%wts(numNzPts) = integrator%xyzw(iPt,4)

                CALL integrator%compAOs( &
                  atmxvec, atmyvec, atmzvec, &
                  iPt, numNzPts, iAtom, nNzSh)

                IF (nNzSh==0) THEN
                    ! the ugly way to clear variable in parallel
                    molGrid%totWts(wtStart+iPt-1,iAtom) = &
                        -(nproc-1)*molGrid%totWts(wtStart+iPt-1,iAtom)
                    needWtSync = 1
                    CYCLE
                END IF

                numNzPts = numNzPts + 1

            END DO
        END DO

        IF (numNzPts==1) THEN
            ! the ugly way to clear variable in parallel
            molGrid%nTotPts(iSlice) = -(nproc-1)*nPtSlc
            CYCLE
        END IF

        ! the ugly way to set variable in parallel
        molGrid%nTotPts(iSlice) = -(nproc-1)*nPtSlc + nproc*numNzPts
        numMolPts = numMolPts + numNzPts
        numNzPts = numNzPts - 1

        ! (Re)map pointers
        CALL integrator%resetXCPointers(numNzPts)

        !write(*,*) 'comp xc'
!       Compute XC and XC atomic gradient
        CALL integrator%compXCgrad(bfGrad_priv, exec, ecorl1, totele, totkin)
        !write(*,*) 'end comp xc'

!!       Normalize BF gradient
!        bfGrad_priv(:,1) = bfGrad_priv(:,1) * nosp_basis%bfnrm
!        bfGrad_priv(:,2) = bfGrad_priv(:,2) * nosp_basis%bfnrm
!        bfGrad_priv(:,3) = bfGrad_priv(:,3) * nosp_basis%bfnrm

!       Update BF gradient
        bfGrad = bfGrad + bfGrad_priv
        !DO j = 1, nosp_basis%nshell
        !    ASSOCIATE(kat  => nosp_basis%katom(j), &
        !              kloc => nosp_basis%kloc(j), &
        !              knbf => nosp_basis%kmax(j)-nosp_basis%kmin(j)+1 )
        !    IF (kat/=iAtom) &
        !    bfGrad(kloc:kloc+knbf-1,:) = bfGrad(kloc:kloc+knbf-1,:) &
        !            + bfGrad_priv(kloc:kloc+knbf-1,:)
        !    END ASSOCIATE
        !END DO

!       !Use translation invariance
        !IF (dft_wt_der) THEN
        !dedft(1,iAtom) = dedft(1,iAtom) + sum(bfGrad_priv(:,1))
        !dedft(2,iAtom) = dedft(2,iAtom) + sum(bfGrad_priv(:,2))
        !dedft(3,iAtom) = dedft(3,iAtom) + sum(bfGrad_priv(:,3))
        !END IF

!       Clear temp gradient storage
        bfGrad_priv = 0.0
        !write(*,*) 'upd ks'

        !write(*,*) 'end'
    END DO slc
!$omp end do nowait
    END DO
    DEALLOCATE(integrator)
    DEALLOCATE(bfGrad_priv)
!$omp end parallel

    IF (dlb) CALL ddi_dlbreset

!    write(*,*)
!    DO i = 1, nat
!        write(*,'(I8,*(ES16.5))') i, dedft(1:3,i)
!    END DO
!    write(*,*)
    DO j = 1, nosp_basis%nshell
        ASSOCIATE(kat  => nosp_basis%katom(j), &
                  kloc => nosp_basis%kloc(j), &
                  knbf => nosp_basis%kmax(j)-nosp_basis%kmin(j)+1 )
        !write(*,'(4I8,*(ES16.5))') i, kloc, knbf, kat, -bfGrad(i,1:3)
        !write(*,'(4I8,*(ES16.5))') j, kloc, kloc+knbf-1, kat, -sum(bfGrad(kloc:kloc+knbf-1,1)), -sum(bfGrad(kloc:kloc+knbf-1,2)), -sum(bfGrad(kloc:kloc+knbf-1,3))
        dedft(1,kat) = dedft(1,kat) - sum(bfGrad(kloc:kloc+knbf-1,1))
        dedft(2,kat) = dedft(2,kat) - sum(bfGrad(kloc:kloc+knbf-1,2))
        dedft(3,kat) = dedft(3,kat) - sum(bfGrad(kloc:kloc+knbf-1,3))
        END ASSOCIATE
    END DO
!    write(*,'(/A)') 'xc deriv contribution:'
!    DO i = 1, nat
!        write(*,'(*(F20.10))') dedft(1:3,i)
!    END DO
!    write(*,*)

    CALL molGrid%syncNumPts(numMolPts)
    CALL ddi_gsumi(FLAG_SYNC,needWtSync,1)
    IF (needWtSync>0) CALL molGrid%syncWts

    eexc = eexc + exec
    ecorl = ecorl + ecorl1

    dftthr = dftthrs

!    write(*,*) 'new totele=', totele
!    write(*,*) 'new totexc=', exec
!    write(*,*) 'new ecorl =', ecorl
!    write(*,*) 'new totkin=', totkin

    DEALLOCATE(da2)
    IF (urohf) DEALLOCATE(db2)
    DEALLOCATE(bfGrad)


 END SUBROUTINE

 END MODULE mod_dft_gridint
