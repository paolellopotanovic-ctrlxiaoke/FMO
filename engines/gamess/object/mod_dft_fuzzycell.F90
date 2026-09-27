 MODULE mod_dft_fuzzycell

    USE prec, ONLY: fp

    USE mx_limits, ONLY: &
            MXGRID, MXGRIDTYP, &
            MXATM, MXGTOT, MXGTOT, MXSH, MXGSH

    USE mod_nosp_basis, ONLY: basis_set

    USE params, ONLY: dft_partfun

    USE mod_dft_partfunc, ONLY: partition_function, PTYPE_BECKE4
    USE mod_grid_storage, ONLY: init_gridlist, get_grid, set_grid
    USE mod_dft_molgrid, ONLY: molGrid

    IMPLICIT NONE

    REAL(KIND=fp), PARAMETER :: HUGEFP = huge(1.0_fp)
    REAL(KIND=fp), PARAMETER :: PI = 3.141592653589793238463_fp
    REAL(KIND=fp), PARAMETER :: FOUR_PI = 4.0_fp*PI

    REAL(KIND=fp), ALLOCATABLE :: prim_mx_dist2(:)
    REAL(KIND=fp), ALLOCATABLE :: shell_mx_dist2(:)

    PRIVATE
    PUBLIC prune_basis
    PUBLIC comp_basis_mxdists
    PUBLIC prim_mx_dist2
    PUBLIC shell_mx_dist2

    PUBLIC dft_fc_blk
!    PUBLIC dft_fc_becke
!    PUBLIC dft_fc_ssf
!    PUBLIC dft_fc_ssf2

 CONTAINS

!-------------------------------------------------------------------------------

!> @brief Find shells and primitives which are significant in a given set
!>  of 3D coordinates
!  TODO: move to basis set related source file
!> @author Vladimir Mironov
 SUBROUTINE prune_basis(inBas, xyzv, xyzat, nSh, nPrim, nBf, &
                 outSh, outShNG, outPrim)
    USE mx_limits, ONLY: MXANG
    TYPE(basis_set), INTENT(IN)           :: inBas
    REAL(KIND=fp), INTENT(IN)             :: xyzv(:,:), xyzat(:)
    INTEGER, INTENT(OUT)                  :: nSh, nPrim, nBf
    INTEGER, CONTIGUOUS, INTENT(OUT)      :: outSh(:), outShNG(:), outPrim(:)

    COMMON /INFOA / nat,ich,mul,num,nqmt,ne,na,nb, &
                    zan(MXATM),c(3,MXATM),ian(MXATM)
        INTEGER :: nat,ich,mul,num,nqmt,ne,na,nb,ian
        REAL(KIND=fp) :: zan,c

    INTEGER :: ish, ig, nCur

    INTEGER, PARAMETER :: NANGBF(MXANG) = [((ish+1)*ish/2, ish=1,MXANG)]

    nSh   = 0
    nPrim = 0
    nBf   = 0
    DO ish = 1, inBas%nshell
        ASSOCIATE ( kng    => inBas%kng(ish), &
                    kstart => inBas%kstart(ish), &
                    ktype  => inBas%ktype(ish), &
                    xyz    => c(1:3,inBas%katom(ish)) - xyzat(1:3))

            nCur = 0
            DO ig = kstart, kstart+kng-1
                IF (bfnz(xyz, xyzv, ig)) THEN
                    nPrim = nPrim + 1
                    nCur  = nCur  + 1
                    outPrim(nPrim) = ig
                END IF
            END DO
            IF (nCur>0) THEN
                nSh = nSh + 1
                nBf = nBf + NANGBF(ktype)
                outSh(nSh)   = ish
                outShNG(nSh) = nCur
            END IF
        END ASSOCIATE
    END DO

 CONTAINS

    LOGICAL FUNCTION bfnz(xyz, xyzv, iPrim)
        REAL(KIND=fp), INTENT(IN) :: xyz(:), xyzv(:,:)
        INTEGER, INTENT(IN) :: iPrim
        INTEGER :: i
        REAL :: r2
        bfnz = .FALSE.
        DO i = 1, ubound(xyzv,1)
            bfnz = sum((xyz(1:3)-xyzv(i,1:3))**2) < prim_mx_dist2(iPrim)
            IF ( bfnz ) EXIT
        END DO
    END FUNCTION

 END SUBROUTINE

!-------------------------------------------------------------------------------

!> @brief Compute maximum extent of basis set primitives up to a given
!>  tolerance.
!> @details Find the largest root of the eqn.: r**n * exp(-a*r**2) = tol,
!>  and fill prim_mx_dist2 array with correspoinding r**2 values
!>  For n==0 (S shells) the solution is trivial.
!>  For n>0 (P,D,F... shells) the equivalent equation is used:
!>      ln(q)/2 - a*q/n - ln(tol)/n == 0, where q = r**2
!>  The solution of this equation is:
!>      q = -n/(2*a) * W_{-1} (-(2*a/n)*tol**(2/n))
!>  where W_{k} (x) - k-th branch of Lambert W function
!>  Assuming that a << 0.5*tol**(-2/n), W_{-1} (x) can be approximated:
!>      W_{-1} (-x) =  log(x) - log(-log(x))
!>  The assumption holds for reasonable basis sets.
!>  Next, the approximated result is then refined by making 1-2
!>  Newton-Raphson steps.
!>  The error of this approximation is (much) less than 10**(-4) Bohr**2
!>  for typical cases:
!>  a < 10^6, n = (1 to 3) (P to F shells), tol = 10**(-10)
!>  a < 10^3, n = (4 to 6) (G to I shells), tol = 10**(-10)
!  TODO: move to basis set related source file
!> @param[in] basis     basis set variable
!> @param[in] mlogtol   -ln(tol) value
!> @author Vladimir Mironov
 SUBROUTINE comp_basis_mxdists(basis, mLogTol)
    USE mx_limits, ONLY: MXANG
    TYPE(basis_set), INTENT(IN) :: basis
    REAL(KIND=fp), INTENT(IN) :: mLogTol

    REAL(KIND=fp) :: tmpLogs(1:MXANG-1), logLogTol
    INTEGER :: ish, i

    IF (allocated(prim_mx_dist2)) DEALLOCATE(prim_mx_dist2)
    ALLOCATE(prim_mx_dist2(basis%nPrim))

    IF (allocated(shell_mx_dist2)) DEALLOCATE(shell_mx_dist2)
    ALLOCATE(shell_mx_dist2(basis%nShell))

    logLogTol = log(mLogTol)
    tmpLogs = [(logLogTol+2.0*mLogTol/i, i = 1, MXANG-1)]

    DO ish = 1, basis%nShell
        shell_mx_dist2(ish) = 0.0
        DO i = basis%kstart(ish), basis%kstart(ish)+basis%kng(ish)-1
            ASSOCIATE(  n => basis%ktype(ish)-1, &
                        a => basis%ex(i), &
                        r2 => prim_mx_dist2(i), &
                        r2sh => shell_mx_dist2(ish) )
                IF (n==0) THEN
!                   Explicit solution:
                    r2 = mLogTol/a
                ELSE IF (n<5) THEN
!                   Approximate result:
                    r2 = 0.5*n/a*(tmpLogs(n)-log(a))
!                   One NR step:
                    r2 = r2*(1 - 2*(0.5*n*log(r2)-a*r2+mLogTol)/(n-2*a*r2))
                ELSE
!                   Approximate result:
                    r2 = 0.5*n/a*(tmplogs(n)-log(a))
!                   Two NR steps:
                    r2 = r2*(1 - 2*(0.5*n*log(r2)-a*r2+mLogTol)/(n-2*a*r2))
                    r2 = r2*(1 - 2*(0.5*n*log(r2)-a*r2+mLogTol)/(n-2*a*r2))
                END IF
                r2sh = max(r2,r2sh)
            END ASSOCIATE
        END DO
    END DO
 END SUBROUTINE

!-------------------------------------------------------------------------------

!> @brief Assemble numerical atomic DFT grids to a molecular grid
!> @param[in]    atmxvec  array of atomic X coordinates
!> @param[in]    atmyvec  array of atomic Y coordinates
!> @param[in]    atmzvec  array of atomic Z coordinates
!> @param[in]    rij      interatomic distances
!> @param[in]    nat      number of atoms
!> @param[in]    curAt    index of current atom
!> @param[in]    rad      effective (e.g. Bragg-Slater) radius of current atom
!> @param[inout] wtab     normalized cell function values for LRD
!> @param[in]    aij      surface shifting factors for Becke's method
!> @author Vladimir Mironov
 SUBROUTINE dft_fc_blk(atmxvec,atmyvec,atmzvec, &
                rij,nat,wtab,aij)

    USE mod_grid_storage, ONLY: &
            get_grid_ptr, grid_3d_t
    USE omp_lib, ONLY: omp_get_num_threads, omp_get_thread_num
    REAL(KIND=fp), PARAMETER :: CHECK = transfer('CHECK   ',1.0_fp)

    INTEGER, INTENT(IN) :: nat
    REAL(KIND=fp), INTENT(IN) :: rij(nat,nat), &
            atmxvec(nat,nat), atmyvec(nat,nat), atmzvec(nat,nat)
    REAL(KIND=fp), INTENT(INOUT) :: wtab(nat,nat,*)
    REAL(KIND=fp), INTENT(IN), OPTIONAL :: aij(nat,nat)

    COMMON /PAR   / me,master,nproc,ibtyp,iptim,goparr,dskwrk,maswrk
      INTEGER :: me,master,nproc,ibtyp,iptim
      LOGICAL :: dskwrk, maswrk, goparr

    COMMON /DFPRUN/ prunerads(MXGRID,MXGRIDTYP),pruneatoms(2,MXGRIDTYP), &
                    iprunecuts(MXATM),ntotgridpoints(MXATM), &
                    ngrids,maxang,ngridtyps
      REAL(KIND=fp) :: prunerads, pruneatoms
      INTEGER :: iprunecuts,ntotgridpoints,ngrids,maxang,ngridtyps

    REAL(KIND=fp) :: braggrad
    EXTERNAL braggrad

    REAL(KIND=fp), ALLOCATABLE :: rijInv(:,:), ri(:), wtintr(:)

    INTEGER :: iChunk, iSlice
    INTEGER :: maxSize, maxPts, nMolPts

    TYPE(partition_function) :: partfunc

!    LOGICAL :: dlb
!    dlb = ibtyp==1

!   Set up selected partition function
    CALL partfunc%set(dft_partfun)

!   Initialize temporary space for atomic cell functions
!   and point-to-atom distances
    ALLOCATE(wtintr(nat), ri(nat))

!   Compute inverse interatomic distances
    ALLOCATE(rijInv(nat,nat))
    WHERE (rij/=0.0_fp)
        rijInv = 1.0_fp/rij
    ELSEWHERE
        rijInv = 0.0_fp
    END WHERE

!   Clear point counters
    maxSize = 0
    maxPts  = 0
    nMolPts = 0
    ntotgridpoints = 0

!   Compute total weights for grid points
!   MPI - static LB, OpenMP - dynamic LB
!$omp parallel do &
!$omp   private(iSlice, iChunk, ri, wtintr) &
!$omp   reduction(max:maxPts, maxSize) &
!$omp   reduction(+:ntotgridpoints, nMolPts) &
!$omp   schedule(dynamic)
    DO iChunk = 1, molGrid%nSlices, nProc
        iSlice = iChunk + me

        IF (iSlice<=molGrid%nSlices) THEN
!         Apply selected algorithm on the current slice
          IF (present(aij)) THEN
            CALL do_bfcSlice(iSlice, partfunc, nat, &
                     atmxvec, atmyvec, atmzvec, &
                     ri, rij, rijInv, wtintr, wtab, aij)
          ELSE
            CALL do_bfcSlice(iSlice, partfunc, nat, &
                     atmxvec, atmyvec, atmzvec, &
                     ri, rij, rijInv, wtintr, wtab)
          END IF

!         Count points:
          nMolPts = nMolPts + molGrid%nTotPts(iSlice)
          ntotgridpoints(molGrid%idOrigin(iSlice)) = &
                  ntotgridpoints(molGrid%idOrigin(iSlice)) &
                + molGrid%nTotPts(iSlice)
          maxPts  = max(maxPts, molGrid%nTotPts(iSlice))
          maxSize = max(maxSize,  molGrid%nAngPts(iSlice) &
                                * molGrid%nRadPts(iSlice))
        END IF

    END DO
!$omp end parallel do

!   Finalize grid info
    molGrid%maxNRadTimesNAng = maxSize
    molGrid%maxSlicePts      = maxPts
    molGrid%nMolPts          = nMolPts

 END SUBROUTINE

!> @brief Compute total weights for points in a slice
!> @param[in]    iSlice    index of current slice
!> @param[in]    partfunc  partition function
!> @param[in]    nat       number of atoms
!> @param[in]    atmxvec   array of atomic X coordinates
!> @param[in]    atmyvec   array of atomic Y coordinates
!> @param[in]    atmzvec   array of atomic Z coordinates
!> @param[inout] ri        tmp array to store point to atoms distances
!> @param[in]    rij       interatomic distances
!> @param[in]    rijInv    inverse interatomic distances
!> @param[inout] wtintr    tmp array to store cell function values
!> @param[inout] wtab      normalized cell function values for LRD
!> @param[in]    aij       surface shifting factors for Becke's method
!> @author Vladimir Mironov
 SUBROUTINE do_bfcSlice(iSlice, partfunc, nAt, &
                 atmxvec, atmyvec, atmzvec, &
                 ri, rij, rijInv, wtintr, wtab, aij)
    USE mod_grid_storage, ONLY: &
            ptrad => rad_grid, wtrad => rad_wts, &
            get_grid_ptr, grid_3d_t

    COMMON /LRDISP/ elrd6,elrd8,elrd10,emult,lrdflg,mltint,dolrd
      REAL(KIND=fp) :: elrd6,elrd8,elrd10,emult
      LOGICAL :: lrdflg,mltint,dolrd

    INTEGER, INTENT(IN) :: iSlice, nAt
    TYPE(partition_function), INTENT(IN) :: partfunc
    REAL(KIND=fp), CONTIGUOUS :: atmxvec(:,:), atmyvec(:,:), atmzvec(:,:), &
             ri(:), rij(:,:), rijInv(:,:), wtintr(:)
    REAL(KIND=fp) :: wtab(nAt,nAt,*)
    REAL(KIND=fp), INTENT(IN), CONTIGUOUS, OPTIONAL :: aij(:,:)

    TYPE(grid_3d_t), POINTER :: curGrid
    INTEGER :: iAng, iRad, iPt, iAtm
    REAL(KIND=fp) :: radwt, r1, xd,yd,zd, xcdnt,ycdnt,zcdnt, wtAngRad
    REAL(KIND=fp) :: wtnrm

    ASSOCIATE ( &
        dummyAtom => molGrid%dummyAtom, &
        rInner    => molGrid%rInner, &
        iAngStart => molGrid%iAngStart(iSlice), &
        iRadStart => molGrid%iRadStart(iSlice), &
        nAngPts   => molGrid%nAngPts(iSlice), &
        nRadPts   => molGrid%nRadPts(iSlice), &
        wtStart   => molGrid%wtStart(iSlice)-1, &
        totWts    => molGrid%totWts, &
        ntp       => molGrid%nTotPts(iSlice), &
        isInner   => molGrid%isInner(iSlice), &
        curAt     => molGrid%idOrigin(iSlice), &
        rad       => molGrid%rAtm(iSlice) )

        ntp = 0

        isInner = 0
        IF (ptRad(iRadStart+nRadPts-1)*rad < rInner(curAt)) THEN
!           Weights of the whole slice are unchanged
            isInner = 1
            ntp = nAngPts * nRadPts
!           Nothing left to do for inner slice
            RETURN
        END IF

        curGrid => get_grid_ptr(molGrid%idAng(iSlice))
        ASSOCIATE ( &
            xAng => curGrid%x(iAngStart:iAngStart+nAngPts-1), &
            yAng => curGrid%y(iAngStart:iAngStart+nAngPts-1), &
            zAng => curGrid%z(iAngStart:iAngStart+nAngPts-1), &
            wAng => curGrid%w(iAngStart:iAngStart+nAngPts-1) )

        DO iAng = 1, nAngPts
rloop:      DO iRad = 1, nRadPts

                r1       =         rad*ptRad(iRadStart+iRad-1)
                radWt    = rad*rad*rad*wtRad(iRadStart+iRad-1)
                wtAngRad = FOUR_PI*radWt*wAng(iAng)

                iPt = (iAng-1)*nRadPts + iRad

!               Quick check for inner quadrature points
                IF (r1 < rInner(curAt)) THEN
!                   Point weight is unchanged
                    totWts(wtStart+iPt,curAt) = wtAngRad
                    ntp = ntp + 1
!                   Next point
                    CYCLE rloop
                END IF

                xd = r1*xAng(iAng)
                yd = r1*yAng(iAng)
                zd = r1*zAng(iAng)

                DO iAtm = 1, nAt
                    IF (dummyAtom(iAtm)) CYCLE
                    xcdnt = xd-atmxvec(iAtm,curAt)
                    ycdnt = yd-atmyvec(iAtm,curAt)
                    zcdnt = zd-atmzvec(iAtm,curAt)
                    ri(iAtm) = sqrt(xcdnt**2+ycdnt**2+zcdnt**2)

!                   Check if the point belongs to any other atom
                    IF (iAtm/=curAt.AND.&
                        (r1-ri(iAtm) > rij(iAtm,curAt) * &
                                       0.5D0*(1.0D0+partfunc%limit))) THEN
!                       This and all next points along the same direction
!                       belong to another atom. Switch to the next angular
!                       point. It is never happen when Becke's function is
!                       selected.
                        EXIT rloop
                    END IF
                END DO

!               If screening fails, follow regular BFC procedure and check
!               all atom pairs
                IF (present(aij)) THEN
                    CALL do_bfc(totWts(wtStart+iPt,curAt), wtnrm, wtintr, wtAngRad, &
                           ri, rijInv, nAt, curAt, dummyAtom, partfunc, aij)
                ELSE
                    CALL do_bfc(totWts(wtStart+iPt,curAt), wtnrm, wtintr, wtAngRad, &
                            ri, rijInv, nAt, curAt, dummyAtom, partfunc)
                END IF

                IF (lrdflg) wtab(1:nat,curAt,iPt) = wtintr(1:nat)*wtnrm

!               Same as before, if weight is zero, the current and
!               all next points along the same direction belong
!               to another atom. Switch to the next angular point.
!               It is never happen when Becke's function is selected.
                IF (totWts(wtStart+iPt,curAt) == 0.0D0) EXIT rloop

                ntp = ntp + 1

            END DO rloop
        END DO
        END ASSOCIATE

    END ASSOCIATE

 END SUBROUTINE

!-------------------------------------------------------------------------------

!> @brief Compute total weight of a grid point
!> @param[out]   wt        total weight
!> @param[out]   wtNorm    sum of cell function values
!> @param[out]   cells     array of cell function values
!> @param[in]    wtAngRad  unmodified weight
!> @param[inout] ri        tmp array to store point to atoms distances
!> @param[in]    dummyAtom array to indicate which atoms are "dummy"
!> @param[in]    rijInv    inverse interatomic distances
!> @param[in]    numAt     number of atoms
!> @param[in]    curAtom   current atom
!> @param[in]    partfunc  partition function
!> @param[in]    aij       surface shifting factors for Becke's method
!> @author Vladimir Mironov
 SUBROUTINE do_bfc(wt, wtNorm, cells, wtAngRad, &
                 ri, rijInv, numAt, curAtom, dummyAtom, partfunc, aij)
    LOGICAL,       CONTIGUOUS, INTENT(IN)  :: dummyAtom(:)
    REAL(KIND=fp), CONTIGUOUS, INTENT(IN)  :: ri(:), rijInv(:,:)
    INTEGER,                   INTENT(IN)  :: numAt, curAtom
    REAL(KIND=fp),             INTENT(IN)  :: wtAngRad
    TYPE(partition_function),  INTENT(IN)  :: partfunc
    REAL(KIND=fp),             INTENT(OUT) :: wt, wtNorm
    REAL(KIND=fp), CONTIGUOUS, INTENT(OUT) :: cells(:)
    REAL(KIND=fp), CONTIGUOUS, INTENT(IN), OPTIONAL :: aij(:,:)

    INTEGER       :: i, j
    REAL(KIND=fp) :: f, mu

    WHERE (dummyAtom)
        cells = 0.0D0
    ELSEWHERE
        cells = 1.0D0
    END WHERE

    DO i = 2, numAt
        IF (dummyAtom(i)) CYCLE
        DO j = 1, i-1
            IF (dummyAtom(j)) CYCLE
            mu = (ri(i)-ri(j))*rijInv(j,i)
            IF (present(aij)) mu = mu + aij(j,i)*(1.0D0-mu*mu)
            f = partfunc%eval(mu)
            cells(i) = cells(i)*abs(f)
            cells(j) = cells(j)*abs(1.0D0-f)
        END DO
    END DO

    wtNorm = 1.0D0/sum(cells(1:numAt))

    wt = cells(curAtom) * wtNorm * wtAngRad

 END SUBROUTINE

!-------------------------------------------------------------------------------

!!> @brief Compute total weight and its derivatives of a grid point
!!> Currently not functioning
!!> @author Vladimir Mironov
! SUBROUTINE do_bfc_deriv(dummyAtom, ri, rij, riInv, rijInv, riv, &
!                 atx, aty, atz, &
!                 numAt, curAtom, wtAngRad, partfunc, &
!                 wt,  cells, dp)
!    LOGICAL,       CONTIGUOUS, INTENT(IN)  :: dummyAtom(:)
!    REAL(KIND=fp), CONTIGUOUS, INTENT(IN)  :: ri(:), riv(:,:), rij(:,:)
!    REAL(KIND=fp), CONTIGUOUS, INTENT(IN)  :: riInv(:), rijInv(:,:)
!    REAL(KIND=fp), CONTIGUOUS, INTENT(IN)  :: atx(:), aty(:), atz(:)
!    INTEGER,                   INTENT(IN)  :: numAt, curAtom
!    REAL(KIND=fp),             INTENT(IN)  :: wtAngRad
!    TYPE(partition_function),  INTENT(IN)  :: partfunc
!    REAL(KIND=fp),             INTENT(OUT) :: wt
!    REAL(KIND=fp), CONTIGUOUS, INTENT(OUT) :: cells(:), dp(:,:)
!
!    INTEGER       :: i, j
!    REAL(KIND=fp) :: f, df, t, mu, dmudi(3), dmudj(3), vi(3), vj(3), wtNorm
!
!    WHERE (dummyAtom)
!        cells = 0.0
!    ELSEWHERE
!        cells = 1.0
!    END WHERE
!    dp = 0.0
!
!    DO i = 2, numAt
!        IF (dummyAtom(i)) CYCLE
!        vi(1) = atx(i)
!        vi(2) = aty(i)
!        vi(3) = atz(i)
!        DO j = 1, i-1
!            IF (dummyAtom(j)) CYCLE
!            vj(1) = atx(j)
!            vj(2) = aty(j)
!            vj(3) = atz(j)
!
!            mu = (ri(i)-ri(j)) * rijInv(j,i)
!
!            dmudi = rij(j,i)*riv(i,:)*riInv(i) - mu*(vi-vj)
!            dmudi = dmudi*rijInv(j,i)*rijInv(j,i)
!
!            dmudj = rij(j,i)*riv(j,:)*riInv(j) + mu*(vi-vj)
!            dmudj = dmudj*rijInv(j,i)*rijInv(j,i)
!
!            f  = partfunc%eval(mu)
!            df = partfunc%deriv(mu)
!
!            t = df/f
!
!            cells(i) = cells(i)*abs(f)
!            cells(j) = cells(j)*abs(1.0-f)
!
!            dp(:,i) = dp(:,i) + t*dmudi
!
!        END DO
!    END DO
!
!    wtNorm = 1.0/sum(cells(1:numAt))
!
!    wt = cells(curAtom) * wtNorm * wtAngRad
!
! END SUBROUTINE

!!-------------------------------------------------------------------------------
!
! SUBROUTINE dft_fc_ssf_radang(nrad,xdat,ydat,zdat,atmxvec,atmyvec,atmzvec,ri, &
!                rij,nat,wtintr,totwt,wght,iangn,ncntr,rad, &
!                maxpts,wtab)
!
!    USE mod_grid_storage, ONLY: ptrad => rad_grid, wtrad => rad_wts
!    REAL(KIND=fp), PARAMETER :: CHECK = transfer('CHECK   ',1.0_fp)
!
!    COMMON /PAR   / me,master,nproc,ibtyp,iptim,goparr,dskwrk,maswrk
!      INTEGER :: me,master,nproc,ibtyp,iptim
!      LOGICAL :: dskwrk, maswrk, goparr
!    COMMON /RUNOPT/ runtyp,exetyp,nevals,nglevl,nhlevl
!      REAL(KIND=fp) :: runtyp,exetyp
!      INTEGER :: nevals,nglevl,nhlevl
!
!    COMMON /DFPRUN/ prunerads(MXGRID,MXGRIDTYP),pruneatoms(2,MXGRIDTYP), &
!                    iprunecuts(MXATM),ntotgridpoints(MXATM), &
!                    ngrids,maxang,ngridtyps
!      REAL(KIND=fp) :: prunerads, pruneatoms
!      INTEGER :: iprunecuts,ntotgridpoints,ngrids,maxang,ngridtyps
!
!    COMMON /LRDISP/ elrd6,elrd8,elrd10,emult,lrdflg,mltint,dolrd
!      REAL(KIND=fp) :: elrd6,elrd8,elrd10,emult
!      LOGICAL :: lrdflg,mltint,dolrd
!
!    INTEGER :: nrad, ncntr, nat
!    REAL(KIND=fp) :: rad
!    REAL(KIND=fp) :: wght(maxpts,nat,*), &
!             xdat(maxpts,nat,*), ydat(maxpts,nat,*), zdat(maxpts,nat,*), &
!             atmxvec(nat,nat),atmyvec(nat,nat),atmzvec(nat,nat), &
!             ri(nat),rij(nat,nat),wtintr(nat),totwt(nat,*)
!    INTEGER :: iangn(nat,2,*),maxpts
!
!    REAL(KIND=fp) :: wtab(nat,nat,*)
!    REAL(KIND=fp) :: braggrad
!    EXTERNAL braggrad
!
!    LOGICAL :: dlb
!    INTEGER :: iang, irad, i, ipt, iptme, iatm, nn
!    INTEGER :: ntotpts, igrid
!    REAL(KIND=fp) :: radwt, r1, xd,yd,zd, xcdnt,ycdnt,zcdnt
!    REAL(KIND=fp) :: wt0
!    REAL(KIND=fp) :: wtnrm
!    REAL(KIND=fp) :: distnn
!    LOGICAL, ALLOCATABLE :: dummyatom(:)
!
!    TYPE(partition_function) :: partfunc
!
!    REAL(KIND=fp) :: rinner
!
!    IF (exetyp==check) THEN
!        totwt(ncntr,1:nrad*maxang) = 0.1
!        RETURN
!    END IF
!
!    CALL partfunc%set(dft_partfun)
!
!    dlb = ibtyp==1
!
!    distnn = HUGEFP
!    nn = nat+1
!    ALLOCATE(dummyatom(nat))
!!   Find nearest neghbour of current atom and the corresponding distance
!!   Also, tag non-real atoms present in the system
!    DO i = 1, nat
!        IF (braggrad(i)==0.0_fp) THEN
!            dummyatom(i) = .TRUE.
!        ELSE
!            dummyatom(i) = .FALSE.
!            IF (i==ncntr) CYCLE
!            IF (rij(i,ncntr)<distnn) THEN
!                nn = i
!                distnn = rij(i,ncntr)
!            END IF
!        END IF
!    END DO
!
!    rinner = 0.5*distnn*(1.0-partfunc%limit)
!    totwt(ncntr,1:nrad:maxang) = 0.0
!
!gl: DO igrid = 1, ngrids
!
!angl:   DO iang = iangn(ncntr,1,igrid), iangn(ncntr,2,igrid)
!
!radl:   DO irad = 1, nrad
!
!            IF (igrid>1.AND.&
!                ptrad(irad)<prunerads(igrid-1,iprunecuts(ncntr))) THEN
!                CYCLE radl
!            END IF
!            IF (ptrad(irad)>=prunerads(igrid,iprunecuts(ncntr))) EXIT radl
!
!            r1 = rad*ptrad(irad)
!
!            ipt = (irad-1)*maxpts+iang
!            !ipt = (iang-1)*nrad+irad
!
!            IF (mod(ipt,nproc)/=me) CYCLE
!
!            iptme = (ipt-1)/nproc+1
!            IF (dlb) iptme = ipt
!
!            radwt = rad*rad*rad*wtrad(irad)
!            wt0 = radwt*wght(iang,ncntr,igrid)
!
!!           Quick check for inner quadrature points
!            IF (r1 < rinner) THEN
!                totwt(ncntr,iptme) = wt0
!                CYCLE
!            END IF
!
!            xd = r1*xdat(iang,ncntr,igrid)
!            yd = r1*ydat(iang,ncntr,igrid)
!            zd = r1*zdat(iang,ncntr,igrid)
!
!            DO iatm = 1, nat
!                IF (dummyAtom(iAtm)) CYCLE
!                xcdnt = xd-atmxvec(iatm,ncntr)
!                ycdnt = yd-atmyvec(iatm,ncntr)
!                zcdnt = zd-atmzvec(iatm,ncntr)
!                ri(iatm) = sqrt(xcdnt**2+ycdnt**2+zcdnt**2)
!!               Quick check if the point belongs to any other atom
!                IF (iatm/=ncntr.AND.&
!                    (r1-ri(iatm) > 0.5*rij(iatm,ncntr) * &
!                                  (1.0+partfunc%limit))) THEN
!                    !totwt(ncntr,iptme) = 0.0
!                    !CYCLE radl
!                    EXIT radl
!                END IF
!            END DO
!
!!           If screening fails, follow regular BFC procedure and check
!!           all atom pairs
!            CALL do_bfc_old(dummyAtom,ri,rij,nAt,ncntr,wt0,partfunc, &
!               totwt(ncntr,iptme), wtnrm, wtintr)
!            IF (lrdflg) wtab(1:nat,ncntr,iPt) = wtintr(1:nat)*wtnrm
!            IF (totwt(ncntr,iptme) == 0.0) EXIT radl
!
!        END DO radl
!        END DO angl
!
!        ntotpts = 0
!        DO irad = 1, nrad
!            ntotpts = ntotpts + iangn(ncntr,2,igrid)-iangn(ncntr,1,igrid) + 1
!        END DO
!        ntotgridpoints(ncntr) = ntotpts
!    END DO gl
!    DEALLOCATE(dummyatom)
!
! END SUBROUTINE
!
!!-------------------------------------------------------------------------------
!
! SUBROUTINE dft_fc_ssf2(nrad,xdat,ydat,zdat,atmxvec,atmyvec,atmzvec,ri,&
!                rij,nat,wtintr,totwt,wght,iangn,ncntr,rad, &
!                maxpts,ptrad,wtrad,wtab)
!
!      REAL(KIND=fp), PARAMETER :: check = transfer('check   ',1.0_fp)
!
!      COMMON /PAR   / me,master,nproc,ibtyp,iptim,goparr,dskwrk,maswrk
!        INTEGER :: me,master,nproc,ibtyp,iptim
!        LOGICAL :: dskwrk, maswrk, goparr
!      COMMON /RUNOPT/ runtyp,exetyp,nevals,nglevl,nhlevl
!        REAL(KIND=fp) :: runtyp,exetyp
!        INTEGER :: nevals,nglevl,nhlevl
!      COMMON /DFPRUN/ prunerads(MXGRID,MXGRIDTYP),pruneatoms(2,MXGRIDTYP), &
!                      iprunecuts(MXATM),ntotgridpoints(MXATM), &
!                      ngrids,maxang,ngridtyps
!        REAL(KIND=fp) :: prunerads, pruneatoms
!        INTEGER :: iprunecuts,ntotgridpoints,ngrids,maxang,ngridtyps
!
!      INTEGER :: nrad, ncntr, nat
!      REAL(KIND=fp) :: rad
!      REAL(KIND=fp) :: wght(maxpts,nat,*), &
!               xdat(maxpts,nat,*), ydat(maxpts,nat,*), zdat(maxpts,nat,*), &
!               atmxvec(nat,nat),atmyvec(nat,nat),atmzvec(nat,nat), &
!               ri(nat),rij(nat,nat),wtintr(nat),totwt(nat,*), &
!               ptrad(nrad),wtrad(nrad)
!      INTEGER :: iangn(nat,2,*),maxpts
!
!      COMMON /LRDISP/ elrd6,elrd8,elrd10,emult,lrdflg,mltint,dolrd
!        REAL(KIND=fp) :: elrd6,elrd8,elrd10,emult
!        LOGICAL :: lrdflg,mltint,dolrd
!
!      REAL(KIND=fp) :: wtab(nat,nat,*)
!
!      REAL(KIND=fp) :: braggrad
!      EXTERNAL braggrad
!
!      LOGICAL :: dlb
!      INTEGER :: irad, i, ipt, iptme, iatm
!      REAL(KIND=fp) :: radwt, r1, xd,yd,zd, xcdnt,ycdnt,zcdnt
!      REAL(KIND=fp) :: wt0
!    REAL(KIND=fp) :: wtnrm
!      REAL(KIND=fp) :: distnn, rptmin
!      LOGICAL, ALLOCATABLE :: dummyatom(:)
!      REAL(KIND=fp) :: nns(0:8)
!      INTEGER :: id_nns(0:8)
!      INTEGER :: m1, m2, m3, moct, nnpt
!      INTEGER :: ntotpts, igrid
!      REAL(KIND=fp) :: dist
!
!      TYPE(partition_function) :: partfunc
!
!      IF (exetyp==check) THEN
!          totwt(ncntr,1:nrad*maxang) = 0.1
!          RETURN
!      END IF
!
!      CALL partfunc%set(dft_partfun)
!
!      dlb = ibtyp==1
!
!      ALLOCATE(dummyatom(nat))
!      DO i = 1, nat
!          dummyatom(i) = braggrad(i)==0.0_fp
!      END DO
!
!!     Find nearest neghbours in each octant of current atom coordinate system
!      nns = HUGEFP
!      nns(0) = 0.0_fp
!      id_nns = nat+1
!      id_nns(0) = ncntr
!      DO i = 1, nat
!          IF (braggrad(i)==0.0_fp) CYCLE
!          IF (i==ncntr) CYCLE
!!         Find which octant w.r.t. NCNTR the atom belongs to
!          m1 = int(sign(0.5_fp,atmxvec(i,ncntr))+0.5)
!          m2 = int(sign(0.5_fp,atmyvec(i,ncntr))+0.5)
!          m3 = int(sign(0.5_fp,atmzvec(i,ncntr))+0.5)
!          moct = 8 - (m1 + 2*m2 + 4*m3)
!          dist = rij(i,ncntr)
!          IF (dist<nns(moct)) THEN
!            nns(moct) = dist
!            id_nns(moct) = i
!          END IF
!      END DO
!      distnn = minval(nns(1:))
!
!      ntotpts = 0
!      igrid = 1
!
!      DO irad = 1, nrad
!
!          radwt = rad*rad*rad*wtrad(irad)
!          r1 = rad*ptrad(irad)
!
!          IF (r1>=prunerads(igrid,iprunecuts(ncntr))*rad) THEN
!             igrid = igrid + 1
!          END IF
!          ntotpts = ntotpts + iangn(ncntr,2,igrid)-iangn(ncntr,1,igrid) + 1
!
!ptloop:   DO i = iangn(ncntr,1,igrid), iangn(ncntr,2,igrid)
!              ipt = (irad-1)*maxpts+i
!              IF (mod(ipt,nproc)/=me) CYCLE
!              iptme = (ipt-1)/nproc+1
!              IF (dlb) iptme = ipt
!
!              wt0 = radwt*wght(i,ncntr,igrid)
!
!!             Quick check for inner quadrature points
!              IF (r1<0.5*distnn*(1.0-partfunc%limit)) THEN
!                  totwt(ncntr,iptme) = wt0
!                  CYCLE
!              END IF
!
!              xd = r1*xdat(i,ncntr,igrid)
!              yd = r1*ydat(i,ncntr,igrid)
!              zd = r1*zdat(i,ncntr,igrid)
!
!!             Neighbouring screening
!!             Find distances from point to all atom neighbours
!              DO iatm = 1, 8
!                  IF (id_nns(iatm)>nat) CYCLE
!                  xcdnt = xd-atmxvec(id_nns(iatm),ncntr)
!                  ycdnt = yd-atmyvec(id_nns(iatm),ncntr)
!                  zcdnt = zd-atmzvec(id_nns(iatm),ncntr)
!                  dist = sqrt(xcdnt**2+ycdnt**2+zcdnt**2)
!
!                  IF (r1-dist>0.5*nns(iatm)*(1.0+partfunc%limit)) THEN
!                      totwt(ncntr,iptme) = 0.0
!                      CYCLE ptloop
!                  END IF
!              END DO
!
!              nnpt = ncntr
!              rptmin = r1
!              DO iatm = 1,nat
!                  IF (dummyAtom(iAtm)) CYCLE
!                  xcdnt = xd-atmxvec(iatm,ncntr)
!                  ycdnt = yd-atmyvec(iatm,ncntr)
!                  zcdnt = zd-atmzvec(iatm,ncntr)
!                  ri(iatm) = sqrt(xcdnt**2+ycdnt**2+zcdnt**2)
!                  IF (ri(iatm)<rptmin) THEN
!                      rptmin = ri(iatm)
!                      nnpt = iatm
!                  END IF
!              END DO
!
!!             Quick check if the point belongs to any neighbour
!              IF (nnpt/=ncntr.AND.&
!                  (r1-rptmin>0.5*rij(ncntr,nnpt)*(1.0+partfunc%limit))) THEN
!                  totwt(ncntr,iptme) = 0.0
!                  CYCLE
!              END IF
!
!!             If screening fails, follow regular BFC procedure and check
!!             all atom pairs
!              CALL do_bfc_old(dummyAtom,ri,rij,nAt,ncntr,wt0,partfunc, &
!                 totwt(ncntr,iptme), wtnrm, wtintr)
!              IF (lrdflg) wtab(1:nat,ncntr,iPt) = wtintr(1:nat)*wtnrm
!
!          END DO ptloop
!      END DO
!      DEALLOCATE(dummyatom)
!      ntotgridpoints(ncntr) = ntotpts
! END SUBROUTINE
!
!!-------------------------------------------------------------------------------
!
! SUBROUTINE dft_fc_ssf(nrad,xdat,ydat,zdat,atmxvec,atmyvec,atmzvec,ri, &
!                rij,nat,wtintr,totwt,wght,iangn,ncntr,rad, &
!                maxpts,ptrad,wtrad,wtab)
!
!      REAL(KIND=fp), PARAMETER :: CHECK = transfer('CHECK   ',1.0_fp)
!
!      COMMON /PAR   / me,master,nproc,ibtyp,iptim,goparr,dskwrk,maswrk
!        INTEGER :: me,master,nproc,ibtyp,iptim
!        LOGICAL :: dskwrk, maswrk, goparr
!      COMMON /RUNOPT/ runtyp,exetyp,nevals,nglevl,nhlevl
!        REAL(KIND=fp) :: runtyp,exetyp
!        INTEGER :: nevals,nglevl,nhlevl
!
!      COMMON /DFPRUN/ prunerads(MXGRID,MXGRIDTYP),pruneatoms(2,MXGRIDTYP), &
!                      iprunecuts(MXATM),ntotgridpoints(MXATM), &
!                      ngrids,maxang,ngridtyps
!        REAL(KIND=fp) :: prunerads, pruneatoms
!        INTEGER :: iprunecuts,ntotgridpoints,ngrids,maxang,ngridtyps
!
!      INTEGER :: nrad, ncntr, nat
!      REAL(KIND=fp) :: rad
!      REAL(KIND=fp) :: wght(maxpts,nat,*), &
!               xdat(maxpts,nat,*), ydat(maxpts,nat,*), zdat(maxpts,nat,*), &
!               atmxvec(nat,nat),atmyvec(nat,nat),atmzvec(nat,nat), &
!               ri(nat),rij(nat,nat),wtintr(nat),totwt(nat,*), &
!               ptrad(nrad),wtrad(nrad)
!      INTEGER :: iangn(nat,2,*),maxpts
!
!      COMMON /LRDISP/ elrd6,elrd8,elrd10,emult,lrdflg,mltint,dolrd
!        REAL(KIND=fp) :: elrd6,elrd8,elrd10,emult
!        LOGICAL :: lrdflg,mltint,dolrd
!
!      REAL(KIND=fp) :: wtab(nat,nat,*)
!      REAL(KIND=fp) :: braggrad
!      EXTERNAL braggrad
!
!      LOGICAL :: dlb
!      INTEGER :: irad, i, ipt, iptme, iatm, nn
!      INTEGER :: ntotpts, igrid
!      REAL(KIND=fp) :: radwt, r1, xd,yd,zd, xcdnt,ycdnt,zcdnt
!      REAL(KIND=fp) :: wt0
!      REAL(KIND=fp) :: wtnrm
!      REAL(KIND=fp) :: distnn
!      LOGICAL, ALLOCATABLE :: dummyatom(:)
!
!      TYPE(partition_function) :: partfunc
!
!      REAL(KIND=fp) :: rinner, rlimit, rmax
!
!      IF (exetyp==check) THEN
!          totwt(ncntr,1:nrad*maxang) = 0.1
!          RETURN
!      END IF
!
!      CALL partfunc%set(dft_partfun)
!
!      dlb = ibtyp==1
!
!      distnn = HUGEFP
!      nn = nat+1
!      ALLOCATE(dummyatom(nat))
!!     Find nearest neghbour of current atom and the corresponding distance
!!     Also, tag non-real atoms present in the system
!      DO i = 1, nat
!          IF (braggrad(i)==0.0_fp) THEN
!              dummyatom(i) = .TRUE.
!          ELSE
!              dummyatom(i) = .FALSE.
!              IF (i==ncntr) CYCLE
!              IF (rij(i,ncntr)<distnn) THEN
!                  nn = i
!                  distnn = rij(i,ncntr)
!              END IF
!          END IF
!      END DO
!
!      rinner = 0.5*distnn*(1.0-partfunc%limit)
!      rlimit = 0.5*distnn*(1.0+partfunc%limit)
!
!      ntotpts = 0
!      igrid = 1
!
!      DO irad = 1, nrad
!
!          radwt = rad*rad*rad*wtrad(irad)
!          r1 = rad*ptrad(irad)
!
!          IF (r1>=prunerads(igrid,iprunecuts(ncntr))*rad) THEN
!             igrid = igrid + 1
!          END IF
!          ntotpts = ntotpts + iangn(ncntr,2,igrid)-iangn(ncntr,1,igrid) + 1
!
!angloop:  DO i = iangn(ncntr,1,igrid), iangn(ncntr,2,igrid)
!              ipt = (irad-1)*maxpts+i
!              IF (mod(ipt,nproc)/=me) CYCLE
!              iptme = (ipt-1)/nproc+1
!              IF (dlb) iptme = ipt
!
!              wt0 = radwt*wght(i,ncntr,igrid)
!
!!             Quick check for inner quadrature points
!              IF (r1 < rinner) THEN
!                  totwt(ncntr,iptme) = wt0
!                  rmax = max(rmax,r1)
!                  CYCLE
!              END IF
!
!              xd = r1*xdat(i,ncntr,igrid)
!              yd = r1*ydat(i,ncntr,igrid)
!              zd = r1*zdat(i,ncntr,igrid)
!
!!             Quick check if the point belongs to the nearest neighbour
!              xcdnt = xd-atmxvec(nn,ncntr)
!              ycdnt = yd-atmyvec(nn,ncntr)
!              zcdnt = zd-atmzvec(nn,ncntr)
!
!              !IF (r1-sqrt(xcdnt**2+ycdnt**2+zcdnt**2) > rlimit) THEN
!              IF (xcdnt**2+ycdnt**2+zcdnt**2 < rinner*rinner) THEN
!                  totwt(ncntr,iptme) = 0.0
!                  CYCLE
!              END IF
!
!              DO iatm = 1, nat
!                  IF (dummyatom(iatm)) CYCLE
!                  xcdnt = xd-atmxvec(iatm,ncntr)
!                  ycdnt = yd-atmyvec(iatm,ncntr)
!                  zcdnt = zd-atmzvec(iatm,ncntr)
!                  ri(iatm) = sqrt(xcdnt**2+ycdnt**2+zcdnt**2)
!!                 Check if the point belongs to any other atom
!                  IF (iatm/=ncntr.AND.&
!                      (r1-ri(iatm)>0.5*rij(iatm,ncntr)*(1.0+partfunc%limit))) THEN
!                      totwt(ncntr,iptme) = 0.0
!                      CYCLE angloop
!                  END IF
!              END DO
!
!!             If screening fails, follow regular BFC procedure and check
!!             all atom pairs
!              CALL do_bfc_old(dummyAtom,ri,rij,nAt,ncntr,wt0,partfunc, &
!                 totwt(ncntr,iptme), wtnrm, wtintr)
!              IF (lrdflg) wtab(1:nat,ncntr,iPt) = wtintr(1:nat)*wtnrm
!
!          END DO angloop
!      END DO
!      DEALLOCATE(dummyatom)
!      ntotgridpoints(ncntr) = ntotpts
! END SUBROUTINE
!
!!-------------------------------------------------------------------------------
!
!!> @brief Compute total weights for points of current atom using original
!!>  Becke's fuzzy cell method
!!> @author Vladimir Mironov
! SUBROUTINE dft_fc_becke(nrad,xdat,ydat,zdat,atmxvec,atmyvec,atmzvec,ri, &
!                rij,nat,aij,wtintr,totwt,wght,iangn,ncntr,rad, &
!                maxpts,ptrad,wtrad,wtab)
!
!      IMPLICIT NONE
!
!      REAL(KIND=fp), PARAMETER :: CHECK = transfer('CHECK   ',1.0_fp)
!
!      COMMON /PAR   / me,master,nproc,ibtyp,iptim,goparr,dskwrk,maswrk
!        INTEGER :: me,master,nproc,ibtyp,iptim
!        LOGICAL :: dskwrk, maswrk, goparr
!      COMMON /RUNOPT/ runtyp,exetyp,nevals,nglevl,nhlevl
!        REAL(KIND=fp) :: runtyp,exetyp
!        INTEGER :: nevals,nglevl,nhlevl
!
!      COMMON /DFPRUN/ prunerads(MXGRID,MXGRIDTYP),pruneatoms(2,MXGRIDTYP), &
!                      iprunecuts(MXATM),ntotgridpoints(MXATM), &
!                      ngrids,maxang,ngridtyps
!        REAL(KIND=fp) :: prunerads, pruneatoms
!        INTEGER :: iprunecuts,ntotgridpoints,ngrids,maxang,ngridtyps
!
!      INTEGER :: nrad, ncntr, nat
!      REAL(KIND=fp) :: rad
!      REAL(KIND=fp) :: wght(maxang,nat,*), &
!               xdat(maxang,nat,*), ydat(maxang,nat,*), zdat(maxang,nat,*), &
!               atmxvec(nat,nat),atmyvec(nat,nat),atmzvec(nat,nat), &
!               ri(nat),rij(nat,nat),wtintr(nat),totwt(nat,*), &
!               aij(nat,nat),ptrad(nrad),wtrad(nrad)
!      INTEGER :: iangn(nat,2,*),maxpts
!
!      COMMON /LRDISP/ elrd6,elrd8,elrd10,emult,lrdflg,mltint,dolrd
!        REAL(KIND=fp) :: elrd6,elrd8,elrd10,emult
!        LOGICAL :: lrdflg,mltint,dolrd
!
!      REAL(KIND=fp) :: wtab(nat,nat,*)
!
!      LOGICAL :: dlb
!      INTEGER :: irad, i, igrid, ipt, iptme, iatm, jatm, ntotpts
!      REAL(KIND=fp) :: radwt, r1, xd,yd,zd, xcdnt,ycdnt,zcdnt
!      REAL(KIND=fp) :: wttot, cutij,cutji, zmuij, xmuij
!
!      TYPE(partition_function) :: partfunc
!
!      IF (exetyp==CHECK) THEN
!         totwt(ncntr,1:nrad*maxang) = 0.1
!         RETURN
!      END IF
!
!      CALL partfunc%set(PTYPE_BECKE4)
!
!      dlb = ibtyp==1
!
!      ntotpts = 0
!      igrid = 1
!
!      DO irad = 1, nrad
!
!         radwt = rad*rad*rad*wtrad(irad)
!         r1 = rad*ptrad(irad)
!
!         IF (r1>=prunerads(igrid,iprunecuts(ncntr))*rad) THEN
!            igrid = igrid + 1
!         END IF
!         ntotpts = ntotpts + iangn(ncntr,2,igrid)-iangn(ncntr,1,igrid) + 1
!
!         DO i = iangn(ncntr,1,igrid), iangn(ncntr,2,igrid)
!            ipt = (irad-1)*maxpts+i
!            IF (mod(ipt,nproc)/=me) CYCLE
!            iptme = (ipt-1)/nproc+1
!            IF (dlb) iptme = ipt
!!
!            xd = r1*xdat(i,ncntr,igrid)
!            yd = r1*ydat(i,ncntr,igrid)
!            zd = r1*zdat(i,ncntr,igrid)
!            DO iatm = 1,nat
!               xcdnt = xd-atmxvec(iatm,ncntr)
!               ycdnt = yd-atmyvec(iatm,ncntr)
!               zcdnt = zd-atmzvec(iatm,ncntr)
!               ri(iatm) = sqrt(xcdnt**2+ycdnt**2+zcdnt**2)
!            END DO
!            wttot = 0.0
!            wtintr(1:nat) = 1.0
!
!            DO iatm = 1,nat
!               IF (abs(aij(1,iatm)+1.0d+00)<1.0d-05) THEN
!                 wtintr(iatm) = 0.0
!                 CYCLE
!               END IF
!
!               DO jatm = 1, iatm-1
!                  IF (abs(aij(jatm,iatm)-1.0d+00)>=1.0d-05) THEN
!                      zmuij = (ri(iatm)-ri(jatm))/rij(iatm,jatm)
!                      xmuij = zmuij + aij(jatm,iatm)*(1.0-zmuij*zmuij)
!                      cutij = partfunc%eval(xmuij)
!                      cutji = 1.0 - cutij
!                      wtintr(iatm) = wtintr(iatm)*cutij
!                      wtintr(jatm) = wtintr(jatm)*cutji
!                  END IF
!               END DO
!            END DO
!
!            wttot = sum(wtintr(1:nat))
!
!            IF (wttot>1.0d-10) THEN
!                totwt(ncntr,iptme) = wtintr(ncntr)/wttot * &
!                                     radwt*wght(i,ncntr,igrid)
!            ELSE
!                totwt(ncntr,iptme) = 0.0
!            END IF
!
!            IF (lrdflg) wtab(1:nat,ncntr,iptme) = wtintr(1:nat)/wttot
!
!         END DO
!      END DO
!      ntotgridpoints(ncntr) = ntotpts
! END SUBROUTINE
!
!!-------------------------------------------------------------------------------
!
! SUBROUTINE do_bfc_old(dummyAtom, ri, rij, numAt, curAtom, wtAngRad, partfunc, &
!                 wt, wtNorm, cells)
!    LOGICAL,       CONTIGUOUS, INTENT(IN)  :: dummyAtom(:)
!    REAL(KIND=fp), CONTIGUOUS, INTENT(IN)  :: ri(:), rij(:,:)
!    INTEGER,                   INTENT(IN)  :: numAt, curAtom
!    REAL(KIND=fp),             INTENT(IN)  :: wtAngRad
!    TYPE(partition_function),  INTENT(IN)  :: partfunc
!    REAL(KIND=fp),             INTENT(OUT) :: wt, wtNorm
!    REAL(KIND=fp), CONTIGUOUS, INTENT(OUT) :: cells(:)
!
!    INTEGER       :: i, j
!    REAL(KIND=fp) :: f, mu
!
!    WHERE (dummyAtom)
!        cells = 0.0
!    ELSEWHERE
!        cells = 1.0
!    END WHERE
!
!    DO i = 2, numAt
!        IF (dummyAtom(i)) CYCLE
!        DO j = 1, i-1
!            IF (dummyAtom(j)) CYCLE
!            mu = (ri(i)-ri(j))/rij(j,i)
!            f = partfunc%eval(mu)
!            cells(i) = cells(i)*abs(f)
!            cells(j) = cells(j)*abs(1.0-f)
!        END DO
!    END DO
!
!    wtNorm = 1.0/sum(cells(1:numAt))
!
!    wt = cells(curAtom) * wtNorm * wtAngRad
!
! END SUBROUTINE

 END MODULE mod_dft_fuzzycell
