!>  @brief This module contains types and subroutines to manipulate
!>   basis set
!>  @details The main goal of this module is to split SP(L) type shells
!>   onto pair of S and P shells. It significantly simplifies code for
!>   one- and two-electron integrals, whereas does not impact on other
!>   parts of GAMESS. The SP-free basis set will be stored in `NOSP_BASIS`
!>   variable.
!   REVISION HISTORY:
!>  @date -Sep, 2018- Initial release
!>  @author  Vladimir Mironov
MODULE mod_nosp_basis
    USE prec, ONLY: fp
    USE mx_limits, ONLY: mxsh, mxgtot, mxgsh, mxg2

    IMPLICIT NONE

    TYPE basis_set
        REAL(KIND=fp), DIMENSION(:), ALLOCATABLE :: &
            ex, &   !< Array of primitive Gaussian exponents
            cc, &   !< Array of contraction coefficients
            bfnrm   !< Array of normalization constants
        INTEGER, DIMENSION(:), ALLOCATABLE :: &
            kstart, & !< Locations of the first Gaussian in shells
            katom, &  !< Tells which atom the shell is centered on
            ktype, &  !< Array of shell types, is 1,2,3,4,5,6,7 for S,P,D,F,G,H,I
            kng, &    !< Array of contraction degrees
            kloc, &   !< Indices of shells in the total AO basis
            kmin, &   !< Starting indices of shells
            kmax      !< Ending indices of shells
        INTEGER :: &
            nshell, & !< Number of shells in the basis set
            nprim,  & !< Number of primitive Gaussians in the basis set
            nbf, &    !< Number of basis set functions
            mxcontr   !< Max. contraction degree
    END TYPE


    TYPE(basis_set), SAVE :: nosp_basis !< SP-free basis set

    CHARACTER(LEN=*), PARAMETER :: &

        fmt_split_stats = '(/,&
                   &" L-TYPE SHELLS WERE SPLIT INTO S- AND P-TYPE SHELLS",/,&
                   &" OLD NUMBER OF BASIS SET SHELLS = ",I4,/,&
                   &" NEW NUMBER OF BASIS SET SHELLS = ",I4&
        &)', &

        fmt_split_notfound = '(/,&
                   &" NO L-TYPE SHELLS FOUND"&
        &)'

    PRIVATE
    PUBLIC split_sp_basis
!    PUBLIC set_bfnorms
    PUBLIC basis_set
    PUBLIC nosp_basis
    PUBLIC bas_norm_matrix
    PUBLIC bas_denorm_matrix

CONTAINS

!--------------------------------------------------------------------------------

!>  @brief Main subroutine to initialize `NOSP_BASIS` variable
!
!   PARAMETERS:
!   This subroutine has no parameters
!
!>  @author  Vladimir Mironov
!
!   REVISION HISTORY:
!>  @date -Sep, 2018- Initial release
 SUBROUTINE split_sp_basis

!   GAMESS common blocks specification
    COMMON /NSHEL / ex(mxgtot),cs(mxgtot),cp(mxgtot),cd(mxgtot),     &
                    cf(mxgtot),cg(mxgtot),ch(mxgtot),ci(mxgtot),     &
                    kstart(mxsh),katom(mxsh),ktype(mxsh),kng(mxsh),  &
                    kloc(mxsh),kmin(mxsh),kmax(mxsh),nshell
    REAL(KIND=fp) :: ex, cs, cp, cd, cf, cg, ch, ci
    INTEGER :: kstart, katom, ktype, kng, kloc, kmin, kmax, nshell

    COMMON /PAR   / me,master,nproc,ibtyp,iptim,goparr,dskwrk,maswrk
    INTEGER :: me,master,nproc,ibtyp,iptim
    LOGICAL :: goparr, dskwrk, maswrk


!   Internal variables
    INTEGER :: &
        iprcnt, ind1, ind2, &
        iaocnt, ishel, ishel1

    LOGICAL :: &
        found

    iprcnt  = 0
    iaocnt  = 1

    CALL omp_sp_alloc(nosp_basis)

    nosp_basis%kloc(1)   = 1

    ishel1 = 0

    DO ishel = 1, nshell

!       Check for L-shell
!       If not L-shell:
        IF ((kmax(ishel)-kmin(ishel)+1).NE.4) THEN

            ishel1 = ishel1+1
            nosp_basis%katom(ishel1) = katom(ishel)
            nosp_basis%kng(ishel1)   = kng(ishel)
            nosp_basis%ktype(ishel1) = ktype(ishel)
            nosp_basis%kmin(ishel1)  = kmin(ishel)
            nosp_basis%kmax(ishel1)  = kmax(ishel)

            nosp_basis%kstart(ishel1) = iprcnt + 1

            nosp_basis%kloc(ishel1)   = iaocnt

            iaocnt = iaocnt + kmax(ishel) - kmin(ishel) + 1

            ind1 = kstart(ishel)
            ind2 = kng(ishel)

            nosp_basis%ex(iprcnt+1:iprcnt+ind2) = ex(ind1:ind1+ind2-1)
            nosp_basis%cc(iprcnt+1:iprcnt+ind2) = &
                cs(ind1:ind1+ind2-1) + &
                cp(ind1:ind1+ind2-1) + &
                cd(ind1:ind1+ind2-1) + &
                cf(ind1:ind1+ind2-1) + &
                cg(ind1:ind1+ind2-1) + &
                ch(ind1:ind1+ind2-1) + &
                ci(ind1:ind1+ind2-1)

            iprcnt = iprcnt + kng(ishel)

!       If L-shell was found:
        ELSE

            found = .TRUE.

            ind1 = kstart(ishel)
            ind2 = kng(ishel)

!           S part
            ishel1 = ishel1+1

            nosp_basis%katom(ishel1)  = katom(ishel)
            nosp_basis%kng(ishel1)    = kng(ishel)
            nosp_basis%ktype(ishel1)  = 1
            nosp_basis%kmin(ishel1)   = 1
            nosp_basis%kmax(ishel1)   = 1
            nosp_basis%kstart(ishel1) = iprcnt + 1
            nosp_basis%kloc(ishel1)   = iaocnt

            nosp_basis%ex(iprcnt+1:iprcnt+ind2) = ex(ind1:ind1+ind2-1)
            nosp_basis%cc(iprcnt+1:iprcnt+ind2) = cs(ind1:ind1+ind2-1)

            iaocnt = iaocnt + 1
            iprcnt = iprcnt + ind2

!           P part
            ishel1 = ishel1+1

            nosp_basis%katom(ishel1)  = katom(ishel)
            nosp_basis%kng(ishel1)    = kng(ishel)
            nosp_basis%ktype(ishel1)  = 2
            nosp_basis%kmin(ishel1)   = 2
            nosp_basis%kmax(ishel1)   = 4
            nosp_basis%kstart(ishel1) = iprcnt + 1
            nosp_basis%kloc(ishel1)   = iaocnt

            nosp_basis%ex(iprcnt+1:iprcnt+ind2) = ex(ind1:ind1+ind2-1)
            nosp_basis%cc(iprcnt+1:iprcnt+ind2) = cp(ind1:ind1+ind2-1)

            iaocnt = iaocnt + 3
            iprcnt = iprcnt + ind2

        END IF
    END DO

    nosp_basis%mxcontr = maxval(kng(1:nshell))

!    IF (found) THEN
!        WRITE(*,fmt_split_stats) nshell, nosp_basis%nshell
!    ELSE
!        WRITE(*,fmt_split_notfound)
!    END IF

    CALL set_bfnorms(nosp_basis%bfnrm)

 END SUBROUTINE

!--------------------------------------------------------------------------------

!>  @brief Allocate arrays in `BASIS_SET` type variable
!
!   PARAMETERS:
!   @param[inout]  basis    basis_set type variable
!
!>  @author  Vladimir Mironov
!
!   REVISION HISTORY:
!>  @date -Sep, 2018- Initial release
 SUBROUTINE omp_sp_alloc(basis)
    TYPE(basis_set), INTENT(INOUT) :: basis

!   GAMESS common blocks specification
    COMMON /NSHEL / ex(mxgtot),cs(mxgtot),cp(mxgtot),cd(mxgtot),     &
                    cf(mxgtot),cg(mxgtot),ch(mxgtot),ci(mxgtot),     &
                    kstart(mxsh),katom(mxsh),ktype(mxsh),kng(mxsh),  &
                    kloc(mxsh),kmin(mxsh),kmax(mxsh),nshell
    REAL(KIND=fp) :: ex, cs, cp, cd, cf, cg, ch, ci
    INTEGER :: kstart, katom, ktype, kng, kloc, kmin, kmax, nshell

    INTEGER :: i, num_shell, num_gauss, num_bf

    num_shell = 0
    num_gauss = 0
    num_bf    = 0

    IF (allocated(basis%ex)) CALL omp_sp_destroy(basis) ! cleanup

    DO i = 1, nshell
        IF ((kmax(i)-kmin(i)+1).NE.4) THEN
            num_shell = num_shell + 1
            num_gauss = num_gauss + kng(i)
        ELSE
            num_shell = num_shell + 2
            num_gauss = num_gauss + 2*kng(i)
        END IF
    END DO

    num_bf = kloc(nshell) + kmax(nshell) - kmin(nshell)

    basis%nshell = num_shell
    basis%nprim  = num_gauss
    basis%nbf    = num_bf

    ALLOCATE(basis%ex(num_gauss))
    ALLOCATE(basis%cc(num_gauss))
    ALLOCATE(basis%bfnrm(num_bf))

    ALLOCATE(basis%kstart(num_shell))
    ALLOCATE(basis%katom(num_shell))
    ALLOCATE(basis%ktype(num_shell))
    ALLOCATE(basis%kng(num_shell))
    ALLOCATE(basis%kloc(num_shell))
    ALLOCATE(basis%kmin(num_shell))
    ALLOCATE(basis%kmax(num_shell))

 END SUBROUTINE

!--------------------------------------------------------------------------------

!>  @brief Dellocate arrays in `BASIS_SET` type variable
!
!   PARAMETERS:
!   @param[inout]  basis    basis_set type variable
!
!>  @author  Vladimir Mironov
!
!   REVISION HISTORY:
!>  @date -Sep, 2018- Initial release
 SUBROUTINE omp_sp_destroy(basis)
    TYPE(basis_set), INTENT(INOUT) :: basis

    DEALLOCATE(basis%ex)
    DEALLOCATE(basis%cc)
    DEALLOCATE(basis%bfnrm)

    DEALLOCATE(basis%kstart)
    DEALLOCATE(basis%katom)
    DEALLOCATE(basis%ktype)
    DEALLOCATE(basis%kng)
    DEALLOCATE(basis%kloc)
    DEALLOCATE(basis%kmin)
    DEALLOCATE(basis%kmax)

 END SUBROUTINE

!--------------------------------------------------------------------------------

!>  @brief Initialize array of basis function normalization factors
!
!   PARAMETERS:
!   @param[inout]  p(:)    array of normalization factors
!
!>  @author  Vladimir Mironov
!
!   REVISION HISTORY:
!>  @date -Sep, 2018- Initial release
 SUBROUTINE set_bfnorms(p)
    REAL(KIND=fp), ALLOCATABLE, INTENT(INOUT) :: p(:)

    COMMON /NSHEL / ex(MXGTOT),cs(MXGTOT),cp(MXGTOT),cd(MXGTOT), &
                    cf(MXGTOT),cg(MXGTOT),ch(MXGTOT),ci(MXGTOT), &
                    kstart(MXSH),katom(MXSH),ktype(MXSH),kng(MXSH), &
                    kloc(MXSH),kmin(MXSH),kmax(MXSH),nshell
      REAL(KIND=fp) :: ex,cs,cp,cd,cf,cg,ch,ci
      INTEGER :: kstart,katom,ktype,kng,kloc,kmin,kmax,nshell

    COMMON /SHLNRM/ pnrm(84)
      REAL(KIND=fp) :: pnrm

    INTEGER :: mini, maxi, n, i

    DO i = 1, nshell
        mini = kmin(i)
        maxi = kmax(i)
        n = kloc(i)
        p(n:n+maxi-mini) = pnrm(mini:maxi)
    END DO

 END SUBROUTINE

!--------------------------------------------------------------------------------

!>  @brief Scale matrix `A` with matrix \f$ P \cdot P^T \f$
!>  @details `A` is a packed square matrix, `P` is a column vector
!
!   PARAMETERS:
!   @param[inout]  a(:)    triangular matrix (dimension `LD`*(`LD`+1)/2)
!   @param[in]     p(:)    vector (dimension `LD`)
!   @param[in]     ld      leading dimension of P and matrix A in unpacked form
!
!>  @author  Vladimir Mironov
!
!   REVISION HISTORY:
!>  @date -Sep, 2018- Initial release
 SUBROUTINE bas_norm_matrix(a, p, ld)

    REAL(KIND=fp), INTENT(INOUT) :: a(:)
    INTEGER, INTENT(IN) :: ld
    REAL(KIND=fp), ALLOCATABLE, INTENT(IN) :: p(:)

    INTEGER :: i, n

!$omp parallel do private(i,n)
    DO i = 1, ld
        n = i*(i-1)/2
        a(n+1:n+i) = a(n+1:n+i) * p(i) * p(1:i)
    END DO
!$omp end parallel do

 END SUBROUTINE

!--------------------------------------------------------------------------------

!>  @brief Scale matrix `A` with matrix \f$ 1/P \cdot 1/P^T \f$
!>  @details `A` is a packed square matrix, `P` is a column vector
!
!   PARAMETERS:
!   @param[inout]  a(:)    triangular matrix (dimension `LD`*(`LD`+1)/2)
!   @param[in]     p(:)    vector (dimension `LD`)
!   @param[in]     ld      leading dimension of P and matrix A in unpacked form
!
!>  @author  Vladimir Mironov
!
!   REVISION HISTORY:
!>  @date -Sep, 2018- Initial release
 SUBROUTINE bas_denorm_matrix(a, p, ld)

    REAL(KIND=fp), INTENT(INOUT) :: a(:)
    INTEGER, INTENT(IN) :: ld
    REAL(KIND=fp), ALLOCATABLE, INTENT(INOUT) :: p(:)

    p = 1/p
    CALL bas_norm_matrix(a, p, ld)
    p = 1/p

 END SUBROUTINE

END MODULE
