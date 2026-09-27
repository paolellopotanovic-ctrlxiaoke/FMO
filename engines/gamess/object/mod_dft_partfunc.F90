 MODULE mod_dft_partfunc

    USE prec, ONLY: fp
    IMPLICIT NONE

    !> @brief Interface for func: double -> double
    ABSTRACT INTERFACE
        PURE REAL(KIND=fp) FUNCTION func_d_d(x)
            IMPORT
            REAL(KIND=fp), INTENT(IN) :: x
        END FUNCTION func_d_d
    END INTERFACE

!> @brief Type for partition function calculation
    TYPE partition_function
        !< Values beyond (-limit,+limit) interval considered as 0.0 or 1.0
        REAL(KIND=fp) :: limit = 1.0_fp
        !< Compute partition function value
        PROCEDURE(func_d_d), NOPASS, POINTER :: eval
        !< Compute partition function derivative
        PROCEDURE(func_d_d), NOPASS, POINTER :: deriv
    CONTAINS
        PROCEDURE :: set => set_partition_function
    END TYPE

    INTEGER, PARAMETER :: PTYPE_BECKE4 = 0
    INTEGER, PARAMETER :: PTYPE_SSF    = 1
    INTEGER, PARAMETER :: PTYPE_ERF    = 2
    INTEGER, PARAMETER :: PTYPE_SMSTP2 = -2
    INTEGER, PARAMETER :: PTYPE_SMSTP3 = 3
    INTEGER, PARAMETER :: PTYPE_SMSTP4 = 4
    INTEGER, PARAMETER :: PTYPE_SMSTP5 = 5

    PRIVATE
    PUBLIC partition_function
    PUBLIC PTYPE_BECKE4
    PUBLIC PTYPE_SSF
    PUBLIC PTYPE_ERF
    PUBLIC PTYPE_SMSTP2
    PUBLIC PTYPE_SMSTP3
    PUBLIC PTYPE_SMSTP4
    PUBLIC PTYPE_SMSTP5


!*******************************************************************************
 CONTAINS

!> @brief Becke's partition function (4th order)
 PURE FUNCTION partf_eval_becke4(x) RESULT(f)
    REAL(KIND=fp) :: f
    REAL(KIND=fp), INTENT(IN) :: x
!    REAL(KIND=fp), PARAMETER :: LIMIT = 1.0, SCALEF = 1.0
    INTEGER :: i
     f = x
     DO i = 1, 4
         f = 0.5_fp*f*(3.0_fp-f*f)
     END DO
     f = 0.5_fp-0.5_fp*f
 END FUNCTION

!> @brief Becke's partition function (4th order) derivative
 PURE FUNCTION partf_diff_becke4(x) RESULT(df)
    REAL(KIND=fp) :: df
    REAL(KIND=fp), INTENT(IN) :: x
!    REAL(KIND=fp), PARAMETER :: LIMIT = 1.0, SCALEF = 1.0
    REAL(KIND=fp), PARAMETER :: FACTOR = 81.0_fp/32.0_fp
    REAL(KIND=fp) :: f
    INTEGER :: i
     f = x
     df = 1.0_fp
     DO i = 1, 4
        df = df*(1.0_fp-f*f)
        f = 0.5_fp*f*(3.0_fp-f*f)
     END DO
     df = -FACTOR*df
 END FUNCTION

!-------------------------------------------------------------------------------

!> @brief SSF partition function
 PURE FUNCTION partf_eval_ssf(x) RESULT(f)
    REAL(KIND=fp) :: f
    REAL(KIND=fp), INTENT(IN) :: x
    REAL(KIND=fp), PARAMETER :: LIMIT = 0.64_fp, SCALEF = 1.0_fp/LIMIT
    REAL(KIND=fp) :: f1, f2, f4

     IF (abs(x)>LIMIT) THEN
        f = 0.5_fp-sign(0.5_fp,x)
        RETURN
     END IF
     f1 = x*SCALEF
     f2 = f1*f1
     f4 = f2*f2
     f = 0.0625_fp * f1 * ((35.0_fp - 35.0_fp*f2) + f4*(21.0_fp - 5.0_fp*f2))
     f = 0.5_fp-0.5_fp*f
 END FUNCTION

!> @brief SSF partition function derivative
 PURE FUNCTION partf_diff_ssf(x) RESULT(df)
    REAL(KIND=fp) :: df
    REAL(KIND=fp), INTENT(IN) :: x
    REAL(KIND=fp), PARAMETER :: LIMIT = 0.64_fp, SCALEF = 1.0_fp/LIMIT
    REAL(KIND=fp) :: f1, f2, f4

     IF (abs(x)>LIMIT) THEN
        df = 0.0_fp
        RETURN
     END IF
     f1 = x*SCALEF
     f2 = f1*f1
     f4 = f2*f2
     df = 0.0625_fp * (   ( 35.0_fp - 105.0_fp*f2) + &
                       f4*(105.0_fp -  35.0_fp*f2) )
     df = -0.5_fp*SCALEF*df
 END FUNCTION

!-------------------------------------------------------------------------------

!> @brief Modified SSF with erf(a*x/(1-x^2)) scaling (like in NWChem)
 PURE FUNCTION partf_eval_erf(x) RESULT(f)
    REAL(KIND=fp) :: f
    REAL(KIND=fp), INTENT(IN) :: x
    REAL(KIND=fp), PARAMETER :: LIMIT = 0.725_fp, SCALEF = 1.0_fp/0.3_fp

     IF (abs(x)>LIMIT) THEN
        f = 0.5_fp-sign(0.5_fp,x)
        RETURN
     END IF
     f = erf(x/(1.0_fp-x**2)*SCALEF)
     f = 0.5_fp-0.5_fp*f
 END FUNCTION

!> @brief Modified SSF with erf(a*x/(1-x^2)) scaling (like in NWChem) derivative
 PURE FUNCTION partf_diff_erf(x) RESULT(df)
    REAL(KIND=fp) :: df
    REAL(KIND=fp), INTENT(IN) :: x
    REAL(KIND=fp), PARAMETER :: LIMIT = 0.725_fp, SCALEF = 1.0_fp/0.3_fp
    REAL(KIND=fp), PARAMETER :: PI = 3.141592653589793D0
    REAL(KIND=fp), PARAMETER :: FACTOR = SCALEF/SQRT(PI)
    REAL(KIND=fp) :: ex, frac

     IF (abs(x)>LIMIT) THEN
        df = 0.0_fp
        RETURN
     END IF
     frac = 1.0_fp/(1.0_fp-x*x)
     ex = exp(-(SCALEF*SCALEF*x*x*frac*frac))
     df = -FACTOR*ex*(1.0_fp+x*x)*frac
 END FUNCTION

!-------------------------------------------------------------------------------

!> @brief Smoothstep function (2nd order)
 PURE FUNCTION partf_eval_smoothstep2(x) RESULT(f)
    REAL(KIND=fp) :: f
    REAL(KIND=fp), INTENT(IN) :: x
    REAL(KIND=fp) :: f1, f2
    REAL(KIND=fp), PARAMETER :: LIMIT = 0.55_fp, SCALEF = 0.5_fp/LIMIT
      IF (abs(x)>LIMIT) THEN
          f = 0.5_fp-sign(0.5_fp,x)
          RETURN
      ELSE
          f1 = 0.5_fp - x*SCALEF
          f2 = f1*f1
          f = f2*f1 * ( 10.0_fp - 15.0_fp*f1 + 6.0_fp*f2 )
      END IF
 END FUNCTION

!> @brief Smoothstep function (2nd order) derivative
 PURE FUNCTION partf_diff_smoothstep2(x) RESULT(df)
    REAL(KIND=fp) :: df
    REAL(KIND=fp), INTENT(IN) :: x
    REAL(KIND=fp) :: f1, f2
    REAL(KIND=fp), PARAMETER :: LIMIT = 0.55_fp, SCALEF = 0.5_fp/LIMIT
      IF (abs(x)>LIMIT) THEN
          df = 0.0_fp
          RETURN
      ELSE
          f1 = 0.5_fp - x*SCALEF
          f2 = f1*f1
          df = f2 * ( 30.0_fp - 60.0_fp*f1 + 30.0_fp*f2 )
      END IF
 END FUNCTION

!-------------------------------------------------------------------------------

!> @brief Smoothstep function (3th order)
 PURE FUNCTION partf_eval_smoothstep3(x) RESULT(f)
    REAL(KIND=fp) :: f
    REAL(KIND=fp), INTENT(IN) :: x
    REAL(KIND=fp) :: f1, f2, f4
    REAL(KIND=fp), PARAMETER :: LIMIT = 0.62_fp, SCALEF = 0.5_fp/LIMIT
      IF (abs(x)>LIMIT) THEN
          f = 0.5_fp-sign(0.5_fp,x)
          RETURN
      ELSE
          f1 = 0.5_fp - x*SCALEF
          f2 = f1*f1
          f4 = f2*f2
          f = f4 * (   (35.0_fp - 84.0_fp*f1) + &
                    f2*(70.0_fp - 20.0_fp*f1) )
      END IF
 END FUNCTION

!> @brief Smoothstep function (3th order) derivative
 PURE FUNCTION partf_diff_smoothstep3(x) RESULT(df)
    REAL(KIND=fp) :: df
    REAL(KIND=fp), INTENT(IN) :: x
    REAL(KIND=fp) :: f1, f2, f4
    REAL(KIND=fp), PARAMETER :: LIMIT = 0.62_fp, SCALEF = 0.5_fp/LIMIT
      IF (abs(x)>LIMIT) THEN
          df = 0.0_fp
          RETURN
      ELSE
          f1 = 0.5_fp - x*SCALEF
          f2 = f1*f1
          f4 = f2*f2
          df = f2*f1 * (   (140.0_fp - 420.0_fp*f1) + &
                        f2*(420.0_fp - 140.0_fp*f1) )
      END IF
 END FUNCTION

!-------------------------------------------------------------------------------

!> @brief Smoothstep function (4th order)
 PURE FUNCTION partf_eval_smoothstep4(x) RESULT(f)
    REAL(KIND=fp) :: f
    REAL(KIND=fp), INTENT(IN) :: x
    REAL(KIND=fp) :: f1, f2, f4
    REAL(KIND=fp), PARAMETER :: LIMIT = 0.69_fp, SCALEF = 0.5_fp/LIMIT
      IF (abs(x)>LIMIT) THEN
          f = 0.5_fp-sign(0.5_fp,x)
          RETURN
      ELSE
          f1 = 0.5_fp - x*SCALEF
          f2 = f1*f1
          f4 = f2*f2
          f = f4*f1 * (   (126.0_fp - 420.0_fp*f1) + &
                       f2*(540.0_fp - 315.0_fp*f1) + &
                       f4*70.0_fp)
      END IF
 END FUNCTION

!> @brief Smoothstep function (4th order) derivative
 PURE FUNCTION partf_diff_smoothstep4(x) RESULT(df)
    REAL(KIND=fp) :: df
    REAL(KIND=fp), INTENT(IN) :: x
    REAL(KIND=fp) :: f1, f2, f4
    REAL(KIND=fp), PARAMETER :: LIMIT = 0.69_fp, SCALEF = 0.5_fp/LIMIT
      IF (abs(x)>LIMIT) THEN
          df = 0.0_fp
          RETURN
      ELSE
          f1 = 0.5_fp - x*SCALEF
          f2 = f1*f1
          f4 = f2*f2
          df = f4 * (   ( 630.0_fp - 2520.0_fp*f1) + &
                     f2*(3780.0_fp - 2520.0_fp*f1) + &
                     f4*630.0_fp)
      END IF
 END FUNCTION

!-------------------------------------------------------------------------------

!> @brief Smoothstep function (5th order)
 PURE FUNCTION partf_eval_smoothstep5(x) RESULT(f)
    REAL(KIND=fp) :: f
    REAL(KIND=fp), INTENT(IN) :: x
    REAL(KIND=fp) :: f1, f2, f4
    REAL(KIND=fp), PARAMETER :: LIMIT = 0.73_fp, SCALEF = 0.5_fp/LIMIT
      IF (abs(x)>LIMIT) THEN
          f = 0.5_fp-sign(0.5_fp,x)
          RETURN
      ELSE
          f1 = 0.5_fp - x*SCALEF
          f2 = f1*f1
          f4 = f2*f2
          f =  f4*f2 * (   ( 462.0_fp - 1980.0_fp*f1) + &
                        f2*(3465.0_fp - 3080.0_fp*f1) + &
                        f4*(1386.0_fp -  252.0_fp*f1))
      END IF
 END FUNCTION

!> @brief Smoothstep function (5th order) derivative
 PURE FUNCTION partf_diff_smoothstep5(x) RESULT(df)
    REAL(KIND=fp) :: df
    REAL(KIND=fp), INTENT(IN) :: x
    REAL(KIND=fp) :: f1, f2, f4
    REAL(KIND=fp), PARAMETER :: LIMIT = 0.73_fp, SCALEF = 0.5_fp/LIMIT
      IF (abs(x)>LIMIT) THEN
          df = 0.0_fp
          RETURN
      ELSE
          f1 = 0.5_fp - x*SCALEF
          f2 = f1*f1
          f4 = f2*f2
          df = f4*f1 *(   ( 2772.0_fp - 13860.0_fp*f1) + &
                       f2*(27720.0_fp - 27720.0_fp*f1) + &
                       f4*(13860.0_fp -  2772.0_fp*f1))
      END IF
 END FUNCTION

!*******************************************************************************

!> @brief Set up partition function parameters
 SUBROUTINE set_partition_function(partfunc, ptype)
    CLASS(partition_function), INTENT(INOUT) :: partfunc
    INTEGER, INTENT(IN) :: ptype

    SELECT CASE (ptype)
    CASE (PTYPE_BECKE4)
      partfunc%limit = 1.0_fp
      partfunc%eval  => partf_eval_becke4
      partfunc%deriv => partf_diff_becke4

    CASE (PTYPE_SSF   )
      partfunc%limit = 0.64_fp
      partfunc%eval  => partf_eval_ssf
      partfunc%deriv => partf_diff_ssf

    CASE (PTYPE_ERF)
      partfunc%limit = 0.725_fp
      partfunc%eval  => partf_eval_erf
      partfunc%deriv => partf_diff_erf

    CASE (PTYPE_SMSTP2)
      partfunc%limit = 0.55_fp
      partfunc%eval  => partf_eval_smoothstep2
      partfunc%deriv => partf_diff_smoothstep2

    CASE (PTYPE_SMSTP3)
      partfunc%limit = 0.62_fp
      partfunc%eval  => partf_eval_smoothstep3
      partfunc%deriv => partf_diff_smoothstep3

    CASE (PTYPE_SMSTP4)
      partfunc%limit = 0.69_fp
      partfunc%eval  => partf_eval_smoothstep4
      partfunc%deriv => partf_diff_smoothstep4

    CASE (PTYPE_SMSTP5)
      partfunc%limit = 0.74_fp
      partfunc%eval  => partf_eval_smoothstep5
      partfunc%deriv => partf_diff_smoothstep5

    CASE DEFAULT
      partfunc%limit = 0.62_fp
      partfunc%eval  => partf_eval_smoothstep3
      partfunc%deriv => partf_diff_smoothstep3
    END SELECT

 END SUBROUTINE

!-------------------------------------------------------------------------------

! SUBROUTINE toupper(s, u)
!     CHARACTER(LEN=*), INTENT(IN)  :: s
!     CHARACTER(LEN=*), INTENT(OUT) :: u
!     INTEGER, PARAMETER :: ISHIFT = (iachar("A")-iachar("a"))
!     INTEGER :: i, code
!
!     DO i = 1, min(len(s), len(u))
!         SELECT CASE (s(i:i))
!         CASE ( "a" : "z" )
!             code = iachar(s(i:i))
!             u(i:i) = achar(code+ISHIFT)
!         CASE DEFAULT
!             u(i:i) = s(i:i)
!         END SELECT
!     END DO
!
! END SUBROUTINE

 END MODULE mod_dft_partfunc
