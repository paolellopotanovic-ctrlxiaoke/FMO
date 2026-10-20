! MODULE PREC
!>    @author  Vladimir Mironov
!>
!>    @brief   Contains constants for floating
!>             point number precision
!
!     REVISION HISTORY:
!>    @date _May, 2016_ Initial release
!
MODULE prec

    IMPLICIT NONE

    PUBLIC

    INTEGER, PARAMETER :: &

!       Single precision:
        sp = selected_real_kind(6, 37),     &

!       Double precision:
!       dp = 8, &
        dp = selected_real_kind(15, 307),   &

!       Quad precision:
!            If compiler lacks support of Quad precision,
!               Douple precision is used.
!            `selected_real_kind(33, 4931)` returns `-1` if
!               compiler does not support such precision.
!               Probably, it is actual only for PGI.
        qp = max(selected_real_kind(15, 307),   &
                 selected_real_kind(33, 4931)), &

!       Default floating point precision:
        fp = dp

END MODULE prec
