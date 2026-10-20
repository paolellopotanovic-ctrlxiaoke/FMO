! MODULE MODMNFUN
!>    @author  Sarom Sok Leang
!
!>    @brief   Contains parameters for Minnesota Density Functionals
!
!>    @detail  Replaces the following common blocks:
!>    COMMON /CM12  / IM12
!>    COMMON /CM05  / CM05XF, IM05
!>    COMMON /CM06  / IM06
!>    COMMON /CM08  / IM08
!>    COMMON /CSOGGA/ ISOGGA
!
!     REVISION HISTORY:
!>    @date April, 2018 - Initial release
!>    @date Aug,   2019 - Adding CM05,CM06,CM08,CSOGGA common blocks
!
!>    @param IM05   Toggles M05 parameters
!>    @param IM06   Toggles M06 parameters
!>    @param IM08   Toggles M08/M11 parameters
!>    @param IM12   Toggles MN12/MN15 parameters
!>    @param ISOGGA Toggles SOGGA parameters
!>    @param CM05XF Scale factors for local and non-local contributions
MODULE modmnfun

    IMPLICIT NONE

    INTEGER :: IM05
    INTEGER :: IM06
    INTEGER :: IM08
    INTEGER :: IM12
    INTEGER :: ISOGGA
    DOUBLE PRECISION :: CM05XF

END MODULE modmnfun