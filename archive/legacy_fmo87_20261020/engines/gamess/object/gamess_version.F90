!*MODULE gamess_version
!> @brief  Provides information about GAMESS build
!> @author Igor S. Gerasimov
!> @date   April, 2022
!> @note   This file should be configured
module gamess_version
  use comm_PAR, only: maswrk
  use comm_IOFILE, only: iw
  implicit none
contains
  !> @brief   Show git revision if available
  !> @author  Igor S. Gerasimov
  !> @date    April, 2022
  !> @details Show nothing if git is not available
  subroutine git_commit()
#ifdef __GIT_COMMIT_HASH__
    if (maswrk) then
      write (iw, *) ""
      write (iw, *) "          git commit: " // __GIT_COMMIT_HASH__
      write (iw, *) ""
    end if
#endif
  end subroutine git_commit
end module gamess_version
