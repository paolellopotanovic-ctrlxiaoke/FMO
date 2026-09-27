module strings
  implicit none
  public
contains
  !*MODULE strings DECK to_upper
  !> @brief  set all symbol in line to upper case
  !> @author Igor S. Gerasimov
  !> @date   Sep, 2019 --Initial release--
  !> @date   May, 2021 Moved to strings
  !> @param  line - (inout)
  subroutine to_upper(line)
     integer, parameter ::       &
       ial   = iachar('a'),      &
       izl   = iachar('z'),      &
       idiff = ial - iachar('A')
    character(len=:), allocatable, intent(inout) :: line
    !internal variables
    integer :: i, c
    do i = 1, len(line)
      c = iachar(line(i:i))
      if(c .ge. ial .and. c .le. izl) then
        line(i:i) = achar(c-idiff)
      end if
    end do
  end subroutine to_upper
  !*MODULE strings DECK count_substring
  !> @brief  This function return count of substring in string
  !> @author Igor S. Gerasimov
  !> @date   Sep, 2019 --Initial release--
  !> @date   May, 2021 Moved to strings
  !> @param  substring - (in)
  !> @param  string    - (in)
  integer function count_substring(substring, string) result(res)
    character(len=*), intent(in) :: substring, string
    ! internal variables
    character(len=:), allocatable :: tmp_string
    res = 0
    tmp_string = string
    do
      if(INDEX(tmp_string, substring) .eq. 0) exit
      res = res + 1
      tmp_string = tmp_string(INDEX(tmp_string, substring) + 1:)
    end do
    return
  end function count_substring
  !*MODULE strings DECK remove_spaces
  !> @brief  This routine remove not needed spaces from section line
  !> @detail This routine used revert reading of lines.
  !>         Example, that this routine do:
  !>       > DFTTYP = PBE0 BASNAM=APC4 , ACC5
  !>       < DFTTYP=PBE0 BASNAM=APC4,ACC5
  !> @author Igor S. Gerasimov
  !> @date   Sep, 2019 --Initial release--
  !> @date   May, 2021 Moved to strings
  !> @param  line - (inout) worked line
  subroutine remove_spaces(line)
    character(len=:), allocatable, intent(inout) :: line
    ! internal variables
    character(len=:), allocatable :: tmp_line, res_line
    integer :: i, ind, ind_end
    logical :: skip, first
    tmp_line = line
    line     = ""
    do
      res_line = ""
      ind_end = index(tmp_line,"=")
      first = .true.
      skip = .false.
      ind = ind_end
      if(ind .eq. 0) ind = len(tmp_line)
      do i = ind, 1, -1
        if(tmp_line(i:i) /= " ") then
          res_line = tmp_line(i:i) // res_line
          if(tmp_line(i:i) .ne. "=" .and. first .and. ind_end .ne. 0) then
            skip = .true.
            first = .false.
          end if
        else if(skip) then
          res_line = " " // res_line
          skip = .false.
        end if
      end do
      line = line // trim(adjustl(res_line))
      if(ind_end .eq. 0) exit
      tmp_line = tmp_line(ind+1:)
    end do
  end subroutine remove_spaces
  !*MODULE strings DECK index_ith
  !> @brief   This function return index of i'th substring in string
  !> @details if ind is negative, backward search will be
  !>          result is equal 0 if search was failed
  !> @author  Igor S. Gerasimov
  !> @date    May, 2021 --Initial release--
  !> @param   substring - (in)
  !> @param   string    - (in)
  !> @param   ith       - (in) index of needed substring
  integer function index_ith(substring, string, ith) result(res)
    character(len=*), intent(in) :: substring, string
    integer,          intent(in) :: ith
    ! internal variables
    character(len=:), allocatable :: tmp_string
    integer :: i, diff
    res = 0
    tmp_string = string
    if(ith .gt. 0) then
      do i = 1, ith
        diff = INDEX(tmp_string, substring)
        res = res + diff
        if(diff .eq. 0) then
          res = 0
          exit
        end if
        tmp_string = tmp_string(INDEX(tmp_string, substring) + 1:)
      end do
    else if(ith .lt. 0) then
      do i = -1, ith, -1
        res = INDEX(tmp_string, substring, back=.true.)
        if(res .eq. 0) exit
        tmp_string = tmp_string(:INDEX(tmp_string, substring, back=.true.) - 1)
      end do
    else
      res = 0
    end if
  end function index_ith
end module strings
