program test_dgesv
  implicit none
  integer :: n = 1
  integer :: nrhs = 1
  integer :: lda = 1
  integer :: ldb = 1
  integer :: info
  double precision, dimension(1) :: a = (/ 1.0d0 /)
  double precision, dimension(1) :: b = (/ 3.0d0 /)
  integer, dimension(1) :: ipiv

  call dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)

  if (info.ne.0) then
    write (*,*) "Check failed."
    call exit(1)
  end if

  write (*,*) "Check was successful."
  stop
end program

