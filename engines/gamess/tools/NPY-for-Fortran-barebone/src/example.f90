program main
    use npy, only: save_npy
    implicit none
    real(8) :: a(10,20), b(10)
    integer(8) :: i, j

    do i = 1,size(a,1)
        do j = 1,size(a,2)
            a(i,j) = i * j
        enddo
    enddo
    call save_npy("a.npy", a)

    do i = 1,size(b,1)
        b(i) = 2 * i
    enddo
    call save_npy("b.npy", b)
end program main
