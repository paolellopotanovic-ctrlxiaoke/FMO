module npy
    implicit none
    interface save_npy
        module procedure write_int64_vec, write_int64_mtx, write_dbl_vec, write_dbl_mtx
    end interface save_npy
contains
#define GENSUB(SNAME,PTYPE,PNAME,PSHAPE,PSTR) \
    subroutine SNAME(filename, PNAME);\
        character(len=*), intent(in) :: filename;\
        PTYPE, intent(in) :: PNAME PSHAPE;\
        character(len=:), allocatable :: header;\
        character(len=13):: num_str;\
        integer :: i;\
        header = "{'descr':'" // PSTR // "','fortran_order':True,'shape':(";\
        do i = 1, size(shape(PNAME));\
            write (num_str, "(I13)") size(PNAME, i);\
            header = header // trim(adjustl(num_str));\
            if (size(shape(PNAME)).eq.1 .or. i.lt.size(shape(PNAME))) header = header // ",";\
        enddo;\
        header = header // ")}"  // repeat(" ", mod(16-mod(6+1+1+4+len(header)+2+1, 16), 16)) // achar(10);\
        open(NPY_UNIT, file=filename, form="unformatted", access="stream");\
        write(NPY_UNIT) achar(int(Z'93')) // "NUMPY", achar(2) // achar(0), len(header, 4), header, PNAME ;\
        close(NPY_UNIT);\
    end subroutine SNAME
    GENSUB(write_int64_mtx,integer(8),int64_mtx,(:,:),"<i8")
    GENSUB(write_int64_vec,integer(8),int64_vec,(:),"<i8")
    GENSUB(write_dbl_mtx,real(8),dbl_mtx,(:,:),"<f8")
    GENSUB(write_dbl_vec,real(8),dbl_vec,(:),"<f8")
end module npy
