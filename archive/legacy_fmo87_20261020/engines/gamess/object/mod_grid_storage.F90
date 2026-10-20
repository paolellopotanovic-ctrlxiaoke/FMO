!> @brief Module to store data of DFT atomic quadratures
!> @author Vladimir Mironov
 MODULE mod_grid_storage

    USE mx_limits, ONLY: MXGRID
    USE prec, ONLY: fp
    IMPLICIT NONE

!*******************************************************************************

!> @brief Type to store 3D quadrature grid
!> @author Vladimir Mironov
    TYPE :: grid_3d_t
        !< X coordinates
        REAL(KIND=fp), ALLOCATABLE :: x(:)
        !< Y coordinates
        REAL(KIND=fp), ALLOCATABLE :: y(:)
        !< Z coordinates
        REAL(KIND=fp), ALLOCATABLE :: z(:)
        !< quadrature weights
        REAL(KIND=fp), ALLOCATABLE :: w(:)
        !< number of points
        INTEGER :: nPts   = 0
        !< grid array index, used for reverse search
        INTEGER :: idGrid = 0
        !< tesselation data: (number or points per zone, tesselation degree)
        INTEGER(KIND=2), ALLOCATABLE :: izones(:,:)
    CONTAINS
        PROCEDURE, PASS :: set   => setGrid
        PROCEDURE, PASS :: get   => getGrid
        PROCEDURE, PASS :: check => checkGrid
        PROCEDURE, PASS :: clear => clearGrid
    END TYPE

!*******************************************************************************

!> @brief Pointer to 3D grid container (to be used in arrays)
!> @author Vladimir Mironov
    TYPE :: grid_3d_pt
        TYPE(grid_3d_t),  POINTER :: p
    END TYPE

!*******************************************************************************

    ! Dummy type, to be changed in future
    TYPE :: grid_3d_stack_t
        TYPE(grid_3d_pt) :: grd(MXGRID)
        INTEGER :: ngrd
    END TYPE

!*******************************************************************************

!> @brief Basic array list type to store 3d grids
!> @note DO not use any method except `p => list%get`
!>  in a performance-critical loop!
!> @author Vladimir Mironov
    TYPE :: list_grid_3d_t
        PRIVATE
        !< Current number of grids stored
        INTEGER, PUBLIC :: nGrids   = 0
        !< Maximum number of grids
        INTEGER         :: maxGrids = 0
        !< Grid data
        TYPE(grid_3d_pt),  ALLOCATABLE :: elem(:)
    CONTAINS
        PROCEDURE :: get     => getListGrid
        PROCEDURE :: findID  => findIDListGrid
        PROCEDURE :: getByID => getByIDListGrid
        PROCEDURE :: set     => setListGrid
        PROCEDURE :: push    => pushListGrid
        PROCEDURE :: pop     => popListGrid
        PROCEDURE :: init    => initListGrid
        PROCEDURE :: clear   => clearListGrid
        PROCEDURE :: delete  => deleteListGrid
        PROCEDURE, PRIVATE :: extend => extendListGrid
    END TYPE

    INTEGER, PARAMETER :: DEFAULT_GRID_CHUNK = 32

!*******************************************************************************

    TYPE(list_grid_3d_t), TARGET :: spherical_grids
    TYPE(grid_3d_t), POINTER :: sph_last => NULL()
!$omp threadprivate(sph_last)
    TYPE(grid_3d_stack_t), TARGET :: grid_stack

!*******************************************************************************

    REAL(KIND=fp), ALLOCATABLE, TARGET :: rad_grid(:), rad_wts(:)

!*******************************************************************************
    PRIVATE

    PUBLIC rad_grid, rad_wts
    PUBLIC save_rad_grid

    PUBLIC grid_3d_t
    PUBLIC grid_3d_pt
    PUBLIC list_grid_3d_t

    PUBLIC init_gridlist
    PUBLIC delete_gridlist
    PUBLIC get_grid
    PUBLIC get_grid_ptr
    PUBLIC set_grid

    PUBLIC grid_stack
    PUBLIC reset_leb_grid_stack
    PUBLIC push_last_leb_grid
CONTAINS

!*******************************************************************************

! Wrappers for legacy Fortran

!> @brief Store radial quadrature
!> @param[in]  npts number of points
!> @param[in]  r    grid abscissae
!> @param[in]  w    grid weights
!> @author Vladimir Mironov
 SUBROUTINE save_rad_grid(npts, r, w)
    INTEGER, INTENT(IN) :: npts
    REAL(KIND=fp), INTENT(IN) :: r(:), w(:)
    rad_grid = r(1:npts)
    rad_wts  = w(1:npts)
 END SUBROUTINE

!-------------------------------------------------------------------------------

!> @brief Reset stack of angular grids
!> @author Vladimir Mironov
 SUBROUTINE reset_leb_grid_stack
    grid_stack%ngrd = 0
 END SUBROUTINE

!-------------------------------------------------------------------------------

!> @brief Push angular grid on the top of stack
!> @author Vladimir Mironov
 SUBROUTINE push_last_leb_grid
    IF (grid_stack%ngrd==MXGRID) CALL abrt
    grid_stack%ngrd = grid_stack%ngrd + 1
    grid_stack%grd(grid_stack%ngrd)%p => sph_last
 END SUBROUTINE

!-------------------------------------------------------------------------------

!> @brief Initialize angular grid storage
!> @author Vladimir Mironov
 SUBROUTINE init_gridlist
    CALL spherical_grids%init
 END SUBROUTINE

!-------------------------------------------------------------------------------

!> @brief Clear angular grid storage
!> @author Vladimir Mironov
 SUBROUTINE delete_gridlist
    CALL spherical_grids%delete
 END SUBROUTINE

!-------------------------------------------------------------------------------

!> @brief Get a pointer to angular grid by specifying its number of points
!> @param[in]  npts  number of points
!> @author Vladimir Mironov
 FUNCTION get_grid_ptr(npts) RESULT(ptr)
    INTEGER, INTENT(IN) :: npts
    TYPE(grid_3d_t), POINTER :: ptr

    ptr => spherical_grids%getByID(npts)

 END FUNCTION

!-------------------------------------------------------------------------------

!> @brief Get angular grid data
!> @param[in]   npts  number of points
!> @param[out]  x        x coordinates
!> @param[out]  y        y coordinates
!> @param[out]  z        z coordinates
!> @param[out]  w        weights
!> @param[out]  izones   tesselation data
!> @param[out]  found    whether the grid was found in the dataset
!> @author Vladimir Mironov
 SUBROUTINE get_grid(npts,x,y,z,w,izones,found)
    INTEGER, INTENT(OUT) :: izones(:,:)
    REAL(KIND=fp), INTENT(OUT) :: x(:), y(:), z(:), w(:)
    INTEGER, INTENT(IN) :: npts
    LOGICAL, INTENT(OUT) :: found
    TYPE(grid_3d_t), POINTER :: tmp

    tmp => spherical_grids%get(npts)
    found = associated(tmp)

    IF (found) THEN
        x = tmp%x
        y = tmp%y
        z = tmp%z
        w = tmp%w
        izones = tmp%izones
        sph_last => tmp
    END IF

 END SUBROUTINE

!-------------------------------------------------------------------------------

!> @brief Set angular grid data
!> @param[in]  npts  number of points
!> @param[in]  x        x coordinates
!> @param[in]  y        y coordinates
!> @param[in]  z        z coordinates
!> @param[in]  w        weights
!> @param[in]  izones   tesselation data
!> @author Vladimir Mironov
 SUBROUTINE set_grid(npts,x,y,z,w,izones)
    REAL(KIND=fp), INTENT(IN) :: x(:), y(:), z(:), w(:)
    INTEGER, INTENT(IN) :: izones(:,:)
    INTEGER, INTENT(IN) :: npts
    !TYPE(grid_3d_t) :: grid
    !CALL grid%set(npts,x,y,z,w,izones)
    !sph_last => spherical_grids%set(grid)
    sph_last => spherical_grids%set(grid_3d_t(x=x, y=y, z=z, w=w, &
                                              npts=npts, izones=izones))
 END SUBROUTINE

!*******************************************************************************

!> @brief Get a pointer to a grid stored in a list
!> @param[in]  npts  number of points
!> @author Vladimir Mironov
 FUNCTION getListGrid(list,npts) RESULT(res)
    CLASS(list_grid_3d_t), TARGET, INTENT(IN) :: list
    TYPE(grid_3d_t), POINTER :: res

    INTEGER, INTENT(IN) :: npts
    INTEGER :: i
    LOGICAL :: found

#ifdef __ibmxl__
    flush(6)
#endif

    found = .FALSE.

    DO i = 1, list%nGrids
        res => list%elem(i)%p
        found = res%check(npts)
        IF (found) EXIT
    END DO

    IF (.NOT.found) res => NULL()

 END FUNCTION getListGrid

!-------------------------------------------------------------------------------

!> @brief Get index of a grid stored in a list
!> @param[in]  npts  number of points
!> @author Vladimir Mironov
 FUNCTION findIDListGrid(list,npts) RESULT(id)
    CLASS(list_grid_3d_t), TARGET, INTENT(IN) :: list
    INTEGER :: id

    INTEGER, INTENT(IN) :: npts
    LOGICAL :: found

    found = .FALSE.

    DO id = 1, list%nGrids
        found = list%elem(id)%p%check(npts)
        IF (found) EXIT
    END DO

    IF (.NOT.found) id = 0

 END FUNCTION findIDListGrid

!-------------------------------------------------------------------------------

!> @brief Get a pointer to a grid stored in a list by its index
!> @param[in]  id    grid index
!> @author Vladimir Mironov
 FUNCTION getByIDListGrid(list,id) RESULT(res)
    CLASS(list_grid_3d_t), TARGET, INTENT(IN) :: list
    TYPE(grid_3d_t), POINTER :: res
    INTEGER, INTENT(IN) :: id

    res => NULL()
    IF ( id>=1 .AND. id<=list%nGrids ) res => list%elem(id)%p
#ifdef __ibmxl__
    IF (.NOT.associated(res)) THEN
        CALL ABRTX("This is IBM XL workaround")
    END IF
#endif

 END FUNCTION getByIDListGrid

!-------------------------------------------------------------------------------

!> @brief Set a grid in the list. Return a pointer to the list entry.
!> @param[in]  grid  grid data
!> @author Vladimir Mironov
 FUNCTION setListGrid(list,grid) RESULT(ptr)
    CLASS(list_grid_3d_t), INTENT(INOUT) :: list
    TYPE(grid_3d_t), INTENT(IN) :: grid
    TYPE(grid_3d_t), POINTER :: ptr

    ptr => list%get(grid%npts)

    IF (associated(ptr)) THEN
        ptr = grid
    ELSE
        ptr => list%push(grid)
    END IF

END FUNCTION setListGrid

!-------------------------------------------------------------------------------

!> @brief Add a new grid to the list. Return a pointer to the list entry.
!> @param[in]  grid  grid data
!> @author Vladimir Mironov
 FUNCTION pushListGrid(list,grid) RESULT(ptr)
    CLASS(list_grid_3d_t), INTENT(INOUT) :: list
    TYPE(grid_3d_t), INTENT(IN) :: grid
    TYPE(grid_3d_t), POINTER :: ptr

    IF (list%nGrids == list%maxGrids) CALL list%extend

    list%nGrids = list%nGrids + 1
    ALLOCATE(list%elem(list%nGrids)%p)!, source=grid)
    list%elem(list%nGrids)%p = grid

    ptr => list%elem(list%nGrids)%p

    ptr%idGrid = list%nGrids

 END FUNCTION pushListGrid

!-------------------------------------------------------------------------------

!> @brief Pop a grid from the top of the list.
!> @author Vladimir Mironov
 FUNCTION popListGrid(list) RESULT(grid)
    CLASS(list_grid_3d_t), TARGET, INTENT(INOUT) :: list
    TYPE(grid_3d_t) :: grid

    IF (list%nGrids/=0) THEN
        grid = list%elem(list%nGrids)%p
        DEALLOCATE(list%elem(list%nGrids)%p)
        list%nGrids = list%nGrids - 1
    END IF

 END FUNCTION popListGrid


!-------------------------------------------------------------------------------

!> @brief Initialize grid list, possibly specifying its size
!> @param[in]  iSize   desired list size
!> @author Vladimir Mironov
 SUBROUTINE initListGrid(list,iSize)
    CLASS(list_grid_3d_t) :: list
    INTEGER, OPTIONAL, INTENT(IN) :: iSize
    INTEGER :: isz

    isz = DEFAULT_GRID_CHUNK
    IF (present(iSize)) isz = iSize

    IF (allocated(list%elem)) THEN
        CALL list%clear
        IF (list%maxGrids<isz) THEN
            DEALLOCATE(list%elem)
        END IF
    END IF

    IF (.NOT. allocated(list%elem)) THEN
        ALLOCATE(list%elem(isz))
    END IF

    list%nGrids = 0
    list%maxGrids = max(isz,list%maxGrids)

 END SUBROUTINE initListGrid

!-------------------------------------------------------------------------------

!> @brief Clear grid list
!> @author Vladimir Mironov
 SUBROUTINE clearListGrid(list)
    CLASS(list_grid_3d_t) :: list
    INTEGER :: i
    DO i = 1, list%nGrids
        IF (associated(list%elem(i)%p)) DEALLOCATE(list%elem(i)%p)
    END DO
    list%nGrids = 0
 END SUBROUTINE clearListGrid

!-------------------------------------------------------------------------------

!> @brief Destroy grid list
!> @author Vladimir Mironov
 SUBROUTINE deleteListGrid(list)
    CLASS(list_grid_3d_t) :: list
    CALL list%clear
    DEALLOCATE(list%elem)
    list%maxGrids = 0
 END SUBROUTINE deleteListGrid

!-------------------------------------------------------------------------------

!> @brief Extend grid list to a new size
!> @param[in] iSize  size to extend the list
!> @author Vladimir Mironov
 SUBROUTINE extendListGrid(list, iSize)
    CLASS(list_grid_3d_t), INTENT(INOUT) :: list
    TYPE(grid_3d_pt), ALLOCATABLE :: elem_new(:)

    INTEGER, OPTIONAL, INTENT(IN) :: iSize
    INTEGER :: isz

    isz = DEFAULT_GRID_CHUNK
    IF (present(iSize)) isz = iSize

    ALLOCATE(elem_new(list%maxGrids+isz), source=list%elem)
    CALL move_alloc(from=elem_new, to=list%elem)
    list%maxGrids = list%maxGrids + isz

 END SUBROUTINE

!*******************************************************************************

!> @brief Set 3D grid data
!> @param[in]  npts  number of points
!> @param[in]  x        x coordinates
!> @param[in]  y        y coordinates
!> @param[in]  z        z coordinates
!> @param[in]  w        weights
!> @param[in]  izones   tesselation data
!> @author Vladimir Mironov
 SUBROUTINE setGrid(grid, npts, x, y, z, w, izones)
    CLASS(grid_3d_t), INTENT(INOUT) :: grid
    REAL(KIND=fp), INTENT(IN) :: x(:), y(:), z(:), w(:)
    INTEGER, INTENT(IN) :: izones(:,:)
    INTEGER, INTENT(IN) :: npts

    grid%npts   = npts
    grid%x      = x
    grid%y      = y
    grid%z      = z
    grid%w      = w
    grid%izones = int(izones,kind=2)

 END SUBROUTINE setGrid

!-------------------------------------------------------------------------------

!> @brief Get 3D grid data
!> @param[in]   npts  number of points
!> @param[out]  x        x coordinates
!> @param[out]  y        y coordinates
!> @param[out]  z        z coordinates
!> @param[out]  w        weights
!> @param[out]  izones   tesselation data
!> @author Vladimir Mironov
 SUBROUTINE getGrid(grid, npts, x, y, z, w, izones)
    CLASS(grid_3d_t), INTENT(IN) :: grid
    REAL(KIND=fp), INTENT(OUT) :: x(:), y(:), z(:), w(:)
    INTEGER, INTENT(OUT) :: izones(:,:)
    INTEGER, INTENT(OUT) :: npts

    npts   = grid%npts
    x      = grid%x
    y      = grid%y
    z      = grid%z
    w      = grid%w
    izones = grid%izones

 END SUBROUTINE getGrid

!-------------------------------------------------------------------------------

!> @brief Check if it is the grid you are looking for.
!> @param[in]   npts  number of points
!> @author Vladimir Mironov
 FUNCTION checkGrid(grid, npts) RESULT(ok)
    CLASS(grid_3d_t) :: grid
    INTEGER, INTENT(IN) :: npts
    LOGICAL :: ok
    ok = .FALSE.
    IF (grid%npts==npts) ok = .TRUE.
 END FUNCTION checkGrid

!-------------------------------------------------------------------------------

!> @brief Destroy grid
!> @author Vladimir Mironov
 SUBROUTINE clearGrid(grid)
    CLASS(grid_3d_t), INTENT(INOUT) :: grid
    grid%npts = 0
    DEALLOCATE(grid%x, grid%y, grid%z, grid%w)
    DEALLOCATE(grid%izones)
 END SUBROUTINE clearGrid

!-------------------------------------------------------------------------------

 END MODULE mod_grid_storage
