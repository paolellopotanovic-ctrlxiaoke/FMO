MODULE mod_dft_molgrid

    USE prec, ONLY: fp
    USE mx_limits, ONLY: MXGRID, MXGRIDTYP, MXATM
    USE mod_nosp_basis, ONLY: nosp_basis, basis_set
    USE mod_grid_storage

    IMPLICIT NONE

    INTEGER, PARAMETER :: MAXDEPTH =  2
    REAL(KIND=fp), PARAMETER :: HUGEFP = huge(1.0_fp)

!*******************************************************************************

!    COMMON /DFPRUN/ pruneRads(MXGRID,MXGRIDTYP),pruneAtoms(2,MXGRIDTYP), &
!                   iPruneCuts(MXATM),nTotGridPoints(MXATM),nGrids,maxAng,nGridTyps
!        REAL(KIND=fp) :: pruneRads, pruneAtoms
!        INTEGER :: iPruneCuts,nTotGridPoints,nGrids,maxAng,nGridTyps
!
    COMMON /INFOA / nat,ich,mul,num,nqmt,ne,na,nb,zan(MXATM),c(3,MXATM),ian(MXATM)
        INTEGER :: nat,ich,mul,num,nqmt,ne,na,nb,ian
        REAL(KIND=fp) :: zan,c


!*******************************************************************************

!> @brief Type to store molecular grid information
    TYPE :: dft_grid_t
        !< current number of slices
        INTEGER :: nSlices = 0
        !< maximum number of slices
        INTEGER :: maxSlices = 0
        !< maximum number of grid points per atom
        INTEGER :: maxAtomPts = 0
        !< maximum number of nonzero grid points per slice
        INTEGER :: maxSlicePts = 0
        !< maximum possible number of grid points per slice
        INTEGER :: maxNRadTimesNAng = 0
        !< total number of nonzero grid points
        INTEGER :: nMolPts = 0

        ! Every slice is a spacially localized subgrid of molecular grid:
        ! slice = (few radial points)x(few angular points)
        ! The following arrays contains data for each slice

        !< index of angular grid
        INTEGER, ALLOCATABLE       :: idAng(:)
        !< index of the first point in angular grid array
        INTEGER, ALLOCATABLE       :: iAngStart(:)
        !< number of angular grid points
        INTEGER, ALLOCATABLE       :: nAngPts(:)
        !< index of the first point in radial grid array
        INTEGER, ALLOCATABLE       :: iRadStart(:)
        !< number of radial grid points
        INTEGER, ALLOCATABLE       :: nRadPts(:)
        !< number of grid points with nonzero weight
        INTEGER, ALLOCATABLE       :: nTotPts(:)
        !< index of atom which owns the slice
        INTEGER, ALLOCATABLE       :: idOrigin(:)
        !< type of chunk -- used for different processing
        INTEGER, ALLOCATABLE       :: chunkType(:)
        !< index of the first point in weights array
        INTEGER, ALLOCATABLE       :: wtStart(:)
        !< isInner == 1 means that the slice is 'inner' -- i.e. its weight is not modified and is equal to wtRad*wtAng
        INTEGER, ALLOCATABLE       :: isInner(:)

        ! The following arrays contains information about atoms
        !< effective radii of atoms
        REAL(KIND=fp), ALLOCATABLE :: rAtm(:)
        !< max radii for 'inner' points
        REAL(KIND=fp), ALLOCATABLE :: rInner(:)
        !< .TRUE. for dummy atoms
        LOGICAL, ALLOCATABLE       :: dummyAtom(:)

        !< array to store grid point weights
        REAL(KIND=fp), ALLOCATABLE :: totWts(:,:)
        !< number of nonzero grid points per atom
        INTEGER, ALLOCATABLE, PRIVATE :: wt_top(:)
    CONTAINS
        PROCEDURE, PASS :: getSliceData    => getSliceData
        PROCEDURE, PASS :: getSliceNonZero => getSliceNonZero
        PROCEDURE, PASS :: exportGrid      => exportGrid
        PROCEDURE, PASS :: setSlice        => setSlice
        PROCEDURE, PASS :: reset           => reset_dft_grid_t
        PROCEDURE, PASS :: sync            => sync_dft_grid_t
        PROCEDURE, PASS :: syncNumPts      => sync_npts_dft_grid_t
        PROCEDURE, PASS :: syncWts         => sync_wts_dft_grid_t
        PROCEDURE, PASS :: compress        => compress_dft_grid_t
        PROCEDURE, PASS :: extend          => extend_dft_grid_t
    END TYPE dft_grid_t

!*******************************************************************************

!< Molecular grid information
    TYPE(dft_grid_t) :: molGrid

!< Array to store vertices for triangulation of sphere, for internal use
    REAL(KIND=fp), TARGET, ALLOCATABLE :: triangles(:,:,:,:)

    PRIVATE
    PUBLIC molGrid
    PUBLIC get_sorted_lebedev_pts
    PUBLIC sync_slices
    PUBLIC init_slices
    PUBLIC split_grid
    PUBLIC export_dft_grid
    PUBLIC tesselate_layer
    PUBLIC get_dft_grid_num_pts
    PUBLIC find_neighbours
    PUBLIC MAXDEPTH

 CONTAINS

!*******************************************************************************
! Legacy Fortran wrappers

!> @brief Export grid for use in legacy TD-DFT code in GAMESS
!> @param[out]   xyz       coordinates of grid points
!> @param[out]   w         weights of grid points
!> @param[out]   kcp       atoms, which the points belongs to
!> @param[out]   npts      number of nonzero point
!> @author Vladimir Mironov
 SUBROUTINE export_dft_grid(xyz, w, kcp, npts)
    REAL(KIND=fp), INTENT(OUT) :: xyz(3,*), w(*)
    INTEGER, INTENT(OUT) :: kcp(*)
    INTEGER, INTENT(OUT) :: npts
    REAL(KIND=fp) :: cutoff
    cutoff=1.0d-08/molGrid%nMolPts
    CALL molGrid%exportGrid(xyz, w, kcp, npts, cutoff)
 END SUBROUTINE export_dft_grid

!-------------------------------------------------------------------------------
!> @brief Get total number of grid points in DFT quadrature
!> @note Written for legacy TD-DFT code in GAMESS
!> @param[out]   npts      number of nonzero point
!> @author Vladimir Mironov
 SUBROUTINE get_dft_grid_num_pts(npts)
    INTEGER, INTENT(OUT) :: npts
    npts = molGrid%nMolPts
 END SUBROUTINE get_dft_grid_num_pts

!-------------------------------------------------------------------------------

!> @brief Syncrhonize slice information over DDI processes
!> @author Vladimir Mironov
 SUBROUTINE sync_slices
    CALL molGrid%sync
 END SUBROUTINE sync_slices

!-------------------------------------------------------------------------------

!> @brief Initialize/reset molecular grid
!> @param[in]   maxSlices   guess for maximum number of slices
!> @param[in]   nAt         number of atoms in a system
!> @param[in]   maxPtPerAt  maximum number of points per atom
!> @note  `maxSlices` need not to be precise. The arrays in `molGrid` will be
!>  automatically grown up if the actual number of slices exceeds it.
!> @author Vladimir Mironov
 SUBROUTINE init_slices(maxSlices, nAt, maxPtPerAt)
    INTEGER, INTENT(IN) :: maxSlices, nAt, maxPtPerAt
    CALL molGrid%reset(maxSlices,nAt,maxPtPerAt)
 END SUBROUTINE init_slices

!-------------------------------------------------------------------------------

!> @brief Find nearest neghbours of every atom and the corresponding distance
!>  The results are stored into `molGrid` dataset
!> @details This is needed for screening inner points in SSF variant of
!>  Becke's fuzzy cell method
!> @param[in]    rij       array of interatomic distances
!> @param[in]    nAt       number of atoms
!> @author Vladimir Mironov
 SUBROUTINE find_neighbours(rij, nAt)

    USE params, ONLY: dft_partfun
    USE mod_dft_partfunc, ONLY: partition_function
    REAL(KIND=fp), EXTERNAL :: braggRad
    INTEGER, INTENT(IN) :: nAt
    REAL(KIND=fp), INTENT(IN) :: rij(nAt,nAt)
    TYPE(partition_function) :: partFunc
    INTEGER :: i, j
    REAL(KIND=fp) :: distNN

    CALL partFunc%set(dft_partfun)

    DO j = 1, nAt
!       Find distance to the nearest neghbour of current atom
!       Also, tag non-real atoms present in the system
        distNN = HUGEFP
        DO i = 1, nAt
            molGrid%dummyAtom(i) = braggRad(i)==0.0_fp
            IF (molGrid%dummyAtom(i) .OR. i==j) CYCLE
            distNN = min(distNN, rij(i,j))
        END DO
        molGrid%rInner(j) = 0.5*distNN*(1.0-partFunc%limit)
    END DO

    !write(*,'(*(f16.10))') molGrid%rInner

 END SUBROUTINE


!-------------------------------------------------------------------------------

!> @brief Get spherical Lebedev grid of nReq number of points
!> @details If the grid is not yet computed, then compute it,
!>  sort according to pre-defined triangulation scheme and save
!>  the result for later use. If the same grid was already computed,
!>  fetch it from the in-memory storage
!> @author Vladimir Mironov
 SUBROUTINE get_sorted_lebedev_pts(xyz,w,nreq,mxdim)
    USE params, ONLY: dft_bfc_algo
    REAL(KIND=fp), TARGET, INTENT(OUT) :: xyz(*)
    REAL(KIND=fp), CONTIGUOUS, INTENT(OUT) :: w(:)
    INTEGER, INTENT(IN) :: nreq, mxdim
    LOGICAL, SAVE :: first = .TRUE.
    LOGICAL :: found
    REAL(KIND=fp), POINTER :: pxyz(:,:)
    INTEGER, ALLOCATABLE :: icnt(:,:)

    IF (first) THEN
        CALL init_gridlist
        first = .FALSE.
    END IF

    ALLOCATE(icnt(8*4**MAXDEPTH,MAXDEPTH+2))

    pxyz(1:nreq,1:3) => xyz(1:3*nreq)

    CALL get_grid(nreq,pxyz(:,1), pxyz(:,2), pxyz(:,3),w,icnt,found)

    IF (.NOT.found) THEN
        CALL lebptw(xyz, xyz, w, nreq, mxdim, .FALSE.)
        IF (dft_bfc_algo>=0) THEN
            CALL tesselate_layer(pxyz(:,1),pxyz(:,2),pxyz(:,3),w,nreq,icnt)
        END IF
        CALL set_grid(nreq,pxyz(:,1),pxyz(:,2),pxyz(:,3),w,icnt)
    END IF

 END SUBROUTINE

!-------------------------------------------------------------------------------

!> @brief Append slice data of an atom to `molGrid` arrays
!> @param[in]   idAtom      index of current atom
!> @param[in]   layers      meta-data for atom grid layers: radial grid, angular grid, and angular grid separation depth
!> @param[in]   nLayers     number of atom grid layers
!> @param[in]   rAtm        effective radius of current atom
!> @author Vladimir Mironov
 SUBROUTINE add_slices(idAtm, layers, nLayers, rAtm)
    USE mod_grid_storage, ONLY: ptrad => rad_grid
    INTEGER, INTENT(IN) :: idAtm, layers(:,:), nLayers
    REAL(KIND=fp), INTENT(IN) :: rAtm
    INTEGER :: i, j, ilen, idrad0, idrad1, idang0, idang1, nAngTr, nPtSlc, depth

    idrad1 = 0

    DO i = 1, nLayers
        ASSOCIATE ( &
                nCur      => molGrid%nSlices, &
                maxSlices => molGrid%maxSlices, &
                curAng    => grid_stack%grd(layers(2,i))%p, &
                nextRad   => layers(1,i), &
                depth0    => layers(3,i), &
                wtTopAtm  => molGrid%wt_top(idAtm) )

            idrad0 = idrad1 + 1
            idrad1 = nextRad
            depth = depth0
            IF ( ptrad(idrad0)<2.0 .AND. depth==-1) THEN
                IF (curAng%npts*(idrad1-idrad0+1)>=10*32) THEN
                    depth = 1
                ELSE IF ( curAng%npts*(idrad1-idrad0+1)>=10*8 ) THEN
                    depth = 0
                END IF
            END IF

            SELECT CASE (depth)
            CASE (-1)

                nPtSlc = (idrad1-idrad0+1)*curAng%nPts

                IF (nPtSlc == 0) CYCLE
                IF (nCur==maxSlices) CALL molGrid%extend

                nCur = nCur + 1

                CALL molGrid%setSlice( nCur, &
                        idrad0, idrad1-idrad0+1, &
                        curAng%idGrid, 1, curAng%nPts, &
                        wtTopAtm+1, &
                        idAtm, rAtm, 0)

                wtTopAtm = wtTopAtm + nPtSlc

            CASE (0:3)
                ilen = 4**depth
                idang1 = 0
                DO j = 1, 8*ilen
                    nAngTr = curAng%izones(j,depth+2)
                    IF (nAngTr==0) CYCLE
                    IF (nCur==maxSlices) CALL molGrid%extend

                    nCur = nCur + 1

                    idang0 = idang1 + 1
                    idang1 = idang0 + nAngTr - 1

                    nPtSlc = nAngTr*(idrad1-idrad0+1)

                CALL molGrid%setSlice( nCur, &
                         idrad0, idrad1-idrad0+1, &
                         curAng%idGrid, idAng0, nAngTr, &
                         wtTopAtm+1, &
                         idAtm, rAtm, 1)

                    wtTopAtm = wtTopAtm + nPtSlc
                END DO

            CASE DEFAULT
                nPtSlc = (idRad1-idRad0+1)*curAng%nPts
                IF (nPtSlc == 0) CYCLE
                IF (nCur==maxSlices) CALL molGrid%extend


                nCur = nCur + 1

                CALL molGrid%setSlice( nCur, &
                        idrad0, idrad1-idrad0+1, &
                        curAng%idGrid, 1, curAng%nPts, &
                        wtTopAtm+1, &
                        idAtm, rAtm, 2)

                wtTopAtm = wtTopAtm + nPtSlc

            END SELECT

        END ASSOCIATE
    END DO

 END SUBROUTINE

!-------------------------------------------------------------------------------

!> @brief Split atomic 3D grid on spacially localized clusters ("slices")
!>  and appended the data to `molGrid` dataset
!> @details First, decompose atomic grid on layers, depending on the pruning
!>  scheme and the distance from the nuclei. It is done in this procedure.
!>  Then, split layers on slices and append their data to the `molGrid`
!>  dataset. It is done by calling subroutine `addSlices`.
!> @param[in]   idAtm       index of current atom
!> @param[in]   nGrids      number of grids in pruning scheme
!> @param[in]   rAtm        effective radius of current atom
!> @param[in]   pruneRads   radii of the pruning scheme
!> @author Vladimir Mironov
 SUBROUTINE split_grid(idAtm, nGrids, rAtm, pruneRads)
    USE mod_grid_storage, ONLY: ptrad => rad_grid
    INTEGER,       INTENT(IN) :: idAtm, nGrids
    REAL(KIND=fp), INTENT(IN) :: rAtm, pruneRads(:)
    REAL(KIND=fp), PARAMETER :: EPS = tiny(0.0_fp)

    INTEGER,       PARAMETER :: depths(5) = [ -1,   0,   0,    1,   -2]
    REAL(KIND=fp), PARAMETER :: stepSz(5) = [1.0, 1.5, 4.0,  8.0,1.0e9]
    REAL(KIND=fp) :: limits(5)

    INTEGER :: i, nLayers, iSpl, iRMin, iRMax, iRNext, nRad!, k
    REAL(KIND=fp) :: rMin, rMax, rNext

    INTEGER :: layers(3,MXGRID*5)

    limits = [1.0, 2.0, 8.0, 12.0, 20.0]
    IF (molGrid%rInner(idAtm)/=0.0_fp) THEN
        limits(1) = molGrid%rInner(idAtm) + EPS
    END IF

    nRad = ubound(ptrad,1)
    nLayers = 0
    iRMax = 0
    iRNext = 0
    DO i = 1, nGrids
        iRMin = iRMax + 1
        iRMax = findl(pruneRads(i),ptRad,iRMin)

        IF (nGrids==1) iRMax = nRad

        rMin = rAtm*ptRad(iRMin)
        rMax = rAtm*ptRad(iRMax)

        iSpl = findl(rAtm*rMin,limits,1)

        rNext = rMin

        DO
            rNext  = min(rMax, rNext + stepSz(iSpl))
            iRNext = min(iRMax, findl(rNext/rAtm,ptRad,iRNext+1))

            nLayers = nLayers+1
            layers(1,nLayers) = iRNext
            layers(2,nLayers) = i
            layers(3,nLayers) = depths(iSpl)

            IF (rNext>=limits(iSpl)) iSpl = iSpl+1
            IF (iRNext==iRMax) EXIT
        END DO

    END DO

    CALL add_slices(idAtm, layers, nLayers, rAtm)

 CONTAINS

!> @brief Find the location of the largest element in
!>  sorted (in ascending order) real-valued array,
!>  which is smaller than the specified value. Return last array index
!>  if `ALL(X>=V)`
!> @note Linear search is used due to problem specifics:
!>  tiny arrays and sequential access
!> @param[in]   x       value to found
!> @param[in]   v       sorted real array
!> @param[in]   hint    the index from which the lookup is started
!> @author Vladimir Mironov
    FUNCTION findl(x,v,hint) RESULT(id)
        INTEGER :: id
        INTEGER :: hint
        REAL(KIND=fp) :: x, v(:)
        DO id = max(lbound(v,1),hint), ubound(v,1)
            IF (x<v(id)) RETURN
        END DO
        id = ubound(v,1)
    END FUNCTION

 END SUBROUTINE

!-------------------------------------------------------------------------------

!> @brief Split spherical grid on localized clusters of points
!> @param[inout] xs     grid point X coordinates
!> @param[inout] ys     grid point Y coordinates
!> @param[inout] zs     grid point Z coordinates
!> @param[inout] wt     grid point weights
!> @param[in]    npts   number of points
!> @param[out]   icnts  sizes of clusters
!> @author Vladimir Mironov
 SUBROUTINE tesselate_layer(xs,ys,zs,wt,npts,icnts)
    REAL(KIND=fp), INTENT(INOUT) :: xs(:), ys(:), zs(:)
    REAL(KIND=fp), INTENT(INOUT) :: wt(:)
    INTEGER, INTENT(IN) :: npts
    INTEGER, INTENT(OUT) :: icnts(:,-1:)
    REAL(KIND=fp), ALLOCATABLE :: tmp(:,:,:)
    LOGICAL, SAVE :: first = .TRUE.
    INTEGER :: i, j, itrng, itrngoct, ntrng, ntrngoct, m1, m2, m3, i0, i1, ilen
    REAL(KIND=fp), POINTER :: trng(:,:,:)

    ntrngoct = 4**MAXDEPTH
    ntrng = 8*ntrngoct

    IF (first) THEN
        ALLOCATE(triangles(3,3,ntrng,1:3))
        DO i = 1, MAXDEPTH+1
            CALL triangulate_sphere(triangles(:,:,:,i), i-1)
        END DO
        first = .FALSE.
    END IF

    trng => triangles(:,:,:,MAXDEPTH+1)

    ALLOCATE(tmp(size(xs,dim=1),4,ntrng))

    icnts(1:ntrng,:) = 0

    ASSOCIATE (icnt => icnts(:,MAXDEPTH) )
        DO i = 1, npts
            ASSOCIATE (x=>xs(i), y=>ys(i), z=>zs(i), w=>wt(i))
!               Find the octant of this point
                m1 = int(sign(0.5_fp,x)+0.5)
                m2 = int(sign(0.5_fp,y)+0.5)
                m3 = int(sign(0.5_fp,z)+0.5)
                itrngoct = 8 - (m1 + 2*m2 + 4*m3)
                DO itrng = ntrngoct*(itrngoct-1)+1, ntrngoct*itrngoct
!                   If point is on the boundary between two cells - first
!                   one gets it
                    IF (is_vector_inside_shape([x,y,z],trng(:,:,itrng))) EXIT
                END DO
                icnt(itrng) = icnt(itrng) + 1
                tmp(icnt(itrng),1,itrng) = x
                tmp(icnt(itrng),2,itrng) = y
                tmp(icnt(itrng),3,itrng) = z
                tmp(icnt(itrng),4,itrng) = w
            END ASSOCIATE
        END DO

        i0 = 0
        DO i = 1, ntrng
            i1 = i0 + icnt(i)
            xs(i0+1:i1) = tmp(1:icnt(i),1,i)
            ys(i0+1:i1) = tmp(1:icnt(i),2,i)
            zs(i0+1:i1) = tmp(1:icnt(i),3,i)
            wt(i0+1:i1) = tmp(1:icnt(i),4,i)
            i0 = i1
        END DO

    END ASSOCIATE

    DO i = MAXDEPTH-1, 0, -1
        ilen = 4**i
        DO j = 1, 8*ilen
            icnts(j,i) = sum(icnts((j-1)*4+1:j*4,i+1))
        END DO
    END DO

    icnts(1,-1) = npts

    DEALLOCATE(tmp)

 END SUBROUTINE

!-------------------------------------------------------------------------------

!> @brief Compute triangulation map of the sphere
!> @details This subroutine returns inscribed polyhedron with triangular
!>  faces and octahedral symmetry. Polyhedron faces are computed by recursive
!>  subdivision of the octahedron faces on four spherical
!>  equilateral triangles.
!> @param[in]   depth   depth of recursion
!> @param[out]  vect    vertices of triangles
!> @author Vladimir Mironov
 SUBROUTINE triangulate_sphere(vect, depth)
    INTEGER, INTENT(IN) :: depth
    REAL(KIND=fp), INTENT(OUT) :: vect(3,3,*)
    REAL(KIND=fp), DIMENSION(3) :: v1, v2, v3
    REAL(KIND=fp), PARAMETER :: xf(8) = [ 1,-1, 1,-1, 1,-1, 1,-1]
    REAL(KIND=fp), PARAMETER :: yf(8) = [ 1, 1,-1,-1, 1, 1,-1,-1]
    REAL(KIND=fp), PARAMETER :: zf(8) = [ 1, 1, 1, 1,-1,-1,-1,-1]
    INTEGER :: ip
    INTEGER :: i, ip0, ip1

    ip = 0

!   Fill in first octant
    v1 = [ 1.0,  0.0,  0.0]
    v2 = [ 0.0,  1.0,  0.0]
    v3 = [ 0.0,  0.0,  1.0]
    CALL subdivide_triangle(v1, v2, v3, vect, ip, depth)

!   Fill in other octants - use octaherdal symmetry of the polyhedron
    DO i = 2, 8
        ip0 = (i-1)*ip+1
        ip1 = i*ip
        vect(1,1:3,ip0:ip1) = xf(i)*vect(1,1:3,1:ip)
        vect(2,1:3,ip0:ip1) = yf(i)*vect(2,1:3,1:ip)
        vect(3,1:3,ip0:ip1) = zf(i)*vect(3,1:3,1:ip)
    END DO

 END SUBROUTINE

!> @brief Recursively subdivide spherical triangle on four equal triangles
!> @param[in]       v1      vertex of initial triangle
!> @param[in]       v2      vertex of initial triangle
!> @param[in]       v3      vertex of initial triangle
!> @param[inout]    res     stack of triangles; dim: (xyz,vertices,triangles)
!> @param[in]       ip      index of the last triangle on stack
!> @param[in]       depth   current depth of recursion
!> @author Vladimir Mironov
 RECURSIVE SUBROUTINE subdivide_triangle(v1, v2, v3, res, ip, depth)
   INTEGER, INTENT(IN) :: depth
   INTEGER, INTENT(INOUT) :: ip
   REAL(KIND=8), INTENT(IN) :: v1(3), v2(3), v3(3)
   REAL(KIND=8), INTENT(INOUT) :: res(3,3,*)
   REAL(KIND=8) :: c1(3), c2(3), c3(3)

   IF (depth<0.OR.depth>MAXDEPTH) THEN
       WRITE(*,*) 'DEPTH=', depth, ' IS .GT. MAXDEPTH=', MAXDEPTH
       CALL ABRT
   END IF

   IF (depth==0) THEN
!      add lowest level triangle to the result array
       ip = ip + 1
       res(:,1,ip) = v1
       res(:,2,ip) = v2
       res(:,3,ip) = v3
       RETURN
   END IF

!  find centers of arcs between input vectors
   c1 = v1 + v2
   c2 = v2 + v3
   c3 = v1 + v3
   c1 = c1/sqrt(sum(c1*c1))
   c2 = c2/sqrt(sum(c2*c2))
   c3 = c3/sqrt(sum(c3*c3))

!  subdivide the resulting four triangles
   CALL subdivide_triangle(v1, c1, c3, res, ip, depth-1)
   CALL subdivide_triangle(c1, v2, c2, res, ip, depth-1)
   CALL subdivide_triangle(c3, c2, v3, res, ip, depth-1)
   CALL subdivide_triangle(c1, c2, c3, res, ip, depth-1)

 END SUBROUTINE

!*******************************************************************************

!> @brief Clean up `molGrid` dataset and reallocate arrays if needed
!> @param[in]   maxSlices   guess for maximum number of slices
!> @param[in]   nAt         number of atoms in a system
!> @param[in]   maxPtPerAt  maximum number of points per atom
!> @author Vladimir Mironov
 SUBROUTINE reset_dft_grid_t(grid, maxSlices, nAt, maxPtPerAt)
    CLASS(dft_grid_t), INTENT(INOUT) :: grid
    INTEGER, INTENT(IN) :: maxSlices, nAt, maxPtPerAt

    IF (maxSlices<=0 .OR. nAt<=0 .OR. maxPtPerAt<=0) CALL abrt

    IF (grid%maxSlices<maxSlices) THEN
        grid%maxSlices = maxSlices
        IF (ALLOCATED(grid%idAng)) THEN
            DEALLOCATE(grid%idAng, grid%iAngStart, grid%nAngPts, &
                grid%iRadStart, grid%nRadPts, grid%nTotPts, &
                grid%idOrigin, grid%chunkType, grid%rAtm, &
                grid%wtStart, grid%isInner)
        END IF

        ALLOCATE(grid%idAng(maxSlices))
        ALLOCATE(grid%iAngStart(maxSlices))
        ALLOCATE(grid%nAngPts(maxSlices))

        ALLOCATE(grid%iRadStart(maxSlices))
        ALLOCATE(grid%nRadPts(maxSlices))

        ALLOCATE(grid%nTotPts(maxSlices))

        ALLOCATE(grid%idOrigin(maxSlices))
        ALLOCATE(grid%chunkType(maxSlices))
        ALLOCATE(grid%rAtm(maxSlices))
        ALLOCATE(grid%wtStart(maxSlices))
        ALLOCATE(grid%isInner(maxSlices))

    END IF

    IF (ALLOCATED(grid%totWts)) DEALLOCATE(grid%totWts)
    IF (ALLOCATED(grid%wt_top)) DEALLOCATE(grid%wt_top)
    ALLOCATE(grid%totWts(maxPtPerAt, nAt))
    ALLOCATE(grid%wt_top(nAt))

    IF (ALLOCATED(grid%rInner)) DEALLOCATE(grid%rInner)
    IF (ALLOCATED(grid%dummyAtom)) DEALLOCATE(grid%dummyAtom)
    ALLOCATE(grid%rInner(nAt))
    ALLOCATE(grid%dummyAtom(nAt))

    grid%idAng     = 0
    grid%iAngStart = 0
    grid%nAngPts   = 0
    grid%iRadStart = 0
    grid%nRadPts   = 0
    grid%nTotPts   = 0
    grid%idOrigin  = 0
    grid%chunkType = 0
    grid%wtStart   = 0
    grid%rAtm      = 0.0_fp
    grid%isInner   = 0
    grid%totWts    = 0.0_fp
    grid%wt_top    = 0

    grid%rInner = 0.0_fp
    grid%dummyAtom = .FALSE.

    grid%maxAtomPts = maxPtPerAt
    grid%maxSlicePts =  0
    grid%maxNRadTimesNAng = 0
    grid%nSlices = 0
    grid%nMolPts = 0

 END SUBROUTINE reset_dft_grid_t

!-------------------------------------------------------------------------------

!> @brief Syncrhonize grid information over DDI processes
!> @author Vladimir Mironov
 SUBROUTINE sync_dft_grid_t(grid)
    CLASS(dft_grid_t), INTENT(INOUT) :: grid
    COMMON /PAR   / me,master,nproc,ibtyp,iptim,goparr,dskwrk,maswrk
    INTEGER :: me,master,nproc,ibtyp,iptim
    LOGICAL :: dskwrk, maswrk, goparr

    INTEGER, PARAMETER :: TOTWT_SYNC    = 3515
    INTEGER, PARAMETER :: TOTPT_SYNC    = 3516
    INTEGER, PARAMETER :: INNERPT_SYNC  = 3517
    INTEGER, PARAMETER :: TOTMOLPT_SYNC = 3518
    INTEGER, PARAMETER :: MAXCHPT_SYNC  = 3519
    INTEGER, PARAMETER :: MAXCHPT2_SYNC = 3520

    INTEGER :: nMolPts(1)

    IF (allocated(grid%totWts) .AND.  goparr) THEN
        CALL ddi_gsumf(TOTWT_SYNC,grid%totWts,size(grid%totWts))
        CALL ddi_gsumi(TOTPT_SYNC,grid%nTotPts,grid%nSlices)
        CALL ddi_gsumi(INNERPT_SYNC,grid%isInner,grid%nSlices)
        nMolPts(1) = grid%nMolPts
        CALL ddi_gsumi(TOTMOLPT_SYNC,nMolPts(1),1)
        grid%nMolPts = nMolPts(1)
    END IF

    IF (goparr) THEN
        CALL ddi_gmaxi_scalar(grid%maxSlicePts)
        CALL ddi_gmaxi_scalar(grid%maxNRadTimesNAng)
    END IF

    CALL grid%compress

 END SUBROUTINE sync_dft_grid_t

!-------------------------------------------------------------------------------

!> @brief Syncrhonize number of points for every slices over DDI processes
!> @param[in]   nptCheck  number of grid points processed by current process
!> @note If the sum of `nptCheck` over processes is equal to current
!>  total number of molecular grid points, then no reduction will be done.
!> @author Vladimir Mironov
 SUBROUTINE sync_npts_dft_grid_t(grid, nptCheck)
    CLASS(dft_grid_t), INTENT(INOUT) :: grid
    INTEGER, INTENT(INOUT) :: nptCheck
    COMMON /PAR   / me,master,nproc,ibtyp,iptim,goparr,dskwrk,maswrk
    INTEGER :: me,master,nproc,ibtyp,iptim
    LOGICAL :: dskwrk, maswrk, goparr
    INTEGER, PARAMETER :: TOTPT_SYNC    = 3516
    INTEGER, PARAMETER :: TOTMOLPT_SYNC = 3518
    INTEGER, PARAMETER :: MAXCHPT_SYNC  = 3519

    INTEGER :: a_nptCheck(1)

    IF ( goparr ) THEN
        a_nptCheck(1) = nptCheck
        CALL ddi_gsumi(TOTMOLPT_SYNC,a_nptCheck(1),1)
        nptCheck = a_nptCheck(1)
    END IF

    IF (nptCheck/=grid%nMolPts) THEN
        IF (allocated(grid%nTotPts).AND.goparr) THEN
            CALL ddi_gsumi(TOTPT_SYNC,grid%nTotPts,grid%nSlices)
            grid%nTotPts = grid%nTotPts/nproc
        END IF

        CALL grid%compress

        grid%maxSlicePts = maxval(grid%nTotPts(1:grid%nSlices))
        grid%nMolPts = nptCheck
    END IF

 END SUBROUTINE sync_npts_dft_grid_t

!-------------------------------------------------------------------------------

!> @brief Syncrhonize quadrature weights over DDI processes
!> @author Vladimir Mironov
 SUBROUTINE sync_wts_dft_grid_t(grid)
    CLASS(dft_grid_t), INTENT(INOUT) :: grid
    COMMON /PAR   / me,master,nproc,ibtyp,iptim,goparr,dskwrk,maswrk
    INTEGER :: me,master,nproc,ibtyp,iptim
    LOGICAL :: dskwrk, maswrk, goparr
    INTEGER, PARAMETER :: WTS_SYNC = 3511
    REAL(KIND=fp) :: revNumProc

    IF (allocated(grid%totWts).AND.goparr) THEN
        revNumProc = 1.0_fp/nproc
        grid%totWts = grid%totWts*revNumProc
        CALL ddi_gsumf(WTS_SYNC,grid%totWts,size(grid%totWts))
    END IF

 END SUBROUTINE sync_wts_dft_grid_t

!-------------------------------------------------------------------------------

!> @brief Compress the grid for eaach atom
!> @note This is a very basic implementation and it will be changed in future
!> @author Vladimir Mironov
 SUBROUTINE compress_dft_grid_t(grid)
    CLASS(dft_grid_t), INTENT(INOUT) :: grid
    COMMON /PAR   / me,master,nproc,ibtyp,iptim,goparr,dskwrk,maswrk
    INTEGER :: me,master,nproc,ibtyp,iptim
    LOGICAL :: dskwrk, maswrk, goparr

    INTEGER :: i
    INTEGER :: iCur, nPts
    INTEGER, ALLOCATABLE :: iCurWt(:)

    ALLOCATE(iCurWt(ubound(grid%totWts,2)))

    iCur = 0
    iCurWt = 1
    DO i = 1, grid%nSlices
        IF (grid%nTotPts(i)>0) THEN
            iCur = iCur + 1
            grid%idAng(iCur)      = grid%idAng(i)
            grid%iAngStart(iCur)  = grid%iAngStart(i)
            grid%nAngPts(iCur)    = grid%nAngPts(i)
            grid%iRadStart(iCur)  = grid%iRadStart(i)
            grid%nRadPts(iCur)    = grid%nRadPts(i)
            grid%idOrigin(iCur)   = grid%idOrigin(i)
            grid%chunkType(iCur)  = grid%chunkType(i)
            grid%rAtm(iCur)       = grid%rAtm(i)
            grid%isInner(iCur)    = grid%isInner(i)

            grid%nTotPts(iCur)    = grid%nTotPts(i)

            ASSOCIATE ( totWts => grid%totWts, &
                        iAt => grid%idOrigin(i), wt0 => grid%wtStart(i), &
                        nAng => grid%nAngPts(i), nRad => grid%nRadPts(i) )

                nPts = nAng*nRad

                totWts(iCurWt(iAt):iCurWt(iAt)+nPts-1,iAt) = &
                    totWts(wt0:wt0+nPts-1,iAt)

                grid%wtStart(iCur) = iCurWt(iAt)
                iCurWt(iAt) = iCurWt(iAt) + nPts
            END ASSOCIATE
        END IF
    END DO

!   IF (maswrk) write(*,'("Removed",I8," empty slices out of", I8," total")') &
!            grid%nSlices-iCur, grid%nSlices

    grid%nSlices = iCur

    DEALLOCATE(iCurWt)

 END SUBROUTINE compress_dft_grid_t

!-------------------------------------------------------------------------------

!> @brief Extend arrays in molecular grid type
!> @author Vladimir Mironov
 SUBROUTINE extend_dft_grid_t(grid)
    CLASS(dft_grid_t), INTENT(INOUT) :: grid
    COMMON /PAR   / me,master,nproc,ibtyp,iptim,goparr,dskwrk,maswrk
    INTEGER :: me,master,nproc,ibtyp,iptim
    LOGICAL :: dskwrk, maswrk, goparr
    INTEGER :: maxSlices

    maxSlices = grid%maxSlices * 2

    CALL reallocate_int(grid%idAng,maxSlices)
    CALL reallocate_int(grid%iAngStart,maxSlices)
    CALL reallocate_int(grid%nAngPts,maxSlices)
    CALL reallocate_int(grid%iRadStart,maxSlices)
    CALL reallocate_int(grid%nRadPts,maxSlices)
    CALL reallocate_int(grid%nTotPts,maxSlices)
    CALL reallocate_int(grid%idOrigin,maxSlices)
    CALL reallocate_int(grid%chunkType,maxSlices)
    CALL reallocate_int(grid%wtStart,maxSlices)
    CALL reallocate_int(grid%isInner,maxSlices)

    CALL reallocate_real(grid%rAtm,maxSlices)

    grid%maxSlices = grid%maxSlices * 2

 CONTAINS
!> @brief Extend allocatable array of integers to a new size
!> @param[inout]  v         the array
!> @param[in]     newsz     new size of the array
!> @author Vladimir Mironov
     SUBROUTINE reallocate_int(v, newsz)
        INTEGER, ALLOCATABLE, INTENT(INOUT) :: v(:)
        INTEGER, INTENT(IN) :: newsz
        INTEGER, ALLOCATABLE :: nv(:)
        IF (allocated(v)) THEN
            ALLOCATE(nv(1:newsz), source=v)
            CALL move_alloc(from=nv, to=v)
        END IF
     END SUBROUTINE reallocate_int

!> @brief Extend allocatable array of integers to a new size
!> @param[inout]  v         the array
!> @param[in]     newsz     new size of the array
!> @author Vladimir Mironov
     SUBROUTINE reallocate_real(v, newsz)
        REAL(kind=fp), ALLOCATABLE, INTENT(INOUT) :: v(:)
        INTEGER, INTENT(IN) :: newsz
        REAL(kind=fp), ALLOCATABLE :: nv(:)
        IF (allocated(v)) THEN
            ALLOCATE(nv(1:newsz), source=v)
            CALL move_alloc(from=nv, to=v)
        END IF
     END SUBROUTINE reallocate_real

 END SUBROUTINE extend_dft_grid_t

!-------------------------------------------------------------------------------

!> @brief Get coordinates and weight for quadrature points, which
!>   belong to a slice
!> @param[in]     iSlice    index of a slice
!> @param[out]    xyzw      coordinates and weights
!> @author Vladimir Mironov
 SUBROUTINE getSliceData(grid, iSlice, xyzw)
    USE mod_grid_storage, ONLY: ptrad => rad_grid, wtRad => rad_wts
    CLASS(dft_grid_t), INTENT(IN) :: grid
    INTEGER, INTENT(IN) :: iSlice
    REAL(KIND=fp), CONTIGUOUS, INTENT(OUT) :: xyzw(:,:)

    REAL(KIND=fp), PARAMETER :: FOUR_PI = 4.0_fp*3.141592653589793238463_fp

    TYPE(grid_3d_t), POINTER :: curGrid
    INTEGER :: iAng, iRad, iPt
    REAL(KIND=fp) :: r1

    curGrid => get_grid_ptr(grid%idAng(iSlice))
    ASSOCIATE ( &
        rad       => grid%rAtm(iSlice), &
        isInner   => grid%isInner(iSlice), &
        wtStart   => grid%wtStart(iSlice)-1, &
        curAt     => grid%idOrigin(iSlice), &
        iAngStart => grid%iAngStart(iSlice), &
        iRadStart => grid%iRadStart(iSlice), &
        nAngPts   => grid%nAngPts(iSlice), &
        nRadPts   => grid%nRadPts(iSlice))

    ASSOCIATE ( &
        xAng      => curGrid%x(iAngStart:iAngStart+nAngPts-1), &
        yAng      => curGrid%y(iAngStart:iAngStart+nAngPts-1), &
        zAng      => curGrid%z(iAngStart:iAngStart+nAngPts-1), &
        wAng      => curGrid%w(iAngStart:iAngStart+nAngPts-1))

        DO iAng = 1, nAngPts
        DO iRad = 1, nRadPts
            r1 = rad*ptRad(iRadStart+iRad-1)
            iPt = (iAng-1)*nRadPts+iRad
            xyzw(iPt,1) = r1*xAng(iAng)
            xyzw(iPt,2) = r1*yAng(iAng)
            xyzw(iPt,3) = r1*zAng(iAng)
            IF ( isInner==0 ) THEN
                xyzw(iPt,4) = grid%totWts(wtStart+iPt,curAt)
            ELSE
                xyzw(iPt,4) = FOUR_PI*rad*rad*rad* &
                        wtRad(iRadStart+iRad-1)*wAng(iAng)
            END IF
        END DO
        END DO

    END ASSOCIATE
    END ASSOCIATE

 END SUBROUTINE

!-------------------------------------------------------------------------------

!> @brief Export grid for use in legacy TD-DFT code in GAMESS
!> @param[out]   xyz       coordinates of grid points
!> @param[out]   w         weights of grid points
!> @param[out]   kcp       atoms, which the points belongs to
!> @param[out]   npts      number of nonzero point
!> @param[in]    cutoff    cutoff to skip small weights
!> @author Vladimir Mironov
 SUBROUTINE exportGrid(grid, xyz, w, kcp, npts, cutoff)
    CLASS(dft_grid_t), INTENT(IN) :: grid
    REAL(KIND=fp), INTENT(OUT) :: xyz(3,*), w(*)
    INTEGER, INTENT(OUT) :: kcp(*)
    INTEGER, INTENT(OUT) :: npts
    REAL(KIND=fp), INTENT(IN) :: cutoff

    INTEGER :: iSlice, iPt, iPtOld, nPt, iAt

    REAL(KIND=fp), ALLOCATABLE :: xyzw(:,:)

    iPt = 0

!$omp parallel private(xyzw, iSlice, nPt, iPtOld, iAt)

    ALLOCATE(xyzw(molGrid%maxSlicePts,4))

!$omp do schedule(dynamic)
    DO iSlice = 1, molGrid%nSlices

        CALL grid%getSliceNonZero(cutoff, iSlice, xyzw, nPt)

        IF (nPt == 0 ) CYCLE

        !$omp atomic capture
            iPtOld = iPt
            iPt    = iPt + nPt
        !$omp end atomic

        iAt = molGrid%idOrigin(iSlice)

        xyz(1,iPtOld+1:iPtOld+nPt) = c(1,iAt) + xyzw(1:nPt,1)
        xyz(2,iPtOld+1:iPtOld+nPt) = c(2,iAt) + xyzw(1:nPt,2)
        xyz(3,iPtOld+1:iPtOld+nPt) = c(3,iAt) + xyzw(1:nPt,3)
        w(iPtOld+1:iPtOld+nPt)     =            xyzw(1:nPt,4)
        kcp(iPtOld+1:iPtOld+nPt)   = iAt

    END DO
!$omp end do
!$omp end parallel

    npts = iPt

 END SUBROUTINE

!-------------------------------------------------------------------------------

!> @brief Get grid points from a slice, which weights are larger
!>  than a cutoff.
!> @param[in]    cutoff    cutoff to skip small weights
!> @param[in]    iSlice    index of current slice
!> @param[out]   xyzw      coordinates and weights of grid points
!> @param[out]   nPt       number of nonzero point for slice
!> @author Vladimir Mironov
 SUBROUTINE getSliceNonZero(grid, cutoff, iSlice, xyzw, nPt)
    USE mod_grid_storage, ONLY: ptrad => rad_grid, wtRad => rad_wts
    CLASS(dft_grid_t), INTENT(IN) :: grid
    REAL(KIND=fp), INTENT(IN) :: cutoff
    INTEGER, INTENT(IN) :: iSlice
    REAL(KIND=fp), CONTIGUOUS, INTENT(OUT) :: xyzw(:,:)
    INTEGER, INTENT(OUT) :: nPt

    REAL(KIND=fp), PARAMETER :: FOUR_PI = 4.0_fp*3.141592653589793238463_fp

    TYPE(grid_3d_t), POINTER :: curGrid
    INTEGER :: iAng, iRad, iPt
    REAL(KIND=fp) :: r1, wtCur

    nPt = 0
    curGrid => get_grid_ptr(grid%idAng(iSlice))

    ASSOCIATE ( &
        rad       => grid%rAtm(iSlice), &
        isInner   => grid%isInner(iSlice), &
        wtStart   => grid%wtStart(iSlice)-1, &
        totWts    => grid%totWts, &
        curAt     => grid%idOrigin(iSlice), &
        iAngStart => grid%iAngStart(iSlice), &
        iRadStart => grid%iRadStart(iSlice), &
        nAngPts   => grid%nAngPts(iSlice), &
        nRadPts   => grid%nRadPts(iSlice))

    ASSOCIATE ( &
        xAng      => curGrid%x(iAngStart:iAngStart+nAngPts-1), &
        yAng      => curGrid%y(iAngStart:iAngStart+nAngPts-1), &
        zAng      => curGrid%z(iAngStart:iAngStart+nAngPts-1), &
        wAng      => curGrid%w(iAngStart:iAngStart+nAngPts-1))

        DO iAng = 1, nAngPts
        DO iRad = 1, nRadPts
            iPt = (iAng-1)*nRadPts+iRad
            IF ( isInner==0 ) THEN
                wtCur = totWts(wtStart+iPt,curAt)
            ELSE
                wtCur = FOUR_PI*rad*rad*rad* &
                        wtRad(iRadStart+iRad-1)*wAng(iAng)
            END IF
            IF (wtCur==0.0_fp) EXIT
            IF (wtCur<cutoff) CYCLE
            nPt = nPt + 1
            r1 = rad*ptRad(iRadStart+iRad-1)
            xyzw(nPt,1) = r1*xAng(iAng)
            xyzw(nPt,2) = r1*yAng(iAng)
            xyzw(nPt,3) = r1*zAng(iAng)
            xyzw(nPt,4) = wtCur
        END DO
        END DO
    END ASSOCIATE
    END ASSOCIATE

 END SUBROUTINE



!*******************************************************************************


!> @brief Set slice data
!> @param[in]    iSlice     index of current slice
!> @param[in]    iRadStart  index of the first point in radial grid array
!> @param[in]    nRadPts    number of radial grid points
!> @param[in]    idAng      index of angular grid
!> @param[in]    idAngStart index of the first point in angular grid array
!> @param[in]    nAngPts    number of angular grid points
!> @param[in]    wtStart    index of the first point in weights array
!> @param[in]    idAtm      index of atom which owns the slice
!> @param[in]    rAtm       effective radius of current atom
!> @param[in]    chunkType  type of chunk -- used for different processing
!> @author Vladimir Mironov
 SUBROUTINE setSlice(grid, iSlice, iRadStart, nRadPts, &
                 idAng, iAngStart, nAngPts, &
                 wtStart, idAtm, rAtm, chunkType)
    CLASS(dft_grid_t), INTENT(INOUT) :: grid

    INTEGER, INTENT(IN) :: iSlice
    INTEGER, INTENT(IN) :: iRadStart, nRadPts
    INTEGER, INTENT(IN) :: idAng, iAngStart, nAngPts
    INTEGER, INTENT(IN) :: wtStart
    INTEGER, INTENT(IN) :: idAtm, chunkType
    REAL(KIND=fp), INTENT(IN) :: rAtm

    grid%idAng(iSlice) = idAng

    grid%iRadStart(iSlice) = iRadStart
    grid%iAngStart(iSlice) = iAngStart

    grid%nRadPts(iSlice) = nRadPts
    grid%nAngPts(iSlice) = nAngPts
    grid%nTotPts(iSlice) = 0

    grid%idOrigin(iSlice) = idAtm
    grid%rAtm(iSlice) = rAtm

    grid%chunkType(iSlice) = chunkType
    grid%wtStart(iSlice) = wtStart

 END SUBROUTINE setSlice

!-------------------------------------------------------------------------------

!> @brief Compute average vector of two 3D vectors
!> @param[in]    vectors    input 3D vectors
!> @author Vladimir Mironov
 PURE FUNCTION vector_average(vectors) RESULT(avg)
    REAL(KIND=fp), INTENT(IN) :: vectors(:,:)
    REAL(KIND=fp) :: avg(3)

    REAL(KIND=fp) :: norm

    avg = sum(vectors(1:3,:), dim=2)

    norm = 1.0/sqrt(sum(avg**2))
    avg = avg * norm

 END FUNCTION vector_average

!-------------------------------------------------------------------------------

!> @brief Compute cross product of two 3D vectors
!> @param[in]    a    first vector
!> @param[in]    b    second vector
!> @author Vladimir Mironov
 PURE FUNCTION cross_product(a, b) RESULT(c)
     REAL(KIND=fp), INTENT(IN) :: a(:), b(:)
     REAL(KIND=fp) :: c(3)
     c(1) = a(2)*b(3)-a(3)*b(2)
     c(2) = a(3)*b(1)-a(1)*b(3)
     c(3) = a(1)*b(2)-a(2)*b(1)
 END FUNCTION

!-------------------------------------------------------------------------------

!> @brief Check if test vector crosses spherical polygon
!> @param[in]    test  test vector
!> @param[in]    vecs  vectors, defining the polygon
!> @author Vladimir Mironov
 PURE FUNCTION is_vector_inside_shape(test, vecs) RESULT(inside)
     REAL(KIND=fp), INTENT(IN) :: test(:), vecs(:,:)
     LOGICAL :: inside

     INTEGER :: i
     REAL(KIND=fp) :: sgn, sgnold, vectmp(3)

!!    Test vector must at least has same direction, as average vector
!     vectmp(:) = sum(vecs(:,:),dim=1)
!     IF (sum(test*vectmp)<0.0) THEN
!         inside = .FALSE.
!         RETURN
!     END IF

     inside = .TRUE.
     vectmp = cross_product(vecs(:,ubound(vecs,2)),vecs(:,lbound(vecs,2)))
     sgnold = sign(1.0_fp, sum(test(1:3)*vectmp(1:3)))

     DO i = lbound(vecs,2), ubound(vecs,2)-1
         vectmp = cross_product(vecs(:,i),vecs(:,i+1))
         sgn = sign(1.0_fp, sum(test(1:3)*vectmp(1:3)))
         IF (sgn/=sgnold) THEN
             inside = .FALSE.
             RETURN
         END IF
         sgnold = sgn
     END DO

 END FUNCTION

!*******************************************************************************

!> @brief Find maximum value of an integer over DDI processes
!> @param[inout]  ivar  value from current process
!> @author Vladimir Mironov
 SUBROUTINE ddi_gmaxi_scalar(ivar)
    INTEGER, INTENT(INOUT) :: ivar
    COMMON /PAR   / me,master,nproc,ibtyp,iptim,goparr,dskwrk,maswrk
    INTEGER :: me,master,nproc,ibtyp,iptim
    LOGICAL :: dskwrk, maswrk, goparr

    INTEGER, ALLOCATABLE, SAVE :: buf(:)
!   Reserve DDI descriptors 31200-31231
    INTEGER, PARAMETER :: DDI_GMAX_MAGIC = 31200
    INTEGER, PARAMETER :: DDI_GMAX_NDESC = 32
    INTEGER, SAVE :: isync = 0

    IF (goparr) THEN
!       Init buffer
        IF (allocated(buf)) THEN
            IF (size(buf)<nproc) THEN
                DEALLOCATE(buf)
            END IF
        END IF
        IF (.NOT.allocated(buf)) ALLOCATE(buf(nproc))
        buf = 0

!       Accumulate the data from all processes
        buf(me+1) = ivar
        CALL ddi_gsumi(DDI_GMAX_MAGIC+isync, buf, nproc)

!       Find max value
        ivar  = maxval(buf(:nproc))

!       Adjust DDI descriptor for future use
        isync = mod(isync+1, DDI_GMAX_NDESC)
    END IF

 END SUBROUTINE


END MODULE mod_dft_molgrid
