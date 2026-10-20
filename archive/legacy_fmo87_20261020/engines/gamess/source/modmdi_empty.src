! MODULE modmdi
!>    @author  Taylor Barnes (MolSSI)
!
!>    @brief   Dummy implementation of an MDI interface in GAMESS
!
!>    @detail
!              This module enables compilation of GAMESS without MDI.
!
!     REVISION HISTORY:
!>    @date Aug, 2021 - Initial release
!
!>    @param   USE_MDI - Flag whether this is an MDI run
MODULE modmdi
    use, intrinsic::iso_c_binding, only : C_NULL_PTR, C_PTR, C_INT

    IMPLICIT NONE

    CHARACTER(len=:), ALLOCATABLE :: MDI_OPTIONS

    ! Flag whether to use MDI
    LOGICAL :: USE_MDI = .FALSE.

    ! Flag whether the latest energy evaluation is still valid
    LOGICAL :: energy_current = .FALSE.

    ! Flag whether the latest force evaluation is still valid
    LOGICAL :: forces_current = .FALSE.

    ! Maximum length of the MDI options
    INTEGER, PARAMETER :: MDI_OPTIONS_LEN = 1024

    ! MDI communicator for inter-code communication
    INTEGER(C_INT) :: mdi_comm

    PRIVATE
    PUBLIC :: USE_MDI, MDI_RUN

CONTAINS


!*MODULE MODMDI  *DECK MDI_READ_OPTIONS
!>    @brief   Read options used to initialze MDI
!>    @detail
!          Dummy version of the MDI_READ_OPTIONS subroutine in modmdi.src.
!>    @author Taylor A. Barnes
!>    @date   Aug, 2021 - Initial release
!
    subroutine MDI_READ_OPTIONS()
        ! Variables for getting the command-line options
        INTEGER :: istatus

        ALLOCATE(character(len=MDI_OPTIONS_LEN) :: mdi_options)

        call GET_ENVIRONMENT_VARIABLE("GAMESS_MDI_OPTIONS",mdi_options, STATUS=istatus)

        ! Return if GAMESS_MDI_OPTIONS is not set
        IF ( istatus .ne. 0 ) RETURN

        ! Return if GAMESS_MDI_OPTIONS is empty
        IF ( TRIM(mdi_options) .eq. "" ) RETURN

        USE_MDI = .TRUE.

    end subroutine MDI_READ_OPTIONS

!*MODULE MODMDI_EMPTY  *DECK MDI_RUN
!>    @brief   Run as an MDI engine, if the -mdi command-line option is present
!>    @detail
!          Dummy version of the MDI_RUN subroutine in modmdi.src.
!>    @author Taylor A. Barnes
!>    @date   Aug, 2021 - Initial release
!
    subroutine MDI_RUN()
        ! GET VARIABLES FROM OLD CODE
        integer :: IR, IW, IP, IJK, IJKT, IDAF, NAV, IODA(950)
        COMMON /IOFILE/ IR, IW, IP, IJK, IJKT, IDAF, NAV, IODA
        logical :: GOPARR,DSKWRK,MASWRK
        integer :: ME,MASTER,NPROC,IBTYP,IPTIM
        common /par   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK

        CALL MDI_READ_OPTIONS()

        IF ( USE_MDI ) THEN
            IF (MASWRK) THEN
                WRITE(IW,"(A)")'Found MDI option, but the code was not compiled without MDI support'
                CALL ABRT()
            END IF
        END IF
    end subroutine MDI_RUN


END MODULE modmdi
