!*MODULE libxc
!> @brief  The head of libxc driver
!> @author Igor S. Gerasimov
!> @date   July, 2019 - Initial release -
!> @date   July, 2021 Making internal subroutines private
module libxc
  use prec, only: fp
  implicit none
  private
  logical :: use_libxc = .FALSE.
  public :: use_libxc
  public :: libxc_input, libxc_calc, libxc_TD_2D, libxc_TD_3D, libxc_destroy
contains
  !
  !*MODULE libxc DECK libxc_input
  !> @brief  setting up of using libxc functionals
  !>         (Analog of INPGDFT)
  !> @author Igor S. Gerasimov
  !> @date   July, 2019 - Initial release -
  !> @date   July, 2021 Using messages module
  !>                    Adding optional arguments
  !> @params read_unit  (in,optional) logical unit of input file
  !> @params functional (in,optional) functional's name
  subroutine libxc_input(read_unit, functional)
    use messages, only: show_message, WITH_ABORT
    integer, intent(in), optional :: read_unit
    character(len=*),    optional :: functional
    call show_message("Please, compile GAMESS-US with LibXC", WITH_ABORT)
  end subroutine libxc_input
  !
  !> @author Igor S. Gerasimov
  !> @date   July, 2021 Using messages module
  !>                    double precision -> real(kind=fp)
  !> @params WEIGHT - weight of point (FTOTWT)
  !> @params RHOA   - alpha-density at point
  !> @params RHOB   - beta-density at point
  !> @params GRD__  - gradients of electron density
  !> @params TAU__  - local kinetic energies of electron density
  !> @params XALPHA - LDA exchange energy (full energy)
  !> @params XGRD   - GGA exchange energy (always is zero)
  !> @params ECF    - correlation energy  (always is zero)
  !> @params VXC_1  - gradient of energy by spin-density
  !> @params DUM__  - gradient of energy by normed gradient of  electron density
  !> @params DMGG_  - gradient of energy by local kinetic energies
  subroutine libxc_calc(WEIGHT,                                    &
                        RHOA, RHOB,                                &
                        GRDAA, GRDAB, GRDBB,                       &
                        GRADXA,GRADYA,GRADZA,GRADXB,GRADYB,GRADZB, &
                        TAUXA,TAUYA,TAUZA,TAUXB,TAUYB,TAUZB,       &
                        XALPHA,XGRD,ECF,                           &
                        VXCA1,VXCB1,                               &
                        DUMAX,DUMAY,DUMAZ,DUMBX,DUMBY,DUMBZ,       &
                        DMGGA, DMGGB)
    use messages, only: show_message, WITH_ABORT
    real(kind=fp), intent(in)  :: WEIGHT
    real(kind=fp), intent(in)  :: RHOA, RHOB
    real(kind=fp), intent(in)  :: GRDAA, GRDAB, GRDBB
    real(kind=fp), intent(in)  :: GRADXA,GRADYA,GRADZA,GRADXB,GRADYB,GRADZB
    real(kind=fp), intent(in)  :: TAUXA,TAUYA,TAUZA,TAUXB,TAUYB,TAUZB
    real(kind=fp), intent(out) :: XALPHA,XGRD,ECF
    real(kind=fp), intent(out) :: VXCA1,VXCB1
    real(kind=fp), intent(out) :: DUMAX,DUMAY,DUMAZ,DUMBX,DUMBY,DUMBZ
    real(kind=fp), intent(out) :: DMGGA, DMGGB
    call show_message("Please, compile GAMESS-US with LibXC", WITH_ABORT)
    XALPHA=0
    XGRD=0
    ECF=0
    VXCA1=0
    VXCB1=0
    DUMAX=0
    DUMAY=0
    DUMAZ=0
    DUMBX=0
    DUMBY=0
    DUMBZ=0
    DMGGA=0
    DMGGB=0
  end subroutine libxc_calc
  !
  !*MODULE libxc DECK libxc_input
  !> @brief  Destroy internal variables of functional
  !> @author Igor S. Gerasimov
  !> @date   Dec, 2020 - Initial release -
  !> @date   July, 2021 Using messages module
  subroutine libxc_destroy()
    use messages, only: show_message, WITH_ABORT
    call show_message("Please, compile GAMESS-US with LibXC", WITH_ABORT)
  end subroutine libxc_destroy
  !
  !*MODULE libxc DECK libxc_calc_2D
  !> @brief  Performing TD-DFT calculations using libxc
  !>         (Analog of any TDFUNC subroutine)
  !> @author Igor S. Gerasimov
  !> @date   Sep, 2019 - Initial release -
  !> @date   Jul, 2021 Using messages module
  !>                   double precision -> real(kind=fp)
  !> @params WaveType    - Type of wavefunction: restricted (1) or unrestricted (2). integer value
  !> @params NGridPoints - count of grid points (also maxgrd)
  !> @params NPoints     - count of points
  !> @params DerivLength - Size of array of derivatives energy by rho/grad/tau
  !> @params EXC         - total energy
  !> @params EnergyPoint - energies at each point
  !> @params DerivsPoint - derivatives energy at each point
  !> @params dRho        - rho and gradients rho at each point
  !> @params dTau        - tau of rho at each point
  !> @params weight      - weight of point at each point
  !> @params MinRho      - except point where rho less that MinRho
  !> @params StartPoint  - The first point, which needed to work (Needs for parallel runs, I think)
  !> @params EndPoint    - The last point, which needed to work (Needs for parallel runs, I think)
  subroutine libxc_TD_2D(WaveType, NGridPoints, NPoints, DerivLength, EXC, EnergyPoint, DerivsPoint, dRho, dTau, weight, MinRho, StartPoint, EndPoint)
    use messages, only: show_message, WITH_ABORT
    integer, intent(in) :: WaveType
    integer, intent(in) :: NGridPoints
    integer, intent(in) :: NPoints
    integer, intent(in) :: DerivLength
    real(kind=fp), intent(out) :: EXC
    real(kind=fp), dimension(NPoints),                intent(out) :: EnergyPoint
    real(kind=fp), dimension(NPoints,DerivLength),    intent(out) :: DerivsPoint
    real(kind=fp), dimension(NGridPoints,4,WaveType), intent(in)  :: dRho
    real(kind=fp), dimension(NGridPoints,WaveType),   intent(in)  :: dTau
    real(kind=fp), dimension(NGridPoints),            intent(in)  :: weight
    real(kind=fp), intent(in) :: MinRho
    integer, intent(in) :: StartPoint
    integer, intent(in) :: EndPoint
    call show_message("Please, compile GAMESS-US with LibXC", WITH_ABORT)
    EnergyPoint=0
    DerivsPoint=0
    EXC=0
  end subroutine libxc_TD_2D
  !
  !*MODULE libxc DECK libxc_TD_3D
  !> @brief  Performing TD-DFT calculations using libxc
  !>         (Analog of any TDFUNC subroutine)
  !> @detail For determining a shifting like (j * 3 + 3) see libxc/scripts/maple2c_new.pl @derivatives.
  !> @author Igor S. Gerasimov
  !> @date   Sep, 2019 - Initial release -
  !> @date   Jul, 2021 Using messages module
  !>                   double precision -> real(kind=fp)
  !> @params WaveType    - Type of wavefunction: restricted (1) or unrestricted (2). integer value
  !> @params NGridPoints - count of grid points (also maxgrd)
  !> @params NPoints     - count of points
  !> @params DerivLength - Size of array of derivatives energy by rho/grad/tau
  !> @params EXC         - total energy
  !> @params EnergyPoint - energies at each point
  !> @params DerivsPoint - derivatives energy at each point
  !> @params dRho        - rho and gradients rho at each point
  !> @params dTau        - tau of rho at each point
  !> @params weight      - weight of point at each point
  !> @params MinRho      - except point where rho less that MinRho
  !> @params StartPoint  - The first point, which needed to work (Needs for parallel runs, I think)
  !> @params EndPoint    - The last point, which needed to work (Needs for parallel runs, I think)
  subroutine libxc_TD_3D(WaveType, NGridPoints, NPoints, DerivLength, EXC, EnergyPoint, DerivsPoint, dRho, dTau, weight, MinRho, StartPoint, EndPoint)
    use messages, only: show_message, WITH_ABORT
    integer, intent(in) :: WaveType
    integer, intent(in) :: NGridPoints
    integer, intent(in) :: NPoints
    integer, intent(in) :: DerivLength
    real(kind=fp), intent(out) :: EXC
    real(kind=fp), dimension(NPoints),                intent(out) :: EnergyPoint
    real(kind=fp), dimension(NPoints,DerivLength),    intent(out) :: DerivsPoint
    real(kind=fp), dimension(NGridPoints,4,WaveType), intent(in)  :: dRho
    real(kind=fp), dimension(NGridPoints,WaveType),   intent(in)  :: dTau
    real(kind=fp), dimension(NGridPoints),            intent(in)  :: weight
    real(kind=fp), intent(in) :: MinRho
    integer, intent(in) :: StartPoint
    integer, intent(in) :: EndPoint
    call show_message("Please, compile GAMESS-US with LibXC", WITH_ABORT)
    EXC=0
    EnergyPoint=0
    DerivsPoint=0
  end subroutine libxc_TD_3D
end module libxc
