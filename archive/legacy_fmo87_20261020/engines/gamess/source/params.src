!*MODULE PARAMS
!>    @author  Vladimir Mironov
!
!>    @brief   (Dummy) module to store input parameters
!
!     REVISION HISTORY:
!>    @date _May, 2016_ Initial release
!>    @date _Jun, 2017_ Update for shared Fock matrix HF method
!
MODULE params

    IMPLICIT NONE

    INTEGER :: &
        intomp      !< Flag to select algorithm for calculating HF 2-electron contribution.
                    !< Possible values:
                    !< - 0 - MPI-only algorithm,
                    !< - 1 - MPI over I loop, OpenMP over JK loops,
                    !< - 2 - MPI over IJ loops, OpenMP over KL loops.
                    !<
                    !< Default:
                    !< - 0 - for MPI build,
                    !< - 1 - for MPI/OpenMP build.

    LOGICAL :: &
        shfock      !< .TRUE. for shared Fock matrix algorithm
                    !< Significantly reduces memory footprint for
                    !< highly parallel machines, but may affect scaling. Use with care. \n
                    !< Default: .FALSE.

    INTEGER :: &
        rad_grid_type !< type of the radial grid in DFT:
                      !< - 0 (default) - Euler-Maclaurin grid (Murray et al.)
                      !< - 1           - Log3 grid (Mura and Knowles)

    INTEGER :: &
        dft_bfc_algo  !< type of the Becke's fuzzy cell method
                      !< - 0 (default) - original BFC algorithm
                      !< - 1           - SSF-like

    LOGICAL :: &
        dft_wt_der    !< .TRUE. if quadrature weights derivative
                      !<  contribution to the nuclear gradient is needed
                      !< Weight derivatives are not always required,
                      !< especially if the fine grid is used

    CHARACTER(LEN=8) :: &
        dft_bfc_algo_code

    INTEGER :: &
        dft_partfun   !< partition function type in grid-based DFT
                      !< -  0 (default) - Becke's 4th degree polynomial
                      !< -  1           - SSF original polynomial
                      !< -  2           - Modified SSF ( erf(x/(1-x**2)) )
                      !< - -2           - Modified SSF (smoothstep-2, 5th order)
                      !<                  Warning: -2 is experimental option!
                      !< -  3           - Modified SSF (smoothstep-3, 7rd order)
                      !< -  4           - Modified SSF (smoothstep-4, 9th order)
                      !< -  5           - Modified SSF (smoothstep-5,11th order)
                      !< Note: Becke's polynomial with 1 iteration is actually
                      !< a smoothstep-1 polynomial (3*x**2 - 2*x**3)

    CHARACTER(LEN=8) :: &
        dft_partfun_code

    PUBLIC

END MODULE params
