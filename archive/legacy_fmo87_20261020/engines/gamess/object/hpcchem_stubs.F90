
module hpc_chem
    implicit none
    logical :: enabled
    integer :: mode, print_mode
    logical :: use_scf, use_fock_build, use_cphf, use_tdhf, use_all
    double precision :: convergence_threshold, dynamic_threshold
    
    ! default values for user defined values 
    integer :: print_mode_d = 1
    double precision :: convergence_threshold_d = 1.0d-06
    double precision :: dynamic_threshold_d = 1.0d-10
    logical :: use_scf_d = .false.
    logical :: use_fock_build_d = .false.
    logical :: use_cphf_d = .false. 
    logical :: use_tdhf_d = .false.
    logical :: use_all_d = .false.
contains
    subroutine gms_hpcchem_initialize()
        implicit none
        character*256 env
        mode = 0
        call gms_getenv("GMS_HPCCHEM ", env)
        if (env.eq."true") then
            enabled = .true.
        end if
    end subroutine

    function gms_hpcchem() result(hpc_enabled)
        implicit none
        logical ::  hpc_enabled
        hpc_enabled = enabled
    end function

    function gms_hpc_mode() result(hpc_mode)
        implicit none
        integer :: hpc_mode
        hpc_mode = mode
    end function


    subroutine print_banner_libcchem()
        implicit none

    end subroutine print_banner_libcchem

    
end module hpc_chem

subroutine cchem_rhf_driver(density_matrix_gms,fock_matrix_gms, gafo_mo_coefficients, &
                              number_of_occupied_frozen_orbitals, &
                              number_of_frozen_orbitals, orbital_shift, total_energy, do_afo, number_of_electrons)
implicit none


double precision, target::  density_matrix_gms(*),fock_matrix_gms(*), gafo_mo_coefficients(*) 
double precision total_energy !resulting energy 
integer number_of_occupied_frozen_orbitals,  number_of_frozen_orbitals
integer number_of_electrons
real orbital_shift
logical :: do_afo

end subroutine cchem_rhf_driver


subroutine cchem_fock_build_driver(density_matrix_gms,fock_matrix_gms, number_of_electrons, do_tdhf_digestion)
implicit none


double precision, target::  density_matrix_gms(*),fock_matrix_gms(*) 
double precision total_energy !resulting energy 
integer number_of_electrons
integer do_tdhf_digestion



return
end subroutine cchem_fock_build_driver


