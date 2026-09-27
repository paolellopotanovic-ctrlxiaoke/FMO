! MODULE MX_LIMITS
!>    @author  Vladimir Mironov
!
!>    @brief   Contains parameters scattered throughout all
!>             of the code that define sizes of static
!>             memory arrays
!
!     REVISION HISTORY:
!>    @date _Jan, 2017_ Initial release
!
!>    @date Jun 2022, Christian Friedl
!>    - Added parameter tinker_maxatm
!
MODULE mx_limits

    IMPLICIT NONE

    INTEGER, PARAMETER :: &
        mxatm=2000, &               !< Max. number of atoms in a job
        mxgtot=20000, &             !< Max. number of primitive Gaussians in basis set
        mxsh=5000, &                !< Max. number of shells in basis set
        mxao=8192, &                !< Max. number of atomic orbitals in basis set
        mxrt=100, &                 !< Max. number of states (roots of Hamiltonian matrix) for post-HF methods
        mxnoro=250, &               !< Max. number of orbital rotation pairs, omitted from MCSCF optimization
        mxgrid=10, &                !< Max. number of grid-based DFT grids
        mxgridtyp=10, &             !< Max. number of atom grid types for DFT grids
        mxang=7, &                  !< Max. angular momentum
        mxang2=mxang*mxang, &
        mxang3=mxang2*mxang, &
        mxcart=28,          &       !< Max. Cartesian elements of an i function
        mxdim=mxcart*mxcart, &
        mxfrz=15, &                 !< Max. number of frozen MCSCF orbitals
        mxgsh=30, &                 !< Max. degree of contraction
        mxg2=mxgsh*mxgsh, &
        maxsh=mxang*(mxang+1)*(mxang+2)/6, & !< Total number of Cartesian components for angular momentums from 0 to 6 (S to I)
        mxabc=6000, &               !< Max. number of COSMO cavity surface segments
        maxden=25*mxatm, &          !< Max. dimension of COSMO array for multipoles
! AAMBS-oriented values
        lenilist=19, &              !< Max. number of subshells
! EFP-oriented values
        mxfrg=1050, &               !< Max. number or EFP framgents in a job
        mxfgpt=12000, &             !< Max. number of EFP expansion points
        !mxfgpt=23000, &             !< Max. number of EFP expansion points
        mxdfg=5, &                  !< Max. number of EFP fragment types
        mxdppt=mxfrg*mxdfg*12, &    !< Max. number of EFP dynamic polarizable points
        mxpt=2000,&                 !< Max. number of EFP polarizable points (?)
        mxifrq=12,&                 !< Max. number of EFP imaginary frequencies
        mxpairs=mxdfg*(mxdfg+1)/2,& !< Max. number of type pairs of EFP potentials
        mxgefp=4000,&               !< Max. number of gaussian basis functions for EFP
        mxshef=1000,&               !< Max. number of shells for EFP
        mxcpuefp=2048,&             !< Max. number of cpus for EFP (?)
! EFMO-oriented values
       mxefmopts=100,&              !< Max. number of multipole expansion points per fragment (multipole expansion points are typically centered on the atoms and possibly on the bond centroids.)
       mxefmoppts=305,&             !< Max. number of polarizable points per fragment (polarizable points are points where the polarization expansion is done, so typically this is the number of LMOs per fragment, which is <= the number of occupied orbitals per fragment)
       mxnefmopts=27,&              !< Max. number of stored points per fragment per multipole expansion point (that is, at each point, there are monopole, dipole, quadrupole, etc.)
       mxnefmoppts=130,&             !< Max. number of stored points per fragment per polarizable point (that is, at each point, there are elements of the polarizability tensor, etc.)
! NEO-oriented values
       mxneo=30, &                  !< Max. number of quantum nuclei for a NEO calculation
! RELATIVISTIC MODEL CORE POTENTIALS (MCP)
       mxcs=1000, &                 !< Max. number of valence shells
       mxnt=40, &                   !< Max. number of ???
       mxc=21, &                    !< Max. number of ???
       mxt1=10, &                   !< Max. number of ???
       mxio=202, &                  !< Max. number of ???
       mxelm=30, &                  !< Max. number of ???
       mxpos=100, &                 !< Max. number of ???
       mxm=12, &                    !< Max. number of ???
       mxt=10, &                    !< Max. number of ???
! TINKER
       tinker_maxatm=12000          !< Contains parameter MAXATM which controls the dimensions in /atmtyp/ and /atmst/ used to interface between GAMESS and Tinker. Note that the value to be used must match exactly the problem size selected in Tinker's sizes.i include file.

END MODULE mx_limits
