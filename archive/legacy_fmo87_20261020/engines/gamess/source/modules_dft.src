! This file contains modules needed by DFT code
! MODULE metaGGA
!> @brief
!> @details Replaces the following common block:
!>          COMMON /METGGA/ NEEDTAU,PRTTAU
!> @author  Igor S. Gerasimov
!> @date    Aug, 2019 - Initial release
!
!> @params  needtau   - flag that toggles calculation of tau (\sum dot_product(\nabla \phi, \nabla \phi))
!> @params  needlapl  - flag that toggles calculation of lapl (\nabla^2 \rho)
!> @params  printtau  - flag that toggles debug output of tau
!> @params  printlapl - flag that toggles debug output of lapl
module metaGGA
    implicit none
    logical :: needtau
    logical :: printtau
    logical :: needlapl
    logical :: printlapl
end module metaGGA
! MODULE FUNCLIB
!> @brief
!> @details Replaces the following common block:
!>          COMMON /FUNLIB/ FUNCL,FUNFL
!> @author  Igor S. Gerasimov
!> @date    Aug, 2019 - Initial release
!
!> @param   funcl -
!> @param   funfl -
!
!> @todo add description for params
module funclib
    implicit none
    logical :: funcl
    logical :: funfl
end module funclib
! MODULE DFT_COEFFICIENTS
!> @brief
!> @details Replaces the following common blocks:
!>          COMMON /CBECKE/ CB88
!>          COMMON /CSTAX / CSLT
!>          COMMON /CLYPC / CLYP
!>          COMMON /CVWNC / CVWN1, CVWN3, CVWN5
!>          COMMON /CP86CF/ CP86
!>          COMMON /CPBE0 / CPBEX
!>          COMMON /CPW91 / CPW91C
!>          COMMON /CX3LYP/ CPW91
!>          COMMON /CPW91L/ CPW91LDA
!>          COMMON /CTPSSH/ CTPSSX, CTPSSMX, CTPSSC, CRTPSSX, CRTPSSC
!>          COMMON /SCPBES/ CPBESX, CPBESC
!>          COMMON /CPKZB / CPKZBX, CPKZBC
!>          COMMON /EXFACT/ COPTX, SCAOPTX
!>          COMMON /B97TYP/ B97TYP
!>          COMMON /SCAEDF/ CEDF1
!>          COMMON /SCAOPC/ COPC
!>          COMMON /SCLGIL/ CGILL
!>          COMMON /SCPFRE/ CPFREE
!>          COMMON /SCPZ81/ CPZ81C
!>          COMMON /SCRPBE/ CRPBE
!>          COMMON /SLPBEC/ CPBEC
!> @author  Igor S. Gerasimov
!> @date    Aug, 2019 - Initial release
!
!> @params  type_B97 - Toggles B97 parameters
module dft_coefficients
    implicit none
    integer          :: type_B97            ! type of B97 (selecting between any HCTHs, B97s, wB97s...
    double precision :: Becke88             ! Becke exchange coefficients in hybrid functionals
    double precision :: Slater              ! Slater exchange coefficients in hybrid functionals
    double precision :: LYP                 ! Lee-Yang-Parr correlation coefficients in hybrid functionals
    double precision :: VWN1                ! VWN1 correlation coefficients in hybrid functionals
    double precision :: VWN3                ! VWN3 correlation coefficients in hybrid functionals
    double precision :: VWN5                ! VWN5 correlation coefficients in hybrid functionals
    double precision :: Perdew86            ! P86 correlation coefficients in hybrid functionals
    double precision :: PBE_Exchange        ! PBE exchange coefficients in hybrid functionals
    double precision :: PBE_Correlation     ! PBE correlation coefficients in hybrid functionals
    double precision :: PW91_Exchange       ! PW91 exchange coefficients in hybrid functionals
    double precision :: PW91_Correlation    ! PW91 correlation coefficients in hybrid functionals
    double precision :: PW91_Local          ! PW91 LDA (PW92) correlation coefficients in hybrid functionals
    double precision :: TPSS_Exchange       ! TPSS exchange coefficients in hybrid functionals
    double precision :: TPSSmod_Exchange    ! modified TPSS exchange coefficients in hybrid functionals
    double precision :: TPSS_Correlation    ! TPSS correlation coefficients in hybrid functionals
    double precision :: revTPSS_Exchange    ! revised TPSS exchange coefficients in hybrid functionals
    double precision :: revTPSS_Correlation ! revised TPSS correlation coefficients in hybrid functionals
    double precision :: PBEsol_Exchange     ! PBEsol exchange coefficients in hybrid functionals
    double precision :: PBEsol_Correlation  ! PBEsol correlation coefficients in hybrid functionals
    double precision :: PKZB_Exchange       ! PKZB exchange coefficients in hybrid functionals
    double precision :: PKZB_Correlation    ! PKZB correlation coefficients in hybrid functionals
    double precision :: OPT                 ! OPT exchange coefficients in hybrid functionals
    double precision :: Scale_OPT           ! Internal OPT scale coefficient
    double precision :: EDF_1               ! EDF1 coefficients in hybrid functionals
    double precision :: OP                  !
    double precision :: Gill                ! GILL exchange coefficients in hybrid functionals
    double precision :: ParamFree           !
    double precision :: PZ81                ! PZ81 correlation coefficients in hybrid functionals
    double precision :: revPBE              ! revised PBE exchange coefficients in hybrid functionals
    double precision :: RPBE                ! RPBE exchange coefficients in hybrid functionals
end module dft_coefficients
!MODULE DFTEXC
!> @brief
!> @details Replaces the following common block:
!>          COMMON /DFTEXC/ PI,QOP,NEXFG,NCORFG,NPFFG,NXCFG
!> @author Igor S. Gerasimov
!> @date   Aug, 2019 - Initial release
!>
!> @param  QOP    -
!> @param  NEXFG  - ID of exchange functional
!> @param  NCORFG - ID of correlation functional
!> @param  NPFFG  - looks like ID of parameter free functional
!> @param  NXCFG  - looks like ID of exchange-correlation functional
module dftexc
    implicit none
    double precision :: QOP
    integer :: NEXFG
    integer :: NCORFG
    integer :: NPFFG
    integer :: NXCFG
end module dftexc
!MODULE LRCDFT
!> @brief
!> @details Replaces the following common block:
!>          COMMON /NLRC  / LCFLAG,EMU,EMU2,LRFILE
!> @author Igor S. Gerasimov
!> @date   Aug, 2019 - Initial release
!>
!> @param  lcflag - switch long-range corrected scheme for DFT exchange functional
!> @param  emu    - mu coefficient
!> @param  emu2   - usually mu^2 (not sure that it is true everywhere)
!> @param  lrfile -
module lrcdft
     implicit none
     logical :: lcflag
     double precision :: emu
     double precision :: emu2
     integer :: lrfile
end module lrcdft
!MODULE CAMDFT
!> @brief
!> @details Replaces the following common block:
!>          COMMON /DFTCAM/ ALPHAC,BETAC,CAMMU,CAMVWN,CAMLYP,CAMFLAG
!> @author Igor S. Gerasimov
!> @date   Aug, 2019 - Initial release
!>
!> @param  cam_alpha - alpha coefficient
!> @param  cam_beta  - beta coefficient
!> @param  cam_mu    - mu coefficient
!> @param  cam_vwn5  - coefficient of VWN5 correlation functional (for CAM-B3LYP and similar)
!> @param  cam_lyp   - coefficient of LYP  correlation functional (for CAM-B3LYP and similar)
!> @param  camflag   - switch Coulomb-Attenuating Method (CAM-) \frac{1}{r_{12}} = \frac{1-[cam_{alpha}+cam_{beta}*\erf(cam_{mu}*r_{12})]}{r_{12}} + \frac{cam_{alpha}+cam_{beta}*\erf(cam_{mu}*r_{12})}{r_{12}}
module camdft
    implicit none
    double precision :: cam_alpha
    double precision :: cam_beta
    double precision :: cam_mu
    double precision :: cam_vwn5
    double precision :: cam_lyp
    logical :: camflag
end module camdft
!MODULE DHDFT
!> @brief
!> @details Replaces the following common block:
!>          COMMON /DFTDH / CHF,CMP2,C2S,C2T,DHFUNC
!> @author Igor S. Gerasimov
!> @date   Aug, 2019 - Initial release
!
! @params CHF    - coefficient of HF exact exchange
! @params CMP2   - coefficient of non-local MP2 correlation
! @params C2S    - coefficient of non-local MP2 same-spin correlation
! @params C2T    - coefficient of non-local MP2 opposite-spin correlation
! @params DHFUNC - logical flag of using DH-DFT functionals
module dhdft
    implicit none
    double precision :: CHF !this variable equivalent to DFTTYP(3)
    double precision :: CMP2
    double precision :: C2S
    double precision :: C2T
    logical :: DHFUNC
end module dhdft
!MODULE XDERIX
!> @brief
!> @details Replaces the following common block:
!>          COMMON /XDERIX/ ...
!> @author Igor S. Gerasimov
!> @date   Aug, 2019 - Initial release
!
! @params
module XDERIX
    IMPLICIT NONE
    INTEGER, DIMENSION(3), PARAMETER :: NXDIM = [6, 18, 38] !The count of derivatives for functional. NCDIM(1) means first derivatives, NCDIM(2) means second...
    INTEGER, PARAMETER :: KRA     =  1                      ! dE/dR_A
    INTEGER, PARAMETER :: KRB     =  2                      ! dE/dR_B
    INTEGER, PARAMETER :: KGA     =  3                      ! dE/dG_AA
    INTEGER, PARAMETER :: KGB     =  4                      ! dE/dG_BB
    INTEGER, PARAMETER :: KTA     =  5                      ! dE/dT_A
    INTEGER, PARAMETER :: KTB     =  6                      ! dE/dT_B
    INTEGER, PARAMETER :: KRARA   =  7                      ! d2E/(dR_A  dR_A )
    INTEGER, PARAMETER :: KRAGA   =  8                      ! d2E/(dR_A  dG_AA)
    INTEGER, PARAMETER :: KRATA   =  9                      ! d2E/(dR_A  dT_A )
    INTEGER, PARAMETER :: KGAGA   = 10                      ! d2E/(dG_AA dG_AA)
    INTEGER, PARAMETER :: KGATA   = 11                      ! d2E/(dG_AA dT_A )
    INTEGER, PARAMETER :: KTATA   = 12                      ! d2E/(dT_A  dT_A )
    INTEGER, PARAMETER :: KRBRB   = 13                      ! d2E/(dR_B  dR_B )
    INTEGER, PARAMETER :: KRBGB   = 14                      ! d2E/(dR_B  dG_BB)
    INTEGER, PARAMETER :: KRBTB   = 15                      ! d2E/(dR_B  dT_B )
    INTEGER, PARAMETER :: KGBGB   = 16                      ! d2E/(dG_BB dG_BB)
    INTEGER, PARAMETER :: KGBTB   = 17                      ! d2E/(dG_BB dT_B )
    INTEGER, PARAMETER :: KTBTB   = 18                      ! d2E/(dT_B  dT_B )
    INTEGER, PARAMETER :: KRARARA = 19                      ! d3E/(dR_A  dR_A  dR_A )
    INTEGER, PARAMETER :: KRARAGA = 20                      ! d3E/(dR_A  dR_A  dG_AA)
    INTEGER, PARAMETER :: KRARATA = 21                      ! d3E/(dR_A  dR_A  dT_A )
    INTEGER, PARAMETER :: KRAGAGA = 22                      ! d3E/(dR_A  dG_AA dG_AA)
    INTEGER, PARAMETER :: KRAGATA = 23                      ! d3E/(dR_A  dG_AA dT_A )
    INTEGER, PARAMETER :: KRATATA = 24                      ! d3E/(dR_A  dT_A  dT_A )
    INTEGER, PARAMETER :: KGAGAGA = 25                      ! d3E/(dG_AA dG_AA dG_AA)
    INTEGER, PARAMETER :: KGAGATA = 26                      ! d3E/(dG_AA dG_AA dT_A )
    INTEGER, PARAMETER :: KGATATA = 27                      ! d3E/(dG_AA dT_A  dT_A )
    INTEGER, PARAMETER :: KTATATA = 28                      ! d3E/(dT_A  dT_A  dT_A )
    INTEGER, PARAMETER :: KRBRBRB = 29                      ! d3E/(dR_B  dR_B  dR_B )
    INTEGER, PARAMETER :: KRBRBGB = 30                      ! d3E/(dR_B  dR_B  dG_BB)
    INTEGER, PARAMETER :: KRBRBTB = 31                      ! d3E/(dR_B  dR_B  dT_B )
    INTEGER, PARAMETER :: KRBGBGB = 32                      ! d3E/(dR_B  dG_BB dG_BB)
    INTEGER, PARAMETER :: KRBGBTB = 33                      ! d3E/(dR_B  dG_BB dT_B )
    INTEGER, PARAMETER :: KRBTBTB = 34                      ! d3E/(dR_B  dT_B  dT_B )
    INTEGER, PARAMETER :: KGBGBGB = 35                      ! d3E/(dG_BB dG_BB dG_BB)
    INTEGER, PARAMETER :: KGBGBTB = 36                      ! d3E/(dG_BB dG_BB dT_B )
    INTEGER, PARAMETER :: KGBTBTB = 37                      ! d3E/(dG_BB dT_B  dT_B )
    INTEGER, PARAMETER :: KTBTBTB = 38                      ! d3E/(dT_B  dT_B  dT_B )
end module XDERIX
!MODULE CDERIX
!> @brief
!> @details Replaces the following common block:
!>          COMMON /CDERIX/ ...
!> @author Igor S. Gerasimov
!> @date   Aug, 2019 - Initial release
!
! @params
module CDERIX
    IMPLICIT NONE
    INTEGER, DIMENSION(3), PARAMETER :: NCDIM = [7, 35, 119] !The count of derivatives for functional. NCDIM(1) means first derivatives, NCDIM(2) means second...
    INTEGER, PARAMETER :: IRA     =   1                      ! dE/dR_A
    INTEGER, PARAMETER :: IRB     =   2                      ! dE/dR_B
    INTEGER, PARAMETER :: IGA     =   3                      ! dE/dG_AA
    INTEGER, PARAMETER :: IGB     =   4                      ! dE/dG_BB
    INTEGER, PARAMETER :: IGC     =   5                      ! dE/dG_AB
    INTEGER, PARAMETER :: ITA     =   6                      ! dE/dT_A
    INTEGER, PARAMETER :: ITB     =   7                      ! dE/dT_B
    INTEGER, PARAMETER :: IRARA   =   8                      ! d2E/(dR_A  dR_A )
    INTEGER, PARAMETER :: IRARB   =   9                      ! d2E/(dR_A  dR_B )
    INTEGER, PARAMETER :: IRAGA   =  10                      ! d2E/(dR_A  dG_AA)
    INTEGER, PARAMETER :: IRAGB   =  11                      ! d2E/(dR_A  dG_BB)
    INTEGER, PARAMETER :: IRAGC   =  12                      ! d2E/(dR_A  dG_AB)
    INTEGER, PARAMETER :: IRATA   =  13                      ! d2E/(dR_A  dT_A )
    INTEGER, PARAMETER :: IRATB   =  14                      ! d2E/(dR_A  dT_B )
    INTEGER, PARAMETER :: IRBRB   =  15                      ! d2E/(dR_B  dR_B )
    INTEGER, PARAMETER :: IRBGA   =  16                      ! d2E/(dR_B  dG_AA)
    INTEGER, PARAMETER :: IRBGB   =  17                      ! d2E/(dR_B  dG_BB)
    INTEGER, PARAMETER :: IRBGC   =  18                      ! d2E/(dR_B  dG_AB)
    INTEGER, PARAMETER :: IRBTA   =  19                      ! d2E/(dR_B  dT_A )
    INTEGER, PARAMETER :: IRBTB   =  20                      ! d2E/(dR_B  dT_B )
    INTEGER, PARAMETER :: IGAGA   =  21                      ! d2E/(dG_AA dG_AA)
    INTEGER, PARAMETER :: IGAGB   =  22                      ! d2E/(dG_AA dG_BB)
    INTEGER, PARAMETER :: IGAGC   =  23                      ! d2E/(dG_AA dG_AB)
    INTEGER, PARAMETER :: IGATA   =  24                      ! d2E/(dG_AA dT_A )
    INTEGER, PARAMETER :: IGATB   =  25                      ! d2E/(dG_AA dT_B )
    INTEGER, PARAMETER :: IGBGB   =  26                      ! d2E/(dG_BB dG_BB)
    INTEGER, PARAMETER :: IGBGC   =  27                      ! d2E/(dG_BB dG_AB)
    INTEGER, PARAMETER :: IGBTA   =  28                      ! d2E/(dG_BB dT_A )
    INTEGER, PARAMETER :: IGBTB   =  29                      ! d2E/(dG_BB dT_B )
    INTEGER, PARAMETER :: IGCGC   =  30                      ! d2E/(dG_AB dG_AB)
    INTEGER, PARAMETER :: IGCTA   =  31                      ! d2E/(dG_AB dT_A )
    INTEGER, PARAMETER :: IGCTB   =  32                      ! d2E/(dG_AB dT_B )
    INTEGER, PARAMETER :: ITATA   =  33                      ! d2E/(dT_A  dT_A )
    INTEGER, PARAMETER :: ITATB   =  34                      ! d2E/(dT_A  dT_B )
    INTEGER, PARAMETER :: ITBTB   =  35                      ! d2E/(dT_B  dT_B )
    INTEGER, PARAMETER :: IRARARA =  36                      ! d3E/(dR_A  dR_A  dR_A )
    INTEGER, PARAMETER :: IRARARB =  37
    INTEGER, PARAMETER :: IRARAGA =  38
    INTEGER, PARAMETER :: IRARAGB =  39
    INTEGER, PARAMETER :: IRARAGC =  40
    INTEGER, PARAMETER :: IRARATA =  41
    INTEGER, PARAMETER :: IRARATB =  42
    INTEGER, PARAMETER :: IRARBRB =  43
    INTEGER, PARAMETER :: IRARBGA =  44
    INTEGER, PARAMETER :: IRARBGB =  45
    INTEGER, PARAMETER :: IRARBGC =  46
    INTEGER, PARAMETER :: IRARBTA =  47
    INTEGER, PARAMETER :: IRARBTB =  48
    INTEGER, PARAMETER :: IRAGAGA =  49
    INTEGER, PARAMETER :: IRAGAGB =  50
    INTEGER, PARAMETER :: IRAGAGC =  51
    INTEGER, PARAMETER :: IRAGATA =  52
    INTEGER, PARAMETER :: IRAGATB =  53
    INTEGER, PARAMETER :: IRAGBGB =  54
    INTEGER, PARAMETER :: IRAGBGC =  55
    INTEGER, PARAMETER :: IRAGBTA =  56
    INTEGER, PARAMETER :: IRAGBTB =  57
    INTEGER, PARAMETER :: IRAGCGC =  58
    INTEGER, PARAMETER :: IRAGCTA =  59
    INTEGER, PARAMETER :: IRAGCTB =  60
    INTEGER, PARAMETER :: IRATATA =  61
    INTEGER, PARAMETER :: IRATATB =  62
    INTEGER, PARAMETER :: IRATBTB =  63
    INTEGER, PARAMETER :: IRBRBRB =  64
    INTEGER, PARAMETER :: IRBRBGA =  65
    INTEGER, PARAMETER :: IRBRBGB =  66
    INTEGER, PARAMETER :: IRBRBGC =  67
    INTEGER, PARAMETER :: IRBRBTA =  68
    INTEGER, PARAMETER :: IRBRBTB =  69
    INTEGER, PARAMETER :: IRBGAGA =  70
    INTEGER, PARAMETER :: IRBGAGB =  71
    INTEGER, PARAMETER :: IRBGAGC =  72
    INTEGER, PARAMETER :: IRBGATA =  73
    INTEGER, PARAMETER :: IRBGATB =  74
    INTEGER, PARAMETER :: IRBGBGB =  75
    INTEGER, PARAMETER :: IRBGBGC =  76
    INTEGER, PARAMETER :: IRBGBTA =  77
    INTEGER, PARAMETER :: IRBGBTB =  78
    INTEGER, PARAMETER :: IRBGCGC =  79
    INTEGER, PARAMETER :: IRBGCTA =  80
    INTEGER, PARAMETER :: IRBGCTB =  81
    INTEGER, PARAMETER :: IRBTATA =  82
    INTEGER, PARAMETER :: IRBTATB =  83
    INTEGER, PARAMETER :: IRBTBTB =  84
    INTEGER, PARAMETER :: IGAGAGA =  85                      ! d3E/(dG_AA dG_AA dG_AA)
    INTEGER, PARAMETER :: IGAGAGB =  86                      ! d3E/(dG_AA dG_AA dG_BB)
    INTEGER, PARAMETER :: IGAGAGC =  87                      ! d3E/(dG_AA dG_AA dG_AB)
    INTEGER, PARAMETER :: IGAGATA =  88
    INTEGER, PARAMETER :: IGAGATB =  89
    INTEGER, PARAMETER :: IGAGBGB =  90                      ! d3E/(dG_AA dG_BB dG_BB)
    INTEGER, PARAMETER :: IGAGBGC =  91                      ! d3E/(dG_AA dG_BB dG_AB)
    INTEGER, PARAMETER :: IGAGBTA =  92
    INTEGER, PARAMETER :: IGAGBTB =  93
    INTEGER, PARAMETER :: IGAGCGC =  94                      ! d3E/(dG_AA dG_AB dG_AB)
    INTEGER, PARAMETER :: IGAGCTA =  95
    INTEGER, PARAMETER :: IGAGCTB =  96
    INTEGER, PARAMETER :: IGATATA =  97
    INTEGER, PARAMETER :: IGATATB =  98
    INTEGER, PARAMETER :: IGATBTB =  99
    INTEGER, PARAMETER :: IGBGBGB = 100                      ! d3E/(dG_BB dG_BB dG_BB)
    INTEGER, PARAMETER :: IGBGBGC = 101                      ! d3E/(dG_BB dG_BB dG_AB)
    INTEGER, PARAMETER :: IGBGBTA = 102
    INTEGER, PARAMETER :: IGBGBTB = 103
    INTEGER, PARAMETER :: IGBGCGC = 104                      ! d3E/(dG_BB dG_AB dG_AB)
    INTEGER, PARAMETER :: IGBGCTA = 105
    INTEGER, PARAMETER :: IGBGCTB = 106
    INTEGER, PARAMETER :: IGBTATA = 107
    INTEGER, PARAMETER :: IGBTATB = 108
    INTEGER, PARAMETER :: IGBTBTB = 109
    INTEGER, PARAMETER :: IGCGCGC = 110                      ! d3E/(dG_AB dG_AB dG_AB)
    INTEGER, PARAMETER :: IGCGCTA = 111
    INTEGER, PARAMETER :: IGCGCTB = 112
    INTEGER, PARAMETER :: IGCTATA = 113
    INTEGER, PARAMETER :: IGCTATB = 114
    INTEGER, PARAMETER :: IGCTBTB = 115
    INTEGER, PARAMETER :: ITATATA = 116                      ! d3E/(dT_A  dT_A  dT_A )
    INTEGER, PARAMETER :: ITATATB = 117                      ! d3E/(dT_A  dT_A  dT_B )
    INTEGER, PARAMETER :: ITATBTB = 118                      ! d3E/(dT_A  dT_B  dT_B )
    INTEGER, PARAMETER :: ITBTBTB = 119                      ! d3E/(dT_B  dT_B  dT_B )
end module CDERIX
!> @brief  Stores the variables EXCHR,CORCHR
!>         in subroutine INPGDFT
!>
!> @author Christian Friedl
!>
!> @date   Jan, 2023 - Initial release
!>
module modinpgdft
    IMPLICIT NONE
    CHARACTER(LEN=10), SAVE :: EXCHR
    CHARACTER(LEN=10), SAVE :: CORCHR
end module
!> @brief  Stores the variable PFKIN
!>         in subroutine INPPFREE
!>
!> @author Christian Friedl
!>
!> @date   Jan, 2023 - Initial release
!>
module modinppfree
    IMPLICIT NONE
    CHARACTER(LEN=10), SAVE :: PFKIN
end module
