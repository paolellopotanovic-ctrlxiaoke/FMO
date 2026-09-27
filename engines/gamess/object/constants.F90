! MODULE constants
!> @brief   This module contains popular constants
!> @author  Igor S. Gerasimov
!> @date    Aug, 2019 - Initial release
!> @date    Jan, 2020 - Add sqrt_PI
!> @date    Nov, 2020 - George Schoendorff - added many more
!> @date    Oct, 2021 - Add quarter, half, two, three, five, six, seven,
!>                      eight, ten
!> @date    Aug, 2023 - George Schoendorff
!>                      - Added nine, eleven, thirteen, fifteen, sqrt3,
!>                        sqrt5, sqrt7, sqrt9, sqrt11, jx, ix, jy, iy,
!>                        jz, iz, ijz, ijy, ijz
!> @date    Dec, 2023 - Peng Xu
!>                      - add twelve
!> @date    May, 2023 - George Schoendorff
!>                      - Added the periodic table as ATOM array
!>
!> @date    Mar, 2024 - George Schoendorff
!>                      - Added ATMLAB (array of atomic symbols)
!
!> @params  zero    ! 0.0
!> @papams  quarter ! 0.25
!> @params  half    ! 0.5
!> @params  one     ! 1.0
!> @params  two     ! 2.0
!> @params  three   ! 3.0
!> @params  four    ! 4.0
!> @params  five    ! 5.0
!> @params  six     ! 6.0
!> @params  seven   ! 7.0
!> @params  eight   ! 8.0
!> @params  nine    ! 9.0
!> @params  ten     ! 10.0
!> @params  eleven  ! 11.0
!> @params  twelve  ! 12.0
!> @params  thirteen! 13.0
!> @params  fifteen ! 15.0
!> @params  sqrt3   ! 1.73205080756888 ...
!> @params  sqrt5   ! 2.23606797749979 ...
!> @params  sqrt7   ! 2.64575131106459 ...
!> @params  sqrt9   ! 3.0
!> @params  sqrt11  ! 3.31662479035540 ...
!> @params  pi      ! 3.1415926...
!> @params  sqrt_pi ! 1.7724538...
!> @params  pi32    ! 5.56832799683170...
!> @params  pi212   ! 1.1283791670955
!> @params  rln10   ! 2.30258509299405...
!> @params  pt5     ! 0.5
!> @params  pt75    ! 0.75
!> @params  tm1     ! 0.1
!> @params  tm3     ! 0.001
!> @params  tm5     ! 0.00001
!> @params  tm6     ! 0.000001
!> @params  tm8     ! 0.00000001
!> @params  tm10    ! 0.0000000001
!> @params  tm13    ! 0.0000000000001
!> @params  tm14    ! 0.00000000000001
!> @params  debye   ! 2.541766
!> @params  tokcal  ! 627.509541D+00
!> @params  ix, jx, iy, jy, iz, yz ! Tensors of i and j indices used
!>                                   when computing one electron
!>                                   integrals
!> @params  ijx     ! 1 + power of X in Cartesian basis functions
!> @params  ijy     ! 1 + power of Y in Cartesian basis functions
!> @params  ijz     ! 1 + power of Z in Cartesian basis functions
!
module constants
    implicit none
    double precision, parameter :: ZERO = 0.0D+00
    double precision, parameter :: QUARTER = 0.25D+00
    double precision, parameter :: HALF = 0.5D+00
    double precision, parameter :: ONE = 1.0D+00
    double precision, parameter :: TWO = 2.0D+00
    double precision, parameter :: THREE = 3.0D+00
    double precision, parameter :: FOUR = 4.0D+00
    double precision, parameter :: FIVE = 5.0D+00
    double precision, parameter :: SIX = 6.0D+00
    double precision, parameter :: SEVEN = 7.0D+00
    double precision, parameter :: NINE = 9.0D+00
    double precision, parameter :: EIGHT = 8.0D+00
    double precision, parameter :: TEN = 1.0D+01
    double precision, parameter :: ELEVEN = 1.1D+01
    double precision, parameter :: TWELVE = 1.2D+01
    double precision, parameter :: THIRTEEN = 1.3D+01
    double precision, parameter :: FIFTEEN = 1.5D+01
    double precision, parameter :: SQRT3 = sqrt(three)
    double precision, parameter :: SQRT5 = sqrt(five)
    double precision, parameter :: SQRT7 = sqrt(seven)
    double precision, parameter :: SQRT9 = 3.0D+00
    double precision, parameter :: SQRT11 = (eleven)
    double precision, parameter :: PI = FOUR*ATAN(ONE)
    double precision, parameter :: sqrt_PI = SQRT(PI)
    double precision, parameter :: PI32 = PI*SQRT_PI
    double precision, parameter :: pi212 = 1.1283791670955D+00
    double precision, parameter :: rln10 = log(ten)
    double precision, parameter :: PT5  = 5.0D-01
    double precision, parameter :: PT75 = 7.5D-01
    double precision, parameter :: TM1  = 1.0D-01
    double precision, parameter :: TM3  = 1.0D-03
    double precision, parameter :: TM5  = 1.0D-05
    double precision, parameter :: TM6  = 1.0D-06
    double precision, parameter :: TM8  = 1.0D-08
    double precision, parameter :: TM10 = 1.0D-10
    double precision, parameter :: TM13 = 1.0D-13
    double precision, parameter :: TM14 = 1.0D-14
    double precision, parameter :: DEBYE = 2.541766
    double precision, parameter :: TOKCAL = 627.509541D+00

    integer, parameter :: jx(84) = [ &
     0, 1, 0, 0, 2, 0, 0, 1, 1, 0, &
     3, 0, 0, 2, 2, 1, 0, 1, 0, 1, &
     4, 0, 0, 3, 3, 1, 0, 1, 0, 2, &
     2, 0, 2, 1, 1, &
     5, 0, 0, 4, 4, 1, 0, 1, 0, 3, & 
     3, 2, 0, 2, 0, 3, 1, 1, 2, 2, 1, &
     6, 0, 0, 5, 5, 1, 0, 1, 0, 4, 4, &
     2, 0, 2, 0, 4, 1, 1, 3, 3, &
     0, 3, 3, 2, 1, 2, 1, 2 ]
    integer, parameter :: ix(84) = [ &
     1, 8, 1, 1,15, 1, 1, 8, 8, 1, &
     22, 1, 1,15,15, 8, 1, 8, 1, 8, &
     29, 1, 1,22,22, 8, 1, 8, 1,15, &
     15, 1,15, 8, 8, &
     36, 1, 1,29,29, 8, 1, 8, 1,22, &
     22,15, 1,15, 1,22, 8, 8,15,15, 8, &
     43, 1, 1,36,36, 8, 1, 8, 1,29, &
     29,15, 1,15, 1,29, 8, 8,22,22, &
     1,22,22,15, 8,15, 8,15 ]
    integer, parameter :: jy(84) = [ &
     0, 0, 1, 0, 0, 2, 0, 1, 0, 1, &
     0, 3, 0, 1, 0, 2, 2, 0, 1, 1, &
     0, 4, 0, 1, 0, 3, 3, 0, 1, 2, &
     0, 2, 1, 2, 1, &
     0, 5, 0, 1, 0, 4, 4, 0, 1, 2, &
     0, 3, 3, 0, 2, 1, 3, 1, 2, 1, 2, &
     0, 6, 0, 1, 0, 5, 5, 0, 1, 2, &
     0, 4, 4, 0, 2, 1, 4, 1, 3, 0, &
     3, 2, 1, 3, 3, 1, 2, 2 ]
    integer, parameter :: iy(84) = [ &
     1, 1, 8, 1, 1,15, 1, 8, 1, 8, &
     1,22, 1, 8, 1,15,15, 1, 8, 8, &
     1,29, 1, 8, 1,22,22, 1, 8,15, &
     1,15, 8,15, 8, &
     1,36, 1, 8, 1,29,29, 1, 8,15, &
     1,22,22, 1,15, 8,22, 8,15, 8, 15, &
     1,43, 1, 8, 1,36,36, 1, 8,15, &
     1,29,29, 1,15, 8,29, 8,22, 1, &
     22,15, 8,22,22, 8,15,15 ]
    integer, parameter :: jz(84) = [ &
     0, 0, 0, 1, 0, 0, 2, 0, 1, 1, &
     0, 0, 3, 0, 1, 0, 1, 2, 2, 1, &
     0, 0, 4, 0, 1, 0, 1, 3, 3, 0, &
     2, 2, 1, 1, 2, &
     0, 0, 5, 0, 1, 0, 1, 4, 4, 0, &
     2, 0, 2, 3, 3, 1, 1, 3, 1, 2, 2, &
     0, 0, 6, 0, 1, 0, 1, 5, 5, 0, &
     2, 0, 2, 4, 4, 1, 1, 4, 0, 3, &
     3, 1, 2, 1, 2, 3, 3, 2 ]
    integer, parameter :: iz(84) = [ &
     1, 1, 1, 8, 1, 1,15, 1, 8, 8, &
     1, 1,22, 1, 8, 1, 8,15,15, 8, &
     1, 1,29, 1, 8, 1, 8,22,22, 1, &
     15,15, 8, 8,15, &
     1, 1,36, 1, 8, 1, 8,29,29, 1, &
     15, 1,15,22,22, 8, 8,22, 8,15, 15, &
     1, 1,43, 1, 8, 1, 8,36,36, 1, &
     15, 1,15,29,29, 8, 8,29, 1,22, &
     22, 8,15, 8,15,22,22,15 ]
    integer, parameter :: ijx(84) = [ &
      1, 2, 1, 1, 3, 1, 1, 2, 2, 1, &
      4, 1, 1, 3, 3, 2, 1, 2, 1, 2, &
      5, 1, 1, 4, 4, 2, 1, 2, 1, 3, &
      3, 1, 3, 2, 2, 6, 1, 1, 5, 5, &
      2, 1, 2, 1, 4, 4, 3, 1, 3, 1, &
      4, 2, 2, 3, 3, 2, 7, 1, 1, 6, &
      6, 2, 1, 2, 1, 5, 5, 3, 1, 3, &
      1, 5, 2, 2, 4, 4, 1, 4, 4, 3, &
      2, 3, 2, 3 ]
    integer, parameter :: ijy(84) = [ &
      1, 1, 2, 1, 1, 3, 1, 2, 1, 2, &
      1, 4, 1, 2, 1, 3, 3, 1, 2, 2, &
      1, 5, 1, 2, 1, 4, 4, 1, 2, 3, &
      1, 3, 2, 3, 2, 1, 6, 1, 2, 1, &
      5, 5, 1, 2, 3, 1, 4, 4, 1, 3, &
      2, 4, 2, 3, 2, 3, 1, 7, 1, 2, &
      1, 6, 6, 1, 2, 3, 1, 5, 5, 1, &
      3, 2, 5, 2, 4, 1, 4, 3, 2, 4, &
      4, 2, 3, 3 ]
    integer, parameter :: ijz(84) = [ &
      1, 1, 1, 2, 1, 1, 3, 1, 2, 2, &
      1, 1, 4, 1, 2, 1, 2, 3, 3, 2, &
      1, 1, 5, 1, 2, 1, 2, 4, 4, 1, &
      3, 3, 2, 2, 3, 1, 1, 6, 1, 2, &
      1, 2, 5, 5, 1, 3, 1, 3, 4, 4, &
      2, 2, 4, 2, 3, 3, 1, 1, 7, 1, &
      2, 1, 2, 6, 6, 1, 3, 1, 3, 5, &
      5, 2, 2, 5, 1, 4, 4, 2, 3, 2, &
      3, 4, 4, 3 ]

    character(len=2), parameter :: atom(118) = [ &
      'H ','He','Li','Be','B ','C ','N ','O ','F ','Ne', &
      'Na','Mg','Al','Si','P ','S ','Cl','Ar','K ','Ca', &
      'Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn', &
      'Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y ','Zr', &
      'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn', &
      'Sb','Te','I ','Xe','Cs','Ba','La','Ce','Pr','Nd', &
      'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb', &
      'Lu','Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg', &
      'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th', &
      'Pa','U ','Np','Pu','Am','Cm','Bk','Cf','Es','Fm', &
      'Md','No','Lr','Rf','Db','Sg','Bh','Hs','Mt','Ds', &
      'Rg','Cn','Nh','Fl','Mc','Lv','Ts','Og']

    character(len=4), parameter :: atmlab(120) = [ &
      'H   ','HE  ','LI  ','BE  ','B   ','C   ', &
      'N   ','O   ','F   ','NE  ','NA  ','MG  ', &
      'AL  ','SI  ','P   ','S   ','CL  ','AR  ', &
      'K   ','CA  ','SC  ','TI  ','V   ','CR  ', &
      'MN  ','FE  ','CO  ','NI  ','CU  ','ZN  ', &
      'GA  ','GE  ','AS  ','SE  ','BR  ','KR  ', &
      'RB  ','SR  ','Y   ','ZR  ','NB  ','MO  ', &
      'TC  ','RU  ','RH  ','PD  ','AG  ','CD  ', &
      'IN  ','SN  ','SB  ','TE  ','I   ','XE  ', &
      'CS  ','BA  ','LA  ','CE  ','PR  ','ND  ', &
      'PM  ','SM  ','EU  ','GD  ','TB  ','DY  ', &
      'HO  ','ER  ','TM  ','YB  ','LU  ','HF  ', &
      'TA  ','W   ','RE  ','OS  ','IR  ','PT  ', &
      'AU  ','HG  ','TL  ','PB  ','BI  ','PO  ', &
      'AT  ','RN  ','FR  ','RA  ','AC  ','TH  ', &
      'PA  ','U   ','NP  ','PU  ','AM  ','CM  ', &
      'BK  ','CF  ','ES  ','FM  ','MD  ','NO  ', &
      'LR  ','RF  ','DB  ','SG  ','BH  ','HS  ', &
      'MT  ','DS  ','RG  ','CN  ','NH  ','FL  ', &
      'MC  ','LV  ','TS  ','OG  ','X   ','BQ  ' ]

end module constants
