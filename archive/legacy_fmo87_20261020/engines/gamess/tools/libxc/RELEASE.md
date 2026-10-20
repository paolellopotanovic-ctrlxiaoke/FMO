# LibXC Interface - Release Notes

## June 30, 2022 R1 Public Release

LibXC 5.1.7 is used.


### Functionals

1.  These functionals were added as pre-defined:
    - CORRKSDT
    - CASE21
    - R++SCAN
    - R2SCAN01
    - R4SCAN
    - regTM
    - rregTM
    - MGGAC
    - RMGGAC


### Program-related changes

1.  In the interface, all common statements were replaced with modules
2.  Unifying interfaces between `libxc.src` and `libxc_empty.src`
3.  `double precision` type declarations in `libxc_empty.src` dummy file were changed into `real(kind=fp)`
4.  New arguments for `libxc_input` and `libxc_choose_functionals` routines
5.  All outputs are provided using `messages` module (see `messages.src`)



## June 30, 2021 R1 Public Release

LibXC 5.1.5 is used.


### Functionals

1.  These functionals were added as pre-defined:
    - RPBE
    - REVPBE
    - R2SCAN
    - TPSS0
    - M06-SX
    - wHPBE0
    - CAMH-B3LYP
    - TASK
    - HLTA
2.  These functionals were fixed in the LibXC interface:
    - SOGGA11
    - tau-HCTH-hyb
    - HCTH-A
    - HSEsol
3.  These functionals were fixed in the LibXC library:
    - TH1
    - KT3
    - HCTH p=7/6
    - CAP0
4.  Functional mixing can be enabled in two ways:
    - typing `FUNCTIONAL=MIXED` in $LIBXC group (preferred way);
    - skipping `FUNCTIONAL` keyword in $LIBXC group (old way)
5.  Kinetic functionals were removed: GAMESS(US) code should be more flexible before they can be enabled
6.  Double-hybrid functionals will fail when trying to calculate gradients
7.  Discrepancies between GAMESS(US) and LibXC energies of functionals (more than 1e-6 Hartree):

    |   Functional   | GAMESS(US) Energy |   LibXC Energy   | Difference |  Status            |
    |      :---:     |       :---:       |       :---:      |    :---:   |   ---              |
    |  rdft-HCTH120  |   -154.101432826  |  -154.101426930  |   5.9e-06  |  Unknown           |
    |  rdft-HCTH407  |   -154.108748553  |  -154.108651459  |   9.7e-05  |  Unknown           |
    |  rdft-MOHLYP2  |   -153.934665050  |  -153.934632077  |   3.3e-05  |  Unknown           |
    |  rdft-REVTPSS  |   -154.225528646  |  -154.085255196  |   1.4e-01  |  Broken in GAMESS  |
    |  rdft-X3LYP    |   -154.041899667  |  -154.041917751  |   1.8e-05  |  Unknown           |
    |  udft-HCTH120  |   -192.591417621  |  -192.591410532  |   7.1e-06  |  Unknown           |
    |  udft-HCTH407  |   -192.582774257  |  -192.582670977  |   1.0e-04  |  Unknown           |
    |  udft-MOHLYP2  |   -192.474560872  |  -192.474526664  |   3.4e-05  |  Unknown           |
    |  udft-REVTPSS  |   -192.746790511  |  -192.577317509  |   1.7e-01  |  Broken in GAMESS  |
    |  udft-X3LYP    |   -192.551370381  |  -192.551389248  |   1.9e-05  |  Unknown           |


### Documentation

1.  Functionals sorted by their place in Jacob's ladder and their family
2.  Fix: MGGA functionals are not included in GGA sections


### Generators

1.  Generators were merged into one file
2.  The Generators now use the Python API to LibXC (`pylibxc`) instead of parsing headers


### Compilation and installation

1.  `download-libxc.csh` and `download-libxc-tests.csh` were expanded: loading of LibXC library and LibXC tests archives does not necessarily require an active internet connection.
    If no arguments are given, the scripts will try to get the archives from `${GMS_PATH}`. If this procedure is not successful, the scripts will try to download from the Internet.
    If one argument is given, the scripts will try to use the argument as the path where the library or tests are.
2.  For Windows, `download-libxc.csh` and `download-libxc-tests.csh` were updated.


### Program-related changes

1.  `mod_nameio`'s interface was changed: use only `group_input` -- other needed routines for reading input are in this class
2.  `mod_nameio` can check allowed keys
3.  `mod_nameio` can show the location where an input error happened.
    It shows the key where the processing stopped.
4.  Fix: `mod_nameio` can correctly process groups with similar names (for example, `GGA_XC` and `MGGA_XC`)
5.  Bibliography indexes in `.log` are now sorted correctly
6.  `functional.src` provides a `functional` class so that, theoretically, several functionals can be constructed in one GAMESS(US) run
7.  Functional module only supports routines with arrays
8.  Into the `add_functional` routine of `functional` class, the optional arguments `hfex`, `alpha`, `beta`, and `omega` were added for getting coefficients of hybrid functionals from LibXC library
9.  Fixed memory problem in the interface when TD calculations are
10. Added `utils_strings` with some routines which work with strings
11. Arrays in `functionals.src` are deallocated when GAMESS(US) ends
12. Fix: GAMESS(US) terminates normally if there are errors produced by the LibXC interface
13. All `double precision` type declarations were changed into `real(kind=fp)`
14. Since FPE errors in LibXC have been fixed, `NaN`s checks were removed


### Known issues

1.  LibXC interface works incorrectly if it is built using IBM XL compilers



## June 30, 2020 R1 Public Release

LibXC 5.0.0 is used.


### Functionals

Set `DFTTYP=USELIBXC` in the `$CONTRL` group in the input file for enabling the LibXC interface.
The choosing of functionals is available in two ways:
1.  Pre-defined functionals with fixed parameters.
    The full list of available functionals can be found in the manual.
    For example, `SCAN0` (`0.25 HF + 0.75 SCAN_X + 1.00 SCAN`) will be selected:
```
 $LIBXC FUNCTIONAL=SCAN0 $END
```

2.  User-defined functionals. Do not set a `FUNCTIONAL` keyword!
```
 $LIBXC HFEX=0.25 MP2OS=0.1 MP2SS=0.05 $END
 $MGGA_X SCAN=0.65 $END
 $GGA_X PBE=0.1 $END
 $MGGA_C SCAN=0.95 $END
```
Here, the functional has the formula `0.25 HF + 0.65 SCAN_X + 0.1 PBE_X + 0.95 SCAN_C + 0.1 MP2OS + 0.05 MP2SS`,
      where `HFEX`  is a fraction of exact exchange (100% is 1.00),
            `MP2OS` is a fraction of non-local correlation energy for opposite spin,
            `MP2SS` is a fraction of non-local correlation energy for same spin.
    The `MP2` (it sets `MP2SS=MP2` and `MP2OS=MP2`) keyword also exists, but `MP2OS` and `MP2SS` have priority.
    Be careful when mixing `MP2` and `MP2SS` or `MP2OS`!
    The full list of available groups and functionals is available in the `docs-input.txt`. Happy mixing!

The total energies of both functional for He atom (3-21G/300.302) are:
1.  `-2.8763983686` (`SCAN0`)
2.  `-2.8726429406` (`MIXED`)


### Compilation and installation

For using `LibXC`, set `GMS_LIBXC` in `install.info` and `Makefile` to `true`.
The total compilation process is the following:
```
./config
./tools/libxc/download-libxc.csh
make ddi
make libxc -j$(nproc)
make modules
make -j$(nproc)
```

For Windows, see instruction here: https://github.com/gms-bbg/gamess-issues/issues/43


### Program-related changes

1.  A new input file parser (`mod_nameio.src`) that is used by the LibXC interface.
    The input file parser uses Fortran 2003 features that lead to a broken compilation when using GCC 4.9.
    The smallest checked version is GCC 5.5, and I hope that GCC 5.0 will be enough.
    The new module supports reading of the fortran data types logical, integer, double, character
    and arrays composed of these data types (single array fields are separated by comma) in free format.
    The XMARKR is not implemented.
    The maximum length of keywords and the values assigned to keywords (arrays can be entered) are 256 and 256*16 characters, respectively.
2.  The LibXC interface.
    It is composed of two files: `functionals.src` and `libxc.src`.
    The file `functionals.src` holds variables that contain information about chosen functionals.
    This file is the glue between `libxc.src` and the LibXC library.
    It calls the LibXC Fortran 2003 API.
    The file `libxc.src` contains subroutines for reading and converting arrays between GAMESS(US) and LibXC.
    Inside `libxc.src` the code-generated subroutine `libxc_choose_functionals` exists.
    Do not edit this routine by hand! The script for generating (`gen_function.py`) is in `tools/libxc/`.


### Known issues

1.  There functionals are broken (Fixed in 2021R1):
    - SOGGA11
    - tau-HCTH-hyb
    - HCTH-A
    - HSEsol
    - TH1
    - KT3
    - HCTH p=7/6
    - CAP0

2.  Double-Hybrid functional's gradients are not stopped. Results, gotten by this way, are not correct! (Fixed in 2021R1)

3.  Discrepancies between GAMESS(US) and LibXC energies of functionals (more than 1e-7 Hartree):

    |  Functional     |   LibXC Energy   | GAMESS(US) Energy |   Difference   |
    |    :---:        |       :---:      |       :---:       |      :---:     |
    |  rdft-APF       |  -153.985351865  |   -153.985352102  |  -0.000000236  |
    |  rdft-B3LYP     |  -154.116250921  |   -154.020768865  |   0.095482056  |
    |  rdft-B3LYPV5   |  -153.706080859  |   -154.020768865  |  -0.314688006  |
    |  rdft-B3PW91    |  -154.061351344  |   -154.061351912  |  -0.000000567  |
    |  rdft-B97-1     |  -153.147441042  |   -154.060554422  |  -0.913113380  |
    |  rdft-CAMQTP00  |  -154.056370566  |   -154.056370142  |   0.000000423  |
    |  rdft-HCTH120   |  -154.101426930  |   -154.101432826  |  -0.000005895  |
    |  rdft-HCTH147   |  -154.114536741  |   -154.092996196  |   0.021540544  |
    |  rdft-HCTH407   |  -154.108651459  |   -154.108748553  |  -0.000097093  |
    |  rdft-HCTH76    |  -155.314847211  |   -150.731191405  |   4.583655805  |
    |  rdft-HLYP2     |  -153.572621732  |   -153.509483856  |   0.063137875  |
    |  rdft-HLYP      |  -153.998025496  |   -154.400683586  |  -0.402658089  |
    |  rdft-M11-L     |  -153.973653255  |   -153.973653069  |   0.000000186  |
    |  rdft-MN12-SX   |  -153.703301863  |   -151.026324056  |   2.676977806  |
    |  rdft-MOHLYP2   |  -153.934632077  |   -153.934665050  |  -0.000032972  |
    |  rdft-OP        |  -154.006461178  |   -154.000471655  |   0.005989522  |
    |  rdft-PBE0M     |  -153.988704195  |   -153.988704901  |  -0.000000705  |
    |  rdft-PW91      |  -154.045812806  |   -154.045813492  |  -0.000000686  |
    |  rdft-REVTPSS   |  -154.085255196  |   -154.225528646  |  -0.140273449  |
    |  rdft-SOGGA11   |  -154.407295817  |   -154.085024608  |   0.322271208  |
    |  rdft-THCTHHYB  |  -153.620194439  |   -154.093062019  |  -0.472867579  |
    |  rdft-WB97X     |  -154.071539470  |   -154.071539351  |   0.000000118  |
    |  rdft-X3LYP     |  -154.041917751  |   -154.041899667  |   0.000018083  |
    |  udft-B3LYP     |  -192.622685369  |   -192.533460315  |   0.089225054  |
    |  udft-B3LYPV5   |  -192.227420564  |   -192.533460315  |  -0.306039751  |
    |  udft-B3PW91    |  -192.551658365  |   -192.551658596  |  -0.000000230  |
    |  udft-B97-1     |  -191.674241984  |   -192.564035505  |  -0.889793520  |
    |  udft-CAMQTP00  |  -192.560219070  |   -192.560218607  |   0.000000462  |
    |  udft-HCTH120   |  -192.591410532  |   -192.591417621  |  -0.000007089  |
    |  udft-HCTH147   |  -192.596962628  |   -192.575529088  |   0.021433540  |
    |  udft-HCTH407   |  -192.582670977  |   -192.582774257  |  -0.000103279  |
    |  udft-HCTH76    | UNCONVERGED SCF  |   -189.422648754  |     ---        |
    |  udft-HLYP2     |  -192.085115537  |   -192.062433628  |   0.022681908  |
    |  udft-HLYP      |  -192.496262805  |   -192.889962762  |  -0.393699957  |
    |  udft-MN12-SX   |  -192.198980680  |   -189.226130389  |   2.972850290  |
    |  udft-MOHLYP2   |  -192.474526664  |   -192.474560872  |  -0.000034207  |
    |  udft-OP        |  -192.480853894  |   -192.488800904  |  -0.007947009  |
    |  udft-PBE0M     |  -192.479410006  |   -192.479410292  |  -0.000000286  |
    |  udft-PW91      |  -192.555722859  |   -192.555723140  |  -0.000000280  |
    |  udft-REVTPSS   |  -192.577317509  |   -192.746790511  |  -0.169473002  |
    |  udft-SOGGA11   |  -192.893402626  |   -192.590418386  |   0.302984240  |
    |  udft-THCTHHYB  |  -192.155470213  |   -192.577752629  |  -0.422282416  |
    |  udft-WB97X     |  -192.573067117  |   -192.573067003  |   0.000000113  |
    |  udft-X3LYP     |  -192.551389248  |   -192.551370381  |   0.000018866  |

4.  mGGA functionals can not be used for DFT hessian and TD-DFT gradient calculation
5.  Some of the mGGAs require the Laplacian of the electron density, but GAMESS(US) can not provide it
6.  `mod_nameio.src` does not show the location of input errors (Partially fixed in 2021R1)
7.  It is not possible to use double-hybrid functionals with different fractions of same-spin and opposite-spin energy components for the U-SCS-MP2 scheme
8.  Some of the flags in `-ffpe-trap=invalid,zero,overflow` are incompatible with the LibXC library (Fixed in 2021R1)
9.  Some of the LibXC functionals can produce `NaN`s when energies or their derivatives are calculated.
    For energies, the `NaN` is checked, and if, in the point, the `NaN` is found this point will be skipped.
    The derivatives are not checked. (Fixed in 2021R1)
10. Sometimes GAMESS(US) produces a local kinetic energy (`tau`) that is less than zero. This leads to NaNs for some mGGA functionals.
    In order to fix that, the absolute value of `tau` is taken.
11. The LibXC interface can not be used together with FMO
12. The LibXC interface works incorrectly if it is built using PGI compilers (Fixed in 2021R1)



## Citation

 Igor S. Gerasimov, Federico Zahariev, Sarom S. Leang, Anton Tesliuk, Mark S. Gordon, Michael G. Medvedev,  
 Introducing LibXC into GAMESS (US),  
 Mendeleev Commun., 2021, 31, 302–305  
 DOI: 10.1016/j.mencom.2021.05.008  
 https://doi.org/10.1016/j.mencom.2021.04.008
