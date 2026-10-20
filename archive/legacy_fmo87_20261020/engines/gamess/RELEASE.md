# GAMESS - Release Notes

## July 15, 2024 R2 Patch 1 Public Release

GAMESS

Along with a few minor fixes:

- Add missing common block in RIMP2_VDX_I32_cublas (#928) (C. Bertoni of Argonne National Laboratory)
    - This patch allows the use of RI-MP2 energies computed using GPU offloading (e.g., EFMO)
    - Test cases impacted prior to this patch:
        - $GMSDIR/tests/efmo-gpu-offload/parallel/msn_5frag_efmo_rimp2_merged-0005N-0005.inp
        - $GMSDIR/tests/efmo-gpu-offload/parallel/msn_5frag_efmo_rimp2_merged-0010N-0005.inp
        - $GMSDIR/tests/efmo-gpu-offload/parallel/msn_5frag_efmo_rimp2_merged-0020N-0005.inp

## July 15, 2024 R2 Public Release

GAMESS

- Optimized initial guess procedure for DEA- and DIP-EOMCC calculations (#916) (J. Shen of Michigan State University)
    - Both NOACT and NUACT in the generating initial guess procedure are used for the DEA- and DIP-EOMCC calculations, instead of only using NUACT for DEA-EOMCC and NOACT for DIP-EOMCC calculations respectively in the previous version. This will reduce the computing cost of the initial guess procedure in the DEA/DIP-EOMCC calculations.
- QUAO updates (#906) (G. Schoendorff of the University of South Dakota)
    - Prints atomic charges based on QUAO populations and prints a warning when using DFT with the QUAO method
- R-8 dispersion for EFP-EFP and EFMO (#889) (P. Xu and T. Sattasathuchana of Iowa State University)
    - MAKEFP:
        - The data structure of dynamic polarizabilities are changed and made into allocatable arrays
        - Generate dipole-octopole and quadrupole-quadrupole polarizabilities for disp8
    - EFP-EFP
        - R-8 dispersion energies for EFP-EFP (both isotropic and anisotropic)
        - Analytic gradient of isotropic R-8 dispersion
        - Relevant changes for MD and MC routines
        - Modularize the main EFPDYN subroutine
    - EFMO
        - Enable R-7 and R-8 dispersion energy for EFMO
    - See $GMSDIR/tests/R-8dispersion and $GMSDIR/tests/R-7gradient for new test cases
- Implementation of Multi-FED-FCD (#888) (C. Friedl of Johannes Kepler University Linz)
    - FCD as described in J. Chem. Phys. 117, 5607–5616 (2002) (https://doi.org/10.1063/1.1502255)
    - Multi-FED-FCD as described in Photosynth Res 137, 215–226 (2018) (https://doi.org/10.1007/s11120-018-0492-1)

Along with a few minor fixes:

- Path for DEA-EOMCC(4p-2h){Nu} (CCTYP=DEAEOM3B) and DIP-EOMCC(4h-2p){No} (CCTYP=DIPEOM3A) calculation (#916) (J. Shen of Michigan State University)
- Patch set 1 (#902) (S. S. Leang of EP Analytics)
    - The exchange contribution in LCMOX had a small numerical issue with the Schwarz screening
    - Directory file had a conflict for some types of internally split parallel jobs (NSUBGR=-1)
    - A small numerical issue with FMO/PCM gradient related to RESDIM dimers
    - $PDB with multiple chains was not read correctly
    - PIEDA for DFTB/AFO had a glitch
    - Printed basis set labels in DFTB/LCAO were wrong
    - NSUBGR=-1 did not work if there were some idle groups in GDDI
    - Address issue with multi-node runs where RISM was not behaving properly
    - Address a bug with PCM
    - Address issues with primitive indexing in basg3x and basg3l
    - Address double precision numbers used as array indices in efdrvr.src
- FMO-FRET transition dipole moments are now printed correctly (#888) (C. Friedl of Johannes Kepler University Linz)

Software related changes:

- Build support updates (#924) (S. S. Leang of EP Analytics) (#900) (J. L. Galvez Vallejo of Australia National University) (G. Schoendorff of the University of South Dakota)
    - NVIDIA HPCX communication library support (see https://developer.nvidia.com/networking/hpc-x)
    - Arm Compiler for Linux (acfl) support (see https://developer.arm.com/)
    - Arm Performance Library (armpl) (see https://developer.arm.com/), NVIDIA Performance Library (nvpl) (see https://developer.nvidia.com/nvpl), AMD Optimizing CPU Libraries (aoclf) (see https://www.amd.com/en/developer/aocl.html) support
    - NEBPX is stubbed for NVHPC versions > 23 (e.g., 24.01) due to an ICE
- AAMBS changes in preparation of planned development (#920) (G. Schoendorff of the University of South Dakota)
    - Load non-relativistic AAMBS by orbital (similar to relativistic AAMBS)
    - File rename: source/modf77_vvorel.src -> source/modf77_aambs.src
- Initial stage of adding flang-new compiler support (see https://github.com/llvm/llvm-project) (#914) (J. Hammond of NVIDIA)
    - Build script for flang-new: https://github.com/jeffhammond/HPCInfo/blob/master/buildscripts/llvm-git.sh
- Updates to GAMESS continuous integration workflow (#904) (C. Friedl of Johannes Kepler University Linz)

Announcements: Minimum GCC (GFortran) support will move from 5.5.0 to a newer version in a future public release.

## March 31, 2024 R1 Public Release

GAMESS

- OpenMP target offload capabilities (#897)
    - RIMP2 energy (B. Pham, M. Alkan, D. Datta, P. Xu, and T. Sattasathuchana of Iowa State University) (C. Bertoni of Argonne National Laboratory) (L. Carrington, A. Tiwari, and S. S. Leang of EP Analytics)
        - Must specify CODE=GPURIMP2 in the $MP2 group.
        - Reference(s): J. Chem. Phys. 158, 164115 (2023)
        - RIMP2 OpenMP target offload test cases can be found in $GMSDIR/tests/rimp2-gpu-offload
    - EFMO (P. Xu, T. Sattasathuchana, and B. Pham of Iowa State University) (C. Bertoni and Y. L. Kim of Argonne National Laboratory) (S. S. Leang of EP Analytics)
        - Must specify MAKEFP=.T. in $OFFLOAD group to offload MAKEFP
        - Must specify RHF=.T. in $OFFLOAD group to offload closed-shell Hartree-Fock
        - Offload RHF and MAKEFP codes are for s and p basis functions currently
        - Reference(s): doi.org/10.1021/acs.jctc.3c01309 (memory-based, parallel CPU implementation for MAKEFP and its application for EFMO)
        - EFMP OpenMP target offload test cases can be found in $GMSDIR/tests/efmo-gpu-offload
    - ECP test cases can be found at https://github.com/gms-bbg/ecp-kpp2
    - Build instructions for the OpenMP target offload code are provided for Frontier, Perlmutter, and Polaris. These instructions can be found in $GMSDIR/machines/{olcf,nersc,alcf}
- Model core potentials for H and I functions (#885) (G. Schoendorff of Air Force Research Laboratory)

Along with a few minor fixes:

GAMESS

- Remove print statements in DDI (#894) (J. L. Galvez Vallejo of Australia National University)
- Patch to relativistic AAMBS to correctly read H and He (#884) (G. Schoendorff of Air Force Research Laboratory)
- Patch to address mismatch in the length of dictionary file DAF386 (#877) (G. Schoendorff of Air Force Research Laboratory)
- Tighter grids for non-converging TDDFT test cases (#876) (C. Friedl of Johannes Kepler University Linz)
- Address issues identified by valgrind (#875) (#878) (C. Friedl of Johannes Kepler University Linz)

Software related changes:

GAMESS

- Update LAPACK support for v3.12.0 (#893) (C. Friedl of Johannes Kepler University Linz)
- Updates to checkgms.py (#886) (G. Schoendorff of Air Force Research Laboratory)
- Documentation updates (#887) (J. Shen of Michigan State University)
- Update OSX build support (#883) (G. Schoendorff of Air Force Research Laboratory)
- Updates to runtest.py (#873) (C. Friedl of Johannes Kepler University Linz)

LIBCCHEM

- Remove of LIBCCHEM (#882) (J. L. Galvez Vallejo of Australia National University)
    - Will be replaced with CCHEM 2.0 in future public release

## September 30, 2023 R2 Public Release

GAMESS

- H and I functions for relativistic integrals (#853) (G. Schoendorff of Air Force Research Laboratory)
- RISM-SCF-cSED (#772) (I. Kosuke, and H. Sato of Kyoto University) (D. Yokogawa of the University of Tokyo)
    - Reference interaction site model (RISM) is an integral equation theory of liquids. The combination of this with quantum chemical methods enables us to evaluate solvent effects well. RISM-SCF-cSED calculations are supported for single point energy, structural optimization, frequency, IRC, and NMR calculations. The electronic structure calculation can be performed with HF, DFT, MCSCF, MP2, and CCSD. TDDFT and MCQDPT are available for excited state calculations, and conical intersections can be searched.
    - RISM-SCF-cSED is a build option when running config.
    - Input documentation available in rsmcsed/docs-input-rism.txt

LIBCCHEM

Along with a few minor fixes:

GAMESS

- Patch an issue with wrong alpha-beta spin two-electron integrals having two-hole-two-particles' indices (VBHHPP) were used when running ACP calculations for open-shell systems. (#861) (J. Shen of Michigan State University)
- Patch an issue with MRSF-TDDFT gradient code for triplets (#860) (K. Komarov of Kyungpook National University)
- Remove stray `MXGSH` and `MXG2` definitions in elgsrc.src (#842) (I. Rostov of National Computational Infrastructure at the Australian National University)
- CCSD3A and CCT3 patch to address double subtracting the number of frozen virtual orbitals from the total number of orbitals. (#835) (J. Shen of Michigan State University)

LIBCCHEM

Software related changes:

GAMESS

- Updated documentation for VVO and SYMLOC (#858) (G. Schoendorff of Air Force Research Laboratory)
- Add missing test validation files (#848) (S. Leang of EP Analytics)
- Updated documentation for mung (#846) (I. Rostov of National Computational Infrastructure at the Australian National University)
- `make cleaned_gamess` (#836) (C. Friedl of Johannes Kepler University Linz)
    - This new make target is a shortcut for `make clean_ddi clean && make`
- Windows build updates. (#834) (S. Leang of EP Analytics)
- Updates to runtest.py (#831) (#840) (C. Friedl of Johannes Kepler University Linz)
- OSX 12.4-12.6 and 13.30 build support. (#828) (#845) (G. Schoendorff of Air Force Research Laboratory)

LIBCCHEM

- Preparations for LIBCCHEM 2.0 (#851) (J. L. Galvez Vallejo of Australian National University)

## June 30, 2023 R1 Public Release

GAMESS
- Adoption of Contributor Covenant for GAMESS Code of Conduct v2.0 (#771) (L. Carrington, A. Tiwari, and S. Leang of EP Analytics) (M. Gordon, F. Zahariev, P. Xu, B. Pham, T. Sattasathuchana, D. Datta, B. Westeimer, M. Alkan, T. Harville, K. Ferreras, M. Schlinsog, G. Elliot, S. Leonard, D. Del Angel, and J. Hayes of Iowa State University) (J. L. Galvez Vallejo of Australian National University) (J. Mato of Pacific Northwest National Laboratory)
    - https://www.contributor-covenant.org/version/2/0/code_of_conduct.html
- New Michigan State University Coupled-Cluster Theory Methods: ACCSD, ACC23, ACCSD3A, ACCT3, DEAEOM2, DEAEOM2A, DEAEOM3, DEAEOM3A, DEAEOM3B, DIPEOM2, DIPEOM3, DIPEOM3A (#805) (J. Shen of Michigan State University)
    - Addition of MSUAUTO build option for GAMESS (formerly MSUCC)
    - Update to CCSDt code to enable CCSD and CCSDT calculations when corresponding numbers of active orbitals are provided.
    - In the manual docs-input.txt file, several new options are added as follows,
        - In $CONTRL group, four ACP options (CCTYP=ACCSD, ACC23, ACCSD3A, ACCT3) and seven DEA/DIP-EOMCC options (CCTYP=DEAEOM2, DEAEOM2A, DEAEOM3, DEAEOM3A, DEAEOM3B, DIPEOM2, DIPEOM3, DIPEOM3A) are added,
        - In $CCINP group, two new keywords (SHIFT and ACP) are added,
        - In $EOMINP group, two new keywords (NREF and NSAVE) are added.
    - Examples for ACP and DEA/DIP-EOMCC calculations are available in the tests/cc directory.
- LUT-IOTC Updates (#798) (C. Takashima of Waseda University)
    - LUT-IOTC method for two-electron term (RELWFN=LUTIOTC2 in $CONTRL group) becomes available with divide-and-conquer method (DCFLG=.T. in $DANDC group).
    - LUT-IOTC method for two-electron term (RELWFN=LUTIOTC2 in $CONTRL group) becomes available with direct SCF method (DIRSCF in $SCF group). In this implementation, one-center IOTC two-electron integrals are stored on disk storage, and multi-center non-relativistic two-electron integrals are computed in each SCF iteration.
    - Since RELWFN=LUTIOTC2 is not parallelized at present, direct SCF has a longer computational time than the conventional SCF calculation.
    - At present, (LUT-)IOTC methods for two-electron term (RELWFN=IOTC2E and LUTIOTC2 in $CONTRL group) are not parallelized. So, it is set up to stop the job with an error message when parallel calculation will be run.
    - In MP2 calculation, (LUT-)IOTC methods for two-electron term (RELWFN=IOTC2E and LUTIOTC2 in $CONTRL group) are only available with CODE=SERIAL in $MP2 group. The code is modified to choose CODE=SERIAL in the default calculation of RELWFN=IOTC2E and LUTIOTC2. And it is set up to stop the job with an error message when other keywords of CODE are chosen.
- RI-G3 (PDG) Method (#756) (B. Pham and D. Datta of Iowa State University)
    - A modification of the G3 method that utilizes RI-MP2 and RI-CC
    - See test case: $GMSDIR/tests/comp-mpiomp/rig3_exam43.inp
    - Phys. Chem. A 2021, 125, 42, 9421-9429
- VB2000 v3.0.1 update (#742) (D. de Sousa of the Federal University of Rio de Janeiro)
- Accurate Atomic Minimal Basis Sets (AAMBS) (#726) (G. Schoendorff of Air Force Research Laboratory)
    - Implementation of the relativistic AAMBS for Li-Rn.
    - This implementation allows for the formation of VVOs and SVD localization for methods using all-electron basis sets with relativistic methods and for the formation of VVOs when using model core potentials.
    - New options in $LOCAL are the following:
        - VVTYP selects the AAMBS version where IVVTYP=0 selects the original non-relativistic AAMBS and IVVTYP=1 selects the relativistic AAMBS implemented here.
        - PADSVD reorders the outer core if a valence orbital in a molecule is lower in energy than a core orbital, e.g., the O 2s in CeO and CeO2 is lower in energy than the Ce 5s an 5p.
        - GUESS=AAMBS - Huckel guess that uses the AAMBS rather than the MIDI sets) expanded valence - nd orbitals for the p-block that allow for sp3d and sp3d2 hybridization.
        - SUBVAL option to include the 5s and 5p orbitals in the valence space for the 4f elements.

LIBCCHEM
- Updates (#760) (J. L. Galvez Vallejo and G. Barca of Australia National University)
    - Changed the NVIDIA intrinsic boys function back to the interpolation to ensure portability among new GPU.
    - Brought back the routines to copy the arrays on to the GPU.
    - Resolve memory issues.

Along with a few minor fixes:

GAMESS
- LUT-IOTC Patch (#798) (C. Takashima of Waseda University)
    - Address out of bound access
- FMO Patch (#766) (#813)
    - Numerical gradients did not work for FMO with GDDI.
    - Hessian reading did not work for RUNTYP=FMOHESS.
    - Here was a technical problem in reading MOs for SCZV in FMO.
    - Manual and aux file update.
- Fix active orbital space in CCSD3A and CCT3 calculations. (#797) (Piotr group of Michigan State University)
    - A bug with how very small active spaces (e.g. 0-0) were treated when open shell CCSD3A and CCT3 calculations were performed, was found and corrected in this PR. Essentially, some diagrams which should not have been computed, had been included in the calculation due to faulty conditional statements related to the size of the active orbital space.
    - This issue only impacted calculations with very small active spaces and very few electrons, which would result in very small active occupied and active unoccupied orbital spaces.
- Synchronize common blocks in demrpt.src and fmolib.src (#790) (C. Friedl of Johannes Kepler University Linz)
- Address issues identified by valgrind (#768) (#769) (#770) (#781) (C. Friedl of Johannes Kepler University Linz)

Software related changes:

GAMESS
- LAPACK v3.11.0 update (#774) (S. Leang of EP Analytics) (C. Friedl of Johannes Kepler University Linz)
- Build updates for 2023 R1 public release
    - Fix msuauto compilation (#819) (S. Leang of EP Analytics)
    - Address ibm64 build issues (#817) (S. Leang of EP Analytics)
    - OLCF Frontier build support (#799) (S. Leang of EP Analytics)
    - DoD HPCMP Navy DSRC Nautilus build support (#794) (S. Leang of EP Analytics)
    - OSX Ventura 13.2.1 build support (#787) (G. Schoendorff of Air Force Research Laboratory)
    - fj-a64fx updates (#777) (I. Rostov of National Computational Infrastructure at the Australian National University)
    - Enable parallel build of GAMESS (#745) (#727) (C. Friedl of Johannes Kepler University Linz)
        - Prior to this PR, DDI and GAMESS modules had to invoke a separate `make`. After this PR, GAMESS can be compiled with just `make -j N` where N is the number of parallel make commands to invoke.
    - ALCF Polaris build support (#757) (S. Leang of EP Analytics)
    - HPE Cray compiler support (#775) (#754) (L. Roskop of HPE Cray)
- DDI
    - compddi updates impacting DDI large message factor builds (#795) (S. Leang of EP Analytics)
    - Improve output (#762) (C. Friedl of Johannes Kepler University Linz)
- Code modernization
    - Remove ENTRY statements:
        - dftxca.src (#764) (C. Friedl of Johannes Kepler University Linz)
        - ecp.src (#763) (C. Friedl of Johannes Kepler University Linz)
    - Reduce legacy preprocessing in unport.src (#739) (C. Friedl of Johannes Kepler University Linz)
    - Pull MXFRZ parameter into mx_limits module (#727) (G. Schoendorff of Air Force Research Laboratory)
- Continuous integration and test suite
    - Legacy test suite updates (#778) (I. Rostov of National Computational Infrastructure at the Australian National University)
    - Test suite updates (#765) (#779) (C. Friedl of Johannes Kepler University Linz)
    - Update test case:
        - utddft-o2-M06-2X-energy.inp (#784) (C. Friedl of Johannes Kepler University Linz)
        - mc-ormpt2-lif-protect.inp (#783) (C. Friedl of Johannes Kepler University Linz)
    - create-install-info.py updates (#776) (#818) (C. Friedl of Johannes Kepler University Linz)
    - TOTAL WALLTIME parsing (#753) (C. Friedl of Johannes Kepler University Linz)
    - NPY-for-Fortran update (#749) (C. Friedl of Johannes Kepler University Linz)
    - Enable debugger attachment during GAMESS execution (#747) (C. Friedl of Johannes Kepler University Linz)
    - Jenkinsfile update and introduction of bwrap (#746) (#773) (#791) (#793) (#800) (C. Friedl of Johannes Kepler University Linz) (#796) (S. Leang of EP Analytics)

## September 30, 2022 R2 Public Release

GAMESS
- Implementation of Fragment Excitation Difference (FED) for the calculation of EET couplings (#688) (C. Friedl of Johannes Kepler University Linz)
    - Capabilities documented under $TDDFT (see docs-input.txt)
    - Test cases under $GMSDIR/tests/tddft/parallel:
        - dn2\_full\_fed.inp
        - dn2\_no\_no\_fed.inp
        - ethylene\_benzene\_fed.inp
        - stacked\_naphtalene\_3\_5\_a\_fed.inp
        - stacked\_naphtalene\_4\_0\_a\_fed.inp
    - See J. Phys. Chem. C 2008, 112, 1204-1212 https://doi.org/10.1021/jp076512i
- ORMAS-PT2 INTRNL keyword (#725) (J. Ivanic)
    - Adds an option to include internal (active -> active) PT2 energy contributions in ORMAS-PT2 energies
    - See $DETPT group in docs-input.txt for more enlightening comments from the ORMAS developers
    - Default behavior:
      - Single state runs (NPTST=1): `INTRNL=.TRUE.`
      - Multistate runs (NPTST>1): `INTRNL=.FALSE.`
    - See test case: $GMSDIR/tests/mcscf/mrpt/parallel/h2co-or1210-pt2_n-pi-star-excite.inp

Along with a few minor fixes:

GAMESS
- FMO Patches (#740) (D. G. Fedorov of the National Institute of Advanced Industrial Science and Technology Japan) (C. Friedl of Johannes Kepler University Linz)
   - FMO-CC did not set the number of virtual orbitals correctly
   - CAM-B3LYP gradient was incorrect for FMO
   - RI methods did not work properly for FMO
   - Address MD simulations except NAMD which had a memory-related bug (a memory leak)
   - Adjust validation files to reflect updated output generated by MD simulations
   - New MD test case (exam49.inp) added to standard (and travis-ci)
- Address varying common block definitions (#735) (C. Friedl of Johannes Kepler University Linz)
   - neodft.src, neohss.src, dftbx.src, dftgrd.src, fcidump.src, fmo.src, fmoprp.src, int2a.src, prpel.src, prplib.src
- Patches for qrel.src (#698) (G. Schoendorff of Iowa State University)
   - Initialization of `DUM1` and `DUM2` in subroutine `PVPINT`
   - Define `C` as `CLIG` in `QRDYY`
- Correct `RS6` value for PBE0 in dftdis.src (#731) (J. Zou Nanjing University)
- Local Response Dispersion (LRD) patch (#696) (Y. Ikabata of Toyohashi University of Technology)
   - Atomic pseudo-polarizabilities of symmetry unique atoms were copied into symmetrically equivalent atoms. Such a treatment properly works in the case of the lowest order (Eq. (27) with l=l'=1 in J. Chem. Phys. 131, 224104 (2009) https://doi.org/10.1063/1.3269802) because polarizability density is isotropic, but the higher order is not so.
The problem has not been solved, but tentatively symmetry is disabled.
   - Three-center and four-center dispersion terms are corrected as reported in the original paper: J. Chem. Phys. 133, 194101 (2010) https://doi.org/10.1063/1.3503040
   - Seminumerical Hessian is set as default

Software related changes:

GAMESS
- Build updates for 2022 R2 public release (#748) (S. Leang of EP Analytics)
    - LAPACK build support under Intel OneAPI on Windows
    - Isolate and defensive compilation for routine identified by EP Analytics COSE tool in mctwo.src when compiling with AOCC
    - Cross-compatibility changes of test case names containing '+'
    - Update special handling of test cases under Windows
    - Bump GCC compiler versions for 10 and 11
    - Add missing basis used by LIBCCHEM test case into auxdata/AUXBAS/auxdef21
- GCC 12 build support (#730) (S. Leang of EP Analytics) (C. Friedl of Johannes Kepler University Linz)
- Test suite update (#719) (I. Rostov of National Computational Infrastructure at the Australian National University)
- Replace sleep with POSIX fsync in subroutine ENDING and ABRT (#713) (C. Friedl of Johannes Kepler University Linz)
- Remove legacy preprocessing of vb2gms.src (#708) (C. Friedl of Johannes Kepler University Linz)  
- Update AMD Optimizing C/C++ and Fortran Compilers (AOCC) build support (#706) (S. Leang of EP Analytics) (I. S. Gerasimov of Kyungpook National University)
- Remove Tinker atom count hack (#703) (C. Friedl of Johannes Kepler University Linz)
   - See `tinker_maxatom` variable in `module mx_limits`
- Merge `*I64` machine type activation lines `*W64` (#694) (C. Friedl of Johannes Kepler University Linz)
- Update Windows build support (#693) (S. Leang of EP Analytics) (I. S. Gerasimov of Kyungpook National University)
- Remove File-Get-Environment from GAMESS (#692) (I. S. Gerasimov of Kyungpook National University)
- Merge `*INT` machine type activation lines `*L64` (#687) (C. Friedl of Johannes Kepler University Linz)
- Move portable code from unport.src to util.src (#686) (C. Friedl of Johannes Kepler University Linz)
   - Build patches (#736) (C. Friedl of Johannes Kepler University Linz)


## July 31, 2022 R1 Public Release

### Notice

GAMESS now requires a LAPACK library for compilation. Your math library may already provide this library. However, if you receive an error during linking about undefined references to dsygst and dspev then you will need to link to a LAPACK library.

You can try compiling the Netlib LAPACK library from scratch by one of the following two options below:

 1. (Requires CMake)

```
./tools/lapack/download-lapack.csh
make -j lapack
```

 2. (Does not require CMake)

```
make -j lapack USE_CMAKE=false
```

Or provide a valid GMS\_LAPACK\_LINK\_LINE in install.info to an installation of LAPACK (compiled with 8-byte integer) already on your 64-bit system.

For example:

```
setenv GMS_LAPACK_LINK_LINE "-L/path/to/lapack/library -llapack"
```

You do not need to specify GMS\_LAPACK\_LINK\_LINE if you decide to do option 1. or 2.

----

GAMESS
- Crusher build support (Frontier TDS at OLCF) (#716) (S. Leang of EP Analytics)
    - Load the appropriate environment (e.g., PrgEnv-cray, PrgEnv-amd, PrgEnv-gnu) and then run ./config
    - Specify `hpe-cray-ex` target and `crusher` HPC system target
- Enable the substitution of virtual canonical orbitals for VVOs in energy calculations for SF-ORMAS and SF-ORMAS-PDFT (#676) (K. Ferreras of Iowa State University)
    - The SF-ORMAS-PDFT approach combines the SF-ORMAS-CI (Phys.Chem.Chem.Phys., 2018, 20, 2615) method with MC-PDFT (J. Chem. Theory Comput., 2014, 10, 3669) method to treat static and dynamical correlation in multiconfigurational systems. The static correlation description is generated via the spin-flip approach (Chem. Phys. Lett., 2001, 350, 522) a well-known and inexpensive method that treats excited states with multiconfigurational character using a high-spin single reference determinant. The dynamical correlation is treated with the on-top pair-density functional theory which translates GGA density functionals, successfully escaping the double counting of electron correlation and symmetry dilemma. The SF-ORMAS-PDFT method requires input for both calculations, SF-ORMAS-CI and PDFT, and only runs in serial. It can also be combined with virtual valence orbitals (J. Chem. Phys., 2004, 120, 2629) to decrease the number of virtual orbitals used in the calculation with keyword SFOVVO=.TRUE. in the $SCF group.
    - Example inputs are available in $GMSDIR/tests/sformpdft
- QM/EFP2 parallelization improvements of two-electron integrals for calculations with large number of EFP fragments (#666) (P. Xu and T. Sattasathuchana of Iowa State University)
- Enabled Hessian calculations when using ECPs with an h-ul potential and coefficient of zero by increasing the value of NCOEF2. (#664) (G. Schoendorff of Iowa State University)
    - Example inputs are available in $GMSDIR/tests/ecp
- Implementation of LUT-IOTC method for two-electron term (#647) (C. Takashima of Waseda University)
    - Example inputs are available in $GMSDIR/tests/lut-iotc-2e
- Implementation of picture change corrected DFT (#647) (C. Takashima of Waseda University)
    - See `PCCDFT` option in $RELWFN group
- Implementation of REKS and MRSF with NAMD support (#635) (K. Komarov, S. Lee, M. Filatov, and C. H. Choi of Kyungpook National University)
    - MRSF-TDDFT (Mixed-Reference Spin-Flip)-TDDFT is a new linear response time-dependent density functional theory for both ground and excited electronic states.(J. Chem. Phys. 2018, 149, 104101; J. Chem. Phys. 2019, 150, 184111) Currently, the energy and gradient of MRSF-TDDFT are available for ROHF reference. From its high spin triplet reference, MRSF-TDDFT generates singlets including ground state, triplets as well as quintets. MRSF-TDDFT eliminates the bulk of notorious spin contamination of SF-TDDFT. As compared to LR-TDDFT, MRSF-TDDFT includes the important doubly excited configurations providing the correct conical intersection topology between electronic states.(J. Phys. Chem. A 2019, 123, 6455; J. Phys. Chem. Lett. 2021, 12, 9720) The MRSF-TDDFT also produces ground singlet state allowing to study diradical systems.(J. Chem. Theory Comput. 2021, 17, 848) Both ionization potentials and electron affinities of both ground and excited states can be calculated by quasiparticle state calculations.(J. Phys. Chem. Lett. 2021, 12, 9963) A new searching method for three or higher state conical intersection is implemented for MRSF-TDDFT.(J. Phys. Chem. A 2021, 125, 1994) Fast computational method is developed for NACT (Nonadiabatic Coupling Term) and NACME (Nonadiabatic Matrix Elements).(J. Phys. Chem. Lett. 2021, 12, 4722; J. Chem. Theory Comput., 2019, 15, 882) NAMD(Nonadiabatic Molecular Dynamics) is also implemented for MRSF-TDDFT.(J. Phys. Chem. Lett. 2021, 12, 4339) Now, the LC (Long Range-Corrected) XC functionals can be used for both SF- and MRSF-TDDFT.
        - Example inputs are available in $GMSDIR/tests/mrsf
    - Spin-Restricted Ensemble-Referenced Kohn-Sham (REKS) theory is a variational time-independent method for both ground and excited electronic states utilizing an ensemble density, which overcomes the limitations of traditional KS-DFT. (Wiley Interdisciplinary Reviews: Computational Molecular Science 2015, 5, 146). The energy and gradient of the particular State-average REKS(2,2) and its variants are implemented. The theory can produce correct conical intersection topology between ground and excited states. Open shell singlet and doubly excited configurations are included. Ionization and electron affinity by EKT for both ground and excited electronic states can be obtained.
        - Example inputs are available in $GMSDIR/tests/reks
- Removed common blocks from LibXC interface (#559) (I. S. Gerasimov of Kyungpook National University)

Along with a few minor fixes:

GAMESS
- FMO patches (#705) (D. G. Fedorov of the National Institute of Advanced Industrial Science and Technology Japan)
    - Some PIEDA term labels were incorrect for DFT
    - AFO erroneously reset NPRINT to -5.
    - FMO-MCSCF numerical gradients crashed (a file conflict).
    - PA had major problems with processing some types of PDB data.
    - MCSCF did not process ESP derivatives properly.
    - MD was erroneously set to ignore NGUESS=8.
    - PCM[1(2)] did not work with DFTB.
    - PCM[1(2)] nor with the partial screening.
- QM/EFP2 patches (#666) (P. Xu and T. Sattasathuchana of Iowa State University)
- OpenMP patch for IBM64 targets (#667) (I. S. Gerasimov of Kyungpook National University)
- Print control for semi-empirical runs (#710) (C. Friedl of Johannes Kepler University Linz)
- Variable initializations in Coupled Cluster (#672) (C. Friedl of Johannes Kepler University Linz)

Software related changes:

GAMESS
- Run support updates for Fujitsu A64FX (I. Rostov of National Computational Infrastructure at the Australian National University)
    - see $GMSDIR/machines/fj-a64fx/rungms
- Improvements to $GMSDIR/tests/checkgms.py (#707) (C. Friedl of Johannes Kepler University Linz)
- Refactor of $GMSDIR/bin/create-install-info.py (#684) (E. Woodruff and S. Leang of EP Analytics)
- GAMESS internal BLAS and LAPACK routines have been replaced with Netlib (#616) (I. S. Gerasimov of Kyungpook National University) (#701) (S. Leang of EP Analytics) (C. Friedl of Johannes Kepler University Linz) (#717) (I. Rostov of National Computational Infrastructure at the Australian National University)
- NEO is now part of GAMESS continuous integration testing (#642) (I. S. Gerasimov of Kyungpook National University)
- Build improvements
    - (#699) (S. Leang of EP Analytics)
    - (#348) (I. S. Gerasimov of Kyungpook National University)
- Drop build support of GCC compilers < 5.0 (#605) (I. S. Gerasimov of Kyungpook National University)
- Per-testcase coverage testing with gcov (#547) (C. Friedl of Johannes Kepler University Linz)
- Remove 32-bit Linux support from GAMESS (#568) (I. S. Gerasimov of Kyungpook National University)

LIBCCHEM
- Increase maxreg count (#678) (J. L. Galvez Vallejo of Iowa State University)

## September 30, 2021 R2 Patch 2 Public Release

Please note that the latest GAMESS build instructions are available in $GMSDIR/README.md (same directory level as this document).

GAMESS
- Build support updates (#674) (#629) (#626) (S. Leang of EP Analytics) (C. Bertoni of Argonne National Laboratory) (M. M. Lazaretti Tormena)
    - NERSC Perlmutter (Cray EX)
    - NERSC Cori (Cray XC 40)
    - Intel OneAPI 2022.1.2
    - NVIDIA HPC SDK 22.1
    - mac64 target updates
    - additional build options added to create-install-info.py
- NEBPATH documentation update (#644) (S. Koseki of Osaka Prefecture University, Japan)

LIBCCHEM
- New accelerated Hartree-Fock (AHFOCK) code for NVIDIA GPUs (#590) (J. L. Galvez Vallejo of Iowa State University) (D. Poole of Iowa State University) (M. Alkan of Iowa State University) (G. M. J. Barca of the Australian National University)
    - new GPU accelerated HF to libcchem with capacity up to F functions, default up until D
    - new GPU code uses Head-Gordon-Pople as its ERI engine
    - new code works with multiple version of g++, tested at least 4.8, 6.4, 7.4, 8.5 and 9.3
    - boys function is evaluated using an nvidia intrinsic, the boys interpolation tables are left there in case they are needed later they are NOT used
    - correction for higher angular momentum is done on the GPUs
    - see $GMSDIR/tests/libcchem/ahf and $GMSDIR/tests/libcchem/genfock for test cases
    - https://pubs.acs.org/doi/abs/10.1021/acs.jctc.0c00768
    - https://ieeexplore.ieee.org/abstract/document/9355281
    - https://pubs.acs.org/doi/abs/10.1021/acs.jctc.1c00720
    - https://dl.acm.org/doi/abs/10.1145/3458817.3476222

Along with a few minor fixes:

GAMESS
- FMO patches (#657) (#637) (D. G. Fedorov of the National Institute of Advanced Industrial Science and Technology Japan)
    - In DFTB Raman, symmetry labels were undefined resulting in NULL symbols printed
    - IR intensity was zero for FMO-MP2
    - GLOBOP with FMO used unallocated memory in some cases
    - Improve error diagnostics
    - Documentation update
    - Spherical constraints in $STATPT did not work when more than 9
were specified
    - Gradient projection issue was patched for DFTB/PBC
    - Patch issue with atoms outside of DFTB/PBC cell
    - Numerical gradient did not work for FMO with GDDI
    - Empty trajectory files were created in FMO/MD
    - Umbrella sampling constraints were not working for some instances in FMO/MD
    - Load balancing was not proper for some special runs of FMO/PCM
- NBO configuration patch in config and lked (#675) (I. Rostov of National Computational Infrastructure at the Australian National University)
- rimp2grd-mpiomp test patches (#659) (C. Friedl of Johannes Kepler University Linz) (E. Guidez of University of Colorado Denver)
- MCP ZFK* patches (#655) (S. Leang of EP Analytics) (T. Zeng of York University) (D. Small of the Molecular Graphics and Computation Facility at UC Berkeley)
    - Corrections for ZFKN  MCP files for missing parameters
    - MCP test cases
- VB2000 patches
    - fj-a64fx build support(#653) (I. Rostov of National Computational Infrastructure at the Australian National University)
    - compall update (#632) (S. Leang of EP Analytics)
- Address FPEs and variable initializations
    - FFMD2 (#650) (C. Friedl of Johannes Kepler University Linz)
    - rxncrd.src (#646) (C. Friedl of Johannes Kepler University Linz)
    - uninitialized variables identified by valgrind (#649) (C. Friedl of Johannes Kepler University Linz)
- Stop GAMESS from proceeding during NO GRADIENT conditional check (#639) (I. S. Gerasimov of Kyungpook National University)
- LOCAL=SVD with SCFTYP=RHF/ROHF when EXTLOC=NONE/SPLITQ3 patch (#531) (G. Schoendorff of Iowa State University)
- Fix OpenBLAS support with NVIDIA HPC SDK (#656) (S. Leang of EP Analytics)
- Fix C print format specifier for INT_MAX (see DDI Large) (#611) (I. Rostov of National Computational Infrastructure at the Australian National University)

Software related changes:

GAMESS
- Bubblewrap sandboxing in runtest.py (#638) (C. Friedl of Johannes Kepler University Linz)
    - see --bwrap argument in runtest.py
- Replace statement function with F90 style function in mccas.src (#534) (I. S. Gerasimov of A.N.Nesmeyanov Institute of Organoelement Compounds of Russian Academy of Sciences)
- Change quadruple precision handling (#426) (I. S. Gerasimov of A.N.Nesmeyanov Institute of Organoelement Compounds of Russian Academy of Sciences)
- Update documentation (#626) (S. Leang of EP Analytics)

## September 30, 2021 R2 Patch 1 Public Release

Software related changes:

GAMESS
- Address build issues. No source files were changed.
    - config updates (#625) (S. Leang of EP Analytics)
        - TINKER, VB2000 XMVB, NEO, NBO build options added to **config** and **create-install-info.py**
        - Not specifying a `yes`/`no` response during the build option stage of **config** will result in the default response of `no`
    - Address linking issues under GCC 7.x and 10.x (#624) (S. Leang of EP Analytics)
    - Serial build option update (#614) (I. S. Gerasimov of Kyungpook National University)
    - MKL handling under Fedora (#609) (#617) (S. Leang of EP Analytics) (I. S. Gerasimov of Kyungpook National University)

## September 30, 2021 R2 Public Release

GAMESS
- VB2000 updates (#591) (#601) (D. de Sousa of the Federal University of Rio de Janeiro)
    - See $GMSDIR/vb2000/DOC/README for more information
- FMO updates and enhancements (#587) (#598) (S. Leang of EP Analytics)
    - Improve MCSCF convergence for FMO runs
    - Add LOGN (logical node count) argument to rungms*
    - new sample inputs
        - tools/fmo/samples/2021/DFTB.dens.inp
        - tools/fmo/samples/2021/LCMOX.inp
        - tools/fmo/samples/2021/PAVE.inp
        - tools/fmo/samples/2021/dual\_basis.inp
        - tools/fmo/samples/2021/DFTB.dens.inp
- The Molecular Sciences Software Institute (MolSSI) Driver Interface (MDI) support (#583) (T. Barnes of the Molecular Sciences Software Institute)
    - See https://molssi-mdi.github.io/MDI_Library/html/index.html for more information
    - https://doi.org/10.1016/j.cpc.2020.107688
- Build support updates
    - Add support for CCE/12.0.1 and CCE/12.0.3 under `hpe-cray-ex/cs`, NVIDIA HPC SDK under `linux64`, OneAPI Beta under `linux64` (#599) (S. Leang, A. Tiwari, and E. Woodruff of EP Analytics) (C. Bertoni of Argonne National Laboratory) (L. Roskop of HPE Cray)
    - `fj-a64fx` support update (I. Rostov of National Computational Infrastructure at the Australian National University)
    - OpenBLAS support update for 64-bit integers (#582) (C. Friedl of Johannes Kepler University Linz)
    - Clean up build output (#581) (C. Friedl of Johannes Kepler University Linz)
    - GNU compilers (#567) (S. Leang of EP Analytics) (#596) (I. S. Gerasimov of Kyungpook National University)
- VeraChem VM2 interface. VM2 is a commercial statistical mechanics based MM package (#546) (P. Xu and T. Sattasathuchana of Iowa State University) (S. Webb of VeraChem)
    - Visit https://www.verachem.com for more information
- NEBPATH updates (#541) (S. Koseki of Osaka Prefecture University, Japan)
- High angular momentum basis functions (H,I) support (#530)
    - in RI-MP2 support (E. Guidez of University of Colorado Denver)
        - For examples see:
            - tests/rimp2-mpiomp/Br2.inp
            - tests/rimp2-mpiomp/I2.inp
            - tests/rimp2grd-mpiomp/Br2.inp
            - tests/rimp2grd-mpiomp/I2.inp
            - tests/rimp2grd-mpiomp/I2NCH.inp
            - tests/rimp2grd-mpiomp/w2o2.hess.ecp.inp
    - in RI-CC support (D. Datta of Iowa State University)

Along with a few minor fixes:

GAMESS
- FMO related patches (#587) (#598) (S. Leang of EP Analytics)
    - Address variable declaration in **fmoio.src** and **fmoprp.src**
    - Address file handling in **dftbpb.src**
    - Address program logic in **dftbpb.src**
    - Address SMD parallelization in **mnsol.src**
    - Address FPE in **fmoio.src**
    - Address issue with reading in MP2 group in **fmolib.src**
    - Update tools/fmo/samples/2021/PA.inp
    - Address issue in **fmoio.src** causing **gfortran** to crash
    - Fix variable initialization in **dftbpb.src**
    - Address convergence issues with DFTB calculations
    - Address memory allocation in SMD
- Guard allocates and deallocates in **rimp2grd.src** which posed issues when performing semi-analytic Hessian calculations with OpenMP-threaded RI-MP2 code using ECPs (#586) (S. Leang of EP Analytics)
- Fix variable initialization in **efdrvr.src** and **efmogrd.src** (#584) (C. Friedl of Johannes Kepler University Linz) (P. Xu and T. Sattasathuchana of Iowa State University)
- Address "Unknown element" issue on `ibm64` (#570) (S. Leang of EP Analytics)
- Address 32-bit issues (#550) (C. Friedl of Johannes Kepler University Linz)

Software related changes:

GAMESS
- Jenkins-CI small and large testing use GNU compiler version 5.5.0 (#602) (I. S. Gerasimov of Kyungpook National University)
- Test suite updates (#575) (S. Leang of EP Analytics)
    - Add `--filepath` flag to **runtest.py** and **checkgms.py** to allow for absolute or relative filepath strings that can be used for filtering
- Common block AND module support:
    - IOFILE, PAR, DFTPAR (#532) (I. S. Gerasimov of A.N.Nesmeyanov Institute of Organoelement Compounds of Russian Academy of Sciences)
        - see $GMSDIR/source/modules\_common.src for example usage
- Replace LOGAND, LOGOR, LOGXOR, LOGSHF in unport.src with standard routines IAND, IOR, IEOR, ISHFT (#427) (I. S. Gerasimov of A.N.Nesmeyanov Institute of Organoelement Compounds of Russian Academy of Sciences)
- Legacy code removal:
    - DDI v1 and DDI v2 (#399) (I. S. Gerasimov of A.N.Nesmeyanov Institute of Organoelement Compounds of Russian Academy of Sciences)
- Common block to module conversions:
    - MCPDAT, MCPITM, MCPITMG, MCPOUT, MCPOUTG, MCPWRKC, MCPWRKCG, MCPWRKD, MCPWRKDG (#306) (I. S. Gerasimov of A.N.Nesmeyanov Institute of Organoelement Compounds of Russian Academy of Sciences)

Notice

- GAMESS minimum compiler version support is GCC-5.
- 32-bit build support will be removed from GAMESS after the September 30, 2021 R2 public release.

## June 30, 2021 R1 Public Release

GAMESS
- Fujitsu ARM64 platform (fj-a64fx) support (#481)(#506) (I. Rostov of National Computational Infrastructure at the Australian National University)
    - See $GMSDIR/machines/fj-a64fx/readme.txt for more information.
- Build support for HPE Cray EX series (#558) (S. Leang of EP Analytics).
    - `GMS_HPC_SYSTEM_TARGET` support for DoD Navy DSRC `narwhal` and DoE OLCF `spock`.
    - See $GMSDIR/machines/hpcmp/narwhal for PBS submission script.
- Build support for Apple M1 with GNU GCC 11 (#558) (G. Schoendorff of Iowa State University).
- LibXC interface update to support LibXC 5.1.5 (#476) (I. S. Gerasimov of A.N.Nesmeyanov Institute of Organoelement Compounds of Russian Academy of Sciences). More information available in $GMSDIR/tools/libxc/RELEASE.md
    - LibXC interface can be used as well as the internal DFT library except FMO runs
    - New functionals were provided; some old ones were fixed
    - Kinetic functionals were removed from the GAMESS LibXC interface
    - It is not possible to calculate gradients with double hybrid functionals
- FMO Updates (#513) (D. G. Fedorov of the National Institute of Advanced Industrial Science and Technology Japan)
    - FMO-DFTB/PBC
    - Intranode GDDI and core affinity for sockets
    - Partition analysis for DFTB
    - FMO/EDA3
- Enable larger MPI message sizes in DDI (going beyond 2GB limit) (#512) (L. Carrington of EP Analytics)
    - This use case becomes necessary for the OpenMP threaded version of GAMESS where the number of GAMESS compute process (and equivalent number of DDI processes) do not saturate the core count per node.
    - To use DDI large messages (e.g., >2GB) you need to set the flag to be `DDI_MSG_FACTOR` > 1.
    - Common usage:
        - `set EXTRA_FLAGS=' -DDDI_MSG_FACTOR=8'`
        - additionally a threshold value is set in $GMSDIR/ddi/src/common.h : `DDI_LG_MSG_THRESHOLD` so that only messages greater than the threshold are set up to transfer with the `DDI_MSG_FACTOR`.
The value of the flag should be an integer value and it is used to divide the current message by that value to ensure it fits into the "count" limit of MPI.

- All for the use of ECPs with h-ul potentials that have coefficients of zero (#511) (G. Schoendorff of Iowa State University). The h-ul potential provides a reference point for the computation of difference potentials, hence the coefficient of zero. This modification allows for the use of many of the ECPs available for heavy elements, particularly the cc-pVnZ-PP sets and the Stuttgart ECPs. This modification also ensures that GAMESS can reproduce the results obtained using these ECPs with NWChem and Molpro. 
- Change keyword to use threaded RIMP2 gradient code (code=omprimp2). See $MP2 in docs-input.txt (#445) (B. Q. Pham of Iowa State University)
- QuanPol update (#417) (H. Li) (University of Nebraska-Lincoln)
    - New Features
        - DFTB, TDDFTB, AM1 and PM3 can be used as QM methods in QuanPol QM/MM.
        - MD simulation ($QUANPO DIFFUSE=1 COM) of ion mobility and collision cross section (CCS) in neutral buffer gases.
        - Force field global optimization method ($QUANPO keyword MDOPT, coded by Rui Lai).
        - Pairwise approximation method ($QUANPO IDIMER=3) for larger QM region in QM/MM.
    - Improvements
        - Fixed a few FPE errors.
        - Fixed a few bugs in free energy simulation.
        - Automatic selection of QM atoms for QM/MM.
        - Added/modified a few tests in $GMSDIR/tests/quanpol
        - Added new force field files to $GMSDIR/auxdata/QUANPOL
- Enable 1e- and 2e- integrals generated by ROHF code in to be written to a file (FCIDUMP) for external use (#407) (J. E. Deustua)

Along with a few minor fixes:

GAMESS
- Rover Holleriths character data from write statements (#558) (N. Prewitt of SAIC).
- Improve multi-node support on fj-a64fx machines for OpenMP threaded RIMP2 gradients (#529) (I. Rostov of National Computational Infrastructure at the Australian National University)
- Synchronize with subroutine and caller parameter list for `DFTB_CALC_DWDER` in dftbhs.src (#545) (C. Friedl of Johannes Kepler University Linz)
- FMO Updates (#513) (D. G. Fedorov of the National Institute of Advanced Industrial Science and Technology Japan)
    - Corrected PIEDA summary output for RC and BS terms
    - CC/PCM calculations were incorrect
    - PCM Hessians had a parallel bug
    - OPTFMO had several bugs
    - FMO-DFTB Hessian had a bug
    - FMO/FRET had bugs when excited fragments were not initial
    - FMO gradient had a bug when some linear dependent MOs are removed
    - FMO seminumeric Hessians did not work
    - FMO numerical gradients did not work
    - Several major bugs for FMO/ECP
    - HF-3c: added a few more supported elements (not a bug)
- Fix MN12-SX and HCTH/147 functionals (#451) (I. S. Gerasimov of A.N.Nesmeyanov Institute of Organoelement Compounds of Russian Academy of Sciences)
- Abort MRPT calculation if E(MCSCF)=0.0 and MPLEVL=2 (#509) (G. Schoendorff of Iowa State University)
- Fix GMRES normed coefficient in mp2grd2 when vector is zeroed (#508) (I. S. Gerasimov of A.N.Nesmeyanov Institute of Organoelement Compounds of Russian Academy of Sciences)
- Correct NOROT option with FULLNR (#495) (G. Schoendorff of Iowa State University)
- Address issues in subroutine efpdyn in efdrvr.src (#498) (P. Xu and T. Sattasathuchana of Iowa State University)
- NEBPATH attribution update (#464) (S. Koseki of Osaka Prefecture University, Japan)
- Fix type conversions in rimp2grd.src (#439) (I. S. Gerasimov of A.N.Nesmeyanov Institute of Organoelement Compounds of Russian Academy of Sciences)
- Uninitialized variable DFDR in subroutine B97X (#474) (C. Friedl of Johannes Kepler University Linz)
- Address issues found with GNU AddressSanitizer (#397) (#518) (I. S. Gerasimov of A.N.Nesmeyanov Institute of Organoelement Compounds of Russian Academy of Sciences)

Software related changes:

GAMESS
- Common block to modules conversion for CCDDI (#456) (T. Harville of Iowa State University) (#552) (S. Leang of EP Analytics)
- Testing improvements
    - (#523) (#554) (#558) (S. Leang of EP Analytics)
    - (#442) (#524) (#527) (I. S. Gerasimov of A.N.Nesmeyanov Institute of Organoelement Compounds of Russian Academy of Sciences)
- Add -finit-real and -finit-integer to default FPE flags (#482) (C. Friedl of Johannes Kepler University Linz)
- Remove FLUSHONWRITE for AIX (#425) (I. S. Gerasimov of A.N.Nesmeyanov Institute of Organoelement Compounds of Russian Academy of Sciences)
- Replace "rm -rf" with explicit commands (#423) (#424) (C. Friedl of Johannes Kepler University Linz)
- Switch from Travis-CI to GitHub actions (#457) (#473) (#485) (I. S. Gerasimov of A.N.Nesmeyanov Institute of Organoelement Compounds of Russian Academy of Sciences)
- Remove un-used architectures from unport.src (#409) (I. S. Gerasimov of A.N.Nesmeyanov Institute of Organoelement Compounds of Russian Academy of Sciences)
- Software infrastructure modernization and improvements
    - (#510) (G. Schoendorff of Iowa State University)
    - (#489) (#497) (#503) (S. Leang of EP Analytics)
    - (#470) (#471) (#493) (C. Friedl of Johannes Kepler University Linz)
    - (#405) (#430) (#435) (#436) (#455) (#460) (#486) (#492) (#520) (#522) (#533) (#535) (I. S. Gerasimov of A.N.Nesmeyanov Institute of Organoelement Compounds of Russian Academy of Sciences)
- Mark Fortran-Free-Form in .gitattributes (#428) (C. Friedl of Johannes Kepler University Linz)

LIBCCHEM
- Software infrastructure modernization and improvements
    - (#528) (J. L. Galvez Vallejo of Iowa State University)

Notice: 32-bit build support will be removed from GAMESS after the September 30, 2021 R2 public release.

## September 30, 2020 R2 Public Release

GAMESS
- OpenMP-threaded RI-CC method (#413) (D. Datta of Iowa State University)
    - This method is currently only available for the OpenMP-threaded version of GAMESS.
    - Only CCSD and CCSD(T) calculations can run using the RI approximation for two-electron repulsion integrals and a hybrid MPI/OpenMP parallelization model
    - See new CCERI flag in $CCINP
    - See new $RICC input group
        - This group controls the calculation of ground state  energies for closed-shell molecules (SCFTYP=RHF) using the coupled-cluster method within the singles and doubles (CCSD) truncation scheme or CCSD augmented with an approximate triples correction: CCSD(T). These calculations run in parallel using a hybrid MPI/OpenMP parallelization model and employ the resolution-of-the-identity (RI) approximation for two-electron repulsion integrals (2ERIs).
- ENTRY statement replacement:
    - in mcqud.src, mcqdpt.src, mcpgrd.src, gugdm.src, gugdm2.src, and gugem.src (#335) (L. Carrington of EP Analytics)

Along with a few minor fixes:

GAMESS
- Set RORB to zero when NROT=0 (#400) (C. Friedl of Johannes Kepler University Linz)
- Fix $GMSDIR/tests/sformas/eth-sfcepa.inp (#391) (J. Mato of University of Colorado Denver)
- Address FPE in standard exams related to optimization level of dgeev.src and dftgrd.src when compiled with PGI 19.10 (#393) (C. Friedl of Johannes Kepler University Linz)

Software related changes:

GAMESS
- Performance regression metrics (#415)  (S. Leang of EP Analytics)
    - Updated parser definition
        - as a result several validation files had to be regenerated
    - New flag for checkgms.py `--test_path TEST_PATH` to specify a directory path for checking
    - New flag `-m {validation, regression, both}, --mode {validation,regression,both}` to allow users to select between performing validation testing, regression testing, or both
    - Flags for `-m regression, --mode regression`
        -  `-a ARCH, --arch ARCH  - architecture string used for performance regression`
        - `-s SYSTEM, --system SYSTEM - HPC system string used for performance regression`
        - `-n NCPUS, --ncpus NCPUS - number of GAMESS compute processes for performance regression`
        - `-g NACCELERATORS, --naccelerators NACCELERATORS - number of accelerators for performance regression`
    - Performance data is stored with .perf extension
    - Example usage: `./checkgms.py --mode regression --system "Summit" --arch "Power9 with V100" --ncpus 40 --naccelerators 1`
    - Running `./checkgms.py --mode regression` without any additional flags will show you what baselines are available and what flags need to be passed.
    - Performance tolerance is currently set at +/- 5%

- Production code does not set `GMS_FPE_FLAGS` (#393) (S. Leang of EP Analytics)
- serial-debug support in misc/automation/rungms and tests/runtest.py (#388) (C. Friedl of Johannes Kepler University Linz)
    - tests/runtest.py now has the ability to print the output of the test runs directly to console via the new argument -i or --print_to_stdout.
    - For a build with the option `GMS_DDI_COMM=serial-debug`: When using tests/runtest.py -i misc/automation/rungms now automatically runs the tests under the control of gdb. In the case of using tests/runtest.py without -i, the tests are executed normally (without gdb) and the output is written to the .log files in their standard location.
- Removal of actvte.x (#383) (I. S. Gerasimov of A.N.Nesmeyanov Institute of Organoelement Compounds of Russian Academy of Sciences)
- Use standard functions in place of compiler extensions (#334) (I. S. Gerasimov of A.N.Nesmeyanov Institute of Organoelement Compounds of Russian Academy of Sciences)
    - IAND, IOR, IEOR
    - DATE_AND_TIME, SYSTEM_CLOCK, CPU_TIME
    - GET_ENVINROMENT_VARIABLE
    - FLUSH
    - Instead of C calls MEMGET and MEMREL the static allocatable array is used (defined in FAST_MEMORY in unport.src) with the hack (function get_ptr_addr in unport.src) for pointer arithmetic (Fortran 2003).
- Add `--pre` and `--post` flags to runtest.py (#389) (S. Leang of EP Analytics)
- IBM64 Summit modifications (#406) (S. Leang of EP Analytics)
    - 64-bit ESSL support
    - Force auto deallocation
- Build improvements (#393) (S. Leang of EP Analytics) (F. P. Rotzinger):
    - Floating point exception flags are not set in production.
    - Add "-Wno-maybe-uninitialized" flag for GNU compiler builds in production.
- Array handling mod_shell_tools for IBM XL compilers (#373) (I. S. Gerasimov of A.N.Nesmeyanov Institute of Organoelement Compounds of Russian Academy of Sciences)

## June 30, 2020 R1 Public Release

GAMESS
- The nuclear-electronic orbital (NEO) method package in GAMESS has been updated to include NEO-DFT capability with electron-proton correlation methods [1,2]. The epc17 functional is now available as a treatment for electron-proton correlation in a NEO-DFT calculation. The capability to calculate proton density on a set of grid points provided by the user has also been added. The documents NEO.DOC and NEOINP.DOC in the qmnuc/ folder have been updated to reflect these changes, but a brief summary follows. Electron-proton treatment is controlled by the `NEODFT` flag. Currently, the only available treatment is `NEODFT=CS4` for epc17. Parameters for the functional [1,2] are specified by `YANGA`, `YANGB`, and `YANGC`. The flag `IPRDEN=1` designates a run that will calculate proton density at specific Cartesian points (given by the user in `$POINTS`, see GAMESS input manual). All new input flags belong to the `$NEO` group. (#183) (P. Schneider of Yale University / Hammes-Schiffer Group)
    - [1] Y. Yang, K. R. Brorsen, T. Culpitt, M. V. Pak, S. Hammes-Schiffer, Development of a practical multicomponent density functional for electron-proton correlation to produce accurate proton densities, J. Chem. Phys., 147, 114113 (2017).
    - [2] K. R. Brorsen, Y. Yang, S. Hammes-Schiffer, Multicomponent density functional theory: impact of nuclear quantum effects on proton affinities and geometries, J. Phys. Chem. Lett., 8, 3488-3493 (2017).
- MPI/OpenMP implementation for QM-EFP2 energy calculation. More specifically, the implementation is for QM-EFP2 exchange repulsion electron repulsion integrals (ERI) and corresponding Fock contribution construction, as well as QM-EFP2 Schwartz inequality screening subroutine.  A QM-EFP2 test folder is created inside `$GMSDIR/tests/solvent` (look in `$GMSDIR/tests/solvent/qm-efp2/exrep/mpiomp`).  The pre-existing `$GMSDIR/tests/exrep` folder was moved inside `$GMSDIR/tests/solvent`.  (#207) (P. Xu and T. Sattasathuchana of Iowa State University)
- QM-EFP improvements (#216) (P. Xu and T. Sattasathuchana of Iowa State University)
    - Use dynamic memory allocation for QM-EFP2 dispersion (EFAIDISP).
    - Post-HF methods, MP2 and Coupled-Cluster, now are interfaced with EFP2
    - New tests added to `$GMSDIR/tests/exrep`
- FMO 5.4 (#233) (D. G. Fedorov of the National Institute of Advanced Industrial Science and Technology Japan)
    - New Features
        - Multiple-level GDDI parallelization (A. Gaenko, S. Pruitt).
        - Analytic FMO-PCM Hessian (H. Nakata).
        - Long-range corrected DFTB (V. Q. Vuong and Y. Nishimoto).
        - New (so called partial) screening model for PIEDA/PCM.
        - Minimum energy crossing with FMO (D. Kaliakin).
        - Molecular electrostatic potential (MEP) with PCM (V. Mironov).
    - Improvements
        - Improved error diagnostic for host lookup in DDI.
        - Acceleration improvements for FMO3.
        - Dimer data stored in FMO for later reuse.
        - Some improvements in SCC of DFTB (Y. Nishimoto).
        - Some non-minimal basis sets are enabled for HF-3c (K and Ca).
        - Acceleration for the old (local) screening in FMO/PCM.
        - FMO-MP2 runs enabled for fragments without valence electrons (Ca2+ etc).
        - RUNTYP=OPTIMIZE for FMO now does not stop when PROJCT fails.
        - Parallelized nuclear contributions to MEP.
        - FMO-DFT and HF-3c with PIEDA now print all components clearly.
        - Cleanup of FMO output.
    - Fixes
        - Array declaration issue in FMO (a compiler error for some compilers).
        - Trap undeclared SCFTYP in open-shell FMO.
        - FMO/FDD/PCM had a bug with RESPPC.
        - MODGAM=8 did not work with FMO-DFTB/PCM (Y. Nishimoto).
        - Coordinate units were wrong for printing tesserae data in FMO/PCM.
        - RI-MP2 ignored SCSPT option.
        - FMO/AP/PCM has a major memory bug.
        - Total binding for modmol options did not add solvent screening in FMO/PCM.
        - FMO-MP2/AP had a major bug for dimer calculations.
        - Divisions by zero in DFTB are removed. (Y. Nishimoto).
        - FMO-UDFT incorrectly calculated the electron count.
        - FMO-TDDFT had a major bug for iexcit(2)=1.
        - TDDFT overwrote exchange fraction.
        - Pople integrals used undeclared variables.
        - LAPACK KDIAG=5 was fixed.
- FMO 5.4 patch (#268) (D. G. Fedorov of the National Institute of Advanced Industrial Science and Technology Japan)
    - Patch DDI to stop warnings during compilation (define prototypes).
    - Enable the first row of transition metals in HF-3c.
    - Describe FMO/FRET in the manual.
    - Remove symlink for DFTB parameters folder
    - Rename DFTB parameters folder
    - Adjust DFTB test inputs
- Additional of revised M11 functional for DFT energy and gradients `DFTTYP=REVM11`. (#244) (I. S. Gerasimov of A.N.Nesmeyanov Institute of Organoelement Compounds of Russian Academy of Sciences)
    - [REVM11 reference](https://pubs.acs.org/doi/10.1021/acs.jpca.8b11499)
- Interface to LIBXC for ground-state DFT  (#274) (I. S. Gerasimov of A.N.Nesmeyanov Institute of Organoelement Compounds of Russian Academy of Sciences)
    - Requires GCC compiler > 5.0
    - Build process
        - `./config` (answer yes to LIBXC build option)
        - `./tools/libxc/download-libxc.csh`
        - `make ddi`
        - `make libxc -j$(nproc)`
        - `make modules`
        - `make -j$(nproc)`
    - See `$LIBXC` for usage in docs-input.txt
- MPI and MPI+OpenMP implementation of QM-EFP2 electron repulsion integral schemes (Rotated axis (sp and spd) and ERIC) (#310) (P. Xu and T. Sattasathuchana of Iowa State University).
    - QM-EFP ERI scheme can be controlled by `INTTYP` in `$CONTRL`:
        - `INTTYP=best (default)`
        - `INTTYP=rotaxis`
        - `INTTYP=Eric`
        - `INTTYP=rysquad`
    - Tests
        - MPI implementation: `$GMSDIR/tests/solvent/qm-efp2/exrep/`
            - `h2o_acetone_efpqm.rys.inp`
            - `h2o_acetone_efpqm.rotaxis.inp`
            - `h2o_acetone_efpqm.eric.inp`
        - MPI+OpenMP implementation: `$GMSDIR/tests/solvent/qm-efp2/exrep/mpiomp/`
            - `h2o_meoh_qmefp.rys.inp`
            - `h2o_meoh_qmefp.rotaxis.inp`
            - `h2o_meoh_qmefp.eric.inp`
- Spin-flip Occupation Restricted Multiple Active Space (ORMAS) energy and gradients (#316) (J. Mato of Iowa State University)
    - `CITYP=SFORM` or `CITYP=SFDET` when ALDET code is used.
    - Test inputs in `$GMSDIR/tests/sformas`
    - J. Mato, M. S. Gordon, A general spin-complete spin-flip configuration interaction method, Phys. Chem. Chem. Phys., 20, 2615-2626 (2018)
    - J. Mato, M. S. Gordon, Analytic Gradients For The Spin-flip ORMAS-CI Method: Optimizing Minima, Saddle Points, And Conical Intersections, J. Phys. Chem. A, 123, 1260-1272 (2019)
- Nudge Elastic Band Release 3 (#319) (S. Koseki of Osaka Prefecture University, Japan)
    - See `tools/neb/README.NEB.md` for more information
- Updates to OpenMP threaded RI-MP2 code (#330) (#361) (B. Q. Pham of Iowa State University)
    - Converts energy code from F77 to F90
    - Optimize (E)FMO/RIMP2 performance for systems there the biggest dimer fits into single process memory
- EFMO dispersion tweaks (#372) (A. Gunina of Iowa State University)
- Initial interface for LibAccInt (#378) (B. Westheimer of Iowa State University)

LIBCCHEM
- Decoupling of LIBCCHEM RI-MP2 and DF-HF. The LibCChem RI-MP2 code can be run using different GAMESS HF drivers. (#371) (D. Poole of Iowa State University)

Along with a few minor fixes:

GAMESS
- Resolve floating point exceptions (FPE) in GAMESS:
    - Numerical exception in DFTB: tests/dftb/parallel/exam04.inp (#192) (Y. Nishimoto of Kyoto University)
    - Overflow when running `$GMSDIR/tests/dftb/raman.inp` (#220) (C. Friedl of Johannes Kepler University Linz)
    - Overflow in `RPBE` and `RPBEX` (#224) (C. Friedl of Johannes Kepler University Linz)
    - FPE in quanpo.src by clearing `ENALL(1:42)` before `GSUMF(2410,ENALL,42)` (#228) (C. Friedl of Johannes Kepler University Linz)
    - Issues uncovered with `-finit-real=snan` (#245) (#248) (#249) (#250) (#278) (C. Friedl of Johannes Kepler University Linz)
    - FPE in H0 (impacts:  M08-HXC, M08-HX, M08-SOC, M08-SO, M11-C, M11, M11-LC, M11-L, MN12-L, MN12-SX, MN15, MN15-L, REVM11) (#285) (C. Friedl of Johannes Kepler University Linz)
    - Issues associated with un-initialized variables (#347) (#342) (C. Friedl of Johannes Kepler University Linz)
- Fix scaling parameters in tddfun.src (#236) (S. Leang of EP Analytics)
- Fix segfaults for Intel 18/19 in efmo.src (#201) (I. Rostov of National Computational Infrastructure at the Australian National University)
- Additional TD-DFT tests added to `$GMSDIR/tests/tddft/parallel` (#237) (I. S. Gerasimov of A.N.Nesmeyanov Institute of Organoelement Compounds of Russian Academy of Sciences)
- Fix H, I orbital output if `AIMPAC=.T.` (#211) (I. S. Gerasimov of A.N.Nesmeyanov Institute of Organoelement Compounds of Russian Academy of Sciences)
- Adjust variable array sizing found in basg3x.src, dftbhs.src, ecp.src, efmogrd.src, and qfmm.src (#253) (I. S. Gerasimov of A.N.Nesmeyanov Institute of Organoelement Compounds of Russian Academy of Sciences)
- Common block alignments (#246) (I. S. Gerasimov of A.N.Nesmeyanov Institute of Organoelement Compounds of Russian Academy of Sciences)
- Increase parameter limits in mx_limits.src: `mxcpuefp`, `mxefmopts`, `mxefmoppts` (#256) (S. Leang of EP Analytics)
- Resolve issue with FMO Hessian calculations (#267) (V. Mironov of Lomonosov Moscow State University)
- Resolve issue with DC-RHF when OpenMP is enabled. (#321) (C. Bertoni of Argonne National Laboratory)
- Resolve OOB issues in dftbsk.src, efinp.src, and gmcpt.src (#327) (T. Sattasathuchana of Iowa State University)
- Resolve issues with MC-PDFT code (#312) (#356) (A. O. Lykhin of University of Minnesota Twin Cities)
- Resolve UHF/UDFT wavefunction output (#331) (I. S. Gerasimov of A.N.Nesmeyanov Institute of Organoelement Compounds of Russian Academy of Sciences)
- Resolve ATTRIBUTES directive discrepancy in mod_1e_primitives.src (#328) and omp_int1.src (#345) (C. Bertoni of Argonne National Laboratory)
- Fix memory allocation issue in inputb.src (#351) (I. S. Gerasimov of A.N.Nesmeyanov Institute of Organoelement Compounds of Russian Academy of Sciences)
- Re-enabled kinetic bond order and bond order printing during a QUAO procedure (#376) (J. L. Galvez Vallejo of Iowa State University)

LIBCCHEM
- No fixes to report.

Software related changes:

GAMESS
- Ability to write arrays as numpy files (#193) (C. Friedl of Johannes Kepler University Linz)
- Update **bin/create-install-info.py** to support `--ddi_comp==serial-debug` (#209) (C. Friedl of Johannes Kepler University Linz)
- Update rungms to check of `LD_LIBRARY_PATH` exists in the user's environment (#212) (C. Friedl of Johannes Kepler University Linz)
- ENTRY statement replacement:
   - in dftgrd.src (#215) (L. Carrington of EP Analytics)
   - in gmcpt.src (#370) (C. Bertoni of Argonne National Laboratory)
- Define `MATHLIBS` for `cray-xc` targets (#227) (C. Bertoni of Argonne National Laboratory)
- Update documentation for DAREAD and IODA (#221) (C. Bertoni of Argonne National Laboratory)
- Add .gitattributes for proper syntax highlighting in .src files (#230) and C-shell scripts (#282) (#311) (C. Friedl of Johannes Kepler University Linz)
- Remove CR LF endings in tddxce.src (#238) (I. S. Gerasimov of A.N.Nesmeyanov Institute of Organoelement Compounds of Russian Academy of Sciences)
- Intel 2019 support:
    - Compilers (#255) (S. Leang of EP Analytics)
    - Intel MPI (#279) (C. Bertoni of Argonne National Laboratory)
- Add DDI stubs for serial and serial-debug compilations. (#263) (S. Leang of EP Analytics)
- Address linking issue with OpenMP builds of GAMESS (#290)  (C. Bertoni of Argonne National Laboratory)
- Adjust line-endings (dos2unix) for several files (#297) (S. Leang of EP Analytics)
- path/to/file.inp support for rungms (#307) (I. S. Gerasimov of A.N.Nesmeyanov Institute of Organoelement Compounds of Russian Academy of Sciences)
- Remove EQUIVALENCE statements in GAMESS source files (#308) (I. S. Gerasimov of A.N.Nesmeyanov Institute of Organoelement Compounds of Russian Academy of Sciences)
- Common block to modules for HERMIT, WERMIT, and XYZCHI (#284) (I. S. Gerasimov of A.N.Nesmeyanov Institute of Organoelement Compounds of Russian Academy of Sciences)
- Common block to modules to non-threaded DFT code (#318) (I. S. Gerasimov of A.N.Nesmeyanov Institute of Organoelement Compounds of Russian Academy of Sciences):
     - CM05; CM06; CM08; CSOGGA; FUNLIB; METGGA; CBECKE; SCTAX; CLYPC; CVWNC; CP86CF; CPBE0; CPW91C; CX3LYP; CPW91L; CTPSSH; SCPBES; CPKZB; SCAEDF; SCAOPC; SCLGIL; SCPFRE; SCPZ81; SCRPBE; B97TYP; SLPBEC, NLRC; DFTDH; DFTCAM.
     - DFTEXC was split to two modules: XDERIX and CDERIX.
- Patches from NCI ANU for issues encountered using Fujitsu fortran compiler (#258) (I. Rostov of National Computational Infrastructure at the Australian National University)
- Makefile changes (#322) (I. S. Gerasimov of A.N.Nesmeyanov Institute of Organoelement Compounds of Russian Academy of Sciences)
- Resolve build issues for mac64 targets (#323) (E. Guidez of University of Colorado Denver)
- Update Travis-CI to used GNU 4.9.4 as baseline GNU compiler version for GAMESS testing (#324) (S. Leang of EP Analytics)
- Update GCC build support (#338) (I. S. Gerasimov of A.N.Nesmeyanov Institute of Organoelement Compounds of Russian Academy of Sciences)
- Enable the running of serial tests calculations in parallel using multi-threading (#341) (I. S. Gerasimov of A.N.Nesmeyanov Institute of Organoelement Compounds of Russian Academy of Sciences)
- Clean up config script (#357) (I. S. Gerasimov of A.N.Nesmeyanov Institute of Organoelement Compounds of Russian Academy of Sciences)
- Build fixes/enhancements for ibm64 target (Summit)
    - Use of environmental variables in config (#360) (S. Leang of EP Analytics)
    - Resolve compiling issues (#333) (#326) (S. Leang of EP Analytics) (#374) (I. S. Gerasimov of A.N.Nesmeyanov Institute of Organoelement Compounds of Russian Academy of Sciences)
    - Address field separator (#301) (S. Leang of EP Analytics)
- Build enhancements (#346) (S. Leang of EP Analytics)
    - Added `GMS_HPC_SYSTEM_TARGET` to enable conditional paths in GAMESS build and launch scripts for specific systems (e.g., generic (default), onyx, theta, summit, hokulea, qcengine)
    - Added `--system_target` flag to `create-install-info.py`
    - Added `my_ipcrm` script to `$GMSDIR/bin` to help clean up user's semaphores (do not use when running multiple GAMESS calculations on same compute nodes)
- Software documentation changes (#379) (S. Leang of EP Analytics):
    - INTRO.DOC -> docs-intro.txt
    - INPUT.DOC -> docs-input.txt
    - TESTS.DOC -> docs-tests.txt
    - REFS.DOC -> docs-references.txt
    - PROG.DOC -> docs-prog.txt
    - IRON.DOC -> docs-hardware.txt
    - Attribution format has changed.  Split list for individuals and their affiliated organization listed in alphabetical order.  GAMESS banner has also slightly changed.

LIBCCHEM
- Link C++ bindings to LIBCCHEM builds using OpenMPI and removed `-lcublas_device` flag from linking. (#367) (K. Keipert of NVIDIA)
- LIBCCHEM build improvements (#368) (D. Poole of Iowa State University)
    - Remove `-march=native` flag in CMake build script
    - Stub additional routines in ga.src

## September 30, 2019 R2 Public Release

GAMESS
- Multiconfiguration Pair-Density Functional Theory (MCPDFT). MCPDFT combines MCSCF and DFT approaches [1,2]. In MC-PDFT, an MCSCF wave function is used to calculate the kinetic energy, electron–nuclear attraction the Coulomb classical electron-electron repulsion contributions to the total energy; the remaining exchange−correlation energy comes from an on-top density functional. Various attempts have been made to combine MCSCF and DFT methodology to get an affordable yet accurate method to treat strongly correlated systems. Most of the previous methods encountered two problems: 1) double counting of electron correlation, as MCSCF contains some amount of dynamic correlation which is counted again if one simply adds a correlation energy computed from a correlation functional to the MCSCF energy, and 2) the incompatibility of available functionals developed for Kohn-Sham theory with the spin densities of a multiconfigurational wave function. An attempt has been made in MC-PDFT method to avoid these two problems. The first problem is avoided by using only the kinetic, electron-nuclear and Coulomb interaction terms of the MCSCF energy, with the remaining part of the energy calculated through on-top energy functional. The second problem is solved by introducing translated Kohn-Sham density functionals that depend on the on-top pair density and the total density in contrast to the alpha and the beta density dependence of  KS-DFT functionals [1,3]. See `PDFTYP` in the `$CONTRL` group. (#182) (Andrew Sand and Prachi Sharma of University of Minnesota)
    - [1] Li Manni, G.; Carlson, R. K.; Luo, S.; Ma, D.; Olsen, J.; Truhlar, D. G.; and Gagliardi, L.  J. Chem. Theory Comput. 2014, 10(9), 3669-3680.
    - [2] Gagliardi, L.; Truhlar, D. G.; Li Manni, G; Carlson, R. K.; Hoyer, C. E.; Bao, J. L. Acc. Chem. Res. 2016, 50(1), 66-73
    - [3] Carlson, R. K.; Truhlar, D. G.; Gagliardi, L. J. Chem. Theory  Comput. 2015, 11( 9), 4077-4085.
- New keyword TRNSD in `$TDDFT` group that allows property computations on the transition density for the transition from the ground state to the state `IROOT`.  This allows one - for example - to calculate PDCs from the electric potential of the transitions density as proposed by the TrESP method [1].  (#174) (C. Friedl of Johannes Kepler University Linz)
    - [1] Madjet, M. E.; Abdurahman, A.; Renger, J. Phys. Chem. B 2006, 110 (34), 17268–17281.

Along with a few minor fixes:

GAMESS
- Fix parallel FMO gradient issue (#202) (C. Bertoni of Argonne National Laboratory)
- Fix issue in shared-memory Fock matrix code (#198) (#204) (V. Mironov of Lomonosov Moscow State University)
- Fix SLB bug in FMO OMP (#197) (V. Mironov of Lomonosov Moscow State University)
- Fix floating point execption trap in EFMO for ./tests/efmo/exam03-serial (#195) (C. Bertoni of Argonne National Laboratory)
- Change `TDTYP` to `TDDFT` in `$FMO` group. (#180) (C. Bertoni of Argonne National Laboratory)

Software related changes:

GAMESS
- Updated Windows 64 support (#190) (S. Leang of EP Analytics)
- New `GMS_DDI_COMM` options: `serial` and `serial-debug`.  Useful for serial calculations and/or debugging. (#189) (#174) (C. Friedl of Johannes Kepler University Linz)
- Aggregation of array limits into Fortran module `mx_limits`. (#178) (C. Bertoni of Argonne National Laboratory)

LIBCCHEM
- Improve IBM Power 64 (Summit) support in LIBCCHEM. (#172) (D. Poole of Iowa State University)
- Modularization of LIBCCHEM build process. (#169) (D. Poole of Iowa State University)

## June 30, 2019 R1 Patch 1 Release

Along with a few minor fixes:
GAMESS
- Fixed GAMESS release version output
- Fixed memory leak when using Intel 2017 & 2018 and running large FMO jobs (#173) (V. Mironov of Lomonosov Moscow State University)
- Fixed PGI 19.4 Community Edition compiler support with new DFT code (#170) (V. Mironov of Lomonosov Moscow State University)
- Fixed math library issue with efpmodule when compiling on IBM systems (Summit and Ascent) (#167) (D. Poole of Iowa State University)
- Fixed numerical discrepancies with GGA energies in new DFT code (#165) (V. Mironov of Lomonosov Moscow State University)

## June 30, 2019 R1 Public Release

GAMESS
- Conversion of EFP related common blocks to modules (#155) (P. Xu and of T. Sattasathuchana of Iowa State University)
- PGI toolchain support (#152) (D. Poole of Iowa State University)
    - The -tp=px flag was added to all files when pgfortran is requested, to remove architecture-specific PGI-compiled binaries
    - DDI is now compiled with the PGI C compiler when pgfortran is requested
    - PGI LAPACK support was added to lked
- HF, DFT, FMO OpenMP Update (#151) (V. Mironov of Lomonosov Moscow State University)
    - Multithreaded FMO
    - Multithreaded one-electron integrals
    - Multithreaded PCM (partial)
    - DFT code improvements
    - DFT code improvements:
        - new driver for grid DFT XC integration with multithreading support;
        - support for SSF (by Stratmann, Scuseria and Frisch) molecular grid partitioning scheme (in addition to Becke's one);
        - various partitioning functions: Becke's, SSF, Erf1 (from NWChem), smoothstep;
        - use BLAS to update Kohn-Sham matrix;
        - option to ignore grid weight contribution to the DFT gradient;
        - two new types of radial quadrature: Mura-Knowles and Treutler-Ahlrichs;
        - few new spherical grids of various symmetry (octahedral with and without inversion, icosahedral) for tests;
        - disable grid switching when $DFT switch=0.0 instead of crashing;
        - cleanup and bug fixes in legacy DFT code, change layout of some arrays from C to Fortran style;
    - Various other improvements:
        - thread-safe ERIC integrals;
        - cleanup and bug fixes in older OpenMP code;
        - remove unnecessary matrix copying in eigen.src when using LAPACK;
- ARM Clang compiler support for linux32 and linux64 targets (#149) (C. Bertoni of Argonne National Laboratory)
- New DFT functionals: APF, HCTH, HCTH14, HCTH76, HCTH407P, HCTH407P, CAMQTP00, CAMQTP01, RM06X, RM06C, REVM06, RM06-LX, RM06-LC, REVM06-L. See $DFT in INPUT.DOC (#148) (I. S. Gerasimov of A.N.Nesmeyanov Institute of Organoelement Compounds of Russian Academy of Sciences)
- AMD software suite support: BLIS BLAS library, LibFLAME LAPACK library, AMD Optimizing C Compiler (AOCC), LibM mathematical functions library. (#143) (D. Poole of Iowa State University)
- Nudged Elastic Band (NEB) Methods. See $NEB in INPUT.DOC (#145) (N Sawada, M. Haruta, T. Asada, S. Koseki of Osaka Prefecture University, Japan)
- QM/EFP2 update (coulomb, exchange-repulsion and Dispersion terms) (#112) (P. Xu of Iowa State University)
- Improve CC file handling through 2GB file chunking. (#111) (J. L. Galvez Vallejo of Iowa State University) (J. Boatz of Air Force Research Laboratory)
- Enable EFMO/CCSD(T) calculations (#108) (C. Bertoni of Argonne National Laboratory)

LIBCCHEM
- PGI toolchain support (#152) (D. Poole of Iowa State University)
    - PGI compiler can now be used to compile LIBCCHEM
- Closed- and Open-shell RI-energy improvements for CPU and GPU. (#122) (L. Roskop of Iowa State University)
- Update LIBCCHEM to use C++11 standard (#121) (J. L. Galvez Vallejo of Iowa State University)
- CMake capabililities for LIBCCHEM (#110) (L. Roskop of Iowa State University)

GAMESS + LIBCCHEM
- Singularity support added to launch Singularity container builds of GAMESS/GAMESS+LIBCCHEM. Run ./config and select singularity as target. (S. Leang of EP Analytics)

Software related changes:
- Improve GAMESS build process.  Users can now perform the following commands to build GAMESS after running ./config:
    - make ddi
    - make modules
    - make -j <ntasks>


## Sept 30, 2018 R3 Public Release

New public release with:

GAMESS

- New general features
    - Analytic Raman intensities for DFTB (Y. Nishimoto).
    - Replica-exchange umbrella sampling (S. Ito).
    - Subgroups in GDDI (A. Gaenko, S. Pruitt)
    - Density of states (DOS).
    - Assign blocks of tasks for DLB (instead of individual tasks).
    - DFTB input simplified (can omit parameter file names).
    - CMIRS Version 1.1 (D. Chipman)
    - New version of QMCPACK-GAMESS interface that includes AFO-type extension of the QMC-EFMO method to fragmentations across covalent bonds
    - OpenMP threaded RI-MP2 gradients (B. Pham)
        - Set `setenv GMS_OPENMP true` in install.info
        - Set `GMS_OPENMP = true` in Makefile
        - Set `$INTGRL INTOMP=1 $END` in input file
        - Set `$MP2 CODE=RIMP2 $END` in input file
    - Improved harmonic restraints
        - Allow angle and dihedral restraint types
        - See IHMCON and SHMCON in INPUT documentation

- New FMO features
    - FMO-DFTB/AFO (Y. Nishimoto)
    - ESP-PC approximation for FMO Hessians (H. Nakata).
    - FMO/LCMO(X) and FMO3/F(X).
    - PIEDA/DFT and PIEDA/DFTB.
    - Fluctuation analysis for FMO/MD.
    - FMO/SMD solvation model.
    - Numerical gradient for FMO.
    - Any number of atoms in FMO/MD.
    - Solvent dipole contributions for FMO/PCM.
    - Direct matrix inversion in FMO/PCM.
    - Charge difference printed in FMO-DFTB showing convergence.

- Changes of option defaults
    - The defaults of DAMPEX and HUBDER for DFTB3 were changed for 3ob-3-1 parameters.
    - Umbrella sampling MD's option IUSTYP=1 was changed to agree with other programs; the old way is available with IUSTYP=-1.

LIBCCHEM

- New features (L. Roskop)
    - Kinetic energy and electron nuclear attraction derivatives added
    - RYSQ library: libint is no longer required for density-fitting-perturbation-theory (DF-PT) gradients
    - DF-ZAPT added (multi-GPU and CPU)
    - DF-MP2 now has multi-GPU support
    - closed- and open-shell DF-SCF added (multi-GPU and CPU)
    - New $LCCHEM subgroup added to for LIBCCHEM control

Along with a few minor fixes:

GAMESS

- For systems without electrons (Na+ etc), DFTB gradient had a bug.
- DFTB did not use the origin in computing dipoles.
- FMO1/PCM had a bug with IDISP=1.
- FMO/FDD/PCM did not work with OPTFMO.
- FMO/F did not work with PCM.
- In FMO/PCM[1(n)], virial components were not computed properly.
- Quadrupole and octupole moments had a parallel bug.
- FMO0 generated wrong F30 file.
- Stability analysis did not work in multilayer FMO.
- FMO/PCM gradient had a bug with ESP-PC.
- Charges were not reused in FMO-DFTB/MD.
- F40, stored in DDI memory, had a bug with Raman.
- In DFTB, d-orbitals resulted in NaN.
- FMO/FDD/PCM had a memory bug in the gradient.
- MD did not print incremental energy for JEVERY>10.
- FMO/SA did not work properly for correlated methods.
- SCS-MP2 energies were wrong for FMO-MP2/AB runs.
- RUNTYP=FMOHESS did not work with DDI-stored F40.
- DFTB Hessian did not work with DISP=SK(HP).
- In the Jan, 2017 changes, FMO gradient may have been affected.
- EOM/CCSD3A was not checked properly.
- FMO-MP2 erroneously aborted by checking consistency in $DATA.
- HW ECP did not work with FMO.
- Better memory management for ORMAS runs (J. Ivanic).

Software related changes:

- VB2000 is not longer compiled by default. To enable VB2000 modify compall and lked.

## Aug 02, 2018 R2 Public Release

New public release with:

- New Michigan State University Coupled-Cluster Theory Methods: CCSD3A/CCT3
  - You must run config and answer "yes" to the CCSD3A/CCT3 build question.
  - CCTYP (works for both RHF and ROHF references)
    - CCSD3A
      - CC calculation with singles, doubles, and active-space triples, designated as CCSDt. When all orbitals in the MO basis are set to be active, CCSDt becomes full CCSDT.  See $CCINP.
    - CCT3
      - In addition to CCSDt (see CCTYP=CCSD3A), compute the CC(t;3) noniterative energy correction due to the missing triples not captured by CCSDt.
  - See $GMSDIR/tests/cc/cc-ccsd3a-\*.inp and $GMSDIR/tests/cc/cc-cct3-\*.inp for example input files.
  - More information is provided in INPUTS.DOC
- New Minnesota Density Functionals: MN12/MN15
  - Available for ground-state energy and gradients.
  - Available for excited-state energy.
  - DFTTYP
    - MN12-L
    - MN12-SX
    - MN15
    - MN15-L

Along with a few minor fixes:

- Chunking I/O was added to locpol.src to address large DAF files during MAKEFP runs. (#61)
- Minimize OMP warnings (#67)
- Fixed DGEMM handing in efpaul.src

Software related changes:

- Added doxygen config file for developers (requires doxygen v.1.8.15+).
  - cd $GMSDIR
  - doxygen Doxyfile
- Cavium ThunderX2 (ARM64) support through aarch64 conditionals during compilation.
  - Uses 'linux64' target.  Conditionals within build scripts will pick up on prescence of aarch64 architecture.
- Running 'make' after successfully creaing an 'install.info' will now create the 'object' directory if it does not already exist.
- Added gfortran 5.5 support in config.
- Changed the maximum number of SMP procs (affects 'sockets' builds of GAMESS)
- Build alternative to 'config' or 'compall' is provided in $GMSDIR/bin/create-install-info.py (requires Python 3 and jinja2 module).
  - cd $GMSDIR
  - bin/create-install.info.py --help
    - Simplest build example: linux64 target, gfortran 4.8.5, use built in GAMESS BLAS library, use sockets for DDI communication:
      - bin/create-install-info.py --fortran_version=4.8.5
      - make ddi
      - make -j `nproc`
- Alternative to $GMSDIR/rungms provided in $GMSDIR/machines/xeon_phi/rungms.interactive
  - Copy into your $GMSDIR
  - Edit 'target' in rungms.interactive
  - Use rungms.interactive as you would rungms
- Improved BLAS handling in LIBCCHEM

## Feb 14, 2018 R1 Public Release

New public release with:

- This adds code for the analytical gradient for EFMO.
  - It is turned on with the IEFMOG flag in $FMO in the input file (see INPUT.DOC).
  - Described in C. Bertoni and M. S. Gordon,
    "Analytic Gradients for the Effective Molecular Orbital Method"
    Journal of Chemical Theory and Computation 12 (10), 4743-4767 (2016).
- New target `cray-xc` added.
  - `GMS_PHI` flag in **install.info** changed from `true|false` to `none|knl|knc`

Along with a few minor fixes:

- Modifications to comp, lked and config (#51)(#50):
  - Removed flags EFP_OPENMP, MAKEFP_OPENMP and RIMP2_OPENMP from install.info
  - Threaded codes are now controlled by the flag GMS_OPENMP in install.info

## Nov 11, 2017 R3 Public Release

- OpenMP two-electron code (#28)
  - To build the OpenMP threaded two-electron code set GMS_OPENMP flag to true in install.info and Makefile after running config.
  - See the INTOMP flag within the $INTGRL group for more information.
- OpenMP HF algorithm with Fock matrix shared among threads (#28)
  - See the SHFOCK flag within the $INTGRL group for more information.
- See repository https://github.com/gms-bbg/gamess-papers/tree/mv/SC17 for more information on the above changes

- GAMESS Release Notes:
  - https://gamess.gitbooks.io/gamess-release-notes/content

## Oct 7, 2017 Hot Fix (A)

Along with a few minor fixes:

- Address issue with using 'compall' for building GAMESS (#39)

## Sept 30, 2017 R2 Public Release

New public release with:

- OpenMP threaded RI-MP2 energy (#31)
- OpenMP threaded EFP Charge Transfer (#2)
- OpenMP threaded MAKEFP (#35)
- QMC-EFMO interface between QMCPACK and GAMESS (#36)

Along with a few minor fixes:

- Address issue with SP shells for 2- and 3- center ERI in LIBCCHEM (#29)
- Address issue with multi-pole moment calculations when running in parallel (#25)
- Address issue with using 'compall' for building GAMESS (#39)

Platform updates:

- Added support for GNU compilers up to version 7.2 (Linux, OSX)
