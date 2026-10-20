#!/usr/bin/env python3

import argparse
import multiprocessing
import os
import re
from jinja2 import Environment, FileSystemLoader
import sys


PATH = None
TEMPLATE_ENVIRONMENT = None

# classes


class internal_error(SystemExit):
    # we can treat some internal exceptions differently
    pass


class handle_shortcut(argparse.Action):
    def __init__(self, long_name, nargs=0, **kwargs):
        self.long_name = long_name
        super().__init__(
            nargs=nargs,
            help='set --' +
            self.long_name +
            '=' +
            kwargs['dest'],
            **kwargs)

    def __call__(self, parser, namespace, values, option_string=None):
        parser.parse_args(['--' + self.long_name, self.dest], namespace)

# Helper functions


def update_context_args(context, name, args, relevant_names):
    args = vars(args)
    relevant_names = list(set(relevant_names))
    for n in relevant_names:
        field = args[n]
        if isinstance(field, str):
            if 'path' not in n:
                field = field.replace('oneapi_', 'oneapi-')
                field = field.replace('grace_', 'grace-')
        try:
            check_dir(field, '')
            context[name][n] = field
        except internal_error:
            if n in context[name]:
                # accumulate things like compiler debug flags, etc
                context[name][n] += str(field).lower()
            else:
                context[name][n] = str(field).lower()


def err(msg):
    raise internal_error("Error : {}".format(msg),)


def set_shortcuts(shortcuts_dict, action, parser):
    for long_name, shortcuts in shortcuts_dict.items():
        for short in shortcuts:
            parser.add_argument(
                '--' + short,
                action=action,
                long_name=long_name)


def set_conditional_checks_helper(shortcuts_dict, args):
    args = vars(args)  # get dictionary
    for long_name, shortcuts in shortcuts_dict.items():
        long_name_res = args[long_name]
        for short in shortcuts:
            # don't know why, but the '-' causes a key error
            short = short.replace('oneapi-', 'oneapi_')
            # TODO  see if we can replace None with False in checks
            args[short] = True if long_name_res == short else None


def check_dir(path, name):
    if path is None or not os.path.isdir(path):
        err("Unable to find {} path: '{}'".format(name, path))


def auto_detect_mkl():
    mkl_lib, mkl_include = (False, False)
    try:
        mkl_path = condition_get_from_env(None, 'MKLROOT')
    except internal_error:
        return (mkl_lib, mkl_include)
    try:
        mkl_lib = os.path.join(mkl_path, "lib/intel64")
        check_dir(mkl_lib, '')
    except internal_error:
        mkl_lib = os.path.join(mkl_path, "lib")
        try:
            check_dir(mkl_lib, '')
        except internal_error:
            mkl_lib = False
    try:
        mkl_include = os.path.join(mkl_path, "include")
        check_dir(mkl_include, '')
    except internal_error:
        mkl_include = False
    return (mkl_lib, mkl_include)


def auto_detect_intel_impi():
    mpi_lib = False
    impi_path = None
    try:
        impi_path = condition_get_from_env(impi_path, "I_MPI_ROOT")
        check_dir(impi_path, "impi_path")
    except internal_error:  # auto-detect failed
        return (mpi_lib, impi_path)
    mpi_lib = "impi"
    return (mpi_lib, impi_path)


def auto_detect_fpe(compiler):
    if is_intel(compiler):
        return "-fpe0"
    elif compiler == "gfortran":
        return "-ffpe-trap=invalid,zero,overflow -finit-real=snan \
-finit-integer=-9223372036854775808"
    elif compiler == "nvfortran" or compiler == "pgfortran":
        return "-Ktrap=fp"
    else:
        # TODO: do we err out  or assume the compiler does not have a fpe flag?
        return ""

def is_intel(compiler):
    intel = [
        "ifort",
        "oneapi-ifort",
        "oneapi_ifort",
        "oneapi-ifx",
        "oneapi_ifx"]
    return compiler in intel

# Parse helper functions


def parse_ASAN(args):
    if args.fortran != "gfortran":
        err("ASAN may be applied only with GFortran")
    args.fpe = (args.fpe + " -fsanitize=address").strip()
    if args.debug_link_flags is None:
        args.debug_link_flags = "-fsanitize=address"
    else:
        args.debug_link_flags += "-fsanitize=address"


def parse_coverage(args):
    flags = ""
    if not args.coverage:
        return
    # not in acceptable:
    if not is_intel(args.fortran) and args.fortran != "gfortran":
        err("You can only use the coverage option with intel or gnu fortran")
    elif is_intel(args.fortran):
        flags = "-prof-gen=srcpos"
    else:
        flags = "-O0 -g -fprofile-arcs -ftest-coverage"

    if args.debug_flags:
        args.debug_flags = args.debug_flags + " " + flags
        args.debug_link_flags = args.debug_link_flags + " " + flags
    else:
        args.debug_flags = flags
        args.debug_link_flags = flags


def parse_mkl_info(args):
    if args.mathlib_path is None:
        mathlib = auto_detect_mkl()[0]
    else:
        mathlib = args.mathlib_path
    if mathlib:
        args.mathlib_path = mathlib
        args.mkl = "mkl"
        if args.mkl_verno is None:
            # default
            args.mkl_verno = "12"
            # parse year from args.mathlib_path
            year = re.search(r"201\d", args.mathlib_path)
            if year:
                year = int(year.group(0))
                if year == 2011:
                    args.mkl_verno = "10"
                elif year >= 2013 and year <= 2016:
                    args.mkl_verno = "11"
    else:
        args.mkl = False
        args.netlib = True

def parse_ddi(args):
    if args.ddi_comm is None and args.mpi_lib:
        # User specified a mpi lib, but no ddi_comm
        # set ddi_comm from default to mpi
        args.ddi_comm = "mpi"
    elif args.ddi_comm is None:
        # User did not specify ddi_comm or mpi_lib, default to sockets
        args.ddi_comm = "sockets"
    if "serial" in args.ddi_comm:
        return
    elif args.ddi_comm == "sockets":
        # Sockets check
        if args.ibm64:
            err("--mpi or --ddi_comm=mpi is required with ibm64 target")
        args.mpi_lib = args.mpi_path = ""
        return
    # ddi_comm is either mpi or mixed
    mpi_lib, mpi_path = args.mpi_lib, args.mpi_path
    if mpi_lib is None:
        args.mpi_lib = auto_detect_intel_impi()[0]
    if mpi_path is None:
        args.mpi_path = auto_detect_intel_impi()[1]
    if not args.mpi_lib and args.ddi_comm == "mpi":
        err("auto detection of mpi lib failed")
    if not args.mpi_path and args.ddi_comm == "mpi" and args.mpi_lib != "CrayCS":
        err("auto detection of mpi_path failed")

def parse_openmp(args):
    if args.hipfort:
        args.openmp=True
        args.openmp_offload=True
    if args.openmp_offload and not args.openmp:
        args.openmp=True
    if args.openmp_offload == True:
       if args.system_target == "summit":
          args.hipblas=False
          args.cublas=True
       elif args.system_target == "crusher" or args.system_target == "frontier":
          args.hipblas=True
          args.cublas=False
       elif args.system_target == "sunspot" or args.system_target == "aurora":
          args.hipblas=False
          args.cublas=False

def parse_fortran_version(args):
    # parse fortran version for certain compilers that is needed by GAMESS
    # build scripts
    if args.fortran == "gfortran":
        major_version = args.fortran_version.split(".")[0]
        minor_version = args.fortran_version.split(".")[1]
        return ("gfortran_verno", '{}.{}'.format(major_version, minor_version))
    # args.fortran == "ifort" or args.fortran == "oneapi-ifort" or
    # args.fortran == "oneapi-ifx":
    elif is_intel(args.fortran):
        try:
            major_version = args.fortran_version.split(".")[0]
        except BaseException:
            major_version = args.fortran_version
        return ("ifort_verno", major_version)
    else:
        return (False, False)


def parse_libcchem(args):
    # LIBCCHEM requires the math library's INCLUDE path
    if args.mathlib_include_path is None:
        if args.math == "mkl":
            mathlib_include_path = auto_detect_mkl()[1]
        else:
            math_path = os.path.abspath(
                os.path.join(args.mathlib_path, os.pardir))
            mathlib_include_path = os.path.join(math_path, "include")
        if os.path.isdir(mathlib_include_path):
            args.mathlib_include_path = mathlib_include_path
        else:
            args.mathlib_include_path = ""
        # else: #TODO add warning
        #    sys.exit("{} : {}".format("Error : Unable to find math \
        #             library INCLUDE path",mathlib_include_path))

# Check conditions of certain arguments after parsing arguments


def condition_dir_exists(condition, path, description):
    if condition:
        check_dir(path, description)


def condition_not_none(obj, description):
    if obj is None:
        err(description)


def condition_get_from_env(arg, env_name):
    if arg is not None:
        return arg
    else:
        # user did not give us this argument, try to read it in from
        # environment
        try:
            arg = os.environ[env_name]
        except KeyError:
            err("{} not defined".format(env_name))
        return arg


def condition_fortran_verno(fortran_version, example):
    condition_not_none(
        fortran_version,
        "--fortran_version is required: e.g., {}".format(example))


def add_arguments(parser):
    # Add arguments to argparse.Parser object, return a dictionary to help
    # handle short cut flags in parsing step
    shortcuts = {}
    parser.add_argument(
        '--path',
        help='path to $GMS_DIR, default: current directory',
        default=os.getcwd(),
        metavar='PATH')
    parser.add_argument(
        '--build_path',
        help='path to build directory, default: current directory',
        default=os.getcwd(),
        metavar='PATH')

    parser.add_argument(
        '--target',
        help='build target, default: linux64',
        default="linux64",
        choices=[
            'linux64',
            'hpe-apollo',
            'hpe-cray-ex',
            'ibm64'])
    # Build target short-cuts
    shortcuts['target'] = [
        'linux64',
        'hpe-apollo',
        'hpe-cray-ex',
        'ibm64']

    parser.add_argument(
        '--system_target',
        help='HPC system target, default: generic',
        default="generic",
        choices=[
            'generic',
            'crusher',
            'frontier',
            'hokulea',
            'nautilus',
            'narwhal',
            'onyx',
            'perlmutter',
            'polaris',
            'summit',
            'sunspot',
            'theta',
            'grace_hopper',
            'qcengine'])
    # HPC system target short-cuts
    shortcuts['system_target'] = [
        'generic',
        'crusher',
        'frontier',
        'hokulea',
        'nautilus',
        'narwhal',
        'onyx',
        'perlmutter',
        'polaris',
        'summit',
        'sunspot',
        'theta',
        'grace_hopper',
        'qcengine']
    parser.add_argument(
        '--skip_process_system_target',
        help='Skip call to process_system_target default: false',
        action="store_true")

    parser.add_argument(
        '--fortran',
        help='Fortran compiler, default: gfortran',
        default="gfortran",
        choices=[
            'acfl',
            'aocc',
            'armflang',
            'gfortran',
            'ifort',
            'oneapi-ifort',
            'oneapi-ifx',
            'oneapi_ifort',
            'oneapi_ifx',
            'nvfortran',
            'pgfortran',
            'xlf'])
    # Fortran compiler short-cuts
    shortcuts['fortran'] = [
        'acfl',
        'aocc',
        'armflang',
        'gfortran',
        'ifort',
        'oneapi-ifort',
        'oneapi-ifx',
        'nvfortran',
        'pgfortran',
        'xlf']

    parser.add_argument(
        '--fortran_version',
        help='Fortran version, default: undefined',
        default=None,
        metavar='VERSION',
        required=False)
    parser.add_argument(
        '--xl_path',
        help='IBM XL path, default: undefined',
        default=None,
        metavar='PATH')
    parser.add_argument(
        '--math',
        help='math library, default: auto-detect MKL or undefined',
        default="netlib",
        choices=[
            'netlib',
            'mkl',
            'openblas',
            'nvblas',
            'essl',
            'aocl',
            'armpl',
            'nvpl'])
    # Math library short-cuts
    shortcuts['math'] = [
        'netlib',
        'mkl',
        'openblas',
        'nvblas',
        'essl',
        'aocl',
        'armpl',
        'nvpl']

    parser.add_argument(
        '--mathlib_path',
        help='math library LIB path, default: auto-detect MKL or undefined',
        default=None,
        metavar='PATH')
    parser.add_argument(
        '--mathlib_include_path',
        help='math library INCLUDE path, default: auto-detect MKL or \
        undefined',
        default=None,
        metavar='PATH')
    parser.add_argument(
        '--libflame_path',
        help='AMD libflame library LIB path, default: use --mathlib_path',
        default=None,
        metavar='PATH')
    parser.add_argument(
        '--mkl_verno',
        help='MKL version, default: auto-detect MKL or undefined',
        default=None)
    parser.add_argument(
        '--threaded-blas',
        help='Link with threaded BLAS',
        action="store_true")

    parser.add_argument(
        '--ddi_comm',
        help='DDI communication choice, default: sockets',
        default=None,
        choices=[
            'sockets',
            'mpi',
            'mixed',
            'serial',
            'serial-debug'])
    # DDI communication short-cuts
    shortcuts['ddi_comm'] = [
        'sockets',
        'mpi',
        'mixed',
        'serial',
        'serial-debug']

    parser.add_argument(
        '--mpi_lib',
        help='MPI library type, default: auto-detected or undefined',
        default=None,
        choices=[
            'impi',
            'hpcx',
            'mpich',
            'mvapich',
            'mvapich2',
            'openmpi',
            'spectrum'])
    # MPI library short-cuts
    shortcuts['mpi_lib'] = [
        'impi',
        'hpcx',
        'mpich',
        'mvapich',
        'mvapich2',
        'openmpi',
        'spectrum']

    parser.add_argument(
        '--mpi_path',
        help='MPI library path, default: auto-detected or undefined',
        default=None,
        metavar='PATH')
    parser.add_argument(
        '--mpi_version',
        help='MPI software version, default: undefined',
        default=None,
        metavar='VERSION',
        required=False)

    parser.add_argument(
        '--libcchem2',
        help='Use LibCChem 2.0, default: false',
        action="store_true"
    )
    parser.add_argument(
        '--libcchem2_path',
        help='path to $GMS_DIR/cchem_2.0, default: current directory/cchem_2.0',
        default=os.getcwd() + '/cchem_2.0',
        metavar='PATH')
    parser.add_argument(
        '--libcchem2_libs',
        help='Additional LIBCCHEM 2.0 libraries default: -lz -lstdc++',
        default='-lz -lstdc++',
        metavar='LIBS')

    parser.add_argument(
        '--libcchem2_use_cuda',
        help='Use cuda, default: false',
        action="store_true"
    )
    parser.add_argument(
        '--libcchem2_use_rocm',
        help='Use rocm, default: false',
        action="store_true"
    )
    parser.add_argument(
        '--cuda_path',
        help='CUDA library path, default: undefined',
        default=None,
        metavar='PATH')
    parser.add_argument(
        '--rocm_path',
        help='ROCM library path, default: undefined',
        default=None,
        metavar='PATH')
    parser.add_argument(
        '--libcchem2_magma_support',
        help='Add GPU support, default: true',
        default=True,
        action="store_true")
    parser.add_argument(
        '--magma_path',
        help='MAGMA library path, default: undefined',
        default=None,
        metavar='PATH')
    parser.add_argument(
        '--libcchem2_use_data_servers',
        help='Use DDI data servers, default: true',
        default=True,
        action="store_true")

    
    parser.add_argument(
        '--libcchem',
        help='Build LIBCCHEM default: false',
        action="store_true")
    parser.add_argument(
        '--libcchem_libs',
        help='Additional LIBCCHEM libraries default: -ldl',
        default='-ldl',
        metavar='LIBS')

    parser.add_argument(
        '--libcchem_gpu_support',
        help='Add GPU support, default: false',
        action="store_true")
    parser.add_argument(
        '--cuda_board',
        help='CUDA board option, default: undefined',
        default=None,
        choices=[
            'tesla',
            'fermi',
            'kepler',
            'pascal',
            'volta'])

    # CUDA board short-cuts
    shortcuts['cuda_board'] = ['tesla', 'fermi', 'kepler', 'pascal', 'volta']

    # LIBCCHEM paths
    parser.add_argument(
        '--hdf5_path',
        help='HDF5 library path, default: undefined',
        default=None,
        metavar='PATH')
    parser.add_argument(
        '--eigen_path',
        help='EIGEN root path, default: undefined',
        default=None,
        metavar='PATH')
    parser.add_argument(
        '--boost_path',
        help='Boost root path, default: undefined',
        default=None,
        metavar='PATH')
    parser.add_argument(
        '--ga_path',
        help='Global Array root path, default: undefined',
        default=None,
        metavar='PATH')

    # LIBCCHEM options
    parser.add_argument(
        '--libcchem_hf_off',
        help='Disable Hartree-Fock (HF)',
        action="store_true")
    parser.add_argument(
        '--libcchem_mp2_off',
        help='Disable second order perturbation theory (MP2)',
        action="store_true")
    parser.add_argument(
        '--libcchem_ri_off',
        help='Disableh resolution-of-the-identify MP2 (RI)',
        action="store_true")
    parser.add_argument(
        '--libcchem_cc_off',
        help='Disable coupled-cluster (CC)',
        action="store_true")
    parser.add_argument(
        '--libcchem_ahfock_off',
        help='Disable accelerated Hartree-Fock (AHFOCK)',
        action="store_true")

    parser.add_argument(
        '--phi',
        help='Intel Xeon Phi option, default: none',
        default="none",
        choices=[
            'knl',
            'knc',
            'none'])
    # Intel Xeon Phi short-cuts
    shortcuts['phi'] = ['knl', 'knc']

    parser.add_argument(
        '--shmtype',
        help='Shared-Memory type, default: sysv',
        default="sysv",
        choices=[
            'sysv',
            'posix'])
    # Shared-memory short-cuts
    shortcuts['shmtype'] = ['sysv', 'posix']

    # Build options
    parser.add_argument(
        '--openmp',
        help='use OpenMP threaded code, default: false',
        action="store_true")

    parser.add_argument(
        '--openmp-offload',
        help='enable OpenMP device offloading, default: false',
        action="store_true")

    parser.add_argument(
        '--hipfort',
        help='use HIPFORT, default: false',
        action="store_true")
    parser.add_argument(
        '--hipfort_path',
        help='HIPFORT library path, default: undefined',
        default=None,
        metavar='PATH')

    parser.add_argument(
        '--cublas',
        help='use CUBLAS during OpenMP Offloading, default: false',
        action="store_true")

    parser.add_argument(
        '--hipblas',
        help='use HIPBLAS during OpenMP Offloading, default: false',
        action="store_true")

    parser.add_argument(
        '--msucc',
        help='build Michigan State University code, default: false',
        action="store_true")

    parser.add_argument(
        '--libxc',
        help='enable library of exchange-correlation functionals \
            for density-functional theory, default: false',
        action="store_true")

    parser.add_argument(
        '--mdi',
        help='enable support for MDI, default: false',
        action="store_true")

    parser.add_argument(
        '--verachem',
        help='build VeraChem VM2 library, default: false',
        action="store_true")
    parser.add_argument(
        '--vm2_path',
        help='full path for the VM2 executable (VC_CompChemPackage_mpi.a)',
        default=None,
        metavar='PATH')

    parser.add_argument(
        '--nbo',
        help='build NBO plug-in, default: false',
        action="store_true")
    parser.add_argument(
        '--nbolib_path',
        help='full path for the NBO library (gmsnbo.i8.a)',
        default=None,
        metavar='PATH')

    parser.add_argument(
        '--neo',
        help='build NEO plug-in, default: false',
        action="store_true")
    parser.add_argument(
        '--tinker',
        help='build TINKER plug-in, default: false',
        action="store_true")
    parser.add_argument(
        '--vb2000',
        help='build VB2000 plug-in, default: false',
        action="store_true")
    parser.add_argument(
        '--xmvb',
        help='build XMVB plug-in, default: false',
        action="store_true")
    parser.add_argument(
        '--rism',
        help='build RISM-SCF-cSED plug-in, default: false',
        action="store_true")

    parser.add_argument(
        '--version',
        help='set GAMESS version name, default: 00',
        default="00",
        metavar='VERSION')

    # Compiler options
    parser.add_argument(
        '--fpe',
        help='build GAMESS binary with floating point exception \
              traps, default: false',
        action="store_true")
    parser.add_argument(
        '--ASAN',
        help='build GAMESS with AddressSanitizer for detecting \
             memory problems (Available only for >GCC 5.0), default: false',
        action="store_true")
    parser.add_argument(
        '--debug_flags',
        help='specify additional flags during compilation, default: none',
        default=None,
        metavar='FLAGS')
    parser.add_argument(
        '--debug_link_flags',
        help='specify additional flags during linking, default: none',
        default=None,
        metavar='FLAGS')

    # Misc. options
    parser.add_argument(
        '--rungms',
        help='render template instantiation of rungms, default: false',
        action="store_true")
    parser.add_argument(
        '--condense',
        help='render condensed template instantiation of rungms, \
             default: false',
        action="store_true")
    parser.add_argument(
        '--scratch',
        help='define scratch directory in rungms',
        default=None,
        metavar='PATH')
    parser.add_argument(
        '--restart',
        help='define restart directory in rungms',
        default=None,
        metavar='PATH')
    parser.add_argument(
        '--coverage',
        help='add GNU or Intel coverage flags for compilation and linking',
        action="store_true")
    # Choose shell output
    parser.add_argument(
        '--shell',
        help='define shell output: tcsh, bash, zsh,csh',
        choices=[
            'tcsh',
            'bash',
            'zsh',
            'csh'],
        default='csh')
    shortcuts['shell'] = ['tcsh', 'bash', 'zsh', 'csh']

    set_shortcuts(shortcuts, handle_shortcut, parser)
    return shortcuts


def parse_arguments(parser, shortcuts):
    args = parser.parse_args()

    # if we used the long_name instead of a shortcut, set that shortcut to
    # 'True' i.e. if we use --shmtype=posix, set args.posix=True
    set_conditional_checks_helper(shortcuts, args)

    # expand flags based on args.system_target
    if not args.skip_process_system_target:
        process_system_target(args)
        set_conditional_checks_helper(shortcuts, args)

    # condense is not a true shortcut, since we're using a different template,
    # but we do want to set rungms to true if it's not already
    if args.condense:
        args.rungms = args.condense

    # some compilers require us to know the Fortran version, so create a list
    # of tuples with the relevant compiler and example version to make it
    # easier to enforce
    fortran_verno_required = [
        (args.ifort,
         "2019.5"),
        (args.oneapi_ifort,
         "2022.1"),
        (args.oneapi_ifx,
         "2022.1"),
        (args.gfortran,
         "4.8.5")]
    for comp, example in fortran_verno_required:
        if comp:
            condition_fortran_verno(args.fortran_version, example)
            break
    # Make sure XLF path is set
    condition_dir_exists(args.xlf, args.xl_path, "IBM XLF compiler")
    if args.mkl:
        # parse and check conditions of mkl option
        parse_mkl_info(args)
    # If we're not using netlib, we need a path to a mathlib
    condition_dir_exists(not args.netlib, args.mathlib_path, "mathlib_path")
    parse_ddi(args)
    parse_openmp(args)
    if args.libcchem:
        parse_libcchem(args)
        if args.libcchem_gpu_support:
            condition_not_none(
                args.cuda_board,
                "CUDA board type is not defined")
            # CUDA root
            args.cuda_path = condition_get_from_env(
                args.cuda_path, 'CUDA_ROOT')
            check_dir(args.cuda_path, "CUDA")
        check_dir(args.hdf5_path, "HDF5")
        check_dir(args.eigen_path, "EIGEN")
        check_dir(args.boost_path, "Boost")
        check_dir(args.ga_path, "Global Arrays")
    if args.libcchem2:
        if args.libcchem2_use_cuda:
            args.cuda_path = condition_get_from_env(
                args.cuda_path, 'CUDA_ROOT'
            )
            check_dir(args.cuda_path, "CUDA")
        elif args.libcchem2_use_rocm:
            args.rocm_path = condition_get_from_env(
                args.rocm_path, 'ROCM_ROOT'
            )
            check_dir(args.rocm_path, "ROCM")
            if args.libcchem2_magma_support:
                args.magma_path = condition_get_from_env(
                    args.magma_path, 'MAGMA_ROOT'
                )
                check_dir(args.magma_path, "MAGMA")
    # VeraChem option requires VM2 executable path
    condition_dir_exists(
        args.verachem,
        args.vm2_path,
        "VM2 executable path variable")
    # NBO requires library path
    condition_dir_exists(args.nbo, args.nbolib_path, "NBO library")
    # coverage check
    parse_coverage(args)
    # fpe option
    if args.fpe:
        args.fpe = auto_detect_fpe(args.fortran)
    else:
        args.fpe = ""
    # ASAN (AddressSanitizer options)
    if args.ASAN:
        parse_ASAN(args)
    return args


def populate_install_info(args, context):
    relevant_fields = [
        "version",
        "path",
        "build_path",
        "target",
        "system_target",
        "skip_process_system_target",
        "fortran",
        "xl_path",
        "ddi_comm",
        "mpi_lib",
        "mpi_path",
        "phi",
        "shmtype",
        "openmp",
        "openmp_offload",
        "hipfort",
        "hipfort_path",
        "cublas",
        "hipblas",
        "msucc",
        "libxc",
        "mdi",
        "shell",
        "threaded_blas"]
    context["install.info"]["cmdline"] = " ".join(sys.argv)
    context["install.info"]["mathlib"] = args.math
    verno, version = parse_fortran_version(args)
    if (verno):
        context["install.info"][verno] = version
    if args.mkl == "mkl":
        relevant_fields.append("mkl_verno")
    if not args.netlib:
        relevant_fields.append("mathlib_path")
    else:
        context["install.info"]["mathlib"] = "none"
    if args.aocl:
        relevant_fields.append("libflame_path")
        process_aocl(args)
    if "mvapich" in args.mpi_lib:
        relevant_fields.append("mpi_version")
    # LIBCCHEM processing
    relevant_fields.append("libcchem")
    if args.libcchem:
        relevant_fields.append("mathlib_include_path")
        relevant_fields.append("libcchem_gpu_support")

        if args.libcchem_gpu_support:
            relevant_fields += ["cuda_board", "cuda_path"]
        relevant_fields += ["hdf5_path", "eigen_path", "boost_path", "ga_path"]

        # Build options
        relevant_fields += ["libcchem_hf_off",
                            "libcchem_mp2_off",
                            "libcchem_ri_off",
                            "libcchem_cc_off",
                            "libcchem_ahfock_off"]
    # LIBCCHEM 2.0 processing
    relevant_fields.append("libcchem2")
    if args.libcchem2:
        relevant_fields.append("libcchem2_use_data_servers")
        relevant_fields.append("libcchem2_path")

        if args.libcchem2_use_cuda:
            relevant_fields += ["cuda_path", "libcchem2_use_cuda"]
        elif args.libcchem2_use_rocm:
            relevant_fields += ["rocm_path", "magma_path", "libcchem2_use_rocm"]
            if args.libcchem2_magma_support:
                relevant_fields += ["libcchem2_magma_support"]

    # OMP_NUM_THREADS
    context["install.info"]["omp_num_threads"] = \
        multiprocessing.cpu_count() - 2
    # OMP_STACKSIZE
    context["install.info"]["omp_stacksize"] = "1G"
    # LIBCCHEM LIBRARIES
    relevant_fields.append("libcchem_libs")
    # VeraChem Option
    relevant_fields.append("verachem")
    if args.verachem:
        relevant_fields.append("vm2_path")
    # NBO Option
    relevant_fields.append("nbo")
    if args.nbo:
        relevant_fields.append("nbolib_path")

    # Plug-in options
    relevant_fields += ["tinker", "vb2000", "xmvb", "neo", "nbo", "rism"]
    # fpe option
    relevant_fields.append("fpe")
    # ASAN (AddressSanitizer options)
    if args.ASAN:
        relevant_fields.append("ASAN")
    if args.coverage:
        relevant_fields.append("coverage")
    # compiler flags
    context["install.info"]["gms_debug_flag"] = str(args.debug_flags)
    context["install.info"]["gms_debug_link_flag"] = str(args.debug_link_flags)
    update_context_args(context, "install.info", args, relevant_fields)


def populate_rungms(args, context):
    relevant_fields = ["system_target", "rungms", "shell", "condense"]
    if args.rungms:
        relevant_fields += ["build_path", "ddi_comm", "scratch", "restart"]
        if args.condense:
            condense = ["libcchem", "cuda_path", 'ddi_comm', 'vb2000']
            relevant_fields += condense
    update_context_args(context, "rungms", args, relevant_fields)


def populate_context(args, context):
    # install.info template instantiation
    populate_install_info(args, context)
    if context["install.info"] == {}:
        err("Unable to build install.info")

    # rungms template instantiation
    populate_rungms(args, context)

    # if args.shell != "csh": #if not default, need to change scripts #NOPE!
    # what if you want to change back?
    files = ["compddi", "compall", "comp", "lked", "Makefile"]
    for f in files:
        context[f]['shell'] = args.shell


def process_system_target(args):
    if args.system_target == "narwhal":
        process_hpcmp_narwhal(args)
    elif args.system_target == "nautilus":
        process_hpcmp_nautilus(args)
    elif args.system_target == "summit":
        process_olcf_summit(args)
    elif args.system_target == "crusher":
        process_olcf_hpe_cray_ex(args)
    elif args.system_target == "frontier":
        process_olcf_hpe_cray_ex(args)
    elif args.system_target == "perlmutter":
        process_nersc_hpe_cray_ex(args)
    elif args.system_target == "polaris":
        process_alcf_hpe_cray_apollo(args)
    elif args.system_target == "sunspot":
        process_alcf_sunspot(args)

def process_hpcmp_narwhal(args):
    args.target="hpe-cray-ex"
    args.fortran="ftn"
    args.mpi_lib="CrayCS"
    args.math="LibSci"
    args.mathlib_path=condition_get_from_env(None,"CRAY_LIBSCI_PREFIX_DIR")


def process_hpcmp_nautilus(args):
    """
    TODO
    """

def process_olcf_summit(args):
    args.target="ibm64"
    if args.gfortran:
        args.fortran_version=condition_get_from_env(None,"OLCF_FAMILY_COMPILER_VERSION")
    try:
       args.xl_path=condition_get_from_env(None,"OLCF_XL_ROOT")
    except internal_error:
       print("Must define environmental variable OLCF_XL_ROOT since XL compiler is not loaded")
    args.mpi_lib="spectrum"
    args.mpi_path=condition_get_from_env(None,"OLCF_SPECTRUM_MPI_ROOT")
    args.math="essl"
    args.mathlib_path=condition_get_from_env(None,"OLCF_ESSL_ROOT")


def process_olcf_hpe_cray_ex(args):
   args.target="hpe-cray-ex"
   args.fortran="ftn"
   args.mpi_lib="CrayCS"
   args.math="LibSci"
   args.mathlib_path=condition_get_from_env(None,"CRAY_LIBSCI_PREFIX_DIR")

def process_nersc_hpe_cray_ex(args):
   args.target="hpe-cray-ex"
   args.fortran="ftn"
   args.mpi_lib="CrayCS"
   args.math="LibSci"
   if "CRAY_PE_LIBSCI_PREFIX_DIR" in os.environ:
       args.mathlib_path=condition_get_from_env(None,"CRAY_PE_LIBSCI_PREFIX_DIR")
   elif "CRAY_LIBSCI_PREFIX_DIR" in os.environ:
       args.mathlib_path=condition_get_from_env(None,"CRAY_LIBSCI_PREFIX_DIR")

def process_alcf_hpe_cray_apollo(args):
   args.target="hpe-apollo"
   args.fortran="ftn"
   args.mpi_lib="CrayCS"
   args.math="LibSci"
   if "CRAY_PE_LIBSCI_PREFIX_DIR" in os.environ:
       args.mathlib_path=condition_get_from_env(None,"CRAY_PE_LIBSCI_PREFIX_DIR")
   elif "CRAY_LIBSCI_PREFIX_DIR" in os.environ:
       args.mathlib_path=condition_get_from_env(None,"CRAY_LIBSCI_PREFIX_DIR")

def process_alcf_sunspot(args):
   args.target="linux64"
   if args.fortran == "gfortran":
       args.fortran="oneapi-ifx"
       args.fortran_version=condition_get_from_env(None,"AURORA_BASE_ENV")
   args.mpi_lib="mpich"
   args.mpi_path=condition_get_from_env(None,"MPI_ROOT")
   args.math="mkl"
   mkl_path=condition_get_from_env(None,"MKLROOT")
   args.mathlib_path=os.path.join(mkl_path, "lib/intel64")


def process_aocl(args):
   if args.libflame_path is None:
      args.libflame_path=args.mathlib_path


def create(context):

    parser = argparse.ArgumentParser(
        description='GAMESS Build Template Instantiation')
    shortcuts = add_arguments(parser)
    # pass a dict of short cuts and long names, in order to set arguments so
    # that we can do simplier conditional checks
    args = parse_arguments(parser, shortcuts)
    populate_context(args, context)


def main():
    context = {
        "install.info": {},
        "rungms": {},
        "compddi": {},
        "compall": {},
        "comp": {},
        "lked": {},
        "Makefile": {}}
    create(context)
    gms_path = context["install.info"]["path"]

    TEMPLATE_ENVIRONMENT = Environment(
        autoescape=False,
        loader=FileSystemLoader(os.path.join(gms_path, 'misc/automation')),
        trim_blocks=False)

    write_files = ["install.info"]

    condense = False
    if context["rungms"]["rungms"] == 'true':
        write_files.append("rungms")
        condense = context["rungms"]["condense"] == "true"
    for filename in write_files:
        template_name = filename + ".template"
        if condense and filename == "rungms":
            print("Condensing rungms")
            template_name = "condensed_" + template_name
        with open(filename, "w") as f:
            content = TEMPLATE_ENVIRONMENT.get_template(
                template_name).render(context[filename])
            f.write(content + "\n")
    shell_files = [
        "compddi",
        "compall",
        "comp",
        "lked",
        "Makefile",
        "rungms"]
    for filename in shell_files:
        # for now, it makes more sense to modify existing files than to create
        # templates for one parameter.
        old = "bin/csh"  # default
        new = "bin/" + context[filename]["shell"]
        if filename == "compddi":
            filename = os.path.join(gms_path, "ddi/" + filename)
        elif filename == "rungms":  # for runtest.py
            filename = os.path.join(gms_path, "misc/automation/rungms")
        assert os.path.exists(filename)
        content = ""
        with open(filename, "r") as f:
            content = f.read()  # assume we can fit entire file in memory
            old = re.search(r"bin\/(z|c|ba|tc)sh", content).group(0)
            content = content.replace(old, new)
        with open(filename, "w") as f:
            f.write(content)


if __name__ == "__main__":
    main()
