#!/bin/bash -le

setup_environment () {
    export PYTHONUNBUFFERED=1
    export GCCVERSION=5.5.0
    export GNU_FORTRAN_VERSION=5.5.0
    export PATH="/shared/compilers/gcc/${GCCVERSION}/bin:${PATH}"
    PATH="/shared/make/bin:${PATH}"
    PATH="/shared/cmake/3.13.4/bin:${PATH}"
    export LD_LIBRARY_PATH="/shared/compilers/gcc/${GCCVERSION}/lib64:${LD_LIBRARY_PATH}"
    LD_LIBRARY_PATH="/shared/make/lib:${LD_LIBRARY_PATH}"
    export LD_RUN_PATH="/shared/compilers/gcc/${GCCVERSION}/lib64:${LD_RUN_PATH}"
    export CPATH="/shared/compilers/gcc/${GCCVERSION}/include:${CPATH}"
    export INCLUDEPATH="/shared/compilers/gcc/${GCCVERSION}/include:${INCLUDEPATH}"
    export INTEL_LICENSE_FILE=/opt/intel/licenses
    INTEL_LICENSE_FILE+=:/export/home/jenkins/intel/licenses
    INTEL_LICENSE_FILE+=:/shared/compilers/intel/oneapi/clck/2021.6.0/licensing
    INTEL_LICENSE_FILE+=:/opt/intel/licenses
    INTEL_LICENSE_FILE+=:/export/home/jenkins/intel/licenses
    INTEL_LICENSE_FILE+=':/Users/Shared/Library/Application Support/Intel/Licenses'
    export INTEL_COMPILER_VERSION=2022.1
    export INTEL_COMPILER_SOURCE_PATH=/shared/compilers/intel/oneapi/setvars.sh
    export INTEL_MKL_SOURCE_PATH=/shared/compilers/intel/oneapi/mkl/latest/env/vars.sh
    export INTEL_MPI_SOURCE_PATH=/shared/compilers/intel/oneapi/mpi/latest/env/vars.sh
    export ATLAS_BASE_PATH="/shared/math/atlas/3.10.3-gnu-${GCCVERSION}"
    export OPENBLAS_GCCVERSION=5.5.0
    export OPENBLAS_VERSION=0.3.15
    export OPENBLAS_BASE_PATH="/shared/math/openblas/${OPENBLAS_VERSION}"
    export OPENMPI_VERSION=1.10.7
    export OPENMPI_BASE_PATH="/shared/mpi/openmpi-${OPENMPI_VERSION}-gnu-${GCCVERSION}"
    export NUM_CPU_CORES=6
    export MAX_GMS_CORES=3
}

setup_environment_omp_gnu_impi_openblas () {
    export OMP_NUM_THREADS=3
}

create_install_info () {
    local stage=$1
    local fortran=$2
    local comm=$3
    local mathlib=$4
    local openblas_path_threaded=$OPENBLAS_BASE_PATH
    local openblas_path_serial=$OPENBLAS_BASE_PATH
    openblas_path_threaded+="/sandybridge-gnu-${OPENBLAS_GCCVERSION}-threaded/lib"
    openblas_path_serial+="/sandybridge-gnu-${OPENBLAS_GCCVERSION}-serial/lib"
    local cii_args
    declare -a cii_args=()
    if [[ $fortran == gnu ]]; then
        cii_args+=(--fortran_version="${GNU_FORTRAN_VERSION}")
        if [[ $stage == Pre-check ]]; then
            cii_args+=(--fpe)
        elif [[ $stage == GAMESS ]]; then
            cii_args+=(--fpe)
        elif [[ $stage == GAMESS_OMP ]]; then
            cii_args+=(--fpe)
	fi
    else
        echo >&2 'Error: Unexpected fortran argument.'
        exit 1
    fi
    if [[ $mathlib == mkl ]]; then
        cii_args+=(--mkl)
    elif [[ $mathlib == openblas ]]; then
        cii_args+=(--openblas)
        if [[ $stage == GAMESS_OMP ]]; then
            cii_args+=(--mathlib_path="${openblas_path_threaded}")
        else
            cii_args+=(--mathlib_path="${openblas_path_serial}")
        fi
    else
        echo >&2 'Error: Unexpected mathlib argument.'
        exit 1
    fi
    if [[ $comm == sockets ]]; then
        cii_args+=(--sockets)
    elif [[ $comm == impi ]]; then
        cii_args+=(--impi)
    else
        echo >&2 'Error: Unexpected comm argument.'
        exit 1
    fi
    if [[ $stage == Pre-check ]]; then
        cii_args+=(--rungms)
    elif [[ $stage == GAMESS ]]; then
        cii_args+=(--rungms)
    elif [[ $stage == GAMESS_OMP ]]; then
        cii_args+=(--openmp --rungms)
    elif [[ $stage == NEO ]]; then
        cii_args+=(--rungms --neo)
    elif [[ $stage == VB2000 ]]; then
        cii_args+=(--rungms --vb2000)
    else
        echo >&2 'Error: Unexpected stage argument.'
        exit 1
    fi
    set -x
    if [[ $mathlib == mkl ]]; then
        # shellcheck source=/dev/null
        source "${INTEL_MKL_SOURCE_PATH}" intel64 > /dev/null 2>&1
    fi
    if [[ $comm == impi ]]; then
        # shellcheck source=/dev/null
        source "${INTEL_MPI_SOURCE_PATH}" > /dev/null 2>&1
    fi
    bin/create-install-info.py "${cii_args[@]}"
    cat install.info
    { set +x; } 2>/dev/null
}

build () {
    set -x
    if [[ $1 == clean_ddi ]]; then
        make clean_ddi
    fi
    if [[ $2 == rm_gamess ]]; then
        rm gamess.00.x
    fi
    if [[ $3 == rm_object ]]; then
        find object -mindepth 1 -delete
    fi
    if [[ $4 == source_mkl ]]; then
        # shellcheck source=/dev/null
        source "${INTEL_MKL_SOURCE_PATH}" intel64 > /dev/null 2>&1
    fi
    if [[ $5 == make_parallel ]]; then
        make -j "$(nproc)" -Otarget
    else
        make
    fi
    { set +x; } 2>/dev/null
}

run_openmp_tests () {
    local test_type=$1
    local file_path=efp-mpiomp/makeefp/local_E-R
    file_path+=,efp-mpiomp/makeefp/makefp_scr02
    file_path+=,efp-mpiomp/makeefp/pol_poldyn
    file_path+=,efp-mpiomp/makeefp/stone-dma
    file_path+=,rimp2grd-mpiomp/Br2
    file_path+=,rimp2grd-mpiomp/I2.inp
    file_path+=,rimp2grd-mpiomp/I2NCH
    file_path+=,rimp2grd-mpiomp/benz.inp
    file_path+=,rimp2grd-mpiomp/benz_optimize
    file_path+=,ricc-mpiomp/benzene
    file_path+=,comp-mpiomp/rig3_exam43
    if [[ $test_type == serial ]]; then
        set -x
        find scratch restart -mindepth 1 -delete
        bin/my_ipcrm
        tests/runtest.py --bwrap --folder=rhf-mpiomp \
                         --skip_folder=hpc --ncpus=1 -t 2
        tests/runtest.py --bwrap --filepath="${file_path}" --ncpus=1 -t 2
        { set +x; } 2>/dev/null
    elif [[ $test_type == parallel ]]; then
        set -x
        find scratch restart -mindepth 1 -delete
        bin/my_ipcrm
        tests/runtest.py --bwrap --folder=rhf-mpiomp/parallel \
                         --skip_folder=hpc --ncpus=2
        tests/runtest.py --bwrap --filepath="${file_path}" --ncpus=2
        { set +x; } 2>/dev/null
    else
        echo >&2 'Error: Unexpected test_type.'
        exit 1
    fi
}

run_travisci_tests () {
    local test_type=$1
    if [[ $test_type == serial ]]; then
        set -x
        find scratch restart -mindepth 1 -delete
        bin/my_ipcrm
        tests/runtest.py --bwrap --folder=travis-ci -t 6
        { set +x; } 2>/dev/null
    elif [[ $test_type == parallel ]]; then
        set -x
        find scratch restart -mindepth 1 -delete
        bin/my_ipcrm
        tests/runtest.py --bwrap --folder=travis-ci/parallel -n 3 -t 2
        { set +x; } 2>/dev/null
    else
        echo >&2 'Error: Unexpected test_type.'
        exit 1
    fi
}

run_additional_tests () {
    local folder=$1
    set -x
    find scratch restart -mindepth 1 -delete
    bin/my_ipcrm
    tests/runtest.py --bwrap --folder="${folder}" -t 6
    { set +x; } 2>/dev/null
}

validate_tests () {
    local test_path=$1
    local checkgms_cmd='tests/checkgms.py -p -v --skip_json_create'
    checkgms_cmd+=" --pass_delete --test_path=${test_path}"
    # shellcheck disable=SC2016
    checkgms_cmd+=' --filepath="$0.log" >"$0.checkgms"'
    local sed_script='1{h;s/ .*$//;s/^/ /;x};G;s/\n//;p'
    # shellcheck disable=SC2016
    local sed_cmd='sed "s:1 / 1:$0 / $2:" $1'
    set -x
    find "${test_path}" -name '*.log' | sort | sed 's/.log$//' | \
        xargs -r -n 1 -P 6 sh -c "${checkgms_cmd}"
    find "${test_path}" -name '*.checkgms' | sort | \
        nl -s' ' -n'ln' -w1 | tac | sed -ne "${sed_script}" | \
        tac | xargs -n 3 sh -c "${sed_cmd}"
    find "${test_path}" -name '*.checkgms' -delete
    tests/checkgms.py -p --skip_json_create --test_path="${test_path}" \
                      --exit_on_fail
    find tests -name '*.log' -print -delete
    { set +x; } 2>/dev/null
}

pipeline () {
    setup_environment
    if [[ $1 == 'Environment; Setup clean environment' ]]; then
        set -x
        hostname
        module purge || true
        env
        { set +x; } 2>/dev/null
    elif [[ $1 == 'Pre-check; Create install.info' ]]; then
        create_install_info Pre-check gnu sockets mkl
    elif [[ $1 == 'Pre-check; Checking for merge conflicts' ]]; then
        set -x
        if [[ -n $(git ls-files -u) ]]; then
            echo 'There is a merge conflict that has not been resolved. Aborting'
            exit 1
        fi
        { set +x; } 2>/dev/null
    elif [[ $1 == 'Pre-check; Build FTNCHEK' ]]; then
        set -x
        bin/travis-ci-build-ftnchek.sh
        { set +x; } 2>/dev/null
    elif [[ $1 == 'Pre-check; Run FTNCHEK' ]]; then
        set -x
        bin/travis-ci-run-ftnchek.sh
        { set +x; } 2>/dev/null
    elif [[ $1 == 'Pre-check; Download LAPACK' ]]; then
        set -x
        tools/lapack/download-lapack.csh
        { set +x; } 2>/dev/null
    elif [[ $1 == 'GAMESS gnu-sockets-mkl; Create install.info' ]]; then
        create_install_info GAMESS gnu sockets mkl
    elif [[ $1 == 'GAMESS gnu-sockets-mkl; Build' ]]; then
        set -x
        set -o pipefail
        rm -rf object
        mkdir object
        make clean_exams
        { set +x; } 2>/dev/null
        build clean_ddi no_rm_gamess no_rm_object source_mkl make_parallel
    elif [[ $1 == 'GAMESS gnu-sockets-mkl; Test serial' ]]; then
        set -x
        rm -rf scratch restart
        mkdir scratch restart
        { set +x; } 2>/dev/null
        run_travisci_tests serial
    elif [[ $1 == 'GAMESS gnu-sockets-mkl; Validate serial' ]]; then
        validate_tests tests/travis-ci
    elif [[ $1 == 'GAMESS gnu-sockets-mkl; Test parallel' ]]; then
        run_travisci_tests parallel
    elif [[ $1 == 'GAMESS gnu-sockets-mkl; Validate parallel' ]]; then
        validate_tests tests/travis-ci/parallel
    elif [[ $1 == 'GAMESS gnu-impi-mkl; Create install.info' ]]; then
        create_install_info GAMESS gnu impi mkl
    elif [[ $1 == 'GAMESS gnu-impi-mkl; Build' ]]; then
        build clean_ddi rm_gamess no_rm_object source_mkl make_parallel
    elif [[ $1 == 'GAMESS gnu-impi-mkl; Test serial' ]]; then
        run_travisci_tests serial
    elif [[ $1 == 'GAMESS gnu-impi-mkl; Validate serial' ]]; then
        validate_tests tests/travis-ci
    elif [[ $1 == 'GAMESS gnu-impi-mkl; Test parallel' ]]; then
        run_travisci_tests parallel
    elif [[ $1 == 'GAMESS gnu-impi-mkl; Validate parallel' ]]; then
        validate_tests tests/travis-ci parallel
    elif [[ $1 == 'GAMESS gnu-impi-openblas; Create install.info' ]]; then
        create_install_info GAMESS gnu impi openblas
    elif [[ $1 == 'GAMESS gnu-impi-openblas; Build' ]]; then
        build no_clean_ddi rm_gamess no_rm_object no_source_mkl make_parallel
    elif [[ $1 == 'GAMESS gnu-impi-openblas; Test serial' ]]; then
        run_travisci_tests serial
    elif [[ $1 == 'GAMESS gnu-impi-openblas; Validate serial' ]]; then
        validate_tests tests/travis-ci
    elif [[ $1 == 'GAMESS gnu-impi-openblas; Test parallel' ]]; then
        run_travisci_tests parallel
    elif [[ $1 == 'GAMESS gnu-impi-openblas; Validate parallel' ]]; then
        validate_tests tests/travis-ci parallel
    elif [[ $1 == 'GAMESS OMP gnu-impi-openblas; Create install.info' ]]; then
        setup_environment_omp_gnu_impi_openblas
        create_install_info GAMESS_OMP gnu impi openblas
    elif [[ $1 == 'GAMESS OMP gnu-impi-openblas; Build' ]]; then
        setup_environment_omp_gnu_impi_openblas
        build no_clean_ddi rm_gamess rm_object no_source_mkl make_parallel
    elif [[ $1 == 'GAMESS OMP gnu-impi-openblas; Test serial' ]]; then
        setup_environment_omp_gnu_impi_openblas
        run_openmp_tests serial
    elif [[ $1 == 'GAMESS OMP gnu-impi-openblas; Validate serial' ]]; then
        setup_environment_omp_gnu_impi_openblas
        validate_tests tests
    elif [[ $1 == 'GAMESS OMP gnu-impi-openblas; Test parallel' ]]; then
        setup_environment_omp_gnu_impi_openblas
        run_openmp_tests parallel
    elif [[ $1 == 'GAMESS OMP gnu-impi-openblas; Validate parallel' ]]; then
        setup_environment_omp_gnu_impi_openblas
        validate_tests tests
    elif [[ $1 == 'NEO gnu-impi-mkl; Create install.info' ]]; then
        create_install_info NEO gnu impi mkl
    elif [[ $1 == 'NEO gnu-impi-mkl; Build' ]]; then
        build no_clean_ddi rm_gamess rm_object source_mkl make_parallel
    elif [[ $1 == 'NEO gnu-impi-mkl; Test serial' ]]; then
        run_additional_tests neo
    elif [[ $1 == 'NEO gnu-impi-mkl; Validate serial' ]]; then
        validate_tests tests/neo
    elif [[ $1 == 'VB2000 gnu-impi-mkl; Create install.info' ]]; then
        create_install_info VB2000 gnu impi mkl
    elif [[ $1 == 'VB2000 gnu-impi-mkl; Build' ]]; then
        build no_clean_ddi rm_gamess no_rm_object source_mkl no_make_parallel
    elif [[ $1 == 'VB2000 gnu-impi-mkl; Test serial' ]]; then
        run_additional_tests vb2000
    elif [[ $1 == 'VB2000 gnu-impi-mkl; Validate serial' ]]; then
        validate_tests tests/vb2000
    else
        echo >&2 'Error: Unexpected combination of command-line arguments.'
        exit 1
    fi
}

if (( $# == 2 )); then
    pipeline "${1}; ${2}"
else
    echo >&2 'Error: Exactly two command-line arguments are needed.'
    exit 1
fi
