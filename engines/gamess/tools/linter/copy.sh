#!/bin/sh

SRCDIR=tools/linter/sources
MAP=tools/linter/map.txt

pwd
mkdir -p $SRCDIR

for i in vb2000/SRC/*.src
do
  name=`basename $i`
  name=${name%.*}
  case $name in
    mod_vb2000 )
      continue
      ;;
    *)
      name=$name.f
      ;;
  esac
  cp $i $SRCDIR/$name
  echo $name $i >> $MAP
done

for i in qmnuc/neo/*.src
do
  name=`basename $i`
  name=${name%.*}.f
  cp $i $SRCDIR/$name
  echo $name $i >> $MAP
done

for i in misc/*.src
do
  name=`basename $i`
  name=${name%.*}.f
  cp $i $SRCDIR/$name
  echo $name $i >> $MAP
done

for i in source/*.src
do
  name=`basename $i`
  name=${name%.*}
  case $name in
    ccsd3aacgreorder | ccsd3aacgsum | cublasf | functionals | hpcchem_driver | hpcchem_module | hpcchem_molecule_module | hpcchem_runner_module | hpcchem_stubs | libxc | modmdi )
      continue
      ;;
    blkint | constants | cphf_tddnlr_omp | dftbpb | efpmodule | efteiomp | gamess_version | gausshermite | gpu_* | grd2_consts | hrmrst | i00 | i01 | i02 | i11 | i12 | i22 | libxc | libxc_empty | messages | mod_1e_primitives | mod_dft_fuzzycell | mod_dft_gridint | mod_dft_molgrid | mod_dft_partfunc | mod_gauss_hermite | mod_grid_storage | mod_lutiotc | mod_nameio | mod_nosp_basis | mod_offload | mod_sformas | mod_shell_tools | mod_vvos | modcc | modmcpdft | modmdi | modmdi_empty | modmnfun | modules_common | modules_dft | modules_rimp2gpu | mpcdatpm6 | mpchbond | mx_limits | norm2 | omp* | params | prec | riccdiag | riccints | riccrhf | riccutils | rimp2gpu | rimp2grd | rimp2int | rimp2omp | rimp2subs | rt1 | rt2 | rt3 | secor | shellpairs | shellprocessing | tgpu_* | utils_strings | vxx | work )
      name=$name.F90
      ;;
    ccsdt | gamess | grd2a | int2a | rimp2 | rimp2omp | ryspol | unport | vb2000 )
      name=$name.F
      ;;
    *)
      name=$name.f
      ;;
  esac
  cp $i $SRCDIR/$name
  echo $name $i >> $MAP
done

if [ $(uname) = "Darwin" ]; then
  sed -i '' "s/^*I64/CI64/" $SRCDIR/unport.f
  sed -i '' "s/^*L64/    /" $SRCDIR/unport.f
  sed -i '' "s/^*I64/    /" $SRCDIR/*
else
  sed -i "s/^*I64/CI64/" $SRCDIR/unport.f
  sed -i "s/^*L64/    /" $SRCDIR/unport.f
  sed -i "s/^*I64/    /" $SRCDIR/*
fi
