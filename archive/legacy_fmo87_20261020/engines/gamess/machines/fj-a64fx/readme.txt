Files included in this directory are adapted for fj-a64fx platform (Fujitsu
FX1000 and FX700 ARM64 machines).

Please ask your admins if unsure what your platform is (FX1000 or FX700) and 
update value of PLATFORM in scripts './rungms' and './rungms.auto' if necessary.

Presently, if your Fujitsu HPC cluster runs jobs under PJM queueing system, your 
platform is FX1000. If the cluster runs PBS queueing system, it ought to be FX700.

After that, files in this drectory are to replace the following scripts 
(backup originals if necessary): 

cp rungms ../../
cp runall ../..
cp rungms.auto ../../tools/localgms
cp rungms.auto ../../misc/automation/rungms

Also in case of building GAMESS with VB2000 do

cp runallvb ../../vb2000
