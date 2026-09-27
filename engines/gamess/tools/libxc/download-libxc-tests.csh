#!/usr/bin/env csh

if (-e install.info) then
   source install.info
else if (-e ../install.info) then
   source ../install.info
else
   echo "Please run 'config' first, to set up GAMESS compiling information"
   exit 4
endif

if ($GMS_TARGET == win64) then
  set SSL_REVOKE="--ssl-no-revoke"
else
  set SSL_REVOKE=""
endif

set archive=libxc-tests.zip

if ( $#argv == 1 ) then
  if ( -f $argv[1]/$archive ) then
    test $argv[1]/$archive -ef ${GMS_PATH}/$archive
    if ( ! $status ) then
      echo "$argv[1]/$archive and ${GMS_PATH}/$archive files are same."
      echo "They will not be copied!"
    else
      echo "Copying $argv[1]/$archive to ${GMS_PATH}/$archive..."
      cp $argv[1]/$archive ${GMS_PATH}/$archive
    endif
  else
    echo "$argv[1]/$archive not found!"
    exit 1
  endif
else if ( $#argv == 0 ) then
  if(! -e ${GMS_PATH}/$archive) then
    curl ${SSL_REVOKE} -Lk https://github.com/gms-bbg/libxc-tests/archive/master.zip --output ${GMS_PATH}/$archive
  endif
else
  echo "Too many arguments..."
  exit 1
endif

unzip ${GMS_PATH}/$archive -d ${GMS_PATH}/tests
rm -rf ${GMS_PATH}/tests/dft-libxc
mv ${GMS_PATH}/tests/libxc-tests-master/dft-libxc ${GMS_PATH}/tests/dft-libxc
rm -rf ${GMS_PATH}/tests/libxc-tests-master
