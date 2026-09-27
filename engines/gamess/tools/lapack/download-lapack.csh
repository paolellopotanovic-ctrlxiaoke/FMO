#!/bin/csh -f

if (-e install.info) then
   source install.info
else if (-e ../install.info) then
   source ../install.info
else
   echo "Please run 'config' first, to set up GAMESS compiling information"
   exit 4
endif

if (${GMS_TARGET} == win64) then
  set SSL_REVOKE="--ssl-no-revoke"
else
  set SSL_REVOKE=""
endif

set folder=lapack

set MAJOR=3
set MINOR=12
set PATCH=0

mkdir -p ${GMS_PATH}/3rd-party/${folder}

set archive=lapack-${MAJOR}.${MINOR}.${PATCH}.tar.gz

if ( ${#argv} == 1 ) then
  if ( -f ${argv}[1]/${archive} ) then
    test ${argv}[1]/${archive} -ef ${GMS_PATH}/${archive}
    if ( ! ${status} ) then
      echo "${argv}[1]/${archive} and ${GMS_PATH}/${archive} files are same."
      echo "They will not be copied!"
    else
      echo "Copying ${argv}[1]/${archive} to ${GMS_PATH}/${archive}..."
      cp ${argv}[1]/${archive} ${GMS_PATH}/${archive}
    endif
  else
    echo "${argv}[1]/${archive} not found!"
    exit 1
  endif
else if ( ${#argv} == 0 ) then
  if(! -e ${GMS_PATH}/${archive}) then
    curl -LJO https://github.com/Reference-LAPACK/lapack/archive/refs/tags/v${MAJOR}.${MINOR}.${PATCH}.tar.gz ${SSL_REVOKE} -o ${GMS_PATH}/$archive
  endif
else
  echo "Too many arguments..."
  exit 1
endif

tar -xzf ${GMS_PATH}/${archive} -C ${GMS_PATH}/3rd-party/${folder} --strip-components=1

patch -d ${GMS_PATH} -p1 < ${GMS_PATH}/tools/lapack/CheckLAPACKCompilerFlags.patch

