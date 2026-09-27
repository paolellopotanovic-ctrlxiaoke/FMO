#!/bin/bash -x

#first try pulling from netlib. if that fails, try www.dsm.fordham.edu
wget --no-check-certificate netlib.org/fortran/ftnchek.tgz
output=$?
# Abort on Error
set -e
if [ $output == 0 ]
then
    gunzip ftnchek.tgz
    tar -xf ftnchek.tar
else
    wget --no-check-certificate "https://www.dsm.fordham.edu/~ftnchek/download/ftnchek-3.3.1.tar.gz"
    gunzip ftnchek-3.3.1.tar.gz
    tar -xf ftnchek-3.3.1.tar
fi
cd ftnchek-3.3.1/
# this is needed to redimension a table in FTNCHEK,
# changing HASHSZ from 20930 to 80930 in ftnchek.h, since GAMESS
# is currently too big to fit in the default space.
# if ftnchek starts hanging in symtab.c, this might need to be increased more 
sed -i s/"#define HASHSZ 20930"/"#define HASHSZ 80930"/ ftnchek.h
./configure --prefix=$HOME
make
make check
make install
