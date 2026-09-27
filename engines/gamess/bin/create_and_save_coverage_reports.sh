#!/bin/bash
pushd object >/dev/null

find . -maxdepth 1 -type f -writable \( -name '*.f' -o -name '*.F90' -o -name '*.F' \) -print0 | xargs -0 -r gcov -m -c 1>gcov_output.txt 2>gcov_errors.txt

num=1
while [[ -d $(printf "%04d" $num)  && $num -lt 10000 ]]; do
    ((num++))
done

nam=$(printf "%04d" $num)
mkdir $nam

find . -maxdepth 1 -type f -writable -name '*.gcov' -print0 | xargs -0 -r mv -t $nam
mv -t $nam gcov_output.txt gcov_errors.txt

echo "$1" >$nam/coverage.info

popd >/dev/null
