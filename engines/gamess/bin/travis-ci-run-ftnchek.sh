#!/bin/bash -x
#
export PATH=$HOME/bin/:$PATH
which ftnchek
cd source && USERSCR=../scratch/ ../tools/checkgms all 

# hold error codes
RESULT=0

grep -A 2 "Common block" ../scratch/gamess.lis | grep -A 2 "varying" > common_varying
if [ -s common_varying ]
then
    echo "There are common blocks with varying length: "
    cat common_varying
    RESULT=$(($RESULT+1))
else
    echo "No common blocks with varying lengths."
fi

# check for data type mismatch
grep -A 2 "Common block" ../scratch/gamess.lis | grep -A 2 "type mismatch" > common_mismatch
if [ -s common_mismatch ]
then
    echo "There are common blocks with data types that don't match: "
    cat common_mismatch
    RESULT=$(($RESULT+1))
else
    echo "No common blocks with mismatched types."
fi


# check for > 72 lines
grep -A 1 "characters" ../scratch/gamess.lis | grep -B 1 "columns" > characters
if [ -s characters ]
then
    echo "There are code lines where the length goes past 72 columns"
    cat characters
    RESULT=$(($RESULT+1))
else
    echo "All lines <= 72 columns."
fi

# returns 0 if no errors, returns non=zero (error) if there are errors
exit $RESULT
