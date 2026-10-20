#!/bin/bash -e
# Written by Mamoru Haruta

input="${1}"

#-- Setup
tmp=".tmp.$$"
infolines=4 #How many lines do the information have
grepword="COORDINATES OF EACH IMAGES \[angstrom\]"
grep "${grepword}" -A ${infolines} -n ${input} > ${tmp}

#-- Get information
nimages=$(grep "Number of images" -m 1 ${tmp} | grep -o "[0-9]*$")
natoms=$(grep "Number of atoms" -m 1 ${tmp} | grep -o "[0-9]*$")
echo "nimages=${nimages}"
echo "natoms=${natoms}"

#-- Get first cycle
firstcycle=$(grep "Current cycle" ${tmp} | grep -o "[0-9]*$" | head -n 1)
echo "firstcycle=${firstcycle}"

#-- Get last cycle
lastcycle=$(grep "Current cycle" ${tmp} | grep -o "[0-9]*$" | tail -n 1)
if [ "${lastcycle}" = "" ]; then
    echo "ERROR: Something wrong with your file."
    echo "Aborted."
    exit 1
fi
echo "lastcycle=${lastcycle}"
echo "This job was finished in ${lastcycle} cycles."

rm -f ${tmp}
exit
