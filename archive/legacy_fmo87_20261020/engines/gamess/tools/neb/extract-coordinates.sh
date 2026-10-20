#!/bin/bash -e
# Written by Mamoru Haruta

input="${1}"
cycle="${2}"
output="${1}.${2}"

#-- Usage out
if [ $# -ne 2 ]; then
    echo "-- Extract coordinate file from output --"
    echo ""
    echo "usage: `basename $0` [input] [cycle]"
    echo ""
    echo "input = Your GAMESS output log or trj file name"
    echo "cycle = The total cycle number you want to extract (number or 'last')"
    exit 1
fi

#-- Check
if echo "${input}" | grep -i ".dat" ; then
    echo "Input log or trj file, not dat file."
    exit 1
fi

#-- File type
xyzsingle=1
xyzoverlap=1
cc1single=1
cc1overlap=1

echo "Running..."

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

#-- Set cycle number to last cycle
if [ "${cycle}" = "last" ]; then
    cycle="${lastcycle}"
fi

#-- Check cycle
if [ "${cycle}" -lt "${firstcycle}" ]; then
    echo "ERROR: Your specified cycle ${cycle} is less than ${firstcycle}."
    echo "Aborted."
    exit 1
elif [ "${cycle}" -gt "${lastcycle}" ]; then
    echo "ERROR: Your specified cycle ${cycle} is grater than ${lastcycle}."
    echo "Aborted."
    exit 1
fi

#-- Get cycle number
line=$(grep "Current cycle.* ${cycle}$" ${tmp} | grep -o "^[0-9]\+")
echo "line=${line}"
line=$((${line} + 4))
echo "line=${line}"

#-- Extract
rm -f ${output}.xyz ${output}.ol.xyz ${output}.cc1 ${output}.ol.cc1

if [ ${xyzoverlap} -eq 1 ]; then
    echo "$((${natoms} * ${nimages}))" >> ${output}.ol.xyz
    echo "overlap" >> ${output}.ol.xyz
fi
if [ ${cc1overlap} -eq 1 ]; then
    echo "$((${natoms} * ${nimages})) #image = ${iimage}" >> ${output}.ol.cc1
fi
continuous=0

for iimage in $(seq 1 ${nimages})
do
    if [ ${xyzsingle} -eq 1 ]; then
        echo "${natoms}" >> ${output}.xyz
        echo "#image = ${iimage}" >> ${output}.xyz
    fi
    if [ ${cc1single} -eq 1 ]; then
        echo "${natoms} #image = ${iimage}" >> ${output}.cc1
    fi
    startline=$((${line} + 1))
    endline=$((${line} + ${natoms}))
    if [ ${xyzsingle} -eq 1 ]; then
        #sed -n ${startline},${endline}p ${input} | awk '{printf "%-8s %18.10f%18.10f%18.10f\n",$1,$3,$4,$5}' >> ${output}.xyz
        sed -n ${startline},${endline}p ${input} | awk '{printf "%-5.1f %18.10f%18.10f%18.10f\n",$2,$3,$4,$5}' >> ${output}.xyz
    fi
    if [ ${xyzoverlap} -eq 1 ]; then
        #sed -n ${startline},${endline}p ${input} | awk '{printf "%-8s %18.10f%18.10f%18.10f\n",$1,$3,$4,$5}' >> ${output}.ol.xyz
        sed -n ${startline},${endline}p ${input} | awk '{printf "%-5.1f %18.10f%18.10f%18.10f\n",$2,$3,$4,$5}' >> ${output}.ol.xyz
    fi
    if [ ${cc1single} -eq 1 ]; then
        sed -n ${startline},${endline}p ${input} | awk '{printf "%-8s %4d %18.10f%18.10f%18.10f\n",$1,NR,$3,$4,$5}' >> ${output}.cc1
    fi
    if [ ${cc1overlap} -eq 1 ]; then
        sed -n ${startline},${endline}p ${input} | awk -v n=${continuous} '{printf "%-8s %4d %18.10f%18.10f%18.10f\n",$1,NR+n,$3,$4,$5}' >> ${output}.ol.cc1
    fi
    line=$((${line} + 1 + ${natoms}))
    continuous=$((${continuous} + ${natoms}))
done

rm -f ${tmp}

echo "Done."
exit
