#!/bin/csh -f
#
while ($#argv > 0)
   set val=$argv[1]
   shift
   switch ($val)
      case -y:
        set OVERWRITE=true
        breaksw
      case -p:
        set PARTITION=$argv[1]
        shift
        breaksw
      case -l:
        set LOGFILE=$argv[1]
        shift
        breaksw
      case -n:
        set NCPUS=$argv[1]
        shift
        breaksw
      case -w:
        set WALL=$argv[1]
        shift
        breaksw
      case -ppn:
        set PPN=$argv[1]
        shift
        breaksw
      case -v:
        set VERNO=$argv[1]
        shift
        breaksw
      case -exepath:
        set XPATH=$argv[1]
        shift
        breaksw
      case -restart:
        set RESTARTPATH=$argv[1]
        shift
        breaksw
      default:
        if ($?INPUT == 1) then
           echo "You have given too many input file names, $INPUT and $val."
           exit 4
        else
           set INPUT=$val
           if ($INPUT:r.inp != $INPUT) set INPUT=$1.inp
           if (! -f $INPUT) then
              echo "Input file $INPUT does not exist."
              exit 4
           endif
        endif
        breaksw
   endsw
end
#
if ($?INPUT == 0)       set INPUT=help
if ($?VERNO == 0)       set VERNO=00
if ($?LOGFILE == 0)     set LOGFILE=default
if ($?XPATH == 0)       set XPATH=$GMSPATH
if ($?RESTARTPATH == 0) set RESTARTPATH=/home/$USER/restart
if ($?NCPUS == 0)       set NCPUS=1
if ($?WALL == 0)        set WALL=default
if ($?PPN == 0)         set PPN=$NCPUS
if ($?PARTITION == 0)   set PARTITION=compute
if ($?OVERWRITE == 0)   set OVERWRITE=false
#
set JOB=`basename $INPUT`
set FULL_PATH=`readlink -f $INPUT`
set JOB_PATH=`dirname $FULL_PATH`
if ($JOB:r.inp == $JOB) set JOB=$JOB:r
#
if ($JOB == help) then
   tput clear
   echo "       ::::::::      :::       :::   :::   :::::::::: ::::::::   :::::::: "
   echo "     :+:    :+:   :+: :+:    :+:+: :+:+:  :+:       :+:    :+: :+:    :+: "
   echo "    +:+         +:+   +:+  +:+ +:+:+ +:+ +:+       +:+        +:+         "
   echo "   :#:        +#++:++#++: +#+  +:+  +#+ +#++:++#  +#++:++#++ +#++:++#++   "
   echo "  +#+   +#+# +#+     +#+ +#+       +#+ +#+              +#+        +#+    "
   echo " #+#    #+# #+#     #+# #+#       #+# #+#       #+#    #+# #+#    #+#     "
   echo " ########  ###     ### ###       ### ########## ########   ########       "
   echo " "
   echo "The syntax to execute GAMESS is"
   echo " "
   echo "   gms      [-l logfile] [-n CPUS] [-w dd:hh:mm:ss] [input]"
   echo " "
   echo " Required arguments"
   echo " "
   echo "   -l        specify log file name or full absolute path to the log file"
   echo "   -n CPUS   specify total number of GAMESS compute processes for this run"
   echo "   -w        specify wall clock limit as dd:hh:mm:ss (default=24:00:00=24 hrs)"
   echo " "
   echo " Optional arguments"
   echo " "
   echo "   -gmspath  specify full absolute path to the folder containing GAMESS executable,"
   echo "             install.info, and rungms-bolt script"
   echo "   -restart  specify full absolute path to use for the restart folder (default=\$HOME/restart)"
   echo "   -ppn CPN  specify the number of GAMESS compute processes in each node"
   echo "   -p        specify partition to use (default=compute)"
   echo "             ------------------------------------------------"
   echo "             Partition    - Description"
   echo "             ------------------------------------------------"
   echo "             compute      - AMD large memory (768 GB) node   (4 available)"
   echo "             grace-hopper - NVIDIA GH200 Grace + Hopper node (3 available)"
   echo "             gpu          - 4-way H100 GPU node              (1 available)"
   echo " "
   echo "   -v        specify GAMESS version to use (default=00)"
   echo "   -y        overwrite previous log file if it found"
   echo " "
   exit
endif
#
if ($LOGFILE == default) then
   set LOGFILE=$JOB.log
   echo -n "Output file name? [$LOGFILE] "
   set LOGFILE=`pwd`/$JOB.log
   set ans=$<
   if (null$ans != null) set LOGFILE=$ans
else
   if ($LOGFILE == $JOB.log) set LOGFILE=`pwd`/$JOB.log
endif
#
if (-e $LOGFILE) then
   if ("$OVERWRITE" == "false") then
     echo -n "$LOGFILE already exists.  OK to delete the old one? [y] "
     set ans=$<
     if (null$ans == null) set ans=y
     if ($ans == y) then
       rm $LOGFILE
     else
       echo "Exiting, so you can think about your old log file's value."
       exit
     endif
   else
     echo "OVERWRITE MODE! Removing existing log file $LOGFILE"
     rm $LOGFILE
   endif
endif
#
if ($NCPUS == 0) then
   echo -n "Total number of GAMESS compute processes for ths run? [$NCPUS] "
   set ans=$<
   if (null$ans != null) set NCPUS=$ans
endif
#
if ($NCPUS < $PPN) then
   set NNODES=1
else
   @ xx = $NCPUS / $PPN
   @ yy = $PPN * $xx
   set NCPUS=$yy
   set NNODES=$xx
endif
#
@ zz = $PPN * 2
set PPN2 = $zz
unset xx
unset yy
unset zz
#
if ($WALL == default) then
   set WALL=24:00:00
endif
#
if ($NCPUS < $PPN) set PPN=$NCPUS
#
set USERSCR=$RESTARTPATH
if (-d $USERSCR) then
   echo "Using $USERSCR to store restart files for this run."
else
  echo "An empty $USERSCR directory is being created for you ..."
  mkdir -p $USERSCR
  if (-d $USERSCR) then
    echo "This directory will receive .dat, .trj, .rst supplemental outputs."
  else
    echo "Problem creating $USERSCR, no job can be submitted."
    exit
  endif
endif
#
set DATESTRING=`date "+%Y.%m.%d_%H.%M.%S.%N"`
set LONGJOBNAME=`echo ${DATESTRING}_${JOB} | sed "s/\//_/g"`
cp $XPATH/machines/internal/rungms-bolt $USERSCR/$LONGJOBNAME.csh
#
if ($XPATH != none) then
   setenv GMSPATH $XPATH
endif
#
echo " "
echo "-----------"
echo "JOB SUMMARY"
echo "==========="
echo "JOB NAME   : $JOB"
echo "INPUT      : $FULL_PATH"
echo "LOG FILE   : $LOGFILE"
echo "GAMESS PATH: $GMSPATH"
echo "BINARY     : gamess.$VERNO.x"
echo "RESTART    : $RESTARTPATH"
echo "# NODES    : $NNODES"
echo "TASKS/NODE : $PPN2"
echo "TIME LIMIT : $WALL"
echo " "
echo "Submitting GAMESS job $JOB.inp"
set EXTRAOPTIONS="--chdir=$XPATH"
set EXTRAOPTIONS="$EXTRAOPTIONS --exclusive"
if ($PARTITION != none) set EXTRAOPTIONS="$EXTRAOPTIONS --partition=$PARTITION"
chmod 755 $USERSCR/$LONGJOBNAME.csh
echo "#\!/bin/csh"                                                             > $USERSCR/$LONGJOBNAME.slurmjob
echo "setenv GMSPATH $GMSPATH"                                                >> $USERSCR/$LONGJOBNAME.slurmjob
echo "setenv RESTARTPATH $RESTARTPATH"                                        >> $USERSCR/$LONGJOBNAME.slurmjob
echo "setenv AUXPATH /shared/gamess/auxdata"                                  >> $USERSCR/$LONGJOBNAME.slurmjob
echo "$USERSCR/$LONGJOBNAME.csh $FULL_PATH $VERNO $NCPUS $PPN >& $LOGFILE"    >> $USERSCR/$LONGJOBNAME.slurmjob
echo "rm $USERSCR/$LONGJOBNAME.csh"                                           >> $USERSCR/$LONGJOBNAME.slurmjob
(set echo; sbatch --export=PATH -o $LOGFILE --no-requeue --job-name=$JOB --nodes=$NNODES --time=$WALL --ntasks-per-node=$PPN2 $EXTRAOPTIONS $USERSCR/$LONGJOBNAME.slurmjob)
rm $USERSCR/$LONGJOBNAME.slurmjob
exit
