#!/bin/bash
# Wooyoung Jang
# Kyungpook National University
# Date : 2015-12-28
# Last update : 2016-01-13
# Purpose : 
#       This script helps to submit large amount of jobs to Portable Batch System.

# Print out current time
date -u

USRNAME=`whoami`
WORKDIR=`pwd`
#RUNSCRIPT="/u/user/$USRNAME/public/pass6/run.sh"
RUNSCRIPT="/u/user/$USRNAME/batch_job/MC_d/run.sh"
#OUTDIR="/u/user/$USRNAME/public/pass6/output"
OUTDIR="/u/user/$USRNAME/batch_job/MC_d/output"
LOGDIR=$OUTDIR/$1/log
ERRDIR=$OUTDIR/$1/err
#TARGETDIR="/u/user/wyjang/public/pass6"
TARGETDIR="/u/user/sinchul/background_study/Deuteron_MC/Data"

# Check existence of the RUNSCRIPT.
if [ -e $RUNSCRIPT ]; then
  echo "Pre-run script file is found ... $RUNSCRIPT"
else
  echo "Can not find $RUNSCRIPT"
fi

# Check the existence of the output directory.
if [ -d $OUTDIR ]; then
  echo "Output files will be located in $OUTDIR."
else
  mkdir $OUTDIR
  echo "Warning: $OUTDIR is made to store output files."
fi

# Check the existence of the output sub-directory.
if [ -d $OUTDIR/$1 ]; then
  echo "Submit for $1 in $OUTDIR/$1"
else
  mkdir $OUTDIR/$1
  echo "Warning: $OUTDIR/$1 is made to store output files."
fi

# Check the existence of the log-file directory.
if [ -d $LOGDIR ]; then
    echo "Log files will be stored at $LOGDIR"
else
  mkdir $LOGDIR
  echo "Warning: $LOGDIR is made to store log files."
fi

# Check the existence of the error-file directory.
if [ -d $ERRDIR ]; then
  echo "Error logs will be stored at $ERRDIR"
else
  mkdir $ERRDIR
  echo "Warning: $ERRDIR is made to store error logs."
fi

# Make a temporary file list with formating like
# /u/user/wyjang/public/pass6/set_1105/1305853512.root (full path of the target file)
find $TARGETDIR/$1/*.root >> pathlist_$1.tmp

# Loop over the whole target files to submit job
while read irun
do
  # Define the name of log, error and output file names.
  LOGFILE=$LOGDIR/${irun#$TARGETDIR/$1/}.log    # The symbol sharp '#' deletes the next following pattern.
  ERRFILE=$ERRDIR/${irun#$TARGETDIR/$1/}.err
  OUTFILE=$OUTDIR${irun#$TARGETDIR}

  # Execute job submission command.
  #echo /usr/bin/qsub -q knu -o $LOGFILE -e $ERRFILE -l walltime=24:00:00,cput=24:00:00 -N ${irun#$OUTDIR} -F \"$irun $OUTFILE\" $RUNSCRIPT
  /usr/bin/qsub -q knu -o $LOGFILE -e $ERRFILE -l walltime=24:00:00,cput=24:00:00 -N $irun -F "$irun $OUTFILE" $RUNSCRIPT
done < "pathlist_$1.tmp"

# Delete the temporary file list
rm pathlist_$1.tmp
