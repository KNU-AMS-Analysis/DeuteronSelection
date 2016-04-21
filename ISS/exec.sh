#!/bin/bash
export USERNAME=`whoami`
export PWD=`pwd`
source /afs/cern.ch/user/${USERNAME:0:1}/$USERNAME/amsvar_gcc.sh
source /afs/cern.ch/user/${USERNAME:0:1}/$USERNAME/AMS-ACsoft/scripts/thisacsoft.sh

$PWD/bin/main $1 $2
