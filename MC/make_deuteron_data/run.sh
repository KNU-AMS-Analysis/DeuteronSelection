#!/bin/sh

#export ROOTSYS=/cvmfs/ams.cern.ch/Offline/AMSsoft/linux_slc6_icc64/root534
#export LD_LIBRARY_PATH=/lib:/cvmfs/ams.cern.ch/Offline/intel/composer_xe_2013_sp1.3.174/compiler/lib/intel64:/cvmfs/ams.cern.ch/Offline/AMSsoft/linux_slc6_icc64/xrootd/lib64:/cvmfs/ams.cern.ch/Offline/AMSsoft/linux_slc6_icc64/root534/lib:/cvmfs/ams.cern.ch/Offline/AMSsoft/linux_slc6_icc64/geant4_ams/lib/geant4/Linux-icc:/u/user/sinchul/usr/local/lib
#export PATH=$ROOTSYS/bin:$PATH
source /u/user/sinchul/.bash_profile

echo "Executing commmand ... root -b -q -l 'test.cpp(\"$1\", \"$2\")'"
#/cvmfs/ams.cern.ch/Offline/AMSsoft/linux_slc6_icc64/root534/bin/root -b -q -l '/u/user/wyjang/public/pass6/test.cpp('\"$1\"', '\"$2\"')'
#/cvmfs/ams.cern.ch/Offline/AMSsoft/linux_slc6_icc64/root534/bin/root -b -q -l '/u/user/sinchul/batch_job/test.cpp('\"$1\"', '\"$2\"')'
#/cvmfs/ams.cern.ch/Offline/AMSsoft/linux_slc6_icc64/root534/bin/root -b -q -l '/u/user/sinchul/batch_job/MC_d/d_BDT.C('\"$1\"', '\"$2\"')'
#/cvmfs/ams.cern.ch/Offline/AMSsoft/linux_slc6_icc64/root534/bin/root -b -q -l '/u/user/sinchul/batch_job/MC_d/make_data.C('\"$1\"', '\"$2\"')'
#/cvmfs/ams.cern.ch/Offline/AMSsoft/linux_slc6_icc64/root534/bin/root -b -q -l '/u/user/sinchul/batch_job/MC_d/iss_BDT.C('\"$1\"', '\"$2\"')'
#/cvmfs/ams.cern.ch/Offline/AMSsoft/linux_slc6_icc64/root534/bin/root -b -q -l '/u/user/sinchul/batch_job/MC_d/bdt_eachbin.C('\"$1\"', '\"$2\"')'
#/cvmfs/ams.cern.ch/Offline/AMSsoft/linux_slc6_icc64/root534/bin/root -b -q -l '/u/user/sinchul/batch_job/MC_d/make_data_raw.C('\"$1\"', '\"$2\"')'
/storage/Software/root/6.06.02/x86_64-sl6-gcc/bin/root -b -q -l '/u/user/sinchul/batch_job/MC_d/make_data.C('\"$1\"', '\"$2\"')'
#/storage/Software/root/6.06.02/bin/root -b -q -l '/u/user/sinchul/batch_job/MC_d/make_accept_data.C('\"$1\"', '\"$2\"')'
