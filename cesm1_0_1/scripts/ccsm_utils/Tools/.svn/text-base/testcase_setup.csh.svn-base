#!/bin/csh -f

# Need the following environment variables to be available
# $SCRIPTSROOT $CASEROOT $CASE $TESTCASE $MACH 
# $GENERATE_BASELINE, $COMPARE_BASELINE, $STOP_OPTION, $STOP_N

set currdir = `pwd`

cd ${CASEROOT} || exit -1

source ./Tools/ccsm_getenv || exit -2

#===================================================================
# Modify env_* in test directory 
#===================================================================

if ( $?GENERATE_BASELINE || $?COMPARE_BASELINE ) then
  ./xmlchange -file env_run.xml -id HIST_OPTION -val '$STOP_OPTION'
  ./xmlchange -file env_run.xml -id HIST_N      -val '$STOP_N'
endif

#===================================================================
# Create testing script and run script
#===================================================================

# create run script (that testing script will call)
#-------------------------------------------------------------
./configure -case4test 
if ($status != 0) then
  echo "ERROR: testcase_setup.csh invokcation of configure -case4test failed"
  exit -1
endif

if (-d $CASEBUILD) then
  # do nothing
else
  ./configure -namelist 
  if ($status != 0) then
    echo "ERROR: testcase_setup.csh invocation of configure -namelist failed"
    exit -1
  endif
endif

# create testing script
#-------------------------------------------------------------
env PHASE=set_batch ${UTILROOT}/Machines/mkbatch.${MACH} || exit -1  

cd $CASEROOT || exit -1
mv $CASE.$MACH.run $CASE.test
# set environment variables for testing...
source ${SCRIPTSROOT}/ccsm_utils/Tools/testcase_env.csh

# fill in the rest of the test script...
cat >> ${CASE}.test << EOF

cd $CASEROOT
./Tools/ccsm_check_lockedfiles || exit -2
source ./Tools/ccsm_getenv || exit -3

EOF

cat ${SCRIPTSROOT}/ccsm_utils/Tools/testcase_begin >> $CASE.test || exit -1
cat ${SCRIPTSROOT}/ccsm_utils/Testcases/${TESTCASE}_script >> $CASE.test || exit -1
cat ${SCRIPTSROOT}/ccsm_utils/Tools/testcase_end >> $CASE.test || exit -1

# copy over compare script
#-------------------------------------------------------------
cp ${SCRIPTSROOT}/ccsm_utils/Tools/check_exactrestart.pl ./Tools/ || exit -1
cp ${SCRIPTSROOT}/ccsm_utils/Tools/hist_compare.csh ./Tools/ || exit -1

chmod 755 $CASE* *pl

cd $currdir

