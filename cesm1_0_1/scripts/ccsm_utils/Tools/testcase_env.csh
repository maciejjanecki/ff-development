#!/bin/csh -f

#======================================================================
# All environment variables common to all tests that are only used by 
# the test cases should be put here.  
#======================================================================

# regression test baseline directory and coupler log files
# TBD  TBH:  add names of coupler history files here later...

if ($?casebaseid) then
  echo "setenv CASEBASEID $casebaseid" >> $CASE.test || exit -1
else
  echo "setenv CASEBASEID $CASE" >> $CASE.test || exit -1
endif

if ( $?BASELINE_ROOT ) then
   # continue
else if ( $?CCSM_BASELINE) then
   setenv BASELINE_ROOT $CCSM_BASELINE
else
   setenv BASELINE_ROOT "/UNSET"
endif
echo "setenv BASELINE_ROOT       $BASELINE_ROOT"       >> $CASE.test || exit -1

# generate baseline test flag
if ( $?GENERATE_BASELINE ) then
  echo "setenv GENERATE_BASELINE"                        >> $CASE.test || exit -1
  echo "setenv BASEGEN_NAME        $baseline_name_gen"   >> $CASE.test || exit -1
  echo "setenv BASEGEN_DIR         $BASELINE_ROOT/$BASEGEN_CASE"           >> $CASE.test || exit -1
  echo "setenv BASEGEN_CPLLOGFILE  $BASELINE_ROOT/$BASEGEN_CASE/cpl.log"   >> $CASE.test || exit -1
  echo "setenv BASEGEN_CPLHISTFILE $BASELINE_ROOT/$BASEGEN_CASE/cpl.hi.nc" >> $CASE.test || exit -1
  echo "setenv BASEGEN_CPLPROFFILE $BASELINE_ROOT/$BASEGEN_CASE/timing_summary" >> $CASE.test || exit -1
else
  echo "unsetenv GENERATE_BASELINE"   >> $CASE.test || exit -1
  echo "unsetenv BASEGEN_DIR"         >> $CASE.test || exit -1
  echo "unsetenv BASEGEN_CPLLOGFILE"  >> $CASE.test || exit -1
  echo "unsetenv BASEGEN_CPLHISTFILE" >> $CASE.test || exit -1
  echo "unsetenv BASEGEN_CPLPROFFILE" >> $CASE.test || exit -1
endif

# regression test comparison flag
if ( $?COMPARE_BASELINE ) then
  echo "setenv COMPARE_BASELINE"                         >> $CASE.test || exit -1
  echo "setenv BASECMP_NAME        $baseline_name_cmp"   >> $CASE.test || exit -1
  echo "setenv BASECMP_DIR         $BASELINE_ROOT/$BASECMP_CASE"           >> $CASE.test || exit -1
  echo "setenv BASECMP_CPLLOGFILE  $BASELINE_ROOT/$BASECMP_CASE/cpl.log"   >> $CASE.test || exit -1
  echo "setenv BASECMP_CPLHISTFILE $BASELINE_ROOT/$BASECMP_CASE/cpl.hi.nc" >> $CASE.test || exit -1
else
  echo "unsetenv COMPARE_BASELINE"    >> $CASE.test || exit -1
  echo "unsetenv BASECMP_DIR"         >> $CASE.test || exit -1
  echo "unsetenv BASECMP_CPLLOGFILE"  >> $CASE.test || exit -1
  echo "unsetenv BASECMP_CPLHISTFILE" >> $CASE.test || exit -1
endif

# cleanup option
if ( $?CLEANUP ) then
  echo "setenv CLEANUP" >> $CASE.test || exit -1
else
  echo "unsetenv CLEANUP" >> $CASE.test || exit -1
endif

