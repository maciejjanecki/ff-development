#! /bin/csh -f

#==============================================================================
# This script generates resolved run and long term archiving batch scripts 
# it adds the correct batch settings for the requested machine given the tasks and
# threads required for the specified machine
#==============================================================================

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# Generate clean_build script
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------

echo  "Generating clean_build script "

cp -f ./Tools/clean_build ./${CASE}.${MACH}.clean_build
chmod 775 ${CASEROOT}/${CASE}.${MACH}.clean_build

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# Generate submit script
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------

echo  "Generating submit script "

cat >> ${CASEROOT}/${CASE}.${MACH}.submit << EOF1
#! /bin/csh -f

./check_case || echo "check_case failed" && exit -99

${BATCHSUBMIT} ${CASE}.${MACH}.run

set sdate = \`date +"%Y-%m-%d %H:%M:%S"\`
echo "run submitted \$sdate" >>& CaseStatus

EOF1
chmod 775 ${CASEROOT}/${CASE}.${MACH}.submit

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# Generate build script
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------

echo  "Generating build script "

# Create the build file
if    (-e ${CASEROOT}/${CASE}.${MACH}.build) then
   rm     ${CASEROOT}/${CASE}.${MACH}.build
endif
touch     ${CASEROOT}/${CASE}.${MACH}.build
chmod 775 ${CASEROOT}/${CASE}.${MACH}.build

cat >> ${CASEROOT}/${CASE}.${MACH}.build << EOF1
#! /bin/csh -f

./Tools/ccsm_check_lockedfiles || exit -1

source ./Tools/ccsm_getenv || exit -2

set SNAME = \$0 ; set SNAME = \(\$SNAME:t\)
if (\$USE_MPISERIAL == "TRUE" && \$MPISERIAL_SUPPORT == "FALSE") then
   echo "\$SNAME ERROR: USE_MPISERIAL == TRUE not supported on ${MACH}"
   echo "\$SNAME set USE_MPISERIAL = FALSE in env_conf.xml"
   exit -1
endif

setenv LID "\`date +%y%m%d-%H%M%S\`"

cd \$CASEROOT
source \$CASETOOLS/ccsm_buildnml.csh || exit -3
cd \$CASEROOT
source \$CASETOOLS/ccsm_prestage.csh || exit -3
cd \$CASEROOT
source \$CASETOOLS/ccsm_buildexe.csh || exit -3

set sdate = \`date +"%Y-%m-%d %H:%M:%S"\`
echo "build complete \$sdate" >>& \$CASEROOT/CaseStatus

EOF1

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# Generate batch run script
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------

echo  "Generating run script "

# Create the run file
touch ${CASEROOT}/${CASE}.${MACH}.run
chmod 775 ${CASEROOT}/${CASE}.${MACH}.run

# batch stuff which needs to be substituted now
setenv LBQUERY "TRUE"
if !($?BATCHQUERY) then
  setenv LBQUERY "FALSE"
  setenv BATCHQUERY "undefined"
else if ( "$BATCHQUERY" == 'UNSET' ) then
  setenv LBQUERY "FALSE"
  setenv BATCHQUERY "undefined"
endif

setenv LBSUBMIT "TRUE"
if !($?BATCHSUBMIT) then
  setenv LBSUBMIT "FALSE"
  setenv BATCHSUBMIT "undefined"
else if ( "$BATCHSUBMIT" == 'UNSET' ) then
  setenv LBSUBMIT "FALSE"
  setenv BATCHSUBMIT "undefined"
endif
   
# -------------------------------------------------------------------------

env PHASE=set_batch ${UTILROOT}/Machines/mkbatch.${MACH}
${UTILROOT}/Tools/taskmaker.pl -document >> ${CASEROOT}/${CASE}.${MACH}.run

# -------------------------------------------------------------------------

cat >> ${CASEROOT}/${CASE}.${MACH}.run << EOF1
#-----------------------------------------------------------------------
# Determine necessary environment variables
#-----------------------------------------------------------------------

cd $CASEROOT

./Tools/ccsm_check_lockedfiles || exit -1
source ./Tools/ccsm_getenv || exit -2

if (\$BUILD_COMPLETE != "TRUE") then
  echo "BUILD_COMPLETE is not TRUE"
  echo "Please rebuild the model interactively via"
  echo "   ./\${CASE}.\${MACH}.build"
  exit -2
endif

setenv LBQUERY  $LBQUERY
setenv LBSUBMIT $LBSUBMIT

#-----------------------------------------------------------------------
# Determine time-stamp/file-ID string
# Clean up previous run timing files
#-----------------------------------------------------------------------

setenv LID "\`date +%y%m%d-%H%M%S\`"
env | egrep '(MP_|LOADL|XLS|FPE|DSM|OMP|MPC)' # document env vars

# -------------------------------------------------------------------------
# Build the namelists and check prestage
# -------------------------------------------------------------------------

cd \$CASEROOT
source \$CASETOOLS/ccsm_buildnml.csh || exit -3
cd \$CASEROOT
source \$CASETOOLS/ccsm_prestage.csh || exit -3

# -------------------------------------------------------------------------
# Create and cleanup the timing directories
# -------------------------------------------------------------------------

if !(-d \$RUNDIR/timing) mkdir \$RUNDIR/timing
if !(-d \$RUNDIR/timing/checkpoints) mkdir \$RUNDIR/timing/checkpoints
rm -f \$RUNDIR/timing/ccsm_timing*

set sdate = \`date +"%Y-%m-%d %H:%M:%S"\`
echo "run started \$sdate" >>& \$CASEROOT/CaseStatus

EOF1

# -------------------------------------------------------------------------

env PHASE=set_exe ${UTILROOT}/Machines/mkbatch.${MACH}

cat >> ${CASEROOT}/${CASE}.${MACH}.run << EOF1

cd \$CASEROOT
./Tools/ccsm_postrun.csh || exit 1

EOF1

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# Generate batch long term archiving script if appropriate
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------

env PHASE=set_larch ${UTILROOT}/Machines/mkbatch.${MACH}

