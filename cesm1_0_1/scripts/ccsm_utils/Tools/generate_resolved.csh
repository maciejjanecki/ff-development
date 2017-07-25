#! /bin/csh -f

#==============================================================================
# This scripts generates 
# a) resolved namelist-prestage scripts in the CASEBUILD directory
#    for the selected model components
# b) resolved scripts to generate component executables in the CASEBUILD
#    directory
# c) scripts to generate required ccsm libraries in the CASEBUILD
#    directory
#==============================================================================

# Note: the environment variables needed are passed to this script from the
# calling script

# -------------------------------------------------------------------------
# Ensure that resolved directories do not already exist
# -------------------------------------------------------------------------

if      (-d ${CASEBUILD}) then
   echo ERROR: directory $CASEBUILD already exists 
   echo Must remove this directory in order to generate resolved scripts
   exit -1
endif

# -------------------------------------------------------------------------
# Error checks
# -------------------------------------------------------------------------

if ($RUN_TYPE != 'startup' && $RUN_TYPE != 'hybrid' && $RUN_TYPE != 'branch') then
   echo "ERROR: RUN_TYPE setting of $RUN_TYPE not recognized"
   echo "must be set to [startup, hybrid or branch]"
   exit 1
endif

# -------------------------------------------------------------------------
# Generating resolved setup namelist-prestage and build scripts
# -------------------------------------------------------------------------

cd $CASEROOT
mkdir $CASEBUILD

echo "Generating resolved namelist, prestage, and build scripts"
#v echo -------------------------------------------------------------------
#v echo " Generating resolved setup namelist-prestage and build scripts"
#v echo " and installing build scripts for libraries  "
#v echo "    See directory ${CASEBUILD} "
#v echo -------------------------------------------------------------------

cd $CASETOOLS/Templates
./$COMP_ATM.*template || echo "ERROR: generate_resolved.csh error for atm template"  && exit -1
./$COMP_LND.*template || echo "ERROR: generate_resolved.csh error for lnd template"  && exit -2
./$COMP_ICE.*template || echo "ERROR: generate_resolved.csh error for ice template"  && exit -3
./$COMP_OCN.*template || echo "ERROR: generate_resolved.csh error for ocn template"  && exit -4
./$COMP_GLC.*template || echo "ERROR: generate_resolved.csh error for glc template"  && exit -3
./$COMP_CPL.template  || echo "ERROR: generate_resolved.csh error for cpl template"  && exit -5
./ccsm.template       || echo "ERROR: generate_resolved.csh error for ccsm template" && exit -6

# -------------------------------------------------------------------------
# Installing build scripts for libraries  
# -------------------------------------------------------------------------

cd ${CASEROOT}
foreach LIB (mct csm_share pio)
   cp ${UTILROOT}/Components/$LIB.buildlib ${CASEBUILD}/. || exit 1
end
chmod 755 ${CASEBUILD}/* 

echo "Successfully generated resolved namelist, prestage, and build scripts"

