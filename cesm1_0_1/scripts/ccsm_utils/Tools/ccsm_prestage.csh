#! /bin/csh -fv

# This script prestages as follows
# DIN_LOC_ROOT is the local inputdata area, check it exists
# check whether all the data is in DIN_LOC_ROOT
# if not, check that DIN_LOC_ROOT_CSMDATA exists
#    export any missing data to DIN_LOC_ROOT_CSMDATA
#    prestage the data to DIN_LOC_ROOT
# prestage the REFCASE data if needed
#
# in many cases, DIN_LOC_ROOT equals DIN_LOC_ROOT_CSMDATA
# so a few of the steps are redundant.

if ($?CCSMBUILDONLY) then
   if ($CCSMBUILDONLY == 'TRUE') then
echo " CCSM PRESTAGE NOT BEING RUN, CCSMBUILDONLY flag is TRUE"
      exit 0
   endif
endif

echo "-------------------------------------------------------------------------"
echo " CCSM PRESTAGE SCRIPT STARTING"

# -------------------------------------------------------------------------
# Determining DIN_LOC_ROOT
# -------------------------------------------------------------------------

echo " - CCSM input data directory, DIN_LOC_ROOT_CSMDATA, is $DIN_LOC_ROOT_CSMDATA"
echo " - Case input data directory, DIN_LOC_ROOT, is $DIN_LOC_ROOT"
echo " - Checking the existence of input datasets in DIN_LOC_ROOT"

if !(-d $DIN_LOC_ROOT) then
  echo " "
  echo "  ERROR DIN_LOC_ROOT $DIN_LOC_ROOT does not exist"
  echo " "
  exit -20
endif

if (`./check_input_data -inputdata $DIN_LOC_ROOT -check | grep "unknown" | wc -l` > 0) then
   echo " "
   echo "The following files were not found, this is informational only"
   ./check_input_data -inputdata $DIN_LOC_ROOT -check
   echo " "
endif

if (`./check_input_data -inputdata $DIN_LOC_ROOT -check | grep "missing" | wc -l` > 0) then

   if !(-d $DIN_LOC_ROOT_CSMDATA) then
     echo " "
     echo "  ERROR DIN_LOC_ROOT_CSMDATA $DIN_LOC_ROOT_CSMDATA does not exist"
     echo " "
     exit -20
   endif

   if (`./check_input_data -inputdata $DIN_LOC_ROOT_CSMDATA -check | grep "missing" | wc -l` > 0) then
      echo "Attempting to download missing data:"
      ./check_input_data -inputdata $DIN_LOC_ROOT_CSMDATA -export
   endif

   if (`./check_input_data -inputdata $DIN_LOC_ROOT_CSMDATA -check | grep "missing" | wc -l` > 0) then
      echo " "
      echo "The following files were not found, they are required"
      ./check_input_data -inputdata $DIN_LOC_ROOT_CSMDATA -check
      echo "Invoke the following command to obtain them"
      echo "   ./check_input_data -inputdata $DIN_LOC_ROOT_CSMDATA -export"
      echo " "
      exit -1
   endif

   if (($DIN_LOC_ROOT != $DIN_LOC_ROOT_CSMDATA)) then
      ./check_input_data -inputdata $DIN_LOC_ROOT_CSMDATA -prestage $DIN_LOC_ROOT
      if ($status != 0) then
          echo " "
          echo  "ERROR ccsm_prestage.csh: prestaging failed"
          echo " "
          exit -2
      endif
   endif

endif

if (`./check_input_data -inputdata $DIN_LOC_ROOT -check | grep "missing" | wc -l` > 0) then
   echo " "
   echo "The following files were not found, they are required"
   ./check_input_data -inputdata $DIN_LOC_ROOT -check
   echo "Invoke the following command to obtain them"
   echo "   ./check_input_data -inputdata $DIN_LOC_ROOT_CSMDATA -export"
   echo "AND review the DIN_LOC_ROOT and DIN_LOC_ROOT_CSMDATA env variable settings"
   echo " "
   exit -1
endif


if (($GET_REFCASE == 'TRUE') && ($RUN_TYPE != 'startup') && ($CONTINUE_RUN == 'FALSE')) then
  set refdir = "ccsm4_init/$RUN_REFCASE/$RUN_REFDATE"

### auto-export of refcase turned off
#  if !(-d $DIN_LOC_ROOT_CSMDATA/$refdir) then
#    echo "$DIN_LOC_ROOT_CSMDATA/$refdir is not on local disk, trying to get from data repository"
#    set refdir_head = $refdir:h
#    if !(-d $DIN_LOC_ROOT_CSMDATA/${refdir_head}) then
#       mkdir -p $DIN_LOC_ROOT_CSMDATA/${refdir_head} >& /dev/null
#    endif
#    if (-d $DIN_LOC_ROOT_CSMDATA/${refdir_head}) then
#       cd $DIN_LOC_ROOT_CSMDATA/${refdir_head}
#       svn export https://svn-ccsm-inputdata.cgd.ucar.edu/trunk/inputdata/$refdir
#    endif
#  endif

  if !(-d $DIN_LOC_ROOT_CSMDATA/$refdir) then
    echo "*****************************************************************"
    echo "configure ERROR: $DIN_LOC_ROOT_CSMDATA/$refdir is not on local disk"
    echo "obtain this data from the svn input data repository:"
    echo "  > mkdir -p $DIN_LOC_ROOT_CSMDATA/$refdir"
    echo "  > cd $DIN_LOC_ROOT_CSMDATA/$refdir"
    echo "  > cd .."
    echo "  > svn export --force https://svn-ccsm-inputdata.cgd.ucar.edu/trunk/inputdata/$refdir"
    echo "or set GET_REFCASE to FALSE in env_conf.xml, reconfigure,"
    echo "   and prestage the restart data to $RUNDIR manually"
    echo "*****************************************************************"
    exit -1
  endif 

  echo " - Prestaging REFCASE ($refdir) to $RUNDIR"
  if !(-d $RUNDIR) mkdir -p $RUNDIR || "cannot make $RUNDIR" && exit -1
  cp $DIN_LOC_ROOT_CSMDATA/$refdir/* $RUNDIR || \
        "cannot prestage $DIN_LOC_ROOT_CSMDATA/$refdir to $RUNDIR" && exit -1
  chmod u+w $RUNDIR/*
endif

echo " CCSM PRESTAGE SCRIPT HAS FINISHED SUCCESSFULLY"
