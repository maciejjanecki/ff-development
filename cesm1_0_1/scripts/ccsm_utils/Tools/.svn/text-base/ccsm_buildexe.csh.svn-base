#! /bin/csh -fv

echo "-------------------------------------------------------------------------"
echo " CCSM BUILDEXE SCRIPT STARTING"

set atmstr = 0
if ($NTHRDS_ATM > 1 || $BUILD_THREADED == "TRUE") set atmstr = 1
set lndstr = 0
if ($NTHRDS_LND > 1 || $BUILD_THREADED == "TRUE") set lndstr = 1
set icestr = 0
if ($NTHRDS_ICE > 1 || $BUILD_THREADED == "TRUE") set icestr = 1
set ocnstr = 0
if ($NTHRDS_OCN > 1 || $BUILD_THREADED == "TRUE") set ocnstr = 1
set glcstr = 0
if ($NTHRDS_GLC > 1 || $BUILD_THREADED == "TRUE") set glcstr = 1
set cplstr = 0
if ($NTHRDS_CPL > 1 || $BUILD_THREADED == "TRUE") set cplstr = 1

set smpstr = a${atmstr}l${lndstr}i${icestr}o${ocnstr}g${glcstr}c${cplstr}
./xmlchange -file env_build.xml -id SMP_VALUE -val ${smpstr}
setenv SMP_VALUE ${smpstr}

if ($SMP_BUILD != $SMP_VALUE && $SMP_BUILD != "0") then
  echo " "
  echo "  ERROR SMP STATUS HAS CHANGED"
  echo "    SMP_BUILD = $SMP_BUILD"
  echo "    SMP_VALUE = $SMP_VALUE"
  echo "    A manual clean of your obj directories is strongly recommendend"
  echo "    You should execute the following:"
  echo "      ./$CASE.$MACH.clean_build"
  echo "    Then rerun the build script interactively"
  echo "    ---- OR ----"
  echo "    You can override this error message at your own risk by executing"
  echo "      ./xmlchange -file env_build.xml -id SMP_BUILD -val 0 "
  echo "    Then rerun the build script interactively"
  echo " "
  exit -20
endif

if ($BUILD_STATUS != "0") then
  echo " "
  echo "  ERROR env_build or Macros FILE HAS CHANGED"
  echo "    A manual clean of your obj directories is strongly recommendend"
  echo "    You should execute the following:"
  echo "      ./$CASE.$MACH.clean_build"
  echo "    Then rerun the build script interactively"
  echo "    ---- OR ----"
  echo "    You can override this error message at your own risk by executing"
  echo "      rm LockedFiles/Macros.$MACH*"
  echo "      rm LockedFiles/env_build*"
  echo "    Then rerun the build script interactively"
  echo " "
  exit -20
endif

if ($COMP_INTERFACE == "ESMF" && $USE_ESMF_LIB != "TRUE") then
  echo " "
  echo "  ERROR COMP_INTERFACE IS ESMF BUT USE_ESMF_LIB IS NOT TRUE"
  echo "    Set USE_ESMF_LIB to TRUE with"
  echo "      ./xmlchange -file env_build.xml -id USE_ESMF_LIB -val TRUE "
  echo " "
  exit -20
endif

# -------------------------------------------------------------------------
# Determine if models must be rebuilt
# -------------------------------------------------------------------------

if !($?LID) setenv LID "`date +%y%m%d-%H%M%S`"

# -------------------------------------------------------------------------
# Create main execution directories"
# -------------------------------------------------------------------------

foreach DIR ( $EXEROOT $LIBROOT $INCROOT $OBJROOT $RUNDIR)
  if !(-d $DIR) mkdir -p $DIR || "cannot make $DIR" && exit -1
end
set rundir = $RUNDIR
set docdir = $CASEROOT/CaseDocs
if !(-d $docdir) then
  mkdir -p $docdir
  echo "  CCSM Resolved Namelist Files" >& $docdir/README
  echo "    For documentation only" >>& $docdir/README
  echo "    DO NOT MODIFY" >>& $docdir/README
endif
cd $EXEROOT

# -------------------------------------------------------------------------
# Build Libraries:
# -------------------------------------------------------------------------

  if ($USE_MPISERIAL == TRUE) then
    cp -p -f $CODEROOT/utils/mct/mpi-serial/*.h   $LIBROOT/include/.
  endif

  set blibs = "mct pio csm_share"
  echo " - Build Libraries: $blibs "
  setenv SMP "FALSE"
  if ($NTHRDS_ATM > 1 || $NTHRDS_CPL > 1 || $NTHRDS_ICE > 1 || \
      $NTHRDS_LND > 1 || $NTHRDS_OCN > 1 || $NTHRDS_GLC > 1) setenv SMP "TRUE"
  foreach lib ($blibs)
     set libdir = $EXEROOT/$lib; if !(-d $libdir) mkdir -p $libdir
     set file_build = ${libdir}/$lib.bldlog.$LID
     cd $libdir
     touch ${file_build}
     echo `date` $file_build | tee -a ${file_build}
     $CASEBUILD/$lib.buildlib >>& ${file_build}
     if ($status != 0) then
        echo  ERROR: $lib.buildlib failed, see ${file_build}
        echo  ERROR: cat ${file_build}
        exit  99
     endif
  end

# -------------------------------------------------------------------------

setenv COMPLIBS ""
@ n = 0
foreach model ($MODELS)
  @ n = $n + 1

  #--- Set component specific variables

  set comp =  $COMPONENTS[$n]
  setenv SMP FALSE
  if ($NTHRDS[$n] > 1) setenv SMP TRUE

  if ($NTASKS[$n] < 1) then
     echo "$model $comp has <1 NTASKS, must be greater than 0"
     exit 99
  endif
  if ($NTHRDS[$n] < 1) then
     echo "$model $comp has <1 NTHRDS, must be greater than 0"
     exit 99
  endif

  #--- Make necessary directories

  set objdir = $OBJROOT/$model/obj ; if !(-d $objdir) mkdir -p $objdir
  set libdir = $EXEROOT/$model     ; if !(-d $libdir) mkdir -p $libdir

  cd   $libdir
  set file_build = ${rundir}/$model.bldlog.$LID
  touch ${file_build}
  echo `date` ${file_build} | tee -a ${file_build}

  #--- Build Model Executables
  $CASEBUILD/$comp.buildexe.csh  >>& ${file_build}
  if ($status != 0) then
    echo  ERROR: $comp.buildexe.csh failed, see ${file_build}
    echo  ERROR: cat ${file_build}
    exit  99
  endif

  #--- Copy .mod files and add library name to env var COMPLIBS
  #--- Need both cases in comp because of lightning
  cp -p $objdir/*_comp_*.mod $INCROOT/ >& /dev/null
  cp -p $objdir/*_COMP_*.mod $INCROOT/ >& /dev/null

end

  #--- Go into executable directory
  set model = ccsm
  set comp  = ccsm
  setenv SMP "FALSE"
  if ($NTHRDS_ATM > 1 || $NTHRDS_CPL > 1 || $NTHRDS_ICE > 1 || \
      $NTHRDS_LND > 1 || $NTHRDS_OCN > 1 || $NTHRDS_GLC > 1) setenv SMP "TRUE"

  set objdir = $OBJROOT/$model/obj ; if !(-d $objdir) mkdir -p $objdir
  set libdir = $EXEROOT/$model     ; if !(-d $libdir) mkdir -p $libdir

  #--- Go into executable directory
  cd   $libdir
  set file_build = ${rundir}/$model.bldlog.$LID
  touch ${file_build}
  echo `date` ${file_build} | tee -a ${file_build}

  #--- Build Model Executables or Libraries
  $CASEBUILD/$comp.buildexe.csh  >>& ${file_build}
  if ($status != 0) then
    echo  ERROR: $comp.buildexe.csh failed, see ${file_build}
    echo  ERROR: cat ${file_build}
    exit  99
  endif

  #--- rename and save the latest executable
  cp  $rundir/ccsm.exe $EXEROOT/$CASE.ccsm.exe.$LID
  cmp -s $EXEROOT/$CASE.ccsm.exe $EXEROOT/$CASE.ccsm.exe.$LID
  if ($status == 0) then
    rm -f $EXEROOT/$CASE.ccsm.exe.$LID
  else
    cp $EXEROOT/$CASE.ccsm.exe.$LID $EXEROOT/$CASE.ccsm.exe
  endif

endif

# -------------------------------------------------------------------------
# Save model output stdout and stderr 
# -------------------------------------------------------------------------

if ($LOGDIR != "") then
  cd $EXEROOT
  gzip */*.bldlog.$LID
  if (! -d $LOGDIR) mkdir -p $LOGDIR
  if (! -d $LOGDIR/bld) mkdir -p $LOGDIR/bld
  cp -p */*bldlog.$LID.* $LOGDIR/bld/ || echo "Error in copy of build logs "
endif

# -------------------------------------------------------------------------
# Go into case directory and set BUILD_COMPLETE flag to TRUE
# -------------------------------------------------------------------------

   cd $CASEROOT
   ./xmlchange -file env_build.xml -id BUILD_COMPLETE -val TRUE
   ./xmlchange -file env_build.xml -id BUILD_STATUS -val 0
   if ($status != 0) then
     echo  ERROR: ccsm_buildexe.csh call to xmlchange failed
     exit  100
   endif
   ./xmlchange -file env_build.xml -id SMP_BUILD -val $SMP_VALUE
   if ($status != 0) then
     echo  ERROR: ccsm_buildexe.csh call to xmlchange failed
     exit  100
   endif

   cd $CASEROOT
   rm LockedFiles/env_build*  >& /dev/null
   rm LockedFiles/Macros.$MACH*  >& /dev/null
   foreach file (env_build.xml Macros.$MACH)    
      if (! -e LockedFiles/$file.locked) then
         cp $file LockedFiles/$file.locked
         if ($status != 0) then
            echo  ERROR: ccsm_buildexe.csh problem in locking $file
            exit  100
         endif
         echo " - Locking file $file"
      endif
   end

echo " CCSM BUILDEXE SCRIPT HAS FINISHED SUCCESSFULLY"



