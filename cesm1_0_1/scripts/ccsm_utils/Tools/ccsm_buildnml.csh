#! /bin/csh -fv

echo "-------------------------------------------------------------------------"
#--- Build Namelists 
set buildnml = TRUE
if ($?CCSMBUILDONLY) then
   if ($CCSMBUILDONLY == 'TRUE') then
echo " CCSM BUILDNML NOT BEING RUN, CCSMBUILDONLY flag is TRUE"
      exit 0
   endif
endif

echo " CCSM BUILDNML SCRIPT STARTING"

#echo -------------------------------------------------------------------------
#echo " Preparing $GRID component models for execution "
#echo -------------------------------------------------------------------------
echo " - To prestage restarts, untar a restart.tar file into $RUNDIR"

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

setenv COMPLIBS ""
@ n = 0
foreach model ($MODELS)
  @ n = $n + 1

  #--- Set component specific variables

  set comp =  $COMPONENTS[$n]

  #--- Make necessary directories

  set objdir = $OBJROOT/$model/obj ; if !(-d $objdir) mkdir -p $objdir
  set libdir = $EXEROOT/$model     ; if !(-d $libdir) mkdir -p $libdir

  cd   $libdir

  #--- Build Namelist
  $CASEBUILD/$comp.buildnml.csh
  if ($status != 0) then
    echo  ERROR: $comp.buildnml.csh failed
    exit  99
  endif
end

# -------------------------------------------------------------------------
# Save namelist to docdir
# -------------------------------------------------------------------------

cd $rundir
chmod +w $docdir/*
cp -p *_in $docdir/  >& /dev/null
cp -p *.stream.txt $docdir/ >& /dev/null
cp -p *.stxt $docdir/ >& /dev/null
cp -p *maps.rc $docdir/ >& /dev/null
chmod 444 $docdir/*

echo " CCSM BUILDNML SCRIPT HAS FINISHED SUCCESSFULLY"



