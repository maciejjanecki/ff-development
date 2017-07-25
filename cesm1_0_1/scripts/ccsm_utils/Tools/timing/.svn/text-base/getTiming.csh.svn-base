#! /bin/csh -f
# under your ccsm running script directory run the following command
# usage:  getTiming.csh -mach machine_name -lid yymmdd-hhmmss
# The output will be a local file

if ($#argv < 4)  then
  echo "Usage:"
  echo '       ./Tools/getTiming.csh -mach machine_name -lid yymmdd-hhmmss'
  echo '         Run from the $CASEROOT'
  echo '         Assumes timing summary file exists in $CASEROOT/timing directory'
  exit
endif 

set dash = "-"
set ldir = `pwd`
setenv caseroot $ldir
setenv timeroot "unset"
setenv costscalenum 1
setenv costscaleden 1
set search1 = ""
set search2 = ""
set search3 = ""

while (1)
  if($#argv < 1) break
  set i = $argv[1]
  shift argv
  if($#argv < 1) then 
      echo "Error: wrong argument for $i."
      exit -1
  endif
  if( $argv[1] =~ $dash* ) then
          echo "Error: wrong argument for $i.";
          exit -1
  endif

  switch ( $i )
    case "-mach"
         setenv mach $argv[1]
         breaksw
    case "-lid"
         setenv lid $argv[1]
         breaksw
  default:
         echo "Unknown option."
         exit -1
         breaksw
  endsw
  shift argv
end  #---end of while loop.

source $caseroot/Tools/ccsm_getenv

set lookfor = `env | grep MAX_TASKS_PER_NODE | wc -l`
if !($lookfor == "0") then
  set lookfor = `env | grep PES_PER_NODE | wc -l`
  if !($lookfor == "0") then
    setenv costscalenum $PES_PER_NODE
    setenv costscaleden $MAX_TASKS_PER_NODE
  endif
endif

if ($timeroot == "unset") setenv timeroot $CASEROOT/Tools
set date = `date`
setenv date "$date"
set inittype = "FALSE"
if ($CONTINUE_RUN == "FALSE" && $RUN_TYPE == "startup") set inittype = "TRUE"
if ($CONTINUE_RUN == "FALSE" && $RUN_TYPE == "hybrid" ) set inittype = "TRUE"
setenv ccsmuser "unknown"
set lookfor = `env | grep CCSMUSER | wc -l`
if !($lookfor == "0") then
  setenv ccsmuser $CCSMUSER
endif

## --- get log/timing files ---
#foreach lpath ($RUNDIR $caseroot/timing $search1 $search2 $search3)
#  if (-f $lpath/ccsm_timing_all.$lid) then
#    cp $lpath/ccsm_timing_all.$lid $dout/
#    set fin = $dout/ccsm_timing_all.$lid
#  endif
#  if (-f $lpath/ccsm_timing_all.$lid.gz) then
#    cp $lpath/ccsm_timing_all.$lid.gz $dout/
#    gunzip $dout/*.gz
#    set fin = $dout/ccsm_timing_all.$lid
#  endif
#end

set fin  = $caseroot/timing/ccsm_timing_summary.$lid
set fout = $caseroot/timing/ccsm_timing.$CASE.${lid}
if !(-e $fin) then
  echo "$caseroot/timing/ccsm_timing_summary.$lid not found, exit"
  exit -1
endif
mv -f ${fout} ${fout}.old

$timeroot/getTiming.pl -fin $fin >! $fout
