#! /usr/bin/csh -f

# -------------------------------------------------------------------------
# CCSM automatic job submitter for Earth Simulator
# -------------------------------------------------------------------------
#
# Usage: ccsm_auto.csh -interval seconds -njob #_of_jobs
#
# Note:
#
#  0) ES batch job nodes do not have NQS-II system and batch jobs can
#     be submitted only from the frontend server, "moon". 
#  1) This script submits "#_of_jobs" batch jobs. Each job is
#     submitted after the completion of its preceding job.
#  2) To do the above, this script refers $CNTLDIR/CURRENT_JOB and
#     $CNTLDIR/status.$LID, the former of which contains $LID for the
#     # last job, and the latter contains its status. Currently, job
#     status # is classified into BUILDING, SUBMITTED, COMPLETED and
#     FAILED.
#  3) If the last job is COMPLETED, then the script will submit its
#     successor job.
#  4) $CNTLDIR/CURRENT_JOB and $CNTLDIR/status.$LID are updated by the
#     $CASE.moon.run script.
#  5) A default of environment variable CNTLDIR is $cwd/cntl. Or you
#     can setenv CNTLDIR.
#  6) This script must be run on $CASEROOT, and is valid only when 
#     $CONTINUE_RUN is TRUE.
#
# -------------------------------------------------------------------------

./Tools/ccsm_check_lockedfiles || exit -1

source ./Tools/ccsm_getenv || exit -1

if ($CONTINUE_RUN != TRUE) then
  echo "You can use $0 only when CONTINUE_RUN is TRUE"
  exit -1
endif

if ($#argv != 4) then
  echo "usage: $0 -interval seconds -njob #_of_jobs"
  exit -1
endif

while (1)
  if ( $#argv < 1 ) break
  set i = $argv[1];
  switch ( $i )
    case -interval:
      shift argv
      set interval = $argv[1]
      shift argv
      breaksw 
    case -njob:
      shift argv
      set njob = $argv[1]
      shift argv
      breaksw
    default:
      echo "usage: $0 -interval seconds -njob #_of_jobs"
      exit -1
      breaksw
  endsw
end

if (! $?CNTLDIR) setenv CNTLDIR ./cntl

@ n = 0
while ($n < $njob)
  set waiting = 1
  while (1)
    date
    echo 'checking status...'
    if (! -f $CNTLDIR/CURRENT_JOB) then
      date
      echo 'no preceding job'
      break
    endif
    set reqid  = `cat $CNTLDIR/CURRENT_JOB`
    set status = `cat $CNTLDIR/status.$reqid`
    if ($status == 'COMPLETED') break
    if ($status == 'FAILED') then
      date
      echo 'FATAL ERROR: previous job was failed'
      exit -1
    endif
    echo 'sleeping...'
    sleep $interval
  end
  date
  echo 'staring build-and-run script...'
  ./$CASE.moon.run || (echo "$CASE.moon.run failed"; exit -1)
  date
  echo 'build-and-run script completed...'
  sleep 30
  @ n++
end
date
echo "$0 completed."
