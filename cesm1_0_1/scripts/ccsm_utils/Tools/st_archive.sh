#!/bin/sh

#short-term archive script - move model output out of run directory
#to free disc space for next invocation of model run
#must be executed from run directory

#function dispose:
#moves output files to specified area of st archive and will
#process interim files along the way: 
#arg1 => interim files flag
#arg2 => destination
#remaining args => actual files to be processed 
dispose() {
    if [ "$1" == "ifiles_y" ] && [ "$DOUT_S_SAVE_INT_REST_FILES" != "TRUE" ]; then
	shift
	shift
	rm $*       2> /dev/null
    else
	shift
	dest=$1
	mkdir -p $dest
	shift
	mv $* $dest 2> /dev/null
    fi
}

echo ""
echo "st_archive.sh: start of short-term archiving"

#validate required env var settings
if [ -z "$DOUT_S_ROOT" ]; then
    echo "st_archive.sh: error, environment variable DOUT_S_ROOT is required "
    echo "               for root location of short-term archive"
    echo "st_archive.sh: exiting"
    exit 1
fi

sta="${DOUT_S_ROOT}"
mkdir -p ${sta} 2> /dev/null
if [ $? -ne 0 ]; then
    echo "st_archive.sh: error, unable to create short-term archive directory"
    echo "st_archive.sh: exiting"
    exit 1
fi

if [ -z "$DOUT_S_SAVE_INT_REST_FILES" ]; then
    echo "st_archive.sh: warning, environment variable DOUT_S_SAVE_INT_REST_FILES is not "
    echo "               set - using "FALSE" as default for saving interim restart files"
    export DOUT_S_SAVE_INT_REST_FILES=FALSE
fi

if [ "$DOUT_S_SAVE_INT_REST_FILES" == "FALSE" ]; then
    echo "st_archive.sh: restart files from end of run will be saved, "
    echo "               interim restart files will be deleted"
fi

#create directory for restart files
set ${CASE}.cpl.r.*
cplfile=`ls -rt $* 2> /dev/null | tail -1`
dname=`echo $cplfile | sed "s/\.nc//; s/^.*\.r\.//;"`
if [ -d ${sta}/rest/${dname} ]; then
    rm -rf ${sta}/rest/${dname}
fi
mkdir -p ${sta}/rest/${dname}
if [ $? -ne 0 ]; then
    echo "st_archive.sh: error, unable to create rest directory"
    echo "st_archive.sh: exiting"
    exit 1
fi

#populate temp directory with pointer files
set rpointer.*
if [ $# -le 0 ]; then
    echo "st_archive.sh: error, script should be invoked from run directory..."
    echo "               expecting restart pointer files of the form 'rpointer.<component>'"
    echo "               but did not find any: $*"
    echo "st_archive.sh: exiting"
    exit 1
fi
mv $* ${sta}/rest/${dname}

set atm.log.*;                                                                                                        dispose ifiles_n ${sta}/atm/logs $*
set lnd.log.*;                                                                                                        dispose ifiles_n ${sta}/lnd/logs $*
set ice.log.*;                                                                                                        dispose ifiles_n ${sta}/ice/logs $*
set ocn.log.*;                                                                                                        dispose ifiles_n ${sta}/ocn/logs $*
set glc.log.*;                                                                                                        dispose ifiles_n ${sta}/glc/logs $*
set cpl.log.*;                                                                                                        dispose ifiles_n ${sta}/cpl/logs $*
set ccsm*.log.*;                                                                                                      dispose ifiles_n ${sta}/cpl/logs $*

#still needs tweaking - remove assimilate_dir.* directories?, dart_log.nml's?
set assimilate_dir.*/dart_log.*;                                                                                      dispose ifiles_n ${sta}/dart/logs $*

set ${CASE}.cam?.r.*;        latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_y ${sta}/atm/rest $*
set ${CASE}.cam?.rs.*;       latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_y ${sta}/atm/rest $*
set ${CASE}.cam?.ra.*;       latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_y ${sta}/atm/rest $*
set ${CASE}.cam?.rh0.*;      latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_y ${sta}/atm/rest $*
set ${CASE}.cam?.rh1.*;      latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_y ${sta}/atm/rest $*
set ${CASE}.cam?.rh2.*;      latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_y ${sta}/atm/rest $*
set ${CASE}.cam?.rh3.*;      latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_y ${sta}/atm/rest $*
set ${CASE}.cam?.rh4.*;      latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_y ${sta}/atm/rest $*
set ${CASE}.cam?.rh5.*;      latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_y ${sta}/atm/rest $*
set ${CASE}.cam?.h0.*;       latest=`ls -rt $* 2> /dev/null | tail -1`; cp $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_n ${sta}/atm/hist $*
set ${CASE}.cam?.h1.*;       latest=`ls -rt $* 2> /dev/null | tail -1`; cp $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_n ${sta}/atm/hist $*
set ${CASE}.cam?.h2.*;       latest=`ls -rt $* 2> /dev/null | tail -1`; cp $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_n ${sta}/atm/hist $*
set ${CASE}.cam?.h3.*;       latest=`ls -rt $* 2> /dev/null | tail -1`; cp $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_n ${sta}/atm/hist $*
set ${CASE}.cam?.h4.*;       latest=`ls -rt $* 2> /dev/null | tail -1`; cp $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_n ${sta}/atm/hist $*
set ${CASE}.cam?.h5.*;       latest=`ls -rt $* 2> /dev/null | tail -1`; cp $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_n ${sta}/atm/hist $*
set ${CASE}.cam?.i.*;        latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_y ${sta}/atm/init $*

set ${CASE}.clm?.r.*;        latest=`ls -rt $* 2> /dev/null | tail -2`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_y ${sta}/lnd/rest $*
set ${CASE}.clm?.h0.*;       latest=`ls -rt $* 2> /dev/null | tail -1`; cp $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_n ${sta}/lnd/hist $*
set ${CASE}.clm?.h1.*;       latest=`ls -rt $* 2> /dev/null | tail -1`; cp $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_n ${sta}/lnd/hist $*
set ${CASE}.clm?.h2.*;       latest=`ls -rt $* 2> /dev/null | tail -1`; cp $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_n ${sta}/lnd/hist $*
set ${CASE}.clm?.h3.*;       latest=`ls -rt $* 2> /dev/null | tail -1`; cp $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_n ${sta}/lnd/hist $*
set ${CASE}.clm?.h4.*;       latest=`ls -rt $* 2> /dev/null | tail -1`; cp $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_n ${sta}/lnd/hist $*
set ${CASE}.clm?.h5.*;       latest=`ls -rt $* 2> /dev/null | tail -1`; cp $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_n ${sta}/lnd/hist $*
set ${CASE}.clm?.hv.*;       latest=`ls -rt $* 2> /dev/null | tail -1`; cp $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_n ${sta}/lnd/hist $*
set ${CASE}.clm?.i.*;                                                                                                 dispose ifiles_y ${sta}/lnd/init $*

set ${CASE}.cice.r.[0-9]*;   latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_y ${sta}/ice/rest $*
set ${CASE}.cice.r.volpn*;   latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_y ${sta}/ice/rest $*
set ${CASE}.cice.r.dEdd*;    latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_y ${sta}/ice/rest $*
set ${CASE}.cice.r.age*;     latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_y ${sta}/ice/rest $*
set ${CASE}.cice.r.aero*;    latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_y ${sta}/ice/rest $*
set ${CASE}.cice.h*;                                                                                                  dispose ifiles_n ${sta}/ice/hist $*
set ${CASE}.cice.i.*;                                                                                                 dispose ifiles_y ${sta}/ice/init $*

set ${CASE}.cism.r.[0-9]*;   latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_y ${sta}/glc/rest $*
set ${CASE}.cism.r.volpn*;   latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_y ${sta}/glc/rest $*
set ${CASE}.cism.r.dEdd*;    latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_y ${sta}/glc/rest $*
set ${CASE}.cism.r.age*;     latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_y ${sta}/glc/rest $*
set ${CASE}.cism.r.aero*;    latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_y ${sta}/glc/rest $*
set ${CASE}.cism.h*;                                                                                                  dispose ifiles_n ${sta}/glc/hist $*
set ${CASE}.cism.i.*;                                                                                                 dispose ifiles_y ${sta}/glc/init $*

set ${CASE}.pop.r.*.hdr;     latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_y ${sta}/ocn/rest $*
set ${CASE}.pop.r.*0000;     latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_y ${sta}/ocn/rest $*
set ${CASE}.pop.r.*0000.nc;  latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_y ${sta}/ocn/rest $*
set ${CASE}.pop.rh.ecosys.*.hdr;    latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_y ${sta}/ocn/rest $*
set ${CASE}.pop.rh.ecosys.*0000;    latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_y ${sta}/ocn/rest $*
set ${CASE}.pop.rh.ecosys.*0000.nc; latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_y ${sta}/ocn/rest $*
set ${CASE}.pop.rh.*.hdr;    latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_y ${sta}/ocn/rest $*
set ${CASE}.pop.rh.*0000;    latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_y ${sta}/ocn/rest $*
set ${CASE}.pop.rh.*0000.nc; latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_y ${sta}/ocn/rest $*
set ${CASE}.pop.ro.*;        latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_y ${sta}/ocn/rest $*
set ${CASE}.pop.d?*;                                                                                                  dispose ifiles_n ${sta}/ocn/hist $*
set ${CASE}.pop.h*;                                                                                                   dispose ifiles_n ${sta}/ocn/hist $*

set ${CASE}.camice.r.*;      latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_y ${sta}/ice/rest $*
set ${CASE}.camdom.r.*;      latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_y ${sta}/ocn/rest $*
set ${CASE}.camsom.r.*;      latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_y ${sta}/ocn/rest $*

set ${CASE}.datm.r.* ;      latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_y ${sta}/atm/rest $*
set ${CASE}.datm.rs* ;      latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_y ${sta}/atm/rest $*
set ${CASE}.datm.h.* ;                                                                                               dispose ifiles_n ${sta}/atm/hist $*

set ${CASE}.dlnd.r.* ;      latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_y ${sta}/lnd/rest $*
set ${CASE}.dlnd.rs* ;      latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_y ${sta}/lnd/rest $*
set ${CASE}.dlnd.h.* ;                                                                                               dispose ifiles_n ${sta}/lnd/hist $*

set ${CASE}.docn.r.* ;      latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_y ${sta}/ocn/rest $*
set ${CASE}.docn.rs* ;      latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_y ${sta}/ocn/rest $*
set ${CASE}.docn.h.* ;                                                                                               dispose ifiles_n ${sta}/ocn/hist $*

set ${CASE}.dice.r.* ;      latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_y ${sta}/ice/rest $*
set ${CASE}.dice.rs* ;      latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_y ${sta}/ice/rest $*
set ${CASE}.dice.h.* ;                                                                                               dispose ifiles_n ${sta}/ice/hist $*

set ${CASE}.cpl.r.*;         latest=`ls -rt $* 2> /dev/null | tail -1`; mv $latest ${sta}/rest/${dname} 2> /dev/null; dispose ifiles_y ${sta}/cpl/rest $*
set ${CASE}.cpl.h* ;                                                                                                  dispose ifiles_n ${sta}/cpl/hist $*

set Prior_Diag.*.nc;                                                                                                  dispose ifiles_n ${sta}/dart/hist $*
set Posterior_Diag.*.nc;                                                                                              dispose ifiles_n ${sta}/dart/hist $*
set obs_seq.*.final;                                                                                                  dispose ifiles_n ${sta}/dart/hist $*

#copy back the required files for next restart
cp ${sta}/rest/${dname}/* .

echo "st_archive.sh: short-term archiving completed successfully"

exit 0
