#!/bin/sh 
#

if [ $# -ne 1 ]; then
    echo "TCB.sh: incorrect number of input arguments" 
    exit 1
fi

test_name=TCB.$1

if [ -f ${CAM_TESTDIR}/${test_name}/TestStatus ]; then
    if grep -c PASS ${CAM_TESTDIR}/${test_name}/TestStatus > /dev/null; then
        echo "TCB.sh: configure and build test has already passed; results are in "
	echo "        ${CAM_TESTDIR}/${test_name}" 
        exit 0
    else
	read fail_msg < ${CAM_TESTDIR}/${test_name}/TestStatus
        prev_jobid=${fail_msg#*job}

	if [ $JOBID = $prev_jobid ]; then
            echo "TCB.sh: configure and build test has already failed for this job - will not reattempt; "
	    echo "        results are in: ${CAM_TESTDIR}/${test_name}" 
	    exit 2
	else
	    echo "TCB.sh: this configure and build test failed under job ${prev_jobid} - moving those results to "
	    echo "        ${CAM_TESTDIR}/${test_name}_FAIL.job$prev_jobid and trying again"
            cp -rp ${CAM_TESTDIR}/${test_name} ${CAM_TESTDIR}/${test_name}_FAIL.job$prev_jobid
        fi
    fi
fi

cfgdir=${CAM_SCRIPTDIR}/../../bld
blddir=${CAM_TESTDIR}/${test_name}
if [ -d ${blddir} ]; then
    rm -r ${blddir}
fi
mkdir -p ${blddir} 
if [ $? -ne 0 ]; then
    echo "TCB.sh: error, unable to create work subdirectory" 
    exit 3
fi
cd ${blddir}

if [ ! -f ${CAM_SCRIPTDIR}/config_files/$1 ]; then
    echo "TCB.sh: configure options file ${CAM_SCRIPTDIR}/config_files/$1 not found" 
    echo "FAIL.job${JOBID}" > TestStatus
    exit 4
fi

## task and thread_flags needed so scam test (nospmd, nosmp) can get proper cice decomposition 
if grep -c nospmd ${CAM_SCRIPTDIR}/config_files/$1 > /dev/null; then
    task_flag=1
    thrd_flag=1
else
    if grep -c nosmp ${CAM_SCRIPTDIR}/config_files/$1 > /dev/null; then
	task_flag=`expr $CAM_TASKS "*" $CAM_THREADS`
        thrd_flag=1
    else
        task_flag=$CAM_TASKS
        thrd_flag=$CAM_THREADS
    fi
fi


##construct string of args to configure
config_string=$CFG_STRING 
while read config_arg; do
    config_string="${config_string}${config_arg} "
done < ${CAM_SCRIPTDIR}/config_files/$1

##determine cice decomposition
if grep -w csim4 ${CAM_SCRIPTDIR}/config_files/$1; then
    echo "TCB.sh: using csim4 ice model"
elif grep -w sice ${CAM_SCRIPTDIR}/config_files/$1; then
    echo "TCB.sh: using sice ice model"
else
    res_string=${config_string#*-res }
    res_string=${res_string%% *}
    cicebld="${CAM_SCRIPTDIR}/../../../../ice/cice/bld"
    cice_config=`$cicebld/generate_cice_decomp.pl -res $res_string -nproc $task_flag -thrds $thrd_flag -output all`
    rc=${cice_config%% *}
    if [ $rc -ge 0 ]; then
	echo "TCB.sh: generate_cice_decomp was successful" 
        echo "TCB.sh: cice_config= $cice_config"
    else
	echo "TCB.sh: generate_cice_decomp.pl failed, error = $cice_config" 
	echo "FAIL.job${JOBID}" > TestStatus
	exit 7
    fi
    cice_decomp=${cice_config##* }
    cice_config=${cice_config% *}
    cice_maxblocks=${cice_config##* }
    cice_config=${cice_config% *}
    cice_bsizey=${cice_config##* }
    cice_config=${cice_config% *}
    cice_bsizex=${cice_config##* }
    cice_config=${cice_config% *}
    if [ $cice_decomp = "cartesian" ]; then
	cice_maxblocks=`expr $cice_maxblocks "*" 2`
    else
	echo "test scripts only support cartesian cice decomposition" > TestStatus
	exit 5
    fi
    config_string="$config_string -cice_bsizex $cice_bsizex -cice_bsizey $cice_bsizey -cice_maxblocks $cice_maxblocks -cice_decomptype $cice_decomp" 
cat >> cice_namelist << EOF
 &domain_nml
 distribution_type = '$cice_decomp'
 / 
EOF
fi

echo "TCB.sh: building cam executable; output in ${CAM_TESTDIR}/${test_name}/test.log" 

attempt=1
still_compiling="TRUE"
while [ $still_compiling = "TRUE" ]; do

    echo "TCB.sh: call to configure:" 
    echo "        ${cfgdir}/configure ${config_string}" 

    ${cfgdir}/configure ${config_string} > test.log 2>&1
    rc=$?
    if [ $rc -eq 0 ]; then
	echo "TCB.sh: configure was successful" 
    else
	echo "TCB.sh: cam configure failed, error from configure= $rc" 
	echo "TCB.sh: see ${CAM_TESTDIR}/${test_name}/test.log for details"
	echo "FAIL.job${JOBID}" > TestStatus
	exit 5
    fi

    echo "TCB.sh: call to make:" 
    echo "        ${MAKE_CMD}" 
    ${MAKE_CMD} >> test.log 2>&1
    rc=$?
    if [ $rc -eq 0 ]; then
	echo "TCB.sh: make was successful" 
	echo "TCB.sh: configure and build test passed"
	echo "PASS" > TestStatus
	if [ $CAM_RETAIN_FILES != "TRUE" ]; then
	    echo "TCB.sh: removing some unneeded files to save disc space" 
	    rm *.o
	    rm *.mod
	fi
	still_compiling="FALSE"
    elif [ $attempt -lt 10 ] && \
        grep -c "LICENSE MANAGER PROBLEM" test.log > /dev/null; then
        attempt=`expr $attempt + 1`
        echo "TCB.sh: encountered License Manager Problem; launching attempt #$attempt"
    else
	echo "TCB.sh: cam build failed, error from make= $rc" 
	echo "TCB.sh: see ${CAM_TESTDIR}/${test_name}/test.log for details"
	echo "FAIL.job${JOBID}" > TestStatus
	exit 6
    fi
done

exit 0
