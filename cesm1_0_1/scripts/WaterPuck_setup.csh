#! /bin/csh -f
#foreach LOC (SF GT GD GB BM)
#foreach MM (01 02 03 04 05 06 07 08 09 10 11 12)
#    set wdir=CH_$LOC$MM
#    echo $wdir
#end
#end 
setenv my_dir $PWD
setenv wdir /scratch/lustre/plgjjakacki/LD/cases/WP05X

#env_run.xml
setenv CONTINUE_RUN FALSE
setenv INFO_DBUG 1
setenv HISTINIT TRUE
cd $wdir
./xmlchange -file env_run.xml -id CONTINUE_RUN -val $CONTINUE_RUN
./xmlchange -file env_run.xml -id INFO_DBUG -val $INFO_DBUG
./xmlchange -file env_run.xml -id HISTINIT -val $HISTINIT
#env_conf.xml
#setenv RUN_STARTDATE "2010-01-02"
setenv RUN_STARTDATE "2010-07-02"
setenv RUN_REFDATE $RUN_STARTDATE
./xmlchange -file env_conf.xml -id RUN_STARTDATE -val $RUN_STARTDATE
./xmlchange -file env_conf.xml -id RUN_REFDATE -val $RUN_REFDATE
#env_mach_pes.xml
setenv ntasks 160
setenv ntasks_ice 1
./xmlchange -file env_mach_pes.xml -id NTASKS_ATM -val $ntasks
./xmlchange -file env_mach_pes.xml -id NTASKS_LND -val $ntasks
./xmlchange -file env_mach_pes.xml -id NTASKS_ICE -val $ntasks_ice
./xmlchange -file env_mach_pes.xml -id NTASKS_OCN -val $ntasks
./xmlchange -file env_mach_pes.xml -id NTASKS_CPL -val $ntasks
./xmlchange -file env_mach_pes.xml -id NTASKS_GLC -val $ntasks
./xmlchange -file env_mach_pes.xml -id TOTALPES  -val $ntasks

# copy cice_decomp.xml - i am not sure that it is needed, but it is requested
cp /scratch/lustre/plgjjakacki/LD/FindFish_Model/cesm1_0_1/models/ice/cice/bld/cice_decomp.xml $wdir/Tools/Templates
#"${n_filezip_name}" "${Dest}/$1_filezip_name"
cd $my_dir
