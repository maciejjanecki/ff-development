#!/bin/bash

### please set number of months (12 o 10)
nmonths=6


my_pwd=$PWD
cyy=2013
syy=0001
cmm=01
cdd=01

pathIn=/scratch/lustre/plgjjakacki/LD/tmp_data/ARTUR/ICM/ICM_${cyy}_115m
pathOut=/scratch/lustre/plgjjakacki/LD/cesm_input_data/atm/datm7/bs01v1/data_v2
#pathOut=/scratch/lustre/plgjjakacki/LD/cesm_input_data/atm/datm7/bs01v1/data
echo $pathIn
echo $pathOut

#exit

mc=( 31 28 31 30 31 30 31 31 30 31 30 31 )

cd $pathOut 
while [[ $m -lt $nmonths ]]; do #### 12 ]]; do
  md=${mc[m]}
  let m=m+1
  d=1
  while [[ $d -le $md ]]; do
    mm=$m
    if [[ $m -le 9 ]]
         then mm=0${m}
    fi
    dd=$d
    if [[ $d -le 9 ]]
         then dd=0${d}
    fi  
    fin=${pathIn}/${cyy}${mm}${dd}_115m.nc
    fout=${syy}${mm}${dd}_115m.nc
    rm $fout
    ln -s $fin $fout    
    echo $fin $fout
    let d=d+1
  done
done 
pwd
cd $my_pwd


#for MC in {1..12}; 
#   do 
#for i in {1..10}; do echo $i; done
#     for D in {1..$DIM{$MC}}
#       do
#         echo $MC $D
#       done
#   done

#exit

##for LOC in SF GT GD GB BM
##  do
#    for MC in 01 02 03 04 05 06 07 08 09 10 11
#      do
##         cd $my_pwd
#         export EX=CH_"$LOC""$MC"
#         echo $EX
#         echo 'input filename'
#         infile=$my_pwd/$EX/$path_in/$EX.$fname.$cyy-$MC-01.nc
#         echo $infile
#         echo 'output filename'
#         export outfile=$out_dir/$LOC/$EX.$fname.$cyy-$MC-01.nc
#         echo $outfile
#         ln -s $infile $outfile
#      done
##  done

#echo $my_pwd

