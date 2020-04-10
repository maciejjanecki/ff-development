#!/bin/bash
#SBATCH -J TS
####SBATCH -p batch
#SBATCH -N 1
#SBATCH --ntasks-per-node 24
#SBATCH --mem 10000
#SBATCH --time 24:00:00
##SBATCH --mail-type=END
#SBATCH --mail-user=jjakacki@iopan.pl



./TS_boundaries_data > TS_boundaries_data.out

