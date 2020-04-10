#!/bin/bash
#SBATCH -J gyre_io_01
####SBATCH -p batch
#SBATCH -N 1
#SBATCH --ntasks-per-node 24
#SBATCH --mem 10000
#SBATCH --time 24:00:00
##SBATCH --mail-type=END
##SBATCH --mail-user=user5@email.domena.pl



./get_data > get_data.out

