#!/bin/bash
#SBATCH --account=def-pauldel
#SBATCH --time=5-00:00
#SBATCH --mem=30G #10G
#SBATCH --cpus-per-task=24 #48, max 48
#SBATCH --ntasks=1
#SBATCH --job-name="trav.time_15-16"
#SBATCH --mail-user=m.stadler.jp.at@gmail.com
#SBATCH --mail-type=ALL
#---------------------------------------------

#load the executables
module load StdEnv/2020 gcc/9.3.0 r/4.2.2 openmpi/4.0.3

#command to execute
Rscript ~/projects/def-pauldel/mstadler/Jobs/3_travel.time_2015-2016.R
#---------------------------------------------
