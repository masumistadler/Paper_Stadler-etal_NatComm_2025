#!/bin/bash
#SBATCH --account=def-pauldel
#SBATCH --time=5-00:00
#SBATCH --mem=100G #10G
#SBATCH --cpus-per-task=32 #48, max 48
#SBATCH --ntasks=1
#SBATCH --job-name="MO_modelling"
#SBATCH --mail-user=m.stadler.jp.at@gmail.com
#SBATCH --mail-type=ALL
#---------------------------------------------

#load the executables
module load StdEnv/2020 gcc/9.3.0 r/4.2.2 openmpi/4.0.3

#command to execute
Rscript ~/projects/def-pauldel/mstadler/Jobs/5.2_MO_modelling.R
#---------------------------------------------
