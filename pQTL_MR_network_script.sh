#!/bin/bash
#
#
#PBS -l nodes=1:ppn=1,walltime=12:00:00


WORK_DIR="/mnt/storage/home/hw15842"
module add languages/r/3.6.0

R CMD BATCH --no-save --no-restore '--args outcome_protein ' script.R script.out 


make sure to do sleep command for a job array



#!/bin/bash
#PBS -N my_job_array
#PBS -l nodes=1:ppn=1,walltime=12:00:00
#PBS -t 1-3282%10
cd $PBS_O_WORKDIR
echo "The Array ID is: ${PBS_ARRAYID}"
let START_VAL=${PBS_ARRAYID}*100
echo "START_VAL is set to: ${START_VAL}"
## and perhaps pass this value to your executable as a command-line argument:
## ./my_prog.exe ${START_VAL}