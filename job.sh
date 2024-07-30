#!/bin/bash -l
#$ -l h_rt=48:00:00
#$ -l mem=75G
#$ -pe mpi 36
#$ -N psam_10_10
#$ -A KCL_Lorenz
#$ -P Gold
#$ -ac allow=Z
#$ -cwd 

module purge
module load compilers/gnu/4.9.2
module load orca/5.0.4-sbindist
module load gerun

export OMP_NUM_THREADS=1
ORCA_EXEC=$(which orca)
$ORCA_EXEC input.inp > output.out
