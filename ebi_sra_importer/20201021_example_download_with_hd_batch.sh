#!/bin/bash -l

### TORQUE stuff here ####

### To send email when the job is completed:
#PBS -m ae
#PBS -M adswafford@eng.ucsd.edu

# Tell Sun grid engine that this is a Bash script
#PBS -S /bin/bash
# Write errors to this file - make sure the directory exists
#PBS -e /home/adswafford/cluster/logs

# Log output to this file - make sure the directory exists
#PBS -o /home/adswafford/cluster/logs

#PBS -l walltime=1:00:00
#PBS -l nodes=1:ppn=16

# maximum amount of memory used by any single process
#PBS -l mem=50gb
# ALTERNATIVELY: -l mem=32gb # amount of physical memory used by the job

### Run in the queue named
### #PBS -q med4gb

# name of the job
#PBS -N blood_topup

#set up as array job, use ${PBS_ARRAYID} as key in dictionary for lookup
#PBS -t 5-6

# BASH stuff here
### Switch to the working directory; by default TORQUE launches processes
### from your home directory.
cd $PBS_O_WORKDIR
echo Working directory is $PBS_O_WORKDIR

set -e
cpus=$PBS_NUM_PPN

export TMPDIR=/panfs/panfs1.ucsd.edu/panscratch/$USER
tmp=$(mktemp -d --tmpdir)
export TMPDIR=$tmp
trap "rm -r $tmp; unset TMPDIR" EXIT

source ~/.bashrc
declare -A projArray=([3]=PRJNA491863   [6]=PRJNA634938 [2]=PRJNA439269 [5]=PRJNA600846 [4]=PRJNA543497 [1]=PRJNA385815)
echo ${projArray[${PBS_ARRAYID}]}
proj_id=${projArray[${PBS_ARRAYID}]}

downloader_env=ebi_sra_importer
downloader=/projects/tcga-data/shogun/2020_blood_greg/EBI_SRA_Downloader.py
out=/projects/tcga-data/shogun/2020_blood_greg/$proj_id
force_yaml=/projects/tcga-data/shogun/2020_blood_greg/human_blood.yml

cd $tmp

#ugly code but should do the job..
conda activate $downloader_env
cd $tmp

python $downloader -project $proj_id -o $out -v -f -yaml $force_yaml -strat WGS RNA-Seq Other -hd -method minimap2 -qc --prep_max 250 -p $cpus

conda deactivate
