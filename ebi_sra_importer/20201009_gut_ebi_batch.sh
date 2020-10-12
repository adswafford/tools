#!/bin/bash -l

### TORQUE stuff here ####

### To send email when the job is completed:
#PBS -m ae
#PBS -M adswafford@eng.ucsd.edu

# Tell Sun grid engine that this is a Bash script
#PBS -S /bin/bash
# Write errors to this file - make sure the directory exists
#PBS -e /projects/cmi_proj/seed_grants/sd_zoo/cheetahs/gotu_errors

# Log output to this file - make sure the directory exists
#PBS -o /projects/cmi_proj/seed_grants/sd_zoo/cheetahs/gotu_output

#PBS -l walltime=6:00:00
#PBS -l nodes=1:ppn=8

# maximum amount of memory used by any single process
#PBS -l mem=50gb
# ALTERNATIVELY: -l mem=32gb # amount of physical memory used by the job

### Run in the queue named
### #PBS -q med4gb

# name of the job
#PBS -N ebi_array

#set up as array job, use ${PBS_ARRAYID} as key in dictionary for lookup
#PBS -t 2-5

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

declare -A projArray=([1]=PRJEB21528   [2]=PRJEB7774   [3]=PRJEB12449  [4]=PRJEB1220 [5]=PRJNA373879)

echo ${projArray[${PBS_ARRAYID}]}
proj_id=${projArray[${PBS_ARRAYID}]}

downloader_env=ebi_sra_importer
downloader=/projects/covid/EBI_auto_test/EBI_SRA_Downloader.py
out=/projects/covid/EBI_auto_test/$proj_id

cd $tmp

#ugly code but should do the job..
conda activate $downloader_env

python $downloader -project $proj_id -o $out --prep_max 250 -v

conda deactivate
