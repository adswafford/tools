#!/bin/bash -l

### TORQUE stuff here ####

### To send email when the job is completed:
#PBS -m ae
#PBS -M [your email]

# Tell Sun grid engine that this is a Bash script
#PBS -S /bin/bash
# Write errors to this file - make sure the directory exists
#PBS -e [error log path]

# Log output to this file - make sure the directory exists
#PBS -o [output log path]

#PBS -l walltime=HH:MM:SS
#PBS -l nodes=1:ppn=##

# maximum amount of memory used by any single process
#PBS -l mem=###gb

# name of the job
#PBS -N [name]_array

#set up as array job, use ${PBS_ARRAYID} as key in dictionary for lookup
#PBS -t [start num]-[end num]%[num jobs at once]

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

studies=([1]=[study ID1] [2]=[study ID2]) #add as many entries as needed
out_dir=/path/to/your/output/folder
EBI_script=/path/to/your/EBI_SRA_Downloader.py
mode=-ebi

cd $tmp

#first download the study
conda activate ebi_sra_importer

echo ${studies[${PBS_ARRAYID}]}
study=${studies[${PBS_ARRAYID}]}

cd $tmp
mkdir $tmp/$study

cd $tmp/$study
python $EBI_script $mode $study -all-seqs -all-platforms

mv $tmp/$study $out_dir/$study
conda deactivate
