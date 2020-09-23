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
#PBS -t 3-34

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

declare -A projArray=([1]=PRJEB21528   [2]=PRJEB7774   [3]=PRJEB6070   [4]=PRJEB12449  [5]=PRJEB10878  [6]=PRJEB1220   [7]=PRJNA385949 [8]=PRJNA389280 [9]=PRJEB15371  [10]=PRJNA278393    [11]=PRJNA422434    [12]=PRJEB1786  [13]=PRJEB6337  [14]=PRJEB12123 [15]=PRJEB19090 [16]=PRJNA305507    [17]=PRJNA268964    [18]=PRJEB1690  [19]=PRJNA319574    [20]=PRJEB11532 [21]=PRJEB6997  [22]=PRJNA299502    [23]=PRJDB3601  [24]=PRJEB8094  [25]=PRJEB13870 [26]=PRJEB6456  [27]=PRJEB4336  [28]=PRJNA177201    [29]=PRJNA373879    [30]=PRJNA328899 [31]=PRJNA48479 [32]=PRJNA290729 [33]=PRJEB12947 [34]=PRJNA598446)

echo ${projArray[${PBS_ARRAYID}]}
proj_id=${projArray[${PBS_ARRAYID}]}

downloader_env=ebi_sra_importer
downloader=/projects/covid/tmp_ghmi/EBI_SRA_Downloader.py
out=/projects/covid/tmp_ghmi/$proj_id

cd $tmp

#ugly code but should do the job..
conda activate $downloader_env
cd $tmp
mkdir $tmp/$proj_id
cd $tmp/$proj_id
python $downloader -project $proj_id -v -f -yaml /projects/covid/tmp_ghmi/human_feces.yml --no_seqs
mv $tmp/$proj_id $out

conda deactivate
