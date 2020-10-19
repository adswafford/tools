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

#PBS -l walltime=24:00:00
#PBS -l nodes=1:ppn=16

# maximum amount of memory used by any single process
#PBS -l mem=50gb
# ALTERNATIVELY: -l mem=32gb # amount of physical memory used by the job

### Run in the queue named
### #PBS -q med4gb

# name of the job
#PBS -N ebi_array

#set up as array job, use ${PBS_ARRAYID} as key in dictionary for lookup
#PBS -t 1-28%5

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

declare -A projArray=([1]=PRJDB6396 [2]=PRJEB13247  [3]=PRJEB30048  [4]=PRJEB30449  [5]=PRJEB30958  [6]=PRJNA257350 [7]=PRJNA306662 [8]=PRJNA328281 [9]=PRJNA385009 [10]=PRJNA385180    [11]=PRJNA414721    [12]=PRJNA471637    [13]=PRJNA484074    [14]=PRJNA507824    [15]=PRJNA513060    [16]=PRJNA517159    [17]=PRJNA544518    [18]=PRJNA554271    [19]=PRJNA562379    [20]=PRJNA578933    [21]=PRJNA596372    [22]=PRJNA601241    [23]=PRJNA605079    [24]=PRJNA607097    [25]=PRJNA632511    [26]=PRJNA633741    [27]=PRJNA639958    [28]=PRJNA642002)
echo ${projArray[${PBS_ARRAYID}]}
proj_id=${projArray[${PBS_ARRAYID}]}

downloader_env=ebi_sra_importer
downloader=/projects/cmi_proj/blood_microbiome/public_cfdna/EBI_SRA_Downloader.py
out=/projects/cmi_proj/blood_microbiome/public_cfdna/$proj_id
force_yaml=/projects/cmi_proj/blood_microbiome/public_cfdna/human_plasma.yml

cd $tmp

#ugly code but should do the job..
conda activate $downloader_env
cd $tmp

python $downloader -project $proj_id -o $out -v -f -hd -yaml $force_yaml -strat WGS WXS AMPLICON RNA-Seq Other

conda deactivate
