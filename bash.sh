#!/bin/sh
### Every line in this header section should start with ### for a comment
### or #PBS for an option for qsub
### Note: No unix commands may be executed until after the last #PBS line
###
### Account information
#PBS -W group_list=srhgroup -A srhgroup
##
### Send mail when job is aborted or terminates abnormally
#PBS -M s174596@student.dtu.dk
#PBS -m abe
###
### Compute resources, here 1 core on 1 node
#PBS -l nodes=1:ppn=1:thinnode
###
### Required RAM in GB
#PBS -l mem=100GB
###
### How long (max) will the job take, here 24h:
#PBS -l walltime=24:00:00
###
### Output files - not required to be specified
### Comment out the next 2 lines to use the job id instead in the file names
#PBS -e /Volumes/T-cells-and-cancer/SRH\ group/Group\ members/Simone/results
#PBS -o /Volumes/T-cells-and-cancer/SRH\ group/Group\ members/Simone/results
###
### Job name - not required to be specified
### It is often easier just to use the job id instead for recognition
### PBS -N Simone_CellRanger
###
### More qsub options can be added here
#PBS -d /home/projects/SRHgroup/projects/SingleCell_KIR/......



# This part is the real job script
# Here follows the user commands:
cd /Volumes/T-cells-and-cancer/SRH group/Group members/Simone/data

# Load all required modules for the job
module load tools cellranger/7.0.0

# Change directory to the result paths
cd /Volumes/T-cells-and-cancer/SRH\ group/Group\ members/Simone/results

...../cellranger-7.0.0/cellranger multi --id=results --csv=/Volumes/T-cells-and-cancer/SRH\ group/Group\ members/Simone/scripts/config.csv

