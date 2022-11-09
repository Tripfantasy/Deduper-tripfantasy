#!/bin/bash
#SBATCH --account=bgmp                  ### SLURM account which will be charged for the job
#SBATCH --job-name=Deduper              ### Job Name
#SBATCH --output=Deduper_%j.out         ### File in which to store job output
#SBATCH --error=Deduper-%j.err          ### File in which to store job error messages
#SBATCH --time=0-24:00:00               ### Wall clock time limit in Days-HH:MM:SS
#SBATCH --nodes=1                       ### Node count required for the job
#SBATCH --cpus-per-task=4               ### Number of cpus (cores) per task
#SBATCH --partition=bgmp                ### partition to run things
#SBATCH --mem=25G
samfile="/projects/bgmp/dmarro/bioinfo/Bi624/Deduper/Deduper-tripfantasy/C1_SE_uniqAlign_sorted.sam"
umifile="/projects/bgmp/dmarro/bioinfo/Bi624/Deduper/Deduper-tripfantasy/STL96.txt"
outsam="/projects/bgmp/dmarro/bioinfo/Bi624/Deduper/Deduper-tripfantasy/Deduped_C1_SE_uniqAlign_sorted.sam"

conda activate bgmp_py310
/usr/bin/time -v ./marro_deduper.py -u $umifile -f $samfile -o $outsam