#!/bin/bash
#SBATCH -p queue
#SBATCH -J jobname
#SBATCH -N 1
#SBATCH -n 24
#SBATCH -t 168:00:00

module load gcc/7.2.0
module load python/3
thr=${SLURM_NTASKS}
npr=$(($thr/4))

# Run RepeatMasker on a specific target sequence using the
# icefish repeat library from Rivera-Colon et al. 2023 MBE

bin=/path/to/repeatmasker-4.1.2-p1/RepeatMasker
work=/path/to/working/dir

# Define target species
spp=chaEso

# Set the data directories and input FASTA
sp_dir=${work}/species_database/${spp}
in_fasta=${sp_dir}/genome/${spp}.fa

# Repeat library
# ARC MBE 2023
rep_db=$work/rep_database/arc_2023_mbe_db.fa

# Make output directory for the target run
outp=${sp_dir}/repeat_masker.out/
mkdir -p $outp
cd $outp

# Repeat Masker command
cmd=(
    ${bin}/RepeatMasker
    -lib $rep_db
    -pa $npr
    -a
    -dir $outp
    -xsmall
    -gff
    $in_fasta
)

echo "${cmd[@]}"
"${cmd[@]}"
