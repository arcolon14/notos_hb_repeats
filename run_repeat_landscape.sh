#!/bin/bash
#SBATCH -p queue
#SBATCH -J jobid
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -t 4:00:00

module load gcc/7.2.0
module load python/3
thr=${SLURM_NTASKS}

src=/path/to/repeatmasker-4.1.2-p1/RepeatMasker/util
work=/path/to/working/dir

# Species
spp=chaEso

# Repeat Masker directory
rep_mask=$work/${spp}/repeat_masker.out
align_f=$rep_mask/${spp}.fa.align
sum_f=$rep_mask/${spp}.divsum

cmd=(
    $src/calcDivergenceFromAlign.pl
    -s $sum_f
    $align_f
)
echo "${cmd[@]}"
"${cmd[@]}"

# Run the summary table
src=${work}/scripts
cross_match=$rep_mask/${spp}.fa.out
cmd=(
    $src/parse_repeat_masker_out.py
    --cross-match $cross_match
    --divsum $sum_f
    --outdir $rep_mask
    --min-length 10
)
echo "${cmd[@]}"
"${cmd[@]}"

# Tabulate the DIVSUM file
cat $sum_f | \
    tail -n +5 | \
    grep -v '\--' | \
    awk 'NF==5 {print $0}' | \
    tr -d '%' > ${sum_f}.tsv

# compress the alignment file
gzip $align_f
