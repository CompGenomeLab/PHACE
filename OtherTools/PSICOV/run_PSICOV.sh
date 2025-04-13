#!/bin/bash
#SBATCH --job-name=PSICOV
#SBATCH --account=investor
#SBATCH --nodes=1
#SBATCH --array=1-652
# #SBATCH --mem=120G
#SBATCH --cpus-per-task=8
#SBATCH --qos=long_investor
#SBATCH --partition=genomics
#SBATCH --time=7-00:00:00
#SBATCH --output=./output/%A_%a-slurm.out
#SBATCH --mail-type=ALL

i=$SLURM_ARRAY_TASK_ID
echo $i


declare -a FILELIST=(PROTEIN_IDs_LIST)

echo ${FILELIST[i]}
id=${FILELIST[i]}

./psicov -o "MSAs/${id}.fasta" -j 0 -z 8  > "/cta/users/nkuru/psicov/Results/Option_o/${id}.txt"

