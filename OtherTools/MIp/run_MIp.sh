#!/bin/bash
#SBATCH --job-name=MIp
#SBATCH --account=investor
#SBATCH --nodes=1
#SBATCH --array=1-652
#SBATCH --qos=long_investor
#SBATCH --partition=genomics
#SBATCH --time=0-08:00:00
#SBATCH --output=./output/%A_%a-slurm.out
#SBATCH --mail-type=ALL
# #SBATCH --mail-user=foo@bar.com

i=$SLURM_ARRAY_TASK_ID
echo $i

ResultsFolder="Results"

declare -a FILELIST=(PROTEIN_IDs_LIST)

id=${FILELIST[i]}

echo $id

# Set the path to your Julia executable
julia_exec="/cta/users/nkuru/julia-1.9.4/bin/julia"  # Adjust this path accordingly

# Run Julia script with the protein ID as a command-line argument
LD_LIBRARY_PATH="" $julia_exec julia_code.jl $id
