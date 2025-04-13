#!/bin/bash
#SBATCH --job-name=CAPS
#SBATCH --array=1-652

k=1
i=$(($SLURM_ARRAY_TASK_ID + (($k-1))*1000))
echo $i

declare -a FILELIST=(Protein_IDs_List)

id=${FILELIST[i]}
echo $id

mkdir MSA/${id}
mkdir Tree/${id}

scp "/cta/groups/adebali/phact_data_nurdan/MaskedMSA_RAxMLProteins/${id}_MaskedMSA.fasta"  MSA/${id}
scp "/cta/groups/adebali/phact_data_nurdan/MaskedASR_RAxMLProteins/${id}/${id}.treefile" Tree/${id}

./caps -F MSA/${id} -T Tree/${id}
