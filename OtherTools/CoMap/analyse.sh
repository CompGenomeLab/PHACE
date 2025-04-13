#!/bin/sh
#SBATCH --job-name=CoMap
#SBATCH --account=investor
#SBATCH --nodes=1
#SBATCH --mem=250GB
#SBATCH --array=1-652
#SBATCH --qos=long_investor
#SBATCH --partition=genomics
#SBATCH --time=14-00:00:00
#SBATCH --output=./output/%A_%a-slurm.out
#SBATCH --mail-type=ALL
# #SBATCH --mail-user=foo@bar.com

k=1
i=$(($SLURM_ARRAY_TASK_ID + (($k-1))*1000))
echo $i

declare -a FILELIST=(PROTEIN_IDs_LIST)

echo ${FILELIST[i]}

id=${FILELIST[i]}

mkdir Results/New/${id}

DATA=$id
DATA_FILE="/cta/groups/adebali/phact_data_nurdan/MaskedMSA_RAxMLProteins/${id}_MaskedMSA.fasta"
TREE_FILE="/cta/groups/adebali/PHACT_Dataset/${id}/3_mltree/${id}.raxml.bestTree"
OUTPUT_INFOS_FILE=Results/New/${id}/${id}.infos
OUTPUT_TAGS_FILE=Results/New/${id}/${id}.tags.dnd
OUTPUT_TAGS_TRANSLATION_FILE=Results/New/${id}/${id}.tags.tln

echo ${TREE_FILE}
echo "Analysing dataset $DATA."

echo "Running CoMap with unweighted mapping:"
/cta/users/nkuru/Files_NotUsed_Freq/SENT/CoMap/comap param=comap.bpp \
      nijt=Uniformization \
      input.sequence.file=$DATA_FILE \
      input.tree.file=$TREE_FILE \
      output.infos=$OUTPUT_INFOS_FILE \
      output.tags.file=$OUTPUT_TAGS_FILE \
      output.tags.translation=$OUTPUT_TAGS_TRANSLATION_FILE \
      clustering.output.tree.file=Results/New/${id}/${DATA}_clust.dnd \
      clustering.output.groups.file=Results/${DATA}_groups.csv \
      clustering.null.output.file=Results/New/${id}/${DATA}_simulations.csv

echo "Done"
