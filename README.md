
# PHACE - Phylogeny-Aware Co-Evolution Algorithm


The co-evolution trends of amino acids within or between genes offer valuable insights into protein structure and function. Existing tools for uncovering
co-evolutionary signals primarily rely on multiple sequence alignments (MSAs), often neglecting considerations of phylogenetic relatedness and shared 
evolutionary history. 

Here, we present a novel approach based on the substitution mapping of amino acid changes onto the phylogenetic tree. We categorize 
amino acids into two groups: 'tolerable' and 'intolerable,' assigned to each position based on the position dynamics concerning the observed amino acids. 
Amino acids deemed 'tolerable' are those observed phylogenetically independently and multiple times at a specific position, signifying the position's 
tolerance to that alteration. Gaps are regarded as a third character type, and we only take phylogenetically independent altered gap characters into 
consideration. 

Our algorithm is based on a tree traversal process through the nodes and computes the total amount of substitution per branch based on 
the probability differences of two groups of amino acids and gaps between neighboring nodes. To mitigate false co-evolution signals from unaligned regions, 
we employ an MSA-masking approach. 

When compared to tools utilizing phylogeny (e.g., CAPS and CoMap) and state-of-the-art MSA-based approaches (DCA, GaussDCA, 
PSICOV, and MIp), our method exhibits significantly superior accuracy in identifying co-evolving position pairs, as measured by statistical metrics including 
MCC, AUC, and F1 score. The success of PHACE stems from our capacity to account for the often-overlooked phylogenetic dependency.

![Outline of the PHACE algorithm](https://github.com/nurdannkuru/PHACE/raw/main/Outline.png)
                                                    **Figure 1. Outline of the PHACE algorithm**
PHACE utilizes the original MSA and ML phylogenetic tree to cluster amino acids into "tolerable" and "intolerable" groups, resulting in MSA1. To address issues with gapped leaves and obtain accurate co-evolution signals, MSA2 is created to distinguish amino acids from gaps. This information is used to update substitution rates per branch from MSA1. The final MSA is used to construct a matrix detailing changes per branch per position and branch diversity. PHACE score is calculated using a weighted concordance correlation coefficient. (Pos. 126-130, distance: 6.54)



# How to Obtain PHACE Results

Here, we assume you have MSA, a phylogenetic tree, and ASR outputs for the protein of interest. These data should be produced similarly to our previous approach. If you need assistance with this part, please refer to the [PHACT Repository](https://github.com/CompGenomeLab/PHACT).

#### MSA1

1. Calculate the tolerance score per amino acid per position using [ToleranceScore.R](https://github.com/nurdannkuru/PHACE/blob/main/PHACE_Codes/ToleranceScore.R).
2. Generate MSA1, which comprises three characters: C (dominant amino acids), A (alternate amino acids), and - (gap), based on the tolerance scores computed in the previous step. Refer to the [MSA1_Code](https://github.com/nurdannkuru/PHACE/blob/main/PHACE_Codes/Part1_MSA1.R) for implementation.
3. Perform Ancestral Sequence Reconstruction (ASR) with the following command:

                                 iqtree2 -s ${file_fasta} -te ${file_nwk} -blfix -m Data/vals_MSA1.txt -asr --prefix ${id}_MSA1 --safe
   

4. Construct the initial matrix, accounting for total changes per branch over MSA1. Refer to [Part1_MSA1](https://github.com/nurdannkuru/PHACE/blob/main/PHACE_Codes/Part1_MSA1.R) for details.

#### MSA2

1. Formulate MSA2, which includes two characters: C (all amino acids) and G (gap).
2. Execute ASR with the following command:

                                 iqtree2 -s ${file_fasta} -te ${file_nwk} -blfix -m Data/vals_MSA2.txt-asr --prefix ${id}_MSA2 --safe


3. Develop the secondary matrix to identify independent gap alterations. Refer to [Part1_MSA2](https://github.com/nurdannkuru/PHACE/blob/main/PHACE_Codes/Part1_MSA2.R) for the procedure.

#### Final Step

1. Merge the matrices obtained from MSA1 and MSA2. Refer to [GetTotalChangeMatrix.R](https://github.com/nurdannkuru/PHACE/blob/main/PHACE_Codes/GetTotalChangeMatrix.R) for implementation details.
2. Execute the final code to acquire PHACE results. Refer to [PHACE.R](https://github.com/nurdannkuru/PHACE/blob/main/PHACE_Codes/PHACE.R) for the code.

# Results

Result for 652 proteins is provided in Figure 2.

![Result](https://github.com/nurdannkuru/PHACE/raw/main/Result.png)
                              **Figure 2. Comparison of all tools over a common set in terms of AUC**

# Data Availability

All data generated in this study and all benchmark analysis scripts and source codes for PHACE are available at https://github.com/CompGenomeLab/PHACE. The PHACE predictions for the 652 proteins used in this manuscript are provided at [https://zenodo.org/records/14038143](https://zenodo.org/records/14043199).

# Citing this work

Kuru N., Adebali O. (2024). PHACE: Phylogeny-Aware Detection of Molecular Co-Evolution. https://doi.org/10.1093/molbev/msaf150

