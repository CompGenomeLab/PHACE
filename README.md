
#  PHACE: Phylogeny-Aware Detection of Molecular Coevolution


The coevolution trends of amino acids within or between genes offer valuable insights into protein structure and function. Existing tools for uncovering
coevolutionary signals primarily rely on multiple sequence alignments (MSAs), often neglecting considerations of phylogenetic relatedness and shared 
evolutionary history. 

Here, we present a novel approach based on the substitution mapping of amino acid changes onto the phylogenetic tree. We categorize 
amino acids into two groups: 'tolerable' and 'intolerable,' assigned to each position based on the position dynamics concerning the observed amino acids. 
Amino acids deemed 'tolerable' are those observed phylogenetically independently and multiple times at a specific position, signifying the position's 
tolerance to that alteration. Gaps are regarded as a third character type, and we only take phylogenetically independent altered gap characters into 
consideration. 

Our algorithm is based on a tree traversal process through the nodes and computes the total amount of substitution per branch based on 
the probability differences of two groups of amino acids and gaps between neighboring nodes. To mitigate false coevolution signals from unaligned regions, 
we employ an MSA-masking approach. 

When compared to tools utilizing phylogeny (e.g., CAPS and CoMap) and state-of-the-art MSA-based approaches (DCA, GaussDCA, 
PSICOV, and MIp), our method exhibits significantly superior accuracy in identifying coevolving position pairs, as measured by statistical metrics including 
MCC, AUC, and F1 score. The success of PHACE stems from our capacity to account for the often-overlooked phylogenetic dependency.

![Outline of the PHACE algorithm](https://github.com/nurdannkuru/PHACE/raw/main/Outline.png)
                                                    **Figure 1. Outline of the PHACE algorithm**
PHACE utilizes the original MSA and ML phylogenetic tree to cluster amino acids into "tolerable" and "intolerable" groups, resulting in MSA1. To address issues with gapped leaves and obtain accurate coevolution signals, MSA2 is created to distinguish amino acids from gaps. This information is used to update substitution rates per branch from MSA1. The final MSA is used to construct a matrix detailing changes per branch per position and branch diversity. PHACE score is calculated using a weighted concordance correlation coefficient. (Pos. 126-130, distance: 6.54)

## System Requirements

- **R** (version 3.6 or higher)
- **IQ-TREE2** for Ancestral Sequence Reconstruction (ASR)
- **Linux/Unix environment** (recommended for optimal performance)

## Installation

### 1. Install Required R Packages

```r
# Core PHACE dependencies
install.packages(c("ape", "Biostrings", "tidytree", "stringr", "dplyr", "bio3d", "mltools", "irr"))

# Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("Biostrings"))

# Install additional packages for analysis
install.packages(c("ggplot2", "ggsignif", "cowplot", "AUC", "PRROC", "caret", "spaa", "tidyr", "broom", "RColorBrewer", "gridExtra", "Peptides"))
```

### 2. Install IQ-TREE2

Download and install IQ-TREE2 from: https://www.iqtree.org/

### 3. Clone the Repository

```bash
git clone https://github.com/CompGenomeLab/PHACE.git
cd PHACE
```

## Repository Structure

```
PHACE/
├── PHACE_Codes/           # Main PHACE implementation scripts
│   ├── ToleranceScore.R   # Calculate tolerance scores
│   ├── MSA1.R            # Generate MSA1
│   ├── MSA2.R            # Generate MSA2
│   ├── Part1_MSA1.R      # Process MSA1 results
│   ├── Part1_MSA2.R      # Process MSA2 results
│   ├── GetTotalChangeMatrix.R # Merge matrices
│   └── PHACE.R           # Main PHACE algorithm
├── Data/                  # Configuration files
│   ├── vals_MSA1.txt     # MSA1 substitution model
│   └── vals_MSA2.txt     # MSA2 substitution model
├── AnalysisCodes/         # Performance evaluation scripts
├── ExtraAnalyses/         # Additional analyses and comparisons
├── ManuscriptFigures/     # Figure generation scripts
├── OtherTools/           # Comparison tools (CAPS, CoMap, etc.)
├── PDB/                  # PDB structure analysis
└── README.md
```

## Input Requirements

To run PHACE, you need:
1. **Multiple Sequence Alignment (MSA)** in FASTA format
2. **Phylogenetic tree** in Newick format
3. **Ancestral Sequence Reconstruction (ASR)** outputs from IQ-TREE2

For assistance with generating these inputs, please refer to the [PHACT Repository](https://github.com/CompGenomeLab/PHACT).

## How to Obtain PHACE Results

Here, we assume you have MSA, a phylogenetic tree, and ASR outputs for the protein of interest. These data should be produced similarly to our previous approach. If you need assistance with this part, please refer to the [PHACT Repository](https://github.com/CompGenomeLab/PHACT).

### Step-by-Step Workflow

#### Step 1: MSA1 Processing

1. **Calculate tolerance scores** per amino acid per position using [ToleranceScore.R](https://github.com/nurdannkuru/PHACE/blob/main/PHACE_Codes/ToleranceScore.R).

2. **Generate MSA1** using [MSA1.R](https://github.com/nurdannkuru/PHACE/blob/main/PHACE_Codes/MSA1.R), which comprises three characters:
   - C (dominant amino acids)
   - A (alternate amino acids) 
   - - (gap)

3. **Perform Ancestral Sequence Reconstruction (ASR)** with IQ-TREE2:
   ```bash
   iqtree2 -s ${file_fasta} -te ${file_nwk} -blfix -m Data/vals_MSA1.txt -asr --prefix ${id}_MSA1 --safe
   ```

4. **Construct the initial matrix** using [Part1_MSA1.R](https://github.com/nurdannkuru/PHACE/blob/main/PHACE_Codes/Part1_MSA1.R) to account for total changes per branch over MSA1.

#### Step 2: MSA2 Processing

1. **Generate MSA2** using [MSA2.R](https://github.com/nurdannkuru/PHACE/blob/main/PHACE_Codes/MSA2.R), which includes two characters:
   - C (all amino acids)
   - G (gap)

2. **Execute ASR** for MSA2:
   ```bash
   iqtree2 -s ${file_fasta} -te ${file_nwk} -blfix -m Data/vals_MSA2.txt -asr --prefix ${id}_MSA2 --safe
   ```

3. **Develop the secondary matrix** using [Part1_MSA2.R](https://github.com/nurdannkuru/PHACE/blob/main/PHACE_Codes/Part1_MSA2.R) to identify independent gap alterations.

#### Step 3: Final PHACE Calculation

1. **Merge the matrices** obtained from MSA1 and MSA2 using [GetTotalChangeMatrix.R](https://github.com/nurdannkuru/PHACE/blob/main/PHACE_Codes/GetTotalChangeMatrix.R).

2. **Execute the final PHACE algorithm** using [PHACE.R](https://github.com/nurdannkuru/PHACE/blob/main/PHACE_Codes/PHACE.R) to obtain PHACE results.

### Example Usage

```r
# Set your protein ID
id <- "your_protein_id"

# Run the complete PHACE pipeline
source("PHACE_Codes/ToleranceScore.R")
source("PHACE_Codes/MSA1.R")
source("PHACE_Codes/Part1_MSA1.R")
source("PHACE_Codes/MSA2.R") 
source("PHACE_Codes/Part1_MSA2.R")
source("PHACE_Codes/GetTotalChangeMatrix.R")
source("PHACE_Codes/PHACE.R")
```

## Additional Analyses

### Performance Evaluation
- **ROC Comparisons**: [AnalysisCodes/ROC_Comparisons.R](https://github.com/nurdannkuru/PHACE/blob/main/AnalysisCodes/ROC_Comparisons.R)
- **MCC/F1 Score Comparisons**: [AnalysisCodes/MCC_F1Score_Comparisons.R](https://github.com/nurdannkuru/PHACE/blob/main/AnalysisCodes/MCC_F1Score_Comparisons.R)

### Manuscript Figures
- **Figure 3**: [ManuscriptFigures/Figure3.R](https://github.com/nurdannkuru/PHACE/blob/main/ManuscriptFigures/Figure3.R)
- **Figure 4**: [ManuscriptFigures/Figure4.R](https://github.com/nurdannkuru/PHACE/blob/main/ManuscriptFigures/Figure4.R)
- **Figure 5**: [ManuscriptFigures/Figure5.R](https://github.com/nurdannkuru/PHACE/blob/main/ManuscriptFigures/Figure5.R)

### Comparison Tools
The `OtherTools/` directory contains implementations of comparison methods:
- **CAPS**: Coevolution analysis using protein sequences
- **CoMap**: Detecting groups of coevolving positions
- **DCA**: Direct-coupling analysis
- **GaussDCA**: Multivariate Gaussian modeling
- **MIp**: Mutual information without phylogeny influence
- **PSICOV**: Precise structural contact prediction

See [OtherTools/README.md](https://github.com/nurdannkuru/PHACE/blob/main/OtherTools/README.md) for detailed information.

### Extra Analyses
Additional analyses conducted in response to reviewer suggestions are available in the `ExtraAnalyses/` directory:
- **AUPR Comparisons**: Precision-recall analysis
- **MSA Categorization**: Performance across different MSA characteristics
- **CoMap Pairwise Analysis**: Comparison of CoMap versions

See [ExtraAnalyses/README.md](https://github.com/nurdannkuru/PHACE/blob/main/ExtraAnalyses/README.md) for details.

## Results

Result for 652 proteins is provided in Figure 2.

![Result](https://github.com/nurdannkuru/PHACE/raw/main/Result.png)
                              **Figure 2. Comparison of all tools over a common set in terms of AUC**

## Troubleshooting

### Common Issues

1. **Missing R packages**: Install all required packages listed in the Installation section
2. **IQ-TREE2 not found**: Ensure IQ-TREE2 is installed and accessible in your PATH
3. **Memory issues**: For large datasets, consider increasing R memory limits
4. **File format errors**: Ensure input files are in correct FASTA and Newick formats

### Performance Tips

- Use a Linux/Unix environment for optimal performance
- Ensure sufficient disk space for temporary files
- For large datasets, consider running on a high-performance computing cluster

## Data Availability

All data generated in this study and all benchmark analysis scripts and source codes for PHACE are available at https://github.com/CompGenomeLab/PHACE. The PHACE predictions for the 652 proteins used in this manuscript are provided at [https://zenodo.org/records/14038143](https://zenodo.org/records/14043199).

## Contributing

We welcome contributions! Please feel free to submit issues, feature requests, or pull requests.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Citing this work

Kuru N., Adebali O. (2025). PHACE: Phylogeny-Aware Detection of Molecular Coevolution. Molecular Biology and Evolution, 42(7), msaf150. 
https://doi.org/10.1093/molbev/msaf150

