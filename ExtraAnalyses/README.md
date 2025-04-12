# Additional Analyses

This repository contains supplementary analyses conducted in response to reviewer suggestions for our study. Below is a detailed description of each analysis:

## Analysis 1: AUPR Comparison of All Tools

**Objective:** Comparing and reporting the precision and recall values for all methods as suggested by the reviewer.

**Contents:** All scripts and data for this analysis are located in the [AUPR Plots folder](https://github.com/CompGenomeLab/PHACE/tree/main/ExtraAnalyses/AUPR_Plots).

## Analysis 2: Categorization of MSAs

**Objective:** To categorize proteins based on the number of sequences, number of positions, and total tree length as suggested by the reviewer.
We conducted analyses in three categories: number of sequences, number of taxa, and total tree length. We categorized each parameter into three groups—small, medium, and large. The thresholds for the 'small' group are set at 250 for number of sequences, 500 for number of positions, and 100 for total tree length. Conversely, the 'large' groups are defined by values greater than 750 for number of sequences, 1000 for number of positions, and 200 for total tree length. Values falling between these thresholds are classified as 'medium'.

**Contents:** All scripts and data for this analysis are located in the [MSACategories folder](https://github.com/CompGenomeLab/PHACE/tree/main/ExtraAnalyses/MSACategories).

## Analysis 3: CoMap Pairwise Version Execution

**Objective:** To assess performance differences between the pairwise and clustering versions of the CoMap algorithm.

**Procedure:**
We ran the pairwise version from Dutheil 2005, but results were obtained for only 161 out of 653 proteins in a two-step running process. The list of proteins with CoMap results are given as [CoMap_PairwiseResult.txt](https://github.com/CompGenomeLab/PHACE/blob/main/ExtraAnalyses/CoMapPairwise/CoMap_PairwiseResult.txt).

### Step-by-Step Results

#### Initial Run:
- **Script Used:** [comap.bpp](https://github.com/CompGenomeLab/PHACE/blob/main/ExtraAnalyses/CoMapPairwise/comap.bpp)
- **Outcomes:**
  - 156 proteins completed without issues.
  - 459 proteins encountered an error suggesting the removal of saturated sites: "ERROR!!! You may want to try `input.sequence.remove_saturated_sites = yes` to ignore positions with likelihood 0."
  - 8 proteins experienced memory errors: "Detected 1 oom-kill event(s) in StepId=117299.batch."
  - 29 proteins had unexpected errors.

#### Second Run (Adjusted based on CoMap's suggestion):
- **Script Used:** [comap2.bpp](https://github.com/CompGenomeLab/PHACE/blob/main/ExtraAnalyses/CoMapPairwise/comap2.bpp)
- **Adjustment:** Added the command `input.sequence.remove_saturated_sites = yes`.
- **Outcomes:**
  - 148 proteins completed, with 143 overlapping results from the initial run, matching the first outcomes perfectly.
  - Remaining issues included unexpected errors for 383 proteins, premature termination of 90 jobs after 11 days without conclusion, and persistent memory errors for 32 proteins despite increased memory allocation.

**Contents:** The codes and data required to reproduce the results are available in the [CoMapPairwise folder](https://github.com/CompGenomeLab/PHACE/tree/main/ExtraAnalyses/CoMapPairwise).

For more details on the analyses and to view the results, please refer to the respective folders and Supplementary Material of our paper.
