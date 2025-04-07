# Additional Analyses

This repository contains supplementary analyses conducted in response to reviewer suggestions for our study. Below is a detailed description of each analysis:

## Analysis 1: Categorization of MSAs

**Objective:** To categorize proteins based on the number of sequences, number of positions, and total tree length as suggested by the reviewer.

**Contents:** All scripts and outputs for this analysis are located in the [MSACategories folder](https://github.com/CompGenomeLab/PHACE/tree/main/ExtraAnalyses/MSACategories).

## Analysis 2: CoMap Pairwise Version Execution

**Objective:** To assess performance differences between the pairwise and clustering versions of the CoMap algorithm.

**Procedure:**
We ran the pairwise version from Dutheil 2005, but results were obtained for only 149 out of 653 proteins in a two-step running process.

### Step-by-Step Results

#### Initial Run:
- **Script Used:** [script1.txt](https://github.com/CompGenomeLab/PHACE/blob/main/ExtraAnalyses/CoMapPairwise/script1.txt)
- **Outcomes:**
  - 156 proteins completed without issues.
  - 459 proteins encountered an error suggesting the removal of saturated sites: "ERROR!!! You may want to try `input.sequence.remove_saturated_sites = yes` to ignore positions with likelihood 0."
  - 8 proteins experienced memory errors: "Detected 1 oom-kill event(s) in StepId=117299.batch."
  - 29 proteins had unexpected errors.

#### Second Run (Adjusted based on CoMap's suggestion):
- **Script Used:** [script2.txt](https://github.com/CompGenomeLab/PHACE/blob/main/ExtraAnalyses/CoMapPairwise/script2.txt)
- **Adjustment:** Added the command `input.sequence.remove_saturated_sites = yes`.
- **Outcomes:**
  - 148 proteins completed, with 143 overlapping results from the initial run, matching the first outcomes perfectly.
  - Remaining issues included unexpected errors for 359 proteins, premature termination of 90 jobs after 11 days without conclusion, and persistent memory errors for 32 proteins despite increased memory allocation.

**Contents:** The codes and detailed procedure documents are available in the [CoMapPairwise folder](https://github.com/CompGenomeLab/PHACE/tree/main/ExtraAnalyses/CoMapPairwise).

For more details on the analyses and to view the results, please refer to the respective folders. Each folder contains all necessary scripts and documentation to understand and reproduce the analyses.
