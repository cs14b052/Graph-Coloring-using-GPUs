Graph Coloring using GPUs

- git clone this repository and follow the instructions in Overview Document for step-by-step procedure on running various graph coloring programs.

The repository contains the following programs
  - sirg (our implementation)
  - chenGC (work of Chen et al., as mentioned in the paper)
  - csrcolor (using NVIDIA's cuSPARSE library, as mentioned in the paper)
  - Our Baseline algorithm (extension of Rokos et al. to GPUs with degree heuristic) - otherSIRGVariants/ folder
  - Our Baseline algorithm + Opt1 - otherSIRGVariants/ folder
  - sirg + adjColors array per vertex - otherSIRGVariants/ folder

The instructions for running each of them are present in the corresponding folders.


