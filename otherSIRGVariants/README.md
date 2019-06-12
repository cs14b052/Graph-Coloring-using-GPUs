SIRG (Scalable and Improved RIG algorithm for GPUs) Variants

The folder contains programs for the following:
  - Our Baseline algorithm (extension of Rokos et al. to GPUs with degree heuristic) 
  - Our Baseline algorithm + Opt1 
  - sirg + adjColors array per vertex 

Compiling: run 'make' (in src folder) to compile all the three programs.
To run individual programs, run the following:
  - 'make baseline' -> creates 'baseline' executable in bin/ folder
  - 'make baselineOpt1' -> creates 'baselineOpt1' executable in bin/ folder
  - 'make sirgpervertex' -> creates 'sirgpervertex' executable in bin/ folder

Running: All the programs take 2 arguments
  - filename of the input graph
  - an integer - 0/1
    - 0 for degree based heuristic (in conflict resolution)
    - 1 for vertex based heuristic (in conflict resolution)
