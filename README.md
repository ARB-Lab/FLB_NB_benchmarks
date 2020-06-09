# FLB_NB_benchmarks
Functions and data for running the benchmarks. Requires MATLAB and CellNetAnalyzer (http://www2.mpi-magdeburg.mpg.de/projects/cna/cna.html, version 2020.2 or later).
The main function is FLB_NB_benchmarks which executes the enumerations for the problems provided by an input file. There are six input files:

iJM658_mcs_input.mat - growth-coupled production for 94 metabolites in the iJM658 model
iJM658_synthetic_lethals_input.mat - synthetic lethals for 11 mue thresholds in the iJM658 model
iJO1366_mcs_input.mat - growth-coupled production for 49 metabolites in the iJO1366 model
iJO1366_synthetic_lethals_input.mat - synthetic lethals for 11 mue thresholds in the iJO1366 model
iMM904_mcs_input.mat - growth-coupled production for 45 metabolites in the iMM904 model
iMM904_synthetic_lethals_input.mat - synthetic lethals for 11 mue thresholds in the iMM904 model

Some examples how to execute the enumerations are given in FLB_NB_benchmark_examples.m.
