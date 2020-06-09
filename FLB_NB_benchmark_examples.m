% the FLB_NB_benchmarks function saves the results as cell arrays in a file whose name
% is composed with this pattern: [input_file]_s[seed]_[idx(first)]_[idx(last)].mat 
% where idx is a variable in the input_file and input_file, seed, first,
% last are parameters of the function; the cell arrays saved are:
% comp_time: computation time for FLB stoichmat desired
% comp_time2: computation time for FLB kernel desired
% comp_time3: computation time for NB stoichmat desired
% comp_time4: computation time for NB kernel desired
% mcs: MCS for FLB stoichmat desired
% mcs2: MCS for FLB kernel desired
% mcs3: MCS for NB stoichmat desired
% mcs4: MCS for NB kernel desired
% the indices of the cells with the results are idx(first:last)
% the examples below are for the C. glutamicum iJM658 model
%%
% enumerate synthetic lethals up to size 5 for mue threshold subselection 1,
% indicators, split v, combined z and irrev_geq
% runs only the FLB approach because irrev_geq is not supported by NB
input_file= 'iJM658_synthetic_lethals_input';
first= 1;
last= 1;
seed= 5;
FLB_NB_benchmarks(input_file, first, last, seed, 5, 0, 1, 1, true, true, true, 1e-6, 1e-6, 1e-7, 72000)
% results
load(input_file, 'idx');
load(sprintf('%s_s%d_%d_%d.mat', input_file, seed, idx(first), idx(last)))
comp_time(:, idx(first:last))
mcs(:, idx(first:last))

%%
% enumerate synthetic lethals up to size 5 for mue threshold subselection 2,
% indicators, split v, separate z and without irrev_geq
% runs the FLB and NB approach
input_file= 'iJM658_synthetic_lethals_input';
first= 2;
last= 2;
seed= 5;
FLB_NB_benchmarks(input_file, first, last, seed, 5, 0, 1, 1, true, false, false, 1e-6, 1e-6, 1e-7, 72000)
% results
load(input_file, 'idx');
load(sprintf('%s_s%d_%d_%d.mat', input_file, seed, idx(first), idx(last)))
comp_time(:, idx(first:last))
comp_time3(:, idx(first:last))
mcs(:, idx(first:last))
mcs3(:, idx(first:last))

%%
% enumerate MCS up to size 8 for species index subselection 38 with seed 5,
% indicators, split v, combined z and irrev_geq
% runs only the FLB approach (stoichmat and kernel desired) because
% irrev_geq is not supported by NB
input_file= 'iJM658_mcs_input';
first= 38;
last= 38;
seed= 5;
FLB_NB_benchmarks(input_file, first, last, seed, 8, 0, 1, 1, true, true, true, 1e-6, 1e-6, 1e-7, 72000)
% results
load(input_file, 'idx');
load(sprintf('%s_s%d_%d_%d.mat', input_file, seed, idx(first), idx(last)))
comp_time(:, idx(first:last))
comp_time2(:, idx(first:last))
mcs(:, idx(first:last))
mcs2(:, idx(first:last))
