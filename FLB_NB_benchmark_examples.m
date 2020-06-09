% the examples below are for the C. glutamicum iJM658 model
%%
% enumerate synthetic lethals up to size 5 for mue threshold subselection 1,
% indicators, split v, combined z and irrev_geq
% runs only the FLB approach because irrev_geq is not supported by NB
FLB_NB_benchmarks('iJM658_synthetic_lethals_input', 1, 1, 5, 5, 0, 1, 1, true, true, true, 1e-6, 1e-6, 1e-7, 72000)

%%
% enumerate synthetic lethals up to size 5 for mue threshold subselection 2,
% indicators, split v, separate z and without irrev_geq
% runs the FLB and NB approach
FLB_NB_benchmarks('iJM658_synthetic_lethals_input', 2, 2, 5, 5, 0, 1, 1, true, false, false, 1e-6, 1e-6, 1e-7, 72000)

%%
% enumerate MCS up to size 8 for species index subselection 38 with seed 5,
% indicators, split v, combined z and irrev_geq
% runs only the FLB approach (stoichmat and kernel desired) because
% irrev_geq is not supported by NB
FLB_NB_benchmarks('iJM658_mcs_input', 38, 38, 5, 8, 0, 1, 1, true, true, true, 1e-6, 1e-6, 1e-7, 72000)
