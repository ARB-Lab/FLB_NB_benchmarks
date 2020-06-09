function FLB_NB_benchmarks(input_file, first, last, seed, max_sz, bigM, threshold, split_level,...
  reduce_constraints, combined_z, irrev_geq, opt_tol, feas_tol, int_tol, time_limit)
% this function runs the MCS enumerations with the main MILP inputs
% together with the selected indices being read from input_file
% first:last specifies the subselection of indices for which the benchmarks are run
% seed: CPLEX random seed
% max_sz: enumeration is performed with successive poulate calls up to MCS size max_sz
% bigM, threshold, split_level, reduce_constraints, combined_z, irrev_geq:
% paramaters for ConstrainedMinimalCutSetsEnumeratorMILP
% opt_tol, feas_tol, int_tol: CPLEX tolerance settings
% time_limit: time limit in seconds for each MCS size level
% the results are saved in a file with the name
% [input_file]_s[seed]_[idx(first)]_[idx(last)].mat where idx is a variable
% in the input_file

solution_limit= Inf;
working_memory= 80000; % GB
load([input_file, '.mat'], 'rd_rat', 'irrev_rd_rat', 'flux_lb', 'flux_ub', 'cuts', 'kn',...
  'idx', 'inh', 'ub', 'des', 'db');
q= length(rd_rat);
comp_time= zeros(max_sz, q); % computation time for FLB stoichmat desired
comp_time2= zeros(max_sz, q); % computation time for FLB kernel desired
comp_time3= zeros(max_sz, q); % computation time for NB stoichmat desired
comp_time4= zeros(max_sz, q); % computation time for NB kernel desired
mcs= cell(1, q); % MCS for FLB stoichmat desired
mcs2= cell(1, q); % MCS for FLB kernel desired
mcs3= cell(1, q); % MCS for NB stoichmat desired
mcs4= cell(1, q); % MCS for NB kernel desired
matching_mcs= false(1, q);

%%
for i= idx(first:last)
  i
  rd= rd_rat{i};
  irrev_rd= irrev_rd_rat{i};
  
  % FLB stoichmat desired
  obj= ConstrainedMinimalCutSetsEnumeratorMILP(rd, irrev_rd, [], inh{i}, ub{i}, cuts{i},...
    flux_lb{i}, flux_ub{i}, des{i}, db{i}, [], bigM, threshold, split_level,...
    reduce_constraints, combined_z, irrev_geq);
  obj.cpx.Param.emphasis.numerical.Cur= 1;
  obj.cpx.Param.workmem.Cur= working_memory;
  obj.cpx.Param.simplex.tolerances.optimality.Cur= opt_tol;
  obj.cpx.Param.simplex.tolerances.feasibility.Cur= feas_tol;
  obj.cpx.Param.mip.tolerances.integrality.Cur= int_tol;
  obj.cpx.Param.randomseed.Cur= seed;
  obj.cpx.Param.output.clonelog.Cur= -1;
  obj.ev_size_lb= 1;
  for sz = 1:max_sz
    [obj, mcs_sz, ~, comp_time(sz, i)]= shortest_minimal_cut_sets(obj, sz, solution_limit, time_limit, 2);
    mcs{i}= [mcs{i}, mcs_sz];
  end
  mcs{i}= sparse(sortrows(logical(mcs{i}'))');

  if ~isempty(des{i})
    % FLB kernel desired
    obj= ConstrainedMinimalCutSetsEnumeratorMILP(rd, irrev_rd, {kn{i}, 0, 1}, inh{i}, ub{i}, cuts{i},...
      flux_lb{i}, flux_ub{i}, des{i}, db{i}, [], bigM, threshold, split_level,...
      reduce_constraints, combined_z, irrev_geq);
    obj.cpx.Param.emphasis.numerical.Cur= 1;
    obj.cpx.Param.workmem.Cur= working_memory;
    obj.cpx.Param.simplex.tolerances.optimality.Cur= opt_tol;
    obj.cpx.Param.simplex.tolerances.feasibility.Cur= feas_tol;
    obj.cpx.Param.mip.tolerances.integrality.Cur= int_tol;
    obj.cpx.Param.randomseed.Cur= seed;
    obj.cpx.Param.output.clonelog.Cur= -1;
    obj.ev_size_lb= 1;
    for sz = 1:max_sz
      [obj, mcs_sz, ~, comp_time2(sz, i)]= shortest_minimal_cut_sets(obj, sz, solution_limit, time_limit, 2);
      mcs2{i}= [mcs2{i}, mcs_sz];
    end
    mcs2{i}= sparse(sortrows(logical(mcs2{i}'))');
  end

  if ~irrev_geq % not supported by NB approach
    % NB stoichmat desired
    obj= ConstrainedMinimalCutSetsEnumeratorMILP(rd, irrev_rd, {kn{i}, 1, 0}, inh{i}, ub{i}, cuts{i},...
      flux_lb{i}, flux_ub{i}, des{i}, db{i}, [], bigM, threshold, split_level,...
      reduce_constraints, combined_z, irrev_geq);
    obj.cpx.Param.emphasis.numerical.Cur= 1;
    obj.cpx.Param.workmem.Cur= working_memory;
    obj.cpx.Param.simplex.tolerances.optimality.Cur= opt_tol;
    obj.cpx.Param.simplex.tolerances.feasibility.Cur= feas_tol;
    obj.cpx.Param.mip.tolerances.integrality.Cur= int_tol;
    obj.cpx.Param.randomseed.Cur= seed;
    obj.cpx.Param.output.clonelog.Cur= -1;
    obj.ev_size_lb= 1;
    for sz = 1:max_sz
      [obj, mcs_sz, ~, comp_time3(sz, i)]= shortest_minimal_cut_sets(obj, sz, solution_limit, time_limit, 2);
      mcs3{i}= [mcs3{i}, mcs_sz];
    end
    mcs3{i}= sparse(sortrows(logical(mcs3{i}'))');
    
    if ~isempty(des{i})
      % NB kernel desired
      obj= ConstrainedMinimalCutSetsEnumeratorMILP(rd, irrev_rd, {kn{i}, 1, 1}, inh{i}, ub{i}, cuts{i},...
        flux_lb{i}, flux_ub{i}, des{i}, db{i}, [], bigM, threshold, split_level,...
        reduce_constraints, combined_z, irrev_geq);
      obj.cpx.Param.emphasis.numerical.Cur= 1;
      obj.cpx.Param.workmem.Cur= working_memory;
      obj.cpx.Param.simplex.tolerances.optimality.Cur= opt_tol;
      obj.cpx.Param.simplex.tolerances.feasibility.Cur= feas_tol;
      obj.cpx.Param.mip.tolerances.integrality.Cur= int_tol;
      obj.cpx.Param.randomseed.Cur= seed;
      obj.cpx.Param.output.clonelog.Cur= -1;
      obj.ev_size_lb= 1;
      for sz = 1:max_sz
        [obj, mcs_sz, ~, comp_time4(sz, i)]= shortest_minimal_cut_sets(obj, sz, solution_limit, time_limit, 2);
        mcs4{i}= [mcs4{i}, mcs_sz];
      end
      mcs4{i}= sparse(sortrows(logical(mcs4{i}'))');
    end
  end

  if isempty(des{i})
    if ~irrev_geq
      matching_mcs(i)= isequal(mcs{i}, mcs3{i});
    end
  else
    if irrev_geq
      matching_mcs(i)= isequal(mcs{i}, mcs2{i});
    else
      matching_mcs(i)= isequal(mcs{i}, mcs2{i}) && isequal(mcs{i}, mcs3{i}) && isequal(mcs{i}, mcs4{i});
    end
  end
  
  save(sprintf('%s_s%d_%d_%d.mat', input_file, seed, idx(first), idx(last)), 'idx', 'matching_mcs',...
    'mcs', 'comp_time', 'mcs2', 'comp_time2', 'mcs3', 'comp_time3', 'mcs4', 'comp_time4');
end
disp('ALL FINISHED');
