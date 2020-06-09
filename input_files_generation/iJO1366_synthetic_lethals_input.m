function iJO1366_synthetic_lethals_input

cplex_inner= setup_cplex_inner_class_access();

emphasizeAccuracy= cplex_inner.ParameterSet.constructor.newInstance([]);
emphasizeAccuracy.setParam(cplex_inner.BooleanParam.NumericalEmphasis, true);
emphasizeAccuracy.setParam(cplex_inner.DoubleParam.EpOpt, 1e-9);
emphasizeAccuracy.setParam(cplex_inner.DoubleParam.EpRHS, 1e-9);
emphasizeAccuracy.setParam(cplex_inner.IntParam.RootAlg, cplex_inner.Algorithm.Dual);

%% base network
load orth1366D.mat orth1366Bs noKO atpm glucose mue o2in reac_off del_exchanges
st= orth1366Bs.stoichMat; % no external metabolites in orth1366Bs
irr= orth1366Bs.reacMin >= 0;
reac_off= union(reac_off, del_exchanges);
separate_reac= unique([noKO, atpm, glucose, mue, o2in]);

%% MCS parameters
glucose_uptake_limit= 15;
min_atpm= 3.15;

%% maximal growth rate
load orthD_10_dx.mat max_mue
max_mue= mode(max_mue);

%% network compression
coli_full= struct();
coli_full.stoichMat= st;
coli_full= CNAgenerateMFNetwork(coli_full);
coli_full.reacMin(~irr)= -Inf;
coli_full.reacMin(reac_off)= 0;
coli_full.reacMax(reac_off)= 0;
lpub= coli_full.reacMax;
lplb= coli_full.reacMin;
lhs= zeros(coli_full.nums, 1);
rhs= zeros(coli_full.nums, 1);
[ind_i, ind_j, val]= find(st);
fvaj= CplexFVA.fva(ind_i-1, ind_j-1, val, lhs, rhs, lplb, lpub, emphasizeAccuracy);
clear ind_i ind_j val;
fvalb= fvaj(1);
fvaub= fvaj(2);
fva_tol= 1e-9;
fvalb(abs(fvalb) < fva_tol)= 0;
fvaub(abs(fvaub) < fva_tol)= 0;
same= fvalb == fvaub;
blocked_fva= same & (fvalb == 0);
separate_reac(blocked_fva(separate_reac))= []; % do not keep blocked reactions separately
[rd, irrev_rd, sub]= CNAcompressMFNetwork(coli_full, separate_reac, [], 1, 0, 1, find(blocked_fva), 1);
sub= sparse(sub');
rd= sparse(rd);

%% compressed network
muer= find(sub(:, mue));
glucr= find(sub(:, glucose));
atpmr= find(sub(:, atpm));

cgp= Cplex();
q= size(rd, 2);
cgp.Param.emphasis.numerical.Cur= 1;
cgp.Param.advance.Cur= 0;
cgp.Param.preprocessing.presolve.Cur= 0;
cgp.Param.simplex.tolerances.optimality.Cur= 1e-9;
cgp.Param.simplex.tolerances.feasibility.Cur= 1e-9;
cgp.Model.A= rd;
cgp.Model.ub= Inf(q, 1);
cgp.Model.lb= -Inf(q, 1);
cgp.Model.lb(irrev_rd ~= 0)= 0;
cgp.Model.lb(glucr)= -glucose_uptake_limit;
cgp.Model.lb(atpmr)= min_atpm;
cgp.Model.lhs= zeros(size(rd, 1), 1);
cgp.Model.rhs= zeros(size(rd, 1), 1);
cgp.Model.sense= 'maximize';
cgp.DisplayFunc= [];
cgp.Model.obj= zeros(q, 1);
cgp.Model.obj(muer)= 1;
res= cgp.solve();
res.objval-max_mue
kn= null_rat_efmtool(rd);
cuts= ~any(sub(:, noKO), 2);

%% setup of mue thresholds and main MILP inputs 
pc_max_mue_levels= [10 20 30 40 50 60 70 80 90 1 0.1];
q= length(pc_max_mue_levels);
inh= cell(1, q);
ub= cell(1, q);
rd_rat= cell(1, q);
[rd_rat{:}]= deal(rd);
irrev_rd_rat= cell(1, q);
[irrev_rd_rat{:}]= deal(irrev_rd);
des= cell(1, q); % remain empty
db= cell(1, q); % remain empty
flux_lb= cell(1, q); % remain empty
flux_ub= cell(1, q); % remain empty
cuts= repmat({cuts}, 1, q);
kn= repmat({kn}, 1, q);
idx= 1:q;
for i= idx
  pc_max_mue= pc_max_mue_levels(i);
  inh{i}= zeros(0, size(rd, 2));
  inh{i}(1, glucr)= -1;
  inh{i}(2, muer)= -1;
  inh{i}(3, atpmr)= -1;
  ub{i}= [glucose_uptake_limit; -pc_max_mue/100*max_mue; -min_atpm];
end

%%
save iJO1366_synthetic_lethals_input.mat rd_rat irrev_rd_rat flux_lb flux_ub kn cuts inh ub des db idx
