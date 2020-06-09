function iJM658_mcs_input

cplex_inner= setup_cplex_inner_class_access();

emphasizeAccuracy= cplex_inner.ParameterSet.constructor.newInstance([]);
emphasizeAccuracy.setParam(cplex_inner.BooleanParam.NumericalEmphasis, true);
emphasizeAccuracy.setParam(cplex_inner.DoubleParam.EpOpt, 1e-9);
emphasizeAccuracy.setParam(cplex_inner.DoubleParam.EpRHS, 1e-9);
emphasizeAccuracy.setParam(cplex_inner.IntParam.AdvInd, 0);
emphasizeAccuracy.setParam(cplex_inner.IntParam.RootAlg, cplex_inner.Algorithm.Dual);

emphasizeFeasibility= cplex_inner.ParameterSet.constructor.newInstance([]);
emphasizeFeasibility.setParam(cplex_inner.BooleanParam.NumericalEmphasis, true);
emphasizeFeasibility.setParam(cplex_inner.DoubleParam.EpRHS, 1e-9);
emphasizeFeasibility.setParam(cplex_inner.IntParam.RootAlg, cplex_inner.Algorithm.Dual);

%%  base network
load cglutJM_B.mat cglutA reac_off del_exch exchange_ind noKO3 atpm glucose mue
irr= cglutA.reacMin >= 0;
st= cglutA.stoichMat;
reac_off= union(reac_off, del_exch);
exchange_ind= setdiff(exchange_ind, del_exch);
separate_reac= unique([noKO3, atpm, glucose, mue]);
separate_reac(end+1)= 0; % placeholder for index of the export reaction

%% MCS parameters
glucose_uptake_limit= 5.4;
min_yield= 0.1; % fraction of glucose that has to be converted to the product
min_mue= 0.01;
min_atpm= 3.2;

%% initialize/load variables
q= size(st, 1);
max_prod= NaN(1, q);
max_mue= NaN(1, q);
max_prod_rd= NaN(1, q);
max_mue_rd= NaN(1, q);
cuts= cell(1, q);
kn= cell(1, q);
inh= cell(1, q);
ub= cell(1, q);
des= cell(1, q);
db= cell(1, q);
gluc_up= NaN(1, q);
flux_lb= cell(1, q);
flux_ub= cell(1, q);
fvalb= cell(1, q);
fvaub= cell(1, q);
sub_rat= cell(1, q);
val_red= NaN(1, q);
load cglutJM_B_10_KO3_dxS.mat rd_rat irrev_rd_rat sub_rat val_red flux_lb flux_ub proper_mcs

%% select species indices for which the MCS will be calculated
sz= cellfun(@sum, proper_mcs);
idx= find(sz <= 8 & sz > 0);

%% network compression (unless compressed networks have been loaded) and setup of main MILP inputs
for i= idx
  i
  production= find(st(i, exchange_ind));
  if isempty(production)
    stex= [st, zeros(size(st, 1), 1)]; % export reaction is the last reaction
    irrex= [irr; true];
    stex(i, end)= -1;
    production= size(st, 2) + 1;
  else
    stex= st;
    irrex= irr;
    production= exchange_ind(production);
  end
  separate_reac(end)= production;
  
  q= size(stex, 2);
  cgp= Cplex();
  cgp.Param.emphasis.numerical.Cur= 1;
  cgp.Param.preprocessing.presolve.Cur= 0;
  cgp.Param.simplex.tolerances.optimality.Cur= 1e-9;
  cgp.Param.simplex.tolerances.feasibility.Cur= 1e-9;
  cgp.Model.A= stex;
  cgp.Model.ub= Inf(q, 1);
  cgp.Model.lb= -Inf(q, 1);
  cgp.Model.lb(irrex)= 0;
  cgp.Model.lhs= zeros(size(st, 1), 1);
  cgp.Model.rhs= zeros(size(st, 1), 1);
  cgp.Model.lb(glucose)= -glucose_uptake_limit;
  cgp.Model.lb(atpm)= min_atpm;
  cgp.Model.lb(reac_off)= 0;
  cgp.Model.ub(reac_off)= 0;
  cgp.DisplayFunc= [];
  cgp.Model.sense= 'maximize';
  cgp.Model.obj= zeros(q, 1);
  cgp.Model.obj(production)= 1;
  res= cgp.solve();
  if res.status == 1
    max_prod(i)= res.objval;
    gluc_up(i)= -res.x(glucose);
  elseif res.status == 2
    max_prod(i)= Inf;
    gluc_up(i)= -res.x(glucose);
  else
    error('Unexpexted solution status');
  end
  cgp.Model.obj= zeros(q, 1);
  cgp.Model.obj(mue)= 1;
  res= cgp.solve();
  if res.status == 1
    max_mue(i)= res.objval;
  else
    error('Unexpexted solution status');
  end
  
  if isinf(max_prod(i)) || (max_prod(i) < 1e-9)
    continue; % metabolite i is not limited by glucose or not producible from glucose
  end
  
  prod_min_yield= max_prod(i)/gluc_up(i)*min_yield;
  if isempty(sub_rat{i})
    desf= zeros(0, size(stex, 2));
    desf(1, [mue, glucose])= [-1, -min_mue]; % growth yield
    desf(2, glucose)= -1;
    desf(3, [production, glucose])= [-1, -prod_min_yield];
    desf(4, atpm)= -1;
    dbf= [0; glucose_uptake_limit; 0; -min_atpm]; % growth yield

    cnap= struct();
    cnap.stoichMat= stex;
    cnap= CNAgenerateMFNetwork(cnap);
    cnap.reacMin(~irrex)= -Inf;
    cnap.reacMin(reac_off)= 0;
    cnap.reacMax(reac_off)= 0;
    lpub= cnap.reacMax;
    lplb= cnap.reacMin;
    lhs= [zeros(cnap.nums, 1); -Inf(length(dbf), 1)];
    rhs= [zeros(cnap.nums, 1); dbf];
    [ind_i, ind_j, val]= find([stex; desf]);
    % do FVA with des so that the results can be used for validation later
    fvaj= CplexFVA.fva(ind_i-1, ind_j-1, val, lhs, rhs, lplb, lpub, emphasizeAccuracy);
    clear ind_i ind_j val;
    fvalb{i}= fvaj(1);
    fvaub{i}= fvaj(2);
    fva_tol= 1e-9;
    fvalb{i}(abs(fvalb{i}) < fva_tol)= 0;
    fvaub{i}(abs(fvaub{i}) < fva_tol)= 0;
    same= fvalb{i} == fvaub{i};
    blocked_fva= same & (fvalb{i} == 0);
 
    separate_reac(blocked_fva(separate_reac))= []; % do not keep blocked reactions separately
    [rd, irrev_rd, sub]= CNAcompressMFNetwork(cnap, separate_reac, [], 1, 0, 1, find(blocked_fva), 1);
    sub= sparse(sub');
    rd= sparse(rd);
    rd_rat{i}= rd;
    irrev_rd= logical(irrev_rd);
    irrev_rd_rat{i}= irrev_rd;
    sub_rat{i}= sub;
  else
    rd= rd_rat{i};
    irrev_rd= irrev_rd_rat{i};
    sub= sub_rat{i};
  end

  muer= find(sub(:, mue));
  glucr= find(sub(:, glucose));
  prodr= find(sub(:, production));
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
  if res.status == 1
    max_mue_rd(i)= res.objval;
  else
    error('LP problem');
  end
  cgp.Model.obj= zeros(q, 1);
  cgp.Model.obj(prodr)= 1;
  res= cgp.solve();
  if res.status == 1
    max_prod_rd(i)= res.objval;
%     gluc_up(i)= -res.x(glucr);
  else
    error('LP problem');
  end
  
  inh{i}= zeros(0, size(rd, 2));
  inh{i}(1, glucr)= -1;
  inh{i}(2, [prodr, glucr])= [1, prod_min_yield]; % +prod_min_yield because glucose influx has a negative sign!
  inh{i}(3, atpmr)= -1;
  ub{i}= [glucose_uptake_limit; 0; -min_atpm];
  des{i}= zeros(0, size(rd, 2));
  des{i}(1, [muer, glucr])= [-1, -min_mue]; % growth yield
  des{i}(2, glucr)= -1;
  des{i}(3, [prodr, glucr])= [-1, -prod_min_yield];
  des{i}(4, atpmr)= -1;
  db{i}= [0; glucose_uptake_limit; 0; -min_atpm]; % growth yield

  [m, n]= size(rd);
  lplb= -Inf(n, 1);
  lplb(irrev_rd ~= 0)= 0;
  lpub= Inf(n, 1);
  lhs= [zeros(m, 1); -Inf(length(db{i}), 1)];
  rhs= [zeros(m, 1); db{i}];
  [ind_i, ind_j, val]= find([rd; des{i}]);
  valj= CplexValidateMCS.validate(ind_i-1, ind_j-1, val, lhs, rhs, lplb, lpub, 0:n-1, 0:n-1, emphasizeFeasibility);
  fvajrd= CplexFVA.fva(ind_i-1, ind_j-1, val, lhs, rhs, lplb, lpub, emphasizeAccuracy);
  clear ind_i ind_j val;
  cuts{i}= valj(1)' == 0;
  ind= find(cuts{i} == true);
  ind= ind(flux_lb{i}(ind) > 0);
  if ~isempty(ind)
    if any(flux_lb{i}(ind) > 1e-9)
      disp('ERROR: Encountered non-essential reaction with lower bound > 1e-9; skipping.');
      continue;
    end
    disp('Warning: Encountered non-essential reaction with positive lower bound <= 1e-9; setting bound to 0.');
    flux_lb{i}(ind)= 0;
  end
  if isnan(val_red(i))
    val_red(i)= max(max(abs(validate_reduction_fva_bounds(sub, fvaj(1), fvaj(2), fvajrd(1),fvajrd(2)))));
  end
  
  flux_lb{i}= fvajrd(1);
  flux_ub{i}= fvajrd(2);
  flux_lb{i}(isinf(flux_lb{i}))= -2000;
  flux_ub{i}(isinf(flux_ub{i}))= 2000;
  fvalb{i}(isinf(fvalb{i}))= -2000;
  fvaub{i}(isinf(fvaub{i}))= 2000;
  cuts{i}(any(sub(:, noKO3), 2))= false;
  kn{i}= null_rat_efmtool(rd);
  kn{i}(abs(kn{i}) < 1e-10)= 0;
end

%%
save iJM658_mcs_input.mat rd_rat irrev_rd_rat sub_rat flux_lb flux_ub gluc_up max_prod cuts kn idx inh ub des db
