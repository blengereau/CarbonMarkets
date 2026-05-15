% 9. Compute direct and indirect effects of carbon tax shock
shock = 'eta_tau'; 

%% Pull steady states
endo = cellstr(M_.endo_names);
getSS = @(x) oo_.steady_state(strcmp(endo,x));
c1_ss = getSS('c_1');    c2_ss = getSS('c_2');
p1_ss = getSS('p_1');    p2_ss = getSS('p_2');
pE_ss = getSS('p_e');    e1_ss = getSS('e_1');   e2_ss = getSS('e_2');

%%% Pull % IRFs
c1h = oo_.irfs.(['c1_pct_'  shock])/100;
c2h = oo_.irfs.(['c2_pct_'  shock])/100;
p1h = oo_.irfs.(['p1_pct_'  shock])/100;
p2h = oo_.irfs.(['p2_pct_'  shock])/100;
r1h = oo_.irfs.(['res_1_pct_' shock])/100;   % assumes you computed res1_pct in the .mod
r2h = oo_.irfs.(['res_2_pct_' shock])/100;   % assumes you computed res2_pct in the .mod
peh = oo_.irfs.(['pe_pct_'  shock])/100;
e1h = oo_.irfs.(['e1_pct_'  shock])/100;
e2h = oo_.irfs.(['e2_pct_'  shock])/100;

chk1 = max(abs(c1h - (r1h - p1h)));   % should be ≪ 1e-6
chk2 = max(abs(c2h - (r2h - p2h)));
format long e
disp([chk1 chk2])
format short

%%% PDV weights
% beta = M_.params(strcmp(cellstr(M_.param_names),'beta'));   % or use 1/(1+rF_ss)
T    = length(c1h);
disc = beta.^(0:T-1)';

PDV = @(x) sum(disc.*x);

%%% Decomposition: Total = Direct (price) + Indirect (resources)
Total1    =  c1_ss * PDV(c1h);
Direct1   =  c1_ss * PDV(-p1h);
Indirect1 =  c1_ss * PDV( r1h);

Total2    =  c2_ss * PDV(c2h);
Direct2   =  c2_ss * PDV(-p2h);
Indirect2 =  c2_ss * PDV( r2h);

gap1 = abs(Total1 - (Direct1 + Indirect1));
gap2 = abs(Total2 - (Direct2 + Indirect2));
fprintf('gaps:', gap1, gap2);

ShareDir1 = Direct1/(Total1);   ShareInd1 = Indirect1/(Total1);
ShareDir2 = Direct2/(Total2);   ShareInd2 = Indirect2/(Total2);
fprintf('Direct share of total decline: Spenders %.1f%% | Savers %.1f%%\n',100*ShareDir1,100*ShareDir2);

% 10. Compute energy and non-energy share of nominal consumption drop

EE1 = pE_ss*e1_ss * PDV(peh + e1h);             % PDV Δ(p_e*e_1)
EE2 = pE_ss*e2_ss * PDV(peh + e2h);             % PDV Δ(p_e*e_2)
Cnom1 = p1_ss*c1_ss * PDV(p1h + c1h);           % PDV Δ(p_1*c_1)
Cnom2 = p2_ss*c2_ss * PDV(p2h + c2h);           % PDV Δ(p_2*c_2)
ShareEE1 = EE1/(-Cnom1);  ShareEE2 = EE2/(-Cnom2);
fprintf('Energy-spending share of nominal drop: Spenders %.1f%% | Savers %.1f%%\n',100*ShareEE1,100*ShareEE2);

% Pull x-IRFs
x1h = oo_.irfs.(['x1_pct_' shock])/100;
x2h = oo_.irfs.(['x2_pct_' shock])/100;

% Steady states
x1_ss = oo_.steady_state(strcmp(cellstr(M_.endo_names),'x_1'));
x2_ss = oo_.steady_state(strcmp(cellstr(M_.endo_names),'x_2'));

% PDV change in non-energy nominal spending (price of x is 1 in your model)
XX1 = x1_ss * PDV(x1h);   % PDV Δx_1
XX2 = x2_ss * PDV(x2h);   % PDV Δx_2

% Shares of the nominal DROP in p_i c_i
ShareXX1 = XX1/(-Cnom1);
ShareXX2 = XX2/(-Cnom2);

fprintf('Non-energy (x) share of nominal drop: Spenders %.1f%% | Savers %.1f%%\n', ...
        100*ShareXX1, 100*ShareXX2);

% 11. Compute uncompensated elasticity

% --- Pull % IRFs and convert to FRACTIONS
c1h = oo_.irfs.(['c1_pct_'  shock]) / 100;   % HtM consumption (fraction)
c2h = oo_.irfs.(['c2_pct_'  shock]) / 100;   % Savers consumption
e1h = oo_.irfs.(['e1_pct_'  shock])  / 100;  % HtM energy quantity
e2h = oo_.irfs.(['e2_pct_'  shock])  / 100;  % Savers energy quantity
p1h = oo_.irfs.(['p1_pct_'  shock])  / 100;  % HtM bundle price index
p2h = oo_.irfs.(['p2_pct_'  shock])  / 100;  % Savers bundle price index
peh = oo_.irfs.(['pe_pct_'  shock])  / 100;  % Energy price

gap = abs(c1h(1) - (r1h(1) - p1h(1)))
fprintf('gaps:', gap);

t0 = 1;  
etaM1_imp = e1h(t0) / peh(t0);   % HtM
etaM2_imp = e2h(t0) / peh(t0);   % Savers

den1 = peh(t0) - p1h(t0);      % Δ log(p_e/p_1) approx
den2 = peh(t0) - p2h(t0);

epsH1_imp = (c1h(t0) - e1h(t0)) / den1;   % HtM Hicksian ≈ eps_x
epsH2_imp = (c2h(t0) - e2h(t0)) / den2;   % Savers Hicksian ≈ eps_x

fprintf('Marshallian (uncompensated) elasticity: HtM = %+.4f | Savers = %+.4f\n', ...
        etaM1_imp, etaM2_imp);

fprintf('Impact Hicksian (compensated / substitution) elasticity: HtM = %+.4f | Savers = %+.4f\n', ...
        epsH1_imp, epsH2_imp);

% 12.  Compute energy and non-energy share of nominal consumption drop on impact

% --- Impact level changes (first-order) ---
EE1_imp   = pE_ss*e1_ss * (peh(t0) + e1h(t0));   % Δ(p_e e_1) on impact
EE2_imp   = pE_ss*e2_ss * (peh(t0) + e2h(t0));   % Δ(p_e e_2)
Cnom1_imp = p1_ss*c1_ss * (p1h(t0) + c1h(t0));   % Δ(p_1 c_1)
Cnom2_imp = p2_ss*c2_ss * (p2h(t0) + c2h(t0));   % Δ(p_2 c_2)
XX1_imp   = Cnom1_imp - EE1_imp;                 % Δx_1  (residual)
XX2_imp   = Cnom2_imp - EE2_imp;                 % Δx_2

% --- Shares on impact (of the DROP only) ---
if Cnom1_imp < 0
  ShareEE1_imp = EE1_imp/(-Cnom1_imp);
  ShareXX1_imp = XX1_imp/(-Cnom1_imp);
else
  ShareEE1_imp = NaN; ShareXX1_imp = NaN;  % nominal didn’t drop; “share of drop” N/A
end

if Cnom2_imp < 0
  ShareEE2_imp = EE2_imp/(-Cnom2_imp);
  ShareXX2_imp = XX2_imp/(-Cnom2_imp);
else
  ShareEE2_imp = NaN; ShareXX2_imp = NaN;
end

fprintf('Impact energy share of nominal drop:  HtM %.1f%% | Savers %.1f%%\n', ...
        100*ShareEE1_imp, 100*ShareEE2_imp);
fprintf('Impact non-energy share of nominal drop: HtM %.1f%% | Savers %.1f%%\n', ...
        100*ShareXX1_imp, 100*ShareXX2_imp);

% 13.  Compute the effect on energy and non-energy bills on impact

% Pull % IRFs and convert to FRACTIONS (divide by 100 exactly once)
c1h = oo_.irfs.(['c1_pct_'  shock])/100;   
c2h = oo_.irfs.(['c2_pct_'  shock])/100;
p1h = oo_.irfs.(['p1_pct_'  shock])/100;
p2h = oo_.irfs.(['p2_pct_'  shock])/100;   
e1h = oo_.irfs.(['e1_pct_'  shock])/100;  
e2h = oo_.irfs.(['e2_pct_'  shock])/100;  
peh = oo_.irfs.(['pe_pct_'  shock])/100;   

% Steady states
names = cellstr(M_.endo_names);
getSS = @(v) oo_.steady_state(strcmp(names,v));
pE_ss = getSS('p_e');
e1_ss = getSS('e_1');  e2_ss = getSS('e_2');
p1_ss = getSS('p_1');  p2_ss = getSS('p_2');
c1_ss = getSS('c_1');  c2_ss = getSS('c_2');

% Baseline nominal consumption (p_i*c_i) and non-energy x_i levels
Cnom1_ss = p1_ss*c1_ss;             Cnom2_ss = p2_ss*c2_ss;
EE1_ss   = pE_ss*e1_ss;             EE2_ss   = pE_ss*e2_ss;
x1_ss    = Cnom1_ss - EE1_ss;       x2_ss    = Cnom2_ss - EE2_ss;

EE1_imp   = EE1_ss * (peh(t0) + e1h(t0));
EE2_imp   = EE2_ss * (peh(t0) + e2h(t0));

% Δ(p_i c_i) ≈ (p_i*c_i)_SS * [ Δlog p_i + Δlog c_i ]
Cnom1_imp = Cnom1_ss * (p1h(t0) + c1h(t0));
Cnom2_imp = Cnom2_ss * (p2h(t0) + c2h(t0));

% Δx_i = Δ(p_i c_i) − Δ(p_e e_i)  (value identity)
X1_imp    = Cnom1_imp - EE1_imp;
X2_imp    = Cnom2_imp - EE2_imp;

% ---------- Also express them as % of their OWN baselines ----------
dEE1_pct = 100*(peh(t0) + e1h(t0));                       % %Δ energy bill (HtM)
dEE2_pct = 100*(peh(t0) + e2h(t0));                       % %Δ energy bill (Savers)
dX1_pct  = 100*( X1_imp / max(x1_ss,eps) );               % %Δ x_1  (uses baseline x_1)
dX2_pct  = 100*( X2_imp / max(x2_ss,eps) );               % %Δ x_2

% Print nicely (model units are in the numeraire; signs tell rise/fall)
fprintf('Impact Δ energy bill:   HtM = %+.4e (%.2f%% of baseline) | Savers = %+.4e (%.2f%%)\n', ...
        EE1_imp, dEE1_pct, EE2_imp, dEE2_pct);
fprintf('Impact Δ non-energy x:  HtM = %+.4e (%.2f%% of baseline) | Savers = %+.4e (%.2f%%)\n', ...
        X1_imp,  dX1_pct,  X2_imp,  dX2_pct);

% (Optional) quick check of the identity on impact:
chk1 = Cnom1_imp - (EE1_imp + X1_imp);
chk2 = Cnom2_imp - (EE2_imp + X2_imp);
fprintf('Impact identity gaps (should be ~0): HtM %.2e | Savers %.2e\n', chk1, chk2);

% === Impact index ===
t0 = 1;  % Dynare IRFs: first element is the impact period

% === Impact decomposition: Total = Direct (price) + Indirect (resources) ===
Total1_imp    =  c1_ss * c1h(t0);
Direct1_imp   =  c1_ss * (-p1h(t0));
Indirect1_imp =  c1_ss * ( r1h(t0));

Total2_imp    =  c2_ss * c2h(t0);
Direct2_imp   =  c2_ss * (-p2h(t0));
Indirect2_imp =  c2_ss * ( r2h(t0));

% Identity check (should be ~0 up to rounding)
gap1 = abs(Total1_imp - (Direct1_imp + Indirect1_imp));
gap2 = abs(Total2_imp - (Direct2_imp + Indirect2_imp));
fprintf('Impact gaps: type1 %.2e | type2 %.2e\n', gap1, gap2);

% === Shares on impact ===
% "Share of total change" (can be >100% or negative if Total is small or flips sign)
ShareDir1_imp = Direct1_imp / Total1_imp;
ShareDir2_imp = Direct2_imp / Total2_imp;

% If you specifically want "share of the DROP", use -Total when Total<0, else mark N/A
ShareDir1_drop = NaN; ShareDir2_drop = NaN;
if Total1_imp < 0, ShareDir1_drop = Direct1_imp / (-Total1_imp); end
if Total2_imp < 0, ShareDir2_drop = Direct2_imp / (-Total2_imp); end

% === Print nicely ===
fprintf('IMPACT (levels):\n');
fprintf('  Spenders:  Total=%+.3e | Direct=%+.3e | Indirect=%+.3e\n', ...
        Total1_imp, Direct1_imp, Indirect1_imp);
fprintf('  Savers:    Total=%+.3e | Direct=%+.3e | Indirect=%+.3e\n', ...
        Total2_imp, Direct2_imp, Indirect2_imp);

fprintf('IMPACT direct share of total change: Spenders %+.1f%% | Savers %+.1f%%\n', ...
        100*ShareDir1_imp, 100*ShareDir2_imp);