% --- 1) Run the model ---
dynare model_thesis.mod noclearall
%% === 1x2 IRF plot: consumption vs income ===
shock = 'eta_tau';   % choose the carbon tax shock
H = options_.irf;
tt = 0:H-1;

% --- Pull IRFs (% deviations already defined in the model) ---
c_pct   = oo_.irfs.(['c_pct_' shock]);
c1_pct  = oo_.irfs.(['c1_pct_' shock]);
c2_pct  = oo_.irfs.(['c2_pct_' shock]);

y_pct    = oo_.irfs.(['Y_priv_pct_' shock]);
y_inc_pct = oo_.irfs.(['Y_inc_priv_pct_' shock]);
inc1_pct = oo_.irfs.(['inc1_pct_' shock]);
inc2_pct = oo_.irfs.(['inc2_pct_' shock]);

% --- Figure layout ---
F2 = figure('Color','w','Name','Consumption and Income IRFs');
t = tiledlayout(F2,1,2,'TileSpacing','compact','Padding','compact');

% --- Left panel: Consumption ---
ax1 = nexttile; hold(ax1,'on'); grid(ax1,'on'); box(ax1,'on');
plot(tt, c_pct(1:H),  'LineWidth',1);
plot(tt, c1_pct(1:H), 'LineWidth',1);
plot(tt, c2_pct(1:H), 'LineWidth',1);
yline(ax1,0,'k-','HandleVisibility','off');
title(ax1,'Consumption');
xlabel(ax1,'Quarters'); ylabel(ax1,'% deviation');
legend(ax1,{'Aggregate C','Spenders C_1','Savers C_2'}, ...
    'Location','best','Box','off');
set(ax1,'FontName','Helvetica','FontSize',11,'LineWidth',0.8);

% --- Right panel: Income ---
ax2 = nexttile; hold(ax2,'on'); grid(ax2,'on'); box(ax2,'on');
% plot(tt, y_pct(1:H),    'LineWidth',1);
plot(tt, y_inc_pct(1:H),    'LineWidth',1);
plot(tt, inc1_pct(1:H), 'LineWidth',1);
plot(tt, inc2_pct(1:H), 'LineWidth',1);
yline(ax2,0,'k-','HandleVisibility','off');
title(ax2,'Income');
xlabel(ax2,'Quarters'); ylabel(ax2,'% deviation');
legend(ax2,{ ...
    % 'Aggregate Y', ...
    'Aggregate private income','Spenders income','Savers income'}, ...
    'Location','best','Box','off');
set(ax2,'FontName','Helvetica','FontSize',11,'LineWidth',0.8);
