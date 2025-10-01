% --- 1) Run the model ---
dynare model_thesis noclearall

%% === Plot consumption & income ===
shock = 'eta_tau';   % choose the carbon tax shock
H = options_.irf;
tt = 0:H-1;

% Pull IRFs
c_pct   = oo_.irfs.(['c_pct_' shock]);
c1_pct  = oo_.irfs.(['c1_pct_' shock]);
c2_pct  = oo_.irfs.(['c2_pct_' shock]);

p1_pct  = oo_.irfs.(['p1_pct_' shock]);
p2_pct  = oo_.irfs.(['p2_pct_' shock]);

y_inc_pct = oo_.irfs.(['Y_inc_priv_pct_' shock]);
inc1_pct = oo_.irfs.(['inc1_pct_' shock]);
inc2_pct = oo_.irfs.(['inc2_pct_' shock]);

x_pct  = oo_.irfs.(['x_pct_'  shock]);   % aggregate consumption goods
x1_pct = oo_.irfs.(['x1_pct_' shock]);   % HtM goods
x2_pct = oo_.irfs.(['x2_pct_' shock]);   % Savers goods

ec_pct = oo_.irfs.(['ec_pct_' shock]);   % aggregate energy consumption
e1_pct = oo_.irfs.(['e1_pct_' shock]);   % HtM energy
e2_pct = oo_.irfs.(['e2_pct_' shock]);   % Savers energy

% --- Figure layout ---
F = figure('Color','w','Name','Consumption and Income IRFs');
% t = tiledlayout(F,2,3,'TileSpacing','compact','Padding','compact');
t = tiledlayout(F,2,6,'TileSpacing','compact','Padding','compact');

% Colors
col_htm = [1 0 0];          
col_sav = [1 0.84 0];     
col_agg = [0 0 1]; 

% Consumption
% ax1 = nexttile; hold(ax1,'on'); grid(ax1,'on'); box(ax1,'on');
ax1 = nexttile(1,[1 2]); hold(ax1,'on'); grid(ax1,'on'); box(ax1,'on');

h_agg = plot(tt, c_pct, 'LineWidth',1);
h_htm = plot(tt, c1_pct, 'LineWidth',1);
h_sav = plot(tt, c2_pct, 'LineWidth',1);

yline(ax1,0,'k-','HandleVisibility','off');
title(ax1,'Consumption'); xlabel(ax1,'Quarters'); ylabel(ax1,'% deviation');
set(ax1,'FontName','Helvetica','FontSize',11,'LineWidth',0.8);
% set_ylim0(ax1, [c2_pct; c1_pct; c_pct]);

y_all = [c1_pct(:); c2_pct(:); c_pct(:)];
ymin = min(y_all,[],'omitnan'); ymax = max(y_all,[],'omitnan');
if ymin > 0, ymin = 0; end
if ymax < 0, ymax = 0; end
pad = 0.03 * (ymax - ymin + eps);
ylim(ax1,[ymin - pad, ymax + 1.5*pad]); % give more space above 0

% Legend
legend(ax1, [h_htm h_sav h_agg], ...
    {'Hand-to-Mouth','Savers','Total'}, ...
    'Location','best','Box','on');  % Box ON and anchored to subplot

% Price Index
ax2 = nexttile(3,[1 2]); hold(ax2,'on'); grid(ax2,'on'); box(ax2,'on');

xNAN = NaN(size(tt));
plot(tt, xNAN, 'LineWidth',1); 
plot(tt, p1_pct, 'LineWidth',1);   % HtM price index
plot(tt, p2_pct, 'LineWidth',1);   % Savers' price index

yline(ax2,0,'k-','HandleVisibility','off');
title(ax2,'Price indices'); 
xlabel(ax2,'Quarters'); ylabel(ax2,'% deviation');
set(ax2,'FontName','Helvetica','FontSize',11,'LineWidth',0.8);
% set_ylim0(ax2, [p2_pct; p1_pct]);

y_all = [p1_pct(:); p2_pct(:)];
ymin = min(y_all,[],'omitnan'); ymax = max(y_all,[],'omitnan');
if ymin > 0, ymin = 0; end
if ymax < 0, ymax = 0; end
pad = 0.03 * (ymax - ymin + eps);
ylim(ax2,[ymin - 1.5*pad, ymax + pad]); % more room below 0

% Income
% ax3 = nexttile; hold(ax3,'on'); grid(ax3,'on'); box(ax3,'on');
ax3 = nexttile(5,[1 2]); hold(ax3,'on'); grid(ax3,'on'); box(ax3,'on');

plot(tt, y_inc_pct, 'LineWidth',1);  % Aggregate private income
plot(tt, inc1_pct, 'LineWidth',1);  % HtM income
plot(tt, inc2_pct, 'LineWidth',1);  % Savers income

yline(ax3,0,'k-','HandleVisibility','off');
title(ax3,'Income'); xlabel(ax3,'Quarters'); ylabel(ax3,'% deviation');
set(ax3,'FontName','Helvetica','FontSize',11,'LineWidth',0.8);

y_all = [c1_pct(:); c2_pct(:); c_pct(:)];
ymin = min(y_all,[],'omitnan'); ymax = max(y_all,[],'omitnan');
if ymin > 0, ymin = 0; end
if ymax < 0, ymax = 0; end
pad = 0.03 * (ymax - ymin + eps);
ylim(ax3,[ymin - pad, ymax + 1.5*pad]); % give more space above 0


   
% Goods consumption
% ax4 = nexttile(5); hold(ax4,'on'); grid(ax4,'on'); box(ax4,'on');
ax4 = nexttile(8,[1 2]);  % goods consumption (x_*)
hold(ax4,'on'); grid(ax4,'on'); box(ax4,'on');

hx_agg = plot(tt, x_pct,  'LineWidth',1);
hx_htm = plot(tt, x1_pct, 'LineWidth',1);
hx_sav = plot(tt, x2_pct, 'LineWidth',1);
yline(ax4,0,'k-','HandleVisibility','off');
title(ax4,'Goods consumption'); xlabel(ax4,'Quarters'); ylabel(ax4,'% deviation');
set(ax4,'FontName','Helvetica','FontSize',11,'LineWidth',0.8);
y_all = [x1_pct(:); x2_pct(:); x_pct(:)];
ymin = min(y_all,[],'omitnan'); ymax = max(y_all,[],'omitnan');
if ymin > 0, ymin = 0; end
if ymax < 0, ymax = 0; end
pad = 0.03 * (ymax - ymin + eps);
ylim(ax4,[ymin - pad, ymax + 1.5*pad]);

% Energy consumption
ax5 = nexttile(10,[1 2]); % energy consumption (e_*)
hold(ax5,'on'); grid(ax5,'on'); box(ax5,'on');

he_agg = plot(tt, ec_pct, 'LineWidth',1);
he_htm = plot(tt, e1_pct, 'LineWidth',1);
he_sav = plot(tt, e2_pct, 'LineWidth',1);
yline(ax5,0,'k-','HandleVisibility','off');
title(ax5,'Energy consumption'); xlabel(ax5,'Quarters'); ylabel(ax5,'% deviation');
set(ax5,'FontName','Helvetica','FontSize',11,'LineWidth',0.8);
y_all = [e1_pct(:); e2_pct(:); ec_pct(:)];
ymin = min(y_all,[],'omitnan'); ymax = max(y_all,[],'omitnan');
if ymin > 0, ymin = 0; end
if ymax < 0, ymax = 0; end
pad = 0.03 * (ymax - ymin + eps);
ylim(ax5,[ymin - pad, ymax + 1.5*pad]);