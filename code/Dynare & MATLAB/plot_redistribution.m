% Graph of the role of redistribution carbon market revenues

%% Pull IRFs
shock = 'eta_tau';
file_baseline    = 'baseline.mat';
file_redistribution   = 'redistribution.mat';
file_redistribution_abatement   = 'redistribution_abatement_rigidity.mat';

sty_base  = {'-','LineWidth',1};
sty_redistribution = {'--','LineWidth',1};
sty_redistribution_abatement  = {':','LineWidth',1.5};

S = struct();
S.base   = load(file_baseline,   'oo_','M_','options_');
S.red    = load(file_redistribution,  'oo_','M_','options_');
S.abred  = load(file_redistribution_abatement,  'oo_','M_','options_');

% variable names in oo_.irfs
v_c   = ['c_pct_' shock];       % total consumption
v_c1  = ['c1_pct_' shock];      % HtM consumption
v_co2 = ['co2_pct_' shock];     % emissions
v_y   = ['y_pct_' shock];       % output (demand)
v_e   = ['e2_pct_' shock];       % energy production (back in)
v_i   = ['l_2_' shock];       % investment

% pull IRFs (% deviations)
c_base    = S.base.oo_.irfs.(v_c)(:);
c_red     = S.red.oo_.irfs.(v_c)(:);
c_abred   = S.abred.oo_.irfs.(v_c)(:);

c1_base   = S.base.oo_.irfs.(v_c1)(:);
c1_red    = S.red.oo_.irfs.(v_c1)(:);
c1_abred  = S.abred.oo_.irfs.(v_c1)(:);

co2_base  = S.base.oo_.irfs.(v_co2)(:);
co2_red   = S.red.oo_.irfs.(v_co2)(:);
co2_abred = S.abred.oo_.irfs.(v_co2)(:);

y_base    = S.base.oo_.irfs.(v_y)(:);
y_red     = S.red.oo_.irfs.(v_y)(:);
y_abred   = S.abred.oo_.irfs.(v_y)(:);

e_base    = S.base.oo_.irfs.(v_e)(:);
e_red     = S.red.oo_.irfs.(v_e)(:);
e_abred   = S.abred.oo_.irfs.(v_e)(:);

i_base    = S.base.oo_.irfs.(v_i)(:);
i_red     = S.red.oo_.irfs.(v_i)(:);
i_abred   = S.abred.oo_.irfs.(v_i)(:);

% align horizons if different
T = numel(c_base);
tt = (0:T-1)';

% PLOT
F = figure('Color','w','Name','Redistribution IRFs (3x3)');
t = tiledlayout(F,3,3,'TileSpacing','compact','Padding','compact');

% Aggregate Consumption
ax1 = nexttile; hold(ax1,'on'); grid(ax1,'on'); box(ax1,'on');
plot(ax1, tt, c_base,  sty_base{:},  'DisplayName','Baseline');
plot(ax1, tt, c_red,   sty_redistribution{:},   'DisplayName','Redistribution');
plot(ax1, tt, c_abred, sty_redistribution_abatement{:}, 'DisplayName','Redistribution & abatement rigidity');
yline(ax1,0,'-','HandleVisibility','off');
title(ax1,'Total consumption'); xlabel(ax1,'Quarters'); ylabel(ax1,'% deviation');
set(ax1,'FontName','Helvetica','FontSize',11,'LineWidth',0.8);
auto_ylim(ax1, [c_base; c_red; c_abred]);
legend(ax1,'Location','best','Box','on');   % legend only here

% Hand-to-Mouth Consumption
ax2 = nexttile; hold(ax2,'on'); grid(ax2,'on'); box(ax2,'on');
plot(ax2, tt, c1_base,   sty_base{:},  'HandleVisibility','off');
plot(ax2, tt, c1_red,    sty_redistribution{:},   'HandleVisibility','off');
plot(ax2, tt, c1_abred,  sty_redistribution_abatement{:}, 'HandleVisibility','off');
yline(ax2,0,'-','HandleVisibility','off');
title(ax2,'Hand-to-Mouth consumption'); xlabel(ax2,'Quarters'); ylabel(ax2,'% deviation');
set(ax2,'FontName','Helvetica','FontSize',11,'LineWidth',0.8);
auto_ylim(ax2, [c1_base; c1_red; c1_abred]);

% Aggregate Emissions
ax3 = nexttile; hold(ax3,'on'); grid(ax3,'on'); box(ax3,'on');
plot(ax3, tt, co2_base,   sty_base{:},  'HandleVisibility','off');
plot(ax3, tt, co2_red,    sty_redistribution{:},   'HandleVisibility','off');
plot(ax3, tt, co2_abred,  sty_redistribution_abatement{:}, 'HandleVisibility','off');
yline(ax3,0,'-','HandleVisibility','off');
title(ax3,'Emissions'); xlabel(ax3,'Quarters'); ylabel(ax3,'% deviation');
set(ax3,'FontName','Helvetica','FontSize',11,'LineWidth',0.8);
auto_ylim(ax3, [co2_base; co2_red; co2_abred]);

% Output (Demand)
ax4 = nexttile; hold(ax4,'on'); grid(ax4,'on'); box(ax4,'on');
plot(ax4, tt, y_base,   sty_base{:},  'HandleVisibility','off');
plot(ax4, tt, y_red,    sty_redistribution{:},   'HandleVisibility','off');
plot(ax4, tt, y_abred,  sty_redistribution_abatement{:}, 'HandleVisibility','off');
yline(ax4,0,'-','HandleVisibility','off');
title(ax4,'Output (demand)'); xlabel(ax4,'Quarters'); ylabel(ax4,'% deviation');
set(ax4,'FontName','Helvetica','FontSize',11,'LineWidth',0.8);
auto_ylim(ax4, [y_base; y_red; y_abred]);

% Energy Production
ax5 = nexttile; hold(ax5,'on'); grid(ax5,'on'); box(ax5,'on');
plot(ax5, tt, e_base,   sty_base{:},  'HandleVisibility','off');
plot(ax5, tt, e_red,    sty_redistribution{:},   'HandleVisibility','off');
plot(ax5, tt, e_abred,  sty_redistribution_abatement{:}, 'HandleVisibility','off');
yline(ax5,0,'-','HandleVisibility','off');
title(ax5,'Energy production'); xlabel(ax5,'Quarters'); ylabel(ax5,'% deviation');
set(ax5,'FontName','Helvetica','FontSize',11,'LineWidth',0.8);
auto_ylim(ax5, [e_base; e_red; e_abred]);

% Investment
ax6 = nexttile; hold(ax6,'on'); grid(ax6,'on'); box(ax6,'on');
plot(ax6, tt, i_base,   sty_base{:},  'HandleVisibility','off');
plot(ax6, tt, i_red,    sty_redistribution{:},   'HandleVisibility','off');
plot(ax6, tt, i_abred,  sty_redistribution_abatement{:}, 'HandleVisibility','off');
yline(ax6,0,'-','HandleVisibility','off');
title(ax6,'Investment'); xlabel(ax6,'Quarters'); ylabel(ax6,'% deviation');
set(ax6,'FontName','Helvetica','FontSize',11,'LineWidth',0.8);
auto_ylim(ax6, [i_base; i_red; i_abred]);


%% Function
function auto_ylim(ax, arr)
    ymin = min(arr,[],'omitnan');
    ymax = max(arr,[],'omitnan');
    if isnan(ymin) || isnan(ymax), return; end
    if ymin > 0, ymin = 0; end
    if ymax < 0, ymax = 0; end
    pad = 0.03 * (ymax - ymin + eps);
    ylim(ax, [ymin - pad, ymax + 1.5*pad]);
end
