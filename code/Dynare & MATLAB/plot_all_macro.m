clearvars; close all; clc;

%% --- Run model ---
dynare model_thesis.mod noclearall

%% Pull IRF series
shock = 'eta_tau';
H = 15; if exist('options_','var') && isfield(options_,'irf') && ~isempty(options_.irf), H = options_.irf; end
tt = 0:H-1;

vars6   = {'pe_pct','w_pct','rK_pct','e_pct','k_agg_pct','y_pct'};
titles6 = {'Energy price','Wages','Interest rate','Energy production','Aggregate capital','Output'};

S6 = cell(size(vars6));
for i = 1:numel(vars6)
    fld = [vars6{i} '_' shock];
    if isfield(oo_.irfs, fld)
        x = oo_.irfs.(fld)(:); x = x(1:min(H,numel(x))); if numel(x)<H, x(end+1:H,1)=NaN; end
    else
        x = NaN(H,1); warning('IRF %s missing (baseline).', fld);
    end
    S6{i} = x;
end

muY0 = []; if isfield(oo_.irfs, ['mu_Y_' shock]),  muY0 = oo_.irfs.(['mu_Y_'  shock])(:);  end
muE0 = []; if isfield(oo_.irfs, ['mu_E_' shock]),  muE0 = oo_.irfs.(['mu_E_'  shock])(:);  end
co20 = []; if isfield(oo_.irfs, ['co2_pct_' shock]), co20 = oo_.irfs.(['co2_pct_' shock])(:); end

%% Run model with rigidity
dynare copy_model_thesis_abatement_1.mod noclearall

% Pull rigidity IRF  
muY1 = []; if isfield(oo_.irfs, ['mu_Y_' shock]),  muY1 = oo_.irfs.(['mu_Y_'  shock])(:);  end
muE1 = []; if isfield(oo_.irfs, ['mu_E_' shock]),  muE1 = oo_.irfs.(['mu_E_'  shock])(:);  end
co21 = []; if isfield(oo_.irfs, ['co2_pct_' shock]), co21 = oo_.irfs.(['co2_pct_' shock])(:); end
if isempty(muY1), muY1 = NaN(H,1); else, muY1 = muY1(1:min(H,numel(muY1))); if numel(muY1)<H, muY1(end+1:H,1)=NaN; end, end
if isempty(muE1), muE1 = NaN(H,1); else, muE1 = muE1(1:min(H,numel(muE1))); if numel(muE1)<H, muE1(end+1:H,1)=NaN; end, end
if isempty(co21), co21 = NaN(H,1); else, co21 = co21(1:min(H,numel(co21))); if numel(co21)<H, co21(end+1:H,1)=NaN; end, end

Sr6 = cell(size(vars6));
for i = 1:numel(vars6)
    fld = [vars6{i} '_' shock];
    if isfield(oo_.irfs, fld)
        x = oo_.irfs.(fld)(:); x = x(1:min(H,numel(x))); if numel(x)<H, x(end+1:H,1)=NaN; end
    else
        x = NaN(H,1); warning('IRF %s missing (rigidity).', fld);
    end
    Sr6{i} = x;
end

%% --- Plot % deviation (price, wage, rate ; energy prod, agg capital, output)
F1 = figure('Color','w','Name','(1) IRFs: 2x3 baseline');
t1 = tiledlayout(F1,2,3,'TileSpacing','compact','Padding','compact');
for i = 1:6
    ax = nexttile; hold(ax,'on'); grid(ax,'off'); box(ax,'on');
    xi = S6{i}; 
    plot(tt, xi, 'LineWidth',1); 
    yline(ax,0,'LineWidth',0.8,'HandleVisibility','off');

    ymin = min(xi,[],'omitnan'); if isempty(ymin) || isinf(ymin), ymin=0; end
    ymax = max(xi,[],'omitnan'); if isempty(ymax) || isinf(ymax), ymax=0; end
    if ymin>0, ymin=0; end
    if ymax<0, ymax=0; end
    pad = 0.1*(ymax - ymin + eps); 
    ylim(ax,[ymin - pad, ymax + pad]);

    title(ax, titles6{i}); 
    ylabel(ax,'% deviation'); 
    if i>3, xlabel(ax,'Quarters'); end
    set(ax,'FontName','Helvetica','FontSize',11,'LineWidth',0.8);
end

%% --- Plot mu_Y level, mu_E level, CO2 %, with rigidity overlay
F2 = figure('Color','w','Name','(2) Abatement (levels) & CO2 (%): baseline vs rigidity');
t2 = tiledlayout(F2,1,1,'TileSpacing','compact','Padding','compact');

% Panel 3: CO2 (% deviation)
ax = nexttile; hold(ax,'on'); grid(ax,'off'); box(ax,'on');
plot(tt, co20, 'LineWidth',1); plot(tt, co21, 'r--','LineWidth',1.2);
yline(ax,0,'LineWidth',0.8,'HandleVisibility','off');
y = [co20; co21; 0]; y = y(isfinite(y)); if isempty(y), y=0; end
pad = 0.1*(max(y)-min(y)+eps); ylim(ax,[min(y)-pad, max(y)+pad]);
title(ax,'Emissions'); ylabel(ax,'% deviation'); xlabel(ax,'Quarters');
set(ax,'FontName','Helvetica','FontSize',11,'LineWidth',0.8);
legend(ax,{'Baseline','Abatement rigidity'},'Location','best');

%% --- Plot baseline vs rigidity
F3 = figure('Color','w','Name','(3) IRFs: 2x3 baseline vs rigidity');
t3 = tiledlayout(F3,2,3,'TileSpacing','compact','Padding','compact');
for i = 1:6
    ax = nexttile; hold(ax,'on'); grid(ax,'off'); box(ax,'on');

    % Correct cell indexing with {}
    xi = S6{i};          % baseline series (vector)
    xj = Sr6{i};         % rigidity series (vector)

    plot(tt, xi, 'LineWidth',1);
    if any(isfinite(xj)), plot(tt, xj, 'r--', 'LineWidth',1.2); end
    yline(ax,0,'LineWidth',0.8,'HandleVisibility','off');

    y = [xi; xj; 0]; y = y(isfinite(y)); if isempty(y), y = 0; end
    pad = 0.1*(max(y)-min(y)+eps); ylim(ax, [min(y)-pad, max(y)+pad]);

    title(ax, titles6{i}); ylabel(ax,'% deviation'); if i>3, xlabel(ax,'Quarters'); end
    set(ax,'FontName','Helvetica','FontSize',11,'LineWidth',0.8);

    if i==1, legend(ax,{'Baseline','Abatement rigidity'},'Location','best'); end
end


