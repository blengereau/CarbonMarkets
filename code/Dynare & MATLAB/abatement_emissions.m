clearvars; close all; clc;

% --- 1) Run the model ---
dynare model_thesis.mod noclearall

% --- 2) Settings ---
shock = 'eta_tau';            % which shock's IRFs to read
H = 15;                       % fallback horizon
if exist('options_', 'var') && isfield(options_, 'irf') && ~isempty(options_.irf)
    H = options_.irf;
end
tt = 0:H-1;

getSS = @(nm) oo_.steady_state(strcmp(cellstr(M_.endo_names), nm));
getIRF = @(nm) (isfield(oo_.irfs, [nm '_' shock])) .* oo_.irfs.([nm '_' shock]);

% --- 3) Pull series ---
mu_Y = getIRF('mu_Y');  
mu_E = getIRF('mu_E');   
co2_pct    = getIRF('co2_pct');

%% --- Abatement & CO2 IRFs (baseline vs. abatement rigidity) ---
% Baseline series (truncate to H if needed)
muY0 = mu_Y(1:min(H,numel(mu_Y)));
muE0 = mu_E(1:min(H,numel(mu_E)));
co20 = co2_pct(1:min(H,numel(co2_pct)));

% Run abatement-rigidity version of the model
dynare copy_model_thesis_abatement_1.mod

shock = 'eta_tau';            % which shock's IRFs to read
H = 15;                       % fallback horizon
if exist('options_', 'var') && isfield(options_, 'irf') && ~isempty(options_.irf)
    H = options_.irf;
end
tt = 0:H-1;

getSS = @(nm) oo_.steady_state(strcmp(cellstr(M_.endo_names), nm));
getIRF = @(nm) (isfield(oo_.irfs, [nm '_' shock])) .* oo_.irfs.([nm '_' shock]);


muY1 = getIRF('mu_Y');    muY1 = muY1(1:min(H,numel(muY1)));
muE1 = getIRF('mu_E');    muE1 = muE1(1:min(H,numel(muE1)));
co21 = getIRF('co2_pct'); co21 = co21(1:min(H,numel(co21)));

% Plot (1x3): solid = baseline, dashed red = rigidity
F = figure('Color','w','Name','IRFs: Abatement & CO2');
tiledlayout(F,1,3,'TileSpacing','compact','Padding','compact');

names = {'Abatement (Y sector)','Abatement (Energy sector)','CO_2 (aggregate, % dev)'};
A0 = {muY0, muE0, co20}; 
A1 = {muY1, muE1, co21};

for i = 1:3
    ax = nexttile; hold(ax,'on'); grid(ax,'off'); box(ax,'on');
    plot(tt(1:numel(A0{i})), A0{i}, 'LineWidth',1);
    if any(isfinite(A1{i})), plot(tt(1:numel(A1{i})), A1{i}, 'r--','LineWidth',1.2); end
    yline(ax,0,'LineWidth',0.8,'HandleVisibility','off');

    title(ax,names{i});
    if i==3          % only CO2 is in % deviation
        ylabel(ax,'% deviation');
    else
        ylabel(ax,'');  % no label for abatement levels
    end
    xlabel(ax,'Quarters');

    set(ax,'FontName','Helvetica','FontSize',11,'LineWidth',0.8);

    y = [A0{i}(:); A1{i}(:); 0]; y = y(isfinite(y)); if isempty(y), y = 0; end
    pad = 0.05*(max(y)-min(y)+eps); ylim(ax,[min(y)-pad, max(y)+pad]);

    if i==1, legend(ax,{'Baseline','Abatement rigidity'},'Location','best'); end
end