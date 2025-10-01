% Plot Role of heterogeneity graph


shock = 'eta_tau';             
file_baseline    = 'baseline.mat';
file_no_income   = 'no_income_incidence.mat';
file_no_energy   = 'no_energy_share.mat';

% line styles 
sty_base_htm   = {'-','LineWidth',1.3};
sty_base_saver = {'-','LineWidth',1.3};
sty_no_income  = {':','LineWidth',1.5};
sty_no_energy  = {'-.' ,'LineWidth',1.3};

%% Pull IRFs
S = struct();
S.base   = load(file_baseline,   'oo_','M_','options_');
S.ninc   = load(file_no_income,  'oo_','M_','options_');
S.neng   = load(file_no_energy,  'oo_','M_','options_');

% variable names 
v_c1 = ['c1_pct_' shock];
v_c2 = ['c2_pct_' shock];

must = @(cond,msg) assert(cond, msg);
must(isfield(S.base.oo_.irfs, v_c1), ['Missing ' v_c1 ' in baseline run.']);
must(isfield(S.base.oo_.irfs, v_c2), ['Missing ' v_c2 ' in baseline run.']);
must(isfield(S.ninc.oo_.irfs, v_c1), ['Missing ' v_c1 ' in no-income-incidence run.']);
must(isfield(S.neng.oo_.irfs, v_c1), ['Missing ' v_c1 ' in no-energy-share run.']);

c1_base = S.base.oo_.irfs.(v_c1)(:);   % HTM baseline
c2_base = S.base.oo_.irfs.(v_c2)(:);   % Saver baseline
c1_ninc = S.ninc.oo_.irfs.(v_c1)(:);   % HtM, no income incidence
c1_neng = S.neng.oo_.irfs.(v_c1)(:);   % HtM, no energy-share

% align horizons if different
T = min([numel(c1_base), numel(c2_base), numel(c1_ninc), numel(c1_neng)]);
t = (0:T-1)';   % periods

%% PLOT
figure('Color','w'); hold on; grid on; box on;
plot(t, c2_base(1:T), sty_base_saver{:});                   % Saver baseline
plot(t, c1_base(1:T), sty_base_htm{:});                     % HTM baseline
plot(t, c1_ninc(1:T), sty_no_income{:});                    % Saver, no income-incidence hetero
plot(t, c1_neng(1:T), sty_no_energy{:});                    % Saver, no energy-share hetero

yline(0, '-', 'HandleVisibility','off');  % zero line, keep out of legend
y = [c2_base; c1_base; 0]; y = y(isfinite(y)); if isempty(y), y=0; end
pad = 0.1*(max(y)-min(y)+eps); ylim([min(y)-pad, max(y)+pad]);;

xlabel('Quarters');
ylabel('% deviation');
% title('Consumption');

legend({ ...
        'Saver', ...
        'HtM - baseline', ...
        'HtM - heterogeneity in energy share only', ...
        'HtM - heterogeneity in income incidence only'}, ...
        'Location','best');
