clear; clc;

% 1) Parse the .mod file BUT keep workspace variables
dynare model_thesis noclearall nolog;

% 2) Calibrate tau_bar to hit the 1.5% SS share
tau_bar_val = calibrate_tau_bar(0.0087);
set_param_value('tau_bar', tau_bar_val);

% 3) Re-compute SS and checks with the calibrated tau_bar
steady; resid; check;

% 4) Set shocks (e.g., 1% price shock if you like)
% shocks; var eta_tau; stderr log(1.01); end;

% % 5) Simulate
% stoch_simul(order=1, pruning, irf=12, tex) ...
%     y k_Y l_Y e_Y k_agg l_E k_E i_agg g c_pct c1_pct c2_pct y_pct Y_pct ...
%     Y_priv_pct Y_inc_priv_pct w_pct co2_pct p1_pct p2_pct pe_pct inc1_pct ...
%     inc2_pct e_pct rK_pct co2_agg s tau_wedge iota_model;
