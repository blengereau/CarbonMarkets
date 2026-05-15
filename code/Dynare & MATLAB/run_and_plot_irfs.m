% run_and_plot_irfs.m
clear; close all;

% 1) Run baseline
dynare model_thesis.mod noclearall
save('baseline.mat','oo_','M_','options_');

% 2) Run “no income-incidence heterogeneity”
dynare copy_model_thesis_income_incidence_0.mod noclearall
save('no_income_incidence.mat','oo_','M_','options_');

% 3) Run “no energy-share heterogeneity”
dynare copy_model_thesis_share_hete_0.mod noclearall
save('no_energy_share.mat','oo_','M_','options_');

% 4) Run baseline with redistribution
dynare copy_model_thesis_redistribution.mod noclearall
save('redistribution.mat','oo_','M_','options_');

% 3) Run “abatement with redistribution”
dynare copy_model_thesis_abatement_redis.mod noclearall
save('redistribution_abatement_rigidity.mat','oo_','M_','options_');

% 4) Plot 
plot_heterogeneity;
