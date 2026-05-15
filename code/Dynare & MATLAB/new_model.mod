
// Code by Ghassane Benmir for a Two Country Production Network
// This code is done with Dynare 5.5 under Matlab 2023a. (It is also compatible with recent versions of Dynare and Matlab)
// Please cite the author when using.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  0. Housekeeping (close all graphic windows) and declaring options:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
warning off;
@#define policy = 1
            // 0: Laissez-faire
            // 1: Carbon Policy (fixed tax)
            // 2: Fixed cap


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Endogenous variables: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

var 

    mu_Y $\mu_t$ (long_name='Abatement')  
    mu_E $\mu_{E,t}$ (long_name='Abatement Energy producer')
    tau $\tau_t$ (long_name='Carbon Tax')  
    co2_agg ${CO_2}_t$ (long_name='Aggregate emisisons')
    co2_Y ${CO_2}_t^Y$ (long_name='Emissions')  
    co2_E ${CO_2}_t^E$ (long_name='Emissions energy sector')
    s     $$ (long_name='Carbon concentration')
    c     $c_t$ (long_name='Consumption')  
    c_1   $c_t^{1}$ (long_name='Consumption of spenders')
    c_2   $c_t^{2}$ (long_name='Consumption of savers')
    x_1   $x_{1,t}$ (long_name='Goods consumption of spenders') 
    x_2   $x_{2,t}$ (long_name='Goods consumption of savers') 
    e_1   $e_{1,t}$ (long_name='Energy consumption of spenders') 
    e_2   $e_{2,t}$ (long_name='Energy consumption of savers') 
    l     $l_t$ (long_name='Labour')
    l_1   $l_t^{1}$ (long_name='Labour of spenders')
    l_2   $l_t^{2}$ (long_name='Labour of savers')
    l_Y   $l_t^{Y}$ (long_name='Labor demand in production sector')   
    l_E   $l_t^{E}$ (long_name='Labor demand in energy production sector')
    k     $k_t$ (long_name='Capital')
    k_agg $k_t^{T}$ (long_name='Aggregate Capital')
    k_Y   $k_t^{Y}$ (long_name='Capital demand in production sector')
    k_E   $k_t^{E}$ (long_name='Capital demand in energy sector')
    e     $e_t$ (long_name='Energy production')
    e_Y   $e_t^{Y}$ (long_name='Energy demand in energy sector')
    i     $i_t$ (long_name='Investment')
    i_agg $i_t^{T}$ (long_name='Aggregate Investment')
    d     $d_t$ (long_name='Dividends')
    g     $g_t$ (long_name='Government Spending')
    t     $t_t$ (long_name='Transfers')
    t_1   $t_t^{1}$ (long_name='Transfers of spenders')
    t_2   $t_t^{2}$ (long_name='Transfers of savers')
    w     $w_t$ (long_name='Wages')
    rK    $r_t^{K}$ (long_name='Capital returns')
    y     $y_t$ (long_name='Output') 
    rF    $r_t^{F}$ (long_name='Risk free rate') 
    b_h   $b_t^{h}$ (long_name='Bonds demand')
    b_g   $b_t^{g}$ (long_name='Bonds supply')
    q     $q_t$ (long_name='The lagrangian with respect to capital and investment')
    lambda_1 $\lambda_t^{1}$ (long_name='The lagrangian with respect to consumption and savings (bonds)') 
    lambda_2 $\lambda_t^{2}$ (long_name='The lagrangian with respect to consumption and savings (bonds)')
    varrho_Y ${varrho}_t^Y$ (long_name='Marginal cost (input shadow price)')
    mc_Y  $mc_t^{Y}$ (long_name='Total marginal cot of production')
    varrho_E ${varrho}_t^E$ (long_name='Marginal cost Energy (input shadow price)')
    p_e   $p_e$ (long_name='Energy price')
    p_1   $p_{c,1}$ (long_name='Price index of Spenders consumption bundle')
    p_2   $p_{c,2}$ (long_name='Price index of Savers consumption bundle')
    e_b   $e_t^{b}$ (long_name='The shock variable within the law of motion to consumption (preference shock)')
    e_a   $e_t^{a}$ (long_name='The shock variable within the law of motion to output')
    e_e   $e_t^{e}$ (long_name='The shock variable within the law of motion to energy')
    e_t   $e_t$ (long_name='The shock variable within the law of motion to government transfers')
    % e_c   $e_t^{c}$ (long_name='The shock variable within the law of motion to the carbon tax')
;  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. Exogenous variables: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

varexo 
    eta_b  $eta_t^{b}$ (long_name='The exogenous consumption shock')
    eta_e  $eta_t^{e}$ (long_name='The exogenous consumption shock energy')
    eta_a  $eta_t^{a}$ (long_name='The exogenous output shock')
    eta_t  $eta_t^{t}$ (long_name='The exogenous transfers spending shock')
    eta_tau  $eta_t^{c}$ (long_name='The exogenous carbon tax shock')
;                          

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. Parameters: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

parameters 
    tau_bar $\bar{\tau}$ (long_name='Steady state carbon tax')
    theta_1 $\theta_1$ (long_name='Abatement cost param')
    theta_2 $\theta_2$ (long_name='Abatement cost param')
    phi_Y   $\phi_Y$ (long_name='Emission intensity non-energy sector')
    phi_E   $\phi_E$ (long_name='Emission intensity energy sector')
    beta    $\beta$ (long_name='Time preference')
    sigma   $\sigma$ (long_name='Relative risk aversion')
    alpha   $\alpha$ (long_name='Share of capital income')
    nu      $\nu$ (long_name='Share of energy income')
    alpha_E $\alpha$ (long_name='Share of capital income - energy sector')
    delta   $\delta$ (long_name='Depreciation capital rate')
    chi     $\chi$ (long_name='Level of disutility')
    phi     $\phi$ (long_name='Frisch elasticity of labour disutility')
    l0      $l_0$ (long_name='Hours worked')
    b0      $b_0$ (long_name='Fixed supply of bonds')
    t0      $t_0$ (long_name='Gov spending as % of GDP')
    sdev_a  $sdev_a$ (long_name='Standard deviation of the exogenous shock to income')
    sig_a   $\sigma_a$ (long_name='The sign of the shock to income')
    rho_a   $\rho_a$ (long_name='The shock persistence level for income shock')
    sdev_e  $sdev_a$ (long_name='Standard deviation of the exogenous shock to income')
    sig_e   $\sigma_a$ (long_name='The sign of the shock to income')
    rho_e   $\rho_a$ (long_name='The shock persistence level for income shock')
    sdev_b  $sdev_b$ (long_name='Standard deviation of the exogenous shock to preferences')
    sig_b   $\sigma_b$ (long_name='The sign of the shock to preferences')
    rho_b   $\rho_b$ (long_name='The shock persistence level for preference shock')
    sdev_t  $sdev_t$ (long_name='Standard deviation of the exogenous shock to preferences')
    rho_t   $\rho_t$ (long_name='The shock persistence level for transfers shock')
    sig_t   $\sigma_t$ (long_name='The sign of the shock to transfers')
    omega_s $\omega_s$ (long_name='Fraction of spenders')
    tau_i   $\tau_i$ (long_name='Income tax rate')
    phi_1   $\phi_1$ (long_name='Emissions decay parameter')
    phi_0   $\phi_0$ (long_name='Emissions staying in atmosphere')
    trans_adjust $tr_{adjust}$ (long_name='Transfer adjustment parameter')
    transfer_fix $tr_{fix}$ (long_name='Permanent increase in transfers')
    a_1x    $a_{1,x}$ (long_name='Spenders distribution parameter for goods')
    a_1e    $a_{1,e}$ (long_name='Spenders distribution parameter for energy')
    a_2x    $a_{2,x}$ (long_name='Savers distribution parameter for goods')
    a_2e    $a_{2,e}$ (long_name='Savers distribution parameter for energy')
    eps_x   $\eps_x$ (long_name='Elasticity of substitution energy/goods')
    eps_p   $\eps_y$ (long_name='Elasticity of substitution between the intermediate goods')
    sdev_tau $sdev_a$ (long_name='Standard deviation of the exogenous shock to the carbon tax')
    sig_tau $\sigma_a$ (long_name='The sign of the shock to the carbon tax')
    rho_tau $\rho_a$ (long_name='The shock persistence level for the carbon shock')
    gamma   $\gamma$ (long_name='Climate damage parameter')
    mu     $\mu$ (long_name='Share of carbon revenues to Spenders')
;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Third (bis) calibration for all the params 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Macro Params
beta        = 0.99;   % Time preference
delta       = .0275;    % Depreciation capital rate (2.5% quarterly) 
alpha       = .275;     % Share of capital income
alpha_E     = .66;      % Share of capital income in the energy sector
nu          = .085;     % Share of energy income
sigma       = 0.99;     % Coefficient of relative risk-aversion        
chi         = 8;        % Level of disutility (we calibrate it such that we have l=1/3 at the Steady state)
phi         = 1;        % Frisch elasticity of labour disutility
omega_s     = .25;       % of spenders (Hand to mouth) in the EU (30%)
tau_i       = 0;      % % of income tax  
l0          = 1/3;      % Hours worked (1/3 of 24h a day)
b0          = -.92;     % Fixed supply of bonds
t0          = .075;     % Gov spending as % of GDP
trans_adjust = 0;       % Transfer adjustment
transfer_fix = 0;       % Permanent increase in transfers (make sure it is 0 for the case of the shock!!)
a_1e        = .099;     % Spenders' distribution parameter for energy
a_1x        = 1 - a_1e; % Spenders' distribution parameter for goods
a_2e        = .068;     % Savers' distribution parameter for energy
a_2x        = 1 - a_2e; % Savers' distribution parameter for energy
eps_x       = 0.2;      % Elasticity of substitution energy/goods
eps_p       = 6;        % Elasticity  of substitution between goods
tau_bar     = 0.039;    % Steady-state carbon tax
phi_1       = 1 - 0.9994; % Emissions decay parameter 
phi_0       = 0.5359;   % Emissions staying in atmosphere
gamma       = 5.3*10^(-5); % Climate damage parameter
mu         = 1;        % Share of carbon revenues to spenders

%% Environmental Params
% @#if policy == 0
% carbon_price = 0.01;
% @#elseif policy == 1 || policy == 2
% carbon_price = 0.05;
% @#endif
theta_1 = .1;
theta_2 = 2.6;
phi_Y = 2;
phi_E = 0.3;

%% Shocks Params
sdev_a     = .01;      % Standard deviation of the exogenous shock to income
rho_a      = .9;       % The shock persistence level
sig_a      = 1;        % The sign of the shock (1 if positive, -1 if negative)
sdev_e     = .01;      % Standard deviation of the exogenous shock to income
rho_e      = .9;       % The shock persistence level
sig_e     = 1;         % The sign of the shock (1 if positive, -1 if negative)
sdev_b     = .01;      % Standard deviation of the exogenous shock to preferences
rho_b      = .9;       % The shock persistence level
sig_b      = 1;        % The sign of the shock (1 if positive, -1 if negative)
sdev_t     = .8;       % Standard deviation of the exogenous shock to preferences
rho_t      = .9;       % The shock persistence level
sig_t      = 1;        % The sign of the shock (1 if positive, -1 if negative)
sdev_tau   = .01;      % Standard deviation of the exogenous shock to the carbon tax
sig_tau    = 1;        % The sign of the shock (1 if positive, -1 if negative)
rho_tau    = 0.85;     % Persistence of the carbon tax shock


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fourth declare your model (i.e. the FOCs and Equations of the model) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

model;  

%%%% Households Side (Spenders) %%%% 

[name='Spenders: Marginal utility of consumption']
lambda_1  = c_1^(-sigma)/p_1;

[name='Spenders: Marginal disutility of labour']
(1-tau_i)*w = chi*(l_1^phi)/lambda_1;

[name='Spenders: Budget constraint']
p_1*c_1 = (1-tau_i)*w*l_1 + t_1 ; 

[name='Spenders: Demand function for goods']
x_1 = a_1x * (1/p_1)^(-eps_x) * c_1;

[name='Spenders: Demand function for energy']
e_1 = a_1e * (p_e/p_1)^(-eps_x) * c_1;

[name='Spenders: Price of consumption bundle']
p_1 = (a_1x + a_1e * (p_e)^(1-eps_x))^(1/(1-eps_x));

%%%% Households Side (Savers) %%%%  

[name='Savers: Bonds Euler Equation']
lambda_2 = e_b(+1)*beta*rF*lambda_2(+1)/e_b;

[name='Savers: Marginal utility of consumption']
lambda_2 = c_2^(-sigma)/p_2;

[name='Savers: Marginal disutility of labour']
(1-tau_i)*w = chi*l_2^phi/lambda_2;

[name='FOC w.r.t k: Enveloppe Condition']
q = e_b(+1)*beta*(lambda_2(+1)*rK(+1) + q(+1)*(1-delta));

[name='FOC w.r.t i']
lambda_2 = q ;

[name='Capital law of motion']
k = (1-delta)*k(-1) + i;

[name='Savers: Consumer Budget Constraint']
b_h + p_2*c_2 + i = rF*b_h(-1) + rK*k(-1) + (1-tau_i)*w*l_2 + d/(1-omega_s) + t_2; 

[name='Savers: Demand function for goods']
x_2 = a_2x * (1/p_2)^(-eps_x) * c_2;

[name='Savers: Demand function for energy']
e_2 = a_2e * (p_e/p_2)^(-eps_x) * c_2;

[name='Savers: Price of consumption bundle']
p_2 = (a_2x + a_2e * (p_e)^(1-eps_x))^(1/(1-eps_x));

%%%% Different Aggregation %%%% 

[name='Aggreg of Consumption']
c = omega_s*c_1 + (1-omega_s)*c_2;

[name='Aggreg Labour']
l = omega_s*l_1 + (1-omega_s)*l_2;

[name='Agg Capital']
k_agg = (1-omega_s)*k ; 

[name='Agg Capital Law']
i_agg = k_agg - (1-delta)*k_agg(-1) ; 

[name='Agg emissions']
co2_agg = co2_Y + co2_E;

[name='Carbon concentration']
s = (1-phi_1)*s(-1) + phi_0*co2_agg;

%%%% Non-energy producers %%%% 

[name='Production function (Cobb-Douglas)']
% y = e_a*k_Y(-1)^(alpha)*e^(nu)*l_Y^(1-alpha-nu);
y = e_a*exp(-gamma*s)*k_Y^(alpha)*e_Y^(nu)*l_Y^(1-alpha-nu);

[name='FOC w.r.t l']
w = varrho_Y*(1-alpha-nu)*y/l_Y;

[name='FOC w.r.t k(-1)']
% rK = mc_Y*alpha*y/k_Y(-1);
rK = varrho_Y*alpha*y/k_Y;

[name='FOC w.r.t e']
p_e = varrho_Y*nu*y/e_Y;

[name='Marginal cost of production (FOC w.r.t y)']
varrho_Y = mc_Y - (1-mu_Y)*phi_Y*tau - theta_1*mu_Y^theta_2; % add markup of 17%

[name='Total marginal cost of production']
mc_Y = (eps_p-1)/eps_p;

[name='FOC w.r.t mu']
tau = theta_1*theta_2*mu_Y^(theta_2-1)/phi_Y;

[name='CO2 emissions of non-energy producers']
co2_Y = (1-mu_Y)*phi_Y*y;

[name='Dividends']
d = y*(1-mc_Y);

%%%% Energy producers %%%%

[name='Energy production function']
% e = e_e*k_E(-1)^(alpha)*l_E^(1-alpha);
e = e_e*k_E^(alpha_E)*l_E^(1-alpha_E);

[name='Energy prod.: FOC w.r.t l']
w = varrho_E*(1-alpha_E)*e/l_E;

[name='Energy prod.: FOC w.r.t k(-1)']
% rK = mc_E*alpha*e/k_E(-1);
rK = varrho_E*alpha_E*e/k_E;

[name='Marginal cost of energy production (FOC w.r.t e)']
varrho_E = p_e - (1-mu_E)*phi_E*tau - theta_1*mu_E^theta_2;

[name='FOC w.r.t mu_E']
tau = theta_1*theta_2*mu_E^(theta_2-1)/phi_E;

[name='CO2 emissions of energy producers']
co2_E = (1-mu_E)*phi_E*e;


%%%% Government Spending %%%% 

% [name='Transfers as % of GDP']
% t = t0*y*e_t; 
[name='Transfers as % of Carbon revenues']
t = tau*co2_agg*e_t; 

[name='Gov Budget Constraint']
% b_g + g = rF*b_g(-1) + tau_i*w*l + tau*co2_agg - t; 
rF*b_g(-1) + g + t = b_g + tau_i*w*l + tau*co2_agg;

[name='Carbon tax rule']
tau = (1 - rho_tau)*tau_bar + rho_tau*tau(-1) + sig_tau*eta_tau;

[name='Total Transfers']
t = omega_s*t_1 + (1-omega_s)*t_2;

[name='Transfers to Spenders']
t_1 = (mu/omega_s)*tau*co2_agg;

% [name='Transfers to Savers']
% t_2 = ((1-mu)/(1-omega_s))*tau*co2_agg;

%%%% Environmental Regulator %%%% 
@#if policy == 0
tau = 0.0001;
@#elseif policy == 1
% tau = 0.05;
@#elseif policy == 2
co2_agg = steady_state(co2_agg);
@#endif


%%%% Market Clearing conditions (to close the model) 

[name='Supply of bonds']
b_g = b0*y ;

[name='No deficit Spending']
(1-omega_s)*b_h + b_g = 0;

[name='Capital market clearing']
k_Y + k_E = k_agg(-1);

[name='Labor market clearing']
l_Y + l_E = l;

[name='Energy market clearing']
e = omega_s*e_1 + (1-omega_s)*e_2 + e_Y;



%%%% Exogenous shock law of motion %%%% 

[name='Law of motion of the output shock']
log(e_a) = rho_a*log(e_a(-1)) +sig_a*eta_a;

[name='Law of motion of the output shock energy']
log(e_e) = rho_e*log(e_e(-1)) +sig_e*eta_e;

[name='Law of motion of the preference shock']
log(e_b) = rho_b*log(e_b(-1)) +sig_b*eta_b;

[name='Law of motion of the transfer shock']
log(e_t) = rho_t*log(e_t(-1)) +sig_t*eta_t;

% [name='Law of motion of the carbon tax shock']
% e_c = rho_c*e_c(-1) + sig_c*eta_c;

end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fifth declare your steady state 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

steady_state_model;

e_a = 1 ;
e_e = 1 ;
e_b = 1 ;
e_t = 1 ;
% mc_Y = (eps_p - 1)/eps_p; % Markup level in the production sector

%% Using a solver
@#if policy == 0 || policy == 1 || policy == 2
[lambda_2,rF,c_2,w,l_2,q,rK,i,b_h,t_2,y,k,g,lambda_1,c_1,l_1,t_1,c,t,l,b_g,k_agg,i_agg,tau,mu_Y,co2_Y,varrho_Y, mu_E, co2_E, varrho_E, k_Y, k_E, l_Y, l_E, e, co2_agg, p_e, x_1, e_1, p_1, x_2, e_2, p_2, e_Y, d, mc_Y, s] = get_steady_state_new_model(beta,sigma,chi,phi,delta,alpha,alpha_E,nu,t0,b0,tau_i,omega_s,trans_adjust,transfer_fix,theta_1,theta_2,phi_Y,phi_E,tau_bar,eps_x,a_1x,a_1e,a_2x,a_2e, eps_p, phi_1, phi_0, gamma, mu);
@#endif
   

end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sixth conduct all the checks (Residual, Blanchard-Khan ...) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

steady;     % Provides you with the steady states values
resid;      % Allow you to check if you computed correctly the steady states equations
check;      % Checks the Blanchard Khan conditions to see if the model is well set

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Seventh declare the shocks 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

shocks;

% var eta_a; stderr sdev_a ; 
% var eta_b; stderr sdev_b ; 
% var eta_t; stderr sdev_t ; 
var eta_tau; stderr sdev_tau ; 

end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Eight performe the simulation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

stoch_simul (order=1, pruning, irf=30, tex) c c_1 c_2 y co2_agg g t b_g b_h rF tau; 
