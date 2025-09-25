% Model of The Uneven Macroeconomic Effects of Carbon Markets
% Benjamin Lengereau, September 2025

%% Declare model specifications
close all;
warning off;

@#define inelastic_labor = 0   // set to 0 for endogenous labour
@#define dividends = 1         // set to 0 for no dividends 
@#define damage_function =  1  // set to 0 for no damage function
@#define abatment_rigidity = 0 // set to 0 for no abatement rigidity
@#define bonds = 0             // set to 0 for no bonds in the model
@#define income_incidence = 1  // set to 0 for no heterogneity in income incidence
@#define energy_share_heterogeneity = 1 // set to 1 for energy share heterogeneity

% 1. Endogeneous variables

var 
    mu_Y $\mu_t$ (long_name='Abatement')  
    mu_E $\mu_{E,t}$ (long_name='Abatement Energy producer')
    tau $\tau_t$ (long_name='Carbon Tax')  
    co2_agg ${CO_2}_t$ (long_name='Aggregate emisisons')
    co2_Y ${CO_2}_t^Y$ (long_name='Emissions')  
    co2_E ${CO_2}_t^E$ (long_name='Emissions energy sector')
    s     $$ (long_name='Carbon concentration')
    c     $c_t$ (long_name='Aggregate consumption')  
    c_1   $c_t^{1}$ (long_name='Consumption of spenders')
    c_2   $c_t^{2}$ (long_name='Consumption of savers')
    x     $x_t$ (long_name='Aggregate non-energy consumption')
    x_1   $x_{1,t}$ (long_name='Goods consumption of spenders') 
    x_2   $x_{2,t}$ (long_name='Goods consumption of savers') 
    e_c   $e_t^C$ (long_name='Aggregate energy consumption')
    e_1   $e_{1,t}$ (long_name='Energy consumption of spenders') 
    e_2   $e_{2,t}$ (long_name='Energy consumption of savers') 
    l     $l_t$ (long_name='Labour')
    l_1   $l_t^{1}$ (long_name='Labour of spenders')
    l_2   $l_t^{2}$ (long_name='Labour of savers')
    l_Y   $l_t^{Y}$ (long_name='Labor demand in production sector')   
    l_E   $l_t^{E}$ (long_name='Labor demand in energy production sector')
    a     $k_t$ (long_name='Capital')
    k_agg $k_t^{agg}$ (long_name='Aggregate Capital')
    k_Y   $k_t^{Y}$ (long_name='Capital demand in production sector')
    k_E   $k_t^{E}$ (long_name='Capital demand in energy sector')
    e     $e_t$ (long_name='Energy production')
    e_Y   $e_t^{Y}$ (long_name='Energy input in non-energy sector')
    i     $i_t$ (long_name='Investment')
    i_agg $i_t^{agg}$ (long_name='Aggregate Investment')
    d     $d_t$ (long_name='Dividends')
    g     $g_t$ (long_name='Government Spending')
    t     $t_t$ (long_name='Transfers')
    t_1   $t_t^{1}$ (long_name='Transfers of spenders')
    t_2   $t_t^{2}$ (long_name='Transfers of savers')
    w     $w_t$ (long_name='Wages')
    rK    $r_t^{K}$ (long_name='Capital returns')
    y     $y_t$ (long_name='Output') 
    Y_priv  $Y_t^{exp,priv}$ (long_name='Aggregate expenses - private sector')
    Y_inc_priv $Y_t^{inc,priv}$ (long_name='Aggregate Income - private sector')
    Y  $Y_t^{exp,pu}$ (long_name='Aggregate expenses')
    Y_inc  $Y_t^{inc,pu}$ (long_name='Aggregate Income')
    rF    $r_t^{F}$ (long_name='Risk free rate') 
    lambda_1 $\lambda_t^{1}$ (long_name='The lagrangian with respect to consumption and savings (bonds)') 
    lambda_2 $\lambda_t^{2}$ (long_name='The lagrangian with respect to consumption and savings (bonds)')
    varrho_Y ${varrho}_t^Y$ (long_name='Marginal cost (input shadow price)')
    mc_Y  $mc_t^{Y}$ (long_name='Total marginal cot of production')
    varrho_E ${varrho}_t^E$ (long_name='Marginal cost Energy (input shadow price)')
    p_e   $p_e$ (long_name='Relative energy price')
    p_1   $p_{c,1}$ (long_name='Price index of Spenders consumption bundle')
    p_2   $p_{c,2}$ (long_name='Price index of Savers consumption bundle')
    e_b   $e_t^{b}$ (long_name='The shock variable within the law of motion to consumption (preference shock)')
    e_a   $e_t^{a}$ (long_name='The shock variable within the law of motion to output')
    e_e   $e_t^{e}$ (long_name='The shock variable within the law of motion to energy')
    e_t   $e_t$ (long_name='The shock variable within the law of motion to government transfers')
    e_tau $e_{\tau}$ (long_name='The shock variable within the carbon price')
    inc_1 $Y_{1,t}$ (long_name='Total income of spenders')
    inc_2 $Y_{2,t}$ (long_name='Total income of savers')
    tau_wedge 
    emission_intensity
    share_energy_1
    share_energy_2
    @# if abatment_rigidity == 1
        muE_star
        muY_star
        tau_tilda
        rho_h
    @# endif
    @# if bonds == 1
        rB    $r_t^{F}$ (long_name='Risk free rate') 
        b_h   $b_t^{h}$ (long_name='Bonds demand')
        b_g   $b_t^{g}$ (long_name='Bonds supply')
    @# endif
    res_1
    res_2
    @#if income_incidence == 0
        r_1 
        r_2 
    @# endif
;  

%% Reporting variables
var c_pct c1_pct c2_pct y_pct Y_pct Y_priv_pct Y_inc_priv_pct w_pct co2_pct inc1_pct inc2_pct p1_pct p2_pct pe_pct e_pct rK_pct l_Y_pct k_Y_pct e_Y_pct l_E_pct k_E_pct k_agg_pct res_1_pct res_2_pct e1_pct e2_pct x1_pct x2_pct ec_pct x_pct;


predetermined_variables a;

% 2. Exogenous variables

varexo 
    eta_b  $eta_t^{b}$ (long_name='The exogenous consumption shock')
    eta_e  $eta_t^{e}$ (long_name='The exogenous consumption shock energy')
    eta_a  $eta_t^{a}$ (long_name='The exogenous output shock')
    eta_t  $eta_t^{t}$ (long_name='The exogenous transfers spending shock')
    eta_tau  $eta_t^{c}$ (long_name='The exogenous carbon tax shock')
;                          

% 3. Parameters

parameters 
    z       $TFP_t$ (long_name='TFP')
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
    phi     $\phi$ (long_name='Frisch elasticity of labour disutility')
    l0      $l_0$ (long_name='Hours worked')
    t0      $t_0$ (long_name='Transfers as percentage of carbon revenues')
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
    a_1x    $a_{1,x}$ (long_name='Spenders distribution parameter for goods')
    a_1e    $a_{1,e}$ (long_name='Spenders distribution parameter for energy')
    a_2x    $a_{2,x}$ (long_name='Savers distribution parameter for goods')
    a_2e    $a_{2,e}$ (long_name='Savers distribution parameter for energy')
    eps_x   $\eps_x$ (long_name='Elasticity of substitution energy/goods')
    eps_p   $\eps_y$ (long_name='Elasticity of substitution between the intermediate goods')
    sdev_tau $sdev_a$ (long_name='Standard deviation of the exogenous shock to the carbon tax')
    rho_tau $\rho_a$ (long_name='The shock persistence level for the carbon shock')
    gamma   $\gamma$ (long_name='Climate damage parameter')
    mu     $\mu$ (long_name='Share of carbon revenues to Spenders')
    chi_1    $\chi$ (long_name='Level of disutility')
    chi_2    $\chi$ (long_name='Level of disutility')
    lambda_mu
    b0      $b_0$ (long_name='Fixed supply of bonds')

;

% 3. Calibration

%% Macro parameters

z           = 3.7;      % TFP parameter
beta        = 0.99;     % Time preference
delta       = .025;     % Depreciation capital rate (2.5% quarterly) - Smets and Wouters (2003)
alpha       = .3;       % Share of capital income - GHKT (2014)
alpha_E     = .6;       % Share of capital income in the energy sector (Barrage, 2019)
nu          = .03;      % Share of energy income (Barrage, 2019)
sigma       = 1;      % Inverse of EIS (Barrage, 2019)
chi_1       = 11.6;     % Spenders - level of disutility (calibrated for l_1=1/3 at SS)
chi_2       = 7.2;      % Savers - level of disutility (calibrated for l_1=1/3 at SS)
% chi         = 4.4;      % Level of disutility (we calibrate it such that we have l=1/3 at the Steady state)
phi         = 1.28;     % inverse of Frisch elasticity of labour disutility - 
omega_s     = .2;       % Fraction of spenders
tau_i       = 0;        % Income tax rate
l0          = 1/3;      % Hours worked (1/3 of 24h a day)
b0          = -.92;     % Fixed supply of bonds
@# if energy_share_heterogeneity == 1
    a_1e        = .1618;     % Spenders' distribution parameter for energy
@# else
    a_1e        = .1031;
@# endif
a_1x        = 1 - a_1e; % Spenders' distribution parameter for goods
a_2e        = .1031;     % Savers' distribution parameter for energy
a_2x        = 1 - a_2e; % Savers' distribution parameter for energy
eps_x       = 0.2;     % Elasticity of substitution energy/goods - Känzig (2023)
eps_p       = 5;        % Elasticity  of substitution between goods -  Christopoulou and Vermeulen (2012)
% tau_bar     = 112*10^(-3);    % Steady-state carbon tax
tau_bar     = 0.0587;
phi_1       = 1 - 0.9994; % Emissions decay parameter - Känzig (2023)
phi_0       = 0.5359;   % Emissions staying in atmosphere - Känzig (2023)
@# if damage_function == 1
    gamma   = 5.3*10^(-5); % Climate damage parameter - Känzig (2023)
@# else
    gamma   = 0;
@# endif
omega2_e = .095;
omega1_e = .06;

% Redistribution policy parameters
mu         = 1;         % Share of carbon revenues to spenders
t0         = 0;         % Transfers as % of carbon revenues

%% Environmental parameters
theta_1 = .1;           % Abatement cost parameter
theta_2 = 2.7;          % Abatement cost parameter
phi_Y = .084;           % Emission intensity - Production sector
phi_E = 0.95;           % Emission intensity - Energy sector

%% Shocks paameters
sdev_a     = .01;      % Standard deviation of the exogenous shock to income
rho_a      = .9;       % The shock persistence level
sig_a      = 1;        % The sign of the shock (1 if positive, -1 if negative)
sdev_e     = .01;      % Standard deviation of the exogenous shock to income
rho_e      = .9;       % The shock persistence level
sig_e      = 1;        % The sign of the shock (1 if positive, -1 if negative)
sdev_b     = .01;      % Standard deviation of the exogenous shock to preferences
rho_b      = .9;       % The shock persistence level
sig_b      = 1;        % The sign of the shock (1 if positive, -1 if negative)
sdev_t     = .8;       % Standard deviation of the exogenous shock to preferences
rho_t      = .9;       % The shock persistence level
sig_t      = 1;        % The sign of the shock (1 if positive, -1 if negative)
sdev_tau   = .1;      % Standard deviation of the exogenous shock to the carbon tax
rho_tau    = 0.6;     % Persistence of the carbon tax shock
lambda_mu  = 0.8;


% 4. Model

model;  

%%%% Hand-to-Mouth households %%%% 

[name='Spenders: Marginal utility of consumption']
lambda_1  = c_1^(-sigma)/p_1;

[name='Spenders: Marginal disutility of labour']
@#if inelastic_labor == 1
    l_1 = l0;
@# else 
    (1-tau_i)*w = chi_1*(l_1^phi)/lambda_1;  
@#endif

[name='Spenders: Budget constraint']
@#if income_incidence == 0
    p_1*c_1 = (1-tau_i)*w*l_1 + t_1 + r_1;
@# else
    p_1*c_1 = (1-tau_i)*w*l_1 + t_1 ; 
@# endif

[name='Spenders: Demand function for goods']
x_1 = a_1x * (1/p_1)^(-eps_x) * c_1;

[name='Spenders: Demand function for energy']
e_1 = a_1e * (p_e/p_1)^(-eps_x) * c_1;

[name='Spenders: Price of consumption bundle']
p_1 = (a_1x + a_1e * (p_e)^(1-eps_x))^(1/(1-eps_x));

[name='Spenders: Total income']
inc_1 = w*l_1 + t_1;

[name='Spenders: Energy share']
share_energy_1 = (p_e*e_1)/(p_1*c_1);

%%%% Savers %%%%  

@# if bonds == 1
    [name='Savers: Bonds Euler Equation']
    lambda_2 = beta*rB*lambda_2(+1);
@# endif

[name='Savers: Euler Equation']
lambda_2 = beta*(1+rF(+1))*lambda_2(+1);

[name='Savers: Marginal utility of consumption']
lambda_2 = c_2^(-sigma)/p_2;

[name='Savers: Marginal disutility of labour']
@#if inelastic_labor == 1
    l_2 = l0;
@#else
    (1-tau_i)*w = chi_2*(l_2^phi)/lambda_2;
@#endif

[name='Capital Law of Motion']
a(+1) = (1-delta)*a + i;

[name='Savers: Consumer Budget Constraint']
@# if bonds == 1
    @# if income_incidence == 0
        b_h + p_2*c_2 + i = rB*b_h(-1) + rK*a + (1-tau_i)*w*l_2 + d/(1-omega_s) + t_2 + r_2; 
    @# else
        b_h + p_2*c_2 + i = rB*b_h(-1) + rK*a + (1-tau_i)*w*l_2 + d/(1-omega_s) + t_2; 
    @# endif
@# else
    @# if income_incidence == 0
        p_2*c_2 + i = rK*a + (1-tau_i)*w*l_2 + d/(1-omega_s) + t_2 + r_2;
    @# else
        // p_2*c_2 + i = rK*a + (1-tau_i)*w*l_2 + d/(1-omega_s) + t_2;
        y = x + i_agg + g + theta_1*mu_Y^theta_2*y + theta_1*mu_E^theta_2*e; % Consumption good resource constraint
    @# endif
@# endif

[name='Risk-free rate']
rF = rK - delta;

[name='Savers: Demand function for goods']
x_2 = a_2x * (1/p_2)^(-eps_x) * c_2;

[name='Savers: Demand function for energy']
e_2 = a_2e * (p_e/p_2)^(-eps_x) * c_2;

[name='Savers: Price of consumption bundle']
p_2 = (a_2x + a_2e * (p_e)^(1-eps_x))^(1/(1-eps_x));

[name='Savers: Total income']
inc_2 = rK*a + w*l_2 + d/(1-omega_s) + t_2;

[name='Savers: Energy share']
share_energy_2 = (p_e*e_2)/(p_2*c_2);

@#if income_incidence == 0  
    [name='Lump sum redistribution 1']
    r_1 = (1-omega_s)*(res_2 - res_1);
    [name='Lump sum redistribution 2']
    r_2 = - (omega_s/(1-omega_s)) * r_1;

@#endif
[name='Disposable resources 1']
res_1 = (1 - tau_i)*w*l_1 + t_1; 
[name='Disposable resources 2']
res_2 = rK*a + (1 - tau_i)*w*l_2 + d/(1-omega_s) + t_2 - i; 

%%%% Non-energy producers %%%% 
[name='Production function (Cobb-Douglas)']
y = e_a*exp(-gamma*s)*z*k_Y^(alpha)*e_Y^(nu)*l_Y^(1-alpha-nu);

[name='FOC w.r.t l']
w = varrho_Y*(1-alpha-nu)*y/l_Y;

[name='FOC w.r.t k(-1)']
rK = varrho_Y*alpha*y/k_Y;

[name='FOC w.r.t e']
p_e = varrho_Y*nu*y/e_Y;

[name='Marginal cost of production (FOC w.r.t y)']
varrho_Y = mc_Y - (1-mu_Y)*phi_Y*tau - theta_1*mu_Y^theta_2; 

[name='Total marginal cost of production']
@# if dividends == 1
    mc_Y = (eps_p-1)/eps_p;
@# else
    mc_Y = 1;
@# endif

[name='FOC w.r.t mu']
@# if abatment_rigidity == 1
    rho_h = beta*(1 - lambda_mu);
    tau_tilda = (1 - rho_h)*tau_tilda + rho_h*tau(+1);
    muY_star = ((phi_Y/(theta_1*theta_2))*tau_tilda)^(1/(theta_2-1));
@# else
    tau = theta_1*theta_2*mu_Y^(theta_2-1)/phi_Y;
@# endif

[name='CO2 emissions of non-energy producers']
co2_Y = (1-mu_Y)*phi_Y*y;

[name='Dividends']
d = y*(1-mc_Y);

%%%% Energy producers %%%%

[name='Energy production function']
e = e_e*k_E^(alpha_E)*l_E^(1-alpha_E);

[name='Energy prod.: FOC w.r.t l']
w = varrho_E*(1-alpha_E)*e/l_E;

[name='Energy prod.: FOC w.r.t k(-1)']
rK = varrho_E*alpha_E*e/k_E;

[name='Marginal cost of energy production (FOC w.r.t e)']
varrho_E = p_e - (1-mu_E)*phi_E*tau - theta_1*mu_E^theta_2;

[name='FOC w.r.t mu_E']
@# if abatment_rigidity == 1
    muE_star = ((phi_E/(theta_1*theta_2))*tau_tilda)^(1/(theta_2-1));
@# else
    tau = theta_1*theta_2*mu_E^(theta_2-1)/phi_E;
@# endif

[name='CO2 emissions of energy producers']
co2_E = (1-mu_E)*phi_E*e;

[name='Emission intensity (level)']
emission_intensity = co2_agg / Y;

@#if abatment_rigidity == 1
    [name='Abatement smoothing (non-energy)']
    mu_Y = lambda_mu*muY_star + (1 - lambda_mu)*mu_Y(-1);

    [name='Abatement smoothing (energy)']
    mu_E = lambda_mu*muE_star + (1 - lambda_mu)*mu_E(-1);
@# else 
@#endif

%%%% Government Spending %%%% 

[name='Transfers as % of Carbon revenues']
t = t0*tau*co2_agg*e_t; 

[name='Gov Budget Constraint']
@# if bonds == 1
    rB*b_g(-1) + g + t = b_g + tau_i*w*l + tau*co2_agg;
@# else
    g + t = tau_i*w*l + tau*co2_agg;
@#endif

[name='Carbon tax rule']
tau = tau_bar*e_tau;

[name='Total Transfers']
t = omega_s*t_1 + (1-omega_s)*t_2;

[name='Transfers to Spenders']
t_1 = (mu/omega_s)*t;

[name='Share of carbon revenues']
tau_wedge = (co2_agg*tau)/Y;

@# if bonds == 1
    [name='Supply of bonds']
    b_g = b0*y ;

    [name='No deficit Spending']
    (1-omega_s)*b_h + b_g = 0;
@# endif

%%%% Aggregation %%%% 

[name='Aggregate consumption of non-energy']
x = omega_s*x_1 + (1-omega_s)*x_2;

[name='Aggregate consumption of energy']
e_c = omega_s*e_1 + (1-omega_s)*e_2;

[name='Aggregate Consumption']
c = omega_s*c_1 + (1-omega_s)*c_2;

[name='Aggregate Labour']
l = omega_s*l_1 + (1-omega_s)*l_2;

[name='Aggregate Capital']
k_agg = (1-omega_s)*a ; 

[name='Aggregate Investment']
i_agg = k_agg(+1) - (1-delta)*k_agg; 

[name='Aggregate emissions']
co2_agg = co2_Y + co2_E;

[name='Carbon concentration']
s = (1-phi_1)*s(-1) + phi_0*co2_agg;

[name='Aggregate income - private']
% Y_inc_priv = (rK*k_agg + (1-tau_i)*w*l + d + t);
Y_inc_priv = omega_s*inc_1 + (1-omega_s)*inc_2;

[name='Aggregate income']
Y_inc  = (rK*k_agg + w*l + d + tau*co2_agg);

[name='Aggregate expenses - private']
Y_priv = x + i_agg + p_e*e_c;

[name='Aggregate expenses']
Y = x + p_e*e_c + i_agg + g;

%%%% Market Clearing conditions %%%%

[name='Capital market clearing']
k_Y + k_E = k_agg;

[name='Labor market clearing']
l_Y + l_E = l;

[name='Energy market clearing']
e = e_c + e_Y;


%%%% Exogenous shock law of motion %%%% 

[name='Law of motion of the carbon price']
log(e_tau) = rho_tau*log(e_tau(-1)) +eta_tau;

[name='Law of motion of the output shock']
log(e_a) = rho_a*log(e_a(-1)) +sig_a*eta_a;

[name='Law of motion of the output shock energy']
log(e_e) = rho_e*log(e_e(-1)) +sig_e*eta_e;

[name='Law of motion of the preference shock']
log(e_b) = rho_b*log(e_b(-1)) +sig_b*eta_b;

[name='Law of motion of the transfer shock']
log(e_t) = rho_t*log(e_t(-1)) +sig_t*eta_t;

% Reporting variables - Percent deviations (100 * (x/SS - 1))
c_pct    = 100*(c   /steady_state(c)    - 1);
c1_pct   = 100*(c_1 /steady_state(c_1)  - 1);
c2_pct   = 100*(c_2 /steady_state(c_2)  - 1);
y_pct    = 100*(y   /steady_state(y)    - 1);
w_pct    = 100*(w   /steady_state(w)    - 1);
co2_pct  = 100*(co2_agg/steady_state(co2_agg) - 1);
p1_pct   = 100*(p_1 /steady_state(p_1)  - 1);
p2_pct   = 100*(p_2 /steady_state(p_2)  - 1);
pe_pct   = 100*(p_e /steady_state(p_e)  - 1);
inc1_pct = 100*(inc_1/steady_state(inc_1) - 1);
inc2_pct = 100*(inc_2/steady_state(inc_2) - 1);
e_pct    = 100*(e/steady_state(e) - 1);
rK_pct   = 100*(rK/steady_state(rK) - 1);
Y_pct    = 100*(Y/steady_state(Y) - 1);
Y_priv_pct = 100*(Y_priv/steady_state(Y_priv) - 1);
Y_inc_priv_pct = 100*(Y_inc_priv/steady_state(Y_inc_priv) - 1);
k_Y_pct = 100*(k_Y/steady_state(k_Y) - 1);
l_Y_pct = 100*(l_Y/steady_state(l_Y) - 1);
e_Y_pct = 100*(e_Y/steady_state(e_Y) - 1);
k_E_pct = 100*(k_E/steady_state(k_E) - 1);
l_E_pct = 100*(l_E/steady_state(l_E) - 1);
k_agg_pct = 100*(k_agg/steady_state(k_agg) - 1);
res_1_pct = 100*(res_1/steady_state(res_1) - 1);
res_2_pct = 100*(res_2/steady_state(res_2) - 1);
e1_pct   = 100*(e_1 /steady_state(e_1)  - 1);
e2_pct   = 100*(e_2 /steady_state(e_2)  - 1);
x1_pct   = 100*(x_1 /steady_state(x_1)  - 1);
x2_pct   = 100*(x_2 /steady_state(x_2)  - 1);
ec_pct   = 100*(e_c /steady_state(e_c)  - 1);
x_pct   = 100*(x /steady_state(x)  - 1);


end;

% 5. Steady state

steady_state_model;
e_tau = 1;
e_a = 1 ;
e_e = 1 ;
e_b = 1 ;
e_t = 1 ;
@#if dividends == 0
    mc_Y = 1;
@#elseif dividends == 1
    mc_Y = (eps_p-1)/eps_p;
@#endif
tau = tau_bar;
@#if abatment_rigidity == 0
    mu_E = (tau*phi_E/(theta_2*theta_1))^(1/(theta_2-1));
    mu_Y = (tau*phi_Y/(theta_2*theta_1))^(1/(theta_2-1));
@#elseif abatment_rigidity == 1
    tau_tilda =tau;
    rho_h = beta*(1 - lambda_mu);
    muE_star = (tau*phi_E/(theta_2*theta_1))^(1/(theta_2-1));
    muY_star = (tau*phi_Y/(theta_2*theta_1))^(1/(theta_2-1));
    mu_E = muE_star;
    mu_Y = muY_star;
@#endif
@#if bonds == 1
    rB = 1/beta;
@# endif

varrho_Y = mc_Y - theta_1*mu_Y^theta_2 - phi_Y*tau*(1-mu_Y);
rF = 1/beta - 1;
rK = rF + delta;

%% Solver
@# if bonds == 1
    [lambda_2,c_2,w,l_2,i,t_2,y,a,g,lambda_1,c_1,l_1,t_1,t,l,k_agg,i_agg,co2_Y, co2_E, varrho_E, k_Y, k_E, l_Y, l_E, e, co2_agg, e_Y, d, s, p_1, e_1,p_2, e_2,p_e,b_g,b_h] = get_steady_state_model_thesis(beta,sigma,chi_1,chi_2,phi,delta,alpha,alpha_E,nu,tau_i,omega_s,theta_1,theta_2,phi_Y,phi_E,tau_bar, eps_p, phi_1, phi_0, gamma, mu, eps_x,a_1x,a_1e,a_2x,a_2e,l0,t0,z,mu_E,mu_Y,mc_Y,tau,varrho_Y,rK, rB,b0);
@# elseif income_incidence == 0
    [lambda_2,c_2,w,l_2,i,t_2,y,a,g,lambda_1,c_1,l_1,t_1,t,l,k_agg,i_agg,co2_Y, co2_E, varrho_E, k_Y, k_E, l_Y, l_E, e, co2_agg, e_Y, d, s, p_1, e_1,p_2, e_2,p_e,r_1,r_2,res_1,res_2] = get_steady_state_model_thesis_income_incidence(beta,sigma,chi_1,chi_2,phi,delta,alpha,alpha_E,nu,tau_i,omega_s,theta_1,theta_2,phi_Y,phi_E,tau_bar, eps_p, phi_1, phi_0, gamma, mu, eps_x,a_1x,a_1e,a_2x,a_2e,l0,t0,z,mu_E,mu_Y,mc_Y,tau,varrho_Y,rK);
@# else
    [lambda_2,c_2,w,l_2,i,t_2,y,a,g,lambda_1,c_1,l_1,t_1,t,l,k_agg,i_agg,co2_Y, co2_E, varrho_E, k_Y, k_E, l_Y, l_E, e, co2_agg, e_Y, d, s, p_1, e_1,p_2, e_2,p_e] = get_steady_state_model_thesis(beta,sigma,chi_1,chi_2,phi,delta,alpha,alpha_E,nu,tau_i,omega_s,theta_1,theta_2,phi_Y,phi_E,tau_bar, eps_p, phi_1, phi_0, gamma, mu, eps_x,a_1x,a_1e,a_2x,a_2e,l0,t0,z,mu_E,mu_Y,mc_Y,tau,varrho_Y,rK);
@# endif

x_1 = a_1x * (1/p_1)^(-eps_x) * c_1;
x_2 = a_2x * (1/p_2)^(-eps_x) * c_2;

inc_1 = w*l_1 + t_1;
inc_2 = rK*a + w*l_2 + d/(1-omega_s) + t_2;
@# if income_incidence == 1
    res_1 = (1 - tau_i)*w*l_1 + t_1; 
    res_2 = rK*a + (1 - tau_i)*w*l_2 + d/(1-omega_s) + t_2 - i; 
@# endif

c = omega_s*c_1 + (1-omega_s)*c_2;
x = omega_s*x_1 + (1-omega_s)*x_2;
e_c = omega_s*e_1 + (1-omega_s)*e_2;
Y = x + p_e*e_c + i_agg + g;
Y_inc = rK*k_agg + w*l + d + tau*co2_agg;
Y_priv = x + i_agg + p_e*e_c;
Y_inc_priv = rK*k_agg + (1-tau_i)*w*l + d + t;

tau_wedge = (co2_agg*tau)/Y;
emission_intensity = co2_agg / Y;
share_energy_1 = (p_e*e_1)/(p_1*c_1);
share_energy_2 = (p_e*e_2)/(p_2*c_2);

c_pct = 0;  c1_pct = 0;  c2_pct = 0;  y_pct = 0;  w_pct = 0; co2_pct = 0; 
p1_pct = 0; p2_pct = 0; pe_pct = 0; inc1_pct = 0; inc2_pct = 0; e_pct = 0; 
rK_pct = 0; Y_pct= 0; Y_priv_pct=0; Y_inc_priv_pct = 0; l_Y_pct = 0; k_Y_pct = 0;
e_Y_pct = 0;  l_E_pct=0; k_E_pct=0; k_agg_pct=0; res_1_pct = 0; res_2_pct = 0; x1_pct = 0; x2_pct = 0; e1_pct = 0; e2_pct = 0;
ec_pct = 0;x_pct = 0;

end;

% 6. Compute steady state, check resisuals and BK conditions

steady;     % Provides you with the steady states values
resid;      % Allow you to check if you computed correctly the steady states equations
check;      % Checks the Blanchard Khan conditions to see if the model is well set

% 7. Report variales to calibrate the model

%% Frish labor supply elasticity
l0 = 1/3;
w_ss  = oo_.steady_state(strmatch('w',        M_.endo_names, 'exact'));
lam1  = oo_.steady_state(strmatch('lambda_1', M_.endo_names, 'exact'));
lam2  = oo_.steady_state(strmatch('lambda_2', M_.endo_names, 'exact'));
chi1_hat = (1 - tau_i) * w_ss * lam1 / (l0^phi);
chi2_hat = (1 - tau_i) * w_ss * lam2 / (l0^phi);
disp(['implied chi (spenders) = ' num2str(chi1_hat)]);
disp(['implied chi (savers)   = ' num2str(chi2_hat)]);

%% Consumption distribution parameters 
s0_1 = 0.1158;
s0_2 = 0.0724;
pe_ss  = oo_.steady_state(strmatch('p_e',        M_.endo_names, 'exact'));
p_1_ss  = oo_.steady_state(strmatch('p_1', M_.endo_names, 'exact'));
p_2_ss  = oo_.steady_state(strmatch('p_2', M_.endo_names, 'exact'));
a1_e_hat = s0_1*(pe_ss/p_1_ss)^(eps_x-1);
a2_e_hat = s0_2*(pe_ss/p_2_ss)^(eps_x-1);
disp(['implied a_e (spenders) = ' num2str(a1_e_hat)]);
disp(['implied a_e (savers)   = ' num2str(a2_e_hat)]);

% 7. Declare shock variable

shocks;

var eta_tau; stderr sdev_tau ; 

end;

% 8. Perform the simulations
stoch_simul(order=1, pruning, irf=13, tex)
    y k_Y l_Y e_Y k_agg l_E k_E i_agg t t_1 t_2 d c_pct c1_pct c2_pct y_pct Y_pct Y_priv_pct Y_inc_priv_pct w_pct co2_pct p1_pct p2_pct pe_pct l_1 l_2 inc1_pct inc2_pct e_pct rK_pct co2_agg s tau mu_Y mu_E varrho_Y varrho_E k_Y_pct l_Y_pct e_Y_pct k_E_pct l_E_pct k_agg_pct a res_1_pct res_2_pct e1_pct e2_pct ec_pct e_c g x x_1 x_2 x1_pct x2_pct x_pct;
