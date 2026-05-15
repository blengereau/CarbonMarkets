
%% Two Agent RBC/NK model 
%% By Ghassane Benmir, January 2025

close all;close all;
warning off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Zero declare all the options:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
@#define model_type = 0
            // 0: RBC
            // 1: NK

@#define policy = 1
            // 0: Gov Spending adjust
            // 1: Transfers adjust


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First declare all the endogenous variables: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

var 
    c       % Consumption  
    c_1     % Consumption of spenders
    c_2     % Consumption  of savers
    l       % Labour
    l_1     % Labour of spenders
    l_2     % Labour of savers
    k       % Capital
    k_agg   % Aggregate Capital
    i       % Investment
    i_agg   % Aggregate Investment
    g       % Government Spending
    t       % Transfers
    t_1     % Transfers of spenders
    t_2     % Transfers of savers
    w       % Wages
    rK      % Capital returns
    y       % Output 
    r_real  % Real interest rate
    b_h     % Bonds demand
    b_g     % Bonds supply
    q       % The lagrangian with respect to capital and investment  - also represents the cost of capital (the Tobin Q)
    lambda_1% The lagrangian with respect to consumption and savings (bonds) - also represents the marginal utility of consumption
    lambda_2% The lagrangian with respect to consumption and savings (bonds) - also represents the marginal utility of consumption
    mc      % Marginal cost (input shadow price)
    e_b     % The shock variable within the law of motion to consumption (preference shock)
    e_a     % The shock variable within the law of motion to output
    e_gt     % The shock variable within the law of motion to government transfers or public spending
    @#if model_type == 1
    pi      % Inflation rate
    rF      % Nominal Risk free rate 
    e_r
    @#endif
    c_log 
    c1_log 
    c2_log
    w_log
;  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Second declare all the exogenous variables: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

varexo 
    eta_b     % The exogenous consumption shock
    eta_a     % The exogenous output shock
    eta_gt     % The exogenous tranfers spending shock
    @#if model_type == 1
    eta_r     % The exogenous monetary shock
    @#endif
;                          

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Third declare all the parameters: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

parameters 
    beta
    alpha
    delta
    chi_1
    chi_2
    phi
    l0 
    b0 
    t0
    g0
    sdev_a 
    sig_a 
    rho_a 
    sdev_b
    sig_b
    rho_b 
    sdev_gt
    rho_gt
    sig_gt
    omega_s
    tau_i
    trans_adjust
    transfer_fix
    @#if model_type == 1
    theta         
    kappa         
    phi_pi        
    phi_y        
    sdev_r        
    rho_r         
    sig_r        
    @#endif
;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Third (bis) calibration for all the params 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Macro Params
beta       = 0.9959;   % Time preference
delta      = .025;    % Depreciation capital rate (2.5% quarterly) 
alpha      = 1/3;      % Share of capital income
chi_1        = 8;        % Level of disutility (we calibrate it such that we have l_1=1/3 at the Steady state)
chi_2        = 8;        % Level of disutility (we calibrate it such that we have l_2=1/3 at the Steady state)
phi        = 1;        % Frisch elasticity of labour
omega_s    = 0.9;//.3;     % % of spenders (Hand to mouth) in the EU (30%)
tau_i      = 0.25;     % 25% of income tax  
l0         = 1/3;      % Hours worked (1/3 of 24h a day)
b0         = -.92;     % Fixed supply of bonds 92% of Dept to GDP ratio (annual-i.e. 0.23 quarterly b_g=.23*4*y => b_g = .92*y)
t0         = 0.07;     % Gog spending as % of GDP (10% annual here)
g0         = .2;       % Public spending to gdp ration (20%)
trans_adjust = 0;      % Transfer adjustment
transfer_fix = 0;      % Permant increase in transfers (make sure it is 0 for the case of the shock!!)
    
@#if model_type == 1
%% New Keynesian Parameters
theta      = 9;         % Elasticty of substitution between intermediate outputs
kappa      = 45;        % Adjusment cost coefficient for the Phillips Curve 
phi_pi     = 1.25;       % Taylor rule response to inflation
phi_y      = 0;//0.3;       % Taylor rule response to output gap
@#endif

%% Shocks Params
sdev_a     = .01;      % Standard deviation of the exogenous shock to income
rho_a      = .9;       % The shock persistence level
sig_a      = 1;        % The signe of the shock (1 if positive, -1 is negative)
sdev_b     = .01;      % Standard deviation of the exogenous shock to preferences
rho_b      = .9;       % The shock persistence level
sig_b      = 1;        % The signe of the shock (1 if positive, -1 is negative)
sdev_gt     = .06;       % Standard deviation of the exogenous shock to preferences
rho_gt      = .9;       % The shock persistence level
sig_gt      = 1;        % The signe of the shock (1 if positive, -1 is negative)
sdev_r     = .01;      % Standard deviation of the exogenous shock
rho_r      = .5;       % The shock persistence level
sig_r      = 1;        % The signe of the shock (1 if positive, -1 is negative)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fourth delcare your model (i.e. the FOCs and Equations of the model) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

model;  

%%%% Households Side (Spenders) %%%%

[name='Marginal utility of consumption']
lambda_1  = 1/c_1;

[name='Labor market equilibrium']
(1-tau_i)*w = chi_1*(l_1^phi)/lambda_1;

[name='Consumer Budget Constraint']
@#if model_type == 0
c_1 = (1-tau_i)*w*l_1 + t_1 ; 
@#elseif model_type == 1
c_1 = (1-tau_i)*w*l_1 + t_1 ;//+ y*(1-mc -kappa/2*(pi-steady_state(pi))^2); 
@#endif

%%%% Households Side (Savers) %%%%

[name='FOC w.r.t to b: Euler Equation']
lambda_2 = e_b(+1)*beta*r_real*lambda_2(+1)/e_b;

[name='FOC w.r.t to c: Marginal Utility of Consumption']
c_2 = 1/lambda_2;

[name='FOC w.r.t to l: Disutility of labour']
(1-tau_i)*w = chi_2*l_2^phi/lambda_2;

[name='FOC w.r.t k']
q = e_b(+1)*beta*lambda_2(+1)*rK(+1) + e_b(+1)*beta*q(+1)*(1-delta);

[name='FOC w.r.t i']
lambda_2 = q ;

[name='Capital Law of Motion']
k = (1-delta)*k(-1) + i;

[name='Consumer Budget Constraint']
@#if model_type == 0
//b_h = r_real*b_h(-1) + rK*k(-1) + (1-tau_i)*w*l_2 + t_2 - c_2 - i   ; 
y = c + i_agg + g ; //(RBC) This is the same equation (always good to try to make sure all is good) 
@#elseif model_type == 1
b_h = r_real*b_h(-1) + rK*k(-1) + (1-tau_i)*w*l_2 + t_2 - c_2 - i + y*(1-mc -kappa/2*(pi-steady_state(pi))^2); 
//(1- omega_s*(1 -mc))*y  = c + i_agg + g + (1-omega_s)*(kappa/2*(pi-steady_state(pi))^2)*y ; //(NK) This is the same equation (always good to try to make sure all is good) 
@#endif
          
%%%% Different Aggregation %%%%

[name='Aggreg of Consumption']
c = omega_s*c_1 + (1-omega_s)*c_2;

[name='Aggreg Labour']
l = omega_s*l_1 + (1-omega_s)*l_2;

[name='Total Transfers']
t = (trans_adjust + omega_s)*t_1 + (1-(omega_s+trans_adjust))*t_2;

[name='Agg Capital']
k_agg = (1-omega_s)*k ; 

[name='Agg Capital Law']
i_agg = (1-omega_s)*i; 
//i_agg = k_agg - (1-delta)*k_agg(-1) ; // Good to double check as well. Here I checked already and all is good.

%%%% Firms Side %%%%

[name='Production function (Cobb-Douglas)']
y = e_a*k_agg(-1)^(alpha)*l^(1-alpha);

[name='FOC w.r.t l']
w = mc*(1-alpha)*y/l;

[name='FOC w.r.t k(-1)']
rK = mc*alpha*y/k_agg(-1);

@#if model_type == 0
[name='FOC w.r.t y (this is the Marginal cost)']
mc = 1;
@#endif

%%%% New Keynsian Part %%%%

@#if model_type == 1
[name='The New Phillips Curve']
(pi-steady_state(pi))*pi=beta*(lambda_2(+1)/lambda_2*y(+1)/y*pi(+1)*(pi(+1)-steady_state(pi)))+theta/kappa*(mc-(theta-1)/theta);

[name='Fisher Equation']
rF = r_real*pi(+1);

[name='Taylor Rule']
rF = rF(-1) + phi_pi*(pi-pi(-1))  + phi_y*(y-y(-1)) + log(e_r);

@#endif

%%%% Government Spending %%%%

@#if policy == 0
[name='public spending as % of GDP']
g = g0*y*e_gt; 
@#elseif policy == 1
[name='Transfers as % of GDP']
t = t0*y*e_gt; % This is a shock to transfers
@#endif

[name='Supply of bonds']
b_g = b0*y ;

[name='Share of Transfers to Sp']
t_1 = t  ;

[name='Gov Budget Constraint']
b_g + g*e_gt = r_real*b_g(-1) + tau_i*w*l - t; 


%%%% Market Clearing conditions (to close the model)

[name='No deficit Spending']
(1-omega_s)*b_h + b_g = 0;

%%%% Exogenous shock law of motion %%%%

[name='Law of motion of the output shock']
log(e_a) = rho_a*log(e_a(-1)) +sig_a*eta_a;

[name='Law of motion of the preference shock']
log(e_b) = rho_b*log(e_b(-1)) +sig_b*eta_b;

[name='Law of motion of the transfer shock']
log(e_gt) = rho_gt*log(e_gt(-1)) +sig_gt*eta_gt;

@#if model_type == 1
[name='Law of motion of the monetary shock']
log(e_r) = rho_r*log(e_r(-1)) +sig_r*eta_r;
@#endif

// Some log definitions to help reading the IRFs as % deviation from the SS.
c_log=log(c) ;
c1_log=log(c_1) ;
c2_log=log(c_2) ;
w_log=log(w);

end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fifth delcare your steady state 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

steady_state_model;

// e_a = 1 ;
// e_b = 1 ;
// e_gt = 1 ;
// e_r = 1 ;
// 
// @#if model_type == 0
// mc = 1;
// @#elseif model_type == 1
// mc = (theta-1)/theta;
// pi = 1;
// @#endif
// 
// l = 1/3; // Then find chi that allows for this! 
// chi = chi_1;
// chi_2 = chi_1;
// %% Using a solver
// [lambda_2,r_real,c_2,w,l_2,q,rK,i,b_h,t_2,y,k,g,lambda_1,c_1,l_1,t_1,c,t,chi,b_g,k_agg,i_agg] = get_steady_state_two_agents(g0,mc,beta,l,phi,delta,alpha,t0,b0,tau_i,omega_s,trans_adjust,transfer_fix);
// 
// @#if model_type == 1
// rF = r_real ;
// @#endif

// c_log=log(c) ;
// c1_log=log(c_1) ;
// c2_log=log(c_2) ;
// w_log=log(w);

e_a = 1 ;
e_b = 1 ;
e_gt = 1 ;
e_r = 1 ;

@#if model_type == 0
mc = 1;
r_real = 1/beta;
@#elseif model_type == 1
mc = (theta-1)/theta;
pi = 1;
r_real = 1/beta;
rF = r_real*pi;
@#endif


l = 1/3;
rK = 1/beta - (1-delta) ;
k_agg_y = mc*alpha/rK; // Here I use the FOC k of firms to get k/y which I call k_y
y = (k_agg_y^(alpha)*l^(1-alpha))^(1/(1-alpha));
k_agg = k_agg_y*y;
i_agg = delta*k_agg ;
k = (1-omega_s)^-1*k_agg;
i = delta*k;
w = mc*(1-alpha)*y/l;
b_g = b0*y ;
@#if policy == 0
g = g0*y; 
t = r_real*b_g + tau_i*w*l - b_g - g; 
@#elseif policy == 1
t = t0*y*e_gt; 
g  = r_real*b_g + tau_i*w*l - b_g - t ; 
@#endif
b_h = -(1-omega_s)^-1*b_g ;
t_1 = t  ;
t_2 = 1/(1-(omega_s+trans_adjust))*t -(trans_adjust + omega_s)/(1-(omega_s+trans_adjust))*t_1 ;
l_1 = 1/3;
c_1 = (1-tau_i)*w*l_1 + t_1 ; 
l_2 = 1/(1-omega_s)*l - omega_s/(1-omega_s)*l_1;
@#if model_type == 0
c_2 = r_real*b_h + rK*k+ (1-tau_i)*w*l_2 + t_2  - i - b_h   ; 
@#elseif model_type == 1
c_2 = r_real*b_h + rK*k + (1-tau_i)*w*l_2 + t_2 - b_h - i + y*(1-mc); 
@#endif
c = omega_s*c_1 + (1-omega_s)*c_2;
lambda_1 = 1/c_1;
lambda_2 = 1/c_2;
chi_1 = (1-tau_i)*w/(l_1^phi/lambda_1) ;
chi_2 = (1-tau_i)*w/(l_2^phi/lambda_2) ;
q = lambda_2 ;
c_log=log(c) ;
c1_log=log(c_1) ;
c2_log=log(c_2) ;
w_log=log(w);

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

//var eta_a; stderr sdev_a ; 
//var eta_b; stderr sdev_b ; 
var eta_gt; stderr sdev_gt ; 
@#if model_type == 1
//var eta_r; stderr sdev_r ; 
@#endif

end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Eight performe the simulation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

@#if model_type == 0
stoch_simul (order=1, pruning, irf=60) g y r_real c w c1_log c2_log c_log w_log;//c_1 c_2 l_1 l_2 c t g y i k l q rK rF w b_h b_g; 
@#elseif model_type == 1
stoch_simul (order=1, pruning, irf=60) g y pi rF r_real c w c1_log c2_log c_log w_log;//c_1 c_2 l_1 l_2 c t g y i k l q rK rF w b_h b_g; 
@#endif


