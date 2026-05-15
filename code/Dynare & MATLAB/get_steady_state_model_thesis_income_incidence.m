function [lambda_2,c_2,w,l_2,i,t_2,y,a,g,...
    lambda_1,c_1,l_1,t_1,t,l,k_agg,i_agg,...
    co2_Y, co2_E, varrho_E,... 
    k_Y, k_E, l_Y, l_E, e, co2_agg, e_Y, d,...
    s, p_1, e_1,p_2, e_2, p_e,r_1,r_2,res_1,res_2] = ...
    get_steady_state_model_thesis_income_incidence( ...
    beta,sigma, chi_1, chi_2, phi,delta,alpha,alpha_E,nu,tau_i,omega_s, ...
    theta_1,theta_2,phi_Y,phi_E,tau_bar, eps_p, phi_1, phi_0, ...
    gamma, mu, eps_x,a_1x,a_1e,a_2x,a_2e,l0,t0,z, ...
    mu_E,mu_Y,mc_Y,tau,varrho_Y,rK)
    % rB, b0...
    % b_g,b_h


    % rB, b0...
func = @(x)([

% HOUSEHOLDS

%% Spenders
x(10)  - x(11)^(-sigma)/x(30);                      % 'Spenders: Marginal utility of consumption'
% chi - x(3)*x(10)*(1/l1^phi);
x(3) - (chi_1*(x(12)^phi)/x(10))/(1-tau_i);         % 'Spenders: Marginal disutility of labor'
% x(12) - l0;   
x(11) - ((1-tau_i)*x(3)*x(12) + x(13) + x(35))/x(30) ;
% x(11) - ((1-tau_i)*x(3)*x(12) + x(13))/x(30) ;     % 'Spenders: Budget constraint'
x(31) - a_1e*(x(34)/x(30))^(-eps_x)*x(11);           % 'Spenders: Demand function for energy'
x(30) - (a_1x + a_1e*(x(34))^(1-eps_x))^(1/(1-eps_x)); % 'Spenders: Price of consumption bundle'

%% Savers
x(1) - x(2)^(-sigma)/x(32);                         % 'Savers: Marginal utility of consumption'
x(3) - (chi_2*(x(4)^phi)/x(1))/(1-tau_i);               % 'Savers: Marginal disutility of labor'
% x(4) - l0;             
x(8) - ((1-delta)*x(8) + x(5));                   % Capital Law of motion
% x(5) - (rK*x(8) + (1-tau_i)*x(3)*l1 + x(28)/(1-omega_s) + x(6) - x(32)*x(2));
x(5) - (rK*x(8) + (1-tau_i)*x(3)*x(4) + x(28)/(1-omega_s) + x(6) + x(36) - x(32)*x(2)); % Savers budget constraint
% x(5) - (rK*x(8) + (1-tau_i)*x(3)*x(4) + x(28)/(1-omega_s) + x(6) - x(32)*x(2) + x(36)*(rB-1)); % Savers budget constraint
x(33) - a_2e*(x(34)/x(32))^(-eps_x)*x(2);           % 'Savers: Demand function for energy'
x(32) - (a_2x + a_2e*(x(34))^(1-eps_x))^(1/(1-eps_x)); % 'Savers: Price of consumption bundle'

% FIRMS

%% Non-energy 
x(7) - exp(-gamma*x(29))*z*x(21)^(alpha)*x(27)^(nu)*x(23)^(1-alpha-nu);   % Y
x(3) - varrho_Y*(1-alpha-nu)*x(7)/x(23);                                   % w
rK - varrho_Y*alpha*x(7)/x(21);                                          % rK
x(34) - varrho_Y*nu*x(7)/x(27);                                            % p_e
x(18) - (1-mu_Y)*phi_Y*x(7);                                           % co2_Y

%% Energy
x(25) - x(22)^alpha_E*x(24)^(1-alpha_E);                                % e
x(22) - x(20)*(alpha_E)*x(25)/rK;                                     % k_E
x(24) - x(20)*(1-alpha_E)*x(25)/x(3);                                   % l_E
x(19) - (1-mu_E)*phi_E*x(25);                                          % co2_E
x(20) - ( x(34) - (1-mu_E)*phi_E*tau - theta_1*mu_E^theta_2);       % varrho_E

%% Markup and profits
x(28) - x(7)*(1-mc_Y);                             % Dividends / Profits

% CLIMATE BLOCK

x(29) - ((1-phi_1)*x(29) + phi_0*x(26));           % Carbon concentration

% GOVERNMENT

%% Redistribution and spending
x(14) - t0*tau*x(26);                               % Transfers as % of Carbon revenues
x(9) - (tau*x(26)- x(14)) ;                         % Government Budget Constraint
% x(9) - (tau*x(26)- x(14) + x(35)*(1-rB)) ; 
x(13) - (mu/omega_s)*x(14);                         % Transfers to spenders
% x(35) - b0*x(7)                                     % bonds supply
% x(35) + x(36)*(1-omega_s)                           % No deficit spending

% AGREGATION

x(14) - (omega_s*x(13) + (1-omega_s)*x(6));         % Aggregate transfers
x(15) - (omega_s*x(12) + (1-omega_s)*x(4));         % Aggregate labor
x(16) - (1-omega_s)*x(8);                          % Aggregate Asset level
x(17) - (x(16) - (1-delta)*x(16)) ;                 % Aggregate Investment
x(26) - x(18) - x(19);                              % Aggregate carbon level

% MARKET CLEARING
x(21) - x(16) + x(22);                              % k_Y + k_E = k_agg;
% x(23) - l + x(24); 
x(23) - x(15) + x(24);                              % l_Y + l_E = l;
x(27) - x(25) + (omega_s*x(31) + (1-omega_s)*x(33));% Energy market clearing

% LUMP SUM REDISTRIBUTION
x(35) - ((1-omega_s)*(x(38)-x(37)));
x(36) + (omega_s/(1-omega_s))*x(35);
x(37) - (x(3)*x(12) + x(13));
x(38) - (rK*x(8) + (1-tau_i)*x(3)*x(4) + x(28)/(1-omega_s) + x(6) - x(5));

    ]);

% INITIAL VALUE
x0 = ones(38,1);

%Set some initial values for the solver
x0(7)  = 5;   % y


options=optimset('Algorithm','trust-region-dogleg','display','off','MaxFunEvals',3000,'MaxIter',3000);
rez=fsolve(func,x0,options);

lambda_2 = rez(1);
c_2 = rez(2);
w= rez(3);
l_2= rez(4);
i= rez(5);
t_2= rez(6);
y= rez(7);
a= rez(8);
g= rez(9);
lambda_1= rez(10);
c_1= rez(11);
l_1= rez(12);
t_1= rez(13);
t= rez(14);
l= rez(15);
k_agg= rez(16);
i_agg = rez(17);
co2_Y = rez(18);
co2_E = rez(19);
varrho_E = rez(20);
k_Y = rez(21);
k_E = rez(22);
l_Y = rez(23);
l_E = rez(24);
e = rez(25);
co2_agg = rez(26);
e_Y = rez(27);
d   = rez(28);
s   = rez(29);
p_1 = rez(30);
e_1 = rez(31);
p_2 = rez(32);
e_2 = rez(33);
p_e = rez(34);
r_1 = rez(35);
r_2 = rez(36);
res_1 = rez(37);
res_2 = rez(38);