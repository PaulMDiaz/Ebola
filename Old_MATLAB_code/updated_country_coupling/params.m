function params()
% parameter initialization
alpha  = 0.024/365; % population growth constant (known empircally) 
beta1  = IG(1);
beta2  = IG(2);
beta3  = IG(3);
delta  = 1/21;%1/incubation period for Ebola Virus
gamma1 = 1/10;
gamma2 = IG(4);
psi    = IG(5);
rho1   = 0.7;
rho2   = 0.6;
omega  = 1;  %population death constant 
