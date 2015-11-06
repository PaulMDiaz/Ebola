%SIERRALEONE
clear;clc;
tic

%Manner in which to order variables
%alpha, beta1, beta2, beta3, delta, gamma1, gamma2, psi, rho1, rho2, omega

Pop = 6.092e6;
%Initial guesses for parameters
beta1 = 0.238/Pop;%0.25/Pop;
beta2 = 0.28/Pop;%0.192/Pop;
beta3 = 0.132/Pop;%0.126/Pop;
gamma2 = 0.139;%0.135;
psi = 0.4;%0.4;
omega = 0.263;%0.315;
gamma1 = 0.047;%0.055;
rho1 = 0.847;%0.9;
rho2 = 0.714;%0.85;

initR0 = (1/(gamma1+psi))*(beta1*Pop + (beta3*Pop*psi/gamma2) + (beta2*Pop*rho1*gamma1/omega))

beta1*Pop < 0.08 || beta1*Pop > 0.4 ||beta2*Pop < 0.08 || beta2*Pop > 0.9 || beta3*Pop < 1e-2 ||  beta3*Pop > 0.2...
    || gamma1 < 0.04 || gamma1 > 0.15 || gamma2 < 0.01 || gamma2 > 0.4...
    || rho1 < 0.6 || rho1 > 0.98 || rho2 < 0.7 || rho2 > 0.9...
    || psi < 0.01 || psi > 0.45|| omega < 0.25 || omega > 0.5 || beta3*Pop > 0.95*gamma2

IG = [beta1; beta2 ; beta3; gamma1; gamma2; rho1; rho2; psi; omega];
options = optimset('MaxFunEvals', 1200*length(IG),'MaxIter', 1200*length(IG), 'TolFun',5e-16,'TolX',5e-16);
[k_val,f_val]=fminsearch(@seir_parameter_fit,IG,options)

toc
%% 
figure
GraphIt(k_val)

delta  = 1/9; 

beta1  = k_val(1);
beta2  = k_val(2);
beta3  = k_val(3);
gamma1 = k_val(4);
gamma2 = k_val(5);
rho1   = k_val(6);
rho2   = k_val(7);
psi    = k_val(8);
omega  = k_val(9);

R0 = (1/(gamma1+psi))*(beta1*Pop + (beta3*Pop*psi/gamma2) + (beta2*Pop*rho1*gamma1/omega))
