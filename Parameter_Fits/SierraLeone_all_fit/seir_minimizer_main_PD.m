%SIERRALEONE
clear;clc;
tic

%Manner in which to order variables
%alpha, beta1, beta2, beta3, delta, gamma1, gamma2, psi, rho1, rho2, omega

Pop = 6.092e6;
%Initial guesses for parameters
beta1 = 0.133/Pop; %0.128/Pop;
beta2 = 0.176/Pop; %0.111/Pop;
beta3 = 0.156/Pop;
gamma2 = 0.1706; %0.3011; %0.01355;
psi = 1; %0.106; %0.0478;
omega = 0.5; %0.222; 0.333;

gamma1 = 0.055;%0.0585; %0.068;
rho2 = 0.7985;%0.7084; %0.884;
rho1 = 0.7097;%0.7265; %0.853;

initR0 = (1/(gamma1+psi))*(beta1*Pop + (beta3*Pop*psi/gamma2) + (beta2*Pop*rho1*gamma1/omega))

IG = [beta1; beta2; beta3; gamma2; psi; omega];
options = optimset('MaxFunEvals', 800*length(IG),'MaxIter', 800*length(IG), 'TolFun',1e-15,'TolX',1e-15);
[k_val,f_val]=fminsearch(@seir_parameter_fit,IG,options)
toc
%% 
figure
GraphIt(k_val)

beta1 = k_val(1);
beta2 = k_val(2);
beta3 = k_val(3);
delta  = 1/9; %1/21;%1/incubation period for Ebola Virus
gamma1 = 0.055;%0.0585; %0.068;
rho2 = 0.7985;%0.7084; %0.884;
rho1 = 0.7097;%0.7265; %0.853;
gamma2 = k_val(4);
psi    = k_val(5);
omega = k_val(6);

R0 = (1/(gamma1+psi))*(beta1*Pop + (beta3*Pop*psi/gamma2) + (beta2*Pop*rho1*gamma1/omega))
