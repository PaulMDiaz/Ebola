%LIBERIA - parameter fitting with all 
clear;clc;
tic

%Manner in which to order variables
%alpha, beta1, beta2, beta3, delta, gamma1, gamma2, psi, rho1, rho2, omega

% beta1 = 1.0e-06;
% beta2 = 2.0e-07; 
% beta3 = 1.0e-8;
% gamma2 = 1/20;
% psi = 1/2;

Pop = 4.294e6;
%Initial guesses for parameters
beta1 = 0.31/Pop;%0.2/Pop;
beta2 = 0.2/Pop;%0.2/Pop; 
beta3 = 0.15/Pop;%0.08/Pop;
gamma2 = 0.171;%0.1747;%0.0513;
psi = 0.44;
omega = 0.3;
gamma1 = 0.0689;%0.0875;%0.0585; %0.068;
rho1 = 0.96;%0.7268; %0.853;
rho2 = 0.93;%0.7084; %0.884;

initR0 = (1/(gamma1+psi))*(beta1*Pop + (beta3*Pop*psi/gamma2) + (beta2*Pop*rho1*gamma1/omega))

beta1*Pop < 0.08 || beta1*Pop > 0.4 ||beta2*Pop < 0.08 || beta2*Pop > 0.5 || beta3*Pop < 1e-2 ||  beta3*Pop > 0.2...
|| gamma1 < 0.03 || gamma1 > 0.2 || gamma2 < 0.15 || gamma2 > 0.4...
|| rho1 < 0.75 || rho1 > 0.98 || rho2 < 0.7 || rho2 > 0.98...
|| psi < 0.04 || psi > 0.5|| omega < 0.25 || omega > 0.5 || beta3*Pop > 0.95*gamma2

IG = [beta1; beta2 ; beta3; gamma1; gamma2; rho1; rho2; psi; omega];
options = optimset('MaxFunEvals', 800*length(IG),'MaxIter', 800*length(IG), 'TolFun',5e-16,'TolX',5e-16);
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
