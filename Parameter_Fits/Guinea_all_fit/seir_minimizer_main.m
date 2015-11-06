%GUINEA - parameter fitting with all 
clear;clc;
tic


Pop = 11.76e6;
%Initial guesses for parameters

beta1 = 0.3152/Pop;%0.2/Pop;
beta2 = 0.16/Pop;%0.2/Pop; 
beta3 = 0.0165/Pop;%0.08/Pop;
gamma2 = 0.016;%0.1747;%0.0513;
psi = 0.4999;
omega = 0.3;
gamma1 = 0.295;
rho1 = 0.98;
rho2 = 0.93;


initR0 = (1/(gamma1+psi))*(beta1*Pop + (beta3*Pop*psi/gamma2) + (beta2*Pop*rho1*gamma1/omega))

beta1*Pop < 0.08 || beta1*Pop > 0.4 ||beta2*Pop < 0.08 || beta2*Pop > 0.9 || beta3*Pop < 1e-2 ||  beta3*Pop > 0.2...
|| gamma1 < 0.05 || gamma1 > 0.32 || gamma2 < 0.01 || gamma2 > 0.4...
|| rho1 < 0.57 || rho1 > 0.98 || rho2 < 0.57 || rho2 > 0.93...
|| psi < 0.01 || psi > 0.5|| omega < 0.25 || omega > 0.5...
|| beta1*Pop + (beta2*Pop*rho1*gamma1/omega) - (gamma1*beta3*Pop/gamma2) < 0 || beta3*Pop*gamma1/gamma2 - beta1*Pop > 0

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
