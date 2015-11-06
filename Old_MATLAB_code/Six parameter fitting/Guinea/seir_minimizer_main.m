clear;clc;
tic
%GUINEA - parameter fitting with omega

Pop = 11.76e6;
% beta1 = 1.0e-07;
% beta2 = 2.0e-07; 
% beta3 = 1.0e-08;
% gamma2 = 1/30;
% psi = 1/2;

%Initial guesses for parameters
beta1 = 0.316/Pop;
beta2 = 0.446/Pop; 
beta3 = 0.0325/Pop;
gamma2 = 0.0269;
psi = 1.11;

omega = 1;


IG = [beta1; beta2 ; beta3; gamma2; psi; omega];
options = optimset('MaxFunEvals', 400*length(IG), 'TolFun',1e-6,'TolX',1e-6);
k_val=fminsearch(@seir_parameter_fit,IG,options);

toc
%% 
figure
GraphIt(k_val)

beta1 = k_val(1);
beta2 = k_val(2);
beta3 = k_val(3);
delta  = 1/9; %1/21;%1/incubation period for Ebola Virus
gamma1 = 1/10;
gamma2 = k_val(4);
psi    = k_val(5);
rho2   = 0.59;
rho1   = 1.1*rho2;
omega = k_val(6);

R0 = (1/(gamma1+psi))*(beta1*Pop + (beta3*Pop*psi/gamma2) + (beta2*Pop*rho1*gamma1/omega))
