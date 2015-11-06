%LIBERIA - parameter fitting with omega
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
beta1 = 0.133/Pop;%0.5/Pop;
beta2 = 0.176/Pop;%0.3/Pop; 
beta3 = 0.156/Pop;%0.047/Pop;
gamma2 = 0.1747;%0.0513;
psi = 0.25;
omega = 0.3;%0.832;

gamma1 = 0.0875;%0.0585; %0.068;
rho2 = 0.8523;%0.7084; %0.884;
rho1 = 0.9111;%0.7265; %0.853;

initR0 = (1/(gamma1+psi))*(beta1*Pop + (beta3*Pop*psi/gamma2) + (beta2*Pop*rho1*gamma1/omega))


IG = [beta1; beta2 ; beta3; gamma2; psi; omega];
options = optimset('MaxFunEvals', 800*length(IG),'MaxIter', 800*length(IG), 'TolFun',1e-15,'TolX',1e-15);
[k_val,f_val]=fminsearch(@seir_parameter_fit,IG,options)

toc
%% 
figure
GraphIt(k_val)

beta1 = k_val(1);
beta2 = k_val(2);
beta3 = k_val(3);
delta  = 1/9; 
gamma1 = 0.0875;%0.0585; %0.068;
rho2 = 0.8523;%0.7084; %0.884;
rho1 = 0.9111;%0.7268; %0.853;
gamma2 = k_val(4);
psi    = k_val(5);
omega = k_val(6);

R0 = (1/(gamma1+psi))*(beta1*Pop + (beta3*Pop*psi/gamma2) + (beta2*Pop*rho1*gamma1/omega))
