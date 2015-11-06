clear;clc;
tic
%LIBERIA

%Manner in which to order variables
%alpha, beta1, beta2, beta3, delta, gamma1, gamma2, psi, rho1, rho2, omega

beta1 = 1.0e-06;
beta2 = 2.0e-07; 
beta3 = 1.0e-8;
gamma2 = 1/20;
psi = 1/2;

Pop = 4.294e6;

IG = [beta1; beta2 ; beta3; gamma2; psi];
options = optimset('MaxFunEvals', 400*length(IG), 'TolFun',1e-6,'TolX',1e-6);
%IG = [delta,gamma1,psi];
k_val=fminsearch(@seir_parameter_fitSP,IG,options);
scaled_k_val  = [k_val(1:3)*Pop; k_val(4:5)];
%GraphItSP(IG)
toc
%% 
scaled_k_val
figure
GraphItSP(k_val)
