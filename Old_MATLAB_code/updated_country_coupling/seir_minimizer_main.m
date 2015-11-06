clear;clc;
tic
%To run the parameter fit execute this block of code, but first you must
%change the filename in the seir_parameter_fit_3_30.m 
%This returns fitted parameter values for everything in our model
%The output very much so depends on the initial k0, I have chosen the
%fitted k1,k2 values and then our theoretical values for 

%Manner in which to order variables
%alpha, beta1, beta2, beta3, delta, gamma1, gamma2, psi, rho1, rho2, omega

%beta1 = 1.0e-07;
%beta2 = 2.0e-07; 
%beta3 = 1.0e-8;
%gamma2 = 1/3;
%psi = 0.1;


%IG = [beta1; beta2 ; beta3;gamma2;  psi];
IG = [.01 .01 .01 .01 .01 .01] % in the order psi12 13 21 23 31 32
%IG = [0 0 0 0 0 0]
%IG = [delta,gamma1,psi];
k_val=fminsearch(@seir_parameter_fit,IG,optimset('Display','iter'));
%GraphIt(IG)
toc
%% 
figure
%GraphIt(k_val)
k_val

%%
final_graphing(k_val)
