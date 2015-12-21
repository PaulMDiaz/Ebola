close all; clear all; clear;clc;

m = 8;
M = 1000;
X = 2*rand(M,m)-1;
G = zeros(m,M);
for i=1:M
    G(:,i) = dR0(X(i,:));
end
[W,Sig,~] = svd(G,'econ');
evals = (1/M)*diag(Sig).^2;

save('liberia_active_subspace_analysis')