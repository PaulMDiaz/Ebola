function [alpha,W,Lambda] = activity_score(dfun, m, n, dim , a , b)
% dfun is anonymous function that returns a gradient
% m is number of input parameters
% n number of points per dimension
% dim is dimension of active subspace
% alpha is the activitiy score 
% a is an mx1 vector containing the lower bound on a parameter's domain of
% definition
% b is an mx1 vector containing the upper bound on a parameter's domain of
% definition

s = [];
for i=1:m
    s = [s; parameter()];
end

order = n*ones(1,m);
[p,w] = gaussian_quadrature(s, order);
N = size(w,1);
C = zeros(m,m);
for i=1:N
    g = dfun(p(i,:),a,b);
    C = C + (g*g')*w(i);
end

% eigenvalues and eigenvectors
[W,Lambda] = eig(C);
[lambda, ind] = sort(diag(Lambda), 'descend');
W = W(:,ind);

% compute activity scores
alpha = W(:,1:dim).^2*lambda(1:dim);
alpha = alpha./sum(alpha);
