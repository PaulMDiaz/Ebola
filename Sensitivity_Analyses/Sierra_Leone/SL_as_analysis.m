
%%
clear all;close all;clear;clc;
M = 5000;
m = 8;
X = 2*rand(M,m)-1;
G = zeros(m,M);
F = zeros(M,1);
a = [0.1,0.1,0.05,.41,0.0275,0.1236,0.25,1/12];
b = [0.4,0.4,0.20, 1,0.1569,0.3840, 0.5,0.7];

for i=1:M
    G(:,i) = dR0(X(i,:), a, b);
    F(i) = R0(X(i,:), a , b);
end
[W,Sig,~] = svd(G,'econ');
evals = (1/M)*diag(Sig).^2;

save('sl_as_analysis')

%% Active Subspace Eigenvalues
clear all; close all; clear; clc;
load('activity_scores_gq')
evals = flipud(diag(Lambda));
f1 =figure;
hold on;
box on;
% C = [ 0.25 0.75 0.25;
%       0 0 1;
%       0 1 0;
%       0 1 1;
%       1 0 0;
%       1 0 1;
%       1 1 0;
%       1 0.5 0;
%     ];

for ind=1:m
plot(ind,log(evals(ind)),'b--*','LineWidth',2,'MarkerSize',20);
end
%plot(1:m,log(evals),'--*','LineWidth',2,'MarkerSize',20);
set(gca,'FontSize',20);
xlabel('i=1, \ldots, 8','Interpreter','latex'); xlim([1 m]);
ylabel('$log(\lambda_i)$','Interpreter','latex');

%%
print('sl_as_eigs','-dpdf','-r300');
%%
whitebg(f1)
set(gcf, 'InvertHardCopy','off');
print('sl_as_eigs_inv','-dpdf','-r300');
%% Active Subspace Eigenvectors
clear all;close all; clear; clc;
load('activity_scores_gq');
w1 = W(:,1);
f1 = figure;
ylim([-1 1]);
box on;
hold on;
for ind = 1:m
plot(ind,(-1)*w1(ind),'b*','LineWidth',2,'MarkerSize',20);
end
ax = gca;
xlim([1 8])
ax.XTick = [1 2 3 4 5 6 7 8];
ax.XTickLabel = {'\beta_1','\beta_2','\beta_3','\rho_1','\gamma_1','\gamma_2','\omega','\psi'};
set(ax,'Fontsize',20)
ylabel('$\mathbf{w}_1$','Interpreter','latex');

%% Active Subspace Eigenvectors
% clear all;close all; clear; clc;
% load('activity_scores_gq');  
% w1 = W(:,1);
% f1 = figure;
% box on;
% hold on;
% set(gca,'FontSize',20);
% C = [ 0.25 0.75 0.25;
%       0 0 1;
%       0 1 0;
%       0 1 1;
%       1 0 0;
%       1 0 1;
%       1 1 0;
%       1 0.5 0;
%     ];
% hold on;
% for ind = 1:m
% plot(ind,(-1)*w1(ind),'*','LineWidth',2,'MarkerSize',20,'MarkerFaceColor',C(ind,:),'markeredgecolor',C(ind,:));
% end
% xlabel('$i$th component','Interpreter','latex'); xlim([1 m]);
% ylabel('$\mathbf{w}_1$','Interpreter','latex');
% h = legend('$\beta_1$','$\beta_2$','$\beta_3$','$\rho_1$','$\gamma_1$', '$\gamma_2$','$\omega$','$\psi$','location','bestoutside');
% set(h,'Interpreter','latex');
% %title('AS 1st eigenvector, sl');
% ylim([-1,1]);

%%
whitebg(f1)
set(gcf, 'InvertHardCopy','off');
print('sl_as_evec_inv','-dpdf','-r300');
%%
print('sl_as_evec','-dpdf','-r300');


%% Sufficient Summary plot 
clear all;close all; clear; clc;
load('activity_scores_gq')
f1 = figure;
box on;
%hold on;
%plot(X*(-1)*W(:,1),F,'bo','MarkerFaceColor','r','MarkerSize',12);
%for i =1:length(X)
%scatter(X(i,:)*(-1)*W(:,1),X(i,:)*(-1)*W(:,2),F(i));
scatter(X*(-1)*W(:,1),X*(-1)*W(:,2),80,F,'filled'),colorbar;
%C = repmat([0.25,0.5,1],numel(X*(-1)*W(:,1)),1);
%c = C(:);
%surf(X*(-1)*W(:,1),X*(-1)*W(:,2),F);
%scatter3(X*(-1)*W(:,1),X*(-1)*W(:,2),F,80,F,'filled');
%axis square; grid on;
%end
set(gca,'FontSize',20);
t = title('$R_0(\mathbf{x})$');
x = xlabel('$\mathbf{w}_1^T\mathbf{x}$');
y = ylabel('$\mathbf{w}_2^T\mathbf{x}$');
%z = zlabel('$R_0(\mathbf{x})$');
set(x,'Interpreter','latex');
set(y,'Interpreter','latex');
%set(z,'Interpreter','latex');
set(t,'Interpreter','latex');

%%
whitebg(f1)
set(gcf, 'InvertHardCopy','off');
print('sl_as_ss_inv','-dpdf','-r300');
%%
print('sl_as_ss','-dpdf','-r300');

%% Activity Scores 
clear all;close all; clear; clc; 
load('activity_scores_gq');  
as = alpha;
%load('sl_as_analysis')
f1 = figure;
% as_dim = 2;
% as = zeros(m,1);
% for i = 1:as_dim
%    as = as + W(:,i).^2*evals(i);
% end
%as = W(:,2).^2*evals();
hold on;
box on;

for ind=1:m
plot(ind,as(ind),'b*','LineWidth',2,'MarkerSize',20);
end
ax = gca;
ylim([-0.05 1])
xlim([1 8])
ax.XTick = [1 2 3 4 5 6 7 8];
ax.XTickLabel = {'\beta_1','\beta_2','\beta_3','\rho_1','\gamma_1','\gamma_2','\omega','\psi'};
set(ax,'Fontsize',20)
ylabel('$\alpha(2)$','Interpreter','latex');

%%
set(gca,'FontSize',20);
xlabel('Parameter','Interpreter','latex');
ylabel('$\alpha(2)$','Interpreter','latex');
h = legend('$\beta_1$','$\beta_2$','$\beta_3$','$\rho_1$','$\gamma_1$', '$\gamma_2$','$\omega$','$\psi$','location','bestoutside');
set(h,'Interpreter','latex');
%%
whitebg(f1)
set(gcf, 'InvertHardCopy','off');
print('sl_as_act_inv','-dpdf','-r300');
%%
print('sl_as_act','-dpdf','-r300');


%% Gaussian Tensor Quadrature activity scores .
clear; close all; clc; clear all;
dfun = @dR0;
 as_dim = 2;
 m = 8; 
 NN = 8; % # of quadrature points per dimension
a = [0.1,0.1,0.05,.41,0.0275,0.1236,0.25,1/12];
b = [0.4,0.4,0.20, 1,0.1569,0.3840, 0.5,0.7];
 %%
 [alpha,W,Lambda] = activity_score(dfun, m, NN, as_dim,a,b);
 
 %%
 M = 5000;
 X = 2*rand(M,8)-1;
 F = zeros(M,1); 
 for i=1:M
    G(:,i) = dR0(X(i,:), a, b);
    F(i) = R0(X(i,:), a , b);
 end

 %%
save('activity_scores_gq');       
