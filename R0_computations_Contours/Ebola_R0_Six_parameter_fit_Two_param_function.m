%Liberia
%Pop = 4.294e6;
clear;clc;

beta1 = 0.376;
beta2 = 0.135;
beta3 = 0.163;
gamma2 = 0.174;
psi = 0.5;
omega = 0.325;
gamma1 = 0.0542;
rho1 = 0.98;
rho2 = 0.88;  %0.98

R0 = (1/(gamma1+psi))*(beta1 + (beta3*psi/gamma2) + (beta2*rho1*gamma1/omega))

R0pconst = beta1 + (beta2*rho1*gamma1/omega) - (gamma1*beta3/gamma2);
R0primepsi = - R0pconst/(gamma1 + psi)^2
R0primeomega = -(1/(gamma1 + psi))*(beta2*rho1*gamma1)/omega^2
gradR0 = sqrt(R0primepsi^2 + R0primeomega^2)
R0at0 = (beta1/gamma1) + (beta2*rho1/omega)
R0atinf = beta3/gamma2

x = psi*linspace(0.5,4,100);
y = omega*linspace(0.5,4,100);
for i = 1:length(x)
    for j = 1:length(y)
        yL(i,j) = (1/(gamma1+x(i)))*(beta1 + (beta3*x(i)/gamma2) + (beta2*rho1*gamma1/y(j)));
    end
end

figure;
C = contourf(x,y,yL'); colorbar; clabel(C)
%surf(x,y,yL)
xlabel('$\psi$','Interpreter','LaTex','FontSize',20), 
ylabel('$\omega$','Interpreter','LaTex','FontSize',20)
h = title(['\makebox[4in][c]{\textbf{Liberia}}', sprintf('\n'), '\makebox[4in][c]{$R_0(\psi,\omega)$}']);
set(h,'Interpreter','latex', 'FontSize', 20);
%zlabel('$R_0(\psi,\omega)$','Interpreter','LaTex','FontSize',20)
%h = title('\makebox[4in][c]{\textbf{Liberia}}');
%set(h,'Interpreter','latex', 'FontSize', 20);

hold on
plot(psi,omega,'ks','MarkerSize',12)
%plot3(psi,omega,R0,'ks','MarkerSize',20,'MarkerEdgeColor','none','MarkerFaceColor','k')
hold off

% numarrows = 20;
% [x,y] = meshgrid(psi*linspace(0.5,4,numarrows),omega*linspace(0.5,4,numarrows));
% for i = 1:numarrows
%     yLpsi(i) = - R0pconst/(g1 + x(i,i))^2;
%     yLomega(i) = -(1/(g1 + psi))*(b2*r1*g1)/y(i,i)^2;
% end
% 
% [yLpsi,yLomega] = meshgrid(yLpsi,yLomega);
% 
% figure;
% quiver(x,y,yLpsi,yLomega)
% xlabel('$\psi$','Interpreter','LaTex','FontSize',20), 
% ylabel('$\omega$','Interpreter','LaTex','FontSize',20)
% h = title(['\makebox[4in][c]{\textbf{Liberia}}', sprintf('\n'), '\makebox[4in][c]{$\nabla R_0(\psi,\omega)$}']);
% set(h,'Interpreter','latex', 'FontSize', 20);
% hold on
% plot(psi,omega,'ks','MarkerSize',12)
% hold off

% psiL = psi;
% omegaL = omega;
% R0L = R0;
% R0pL = R0primepsi;

%%
%Guinea
%Pop = 11.76e6;
clear;clc;

beta1 = 0.315;
beta2 = 0.16;
beta3 = 0.0165;
gamma1 = 0.295;
gamma2 = 0.016;
rho1 = 0.98;
rho2 = 0.93;
psi = 0.5;
omega = 0.3;

R0 = (1/(gamma1+psi))*(beta1 + (beta3*psi/gamma2) + (beta2*rho1*gamma1/omega))

R0pconst = beta1 + (beta2*rho1*gamma1/omega) - (gamma1*beta3/gamma2);
R0primepsi = - R0pconst/(gamma1 + psi)^2
R0primeomega = -(1/(gamma1 + psi))*(beta2*rho1*gamma1)/omega^2
gradR0 = sqrt(R0primepsi^2 + R0primeomega^2)
R0atpsi0 = (beta1/gamma1) + (beta2*rho1/omega)
R0atpsiinf = beta3/gamma2
R0atomegainf = (1/(gamma1+psi))*(beta1 + (beta3*psi/gamma2))

x = psi*linspace(0.5,4,100);
y = omega*linspace(0.5,4,100);
for i = 1:length(x)
    for j = 1:length(y)
        yG(i,j) = (1/(gamma1+x(i)))*(beta1 + (beta3*x(i)/gamma2) + (beta2*rho1*gamma1/y(j)));
        %R0pnum(j) = b1 + (b2*r1*g1/y(j)) - (g1*b3/g2);
        %yLp(i,j) = - R0pnum(j)/(g1 + x(i))^2;
    end
end

figure;
C = contourf(x,y,yG'); colorbar; clabel(C)
%surf(x,y,yG')
xlabel('$\psi$','Interpreter','LaTex','FontSize',20), 
ylabel('$\omega$','Interpreter','LaTex','FontSize',20)
h = title(['\makebox[4in][c]{\textbf{Guinea}}', sprintf('\n'), '\makebox[4in][c]{$R_0(\psi,\omega)$}']);
set(h,'Interpreter','latex', 'FontSize', 20);
%zlabel('$R_0(\psi,\omega)$','Interpreter','LaTex','FontSize',20)
%h = title('\makebox[4in][c]{\textbf{Guinea}}');
%set(h,'Interpreter','latex', 'FontSize', 20);

hold on
plot(psi,omega,'ks','MarkerSize',12)
%plot3(psi,omega,R0,'ks','MarkerSize',20,'MarkerEdgeColor','none','MarkerFaceColor','k')
hold off

% numarrows = 20;
% [x,y] = meshgrid(psi*linspace(0.5,4,numarrows),omega*linspace(0.5,4,numarrows));
% for i = 1:numarrows
%     yGpsi(i) = - R0pconst/(g1 + x(i,i))^2;
%     yGomega(i) = -(1/(g1 + psi))*(b2*r1*g1)/y(i,i)^2;
% end
% 
% [yGpsi,yGomega] = meshgrid(yGpsi,yGomega);
% 
% figure;
% quiver(x,y,yGpsi,yGomega)
% xlabel('$\psi$','Interpreter','LaTex','FontSize',20), 
% ylabel('$\omega$','Interpreter','LaTex','FontSize',20)
% h = title(['\makebox[4in][c]{\textbf{Guinea}}', sprintf('\n'), '\makebox[4in][c]{$\nabla R_0(\psi,\omega)$}']);
% set(h,'Interpreter','latex', 'FontSize', 20);
% hold on
% plot(psi,omega,'ks','MarkerSize',12)
% hold off

% psiG = psi;
% omegaG = omega;
% R0G = R0;
% R0pG = R0primepsi;

%%
%Sierra Leone
%Pop = 6.092e6;
clear;clc;

beta1 = 0.251;
beta2 = 0.395;
beta3 = 0.0791;
gamma1 = 0.051;
gamma2 = 0.0833;
rho1 = 0.766;
rho2 = 0.74;    %0.885;
psi = 0.442;
omega = 0.37;

R0 = (1/(gamma1+psi))*(beta1 + (beta3*psi/gamma2) + (beta2*rho1*gamma1/omega))

R0pconst = beta1 + (beta2*rho1*gamma1/omega) - (gamma1*beta3/gamma2);
R0primepsi = - R0pconst/(gamma1 + psi)^2
R0primeomega = -(1/(gamma1 + psi))*(beta2*rho1*gamma1)/omega^2
gradR0 = sqrt(R0primepsi^2 + R0primeomega^2)
R0at0 = (beta1/gamma1) + (beta2*rho1/omega)
R0atinf = beta3/gamma2

x = psi*linspace(0.5,4,100);
y = omega*linspace(0.5,4,100);
for i = 1:length(x)
    for j = 1:length(y)
        ySL(i,j) = (1/(gamma1+x(i)))*(beta1 + (beta3*x(i)/gamma2) + (beta2*rho1*gamma1/y(j)));
        %R0pnum(j) = b1 + (b2*r1*g1/y(j)) - (g1*b3/g2);
        %yLp(i,j) = - R0pnum(j)/(g1 + x(i))^2;
    end
end

figure;
C = contourf(x,y,ySL'); colorbar; clabel(C)
%surf(x,y,ySL')
xlabel('$\psi$','Interpreter','LaTex','FontSize',14), 
ylabel('$\omega$','Interpreter','LaTex','FontSize',20)
h = title(['\makebox[4in][c]{\textbf{Sierra Leone}}', sprintf('\n'), '\makebox[4in][c]{$R_0(\psi,\omega)$}']);
set(h,'Interpreter','latex', 'FontSize', 20);
%zlabel('$R_0(\psi,\omega)$','Interpreter','LaTex','FontSize',20)
%h = title('\makebox[4in][c]{\textbf{Sierra Leone}}');
%set(h,'Interpreter','latex', 'FontSize', 20);

hold on
plot(psi,omega,'ks','MarkerSize',12)
%plot3(psi,omega,R0,'ks','MarkerSize',20,'MarkerEdgeColor','none','MarkerFaceColor','k')
hold off


% numarrows = 20;
% [x,y] = meshgrid(psi*linspace(0.5,4,numarrows),omega*linspace(0.5,4,numarrows));
% for i = 1:numarrows
%     ySLpsi(i) = - R0pconst/(g1 + x(i,i))^2;
%     ySLomega(i) = -(1/(g1 + psi))*(b2*r1*g1)/y(i,i)^2;
% end
% 
% [ySLpsi,ySLomega] = meshgrid(ySLpsi,ySLomega);
% 
% figure;
% quiver(x,y,ySLpsi,ySLomega)
% xlabel('$\psi$','Interpreter','LaTex','FontSize',20), 
% ylabel('$\omega$','Interpreter','LaTex','FontSize',20)
% h = title(['\makebox[4in][c]{\textbf{Sierra Leone}}', sprintf('\n'), '\makebox[4in][c]{$\nabla R_0(\psi,\omega)$}']);
% set(h,'Interpreter','latex', 'FontSize', 20);
% hold on
% plot(psi,omega,'ks','MarkerSize',12)
% hold off

% psiSL = psi;
% omegaSL = omega;
% R0SL = R0;
% R0pSL = R0primepsi;

