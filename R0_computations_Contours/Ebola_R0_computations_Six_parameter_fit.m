%Liberia
%Pop = 4.294e6;
clear;clc;

b1 = 0.376;
b2 = 0.135;
b3 = 0.163;
g1 = 0.0542;
g2 = 0.174;
r1 = 0.98;
r2 = 0.88;  %0.98
psi = 0.5;
omega = 0.325;

R0 = (1/(g1+psi))*(b1 + (b3*psi/g2) + (b2*r1*g1/omega))

R0pconst = b1 + (b2*r1*g1/omega) - (g1*b3/g2);
R0primepsi = - R0pconst/(g1 + psi)^2
R0at0 = (b1/g1) + (b2*r1/omega)
R0atinf = b3/g2

%R0tilde = (1/(g1+psi))*(b1 + (b3*psi/g2))

x = linspace(0.3,3,100);
for i = 1:length(x)
    yL(i) = (1/(g1+x(i)))*(b1 + (b3*x(i)/g2) + (b2*r1*g1/omega));
    yLp(i) = - R0pconst/(g1 + x(i))^2;
end
figure;
plot(x,yL, x, yLp,'--');
legend('R_0(\psi)', 'R_0\prime (\psi)');

psiL = psi;
R0L = R0;
R0pL = R0primepsi;

%%
%Guinea
%Pop = 11.76e6;
%clear;clc;

b1 = 0.315;
b2 = 0.16;
b3 = 0.0165;
g1 = 0.295;
g2 = 0.016;
r1 = 0.98;
r2 = 0.93;
psi = 0.5;
omega = 0.3;

R0 = (1/(g1+psi))*(b1 + (b3*psi/g2) + (b2*r1*g1/omega))

R0pconst = b1 + (b2*r1*g1/omega) - (g1*b3/g2);
R0primepsi = - R0pconst/(g1 + psi)^2
R0at0 = (b1/g1) + (b2*r1/omega)
R0atinf = b3/g2

%R0tilde = (1/(g1+psi))*(b1 + (b3*psi/g2))

x = linspace(0.3,3,100);
for i = 1:length(x)
    yG(i) = (1/(g1+x(i)))*(b1 + (b3*x(i)/g2) + (b2*r1*g1/omega));
    yGp(i) = - R0pconst/(g1 + x(i))^2;
end
figure;
plot(x,yG, x, yGp,'--');
legend('R_0(\psi)', 'R_0\prime (\psi)');

psiG = psi;
R0G = R0;
R0pG = R0primepsi;

%%
%Sierra Leone
%Pop = 6.092e6;
%clear;clc;

b1 = 0.251;
b2 = 0.395;
b3 = 0.079;
g1 = 0.051;
g2 = 0.0833;
r1 = 0.766;
r2 = 0.74;  %0.885;
psi = 0.442;
omega = 0.370;

R0 = (1/(g1+psi))*(b1 + (b3*psi/g2) + (b2*r1*g1/omega))

R0pconst = b1 + (b2*r1*g1/omega) - (g1*b3/g2);
R0primepsi = - R0pconst/(g1 + psi)^2
R0at0 = (b1/g1) + (b2*r1/omega)
R0atinf = b3/g2

%R0tilde = (1/(g1+psi))*(b1 + (b3*psi/g2))

x = linspace(0.3,3,100);
for i = 1:length(x)
    ySL(i) = (1/(g1+x(i)))*(b1 + (b3*x(i)/g2) + (b2*r1*g1/omega));
    ySLp(i) = - R0pconst/(g1 + x(i))^2;
end
figure;
plot(x,ySL, x, ySLp,'--');
legend('R_0(\psi)', 'R_0\prime (\psi)');

psiSL = psi;
R0SL = R0;
R0pSL = R0primepsi;

%%
plot(x, yL, x, yG,'--', x, ySL,'-.');
hold on
plot(psiL,R0L,'o','MarkerSize',12)
plot(psiG,R0G,'og','MarkerSize',12)
plot(psiSL,R0SL,'or','MarkerSize',12)
%title('Basic Reproduction Number as a function of hospitalization rate')
legend('Liberia', 'Guinea', 'Sierra Leone');
axis([0.3,3,0.9,2])
xlabel('$\psi$','Interpreter','LaTex','FontSize',16), 
ylabel('$R_0$','Interpreter','LaTex','FontSize',16)
%%
figure;
a = 5:100;
plot(x(a), yLp(a), x(a), yGp(a), '--', x(a), ySLp(a),'-.');
hold on
plot(psiL,R0pL,'o','MarkerSize',12)
plot(psiG,R0pG,'og','MarkerSize',12)
plot(psiSL,R0pSL,'or','MarkerSize',12)

%title('Derivative of Basic Reproduction Number')
legend('Liberia', 'Guinea', 'Sierra Leone');
xlabel('$\psi$','Interpreter','LaTex','FontSize',14), ylabel('$\frac{\partial R_0}{\partial \psi}$','Interpreter','LaTex','FontSize',20)
