%LIBERIA - Global Sensitivity
%clear;clc;
clearvars -except totsumG tG densG;clc;

% import .csv value
filename= 'matlab_data.csv';
%%%%%%%%%%%%%%% MODIFY FILE PATH TO DATA DESIRED ABOVE %%%%%%%%%%%%%%% 
fid = fopen(filename);
raw_data=textscan(fid, '%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f','delimiter',',');
fclose(fid);

% data of the form [day of outbreak, cases (i.e. infected) on this day, deaths
% (i.e. removed) by this day]
lib_data2 = [raw_data{2}, raw_data{4}, raw_data{12}];
temp_matrix = zeros(size(lib_data2));

count=1;
for i=1:length(lib_data2)
    lib_data2(i,1) = lib_data2(i,1); %- 102; % shift the data to start on this
                                           % (arbitrary!) date
    if (~isnan(lib_data2(i,2)) & ~isnan(lib_data2(i,3)))
        if (lib_data2(i,1) >= 0)
            temp_matrix(count,:) = lib_data2(i,:);
            count=count+1;
        end
    end
    
end
lib_data=flipud(temp_matrix(find(temp_matrix(:,1),1,'first'):find(temp_matrix(:,1),1,'last')+1,:));
lib_data = lib_data(25:57,:);
%tSpan = lib_data(:,1); %Final time is 193
tSpan = linspace(193, 400, 40); %tSpan = linspace(193,243, 20);

% parameter initialization
Pop = 4.294e6;
alpha  = 0;%0.024/365; % population growth constant (known empircally) 
beta1  = 0.376/Pop;
beta2  = 0.135/Pop;
beta3  = 0.163/Pop;
delta  = 1/9;
gamma1 = 0.0542;
gamma2 = 0.174;
rho1   = 0.98;
rho2   = 0.88; %0.98
omega  = 0.325;

% using initial condition from Liberia on Last Day of Data
y0 = [4.294e6, 0, lib_data(end,2), 0, 0, lib_data(end,3), 0];

N = 1000; %Number of trials
h = 1e-6; %Finite difference step size

%psi mean is 0.45; psi spread is 10 percent
psi0 = 0.5;
%ll = 0.9*psi0;
%ul = 1.1*psi0;

mn = psi0;
stddev = 0.34*mn;
v = stddev^2;

%mn = exp(loc + 0.5*scale^2)
%stddev = sqrt( (exp(scale^2) - 1)*exp(2*loc+scale^2)) 
con = 1 + (v/mn^2);
loc = log(mn/sqrt(con));
scale = sqrt(log(con));

pdflogN = @(z) (1./z)*(1/scale/sqrt(2*pi)).*exp(-(log(z)-loc).^2/(2*scale^2));


x = 0:1e-3:2;
densL = pdflogN(x);
figure;
plot(x, densL);

totsum = 0;

for trial = 1:N
%     draw = rand;
%     psi(trial) = ll + draw*(ul-ll); 
%     
%    drawplus = draw + h;
%    psiplus(trial) = ll + drawplus*(ul-ll); 

    draw = randn;
    psi(trial) = exp(loc + scale*draw);
    
    drawplus = draw + h;
    psiplus(trial) = exp(loc + scale*drawplus);
     
    %R0(ind) = (1/(gamma1+psi))*(beta1*Pop + (beta3*Pop*psi/gamma2) + (beta2*Pop*rho1*gamma1));
    
    %Solve ODE
    odejac = @(t,u,up) jac(u, alpha, beta1, beta2, beta3, delta, gamma1, gamma2, psi(trial), rho1, rho2, omega); 
    odefun = @(t,u) SEIHRRR(t, u, alpha, beta1, beta2, beta3, delta, gamma1, gamma2, psi(trial), rho1, rho2, omega);
    opts = odeset('Jacobian', odejac);
    [t,y] = ode15s(odefun, tSpan, y0, opts);
    Infected(:,trial) = y(:,3);
    
    %Solve Perturbed ODE
    odejac2 = @(t,u,up) jac(u, alpha, beta1, beta2, beta3, delta, gamma1, gamma2, psiplus(trial), rho1, rho2, omega); 
    odefun2 = @(t,u) SEIHRRR(t, u, alpha, beta1, beta2, beta3, delta, gamma1, gamma2, psiplus(trial), rho1, rho2, omega);
    opts2 = odeset('Jacobian', odejac2);
    [t,y2] = ode15s(odefun2, tSpan, y0, opts2);
    Infectedplus(:,trial) = y2(:,3);
    
    %Finite Difference
    DiffI = (Infected(:,trial)-Infectedplus(:,trial))/h;
    totsum = totsum + DiffI.^2*pdflogN(psi(trial));
    %plot(t, Diff,'--') 
    %hold on
end
    totsum = totsum/N;
    plot(t, totsum,'--'); 
    
    totsumL = totsum;
    tL = t;
    
xlabel('t (days)','Interpreter','LaTex','FontSize',14), 
%ylabel('Change in $I(t)$','Interpreter','LaTex','FontSize',14)
ylabel('$C(t)$','Interpreter','LaTex','FontSize',14)
title('Liberia','Interpreter','LaTex','FontSize',14)