%GUINEA
clear;clc;

% import .csv value
filename= 'matlab_data.csv';
%%%%%%%%%%%%%%% MODIFY FILE PATH TO DATA DESIRED ABOVE %%%%%%%%%%%%%%% 
fid = fopen(filename);
raw_data=textscan(fid, '%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f','delimiter',',');
fclose(fid);

% data of the form [day of outbreak, cases (i.e. infected) on this day, deaths
% (i.e. removed) by this day]
lib_data2 = [raw_data{2}, raw_data{3}, raw_data{11}];
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

%Adjusting from cumulative case counts to current
lib_data(:,2) = lib_data(:,2) - lib_data(:,3);
A = lib_data(2:72,:);
B=lib_data(74:end,:);
lib_data=[A;B];
tSpan = lib_data(:,1); %Final time is 289
tSpan = linspace(289, 349, 120);

% parameter initialization
Pop = 11.76e6;
alpha  = 0;%0.024/365; % population growth constant (known empircally) 
beta1  = 0.307/Pop;
beta2  = 0.477/Pop;
beta3  = 0.0325/Pop;
delta  = 1/9;%1/21;%1/incubation period for Ebola Virus
gamma1 = 1/10;
gamma2 = 0.0268;
rho2   = 0.59;
rho1   = 1.1*rho2;
omega= 1.07;  %population death constant 

% using initial condition from Guinea on Last Day of Data
y0 = [11.76e6, lib_data(end,2), lib_data(end,2), 0, lib_data(end,3), 0, 0 ];

figure;

ind = 1;
for i = 0.5:0.1:2
    psi = 1.08*i;    
    R0(ind) = (1/(gamma1+psi))*(beta1*Pop + (beta3*Pop*psi/gamma2) + (beta2*Pop*rho1*gamma1));

    odejac = @(t,u,up) jac(u, alpha, beta1, beta2, beta3, delta, gamma1, gamma2, psi, rho1, rho2, omega); 
    odefun =@(t,u) SEIHRRR(t, u, alpha, beta1, beta2, beta3, delta, gamma1, gamma2, psi, rho1, rho2, omega);
    opts = odeset('Jacobian', odejac);
    %[t,y] = ode23s(odefun, tSpan, y0, opts);
    [t,y] = ode15s(odefun, tSpan, y0,opts);

    %plot(t, y(:,3)+y(:,4))

    Infected(:,ind) = y(:,3)+y(:,4);
    plot(t, log10(Infected(:,ind)),'--')
    hold on
    ind = ind+1;
end
plot(t, log10(Infected(:,6)),'r')
axis([289,349,2.8,3.8])
xlabel('$t$ (days)','Interpreter','LaTex','FontSize',14), ylabel('$\log_{\ 10}I(t)$','Interpreter','LaTex','FontSize',14)
title('Guinea','Interpreter','LaTex','FontSize',14)