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

A = lib_data(2:72,:);
B=lib_data(74:end,:);
lib_data=[A;B];
tSpan = lib_data(:,1); %Final time is 289
tSpan = linspace(289, 439, 50);

% parameter initialization
Pop = 11.76e6;
alpha = 0;
beta1  = 0.315/Pop;
beta2  = 0.16/Pop;
beta3  = 0.0165/Pop;
delta  = 1/9;%1/21;%1/incubation period for Ebola Virus
gamma1 = 0.295;
gamma2 = 0.016;
rho1   = 0.98;
rho2   = 0.93;
omega  = 0.3;  %population death constant 

% using initial condition from Guinea on Last Day of Data
y0 = [11.76e6, 0, lib_data(end,2), 0, 0, lib_data(end,3), 0 ];

ind = 1;
psi0 = 0.5;
for i = 1:0.04:2
    psi = psi0*i;    
    R0(ind) = (1/(gamma1+psi))*(beta1*Pop + (beta3*Pop*psi/gamma2) + (beta2*Pop*rho1*gamma1));

    odejac = @(t,u,up) jac(u, alpha, beta1, beta2, beta3, delta, gamma1, gamma2, psi, rho1, rho2, omega); 
    odefun =@(t,u) SEIHRRR(t, u, alpha, beta1, beta2, beta3, delta, gamma1, gamma2, psi, rho1, rho2, omega);
    opts = odeset('Jacobian', odejac);
    %[t,y] = ode23s(odefun, tSpan, y0, opts);
    [t,y] = ode15s(odefun, tSpan, y0,opts);

    %plot(t, y(:,3)+y(:,4))

    Infected(:,ind) = y(:,3);
        %Percentage change
        %plot(t, Infected(:,ind),'--') 
        
        %Finite Difference
        Diff = (Infected(:,ind)-Infected(:,1))/(psi-psi0);
        plot(t, log10(-Diff),'--b') 
        if ind == 2
            Diff2 = Diff;
        end
    hold on
    ind = ind+1;
end
    fill([t,t], log10(-[Diff2,Diff]), '-r')
    %Percentage change
    %plot(t, Infected(:,1)/lib_data(end,2),'r')
axis([290,340, 2, 3])    
xlabel('t (days)','Interpreter','LaTex','FontSize',14), 
%ylabel('Change in $I(t)$','Interpreter','LaTex','FontSize',14)
ylabel('$\log_{\ 10} \left (\frac{\partial I}{\partial \psi} \right )$','Interpreter','LaTex','FontSize',14)
title('Guinea','Interpreter','LaTex','FontSize',14)