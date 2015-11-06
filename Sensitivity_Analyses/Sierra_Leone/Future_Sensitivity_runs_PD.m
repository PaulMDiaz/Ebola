%SIERRA LEONE
clear;clc;
% import .csv value
filename= 'matlab_data.csv';
%%%%%%%%%%%%%%% MODIFY FILE PATH TO DATA DESIRED ABOVE %%%%%%%%%%%%%%% 
fid = fopen(filename);
raw_data=textscan(fid, '%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f','delimiter',',');
fclose(fid);

% data of the form [day of outbreak, cases (i.e. infected) on this day, deaths
% (i.e. removed) by this day]
lib_data2 = [raw_data{2}, raw_data{5}, raw_data{13}];
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
lib_data = lib_data(21:66,:);
tSpan = lib_data(:,1); %Final time is 214
tSpan = linspace(214,264, 20);


% parameter initialization
Pop = 6.09e6;
alpha  = 0;%0.024/365; % population growth constant (known empircally) 
beta1  = 0.251/Pop;%0.803/Pop;
beta2  = 0.395/Pop;%5.35/Pop;
beta3  = 0.0791/Pop;%0.000169/Pop;
delta  = 1/9;%1/21;%1/incubation period for Ebola Virus
gamma1 = 0.0510;%/10;
gamma2 = 0.0833;%0.0271;
rho2   = 0.885;%0.41;
rho1   = 0.766;%1.1*rho2;
omega  = 0.37;%1;  %population death constant 

% using initial condition from Sierra Leone on Last Day of Data
y0 = [6.09e6, 280.67, 60.49, 252.49, 6.02, 765.88, 264.79 ];
imin = 0.5;
istep = 0.1;
imax = 2;
ind = 1;
psi0 = 0.442;
for i = imin:istep:imax
    psi = psi0*i;    
    R0(ind) = (1/(gamma1+psi))*(beta1*Pop + (beta3*Pop*psi/gamma2) + (beta2*Pop*rho1*gamma1));

    odejac = @(t,u,up) jac(u, alpha, beta1, beta2, beta3, delta, gamma1, gamma2, psi, rho1, rho2, omega); 
    odefun =@(t,u) SEIHRRR(t, u, alpha, beta1, beta2, beta3, delta, gamma1, gamma2, psi, rho1, rho2, omega);
    opts = odeset('Jacobian', odejac);
    [t,y] = ode15s(odefun, tSpan, y0,opts);
    intI(1) = y(1,3);
    for k = 2:length(y(:,2))
        intI(k) = trapz(t(1:k),y(1:k,3));
    end
    intI = intI';
    Infected(:,ind) = y(:,3)+y(:,4);
    ind = ind+1;
end
%%
close all; 
I = imin:istep:imax;
hold on;
box on;
for l = 1:length(I)
   if(I(l) ~= 1) 
    plot(t,Infected(:,l),'b--','LineWidth',1); 
   else 
    plot(t,Infected(:,l),'r-','LineWidth',2); 
   end
end 
xlabel('t (Days)','Interpreter','LaTex','FontSize',14), 
ylabel('I(t)','Interpreter','LaTex','FontSize',14)

% for k = 1:length(I)
%     if( I(k) ~= 1 && I(k) ~= 2 && I(k) ~= 3 && I(k) ~= 4 && I(k) ~= 5 )
%         legendary(k,:) = sprintf('%0.2g\\psi', I(k));
%     end     
%     if I(k) == 1 
%         legendary(k,:) = '1.0\psi';
%     end
%     if  I(k) == 2
%         legendary(k,:) = '2.0\psi';
%     end
%     if I(k) == 3
%         legendary(k,:) = '3.0\psi';
%     end
%     if I(k) == 4
%         legendary(k,:) = '4.0\psi';
%     end
%     if I(k) == 5
%         legendary(k,:) = '5.0\psi';
%     end
%         
%     
% end
%legend(legendary,'Location','EastOutside')
title('Sierra Leone','Interpreter','LaTex','FontSize',14)
    
%%
FigHandle = figure('Position', [100, 100, 1049, 895]);
surf([1:0.1:2]',t,Infected/Infected(end,1))
colormap jet
colorbar('Eastoutside')
view(45,30)
xlabel('Multiplicative increase in $\psi$','Interpreter','LaTex','FontSize',14), 
ylabel('Day','Interpreter','LaTex','FontSize',14)
zlabel('Infected Proportion','Interpreter','LaTex','FontSize',14)
title(' Change in Sierra Leone Infected Population Based on Change in $\psi$','Interpreter','LaTex','FontSize',14)
%%
print('SL_future_SA-11-04','-dpdf','-r300');