%LIBERIA
clear;clc;

% import .csv value
%filename= '/Users/redpatrioteq2/Dropbox/PIC_Math_group/ebola_modeling_484/ebola_data/matlab_data.csv';
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
tSpan = lib_data(:,1); %Final time is 193
tSpan = linspace(193,243, 20);

% parameter initialization
Pop = 4.28e6;
alpha  = 0;%0.024/365; % population growth constant (known empircally) 
beta1  = 0.376/Pop;%0.718/Pop;
beta2  = 0.135/Pop;%1.34/Pop;
beta3  = 0.163/Pop;%1.05e-5/Pop;
delta  = 1/9;%1/21;%1/incubation period for Ebola Virus
gamma1 = 0.0542;%1/10;
gamma2 = 0.174;%0.0524;
rho2   = 0.98;%0.7;
rho1   = 0.98;%1.1*rho2;
omega= 0.325;%1;  %population death constant 

% using initial condition from Liberia on Last Day of Data
y0 = [4.28e6, 1322.03, 248.80, 591.83, 36.59, 2799.88, 341.88 ];
imin = 0.5;
istep = 0.1;
imax = 2;
ind = 1;
psi0 = 0.5;    
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
I = imin:istep:imax;
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
% legend(legendary,'Location','EastOutside')
title('Liberia','Interpreter','LaTex','FontSize',14)
%%
print('Liberia_future_SA-11-04','-dpdf','-r300');
