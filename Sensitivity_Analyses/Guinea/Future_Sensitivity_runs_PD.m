%GUINEA
clear;clc;
format short e
% import .csv value
filename = 'matlab_data.csv';
%filename= '/Users/redpatrioteq2/Dropbox/PIC_Math_group/ebola_modeling_484/ebola_data/matlab_data.csv';
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
A = lib_data(2:72,:);
B=lib_data(74:end,:);
lib_data=[A;B];
tSpan = lib_data(:,1); %Final time is 289
tSpan = linspace(289, 339, 20);

% parameter initialization
Pop = 1.18e7;
alpha  = 0;%0.024/365; % population growth constant (known empircally) 
beta1  = 0.315/Pop;%0.316/Pop;
beta2  = 0.16/Pop;%0.446/Pop;
beta3  = 0.0165/Pop;%0.0325/Pop;
delta  = 1/9;%1/21;%1/incubation period for Ebola Virus
gamma1 = 0.295;%1/10;
gamma2 = 0.016;%0.0269;
rho2   = 0.93;%0.59;
rho1   = 0.98;%1.1*rho2;
omega= 0.3;%1;  %population death constant 

% using initial condition from Guinea on Last Day of Data
y0 = [1.18e7, 135.73, 18.83, 430.42, 17.80, 1715.3, 82.38 ];
imin = 0.5;
istep = 0.1;
imax = 2;
ind = 1;
psi0 = 0.5;
Infected = zeros(length(tSpan),length([imin:istep:imax]));
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
    if(i == 1)
       relative_infected = Infected(end,ind);
    end
    ind = ind+1;
    
end
%%
%close all; 
%FigHandle = figure('Position', [100, 100, 1049, 895]);
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
%plot(t,Infected,'--','LineWidth',1)
%axis square;
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
title('Guinea','Interpreter','LaTex','FontSize',14)
%%
%surf(Infected/Infected(end,1));    
close all;
FigHandle = figure('Position', [100, 100, 1049, 895]);
surf([1:0.1:2]',t,Infected/Infected(end,1))
colormap jet
view(45,30)
xlabel('Multiplicative increase in $\psi$','Interpreter','LaTex','FontSize',14), 
ylabel('Day','Interpreter','LaTex','FontSize',14)
zlabel('Infected Proportion','Interpreter','LaTex','FontSize',14)
title('Change in Guinea Infected Population Based on Change in $\psi$','Interpreter','LaTex','FontSize',14)
n=get(gca,'Ztick');
colorbar('Eastoutside')
%set(gca,'zticklabel',sprintf('%.e |',n'));
%colorbar('Ticks',[-5,-2,1,4,7],...
%         'TickLabels',{'Cold','Cool','Neutral','Warm','Hot'})
%%
print('Guinea_future_SA-11-04','-dpdf','-r300');