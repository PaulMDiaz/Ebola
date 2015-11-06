function err = final_graphing(IG)

% import .csv value
%filename= '/Users/redpatrioteq2/Dropbox/PIC_Math_group/ebola_data/matlab_data.csv';
filename='/home/eric/ebola/ebola_modeling_484/ebola_data/matlab_data.csv';
%%%%%%%%%%%%%%% MODIFY FILE PATH TO DATA DESIRED ABOVE %%%%%%%%%%%%%%% 
fid = fopen(filename);
raw_data=textscan(fid, '%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f','delimiter',',');
fclose(fid);

% data of the form [day of outbreak, cases (i.e. infected) on this day, deaths
% (i.e. removed) by this day]
lib_data2 = [raw_data{2}, raw_data{4}, raw_data{12}];
guin_data2 = [raw_data{2}, raw_data{3}, raw_data{11}];
sl_data2 = [raw_data{2}, raw_data{5}, raw_data{13}];

temp_matrix1 = zeros(size(lib_data2));
temp_matrix2 = zeros(size(guin_data2));
temp_matrix3 = zeros(size(sl_data2));

count1=1;
count2=1;
count3=1;
for i=1:length(lib_data2)
    lib_data2(i,1) = lib_data2(i,1); %- 102; % shift the data to start on this
                                           % (arbitrary!) date
    if (~isnan(lib_data2(i,2)) & ~isnan(lib_data2(i,3)))
        if (lib_data2(i,1) >= 0)
            temp_matrix1(count1,:) = lib_data2(i,:);
            count1=count1+1;
        end
    end
    if (~isnan(guin_data2(i,2)) & ~isnan(guin_data2(i,3)))
        if (guin_data2(i,1) >= 0)
            temp_matrix2(count2,:) = guin_data2(i,:);
            count2=count2+1;
        end
    end
    if (~isnan(sl_data2(i,2)) & ~isnan(sl_data2(i,3)))
        if (sl_data2(i,1) >= 0)
            temp_matrix3(count3,:) = sl_data2(i,:);
            count3=count3+1;
        end
    end
end

lib_data=flipud(temp_matrix1(find(temp_matrix1(:,1),1,'first'):find(temp_matrix1(:,1),1,'last')+1,:));
guin_data=flipud(temp_matrix2(find(temp_matrix2(:,1),1,'first'):find(temp_matrix2(:,1),1,'last')+1,:));
sl_data=flipud(temp_matrix3(find(temp_matrix3(:,1),1,'first'):find(temp_matrix3(:,1),1,'last')+1,:));

%size(lib_data)
%size(guin_data)
%size(sl_data)

%Adjusting from cumulative case counts to current
lib_data(:,2) = lib_data(:,2) - lib_data(:,3);
guin_data(:,2) = guin_data(:,2) - guin_data(:,3);
sl_data(:,2) = sl_data(:,2) - sl_data(:,3);

lib_data = lib_data(25:57,:);
guin_data_a = guin_data([2:72],:);
guin_data_b = guin_data([74:end],:);
guin_data = [guin_data_a; guin_data_b];
sl_data = sl_data(21:66,:);

tSpan = unique(sort([lib_data(:,1); guin_data(:,1); sl_data(:,1)]));
tSpan1 = lib_data(:,1);
tSpan2 = guin_data(:,1);
tSpan3 = sl_data(:,1);

first_day_all = max([tSpan1(1), tSpan2(1), tSpan3(1)]); % it's in liberia
temp = tSpan;
tSpan = [];
for i=1:length(temp)
    if (temp(i) >= first_day_all)
        tSpan = [tSpan; temp(i)];
    end
end
tSpan;

length(tSpan);
%length(tSpan1)

%lib_data(1,:)
%guin_data
%sl_data(1,:)

% parameter initialization
% order: [liberia, guinea, sl]
alpha  = [0; 0; 0];  % population growth constant (known empirically) 
pop1 = 4.294e6;
pop2 = 11.76e6;
pop3 = 6.092e6;
beta1 = [.718/pop1; .316/pop2; .803/pop3];
beta2 = [1.34/pop1; .446/pop2; 5.35/pop3];
beta3 = [1.05e-5/pop1; .0325/pop2; 1.69e-4/pop3];
delta  = [1/9;      1/9;       1/9     ];  %1/incubation period for Ebola Virus
gamma1 = [1/10;      1/10;       1/10     ];
gamma2 = [.0524; .0269; .0271];
psi    = [.486, 1.11, .739];
rho1   = [.77; .649; .451];
rho2   = [.7; .59; .41];
omega  = [1;         1;          1        ]; %population death constant 
phi = IG;


first_day_guin = tSpan2(1);
first_day_sl = tSpan3(1);
first_day_all;
guin_ic_0 = [11.76e6, guin_data(1,2), guin_data(1,2), 0, guin_data(1,3), 0, 0];
sl_ic_0 = [6.092e6, sl_data(1,2), sl_data(1,2), 0, sl_data(1,3), 0, 0]; 
guin_ic = find_coupled_ic(first_day_guin, first_day_all, guin_ic_0, alpha(2), beta1(2), beta2(2), ...
          beta3(2), delta(2), gamma1(2), gamma2(2), psi(2), rho1(2), rho2(2), omega(2));
sl_ic = find_coupled_ic(first_day_sl, first_day_all, sl_ic_0, alpha(3), beta1(3), beta2(3), ...
          beta3(3), delta(3), gamma1(3), gamma2(3), psi(3), rho1(3), rho2(3), omega(3));

          %sl_ic = find_coupled_ic()



% using initial condition from Liberia on 7/2/14
y0 = [4.294e6, lib_data(1,2), lib_data(1,2), 0, lib_data(1,3), 0, 0, ...
      guin_ic, sl_ic];
%    11.75e6, guin_data(1,2), guin_data(1,2), 0, guin_data(1,3), 0, 0, ...
%    6.092e6, sl_data(1,2), sl_data(1,2), 0, sl_data(1,3), 0, 0]; 
%odejac = @(t,u,up) jac(u, alpha, beta1, beta2, beta3, delta, gamma1, gamma2, psi, rho1, rho2, omega); 
odefun =@(t,u) SEIHRRR(t, u, alpha, beta1, beta2, beta3, delta, gamma1, gamma2, ...
                       psi, rho1, rho2, omega, phi);
%opts = odeset('Jacobian', odejac);
[t,y] = ode15s(odefun, tSpan, y0);


clf;
figure;
plot(t,y(:,3)+y(:,4),'b-',t,y(:,5)+y(:,6),'r--',lib_data(:,1), ...
lib_data(:,2),'b.',lib_data(:,1), lib_data(:,3), 'rs','MarkerSize',10)
legend('I+H','R_I + R_B','Number of cases data','Deaths data', 'Location', 'Northwest')
axis([tSpan(1) (max(tSpan1)+10) 0 max(lib_data(:,2))*1.5])
title(sprintf('6-parameter Liberia Coupled Cases & Death Fit \n phi12 = %0.3g, phi13 = %0.3g, phi21 = %0.3g, phi23 = %0.3g, phi31= %0.3g, phi32= %0.3g',IG(1),IG(2),IG(3),IG(4),IG(5),IG(6)))
xlabel('Days')
ylabel('Population')
print('-painters', '-depsc', '/home/eric/ebola_new/liberia.eps')
saveas(gcf,'/home/eric/ebola_new/liberia.fig','fig')

clf;
figure;
plot(t,y(:,10)+y(:,11),'b-',t,y(:,12)+y(:,13),'r--',guin_data(:,1), ...
guin_data(:,2),'b.',guin_data(:,1), guin_data(:,3), 'rs','MarkerSize',10)
legend('I+H','R_I + R_B','Number of cases data','Deaths data', 'Location', 'Northwest')
axis([tSpan(1) (max(tSpan2)+10) 0 max(guin_data(:,2))*1.5])
title(sprintf('6-parameter Guinea Coupled Cases & Death Fit \n phi12 = %0.3g, phi13 = %0.3g, phi21 = %0.3g, phi23 = %0.3g, phi31= %0.3g, phi32= %0.3g',IG(1),IG(2),IG(3),IG(4),IG(5),IG(6)))
xlabel('Days')
ylabel('Population')
print('-painters', '-depsc', '/home/eric/ebola_new/guinea.eps')
saveas(gcf,'/home/eric/ebola_new/guinea.fig','fig')

clf;
figure;
plot(t,y(:,17)+y(:,18),'b-',t,y(:,19)+y(:,20),'r--',sl_data(:,1), ...
sl_data(:,2),'b.',sl_data(:,1), sl_data(:,3), 'rs','MarkerSize',10)
legend('I+H','R_I + R_B','Number of cases data','Deaths data', 'Location', 'Northwest')
axis([tSpan(1) (max(tSpan3)+10) 0 max(sl_data(:,2))*1.5])
title(sprintf('6-parameter Sierra Leone Coupled Cases & Death Fit \n phi12 = %0.3g, phi13 = %0.3g, phi21 = %0.3g, phi23 = %0.3g, phi31= %0.3g, phi32= %0.3g',IG(1),IG(2),IG(3),IG(4),IG(5),IG(6)))
xlabel('Days')
ylabel('Population')
print('-painters', '-depsc', '/home/eric/ebola_new/sl.eps')
saveas(gcf,'/home/eric/ebola_new/sl.fig','fig')



err = 0;
for i=1:length(lib_data)
    temp_t = lib_data(i,1);
    if (temp_t > first_day_all)
        index = find(t == temp_t);
        err = err + sqrt((lib_data(i, 3) - (y(index,5) + y(index, 6)))^2 + (lib_data(i,2) - ...
        (y(index,3) + y(index,4)))^2);
    end
end
for i=1:length(guin_data)
    temp_t = guin_data(i,1);
    if (temp_t > first_day_all)
        index = find(t == temp_t);
        err = err + sqrt((guin_data(i, 3) - (y(index,12) + y(index, 13))).^2 + (guin_data(i,2) - ...
        (y(index,10) + y(index,11)))^2);
    end
end
for i=1:length(sl_data)
    temp_t = sl_data(i,1);
    if (temp_t > first_day_all)
        index = find(t == temp_t);
        err = err + sqrt((sl_data(i, 3) - (y(index,19) + y(index, 20)))^2 + (sl_data(i,2) - ...
        (y(index,17) + y(index,18)))^2);
    end
end

if (sum(phi < 0) > 0)
    err = inf;
end

