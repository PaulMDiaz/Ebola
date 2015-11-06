%SIERRA LEONE

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

lib_data = lib_data(2:57,:);
tSpan = lib_data(:,1);

tSpan = [tSpan; tSpan(end) + (5:5:30)'];

% parameter initialization
Pop = 4.294e6;
alpha  = 0;%0.024/365; % population growth constant (known empircally) 
beta1  = 6.59/Pop; 
beta2  = 5.01/Pop; 
beta3  = 2.43e-4/Pop; 
delta  = 1/9;%1/21;%1/incubation period for Ebola Virus
gamma1 = 1/10;
gamma2 = 0.117;
psi    = 5.11;
rho2   = 0.41;
rho1   = 1.1*rho2;
omega  = 2.73;

R0 = (1/(gamma1+psi))*(beta1*Pop + (beta3*Pop*psi/gamma2) + (beta2*Pop*rho1*gamma1/omega))

init_inf = lib_data(1,2);% - lib_data(1,3);

% using initial condition from Liberia on 7/2/14
y0 = [4.294e6, 0, init_inf, 0, lib_data(1,3), 0, 0 ];


odejac = @(t,u,up) jac(u, alpha, beta1, beta2, beta3, delta, gamma1, gamma2, psi, rho1, rho2, omega); 
odefun =@(t,u) SEIHRRR(t, u, alpha, beta1, beta2, beta3, delta, gamma1, gamma2, psi, rho1, rho2, omega);
opts = odeset('Jacobian', odejac);
[t,y] = ode23s(odefun, tSpan, y0, opts);

figure;
plot(t, y(:,3))
