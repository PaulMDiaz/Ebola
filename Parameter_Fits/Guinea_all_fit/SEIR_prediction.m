%GUINEA

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
tSpan = lib_data(:,1);

origspan = length(tSpan);
future = (5:5:30)';
tSpan = [tSpan; tSpan(end) + future];

% parameter initialization
Pop = 11.76e6;
beta1 = 0.315/Pop;
beta2 = 0.16/Pop;
beta3 = 0.0165/Pop;
gamma1 = 0.295;
gamma2 = 0.016;
rho1 = 0.98;
rho2 = 0.93;
psi = 0.5;
omega = 0.3;

R0 = (1/(gamma1+psi))*(beta1*Pop + (beta3*Pop*psi/gamma2) + (beta2*Pop*rho1*gamma1/omega))

% using initial condition from Guinea on 7/2/14
y0 = [11.76e6, 0, lib_data(1,2), 0, 0, lib_data(1,3), 0 ];


odejac = @(t,u,up) jac(u, beta1, beta2, beta3, delta, gamma1, gamma2, psi, rho1, rho2, omega); 
odefun =@(t,u) SEIHRRR(t, u, beta1, beta2, beta3, delta, gamma1, gamma2, psi, rho1, rho2, omega);
opts = odeset('Jacobian', odejac);
%[t,y] = ode23s(odefun, tSpan, y0, opts);
[t,y] = ode15s(odefun, tSpan, y0,opts);

%figure;
%plot(t, y(:,3))

newInf = y(origspan+1:end, 3);
newt = t(origspan+1:end);
intI = zeros(length(future),1);

for k = 2:length(future)
    intI(k) = trapz(newt(1:k),newInf(1:k));
end

figure;
plot(future, intI)

TotNewInf = intI(end)

DailyInf = y(end,3)
