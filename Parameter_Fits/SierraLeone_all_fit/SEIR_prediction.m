%SIERRA LEONE

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
tSpan = lib_data(:,1);

origspan = length(tSpan);
future = (5:5:30)';
tSpan = [tSpan; tSpan(end) + future];

% parameter initialization
Pop = 6.092e6;
delta = 1/9;

beta1 = 0.251/Pop;
beta2 = 0.395/Pop;
beta3 = 0.0791/Pop;
gamma1 = 0.051;
gamma2 = 0.0833;
rho1 = 0.766;
rho2 = 0.74; %0.885
psi = 0.442;
omega = 0.37;


R0 = (1/(gamma1+psi))*(beta1*Pop + (beta3*Pop*psi/gamma2) + (beta2*Pop*rho1*gamma1/omega))

% using initial condition from Liberia on 7/2/14
y0 = [Pop, 0, lib_data(1,2), 0, 0, lib_data(1,3), 0 ];


odejac = @(t,u,up) jac(u, beta1, beta2, beta3, delta, gamma1, gamma2, psi, rho1, rho2, omega); 
odefun =@(t,u) SEIHRRR(t, u, beta1, beta2, beta3, delta, gamma1, gamma2, psi, rho1, rho2, omega);
opts = odeset('Jacobian', odejac);
[t,y] = ode23s(odefun, tSpan, y0, opts);

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
