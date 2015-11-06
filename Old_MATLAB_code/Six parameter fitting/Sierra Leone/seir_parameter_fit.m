function err = seir_parameter_fit(IG)
%SIERRALEONE

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

%Adjusting from cumulative case counts to current
%lib_data(:,2) = lib_data(:,2) - lib_data(:,3);
%lib_data(:,2) = [lib_data(1,2); diff(lib_data(:,2))];
    
lib_data = lib_data(2:62,:);
tSpan = lib_data(:,1);

% parameter initialization
beta1  = IG(1);
beta2  = IG(2);
beta3  = IG(3);
delta  = 1/9;%1/21;%1/incubation period for Ebola Virus
gamma2 = IG(4);
psi    = IG(5);

gamma1 = 0.055;%0.0585; %0.068;
rho2 = 0.7985;%0.7084; %0.884;
rho1 = 0.7097;%0.7265; %0.853;

omega  = IG(6);  %population death constant 


y0 = [6.092e6, 0, lib_data(1,2), 0, 0, lib_data(1,3), 0];


odejac = @(t,u,up) jac(u, beta1, beta2, beta3, delta, gamma1, gamma2, psi, rho1, rho2, omega); 
odefun =@(t,u) SEIHRRR(t, u, beta1, beta2, beta3, delta, gamma1, gamma2, psi, rho1, rho2, omega);
opts = odeset('Jacobian', odejac);
[t,y] = ode15s(odefun, tSpan, y0, opts);

intI(1) = y(1,3);
for k = 2:length(y(:,3))
    intI(k) = trapz(t(1:k),y(1:k,3));
end
intI = intI';


 if beta1 < 0 || beta2 < 0 || beta3 < 0 || gamma2 < 0 || psi < 0 || omega < 1e-1 || gamma2 > (1/4.4) + (1/12.4) || beta2 > 2 || omega > 0.5
    err = inf;
 else
    err = sum(sqrt( (lib_data(:,3) - y(:,6) ).^2 + (lib_data(:,2) - intI).^2));
 end