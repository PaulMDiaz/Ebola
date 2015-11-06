function err = seir_parameter_fit(IG)
%LIBERIA

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

%Adjusting total cases to non-deaths only
%lib_data(:,2) = lib_data(:,2) - lib_data(:,3);

%Adjusting from cumulative case counts to current
% for i = 2:length(lib_data(:,2))
%     temp(i-1) = max(0,lib_data(i,2) - lib_data(i-1,2));
% end
% lib_data(2:end,2) = temp;

lib_data = lib_data(2:57,:); %25:57
tSpan = lib_data(:,1);

% parameter initialization
beta1  = IG(1);
beta2  = IG(2);
beta3  = IG(3);
delta  = 1/9;%1/21;%1/incubation period for Ebola Virus
gamma2 = IG(4);
psi    = IG(5);

gamma1 = 0.0585;%0.0585; %0.068;
rho2 = 0.8523;%0.7084; %0.884;
rho1 = 0.7268;%0.7265; %0.853;

omega  = IG(6);  %population death constant 

% using initial condition from Liberia on 7/2/14
y0 = [4.294e6, 0, lib_data(1,2), 0, 0, lib_data(1,3), 0 ];


odejac = @(t,u,up) jac(u, beta1, beta2, beta3, delta, gamma1, gamma2, psi, rho1, rho2, omega); 
odefun =@(t,u) SEIHRRR(t, u, beta1, beta2, beta3, delta, gamma1, gamma2, psi, rho1, rho2, omega);
opts = odeset('Jacobian', odejac);
[t,y] = ode15s(odefun, tSpan, y0, opts);

intI(1) = y(1,3);
for k = 2:length(y(:,2))
    intI(k) = trapz(t(1:k),y(1:k,3));
end
intI = intI';

 if beta1 < 0 || beta2 < 0 || beta3 < 0 || gamma2 < 0 || psi < 0 || omega < 1e-1 || gamma2 > (1/4.5) + (1/12.8) || beta2 > 2 || omega > 1
    err = inf;
 else
    err = sum(sqrt( 1.5*(lib_data(:,3) - (y(:,6)) ).^2 +  (lib_data(:,2) - intI).^2));
 end