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
lib_data = lib_data(2:62,:);
tSpan = lib_data(:,1);

Pop = 6.092e6;
delta  = 1/9;
%Initial guesses for parameters
beta1 = 0.133/Pop;%0.2/Pop;
beta2 = 0.192/Pop;%0.2/Pop; 
beta3 = 0.156/Pop;%0.08/Pop;
gamma2 = 0.12359;%0.1747;%0.0513;
psi = 0.445;
omega = 0.315;
gamma1 = 0.055;%0.0875;%0.0585; %0.068;
rho1 = 0.7097;%0.7268; %0.853;
rho2 = 0.8177;%0.7084; %0.884;

% using initial condition from Guinea
y0 = [6.092e6, 0, lib_data(1,2), 0, 0, lib_data(1,3), 0];


odejac = @(t,u,up) jac(u, beta1, beta2, beta3, delta, gamma1, gamma2, psi, rho1, rho2, omega); 
odefun = @(t,u) SEIHRRR(t, u, beta1, beta2, beta3, delta, gamma1, gamma2, psi, rho1, rho2, omega);
opts = odeset('Jacobian', odejac);
[t,y] = ode15s(odefun, tSpan, y0, opts);

intI = zeros(61,1);
intI(1) = y(1,3);
for k = 2:length(y(:,3))
    intI(k) = trapz(t(1:k),y(1:k,3));
end

err = sum(sqrt( (lib_data(:,3) - (y(:,6)) ).^2 +  (lib_data(:,2) - intI).^2))

R0 = (1/(gamma1+psi))*(beta1*Pop + (beta3*Pop*psi/gamma2) + (beta2*Pop*rho1*gamma1/omega))

Final_Infected = y(end,3)

figure;
plot(t,intI,'b',t,y(:,6),'r--',lib_data(:,1), ...
lib_data(:,2),'b.',lib_data(:,1), lib_data(:,3), 'rs','MarkerSize',10)
legend('Cumulative I(t)','R_B(t)','Cumulative cases data','Deaths data','Location','Northwest')
axis([tSpan(1) max(tSpan) 0 max(lib_data(:,2))*1.25])

