
% import .csv value
filename= 'matlab_data.csv';
%filename= '/Users/redpatrioteq2/Dropbox/PIC_Math_group/Optimal Treatment_Ebola_Outbreak_Western_Africa/ebola_modeling_484/ebola_data/matlab_data.csv';

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

lib_data = lib_data(22:57,:); %25:57
tSpan = lib_data(:,1);

Pop = 4.294e6;
delta  = 1/9;
%Initial guesses for parameters
beta1 = 0.376/Pop;
beta2 = 0.135/Pop;
beta3 = 0.163/Pop;
gamma2 = 0.174;
psi = 0.5;
omega = 0.325;
gamma1 = 0.0542;
rho1 = 0.98;
rho2 = 0.98;

% using initial condition from Liberia on 7/2/14
y0 = [4.294e6, 0, lib_data(1,2), 0, 0, lib_data(1,3), 0 ];


odejac = @(t,u,up) jac(u, beta1, beta2, beta3, delta, gamma1, gamma2, psi, rho1, rho2, omega); 
odefun =@(t,u) SEIHRRR(t, u, beta1, beta2, beta3, delta, gamma1, gamma2, psi, rho1, rho2, omega);
opts = odeset('Jacobian', odejac);
[t,y] = ode15s(odefun, tSpan, y0, opts);

intI = zeros(length(y(:,2)),1);
intI(1) = y(1,3);
for k = 2:length(y(:,2))
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

