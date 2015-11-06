function IC = find_coupled_ic(first_day_single, first_day_all, given_ic, alpha, beta1, beta2, ...
beta3, delta, gamma1, gamma2, psi, rho1, rho2, omega)

y0 = given_ic;
odefun = @(t,u) SEIHRRR_single(t, u, alpha, beta1, beta2, beta3, delta, gamma1,...
gamma2, psi, rho1, rho2, omega);
tSpan = [first_day_single, first_day_all];

[t, y] = ode15s(odefun, tSpan, y0);
plot(t,y(:,2),'y-',t,y(:,3)+y(:,4),'r-.',t,y(:,5)+y(:,6),'k-',t,y(:,7),'b-')
%print('-painters', '-depsc', '/home/eric/ebola/Optimal Treatment_Ebola_Outbreak_Western_Africa/Submission_JMB/MATLAB_runs/updated_country_coupling/temp.eps')
IC = y(end,:);

