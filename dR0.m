function [ dr0 ] = dR0( xx )
%DR0 This function takes as input xx in [-1,1]^8 and returns a scaled value
%of the gradient function for R0 
%   Detailed explanation goes here


% Libera
pop = 4.294e6;
b1 = (xx(1)+1)*(0.4-0.1)/2 + 0.1; 
b1 = b1/pop;
b2 = (xx(2)+1)*(0.4-0.1)/2 + 0.1;
b2 = b2/pop;
b3 = (xx(3)+1)*(0.2-0.05)/2 + 0.05;
b3 = b3/pop;
rho1 = (xx(4)+1)*(1-0.5)/2 + 0.5;
g1 = (xx(5)+1)*(0.1702-0.0276)/2 + 0.0276;
g2 = (xx(6)+1)*(.21 -0.081)/2 + 0.081;
o = (xx(7)+1)*(0.5-0.25)/2 + 0.25;
%psi = (xx(8)+1)*(0.7-0.3)/2+0.3;
psi = (xx(8)+1)*(1-0.04)/2+0.3;


dr0 = [
%beta1 
g1+psi).^(-1);
%beta2
1.*o.^(-1).*(g1+psi).^(-1).*rho1;
%beta3
2.^(-1).*psi.*(g1+psi).^(-1);
%rho1
b2.*g1.*o.^(-1).*(g1+psi).^(-1);
%gamma1
b2.*o.^(-1).*(g1+psi).^(-1).*rho1+(-1).*(g1+psi).^(-2).*(b1+b3.*g2.^(-1).*psi+b2.*g1.*o.^(-1).*rho1);
%gamma2
(-1).*b3.*g2.^(-2).*psi.*(g1+psi).^(-1);
%omega
(-1).*b2.*g1.*o.^(-2).*(g1+psi).^(-1).*rho1;
%psi
b3.*g2.^(-1).*(g1+psi).^(-1)+(-1).*(g1+psi).^(-2).*(b1+b3.*g2.^(-1).*psi+b2.*g1.*o.^(-1).*rho1);

];

dr0 = dr0.*[0.3/pop;0.3/pop;0.15/pop;0.5;(0.1702-0.0276);(.21 -0.081);0.25;1-0.04]*0.5;  


end

