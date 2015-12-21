function [ dr0 ] = dR0( xx , a , b )
%DR0 This function takes as input xx in [-1,1]^8 and returns a scaled value
%of the gradient function for R0 in vector format
%   Detailed explanation goes here



b1 = (xx(1)+1)*(b(1)-a(1))/2 + a(1); 
b2 = (xx(2)+1)*(b(2)-a(2))/2 + a(2);
b3 = (xx(3)+1)*(b(3)-a(3))/2 + a(3);
rho1=(xx(4)+1)*(b(4)-a(4))/2 + a(4);
g1 = (xx(5)+1)*(b(5)-a(5))/2 + a(5);
g2 = (xx(6)+1)*(b(6)-a(6))/2 + a(6);
o =  (xx(7)+1)*(b(7)-a(7))/2 + a(7);
psi= (xx(8)+1)*(b(8)-a(8))/2 + a(8);




dr0 = [
%beta1 
(g1+psi).^(-1);
%beta2
g1.*o.^(-1).*(g1+psi).^(-1).*rho1;
%beta3
g2.^(-1).*psi.*(g1+psi).^(-1);
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

%This step applies the chain rule to correctly scale gradient values.
c = b - a;
dr0 = dr0.*c';

end

