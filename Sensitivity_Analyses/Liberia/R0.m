function [ r0 ] = R0( xx , a , b )
%R0 takes as input values xx uniformly sampled from the interval [-1,1]^8
%scales these random samples to be within the defined intervals [a,b] then
%computes the reproductive ratio

b1 = (xx(1)+1)*(b(1)-a(1))/2 + a(1); 
b2 = (xx(2)+1)*(b(2)-a(2))/2 + a(2);
b3 = (xx(3)+1)*(b(3)-a(3))/2 + a(3);
rho1=(xx(4)+1)*(b(4)-a(4))/2 + a(4);
g1 = (xx(5)+1)*(b(5)-a(5))/2 + a(5);
g2 = (xx(6)+1)*(b(6)-a(6))/2 + a(6);
o =  (xx(7)+1)*(b(7)-a(7))/2 + a(7);
psi= (xx(8)+1)*(b(8)-a(8))/2 + a(8);


r0 =  (b1 + (b2*rho1*g1)*(o)^(-1) + b3*psi*(g2)^(-1) )*(g1 + psi)^(-1);

end

