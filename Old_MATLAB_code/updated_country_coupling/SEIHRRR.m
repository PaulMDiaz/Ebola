function f=SEIHRRR(t, u, alpha1, beta1, beta2, beta3, delta, gamma1, gamma2, ...
psi, rho1, rho2, omega, phi)
S = [u(1), u(8), u(15)]; 
E= [u(2),u(9),u(16)]; 
I= [u(3),u(10),u(17)]; 
H= [u(4),u(11),u(18)]; 
RI= [u(5),u(12),u(19)]; 
RB=[u(6),u(13),u(20)]; 
RR=[u(7),u(14),u(21)];
%phi is [phi12 13 21 23 31 32]

f =[ 
    (alpha1(1)*S(1) - beta1(1)*S(1)*I(1) - beta2(1)*S(1)*RI(1) - beta3(1)*S(1)*H(1)) ; 
    
     (beta1(1)*S(1)*I(1)-delta(1)*E(1) + beta2(1)*S(1)*RI(1) + beta3(1)*S(1)*H(1)) ... 
     -  phi(1)*E(1) - phi(2)*E(1) + phi(3)*E(2) + phi(5)*E(3);
     
     (delta(1)*E(1)-gamma1(1)*I(1) - psi(1)*I(1) );
     
     psi(1)*I(1)-gamma2(1)*H(1) ;  
     
     rho1(1)*gamma1(1)*I(1) - omega(1)*RI(1) ; 
     
     omega(1)*RI(1)+rho2(1)*gamma2(1)*H(1) ; 
     
     (1-rho1(1))*gamma1(1)*I(1) + (1-rho2(1))*gamma2(1)*H(1) ;
     alpha1(2)*S(2) - beta1(2)*S(2)*I(2)-beta2(2)*S(2)*RI(2) - beta3(2)*S(2)*H(2) ;

    (beta1(2)*S(2)*I(2)-delta(2)*E(2) + beta2(2)*S(2)*RI(2) + beta3(2)*S(2)*H(2)) ... 
     -  phi(3)*E(2) - phi(4)*E(2) + phi(1)*E(1) + phi(6)*E(3);
   
    (delta(2)*E(2)-gamma1(2)*I(2) - psi(2)*I(2) );
   
    psi(2)*I(2)-gamma2(2)*H(2) ;  
   
    rho1(2)*gamma1(2)*I(2) - omega(2)*RI(2) ; 
    
    omega(2)*RI(2)+rho2(2)*gamma2(2)*H(2) ; 
    
    (1-rho1(2))*gamma1(2)*I(2) + (1-rho2(2))*gamma2(2)*H(2) ;

     alpha1(3)*S(3) - beta1(3)*S(3)*I(3) - beta2(3)*S(3)*RI(3) - beta3(3)*S(3)*H(3) ; 
    
    (beta1(3)*S(3)*I(3)-delta(3)*E(3) + beta2(3)*S(3)*RI(3) +beta3(3)*S(3)*H(3)) ...
     -  phi(5)*E(3) - phi(6)*E(3) + phi(2)*E(1) + phi(4)*E(2);
    
    (delta(3)*E(3)-gamma1(3)*I(3) - psi(3)*I(3) );
    
    psi(3)*I(3)-gamma2(3)*H(3) ;  
    
    rho1(3)*gamma1(3)*I(3) - omega(3)*RI(3) ; 
    
    omega(3)*RI(3)+rho2(3)*gamma2(3)*H(3) ; 
    
    (1-rho1(3))*gamma1(3)*I(3) + (1-rho2(3))*gamma2(3)*H(3)
                                               ] ;

  
