function f=SEIHRRR(t, u, alpha, beta1, beta2, beta3, delta, gamma1, gamma2, psi, rho1, rho2, omega)
S = u(1); E=u(2); I=u(3); H=u(4); RI=u(5); RB=u(6); RR = u(7);
f =[ 
    (alpha*S -beta1*S*I -beta2*S*RI -beta3*S*H) ; 
    
     (beta1*S*I-delta*E + beta2*S*RI +beta3*S*H)  ;
     
     (delta*E-gamma1*I -psi*I );
     
     psi*I-gamma2*H ;  
     
     rho1*gamma1*I - omega*RI ; 
     
     omega*RI+rho2*gamma2*H ; 
     
     (1-rho1)*gamma1*I + (1-rho2)*gamma2*H 
                                                ] ;
 
   