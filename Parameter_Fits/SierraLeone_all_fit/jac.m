function dF=jac(u, beta1, beta2, beta3, delta, gamma1, gamma2, psi, rho1, rho2, omega)
S = u(1); E=u(2); I=u(3); H=u(4); RI=u(5); RB=u(6); RR = u(7);


dF = [     -beta1*I-beta2*RI-beta3*H,            0,      -beta1*S,         -beta3*S,  -beta2*S, 0, 0;
            beta1*I+beta2*RI+beta3*H,       -delta,       beta1*S,          beta3*S,   beta2*S, 0, 0;
                                   0,        delta,   -gamma1-psi,                0,         0, 0, 0;
                                   0,            0,           psi,          -gamma2,         0, 0, 0;
                                   0,            0,   rho1*gamma1,                0,    -omega, 0, 0;
                                   0,            0,             0,      rho2*gamma2,     omega, 0, 0;
                                   0,            0, (1-rho1)*gamma1,(1-rho2)*gamma2,         0, 0, 0];
                               
                               