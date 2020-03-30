function out = f(x,p)

    lam = p(1);
    ip = p(2);
    im = p(3);
    d = p(4);
    rho = p(5);
    lamE = p(6);
    rhoI = p(7);

    S = x(1);
    I = x(2);
    R = x(3);
    Si = x(4);
    Ii = x(5);
    Ri = x(6);
    D = x(7);
    
    Sdot = -lam*S*I -lamE*S*Ii - ip*S + im*Si;
    Idot = lam*S*I + lamE*S*Ii - ip*I + im*Ii - (d+rho)*I;
    Rdot = rho*I - ip*R + im*Ri;

    Sidot = ip*S - im*Si - lamE*(I+Ii)*Si;
    Iidot = ip*I - im*Ii + lamE*(I+Ii)*Si - (d+rhoI)*Ii;
    Ridot = ip*R - im*Ri + rhoI*Ii;
    
    Ddot = d*(Ii + I);
    
    out = [Sdot,Idot,Rdot,Sidot,Iidot,Ridot,Ddot]';

end