function out = isolationSIRmodel(t,x,p,isoTime)

    lam = abs(p(1));
    %this limits the policy-on rate:
    ip = 5*exp(-abs(p(2)));
    
    %we are assuming de-isolation does not occur:
    %im = abs(p(3));
    im = 0;
    %if isoTime is NaN, the following becomes 0 and no-one can be placed in
    %isolation:
    ip = ip*(t >= isoTime);
    
    d = abs(p(3));
    
    %to ensure lamE < lam we do this:
    lamE = exp(-abs(p(4)))*lam;
    rho = abs(p(5));

    S = x(1);
    I = x(2);
    R = x(3);
    Si = x(4);
    Ii = x(5);
    Ri = x(6);
    D = x(7);
        
    Sdot = -lam*S*I -lamE*S*Ii - ip*S + im*Si;
    Idot = lam*S*I + lamE*S*Ii - ip*I + im*Ii - (rho+d)*I;
    Rdot = rho*I - ip*R + im*Ri;

    Sidot = ip*S - im*Si - lamE*(I+Ii)*Si;
    Iidot = ip*I - im*Ii + lamE*(I+Ii)*Si - (rho+d)*Ii;
    Ridot = ip*R - im*Ri + rho*Ii;
    
    Ddot = d*(Ii + I);
    
    out = [Sdot,Idot,Rdot,Sidot,Iidot,Ridot,Ddot]';

end