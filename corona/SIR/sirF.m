function Ddot = sirF(x,p)

    A = abs(p(1));
    B = abs(p(2));
    C = abs(p(3));

    D = x(1);    
    Ddot = A - B*exp(-abs(D)) - C*D;

end

