function rhs = SIRmodel(x,p)

    S = x(1);
    I = x(2);
    R = x(3);
    D = x(4);
    
    lambda = p(1);
    rho = p(2);
    d = p(3);
    
    dotS = -lambda*S*I;
    dotI = lambda*S*I - (d+rho)*I;
    dotR = rho*I;
    dotD = d*I;
    
    rhs = [dotS ; dotI ; dotR ; dotD];

end
