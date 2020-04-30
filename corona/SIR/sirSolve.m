function solution = sirSolve(p,T)

    IC = 0;
    model = @(t,x)sirF(x,p);
    options = odeset('NonNegative',[1]);
    [soltimes,output] = ode113(model,[1,T(end)],IC,options);
    solution = interp1(soltimes,output,T,'pchip');
    
end