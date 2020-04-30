function solution = dailySolve(p,T,IC)

    %IC = 0;
    model = @(t,x)dailyF(x,p);
    %options = odeset('NonNegative',[1]);
    [soltimes,output] = ode113(model,[T(1),T(end)],IC);
    solution = interp1(soltimes,output,T,'pchip');
    
end