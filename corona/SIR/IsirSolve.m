function solution = IsirSolve(p,T,isoTime,initialIsolationGuess)

    S0 = initialIsolationGuess.S0;
    I0 = initialIsolationGuess.I0;
    %delta = initialIsolationGuess.delta;
    
    IC = [S0, I0, 0, 0, 0, 0, 0]';
    model = @(t,x)isolationSIRmodel(t,x,p,isoTime);
    options = odeset('NonNegative',[1 1 1 1 1 1 1]');
    
    [soltimes,output] = ode113(model,[1,T(end)],IC,options);
    D = output(:,7);
    solution = interp1(soltimes,D,T,'pchip');
    
end