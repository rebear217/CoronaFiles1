function dailyDeathsAnalysis(times,dataToFit,pg,infit)
    
    figure(111)
    plot(times(2:end),diff(infit.feval(times)));
    hold on
    plot(times(2:end),diff(dataToFit));
    
    dDeaths = diff(infit.feval(times));
    
    newtimes = times(2:end);
    dailyDeaths = diff(dataToFit);
    
    opts = statset('MaxIter',2000,'TolX',1e-30,'TolFun',1e-30,'Display','iter');
    fit = fitnlm(newtimes,dailyDeaths,@(p,T)dailySolve(p,T,dDeaths(1)),pg,'Options',opts)
    betterGuess = abs(fit.Coefficients.Estimate);

    plot(newtimes,fit.feval(newtimes));
     
end