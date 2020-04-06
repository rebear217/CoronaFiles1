function fit = fitToTimeseries(data,extendTimeLength)

    if nargin < 2
        extendTimeLength = 7;
    end

    myRobustOptions = {[],'andrews','bisquare','cauchy','fair',...
        'huber','logistic','talwar','welsch'};    
    RobustOption = 1;
    
    totalTime = length(data);
    times = 1:(totalTime);

    baseFun0 = @(t,p) abs(p(3))./(1+exp(p(1))*exp(-abs(p(2))*t));
    baseFun1 = @(t,p) 1./(1+exp(p(1))*exp(-abs(p(2))*t));
    baseFun2 = @(t,p) 1./(1+exp(p(1))*exp(-abs(p(2))*(t-abs(p(3)))));
    
    fullFun = @(t,p) p(7)*(exp(-abs(p(6)))*baseFun1(t,p(1:2)) + ...
        (1-exp(-abs(p(6))))*baseFun2(t,p(3:5)));

	doubleDelayFun = @(t,p) p(8)*(exp(-abs(p(7)))*baseFun1(t,p(1:3)) + ...
        (1-exp(-abs(p(7))))*baseFun2(t,p(4:6)));

    %%
    
    deathsSoFar = data(end);
    
    fitdata = data;
    fittimes = times;

    p1guess = (1/deathsSoFar + 1);
    p3guess = totalTime/2;
    p2guess = 2*log(p1guess);

    %guess = [p1guess , p2guess];
    guess = [0 , 1];

    opts = statset('MaxIter',1000,'TolX',1e-10,'TolFun',1e-10);
    
    
    try
    	fit0 = fitnlm(fittimes,fitdata/deathsSoFar,@(b,x)baseFun1(x,b),guess);
    catch
        J = length(fittimes) - 1;
        worked = false;
        while J <= length(fittimes)
            try
                fit0 = fitnlm(fittimes(1:J),fitdata(1:J)/fitdata(J),@(b,x)baseFun1(x,b),guess);
                guess = fit0.Coefficients.Estimate;
                worked = true;
                J = J + 1;
            catch
                J = J + 1;
            end
        end        
    end
    
    guess0 = fit0.Coefficients.Estimate;
    
    guess = [guess0 ; deathsSoFar];

    %NLM = NonLinearModel.fit(fittimes,fitdata,@(b,x)baseFun(x,b),guess);
    NLM = fitnlm(fittimes,fitdata,@(b,x)baseFun0(x,b),guess);
    betterGuess = NLM.Coefficients.Estimate;
    SE = NLM.Coefficients.SE;
    pValues = NLM.Coefficients.pValue;
    fitFunction = @(t)NLM.feval(t);
    AIC = NLM.ModelCriterion.AIC;

    fullGuess = [betterGuess(1) betterGuess(2) 1 1 totalTime/2 0.1 deathsSoFar];

    %fullNLM = NonLinearModel.fit(fittimes,fitdata,@(b,x)fullFun(x,b),fullGuess);    
    
    
	opts = statset('TolFun',1e-10,'MaxIter',1000); 
    opts.RobustWgtFun = myRobustOptions{RobustOption};    
    try
        fullNLM = fitnlm(fittimes,fitdata,@(b,x)fullFun(x,b),fullGuess,'Options',opts);
        fullGuess = fullNLM.Coefficients.Estimate;
    catch
        %the above method can fail, so we can define a new grid to solve on
        %of reduced length, by J points and try again:
        J = length(fittimes) - 4;
        while J < length(fittimes)
            try
                tempNLM = fitnlm(fittimes(1:J),fitdata(1:J),@(b,x)fullFun(x,b),fullGuess);
                fullNLM = fitnlm(fittimes(1:J+1),fitdata(1:J+1),@(b,x)fullFun(x,b),tempNLM.Coefficients.Estimate);
                fullGuess = fullNLM.Coefficients.Estimate;
                %should this update 'worked'?
            end
            J = J + 1;
        end
    end
    
    ddGuess = [fullGuess(1:2) ; fullGuess(6) ; fullGuess(3:7)];
    ddNLM = fitnlm(fittimes,fitdata,@(b,x)doubleDelayFun(x,b),ddGuess,'Options',opts);
    
    fulllbetterGuess = fullNLM.Coefficients.Estimate;
    fullSE = fullNLM.Coefficients.SE;
    fullpValues = fullNLM.Coefficients.pValue;
    fullFitFunction = @(t)fullNLM.feval(t);
    fullAIC = fullNLM.ModelCriterion.AIC;    
    
    fit.NLM = NLM;
    fit.fitFunction = fitFunction;
    fit.solution = betterGuess;
    fit.SE = SE;
    fit.pValues = pValues;
    fit.AIC = AIC;
    fit.adjR2 = NLM.Rsquared.Adjusted;
    
    fit.fullNLM = fullNLM;
    fit.fullFitFunction = fullFitFunction;
    fit.fullsolution = fulllbetterGuess;
    fit.fullSE = fullSE;
    fit.fullpValues = fullpValues;
    fit.fulLAIC = fullAIC;
    fit.fulladjR2 = fullNLM.Rsquared.Adjusted;

    ddbetterGuess = ddNLM.Coefficients.Estimate;
    ddSE = ddNLM.Coefficients.SE;
    ddpValues = ddNLM.Coefficients.pValue;
    ddFitFunction = @(t)ddNLM.feval(t);
    ddAIC = ddNLM.ModelCriterion.AIC;    

    fit.ddNLM = ddNLM;
    fit.ddFitFunction = ddFitFunction;
    fit.ddsolution = ddbetterGuess;
    fit.ddSE = ddSE;
    fit.ddpValues = ddpValues;
    fit.ddAIC = ddAIC;
    fit.ddadjR2 = ddNLM.Rsquared.Adjusted;    
    
    fit.times = times;
    
    %now do error bars out of the way of the above code for clarity:
    [beta,resid,J,sigma] = nlinfit(fittimes,fitdata,@(b,x)baseFun0(x,b),betterGuess);
    extendTime = [fittimes fittimes(end)+(1:extendTimeLength)];
    [deltaFit, delta] = nlpredci(@(b,x)baseFun0(x,b),extendTime,beta,resid,'Covar',sigma);

    fit.upperCI = deltaFit + delta;
    fit.lowerCI = deltaFit - delta;
    fit.midCI = deltaFit;
    fit.extendTime = extendTime;
    
    %policy model error bars:    
    [beta,resid,J,sigma] = nlinfit(fittimes,fitdata,@(b,x)fullFun(x,b),fulllbetterGuess);
    [deltaFit, delta] = nlpredci(@(b,x)fullFun(x,b),extendTime,beta,resid,'Covar',sigma);

    fit.fullupperCI = deltaFit + delta;
    fit.fulllowerCI = deltaFit - delta;
    fit.fullmidCI = deltaFit;
    
    gprMdl = fitrgp(fittimes',fitdata','Basis','linear','FitMethod','exact','PredictMethod','exact');
    deathGPprediction = resubPredict(gprMdl)';
    
    fit.GRP = deathGPprediction;
end
