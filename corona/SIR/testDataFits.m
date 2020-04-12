function testDataFits()

    % this script generate noisy synthetic death timeseries
    % data and compares the fit of 2 model epidemics to data
    % We are seeking whether the model fits can distinguish
    % between controlled and non-controlled epidemics and these
    % synthetic examples work with very high likelihood:

    clear variables
    close all
    clc

    figure(1)
    set(1,'pos',[0           1        1249         704])

    parameters = defaulParameters();

    for j = 1:4

        load(['./synthMatFiles/syntheticDeaths_',num2str(j),'.mat'])

        baseoutput = deathNumbers;
        basetimes = allEpitimes;

        sigma = 0.1;
        R = (rand(size(basetimes)));
        simulatedData = deathNumbers + sigma*R;

        %%

        %[basestate,baseoutput,basetimes] = moveForwardOne(parameters,[],[],parameters.totalTime);
        %deathFromNoPolicy = basestate(end,7);
        %sigma = 0.02;
        %R = cumsum(rand(size(basetimes)));
        %simulatedData = baseoutput(:,7) + sigma*R;

        %%

        totalTime = 35;
        times = 0:0.1:totalTime;

        baseFun0 = @(t,p) abs(p(3))./(1+exp(p(1))*exp(-abs(p(2))*t));
        baseFun1 = @(t,p) 1./(1+exp(p(1))*exp(-abs(p(2))*t));
        baseFun2 = @(t,p) 1./(1+exp(p(1))*exp(-abs(p(2))*(t-abs(p(3)))));

        fullFun = @(t,p) p(7)*(exp(-abs(p(6)))*baseFun1(t,p(1:2)) + ...
            (1-exp(-abs(p(6))))*baseFun2(t,p(3:5)));

        %%

        fitdata = simulatedData(1:4:end);
        fittimes = basetimes(1:4:end);

        guess = [0 , 1 , simulatedData(end)];

        NLM = NonLinearModel.fit(fittimes,fitdata,@(b,x)baseFun0(x,b),guess)
        betterGuess = NLM.Coefficients.Estimate;
        SE = NLM.Coefficients.SE;
        pValues = NLM.Coefficients.pValue;
        fitFunction = @(t)NLM.feval(t);
        adjR2 = NLM.Rsquared.Adjusted;
        AIC = NLM.ModelCriterion.AIC;

        fullGuess = [betterGuess(1) betterGuess(2) 1 1 totalTime/2 0.1 deathNumbers(end)];

        fullNLM = NonLinearModel.fit(fittimes,fitdata,@(b,x)fullFun(x,b),fullGuess)
        fulllbetterGuess = fullNLM.Coefficients.Estimate;
        fullSE = fullNLM.Coefficients.SE;
        fullpValues = fullNLM.Coefficients.pValue;
        fullFitFunction = @(t)fullNLM.feval(t);
        fulladjR2 = fullNLM.Rsquared.Adjusted;
        fullAIC = fullNLM.ModelCriterion.AIC;

        %%

        subplot(2,2,j)
        semilogx(fittimes,fitdata,'.k','markersize',20)
        set(gca,'Ytick',[0 1 2 3]);
        set(gca,'Yticklabel',{'0%','1%','2%','3%'});
        set(gca,'Xtick',[0 1 2 5 10 50 100 200]);
        hold on
        plot(basetimes,fitFunction(basetimes),'-k','linewidth',1)
        plot(basetimes,fullFitFunction(basetimes),'-k','linewidth',2)
        xlabel('time')
        ylabel('freq death (%S(0))')

        thresh = parameters.OFFcontrolThresholds(j);

        axis tight
        legend({'synthetic controlled epidemic data',...
            ['no iso epidemic datafit (adj R^2 \approx ',num2str(adjR2,3),')'],...
            ['iso policy epidemic datafit (adj R^2 \approx ',num2str(fulladjR2,3),', \theta=',num2str(thresh),')']},...
            'Location','south')
        %legend('boxoff')

        RL = exp(-abs(AIC - fullAIC)/2);
        YL = ylim;
        text(0.6,0.95*YL(2),['relative likelihood \approx ',num2str(RL,3)]);

    end

    %%

    my_export_fig('syntheticTestTest.pdf')
    
end
