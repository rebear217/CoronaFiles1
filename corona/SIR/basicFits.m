function gprMdl = basicFits(MATdata,savePDF)

    if nargin == 0 || isempty(MATdata)
        load('../data/deathData.mat')
    end
    if nargin < 2
        savePDF = false;
    end

    clc
    close all
    
    parameters = defaulParameters();
    
    extendDays = 2;
    
    countryStrings = {'China', ...
        'Denmark', ...
        'France', ...
        'Germany', ...
        'Italy', ...
        'Switzerland', ...
        'Spain', ...
        'United Kingdom', ...
        'US'};

    turningPointCountries = [1 3 5 6 7];
    PL = zeros(1,length(turningPointCountries)+1);
    
	fitFunctions{1} = @(p,t)exp(p(1) + p(2)*t);
    guesses{1} = [1 0];
	dfitFunctions{1} = @(p,t)fitFunctions{1}(p,t)*p(2);

	fitFunctions{2} = @(p,t)exp(p(1) + p(2)*t + p(3)*t.^2);
	dfitFunctions{2} = @(p,t)fitFunctions{2}(p,t).*(p(2) + 2*p(3)*t);
    guesses{2} = [1 0 0];

	fitFunctions{3} = @(p,t)(p(1) + p(2)*t + p(3)*t.^2).^2;
	dfitFunctions{3} = @(p,t)2*fitFunctions{3}(p,t).*(p(2) + 2*p(3)*t);
    guesses{3} = [1 1 1];
    
	fitFunctions{4} = @(p,t)(p(1) + p(2)*t + ...
        p(3)*t.^2 + p(4).*t.^3).^2;
	dfitFunctions{4} = @(p,t)2*fitFunctions{4}(p,t).*(p(2) + 2*p(3)*t + 3*p(4).*t.^2);
    guesses{4} = [1 1 1 1];
    
    fitFunctions{5} = @(p,t)(p(1) + p(2)*t + ...
        p(3)*t.^2 + p(4).*t.^3 + p(5).*t.^4).^2;
	dfitFunctions{5} = @(p,t)2*fitFunctions{5}(p,t).*(p(2) + 2*p(3)*t + 3*p(4).*t.^2 + 4*p(5).*t.^3);
    guesses{5} = [1 1 1 1 1];
    
    figure(1)
	set(1,'pos',[38 1 969 704]);
    figure(2)
	set(2,'pos',[272   147   846   474]);

    G = length(guesses);
    qj = 1;
	for ctry = 1:length(countryStrings)
        figure(1)

        countryStr = countryStrings{ctry};

        cj = find(strcmp(countryStr,MATdata.country));

        cdata = MATdata.deathData{cj};
        scdata = sum(cdata,1);
        deathdayData = diff(scdata);

        firstDeath = find(deathdayData,1,'first');
        dataToFit = scdata(firstDeath:end);
        L = length(dataToFit);        
        fittimes = 1:L;

        daydeathdataToFit = deathdayData(firstDeath:end);

        FT = fittimes(2:end);
        extendTime = FT(1):0.5:FT(end)+extendDays;

        AICs = zeros(1,G);
        
        gprMdl = fitrgp(FT',daydeathdataToFit','Basis','linear',...
            'FitMethod','exact','PredictMethod','exact');
        deathGPprediction = resubPredict(gprMdl)';
        
        for fj = 1:G
            
            fitFun = fitFunctions{fj};
            guess = guesses{fj};
            fit = fitnlm(FT,daydeathdataToFit,fitFun,guess);
            SE = fit.Coefficients.SE;
            pValues = fit.Coefficients.pValue;
            fitFunction = @(t)fit.feval(t);
            AIC = fit.ModelCriterion.AIC;
            adjR2 = fit.Rsquared.Adjusted;
            betterGuess{fj} = fit.Coefficients.Estimate;
            
            AICs(fj) = AIC;
            
            subplot(3,3,ctry)
            pl{fj} = plot(fittimes(2:end),daydeathdataToFit,'.k','markersize',24);
            hold on
            s = (fj - 1)/(G-1);
            col = [s 0 1-s];
            pl{fj+G} = plot(extendTime,fitFunction(extendTime),'--k','color',col,'linewidth',1);
            axis tight
            
        end
        
        [Amax,I] = min(AICs);
        fitderiv = dfitFunctions{I}(betterGuess{I},extendTime);
        disp([countryStr,': best fit function ',num2str(I)]);
        set(pl{I},'linewidth',2,'color','k','linestyle','-');
        
        for fj = 1:2*G
            delete(pl{fj});
        end
        
        fitFun = fitFunctions{I};
        [beta,resid,J,sigma] = nlinfit(FT,daydeathdataToFit,fitFun,betterGuess{I});
        [deltaFit, delta] = nlpredci(fitFun,extendTime,beta,resid,'Covar',sigma);
        
        upperCI = deltaFit + delta;
        lowerCI = deltaFit - delta;
        
        plotshaded(extendTime, [upperCI ; lowerCI],parameters.grey);
        plot(extendTime,fitFunction(extendTime),'-','color',[0 0 0],'linewidth',2);
        plot(FT,deathGPprediction,'-r','linewidth',1);
        plot(fittimes(2:end),daydeathdataToFit,'.k','markersize',24)

        legend({'95% CIs','regression','GPR',countryStr},'Location','northwest')
        legend('boxoff')

        axis tight
        YL = ylim;
        ylim([0 YL(2)])
        %xlim([1 extendTime(end)])
        
        plot([extendTime(end-extendDays) extendTime(end-extendDays)],[0 YL(2)],'-k','linewidth',1)
        
        if ismember(ctry,[7 8 9])
            xlabel('day')
        end
        if ismember(ctry,[1 4 7])
            ylabel('deaths each day')
        end

        figure(2)
        s = (ctry - 1) / (length(countryStrings) - 1);
        col = [s 0.5 1-s];
        sfitderiv = 2*(-0.5 + (fitderiv - min(fitderiv)) ./ (max(fitderiv) - min(fitderiv)));
        if ismember(ctry,turningPointCountries)
            PL(qj) = plot(extendTime,sfitderiv,'-r','linewidth',2,'color',col);
            hold on
            extP = plot(extendTime(end-extendDays+1:end),sfitderiv(end-extendDays+1:end),...
                'o','color',col);
            sc = sfitderiv(1:end-1).*sfitderiv(2:end);
            F = find(sc < 0,1,'last');
            ZT = (extendTime(F)+extendTime(F+1))/2;
            if sfitderiv(F-1) > 0
                PL(length(turningPointCountries)+1) = plot(ZT,0,'.k','markersize',40);
            end
            qj = qj + 1;
        end
        
    end

    figure(2)
    axis tight
    XL = xlim;
    plot([2 XL(2)],[0 0],'-k');
    grid on
        
    legend([PL extP],{countryStrings{turningPointCountries},...
        'days of maximal death','circles: future extrapolations'},'Location','northeast');
    legend('boxoff');
    ylabel('normalised daily death derivatives');
    xlabel('days')
    
    if savePDF
        my_export_fig('derivativePlotOfDeathData.pdf')

        figure(1)
        my_export_fig('basicRegressionsToDeathData.pdf')
    end
    
end
