function [outsideFilename,...
    UKestimate,UKestimateMax99,latexStr] = fitToCountryData(MATdata)

    warning('off')
    
    close all
    clc
    
    pickOutctry = 10;
	extendTimeLength = 14;
    
    parameters = defaulParameters();
	
    %'Austria',...

    countryStrings = parameters.countryStrings;
	outsideFilename = [countryStrings{pickOutctry},'additionalDatafit.pdf'];

    %%

    %for ctry = 1:2

    latexStr = ['-------------------' char(13)];
    latexStr = [latexStr 'latex format table:' char(13)];
    latexStr = [latexStr '-------------------' char(13)];
    latexStr = [latexStr '\hline' char(13)];
    latexStr = [latexStr 'Country & $w$ & $\sigma_w$ & p($H_0:w=1$) & $w_{{\rm delay}}$ & $\sigma_{w-delay}$ & p(same $H_0$) \\' char(13)];
    latexStr = [latexStr '\hline' char(13)];
    latexStr = [latexStr '\hline' char(13)];
    
    for ctry = 1:length(countryStrings)

        countryStr = countryStrings{ctry};

        cj = find(strcmp(countryStr,MATdata.country));

        cdata = MATdata.deathData{cj};
        scdata = sum(cdata,1);

        firstDeath = find(scdata,1,'first');

        dataToFit = scdata(firstDeath:end);
        
        try
            fit = fitToTimeseries(dataToFit,extendTimeLength);

            StErr = exp(-abs(fit.fullSE(6)))*abs(fit.fullSE(6));
            ddStErr = exp(-abs(fit.ddSE(6)))*abs(fit.ddSE(6));

            dp = 4;

            latexStr = [latexStr [countryStr,' & ',...
                num2str(exp(-abs(fit.fullsolution(6))),dp),' & ',...
                num2str(StErr,dp),' & ',...
                num2str(fit.fullpValues(6),dp),' & ',...
                num2str(exp(-abs(fit.ddsolution(7))),dp),' & ',...
                num2str(ddStErr,dp),' & ',...
                num2str(fit.ddpValues(7),dp),...
                ' \\'] char(13)];

            times = fit.times;
            compositeFit = fit.fullFitFunction(times);
            basicFit = fit.fitFunction(times);

            if ctry == pickOutctry
                figure(1)
                set(1,'pos',[77         291        1109         414])

                subplot(1,2,1)
                semilogPlotType = false;
                countryOutplot = true;
                plotMain()
                %xlim([35 65])
                %ylim([2600 3400])

                subplot(1,2,2)
                plot(times,dataToFit - basicFit)
                hold on
                plot(times,dataToFit - compositeFit)

                xlabel('days from first recorded case')
                ylabel('death differential')
                axis tight
                legend({'data - SIR model','data - isolation model'})
                legend('boxoff')

            end

            if ctry <= 9
                figure(2)
                set(2,'pos',[32 1 1229 704])
                subplot(3,3,ctry)
                semilogPlotType = false;
                countryOutplot = false;
                plotMain()

                figure(3)
                set(3,'pos',[32 1 1229 704])
                subplot(3,3,ctry)
                plotDeriv()

            end
        catch
            disp('****************************************')
            disp([countryStr,' double logistic fit failed']);
            disp('****************************************')
        end
    end

    function plotDeriv()
        
        diffData = diff(dataToFit);
        yLimit = max([diff(fit.upperCI) diff(fit.fullupperCI) max(diffData)]);
            
        if semilogPlotType
            s = semilogy(dataToFit,0*dataToFit);
            axis tight
            delete(s)
            set(gca,'Ytick',[1 10 100 1000]);
            %ylim([10 5000])
        end
        
        rectangle('Position',[times(end) 0 extendTimeLength yLimit],...
            'facecolor',parameters.grey,'edgecolor','none')
        hold on
        
        plot(diffData,'.','markersize',26)

        axis tight
        YL = ylim;
        ylim([YL(1) yLimit]);
        RL = exp(-abs(fit.AIC - fit.fulLAIC)/2);
        XL = xlim;
        dX = XL(2) - XL(1);

        if ctry == 1
            text(times(end)-11.5,0.15*YL(2),...
                [{[num2str(extendTimeLength),'day']},{'projection'}]);
        end
        
        if ismember(ctry , [1,4,7]) || countryOutplot
            ylabel('daily deaths')
        end
        if ismember(ctry , [7,8,9]) || countryOutplot
            xlabel('days from first recorded case')
        end

        RLscale = 0.35;
        if countryOutplot
            RLscale = 0.7;
        end
        text(XL(1) + 0.075*dX,RLscale*YL(2),['RL ',num2str(RL,3)],'fontsize',12);

        plot(fit.extendTime(2:end),diff(fit.midCI),'-b','linewidth',1)
        plot(fit.extendTime(2:end),diff(fit.fullmidCI),'-k','linewidth',1)       
        
        plot(fit.extendTime(2:end),diff(fit.upperCI),':b')
        plot(fit.extendTime(2:end),diff(fit.fullupperCI),':k')    
        
        %this puts Gaussian process regression on there too:
        %plot(fit.times(2:end),diff(fit.GRP),'-.r')

        %this puts the double delay regression on there too:
        %plot(fit.times(2:end),diff(fit.ddFitFunction(fit.times)),'-.r')
        
        axis tight
        YL = ylim;
        finalYLimit = min([YL(2) 5*max(diffData)]);

        ylim([0 finalYLimit])
        
        legend({[countryStr,' data'],...
            'SIR model','SIR model upper 99% CI',...            
            'iso model', 'iso model 99% CI'},'Location','northwest')
        legend('boxoff')        
        
    end
    
    latexStr = [latexStr '\hline' char(13)];
    disp(latexStr);

    function plotMain()
    
        yLimit = max([fit.upperCI fit.fullupperCI max(dataToFit)]);
            
        if semilogPlotType
            s = semilogy(dataToFit,0*dataToFit);
            axis tight
            delete(s)
            set(gca,'Ytick',[10 100 500 1000 2000 3000 4000 5000]);
            %ylim([10 5000])
        end
        
        rectangle('Position',[times(end) 0 extendTimeLength yLimit],...
            'facecolor',parameters.grey,'edgecolor','none')
        hold on
        
        plot(dataToFit,'.','markersize',26)

        axis tight
        YL = ylim;
        ylim([YL(1) yLimit]);
        RL = exp(-abs(fit.AIC - fit.fulLAIC)/2);
        XL = xlim;
        dX = XL(2) - XL(1);

        if ctry == 1
            text(times(end)-11.5,0.15*YL(2),...
                [{[num2str(extendTimeLength),'day']},{'projection'}]);
        end
        
        if ismember(ctry , [1,4,7]) || countryOutplot
            ylabel('cumulative deaths')
        end
        if ismember(ctry , [7,8,9]) || countryOutplot
            xlabel('days from first recorded case')
        end

        RLscale = 0.35;
        if countryOutplot
            RLscale = 0.7;
        end
        text(XL(1) + 0.075*dX,RLscale*YL(2),['RL ',num2str(RL,3)],'fontsize',12);

        plot(fit.extendTime,fit.midCI,'-b','linewidth',1)
        %plot(fit.extendTime,fit.lowerCI,':b')
        plot(fit.extendTime,fit.upperCI,':b')

        plot(fit.extendTime,fit.fullmidCI,'-k','linewidth',1)
        %plot(fit.extendTime,fit.fulllowerCI,':k')
        plot(fit.extendTime,fit.fullupperCI,':k')    
        
        %this puts Gaussian process regression on there too:
        %plot(fit.times,fit.GRP,'-.r')

        %this puts the double delay regression on there too:
        %plot(fit.times,fit.ddFitFunction(fit.times),'-.r')
        
        axis tight
        YL = ylim;
        finalYLimit = min([YL(2) 5*max(dataToFit)]);

        ylim([0 finalYLimit])
        
        legend({[countryStr,' data'],...
            ['SIR model (adj R^2 \approx',num2str(fit.adjR2,4),')'],...            
            'SIR model 99% CI',...            
            ['iso model (',num2str(fit.fulladjR2,4),')']...
            'iso model 99% CI'},'Location','northwest')
        legend('boxoff')
        
        if strcmp(countryStr,'United Kingdom')
            UKestimate = [fit.midCI(end), fit.fullmidCI(end)];
            UKestimateMax99 = [fit.upperCI(end), fit.fullupperCI(end)];
        end
    end

end
