function outsideFilename = fitToCountryData(MATdata)

    warning('off')
    
    close all
    clc
    
    pickOutctry = 10;
	extendTimeLength = 7;
    
    parameters = defaulParameters();
	
    %'Austria',...

    countryStrings = {'China', ...
        'Denmark', ...
        'France', ...
        'Germany', ...
        'Italy', ...
        'Switzerland', ...
        'Spain', ...
        'United Kingdom', ...
        'US', ...
        'Iran'};

        outsideFilename = [countryStrings{pickOutctry},'additionalDatafit.pdf'];

    %%

    %for ctry = 1:2

    disp('latex format table:')
    disp('-------------------')
    disp('\hline')
	disp('Country & $w$ & $\sigma_w$ & p($H_0:w=1$) & $w_{{\rm delay}}$ & $\sigma_{w-delay}$ & p(same $H_0$) \\')
    disp('\hline')
    disp('\hline')
    for ctry = 1:length(countryStrings)

        countryStr = countryStrings{ctry};

        cj = find(strcmp(countryStr,MATdata.country));

        cdata = MATdata.deathData{cj};
        scdata = sum(cdata,1);

        firstDeath = find(scdata,1,'first');

        dataToFit = scdata(firstDeath:end);
        fit = fitToTimeseries(dataToFit,extendTimeLength);
        
        StErr = exp(-abs(fit.fullSE(6)))*abs(fit.fullSE(6));
        ddStErr = exp(-abs(fit.ddSE(6)))*abs(fit.ddSE(6));
        
        dp = 4;

        disp([countryStr,' & ',...
            num2str(exp(-abs(fit.fullsolution(6))),dp),' & ',...
            num2str(StErr,dp),' & ',...
            num2str(fit.fullpValues(6),dp),' & ',...
            num2str(exp(-abs(fit.ddsolution(7))),dp),' & ',...
            num2str(ddStErr,dp),' & ',...
            num2str(fit.ddpValues(7),dp),...
            ' \\']);
        
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
        end
    end

	disp('\hline')

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
        ylim([0 YL(2)])
        
        legend({[countryStr,' data'],...
            ['SIR model (adj R^2 \approx',num2str(fit.adjR2,4),')'],...            
            'SIR model 95% CI',...            
            ['iso model (',num2str(fit.fulladjR2,4),')']...
            'iso model 95% CI'},'Location','northwest')
        legend('boxoff')
        
    end

end
