function [AIC,deathRate] = OneCountryProjectionAnalysis_SIR(MATdata,countryStr,delta,iGuessStrategy)

    close all
    clc
    
    if nargin < 3
        delta = 0.02;
    end
    if nargin < 4
        iGuessStrategy = 1;
    end
    
    parameters = defaulParameters();
    countryStrings = parameters.countryStrings;
    Pguesses = parameters.Pguesses;
    extendtime = 40;

	cj = find(strcmp(countryStr,MATdata.country));
    cdata = MATdata.deathData{cj};
    scdata = sum(cdata,1);
    firstDeath = find(scdata,1,'first');

    daycaseData = diff(scdata);
    dailyDeaths = daycaseData(firstDeath:end);    
    
    fc = find(strcmp(countryStr,countryStrings));
    if ~isempty(fc)
        absolutelockdownStart = parameters.LockdownStart(fc);
    else
        absolutelockdownStart = NaN;
        disp(['Press any key: isolation time not known for ',countryStr]);
        pause
    end
    isoTime = absolutelockdownStart - firstDeath + 1;
        
    firstTime = 14;
    stepTime = 7;
    
    dataToFit = scdata(firstDeath:end);
    deaths = dataToFit(end);
    
    T = length(dataToFit);
    times = 1:T;

    fguess = find(strcmp(countryStr,countryStrings));

    if not(isempty(fguess))
        betterGuess = Pguesses{fguess};
    else
        betterGuess = [1 1 1]*0.5;
    end

    figure(1)
    set(1,'pos',[70         254        1186         450])
	subplot(1,2,1)
    
    slp = semilogy(times,dataToFit,'.k','markersize',24);
    hold on
    xlabel('days after first recorded death');
    ylabel('deaths (people)');
    
	extendtimes = [times times(end)+(1:extendtime)];
    
    allFitDays = firstTime:stepTime:T;
    %allFitDays = fliplr(T:-stepTime:firstTime);    
    allFitDays(end) = T;
    
    for fd = 1:length(allFitDays)
        
        opts = statset('MaxIter',140,'TolX',1e-20,'TolFun',1e-20);
        opts2 = statset(opts,'Display','iter');
        
        figure(1)
        subplot(1,2,1)

        s = (fd-1) / (length(allFitDays) - 1);
        col = s*parameters.grey/2 + (1-s)*[1 0 0];
        
        %fitDays = allFitDays(length(allFitDays) - fd + 1);
        fitDays = allFitDays(fd);
        
        scdataToFit = dataToFit / deaths;
        
        shorterData = scdataToFit(1:fitDays);
        shorterTimes = times(1:fitDays);
        fit = fitnlm(shorterTimes,shorterData,@(p,t)sirSolve(p,t),betterGuess,'Options',opts)
        betterGuess = abs(fit.Coefficients.Estimate);

        projectedDeaths = deaths*fit.feval(extendtimes);        

        P = plot(extendtimes,projectedDeaths,'-','linewidth',1,...
            'color',col);
        p(fd) = plot(shorterTimes(end),deaths*shorterData(end),'o','markersize',14,...
            'color',col,'linewidth',3);
        text(shorterTimes(end)+1,0.95*deaths*shorterData(end),num2str(fitDays));
        error = floor(100*(projectedDeaths(T))/dataToFit(T));
        l{fd} = [num2str(fitDays),' days : \Delta_r \approx ',num2str(error),'%'];
        
        Y = projectedDeaths;
        dY = diff(Y);
        a = abs(betterGuess(1));
        b = abs(betterGuess(2));
        c = abs(betterGuess(3));
        
    	[~,~,~,~,~,~,~,~,~,iIG] = parametersFromDelta(delta);

        legend([slp P p],{ [countryStr,' cumulative deaths'], ['projections to ',num2str(T),' days'] , l{:}},'location','southeast');
        legend('boxoff')
        axis tight
        xlim([2 T+7])
        YL = ylim;
        ylim([YL(1) 1.1*YL(2)]);

        %%

        subplot(1,2,2)

        sslp = semilogy(times,dataToFit,'.k','markersize',24);
        hold on
        xlabel('days after first recorded death');
        ylabel('deaths (people)');

        if fd == 1 || iGuessStrategy == 2
            standardGuess = [delta , iIG.lambda , 0.5 , iIG.d , 0.5 , iIG.rho];
        end
            
        %s = (fd-1) / (length(allFitDays) - 1);
        %col = (1-s)*parameters.grey/2 + s*[1 0 0];        
        %fitDays = allFitDays(fd);  
        
        shorterData = dataToFit(1:fitDays);
        shorterTimes = times(1:fitDays);
        
        fit = fitnlm(shorterTimes,shorterData,@(p,t)IsirSolveWithDelta(p,t,isoTime),...
            standardGuess,'Options',opts2)
        
        betterdelta = abs(fit.Coefficients.Estimate(1));
        standardGuess = abs(fit.Coefficients.Estimate);
        selectionCriterion = fit.ModelCriterion.AIC;
                
        projectedDeaths = fit.feval(extendtimes);        

        PP = plot(extendtimes,projectedDeaths,'-','linewidth',1,'color',col);
        
        pp(fd) = plot(shorterTimes(end),shorterData(end),'o','markersize',14,...
            'color',col,'linewidth',3);
        text(shorterTimes(end)+1,0.95*shorterData(end),num2str(fitDays));
        pcerrorEstimates = 100*mean(projectedDeaths(T-7:T)./dataToFit(T-7:T));
        error = floor(pcerrorEstimates);
        if abs(error) > 1000
            errorStr = 'too large';
        else
            errorStr = [' \approx ',num2str(error),'%'];
        end
        ll{fd} = [num2str(fitDays),' days : \Delta_r ',errorStr];
            
        legend([sslp PP pp],{ [countryStr,' cumulative deaths'], ['projections to ',num2str(T),' days'] , ll{:}},'location','southeast');
        legend('boxoff')
        axis tight
        xlim([2 T+7])
        YL = ylim;
        ylim([YL(1) 1.1*YL(2)]);
        
        figure(2)
        set(2,'pos',[40 1 1100 704]);
        
        [~,~,~,~,~,~,~,~,~,betteriIG] = parametersFromDelta(betterdelta);        
        
        IC = [betteriIG.S0, betteriIG.I0, 0, 0, 0, 0, 0]';
        options = odeset('NonNegative',[1 1 1 1 1 1 1]);            
        params = standardGuess(2:end);     

        [SIRtimes,SIRoutput] = ode113(@(t,x)isolationSIRmodel(t,x,params,isoTime),...
            [1,extendtimes(end)],IC,options);
        
        S = SIRoutput(:,1);
        I = SIRoutput(:,2);
        R = SIRoutput(:,3);
        Si = SIRoutput(:,4);
        Ii = SIRoutput(:,5);
        Ri = SIRoutput(:,6);        
        D = SIRoutput(:,7);
        
        lw = 1;
        lw2 = 0.5;
        
        if (fd == length(allFitDays))
            lw2 = 2;
            lw = 2;
        end
        
        P1 = semilogy(SIRtimes,S,'-','color',parameters.blue,'linewidth',lw2);
        hold on
        P2 = plot(SIRtimes,I,'-','color',parameters.red,'linewidth',lw2);
        P3 = plot(SIRtimes,R,'-','color',parameters.green,'linewidth',lw2);
        P4 = plot(SIRtimes,Si,'--','color',parameters.blue,'linewidth',lw2);
        P5 = plot(SIRtimes,Ii,'--','color',parameters.red,'linewidth',lw2);
        P6 = plot(SIRtimes,Ri,'--','color',parameters.green,'linewidth',lw2);
        P7 = plot(SIRtimes,D,'-','color',parameters.grey/2,'linewidth',lw2);

        P8 = plot(times,dataToFit,'.','markersize',26,'color',parameters.grey/2);

        totR = Ri(end)+R(end);
        if (fd == length(allFitDays))
            text(SIRtimes(end)-20,1.3*D(end),['D_f \approx ',num2str(floor(D(end)))]);
            text(SIRtimes(end)-20,1.3*totR,['R_f \approx ',num2str(floor(totR(end)))]);
        end
        
        plotList = [P1,P2,P3,P4,P5,P6,P7,P8];
        legendList = {['susceptible (S - for \delta = ',num2str(100*betterdelta,3),'%)'],...
            'infected (I)','recovered (R)','i-susceptible (S_i)',...
            'i-infected (I_i)','i-recovered (R_i)',...
            ['deaths (D : i-SIR AIC \approx ',num2str(selectionCriterion,3),')'],...
            [countryStr,' deaths']};

        %in case this is country data, plot I cases too:
        try
            IcaseData = sum(MATdata.ICaseData{cj},1);
            dailyIcases = diff(IcaseData(firstDeath:end));
            P9 = plot(times(2:end),dailyIcases,'o','color',parameters.red);
            plotList = [plotList P9];
            legendList = { legendList{:} , '+ve Cov2 cases' };

            figure(3)
            set(3,'pos',[40   237   806   468]);
            Mlist = {'.','*','s','o','x','p','h'};
            
            MS = 20;
            if (fd == length(allFitDays))
                MS = 30;
            end
            %S = (fd-1)/(length(allFitDays)-1);
            %rColor = [S 0.5 (1-S)];
            totalI = I + Ii;
            esttotalIdata = max(totalI) * dailyIcases/max(dailyIcases);
            IP1 = plot(SIRtimes,totalI,'-','color',col,'linewidth',lw);
            hold on
            IP2 = plot(times(2:end),esttotalIdata,'.','Marker',Mlist{1},...
                'markersize',MS,'color',col);
            legend([IP1 IP2],{'total I model predictions','scaled data I'});
            legend('boxoff')
            xlabel('days from first recorded death')
            ylabel('predicted total infecteds')
            axis tight
            YL = ylim;
            ylim([0 YL(2)*1.1]);
            
        catch
            disp('Could not plot infective case data')
        end
        
        figure(4)
        set(4,'pos',[226    27   750   486]);
        
        dDS = diff(D)./diff(SIRtimes);
        q2(fd) = plot(SIRtimes(2:end),dDS,'-','color',col,'linewidth',lw);
        hold on
        lq2{fd} = ['i-SIR day ',num2str(fitDays)];
        q3(fd) = plot(extendtimes(2:end),diff(projectedDeaths),'--','color',col,'linewidth',lw);
        lq3{fd} = ['SIR day ',num2str(fitDays)];
        
        %[mdds,idds] = max(dDS);
        %text(SIRtimes(idds+1),1.1*mdds,num2str(fitDays));
        
        if fd == length(allFitDays)
            q1 = plot(times(2:end),dailyDeaths,'.','markersize',30,'color',col);
            legend([q1 q2 q3],...
                { countryStr , lq2{:} , lq3{:} });
            legend('boxoff')
        end
        xlabel('days from first recorded death')
        ylabel('daily deaths')
        axis tight
        YL = ylim;
        ylim([0 YL(2)*1.1]);
        %ylim([0 max(dailyDeaths)*2]);
        
    end
    
    AIC = selectionCriterion;
    deathRate = betterdelta;
    
    figure(1)
    subplot(1,2,1)
    ylim([1 10*max(dataToFit)])
    subplot(1,2,2)
    ylim([1 10*max(dataToFit)])

	figure(4)
    YL = ylim;
    ylim([YL(1) min([YL(2),2*max(dDS)]) ])

    figure(2)
    grid on
    axis tight
    set(gca,'Ytick',[1e0 1e1 1e2 1e3 1e4 1e5 1e6 1e7 1e8])
    YL = ylim;
    ylim([1e0 5*YL(2)])
    YL = ylim;
	ylim([YL(1) min([1e8 YL(2)])]);
    
    xlabel('days from first recorded death')
    ylabel('people')
    legend(plotList,legendList,'Location','northwest');
    legend('boxoff')    

    figure(1)
	my_export_fig(['countryODEanalysis/',countryStr,'SIRHindcasting.pdf'])
    figure(2)
	my_export_fig(['countryODEanalysis/',countryStr,'SIRProjectionAnalysis.pdf'])
    figure(4)
	my_export_fig(['countryODEanalysis/',countryStr,'SIRProjectionDailyDeaths.pdf'])
    
    function [lambda,d,S0,P0,I0,rho,Rend,Iend,Send,iIG] = parametersFromDelta(delta)
        
        lambda = delta*c/deaths;
        d = delta*c;
        S0 = b./lambda;
        P0 = S0 - (b-a)*deaths./(c*delta);
        I0 = P0 - S0;
        rho = c-d;
        Rend = rho*deaths./d;
        Iend = d*dY(end);
        Send = P0 - Rend - deaths - Iend;
        
        iIG.S0 = S0;
        iIG.I0 = I0;
        iIG.P0 = P0;
        iIG.lambda = lambda;
        iIG.d = d;
        iIG.rho = rho;
        iIG.delta = delta;    
        
    end    
    
    function output = IsirSolveWithDelta(p,t,isoTime)
        delta = abs(p(1));
        [~,~,~,~,~,~,~,~,~,iIG] = parametersFromDelta(delta);        
        output = IsirSolve(p(2:end),t,isoTime,iIG);
    end

end

