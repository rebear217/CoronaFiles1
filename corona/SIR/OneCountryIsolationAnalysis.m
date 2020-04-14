function Deathsestimate = OneCountryIsolationAnalysis(MATdata,countryStr,...
    initialIsolationGuess,lockdownDate,iguessStrategy)
    
    close all
    clc

    figure(1)
    set(1,'pos',[60     1   885   704])
    
    parameters = defaulParameters();
    countryStrings = parameters.countryStrings;
    extendtime = 50;
    Nfits = length(initialIsolationGuess.delta);

	%MM = 2;

    fc = find(strcmp(countryStr,countryStrings));
    if ~isempty(fc)
        absolutelockdownStart = parameters.LockdownStart(fc);
    else
        absolutelockdownStart = NaN;
    end
    if nargin >= 4
        if ~isempty(lockdownDate)
            absolutelockdownStart = lockdownDate;
        end
    end
    
    if nargin < 5
        iguessStrategy = ones(Nfits,1);
    end
    if length(iguessStrategy) < Nfits
        error('Error - initial guess strategy has wrong number of entries');
        %iguessStrategy = ones(Nfits,1);        
    end
    if ~(iguessStrategy(Nfits) == 1)
        error('Error - first initial guess assumes a prior solution is known when it is not');        
        %iguessStrategy = ones(Nfits,1);        
    end
    
    
    %this is the default in case a changepoint isn't detected:
	defaultChangePoint = 20;
    
    cj = find(strcmp(countryStr,MATdata.country));
    cdata = MATdata.deathData{cj};
    
    scdata = sum(cdata,1);
    %scdata = [scdata , scdata(end) + 778];
    
    firstDeath = find(scdata,1,'first');
    isoTime = absolutelockdownStart - firstDeath + 1;
    
    dataToFit = scdata(firstDeath:end);
    deaths = dataToFit(end);
    scdataToFit = dataToFit / deaths;
    
    T = length(dataToFit);
    times = 1:T;
        
    for ig = Nfits:-1:1
        iIG = unpack(initialIsolationGuess,ig);
        if (ig == Nfits)
            betterGuess = [iIG.lambda , 0 , iIG.d , 0 , iIG.rho];
        else
            %2 initial guess strategies:
            if iguessStrategy(ig) == 1
                betterGuess = [iIG.lambda , 0 , iIG.d , 0 , iIG.rho];
            else
                betterGuess = abs(thisFit.Coefficients.Estimate);
            end
        end

        opts = statset('MaxIter',100,'TolX',1e-20,'TolFun',1e-20,'Display','iter');
        thisFit = fitnlm(times,scdataToFit,@(p,T)Solve(p,T,isoTime,iIG),betterGuess,'Options',opts)
        Fit{ig} = thisFit;         
    end
    
    Deathsestimate = zeros(1,Nfits);
    plotData = true;
    ylimit = zeros(1,Nfits);
    for ig = Nfits:-1:1
        fit = Fit{ig};
        %iIG = unpack(initialIsolationGuess,ig);

        opts = statset('MaxIter',100,'TolX',1e-10,'TolFun',1e-10);

        extendtimes = [times times(end)+(1:extendtime)];
        Y = deaths*fit.feval(extendtimes);

        %policy model error bars:    
        %[beta,resid,J,sigma] = nlinfit(times,scdataToFit,@(p,T)Solve(p,T,isoTime,iIG),betterGuess);
        %[deltaFit, idelta] = nlpredci(@(p,T)Solve(p,T,isoTime,iIG),extendtimes,beta,resid,'Covar',sigma,'alpha',0.01);

        %upperCI = deaths*(deltaFit + idelta);
        %lowerCI = deaths*(deltaFit - idelta);
        %midCI = deaths*deltaFit;    

        ylimit(ig) = max([1.5*Y dataToFit]);
        Ylimit = max(ylimit);
        
        figure(1)
        subplot(2,2,1)

        if plotData
            R = rectangle('position',[times(end) 0 extendtime Ylimit],...
                'facecolor',parameters.grey,'edgecolor','none');
            hold on
        else
            R.Position = [times(end) 0 extendtime Ylimit];
        end

        lw = 0.5*ig;
        
        p1 = plot(extendtimes,Y,'-k','linewidth',lw);
        hold on
        
        if plotData
            p2 = plot(times,dataToFit,'.','markersize',26);
        end
        
        axis tight
        ylim([0 Ylimit])
        %plot(extendtimes,deaths*midCI,'-k','linewidth',1)
        %p3 = plot(extendtimes,lowerCI,':k');
        %plot(extendtimes,upperCI,':k')

        legend([p1 p2],{['i-SIR fit (adj R^2 \approx ',num2str(fit.Rsquared.Adjusted,4),...
            ', \delta = ',initialIsolationGuess.deltaStr,'%)'],...
            [countryStr,' data']},'Location','northwest')
        legend('boxoff')

        ylabel('cumulative deaths')
        xlabel('days from first recorded death')

        dDTF = diff(dataToFit);

        subplot(2,2,2)
        if plotData
            R2 = rectangle('position',[times(end) 0 extendtime Ylimit],...
                'facecolor',parameters.grey,'edgecolor','none');
            hold on
        else
            R2.Position = [times(end) 0 extendtime Ylimit];
        end
        
        dY = diff(Y);
        p1 = plot(extendtimes(2:end),dY,'-k','linewidth',lw);
        hold on
        if plotData
            p2 = plot(times(2:end),dDTF,'.','markersize',26);
        end
        
        %p3 = plot(extendtimes(2:end),movmean(diff(upperCI),MM),':k');
        %plot(extendtimes(2:end),movmean(diff(lowerCI),MM),':k')

        axis tight
        ylim([0 1.2*max([dY dDTF])])

        ylabel('daily deaths')
        xlabel('days from first recorded death')

        legend([p1 p2],{'i-SIR fit',countryStr},'Location','northeast')
        legend('boxoff')

        subplot(2,2,3)
        if plotData
            R3 = rectangle('position',[times(end) 0 extendtime Ylimit],...
                'facecolor',parameters.grey,'edgecolor','none');
            hold on
        else
            R3.Position = [times(end) 0 extendtime Ylimit];
        end
        
        PCD = dDTF./dataToFit(2:end);
        movPCD = movmean(PCD,parameters.movingAverageDays);

        changePoint = findchangepts(movPCD,'Statistic','linear');
        changePoint = max([changePoint wvarchg(movPCD)]);
        if isempty(changePoint) || (changePoint == 1)
            changePoint = defaultChangePoint;
        end

        p1 = plot(extendtimes(2:end),dY./Y(2:end),'-k','linewidth',1);
        hold on
        if plotData
            p2 = plot(times(2:end),PCD,'.','markersize',26);
        end
        %plot(extendtimes(2:end),movmean(diff(upperCI),MM)./upperCI(2:end),':k')
        %plot(extendtimes(2:end),movmean(diff(lowerCI),MM)./lowerCI(2:end),':k')
        axis tight
        ylim([0 max(PCD)*1.05]);

        halfchangePoint = floor(changePoint/2);

        plot(extendtimes(1:halfchangePoint),1./extendtimes(1:halfchangePoint),...
            '--r','linewidth',2);
        plot([extendtimes(2) changePoint],[1 1]*mean(PCD(halfchangePoint:changePoint)),...
            '--r','linewidth',2);

        if ig == Nfits
            expTail = fitnlm(changePoint:length(PCD),PCD(changePoint:end),...
                @(p,t)abs(p(1))*exp(-abs(p(2)*t)),[0.5 0],'Options',opts);

            decayTail = abs(expTail.Coefficients.Estimate);

            disp(' ')
            disp([countryStr,' decay tail : ',num2str(decayTail(2))]);
            disp(' ')
        end

        p3 = plot(extendtimes(changePoint:end),expTail.feval(extendtimes(changePoint:end)),...
            '--r','linewidth',2);

        ylabel('per capita daily deaths')
        xlabel('days from first recorded death')

        legend([p1 p2 p3],{'i-SIR fit',countryStr,'asymptotes'},'Location','northeast')
        legend('boxoff')

        subplot(2,2,4)

        p1 = semilogy(extendtimes(2:end),dY./Y(2:end),'-k','linewidth',1);
        hold on
        if plotData
            p2 = plot(times(2:end),PCD,'.','markersize',26);
        end
        p3 = plot(extendtimes(changePoint:end),expTail.feval(extendtimes(changePoint:end)),...
            '--r','linewidth',2);
        axis tight
        ylim([0 max(PCD)*1.05]);

        ylabel('log per capita daily deaths')
        xlabel('days from first recorded death')
        
        if ig == Nfits
            Ccj = find(strcmp('China',MATdata.country));
            Ccdata = MATdata.deathData{Ccj};    
            Cscdata = sum(Ccdata,1);    
            CfirstDeath = find(Cscdata,1,'first');    
            CdataToFit = Cscdata(CfirstDeath:end);

            CPCD = diff(CdataToFit)./CdataToFit(2:end);
            CexpTail = fitnlm(changePoint:length(CPCD),CPCD(changePoint:end),...
                @(p,t)abs(p(1))*exp(-abs(p(2)*t)),[0.367 0.0746],'Options',opts);
            CdecayTail = abs(CexpTail.Coefficients.Estimate);
        
            p4 = plot(extendtimes(1:end),CexpTail.feval(extendtimes((1:end))),...
                '--b','linewidth',2);

            disp(['China decay tail : ',num2str(CdecayTail(2))]);
            disp(' ')
        end

        legend([p1 p2 p3 p4], {'SIR fit',countryStr,['asym decay -',num2str(decayTail(2),3)],...
            ['China (-',num2str(CdecayTail(2),3),')']},...
            'Location','southwest')
        legend('boxoff')

        Deathsestimate(ig) = Y(end);
        
        plotData = false;
    end
    
    [MD,MDi] = max(Deathsestimate);
    subplot(2,2,1);
    text(0.7*extendtimes(end),1.1*MD,{[num2str(floor(MD)),' deaths']...
        ['\delta = ',num2str(100*initialIsolationGuess.delta(MDi)),'%']},'FontSize',10);
    
    figure(1)
	my_export_fig(['countryODEanalysis/',countryStr,'IsolationODEAnalysis.pdf'])    
    
end

function solution = Solve(p,T,isoTime,initialIsolationGuess)

    S0 = initialIsolationGuess.S0;
    I0 = initialIsolationGuess.I0;
    %delta = initialIsolationGuess.delta;
    
    IC = [S0, I0, 0, 0, 0, 0, 0]';
    model = @(t,x)isolationSIR(t,x,p,isoTime);
    options = odeset('NonNegative',[1 1 1 1 1 1 1]');
    [soltimes,output] = ode113(model,[1,T(end)],IC,options);
    D = output(:,7);
    solution = interp1(soltimes,D,T,'pchip');
    
end

function out = isolationSIR(t,x,p,isoTime)

    lam = abs(p(1));
    ip = abs(p(2));
    
    %we are assuming de-isolation does not occur:
    %im = abs(p(3));
    im = 0;
    
    d = abs(p(3));
    lamE = abs(p(4));
    rho = abs(p(5));

    S = x(1);
    I = x(2);
    R = x(3);
    Si = x(4);
    Ii = x(5);
    Ri = x(6);
    D = x(7);
    
    ip = ip*(t > isoTime);
    im = im*(t <= isoTime);
    
    Sdot = -lam*S*I -lamE*S*Ii - ip*S + im*Si;
    Idot = lam*S*I + lamE*S*Ii - ip*I + im*Ii - (rho+d)*I;
    Rdot = rho*I - ip*R + im*Ri;

    Sidot = ip*S - im*Si - lamE*(I+Ii)*Si;
    Iidot = ip*I - im*Ii + lamE*(I+Ii)*Si - (rho+d)*Ii;
    Ridot = ip*R - im*Ri + rho*Ii;
    
    Ddot = d*(Ii + I);
    
    out = [Sdot,Idot,Rdot,Sidot,Iidot,Ridot,Ddot]';

end

function iIG = unpack(initialIsolationGuess,j)

    iIG.S0 = initialIsolationGuess.S0(j);
    iIG.I0 = initialIsolationGuess.I0(j);
    iIG.P0 = initialIsolationGuess.P0(j);
    iIG.lambda = initialIsolationGuess.lambda(j);
    iIG.d = initialIsolationGuess.d(j);
    iIG.rho = initialIsolationGuess.rho(j);
    iIG.delta = initialIsolationGuess.delta(j);
    
end