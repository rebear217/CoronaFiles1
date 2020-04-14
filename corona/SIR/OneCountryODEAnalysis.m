function [Deathsestimate,DeathsestimateMaxCI99,initialIsolationGuess] = OneCountryODEAnalysis(MATdata,countryStr,rectangleOn)

    if nargin < 3
        rectangleOn = true;
    end
    
    if rectangleOn
        close all
        clc        
    end
    
    figure(1)
    set(1,'pos',[60     1   885   704])
    
    MM = 2;
    tryThisDeltaRange = [0.01 0.02 0.03 0.04 0.05];

    %this is the default in case a changepoint isn't detected:
	defaultChangePoint = 20;

    parameters = defaulParameters();
    countryStrings = parameters.countryStrings;
    extendtime = 40;
    
    cj = find(strcmp(countryStr,MATdata.country));
    cdata = MATdata.deathData{cj};
    
    scdata = sum(cdata,1);
    %scdata = [scdata , scdata(end) + 717];
    
    firstDeath = find(scdata,1,'first');

    dataToFit = scdata(firstDeath:end);
    deaths = dataToFit(end);
    scdataToFit = dataToFit / deaths;
    
    T = length(dataToFit);
    times = 1:T;

    fguess = find(strcmp(countryStr,countryStrings));
    Pguesses = {[0.44775 0.44341 0.28533 ],...
        [0.49097 0.48447 0.21848 ],...
        [0.40009 0.40009 0.13918 ],...
        [0.31967 0.31884 0.03864 ],...
        [0.36294 0.36257 0.17064 ],...
        [0.23272 0.23168 0.012949 ],...
        [0.41816 0.41768 0.162 ],...
        [0.41445 0.4143 0.10843 ],...
        [0.41438 0.41436 0.093819 ],...
        };

    opts = statset('MaxIter',2000,'TolX',1e-12,'TolFun',1e-12);

    if not(isempty(fguess))
        pg = Pguesses{fguess};
    else
        pg = [1 1 1]*0.5;
    end

    fit = fitnlm(times,scdataToFit,@(p,T)Solve(p,T),pg,'Options',opts)
    betterGuess = abs(fit.Coefficients.Estimate);
    
    a = abs(betterGuess(1));
    b = abs(betterGuess(2));
    c = abs(betterGuess(3));

    figure(1)
    subplot(2,2,1)
    extendtimes = [times times(end)+(1:extendtime)];
    Y = deaths*fit.feval(extendtimes);

    %policy model error bars:    
    [beta,resid,J,sigma] = nlinfit(times,scdataToFit,@(p,T)Solve(p,T),betterGuess);
    [deltaFit, delta] = nlpredci(@(p,t)Solve(p,t),extendtimes,...
        beta,resid,'Covar',sigma,'alpha',0.01);

    upperCI = deaths*(deltaFit + delta);
    lowerCI = deaths*(deltaFit - delta);
    midCI = deaths*deltaFit;

    ylimit = max([1.5*Y dataToFit]);

    if rectangleOn
        rectangle('position',[times(end) 0 extendtime ylimit],...
            'facecolor',parameters.grey,'edgecolor','none');
        hold on
    end
    
    plot(extendtimes,Y,'-k','linewidth',1);
    hold on
    %plot(times,dataToFit,'o','markersize',6,'linewidth',1);
    plot(times,dataToFit,'.','markersize',26)        
    axis tight
    ylim([0 ylimit])

    %plot(extendtimes,deaths*midCI,'-k','linewidth',1)
    plot(extendtimes,lowerCI,':k')
    plot(extendtimes,upperCI,':k')

    legend({['SIR fit (adj R^2 \approx ',num2str(fit.Rsquared.Adjusted,4),')'],...
        [countryStr,' data'],'99% CI'},'Location','northwest')
    legend('boxoff')

    ylabel('cumulative deaths')
    xlabel('days from first recorded death')

    text(5,ylimit*0.65,'$\frac{dD}{dt} = \alpha - \beta D - \gamma e^{-\mu D}$',...
        'interpreter','latex');
    text(times(end) + 1,ylimit*0.2,{[num2str(extendtime),' day'],'projection'});

    dDTF = diff(dataToFit);

    subplot(2,2,2)
    if rectangleOn
        rectangle('position',[times(end) 0 extendtime ylimit],...
            'facecolor',parameters.grey,'edgecolor','none');
        hold on
    end
    dY = diff(Y);
    plot(extendtimes(2:end),dY,'-k','linewidth',1);
    hold on
    plot(times(2:end),dDTF,'.','markersize',26)   
    plot(extendtimes(2:end),movmean(diff(upperCI),MM),':k')
    plot(extendtimes(2:end),movmean(diff(lowerCI),MM),':k')
    
    axis tight
    ylim([0 1.2*max([dY dDTF])])

    ylabel('daily deaths')
    xlabel('days from first recorded death')

    legend({'SIR fit',countryStr,'99% CI'},'Location','northeast')
    legend('boxoff')

    subplot(2,2,3)
    if rectangleOn
        rectangle('position',[times(end) 0 extendtime ylimit],...
            'facecolor',parameters.grey,'edgecolor','none');
        hold on
    end
    PCD = dDTF./dataToFit(2:end);
    
    movPCD = movmean(PCD,parameters.movingAverageDays);
    
    changePoint = findchangepts(movPCD,'Statistic','linear');
    changePoint = max([changePoint wvarchg(movPCD)]);
    if isempty(changePoint) || (changePoint == 1)
        changePoint = defaultChangePoint;
    end
    
    plot(extendtimes(2:end),dY./Y(2:end),'-k','linewidth',1);
    hold on
    plot(times(2:end),PCD,'.','markersize',26)
    %plot(extendtimes(2:end),movmean(diff(upperCI),MM)./upperCI(2:end),':k')
    %plot(extendtimes(2:end),movmean(diff(lowerCI),MM)./lowerCI(2:end),':k')
    axis tight
    ylim([0 max(PCD)*1.05]);
    
    halfchangePoint = floor(changePoint/2);
    
	plot(extendtimes(1:halfchangePoint),1./extendtimes(1:halfchangePoint),...
        '--r','linewidth',2);
    plot([extendtimes(2) changePoint],[1 1]*mean(PCD(halfchangePoint:changePoint)),...
        '--r','linewidth',2);
    
    expTail = fitnlm(changePoint:length(PCD),PCD(changePoint:end),...
        @(p,t)abs(p(1))*exp(-abs(p(2)*t)),[0.5 0],'Options',opts);
    
    decayTail = abs(expTail.Coefficients.Estimate);
    
    disp(' ')
    disp([countryStr,' decay tail : ',num2str(decayTail(2))]);
    disp(' ')
    
    plot(extendtimes(changePoint:end),expTail.feval(extendtimes(changePoint:end)),...
        '--r','linewidth',2);
    
    ylabel('per capita daily deaths')
    xlabel('days from first recorded death')

    legend({'SIR fit',countryStr,'asymptotes'},'Location','northeast')
    legend('boxoff')

    subplot(2,2,4)
    %rectangle('position',[times(end) 0 extendtime ylimit],'facecolor',parameters.grey,'edgecolor','none');
    semilogy(extendtimes(2:end),dY./Y(2:end),'-k','linewidth',1);
    hold on
    plot(times(2:end),PCD,'.','markersize',26)   
    plot(extendtimes(changePoint:end),expTail.feval(extendtimes(changePoint:end)),...
        '--r','linewidth',2);
    axis tight
    ylim([0 max(PCD)*1.05]);

    ylabel('log per capita daily deaths')
    xlabel('days from first recorded death')

    Ccj = find(strcmp('China',MATdata.country));
    Ccdata = MATdata.deathData{Ccj};    
    Cscdata = sum(Ccdata,1);    
    CfirstDeath = find(Cscdata,1,'first');    
    CdataToFit = Cscdata(CfirstDeath:end);
    
    CPCD = diff(CdataToFit)./CdataToFit(2:end);
    CexpTail = fitnlm(changePoint:length(CPCD),CPCD(changePoint:end),...
        @(p,t)abs(p(1))*exp(-abs(p(2)*t)),[0.367 0.0746],'Options',opts);
    CdecayTail = abs(CexpTail.Coefficients.Estimate);
    plot(extendtimes(1:end),CexpTail.feval(extendtimes((1:end))),...
        '--b','linewidth',2);

    disp(['China decay tail : ',num2str(CdecayTail(2))]);
    disp(' ')
    
    legend({'SIR fit',countryStr,['asym decay -',num2str(decayTail(2),3)],...
        ['China (-',num2str(CdecayTail(2),3),')']},...
        'Location','southwest')
    legend('boxoff')
    
    Deathsestimate = Y(end);
    DeathsestimateMaxCI99 = deaths*upperCI(end);            

    deltaRange = 0.001:0.001:0.2;
    
    delta = deltaRange;
    [lambda,d,S0,P0,I0,rho,Rend,Iend,Send] = parametersFromDelta(delta);
    
    delta = 100*delta;
    
    figure(2)
    loglog(delta,P0,'-k');
    hold on
    %plot(delta,S0,':k');
    plot(delta,I0);
    plot(delta,Rend,'-r');
    plot(delta,Iend);
    plot(delta,Send);
    plot([delta(1) delta(end)],[deaths deaths],'--k');
    grid on
    
    %plot(delta,lambda);
    %plot(delta,d);
    %plot(delta,rho);
    
    Df = [0 0];
	deltas = [min(deltaRange) max(deltaRange)];
    options = optimoptions('fsolve','Display','none');
    
    for j = 1:2
        delta = deltas(j);
        [lambda,d,S0,P0,I0,rho,Rend,Iend,Send] = parametersFromDelta(delta);
        fun = @(Df)(-P0 + S0 * exp(-lambda * Df./d) + Df*(d+rho)./d);

        Dfguess = fsolve(fun,Y(end),options);
        Df(j) = Dfguess;

        xlabel('deaths per infection ($\delta\%$)','interpreter','latex')
        ylabel('modelled variable (people)')
        set(gca,'Xtick',[0.1 1 10])
        set(gca,'Xticklabel',{'0.1%','1%','10%'})
    end
    
    plot(100*deltas,Df,'-b');
    
    legend({[countryStr,' P_0'],'I(0)','R_f','I_f','S_f','D_T',...
        ['D_f (',num2str(floor(Dfguess)),')']},'Location','southwest');
    legend('boxoff')
    axis tight
    
    figure(3)
    set(3,'pos',[141   141   739   557]);
	deltaRange = tryThisDeltaRange;
    deltaStr = '';
    for j = 1:length(deltaRange)
        deltaStr = [deltaStr,num2str(100*deltaRange(j)),','];
    end
    deltaStr = deltaStr(1:end-1);
    
    options = odeset('NonNegative',[1 1 1 1]);    
    
    for j = 1:length(deltaRange)
        lw = 3;
        if j < length(deltaRange)
            lw = 1;
        end
        delta = deltaRange(j);
    	[lambda,d,S0,P0,I0,rho,Rend,Iend,Send] = parametersFromDelta(delta);    
        IC = [S0 I0 0 0];
        p = [lambda rho d];
        [SIRtimes,SIRoutput] = ode113(@(t,x)modelSIR(x,p),[1,extendtimes(end)],IC,options);
        S = SIRoutput(:,1);
        I = SIRoutput(:,2);
        R = SIRoutput(:,3);
        D = SIRoutput(:,4);
        
        p1 = semilogy(SIRtimes,S,'-','color',parameters.blue,'linewidth',lw);
        hold on
        p2 = plot(SIRtimes,I,'-','color',parameters.red,'linewidth',lw);
        p3 = plot(SIRtimes,R,'-','color',parameters.green,'linewidth',lw);
        p4 = plot(SIRtimes,D,'-','color',parameters.grey/2,'linewidth',lw);
        
    end
    p5 = plot(times,dataToFit,'.','markersize',26,'color',parameters.grey/2);
    text(floor(4*extendtimes(end)/5),1.4*Df(end),['D_f ~ ',num2str(floor(Df(end)))]);

    plotList = [p1,p2,p3,p4,p5];
    legendList = {['susceptible (S - for \delta = ',deltaStr,'% - thick line last value)'],...
        'infected (I)','recovered (R)','deaths (D)',[countryStr,' deaths']};
    
    %in case this is country data, plot I cases too:
    try
        IcaseData = sum(MATdata.ICaseData{cj},1);
        p6 = plot(times(2:end),diff(IcaseData(firstDeath:end)),'o','color',parameters.red);
        plotList = [plotList p6];
        legendList = { legendList{:} , '+ve Cov2 cases' };
    catch
        disp('Could not plot infective case data')
    end
    
    grid on
    axis tight
    set(gca,'Ytick',[1e2 1e3 1e4 1e5 1e6 1e7])
    YL = ylim;
    ylim([1e2 5*YL(2)])
    xlabel('days from first recorded death')
    ylabel('people')    
    legend(plotList,legendList,'Location','northwest');
    legend('boxoff')    
    
    figure(1)
	my_export_fig(['countryODEanalysis/',countryStr,'ODEAnalysis.pdf'])

    figure(2)
	my_export_fig(['countryODEanalysis/',countryStr,'P0analysis.pdf'])

    figure(3)
	my_export_fig(['countryODEanalysis/',countryStr,'SIRdynamics.pdf'])
    
    for j = 1:length(deltaRange)
        delta = deltaRange(j);
    	[lambda,d,S0,P0,I0,rho,Rend,Iend,Send] = parametersFromDelta(delta);    

        initialIsolationGuess.S0(j) = S0;
        initialIsolationGuess.I0(j) = I0;
        initialIsolationGuess.P0(j) = P0;
        initialIsolationGuess.lambda(j) = lambda;
        initialIsolationGuess.d(j) = d;
        initialIsolationGuess.rho(j) = rho;
        initialIsolationGuess.delta(j) = delta;
        
    end
    
    initialIsolationGuess.deltaStr = deltaStr;
    
    function [lambda,d,S0,P0,I0,rho,Rend,Iend,Send] = parametersFromDelta(delta)
        lambda = delta*c/deaths;
        d = delta*c;
        S0 = b./lambda;
        P0 = S0 - (b-a)*deaths./(c*delta);
        I0 = P0 - S0;
        rho = c-d;
        Rend = rho*deaths./d;
        Iend = d*dY(end);
        Send = P0 - Rend - deaths - Iend;
    end
        
end

function solution = Solve(p,T)

    IC = 0;
    model = @(t,x)F(x,p);
    options = odeset('NonNegative',[1]);
    [soltimes,output] = ode113(model,[1,T(end)],IC,options);
    solution = interp1(soltimes,output,T,'pchip');
    
end

function Ddot = F(x,p)

    A = abs(p(1));
    B = abs(p(2));
    C = abs(p(3));

    D = x(1);    
    Ddot = A - B*exp(-abs(D)) - C*D;

end

function rhs = modelSIR(x,p)

    S = x(1);
    I = x(2);
    R = x(3);
    D = x(4);
    
    lambda = p(1);
    rho = p(2);
    d = p(3);
    
    dotS = -lambda*S*I;
    dotI = lambda*S*I - (d+rho)*I;
    dotR = rho*I;
    dotD = d*I;
    
    rhs = [dotS ; dotI ; dotR ; dotD];

end

