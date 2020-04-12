function [Deathsestimate,DeathsestimateMaxCI99] = OneCountryODEAnalysis(MATdata,countryStr,rectangleOn)

    close all
    clc
    
    if nargin < 3
        rectangleOn = true;
    end
    
    figure(1)
    set(1,'pos',[60     1   885   704])
    
    MM = 2;
    
	changePoint = 30;
    if strcmp(countryStr , 'Spain')
        changePoint = 23;
    end
    if strcmp(countryStr , 'United Kingdom')
        changePoint = 33;
    end

    parameters = defaulParameters();
    countryStrings = parameters.countryStrings;
    extendtime = 35;
    
    cj = find(strcmp(countryStr,MATdata.country));
    cdata = MATdata.deathData{cj};
    
    scdata = sum(cdata,1);
    
    firstDeath = find(scdata,1,'first');

    dataToFit = scdata(firstDeath:end);
    deaths = dataToFit(end);
    scdataToFit = dataToFit / deaths;

    changePoint = min([changePoint length(dataToFit)-1]);
    
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

    fit = fitnlm(times,scdataToFit,@(p,t)Solve(p,t),pg,'Options',opts)
    betterGuess = abs(fit.Coefficients.Estimate);
    
    a = abs(betterGuess(1));
    b = abs(betterGuess(2));
    c = abs(betterGuess(3));

    figure(1)
    subplot(2,2,1)
    extendtimes = [times times(end)+(1:extendtime)];
    Y = deaths*fit.feval(extendtimes);

    %policy model error bars:    
    [beta,resid,J,sigma] = nlinfit(times,scdataToFit,@(p,t)Solve(p,t),betterGuess);
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

    legend({['SIR fit (adj R^2 \approx ',num2str(fit.Rsquared.Adjusted,4),' )'],...
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
    
    plot(extendtimes(changePoint-5:end),expTail.feval(extendtimes(changePoint-5:end)),...
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
    plot(extendtimes(changePoint-3:end),expTail.feval(extendtimes(changePoint-3:end)),...
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
        @(p,t)abs(p(1))*exp(-abs(p(2)*t)),[0.5 0],'Options',opts);
    CdecayTail = abs(expTail.Coefficients.Estimate);
    plot(extendtimes(changePoint:end),CexpTail.feval(extendtimes((changePoint:end))),...
        '--b','linewidth',2);

    legend({'SIR fit',countryStr,['asym decay -',num2str(decayTail(2),3)],'China'},...
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
	deltaRange = [0.01 0.03 0.02];
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
        
        p1 = semilogy(SIRtimes,S,'-','color',0.75*[0.5 0.5 1],'linewidth',lw);
        hold on
        p2 = plot(SIRtimes,I,'-','color',0.75*[1 0.5 0.5],'linewidth',lw);
        p3 = plot(SIRtimes,R,'-','color',0.75*[0.5 1 0.5],'linewidth',lw);
        p4 = plot(SIRtimes,D,'-','color',[1 1 1]*0.5,'linewidth',lw);
        
    end
    p5 = plot(times,dataToFit,'.','markersize',26);
    text(floor(4*extendtimes(end)/5),Df(end),['D_f ~ ',num2str(floor(Df(end)))]);
    
    grid on
    axis tight
    set(gca,'Ytick',[1e2 1e3 1e4 1e5 1e6 1e7])
    YL = ylim;
    ylim([1e2 5*YL(2)])
    xlabel('days from first recorded death')
    ylabel('people')    
    legend([p1,p2,p3,p4,p5],{'S','I','R','D',[countryStr,' data']},'Location','southeast');
    
    figure(1)
	my_export_fig(['countryODEanalysis/',countryStr,'ODEAnalysis.pdf'])

    figure(2)
	my_export_fig(['countryODEanalysis/',countryStr,'P0analysis.pdf'])

    figure(3)
	my_export_fig(['countryODEanalysis/',countryStr,'SIRdynamics.pdf'])
    
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
    
%{
function [Deathsestimate,DeathsestimateMaxCI99] = oldOneCountryODEAnalysis(MATdata,countryStr,rectangleOn)

    %close all
    %clc
    
    if nargin < 3
        rectangleOn = true;
    end
    
    figure(1)
    set(1,'pos',[60     1   885   704])
    
    MM = 2;
    changePoint = 33;

    parameters = defaulParameters();
    countryStrings = parameters.countryStrings;
    extendtime = 35;
    
    cj = find(strcmp(countryStr,MATdata.country));
    cdata = MATdata.deathData{cj};
    
    scdata = sum(cdata,1);
    
    firstDeath = find(scdata,1,'first');

    dataToFit = scdata(firstDeath:end);
    deaths = dataToFit(end);
    scdataToFit = dataToFit / deaths;

    changePoint = min([changePoint length(dataToFit)-1]);
    
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

    fit = fitnlm(times,scdataToFit,@(p,t)Solve(p,t),pg,...
        'Options',opts)
    betterGuess = abs(fit.Coefficients.Estimate);
    
    a = abs(betterGuess(1));
    b = abs(betterGuess(2));
    c = abs(betterGuess(3));

    d = deaths*a;
    lamP = a/d;
    S0 = b/a;
    rhoP = c;
    
	fit = fitnlm(times,dataToFit,@(p,t)Solve2(p,t),[d rhoP S0*d lamP],'Options',opts)
    betterGuess = abs(fit.Coefficients.Estimate);

    figure(1)
    subplot(2,2,1)
    extendtimes = [times times(end)+(1:extendtime)];
    Y = fit.feval(extendtimes);

    %policy model error bars:    
    [beta,resid,J,sigma] = nlinfit(times,dataToFit,@(p,t)Solve2(p,t),betterGuess);
    [deltaFit, delta] = nlpredci(@(p,t)Solve2(p,t),extendtimes,...
        beta,resid,'Covar',sigma,'alpha',0.01);

    upperCI = deltaFit + delta;
    lowerCI = deltaFit - delta;
    midCI = deltaFit;

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

    legend({['SIR fit (adj R^2 \approx ',num2str(fit.Rsquared.Adjusted,4),' )'],...
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
    
    plot(extendtimes(changePoint-5:end),expTail.feval(extendtimes(changePoint-5:end)),...
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
    plot(extendtimes(changePoint-3:end),expTail.feval(extendtimes(changePoint-3:end)),...
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
        @(p,t)abs(p(1))*exp(-abs(p(2)*t)),[0.5 0],'Options',opts);
    CdecayTail = abs(expTail.Coefficients.Estimate);
    plot(extendtimes(changePoint:end),CexpTail.feval(extendtimes((changePoint:end))),...
        '--b','linewidth',2);

    legend({'SIR fit',countryStr,['asym decay -',num2str(decayTail(2),3)],'China'},...
        'Location','southwest')
    legend('boxoff')
    
    Deathsestimate = Y(end);
    DeathsestimateMaxCI99 = deaths*upperCI(end);            

    delta = 0.1;

    lambda = delta*c/deaths;
    d = delta*c;
    S0 = b/lambda;
    P0 = S0 - (b-a)*deaths/(c*delta);
    rho = c-d;
    
    disp(['P0 : ',num2str(floor(P0))]);
    disp(['S0 : ',num2str(floor(S0))]);
    
    figure(1)
	my_export_fig(['countryODEanalysis/',countryStr,'ODEAnalysis.pdf'])

end
%}

function solution = Solve(p,T)

    IC = 0;
    model = @(t,x)F(x,p);
    options = odeset('NonNegative',[1]);
    [soltimes,output] = ode113(model,[1,T(end)],IC,options);
    solution = interp1(soltimes,output,T,'pchip');
    
end

function Ddot = F(x,p)

    %this has flipped B (beta) and C (gamma) around wrt the latex

    A = abs(p(1));
    B = abs(p(2));
    C = abs(p(3));

    D = x(1);    
    Ddot = A - B*exp(-abs(D)) - C*D;

end

%{
function solution = Solve2(p,T)

    IC = 0;
    model = @(t,x)F2(x,p);
    options = odeset('NonNegative',[1]);
    [soltimes,output] = ode113(model,[1,T(end)],IC,options);
    solution = interp1(soltimes,output,T,'pchip');
    
end

function Ddot = F2(x,p)

    dP0 = abs(p(1));
    rhoP = abs(p(2));
    dS0 = abs(p(3));
    lamP = abs(p(4));

    D = x(1);    
    Ddot = dP0 - rhoP*D - dS0*exp(-lamP*D);
    %Ddot = d - (rho + d)*D - S0*d*exp(-lam*D/d);

end
%}

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