function [Deathsestimate,AIC,RLMatrix] = OneCountryIsolationAnalysis(MATdata,countryStr,...
    initialIsolationGuess,lockdownDate,iguessStrategy)
    
    close all
    clc

    figure(1)
    set(1,'pos',[60     1   885   704])
    
    parameters = defaulParameters();
    countryStrings = parameters.countryStrings;
    deltaStr = initialIsolationGuess.deltaStr;
    
    extendtime = 50;
    Nfits = length(initialIsolationGuess.delta);

    plottedSomething = false;
    
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
    %scdata = [scdata , scdata(end) + 596];
    
    firstDeath = find(scdata,1,'first');
    isoTime = absolutelockdownStart - firstDeath + 1;
    
    dataToFit = scdata(firstDeath:end);
    %deaths = dataToFit(end);
    
    scdataToFit = dataToFit;
    %scdataToFit = dataToFit / deaths;
    
    T = length(dataToFit);
    times = 1:T;

	PDFselectionCriterion = zeros(1,Nfits);

    for ig = Nfits:-1:1
        iIG = unpack(initialIsolationGuess,ig);
        standardGuess = [iIG.lambda , 0.5 , iIG.d , 0.5 , iIG.rho];
        if (ig == Nfits)
            betterGuess = standardGuess;
        else
            %2 initial guess strategies:
            if iguessStrategy(ig) == 1
                betterGuess = standardGuess;
            else
                betterGuess = abs(thisFit.Coefficients.Estimate);
            end
        end

        opts = statset('MaxIter',120,'TolX',1e-20,'TolFun',1e-20,'Display','iter');
        thisFit = fitnlm(times,scdataToFit,@(p,T)IsirSolve(p,T,isoTime,iIG),betterGuess,'Options',opts)
        Fit{ig} = thisFit;         
        PDFselectionCriterion(ig) = thisFit.ModelCriterion.AIC;
    end
    
    AIC = [];
    Deathsestimate = zeros(1,Nfits);
    RLMatrix = zeros(Nfits);
    plotData = true;
    ylimit = zeros(1,Nfits);
    
    for ig = Nfits:-1:1
        fit = Fit{ig};
        if (fit.Rsquared.Adjusted > 0.97)
            plottedSomething = true;
            iIG = unpack(initialIsolationGuess,ig);

            opts = statset('MaxIter',100,'TolX',1e-10,'TolFun',1e-10);

            extendtimes = [times times(end)+(1:extendtime)];
            Y = fit.feval(extendtimes);

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
                R1 = rectangle('position',[times(end) 0 extendtime Ylimit],...
                    'facecolor',parameters.grey,'edgecolor','none');
                hold on
            else
                R1.Position = [times(end) 0 extendtime Ylimit];
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
                ', \delta = ',deltaStr,'%)'],...
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
            YL = ylim;
            p11 = plot([isoTime isoTime],[YL(1) YL(2)],':r','linewidth',1);
            
            %plot(extendtimes(2:end),movmean(diff(upperCI),MM)./upperCI(2:end),':k')
            %plot(extendtimes(2:end),movmean(diff(lowerCI),MM)./lowerCI(2:end),':k')
            axis tight
            ylim([0 max(PCD)*1.05]);

            halfchangePoint = floor(changePoint/2);

            plot(extendtimes(1:halfchangePoint),1./extendtimes(1:halfchangePoint),...
                '--r','linewidth',2);
            plot([extendtimes(2) changePoint],[1 1]*mean(PCD(halfchangePoint:changePoint)),...
                '--r','linewidth',2);

            if ~exist('expTail','var')
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

            legend([p1 p2 p3 p11],...
                {'i-SIR fit',countryStr,'asymptotes','lockdown'},'Location','northeast')
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

            if ig == Nfits || ~exist('p4','var')
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

            figure(ig+1)
            set(ig+1,'pos',[40 1 1100 704]);

            lw = 2;

            IC = [iIG.S0, iIG.I0, 0, 0, 0, 0, 0]';
            options = odeset('NonNegative',[1 1 1 1 1 1 1]);            
            params = abs(fit.Coefficients.Estimate);     

            [SIRtimes,SIRoutput] = ode113(@(t,x)isolationSIRmodel(t,x,params,isoTime),...
                [1,extendtimes(end)],IC,options);

            S = SIRoutput(:,1);
            I = SIRoutput(:,2);
            R = SIRoutput(:,3);
            Si = SIRoutput(:,4);
            Ii = SIRoutput(:,5);
            Ri = SIRoutput(:,6);        
            D = SIRoutput(:,7);

            P1 = semilogy(SIRtimes,S,'-','color',parameters.blue,'linewidth',lw);
            hold on
            P2 = plot(SIRtimes,I,'-','color',parameters.red,'linewidth',lw);
            P3 = plot(SIRtimes,R,'-','color',parameters.green,'linewidth',lw);
            P4 = plot(SIRtimes,Si,'--','color',parameters.blue,'linewidth',lw);
            P5 = plot(SIRtimes,Ii,'--','color',parameters.red,'linewidth',lw);
            P6 = plot(SIRtimes,Ri,'--','color',parameters.green,'linewidth',lw);
            P7 = plot(SIRtimes,D,'-','color',parameters.grey/2,'linewidth',lw);

            P8 = plot(times,dataToFit,'.','markersize',26,'color',parameters.grey/2);
            
            totR = Ri(end)+R(end);
            text(SIRtimes(end)-20,1.3*D(end),['D_f \approx ',num2str(floor(D(end)))]);
            text(SIRtimes(end)-20,1.3*totR,['R_f \approx ',num2str(floor(totR(end)))]);

            plotList = [P1,P2,P3,P4,P5,P6,P7,P8];
            legendList = {['susceptible (S - for \delta = ',num2str(100*iIG.delta),'%)'],...
                'infected (I)','recovered (R)','i-susceptible (S_i)',...
                'i-infected (I_i)','i-recovered (R_i)',...
                ['deaths (D: i-SIR AIC \approx ',num2str(PDFselectionCriterion(ig),3),')'],...
                [countryStr,' deaths']};

            %in case this is country data, plot I cases too:
            try
                IcaseData = sum(MATdata.ICaseData{cj},1);
                dailyIcases = diff(IcaseData(firstDeath:end));
                P9 = plot(times(2:end),dailyIcases,'o','color',parameters.red);
                plotList = [plotList P9];
                legendList = { legendList{:} , '+ve Cov2 cases' };
                
                figure(20)
                set(20,'pos',[40   237   806   468]);
                Mlist = {'.','*','s','o','x','p','h'};
                MS = [28 10 10 10 10 10 10];
                S = (ig-1)/(Nfits-1);
                rColor = [S 0.5 (1-S)];
                totalI = I + Ii;
                esttotalIdata = max(totalI) * dailyIcases/max(dailyIcases);
                IP1 = plot(SIRtimes,totalI,'-','color',rColor,'linewidth',1);
                hold on
                IP2 = plot(times(2:end),esttotalIdata,'.','Marker',Mlist{ig},...
                    'markersize',MS(ig),'color',rColor);
                legend([IP1 IP2],{'total I model predictions','scaled data I'});
                legend('boxoff')
                
            catch
                disp('Could not plot infective case data')
            end

            figure(ig+1)
            grid on
            axis tight
            set(gca,'Ytick',[1e0 1e1 1e2 1e3 1e4 1e5 1e6 1e7 1e8])
            YL = ylim;
            ylim([1e0 5*YL(2)])
            xlabel('days from first recorded death')
            ylabel('people')
            legend(plotList,legendList,'Location','northwest');
            legend('boxoff')    

            Deathsestimate(ig) = Y(end);        
            plotData = false;
        else
            disp('');
            disp(['Did not include plot of delta = ',...
                num2str(initialIsolationGuess.delta(ig))]);
            disp('');
        end        
    end
    
    
    if plottedSomething
        figure(20)
        xlabel('days from first recorded death')
        ylabel('predicted total infecteds')
        axis tight
        YL = ylim;
        ylim([0 YL(2)*1.1]);
        
        figure(1)
        %[MD,MDi] = max(Deathsestimate);
        [~,MDi] = min(PDFselectionCriterion);
        MD = Deathsestimate(MDi);
        AIC = PDFselectionCriterion;
        AICm = repmat(AIC,5,1)- repmat(AIC',1,5);
        RLMatrix = exp(-abs(AICm)/2);
        
        disp('')
        disp('---------------------------')
        disp('iSIR fit AIC values are ... ');
        disp(num2str(PDFselectionCriterion));
        disp('This figure has the lowest AIC :')
        disp(num2str(MDi+1));
        disp('---------------------------')
        disp('RL matrix is ...')
        disp(RLMatrix)
        disp('---------------------------')
        disp('')        
        
        subplot(2,2,1);
        text(15,1.1*MD,{['AIC optimal ',num2str(floor(MD)),' deaths']...
            ['\delta = ',num2str(100*initialIsolationGuess.delta(MDi)),'%']},'FontSize',10);

        my_export_fig(['countryODEanalysis/',countryStr,'IsolationODEAnalysis.pdf'])    

        figure(MDi+1)
        my_export_fig(['countryODEanalysis/',countryStr,'IsolationODEAnalysisISIR.pdf'])    
    end
    
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