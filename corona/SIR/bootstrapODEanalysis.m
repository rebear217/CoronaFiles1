function JNparameters = bootstrapODEanalysis(MATdata,countryStr)

    close all
    
    parameters = defaulParameters();

    ctry = find(strcmp(countryStr,MATdata.country));
    cdata = MATdata.deathData{ctry};
    scdata = sum(cdata,1);
    firstDeath = find(scdata,1,'first');

    dataToFit = scdata(firstDeath:end);
    deaths = dataToFit(end);
    scdataToFit = dataToFit / deaths;

    T = length(dataToFit);
    times = 1:T;
    
    extendtime = 35;
	extendtimes = [times times(end)+(1:extendtime)];

    Pguesses = parameters.Pguesses;
    Pctry = find(strcmp(countryStr,parameters.countryStrings));
    if ~isempty(Pctry)
        pg = Pguesses{Pctry};
    else
        pg = [0.44      0.43       0.25];
    end

    [~,betterGuess] = mybootStrapper(times,scdataToFit,pg);
    
    M = length(times);
    JNparameters = zeros(M,3);
    timesBS = zeros(M,M-1);
    scdataToFitBS = zeros(M,M-1);
    
    for j = 1:M-1
        timesBS(j,:) = [times(1:j) times(j+2:end)];
        scdataToFitBS(j,:) = [scdataToFit(1:j) scdataToFit(j+2:end)];
    end
    timesBS(M,:) = times(1:M-1);
    scdataToFitBS(M,:) = scdataToFit(1:M-1);
    
    fitmodels = {};
    parfor j = 1:M
        [~,bG,FM] = mybootStrapper(timesBS(j,:),scdataToFitBS(j,:),betterGuess);
        fitmodels{j} = FM;
        JNparameters(j,:) = bG;
        %jackfun = @(X)mybootStrapper(times,X,betterGuess);
        %jackstat = jackknife(jackfun,scdataToFit)
    end
    
    figure(1)
	ylimit = max(1.5*dataToFit);
    rectangle('position',[times(end) 0 extendtime ylimit],...
        'facecolor',parameters.grey,'edgecolor','none');
    hold on
    ylabel('cumulative deaths')
    xlabel('days from first recorded death')

    DiffdataToFit = diff(dataToFit);
    figure(2)
	ylimit2 = max(1.5*DiffdataToFit);
    rectangle('position',[times(end) 0 extendtime ylimit2],...
        'facecolor',parameters.grey,'edgecolor','none');
    hold on
    ylabel('daily deaths')
    xlabel('days from first recorded death')
    
    sigma = 1;
    Nsigma = 2;
    p = zeros(1,Nsigma);
    P = zeros(1,Nsigma);
    
    pLeg = {};
    PLeg = {};
    
    for n = Nsigma:-1:0
        
        pSTD = std(JNparameters,0,1);
        pSTD = repmat(pSTD,M,1);
        pSTD = n*sigma*(rand(size(JNparameters))-0.5).*pSTD;
        
        pLeg{n+1} = ['JK SIR + ' , num2str(n*sigma),'\sigma_p'];
        PLeg{n+1} = ['JK SIR + ' , num2str(n*sigma),'\sigma_p'];        
            
        JNP = JNparameters + pSTD;
        s = n/Nsigma;
        col = [1-s 0 s];
        
        for j = 1:M
            
            %prediction = deaths*fitmodels{j}(extendtimes);
            params = JNP(j,:);
            prediction = deaths*sirSolve(params,extendtimes);
            
            figure(1)
            p(n+1) = plot(extendtimes,prediction,'-','color',col,'linewidth',1);
            
            figure(2)
        	P(n+1) = plot(extendtimes(2:end),diff(prediction),'-','color',col,'linewidth',1);
            
        end
    end
    
    figure(1)
    p1 = plot(times,dataToFit,'.','markersize',32,'color',parameters.blue);
    plot(times,dataToFit,'o','markersize',10,'color',parameters.grey,'linewidth',1);
    axis tight
    ylim([0 ylimit]);
	text(times(end)+3,ylimit*0.9,{[num2str(extendtime),' day'],'projection'});
    pLeg{Nsigma+2} = countryStr;
    legend([p p1],pLeg,'Location','northwest');
    legend('boxoff')

    figure(2)    
    P1 = plot(times(2:end),DiffdataToFit,'.','markersize',32,'color',parameters.blue);    
    plot(times(2:end),DiffdataToFit,'o','markersize',10,'color',parameters.grey,'linewidth',1);
    axis tight
    ylim([0 ylimit2]);
	text(times(end)+3,ylimit2*0.9,{[num2str(extendtime),' day'],'projection'});
    PLeg{Nsigma+2} = countryStr;
    legend([P P1],PLeg,'Location','northwest');
    legend('boxoff')
    
end

function [predictor,betterGuess,fitmodel] = mybootStrapper(inTimes,inscdataToFit,pgguess)
    
    opts = statset('MaxIter',200,'TolX',1e-18,'TolFun',1e-18);
    fit = fitnlm(inTimes,inscdataToFit,@(p,t)sirSolve(p,t),...
        pgguess,'Options',opts);
    predictor = fit.feval(inTimes);
    betterGuess = abs(fit.Coefficients.Estimate);
    fitmodel = @(t)fit.feval(t);
    
end

