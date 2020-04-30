function OneCountryProjectionAnalysis_PCA(MATdata,countryStr)

    close all
    clc

	parameters = defaulParameters();

	cj = find(strcmp(countryStr,MATdata.country));
    cdata = MATdata.deathData{cj};
    scdata = sum(cdata,1);
    firstDeath = find(scdata,1,'first');
    
    dataToFit = scdata(firstDeath:end);
    deaths = dataToFit(end);

    daycaseData = diff(scdata);
    dailyDeaths = daycaseData(firstDeath:end);    
    
    firstTime = 14;
    stepTime = 7;    

	figure(1)
    set(1,'pos',[70         254        550         450])
    
    slp = semilogy(1:length(dataToFit),dataToFit,'.k','markersize',24);
    hold on
    xlabel('days after first recorded death');
    ylabel('deaths (people)');
    
    reduceBy = fliplr(0:7:42);
    %reduceBy = (0:7:42);
    
    pp = [];
    ll = {};
    
    for j = 1:length(reduceBy)
        
        rb = reduceBy(j);

        truncDataToFit = dataToFit(1:end-rb);
        truncatedTimes = 1:length(truncDataToFit);
        DtruncDataToFit = diff(truncDataToFit);

        s = (j-1) / (length(reduceBy) - 1);
        col = s*parameters.grey/2 + (1-s)*[1 0 0];

        for iterate = 1:10
            projectedDeaths = BayesPCA(MATdata,countryStr,rb);

            if ~any(isnan(projectedDeaths))

                T = length(projectedDeaths);
                t = min([length(projectedDeaths) length(dataToFit)]);
                times = 1:T;

                figure(1)
                PP = plot(times,projectedDeaths,'-','color',col,'linewidth',1);

                p = plot(truncatedTimes(end),truncDataToFit(end),'o','markersize',14,...
                    'color',col,'linewidth',3);
                if iterate == 1
                    text(truncatedTimes(end)+2,0.95*truncDataToFit(end),num2str(rb));
                    pp = [pp p];
                    pcerrorEstimates = 100*mean(projectedDeaths(t-7:t)./dataToFit(t-7:t));
                    error = floor(pcerrorEstimates);
                end
                if iterate == 1
                    if abs(error) > 1000
                        errorStr = 'too large';
                    else
                        errorStr = [' \approx ',num2str(error),'%'];
                    end
                    ll = {ll{:},['reduce by ',num2str(rb),' days : \Delta_r ',errorStr]};
                end
            end
        end
        
    end
    
    figure(1)
	legend([slp PP pp],{ [countryStr,' cumulative deaths'], ...
        ['projections'] , ...
        ll{:}},'location','southeast');    
    legend('boxoff')
    axis tight
    XL = xlim;
    xlim([10 min([XL(2),length(dataToFit)+7])]);
    YL = ylim;
    ylim([1e2 YL(2)*1.1]);

    figure(1)
	my_export_fig(['countryODEanalysis/',countryStr,'PCAHindcasting.pdf'])    
    
end