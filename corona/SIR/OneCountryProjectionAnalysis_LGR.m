function OneCountryProjectionAnalysis_LGR(MATdata,countryStr)

    parameters = defaulParameters();

    cj = find(strcmp(countryStr,MATdata.country));

    cdata = MATdata.deathData{cj};
    scdata = sum(cdata,1);

    firstDeath = find(scdata,1,'first');

    dataToFit = scdata(firstDeath:end);
    times = 1:length(dataToFit);
    ET = 7;
    extendTimes = 1:(times(end)+ET);
    
    truncatedTimes = 14:7:length(dataToFit);
    truncatedTimes(end) = length(dataToFit);
    
    T = length(dataToFit);
    
    figure(1)
    set(1,'pos',[70         254        570         450])
    
    sslp = semilogy(times,dataToFit,'.k','markersize',24);
    hold on
    xlabel('days after first recorded death');
    ylabel('cumulative deaths (people)');
    axis tight
    YL = ylim;
    ylim([1 YL(2)*1.5]);
    
    figure(2)
    set(2,'pos',[70           1        1155         703])
    
    Dsslp = semilogy(times(2:end),diff(dataToFit),'.k','markersize',24);
    hold on
    xlabel('days after first recorded death');
    ylabel('daily deaths (people)');
    axis tight
    YL = ylim;
    ylim([1 YL(2)*1.5]);    
    
    for fd = 1:length(truncatedTimes)
        tr = truncatedTimes(fd);
        
        truncDataToFit = dataToFit(1:tr);
        DtruncDataToFit = diff(truncDataToFit);
        
        s = (fd-1) / (length(truncatedTimes) - 1);
        col = s*parameters.grey/2 + (1-s)*[1 0 0];
        
        fit = fitToTimeseries(truncDataToFit,T);

        ftimes = fit.times;
        compositeFit = fit.fullFitFunction(extendTimes);
        basicFit = fit.fitFunction(extendTimes);

        projectedDeaths = compositeFit;
        
        figure(1)
        plot(extendTimes,basicFit,':','color',col,'linewidth',1);
        PP = plot(extendTimes,compositeFit,'-','color',col,'linewidth',1);
        
        pp(fd) = plot(tr(end),truncDataToFit(end),'o','markersize',14,...
            'color',col,'linewidth',3);

        pcerrorEstimates = 100*mean(projectedDeaths(T-7:T)./dataToFit(T-7:T));
        error = floor(pcerrorEstimates);
        if abs(error) > 1000
            errorStr = 'too large';
        else
            errorStr = [' \approx ',num2str(error),'%'];
        end
        ll{fd} = [num2str(tr),' days : \Delta_r ',errorStr];
                    
        figure(2)
        %plot(extendTimes(2:end),diff(basicFit),':','color',col,'linewidth',1);
        DPP = plot(extendTimes(2:end),diff(compositeFit),'-','color',col,'linewidth',1);
        
        Dpp(fd) = plot(tr(end),DtruncDataToFit(end),'o','markersize',14,...
            'color',col,'linewidth',3);

        if abs(error) > 1000
            errorStr = 'too large';
        else
            errorStr = [' \approx ',num2str(error),'%'];
        end
        ll{fd} = [num2str(tr),' days : \Delta_r ',errorStr];
    end
    
    figure(1)
	legend([sslp PP pp],{ [countryStr,' cumulative deaths'], ...
        ['projections to ',num2str(extendTimes(end)),' days'] , ...
        ll{:}},'location','southeast');    
    legend('boxoff')
	my_export_fig(['countryODEanalysis/',countryStr,'LGRHindcasting.pdf'])    

    figure(2)
	legend([Dsslp DPP Dpp],{ [countryStr,' cumulative deaths'], ...
        ['projections to ',num2str(extendTimes(end)),' days'] , ...
        ll{:}},'location','northwest');    
    legend('boxoff')

end