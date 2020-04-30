function [cList,UKestimate] = loglogPlot(MATdata,ignoreUS,highlightCountries)

    clc
    close all
    
    figure(1);
    set(1,'pos',[37         199        1244         506])
    
    if nargin < 3
    	highlightCountries = {};
    end
    
    %places with bad data:
    forceIgnoreList = {'Guinea-Bissau','Chad','Maldives'};
    
    UKestimate = NaN;

    parameters = defaulParameters();
    forceInclude = {parameters.countryStrings{1:6} , highlightCountries{:}};
    
    N = length(MATdata.country);
    simTimeLength = 85;
    PCAdims = 5;
    minDeathsToIncludeCountry = 0*parameters.limitDeaths;
    movingAverage = parameters.movingAverageDays;
    
    legendDeathsMax = 5000;
    
    cList = [];
    j = 1;
    
    compactData = {};
    
    for ctry = 1:N
        if ~strcmp(MATdata.country{ctry},forceIgnoreList)
            Ddata = sum(MATdata.deathData{ctry},1);
            if not(strcmp('Global',MATdata.country{ctry})) && ...
                (sum(Ddata) > minDeathsToIncludeCountry) && ...
                    (not(ignoreUS) || not(strcmp('US',MATdata.country{ctry})))

                if strcmp('US',MATdata.country{ctry})
                    disp('included US');
                end
                try
                    [PCD,y,x] = getPCD(Ddata,MATdata.country{ctry});
                catch
                    error([MATdata.country{ctry},' caused an error.'])
                end
                s = (ctry-1)/(N-1);
                col = [s 0.5 1-s];
                cList = [cList ctry];
                lw = 1 + 1.0*(sum(Ddata) > 1000);

                subplot(1,2,1)
                loglog(x,y,'-','linewidth',lw,'color',col);
                hold on

                subplot(1,2,2)
                plot(PCD,'-','linewidth',lw,'color',col);
                hold on

                %compactData{j} = x;
                compactData{j} = y/max(y);
                j = j + 1;
            end
        end
    end

	subplot(1,2,1)
    legend({MATdata.country{cList}},'Location','northwest');
    legend('boxoff')

    xlabel('cumulative deaths')
    ylabel('deaths each day')
    axis tight

	subplot(1,2,2)
    legend(MATdata.country{cList},'Location','northeast');
    legend('boxoff')

    xlabel('days from first death')
    ylabel('per capita deaths each day')
    axis tight
    
    %%
    
    figure(2)
    set(2,'pos',[66     1 920   704])
    
    lenmx = 1;
    for c = 1:length(cList)
        l = length(compactData{c});
        lenmx = max(l,lenmx);
    end
    
    NN = NaN(length(cList),lenmx);
    smallTimes = 1:min([simTimeLength lenmx]);
    
    for c = 1:length(cList)
        NN(c,1:length(compactData{c})) = compactData{c};
    end
    
    [C, ss, M, X,Ye] = ppca_mv(NN,PCAdims,0);
    
    LegText = {};
    
    pl = zeros(1,length(cList));
    pl2 = zeros(1,length(cList));
    
	iC = 1;
    
    for c = 1:length(cList)
        ctry = cList(c);
        
        s = (c-1)/(length(cList)-1);
        col = [s 0.5 1-s];
        lw = 1;
        if movingAverage > 0
            Y = movmean(Ye(c,smallTimes),movingAverage);
        else
            Y = Ye(c,smallTimes);
        end

        Ddata = sum(MATdata.deathData{ctry},1);
        [PCD,y,x] = getPCD(Ddata,MATdata.country{ctry});
        Y = max(y)*Y;
        CY = cumsum(Y);
        
        %cumsum(sum(Ye,1))        
        
        if any(strcmp(MATdata.country{cList(c)},highlightCountries))
            lw = 2;
            col = [1 1 1]*0;
        end
        
        if strcmp(MATdata.country{ctry},'United Kingdom')
            UKestimate = CY(end);
        end
        
        if CY(end) > legendDeathsMax || any(find(strcmp(MATdata.country{ctry},forceInclude)))
            
            subplot(2,1,1)
            CYCD = max(y)*cumsum(compactData{c});
            if iC == 1
                LegText0{iC} = [MATdata.country{ctry},' ( max \approx ',num2str(floor(max(diff(CY)))),' )'];
            else
                LegText0{iC} = [MATdata.country{ctry},' ( ',num2str(floor(max(diff(CY)))),' )'];
            end
            pl(iC) = plot(smallTimes,Y,'-','color',col,'linewidth',lw);
            hold on
            
            if any(find(strcmp(MATdata.country{ctry},forceInclude)))
                try
                    plot(smallTimes(1:length(CYCD)-1),diff(CYCD),'.','color',col,'markersize',24);
                catch
                    plot(smallTimes(2:end),diff(CYCD(1:length(smallTimes))),'.','color',col,'markersize',24);
                end
            end

            subplot(2,1,2)
            pl2(iC) = plot(smallTimes,CY,'-','color',col,'linewidth',lw);
            hold on
            if iC == 1
                LegText{iC} = [MATdata.country{ctry},' ( total \approx ',num2str(floor(CY(end))),' )'];
            else
               LegText{iC} = [MATdata.country{ctry},' ( ',num2str(floor(CY(end))),' )'];
            end
               
            iC = iC + 1;
            if any(find(strcmp(MATdata.country{ctry},forceInclude)))
                try
                    plot(smallTimes(1:length(CYCD)),CYCD,'.','color',col,'markersize',24);
                catch
                    pl2(iC) = plot(smallTimes,CYCD(1:length(smallTimes)),'.','color',col,'markersize',24);
                end
            end
        end
    end
    
    pl = pl(1:iC-1);
    pl2 = pl2(1:iC-1);
    
    subplot(2,1,1)
    legend(pl,{LegText0{:}},'Location','northwest','fontsize',9);
    legend('boxoff')
    axis tight
    YL = ylim;
    ylim([0 1.1*YL(2)])
    xlabel('days after first death')
    ylabel('deaths per day')
    
    subplot(2,1,2)
    legend(pl2,LegText,'Location','northwest','fontsize',9);
    legend('boxoff')
    xlabel('days after first death')
    ylabel('cumulative deaths')
    axis tight
    YL = ylim;
    ylim([0 1.1*YL(2)])
    
    str = '';
    if ignoreUS
        str = 'IgnoreUS';
    end
    
    %{
    figure(1)
    my_export_fig(['perCapitaLogLogPlot',str,'.png'])
    my_export_fig(['perCapitaLogLogPlot',str,'.pdf'])
    %}

    figure(2)
    %my_export_fig(['PCAPlot',str,'.png'])
    my_export_fig(['PCAPlot',str,'.pdf'])

end

