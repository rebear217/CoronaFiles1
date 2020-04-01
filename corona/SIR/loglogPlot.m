function cList = loglogPlot(MATdata,ignoreChina)

    clc
    close all
    
    figure(1);
    set(1,'pos',[37         199        1244         506])
    
    if nargin < 2
    	ignoreChina = false;
    end

    N = length(MATdata.country);
    simTimeLength = 40;
    PCAdims = 15;
    maxDeathsToIncludeCountry = 50;
    movingAverage = 0;
    highlightUS = false;
    highlightUK = false;
    
    cList = [];
    j = 1;
    compactData = {};
    for ctry = 1:N
        Ddata = sum(MATdata.deathData{ctry},1);
        F = find(Ddata,1,'first');
        if (Ddata(end) > maxDeathsToIncludeCountry) && (not(ignoreChina) || not(strcmp('China',MATdata.country{ctry})))
            s = (ctry-1)/(N-1);
            col = [s 0.5 1-s];
            x = Ddata(F:end);
            y = diff(x);
            if any(y < 0) || any(x < 0)
                disp(['-ve data in ',MATdata.country{ctry}])
            else
                cList = [cList ctry];
                x = Ddata((F+1):end);
                lw = 1 + 1.0*(Ddata(end) > 1000);

                subplot(1,2,1)
                loglog(x,y,'-','linewidth',lw,'color',col);
                hold on
                subplot(1,2,2)
                plot(y./x,'-','linewidth',lw,'color',col);
                hold on
                %compactData{j} = x;
                compactData{j} = y;
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
        ctry = cList(c);
        NN(c,1:length(compactData{c})) = compactData{c};
    end
    
    [C, ss, M, X,Ye] = ppca_mv(NN,PCAdims,0);
    
    LegText = {MATdata.country{cList}};
    
    pl = zeros(1,length(cList));
    pl2 = zeros(1,length(cList));
    
    for c = 1:length(cList)
        s = (c-1)/(length(cList)-1);
        col = [s 0.5 1-s];
        lw = 1;
        if movingAverage > 0
            Y = movmean(Ye(c,smallTimes),movingAverage);
        else
            Y = Ye(c,smallTimes);
        end

        if strcmp(MATdata.country{cList(c)},'United Kingdom') && highlightUK
            lw = 2;
            col = [0 0 0];
        end
        if strcmp(MATdata.country{cList(c)},'US') && highlightUS
            lw = 2;
            col = [0 0 0];
        end
        
        subplot(2,1,1)
        pl(c) = plot(smallTimes,Y,'-','color',col,'linewidth',lw);
        hold on
        
        subplot(2,1,2)
        CY = cumsum(Y);
        plot(smallTimes,CY,'-','color',col,'linewidth',lw);
        hold on
        CYCD = cumsum(compactData{c});
        try
            pl2(c) = plot(smallTimes(1:length(CYCD)),CYCD,'.','color',col,'markersize',24);
        catch
            pl2(c) = plot(smallTimes,CYCD(1:length(smallTimes)),'.','color',col,'markersize',24);
        end
        LegText{c} = [LegText{c},' data ( ',num2str(floor(CY(end))),' )'];
    end
    
    subplot(2,1,1)
    legend(pl,{MATdata.country{cList}},'Location','northwest',...
        'fontsize',7);
    legend('boxoff')
    axis tight
    YL = ylim;
    ylim([0 YL(2)])
    xlabel('days after first death')
    ylabel('deaths per day')
    
    subplot(2,1,2)
    legend(pl2,LegText,'Location','northwest',...
        'fontsize',7);
    legend('boxoff')
    xlabel('days after first death')
    ylabel('cumulative deaths')
    axis tight
    
    str = '';
    if ignoreChina
        str = 'IgnoreChina';
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

