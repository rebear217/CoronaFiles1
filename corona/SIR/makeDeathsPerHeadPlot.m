function makeDeathsPerHeadPlot(MATdata,wormCountry)

    clc
    close all

    %%

    figure(1)
    set(1,'pos',[127   103   894   566])
    figure(2)
    set(2,'pos',[127   103   894   566])

    %%

    timePoint = 0;

    countries = length(MATdata.country);
    deathsPerHead = zeros(1,countries);
    IPerHead = zeros(1,countries);

    for c = 1:countries
        if ~isempty(MATdata.population{c})
            thisData = sum(MATdata.deathData{c},1);
            thisData = thisData(end-timePoint);
            thisIData = sum(MATdata.ICaseData{c},1);            
            thisIData = thisIData(end-timePoint);
            deathsPerHead(c) = thisData/MATdata.population{c};
            IPerHead(c) = thisIData/MATdata.population{c};
        end
    end

    [R,I] = sort(deathsPerHead);

    for c = 1:countries
        i = I(c);
        figure(1)
        MS = 26;
        col = 'k';
        if strcmp(MATdata.country{i},'Sweden')
            MS = 50;
            col = 'r';
        end
        if c >= countries-19
            semilogy(countries-c+1,deathsPerHead(i),'.k','markersize',MS,'color',col)
            text(countries-c+1+0.25,deathsPerHead(i),MATdata.country{i})
            hold on
        end
        figure(2)    
        p1 = loglog(IPerHead(i),deathsPerHead(i),'.k','markersize',26);
        %text(IPerHead(i),deathsPerHead(i),MATdata.country{i})
        hold on
    end

    figure(1)
    axis tight
    xlim([1 20]);

    set(gca,'Xticklabel',1:countries);
    set(gca,'Xtick',1:countries);
    ylabel('covid19 deaths per person since death-1')
    xlabel('country rank')

    figure(2)
    xlabel('covid19 infections per person since death-1')
    ylabel('covid19 deaths per person since death-1')

    X = IPerHead(IPerHead > 0);
    Y = deathsPerHead(IPerHead > 0);

    X = IPerHead(deathsPerHead > 0);
    Y = deathsPerHead(deathsPerHead > 0);

    %%

    x = -6.5:0.1:-1;
    fit = fitlm(log10(X),log10(Y),'RobustOpts','on')
    plot(10.^(x),10.^(fit.feval(x)),'-k');
    
    p11 = plot(10.^x,10.^(x-2),'--k');
    p111 = plot(10.^x,10.^(x-1),'-.k');
    
    [p2,p2st,p2end] = plotTrajectory(wormCountry,'b');
    [p3,p3st,p3end] = plotTrajectory('Sweden','r');
    
    legend([p1 p11 p111 p2st p2end p3st p3end],...
        {'all countries','1% line','10% line',[wormCountry,' worm start'],[wormCountry,' worm end'],...
        'Sweden worm start','Sweden worm end'},'location','northwest');
    legend('boxoff')
    
    function [p,pst,pend] = plotTrajectory(wormCountry,colour)
        F = find(strcmp(wormCountry,MATdata.country));
        deathdata = sum(MATdata.deathData{F},1)/MATdata.population{F};
        Icasedata = sum(MATdata.ICaseData{F},1)/MATdata.population{F};
        j = find(deathdata > 0,1,'first');
        deathdata = deathdata(j:end);
        Icasedata = Icasedata(j:end);
        
        p = plot(Icasedata,deathdata,'-','color',colour);
        pst = plot(Icasedata(1),deathdata(1),'.','markersize',24,'color',colour);
        pend = plot(Icasedata(end),deathdata(end),'o','markersize',14,'color',colour);
        %text(Icasedata(end),deathdata(end),wormCountry)        
        %text(Icasedata(1),deathdata(1),wormCountry)        
    end
    axis tight

    figure(1)
    my_export_fig('countryCovidRank.PDF');
    figure(2)
    my_export_fig([wormCountry,'Covidworm.PDF']);
end

