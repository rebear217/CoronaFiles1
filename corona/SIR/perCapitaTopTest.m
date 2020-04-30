function perCapitaTopTest(MATdata,country)

    clc
    close all
    
    figure(1);
    set(1,'pos',[149   158   731   547])

    N = length(MATdata.country);
    
    ctry = find(strcmp(country,MATdata.country));
    
    Ddata = sum(MATdata.deathData{ctry},1);

    d = getPCD(Ddata,MATdata.country{ctry});
    t = 0:length(d)-1;
    
    lm = fitlm(t,d);
    deltaD = lm.Coefficients.Estimate(2);

    Y = -log(d);
    Y = Y(2:end);
    t = t(2:end);
    l = '-';
    
    findchangepts(Y,'Statistic','linear','MaxNumChanges',1)
    hold on
    
    lw = 2;
    col = 'k';

    plot(t,Y,l,'linewidth',lw,'color',col);
    text(t(end)+1,Y(end),MATdata.country{ctry});
    plot(t(end),Y(end),'.','color',col,'markersize',50);

    MATdata.country{ctry} = [MATdata.country{ctry},...
        ' (\Delta \approx ',num2str(deltaD,2),' d^{-1})'];
    
    axis tight
    legend(country,'Location','northeast');
    legend('boxoff')

    xlabel('days from first recorded death')
    ylabel('deaths per day per capita')
    
end