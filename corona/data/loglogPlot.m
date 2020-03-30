function loglogPlot(MATdata)

    figure(1);
    set(1,'pos',[57     8   894   697])

    N = length(MATdata.country);

    cList = [];
    for ctry = 1:N
        Ddata = sum(MATdata.deathData{ctry},1);
        if (Ddata(end) > 30)
            cList = [cList ctry];
            s = (ctry-1)/(N-1);
            col = [s 0.5 1-s];
            F = find(Ddata,1,'first');
            x = Ddata(F:end);
            y = diff(x);
            x = Ddata((F+1):end);
            lw = 1 + 1.0*(Ddata(end) > 1000);
            loglog(x,y,'-','linewidth',lw,'color',col);
            hold on
        end
    end

    legend(MATdata.country{cList},'Location','northwest');
    legend('boxoff')

    xlabel('cumulative deaths')
    ylabel('deaths each day')

end
