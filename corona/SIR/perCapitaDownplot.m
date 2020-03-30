function perCapitaDownplot(MATdata)

    clc
    
    figure(1);
    set(1,'pos',[149   158   731   547])

    N = length(MATdata.country);

    cList = [];
    plots = [];
    
    for ctry = 1:N
        Ddata = sum(MATdata.deathData{ctry},1);        
        if (Ddata(end) > 50)
            s = (ctry-1)/(N-1);
            col = [s 0.5 1-s];
            F = find(Ddata,1,'first');
            x = Ddata(F:end);
            dx = diff(x);
            x = Ddata((F+1):end);
            if any(dx < 0)
                disp(['-ve ddata in ',MATdata.country{ctry}])
            else
                y = (dx./x);
                t = 0:length(x)-1;
                %lw = 1 + 2.0*(Ddata(end) > 100);
                d = movmean(y,5);
                lm = fitlm(t,d);
                deltaD = lm.Coefficients.Estimate(2);
                hold on
%                if strcmp(MATdata.country{ctry}(1:2),'US') || ...
%                        lm.coefTest < 0.05 ...
%                        && lm.Coefficients.Estimate(2) < 0
                if lm.coefTest < 0.05 && lm.Coefficients.Estimate(2) < 0
                    l = '-';
                    lw = 1;
                    %if strcmp(MATdata.country{ctry}(1:2),'US')
                    %    lw = 4;
                    %end
                    %if strcmp(MATdata.country{ctry}(1:2),'Un')
                    %    col = 'k';
                    %    lw = 2;
                    %end
                    if strcmp(MATdata.country{ctry}(1:2),'Ch')
                        lw = 3;
                    end
                    pl = plot(t,d,l,'linewidth',lw,'color',col);
                    text(t(end)+1,d(end),MATdata.country{ctry});
                    plot(t(end),d(end),'.','color',col,'markersize',50);
                    cList = [cList ctry];
                    plots = [plots pl];
                    
                    MATdata.country{ctry} = [MATdata.country{ctry},...
                        ' (\Delta \approx ',num2str(deltaD,2),' d^{-1})'];
                else
                    l = ':';
                end
            end
        end
    end
    
    axis tight
    legend(plots,MATdata.country{cList},'Location','northeast');
    legend('boxoff')

    xlabel('days from first recorded death')
    ylabel('deaths per day per capita')
    %title('Constant data above zero means exponential growth')

    my_export_fig('perCapitaDeathDeclines.pdf')

end
