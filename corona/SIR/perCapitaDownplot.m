function perCapitaDownplot(MATdata,includeList)

    clc
    
    figure(1);
    set(1,'pos',[149   158   731   547])

    N = length(MATdata.country);

    cList = [];
    plots = [];
    
    if nargin < 2
        includeList = {};
    end
    
    parameters = defaulParameters();
    
    for ctry = 1:N
        Ddata = sum(MATdata.deathData{ctry},1);        
        if (Ddata(end) > parameters.limitDeaths)
            s = (ctry-1)/(N-1);
            col = [s 0.5 1-s];
            d = getPCD(Ddata,MATdata.country{ctry});
            t = 0:length(d)-1;
            lm = fitlm(t,d);
            deltaD = lm.Coefficients.Estimate(2);
            hold on
            if (lm.coefTest < 0.01 && lm.Coefficients.Estimate(2) < 0) || ...
                ismember(MATdata.country{ctry},includeList)
                l = '-';
                lw = 1;
                if ismember(MATdata.country{ctry},includeList)
                    lw = 4;
                    col = 'k';
                end
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
    
    axis tight
    legend(plots,MATdata.country{cList},'Location','northeast');
    legend('boxoff')

    xlabel('days from first recorded death')
    ylabel('deaths per day per capita')
    %title('Constant data above zero means exponential growth')
    
    if isempty(includeList)
        my_export_fig('perCapitaDeathDeclines.pdf')
    else
        my_export_fig('perCapitaDeathDeclines_IL.png')
    end
    
end
