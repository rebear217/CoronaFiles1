function allDeathsPerInfection(MATdata,highlistList)

    clc
    close all

    figure(1)
    set(1,'pos',[43          69        1007         636]);
    figure(2)
    set(2,'pos',[43         264        1238         441]);
    
    %%
    
    parameters = defaulParameters();
    MM = parameters.movingAverageDays;
    MM = 7;
    
    cList = [];
    plots = [];
    N = length(MATdata.country);
    
    if nargin < 2
        highlistList = {};
    end
    
    legends = {};
    deltaDs = [];
    for ctry = 1:N
        deaths = movmean(sum(MATdata.deathData{ctry},1),MM);        
        infecteds = movmean(sum(MATdata.ICaseData{ctry},1),MM);
        
        F = find(deaths,1,'first');
        deaths = deaths(F:end);
        infecteds = infecteds(F:end);
        
        deaths = diff(deaths);
        infecteds = diff(infecteds);
        
        times = 1:length(infecteds);
        
        if ~isempty(deaths)
            if (max(deaths) > 50) || any(strcmp(MATdata.country{ctry},highlistList))
                
                figure(1)
                    
                deathsPerInfectionPercent = 100*deaths ./ infecteds;
                deathsPerInfectionPercent(infecteds == 0) = 0;
                
                lm = fitlm(times,deathsPerInfectionPercent);
                deltaD = abs(lm.Coefficients.Estimate(2));
                pVal = lm.Coefficients.pValue(2);
                
                s = (ctry-1)/(N-1);
                col = [s 0.5 1-s];
                lw = 1;
                if ismember(MATdata.country{ctry},highlistList)
                    lw = 4;
                    col = 'k';
                    plot(times,lm.feval(times),'--','color',[0 0 0]);
                    text(times(end),lm.feval(times(end)),['p \approx ',num2str(deltaD)]);                    
                end
                
                if strcmp(MATdata.country{ctry},'China')
                    lw = 2;
                    col = [1 1 1]/2;
                end
                
                %if (deltaD > -0.0 && pVal < 0.01) || strcmp(MATdata.country{ctry},'China')

                    pl = plot(times,deathsPerInfectionPercent,'Linewidth',lw,'color',col);
                    hold on

                    legends = { legends{:} , [MATdata.country{ctry} , ' (\Delta \approx ',num2str(deltaD,2),' d^{-1})'] };

                    cList = [cList ctry];
                    plots = [plots pl];
                    deltaDs = [deltaDs deltaD];
                    
                %end
                
                figure(2)
                subplot(1,2,1)
                plot(deaths(2:end)./deaths(1:end-1),'color',col,'linewidth',lw);                
                hold on
                subplot(1,2,2)
                plot(infecteds(2:end)./infecteds(1:end-1),'color',col,'linewidth',lw);                
                hold on
                
            end
        end

    end
    
    figure(1)
    [~,R] = sort(deltaDs,'descend');
        
    axis tight
    YL = ylim;
    %ylim([0 YL(2)]);
    ylim([0 30]);

    legend(plots(R),{legends{R}},'Location','northeastoutside','fontsize',8);
    legend('boxoff')

    xlabel('days from first recorded death')
    ylabel('deaths per new infection (%)')
    
    ste = std(deltaDs) / sqrt(length(deltaDs) - 1);
    CI = [mean(deltaDs) - 2*ste, mean(deltaDs) + 2*ste];
    disp(CI)
    
    figure(2)
    subplot(1,2,1)
    xlabel('days from first recorded death')
    ylabel('D(t+1) / D(t)')
    axis tight
    YL = ylim;
    ylim([0 YL(2)]);
    ylim([0 3]);
    
    subplot(1,2,2)
    xlabel('days from first recorded death')
    ylabel('I(t+1) / I(t)')
    axis tight
    YL = ylim;
    ylim([0 YL(2)]);
    ylim([0 3]);
    
end