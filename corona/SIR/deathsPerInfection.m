function deathsPerInfection(MATdata,highlistList)

    clc
    close all

    figure(1)
    set(1,'pos',[43   223   769   482]);
    
    %%
    
    parameters = defaulParameters();
    MM = parameters.movingAverageDays;
    
    cList = [];
    plots = [];
    
    if nargin < 2
        highlistList = {parameters.countryStrings{1:6}};        
    end
    N = length(highlistList);
    
    legends = {};
    deltaDs = [];
    for j = 1:N
        ctry = find(strcmp(highlistList{j},MATdata.country));
        deaths = movmean(sum(MATdata.deathData{ctry},1),MM);        
        infecteds = movmean(sum(MATdata.ICaseData{ctry},1),MM);

        F = find(deaths,1,'first');
        deaths = deaths(F:end);
        infecteds = infecteds(F:end);

        deaths = diff(deaths);
        infecteds = diff(infecteds);

        times = 1:length(infecteds);

        if ~isempty(deaths)

            figure(1)

            deathsPerInfectionPercent = 100*deaths ./ infecteds;
            deathsPerInfectionPercent(infecteds == 0) = 0;

            lm = fitlm(times,deathsPerInfectionPercent);
            deltaD = abs(lm.Coefficients.Estimate(2));
            pVal = lm.Coefficients.pValue(2);
            
            if N == 1
                s = 0.5;
            else
                s = (j-1)/(N-1);
            end
            col = [s 0.5 1-s];
            lw = 3;

            %if ismember(MATdata.country{ctry},highlistList)
            %    lw = 4;
            %    col = 'k';
            %plot(times,lm.feval(times),'--','color',col);
            %hold on
            %text(times(end),lm.feval(times(end)),['\Delta \approx ',num2str(deltaD,2)]);                    
            %end

            %if strcmp(MATdata.country{ctry},'China')
            %    lw = 2;
            %    col = [1 1 1]/2;
            %end

            pl = semilogy(times,deathsPerInfectionPercent,'Linewidth',lw,'color',col);
            hold on
            legends = { legends{:} , [MATdata.country{ctry} , ' (\Delta \approx ',num2str(deltaD,2),' d^{-1})'] };

            cList = [cList ctry];
            plots = [plots pl];
            deltaDs = [deltaDs deltaD];

        end

    end
    
    figure(1)
    [~,R] = sort(deltaDs,'descend');
        
    axis tight
    YL = ylim;
    ylim([0 YL(2)]);
    %ylim([0 30]);
    set(gca,'Ytick',[0.1 0.2 0.3 0.4 0.5 1 2 3 4 5 10 20 30 40 50 100]);
    
    legend(plots(R),{legends{R}},'Location','northwest','fontsize',12);
    legend('boxoff')

    xlabel('days from first recorded death')
    ylabel('deaths per new infection (%)')
    
    my_export_fig('deathsPerInfectionData.PDF')
    
end