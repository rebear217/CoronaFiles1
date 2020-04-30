function clusterPCA(MATdata,timePoints)

    clc
    close all
    
    parameters = defaulParameters();

    N = length(MATdata.country);

    plots = [];
    
    deathMinimum = parameters.limitDeaths;
    
    if nargin < 2
        timePoints = 25;
    end
    
    allData = [];
    allCountries = [];

    figure(1);
    figure(2);
    set(2,'pos',[149   158   731   547])
    
    for ctry = 1:N
        
        Ddata = sum(MATdata.deathData{ctry},1);
        
        if (Ddata(end) > deathMinimum)
            
            s = (ctry-1)/(N-1);
            col = [s 0.5 1-s];
            PCD = getPCD(Ddata,MATdata.country{ctry});
            
            if length(PCD) >= timePoints
                thisData = diff(movmean(PCD,parameters.movingAverageDays));
                thisData = thisData(1:timePoints-1);
                %thisData = (thisData - mean(thisData))/std(thisData);
                thisData = (thisData - mean(thisData));
                if std(thisData) > 0

                    %thisData = diff(thisData);
                    
                    allData = [allData ; thisData];
                    allCountries = [allCountries ; ctry];

                    figure(1);
                    plot(1:timePoints-1,thisData,'-','color',col);
                    hold on
                    
                    %plot(1:timePoints,fit.feval(1:timePoints),':','color',col);
                end
            end            
        end
    end
    
    pcaCoeffs = pca(allData');
    %pcaCoeffs = sign(pcaCoeffs).*log10(1 + abs(pcaCoeffs));

    figure(2);    
    
    plot3(pcaCoeffs(:,1),pcaCoeffs(:,2),pcaCoeffs(:,3),'.','markersize',24);
    for c = 1:length(allCountries)
        text(pcaCoeffs(c,1),pcaCoeffs(c,2),pcaCoeffs(c,3),MATdata.country{allCountries(c)});
    end
    
    %view(2)
    
    %{
    plot(pcaCoeffs(:,1),pcaCoeffs(:,2),'.','markersize',24);
    for c = 1:length(allCountries)
        text(pcaCoeffs(c,1),pcaCoeffs(c,2),MATdata.country{allCountries(c)});
    end
    %}
    
end
