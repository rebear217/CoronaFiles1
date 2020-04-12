
if not(exist('figures','dir'))
    mkdir('figures')
end
if not(exist('synthMatFiles','dir'))
    mkdir('synthMatFiles')    
end

%%

load('../data/deathData.mat')
clc
close all
perCapitaDownplot(MATdata);
%perCapitaDownplot(MATdata,{'US'});

%%

close all
perCapitaDownplotWExpFit(MATdata,{},true);
close all
perCapitaDownplotWExpFit(MATdata,{},false);
%perCapitaDownplotWExpFit(MATdata,{'United Kingdom','Italy'});

%%

clc
close all
ignoreChina = {true,false};
[~,loglogestimateNoChina] = loglogPlot(MATdata,ignoreChina{1});
[~,loglogestimate] = loglogPlot(MATdata,ignoreChina{2});

%%

clc
close all
[ODEestimate,ODEestimate99] = ODEAnalysis(MATdata);

%%

close all
generateExemplars;

%%

close all
generatePlots;

%%

close all
generateHeatmap;

%%

close all
trialOneIntervention;

%%

%this displays a p-values table in latex:
close all
[DSDestimate,DSDestimate99,latexStr] = doSomeDataFitting;

%%

close all
testDataFits;

%%

close all
basicFits([],true);

%%

close all
OneCountryODEAnalysis(MATdata,'Spain');
close all
OneCountryODEAnalysis(MATdata,'United Kingdom');

%%

disp(latexStr);

%%

disp('')
disp('')

disp('------------------------')
disp('Prediction stats...')
disp('------------------------')

UKestimateVector = [loglogestimateNoChina,...
    loglogestimate,...
    DSDestimate(1),...
    DSDestimate(2),...
    DSDestimate99(1),...
    DSDestimate99(2),...    
    ODEestimate,...
    ODEestimate99];

disp(['Vector : ',num2str(UKestimateVector,5)]);
disp(['Mean : ',num2str(mean(UKestimateVector),5)]);
disp(['Std: ',num2str(std(UKestimateVector),5)]);

%%

disp('')
disp('')

disp('------------')
disp('Finished...')
disp('------------')

