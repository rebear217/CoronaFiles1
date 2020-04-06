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
loglogPlot(MATdata,ignoreChina{1});
loglogPlot(MATdata,ignoreChina{2});

%%

clc
close all
ODEAnalysis(MATdata);

%%

close all
clear variables
generateExemplars;

%%

close all
clear variables
generatePlots;

%%

close all
clear variables
generateHeatmap;

%%

close all
clear variables
trialOneIntervention;

%%

%this displays a p-values table in latex:
close all
clear variables
doSomeDataFitting;

%%

close all
clear variables
testDataFits;

%%

close all
clear variables
basicFits([],true);

%%


clear variables
close all
load('../data/deathData.mat')

doSomeDataFitting;

disp('')
disp('------------')
disp('Finished...')
disp('------------')

