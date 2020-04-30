close all
clc

%%

if not(exist('figures','dir'))
    mkdir('figures')
end
if not(exist('synthMatFiles','dir'))
    mkdir('synthMatFiles')    
end

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

close all
testDataFits;


