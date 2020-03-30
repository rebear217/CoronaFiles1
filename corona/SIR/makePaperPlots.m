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

disp('Finished...')

