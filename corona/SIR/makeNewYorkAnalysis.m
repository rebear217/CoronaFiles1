clc
close all

warning('off');
parameters = defaulParameters();

%%

j = 1;

states = {'New York','California'};
deltaStart = [0.005,0.005];
iGuessStrategy = [2,2];

state = states{j};

[~,~,IG_NY] = OneCountryODEAnalysis(USdata,state);
[NYdeaths , NYAIC , NYRLmatrix] = OneCountryIsolationAnalysis(USdata,state,IG_NY);
[NYAIC2,deathRate] = OneCountryProjectionAnalysis_SIR(USdata,state,deltaStart(j),iGuessStrategy(j));
NYAICvector = [NYAIC NYAIC2];
NYfullAIC = repmat(NYAICvector,length(NYAICvector),1) - repmat(NYAICvector',1,length(NYAICvector));
NYfullRLmatrix = sign(NYfullAIC).*exp(-abs(NYfullAIC)/2);

%%

allDeltas = [parameters.tryThisDeltaRange deathRate];
[~,iAIC] = min(NYAICvector);

disp([state,' death rate estimate ~ ',num2str(100*allDeltas(iAIC)),'%']);

