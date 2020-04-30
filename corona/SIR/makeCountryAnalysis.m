clc
close all

warning('off');
parameters = defaulParameters();

%%

j = 6;

countries = {'United Kingdom','Sweden','China','Italy','Spain','Germany'};
deltaStart = [0.0075,0.01,0.001,0.022,0.003,0.001];
iGuessStrategy = [2,2,1,1,2,1];

thisCountry = countries{j};

[~,~,IG] = OneCountryODEAnalysis(MATdata,thisCountry);
if strcmp('Germany',thisCountry)
    [deaths , AIC , RLmatrix] = OneCountryIsolationAnalysis(MATdata,'Germany',IG,[],[1 1 2 2 1]);
else
    [deaths , AIC , RLmatrix] = OneCountryIsolationAnalysis(MATdata,thisCountry,IG);
end
[AIC2,deathRate] = OneCountryProjectionAnalysis_SIR(MATdata,thisCountry,deltaStart(j),iGuessStrategy(j));
AICvector = [AIC AIC2];
fullAIC = repmat(AICvector,length(AICvector),1) - repmat(AICvector',1,length(AICvector));
fullRLmatrix = sign(fullAIC).*exp(-abs(fullAIC)/2) + eye(length(AICvector));

allDeltas = [parameters.tryThisDeltaRange deathRate];
[~,iAIC] = min(AICvector);

%%
disp(' ');
disp('---------------------------------------');
disp([thisCountry,' death rate estimate ~ ',num2str(100*allDeltas(iAIC)),'%']);
disp('---------------------------------------');

