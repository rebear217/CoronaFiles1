
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
[~,~,IG_Italy] = OneCountryODEAnalysis(MATdata,'Italy');
close all
[~,~,IG_Spain] = OneCountryODEAnalysis(MATdata,'Spain');
close all
[~,~,IG_UK] = OneCountryODEAnalysis(MATdata,'United Kingdom');
close all
[~,~,IG_US] = OneCountryODEAnalysis(MATdata,'US');
close all
[~,~,IG_China] = OneCountryODEAnalysis(MATdata,'China');
close all
[~,~,IG_NY] = OneCountryODEAnalysis(USdata,'New York');
close all
[~,~,IG_CA] = OneCountryODEAnalysis(USdata,'California');
close all
[~,~,IG_Denmark] = OneCountryODEAnalysis(MATdata,'Denmark');
close all
[~,~,IG_Germany] = OneCountryODEAnalysis(MATdata,'Germany');
close all
[~,~,IG_France] = OneCountryODEAnalysis(MATdata,'France');
close all
[~,~,IG_Switzerland] = OneCountryODEAnalysis(MATdata,'Switzerland');
close all
[~,~,IG_Sweden] = OneCountryODEAnalysis(MATdata,'Sweden');
close all

Italydeaths = OneCountryIsolationAnalysis(MATdata,'Italy',IG_Italy);
close all
Spaindeaths = OneCountryIsolationAnalysis(MATdata,'Spain',IG_Spain);
close all
UKdeaths = OneCountryIsolationAnalysis(MATdata,'United Kingdom',IG_UK);
close all
USdeaths = OneCountryIsolationAnalysis(MATdata,'US',IG_US,[],[2 2 2 2 1]);
close all
Chinadeaths = OneCountryIsolationAnalysis(MATdata,'China',IG_China);
close all
%CA & NY needs wierd initial guess strategies for each fit seperately:
NYdeaths = OneCountryIsolationAnalysis(USdata,'New York',IG_NY,57,[1 2 1 1 1]);
close all
CAdeaths = OneCountryIsolationAnalysis(USdata,'California',IG_CA,56,[1 2 2 2 1]);
close all
Denmarkdeaths = OneCountryIsolationAnalysis(MATdata,'Denmark',IG_Denmark);
close all
Germanydeaths = OneCountryIsolationAnalysis(MATdata,'Germany',IG_Germany,[],[1 1 2 2 1]);
close all
Francedeaths = OneCountryIsolationAnalysis(MATdata,'France',IG_France);
close all
Switzerlanddeaths = OneCountryIsolationAnalysis(MATdata,'Switzerland',IG_Switzerland);
close all
Swedendeaths = OneCountryIsolationAnalysis(MATdata,'Sweden',IG_Sweden);

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
    ODEestimate99 UKdeaths];

%%

disp(['Estimates : ',num2str(UKestimateVector,5)]);
disp(['Mean : ',num2str(mean(UKestimateVector),5)]);
disp(['95% CI: ',num2str(mean(UKestimateVector) + 2*std(UKestimateVector)/sqrt(length(UKestimateVector)),5)]);

%%

disp('')
disp('')

disp('------------')
disp('Finished...')
disp('------------')

