clc
close all
clear variables

load('deathData.mat')

%%

%csv source
%https://github.com/CSSEGISandData/COVID-19/blob/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_US.csv

URL = 'https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_US.csv';

options = weboptions('RequestMethod','get','ArrayFormat','csv','ContentType','text');
CSVfilename = 'US-jhData.CSV';

outfilename = websave(CSVfilename,URL,options);

%%

%get the data headers by importing text:
data = USimportData(outfilename);
%get the data itself:
dataMatrix = csvread(outfilename,2,14);

%%

allStates = {data{:,7}};
states = unique(allStates);

[n,m] = size(data);
[N,M] = size(dataMatrix);

USdata.country = {};
USdata.deathData = {};
USdata.M = M;

%leave off any trailing zeros:
%M = M - 1;

ss = 1;
for s = 1:length(states)
    try
        thisState = states{s};
        F = find(strcmp(thisState,allStates));

        rowData = zeros(length(F),M);
        for j = 1:length(F)
            rowData(j,:) = rowData(j,:) + dataMatrix(F(j),1:M);
        end

        %data has had -ve entries, so assume they are zero:
        rowData = rowData.*(rowData >= 0) + 0.*(rowData < 0);

        USdata.country{ss} = thisState;
        USdata.deathData{ss} = rowData;

        figure(2)
        semilogy(sum(rowData,1));
        hold on
        xlabel('day')
        ylabel('deaths')
        axis tight
        ss = ss + 1;
    catch
        disp([thisState,' data caused an error - ignored']);        
    end    
end

F = find(strcmp(MATdata.country,'China'));
USdata.country{ss} = 'China';
%Just in case China has a trailing zero, remove it:
%USdata.deathData{ss} = MATdata.deathData{F}(:,1:end-1);
USdata.deathData{ss} = MATdata.deathData{F};

save('USdeathData.mat','USdata')
axis tight
export_fig('USbasicDeathDataPlot.PDF')

%%

clear variables
load('USdeathData.mat')
load('deathData.mat')
