close all
clc
clear variables

%%

%csv source
%https://github.com/CSSEGISandData/COVID-19/blob/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_global.csv

URL = 'https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_global.csv';
options = weboptions('RequestMethod','get','ArrayFormat','csv','ContentType','text');
CSVfilename = ['jhData-',date,'.CSV'];

outfilename = websave(CSVfilename,URL,options);

%%

%get the data headers by importing text:
data = importData(outfilename);
%get the data itself:
dataMatrix = csvread(outfilename,2,5);

%data = importData('time_series_covid19_deaths_global.csv');

%%

clc
close all

[n,m] = size(data);
[N,M] = size(dataMatrix);

%sometimes data is missing from some countries
%for the most recent times - therefore are zero - so remove them:
%M = M - 1;
%M = 64;%worked up to 27 March
%M = 65;%worked up to 28 March
M = 66;%worked up to 29 March

allCountries = {data{:,2}};
countries = unique(allCountries);

MATdata.country = {};
MATdata.deathData = {};

for c = 1:length(countries)
    
    thisCountry = countries{c};
    F = find(strcmp(thisCountry,allCountries));
    rowData = zeros(length(F),M);
    for j = 1:length(F)
    	rowData(j,:) = rowData(j,:) + dataMatrix(F(j),1:M);
    end
    
    %data has had -ve entries, so assume they are zero:
    rowData = rowData.*(rowData >= 0) + 0.*(rowData < 0);

    MATdata.country{c} = thisCountry;
    MATdata.deathData{c} = rowData;
    
    figure(2)
    semilogy(sum(rowData,1));
    hold on
    xlabel('day')
    ylabel('deaths')
    
end

save('deathData.mat','MATdata','M')
axis tight
export_fig('basicDeathDataPlot.PDF')

%%

figure(1)
close all
loglogPlot(MATdata);
export_fig('basicLogLogPlot.PDF')

%%

clear variables
load('deathData.mat')

