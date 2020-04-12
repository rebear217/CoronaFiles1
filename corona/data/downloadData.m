close all
clc
clear variables

%%

%there are inconsistencies in French data: so this
%is recovered by hand from https://dashboard.covid19.data.gouv.fr
%on Thur 9 April
urlFranceData = flipud([7632;7091;6494;5889;5532;5091;4503;4032;3523;3024;2606;...
2314;1995;1696;1331;1100;860;674;562;372;264;175;148;127;91;79;61;48;33;25;19;16;9;7;...
5;3;1])';
%this is merged with the JH data below
%I am still not sure what the real data is having done this...

%China data head has been taken from ECDC and sent to me by Tom Pike

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
%M = 66;%worked up to 29 March
%M = 67;%worked up to 30 March
%M = 68;%worked up to 30 March
%M = 69;%worked up to 1 April
%M = 70;%worked up to 2 April
%M = 71;%worked up to 3 April
%M = 72;%worked up to 4 April
%M = 73;%worked up to 5 April
%M = 74;%worked up to 6 April
%M = 75;%worked up to 7 April
%M = 76;%worked up to 8 April
%M = 77;%worked up to 9 April
%M = 78;%worked up to 10 April
M = 79;%worked up to

allCountries = {data{:,2}};
countries = unique(allCountries);

MATdata.country = {};
MATdata.deathData = {};
MATdata.M = M;

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
    
    if strcmp(thisCountry,'France')
        FranceData = sum(MATdata.deathData{c},1);        
        FranceDataTail = FranceData(74:end);
        MATdata.deathData{c} = [urlFranceData FranceDataTail];        
    end
    
	%if strcmp(thisCountry,'China')
    %    ChinaHead = cumsum([1 0 0 0 1 0 0 0 1 0 3]);
    %	ChinaData = sum(MATdata.deathData{c},1);        
    %    MATdata.deathData{c} = [ChinaHead ChinaHead(end) + ChinaData];
    %end
    
    figure(1)
    semilogy(sum(rowData,1));
    hold on
    xlabel('day')
    ylabel('deaths')
    
end

save('deathData.mat','MATdata')
axis tight
export_fig('basicDeathDataPlot.PDF')

%%

clear variables
load('deathData.mat')

