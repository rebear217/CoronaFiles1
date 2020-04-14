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

%% get death data:

%csv source
%https://github.com/CSSEGISandData/COVID-19/blob/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_global.csv

URL = 'https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_global.csv';
options = weboptions('RequestMethod','get','ArrayFormat','csv','ContentType','text');
CSVfilename = 'jhData.CSV';
outfilename = websave(CSVfilename,URL,options);

%get the data headers by importing text:
data = importData(outfilename);
%get the data itself:
dataMatrix = csvread(outfilename,1,5);

%% get infected case data:

URLI = 'https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv';
optionsI = weboptions('RequestMethod','get','ArrayFormat','csv','ContentType','text');
CSVfilenameI = 'jhDataIcases.CSV';

outfilenameI = websave(CSVfilenameI,URLI,optionsI);

IcaseCountries = importICaseHeadings(outfilenameI);
IcaseDataMatrix = csvread(outfilenameI,1,5);

%%

clc
close all

[n,m] = size(data);
[N,M] = size(dataMatrix);

%JH puts a trailing 0 onto data in readiness
%for the next day so we need this:

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
%M = 79;%worked up to 11 April
%M = 80;%worked up to 12 April
%M = 81;%worked up to 13 April
M = 82;%worked up to 

%or this:
%M = M-1;

allCountries = {data{:,2}};
countries = unique(allCountries);

MATdata.country = {};
MATdata.deathData = {};
MATdata.M = M;

figure(1)
set(1,'pos', [148     1   958   704]);

for c = 1:length(countries)
    s = (c-1) / (length(countries) - 1);
    col = [s 0.5 1-s];
    
    thisCountry = countries{c};
    F = find(strcmp(thisCountry,allCountries));
    Fi = find(strcmp(thisCountry,IcaseCountries));
    
    rowData = zeros(length(F),M);
    rowDataICase = zeros(length(Fi),M);
    for j = 1:length(F)
    	rowData(j,:) = rowData(j,:) + dataMatrix(F(j),1:M);
    	rowDataICase(j,:) = rowDataICase(j,:) + IcaseDataMatrix(Fi(j),1:M);
    end
    
    %data has had -ve entries, so assume they are zero:
    rowData = rowData.*(rowData >= 0) + 0.*(rowData < 0);
    rowDataICase = rowDataICase.*(rowDataICase >= 0) + 0.*(rowDataICase < 0);

    MATdata.country{c} = thisCountry;
    MATdata.deathData{c} = rowData;
    MATdata.ICaseData{c} = rowDataICase;
    
    %if strcmp(thisCountry,'France')
    %    FranceData = sum(MATdata.deathData{c},1);        
    %    FranceDataTail = FranceData(74:end);
    %    MATdata.deathData{c} = [urlFranceData FranceDataTail];        
    %end
    
	%if strcmp(thisCountry,'China')
    %    ChinaHead = cumsum([1 0 0 0 1 0 0 0 1 0 3]);
    %	ChinaData = sum(MATdata.deathData{c},1);        
    %    MATdata.deathData{c} = [ChinaHead ChinaHead(end) + ChinaData];
    %end
    
    subplot(2,2,1)
    semilogy(sum(rowData,1));
    xlabel('day')
    ylabel('deaths')
    hold on

    subplot(2,2,2)
    semilogy(sum(rowDataICase,1));
    xlabel('day')
    ylabel('I cases')
    hold on
    
    subplot(2,2,3)
    loglog(sum(rowDataICase,1),sum(rowData,1));
    xlabel('I cases')
    ylabel('deaths')
    hold on
    
    subplot(2,2,4)
    semilogx(rowData(:),100*rowData(:)./rowDataICase(:),'.','color',col,'markersize',18);
    ylabel('+ve SARS CoV-2 tests per death (%)')
    xlabel('deaths')
    set(gca,'Xtick',[1 10 100 1000 10000]);
    hold on
    
end

save('deathData.mat','MATdata')
axis tight

subplot(2,2,1);
axis tight
subplot(2,2,2);
axis tight
subplot(2,2,3);
axis tight
subplot(2,2,4);
axis tight

export_fig('basicDeathDataPlot.PDF')

%%

clear variables
load('deathData.mat')
load('USdeathData.mat')
