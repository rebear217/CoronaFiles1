UKURL = 'https://raw.githubusercontent.com/tomwhite/covid-19-uk-data/master/data/covid-19-totals-uk.csv';
options = weboptions('RequestMethod','get','ArrayFormat','csv','ContentType','text');
CSVfilename = 'UKOnlyData.CSV';
UKoutfilename = websave(CSVfilename,UKURL,options);

%%

UKdataMatrix = csvread(UKoutfilename,1,1);

%%

%because this data starts on 25/1 we need to pad with
%zeros to be consistent:

Z = zeros(2,3);
finalUKdataMatrix = [Z ; UKdataMatrix];

%%

jUK = find(strcmp(MATdata.country,'United Kingdom'));
MATdata.deathData{jUK} = finalUKdataMatrix(:,3)';
MATdata.ICaseData{jUK} = finalUKdataMatrix(:,2)';
%the following number of tests data is not being used:
%MATdata.ITestData{jUK} = finalUKdataMatrix(:,1)';

