clc
close all
clear variables

load('../data/deathData.mat')
outsideFilename = fitToCountryData(MATdata);

%%

figure(1)
my_export_fig(outsideFilename)

figure(2)
my_export_fig('datafitPanel.pdf')
