function [DSDestimate,DSDestimate99,latexStr] = doSomeDataFitting()

    close all
    load('../data/deathData.mat')

    [outsideFilename,DSDestimate,DSDestimate99,latexStr] = fitToCountryData(MATdata);

    %%

    figure(1)
    my_export_fig(outsideFilename)

    figure(2)
    my_export_fig('datafitPanel.pdf')

    figure(3)
    my_export_fig('datafitPanelDiff.pdf')

end