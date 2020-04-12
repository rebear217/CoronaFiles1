function generateExemplars()

    clear vars
    close all

    %%

    parameters = defaulParameters();
    parameters.plot = true;

    parameters.timeInitial = 2;
    parameters.timeIso = 5;
    parameters.totalTime = 20;

    parameters.positiveComplianceRate = 0;
    parameters.p(6) = parameters.positiveComplianceRate;

    %%

    d = runCorona(parameters);
    disp(['Check: dead with isolation as a % of S(0) is ',num2str(100*d)]);

    my_export_fig('basicModelfullCompliance.pdf')

    %%

    parameters.timeIso = 0;
    d = runCorona(parameters);
    disp(['Check: dead with no isolation as a % of S(0) is ',num2str(100*d)]);

    my_export_fig('basicModelNoIsolation.pdf')

end