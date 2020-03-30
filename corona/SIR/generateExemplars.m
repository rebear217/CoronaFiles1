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

runCorona(parameters)
my_export_fig('basicModelfullCompliance.pdf')

%%

parameters.timeIso = 0;
runCorona(parameters)

my_export_fig('basicModelNoIsolation.pdf')
