function parameters = defaulParameters()

    parameters.isolateRate = 1;
    parameters.deisolateRate = 1;

    %the total time for a simulation
    parameters.totalTime = 50;
    %the default isolation time
    parameters.timeIso = 15;
    %the start time of isolation:
    parameters.timeInitial = 1;
    
    %when in isolation, the rate at which the isolate order is complied
    %with when moving into isolation
    parameters.positiveComplianceRate = 9;
    
    %2 and 3 are the control variables which are zero by default:
    %p(6) is the compliance rate, set to 0 to begin with: perfect compliance
    parameters.p = [3, 0.0, 0.0, 0.02, 0.5, 0, 0.5];

    parameters.plot = true;

    %    lam = p(1);
    %    ip = p(2);
    %    im = p(3);
    %    d = p(4);
    %    rho = p(5);
    %    lamE = p(6);
    %    rhoI = p(7);

    parameters.grey = [1 1 1]*0.9;

    %standard initial conditions:
    I0 = 1e-3;
    S0 = 1-I0;
    R0 = 0;
    Si0 = 0;
    Ii0 = 0;
    Ri0 = 0;
    D0 = 0;
    
    IC = [S0,I0,R0,Si0,Ii0,Ri0,D0];    
    parameters.IC = IC;
    parameters.statesLegendText = {'susceptible','infected',...
        'recovered','i-susceptible',...
        'i-infected','i-recovered','dead'};
    
    parameters.ylabel = 'freq';
    parameters.xlabel = 'time';
    
    parameters.deadlabel = 'dead (as %S(0)';
    
    %the observation time used for a feedback control policy:
    parameters.T = 0.5;
    
    parameters.OFFcontrolThresholds = [1e-1, 1e-2, 1e-3, 1e-4];
    
end