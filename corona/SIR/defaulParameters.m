function parameters = defaulParameters()

    parameters.limitDeaths = 250;

    parameters.isolateRate = 1;
    parameters.deisolateRate = 1;

    %the total time for a simulation
    parameters.totalTime = 60;
    %the default isolation time
    parameters.timeIso = 20;
    %the start time of isolation:
    parameters.timeInitial = 3;
    
    %when in isolation, the rate at which the isolate order is complied
    %with when moving into isolation
    parameters.positiveComplianceRate = 9;
    
    %2 and 3 are the control variables which are zero by default:
    %p(6) is the compliance rate, set to 0 to begin with: perfect compliance
    parameters.p = [3, 0.0, 0.0, 0.02, 0.5, 0, 0.5];

    parameters.plot = true;

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
    
    parameters.deadlabel = 'deaths (%S(0))';
    parameters.dailydeadlabel = 'daily deaths (%S(0) per diem)';
    
    %the observation time used for a feedback control policy:
    parameters.T = 0.5;
    
    parameters.OFFcontrolThresholds = [1e-1, 1e-2, 1e-3, 1e-4];

    parameters.tryThisDeltaRange = [0.005 0.0075 0.01 0.02 0.03];

    parameters.countryStrings = {'China', ...
        'United Kingdom', ...
        'Italy', ...
        'Germany', ...
        'Sweden', ...
        'Spain', ...
        'Denmark', ...
        'France', ...
        'Switzerland', ...
        'Iran',...
        'US','New York','California'};
    
    %the following is needed for OneCountryODEAnalysis:
    parameters.Pguesses = {[0.4399      0.43545       0.2806 ],...
        [0.41452      0.41421      0.24705 ],...
        [0.35112      0.3504      0.21804 ],...
        [0.32079      0.31977      0.16729 ],...
        [0.4157      0.41514      0.23485 ],...
        [0.41885      0.41781      0.25347],...        
        [0.49097 0.48447 0.21848 ],...
        [0.40009 0.40009 0.13918 ],...
        [0.23272 0.23168 0.012949 ],...
        [0.41438 0.41436 0.093819 ],...
        [0.4997 0.49969 0.26055] , [0.4997 0.49969 0.26055] , [0.4997 0.49969 0.26055]};
    
    parameters.LockdownStart = [1 , 61 , 46 , 57 , NaN , 51 , 48 , 54 , 57 , NaN , 57 , 57 , 56];
    
    parameters.movingAverageDays = 3;
    
    parameters.red = 0.75*[1 0.5 0.5];
    parameters.green = 0.75*[0.5 1 0.5];
    parameters.blue = 0.75*[0.5 0.5 1];
    parameters.grey = [1 1 1]*0.95;
    
end