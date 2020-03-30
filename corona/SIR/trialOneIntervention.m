close all
clear vars
clc

%%

controlTestTime = 200;

parameters = defaulParameters();
    parameters.totalTime = controlTestTime;
    thresholds = parameters.OFFcontrolThresholds;

[basestate,baseoutput,basetimes] = moveForwardOne(parameters,[],[],parameters.totalTime);
dnopolicy = basestate(end,7);

figure(1)
semilogx(basetimes,100*baseoutput(:,7),'-k','linewidth',3)
xlim([1 controlTestTime])
hold on

%%

controllegends = {'baseline deaths'};

for j = 1:length(thresholds)
    
    th = thresholds(j);
    controllegends{j+1} = num2str(th);

    %for cleanliness, reset parameters:
    parameters = defaulParameters();
        parameters.totalTime = controlTestTime;

    %there is no point simulating non-compliance, it seems
    %clear that it makes matters worse:
    %parameters.p(6) = parameters.positiveComplianceRate;

    isolateRate = parameters.isolateRate;
    deisolateRate = parameters.deisolateRate;

    dT = parameters.T;
    observationPoints = parameters.totalTime / dT;

    [outstate,allstates,times] = moveForwardOne();
        %infecteds = outstate(:,2);% + outstate(:,5);
        infecteds = allstates(:,2);
        di = diff(infecteds) ./ diff(times);
        d2i = diff(di) ./ diff(times(1:end-1));
        observeON = mean(d2i)/mean(infecteds);
        observeOFF = mean(di)/mean(infecteds);

    startedIsolation = false;
    stoppedIsolation = false;

    control.plusPolicy = 0;
    control.minusPolicy = 0;
    control.useTheControl = true;
    control.forceTheControlON = false;
    control.forceTheControlONforever = false;

    allEpistates = zeros(2+observationPoints,7);
    allEpitimes = zeros(2+observationPoints,1);

    policyOntimes = zeros(2+observationPoints,1);

    %initial conditions:
    allEpitimes(1) = 0;
    allEpistates(1,:) = parameters.IC';
    S0 = parameters.IC(1);

    %epidemic begins unknown to policy makers: controls cannot yet start
    allEpitimes(2) = dT;
    allEpistates(2,:) = outstate';

    %now assume epidemic has been detected by policy makers:

    %%

    for timepoint = 1:observationPoints
        thisTime = (1+timepoint)*dT;

        [outstate,allstates,times] = moveForwardOne(parameters,outstate,control);

        previousObserveON = observeON;
        previousObserveOFF = observeOFF;

        %infecteds = outstate(:,2);% + outstate(:,5);
        infecteds = allstates(:,2);
        dead = allstates(:,7);

        DT = diff(times);
        DT2 = diff(DT);

        di = interp1( times(1:end-1), diff(infecteds) ./ DT , times , 'linear' , 'extrap');    
        d2i = interp1( times(1:end-1) , diff(di) ./ DT , times , 'linear' , 'extrap' );

        dd = interp1( times(1:end-1) , diff(dead) ./ DT , times , 'linear' , 'extrap' );
        d2d = interp1( times(1:end-1) , diff(dd) ./ DT , times , 'linear' , 'extrap' );
        d3d = interp1( times(1:end-1), diff(d2d) ./ DT , times , 'linear' , 'extrap' );

        observeON = mean(d2i)/mean(infecteds);
        nextObserveONGuess = -previousObserveON + 2*observeON;

        observeOFF = mean(dd)/mean(dead) < th*S0;
        %observeOFF = mean(d2d)/mean(dead) < -0.2;

        %disp([mean(dd),mean(d2d)/mean(dead)])

        if control.useTheControl
            if (previousObserveON*observeON < 0)
            %if (nextObserveONGuess*observeON < 0) || (previousObserveON*observeON < 0)
                if not(startedIsolation)
                    startedIsolation = true;
                    policyStartTime = thisTime;
                end
            end

            if observeOFF
                if startedIsolation
                    stoppedIsolation = true;
                    policyEndTime = thisTime;                
                end
            end

            if startedIsolation && not(stoppedIsolation)
                control.plusPolicy = isolateRate;
                control.minusPolicy = 0;
            end

            if startedIsolation && stoppedIsolation            
                control.plusPolicy = 0;
                control.minusPolicy = deisolateRate;
            end

        end

        if control.forceTheControlON 
            policyStartTime = 2*dT;
            if (thisTime < parameters.timeIso)
                control.plusPolicy = isolateRate;
                control.minusPolicy = 0;
            else
                control.plusPolicy = 0;
                control.minusPolicy = deisolateRate;        
            end
        end

        if control.forceTheControlONforever
            policyStartTime = 2*dT;
            control.plusPolicy = isolateRate;
            control.minusPolicy = 0;
        end

        policyOntimes(2+timepoint) = (control.plusPolicy > 0);

        allEpitimes(2+timepoint) = thisTime;
        allEpistates(2+timepoint,:) = outstate';

        %figure(2)
        %plot(thisTime,observeON,'ok');
        %hold on

    end

    %%

    D = outstate(7);

    relD = 100*(1 - D / dnopolicy);

    figure(1)
    
    cp = (j-1)/(length(thresholds)-1);
    colour = cp*[1 0 0] + (1-cp)*[0 0 1];

    deathNumbers = 100*allEpistates(:,7);
    plot(allEpitimes,deathNumbers,'-','color',colour)
    save(['./synthMatFiles/syntheticDeaths_',num2str(j)],'allEpitimes','deathNumbers');
    
    axis tight
    xlim([1 controlTestTime])
    xlabel(parameters.xlabel)
    ylabel(parameters.deadlabel)

    figure(2)
    set(2,'pos',[29 1 1117 704])
    subplot(length(thresholds)/2,2,j)
    
    for timepoint = 1:observationPoints
        if policyOntimes(2+timepoint)
            thisTime = (1+timepoint)*dT;
            rectangle('Position',[thisTime 0 dT S0],'FaceColor',parameters.grey,'EdgeColor','none')
            hold on
        end
    end
    
    plot(allEpitimes,allEpistates,'-')
    axis tight
    xlim([0 30])
    parameters.statesLegendText{7} = ['dead: ',num2str(relD,3),'% reduction'];
        legend(parameters.statesLegendText)
        legend('boxoff')
    xlabel(parameters.xlabel)
    ylabel(parameters.ylabel)
    
end

figure(1)
legend(controllegends,'Location','southeast');
legend('boxoff')

%%

my_export_fig('feedbackControlComparisons.pdf')

%%

figure(2)
text(policyStartTime+1,0.9*S0,'isolation')

my_export_fig('feedbackControlDeaths.pdf')

