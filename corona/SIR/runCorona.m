function [d,s] = runCorona(parameters,fig,sub)
    
    if nargin == 0
        parameters = defaulParameters();

        parameters.timeInitial = 1;
        parameters.timeIso = 5;

        parameters.positiveComplianceRate = 0;

        parameters.p(6) = parameters.positiveComplianceRate;
        parameters.plot = true;
        
    end

    if nargin <= 1
        if parameters.plot
            fig = figure;
            set(fig,'pos',[360 1 592 697])
        end
        sub.p = 3;
        sub.q = 1;
        sub.r = [1 2 3];
        sub.print = true;
        sub.label = 'exemplar case';
        sub.legend = true;
        sub.factor = 0.9;
        sub.colour = [0 0 0];
    else
        if parameters.plot
            figure(fig)
        end
    end
    
    isolateRate = parameters.isolateRate;
    deisolateRate = parameters.deisolateRate;
    timeInitial = parameters.timeInitial;
    timeIso = parameters.timeIso;
    totalTime = parameters.totalTime;
    grey = parameters.grey;

    p = parameters.p;

    timePostIso = totalTime - (timeIso+timeInitial);        

    IC = parameters.IC;
    S0 = IC(1);
    
    model = @(t,x)f(x,p);
    options = odeset('NonNegative',[1 1 1 1 1 1 1]);
    
    [noisot,noisoOutput] = ode45(model,[0,timeInitial],IC,options);
    ICisolate = noisoOutput(end,:);
    
    if timeIso > 0
        %ip = p(2); im = p(3);
        p(2) = isolateRate;
        p(3) = 0;
        isoModel = @(t,x)f(x,p);
        [isot,isoOutput] = ode45(isoModel,[timeInitial,timeIso+timeInitial],ICisolate,options);
    else
        isot = timeInitial;
        isoOutput = ICisolate;
    end

    ICnoisolate = isoOutput(end,:);
    %ip = p(2); im = p(3);
    p(2) = 0;
    p(3) = deisolateRate;
    noisoModel = @(t,x)f(x,p);
    [postt,postOutput] = ode45(noisoModel,[timeIso+timeInitial,...
        timeIso+timeInitial+timePostIso],ICnoisolate,options);
    
    fullt = [ noisot ; isot ; postt ];
    fulloutput = [noisoOutput ; isoOutput ; postOutput];
    
    S = fulloutput(1);
    I = fulloutput(2);
    R = fulloutput(3);
    Si = fulloutput(4);
    Ii = fulloutput(5);
    Ri = fulloutput(6);
    D = fulloutput(7);
    
    alldead = fulloutput(:,7);
    allinfecteds = fulloutput(:,2) + fulloutput(:,5);
    d = alldead(end);
    s = [num2str(100*alldead(end),4),'%'];
    
    if parameters.plot
        subplot(sub.p,sub.q,sub.r(1))
        rectangle('Position',[timeInitial 0 timeIso S0],'FaceColor',grey,'EdgeColor','none')
        hold on
        if parameters.timeIso > 0
            text(timeInitial+2,0.5,'isolation')
        else
            text(timeInitial+2,0.5,'no isolation period')
        end

        plot(fullt,fulloutput(:,1:end-1));
        axis tight
        if sub.legend
            legend(parameters.statesLegendText(1:6))
            legend('boxoff')
        end
        %xlim([ 1e-8 S0]);
        xlabel(parameters.xlabel)
        ylabel(parameters.ylabel)

        %plot([timeInitial timeInitial],[0 S0],'--r','color',grey);
        %plot([timeIso+timeInitial timeIso+timeInitial],[0 S0],'--k','color',grey);

        subplot(sub.p,sub.q,sub.r(2))
        plot(fullt,allinfecteds,'-r');
        hold on
        xlabel(parameters.xlabel)
        ylabel('infecteds')
        axis tight
        if nargin > 1
            ylim([0 S0])
        end
        legend(sub.label)
        legend('boxoff')

        if nargin <= 1
            subplot(sub.p,sub.q,sub.r(3))    
        else
            figure(2)
        end

        I = floor(length(fullt)*sub.factor);
        l = '-';

        if nargin <= 1
            text(fullt(I),0.9*alldead(I),s,'color',sub.colour)
            hold on
        end

        plot(fullt,alldead,l,'color',sub.colour)
        hold on
        xlabel(parameters.xlabel)
        ylabel(parameters.deadlabel)
        axis tight
        if nargin > 1
            ylim([0 0.04])
        end

        if sub.print
            my_export_fig(['latest_',num2str(p(6)),'-',num2str(timeInitial),'.pdf'])
        end
    end

end

