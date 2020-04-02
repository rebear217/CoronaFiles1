function PikeAnalysis(MATdata)

    close all
    figure(1)
    set(1,'pos',[60     1   885   704])
    
    parameters = defaulParameters();
    countryStrings = parameters.countryStrings;
    extendtime = 35;
    
    for ctry = 1:9

        countryStr = countryStrings{ctry};
        
        disp(countryStr);
        cj = find(strcmp(countryStr,MATdata.country));
        cdata = MATdata.deathData{cj};
        scdata = sum(cdata,1);
        firstDeath = find(scdata,1,'first');

        dataToFit = scdata(firstDeath:end);
        deaths = dataToFit(end);
        scdataToFit = dataToFit / deaths;
        
        T = length(dataToFit);
        times = 1:T;
        
        Pguesses = {[0.44775 0.28533 0.44341],...
            [0.49097 0.21848 0.48447],...
            [0.40009 0.13918 0.40009],...
            [0.31967 0.03864 0.31884],...
            [0.36294 0.17064 0.36257],...
            [0.23272 0.012949 0.23168],...
            [0.41816 0.162 0.41768],...
            [0.41445 0.10843 0.4143],...
            [0.41438 0.093819 0.41436],...
            };
        
        opts = statset('MaxIter',1000,'TolX',1e-10,'TolFun',1e-10);
        
        pg = Pguesses{ctry};
        %pg = [1 1 1]*0.5;
        
        fit = fitnlm(times,scdataToFit,@(p,t)Solve(p,t),pg,...
            'Options',opts);
        betterGuess = fit.Coefficients.Estimate;

        figure(1)
        subplot(3,3,ctry)
        extendtimes = [times times(end)+(1:extendtime)];
        Y = deaths*fit.feval(extendtimes);

        %policy model error bars:    
        [beta,resid,J,sigma] = nlinfit(times,scdataToFit,@(p,t)Solve(p,t),betterGuess);
        [deltaFit, delta] = nlpredci(@(p,t)Solve(p,t),extendtimes,beta,resid,'Covar',sigma);

        upperCI = deltaFit + delta;
        lowerCI = deltaFit - delta;
        midCI = deltaFit;
        
        ylimit = max([1.5*Y dataToFit]);
        
        rectangle('position',[times(end) 0 extendtime ylimit],...
            'facecolor',[1 1 1]*0.8,'edgecolor','none');
        hold on
        
        plot(extendtimes,Y,'-k','linewidth',1);
        %plot(times,dataToFit,'o','markersize',6,'linewidth',1);
        plot(times,dataToFit,'.','markersize',26)        
        axis tight
        ylim([0 ylimit])
        
        legend({['SIR fit (adj R^2 \approx ',num2str(fit.Rsquared.Adjusted,4),' )'],...
            [countryStrings{ctry},' data']},'Location','northwest')
        legend('boxoff')

        %plot(extendtimes,deaths*midCI,'-k','linewidth',1)
        %plot(extendtimes,deaths*lowerCI,':k')
        plot(extendtimes,deaths*upperCI,':k')
        
        if ismember(ctry , [1,4,7])
            ylabel('cumulative deaths')
        end
        if ismember(ctry , [7,8,9])
            xlabel('days from first recorded case')
        end
        
        if ctry == 2
            text(10,ylimit*0.5,'$\frac{dD}{dt} = a - bD - \gamma e^{-D}$',...
                'interpreter','latex');
        end
        if ctry == 1
            text(times(end) + 1,ylimit*0.2,{[num2str(extendtime),' day'],'projection'});
        end
    end

    my_export_fig('PikeODEAnalysis.pdf')
    
end
    
function solution = Solve(p,T)

    IC = 0;
    model = @(t,x)F(x,p);
    options = odeset('NonNegative',[1]);
    [soltimes,output] = ode113(model,[1,T(end)],IC,options);
    solution = interp1(soltimes,output,T,'pchip');
    
end

function Ddot = F(x,p)

    A = abs(p(1));
    B = abs(p(2));
    C = abs(p(3));

    D = x(1);    
    %Ddot = d - (rho + d)*D - S0*d*exp(-lam*D/d);
    Ddot = A - B*D - C*exp(-abs(D));

end