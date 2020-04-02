function perCapitaDownplotWExpFit(MATdata,includeList,logtrue)

    clc
    
    figure(1);
    set(1,'pos',[149   158   731   547])

    N = length(MATdata.country);

    cList = [];
    plots = [];
    
    if nargin < 2
        includeList = {};
    end
    if nargin < 3
        logtrue = false;
    end

    fitfun = @(p,t) p(3) * exp(-abs(p(1))*t) ./ (1 + abs(p(2))*exp(-abs(p(1))*t));
    extendTime = 0;
    adjR2limit = 0.6;
    
    for ctry = 1:N
        Ddata = sum(MATdata.deathData{ctry},1);        
        if (Ddata(end) > 50)
            s = (ctry-1)/(N-1);
            col = [s 0.5 1-s];
            F = find(Ddata,1,'first');
            x = Ddata(F:end);
            dx = diff(x);
            x = Ddata((F+1):end);
            if any(dx < 0)
                disp(['-ve ddata in ',MATdata.country{ctry}])
            else
                y = (dx./x);
                t = 0:length(x)-1;
                d = movmean(y,5);
                lm = fitlm(t,d);
                deltaD = lm.Coefficients.Estimate(2);
                if (lm.coefTest < 0.01 && lm.Coefficients.Estimate(2) < 0) || ...
                    ismember(MATdata.country{ctry},includeList)
                    l = '-';
                    lw = 1;
                    if ismember(MATdata.country{ctry},includeList)
                        lw = 4;
                        col = 'k';
                    end
                    if strcmp(MATdata.country{ctry},'China')
                        lw = 3;
                    end
                    
                    if logtrue
                        pl = semilogy(t,d,l,'linewidth',lw,'color',col);
                    else
                        pl = plot(t,d,l,'linewidth',lw,'color',col);
                    end
                    hold on
                    plot(t(end),d(end),'.','color',col,'markersize',50);
                    
                    try
                        expfit = fitnlm(t,d,fitfun,[0.1 0.1 0.1]);
                        if expfit.Rsquared.Adjusted > adjR2limit
                            extendt = [t (t(end) + (1:extendTime))];
                            EF = expfit.feval(extendt);
                            plot(extendt,EF,'--','linewidth',2,'color',col);
                            text(extendt(end),EF(end),MATdata.country{ctry});                            
                        else
                            text(t(end)+1,d(end),MATdata.country{ctry});
                        end
                    catch
                        text(t(end)+1,d(end),MATdata.country{ctry});
                    end
                    
                    cList = [cList ctry];
                    plots = [plots pl];
                    
                    MATdata.country{ctry} = [MATdata.country{ctry},...
                        ' (\Delta \approx ',num2str(deltaD,2),' d^{-1})'];
                    
                else
                    l = ':';
                end
            end
        end
    end
    
    axis tight
    legend(plots,MATdata.country{cList},'Location','northeast');
    legend('boxoff')

    xlabel('days from first recorded death')
    ylabel('deaths per day per capita')
    %title('Constant data above zero means exponential growth')
    
    if logtrue
        str = 'log';
    else
        str = '';
    end
    
    if isempty(includeList)
        my_export_fig([str,'perCapitaDeathDeclinesWexpFit.pdf'])
    else
        my_export_fig([str,'perCapitaDeathDeclinesWexpFit_IL.jpg'])
    end
end
