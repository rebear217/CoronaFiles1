function stupidity(MATdata,extendTime,country)

    if nargin < 2
        extendTime = 10;
    end

	if nargin < 3
        country = 'United Kingdom';
    end

    F = find(strcmp(country,MATdata.country));
    data = sum(MATdata.deathData{F},1);
    f = find(data,1,'first');
    data = data(f:end);
    X = 1:length(data);

    fit = fitnlm(X,data,@(p,t)(p(1)*exp(p(2)*t)),[1 0.1])

    xX = [X (X(end)+(1:extendTime))];

    Y = fit.feval(xX);
    mx = max([data Y]);
    
    %rectangle('Position',[X(end) 0 extendTime mx],...
    %    'facecolor',0.7*[1 1 1],'edgecolor','none');
    %hold on
    semilogy(X,data,'ok');
    hold on
    plot(xX,20000*ones(size(xX)),'-b');
    plot(xX,fit.feval(xX),'-k');
    xlabel('days since 1st death')
    ylabel('cumulative deaths')
    axis tight
    YL = ylim;
    text(X(end)+1.5,0.95*YL(2),'projection')
    legend(country,'location','northwest');
    
    Tcrit = log(20000/fit.Coefficients.Estimate(1))/fit.Coefficients.Estimate(2);
    plot([Tcrit Tcrit],[0 YL(2)],':r');
    T = Tcrit - X(end);
    text(X(end)-1.5,0.2*YL(2),['govt estimate in ',num2str(T,2),'days']);

    my_export_fig('stupid.png')
end

