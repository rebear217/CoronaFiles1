function generatePlots()

    clear vars
    close all

    %%

    parameters = defaulParameters();

    %%
    fig = figure(1);
    set(fig,'pos',[36 1 1191 704])

    s.p = 5;%leave this
    s.q = 2;%leave this
    s.r = [1 2];
    s.print = false;

    col = @(s)0.75*((1-s)*[1 0.5 0.1] + s*[0.5 0.5 1]);

    legends = {'early intervention imperfect isolation',....
        'late intervention imperfect isolation',...
        'early intervention complete isolation',...
        'late intervention complete isolation',...
        'no isolation'};

    fulllegends = {'','','','',''};

    d = [0 0 0 0];
    t = {'','','',''};

    %%

    parameters.timeInitial = 1;
    parameters.p(6) = parameters.positiveComplianceRate;
    s.r = [1 2];
    s.legend = true;
    s.factor = 0.6;
    s.colour = col(0.0);
    s.label = legends{1};

    [d(1),t{1}] = runCorona(parameters,fig,s);

    %%

    parameters.timeInitial = 3;
    parameters.p(6) = parameters.positiveComplianceRate;
    s.r = [3 4];
    s.legend = false;
    s.factor = 0.4;
    s.colour = col(0.33);
    s.label = legends{2};

    [d(2),t{2}] = runCorona(parameters,fig,s);

    %%

    s.r = [5 6];
    parameters.timeInitial = 1;
    parameters.p(6) = 0;
    s.legend = false;
    s.factor = 0.5;
    s.colour = col(0.66);
    s.label = legends{3};

    [d(3),t{3}] = runCorona(parameters,fig,s);

    %%

    s.r = [7 8];
    parameters.timeInitial = 3;
    parameters.p(6) = 0;
    s.legend = false;
    s.factor = 0.3;
    s.colour = col(1);
    s.label = legends{4};

    [d(4),t{4}] = runCorona(parameters,fig,s);

    %%

    parameters.timeIso = 0;

    s.r = [9 10];
    s.factor = 1.0;
    s.colour = [0 0 0];
    s.label = legends{5};

    [d(5),t{5}] = runCorona(parameters,fig,s);

    %%

    for j = 1:5

        if j < 5
            lg = [num2str(100*d(j)/d(5),4),'%'];
            fulllegends{j} = [legends{j},' - ',lg];
        end

        if j == 5
            fulllegends{j} = [legends{j},' - ',num2str(d(5),4)];
        end
    end

    %%

    figure(1)
    my_export_fig('allDynamics.pdf')

    figure(2)
    legend(fulllegends,'Location','east')
    legend('boxoff')
    my_export_fig('dead.pdf')

end