function generateHeatmap()

    for trial = 1:2
        clear vars

        parameters = defaulParameters();

        parameters.timeIso = 0;
        parameters.plot = false;

        dnopolicy = runCorona(parameters);

        if trial == 2
            parameters.p(6) = parameters.positiveComplianceRate;
        end

        %%

        imax = 30;
        jmax = 20;
        isolationDuration = (1:jmax);
        isolationStartTimes = (1:imax)/2;

        D = zeros(imax,jmax);

        for i = 1:imax
            for j = 1:jmax
                parameters.timeInitial = isolationStartTimes(i);
                parameters.timeIso = isolationDuration(j);

                [D(i,j),s] = runCorona(parameters);
            end
        end

        relD = 100*(1 - D / dnopolicy);

        %%

        imagesc(relD);
        colormap(hot)

        colorbar
        xlabel('isolation duration')
        ylabel('isolation start time')

        set(gca,'Ytick',2:2:imax)
        set(gca,'Xtick',1:jmax)

        set(gca,'YtickLabel',isolationStartTimes(2:2:imax))
        set(gca,'XtickLabel',isolationDuration)

        %%

        title('% gain from isolation policy')

        if trial == 1
            my_export_fig('gainHeatMap_completeIso.png')
        else
            my_export_fig('gainHeatMap_partialIso.png')
        end

    end

end