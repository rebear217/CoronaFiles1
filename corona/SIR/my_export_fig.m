function my_export_fig(filename)

    if not(exist('figures','dir'))
        mkdir('figures');
    end

    try
        export_fig(['./figures/',filename])
    catch
        disp(['Cannot export_fig ',filename',' to a figure file ... trying print']);
        print(['./figures/',filename], '.pdf')
        disp('See https://uk.mathworks.com/matlabcentral/fileexchange/23629-export_fig');
        %web('https://uk.mathworks.com/matlabcentral/fileexchange/23629-export_fig')
    end
    
end

