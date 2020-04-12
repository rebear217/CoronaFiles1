function my_export_fig(filename,forceFilename)

if nargin < 2
        forceFilename = 0;
    end
    
    if not(exist('figures','dir'))
        mkdir('figures');
    end
    if not(exist('./figures/countryODEanalysis','dir'))
        mkdir('./figures/countryODEanalysis');
    end
        
    if forceFilename
        Filename = filename;
    else
        Filename = ['./figures/',filename];
    end

    Filename = Filename(~isspace(Filename));
    
    try
        export_fig(Filename);
    catch
        disp(['Cannot export_fig ',filename',' to a figure file ... trying print']);
        print(Filename, '.pdf')
        disp('See https://uk.mathworks.com/matlabcentral/fileexchange/23629-export_fig');
        %web('https://uk.mathworks.com/matlabcentral/fileexchange/23629-export_fig')
    end
    
end

