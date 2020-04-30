function [PCD,dx,x] = getPCD(Ddata,country,DisplayOff)

    if nargin < 3
        DisplayOff = false;
    end

    parameters = defaulParameters();

    if nargin < 2
        country = 'unspecified country';
    end
    
    F = find(Ddata,1,'first');
    x = Ddata(F:end);
    dx = diff(x);
    
    if any(dx < 0) || any(x < 0)
        %remove possible -ve cumulative death data:
        Ddata = Ddata.*(Ddata >= 0) + 0*(Ddata < 0);
        x = Ddata(F:end);
        dx = diff(x);
        
        %remove possible -ve daily death data:
        dx = dx.*(dx >= 0) + 0*(dx < 0);
        x = cumsum(dx);
        
        if ~DisplayOff
            disp(['Assuming 0 for -ve ddata in ',country])
        end
    end
    
    x = Ddata((F+1):end);
    PCD = movmean(dx./x,parameters.movingAverageDays);                

end