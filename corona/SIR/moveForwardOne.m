function [outstate,fulloutput,times] = moveForwardOne(parameters,instate,control,forceTime)

    if nargin == 0
        parameters = defaulParameters();
    end
    
    if nargin < 2 || isempty(instate)
        instate = parameters.IC;
    end
    
    if nargin < 3 || isempty(control)
        control.plusPolicy = 0;
        control.minusPolicy = 0;
    end
    
    if nargin < 4
        T = parameters.T;
    else
        T = forceTime;
    end

    p = parameters.p;
    
    %controls:
    p(2) = control.plusPolicy;
    p(3) = control.minusPolicy;
    
    model = @(t,x)f(x,p);
    options = odeset('NonNegative',[1 1 1 1 1 1 1]);

    [times,fulloutput] = ode45(model,[0,T],instate,options);
    
    outstate = fulloutput(end,:);
    
end