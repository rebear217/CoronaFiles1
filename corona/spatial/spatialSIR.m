function spatialSIR(T)

    close all
    clc
    figure(1)
    set(1,'pos',[96           1        1034         704])
    
    if nargin < 1
        T = 0.000001;
    end
    
    %%
    
    p = 1;
    rho = 1;
    lambda = 0.00001;
    dr = 0.001;
    sigma = 0.001;
    
    %%
    
    Nx = 10;
    x = 0:(1/(Nx-1)):1;
    
    D2x = diag(-2*ones(Nx,1)) + ...
        diag(ones(Nx-1,1),1) + diag(ones(Nx-1,1),-1);
    D2x(2,1) = 2;
    D2x(end-1,end) = 2;
    
    D2x = sigma*D2x * (Nx-1)^2;
    
    %%
    
    Ny = Nx;
    y = x;
    D2y = D2x;

    n = (Nx)^2;
    
    % Grid and initial data:
    [X,Y] = meshgrid(x,y);
    V0 = zeros(size(X));
    S0 = ones(size(X));
    I0 = zeros(size(X));
    R0 = zeros(size(X));
    D0 = zeros(size(X));

    %S0(2:5,2:5) = 1;
    %S0 = 100*exp(-5*(X.^2 + Y.^2) ) + 1;
    S0(5,5) = 1000;
    
    %V0 = exp(-15*(X.^2 + Y.^2) );
    V0(2,2) = 0.1;
    
    IC = flatten(V0,S0,I0,R0,D0);
    [times,states] = ode23(@(t,x)F(x),[0 T],IC);

    figure(1)
    for t = times
        subplot(2,3,1)
        V = states(1:n);
        v = shape(V);
        surf(v);
        title('V');
        axis tight
        subplot(2,3,2)
        S = states(n+1:2*n);
        s = shape(S);
        surf(s);
        title('S');
        axis tight
        subplot(2,3,3)
        I = states(2*n+1:3*n);
        i = shape(I);
        surf(i);
        title('I');
        axis tight
        subplot(2,3,4)
        R = states(3*n+1:4*n);
        r = shape(R);
        surf(r);
        title('R');
        axis tight
        subplot(2,3,5)
        D = states(4*n+1:5*n);
        d = shape(D);
        surf(d);
        title('D');        
        axis tight
    end
    
    function outstate = F(state)

        V = state(1:n);
        v = shape(V);
        S = state(n+1:2*n);
        s = shape(S);
        I = state(2*n+1:3*n);
        i = shape(I);
        R = state(3*n+1:4*n);
        r = shape(R);
        D = state(4*n+1:5*n);
        d = shape(D);

        %v([1 Ny+1],:) = BC*v(2:Ny,:);         
        %w = v';
        %w([1 Ny+1],:) = BC*w(2:Ny,:);        
        %v = w';

        Vt = v*D2x + D2y*v + p*i;
        St = -lambda*s.*v;
        It = -rho*i - dr*i + lambda*s.*v;
        Rt = rho*i;
        Dt = dr*i;

        outstate = flatten(Vt,St,It,Rt,Dt);

    end

    function f = flatten(V,S,I,R,D)
        f = [V(:) ; S(:) ; I(:) ; R(:) ; D(:)];
    end

    function sX = shape(X)
        sX = reshape(X,Nx,Ny);
    end

end