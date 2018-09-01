%% The aim of this script is to establish whether the full RNA localization
%% model is identifiable or not.
%%
%%%%%%%%%%%%%%%%%%
is_symbolic = 1;

if is_symbolic
    syms a b gamma nu t
    B = construct_matrix(nu,gamma,is_symbolic);
    [V,D] = eig(B);
    y0 = [0; 0*ones(15,1)];
    producers = [0; ones(15,1)];
    % compute steady state distn
    k1 = -(a/b) * (B\producers);
    % incorporate initial conditions
    c = V\(y0 - k1);
    y =  V*expm(D*t)*c + k1;
    diff(y,a)
    
    %%%%%%%%%%%%%%%%%%%%%%%
else
    %set or sample model parameters
    params.a = 1;
    params.b = 0.1;
    params.gamma = 0.01;
    params.nu = 0.72;
    params.phi = 0.5;
    params.sigma = 1;
    
    %set initial conditions and producers
    y0 = [0; 0*ones(15,1)];
    producers = [0; ones(15,1)];
    
    
    %get eigenvalues
    B = construct_matrix(params.nu, params.gamma,is_symbolic);
    [V,D] = eig(B);
    % compute steady state distn
    k1 = -(params.a/params.b) * (B\producers);
    % incorporate initial conditions
    c = V\(y0 - k1);
    
    y = @(t) V*expm(D*t)*c + k1;
    
    nt = 20;
    tt = linspace(7.5,9.5,20);
    y1 = zeros(nt,16);
    for i=1:nt
        y1(i,:) = y(tt(i));
    end
    
    figure;
    plot(tt,y1);
    
    y1
    k1
end


