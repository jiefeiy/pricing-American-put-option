% discounted to present first
% CC quadrature 
% + barycentric formula at chebyshev knots (roots of first kind cheb polys)
clear,clc
S0 = 36; r = 0.06; sig = 0.2; T = 1; kappa = 40; K = 50; 
M = 100;     % number of interpolation pts
n = 500;     % number of quadrature pts
dt = T/K;
lb = log(S0) + (r-sig*sig/2)*T - 4*sig*sqrt(T);
ub = log(S0) + (r-sig*sig/2)*T + 4*sig*sqrt(T);

x = (ub-lb)/2*cos( (2*(1:M)-1)/M*pi/2 ) + (ub+lb)/2; % interpolation pts x
C = zeros(size(x)); s = zeros(n+1, M);
%% compute quadrature pts s(:)
for m = 1:M
    li = x(m) + (r-sig*sig/2)*dt - 4*sig*sqrt(dt);
    ui = x(m) + (r-sig*sig/2)*dt + 4*sig*sqrt(dt);
    t = pi*(0:n)'/n;
    s(:,m) = (ui-li)/2*cos(t) + (ui+li)/2;
end
V = max(kappa - exp(s(:)), 0);
V = exp(-r*T)*reshape(V, n+1, M);
%% Backward induction
for k = K-1:-1:1
    for m = 1:M
        li = x(m) + (r-sig*sig/2)*dt - 4*sig*sqrt(dt);
        ui = x(m) + (r-sig*sig/2)*dt + 4*sig*sqrt(dt);
        t = pi*(0:n)'/n;
        ftemp = 1/sig/sqrt(2*pi*dt)*exp( -(s(:,m) - x(m) - (r-sig*sig/2)*dt).^2/(2*dt*sig*sig) );
        C(m) = clenshaw_curtis(s(:,m), V(:,m).*ftemp, li, ui);
    end
    f = barycentric_2(x', C', lb, ub, s(:));
    g = exp(-r*k*dt)*(kappa - exp(s(:)));
    V = max([g, f], [], 2);
    V = reshape(V, n+1, M);
    if k == 47
        plot(exp(s(:)), V(:),'.r'); hold on;
        plot(exp(x), C, 'ok'); 
        
    end
    if k == 2
        plot(exp(s(:)), V(:),'.b'); hold on;
        plot(exp(x), C, 'og'); 
        
    end
end
%% compute the initial value
for m = 1:M
    li = x(m) + (r-sig*sig/2)*dt - 4*sig*sqrt(dt);
    ui = x(m) + (r-sig*sig/2)*dt + 4*sig*sqrt(dt);
    t = pi*(0:n)'/n;
    ftemp = 1/sig/sqrt(2*pi*dt)*exp( -(s(:,m) - x(m) - (r-sig*sig/2)*dt).^2/(2*dt*sig*sig) );
    C(m) = clenshaw_curtis(s(:,m), V(:,m).*ftemp, li, ui);
end
V0 = barycentric_2(x', C', lb, ub, log(S0))
