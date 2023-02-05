function V0 = interp1d
% !! quadrature formula has some bugs!!
% pricing American put option with single asset 
% using interpolation of hold value function 

S0 = 36; r = 0.06; sig = 0.2; T = 1; K = 50; 
M = 100;     % number of interpolation pts
n = 500;    % number of quadrature pts
kappa = 40;
% paths = pathBS1d(S0, r, sig, T, M, K);
% showpath(paths, 100, 50)

dt = T/K;
lb = log(S0) + (r-sig*sig/2)*T - 3*sig*sqrt(T);
ub = log(S0) + (r-sig*sig/2)*T + 3*sig*sqrt(T);

x = (ub-lb)/2*cos( (2*(1:M)-1)/M*pi/2 ) + (ub+lb)/2;
C = zeros(size(x)); s = zeros(n, M); V = zeros(n,M);
for m = 1:M
    li = x(m) + (r-sig*sig/2)*dt - 3*sig*sqrt(dt);
    ui = x(m) + (r-sig*sig/2)*dt + 3*sig*sqrt(dt);
    t = (2*(1:n)'-1)/n*pi/2;
    s(:,m) = (ui-li)/2*cos(t) + (ui+li)/2;
    V(:,m) = max([kappa - exp(s(:,m)), zeros(n,1)], [], 2);
end
for k = K-1:-1:1
    for m = 1:M
        li = x(m) + (r-sig*sig/2)*dt - 3*sig*sqrt(dt);
        ui = x(m) + (r-sig*sig/2)*dt + 3*sig*sqrt(dt);
        t = (2*(1:n)'-1)/n*pi/2;
        ftemp = 1/sig/sqrt(2*pi*dt)*exp( -(s(:,m) - x(m) - (r-sig*sig/2)*dt).^2/(2*dt*sig*sig) );
        C(m) = (ui-li)/2*pi/n*( sum( V(:,m).*ftemp.*sin(t)) );
    end
    f = barycentric_2(x', C', lb, ub, s(:));
    V = max([kappa - exp(s(:)), exp(-r*dt)*f], [], 2);
    V = reshape(V, n, M);
    if k == 47
        plot(exp(s(:)), V(:),'.r'); hold on;
        plot(exp(x), C, 'ok'); 
        
    end
    if k == 1
        plot(exp(s(:)), V(:),'.b'); hold on;
        plot(exp(x), C, 'og'); 
        
    end
end
for m = 1:M
    li = x(m) + (r-sig*sig/2)*dt - 3*sig*sqrt(dt);
    ui = x(m) + (r-sig*sig/2)*dt + 3*sig*sqrt(dt);
    t = (2*(1:n)'-1)/n*pi/2;
    ftemp = 1/sig/sqrt(2*pi*dt)*exp( -(s(:,m) - x(m) - (r-sig*sig/2)*dt).^2/(2*dt*sig*sig) );
    C(m) = (ui-li)/2*pi/n*( sum( V(:,m).*ftemp.*sin(t)) );
end
V0 = exp(-r*dt)*barycentric_2(x', C', lb, ub, log(S0));
legend('value function V_{48} at t_{48}', 'interpolation knots of C at t_{48}', ...
    'value function V_{1} at t_{1}', 'interpolation knots of C at t_{1}');
title('Interpolate hold value C, quadrature of value V, K = 50 time steps')
