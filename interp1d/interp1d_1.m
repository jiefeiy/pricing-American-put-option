% !! quadrature formula has some bugs!!
clear,clc
S0 = 36; r = 0.06; sig = 0.2; T = 1; 
M = 100; K = 50; n = 500;
kappa = 40;
% paths = pathBS1d(S0, r, sig, T, M, K);
% showpath(paths, 100, 50)

dt = T/K;
lb = log(S0) + (r-sig*sig/2)*T - 3*sig*sqrt(T);
ub = log(S0) + (r-sig*sig/2)*T + 3*sig*sqrt(T);

x = (ub-lb)/2*cos( (2*(1:M)-1)/M*pi/2 ) + (ub+lb)/2; % interpolation pts x
C = zeros(size(x)); s = zeros(n, M);
%% compute quadrature pts s(:)
for m = 1:M
    li = x(m) + (r-sig*sig/2)*dt - 3*sig*sqrt(dt);
    ui = x(m) + (r-sig*sig/2)*dt + 3*sig*sqrt(dt);
    t = (2*(1:n)'-1)/n*pi/2;
    s(:,m) = (ui-li)/2*cos(t) + (ui+li)/2;
end
%% first step using BS formula
V = bspde(exp(s(:)),kappa,r,T-dt,T,sig);
V = reshape(V,n,M); 
%% Backward induction
for k = K-2:-1:1
    for m = 1:M
        li = x(m) + (r-sig*sig/2)*dt - 3*sig*sqrt(dt);
        ui = x(m) + (r-sig*sig/2)*dt + 3*sig*sqrt(dt);
        t = (2*(1:n)'-1)/n*pi/2;
        ftemp = 1/sig/sqrt(2*pi*dt)*exp( -(s(:,m) - x(m) - (r-sig*sig/2)*dt).^2/(2*dt*sig*sig) );
        C(m) = exp(-r*dt)*(ui-li)/2*pi/n*( sum( V(:,m).*ftemp.*sin(t)) );
    end
    f = barycentric_2(x', C', lb, ub, s(:));
    V = max([kappa - exp(s(:)), f], [], 2);
    V = reshape(V, n, M);
%     if k == 47
%         plot(exp(s(:)), V(:),'.r'); hold on;
%         plot(exp(x), C, 'ok'); 
%         
%     end
%     if k == 1
%         plot(exp(s(:)), V(:),'.b'); hold on;
%         plot(exp(x), C, 'og'); 
%         
%     end
end
%% compute the initial value
for m = 1:M
    li = x(m) + (r-sig*sig/2)*dt - 3*sig*sqrt(dt);
    ui = x(m) + (r-sig*sig/2)*dt + 3*sig*sqrt(dt);
    t = (2*(1:n)'-1)/n*pi/2;
    ftemp = 1/sig/sqrt(2*pi*dt)*exp( -(s(:,m) - x(m) - (r-sig*sig/2)*dt).^2/(2*dt*sig*sig) );
    C(m) = exp(-r*dt)*(ui-li)/2*pi/n*( sum( V(:,m).*ftemp.*sin(t)) );
end
V0 = barycentric_2(x', C', lb, ub, log(S0))
