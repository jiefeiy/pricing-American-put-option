% discounted to present first
% CC quadrature 
% + barycentric formula at chebyshev-Lobatto pts
clear, close all
S0 = 100; r = 0.05; sig = 0.2; T = 1; kappa = 100; K = 50; DivYield = 0.1;
M = 201;      % number of interpolation pts
n = 200;     % number of quadrature pts
dt = T/K;
NumInterp = 3.89;    % NumInterp standard error of interpolation interval
NumSig = 5;       % NumSig standard error of quadrature interval
lb = log(S0) + (r-DivYield-sig*sig/2)*T - NumInterp*sig*sqrt(T);
ub = log(S0) + (r-DivYield-sig*sig/2)*T + NumInterp*sig*sqrt(T);

x = (ub-lb)/2*cos( pi*(0:M)/M ) + (ub+lb)/2; % interpolation pts x
C = zeros(size(x)); s = zeros(n+1, M+1);
%% compute quadrature pts s(:)
for m = 1:M+1
    li = x(m) + (r-DivYield-sig*sig/2)*dt - NumSig*sig*sqrt(dt);
    ui = x(m) + (r-DivYield-sig*sig/2)*dt + NumSig*sig*sqrt(dt);
    t = pi*(0:n)'/n;
    s(:,m) = (ui-li)/2*cos(t) + (ui+li)/2;
end
V = max(exp(s(:)) - kappa, 0);
V = exp(-r*T)*reshape(V, n+1, M+1);
%% Backward induction
for k = K-1:-1:1
    for m = 1:M+1
        li = x(m) + (r-DivYield-sig*sig/2)*dt - NumSig*sig*sqrt(dt);
        ui = x(m) + (r-DivYield-sig*sig/2)*dt + NumSig*sig*sqrt(dt);
        t = pi*(0:n)'/n;
        ftemp = 1/sig/sqrt(2*pi*dt)*exp( -(s(:,m) - x(m) - (r-DivYield-sig*sig/2)*dt).^2/(2*dt*sig*sig) );
        C(m) = clenshaw_curtis(s(:,m), V(:,m).*ftemp, li, ui);
    end
    f = barycentric_1(x', C', lb, ub, s(:));
    g = exp(-r*k*dt)*(exp(s(:)) - kappa);
    V = max([g, f], [], 2);
    V = reshape(V, n+1, M+1);
%     if k == 48
%         plot(exp(s(:)), V(:),'.r'); hold on;
%         plot(exp(x), C, 'ok'); 
%         
%     end
    if k == 1
        % plot(exp(s(:)), V(:),'.b'); hold on;
        plot(exp(x), C, '-k'); 
        
    end
end
%% compute the initial value
for m = 1:M+1
    li = x(m) + (r-DivYield-sig*sig/2)*dt - NumSig*sig*sqrt(dt);
    ui = x(m) + (r-DivYield-sig*sig/2)*dt + NumSig*sig*sqrt(dt);
    t = pi*(0:n)'/n;
    ftemp = 1/sig/sqrt(2*pi*dt)*exp( -(s(:,m) - x(m) - (r-DivYield-sig*sig/2)*dt).^2/(2*dt*sig*sig) );
    C(m) = clenshaw_curtis(s(:,m), V(:,m).*ftemp, li, ui);
end
V0 = barycentric_1(x', C', lb, ub, log(S0))


%% show European call price
plot(exp(x),C, 'Color',	"#77AC30", 'LineWidth',1, 'DisplayName', 'American call');

V_noDiv = bsCallDividend(exp(x),kappa,r,0,T,sig,DivYield);
hold on; plot(exp(x), max(exp(x) - kappa,0), '-r', 'DisplayName', 'intrinsic value');
hold on; plot(exp(x), V_noDiv, '-b', 'LineWidth',1, 'DisplayName', 'European call');
hold on; plot(exp(x), max(exp(-DivYield*T)*exp(x) - kappa*exp(-r*T)*ones(size(x)), 0),'--b', 'DisplayName', 'European asympototic');
legend('Location','northwest');
title('Comparison of American call and European call with dividends');
