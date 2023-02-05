function f = barycentric_1(xs, ys, lb, ub, x)
% Chebyshev interpolation using barycentric interpolation formula
% Chebyshev pts (2nd kind) / Chebyshev-Lobatto points
% -------------------------------------------------------------
% xs     =   Chebyshev pts                            (m-by-1 vector)
% ys     =   function value at interpolation knots    (m-by-1 vector)
% lb     =   lower bound of the interpolant interval          (float)
% ub     =   upper bound of the interpolant interval          (float)
% x      =   the points needed to evaluate                   (vector)
% --------------------------------------------------------------
% f      =   interpolation value at the points x             (vector)
% ---------------------------------------------------------------
m = size(xs,1) - 1; % number of interpolation points-1
% Barycentric interpolation
c =  [1/2; ones(m-1,1); 1/2].*(-1).^((0:m)');
numer = zeros(size(x));
denom = zeros(size(x));
exact = zeros(size(x));
for j = 1:m+1
    xdiff = x-xs(j);
    temp = c(j)./xdiff;
    numer = numer + temp*ys(j);
    denom = denom + temp;
    exact(xdiff==0) = j;
end
f = numer./denom; 
jj = find(exact); f(jj) = ys(exact(jj));
f(find(x<lb | x>ub)) = 0;                   % set extrapolant to zero

