function f = barycentric_2(xs, ys, lb, ub, x)
% Chebyshev interpolation using barycentric interpolation formula
% Chebyshev pts (roots of the 1st kind Chebyshev polynomials)
% -------------------------------------------------------------
% xs     =   Chebyshev pts                            (m-by-1 vector)
% ys     =   function value at interpolation knots    (m-by-1 vector)
% lb     =   lower bound of the interpolant interval          (float)
% ub     =   upper bound of the interpolant interval          (float)
% x      =   the points needed to evaluate                   (vector)
% --------------------------------------------------------------
% f      =   interpolation value at the points x             (vector)
% ---------------------------------------------------------------
m = size(xs,1); % number of interpolation points
% Barycentric interpolation
t = (2*(1:m)-1)/m*pi/2;
c = (-1).^(0:m-1).*sin(t);
numer = zeros(size(x));
denom = zeros(size(x));
exact = zeros(size(x));
for j = 1:m
    xdiff = x-xs(j);
    temp = c(j)./xdiff;
    numer = numer + temp*ys(j);
    denom = denom + temp;
    exact(xdiff==0) = j;
end
f = numer./denom; 
f(find(x<lb | x>ub)) = 0;                   % set extrapolant to zero
jj = find(exact); y(jj) = ys(exact(jj));
