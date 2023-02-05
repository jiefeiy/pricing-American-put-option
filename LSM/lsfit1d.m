function fx = lsfit1d(coeff, X)
% compute the value of least square approximation
% ---------------------------------------------------------
% coeff =   the coefficients wrt N basis         (N-by-1 column vector)
% X     =   M data points                           (1-by-M row vector)
% ---------------------------------------------------------
% fx    =   function values                         (1-by-M row vector)
% ---------------------------------------------------------
M = size(X,2);
%% weighted Laguerre polynomials
p0 = exp(-X/2)';
p1 = (1-X').*p0;
p2 = (1-2*X' + X'.^2/2).*p0;
% p3 = 1/6*(-X.^3+9*X.^2 - 18*X + 6)'.*p0;
% p4 = 1/24*(X.^4-16*X.^3+72*X.^2-96*X+24)'.*p0;
% p5 = 1/120*(-X.^5+25*X.^4-200*X.^3+600*X.^2-600*X+120)'.*p0;
% p6 = 1/720*(X.^6-36*X.^5+450*X.^4-2400*X.^3+5400*X.^2-4320*X+720)'.*p0;
% p7 = 1/5040*(-X.^7+49*X.^6-882*X.^5+7350*X.^4-29400*X.^3+105840*X.^2-35280*X+5040)'.*p0;
%% compute pseudoinverse 
A = [ones(M,1), p0, p1, p2];
% A = [p0, p1, p2, p3, p4];
fx = A*coeff; fx = fx';
