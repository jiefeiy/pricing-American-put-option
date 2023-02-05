function fx = leastSquare(X,Y)
% ------------------------------------------------------------
% f(X) = E[Y|X]
% compute M samples of the random variable f(X)
% given M samples of R.V.s X and Y
% (or say, compute regression function values of the data points 
% using least square method combine with Monte Carlo)
% polynomial basis of total degree 2 are used as basis
% ------------------------------------------------------------
% X    =   M data points with d dimension          (d-by-M matrix)
% Y    =   M observations wrt the data points      (1-by-M vector)
% ------------------------------------------------------------
% fx   =   M conditional expectation samples       (1-by-M vector)
% ------------------------------------------------------------
d = size(X,1); M = size(X,2);
p0 = ones(M,1);
p1 = 1-X';
p2 = ( (3-X').*p1 - 1 )/2;
p1p1 = [];
for i = 1:d-1
    for j = i+1:d
        p1p1 = [p1p1, p1(:,i).*p1(:,j)];
    end
end
A = [p0, p1, p2, p1p1]; % total degree 2 polynomials
% A is an M-by-N matrix which stores N basis evaluated on M data points
coeff = (A'*A)\(A'*Y');
% error = max( abs(Y' - A*coeff) );
fx = A*coeff; fx = fx';



