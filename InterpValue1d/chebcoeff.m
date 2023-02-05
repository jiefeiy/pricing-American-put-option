function a = chebcoeff(x, fx)
% obtain chebyshev coefficient a_n for discrete chebyshev series
% f(x) = \sum_n a_n T_n(x)
% x and fx should be column vectors
% x is the Chebyshev-Lobatto points 
% eg: -------------------------------------------------------------------
%    n = 10; 
%    range = [0 2];
%    xtest = (range(1,2)-range(1,1))/2*cos( pi*(0:n)'/n ) + (range(1,2)+range(1,1))/2;
% -----------------------------------------------------------------------
n = size(x,1)-1;
fx = fx/(2*n);                            
g = real(fft(fx([1:n+1 n:-1:2])));        % Fast Fourier Transform
a = [g(1); g(2:n)+g(2*n:-1:n+2); g(n+1)]; % Chebyshev coefficients
