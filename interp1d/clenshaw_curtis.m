function I = clenshaw_curtis(x,fx,xmin,xmax) % (n+1)-pt C-C quadrature
n = size(x,1)-1;
scale = (xmax - xmin)/2;
fx = fx/(2*n);                            
g = real(fft(fx([1:n+1 n:-1:2])));        % Fast Fourier Transform
a = [g(1); g(2:n)+g(2*n:-1:n+2); g(n+1)]; % Chebyshev coefficients
w = 0*a'; w(1:2:end) = 2./(1-(0:2:n).^2); % weight vector
I = scale*w*a;                            % the integral

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% li = 2; ui = 5.5;
% n = 50;
% t = pi*(0:n)'/n;
% p = ui-li; q = ui+li;
% x = (p*cos(t) + q)/2;
% fx = sin(x);
% I_exact = integral(@(x) sin(x), li, ui);
% I = clenshaw_curtis(x, fx, li, ui);
% disp(I_exact - I);


