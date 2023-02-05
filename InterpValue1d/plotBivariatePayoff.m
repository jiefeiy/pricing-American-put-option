clear,close all
[X,Y] = meshgrid(50:1:150, 50:1:150);
G = max(X - 100, 0).*(X<140).*(Y<140);
surf(X,Y,G);
xlabel('asset price S_1'); ylabel('asset price S_2');
zlabel('payoff');
title('payoff of bivariate barrier option, strike price = 100, b_1 = b_2 = log(140)');