clear,clc

S0 = 36; r = 0.06; sig = 0.2; T = 1; 
M = 2; K = 50;
% paths = pathBS1d(S0, r, sig, T, M, K);
paths = antipathBS1d(S0, r, sig, T, M, K);
showpath(paths, S0, 2*M)