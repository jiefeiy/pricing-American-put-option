clear, close all
% test pathBS1d
S0 = 36; r = 0.06; sig = 0.2; T = 1; K = 50;
num_paths = 3;
paths = pathBS1d(S0, r, sig, T, num_paths, K);
figure(1);
showpath(paths, S0, num_paths); % show paths

% test antipathBS1d
antipaths = antipathBS1d(S0, r, sig, T, num_paths, K);
figure(2);
showpath(antipaths, S0, 2*num_paths); % show paths and antipaths