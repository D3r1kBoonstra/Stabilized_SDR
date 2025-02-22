function f = F4manifold1D(W,FParameters)

M = FParameters.M;
U = FParameters.U;
% a = log(W' * M * W);
% b = log(W' * inv(M+U) * W);
a = log(W'*M*W);
b = log(W'/(M+U)*W);
f = a + b;