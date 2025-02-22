function f = F4manifoldFG(W,FParameters)
M = FParameters.M;
U = FParameters.U;
a = logdet(W' * M * W);
b = logdet(W'/(M+U) * W);
f = a + b;