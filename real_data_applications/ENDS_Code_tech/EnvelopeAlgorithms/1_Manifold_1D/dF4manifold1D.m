function df = dF4manifold1D(W,FParameters)


M = FParameters.M;
U = FParameters.U;
a = M*W/(W'*M*W);
b = (M+U)\W/(W'/(M+U)*W);
% a = M*W*inv(W'*M*W);
% b = inv(M+U)*W*inv(W'*inv(M+U)*W);
df = 2*(a+b);
