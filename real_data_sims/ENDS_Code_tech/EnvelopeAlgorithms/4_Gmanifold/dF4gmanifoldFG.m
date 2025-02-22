function df = dF4gmanifoldFG(W,FParameters)
Sx = FParameters.Sx;
Sk = FParameters.Sk;
nk = FParameters.nk;
nk = nk/sum(nk);
K = size(nk,1);
df = Sx\W/(W'/Sx*W);
% = inv(Sx)*W*inv(W'*inv(Sx)*W);
for k=1:K
    df = df + nk(k)*Sk(k)*W/(W'*Sk(k)* W);
    %  = df + nk(k)*Sk(k)*W*inv(W'*Sk(k)* W);
end
df = 2*df;
