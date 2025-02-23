function f = F4gmanifoldFG(W,FParameters)
Sx = FParameters.Sx;
Sk = FParameters.Sk;
nk = FParameters.nk;
K = size(nk,1);
f = logdet(W'/Sx*W);
for k=1:K
	f = f + nk(k)*logdet(W'*Sk(:,:,k)* W);
end

