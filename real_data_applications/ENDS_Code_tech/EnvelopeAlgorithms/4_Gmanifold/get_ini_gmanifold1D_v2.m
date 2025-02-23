function W0 = get_ini_gmanifold1D_v2(Sx,Sk,nk)
p = size(Sx,2);
K = size(nk,1);
W0 = orth(rand(p,1));
FW0 = -logdet(W0'*Sx*W0);
for k=1:K
	FW0 = FW0 + nk(k)*logdet(W0'*Sk(:,:,k)* W0);
end
nsim = 100;
for i=1:nsim
    W = orth(rand(p,1));
    FW = -logdet(W'*Sx*W);
    for k=1:K
        FW = FW + nk(k)*logdet(W'*Sk(:,:,k)* W);
    end
    if FW<FW0
        W0 = W;
        FW0 = FW;
    end
end
