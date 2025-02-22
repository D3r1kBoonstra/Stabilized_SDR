function W0 = get_ini_gmanifold1D(Sx,Sk,nk)
p = size(Sx,2);
K = size(nk,1);
[v, d1]=eig(Sx);
for k=1:K
    [v2, d2]=eig(Sk(:,:,k));
    v = [v v2];
end
W0 = v(:,1);
FW0 = logdet(W0'/Sx*W0);
for k=1:K
	FW0 = FW0 + nk(k)*logdet(W0'*Sk(:,:,k)* W0);
end
for i=2:(K*p)
    W = v(:,i);
    FW = logdet(W'/Sx*W);
    for k=1:K
        FW = FW + nk(k)*logdet(W'*Sk(:,:,k)* W);
    end
    if FW<FW0
        W0 = W;
        FW0 = FW;
    end
end