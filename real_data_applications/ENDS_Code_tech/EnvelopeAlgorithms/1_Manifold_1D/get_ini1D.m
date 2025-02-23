function W0 = get_ini1D(M,U)
p = size(U,2);
[v1, d1]=eig(M);
[v2, d2]=eig(M+U);
v = [v1 v2];
W0 = v(:,1);
Fw0 = log(W0' /(M+U) * W0)+log(W0' *(M) * W0);
for i=2:(2*p)
    W = v(:,i);
    Fw = log(W' /(M+U) * W)+log(W' *(M) * W);
    if Fw<Fw0
        W0 = W;
        Fw0 = Fw;
    end
end