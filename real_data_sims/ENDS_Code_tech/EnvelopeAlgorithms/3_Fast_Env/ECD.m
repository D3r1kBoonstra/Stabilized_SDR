function Ghat = ECD(M,U,d)
% estimating M-envelope contains span(U)
% where M>0 and is symmetric
% dimension of the envelope is d
% based on inv(M+U) and (M)
maxiter = 2000;
p = size(M,2);
Mnew = M;
Unew = U;
G = zeros(p,d);
G0 = eye(p);
for k=1:d
    gk = ECD1st(Mnew,Unew,maxiter);
    G(:,k)=G0*gk;
    G0 = null(G(:,1:k)');
    Mnew = G0'*M*G0;
    Unew = G0'*U*G0;
end
Ghat = G;