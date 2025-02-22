function [A, Mnew, Unew] = ECS1(M,U,k)
% pre-screening M-envelope of U, based on invM and M+U
% k is the number of components wanted

p = size(M,2);
rkM = rank(M);

if rkM>(p-1)
    [v, dd] = eig(M+M');
    [dd, idx] = sort(diag(dd),'descend');
    dd = dd/2;
    v = v(:,idx);
    
    q = min(p,rkM);
    f = zeros(q,1);
    
    for i=1:q
        vi = v(:,i);
        ei = dd(i);
        f(i) = log(ei) + log(vi'*pinv(M+U)*vi);
    end
    [fs, idf] = sort(f);
    A = v(:,idf(1:k));
else
    [v, dd] = eig(M+M');
    v = v(:,(p-rkM+1):end);
    M1 = v'*M*v;
    U1 = v'*U*v;
    [A, M2, U2] = ECS1(M1,U1,k);
    A = v*A;
end

Mnew = A'*M*A;
Unew = A'*U*A;
