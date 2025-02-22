function [d, A, Mnew, Unew] = ECS2(M,U,C1)
% pre-screening M-envelope of U, based on invM and M+U
% C1 is a negative number, suggested as -1/n in practice

p = size(M,2);
rkM = rank(M);
if rkM>(p-1)
    [v, dd] = eig(M+M');
    [dd, idx] = sort(diag(dd),'descend');
    dd = dd/2;
    v = v(:,idx);
    f = zeros(p,1);
    for i=1:p
        vi = v(:,i);
        ei = dd(i);
        f(i) = log(ei) + log(vi'/(M+U)*vi);
    end
    [fs, idf] = sort(f);
    v = v(:,idf);
    for i=1:p
        %     vi = v(:,1:i);
        %     fc = logdet(vi'*inv(M)*vi) + logdet(vi'*(M+U)*vi);
        fc = sum(fs(i:end));
        if fc > C1
            break;
        end
    end
    %A = v(:,1:i);
    A = v(:,1:(i-1));
else
    [v, dd] = eig(M+M');
    v = v(:,(p-rkM+1):end);
    M1 = v'*M*v;
    U1 = v'*U*v;
    [d, A, M2, U2] = ECS2(M1,U1,C1);
    A = v*A;
end
Mnew = A'*M*A;
Unew = A'*U*A;
d = size(A,2);

