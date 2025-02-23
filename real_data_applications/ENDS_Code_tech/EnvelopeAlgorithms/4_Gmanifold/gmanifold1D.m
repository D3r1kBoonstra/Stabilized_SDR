function Ghat = gmanifold1D(Sx,Sk,nk,d)
% 1D Estimating basis for ENDS
% Input: 
%   sample covariance of X: Sx, p*p postive definite
%   conditional sample covariances of X in each group: Sk, p*p*K array of 
%       K matrices each has dimension p*p, postive definite
%   groupsize: nk, K*1 vector of sample sizes in each group
%   envelope dimension: u, integer between 1 and p
% Output:
%   ENDS basis matrix: Ghat, semi-orthogonal p*u
p = size(Sx,2);
nk = nk/sum(nk);
K = size(nk,1);

if d<(p-1)
    Sxnew = Sx;
    Sknew = Sk;
    G = zeros(p,d);
    G0 = eye(p);
    for k=1:d
        gk = first_gmanifold1D(Sxnew,Sknew,nk);
        G(:,k)=G0*gk;
        G0 = null(G(:,1:k)');
        Sxnew = G0'*Sx*G0; 
        for j=1:K
            Sknew(1:end-1,1:end-1,j) = G0'*Sk(:,:,j)*G0;
        end
        Sknew(end,:,:) = [];
        Sknew(:,end,:) = [];    
    end
    Ghat = G;
elseif d<p
    Sxnew = Sx;
    Sknew = Sk;
    G = zeros(p,d);
    G0 = eye(p);
    for k=1:(p-2)
        gk = first_gmanifold1D(Sxnew,Sknew,nk);
        G(:,k)=G0*gk;
        G0 = null(G(:,1:k)');
        Sxnew = G0'*Sx*G0; 
        for j=1:K
            Sknew(1:end-1,1:end-1,j) = G0'*Sk(:,:,j)*G0;
        end
        Sknew(end,:,:) = [];
        Sknew(:,end,:) = [];    
    end
    k = p-1;
    gk = get_ini_gmanifold1D(Sxnew,Sknew,nk);
    G(:,k)=G0*gk;
    Ghat = G;
else
    Ghat = eye(p);
end
