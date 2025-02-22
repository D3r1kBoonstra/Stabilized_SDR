function [u, bic_val] = ENDS_bic(Xn,Yn,umax,lambda)
[n,p]=size(Xn);
K = length(unique(Yn)); % number of classes
nk = zeros(K,1);  %number of obs. in each class
for k=1:K
    nk(k) = sum(Yn==k);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Sample means and covariance matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Xbar = mean(Xn);    %overall mean
Xcn = Xn - Xbar(ones(n,1),:);   %centered predictor
Xcnk = Xn;  %groupwise centered predictor
Xkbar = zeros(p,K); %group mean
Nind = 1:n;
Sk = zeros(p,p,K);
for k=1:K
    idx = Nind(Yn==k);
    Xkbar(:,k) = mean(Xn(idx,:));
    Xcnk(idx,:) = Xcnk(idx,:) - repmat(Xkbar(:,k)',[nk(k),1]);
    Sk(:,:,k) = Xcnk(idx,:)'*Xcnk(idx,:)/nk(k);
end
Sx = Xcn'*Xcn/n;
S = Xcnk'*Xcnk/n;
for k=1:K
    Sk(:,:,k) = lambda*S + (1-lambda)*Sk(:,:,k);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  BIC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bic_max = umax;
bic_val = zeros(bic_max,1);
for d=1:bic_max
    Ghat = ENDS(Xn,Yn,d,lambda);
    bic_val(d) = n*logdet(Ghat'/Sx*Ghat) + (K-1)*log(n)*((1-lambda)*d*(d+3)/2+lambda*d);
    for k=1:K
        bic_val(d) =  bic_val(d) + nk(k)*logdet(Ghat'*Sk(:,:,k)* Ghat);
    end
end
[~, u] = min(bic_val);


