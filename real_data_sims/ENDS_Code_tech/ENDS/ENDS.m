function Ghat = ENDS(Xn,Yn,u,lambda)
% Finding a basis matrix of the ENDS
% Input: 
%   predictor: Xn, n*p, continuous
%   class label: Yn, n*1, positive integers 1,...,K
%   envelope dimension: u, integer between 1 and p
%   tuning parameter: lambda, [0,1], where 1 is LDA, 0 is QDA
% Output:
%   ENDS basis matrix: Ghat, semi-orthogonal p*u

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Fit envelope from ECD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ghat = gmanifold1D(Sx,Sk,nk,u);    % Initialize from 1D
if ((p-u)*u)<(n/2)
    if u>1
        Ghat = gmanifoldFG(S,Sk,nk,u,Ghat);  % Fully optimize for MLE
    end
end 
