function Ytest = ENDS_classify(Xtest,Xn,Yn,u,lambdahat)
% Finding a basis matrix of the ENDS
% Input:
%   predictor (testing sample): Xtest, ntest*p, continuous
%   predictor (training sample): Xn, n*p, continuous
%   class label (training sample): Yn, n*1, positive integers 1,...,K
%   envelope dimension: u.lda, u.qda, integer between 1 and p
%   tuning parameter: lambda, [0,1], where 1 is LDA, 0 is QDA
% Output:
%   predicted class label (testing sample): Ytest, ntest*2, first column is
%   based on LDA (lambda=1), second column is based on QDA (lambda=lambdahat)

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Envelope LDA prediction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k=1:K
    Sk(:,:,k) = S;
end
Ghat = gmanifold1D(Sx,Sk,nk,u.lda);    % Initialize from 1D
if ((p-u.lda)*u.lda)<min(nk)
    if u.lda>1
        Ghat = gmanifoldFG(S,Sk,nk,u.lda,Ghat);  % Fully optimize for MLE
    end
end
lda = fitcdiscr(Xn*Ghat,Yn,'DiscrimType','linear');
Ytest = predict(lda,Xtest*Ghat);
Ytest = [Ytest Ytest];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Envelope QDA prediction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k=1:K
    Sk(:,:,k) = lambdahat*S + (1-lambdahat)*Sk(:,:,k);
end
Ghat = gmanifold1D(Sx,Sk,nk,u.qda);    % Initialize from 1D
if ((p-u.qda)*u.qda)<min(nk)
    if u.qda>1
        Ghat = gmanifoldFG(S,Sk,nk,u.qda,Ghat);  % Fully optimize for MLE
    end
end
qda = fitcdiscr(Xn*Ghat,Yn,'DiscrimType','PseudoQuadratic');
Ytest(:,2) = predict(qda,Xtest*Ghat);
