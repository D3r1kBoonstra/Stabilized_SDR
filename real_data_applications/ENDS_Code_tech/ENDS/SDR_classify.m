function Ytest = SDR_classify(Xtest,Xn,Yn,u)
% Input:
%   predictor (testing sample): Xtest, ntest*p, continuous
%   predictor (training sample): Xn, n*p, continuous
%   class label (training sample): Yn, n*1, positive integers 1,...,K
%   SDR dimension: u, integer between 1 and p
% Output:
%   predicted class label (testing sample): Ytest, ntest*2, first column is
%   based on LDA (lambda=1), second column is based on QDA (lambda=lambdahat)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   SDR subspaces
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~,bsir]=SIR(Yn,Xn,'disc',u);
[~,bsave]=SAVE(Yn,Xn,'disc',u);
[~,bdr]=DR(Yn,Xn,'disc',u);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   SIR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lda = fitcdiscr(Xn*bsir,Yn,'DiscrimType','linear');
Ytest = predict(lda,Xtest*bsir);
Ytest = repmat(Ytest,[1,6]);
qda = fitcdiscr(Xn*bsir,Yn,'DiscrimType','PseudoQuadratic');
Ytest(:,2) = predict(qda,Xtest*bsir);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   SAVE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lda = fitcdiscr(Xn*bsave,Yn,'DiscrimType','linear');
Ytest(:,3) = predict(lda,Xtest*bsave);
qda = fitcdiscr(Xn*bsave,Yn,'DiscrimType','PseudoQuadratic');
Ytest(:,4) = predict(qda,Xtest*bsave);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   DR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lda = fitcdiscr(Xn*bdr,Yn,'DiscrimType','linear');
Ytest(:,5) = predict(lda,Xtest*bdr);
qda = fitcdiscr(Xn*bdr,Yn,'DiscrimType','PseudoQuadratic');
Ytest(:,6) = predict(qda,Xtest*bdr);
