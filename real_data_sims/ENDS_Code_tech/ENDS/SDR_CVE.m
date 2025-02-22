function cverr = SDR_CVE(Xn,Yn,u,nfold)
% Input: 
%   predictor: Xn, n*p, continuous
%   class label: Yn, n*1, positive integers 1,...,K
%   SDR dimension: u, integer between 1 and p
%   number of folds in cross-validation: nfold, {2,3,...,n}
% Output:
%   ENDS cross-validation error: cverr.sirlda, cverr.savelda, cverr.drlda
%   and cverr.sirqda, cverr.saveqda, cverr.drqda
CVO = cvpartition(Yn,'k',nfold);
err.sirlda = zeros(CVO.NumTestSets,1); err.sirqda = err.sirlda;
err.savelda = err.sirlda; err.saveqda = err.sirqda;
err.drlda = err.sirlda; err.drqda = err.sirqda;
for i = 1:CVO.NumTestSets
    trIdx = CVO.training(i);
    teIdx = CVO.test(i);
    Ytest = SDR_classify(Xn(teIdx,:),Xn(trIdx,:),Yn(trIdx,:),u);
    err.sirlda(i) = sum(Yn(teIdx,:)~=Ytest(:,1));
    err.sirqda(i) = sum(Yn(teIdx,:)~=Ytest(:,2));
    err.savelda(i) = sum(Yn(teIdx,:)~=Ytest(:,3));
    err.saveqda(i) = sum(Yn(teIdx,:)~=Ytest(:,4));
    err.drlda(i) = sum(Yn(teIdx,:)~=Ytest(:,5));
    err.drqda(i) = sum(Yn(teIdx,:)~=Ytest(:,6));
    
end
cverr.sirlda = sum(err.sirlda)/sum(CVO.TestSize);
cverr.sirqda = sum(err.sirqda)/sum(CVO.TestSize);
cverr.savelda = sum(err.savelda)/sum(CVO.TestSize);
cverr.saveqda = sum(err.saveqda)/sum(CVO.TestSize);
cverr.drlda = sum(err.drlda)/sum(CVO.TestSize);
cverr.drqda = sum(err.drqda)/sum(CVO.TestSize);

cverr.all = [cverr.sirlda cverr.sirqda cverr.savelda cverr.saveqda ...
cverr.drlda cverr.drqda];
