function cverr = ENDS_CVE(Xn,Yn,u,lambdahat,nfold)
% Finding a basis matrix of the ENDS
% Input: 
%   predictor: Xn, n*p, continuous
%   class label: Yn, n*1, positive integers 1,...,K
%   envelope dimension: u.lda, u.qda, integer between 1 and p
%   tuning parameter: lambda, [0,1], where 1 is LDA, 0 is QDA
%   number of folds in cross-validation: nfold, {2,3,...,n}
% Output:
%   ENDS cross-validation error: cverr.lda, cverr.elda, cverr.qda and cverr.eqda
CVO = cvpartition(Yn,'k',nfold);
err.lda = zeros(CVO.NumTestSets,1); err.elda = err.lda;
err.qda = err.lda; err.eqda = err.lda;
err.knn = err.lda; err.nb = err.lda; err.svm = err.lda;
for i = 1:CVO.NumTestSets
    trIdx = CVO.training(i);
    teIdx = CVO.test(i);
    Ytest = ENDS_classify(Xn(teIdx,:),Xn(trIdx,:),Yn(trIdx,:),u,lambdahat);
    err.elda(i) = sum(Yn(teIdx,:)~=Ytest(:,1));
    err.eqda(i) = sum(Yn(teIdx,:)~=Ytest(:,2));
        
    lda = fitcdiscr(Xn(trIdx,:),Yn(trIdx,:),'DiscrimType','linear');
    yhat = predict(lda,Xn(teIdx,:));
    err.lda(i) = sum(Yn(teIdx,:)~=yhat);
    
    qda = fitcdiscr(Xn(trIdx,:),Yn(trIdx,:),'DiscrimType','PseudoQuadratic');
    yhat = predict(qda,Xn(teIdx,:));
    err.qda(i) = sum(Yn(teIdx,:)~=yhat);  
    
    mKNN = fitcknn(Xn(trIdx,:),Yn(trIdx,:),'NumNeighbors',3);
    yhat = predict(mKNN,Xn(teIdx,:));
    err.knn3(i) = sum(Yn(teIdx,:)~=yhat);
    mKNN = fitcknn(Xn(trIdx,:),Yn(trIdx,:),'NumNeighbors',5);
    yhat = predict(mKNN,Xn(teIdx,:));
    err.knn5(i) = sum(Yn(teIdx,:)~=yhat);
    mKNN = fitcknn(Xn(trIdx,:),Yn(trIdx,:),'NumNeighbors',10);
    yhat = predict(mKNN,Xn(teIdx,:));
    err.knn10(i) = sum(Yn(teIdx,:)~=yhat);
    
    mNB = fitcnb(Xn(trIdx,:),Yn(trIdx,:));
    yhat = predict(mNB,Xn(teIdx,:));
    err.nb(i) = sum(Yn(teIdx,:)~=yhat);

    mSVM = fitcecoc(Xn(trIdx,:),Yn(trIdx,:));
    yhat = predict(mSVM,Xn(teIdx,:));
    err.svm(i) = sum(Yn(teIdx,:)~=yhat);


end
cverr.lda = sum(err.lda)/sum(CVO.TestSize);
cverr.qda = sum(err.qda)/sum(CVO.TestSize);
cverr.elda = sum(err.elda)/sum(CVO.TestSize);
cverr.eqda = sum(err.eqda)/sum(CVO.TestSize);
cverr.knn3 = sum(err.knn3)/sum(CVO.TestSize);
cverr.knn5 = sum(err.knn5)/sum(CVO.TestSize);
cverr.knn10 = sum(err.knn10)/sum(CVO.TestSize);
cverr.nb = sum(err.nb)/sum(CVO.TestSize);
cverr.svm = sum(err.svm)/sum(CVO.TestSize);

cverr.all = [cverr.lda cverr.qda cverr.elda cverr.eqda ...
cverr.knn3 cverr.knn5 cverr.knn10 cverr.nb cverr.svm];
