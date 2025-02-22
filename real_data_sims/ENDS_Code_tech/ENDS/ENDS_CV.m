function [lambda, Ghat] = ENDS_CV(Xn,Yn,u,listlambda)
% Finding the lambda and the basis matrix of the ENDS
% based on nfold CV classification error rate of RDA
% Input: 
%   predictor: Xn, n*p, continuous
%   class label: Yn, n*1, positive integers 1,...,K
%   envelope dimension: u, integer between 1 and p
%   number of folds in cross-validation: nfold, integer between 2 and n
%   a list of tuning parameters (listed in increasing order): listlambda, all values between [0,1]
% Output:
%   Tuning parameter: lambda, chosen from listlambda, where 1 is LDA, 0 is QDA
%   ENDS basis matrix: Ghat, semi-orthogonal p*u
%   lambda.lda, Ghat.lda for ENDS_lda
%   lambda.qda, Ghat.qda for ENDS_qda

nsim = length(listlambda);

lambda.lda = listlambda(1);
lambda.qda = listlambda(1);
G = ENDS(Xn,Yn,u,lambda.lda);
y = num2str(Yn);
c = cvpartition(y,'k',5);
classf = @(XTRAIN, ytrain,XTEST)(classify(XTEST,XTRAIN,ytrain,'linear'));
rate_lda1 = crossval('mcr',Xn*G,y,'predfun',classf,'partition',c);
classf = @(XTRAIN, ytrain,XTEST)(classify(XTEST,XTRAIN,ytrain,'quadratic'));
rate_qda1 = crossval('mcr',Xn*G,y,'predfun',classf,'partition',c);
Ghat.lda = G;
Ghat.qda = G;

for k = 2:nsim
    lambdak = listlambda(k);
    G = ENDS(Xn,Yn,u,lambdak);
    y = num2str(Yn);
    c = cvpartition(y,'k',5);
    classf = @(XTRAIN, ytrain,XTEST)(classify(XTEST,XTRAIN,ytrain,'linear'));
    rate_lda = crossval('mcr',Xn*G,y,'predfun',classf,'partition',c);
    if rate_lda<=rate_lda1   %when there's a tie, we choose the larger lambda (i.e., towards LDA)
        rate_lda1 = rate_lda;
        Ghat.lda = G;
        lambda.lda = lambdak;
    end
    classf = @(XTRAIN, ytrain,XTEST)(classify(XTEST,XTRAIN,ytrain,'quadratic'));
    rate_qda = crossval('mcr',Xn*G,y,'predfun',classf,'partition',c);
    if rate_qda<=rate_qda1
        rate_qda1 = rate_qda;
        Ghat.qda = G;
        lambda.qda = lambdak;
    end
end
