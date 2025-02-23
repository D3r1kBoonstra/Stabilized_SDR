addpath(genpath(['.' filesep]));
rng(20170516)

K = 4; % number of classes
r = 15; % number of responses
u = 2; % dimension of the envelope
n = 300; % total number of observations
ntest = 3000; % always test on 1000 sample
nk = n/K*ones(K,1);   % number of observations for each class
nktest = ntest/K*ones(K,1); % number of observations for each class in test set

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   fixed response
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Yn = reshape(repmat(1:K,[n/K,1]),[n,1]);
Ytest = reshape(repmat(1:K,[ntest/K,1]),[ntest,1]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Set parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Gamma = orth(rand(r,u));
Gamma0 = null(Gamma');

Omega0 = rand(r-u,r-u); Omega0 = Omega0*Omega0'/norm(Omega0*Omega0','fro'); 
Sigmak = zeros(r,r,K);
Sigma = zeros(r,r);
for k=1:K
    Omega = rand(u,u); Omega = exp(-k)*Omega*Omega'/norm(Omega*Omega','fro'); 
    Sigk = Gamma*Omega*Gamma' + Gamma0*Omega0*Gamma0';
    Sigmak(:,:,k) = 25*Sigk/norm(Sigk,'fro');
    Sigma = Sigma + nk(k)/n*Sigmak(:,:,k);
end

muk = zeros(r,K);
for k=1:K
    alphak = normrnd(0,1,[u,1]);
    muk(:,k) = Gamma*alphak;
end
mubar = muk*nk/n;
muk = muk - mubar(:,ones(1,K));
betatrue = Sigma\(muk(:,2:K)-repmat(muk(:,1),[1,K-1]));


nsim = 100;
Err_lda = zeros(nsim,1); Err_qda = zeros(nsim,1); Err_bayes = zeros(nsim,1);
Err_elda = zeros(nsim,1); Err_sirqda = zeros(nsim,1);
Err_sirlda = zeros(nsim,1); Err_saveqda = zeros(nsim,1);
Err_savelda = zeros(nsim,1); Err_drqda = zeros(nsim,1);
Err_drlda = zeros(nsim,1); Err_eqda = zeros(nsim,1);
Err_cvelda = zeros(nsim,1); Err_cveqda = zeros(nsim,1);
Err_nb = zeros(nsim,1); Err_svm = zeros(nsim,1);
lambda_lda = zeros(nsim,1); lambda_qda = zeros(nsim,1);
for isim=1:nsim
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Simulate X|Y and Xtest|Y, which are IID
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Xn = zeros(n,r);
    Xtest = zeros(ntest,r);
    for k=1:K
        idx = [1:nk(k)] + sum(nk(1:(k-1)));
        Xn(idx,:) = mvnrnd(muk(:,k),Sigmak(:,:,k),nk(k));
        idx = [1:nktest(k)] + sum(nktest(1:(k-1)));
        Xtest(idx,:) = mvnrnd(muk(:,k),Sigmak(:,:,k),nktest(k));
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Bayes
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    pik = nk/n;
    Ymat = zeros(ntest,K);
    for k=1:K
        muk_diff = muk(:,k);   
        E_betak = Sigma\muk_diff;
        for itest=1:ntest
            Ymat(itest,k) = log(pik(k)) - ecmnobj(Xtest(itest,:),muk(:,k),Sigmak(:,:,k));
        end
    end
    yhat = ones(ntest,1);
    for i=1:ntest
        [lik,yhat(i)] = max(Ymat(i,:));
    end
    Err_bayes(isim) = sum(yhat~=Ytest);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % LDA/QDA
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    lda = fitcdiscr(Xn,Yn,'DiscrimType','linear');
    yhat = predict(lda,Xtest);
    Err_lda(isim) = sum(yhat~=Ytest);
    
    qda = fitcdiscr(Xn,Yn,'DiscrimType','PseudoQuadratic');
    yhat = predict(qda,Xtest);
    Err_qda(isim) = sum(yhat~=Ytest);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ENDS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    lambda=1;
    Ghat = ENDS(Xn,Yn,u,lambda);
    lda = fitcdiscr(Xn*Ghat,Yn,'DiscrimType','linear');
    yhat = predict(lda,Xtest*Ghat);
    Err_elda(isim) = sum(yhat~=Ytest);
    
    lambda=0;
    Ghat = ENDS(Xn,Yn,u,lambda); 
    qda = fitcdiscr(Xn*Ghat,Yn,'DiscrimType','PseudoQuadratic');
    yhat = predict(qda,Xtest*Ghat);
    Err_eqda(isim) = sum(yhat~=Ytest);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ENDS_CV_lambda
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [lambda, Ghatcv] = ENDS_CV(Xn,Yn,u,0:0.1:1);
    
    lda = fitcdiscr(Xn*Ghatcv.lda,Yn,'DiscrimType','linear');
    yhat = predict(lda,Xtest*Ghatcv.lda);
    Err_cvelda(isim) = sum(yhat~=Ytest);
    lambda_lda(isim) = lambda.lda;
    
    qda = fitcdiscr(Xn*Ghatcv.qda,Yn,'DiscrimType','PseudoQuadratic');
    yhat = predict(qda,Xtest*Ghatcv.qda);
    Err_cveqda(isim) = sum(yhat~=Ytest);
    lambda_qda(isim) = lambda.qda;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SDR
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    yhat = SDR_classify(Xtest,Xn,Yn,u);
    Err_sirlda(isim) = sum(yhat(:,1)~=Ytest);
    Err_sirqda(isim) = sum(yhat(:,2)~=Ytest);
    Err_savelda(isim) = sum(yhat(:,3)~=Ytest);
    Err_saveqda(isim) = sum(yhat(:,4)~=Ytest);
    Err_drlda(isim) = sum(yhat(:,5)~=Ytest);
    Err_drqda(isim) = sum(yhat(:,6)~=Ytest);

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % NB and SVM
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    mNB = fitcnb(Xn,Yn);
    yhat = predict(mNB,Xtest);
    Err_nb(isim) = sum(Ytest~=yhat);

    mSVM = fitcecoc(Xn,Yn);
    yhat = predict(mSVM,Xtest);
    Err_svm(isim) = sum(Ytest~=yhat);

end

mean([Err_bayes Err_nb Err_svm Err_lda Err_elda Err_cvelda Err_qda Err_eqda Err_cveqda Err_sirlda Err_sirqda Err_savelda Err_saveqda Err_drlda Err_drqda]/ntest)
std([Err_bayes Err_nb Err_svm Err_lda Err_elda Err_cvelda Err_qda Err_eqda Err_cveqda Err_sirlda Err_sirqda Err_savelda Err_saveqda Err_drlda Err_drqda]/ntest)/sqrt(nsim)
mean([lambda_lda lambda_qda])
std([lambda_lda lambda_qda])/sqrt(nsim)


