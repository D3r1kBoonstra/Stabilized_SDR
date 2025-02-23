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
bic_u = zeros(nsim,1);
aic_u = zeros(nsim,1);
for isim=1:nsim
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Simulate X|Y and Xtest|Y, which are IID
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Xn = zeros(n,r);
    Xtest = zeros(ntest,r);
    for k=1:K
        idx = [1:nk(k)] + sum(nk(1:(k-1)));
        Xn(idx,:) = mvnrnd(muk(:,k),Sigma,nk(k));
        idx = [1:nktest(k)] + sum(nktest(1:(k-1)));
        Xtest(idx,:) = mvnrnd(muk(:,k),Sigma,nktest(k));
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Sample means and covariance matrices
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    umax = 10;
    lambda = 0.5;
    
%     Ghat = ENDS(Xn,Yn,u,lambda);
%     subspace(Gamma,Ghat)

    [biciu, ~] = ENDS_bic(Xn,Yn,umax,lambda);
    bic_u(isim) = biciu;
    [aiciu, ~] = ENDS_aic(Xn,Yn,umax,lambda);
    aic_u(isim) = aiciu;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % subspaces
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp(isim);
end



sum(aic_u==u)
sum(bic_u==u)

