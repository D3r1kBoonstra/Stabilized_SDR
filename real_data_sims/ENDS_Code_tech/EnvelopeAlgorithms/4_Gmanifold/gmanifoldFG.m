function Ghat = gmanifoldFG(Sx,Sk,nk,d,G_ini)
% Estimating basis for ENDS
% Input: 
%   sample covariance of X: Sx, p*p postive definite
%   conditional sample covariances of X in each group: Sk, p*p*K array of 
%       K matrices each has dimension p*p, postive definite
%   groupsize: nk, K*1 vector of sample sizes in each group
%   envelope dimension: u, integer between 1 and p
%   initial guess of ENDS basis matrix: Ghat, semi-orthogonal p*u
% Output:
%   ENDS basis matrix: Ghat, semi-orthogonal p*u
maxitera = 50;
nk = nk/sum(nk);

%--- get sample statistics ................................................
data_parameters.Sx = Sx;
data_parameters.Sk = Sk;
data_parameters.nk = nk;

%--- get handle to objective function and derivative ......................
Fhandle = F(@F4gmanifoldFG,data_parameters);
dFhandle = dF(@dF4gmanifoldFG,data_parameters);
W0 = G_ini;
[~, gamma] = sg_min(Fhandle,dFhandle,W0,'prcg','euclidean',{1:d},'quiet',maxitera);
Ghat = gamma;
