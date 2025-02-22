function Ghat = manifoldFG(M,U,d,G_ini)
% estimating M-envelope contains span(U)
% where M>0 and is symmetric
% dimension of the envelope is d
% based on inv(M+U) and M

maxitera = 50;
if size(U,1)~=size(U,2) %not a squared matrix
    U = U*U';
end
p = size(M,2);
if d==p
    Ghat = eye(p);
    return;
end
%--- get sample statistics ................................................
data_parameters.M = M;
data_parameters.U = U;

%--- get handle to objective function and derivative ......................
Fhandle = F(@F4manifoldFG,data_parameters);
dFhandle = dF(@dF4manifoldFG,data_parameters);
W0 = G_ini;
[~, gamma] = sg_min(Fhandle,dFhandle,W0,'prcg','euclidean',{1:d},'quiet',maxitera);
Ghat = gamma;
