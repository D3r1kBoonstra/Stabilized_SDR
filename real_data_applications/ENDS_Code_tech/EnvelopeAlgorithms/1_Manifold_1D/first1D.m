function gamma = first1D(M,U)
%--- get sample statistics ................................................
data_parameters.M = M;
data_parameters.U = U;

%--- get handle to objective function and derivative ......................
Fhandle = F(@F4manifold1D,data_parameters);
dFhandle = dF(@dF4manifold1D,data_parameters);
W0 = get_ini1D(M,U);
% W0 = orth(rand(p,1));
[~, gamma] = sg_min(Fhandle,dFhandle,W0,'prcg','euclidean',{1:1},'quiet',200);