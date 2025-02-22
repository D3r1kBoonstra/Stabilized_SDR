function gamma = first_gmanifold1D(Sx,Sk,nk)
maxitera = 50;
nk = nk/sum(nk);

%--- get sample statistics ................................................
data_parameters.Sx = Sx;
data_parameters.Sk = Sk;
data_parameters.nk = nk;


%--- get handle to objective function and derivative ......................
Fhandle = F(@F4gmanifold1D,data_parameters);
dFhandle = dF(@dF4gmanifold1D,data_parameters);
W0 =  get_ini_gmanifold1D(Sx,Sk,nk);
% W0 = orth(rand(p,1));
[~, gamma] = sg_min(Fhandle,dFhandle,W0,'prcg','euclidean',{1:1},'quiet',maxitera);