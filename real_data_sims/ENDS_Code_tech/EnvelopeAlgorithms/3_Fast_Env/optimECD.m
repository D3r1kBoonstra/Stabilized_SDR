function w = optimECD(A,B,w0,maxiter)
p = length(w0);
epsilon = 1e-8;
[Gp, dn] = eig(A+A');
[dn, idx] = sort(diag(dn),'descend');
Gp = Gp(:,idx);
dn = diag(dn/2);
v0 = Gp'*w0;
GBG = Gp'*B*Gp;
fk = objF_ECD(dn,GBG,v0);
v = v0;
for itera=1:maxiter
    flg = 0;
    alpha = 1/(v'*dn*v);
    beta = 1/(v'*GBG*v);
    delta = 1/(v'*v);
    A1 = alpha*dn;
    B1 = beta*GBG;
    for j = 1:p
        AB1 = A1(j,j) + B1(j,j);
        if (2*delta - AB1)~=0
            v(j) = (v'*A1(:,j) + v'*B1(:,j) - AB1*v(j) )/(2*delta - AB1);
            flg = flg+1;
        end
        if objF_ECD(dn,GBG,v) > (objF_ECD(dn,GBG,v0) + epsilon)
            v = v0;
            flg = flg-1;
        end
    end
%     if flg==0
%         data_parameters.A = dn;
%         data_parameters.B = GBG;
%         %--- get handle to objective function and derivative ......................
%         Fhandle = F(@F4mani1d,data_parameters);
%         dFhandle = dF(@dFmani1d,data_parameters);
%         W0 = v;
%         % W0 = orth(rand(p,1));
%         [fn v] = sg_min(Fhandle,dFhandle,W0,'prcg','euclidean',{1:1},'quiet',50);
%     end
    
    fk1 = objF_ECD(dn,GBG,v);
    if abs(fk-fk1)<epsilon, break, end
    fk = fk1;
end
w = Gp*v;
w = w/sqrt(w'*w);

