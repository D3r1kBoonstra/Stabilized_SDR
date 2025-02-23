%% Envelope Algorithms
% Here we explain the basics of fitting an envelope with different
% algorithms. We explain the fitting of a generic M-envelope of span(U)
% with a symmetric positive definite matrix M>0 and a symmetric positive
% semi-definite matrix U>=0.
%% Load the package
clear all;
cd D:\EnvelopeComputing\EnvelopeAlgorithms; % change directory to the package folder
setpaths;   % load all the functions in the package directory
rng(2016);  % set random seed

%% Simulate the two matrices M and U with an envelope structure
% M and U both have size p-by-p and the envelope has dimension u
p = 20; 
u = 5;	
% randomly generate a semi-orthogonal p-by-u basis matrix (Gamma) for the
% envelope and its orthogonal completion (Gamma0) of dimension p-by-(p-u)
Gamma = orth(rand(p,u));    
Gamma0 = null(Gamma'); 
% randomly generated symmetric positive definite matrices, M and U, to have
% an exact u-dimensional envelope structure
Phi = rand(u,u); Phi = Phi*Phi';   
Omega = rand(u,u); Omega = Omega*Omega';
Omega0 = rand(p-u,p-u); Omega0 = Omega0*Omega0';
M = Gamma*Omega*Gamma' + Gamma0*Omega0*Gamma0';
U = Gamma*Phi*Gamma';
% randomly generate symmetric positive definite matrices, Mhat and Uhat, as
% root-n consistent sample estimators for M and U
n = 200;
X = mvnrnd(zeros(p,1),M,n);
Y = mvnrnd(zeros(p,1),U,n);
Mhat = X'*X/n;
Uhat = Y'*Y/n;


%% Envelope estimation from the 1D algorithm
% The 1D algorithm (Cook and Zhang 2016 JCGS) is the most reliable envelope
% estimation method in the literature. It requires no initial value input
% and is guarenteed to converge at a root-n consistent solution. However,
% the 1D envelope estimators are typically not the maximum likelihood
% estimator (MLE). In specific envelope models, the MLE should be obtained
% through full Grassmannian (FG) envelope estimation where the 1D algorithm
% estimators are often used as an initial value.
% In the 1D algorithm, we are using 500 as the default maximum number of iterations.
% If you want to change this number, please go to the "first1D.m" file.

% computational time and estimation error
tic
Ghat_1D = manifold1D(Mhat,Uhat,u);
toc
disp(norm(Gamma*Gamma'-Ghat_1D*Ghat_1D','fro'));

%% Envelope estimation from the ECD algorithm
% The Envelope Coordinate Descent (ECD) algorithm (Cook and Zhang 2016+,
% manuscript, "Fast Envelope Algorithms") is developed under the 1D 
% algorithm framework. Thus, the ECD algorithm and the 1D algorithm
% share the same statistical properties (Fisher consistency, root-n
% consistency, etc.) The ECD algorithm is much faster than the 1D
% algorithm. In the ECD algorithm, we are using 2000 as the default
% maximum number of iterations.
% If you want to change this number, please go to the "ECD.m" file.

% computational time and estimation error
tic
Ghat_ECD = ECD(Mhat,Uhat,u);
toc
disp(norm(Gamma*Gamma'-Ghat_ECD*Ghat_ECD','fro'));

%% Envelope estimation from full Grassmannian (FG) optimization
% The full Grassmannian (FG) optimization means directly optimize the
% envelope objective function over p-by-u Grassmannian. The core
% optimization part, which is a form of conjugate gradient optimization, 
% is borrowed from the "sg_min" Matlab package, which is inside the "LDR"
% Matlab package. Because the envelope objective function is non-convex and 
% has multiple local minima, the FG optimization is highly sensitive to the
% starting value. We recommand to use FG optimization only when p and u are
% small. Unlike the 1D and the ECD algorithms, the FG optimization requires
% a good (e.g. root-n consistent) initial value to begin with. In practice,
% either the 1D or the ECD estimator can be a good initial value for FG optimization.
% In the FG optimization, we are using 500 as the default
% maximum number of iterations.
% If you want to change this number, please go to the "manifoldFG.m" file.
 
% Using a random initial value
tic
Ghat_FG = manifoldFG(Mhat,Uhat,u,orth(rand(p,u)));
toc
disp(norm(Gamma*Gamma'-Ghat_FG*Ghat_FG','fro'));
% Using 1D algorithm estimator as initial value
tic
Ghat_FG = manifoldFG(Mhat,Uhat,u,Ghat_1D);
toc
disp(norm(Gamma*Gamma'-Ghat_FG*Ghat_FG','fro'));
% Using ECD algorithm estimator as initial value
tic
Ghat_FG = manifoldFG(Mhat,Uhat,u,Ghat_ECD);
toc
disp(norm(Gamma*Gamma'-Ghat_FG*Ghat_FG','fro'));

%% Envelope component screeing (ECS) algorithms
% All the above envelope estimation methods (1D, ECD, FG, among others) are not directly applicable to
% n<p settings, therefore, the Envelope component screeing (ECS) algorithm (Cook and Zhang 2016+,
% manuscript, "Fast Envelope Algorithms") is developed as a pre-screening method
% for envelope estimation in high-dimensional settings. After applying the
% ECS algorithm, the ECD algorithm or the 1D algorithm can be applied to
% get an estimator for the original envelope.

% Generate data matrices Mhat and Uhat with dimension p larger than sample size n.
p = 2000; 
u = 5;	
n = 100;
Gamma = orth(rand(p,u));
Gamma0 = null(Gamma');
Phi = rand(u,u); Phi = Phi*Phi'; Phi = Phi/norm(Phi,'fro');
Omega = rand(u,u); Omega = Omega*Omega'; Omega = Omega/norm(Omega,'fro');
Omega0 = rand(p-u,p-u); Omega0 = Omega0*Omega0'; Omega0 = Omega0/norm(Omega0,'fro');
M = Gamma*(Omega+eye(u))*Gamma' + Gamma0*Omega0*Gamma0';
U = Gamma*Phi*Gamma';
X = mvnrnd(zeros(p,1),M,n);
Y = mvnrnd(zeros(p,1),U,n);
Mhat = X'*X/n;
Uhat = Y'*Y/n;
%%
% ECS1 is the screening method that keeps a user-specified number of components.
k = 50; % reduce the dimension from p to k
[A, Mnew, Unew] = ECS1(Mhat,Uhat,k);    
% A is p-by-k basis matrix for the selected envelope components
% Mnew and Unew are the k-by-k reduced matrices
Ghat_ECS = A*ECD(Mnew,Unew,u);  % p-by-u basis matrix for the envelope
disp(norm(Gamma*Gamma'-Ghat_ECS*Ghat_ECS','fro'));
%%
% ECS2 is the screening method that select the number of components to keep in a data-driven way.
[k, A, Mnew, Unew] = ECS2(Mhat,Uhat,-1/n); 
disp(k) % display how many components were kept after screening
Ghat_ECS = A*ECD(Mnew,Unew,u);
disp(norm(Gamma*Gamma'-Ghat_ECS*Ghat_ECS','fro'));



