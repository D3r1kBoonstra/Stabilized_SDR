function gamma = ECD1st(M,U,maxiter)
% estimating M-envelope contains span(U)
% where M>0 and is symmetric
% dimension of the envelope is d
% based on inv(M+U) and (M)
gamma = optimECD(M,inv(M+U),ECDini(M,U),maxiter);