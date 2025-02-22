function fk = objF_ECD(A,B,w)
fk = log(w'*A*w) + log(w'*B*w) - 2*log(w'*w);
