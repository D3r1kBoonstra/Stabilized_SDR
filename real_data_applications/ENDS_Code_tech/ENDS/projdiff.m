function dist_ab = projdiff(A,B)
A = orth(A);
B = orth(B);
dist_ab = norm(A*A'-B*B','fro');
end