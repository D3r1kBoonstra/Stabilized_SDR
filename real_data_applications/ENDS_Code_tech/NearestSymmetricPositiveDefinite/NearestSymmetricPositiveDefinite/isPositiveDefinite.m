function isPD = isPositiveDefinite(A)
    % Check if matrix A is positive definite
    [~,p] = chol(A);
    isPD = (p == 0);
end