function A_PD = ensurePositiveDefinite(A)
    if isPositiveDefinite(A)
        A_PD = A;
    else
        A_PD = nearestSPD(A);
    end
end
