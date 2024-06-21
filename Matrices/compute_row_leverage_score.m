function  [probabilities] = compute_row_leverage_score(A, rank)
    
    %%%
    % Inputs : A matrix A with full column rank, the rank of the matrix
    %
    % Output : The leverage scores distribution of this matrix
    %%%

    [U,~,~] = svd(A,'econ');
    U_A = U(:,1:rank);
    probabilities = vecnorm(U_A,2,2)./rank;

end
