function  [probabilities] = compute_row_leverage_score(A, rank)
    
    %%%
    % Inputs : A matrix "A" with full column rank, an integer "rank"
    % representing the rank of the matrix
    %
    % Output : The leverage scores distribution of the matrix "A"
    %%%

    [U_A,~,~] = svd(A,"econ");
    U_A = U_A(:,1:rank);
    probabilities = vecnorm(U_A,2,2)./rank;

end
