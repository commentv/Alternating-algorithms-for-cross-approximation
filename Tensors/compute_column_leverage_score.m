function  [probabilities] = compute_column_leverage_score(B, rank)
    
    %%%
    % Inputs : A matrix "B" with full row rank, an integer "rank"
    % representing the rank of the matrix
    %
    % Output : The leverage scores distribution of the matrix "B"
    %%%

    [~,~,V] = svd(B,'econ');
    V_B = V(:,1:rank);
    probabilities = vecnorm(V_B,2,2)./rank;

end