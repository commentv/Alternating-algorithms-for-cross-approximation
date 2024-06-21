function [A,B,error] = iterative_refinement_fast_transition_diagram(M,A,B,tau,r,l)

    %%%
    % Inputs : A matrix "M" representing the matrix we want to approximate, 
    % 2 matrices "A" and "B" representing the crude intial low rank 
    % approximations, an integer "tau" representing the number of
    % iterations, an integer "r" representing the rank of the approximation
    % and an integer "l" representing the size of the sampling.
    %
    % Outputs : The final matrices "A" and "B" representing a low rank
    % approximation of "M", as well as the error "error" for the last 
    % iteration.
    %%%

    [m,n] = size(M,[1,2]);

    for i = 1:tau
        %Compute the proba from the row leverage score of current A
    
        [proba_A] = compute_row_leverage_score(A,r);
        
        %Sample the rows of M and current A
        
        indexes_sampled = randsample(1:m,l,true,proba_A);
        scaling = 1./sqrt(l*proba_A(indexes_sampled));
        sampled_A = scaling.*A(indexes_sampled,:);
        sampled_M = scaling.*M(indexes_sampled,:);

        %Update current B

        B = lsqminnorm(sampled_A,sampled_M);
   
        %Compute the proba from the column leverage score of current B
        
        [proba_B] = compute_column_leverage_score(B,r);

        %Sample the columns of M and current B

        indexes_sampled = randsample(1:n,l,true,proba_B);
        scaling = 1./sqrt(l*proba_B(indexes_sampled));
        sampled_B = scaling'.*B(:,indexes_sampled);
        sampled_M = scaling'.*M(:,indexes_sampled);

        %Update current B

        A = lsqminnorm(sampled_B',sampled_M')';
        
    end

    % Gather the error only at the end of the iterations to accelerate the
    % code
    error = norm(M - A*B,'fro');

return