function [A_mu] = one_step_iterative_refinement(M_mu,number_columns,B_mu,R,l)
    
    %%%
    % Input : A matrix "M_mu" representing the mu-th matriciation of the
    % tensor M, an integer "number_columns" giving the number of columns of
    % "M_mu" and "B_mu", a matrix "B_mu" representing the matrix from which
    % we sample the leverage scores, an integer "R" giving the rank of the 
    % approximation and an integer "l" denoting the number of samplings.
    %
    % Output : An approximation "A_mu" after one step of the iterative 
    % refinement.
    %%%

    [proba_B_mu] = compute_row_leverage_score(B_mu',R);

    % Sample the columns of M and current B_mu using the column leverage
    % score of current B_mu
    indexes_sampled = randsample(1:number_columns,l,true,proba_B_mu);
    scaling = 1./sqrt(l*proba_B_mu(indexes_sampled));
    sampled_B_mu = scaling'.*B_mu(:,indexes_sampled);
    sampled_M_mu = scaling'.*M_mu(:,indexes_sampled);

    %Update current A_mu

    A_mu = lsqminnorm(sampled_B_mu',sampled_M_mu')';
return