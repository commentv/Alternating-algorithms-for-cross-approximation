function [G,List_Matrix,List_distances] = sublinear_iterative_refinement_rescaled_fast_distance_col(M,G,List_Matrix,tau,R,l,l_tilde)

    %%%
    % Inputs : A tensor "M" representing the tensor we want to approximate, 
    % A core tensor "G" and a list of factor matrices "List_Matrix"
    % representing the crude intial low rank tensor approximation, an 
    % integer "tau" representing the number of an integer "R" representing 
    % the rank of the approximation ,an integer "l" representing the size 
    % of the sampling and an integer "l_tilde representing the size of the 
    % presampling.
    %
    % Outputs : A core tensor "G" and a list of factor matrices "List_Matrix"
    % representing the approximation of "M" in Tucker's format as well as 
    % the error "error" of the subsampled chain of kronecker products for 
    % every iterations.
    %%%

    dims = size(M);
    d = length(dims);
    List_probabilities = {};
    List_Matricization_M = {};
    List_distances = {};

    % Compute the probability distributions for the presampling
    for mu = 1:d
        [List_probabilities{mu}] = compute_row_leverage_score(List_Matrix{mu},R);
        List_Matricization_M{mu} = matricization(M,mu);
        List_distances{mu} = zeros([1,tau]);
        
    end

    for i = 1:tau
        for mu = d:-1:1

            % Subsample the chain of kronecker products
            List_probabilities_temp = List_probabilities;
            List_probabilities_temp(mu) = [];
            List_matrix_temp = List_Matrix;
            List_matrix_temp(mu) = [];
            [subsampled_kron_prod,index_set,scaling_array] = subsample_kron_prod(List_matrix_temp,List_probabilities_temp,l_tilde);
            
            % Compute the error between the subsampled chain of kronecker
            % products and the original one
            List_distances{mu}(i) = Compute_Distance_Kron_Prod(subsampled_kron_prod,List_Matrix,mu);
            
            % Update the tensor we want to approximate from the presampling
            M_mu = List_Matricization_M{mu};
            Subsampled_M_mu = scaling_array'.*M_mu(:,index_set);
            G_mu = matricization(G,mu);
            B_mu = G_mu*subsampled_kron_prod';
            
            % Update the mu-th factor matrix
            List_Matrix{mu} = one_step_iterative_refinement(Subsampled_M_mu,l_tilde,B_mu,R,l);
            [List_probabilities{mu}] = compute_row_leverage_score(List_Matrix{mu},R);
    
        end

        % Obliquely project the initial tensor and update the core tensor 
        % "G" and the list of factor matrices "List_Matrix"
        [G,List_Matrix] = oblique_projection(M,List_Matrix);

    end

return