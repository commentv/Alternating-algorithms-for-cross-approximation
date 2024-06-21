function [G,List_Matrix,error] = sublinear_iterative_refinement_rescaled_fast_infl_l(M,G,List_Matrix,tau,R,l,l_tilde)

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
    % the error "error" at the end of the iterations.
    %%%

    dims = size(M);
    d = length(dims);
    List_probabilities = {};
    List_Matricization_M = {};

    % Compute the probability distributions for the presampling
    for mu = 1:d
        [List_probabilities{mu}] = compute_row_leverage_score(List_Matrix{mu},R);
        List_Matricization_M{mu} = matricization(M,mu);
        
    end

    [A_1,B_1] = matricization_Tucker(G,List_Matrix,1);
    
    % Compute the initial error
    error = [norm(List_Matricization_M{1} - A_1*B_1,'fro')];

    for i = 1:tau
        for mu = d:-1:1
            % Subsample the chain of kronecker products
            List_probabilities_temp = List_probabilities;
            List_probabilities_temp(mu) = [];
            List_matrix_temp = List_Matrix;
            List_matrix_temp(mu) = [];
            [subsampled_kron_prod,index_set,scaling_array] = subsample_kron_prod(List_matrix_temp,List_probabilities_temp,l_tilde);
             
            % Update the tensor we want to approximate from the presampling
            M_mu = List_Matricization_M{mu};
            Subsampled_M_mu = scaling_array'.*M_mu(:,index_set);
            number_columns = size(Subsampled_M_mu,2);
            G_mu = matricization(G,mu);
            B_mu = G_mu*subsampled_kron_prod';
            
            % Update the mu-th factor matrix
            A_mu = one_step_iterative_refinement(Subsampled_M_mu,number_columns,B_mu,R,l);
            List_Matrix{mu} = A_mu;
            [List_probabilities{mu}] = compute_row_leverage_score(List_Matrix{mu},R);
    
        end
        
        % Obliquely project the initial tensor and update the core tensor 
        % "G" and the list of factor matrices "List_Matrix"
        [G,List_Matrix] = oblique_projection(M,List_Matrix);

        [A_1,B_1] = matricization_Tucker(G,List_Matrix,1);
        
        % Gather the error
        error = [error norm(List_Matricization_M{1} - A_1*B_1,'fro')];

    end

return