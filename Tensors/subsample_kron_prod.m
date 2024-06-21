function [subsampled_kron_prod,index_set,scaling_array] = subsample_kron_prod(List_Matrix,List_probabilities,tilde_l)

    %%%
    % Inputs : A list of factor matrices "List_Matrix"
    % representing the current factor matrices, a list of probabilites 
    % "List_probabilities" denoting the distribution derived from the
    % leverage scores of the different factor matrices, and an integer 
    % "l_tilde representing the size of the presampling.
    %
    % Outputs : A matrix "subsampled_kron_prod" denoting the subsampled 
    % chain of kronecker products, an index array "index_set" representing
    % the index sampled and an array "scaling_array" which contains the
    % rescaling of the different columns chosen.
    %%%

    d = length(List_Matrix);
    Size_matrices = [];

    %Gather the number of rows of every factor matrices
    for j = 1:d
        Size_matrices = [Size_matrices size(List_Matrix{j},1)];
    end


    index_set = 1;
    scaling_array = 1/sqrt(tilde_l);
    List_Matrix_subsampled = {};
    for j = d:-1:1
        % Sample one factor matrix
        index_sampled = randsample(1:Size_matrices(j),tilde_l,true,List_probabilities{j});
        % Update the index sampled
        index_set = index_set + (index_sampled-1)*prod(Size_matrices(1:(j-1)));
        % Update the scaling factor
        scaling_array = scaling_array./sqrt(List_probabilities{j}(index_sampled));
        % Get the rows sampled
        List_Matrix_subsampled{j} = List_Matrix{j}(index_sampled,:)';
    end
    % Form the final subsampled chain of kronecker products
    subsampled_kron_prod = scaling_array.*khatri_rao_product(List_Matrix_subsampled)';

return
    
