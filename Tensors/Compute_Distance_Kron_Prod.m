function [distance] = Compute_Distance_Kron_Prod(Subsampled_kron_prod,List_factor_matrix,mu)

    %%%
    % Inputs : A matrix "Subsampled_kron_prod" representing the 
    % subsampled chain of kronecker products, a list "List_factor_matrix" 
    % containing the different matrices, an integer "mu" denoting the mode
    % of the matricization
    %
    % Output : The distance between the subsampled chain of kronecker
    % products and the exact chain of kronecker products.
    %%%

    List_factor_matrix(mu) = [];

    % Form the exact chain of kronecker products
    kron_factor_matrix = kron_prod_list_matrices(List_factor_matrix);
    % Compute the distance
    [~,~,V_subsampled_kron_prod] = svd(Subsampled_kron_prod,'econ');
    distance = norm(kron_factor_matrix' - ...
        V_subsampled_kron_prod*V_subsampled_kron_prod'*kron_factor_matrix');
end