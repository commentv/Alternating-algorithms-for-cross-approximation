function [subsampled_kron_prod] = khatri_rao_product(List_Matrix_subsampled)

    %%%
    % Input : A list of matrices "List_Matrix_subsampled".
    %
    % Output : The subsampled chain of kronecker product "subsampled_kron_prod".
    %%%
    
    % Apply the khatri-rao product to from the subsampled chain of
    % kronecker products from the different columns sampled.

    d = length(List_Matrix_subsampled);
    J = size(List_Matrix_subsampled{1},1);
    R = J;
    K = size(List_Matrix_subsampled{2},2);

    subsampled_kron_prod = reshape(List_Matrix_subsampled{1},[J 1 K]);

    for i = 2:d
        A = reshape(List_Matrix_subsampled{i},[1 R K]);
        subsampled_kron_prod = reshape(bsxfun(@times,A,subsampled_kron_prod),[R*J 1 K]);
        J = J*R;
    end
    subsampled_kron_prod = reshape(subsampled_kron_prod,[J K]);

return