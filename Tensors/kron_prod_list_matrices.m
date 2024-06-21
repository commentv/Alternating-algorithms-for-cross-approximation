function result = kron_prod_list_matrices(List_Matrices)
    
    %%%
    % Input : A  list of factor matrices "Lsit_Matrix"
    %
    % Output : The chain of kronecker products of the factor matrices in 
    % "List_Matrix"
    %%%

    d = length(List_Matrices);
    result = List_Matrices{d};

    for mu = d-1:-1:1
        result = kron(result,List_Matrices{mu});
    end

return