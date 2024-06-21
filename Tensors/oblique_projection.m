function [G,Matrix_List_oblique] = oblique_projection(M,Matrix_List)
    
    %%%
    % Input : A core tensor "M", a list of factor matrices "Matrix_List".
    %
    % Output : A core tensor "G" and a list of factor matrices
    % "Matrix_List_oblique" resulting from the oblique projection of the
    % initial tensor "M" with respect to the list of factor matrices
    % "Matrix_List".
    %%%

    dims = size(M);
    d = length(dims);
    List_Index = {};

    for mu = 1:d
        
        [Matrix_List_oblique{mu},~] = qr(Matrix_List{mu},'econ');
        
        % Compute the indexes of the rows chosen via the DEIM method
        I_mu = DEIM(Matrix_List_oblique{mu});

        List_Index{mu} = I_mu;
        
        % Update the corresponding factor matrix
        Matrix_List_oblique{mu} = (Matrix_List_oblique{mu}(I_mu,:)'\Matrix_List_oblique{mu}')';

    end

    % Update the core tensor
    G = M(List_Index{:});


return