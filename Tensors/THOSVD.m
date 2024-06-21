function [G,Matrix_List] = THOSVD(M,R)
    
    %%%  
    % Input : A tensor "M" and the target rank "R".
    %
    % Output : The core tensor "G" and the list of factor matrices "Matrix_List"
    % Computed from the well known Truncated High-Order SVD.
    %%%


    d = length(size(M));
    Matrix_List = {};
    G = M;
    for mu = 1:d
        M_mu = matricization(M,mu);
        [U_mu,~,~] = svd(M_mu,'econ');
        U_mu = U_mu(:,1:R);
        G = tensorproduct(G,U_mu',mu);
        Matrix_List{mu} = U_mu;
    end
end