function [A,B] = matricization_Tucker(G,Matrix_List,mu)
    
    %%%
    % Inputs : A core tensor "G", a list of factor matrices "Matrix_List" and
    % an integer "mu" representing the mu-th mode.
    %
    % Output : Two matrices "A" and "B" where "A" represent the mu-th factor
    % matrix and "B" represent the rest of the mu-th mode of the initial
    % tensor in Tucker's format.
    %%%

    dims = size(G);
    d = length(dims);

    % Get the mu-th factor matrix.
    A = Matrix_List{mu};
    % Compute the mu-th mode matricization of the core tensor.
    G_mu = matricization(G,mu);

    indexes = flip(setdiff(1:d,mu));

    % Compute the kronecker product of the remaining factor matrices to
    % form "B".
    right_term = 1;
    for i = 1:(d-1)
        index = indexes(i);
        right_term = kron(right_term,Matrix_List{index});
    end

    B = G_mu*right_term';

return

