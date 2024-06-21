function Q = randomized_range_finder(M,rank)

    %%%
    % Inputs : A matrix "M" and the rank "r" of its initialization as well as
    %
    % Output : The orthonormal basis of the range finder algorithm, that is
    % Q is the orthonormal basis of M*Omega where Omega is a n x rank
    % matrix with independent gaussian entries
    %%%

    n = size(M,2);
    
    %Create the matrix with independent gaussian entries
    Omega = randn([n,rank]);
    %Get the orthonormal basis of M*Omega
    [Q,~] = qr(M*Omega,'econ');


return