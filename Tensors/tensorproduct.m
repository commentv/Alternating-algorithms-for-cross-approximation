function Y = tensorproduct(X,A,mu)
    
    %%%
    % Inputs : A tensor "X", a matrix "A" and an integer "mu" representing
    % The mode of the tensor-matrix product.
    %
    % Output : A tensor "Y" which corresponds to the tensor "X" multiplied 
    % by the matrix "A" on its mu-th mode matricization.
    %%%

    % Save important parameters
    dims = size(X);
    
    %Matricize X
    X_mu = matricization(X,mu);
    
    %Product
    Y_mu = A*X_mu;
    
    dims(mu) = size(Y_mu,1);
    %Tensorize Y
    Y = tensorization(Y_mu,mu,dims);
return