function X = tensorization(X_mu,mu,dims)

    %%%
    % Inputs : A matrix "X_mu", an integer "mu" representing the mode of 
    % the tensor and the dimensions of the tensor "dims".
    %
    % Output : A tensor "X" which has been transformed back from its mu-th 
    % matricization "X_mu".
    %%%

d = length(dims);
perm = [mu setdiff(1:d,mu)];
perminv(perm) = 1:d;
dims = dims(perm);
X = reshape(X_mu,dims);
X = permute(X,perminv);
end