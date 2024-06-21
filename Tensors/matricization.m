function X = matricization(X,mu)

%%%
% Input : A tensor "X", an integer "mu" representing the mu-th mode.
%
% Output : A matrix "X" representing the mu-th mode matricization of the 
% input tensor "X".
%%%

dims = size(X);
d = length(dims);

% Matricize X
perm = [mu setdiff(1:d,mu)];
X = permute(X,perm);
X = reshape(X,dims(mu),prod(dims(perm(2:end))));

end