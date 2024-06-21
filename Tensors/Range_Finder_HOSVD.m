function [G,Matrix_List] = Range_Finder_HOSVD(M,R)
    G = M;
    d = length(size(M));
    Matrix_List = {};
    for mu = 1:d
        M_mu = matricization(M,mu);
        M_mu = M_mu*randn([size(M_mu,2),R]);
        [Q_mu,~] = qr(M_mu,'econ');
        G = tensorproduct(G,Q_mu',mu);
        Matrix_List{mu} = Q_mu;
    end
end