function I = DEIM(Q_A_mu)

    %%%
    % Input : A matrix "Q_A_mu" with orthonormal columns.
    %
    % Output : An index set "I" chosen via the discrete empirical
    % interpolation method applied on "Q_A_mu".
    %%%

    R = size(Q_A_mu,2);
    [~,i_max] = max(abs(Q_A_mu(:,1)));
    I = [i_max];

    for i = 2:R
        c = Q_A_mu(I,1:i-1)\Q_A_mu(I,i);
        r = Q_A_mu(:,i)-Q_A_mu(:,1:i-1)*c;
        [~,i_max] = max(abs(r));
        I = [I,i_max];
    end
    
return