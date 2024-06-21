function [G,Matrix_List] = create_initialization(init,M,R)

    %%%
    % Inputs : A string "init" denoting the intialization to create, the
    % tensor M from which we derive the initialization and the parameter R
    % denoting the rank we want to approximate.
    %
    % Outputs : The intializations "G" and "Matrix_List" representing
    % respectively the core tensor and the factor matrices of the
    % initialization in Tucker's format.
    %%%

    dims = size(M);
    d = length(dims);
    Matrix_List = {};
    G = M;
    dims_G = dims;
    switch init

        case 'THOSVD'

            % Compute the Truncated High-order SVD
            [G,Matrix_List] = THOSVD(M,R);

        case 'STHOSVD'

            % Compute the Sequentially Truncated High-order SVD
            [G,Matrix_List] = STHOSVD(M,R);

        case 'Sampling'

            % Form an initialization by sampling uniformly at random
            for mu = 1:d

                M_mu = matricization(M,mu);

                number_columns = size(M_mu,2);
                number_rows = size(M_mu,1);

                Indexes_column = randsample(1:number_columns,R,true);
                Indexes_row = randsample(1:number_rows,R,true);

                Matrix_List{mu} = lsqminnorm(M_mu(Indexes_row,Indexes_column)',M_mu(:,Indexes_column)')';

                G_mu = matricization(G,mu);
                dims_G(mu) = R;

                G = tensorization(G_mu(Indexes_row,:),mu,dims_G);
            end

        case 'Range Finder'

            % Use the Range Finder algorithm in the Sequentially Truncated High-order SVD
            [G,Matrix_List] = Range_Finder_HOSVD(M,R);

    end

end