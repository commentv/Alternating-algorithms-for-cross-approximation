function [A_0,B_0] = create_initialization(init,M,U_r,Sigma_rV_r,r,iter_CUR)
    

    %%%
    % Inputs : A string "init" denoting the intialization to create, the
    % matrix M from which we derive the initialization, the parameters
    % "U_r" and "Sigma_rV_r" representing the truncated SVD (to avoid
    % computing twice the SVD if the truncated initialization is used),
    % the parameter r denoting the rank we want to approximate, and the 
    % parameter "iter_CUR" denoting the number of iterations in the CUR
    % intialization.
    %
    % Outputs : The intializations "A_0" and "B_0" for the matrix iterative
    % refinement algorithm.
    %%%

    [m,n] = size(M,[1,2]);
    
    switch init
        case 'Range Finder'
            % Use the randomized_range_finder.m file to get the orthonormal
            % basis of the range finder algorithm.
            Q = randomized_range_finder(M,r);

            % Form A_0 and B_0
            A_0 = Q;
            B_0 = Q'*M;

        case 'CUR'
            % Use the RRQR_CUR.m file
            [C,UR,~,~] = RRQR_CUR(M,r,iter_CUR);

            % Form A_0 and B_0
            A_0 = C;
            B_0 = UR;

        case 'Truncation'
            % Form A_0 and B_0 via the truncated SVD
            A_0 = U_r;
            B_0 = Sigma_rV_r;

        case 'Sampling'

            % Form A_0 and B_0 by sampling randomly at uniform respectively
            %  r columns and r rows of M.
            indexes_row = randsample(1:m,r,true);
            indexes_column = randsample(1:n,r,true);

            A_0 = M(:,indexes_column);
            B_0 = M(indexes_row,:);
            

    end