function [C,UR] = RRQR_CUR(M,r,iter)

    %%%
    % Inputs : A matrix "M", the rank "r" of its initialization as well as
    % the number of iterations "iter" to perform.
    %
    % Outputs : The intializations "A_0 = C" and "B_0 = UR" for the matrix 
    % iterative refinement algorithm.
    %%%

    [~,n] = size(M,[1,2]);

    %First sample r columns at random
    col_choice = randsample(1:n,r,false);

    for i = 1:iter
        %Perform a qr factorization on the matrix M reduced
        % on the columns chosen and get the rows
        [~,~,row_choice] = qr(M(:,col_choice),"econ","vector");
        row_choice = sort(row_choice(1:r),'ascend');
    
        %Perform a qr factorization on the matrix M reduced
        % on the rows chosen and get the columns
        [~,~,col_choice] = qr(M(row_choice,:),"econ","vector");
        col_choice = sort(col_choice(1:r),'ascend');
    end
    
    %Form the CUR approximation
    C = M(:,col_choice);
    UR = lsqminnorm(M(row_choice,col_choice),M(row_choice,:));
return