function delta = Compute_Subspace_Distance(U,V)
    
    orth_compl_U = null(U');
    delta = norm(orth_compl_U'*V);
    
return